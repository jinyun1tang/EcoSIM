module SnowBalanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL
  use abortutils, only : endrun
  use minimathmod, only : AZMAX1,isclose,AZMIN1  
  USE SnowDataType
  use GridConsts
  use GridDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use SoilHeatDataType
  use ClimForcDataType
!  use TFlxTypeMod
  use EcosimConst
  USE SoilPhysDataType
  USE EcoSimSumDataType
  USE FlagDataType
  use SnowPhysData
  use AqueChemDatatype
  use SoilBGCDataType
  use ChemTranspDataType  
  use UnitMod, only : units
implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: SnowMassUpdate     !called in redist
  public :: SnowpackLayering   !called in redist after SnowMassUpdate
  public :: ZeroSnowArrays
  contains

  subroutine SnowMassUpdate(I,J,NY,NX)
  !
  !Description: update snow mass
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: VOLSWI,ENGYW
  integer :: L
      !
      !     CALCULATE SNOWPACK TEMPERATURE FROM ITS CHANGE
      !     IN HEAT STORAGE
      !

  VOLSWI=0.0_r8
  !update snow layer from top to bottom
  
  D9780: DO L=1,JS

    IF(VLHeatCapSnow_col(L,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) .and. L.EQ.1)THEN
      IF(abs(SnoXfer2SnoLay(L,NY,NX))>0._r8)THEN
        CumSno2SnowLay(L,NY,NX)=CumSno2SnowLay(L,NY,NX)+SnoXfer2SnoLay(L,NY,NX)
        CumWat2SnowLay(L,NY,NX)=CumWat2SnowLay(L,NY,NX)+WatXfer2SnoLay(L,NY,NX)
        CumIce2SnowLay(L,NY,NX)=CumIce2SnowLay(L,NY,NX)+IceXfer2SnoLay(L,NY,NX)
        CumHeat2SnowLay(L,NY,NX)=CumHeat2SnowLay(L,NY,NX)+HeatXfer2SnoLay(L,NY,NX)
      ENDIF  
    ENDIF

    call UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)

    call UpdateSoluteInSnow(L,NY,NX)
  ENDDO D9780

!
!     SNOW RUNOFF from snow laye near soil surface
!
  VLDrySnoWE_col(1,NY,NX)=VLDrySnoWE_col(1,NY,NX)+TDrysnoBySnowRedist(NY,NX)
  VLWatSnow_col(1,NY,NX)=VLWatSnow_col(1,NY,NX)+TWatBySnowRedist(NY,NX)
  VLIceSnow_col(1,NY,NX)=VLIceSnow_col(1,NY,NX)+TIceBySnowRedist(NY,NX)
  ENGYW=VLHeatCapSnow_col(1,NY,NX)*TKSnow(1,NY,NX)
  VLHeatCapSnow_col(1,NY,NX)=cps*VLDrySnoWE_col(1,NY,NX)+cpw*VLWatSnow_col(1,NY,NX)+cpi*VLIceSnow_col(1,NY,NX)

  IF(VLHeatCapSnow_col(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
!    if(I>=138.and.I<=139)print*,I+J/24.,'revi',TKSnow(1,NY,NX),ENGYW,THeatBySnowRedist(NY,NX),VLHeatCapSnow_col(1,NY,NX)    
    TKSnow(1,NY,NX)=(ENGYW+THeatBySnowRedist(NY,NX))/VLHeatCapSnow_col(1,NY,NX)    
  ELSE
    TKSnow(1,NY,NX)=TairK(NY,NX)
  ENDIF
!  if(I>=138.and.I<=139)print*,'tairk',TKSnow(1,NY,NX),TairK(NY,NX)
  if(TKSnow(1,NY,NX)<0._r8)call endrun(trim(mod_filename)//' at line',__LINE__)   
!  write(133,*)I+J/24.,'snow',VcumSnowWE(NY,NX),TKSnow(1,NY,NX),TairK(NY,NX),THeatBySnowRedist(NY,NX)
  VcumDrySnoWE_col(NY,NX)=sum(VLDrySnoWE_col(1:JS,NY,NX))
  VcumWatSnow_col(NY,NX)=sum(VLWatSnow_col(1:JS,NY,NX))
  VcumIceSnow_col(NY,NX)=sum(VLIceSnow_col(1:JS,NY,NX))
  VcumSnoDWI(NY,NX)=sum(VLSnoDWIprev_col(1:JS,NY,NX))
  SnowDepth(NY,NX)=sum(SnowThickL_col(1:JS,NY,NX))
  VcumSnowWE(NY,NX)=VcumDrySnoWE_col(NY,NX)+VcumIceSnow_col(NY,NX)*DENSICE+VcumWatSnow_col(NY,NX) 

  write(*,*) "----- SNOW Parameters -----"
  write(*,*) "Number of Snow Layers: ", JS
  write(*,*) "Snow Depth: ",SnowDepth(NY,NX)
  write(*,*) "TKSnow: ", TKSnow (1,NY,NX)
  write(*,*) "Dry SWE: ", VcumDrySnoWE_col(NY,NX)
  write(*,*) "Ice SWE: ", VcumIceSnow_col(NY,NX)
  write(*,*) "Wat SWE: ", VcumWatSnow_col(NY,NX)
  write(*,*) "---------------------------"
!
! IF SNOWPACK DISAPPEARS

  call SnowpackDisapper(I,J,NY,NX)

  TCSnow(1,NY,NX)=units%Kelvin2Celcius(TKSnow(1,NY,NX))  

  end subroutine SnowMassUpdate

!------------------------------------------------------------------------------------------

  subroutine SnowpackDisapper(I,J,NY,NX)
  implicit none
  integer, intent(in) :: NY,NX,I,J

  integer :: L
  real(r8) :: tksx
  real(r8) :: FLWS,FLWW,FLWI
  real(r8) :: HeatFlo2Surface,ENGYS,ENGY1,ENGY
!     begin_execution
!

  IF(VLHeatCapSnow_col(1,NY,NX).GT.0.0_r8 .AND. VLHeatCapSnow_col(1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) &
    .AND. TairK(NY,NX).GT.TFICE)THEN
    !air temperature above freezing, surface snow layer heat insignificant, so it is merged
    !to the surface layer, and all varaibles are reset
    ENGYS=TKSnow(1,NY,NX)*VLHeatCapSnow_col(1,NY,NX)
    ENGY1=TKS(NUM(NY,NX),NY,NX)*VHeatCapacity(NUM(NY,NX),NY,NX)
    FLWS=VLDrySnoWE_col(1,NY,NX)
    FLWW=VLWatSnow_col(1,NY,NX)
    FLWI=VLIceSnow_col(1,NY,NX)
    HeatFlo2Surface=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TKSnow(1,NY,NX)

    !reset snow layer variables
    VLDrySnoWE_col(1,NY,NX)=0.0_r8
    VLWatSnow_col(1,NY,NX)=0.0_r8
    VLIceSnow_col(1,NY,NX)=0.0_r8
    VLHeatCapSnow_col(1,NY,NX)=0.0_r8
    VcumDrySnoWE_col(NY,NX)=0.0_r8
    VcumWatSnow_col(NY,NX)=0.0_r8
    VcumIceSnow_col(NY,NX)=0.0_r8
    VcumSnoDWI(NY,NX)=0.0_r8
    SnowDepth(NY,NX)=0.0_r8    
    VcumSnowWE(NY,NX)=0._r8
    D9770: DO L=1,JS
      SnoDensL(L,NY,NX)=NewSnowDens(NY,NX)
      if(L/=1)TKSnow(L,NY,NX)=spval
    ENDDO D9770

    !update top soil layer variables
    !maybe should be add to surface residual layer?
    VLWatMicP_vr(NUM(NY,NX),NY,NX)=VLWatMicP_vr(NUM(NY,NX),NY,NX)+FLWW
    VLiceMicP(NUM(NY,NX),NY,NX)=VLiceMicP(NUM(NY,NX),NY,NX)+FLWI+FLWS/DENSICE   

    ENGY=VHeatCapacity(NUM(NY,NX),NY,NX)*TKS(NUM(NY,NX),NY,NX)
    VHeatCapacity(NUM(NY,NX),NY,NX)=VHeatCapacitySoilM(NUM(NY,NX),NY,NX) &
      +cpw*(VLWatMicP_vr(NUM(NY,NX),NY,NX)+VLWatMacP(NUM(NY,NX),NY,NX)) &
      +cpi*(VLiceMicP(NUM(NY,NX),NY,NX)+VLiceMacP(NUM(NY,NX),NY,NX))

    IF(VHeatCapacity(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      TKSX=TKS(NUM(NY,NX),NY,NX)
      TKS(NUM(NY,NX),NY,NX)=(ENGY+HeatFlo2Surface)/VHeatCapacity(NUM(NY,NX),NY,NX)
!      if(abs(TKS(NUM(NY,NX),NY,NX)/tksx-1._r8)>0.025_r8)then
!        TKS(NUM(NY,NX),NY,NX)=TKSX
!      endif
    ELSE
      TKS(NUM(NY,NX),NY,NX)=TairK(NY,NX)
    ENDIF

  ENDIF
  end subroutine SnowpackDisapper
!------------------------------------------------------------------------------------------
  subroutine DealHighTempSnow(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: SnoIce,dHeat,dIce,ENGY,dHPhaseChange
  real(r8) :: dVice,fi,fs,cpold

  SnoIce=VLDrySnoWE_col(L,NY,NX)+VLIceSnow_col(L,NY,NX)*DENSICE
  if(SnoIce<ZEROS(NY,NX))return
!    write(*,*)'not sufficient ice to counter high temperature drysnow ice',I+J/24.,VLDrySnoWE_col(L,NY,NX),VLIceSnow_col(L,NY,NX)
!  endif
  
  !the starting enthalpy
  ENGY=VLHeatCapSnow_col(L,NY,NX)*TKSnow(L,NY,NX)
  dHeat=VLHeatCapSnow_col(L,NY,NX)*(TKSnow(L,NY,NX)-TFICE)
  fs=AZMAX1(VLDrySnoWE_col(L,NY,NX))/SnoIce
  fi=(1._r8-fs)/DENSICE

  dIce=dHeat/(LtHeatIceMelt+(cpw-cpi*fi-cps*fs)*TFICE)        
  
  if(dIce<SnoIce)then
  !some ice melt
  !  print*,VLDrySnoWE_col(L,NY,NX),VLIceSnow_col(L,NY,NX),VLWatSnow_col(L,NY,NX),fs,fi  
    dHPhaseChange=-dIce*LtHeatIceMelt  !thaw, cooling    
    dVice=dIce/SnoIce
    VLDrySnoWE_col(L,NY,NX)=VLDrySnoWE_col(L,NY,NX)-fs*dIce
    VLIceSnow_col(L,NY,NX)=VLIceSnow_col(L,NY,NX)-fi*dIce
    VLWatSnow_col(L,NY,NX)=VLWatSnow_col(L,NY,NX)+dIce
    cpold=VLHeatCapSnow_col(L,NY,NX)
    VLHeatCapSnow_col(L,NY,NX)=cps*VLDrySnoWE_col(L,NY,NX)+cpw*VLWatSnow_col(L,NY,NX)+cpi*VLIceSnow_col(L,NY,NX)
!    print*,cpold,VLHeatCapSnow_col(L,NY,NX),VLDrySnoWE_col(L,NY,NX),VLIceSnow_col(L,NY,NX),VLWatSnow_col(L,NY,NX)
    TKSnow(L,NY,NX)=(ENGY+dHPhaseChange)/VLHeatCapSnow_col(L,NY,NX)
  else
    !all ice melt
    dHPhaseChange=-SnoIce*LtHeatIceMelt  !thaw, cooling
    VLDrySnoWE_col(L,NY,NX)=0._r8
    VLIceSnow_col(L,NY,NX)=0._r8
    VLWatSnow_col(L,NY,NX)=VLWatSnow_col(L,NY,NX)+SnoIce
    VLHeatCapSnow_col(L,NY,NX)=cpw*VLWatSnow_col(L,NY,NX)
  endif
  TKSnow(L,NY,NX)=(ENGY+dHPhaseChange)/VLHeatCapSnow_col(L,NY,NX)  
!  print*,'Hightemp',TKSnow(L,NY,NX)
  end subroutine DealHighTempSnow
!------------------------------------------------------------------------------------------

  subroutine DealNegativeSnowMass(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: dHPhaseChange,VLDrySnoWEtmp
  real(r8) :: SnoWEtot,ENGYW,dVLDrySnoWEtmp
  real(r8) :: tKNew,vlheatnew,dIce,dHeat,fs,fi

  !total snow mass
  SnoWEtot=VLDrySnoWE_col(L,NY,NX)+VLWatSnow_col(L,NY,NX)+VLIceSnow_col(L,NY,NX)*DENSICE  
  if(TKSnow(L,NY,NX)<TFICE)then
!     write(*,*)'negative water',L,TKSnow(L,NY,NX)
!    call endrun('Temeprature too low to resolve negative water '//trim(mod_filename)//' at line',__LINE__)         
  endif
  if(SnoWEtot<0._r8)then
    call endrun('Negative snow mass '//trim(mod_filename)//' at line',__LINE__)         
  endif
  !the starting enthalpy
  ENGYW=VLHeatCapSnow_col(L,NY,NX)*TKSnow(L,NY,NX)  

  !thaw all ice + snow, absorb heat/cooling (<0)
!  if(I>=138.and.I<=139)write(149,*)'neg',I+J/24.,VLWatSnow_col(L,NY,NX),VLIceSnow_col(L,NY,NX),VLDrySnoWE_col(L,NY,NX),&
!    TKSnow(L,NY,NX),(CumHeat2SnowLay(L,NY,NX)+XPhaseChangeHeatL(L,NY,NX)+ENGYW)/(cpi*TFICE)
  dHPhaseChange=-LtHeatIceMelt*(VLDrySnoWE_col(L,NY,NX)+VLIceSnow_col(L,NY,NX)*DENSICE)
  
  vlheatnew=cpw*SnoWEtot
  !compute potential temperature
  tkNew=(ENGYW+CumHeat2SnowLay(L,NY,NX)+XPhaseChangeHeatL(L,NY,NX)+dHPhaseChange)/vlheatnew

  if(tkNew>TFICE)then
    !all ice & snow are melt
    XPhaseChangeHeatL(L,NY,NX)=XPhaseChangeHeatL(L,NY,NX)+dHPhaseChange
    VLWatSnow_col(L,NY,NX)=SnoWEtot
    VLDrySnoWE_col(L,NY,NX)=0._r8
    VLIceSnow_col(L,NY,NX)=0._r8
  else
    !some snow and ice remains, refreeze/release heat
    dHeat=-(tkNew-TFICE)*vlheatnew    !>0, enthalpy excess for freeze thaw
    VLDrySnoWEtmp=AZMAX1(VLDrySnoWE_col(L,NY,NX))+AZMAX1(VLIceSnow_col(L,NY,NX))*DENSICE    
    fs=AZMAX1(VLDrySnoWE_col(L,NY,NX))/VLDrySnoWEtmp
    fi=(1._r8-fs)/DENSICE
    dIce=dHeat/(LtHeatIceMelt+(cpw-cpi*fi-cps*fs)*TFICE)          !>0

    if(dIce <= SnoWEtot)then
      if(VLDrySnoWEtmp>0._r8)then      
        dVLDrySnoWEtmp=dIce/VLDrySnoWEtmp
        VLDrySnoWE_col(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLDrySnoWE_col(L,NY,NX))
        VLIceSnow_col(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLIceSnow_col(L,NY,NX))
        VLWatSnow_col(L,NY,NX)=AZMAX1(SnoWEtot-VLDrySnoWE_col(L,NY,NX)-VLIceSnow_col(L,NY,NX)*DENSICE)
        !freeze releases heat
        dHPhaseChange=(VLDrySnoWE_col(L,NY,NX)-VLIceSnow_col(L,NY,NX)*DENSICE-VLDrySnoWEtmp)*LtHeatIceMelt
      else
        !meaning water is positive, freezing occured, heat release
        VLIceSnow_col(L,NY,NX)=dIce/DENSICE
        VLDrySnoWE_col(L,NY,NX)=0._r8
        VLWatSnow_col(L,NY,NX)=AZMAX1(SnoWEtot-VLDrySnoWE_col(L,NY,NX)-VLIceSnow_col(L,NY,NX)*DENSICE)
        dHPhaseChange=dIce*LtHeatIceMelt
      endif

      XPhaseChangeHeatL(L,NY,NX)=XPhaseChangeHeatL(L,NY,NX)+dHPhaseChange
    else
      !negative water is made by ice, no phase change, but there is conversion
      dVLDrySnoWEtmp=SnoWEtot/VLDrySnoWEtmp
      VLDrySnoWE_col(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLDrySnoWE_col(L,NY,NX))
      VLIceSnow_col(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLIceSnow_col(L,NY,NX))  
      VLWatSnow_col(L,NY,NX)=0._r8    
    endif
  endif
  
  end subroutine DealNegativeSnowMass  

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)
  implicit none
  integer, intent(in) :: L,NY,NX,I,J
  real(r8), intent(inout) :: VOLSWI
  real(r8) :: TKWX,VHCPWZ(JZ,JY,JX)
  real(r8) :: ENGYW
  real(r8) :: VOLSF,TCASF
  real(r8) :: DENSF,DDENS2,DDENS1
  real(r8) :: CVISC,DENSX
  real(r8) :: frcnew,SnowIceMass,VLwatNet
  real(r8) :: vwat,vice,vdry
! begin_execution

  !the line below is a hack, and likely a better snow layering scheme is needed.
  if(L.eq.1.and.isclose(TCSnow(1,NY,NX),spval))then
    TCSnow(1,NY,NX)=units%Kelvin2Celcius(TairK(NY,NX))    
  endif
!
! ADD CHANGES IN SNOW, WATER AND ICE
!
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! CumSno2SnowLay,CumWat2SnowLay,CumIce2SnowLay=net fluxes of snow,water,ice in snowpack
! XSnowThawMassL,XIceThawMassL=freeze-thaw flux from watsub.f
!
!
  vwat=VLWatSnow_col(L,NY,NX);vice=VLIceSnow_col(L,NY,NX);vdry=VLDrySnoWE_col(L,NY,NX)
  VLDrySnoWE_col(L,NY,NX)=VLDrySnoWE_col(L,NY,NX)+CumSno2SnowLay(L,NY,NX)-XSnowThawMassL(L,NY,NX)
  VLWatSnow_col(L,NY,NX)=VLWatSnow_col(L,NY,NX)+CumWat2SnowLay(L,NY,NX)+XSnowThawMassL(L,NY,NX)+XIceThawMassL(L,NY,NX)
  VLIceSnow_col(L,NY,NX)=VLIceSnow_col(L,NY,NX)+CumIce2SnowLay(L,NY,NX)-XIceThawMassL(L,NY,NX)/DENSICE

  if(any((/VLDrySnoWE_col(L,NY,NX),VLWatSnow_col(L,NY,NX),VLIceSnow_col(L,NY,NX)/)<0._r8))then    
    call DealNegativeSnowMass(I,J,L,NY,NX)
  endif
!
! ACCUMULATE SNOW MASS FOR CALCULATING COMPRESSION
!
! VOLWSI=accumulated water equivalent volume in snowpack
! XFLWS=snow transfer from watsub.f
! VOLSF=snowfall volume
! DENSS=snow density in layer
!
  
  IF(L.EQ.1)THEN
    VOLSWI=VOLSWI+(VLDrySnoWE_col(L,NY,NX)+VLWatSnow_col(L,NY,NX)+VLIceSnow_col(L,NY,NX)*DENSICE)
    if(VOLSWI<0._r8)then
      write(*,*)'VOLSWI=',VOLSWI,VLDrySnoWE_col(L,NY,NX),VLWatSnow_col(L,NY,NX),VLIceSnow_col(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
!
!   RESET SNOW SURFACE DENSITY FOR SNOWFALL
!
    IF(SnoXfer2SnoLay(L,NY,NX).GT.0.0_r8)THEN
      DENSX=SnoDensL(L,NY,NX)
      TCASF=AMAX1(-15.0_r8,AMIN1(2.0_r8,TCA(NY,NX)))
      !fresh snow density
      DENSF=0.05_r8+1.7E-03_r8*(TCASF+15.0_r8)**1.5_r8
      VOLSF=AMIN1(SnoXfer2SnoLay(L,NY,NX),VLDrySnoWE_col(L,NY,NX))/DENSF + &
        AZMAX1(VLDrySnoWE_col(L,NY,NX)-SnoXfer2SnoLay(L,NY,NX))/SnoDensL(L,NY,NX)        
      if(VOLSF>0._r8)SnoDensL(L,NY,NX)=VLDrySnoWE_col(L,NY,NX)/VOLSF
      !write(*,*)'xVOLSSL=',VLDrySnoWE_col(L,NY,NX),SnoXfer2SnoLay(L,NY,NX),SnoDensL(L,NY,NX),VOLSF
    ENDIF
  ELSE
    VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE_col(L-1,NY,NX)+VLWatSnow_col(L-1,NY,NX) &
      +VLIceSnow_col(L-1,NY,NX)*DENSICE+VLDrySnoWE_col(L,NY,NX)+VLWatSnow_col(L,NY,NX) &
      +VLIceSnow_col(L,NY,NX)*DENSICE)
      
    if(VOLSWI<0._r8)then
      write(*,*)'iVOLSWI=',VOLSWI,VLDrySnoWE_col(L-1,NY,NX)+VLWatSnow_col(L-1,NY,NX) &
        +VLIceSnow_col(L-1,NY,NX)*DENSICE,VLDrySnoWE_col(L,NY,NX)+VLWatSnow_col(L,NY,NX) &
        +VLIceSnow_col(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDIF
!
! SNOWPACK COMPRESSION
!
! DDENS1,DDENS2=temperature, compression effect on snow density
! DENSS=snow density in layer
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! VLSnoDWIprev_col=snowpack layer volume
! DLYRS=snowpack layer depth
! cumSnowDepz_col=cumulative depth to bottom of snowpack layer
! VHCPW=snowpack layer heat capacity
! TKW,TCSnow=snowpack layer temperature K,oC
! CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
! XPhaseChangeHeatL=latent heat flux from freeze-thaw from watsub.f
! HEATIN=cumulative net surface heat transfer
! VOLSS,VOLWS,VOLIS=total snow water equivalent, water, ice content of snowpack
! VOLS,SnowDepth=total snowpack volume, depth
!
!  if(I>=138.and.I<=139)print*,L,SnoDensL(L,NY,NX),TKSnow(L,NY,NX)
  if(TKSnow(L,NY,NX)<0.)call endrun(trim(mod_filename)//' at line',__LINE__)   
  if(SnoDensL(L,NY,NX)>0._r8)then
    IF(SnoDensL(L,NY,NX).LT.0.25_r8)THEN
      DDENS1=SnoDensL(L,NY,NX)*1.0E-05_r8*EXP(0.04_r8*TCSnow(L,NY,NX))
    ELSE
      DDENS1=0.0_r8
    ENDIF
    CVISC=0.25_r8*EXP(-0.08_r8*TCSnow(L,NY,NX)+23.0_r8*SnoDensL(L,NY,NX))

    DDENS2=SnoDensL(L,NY,NX)*VOLSWI/(AREA(3,NU(NY,NX),NY,NX)*CVISC)
    
    SnoDensL(L,NY,NX)=SnoDensL(L,NY,NX)+DDENS1+DDENS2
    if(SnoDensL(L,NY,NX)<0._r8)then
      write(*,*)'DDENS1=',SnoDensL(L,NY,NX),DDENS1,DDENS2,L
      write(*,*)SnoXfer2SnoLay(L,NY,NX),VLDrySnoWE_col(L,NY,NX)
      call endrun("negative snow dens")
    endif  
  endif
  
  !there is snow in layer L
  VLSnoDWIprev_col(L,NY,NX)=VLDrySnoWE_col(L,NY,NX)/SnoDensL(L,NY,NX)+VLWatSnow_col(L,NY,NX)+VLIceSnow_col(L,NY,NX)

  IF(VLSnoDWIprev_col(L,NY,NX) > ZEROS2(NY,NX))THEN
    SnowThickL_col(L,NY,NX)=AZMAX1(VLSnoDWIprev_col(L,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
    cumSnowDepz_col(L,NY,NX)=cumSnowDepz_col(L-1,NY,NX)+SnowThickL_col(L,NY,NX)
    VHCPWZ(L,NY,NX)=VLHeatCapSnow_col(L,NY,NX)
    TKWX=TKSnow(L,NY,NX)
    ENGYW=VLHeatCapSnow_col(L,NY,NX)*TKWX
    VLHeatCapSnow_col(L,NY,NX)=cps*VLDrySnoWE_col(L,NY,NX)+cpw*VLWatSnow_col(L,NY,NX)+cpi*VLIceSnow_col(L,NY,NX)    
    IF(VHCPWZ(L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX) .AND. VLHeatCapSnow_col(L,NY,NX).GT.ZEROS(NY,NX))THEN
      !there is significant snow layer mass
!      if(I>=138.and.I<=139)write(149,*)'rest',L,ENGYW,CumHeat2SnowLay(L,NY,NX),XPhaseChangeHeatL(L,NY,NX),VLHeatCapSnow_col(L,NY,NX)
      TKSnow(L,NY,NX)=(ENGYW+CumHeat2SnowLay(L,NY,NX)+XPhaseChangeHeatL(L,NY,NX))/VLHeatCapSnow_col(L,NY,NX)

      if(TKSnow(L,NY,NX)>280._r8)call DealHighTempSnow(I,J,L,NY,NX)
    ELSE
      !there is no significant snow mass      
      IF(L.EQ.1)THEN
        !if current layer is top layer
        TKSnow(L,NY,NX)=TairK(NY,NX)
      ELSE
        !if it is not the top layer
        TKSnow(L,NY,NX)=TKSnow(L-1,NY,NX)
      ENDIF
      IF(VHCPWZ(L,NY,NX).GT.ZEROS(NY,NX))THEN
        HEATIN=HEATIN+(TKSnow(L,NY,NX)-TKWX)*VLHeatCapSnow_col(L,NY,NX)
      ENDIF
    ENDIF
  ELSE
    !there is no snow in layer L
    VLDrySnoWE_col(L,NY,NX)=0.0_r8
    VLWatSnow_col(L,NY,NX)=0.0_r8
    VLIceSnow_col(L,NY,NX)=0.0_r8
    VLSnoDWIprev_col(L,NY,NX)=0.0_r8
    SnowThickL_col(L,NY,NX)=0.0_r8
    cumSnowDepz_col(L,NY,NX)=cumSnowDepz_col(L-1,NY,NX)
    IF(L.EQ.1)THEN
      TKSnow(L,NY,NX)=TairK(NY,NX)
    ELSE
      TKSnow(L,NY,NX)=TKSnow(L-1,NY,NX)
    ENDIF
  ENDIF  
  TCSnow(L,NY,NX)=units%Kelvin2Celcius(TKSnow(L,NY,NX))
!  if(I>=138.and.I<=139)print*,'updtesnowl',L,TKSnow(L,NY,NX)
  end subroutine UpdateSnowLayerL
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSnow(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: NTG,NTN,NTSA
  !     begin_execution
  !     SNOWPACK SOLUTE CONTENT
  !
  !     *W2=solute content of snowpack
  !     T*BLS=net solute flux in snowpack
  !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
  !             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
  !
  DO NTG=idg_beg,idg_end-1
    trcg_solsml(NTG,L,NY,NX)=trcg_solsml(NTG,L,NY,NX)+trcg_TBLS(NTG,L,NY,NX)
  ENDDO

  DO NTN =ids_nut_beg,ids_nuts_end
    trcn_solsml(NTN,L,NY,NX)=trcn_solsml(NTN,L,NY,NX)+trcn_TBLS(NTN,L,NY,NX)
  ENDDO
  !
  !     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
  !          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
  !          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
  !          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
  !          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
  !          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
  !          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
  !     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
  !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
  !
  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      trcs_solsml(NTSA,L,NY,NX)=trcs_solsml(NTSA,L,NY,NX)+trcSalt_TBLS(NTSA,L,NY,NX)
    ENDDO

  ENDIF

  end subroutine UpdateSoluteInSnow

!------------------------------------------------------------------------------------------

  subroutine SnowpackLayering(I,J,NY,NX)
  !
  !Description:
  !Relayering snow after mass update 
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: FX,FY
  integer :: L,L1,L0
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
  real(r8) :: DDLYXS,DDLYRS
  real(r8) :: DDLYXX,VOLSLX
  integer :: IFLGLS,NTN,NTG,NTU,NTSA
  integer, parameter :: inochange=0
  integer, parameter :: iexpand=1
  integer, parameter :: ishrink=2
!     begin_execution
! from surface to bottom, and modify the bottom layer
! there is snow

!  if(I>=138.and.I<=139)print*,I+J/24.,'bflay',TKSnow(1:JS,NY,NX)
  IF(VLHeatCapSnow_col(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
    D325: DO L=1,JS-1
!      VOLSLX=VLSnoDWIprev_col(L,NY,NX)
      IF(VLSnoDWIprev_col(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        !compute the free volume/thickness: DDLYXS
        DDLYXS=(VLSnoDWIMax_col(L,NY,NX)-VLDrySnoWE_col(L,NY,NX)/SnoDensL(L,NY,NX) &
          -VLWatSnow_col(L,NY,NX)-VLIceSnow_col(L,NY,NX))/AREA(3,L,NY,NX)
!        DDLYXX=DDLYXS
        IF(DDLYXS.LT.-ZERO .OR. SnowThickL_col(L+1,NY,NX).GT.ZERO)THEN
          !current volume is greater than allowed, or next layer exists   
          !case 1: DDLYXS< 0, layer L extends into layer L+1, DDLYRS<0: amount of expand layer L into layer L+1
          !case 2: 0<DDLYXS<SnowThickL_col(L+1,NY,NX), layer L still has space, can take DDLYRS of layer L+1       
          !case 3: DDLYXS>SnowThickL_col(L+1,NY,NX), layer L still has space, and can even hold DDLYRS of layer L+1
          !
          DDLYRS=AMIN1(DDLYXS,SnowThickL_col(L+1,NY,NX))
          IFLGLS=iexpand         !expand
        ELSE
          !layer L has space and layer L+1 disappears, so layer L+1 combines into layer L
          !volume less than allowed, and no next layer
          !DDLYXS: is the depth change of layer L
          DDLYXS=(VLSnoDWIprev_col(L,NY,NX)-VLDrySnoWE_col(L,NY,NX)/SnoDensL(L,NY,NX) &
            -VLWatSnow_col(L,NY,NX)-VLIceSnow_col(L,NY,NX))/AREA(3,L,NY,NX)
          DDLYRS=DDLYXS
          IFLGLS=ishrink         !shrink
        ENDIF
      ELSE
        !current layer is empty, do nothing
        DDLYRS=0.0_r8      !no change
        IFLGLS=inochange
      ENDIF
      !
      !     RESET SNOW LAYER DEPTHS
      !
      cumSnowDepz_col(L,NY,NX)=cumSnowDepz_col(L,NY,NX)+DDLYRS
      SnowThickL_col(L,NY,NX)=cumSnowDepz_col(L,NY,NX)-cumSnowDepz_col(L-1,NY,NX)
!
      !     TRANSFER STATE VARIABLES BETWEEN LAYERS
      !
      IF(ABS(DDLYRS).GT.ZERO)THEN
        IF(DDLYRS.GT.0.0_r8)THEN
          !incoporating L+1 into L
          L1=L
          L0=L+1
          IF(DDLYRS.LT.DDLYXS)THEN
            !full L0 into L1
            FX=1.0_r8
          ELSE
            !partial L0 into L1
            FX=AMIN1(1.0_r8,DDLYRS*AREA(3,L0,NY,NX)/VLSnoDWIprev_col(L0,NY,NX))
          ENDIF
        ELSE
          !expanding L into L+1
          L1=L+1
          L0=L
          IF(VLSnoDWIprev_col(L0,NY,NX).LT.VLSnoDWIMax_col(L0,NY,NX))THEN
            FX=0.0_r8
          ELSE
            !FX fraction to be donated from L0 to L1
            FX=AMIN1(1.0_r8,-DDLYRS*AREA(3,L0,NY,NX)/VLSnoDWIprev_col(L0,NY,NX))
          ENDIF
        ENDIF
!   donor L0, target L1
        IF(FX.GT.0.0_r8)THEN
          FY=1.0_r8-FX
!          print*,'relay',L0,L1,FX,FY,DDLYRS,DDLYXS,-ZERO,SnowThickL_col(L+1,NY,NX),VLSnoDWIprev_col(L0,NY,NX),VLSnoDWIMax_col(L0,NY,NX)
!
!     TARGET SNOW LAYER
!         volume/mass
          VLDrySnoWE_col(L1,NY,NX)=VLDrySnoWE_col(L1,NY,NX)+FX*VLDrySnoWE_col(L0,NY,NX)
          VLWatSnow_col(L1,NY,NX)=VLWatSnow_col(L1,NY,NX)+FX*VLWatSnow_col(L0,NY,NX)
          VLIceSnow_col(L1,NY,NX)=VLIceSnow_col(L1,NY,NX)+FX*VLIceSnow_col(L0,NY,NX)          
          VLSnoDWIprev_col(L1,NY,NX)=VLDrySnoWE_col(L1,NY,NX)/SnoDensL(L1,NY,NX)+VLWatSnow_col(L1,NY,NX)+VLIceSnow_col(L1,NY,NX)
!         energy
          ENGY1X=VLHeatCapSnow_col(L1,NY,NX)*TKSnow(L1,NY,NX)
          ENGY0X=VLHeatCapSnow_col(L0,NY,NX)*TKSnow(L0,NY,NX)
          ENGY1=ENGY1X+FX*ENGY0X
          VLHeatCapSnow_col(L1,NY,NX)=cps*VLDrySnoWE_col(L1,NY,NX)+cpw*VLWatSnow_col(L1,NY,NX)+cpi*VLIceSnow_col(L1,NY,NX)

          IF(VLHeatCapSnow_col(L1,NY,NX).GT.ZEROS(NY,NX))THEN
            TKSnow(L1,NY,NX)=ENGY1/VLHeatCapSnow_col(L1,NY,NX)
          ELSE
            TKSnow(L1,NY,NX)=TKSnow(L0,NY,NX)
          ENDIF
          TCSnow(L1,NY,NX)=units%Kelvin2Celcius(TKSnow(L1,NY,NX))
!          chemicals
          !gas
          DO NTG=idg_beg,idg_end-1
            trcg_solsml(NTG,L1,NY,NX)=trcg_solsml(NTG,L1,NY,NX)+FX*trcg_solsml(NTG,L0,NY,NX)
          ENDDO

          !nutrients
          DO NTN=ids_nut_beg,ids_nuts_end
            trcn_solsml(NTN,L1,NY,NX)=trcn_solsml(NTN,L1,NY,NX)+FX*trcn_solsml(NTN,L0,NY,NX)
          ENDDO
          !salt
          IF(salt_model)THEN
            DO NTSA=idsalt_beg,idsalt_end
              trcs_solsml(NTSA,L1,NY,NX)=trcs_solsml(NTSA,L1,NY,NX)+FX*trcs_solsml(NTSA,L0,NY,NX)
            ENDDO
          ENDIF
!
!     SOURCE SNOW LAYER
!         volume/mass
          VLDrySnoWE_col(L0,NY,NX)=FY*VLDrySnoWE_col(L0,NY,NX)
          VLWatSnow_col(L0,NY,NX)=FY*VLWatSnow_col(L0,NY,NX)
          VLIceSnow_col(L0,NY,NX)=FY*VLIceSnow_col(L0,NY,NX)
          VLSnoDWIprev_col(L0,NY,NX)=VLDrySnoWE_col(L0,NY,NX)/SnoDensL(L0,NY,NX)+VLWatSnow_col(L0,NY,NX)+VLIceSnow_col(L0,NY,NX)
!     energy 
          ENGY0=FY*ENGY0X
          VLHeatCapSnow_col(L0,NY,NX)=cps*VLDrySnoWE_col(L0,NY,NX)+cpw*VLWatSnow_col(L0,NY,NX)+cpi*VLIceSnow_col(L0,NY,NX)
          IF(VLHeatCapSnow_col(L0,NY,NX).GT.ZEROS(NY,NX))THEN
            TKSnow(L0,NY,NX)=ENGY0/VLHeatCapSnow_col(L0,NY,NX)
          ELSE
            TKSnow(L0,NY,NX)=TKSnow(L1,NY,NX)
          ENDIF
          TCSnow(L0,NY,NX)=units%Kelvin2Celcius(TKSnow(L0,NY,NX))
!     chemicals
          DO NTG=idg_beg,idg_end-1
            trcg_solsml(NTG,L0,NY,NX)=FY*trcg_solsml(NTG,L0,NY,NX)
          ENDDO
          DO NTU=ids_nut_beg,ids_nuts_end
            trcn_solsml(NTU,L0,NY,NX)=FY*trcn_solsml(NTU,L0,NY,NX)
          ENDDO

          IF(salt_model)THEN
            DO NTSA=idsalt_beg,idsalt_end
              trcs_solsml(NTSA,L0,NY,NX)=FY*trcs_solsml(NTSA,L0,NY,NX)
            ENDDO

          ENDIF
        ENDIF
      ENDIF
    ENDDO D325
  ENDIF
  !update volumetric snow heat capacity
  nsnol_col(NY,NX)=0
  DO L=1,JS
!     if(I>=138.and.I<140)print*,I+J/24.,'layering',L,VLDrySnoWE_col(L,NY,NX),VLWatSnow_col(L,NY,NX),VLIceSnow_col(L,NY,NX),TKSnow(L,NY,NX)
    if(TKSnow(L,NY,NX)/=spval)then
      VLHeatCapSnow_col(L,NY,NX)=cps*VLDrySnoWE_col(L,NY,NX)+cpw*VLWatSnow_col(L,NY,NX)+cpi*VLIceSnow_col(L,NY,NX)
      if(VLHeatCapSnow_col(L,NY,NX)>ZEROS(NY,NX))then
        nsnol_col(NY,NX)=nsnol_col(NY,NX)+1
      endif
    else
      VLDrySnoWE_col(L,NY,NX)=0._r8
      VLWatSnow_col(L,NY,NX)=0._r8
      VLIceSnow_col(L,NY,NX)=0._r8
      VLHeatCapSnow_col(L,NY,NX)=0._r8
      SnoDensL(L,NY,NX)=NewSnowDens(NY,NX)
    endif
    if(VLHeatCapSnow_col(L,NY,NX)<ZEROS(NY,NX))then
      if(L>1)then
        TKSnow(L,NY,NX)=spval        
      else
        TKSnow(L,NY,NX)=TairK(NY,NX) 
      endif  
    endif  
    !move up to handel loss of surface layer, this should rarely occur, but round off error may trigger it
    if(L>1 .and. VLHeatCapSnow_col(L,NY,NX)>=ZEROS(NY,NX))then
      if (TKSnow(L-1,NY,NX)==spval)then
        VLDrySnoWE_col(L-1,NY,NX)=VLDrySnoWE_col(L,NY,NX)
        VLWatSnow_col(L-1,NY,NX)=VLWatSnow_col(L,NY,NX)
        VLIceSnow_col(L-1,NY,NX)=VLIceSnow_col(L,NY,NX)
        VLHeatCapSnow_col(L-1,NY,NX)=VLHeatCapSnow_col(L,NY,NX)
        TKSnow(L-1,NY,NX)=TKSnow(L,NY,NX)
        SnoDensL(L-1,NY,NX)=SnoDensL(L,NY,NX)

        VLDrySnoWE_col(L,NY,NX)=0._r8
        VLWatSnow_col(L,NY,NX)=0._r8
        VLIceSnow_col(L,NY,NX)=0._r8
        VLHeatCapSnow_col(L,NY,NX)=0._r8
        TKSnow(L,NY,NX)=spval
        SnoDensL(L,NY,NX)=NewSnowDens(NY,NX)
      endif
    endif 
  ENDDO  
  
  end subroutine SnowpackLayering

!------------------------------------------------------------------------------------------
  subroutine ZeroSnowArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  TDrysnoBySnowRedist(NY,NX)=0.0_r8
  TWatBySnowRedist(NY,NX)=0.0_r8
  TIceBySnowRedist(NY,NX)=0.0_r8
  THeatBySnowRedist(NY,NX)=0.0_r8

  trcn_TFloXSurRunoff(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8
  trcg_QSS(idg_beg:idg_end-1,NY,NX)=0.0_r8  
  trcn_QSS(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8
  trcg_TFloXSurRunoff(idg_beg:idg_end-1,NY,NX)=0.0_r8

  DO  L=1,JS
    trcg_TBLS(idg_beg:idg_end-1,L,NY,NX)=0.0_r8
    trcn_TBLS(ids_nut_beg:ids_nuts_end,L,NY,NX)=0.0_r8

  ENDDO

  IF(salt_model)THEN
!     INITIALIZE NET SOLUTE AND GAS FLUXES FROM SNOWPACK DRIFT
!
    trcSalt_TQR(idsalt_beg:idsalt_end,NY,NX)=0.0_r8
    trcSalt_TQS(idsalt_beg:idsalt_end,NY,NX)=0.0_r8
    DO  L=1,JS
      trcSalt_TBLS(idsalt_beg:idsalt_end,L,NY,NX)=0.0_r8
    ENDDO
  endif
  end subroutine ZeroSnowArrays

end module SnowBalanceMod
