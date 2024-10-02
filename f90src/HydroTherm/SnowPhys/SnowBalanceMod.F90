module SnowBalanceMod
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: spval => DAT_CONST_SPVAL
  use abortutils,     only: endrun
  use EcoSIMCtrlMod,  only: lverb,snowRedist_model
  use minimathmod,    only: AZMAX1, isclose, AZMIN1,AZMAX1d  
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
  public :: DebugSnowPrint
  contains

  subroutine DebugSnowPrint(I,J,NY,NX,info)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), intent(in) :: info
  integer :: L
    if(I==152 .and. J==1 .and. NX==7)then
    L=1
    write(113,*)trim(info)//'snoadd',L,CumSno2SnowL_snvr(L,NY,NX),CumWat2SnowL_snvr(L,NY,NX),CumIce2SnowL_snvr(L,NY,NX),CumHeat2SnowL_snvr(L,NY,NX)
    endif

  end subroutine DebugSnowPrint
!------------------------------------------------------------------------------------------

  subroutine SnowMassUpdate(I,J,NY,NX)
  !
  !Description: update snow mass
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: VOLSWI,ENGYW,TKSnow_pre
  integer :: L
      !
      !     CALCULATE SNOWPACK TEMPERATURE FROM ITS CHANGE
      !     IN HEAT STORAGE
      !
  if(lverb)write(*,*)'SnowMassUpdate'
  !update snow layer from top to bottom  
  VOLSWI     = 0.0_r8
  TKSnow_pre = TKSnow_snvr(1,NY,NX)
!  if(TKSnow_snvr(1,NY,NX)<230._r8)then
!    write(*,*)I,J,TKSnow_pre,TairK_col(NY,NX)
!    call endrun(trim(mod_filename)//' at line',__LINE__)     
!  endif  


  D9780: DO L=1,JS
    IF(VLHeatCapSnow_snvr(L,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) .and. L.EQ.1)THEN
      IF(abs(SnoXfer2SnoLay_snvr(L,NY,NX))>0._r8)THEN
        CumSno2SnowL_snvr(L,NY,NX)  = CumSno2SnowL_snvr(L,NY,NX)+SnoXfer2SnoLay_snvr(L,NY,NX)
        CumWat2SnowL_snvr(L,NY,NX)  = CumWat2SnowL_snvr(L,NY,NX)+WatXfer2SnoLay_snvr(L,NY,NX)
        CumIce2SnowL_snvr(L,NY,NX)  = CumIce2SnowL_snvr(L,NY,NX)+IceXfer2SnoLay_snvr(L,NY,NX)
        CumHeat2SnowL_snvr(L,NY,NX) = CumHeat2SnowL_snvr(L,NY,NX)+HeatXfer2SnoLay_snvr(L,NY,NX)
      ENDIF  
    ENDIF
    
    call UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)

    call UpdateSoluteInSnow(L,NY,NX)
    
  ENDDO D9780

!  if(TKSnow_snvr(1,NY,NX)<230._r8)then
!    write(*,*)I,J,TKSnow_pre,TairK_col(NY,NX)
!    call endrun(trim(mod_filename)//' at line',__LINE__)     
!  endif  

  VcumDrySnoWE_col(NY,NX) = sum(VLDrySnoWE_snvr(1:JS,NY,NX))
  VcumWatSnow_col(NY,NX)  = sum(VLWatSnow_snvr(1:JS,NY,NX))
  VcumIceSnow_col(NY,NX)  = sum(VLIceSnow_snvr(1:JS,NY,NX))
  VcumSnoDWI_col(NY,NX)   = sum(VLSnoDWIprev_snvr(1:JS,NY,NX))
  SnowDepth_col(NY,NX)    = sum(SnowThickL_snvr(1:JS,NY,NX))
  VcumSnowWE_col(NY,NX)   = VcumDrySnoWE_col(NY,NX)+VcumIceSnow_col(NY,NX)*DENSICE+VcumWatSnow_col(NY,NX)
!
! IF SNOWPACK DISAPPEARS

  call SnowpackDisapper(I,J,NY,NX)

  TCSnow_snvr(1,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(1,NY,NX))  

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

  IF(VLHeatCapSnow_snvr(1,NY,NX).GT.0.0_r8 .AND. VLHeatCapSnow_snvr(1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) &
    .AND. TairK_col(NY,NX).GT.TFICE)THEN

    !air temperature above freezing, surface snow layer heat insignificant, so it is merged
    !to the surface layer, and all varaibles are reset
    ENGYS           = TKSnow_snvr(1,NY,NX)*VLHeatCapSnow_snvr(1,NY,NX)
    ENGY1           = TKS_vr(NUM(NY,NX),NY,NX)*VHeatCapacity_vr(NUM(NY,NX),NY,NX)
    FLWS            = VLDrySnoWE_snvr(1,NY,NX)
    FLWW            = VLWatSnow_snvr(1,NY,NX)
    FLWI            = VLIceSnow_snvr(1,NY,NX)
    HeatFlo2Surface = (cpw*FLWW+cps*FLWS+cpi*FLWI)*TKSnow_snvr(1,NY,NX)

    !reset snow layer variables
    VLDrySnoWE_snvr(1,NY,NX)    = 0.0_r8
    VLWatSnow_snvr(1,NY,NX)     = 0.0_r8
    VLIceSnow_snvr(1,NY,NX)     = 0.0_r8
    VLHeatCapSnow_snvr(1,NY,NX) = 0.0_r8
    VcumDrySnoWE_col(NY,NX)     = 0.0_r8
    VcumWatSnow_col(NY,NX)      = 0.0_r8
    VcumIceSnow_col(NY,NX)      = 0.0_r8
    VcumSnoDWI_col(NY,NX)       = 0.0_r8
    SnowDepth_col(NY,NX)        = 0.0_r8
    VcumSnowWE_col(NY,NX)       = 0.0_r8

    D9770: DO L=1,JS
      SnoDens_snvr(L,NY,NX)=NewSnowDens_col(NY,NX)
      if(L/=1)TKSnow_snvr(L,NY,NX)=spval
    ENDDO D9770

    !update top soil layer variables
    !maybe should be add to surface residual layer?
    VLWatMicP_vr(NUM(NY,NX),NY,NX) = VLWatMicP_vr(NUM(NY,NX),NY,NX)+FLWW
    VLiceMicP_vr(NUM(NY,NX),NY,NX) = VLiceMicP_vr(NUM(NY,NX),NY,NX)+FLWI+FLWS/DENSICE

    ENGY                               = VHeatCapacity_vr(NUM(NY,NX),NY,NX)*TKS_vr(NUM(NY,NX),NY,NX)
    VHeatCapacity_vr(NUM(NY,NX),NY,NX) = VHeatCapacitySoilM_vr(NUM(NY,NX),NY,NX) &
      +cpw*(VLWatMicP_vr(NUM(NY,NX),NY,NX)+VLWatMacP_vr(NUM(NY,NX),NY,NX)) &
      +cpi*(VLiceMicP_vr(NUM(NY,NX),NY,NX)+VLiceMacP_vr(NUM(NY,NX),NY,NX))

    IF(VHeatCapacity_vr(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      TKSX                     = TKS_vr(NUM(NY,NX),NY,NX)
      TKS_vr(NUM(NY,NX),NY,NX) = (ENGY+HeatFlo2Surface)/VHeatCapacity_vr(NUM(NY,NX),NY,NX)
!      if(abs(TKS_vr(NUM(NY,NX),NY,NX)/tksx-1._r8)>0.025_r8)then
!        TKS_vr(NUM(NY,NX),NY,NX)=TKSX
!      endif
    ELSE
      TKS_vr(NUM(NY,NX),NY,NX)=TairK_col(NY,NX)
    ENDIF

  ENDIF
  end subroutine SnowpackDisapper
!------------------------------------------------------------------------------------------
  subroutine DealHighTempSnow(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: SnoIce,dHeat,dIce,ENGY,dHPhaseChange
  real(r8) :: dVice,fi,fs,cpold

  SnoIce=VLDrySnoWE_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  if(SnoIce<ZEROS(NY,NX))return
!    write(*,*)'not sufficient ice to counter high temperature drysnow ice',I+J/24.,VLDrySnoWE_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)
!  endif
  
  !the starting enthalpy
  ENGY=VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)
  dHeat=VLHeatCapSnow_snvr(L,NY,NX)*(TKSnow_snvr(L,NY,NX)-TFICE)
  fs=AZMAX1(VLDrySnoWE_snvr(L,NY,NX))/SnoIce
  fi=(1._r8-fs)/DENSICE

  dIce=dHeat/(LtHeatIceMelt+(cpw-cpi*fi-cps*fs)*TFICE)        
  
  if(dIce<SnoIce)then
  !some ice melt
  !  print*,VLDrySnoWE_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),fs,fi  
    dHPhaseChange=-dIce*LtHeatIceMelt  !thaw, cooling    
    dVice=dIce/SnoIce
    VLDrySnoWE_snvr(L,NY,NX)=VLDrySnoWE_snvr(L,NY,NX)-fs*dIce
    VLIceSnow_snvr(L,NY,NX)=VLIceSnow_snvr(L,NY,NX)-fi*dIce
    VLWatSnow_snvr(L,NY,NX)=VLWatSnow_snvr(L,NY,NX)+dIce
    cpold=VLHeatCapSnow_snvr(L,NY,NX)
    VLHeatCapSnow_snvr(L,NY,NX)=cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
!    print*,cpold,VLHeatCapSnow_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX)
    TKSnow_snvr(L,NY,NX)=(ENGY+dHPhaseChange)/VLHeatCapSnow_snvr(L,NY,NX)
  else
    !all ice melt
    dHPhaseChange=-SnoIce*LtHeatIceMelt  !thaw, cooling
    VLDrySnoWE_snvr(L,NY,NX)=0._r8
    VLIceSnow_snvr(L,NY,NX)=0._r8
    VLWatSnow_snvr(L,NY,NX)=VLWatSnow_snvr(L,NY,NX)+SnoIce
    VLHeatCapSnow_snvr(L,NY,NX)=cpw*VLWatSnow_snvr(L,NY,NX)
  endif
  TKSnow_snvr(L,NY,NX)=(ENGY+dHPhaseChange)/VLHeatCapSnow_snvr(L,NY,NX)  
!  print*,'Hightemp',TKSnow_snvr(L,NY,NX)
  end subroutine DealHighTempSnow
!------------------------------------------------------------------------------------------

  subroutine DealNegativeSnowMass(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: dHPhaseChange,VLDrySnoWEtmp
  real(r8) :: SnoWEtot,ENGYW,dVLDrySnoWEtmp
  real(r8) :: tKNew,vlheatnew,dIce,dHeat,fs,fi
  real(r8) :: dwat,cphwat,vcphsnw


  cphwat=cpw*VLWatSnow_snvr(L,NY,NX)
  vcphsnw=cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
  if(abs(cphwat/vcphsnw)<1.e-3_r8)return

  !total snow mass
  SnoWEtot=VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE  
  if(TKSnow_snvr(L,NY,NX)<TFICE)then
!     write(*,*)'negative water',L,TKSnow_snvr(L,NY,NX)
!    call endrun('Temeprature too low to resolve negative water '//trim(mod_filename)//' at line',__LINE__)         
  endif
  if(SnoWEtot<0._r8)then
    call endrun('Negative snow mass '//trim(mod_filename)//' at line',__LINE__)         
  endif
  !the starting enthalpy
  ENGYW=VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)  

  !thaw all ice + snow, absorb heat/cooling (<0)
!  if(I>=138.and.I<=139)write(149,*)'neg',I+J/24.,VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX),&
!    TKSnow_snvr(L,NY,NX),(CumHeat2SnowL_snvr(L,NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX)+ENGYW)/(cpi*TFICE)
  dHPhaseChange=-LtHeatIceMelt*(VLDrySnoWE_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE)
  
  vlheatnew=cpw*SnoWEtot
  !compute potential temperature
  tkNew=(ENGYW+CumHeat2SnowL_snvr(L,NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX)+dHPhaseChange)/vlheatnew

  if(tkNew>TFICE)then
    !all ice & snow are melt
    XPhaseChangeHeatL_snvr(L,NY,NX)=XPhaseChangeHeatL_snvr(L,NY,NX)+dHPhaseChange
    VLWatSnow_snvr(L,NY,NX)=SnoWEtot
    VLDrySnoWE_snvr(L,NY,NX)=0._r8
    VLIceSnow_snvr(L,NY,NX)=0._r8
  else
    !some snow and ice remains, refreeze/release heat
    dHeat=-(tkNew-TFICE)*vlheatnew    !>0, enthalpy excess for freeze thaw
    VLDrySnoWEtmp=AZMAX1(VLDrySnoWE_snvr(L,NY,NX))+AZMAX1(VLIceSnow_snvr(L,NY,NX))*DENSICE    
    fs=AZMAX1(VLDrySnoWE_snvr(L,NY,NX))/VLDrySnoWEtmp
    fi=(1._r8-fs)/DENSICE
    dIce=dHeat/(LtHeatIceMelt+(cpw-cpi*fi-cps*fs)*TFICE)          !>0

    if(dIce <= SnoWEtot)then
      if(VLDrySnoWEtmp>0._r8)then      
        dVLDrySnoWEtmp=dIce/VLDrySnoWEtmp
        VLDrySnoWE_snvr(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLDrySnoWE_snvr(L,NY,NX))
        VLIceSnow_snvr(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLIceSnow_snvr(L,NY,NX))
        VLWatSnow_snvr(L,NY,NX)=AZMAX1(SnoWEtot-VLDrySnoWE_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX)*DENSICE)
        !freeze releases heat
        dHPhaseChange=(VLDrySnoWE_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX)*DENSICE-VLDrySnoWEtmp)*LtHeatIceMelt
      else
        !meaning water is positive, freezing occured, heat release
        VLIceSnow_snvr(L,NY,NX)=dIce/DENSICE
        VLDrySnoWE_snvr(L,NY,NX)=0._r8
        VLWatSnow_snvr(L,NY,NX)=AZMAX1(SnoWEtot-VLDrySnoWE_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX)*DENSICE)
        dHPhaseChange=dIce*LtHeatIceMelt
      endif

      XPhaseChangeHeatL_snvr(L,NY,NX)=XPhaseChangeHeatL_snvr(L,NY,NX)+dHPhaseChange
    else
      !negative water is made by ice, no phase change, but there is conversion
      dVLDrySnoWEtmp=SnoWEtot/VLDrySnoWEtmp
      VLDrySnoWE_snvr(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLDrySnoWE_snvr(L,NY,NX))
      VLIceSnow_snvr(L,NY,NX)=dVLDrySnoWEtmp*AZMAX1(VLIceSnow_snvr(L,NY,NX))  
      VLWatSnow_snvr(L,NY,NX)=0._r8    
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
  logical :: watchk,icechk
  real(r8), parameter :: tinyw=1.e-13_r8
! begin_execution

  !the line below is a hack, and likely a better snow layering scheme is needed.
  if(L.eq.1.and.isclose(TCSnow_snvr(1,NY,NX),spval))then
    TCSnow_snvr(1,NY,NX)=units%Kelvin2Celcius(TairK_col(NY,NX))    
  endif
!
! ADD CHANGES IN SNOW, WATER AND ICE
!
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! CumSno2SnowLay,CumWat2SnowLay,CumIce2SnowLay=net fluxes of snow,water,ice in snowpack
! XSnowThawMassL,XIceThawMassL=freeze-thaw flux from watsub.f
!
!
  vwat                     = VLWatSnow_snvr(L,NY,NX);vice = VLIceSnow_snvr(L,NY,NX);vdry = VLDrySnoWE_snvr(L,NY,NX)
  VLDrySnoWE_snvr(L,NY,NX) = AZMAX1d(VLDrySnoWE_snvr(L,NY,NX)+CumSno2SnowL_snvr(L,NY,NX)-XSnowThawMassL_snvr(L,NY,NX),tinyw)
  VLWatSnow_snvr(L,NY,NX)  = AZMAX1d(VLWatSnow_snvr(L,NY,NX)+CumWat2SnowL_snvr(L,NY,NX)+XSnowThawMassL_snvr(L,NY,NX)+XIceThawMassL_snvr(L,NY,NX),tinyw)
  VLIceSnow_snvr(L,NY,NX)  = AZMAX1d(VLIceSnow_snvr(L,NY,NX)+CumIce2SnowL_snvr(L,NY,NX)-XIceThawMassL_snvr(L,NY,NX)/DENSICE,tinyw)

!  if(I==152 .and. J==1 .and. NX==7)then
!    write(113,*)'pdateSnowLayerL',I+J/24.,L,NY,NX
!    write(113,*)'oldwat',vdry,vwat,vice,SnoXfer2SnoLay_snvr(L,NY,NX),VLHeatCapSnow_snvr(L,NY,NX)  
!    write(113,*)'newat',VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)
!    write(113,*)'heat incom, fase',CumHeat2SnowL_snvr(L,NY,NX),XPhaseChangeHeatL_snvr(L,NY,NX)    
!    write(113,*)'diffsnow',CumSno2SnowL_snvr(L,NY,NX)-XSnowThawMassL_snvr(L,NY,NX),&
!          CumWat2SnowL_snvr(L,NY,NX)+XSnowThawMassL_snvr(L,NY,NX)+XIceThawMassL_snvr(L,NY,NX),&
!          CumIce2SnowL_snvr(L,NY,NX)-XIceThawMassL_snvr(L,NY,NX)/DENSICE    
!    write(113,*)'diffcomp',CumSno2SnowL_snvr(L,NY,NX),XSnowThawMassL_snvr(L,NY,NX),&
!          CumWat2SnowL_snvr(L,NY,NX),XIceThawMassL_snvr(L,NY,NX),&
!          CumIce2SnowL_snvr(L,NY,NX)
!    write(113,*)'---------------------------'      
!  endif

  if(any((/VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)/)<0._r8))then      
    write(113,*)I+J/24.,J,L,NY,NX,VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)
    call endrun(trim(mod_filename)//' at line',__LINE__)    
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
    VOLSWI=VOLSWI+(VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE)
    if(VOLSWI<0._r8)then
      write(*,*)'VOLSWI=',VOLSWI,VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
!
!   RESET SNOW SURFACE DENSITY FOR SNOWFALL
!
    IF(SnoXfer2SnoLay_snvr(L,NY,NX).GT.0.0_r8)THEN
      !DENSF: fresh snow density    
      DENSX = SnoDens_snvr(L,NY,NX)
      TCASF = AMAX1(-15.0_r8,AMIN1(2.0_r8,TCA_col(NY,NX)))
      DENSF = 0.05_r8+1.7E-03_r8*(TCASF+15.0_r8)**1.5_r8
      VOLSF = AMIN1(SnoXfer2SnoLay_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX))/DENSF + &
        AZMAX1(VLDrySnoWE_snvr(L,NY,NX)-SnoXfer2SnoLay_snvr(L,NY,NX))/SnoDens_snvr(L,NY,NX)        
      if(VOLSF>0._r8)SnoDens_snvr(L,NY,NX)=VLDrySnoWE_snvr(L,NY,NX)/VOLSF
      !write(*,*)'xVOLSSL=',VLDrySnoWE_snvr(L,NY,NX),SnoXfer2SnoLay_snvr(L,NY,NX),SnoDens_snvr(L,NY,NX),VOLSF
    ENDIF
  ELSE
    VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE_snvr(L-1,NY,NX)+VLWatSnow_snvr(L-1,NY,NX) &
      +VLIceSnow_snvr(L-1,NY,NX)*DENSICE+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX) &
      +VLIceSnow_snvr(L,NY,NX)*DENSICE)
      
    if(VOLSWI<0._r8)then
      write(*,*)'iVOLSWI=',VOLSWI,VLDrySnoWE_snvr(L-1,NY,NX)+VLWatSnow_snvr(L-1,NY,NX) &
        +VLIceSnow_snvr(L-1,NY,NX)*DENSICE,VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX) &
        +VLIceSnow_snvr(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDIF
!
! SNOWPACK COMPRESSION
!
! DDENS1,DDENS2=temperature, compression effect on snow density
! DENSS=snow density in layer
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! VLSnoDWIprev_snvr=snowpack layer volume
! DLYRS=snowpack layer depth
! cumSnowDepz_col=cumulative depth to bottom of snowpack layer
! VHCPW=snowpack layer heat capacity
! TKW,TCSnow=snowpack layer temperature K,oC
! CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
! XPhaseChangeHeatL=latent heat flux from freeze-thaw from watsub.f
! HEATIN_lnd=cumulative net surface heat transfer
! VOLSS,VOLWS,VOLIS=total snow water equivalent, water, ice content of snowpack
! VOLS,SnowDepth=total snowpack volume, depth
!
!  if(I>=138.and.I<=139)print*,L,SnoDens_snvr(L,NY,NX),TKSnow_snvr(L,NY,NX)
  if(TKSnow_snvr(L,NY,NX)<0.)call endrun('too low snow temp '//trim(mod_filename)//' at line',__LINE__)   
  if(SnoDens_snvr(L,NY,NX)>0._r8)then
    IF(SnoDens_snvr(L,NY,NX).LT.0.25_r8)THEN
      DDENS1=SnoDens_snvr(L,NY,NX)*1.0E-05_r8*EXP(0.04_r8*TCSnow_snvr(L,NY,NX))
    ELSE
      DDENS1=0.0_r8
    ENDIF
    CVISC=0.25_r8*EXP(-0.08_r8*TCSnow_snvr(L,NY,NX)+23.0_r8*SnoDens_snvr(L,NY,NX))

    DDENS2=SnoDens_snvr(L,NY,NX)*VOLSWI/(AREA(3,NU(NY,NX),NY,NX)*CVISC)
    
    SnoDens_snvr(L,NY,NX)=SnoDens_snvr(L,NY,NX)+DDENS1+DDENS2
    if(SnoDens_snvr(L,NY,NX)<0._r8)then
      write(*,*)'DDENS1=',SnoDens_snvr(L,NY,NX),DDENS1,DDENS2,L
      write(*,*)SnoXfer2SnoLay_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX)
      call endrun("negative snow dens")
    endif  
  endif
  
  !there is snow in layer L
  VLSnoDWIprev_snvr(L,NY,NX)=VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)

  IF(VLSnoDWIprev_snvr(L,NY,NX) > ZEROS2(NY,NX))THEN
    SnowThickL_snvr(L,NY,NX)    = AZMAX1(VLSnoDWIprev_snvr(L,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
    cumSnowDepz_col(L,NY,NX)    = cumSnowDepz_col(L-1,NY,NX)+SnowThickL_snvr(L,NY,NX)
    VHCPWZ(L,NY,NX)             = VLHeatCapSnow_snvr(L,NY,NX)
    TKWX                        = TKSnow_snvr(L,NY,NX)
    ENGYW                       = VLHeatCapSnow_snvr(L,NY,NX)*TKWX
    VLHeatCapSnow_snvr(L,NY,NX) = cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)

    IF(VLHeatCapSnow_snvr(L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      !there is significant snow layer mass
      TKSnow_snvr(L,NY,NX)=(ENGYW+CumHeat2SnowL_snvr(L,NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX))/VLHeatCapSnow_snvr(L,NY,NX)
!      if(I==152 .and. J==1 .and. NX==7)then
!        write(113,*)'HHeatcap old new',VHCPWZ(L,NY,NX),VLHeatCapSnow_snvr(L,NY,NX),ENGYW+CumHeat2SnowL_snvr(L,NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX)
!        write(113,*)'sepengy',ENGYW,CumHeat2SnowL_snvr(L,NY,NX),XPhaseChangeHeatL_snvr(L,NY,NX)
!        write(*,*)'LL',L,TKWX,TKSnow_snvr(L,NY,NX),VLHeatCapSnow_snvr(L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)
!      endif
      if((TKSnow_snvr(L,NY,NX)<200._r8 .or. TKSnow_snvr(L,NY,NX)>290._r8) .and. VLHeatCapSnow_snvr(L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))then
        write(*,*)'ijl',I,J,L,NY,NX,TKSnow_snvr(L,NY,NX),TKWX,TairK_col(NY,NX)
        write(*,*)VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)
        icechk=TKSnow_snvr(L,NY,NX)>290._r8 .and. VLDrySnoWE_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)>1.e-10_r8
        watchk=TKSnow_snvr(L,NY,NX)<200._r8 .and. VLWatSnow_snvr(L,NY,NX) > 1.e-10_r8
        if(icechk .or. watchk)call endrun(trim(mod_filename)//' at line',__LINE__)     
      endif
      if(TKSnow_snvr(L,NY,NX)>280._r8)call DealHighTempSnow(I,J,L,NY,NX)
    ELSE
      !there is no significant snow mass      
      IF(L.EQ.1)THEN
        !if current layer is top layer
        TKSnow_snvr(L,NY,NX)=TairK_col(NY,NX)
      ELSE
        !if it is not the top layer
        TKSnow_snvr(L,NY,NX)=TKSnow_snvr(L-1,NY,NX)
      ENDIF
      IF(VHCPWZ(L,NY,NX).GT.ZEROS(NY,NX))THEN
        HEATIN_lnd=HEATIN_lnd+(TKSnow_snvr(L,NY,NX)-TKWX)*VLHeatCapSnow_snvr(L,NY,NX)
      ENDIF
    ENDIF
  ELSE
    !there is no snow in layer L
    VLDrySnoWE_snvr(L,NY,NX)   = 0.0_r8
    VLWatSnow_snvr(L,NY,NX)    = 0.0_r8
    VLIceSnow_snvr(L,NY,NX)    = 0.0_r8
    VLSnoDWIprev_snvr(L,NY,NX) = 0.0_r8
    SnowThickL_snvr(L,NY,NX)   = 0.0_r8
    cumSnowDepz_col(L,NY,NX)   = cumSnowDepz_col(L-1,NY,NX)
    IF(L.EQ.1)THEN
      TKSnow_snvr(L,NY,NX)=TairK_col(NY,NX)
    ELSE
      TKSnow_snvr(L,NY,NX)=TKSnow_snvr(L-1,NY,NX)
    ENDIF
  ENDIF  
  TCSnow_snvr(L,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L,NY,NX))

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
    trcg_solsml_snvr(NTG,L,NY,NX)=trcg_solsml_snvr(NTG,L,NY,NX)+trcg_TBLS(NTG,L,NY,NX)
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

  if(lverb)write(*,*)'SnowpackLayering'
!  if(I>=138.and.I<=139)print*,I+J/24.,'bflay',TKSnow_snvr(1:JS,NY,NX)
  IF(VLHeatCapSnow_snvr(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
    D325: DO L=1,JS-1
!      VOLSLX=VLSnoDWIprev_snvr(L,NY,NX)
      IF(VLSnoDWIprev_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        !compute the free volume/thickness: DDLYXS
        DDLYXS=(VLSnoDWIMax_col(L,NY,NX)-VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX) &
          -VLWatSnow_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX))/AREA(3,L,NY,NX)
!        DDLYXX=DDLYXS
        IF(DDLYXS.LT.-ZERO .OR. SnowThickL_snvr(L+1,NY,NX).GT.ZERO)THEN
          !current volume is greater than allowed, or next layer exists   
          !case 1: DDLYXS< 0, layer L extends into layer L+1, DDLYRS<0: amount of expand layer L into layer L+1
          !case 2: 0<DDLYXS<SnowThickL_snvr(L+1,NY,NX), layer L still has space, can take DDLYRS of layer L+1       
          !case 3: DDLYXS>SnowThickL_snvr(L+1,NY,NX), layer L still has space, and can even hold DDLYRS of layer L+1
          !
          DDLYRS=AMIN1(DDLYXS,SnowThickL_snvr(L+1,NY,NX))
          IFLGLS=iexpand         !expand
        ELSE
          !layer L has space and layer L+1 disappears, so layer L+1 combines into layer L
          !volume less than allowed, and no next layer
          !DDLYXS: is the depth change of layer L
          DDLYXS=(VLSnoDWIprev_snvr(L,NY,NX)-VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX) &
            -VLWatSnow_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX))/AREA(3,L,NY,NX)
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
      SnowThickL_snvr(L,NY,NX)=cumSnowDepz_col(L,NY,NX)-cumSnowDepz_col(L-1,NY,NX)
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
            FX=AMIN1(1.0_r8,DDLYRS*AREA(3,L0,NY,NX)/VLSnoDWIprev_snvr(L0,NY,NX))
          ENDIF
        ELSE
          !expanding L into L+1
          L1=L+1
          L0=L
          IF(VLSnoDWIprev_snvr(L0,NY,NX).LT.VLSnoDWIMax_col(L0,NY,NX))THEN
            FX=0.0_r8
          ELSE
            !FX fraction to be donated from L0 to L1
            FX=AMIN1(1.0_r8,-DDLYRS*AREA(3,L0,NY,NX)/VLSnoDWIprev_snvr(L0,NY,NX))
          ENDIF
        ENDIF
!   donor L0, target L1
        IF(FX.GT.0.0_r8)THEN
          FY=1.0_r8-FX
!          print*,'relay',L0,L1,FX,FY,DDLYRS,DDLYXS,-ZERO,SnowThickL_snvr(L+1,NY,NX),VLSnoDWIprev_snvr(L0,NY,NX),VLSnoDWIMax_col(L0,NY,NX)
!
!     TARGET SNOW LAYER
!         volume/mass
          VLDrySnoWE_snvr(L1,NY,NX)=VLDrySnoWE_snvr(L1,NY,NX)+FX*VLDrySnoWE_snvr(L0,NY,NX)
          VLWatSnow_snvr(L1,NY,NX)=VLWatSnow_snvr(L1,NY,NX)+FX*VLWatSnow_snvr(L0,NY,NX)
          VLIceSnow_snvr(L1,NY,NX)=VLIceSnow_snvr(L1,NY,NX)+FX*VLIceSnow_snvr(L0,NY,NX)          
          VLSnoDWIprev_snvr(L1,NY,NX)=VLDrySnoWE_snvr(L1,NY,NX)/SnoDens_snvr(L1,NY,NX)+VLWatSnow_snvr(L1,NY,NX)+VLIceSnow_snvr(L1,NY,NX)
!         energy
          ENGY1X=VLHeatCapSnow_snvr(L1,NY,NX)*TKSnow_snvr(L1,NY,NX)
          ENGY0X=VLHeatCapSnow_snvr(L0,NY,NX)*TKSnow_snvr(L0,NY,NX)
          ENGY1=ENGY1X+FX*ENGY0X
          VLHeatCapSnow_snvr(L1,NY,NX)=cps*VLDrySnoWE_snvr(L1,NY,NX)+cpw*VLWatSnow_snvr(L1,NY,NX)+cpi*VLIceSnow_snvr(L1,NY,NX)

          IF(VLHeatCapSnow_snvr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
            TKSnow_snvr(L1,NY,NX)=ENGY1/VLHeatCapSnow_snvr(L1,NY,NX)
          ELSE
            TKSnow_snvr(L1,NY,NX)=TKSnow_snvr(L0,NY,NX)
          ENDIF
          TCSnow_snvr(L1,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L1,NY,NX))
!          chemicals
          !gas
          DO NTG=idg_beg,idg_end-1
            trcg_solsml_snvr(NTG,L1,NY,NX)=trcg_solsml_snvr(NTG,L1,NY,NX)+FX*trcg_solsml_snvr(NTG,L0,NY,NX)
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
          VLDrySnoWE_snvr(L0,NY,NX)=FY*VLDrySnoWE_snvr(L0,NY,NX)
          VLWatSnow_snvr(L0,NY,NX)=FY*VLWatSnow_snvr(L0,NY,NX)
          VLIceSnow_snvr(L0,NY,NX)=FY*VLIceSnow_snvr(L0,NY,NX)
          VLSnoDWIprev_snvr(L0,NY,NX)=VLDrySnoWE_snvr(L0,NY,NX)/SnoDens_snvr(L0,NY,NX)+VLWatSnow_snvr(L0,NY,NX)+VLIceSnow_snvr(L0,NY,NX)
!     energy 
          ENGY0=FY*ENGY0X
          VLHeatCapSnow_snvr(L0,NY,NX)=cps*VLDrySnoWE_snvr(L0,NY,NX)+cpw*VLWatSnow_snvr(L0,NY,NX)+cpi*VLIceSnow_snvr(L0,NY,NX)
          IF(VLHeatCapSnow_snvr(L0,NY,NX).GT.ZEROS(NY,NX))THEN
            TKSnow_snvr(L0,NY,NX)=ENGY0/VLHeatCapSnow_snvr(L0,NY,NX)
          ELSE
            TKSnow_snvr(L0,NY,NX)=TKSnow_snvr(L1,NY,NX)
          ENDIF
          TCSnow_snvr(L0,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L0,NY,NX))
!     chemicals
          DO NTG=idg_beg,idg_end-1
            trcg_solsml_snvr(NTG,L0,NY,NX)=FY*trcg_solsml_snvr(NTG,L0,NY,NX)
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
!     if(I>=138.and.I<140)print*,I+J/24.,'layering',L,VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX),TKSnow_snvr(L,NY,NX)
    if(TKSnow_snvr(L,NY,NX)/=spval)then
      VLHeatCapSnow_snvr(L,NY,NX)=cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
      if(VLHeatCapSnow_snvr(L,NY,NX)>ZEROS(NY,NX))then
        nsnol_col(NY,NX)=nsnol_col(NY,NX)+1
      endif
    else
      VLDrySnoWE_snvr(L,NY,NX)=0._r8
      VLWatSnow_snvr(L,NY,NX)=0._r8
      VLIceSnow_snvr(L,NY,NX)=0._r8
      VLHeatCapSnow_snvr(L,NY,NX)=0._r8
      SnoDens_snvr(L,NY,NX)=NewSnowDens_col(NY,NX)
    endif
    if(VLHeatCapSnow_snvr(L,NY,NX)<ZEROS(NY,NX))then
      if(L>1)then
        TKSnow_snvr(L,NY,NX)=spval        
      else
        TKSnow_snvr(L,NY,NX)=TairK_col(NY,NX) 
      endif  
    endif  
    !move up to handel loss of surface layer, this should rarely occur, but round off error may trigger it
    if(L>1 .and. VLHeatCapSnow_snvr(L,NY,NX)>=ZEROS(NY,NX))then
      if (TKSnow_snvr(L-1,NY,NX)==spval)then
        VLDrySnoWE_snvr(L-1,NY,NX)=VLDrySnoWE_snvr(L,NY,NX)
        VLWatSnow_snvr(L-1,NY,NX)=VLWatSnow_snvr(L,NY,NX)
        VLIceSnow_snvr(L-1,NY,NX)=VLIceSnow_snvr(L,NY,NX)
        VLHeatCapSnow_snvr(L-1,NY,NX)=VLHeatCapSnow_snvr(L,NY,NX)
        TKSnow_snvr(L-1,NY,NX)=TKSnow_snvr(L,NY,NX)
        SnoDens_snvr(L-1,NY,NX)=SnoDens_snvr(L,NY,NX)

        VLDrySnoWE_snvr(L,NY,NX)=0._r8
        VLWatSnow_snvr(L,NY,NX)=0._r8
        VLIceSnow_snvr(L,NY,NX)=0._r8
        VLHeatCapSnow_snvr(L,NY,NX)=0._r8
        TKSnow_snvr(L,NY,NX)=spval
        SnoDens_snvr(L,NY,NX)=NewSnowDens_col(NY,NX)
      endif
    endif 
  ENDDO  
  
  end subroutine SnowpackLayering

!------------------------------------------------------------------------------------------
  subroutine ZeroSnowArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  TDrysnoBySnowRedist(NY,NX)   = 0.0_r8
  TWatBySnowRedist(NY,NX)      = 0.0_r8
  TIceBySnowRedist(NY,NX)      = 0.0_r8
  THeatBySnowRedist_col(NY,NX) = 0.0_r8

  trcn_TFloXSurRunoff_2D(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8
  trcg_QSS(idg_beg:idg_end-1,NY,NX)                      = 0.0_r8
  trcn_QSS(ids_nut_beg:ids_nuts_end,NY,NX)               = 0.0_r8
  trcg_TFloXSurRunoff(idg_beg:idg_end-1,NY,NX)           = 0.0_r8

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
