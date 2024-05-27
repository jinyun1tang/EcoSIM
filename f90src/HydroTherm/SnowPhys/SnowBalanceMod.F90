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
  public :: SaltFromRunoffSnowpack
  public :: OverlandSnowFlow
  public :: ChemicalBySnowRedistribution
  public :: MassFluxFromSnowRunoff
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

    call UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)

    call UpdateSoluteInSnow(L,NY,NX)
  ENDDO D9780

!
!     SNOW RUNOFF from snow laye near soil surface
!
  VLDrySnoWE(1,NY,NX)=VLDrySnoWE(1,NY,NX)+TDrysnoBySnowRedist(NY,NX)
  VLWatSnow(1,NY,NX)=VLWatSnow(1,NY,NX)+TWatBySnowRedist(NY,NX)
  VLIceSnow(1,NY,NX)=VLIceSnow(1,NY,NX)+TIceBySnowRedist(NY,NX)
  ENGYW=VLHeatCapSnow_col(1,NY,NX)*TKSnow(1,NY,NX)
  VLHeatCapSnow_col(1,NY,NX)=cps*VLDrySnoWE(1,NY,NX)+cpw*VLWatSnow(1,NY,NX)+cpi*VLIceSnow(1,NY,NX)

  IF(VLHeatCapSnow_col(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
    TKSnow(1,NY,NX)=(ENGYW+THeatBySnowRedist(NY,NX))/VLHeatCapSnow_col(1,NY,NX)
  ELSE
    TKSnow(1,NY,NX)=TairK(NY,NX)
    if(VcumSnowWE(NY,NX)<ZEROS(NY,NX))TKSnow(1,NY,NX)=spval
  ENDIF
  write(133,*)I+J/24.,'snow',VcumSnowWE(NY,NX),TKSnow(1,NY,NX),TairK(NY,NX),THeatBySnowRedist(NY,NX)
  VcumDrySnoWE_col(NY,NX)=sum(VLDrySnoWE(1:JS,NY,NX))
  VcumWatSnow_col(NY,NX)=sum(VLWatSnow(1:JS,NY,NX))
  VcumIceSnow_col(NY,NX)=sum(VLIceSnow(1:JS,NY,NX))
  VcumSnoDWI(NY,NX)=sum(VLSnoDWIprev_col(1:JS,NY,NX))
  SnowDepth(NY,NX)=sum(SnowThickL_col(1:JS,NY,NX))
  VcumSnowWE(NY,NX)=VcumDrySnoWE_col(NY,NX)+VcumIceSnow_col(NY,NX)*DENSICE+VcumWatSnow_col(NY,NX) 
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
    FLWS=VLDrySnoWE(1,NY,NX)
    FLWW=VLWatSnow(1,NY,NX)
    FLWI=VLIceSnow(1,NY,NX)
    HeatFlo2Surface=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TKSnow(1,NY,NX)

    !reset snow layer variables
    VLDrySnoWE(1,NY,NX)=0.0_r8
    VLWatSnow(1,NY,NX)=0.0_r8
    VLIceSnow(1,NY,NX)=0.0_r8
    VLHeatCapSnow_col(1,NY,NX)=0.0_r8
    VcumDrySnoWE_col(NY,NX)=0.0_r8
    VcumWatSnow_col(NY,NX)=0.0_r8
    VcumIceSnow_col(NY,NX)=0.0_r8
    VcumSnoDWI(NY,NX)=0.0_r8
    SnowDepth(NY,NX)=0.0_r8    
    VcumSnowWE(NY,NX)=0._r8
    D9770: DO L=1,JS
      SnoDensL(L,NY,NX)=NewSnowDens(NY,NX)
      TKSnow(L,NY,NX)=spval
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

  subroutine UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)
  implicit none
  integer, intent(in) :: L,NY,NX,I,J
  real(r8), intent(inout) :: VOLSWI
  real(r8) :: TKWX,VHCPWZ(JZ,JY,JX)
  real(r8) :: ENGYW
  real(r8) :: VOLSF,TCASF
  real(r8) :: DENSF,DDENS2,DDENS1
  real(r8) :: CVISC,DENSX
  real(r8) :: dHPhaseChange,VLDrySnoWEtmp
  real(r8) :: frcnew,SnowIceMass,VLwatNet
  
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
  VLDrySnoWE(L,NY,NX)=VLDrySnoWE(L,NY,NX)+CumSno2SnowLay(L,NY,NX)-XSnowThawMassL(L,NY,NX)
  VLWatSnow(L,NY,NX)=VLWatSnow(L,NY,NX)+CumWat2SnowLay(L,NY,NX)+XSnowThawMassL(L,NY,NX)+XIceThawMassL(L,NY,NX)
  VLIceSnow(L,NY,NX)=VLIceSnow(L,NY,NX)+CumIce2SnowLay(L,NY,NX)-XIceThawMassL(L,NY,NX)/DENSICE

  if(VLWatSnow(L,NY,NX)<0._r8)then
    !melt dry snow
    VLDrySnoWEtmp=VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX)
    if(VLDrySnoWEtmp>0._r8)then
      VLDrySnoWE(L,NY,NX)=VLDrySnoWEtmp
      dHPhaseChange=-LtHeatIceMelt*VLWatSnow(L,NY,NX)
      XPhaseChangeHeatL(L,NY,NX)=XPhaseChangeHeatL(L,NY,NX)+dHPhaseChange
      VLWatSnow(L,NY,NX)=0._r8      
    else
      !also check for ice
      SnowIceMass=VLDrySnoWE(L,NY,NX)+VLIceSnow(L,NY,NX)*DENSICE
      VLwatNet=SnowIceMass+VLWatSnow(L,NY,NX)
      if(VLwatNet>0._r8)then
        frcnew=VLwatNet/SnowIceMass
        VLDrySnoWE(L,NY,NX)=VLDrySnoWE(L,NY,NX)*frcnew
        VLIceSnow(L,NY,NX)=VLIceSnow(L,NY,NX)*frcnew
        dHPhaseChange=-LtHeatIceMelt*VLWatSnow(L,NY,NX)
        XPhaseChangeHeatL(L,NY,NX)=XPhaseChangeHeatL(L,NY,NX)+dHPhaseChange        
        VLWatSnow(L,NY,NX)=0._r8
      else
        if(abs(VLWatSnow(L,NY,NX)) < ZEROS2(NY,NX))then
          VLWatSnow(L,NY,NX)=0._r8
        else
          write(*,*)'negative snowmass cannot be fixed',L,VLDrySnoWE(L,NY,NX),&
            VLWatSnow(L,NY,NX),VLIceSnow(L,NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)       
        endif   
      endif
    endif
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
    VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX) &
      +VLIceSnow(L,NY,NX)*DENSICE)
    if(VOLSWI<0._r8)then
      write(*,*)'VOLSWI=',VOLSWI,VLDrySnoWE(L,NY,NX),VLWatSnow(L,NY,NX), &
        VLIceSnow(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
!
!   RESET SNOW SURFACE DENSITY FOR SNOWFALL
!
!    print*,'here'
    IF(SnoXfer2SnoLay(L,NY,NX).GT.0.0_r8)THEN
      DENSX=SnoDensL(L,NY,NX)
      TCASF=AMAX1(-15.0_r8,AMIN1(2.0_r8,TCA(NY,NX)))
      !fresh snow density
      DENSF=0.05_r8+1.7E-03_r8*(TCASF+15.0_r8)**1.5_r8
      VOLSF=AMIN1(SnoXfer2SnoLay(L,NY,NX),VLDrySnoWE(L,NY,NX))/DENSF + &
        AZMAX1(VLDrySnoWE(L,NY,NX)-SnoXfer2SnoLay(L,NY,NX))/SnoDensL(L,NY,NX)
      SnoDensL(L,NY,NX)=VLDrySnoWE(L,NY,NX)/VOLSF
      !write(*,*)'xVOLSSL=',VLDrySnoWE(L,NY,NX),SnoXfer2SnoLay(L,NY,NX),SnoDensL(L,NY,NX),VOLSF
    ENDIF
  ELSE
    VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE(L-1,NY,NX)+VLWatSnow(L-1,NY,NX) &
      +VLIceSnow(L-1,NY,NX)*DENSICE+VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX) &
      +VLIceSnow(L,NY,NX)*DENSICE)
    !write(*,*)'VOLSWI',VOLSWI
    if(VOLSWI<0._r8)then
      write(*,*)'iVOLSWI=',VOLSWI,VLDrySnoWE(L-1,NY,NX)+VLWatSnow(L-1,NY,NX) &
        +VLIceSnow(L-1,NY,NX)*DENSICE,VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX) &
        +VLIceSnow(L,NY,NX)*DENSICE
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
  IF(SnoDensL(L,NY,NX).LT.0.25_r8)THEN
    DDENS1=SnoDensL(L,NY,NX)*1.0E-05_r8*EXP(0.04_r8*TCSnow(L,NY,NX))
  ELSE
    DDENS1=0.0_r8
  ENDIF
  CVISC=0.25_r8*EXP(-0.08_r8*TCSnow(L,NY,NX)+23.0_r8*SnoDensL(L,NY,NX))

  DDENS2=SnoDensL(L,NY,NX)*VOLSWI/(AREA(3,NU(NY,NX),NY,NX)*CVISC)
  
  if(DDENS2<0._r8)write(*,*)'DDENS2=',SnoDensL(L,NY,NX),VOLSWI,CVISC
  SnoDensL(L,NY,NX)=SnoDensL(L,NY,NX)+DDENS1+DDENS2
  if(SnoDensL(L,NY,NX)<0._r8)then
    write(*,*)'DDENS1=',SnoDensL(L,NY,NX),DDENS1,DDENS2,L
    write(*,*)SnoXfer2SnoLay(L,NY,NX),VLDrySnoWE(L,NY,NX)
    call endrun("negative snow dens")
  endif  

  !there is snow in layer L
  IF(VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX)+VLIceSnow(L,NY,NX).GT.ZEROS2(NY,NX) .and. SnoDensL(L,NY,NX) > 0._r8)THEN
    VLSnoDWIprev_col(L,NY,NX)=VLDrySnoWE(L,NY,NX)/SnoDensL(L,NY,NX)+VLWatSnow(L,NY,NX)+VLIceSnow(L,NY,NX)
    SnowThickL_col(L,NY,NX)=AZMAX1(VLSnoDWIprev_col(L,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
    cumSnowDepz_col(L,NY,NX)=cumSnowDepz_col(L-1,NY,NX)+SnowThickL_col(L,NY,NX)
    VHCPWZ(L,NY,NX)=VLHeatCapSnow_col(L,NY,NX)
    TKWX=TKSnow(L,NY,NX)
    ENGYW=VLHeatCapSnow_col(L,NY,NX)*TKSnow(L,NY,NX)
    VLHeatCapSnow_col(L,NY,NX)=cps*VLDrySnoWE(L,NY,NX)+cpw*VLWatSnow(L,NY,NX)+cpi*VLIceSnow(L,NY,NX)    
    IF(VHCPWZ(L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX) .AND. VLHeatCapSnow_col(L,NY,NX).GT.ZEROS(NY,NX))THEN
      !there is significant snow layer mass
      TKSnow(L,NY,NX)=(ENGYW+CumHeat2SnowLay(L,NY,NX)+XPhaseChangeHeatL(L,NY,NX))/VLHeatCapSnow_col(L,NY,NX)
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
    VLDrySnoWE(L,NY,NX)=0.0_r8
    VLWatSnow(L,NY,NX)=0.0_r8
    VLIceSnow(L,NY,NX)=0.0_r8
    VLSnoDWIprev_col(L,NY,NX)=0.0_r8
    SnowThickL_col(L,NY,NX)=0.0_r8
    cumSnowDepz_col(L,NY,NX)=cumSnowDepz_col(L-1,NY,NX)
    VLHeatCapSnow_col(L,NY,NX)=0.0_r8
    IF(L.EQ.1)THEN
      TKSnow(L,NY,NX)=TairK(NY,NX)
    ELSE
      TKSnow(L,NY,NX)=TKSnow(L-1,NY,NX)
    ENDIF
  ENDIF
!  print*,I+J/24.,'snowmass up',L,VLSnoDWIprev_col(L,NY,NX),SnoDensL(L,NY,NX),&
!    VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX)+VLIceSnow(L,NY,NX)
  TCSnow(L,NY,NX)=units%Kelvin2Celcius(TKSnow(L,NY,NX))
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
  IF(VLHeatCapSnow_col(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
    D325: DO L=1,JS-1
!      VOLSLX=VLSnoDWIprev_col(L,NY,NX)
!      print*,I+J/24.,'layer L',L,VLSnoDWIprev_col(L,NY,NX),SnoDensL(L,NY,NX)
      IF(VLSnoDWIprev_col(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        !compute the free volume/thickness: DDLYXS
        DDLYXS=(VLSnoDWIMax_col(L,NY,NX)-VLDrySnoWE(L,NY,NX)/SnoDensL(L,NY,NX) &
          -VLWatSnow(L,NY,NX)-VLIceSnow(L,NY,NX))/AREA(3,L,NY,NX)
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
          DDLYXS=(VLSnoDWIprev_col(L,NY,NX)-VLDrySnoWE(L,NY,NX)/SnoDensL(L,NY,NX) &
            -VLWatSnow(L,NY,NX)-VLIceSnow(L,NY,NX))/AREA(3,L,NY,NX)
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
!        print*,'relay',L0,L1,FX,DDLYRS,DDLYXS,-ZERO,SnowThickL_col(L+1,NY,NX)
        IF(FX.GT.0.0_r8)THEN
          FY=1.0_r8-FX
!
!     TARGET SNOW LAYER
!         volume/mass
          VLDrySnoWE(L1,NY,NX)=VLDrySnoWE(L1,NY,NX)+FX*VLDrySnoWE(L0,NY,NX)
          VLWatSnow(L1,NY,NX)=VLWatSnow(L1,NY,NX)+FX*VLWatSnow(L0,NY,NX)
          VLIceSnow(L1,NY,NX)=VLIceSnow(L1,NY,NX)+FX*VLIceSnow(L0,NY,NX)
          if(SnoDensL(L1,NY,NX)<ZEROS(NY,NX))SnoDensL(L1,NY,NX)=SnoDensL(L0,NY,NX)
          VLSnoDWIprev_col(L1,NY,NX)=VLDrySnoWE(L1,NY,NX)/SnoDensL(L1,NY,NX)+VLWatSnow(L1,NY,NX)+VLIceSnow(L1,NY,NX)
!         energy
          ENGY1X=VLHeatCapSnow_col(L1,NY,NX)*TKSnow(L1,NY,NX)
          ENGY0X=VLHeatCapSnow_col(L0,NY,NX)*TKSnow(L0,NY,NX)
          ENGY1=ENGY1X+FX*ENGY0X
          VLHeatCapSnow_col(L1,NY,NX)=cps*VLDrySnoWE(L1,NY,NX)+cpw*VLWatSnow(L1,NY,NX)+cpi*VLIceSnow(L1,NY,NX)

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
          VLDrySnoWE(L0,NY,NX)=FY*VLDrySnoWE(L0,NY,NX)
          VLWatSnow(L0,NY,NX)=FY*VLWatSnow(L0,NY,NX)
          VLIceSnow(L0,NY,NX)=FY*VLIceSnow(L0,NY,NX)
          VLSnoDWIprev_col(L0,NY,NX)=VLDrySnoWE(L0,NY,NX)/SnoDensL(L0,NY,NX)+VLWatSnow(L0,NY,NX)+VLIceSnow(L0,NY,NX)
!     energy 
          ENGY0=FY*ENGY0X
          VLHeatCapSnow_col(L0,NY,NX)=cps*VLDrySnoWE(L0,NY,NX)+cpw*VLWatSnow(L0,NY,NX)+cpi*VLIceSnow(L0,NY,NX)
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
!     IF(VLWatSnow(L0,NY,NX)+VLIceSnow(L0,NY,NX)
!    2+VLDrySnoWE(L0,NY,NX).LE.ZEROS(NY,NX))THEN
!     cumSnowDepz_col(L1,NY,NX)=cumSnowDepz_col(L0,NY,NX)
!     ENDIF
        ENDIF
      ENDIF
    ENDDO D325
  ENDIF
  !update volumetric snow heat capacity
  nsnol_col(NY,NX)=0
  DO L=1,JS
    if(TKSnow(L,NY,NX)/=spval)then
      VLHeatCapSnow_col(L,NY,NX)=cps*VLDrySnoWE(L,NY,NX)+cpw*VLWatSnow(L,NY,NX)+cpi*VLIceSnow(L,NY,NX)
      nsnol_col(NY,NX)=nsnol_col(NY,NX)+1
    else
      VLDrySnoWE(L,NY,NX)=0._r8
      VLWatSnow(L,NY,NX)=0._r8
      VLIceSnow(L,NY,NX)=0._r8
    endif
    if(VLHeatCapSnow_col(L,NY,NX)<ZEROS(NY,NX))TKSnow(L,NY,NX)=spval  
  ENDDO  
  end subroutine SnowpackLayering

!------------------------------------------------------------------------------------------

  subroutine SaltFromRunoffSnowpack(N1,N2,NY,NX)
  implicit none
  integer, intent(in) :: N1,N2,NY,NX
  integer :: LS, LS2
  integer :: NTG,NTN,NTSA,NTS
!     begin_execution
!     NET WATER AND HEAT FLUXES THROUGH SNOWPACK
!
!     VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
!     CumSno2SnowLay,CumWat2SnowLay,CumIce2SnowLay=net fluxes of snow,water,ice in snowpack
!     CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
!     XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=snow,water,ice transfer from watsub.f
!     HeatXfer2SnoLay=convective heat flux from snow,water,ice transfer from watsub.f
!     WatConvSno2MicP,WatConvSno2MacP,WatConvSno2LitR=water flux from lowest snow layer to soil macropore,micropore,litter
!     HeatConvSno2Soi,HeatConvSno2LitR=heat flux from lowest snow layer to soil,litter

  D1205: DO LS=1,JS
    IF(VLHeatCapSnow_col(LS,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      !id of next snow layer
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS.AND.VLHeatCapSnow_col(LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        !not surface layer, and is heat significant
        CumSno2SnowLay(LS,N2,N1)=CumSno2SnowLay(LS,N2,N1)+SnoXfer2SnoLay(LS,N2,N1)-SnoXfer2SnoLay(LS2,N2,N1)
        CumWat2SnowLay(LS,N2,N1)=CumWat2SnowLay(LS,N2,N1)+WatXfer2SnoLay(LS,N2,N1)-WatXfer2SnoLay(LS2,N2,N1) &
          -WatConvSno2LitR(LS,N2,N1)-WatConvSno2MicP(LS,N2,N1)-WatConvSno2MacP(LS,N2,N1)
        CumIce2SnowLay(LS,N2,N1)=CumIce2SnowLay(LS,N2,N1)+IceXfer2SnoLay(LS,N2,N1)-IceXfer2SnoLay(LS2,N2,N1)
        CumHeat2SnowLay(LS,N2,N1)=CumHeat2SnowLay(LS,N2,N1)+HeatXfer2SnoLay(LS,N2,N1) &
          -HeatXfer2SnoLay(LS2,N2,N1)-HeatConvSno2LitR(LS,N2,N1)-HeatConvSno2Soi(LS,N2,N1)

        !     NET SOLUTE FLUXES THROUGH SNOWPACK
        !
        !     T*BLS=net solute flux in snowpack
        !     X*BLS=solute flux in snowpack from TranspNoSalt.f
        !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
        !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
        !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
        !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
        !
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_Xbndl_flx(NTG,LS,N2,N1) &
            -trcg_Xbndl_flx(NTG,LS2,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_Xbndl_flx(NTN,LS,N2,N1) &
            -trcn_Xbndl_flx(NTN,LS2,N2,N1)
        ENDDO
        !
        !     NET SALT FLUXES THROUGH SNOWPACK
        !
        !     T*BLS=net solute flux in snowpack
        !     X*BLS=solute flux in snowpack from TranspSalt.f
        !     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
        !          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
        !          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
        !          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
        !          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
        !          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
        !          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
        !     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
        !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
        !          :*1=non-band,*B=band
        !
        IF(salt_model)THEN
          DO NTSA=idsalt_beg,idsalt_end
            trcSalt_TBLS(NTSA,LS,N2,N1)=trcSalt_TBLS(NTSA,LS,N2,N1)+trcSaltFlo2SnowLay(NTSA,LS,N2,N1)&
              -trcSaltFlo2SnowLay(NTSA,LS2,N2,N1)
          ENDDO
        ENDIF
        !
        !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
        !
      ELSE
        CumSno2SnowLay(LS,N2,N1)=CumSno2SnowLay(LS,N2,N1)+SnoXfer2SnoLay(LS,N2,N1)
        CumWat2SnowLay(LS,N2,N1)=CumWat2SnowLay(LS,N2,N1)+WatXfer2SnoLay(LS,N2,N1) &
          -WatConvSno2LitR(LS,N2,N1)-WatConvSno2MicP(LS,N2,N1)-WatConvSno2MacP(LS,N2,N1)
        CumIce2SnowLay(LS,N2,N1)=CumIce2SnowLay(LS,N2,N1)+IceXfer2SnoLay(LS,N2,N1)
        CumHeat2SnowLay(LS,N2,N1)=CumHeat2SnowLay(LS,N2,N1)+HeatXfer2SnoLay(LS,N2,N1) &
          -HeatConvSno2LitR(LS,N2,N1)-HeatConvSno2Soi(LS,N2,N1)

        ! and NH3B
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_Xbndl_flx(NTG,LS,N2,N1) &
            -trcs_3DTransp2MicP_vr(NTG,3,0,N2,N1)-trcs_3DTransp2MicP_vr(NTG,3,NUM(N2,N1),N2,N1) &
            -trcs_3DTransp2MacP(NTG,3,NUM(N2,N1),N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_Xbndl_flx(NTN,LS,N2,N1) &
            -trcs_3DTransp2MicP_vr(NTN,3,0,N2,N1)-trcs_3DTransp2MicP_vr(NTN,3,NUM(N2,N1),N2,N1) &
            -trcs_3DTransp2MacP(NTN,3,NUM(N2,N1),N2,N1)
        ENDDO

        !add band flux
        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1) &
          -trcs_3DTransp2MicP_vr(idg_NH3B,3,NUM(N2,N1),N2,N1)-trcs_3DTransp2MacP(idg_NH3B,3,NUM(N2,N1),N2,N1)

        DO NTS=0,ids_nuts
          trcn_TBLS(ids_NH4+NTS,LS,N2,N1)=trcn_TBLS(ids_NH4+NTS,LS,N2,N1) &
            -trcs_3DTransp2MicP_vr(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)-trcs_3DTransp2MacP(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO NTSA=idsalt_beg,idsalt_end
            trcSalt_TBLS(NTSA,LS,NY,NX)=trcSalt_TBLS(NTSA,LS,NY,NX)+trcSaltFlo2SnowLay(NTSA,LS,NY,NX) &
              -trcSalt3DFlo2Cell(NTSA,3,0,N2,N1)-trcSalt3DFlo2Cell(NTSA,3,NUM(N2,N1),N2,N1) &
              -trcSalt_XFHS(NTSA,3,NUM(N2,N1),N2,N1)
          ENDDO

          !add band flux
          DO NTSA=0,idsalt_nuts
            trcSalt_TBLS(idsalt_H0PO4+NTSA,LS,NY,NX)=trcSalt_TBLS(idsalt_H0PO4+NTSA,LS,NY,NX) &
              -trcSalt3DFlo2Cell(idsalt_H0PO4B+NTSA,3,NUM(N2,N1),N2,N1) &
              -trcSalt_XFHS(idsalt_H0PO4B+NTSA,3,NUM(N2,N1),N2,N1)
          ENDDO
        ENDIF
      ENDIF
!
!     WATER,GAS,SOLUTE,SALT FLUXES INTO SNOWPACK SURFACE
!
    ELSEIF(LS.EQ.1)THEN
      IF(abs(SnoXfer2SnoLay(LS,N2,N1))>0._r8)THEN
        CumSno2SnowLay(LS,N2,N1)=CumSno2SnowLay(LS,N2,N1)+SnoXfer2SnoLay(LS,N2,N1)
        CumWat2SnowLay(LS,N2,N1)=CumWat2SnowLay(LS,N2,N1)+WatXfer2SnoLay(LS,N2,N1)
        CumIce2SnowLay(LS,N2,N1)=CumIce2SnowLay(LS,N2,N1)+IceXfer2SnoLay(LS,N2,N1)
        CumHeat2SnowLay(LS,N2,N1)=CumHeat2SnowLay(LS,N2,N1)+HeatXfer2SnoLay(LS,N2,N1)

        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_Xbndl_flx(NTG,LS,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_Xbndl_flx(NTN,LS,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO NTSA=idsalt_beg,idsalt_end
            trcSalt_TBLS(NTSA,LS,N2,N1)=trcSalt_TBLS(NTSA,LS,N2,N1)+trcSaltFlo2SnowLay(NTSA,LS,N2,N1)
          ENDDO

        ENDIF
      ENDIF
    ENDIF
  ENDDO D1205
  end subroutine SaltFromRunoffSnowpack
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
    CumSno2SnowLay(L,NY,NX)=0.0_r8
    CumWat2SnowLay(L,NY,NX)=0.0_r8
    CumIce2SnowLay(L,NY,NX)=0.0_r8
    CumHeat2SnowLay(L,NY,NX)=0.0_r8
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
!------------------------------------------------------------------------------------------
  subroutine MassFluxFromSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)

  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B
  integer :: NN,NTG,NTN

  D1202: DO NN=1,2
    !gaseous tracers
    DO NTG=idg_beg,idg_end-1
      trcg_TFloXSurRunoff(NTG,N2,N1)=trcg_TFloXSurRunoff(NTG,N2,N1)+trcg_2DFloXSurRunoff(NTG,N,NN,N2,N1)
    ENDDO

    !nutrient tracres
    DO NTN=ids_nut_beg,ids_nuts_end
      trcn_TFloXSurRunoff(NTN,N2,N1)=trcn_TFloXSurRunoff(NTN,N2,N1)+trcn_2DFloXSurRunoff(NTN,N,NN,N2,N1)
    ENDDO

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN    

      DO NTG=idg_beg,idg_end-1
        trcg_TFloXSurRunoff(NTG,N2,N1)=trcg_TFloXSurRunoff(NTG,N2,N1)-trcg_2DFloXSurRunoff(NTG,N,NN,N5,N4)
      ENDDO
      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_TFloXSurRunoff(NTN,N2,N1)=trcn_TFloXSurRunoff(NTN,N2,N1)-trcn_2DFloXSurRunoff(NTN,N,NN,N5,N4)
      ENDDO

    ENDIF 
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO NTG=idg_beg,idg_end-1
        trcg_TFloXSurRunoff(NTG,N2,N1)=trcg_TFloXSurRunoff(NTG,N2,N1)-trcg_2DFloXSurRunoff(NTG,N,NN,N5B,N4B)
      ENDDO
      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_TFloXSurRunoff(NTN,N2,N1)=trcn_TFloXSurRunoff(NTN,N2,N1)-trcn_2DFloXSurRunoff(NTN,N,NN,N5B,N4B)
      ENDDO

    ENDIF
  ENDDO D1202

  TDrysnoBySnowRedist(N2,N1)=TDrysnoBySnowRedist(N2,N1)+DrysnoBySnowRedistrib(N,N2,N1)&
    -DrysnoBySnowRedistrib(N,N5,N4)
  TWatBySnowRedist(N2,N1)=TWatBySnowRedist(N2,N1)+WatBySnowRedistrib(N,N2,N1)-WatBySnowRedistrib(N,N5,N4)
  TIceBySnowRedist(N2,N1)=TIceBySnowRedist(N2,N1)+IceBySnowRedistrib(N,N2,N1)-IceBySnowRedistrib(N,N5,N4)
  THeatBySnowRedist(N2,N1)=THeatBySnowRedist(N2,N1)+HeatBySnowRedistribution(N,N2,N1)-HeatBySnowRedistribution(N,N5,N4)
  !
  !     NET GAS AND SOLUTE FLUXES FROM RUNOFF AND SNOWPACK
  !
  !     T*QRS=net overland solute flux from runoff
  !     X*QRS=solute in runoff from TranspNoSalt.f
  !     T*QSS=net overland solute flux from snowpack
  !     X*QSS=solute in snowpack flux from TranspNoSalt.f
  !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
  !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
  !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
  !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
  !
  trcg_QSS(idg_CO2,N2,N1)=trcg_QSS(idg_CO2,N2,N1)+trcg_FloXSnow(idg_CO2,N,N2,N1)-trcg_FloXSnow(idg_CO2,N,N5,N4)
  trcg_QSS(idg_CH4,N2,N1)=trcg_QSS(idg_CH4,N2,N1)+trcg_FloXSnow(idg_CH4,N,N2,N1)-trcg_FloXSnow(idg_CH4,N,N5,N4)
  trcg_QSS(idg_O2,N2,N1)=trcg_QSS(idg_O2,N2,N1)+trcg_FloXSnow(idg_O2,N,N2,N1)-trcg_FloXSnow(idg_O2,N,N5,N4)
  trcg_QSS(idg_N2,N2,N1)=trcg_QSS(idg_N2,N2,N1)+trcg_FloXSnow(idg_N2,N,N2,N1)-trcg_FloXSnow(idg_N2,N,N5,N4)
  trcg_QSS(idg_N2O,N2,N1)=trcg_QSS(idg_N2O,N2,N1)+trcg_FloXSnow(idg_N2O,N,N2,N1)-trcg_FloXSnow(idg_N2O,N,N5,N4)
  trcg_QSS(idg_NH3,N2,N1)=trcg_QSS(idg_NH3,N2,N1)+trcg_FloXSnow(idg_NH3,N,N2,N1)-trcg_FloXSnow(idg_NH3,N,N5,N4)

  trcn_QSS(ids_NH4,N2,N1)=trcn_QSS(ids_NH4,N2,N1)+trcn_FloXSnow(ids_NH4,N,N2,N1)-trcn_FloXSnow(ids_NH4,N,N5,N4)
  trcn_QSS(ids_NO3,N2,N1)=trcn_QSS(ids_NO3,N2,N1)+trcn_FloXSnow(ids_NO3,N,N2,N1)-trcn_FloXSnow(ids_NO3,N,N5,N4)
  trcn_QSS(ids_H1PO4,N2,N1)=trcn_QSS(ids_H1PO4,N2,N1)+trcn_FloXSnow(ids_H1PO4,N,N2,N1)-trcn_FloXSnow(ids_H1PO4,N,N5,N4)
  trcn_QSS(ids_H2PO4,N2,N1)=trcn_QSS(ids_H2PO4,N2,N1)+trcn_FloXSnow(ids_H2PO4,N,N2,N1)-trcn_FloXSnow(ids_H2PO4,N,N5,N4)

  !     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
  call SaltThruFluxRunoffAndSnowpack(N,N1,N2,N4,N5,N4B,N5B)
  !

  end subroutine MassFluxFromSnowRunoff
!------------------------------------------------------------------------------------------

  subroutine ChemicalBySnowRedistribution(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: NTA,NTG,NTS
!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TDrysnoBySnowRedist(NY,NX))>0._r8)THEN
    DO NTG=idg_beg,idg_end-1
      trcg_solsml(NTG,1,NY,NX)=trcg_solsml(NTG,1,NY,NX)+trcg_QSS(NTG,NY,NX)
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trcn_solsml(NTS,1,NY,NX)=trcn_solsml(NTS,1,NY,NX)+trcn_QSS(NTS,NY,NX)
    ENDDO
    IF(salt_model)THEN
      DO NTA=idsalt_beg,idsalt_end
        trcs_solsml(NTA,1,NY,NX)=trcs_solsml(NTA,1,NY,NX)+trcSalt_TQS(NTA,NY,NX)
      ENDDO
    ENDIF
  ENDIF
  end subroutine ChemicalBySnowRedistribution
!------------------------------------------------------------------------------------------

  subroutine SaltThruFluxRunoffAndSnowpack(N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B

  integer :: NN,NTSA
!     begin_execution
!
!     TQR*=net overland solute flux in runoff
!     XQR*=solute in runoff from TranspSalt.f
!     TQS*=net overland solute flux in snow drift
!     XQS*=solute in snow drift from TranspSalt.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  IF(salt_model)THEN
    D1203: DO NN=1,2

      DO NTSA=idsalt_beg,idsalt_end
        trcSalt_TQR(NTSA,N2,N1)=trcSalt_TQR(NTSA,N2,N1)+trc_salt_rof_bounds(NTSA,N,NN,N2,N1)
      ENDDO

      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
! runoff direction
        DO NTSA=idsalt_beg,idsalt_end
          trcSalt_TQR(NTSA,N2,N1)=trcSalt_TQR(NTSA,N2,N1)-trc_salt_rof_bounds(NTSA,N,NN,N5,N4)
        ENDDO
      ENDIF

      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        DO NTSA=idsalt_beg,idsalt_end
          trcSalt_TQR(NTSA,N2,N1)=trcSalt_TQR(NTSA,N2,N1)-trc_salt_rof_bounds(NTSA,N,NN,N5B,N4B)
        ENDDO
      ENDIF
    ENDDO D1203

    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_TQS(NTSA,N2,N1)=trcSalt_TQS(NTSA,N2,N1)+trcSalt_XQS(NTSA,N,N2,N1)-trcSalt_XQS(NTSA,N,N5,N4)
    ENDDO
  ENDIF
  end subroutine SaltThruFluxRunoffAndSnowpack

!------------------------------------------------------------------------------------------

  subroutine OverlandSnowFlow(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: NTU,NTSA,NTG

    !    SOLUTES
!  exclude NH3B
    DO NTG=idg_beg,idg_end-1
      trc_solml_vr(NTG,0,NY,NX)=trc_solml_vr(NTG,0,NY,NX)+trcg_TFloXSurRunoff(NTG,NY,NX)
    ENDDO

  DO NTU=ids_nut_beg,ids_nuts_end
    trc_solml_vr(NTU,0,NY,NX)=trc_solml_vr(NTU,0,NY,NX)+trcn_TFloXSurRunoff(NTU,NY,NX)
  ENDDO

  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_solml(NTSA,0,NY,NX)=trcSalt_solml(NTSA,0,NY,NX)+trcSalt_TQR(NTSA,NY,NX)
    ENDDO
  ENDIF
  end subroutine OverlandSnowFlow    
end module SnowBalanceMod
