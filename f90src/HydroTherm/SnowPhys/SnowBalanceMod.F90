module SnowBalanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : endrun
  use minimathmod, only : AZMAX1
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

  character(len=*), parameter :: mod_filename = __FILE__

  public :: SnowMassUpdate
  public :: SnowpackLayering
  public :: ZeroSnowArrays
  public :: SaltFromRunoffSnowpack
  public :: OverlandSnowFlow
  public :: ChemicalBySnowRedistribution
  public :: FluxFromSnowRunoff
  contains

  subroutine SnowMassUpdate(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: VOLSWI,ENGYW
  integer :: L
      !
      !     CALCULATE SNOWPACK TEMPERATURE FROM ITS CHANGE
      !     IN HEAT STORAGE
      !
  VOLSS(NY,NX)=0.0_r8
  VOLWS(NY,NX)=0.0_r8
  VOLIS(NY,NX)=0.0_r8
  VOLS(NY,NX)=0.0_r8
  DPTHS(NY,NX)=0.0_r8
  VOLSWI=0.0_r8

  D9780: DO L=1,JS
    call UpdateSnowLayers(L,NY,NX,VOLSWI)

    call UpdateSoluteInSnow(L,NY,NX)
  ENDDO D9780
!
!     SNOW RUNOFF
!
  VOLSSL(1,NY,NX)=VOLSSL(1,NY,NX)+TQS(NY,NX)
  VOLWSL(1,NY,NX)=VOLWSL(1,NY,NX)+TQW(NY,NX)
  VOLISL(1,NY,NX)=VOLISL(1,NY,NX)+TQI(NY,NX)
  ENGYW=VHCPW(1,NY,NX)*TKW(1,NY,NX)
  VHCPW(1,NY,NX)=cps*VOLSSL(1,NY,NX)+cpw*VOLWSL(1,NY,NX)+cpi*VOLISL(1,NY,NX)
  IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
    TKW(1,NY,NX)=(ENGYW+THQS(NY,NX))/VHCPW(1,NY,NX)
  ELSE
    TKW(1,NY,NX)=TairK(NY,NX)
  ENDIF
!
! IF SNOWPACK DISAPPEARS
  call SnowpackDisapper(NY,NX)

  end subroutine SnowMassUpdate

!------------------------------------------------------------------------------------------

  subroutine SnowpackDisapper(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
  real(r8) :: tksx
  real(r8) :: FLWS,FLWW,FLWI
  real(r8) :: HFLWS,ENGYS,ENGY2,ENGY1,ENGY
!     begin_execution
!
  IF(VHCPW(1,NY,NX).GT.0.0_r8.AND.VHCPW(1,NY,NX).LE.VHCPWX(NY,NX) &
    .AND.TairK(NY,NX).GT.TFICE)THEN
    ENGYS=TKW(1,NY,NX)*VHCPW(1,NY,NX)
    ENGY1=TKS(NUM(NY,NX),NY,NX)*VHCP(NUM(NY,NX),NY,NX)
    FLWS=VOLSSL(1,NY,NX)
    FLWW=VOLWSL(1,NY,NX)
    FLWI=VOLISL(1,NY,NX)
    HFLWS=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TKW(1,NY,NX)
    VOLSSL(1,NY,NX)=0.0_r8
    VOLWSL(1,NY,NX)=0.0_r8
    VOLISL(1,NY,NX)=0.0_r8
    VHCPW(1,NY,NX)=0.0_r8

    VOLSS(NY,NX)=0.0_r8
    VOLWS(NY,NX)=0.0_r8
    VOLIS(NY,NX)=0.0_r8
    VOLS(NY,NX)=0.0_r8
    DPTHS(NY,NX)=0.0_r8
    D9770: DO L=1,JS
      DENSS(L,NY,NX)=DENS0(NY,NX)
    ENDDO D9770
    VOLW(0,NY,NX)=VOLW(0,NY,NX)+FLWW
    VOLI(0,NY,NX)=VOLI(0,NY,NX)+FLWI+FLWS/DENSI
    ENGY=VHCP(NUM(NY,NX),NY,NX)*TKS(NUM(NY,NX),NY,NX)
    VHCP(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX) &
      +cpw*(VOLW(NUM(NY,NX),NY,NX)+VOLWH(NUM(NY,NX),NY,NX)) &
      +cpi*(VOLI(NUM(NY,NX),NY,NX)+VOLIH(NUM(NY,NX),NY,NX))
    IF(VHCP(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      TKSX=TKS(NUM(NY,NX),NY,NX)
      TKS(NUM(NY,NX),NY,NX)=(ENGY+HFLWS)/VHCP(NUM(NY,NX),NY,NX)
      if(abs(TKS(NUM(NY,NX),NY,NX)/tksx-1._r8)>0.025_r8)then
        TKS(NUM(NY,NX),NY,NX)=TKSX
      endif
    ELSE
      TKS(NUM(NY,NX),NY,NX)=TairK(NY,NX)
    ENDIF
    ENGY2=VHCP(NUM(NY,NX),NY,NX)*TKS(NUM(NY,NX),NY,NX)

  ENDIF
  end subroutine SnowpackDisapper

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowLayers(L,NY,NX,VOLSWI)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(inout) :: VOLSWI
  real(r8) :: TKWX,VHCPWZ(JZ,JY,JX)
  real(r8) :: ENGYW
  real(r8) :: VOLSF,TCASF
  real(r8) :: DENSF,DDENS2,DDENS1
  real(r8) :: CVISC,DENSX
! begin_execution
!
! ADD CHANGES IN SNOW, WATER AND ICE
!
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
! XWFLXS,XWFLXI=freeze-thaw flux from watsub.f
!
!
  VOLSSL(L,NY,NX)=VOLSSL(L,NY,NX)+TFLWS(L,NY,NX)-XWFLXS(L,NY,NX)
  VOLWSL(L,NY,NX)=VOLWSL(L,NY,NX)+TFLWW(L,NY,NX)+XWFLXS(L,NY,NX)+XWFLXI(L,NY,NX)
  if(VOLWSL(L,NY,NX)<0._r8)then
    if(L>1)then
      VOLWSL(L-1,NY,NX)=VOLWSL(L-1,NY,NX)+VOLWSL(L,NY,NX)
    elseif(L==1)then
      VOLW(0,NY,NX)=VOLW(0,NY,NX)+VOLWSL(L,NY,NX)
    endif
    VOLWSL(L,NY,NX)=0._r8
  endif
  VOLISL(L,NY,NX)=VOLISL(L,NY,NX)+TFLWI(L,NY,NX)-XWFLXI(L,NY,NX)/DENSI
  if(VOLWSL(L,NY,NX)<0._r8)then
    if(L>1)then
      VOLISL(L-1,NY,NX)=VOLISL(L-1,NY,NX)+VOLISL(L,NY,NX)
    elseif(L==1)then
      VOLI(0,NY,NX)=VOLI(0,NY,NX)+VOLISL(L,NY,NX)
    endif
    VOLISL(L,NY,NX)=0._r8
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
    VOLSWI=VOLSWI+0.5_r8*(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
      +VOLISL(L,NY,NX)*DENSI)
    if(VOLSWI<0._r8)then
      write(*,*)'VOLSWI=',VOLSWI,VOLSSL(L,NY,NX),VOLWSL(L,NY,NX), &
        VOLISL(L,NY,NX)*DENSI
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
!
!   RESET SNOW SURFACE DENSITY FOR SNOWFALL
!
    IF(XFLWS(L,NY,NX).GT.0.0_r8)THEN
      DENSX=DENSS(L,NY,NX)
      TCASF=AMAX1(-15.0_r8,AMIN1(2.0_r8,TCA(NY,NX)))
      DENSF=0.05_r8+1.7E-03_r8*(TCASF+15.0_r8)**1.5_r8
      VOLSF=XFLWS(L,NY,NX)/DENSF+(VOLSSL(L,NY,NX)-XFLWS(L,NY,NX))/DENSS(L,NY,NX)
      DENSS(L,NY,NX)=VOLSSL(L,NY,NX)/VOLSF
      if(DENSS(L,NY,NX)<0._r8)write(*,*)'VOLSSL=',VOLSSL(L,NY,NX),VOLSF
    ENDIF
  ELSE
    VOLSWI=VOLSWI+0.5_r8*(VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX) &
      +VOLISL(L-1,NY,NX)*DENSI+VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
      +VOLISL(L,NY,NX)*DENSI)
    if(VOLSWI<0._r8)then
      write(*,*)'iVOLSWI=',VOLSWI,VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX) &
        +VOLISL(L-1,NY,NX)*DENSI,VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
        +VOLISL(L,NY,NX)*DENSI
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDIF
!
! SNOWPACK COMPRESSION
!
! DDENS1,DDENS2=temperature, compression effect on snow density
! DENSS=snow density in layer
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! VOLSL=snowpack layer volume
! DLYRS=snowpack layer depth
! CDPTHS=cumulative depth to bottom of snowpack layer
! VHCPW=snowpack layer heat capacity
! TKW,TCW=snowpack layer temperature K,oC
! THFLWW=convective heat fluxes of snow,water,ice in snowpack
! XTHAWW=latent heat flux from freeze-thaw from watsub.f
! HEATIN=cumulative net surface heat transfer
! VOLSS,VOLWS,VOLIS=total snow water equivalent, water, ice content of snowpack
! VOLS,DPTHS=total snowpack volume, depth
!
  IF(DENSS(L,NY,NX).LT.0.25_r8)THEN
    DDENS1=DENSS(L,NY,NX)*1.0E-05_r8*EXP(0.04_r8*TCW(L,NY,NX))
  ELSE
    DDENS1=0.0_r8
  ENDIF
  CVISC=0.25_r8*EXP(-0.08_r8*TCW(L,NY,NX)+23.0_r8*DENSS(L,NY,NX))

  !if(curday>=83)write(*,*)'CVISC=',CVISC,TCW(L,NY,NX),DENSS(L,NY,NX)
  DDENS2=DENSS(L,NY,NX)*VOLSWI/(AREA(3,NU(NY,NX),NY,NX)*CVISC)
  
  if(DDENS2<0._r8)write(*,*)'DDENS2=',DENSS(L,NY,NX),VOLSWI,CVISC
  DENSS(L,NY,NX)=DENSS(L,NY,NX)+DDENS1+DDENS2
  if(DENSS(L,NY,NX)<0._r8)write(*,*)'DDENS1=',DENSS(L,NY,NX),DDENS1,DDENS2

  IF(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLSL(L,NY,NX)=VOLSSL(L,NY,NX)/DENSS(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)
    DLYRS(L,NY,NX)=AZMAX1(VOLSL(L,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
    CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)+DLYRS(L,NY,NX)
    VHCPWZ(L,NY,NX)=VHCPW(L,NY,NX)
    TKWX=TKW(L,NY,NX)
    ENGYW=VHCPW(L,NY,NX)*TKW(L,NY,NX)
    VHCPW(L,NY,NX)=cps*VOLSSL(L,NY,NX)+cpw*VOLWSL(L,NY,NX)+cpi*VOLISL(L,NY,NX)
    IF(VHCPWZ(L,NY,NX).GT.VHCPWX(NY,NX).AND.VHCPW(L,NY,NX).GT.ZEROS(NY,NX))THEN
      TKW(L,NY,NX)=(ENGYW+THFLWW(L,NY,NX)+XTHAWW(L,NY,NX))/VHCPW(L,NY,NX)
    ELSE
      IF(L.EQ.1)THEN
        TKW(L,NY,NX)=TairK(NY,NX)
      ELSE
        TKW(L,NY,NX)=TKW(L-1,NY,NX)
      ENDIF
      IF(VHCPWZ(L,NY,NX).GT.ZEROS(NY,NX))THEN
        HEATIN=HEATIN+(TKW(L,NY,NX)-TKWX)*VHCPW(L,NY,NX)
      ENDIF
    ENDIF
    VOLSS(NY,NX)=VOLSS(NY,NX)+VOLSSL(L,NY,NX)
    VOLWS(NY,NX)=VOLWS(NY,NX)+VOLWSL(L,NY,NX)
    VOLIS(NY,NX)=VOLIS(NY,NX)+VOLISL(L,NY,NX)
    VOLS(NY,NX)=VOLS(NY,NX)+VOLSL(L,NY,NX)
    DPTHS(NY,NX)=DPTHS(NY,NX)+DLYRS(L,NY,NX)
  ELSE
    VOLSSL(L,NY,NX)=0.0_r8
    VOLWSL(L,NY,NX)=0.0_r8
    VOLISL(L,NY,NX)=0.0_r8
    VOLSL(L,NY,NX)=0.0_r8
    DLYRS(L,NY,NX)=0.0_r8
    CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)
    VHCPW(L,NY,NX)=0.0_r8
    IF(L.EQ.1)THEN
      TKW(L,NY,NX)=TairK(NY,NX)
    ELSE
      TKW(L,NY,NX)=TKW(L-1,NY,NX)
    ENDIF
  ENDIF
  TCW(L,NY,NX)=units%Kelvin2Celcius(TKW(L,NY,NX))
  end subroutine UpdateSnowLayers
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
  !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
  !
  IF(salt_model)THEN
    DO NTSA=idsa_beg,idsa_end
      trcs_solsml(NTSA,L,NY,NX)=trcs_solsml(NTSA,L,NY,NX)+trcsa_TBLS(NTSA,L,NY,NX)
    ENDDO

  ENDIF

  end subroutine UpdateSoluteInSnow

!------------------------------------------------------------------------------------------

  subroutine SnowpackLayering(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: FX,FY
  integer :: L,L1,L0
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
  real(r8) :: DDLYXS,DDLYRS
  real(r8) :: DDLYXX,VOLSLX
  integer :: IFLGLS,NTN,NTG,NTU,NTSA

!     begin_execution
! from surface to bottom, and modify the bottom layer
  IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
    D325: DO L=1,JS-1
!      VOLSLX=VOLSL(L,NY,NX)
      IF(VOLSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        !compute the difference in thickness compared to the initial
        DDLYXS=(VOLSI(L,NY,NX)-VOLSSL(L,NY,NX)/DENSS(L,NY,NX) &
          -VOLWSL(L,NY,NX)-VOLISL(L,NY,NX))/AREA(3,L,NY,NX)
!        DDLYXX=DDLYXS
        IF(DDLYXS.LT.-ZERO.OR.DLYRS(L+1,NY,NX).GT.ZERO)THEN
        ! volume is greater than allowed, or next layer exists
          DDLYRS=AMIN1(DDLYXS,DLYRS(L+1,NY,NX))
          IFLGLS=1         !expand
        ELSE
          !volume less than allowed, and no next layer
          DDLYXS=(VOLSL(L,NY,NX)-VOLSSL(L,NY,NX)/DENSS(L,NY,NX) &
            -VOLWSL(L,NY,NX)-VOLISL(L,NY,NX))/AREA(3,L,NY,NX)
          DDLYRS=DDLYXS
          IFLGLS=2         !shrink
        ENDIF
      ELSE
        DDLYRS=0.0_r8      !no change
        IFLGLS=0
      ENDIF
      !
      !     RESET SNOW LAYER DEPTHS
      !
      CDPTHS(L,NY,NX)=CDPTHS(L,NY,NX)+DDLYRS
      DLYRS(L,NY,NX)=CDPTHS(L,NY,NX)-CDPTHS(L-1,NY,NX)
!
      !     TRANSFER STATE VARIABLES BETWEEN LAYERS
      !
      IF(ABS(DDLYRS).GT.ZERO)THEN
        IF(DDLYRS.GT.0.0_r8)THEN
          L1=L
          L0=L+1
          IF(DDLYRS.LT.DDLYXS)THEN
            FX=1.0_r8
          ELSE
            FX=AMIN1(1.0_r8,DDLYRS*AREA(3,L0,NY,NX)/VOLSL(L0,NY,NX))
          ENDIF
        ELSE
          L1=L+1
          L0=L
          IF(VOLSL(L0,NY,NX).LT.VOLSI(L0,NY,NX))THEN
            FX=0.0_r8
          ELSE
            FX=AMIN1(1.0_r8,-DDLYRS*AREA(3,L0,NY,NX)/VOLSL(L0,NY,NX))
          ENDIF
        ENDIF
        IF(FX.GT.0.0_r8)THEN
          FY=1.0_r8-FX

!
!     TARGET SNOW LAYER
!
          VOLSSL(L1,NY,NX)=VOLSSL(L1,NY,NX)+FX*VOLSSL(L0,NY,NX)
          VOLWSL(L1,NY,NX)=VOLWSL(L1,NY,NX)+FX*VOLWSL(L0,NY,NX)
          VOLISL(L1,NY,NX)=VOLISL(L1,NY,NX)+FX*VOLISL(L0,NY,NX)
          VOLSL(L1,NY,NX)=VOLSSL(L1,NY,NX)/DENSS(L1,NY,NX)+VOLWSL(L1,NY,NX)+VOLISL(L1,NY,NX)
          ENGY1X=VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
          ENGY0X=VHCPW(L0,NY,NX)*TKW(L0,NY,NX)
          ENGY1=ENGY1X+FX*ENGY0X
          VHCPW(L1,NY,NX)=cps*VOLSSL(L1,NY,NX)+cpw*VOLWSL(L1,NY,NX)+cpi*VOLISL(L1,NY,NX)
          IF(VHCPW(L1,NY,NX).GT.ZEROS(NY,NX))THEN
            TKW(L1,NY,NX)=ENGY1/VHCPW(L1,NY,NX)
          ELSE
            TKW(L1,NY,NX)=TKW(L0,NY,NX)
          ENDIF
          TCW(L1,NY,NX)=units%Kelvin2Celcius(TKW(L1,NY,NX))

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
            DO NTSA=idsa_beg,idsa_end
              trcs_solsml(NTSA,L1,NY,NX)=trcs_solsml(NTSA,L1,NY,NX)+FX*trcs_solsml(NTSA,L0,NY,NX)
            ENDDO
          ENDIF
!
!     SOURCE SNOW LAYER
!
          VOLSSL(L0,NY,NX)=FY*VOLSSL(L0,NY,NX)
          VOLWSL(L0,NY,NX)=FY*VOLWSL(L0,NY,NX)
          VOLISL(L0,NY,NX)=FY*VOLISL(L0,NY,NX)
          VOLSL(L0,NY,NX)=VOLSSL(L0,NY,NX)/DENSS(L0,NY,NX)+VOLWSL(L0,NY,NX)+VOLISL(L0,NY,NX)
          ENGY0=FY*ENGY0X
          VHCPW(L0,NY,NX)=cps*VOLSSL(L0,NY,NX)+cpw*VOLWSL(L0,NY,NX)+cpi*VOLISL(L0,NY,NX)
          IF(VHCPW(L0,NY,NX).GT.ZEROS(NY,NX))THEN
            TKW(L0,NY,NX)=ENGY0/VHCPW(L0,NY,NX)
          ELSE
            TKW(L0,NY,NX)=TKW(L1,NY,NX)
          ENDIF
          TCW(L0,NY,NX)=units%Kelvin2Celcius(TKW(L0,NY,NX))

          DO NTG=idg_beg,idg_end-1
            trcg_solsml(NTG,L0,NY,NX)=FY*trcg_solsml(NTG,L0,NY,NX)
          ENDDO
          DO NTU=ids_nut_beg,ids_nuts_end
            trcn_solsml(NTU,L0,NY,NX)=FY*trcn_solsml(NTU,L0,NY,NX)
          ENDDO

          IF(salt_model)THEN
            DO NTSA=idsa_beg,idsa_end
              trcs_solsml(NTSA,L0,NY,NX)=FY*trcs_solsml(NTSA,L0,NY,NX)
            ENDDO

          ENDIF
!     IF(VOLWSL(L0,NY,NX)+VOLISL(L0,NY,NX)
!    2+VOLSSL(L0,NY,NX).LE.ZEROS(NY,NX))THEN
!     CDPTHS(L1,NY,NX)=CDPTHS(L0,NY,NX)
!     ENDIF
        ENDIF
      ENDIF
    ENDDO D325
  ENDIF
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
!     VHCPW,VHCPWX=current, minimum snowpack heat capacities
!     TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
!     THFLWW=convective heat fluxes of snow,water,ice in snowpack
!     XFLWS,XFLWW,XFLWI=snow,water,ice transfer from watsub.f
!     XHFLWW=convective heat flux from snow,water,ice transfer from watsub.f
!     FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil macropore,micropore,litter
!     HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter

  D1205: DO LS=1,JS
    IF(VHCPW(LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS.AND.VHCPW(LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
        TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)-XFLWS(LS2,N2,N1)
        TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1)-XFLWW(LS2,N2,N1) &
          -FLSWR(LS,N2,N1)-FLSW(LS,N2,N1)-FLSWH(LS,N2,N1)
        TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)-XFLWI(LS2,N2,N1)
        THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1) &
          -XHFLWW(LS2,N2,N1)-HFLSWR(LS,N2,N1)-HFLSW(LS,N2,N1)

!     NET SOLUTE FLUXES THROUGH SNOWPACK
!
!     T*BLS=net solute flux in snowpack
!     X*BLS=solute flux in snowpack from trnsfr.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_XBLS(NTG,LS,N2,N1) &
            -trcg_XBLS(NTG,LS2,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_XBLS(NTN,LS,N2,N1) &
            -trcn_XBLS(NTN,LS2,N2,N1)
        ENDDO
!
!     NET SALT FLUXES THROUGH SNOWPACK
!
!     T*BLS=net solute flux in snowpack
!     X*BLS=solute flux in snowpack from trnsfrs.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
        IF(salt_model)THEN
          DO NTSA=idsa_beg,idsa_end
            trcsa_TBLS(NTSA,LS,N2,N1)=trcsa_TBLS(NTSA,LS,N2,N1)+trcsa_XBLS(NTSA,LS,N2,N1)-trcsa_XBLS(NTSA,LS2,N2,N1)
          ENDDO
        ENDIF
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
      ELSE
        TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
        TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1) &
          -FLSWR(LS,N2,N1)-FLSW(LS,N2,N1)-FLSWH(LS,N2,N1)
        TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
        THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1) &
          -HFLSWR(LS,N2,N1)-HFLSW(LS,N2,N1)
! and NH3B
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_XBLS(NTG,LS,N2,N1) &
            -trcs_XFLS(NTG,3,0,N2,N1)-trcs_XFLS(NTG,3,NUM(N2,N1),N2,N1) &
            -trcs_XFHS(NTG,3,NUM(N2,N1),N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_XBLS(NTN,LS,N2,N1) &
            -trcs_XFLS(NTN,3,0,N2,N1)-trcs_XFLS(NTN,3,NUM(N2,N1),N2,N1) &
            -trcs_XFHS(NTN,3,NUM(N2,N1),N2,N1)
        ENDDO

!add band flux
        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1) &
          -trcs_XFLS(idg_NH3B,3,NUM(N2,N1),N2,N1)-trcs_XFHS(idg_NH3B,3,NUM(N2,N1),N2,N1)

        DO NTS=0,ids_nuts
          trcn_TBLS(ids_NH4+NTS,LS,N2,N1)=trcn_TBLS(ids_NH4+NTS,LS,N2,N1) &
            -trcs_XFLS(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)-trcs_XFHS(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO NTSA=idsa_beg,idsa_end
            trcsa_TBLS(NTSA,LS,NY,NX)=trcsa_TBLS(NTSA,LS,NY,NX)+trcsa_XBLS(NTSA,LS,NY,NX) &
              -trcsa_XFLS(NTSA,3,0,N2,N1)-trcsa_XFLS(NTSA,3,NUM(N2,N1),N2,N1) &
              -trcsa_XFHS(NTSA,3,NUM(N2,N1),N2,N1)
          ENDDO

!add band flux
          DO NTSA=0,idsa_nuts
            trcsa_TBLS(idsa_H0PO4+NTSA,LS,NY,NX)=trcsa_TBLS(idsa_H0PO4+NTSA,LS,NY,NX) &
              -trcsa_XFLS(idsa_H0PO4B+NTSA,3,NUM(N2,N1),N2,N1) &
              -trcsa_XFHS(idsa_H0PO4B+NTSA,3,NUM(N2,N1),N2,N1)
          ENDDO
        ENDIF
      ENDIF
!
!     WATER,GAS,SOLUTE,SALT FLUXES INTO SNOWPACK SURFACE
!
    ELSEIF(LS.EQ.1)THEN
      IF(abs(XFLWS(LS,N2,N1))>0._r8)THEN
        TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
        TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1)
        TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
        THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1)

        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_XBLS(NTG,LS,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_XBLS(NTN,LS,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO NTSA=idsa_beg,idsa_end
            trcsa_TBLS(NTSA,LS,N2,N1)=trcsa_TBLS(NTSA,LS,N2,N1)+trcsa_XBLS(NTSA,LS,N2,N1)
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

  TQS(NY,NX)=0.0_r8
  TQW(NY,NX)=0.0_r8
  TQI(NY,NX)=0.0_r8
  THQS(NY,NX)=0.0_r8

  trcn_TQR(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8
  trcg_QSS(idg_beg:idg_end-1,NY,NX)=0.0_r8  
  trcn_QSS(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8
  trcg_TQR(idg_beg:idg_end-1,NY,NX)=0.0_r8

  DO  L=1,JS
    trcg_TBLS(idg_beg:idg_end-1,L,NY,NX)=0.0_r8
    trcn_TBLS(ids_nut_beg:ids_nuts_end,L,NY,NX)=0.0_r8
    TFLWS(L,NY,NX)=0.0_r8
    TFLWW(L,NY,NX)=0.0_r8
    TFLWI(L,NY,NX)=0.0_r8
    THFLWW(L,NY,NX)=0.0_r8
  ENDDO

  IF(salt_model)THEN
!     INITIALIZE NET SOLUTE AND GAS FLUXES FROM SNOWPACK DRIFT
!
    trcsa_TQR(idsa_beg:idsa_end,NY,NX)=0.0_r8
    trcsa_TQS(idsa_beg:idsa_end,NY,NX)=0.0_r8
    DO  L=1,JS
      trcsa_TBLS(idsa_beg:idsa_end,L,NY,NX)=0.0_r8
    ENDDO
  endif
  end subroutine ZeroSnowArrays
!------------------------------------------------------------------------------------------
  subroutine FluxFromSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)

  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B
  integer :: NN,NTG,NTN

  D1202: DO NN=1,2
    !gaseous tracers
    DO NTG=idg_beg,idg_end-1
      trcg_TQR(NTG,N2,N1)=trcg_TQR(NTG,N2,N1)+trcg_XRS(NTG,N,NN,N2,N1)
    ENDDO

    !nutrient tracres
    DO NTN=ids_nut_beg,ids_nuts_end
      trcn_TQR(NTN,N2,N1)=trcn_TQR(NTN,N2,N1)+trcn_XRS(NTN,N,NN,N2,N1)
    ENDDO

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN    

      DO NTG=idg_beg,idg_end-1
        trcg_TQR(NTG,N2,N1)=trcg_TQR(NTG,N2,N1)-trcg_XRS(NTG,N,NN,N5,N4)
      ENDDO
      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_TQR(NTN,N2,N1)=trcn_TQR(NTN,N2,N1)-trcn_XRS(NTN,N,NN,N5,N4)
      ENDDO

    ENDIF 
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO NTG=idg_beg,idg_end-1
        trcg_TQR(NTG,N2,N1)=trcg_TQR(NTG,N2,N1)-trcg_XRS(NTG,N,NN,N5B,N4B)
      ENDDO
      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_TQR(NTN,N2,N1)=trcn_TQR(NTN,N2,N1)-trcn_XRS(NTN,N,NN,N5B,N4B)
      ENDDO

    ENDIF
  ENDDO D1202

  TQS(N2,N1)=TQS(N2,N1)+QS(N,N2,N1)-QS(N,N5,N4)
  TQW(N2,N1)=TQW(N2,N1)+QW(N,N2,N1)-QW(N,N5,N4)
  TQI(N2,N1)=TQI(N2,N1)+QI(N,N2,N1)-QI(N,N5,N4)
  THQS(N2,N1)=THQS(N2,N1)+HQS(N,N2,N1)-HQS(N,N5,N4)
  !
  !     NET GAS AND SOLUTE FLUXES FROM RUNOFF AND SNOWPACK
  !
  !     T*QRS=net overland solute flux from runoff
  !     X*QRS=solute in runoff from trnsfr.f
  !     T*QSS=net overland solute flux from snowpack
  !     X*QSS=solute in snowpack flux from trnsfr.f
  !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
  !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
  !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
  !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
  !
  trcg_QSS(idg_CO2,N2,N1)=trcg_QSS(idg_CO2,N2,N1)+XCOQSS(N,N2,N1)-XCOQSS(N,N5,N4)
  trcg_QSS(idg_CH4,N2,N1)=trcg_QSS(idg_CH4,N2,N1)+XCHQSS(N,N2,N1)-XCHQSS(N,N5,N4)
  trcg_QSS(idg_O2,N2,N1)=trcg_QSS(idg_O2,N2,N1)+XOXQSS(N,N2,N1)-XOXQSS(N,N5,N4)
  trcg_QSS(idg_N2,N2,N1)=trcg_QSS(idg_N2,N2,N1)+XNGQSS(N,N2,N1)-XNGQSS(N,N5,N4)
  trcg_QSS(idg_N2O,N2,N1)=trcg_QSS(idg_N2O,N2,N1)+XN2QSS(N,N2,N1)-XN2QSS(N,N5,N4)
  trcg_QSS(idg_NH3,N2,N1)=trcg_QSS(idg_NH3,N2,N1)+XN3QSS(N,N2,N1)-XN3QSS(N,N5,N4)

  trcn_QSS(ids_NH4,N2,N1)=trcn_QSS(ids_NH4,N2,N1)+XN4QSS(N,N2,N1)-XN4QSS(N,N5,N4)
  trcn_QSS(ids_NO3,N2,N1)=trcn_QSS(ids_NO3,N2,N1)+XNOQSS(N,N2,N1)-XNOQSS(N,N5,N4)
  trcn_QSS(ids_H1PO4,N2,N1)=trcn_QSS(ids_H1PO4,N2,N1)+XP1QSS(N,N2,N1)-XP1QSS(N,N5,N4)
  trcn_QSS(ids_H2PO4,N2,N1)=trcn_QSS(ids_H2PO4,N2,N1)+XP4QSS(N,N2,N1)-XP4QSS(N,N5,N4)

  !     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
  call SaltThruFluxRunoffAndSnowpack(N,N1,N2,N4,N5,N4B,N5B)
  !

  end subroutine FluxFromSnowRunoff
!------------------------------------------------------------------------------------------

  subroutine ChemicalBySnowRedistribution(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: NTA,NTG,NTS
!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TQS(NY,NX))>0._r8)THEN
    DO NTG=idg_beg,idg_end-1
      trcg_solsml(NTG,1,NY,NX)=trcg_solsml(NTG,1,NY,NX)+trcg_QSS(NTG,NY,NX)
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trcn_solsml(NTS,1,NY,NX)=trcn_solsml(NTS,1,NY,NX)+trcn_QSS(NTS,NY,NX)
    ENDDO
    IF(salt_model)THEN
      DO NTA=idsa_beg,idsa_end
        trcs_solsml(NTA,1,NY,NX)=trcs_solsml(NTA,1,NY,NX)+trcsa_TQS(NTA,NY,NX)
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
!     XQR*=solute in runoff from trnsfrs.f
!     TQS*=net overland solute flux in snow drift
!     XQS*=solute in snow drift from trnsfrs.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  IF(salt_model)THEN
    D1203: DO NN=1,2

      DO NTSA=idsa_beg,idsa_end
        trcsa_TQR(NTSA,N2,N1)=trcsa_TQR(NTSA,N2,N1)+trcsa_XQR(NTSA,N,NN,N2,N1)
      ENDDO

      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
! runoff direction
        DO NTSA=idsa_beg,idsa_end
          trcsa_TQR(NTSA,N2,N1)=trcsa_TQR(NTSA,N2,N1)-trcsa_XQR(NTSA,N,NN,N5,N4)
        ENDDO
      ENDIF

      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        DO NTSA=idsa_beg,idsa_end
          trcsa_TQR(NTSA,N2,N1)=trcsa_TQR(NTSA,N2,N1)-trcsa_XQR(NTSA,N,NN,N5B,N4B)
        ENDDO
      ENDIF
    ENDDO D1203

    DO NTSA=idsa_beg,idsa_end
      trcsa_TQS(NTSA,N2,N1)=trcsa_TQS(NTSA,N2,N1)+trcsa_XQS(NTSA,N,N2,N1)-trcsa_XQS(NTSA,N,N5,N4)
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
      trc_solml(NTG,0,NY,NX)=trc_solml(NTG,0,NY,NX)+trcg_TQR(NTG,NY,NX)
    ENDDO

  DO NTU=ids_nut_beg,ids_nuts_end
    trc_solml(NTU,0,NY,NX)=trc_solml(NTU,0,NY,NX)+trcn_TQR(NTU,NY,NX)
  ENDDO

  IF(salt_model)THEN
    DO NTSA=idsa_beg,idsa_end
      trcsa_solml(NTSA,0,NY,NX)=trcsa_solml(NTSA,0,NY,NX)+trcsa_TQR(NTSA,NY,NX)
    ENDDO
  ENDIF
  end subroutine OverlandSnowFlow    
end module SnowBalanceMod
