module SoilHydroParaMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  use SoilWaterDataType
  use SoilPropertyDataType
  use SoilHeatDataType
  use GridDataType
  use SoilPhysDataType
  use AqueChemDatatype
  use EcoSIMCtrlDataType
  use SOMDataType
  use FlagDataType
  use LandSurfDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use SurfLitterDataType
  use EcoSIMCtrlMod  
  use PhysPars  
  use EcoSiMParDataMod   , only : micpar
  use minimathmod  , only : test_aeqb,AZMAX1,AZMIN1
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__
  real(r8), parameter :: FORGW=0.25E+06_r8 !threshold for  C concentration in organic soil 	g Mg-1

  public :: GetSoilHydraulicVars
  public :: SoilHydroProperty
  public :: LitterHydroproperty
contains
!------------------------------------------------------------------------------------------

  subroutine GetSoilHydraulicVars(NY,NX)
  !
  !DESCRIPTIONS
  !compute hydraulic properties
  !called in hour1.F90 before doing hydrology 
  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8) :: FCX,FCLX
  real(r8) :: FCDX
  real(r8) :: PSDX
  real(r8) :: THETW1
  real(r8) :: WPX,WPLX
  integer :: K,L

  ! begin_execution
  DO L=NUI(NY,NX),NLI(NY,NX)
  ! WATER POTENTIALS
  !
  ! FC,WP=water contents at field capacity,wilting point,saturation
  ! PSISM,PSISE=matric,saturation water potential
  ! SRP=parameter for deviation from linear log-log water retention
  ! FC,WP=water contents at field capacity,wilting point from soil file
  ! FCL,WPL=log FC,WP
  ! FCD,PSD=FCL-WPL,log(POROS)-FCL
  ! FCI,WPI=FC,WP of ice
  ! THETIX=ice concentration
!
    IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX).AND.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AZMAX1(AMIN1(POROS(L,NY,NX),VOLW(L,NY,NX)/VOLY(L,NY,NX)),1.e-6_r8)
      IF(THETW1.LT.FC(L,NY,NX))THEN
        PSISM(L,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
          +((FCL(L,NY,NX)-LOG(THETW1))/FCD(L,NY,NX)*PSIMD(NY,NX))))
      ELSE IF(THETW1.LT.POROS(L,NY,NX)-DTHETW)THEN
        PSISM(L,NY,NX)=-EXP(PSIMS(NY,NX)+(((PSL(L,NY,NX)-LOG(THETW1)) &
          /PSD(L,NY,NX))**SRP(L,NY,NX)*PSISD(NY,NX)))        
      ELSE
        PSISM(L,NY,NX)=PSISE(L,NY,NX)
      ENDIF
    ELSE IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX).and.THETI(L,NY,NX)>ZEROS2(NY,NX))THEN
      FCX=FCI*THETI(L,NY,NX)
      WPX=WPI*THETI(L,NY,NX)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(L,NY,NX)-FCLX
      FCDX=FCLX-WPLX
      IF(THETW(L,NY,NX).LT.FCX)THEN
        PSISM(L,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
          +((FCLX-LOG(THETW(L,NY,NX)))/FCDX*PSIMD(NY,NX))))
      ELSE IF(THETW(L,NY,NX).LT.POROS(L,NY,NX)-DTHETW)THEN
        PSISM(L,NY,NX)=-EXP(PSIMS(NY,NX)+(((PSL(L,NY,NX)-LOG(THETW(L,NY,NX))) &
          /PSDX)*PSISD(NY,NX)))
      ELSE
        PSISM(L,NY,NX)=PSISE(L,NY,NX)
      ENDIF
    ELSE
      PSISM(L,NY,NX)=PSISE(L,NY,NX)
    ENDIF
!
!     SOIL OSMOTIC, GRAVIMETRIC AND MATRIC WATER POTENTIALS
!
!     PSISM,PSISO,PSISH,PSIST=matric,osmotic,gravimetric,total water potential
!
    PSISO(L,NY,NX)=-RGAS*1.E-6_r8*TKS(L,NY,NX)*CION(L,NY,NX)
    PSISH(L,NY,NX)=0.0098_r8*(ALT(NY,NX)-DPTH(L,NY,NX))
    PSIST(L,NY,NX)=AZMIN1(PSISM(L,NY,NX)+PSISO(L,NY,NX)+PSISH(L,NY,NX))

!
!     SOIL RESISTANCE TO ROOT PENETRATION
!
!     RSCS=soil resistance to root penetration (MPa)
!
!     IF(BKDS(L,NY,NX).GT.ZERO)THEN
!     CCLAYT=CCLAY(L,NY,NX)*1.0E+02
!     CORGCT=CORGC(L,NY,NX)*1.0E-04
!     CC=EXP(-3.6733-0.1447*CCLAYT+0.7653*CORGCT)
!     DD=-0.4805-0.1239*CCLAYT+0.2080*CORGCT
!     EE=3.8521+0.0963*CCLAYT
!     RSCS(L,NY,NX)=CC*THETW(L,NY,NX)**DD*BKDS(L,NY,NX)**EE
!     ELSE
    RSCS(L,NY,NX)=0.0_r8
!     ENDIF
!
!     SOIL HYDRAULIC CONDUCTIVITIES FROM AMBIENT SOIL WATER CONTENTS
!
!     CNDU=soil hydraulic conductivity for root uptake
!
    K=MAX(1,MIN(100,INT(100.0_r8*(POROS(L,NY,NX)-THETW(L,NY,NX))/POROS(L,NY,NX))+1))
    CNDU(L,NY,NX)=0.5_r8*(HCND(1,K,L,NY,NX)+HCND(3,K,L,NY,NX))
  END DO
  end subroutine GetSoilHydraulicVars


!------------------------------------------------------------------------------------------

  subroutine SoilHydroProperty(L,NY,NX,I,J)
  !
  !Set up soil hydraulic property
  implicit none
  integer, intent(in) :: L,NY,NX
  integer, intent(in) :: I,J
  real(r8) :: THETF
  real(r8) :: THETK(100),PSISK(0:100)

  integer :: K,M,N
  real(r8) :: XK,YK,SUM1,SUM2
  real(r8) :: VISCWL

  IF(CORGC(L,NY,NX).GT.FORGC)THEN
    SRP(L,NY,NX)=0.25_r8
  ELSE IF(CORGC(L,NY,NX).GT.0.5_r8*FORGC)THEN
    SRP(L,NY,NX)=0.33_r8
  ELSE
    SRP(L,NY,NX)=1.00_r8
  ENDIF

! double check cold_run() setup
  PSL(L,NY,NX)=LOG(POROS(L,NY,NX))
  IF((ISOIL(isoi_fc,L,NY,NX).EQ.0.AND.ISOIL(isoi_wp,L,NY,NX).EQ.0).OR.(.not.cold_run()))THEN
  ! read from check point file or if soil properties are set with soil file
    FCL(L,NY,NX)=LOG(FC(L,NY,NX))
    WPL(L,NY,NX)=LOG(WP(L,NY,NX))
    PSD(L,NY,NX)=PSL(L,NY,NX)-FCL(L,NY,NX)
    FCD(L,NY,NX)=FCL(L,NY,NX)-WPL(L,NY,NX)
  ELSE
    !
    !     DEFAULT SOIL HYDROLOGIC PPTYS (FIELD CAPACITY, WILTING POINT)
    !     IF ACTUAL VALUES WERE NOT INPUT TO THE SOIL FILE
    !
    !     THW,THI=initial soil water,ice content from soil file
    !
    IF(cold_run())THEN
    ! restart is defined as simulation starting from a previous run
      IF(ISOIL(isoi_fc,L,NY,NX).EQ.1.OR.ISOIL(isoi_wp,L,NY,NX).EQ.1)THEN
        !calculating FC or WP
        IF(CORGC(L,NY,NX).LT.FORGW)THEN
          FC(L,NY,NX)=0.2576_r8-0.20_r8*CSAND(L,NY,NX) &
                +0.36_r8*CCLAY(L,NY,NX)+0.60E-06*CORGC(L,NY,NX)
        ELSE
          IF(BKDS(L,NY,NX).LT.0.075_r8)THEN
            FC(L,NY,NX)=0.27_r8
          ELSEIF(BKDS(L,NY,NX).LT.0.195_r8)THEN
            FC(L,NY,NX)=0.62_r8
          ELSE
            FC(L,NY,NX)=0.71_r8
          ENDIF
        ENDIF
        FC(L,NY,NX)=FC(L,NY,NX)/(1.0_r8-FHOL(L,NY,NX))
        FC(L,NY,NX)=AMIN1(0.75_r8*POROS(L,NY,NX),FC(L,NY,NX))
        IF(CORGC(L,NY,NX).LT.FORGW)THEN
          WP(L,NY,NX)=0.0260_r8+0.50_r8*CCLAY(L,NY,NX)+0.32E-06_r8*CORGC(L,NY,NX)
        ELSE
          IF(BKDS(L,NY,NX).LT.0.075_r8)THEN
            WP(L,NY,NX)=0.04_r8
          ELSEIF(BKDS(L,NY,NX).LT.0.195_r8)THEN
            WP(L,NY,NX)=0.15_r8
          ELSE
            WP(L,NY,NX)=0.22_r8
          ENDIF
        ENDIF
        WP(L,NY,NX)=WP(L,NY,NX)/(1.0_r8-FHOL(L,NY,NX))
        WP(L,NY,NX)=AMIN1(0.75_r8*FC(L,NY,NX),WP(L,NY,NX))
      ENDIF
      FCL(L,NY,NX)=LOG(FC(L,NY,NX))
      WPL(L,NY,NX)=LOG(WP(L,NY,NX))
      PSD(L,NY,NX)=PSL(L,NY,NX)-FCL(L,NY,NX)
      FCD(L,NY,NX)=FCL(L,NY,NX)-WPL(L,NY,NX)
    ENDIF
!   IBEGIN:   start date of model run

    !IDATA(9): start year of model run
    IF(I.EQ.IBEGIN.AND.J.EQ.1.AND.is_first_year)THEN
      !first time step at the beginning year
      !THW=initial soil water content 
      !DPTH=depth to middle of soil layer [m]
      !DTBLZ=external water table depth, [m]
      IF(THW(L,NY,NX).GT.1.0_r8.OR.DPTH(L,NY,NX).GE.DTBLZ(NY,NX))THEN
        !below the water table, thus it is saturated
        THETW(L,NY,NX)=POROS(L,NY,NX)
      ELSEIF(test_aeqb(THW(L,NY,NX),1._r8))THEN
        !at field capacity
        THETW(L,NY,NX)=FC(L,NY,NX)
      ELSEIF(test_aeqb(THW(L,NY,NX),0._r8))THEN
        !at wilting point
        THETW(L,NY,NX)=WP(L,NY,NX)
      ELSEIF(THW(L,NY,NX).LT.0.0_r8)THEN
        !completely dry
        THETW(L,NY,NX)=0.0_r8
      ENDIF

      IF(THI(L,NY,NX).GT.1.0_r8.OR.DPTH(L,NY,NX).GE.DTBLZ(NY,NX))THEN
        THETI(L,NY,NX)=AZMAX1(AMIN1(POROS(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
      ELSEIF(test_aeqb(THI(L,NY,NX),1._r8))THEN
        THETI(L,NY,NX)=AZMAX1(AMIN1(FC(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
      ELSEIF(test_aeqb(THI(L,NY,NX),0._r8))THEN
        THETI(L,NY,NX)=AZMAX1(AMIN1(WP(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
      ELSEIF(THI(L,NY,NX).LT.0.0_r8)THEN
        THETI(L,NY,NX)=0.0_r8
      ENDIF
      
      IF(cold_run())THEN
      !in a cold run, set it 
        VOLW(L,NY,NX)=THETW(L,NY,NX)*VOLX(L,NY,NX)
        VOLWX(L,NY,NX)=VOLW(L,NY,NX)
        VOLWH(L,NY,NX)=THETW(L,NY,NX)*VOLAH(L,NY,NX)
        VOLI(L,NY,NX)=THETI(L,NY,NX)*VOLX(L,NY,NX)
        VOLIH(L,NY,NX)=THETI(L,NY,NX)*VOLAH(L,NY,NX)
        VHCP(L,NY,NX)=VHCM(L,NY,NX)+Cpw*(VOLW(L,NY,NX) &
          +VOLWH(L,NY,NX))+Cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
        THETWZ(L,NY,NX)=THETW(L,NY,NX)
        THETIZ(L,NY,NX)=THETI(L,NY,NX)
      ENDIF
    ENDIF
  ENDIF
  VOLP(L,NY,NX)=AZMAX1(VOLA(L,NY,NX)-VOLW(L,NY,NX)-VOLI(L,NY,NX)) &
    +AZMAX1(VOLAH(L,NY,NX)-VOLWH(L,NY,NX)-VOLIH(L,NY,NX))
  IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    THETP(L,NY,NX)=VOLP(L,NY,NX)/VOLY(L,NY,NX)
  ELSE
    THETP(L,NY,NX)=0.0_r8
  ENDIF
  IF(BKDS(L,NY,NX).GT.ZERO)THEN
    THETY(L,NY,NX)=EXP((PSIMX(NY,NX)-LOG(-PSIHY))*FCD(L,NY,NX)/PSIMD(NY,NX)+FCL(L,NY,NX))
  ELSE
    THETY(L,NY,NX)=ZERO2
  ENDIF
      !
      !     SATURATED HYDRAULIC CONDUCTIVITY FROM SWC AT SATURATION VS.
      !     -0.033 MPA (MINERAL SOILS) IF NOT ENTERED IN SOIL FILE IN 'READS'
      !
      !     SCNV,SCNH=vertical,lateral saturated hydraulic conductivity
  !
  IF(ISOIL(isoi_scnv,L,NY,NX).EQ.1)THEN
    !computing vertical saturated hydraulic conductivity
    IF(CORGC(L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS(L,NY,NX),EXP((PSIMS(NY,NX)-LOG(0.033_r8)) &
        *(PSL(L,NY,NX)-FCL(L,NY,NX))/PSISD(NY,NX)+PSL(L,NY,NX)))
      SCNV(L,NY,NX)=1.54_r8*((POROS(L,NY,NX)-THETF)/THETF)**2
    ELSE
      SCNV(L,NY,NX)=0.10_r8+75.0_r8*1.0E-15_r8**BKDS(L,NY,NX)
      SCNV(L,NY,NX)=SCNV(L,NY,NX)*FMPR(L,NY,NX)
    ENDIF
  ENDIF

  IF(ISOIL(isoi_scnh,L,NY,NX).EQ.1)THEN
    !computing horizontal saturated hydraulic conductivity
    IF(CORGC(L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS(L,NY,NX),EXP((PSIMS(NY,NX)-LOG(0.033_r8)) &
        *(PSL(L,NY,NX)-FCL(L,NY,NX))/PSISD(NY,NX)+PSL(L,NY,NX)))
      SCNH(L,NY,NX)=1.54_r8*((POROS(L,NY,NX)-THETF)/THETF)**2._r8
    ELSE
      SCNH(L,NY,NX)=0.10_r8+75.0_r8*1.0E-15_r8**BKDS(L,NY,NX)
      SCNH(L,NY,NX)=SCNH(L,NY,NX)*FMPR(L,NY,NX)
    ENDIF
  ENDIF

  !
  !     HYDRAULIC CONDUCTIVITY FUNCTION FROM KSAT AND SOIL WATER RELEASE CURVE
  !
  !     THETK,PSISK=micropore class water content,potential
  !     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
  !
  !     IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
  SUM2=0.0_r8
  DO  K=1,100
    XK=K-1
    THETK(K)=POROS(L,NY,NX)-(XK/100.0_r8*POROS(L,NY,NX))
    IF(THETK(K).LT.FC(L,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)+((FCL(L,NY,NX)-LOG(THETK(K))) &
        /FCD(L,NY,NX)*PSIMD(NY,NX))))
    ELSEIF(THETK(K).LT.POROS(L,NY,NX)-DTHETW)THEN
      PSISK(K)=-EXP(PSIMS(NY,NX)+(((PSL(L,NY,NX)-LOG(THETK(K))) &
        /PSD(L,NY,NX))**SRP(L,NY,NX)*PSISD(NY,NX)))
    ELSE
      PSISK(K)=PSISE(L,NY,NX)
    ENDIF
    SUM2=SUM2+(2*K-1)/(PSISK(K)**2_r8)
  ENDDO

  DO  K=1,100
    SUM1=0.0_r8
    XK=K-1
    YK=((100.0_r8-XK)/100.0_r8)**1.33_r8
    DO  M=K,100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2_r8)
    ENDDO
    DO  N=1,3
      IF(N.EQ.3)THEN
        !vertical
        HCND(N,K,L,NY,NX)=SCNV(L,NY,NX)*YK*SUM1/SUM2
        IF(K.GT.1.AND.PSISK(K).LT.PSISA(L,NY,NX).AND.PSISK(K-1).GE.PSISA(L,NY,NX))THEN
          THETS(L,NY,NX)=THETK(K)
        ENDIF
      ELSE
        !horizontal
        HCND(N,K,L,NY,NX)=SCNH(L,NY,NX)*YK*SUM1/SUM2
      ENDIF
    ENDDO
  ENDDO

!     SOIL MACROPORE DIMENSIONS AND CONDUCTIVITY FROM MACROPORE FRACTION
!     ENTERED IN 'READS'
!
!     PHOL,NHOL,HRAD=path length between, number,radius of macropores
!     CNDH=macropore hydraulic conductivity
!
  HRAD(L,NY,NX)=0.5E-03_r8
  NHOL(L,NY,NX)=INT(VOLAH(L,NY,NX)/(PICON*HRAD(L,NY,NX)**2._r8*VOLTI(L,NY,NX)))

  IF(NHOL(L,NY,NX).GT.0.0_r8)THEN
    PHOL(L,NY,NX)=1.0_r8/(SQRT(PICON*NHOL(L,NY,NX)))
  ELSE
    PHOL(L,NY,NX)=1.0_r8
  ENDIF
  VISCWL=VISCW*EXP(0.533_r8-0.0267_r8*TCS(L,NY,NX))
  CNDH(L,NY,NX)=3.6E+03_r8*PICON*NHOL(L,NY,NX)*HRAD(L,NY,NX)**4._r8/(8.0_r8*VISCWL)
  end subroutine SoilHydroProperty

!------------------------------------------------------------------------------------------

  subroutine LitterHydroproperty(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer  :: K,M
  real(r8) :: XK,YK
  real(r8) :: SUM1,SUM2
  integer, parameter :: n100=100
  real(r8) :: THETK(n100),PSISK(0:n100)

  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    BKDS(0,NY,NX)=BKVL(0,NY,NX)/VOLT(0,NY,NX)
  ELSE
    BKDS(0,NY,NX)=BKRS(micpar%k_fine_litr)
  ENDIF
  THETY(0,NY,NX)=EXP((PSIMX(NY,NX)-LOG(-PSIHY))*FCD(0,NY,NX)/PSIMD(NY,NX)+FCL(0,NY,NX))
  SUM2=0.0_r8
  D1220: DO  K=1,n100
    XK=K-1
    THETK(K)=POROS0(NY,NX)-(XK/n100*POROS0(NY,NX))
    IF(THETK(K).LT.FC(0,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)+((FCL(0,NY,NX)-LOG(THETK(K))) &
          /FCD(0,NY,NX)*PSIMD(NY,NX))))
    ELSEIF(THETK(K).LT.POROS0(NY,NX))THEN
      PSISK(K)=-EXP(PSIMS(NY,NX)+(((PSL(0,NY,NX)-LOG(THETK(K))) &
          /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
    ELSE
      PSISK(K)=PSISE(0,NY,NX)
    ENDIF
    SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
  ENDDO D1220

  D1235: DO  K=1,n100
    SUM1=0.0_r8
    XK=K-1
    YK=((n100-XK)/n100)**1.33_r8
    D1230: DO M=K,n100
        SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2._r8)
    ENDDO D1230
    HCND(3,K,0,NY,NX)=SCNV(0,NY,NX)*YK*SUM1/SUM2
    HCND(1,K,0,NY,NX)=0.0_r8
    HCND(2,K,0,NY,NX)=0.0_r8
    IF(K.GT.1.AND.PSISK(K).LT.PSISA(0,NY,NX).AND.PSISK(K-1).GE.PSISA(0,NY,NX))THEN
      THETS(0,NY,NX)=THETK(K)
    ENDIF
  ENDDO D1235
  end subroutine LitterHydroproperty

end module SoilHydroParaMod
