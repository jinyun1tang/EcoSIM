module SoilPhysParaMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  USE SoilWaterDataType
  use SoilPropertyDataType
  USE SoilPhysDataType
  use EcoSIMCtrlDataType
  use HydroThermData
  use PhysPars
  use EcoSiMParDataMod, only : micpar
  use EcoSimConst
  use FlagDataType
  USE AqueChemDatatype
  USE SOMDataType
  use SoilBGCDataType
  use MiniMathMod
  use SoilHeatDataType
  USE GridDataType
implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  PUBLIC :: CalcSoilWatPotential
  public :: SetDeepSoil
  public :: CalcSoilThermConductivity

  contains
!------------------------------------------------------------------------------------------
  subroutine CalcSoilWatPotential(NY,NX,N1,N2,N3,PSISoilMatric,THETA1S)
  implicit none
  integer, intent(in) :: NY,NX,N1,N2,N3
  real(r8), intent(out) :: PSISoilMatric    !< 0
  real(r8), optional, intent(out) :: THETA1S
  real(r8) :: FCDX,FCX,FCLX,WPLX,PSDX,WPX,THETA1
  !     BKVL=soil mass
  !     FC,WP=water contents at field capacity,wilting point
  !     FCL,WPL=log FC,WP
  !     FCD,PSD=FCL-WPL,log(POROS)-FCL
  !     PSISoilMatric,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
  !     PSIMX,PSIMD,LOGPSIAtSat=log water potential at FC,WP,saturation
  !     PSISD=PSIMX-LOGPSIAtSat
  !     SRP=parameter for deviation from linear log-log water retention
  !     PSISO=osmotic potential
  !     DTHETW=minimum water content for numerical purpose
  ! soil matric potential upper layer

  THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1),safe_adb(VLWatMicP1(N3,N2,N1),VLSoilMicP(N3,N2,N1))))

  IF(SoilMicPMassLayer(N3,N2,N1).GT.ZEROS(NY,NX))THEN
    !source layer is active soil
    IF(THETA1.LT.FieldCapacity(N3,N2,N1))THEN
      !water less than field capacity
      !PSIHY is the minimum water potential allowed, where hygroscopic water is < 0
      !held tightly on the surfaces of soil particles and exists as a thin layer of vapor
      PSISoilMatric=AMAX1(PSIHY,-EXP(LOGPSIFLD(N2,N1)+((LOGFldCapacity(N3,N2,N1)-LOG(THETA1))/FCD(N3,N2,N1)*LOGPSIMND(N2,N1))))
    ELSEIF(THETA1.LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISoilMatric=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS(N3,N2,N1)-LOG(THETA1))/PSD(N3,N2,N1))**SRP(N3,N2,N1)*LOGPSIMXD(N2,N1)))
    ELSE
      THETA1=POROS(N3,N2,N1)
      PSISoilMatric=PSISE(N3,N2,N1)
    ENDIF
    !
    !     SUBSURFCE UPPER WATER LAYER
    !
    !     FracSoiPAsIce,FracSoiPAsWat=ice,water concentration
    !     FCI,WPI=ice field capacity,wilting point
    !     PSISoilMatric=matric water potential
    !
  ELSEIF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1).and.FracSoiPAsIce(N3,N2,N1)>ZEROS2(N2,N1))THEN
    FCX=FCI*FracSoiPAsIce(N3,N2,N1)
    WPX=WPI*FracSoiPAsIce(N3,N2,N1)
    FCLX=LOG(FCX)
    WPLX=LOG(WPX)
    PSDX=LOGPOROS(N3,N2,N1)-FCLX
    FCDX=FCLX-WPLX
    IF(FracSoiPAsWat(N3,N2,N1).LT.FCX)THEN
      PSISoilMatric=AMAX1(PSIHY,-EXP(LOGPSIFLD(N2,N1)+((FCLX-LOG(FracSoiPAsWat(N3,N2,N1)))/FCDX*LOGPSIMND(NY,NX))))
    ELSEIF(FracSoiPAsWat(N3,N2,N1).LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISoilMatric=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS(N3,N2,N1)-LOG(FracSoiPAsWat(N3,N2,N1)))/PSDX)*LOGPSIMXD(N2,N1)))
    ELSE
      !saturated
      THETA1=POROS(N3,N2,N1)
      PSISoilMatric=PSISE(N3,N2,N1)
    ENDIF
  ELSE
    THETA1=POROS(N3,N2,N1)
    PSISoilMatric=PSISE(N3,N2,N1)
  ENDIF

  if(present(THETA1S))THETA1S=THETA1
  end subroutine CalcSoilWatPotential
!------------------------------------------------------------------------------------------


  subroutine SetDeepSoil(NY,NX,NM,JZ)

  implicit none
  integer, intent(in) :: NY,NX,NM,JZ
  integer :: L
  associate(                            &
  k_woody_litr => micpar%k_woody_litr , &
  k_fine_litr  => micpar%k_fine_litr  , &
  k_manure     => micpar%k_manure       &
  )

  write(*,*) "In SetDeepSoil, starting loop..."
  DO L=NM+1,JZ
    write(*,*) "On loop: ", L
    write(*,*) "Cumulative depth"
    CumDepth2LayerBottom(L,NY,NX)=2.0_r8*CumDepth2LayerBottom(L-1,NY,NX)-1.0_r8*CumDepth2LayerBottom(L-2,NY,NX)
    SoiBulkDensityt0(L,NY,NX)=SoiBulkDensityt0(L-1,NY,NX)
    FieldCapacity(L,NY,NX)=FieldCapacity(L-1,NY,NX)
    WiltPoint(L,NY,NX)=WiltPoint(L-1,NY,NX)
    SatHydroCondVert(L,NY,NX)=SatHydroCondVert(L-1,NY,NX)
    SatHydroCondHrzn(L,NY,NX)=SatHydroCondHrzn(L-1,NY,NX)
    CSAND(L,NY,NX)=CSAND(L-1,NY,NX)
    CSILT(L,NY,NX)=CSILT(L-1,NY,NX)
    CCLAY(L,NY,NX)=CCLAY(L-1,NY,NX)
    SoilFracAsMacP(L,NY,NX)=SoilFracAsMacP(L-1,NY,NX)
    ROCK(L,NY,NX)=ROCK(L-1,NY,NX)
    PH(L,NY,NX)=PH(L-1,NY,NX)
    CEC(L,NY,NX)=CEC(L-1,NY,NX)
    AEC(L,NY,NX)=AEC(L-1,NY,NX)

!     IF(IDWaterTable(NY,NX).EQ.0)THEN
!       0.25_r8 is the geometric decreasing ratio (tunable)
!  or different scheme can be used
    CORGC(L,NY,NX)=0.25_r8*CORGC(L-1,NY,NX)
    CORGR(L,NY,NX)=0.25_r8*CORGR(L-1,NY,NX)
    CORGN(L,NY,NX)=0.25_r8*CORGN(L-1,NY,NX)
    CORGP(L,NY,NX)=0.25_r8*CORGP(L-1,NY,NX)
!     ELSE
!       CORGC(L,NY,NX)=CORGC(L-1,NY,NX)
!       CORGR(L,NY,NX)=CORGR(L-1,NY,NX)
!       CORGN(L,NY,NX)=CORGN(L-1,NY,NX)
!       CORGP(L,NY,NX)=CORGP(L-1,NY,NX)
!     ENDIF
    write(*,*) "atm abundances"
    CNH4(L,NY,NX)=CNH4(L-1,NY,NX)
    write(*,*) "1"
    CNO3(L,NY,NX)=CNO3(L-1,NY,NX)
    write(*,*) "2"
    CPO4(L,NY,NX)=CPO4(L-1,NY,NX)
    write(*,*) "3"
    !9/12 - error seems to be here
    !CAL(L,NY,NX)=CAL(L-1,NY,NX)
    !write(*,*) "4"
    !CFE(L,NY,NX)=CFE(L-1,NY,NX)
    !write(*,*) "5"
    !CCA(L,NY,NX)=CCA(L-1,NY,NX)
    !write(*,*) "6"
    !CMG(L,NY,NX)=CMG(L-1,NY,NX)
    !write(*,*) "7"
    !CNA(L,NY,NX)=CNA(L-1,NY,NX)
    !write(*,*) "8"
    !CKA(L,NY,NX)=CKA(L-1,NY,NX)
    !write(*,*) "9"
    !CSO4(L,NY,NX)=CSO4(L-1,NY,NX)
    !write(*,*) "10"
    !CCL(L,NY,NX)=CCL(L-1,NY,NX)
    !write(*,*) "11"
    !CALOH(L,NY,NX)=CALOH(L-1,NY,NX)
    !write(*,*) "12"
    !CFEOH(L,NY,NX)=CFEOH(L-1,NY,NX)
    !write(*,*) "13"
    !CCACO(L,NY,NX)=CCACO(L-1,NY,NX)
    !write(*,*) "14"
    !CCASO(L,NY,NX)=CCASO(L-1,NY,NX)
    !write(*,*) "15"
    !CALPO(L,NY,NX)=CALPO(L-1,NY,NX)
    !write(*,*) "16"
    !CFEPO(L,NY,NX)=CFEPO(L-1,NY,NX)
    !write(*,*) "17"
    !CCAPD(L,NY,NX)=CCAPD(L-1,NY,NX)
    !write(*,*) "18"
    !CCAPH(L,NY,NX)=CCAPH(L-1,NY,NX)
    !write(*,*) "19"
    !GKC4(L,NY,NX)=GKC4(L-1,NY,NX)
    !write(*,*) "20"
    !GKCH(L,NY,NX)=GKCH(L-1,NY,NX)
    !write(*,*) "21"
    !GKCA(L,NY,NX)=GKCA(L-1,NY,NX)
    !write(*,*) "22"
    !GKCM(L,NY,NX)=GKCM(L-1,NY,NX)
    !write(*,*) "23"
    !GKCN(L,NY,NX)=GKCN(L-1,NY,NX)
    !write(*,*) "24"
    !GKCK(L,NY,NX)=GKCK(L-1,NY,NX)
    !write(*,*) "25"
    !THW(L,NY,NX)=THW(L-1,NY,NX)
    !write(*,*) "26"
    !THI(L,NY,NX)=THI(L-1,NY,NX)
    !write(*,*) "27"
    !ISOIL(1:4,L,NY,NX)=ISOIL(1:4,L-1,NY,NX)
    !write(*,*) "28"
    !RSC(k_fine_litr,L,NY,NX)=0.0_r8
    !write(*,*) "29"
    !RSN(k_fine_litr,L,NY,NX)=0.0_r8
    !write(*,*) "30"
    !RSP(k_fine_litr,L,NY,NX)=0.0_r8
    !write(*,*) "31"
    !RSC(k_woody_litr,L,NY,NX)=0.0_r8
    !write(*,*) "32"
    !RSN(k_woody_litr,L,NY,NX)=0.0_r8
    !write(*,*) "33"
    !RSP(k_woody_litr,L,NY,NX)=0.0_r8
    !write(*,*) "34"
    !RSC(k_manure,L,NY,NX)=0.0_r8
    !write(*,*) "35"
    !RSN(k_manure,L,NY,NX)=0.0_r8
    !write(*,*) "36"
    !RSP(k_manure,L,NY,NX)=0.0_r8
  ENDDO
  END associate
  end subroutine SetDeepSoil

!------------------------------------------------------------------------------------------
  subroutine CalcSoilThermConductivity(N1,N2,N3,DTKX,TCND1)
  implicit none
  integer , intent(in) :: N1,N2,N3
  real(r8), intent(in) :: DTKX     !absolute temeprature gradient
  real(r8), intent(out):: TCND1

  real(r8) :: HeatDiffusByWat1,HeatDiffusByAir1,RYLXW1,RYLXA1,RYLNW1,RYLNA1
  REAL(R8) :: XNUSW1,XNUSA1,ThermalConducByWater,ThermalConducByAir,WTHET1

  IF(SoiBulkDensity(N3,N2,N1).GT.ZERO.OR.FracSoiPAsWat(N3,N2,N1)+FracSoiPAsIce(N3,N2,N1).GT.ZERO)THEN
    !it is a soil layer or pure water layer
    HeatDiffusByWat1=AZMAX1(FracSoiPAsWat(N3,N2,N1)-TRBW)**3._r8
    HeatDiffusByAir1=AZMAX1(FracSoiPAsAir(N3,N2,N1)-TRBA)**3._r8
    RYLXW1=DTKX*HeatDiffusByWat1
    RYLXA1=DTKX*HeatDiffusByAir1
    RYLNW1=AMIN1(1.0E+04_r8,RYLXW*RYLXW1)
    RYLNA1=AMIN1(1.0E+04_r8,RYLXA*RYLXA1)
    XNUSW1=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW1**0.25_r8/DNUSW)
    XNUSA1=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA1**0.25_r8/DNUSA)
    ThermalConducByWater=2.067E-03_r8*XNUSW1
    ThermalConducByAir=9.050E-05_r8*XNUSA1
    WTHET1=1.467_r8-0.467_r8*FracSoilAsAirt(N3,N2,N1)
    TCND1=(STC(N3,N2,N1)+FracSoiPAsWat(N3,N2,N1)*ThermalConducByWater &
      +0.611_r8*FracSoiPAsIce(N3,N2,N1)*7.844E-03_r8 &
      +WTHET1*FracSoiPAsAir(N3,N2,N1)*ThermalConducByAir) &
      /(DTC(N3,N2,N1)+FracSoiPAsWat(N3,N2,N1)+0.611_r8*FracSoiPAsIce(N3,N2,N1) &
      +WTHET1*FracSoiPAsAir(N3,N2,N1))
  ELSE
    TCND1=0.0_r8
  ENDIF
  end subroutine CalcSoilThermConductivity
end  module SoilPhysParaMod
