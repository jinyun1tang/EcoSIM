module SoilPhysParaMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  USE SoilWaterDataType
  use SoilPropertyDataType
  USE SoilPhysDataType
  use EcoSIMCtrlDataType
  use HydroThermData
  use MiniMathMod, only : safe_adb
  use PhysPars
  use EcoSiMParDataMod, only : micpar  
  use EcoSimConst
  use FlagDataType
  USE AqueChemDatatype
  USE SOMDataType
  use SoilBGCDataType
  USE GridDataType
implicit none
  private
  character(len=*), parameter :: mod_filename=__FILE__
  PUBLIC :: CalcSoilWaterPotential
  public :: SetDeepSoil
  contains
!------------------------------------------------------------------------------------------
  subroutine CalcSoilWaterPotential(NY,NX,N1,N2,N3,PSISA1,THETA1S)
  implicit none
  integer, intent(in) :: NY,NX,N1,N2,N3
  real(r8), intent(out) :: PSISA1
  real(r8), optional, intent(out) :: THETA1S
  real(r8) :: FCDX,FCX,FCLX,WPLX,PSDX,WPX,THETA1
  !     BKVL=soil mass
  !     FC,WP=water contents at field capacity,wilting point
  !     FCL,WPL=log FC,WP
  !     FCD,PSD=FCL-WPL,log(POROS)-FCL
  !     PSISA1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
  !     PSIMX,PSIMD,PSIMS=log water potential at FC,WP,saturation
  !     PSISD=PSIMX-PSIMS
  !     SRP=parameter for deviation from linear log-log water retention
  !     PSISO=osmotic potential
  !     DTHETW=minimum water content for numerical purpose
  ! soil matric potential upper layer
  
  THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1),safe_adb(VWatMicP1(N3,N2,N1),VOLY(N3,N2,N1))))
  IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
    !source layer is active soil  
    IF(THETA1.LT.FC(N3,N2,N1))THEN
      !water less than field capacity
      !PSIHY is the minimum water potential allowed
      PSISA1=AMAX1(PSIHY,-EXP(PSIMX(N2,N1)+((FCL(N3,N2,N1)-LOG(THETA1))/FCD(N3,N2,N1)*PSIMD(N2,N1))))
    ELSEIF(THETA1.LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISA1=-EXP(PSIMS(N2,N1)+(((PSL(N3,N2,N1)-LOG(THETA1))/PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
    ELSE
      THETA1=POROS(N3,N2,N1)    
      PSISA1=PSISE(N3,N2,N1)
    ENDIF
    !
    !     SUBSURFCE UPPER WATER LAYER
    !
    !     THETIX,THETWX=ice,water concentration
    !     FCI,WPI=ice field capacity,wilting point
    !     PSISA1=matric water potential
    !
  ELSEIF(VSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1).and.THETIX(N3,N2,N1)>ZEROS2(N2,N1))THEN
    FCX=FCI*THETIX(N3,N2,N1)
    WPX=WPI*THETIX(N3,N2,N1)
    FCLX=LOG(FCX)
    WPLX=LOG(WPX)
    PSDX=PSL(N3,N2,N1)-FCLX
    FCDX=FCLX-WPLX
    IF(THETWX(N3,N2,N1).LT.FCX)THEN
      PSISA1=AMAX1(PSIHY,-EXP(PSIMX(N2,N1)+((FCLX-LOG(THETWX(N3,N2,N1)))/FCDX*PSIMD(NY,NX))))
    ELSEIF(THETWX(N3,N2,N1).LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISA1=-EXP(PSIMS(N2,N1)+(((PSL(N3,N2,N1)-LOG(THETWX(N3,N2,N1)))/PSDX)*PSISD(N2,N1)))
    ELSE
      !saturated
      THETA1=POROS(N3,N2,N1)
      PSISA1=PSISE(N3,N2,N1)
    ENDIF
  ELSE
    THETA1=POROS(N3,N2,N1)
    PSISA1=PSISE(N3,N2,N1)
  ENDIF

  if(present(THETA1S))THETA1S=THETA1
  end subroutine CalcSoilWaterPotential
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

  DO L=NM+1,JZ
    CumDepth2LayerBottom(L,NY,NX)=2.0_r8*CumDepth2LayerBottom(L-1,NY,NX)-1.0_r8*CumDepth2LayerBottom(L-2,NY,NX)
    BKDSI(L,NY,NX)=BKDSI(L-1,NY,NX)
    FC(L,NY,NX)=FC(L-1,NY,NX)
    WP(L,NY,NX)=WP(L-1,NY,NX)
    SCNV(L,NY,NX)=SCNV(L-1,NY,NX)
    SCNH(L,NY,NX)=SCNH(L-1,NY,NX)
    CSAND(L,NY,NX)=CSAND(L-1,NY,NX)
    CSILT(L,NY,NX)=CSILT(L-1,NY,NX)
    CCLAY(L,NY,NX)=CCLAY(L-1,NY,NX)
    FHOL(L,NY,NX)=FHOL(L-1,NY,NX)
    ROCK(L,NY,NX)=ROCK(L-1,NY,NX)
    PH(L,NY,NX)=PH(L-1,NY,NX)
    CEC(L,NY,NX)=CEC(L-1,NY,NX)
    AEC(L,NY,NX)=AEC(L-1,NY,NX)

!     IF(IDTBL(NY,NX).EQ.0)THEN
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
    CNH4(L,NY,NX)=CNH4(L-1,NY,NX)
    CNO3(L,NY,NX)=CNO3(L-1,NY,NX)
    CPO4(L,NY,NX)=CPO4(L-1,NY,NX)
    CAL(L,NY,NX)=CAL(L-1,NY,NX)
    CFE(L,NY,NX)=CFE(L-1,NY,NX)
    CCA(L,NY,NX)=CCA(L-1,NY,NX)
    CMG(L,NY,NX)=CMG(L-1,NY,NX)
    CNA(L,NY,NX)=CNA(L-1,NY,NX)
    CKA(L,NY,NX)=CKA(L-1,NY,NX)
    CSO4(L,NY,NX)=CSO4(L-1,NY,NX)
    CCL(L,NY,NX)=CCL(L-1,NY,NX)
    CALOH(L,NY,NX)=CALOH(L-1,NY,NX)
    CFEOH(L,NY,NX)=CFEOH(L-1,NY,NX)
    CCACO(L,NY,NX)=CCACO(L-1,NY,NX)
    CCASO(L,NY,NX)=CCASO(L-1,NY,NX)
    CALPO(L,NY,NX)=CALPO(L-1,NY,NX)
    CFEPO(L,NY,NX)=CFEPO(L-1,NY,NX)
    CCAPD(L,NY,NX)=CCAPD(L-1,NY,NX)
    CCAPH(L,NY,NX)=CCAPH(L-1,NY,NX)
    GKC4(L,NY,NX)=GKC4(L-1,NY,NX)
    GKCH(L,NY,NX)=GKCH(L-1,NY,NX)
    GKCA(L,NY,NX)=GKCA(L-1,NY,NX)
    GKCM(L,NY,NX)=GKCM(L-1,NY,NX)
    GKCN(L,NY,NX)=GKCN(L-1,NY,NX)
    GKCK(L,NY,NX)=GKCK(L-1,NY,NX)
    THW(L,NY,NX)=THW(L-1,NY,NX)
    THI(L,NY,NX)=THI(L-1,NY,NX)
    ISOIL(1:4,L,NY,NX)=ISOIL(1:4,L-1,NY,NX)
    RSC(k_fine_litr,L,NY,NX)=0.0_r8
    RSN(k_fine_litr,L,NY,NX)=0.0_r8
    RSP(k_fine_litr,L,NY,NX)=0.0_r8
    RSC(k_woody_litr,L,NY,NX)=0.0_r8
    RSN(k_woody_litr,L,NY,NX)=0.0_r8
    RSP(k_woody_litr,L,NY,NX)=0.0_r8
    RSC(k_manure,L,NY,NX)=0.0_r8
    RSN(k_manure,L,NY,NX)=0.0_r8
    RSP(k_manure,L,NY,NX)=0.0_r8
  ENDDO
  END associate
  end subroutine SetDeepSoil
end  module SoilPhysParaMod
