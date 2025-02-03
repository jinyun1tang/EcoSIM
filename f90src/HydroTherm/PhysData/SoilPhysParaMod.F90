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
  !
  !Description
  !Compute soil matric potential
  implicit none
  integer, intent(in) :: NY,NX,N1,N2,N3
  real(r8), intent(out) :: PSISoilMatric      ! soil matric potential
  real(r8), optional, intent(out) :: THETA1S  ! porosity
  real(r8) :: FCDX,FCX,FCLX,WPLX,PSDX,WPX
  real(r8) :: THETA1                          !volumetric soil moisture
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

  THETA1=AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicP1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))

  IF(VLSoilMicPMass_vr(N3,N2,N1).GT.ZEROS(NY,NX))THEN
    !source layer is active soil
    IF(THETA1.LT.FieldCapacity_vr(N3,N2,N1))THEN
      !water less than field capacity
      !PSIHY is the minimum water potential allowed, where hygroscopic water is < 0
      !held tightly on the surfaces of soil particles and exists as a thin layer of vapor
      PSISoilMatric=AMAX1(PSIHY,-EXP(LOGPSIFLD(N2,N1)+((LOGFldCapacity_vr(N3,N2,N1)-LOG(THETA1))/FCD_vr(N3,N2,N1)*LOGPSIMND(N2,N1))))
    ELSEIF(THETA1.LT.POROS_vr(N3,N2,N1)-DTHETW)THEN
      PSISoilMatric=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS_vr(N3,N2,N1)-LOG(THETA1))/PSD_vr(N3,N2,N1))**SRP_vr(N3,N2,N1)*LOGPSIMXD(N2,N1)))
    ELSE
      THETA1        = POROS_vr(N3,N2,N1)
      PSISoilMatric = PSISE_vr(N3,N2,N1)
    ENDIF
    !
    !     SUBSURFCE UPPER WATER LAYER
    !
    !     FracSoiPAsIce,FracSoiPAsWat=ice,water concentration
    !     FCI,WPI=ice field capacity,wilting point
    !     PSISoilMatric=matric water potential
    !
  ELSEIF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1).and.FracSoiPAsIce_vr(N3,N2,N1)>ZEROS2(N2,N1))THEN
    FCX  = FCI*FracSoiPAsIce_vr(N3,N2,N1)
    WPX  = WPI*FracSoiPAsIce_vr(N3,N2,N1)
    FCLX = LOG(FCX)
    WPLX = LOG(WPX)
    PSDX = LOGPOROS_vr(N3,N2,N1)-FCLX
    FCDX = FCLX-WPLX
    IF(FracSoiPAsWat_vr(N3,N2,N1).LT.FCX)THEN
      PSISoilMatric=AMAX1(PSIHY,-EXP(LOGPSIFLD(N2,N1)+((FCLX-LOG(FracSoiPAsWat_vr(N3,N2,N1)))/FCDX*LOGPSIMND(NY,NX))))
    ELSEIF(FracSoiPAsWat_vr(N3,N2,N1).LT.POROS_vr(N3,N2,N1)-DTHETW)THEN
      PSISoilMatric=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS_vr(N3,N2,N1)-LOG(FracSoiPAsWat_vr(N3,N2,N1)))/PSDX)*LOGPSIMXD(N2,N1)))
    ELSE
      !saturated
      THETA1        = POROS_vr(N3,N2,N1)
      PSISoilMatric = PSISE_vr(N3,N2,N1)
    ENDIF
  ELSE
    THETA1        = POROS_vr(N3,N2,N1)
    PSISoilMatric = PSISE_vr(N3,N2,N1)
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

  DO L=NM+1,JZ
    CumDepz2LayBottom_vr(L,NY,NX) = 2.0_r8*CumDepz2LayBottom_vr(L-1,NY,NX)-1.0_r8*CumDepz2LayBottom_vr(L-2,NY,NX)
    SoiBulkDensityt0_vr(L,NY,NX) = SoiBulkDensityt0_vr(L-1,NY,NX)
    FieldCapacity_vr(L,NY,NX)    = FieldCapacity_vr(L-1,NY,NX)
    WiltPoint_vr(L,NY,NX)        = WiltPoint_vr(L-1,NY,NX)
    SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L-1,NY,NX)
    SatHydroCondHrzn_vr(L,NY,NX) = SatHydroCondHrzn_vr(L-1,NY,NX)
    CSAND_vr(L,NY,NX)               = CSAND_vr(L-1,NY,NX)
    CSILT(L,NY,NX)               = CSILT(L-1,NY,NX)
    CCLAY_vr(L,NY,NX)               = CCLAY_vr(L-1,NY,NX)
    SoilFracAsMacP_vr(L,NY,NX)   = SoilFracAsMacP_vr(L-1,NY,NX)
    ROCK_vr(L,NY,NX)             = ROCK_vr(L-1,NY,NX)
    PH_vr(L,NY,NX)                  = PH_vr(L-1,NY,NX)
    CEC_vr(L,NY,NX)              = CEC_vr(L-1,NY,NX)
    AEC_vr(L,NY,NX)              = AEC_vr(L-1,NY,NX)

!     IF(IDWaterTable_col(NY,NX).EQ.0)THEN
!       0.25_r8 is the geometric decreasing ratio (tunable)
!  or different scheme can be used
    CSoilOrgM_vr(ielmc,L,NY,NX) = 0.25_r8*CSoilOrgM_vr(ielmc,L-1,NY,NX)
    COMLitrC_vr(L,NY,NX)        = 0.25_r8*COMLitrC_vr(L-1,NY,NX)
    CSoilOrgM_vr(ielmn,L,NY,NX) = 0.25_r8*CSoilOrgM_vr(ielmn,L-1,NY,NX)
    CSoilOrgM_vr(ielmp,L,NY,NX) = 0.25_r8*CSoilOrgM_vr(ielmp,L-1,NY,NX)
!     ELSE
!       CSoilOrgM_vr(ielmc,L,NY,NX)=CSoilOrgM_vr(ielmc,L-1,NY,NX)
!       COMLitrC_vr(L,NY,NX)=COMLitrC_vr(L-1,NY,NX)
!       CSoilOrgM_vr(ielmn,L,NY,NX)=CSoilOrgM_vr(ielmn,L-1,NY,NX)
!       CSoilOrgM_vr(ielmp,L,NY,NX)=CSoilOrgM_vr(ielmp,L-1,NY,NX)
!     ENDIF
    CNH4_vr(L,NY,NX)             = CNH4_vr(L-1,NY,NX)
    CNO3_vr(L,NY,NX)             = CNO3_vr(L-1,NY,NX)
    CPO4_vr(L,NY,NX)             = CPO4_vr(L-1,NY,NX)
    CAL_vr(L,NY,NX)              = CAL_vr(L-1,NY,NX)
    CFE_vr(L,NY,NX)              = CFE_vr(L-1,NY,NX)
    CCA_vr(L,NY,NX)              = CCA_vr(L-1,NY,NX)
    CMG_vr(L,NY,NX)              = CMG_vr(L-1,NY,NX)
    CNA_vr(L,NY,NX)              = CNA_vr(L-1,NY,NX)
    CKA_vr(L,NY,NX)              = CKA_vr(L-1,NY,NX)
    CSO4_vr(L,NY,NX)             = CSO4_vr(L-1,NY,NX)
    CCL_vr(L,NY,NX)              = CCL_vr(L-1,NY,NX)
    CALOH_vr(L,NY,NX)            = CALOH_vr(L-1,NY,NX)
    CFEOH_vr(L,NY,NX)            = CFEOH_vr(L-1,NY,NX)
    CCACO_vr(L,NY,NX)            = CCACO_vr(L-1,NY,NX)
    CCASO_vr(L,NY,NX)            = CCASO_vr(L-1,NY,NX)
    CALPO_vr(L,NY,NX)            = CALPO_vr(L-1,NY,NX)
    CFEPO_vr(L,NY,NX)            = CFEPO_vr(L-1,NY,NX)
    CCAPD_vr(L,NY,NX)            = CCAPD_vr(L-1,NY,NX)
    CCAPH_vr(L,NY,NX)            = CCAPH_vr(L-1,NY,NX)
    GKC4_vr(L,NY,NX)             = GKC4_vr(L-1,NY,NX)
    GKCH_vr(L,NY,NX)             = GKCH_vr(L-1,NY,NX)
    GKCA_vr(L,NY,NX)             = GKCA_vr(L-1,NY,NX)
    GKCM_vr(L,NY,NX)             = GKCM_vr(L-1,NY,NX)
    GKCN_vr(L,NY,NX)             = GKCN_vr(L-1,NY,NX)
    GKCK_vr(L,NY,NX)             = GKCK_vr(L-1,NY,NX)
    THW_vr(L,NY,NX)           = THW_vr(L-1,NY,NX)
    THI_vr(L,NY,NX)           = THI_vr(L-1,NY,NX)
    ISOIL(1:4,L,NY,NX)        = ISOIL(1:4,L-1,NY,NX)
    RSC(k_fine_litr,L,NY,NX)  = 0.0_r8
    RSN(k_fine_litr,L,NY,NX)  = 0.0_r8
    RSP(k_fine_litr,L,NY,NX)  = 0.0_r8
    RSC(k_woody_litr,L,NY,NX) = 0.0_r8
    RSN(k_woody_litr,L,NY,NX) = 0.0_r8
    RSP(k_woody_litr,L,NY,NX) = 0.0_r8
    RSC(k_manure,L,NY,NX)     = 0.0_r8
    RSN(k_manure,L,NY,NX)     = 0.0_r8
    RSP(k_manure,L,NY,NX)     = 0.0_r8
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

  IF(SoilBulkDensity_vr(N3,N2,N1).GT.ZERO.OR.FracSoiPAsWat_vr(N3,N2,N1)+FracSoiPAsIce_vr(N3,N2,N1).GT.ZERO)THEN
    !it is a soil layer or pure water layer
    HeatDiffusByWat1     = AZMAX1(FracSoiPAsWat_vr(N3,N2,N1)-TRBW)**3._r8
    HeatDiffusByAir1     = AZMAX1(AirFilledSoilPore_vr(N3,N2,N1)-TRBA)**3._r8
    RYLXW1               = DTKX*HeatDiffusByWat1
    RYLXA1               = DTKX*HeatDiffusByAir1
    RYLNW1               = AMIN1(1.0E+04_r8,RYLXW*RYLXW1)
    RYLNA1               = AMIN1(1.0E+04_r8,RYLXA*RYLXA1)
    XNUSW1               = AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW1**0.25_r8/DNUSW)
    XNUSA1               = AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA1**0.25_r8/DNUSA)
    ThermalConducByWater = 2.067E-03_r8*XNUSW1
    ThermalConducByAir   = 9.050E-05_r8*XNUSA1
    WTHET1               = 1.467_r8-0.467_r8*FracSoilAsAirt(N3,N2,N1)
    TCND1                = (NumerSolidThermCond_vr(N3,N2,N1)+FracSoiPAsWat_vr(N3,N2,N1)*ThermalConducByWater &
      +0.611_r8*FracSoiPAsIce_vr(N3,N2,N1)*7.844E-03_r8+WTHET1*AirFilledSoilPore_vr(N3,N2,N1)*ThermalConducByAir) &
      /(DenomSolidThermCond_vr(N3,N2,N1)+FracSoiPAsWat_vr(N3,N2,N1)+0.611_r8*FracSoiPAsIce_vr(N3,N2,N1) &
      +WTHET1*AirFilledSoilPore_vr(N3,N2,N1))
  ELSE
    TCND1=0.0_r8
  ENDIF
  end subroutine CalcSoilThermConductivity
end  module SoilPhysParaMod
