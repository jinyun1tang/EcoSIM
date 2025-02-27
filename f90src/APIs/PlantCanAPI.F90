module PlantCanAPI

! interface to integrate the plant model
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar, pltpar
  use SoilPhysDataType, only: SurfAlbedo_col
  use MiniMathMod,      only: AZMAX1
  use EcoSIMSolverPar
  use EcoSIMHistMod
  use SnowDataType
  use TracerIDMod
  use SurfLitterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use EcoSimSumDataType
  use SoilHeatDataType
  use SOMDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use GridDataType
  use RootDataType
  use SoilWaterDataType
  use CanopyDataType
  use PlantDataRateType
  use PlantTraitDataType
  use CanopyRadDataType
  use FlagDataType
  use EcosimBGCFluxType
  use FertilizerDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use PlantAPIData
implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantAPICanMSend
  public :: PlantAPICanMRecv

  contains
!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMSend(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,N,M,NN,NZ,K,NB
!  Integers
  plt_site%KoppenClimZone        = KoppenClimZone_col(NY,NX)
  plt_site%NP                    = NP(NY,NX)
  plt_site%NU                    = NU(NY,NX)
  plt_morph%StemArea_col         = StemArea_col(NY,NX)
  plt_morph%CanopyLeafArea_col   = CanopyLeafArea_col(NY,NX)
  plt_site%ZEROS                 = ZEROS(NY,NX)
  plt_site%ZERO                  = ZERO
  plt_ew%SnowDepth               = SnowDepth_col(NY,NX)
  plt_ew%TairK                   = TairK_col(NY,NX)
  plt_morph%CanopyHeight_col     = CanopyHeight_col(NY,NX)
  plt_site%WindMesureHeight_col  = WindMesureHeight_col(NY,NX)
  plt_ew%ZERO4PlantDisplace_col  = ZERO4PlantDisplace_col(NY,NX)
  plt_site%WindSpeedAtm_col      = WindSpeedAtm_col(NY,NX)
  plt_ew%VLHeatCapSnowMin_col    = VLHeatCapSnowMin_col(NY,NX)
  plt_ew%VLHeatCapSurfSnow_col   = VLHeatCapSnow_snvr(1,NY,NX)
  plt_morph%CanopyHeightZ_col(0) = CanopyHeightZ_col(0,NY,NX)
  DO L=1,NumOfCanopyLayers
    plt_morph%CanopyStemAareZ_col(L) = CanopyStemAareZ_col(L,NY,NX)
    plt_morph%CanopyLeafAareZ_col(L) = CanopyLeafAareZ_col(L,NY,NX)
    plt_morph%CanopyHeightZ_col(L)   = CanopyHeightZ_col(L,NY,NX)
    plt_rad%TAU_DirRadTransm(L)      = TAU_DirRadTransm(L,NY,NX)
  ENDDO
  plt_rad%TAU_DirRadTransm(NumOfCanopyLayers+1)=TAU_DirRadTransm(NumOfCanopyLayers+1,NY,NX)

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)                 = AREA(3,L,NY,NX)
    plt_soilchem%VLSoilPoreMicP_vr(L) = VLSoilPoreMicP_vr(L,NY,NX)
    plt_soilchem%VLSoilMicP_vr(L)     = VLSoilMicP_vr(L,NY,NX)
    plt_soilchem%VLWatMicP_vr(L)      = VLWatMicP_vr(L,NY,NX)
  ENDDO
  plt_rad%SineSunInclAngle_col     = SineSunInclAngle_col(NY,NX)
  plt_site%SolarNoonHour_col       = SolarNoonHour_col(NY,NX)
  plt_morph%LeafStalkArea_col      = LeafStalkArea_col(NY,NX)
  plt_rad%GroundSurfAzimuth_col    = GroundSurfAzimuth_col(NY,NX)
  plt_rad%CosineGrndSlope_col      = CosineGrndSlope_col(NY,NX)
  plt_rad%SineGrndSlope_col        = SineGrndSlope_col(NY,NX)
  plt_rad%RadSWDiffus_col          = RadSWDiffus_col(NY,NX)
  plt_rad%RadPARDiffus_col         = RadPARDiffus_col(NY,NX)
  plt_rad%RadSWDirect_col          = RadSWDirect_col(NY,NX)
  plt_rad%RadPARDirect_col         = RadPARDirect_col(NY,NX)
  plt_site%SoilSurfRoughnesst0_col = SoilSurfRoughnesst0_col(NY,NX)
  plt_ew%VcumWatSnow_col           = VcumWatSnow_col(NY,NX)
  plt_ew%VcumIceSnow_col           = VcumIceSnow_col(NY,NX)
  plt_ew%VcumDrySnoWE_col          = VcumDrySnoWE_col(NY,NX)
  plt_rad%TotSineSkyAngles_grd     = TotSineSkyAngles_grd
  plt_rad%SoilAlbedo               = SoilAlbedo_col(NY,NX)
  plt_rad%SurfAlbedo_col           = SurfAlbedo_col(NY,NX)
  plt_site%ZEROS2                  = ZEROS2(NY,NX)
  plt_site%POROS1                  = POROS_vr(NU(NY,NX),NY,NX)
  DO NZ=1,NP(NY,NX)
    plt_morph%CanopyLeafArea_pft(NZ)   = CanopyLeafArea_pft(NZ,NY,NX)
    plt_morph%CanopyHeight_pft(NZ)     = CanopyHeight_pft(NZ,NY,NX)
    plt_morph%ClumpFactorNow_pft(NZ)   = ClumpFactorNow_pft(NZ,NY,NX)
    plt_rad%LeafSWabsorpty_pft(NZ)     = LeafSWabsorpty_pft(NZ,NY,NX)
    plt_rad%LeafPARabsorpty_pft(NZ)    = LeafPARabsorpty_pft(NZ,NY,NX)
    plt_rad%RadSWLeafTransmis_pft(NZ)  = RadSWLeafTransmis_pft(NZ,NY,NX)
    plt_rad%RadSWLeafAlbedo_pft(NZ)    = RadSWLeafAlbedo_pft(NZ,NY,NX)
    plt_rad%RadPARLeafTransmis_pft(NZ) = RadPARLeafTransmis_pft(NZ,NY,NX)
    plt_rad%CanopyPARalbedo_pft(NZ)    = CanopyPARalbedo_pft(NZ,NY,NX)
    plt_morph%NumOfBranches_pft(NZ)    = NumOfBranches_pft(NZ,NY,NX)
    plt_morph%ClumpFactor_pft(NZ)      = ClumpFactor_pft(NZ,NY,NX)

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO K=0,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers
          plt_morph%CanopyLeafArea_lpft(L,K,NB,NZ)=CanopyLeafArea_lpft(L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,NumOfCanopyLayers
        plt_morph%CanopyStalkArea_lbrch(L,NB,NZ)=CanopyStalkArea_lbrch(L,NB,NZ,NY,NX)
      ENDDO
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers
          DO N=1,NumOfLeafZenithSectors
            plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)=LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      DO  L=1,NumOfCanopyLayers
        DO N=1,NumOfLeafZenithSectors
          plt_morph%StemAreaZsec_brch(N,L,NB,NZ)=StemAreaZsec_brch(N,L,NB,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO N=1,NumOfSkyAzimuthSects
    plt_rad%OMEGAG(N)=OMEGAG(N,NY,NX)
  ENDDO
  DO N=1,NumOfLeafZenithSectors
    plt_rad%CosineLeafAngle(N) = CosineLeafAngle(N)
    plt_rad%SineLeafAngle(N)   = SineLeafAngle(N)
  ENDDO
  DO NN=1,NumOfLeafAzimuthSectors
    DO M=1,NumOfLeafZenithSectors
      DO N=1,NumOfSkyAzimuthSects
        plt_rad%OMEGA(N,M,NN)             = OMEGA(N,M,NN)
        plt_rad%OMEGX(N,M,NN)             = OMEGX(N,M,NN)
        plt_rad%iScatteringDiffus(N,M,NN) = iScatteringDiffus(N,M,NN)
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMSend

!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMRecv(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer :: N,M,NN,L,NZ,K,NB

  ZERO4PlantDisplace_col(NY,NX)  = plt_ew%ZERO4PlantDisplace_col
  RoughHeight_col(NY,NX)         = plt_ew%RoughHeight
  AbvCanopyBndlResist_col(NY,NX) = plt_ew%AbvCanopyBndlResist_col
  RIB_col(NY,NX)                     = plt_ew%RIB
  
  CanopyHeight_col(NY,NX)    = plt_morph%CanopyHeight_col
  RadSWDirect_col(NY,NX)     = plt_rad%RadSWDirect_col
  RadSWDiffus_col(NY,NX)     = plt_rad%RadSWDiffus_col
  RadPARDirect_col(NY,NX)    = plt_rad%RadPARDirect_col
  RadPARDiffus_col(NY,NX)    = plt_rad%RadPARDiffus_col
  RadSWGrnd_col(NY,NX)       = plt_rad%RadSWGrnd_col
  FracSWRad2Grnd_col(NY,NX)  = plt_rad%FracSWRad2Grnd_col
  RadSWSolarBeam_col(NY,NX)  = plt_rad%RadSWSolarBeam_col
  RadPARSolarBeam_col(NY,NX) = plt_rad%RadPARSolarBeam_col
  
  DO L=0,NumOfCanopyLayers
    CanopyHeightZ_col(L,NY,NX)=plt_morph%CanopyHeightZ_col(L)
  ENDDO
  DO L=1,NumOfCanopyLayers
    TAU_DirRadTransm(L,NY,NX) = plt_rad%TAU_DirRadTransm(L)
    TAU_RadThru(L,NY,NX)      = plt_rad%TAU_RadThru(L)
  ENDDO
  LeafStalkArea_col(NY,NX)=plt_morph%LeafStalkArea_col
  DO NZ=1,NP(NY,NX)
    LeafStalkArea_pft(NZ,NY,NX)    =plt_morph%LeafStalkArea_pft(NZ)
    RadSWbyCanopy_pft(NZ,NY,NX)    =plt_rad%RadSWbyCanopy_pft(NZ)
    RadPARbyCanopy_pft(NZ,NY,NX)   =plt_rad%RadPARbyCanopy_pft(NZ)
    ClumpFactorNow_pft(NZ,NY,NX)   =plt_morph%ClumpFactorNow_pft(NZ)
    FracPARads2Canopy_pft(NZ,NY,NX)=plt_rad%FracPARads2Canopy_pft(NZ)
    StomatalStress_pft(NZ,NY,NX)   =plt_biom%StomatalStress_pft(NZ)
    Eco_RadSW_col(NY,NX)           =Eco_RadSW_col(NY,NX)+RadSWbyCanopy_pft(NZ,NY,NX)
    DO L=1,NumOfCanopyLayers
      DO M=1,NumOfSkyAzimuthSects
        DO  N=1,NumOfLeafZenithSectors
          RadDifPAR_zsec(N,M,L,NZ,NY,NX)=plt_rad%RadDifPAR_zsec(N,M,L,NZ)
          RadPAR_zsec(N,M,L,NZ,NY,NX)   =plt_rad%RadPAR_zsec(N,M,L,NZ)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMRecv

end module PlantCanAPI