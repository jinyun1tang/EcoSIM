module PlantAPI4Uptake
!
! interface to integrate the plant model
! for prescribed phenology
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar, pltpar
  use SoilPhysDataType, only: SurfAlbedo_col
  use MiniMathMod,      only: AZMAX1
  use NumericalAuxMod
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
  public :: PlantUptakeAPISend
  public :: PlantUPtakeAPIRecv

  contains
!------------------------------------------------------------------------------------------

  subroutine PlantUptakeAPISend(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(in) :: M,L,NZ,N

  !sent variables also modified
  plt_site%NumActivePlants     = NumActivePlants_col(NY,NX)
  plt_site%PlantPopu_col       = PlantPopu_col(NY,NX)
  plt_site%ZERO                = ZERO                      !numerical threshold
  plt_site%ZERO2               = ZERO2                     !numerical threshold
  plt_morph%LeafStalkArea_col  = LeafStalkArea_col(NY,NX)
  plt_morph%CanopyLeafArea_col = CanopyLeafArea_col(NY,NX)
  plt_site%NL                  = NL_col(NY,NX)
  plt_site%NP0                 = NP0_col(NY,NX)
  plt_site%MaxNumRootLays      = NL_col(NY,NX)
  plt_site%NP                  = NP_col(NY,NX)
  plt_site%NU                  = NUM_col(NY,NX)
  plt_site%NK                  = NK_col(NY,NX)

  DO L=1,NK_col(NY,NX)
    plt_site%CumSoilThickness_vr(L)  = CumSoilThickness_vr(L,NY,NX)
    plt_site%AREA3(L)                = AREA_3D(3,L,NY,NX)
    plt_ew%TKS_vr(L)                 = TKS_vr(L,NY,NX)
    plt_soilchem%THETW_vr(L)         = THETW_vr(L,NY,NX)
    plt_soilchem%SoilWatAirDry_vr(L) = SoilWatAirDry_vr(L,NY,NX)
    plt_soilchem%VLSoilMicP_vr(L)    = VLSoilMicP_vr(L,NY,NX)
    plt_soilchem%VLiceMicP_vr(L)     = VLiceMicP_vr(L,NY,NX)
    plt_soilchem%VLWatMicP_vr(L)     = VLWatMicP_vr(L,NY,NX)
    plt_soilchem%VLMicP_vr(L)        = VLMicP_vr(L,NY,NX)
  ENDDO


  DO NZ=1,NP0_col(NY,NX)
    plt_ew%CanopyBiomWater_pft(NZ)    = CanopyBiomWater_pft(NZ,NY,NX)  
    plt_morph%MaxSoiL4Root_pft(NZ)    = NK_col(NY,NX)                     !can be derived from root type, set to maximum for simplicity
    plt_morph%Myco_pft(NZ)            = 1                                 !no mycorrhizae
    plt_biom%CanopyLeafShethC_pft(NZ) = CanopyLeafShethC_pft(NZ,NY,NX)    !need to convert from leaf area
    plt_biom%CanopyStalkC_pft(NZ)     = CanopyStalkC_pft(NZ,NY,NX)        !need to convert from stem area
    plt_ew%WatHeldOnCanopy_pft(NZ)    = WatHeldOnCanopy_pft(NZ,NY,NX)     !water held by canopy surface
    plt_rad%RadPARbyCanopy_pft(NZ)    = RadPARbyCanopy_pft(NZ,NY,NX)      !computed from surface energy module
    plt_rad%RadSWbyCanopy_pft(NZ)     = RadSWbyCanopy_pft(NZ,NY,NX)       !computed from surface energy module
    plt_ew%PrecIntcptByCanopy_pft(NZ) = PrecIntcptByCanopy_pft(NZ,NY,NX)  !computed from hour1, rainfall partition
    plt_pheno%IsPlantActive_pft(NZ)   = IsPlantActive_pft(NZ,NY,NX)       !lai >0, active
    plt_site%PlantPopulation_pft(NZ)  = PlantPopulation_pft(NZ,NY,NX)
    plt_biom%ZERO4LeafVar_pft(NZ)     = ZERO4LeafVar_pft(NZ,NY,NX)
    plt_biom%ZERO4Groth_pft(NZ)       = ZERO4Groth_pft(NZ,NY,NX)
    plt_morph%CanopyHeight_pft(NZ)    = CanopyHeight_pft(NZ,NY,NX)
    plt_ew%PSICanopy_pft(NZ)          = PSICanopy_pft(NZ,NY,NX)
    DO L=1,NumCanopyLayers
      DO  M=1,NumOfSkyAzimuthSects
        DO  N=1,NumLeafZenithSectors
          plt_rad%RadPAR_zsec(N,M,L,NZ)    = RadPAR_zsec(N,M,L,NZ,NY,NX)
          plt_rad%RadDifPAR_zsec(N,M,L,NZ) = RadDifPAR_zsec(N,M,L,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NK_col(NY,NX)
    DO M=1,NPH
      plt_site%VLWatMicPM_vr(M,L)     = VLWatMicPM_vr(M,L,NY,NX)  !micropore soil moisture
    ENDDO
  ENDDO    

  end subroutine PlantUptakeAPISend
!------------------------------------------------------------------------------------------

  subroutine PlantUPtakeAPIRecv(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NY,NX

  
  end subroutine PlantUPtakeAPIRecv

end module PlantAPI4Uptake