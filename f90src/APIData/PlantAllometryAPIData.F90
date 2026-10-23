module PlantAllometryAPIData
  ! Owns the allometry API type, instance and lifecycle procedures.
  use PlantAPICommonData
  implicit none
  private
  public :: plt_allom_init, plt_allom_destroy

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTS_pft(:)                     => null()  !root P:C ratio x root growth yield,                      [-]
  real(r8), pointer :: CNRTS_pft(:)                     => null()  !root N:C ratio x root growth yield,                      [-]
  real(r8), pointer :: rNCNodule_pft(:)                 => null()  !nodule N:C ratio,                                        [gN gC-1]
  real(r8), pointer :: rPCNoduler_pft(:)                => null()  !nodule P:C ratio,                                        [gP gC-1]
  real(r8), pointer :: rNCRoot_pft(:)                   => null()  !root N:C ratio,                                          [gN gC-1]
  real(r8), pointer :: rPCRootr_pft(:)                  => null()  !root P:C ratio,                                          [gP gC-1]
  real(r8), pointer :: rProteinC2RootN_pft(:)           => null()  !Protein C to root N ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2RootP_pft(:)           => null()  !Protein C to root P ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2LeafN_pft(:)          => null()  !Protein C to leaf N ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2LeafP_pft(:)          => null()  !Protein C to leaf P ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: NoduGrowthYield_pft(:)           => null()  !nodule growth yield,                                     [g g-1]
  real(r8), pointer :: RootBiomGrosYld_pft(:)           => null()  !root growth yield,                                       [g g-1]
  real(r8), pointer :: rPCEar_pft(:)                    => null()  !ear P:C ratio,                                           [gP gC-1]
  real(r8), pointer :: PetolShethBiomGrowthYld_pft(:)      => null()  !sheath growth yield,                                     [g g-1]
  real(r8), pointer :: rPCHusk_pft(:)                   => null()  !husk P:C ratio,                                          [gP gC-1]
  real(r8), pointer :: StalkBiomGrowthYld_pft(:)        => null()  !stalk growth yield,                                      [gC gC-1]
  real(r8), pointer :: HuskBiomGrowthYld_pft(:)         => null()  !husk growth yield,                                       [gC gC-1]
  real(r8), pointer :: ReserveBiomGrowthYld_pft(:)      => null()  !reserve growth yield,                                    [gC gC-1]
  real(r8), pointer :: GrainBiomGrowthYld_pft(:)        => null()  !grain growth yield,                                      [gC gC-1]
  real(r8), pointer :: EarBiomGrowthYld_pft(:)          => null()  !ear growth yield,                                        [gC gC-1]
  real(r8), pointer :: rNCHusk_pft(:)                   => null()  !husk N:C ratio,                                          [gN gC-1]
  real(r8), pointer :: rNCReserve_pft(:)                => null()  !reserve N:C ratio,                                       [gN gC-1]
  real(r8), pointer :: rNCEar_pft(:)                    => null()  !ear N:C ratio,                                           [gN gC-1]
  real(r8), pointer :: rPCReserve_pft(:)                => null()  !reserve P:C ratio,                                       [gP gC-1]
  real(r8), pointer :: rPCGrain_pft(:)                      => null()  !grain P:C ratio,                                         [gP gP-1]
  real(r8), pointer :: rNCStalk_pft(:)                  => null()  !stalk N:C ratio,                                         [gN gC-1]
  real(r8), pointer :: rECLiveCRoot_pft(:,:)           => null()   ! element:C ratio of live coarse root, [gE gC-1]
  real(r8), pointer :: rECDeadCRoot_pft(:,:)           => null()   ! element:C ratio of dead coarse root, [gE gC-1]
  real(r8), pointer :: rNCLigRoot_pft(:)               => null()  !NC ratio of lignified root, [gN gC-1]
  real(r8), pointer :: rPCLigRoot_pft(:)               => null()  !PC ratio of lignified root, [gP gC-1]
  real(r8), pointer :: FracLeafShethElmAlloc2Litr(:,:)  => null()  !woody element allocation, [-]
  real(r8), pointer :: FracPetolShethAlloc2Litr(:,:) => null()  !leaf element allocation,[-]
  real(r8), pointer :: FracRootElmAllocm(:,:)       => null()  !C woody fraction in root,[-]
  real(r8), pointer :: FracWoodStalkElmAlloc2Litr(:,:)  => null()  !woody element allocation,[-]
  real(r8), pointer :: LeafBiomGrowthYld_pft(:)         => null()  !leaf growth yield,                                       [g g-1]
  real(r8), pointer :: rNCGrain_pft(:)                      => null()  !grain N:C ratio,                                         [g g-1]
  real(r8), pointer :: rPCLeaf_pft(:)                      => null()  !maximum leaf P:C ratio,                                  [g g-1]
  real(r8), pointer :: rPCSheath_pft(:)                     => null()  !sheath P:C ratio,                                        [g g-1]
  real(r8), pointer :: rNCSheath_pft(:)                     => null()  !sheath N:C ratio,                                        [g g-1]
  real(r8), pointer :: rPCStalk_pft(:)                  => null()  !stalk P:C ratio,                                         [g g-1]
  real(r8), pointer :: rNCLeaf_pft(:)                      => null()  !maximum leaf N:C ratio,                                  [g g-1]
  real(r8), pointer :: SingleGrainMeanBiomC_brch(:,:)     => null()  !potential carbon mass per grain, [gC seed-1]
  real(r8), pointer :: FracGroth2Node_pft(:)            => null()  !parameter for allocation of growth to nodes,             [-]
  real(r8), pointer ::RootProteinCMax_pft(:)      => null()  !reference root protein N, [gN g-1]

  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type


  type(plant_allometry_type), public, target :: plt_allom     !plant allometric parameters

contains

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this


  allocate(this%RootProteinCMax_pft(JP1));this%RootProteinCMax_pft=spval
  allocate(this%FracGroth2Node_pft(JP1));this%FracGroth2Node_pft=spval
  allocate(this%SingleGrainMeanBiomC_brch(MaxNumBranches,JP1));this%SingleGrainMeanBiomC_brch=spval
  allocate(this%NoduGrowthYield_pft(JP1));this%NoduGrowthYield_pft=spval
  allocate(this%RootBiomGrosYld_pft(JP1));this%RootBiomGrosYld_pft=spval
  allocate(this%rPCRootr_pft(JP1));this%rPCRootr_pft=spval
  allocate(this%rProteinC2RootN_pft(JP1));this%rProteinC2RootN_pft=spval
  allocate(this%rProteinC2RootP_pft(JP1));this%rProteinC2RootP_pft=spval
  allocate(this%rProteinC2LeafN_pft(JP1));this%rProteinC2LeafN_pft=spval
  allocate(this%rProteinC2LeafP_pft(JP1));this%rProteinC2LeafP_pft=spval
  allocate(this%CPRTS_pft(JP1));this%CPRTS_pft=spval
  allocate(this%CNRTS_pft(JP1));this%CNRTS_pft=spval
  allocate(this%rNCNodule_pft(JP1));this%rNCNodule_pft=spval
  allocate(this%rPCNoduler_pft(JP1));this%rPCNoduler_pft=spval
  allocate(this%rNCRoot_pft(JP1));this%rNCRoot_pft=spval
  allocate(this%rPCLeaf_pft(JP1));this%rPCLeaf_pft=spval
  allocate(this%rPCSheath_pft(JP1));this%rPCSheath_pft=spval
  allocate(this%rNCLeaf_pft(JP1));this%rNCLeaf_pft=spval
  allocate(this%rNCSheath_pft(JP1));this%rNCSheath_pft=spval
  allocate(this%rNCGrain_pft(JP1));this%rNCGrain_pft=spval
  allocate(this%rPCStalk_pft(JP1));this%rPCStalk_pft=spval
  allocate(this%rNCStalk_pft(JP1));this%rNCStalk_pft=spval
  allocate(this%rECDeadCRoot_pft(NumPlantChemElms,JP1));    this%rECDeadCRoot_pft=spval
  allocate(this%rECLiveCRoot_pft(NumPlantChemElms,JP1)); this%rECLiveCRoot_pft=spval
  allocate(this%rNCLigRoot_pft(JP1));this%rNCLigRoot_pft=spval
  allocate(this%rPCLigRoot_pft(JP1)); this%rPCLigRoot_pft=spval
  allocate(this%rPCGrain_pft(JP1));this%rPCGrain_pft=spval
  allocate(this%rPCEar_pft(JP1));this%rPCEar_pft=spval
  allocate(this%rPCReserve_pft(JP1));this%rPCReserve_pft=spval
  allocate(this%rNCReserve_pft(JP1));this%rNCReserve_pft=spval
  allocate(this%rPCHusk_pft(JP1));this%rPCHusk_pft=spval
  allocate(this%FracWoodStalkElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracWoodStalkElmAlloc2Litr=spval
  allocate(this%FracPetolShethAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracPetolShethAlloc2Litr=spval
  allocate(this%FracRootElmAllocm(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracRootElmAllocm=spval
  allocate(this%FracLeafShethElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracLeafShethElmAlloc2Litr=spval

  allocate(this%PetolShethBiomGrowthYld_pft(JP1));this%PetolShethBiomGrowthYld_pft=spval
  allocate(this%HuskBiomGrowthYld_pft(JP1));this%HuskBiomGrowthYld_pft=spval
  allocate(this%StalkBiomGrowthYld_pft(JP1));this%StalkBiomGrowthYld_pft=spval
  allocate(this%ReserveBiomGrowthYld_pft(JP1));this%ReserveBiomGrowthYld_pft=spval
  allocate(this%EarBiomGrowthYld_pft(JP1));this%EarBiomGrowthYld_pft=spval
  allocate(this%GrainBiomGrowthYld_pft(JP1));this%GrainBiomGrowthYld_pft=spval
  allocate(this%rNCHusk_pft(JP1));this%rNCHusk_pft=spval
  allocate(this%rNCEar_pft(JP1));this%rNCEar_pft=spval
  allocate(this%LeafBiomGrowthYld_pft(JP1));this%LeafBiomGrowthYld_pft=spval

  end subroutine plt_allom_init

  subroutine plt_allom_destroy(this)
  implicit none

  class(plant_allometry_type) :: this

!  if(allocated(FracGroth2Node_pft))deallocate(FracGroth2Node_pft)
!  if(allocated(SingleGrainMeanBiomC_brch))deallocate(SingleGrainMeanBiomC_brch)
!  if(allocated(NoduGrowthYield_pft))deallocate(NoduGrowthYield_pft)
!  if(allocated(RootBiomGrosYld_pft))deallocate(RootBiomGrosYld_pft)
!  if(allocated(rPCRootr_pft))deallocate(rPCRootr_pft)
!  if(allocated(rProteinC2LeafN_pft))deallocate(rProteinC2LeafN_pft)
!  if(allocated(rProteinC2LeafP_pft))deallocate(rProteinC2LeafP_pft)
!  if(allocated(CPRTS_pft))deallocate(CPRTS_pft)
!  if(allocated(CNRTS_pft))deallocate(CNRTS_pft)
!  if(allocated(rNCNodule_pft))deallocate(rNCNodule_pft)
!  if(allocated(rPCNoduler_pft))deallocate(rPCNoduler_pft)
!  if(allocated(rNCRoot_pft))deallocate(rNCRoot_pft)
!  if(allocated(RootProteinCMax_pft))deallocate(RootProteinCMax_pft)
!  if(allocated(CPLF))deallocate(CPLF)
!  if(allocated(CPSHE))deallocate(CPSHE)
!  if(allocated(CNSHE))deallocate(CNSHE)
!  if(allocated(CNLF))deallocate(CNLF)
!  if(allocated(LeafBiomGrowthYld_pft))deallocate(LeafBiomGrowthYld_pft)
!  if(allocated(rPCStalk_pft))deallocate(rPCStalk_pft)
!  if(allocated(CNGR))deallocate(CNGR)
!  if(allocated(rNCStalk_pft))deallocate(rNCStalk_pft)
!  if(allocated(CPGR))deallocate(CPGR)
!  if(allocated(PetolShethBiomGrowthYld_pft))deallocate(PetolShethBiomGrowthYld_pft)
!  if(allocated(StalkBiomGrowthYld_pft))deallocate(StalkBiomGrowthYld_pft)
!  if(allocated(GrainBiomGrowthYld_pft))deallocate(GrainBiomGrowthYld_pft)
!  if(allocated(ReserveBiomGrowthYld_pft))deallocate(ReserveBiomGrowthYld_pft)
!  if(allocated(EarBiomGrowthYld_pft))deallocate(EarBiomGrowthYld_pft)
!  if(allocated(HuskBiomGrowthYld_pft))deallocate(HuskBiomGrowthYld_pft)
!  if(allocated(FracPetolShethAlloc2Litr))deallocate(FracPetolShethAlloc2Litr)
!  if(allocated(FracWoodStalkElmAlloc2Litr))deallocate(FracWoodStalkElmAlloc2Litr)
!  if(allocated(FracRootElmAllocm))deallocate(FracRootElmAllocm)
!  if(allocated(FracLeafShethElmAlloc2Litr))deallocate(FracLeafShethElmAlloc2Litr)
!  if(allocated(rPCEar_pft))deallocate(rPCEar_pft)
!  if(allocated(rPCHusk_pft))deallocate(rPCHusk_pft)
!  if(allocated(rNCHusk_pft))deallocate(rNCHusk_pft)
!  if(allocated(rNCEar_pft))deallocate(rNCEar_pft)
!  if(allocated(rNCReserve_pft))deallocate(rNCReserve_pft)
!  if(allocated(rPCReserve_pft))deallocate(rPCReserve_pft)

  end subroutine plt_allom_destroy
end module PlantAllometryAPIData
