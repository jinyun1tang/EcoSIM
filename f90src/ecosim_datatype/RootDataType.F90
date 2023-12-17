module RootDataType

!
!!
! data types of plant characteristics
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  integer,target,allocatable ::  NumRootAxes_pft(:,:,:)                          !root primary axis number, [-]
  integer,target,allocatable ::  NIXBotRootLayer_rpft(:,:,:,:)                       !maximum soil layer number for root axes, [-]
  integer,target,allocatable ::  iPlantRootState_pft(:,:,:)                        !flag to detect root system death , [-]
  integer,target,allocatable ::  NIXBotRootLayer_pft(:,:,:)              !maximum soil layer number for all root axes, [-]
  integer,target,allocatable ::  NI(:,:,:)                           !maximum soil layer number for all root axes, [-]
  real(r8),target,allocatable ::  RootBiomGrowthYield(:,:,:)                        !root growth yield, [g g-1]
  real(r8),target,allocatable ::  MinNonstructuralC4InitRoot_pft(:,:,:)                          !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(r8),target,allocatable ::  RootFracRemobilizableBiom(:,:,:)                       !fraction of remobilizable nonstructural biomass in root, [-]
  real(r8),target,allocatable ::  RootVolPerMassC_pft(:,:,:,:)                      !root volume:mass ratio, [m3 g-1]
  real(r8),target,allocatable ::  Max1stRootRadius1(:,:,:,:)                    !root diameter primary axes, [m]
  real(r8),target,allocatable ::  Max2ndRootRadius1(:,:,:,:)                    !root diameter secondary axes, [m]
  real(r8),target,allocatable ::  PrimRootXSecArea(:,:,:,:)          !root cross-sectional area primary axes, [m2]
  real(r8),target,allocatable ::  SecndRootXSecArea(:,:,:,:)                    !root  cross-sectional area  secondary axes, [m2]
  real(r8),target,allocatable ::  fTgrowRootP(:,:,:,:)                      !root layer temperature growth functiom, [-]
  real(r8),target,allocatable ::  RootrNC_pft(:,:,:)                        !root N:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootrPC_pft(:,:,:)                        !root P:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootPorosity(:,:,:,:)                      !root porosity, [m3 m-3]
  real(r8),target,allocatable ::  RSRR(:,:,:,:)                      !root radial resistivity, [MPa h m-2]
  real(r8),target,allocatable ::  RSRA(:,:,:,:)                      !root axial resistivity, [MPa h m-4]
  real(r8),target,allocatable ::  ShutRutNonstructElmntConducts_pft(:,:,:)                       !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8),target,allocatable ::  VmaxNH4Root_pft(:,:,:,:)                    !maximum root NH4 uptake rate, [g m-2 h-1]
  real(r8),target,allocatable ::  KmNH4Root_pft(:,:,:,:)                    !Km for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  CMinNH4Root_pft(:,:,:,:)                    !minimum NH4 concentration for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  VmaxNO3Root_pft(:,:,:,:)                    !maximum root NO3 uptake rate, [g m-2 h-1]
  real(r8),target,allocatable ::  KmNO3Root_pft(:,:,:,:)                    !Km for root NO3 uptake, [g m-3]
  real(r8),target,allocatable ::  CminNO3Root_pft(:,:,:,:)                    !minimum NO3 concentration for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  VmaxPO4Root_pft(:,:,:,:)                    !maximum root PO4 uptake rate, [g m-2 h-1]
  real(r8),target,allocatable ::  KmPO4Root_pft(:,:,:,:)                    !Km for root PO4 uptake, [g m-3]
  real(r8),target,allocatable ::  CMinPO4Root_pft(:,:,:,:)                    !minimum PO4 concentration for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  RRADP(:,:,:,:)                     !root internal radius, [m]
  real(r8),target,allocatable ::  CNRTS(:,:,:)                       !root N:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  CPRTS(:,:,:)                       !root P:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  Max1stRootRadius(:,:,:,:)                    !maximum radius of primary roots, [m]
  real(r8),target,allocatable ::  Max2ndRootRadius(:,:,:,:)                    !maximum radius of secondary roots, [m]
  real(r8),target,allocatable ::  RootBranchFreq_pft(:,:,:)                        !root brancing frequency, [m-1]
  real(r8),target,allocatable ::  RootPoreTortu4Gas(:,:,:,:)                     !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8),target,allocatable ::  RootNoduleNonstructElmnt_vr(:,:,:,:,:)                  !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootLenPerPopu_pvr(:,:,:,:,:)                   !root layer length per plant, [m p-1]
  real(r8),target,allocatable ::  PrimRootLen(:,:,:,:,:,:)                 !root layer length primary axes, [m d-2]
  real(r8),target,allocatable ::  SecndRootLen(:,:,:,:,:,:)                 !root layer length secondary axes, [m d-2]
  real(r8),target,allocatable ::  RootLenthDensPerPopu_pvr(:,:,:,:,:)          !root length density in soil layers, [m m-3]
  real(r8),target,allocatable ::  PrimRootXNumL_pvr(:,:,:,:,:)                    !root layer number primary axes, [d-2]
  real(r8),target,allocatable ::  SecndRootXNum_pvr(:,:,:,:,:)                    !root layer number axes, [d-2]
  real(r8),target,allocatable ::  SecndRootXNum_rpvr(:,:,:,:,:,:)                  !root layer number secondary axes, [d-2]
  real(r8),target,allocatable ::  AveSecndRootLen(:,:,:,:,:)                   !root layer average length, [m]
  real(r8),target,allocatable ::  RootAreaPerPlant_vr(:,:,:,:,:)                   !root layer area per plant, [m p-1]
  real(r8),target,allocatable ::  RootVH2O_vr(:,:,:,:,:)                   !root layer volume water, [m2 d-2]
  real(r8),target,allocatable ::  PrimRootRadius_pvr(:,:,:,:,:)                   !root layer diameter primary axes, [m ]
  real(r8),target,allocatable ::  RootVolume_vr(:,:,:,:,:)                   !root layer volume air, [m2 d-2]
  real(r8),target,allocatable ::  PrimRootDepth(:,:,:,:,:)                   !root layer depth, [m]
  real(r8),target,allocatable ::  SecndRootRadius_pvr(:,:,:,:,:)                   !root layer diameter secondary axes, [m ]
  real(r8),target,allocatable ::  PrimRootSpecLen(:,:,:,:)                    !specific root length primary axes, [m g-1]
  real(r8),target,allocatable ::  SecndRootSpecLen(:,:,:,:)                    !specific root length secondary axes, [m g-1]
  real(r8),target,allocatable ::  AllPlantRootH2OUptake_vr(:,:,:,:,:)                   !root water uptake, [m2 d-2 h-1]
  real(r8),target,allocatable ::  PSIRoot_vr(:,:,:,:,:)                   !root total water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootOSMO_vr(:,:,:,:,:)                   !root osmotic water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootTurg_vr(:,:,:,:,:)                   !root turgor water potential , [Mpa]
  real(r8),target,allocatable ::  trcg_rootml_vr(:,:,:,:,:,:)           !root gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  trcs_rootml_vr(:,:,:,:,:,:)           !root dissolved gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  TRootGasLossDisturb_pft(:,:,:)                 !total root gas content, [g d-2]
  real(r8),target,allocatable ::  RootBiomCPerPlant_pft(:,:,:)                       !root C per plant, [g p-1]
  real(r8),target,allocatable ::  RootChemElmnts_pft(:,:,:,:)                     !plant root element, [g d-2]
  real(r8),target,allocatable ::  RootStructChemElmnt_pft(:,:,:,:)                    !plant root structural element, [g d-2]
  real(r8),target,allocatable ::  RootProteinC_pvr(:,:,:,:,:)                   !root layer protein C, [g d-2]
  real(r8),target,allocatable ::  Root1stStructChemElmnt_pvr(:,:,:,:,:,:,:)              !root layer element primary axes, [g d-2]
  real(r8),target,allocatable ::  Root2ndStructChemElmnt_pvr(:,:,:,:,:,:,:)              !root layer element secondary axes, [g d-2]
  real(r8),target,allocatable ::   PopuPlantRootC_vr(:,:,:,:,:)                   !root layer C, [g d-2]
  real(r8),target,allocatable ::  RootNodueChemElmnt_pvr(:,:,:,:,:)                  !root layer nodule element, [g d-2]
  real(r8),target,allocatable ::  NoduleChemElmnts_pft(:,:,:,:)                     !root total nodule mass, [g d-2]
  real(r8),target,allocatable ::  RootStructBiomC_vr(:,:,:,:,:)                   !root layer structural C, [g d-2]
  real(r8),target,allocatable ::   RootMycoNonstructElmnt_vr(:,:,:,:,:,:)                !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootNonstructElementConcpft_vr(:,:,:,:,:,:)                !root  layer nonstructural element concentration, [g g-1]
  real(r8),target,allocatable ::  Root1stChemElmnt(:,:,:,:,:,:)                   !root C primary axes, [g d-2]
  real(r8),target,allocatable ::  RootProteinConc_pftvr(:,:,:,:,:)                  !root layer protein C concentration, [g g-1]
!----------------------------------------------------------------------

contains
  subroutine InitRootData(jroots)

  implicit none
  integer, intent(in) :: jroots

  allocate(NumRootAxes_pft(JP,JY,JX));      NumRootAxes_pft=0
  allocate(NIXBotRootLayer_rpft(MaxNumRootAxes,JP,JY,JX));  NIXBotRootLayer_rpft=1  !set to one to avoid numerical failure
  allocate(iPlantRootState_pft(JP,JY,JX));    iPlantRootState_pft=iDead
  allocate(NIXBotRootLayer_pft(JP,JY,JX));      NIXBotRootLayer_pft=0
  allocate(NI(JP,JY,JX));       NI=0
  allocate(RootBiomGrowthYield(JP,JY,JX));     RootBiomGrowthYield=0._r8
  allocate(MinNonstructuralC4InitRoot_pft(JP,JY,JX));       MinNonstructuralC4InitRoot_pft=0._r8
  allocate(RootFracRemobilizableBiom(JP,JY,JX));    RootFracRemobilizableBiom=0._r8
  allocate(RootVolPerMassC_pft(jroots,JP,JY,JX));   RootVolPerMassC_pft=0._r8
  allocate(Max1stRootRadius1(jroots,JP,JY,JX)); Max1stRootRadius1=0._r8
  allocate(Max2ndRootRadius1(jroots,JP,JY,JX)); Max2ndRootRadius1=0._r8
  allocate(PrimRootXSecArea(jroots,JP,JY,JX)); PrimRootXSecArea=0._r8
  allocate(SecndRootXSecArea(jroots,JP,JY,JX)); SecndRootXSecArea=0._r8
  allocate(fTgrowRootP(JZ,JP,JY,JX));  fTgrowRootP=0._r8
  allocate(RootrNC_pft(JP,JY,JX));     RootrNC_pft=0._r8
  allocate(RootrPC_pft(JP,JY,JX));     RootrPC_pft=0._r8
  allocate(RootPorosity(jroots,JP,JY,JX));   RootPorosity=0._r8
  allocate(RSRR(jroots,JP,JY,JX));   RSRR=0._r8
  allocate(RSRA(jroots,JP,JY,JX));   RSRA=0._r8
  allocate(ShutRutNonstructElmntConducts_pft(JP,JY,JX));    ShutRutNonstructElmntConducts_pft=0._r8
  allocate(VmaxNH4Root_pft(jroots,JP,JY,JX)); VmaxNH4Root_pft=0._r8
  allocate(KmNH4Root_pft(jroots,JP,JY,JX)); KmNH4Root_pft=0._r8
  allocate(CMinNH4Root_pft(jroots,JP,JY,JX)); CMinNH4Root_pft=0._r8
  allocate(VmaxNO3Root_pft(jroots,JP,JY,JX)); VmaxNO3Root_pft=0._r8
  allocate(KmNO3Root_pft(jroots,JP,JY,JX)); KmNO3Root_pft=0._r8
  allocate(CminNO3Root_pft(jroots,JP,JY,JX)); CminNO3Root_pft=0._r8
  allocate(VmaxPO4Root_pft(jroots,JP,JY,JX)); VmaxPO4Root_pft=0._r8
  allocate(KmPO4Root_pft(jroots,JP,JY,JX)); KmPO4Root_pft=0._r8
  allocate(CMinPO4Root_pft(jroots,JP,JY,JX)); CMinPO4Root_pft=0._r8
  allocate(RRADP(jroots,JP,JY,JX));  RRADP=0._r8
  allocate(CNRTS(JP,JY,JX));    CNRTS=0._r8
  allocate(CPRTS(JP,JY,JX));    CPRTS=0._r8
  allocate(Max1stRootRadius(jroots,JP,JY,JX)); Max1stRootRadius=0._r8
  allocate(Max2ndRootRadius(jroots,JP,JY,JX)); Max2ndRootRadius=0._r8
  allocate(RootBranchFreq_pft(JP,JY,JX));     RootBranchFreq_pft=0._r8
  allocate(RootPoreTortu4Gas(jroots,JP,JY,JX));  RootPoreTortu4Gas=0._r8
  allocate(RootNoduleNonstructElmnt_vr(NumOfPlantChemElmnts,JZ,JP,JY,JX));RootNoduleNonstructElmnt_vr=0._r8
  allocate(RootLenPerPopu_pvr(jroots,JZ,JP,JY,JX));RootLenPerPopu_pvr=0._r8
  allocate(PrimRootLen(jroots,JZ,JC,JP,JY,JX));PrimRootLen=0._r8
  allocate(SecndRootLen(jroots,JZ,JC,JP,JY,JX));SecndRootLen=0._r8
  allocate(RootLenthDensPerPopu_pvr(jroots,JZ,JP,JY,JX));RootLenthDensPerPopu_pvr=0._r8
  allocate(PrimRootXNumL_pvr(jroots,JZ,JP,JY,JX));PrimRootXNumL_pvr=0._r8
  allocate(SecndRootXNum_pvr(jroots,JZ,JP,JY,JX));SecndRootXNum_pvr=0._r8
  allocate(SecndRootXNum_rpvr(jroots,JZ,JC,JP,JY,JX));SecndRootXNum_rpvr=0._r8
  allocate(AveSecndRootLen(jroots,JZ,JP,JY,JX));AveSecndRootLen=0._r8
  allocate(RootAreaPerPlant_vr(jroots,JZ,JP,JY,JX));RootAreaPerPlant_vr=0._r8
  allocate(RootVH2O_vr(jroots,JZ,JP,JY,JX));RootVH2O_vr=0._r8
  allocate(PrimRootRadius_pvr(jroots,JZ,JP,JY,JX));PrimRootRadius_pvr=0._r8
  allocate(RootVolume_vr(jroots,JZ,JP,JY,JX));RootVolume_vr=0._r8
  allocate(PrimRootDepth(jroots,JC,JP,JY,JX));PrimRootDepth=0._r8
  allocate(SecndRootRadius_pvr(jroots,JZ,JP,JY,JX));SecndRootRadius_pvr=0._r8
  allocate(PrimRootSpecLen(jroots,JP,JY,JX)); PrimRootSpecLen=0._r8
  allocate(SecndRootSpecLen(jroots,JP,JY,JX)); SecndRootSpecLen=0._r8
  allocate(AllPlantRootH2OUptake_vr(jroots,JZ,JP,JY,JX));AllPlantRootH2OUptake_vr=0._r8
  allocate(PSIRoot_vr(jroots,JZ,JP,JY,JX));PSIRoot_vr=0._r8
  allocate(PSIRootOSMO_vr(jroots,JZ,JP,JY,JX));PSIRootOSMO_vr=0._r8
  allocate(PSIRootTurg_vr(jroots,JZ,JP,JY,JX));PSIRootTurg_vr=0._r8
  allocate(trcg_rootml_vr(idg_beg:idg_end-1,2,JZ,JP,JY,JX)); trcg_rootml_vr =0._r8
  allocate(trcs_rootml_vr(idg_beg:idg_end-1,2,JZ,JP,JY,JX)); trcs_rootml_vr =0._r8
  allocate(TRootGasLossDisturb_pft(idg_beg:idg_end-1,JY,JX));TRootGasLossDisturb_pft=0._r8
  allocate(RootBiomCPerPlant_pft(JP,JY,JX));    RootBiomCPerPlant_pft=0._r8
  allocate(RootChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX)); RootChemElmnts_pft=0._r8
  allocate(RootStructChemElmnt_pft(NumOfPlantChemElmnts,JP,JY,JX));   RootStructChemElmnt_pft=0._r8
  allocate(RootProteinC_pvr(jroots,JZ,JP,JY,JX));RootProteinC_pvr=0._r8
  allocate(Root1stStructChemElmnt_pvr(NumOfPlantChemElmnts,jroots,JZ,JC,JP,JY,JX));Root1stStructChemElmnt_pvr=0._r8
  allocate(Root2ndStructChemElmnt_pvr(NumOfPlantChemElmnts,jroots,JZ,JC,JP,JY,JX));Root2ndStructChemElmnt_pvr=0._r8
  allocate( PopuPlantRootC_vr(jroots,JZ,JP,JY,JX)); PopuPlantRootC_vr=0._r8
  allocate(RootNodueChemElmnt_pvr(NumOfPlantChemElmnts,JZ,JP,JY,JX)); RootNodueChemElmnt_pvr=0._r8
  allocate(NoduleChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));  NoduleChemElmnts_pft=0._r8
  allocate(RootStructBiomC_vr(jroots,JZ,JP,JY,JX));RootStructBiomC_vr=0._r8
  allocate(RootMycoNonstructElmnt_vr(NumOfPlantChemElmnts,jroots,JZ,JP,JY,JX)); RootMycoNonstructElmnt_vr=0._r8
  allocate(RootNonstructElementConcpft_vr(NumOfPlantChemElmnts,jroots,JZ,JP,JY,JX));RootNonstructElementConcpft_vr=0._r8
  allocate(Root1stChemElmnt(NumOfPlantChemElmnts,jroots,MaxNumRootAxes,JP,JY,JX));Root1stChemElmnt=0._r8
  allocate(RootProteinConc_pftvr(jroots,JZ,JP,JY,JX));RootProteinConc_pftvr=0._r8
  end subroutine InitRootData

!----------------------------------------------------------------------
  subroutine DestructRootData
  use abortutils, only : destroy
  implicit none
  call destroy(NumRootAxes_pft)
  call destroy(NIXBotRootLayer_rpft)
  call destroy(iPlantRootState_pft)
  call destroy(NIXBotRootLayer_pft)
  call destroy(NI)
  call destroy(RootBiomGrowthYield)
  call destroy(MinNonstructuralC4InitRoot_pft)
  call destroy(RootFracRemobilizableBiom)
  call destroy(RootVolPerMassC_pft)
  call destroy(Max1stRootRadius1)
  call destroy(Max2ndRootRadius1)
  call destroy(PrimRootXSecArea)
  call destroy(SecndRootXSecArea)
  call destroy(fTgrowRootP)
  call destroy(RootrNC_pft)
  call destroy(RootrPC_pft)
  call destroy(RootPorosity)
  call destroy(RSRR)
  call destroy(RSRA)
  call destroy(ShutRutNonstructElmntConducts_pft)
  call destroy(VmaxNH4Root_pft)
  call destroy(KmNH4Root_pft)
  call destroy(CMinNH4Root_pft)
  call destroy(VmaxNO3Root_pft)
  call destroy(KmNO3Root_pft)
  call destroy(CminNO3Root_pft)
  call destroy(VmaxPO4Root_pft)
  call destroy(KmPO4Root_pft)
  call destroy(CMinPO4Root_pft)
  call destroy(RRADP)
  call destroy(CNRTS)
  call destroy(CPRTS)
  call destroy(Max1stRootRadius)
  call destroy(Max2ndRootRadius)
  call destroy(RootBranchFreq_pft)
  call destroy(RootPoreTortu4Gas)
  call destroy(RootNoduleNonstructElmnt_vr)
  call destroy(RootLenPerPopu_pvr)
  call destroy(PrimRootLen)
  call destroy(SecndRootLen)
  call destroy(RootLenthDensPerPopu_pvr)
  call destroy(PrimRootXNumL_pvr)
  call destroy(SecndRootXNum_pvr)
  call destroy(SecndRootXNum_rpvr)
  call destroy(AveSecndRootLen)
  call destroy(RootAreaPerPlant_vr)
  call destroy(RootVH2O_vr)
  call destroy(PrimRootRadius_pvr)
  call destroy(RootVolume_vr)
  call destroy(PrimRootDepth)
  call destroy(SecndRootRadius_pvr)
  call destroy(PrimRootSpecLen)
  call destroy(SecndRootSpecLen)
  call destroy(AllPlantRootH2OUptake_vr)
  call destroy(PSIRoot_vr)
  call destroy(PSIRootOSMO_vr)
  call destroy(PSIRootTurg_vr)
  call destroy(trcg_rootml_vr)
  call destroy(trcs_rootml_vr)
  call destroy(TRootGasLossDisturb_pft)
  call destroy(RootBiomCPerPlant_pft)
  call destroy(RootChemElmnts_pft)
  call destroy(RootStructChemElmnt_pft)
  call destroy(RootProteinC_pvr)
  call destroy(Root1stStructChemElmnt_pvr)
  call destroy(Root2ndStructChemElmnt_pvr)
  call destroy( PopuPlantRootC_vr)
  call destroy(RootNodueChemElmnt_pvr)
  call destroy(NoduleChemElmnts_pft)
  call destroy(RootStructBiomC_vr)
  call destroy(RootMycoNonstructElmnt_vr)
  call destroy(RootNonstructElementConcpft_vr)
  call destroy(Root1stChemElmnt)
  call destroy(RootProteinConc_pftvr)
  end subroutine DestructRootData

end module RootDataType
