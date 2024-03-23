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
  integer,target,allocatable ::  NIXBotRootLayer_rpft(:,:,:,:)                   !maximum soil layer number for root axes, [-]
  integer,target,allocatable ::  iPlantRootState_pft(:,:,:)                      !flag to detect root system death , [-]
  integer,target,allocatable ::  NIXBotRootLayer_pft(:,:,:)              !maximum soil layer number for all root axes, [-]
  integer,target,allocatable ::  MaxSoiL4Root(:,:,:)                           !maximum soil layer number for all root axes, [-]
  real(r8),target,allocatable ::  RootBiomGrowthYield(:,:,:)                        !root growth yield, [g g-1]
  real(r8),target,allocatable ::  MinNonstructuralC4InitRoot_pft(:,:,:)                          !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(r8),target,allocatable ::  RootFracRemobilizableBiom(:,:,:)                       !fraction of remobilizable nonstructural biomass in root, [-]
  real(r8),target,allocatable ::  RootVolPerMassC_pft(:,:,:,:)                      !root volume:mass ratio, [m3 g-1]
  real(r8),target,allocatable ::  Max1stRootRadius_pft1(:,:,:,:)                    !root diameter primary axes, [m]
  real(r8),target,allocatable ::  Max2ndRootRadius_pft1(:,:,:,:)                    !root diameter secondary axes, [m]
  real(r8),target,allocatable ::  Root1stXSecArea_pft(:,:,:,:)          !root cross-sectional area primary axes, [m2]
  real(r8),target,allocatable ::  Root2ndXSecArea_pft(:,:,:,:)                    !root  cross-sectional area  secondary axes, [m2]
  real(r8),target,allocatable ::  fTgrowRootP_vr(:,:,:,:)                      !root layer temperature growth functiom, [-]
  real(r8),target,allocatable ::  RootrNC_pft(:,:,:)                        !root N:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootrPC_pft(:,:,:)                        !root P:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootPorosity_pft(:,:,:,:)                      !root porosity, [m3 m-3]
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
  real(r8),target,allocatable ::  RootRaidus_rpft(:,:,:,:)                     !root internal radius, [m]
  real(r8),target,allocatable ::  CNRTS(:,:,:)                       !root N:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  CPRTS(:,:,:)                       !root P:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  Max1stRootRadius_pft(:,:,:,:)                    !maximum radius of primary roots, [m]
  real(r8),target,allocatable ::  Max2ndRootRadius_pft(:,:,:,:)                    !maximum radius of secondary roots, [m]
  real(r8),target,allocatable ::  RootBranchFreq_pft(:,:,:)                        !root brancing frequency, [m-1]
  real(r8),target,allocatable ::  RootPoreTortu4Gas(:,:,:,:)                     !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8),target,allocatable ::  RootNodulNonstElm_pvr(:,:,:,:,:)                  !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootLenPerPlant_pvr(:,:,:,:,:)                   !root layer length per plant, [m p-1]
  real(r8),target,allocatable ::  Root1stLen_rpvr(:,:,:,:,:,:)                 !root layer length primary axes, [m d-2]
  real(r8),target,allocatable ::  Root2ndLen_pvr(:,:,:,:,:,:)                 !root layer length secondary axes, [m d-2]
  real(r8),target,allocatable ::  RootLenDensPerPlant_pvr(:,:,:,:,:)          !root length density in soil layers, [m m-3]
  real(r8),target,allocatable ::  Root1stXNumL_pvr(:,:,:,:,:)                    !root layer number primary axes, [d-2]
  real(r8),target,allocatable ::  Root2ndXNum_pvr(:,:,:,:,:)                    !root layer number axes, [d-2]
  real(r8),target,allocatable ::  Root2ndXNum_rpvr(:,:,:,:,:,:)                  !root layer number secondary axes, [d-2]
  real(r8),target,allocatable ::  AveLen2ndRoot_pvr(:,:,:,:,:)                   !root layer average length, [m]
  real(r8),target,allocatable ::  RootAreaPerPlant_pvr(:,:,:,:,:)                   !root layer area per plant, [m p-1]
  real(r8),target,allocatable ::  RootVH2O_pvr(:,:,:,:,:)                   !root layer volume water, [m2 d-2]
  real(r8),target,allocatable ::  Root1stRadius_pvr(:,:,:,:,:)                   !root layer diameter primary axes, [m ]
  real(r8),target,allocatable ::  RootPoreVol_pvr(:,:,:,:,:)                   !root layer volume air, [m2 d-2]
  real(r8),target,allocatable ::  Root1stDepz_pft(:,:,:,:,:)                   !root layer depth, [m]
  real(r8),target,allocatable ::  Radius2ndRoot_pvr(:,:,:,:,:)                   !root layer diameter secondary axes, [m ]
  real(r8),target,allocatable ::  Root1stSpecLen_pft(:,:,:,:)                    !specific root length primary axes, [m g-1]
  real(r8),target,allocatable ::  Root2ndSpecLen_pft(:,:,:,:)                    !specific root length secondary axes, [m g-1]
  real(r8),target,allocatable ::  AllPlantRootH2OUptake_vr(:,:,:,:,:)                   !root water uptake, [m2 d-2 h-1]
  real(r8),target,allocatable ::  PSIRoot_pvr(:,:,:,:,:)                   !root total water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootOSMO_vr(:,:,:,:,:)                   !root osmotic water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootTurg_vr(:,:,:,:,:)                   !root turgor water potential , [Mpa]
  real(r8),target,allocatable ::  trcg_rootml_pvr(:,:,:,:,:,:)           !root gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  trcs_rootml_pvr(:,:,:,:,:,:)           !root dissolved gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  TRootGasLossDisturb_pft(:,:,:)                 !total root gas content, [g d-2]
  real(r8),target,allocatable ::  RootBiomCPerPlant_pft(:,:,:)                       !root C per plant, [g p-1]
  real(r8),target,allocatable ::  RootElms_pft(:,:,:,:)                     !plant root element, [g d-2]
  real(r8),target,allocatable ::  RootStructElmnt_pft(:,:,:,:)                    !plant root structural element, [g d-2]
  real(r8),target,allocatable ::  RootProteinC_pvr(:,:,:,:,:)                   !root layer protein C, [g d-2]
  real(r8),target,allocatable ::  Root1stStructElm_rpvr(:,:,:,:,:,:,:)              !root layer element primary axes, [g d-2]
  real(r8),target,allocatable ::  Root2ndStructElm_pvr(:,:,:,:,:,:,:)              !root layer element secondary axes, [g d-2]
  real(r8),target,allocatable ::   PopuPlantRootC_vr(:,:,:,:,:)                   !root layer C, [g d-2]
  real(r8),target,allocatable ::  RootNodulElm_pvr(:,:,:,:,:)                  !root layer nodule element, [g d-2]
  real(r8),target,allocatable ::  NoduleChemElms_pft(:,:,:,:)                     !root total nodule mass, [g d-2]
  real(r8),target,allocatable ::  RootMycoActiveBiomC_pvr(:,:,:,:,:)                   !root layer structural C, [g d-2]
  real(r8),target,allocatable ::   RootMycoNonstElm_pvr(:,:,:,:,:,:)                !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootNonstructElmConc_pvr(:,:,:,:,:,:)                !root  layer nonstructural element concentration, [g g-1]
  real(r8),target,allocatable ::  Root1stElm_raxs(:,:,:,:,:,:)                   !root C primary axes, [g d-2]
  real(r8),target,allocatable ::  RootProteinConc_pvr(:,:,:,:,:)                  !root layer protein C concentration, [g g-1]
!----------------------------------------------------------------------

contains
  subroutine InitRootData(jroots)

  implicit none
  integer, intent(in) :: jroots

  allocate(NumRootAxes_pft(JP,JY,JX));      NumRootAxes_pft=0
  allocate(NIXBotRootLayer_rpft(MaxNumRootAxes,JP,JY,JX));  NIXBotRootLayer_rpft=1  !set to one to avoid numerical failure
  allocate(iPlantRootState_pft(JP,JY,JX));    iPlantRootState_pft=iDead
  allocate(NIXBotRootLayer_pft(JP,JY,JX));      NIXBotRootLayer_pft=0
  allocate(MaxSoiL4Root(JP,JY,JX));       MaxSoiL4Root=0
  allocate(RootBiomGrowthYield(JP,JY,JX));     RootBiomGrowthYield=0._r8
  allocate(MinNonstructuralC4InitRoot_pft(JP,JY,JX));       MinNonstructuralC4InitRoot_pft=0._r8
  allocate(RootFracRemobilizableBiom(JP,JY,JX));    RootFracRemobilizableBiom=0._r8
  allocate(RootVolPerMassC_pft(jroots,JP,JY,JX));   RootVolPerMassC_pft=0._r8
  allocate(Max1stRootRadius_pft1(jroots,JP,JY,JX)); Max1stRootRadius_pft1=0._r8
  allocate(Max2ndRootRadius_pft1(jroots,JP,JY,JX)); Max2ndRootRadius_pft1=0._r8
  allocate(Root1stXSecArea_pft(jroots,JP,JY,JX)); Root1stXSecArea_pft=0._r8
  allocate(Root2ndXSecArea_pft(jroots,JP,JY,JX)); Root2ndXSecArea_pft=0._r8
  allocate(fTgrowRootP_vr(JZ,JP,JY,JX));  fTgrowRootP_vr=0._r8
  allocate(RootrNC_pft(JP,JY,JX));     RootrNC_pft=0._r8
  allocate(RootrPC_pft(JP,JY,JX));     RootrPC_pft=0._r8
  allocate(RootPorosity_pft(jroots,JP,JY,JX));   RootPorosity_pft=0._r8
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
  allocate(RootRaidus_rpft(jroots,JP,JY,JX));  RootRaidus_rpft=0._r8
  allocate(CNRTS(JP,JY,JX));    CNRTS=0._r8
  allocate(CPRTS(JP,JY,JX));    CPRTS=0._r8
  allocate(Max1stRootRadius_pft(jroots,JP,JY,JX)); Max1stRootRadius_pft=0._r8
  allocate(Max2ndRootRadius_pft(jroots,JP,JY,JX)); Max2ndRootRadius_pft=0._r8
  allocate(RootBranchFreq_pft(JP,JY,JX));     RootBranchFreq_pft=0._r8
  allocate(RootPoreTortu4Gas(jroots,JP,JY,JX));  RootPoreTortu4Gas=0._r8
  allocate(RootNodulNonstElm_pvr(NumPlantChemElms,JZ,JP,JY,JX));RootNodulNonstElm_pvr=0._r8
  allocate(RootLenPerPlant_pvr(jroots,JZ,JP,JY,JX));RootLenPerPlant_pvr=0._r8
  allocate(Root1stLen_rpvr(jroots,JZ,NumOfCanopyLayers,JP,JY,JX));Root1stLen_rpvr=0._r8
  allocate(Root2ndLen_pvr(jroots,JZ,NumOfCanopyLayers,JP,JY,JX));Root2ndLen_pvr=0._r8
  allocate(RootLenDensPerPlant_pvr(jroots,JZ,JP,JY,JX));RootLenDensPerPlant_pvr=0._r8
  allocate(Root1stXNumL_pvr(jroots,JZ,JP,JY,JX));Root1stXNumL_pvr=0._r8
  allocate(Root2ndXNum_pvr(jroots,JZ,JP,JY,JX));Root2ndXNum_pvr=0._r8
  allocate(Root2ndXNum_rpvr(jroots,JZ,NumOfCanopyLayers,JP,JY,JX));Root2ndXNum_rpvr=0._r8
  allocate(AveLen2ndRoot_pvr(jroots,JZ,JP,JY,JX));AveLen2ndRoot_pvr=0._r8
  allocate(RootAreaPerPlant_pvr(jroots,JZ,JP,JY,JX));RootAreaPerPlant_pvr=0._r8
  allocate(RootVH2O_pvr(jroots,JZ,JP,JY,JX));RootVH2O_pvr=0._r8
  allocate(Root1stRadius_pvr(jroots,JZ,JP,JY,JX));Root1stRadius_pvr=0._r8
  allocate(RootPoreVol_pvr(jroots,JZ,JP,JY,JX));RootPoreVol_pvr=0._r8
  allocate(Root1stDepz_pft(jroots,NumOfCanopyLayers,JP,JY,JX));Root1stDepz_pft=0._r8
  allocate(Radius2ndRoot_pvr(jroots,JZ,JP,JY,JX));Radius2ndRoot_pvr=0._r8
  allocate(Root1stSpecLen_pft(jroots,JP,JY,JX)); Root1stSpecLen_pft=0._r8
  allocate(Root2ndSpecLen_pft(jroots,JP,JY,JX)); Root2ndSpecLen_pft=0._r8
  allocate(AllPlantRootH2OUptake_vr(jroots,JZ,JP,JY,JX));AllPlantRootH2OUptake_vr=0._r8
  allocate(PSIRoot_pvr(jroots,JZ,JP,JY,JX));PSIRoot_pvr=0._r8
  allocate(PSIRootOSMO_vr(jroots,JZ,JP,JY,JX));PSIRootOSMO_vr=0._r8
  allocate(PSIRootTurg_vr(jroots,JZ,JP,JY,JX));PSIRootTurg_vr=0._r8
  allocate(trcg_rootml_pvr(idg_beg:idg_end-1,2,JZ,JP,JY,JX)); trcg_rootml_pvr =0._r8
  allocate(trcs_rootml_pvr(idg_beg:idg_end-1,2,JZ,JP,JY,JX)); trcs_rootml_pvr =0._r8
  allocate(TRootGasLossDisturb_pft(idg_beg:idg_end-1,JY,JX));TRootGasLossDisturb_pft=0._r8
  allocate(RootBiomCPerPlant_pft(JP,JY,JX));    RootBiomCPerPlant_pft=0._r8
  allocate(RootElms_pft(NumPlantChemElms,JP,JY,JX)); RootElms_pft=0._r8
  allocate(RootStructElmnt_pft(NumPlantChemElms,JP,JY,JX));   RootStructElmnt_pft=0._r8
  allocate(RootProteinC_pvr(jroots,JZ,JP,JY,JX));RootProteinC_pvr=0._r8
  allocate(Root1stStructElm_rpvr(NumPlantChemElms,jroots,JZ,NumOfCanopyLayers,JP,JY,JX));Root1stStructElm_rpvr=0._r8
  allocate(Root2ndStructElm_pvr(NumPlantChemElms,jroots,JZ,NumOfCanopyLayers,JP,JY,JX));Root2ndStructElm_pvr=0._r8
  allocate( PopuPlantRootC_vr(jroots,JZ,JP,JY,JX)); PopuPlantRootC_vr=0._r8
  allocate(RootNodulElm_pvr(NumPlantChemElms,JZ,JP,JY,JX)); RootNodulElm_pvr=0._r8
  allocate(NoduleChemElms_pft(NumPlantChemElms,JP,JY,JX));  NoduleChemElms_pft=0._r8
  allocate(RootMycoActiveBiomC_pvr(jroots,JZ,JP,JY,JX));RootMycoActiveBiomC_pvr=0._r8
  allocate(RootMycoNonstElm_pvr(NumPlantChemElms,jroots,JZ,JP,JY,JX)); RootMycoNonstElm_pvr=0._r8
  allocate(RootNonstructElmConc_pvr(NumPlantChemElms,jroots,JZ,JP,JY,JX));RootNonstructElmConc_pvr=0._r8
  allocate(Root1stElm_raxs(NumPlantChemElms,jroots,MaxNumRootAxes,JP,JY,JX));Root1stElm_raxs=0._r8
  allocate(RootProteinConc_pvr(jroots,JZ,JP,JY,JX));RootProteinConc_pvr=0._r8
  end subroutine InitRootData

!----------------------------------------------------------------------
  subroutine DestructRootData
  use abortutils, only : destroy
  implicit none
  call destroy(NumRootAxes_pft)
  call destroy(NIXBotRootLayer_rpft)
  call destroy(iPlantRootState_pft)
  call destroy(NIXBotRootLayer_pft)
  call destroy(MaxSoiL4Root)
  call destroy(RootBiomGrowthYield)
  call destroy(MinNonstructuralC4InitRoot_pft)
  call destroy(RootFracRemobilizableBiom)
  call destroy(RootVolPerMassC_pft)
  call destroy(Max1stRootRadius_pft1)
  call destroy(Max2ndRootRadius_pft1)
  call destroy(Root1stXSecArea_pft)
  call destroy(Root2ndXSecArea_pft)
  call destroy(fTgrowRootP_vr)
  call destroy(RootrNC_pft)
  call destroy(RootrPC_pft)
  call destroy(RootPorosity_pft)
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
  call destroy(RootRaidus_rpft)
  call destroy(CNRTS)
  call destroy(CPRTS)
  call destroy(Max1stRootRadius_pft)
  call destroy(Max2ndRootRadius_pft)
  call destroy(RootBranchFreq_pft)
  call destroy(RootPoreTortu4Gas)
  call destroy(RootNodulNonstElm_pvr)
  call destroy(RootLenPerPlant_pvr)
  call destroy(Root1stLen_rpvr)
  call destroy(Root2ndLen_pvr)
  call destroy(RootLenDensPerPlant_pvr)
  call destroy(Root1stXNumL_pvr)
  call destroy(Root2ndXNum_pvr)
  call destroy(Root2ndXNum_rpvr)
  call destroy(AveLen2ndRoot_pvr)
  call destroy(RootAreaPerPlant_pvr)
  call destroy(RootVH2O_pvr)
  call destroy(Root1stRadius_pvr)
  call destroy(RootPoreVol_pvr)
  call destroy(Root1stDepz_pft)
  call destroy(Radius2ndRoot_pvr)
  call destroy(Root1stSpecLen_pft)
  call destroy(Root2ndSpecLen_pft)
  call destroy(AllPlantRootH2OUptake_vr)
  call destroy(PSIRoot_pvr)
  call destroy(PSIRootOSMO_vr)
  call destroy(PSIRootTurg_vr)
  call destroy(trcg_rootml_pvr)
  call destroy(trcs_rootml_pvr)
  call destroy(TRootGasLossDisturb_pft)
  call destroy(RootBiomCPerPlant_pft)
  call destroy(RootElms_pft)
  call destroy(RootStructElmnt_pft)
  call destroy(RootProteinC_pvr)
  call destroy(Root1stStructElm_rpvr)
  call destroy(Root2ndStructElm_pvr)
  call destroy( PopuPlantRootC_vr)
  call destroy(RootNodulElm_pvr)
  call destroy(NoduleChemElms_pft)
  call destroy(RootMycoActiveBiomC_pvr)
  call destroy(RootMycoNonstElm_pvr)
  call destroy(RootNonstructElmConc_pvr)
  call destroy(Root1stElm_raxs)
  call destroy(RootProteinConc_pvr)
  end subroutine DestructRootData

end module RootDataType
