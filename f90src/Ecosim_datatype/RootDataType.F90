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
  integer,target,allocatable ::  MaxSoiL4Root_pft(:,:,:)                           !maximum soil layer number for all root axes, [-]
  real(r8),target,allocatable :: RootElmsbeg_pft(:,:,:,:)            !root biomass per pft
  real(r8),target,allocatable ::  RootBiomGrosYld_pft(:,:,:)                        !root growth yield, [g g-1]
  real(r8),target,allocatable ::  MinNonstC2InitRoot_pft(:,:,:)                          !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(r8),target,allocatable ::  RootFracRemobilizableBiom(:,:,:)                       !fraction of remobilizable nonstructural biomass in root, [-]
  real(r8),target,allocatable ::  RootVolPerMassC_pft(:,:,:,:)                      !root volume:mass ratio, [m3 g-1]
  real(r8),target,allocatable ::  Root1stMaxRadius1_pft(:,:,:,:)                    !root diameter primary axes, [m]
  real(r8),target,allocatable ::  Root2ndMaxRadius1_pft(:,:,:,:)                    !root diameter secondary axes, [m]
  real(r8),target,allocatable ::  Root1stXSecArea_pft(:,:,:,:)          !root cross-sectional area primary axes, [m2]
  real(r8),target,allocatable ::  Root2ndXSecArea_pft(:,:,:,:)                    !root  cross-sectional area  secondary axes, [m2]
  real(r8),target,allocatable ::  fTgrowRootP_vr(:,:,:,:)                      !root layer temperature growth functiom, [-]
  real(r8),target,allocatable ::  RootrNC_pft(:,:,:)                        !root N:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootrPC_pft(:,:,:)                        !root P:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootPorosity_pft(:,:,:,:)                      !root porosity, [m3 m-3]
  real(r8),target,allocatable ::  RootRadialResist_pft(:,:,:,:)                      !root radial resistivity, [MPa h m-2]
  real(r8),target,allocatable ::  RootAxialResist_pft(:,:,:,:)                      !root axial resistivity, [MPa h m-4]
  real(r8),target,allocatable ::  ShutRutNonstElmntConducts_pft(:,:,:)                       !shoot-root rate constant for nonstructural C exchange, [h-1]
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
  real(r8),target,allocatable ::  CNRTS_pft(:,:,:)                       !root N:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  CPRTS_pft(:,:,:)                       !root P:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  RootMycoNonstElms_pft(:,:,:,:,:)
  real(r8),target,allocatable ::  Root1stMaxRadius_pft(:,:,:,:)                    !maximum radius of primary roots, [m]
  real(r8),target,allocatable ::  Root2ndMaxRadius_pft(:,:,:,:)                    !maximum radius of secondary roots, [m]
  real(r8),target,allocatable ::  RootBranchFreq_pft(:,:,:)                        !root brancing frequency, [m-1]
  real(r8),target,allocatable ::  RootNodulElms_pft(:,:,:,:)            
  real(r8),target,allocatable ::  RootPoreTortu4Gas(:,:,:,:)                     !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8),target,allocatable ::  RootNodulNonstElms_rpvr(:,:,:,:,:)                  !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootLenPerPlant_pvr(:,:,:,:,:)                   !root layer length per plant, [m p-1]
  real(r8),target,allocatable ::  Root1stLen_rpvr(:,:,:,:,:,:)                 !root layer length primary axes, [m d-2]
  real(r8),target,allocatable ::  Root2ndLen_rpvr(:,:,:,:,:,:)                 !root layer length secondary axes, [m d-2]
  real(r8),target,allocatable ::  RootLenDensPerPlant_pvr(:,:,:,:,:)          !root length density in soil layers, [m m-3]
  real(r8),target,allocatable ::  Root1stXNumL_pvr(:,:,:,:,:)                    !root layer number primary axes, [d-2]
  real(r8),target,allocatable ::  Root2ndXNum_pvr(:,:,:,:,:)                    !root layer number axes, [d-2]
  real(r8),target,allocatable ::  Root2ndXNum_rpvr(:,:,:,:,:,:)                  !root layer number secondary axes, [d-2]
  real(r8),target,allocatable ::  Root2ndAveLen_pvr(:,:,:,:,:)                   !root layer average length, [m]
  real(r8),target,allocatable ::  RootAreaPerPlant_pvr(:,:,:,:,:)                   !root layer area per plant, [m p-1]
  real(r8),target,allocatable ::  RootVH2O_pvr(:,:,:,:,:)                   !root layer volume water, [m2 d-2]
  real(r8),target,allocatable ::  Root1stRadius_pvr(:,:,:,:,:)                   !root layer diameter primary axes, [m ]
  real(r8),target,allocatable ::  RootPoreVol_pvr(:,:,:,:,:)                   !root layer volume air, [m2 d-2]
  real(r8),target,allocatable ::  Root1stDepz_pft(:,:,:,:,:)                   !root layer depth, [m]
  real(r8),target,allocatable ::  Root2ndRadius_pvr(:,:,:,:,:)                   !root layer diameter secondary axes, [m ]
  real(r8),target,allocatable ::  Root1stSpecLen_pft(:,:,:,:)                    !specific root length primary axes, [m g-1]
  real(r8),target,allocatable ::  Root2ndSpecLen_pft(:,:,:,:)                    !specific root length secondary axes, [m g-1]
  real(r8),target,allocatable ::  AllPlantRootH2OLoss_vr(:,:,:,:,:)      !root water uptake, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PSIRoot_pvr(:,:,:,:,:)                   !root total water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootOSMO_vr(:,:,:,:,:)                   !root osmotic water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootTurg_vr(:,:,:,:,:)                   !root turgor water potential , [Mpa]
  real(r8),target,allocatable ::  trcg_rootml_pvr(:,:,:,:,:,:)           !root gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  trcs_rootml_pvr(:,:,:,:,:,:)           !root dissolved gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  TRootGasLossDisturb_col(:,:,:)                 !total root gas content, [g d-2]
  real(r8),target,allocatable ::  RootBiomCPerPlant_pft(:,:,:)                       !root C per plant, [g p-1]
  real(r8),target,allocatable ::  RootElms_pft(:,:,:,:)                     !plant root element, [g d-2]
  real(r8),target,allocatable ::  RootStrutElms_pft(:,:,:,:)                    !plant root structural element, [g d-2]
  real(r8),target,allocatable ::  RootProteinC_pvr(:,:,:,:,:)                   !root layer protein C, [g d-2]
  real(r8),target,allocatable ::  RootMyco1stStrutElms_rpvr(:,:,:,:,:,:,:)              !root layer element primary axes, [g d-2]
  real(r8),target,allocatable ::  RootMyco2ndStrutElms_rpvr(:,:,:,:,:,:,:)              !root layer element secondary axes, [g d-2]
  real(r8),target,allocatable ::   PopuRootMycoC_pvr(:,:,:,:,:)                   !root layer C, [g d-2]
  real(r8),target,allocatable ::  RootNodulStrutElms_rpvr(:,:,:,:,:)                  !root layer nodule element, [g d-2]
  real(r8),target,allocatable ::  NodulStrutElms_pft(:,:,:,:)                     !root total nodule mass, [g d-2]
  real(r8),target,allocatable ::  RootMycoActiveBiomC_pvr(:,:,:,:,:)                   !root layer structural C, [g d-2]
  real(r8),target,allocatable ::   RootMycoNonstElms_rpvr(:,:,:,:,:,:)                !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootNonstructElmConc_rpvr(:,:,:,:,:,:)                !root  layer nonstructural element concentration, [g g-1]
  real(r8),target,allocatable ::  RootMyco1stElm_raxs(:,:,:,:,:,:)                   !root C primary axes, [g d-2]
  real(r8),target,allocatable ::  RootProteinConc_rpvr(:,:,:,:,:)                  !root layer protein C concentration, [g g-1]
  real(r8),target,allocatable :: RootMassElm_vr(:,:,:,:)
  real(r8),target,allocatable :: RootGasConductance_pvr(:,:,:,:,:,:)       !1/hr
!----------------------------------------------------------------------

contains
  subroutine InitRootData(jroots)

  implicit none
  integer, intent(in) :: jroots

  allocate(RootMassElm_vr(NumPlantChemElms,JZ,JY,JX)); RootMassElm_vr =0._r8
  allocate(NumRootAxes_pft(JP,JY,JX));      NumRootAxes_pft=0
  allocate(NIXBotRootLayer_rpft(MaxNumRootAxes,JP,JY,JX));  NIXBotRootLayer_rpft=1  !set to one to avoid numerical failure
  allocate(iPlantRootState_pft(JP,JY,JX));    iPlantRootState_pft=iDead
  allocate(NIXBotRootLayer_pft(JP,JY,JX));      NIXBotRootLayer_pft=0
  allocate(MaxSoiL4Root_pft(JP,JY,JX));       MaxSoiL4Root_pft=0
  allocate(RootElmsbeg_pft(NumPlantChemElms,JP,JY,JX)); RootElmsbeg_pft=0._r8
  allocate(RootBiomGrosYld_pft(JP,JY,JX));     RootBiomGrosYld_pft=0._r8
  allocate(MinNonstC2InitRoot_pft(JP,JY,JX));       MinNonstC2InitRoot_pft=0._r8
  allocate(RootFracRemobilizableBiom(JP,JY,JX));    RootFracRemobilizableBiom=0._r8
  allocate(RootVolPerMassC_pft(jroots,JP,JY,JX));   RootVolPerMassC_pft=0._r8
  allocate(Root1stMaxRadius1_pft(jroots,JP,JY,JX)); Root1stMaxRadius1_pft=0._r8
  allocate(Root2ndMaxRadius1_pft(jroots,JP,JY,JX)); Root2ndMaxRadius1_pft=0._r8
  allocate(Root1stXSecArea_pft(jroots,JP,JY,JX)); Root1stXSecArea_pft=0._r8
  allocate(Root2ndXSecArea_pft(jroots,JP,JY,JX)); Root2ndXSecArea_pft=0._r8
  allocate(fTgrowRootP_vr(JZ,JP,JY,JX));  fTgrowRootP_vr=0._r8
  allocate(RootrNC_pft(JP,JY,JX));     RootrNC_pft=0._r8
  allocate(RootrPC_pft(JP,JY,JX));     RootrPC_pft=0._r8
  allocate(RootGasConductance_pvr(idg_beg:idg_NH3,jroots,JZ,JP,JY,JX));RootGasConductance_pvr=0._r8
  allocate(RootPorosity_pft(jroots,JP,JY,JX));   RootPorosity_pft=0._r8
  allocate(RootRadialResist_pft(jroots,JP,JY,JX));   RootRadialResist_pft=0._r8
  allocate(RootAxialResist_pft(jroots,JP,JY,JX));   RootAxialResist_pft=0._r8
  allocate(ShutRutNonstElmntConducts_pft(JP,JY,JX));    ShutRutNonstElmntConducts_pft=0._r8
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
  allocate(CNRTS_pft(JP,JY,JX));    CNRTS_pft=0._r8
  allocate(CPRTS_pft(JP,JY,JX));    CPRTS_pft=0._r8
  allocate(RootMycoNonstElms_pft(NumPlantChemElms,jroots,JP,JY,JX));RootMycoNonstElms_pft=0._r8
  allocate(RootNodulElms_pft(NumPlantChemElms,JP,JY,JX));RootNodulElms_pft=0._r8
  allocate(Root1stMaxRadius_pft(jroots,JP,JY,JX)); Root1stMaxRadius_pft=0._r8
  allocate(Root2ndMaxRadius_pft(jroots,JP,JY,JX)); Root2ndMaxRadius_pft=0._r8
  allocate(RootBranchFreq_pft(JP,JY,JX));     RootBranchFreq_pft=0._r8
  allocate(RootPoreTortu4Gas(jroots,JP,JY,JX));  RootPoreTortu4Gas=0._r8
  allocate(RootNodulNonstElms_rpvr(NumPlantChemElms,JZ,JP,JY,JX));RootNodulNonstElms_rpvr=0._r8
  allocate(RootLenPerPlant_pvr(jroots,JZ,JP,JY,JX));RootLenPerPlant_pvr=0._r8
  allocate(Root1stLen_rpvr(jroots,JZ,MaxNumRootAxes,JP,JY,JX));Root1stLen_rpvr=0._r8
  allocate(Root2ndLen_rpvr(jroots,JZ,MaxNumRootAxes,JP,JY,JX));Root2ndLen_rpvr=0._r8
  allocate(RootLenDensPerPlant_pvr(jroots,JZ,JP,JY,JX));RootLenDensPerPlant_pvr=0._r8
  allocate(Root1stXNumL_pvr(jroots,JZ,JP,JY,JX));Root1stXNumL_pvr=0._r8
  allocate(Root2ndXNum_pvr(jroots,JZ,JP,JY,JX));Root2ndXNum_pvr=0._r8
  allocate(Root2ndXNum_rpvr(jroots,JZ,MaxNumRootAxes,JP,JY,JX));Root2ndXNum_rpvr=0._r8
  allocate(Root2ndAveLen_pvr(jroots,JZ,JP,JY,JX));Root2ndAveLen_pvr=0._r8
  allocate(RootAreaPerPlant_pvr(jroots,JZ,JP,JY,JX));RootAreaPerPlant_pvr=0._r8
  allocate(RootVH2O_pvr(jroots,JZ,JP,JY,JX));RootVH2O_pvr=0._r8
  allocate(Root1stRadius_pvr(jroots,JZ,JP,JY,JX));Root1stRadius_pvr=0._r8
  allocate(RootPoreVol_pvr(jroots,JZ,JP,JY,JX));RootPoreVol_pvr=0._r8
  allocate(Root1stDepz_pft(jroots,MaxNumRootAxes,JP,JY,JX));Root1stDepz_pft=0._r8
  allocate(Root2ndRadius_pvr(jroots,JZ,JP,JY,JX));Root2ndRadius_pvr=0._r8
  allocate(Root1stSpecLen_pft(jroots,JP,JY,JX)); Root1stSpecLen_pft=0._r8
  allocate(Root2ndSpecLen_pft(jroots,JP,JY,JX)); Root2ndSpecLen_pft=0._r8
  allocate(AllPlantRootH2OLoss_vr(jroots,JZ,JP,JY,JX));AllPlantRootH2OLoss_vr=0._r8
  allocate(PSIRoot_pvr(jroots,JZ,JP,JY,JX));PSIRoot_pvr=0._r8
  allocate(PSIRootOSMO_vr(jroots,JZ,JP,JY,JX));PSIRootOSMO_vr=0._r8
  allocate(PSIRootTurg_vr(jroots,JZ,JP,JY,JX));PSIRootTurg_vr=0._r8
  allocate(trcg_rootml_pvr(idg_beg:idg_NH3,jroots,JZ,JP,JY,JX)); trcg_rootml_pvr =0._r8
  allocate(trcs_rootml_pvr(idg_beg:idg_NH3,jroots,JZ,JP,JY,JX)); trcs_rootml_pvr =0._r8
  allocate(TRootGasLossDisturb_col(idg_beg:idg_NH3,JY,JX));TRootGasLossDisturb_col=0._r8
  allocate(RootBiomCPerPlant_pft(JP,JY,JX));    RootBiomCPerPlant_pft=0._r8
  allocate(RootElms_pft(NumPlantChemElms,JP,JY,JX)); RootElms_pft=0._r8
  allocate(RootStrutElms_pft(NumPlantChemElms,JP,JY,JX));   RootStrutElms_pft=0._r8
  allocate(RootProteinC_pvr(jroots,JZ,JP,JY,JX));RootProteinC_pvr=0._r8
  allocate(RootMyco1stStrutElms_rpvr(NumPlantChemElms,jroots,JZ,MaxNumRootAxes,JP,JY,JX));RootMyco1stStrutElms_rpvr=0._r8
  allocate(RootMyco2ndStrutElms_rpvr(NumPlantChemElms,jroots,JZ,MaxNumRootAxes,JP,JY,JX));RootMyco2ndStrutElms_rpvr=0._r8
  allocate( PopuRootMycoC_pvr(jroots,JZ,JP,JY,JX)); PopuRootMycoC_pvr=0._r8
  allocate(RootNodulStrutElms_rpvr(NumPlantChemElms,JZ,JP,JY,JX)); RootNodulStrutElms_rpvr=0._r8
  allocate(NodulStrutElms_pft(NumPlantChemElms,JP,JY,JX));  NodulStrutElms_pft=0._r8
  allocate(RootMycoActiveBiomC_pvr(jroots,JZ,JP,JY,JX));RootMycoActiveBiomC_pvr=0._r8
  allocate(RootMycoNonstElms_rpvr(NumPlantChemElms,jroots,JZ,JP,JY,JX)); RootMycoNonstElms_rpvr=0._r8
  allocate(RootNonstructElmConc_rpvr(NumPlantChemElms,jroots,JZ,JP,JY,JX));RootNonstructElmConc_rpvr=0._r8
  allocate(RootMyco1stElm_raxs(NumPlantChemElms,jroots,MaxNumRootAxes,JP,JY,JX));RootMyco1stElm_raxs=0._r8
  allocate(RootProteinConc_rpvr(jroots,JZ,JP,JY,JX));RootProteinConc_rpvr=0._r8
  end subroutine InitRootData

!----------------------------------------------------------------------
  subroutine DestructRootData
  use abortutils, only : destroy
  implicit none

  call destroy(RootGasConductance_pvr)
  call destroy(RootMassElm_vr)
  call destroy(NumRootAxes_pft)
  call destroy(NIXBotRootLayer_rpft)
  call destroy(iPlantRootState_pft)
  call destroy(NIXBotRootLayer_pft)
  call destroy(MaxSoiL4Root_pft)
  call destroy(RootElmsbeg_pft)
  call destroy(RootBiomGrosYld_pft)
  call destroy(MinNonstC2InitRoot_pft)
  call destroy(RootFracRemobilizableBiom)
  call destroy(RootVolPerMassC_pft)
  call destroy(Root1stMaxRadius1_pft)
  call destroy(Root2ndMaxRadius1_pft)
  call destroy(Root1stXSecArea_pft)
  call destroy(Root2ndXSecArea_pft)
  call destroy(fTgrowRootP_vr)
  call destroy(RootrNC_pft)
  call destroy(RootrPC_pft)
  call destroy(RootPorosity_pft)
  call destroy(RootRadialResist_pft)
  call destroy(RootAxialResist_pft)
  call destroy(ShutRutNonstElmntConducts_pft)
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
  call destroy(CNRTS_pft)
  call destroy(CPRTS_pft)
  call destroy(RootMycoNonstElms_pft)
  call destroy(Root1stMaxRadius_pft)
  call destroy(Root2ndMaxRadius_pft)
  call destroy(RootBranchFreq_pft)
  call destroy(RootPoreTortu4Gas)
  call destroy(RootNodulNonstElms_rpvr)
  call destroy(RootLenPerPlant_pvr)
  call destroy(Root1stLen_rpvr)
  call destroy(Root2ndLen_rpvr)
  call destroy(RootLenDensPerPlant_pvr)
  call destroy(Root1stXNumL_pvr)
  call destroy(Root2ndXNum_pvr)
  call destroy(Root2ndXNum_rpvr)
  call destroy(Root2ndAveLen_pvr)
  call destroy(RootAreaPerPlant_pvr)
  call destroy(RootVH2O_pvr)
  call destroy(Root1stRadius_pvr)
  call destroy(RootPoreVol_pvr)
  call destroy(Root1stDepz_pft)
  call destroy(Root2ndRadius_pvr)
  call destroy(Root1stSpecLen_pft)
  call destroy(Root2ndSpecLen_pft)
  call destroy(RootNodulElms_pft)
  call destroy(AllPlantRootH2OLoss_vr)
  call destroy(PSIRoot_pvr)
  call destroy(PSIRootOSMO_vr)
  call destroy(PSIRootTurg_vr)
  call destroy(trcg_rootml_pvr)
  call destroy(trcs_rootml_pvr)
  call destroy(TRootGasLossDisturb_col)
  call destroy(RootBiomCPerPlant_pft)
  call destroy(RootElms_pft)
  call destroy(RootStrutElms_pft)
  call destroy(RootProteinC_pvr)
  call destroy(RootMyco1stStrutElms_rpvr)
  call destroy(RootMyco2ndStrutElms_rpvr)
  call destroy( PopuRootMycoC_pvr)
  call destroy(RootNodulStrutElms_rpvr)
  call destroy(NodulStrutElms_pft)
  call destroy(RootMycoActiveBiomC_pvr)
  call destroy(RootMycoNonstElms_rpvr)
  call destroy(RootNonstructElmConc_rpvr)
  call destroy(RootMyco1stElm_raxs)
  call destroy(RootProteinConc_rpvr)
  end subroutine DestructRootData

end module RootDataType
