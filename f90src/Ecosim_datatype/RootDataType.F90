module RootDataType

!
!!
! data types of plant characteristics
  use data_kind_mod, only : sp => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  integer,target,allocatable ::  NumPrimeRootAxes_pft(:,:,:)                          !root primary axis number, [-]
  integer,target,allocatable ::  NRoot1stTipLay_raxes(:,:,:,:)                   !maximum soil layer number for root axes, [-]
  integer,target,allocatable ::  iPlantRootState_pft(:,:,:)                      !flag to detect root system death , [-]
  integer,target,allocatable ::  NMaxRootBotLayer_pft(:,:,:)                      !maximum soil layer number for all root axes, [-]
  integer,target,allocatable ::  MaxSoiL4Root_pft(:,:,:)                         !maximum soil layer number for all root axes, [-]
  real(sp),target,allocatable :: RootElmsbeg_pft(:,:,:,:)                        !root biomass per pft
  real(sp),target,allocatable ::  RootBiomGrosYld_pft(:,:,:)                     !root growth yield, [g g-1]
  real(sp),target,allocatable ::  MinNonstC2InitRoot_pft(:,:,:)                  !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(sp),target,allocatable ::  RootProteinCMax_pft(:,:,:)                     !reference root protein N, [gN g-1]
  real(sp),target,allocatable ::  FineRootVolPerMassC_pft(:,:,:,:)               !fine root volume:mass ratio, [m3 g-1]
  real(r8),target,allocatable ::  CoarseRootVolPerMassC_pft(:,:,:)               !coarse root volume:mass ratio,[m3 g-1]
  real(sp),target,allocatable ::  Root1stMaxRadius1_pft(:,:,:,:)                 !root diameter primary axes, [m]
  real(sp),target,allocatable ::  Root2ndMaxRadius1_pft(:,:,:,:)                 !root diameter secondary axes, [m]
  real(sp),target,allocatable ::  Root1stXSecArea_pft(:,:,:,:)                   !root cross-sectional area primary axes, [m2]
  real(sp),target,allocatable ::  Root2ndXSecArea_pft(:,:,:,:)                   !root  cross-sectional area  secondary axes, [m2]
  real(sp),target,allocatable ::  fTgrowRootP_vr(:,:,:,:)                        !root layer temperature growth functiom, [-]
  real(sp),target,allocatable ::  rNCRoot_pft(:,:,:)                             !root N:C ratio, [g g-1]
  real(sp),target,allocatable ::  rPCRootr_pft(:,:,:)                             !root P:C ratio, [g g-1]
  real(sp),target,allocatable ::  RootPorosity_pft(:,:,:,:)                      !root porosity, [m3 m-3]
  real(sp),target,allocatable ::  RootRadialResist_pft(:,:,:,:)                  !root radial resistivity, [MPa h m-2]
  real(sp),target,allocatable ::  RootAxialResist_pft(:,:,:,:)                   !root axial resistivity, [MPa h m-4]
  real(sp),target,allocatable ::  ShootRootNonstElmConduts_pft(:,:,:)           !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(sp),target,allocatable ::  VmaxNH4Root_pft(:,:,:,:)                       !maximum root NH4 uptake rate, [g m-2 h-1]
  real(sp),target,allocatable ::  KmNH4Root_pft(:,:,:,:)                         !Km for root NH4 uptake, [g m-3]
  real(sp),target,allocatable ::  CMinNH4Root_pft(:,:,:,:)                       !minimum NH4 concentration for root NH4 uptake, [g m-3]
  real(sp),target,allocatable ::  VmaxNO3Root_pft(:,:,:,:)                       !maximum root NO3 uptake rate, [g m-2 h-1]
  real(sp),target,allocatable ::  KmNO3Root_pft(:,:,:,:)                         !Km for root NO3 uptake, [g m-3]
  real(sp),target,allocatable ::  CminNO3Root_pft(:,:,:,:)                       !minimum NO3 concentration for root NH4 uptake, [g m-3]
  real(sp),target,allocatable ::  VmaxPO4Root_pft(:,:,:,:)                       !maximum root PO4 uptake rate, [g m-2 h-1]
  real(sp),target,allocatable ::  KmPO4Root_pft(:,:,:,:)                         !Km for root PO4 uptake, [g m-3]
  real(sp),target,allocatable ::  CMinPO4Root_pft(:,:,:,:)                       !minimum PO4 concentration for root NH4 uptake, [g m-3]
  real(sp),target,allocatable ::  RootRaidus_rpft(:,:,:,:)                       !root internal radius, [m]
  real(sp),target,allocatable ::  CNRTS_pft(:,:,:)                               !root N:C ratio x root growth yield, [-]
  real(sp),target,allocatable ::  CPRTS_pft(:,:,:)                               !root P:C ratio x root growth yield, [-]
  real(sp),target,allocatable ::  RootMycoNonstElms_pft(:,:,:,:,:)               !non-structural chemical element in roots,[g d-2]
  real(sp),target,allocatable ::  Root1stMaxRadius_pft(:,:,:,:)                  !maximum radius of primary roots, [m]
  real(sp),target,allocatable ::  RootMatureAge_pft(:,:,:)                       !Root age to trigger secondary growth, [h]
  real(sp),target,allocatable ::  Root2ndMaxRadius_pft(:,:,:,:)                  !maximum radius of secondary roots, [m]
  real(sp),target,allocatable ::  RootBranchFreq_pft(:,:,:)                      !root brancing frequency, [m-1]
  real(sp),target,allocatable ::  RootPoreTortu4Gas_pft(:,:,:,:)                 !root tortuosity to calculate root gaseous diffusivity, [-]
  real(sp),target,allocatable ::  RootNodulNonstElms_rpvr(:,:,:,:,:)             !root  layer nonstructural element, [g d-2]
  real(sp),target,allocatable ::  RootTotLenPerPlant_pvr(:,:,:,:,:)              !root layer length per plant, including root hair [m p-1]
  real(sp),target,allocatable ::  RootLenPerPlant_pvr(:,:,:,:,:)                 !root layer length per plant, excluding root hair [m p-1]       
  real(sp),target,allocatable ::  Root1stLen_rpvr(:,:,:,:,:)                   !root layer length primary axes, [m d-2]
  real(sp),target,allocatable ::  RootAge_rpvr(:,:,:,:,:)                        !root age, [h]
  real(sp),target,allocatable ::  Root2ndLen_rpvr(:,:,:,:,:,:)                   !root layer length secondary axes, [m d-2]
  real(sp),target,allocatable ::  RootLenDensPerPlant_pvr(:,:,:,:,:)             !root length density in soil layers, [m m-3]
  real(sp),target,allocatable ::  Root1stXNumL_rpvr(:,:,:,:)                    !root layer number primary axes, [d-2]
  real(sp),target,allocatable ::  Root2ndXNumL_rpvr(:,:,:,:,:)                     !root layer number axes, [d-2]
  real(sp),target,allocatable ::  Root2ndXNum_rpvr(:,:,:,:,:,:)                  !root layer number secondary axes, [d-2]
  real(sp),target,allocatable ::  Root2ndEffLen4uptk_rpvr(:,:,:,:,:)              !Layer effective root length four resource uptake, [m]
  real(sp),target,allocatable ::  RootSAreaPerPlant_pvr(:,:,:,:,:)                !root layer area per plant, [m p-1]
  real(sp),target,allocatable  :: RootArea1stPP_pvr(:,:,:,:,:)                   !layer 1st root area per plant, [m2 plant-1]
  real(sp),target,allocatable  :: RootArea2ndPP_pvr(:,:,:,:,:)                  !layer 2nd root area per plant, [m2 plant-1]
  real(sp),target,allocatable ::  Root1stStructE_buf(:,:,:,:,:)                 !buffer structural element for root secondary growth calculation, [g d-2]
  real(sp),target,allocatable ::  RootVH2O_pvr(:,:,:,:,:)                        !root layer volume water, [m2 d-2]
  real(sp),target,allocatable ::  Root1stRadius_pvr(:,:,:,:,:)                   !root layer diameter primary axes, [m]
  real(sp),target,allocatable ::  Root1stRadius_rpvr(:,:,:,:,:)                  !root layer diameter for each primary axes, [m]
  real(sp),target,allocatable ::  RootCRRadius0_rpvr(:,:,:,:,:)                     !initial radius for root that may undergo secondary growth, [m]
  real(sp),target,allocatable ::  RootPoreVol_pvr(:,:,:,:,:)                     !root layer volume air, [m2 d-2]
  real(sp),target,allocatable ::  Root1stDepz_raxes(:,:,:,:)                     !root layer depth, [m]
  real(sp),target,allocatable ::  Root2ndRadius_rpvr(:,:,:,:,:)                   !root layer diameter secondary axes, [m ]
  real(sp),target,allocatable ::  Root1stSpecLen_pft(:,:,:,:)                    !specific root length primary axes, [m g-1]
  real(sp),target,allocatable ::  Root2ndSpecLen_pft(:,:,:,:)                    !specific root length secondary axes, [m g-1]
  real(sp),target,allocatable ::  RPlantRootH2OUptk_pvr(:,:,:,:,:)               !root water uptake, [m3 d-2 h-1]
  real(sp),target,allocatable ::  RootH2OUptkStress_pvr(:,:,:,:,:)               !root water uptake stress indicated by rate,  [m3 d-2 h-1] 
  real(sp),target,allocatable ::  PSIRoot_pvr(:,:,:,:,:)                         !root total water potential , [Mpa]
  real(sp),target,allocatable ::  PSIRootOSMO_vr(:,:,:,:,:)                      !root osmotic water potential , [Mpa]
  real(sp),target,allocatable ::  PSIRootTurg_vr(:,:,:,:,:)                      !root turgor water potential , [Mpa]
  real(sp),target,allocatable ::  RootSinkWeight_pvr(:,:,:,:)                     !Root nonst element sink profile, [d-2]
  real(sp),target,allocatable ::  trcg_rootml_pvr(:,:,:,:,:,:)                   !root gaseous tracer content [g d-2]
  real(sp),target,allocatable ::  trcs_rootml_pvr(:,:,:,:,:,:)                   !root dissolved gaseous tracer content [g d-2]
  real(sp),target,allocatable ::  TRootGasLossDisturb_col(:,:,:)                 !total root gas content, [g d-2]
  real(sp),target,allocatable ::  RootBiomCPerPlant_pft(:,:,:)                   !root C per plant, [g p-1]
  real(sp),target,allocatable ::  RootElms_pft(:,:,:,:)                          !plant root element, [g d-2]
  real(sp),target,allocatable ::  RootNoduleElms_pft(:,:,:,:)                    !plant root nodule element mass,                     [g d-2]
  real(sp),target,allocatable ::  RootStrutElms_pft(:,:,:,:)                     !plant root structural element, [g d-2]
  real(sp),target,allocatable ::  RootProteinC_pvr(:,:,:,:,:)                    !root layer protein C, [g d-2]
  real(sp),target,allocatable ::  RootMyco1stStrutElms_rpvr(:,:,:,:,:,:)       !root layer element primary axes, [g d-2]
  real(sp),target,allocatable ::  RootMyco2ndStrutElms_rpvr(:,:,:,:,:,:,:)       !root layer element secondary axes, [g d-2]
  real(sp),target,allocatable ::   PopuRootMycoC_pvr(:,:,:,:,:)                  !root layer C, [g d-2]
  real(sp),target,allocatable ::  RootNodulStrutElms_rpvr(:,:,:,:,:)             !root layer nodule element, [g d-2]
  real(sp),target,allocatable ::  RootMycoActiveBiomC_pvr(:,:,:,:,:)             !root layer structural C, [g d-2]
  real(sp),target,allocatable ::   RootMycoNonstElms_rpvr(:,:,:,:,:,:)           !root  layer nonstructural element, [g d-2]
  real(sp),target,allocatable ::  RootNonstructElmConc_rpvr(:,:,:,:,:,:)         !root  layer nonstructural element concentration, [g g-1]
  real(sp),target,allocatable ::  RootMyco1stElm_raxs(:,:,:,:,:)               !root C primary axes, [g d-2]
  real(sp),target,allocatable ::  RootProteinConc_rpvr(:,:,:,:,:)                !root layer protein C concentration, [g g-1]
  real(sp),target,allocatable :: RootMycoMassElm_vr(:,:,:,:,:)                   !root chemical element mass in soil layer, [g d-2]
  real(sp),target,allocatable :: RootMycoMassElm_pvr(:,:,:,:,:,:)                  !pft root chemical elemental mass in soil layer, [g d-2]
  real(sp),target,allocatable :: RootGasConductance_rpvr(:,:,:,:,:,:)            !Root Conductance for gas uptake, [m3 d-2 h-1]
  real(sp),target,allocatable :: Nutruptk_fClim_rpvr(:,:,:,:,:)                  !Carbon limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(sp),target,allocatable :: Nutruptk_fNlim_rpvr(:,:,:,:,:)                  !Nitrogen limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(sp),target,allocatable :: Nutruptk_fPlim_rpvr(:,:,:,:,:)                  !Phosphorus limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(sp),target,allocatable :: Nutruptk_fProtC_rpvr(:,:,:,:,:)                 !transporter scalar indicated by protein for root nutrient uptake, greater value greater capacity, [-] 
  real(sp),target,allocatable :: RootMaintDef_CO2_pvr(:,:,:,:,:)                 !plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]  
  real(sp),target,allocatable :: ROOTNLim_rpvr(:,:,:,:,:)                        !root N-limitation, 0->1 weaker limitation, [-]     
  real(sp),target,allocatable :: ROOTPLim_rpvr(:,:,:,:,:)                        !root P-limitation, 0->1 weaker limitation, [-]         
  real(sp),target,allocatable :: RootResist4H2O_pvr(:,:,:,:,:)                   !total root (axial+radial) resistance for water uptake, [MPa-1 h-1]     
  real(sp),target,allocatable :: RootRadialKond2H2O_pvr(:,:,:,:,:)               !radial root conductance for water uptake, [m3 H2O h-1 MPa-1]
  real(sp),target,allocatable :: RootAxialKond2H2O_pvr(:,:,:,:,:)                !axial root conductance for water uptake, [m3 H2O h-1 MPa-1]
  real(sp),target, allocatable :: RootSegAges_raxes(:,:,:,:,:)                   !age of different active root segments, [h]
  integer ,target, allocatable :: NActiveRootSegs_raxes(:,:,:,:)                 !number of active root segments, [-]
  real(sp),target, allocatable :: RootSegBaseDepth_raxes(:,:,:,:)                !base depth of different root axes, [m]
  real(sp),target, allocatable :: RootSeglengths_raxes(:,:,:,:,:)                !root length in each segment, [m]
  integer ,target, allocatable :: IndRootSegBase_raxes(:,:,:,:)                  !Index for base segment under tracking, [-]
  integer ,target, allocatable :: IndRootSegTip_raxes(:,:,:,:)                   !Index for tip segment under tracking, [-]  

!----------------------------------------------------------------------

contains
  subroutine InitRootData(jroots)

  implicit none
  integer, intent(in) :: jroots

  allocate(ROOTNLim_rpvr(jroots,JZ,JP,JY,JX)); ROOTNLim_rpvr=0._sp
  allocate(ROOTPLim_rpvr(jroots,JZ,JP,JY,JX)); ROOTPLim_rpvr=0._sp
  allocate(RootMaintDef_CO2_pvr(jroots,JZ,JP,JY,JX)); RootMaintDef_CO2_pvr=0._sp
  allocate(Nutruptk_fClim_rpvr(jroots,JZ,JP,JY,JX));Nutruptk_fClim_rpvr=0._sp
  allocate(Nutruptk_fNlim_rpvr(jroots,JZ,JP,JY,JX));Nutruptk_fNlim_rpvr=0._sp
  allocate(Nutruptk_fPlim_rpvr(jroots,JZ,JP,JY,JX));Nutruptk_fPlim_rpvr=0._sp
  allocate(Nutruptk_fProtC_rpvr(jroots,JZ,JP,JY,JX));Nutruptk_fProtC_rpvr=0._sp
  allocate(RootMycoMassElm_vr(NumPlantChemElms,jroots,JZ,JY,JX)); RootMycoMassElm_vr =0._sp
  allocate(RootMycoMassElm_pvr(NumPlantChemElms,jroots,JZ,JP,JY,JX)); RootMycoMassElm_pvr=0._sp
  allocate(NumPrimeRootAxes_pft(JP,JY,JX));      NumPrimeRootAxes_pft=0
  allocate(NRoot1stTipLay_raxes(MaxNumRootAxes,JP,JY,JX));  NRoot1stTipLay_raxes=1  !set to one to avoid numerical failure
  allocate(iPlantRootState_pft(JP,JY,JX));    iPlantRootState_pft=iDead
  allocate(NMaxRootBotLayer_pft(JP,JY,JX));      NMaxRootBotLayer_pft=0
  allocate(MaxSoiL4Root_pft(JP,JY,JX));       MaxSoiL4Root_pft=0
  allocate(RootElmsbeg_pft(NumPlantChemElms,JP,JY,JX)); RootElmsbeg_pft=0._sp
  allocate(RootBiomGrosYld_pft(JP,JY,JX));     RootBiomGrosYld_pft=0._sp
  allocate(MinNonstC2InitRoot_pft(JP,JY,JX));       MinNonstC2InitRoot_pft=0._sp
  allocate(RootProteinCMax_pft(JP,JY,JX));    RootProteinCMax_pft=0._sp
  allocate(FineRootVolPerMassC_pft(jroots,JP,JY,JX));   FineRootVolPerMassC_pft=0._sp
  allocate(CoarseRootVolPerMassC_pft(JP,JY,JX)); CoarseRootVolPerMassC_pft = 0._sp
  allocate(Root1stMaxRadius1_pft(jroots,JP,JY,JX)); Root1stMaxRadius1_pft=0._sp
  allocate(Root2ndMaxRadius1_pft(jroots,JP,JY,JX)); Root2ndMaxRadius1_pft=0._sp
  allocate(Root1stXSecArea_pft(jroots,JP,JY,JX)); Root1stXSecArea_pft=0._sp
  allocate(Root2ndXSecArea_pft(jroots,JP,JY,JX)); Root2ndXSecArea_pft=0._sp
  allocate(fTgrowRootP_vr(JZ,JP,JY,JX));  fTgrowRootP_vr=0._sp
  allocate(rNCRoot_pft(JP,JY,JX));     rNCRoot_pft=0._sp
  allocate(rPCRootr_pft(JP,JY,JX));     rPCRootr_pft=0._sp
  allocate(RootGasConductance_rpvr(idg_beg:idg_NH3,jroots,JZ,JP,JY,JX));RootGasConductance_rpvr=0._sp
  allocate(RootPorosity_pft(jroots,JP,JY,JX));   RootPorosity_pft=0._sp
  allocate(RootRadialResist_pft(jroots,JP,JY,JX));   RootRadialResist_pft=0._sp
  allocate(RootAxialResist_pft(jroots,JP,JY,JX));   RootAxialResist_pft=0._sp
  allocate(ShootRootNonstElmConduts_pft(JP,JY,JX));    ShootRootNonstElmConduts_pft=0._sp
  allocate(VmaxNH4Root_pft(jroots,JP,JY,JX)); VmaxNH4Root_pft=0._sp
  allocate(KmNH4Root_pft(jroots,JP,JY,JX)); KmNH4Root_pft=0._sp
  allocate(CMinNH4Root_pft(jroots,JP,JY,JX)); CMinNH4Root_pft=0._sp
  allocate(VmaxNO3Root_pft(jroots,JP,JY,JX)); VmaxNO3Root_pft=0._sp
  allocate(KmNO3Root_pft(jroots,JP,JY,JX)); KmNO3Root_pft=0._sp
  allocate(CminNO3Root_pft(jroots,JP,JY,JX)); CminNO3Root_pft=0._sp
  allocate(VmaxPO4Root_pft(jroots,JP,JY,JX)); VmaxPO4Root_pft=0._sp
  allocate(KmPO4Root_pft(jroots,JP,JY,JX)); KmPO4Root_pft=0._sp
  allocate(CMinPO4Root_pft(jroots,JP,JY,JX)); CMinPO4Root_pft=0._sp
  allocate(RootRaidus_rpft(jroots,JP,JY,JX));  RootRaidus_rpft=0._sp
  allocate(CNRTS_pft(JP,JY,JX));    CNRTS_pft=0._sp
  allocate(CPRTS_pft(JP,JY,JX));    CPRTS_pft=0._sp
  allocate(RootMatureAge_pft(JP,JY,JX)); RootMatureAge_pft=0._sp
  allocate(RootMycoNonstElms_pft(NumPlantChemElms,jroots,JP,JY,JX));RootMycoNonstElms_pft=0._sp
  allocate(Root1stMaxRadius_pft(jroots,JP,JY,JX)); Root1stMaxRadius_pft=0._sp
  allocate(Root2ndMaxRadius_pft(jroots,JP,JY,JX)); Root2ndMaxRadius_pft=0._sp
  allocate(RootBranchFreq_pft(JP,JY,JX));     RootBranchFreq_pft=0._sp
  allocate(RootPoreTortu4Gas_pft(jroots,JP,JY,JX));  RootPoreTortu4Gas_pft=0._sp
  allocate(RootNodulNonstElms_rpvr(NumPlantChemElms,JZ,JP,JY,JX));RootNodulNonstElms_rpvr=0._sp
  allocate(RootTotLenPerPlant_pvr(jroots,JZ,JP,JY,JX));RootTotLenPerPlant_pvr=0._sp
  allocate(RootLenPerPlant_pvr(jroots,JZ,JP,JY,JX));RootLenPerPlant_pvr=0._sp
  allocate(Root1stLen_rpvr(JZ,MaxNumRootAxes,JP,JY,JX));Root1stLen_rpvr=0._sp
  allocate(RootAge_rpvr(JZ,MaxNumRootAxes,JP,JY,JX)); RootAge_rpvr=0._sp
  allocate(Root2ndLen_rpvr(jroots,JZ,MaxNumRootAxes,JP,JY,JX));Root2ndLen_rpvr=0._sp
  allocate(RootLenDensPerPlant_pvr(jroots,JZ,JP,JY,JX));RootLenDensPerPlant_pvr=0._sp
  allocate(Root1stXNumL_rpvr(JZ,JP,JY,JX));Root1stXNumL_rpvr=0._sp
  allocate(Root2ndXNumL_rpvr(jroots,JZ,JP,JY,JX));Root2ndXNumL_rpvr=0._sp
  allocate(Root2ndXNum_rpvr(jroots,JZ,MaxNumRootAxes,JP,JY,JX));Root2ndXNum_rpvr=0._sp
  allocate(Root2ndEffLen4uptk_rpvr(jroots,JZ,JP,JY,JX));Root2ndEffLen4uptk_rpvr=0._sp
  allocate(RootSAreaPerPlant_pvr(jroots,JZ,JP,JY,JX));RootSAreaPerPlant_pvr=0._sp
  allocate(RootArea1stPP_pvr(jroots,JZ,JP,JY,JX));RootArea1stPP_pvr=0._sp
  allocate(RootArea2ndPP_pvr(jroots,JZ,JP,JY,JX));RootArea2ndPP_pvr=0._sp
  allocate(RootVH2O_pvr(jroots,JZ,JP,JY,JX));RootVH2O_pvr=0._sp
  allocate(Root1stRadius_pvr(jroots,JZ,JP,JY,JX));Root1stRadius_pvr=0._sp
  allocate(Root1stRadius_rpvr(JZ,MaxNumRootAxes,JP,JY,JX)); Root1stRadius_rpvr=0._sp
  allocate(RootCRRadius0_rpvr(JZ,MaxNumRootAxes,JP,JY,JX)); RootCRRadius0_rpvr=0._sp
  allocate(RootPoreVol_pvr(jroots,JZ,JP,JY,JX));RootPoreVol_pvr=0._sp
  allocate(Root1stDepz_raxes(MaxNumRootAxes,JP,JY,JX));Root1stDepz_raxes=0._sp
  allocate(Root2ndRadius_rpvr(jroots,JZ,JP,JY,JX));Root2ndRadius_rpvr=0._sp
  allocate(Root1stSpecLen_pft(jroots,JP,JY,JX)); Root1stSpecLen_pft=0._sp
  allocate(Root2ndSpecLen_pft(jroots,JP,JY,JX)); Root2ndSpecLen_pft=0._sp
  allocate(RPlantRootH2OUptk_pvr(jroots,JZ,JP,JY,JX));RPlantRootH2OUptk_pvr=0._sp
  allocate(RootH2OUptkStress_pvr(jroots,JZ,JP,JY,JX)); RootH2OUptkStress_pvr=0._sp
  allocate(PSIRoot_pvr(jroots,JZ,JP,JY,JX));PSIRoot_pvr=0._sp
  allocate(PSIRootOSMO_vr(jroots,JZ,JP,JY,JX));PSIRootOSMO_vr=0._sp
  allocate(PSIRootTurg_vr(jroots,JZ,JP,JY,JX));PSIRootTurg_vr=0._sp
  allocate(RootSinkWeight_pvr(JZ,JP,JY,JX)); RootSinkWeight_pvr=0._sp
  allocate(trcg_rootml_pvr(idg_beg:idg_NH3,jroots,JZ,JP,JY,JX)); trcg_rootml_pvr =0._sp
  allocate(trcs_rootml_pvr(idg_beg:idg_NH3,jroots,JZ,JP,JY,JX)); trcs_rootml_pvr =0._sp
  allocate(TRootGasLossDisturb_col(idg_beg:idg_NH3,JY,JX));TRootGasLossDisturb_col=0._sp
  allocate(RootBiomCPerPlant_pft(JP,JY,JX));    RootBiomCPerPlant_pft=0._sp
  allocate(RootElms_pft(NumPlantChemElms,JP,JY,JX)); RootElms_pft=0._sp
  allocate(RootNoduleElms_pft(NumPlantChemElms,JP,JY,JX)); RootNoduleElms_pft=0._sp
  allocate(RootStrutElms_pft(NumPlantChemElms,JP,JY,JX));   RootStrutElms_pft=0._sp
  allocate(RootProteinC_pvr(jroots,JZ,JP,JY,JX));RootProteinC_pvr=0._sp
  allocate(RootMyco1stStrutElms_rpvr(NumPlantChemElms,JZ,MaxNumRootAxes,JP,JY,JX));RootMyco1stStrutElms_rpvr=0._sp
  allocate(Root1stStructE_buf(NumPlantChemElms,MaxNumRootAxes,JP,JY,JX)); Root1stStructE_buf=0._sp
  allocate(RootMyco2ndStrutElms_rpvr(NumPlantChemElms,jroots,JZ,MaxNumRootAxes,JP,JY,JX));RootMyco2ndStrutElms_rpvr=0._sp
  allocate( PopuRootMycoC_pvr(jroots,JZ,JP,JY,JX)); PopuRootMycoC_pvr=0._sp
  allocate(RootNodulStrutElms_rpvr(NumPlantChemElms,JZ,JP,JY,JX)); RootNodulStrutElms_rpvr=0._sp
  allocate(RootMycoActiveBiomC_pvr(jroots,JZ,JP,JY,JX));RootMycoActiveBiomC_pvr=0._sp
  allocate(RootMycoNonstElms_rpvr(NumPlantChemElms,jroots,JZ,JP,JY,JX)); RootMycoNonstElms_rpvr=0._sp
  allocate(RootNonstructElmConc_rpvr(NumPlantChemElms,jroots,JZ,JP,JY,JX));RootNonstructElmConc_rpvr=0._sp
  allocate(RootMyco1stElm_raxs(NumPlantChemElms,MaxNumRootAxes,JP,JY,JX));RootMyco1stElm_raxs=0._sp
  allocate(RootProteinConc_rpvr(jroots,JZ,JP,JY,JX));RootProteinConc_rpvr=0._sp
  allocate(RootResist4H2O_pvr(jroots,JZ,JP,JY,JX)); RootResist4H2O_pvr=0._sp
  allocate(RootRadialKond2H2O_pvr(jroots,JZ,JP,JY,JX));RootRadialKond2H2O_pvr=0._sp
  allocate(RootAxialKond2H2O_pvr(jroots,JZ,JP,JY,JX));RootAxialKond2H2O_pvr=0._sp

  allocate(RootSegAges_raxes(1:NMaxRootSegs,1:MaxNumRootAxes,JP,JY,JX));RootSegAges_raxes=0._sp
  allocate(NActiveRootSegs_raxes(1:MaxNumRootAxes,JP,JY,JX));NActiveRootSegs_raxes=0
  allocate(RootSegBaseDepth_raxes(1:MaxNumRootAxes,JP,JY,JX));RootSegBaseDepth_raxes=0._sp     
  allocate(RootSeglengths_raxes(1:NMaxRootSegs,1:MaxNumRootAxes,JP,JY,JX));RootSeglengths_raxes=0._sp       
  allocate(IndRootSegBase_raxes(1:MaxNumRootAxes,JP,JY,JX));IndRootSegBase_raxes=0
  allocate(IndRootSegTip_raxes(1:MaxNumRootAxes,JP,JY,JX));IndRootSegTip_raxes=0             

  end subroutine InitRootData

!----------------------------------------------------------------------
  subroutine DestructRootData
  use abortutils, only : destroy
  implicit none
  call destroy(RootRadialKond2H2O_pvr)
  call destroy(RootAxialKond2H2O_pvr)
  call destroy(RootResist4H2O_pvr)
  call destroy(ROOTNLim_rpvr)
  call destroy(ROOTPLim_rpvr)  
  call destroy(RootMaintDef_CO2_pvr)
  call destroy(RootGasConductance_rpvr)
  call destroy(RootMycoMassElm_pvr)
  call destroy(RootMycoMassElm_vr)
  call destroy(NumPrimeRootAxes_pft)
  call destroy(NRoot1stTipLay_raxes)
  call destroy(iPlantRootState_pft)
  call destroy(NMaxRootBotLayer_pft)
  call destroy(MaxSoiL4Root_pft)
  call destroy(RootElmsbeg_pft)
  call destroy(RootBiomGrosYld_pft)
  call destroy(MinNonstC2InitRoot_pft)
  call destroy(RootProteinCMax_pft)
  call destroy(FineRootVolPerMassC_pft)
  call destroy(CoarseRootVolPerMassC_pft)
  call destroy(Root1stMaxRadius1_pft)
  call destroy(Root2ndMaxRadius1_pft)
  call destroy(Root1stXSecArea_pft)
  call destroy(Root2ndXSecArea_pft)
  call destroy(fTgrowRootP_vr)
  call destroy(rNCRoot_pft)
  call destroy(rPCRootr_pft)
  call destroy(RootPorosity_pft)
  call destroy(RootRadialResist_pft)
  call destroy(RootAxialResist_pft)
  call destroy(ShootRootNonstElmConduts_pft)
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
  call destroy(RootMatureAge_pft)
  call destroy(RootMycoNonstElms_pft)
  call destroy(Root1stMaxRadius_pft)
  call destroy(Root2ndMaxRadius_pft)
  call destroy(RootBranchFreq_pft)
  call destroy(RootPoreTortu4Gas_pft)
  call destroy(RootNodulNonstElms_rpvr)
  call destroy(RootTotLenPerPlant_pvr)
  call destroy(RootLenPerPlant_pvr)
  call destroy(Root1stLen_rpvr)
  call destroy(RootAge_rpvr)
  call destroy(Root2ndLen_rpvr)
  call destroy(RootLenDensPerPlant_pvr)
  call destroy(Root1stXNumL_rpvr)
  call destroy(Root2ndXNumL_rpvr)
  call destroy(Root2ndXNum_rpvr)
  call destroy(Root2ndEffLen4uptk_rpvr)
  call destroy(RootSAreaPerPlant_pvr)
  call destroy(RootArea2ndPP_pvr)
  call destroy(RootArea1stPP_pvr)  
  call destroy(RootVH2O_pvr)
  call destroy(Root1stRadius_pvr)
  call destroy(Root1stRadius_rpvr)
  call destroy(RootCRRadius0_rpvr)
  call destroy(RootPoreVol_pvr)
  call destroy(Root1stDepz_raxes)
  call destroy(Root2ndRadius_rpvr)
  call destroy(Root1stSpecLen_pft)
  call destroy(Root2ndSpecLen_pft)
  call destroy(RPlantRootH2OUptk_pvr)
  call destroy(RootH2OUptkStress_pvr)
  call destroy(PSIRoot_pvr)
  call destroy(PSIRootOSMO_vr)
  call destroy(PSIRootTurg_vr)
  call destroy(RootSinkWeight_pvr)
  call destroy(trcg_rootml_pvr)
  call destroy(trcs_rootml_pvr)
  call destroy(TRootGasLossDisturb_col)
  call destroy(RootBiomCPerPlant_pft)
  call destroy(RootElms_pft)
  call destroy(RootNoduleElms_pft)
  call destroy(RootStrutElms_pft)
  call destroy(RootProteinC_pvr)
  call destroy(Root1stStructE_buf)
  call destroy(RootMyco1stStrutElms_rpvr)
  call destroy(RootMyco2ndStrutElms_rpvr)
  call destroy( PopuRootMycoC_pvr)
  call destroy(RootNodulStrutElms_rpvr)
  call destroy(RootMycoActiveBiomC_pvr)
  call destroy(RootMycoNonstElms_rpvr)
  call destroy(RootNonstructElmConc_rpvr)
  call destroy(RootMyco1stElm_raxs)
  call destroy(RootProteinConc_rpvr)
  call destroy(RootSegAges_raxes)
  call destroy(NActiveRootSegs_raxes)
  call destroy(RootSegBaseDepth_raxes)
  call destroy(RootSeglengths_raxes)
  call destroy(IndRootSegBase_raxes)
  call destroy(IndRootSegTip_raxes)
  call destroy(Nutruptk_fClim_rpvr)
  call destroy(Nutruptk_fNlim_rpvr)
  call destroy(Nutruptk_fPlim_rpvr)
  call destroy(Nutruptk_fProtC_rpvr)

  end subroutine DestructRootData

end module RootDataType
