module PlantBGCRatesAPIData
  ! Owns the bgcrates API type, instance and lifecycle procedures.
  use PlantAPICommonData
  implicit none
  private
  public :: plt_bgcrate_init, plt_bgcrate_destroy

  type, public :: plant_bgcrate_type
  real(r8) :: Eco_NBP_CumYr_col        !total NBP, [g d-2]
  real(r8) :: NetCO2Flx2Canopy_col     !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8) :: ECO_ER_col               !ecosystem respiration, [g d-2 h-1]
  real(r8) :: Eco_AutoR_CumYr_col      !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: TRootH2Flx_col           !total root H2 flux, [g d-2]
  real(r8) :: Canopy_NEE_col           !total net CO2 fixation, [gC d-2]
  real(r8), pointer :: PTSHTR_pft(:)                       => null()  !root-shoot coupling conductance, [h-1]
  real(r8), pointer :: RootShootExch_pvr(:,:,:)            => null()  !Root-shoot nonstrucal element exchange, [g d-2 h-1]
  real(r8), pointer :: Nutruptk_fClim_rpvr(:,:,:)          => null()  !Carbon limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(r8), pointer :: Nutruptk_fNlim_rpvr(:,:,:)          => null()  !Nitrogen limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(r8), pointer :: Nutruptk_fPlim_rpvr(:,:,:)          => null()  !Phosphorus limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(r8), pointer :: Nutruptk_fProtC_rpvr(:,:,:)         => null()  !transporter scalar indicated by protein for root nutrient uptake, greater value greater capacity, [-]
  real(r8), pointer :: LitrFallStrutElms_col(:)            => null()  !total LitrFall structural element mass,     [g d-2 h-1]
  real(r8), pointer :: NetPrimProduct_pft(:)               => null()  !total net primary productivity,             [gC d-2]
  real(r8), pointer :: NH3Dep2Can_pft(:)                   => null()  !canopy NH3 flux,                            [g d-2 h-1]
  real(r8), pointer :: tRootMycoExud2Soil_vr(:,:,:)        => null()  !total root element exchange,                 [g d-2 h-1]
  real(r8), pointer :: RCO2Nodule_pvr(:,:)                 => null()  !layered root nodule respiration, [gC d-2 h-1]
  real(r8), pointer :: RNFixCO2_pft(:)                     => null()  !CO2 respired by nodules, [gC d-2]
  real(r8), pointer :: RootN2Fix_pvr(:,:)                  => null()  !root N2 fixation,                            [gN d-2 h-1]
  real(r8), pointer :: RootO2_TotSink_pvr(:,:,:)           => null()  !root O2 sink for autotrophic respiraiton,     [gC d-2 h-1]
  real(r8), pointer :: RootO2_TotSink_vr(:)                => null()  !all root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
  real(r8), pointer :: RCanMaintDef_CO2_pft(:)             => null()  !canopy maintenance respiraiton deficit as CO2, [gC d-2 h-1]
  real(r8), pointer :: CanopyRespC_CumYr_pft(:)            => null()  !total autotrophic respiration,               [gC d-2 ]
  real(r8), pointer :: CanopyNLimFactor_brch(:,:)          => null()  !Canopy N-limitation factor, [0->1] weaker limitation,[-]
  real(r8), pointer :: CanopyPLimFactor_brch(:,:)          => null()  !Canopy P-limitation factor, [0->1] weaker limitation,[-]
  real(r8), pointer :: LitrfallElms_pvr(:,:,:,:,:)         => null()  !plant LitrFall element,                      [g d-2 h-1]
  real(r8), pointer :: SSXferElms_pft(:,:)                 => null()  !export flux from the seasonal storage, [g h-1 d-2]
  real(r8), pointer :: SSXfer2ShootElms_pft(:,:)           => null()  !flux export from seasonal storage to shoot, [g h-1 d-2]
  real(r8), pointer :: Xfer2RootsC_pft(:)                  => null()  !carbon transfer from other nonstructural source to roots, [gC d-2 h-1]
  real(r8), pointer :: LitrFallElms_brch(:,:,:)            => null()  !litterfall from the branch, [g d-2 h-1]
  real(r8), pointer :: RootMaintDef_CO2_pvr(:,:,:)         => null()  !plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]
  real(r8), pointer :: REcoO2DmndResp_vr(:)                => null()  !total root + microbial O2 uptake,            [g d-2 h-1]
  real(r8), pointer :: REcoNH4DmndBand_vr(:)               => null()   !total root + microbial NH4 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: REcoH1PO4DmndSoil_vr(:)             => null()   !HPO4 demand in non-band by all microbial,   root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: REcoH1PO4DmndBand_vr(:)             => null()   !HPO4 demand in band by all microbial,       root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: REcoNO3DmndSoil_vr(:)               => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: REcoNH4DmndSoil_vr(:)               => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: REcoNO3DmndBand_vr(:)               => null()   !total root + microbial NO3 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: REcoH2PO4DmndSoil_vr(:)             => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: REcoH2PO4DmndBand_vr(:)             => null()   !total root + microbial PO4 uptake band,     [gP d-2 h-1]
  real(r8), pointer :: RH2PO4EcoDmndSoilPrev_vr(:)         => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: RH2PO4EcoDmndBandPrev_vr(:)         => null()   !total root + microbial PO4 uptake band,     [gP d-2 h-1]
  real(r8), pointer :: RH1PO4EcoDmndSoilPrev_vr(:)         => null()   !HPO4 demand in non-band by all microbial,   root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: RH1PO4EcoDmndBandPrev_vr(:)         => null()   !HPO4 demand in band by all microbial,       root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: RNO3EcoDmndSoilPrev_vr(:)           => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4EcoDmndSoilPrev_vr(:)           => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4EcoDmndBandPrev_vr(:)           => null()   !total root + microbial NH4 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: RNO3EcoDmndBandPrev_vr(:)           => null()   !total root + microbial NO3 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: RGasTranspFlxPrev_vr(:,:)           => null()   !net gaseous flux,                           [g d-2 h-1]
  real(r8), pointer :: RO2AquaSourcePrev_vr(:)             => null()   !net aqueous O2 flux,                        [g d-2 h-1]
  real(r8), pointer :: RO2EcoDmndPrev_vr(:)                => null()   !total root + microbial O2 uptake,           [g d-2 h-1]
  real(r8), pointer :: RootCO2Emis2Root_vr(:)              => null()   !total root CO2 flux,                        [gC d-2 h-1]
  real(r8), pointer :: RootCO2Emis2Root_pvr(:,:)           => null()   !root CO2 flux,                        [gC d-2 h-1]
  real(r8), pointer :: RUptkRootO2_vr(:)                   => null()   !total root internal O2 flux,                [g d-2 h-1]
  real(r8), pointer :: LitrfalStrutElms_vr(:,:,:,:)        => null()   !total LitrFall element,                       [g d-2 h-1]
  real(r8), pointer :: REcoDOMProd_vr(:,:,:)               => null()   !net microbial DOC flux,                      [gC d-2 h-1]
  real(r8), pointer :: CO2NetFix_pft(:)                    => null()   !canopy net CO2 exchange,                     [gC d-2 h-1]
  real(r8), pointer :: GrossCO2Fix_pft(:)                  => null()   !total gross CO2 fixation,                    [gC d-2 ]
  real(r8), pointer :: LitrfallElms_pft(:,:)               => null()   !plant element Litrfall,                      [g d-2 h-1]
  real(r8), pointer :: LitrfallAbvgElms_pft(:,:)           => null()   !aboveground litterfall, [g d-2 h-1]
  real(r8), pointer :: SurfLitrfallElms_pft(:,:)           => null()   !surface litterfall, [g d-2 h-1]
  real(r8), pointer :: LitrfallBlgrElms_pft(:,:)           => null()   !belwground litterfall, [g d-2 h-1]
  real(r8), pointer :: RootGasLossDisturb_pft(:,:)         => null()   !gaseous flux fron root disturbance,           [g d-2 h-1]
  real(r8), pointer :: SurfLitrfalStrutElms_CumYr_pft(:,:) => null()   !total surface LitrFall element,              [g d-2]
  real(r8), pointer :: LitrfalStrutElms_CumYr_pft(:,:)     => null()   !total plant element LitrFall,                [g d-2 ]
  real(r8), pointer :: GrossResp_pft(:)                    => null()   !total plant respiration,                     [gC d-2 ]
  real(r8), pointer :: CanopyGrosRCO2_pft(:)               => null()   !canopy plant+nodule autotrophic respiraiton, [gC d-2]
  real(r8), pointer :: CanopyResp_brch(:,:)                => null()   !canopy respiration for a branch, [gC d-2 h-1]
  real(r8), pointer :: RootAutoCO2_pft(:)                  => null()   !root autotrophic respiraiton, [gC d-2]
  real(r8), pointer :: NodulInfectElms_pft(:,:)            => null()   !nodule infection chemical element mass, [g d-2]
  real(r8), pointer :: NH3Emis_CumYr_pft(:)                => null()   !total canopy NH3 flux,                       [gN d-2 ]
  real(r8), pointer :: PlantN2Fix_CumYr_pft(:)             => null()   !total plant N2 fixation,                     [g d-2 ]
  real(r8), pointer :: ShootRootXferElm_pft(:,:)           => null()   !shoot-root nonstructural element transfer, [ g d-2 h-1]
  contains
    procedure, public :: Init  => plt_bgcrate_init
    procedure, public :: Destroy  => plt_bgcrate_destroy
  end type plant_bgcrate_type


  type(plant_bgcrate_type)  , public, target :: plt_bgcr      !bgc reaction

contains

  subroutine plt_bgcrate_init(this)
  implicit none
  class(plant_bgcrate_type) :: this

  allocate(this%RootCO2Emis2Root_vr(JZ1)); this%RootCO2Emis2Root_vr=spval
  allocate(this%RootCO2Emis2Root_pvr(JZ1,JP1)); this%RootCO2Emis2Root_pvr=spval
  allocate(this%RUptkRootO2_vr(JZ1)); this%RUptkRootO2_vr=0._r8
  allocate(this%RH2PO4EcoDmndSoilPrev_vr(0:JZ1)); this%RH2PO4EcoDmndSoilPrev_vr=spval
  allocate(this%RH2PO4EcoDmndBandPrev_vr(0:JZ1)); this%RH2PO4EcoDmndBandPrev_vr=spval
  allocate(this%RH1PO4EcoDmndSoilPrev_vr(0:JZ1)); this%RH1PO4EcoDmndSoilPrev_vr=spval
  allocate(this%RH1PO4EcoDmndBandPrev_vr(0:JZ1)); this%RH1PO4EcoDmndBandPrev_vr=spval
  allocate(this%RNO3EcoDmndSoilPrev_vr(0:JZ1)); this%RNO3EcoDmndSoilPrev_vr=spval
  allocate(this%RNH4EcoDmndSoilPrev_vr(0:JZ1)); this%RNH4EcoDmndSoilPrev_vr =spval
  allocate(this%RNH4EcoDmndBandPrev_vr(0:JZ1)); this%RNH4EcoDmndBandPrev_vr=spval
  allocate(this%RNO3EcoDmndBandPrev_vr(0:JZ1)); this%RNO3EcoDmndBandPrev_vr=spval
  allocate(this%RGasTranspFlxPrev_vr(idg_beg:idg_end,0:JZ1)); this%RGasTranspFlxPrev_vr=spval
  allocate(this%RO2AquaSourcePrev_vr(0:JZ1)); this%RO2AquaSourcePrev_vr=spval
  allocate(this%RO2EcoDmndPrev_vr(0:JZ1)); this%RO2EcoDmndPrev_vr=spval
  allocate(this%LitrfalStrutElms_vr(NumPlantChemElms,jsken,NumOfPlantLitrCmplxs,0:JZ1));this%LitrfalStrutElms_vr=spval
  allocate(this%GrossCO2Fix_pft(JP1));this%GrossCO2Fix_pft=spval
  allocate(this%REcoDOMProd_vr(1:NumPlantChemElms,1:jcplx,0:JZ1));this%REcoDOMProd_vr=spval
  allocate(this%CO2NetFix_pft(JP1));this%CO2NetFix_pft=spval
  allocate(this%RCanMaintDef_CO2_pft(JP1));this%RCanMaintDef_CO2_pft=spval
  allocate(this%RootGasLossDisturb_pft(idg_beg:idg_NH3,JP1));this%RootGasLossDisturb_pft=spval
  allocate(this%GrossResp_pft(JP1));this%GrossResp_pft=spval
  allocate(this%CanopyGrosRCO2_pft(JP1));this%CanopyGrosRCO2_pft=spval
  allocate(this%CanopyResp_brch(MaxNumBranches,JP1));this%CanopyResp_brch=spval
  allocate(this%RootAutoCO2_pft(JP1));this%RootAutoCO2_pft=spval
  allocate(this%PlantN2Fix_CumYr_pft(JP1));this%PlantN2Fix_CumYr_pft=spval
  allocate(this%NH3Emis_CumYr_pft(JP1));this%NH3Emis_CumYr_pft=spval
  allocate(this%NodulInfectElms_pft(NumPlantChemElms,JP1));this%NodulInfectElms_pft=spval
  allocate(this%SurfLitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP1));this%SurfLitrfalStrutElms_CumYr_pft=0._r8
  allocate(this%LitrFallStrutElms_col(NumPlantChemElms));this%LitrFallStrutElms_col=0._r8
  allocate(this%NetPrimProduct_pft(JP1));this%NetPrimProduct_pft=spval
  allocate(this%PTSHTR_pft(JP1)); this%PTSHTR_pft(:)=spval
  allocate(this%RootShootExch_pvr(NumPlantChemElms,JZ1,JP1)); this%RootShootExch_pvr=0._r8
  allocate(this%Nutruptk_fClim_rpvr(jroots,JZ1,JP1));this%Nutruptk_fClim_rpvr=0._r8
  allocate(this%Nutruptk_fNlim_rpvr(jroots,JZ1,JP1));this%Nutruptk_fNlim_rpvr=0._r8
  allocate(this%Nutruptk_fPlim_rpvr(jroots,JZ1,JP1));this%Nutruptk_fPlim_rpvr=0._r8
  allocate(this%Nutruptk_fProtC_rpvr(jroots,JZ1,JP1));this%Nutruptk_fProtC_rpvr=0._r8
  allocate(this%NH3Dep2Can_pft(JP1));this%NH3Dep2Can_pft=0._r8
  allocate(this%tRootMycoExud2Soil_vr(NumPlantChemElms,1:jcplx,JZ1));this%tRootMycoExud2Soil_vr=spval
  allocate(this%RootN2Fix_pvr(JZ1,JP1));this%RootN2Fix_pvr=0._r8
  allocate(this%RCO2Nodule_pvr(JZ1,JP1)); this%RCO2Nodule_pvr=0._r8
  allocate(this%RNFixCO2_pft(JP1)); this%RNFixCO2_pft=0._r8
  allocate(this%RootO2_TotSink_pvr(jroots,JZ1,JP1)); this%RootO2_TotSink_pvr=0._r8
  allocate(this%RootO2_TotSink_vr(JZ1)); this%RootO2_TotSink_vr=0._r8
  allocate(this%CanopyRespC_CumYr_pft(JP1));this%CanopyRespC_CumYr_pft=spval
  allocate(this%CanopyNLimFactor_brch(MaxNumBranches,JP1));this%CanopyNLimFactor_brch=1._r8
  allocate(this%CanopyPLimFactor_brch(MaxNumBranches,JP1));this%CanopyPLimFactor_brch=1._r8
  allocate(this%REcoH1PO4DmndBand_vr(0:JZ1));this%REcoH1PO4DmndBand_vr=spval
  allocate(this%REcoNO3DmndSoil_vr(0:JZ1));this%REcoNO3DmndSoil_vr=spval
  allocate(this%REcoH2PO4DmndSoil_vr(0:JZ1));this%REcoH2PO4DmndSoil_vr=spval
  allocate(this%REcoNH4DmndSoil_vr(0:JZ1));this%REcoNH4DmndSoil_vr=spval
  allocate(this%REcoO2DmndResp_vr(0:JZ1));this%REcoO2DmndResp_vr=spval
  allocate(this%REcoH2PO4DmndBand_vr(0:JZ1));this%REcoH2PO4DmndBand_vr=spval
  allocate(this%REcoNO3DmndBand_vr(0:JZ1));this%REcoNO3DmndBand_vr=spval
  allocate(this%REcoNH4DmndBand_vr(0:JZ1));this%REcoNH4DmndBand_vr=spval
  allocate(this%REcoH1PO4DmndSoil_vr(0:JZ1));this%REcoH1PO4DmndSoil_vr=spval
  allocate(this%LitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP1));this%LitrfalStrutElms_CumYr_pft=0._r8
  allocate(this%LitrfallElms_pft(NumPlantChemElms,JP1));this%LitrfallElms_pft=spval
  allocate(this%SSXfer2ShootElms_pft(NumPlantChemElms,JP1));this%SSXfer2ShootElms_pft=spval
  allocate(this%SSXferElms_pft(NumPlantChemElms,JP1));this%SSXferElms_pft=spval
  allocate(this%Xfer2RootsC_pft(JP1));this%Xfer2RootsC_pft=0._r8
  allocate(this%LitrfallElms_pvr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ1,JP1));this%LitrfallElms_pvr=spval
  allocate(this%LitrFallElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LitrFallElms_brch=spval
  allocate(this%LitrfallAbvgElms_pft(NumPlantChemElms,JP1)); this%LitrfallAbvgElms_pft=spval
  allocate(this%SurfLitrfallElms_pft(NumPlantChemElms,JP1)); this%SurfLitrfallElms_pft=spval
  allocate(this%LitrfallBlgrElms_pft(NumPlantChemElms,JP1)); this%LitrfallBlgrElms_pft=spval
  allocate(this%RootMaintDef_CO2_pvr(jroots,JZ1,JP1));this%RootMaintDef_CO2_pvr=0._r8
  allocate(this%ShootRootXferElm_pft(NumPlantChemElms,JP1)); this%ShootRootXferElm_pft=spval
  end subroutine plt_bgcrate_init

  subroutine plt_bgcrate_destroy(this)

  implicit none
  class(plant_bgcrate_type) :: this


  end subroutine plt_bgcrate_destroy
end module PlantBGCRatesAPIData
