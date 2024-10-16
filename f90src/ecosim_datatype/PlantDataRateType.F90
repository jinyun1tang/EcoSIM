module PlantDataRateType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  Eco_NEE_col(:,:)                         !total canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  NH3Dep2Can_pft(:,:,:)                       !canopy NH3 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  NodulInfectElms_pft(:,:,:,:)             
  real(r8),target,allocatable ::  NodulInfectElmsCum_pft(:,:,:,:)             
  real(r8),target,allocatable ::  NH3Emis_CumYr_pft(:,:,:)                       !total canopy NH3 flux, [g d-2 ]
  real(r8),target,allocatable ::  SurfLitrfalStrutElms_CumYr_pft(:,:,:,:)                     !total surface LitrFall element, [g d-2]
  real(r8),target,allocatable ::  RootMycoExudElm_pvr(:,:,:,:,:,:,:)              !root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNutUptake_pvr(:,:,:,:,:,:)       !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_pvr(:,:,:,:)             !root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPGasSol_vr(:,:,:,:,:,:)          !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH2PO4DmndSoil_pvr(:,:,:,:,:)                  !root uptake of H2PO4 non-band
  real(r8),target,allocatable ::  RootH2PO4DmndBand_pvr(:,:,:,:,:)                  !root uptake of H2PO4 band
  real(r8),target,allocatable ::  RootH1PO4DmndSoil_pvr(:,:,:,:,:)                  !HPO4 demand in non-band by each root population
  real(r8),target,allocatable ::  RootH1PO4DmndBand_pvr(:,:,:,:,:)                  !HPO4 demand in band by each root population
  real(r8),target,allocatable ::  LeafElmntRemobFlx_brch(:,:,:,:,:)                   !element translocated from leaf during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCZSX(:,:,:,:)                     !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCPSX(:,:,:,:)                     !P translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  PetioleChemElmRemobFlx_brch(:,:,:,:,:)                   !element translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  GrossCO2Fix_pft(:,:,:)                       !total gross CO2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  GrossCO2Fix_CumYr_pft(:,:,:)
  real(r8),target,allocatable ::  LitrfalStrutElms_pft(:,:,:,:)                     !total plant element LitrFall , [g d-2 ]
  real(r8),target,allocatable ::  PlantN2Fix_CumYr_pft(:,:,:)                      !total plant N2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  GrossRespC_CumYr_pft(:,:,:)                       !total plant respiration, [g d-2 ]
  real(r8),target,allocatable ::  GrossResp_pft(:,:,:)
  real(r8),target,allocatable ::  ElmBalanceCum_pft(:,:,:,:)                      !plant element balance, [g d-2]
  real(r8),target,allocatable ::  LitrfalStrutElms_CumYr_pft(:,:,:,:)                     !plant element LitrFall, [g d-2 h-1]
  real(r8),target,allocatable ::  LitrfalStrutElms_pvr(:,:,:,:,:,:,:)                !plant LitrFall element, [g d-2 h-1]
  real(r8),target,allocatable ::  NetPrimProduct_pft(:,:,:)                        !total net primary productivity, [g d-2]
  real(r8),target,allocatable ::  ETCanopy_CumYr_pft(:,:,:)                       !total transpiration, [m d-2], <0 into atmosphere
  real(r8),target,allocatable ::  CanopyRespC_CumYr_pft(:,:,:)                       !total autotrophic respiration, [g d-2 ]
  real(r8),target,allocatable ::  EcoHavstElmnt_CumYr_pft(:,:,:,:)                     !plant element harvest, [g d-2 ]
  real(r8),target,allocatable ::  EcoHavstElmntCum_pft(:,:,:,:)                    !total plant harvest, [g d-2 ]
  real(r8),target,allocatable ::  CO2ByFire_CumYr_pft(:,:,:)                       !plant CO2 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  CH4ByFire_CumYr_pft(:,:,:)                       !plant CH4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  O2ByFire_CumYr_pft(:,:,:)                       !plant O2 uptake from fire, [g d-2 ]
  real(r8),target,allocatable ::  NH3byFire_CumYr_pft(:,:,:)                       !plant NH3 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  N2ObyFire_CumYr_pft(:,:,:)                       !plant N2O emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  PO4byFire_CumYr_pft(:,:,:)                       !plant PO4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  RootO2Dmnd4Resp_pvr(:,:,:,:,:)                   !root  O2 demand from respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_air2root_flx__pvr(:,:,:,:,:,:)                  !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_Root_DisEvap_flx_vr(:,:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Emis_pvr(:,:,:,:,:)                   !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RootO2Uptk_pvr(:,:,:,:,:)                  !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RootRespPotent_pvr(:,:,:,:,:)                   !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Autor_pvr(:,:,:,:,:)                   !root respiration constrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootMycoExudElms_pft(:,:,:,:)                     !total root uptake (+ve) - exudation (-ve) of dissovled element, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNH4Uptake_pft(:,:,:)                       !total root uptake of NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3Uptake_pft(:,:,:)                       !total root uptake of NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH2PO4Uptake_pft(:,:,:)                       !total root uptake of PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootHPO4Uptake_pft(:,:,:)                       !total root uptake of HPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_pft(:,:,:)                        !total root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RootGasLossDisturb_pft(:,:,:,:)                !gas flux from root disturbance [g d-2 h-1]
  real(r8),target,allocatable ::  RootOUlmNutUptake_pvr(:,:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCUlmNutUptake_pvr(:,:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2EmisPot_pvr(:,:,:,:,:)                   !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNH4DmndSoil_pvr(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3DmndSoil_pvr(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNH4DmndBand_pvr(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3DmndBand_pvr(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH3Z(:,:,:)                       !gaseous NH3 flux fron root disturbance non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  NH3Dep2Can_brch(:,:,:,:)                     !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8),target,allocatable ::  RAutoRootO2Limter_pvr(:,:,:,:,:)                     !O2 constraint to root respiration, []
  real(r8),target,allocatable ::  PlantRootSoilElmNetX_pft(:,:,:,:)                    !net root element uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  PlantExudElm_CumYr_pft(:,:,:,:)                    !total net root element uptake (+ve) - exudation (-ve), [g d-2 ]
  real(r8),target,allocatable ::  RootUptk_N_CumYr_pft(:,:,:)
  real(r8),target,allocatable ::  RootUptk_P_CumYr_pft(:,:,:)
  real(r8),target,allocatable ::  TPlantRootH2OUptake_vr(:,:,:)                      !total root water uptake, [m3 d-2]
  real(r8),target,allocatable ::  THeatRootUptake_vr(:,:,:)                       !vertically profile of root heat uptake, [MJ d-2]
  real(r8),target,allocatable :: THeatRootUptake_col(:,:)                        !total root heat uptake, [MJ d-2]
  real(r8),target,allocatable ::  trcg_air2root_flx_vr(:,:,:,:)                 !total internal root gas flux , [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_root_vr(:,:,:,:)                  !total root internal gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_plant_uptake_vr(:,:,:,:)      !total root-soil solute flux, [g d-2 h-1]
  real(r8),target,allocatable ::  tRootMycoExud2Soil_vr(:,:,:,:,:)                  !total root element exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  tRootCO2Emis_vr(:,:,:)                       !total root CO2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  tRO2MicrbUptk_vr(:,:,:)                      !total root internal O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  totRootLenDens_vr(:,:,:)                       !total root length density, [m m-3]
  real(r8),target,allocatable ::  REcoO2DmndResp_vr(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RO2EcoDmndPrev_vr(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNH4DmndSoil_vr(:,:,:)                       !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH4EcoDmndSoilPrev_vr(:,:,:)                       !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNO3DmndSoil_vr(:,:,:)                       !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO3EcoDmndSoilPrev_vr(:,:,:)                       !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkSoil_vr(:,:,:)                       !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkSoilPrev_vr(:,:,:)                       !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoH2PO4DmndSoil_vr(:,:,:)                       !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2PO4EcoDmndSoilPrev_vr(:,:,:)                       !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2OEcoUptkSoil_vr(:,:,:)                       !total root + microbial N2O uptake , [g d-2 h-1]
  real(r8),target,allocatable ::  RN2OEcoUptkSoilPrev_vr(:,:,:)                       !total root + microbial N2O uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNH4DmndBand_vr(:,:,:)                       !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH4EcoDmndBandPrev_vr(:,:,:)                       !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNO3DmndBand_vr(:,:,:)                       !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO3EcoDmndBandPrev_vr(:,:,:)                       !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkBand_vr(:,:,:)                       !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkBandPrev_vr(:,:,:)                       !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoH2PO4DmndBand_vr(:,:,:)                       !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2PO4EcoDmndBandPrev_vr(:,:,:)                       !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RDOMEcoDmndK_vr(:,:,:,:)                     !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RDOMEcoDmndPrev_vr(:,:,:,:)                     !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RAcetateEcoDmndK_vr(:,:,:,:)                     !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RAcetateEcoDmndPrev_vr(:,:,:,:)                     !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  TH2GZ(:,:)                         !total root H2 flux, [g d-2]
  private :: InitAllocate
  contains

!----------------------------------------------------------------------
  subroutine InitPlantRates(NumOfPlantLitrCmplxs,jroots)
  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs
  integer, intent(in) :: jroots
  call InitAllocate(NumOfPlantLitrCmplxs,jroots)

  end subroutine InitPlantRates
!----------------------------------------------------------------------

  subroutine InitAllocate(NumOfPlantLitrCmplxs,jroots)

  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs
  integer, intent(in) :: jroots    !number of root types, root,mycos
  allocate(Eco_NEE_col(JY,JX));       Eco_NEE_col=0._r8
  allocate(NH3Dep2Can_pft(JP,JY,JX));    NH3Dep2Can_pft=0._r8
  allocate(NH3Emis_CumYr_pft(JP,JY,JX));    NH3Emis_CumYr_pft=0._r8  
  allocate(NodulInfectElms_pft(NumPlantChemElms,JP,JY,JX));NodulInfectElms_pft=0._r8
  allocate(NodulInfectElmsCum_pft(NumPlantChemElms,JP,JY,JX));NodulInfectElmsCum_pft=0._r8  
  allocate(SurfLitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP,JY,JX));    SurfLitrfalStrutElms_CumYr_pft=0._r8
  allocate(RootMycoExudElm_pvr(NumPlantChemElms,jroots,1:jcplx,JZ,JP,JY,JX));RootMycoExudElm_pvr=0._r8
  allocate(RootNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootNutUptake_pvr=0._r8
  allocate(RootN2Fix_pvr(JZ,JP,JY,JX)); RootN2Fix_pvr=0._r8
  allocate(RootH2PO4DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootH2PO4DmndSoil_pvr=0._r8
  allocate(RootH2PO4DmndBand_pvr(jroots,JZ,JP,JY,JX));RootH2PO4DmndBand_pvr=0._r8
  allocate(RootH1PO4DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootH1PO4DmndSoil_pvr=0._r8
  allocate(RootH1PO4DmndBand_pvr(jroots,JZ,JP,JY,JX));RootH1PO4DmndBand_pvr=0._r8
  allocate(LeafElmntRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); LeafElmntRemobFlx_brch=0._r8
  allocate(PetioleChemElmRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); PetioleChemElmRemobFlx_brch=0._r8
  allocate(GrossCO2Fix_pft(JP,JY,JX));    GrossCO2Fix_pft=0._r8
  allocate(GrossCO2Fix_CumYr_pft(JP,JY,JX)); GrossCO2Fix_CumYr_pft=0._r8
  allocate(LitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP,JY,JX));    LitrfalStrutElms_CumYr_pft=0._r8
  allocate(PlantN2Fix_CumYr_pft(JP,JY,JX));   PlantN2Fix_CumYr_pft=0._r8
  allocate(GrossResp_pft(JP,JY,JX));    GrossResp_pft=0._r8
  allocate(GrossRespC_CumYr_pft(JP,JY,JX)); GrossRespC_CumYr_pft=0._r8
  allocate(ElmBalanceCum_pft(NumPlantChemElms,JP,JY,JX));     ElmBalanceCum_pft=0._r8
  allocate(LitrfalStrutElms_pft(NumPlantChemElms,JP,JY,JX));    LitrfalStrutElms_pft=0._r8
  allocate(LitrfalStrutElms_pvr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JP,JY,JX));LitrfalStrutElms_pvr=0._r8
  allocate(NetPrimProduct_pft(JP,JY,JX));     NetPrimProduct_pft=0._r8
  allocate(ETCanopy_CumYr_pft(JP,JY,JX));    ETCanopy_CumYr_pft=0._r8
  allocate(CanopyRespC_CumYr_pft(JP,JY,JX));    CanopyRespC_CumYr_pft=0._r8
  allocate(EcoHavstElmnt_CumYr_pft(NumPlantChemElms,JP,JY,JX));    EcoHavstElmnt_CumYr_pft=0._r8
  allocate(EcoHavstElmntCum_pft(NumPlantChemElms,JP,JY,JX));   EcoHavstElmntCum_pft=0._r8
  allocate(CO2ByFire_CumYr_pft(JP,JY,JX));    CO2ByFire_CumYr_pft=0._r8
  allocate(CH4ByFire_CumYr_pft(JP,JY,JX));    CH4ByFire_CumYr_pft=0._r8
  allocate(O2ByFire_CumYr_pft(JP,JY,JX));    O2ByFire_CumYr_pft=0._r8
  allocate(NH3byFire_CumYr_pft(JP,JY,JX));    NH3byFire_CumYr_pft=0._r8
  allocate(N2ObyFire_CumYr_pft(JP,JY,JX));    N2ObyFire_CumYr_pft=0._r8
  allocate(PO4byFire_CumYr_pft(JP,JY,JX));    PO4byFire_CumYr_pft=0._r8
  allocate(RootO2Dmnd4Resp_pvr(jroots,JZ,JP,JY,JX));RootO2Dmnd4Resp_pvr=0._r8
  allocate(trcg_air2root_flx__pvr(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_air2root_flx__pvr=0._r8
  allocate(trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_Root_DisEvap_flx_vr=0._r8
  allocate(RUPGasSol_vr(idg_beg:idg_end,jroots,JZ,JP,JY,JX));RUPGasSol_vr=0._r8
  allocate(RootCO2Emis_pvr(jroots,JZ,JP,JY,JX));RootCO2Emis_pvr=0._r8
  allocate(RootO2Uptk_pvr(jroots,JZ,JP,JY,JX));RootO2Uptk_pvr=0._r8
  allocate(RootRespPotent_pvr(jroots,JZ,JP,JY,JX));RootRespPotent_pvr=0._r8
  allocate(RootCO2Autor_pvr(jroots,JZ,JP,JY,JX));RootCO2Autor_pvr=0._r8
  allocate(RootMycoExudElms_pft(1:NumPlantChemElms,JP,JY,JX));    RootMycoExudElms_pft=0._r8
  allocate(RootNH4Uptake_pft(JP,JY,JX));    RootNH4Uptake_pft=0._r8
  allocate(RootNO3Uptake_pft(JP,JY,JX));    RootNO3Uptake_pft=0._r8
  allocate(RootH2PO4Uptake_pft(JP,JY,JX));    RootH2PO4Uptake_pft=0._r8
  allocate(RootHPO4Uptake_pft(JP,JY,JX));    RootHPO4Uptake_pft=0._r8
  allocate(RootN2Fix_pft(JP,JY,JX));     RootN2Fix_pft=0._r8
  allocate(RootGasLossDisturb_pft(idg_beg:idg_end-1,JP,JY,JX)); RootGasLossDisturb_pft=0._r8
  allocate(RootCUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootCUlmNutUptake_pvr=0._r8
  allocate(RootOUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootOUlmNutUptake_pvr=0._r8
  allocate(RootCO2EmisPot_pvr(jroots,JZ,JP,JY,JX));RootCO2EmisPot_pvr=0._r8
  allocate(RootNH4DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootNH4DmndSoil_pvr=0._r8
  allocate(RootNO3DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootNO3DmndSoil_pvr=0._r8
  allocate(RootNH4DmndBand_pvr(jroots,JZ,JP,JY,JX));RootNH4DmndBand_pvr=0._r8
  allocate(RootNO3DmndBand_pvr(jroots,JZ,JP,JY,JX));RootNO3DmndBand_pvr=0._r8
  allocate(RNH3Z(JP,JY,JX));    RNH3Z=0._r8
  allocate(NH3Dep2Can_brch(MaxNumBranches,JP,JY,JX)); NH3Dep2Can_brch=0._r8
  allocate(RAutoRootO2Limter_pvr(jroots,JZ,JP,JY,JX)); RAutoRootO2Limter_pvr=0._r8
  allocate(PlantRootSoilElmNetX_pft(NumPlantChemElms,JP,JY,JX));   PlantRootSoilElmNetX_pft=0._r8
  allocate(PlantExudElm_CumYr_pft(NumPlantChemElms,JP,JY,JX));   PlantExudElm_CumYr_pft=0._r8
  allocate(RootUptk_N_CumYr_pft(JP,JY,JX)); RootUptk_N_CumYr_pft=0._r8
  allocate(RootUptk_P_CumYr_pft(JP,JY,JX)); RootUptk_P_CumYr_pft=0._r8
  allocate(TPlantRootH2OUptake_vr(0:JZ,JY,JX)); TPlantRootH2OUptake_vr=0._r8
  allocate(THeatRootUptake_vr(0:JZ,JY,JX));  THeatRootUptake_vr=0._r8
  allocate(THeatRootUptake_col(JY,JX)); THeatRootUptake_col=0._r8
  allocate(trcg_air2root_flx_vr(idg_beg:idg_end-1,JZ,JY,JX));   trcg_air2root_flx_vr=0._r8
  allocate(trcg_root_vr(idg_beg:idg_end-1,JZ,JY,JX));   trcg_root_vr=0._r8
  allocate(trcs_plant_uptake_vr(ids_beg:ids_end,JZ,JY,JX));    trcs_plant_uptake_vr=0._r8
  allocate(tRootMycoExud2Soil_vr(NumPlantChemElms,1:jcplx,JZ,JY,JX));tRootMycoExud2Soil_vr=0._r8
  allocate(tRootCO2Emis_vr(JZ,JY,JX));    tRootCO2Emis_vr=0._r8
  allocate(tRO2MicrbUptk_vr(JZ,JY,JX));   tRO2MicrbUptk_vr=0._r8
  allocate(totRootLenDens_vr(JZ,JY,JX));    totRootLenDens_vr=0._r8
  allocate(REcoO2DmndResp_vr(0:JZ,JY,JX));  REcoO2DmndResp_vr=0._r8
  allocate(RO2EcoDmndPrev_vr(0:JZ,JY,JX));  RO2EcoDmndPrev_vr=0._r8
  allocate(REcoNH4DmndSoil_vr(0:JZ,JY,JX));  REcoNH4DmndSoil_vr=0._r8
  allocate(RNH4EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RNH4EcoDmndSoilPrev_vr=0._r8
  allocate(REcoNO3DmndSoil_vr(0:JZ,JY,JX));  REcoNO3DmndSoil_vr=0._r8
  allocate(RNO3EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RNO3EcoDmndSoilPrev_vr=0._r8
  allocate(RNO2EcoUptkSoil_vr(0:JZ,JY,JX));  RNO2EcoUptkSoil_vr=0._r8
  allocate(RNO2EcoUptkSoilPrev_vr(0:JZ,JY,JX));  RNO2EcoUptkSoilPrev_vr=0._r8
  allocate(REcoH2PO4DmndSoil_vr(0:JZ,JY,JX));  REcoH2PO4DmndSoil_vr=0._r8
  allocate(RH2PO4EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RH2PO4EcoDmndSoilPrev_vr=0._r8
  allocate(RN2OEcoUptkSoil_vr(0:JZ,JY,JX));  RN2OEcoUptkSoil_vr=0._r8
  allocate(RN2OEcoUptkSoilPrev_vr(0:JZ,JY,JX));  RN2OEcoUptkSoilPrev_vr=0._r8
  allocate(REcoNH4DmndBand_vr(0:JZ,JY,JX));  REcoNH4DmndBand_vr=0._r8
  allocate(RNH4EcoDmndBandPrev_vr(0:JZ,JY,JX));  RNH4EcoDmndBandPrev_vr=0._r8
  allocate(REcoNO3DmndBand_vr(0:JZ,JY,JX));  REcoNO3DmndBand_vr=0._r8
  allocate(RNO3EcoDmndBandPrev_vr(0:JZ,JY,JX));  RNO3EcoDmndBandPrev_vr=0._r8
  allocate(RNO2EcoUptkBand_vr(0:JZ,JY,JX));  RNO2EcoUptkBand_vr=0._r8
  allocate(RNO2EcoUptkBandPrev_vr(0:JZ,JY,JX));  RNO2EcoUptkBandPrev_vr=0._r8
  allocate(REcoH2PO4DmndBand_vr(0:JZ,JY,JX));  REcoH2PO4DmndBand_vr=0._r8
  allocate(RH2PO4EcoDmndBandPrev_vr(0:JZ,JY,JX));  RH2PO4EcoDmndBandPrev_vr=0._r8
  allocate(RDOMEcoDmndK_vr(1:jcplx,0:JZ,JY,JX));RDOMEcoDmndK_vr=0._r8
  allocate(RDOMEcoDmndPrev_vr(1:jcplx,0:JZ,JY,JX));RDOMEcoDmndPrev_vr=0._r8
  allocate(RAcetateEcoDmndK_vr(1:jcplx,0:JZ,JY,JX));RAcetateEcoDmndK_vr=0._r8
  allocate(RAcetateEcoDmndPrev_vr(1:jcplx,0:JZ,JY,JX));RAcetateEcoDmndPrev_vr=0._r8
  allocate(TH2GZ(JY,JX));       TH2GZ=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantRates
  use abortutils, only : destroy
  implicit none
  call destroy(Eco_NEE_col)
  call destroy(NH3Dep2Can_pft)
  call destroy(NH3Emis_CumYr_pft)
  call destroy(SurfLitrfalStrutElms_CumYr_pft)
  call destroy(RootMycoExudElm_pvr)
  call destroy(RootNutUptake_pvr)
  call destroy(RUPGasSol_vr)
  call destroy(RootN2Fix_pvr)
  call destroy(RootH2PO4DmndSoil_pvr)
  call destroy(RootH2PO4DmndBand_pvr)
  call destroy(RootH1PO4DmndSoil_pvr)
  call destroy(RootH1PO4DmndBand_pvr)
  call destroy(LeafElmntRemobFlx_brch)
  call destroy(PetioleChemElmRemobFlx_brch)
  call destroy(GrossCO2Fix_pft)
  call destroy(GrossCO2Fix_CumYr_pft)
  call destroy(LitrfalStrutElms_CumYr_pft)
  call destroy(PlantN2Fix_CumYr_pft)
  call destroy(GrossResp_pft)
  call destroy(GrossRespC_CumYr_pft)
  call destroy(ElmBalanceCum_pft)
  call destroy(LitrfalStrutElms_pft)
  call destroy(LitrfalStrutElms_pvr)
  call destroy(NetPrimProduct_pft)
  call destroy(ETCanopy_CumYr_pft)
  call destroy(CanopyRespC_CumYr_pft)
  call destroy(EcoHavstElmnt_CumYr_pft)
  call destroy(EcoHavstElmntCum_pft)
  call destroy(NodulInfectElms_pft)
  call destroy(NodulInfectElmsCum_pft)
  call destroy(CO2ByFire_CumYr_pft)
  call destroy(CH4ByFire_CumYr_pft)
  call destroy(O2ByFire_CumYr_pft)
  call destroy(NH3byFire_CumYr_pft)
  call destroy(N2ObyFire_CumYr_pft)
  call destroy(PO4byFire_CumYr_pft)
  call destroy(RootO2Dmnd4Resp_pvr)  
  call destroy(RootCO2Emis_pvr)
  call destroy(RootO2Uptk_pvr)
  call destroy(RootRespPotent_pvr)
  call destroy(RootCO2Autor_pvr)
  call destroy(RootMycoExudElms_pft)
  call destroy(RootNH4Uptake_pft)
  call destroy(RootNO3Uptake_pft)
  call destroy(RootH2PO4Uptake_pft)
  call destroy(RootHPO4Uptake_pft)
  call destroy(RootN2Fix_pft)
  call destroy(RootOUlmNutUptake_pvr)
  call destroy(RootCUlmNutUptake_pvr)
  call destroy(RootCO2EmisPot_pvr)
  call destroy(RootNH4DmndSoil_pvr)
  call destroy(RootNO3DmndSoil_pvr)
  call destroy(RootNH4DmndBand_pvr)
  call destroy(RootNO3DmndBand_pvr)
  call destroy(RNH3Z)
  call destroy(NH3Dep2Can_brch)
  call destroy(RAutoRootO2Limter_pvr)
  call destroy(PlantRootSoilElmNetX_pft)
  call destroy(PlantExudElm_CumYr_pft)
  call destroy(RootUptk_N_CumYr_pft)
  call destroy(RootUptk_P_CumYr_pft)
  call destroy(TPlantRootH2OUptake_vr)
  call destroy(THeatRootUptake_vr)
  call destroy(THeatRootUptake_col)
  call destroy(tRootMycoExud2Soil_vr)
  call destroy(tRootCO2Emis_vr)
  call destroy(tRO2MicrbUptk_vr)
  call destroy(totRootLenDens_vr)
  call destroy(REcoO2DmndResp_vr)
  call destroy(RO2EcoDmndPrev_vr)
  call destroy(REcoNH4DmndSoil_vr)
  call destroy(RNH4EcoDmndSoilPrev_vr)
  call destroy(REcoNO3DmndSoil_vr)
  call destroy(RNO3EcoDmndSoilPrev_vr)
  call destroy(RNO2EcoUptkSoil_vr)
  call destroy(RNO2EcoUptkSoilPrev_vr)
  call destroy(REcoH2PO4DmndSoil_vr)
  call destroy(RH2PO4EcoDmndSoilPrev_vr)
  call destroy(RN2OEcoUptkSoil_vr)
  call destroy(RN2OEcoUptkSoilPrev_vr)
  call destroy(REcoNH4DmndBand_vr)
  call destroy(RNH4EcoDmndBandPrev_vr)
  call destroy(REcoNO3DmndBand_vr)
  call destroy(RNO3EcoDmndBandPrev_vr)
  call destroy(RNO2EcoUptkBand_vr)
  call destroy(RNO2EcoUptkBandPrev_vr)
  call destroy(REcoH2PO4DmndBand_vr)
  call destroy(RH2PO4EcoDmndBandPrev_vr)
  call destroy(RDOMEcoDmndK_vr)
  call destroy(RDOMEcoDmndPrev_vr)
  call destroy(RAcetateEcoDmndK_vr)
  call destroy(RAcetateEcoDmndPrev_vr)
  call destroy(TH2GZ)
  call destroy(trcs_plant_uptake_vr)
  end subroutine DestructPlantRates

end module PlantDataRateType
