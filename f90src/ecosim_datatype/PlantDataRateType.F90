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
  real(r8),target,allocatable ::  NH3EmiCum_pft(:,:,:)                       !total canopy NH3 flux, [g d-2 ]
  real(r8),target,allocatable ::  SurfLitrfalStrutElms_pft(:,:,:,:)                     !total surface LitrFall element, [g d-2]
  real(r8),target,allocatable ::  RootMycoExudElm_pvr(:,:,:,:,:,:,:)              !root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNutUptake_pvr(:,:,:,:,:,:)       !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_pvr(:,:,:,:)             !root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPGasSol_vr(:,:,:,:,:,:)          !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPP2P(:,:,:,:,:)                  !root uptake of H2PO4 non-band
  real(r8),target,allocatable ::  RUPP2B(:,:,:,:,:)                  !root uptake of H2PO4 band
  real(r8),target,allocatable ::  RUPP1P(:,:,:,:,:)                  !HPO4 demand in non-band by each root population
  real(r8),target,allocatable ::  RUPP1B(:,:,:,:,:)                  !HPO4 demand in band by each root population
  real(r8),target,allocatable ::  LeafElmntRemobFlx_brch(:,:,:,:,:)                   !element translocated from leaf during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCZSX(:,:,:,:)                     !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCPSX(:,:,:,:)                     !P translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  PetioleChemElmRemobFlx_brch(:,:,:,:,:)                   !element translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  GrossCO2Fix_pft(:,:,:)                       !total gross CO2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  LitrfalStrutElms_pft(:,:,:,:)                     !total plant element LitrFall , [g d-2 ]
  real(r8),target,allocatable ::  PlantN2FixCum_pft(:,:,:)                      !total plant N2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  GrossResp_pft(:,:,:)                       !total plant respiration, [g d-2 ]
  real(r8),target,allocatable ::  ElmBalanceCum_pft(:,:,:,:)                      !plant element balance, [g d-2]
  real(r8),target,allocatable ::  LitrfalStrutElms_pft(:,:,:,:)                     !plant element LitrFall, [g d-2 h-1]
  real(r8),target,allocatable ::  LitfalStrutElms_pvr(:,:,:,:,:,:,:)                !plant LitrFall element, [g d-2 h-1]
  real(r8),target,allocatable ::  NetPrimProduct_pft(:,:,:)                        !total net primary productivity, [g d-2]
  real(r8),target,allocatable ::  ETCanopy_pft(:,:,:)                       !total transpiration, [m d-2], <0 into atmosphere
  real(r8),target,allocatable ::  CanopyPlusNoduRespC_pft(:,:,:)                       !total autotrophic respiration, [g d-2 ]
  real(r8),target,allocatable ::  EcoHavstElmnt_pft(:,:,:,:)                     !plant element harvest, [g d-2 ]
  real(r8),target,allocatable ::  EcoHavstElmntCum_pft(:,:,:,:)                    !total plant harvest, [g d-2 ]
  real(r8),target,allocatable ::  CO2ByFire_pft(:,:,:)                       !plant CO2 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  CH4ByFire_pft(:,:,:)                       !plant CH4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  O2ByFire_pft(:,:,:)                       !plant O2 uptake from fire, [g d-2 ]
  real(r8),target,allocatable ::  NH3byFire_pft(:,:,:)                       !plant NH3 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  N2ObyFire_pft(:,:,:)                       !plant N2O emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  PO4byFire_pft(:,:,:)                       !plant PO4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  ROXYP(:,:,:,:,:)                   !root  O2 demand from respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_air2root_flx__pvr(:,:,:,:,:,:)                  !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_Root_DisEvap_flx_vr(:,:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2P(:,:,:,:,:)                   !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RUPOXP(:,:,:,:,:)                  !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RootRespPotent_pvr(:,:,:,:,:)                   !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2A_pvr(:,:,:,:,:)                   !root respiration constrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootMycoExudElms_pft(:,:,:,:)                     !total root uptake (+ve) - exudation (-ve) of dissovled element, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNH4Uptake_pft(:,:,:)                       !total root uptake of NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3Uptake_pft(:,:,:)                       !total root uptake of NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH2PO4Uptake_pft(:,:,:)                       !total root uptake of PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootHPO4Uptake_pft(:,:,:)                       !total root uptake of HPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_pft(:,:,:)                        !total root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RootGasLossDisturb_pft(:,:,:,:)                !gas flux from root disturbance [g d-2 h-1]
  real(r8),target,allocatable ::  RootOUlmNutUptake_pvr(:,:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCUlmNutUptake_pvr(:,:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2N_pvr(:,:,:,:,:)                   !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNHP(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNOP(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNBP(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNXP(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH3Z(:,:,:)                       !gaseous NH3 flux fron root disturbance non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  NH3Dep2_brch(:,:,:,:)                     !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8),target,allocatable ::  RAutoRootO2Limter_pvr(:,:,:,:,:)                     !O2 constraint to root respiration, []
  real(r8),target,allocatable ::  RH2GZ(:,:,:)                       !gaseous H2 flux fron root disturbance, [g d-2 h-1]
  real(r8),target,allocatable ::  PlantRootSoilElmNetX_pft(:,:,:,:)                    !net root element uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  PlantExudChemElmCum_pft(:,:,:,:)                    !total net root element uptake (+ve) - exudation (-ve), [g d-2 ]
  real(r8),target,allocatable ::  GridPlantRootH2OUptake_vr(:,:,:)                      !total root water uptake, [m3 d-2]
  real(r8),target,allocatable ::  THeatRootUptake(:,:,:)                       !total root heat uptake, [MJ d-2]
  real(r8),target,allocatable ::  trcg_air2root_flx_vr(:,:,:,:)                 !total internal root gas flux , [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_TLP(:,:,:,:)                  !total root internal gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_plant_uptake_vr(:,:,:,:)      !total root-soil solute flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TDFOME(:,:,:,:,:)                  !total root element exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  TCO2P(:,:,:)                       !total root CO2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPOXP(:,:,:)                      !total root internal O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RTDNT(:,:,:)                       !total root length density, [m m-3]
  real(r8),target,allocatable ::  ROXYX(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  ROXYY(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH4X(:,:,:)                       !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH4Y(:,:,:)                       !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO3X(:,:,:)                       !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO3Y(:,:,:)                       !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2X(:,:,:)                       !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2Y(:,:,:)                       !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RPO4X(:,:,:)                       !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RPO4Y(:,:,:)                       !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2OX(:,:,:)                       !total root + microbial N2O uptake , [g d-2 h-1]
  real(r8),target,allocatable ::  RN2OY(:,:,:)                       !total root + microbial N2O uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RNHBX(:,:,:)                       !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNHBY(:,:,:)                       !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN3BX(:,:,:)                       !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN3BY(:,:,:)                       !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2BX(:,:,:)                       !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2BY(:,:,:)                       !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RPOBX(:,:,:)                       !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RPOBY(:,:,:)                       !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  ROQCX(:,:,:,:)                     !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  ROQCY(:,:,:,:)                     !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  ROQAX(:,:,:,:)                     !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  ROQAY(:,:,:,:)                     !total root + microbial acetate uptake, [g d-2 h-1]
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
  allocate(NH3EmiCum_pft(JP,JY,JX));    NH3EmiCum_pft=0._r8
  allocate(SurfLitrfalStrutElms_pft(NumPlantChemElms,JP,JY,JX));    SurfLitrfalStrutElms_pft=0._r8
  allocate(RootMycoExudElm_pvr(NumPlantChemElms,2,1:jcplx,JZ,JP,JY,JX));RootMycoExudElm_pvr=0._r8
  allocate(RootNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootNutUptake_pvr=0._r8
  allocate(RootN2Fix_pvr(JZ,JP,JY,JX)); RootN2Fix_pvr=0._r8
  allocate(RUPP2P(jroots,JZ,JP,JY,JX));RUPP2P=0._r8
  allocate(RUPP2B(jroots,JZ,JP,JY,JX));RUPP2B=0._r8
  allocate(RUPP1P(jroots,JZ,JP,JY,JX));RUPP1P=0._r8
  allocate(RUPP1B(jroots,JZ,JP,JY,JX));RUPP1B=0._r8
  allocate(LeafElmntRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); LeafElmntRemobFlx_brch=0._r8
  allocate(PetioleChemElmRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); PetioleChemElmRemobFlx_brch=0._r8
  allocate(GrossCO2Fix_pft(JP,JY,JX));    GrossCO2Fix_pft=0._r8
  allocate(LitrfalStrutElms_pft(NumPlantChemElms,JP,JY,JX));    LitrfalStrutElms_pft=0._r8
  allocate(PlantN2FixCum_pft(JP,JY,JX));   PlantN2FixCum_pft=0._r8
  allocate(GrossResp_pft(JP,JY,JX));    GrossResp_pft=0._r8
  allocate(ElmBalanceCum_pft(NumPlantChemElms,JP,JY,JX));     ElmBalanceCum_pft=0._r8
  allocate(LitrfalStrutElms_pft(NumPlantChemElms,JP,JY,JX));    LitrfalStrutElms_pft=0._r8
  allocate(LitfalStrutElms_pvr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JP,JY,JX));LitfalStrutElms_pvr=0._r8
  allocate(NetPrimProduct_pft(JP,JY,JX));     NetPrimProduct_pft=0._r8
  allocate(ETCanopy_pft(JP,JY,JX));    ETCanopy_pft=0._r8
  allocate(CanopyPlusNoduRespC_pft(JP,JY,JX));    CanopyPlusNoduRespC_pft=0._r8
  allocate(EcoHavstElmnt_pft(NumPlantChemElms,JP,JY,JX));    EcoHavstElmnt_pft=0._r8
  allocate(EcoHavstElmntCum_pft(NumPlantChemElms,JP,JY,JX));   EcoHavstElmntCum_pft=0._r8
  allocate(CO2ByFire_pft(JP,JY,JX));    CO2ByFire_pft=0._r8
  allocate(CH4ByFire_pft(JP,JY,JX));    CH4ByFire_pft=0._r8
  allocate(O2ByFire_pft(JP,JY,JX));    O2ByFire_pft=0._r8
  allocate(NH3byFire_pft(JP,JY,JX));    NH3byFire_pft=0._r8
  allocate(N2ObyFire_pft(JP,JY,JX));    N2ObyFire_pft=0._r8
  allocate(PO4byFire_pft(JP,JY,JX));    PO4byFire_pft=0._r8
  allocate(ROXYP(jroots,JZ,JP,JY,JX));ROXYP=0._r8
  allocate(trcg_air2root_flx__pvr(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_air2root_flx__pvr=0._r8
  allocate(trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_Root_DisEvap_flx_vr=0._r8
  allocate(RUPGasSol_vr(idg_beg:idg_end,jroots,JZ,JP,JY,JX));RUPGasSol_vr=0._r8
  allocate(RCO2P(jroots,JZ,JP,JY,JX));RCO2P=0._r8
  allocate(RUPOXP(jroots,JZ,JP,JY,JX));RUPOXP=0._r8
  allocate(RootRespPotent_pvr(jroots,JZ,JP,JY,JX));RootRespPotent_pvr=0._r8
  allocate(RCO2A_pvr(jroots,JZ,JP,JY,JX));RCO2A_pvr=0._r8
  allocate(RootMycoExudElms_pft(1:NumPlantChemElms,JP,JY,JX));    RootMycoExudElms_pft=0._r8
  allocate(RootNH4Uptake_pft(JP,JY,JX));    RootNH4Uptake_pft=0._r8
  allocate(RootNO3Uptake_pft(JP,JY,JX));    RootNO3Uptake_pft=0._r8
  allocate(RootH2PO4Uptake_pft(JP,JY,JX));    RootH2PO4Uptake_pft=0._r8
  allocate(RootHPO4Uptake_pft(JP,JY,JX));    RootHPO4Uptake_pft=0._r8
  allocate(RootN2Fix_pft(JP,JY,JX));     RootN2Fix_pft=0._r8
  allocate(RootGasLossDisturb_pft(idg_beg:idg_end-1,JP,JY,JX)); RootGasLossDisturb_pft=0._r8
  allocate(RootCUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootCUlmNutUptake_pvr=0._r8
  allocate(RootOUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootOUlmNutUptake_pvr=0._r8
  allocate(RCO2N_pvr(jroots,JZ,JP,JY,JX));RCO2N_pvr=0._r8
  allocate(RUNNHP(jroots,JZ,JP,JY,JX));RUNNHP=0._r8
  allocate(RUNNOP(jroots,JZ,JP,JY,JX));RUNNOP=0._r8
  allocate(RUNNBP(jroots,JZ,JP,JY,JX));RUNNBP=0._r8
  allocate(RUNNXP(jroots,JZ,JP,JY,JX));RUNNXP=0._r8
  allocate(RNH3Z(JP,JY,JX));    RNH3Z=0._r8
  allocate(NH3Dep2_brch(MaxNumBranches,JP,JY,JX)); NH3Dep2_brch=0._r8
  allocate(RAutoRootO2Limter_pvr(jroots,JZ,JP,JY,JX)); RAutoRootO2Limter_pvr=0._r8
  allocate(RH2GZ(JP,JY,JX));    RH2GZ=0._r8
  allocate(PlantRootSoilElmNetX_pft(NumPlantChemElms,JP,JY,JX));   PlantRootSoilElmNetX_pft=0._r8
  allocate(PlantExudChemElmCum_pft(NumPlantChemElms,JP,JY,JX));   PlantExudChemElmCum_pft=0._r8
  allocate(GridPlantRootH2OUptake_vr(0:JZ,JY,JX)); GridPlantRootH2OUptake_vr=0._r8
  allocate(THeatRootUptake(0:JZ,JY,JX));  THeatRootUptake=0._r8
  allocate(trcg_air2root_flx_vr(idg_beg:idg_end-1,JZ,JY,JX));   trcg_air2root_flx_vr=0._r8
  allocate(trcg_TLP(idg_beg:idg_end-1,JZ,JY,JX));   trcg_TLP=0._r8
  allocate(trcs_plant_uptake_vr(ids_beg:ids_end,JZ,JY,JX));    trcs_plant_uptake_vr=0._r8
  allocate(TDFOME(NumPlantChemElms,1:jcplx,JZ,JY,JX));TDFOME=0._r8
  allocate(TCO2P(JZ,JY,JX));    TCO2P=0._r8
  allocate(TUPOXP(JZ,JY,JX));   TUPOXP=0._r8
  allocate(RTDNT(JZ,JY,JX));    RTDNT=0._r8
  allocate(ROXYX(0:JZ,JY,JX));  ROXYX=0._r8
  allocate(ROXYY(0:JZ,JY,JX));  ROXYY=0._r8
  allocate(RNH4X(0:JZ,JY,JX));  RNH4X=0._r8
  allocate(RNH4Y(0:JZ,JY,JX));  RNH4Y=0._r8
  allocate(RNO3X(0:JZ,JY,JX));  RNO3X=0._r8
  allocate(RNO3Y(0:JZ,JY,JX));  RNO3Y=0._r8
  allocate(RNO2X(0:JZ,JY,JX));  RNO2X=0._r8
  allocate(RNO2Y(0:JZ,JY,JX));  RNO2Y=0._r8
  allocate(RPO4X(0:JZ,JY,JX));  RPO4X=0._r8
  allocate(RPO4Y(0:JZ,JY,JX));  RPO4Y=0._r8
  allocate(RN2OX(0:JZ,JY,JX));  RN2OX=0._r8
  allocate(RN2OY(0:JZ,JY,JX));  RN2OY=0._r8
  allocate(RNHBX(0:JZ,JY,JX));  RNHBX=0._r8
  allocate(RNHBY(0:JZ,JY,JX));  RNHBY=0._r8
  allocate(RN3BX(0:JZ,JY,JX));  RN3BX=0._r8
  allocate(RN3BY(0:JZ,JY,JX));  RN3BY=0._r8
  allocate(RN2BX(0:JZ,JY,JX));  RN2BX=0._r8
  allocate(RN2BY(0:JZ,JY,JX));  RN2BY=0._r8
  allocate(RPOBX(0:JZ,JY,JX));  RPOBX=0._r8
  allocate(RPOBY(0:JZ,JY,JX));  RPOBY=0._r8
  allocate(ROQCX(1:jcplx,0:JZ,JY,JX));ROQCX=0._r8
  allocate(ROQCY(1:jcplx,0:JZ,JY,JX));ROQCY=0._r8
  allocate(ROQAX(1:jcplx,0:JZ,JY,JX));ROQAX=0._r8
  allocate(ROQAY(1:jcplx,0:JZ,JY,JX));ROQAY=0._r8
  allocate(TH2GZ(JY,JX));       TH2GZ=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantRates
  use abortutils, only : destroy
  implicit none
  call destroy(Eco_NEE_col)
  call destroy(NH3Dep2Can_pft)
  call destroy(NH3EmiCum_pft)
  call destroy(SurfLitrfalStrutElms_pft)
  call destroy(RootMycoExudElm_pvr)
  call destroy(RootNutUptake_pvr)
  call destroy(RUPGasSol_vr)
  call destroy(RootN2Fix_pvr)
  call destroy(RUPP2P)
  call destroy(RUPP2B)
  call destroy(RUPP1P)
  call destroy(RUPP1B)
  call destroy(LeafElmntRemobFlx_brch)
  call destroy(PetioleChemElmRemobFlx_brch)
  call destroy(GrossCO2Fix_pft)
  call destroy(LitrfalStrutElms_pft)
  call destroy(PlantN2FixCum_pft)
  call destroy(GrossResp_pft)
  call destroy(ElmBalanceCum_pft)
  call destroy(LitrfalStrutElms_pft)
  call destroy(LitfalStrutElms_pvr)
  call destroy(NetPrimProduct_pft)
  call destroy(ETCanopy_pft)
  call destroy(CanopyPlusNoduRespC_pft)
  call destroy(EcoHavstElmnt_pft)
  call destroy(EcoHavstElmntCum_pft)
  call destroy(CO2ByFire_pft)
  call destroy(CH4ByFire_pft)
  call destroy(O2ByFire_pft)
  call destroy(NH3byFire_pft)
  call destroy(N2ObyFire_pft)
  call destroy(PO4byFire_pft)
  call destroy(ROXYP)  
  call destroy(RCO2P)
  call destroy(RUPOXP)
  call destroy(RootRespPotent_pvr)
  call destroy(RCO2A_pvr)
  call destroy(RootMycoExudElms_pft)
  call destroy(RootNH4Uptake_pft)
  call destroy(RootNO3Uptake_pft)
  call destroy(RootH2PO4Uptake_pft)
  call destroy(RootHPO4Uptake_pft)
  call destroy(RootN2Fix_pft)
  call destroy(RootOUlmNutUptake_pvr)
  call destroy(RootCUlmNutUptake_pvr)
  call destroy(RCO2N_pvr)
  call destroy(RUNNHP)
  call destroy(RUNNOP)
  call destroy(RUNNBP)
  call destroy(RUNNXP)
  call destroy(RNH3Z)
  call destroy(NH3Dep2_brch)
  call destroy(RAutoRootO2Limter_pvr)
  call destroy(RH2GZ)
  call destroy(PlantRootSoilElmNetX_pft)
  call destroy(PlantExudChemElmCum_pft)
  call destroy(GridPlantRootH2OUptake_vr)
  call destroy(THeatRootUptake)
  call destroy(TDFOME)
  call destroy(TCO2P)
  call destroy(TUPOXP)
  call destroy(RTDNT)
  call destroy(ROXYX)
  call destroy(ROXYY)
  call destroy(RNH4X)
  call destroy(RNH4Y)
  call destroy(RNO3X)
  call destroy(RNO3Y)
  call destroy(RNO2X)
  call destroy(RNO2Y)
  call destroy(RPO4X)
  call destroy(RPO4Y)
  call destroy(RN2OX)
  call destroy(RN2OY)
  call destroy(RNHBX)
  call destroy(RNHBY)
  call destroy(RN3BX)
  call destroy(RN3BY)
  call destroy(RN2BX)
  call destroy(RN2BY)
  call destroy(RPOBX)
  call destroy(RPOBY)
  call destroy(ROQCX)
  call destroy(ROQCY)
  call destroy(ROQAX)
  call destroy(ROQAY)
  call destroy(TH2GZ)
  call destroy(trcs_plant_uptake_vr)
  end subroutine DestructPlantRates

end module PlantDataRateType
