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
  real(r8),target,allocatable ::  RNH3C(:,:,:)                       !canopy NH3 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TNH3C(:,:,:)                       !total canopy NH3 flux, [g d-2 ]
  real(r8),target,allocatable ::  SurfLitrfallChemElmnts_pft(:,:,:,:)                     !total surface litterfall element, [g d-2]
  real(r8),target,allocatable ::  RDFOME(:,:,:,:,:,:,:)                !root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPNH4(:,:,:,:,:)                  !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPNHB(:,:,:,:,:)                  !root uptake of NH4 band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPNO3(:,:,:,:,:)                  !root uptake of NO3 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPNOB(:,:,:,:,:)                  !root uptake of NO3 band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPH2P(:,:,:,:,:)                  !root uptake of PO4 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPH2B(:,:,:,:,:)                  !root uptake of PO4 band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPNF(:,:,:,:)                     !root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPHGS(:,:,:,:,:)                  !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPH1P(:,:,:,:,:)                  !root uptake HPO4 non-band
  real(r8),target,allocatable ::  RUPH1B(:,:,:,:,:)                  !root uptake of HPO4 band
  real(r8),target,allocatable ::  RUPP2P(:,:,:,:,:)                  !root uptake of H2PO4 non-band
  real(r8),target,allocatable ::  RUPP2B(:,:,:,:,:)                  !root uptake of H2PO4 band
  real(r8),target,allocatable ::  RUPP1P(:,:,:,:,:)                  !HPO4 demand in non-band by each root population
  real(r8),target,allocatable ::  RUPP1B(:,:,:,:,:)                  !HPO4 demand in band by each root population
  real(r8),target,allocatable ::  RCELX(:,:,:,:,:)                   !element translocated from leaf during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCZSX(:,:,:,:)                     !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCPSX(:,:,:,:)                     !P translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  RCESX(:,:,:,:,:)                   !element translocated from sheath during senescence, [g d-2 h-1]
  real(r8),target,allocatable ::  GrossCO2Fix_pft(:,:,:)                       !total gross CO2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  LitrfallChemElmnts_pft(:,:,:,:)                     !total plant element litterfall , [g d-2 ]
  real(r8),target,allocatable ::  TZUPFX(:,:,:)                      !total plant N2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  GrossResp_pft(:,:,:)                       !total plant respiration, [g d-2 ]
  real(r8),target,allocatable ::  BALE(:,:,:,:)                      !plant element balance, [g d-2]
  real(r8),target,allocatable ::  HESNC(:,:,:,:)                     !plant element litterfall, [g d-2 h-1]
  real(r8),target,allocatable ::  ESNC(:,:,:,:,:,:,:)                !plant litterfall element, [g d-2 h-1]
  real(r8),target,allocatable ::  NetPrimaryProductvity_pft(:,:,:)                        !total net primary productivity, [g d-2]
  real(r8),target,allocatable ::  ETCanP(:,:,:)                       !total transpiration, [m d-2], <0 into atmosphere
  real(r8),target,allocatable ::  TCO2A(:,:,:)                       !total autotrophic respiration, [g d-2 ]
  real(r8),target,allocatable ::  HVSTE(:,:,:,:)                     !plant element harvest, [g d-2 ]
  real(r8),target,allocatable ::  THVSTE(:,:,:,:)                    !total plant harvest, [g d-2 ]
  real(r8),target,allocatable ::  CO2ByFire_pft(:,:,:)                       !plant CO2 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  CH4ByFire_pft(:,:,:)                       !plant CH4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VOXYF(:,:,:)                       !plant O2 uptake from fire, [g d-2 ]
  real(r8),target,allocatable ::  VNH3F(:,:,:)                       !plant NH3 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VN2OF(:,:,:)                       !plant N2O emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VPO4F(:,:,:)                       !plant PO4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  ROXYP(:,:,:,:,:)                   !root  O2 demand from respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_air2root_flx_pft_vr(:,:,:,:,:,:)                  !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_Root_DisEvap_flx_vr(:,:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2S(:,:,:,:,:)                   !aqueous CO2 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPOXS(:,:,:,:,:)                  !aqueous O2 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPCHS(:,:,:,:,:)                  !aqueous CH4 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPN2S(:,:,:,:,:)                  !aqueous N2O flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPN3S(:,:,:,:,:)                  !aqueous NH3 flux from roots to soil water non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPN3B(:,:,:,:,:)                  !aqueous NH3 flux from roots to soil water band, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2P(:,:,:,:,:)                   !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RUPOXP(:,:,:,:,:)                  !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2M(:,:,:,:,:)                   !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2A(:,:,:,:,:)                   !root respiration constrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  UPOME(:,:,:,:)                     !total root uptake (+ve) - exudation (-ve) of dissovled element, [g d-2 h-1]
  real(r8),target,allocatable ::  UPNH4(:,:,:)                       !total root uptake of NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  UPNO3(:,:,:)                       !total root uptake of NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  UPH2P(:,:,:)                       !total root uptake of PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  UPH1P(:,:,:)                       !total root uptake of HPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  UPNF(:,:,:)                        !total root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RootGasLoss_disturb(:,:,:,:)                !gas flux from root disturbance [g d-2 h-1]
  real(r8),target,allocatable ::  RUONH4(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RUONHB(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RUONO3(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RUONOB(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RUOH2P(:,:,:,:,:)                  !root uptake of PO4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RUOH2B(:,:,:,:,:)                  !root uptake of PO4 band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RUCNH4(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUCNHB(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUCNO3(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUCNOB(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUCH2P(:,:,:,:,:)                  !root uptake of PO4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUCH2B(:,:,:,:,:)                  !root uptake of PO4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUOH1P(:,:,:,:,:)                  !root HPO4 uptake in non-band unlimited by O2
  real(r8),target,allocatable ::  RUCH1P(:,:,:,:,:)                  !root HPO4 uptake in non-band unlimited by nonstructural C
  real(r8),target,allocatable ::  RUOH1B(:,:,:,:,:)                  !root HPO4 uptake in band unlimited by O2
  real(r8),target,allocatable ::  RUCH1B(:,:,:,:,:)                  !root HPO4 uptake in band unlimited by nonstructural C
  real(r8),target,allocatable ::  RCO2N(:,:,:,:,:)                   !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNHP(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNOP(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNBP(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RUNNXP(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH3Z(:,:,:)                       !gaseous NH3 flux fron root disturbance non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH3B(:,:,:,:)                     !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8),target,allocatable ::  WFR(:,:,:,:,:)                     !O2 constraint to root respiration, []
  real(r8),target,allocatable ::  RH2GZ(:,:,:)                       !gaseous H2 flux fron root disturbance, [g d-2 h-1]
  real(r8),target,allocatable ::  HEUPTK(:,:,:,:)                    !net root element uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  PlantExudChemElmnts_pft(:,:,:,:)                    !total net root element uptake (+ve) - exudation (-ve), [g d-2 ]
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
  allocate(RNH3C(JP,JY,JX));    RNH3C=0._r8
  allocate(TNH3C(JP,JY,JX));    TNH3C=0._r8
  allocate(SurfLitrfallChemElmnts_pft(NumOfPlantChemElements,JP,JY,JX));    SurfLitrfallChemElmnts_pft=0._r8
  allocate(RDFOME(NumOfPlantChemElements,2,1:jcplx,JZ,JP,JY,JX));RDFOME=0._r8
  allocate(RUPNH4(jroots,JZ,JP,JY,JX));RUPNH4=0._r8
  allocate(RUPNHB(jroots,JZ,JP,JY,JX));RUPNHB=0._r8
  allocate(RUPNO3(jroots,JZ,JP,JY,JX));RUPNO3=0._r8
  allocate(RUPNOB(jroots,JZ,JP,JY,JX));RUPNOB=0._r8
  allocate(RUPH2P(jroots,JZ,JP,JY,JX));RUPH2P=0._r8
  allocate(RUPH2B(jroots,JZ,JP,JY,JX));RUPH2B=0._r8
  allocate(RUPNF(JZ,JP,JY,JX)); RUPNF=0._r8
  allocate(RUPHGS(jroots,JZ,JP,JY,JX));RUPHGS=0._r8
  allocate(RUPH1P(jroots,JZ,JP,JY,JX));RUPH1P=0._r8
  allocate(RUPH1B(jroots,JZ,JP,JY,JX));RUPH1B=0._r8
  allocate(RUPP2P(jroots,JZ,JP,JY,JX));RUPP2P=0._r8
  allocate(RUPP2B(jroots,JZ,JP,JY,JX));RUPP2B=0._r8
  allocate(RUPP1P(jroots,JZ,JP,JY,JX));RUPP1P=0._r8
  allocate(RUPP1B(jroots,JZ,JP,JY,JX));RUPP1B=0._r8
  allocate(RCELX(NumOfPlantChemElements,MaxNumBranches,JP,JY,JX)); RCELX=0._r8
  allocate(RCESX(NumOfPlantChemElements,MaxNumBranches,JP,JY,JX)); RCESX=0._r8
  allocate(GrossCO2Fix_pft(JP,JY,JX));    GrossCO2Fix_pft=0._r8
  allocate(LitrfallChemElmnts_pft(NumOfPlantChemElements,JP,JY,JX));    LitrfallChemElmnts_pft=0._r8
  allocate(TZUPFX(JP,JY,JX));   TZUPFX=0._r8
  allocate(GrossResp_pft(JP,JY,JX));    GrossResp_pft=0._r8
  allocate(BALE(NumOfPlantChemElements,JP,JY,JX));     BALE=0._r8
  allocate(HESNC(NumOfPlantChemElements,JP,JY,JX));    HESNC=0._r8
  allocate(ESNC(NumOfPlantChemElements,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JP,JY,JX));ESNC=0._r8
  allocate(NetPrimaryProductvity_pft(JP,JY,JX));     NetPrimaryProductvity_pft=0._r8
  allocate(ETCanP(JP,JY,JX));    ETCanP=0._r8
  allocate(TCO2A(JP,JY,JX));    TCO2A=0._r8
  allocate(HVSTE(NumOfPlantChemElements,JP,JY,JX));    HVSTE=0._r8
  allocate(THVSTE(NumOfPlantChemElements,JP,JY,JX));   THVSTE=0._r8
  allocate(CO2ByFire_pft(JP,JY,JX));    CO2ByFire_pft=0._r8
  allocate(CH4ByFire_pft(JP,JY,JX));    CH4ByFire_pft=0._r8
  allocate(VOXYF(JP,JY,JX));    VOXYF=0._r8
  allocate(VNH3F(JP,JY,JX));    VNH3F=0._r8
  allocate(VN2OF(JP,JY,JX));    VN2OF=0._r8
  allocate(VPO4F(JP,JY,JX));    VPO4F=0._r8
  allocate(ROXYP(jroots,JZ,JP,JY,JX));ROXYP=0._r8
  allocate(trcg_air2root_flx_pft_vr(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_air2root_flx_pft_vr=0._r8
  allocate(trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_Root_DisEvap_flx_vr=0._r8
  allocate(RCO2S(jroots,JZ,JP,JY,JX));RCO2S=0._r8
  allocate(RUPOXS(jroots,JZ,JP,JY,JX));RUPOXS=0._r8
  allocate(RUPCHS(jroots,JZ,JP,JY,JX));RUPCHS=0._r8
  allocate(RUPN2S(jroots,JZ,JP,JY,JX));RUPN2S=0._r8
  allocate(RUPN3S(jroots,JZ,JP,JY,JX));RUPN3S=0._r8
  allocate(RUPN3B(jroots,JZ,JP,JY,JX));RUPN3B=0._r8
  allocate(RCO2P(jroots,JZ,JP,JY,JX));RCO2P=0._r8
  allocate(RUPOXP(jroots,JZ,JP,JY,JX));RUPOXP=0._r8
  allocate(RCO2M(jroots,JZ,JP,JY,JX));RCO2M=0._r8
  allocate(RCO2A(jroots,JZ,JP,JY,JX));RCO2A=0._r8
  allocate(UPOME(1:NumOfPlantChemElements,JP,JY,JX));    UPOME=0._r8
  allocate(UPNH4(JP,JY,JX));    UPNH4=0._r8
  allocate(UPNO3(JP,JY,JX));    UPNO3=0._r8
  allocate(UPH2P(JP,JY,JX));    UPH2P=0._r8
  allocate(UPH1P(JP,JY,JX));    UPH1P=0._r8
  allocate(UPNF(JP,JY,JX));     UPNF=0._r8
  allocate(RootGasLoss_disturb(idg_beg:idg_end-1,JP,JY,JX)); RootGasLoss_disturb=0._r8
  allocate(RUONH4(jroots,JZ,JP,JY,JX));RUONH4=0._r8
  allocate(RUONHB(jroots,JZ,JP,JY,JX));RUONHB=0._r8
  allocate(RUONO3(jroots,JZ,JP,JY,JX));RUONO3=0._r8
  allocate(RUONOB(jroots,JZ,JP,JY,JX));RUONOB=0._r8
  allocate(RUOH2P(jroots,JZ,JP,JY,JX));RUOH2P=0._r8
  allocate(RUOH2B(jroots,JZ,JP,JY,JX));RUOH2B=0._r8
  allocate(RUCNH4(jroots,JZ,JP,JY,JX));RUCNH4=0._r8
  allocate(RUCNHB(jroots,JZ,JP,JY,JX));RUCNHB=0._r8
  allocate(RUCNO3(jroots,JZ,JP,JY,JX));RUCNO3=0._r8
  allocate(RUCNOB(jroots,JZ,JP,JY,JX));RUCNOB=0._r8
  allocate(RUCH2P(jroots,JZ,JP,JY,JX));RUCH2P=0._r8
  allocate(RUCH2B(jroots,JZ,JP,JY,JX));RUCH2B=0._r8
  allocate(RUOH1P(jroots,JZ,JP,JY,JX));RUOH1P=0._r8
  allocate(RUCH1P(jroots,JZ,JP,JY,JX));RUCH1P=0._r8
  allocate(RUOH1B(jroots,JZ,JP,JY,JX));RUOH1B=0._r8
  allocate(RUCH1B(jroots,JZ,JP,JY,JX));RUCH1B=0._r8
  allocate(RCO2N(jroots,JZ,JP,JY,JX));RCO2N=0._r8
  allocate(RUNNHP(jroots,JZ,JP,JY,JX));RUNNHP=0._r8
  allocate(RUNNOP(jroots,JZ,JP,JY,JX));RUNNOP=0._r8
  allocate(RUNNBP(jroots,JZ,JP,JY,JX));RUNNBP=0._r8
  allocate(RUNNXP(jroots,JZ,JP,JY,JX));RUNNXP=0._r8
  allocate(RNH3Z(JP,JY,JX));    RNH3Z=0._r8
  allocate(RNH3B(MaxNumBranches,JP,JY,JX)); RNH3B=0._r8
  allocate(WFR(2,JZ,JP,JY,JX)); WFR=0._r8
  allocate(RH2GZ(JP,JY,JX));    RH2GZ=0._r8
  allocate(HEUPTK(NumOfPlantChemElements,JP,JY,JX));   HEUPTK=0._r8
  allocate(PlantExudChemElmnts_pft(NumOfPlantChemElements,JP,JY,JX));   PlantExudChemElmnts_pft=0._r8
  allocate(GridPlantRootH2OUptake_vr(0:JZ,JY,JX)); GridPlantRootH2OUptake_vr=0._r8
  allocate(THeatRootUptake(0:JZ,JY,JX));  THeatRootUptake=0._r8
  allocate(trcg_air2root_flx_vr(idg_beg:idg_end-1,JZ,JY,JX));   trcg_air2root_flx_vr=0._r8
  allocate(trcg_TLP(idg_beg:idg_end-1,JZ,JY,JX));   trcg_TLP=0._r8
  allocate(trcs_plant_uptake_vr(ids_beg:ids_end,JZ,JY,JX));    trcs_plant_uptake_vr=0._r8
  allocate(TDFOME(NumOfPlantChemElements,1:jcplx,JZ,JY,JX));TDFOME=0._r8
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
  call destroy(RNH3C)
  call destroy(TNH3C)
  call destroy(SurfLitrfallChemElmnts_pft)
  call destroy(RDFOME)
  call destroy(RUPNH4)
  call destroy(RUPNHB)
  call destroy(RUPNO3)
  call destroy(RUPNOB)
  call destroy(RUPH2P)
  call destroy(RUPH2B)
  call destroy(RUPNF)
  call destroy(RUPHGS)
  call destroy(RUPH1P)
  call destroy(RUPH1B)
  call destroy(RUPP2P)
  call destroy(RUPP2B)
  call destroy(RUPP1P)
  call destroy(RUPP1B)
  call destroy(RCELX)
  call destroy(RCESX)
  call destroy(GrossCO2Fix_pft)
  call destroy(LitrfallChemElmnts_pft)
  call destroy(TZUPFX)
  call destroy(GrossResp_pft)
  call destroy(BALE)
  call destroy(HESNC)
  call destroy(ESNC)
  call destroy(NetPrimaryProductvity_pft)
  call destroy(ETCanP)
  call destroy(TCO2A)
  call destroy(HVSTE)
  call destroy(THVSTE)
  call destroy(CO2ByFire_pft)
  call destroy(CH4ByFire_pft)
  call destroy(VOXYF)
  call destroy(VNH3F)
  call destroy(VN2OF)
  call destroy(VPO4F)
  call destroy(ROXYP)
  call destroy(RCO2S)
  call destroy(RUPOXS)
  call destroy(RUPCHS)
  call destroy(RUPN2S)
  call destroy(RUPN3S)
  call destroy(RUPN3B)
  call destroy(RCO2P)
  call destroy(RUPOXP)
  call destroy(RCO2M)
  call destroy(RCO2A)
  call destroy(UPOME)
  call destroy(UPNH4)
  call destroy(UPNO3)
  call destroy(UPH2P)
  call destroy(UPH1P)
  call destroy(UPNF)
  call destroy(RUONH4)
  call destroy(RUONHB)
  call destroy(RUONO3)
  call destroy(RUONOB)
  call destroy(RUOH2P)
  call destroy(RUOH2B)
  call destroy(RUCNH4)
  call destroy(RUCNHB)
  call destroy(RUCNO3)
  call destroy(RUCNOB)
  call destroy(RUCH2P)
  call destroy(RUCH2B)
  call destroy(RUOH1P)
  call destroy(RUCH1P)
  call destroy(RUOH1B)
  call destroy(RUCH1B)
  call destroy(RCO2N)
  call destroy(RUNNHP)
  call destroy(RUNNOP)
  call destroy(RUNNBP)
  call destroy(RUNNXP)
  call destroy(RNH3Z)
  call destroy(RNH3B)
  call destroy(WFR)
  call destroy(RH2GZ)
  call destroy(HEUPTK)
  call destroy(PlantExudChemElmnts_pft)
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
