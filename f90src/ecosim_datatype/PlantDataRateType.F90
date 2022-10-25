module PlantDataRateType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use ElmIDMod
  use EcoSIMConfig, only : jcplx1=> jcplx1c,jsken=>jskenc
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  TCNET(:,:)                         !total canopy net CO2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  RNH3C(:,:,:)                       !canopy NH3 flux, [g d-2 h-1]
  real(r8),allocatable ::  TNH3C(:,:,:)                       !total canopy NH3 flux, [g d-2 ]
  real(r8),allocatable ::  TESN0(:,:,:,:)                     !total surface litterfall element, [g d-2]
  real(r8),allocatable ::  RDFOMC(:,:,:,:,:,:)                !root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8),allocatable ::  RDFOMN(:,:,:,:,:,:)                !root uptake (+ve) - exudation (-ve) of DON, [g d-2 h-1]
  real(r8),allocatable ::  RDFOMP(:,:,:,:,:,:)                !root uptake (+ve) - exudation (-ve) of DOP, [g d-2 h-1]
  real(r8),allocatable ::  RUPNH4(:,:,:,:,:)                  !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8),allocatable ::  RUPNHB(:,:,:,:,:)                  !root uptake of NH4 band, [g d-2 h-1]
  real(r8),allocatable ::  RUPNO3(:,:,:,:,:)                  !root uptake of NO3 non-band, [g d-2 h-1]
  real(r8),allocatable ::  RUPNOB(:,:,:,:,:)                  !root uptake of NO3 band, [g d-2 h-1]
  real(r8),allocatable ::  RUPH2P(:,:,:,:,:)                  !root uptake of PO4 non-band, [g d-2 h-1]
  real(r8),allocatable ::  RUPH2B(:,:,:,:,:)                  !root uptake of PO4 band, [g d-2 h-1]
  real(r8),allocatable ::  RUPNF(:,:,:,:)                     !root N2 fixation, [g d-2 h-1]
  real(r8),allocatable ::  RUPHGS(:,:,:,:,:)                  !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8),allocatable ::  RUPH1P(:,:,:,:,:)                  !root uptake HPO4 non-band
  real(r8),allocatable ::  RUPH1B(:,:,:,:,:)                  !root uptake of HPO4 band
  real(r8),allocatable ::  RUPP2P(:,:,:,:,:)                  !root uptake of H2PO4 non-band
  real(r8),allocatable ::  RUPP2B(:,:,:,:,:)                  !root uptake of H2PO4 band
  real(r8),allocatable ::  RUPP1P(:,:,:,:,:)                  !HPO4 demand in non-band by each root population
  real(r8),allocatable ::  RUPP1B(:,:,:,:,:)                  !HPO4 demand in band by each root population
  real(r8),allocatable ::  RCZLX(:,:,:,:)                     !N translocated from leaf during senescence, [g d-2 h-1]
  real(r8),allocatable ::  RCPLX(:,:,:,:)                     !P translocated from leaf during senescence, [g d-2 h-1]
  real(r8),allocatable ::  RCCLX(:,:,:,:)                     !C translocated from leaf during senescence, [g d-2 h-1]
  real(r8),allocatable ::  RCZSX(:,:,:,:)                     !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8),allocatable ::  RCPSX(:,:,:,:)                     !P translocated from sheath during senescence, [g d-2 h-1]
  real(r8),allocatable ::  RCCSX(:,:,:,:)                     !C translocated from sheath during senescence, [g d-2 h-1]
  real(r8),allocatable ::  CARBN(:,:,:)                       !total gross CO2 fixation, [g d-2 ]
  real(r8),allocatable ::  TESNC(:,:,:,:)                     !total plant element litterfall , [g d-2 ]
  real(r8),allocatable ::  TZUPFX(:,:,:)                      !total plant N2 fixation, [g d-2 ]
  real(r8),allocatable ::  TCO2T(:,:,:)                       !total plant respiration, [g d-2 ]
  real(r8),allocatable ::  BALE(:,:,:,:)                      !plant element balance, [g d-2]
  real(r8),allocatable ::  HESNC(:,:,:,:)                     !plant element litterfall, [g d-2 h-1]
  real(r8),allocatable ::  ESNC(:,:,:,:,:,:,:)                !plant litterfall element, [g d-2 h-1]
  real(r8),allocatable ::  ZNPP(:,:,:)                        !total net primary productivity, [g d-2]
  real(r8),allocatable ::  CTRAN(:,:,:)                       !total transpiration, [m d-2]
  real(r8),allocatable ::  TCO2A(:,:,:)                       !total autotrophic respiration, [g d-2 ]
  real(r8),allocatable ::  HVSTE(:,:,:,:)                     !plant element harvest, [g d-2 ]
  real(r8),allocatable ::  THVSTE(:,:,:,:)                    !total plant harvest, [g d-2 ]
  real(r8),allocatable ::  VCO2F(:,:,:)                       !plant CO2 emission from fire, [g d-2 ]
  real(r8),allocatable ::  VCH4F(:,:,:)                       !plant CH4 emission from fire, [g d-2 ]
  real(r8),allocatable ::  VOXYF(:,:,:)                       !plant O2 uptake from fire, [g d-2 ]
  real(r8),allocatable ::  VNH3F(:,:,:)                       !plant NH3 emission from fire, [g d-2 ]
  real(r8),allocatable ::  VN2OF(:,:,:)                       !plant N2O emission from fire, [g d-2 ]
  real(r8),allocatable ::  VPO4F(:,:,:)                       !plant PO4 emission from fire, [g d-2 ]
  real(r8),allocatable ::  ROXYP(:,:,:,:,:)                   !root  O2 demand from respiration, [g d-2 h-1]
  real(r8),allocatable ::  RCOFLA(:,:,:,:,:)                  !gaseous CO2 flux through roots, [g d-2 h-1]
  real(r8),allocatable ::  ROXFLA(:,:,:,:,:)                  !gaseous O2 flux through roots, [g d-2 h-1]
  real(r8),allocatable ::  RCHFLA(:,:,:,:,:)                  !gaseous CH4 flux through roots, [g d-2 h-1]
  real(r8),allocatable ::  RN2FLA(:,:,:,:,:)                  !gaseous N2O flux through roots, [g d-2 h-1]
  real(r8),allocatable ::  RNHFLA(:,:,:,:,:)                  !gaseous NH3 flux through roots, [g d-2 h-1]
  real(r8),allocatable ::  RCODFA(:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) CO2 flux in roots, [g d-2 h-1]
  real(r8),allocatable ::  ROXDFA(:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) O2 flux in roots, [g d-2 h-1]
  real(r8),allocatable ::  RCHDFA(:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) CH4 flux in roots, [g d-2 h-1]
  real(r8),allocatable ::  RN2DFA(:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) N2O flux in roots, [g d-2 h-1]
  real(r8),allocatable ::  RNHDFA(:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) NH3 flux in roots, [g d-2 h-1]
  real(r8),allocatable ::  RCO2S(:,:,:,:,:)                   !aqueous CO2 flux from roots to soil water, [g d-2 h-1]
  real(r8),allocatable ::  RUPOXS(:,:,:,:,:)                  !aqueous O2 flux from roots to soil water, [g d-2 h-1]
  real(r8),allocatable ::  RUPCHS(:,:,:,:,:)                  !aqueous CH4 flux from roots to soil water, [g d-2 h-1]
  real(r8),allocatable ::  RUPN2S(:,:,:,:,:)                  !aqueous N2O flux from roots to soil water, [g d-2 h-1]
  real(r8),allocatable ::  RUPN3S(:,:,:,:,:)                  !aqueous NH3 flux from roots to soil water non-band, [g d-2 h-1]
  real(r8),allocatable ::  RUPN3B(:,:,:,:,:)                  !aqueous NH3 flux from roots to soil water band, [g d-2 h-1]
  real(r8),allocatable ::  RCO2P(:,:,:,:,:)                   !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8),allocatable ::  RUPOXP(:,:,:,:,:)                  !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8),allocatable ::  RCO2M(:,:,:,:,:)                   !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RCO2A(:,:,:,:,:)                   !root respiration constrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  UPOME(:,:,:,:)                     !total root uptake (+ve) - exudation (-ve) of dissovled element, [g d-2 h-1]
  real(r8),allocatable ::  UPNH4(:,:,:)                       !total root uptake of NH4, [g d-2 h-1]
  real(r8),allocatable ::  UPNO3(:,:,:)                       !total root uptake of NO3, [g d-2 h-1]
  real(r8),allocatable ::  UPH2P(:,:,:)                       !total root uptake of PO4, [g d-2 h-1]
  real(r8),allocatable ::  UPH1P(:,:,:)                       !total root uptake of HPO4, [g d-2 h-1]
  real(r8),allocatable ::  UPNF(:,:,:)                        !total root N2 fixation, [g d-2 h-1]
  real(r8),allocatable ::  RCO2Z(:,:,:)                       !gaseous CO2 flux fron root disturbance, [g d-2 h-1]
  real(r8),allocatable ::  ROXYZ(:,:,:)                       !gaseous O2 flux fron root disturbance, [g d-2 h-1]
  real(r8),allocatable ::  RCH4Z(:,:,:)                       !gaseous CH4 flux fron root disturbance, [g d-2 h-1]
  real(r8),allocatable ::  RN2OZ(:,:,:)                       !gaseous N2O flux fron root disturbance, [g d-2 h-1]
  real(r8),allocatable ::  RUONH4(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RUONHB(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RUONO3(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RUONOB(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RUOH2P(:,:,:,:,:)                  !root uptake of PO4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RUOH2B(:,:,:,:,:)                  !root uptake of PO4 band unconstrained by O2, [g d-2 h-1]
  real(r8),allocatable ::  RUCNH4(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUCNHB(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUCNO3(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUCNOB(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUCH2P(:,:,:,:,:)                  !root uptake of PO4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUCH2B(:,:,:,:,:)                  !root uptake of PO4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUOH1P(:,:,:,:,:)                  !root HPO4 uptake in non-band unlimited by O2
  real(r8),allocatable ::  RUCH1P(:,:,:,:,:)                  !root HPO4 uptake in non-band unlimited by nonstructural C
  real(r8),allocatable ::  RUOH1B(:,:,:,:,:)                  !root HPO4 uptake in band unlimited by O2
  real(r8),allocatable ::  RUCH1B(:,:,:,:,:)                  !root HPO4 uptake in band unlimited by nonstructural C
  real(r8),allocatable ::  RCO2N(:,:,:,:,:)                   !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),allocatable ::  RUNNHP(:,:,:,:,:)                  !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8),allocatable ::  RUNNOP(:,:,:,:,:)                  !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8),allocatable ::  RUNNBP(:,:,:,:,:)                  !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8),allocatable ::  RUNNXP(:,:,:,:,:)                  !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8),allocatable ::  RNH3Z(:,:,:)                       !gaseous NH3 flux fron root disturbance non-band, [g d-2 h-1]
  real(r8),allocatable ::  RNH3B(:,:,:,:)                     !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8),allocatable ::  WFR(:,:,:,:,:)                     !O2 constraint to root respiration, []
  real(r8),allocatable ::  RHGFLA(:,:,:,:,:)                  !gaseous H2 flux through roots, [g d-2 h-1]
  real(r8),allocatable ::  RHGDFA(:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) H2 flux in roots, [g d-2 h-1]
  real(r8),allocatable ::  RH2GZ(:,:,:)                       !gaseous H2 flux fron root disturbance, [g d-2 h-1]
  real(r8),allocatable ::  HCUPTK(:,:,:)                      !net root C uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),allocatable ::  HZUPTK(:,:,:)                      !net root N uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),allocatable ::  HPUPTK(:,:,:)                      !net root P uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),allocatable ::  TEUPTK(:,:,:,:)                    !total net root element uptake (+ve) - exudation (-ve), [g d-2 ]
  real(r8),allocatable ::  TUPWTR(:,:,:)                      !total root water uptake, [m3 d-2]
  real(r8),allocatable ::  TUPHT(:,:,:)                       !total root heat uptake, [MJ d-2]
  real(r8),allocatable ::  TCOFLA(:,:,:)                      !total internal root CO2 flux , [g d-2 h-1]
  real(r8),allocatable ::  TOXFLA(:,:,:)                      !total internal root O2 flux , [g d-2 h-1]
  real(r8),allocatable ::  TCHFLA(:,:,:)                      !total internal root CH4 flux , [g d-2 h-1]
  real(r8),allocatable ::  TN2FLA(:,:,:)                      !total internal root N2O flux , [g d-2 h-1]
  real(r8),allocatable ::  TNHFLA(:,:,:)                      !total internal root NH3 flux , [g d-2 h-1]
  real(r8),allocatable ::  TLCO2P(:,:,:)                      !total root internal CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  TLOXYP(:,:,:)                      !total root internal O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  TLCH4P(:,:,:)                      !total root internal CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  TLN2OP(:,:,:)                      !total root internal N2O flux, [g d-2 h-1]
  real(r8),allocatable ::  TLNH3P(:,:,:)                      !total root internal NH3 flux, [g d-2 h-1]
  real(r8),allocatable ::  TCO2S(:,:,:)                       !total root-soil CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  TUPOXS(:,:,:)                      !total root-soil O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  TUPCHS(:,:,:)                      !total root-soil CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  TUPN2S(:,:,:)                      !total root-soil N2O flux, [g d-2 h-1]
  real(r8),allocatable ::  TUPN3S(:,:,:)                      !total root-soil NH3 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  TUPNH4(:,:,:)                      !total root-soil NH4 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  TUPNO3(:,:,:)                      !total root-soil NO3 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  TUPH2P(:,:,:)                      !total root-soil PO4 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  TUPN3B(:,:,:)                      !total root-soil NH3 flux band, [g d-2 h-1]
  real(r8),allocatable ::  TUPNHB(:,:,:)                      !total root-soil NH4 flux band, [g d-2 h-1]
  real(r8),allocatable ::  TUPNOB(:,:,:)                      !total root-soil NO3 flux band, [g d-2 h-1]
  real(r8),allocatable ::  TUPH2B(:,:,:)                      !total root-soil PO4 flux band, [g d-2 h-1]
  real(r8),allocatable ::  TUPNF(:,:,:)                       !total root N2 fixation, [g d-2 h-1]
  real(r8),allocatable ::  TDFOMC(:,:,:,:)                    !total root C exchange, [g d-2 h-1]
  real(r8),allocatable ::  TDFOMN(:,:,:,:)                    !total root N exchange, [g d-2 h-1]
  real(r8),allocatable ::  TDFOMP(:,:,:,:)                    !total root P exchange, [g d-2 h-1]
  real(r8),allocatable ::  TCO2P(:,:,:)                       !total root CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  TUPOXP(:,:,:)                      !total root internal O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  RTDNT(:,:,:)                       !total root length density, [m m-3]
  real(r8),allocatable ::  TUPH1P(:,:,:)                      !soil-root exch of HPO4 in non-band
  real(r8),allocatable ::  TUPH1B(:,:,:)                      !soil-root exch of HPO4 in band
  real(r8),allocatable ::  ROXYX(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),allocatable ::  ROXYY(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),allocatable ::  RNH4X(:,:,:)                       !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RNH4Y(:,:,:)                       !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RNO3X(:,:,:)                       !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RNO3Y(:,:,:)                       !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RNO2X(:,:,:)                       !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RNO2Y(:,:,:)                       !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RPO4X(:,:,:)                       !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RPO4Y(:,:,:)                       !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),allocatable ::  RN2OX(:,:,:)                       !total root + microbial N2O uptake , [g d-2 h-1]
  real(r8),allocatable ::  RN2OY(:,:,:)                       !total root + microbial N2O uptake, [g d-2 h-1]
  real(r8),allocatable ::  RNHBX(:,:,:)                       !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RNHBY(:,:,:)                       !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RN3BX(:,:,:)                       !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RN3BY(:,:,:)                       !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RN2BX(:,:,:)                       !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RN2BY(:,:,:)                       !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RPOBX(:,:,:)                       !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  RPOBY(:,:,:)                       !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),allocatable ::  ROQCX(:,:,:,:)                     !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),allocatable ::  ROQCY(:,:,:,:)                     !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),allocatable ::  ROQAX(:,:,:,:)                     !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),allocatable ::  ROQAY(:,:,:,:)                     !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),allocatable ::  TH2GZ(:,:)                         !total root H2 flux, [g d-2]
  real(r8),allocatable ::  TUPHGS(:,:,:)                      !total root-soil H2 flux, [g d-2 h-1]
  real(r8),allocatable ::  THGFLA(:,:,:)                      !total root-atmosphere H2 flux, [g d-2 h-1]
  real(r8),allocatable ::  TLH2GP(:,:,:)                      !total root-soil H2 flux, [g d-2 h-1]
  private :: InitAllocate
  contains

!----------------------------------------------------------------------
  subroutine InitPlantRates
  implicit none

  call InitAllocate

  end subroutine InitPlantRates
!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(TCNET(JY,JX));       TCNET=0._r8
  allocate(RNH3C(JP,JY,JX));    RNH3C=0._r8
  allocate(TNH3C(JP,JY,JX));    TNH3C=0._r8
  allocate(TESN0(npelms,JP,JY,JX));    TESN0=0._r8
  allocate(RDFOMC(2,0:jcplx1,JZ,JP,JY,JX));RDFOMC=0._r8
  allocate(RDFOMN(2,0:jcplx1,JZ,JP,JY,JX));RDFOMN=0._r8
  allocate(RDFOMP(2,0:jcplx1,JZ,JP,JY,JX));RDFOMP=0._r8
  allocate(RUPNH4(2,JZ,JP,JY,JX));RUPNH4=0._r8
  allocate(RUPNHB(2,JZ,JP,JY,JX));RUPNHB=0._r8
  allocate(RUPNO3(2,JZ,JP,JY,JX));RUPNO3=0._r8
  allocate(RUPNOB(2,JZ,JP,JY,JX));RUPNOB=0._r8
  allocate(RUPH2P(2,JZ,JP,JY,JX));RUPH2P=0._r8
  allocate(RUPH2B(2,JZ,JP,JY,JX));RUPH2B=0._r8
  allocate(RUPNF(JZ,JP,JY,JX)); RUPNF=0._r8
  allocate(RUPHGS(2,JZ,JP,JY,JX));RUPHGS=0._r8
  allocate(RUPH1P(2,JZ,JP,JY,JX));RUPH1P=0._r8
  allocate(RUPH1B(2,JZ,JP,JY,JX));RUPH1B=0._r8
  allocate(RUPP2P(2,JZ,JP,JY,JX));RUPP2P=0._r8
  allocate(RUPP2B(2,JZ,JP,JY,JX));RUPP2B=0._r8
  allocate(RUPP1P(2,JZ,JP,JY,JX));RUPP1P=0._r8
  allocate(RUPP1B(2,JZ,JP,JY,JX));RUPP1B=0._r8
  allocate(RCZLX(JC,JP,JY,JX)); RCZLX=0._r8
  allocate(RCPLX(JC,JP,JY,JX)); RCPLX=0._r8
  allocate(RCCLX(JC,JP,JY,JX)); RCCLX=0._r8
  allocate(RCZSX(JC,JP,JY,JX)); RCZSX=0._r8
  allocate(RCPSX(JC,JP,JY,JX)); RCPSX=0._r8
  allocate(RCCSX(JC,JP,JY,JX)); RCCSX=0._r8
  allocate(CARBN(JP,JY,JX));    CARBN=0._r8
  allocate(TESNC(npelms,JP,JY,JX));    TESNC=0._r8
  allocate(TZUPFX(JP,JY,JX));   TZUPFX=0._r8
  allocate(TCO2T(JP,JY,JX));    TCO2T=0._r8
  allocate(BALE(npelms,JP,JY,JX));     BALE=0._r8
  allocate(HESNC(npelms,JP,JY,JX));    HESNC=0._r8
  allocate(ESNC(jsken,npelms,0:1,0:JZ,JP,JY,JX));ESNC=0._r8
  allocate(ZNPP(JP,JY,JX));     ZNPP=0._r8
  allocate(CTRAN(JP,JY,JX));    CTRAN=0._r8
  allocate(TCO2A(JP,JY,JX));    TCO2A=0._r8
  allocate(HVSTE(npelms,JP,JY,JX));    HVSTE=0._r8
  allocate(THVSTE(npelms,JP,JY,JX));   THVSTE=0._r8
  allocate(VCO2F(JP,JY,JX));    VCO2F=0._r8
  allocate(VCH4F(JP,JY,JX));    VCH4F=0._r8
  allocate(VOXYF(JP,JY,JX));    VOXYF=0._r8
  allocate(VNH3F(JP,JY,JX));    VNH3F=0._r8
  allocate(VN2OF(JP,JY,JX));    VN2OF=0._r8
  allocate(VPO4F(JP,JY,JX));    VPO4F=0._r8
  allocate(ROXYP(2,JZ,JP,JY,JX));ROXYP=0._r8
  allocate(RCOFLA(2,JZ,JP,JY,JX));RCOFLA=0._r8
  allocate(ROXFLA(2,JZ,JP,JY,JX));ROXFLA=0._r8
  allocate(RCHFLA(2,JZ,JP,JY,JX));RCHFLA=0._r8
  allocate(RN2FLA(2,JZ,JP,JY,JX));RN2FLA=0._r8
  allocate(RNHFLA(2,JZ,JP,JY,JX));RNHFLA=0._r8
  allocate(RCODFA(2,JZ,JP,JY,JX));RCODFA=0._r8
  allocate(ROXDFA(2,JZ,JP,JY,JX));ROXDFA=0._r8
  allocate(RCHDFA(2,JZ,JP,JY,JX));RCHDFA=0._r8
  allocate(RN2DFA(2,JZ,JP,JY,JX));RN2DFA=0._r8
  allocate(RNHDFA(2,JZ,JP,JY,JX));RNHDFA=0._r8
  allocate(RCO2S(2,JZ,JP,JY,JX));RCO2S=0._r8
  allocate(RUPOXS(2,JZ,JP,JY,JX));RUPOXS=0._r8
  allocate(RUPCHS(2,JZ,JP,JY,JX));RUPCHS=0._r8
  allocate(RUPN2S(2,JZ,JP,JY,JX));RUPN2S=0._r8
  allocate(RUPN3S(2,JZ,JP,JY,JX));RUPN3S=0._r8
  allocate(RUPN3B(2,JZ,JP,JY,JX));RUPN3B=0._r8
  allocate(RCO2P(2,JZ,JP,JY,JX));RCO2P=0._r8
  allocate(RUPOXP(2,JZ,JP,JY,JX));RUPOXP=0._r8
  allocate(RCO2M(2,JZ,JP,JY,JX));RCO2M=0._r8
  allocate(RCO2A(2,JZ,JP,JY,JX));RCO2A=0._r8
  allocate(UPOME(1:npelms,JP,JY,JX));    UPOME=0._r8
  allocate(UPNH4(JP,JY,JX));    UPNH4=0._r8
  allocate(UPNO3(JP,JY,JX));    UPNO3=0._r8
  allocate(UPH2P(JP,JY,JX));    UPH2P=0._r8
  allocate(UPH1P(JP,JY,JX));    UPH1P=0._r8
  allocate(UPNF(JP,JY,JX));     UPNF=0._r8
  allocate(RCO2Z(JP,JY,JX));    RCO2Z=0._r8
  allocate(ROXYZ(JP,JY,JX));    ROXYZ=0._r8
  allocate(RCH4Z(JP,JY,JX));    RCH4Z=0._r8
  allocate(RN2OZ(JP,JY,JX));    RN2OZ=0._r8
  allocate(RUONH4(2,JZ,JP,JY,JX));RUONH4=0._r8
  allocate(RUONHB(2,JZ,JP,JY,JX));RUONHB=0._r8
  allocate(RUONO3(2,JZ,JP,JY,JX));RUONO3=0._r8
  allocate(RUONOB(2,JZ,JP,JY,JX));RUONOB=0._r8
  allocate(RUOH2P(2,JZ,JP,JY,JX));RUOH2P=0._r8
  allocate(RUOH2B(2,JZ,JP,JY,JX));RUOH2B=0._r8
  allocate(RUCNH4(2,JZ,JP,JY,JX));RUCNH4=0._r8
  allocate(RUCNHB(2,JZ,JP,JY,JX));RUCNHB=0._r8
  allocate(RUCNO3(2,JZ,JP,JY,JX));RUCNO3=0._r8
  allocate(RUCNOB(2,JZ,JP,JY,JX));RUCNOB=0._r8
  allocate(RUCH2P(2,JZ,JP,JY,JX));RUCH2P=0._r8
  allocate(RUCH2B(2,JZ,JP,JY,JX));RUCH2B=0._r8
  allocate(RUOH1P(2,JZ,JP,JY,JX));RUOH1P=0._r8
  allocate(RUCH1P(2,JZ,JP,JY,JX));RUCH1P=0._r8
  allocate(RUOH1B(2,JZ,JP,JY,JX));RUOH1B=0._r8
  allocate(RUCH1B(2,JZ,JP,JY,JX));RUCH1B=0._r8
  allocate(RCO2N(2,JZ,JP,JY,JX));RCO2N=0._r8
  allocate(RUNNHP(2,JZ,JP,JY,JX));RUNNHP=0._r8
  allocate(RUNNOP(2,JZ,JP,JY,JX));RUNNOP=0._r8
  allocate(RUNNBP(2,JZ,JP,JY,JX));RUNNBP=0._r8
  allocate(RUNNXP(2,JZ,JP,JY,JX));RUNNXP=0._r8
  allocate(RNH3Z(JP,JY,JX));    RNH3Z=0._r8
  allocate(RNH3B(JC,JP,JY,JX)); RNH3B=0._r8
  allocate(WFR(2,JZ,JP,JY,JX)); WFR=0._r8
  allocate(RHGFLA(2,JZ,JP,JY,JX));RHGFLA=0._r8
  allocate(RHGDFA(2,JZ,JP,JY,JX));RHGDFA=0._r8
  allocate(RH2GZ(JP,JY,JX));    RH2GZ=0._r8
  allocate(HCUPTK(JP,JY,JX));   HCUPTK=0._r8
  allocate(HZUPTK(JP,JY,JX));   HZUPTK=0._r8
  allocate(HPUPTK(JP,JY,JX));   HPUPTK=0._r8
  allocate(TEUPTK(npelms,JP,JY,JX));   TEUPTK=0._r8
  allocate(TUPWTR(0:JZ,JY,JX)); TUPWTR=0._r8
  allocate(TUPHT(0:JZ,JY,JX));  TUPHT=0._r8
  allocate(TCOFLA(JZ,JY,JX));   TCOFLA=0._r8
  allocate(TOXFLA(JZ,JY,JX));   TOXFLA=0._r8
  allocate(TCHFLA(JZ,JY,JX));   TCHFLA=0._r8
  allocate(TN2FLA(JZ,JY,JX));   TN2FLA=0._r8
  allocate(TNHFLA(JZ,JY,JX));   TNHFLA=0._r8
  allocate(TLCO2P(JZ,JY,JX));   TLCO2P=0._r8
  allocate(TLOXYP(JZ,JY,JX));   TLOXYP=0._r8
  allocate(TLCH4P(JZ,JY,JX));   TLCH4P=0._r8
  allocate(TLN2OP(JZ,JY,JX));   TLN2OP=0._r8
  allocate(TLNH3P(JZ,JY,JX));   TLNH3P=0._r8
  allocate(TCO2S(JZ,JY,JX));    TCO2S=0._r8
  allocate(TUPOXS(JZ,JY,JX));   TUPOXS=0._r8
  allocate(TUPCHS(JZ,JY,JX));   TUPCHS=0._r8
  allocate(TUPN2S(JZ,JY,JX));   TUPN2S=0._r8
  allocate(TUPN3S(JZ,JY,JX));   TUPN3S=0._r8
  allocate(TUPNH4(JZ,JY,JX));   TUPNH4=0._r8
  allocate(TUPNO3(JZ,JY,JX));   TUPNO3=0._r8
  allocate(TUPH2P(JZ,JY,JX));   TUPH2P=0._r8
  allocate(TUPN3B(JZ,JY,JX));   TUPN3B=0._r8
  allocate(TUPNHB(JZ,JY,JX));   TUPNHB=0._r8
  allocate(TUPNOB(JZ,JY,JX));   TUPNOB=0._r8
  allocate(TUPH2B(JZ,JY,JX));   TUPH2B=0._r8
  allocate(TUPNF(JZ,JY,JX));    TUPNF=0._r8
  allocate(TDFOMC(0:jcplx1,JZ,JY,JX));TDFOMC=0._r8
  allocate(TDFOMN(0:jcplx1,JZ,JY,JX));TDFOMN=0._r8
  allocate(TDFOMP(0:jcplx1,JZ,JY,JX));TDFOMP=0._r8
  allocate(TCO2P(JZ,JY,JX));    TCO2P=0._r8
  allocate(TUPOXP(JZ,JY,JX));   TUPOXP=0._r8
  allocate(RTDNT(JZ,JY,JX));    RTDNT=0._r8
  allocate(TUPH1P(JZ,JY,JX));   TUPH1P=0._r8
  allocate(TUPH1B(JZ,JY,JX));   TUPH1B=0._r8
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
  allocate(ROQCX(0:jcplx1,0:JZ,JY,JX));ROQCX=0._r8
  allocate(ROQCY(0:jcplx1,0:JZ,JY,JX));ROQCY=0._r8
  allocate(ROQAX(0:jcplx1,0:JZ,JY,JX));ROQAX=0._r8
  allocate(ROQAY(0:jcplx1,0:JZ,JY,JX));ROQAY=0._r8
  allocate(TH2GZ(JY,JX));       TH2GZ=0._r8
  allocate(TUPHGS(JZ,JY,JX));   TUPHGS=0._r8
  allocate(THGFLA(JZ,JY,JX));   THGFLA=0._r8
  allocate(TLH2GP(JZ,JY,JX));   TLH2GP=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantRates
  use abortutils, only : destroy
  implicit none
  call destroy(TCNET)
  call destroy(RNH3C)
  call destroy(TNH3C)
  call destroy(TESN0)
  call destroy(RDFOMC)
  call destroy(RDFOMN)
  call destroy(RDFOMP)
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
  call destroy(RCZLX)
  call destroy(RCPLX)
  call destroy(RCCLX)
  call destroy(RCZSX)
  call destroy(RCPSX)
  call destroy(RCCSX)
  call destroy(CARBN)
  call destroy(TESNC)
  call destroy(TZUPFX)
  call destroy(TCO2T)
  call destroy(BALE)
  call destroy(HESNC)
  call destroy(ESNC)
  call destroy(ZNPP)
  call destroy(CTRAN)
  call destroy(TCO2A)
  call destroy(HVSTE)
  call destroy(THVSTE)
  call destroy(VCO2F)
  call destroy(VCH4F)
  call destroy(VOXYF)
  call destroy(VNH3F)
  call destroy(VN2OF)
  call destroy(VPO4F)
  call destroy(ROXYP)
  call destroy(RCOFLA)
  call destroy(ROXFLA)
  call destroy(RCHFLA)
  call destroy(RN2FLA)
  call destroy(RNHFLA)
  call destroy(RCODFA)
  call destroy(ROXDFA)
  call destroy(RCHDFA)
  call destroy(RN2DFA)
  call destroy(RNHDFA)
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
  call destroy(RCO2Z)
  call destroy(ROXYZ)
  call destroy(RCH4Z)
  call destroy(RN2OZ)
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
  call destroy(RHGFLA)
  call destroy(RHGDFA)
  call destroy(RH2GZ)
  call destroy(HCUPTK)
  call destroy(HZUPTK)
  call destroy(HPUPTK)
  call destroy(TEUPTK)
  call destroy(TUPWTR)
  call destroy(TUPHT)
  call destroy(TCOFLA)
  call destroy(TOXFLA)
  call destroy(TCHFLA)
  call destroy(TN2FLA)
  call destroy(TNHFLA)
  call destroy(TLCO2P)
  call destroy(TLOXYP)
  call destroy(TLCH4P)
  call destroy(TLN2OP)
  call destroy(TLNH3P)
  call destroy(TCO2S)
  call destroy(TUPOXS)
  call destroy(TUPCHS)
  call destroy(TUPN2S)
  call destroy(TUPN3S)
  call destroy(TUPNH4)
  call destroy(TUPNO3)
  call destroy(TUPH2P)
  call destroy(TUPN3B)
  call destroy(TUPNHB)
  call destroy(TUPNOB)
  call destroy(TUPH2B)
  call destroy(TUPNF)
  call destroy(TDFOMC)
  call destroy(TDFOMN)
  call destroy(TDFOMP)
  call destroy(TCO2P)
  call destroy(TUPOXP)
  call destroy(RTDNT)
  call destroy(TUPH1P)
  call destroy(TUPH1B)
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
  call destroy(TUPHGS)
  call destroy(THGFLA)
  call destroy(TLH2GP)
  end subroutine DestructPlantRates

end module PlantDataRateType
