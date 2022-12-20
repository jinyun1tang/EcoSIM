module PlantDataRateType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  TCNET(:,:)                         !total canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH3C(:,:,:)                       !canopy NH3 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TNH3C(:,:,:)                       !total canopy NH3 flux, [g d-2 ]
  real(r8),target,allocatable ::  TESN0(:,:,:,:)                     !total surface litterfall element, [g d-2]
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
  real(r8),target,allocatable ::  CARBN(:,:,:)                       !total gross CO2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  TESNC(:,:,:,:)                     !total plant element litterfall , [g d-2 ]
  real(r8),target,allocatable ::  TZUPFX(:,:,:)                      !total plant N2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  TCO2T(:,:,:)                       !total plant respiration, [g d-2 ]
  real(r8),target,allocatable ::  BALE(:,:,:,:)                      !plant element balance, [g d-2]
  real(r8),target,allocatable ::  HESNC(:,:,:,:)                     !plant element litterfall, [g d-2 h-1]
  real(r8),target,allocatable ::  ESNC(:,:,:,:,:,:,:)                !plant litterfall element, [g d-2 h-1]
  real(r8),target,allocatable ::  ZNPP(:,:,:)                        !total net primary productivity, [g d-2]
  real(r8),target,allocatable ::  CTRAN(:,:,:)                       !total transpiration, [m d-2]
  real(r8),target,allocatable ::  TCO2A(:,:,:)                       !total autotrophic respiration, [g d-2 ]
  real(r8),target,allocatable ::  HVSTE(:,:,:,:)                     !plant element harvest, [g d-2 ]
  real(r8),target,allocatable ::  THVSTE(:,:,:,:)                    !total plant harvest, [g d-2 ]
  real(r8),target,allocatable ::  VCO2F(:,:,:)                       !plant CO2 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VCH4F(:,:,:)                       !plant CH4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VOXYF(:,:,:)                       !plant O2 uptake from fire, [g d-2 ]
  real(r8),target,allocatable ::  VNH3F(:,:,:)                       !plant NH3 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VN2OF(:,:,:)                       !plant N2O emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  VPO4F(:,:,:)                       !plant PO4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  ROXYP(:,:,:,:,:)                   !root  O2 demand from respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_RFLA(:,:,:,:,:,:)                  !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_RDFA(:,:,:,:,:,:)                  !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
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
  real(r8),target,allocatable ::  RFGas_root(:,:,:,:)                !gas flux from root disturbance [g d-2 h-1]
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
  real(r8),target,allocatable ::  TEUPTK(:,:,:,:)                    !total net root element uptake (+ve) - exudation (-ve), [g d-2 ]
  real(r8),target,allocatable ::  TUPWTR(:,:,:)                      !total root water uptake, [m3 d-2]
  real(r8),target,allocatable ::  TUPHT(:,:,:)                       !total root heat uptake, [MJ d-2]
  real(r8),target,allocatable ::  trcg_TFLA(:,:,:,:)                 !total internal root gas flux , [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_TLP(:,:,:,:)                  !total root internal gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TCO2S(:,:,:)                       !total root-soil CO2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPOXS(:,:,:)                      !total root-soil O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPCHS(:,:,:)                      !total root-soil CH4 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPN2S(:,:,:)                      !total root-soil N2O flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPN3S(:,:,:)                      !total root-soil NH3 flux non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPNH4(:,:,:)                      !total root-soil NH4 flux non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPNO3(:,:,:)                      !total root-soil NO3 flux non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPH2P(:,:,:)                      !total root-soil PO4 flux non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPN3B(:,:,:)                      !total root-soil NH3 flux band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPNHB(:,:,:)                      !total root-soil NH4 flux band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPNOB(:,:,:)                      !total root-soil NO3 flux band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPH2B(:,:,:)                      !total root-soil PO4 flux band, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPNF(:,:,:)                       !total root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  TDFOMC(:,:,:,:)                    !total root C exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  TDFOMN(:,:,:,:)                    !total root N exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  TDFOMP(:,:,:,:)                    !total root P exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  TCO2P(:,:,:)                       !total root CO2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TUPOXP(:,:,:)                      !total root internal O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RTDNT(:,:,:)                       !total root length density, [m m-3]
  real(r8),target,allocatable ::  TUPH1P(:,:,:)                      !soil-root exch of HPO4 in non-band
  real(r8),target,allocatable ::  TUPH1B(:,:,:)                      !soil-root exch of HPO4 in band
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
  real(r8),target,allocatable ::  TUPHGS(:,:,:)                      !total root-soil H2 flux, [g d-2 h-1]
  private :: InitAllocate
  contains

!----------------------------------------------------------------------
  subroutine InitPlantRates(n_pltlitrk)
  implicit none
  integer, intent(in) :: n_pltlitrk

  call InitAllocate(n_pltlitrk)

  end subroutine InitPlantRates
!----------------------------------------------------------------------

  subroutine InitAllocate(n_pltlitrk)

  implicit none
  integer, intent(in) :: n_pltlitrk

  allocate(TCNET(JY,JX));       TCNET=0._r8
  allocate(RNH3C(JP,JY,JX));    RNH3C=0._r8
  allocate(TNH3C(JP,JY,JX));    TNH3C=0._r8
  allocate(TESN0(npelms,JP,JY,JX));    TESN0=0._r8
  allocate(RDFOME(npelms,2,1:jcplx,JZ,JP,JY,JX));RDFOME=0._r8
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
  allocate(RCELX(JBR,npelms,JP,JY,JX)); RCELX=0._r8
  allocate(RCESX(JBR,npelms,JP,JY,JX)); RCESX=0._r8
  allocate(CARBN(JP,JY,JX));    CARBN=0._r8
  allocate(TESNC(npelms,JP,JY,JX));    TESNC=0._r8
  allocate(TZUPFX(JP,JY,JX));   TZUPFX=0._r8
  allocate(TCO2T(JP,JY,JX));    TCO2T=0._r8
  allocate(BALE(npelms,JP,JY,JX));     BALE=0._r8
  allocate(HESNC(npelms,JP,JY,JX));    HESNC=0._r8
  allocate(ESNC(jsken,npelms,1:n_pltlitrk,0:JZ,JP,JY,JX));ESNC=0._r8
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
  allocate(trcg_RFLA(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_RFLA=0._r8
  allocate(trcg_RDFA(idg_beg:idg_end-1,2,JZ,JP,JY,JX));trcg_RDFA=0._r8
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
  allocate(RFGas_root(idg_beg:idg_end-1,JP,JY,JX)); RFGas_root=0._r8
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
  allocate(RNH3B(JBR,JP,JY,JX)); RNH3B=0._r8
  allocate(WFR(2,JZ,JP,JY,JX)); WFR=0._r8
  allocate(RH2GZ(JP,JY,JX));    RH2GZ=0._r8
  allocate(HEUPTK(npelms,JP,JY,JX));   HEUPTK=0._r8
  allocate(TEUPTK(npelms,JP,JY,JX));   TEUPTK=0._r8
  allocate(TUPWTR(0:JZ,JY,JX)); TUPWTR=0._r8
  allocate(TUPHT(0:JZ,JY,JX));  TUPHT=0._r8
  allocate(trcg_TFLA(idg_beg:idg_end-1,JZ,JY,JX));   trcg_TFLA=0._r8
  allocate(trcg_TLP(idg_beg:idg_end-1,JZ,JY,JX));   trcg_TLP=0._r8
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
  allocate(TDFOMC(1:jcplx,JZ,JY,JX));TDFOMC=0._r8
  allocate(TDFOMN(1:jcplx,JZ,JY,JX));TDFOMN=0._r8
  allocate(TDFOMP(1:jcplx,JZ,JY,JX));TDFOMP=0._r8
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
  allocate(ROQCX(1:jcplx,0:JZ,JY,JX));ROQCX=0._r8
  allocate(ROQCY(1:jcplx,0:JZ,JY,JX));ROQCY=0._r8
  allocate(ROQAX(1:jcplx,0:JZ,JY,JX));ROQAX=0._r8
  allocate(ROQAY(1:jcplx,0:JZ,JY,JX));ROQAY=0._r8
  allocate(TH2GZ(JY,JX));       TH2GZ=0._r8
  allocate(TUPHGS(JZ,JY,JX));   TUPHGS=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantRates
  use abortutils, only : destroy
  implicit none
  call destroy(TCNET)
  call destroy(RNH3C)
  call destroy(TNH3C)
  call destroy(TESN0)
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
  call destroy(TEUPTK)
  call destroy(TUPWTR)
  call destroy(TUPHT)
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
  end subroutine DestructPlantRates

end module PlantDataRateType
