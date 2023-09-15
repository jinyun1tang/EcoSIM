module CanopyDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use EcoSIMConfig, only : jsken => jskenc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  ALBP(:,:,:)                        !canopy PAR albedo , [-]
  real(r8),target,allocatable ::  TAUP(:,:,:)                        !canopy PAR transmissivity , [-]
  real(r8),target,allocatable ::  ABSR(:,:,:)                        !canopy shortwave absorptivity , [-]
  real(r8),target,allocatable ::  ABSP(:,:,:)                        !canopy PAR absorptivity , [-]
  real(r8),target,allocatable ::  RSMX(:,:,:)                        !maximum stomatal resistance to vapor, [s m-1]
  real(r8),target,allocatable ::  RCMX(:,:,:)                        !maximum stomatal resistance to CO2, [s h-1]
  real(r8),target,allocatable ::  MaxCanPStomaResistH2O(:,:,:)      !maximum stomatal resistance to vapor, [s h-1]
  real(r8),target,allocatable ::  RCS(:,:,:)                         !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8),target,allocatable ::  CanPStomaResistH2O(:,:,:)         !canopy stomatal resistance, [h m-1]
  real(r8),target,allocatable ::  MinCanPStomaResistH2O(:,:,:)      !canopy minimum stomatal resistance, [s m-1]
  real(r8),target,allocatable ::  RAC(:,:)                           !canopy boundary layer resistance, [m h-1]
  real(r8),target,allocatable ::  O2I(:,:,:)                         !leaf gaseous O2 concentration, [umol m-3]
  real(r8),target,allocatable ::  CO2I(:,:,:)                        !leaf gaseous CO2 concentration, [umol m-3]
  real(r8),target,allocatable ::  FMOL(:,:,:)                        !total gas concentration, [mol m-3]
  real(r8),target,allocatable ::  DCO2(:,:,:)                        !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8),target,allocatable ::  CO2Q(:,:,:)                        !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8),target,allocatable ::  CO2L(:,:,:)                        !leaf aqueous CO2 concentration, [uM]
  real(r8),target,allocatable ::  O2L(:,:,:)                         !leaf aqueous O2 concentration, [uM]
  real(r8),target,allocatable ::  SCO2(:,:,:)                        !leaf CO2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  SO2(:,:,:)                         !leaf O2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  XKCO2L(:,:,:)                      !leaf aqueous CO2 Km no O2, [uM]
  real(r8),target,allocatable ::  XKCO2O(:,:,:)                      !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8),target,allocatable ::  CHILL(:,:,:)                       !chilling effect on CO2 fixation, [-]
  real(r8),target,allocatable ::  VCGRO(:,:,:,:,:)                   !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  VGRO(:,:,:,:,:)                    !carboxylation rate, [umol m-2 s-1]
  real(r8),target,allocatable ::  COMPL(:,:,:,:,:)                   !CO2 compensation point, [uM]
  real(r8),target,allocatable ::  ETGRO(:,:,:,:,:)                   !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  CBXN(:,:,:,:,:)                    !carboxylation efficiency, [umol umol-1]
  real(r8),target,allocatable ::  CO2B(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  VCGR4(:,:,:,:,:)                   !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  VGRO4(:,:,:,:,:)                   !C4 carboxylation rate, [umol m-2 s-1]
  real(r8),target,allocatable ::  ETGR4(:,:,:,:,:)                   !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  CBXN4(:,:,:,:,:)                   !C4 carboxylation efficiency, [umol umol-1]
  real(r8),target,allocatable ::  CPOOL4(:,:,:,:,:)                  !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  HCOB(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  FDBK(:,:,:,:)                      !branch down-regulation of CO2 fixation, [-]
  real(r8),target,allocatable ::  FDBK4(:,:,:,:,:)                   !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  FDBKX(:,:,:,:)                     !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  CNETX(:,:)                         !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  VCMX(:,:,:)                        !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VOMX(:,:,:)                        !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VCMX4(:,:,:)                       !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  XKCO2(:,:,:)                       !Km for rubisco carboxylase activity, [uM]
  real(r8),target,allocatable ::  XKO2(:,:,:)                        !Km for rubisco oxygenase activity, [uM]
  real(r8),target,allocatable ::  XKCO24(:,:,:)                      !Km for PEP carboxylase activity, [uM]
  real(r8),target,allocatable ::  RUBP(:,:,:)                        !leaf rubisco content, [g g-1]
  real(r8),target,allocatable ::  PEPC(:,:,:)                        !leaf PEP carboxylase content, [g g-1]
  real(r8),target,allocatable ::  ETMX(:,:,:)                        !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  CHL(:,:,:)                         !leaf C3 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  CHL4(:,:,:)                        !leaf C4 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  FCO2(:,:,:)                        !Ci:Ca ratio, [-]
  real(r8),target,allocatable ::  RadNet2CanP(:,:,:)                 !canopy net radiation , [MJ d-2 h-1] >0
  real(r8),target,allocatable ::  LWRadCanP(:,:,:)                   !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  SWRadByCanP(:,:,:)                 !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  PARByCanP(:,:,:)                   !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),target,allocatable ::  FracPARByCanP(:,:,:)                       !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),target,allocatable ::  TAU0(:,:,:)                        !fraction of radiation transmitted by canopy layer, [-]
  real(r8),target,allocatable ::  TAUS(:,:,:)                        !fraction of radiation intercepted by canopy layer, [-]
  real(r8),target,allocatable ::  FRADG(:,:)                         !fraction of radiation intercepted by ground surface, [-]
  real(r8),target,allocatable ::  RADG(:,:)                          !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  LWRadCanGPrev(:,:)                        !longwave radiation emitted by canopy, [MJ m-2 h-1]
  real(r8),target,allocatable ::  LWRadGrnd(:,:)                        !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  CanH2OHeldVg(:,:)                  !grid canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  TFLWCI(:,:)                        !net ice transfer to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  PrecIntcptByCanG(:,:)              !grid net precipitation water interception to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  EvapTransHeatP(:,:,:)                       !canopy latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatXAir2PCan(:,:,:)               !air to canopy sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatStorCanP(:,:,:)                       !canopy storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ENGYX(:,:,:)                       !canopy heat storage from previous time step, [MJ d-2]
  real(r8),target,allocatable ::  VHeatCapCanP(:,:,:)                       !canopy heat capacity, [MJ d-2 K-1]
  real(r8),target,allocatable ::  PSICanP(:,:,:)                     !plant canopy total water potential , [Mpa]
  real(r8),target,allocatable ::  PSICanPTurg(:,:,:)                       !plant canopy turgor water potential, [Mpa]
  real(r8),target,allocatable ::  PSICanPOsmo(:,:,:)                 !platn canopy osmotic water potential, [Mpa]
  real(r8),target,allocatable ::  CanPbndlResist(:,:,:)                          !canopy boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  PTrans(:,:,:)                          !canopy transpiration, [m2 d-2 h-1]
  real(r8),target,allocatable ::  VapXAir2PCan(:,:,:)                !negative of canopy evaporation, [m2 d-2 h-1]
  real(r8),target,allocatable ::  CanWatP(:,:,:)                       !canopy water content associated with dry matter, [m3 d-2]
  real(r8),target,allocatable ::  TEVAPP(:,:)                        !total canopy evaporation + transpiration, [m3 d-2]
  real(r8),target,allocatable ::  VapXAir2CanG(:,:)                        !total canopy evaporation, [m3 d-2]
  real(r8),target,allocatable ::  TENGYC(:,:)                        !total canopy heat content, [MJ  d-2]
  real(r8),target,allocatable ::  THFLXC(:,:)                        !total canopy heat flux, [MJ  d-2]
  real(r8),target,allocatable ::  CanWatg(:,:)                       !total canopy water content stored in dry matter, [m3 d-2]
  real(r8),target,allocatable ::  LWRadCanG(:,:)                         !total canopy LW emission, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ALBR(:,:,:)                        !canopy shortwave albedo , [-]
  real(r8),target,allocatable ::  TAUR(:,:,:)                        !canopy shortwave transmissivity , [-]
  real(r8),target,allocatable ::  PrecIntcptByCanP(:,:,:)                        !water flux into plant canopy, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WatByPCan(:,:,:)                   !canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  TKC(:,:,:)                         !canopy temperature, [K]
  real(r8),target,allocatable ::  TCC(:,:,:)                         !canopy temperature, [oC]
  real(r8),target,allocatable ::  DTKC(:,:,:)                        !change in canopy temperature, [K]
  real(r8),target,allocatable ::  TKCZ(:,:,:)                        !canopy temperature, [K]
  real(r8),target,allocatable ::  CPOOL3(:,:,:,:,:)                  !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8),target,allocatable ::  RSETE(:,:,:,:)                     !effect of canopy element status on seed set , []
  real(r8),target,allocatable ::  WGLFT(:,:,:)                       !total leaf mass, [g d-2]
  real(r8),target,allocatable ::  CFOPE(:,:,:,:,:,:)                 !litter kinetic fraction, [-]
  real(r8),target,allocatable ::  CanPShootElmMass(:,:,:,:)                    !canopy shoot element, [g d-2]
  real(r8),target,allocatable ::  WTLFE(:,:,:,:)                     !canopy leaf element, [g d-2]
  real(r8),target,allocatable ::  WTSHEE(:,:,:,:)                    !canopy sheath element , [g d-2]
  real(r8),target,allocatable ::  WTSTKE(:,:,:,:)                    !canopy stalk element, [g d-2]
  real(r8),target,allocatable ::  CanPStalkC(:,:,:)                       !canopy active stalk C, [g d-2]
  real(r8),target,allocatable ::  WTRSVE(:,:,:,:)                    !canopy reserve element, [g d-2]
  real(r8),target,allocatable ::  WTHSKE(:,:,:,:)                    !canopy husk element, [g d-2]
  real(r8),target,allocatable ::  WTEARE(:,:,:,:)                    !canopy ear element, [g d-2]
  real(r8),target,allocatable ::  WTGRE(:,:,:,:)                     !canopy grain element, [g d-2]
  real(r8),target,allocatable ::  CanPLeafShethC(:,:,:)              !plant canopy leaf + sheath C, [gC d-2]
  real(r8),target,allocatable ::  ARLFV(:,:,:,:)                     !canopy layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CNET(:,:,:)                        !canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  WGLFV(:,:,:,:)                     !canopy layer leaf C, [g d-2]
  real(r8),target,allocatable ::  EPOOLP(:,:,:,:)                    !canopy nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CEPOLP(:,:,:,:)                    !canopy nonstructural element concentration, [g d-2]
  real(r8),target,allocatable ::  CanPLSA(:,:,:,:)                   !plant canopy layer stem area, [m2 d-2]
  real(r8),target,allocatable ::  EPOLNP(:,:,:,:)                    !canopy nodule nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CanPBStalkC(:,:,:,:)                    !branch active stalk C, [g d-2]
  real(r8),target,allocatable ::  EPOOL(:,:,:,:,:)                   !branch nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CanPBLeafShethC(:,:,:,:)           !plant branch leaf + sheath C, [g d-2]
  real(r8),target,allocatable ::  WTSHTBE(:,:,:,:,:)                 !branch shoot C, [g d-2]
  real(r8),target,allocatable ::  WTLFBE(:,:,:,:,:)                  !branch leaf element, [g d-2]
  real(r8),target,allocatable ::  WTSHEBE(:,:,:,:,:)                 !branch sheath element , [g d-2]
  real(r8),target,allocatable ::  WTSTKBE(:,:,:,:,:)                  !branch stalk element, [g d-2]
  real(r8),target,allocatable ::  WTRSVBE(:,:,:,:,:)                  !branch reserve element, [g d-2]
  real(r8),target,allocatable ::  WTHSKBE(:,:,:,:,:)                  !branch husk element, [g d-2]
  real(r8),target,allocatable ::  WTEARBE(:,:,:,:,:)                 !branch ear element, [g d-2]
  real(r8),target,allocatable ::  WTGRBE(:,:,:,:,:)                  !branch grain element, [g d-2]
  real(r8),target,allocatable ::  CEPOLB(:,:,:,:,:)                    !branch nonstructural C concentration, [g d-2]
  real(r8),target,allocatable ::  EPOLNB(:,:,:,:,:)                  !branch nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  WTNDBE(:,:,:,:,:)                  !branch nodule element, [g d-2]
  real(r8),target,allocatable ::  WGSHEXE(:,:,:,:,:)                  !branch sheath structural element, [g d-2]
  real(r8),target,allocatable ::  WTSTXBE(:,:,:,:,:)                    !branch stalk structural C, [g d-2]
  real(r8),target,allocatable ::  WGLFEX(:,:,:,:,:)                     !branch leaf structural element, [g d-2]
  real(r8),target,allocatable ::  WGLFE(:,:,:,:,:,:)                    !leaf element, [g d-2]
  real(r8),target,allocatable ::  WGSHE(:,:,:,:,:,:)                 !sheath element , [g d-2]
  real(r8),target,allocatable ::  WGNODE(:,:,:,:,:,:)                  !internode element, [g d-2]
  real(r8),target,allocatable ::  WGLFLE(:,:,:,:,:,:,:)                 !layer leaf element, [g d-2]
  real(r8),target,allocatable ::  CanPLNBLA(:,:,:,:,:,:)                 !layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  WSLF(:,:,:,:,:)                    !layer leaf protein C, [g d-2]
  real(r8),target,allocatable ::  WSSHE(:,:,:,:,:)                   !layer sheath protein C, [g d-2]
  real(r8),target,allocatable ::  CCPLNP(:,:,:)                      !nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  GRWTB(:,:,:,:)                     !maximum grain C during grain fill, [g d-2]
  real(r8),target,allocatable ::  WTSTDE(:,:,:,:,:)                  !standing dead element fraction, [g d-2]
  real(r8),target,allocatable ::  WTSTGE(:,:,:,:)                    !standing dead element, [g d-2]
  real(r8),target,allocatable ::  WTRVE(:,:,:,:)                     !plant stored nonstructural element, [g d-2]
  real(r8),target,allocatable ::  WTRVX(:,:,:)                       !plant stored nonstructural C at planting, [g d-2]
  REAL(R8),target,allocatable ::  WTSHTA(:,:,:)                      !landscape average canopy shoot C, [g d-2]
  contains
!----------------------------------------------------------------------

  subroutine InitCanopyData

  implicit none
  allocate(ALBP(JP,JY,JX));     ALBP=0._r8
  allocate(TAUP(JP,JY,JX));     TAUP=0._r8
  allocate(ABSR(JP,JY,JX));     ABSR=0._r8
  allocate(ABSP(JP,JY,JX));     ABSP=0._r8
  allocate(RSMX(JP,JY,JX));     RSMX=0._r8
  allocate(RCMX(JP,JY,JX));     RCMX=0._r8
  allocate(MaxCanPStomaResistH2O(JP,JY,JX));     MaxCanPStomaResistH2O=0._r8
  allocate(RCS(JP,JY,JX));      RCS=0._r8
  allocate(CanPStomaResistH2O(JP,JY,JX));       CanPStomaResistH2O=0._r8
  allocate(MinCanPStomaResistH2O(JP,JY,JX));     MinCanPStomaResistH2O=0._r8
  allocate(RAC(JY,JX));         RAC=0._r8
  allocate(O2I(JP,JY,JX));      O2I=0._r8
  allocate(CO2I(JP,JY,JX));     CO2I=0._r8
  allocate(FMOL(JP,JY,JX));     FMOL=0._r8
  allocate(DCO2(JP,JY,JX));     DCO2=0._r8
  allocate(CO2Q(JP,JY,JX));     CO2Q=0._r8
  allocate(CO2L(JP,JY,JX));     CO2L=0._r8
  allocate(O2L(JP,JY,JX));      O2L=0._r8
  allocate(SCO2(JP,JY,JX));     SCO2=0._r8
  allocate(SO2(JP,JY,JX));      SO2=0._r8
  allocate(XKCO2L(JP,JY,JX));   XKCO2L=0._r8
  allocate(XKCO2O(JP,JY,JX));   XKCO2O=0._r8
  allocate(CHILL(JP,JY,JX));    CHILL=0._r8
  allocate(VCGRO(JNODS,JBR,JP,JY,JX));VCGRO=0._r8
  allocate(VGRO(JNODS,JBR,JP,JY,JX));VGRO=0._r8
  allocate(COMPL(JNODS,JBR,JP,JY,JX));COMPL=0._r8
  allocate(ETGRO(JNODS,JBR,JP,JY,JX));ETGRO=0._r8
  allocate(CBXN(JNODS,JBR,JP,JY,JX));CBXN=0._r8
  allocate(CO2B(JNODS,JBR,JP,JY,JX));CO2B=0._r8
  allocate(VCGR4(JNODS,JBR,JP,JY,JX));VCGR4=0._r8
  allocate(VGRO4(JNODS,JBR,JP,JY,JX));VGRO4=0._r8
  allocate(ETGR4(JNODS,JBR,JP,JY,JX));ETGR4=0._r8
  allocate(CBXN4(JNODS,JBR,JP,JY,JX));CBXN4=0._r8
  allocate(CPOOL4(JNODS,JBR,JP,JY,JX));CPOOL4=0._r8
  allocate(HCOB(JNODS,JBR,JP,JY,JX));HCOB=0._r8
  allocate(FDBK(JBR,JP,JY,JX));  FDBK=0._r8
  allocate(FDBK4(JNODS,JBR,JP,JY,JX));FDBK4=0._r8
  allocate(FDBKX(JBR,JP,JY,JX)); FDBKX=0._r8
  allocate(CNETX(JY,JX));       CNETX=0._r8
  allocate(VCMX(JP,JY,JX));     VCMX=0._r8
  allocate(VOMX(JP,JY,JX));     VOMX=0._r8
  allocate(VCMX4(JP,JY,JX));    VCMX4=0._r8
  allocate(XKCO2(JP,JY,JX));    XKCO2=0._r8
  allocate(XKO2(JP,JY,JX));     XKO2=0._r8
  allocate(XKCO24(JP,JY,JX));   XKCO24=0._r8
  allocate(RUBP(JP,JY,JX));     RUBP=0._r8
  allocate(PEPC(JP,JY,JX));     PEPC=0._r8
  allocate(ETMX(JP,JY,JX));     ETMX=0._r8
  allocate(CHL(JP,JY,JX));      CHL=0._r8
  allocate(CHL4(JP,JY,JX));     CHL4=0._r8
  allocate(FCO2(JP,JY,JX));     FCO2=0._r8
  allocate(RadNet2CanP(JP,JY,JX));     RadNet2CanP=0._r8
  allocate(LWRadCanP(JP,JY,JX));    LWRadCanP=0._r8
  allocate(SWRadByCanP(JP,JY,JX));     SWRadByCanP=0._r8
  allocate(PARByCanP(JP,JY,JX));     PARByCanP=0._r8
  allocate(FracPARByCanP(JP,JY,JX));    FracPARByCanP=0._r8
  allocate(TAU0(JC+1,JY,JX));   TAU0=0._r8
  allocate(TAUS(JC+1,JY,JX));   TAUS=0._r8
  allocate(FRADG(JY,JX));       FRADG=0._r8
  allocate(RADG(JY,JX));        RADG=0._r8
  allocate(LWRadCanGPrev(JY,JX));      LWRadCanGPrev=0._r8
  allocate(LWRadGrnd(JY,JX));      LWRadGrnd=0._r8
  allocate(CanH2OHeldVg(JY,JX));      CanH2OHeldVg=0._r8
  allocate(TFLWCI(JY,JX));      TFLWCI=0._r8
  allocate(PrecIntcptByCanG(JY,JX));       PrecIntcptByCanG=0._r8
  allocate(EvapTransHeatP(JP,JY,JX));    EvapTransHeatP=0._r8
  allocate(HeatXAir2PCan(JP,JY,JX));    HeatXAir2PCan=0._r8
  allocate(HeatStorCanP(JP,JY,JX));    HeatStorCanP=0._r8
  allocate(ENGYX(JP,JY,JX));    ENGYX=0._r8
  allocate(VHeatCapCanP(JP,JY,JX));    VHeatCapCanP=0._r8
  allocate(PSICanP(JP,JY,JX));    PSICanP=0._r8
  allocate(PSICanPTurg(JP,JY,JX));    PSICanPTurg=0._r8
  allocate(PSICanPOsmo(JP,JY,JX));    PSICanPOsmo=0._r8
  allocate(CanPbndlResist(JP,JY,JX));       CanPbndlResist=0._r8
  allocate(PTrans(JP,JY,JX));       PTrans=0._r8
  allocate(VapXAir2PCan(JP,JY,JX));    VapXAir2PCan=0._r8
  allocate(CanWatP(JP,JY,JX));    CanWatP=0._r8
  allocate(TEVAPP(JY,JX));      TEVAPP=0._r8
  allocate(VapXAir2CanG(JY,JX));      VapXAir2CanG=0._r8
  allocate(TENGYC(JY,JX));      TENGYC=0._r8
  allocate(THFLXC(JY,JX));      THFLXC=0._r8
  allocate(CanWatg(JY,JX));      CanWatg=0._r8
  allocate(LWRadCanG(JY,JX));       LWRadCanG=0._r8
  allocate(ALBR(JP,JY,JX));     ALBR=0._r8
  allocate(TAUR(JP,JY,JX));     TAUR=0._r8
  allocate(PrecIntcptByCanP(JP,JY,JX));     PrecIntcptByCanP=0._r8
  allocate(WatByPCan(JP,JY,JX));    WatByPCan=0._r8
  allocate(TKC(JP,JY,JX));      TKC=0._r8
  allocate(TCC(JP,JY,JX));      TCC=0._r8
  allocate(DTKC(JP,JY,JX));     DTKC=0._r8
  allocate(TKCZ(JP,JY,JX));     TKCZ=0._r8
  allocate(CPOOL3(JNODS,JBR,JP,JY,JX));CPOOL3=0._r8
  allocate(RSETE(npelms,JP,JY,JX));    RSETE=0._r8
  allocate(WGLFT(JC,JY,JX));    WGLFT=0._r8
  allocate(CFOPE(npelms,0:Jlitgrp,jsken,JP,JY,JX));CFOPE=0._r8
  allocate(CanPShootElmMass(npelms,JP,JY,JX)); CanPShootElmMass=0._r8
  allocate(WTLFE(npelms,JP,JY,JX));  WTLFE=0._r8
  allocate(WTSHEE(npelms,JP,JY,JX)); WTSHEE=0._r8
  allocate(WTSTKE(npelms,JP,JY,JX)); WTSTKE=0._r8
  allocate(CanPStalkC(JP,JY,JX));    CanPStalkC=0._r8
  allocate(WTRSVE(npelms,JP,JY,JX));    WTRSVE=0._r8
  allocate(WTHSKE(npelms,JP,JY,JX));    WTHSKE=0._r8
  allocate(WTEARE(npelms,JP,JY,JX));    WTEARE=0._r8
  allocate(WTGRE(npelms,JP,JY,JX));     WTGRE=0._r8
  allocate(CanPLeafShethC(JP,JY,JX));     CanPLeafShethC=0._r8
  allocate(ARLFV(JC,JP,JY,JX)); ARLFV=0._r8
  allocate(CNET(JP,JY,JX));     CNET=0._r8
  allocate(WGLFV(JC,JP,JY,JX)); WGLFV=0._r8
  allocate(EPOOLP(npelms,JP,JY,JX));   EPOOLP=0._r8
  allocate(CEPOLP(npelms,JP,JY,JX));   CEPOLP=0._r8
  allocate(CanPLSA(JC,JP,JY,JX)); CanPLSA=0._r8
  allocate(EPOLNP(npelms,JP,JY,JX));   EPOLNP=0._r8
  allocate(CanPBStalkC(JBR,JP,JY,JX));CanPBStalkC=0._r8
  allocate(EPOOL(npelms,JBR,JP,JY,JX)); EPOOL=0._r8
  allocate(CanPBLeafShethC(JBR,JP,JY,JX)); CanPBLeafShethC=0._r8
  allocate(WTSHTBE(npelms,JBR,JP,JY,JX));WTSHTBE=0._r8
  allocate(WTLFBE(npelms,JBR,JP,JY,JX)); WTLFBE=0._r8
  allocate(WTSHEBE(npelms,JBR,JP,JY,JX));WTSHEBE=0._r8
  allocate(WTSTKBE(npelms,JBR,JP,JY,JX));WTSTKBE=0._r8
  allocate(WTRSVBE(npelms,JBR,JP,JY,JX));WTRSVBE=0._r8
  allocate(WTHSKBE(npelms,JBR,JP,JY,JX));WTHSKBE=0._r8
  allocate(WTEARBE(npelms,JBR,JP,JY,JX));WTEARBE=0._r8
  allocate(WTGRBE(npelms,JBR,JP,JY,JX)); WTGRBE=0._r8
  allocate(CEPOLB(npelms,JBR,JP,JY,JX));CEPOLB=0._r8
  allocate(EPOLNB(npelms,JBR,JP,JY,JX));EPOLNB=0._r8
  allocate(WTNDBE(npelms,JBR,JP,JY,JX)); WTNDBE=0._r8
  allocate(WGSHEXE(npelms,JBR,JP,JY,JX));WGSHEXE=0._r8
  allocate(WTSTXBE(npelms,JBR,JP,JY,JX));WTSTXBE=0._r8
  allocate(WGLFEX(npelms,JBR,JP,JY,JX)); WGLFEX=0._r8
  allocate(WGLFE(npelms,0:JNODS,JBR,JP,JY,JX));WGLFE=0._r8
  allocate(WGSHE(npelms,0:JNODS,JBR,JP,JY,JX));WGSHE=0._r8
  allocate(WGNODE(npelms,0:JNODS,JBR,JP,JY,JX));WGNODE=0._r8
  allocate(WGLFLE(npelms,JC,0:JNODS,JBR,JP,JY,JX));WGLFLE=0._r8
  allocate(CanPLNBLA(JC,0:JNODS,JBR,JP,JY,JX));CanPLNBLA=0._r8
  allocate(WSLF(0:JNODS,JBR,JP,JY,JX));WSLF=0._r8
  allocate(WSSHE(0:JNODS,JBR,JP,JY,JX));WSSHE=0._r8
  allocate(CCPLNP(JP,JY,JX));   CCPLNP=0._r8
  allocate(GRWTB(JBR,JP,JY,JX)); GRWTB=0._r8
  allocate(WTSTDE(npelms,jsken,JP,JY,JX)); WTSTDE=0._r8
  allocate(WTSTGE(npelms,JP,JY,JX));    WTSTGE=0._r8
  allocate(WTRVE(npelms,JP,JY,JX));  WTRVE=0._r8
  allocate(WTRVX(JP,JY,JX));    WTRVX=0._r8
  allocate(WTSHTA(JP,JY,JX));   WTSHTA=0._r8
  end subroutine InitCanopyData

!----------------------------------------------------------------------
  subroutine DestructCanopyData
  use abortutils, only : destroy
  implicit none

  call destroy(ALBP)
  call destroy(TAUP)
  call destroy(ABSR)
  call destroy(ABSP)
  call destroy(RSMX)
  call destroy(RCMX)
  call destroy(MaxCanPStomaResistH2O)
  call destroy(RCS)
  call destroy(CanPStomaResistH2O)
  call destroy(MinCanPStomaResistH2O)
  call destroy(RAC)
  call destroy(O2I)
  call destroy(CO2I)
  call destroy(FMOL)
  call destroy(DCO2)
  call destroy(CO2Q)
  call destroy(CO2L)
  call destroy(O2L)
  call destroy(SCO2)
  call destroy(SO2)
  call destroy(XKCO2L)
  call destroy(XKCO2O)
  call destroy(CHILL)
  call destroy(VCGRO)
  call destroy(VGRO)
  call destroy(COMPL)
  call destroy(ETGRO)
  call destroy(CBXN)
  call destroy(CO2B)
  call destroy(VCGR4)
  call destroy(VGRO4)
  call destroy(ETGR4)
  call destroy(CBXN4)
  call destroy(CPOOL4)
  call destroy(HCOB)
  call destroy(FDBK)
  call destroy(FDBK4)
  call destroy(FDBKX)
  call destroy(CNETX)
  call destroy(VCMX)
  call destroy(VOMX)
  call destroy(VCMX4)
  call destroy(XKCO2)
  call destroy(XKO2)
  call destroy(XKCO24)
  call destroy(RUBP)
  call destroy(PEPC)
  call destroy(ETMX)
  call destroy(CHL)
  call destroy(CHL4)
  call destroy(FCO2)
  call destroy(RadNet2CanP)
  call destroy(LWRadCanP)
  call destroy(SWRadByCanP)
  call destroy(PARByCanP)
  call destroy(FracPARByCanP)
  call destroy(TAU0)
  call destroy(TAUS)
  call destroy(FRADG)
  call destroy(RADG)
  call destroy(LWRadCanGPrev)
  call destroy(LWRadGrnd)
  call destroy(CanH2OHeldVg)
  call destroy(TFLWCI)
  call destroy(PrecIntcptByCanG)
  call destroy(EvapTransHeatP)
  call destroy(HeatXAir2PCan)
  call destroy(HeatStorCanP)
  call destroy(ENGYX)
  call destroy(VHeatCapCanP)
  call destroy(PSICanP)
  call destroy(PSICanPTurg)
  call destroy(PSICanPOsmo)
  call destroy(CanPbndlResist)
  call destroy(PTrans)
  call destroy(VapXAir2PCan)
  call destroy(CanWatP)
  call destroy(TEVAPP)
  call destroy(VapXAir2CanG)
  call destroy(TENGYC)
  call destroy(THFLXC)
  call destroy(CanWatg)
  call destroy(LWRadCanG)
  call destroy(ALBR)
  call destroy(TAUR)
  call destroy(PrecIntcptByCanP)
  call destroy(WatByPCan)
  call destroy(TKC)
  call destroy(TCC)
  call destroy(DTKC)
  call destroy(TKCZ)
  call destroy(CPOOL3)
  call destroy(RSETE)
  call destroy(WGLFT)
  call destroy(CFOPE)
  call destroy(CanPShootElmMass)
  call destroy(WTLFE)
  call destroy(WTSHEE)
  call destroy(WTSTKE)
  call destroy(CanPStalkC)
  call destroy(WTRSVE)
  call destroy(WTHSKE)
  call destroy(WTEARE)
  call destroy(WTGRE)
  call destroy(CanPLeafShethC)
  call destroy(ARLFV)
  call destroy(CNET)
  call destroy(WGLFV)
  call destroy(EPOOLP)
  call destroy(CEPOLP)
  call destroy(CanPLSA)
  call destroy(EPOLNP)
  call destroy(CanPBStalkC)
  call destroy(EPOOL)
  call destroy(CanPBLeafShethC)
  call destroy(WTSHTBE)
  call destroy(WTLFBE)
  call destroy(WTSHEBE)
  call destroy(WTSTKBE)
  call destroy(WTRSVBE)
  call destroy(WTHSKBE)
  call destroy(WTEARBE)
  call destroy(WTGRBE)
  call destroy(CEPOLB)
  call destroy(EPOLNB)
  call destroy(WTNDBE)
  call destroy(WGSHEXE)
  call destroy(WTSTXBE)
  call destroy(WGLFEX)
  call destroy(WGLFE)
  call destroy(WGSHE)
  call destroy(WGNODE)
  call destroy(WGLFLE)
  call destroy(CanPLNBLA)
  call destroy(WSLF)
  call destroy(WSSHE)
  call destroy(CCPLNP)
  call destroy(GRWTB)
  call destroy(WTSTDE)
  call destroy(WTSTGE)
  call destroy(WTRVE)
  call destroy(WTRVX)
  call destroy(WTSHTA)
  end subroutine DestructCanopyData

end module CanopyDataType
