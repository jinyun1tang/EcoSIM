module CanopyDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
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
  real(r8),target,allocatable ::  RSMH(:,:,:)                        !maximum stomatal resistance to vapor, [s h-1]
  real(r8),target,allocatable ::  RCS(:,:,:)                         !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8),target,allocatable ::  RC(:,:,:)                          !canopy stomatal resistance, [h m-1]
  real(r8),target,allocatable ::  RSMN(:,:,:)                        !canopy minimum stomatal resistance, [s m-1]
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
  real(r8),target,allocatable ::  RAD1(:,:,:)                        !canopy net radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  THRM1(:,:,:)                       !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RADC(:,:,:)                        !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RADP(:,:,:)                        !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),target,allocatable ::  FRADP(:,:,:)                       !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),target,allocatable ::  TAU0(:,:,:)                        !fraction of radiation transmitted by canopy layer, [-]
  real(r8),target,allocatable ::  TAUS(:,:,:)                        !fraction of radiation intercepted by canopy layer, [-]
  real(r8),target,allocatable ::  FRADG(:,:)                         !fraction of radiation intercepted by ground surface, [-]
  real(r8),target,allocatable ::  RADG(:,:)                          !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  THRMCX(:,:)                        !longwave radiation emitted by canopy, [MJ m-2 h-1]
  real(r8),target,allocatable ::  THRMGX(:,:)                        !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  TVOLWC(:,:)                        !canopy surface water content, [m3 d-2]
  real(r8),target,allocatable ::  TFLWCI(:,:)                        !net ice transfer to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  TFLWC(:,:)                         !net water transfer to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  EFLXC(:,:,:)                       !canopy latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  SFLXC(:,:,:)                       !canopy sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HFLXC(:,:,:)                       !canopy storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ENGYX(:,:,:)                       !canopy heat storage from previous time step, [MJ d-2]
  real(r8),target,allocatable ::  VHCPC(:,:,:)                       !canopy heat capacity, [MJ d-2 K-1]
  real(r8),target,allocatable ::  PSILT(:,:,:)                       !canopy total water potential , [Mpa]
  real(r8),target,allocatable ::  PSILG(:,:,:)                       !canopy turgor water potential, [Mpa]
  real(r8),target,allocatable ::  PSILO(:,:,:)                       !canopy osmotic water potential, [Mpa]
  real(r8),target,allocatable ::  RA(:,:,:)                          !canopy boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  EP(:,:,:)                          !canopy transpiration, [m2 d-2 h-1]
  real(r8),target,allocatable ::  EVAPC(:,:,:)                       !canopy evaporation, [m2 d-2 h-1]
  real(r8),target,allocatable ::  VOLWP(:,:,:)                       !canopy water content, [m3 d-2]
  real(r8),target,allocatable ::  TEVAPP(:,:)                        !total canopy evaporation + transpiration, [m3 d-2]
  real(r8),target,allocatable ::  TEVAPC(:,:)                        !total canopy evaporation, [m3 d-2]
  real(r8),target,allocatable ::  TENGYC(:,:)                        !total canopy heat content, [MJ  d-2]
  real(r8),target,allocatable ::  THFLXC(:,:)                        !total canopy heat flux, [MJ  d-2]
  real(r8),target,allocatable ::  TVOLWP(:,:)                        !total canopy water content, [m3 d-2]
  real(r8),target,allocatable ::  THRMC(:,:)                         !total canopy LW emission, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ALBR(:,:,:)                        !canopy shortwave albedo , [-]
  real(r8),target,allocatable ::  TAUR(:,:,:)                        !canopy shortwave transmissivity , [-]
  real(r8),target,allocatable ::  FLWC(:,:,:)                        !water flux into canopy, [m3 d-2 h-1]
  real(r8),target,allocatable ::  VOLWC(:,:,:)                       !canopy surface water content, [m3 d-2]
  real(r8),target,allocatable ::  TKC(:,:,:)                         !canopy temperature, [K]
  real(r8),target,allocatable ::  TCC(:,:,:)                         !canopy temperature, [oC]
  real(r8),target,allocatable ::  DTKC(:,:,:)                        !change in canopy temperature, [K]
  real(r8),target,allocatable ::  TKCZ(:,:,:)                        !canopy temperature, [K]
  real(r8),target,allocatable ::  CPOOL3(:,:,:,:,:)                  !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8),target,allocatable ::  RSETE(:,:,:,:)                     !effect of canopy element status on seed set , []
  real(r8),target,allocatable ::  WGLFT(:,:,:)                       !total leaf mass, [g d-2]
  real(r8),target,allocatable ::  CFOPE(:,:,:,:,:,:)                 !litter kinetic fraction, [-]
  real(r8),target,allocatable ::  WTSHTE(:,:,:,:)                    !canopy shoot element, [g d-2]
  real(r8),target,allocatable ::  WTLFE(:,:,:,:)                     !canopy leaf element, [g d-2]
  real(r8),target,allocatable ::  WTSHEE(:,:,:,:)                    !canopy sheath element , [g d-2]
  real(r8),target,allocatable ::  WTSTKE(:,:,:,:)                    !canopy stalk element, [g d-2]
  real(r8),target,allocatable ::  WVSTK(:,:,:)                       !canopy active stalk C, [g d-2]
  real(r8),target,allocatable ::  WTRSVE(:,:,:,:)                    !canopy reserve element, [g d-2]
  real(r8),target,allocatable ::  WTHSKE(:,:,:,:)                    !canopy husk element, [g d-2]
  real(r8),target,allocatable ::  WTEARE(:,:,:,:)                    !canopy ear element, [g d-2]
  real(r8),target,allocatable ::  WTGRE(:,:,:,:)                     !canopy grain element, [g d-2]
  real(r8),target,allocatable ::  WTLS(:,:,:)                        !canopy leaf + sheath C, [g d-2]
  real(r8),target,allocatable ::  ARLFV(:,:,:,:)                     !canopy layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CNET(:,:,:)                        !canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  WGLFV(:,:,:,:)                     !canopy layer leaf C, [g d-2]
  real(r8),target,allocatable ::  EPOOLP(:,:,:,:)                    !canopy nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CEPOLP(:,:,:,:)                    !canopy nonstructural element concentration, [g d-2]
  real(r8),target,allocatable ::  ARSTV(:,:,:,:)                     !canopy layer stem area, [m2 d-2]
  real(r8),target,allocatable ::  EPOLNP(:,:,:,:)                    !canopy nodule nonstructural element, [g d-2]
  real(r8),target,allocatable ::  WVSTKB(:,:,:,:)                    !branch active stalk C, [g d-2]
  real(r8),target,allocatable ::  EPOOL(:,:,:,:,:)                   !branch nonstructural element, [g d-2]
  real(r8),target,allocatable ::  WTLSB(:,:,:,:)                     !branch leaf + sheath C, [g d-2]
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
  real(r8),target,allocatable ::  ARLFL(:,:,:,:,:,:)                 !layer leaf area, [m2 d-2]
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
  allocate(RSMH(JP,JY,JX));     RSMH=0._r8
  allocate(RCS(JP,JY,JX));      RCS=0._r8
  allocate(RC(JP,JY,JX));       RC=0._r8
  allocate(RSMN(JP,JY,JX));     RSMN=0._r8
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
  allocate(RAD1(JP,JY,JX));     RAD1=0._r8
  allocate(THRM1(JP,JY,JX));    THRM1=0._r8
  allocate(RADC(JP,JY,JX));     RADC=0._r8
  allocate(RADP(JP,JY,JX));     RADP=0._r8
  allocate(FRADP(JP,JY,JX));    FRADP=0._r8
  allocate(TAU0(JC+1,JY,JX));   TAU0=0._r8
  allocate(TAUS(JC+1,JY,JX));   TAUS=0._r8
  allocate(FRADG(JY,JX));       FRADG=0._r8
  allocate(RADG(JY,JX));        RADG=0._r8
  allocate(THRMCX(JY,JX));      THRMCX=0._r8
  allocate(THRMGX(JY,JX));      THRMGX=0._r8
  allocate(TVOLWC(JY,JX));      TVOLWC=0._r8
  allocate(TFLWCI(JY,JX));      TFLWCI=0._r8
  allocate(TFLWC(JY,JX));       TFLWC=0._r8
  allocate(EFLXC(JP,JY,JX));    EFLXC=0._r8
  allocate(SFLXC(JP,JY,JX));    SFLXC=0._r8
  allocate(HFLXC(JP,JY,JX));    HFLXC=0._r8
  allocate(ENGYX(JP,JY,JX));    ENGYX=0._r8
  allocate(VHCPC(JP,JY,JX));    VHCPC=0._r8
  allocate(PSILT(JP,JY,JX));    PSILT=0._r8
  allocate(PSILG(JP,JY,JX));    PSILG=0._r8
  allocate(PSILO(JP,JY,JX));    PSILO=0._r8
  allocate(RA(JP,JY,JX));       RA=0._r8
  allocate(EP(JP,JY,JX));       EP=0._r8
  allocate(EVAPC(JP,JY,JX));    EVAPC=0._r8
  allocate(VOLWP(JP,JY,JX));    VOLWP=0._r8
  allocate(TEVAPP(JY,JX));      TEVAPP=0._r8
  allocate(TEVAPC(JY,JX));      TEVAPC=0._r8
  allocate(TENGYC(JY,JX));      TENGYC=0._r8
  allocate(THFLXC(JY,JX));      THFLXC=0._r8
  allocate(TVOLWP(JY,JX));      TVOLWP=0._r8
  allocate(THRMC(JY,JX));       THRMC=0._r8
  allocate(ALBR(JP,JY,JX));     ALBR=0._r8
  allocate(TAUR(JP,JY,JX));     TAUR=0._r8
  allocate(FLWC(JP,JY,JX));     FLWC=0._r8
  allocate(VOLWC(JP,JY,JX));    VOLWC=0._r8
  allocate(TKC(JP,JY,JX));      TKC=0._r8
  allocate(TCC(JP,JY,JX));      TCC=0._r8
  allocate(DTKC(JP,JY,JX));     DTKC=0._r8
  allocate(TKCZ(JP,JY,JX));     TKCZ=0._r8
  allocate(CPOOL3(JNODS,JBR,JP,JY,JX));CPOOL3=0._r8
  allocate(RSETE(npelms,JP,JY,JX));    RSETE=0._r8
  allocate(WGLFT(JC,JY,JX));    WGLFT=0._r8
  allocate(CFOPE(0:Jlitgrp,jsken,npelms,JP,JY,JX));CFOPE=0._r8
  allocate(WTSHTE(npelms,JP,JY,JX)); WTSHTE=0._r8
  allocate(WTLFE(npelms,JP,JY,JX));  WTLFE=0._r8
  allocate(WTSHEE(npelms,JP,JY,JX)); WTSHEE=0._r8
  allocate(WTSTKE(npelms,JP,JY,JX)); WTSTKE=0._r8
  allocate(WVSTK(JP,JY,JX));    WVSTK=0._r8
  allocate(WTRSVE(npelms,JP,JY,JX));    WTRSVE=0._r8
  allocate(WTHSKE(npelms,JP,JY,JX));    WTHSKE=0._r8
  allocate(WTEARE(npelms,JP,JY,JX));    WTEARE=0._r8
  allocate(WTGRE(npelms,JP,JY,JX));     WTGRE=0._r8
  allocate(WTLS(JP,JY,JX));     WTLS=0._r8
  allocate(ARLFV(JC,JP,JY,JX)); ARLFV=0._r8
  allocate(CNET(JP,JY,JX));     CNET=0._r8
  allocate(WGLFV(JC,JP,JY,JX)); WGLFV=0._r8
  allocate(EPOOLP(npelms,JP,JY,JX));   EPOOLP=0._r8
  allocate(CEPOLP(npelms,JP,JY,JX));   CEPOLP=0._r8
  allocate(ARSTV(JC,JP,JY,JX)); ARSTV=0._r8
  allocate(EPOLNP(npelms,JP,JY,JX));   EPOLNP=0._r8
  allocate(WVSTKB(JBR,JP,JY,JX));WVSTKB=0._r8
  allocate(EPOOL(JBR,npelms,JP,JY,JX)); EPOOL=0._r8
  allocate(WTLSB(JBR,JP,JY,JX)); WTLSB=0._r8
  allocate(WTSHTBE(JBR,npelms,JP,JY,JX));WTSHTBE=0._r8
  allocate(WTLFBE(JBR,npelms,JP,JY,JX)); WTLFBE=0._r8
  allocate(WTSHEBE(JBR,npelms,JP,JY,JX));WTSHEBE=0._r8
  allocate(WTSTKBE(JBR,npelms,JP,JY,JX));WTSTKBE=0._r8
  allocate(WTRSVBE(JBR,npelms,JP,JY,JX));WTRSVBE=0._r8
  allocate(WTHSKBE(JBR,npelms,JP,JY,JX));WTHSKBE=0._r8
  allocate(WTEARBE(JBR,npelms,JP,JY,JX));WTEARBE=0._r8
  allocate(WTGRBE(JBR,npelms,JP,JY,JX)); WTGRBE=0._r8
  allocate(CEPOLB(JBR,npelms,JP,JY,JX));CEPOLB=0._r8
  allocate(EPOLNB(JBR,npelms,JP,JY,JX));EPOLNB=0._r8
  allocate(WTNDBE(JBR,npelms,JP,JY,JX)); WTNDBE=0._r8
  allocate(WGSHEXE(JBR,npelms,JP,JY,JX));WGSHEXE=0._r8
  allocate(WTSTXBE(JBR,npelms,JP,JY,JX));WTSTXBE=0._r8
  allocate(WGLFEX(JBR,npelms,JP,JY,JX)); WGLFEX=0._r8
  allocate(WGLFE(0:JNODS,JBR,npelms,JP,JY,JX));WGLFE=0._r8
  allocate(WGSHE(0:JNODS,JBR,npelms,JP,JY,JX));WGSHE=0._r8
  allocate(WGNODE(0:JNODS,JBR,npelms,JP,JY,JX));WGNODE=0._r8
  allocate(WGLFLE(JC,0:JNODS,JBR,npelms,JP,JY,JX));WGLFLE=0._r8
  allocate(ARLFL(JC,0:JNODS,JBR,JP,JY,JX));ARLFL=0._r8
  allocate(WSLF(0:JNODS,JBR,JP,JY,JX));WSLF=0._r8
  allocate(WSSHE(0:JNODS,JBR,JP,JY,JX));WSSHE=0._r8
  allocate(CCPLNP(JP,JY,JX));   CCPLNP=0._r8
  allocate(GRWTB(JBR,JP,JY,JX)); GRWTB=0._r8
  allocate(WTSTDE(jsken,npelms,JP,JY,JX)); WTSTDE=0._r8
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
  call destroy(RSMH)
  call destroy(RCS)
  call destroy(RC)
  call destroy(RSMN)
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
  call destroy(RAD1)
  call destroy(THRM1)
  call destroy(RADC)
  call destroy(RADP)
  call destroy(FRADP)
  call destroy(TAU0)
  call destroy(TAUS)
  call destroy(FRADG)
  call destroy(RADG)
  call destroy(THRMCX)
  call destroy(THRMGX)
  call destroy(TVOLWC)
  call destroy(TFLWCI)
  call destroy(TFLWC)
  call destroy(EFLXC)
  call destroy(SFLXC)
  call destroy(HFLXC)
  call destroy(ENGYX)
  call destroy(VHCPC)
  call destroy(PSILT)
  call destroy(PSILG)
  call destroy(PSILO)
  call destroy(RA)
  call destroy(EP)
  call destroy(EVAPC)
  call destroy(VOLWP)
  call destroy(TEVAPP)
  call destroy(TEVAPC)
  call destroy(TENGYC)
  call destroy(THFLXC)
  call destroy(TVOLWP)
  call destroy(THRMC)
  call destroy(ALBR)
  call destroy(TAUR)
  call destroy(FLWC)
  call destroy(VOLWC)
  call destroy(TKC)
  call destroy(TCC)
  call destroy(DTKC)
  call destroy(TKCZ)
  call destroy(CPOOL3)
  call destroy(RSETE)
  call destroy(WGLFT)
  call destroy(CFOPE)
  call destroy(WTSHTE)
  call destroy(WTLFE)
  call destroy(WTSHEE)
  call destroy(WTSTKE)
  call destroy(WVSTK)
  call destroy(WTRSVE)
  call destroy(WTHSKE)
  call destroy(WTEARE)
  call destroy(WTGRE)
  call destroy(WTLS)
  call destroy(ARLFV)
  call destroy(CNET)
  call destroy(WGLFV)
  call destroy(EPOOLP)
  call destroy(CEPOLP)
  call destroy(ARSTV)
  call destroy(EPOLNP)
  call destroy(WVSTKB)
  call destroy(EPOOL)
  call destroy(WTLSB)
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
  call destroy(ARLFL)
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
