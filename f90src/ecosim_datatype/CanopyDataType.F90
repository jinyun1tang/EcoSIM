module CanopyDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use ElmIDMod
  use EcoSIMConfig, only : jsken => jskenc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  ALBP(:,:,:)                        !canopy PAR albedo , [-]
  real(r8),allocatable ::  TAUP(:,:,:)                        !canopy PAR transmissivity , [-]
  real(r8),allocatable ::  ABSR(:,:,:)                        !canopy shortwave absorptivity , [-]
  real(r8),allocatable ::  ABSP(:,:,:)                        !canopy PAR absorptivity , [-]
  real(r8),allocatable ::  RSMX(:,:,:)                        !maximum stomatal resistance to vapor, [s m-1]
  real(r8),allocatable ::  RCMX(:,:,:)                        !maximum stomatal resistance to CO2, [s h-1]
  real(r8),allocatable ::  RSMH(:,:,:)                        !maximum stomatal resistance to vapor, [s h-1]
  real(r8),allocatable ::  RCS(:,:,:)                         !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8),allocatable ::  RC(:,:,:)                          !canopy stomatal resistance, [h m-1]
  real(r8),allocatable ::  RSMN(:,:,:)                        !canopy minimum stomatal resistance, [s m-1]
  real(r8),allocatable ::  RAC(:,:)                           !canopy boundary layer resistance, [m h-1]
  real(r8),allocatable ::  O2I(:,:,:)                         !leaf gaseous O2 concentration, [umol m-3]
  real(r8),allocatable ::  CO2I(:,:,:)                        !leaf gaseous CO2 concentration, [umol m-3]
  real(r8),allocatable ::  FMOL(:,:,:)                        !total gas concentration, [mol m-3]
  real(r8),allocatable ::  DCO2(:,:,:)                        !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8),allocatable ::  CO2Q(:,:,:)                        !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8),allocatable ::  CO2L(:,:,:)                        !leaf aqueous CO2 concentration, [uM]
  real(r8),allocatable ::  O2L(:,:,:)                         !leaf aqueous O2 concentration, [uM]
  real(r8),allocatable ::  SCO2(:,:,:)                        !leaf CO2 solubility, [uM /umol mol-1]
  real(r8),allocatable ::  SO2(:,:,:)                         !leaf O2 solubility, [uM /umol mol-1]
  real(r8),allocatable ::  XKCO2L(:,:,:)                      !leaf aqueous CO2 Km no O2, [uM]
  real(r8),allocatable ::  XKCO2O(:,:,:)                      !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8),allocatable ::  CHILL(:,:,:)                       !chilling effect on CO2 fixation, [-]
  real(r8),allocatable ::  VCGRO(:,:,:,:,:)                   !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),allocatable ::  VGRO(:,:,:,:,:)                    !carboxylation rate, [umol m-2 s-1]
  real(r8),allocatable ::  COMPL(:,:,:,:,:)                   !CO2 compensation point, [uM]
  real(r8),allocatable ::  ETGRO(:,:,:,:,:)                   !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),allocatable ::  CBXN(:,:,:,:,:)                    !carboxylation efficiency, [umol umol-1]
  real(r8),allocatable ::  CO2B(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),allocatable ::  VCGR4(:,:,:,:,:)                   !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),allocatable ::  VGRO4(:,:,:,:,:)                   !C4 carboxylation rate, [umol m-2 s-1]
  real(r8),allocatable ::  ETGR4(:,:,:,:,:)                   !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),allocatable ::  CBXN4(:,:,:,:,:)                   !C4 carboxylation efficiency, [umol umol-1]
  real(r8),allocatable ::  CPOOL4(:,:,:,:,:)                  !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8),allocatable ::  HCOB(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),allocatable ::  FDBK(:,:,:,:)                      !branch down-regulation of CO2 fixation, [-]
  real(r8),allocatable ::  FDBK4(:,:,:,:,:)                   !down-regulation of C4 photosynthesis, [-]
  real(r8),allocatable ::  FDBKX(:,:,:,:)                     !down-regulation of C4 photosynthesis, [-]
  real(r8),allocatable ::  CNETX(:,:)                         !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  VCMX(:,:,:)                        !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),allocatable ::  VOMX(:,:,:)                        !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8),allocatable ::  VCMX4(:,:,:)                       !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),allocatable ::  XKCO2(:,:,:)                       !Km for rubisco carboxylase activity, [uM]
  real(r8),allocatable ::  XKO2(:,:,:)                        !Km for rubisco oxygenase activity, [uM]
  real(r8),allocatable ::  XKCO24(:,:,:)                      !Km for PEP carboxylase activity, [uM]
  real(r8),allocatable ::  RUBP(:,:,:)                        !leaf rubisco content, [g g-1]
  real(r8),allocatable ::  PEPC(:,:,:)                        !leaf PEP carboxylase content, [g g-1]
  real(r8),allocatable ::  ETMX(:,:,:)                        !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8),allocatable ::  CHL(:,:,:)                         !leaf C3 chlorophyll content, [g g-1]
  real(r8),allocatable ::  CHL4(:,:,:)                        !leaf C4 chlorophyll content, [g g-1]
  real(r8),allocatable ::  FCO2(:,:,:)                        !Ci:Ca ratio, [-]
  real(r8),allocatable ::  RAD1(:,:,:)                        !canopy net radiation , [MJ d-2 h-1]
  real(r8),allocatable ::  THRM1(:,:,:)                       !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),allocatable ::  RADC(:,:,:)                        !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),allocatable ::  RADP(:,:,:)                        !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),allocatable ::  FRADP(:,:,:)                       !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),allocatable ::  TAU0(:,:,:)                        !fraction of radiation transmitted by canopy layer, [-]
  real(r8),allocatable ::  TAUS(:,:,:)                        !fraction of radiation intercepted by canopy layer, [-]
  real(r8),allocatable ::  FRADG(:,:)                         !fraction of radiation intercepted by ground surface, [-]
  real(r8),allocatable ::  RADG(:,:)                          !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8),allocatable ::  THRMCX(:,:)                        !longwave radiation emitted by canopy, [MJ m-2 h-1]
  real(r8),allocatable ::  THRMGX(:,:)                        !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8),allocatable ::  TVOLWC(:,:)                        !canopy surface water content, [m3 d-2]
  real(r8),allocatable ::  TFLWCI(:,:)                        !net ice transfer to canopy, [MJ d-2 t-1]
  real(r8),allocatable ::  TFLWC(:,:)                         !net water transfer to canopy, [MJ d-2 t-1]
  real(r8),allocatable ::  EFLXC(:,:,:)                       !canopy latent heat flux, [MJ d-2 h-1]
  real(r8),allocatable ::  SFLXC(:,:,:)                       !canopy sensible heat flux, [MJ d-2 h-1]
  real(r8),allocatable ::  HFLXC(:,:,:)                       !canopy storage heat flux, [MJ d-2 h-1]
  real(r8),allocatable ::  ENGYX(:,:,:)                       !canopy heat storage from previous time step, [MJ d-2]
  real(r8),allocatable ::  VHCPC(:,:,:)                       !canopy heat capacity, [MJ d-2 K-1]
  real(r8),allocatable ::  PSILT(:,:,:)                       !canopy total water potential , [Mpa]
  real(r8),allocatable ::  PSILG(:,:,:)                       !canopy turgor water potential, [Mpa]
  real(r8),allocatable ::  PSILO(:,:,:)                       !canopy osmotic water potential, [Mpa]
  real(r8),allocatable ::  RA(:,:,:)                          !canopy boundary layer resistance, [h m-1]
  real(r8),allocatable ::  EP(:,:,:)                          !canopy transpiration, [m2 d-2 h-1]
  real(r8),allocatable ::  EVAPC(:,:,:)                       !canopy evaporation, [m2 d-2 h-1]
  real(r8),allocatable ::  VOLWP(:,:,:)                       !canopy water content, [m3 d-2]
  real(r8),allocatable ::  TEVAPP(:,:)                        !total canopy evaporation + transpiration, [m3 d-2]
  real(r8),allocatable ::  TEVAPC(:,:)                        !total canopy evaporation, [m3 d-2]
  real(r8),allocatable ::  TENGYC(:,:)                        !total canopy heat content, [MJ  d-2]
  real(r8),allocatable ::  THFLXC(:,:)                        !total canopy heat flux, [MJ  d-2]
  real(r8),allocatable ::  TVOLWP(:,:)                        !total canopy water content, [m3 d-2]
  real(r8),allocatable ::  THRMC(:,:)                         !total canopy LW emission, [MJ d-2 h-1]
  real(r8),allocatable ::  ALBR(:,:,:)                        !canopy shortwave albedo , [-]
  real(r8),allocatable ::  TAUR(:,:,:)                        !canopy shortwave transmissivity , [-]
  real(r8),allocatable ::  FLWC(:,:,:)                        !water flux into canopy, [m3 d-2 h-1]
  real(r8),allocatable ::  VOLWC(:,:,:)                       !canopy surface water content, [m3 d-2]
  real(r8),allocatable ::  TKC(:,:,:)                         !canopy temperature, [K]
  real(r8),allocatable ::  TCC(:,:,:)                         !canopy temperature, [oC]
  real(r8),allocatable ::  DTKC(:,:,:)                        !change in canopy temperature, [K]
  real(r8),allocatable ::  TKCZ(:,:,:)                        !canopy temperature, [K]
  real(r8),allocatable ::  CPOOL3(:,:,:,:,:)                  !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8),allocatable ::  RSETE(:,:,:,:)                     !effect of canopy element status on seed set , []
  real(r8),allocatable ::  WGLFT(:,:,:)                       !total leaf mass, [g d-2]
  real(r8),allocatable ::  CFOPC(:,:,:,:,:)                   !litter kinetic fraction, [-]
  real(r8),allocatable ::  CFOPN(:,:,:,:,:)                   !litterfall kinetic N fraction, [-]
  real(r8),allocatable ::  CFOPP(:,:,:,:,:)                   !litter P kinetic fraction, [-]
  real(r8),allocatable ::  WTSHTE(:,:,:,:)                    !canopy shoot element, [g d-2]
  real(r8),allocatable ::  WTLFE(:,:,:,:)                     !canopy leaf element, [g d-2]
  real(r8),allocatable ::  WTSHEE(:,:,:,:)                    !canopy sheath element , [g d-2]
  real(r8),allocatable ::  WTSTKE(:,:,:,:)                    !canopy stalk element, [g d-2]
  real(r8),allocatable ::  WVSTK(:,:,:)                       !canopy active stalk C, [g d-2]
  real(r8),allocatable ::  WTRSVE(:,:,:,:)                    !canopy reserve element, [g d-2]
  real(r8),allocatable ::  WTHSKE(:,:,:,:)                    !canopy husk element, [g d-2]
  real(r8),allocatable ::  WTEARE(:,:,:,:)                    !canopy ear element, [g d-2]
  real(r8),allocatable ::  WTGRE(:,:,:,:)                     !canopy grain element, [g d-2]
  real(r8),allocatable ::  WTLS(:,:,:)                        !canopy leaf + sheath C, [g d-2]
  real(r8),allocatable ::  ARLFV(:,:,:,:)                     !canopy layer leaf area, [m2 d-2]
  real(r8),allocatable ::  CNET(:,:,:)                        !canopy net CO2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  WGLFV(:,:,:,:)                     !canopy layer leaf C, [g d-2]
  real(r8),allocatable ::  EPOOLP(:,:,:,:)                    !canopy nonstructural element, [g d-2]
  real(r8),allocatable ::  CEPOLP(:,:,:,:)                    !canopy nonstructural element concentration, [g d-2]
  real(r8),allocatable ::  ARSTV(:,:,:,:)                     !canopy layer stem area, [m2 d-2]
  real(r8),allocatable ::  EPOLNP(:,:,:,:)                    !canopy nodule nonstructural element, [g d-2]
  real(r8),allocatable ::  WVSTKB(:,:,:,:)                    !branch active stalk C, [g d-2]
  real(r8),allocatable ::  EPOOL(:,:,:,:,:)                   !branch nonstructural element, [g d-2]
  real(r8),allocatable ::  WTLSB(:,:,:,:)                     !branch leaf + sheath C, [g d-2]
  real(r8),allocatable ::  WTSHTBE(:,:,:,:,:)                 !branch shoot C, [g d-2]
  real(r8),allocatable ::  WTLFB(:,:,:,:)                     !branch leaf C, [g d-2]
  real(r8),allocatable ::  WTSHEBE(:,:,:,:,:)                 !branch sheath element , [g d-2]
  real(r8),allocatable ::  WTSTKB(:,:,:,:)                    !branch stalk C, [g d-2]
  real(r8),allocatable ::  WTRSVB(:,:,:,:)                    !branch reserve C, [g d-2]
  real(r8),allocatable ::  WTHSKB(:,:,:,:)                    !branch husk C, [g d-2]
  real(r8),allocatable ::  WTEARB(:,:,:,:)                    !branch ear C, [g d-2]
  real(r8),allocatable ::  WTGRB(:,:,:,:)                     !branch grain C, [g d-2]
  real(r8),allocatable ::  CEPOLB(:,:,:,:,:)                    !branch nonstructural C concentration, [g d-2]
  real(r8),allocatable ::  EPOLNB(:,:,:,:,:)                  !branch nodule nonstructural C, [g d-2]
  real(r8),allocatable ::  WTNDB(:,:,:,:)                     !branch nodule C, [g d-2]
  real(r8),allocatable ::  WGSHEX(:,:,:,:)                    !branch sheath structural C, [g d-2]
  real(r8),allocatable ::  WTSTXB(:,:,:,:)                    !branch stalk structural C, [g d-2]
  real(r8),allocatable ::  WGLFX(:,:,:,:)                     !branch leaf structural C, [g d-2]
  real(r8),allocatable ::  WGLFPX(:,:,:,:)                    !branch leaf structural P, [g d-2]
  real(r8),allocatable ::  WGSHPX(:,:,:,:)                    !branch sheath structural P, [g d-2]
  real(r8),allocatable ::  WTSTXP(:,:,:,:)                    !branch stalk structural P, [g d-2]
  real(r8),allocatable ::  WTLFBN(:,:,:,:)                    !branch leaf N, [g d-2]
  real(r8),allocatable ::  WTSTBN(:,:,:,:)                    !branch stalk N, [g d-2]
  real(r8),allocatable ::  WTRSBN(:,:,:,:)                    !branch reserve N, [g d-2]
  real(r8),allocatable ::  WTHSBN(:,:,:,:)                    !branch husk N, [g d-2]
  real(r8),allocatable ::  WTEABN(:,:,:,:)                    !branch ear N, [g d-2]
  real(r8),allocatable ::  WTGRBN(:,:,:,:)                    !branch grain N, [g d-2]
  real(r8),allocatable ::  WTNDBN(:,:,:,:)                    !branch nodule N, [g d-2]
  real(r8),allocatable ::  WGSHNX(:,:,:,:)                    !branch sheath structural N, [g d-2]
  real(r8),allocatable ::  WTSTXN(:,:,:,:)                    !branch stalk structural N, [g d-2]
  real(r8),allocatable ::  WGLFNX(:,:,:,:)                    !branch leaf structural N, [g d-2]
  real(r8),allocatable ::  WTLFBP(:,:,:,:)                    !branch leaf P, [g d-2]
  real(r8),allocatable ::  WTSTBP(:,:,:,:)                    !branch stalk P, [g d-2]
  real(r8),allocatable ::  WTRSBP(:,:,:,:)                    !branch reserve P, [g d-2]
  real(r8),allocatable ::  WTHSBP(:,:,:,:)                    !branch husk P, [g d-2]
  real(r8),allocatable ::  WTEABP(:,:,:,:)                    !branch ear P, [g d-2]
  real(r8),allocatable ::  WTGRBP(:,:,:,:)                    !branch grain P, [g d-2]
  real(r8),allocatable ::  WTNDBP(:,:,:,:)                    !branch nodule P, [g d-2]
  real(r8),allocatable ::  WGLF(:,:,:,:,:)                    !leaf C, [g d-2]
  real(r8),allocatable ::  WGSHE(:,:,:,:,:)                   !sheath C , [g d-2]
  real(r8),allocatable ::  WGNODE(:,:,:,:,:)                  !internode C, [g d-2]
  real(r8),allocatable ::  WGLFL(:,:,:,:,:,:)                 !layer leaf C, [g d-2]
  real(r8),allocatable ::  ARLFL(:,:,:,:,:,:)                 !layer leaf area, [m2 d-2]
  real(r8),allocatable ::  WGLFLN(:,:,:,:,:,:)                !layer leaf N, [g d-2]
  real(r8),allocatable ::  WSLF(:,:,:,:,:)                    !layer leaf protein C, [g d-2]
  real(r8),allocatable ::  WSSHE(:,:,:,:,:)                   !layer sheath protein C, [g d-2]
  real(r8),allocatable ::  WGLFN(:,:,:,:,:)                   !leaf N, [g d-2]
  real(r8),allocatable ::  WGSHN(:,:,:,:,:)                   !sheath N, [g d-2]
  real(r8),allocatable ::  WGNODN(:,:,:,:,:)                  !internode N, [g d-2]
  real(r8),allocatable ::  CCPLNP(:,:,:)                      !nodule nonstructural C, [g d-2]
  real(r8),allocatable ::  GRWTB(:,:,:,:)                     !maximum grain C during grain fill, [g d-2]
  real(r8),allocatable ::  WTSTDE(:,:,:,:,:)                  !standing dead element fraction, [g d-2]
  real(r8),allocatable ::  WTSTGE(:,:,:,:)                    !standing dead element, [g d-2]
  real(r8),allocatable ::  WGLFP(:,:,:,:,:)                   !leaf P, [g d-2]
  real(r8),allocatable ::  WGSHP(:,:,:,:,:)                   !sheath P, [g d-2]
  real(r8),allocatable ::  WGNODP(:,:,:,:,:)                  !nodule P, [g d-2]
  real(r8),allocatable ::  WGLFLP(:,:,:,:,:,:)                !leaf layer P, [g d-2]
  real(r8),allocatable ::  WTRVE(:,:,:,:)                     !plant stored nonstructural element, [g d-2]
  real(r8),allocatable ::  WTRVX(:,:,:)                       !plant stored nonstructural C at planting, [g d-2]
  REAL(R8),allocatable ::  WTSHTA(:,:,:)                      !landscape average canopy shoot C, [g d-2]
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
  allocate(VCGRO(JNODS,JC,JP,JY,JX));VCGRO=0._r8
  allocate(VGRO(JNODS,JC,JP,JY,JX));VGRO=0._r8
  allocate(COMPL(JNODS,JC,JP,JY,JX));COMPL=0._r8
  allocate(ETGRO(JNODS,JC,JP,JY,JX));ETGRO=0._r8
  allocate(CBXN(JNODS,JC,JP,JY,JX));CBXN=0._r8
  allocate(CO2B(JNODS,JC,JP,JY,JX));CO2B=0._r8
  allocate(VCGR4(JNODS,JC,JP,JY,JX));VCGR4=0._r8
  allocate(VGRO4(JNODS,JC,JP,JY,JX));VGRO4=0._r8
  allocate(ETGR4(JNODS,JC,JP,JY,JX));ETGR4=0._r8
  allocate(CBXN4(JNODS,JC,JP,JY,JX));CBXN4=0._r8
  allocate(CPOOL4(JNODS,JC,JP,JY,JX));CPOOL4=0._r8
  allocate(HCOB(JNODS,JC,JP,JY,JX));HCOB=0._r8
  allocate(FDBK(JC,JP,JY,JX));  FDBK=0._r8
  allocate(FDBK4(JNODS,JC,JP,JY,JX));FDBK4=0._r8
  allocate(FDBKX(JC,JP,JY,JX)); FDBKX=0._r8
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
  allocate(CPOOL3(JNODS,JC,JP,JY,JX));CPOOL3=0._r8
  allocate(RSETE(npelms,JP,JY,JX));    RSETE=0._r8
  allocate(WGLFT(JC,JY,JX));    WGLFT=0._r8
  allocate(CFOPC(0:5,4,JP,JY,JX));CFOPC=0._r8
  allocate(CFOPN(0:5,4,JP,JY,JX));CFOPN=0._r8
  allocate(CFOPP(0:5,4,JP,JY,JX));CFOPP=0._r8
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
  allocate(WVSTKB(JC,JP,JY,JX));WVSTKB=0._r8
  allocate(EPOOL(JC,npelms,JP,JY,JX)); EPOOL=0._r8
  allocate(WTLSB(JC,JP,JY,JX)); WTLSB=0._r8
  allocate(WTSHTBE(JC,npelms,JP,JY,JX));WTSHTBE=0._r8
  allocate(WTLFB(JC,JP,JY,JX)); WTLFB=0._r8
  allocate(WTSHEBE(JC,npelms,JP,JY,JX));WTSHEBE=0._r8
  allocate(WTSTKB(JC,JP,JY,JX));WTSTKB=0._r8
  allocate(WTRSVB(JC,JP,JY,JX));WTRSVB=0._r8
  allocate(WTHSKB(JC,JP,JY,JX));WTHSKB=0._r8
  allocate(WTEARB(JC,JP,JY,JX));WTEARB=0._r8
  allocate(WTGRB(JC,JP,JY,JX)); WTGRB=0._r8
  allocate(CEPOLB(JC,npelms,JP,JY,JX));CEPOLB=0._r8
  allocate(EPOLNB(JC,npelms,JP,JY,JX));EPOLNB=0._r8
  allocate(WTNDB(JC,JP,JY,JX)); WTNDB=0._r8
  allocate(WGSHEX(JC,JP,JY,JX));WGSHEX=0._r8
  allocate(WTSTXB(JC,JP,JY,JX));WTSTXB=0._r8
  allocate(WGLFX(JC,JP,JY,JX)); WGLFX=0._r8
  allocate(WGLFPX(JC,JP,JY,JX));WGLFPX=0._r8
  allocate(WGSHPX(JC,JP,JY,JX));WGSHPX=0._r8
  allocate(WTSTXP(JC,JP,JY,JX));WTSTXP=0._r8
  allocate(WTLFBN(JC,JP,JY,JX));WTLFBN=0._r8
  allocate(WTSTBN(JC,JP,JY,JX));WTSTBN=0._r8
  allocate(WTRSBN(JC,JP,JY,JX));WTRSBN=0._r8
  allocate(WTHSBN(JC,JP,JY,JX));WTHSBN=0._r8
  allocate(WTEABN(JC,JP,JY,JX));WTEABN=0._r8
  allocate(WTGRBN(JC,JP,JY,JX));WTGRBN=0._r8
  allocate(WTNDBN(JC,JP,JY,JX));WTNDBN=0._r8
  allocate(WGSHNX(JC,JP,JY,JX));WGSHNX=0._r8
  allocate(WTSTXN(JC,JP,JY,JX));WTSTXN=0._r8
  allocate(WGLFNX(JC,JP,JY,JX));WGLFNX=0._r8
  allocate(WTLFBP(JC,JP,JY,JX));WTLFBP=0._r8
  allocate(WTSTBP(JC,JP,JY,JX));WTSTBP=0._r8
  allocate(WTRSBP(JC,JP,JY,JX));WTRSBP=0._r8
  allocate(WTHSBP(JC,JP,JY,JX));WTHSBP=0._r8
  allocate(WTEABP(JC,JP,JY,JX));WTEABP=0._r8
  allocate(WTGRBP(JC,JP,JY,JX));WTGRBP=0._r8
  allocate(WTNDBP(JC,JP,JY,JX));WTNDBP=0._r8
  allocate(WGLF(0:JNODS,JC,JP,JY,JX));WGLF=0._r8
  allocate(WGSHE(0:JNODS,JC,JP,JY,JX));WGSHE=0._r8
  allocate(WGNODE(0:JNODS,JC,JP,JY,JX));WGNODE=0._r8
  allocate(WGLFL(JC,0:JNODS,JC,JP,JY,JX));WGLFL=0._r8
  allocate(ARLFL(JC,0:JNODS,JC,JP,JY,JX));ARLFL=0._r8
  allocate(WGLFLN(JC,0:JNODS,JC,JP,JY,JX));WGLFLN=0._r8
  allocate(WSLF(0:JNODS,JC,JP,JY,JX));WSLF=0._r8
  allocate(WSSHE(0:JNODS,JC,JP,JY,JX));WSSHE=0._r8
  allocate(WGLFN(0:JNODS,JC,JP,JY,JX));WGLFN=0._r8
  allocate(WGSHN(0:JNODS,JC,JP,JY,JX));WGSHN=0._r8
  allocate(WGNODN(0:JNODS,JC,JP,JY,JX));WGNODN=0._r8
  allocate(CCPLNP(JP,JY,JX));   CCPLNP=0._r8
  allocate(GRWTB(JC,JP,JY,JX)); GRWTB=0._r8
  allocate(WTSTDE(jsken,npelms,JP,JY,JX)); WTSTDE=0._r8
  allocate(WTSTGE(npelms,JP,JY,JX));    WTSTGE=0._r8
  allocate(WGLFP(0:JNODS,JC,JP,JY,JX));WGLFP=0._r8
  allocate(WGSHP(0:JNODS,JC,JP,JY,JX));WGSHP=0._r8
  allocate(WGNODP(0:JNODS,JC,JP,JY,JX));WGNODP=0._r8
  allocate(WGLFLP(JC,0:JNODS,JC,JP,JY,JX));WGLFLP=0._r8
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
  call destroy(CFOPC)
  call destroy(CFOPN)
  call destroy(CFOPP)
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
  call destroy(WTLFB)
  call destroy(WTSHEBE)
  call destroy(WTSTKB)
  call destroy(WTRSVB)
  call destroy(WTHSKB)
  call destroy(WTEARB)
  call destroy(WTGRB)
  call destroy(CEPOLB)
  call destroy(EPOLNB)
  call destroy(WTNDB)
  call destroy(WGSHEX)
  call destroy(WTSTXB)
  call destroy(WGLFX)
  call destroy(WGLFPX)
  call destroy(WGSHPX)
  call destroy(WTSTXP)
  call destroy(WTLFBN)
  call destroy(WTSTBN)
  call destroy(WTRSBN)
  call destroy(WTHSBN)
  call destroy(WTEABN)
  call destroy(WTGRBN)
  call destroy(WTNDBN)
  call destroy(WGSHNX)
  call destroy(WTSTXN)
  call destroy(WGLFNX)
  call destroy(WTLFBP)
  call destroy(WTSTBP)
  call destroy(WTRSBP)
  call destroy(WTHSBP)
  call destroy(WTEABP)
  call destroy(WTGRBP)
  call destroy(WTNDBP)
  call destroy(WGLF)
  call destroy(WGSHE)
  call destroy(WGNODE)
  call destroy(WGLFL)
  call destroy(ARLFL)
  call destroy(WGLFLN)
  call destroy(WSLF)
  call destroy(WSSHE)
  call destroy(WGLFN)
  call destroy(WGSHN)
  call destroy(WGNODN)
  call destroy(CCPLNP)
  call destroy(GRWTB)
  call destroy(WTSTDE)
  call destroy(WTSTGE)
  call destroy(WGLFP)
  call destroy(WGSHP)
  call destroy(WGNODP)
  call destroy(WGLFLP)
  call destroy(WTRVE)
  call destroy(WTRVX)
  call destroy(WTSHTA)
  end subroutine DestructCanopyData

end module CanopyDataType
