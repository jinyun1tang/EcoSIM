module RootDataType

!
!!
! data types of plant characteristics
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  integer,allocatable ::  NRT(:,:,:)                          !root primary axis number, [-]
  integer,allocatable ::  NINR(:,:,:,:)                       !maximum soil layer number for root axes, [-]
  integer,allocatable ::  IDTHR(:,:,:)                        !flag to detect root system death , [-]
  integer,allocatable ::  NIX(:,:,:)                          !maximum soil layer number for all root axes, [-]
  integer,allocatable ::  NI(:,:,:)                           !maximum soil layer number for all root axes, [-]
  real(r8),allocatable ::  DMRT(:,:,:)                        !root growth yield, [g g-1]
  real(r8),allocatable ::  PR(:,:,:)                          !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(r8),allocatable ::  CWSRT(:,:,:)                       !fraction of remobilizable nonstructural biomass in root, [-]
  real(r8),allocatable ::  DMVL(:,:,:,:)                      !root volume:mass ratio, [m3 g-1]
  real(r8),allocatable ::  RRAD1X(:,:,:,:)                    !root diameter primary axes, [m]
  real(r8),allocatable ::  RRAD2X(:,:,:,:)                    !root diameter secondary axes, [m]
  real(r8),allocatable ::  RTAR1X(:,:,:,:)                    !root cross-sectional area primary axes, [m2]
  real(r8),allocatable ::  RTAR2X(:,:,:,:)                    !root  cross-sectional area  secondary axes, [m2]
  real(r8),allocatable ::  TFN4(:,:,:,:)                      !root layer temperature growth functiom, [-]
  real(r8),allocatable ::  CNRT(:,:,:)                        !root N:C ratio, [g g-1]
  real(r8),allocatable ::  CPRT(:,:,:)                        !root P:C ratio, [g g-1]
  real(r8),allocatable ::  PORT(:,:,:,:)                      !root porosity, [m3 m-3]
  real(r8),allocatable ::  RSRR(:,:,:,:)                      !root radial resistivity, [MPa h m-2]
  real(r8),allocatable ::  RSRA(:,:,:,:)                      !root axial resistivity, [MPa h m-4]
  real(r8),allocatable ::  PTSHT(:,:,:)                       !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8),allocatable ::  UPMXZH(:,:,:,:)                    !maximum root NH4 uptake rate, [g m-2 h-1]
  real(r8),allocatable ::  UPKMZH(:,:,:,:)                    !Km for root NH4 uptake, [g m-3]
  real(r8),allocatable ::  UPMNZH(:,:,:,:)                    !minimum NH4 concentration for root NH4 uptake, [g m-3]
  real(r8),allocatable ::  UPMXZO(:,:,:,:)                    !maximum root NO3 uptake rate, [g m-2 h-1]
  real(r8),allocatable ::  UPKMZO(:,:,:,:)                    !Km for root NO3 uptake, [g m-3]
  real(r8),allocatable ::  UPMNZO(:,:,:,:)                    !minimum NO3 concentration for root NH4 uptake, [g m-3]
  real(r8),allocatable ::  UPMXPO(:,:,:,:)                    !maximum root PO4 uptake rate, [g m-2 h-1]
  real(r8),allocatable ::  UPKMPO(:,:,:,:)                    !Km for root PO4 uptake, [g m-3]
  real(r8),allocatable ::  UPMNPO(:,:,:,:)                    !minimum PO4 concentration for root NH4 uptake, [g m-3]
  real(r8),allocatable ::  RRADP(:,:,:,:)                     !root internal radius, [m]
  real(r8),allocatable ::  CNRTS(:,:,:)                       !root N:C ratio x root growth yield, [-]
  real(r8),allocatable ::  CPRTS(:,:,:)                       !root P:C ratio x root growth yield, [-]
  real(r8),allocatable ::  RRAD1M(:,:,:,:)                    !maximum radius of primary roots, [m]
  real(r8),allocatable ::  RRAD2M(:,:,:,:)                    !maximum radius of secondary roots, [m]
  real(r8),allocatable ::  RTFQ(:,:,:)                        !root brancing frequency, [m-1]
  real(r8),allocatable ::  PORTX(:,:,:,:)                     !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8),allocatable ::  CPPOLR(:,:,:,:,:)                  !root layer nonstructural P concentration, [g g-1]
  real(r8),allocatable ::  WTNDLP(:,:,:,:)                    !root layer nodule P, [g d-2]
  real(r8),allocatable ::  PPOOLN(:,:,:,:)                    !nodule layer nonstructural P, [g d-2]
  real(r8),allocatable ::  ZPOOLN(:,:,:,:)                    !root nodule nonstructural N, [g d-2]
  real(r8),allocatable ::  RTWT1P(:,:,:,:,:)                  !root P primary axes, [g d-2]
  real(r8),allocatable ::  CZPOLR(:,:,:,:,:)                  !root layer nonstructural N concentration, [g g-1]
  real(r8),allocatable ::  WTNDLN(:,:,:,:)                    !root layer nodule N, [g d-2]
  real(r8),allocatable ::  RTWT1N(:,:,:,:,:)                  !root N primary axes, [g d-2]
  real(r8),allocatable ::  CPOOLN(:,:,:,:)                    !root  layer nonstructural N, [g d-2]
  real(r8),allocatable ::  RTLGP(:,:,:,:,:)                   !root layer length per plant, [m p-1]
  real(r8),allocatable ::  RTLG1(:,:,:,:,:,:)                 !root layer length primary axes, [m d-2]
  real(r8),allocatable ::  RTLG2(:,:,:,:,:,:)                 !root layer length secondary axes, [m d-2]
  real(r8),allocatable ::  RTDNP(:,:,:,:,:)                   !root layer length density, [m m-3]
  real(r8),allocatable ::  RTN1(:,:,:,:,:)                    !root layer number primary axes, [d-2]
  real(r8),allocatable ::  RTNL(:,:,:,:,:)                    !root layer number axes, [d-2]
  real(r8),allocatable ::  RTN2(:,:,:,:,:,:)                  !root layer number secondary axes, [d-2]
  real(r8),allocatable ::  RTLGA(:,:,:,:,:)                   !root layer average length, [m]
  real(r8),allocatable ::  RTARP(:,:,:,:,:)                   !root layer area per plant, [m p-1]
  real(r8),allocatable ::  RTVLW(:,:,:,:,:)                   !root layer volume water, [m2 d-2]
  real(r8),allocatable ::  RRAD1(:,:,:,:,:)                   !root layer diameter primary axes, [m ]
  real(r8),allocatable ::  RTVLP(:,:,:,:,:)                   !root layer volume air, [m2 d-2]
  real(r8),allocatable ::  RTDP1(:,:,:,:,:)                   !root layer depth, [m]
  real(r8),allocatable ::  RRAD2(:,:,:,:,:)                   !root layer diameter secondary axes, [m ]
  real(r8),allocatable ::  RTLG1X(:,:,:,:)                    !specific root length primary axes, [m g-1]
  real(r8),allocatable ::  RTLG2X(:,:,:,:)                    !specific root length secondary axes, [m g-1]
  real(r8),allocatable ::  UPWTR(:,:,:,:,:)                   !root water uptake, [m2 d-2 h-1]
  real(r8),allocatable ::  PSIRT(:,:,:,:,:)                   !root total water potential , [Mpa]
  real(r8),allocatable ::  PSIRO(:,:,:,:,:)                   !root osmotic water potential , [Mpa]
  real(r8),allocatable ::  PSIRG(:,:,:,:,:)                   !root turgor water potential , [Mpa]
  real(r8),allocatable ::  CO2A(:,:,:,:,:)                    !root gaseous CO2 content, [g d-2 ]
  real(r8),allocatable ::  OXYA(:,:,:,:,:)                    !root gaseous O2 content, [g d-2 ]
  real(r8),allocatable ::  CH4A(:,:,:,:,:)                    !root gaseous CH4 content, [g d-2 ]
  real(r8),allocatable ::  Z2OA(:,:,:,:,:)                    !root gaseous N2O content, [g d-2 ]
  real(r8),allocatable ::  ZH3A(:,:,:,:,:)                    !root gaseous NH3 content, [g d-2 ]
  real(r8),allocatable ::  CO2P(:,:,:,:,:)                    !root aqueous CO2 content, [g d-2 ]
  real(r8),allocatable ::  OXYP(:,:,:,:,:)                    !root aqueous O2 content, [g d-2 ]
  real(r8),allocatable ::  CH4P(:,:,:,:,:)                    !root aqueous CH4 content, [g d-2 ]
  real(r8),allocatable ::  Z2OP(:,:,:,:,:)                    !root aqueous N2O content, [g d-2 ]
  real(r8),allocatable ::  ZH3P(:,:,:,:,:)                    !root aqueous NH3 content, [g d-2 ]
  real(r8),allocatable ::  TCO2Z(:,:)                         !total root CO2 content, [g d-2]
  real(r8),allocatable ::  TOXYZ(:,:)                         !total root O2 content, [g d-2]
  real(r8),allocatable ::  TCH4Z(:,:)                         !total root CH4 content, [g d-2]
  real(r8),allocatable ::  TN2OZ(:,:)                         !total root N2O content, [g d-2]
  real(r8),allocatable ::  TNH3Z(:,:)                         !total root NH3 content, [g d-2]
  real(r8),allocatable ::  H2GA(:,:,:,:,:)                    !gaseous H2 content of roots, [g d-2]
  real(r8),allocatable ::  H2GP(:,:,:,:,:)                    !aqueous H2 content of roots, [g d-2]
  real(r8),allocatable ::  WTRTA(:,:,:)                       !root C per plant, [g p-1]
  real(r8),allocatable ::  WTRTE(:,:,:,:)                     !plant root element, [g d-2]
  real(r8),allocatable ::  WTRTSE(:,:,:,:)                    !plant root structural element, [g d-2]
  real(r8),allocatable ::  WSRTL(:,:,:,:,:)                   !root layer protein C, [g d-2]
  real(r8),allocatable ::  WTRT1E(:,:,:,:,:,:,:)              !root layer element primary axes, [g d-2]
  real(r8),allocatable ::  WTRT2E(:,:,:,:,:,:,:)              !root layer element secondary axes, [g d-2]
  real(r8),allocatable ::  WTRTD(:,:,:,:,:)                   !root layer C, [g d-2]
  real(r8),allocatable ::  WTNDL(:,:,:,:)                     !root layer nodule mass, [g d-2]
  real(r8),allocatable ::  WTNDE(:,:,:,:)                     !root total nodule mass, [g d-2]
  real(r8),allocatable ::  WTRTL(:,:,:,:,:)                   !root layer structural C, [g d-2]
  real(r8),allocatable ::  EPOOLR(:,:,:,:,:,:)                !root  layer nonstructural element, [g d-2]
  real(r8),allocatable ::  CCPOLR(:,:,:,:,:)                  !root  layer nonstructural C concentration, [g g-1]
  real(r8),allocatable ::  RTWT1(:,:,:,:,:)                   !root C primary axes, [g d-2]
  real(r8),allocatable ::  CWSRTL(:,:,:,:,:)                  !root layer protein C concentration, [g g-1]
!----------------------------------------------------------------------

contains
  subroutine InitRootData

  implicit none
  allocate(NRT(JP,JY,JX));      NRT=0
  allocate(NINR(JC,JP,JY,JX));  NINR=1  !set to one to avoid numerical failure
  allocate(IDTHR(JP,JY,JX));    IDTHR=0
  allocate(NIX(JP,JY,JX));      NIX=0
  allocate(NI(JP,JY,JX));       NI=0
  allocate(DMRT(JP,JY,JX));     DMRT=0._r8
  allocate(PR(JP,JY,JX));       PR=0._r8
  allocate(CWSRT(JP,JY,JX));    CWSRT=0._r8
  allocate(DMVL(2,JP,JY,JX));   DMVL=0._r8
  allocate(RRAD1X(2,JP,JY,JX)); RRAD1X=0._r8
  allocate(RRAD2X(2,JP,JY,JX)); RRAD2X=0._r8
  allocate(RTAR1X(2,JP,JY,JX)); RTAR1X=0._r8
  allocate(RTAR2X(2,JP,JY,JX)); RTAR2X=0._r8
  allocate(TFN4(JZ,JP,JY,JX));  TFN4=0._r8
  allocate(CNRT(JP,JY,JX));     CNRT=0._r8
  allocate(CPRT(JP,JY,JX));     CPRT=0._r8
  allocate(PORT(2,JP,JY,JX));   PORT=0._r8
  allocate(RSRR(2,JP,JY,JX));   RSRR=0._r8
  allocate(RSRA(2,JP,JY,JX));   RSRA=0._r8
  allocate(PTSHT(JP,JY,JX));    PTSHT=0._r8
  allocate(UPMXZH(2,JP,JY,JX)); UPMXZH=0._r8
  allocate(UPKMZH(2,JP,JY,JX)); UPKMZH=0._r8
  allocate(UPMNZH(2,JP,JY,JX)); UPMNZH=0._r8
  allocate(UPMXZO(2,JP,JY,JX)); UPMXZO=0._r8
  allocate(UPKMZO(2,JP,JY,JX)); UPKMZO=0._r8
  allocate(UPMNZO(2,JP,JY,JX)); UPMNZO=0._r8
  allocate(UPMXPO(2,JP,JY,JX)); UPMXPO=0._r8
  allocate(UPKMPO(2,JP,JY,JX)); UPKMPO=0._r8
  allocate(UPMNPO(2,JP,JY,JX)); UPMNPO=0._r8
  allocate(RRADP(2,JP,JY,JX));  RRADP=0._r8
  allocate(CNRTS(JP,JY,JX));    CNRTS=0._r8
  allocate(CPRTS(JP,JY,JX));    CPRTS=0._r8
  allocate(RRAD1M(2,JP,JY,JX)); RRAD1M=0._r8
  allocate(RRAD2M(2,JP,JY,JX)); RRAD2M=0._r8
  allocate(RTFQ(JP,JY,JX));     RTFQ=0._r8
  allocate(PORTX(2,JP,JY,JX));  PORTX=0._r8
  allocate(CPPOLR(2,JZ,JP,JY,JX));CPPOLR=0._r8
  allocate(WTNDLP(JZ,JP,JY,JX));WTNDLP=0._r8
  allocate(PPOOLN(JZ,JP,JY,JX));PPOOLN=0._r8
  allocate(ZPOOLN(JZ,JP,JY,JX));ZPOOLN=0._r8
  allocate(RTWT1P(2,JC,JP,JY,JX));RTWT1P=0._r8
  allocate(CZPOLR(2,JZ,JP,JY,JX));CZPOLR=0._r8
  allocate(WTNDLN(JZ,JP,JY,JX));WTNDLN=0._r8
  allocate(RTWT1N(2,JC,JP,JY,JX));RTWT1N=0._r8
  allocate(CPOOLN(JZ,JP,JY,JX));CPOOLN=0._r8
  allocate(RTLGP(2,JZ,JP,JY,JX));RTLGP=0._r8
  allocate(RTLG1(2,JZ,JC,JP,JY,JX));RTLG1=0._r8
  allocate(RTLG2(2,JZ,JC,JP,JY,JX));RTLG2=0._r8
  allocate(RTDNP(2,JZ,JP,JY,JX));RTDNP=0._r8
  allocate(RTN1(2,JZ,JP,JY,JX));RTN1=0._r8
  allocate(RTNL(2,JZ,JP,JY,JX));RTNL=0._r8
  allocate(RTN2(2,JZ,JC,JP,JY,JX));RTN2=0._r8
  allocate(RTLGA(2,JZ,JP,JY,JX));RTLGA=0._r8
  allocate(RTARP(2,JZ,JP,JY,JX));RTARP=0._r8
  allocate(RTVLW(2,JZ,JP,JY,JX));RTVLW=0._r8
  allocate(RRAD1(2,JZ,JP,JY,JX));RRAD1=0._r8
  allocate(RTVLP(2,JZ,JP,JY,JX));RTVLP=0._r8
  allocate(RTDP1(2,JC,JP,JY,JX));RTDP1=0._r8
  allocate(RRAD2(2,JZ,JP,JY,JX));RRAD2=0._r8
  allocate(RTLG1X(2,JP,JY,JX)); RTLG1X=0._r8
  allocate(RTLG2X(2,JP,JY,JX)); RTLG2X=0._r8
  allocate(UPWTR(2,JZ,JP,JY,JX));UPWTR=0._r8
  allocate(PSIRT(2,JZ,JP,JY,JX));PSIRT=0._r8
  allocate(PSIRO(2,JZ,JP,JY,JX));PSIRO=0._r8
  allocate(PSIRG(2,JZ,JP,JY,JX));PSIRG=0._r8
  allocate(CO2A(2,JZ,JP,JY,JX));CO2A=0._r8
  allocate(OXYA(2,JZ,JP,JY,JX));OXYA=0._r8
  allocate(CH4A(2,JZ,JP,JY,JX));CH4A=0._r8
  allocate(Z2OA(2,JZ,JP,JY,JX));Z2OA=0._r8
  allocate(ZH3A(2,JZ,JP,JY,JX));ZH3A=0._r8
  allocate(CO2P(2,JZ,JP,JY,JX));CO2P=0._r8
  allocate(OXYP(2,JZ,JP,JY,JX));OXYP=0._r8
  allocate(CH4P(2,JZ,JP,JY,JX));CH4P=0._r8
  allocate(Z2OP(2,JZ,JP,JY,JX));Z2OP=0._r8
  allocate(ZH3P(2,JZ,JP,JY,JX));ZH3P=0._r8
  allocate(TCO2Z(JY,JX));       TCO2Z=0._r8
  allocate(TOXYZ(JY,JX));       TOXYZ=0._r8
  allocate(TCH4Z(JY,JX));       TCH4Z=0._r8
  allocate(TN2OZ(JY,JX));       TN2OZ=0._r8
  allocate(TNH3Z(JY,JX));       TNH3Z=0._r8
  allocate(H2GA(2,JZ,JP,JY,JX));H2GA=0._r8
  allocate(H2GP(2,JZ,JP,JY,JX));H2GP=0._r8
  allocate(WTRTA(JP,JY,JX));    WTRTA=0._r8
  allocate(WTRTE(npelms,JP,JY,JX)); WTRTE=0._r8
  allocate(WTRTSE(npelms,JP,JY,JX));   WTRTSE=0._r8
  allocate(WSRTL(2,JZ,JP,JY,JX));WSRTL=0._r8
  allocate(WTRT1E(npelms,2,JZ,JC,JP,JY,JX));WTRT1E=0._r8
  allocate(WTRT2E(npelms,2,JZ,JC,JP,JY,JX));WTRT2E=0._r8
  allocate(WTRTD(2,JZ,JP,JY,JX));WTRTD=0._r8
  allocate(WTNDL(JZ,JP,JY,JX)); WTNDL=0._r8
  allocate(WTNDE(npelms,JP,JY,JX));  WTNDE=0._r8
  allocate(WTRTL(2,JZ,JP,JY,JX));WTRTL=0._r8
  allocate(EPOOLR(npelms,2,JZ,JP,JY,JX));EPOOLR=0._r8
  allocate(CCPOLR(2,JZ,JP,JY,JX));CCPOLR=0._r8
  allocate(RTWT1(2,JC,JP,JY,JX));RTWT1=0._r8
  allocate(CWSRTL(2,JZ,JP,JY,JX));CWSRTL=0._r8
  end subroutine InitRootData

!----------------------------------------------------------------------
  subroutine DestructRootData
  use abortutils, only : destroy
  implicit none
  call destroy(NRT)
  call destroy(NINR)
  call destroy(IDTHR)
  call destroy(NIX)
  call destroy(NI)
  call destroy(DMRT)
  call destroy(PR)
  call destroy(CWSRT)
  call destroy(DMVL)
  call destroy(RRAD1X)
  call destroy(RRAD2X)
  call destroy(RTAR1X)
  call destroy(RTAR2X)
  call destroy(TFN4)
  call destroy(CNRT)
  call destroy(CPRT)
  call destroy(PORT)
  call destroy(RSRR)
  call destroy(RSRA)
  call destroy(PTSHT)
  call destroy(UPMXZH)
  call destroy(UPKMZH)
  call destroy(UPMNZH)
  call destroy(UPMXZO)
  call destroy(UPKMZO)
  call destroy(UPMNZO)
  call destroy(UPMXPO)
  call destroy(UPKMPO)
  call destroy(UPMNPO)
  call destroy(RRADP)
  call destroy(CNRTS)
  call destroy(CPRTS)
  call destroy(RRAD1M)
  call destroy(RRAD2M)
  call destroy(RTFQ)
  call destroy(PORTX)
  call destroy(CPPOLR)
  call destroy(WTNDLP)
  call destroy(PPOOLN)
  call destroy(ZPOOLN)
  call destroy(RTWT1P)
  call destroy(CZPOLR)
  call destroy(WTNDLN)
  call destroy(RTWT1N)
  call destroy(CPOOLN)
  call destroy(RTLGP)
  call destroy(RTLG1)
  call destroy(RTLG2)
  call destroy(RTDNP)
  call destroy(RTN1)
  call destroy(RTNL)
  call destroy(RTN2)
  call destroy(RTLGA)
  call destroy(RTARP)
  call destroy(RTVLW)
  call destroy(RRAD1)
  call destroy(RTVLP)
  call destroy(RTDP1)
  call destroy(RRAD2)
  call destroy(RTLG1X)
  call destroy(RTLG2X)
  call destroy(UPWTR)
  call destroy(PSIRT)
  call destroy(PSIRO)
  call destroy(PSIRG)
  call destroy(CO2A)
  call destroy(OXYA)
  call destroy(CH4A)
  call destroy(Z2OA)
  call destroy(ZH3A)
  call destroy(CO2P)
  call destroy(OXYP)
  call destroy(CH4P)
  call destroy(Z2OP)
  call destroy(ZH3P)
  call destroy(TCO2Z)
  call destroy(TOXYZ)
  call destroy(TCH4Z)
  call destroy(TN2OZ)
  call destroy(TNH3Z)
  call destroy(H2GA)
  call destroy(H2GP)
  call destroy(WTRTA)
  call destroy(WTRTE)
  call destroy(WTRTSE)
  call destroy(WSRTL)
  call destroy(WTRT1E)
  call destroy(WTRT2E)
  call destroy(WTRTD)
  call destroy(WTNDL)
  call destroy(WTNDE)
  call destroy(WTRTL)
  call destroy(EPOOLR)
  call destroy(CCPOLR)
  call destroy(RTWT1)
  call destroy(CWSRTL)
  end subroutine DestructRootData

end module RootDataType
