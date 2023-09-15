module RootDataType

!
!!
! data types of plant characteristics
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  integer,target,allocatable ::  NRT(:,:,:)                          !root primary axis number, [-]
  integer,target,allocatable ::  NINR(:,:,:,:)                       !maximum soil layer number for root axes, [-]
  integer,target,allocatable ::  IDTHR(:,:,:)                        !flag to detect root system death , [-]
  integer,target,allocatable ::  NIXBotRootLayer(:,:,:)              !maximum soil layer number for all root axes, [-]
  integer,target,allocatable ::  NI(:,:,:)                           !maximum soil layer number for all root axes, [-]
  real(r8),target,allocatable ::  DMRT(:,:,:)                        !root growth yield, [g g-1]
  real(r8),target,allocatable ::  PR(:,:,:)                          !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(r8),target,allocatable ::  CWSRT(:,:,:)                       !fraction of remobilizable nonstructural biomass in root, [-]
  real(r8),target,allocatable ::  DMVL(:,:,:,:)                      !root volume:mass ratio, [m3 g-1]
  real(r8),target,allocatable ::  MaxPrimRootRadius1(:,:,:,:)                    !root diameter primary axes, [m]
  real(r8),target,allocatable ::  MaxSecndRootRadius1(:,:,:,:)                    !root diameter secondary axes, [m]
  real(r8),target,allocatable ::  PrimRootXSecArea(:,:,:,:)          !root cross-sectional area primary axes, [m2]
  real(r8),target,allocatable ::  SecndRootXSecArea(:,:,:,:)                    !root  cross-sectional area  secondary axes, [m2]
  real(r8),target,allocatable ::  fTgrowRootP(:,:,:,:)                      !root layer temperature growth functiom, [-]
  real(r8),target,allocatable ::  CNRT(:,:,:)                        !root N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPRT(:,:,:)                        !root P:C ratio, [g g-1]
  real(r8),target,allocatable ::  RootPorosity(:,:,:,:)                      !root porosity, [m3 m-3]
  real(r8),target,allocatable ::  RSRR(:,:,:,:)                      !root radial resistivity, [MPa h m-2]
  real(r8),target,allocatable ::  RSRA(:,:,:,:)                      !root axial resistivity, [MPa h m-4]
  real(r8),target,allocatable ::  PTSHT(:,:,:)                       !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8),target,allocatable ::  UPMXZH(:,:,:,:)                    !maximum root NH4 uptake rate, [g m-2 h-1]
  real(r8),target,allocatable ::  UPKMZH(:,:,:,:)                    !Km for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  UPMNZH(:,:,:,:)                    !minimum NH4 concentration for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  UPMXZO(:,:,:,:)                    !maximum root NO3 uptake rate, [g m-2 h-1]
  real(r8),target,allocatable ::  UPKMZO(:,:,:,:)                    !Km for root NO3 uptake, [g m-3]
  real(r8),target,allocatable ::  UPMNZO(:,:,:,:)                    !minimum NO3 concentration for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  UPMXPO(:,:,:,:)                    !maximum root PO4 uptake rate, [g m-2 h-1]
  real(r8),target,allocatable ::  UPKMPO(:,:,:,:)                    !Km for root PO4 uptake, [g m-3]
  real(r8),target,allocatable ::  UPMNPO(:,:,:,:)                    !minimum PO4 concentration for root NH4 uptake, [g m-3]
  real(r8),target,allocatable ::  RRADP(:,:,:,:)                     !root internal radius, [m]
  real(r8),target,allocatable ::  CNRTS(:,:,:)                       !root N:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  CPRTS(:,:,:)                       !root P:C ratio x root growth yield, [-]
  real(r8),target,allocatable ::  MaxPrimRootRadius(:,:,:,:)                    !maximum radius of primary roots, [m]
  real(r8),target,allocatable ::  MaxSecndRootRadius(:,:,:,:)                    !maximum radius of secondary roots, [m]
  real(r8),target,allocatable ::  RTFQ(:,:,:)                        !root brancing frequency, [m-1]
  real(r8),target,allocatable ::  PORTX(:,:,:,:)                     !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8),target,allocatable ::  EPOOLN(:,:,:,:,:)                  !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  RootLenPerP(:,:,:,:,:)                   !root layer length per plant, [m p-1]
  real(r8),target,allocatable ::  PrimRootLen(:,:,:,:,:,:)                 !root layer length primary axes, [m d-2]
  real(r8),target,allocatable ::  SecndRootLen(:,:,:,:,:,:)                 !root layer length secondary axes, [m d-2]
  real(r8),target,allocatable ::  RootLenDensNLP(:,:,:,:,:)          !root length density in soil layers, [m m-3]
  real(r8),target,allocatable ::  PrimRootXNumL(:,:,:,:,:)                    !root layer number primary axes, [d-2]
  real(r8),target,allocatable ::  SecndRootXNumL(:,:,:,:,:)                    !root layer number axes, [d-2]
  real(r8),target,allocatable ::  RTN2(:,:,:,:,:,:)                  !root layer number secondary axes, [d-2]
  real(r8),target,allocatable ::  AveSecndRootLen(:,:,:,:,:)                   !root layer average length, [m]
  real(r8),target,allocatable ::  RTARP(:,:,:,:,:)                   !root layer area per plant, [m p-1]
  real(r8),target,allocatable ::  RTVLW(:,:,:,:,:)                   !root layer volume water, [m2 d-2]
  real(r8),target,allocatable ::  PrimRootRadius(:,:,:,:,:)                   !root layer diameter primary axes, [m ]
  real(r8),target,allocatable ::  RTVLP(:,:,:,:,:)                   !root layer volume air, [m2 d-2]
  real(r8),target,allocatable ::  PrimRootDepth(:,:,:,:,:)                   !root layer depth, [m]
  real(r8),target,allocatable ::  SecndRootRadius(:,:,:,:,:)                   !root layer diameter secondary axes, [m ]
  real(r8),target,allocatable ::  PrimRootSpecLen(:,:,:,:)                    !specific root length primary axes, [m g-1]
  real(r8),target,allocatable ::  SecndRootSpecLen(:,:,:,:)                    !specific root length secondary axes, [m g-1]
  real(r8),target,allocatable ::  PopPlantRootH2OUptake_vr(:,:,:,:,:)                   !root water uptake, [m2 d-2 h-1]
  real(r8),target,allocatable ::  PSIRoot(:,:,:,:,:)                   !root total water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootOSMO(:,:,:,:,:)                   !root osmotic water potential , [Mpa]
  real(r8),target,allocatable ::  PSIRootTurg(:,:,:,:,:)                   !root turgor water potential , [Mpa]
  real(r8),target,allocatable ::  trcg_rootml(:,:,:,:,:,:)           !root gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  trcs_rootml(:,:,:,:,:,:)           !root dissolved gaseous tracer content [g d-2]
  real(r8),target,allocatable ::  TRootGasLoss_disturb(:,:,:)                 !total root gas content, [g d-2]
  real(r8),target,allocatable ::  WTRTA(:,:,:)                       !root C per plant, [g p-1]
  real(r8),target,allocatable ::  WTRTE(:,:,:,:)                     !plant root element, [g d-2]
  real(r8),target,allocatable ::  WTRTSE(:,:,:,:)                    !plant root structural element, [g d-2]
  real(r8),target,allocatable ::  WSRTL(:,:,:,:,:)                   !root layer protein C, [g d-2]
  real(r8),target,allocatable ::  WTRT1E(:,:,:,:,:,:,:)              !root layer element primary axes, [g d-2]
  real(r8),target,allocatable ::  WTRT2E(:,:,:,:,:,:,:)              !root layer element secondary axes, [g d-2]
  real(r8),target,allocatable ::  RootCPZR(:,:,:,:,:)                   !root layer C, [g d-2]
  real(r8),target,allocatable ::  WTNDLE(:,:,:,:,:)                  !root layer nodule element, [g d-2]
  real(r8),target,allocatable ::  WTNDE(:,:,:,:)                     !root total nodule mass, [g d-2]
  real(r8),target,allocatable ::  WTRTL(:,:,:,:,:)                   !root layer structural C, [g d-2]
  real(r8),target,allocatable ::  EPOOLR(:,:,:,:,:,:)                !root  layer nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CEPOLR(:,:,:,:,:,:)                !root  layer nonstructural element concentration, [g g-1]
  real(r8),target,allocatable ::  RTWT1E(:,:,:,:,:,:)                   !root C primary axes, [g d-2]
  real(r8),target,allocatable ::  CWSRTL(:,:,:,:,:)                  !root layer protein C concentration, [g g-1]
!----------------------------------------------------------------------

contains
  subroutine InitRootData(jroots)

  implicit none
  integer, intent(in) :: jroots

  allocate(NRT(JP,JY,JX));      NRT=0
  allocate(NINR(JRS,JP,JY,JX));  NINR=1  !set to one to avoid numerical failure
  allocate(IDTHR(JP,JY,JX));    IDTHR=0
  allocate(NIXBotRootLayer(JP,JY,JX));      NIXBotRootLayer=0
  allocate(NI(JP,JY,JX));       NI=0
  allocate(DMRT(JP,JY,JX));     DMRT=0._r8
  allocate(PR(JP,JY,JX));       PR=0._r8
  allocate(CWSRT(JP,JY,JX));    CWSRT=0._r8
  allocate(DMVL(jroots,JP,JY,JX));   DMVL=0._r8
  allocate(MaxPrimRootRadius1(jroots,JP,JY,JX)); MaxPrimRootRadius1=0._r8
  allocate(MaxSecndRootRadius1(jroots,JP,JY,JX)); MaxSecndRootRadius1=0._r8
  allocate(PrimRootXSecArea(jroots,JP,JY,JX)); PrimRootXSecArea=0._r8
  allocate(SecndRootXSecArea(jroots,JP,JY,JX)); SecndRootXSecArea=0._r8
  allocate(fTgrowRootP(JZ,JP,JY,JX));  fTgrowRootP=0._r8
  allocate(CNRT(JP,JY,JX));     CNRT=0._r8
  allocate(CPRT(JP,JY,JX));     CPRT=0._r8
  allocate(RootPorosity(jroots,JP,JY,JX));   RootPorosity=0._r8
  allocate(RSRR(jroots,JP,JY,JX));   RSRR=0._r8
  allocate(RSRA(jroots,JP,JY,JX));   RSRA=0._r8
  allocate(PTSHT(JP,JY,JX));    PTSHT=0._r8
  allocate(UPMXZH(jroots,JP,JY,JX)); UPMXZH=0._r8
  allocate(UPKMZH(jroots,JP,JY,JX)); UPKMZH=0._r8
  allocate(UPMNZH(jroots,JP,JY,JX)); UPMNZH=0._r8
  allocate(UPMXZO(jroots,JP,JY,JX)); UPMXZO=0._r8
  allocate(UPKMZO(jroots,JP,JY,JX)); UPKMZO=0._r8
  allocate(UPMNZO(jroots,JP,JY,JX)); UPMNZO=0._r8
  allocate(UPMXPO(jroots,JP,JY,JX)); UPMXPO=0._r8
  allocate(UPKMPO(jroots,JP,JY,JX)); UPKMPO=0._r8
  allocate(UPMNPO(jroots,JP,JY,JX)); UPMNPO=0._r8
  allocate(RRADP(jroots,JP,JY,JX));  RRADP=0._r8
  allocate(CNRTS(JP,JY,JX));    CNRTS=0._r8
  allocate(CPRTS(JP,JY,JX));    CPRTS=0._r8
  allocate(MaxPrimRootRadius(jroots,JP,JY,JX)); MaxPrimRootRadius=0._r8
  allocate(MaxSecndRootRadius(jroots,JP,JY,JX)); MaxSecndRootRadius=0._r8
  allocate(RTFQ(JP,JY,JX));     RTFQ=0._r8
  allocate(PORTX(jroots,JP,JY,JX));  PORTX=0._r8
  allocate(EPOOLN(npelms,JZ,JP,JY,JX));EPOOLN=0._r8
  allocate(RootLenPerP(jroots,JZ,JP,JY,JX));RootLenPerP=0._r8
  allocate(PrimRootLen(jroots,JZ,JC,JP,JY,JX));PrimRootLen=0._r8
  allocate(SecndRootLen(jroots,JZ,JC,JP,JY,JX));SecndRootLen=0._r8
  allocate(RootLenDensNLP(jroots,JZ,JP,JY,JX));RootLenDensNLP=0._r8
  allocate(PrimRootXNumL(jroots,JZ,JP,JY,JX));PrimRootXNumL=0._r8
  allocate(SecndRootXNumL(jroots,JZ,JP,JY,JX));SecndRootXNumL=0._r8
  allocate(RTN2(jroots,JZ,JC,JP,JY,JX));RTN2=0._r8
  allocate(AveSecndRootLen(jroots,JZ,JP,JY,JX));AveSecndRootLen=0._r8
  allocate(RTARP(jroots,JZ,JP,JY,JX));RTARP=0._r8
  allocate(RTVLW(jroots,JZ,JP,JY,JX));RTVLW=0._r8
  allocate(PrimRootRadius(jroots,JZ,JP,JY,JX));PrimRootRadius=0._r8
  allocate(RTVLP(jroots,JZ,JP,JY,JX));RTVLP=0._r8
  allocate(PrimRootDepth(jroots,JC,JP,JY,JX));PrimRootDepth=0._r8
  allocate(SecndRootRadius(jroots,JZ,JP,JY,JX));SecndRootRadius=0._r8
  allocate(PrimRootSpecLen(jroots,JP,JY,JX)); PrimRootSpecLen=0._r8
  allocate(SecndRootSpecLen(jroots,JP,JY,JX)); SecndRootSpecLen=0._r8
  allocate(PopPlantRootH2OUptake_vr(jroots,JZ,JP,JY,JX));PopPlantRootH2OUptake_vr=0._r8
  allocate(PSIRoot(jroots,JZ,JP,JY,JX));PSIRoot=0._r8
  allocate(PSIRootOSMO(jroots,JZ,JP,JY,JX));PSIRootOSMO=0._r8
  allocate(PSIRootTurg(jroots,JZ,JP,JY,JX));PSIRootTurg=0._r8
  allocate(trcg_rootml(idg_beg:idg_end-1,2,JZ,JP,JY,JX)); trcg_rootml =0._r8
  allocate(trcs_rootml(idg_beg:idg_end-1,2,JZ,JP,JY,JX)); trcs_rootml =0._r8
  allocate(TRootGasLoss_disturb(idg_beg:idg_end-1,JY,JX));TRootGasLoss_disturb=0._r8
  allocate(WTRTA(JP,JY,JX));    WTRTA=0._r8
  allocate(WTRTE(npelms,JP,JY,JX)); WTRTE=0._r8
  allocate(WTRTSE(npelms,JP,JY,JX));   WTRTSE=0._r8
  allocate(WSRTL(jroots,JZ,JP,JY,JX));WSRTL=0._r8
  allocate(WTRT1E(npelms,jroots,JZ,JC,JP,JY,JX));WTRT1E=0._r8
  allocate(WTRT2E(npelms,jroots,JZ,JC,JP,JY,JX));WTRT2E=0._r8
  allocate(RootCPZR(jroots,JZ,JP,JY,JX));RootCPZR=0._r8
  allocate(WTNDLE(npelms,JZ,JP,JY,JX)); WTNDLE=0._r8
  allocate(WTNDE(npelms,JP,JY,JX));  WTNDE=0._r8
  allocate(WTRTL(jroots,JZ,JP,JY,JX));WTRTL=0._r8
  allocate(EPOOLR(npelms,jroots,JZ,JP,JY,JX));EPOOLR=0._r8
  allocate(CEPOLR(npelms,jroots,JZ,JP,JY,JX));CEPOLR=0._r8
  allocate(RTWT1E(npelms,jroots,JRS,JP,JY,JX));RTWT1E=0._r8
  allocate(CWSRTL(jroots,JZ,JP,JY,JX));CWSRTL=0._r8
  end subroutine InitRootData

!----------------------------------------------------------------------
  subroutine DestructRootData
  use abortutils, only : destroy
  implicit none
  call destroy(NRT)
  call destroy(NINR)
  call destroy(IDTHR)
  call destroy(NIXBotRootLayer)
  call destroy(NI)
  call destroy(DMRT)
  call destroy(PR)
  call destroy(CWSRT)
  call destroy(DMVL)
  call destroy(MaxPrimRootRadius1)
  call destroy(MaxSecndRootRadius1)
  call destroy(PrimRootXSecArea)
  call destroy(SecndRootXSecArea)
  call destroy(fTgrowRootP)
  call destroy(CNRT)
  call destroy(CPRT)
  call destroy(RootPorosity)
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
  call destroy(MaxPrimRootRadius)
  call destroy(MaxSecndRootRadius)
  call destroy(RTFQ)
  call destroy(PORTX)
  call destroy(EPOOLN)
  call destroy(RootLenPerP)
  call destroy(PrimRootLen)
  call destroy(SecndRootLen)
  call destroy(RootLenDensNLP)
  call destroy(PrimRootXNumL)
  call destroy(SecndRootXNumL)
  call destroy(RTN2)
  call destroy(AveSecndRootLen)
  call destroy(RTARP)
  call destroy(RTVLW)
  call destroy(PrimRootRadius)
  call destroy(RTVLP)
  call destroy(PrimRootDepth)
  call destroy(SecndRootRadius)
  call destroy(PrimRootSpecLen)
  call destroy(SecndRootSpecLen)
  call destroy(PopPlantRootH2OUptake_vr)
  call destroy(PSIRoot)
  call destroy(PSIRootOSMO)
  call destroy(PSIRootTurg)
  call destroy(trcg_rootml)
  call destroy(trcs_rootml)
  call destroy(TRootGasLoss_disturb)
  call destroy(WTRTA)
  call destroy(WTRTE)
  call destroy(WTRTSE)
  call destroy(WSRTL)
  call destroy(WTRT1E)
  call destroy(WTRT2E)
  call destroy(RootCPZR)
  call destroy(WTNDLE)
  call destroy(WTNDE)
  call destroy(WTRTL)
  call destroy(EPOOLR)
  call destroy(CEPOLR)
  call destroy(RTWT1E)
  call destroy(CWSRTL)
  end subroutine DestructRootData

end module RootDataType
