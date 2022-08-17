module HfuncsMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use PlantAPIData
  use StartqsMod, only : startqs
  implicit none

  private


  character(len=*), parameter :: mod_filename = __FILE__

  real(r8) :: RLA,STK,TKCO,TFNP
  real(r8) :: WFNG
  integer :: KVSTGX,NB,NBTX,N

!
! PSILM=minimum canopy turgor potential for leaf expansion (MPa)
! PSILX=minimum canopy water potential for leafout of drought-deciduous PFT (MPa)
! PSILY=minimum canopy water potential for leafoff of drought-deciduous PFT (MPa
! GSTGG,GSTGR=normalized growth stage durations for vegetative,reproductive phenology
! NBX=maximum branch number for PFT defined by IBTYP in PFT file
! VRNE=maximum hours for leafout,leafoff
!
  real(r8), PARAMETER :: PSILM=0.1_r8
  real(r8), parameter :: PSILX=-0.2_r8
  real(r8) ,PARAMETER :: GSTGG=2.00_r8
  real(r8), PARAMETER :: GSTGR=0.667_r8
  real(r8), PARAMETER :: VRNE=3600.0_r8
  real(r8), parameter :: PSILY(0:3)=real((/-200.0,-2.0,-2.0,-2.0/),r8)
  integer, parameter :: NBX(0:3)=(/5,1,1,1/)

  public :: hfuncs
  contains

  subroutine hfuncs(I,J)
!
!     THIS subroutine CALCULATES PLANT PHENOLOGY
!
  implicit none

  integer, intent(in) :: I,J
  INTEGER :: NZ

! begin_execution
  associate(                          &
    NBRs1   =>  plt_morph%NBRs1     , &
    NB1s1   =>  plt_morph%NB1s1       &
  )
  DO 9985 NZ=1,NPs1

    IF(DATAPs1(NZ).NE.'NO')THEN
!
!         PPT=total biome population
!
      PPTs1=PPTs1+PPs1(NZ)
!
!         SET CROP FLAG ACCORDING TO PLANTING, HARVEST DATES, DEATH,
!         1 = ALIVE, 0 = NOT ALIVE
!
!         IDAY0,IYR0,IDAYH,IYRH=day,year of planting,harvesting
!         IYRC=current year
!         IDATAs1(3)=start year of current scenario
!         IDTH=PFT flag:0=alive,1=dead
!         IFLGC=PFT flag:0=not active,1=active
!         IFLGT=number of active PFT
!         DATAP=PFT file name
!
      call set_flags(I,J,NZ)
!
!         INITIALIZE VARIABLES IN ACTIVE PFT
!
      IF(IFLGCs1(NZ).EQ.1)THEN

        call stage_phenology_vars(I,J,NZ)

        call root_shoot_branching(I,J,NZ)
!
!           THE REST OF THE subroutine MODELS THE PHENOLOGY OF EACH BRANCH
!
!           IFLGA,IFLGE=flags for initializing leafout,leafoff
!           VRNS=leafout hours
!
        IF(IDAYs1(1,NB1s1(NZ),NZ).NE.0.OR.IFLGIs1(NZ).EQ.1)THEN
          DO 2010 NB=1,NBRs1(NZ)
            IF(IDTHBs1(NB,NZ).EQ.0)THEN
              call living_branch_phenology(I,J,NB,nz)
            ENDIF
!
!               KVSTG=integer of most recent leaf number currently growing
!               IFLGP=flag for remobilization
!
            KVSTGX=KVSTGs1(NB,NZ)
            IF(VSTGXs1(NB,NZ).LE.1.0E-06)THEN
              KVSTGs1(NB,NZ)=INT(VSTGs1(NB,NZ))+1
            ELSE
              KVSTGs1(NB,NZ)=INT(AMIN1(VSTGs1(NB,NZ),VSTGXs1(NB,NZ)))+1
            ENDIF
            KLEAFs1(NB,NZ)=MIN(24,KVSTGs1(NB,NZ))
            IF(KVSTGs1(NB,NZ).GT.KVSTGX)THEN
              IFLGPs1(NB,NZ)=1
            ELSE
              IFLGPs1(NB,NZ)=0
            ENDIF
!
!               PHENOLOGY
!
!               DYLX,DLYN=daylength of previous,current day
!               VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!
            IF(IDTHBs1(NB,NZ).EQ.0.OR.IFLGIs1(NZ).EQ.1)THEN
              IF(DYLNs1.GE.DYLXs1)THEN
                VRNYs1(NB,NZ)=VRNYs1(NB,NZ)+1.0
                VRNZs1(NB,NZ)=0.0
              ELSE
                VRNYs1(NB,NZ)=0.0
                VRNZs1(NB,NZ)=VRNZs1(NB,NZ)+1.0
              ENDIF

              call pft_specific_phenology(I,J,NZ)

            ENDIF
2010      CONTINUE
!
!             WATER STRESS INDICATOR
!
!             PSILT=canopy total water potential
!             PSILY=minimum canopy water potential for leafoff
!             WSTR=number of hours PSILT < PSILY (for output only)
!
            IF(PSILTs1(NZ).LT.PSILY(IGTYPs1(NZ)))THEN
              WSTRs1(NZ)=WSTRs1(NZ)+1.0
            ENDIF
          ENDIF
        ENDIF
      ENDIF
9985  CONTINUE
  RETURN
  end associate
  END subroutine hfuncs
!------------------------------------------------------------------------------------------

  subroutine set_flags(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ

  INTEGER :: L

! begin_execution
  IF(J.EQ.1)THEN
    IF(IDAY0s1(NZ).LE.IDAYHs1(NZ).OR.IYR0s1(NZ).LT.IYRHs1(NZ))THEN
      IF(I.GE.IDAY0s1(NZ).OR.IDATAs1(3).GT.IYR0s1(NZ))THEN
        IF(I.GT.IDAYHs1(NZ).AND.IYRCs1.GE.IYRHs1(NZ).AND.IDTHs1(NZ).EQ.1)THEN
          IFLGCs1(NZ)=0
        ELSE
          IF(I.EQ.IDAY0s1(NZ).AND.IDATAs1(3).EQ.IYR0s1(NZ))THEN
            IFLGCs1(NZ)=0
            IDTHs1(NZ)=0
            CALL STARTQs(NZ,NZ)
            TNBPs1=TNBPs1+WTRVXs1(NZ)
          ENDIF

          IF(DATAPs1(NZ).NE.'NO'.AND.IDTHs1(NZ).EQ.0)IFLGCs1(NZ)=1
        ENDIF
      ELSE
        IFLGCs1(NZ)=0
      ENDIF
    ELSE
      IF((I.LT.IDAY0s1(NZ).AND.I.GT.IDAYHs1(NZ) &
        .AND.IYRCs1.GE.IYRHs1(NZ).AND.IDTHs1(NZ).EQ.1) &
        .OR.(I.LT.IDAY0s1(NZ).AND.IYR0s1(NZ) &
        .GT.IYRHs1(NZ)))THEN
        IFLGCs1(NZ)=0
      ELSE
        IF(I.EQ.IDAY0s1(NZ).AND.IDATAs1(3).EQ.IYR0s1(NZ))THEN
          IFLGCs1(NZ)=0
          IDTHs1(NZ)=0
          CALL STARTQs(NZ,NZ)
          TNBPs1=TNBPs1+WTRVXs1(NZ)
        ENDIF
        IF(DATAPs1(NZ).NE.'NO'.AND.IDTHs1(NZ).EQ.0)IFLGCs1(NZ)=1
      ENDIF
    ENDIF
    IFLGTs1=IFLGTs1+IFLGCs1(NZ)
  ENDIF
  end subroutine set_flags
!------------------------------------------------------------------------------------------

  subroutine root_shoot_branching(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ

! begin_execution
  associate(                            &
    NRTs1     =>   plt_morph%NRTs1    , &
    NB1s1     =>   plt_morph%NB1s1    , &
    NBRs1     =>   plt_morph%NBRs1    , &
    NNODs1    =>   plt_morph%NNODs1   , &
    NBTs1     =>   plt_morph%NBTs1    , &
    NBTBs1    =>   plt_morph%NBTBs1   , &
    ARLFLs1   =>   plt_morph%ARLFLs1  , &
    PSTGs1    =>   plt_morph%PSTGs1     &
  )

!
! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
! CONCENTRATION PERMIT
!
! IFLGI=PFT initialization flag:0=no,1=yes
! PSIRG=root turgor potential
! ISTYP=growth habit from PFT file
! IDAYs1(2,=floral initiation date
! NBR=primary root axis number
! WTRVC=nonstructural C storage
! PB=nonstructural C concentration needed for branching
! IDTHB=branch life flag:0=living,1=dead
! PSTG=node number
! FNOD=scales node number for perennial vegetation (e.g. trees)
! NNOD=number of concurrently growing nodes
! XTLI,GROUP=node number at planting,floral initiation
!


  IF(IFLGIs1(NZ).EQ.0)THEN
    IF(J.EQ.1.AND.PPs1(NZ).GT.0.0)THEN
      IF(PSIRGs1(1,NGs1(NZ),NZ).GT.PSILM)THEN
        IF(ISTYPs1(NZ).NE.0.OR.IDAYs1(2,NB1s1(NZ),NZ).EQ.0)THEN
          IF((NBRs1(NZ).EQ.0.AND.WTRVCs1(NZ).GT.0.0) &
            .OR.(CCPOLPs1(NZ).GT.PBs1(NZ).AND.PBs1(NZ).GT.0.0))THEN
            DO 120 NB=1,JC1
              IF(IDTHBs1(NB,NZ).EQ.1)THEN
                IF(NB.EQ.NB1s1(NZ) &
                  .OR.PSTGs1(NB1s1(NZ),NZ).GT.NBTs1(NZ) &
                  +NNODs1(NZ)/FNODs1(NZ)+XTLIs1(NZ))THEN
                  NBTs1(NZ)=NBTs1(NZ)+1
                  NBRs1(NZ)=MIN(NBX(IBTYPs1(NZ)),MAX(NB,NBRs1(NZ)))
                  NBTBs1(NB,NZ)=NBTs1(NZ)-1
                  IDTHPs1(NZ)=0
                  IDTHBs1(NB,NZ)=0
                  VRNSs1(NB,NZ)=0.0
                  IF(ISTYPs1(NZ).EQ.0)THEN
                    GROUPs1(NB,NZ)=AMAX1(0.0,GROUPIs1(NZ)-NBTBs1(NB,NZ))
                  ELSE
                    GROUPs1(NB,NZ)=GROUPIs1(NZ)
                  ENDIF
                  exit
                ENDIF
              ENDIF
120         CONTINUE
          ENDIF
        ENDIF
      ENDIF
!
!     ADD AXIS TO ROOT IF PLANT GROWTH STAGE, ROOT NON-STRUCTURAL C
!     CONCENTRATION PERMIT
!
!     PR=nonstructural C concentration needed for root branching
!
      IF(PSIRGs1(1,NGs1(NZ),NZ).GT.PSILM)THEN
        IF(NRTs1(NZ).EQ.0.OR.PSTGs1(NB1s1(NZ),NZ) &
          .GT.NRTs1(NZ)/FNODs1(NZ)+XTLIs1(NZ))THEN
          IF((NRTs1(NZ).EQ.0.AND.WTRVCs1(NZ).GT.0.0) &
            .OR.(CCPOLPs1(NZ).GT.PRs1(NZ).AND.PRs1(NZ).GT.0.0))THEN
            NRTs1(NZ)=MIN(JC1,NRTs1(NZ)+1)
            IDTHRs1(NZ)=0
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine root_shoot_branching
!------------------------------------------------------------------------------------------

  subroutine stage_phenology_vars(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NB,N,L
  real(r8):: ARLSP
  associate(                          &
    NB1s1   =>  plt_morph%NB1s1     , &
    ARLFPs1 =>  plt_morph%ARLFPs1   , &
    NIXs1   =>  plt_morph%NIXs1     , &
    NBRs1   =>  plt_morph%NBRs1     , &
    NBTBs1  =>  plt_morph%NBTBs1    , &
    SDPTHs1 =>  plt_morph%SDPTHs1   , &
    ARSTPs1 =>  plt_morph%ARSTPs1   , &
    NIs1    =>  plt_morph%NIs1        &
  )
  RCO2Zs1(NZ)=0.0_r8
  ROXYZs1(NZ)=0.0_r8
  RCH4Zs1(NZ)=0.0_r8
  RN2OZs1(NZ)=0.0_r8
  RNH3Zs1(NZ)=0.0_r8
  RH2GZs1(NZ)=0.0_r8
  CPOOLPs1(NZ)=0.0_r8
  ZPOOLPs1(NZ)=0.0_r8
  PPOOLPs1(NZ)=0.0_r8
  NIs1(NZ)=NIXs1(NZ)
  NGs1(NZ)=MIN(NIs1(NZ),MAX(NGs1(NZ),NUs1))
  NB1s1(NZ)=1
  NBTX=1.0E+06
!
! TOTAL PLANT NON-STRUCTURAL C, N, P
!
! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
! NB1=main branch number
!
  DO 140 NB=1,NBRs1(NZ)
    IF(IDTHBs1(NB,NZ).EQ.0)THEN
      CPOOLPs1(NZ)=CPOOLPs1(NZ)+CPOOLs1(NB,NZ)
      ZPOOLPs1(NZ)=ZPOOLPs1(NZ)+ZPOOLs1(NB,NZ)
      PPOOLPs1(NZ)=PPOOLPs1(NZ)+PPOOLs1(NB,NZ)
      CPOLNPs1(NZ)=CPOLNPs1(NZ)+CPOLNBs1(NB,NZ)
      ZPOLNPs1(NZ)=ZPOLNPs1(NZ)+ZPOLNBs1(NB,NZ)
      PPOLNPs1(NZ)=PPOLNPs1(NZ)+PPOLNBs1(NB,NZ)
      IF(NBTBs1(NB,NZ).LT.NBTX)THEN
        NB1s1(NZ)=NB
        NBTX=NBTBs1(NB,NZ)
      ENDIF
    ENDIF
140   CONTINUE
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN ROOT
!
! WTRTL=root mass1(g)
! CPOOLR,ZPOOLR,PPOOLR=non-structl C,N,P in root(1),myco(2)(g)
! CCPOLR,CZPOLR,CPPOLR=non-structl C,N,P concn in root(1),myco(2)(g g-1)
!
  DO 180 N=1,MYs1(NZ)
    DO 160 L=NUs1,NIs1(NZ)
      IF(WTRTLs1(N,L,NZ).GT.ZEROLs1(NZ))THEN
        CCPOLRs1(N,L,NZ)=AMAX1(0.0,CPOOLRs1(N,L,NZ)/WTRTLs1(N,L,NZ))
        CZPOLRs1(N,L,NZ)=AMAX1(0.0,ZPOOLRs1(N,L,NZ)/WTRTLs1(N,L,NZ))
        CPPOLRs1(N,L,NZ)=AMAX1(0.0,PPOOLRs1(N,L,NZ)/WTRTLs1(N,L,NZ))
!       CCPOLRs1(N,L,NZ)=AMIN1(1.0,CCPOLRs1(N,L,NZ))
      ELSE
        CCPOLRs1(N,L,NZ)=1.0
        CZPOLRs1(N,L,NZ)=1.0
        CPPOLRs1(N,L,NZ)=1.0
      ENDIF
160 CONTINUE
180 CONTINUE
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN SHOOT
!
! CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy(g g-1)
! CCPLNP=nonstructural C concentration in canopy nodules
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!
  IF(WTLSs1(NZ).GT.ZEROLs1(NZ))THEN
    CCPOLPs1(NZ)=AMAX1(0.0,AMIN1(1.0,CPOOLPs1(NZ)/WTLSs1(NZ)))
    CCPLNPs1(NZ)=AMAX1(0.0,AMIN1(1.0,CPOLNPs1(NZ)/WTLSs1(NZ)))
    CZPOLPs1(NZ)=AMAX1(0.0,AMIN1(1.0,ZPOOLPs1(NZ)/WTLSs1(NZ)))
    CPPOLPs1(NZ)=AMAX1(0.0,AMIN1(1.0,PPOOLPs1(NZ)/WTLSs1(NZ)))
  ELSE
    CCPOLPs1(NZ)=1.0
    CCPLNPs1(NZ)=1.0
    CZPOLPs1(NZ)=1.0
    CPPOLPs1(NZ)=1.0
  ENDIF
  DO 190 NB=1,NBRs1(NZ)
    IF(WTLSBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      CCPOLBs1(NB,NZ)=AMAX1(0.0,CPOOLs1(NB,NZ)/WTLSBs1(NB,NZ))
      CZPOLBs1(NB,NZ)=AMAX1(0.0,ZPOOLs1(NB,NZ)/WTLSBs1(NB,NZ))
      CPPOLBs1(NB,NZ)=AMAX1(0.0,PPOOLs1(NB,NZ)/WTLSBs1(NB,NZ))
    ELSE
      CCPOLBs1(NB,NZ)=1.0
      CZPOLBs1(NB,NZ)=1.0
      CPPOLBs1(NB,NZ)=1.0
    ENDIF
190 CONTINUE
!
! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
!
! IDAYs1(1,=emergence date
! ARLFP,ARSTP=leaf,stalk areas
! HTCTL=hypocotyledon height
! SDPTH=seeding depth
! RTDP1=primary root depth
! VHCPC,WTSHT,VOLWC=canopy heat capacity,mass,water content
!
  IF(IDAYs1(1,NB1s1(NZ),NZ).EQ.0)THEN
    ARLSP=ARLFPs1(NZ)+ARSTPs1(NZ)
    IF((HTCTLs1(NZ).GT.SDPTHs1(NZ)) &
      .AND.(ARLSP.GT.ZEROLs1(NZ)) &
      .AND.(RTDP1s1(1,1,NZ).GT.SDPTHs1(NZ)+1.0E-06))THEN
      IDAYs1(1,NB1s1(NZ),NZ)=I
      VHCPCs1(NZ)=cpw*(WTSHTs1(NZ)*10.0E-06+VOLWCs1(NZ))
    ENDIF
  ENDIF
  end associate
  end subroutine stage_phenology_vars
!------------------------------------------------------------------------------------------

  subroutine pft_specific_phenology(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

!
! CALCULATE EVERGREEN PHENOLOGY DURING LENGTHENING PHOTOPERIODS
!
! IWTYP=phenology type from PFT file
! DYLX,DLYN=daylength of previous,current day
! VRNS,VRNF=leafout,leafoff hours
! VRNY=hourly counter for lengthening photoperiods
! IFLGF=flag for enabling leafoff:0=enable,1=disable
! ALAT=latitude
!
  IF(IWTYPs1(NZ).EQ.0)THEN
    IF(DYLNs1.GE.DYLXs1)THEN
      VRNSs1(NB,NZ)=VRNYs1(NB,NZ)
      IF(VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ) &
        .OR.(ALATs1.GT.0.0.AND.I.EQ.173) &
        .OR.(ALATs1.LT.0.0.AND.I.EQ.355))THEN
        VRNFs1(NB,NZ)=0.0
        IFLGFs1(NB,NZ)=0
      ENDIF
    ENDIF
!
!   CALCULATE EVERGREEN PHENOLOGY DURING SHORTENING PHOTOPERIODS
!
!   VRNS,VRNF=leafout,leafoff hours
!   VRNZ=hourly counter for shortening photoperiods
!   IFLGE=flag for enabling leafout:0=enable,1=disable
!   ALAT=latitude
!
    IF(DYLNs1.LT.DYLXs1)THEN
      VRNFs1(NB,NZ)=VRNZs1(NB,NZ)
      IF(VRNFs1(NB,NZ).GE.VRNXs1(NB,NZ) &
        .OR.(ALATs1.GT.0.0.AND.I.EQ.355) &
        .OR.(ALATs1.LT.0.0.AND.I.EQ.173))THEN
        VRNSs1(NB,NZ)=0.0
        IFLGEs1(NB,NZ)=0
      ENDIF
    ENDIF
!
!   CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATING HOURS ABOVE
!   SPECIFIED TEMPERATURE DURING LENGTHENING PHOTOPERIODS
!
!   IWTYP=phenology type from PFT file
!   DYLX,DLYN=daylength of previous,current day
!   VRNS,VRNL=leafout hours,hours required for leafout
!   VRNF,VRNX=leafoff hours,hours required for leafoff
!   IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!   TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!   ALAT=latitude
!   IDAYs1(2,=date of floral initiation
!
  ELSEIF(IWTYPs1(NZ).EQ.1)THEN
    IF((DYLNs1.GE.DYLXs1.OR.(DYLNs1.LT.DYLXs1.AND.VRNFs1(NB,NZ).LT.VRNXs1(NB,NZ))) &
      .AND.IFLGEs1(NB,NZ).EQ.0)THEN
        IF(TCGs1(NZ).GE.TCZs1(NZ))THEN
          VRNSs1(NB,NZ)=VRNSs1(NB,NZ)+1.0
        ENDIF
        IF(VRNSs1(NB,NZ).LT.VRNLs1(NB,NZ))THEN
          IF(TCGs1(NZ).LT.CTCs1(NZ))THEN
            VRNSs1(NB,NZ)=AMAX1(0.0,VRNSs1(NB,NZ)-1.0)
          ENDIF
        ENDIF
        IF(VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ) &
          .OR.(ALATs1.GT.0.0.AND.I.EQ.173) &
          .OR.(ALATs1.LT.0.0.AND.I.EQ.355))THEN
          VRNFs1(NB,NZ)=0.0
        ENDIF
      ENDIF
      IF(IDAYs1(2,NB,NZ).NE.0.OR.(DYLNs1.LT.DYLXs1.AND.DYLNs1.LT.12.0))THEN
        IFLGFs1(NB,NZ)=0
      ENDIF
!
!     CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATING HOURS BELOW
!     SPECIFIED TEMPERATURE DURING SHORTENING PHOTOPERIODS
!
!     DYLX,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     IDAYs1(2,=date of floral initiation
!
      IF(DYLNs1.LT.DYLXs1.AND.IFLGFs1(NB,NZ).EQ.0 &
        .AND.IDAYs1(2,NB,NZ).NE.0)THEN
        IF(TCGs1(NZ).LE.TCXs1(NZ))THEN
          VRNFs1(NB,NZ)=VRNFs1(NB,NZ)+1.0
        ENDIF
        IF(VRNFs1(NB,NZ).GE.VRNXs1(NB,NZ).AND.IFLGEs1(NB,NZ).EQ.1)THEN
          VRNSs1(NB,NZ)=0.0
          IFLGEs1(NB,NZ)=0
        ENDIF
      ENDIF

!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING HOURS
!     ABOVE SPECIFIED WATER POTENTIAL DURING DORMANCY
!
!     IWTYP=phenology type from PFT file
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF=leafoff hours
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSILT=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     IDAYs1(2,=date of floral initiation
!
    ELSEIF(IWTYPs1(NZ).EQ.2.OR.IWTYPs1(NZ).EQ.4 &
      .OR.IWTYPs1(NZ).EQ.5)THEN
      IF(IFLGEs1(NB,NZ).EQ.0)THEN
        IF(PSILTs1(NZ).GE.PSILX)THEN
          VRNSs1(NB,NZ)=VRNSs1(NB,NZ)+1.0
        ENDIF
        IF(VRNSs1(NB,NZ).LT.VRNLs1(NB,NZ))THEN
          IF(PSILTs1(NZ).LT.PSILY(IGTYPs1(NZ)))THEN
            VRNSs1(NB,NZ)=AMAX1(0.0,VRNSs1(NB,NZ)-12.0)
          ENDIF
        ENDIF
        IF(VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ))THEN
          VRNFs1(NB,NZ)=0.0
          IF(IDAYs1(2,NB,NZ).NE.0)IFLGFs1(NB,NZ)=0
        ENDIF
      ENDIF
      IF(IDAYs1(2,NB,NZ).NE.0)IFLGFs1(NB,NZ)=0
!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING HOURS
!     BELOW SPECIFIED WATER POTENTIAL DURING GROWING SEASON
!
!     VRNS=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSILT=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     IDAYs1(2,=date of floral initiation
!     VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!     VRNE=maximum hours for leafout,leafoff
!
      IF(IFLGEs1(NB,NZ).EQ.1 &
        .AND.IFLGFs1(NB,NZ).EQ.0)THEN
        IF(PSILTs1(NZ).LT.PSILY(IGTYPs1(NZ)))THEN
          VRNFs1(NB,NZ)=VRNFs1(NB,NZ)+1.0
        ENDIF
        IF(IWTYPs1(NZ).EQ.4)THEN
          IF(VRNZs1(NB,NZ).GT.VRNE)THEN
            VRNFs1(NB,NZ)=VRNZs1(NB,NZ)
          ENDIF
        ELSEIF(IWTYPs1(NZ).EQ.5)THEN
          IF(VRNYs1(NB,NZ).GT.VRNE)THEN
            VRNFs1(NB,NZ)=VRNYs1(NB,NZ)
          ENDIF
        ENDIF
        IF(VRNFs1(NB,NZ).GE.VRNXs1(NB,NZ) &
          .AND.IFLGEs1(NB,NZ).EQ.1)THEN
          VRNSs1(NB,NZ)=0.0
          IFLGEs1(NB,NZ)=0
        ENDIF
      ENDIF
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS ABOVE SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     LENGTHENING PHOTOPERIODS
!
!     IWTYP=phenology type from PFT file
!     DYLX,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     PSILT=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     IDAYs1(2,=date of floral initiation
!
    ELSEIF(IWTYPs1(NZ).EQ.3)THEN
      IF((DYLNs1.GE.DYLXs1.OR.DYLNs1.GE.DYLMs1-2.0_r8).AND.IFLGEs1(NB,NZ).EQ.0)THEN
        IF(TCGs1(NZ).GE.TCZs1(NZ).AND.PSILGs1(NZ).GT.PSILM)THEN
          VRNSs1(NB,NZ)=VRNSs1(NB,NZ)+1.0
        ENDIF
        IF(VRNSs1(NB,NZ).LT.VRNLs1(NB,NZ))THEN
          IF(TCGs1(NZ).LT.CTCs1(NZ) &
            .OR.PSILGs1(NZ).LT.PSILM)THEN
            VRNSs1(NB,NZ)=AMAX1(0.0,VRNSs1(NB,NZ)-1.5)
          ENDIF
        ENDIF
        IF(VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ))THEN
          VRNFs1(NB,NZ)=0.0
          IF(IDAYs1(2,NB,NZ).NE.0)IFLGFs1(NB,NZ)=0
        ENDIF
      ENDIF
      IF(IDAYs1(2,NB,NZ).NE.0)IFLGFs1(NB,NZ)=0
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS BELOW SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     SHORTENING PHOTOPERIODS
!
!     DYLX,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     PSILT=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     IDAYs1(2,=date of floral initiation
!
      IF((DYLNs1.LT.DYLXs1.OR.DYLNs1.LT.24.0-DYLMs1+2.0_r8).AND.IFLGFs1(NB,NZ).EQ.0)THEN
        IF(TCGs1(NZ).LE.TCXs1(NZ).OR.PSILTs1(NZ).LT.PSILY(IGTYPs1(NZ)))THEN
          VRNFs1(NB,NZ)=VRNFs1(NB,NZ)+1.0
        ENDIF
        IF(VRNFs1(NB,NZ).GE.VRNXs1(NB,NZ) &
          .AND.IFLGEs1(NB,NZ).EQ.1)THEN
        VRNSs1(NB,NZ)=0.0
        IFLGEs1(NB,NZ)=0
      ENDIF
    ENDIF
  ENDIF
  end subroutine pft_specific_phenology
!------------------------------------------------------------------------------------------

  subroutine living_branch_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J,NB,NZ

  real(r8) :: ACTV,OFNG
  real(r8) :: PPD
  real(r8) :: RTK
  real(r8) :: RNI

! begin_execution
  associate(                            &
    PSTGs1    =>   plt_morph%PSTGs1   , &
    PSTGIs1   =>   plt_morph%PSTGIs1  , &
    PSTGFs1   =>   plt_morph%PSTGFs1  , &
    NB1s1     =>   plt_morph%NB1s1      &
  )
  IF(IDAYs1(1,NB,NZ).EQ.0)THEN
    IDAYs1(1,NB,NZ)=I
    IFLGAs1(NB,NZ)=1
    IFLGEs1(NB,NZ)=0
    VRNSs1(NB,NZ)=0.5_r8*VRNSs1(NB1s1(NZ),NZ)
  ENDIF
!
! CALCULATE NODE INITIATION AND LEAF APPEARANCE RATES
! FROM TEMPERATURE FUNCTION CALCULATED IN 'UPTAKE' AND
! RATES AT 25C ENTERED IN 'READQ' EXCEPT WHEN DORMANT
!
! IWTYP=phenology type from PFT file
! VRNF,VRNX=leafoff hours,hours required for leafoff
! TKG,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
! OFFST=shift in Arrhenius curve for thermal adaptation
! TFNP=temperature function for phenology (25 oC =1 )
! 8.313,710.0=gas constant,enthalpy
! 60000,197500,218500=energy of activn,high,low temp inactivn(KJ mol-1)
! RNI,RLA=rates of node initiation,leaf appearance
! XRNI,XRLA=rate of node initiation,leaf appearance at 25 oC (h-1)
!
  IF(IWTYPs1(NZ).EQ.0.OR.VRNFs1(NB,NZ).LT.VRNXs1(NB,NZ))THEN
    TKCO=TKGs1(NZ)+OFFSTs1(NZ)
    RTK=8.3143_r8*TKCO
    STK=710.0_r8*TKCO
    ACTV=1+EXP((197500_r8-STK)/RTK)+EXP((STK-218500._r8)/RTK)
    TFNP=EXP(24.269_r8-60000._r8/RTK)/ACTV
    RNI=AMAX1(0.0_r8,XRNIs1(NZ)*TFNP)
    RLA=AMAX1(0.0_r8,XRLAs1(NZ)*TFNP)
!
!   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
!
!   PSILG=leaf turgor potential
!   WFNG=water stress effect on phenology
!
    IF(ISTYPs1(NZ).EQ.0)THEN
      IF(IDAYs1(6,NB,NZ).EQ.0)THEN
        WFNG=EXP(0.025_r8*PSILTs1(NZ))
        RNI=RNI*WFNG
        RLA=RLA*WFNG
      ENDIF
      IF(IDAYs1(2,NB,NZ).EQ.0)THEN
        OFNG=SQRT(OSTRs1(NZ))
        RNI=RNI*OFNG
        RLA=RLA*OFNG
      ENDIF
    ENDIF
!
!   ACCUMULATE NODE INITIATION AND LEAF APPEARANCE RATES
!   INTO TOTAL NUMBER OF NODES AND LEAVES
!
!   PSTG,VSTG=number of nodes initiated,leaves appeared
!
    PSTGs1(NB,NZ)=PSTGs1(NB,NZ)+RNI
    VSTGs1(NB,NZ)=VSTGs1(NB,NZ)+RLA
!
!   USE TOTAL NUMBER OF NODES TO CALCULATE PROGRESSION THROUGH
!   VEGETATIVE AND REPRODUCTIVE GROWTH STAGES. THIS PROGRESSION
!   IS USED TO SET START AND END DATES FOR GROWTH STAGES BELOW
!
!   GSTGI=vegetative node number normalized for maturity group
!   PSTGI=node number at floral initiation
!   GROUPI=node number required for floral initiation
!   DGSTGI,TGSTGI=hourly,total change in GSTGI
!   GSTGF=reproductive node number normalized for maturity group
!   PSTGF=node number at flowering
!   DGSTGF,TGSTGF=hourly,total change in GSTGF
!   IFLGG=PFT senescence flag
!
    IF(IDAYs1(2,NB,NZ).NE.0)THEN
      GSTGIs1(NB,NZ)=(PSTGs1(NB,NZ)-PSTGIs1(NB,NZ))/GROUPIs1(NZ)
      DGSTGIs1(NB,NZ)=RNI/(GROUPIs1(NZ)*GSTGG)
      TGSTGIs1(NB,NZ)=TGSTGIs1(NB,NZ)+DGSTGIs1(NB,NZ)
    ENDIF
    IF(IDAYs1(6,NB,NZ).NE.0)THEN
      GSTGFs1(NB,NZ)=(PSTGs1(NB,NZ)-PSTGFs1(NB,NZ))/GROUPIs1(NZ)
      DGSTGFs1(NB,NZ)=RNI/(GROUPIs1(NZ)*GSTGR)
      TGSTGFs1(NB,NZ)=TGSTGFs1(NB,NZ)+DGSTGFs1(NB,NZ)
    ENDIF
    IFLGGs1(NB,NZ)=1
  ELSE
    IFLGGs1(NB,NZ)=0
  ENDIF
!
! REPRODUCTIVE GROWTH STAGES ADVANCE WHEN THRESHOLD NUMBER
! OF NODES HAVE BEEN INITIATED. FIRST DETERMINE PHOTOPERIOD
! AND TEMPERATURE EFFECTS ON FINAL VEG NODE NUMBER FROM
! NUMBER OF INITIATED NODES
!
! PSTG=number of nodes initiated
! PSTGI=node number at floral initiation
! GROUP=node number required for floral initiation
! VRNS,VRNL=leafout hours,hours required for leafout
! DYLX,DLYN=daylength of previous,current day
! ISTYP=growth habit from PFT file
! IWTYP=phenology type from PFT file
! ZC,DPTHS=canopy height,snowpack depth
!
  IF(IDAYs1(2,NB,NZ).EQ.0)THEN
    IF(PSTGs1(NB,NZ).GT.GROUPs1(NB,NZ)+PSTGIs1(NB,NZ) &
      .AND.((VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ)) &
      .OR.(I.GE.IDAY0s1(NZ).AND.IYRCs1.EQ.IYR0s1(NZ) &
      .AND.DYLNs1.GT.DYLXs1)) &
      .OR.(((ISTYPs1(NZ).EQ.1.AND.(IWTYPs1(NZ).EQ.1 &
      .OR.IWTYPs1(NZ).EQ.3)) &
      .OR.(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).EQ.0)) &
      .AND.ZCs1(NZ).GE.DPTHSs1-ZEROs1 &
      .AND.DYLNs1.LT.DYLXs1))THEN
!
!     FINAL VEGETATIVE NODE NUMBER DEPENDS ON PHOTOPERIOD FROM 'DAY'
!     AND ON MATURITY GROUP, CRITICAL PHOTOPERIOD AND PHOTOPERIOD
!     SENSITIVITY ENTERED IN 'READQ'
!
!     IPTYP=photoperiod type from PFT file
!     PPD=photoperiod sensitivity
!     XDL=critical photoperiod from PFT file
!     IDAYs1(2,=date of floral initiation
!     VSTGX=node number on date of floral initiation
!
      IF(IPTYPs1(NZ).EQ.0)THEN
        PPD=0.0
      ELSE
        PPD=AMAX1(0.0,XDLs1(NZ)-DYLNs1)
        IF(IPTYPs1(NZ).EQ.1.AND.DYLNs1.GE.DYLXs1)PPD=0.0
      ENDIF

      IF(IPTYPs1(NZ).EQ.0 &
        .OR.(IPTYPs1(NZ).EQ.1.AND.PPD.GT.XPPDs1(NZ)) &
        .OR.(IPTYPs1(NZ).EQ.2.AND.PPD.LT.XPPDs1(NZ)) &
        .OR.(((ISTYPs1(NZ).EQ.1.AND.(IWTYPs1(NZ).EQ.1 &
        .OR.IWTYPs1(NZ).EQ.3)) &
        .OR.(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).EQ.0)) &
        .AND.ZCs1(NZ).GE.DPTHSs1-ZEROs1 &
        .AND.DYLNs1.LT.DYLXs1))THEN
        IDAYs1(2,NB,NZ)=I
        PSTGIs1(NB,NZ)=PSTGs1(NB,NZ)
        IF(ISTYPs1(NZ).EQ.0.AND.IDTYPs1(NZ).EQ.0)THEN
          VSTGXs1(NB,NZ)=PSTGs1(NB,NZ)
        ENDIF
      ENDIF
    ENDIF
!
!   STEM ELONGATION
!
!   GSTGI=vegetative node number normalized for maturity group
!   GSTGG=normalized growth stage durations for vegetative phenology
!   IDAYs1(3,=start of stem elongation and setting max seed number
!
  ELSEIF(IDAYs1(3,NB,NZ).EQ.0)THEN
    IF(GSTGIs1(NB,NZ).GT.0.25_r8*GSTGG &
      .OR.((IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3) &
      .AND.ISTYPs1(NZ).NE.0.AND.IPTYPs1(NZ).NE.1 &
      .AND.DYLNs1.LT.DYLXs1.AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ)) &
      .OR.(IWTYPs1(NZ).EQ.2.AND.ISTYPs1(NZ).EQ.0) &
      .AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ))THEN
      IDAYs1(3,NB,NZ)=I
    ENDIF
!
!   IDAYs1(4,=mid stem elongation
!
  ELSEIF(IDAYs1(4,NB,NZ).EQ.0)THEN
    IF(GSTGIs1(NB,NZ).GT.0.50*GSTGG &
      .OR.((IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3) &
      .AND.ISTYPs1(NZ).NE.0.AND.IPTYPs1(NZ).NE.1 &
      .AND.DYLNs1.LT.DYLXs1.AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ)) &
      .OR.(IWTYPs1(NZ).EQ.2.AND.ISTYPs1(NZ).EQ.0) &
      .AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ))THEN
      IDAYs1(4,NB,NZ)=I
      IF(ISTYPs1(NZ).EQ.0.AND.IDTYPs1(NZ).NE.0)THEN
        VSTGXs1(NB,NZ)=PSTGs1(NB,NZ)
      ENDIF
    ENDIF
!
!   IDAYs1(5,=end of stem elongation and setting max seed number
!
  ELSEIF(IDAYs1(5,NB,NZ).EQ.0)THEN
    IF(GSTGIs1(NB,NZ).GT.1.00*GSTGG &
      .OR.((IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3) &
      .AND.ISTYPs1(NZ).NE.0.AND.IPTYPs1(NZ).NE.1 &
      .AND.DYLNs1.LT.DYLXs1.AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ)) &
      .OR.(IWTYPs1(NZ).EQ.2.AND.ISTYPs1(NZ).EQ.0) &
      .AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ))THEN
      IDAYs1(5,NB,NZ)=I
    ENDIF
!
!   ANTHESIS OCCURS WHEN THE NUMBER OF LEAVES THAT HAVE APPEARED
!   EQUALS THE NUMBER OF NODES INITIATED WHEN THE FINAL VEGETATIVE
!   NODE NUMBER WAS SET ABOVE
!
!   IDAYs1(6,=start of anthesis and setting final seed number
!   VSTG=number of leaves appeared
!   PSTGI=node number at floral initiation
!   ISTYP,IWTYP,IPTYP=growth habit,phenology,photoperiod type from PFT file
!   IFLGE=flag for enabling leafout:0=enable,1=disable
!   VRNF,VRNX=leafoff hours,hours required for leafoff
!   DYLX,DLYN=daylength of previous,current day
!   PSTGF=number of nodes at anthesis
!
  ELSEIF(IDAYs1(6,NB,NZ).EQ.0)THEN
    IF((VSTGs1(NB,NZ).GT.PSTGIs1(NB,NZ)) &
      .OR.(ISTYPs1(NZ).NE.0.AND.IDAYs1(5,NB,NZ).NE.0) &
      .OR.((IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3) &
      .AND.ISTYPs1(NZ).NE.0.AND.IPTYPs1(NZ).NE.1 &
      .AND.DYLNs1.LT.DYLXs1.AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ)) &
      .OR.(IWTYPs1(NZ).EQ.2.AND.ISTYPs1(NZ).EQ.0) &
      .AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ))THEN
      IF(NB.EQ.NB1s1(NZ).OR.IDAYs1(6,NB1s1(NZ),NZ).NE.0)THEN
        IDAYs1(6,NB,NZ)=I
        PSTGFs1(NB,NZ)=PSTGs1(NB,NZ)
      ENDIF
    ENDIF
!
!   START GRAIN FILL PERIOD
!
!   IDAYs1(7,=start of grain filling and setting max seed size
!   GSTGF=reproductive node number normalized for maturity group
!   GSTGR=normalized growth stage durations for reproductive phenology
!
!
  ELSEIF(IDAYs1(7,NB,NZ).EQ.0)THEN
    IF(GSTGFs1(NB,NZ).GT.0.50*GSTGR &
      .OR.((IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3) &
      .AND.ISTYPs1(NZ).NE.0.AND.IPTYPs1(NZ).NE.1 &
      .AND.DYLNs1.LT.DYLXs1.AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ)) &
      .OR.(IWTYPs1(NZ).EQ.2.AND.ISTYPs1(NZ).EQ.0) &
      .AND.IFLGEs1(NB,NZ).EQ.1 &
      .AND.VRNFs1(NB,NZ).GT.VRNXs1(NB,NZ))THEN
        IDAYs1(7,NB,NZ)=I
    ENDIF
!
!   END SEED NUMBER SET PERIOD
!
!   IDAYs1(8,=end date setting for final seed number
!
  ELSEIF(IDAYs1(8,NB,NZ).EQ.0)THEN
    IF(GSTGFs1(NB,NZ).GT.1.00*GSTGR)THEN
      IDAYs1(8,NB,NZ)=I
    ENDIF
!
!   END SEED SIZE SET PERIOD
!
!   IDAYs1(9,=end of setting max seed size
!
  ELSEIF(IDAYs1(9,NB,NZ).EQ.0)THEN
    IF(GSTGFs1(NB,NZ).GT.1.50*GSTGR)THEN
      IDAYs1(9,NB,NZ)=I
    ENDIF
  ENDIF
  end associate
  end subroutine living_branch_phenology

end module HfuncsMod
