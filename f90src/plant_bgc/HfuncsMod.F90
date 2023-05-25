module HfuncsMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use PlantAPIData
  use minimathmod, only : AZMAX1
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
  integer , parameter :: NBX(0:3)=(/5,1,1,1/)

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
  associate(                           &
    PSICanP   =>  plt_ew%PSICanP     , &
    VRNZ    =>  plt_pheno%VRNZ   , &
    IFLGI   =>  plt_pheno%IFLGI  , &
    IFLGP   =>  plt_pheno%IFLGP  , &
    VSTGX   =>  plt_pheno%VSTGX  , &
    VRNY    =>  plt_pheno%VRNY   , &
    IFLGC   =>  plt_pheno%IFLGC  , &
    IDTHB   =>  plt_pheno%IDTHB  , &
    IDAY    =>  plt_pheno%IDAY   , &
    WSTR    =>  plt_pheno%WSTR   , &
    KVSTG   =>  plt_pheno%KVSTG  , &
    IGTYP   =>  plt_pheno%IGTYP  , &
    DYLN    =>  plt_site%DYLN    , &
    DATAP   =>  plt_site%DATAP   , &
    PPT     =>  plt_site%PPT     , &
    DYLX    =>  plt_site%DYLX    , &
    NP      =>  plt_site%NP      , &
    pftPlantPopulation      =>  plt_site%pftPlantPopulation      , &
    KLEAF   =>  plt_morph%KLEAF  , &
    VSTG    =>  plt_morph%VSTG   , &
    NBR     =>  plt_morph%NBR    , &
    NB1     =>  plt_morph%NB1      &
  )
  D9985: DO NZ=1,NP

    IF(DATAP(NZ).NE.'NO')THEN
!
!         PPT=total biome population
!
      PPT=PPT+pftPlantPopulation(NZ)
!
!         SET CROP FLAG ACCORDINGTopRootLayerTO PLANTING, HARVEST DATES, DEATH,
!         1 = ALIVE, 0 = NOT ALIVE
!
!         IDAY0,IYR0,IDAYH,IYRH=day,year of planting,harvesting
!         IYRC=current year
!         IDATA(3)=start year of current scenario
!         IDTH=PFT flag:0=alive,1=dead
!         IFLGC=PFT flag:0=not active,1=active
!         IFLGT=number of active PFT
!         DATAP=PFT file name
!
      call set_flags(I,J,NZ)
!
!         INITIALIZE VARIABLES IN ACTIVE PFT
!
      IF(IFLGC(NZ).EQ.PlantIsActive)THEN

        call stage_phenology_vars(I,J,NZ)

        call root_shoot_branching(I,J,NZ)
!
!           THE REST OF THE subroutine MODELS THE PHENOLOGY OF EACH BRANCH
!
!           IFLGA,IFLGE=flags for initializing leafout,leafoff
!           VRNS=leafout hours
!
        IF(IDAY(1,NB1(NZ),NZ).NE.0.OR.IFLGI(NZ).EQ.itrue)THEN
          D2010: DO NB=1,NBR(NZ)
            IF(IDTHB(NB,NZ).EQ.iliving_branch)THEN
              call living_branch_phenology(I,J,NB,nz)
            ENDIF
!
!               KVSTG=integer of most recent leaf number currently growing
!               IFLGP=flag for remobilization
!
            KVSTGX=KVSTG(NB,NZ)
            IF(VSTGX(NB,NZ).LE.ppmc)THEN
              KVSTG(NB,NZ)=INT(VSTG(NB,NZ))+1
            ELSE
              KVSTG(NB,NZ)=INT(AMIN1(VSTG(NB,NZ),VSTGX(NB,NZ)))+1
            ENDIF
            KLEAF(NB,NZ)=MIN(24,KVSTG(NB,NZ))
            IF(KVSTG(NB,NZ).GT.KVSTGX)THEN
              IFLGP(NB,NZ)=1
            ELSE
              IFLGP(NB,NZ)=0
            ENDIF
!
!               PHENOLOGY
!
!               DYLX,DLYN=daylength of previous,current day
!               VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!
            IF(IDTHB(NB,NZ).EQ.iliving_branch.OR.IFLGI(NZ).EQ.itrue)THEN
              IF(DYLN.GE.DYLX)THEN
                VRNY(NB,NZ)=VRNY(NB,NZ)+1.0_r8
                VRNZ(NB,NZ)=0.0_r8
              ELSE
                VRNY(NB,NZ)=0.0_r8
                VRNZ(NB,NZ)=VRNZ(NB,NZ)+1.0_r8
              ENDIF

              call pft_specific_phenology(I,J,NZ)

            ENDIF
          ENDDO D2010
!
!             WATER STRESS INDICATOR
!
!             PSICanP=canopy total water potential
!             PSILY=minimum canopy water potential for leafoff
!             WSTR=number of hours PSICanP < PSILY (for output only)
!
          IF(PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
            WSTR(NZ)=WSTR(NZ)+1.0_r8
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO D9985
  RETURN
  end associate
  END subroutine hfuncs
!------------------------------------------------------------------------------------------

  subroutine set_flags(I,J,NZ)
  use EcoSIMCtrlDataType, only : iyear_cur

  implicit none
  integer, intent(in) :: I,J,NZ

  INTEGER :: L

! begin_execution
  associate(                        &
    IYR0    =>  plt_distb%IYR0    , &
    IYRH    =>  plt_distb%IYRH    , &  !year of harvest
    IDAY0   =>  plt_distb%IDAY0   , &  !day of planting
    IDAYH   =>  plt_distb%IDAYH   , &  !day of harvest
    DATAP   =>  plt_site%DATAP    , &
    IDATA   =>  plt_site%IDATA    , &
    IYRC    =>  plt_site%IYRC     , &
    IFLGT   =>  plt_site%IFLGT    , &
    TNBP    =>  plt_bgcr%TNBP     , &
    WTRVX   =>  plt_biom%WTRVX    , &
    IDTH    =>  plt_pheno%IDTH    , &
    IFLGC   =>  plt_pheno%IFLGC     &
  )
  IF(J.EQ.1)THEN
    IF(IDAY0(NZ).LE.IDAYH(NZ).OR.IYR0(NZ).LT.IYRH(NZ))THEN
      IF(I.GE.IDAY0(NZ).OR.iyear_cur.GT.IYR0(NZ))THEN
        IF(I.GT.IDAYH(NZ).AND.IYRC.GE.IYRH(NZ).AND.IDTH(NZ).EQ.1)THEN
          IFLGC(NZ)=ipltdorm
        ELSE
          IF(I.EQ.IDAY0(NZ).AND.iyear_cur.EQ.IYR0(NZ))THEN
            IFLGC(NZ)=ipltdorm
            IDTH(NZ)=0
            CALL STARTQs(NZ,NZ)
            TNBP=TNBP+WTRVX(NZ)
          ENDIF

          IF(DATAP(NZ).NE.'NO'.AND.IDTH(NZ).EQ.0)then
            IFLGC(NZ)=PlantIsActive
          endif
        ENDIF
      ELSE
        IFLGC(NZ)=ipltdorm
      ENDIF
    ELSE
      IF((I.LT.IDAY0(NZ).AND.I.GT.IDAYH(NZ) &
        .AND.IYRC.GE.IYRH(NZ).AND.IDTH(NZ).EQ.1) &
        .OR.(I.LT.IDAY0(NZ).AND.IYR0(NZ) &
        .GT.IYRH(NZ)))THEN
        IFLGC(NZ)=ipltdorm
      ELSE
        IF(I.EQ.IDAY0(NZ).AND.iyear_cur.EQ.IYR0(NZ))THEN
          IFLGC(NZ)=ipltdorm
          IDTH(NZ)=0
          CALL STARTQs(NZ,NZ)
          TNBP=TNBP+WTRVX(NZ)
        ENDIF
        IF(DATAP(NZ).NE.'NO'.AND.IDTH(NZ).EQ.0)then
          IFLGC(NZ)=PlantIsActive
        endif
      ENDIF
    ENDIF
    IFLGT=IFLGT+IFLGC(NZ)
  ENDIF
  end associate
  end subroutine set_flags
!------------------------------------------------------------------------------------------

  subroutine root_shoot_branching(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ

! begin_execution
  associate(                            &
    CEPOLP  =>   plt_biom%CEPOLP  , &
    WTRVE   =>   plt_biom%WTRVE   , &
    GROUP   =>   plt_pheno%GROUP  , &
    IDAY    =>   plt_pheno%IDAY   , &
    IFLGI   =>   plt_pheno%IFLGI  , &
    IDTHR   =>   plt_pheno%IDTHR  , &
    ISTYP   =>   plt_pheno%ISTYP  , &
    IDTHB   =>   plt_pheno%IDTHB  , &
    IDTHP   =>   plt_pheno%IDTHP  , &
    IBTYP   =>   plt_pheno%IBTYP  , &
    PR      =>   plt_pheno%PR     , &
    GROUPI  =>   plt_pheno%GROUPI , &
    PB      =>   plt_pheno%PB     , &
    pftPlantPopulation      =>   plt_site%pftPlantPopulation      , &
    VRNS    =>   plt_pheno%VRNS   , &
    PSIRootTurg   =>   plt_ew%PSIRootTurg     , &
    FNOD    =>   plt_allom%FNOD   , &
    NRT     =>   plt_morph%NRT    , &
    NB1     =>   plt_morph%NB1    , &
    NBR     =>   plt_morph%NBR    , &
    NNOD    =>   plt_morph%NNOD   , &
    NBT     =>   plt_morph%NBT    , &
    NBTB    =>   plt_morph%NBTB   , &
    NGTopRootLayer     =>   plt_morph%NGTopRootLayer    , &
    XTLI    =>   plt_morph%XTLI   , &
    PSTG    =>   plt_morph%PSTG     &
  )

!
! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
! CONCENTRATION PERMIT
!
! IFLGI=PFT initialization flag:0=no,1=yes
! PSIRootTurg=root turgor potential
! ISTYP=growth habit from PFT file
! IDAY(2,=floral initiation date
! NBR=primary root axis number
! WTRVC=nonstructural C storage
! PB=nonstructural C concentration needed for branching
! IDTHB=branch life flag:0=living,1=dead
! PSTG=node number
! FNOD=scales node number for perennial vegetation (e.g. trees)
! NNOD=number of concurrently growing nodes
! XTLI,GROUP=node number at planting,floral initiation
! IBTYP: setup for phenologically-driven above-ground turnover


  IF(IFLGI(NZ).EQ.0)THEN
    IF(J.EQ.1.AND.pftPlantPopulation(NZ).GT.0.0_r8)THEN
      IF(PSIRootTurg(ipltroot,NGTopRootLayer(NZ),NZ).GT.PSILM)THEN
        IF(ISTYP(NZ).NE.iplt_annual.OR.IDAY(2,NB1(NZ),NZ).EQ.0)THEN
          IF((NBR(NZ).EQ.0.AND.WTRVE(ielmc,NZ).GT.0.0_r8) &
            .OR.(CEPOLP(ielmc,NZ).GT.PB(NZ).AND.PB(NZ).GT.0.0_r8))THEN
            D120: DO NB=1,JC1
              IF(IDTHB(NB,NZ).EQ.ibrdead)THEN
                IF(NB.EQ.NB1(NZ).OR.PSTG(NB1(NZ),NZ).GT.NBT(NZ) &
                  +NNOD(NZ)/FNOD(NZ)+XTLI(NZ))THEN
                  NBT(NZ)=NBT(NZ)+1
                  NBR(NZ)=MIN(NBX(IBTYP(NZ)),MAX(NB,NBR(NZ)))
                  NBTB(NB,NZ)=NBT(NZ)-1
                  IDTHP(NZ)=ibralive
                  IDTHB(NB,NZ)=ibralive
                  VRNS(NB,NZ)=0.0_r8
                  IF(ISTYP(NZ).EQ.iplt_annual)THEN
                    GROUP(NB,NZ)=AZMAX1(GROUPI(NZ)-NBTB(NB,NZ))
                  ELSE
                    GROUP(NB,NZ)=GROUPI(NZ)
                  ENDIF
                  exit
                ENDIF
              ENDIF
            ENDDO D120
          ENDIF
        ENDIF
      ENDIF
!
!     ADD AXIS TO ROOT IF PLANT GROWTH STAGE, ROOT NON-STRUCTURAL C
!     CONCENTRATION PERMIT
!
!     PR=nonstructural C concentration needed for root branching
!     FNOD: parameter for allocation of growth to nodes
!     XLI: number of nodes in seed
!     PSTG: node number
!     NB1: number of main branch
!     CEPOLP: canopy nonstructural element concentration
!     PSIRootTurg: root turgor pressure
!     WTRVE: non-structural carbon
      IF(PSIRootTurg(ipltroot,NGTopRootLayer(NZ),NZ).GT.PSILM)THEN
        IF(NRT(NZ).EQ.0 .OR. PSTG(NB1(NZ),NZ).GT.NRT(NZ)/FNOD(NZ)+XTLI(NZ))THEN
          IF((NRT(NZ).EQ.0 .AND. WTRVE(ielmc,NZ).GT.0.0_r8) &
            .OR.(CEPOLP(ielmc,NZ).GT.PR(NZ) .AND. PR(NZ).GT.0.0_r8))THEN
            NRT(NZ)=MIN(JC1,NRT(NZ)+1)
            IDTHR(NZ)=0
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

  integer :: NB,N,L,NE
  real(r8):: ARLSP
  associate(                           &
    CanPLeafShethC   =>  plt_biom%CanPLeafShethC     , &
    CanPBLeafShethC  =>  plt_biom%CanPBLeafShethC    , &
    CanPShootElmMass =>  plt_biom%CanPShootElmMass   , &
    CCPLNP =>  plt_biom%CCPLNP   , &
    CEPOLB =>  plt_biom%CEPOLB   , &
    CEPOLP =>  plt_biom%CEPOLP   , &
    EPOOLP =>  plt_biom%EPOOLP   , &
    EPOLNB =>  plt_biom%EPOLNB   , &
    EPOOLR =>  plt_biom%EPOOLR   , &
    EPOOL  =>  plt_biom%EPOOL    , &
    CEPOLR =>  plt_biom%CEPOLR   , &
    ZEROL  =>  plt_biom%ZEROL    , &
    ZEROP  =>  plt_biom%ZEROP    , &
    WTRTL  =>  plt_biom%WTRTL    , &
    EPOLNP =>  plt_biom%EPOLNP   , &
    WatByPCan  =>  plt_ew%WatByPCan      , &
    VHeatCapCanP  =>  plt_ew%VHeatCapCanP      , &
    NU     =>  plt_site%NU       , &
    IDTHB  =>  plt_pheno%IDTHB   , &
    IDAY   =>  plt_pheno%IDAY    , &
    NB1    =>  plt_morph%NB1     , &
    PrimRootDepth  =>  plt_morph%PrimRootDepth   , &
    MY     =>  plt_morph%MY      , &
    CanPLA  =>  plt_morph%CanPLA   , &
    NGTopRootLayer    =>  plt_morph%NGTopRootLayer     , &
    NIXBotRootLayer   =>  plt_morph%NIXBotRootLayer     , &
    NBR    =>  plt_morph%NBR     , &
    NBTB   =>  plt_morph%NBTB    , &
    HypoctoylHeight  =>  plt_morph%HypoctoylHeight   , &
    SeedinDepth  =>  plt_morph%SeedinDepth   , &
    CanPSA  =>  plt_morph%CanPSA   , &
    NI     =>  plt_morph%NI        &
  )
  plt_bgcr%RootGasLoss_disturb(idg_beg:idg_end-1,NZ)=0.0_r8
  EPOOLP(1:npelms,NZ)=0.0_r8
  NI(NZ)=NIXBotRootLayer(NZ)
  NGTopRootLayer(NZ)=MIN(NI(NZ),MAX(NGTopRootLayer(NZ),NU))
  NB1(NZ)=1
  NBTX=1.0E+06_r8
!
! TOTAL PLANT NON-STRUCTURAL C, N, P
!
! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
! NB1=main branch number
!
  DO NE=1,npelms
    D140: DO NB=1,NBR(NZ)
      IF(IDTHB(NB,NZ).EQ.ibralive)THEN
        EPOOLP(NE,NZ)=EPOOLP(NE,NZ)+EPOOL(NE,NB,NZ)
        EPOLNP(NE,NZ)=EPOLNP(NE,NZ)+EPOLNB(NE,NB,NZ)
      ENDIF
    ENDDO D140
  ENDDO

  DO NB=1,NBR(NZ)
    IF(IDTHB(NB,NZ).EQ.ibralive)THEN
      IF(NBTB(NB,NZ).LT.NBTX)THEN
        NB1(NZ)=NB
        NBTX=NBTB(NB,NZ)
      ENDIF
    ENDIF
  ENDDO
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN ROOT
!
! WTRTL=root mas(g)
! CPOOLR,ZPOOLR,PPOOLR=non-structl C,N,P in root(1),myco(2)(g)
! CCPOLR,CZPOLR,CPPOLR=non-structl C,N,P concn in root(1),myco(2)(g g-1)
!
  D180: DO N=1,MY(NZ)
    D160: DO L=NU,NI(NZ)
      IF(WTRTL(N,L,NZ).GT.ZEROL(NZ))THEN
        DO NE=1,npelms
          CEPOLR(NE,N,L,NZ)=AZMAX1(EPOOLR(NE,N,L,NZ)/WTRTL(N,L,NZ))
        ENDDO
      ELSE
        DO NE=1,npelms
          CEPOLR(NE,N,L,NZ)=1.0_r8
        ENDDO
      ENDIF
    ENDDO D160
  ENDDO D180
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN SHOOT
!
! CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy(g g-1)
! CCPLNP=nonstructural C concentration in canopy nodules
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!
  IF(CanPLeafShethC(NZ).GT.ZEROL(NZ))THEN
    DO NE=1,npelms
      CEPOLP(NE,NZ)=AZMAX1(AMIN1(1.0_r8,EPOOLP(NE,NZ)/CanPLeafShethC(NZ)))
    ENDDO
    CCPLNP(NZ)=AZMAX1(AMIN1(1.0_r8,EPOLNP(ielmc,NZ)/CanPLeafShethC(NZ)))
  ELSE
    CEPOLP(1:npelms,NZ)=1.0_r8
    CCPLNP(NZ)=1.0_r8
  ENDIF
  DO NE=1,npelms
    D190: DO NB=1,NBR(NZ)
      IF(CanPBLeafShethC(NB,NZ).GT.ZEROP(NZ))THEN
        CEPOLB(NE,NB,NZ)=AZMAX1(EPOOL(NE,NB,NZ)/CanPBLeafShethC(NB,NZ))
      ELSE
        CEPOLB(NE,NB,NZ)=1.0_r8
      ENDIF
    ENDDO D190
  ENDDO
!
! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
!
! IDAY(1,=emergence date
! CanPLA,CanPSA=leaf,stalk areas
! HypoctoylHeight=hypocotyledon height
! SeedinDepth=seeding depth
! PrimRootDepth=primary root depth
! VHeatCapCanP,WTSHT,WatByPCan=canopy heat capacity,mass,water content
!
  IF(IDAY(1,NB1(NZ),NZ).EQ.0)THEN
    ARLSP=CanPLA(NZ)+CanPSA(NZ)
    IF((HypoctoylHeight(NZ).GT.SeedinDepth(NZ)).AND.(ARLSP.GT.ZEROL(NZ)) &
      .AND.(PrimRootDepth(1,1,NZ).GT.SeedinDepth(NZ)+ppmc))THEN
      IDAY(1,NB1(NZ),NZ)=I
      VHeatCapCanP(NZ)=cpw*(CanPShootElmMass(ielmc,NZ)*10.0E-06_r8+WatByPCan(NZ))
    ENDIF
  ENDIF
  end associate
  end subroutine stage_phenology_vars
!------------------------------------------------------------------------------------------

  subroutine pft_specific_phenology(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  associate(                           &
    PSICanP   =>  plt_ew%PSICanP     , &
    PSICanPTurg   =>  plt_ew%PSICanPTurg     , &
    DYLN    =>  plt_site%DYLN    , &
    DYLX    =>  plt_site%DYLX    , &
    DYLM    =>  plt_site%DYLM    , &
    ALAT    =>  plt_site%ALAT    , &
    IFLGE   =>  plt_pheno%IFLGE  , &
    IDAY    =>  plt_pheno%IDAY   , &
    IGTYP   =>  plt_pheno%IGTYP  , &
    VRNZ    =>  plt_pheno%VRNZ   , &
    VRNS    =>  plt_pheno%VRNS   , &
    VRNL    =>  plt_pheno%VRNL   , &
    IWTYP   =>  plt_pheno%IWTYP  , &
    TCZ     =>  plt_pheno%TCZ    , &
    VRNY    =>  plt_pheno%VRNY   , &
    TCG     =>  plt_pheno%TCG    , &
    CTC     =>  plt_pheno%CTC    , &
    TCX     =>  plt_pheno%TCX    , &
    VRNF    =>  plt_pheno%VRNF   , &
    VRNX    =>  plt_pheno%VRNX   , &
    IFLGF   =>  plt_pheno%IFLGF    &
  )
!
! CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayerLENGTHENINGTopRootLayerPHOTOPERIODS
!
! IWTYP=phenology type from PFT file
! DYLX,DLYN=daylength of previous,current day
! VRNS,VRNF=leafout,leafoff hours
! VRNY=hourly counter for lengthening photoperiods
! IFLGF=flag for enabling leafoff:0=enable,1=disable
! ALAT=latitude
!
  IF(IWTYP(NZ).EQ.0)THEN
    IF(DYLN.GE.DYLX)THEN
      VRNS(NB,NZ)=VRNY(NB,NZ)
      IF(VRNS(NB,NZ).GE.VRNL(NB,NZ) &
        .OR.(ALAT.GT.0.0_r8.AND.I.EQ.173) &
        .OR.(ALAT.LT.0.0_r8.AND.I.EQ.355))THEN
        VRNF(NB,NZ)=0.0
        IFLGF(NB,NZ)=0
      ENDIF
    ENDIF
!
!   CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayerSHORTENINGTopRootLayerPHOTOPERIODS
!
!   VRNS,VRNF=leafout,leafoff hours
!   VRNZ=hourly counter for shortening photoperiods
!   IFLGE=flag for enabling leafout:0=enable,1=disable
!   ALAT=latitude
!
    IF(DYLN.LT.DYLX)THEN
      VRNF(NB,NZ)=VRNZ(NB,NZ)
      IF(VRNF(NB,NZ).GE.VRNX(NB,NZ) &
        .OR.(ALAT.GT.0.0_r8.AND.I.EQ.355) &
        .OR.(ALAT.LT.0.0_r8.AND.I.EQ.173))THEN
        VRNS(NB,NZ)=0.0
        IFLGE(NB,NZ)=0
      ENDIF
    ENDIF
!
!   CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS ABOVE
!   SPECIFIED TEMPERATURE DURINGTopRootLayerLENGTHENINGTopRootLayerPHOTOPERIODS
!
!   IWTYP=phenology type from PFT file
!   DYLX,DLYN=daylength of previous,current day
!   VRNS,VRNL=leafout hours,hours required for leafout
!   VRNF,VRNX=leafoff hours,hours required for leafoff
!   IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!   TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!   ALAT=latitude
!   IDAY(2,=date of floral initiation
!
  ELSEIF(IWTYP(NZ).EQ.1)THEN
    IF((DYLN.GE.DYLX.OR.(DYLN.LT.DYLX.AND.VRNF(NB,NZ).LT.VRNX(NB,NZ))) &
      .AND.IFLGE(NB,NZ).EQ.0)THEN
        IF(TCG(NZ).GE.TCZ(NZ))THEN
          VRNS(NB,NZ)=VRNS(NB,NZ)+1.0
        ENDIF
        IF(VRNS(NB,NZ).LT.VRNL(NB,NZ))THEN
          IF(TCG(NZ).LT.CTC(NZ))THEN
            VRNS(NB,NZ)=AZMAX1(VRNS(NB,NZ)-1.0)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ).GE.VRNL(NB,NZ) &
          .OR.(ALAT.GT.0.0.AND.I.EQ.173) &
          .OR.(ALAT.LT.0.0.AND.I.EQ.355))THEN
          VRNF(NB,NZ)=0.0
        ENDIF
      ENDIF
      IF(IDAY(2,NB,NZ).NE.0.OR.(DYLN.LT.DYLX.AND.DYLN.LT.12.0))THEN
        IFLGF(NB,NZ)=0
      ENDIF
!
!     CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS BELOW
!     SPECIFIED TEMPERATURE DURINGTopRootLayerSHORTENINGTopRootLayerPHOTOPERIODS
!
!     DYLX,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     IDAY(2,=date of floral initiation
!
      IF(DYLN.LT.DYLX.AND.IFLGF(NB,NZ).EQ.0 &
        .AND.IDAY(2,NB,NZ).NE.0)THEN
        IF(TCG(NZ).LE.TCX(NZ))THEN
          VRNF(NB,NZ)=VRNF(NB,NZ)+1.0
        ENDIF
        IF(VRNF(NB,NZ).GE.VRNX(NB,NZ).AND.IFLGE(NB,NZ).EQ.1)THEN
          VRNS(NB,NZ)=0.0
          IFLGE(NB,NZ)=0
        ENDIF
      ENDIF

!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS
!     ABOVE SPECIFIED WATER POTENTIAL DURINGTopRootLayerDORMANCY
!
!     IWTYP=phenology type from PFT file
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF=leafoff hours
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     IDAY(2,=date of floral initiation
!
    ELSEIF(IWTYP(NZ).EQ.2.OR.IWTYP(NZ).EQ.4 &
      .OR.IWTYP(NZ).EQ.5)THEN
      IF(IFLGE(NB,NZ).EQ.0)THEN
        IF(PSICanP(NZ).GE.PSILX)THEN
          VRNS(NB,NZ)=VRNS(NB,NZ)+1.0
        ENDIF
        IF(VRNS(NB,NZ).LT.VRNL(NB,NZ))THEN
          IF(PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
            VRNS(NB,NZ)=AZMAX1(VRNS(NB,NZ)-12.0)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
          VRNF(NB,NZ)=0.0
          IF(IDAY(2,NB,NZ).NE.0)IFLGF(NB,NZ)=0
        ENDIF
      ENDIF
      IF(IDAY(2,NB,NZ).NE.0)IFLGF(NB,NZ)=0
!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS
!     BELOW SPECIFIED WATER POTENTIAL DURINGTopRootLayerGROWINGTopRootLayerSEASON
!
!     VRNS=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     IDAY(2,=date of floral initiation
!     VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!     VRNE=maximum hours for leafout,leafoff
!
      IF(IFLGE(NB,NZ).EQ.1.AND.IFLGF(NB,NZ).EQ.0)THEN
        IF(PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
          VRNF(NB,NZ)=VRNF(NB,NZ)+1.0
        ENDIF
        IF(IWTYP(NZ).EQ.4)THEN
          IF(VRNZ(NB,NZ).GT.VRNE)THEN
            VRNF(NB,NZ)=VRNZ(NB,NZ)
          ENDIF
        ELSEIF(IWTYP(NZ).EQ.5)THEN
          IF(VRNY(NB,NZ).GT.VRNE)THEN
            VRNF(NB,NZ)=VRNY(NB,NZ)
          ENDIF
        ENDIF
        IF(VRNF(NB,NZ).GE.VRNX(NB,NZ) &
          .AND.IFLGE(NB,NZ).EQ.1)THEN
          VRNS(NB,NZ)=0.0
          IFLGE(NB,NZ)=0
        ENDIF
      ENDIF
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS ABOVE SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     LENGTHENINGTopRootLayerPHOTOPERIODS
!
!     IWTYP=phenology type from PFT file
!     DYLX,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     IDAY(2,=date of floral initiation
!
    ELSEIF(IWTYP(NZ).EQ.3)THEN
      IF((DYLN.GE.DYLX.OR.DYLN.GE.DYLM-2.0_r8).AND.IFLGE(NB,NZ).EQ.0)THEN
        IF(TCG(NZ).GE.TCZ(NZ).AND.PSICanPTurg(NZ).GT.PSILM)THEN
          VRNS(NB,NZ)=VRNS(NB,NZ)+1.0
        ENDIF
        IF(VRNS(NB,NZ).LT.VRNL(NB,NZ))THEN
          IF(TCG(NZ).LT.CTC(NZ) &
            .OR.PSICanPTurg(NZ).LT.PSILM)THEN
            VRNS(NB,NZ)=AZMAX1(VRNS(NB,NZ)-1.5)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
          VRNF(NB,NZ)=0.0
          IF(IDAY(2,NB,NZ).NE.0)IFLGF(NB,NZ)=0
        ENDIF
      ENDIF
      IF(IDAY(2,NB,NZ).NE.0)IFLGF(NB,NZ)=0
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS BELOW SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     SHORTENINGTopRootLayerPHOTOPERIODS
!
!     DYLX,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IFLGE,IFLGF=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     IDAY(2,=date of floral initiation
!
      IF((DYLN.LT.DYLX.OR.DYLN.LT.24.0-DYLM+2.0_r8).AND.IFLGF(NB,NZ).EQ.0)THEN
        IF(TCG(NZ).LE.TCX(NZ).OR.PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
          VRNF(NB,NZ)=VRNF(NB,NZ)+1.0
        ENDIF
        IF(VRNF(NB,NZ).GE.VRNX(NB,NZ) &
          .AND.IFLGE(NB,NZ).EQ.1)THEN
        VRNS(NB,NZ)=0.0
        IFLGE(NB,NZ)=0
      ENDIF
    ENDIF
  ENDIF
  end associate
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
  associate(                        &
    IYR0    =>  plt_distb%IYR0    , &
    IDAY0   =>  plt_distb%IDAY0   , &
    SnowDepth   =>  plt_ew%SnowDepth      , &
    PSICanP   =>  plt_ew%PSICanP      , &
    ZERO    =>  plt_site%ZERO     , &
    DYLN    =>  plt_site%DYLN     , &
    IYRC    =>  plt_site%IYRC     , &
    DYLX    =>  plt_site%DYLX     , &
    TKG     =>  plt_pheno%TKG     , &
    IDAY    =>  plt_pheno%IDAY    , &
    XRLA    =>  plt_pheno%XRLA    , &
    GSTGF   =>  plt_pheno%GSTGF   , &
    OFFST   =>  plt_pheno%OFFST   , &
    ISTYP   =>  plt_pheno%ISTYP   , &
    IFLGE   =>  plt_pheno%IFLGE   , &
    VSTG    =>  plt_morph%VSTG    , &
    GSTGI   =>  plt_pheno%GSTGI   , &
    DGSTGF  =>  plt_pheno%DGSTGF  , &
    IFLGG   =>  plt_pheno%IFLGG   , &
    VSTGX   =>  plt_pheno%VSTGX   , &
    IDTYP   =>  plt_pheno%IDTYP   , &
    TGSTGF  =>  plt_pheno%TGSTGF  , &
    XDL     =>  plt_pheno%XDL     , &
    DGSTGI  =>  plt_pheno%DGSTGI  , &
    XPPD    =>  plt_pheno%XPPD    , &
    TGSTGI  =>  plt_pheno%TGSTGI  , &
    GROUP   =>  plt_pheno%GROUP   , &
    IWTYP   =>  plt_pheno%IWTYP   , &
    VRNL    =>  plt_pheno%VRNL    , &
    IPTYP   =>  plt_pheno%IPTYP   , &
    IFLGA   =>  plt_pheno%IFLGA   , &
    GROUPI  =>  plt_pheno%GROUPI  , &
    PlantO2Stress    =>  plt_pheno%PlantO2Stress    , &
    VRNS    =>  plt_pheno%VRNS    , &
    VRNF    =>  plt_pheno%VRNF    , &
    VRNX    =>  plt_pheno%VRNX    , &
    XRNI    =>  plt_pheno%XRNI    , &
    CanopyHeight      =>  plt_morph%CanopyHeight      , &
    PSTG    =>   plt_morph%PSTG   , &
    PSTGI   =>   plt_morph%PSTGI  , &
    PSTGF   =>   plt_morph%PSTGF  , &
    NB1     =>   plt_morph%NB1      &
  )
  IF(IDAY(1,NB,NZ).EQ.0)THEN
    IDAY(1,NB,NZ)=I
    IFLGA(NB,NZ)=1
    IFLGE(NB,NZ)=0
    VRNS(NB,NZ)=0.5_r8*VRNS(NB1(NZ),NZ)
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
  IF(IWTYP(NZ).EQ.0.OR.VRNF(NB,NZ).LT.VRNX(NB,NZ))THEN
    TKCO=TKG(NZ)+OFFST(NZ)
    RTK=RGAS*TKCO
    STK=710.0_r8*TKCO
    ACTV=1+EXP((197500_r8-STK)/RTK)+EXP((STK-218500._r8)/RTK)
    TFNP=EXP(24.269_r8-60000._r8/RTK)/ACTV
    RNI=AZMAX1(XRNI(NZ)*TFNP)
    RLA=AZMAX1(XRLA(NZ)*TFNP)
!
!   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
!
!   PSICanPTurg=leaf turgor potential
!   WFNG=water stress effect on phenology
!
    IF(ISTYP(NZ).EQ.iplt_annual)THEN
      IF(IDAY(6,NB,NZ).EQ.0)THEN
        WFNG=EXP(0.025_r8*PSICanP(NZ))
        RNI=RNI*WFNG
        RLA=RLA*WFNG
      ENDIF
      IF(IDAY(2,NB,NZ).EQ.0)THEN
        OFNG=SQRT(PlantO2Stress(NZ))
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
    PSTG(NB,NZ)=PSTG(NB,NZ)+RNI
    VSTG(NB,NZ)=VSTG(NB,NZ)+RLA
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
    IF(IDAY(2,NB,NZ).NE.0)THEN
      GSTGI(NB,NZ)=(PSTG(NB,NZ)-PSTGI(NB,NZ))/GROUPI(NZ)
      DGSTGI(NB,NZ)=RNI/(GROUPI(NZ)*GSTGG)
      TGSTGI(NB,NZ)=TGSTGI(NB,NZ)+DGSTGI(NB,NZ)
    ENDIF
    IF(IDAY(6,NB,NZ).NE.0)THEN
      GSTGF(NB,NZ)=(PSTG(NB,NZ)-PSTGF(NB,NZ))/GROUPI(NZ)
      DGSTGF(NB,NZ)=RNI/(GROUPI(NZ)*GSTGR)
      TGSTGF(NB,NZ)=TGSTGF(NB,NZ)+DGSTGF(NB,NZ)
    ENDIF
    IFLGG(NB,NZ)=1
  ELSE
    IFLGG(NB,NZ)=0
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
! ZC,SnowDepth=canopy height,snowpack depth
!
  IF(IDAY(2,NB,NZ).EQ.0)THEN
    IF(PSTG(NB,NZ).GT.GROUP(NB,NZ)+PSTGI(NB,NZ) &
      .AND.((VRNS(NB,NZ).GE.VRNL(NB,NZ)) &
      .OR.(I.GE.IDAY0(NZ).AND.IYRC.EQ.IYR0(NZ) &
      .AND.DYLN.GT.DYLX)) &
      .OR.(((ISTYP(NZ).EQ.iplt_preanu.AND.(IWTYP(NZ).EQ.1 &
      .OR.IWTYP(NZ).EQ.3)) &
      .OR.(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).EQ.0)) &
      .AND.CanopyHeight(NZ).GE.SnowDepth-ZERO &
      .AND.DYLN.LT.DYLX))THEN
!
!     FINAL VEGETATIVE NODE NUMBER DEPENDS ON PHOTOPERIOD FROM 'DAY'
!     AND ON MATURITY GROUP, CRITICAL PHOTOPERIOD AND PHOTOPERIOD
!     SENSITIVITY ENTERED IN 'READQ'
!
!     IPTYP=photoperiod type from PFT file
!     PPD=photoperiod sensitivity
!     XDL=critical photoperiod from PFT file
!     IDAY(2,=date of floral initiation
!     VSTGX=node number on date of floral initiation
!
      IF(IPTYP(NZ).EQ.0)THEN
        PPD=0.0_r8
      ELSE
        PPD=AZMAX1(XDL(NZ)-DYLN)
        IF(IPTYP(NZ).EQ.1.AND.DYLN.GE.DYLX)PPD=0.0_r8
      ENDIF

      IF(IPTYP(NZ).EQ.0 &
        .OR.(IPTYP(NZ).EQ.1.AND.PPD.GT.XPPD(NZ)) &
        .OR.(IPTYP(NZ).EQ.2.AND.PPD.LT.XPPD(NZ)) &
        .OR.(((ISTYP(NZ).EQ.iplt_preanu.AND.(IWTYP(NZ).EQ.1 &
        .OR.IWTYP(NZ).EQ.3)) &
        .OR.(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).EQ.0)) &
        .AND.CanopyHeight(NZ).GE.SnowDepth-ZERO &
        .AND.DYLN.LT.DYLX))THEN
        IDAY(2,NB,NZ)=I
        PSTGI(NB,NZ)=PSTG(NB,NZ)
        IF(ISTYP(NZ).EQ.iplt_annual.AND.IDTYP(NZ).EQ.0)THEN
          VSTGX(NB,NZ)=PSTG(NB,NZ)
        ENDIF
      ENDIF
    ENDIF
!
!   STEM ELONGATION
!
!   GSTGI=vegetative node number normalized for maturity group
!   GSTGG=normalized growth stage durations for vegetative phenology
!   IDAY(3,=start of stem elongation and setting max seed number
!
  ELSEIF(IDAY(3,NB,NZ).EQ.0)THEN
    IF(GSTGI(NB,NZ).GT.0.25_r8*GSTGG &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DYLN.LT.DYLX.AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ))THEN
      IDAY(3,NB,NZ)=I
    ENDIF
!
!   IDAY(4,=mid stem elongation
!
  ELSEIF(IDAY(4,NB,NZ).EQ.0)THEN
    IF(GSTGI(NB,NZ).GT.0.50*GSTGG &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DYLN.LT.DYLX.AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ))THEN
      IDAY(4,NB,NZ)=I
      IF(ISTYP(NZ).EQ.iplt_annual.AND.IDTYP(NZ).NE.0)THEN
        VSTGX(NB,NZ)=PSTG(NB,NZ)
      ENDIF
    ENDIF
!
!   IDAY(5,=end of stem elongation and setting max seed number
!
  ELSEIF(IDAY(5,NB,NZ).EQ.0)THEN
    IF(GSTGI(NB,NZ).GT.1.00*GSTGG &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DYLN.LT.DYLX.AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ))THEN
      IDAY(5,NB,NZ)=I
    ENDIF
!
!   ANTHESIS OCCURS WHEN THE NUMBER OF LEAVES THAT HAVE APPEARED
!   EQUALS THE NUMBER OF NODES INITIATED WHEN THE FINAL VEGETATIVE
!   NODE NUMBER WAS SET ABOVE
!
!   IDAY(6,=start of anthesis and setting final seed number
!   VSTG=number of leaves appeared
!   PSTGI=node number at floral initiation
!   ISTYP,IWTYP,IPTYP=growth habit,phenology,photoperiod type from PFT file
!   IFLGE=flag for enabling leafout:0=enable,1=disable
!   VRNF,VRNX=leafoff hours,hours required for leafoff
!   DYLX,DLYN=daylength of previous,current day
!   PSTGF=number of nodes at anthesis
!
  ELSEIF(IDAY(6,NB,NZ).EQ.0)THEN
    IF((VSTG(NB,NZ).GT.PSTGI(NB,NZ)) &
      .OR.(ISTYP(NZ).NE.iplt_annual.AND.IDAY(5,NB,NZ).NE.0) &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DYLN.LT.DYLX.AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ))THEN
      IF(NB.EQ.NB1(NZ).OR.IDAY(6,NB1(NZ),NZ).NE.0)THEN
        IDAY(6,NB,NZ)=I
        PSTGF(NB,NZ)=PSTG(NB,NZ)
      ENDIF
    ENDIF
!
!   START GRAIN FILL PERIOD
!
!   IDAY(7,=start of grain filling and setting max seed size
!   GSTGF=reproductive node number normalized for maturity group
!   GSTGR=normalized growth stage durations for reproductive phenology
!
!
  ELSEIF(IDAY(7,NB,NZ).EQ.0)THEN
    IF(GSTGF(NB,NZ).GT.0.50*GSTGR &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DYLN.LT.DYLX.AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.IFLGE(NB,NZ).EQ.1 &
      .AND.VRNF(NB,NZ).GT.VRNX(NB,NZ))THEN
        IDAY(7,NB,NZ)=I
    ENDIF
!
!   END SEED NUMBER SET PERIOD
!
!   IDAY(8,=end date setting for final seed number
!
  ELSEIF(IDAY(8,NB,NZ).EQ.0)THEN
    IF(GSTGF(NB,NZ).GT.1.00*GSTGR)THEN
      IDAY(8,NB,NZ)=I
    ENDIF
!
!   END SEED SIZE SET PERIOD
!
!   IDAY(9,=end of setting max seed size
!
  ELSEIF(IDAY(9,NB,NZ).EQ.0)THEN
    IF(GSTGF(NB,NZ).GT.1.50*GSTGR)THEN
      IDAY(9,NB,NZ)=I
    ENDIF
  ENDIF
  end associate
  end subroutine living_branch_phenology

end module HfuncsMod
