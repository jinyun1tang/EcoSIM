module HfuncMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StartqMod    , only : startq
  use EcosimConst
  use GridDataType
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use PlantTraitDataType
  use PlantMngmtDataType
  use PlantDataRateType
  use SnowDataType
  use CanopyDataType
  use RootDataType
  use SOMDataType
  use EcoSIMHistMod
  use EcosysBGCFluxType
  implicit none

  private


  character(len=*), parameter :: mod_filename = __FILE__
  real(r8) :: PSILY(0:3)
  real(r8) :: ARLSP,ACTV,OFNG,PPD,RTK,RNI,RLA,STK,TKCO,TFNP
  real(r8) :: WFNG
  integer :: KVSTGX,NB,NBTX,N
  integer :: NBX(0:3)
!
! PSILM=minimum canopy turgor potential for leaf expansion (MPa)
! PSILX=minimum canopy water potential for leafout of drought-deciduous PFT (MPa)
! PSILY=minimum canopy water potential for leafoff of drought-deciduous PFT (MPa
! GSTGG,GSTGR=normalized growth stage durations for vegetative,reproductive phenology
! NBX=maximum branch number for PFT defined by IBTYP in PFT file
! VRNE=maximum hours for leafout,leafoff
!
  real(r8), PARAMETER :: PSILM=0.1,PSILX=-0.2
  real(r8) ,PARAMETER :: GSTGG=2.00,GSTGR=0.667,VRNE=3600.0
  DATA PSILY/-200.0,-2.0,-2.0,-2.0/
  DATA NBX /5,1,1,1/

  public :: hfunc
  contains

  SUBROUTINE hfunc(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES PLANT PHENOLOGY
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  INTEGER :: NZ,NY,NX

! begin_execution
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      DO 9985 NZ=1,NP(NY,NX)
!       WRITE(*,4444)'IFLGC',I,J,NX,NY,NZ,DATAP(NZ,NY,NX),IYRC
!      2,IDAY0(NZ,NY,NX),IDAYH(NZ,NY,NX),IYR0(NZ,NY,NX),IYRH(NZ,NY,NX)
!      3,IDTH(NZ,NY,NX),IFLGC(NZ,NY,NX),IFLGT(NY,NX)
!4444  FORMAT(A8,5I8,A16,20I8)
        IF(DATAP(NZ,NY,NX).NE.'NO')THEN
!
!         PPT=total biome population
!
          PPT(NY,NX)=PPT(NY,NX)+PP(NZ,NY,NX)
!
!         SET CROP FLAG ACCORDING TO PLANTING, HARVEST DATES, DEATH,
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
          call set_flags(I,J,NZ,NY,NX)
!
!         INITIALIZE VARIABLES IN ACTIVE PFT
!
          IF(IFLGC(NZ,NY,NX).EQ.1)THEN

            call stage_phenology_vars(I,J,NZ,NY,NX)

            call root_shoot_branching(I,J,NZ,NY,NX)
!
!           THE REST OF THE SUBROUTINE MODELS THE PHENOLOGY OF EACH BRANCH
!
!           IFLGA,IFLGE=flags for initializing leafout,leafoff
!           VRNS=leafout hours
!
            IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0 &
              .OR.IFLGI(NZ,NY,NX).EQ.1)THEN
              DO 2010 NB=1,NBR(NZ,NY,NX)
                IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
                  call living_branch_phenology(I,J,NB,nz,ny,nx)
                ENDIF
!
!               KVSTG=integer of most recent leaf number currently growing
!               IFLGP=flag for remobilization
!
                KVSTGX=KVSTG(NB,NZ,NY,NX)
                IF(VSTGX(NB,NZ,NY,NX).LE.1.0E-06)THEN
                  KVSTG(NB,NZ,NY,NX)=INT(VSTG(NB,NZ,NY,NX))+1
                ELSE
                  KVSTG(NB,NZ,NY,NX)=INT(AMIN1(VSTG(NB,NZ,NY,NX) &
                    ,VSTGX(NB,NZ,NY,NX)))+1
                ENDIF
                KLEAF(NB,NZ,NY,NX)=MIN(24,KVSTG(NB,NZ,NY,NX))
                IF(KVSTG(NB,NZ,NY,NX).GT.KVSTGX)THEN
                  IFLGP(NB,NZ,NY,NX)=1
                ELSE
                  IFLGP(NB,NZ,NY,NX)=0
                ENDIF
!
!               PHENOLOGY
!
!               DYLX,DLYN=daylength of previous,current day
!               VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!
                IF(IDTHB(NB,NZ,NY,NX).EQ.0.OR.IFLGI(NZ,NY,NX).EQ.1)THEN
                  IF(DYLN(NY,NX).GE.DYLX(NY,NX))THEN
                    VRNY(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)+1.0
                    VRNZ(NB,NZ,NY,NX)=0.0
                  ELSE
                    VRNY(NB,NZ,NY,NX)=0.0
                    VRNZ(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)+1.0
                  ENDIF

                  call pft_specific_phenology(I,J,nz,ny,nx)

                ENDIF
2010          CONTINUE
!
!             WATER STRESS INDICATOR
!
!             PSILT=canopy total water potential
!             PSILY=minimum canopy water potential for leafoff
!             WSTR=number of hours PSILT < PSILY (for output only)
!
              IF(PSILT(NZ,NY,NX).LT.PSILY(IGTYP(NZ,NY,NX)))THEN
                WSTR(NZ,NY,NX)=WSTR(NZ,NY,NX)+1.0
              ENDIF
            ENDIF
          ENDIF
        ENDIF
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
  RETURN
  END subroutine hfunc
!------------------------------------------------------------------------------------------

  subroutine set_flags(I,J,NZ,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX

  INTEGER :: L

! begin_execution
  IF(J.EQ.1)THEN
    IF(IDAY0(NZ,NY,NX).LE.IDAYH(NZ,NY,NX).OR.IYR0(NZ,NY,NX).LT.IYRH(NZ,NY,NX))THEN
      IF(I.GE.IDAY0(NZ,NY,NX).OR.IDATA(3).GT.IYR0(NZ,NY,NX))THEN
        IF(I.GT.IDAYH(NZ,NY,NX).AND.IYRC.GE.IYRH(NZ,NY,NX).AND.IDTH(NZ,NY,NX).EQ.1)THEN
          IFLGC(NZ,NY,NX)=0
        ELSE
          IF(I.EQ.IDAY0(NZ,NY,NX).AND.IDATA(3).EQ.IYR0(NZ,NY,NX))THEN
            IFLGC(NZ,NY,NX)=0
            IDTH(NZ,NY,NX)=0
            CALL STARTQ(NX,NX,NY,NY,NZ,NZ)
            TNBP(NY,NX)=TNBP(NY,NX)+WTRVX(NZ,NY,NX)
          ENDIF

          IF(DATAP(NZ,NY,NX).NE.'NO'.AND.IDTH(NZ,NY,NX).EQ.0)IFLGC(NZ,NY,NX)=1
        ENDIF
      ELSE
        IFLGC(NZ,NY,NX)=0
      ENDIF
    ELSE
      IF((I.LT.IDAY0(NZ,NY,NX).AND.I.GT.IDAYH(NZ,NY,NX) &
        .AND.IYRC.GE.IYRH(NZ,NY,NX).AND.IDTH(NZ,NY,NX).EQ.1) &
        .OR.(I.LT.IDAY0(NZ,NY,NX).AND.IYR0(NZ,NY,NX) &
        .GT.IYRH(NZ,NY,NX)))THEN
        IFLGC(NZ,NY,NX)=0
      ELSE
        IF(I.EQ.IDAY0(NZ,NY,NX).AND.IDATA(3).EQ.IYR0(NZ,NY,NX))THEN
          IFLGC(NZ,NY,NX)=0
          IDTH(NZ,NY,NX)=0
          CALL STARTQ(NX,NX,NY,NY,NZ,NZ)
          TNBP(NY,NX)=TNBP(NY,NX)+WTRVX(NZ,NY,NX)
        ENDIF
        IF(DATAP(NZ,NY,NX).NE.'NO'.AND.IDTH(NZ,NY,NX).EQ.0)IFLGC(NZ,NY,NX)=1
      ENDIF
    ENDIF
    IFLGT(NY,NX)=IFLGT(NY,NX)+IFLGC(NZ,NY,NX)
  ENDIF
  end subroutine set_flags
!------------------------------------------------------------------------------------------

  subroutine root_shoot_branching(I,J,NZ,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX

! begin_execution

!
! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
! CONCENTRATION PERMIT
!
! IFLGI=PFT initialization flag:0=no,1=yes
! PSIRG=root turgor potential
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
!
! WRITE(*,224)'HFUNC',I,J,IFLGI(NZ,NY,NX),PP(NZ,NY,NX)
!    2,TCG(NZ,NY,NX),PSIRG(1,NG(NZ,NY,NX),NZ,NY,NX)
!    3,PSILM,ISTYP(NZ,NY,NX),IDAY(2,NB1(NZ,NY,NX),NZ,NY,NX)
!    4,NBR(NZ,NY,NX),WTRVC(NZ,NY,NX),CCPOLP(NZ,NY,NX)
!    5,PB(NZ,NY,NX),IDTHB(NB,NZ,NY,NX),NB1(NZ,NY,NX)
!    6,PSTG(NB1(NZ,NY,NX),NZ,NY,NX),NBT(NZ,NY,NX)
!    7,NNOD(NZ,NY,NX),FNOD(NZ,NY,NX),XTLI(NZ,NY,NX)
!224   FORMAT(A8,3I6,5E12.4,3I6,3E12.4,2I6,1E12.4,2I6,2E12.4)
  IF(IFLGI(NZ,NY,NX).EQ.0)THEN
    IF(J.EQ.1.AND.PP(NZ,NY,NX).GT.0.0)THEN
      IF(PSIRG(1,NG(NZ,NY,NX),NZ,NY,NX).GT.PSILM)THEN
        IF(ISTYP(NZ,NY,NX).NE.0.OR.IDAY(2,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
          IF((NBR(NZ,NY,NX).EQ.0.AND.WTRVC(NZ,NY,NX).GT.0.0) &
            .OR.(CCPOLP(NZ,NY,NX).GT.PB(NZ,NY,NX).AND.PB(NZ,NY,NX).GT.0.0))THEN
            DO 120 NB=1,10
              IF(IDTHB(NB,NZ,NY,NX).EQ.1)THEN
                IF(NB.EQ.NB1(NZ,NY,NX) &
                  .OR.PSTG(NB1(NZ,NY,NX),NZ,NY,NX).GT.NBT(NZ,NY,NX) &
                  +NNOD(NZ,NY,NX)/FNOD(NZ,NY,NX)+XTLI(NZ,NY,NX))THEN
                  NBT(NZ,NY,NX)=NBT(NZ,NY,NX)+1
                  NBR(NZ,NY,NX)=MIN(NBX(IBTYP(NZ,NY,NX)),MAX(NB,NBR(NZ,NY,NX)))
                  NBTB(NB,NZ,NY,NX)=NBT(NZ,NY,NX)-1
                  IDTHP(NZ,NY,NX)=0
                  IDTHB(NB,NZ,NY,NX)=0
                  VRNS(NB,NZ,NY,NX)=0.0
                  IF(ISTYP(NZ,NY,NX).EQ.0)THEN
                    GROUP(NB,NZ,NY,NX)=AMAX1(0.0,GROUPI(NZ,NY,NX)-NBTB(NB,NZ,NY,NX))
                  ELSE
                    GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
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
      IF(PSIRG(1,NG(NZ,NY,NX),NZ,NY,NX).GT.PSILM)THEN
        IF(NRT(NZ,NY,NX).EQ.0.OR.PSTG(NB1(NZ,NY,NX),NZ,NY,NX) &
          .GT.NRT(NZ,NY,NX)/FNOD(NZ,NY,NX)+XTLI(NZ,NY,NX))THEN
          IF((NRT(NZ,NY,NX).EQ.0.AND.WTRVC(NZ,NY,NX).GT.0.0) &
            .OR.(CCPOLP(NZ,NY,NX).GT.PR(NZ,NY,NX).AND.PR(NZ,NY,NX).GT.0.0))THEN
            NRT(NZ,NY,NX)=MIN(10,NRT(NZ,NY,NX)+1)
            IDTHR(NZ,NY,NX)=0
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
!2224  FORMAT(A8,6I4)
  end subroutine root_shoot_branching
!------------------------------------------------------------------------------------------

  subroutine stage_phenology_vars(I,J,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX

  integer :: NB,N,L

  RCO2Z(NZ,NY,NX)=0.0
  ROXYZ(NZ,NY,NX)=0.0
  RCH4Z(NZ,NY,NX)=0.0
  RN2OZ(NZ,NY,NX)=0.0
  RNH3Z(NZ,NY,NX)=0.0
  RH2GZ(NZ,NY,NX)=0.0
  CPOOLP(NZ,NY,NX)=0.0
  ZPOOLP(NZ,NY,NX)=0.0
  PPOOLP(NZ,NY,NX)=0.0
  NI(NZ,NY,NX)=NIX(NZ,NY,NX)
  NG(NZ,NY,NX)=MIN(NI(NZ,NY,NX),MAX(NG(NZ,NY,NX),NU(NY,NX)))
  NB1(NZ,NY,NX)=1
  NBTX=1.0E+06
!
! TOTAL PLANT NON-STRUCTURAL C, N, P
!
! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
! NB1=main branch number
!
  DO 140 NB=1,NBR(NZ,NY,NX)
    IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      CPOOLP(NZ,NY,NX)=CPOOLP(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX)
      ZPOOLP(NZ,NY,NX)=ZPOOLP(NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX)
      PPOOLP(NZ,NY,NX)=PPOOLP(NZ,NY,NX)+PPOOL(NB,NZ,NY,NX)
      CPOLNP(NZ,NY,NX)=CPOLNP(NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
      ZPOLNP(NZ,NY,NX)=ZPOLNP(NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX)
      PPOLNP(NZ,NY,NX)=PPOLNP(NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX)
      IF(NBTB(NB,NZ,NY,NX).LT.NBTX)THEN
        NB1(NZ,NY,NX)=NB
        NBTX=NBTB(NB,NZ,NY,NX)
      ENDIF
    ENDIF
140   CONTINUE
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN ROOT
!
! WTRTL=root mass(g)
! CPOOLR,ZPOOLR,PPOOLR=non-structl C,N,P in root(1),myco(2)(g)
! CCPOLR,CZPOLR,CPPOLR=non-structl C,N,P concn in root(1),myco(2)(g g-1)
!
  DO 180 N=1,MY(NZ,NY,NX)
    DO 160 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(WTRTL(N,L,NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
        CCPOLR(N,L,NZ,NY,NX)=AMAX1(0.0,CPOOLR(N,L,NZ,NY,NX)/WTRTL(N,L,NZ,NY,NX))
        CZPOLR(N,L,NZ,NY,NX)=AMAX1(0.0,ZPOOLR(N,L,NZ,NY,NX)/WTRTL(N,L,NZ,NY,NX))
        CPPOLR(N,L,NZ,NY,NX)=AMAX1(0.0,PPOOLR(N,L,NZ,NY,NX)/WTRTL(N,L,NZ,NY,NX))
!       CCPOLR(N,L,NZ,NY,NX)=AMIN1(1.0,CCPOLR(N,L,NZ,NY,NX))
      ELSE
        CCPOLR(N,L,NZ,NY,NX)=1.0
        CZPOLR(N,L,NZ,NY,NX)=1.0
        CPPOLR(N,L,NZ,NY,NX)=1.0
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
  IF(WTLS(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
    CCPOLP(NZ,NY,NX)=AMAX1(0.0,AMIN1(1.0,CPOOLP(NZ,NY,NX)/WTLS(NZ,NY,NX)))
    CCPLNP(NZ,NY,NX)=AMAX1(0.0,AMIN1(1.0,CPOLNP(NZ,NY,NX)/WTLS(NZ,NY,NX)))
    CZPOLP(NZ,NY,NX)=AMAX1(0.0,AMIN1(1.0,ZPOOLP(NZ,NY,NX)/WTLS(NZ,NY,NX)))
    CPPOLP(NZ,NY,NX)=AMAX1(0.0,AMIN1(1.0,PPOOLP(NZ,NY,NX)/WTLS(NZ,NY,NX)))
  ELSE
    CCPOLP(NZ,NY,NX)=1.0
    CCPLNP(NZ,NY,NX)=1.0
    CZPOLP(NZ,NY,NX)=1.0
    CPPOLP(NZ,NY,NX)=1.0
  ENDIF
  DO 190 NB=1,NBR(NZ,NY,NX)
    IF(WTLSB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CCPOLB(NB,NZ,NY,NX)=AMAX1(0.0,CPOOL(NB,NZ,NY,NX)/WTLSB(NB,NZ,NY,NX))
      CZPOLB(NB,NZ,NY,NX)=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX)/WTLSB(NB,NZ,NY,NX))
      CPPOLB(NB,NZ,NY,NX)=AMAX1(0.0,PPOOL(NB,NZ,NY,NX)/WTLSB(NB,NZ,NY,NX))
    ELSE
      CCPOLB(NB,NZ,NY,NX)=1.0
      CZPOLB(NB,NZ,NY,NX)=1.0
      CPPOLB(NB,NZ,NY,NX)=1.0
    ENDIF
190 CONTINUE
!
! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
!
! IDAY(1,=emergence date
! ARLFP,ARSTP=leaf,stalk areas
! HTCTL=hypocotyledon height
! SDPTH=seeding depth
! RTDP1=primary root depth
! VHCPC,WTSHT,VOLWC=canopy heat capacity,mass,water content
!
! WRITE(*,223)'EMERG',I,J,NZ,NB1(NZ,NY,NX)
!   2,IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX),HTCTL(NZ,NY,NX),SDPTH(NZ,NY,NX)
!   3,ARLSP,RTDP1(1,1,NZ,NY,NX)
!223   FORMAT(A8,5I4,12E12.4)
  IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
    ARLSP=ARLFP(NZ,NY,NX)+ARSTP(NZ,NY,NX)
    IF((HTCTL(NZ,NY,NX).GT.SDPTH(NZ,NY,NX)) &
      .AND.(ARLSP.GT.ZEROL(NZ,NY,NX)) &
      .AND.(RTDP1(1,1,NZ,NY,NX).GT.SDPTH(NZ,NY,NX)+1.0E-06))THEN
      IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX)=I
      VHCPC(NZ,NY,NX)=cpw*(WTSHT(NZ,NY,NX)*10.0E-06+VOLWC(NZ,NY,NX))
    ENDIF
  ENDIF
  end subroutine stage_phenology_vars
!------------------------------------------------------------------------------------------

  subroutine pft_specific_phenology(I,J,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX

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
  IF(IWTYP(NZ,NY,NX).EQ.0)THEN
    IF(DYLN(NY,NX).GE.DYLX(NY,NX))THEN
      VRNS(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)
      IF(VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX) &
        .OR.(ALAT(NY,NX).GT.0.0.AND.I.EQ.173) &
        .OR.(ALAT(NY,NX).LT.0.0.AND.I.EQ.355))THEN
        VRNF(NB,NZ,NY,NX)=0.0
        IFLGF(NB,NZ,NY,NX)=0
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
    IF(DYLN(NY,NX).LT.DYLX(NY,NX))THEN
      VRNF(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)
      IF(VRNF(NB,NZ,NY,NX).GE.VRNX(NB,NZ,NY,NX) &
        .OR.(ALAT(NY,NX).GT.0.0.AND.I.EQ.355) &
        .OR.(ALAT(NY,NX).LT.0.0.AND.I.EQ.173))THEN
        VRNS(NB,NZ,NY,NX)=0.0
        IFLGE(NB,NZ,NY,NX)=0
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
!   IDAY(2,=date of floral initiation
!
  ELSEIF(IWTYP(NZ,NY,NX).EQ.1)THEN
    IF((DYLN(NY,NX).GE.DYLX(NY,NX) &
      .OR.(DYLN(NY,NX).LT.DYLX(NY,NX) &
      .AND.VRNF(NB,NZ,NY,NX).LT.VRNX(NB,NZ,NY,NX))) &
      .AND.IFLGE(NB,NZ,NY,NX).EQ.0)THEN
        IF(TCG(NZ,NY,NX).GE.TCZ(NZ,NY,NX))THEN
          VRNS(NB,NZ,NY,NX)=VRNS(NB,NZ,NY,NX)+1.0
        ENDIF
        IF(VRNS(NB,NZ,NY,NX).LT.VRNL(NB,NZ,NY,NX))THEN
          IF(TCG(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
            VRNS(NB,NZ,NY,NX)=AMAX1(0.0,VRNS(NB,NZ,NY,NX)-1.0)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX) &
          .OR.(ALAT(NY,NX).GT.0.0.AND.I.EQ.173) &
          .OR.(ALAT(NY,NX).LT.0.0.AND.I.EQ.355))THEN
          VRNF(NB,NZ,NY,NX)=0.0
        ENDIF
      ENDIF
      IF(IDAY(2,NB,NZ,NY,NX).NE.0.OR.(DYLN(NY,NX).LT.DYLX(NY,NX) &
        .AND.DYLN(NY,NX).LT.12.0))THEN
        IFLGF(NB,NZ,NY,NX)=0
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
!     IDAY(2,=date of floral initiation
!
      IF(DYLN(NY,NX).LT.DYLX(NY,NX) &
        .AND.IFLGF(NB,NZ,NY,NX).EQ.0 &
        .AND.IDAY(2,NB,NZ,NY,NX).NE.0)THEN
        IF(TCG(NZ,NY,NX).LE.TCX(NZ,NY,NX))THEN
          VRNF(NB,NZ,NY,NX)=VRNF(NB,NZ,NY,NX)+1.0
        ENDIF
        IF(VRNF(NB,NZ,NY,NX).GE.VRNX(NB,NZ,NY,NX) &
          .AND.IFLGE(NB,NZ,NY,NX).EQ.1)THEN
          VRNS(NB,NZ,NY,NX)=0.0
          IFLGE(NB,NZ,NY,NX)=0
        ENDIF
      ENDIF
!     WRITE(*,4646)'VRNS',I,J,NZ,NB,IDAY(2,NB,NZ,NY,NX)
!    2,IFLGE(NB,NZ,NY,NX),IFLGF(NB,NZ,NY,NX),VRNS(NB,NZ,NY,NX)
!    2,TCG(NZ,NY,NX),TCZ(NZ,NY,NX),TCX(NZ,NY,NX),PSILG(NZ,NY,NX)
!    3,DYLN(NY,NX),DYLX(NY,NX),DYLM(NY,NX),VRNF(NB,NZ,NY,NX)
!    4,VRNL(NB,NZ,NY,NX),VRNX(NB,NZ,NY,NX)
!4646  FORMAT(A8,7I4,20E12.4)
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
!     IDAY(2,=date of floral initiation
!
    ELSEIF(IWTYP(NZ,NY,NX).EQ.2.OR.IWTYP(NZ,NY,NX).EQ.4 &
      .OR.IWTYP(NZ,NY,NX).EQ.5)THEN
      IF(IFLGE(NB,NZ,NY,NX).EQ.0)THEN
        IF(PSILT(NZ,NY,NX).GE.PSILX)THEN
          VRNS(NB,NZ,NY,NX)=VRNS(NB,NZ,NY,NX)+1.0
        ENDIF
        IF(VRNS(NB,NZ,NY,NX).LT.VRNL(NB,NZ,NY,NX))THEN
          IF(PSILT(NZ,NY,NX).LT.PSILY(IGTYP(NZ,NY,NX)))THEN
            VRNS(NB,NZ,NY,NX)=AMAX1(0.0,VRNS(NB,NZ,NY,NX)-12.0)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
          VRNF(NB,NZ,NY,NX)=0.0
          IF(IDAY(2,NB,NZ,NY,NX).NE.0)IFLGF(NB,NZ,NY,NX)=0
        ENDIF
      ENDIF
      IF(IDAY(2,NB,NZ,NY,NX).NE.0)IFLGF(NB,NZ,NY,NX)=0
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
!     IDAY(2,=date of floral initiation
!     VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!     VRNE=maximum hours for leafout,leafoff
!
      IF(IFLGE(NB,NZ,NY,NX).EQ.1 &
        .AND.IFLGF(NB,NZ,NY,NX).EQ.0)THEN
        IF(PSILT(NZ,NY,NX).LT.PSILY(IGTYP(NZ,NY,NX)))THEN
          VRNF(NB,NZ,NY,NX)=VRNF(NB,NZ,NY,NX)+1.0
        ENDIF
        IF(IWTYP(NZ,NY,NX).EQ.4)THEN
          IF(VRNZ(NB,NZ,NY,NX).GT.VRNE)THEN
            VRNF(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)
          ENDIF
        ELSEIF(IWTYP(NZ,NY,NX).EQ.5)THEN
          IF(VRNY(NB,NZ,NY,NX).GT.VRNE)THEN
            VRNF(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)
          ENDIF
        ENDIF
        IF(VRNF(NB,NZ,NY,NX).GE.VRNX(NB,NZ,NY,NX) &
          .AND.IFLGE(NB,NZ,NY,NX).EQ.1)THEN
          VRNS(NB,NZ,NY,NX)=0.0
          IFLGE(NB,NZ,NY,NX)=0
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
!     IDAY(2,=date of floral initiation
!
    ELSEIF(IWTYP(NZ,NY,NX).EQ.3)THEN
      IF((DYLN(NY,NX).GE.DYLX(NY,NX).OR.DYLN(NY,NX).GE.DYLM(NY,NX)-2.0) &
        .AND.IFLGE(NB,NZ,NY,NX).EQ.0)THEN
        IF(TCG(NZ,NY,NX).GE.TCZ(NZ,NY,NX).AND.PSILG(NZ,NY,NX).GT.PSILM)THEN
          VRNS(NB,NZ,NY,NX)=VRNS(NB,NZ,NY,NX)+1.0
        ENDIF
        IF(VRNS(NB,NZ,NY,NX).LT.VRNL(NB,NZ,NY,NX))THEN
          IF(TCG(NZ,NY,NX).LT.CTC(NZ,NY,NX) &
            .OR.PSILG(NZ,NY,NX).LT.PSILM)THEN
            VRNS(NB,NZ,NY,NX)=AMAX1(0.0,VRNS(NB,NZ,NY,NX)-1.5)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
          VRNF(NB,NZ,NY,NX)=0.0
          IF(IDAY(2,NB,NZ,NY,NX).NE.0)IFLGF(NB,NZ,NY,NX)=0
        ENDIF
      ENDIF
      IF(IDAY(2,NB,NZ,NY,NX).NE.0)IFLGF(NB,NZ,NY,NX)=0
!       WRITE(*,4647)'VRNS',I,J,NZ,NB,VRNS(NB,NZ,NY,NX),TCG(NZ,NY,NX)
!     2,TCZ(NZ,NY,NX),PSILG(NZ,NY,NX),PSILM,CTC(NZ,NY,NX)
!     3,DYLN(NY,NX),DYLX(NY,NX),DYLM(NY,NX),VRNL(NB,NZ,NY,NX)
!4647 FORMAT(A8,4I4,20E12.4)
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
!     IDAY(2,=date of floral initiation
!
      IF((DYLN(NY,NX).LT.DYLX(NY,NX).OR.DYLN(NY,NX) &
        .LT.24.0-DYLM(NY,NX)+2.0).AND.IFLGF(NB,NZ,NY,NX).EQ.0)THEN
        IF(TCG(NZ,NY,NX).LE.TCX(NZ,NY,NX) &
          .OR.PSILT(NZ,NY,NX).LT.PSILY(IGTYP(NZ,NY,NX)))THEN
          VRNF(NB,NZ,NY,NX)=VRNF(NB,NZ,NY,NX)+1.0
        ENDIF
        IF(VRNF(NB,NZ,NY,NX).GE.VRNX(NB,NZ,NY,NX) &
          .AND.IFLGE(NB,NZ,NY,NX).EQ.1)THEN
        VRNS(NB,NZ,NY,NX)=0.0
        IFLGE(NB,NZ,NY,NX)=0
      ENDIF
    ENDIF
  ENDIF
  end subroutine pft_specific_phenology
!------------------------------------------------------------------------------------------

  subroutine living_branch_phenology(I,J,NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NB,NZ,NY,NX

  IF(IDAY(1,NB,NZ,NY,NX).EQ.0)THEN
    IDAY(1,NB,NZ,NY,NX)=I
    IFLGA(NB,NZ,NY,NX)=1
    IFLGE(NB,NZ,NY,NX)=0
    VRNS(NB,NZ,NY,NX)=0.5*VRNS(NB1(NZ,NY,NX),NZ,NY,NX)
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
  IF(IWTYP(NZ,NY,NX).EQ.0 &
    .OR.VRNF(NB,NZ,NY,NX).LT.VRNX(NB,NZ,NY,NX))THEN
    TKCO=TKG(NZ,NY,NX)+OFFST(NZ,NY,NX)
    RTK=8.3143*TKCO
    STK=710.0*TKCO
    ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-218500)/RTK)
    TFNP=EXP(24.269-60000/RTK)/ACTV
    RNI=AMAX1(0.0,XRNI(NZ,NY,NX)*TFNP)
    RLA=AMAX1(0.0,XRLA(NZ,NY,NX)*TFNP)
!
!   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
!
!   PSILG=leaf turgor potential
!   WFNG=water stress effect on phenology
!
    IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      IF(IDAY(6,NB,NZ,NY,NX).EQ.0)THEN
        WFNG=EXP(0.025*PSILT(NZ,NY,NX))
        RNI=RNI*WFNG
        RLA=RLA*WFNG
      ENDIF
      IF(IDAY(2,NB,NZ,NY,NX).EQ.0)THEN
        OFNG=SQRT(OSTR(NZ,NY,NX))
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
    PSTG(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)+RNI
    VSTG(NB,NZ,NY,NX)=VSTG(NB,NZ,NY,NX)+RLA
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
    IF(IDAY(2,NB,NZ,NY,NX).NE.0)THEN
      GSTGI(NB,NZ,NY,NX)=(PSTG(NB,NZ,NY,NX)-PSTGI(NB,NZ,NY,NX))/GROUPI(NZ,NY,NX)
      DGSTGI(NB,NZ,NY,NX)=RNI/(GROUPI(NZ,NY,NX)*GSTGG)
      TGSTGI(NB,NZ,NY,NX)=TGSTGI(NB,NZ,NY,NX)+DGSTGI(NB,NZ,NY,NX)
    ENDIF
    IF(IDAY(6,NB,NZ,NY,NX).NE.0)THEN
      GSTGF(NB,NZ,NY,NX)=(PSTG(NB,NZ,NY,NX)-PSTGF(NB,NZ,NY,NX))/GROUPI(NZ,NY,NX)
      DGSTGF(NB,NZ,NY,NX)=RNI/(GROUPI(NZ,NY,NX)*GSTGR)
      TGSTGF(NB,NZ,NY,NX)=TGSTGF(NB,NZ,NY,NX)+DGSTGF(NB,NZ,NY,NX)
    ENDIF
    IFLGG(NB,NZ,NY,NX)=1
  ELSE
    IFLGG(NB,NZ,NY,NX)=0
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
  IF(IDAY(2,NB,NZ,NY,NX).EQ.0)THEN
    IF(PSTG(NB,NZ,NY,NX).GT.GROUP(NB,NZ,NY,NX)+PSTGI(NB,NZ,NY,NX) &
      .AND.((VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX)) &
      .OR.(I.GE.IDAY0(NZ,NY,NX).AND.IYRC.EQ.IYR0(NZ,NY,NX) &
      .AND.DYLN(NY,NX).GT.DYLX(NY,NX))) &
      .OR.(((ISTYP(NZ,NY,NX).EQ.1.AND.(IWTYP(NZ,NY,NX).EQ.1 &
      .OR.IWTYP(NZ,NY,NX).EQ.3)) &
      .OR.(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).EQ.0)) &
      .AND.ZC(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO &
      .AND.DYLN(NY,NX).LT.DYLX(NY,NX)))THEN
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
      IF(IPTYP(NZ,NY,NX).EQ.0)THEN
        PPD=0.0
      ELSE
        PPD=AMAX1(0.0,XDL(NZ,NY,NX)-DYLN(NY,NX))
        IF(IPTYP(NZ,NY,NX).EQ.1.AND.DYLN(NY,NX).GE.DYLX(NY,NX))PPD=0.0
      ENDIF
!     IF(NZ.EQ.1)THEN
!     WRITE(*,333)'IDAY2',I,J,NZ,NB,IDAY(2,NB,NZ,NY,NX),IDAY0(NZ,NY,NX)
!    2,IYR0(NZ,NY,NX),IPTYP(NZ,NY,NX)
!    2,PPD,XDL(NZ,NY,NX),DYLN(NY,NX),DYLX(NY,NX),VRNS(NB,NZ,NY,NX)
!    3,VRNL(NB,NZ,NY,NX),PSTG(NB,NZ,NY,NX),GROUP(NB,NZ,NY,NX)
!    4,PSTGI(NB,NZ,NY,NX),XPPD(NZ,NY,NX)
!333   FORMAT(A8,8I4,20E12.4)
!     ENDIF
      IF(IPTYP(NZ,NY,NX).EQ.0 &
        .OR.(IPTYP(NZ,NY,NX).EQ.1.AND.PPD.GT.XPPD(NZ,NY,NX)) &
        .OR.(IPTYP(NZ,NY,NX).EQ.2.AND.PPD.LT.XPPD(NZ,NY,NX)) &
        .OR.(((ISTYP(NZ,NY,NX).EQ.1.AND.(IWTYP(NZ,NY,NX).EQ.1 &
        .OR.IWTYP(NZ,NY,NX).EQ.3)) &
        .OR.(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).EQ.0)) &
        .AND.ZC(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO &
        .AND.DYLN(NY,NX).LT.DYLX(NY,NX)))THEN
        IDAY(2,NB,NZ,NY,NX)=I
        PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
        IF(ISTYP(NZ,NY,NX).EQ.0.AND.IDTYP(NZ,NY,NX).EQ.0)THEN
          VSTGX(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
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
  ELSEIF(IDAY(3,NB,NZ,NY,NX).EQ.0)THEN
    IF(GSTGI(NB,NZ,NY,NX).GT.0.25*GSTGG &
      .OR.((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3) &
      .AND.ISTYP(NZ,NY,NX).NE.0.AND.IPTYP(NZ,NY,NX).NE.1 &
      .AND.DYLN(NY,NX).LT.DYLX(NY,NX).AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX)) &
      .OR.(IWTYP(NZ,NY,NX).EQ.2.AND.ISTYP(NZ,NY,NX).EQ.0) &
      .AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX))THEN
      IDAY(3,NB,NZ,NY,NX)=I
    ENDIF
!
!   IDAY(4,=mid stem elongation
!
  ELSEIF(IDAY(4,NB,NZ,NY,NX).EQ.0)THEN
    IF(GSTGI(NB,NZ,NY,NX).GT.0.50*GSTGG &
      .OR.((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3) &
      .AND.ISTYP(NZ,NY,NX).NE.0.AND.IPTYP(NZ,NY,NX).NE.1 &
      .AND.DYLN(NY,NX).LT.DYLX(NY,NX).AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX)) &
      .OR.(IWTYP(NZ,NY,NX).EQ.2.AND.ISTYP(NZ,NY,NX).EQ.0) &
      .AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX))THEN
      IDAY(4,NB,NZ,NY,NX)=I
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IDTYP(NZ,NY,NX).NE.0)THEN
        VSTGX(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      ENDIF
    ENDIF
!
!   IDAY(5,=end of stem elongation and setting max seed number
!
  ELSEIF(IDAY(5,NB,NZ,NY,NX).EQ.0)THEN
    IF(GSTGI(NB,NZ,NY,NX).GT.1.00*GSTGG &
      .OR.((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3) &
      .AND.ISTYP(NZ,NY,NX).NE.0.AND.IPTYP(NZ,NY,NX).NE.1 &
      .AND.DYLN(NY,NX).LT.DYLX(NY,NX).AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX)) &
      .OR.(IWTYP(NZ,NY,NX).EQ.2.AND.ISTYP(NZ,NY,NX).EQ.0) &
      .AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX))THEN
      IDAY(5,NB,NZ,NY,NX)=I
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
  ELSEIF(IDAY(6,NB,NZ,NY,NX).EQ.0)THEN
    IF((VSTG(NB,NZ,NY,NX).GT.PSTGI(NB,NZ,NY,NX)) &
      .OR.(ISTYP(NZ,NY,NX).NE.0.AND.IDAY(5,NB,NZ,NY,NX).NE.0) &
      .OR.((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3) &
      .AND.ISTYP(NZ,NY,NX).NE.0.AND.IPTYP(NZ,NY,NX).NE.1 &
      .AND.DYLN(NY,NX).LT.DYLX(NY,NX).AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX)) &
      .OR.(IWTYP(NZ,NY,NX).EQ.2.AND.ISTYP(NZ,NY,NX).EQ.0) &
      .AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX))THEN
      IF(NB.EQ.NB1(NZ,NY,NX) &
        .OR.IDAY(6,NB1(NZ,NY,NX),NZ,NY,NX).NE.0)THEN
        IDAY(6,NB,NZ,NY,NX)=I
        PSTGF(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
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
  ELSEIF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
    IF(GSTGF(NB,NZ,NY,NX).GT.0.50*GSTGR &
      .OR.((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3) &
      .AND.ISTYP(NZ,NY,NX).NE.0.AND.IPTYP(NZ,NY,NX).NE.1 &
      .AND.DYLN(NY,NX).LT.DYLX(NY,NX).AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX)) &
      .OR.(IWTYP(NZ,NY,NX).EQ.2.AND.ISTYP(NZ,NY,NX).EQ.0) &
      .AND.IFLGE(NB,NZ,NY,NX).EQ.1 &
      .AND.VRNF(NB,NZ,NY,NX).GT.VRNX(NB,NZ,NY,NX))THEN
        IDAY(7,NB,NZ,NY,NX)=I
!       IF(IWTYP(NZ,NY,NX).NE.0.AND.NB.EQ.NB1(NZ,NY,NX))THEN
!       DO 1500 NBB=1,NBR(NZ,NY,NX)
!       IF(NBB.NE.NB.AND.IDAY(5,NBB,NZ,NY,NX).EQ.0)THEN
!       IDAY(5,NBB,NZ,NY,NX)=I
!       PSTGF(NBB,NZ,NY,NX)=PSTG(NBB,NZ,NY,NX)
!       ENDIF
!1500  CONTINUE
!       ENDIF
    ENDIF
!
!   END SEED NUMBER SET PERIOD
!
!   IDAY(8,=end date setting for final seed number
!
  ELSEIF(IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
    IF(GSTGF(NB,NZ,NY,NX).GT.1.00*GSTGR)THEN
      IDAY(8,NB,NZ,NY,NX)=I
!     IF(IWTYP(NZ,NY,NX).NE.0.AND.NB.EQ.NB1(NZ,NY,NX))THEN
!     DO 1495 NBB=1,NBR(NZ,NY,NX)
!     IF(NBB.NE.NB.AND.IDAY(6,NBB,NZ,NY,NX).EQ.0)THEN
!     IDAY(6,NBB,NZ,NY,NX)=I
!     ENDIF
!1495  CONTINUE
!     ENDIF
    ENDIF
!
!   END SEED SIZE SET PERIOD
!
!   IDAY(9,=end of setting max seed size
!
  ELSEIF(IDAY(9,NB,NZ,NY,NX).EQ.0)THEN
    IF(GSTGF(NB,NZ,NY,NX).GT.1.50*GSTGR)THEN
      IDAY(9,NB,NZ,NY,NX)=I
    ENDIF
  ENDIF
  end subroutine living_branch_phenology

end module HfuncMod
