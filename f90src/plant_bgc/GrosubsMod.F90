module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : test_aeqb,safe_adb
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  implicit none

  private

! end_include_section

  character(len=*), private, parameter :: mod_filename = __FILE__
! DIMENSION VCO2(400,366,05)
!
!
!     RTSK=relative primary root sink strength 0.25=shallow,4.0=deep root profile
!     FXRN=rate constant for plant-bacteria nonstructl C,N,P exchange (h-1)
!     FXFB=rate constant for leaf-storage nonstructl C,N,P exchange (h-1)
!     FXFR=rate constant for root-storage nonstructl C,N,P exchange (h-1)
!     FPART=parameter for allocating nonstructural C to shoot organs
!     FXSH,FXRT=shoot-root partitioning of storage C during leafout
!     FRSV=rate constant for remobiln of storage C,N,P during leafout (h-1)
!     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P exchange (h-1)
!     EFIRE=combustion  of N,P relative to C
!     PSILY=canopy water potential below which leafoff is induced (MPa)
!     FLG4Y=number of hours after physiol maturity required for senescence
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!     FVRN=fraction of hours required for leafoff to initiate remobilization
!

  integer,private :: curday,curhour
  public :: grosubs
  public :: InitGrosub
  contains

  subroutine InitGrosub

  implicit none

  call InitVegPars


  end subroutine InitGrosub
!------------------------------------------------------------------------------------------

  subroutine grosubs(I,J)
!
!     THIS subroutine CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  use PlantDisturbsMod, only : RemoveBiomassByDisturbance
  implicit none
  integer, intent(in) :: I, J

  real(r8) :: ZCX(JP1)
  integer :: L,K,M
  integer :: NZ
  real(r8) :: CPOOLK(JC1,JP1)
! begin_execution
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
!

      DO 9980 NZ=1,NP0s1
        DO 1 L=0,NJs1
          DO K=0,1
            DO M=1,4
              CSNCs1(M,K,L,NZ)=0._r8
              ZSNCs1(M,K,L,NZ)=0._r8
              PSNCs1(M,K,L,NZ)=0._r8
            ENDDO
          ENDDO
1       CONTINUE
        HCSNCs1(NZ)=0._r8
        HZSNCs1(NZ)=0._r8
        HPSNCs1(NZ)=0._r8
        CNETs1(NZ)=0._r8
        ZCX(NZ)=ZCs1(NZ)
        ZCs1(NZ)=0._r8
9980  CONTINUE
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
      DO 9985 NZ=1,NPs1

! IFLGC= flag for living pft
        IF(IFLGCs1(NZ).EQ.1)THEN
          call GrowPlant(I,J,NZ,ZCX,CPOOLK)
        ENDIF

!     HARVEST STANDING DEAD

        call RemoveBiomassByDisturbance(I,J,NZ,CPOOLK)
9985  CONTINUE
!
!     TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS
      call LiveDeadTransformation(I,J)
  END subroutine grosubs

!------------------------------------------------------------------------------------------

  subroutine LiveDeadTransformation(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: L,K,NZ,M,NB
  real(r8) :: XFRC,XFRN,XFRP
!     begin_execution

  DO 9975 NZ=1,NP0s1
!
!     ACTIVATE DORMANT SEEDS
!
    DO 205 NB=1,NBRs1(NZ)
      IF(IFLGIs1(NZ).EQ.1)THEN
        IF(IFLGEs1(NB,NZ).EQ.0 &
          .AND.VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ))THEN
          IDAY0s1(NZ)=I
          IYR0s1(NZ)=IYRCs1
          SDPTHIs1(NZ)=0.005_r8+CDPTHZs1(0)
          IFLGIs1(NZ)=0
        ENDIF
      ENDIF
205 CONTINUE
!
!     LITTERFALL FROM STANDING DEAD
!
!     XFRC,XFRN,XFRP=litterfall from standing dead
!     TFN3=temperature function for canopy growth
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall
!
    DO 6235 M=1,4
      XFRC=1.5814E-05*TFN3s1(NZ)*WTSTDGs1(M,NZ)
      XFRN=1.5814E-05*TFN3s1(NZ)*WTSTDNs1(M,NZ)
      XFRP=1.5814E-05*TFN3s1(NZ)*WTSTDPs1(M,NZ)
      IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+XFRC
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+XFRN
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+XFRP
      ELSE
        CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+XFRC
        ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+XFRN
        PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+XFRP
      ENDIF
      WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ)-XFRC
      WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ)-XFRN
      WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ)-XFRP
6235  CONTINUE
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LITTERFALL
!
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
!
    DO 6430 M=1,4
      DO K=0,1
        TCSN0s1(NZ)=TCSN0s1(NZ)+CSNCs1(M,K,0,NZ)
        TZSN0s1(NZ)=TZSN0s1(NZ)+ZSNCs1(M,K,0,NZ)
        TPSN0s1(NZ)=TPSN0s1(NZ)+PSNCs1(M,K,0,NZ)
        DO 8955 L=0,NJs1
          HCSNCs1(NZ)=HCSNCs1(NZ)+CSNCs1(M,K,L,NZ)
          HZSNCs1(NZ)=HZSNCs1(NZ)+ZSNCs1(M,K,L,NZ)
          HPSNCs1(NZ)=HPSNCs1(NZ)+PSNCs1(M,K,L,NZ)
          TCSNCs1(NZ)=TCSNCs1(NZ)+CSNCs1(M,K,L,NZ)
          TZSNCs1(NZ)=TZSNCs1(NZ)+ZSNCs1(M,K,L,NZ)
          TPSNCs1(NZ)=TPSNCs1(NZ)+PSNCs1(M,K,L,NZ)
8955    CONTINUE
      enddo
6430  CONTINUE
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    WTSTGs1(NZ)=WTSTDGs1(1,NZ)+WTSTDGs1(2,NZ) &
      +WTSTDGs1(3,NZ)+WTSTDGs1(4,NZ)
    WTSTGNs1(NZ)=WTSTDNs1(1,NZ)+WTSTDNs1(2,NZ) &
      +WTSTDNs1(3,NZ)+WTSTDNs1(4,NZ)
    WTSTGPs1(NZ)=WTSTDPs1(1,NZ)+WTSTDPs1(2,NZ) &
      +WTSTDPs1(3,NZ)+WTSTDPs1(4,NZ)
!
!     PLANT C BALANCE = TOTAL C STATE VARIABLES + TOTAL
!     AUTOTROPHIC RESPIRATION + TOTAL LITTERFALL - TOTAL EXUDATION
!     - TOTAL CO2 FIXATION
!
!     BALC=PFT C balance
!     WTSHT,WTRT,WTND,WTRVC,WTSTG=PFT shoot,root,bacteria,storage,standing dead C
!     ZNPP=cumulative PFT NPP
!     TCSNC=cumulative PFT C litterfall
!     TCUPTK=cumulative PFT root-soil C exchange
!     RSETC=cumulative C balance from previous year
!     THVSTC=cumulative PFT C removed from ecosystem from previous year
!     HVSTC=total PFT C removed from ecosystem in current year
!     VCO2F,VCH4F=CO2,CH4 emission from disturbance
!
    ZNPPs1(NZ)=CARBNs1(NZ)+TCO2Ts1(NZ)
    IF(IFLGCs1(NZ).EQ.1)THEN
      BALCs1(NZ)=WTSHTs1(NZ)+WTRTs1(NZ)+WTNDs1(NZ) &
        +WTRVCs1(NZ)-ZNPPs1(NZ)+TCSNCs1(NZ)-TCUPTKs1(NZ) &
        -RSETCs1(NZ)+WTSTGs1(NZ)+THVSTCs1(NZ) &
        +HVSTCs1(NZ)-VCO2Fs1(NZ)-VCH4Fs1(NZ)
!
!     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LITTERFALL
!     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
!
!     BALN=PFT N balance
!     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT shoot,root,bacteria,storage,standing dead N
!     TZSNC=cumulative PFT N litterfall
!     TZUPTK=cumulative PFT root-soil N exchange
!     TNH3C=cumulative NH3 exchange
!     RSETN=cumulative N balance from previous year
!     THVSTN=cumulative PFT N removed from ecosystem from previous year
!     HVSTN=total PFT N removed from ecosystem in current year
!     VNH3F,VN2OF=NH3,N2O emission from disturbance
!     TZUPFX=cumulative PFT N2 fixation
!
      BALNs1(NZ)=WTSHNs1(NZ)+WTRTNs1(NZ)+WTNDNs1(NZ) &
        +WTRVNs1(NZ)+TZSNCs1(NZ)-TZUPTKs1(NZ)-TNH3Cs1(NZ) &
        -RSETNs1(NZ)+WTSTGNs1(NZ)+HVSTNs1(NZ)+THVSTNs1(NZ) &
        -VNH3Fs1(NZ)-VN2OFs1(NZ)-TZUPFXs1(NZ)
!
!     PLANT P BALANCE = TOTAL P STATE VARIABLES + TOTAL P LITTERFALL
!     - TOTAL P UPTAKE FROM SOIL
!
!     BALP=PFT N balance
!     WTSHP,WTRTP,WTNDP,WTRVP,WTSTGP=PFT shoot,root,bacteria,storage,standing dead P
!     TPSNC=cumulative PFT P litterfall
!     TPUPTK=cumulative PFT root-soil P exchange
!     RSETP=cumulative P balance from previous year
!     THVSTP=cumulative PFT P removed from ecosystem from previous year
!     HVSTP=total PFT P removed from ecosystem in current year
!     VPO4F=PO4 emission from disturbance
!
      BALPs1(NZ)=WTSHPs1(NZ)+WTRTPs1(NZ)+WTNDPs1(NZ) &
        +WTRVPs1(NZ)+TPSNCs1(NZ)-TPUPTKs1(NZ) &
        -RSETPs1(NZ)+WTSTDPs1(1,NZ)+WTSTGPs1(NZ) &
        +HVSTPs1(NZ)+THVSTPs1(NZ)-VPO4Fs1(NZ)
!     IF(NZ.EQ.4)THEN
!     WRITE(*,1112)'BALP',I,J,NX,NY,NZ,BALPs1(NZ),WTSHPs1(NZ)
!    2,WTRTPs1(NZ),WTNDPs1(NZ),WTRVPs1(NZ),TPSNCs1(NZ)
!    3,TPUPTKs1(NZ),RSETPs1(NZ)
!    4,WTSTDPs1(1,NZ),WTSTGPs1(NZ),HVSTPs1(NZ)
!    5,THVSTPs1(NZ),VPO4Fs1(NZ)
!     ENDIF
    ENDIF
9975  CONTINUE
  end subroutine LiveDeadTransformation
!------------------------------------------------------------------------------------------

  subroutine GrowPlant(I,J,NZ,ZCX,CPOOLK)
  use PlantDisturbsMod, only : RemoveBiomByManagement
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: ZCX(JP1)
  real(r8), intent(out) :: CPOOLK(JC1,JP1)

  real(r8)  :: UPNFC(JP1)
  integer  :: ICHK1(2,JZ1),IDTHRN,NB
  integer  :: NRX(2,JZ1),IFLGZ
  real(r8) :: TFN6(JZ1)
  real(r8) :: WFNGR(2,JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT
  real(r8) :: RTVL
  real(r8) :: RLNT(2,JZ1)
  real(r8) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  real(r8) :: RTNT(2)
  real(r8) :: XRTN1
  real(r8) :: TFN5
  real(r8) :: WFNG
  real(r8) :: WFNC
  real(r8) :: WFNS,WFNSG
! begin_execution

  IF(IDTHPs1(NZ).EQ.0.OR.IDTHRs1(NZ).EQ.0)THEN
    UPNFC(NZ)=0._r8
    IFLGZ = 0
    call StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,WFNC,WFNS,WFNSG)
!
!     CALCULATE GROWTH OF EACH BRANCH
!
!     WTLFB,WTSHEB,WTLSB=leaf,petiole,leaf+petiole mass
!     IDTHB=branch living flag: 0=alive,1=dead
!
    DO 105 NB=1,NBRs1(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6,ZCX,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WFNG,WFNC,WFNS,WFNSG,PTRT,UPNFC,IFLGZ)
105 CONTINUE
!
!     ROOT GROWTH
!
    call RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,XRTN1,WFNGR,RLNT,RTSK1,RTSK2,RTNT)
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
    RTLGPs1(1,NGs1(NZ),NZ)=RTLGPs1(1,NGs1(NZ),NZ)+SDLGs1(NZ)
    IF(DLYR3s1(NGs1(NZ)).GT.ZEROs1)THEN
      RTDNPs1(1,NGs1(NZ),NZ)=RTLGPs1(1,NGs1(NZ),NZ)/DLYR3s1(NGs1(NZ))
    ELSE
      RTDNPs1(1,NGs1(NZ),NZ)=0._r8
    ENDIF
    RTVL=RTVLPs1(1,NGs1(NZ),NZ)+RTVLWs1(1,NGs1(NZ),NZ)+SDVLs1(NZ)*PPs1(NZ)
    RTVLPs1(1,NGs1(NZ),NZ)=PORTs1(1,NZ)*RTVL
    RTVLWs1(1,NGs1(NZ),NZ)=(1.0-PORTs1(1,NZ))*RTVL
    RTARPs1(1,NGs1(NZ),NZ)=RTARPs1(1,NGs1(NZ),NZ)+SDARs1(NZ)
    IF(IDTHRN.EQ.NRTs1(NZ).OR.(WTRVCs1(NZ).LE.ZEROLs1(NZ).AND.ISTYPs1(NZ).NE.0))THEN
      IDTHRs1(NZ)=1
      IDTHPs1(NZ)=1
    ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
!
    call RootNoduleBiomchemistry(I,J,NZ,TFN6,WFNGR)

    call NonstructlBiomTransfer(I,J,NZ,PTRT,RLNT,RTSK1,RTSK2,RTNT,IFLGZ)
!
    call ComputeTotalBiom(NZ,CPOOLK)
  ELSE
    HCUPTKs1(NZ)=UPOMCs1(NZ)
    HZUPTKs1(NZ)=UPOMNs1(NZ)+UPNH4s1(NZ)+UPNO3s1(NZ)+UPNFs1(NZ)
    HPUPTKs1(NZ)=UPOMPs1(NZ)+UPH2Ps1(NZ)+UPH1Ps1(NZ)
  ENDIF
!
  call RemoveBiomByManagement(I,J,NZ,CPOOLK)
!
!     RESET DEAD BRANCHES
  call ResetDeadBranch(I,J,NZ,CPOOLK)
!
  call AccumulateStates(I,J,NZ,UPNFC)

  end subroutine GrowPlant
!------------------------------------------------------------------------------------------

  subroutine ResetDeadBranch(I,J,NZ,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8),intent(inout) :: CPOOLK(JC1,JP1)
  integer :: IDTHY

!     begin_execution
!
!     ZNOON=hour of solar noon
!     IDAYs1(1,=emergence date
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAYH,IYRH=day,year of harvesting
!     IYRCs1=current year
!     IDTHB=branch living flag: 0=alive,1=dead
!     GROUP=node number required for floral initiation
!     PSTGI=node number at floral initiation
!     PSTGF=node number at flowering
!     VSTG=number of leaves appeared
!     KVSTG=integer of most recent leaf number currently growing
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     IFLGA=flag for initializing leafout
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     ATRP=hourly leafout counter
!     FDBK,FDBKX=N,P feedback inhibition on C3 CO2 fixation
!     IFLGA,IFLGE=flag for initializing,enabling leafout
!     IFLGF=flag for enabling leafoff:0=enable,1=disable
!     IFLGQ=current hours after physl maturity until start of litterfall
!
  IF(J.EQ.INT(ZNOONs1) &
    .AND.IDAYs1(1,NB1s1(NZ),NZ).NE.0 &
    .AND.(ISTYPs1(NZ).NE.0.OR.(I.GE.IDAYHs1(NZ) &
    .AND.IYRCs1.GE.IYRHs1(NZ))))THEN
    IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)

    IF(IDTHY.EQ.NBRs1(NZ))THEN
      IDTHPs1(NZ)=1
      NBTs1(NZ)=0
      WSTRs1(NZ)=0._r8
      IF(IFLGIs1(NZ).EQ.1)THEN
        NBRs1(NZ)=1
      ELSE
        NBRs1(NZ)=0
      ENDIF
      HTCTLs1(NZ)=0._r8
      VOLWOUs1=VOLWOUs1+VOLWPs1(NZ)
      UVOLOs1=UVOLOs1+VOLWPs1(NZ)
      VOLWPs1(NZ)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     ISTYP=growth habit:0=annual,1=perennial
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(WTRVCs1(NZ).LT.1.0E-04*WTRTs1(NZ) &
        .AND.ISTYPs1(NZ).NE.0)IDTHRs1(NZ)=1
      IF(ISTYPs1(NZ).EQ.0)IDTHRs1(NZ)=1
      IF(JHVSTs1(NZ).NE.0)IDTHRs1(NZ)=1
      IF(PPs1(NZ).LE.0.0)IDTHRs1(NZ)=1
      IF(IDTHRs1(NZ).EQ.1)IDTHPs1(NZ)=1
    ENDIF
!
!     DEAD ROOTS
!
!
!     LITTERFALL FROM DEAD ROOTS
!
    call LiterfallFromDeadRoots(I,J,NZ)
!
    call LiterfallFromRootShootStorage(I,J,NZ,CPOOLK)
  ENDIF
  end subroutine ResetDeadBranch
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromRootShootStorage(I,J,NZ,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: CPOOLK(JC1,JP1)
  integer :: L,M,NR,NB,N
!     begin_execution
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM SHOOT AT DEATH
!
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!     IFLGI=PFT initialization flag:0=no,1=yes
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
  IF(IDTHPs1(NZ).EQ.1.AND.IDTHRs1(NZ).EQ.1)THEN
    IF(IFLGIs1(NZ).EQ.0)THEN
      DO 6425 M=1,4
        DO 8825 NB=1,NBRs1(NZ)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(0,M,NZ)*(CPOOLs1(NB,NZ)+CPOLNBs1(NB,NZ) &
            +CPOOLK(NB,NZ)+WTRSVBs1(NB,NZ)) &
            +CFOPCs1(1,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(1) &
            +WTNDBs1(NB,NZ)) &
            +CFOPCs1(2,M,NZ)*(WTSHEBs1(NB,NZ)*FWODBs1(1) &
            +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ))
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ) &
            +CFOPCs1(5,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(0) &
            +WTSHEBs1(NB,NZ)*FWODBs1(0))
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(0,M,NZ)*(ZPOOLs1(NB,NZ)+ZPOLNBs1(NB,NZ) &
            +WTRSBNs1(NB,NZ)) &
            +CFOPNs1(1,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(1) &
            +WTNDBNs1(NB,NZ)) &
            +CFOPNs1(2,M,NZ)*(WTSHBNs1(NB,NZ)*FWODSNs1(1) &
            +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ))
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ) &
            +CFOPNs1(5,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(0) &
            +WTSHBNs1(NB,NZ)*FWODSNs1(0))
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(0,M,NZ)*(PPOOLs1(NB,NZ)+PPOLNBs1(NB,NZ) &
            +WTRSBPs1(NB,NZ)) &
            +CFOPPs1(1,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(1) &
            +WTNDBPs1(NB,NZ)) &
            +CFOPPs1(2,M,NZ)*(WTSHBPs1(NB,NZ)*FWODSPs1(1) &
            +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ))
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ) &
            +CFOPPs1(5,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(0) &
            +WTSHBPs1(NB,NZ)*FWODSPs1(0))
          IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
            WTRVCs1(NZ)=WTRVCs1(NZ) &
              +CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
            WTRVNs1(NZ)=WTRVNs1(NZ) &
              +CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
            WTRVPs1(NZ)=WTRVPs1(NZ) &
              +CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
          ELSE
            CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
              +CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
              +CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
            PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
              +CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
          ENDIF
          IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
            CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
              +CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
              +CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
            PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
              +CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
          ELSE
            WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ) &
              +CFOPCs1(5,M,NZ)*WTSTKBs1(NB,NZ)
            WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ) &
              +CFOPNs1(5,M,NZ)*WTSTBNs1(NB,NZ)
            WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ) &
              +CFOPPs1(5,M,NZ)*WTSTBPs1(NB,NZ)
          ENDIF
8825    CONTINUE
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM ROOT AND STORGE AT DEATH
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        DO 6415 L=NUs1,NJs1
          DO N=1,MYs1(NZ)
            CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(0,M,NZ) &
              *CPOOLRs1(N,L,NZ)
            ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(0,M,NZ) &
              *ZPOOLRs1(N,L,NZ)
            PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(0,M,NZ) &
              *PPOOLRs1(N,L,NZ)
            DO NR=1,NRTs1(NZ)
              CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
                *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(0)
              ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
                *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(0)
              PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
                *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(0)
              CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
                *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(1)
              ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
                *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(1)
              PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
                *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(1)
            ENDDO
          ENDDO
6415    CONTINUE
        CSNCs1(M,0,NGs1(NZ),NZ)=CSNCs1(M,0,NGs1(NZ),NZ) &
          +CFOPCs1(0,M,NZ)*WTRVCs1(NZ)*FWOODs1(0)
        ZSNCs1(M,0,NGs1(NZ),NZ)=ZSNCs1(M,0,NGs1(NZ),NZ) &
          +CFOPNs1(0,M,NZ)*WTRVNs1(NZ)*FWOODNs1(0)
        PSNCs1(M,0,NGs1(NZ),NZ)=PSNCs1(M,0,NGs1(NZ),NZ) &
          +CFOPPs1(0,M,NZ)*WTRVPs1(NZ)*FWOODPs1(0)
        CSNCs1(M,1,NGs1(NZ),NZ)=CSNCs1(M,1,NGs1(NZ),NZ) &
          +CFOPCs1(0,M,NZ)*WTRVCs1(NZ)*FWOODs1(1)
        ZSNCs1(M,1,NGs1(NZ),NZ)=ZSNCs1(M,1,NGs1(NZ),NZ) &
          +CFOPNs1(0,M,NZ)*WTRVNs1(NZ)*FWOODNs1(1)
        PSNCs1(M,1,NGs1(NZ),NZ)=PSNCs1(M,1,NGs1(NZ),NZ) &
          +CFOPPs1(0,M,NZ)*WTRVPs1(NZ)*FWOODPs1(1)
6425  CONTINUE
!
      call ResetBranchRootStates(NZ,CPOOLK)
    ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     LYRC=number of days in current year
!     IDAY0,IYR0=day,year of planting
!
    IF(ISTYPs1(NZ).NE.0.AND.JHVSTs1(NZ).EQ.0)THEN
      IF(I.LT.LYRCs1)THEN
        IDAY0s1(NZ)=I+1
        IYR0s1(NZ)=IDATAs1(3)
      ELSE
        IDAY0s1(NZ)=1
        IYR0s1(NZ)=IDATAs1(3)+1
      ENDIF
    ENDIF
  ENDIF
  end subroutine LiterfallFromRootShootStorage
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadRoots(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,M,NR,N
!     begin_execution
!     IDTHR=PFT root living flag: 0=alive,1=dead
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
  IF(IDTHRs1(NZ).EQ.1)THEN
    DO 8900 N=1,MYs1(NZ)
      DO 8895 L=NUs1,NJs1
        DO 6410 M=1,4
          CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(0,M,NZ)*CPOOLRs1(N,L,NZ)
          ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(0,M,NZ)*ZPOOLRs1(N,L,NZ)
          PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(0,M,NZ)*PPOOLRs1(N,L,NZ)
          DO  NR=1,NRTs1(NZ)
            CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
              *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(0)
            ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
              *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(0)
            PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
              *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(0)
            CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
              *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(1)
            ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
              *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(1)
            PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
              *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(1)
          enddo
6410    CONTINUE
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        RCO2Zs1(NZ)=RCO2Zs1(NZ)-CO2As1(N,L,NZ)-CO2Ps1(N,L,NZ)
        ROXYZs1(NZ)=ROXYZs1(NZ)-OXYAs1(N,L,NZ)-OXYPs1(N,L,NZ)
        RCH4Zs1(NZ)=RCH4Zs1(NZ)-CH4As1(N,L,NZ)-CH4Ps1(N,L,NZ)
        RN2OZs1(NZ)=RN2OZs1(NZ)-Z2OAs1(N,L,NZ)-Z2OPs1(N,L,NZ)
        RNH3Zs1(NZ)=RNH3Zs1(NZ)-ZH3As1(N,L,NZ)-ZH3Ps1(N,L,NZ)
        RH2GZs1(NZ)=RH2GZs1(NZ)-H2GAs1(N,L,NZ)-H2GPs1(N,L,NZ)
        CO2As1(N,L,NZ)=0._r8
        OXYAs1(N,L,NZ)=0._r8
        CH4As1(N,L,NZ)=0._r8
        Z2OAs1(N,L,NZ)=0._r8
        ZH3As1(N,L,NZ)=0._r8
        H2GAs1(N,L,NZ)=0._r8
        CO2Ps1(N,L,NZ)=0._r8
        OXYPs1(N,L,NZ)=0._r8
        CH4Ps1(N,L,NZ)=0._r8
        Z2OPs1(N,L,NZ)=0._r8
        ZH3Ps1(N,L,NZ)=0._r8
        H2GPs1(N,L,NZ)=0._r8
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RTLGA=average secondary root length
!
        DO 8870 NR=1,NRTs1(NZ)
          WTRT1s1(N,L,NR,NZ)=0._r8
          WTRT1Ns1(N,L,NR,NZ)=0._r8
          WTRT1Ps1(N,L,NR,NZ)=0._r8
          WTRT2s1(N,L,NR,NZ)=0._r8
          WTRT2Ns1(N,L,NR,NZ)=0._r8
          WTRT2Ps1(N,L,NR,NZ)=0._r8
          RTWT1s1(N,NR,NZ)=0._r8
          RTWT1Ns1(N,NR,NZ)=0._r8
          RTWT1Ps1(N,NR,NZ)=0._r8
          RTLG1s1(N,L,NR,NZ)=0._r8
          RTLG2s1(N,L,NR,NZ)=0._r8
          RTN2s1(N,L,NR,NZ)=0._r8
8870    CONTINUE
        CPOOLRs1(N,L,NZ)=0._r8
        ZPOOLRs1(N,L,NZ)=0._r8
        PPOOLRs1(N,L,NZ)=0._r8
        WTRTLs1(N,L,NZ)=0._r8
        WTRTDs1(N,L,NZ)=0._r8
        WSRTLs1(N,L,NZ)=0._r8
        RTN1s1(N,L,NZ)=0._r8
        RTNLs1(N,L,NZ)=0._r8
        RTLGPs1(N,L,NZ)=0._r8
        RTDNPs1(N,L,NZ)=0._r8
        RTVLPs1(N,L,NZ)=0._r8
        RTVLWs1(N,L,NZ)=0._r8
        RRAD1s1(N,L,NZ)=RRAD1Ms1(N,NZ)
        RRAD2s1(N,L,NZ)=RRAD2Ms1(N,NZ)
        RTARPs1(N,L,NZ)=0._r8
        RTLGAs1(N,L,NZ)=RTLGAX
!
!     LITTERFALL AND STATE VARIABLES FROM DEAD NODULES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(INTYPs1(NZ).NE.0.AND.N.EQ.1)THEN
          DO 6420 M=1,4
            CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
              *WTNDLs1(L,NZ)+CFOPCs1(0,M,NZ)*CPOOLNs1(L,NZ)
            ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
              *WTNDLNs1(L,NZ)+CFOPNs1(0,M,NZ)*ZPOOLNs1(L,NZ)
            PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
              *WTNDLPs1(L,NZ)+CFOPPs1(0,M,NZ)*PPOOLNs1(L,NZ)
6420      CONTINUE
          WTNDLs1(L,NZ)=0._r8
          WTNDLNs1(L,NZ)=0._r8
          WTNDLPs1(L,NZ)=0._r8
          CPOOLNs1(L,NZ)=0._r8
          ZPOOLNs1(L,NZ)=0._r8
          PPOOLNs1(L,NZ)=0._r8
        ENDIF
8895  CONTINUE
8900  CONTINUE
!
!     RESET DEPTH VARIABLES OF DEAD ROOTS
!
!     NINR=deepest root layer
!     RTDP1=primary root depth from soil surface
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!
    DO 8795 NR=1,NRTs1(NZ)
      NINRs1(NR,NZ)=NGs1(NZ)
      DO 8790 N=1,MYs1(NZ)
        RTDP1s1(N,NR,NZ)=SDPTHs1(NZ)
        RTWT1s1(N,NR,NZ)=0._r8
        RTWT1Ns1(N,NR,NZ)=0._r8
        RTWT1Ps1(N,NR,NZ)=0._r8
8790  CONTINUE
8795  CONTINUE
    NIXs1(NZ)=NGs1(NZ)
    NRTs1(NZ)=0
  ENDIF
  end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: IDTHY
  real(r8), intent(inout) :: CPOOLK(JC1,JP1)
  integer :: M,NB
!     begin_execution

  DO 8845 NB=1,NBRs1(NZ)
    IF(IDTHBs1(NB,NZ).EQ.1)THEN
      GROUPs1(NB,NZ)=GROUPIs1(NZ)
      PSTGs1(NB,NZ)=XTLIs1(NZ)
      PSTGIs1(NB,NZ)=PSTGs1(NB,NZ)
      PSTGFs1(NB,NZ)=0._r8
      VSTGs1(NB,NZ)=0._r8
      VSTGXs1(NB,NZ)=0._r8
      KLEAFs1(NB,NZ)=1
      KVSTGs1(NB,NZ)=1
      TGSTGIs1(NB,NZ)=0._r8
      TGSTGFs1(NB,NZ)=0._r8
      VRNSs1(NB,NZ)=0._r8
      VRNFs1(NB,NZ)=0._r8
      VRNYs1(NB,NZ)=0._r8
      VRNZs1(NB,NZ)=0._r8
      ATRPs1(NB,NZ)=0._r8
      FLG4s1(NB,NZ)=0._r8
      FDBKs1(NB,NZ)=1.0_r8
      FDBKXs1(NB,NZ)=1.0_r8
      IFLGAs1(NB,NZ)=0
      IFLGEs1(NB,NZ)=1
      IFLGFs1(NB,NZ)=0
      IFLGRs1(NB,NZ)=0
      IFLGQs1(NB,NZ)=0
      NBTBs1(NB,NZ)=0
      DO 8850 M=1,10
        IDAYs1(M,NB,NZ)=0
8850  CONTINUE
!
!     LITTERFALL FROM DEAD BRANCHES
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      DO 6405 M=1,4
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
          +CFOPCs1(0,M,NZ)*CPOLNBs1(NB,NZ) &
          +CFOPCs1(1,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(1) &
          +WTNDBs1(NB,NZ)) &
          +CFOPCs1(2,M,NZ)*(WTSHEBs1(NB,NZ)*FWODBs1(1) &
          +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ))
        CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ) &
          +CFOPCs1(5,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(0) &
          +WTSHEBs1(NB,NZ)*FWODBs1(0))
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
          +CFOPNs1(0,M,NZ)*ZPOLNBs1(NB,NZ) &
          +CFOPNs1(1,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(1) &
          +WTNDBNs1(NB,NZ)) &
          +CFOPNs1(2,M,NZ)*(WTSHBNs1(NB,NZ)*FWODSNs1(1) &
          +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ))
        ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ) &
          +CFOPNs1(5,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(0) &
          +WTSHBNs1(NB,NZ)*FWODSNs1(0))
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
          +CFOPPs1(0,M,NZ)*PPOLNBs1(NB,NZ) &
          +CFOPPs1(1,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(1) &
          +WTNDBPs1(NB,NZ)) &
          +CFOPPs1(2,M,NZ)*(WTSHBPs1(NB,NZ)*FWODSPs1(1) &
          +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ))
        PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ) &
          +CFOPPs1(5,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(0) &
          +WTSHBPs1(NB,NZ)*FWODSPs1(0))
        IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
          WTRVCs1(NZ)=WTRVCs1(NZ) &
            +CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
          WTRVNs1(NZ)=WTRVNs1(NZ) &
            +CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
          WTRVPs1(NZ)=WTRVPs1(NZ) &
            +CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
        ELSE
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
        ENDIF
        IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
        ELSE
          WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ) &
            +CFOPCs1(5,M,NZ)*WTSTKBs1(NB,NZ)
          WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ) &
            +CFOPNs1(5,M,NZ)*WTSTBNs1(NB,NZ)
          WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ) &
            +CFOPPs1(5,M,NZ)*WTSTBPs1(NB,NZ)
        ENDIF
6405  CONTINUE
!
!     RECOVER NON-STRUCTURAL C,N,P FROM BRANCH TO
!     SEASONAL STORAGE RESERVES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
      WTRVCs1(NZ)=WTRVCs1(NZ)+CPOOLs1(NB,NZ)+CPOOLK(NB,NZ)
      WTRVNs1(NZ)=WTRVNs1(NZ)+ZPOOLs1(NB,NZ)
      WTRVPs1(NZ)=WTRVPs1(NZ)+PPOOLs1(NB,NZ)
      IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
        DO 6406 M=1,4
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(0,M,NZ)*WTRSVBs1(NB,NZ)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(0,M,NZ)*WTRSBNs1(NB,NZ)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(0,M,NZ)*WTRSBPs1(NB,NZ)
6406    CONTINUE
      ELSE
        WTRVCs1(NZ)=WTRVCs1(NZ)+WTRSVBs1(NB,NZ)
        WTRVNs1(NZ)=WTRVNs1(NZ)+WTRSBNs1(NB,NZ)
        WTRVPs1(NZ)=WTRVPs1(NZ)+WTRSBPs1(NB,NZ)
      ENDIF
!
      call ResetDeadRootStates(NB,NZ,CPOOLK)

      IDTHY=IDTHY+1
    ENDIF
8845  CONTINUE
  end subroutine LiterfallFromDeadBranches
!------------------------------------------------------------------------------------------

  subroutine RootNoduleBiomchemistry(I,J,NZ,TFN6,WFNGR)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6(JZ1)
  real(r8), intent(in) :: WFNGR(2,JZ1)
  integer :: L,M
  real(r8) :: ZADDN
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLT
  real(r8) :: CPOOLD
  real(r8) :: CCPOLN,CZPOLN
  real(r8) :: CPPOLN
  real(r8) :: CGNDL,CPOOLNX
  real(r8) :: CCNDLR
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: GRNDG
  real(r8) :: PPOOLD
  real(r8) :: PADDN
  real(r8) :: RCO2T
  real(r8) :: RCO2TM
  real(r8) :: RCNDL,RMNDL,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,RGN2F
  real(r8) :: RXNDLC,RXNDLN,RXNDLP
  real(r8) :: RDNDLC,RDNDLN,RDNDLP
  real(r8) :: RCNDLC,RCNDLN,RCNDLP
  real(r8) :: RGNDG
  real(r8) :: RXNSNC,RXNSNN,RXNSNP
  real(r8) :: RDNSNC,RDNSNN,RDNSNP
  real(r8) :: RCNSNC,RCNSNN,RCNSNP
  real(r8) :: RCNDLM,RXNDLM,RGNDLM
  real(r8) :: RSNDLM
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTRTD1,WTNDL1,WTRTDT
  real(r8) :: XFRC,XFRN,XFRP
  real(r8) :: RCCC,RCCN,RCCP
!     begin_execution
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
  IF(INTYPs1(NZ).GE.1.AND.INTYPs1(NZ).LE.3)THEN
    DO 5400 L=NUs1,NIXs1(NZ)
      IF(WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
!
!     INITIAL INFECTION
!
        IF(WTNDLs1(L,NZ).LE.0.0)THEN
          WTNDLs1(L,NZ)=WTNDLs1(L,NZ)+WTNDI*AREA3s1(NUs1)
          WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)+WTNDI*AREA3s1(NUs1)*CNNDs1(NZ)
          WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)+WTNDI*AREA3s1(NUs1)*CPNDs1(NZ)
        ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
        IF(WTNDLs1(L,NZ).GT.ZEROPs1(NZ))THEN
          CCPOLN=AMAX1(0.0,CPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
          CZPOLN=AMAX1(0.0,ZPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
          CPPOLN=AMAX1(0.0,PPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
        ELSE
          CCPOLN=1.0_r8
          CZPOLN=1.0_r8
          CPPOLN=1.0_r8
        ENDIF
        IF(CCPOLN.GT.ZEROs1)THEN
          CCC=AMAX1(0.0,AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
            ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
!          if(curday==73)write(*,*)CCPOLN,CCPOLN,CZPOLN,CNKI
          CNC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
          CPC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        IF(WTNDLs1(L,NZ).GT.ZEROPs1(NZ))THEN
          FCNPF=AMIN1(1.0 &
            ,SQRT(WTNDLNs1(L,NZ)/(WTNDLs1(L,NZ)*CNNDs1(NZ))) &
            ,SQRT(WTNDLPs1(L,NZ)/(WTNDLs1(L,NZ)*CPNDs1(NZ))))
        ELSE
          FCNPF=1.0_r8
        ENDIF
        SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDLM=respiration from non-structural C unltd by O2
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDL=bacterial C mass
!     TFN4=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNGR=growth function of root water potential
!
        RCNDLM=AMAX1(0.0,AMIN1(CPOOLNs1(L,NZ) &
          ,VMXO*WTNDLs1(L,NZ))*FCNPF*TFN4s1(L,NZ)*WFNGR(1,L))
        CPOOLNX=CPOOLNs1(L,NZ)
!
!     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCNDL=respiration from non-structural C ltd by O2
!     WFR=constraint by O2 consumption on all root processes
!
        RCNDL=RCNDLM*WFRs1(1,L,NZ)
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN6=temperature function for root maintenance respiration
!     WTNDLN=bacterial N mass
!
        RMNDL=AMAX1(0.0,RMPLT*TFN6(L)*WTNDLNs1(L,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDLM,RXNDL=difference between non-structural C respn and mntc respn unltd,ltd by O2
!     RGNDLM,RGNDL=growth respiration unlimited by N,P and unltd,ltd by O2
!     RSNDLM,RSNDL=excess maintenance respiration unltd,ltd by O2
!
        RXNDLM=RCNDLM-RMNDL
        RXNDL=RCNDL-RMNDL
        RGNDLM=AMAX1(0.0,RXNDLM)
        RGNDL=AMAX1(0.0,RXNDL)
        RSNDLM=AMAX1(0.0,-RXNDLM)
        RSNDL=AMAX1(0.0,-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDL,WTNDLN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNF,UPNF=layer,total root N2 fixation
!
        RGN2P=AMAX1(0.0,WTNDLs1(L,NZ)*CNNDs1(NZ)-WTNDLNs1(L,NZ))/EN2F
        IF(RGNDL.GT.ZEROPs1(NZ))THEN
          RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
        ELSE
          RGN2F=0._r8
        ENDIF
        RUPNFs1(L,NZ)=RGN2F*EN2F
        UPNFs1(NZ)=UPNFs1(NZ)+RUPNFs1(L,NZ)
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTRTD=root C mass
!     CCNDLR=bacteria:root ratio
!     RDNDLX=effect of CCNDLR on bacteria decomposition rate
!     CCNKR=Km for bacterial vs root mass in decomposition
!     SPNDX=specific bacterial decomposition rate at current CCNDLR
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
        RCCC=RCCZN+CCC*RCCYN
        RCCN=CNC*RCCXN
        RCCP=CPC*RCCQN
        SPNDX=SPNDL*SQRT(TFN4s1(L,NZ)*WFNGR(1,L))
        RXNDLC=SPNDX*WTNDLs1(L,NZ)
        RXNDLN=SPNDX*WTNDLNs1(L,NZ)
        RXNDLP=SPNDX*WTNDLPs1(L,NZ)
        RDNDLC=RXNDLC*(1.0-RCCC)
        RDNDLN=RXNDLN*(1.0-RCCC)*(1.0-RCCN)
        RDNDLP=RXNDLP*(1.0-RCCC)*(1.0-RCCP)
        RCNDLC=RXNDLC-RDNDLC
        RCNDLN=RXNDLN-RDNDLN
        RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
        CGNDL=AMIN1(CPOOLNs1(L,NZ)-AMIN1(RMNDL,RCNDL) &
          -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0-DMNDs1(NZ)))
        GRNDG=CGNDL*DMNDs1(NZ)
        RGNDG=RGN2F+CGNDL*(1.0-DMNDs1(NZ))
        ZADDN=AMAX1(0.0,AMIN1(ZPOOLNs1(L,NZ),GRNDG*CNNDs1(NZ)))*CZPOLN/(CZPOLN+CZKM)
        PADDN=AMAX1(0.0,AMIN1(PPOOLNs1(L,NZ),GRNDG*CPNDs1(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
        IF(RSNDL.GT.0.0.AND.WTNDLs1(L,NZ).GT.ZEROPs1(NZ).AND.RCCC.GT.ZEROs1)THEN
          RXNSNC=RSNDL/RCCC
          RXNSNN=RXNSNC*WTNDLNs1(L,NZ)/WTNDLs1(L,NZ)
          RXNSNP=RXNSNC*WTNDLPs1(L,NZ)/WTNDLs1(L,NZ)
          RDNSNC=RXNSNC*(1.0-RCCC)
          RDNSNN=RXNSNN*(1.0-RCCC)*(1.0-RCCN)
          RDNSNP=RXNSNP*(1.0-RCCC)*(1.0-RCCP)
          RCNSNC=RXNSNC-RDNSNC
          RCNSNN=RXNSNN-RDNSNN
          RCNSNP=RXNSNP-RDNSNP
        ELSE
          RXNSNC=0._r8
          RXNSNN=0._r8
          RXNSNP=0._r8
          RDNSNC=0._r8
          RDNSNN=0._r8
          RDNSNP=0._r8
          RCNSNC=0._r8
          RCNSNN=0._r8
          RCNSNP=0._r8
        ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     RCO2A=total root respiration
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
        RCO2TM=AMIN1(RMNDL,RCNDLM)+RGNDLM+RCNSNC
        RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
        RCO2Ms1(1,L,NZ)=RCO2Ms1(1,L,NZ)+RCO2TM
        RCO2Ns1(1,L,NZ)=RCO2Ns1(1,L,NZ)+RCO2T
        RCO2As1(1,L,NZ)=RCO2As1(1,L,NZ)-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
        DO 6370 M=1,4
          CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ)*(RDNDLC+RDNSNC)
          ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ)*(RDNDLN+RDNSNN)
          PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ)*(RDNDLP+RDNSNP)
6370    CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNF=root N2 fixation
!
        CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLC
        ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)-ZADDN+RCNDLN+RCNSNN+RUPNFs1(L,NZ)
        PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
        WTNDLs1(L,NZ)=WTNDLs1(L,NZ)+GRNDG-RXNDLC-RXNSNC
        WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)+ZADDN-RXNDLN-RXNSNN
        WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)+PADDN-RXNDLP-RXNSNP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2122)'NODGR',I,J,NZ,L,RCNDLM,RCNDL,RMNDL,RGNDL,RGN2P
!    2,RGN2F,CGNDL,RSNDL,GRNDG,ZADDN,PADDN,RCCC,RCCN,RCCP
!    8,RDNDLC,RDNDLN,RDNDLP,RCNDLC,RDNDLX,WFRs1(1,L,NZ)
!    3,WTNDLs1(L,NZ),WTNDLNs1(L,NZ),WTNDLPs1(L,NZ)
!    2,CPOOLNs1(L,NZ),ZPOOLNs1(L,NZ),PPOOLNs1(L,NZ)
!    5,FCNPF,TFN4s1(L,NZ),WFNGR(1,L),PSIRTs1(1,L,NZ)
!    5,CCPOLN,CZPOLN,CPPOLN,CPOOLNX
!    6,VMXO*WTNDLs1(L,NZ)*TFN4s1(L,NZ)*FCNPF*WFNGR(1,L)
!2122  FORMAT(A8,4I4,60E14.6)
!     ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
!     WTRTD=root C mass
!     WTNDL=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGR=parameter to calculate nonstructural C,N,P exchange
!     CCNDLR=bacteria:root ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(CPOOLRs1(1,L,NZ).GT.ZEROPs1(NZ) &
          .AND.WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
          CCNDLR=WTNDLs1(L,NZ)/WTRTDs1(1,L,NZ)
          WTRTD1=WTRTDs1(1,L,NZ)
          WTNDL1=AMIN1(WTRTDs1(1,L,NZ) &
          ,AMAX1(WTNDI*AREA3s1(NUs1),WTNDLs1(L,NZ)))
          WTRTDT=WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZEROPs1(NZ))THEN
            FXRNX=FXRN(INTYPs1(NZ))/(1.0+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(INTYPs1(NZ))))
            CPOOLD=(CPOOLRs1(1,L,NZ)*WTNDL1-CPOOLNs1(L,NZ)*WTRTD1)/WTRTDT
            XFRC=FXRNX*CPOOLD
            CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)-XFRC
            CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)+XFRC
            CPOOLT=CPOOLRs1(1,L,NZ)+CPOOLNs1(L,NZ)
            IF(CPOOLT.GT.ZEROPs1(NZ))THEN
              ZPOOLD=(ZPOOLRs1(1,L,NZ)*CPOOLNs1(L,NZ) &
                -ZPOOLNs1(L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRN=FXRNX*ZPOOLD
              PPOOLD=(PPOOLRs1(1,L,NZ)*CPOOLNs1(L,NZ) &
                -PPOOLNs1(L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRP=FXRNX*PPOOLD
              ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)-XFRN
              PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)-XFRP
              ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)+XFRN
              PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)+XFRP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2122)'NODEX',I,J,NZ,L,XFRC,XFRN,XFRP
!    3,WTRTDs1(1,L,NZ),WTNDL1,CPOOLT,CCNDLR,FXRNX
!    4,WTNDLs1(L,NZ),WTNDLNs1(L,NZ),WTNDLPs1(L,NZ)
!    2,CPOOLNs1(L,NZ),ZPOOLNs1(L,NZ),PPOOLNs1(L,NZ)
!    3,CPOOLRs1(1,L,NZ),ZPOOLRs1(1,L,NZ),PPOOLRs1(1,L,NZ)
!     ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
5400  CONTINUE
  ENDIF
  end subroutine RootNoduleBiomchemistry
!------------------------------------------------------------------------------------------
  subroutine SummarizeRootSink(NZ,XRTN1,RLNT,RTSK1,RTSK2,RTNT)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: XRTN1
  real(r8),INTENT(OUT) :: RLNT(2,JZ1)
  real(r8),intent(out) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  real(r8),INTENT(OUT) :: RTNT(2)
  integer :: N,L,K,NR
  REAL(R8) :: RTDPL(10,JZ1)
  real(r8) :: CUPRL,CUPRO,CUPRC
  real(r8) :: RTDPP,RTDPS,RTSKP
  real(r8) :: RTSKS
!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER

  RLNT=0._R8
  RTSK1=0._r8
  RTSK2=0._r8
  RTNT=0._r8
  DO 4995 N=1,MYs1(NZ)
    DO 4990 L=NUs1,NIs1(NZ)
!
!     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
!     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
!
!     VOLX=soil layer volume excluding macropore, rocks
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!     RUONH4,RUONHB,RUON03,RUONOB=uptake from non-band,band of NH4,NO3 unlimited by O2
!     RUOH2P,RUOH2B,RUOH1P,RUOH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by O2
!     RUCNH4,RUCNHB,RUCN03,RUCNOB=uptake from non-band,band of NH4,NO3 unlimited by nonstructural C
!     RUCH2P,RUCH2B,RUCH1P,RUCH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by nonstructural C
!
      IF(VOLXs1(L).GT.ZEROS2s1)THEN
        CUPRL=0.86*(RUPNH4s1(N,L,NZ)+RUPNHBs1(N,L,NZ) &
          +RUPNO3s1(N,L,NZ)+RUPNOBs1(N,L,NZ)+RUPH2Ps1(N,L,NZ) &
          +RUPH2Bs1(N,L,NZ)+RUPH1Ps1(N,L,NZ)+RUPH1Bs1(N,L,NZ))
        CUPRO=0.86*(RUONH4s1(N,L,NZ)+RUONHBs1(N,L,NZ) &
          +RUONO3s1(N,L,NZ)+RUONOBs1(N,L,NZ)+RUOH2Ps1(N,L,NZ) &
          +RUOH2Bs1(N,L,NZ)+RUOH1Ps1(N,L,NZ)+RUOH1Bs1(N,L,NZ))
        CUPRC=0.86*(RUCNH4s1(N,L,NZ)+RUCNHBs1(N,L,NZ) &
          +RUCNO3s1(N,L,NZ)+RUCNOBs1(N,L,NZ)+RUCH2Ps1(N,L,NZ) &
          +RUCH2Bs1(N,L,NZ)+RUCH1Ps1(N,L,NZ)+RUCH1Bs1(N,L,NZ))
!
!     ACCUMULATE RESPIRATION IN FLUX ARRAYS
!
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
        RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)+CUPRO
        RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)+CUPRC
        RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)-CUPRL
        CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-CUPRL
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!
        DO 195 K=0,jcplx11
          CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+RDFOMCs1(N,K,L,NZ)
          ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)+RDFOMNs1(N,K,L,NZ)
          PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)+RDFOMPs1(N,K,L,NZ)
195     CONTINUE
        ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ) &
          +RUPNH4s1(N,L,NZ)+RUPNHBs1(N,L,NZ) &
          +RUPNO3s1(N,L,NZ)+RUPNOBs1(N,L,NZ)
        PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ) &
          +RUPH2Ps1(N,L,NZ)+RUPH2Bs1(N,L,NZ) &
          +RUPH1Ps1(N,L,NZ)+RUPH1Bs1(N,L,NZ)
!     IF(L.EQ.1)THEN
!     WRITE(*,9881)'CUPNH4',I,J,NZ,L,N,CPOOLRs1(N,L,NZ)
!    2,ZPOOLRs1(N,L,NZ),PPOOLRs1(N,L,NZ),CUPRL
!    2,RDFOMCs1(N,L,NZ),RDFOMNs1(N,L,NZ),RDFOMPs1(N,L,NZ)
!    2,RUPNH4s1(N,L,NZ),RUPNHBs1(N,L,NZ),RUPNO3s1(N,L,NZ)
!    2,RUPNOBs1(N,L,NZ),RUPH2Ps1(N,L,NZ),RUPH2Bs1(N,L,NZ)
!    3,RUPH12Ps1(N,L,NZ),RUPH1Bs1(N,L,NZ),WFRs1(N,L,NZ)
!9881  FORMAT(A8,5I4,30E24.16)
!     ENDIF
!
!     GROWTH OF EACH ROOT AXIS
!
        DO 4985 NR=1,NRTs1(NZ)
!
!     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
!
!     RTDP1=primary root depth from soil surface
!     RTDPP=primary root depth from canopy
!     CDPTHZ=depth from soil surface to layer bottom
!     HTSTZ=canopy height for water uptake
!     RTSK=relative primary root sink strength
!     RTSK1=primary root sink strength
!     XRTN1=number of primary root axes
!     RRAD1,RRAD2=primary, secondary root radius
!     RTNT,RLNT=total root sink strength
!
          IF(N.EQ.1)THEN
            IF(RTDP1s1(N,NR,NZ).GT.CDPTHZs1(L-1))THEN
              IF(RTDP1s1(N,NR,NZ).LE.CDPTHZs1(L))THEN
                RTDPP=RTDP1s1(N,NR,NZ)+HTSTZs1(NZ)
                RTSK1(N,L,NR)=RTSK(IGTYPs1(NZ))*XRTN1*RRAD1s1(N,L,NZ)**2/RTDPP
                RTNT(N)=RTNT(N)+RTSK1(N,L,NR)
                RLNT(N,L)=RLNT(N,L)+RTSK1(N,L,NR)
              ENDIF
            ENDIF
          ENDIF
!
!     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
!     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
!
!     RTDPL=depth of primary root axis in layer
!     RTDP1=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     RTDPX=distance behind growing point for secondary roots
!     DLYR=layer thickness
!     SDPTH=seeding depth
!     HTCTL=hypocotyledon height
!     HTSTZ=canopy height for water uptake
!     RTDPS=secondary root depth from canopy
!     RTSKP,RTSKS=primary,secondary root sink strength
!     RTN2=number of secondary root axes
!     RTSK2=total secondary root sink strength
!     RTLGA=average secondary root length
!     RTNT,RLNT=total root sink strength
!
          IF(N.EQ.1)THEN
            RTDPL(NR,L)=AMAX1(0.0,RTDP1s1(1,NR,NZ)-CDPTHZs1(L-1)-RTDPX)
            RTDPL(NR,L)=AMAX1(0.0,AMIN1(DLYR3s1(L),RTDPL(NR,L)) &
              -AMAX1(0.0,SDPTHs1(NZ)-CDPTHZs1(L-1)-HTCTLs1(NZ)))
            RTDPS=AMAX1(SDPTHs1(NZ),CDPTHZs1(L-1))+0.5*RTDPL(NR,L)+HTSTZs1(NZ)
            IF(RTDPS.GT.ZEROs1)THEN
              RTSKP=XRTN1*RRAD1s1(N,L,NZ)**2/RTDPS
              RTSKS=safe_adb(RTN2s1(N,L,NR,NZ)*RRAD2s1(N,L,NZ)**2,RTLGAs1(N,L,NZ))
              IF(RTSKP+RTSKS.GT.ZEROPs1(NZ))THEN
                RTSK2(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                RTSK2(N,L,NR)=0._r8
              ENDIF
            ELSE
              RTSK2(N,L,NR)=0._r8
            ENDIF
          ELSE
            RTSK2(N,L,NR)=safe_adb(RTN2s1(N,L,NR,NZ)*RRAD2s1(N,L,NZ)**2,RTLGAs1(N,L,NZ))
          ENDIF
          RTNT(N)=RTNT(N)+RTSK2(N,L,NR)
          RLNT(N,L)=RLNT(N,L)+RTSK2(N,L,NR)
!     IF(IYRCs1.EQ.2000.AND.I.LE.160)THEN
!     WRITE(*,3341)'SINK',I,J,NX,NY,NZ,L,NR,N,RTDP1s1(N,NR,NZ)
!    2,HTCTLs1(NZ),RTSK1(N,L,NR),RTSK2(N,L,NR),RLNT(N,L),RTNT(N)
!    3,XRTN1,PPs1(NZ),RRAD1s1(N,L,NZ),RTDPS,RTDPP
!    4,RTDPL(NR,L),RTN2s1(N,L,NR,NZ),RRAD2s1(N,L,NZ)
!    2,RTLGAs1(N,L,NZ),CDPTHZs1(L-1),CDPTHZs1(L)
!3341  FORMAT(A8,8I4,30E12.4)
!     ENDIF
4985    CONTINUE
      ENDIF
4990  CONTINUE
4995  CONTINUE
  end subroutine SummarizeRootSink
!------------------------------------------------------------------------------------------

  subroutine GrowRootAxes(N,L,L1,NZ,NRX,WFNGR,ICHK1,WFNR,WFNRG,TFN6,XRTN1,DMRTD,RLNT,&
    RTSK1,RTSK2,CNRTW,CPRTW,RTLGL,RTLGZ,WTRTX,WTRTZ)
  implicit none
  INTEGER, INTENT(IN) :: N,L,L1,NZ
  integer, intent(inout) :: NRX(2,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),XRTN1
  real(r8), intent(in) :: WFNGR(2,JZ1)
  real(r8), intent(in) :: DMRTD
  REAL(R8), INTENT(IN) :: RLNT(2,JZ1)
  real(r8), intent(in) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10),CNRTW,CPRTW
  integer, intent(inout) :: ICHK1(2,JZ1)
  real(r8), intent(inout):: WFNR,WFNRG
  real(r8), intent(out) :: RTLGL,RTLGZ,WTRTX,WTRTZ
  real(r8) :: CNRDA,CNRDM
  real(r8) :: CNPG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: FSNC1,FSNC2
  real(r8) :: ZADD1,ZADD1M,ZADD2,ZADD2M
  real(r8) :: ZPOOLB
  real(r8) :: CGRORM
  real(r8) :: CGROR
  real(r8) :: DMRTR
  real(r8) :: FNP
  real(r8) :: FRTN
  real(r8) :: FRCO2
  real(r8) :: FSNCM
  real(r8) :: FSNCP
  real(r8) :: GRTWGM
  real(r8) :: GRTLGL
  real(r8) :: GRTWTG
  real(r8) :: GRTWTP,GRTWTN,GRTWTL
  real(r8) :: GRTWTM
  real(r8) :: PPOOLB
  real(r8) :: PADD2,PADD1
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RCO2G,RCO2T
  real(r8) :: RCO2XM
  real(r8) :: RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: RMNCR,RCO2RM,RCO2R
  real(r8) :: RCCR,RCZR,RCPR
  real(r8) :: RTN2X,RTN2Y
  real(r8) :: RTDP1X,RSCS1
  REAL(R8) :: SNCR,SNCRM
  real(r8) :: TFRCO2
  real(r8) :: RCCC,RCCN,RCCP
  integer :: NR,M,LL,LX

!begin_execution

  RTLGL=0._r8
  RTLGZ=0._r8
  WTRTX=0._r8
  WTRTZ=0._r8
  DO 5050 NR=1,NRTs1(NZ)
!
!     SECONDARY ROOT EXTENSION
!
    IF(L.LE.NINRs1(NR,NZ).AND.NRX(N,NR).EQ.0)THEN
!
!     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
!     TO CURRENT AXIS
!
!     RTSK2=total secondary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of secondary root sink strength in axis
!
      IF(RLNT(N,L).GT.ZEROPs1(NZ))THEN
        FRTN=RTSK2(N,L,NR)/RLNT(N,L)
      ELSE
        FRTN=1.0_r8
      ENDIF
!
!     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
      IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
        CNPG=AMIN1(CZPOLRs1(N,L,NZ)/(CZPOLRs1(N,L,NZ) &
          +CCPOLRs1(N,L,NZ)*CNKI),CPPOLRs1(N,L,NZ) &
          /(CPPOLRs1(N,L,NZ)+CCPOLRs1(N,L,NZ)*CPKI))
      ELSE
        CNPG=1.0_r8
      ENDIF
!
!     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT2N=secondary root N mass
!     TFN6=temperature function for root maintenance respiration
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WFNGR=growth function of root water potential
!
      RMNCR=AMAX1(0.0,RMPLT*WTRT2Ns1(N,L,NR,NZ))*TFN6(L)
      IF(IGTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.2)THEN
        RMNCR=RMNCR*WFNGR(N,L)
      ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     TFN4=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     FDBKX=termination feedback inhibition on C3 CO2
!     WFNGR=growth function of root water potential
!
      RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLRs1(N,L,NZ) &
        *TFN4s1(L,NZ))*CNPG*FDBKXs1(NB1s1(NZ),NZ)*WFNGR(N,L)
!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RCO2R=RCO2RM*WFRs1(N,L,NZ)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AMAX1(0.0,RCO2XM)*WFNRG
      RCO2Y=AMAX1(0.0,RCO2X)*WFNRG
!
!     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AMAX1(0.0,ZPOOLRs1(N,L,NZ))
      PPOOLB=AMAX1(0.0,PPOOLRs1(N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTSs1(NZ) &
        ,PPOOLB*DMRTR/CPRTSs1(NZ))
      IF(RCO2YM.GT.0.0)THEN
        RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
        RCO2GM=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
        RCO2G=AMIN1(RCO2Y,FNP*WFRs1(N,L,NZ))
      ELSE
        RCO2G=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,GRTWTG=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD2M,ZADD2,PADD2=nonstructural N,P unlimited,limited by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*DMRTs1(NZ)
      GRTWTG=CGROR*DMRTs1(NZ)
      ZADD2M=AMAX1(0.0,GRTWGM*CNRTW)
      ZADD2=AMAX1(0.0,AMIN1(FRTN*ZPOOLRs1(N,L,NZ),GRTWTG*CNRTW))
      PADD2=AMAX1(0.0,AMIN1(FRTN*PPOOLRs1(N,L,NZ),GRTWTG*CPRTW))
      CNRDM=AMAX1(0.0,1.70*ZADD2M)
      CNRDA=AMAX1(0.0,1.70*ZADD2)
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAYs1(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(IDAYs1(1,NB1s1(NZ),NZ).NE.0.AND.CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
        CCC=AMAX1(0.0,AMIN1(1.0,safe_adb(CZPOLRs1(N,L,NZ),CZPOLRs1(N,L,NZ) &
          +CCPOLRs1(N,L,NZ)*CNKI) &
          ,safe_adb(CPPOLRs1(N,L,NZ),CPPOLRs1(N,L,NZ) &
          +CCPOLRs1(N,L,NZ)*CPKI)))
        CNC=AMAX1(0.0,AMIN1(1.0 &
          ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ) &
          +CZPOLRs1(N,L,NZ)/CNKI)))
        CPC=AMAX1(0.0,AMIN1(1.0 &
          ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ) &
          +CPPOLRs1(N,L,NZ)/CPKI)))
      ELSE
        CCC=0._r8
        CNC=0._r8
        CPC=0._r8
      ENDIF
      RCCC=RCCZR+CCC*RCCYR
      RCCN=CNC*RCCXR
      RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
!     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     WFR=constraint by O2 consumption on all root processes
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC2=fraction of secondary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0)THEN
        IF(-RCO2XM.LT.WTRT2s1(N,L,NR,NZ)*RCCC)THEN
          SNCRM=-RCO2XM
        ELSE
          SNCRM=AMAX1(0.0,WTRT2s1(N,L,NR,NZ)*RCCC)
        ENDIF
      ELSE
        SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0)THEN
        IF(-RCO2X.LT.WTRT2s1(N,L,NR,NZ)*RCCC)THEN
          SNCR=-RCO2X
        ELSE
          SNCR=AMAX1(0.0,WTRT2s1(N,L,NR,NZ)*RCCC)*WFRs1(N,L,NZ)
        ENDIF
      ELSE
        SNCR=0._r8
      ENDIF
      IF(SNCR.GT.0.0.AND.WTRT2s1(N,L,NR,NZ).GT.ZEROPs1(NZ))THEN
        RCCR=RCCC*WTRT2s1(N,L,NR,NZ)
        RCZR=WTRT2Ns1(N,L,NR,NZ)*(RCCN+(1.0-RCCN)*RCCR/WTRT2s1(N,L,NR,NZ))
        RCPR=WTRT2Ps1(N,L,NR,NZ)*(RCCP+(1.0-RCCP)*RCCR/WTRT2s1(N,L,NR,NZ))
        IF(RCCR.GT.ZEROPs1(NZ))THEN
          FSNC2=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
        ELSE
          FSNC2=1.0_r8
        ENDIF
      ELSE
        RCCR=0._r8
        RCZR=0._r8
        RCPR=0._r8
        FSNC2=0._r8
      ENDIF
!
!     SECONDARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC2=fraction of secondary root C to be remobilized
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      DO 6350 M=1,4
        CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
          *FSNC2*(WTRT2s1(N,L,NR,NZ)-RCCR)*FWODRs1(0)
        ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
          *FSNC2*(WTRT2Ns1(N,L,NR,NZ)-RCZR)*FWODRNs1(0)
        PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
          *FSNC2*(WTRT2Ps1(N,L,NR,NZ)-RCPR)*FWODRPs1(0)
        CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
          *FSNC2*(WTRT2s1(N,L,NR,NZ)-RCCR)*FWODRs1(1)
        ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
          *FSNC2*(WTRT2Ns1(N,L,NR,NZ)-RCZR)*FWODRNs1(1)
        PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
          *FSNC2*(WTRT2Ps1(N,L,NR,NZ)-RCPR)*FWODRPs1(1)
6350  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC2=fraction of secondary root C to be remobilized
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-AMIN1(RMNCR,RCO2R) &
        -CGROR-CNRDA-SNCR+FSNC2*RCCR
      ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-ZADD2+FSNC2*RCZR
      PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-PADD2+FSNC2*RCPR
!
!     TOTAL SECONDARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
      RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)+RCO2TM
      RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)+RCO2T
      RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)-RCO2T
!
!     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=secondary root length extension
!     GRTWTG=secondary root C growth ltd by O2
!     RTLG2X=specific secondary root length from startq.f
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     FSNC2=fraction of secondary root C to be remobilized
!     RTLG2=secondary root length
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      GRTLGL=GRTWTG*RTLG2Xs1(N,NZ)*WFNR*FWODRs1(1) &
        -FSNC2*RTLG2s1(N,L,NR,NZ)
      GRTWTL=GRTWTG-FSNC2*WTRT2s1(N,L,NR,NZ)
      GRTWTN=ZADD2-FSNC2*WTRT2Ns1(N,L,NR,NZ)
      GRTWTP=PADD2-FSNC2*WTRT2Ps1(N,L,NR,NZ)
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     RTLG2=secondary root length
!     GRTLGL=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTFQ=root branching frequency from PFT file
!     RTN2,RTNL=number of secondary root axes
!
      RTLG2s1(N,L,NR,NZ)=RTLG2s1(N,L,NR,NZ)+GRTLGL
      WTRT2s1(N,L,NR,NZ)=WTRT2s1(N,L,NR,NZ)+GRTWTL
      WTRT2Ns1(N,L,NR,NZ)=WTRT2Ns1(N,L,NR,NZ)+GRTWTN
      WTRT2Ps1(N,L,NR,NZ)=WTRT2Ps1(N,L,NR,NZ)+GRTWTP
      WSRTLs1(N,L,NZ)=WSRTLs1(N,L,NZ) &
        +AMIN1(CNWSs1(NZ)*WTRT2Ns1(N,L,NR,NZ) &
        ,CPWSs1(NZ)*WTRT2Ps1(N,L,NR,NZ))
      RTLGL=RTLGL+RTLG2s1(N,L,NR,NZ)
      WTRTX=WTRTX+WTRT2s1(N,L,NR,NZ)
      RTN2X=RTFQs1(NZ)*XRTN1
      RTN2Y=RTFQs1(NZ)*RTN2X
      RTN2s1(N,L,NR,NZ)=(RTN2X+RTN2Y)*DLYR3s1(L)
      RTNLs1(N,L,NZ)=RTNLs1(N,L,NZ)+RTN2s1(N,L,NR,NZ)
!
!     PRIMARY ROOT EXTENSION
!
!     BKDS=soil bulk density
!     RTDP1,RTDP1X=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     ICHKL=flag for identifying layer with primary root tip
!     RTN1=number of primary root axes
!     XRTN1=multiplier for number of primary root axes
!
      IF(N.EQ.1)THEN
        IF(BKDSs1(L).GT.ZEROs1)THEN
          RTDP1X=RTDP1s1(N,NR,NZ)-CDPTHZs1(0)
        ELSE
          RTDP1X=RTDP1s1(N,NR,NZ)
        ENDIF
        IF(RTDP1X.GT.CDPTHZs1(L-1).AND.ICHK1(N,NR).EQ.0)THEN
            RTN1s1(N,L,NZ)=RTN1s1(N,L,NZ)+XRTN1
            IF(RTDP1X.LE.CDPTHZs1(L).OR.L.EQ.NJs1)THEN
              ICHK1(N,NR)=1
!
!     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
!     ATTRIBUTED TO CURRENT AXIS
!
!     RTSK1=primary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of primary root sink strength in axis
!
              IF(RLNT(N,L).GT.ZEROPs1(NZ))THEN
                FRTN=RTSK1(N,L,NR)/RLNT(N,L)
              ELSE
                FRTN=1.0_r8
              ENDIF
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     RSCS,RSCS1=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
!     WFNR=water function for root extension
!     WFNRG=respiration function of root water potential
!
              RSCS1=RSCSs1(L)*RRAD1s1(N,L,NZ)/1.0E-03
              WFNR=AMIN1(1.0,AMAX1(0.0,PSIRGs1(N,L,NZ)-PSILM-RSCS1))
              IF(IGTYPs1(NZ).EQ.0)THEN
                WFNRG=WFNR**0.10
              ELSE
                WFNRG=WFNR**0.25
              ENDIF
!
!     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
              IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
                CNPG=AMIN1(CZPOLRs1(N,L,NZ)/(CZPOLRs1(N,L,NZ) &
                  +CCPOLRs1(N,L,NZ)*CNKI),CPPOLRs1(N,L,NZ) &
                  /(CPPOLRs1(N,L,NZ)+CCPOLRs1(N,L,NZ)*CPKI))
              ELSE
                CNPG=1.0_r8
              ENDIF
!
!     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT1N=primary root N mass
!     TFN6=temperature function for root maintenance respiration
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WFNGR=growth function of root water potential
!
              RMNCR=AMAX1(0.0,RMPLT*RTWT1Ns1(N,NR,NZ))*TFN6(L)
              IF(IGTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.2)THEN
                RMNCR=RMNCR*WFNGR(N,L)
              ENDIF
!
!     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of primary root sink strength in axis
!     CPOOL=non-structural C mass
!     TFN4=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     FDBKX=termination feedback inhibition on C3 CO2
!     WFNGR=growth function of root water potential
!
              RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLRs1(N,L,NZ) &
                *TFN4s1(L,NZ))*CNPG*FDBKXs1(NB1s1(NZ),NZ)*WFNGR(N,L)
              IF(RTDP1X.GE.CDPTHZs1(NJs1))THEN
                RCO2RM=AMIN1(RMNCR,RCO2RM)
              ENDIF
!
!     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
              RCO2R=RCO2RM*WFRs1(N,L,NZ)
              RCO2XM=RCO2RM-RMNCR
              RCO2X=RCO2R-RMNCR
              RCO2YM=AMAX1(0.0,RCO2XM)*WFNRG
              RCO2Y=AMAX1(0.0,RCO2X)*WFNRG
!
!     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
              DMRTR=DMRTD*FRTN
              ZPOOLB=AMAX1(0.0,ZPOOLRs1(N,L,NZ))
              PPOOLB=AMAX1(0.0,PPOOLRs1(N,L,NZ))
              FNP=AMIN1(ZPOOLB*DMRTR/CNRTSs1(NZ),PPOOLB*DMRTR/CPRTSs1(NZ))
              IF(RCO2YM.GT.0.0)THEN
                RCO2GM=AMIN1(RCO2YM,FNP)
              ELSE
                RCO2GM=0._r8
              ENDIF
              IF(RCO2Y.GT.0.0)THEN
                RCO2G=AMIN1(RCO2Y,FNP*WFRs1(N,L,NZ))
              ELSE
                RCO2G=0._r8
              ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,GRTWTG=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD1M,ZADD1,PADD1=nonstructural N,P unltd,ltd by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
              CGRORM=RCO2GM/DMRTD
              CGROR=RCO2G/DMRTD
              GRTWGM=CGRORM*DMRTs1(NZ)
              GRTWTG=CGROR*DMRTs1(NZ)
              ZADD1M=AMAX1(0.0,GRTWGM*CNRTW)
              ZADD1=AMAX1(0.0,AMIN1(FRTN*ZPOOLRs1(N,L,NZ),GRTWTG*CNRTW))
              PADD1=AMAX1(0.0,AMIN1(FRTN*PPOOLRs1(N,L,NZ),GRTWTG*CPRTW))
              CNRDM=AMAX1(0.0,1.70*ZADD1M)
              CNRDA=AMAX1(0.0,1.70*ZADD1)

              call PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC1=fraction of primary root C to be remobilized
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!
              CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-AMIN1(RMNCR,RCO2R)-CGROR-CNRDA-SNCR+FSNC1*RCCR
              ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-ZADD1+FSNC1*RCZR
              PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-PADD1+FSNC1*RCPR
!
!     TOTAL PRIMARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
              RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
              RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
!
!     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
!     THROUGH WHICH PRIMARY ROOTS GROW
!
!     RTDP1=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     RTLG1=primary root length
!     SDPTH=seeding depth
!     FRCO2=fraction of primary root respiration attributed to layer
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!
              IF(RTDP1s1(N,NR,NZ).GT.CDPTHZs1(NGs1(NZ)))THEN
                TFRCO2=0._r8
                DO 5100 LL=NGs1(NZ),NINRs1(NR,NZ)
                  IF(LL.LT.NINRs1(NR,NZ))THEN
                    FRCO2=AMIN1(1.0,RTLG1s1(N,LL,NR,NZ) &
                      /(RTDP1s1(N,NR,NZ)-SDPTHs1(NZ)))
                  ELSE
                    FRCO2=1.0_r8-TFRCO2
                  ENDIF
                  TFRCO2=TFRCO2+FRCO2
                  RCO2Ms1(N,LL,NZ)=RCO2Ms1(N,LL,NZ)+RCO2TM*FRCO2
                  RCO2Ns1(N,LL,NZ)=RCO2Ns1(N,LL,NZ)+RCO2T*FRCO2
                  RCO2As1(N,LL,NZ)=RCO2As1(N,LL,NZ)-RCO2T*FRCO2
5100            CONTINUE
              ELSE
                RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)+RCO2TM
                RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)+RCO2T
                RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)-RCO2T
              ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     GRTWTG=primary root C growth ltd by O2
!     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RTLG2=secondary root length
!
              GRTWTL=GRTWTG-FSNC1*RTWT1s1(N,NR,NZ)
              GRTWTN=ZADD1-FSNC1*RTWT1Ns1(N,NR,NZ)
              GRTWTP=PADD1-FSNC1*RTWT1Ps1(N,NR,NZ)
              IF(GRTWTL.LT.0.0)THEN
                LX=MAX(1,L-1)
                DO 5105 LL=L,LX,-1
                  GRTWTM=GRTWTL
                  IF(GRTWTL.LT.0.0)THEN
                    IF(GRTWTL.GT.-WTRT2s1(N,LL,NR,NZ))THEN
                      RTLG2s1(N,LL,NR,NZ)=RTLG2s1(N,LL,NR,NZ)+GRTWTL &
                        *RTLG2s1(N,LL,NR,NZ)/WTRT2s1(N,LL,NR,NZ)
                      WTRT2s1(N,LL,NR,NZ)=WTRT2s1(N,LL,NR,NZ)+GRTWTL
                      GRTWTL=0._r8
                    ELSE
                      GRTWTL=GRTWTL+WTRT2s1(N,LL,NR,NZ)
                      RTLG2s1(N,LL,NR,NZ)=0._r8
                      WTRT2s1(N,LL,NR,NZ)=0._r8
                    ENDIF
                  ENDIF
                  IF(GRTWTN.LT.0.0)THEN
                    IF(GRTWTN.GT.-WTRT2Ns1(N,LL,NR,NZ))THEN
                      WTRT2Ns1(N,LL,NR,NZ)=WTRT2Ns1(N,LL,NR,NZ)+GRTWTN
                      GRTWTN=0._r8
                    ELSE
                      GRTWTN=GRTWTN+WTRT2Ns1(N,LL,NR,NZ)
                      WTRT2Ns1(N,LL,NR,NZ)=0._r8
                    ENDIF
                  ENDIF
                  IF(GRTWTP.LT.0.0)THEN
                    IF(GRTWTP.GT.-WTRT2Ps1(N,LL,NR,NZ))THEN
                      WTRT2Ps1(N,LL,NR,NZ)=WTRT2Ps1(N,LL,NR,NZ)+GRTWTP
                      GRTWTP=0._r8
                    ELSE
                      GRTWTP=GRTWTP+WTRT2Ps1(N,LL,NR,NZ)
                      WTRT2Ps1(N,LL,NR,NZ)=0._r8
                    ENDIF
                  ENDIF
!
!     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
!     OF SECONDARY ROOTS
!
!     GRTWTM=negative primary root C growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     FSNCM,FSNCP=fraction of mycorrhizal structural,nonstructural C to be remobilized
!     WTRTL=active root C mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     RTLG2=mycorrhizal length
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
!
                  IF(GRTWTM.LT.0.0)THEN
                    IF(WTRT2s1(1,LL,NR,NZ).GT.ZEROPs1(NZ))THEN
                      FSNCM=AMIN1(1.0,ABS(GRTWTM)/WTRT2s1(1,LL,NR,NZ))
                    ELSE
                      FSNCM=1.0_r8
                    ENDIF
                    IF(WTRTLs1(1,LL,NZ).GT.ZEROPs1(NZ))THEN
                      FSNCP=AMIN1(1.0,ABS(GRTWTM)/WTRTLs1(1,LL,NZ))
                    ELSE
                      FSNCP=1.0_r8
                    ENDIF
                    DO 6450 M=1,4
                      CSNCs1(M,0,LL,NZ)=CSNCs1(M,0,LL,NZ)+CFOPCs1(5,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2s1(2,LL,NR,NZ))*FWODRs1(0)
                      ZSNCs1(M,0,LL,NZ)=ZSNCs1(M,0,LL,NZ)+CFOPNs1(5,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ns1(2,LL,NR,NZ))*FWODRNs1(0)
                      PSNCs1(M,0,LL,NZ)=PSNCs1(M,0,LL,NZ)+CFOPPs1(5,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ps1(2,LL,NR,NZ))*FWODRPs1(0)
                      CSNCs1(M,1,LL,NZ)=CSNCs1(M,1,LL,NZ)+CFOPCs1(4,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2s1(2,LL,NR,NZ))*FWODRs1(1)
                      ZSNCs1(M,1,LL,NZ)=ZSNCs1(M,1,LL,NZ)+CFOPNs1(4,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ns1(2,LL,NR,NZ))*FWODRNs1(1)
                      PSNCs1(M,1,LL,NZ)=PSNCs1(M,1,LL,NZ)+CFOPPs1(4,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ps1(2,LL,NR,NZ))*FWODRPs1(1)
                      CSNCs1(M,1,LL,NZ)=CSNCs1(M,1,LL,NZ)+CFOPCs1(0,M,NZ) &
                        *FSNCP*AMAX1(0.0,CPOOLRs1(2,LL,NZ))
                      ZSNCs1(M,1,LL,NZ)=ZSNCs1(M,1,LL,NZ)+CFOPNs1(0,M,NZ) &
                        *FSNCP*AMAX1(0.0,ZPOOLRs1(2,LL,NZ))
                      PSNCs1(M,1,LL,NZ)=PSNCs1(M,1,LL,NZ)+CFOPPs1(0,M,NZ) &
                        *FSNCP*AMAX1(0.0,PPOOLRs1(2,LL,NZ))
6450                CONTINUE
                    RTLG2s1(2,LL,NR,NZ)=AMAX1(0.0,RTLG2s1(2,LL,NR,NZ))*(1.0-FSNCM)
                    WTRT2s1(2,LL,NR,NZ)=AMAX1(0.0,WTRT2s1(2,LL,NR,NZ))*(1.0-FSNCM)
                    WTRT2Ns1(2,LL,NR,NZ)=AMAX1(0.0,WTRT2Ns1(2,LL,NR,NZ))*(1.0-FSNCM)
                    WTRT2Ps1(2,LL,NR,NZ)=AMAX1(0.0,WTRT2Ps1(2,LL,NR,NZ))*(1.0-FSNCM)
                    CPOOLRs1(2,LL,NZ)=AMAX1(0.0,CPOOLRs1(2,LL,NZ))*(1.0-FSNCP)
                    ZPOOLRs1(2,LL,NZ)=AMAX1(0.0,ZPOOLRs1(2,LL,NZ))*(1.0-FSNCP)
                    PPOOLRs1(2,LL,NZ)=AMAX1(0.0,PPOOLRs1(2,LL,NZ))*(1.0-FSNCP)
                  ENDIF
5105            CONTINUE
              ENDIF
!
              call PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,GRTWTG,GRTWTL,GRTWTN,GRTWTP,&
                GRTLGL,RTLGZ,WTRTZ)
            ENDIF
!
!
            IF(L.EQ.NINRs1(NR,NZ))THEN
              call WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RLNT,RTSK1,RTSK2,XRTN1)
            ENDIF

!
!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!
            IF(WTRT1s1(N,L,NR,NZ).LT.0.0)THEN
              CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+WTRT1s1(N,L,NR,NZ)
              WTRT1s1(N,L,NR,NZ)=0._r8
            ENDIF
            IF(WTRT2s1(N,L,NR,NZ).LT.0.0)THEN
              CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+WTRT2s1(N,L,NR,NZ)
              WTRT2s1(N,L,NR,NZ)=0._r8
            ENDIF
!
!     TOTAL PRIMARY ROOT LENGTH AND MASS
!
!     RTLGZ=total primary root length
!     WTRTZ=total primary root C mass
!     RTLG1=primary root length in soil layer
!     WTRT1=primary root C mass in soil layer
!     NINR=deepest root layer
!
            RTLGZ=RTLGZ+RTLG1s1(N,L,NR,NZ)
            WTRTZ=WTRTZ+WTRT1s1(N,L,NR,NZ)
            NINRs1(NR,NZ)=MIN(NINRs1(NR,NZ),NJs1)
            IF(L.EQ.NINRs1(NR,NZ))NRX(N,NR)=1
          ENDIF
        ENDIF
        RTLGZ=RTLGZ+RTLG1s1(N,L,NR,NZ)
        WTRTZ=WTRTZ+WTRT1s1(N,L,NR,NZ)
!     ENDIF
      ENDIF
      NIXs1(NZ)=MAX(NIXs1(NZ),NINRs1(NR,NZ))
5050  CONTINUE
  end subroutine GrowRootAxes

!------------------------------------------------------------------------------------------

  subroutine RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,XRTN1,WFNGR,RLNT,RTSK1,RTSK2,RTNT)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: ICHK1(2,JZ1)
  integer, intent(out) :: IDTHRN
  INTEGER, Intent(inout) :: NRX(2,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),CNRTW,CPRTW,XRTN1
  real(r8), intent(out) :: WFNGR(2,JZ1)
  REAL(R8), INTENT(OUT)  :: RLNT(2,JZ1)
  real(r8), INTENT(OUT) :: RTNT(2)
  real(r8), INTENT(OUT) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  integer :: LL,LZ,L1,L,K,lx,M,NR,N
  real(r8) :: WFNR
  real(r8) :: WFNRG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLD
  real(r8) :: CPOOLX
  real(r8) :: DMRTD
  real(r8) :: FWTRT
  real(r8) :: RSCS2
  real(r8) :: RTLGL,RTLGZ
  real(r8) :: RTLGX
  real(r8) :: RTLGT
  real(r8) :: RTVL
  real(r8) :: RTAR
  real(r8) :: WTRTTX
  real(r8) :: WTRVCX
  real(r8) :: WTRTX,WTRTZ
  real(r8) :: WTRTLX
  real(r8) :: WTRTTT
  real(r8) :: WTRTT
  real(r8) :: XFRC

!     begin_execution

  NIXs1(NZ)=NGs1(NZ)
  IDTHRN=0
!
  call SummarizeRootSink(NZ,XRTN1,RLNT,RTSK1,RTSK2,RTNT)
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
  DO 5010 N=1,MYs1(NZ)
    DO 5000 L=NUs1,NIs1(NZ)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VOLX=soil layer volume excluding macropore, rocks
!
      IF(VOLXs1(L).GT.ZEROS2s1)THEN
!     WRITE(*,4994)'5004',I,J,NZ,N,L,NIs1(NZ)
!    2,NLs1,VOLXs1(L),CDPTHZs1(L-1)
        DO 5003 LZ=L+1,NLs1
!     WRITE(*,4994)'5003',I,J,NZ,N,L,LZ
!    2,LZ,VOLXs1(L),CDPTHZs1(LZ)
          IF(VOLXs1(LZ).GT.ZEROS2s1 &
            .OR.LZ.EQ.NLs1)THEN
            L1=LZ
            EXIT
          ENDIF
5003    CONTINUE
!     WRITE(*,4994)'5005',I,J,NZ,N,L,LZ
!    2,L1,VOLXs1(L),CDPTHZs1(L1)
!4994  FORMAT(A8,7I4,12E12.4)
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     RSCS,RSCS2=soil resistance to secondary root penetration (MPa)
!     RRAD2=secondary root radius
!     WFNR=water function for root extension
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WFNGR,WFNRG=growth,respiration function of root water potential
!     PSIRT,PSIRG=root total,turgor water potential
!     DMRT=root growth yield
!
        RSCS2=RSCSs1(L)*RRAD2s1(N,L,NZ)/1.0E-03
        WFNR=AMIN1(1.0,AMAX1(0.0,PSIRGs1(N,L,NZ)-PSILM-RSCS2))
        IF(IGTYPs1(NZ).EQ.0)THEN
          WFNGR(N,L)=EXP(0.05*PSIRTs1(N,L,NZ))
          WFNRG=WFNR**0.10
        ELSE
          WFNGR(N,L)=EXP(0.10*PSIRTs1(N,L,NZ))
          WFNRG=WFNR**0.25
        ENDIF
        DMRTD=1.0_r8-DMRTs1(NZ)
!
!     FOR EACH ROOT AXIS
!
        call GrowRootAxes(N,L,L1,NZ,NRX,WFNGR,ICHK1,WFNR,WFNRG,TFN6,XRTN1,DMRTD,&
          RLNT,RTSK1,RTSK2,CNRTW,CPRTW,RTLGZ,WTRTX,WTRTZ,RTLGL)

!
!     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
!     SEASONAL STORAGE POOL IS DEPLETED
!
!     WTRTL,WTRT=total root C mass
!     WTRVC=storage C
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!     CPOOLR=non-structural C mass in root
!
        IF(L.LE.NIXs1(NZ))THEN
          IF(WTRTLs1(N,L,NZ).GT.ZEROPs1(NZ) &
            .AND.WTRTs1(NZ).GT.ZEROPs1(NZ) &
            .AND.WTRVCs1(NZ).LT.XFRX*WTRTs1(NZ))THEN
            FWTRT=WTRTLs1(N,L,NZ)/WTRTs1(NZ)
            WTRTLX=WTRTLs1(N,L,NZ)
            WTRTTX=WTRTs1(NZ)*FWTRT
            WTRTTT=WTRTLX+WTRTTX
            CPOOLX=AMAX1(0.0,CPOOLRs1(N,L,NZ))
            WTRVCX=AMAX1(0.0,WTRVCs1(NZ)*FWTRT)
            CPOOLD=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
            XFRC=AMIN1(0.0,XFRY*CPOOLD)
            CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+XFRC
            WTRVCs1(NZ)=WTRVCs1(NZ)-XFRC
!     WRITE(*,3471)'RVC',I,J,NX,NY,NZ,L
!    2,XFRC,CPOOLRs1(N,L,NZ),WTRTDs1(N,L,NZ)
!    3,WTRVCs1(NZ),WTRTs1(NZ),FWTRT
!3471  FORMAT(A8,6I4,12E12.4)
          ENDIF
        ENDIF
!
!     ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
!     TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
!
!     RTLGZ=total primary root length
!     WTRTZ=total primary root C mass
!     RTLGL=total secondary root length
!     WTRTX=total secondary root C mass
!     RTLGT=total root length
!     WTRTT=total root C mass
!     FWOOD=C woody fraction in root:0=woody,1=non-woody
!     PP=PFT population
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVL,RTVLW,RTVLP=root or myco total,aqueous,gaseous volume
!     RRAD1,RRAD2=primary,secondary root radius
!     RTARP=root surface area per plant
!     RTLGA=average secondary root length
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!
        IF(N.EQ.1)THEN
          RTLGZ=RTLGZ*FWODRs1(1)
          RTLGL=RTLGL*FWODRs1(1)
        ENDIF
        RTLGX=RTLGZ*PPs1(NZ)
        RTLGT=RTLGL+RTLGX
        WTRTT=WTRTX+WTRTZ
        IF(RTLGT.GT.ZEROPs1(NZ).AND.WTRTT.GT.ZEROPs1(NZ) &
          .AND.PPs1(NZ).GT.ZEROPs1(NZ))THEN
          RTLGPs1(N,L,NZ)=RTLGT/PPs1(NZ)
          IF(DLYR3s1(L).GT.ZEROs1)THEN
            RTDNPs1(N,L,NZ)=RTLGPs1(N,L,NZ)/DLYR3s1(L)
          ELSE
            RTDNPs1(N,L,NZ)=0._r8
          ENDIF
          RTVL=AMAX1(RTAR1Xs1(N,NZ)*RTLGX+RTAR2Xs1(N,NZ)*RTLGL &
            ,WTRTT*DMVLs1(N,NZ)*PSIRGs1(N,L,NZ))
          RTVLPs1(N,L,NZ)=PORTs1(N,NZ)*RTVL
          RTVLWs1(N,L,NZ)=(1.0-PORTs1(N,NZ))*RTVL
          RRAD1s1(N,L,NZ)=AMAX1(RRAD1Xs1(N,NZ) &
            ,(1.0+PSIRTs1(N,L,NZ)/EMODR)*RRAD1Ms1(N,NZ))
          RRAD2s1(N,L,NZ)=AMAX1(RRAD2Xs1(N,NZ) &
            ,(1.0+PSIRTs1(N,L,NZ)/EMODR)*RRAD2Ms1(N,NZ))
          RTAR=6.283*RRAD1s1(N,L,NZ)*RTLGX &
            +6.283*RRAD2s1(N,L,NZ)*RTLGL
          IF(RTNLs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
            RTLGAs1(N,L,NZ)=AMAX1(RTLGAX,RTLGL/RTNLs1(N,L,NZ))
          ELSE
            RTLGAs1(N,L,NZ)=RTLGAX
          ENDIF
          RTARPs1(N,L,NZ)=RTAR/PPs1(NZ)
!     IF(N.EQ.1)THEN
!     RTARPs1(N,L,NZ)=RTARPs1(N,L,NZ)*RTLGAX/RTLGAs1(N,L,NZ)
!     ENDIF
!     IF((I/10)*10.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2124)'RTLGA',I,J,NZ,L,N
!    2,RTLGAX,RTLGAs1(N,L,NZ),RTLGPs1(N,L,NZ),RTLGT,PPs1(NZ)
!    3,RTLGL,CPOOLRs1(N,L,NZ),WTRTDs1(N,L,NZ)
!2124  FORMAT(A8,5I4,12E12.4)
!     ENDIF
        ELSE
          RTLGPs1(N,L,NZ)=0._r8
          RTDNPs1(N,L,NZ)=0._r8
          RTVLPs1(N,L,NZ)=0._r8
          RTVLWs1(N,L,NZ)=0._r8
          RRAD1s1(N,L,NZ)=RRAD1Ms1(N,NZ)
          RRAD2s1(N,L,NZ)=RRAD2Ms1(N,NZ)
          RTARPs1(N,L,NZ)=0._r8
          RTLGAs1(N,L,NZ)=RTLGAX
          RCO2Zs1(NZ)=RCO2Zs1(NZ)-(CO2As1(N,L,NZ) &
            +CO2Ps1(N,L,NZ))
          ROXYZs1(NZ)=ROXYZs1(NZ)-(OXYAs1(N,L,NZ) &
            +OXYPs1(N,L,NZ))
          RCH4Zs1(NZ)=RCH4Zs1(NZ)-(CH4As1(N,L,NZ) &
            +CH4Ps1(N,L,NZ))
          RN2OZs1(NZ)=RN2OZs1(NZ)-(Z2OAs1(N,L,NZ) &
            +Z2OPs1(N,L,NZ))
          RNH3Zs1(NZ)=RNH3Zs1(NZ)-(ZH3As1(N,L,NZ) &
            +ZH3Ps1(N,L,NZ))
          RH2GZs1(NZ)=RH2GZs1(NZ)-(H2GAs1(N,L,NZ) &
            +H2GPs1(N,L,NZ))
          CO2As1(N,L,NZ)=0._r8
          OXYAs1(N,L,NZ)=0._r8
          CH4As1(N,L,NZ)=0._r8
          Z2OAs1(N,L,NZ)=0._r8
          ZH3As1(N,L,NZ)=0._r8
          H2GAs1(N,L,NZ)=0._r8
          CO2Ps1(N,L,NZ)=0._r8
          OXYPs1(N,L,NZ)=0._r8
          CH4Ps1(N,L,NZ)=0._r8
          Z2OPs1(N,L,NZ)=0._r8
          ZH3Ps1(N,L,NZ)=0._r8
          H2GPs1(N,L,NZ)=0._r8
        ENDIF
      ENDIF
5000  CONTINUE
5010  CONTINUE
  end subroutine RootBiochemistry

!------------------------------------------------------------------------------------------
  subroutine PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,GRTWTG,GRTWTL,GRTWTN,GRTWTP,GRTLGL,RTLGZ,WTRTZ)
  implicit none
  integer, intent(in) :: L,L1,N,NR,NZ
  real(r8), intent(in):: WFNR,FRTN,GRTWTG,GRTWTL,GRTWTN,GRTWTP
  real(r8), intent(inout) :: RTLGZ,WTRTZ
  real(r8), intent(out):: GRTLGL
  real(r8) :: FGROL,FGROZ
  REAL(R8) :: XFRC,XFRN,XFRP
! begin_execution

!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=primary root length extension
!     GRTWTG=primary root C growth ltd by O2
!     RTLG1X=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SDPTH=seeding depth
!     FSNC1=fraction of primary root C to be remobilized
!     RTLG1=primary root length
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
  IF(GRTWTL.LT.0.0.AND.RTWT1s1(N,NR,NZ) &
    .GT.ZEROPs1(NZ))THEN
    GRTLGL=GRTWTG*RTLG1Xs1(N,NZ)/PPs1(NZ)*WFNR*FWODRs1(1) &
      +GRTWTL*(RTDP1s1(N,NR,NZ)-SDPTHs1(NZ)) &
      /RTWT1s1(N,NR,NZ)
  ELSE
    GRTLGL=GRTWTG*RTLG1Xs1(N,NZ)/PPs1(NZ)*WFNR*FWODRs1(1)
  ENDIF
  IF(L.LT.NJs1)THEN
    GRTLGL=AMIN1(DLYR3s1(L1),GRTLGL)
  ENDIF
!
!     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
!     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
!     BOUNDARY OF CURRENT LAYER
!
!     GRTLGL=primary root length extension
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!
  IF(GRTLGL.GT.ZEROPs1(NZ).AND.L.LT.NJs1)THEN
    FGROL=AMAX1(0.0,AMIN1(1.0,(CDPTHZs1(L) &
      -RTDP1s1(N,NR,NZ))/GRTLGL))
    IF(FGROL.LT.1.0)FGROL=0._r8
    FGROZ=AMAX1(0.0,1.0-FGROL)
  ELSE
    FGROL=1.0_r8
    FGROZ=0._r8
  ENDIF
!
!     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
!     AND AXIS NUMBER
!
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     GRTLGL=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTLG1=primary root length
!
  RTWT1s1(N,NR,NZ)=RTWT1s1(N,NR,NZ)+GRTWTL
  RTWT1Ns1(N,NR,NZ)=RTWT1Ns1(N,NR,NZ)+GRTWTN
  RTWT1Ps1(N,NR,NZ)=RTWT1Ps1(N,NR,NZ)+GRTWTP
  RTDP1s1(N,NR,NZ)=RTDP1s1(N,NR,NZ)+GRTLGL
  WTRT1s1(N,L,NR,NZ)=WTRT1s1(N,L,NR,NZ)+GRTWTL*FGROL
  WTRT1Ns1(N,L,NR,NZ)=WTRT1Ns1(N,L,NR,NZ)+GRTWTN*FGROL
  WTRT1Ps1(N,L,NR,NZ)=WTRT1Ps1(N,L,NR,NZ)+GRTWTP*FGROL
  WSRTLs1(N,L,NZ)=WSRTLs1(N,L,NZ) &
    +AMIN1(CNWSs1(NZ)*WTRT1Ns1(N,L,NR,NZ) &
    ,CPWSs1(NZ)*WTRT1Ps1(N,L,NR,NZ))
  RTLG1s1(N,L,NR,NZ)=RTLG1s1(N,L,NR,NZ)+GRTLGL*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of GRTLGL in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     WTRTD=root C mass
!     RTLG1=primary root length
!     GRTLGL=primary root length extension
!     FRTN=fraction of primary root sink strength in axis
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     PSIRT,PSIRG,PSIRO=root total,turgor,osmotic water potential
!     NINR=deepest root layer
!
  IF(FGROZ.GT.0.0)THEN
    WTRT1s1(N,L1,NR,NZ)=WTRT1s1(N,L1,NR,NZ)+GRTWTL*FGROZ
    WTRT1Ns1(N,L1,NR,NZ)=WTRT1Ns1(N,L1,NR,NZ)+GRTWTN*FGROZ
    WTRT1Ps1(N,L1,NR,NZ)=WTRT1Ps1(N,L1,NR,NZ)+GRTWTP*FGROZ
    WSRTLs1(N,L1,NZ)=WSRTLs1(N,L1,NZ)+AMIN1(CNWSs1(NZ)*WTRT1Ns1(N,L1,NR,NZ) &
      ,CPWSs1(NZ)*WTRT1Ps1(N,L1,NR,NZ))
    WTRTDs1(N,L1,NZ)=WTRTDs1(N,L1,NZ)+WTRT1s1(N,L1,NR,NZ)
    RTLG1s1(N,L1,NR,NZ)=RTLG1s1(N,L1,NR,NZ)+GRTLGL*FGROZ
    RRAD1s1(N,L1,NZ)=RRAD1s1(N,L,NZ)
    RTLGZ=RTLGZ+RTLG1s1(N,L1,NR,NZ)
    WTRTZ=WTRTZ+WTRT1s1(N,L1,NR,NZ)
    XFRC=FRTN*CPOOLRs1(N,L,NZ)
    XFRN=FRTN*ZPOOLRs1(N,L,NZ)
    XFRP=FRTN*PPOOLRs1(N,L,NZ)
    CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-XFRC
    ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-XFRN
    PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-XFRP
    CPOOLRs1(N,L1,NZ)=CPOOLRs1(N,L1,NZ)+XFRC
    ZPOOLRs1(N,L1,NZ)=ZPOOLRs1(N,L1,NZ)+XFRN
    PPOOLRs1(N,L1,NZ)=PPOOLRs1(N,L1,NZ)+XFRP
    PSIRTs1(N,L1,NZ)=PSIRTs1(N,L,NZ)
    PSIROs1(N,L1,NZ)=PSIROs1(N,L,NZ)
    PSIRGs1(N,L1,NZ)=PSIRGs1(N,L,NZ)
    NINRs1(NR,NZ)=MAX(NGs1(NZ),L+1)
!     IF(NZ.EQ.2)THEN
!     WRITE(*,9877)'INFIL',I,J,NZ,NR,L,N,NINRs1(NR,NZ)
!    2,FRTN,WTRTDs1(N,L1,NZ),CPOOLRs1(N,L1,NZ)
!    2,FGROZ,RTDP1s1(N,NR,NZ),GRTLGL,CDPTHZs1(L)
!     ENDIF
  ENDIF

  end subroutine PrimRootExtension
!------------------------------------------------------------------------------------------

  subroutine PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

  implicit none
  integer, intent(in) :: N,L,NZ,NR
  real(r8), intent(in) :: RCO2X,RCO2XM
  real(r8), intent(out) :: FSNC1
  integer :: M
  real(r8) :: CCC,CNC,CPC
  real(r8) :: RCCR,RCZR,RCPR
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SNCR,SNCRM
! begin_execution
!
!     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAYs1(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
  IF(IDAYs1(1,NB1s1(NZ),NZ).NE.0 &
    .AND.CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
    CCC=AMAX1(0.0,AMIN1(1.0 &
      ,safe_adb(CZPOLRs1(N,L,NZ),CZPOLRs1(N,L,NZ) &
      +CCPOLRs1(N,L,NZ)*CNKI) &
      ,safe_adb(CPPOLRs1(N,L,NZ),CPPOLRs1(N,L,NZ) &
      +CCPOLRs1(N,L,NZ)*CPKI)))
    CNC=AMAX1(0.0,AMIN1(1.0 &
      ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ) &
      +CZPOLRs1(N,L,NZ)/CNKI)))
    CPC=AMAX1(0.0,AMIN1(1.0 &
      ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ) &
      +CPPOLRs1(N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8
    CNC=0._r8
    CPC=0._r8
  ENDIF
  RCCC=RCCZR+CCC*RCCYR
  RCCN=CNC*RCCXR
  RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
!     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     WFR=constraint by O2 consumption on all root processes
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC1=fraction of primary root C to be remobilized
!
  IF(-RCO2XM.GT.0.0)THEN
    IF(-RCO2XM.LT.RTWT1s1(N,NR,NZ)*RCCC)THEN
      SNCRM=-RCO2XM
    ELSE
      SNCRM=AMAX1(0.0,RTWT1s1(N,NR,NZ)*RCCC)
    ENDIF
  ELSE
    SNCRM=0._r8
  ENDIF
  IF(-RCO2X.GT.0.0)THEN
    IF(-RCO2X.LT.RTWT1s1(N,NR,NZ)*RCCC)THEN
      SNCR=-RCO2X
    ELSE
      SNCR=AMAX1(0.0,RTWT1s1(N,NR,NZ)*RCCC)*WFRs1(N,L,NZ)
    ENDIF
  ELSE
    SNCR=0._r8
  ENDIF
  IF(SNCR.GT.0.0.AND.RTWT1s1(N,NR,NZ).GT.ZEROPs1(NZ))THEN
    RCCR=RCCC*RTWT1s1(N,NR,NZ)
    RCZR=RTWT1Ns1(N,NR,NZ)*(RCCN+(1.0-RCCN)*RCCR/RTWT1s1(N,NR,NZ))
    RCPR=RTWT1Ps1(N,NR,NZ)*(RCCP+(1.0-RCCP)*RCCR/RTWT1s1(N,NR,NZ))
    IF(RCCR.GT.ZEROPs1(NZ))THEN
      FSNC1=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
    ELSE
      FSNC1=1.0_r8
    ENDIF
  ELSE
    RCCR=0._r8
    RCZR=0._r8
    RCPR=0._r8
    FSNC1=0._r8
  ENDIF
!
!     PRIMARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
  DO 6355 M=1,4
    CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
      *FSNC1*(RTWT1s1(N,NR,NZ)-RCCR)*FWODRs1(0)
    ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
      *FSNC1*(RTWT1Ns1(N,NR,NZ)-RCZR)*FWODRNs1(0)
    PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
      *FSNC1*(RTWT1Ps1(N,NR,NZ)-RCPR)*FWODRPs1(0)
    CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
      *FSNC1*(RTWT1s1(N,NR,NZ)-RCCR)*FWODRs1(1)
    ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
      *FSNC1*(RTWT1Ns1(N,NR,NZ)-RCZR)*FWODRNs1(1)
    PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
      *FSNC1*(RTWT1Ps1(N,NR,NZ)-RCPR)*FWODRPs1(1)
6355  CONTINUE
end subroutine PrimRootRemobilization
!------------------------------------------------------------------------------------------

  subroutine WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RLNT,RTSK1,RTSK2,XRTN1)
  implicit none
  integer, intent(in) :: L,NR,NZ,N
  real(r8), intent(in):: RTDP1X
  real(r8), INTENT(IN) :: RLNT(2,JZ1)
  real(r8),intent(in) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10),XRTN1
  integer :: LL,NN
  real(r8) :: XFRD,XFRW,FRTN
  real(r8) :: XFRC,XFRN,XFRP


! begin_execution
!     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
!     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
!     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
!     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
!
!     NINR=deepest root layer
!     VOLX=soil layer volume excluding macropore, rocks
!     RTDP1X=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!     FRTN=fraction of primary root sink strength in axis
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTLG1=primary root length
!     WSRTL=root protein C mass
!     WTRTD=root C mass
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  DO 5115 LL=L,NGs1(NZ)+1,-1
    IF(VOLXs1(LL-1).GT.ZEROS2s1 &
      .AND.(RTDP1X.LT.CDPTHZs1(LL-1) &
      .OR.RTDP1X.LT.SDPTHs1(NZ)))THEN
      IF(RLNT(N,LL).GT.ZEROPs1(NZ))THEN
        FRTN=(RTSK1(N,LL,NR)+RTSK2(N,LL,NR))/RLNT(N,LL)
      ELSE
        FRTN=1.0_r8
      ENDIF
      DO 5110 NN=1,MYs1(NZ)
        WTRT1s1(NN,LL-1,NR,NZ)=WTRT1s1(NN,LL-1,NR,NZ)+WTRT1s1(NN,LL,NR,NZ)
        WTRT1Ns1(NN,LL-1,NR,NZ)=WTRT1Ns1(NN,LL-1,NR,NZ)+WTRT1Ns1(NN,LL,NR,NZ)
        WTRT1Ps1(NN,LL-1,NR,NZ)=WTRT1Ps1(NN,LL-1,NR,NZ)+WTRT1Ps1(NN,LL,NR,NZ)
        WTRT2s1(NN,LL-1,NR,NZ)=WTRT2s1(NN,LL-1,NR,NZ)+WTRT2s1(NN,LL,NR,NZ)
        WTRT2Ns1(NN,LL-1,NR,NZ)=WTRT2Ns1(NN,LL-1,NR,NZ)+WTRT2Ns1(NN,LL,NR,NZ)
        WTRT2Ps1(NN,LL-1,NR,NZ)=WTRT2Ps1(NN,LL-1,NR,NZ)+WTRT2Ps1(NN,LL,NR,NZ)
        RTLG1s1(NN,LL-1,NR,NZ)=RTLG1s1(NN,LL-1,NR,NZ)+RTLG1s1(NN,LL,NR,NZ)
        WTRT1s1(NN,LL,NR,NZ)=0._r8
        WTRT1Ns1(NN,LL,NR,NZ)=0._r8
        WTRT1Ps1(NN,LL,NR,NZ)=0._r8
        WTRT2s1(NN,LL,NR,NZ)=0._r8
        WTRT2Ns1(NN,LL,NR,NZ)=0._r8
        WTRT2Ps1(NN,LL,NR,NZ)=0._r8
        RTLG1s1(NN,LL,NR,NZ)=0._r8
        XFRC=FRTN*CPOOLRs1(NN,LL,NZ)
        XFRN=FRTN*ZPOOLRs1(NN,LL,NZ)
        XFRP=FRTN*PPOOLRs1(NN,LL,NZ)
        XFRW=FRTN*WSRTLs1(NN,L,NZ)
        XFRD=FRTN*WTRTDs1(NN,LL,NZ)
        CPOOLRs1(NN,LL,NZ)=CPOOLRs1(NN,LL,NZ)-XFRC
        ZPOOLRs1(NN,LL,NZ)=ZPOOLRs1(NN,LL,NZ)-XFRN
        PPOOLRs1(NN,LL,NZ)=PPOOLRs1(NN,LL,NZ)-XFRP
        WSRTLs1(NN,LL,NZ)=WSRTLs1(NN,LL,NZ)-XFRW
        WTRTDs1(NN,LL,NZ)=WTRTDs1(NN,LL,NZ)-XFRD
        CPOOLRs1(NN,LL-1,NZ)=CPOOLRs1(NN,LL-1,NZ)+XFRC
        ZPOOLRs1(NN,LL-1,NZ)=ZPOOLRs1(NN,LL-1,NZ)+XFRN
        PPOOLRs1(NN,LL-1,NZ)=PPOOLRs1(NN,LL-1,NZ)+XFRP
        WSRTLs1(NN,LL-1,NZ)=WSRTLs1(NN,LL-1,NZ)+XFRW
        WTRTDs1(NN,LL-1,NZ)=WTRTDs1(NN,LL-1,NZ)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
        RCO2Zs1(NZ)=RCO2Zs1(NZ)-FRTN*(CO2As1(NN,LL,NZ)+CO2Ps1(NN,LL,NZ))
        ROXYZs1(NZ)=ROXYZs1(NZ)-FRTN*(OXYAs1(NN,LL,NZ)+OXYPs1(NN,LL,NZ))
        RCH4Zs1(NZ)=RCH4Zs1(NZ)-FRTN*(CH4As1(NN,LL,NZ)+CH4Ps1(NN,LL,NZ))
        RN2OZs1(NZ)=RN2OZs1(NZ)-FRTN*(Z2OAs1(NN,LL,NZ)+Z2OPs1(NN,LL,NZ))
        RNH3Zs1(NZ)=RNH3Zs1(NZ)-FRTN*(ZH3As1(NN,LL,NZ)+ZH3Ps1(NN,LL,NZ))
        RH2GZs1(NZ)=RH2GZs1(NZ)-FRTN*(H2GAs1(NN,LL,NZ)+H2GPs1(NN,LL,NZ))
        CO2As1(NN,LL,NZ)=(1.0-FRTN)*CO2As1(NN,LL,NZ)
        OXYAs1(NN,LL,NZ)=(1.0-FRTN)*OXYAs1(NN,LL,NZ)
        CH4As1(NN,LL,NZ)=(1.0-FRTN)*CH4As1(NN,LL,NZ)
        Z2OAs1(NN,LL,NZ)=(1.0-FRTN)*Z2OAs1(NN,LL,NZ)
        ZH3As1(NN,LL,NZ)=(1.0-FRTN)*ZH3As1(NN,LL,NZ)
        H2GAs1(NN,LL,NZ)=(1.0-FRTN)*H2GAs1(NN,LL,NZ)
        CO2Ps1(NN,LL,NZ)=(1.0-FRTN)*CO2Ps1(NN,LL,NZ)
        OXYPs1(NN,LL,NZ)=(1.0-FRTN)*OXYPs1(NN,LL,NZ)
        CH4Ps1(NN,LL,NZ)=(1.0-FRTN)*CH4Ps1(NN,LL,NZ)
        Z2OPs1(NN,LL,NZ)=(1.0-FRTN)*Z2OPs1(NN,LL,NZ)
        ZH3Ps1(NN,LL,NZ)=(1.0-FRTN)*ZH3Ps1(NN,LL,NZ)
        H2GPs1(NN,LL,NZ)=(1.0-FRTN)*H2GPs1(NN,LL,NZ)
!     IF(NZ.EQ.2)THEN
!     WRITE(*,9868)'WITHDR',I,J,NZ,NR,LL,NN,NINRs1(NR,NZ)
!    2,FRTN,RTSK1(N,LL,NR),RTSK2(N,LL,NR),RLNT(N,LL)
!    2,WTRTDs1(NN,LL-1,NZ),WTRTDs1(NN,LL,NZ)
!    2,RTLG1s1(NN,LL-1,NR,NZ),RTLG1s1(NN,LL,NR,NZ)
!    2,RTLG2s1(NN,LL-1,NR,NZ),RTLG2s1(NN,LL,NR,NZ)
!    3,RTDP1s1(N,NR,NZ),RTDP1s1(NN,NR,NZ)
!    4,CPOOLRs1(NN,LL-1,NZ),CPOOLRs1(NN,LL,NZ)
!    4,WTRT1s1(NN,LL-1,NR,NZ),WTRT1s1(NN,LL,NR,NZ)
!    4,WTRT2s1(NN,LL-1,NR,NZ),WTRT2s1(NN,LL,NR,NZ)
!9868  FORMAT(A8,7I4,100E24.16)
!      ENDIF
5110  CONTINUE
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     RTN2,RTNL=number of secondary root axes
!     RTN1=number of primary root axes
!     RTLG1=primary root length
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!
      RTNLs1(N,LL,NZ)=RTNLs1(N,LL,NZ)-RTN2s1(N,LL,NR,NZ)
      RTNLs1(N,LL-1,NZ)=RTNLs1(N,LL-1,NZ)+RTN2s1(N,LL,NR,NZ)
      RTN2s1(N,LL,NR,NZ)=0._r8
      RTN1s1(N,LL,NZ)=RTN1s1(N,LL,NZ)-XRTN1
      IF(LL-1.GT.NGs1(NZ))THEN
        RTLG1s1(N,LL-1,NR,NZ)=DLYR3s1(LL-1)-(CDPTHZs1(LL-1)-RTDP1s1(N,NR,NZ))
      ELSE
        RTLG1s1(N,LL-1,NR,NZ)=DLYR3s1(LL-1)-(CDPTHZs1(LL-1)-RTDP1s1(N,NR,NZ)) &
          -(SDPTHs1(NZ)-CDPTHZs1(LL-2))
      ENDIF
!
!     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     FRTN=fraction of primary root sink strength in axis
!     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
!
      IF(INTYPs1(NZ).GE.1.AND.INTYPs1(NZ).LE.3)THEN
        XFRC=FRTN*WTNDLs1(LL,NZ)
        XFRN=FRTN*WTNDLNs1(LL,NZ)
        XFRP=FRTN*WTNDLPs1(LL,NZ)
        WTNDLs1(LL,NZ)=WTNDLs1(LL,NZ)-XFRC
        WTNDLNs1(LL,NZ)=WTNDLNs1(LL,NZ)-XFRN
        WTNDLPs1(LL,NZ)=WTNDLPs1(LL,NZ)-XFRP
        WTNDLs1(LL-1,NZ)=WTNDLs1(LL-1,NZ)+XFRC
        WTNDLNs1(LL-1,NZ)=WTNDLNs1(LL-1,NZ)+XFRN
        WTNDLPs1(LL-1,NZ)=WTNDLPs1(LL-1,NZ)+XFRP
        XFRC=FRTN*CPOOLNs1(LL,NZ)
        XFRN=FRTN*ZPOOLNs1(LL,NZ)
        XFRP=FRTN*PPOOLNs1(LL,NZ)
        CPOOLNs1(LL,NZ)=CPOOLNs1(LL,NZ)-XFRC
        ZPOOLNs1(LL,NZ)=ZPOOLNs1(LL,NZ)-XFRN
        PPOOLNs1(LL,NZ)=PPOOLNs1(LL,NZ)-XFRP
        CPOOLNs1(LL-1,NZ)=CPOOLNs1(LL-1,NZ)+XFRC
        ZPOOLNs1(LL-1,NZ)=ZPOOLNs1(LL-1,NZ)+XFRN
        PPOOLNs1(LL-1,NZ)=PPOOLNs1(LL-1,NZ)+XFRP
!     WRITE(*,9868)'WITHDRN',I,J,NZ,NR,LL,NN,NINRs1(NR,NZ)
!    2,WTNDLs1(LL,NZ),CPOOLNs1(LL,NZ),RTDP1s1(N,NR,NZ)
      ENDIF
      NINRs1(NR,NZ)=MAX(NGs1(NZ),LL-1)
    ELSE
      EXIT
    ENDIF
5115  CONTINUE
  end subroutine WithdrawPrimRoot
!------------------------------------------------------------------------------------------

  subroutine CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WFNG,UPNFC)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: TFN5,WFNG
  real(r8), intent(inout) :: UPNFC(JP1)
  integer :: M
  real(r8) :: ZADDN,ZPOOLD
  real(r8) :: CCC,CNC,CPC
  REAL(R8) :: cpoolt
  real(r8) :: CPOOLD
  real(r8) :: CCPOLN,CZPOLN
  real(r8) :: CPPOLN,CGNDL
  real(r8) :: CCNDLB
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: GRNDG
  real(r8) :: PPOOLD
  real(r8) :: PADDN,RCO2T
  real(r8) :: RCNDL,RMNDL,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,RGN2F
  real(r8) :: RUPNFB
  real(r8) :: RXNDLC,RXNDLN,RXNDLP
  real(r8) :: RDNDLC,RDNDLN,RDNDLP
  real(r8) :: RCNDLC,RCNDLN,RCNDLP
  real(r8) :: RGNDG
  real(r8) :: RXNSNC,RXNSNN,RXNSNP
  real(r8) :: RDNSNC,RDNSNN,RDNSNP
  real(r8) :: RCNSNC,RCNSNN,RCNSNP
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTLSB1,WTNDB1,WTLSBT
  real(r8) :: XFRC,XFRN,XFRP
  REAL(R8) :: RCCC,RCCN,RCCP
!     begin_execution
!     INTYP=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
  IF(INTYPs1(NZ).GE.4)THEN
!
!     INITIAL INFECTION
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
    IF(WTNDBs1(NB,NZ).LE.0.0)THEN
      WTNDBs1(NB,NZ)=WTNDBs1(NB,NZ) &
        +WTNDI*AREA3s1(NUs1)
      WTNDBNs1(NB,NZ)=WTNDBNs1(NB,NZ) &
        +WTNDI*AREA3s1(NUs1)*CNNDs1(NZ)
      WTNDBPs1(NB,NZ)=WTNDBPs1(NB,NZ) &
        +WTNDI*AREA3s1(NUs1)*CPNDs1(NZ)
    ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
    IF(WTNDBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      CCPOLN=AMAX1(0.0,CPOLNBs1(NB,NZ)/WTNDBs1(NB,NZ))
      CZPOLN=AMAX1(0.0,ZPOLNBs1(NB,NZ)/WTNDBs1(NB,NZ))
      CPPOLN=AMAX1(0.0,PPOLNBs1(NB,NZ)/WTNDBs1(NB,NZ))
    ELSE
      CCPOLN=1.0_r8
      CZPOLN=1.0_r8
      CPPOLN=1.0_r8
    ENDIF
    IF(CCPOLN.GT.ZEROs1)THEN
      CCC=AMAX1(0.0,AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
        ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    IF(WTNDBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      FCNPF=AMIN1(1.0 &
        ,SQRT(WTNDBNs1(NB,NZ)/(WTNDBs1(NB,NZ)*CNNDs1(NZ))) &
        ,SQRT(WTNDBPs1(NB,NZ)/(WTNDBs1(NB,NZ)*CPNDs1(NZ))))
    ELSE
      FCNPF=1.0_r8
    ENDIF
    SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDL=respiration from non-structural C
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDB=bacterial C mass
!     TFN3=temperature function for canopy growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNG=growth function of canopy water potential
!
    RCNDL=AMAX1(0.0,AMIN1(CPOLNBs1(NB,NZ) &
      ,VMXO*WTNDBs1(NB,NZ))*FCNPF*TFN3s1(NZ)*WFNG)
!     CPOOLNX=CPOLNBs1(NB,NZ)
!     VMXOX=VMXO*WTNDBs1(NB,NZ)*FCNPF*TFN3s1(NZ)*WFNG
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN5=temperature function for canopy maintenance respiration
!     WTNDBN=bacterial N mass
!
    RMNDL=AMAX1(0.0,RMPLT*TFN5*WTNDBNs1(NB,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDL=difference between non-structural C respn and mntc respn
!     RGNDL=growth respiration unlimited by N,P
!     RSNDL=excess maintenance respiration
!
    RXNDL=RCNDL-RMNDL
    RGNDL=AMAX1(0.0,RXNDL)
    RSNDL=AMAX1(0.0,-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDB,WTNDBN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNFB,UPNFC=branch,total N2 fixation
!
    RGN2P=AMAX1(0.0,WTNDBs1(NB,NZ)*CNNDs1(NZ)-WTNDBNs1(NB,NZ))/EN2F
    IF(RGNDL.GT.ZEROPs1(NZ))THEN
      RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
    ELSE
      RGN2F=0._r8
    ENDIF
    RUPNFB=RGN2F*EN2F
    UPNFC(NZ)=UPNFC(NZ)+RUPNFB
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION AND LEAKAGE
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTLSB=leaf+petiole mass
!     CCNDLB=bacteria:leaf+petiole ratio
!     RDNDBX=effect of CCNDLB on bacteria decomposition rate
!     SPNDX=specific bacterial decomposition rate at current CCNDLB
!     SPNDL=specific decomposition rate by bacterial N2 fixers
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
    RCCC=RCCZN+CCC*RCCYN
    RCCN=CNC*RCCXN
    RCCP=CPC*RCCQN
    SPNDX=SPNDL*SQRT(TFN3s1(NZ)*WFNG)
    RXNDLC=SPNDX*WTNDBs1(NB,NZ)
    RXNDLN=SPNDX*WTNDBNs1(NB,NZ)
    RXNDLP=SPNDX*WTNDBPs1(NB,NZ)
    RDNDLC=RXNDLC*(1.0-RCCC)
    RDNDLN=RXNDLN*(1.0-RCCC)*(1.0-RCCN)
    RDNDLP=RXNDLP*(1.0-RCCC)*(1.0-RCCP)
    RCNDLC=RXNDLC-RDNDLC
    RCNDLN=RXNDLN-RDNDLN
    RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
    CGNDL=AMIN1(CPOLNBs1(NB,NZ)-AMIN1(RMNDL,RCNDL) &
      -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0-DMNDs1(NZ)))
    GRNDG=CGNDL*DMNDs1(NZ)
    RGNDG=RGN2F+CGNDL*(1.0-DMNDs1(NZ))
    ZADDN=AMAX1(0.0,AMIN1(ZPOLNBs1(NB,NZ) &
      ,GRNDG*CNNDs1(NZ)))*CZPOLN/(CZPOLN+CZKM)
    PADDN=AMAX1(0.0,AMIN1(PPOLNBs1(NB,NZ) &
      ,GRNDG*CPNDs1(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
    IF(RSNDL.GT.0.0.AND.WTNDBs1(NB,NZ).GT.ZEROPs1(NZ).AND.RCCC.GT.ZEROs1)THEN
      RXNSNC=RSNDL/RCCC
      RXNSNN=RXNSNC*WTNDBNs1(NB,NZ)/WTNDBs1(NB,NZ)
      RXNSNP=RXNSNC*WTNDBPs1(NB,NZ)/WTNDBs1(NB,NZ)
      RDNSNC=RXNSNC*(1.0-RCCC)
      RDNSNN=RXNSNN*(1.0-RCCC)*(1.0-RCCN)
      RDNSNP=RXNSNP*(1.0-RCCC)*(1.0-RCCP)
      RCNSNC=RXNSNC-RDNSNC
      RCNSNN=RXNSNN-RDNSNN
      RCNSNP=RXNSNP-RDNSNP
    ELSE
      RXNSNC=0._r8
      RXNSNN=0._r8
      RXNSNP=0._r8
      RDNSNC=0._r8
      RDNSNN=0._r8
      RDNSNP=0._r8
      RCNSNC=0._r8
      RCNSNN=0._r8
      RCNSNP=0._r8
    ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2T=total C respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
    RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
    TCO2Ts1(NZ)=TCO2Ts1(NZ)-RCO2T
    TCO2As1(NZ)=TCO2As1(NZ)-RCO2T
    CNETs1(NZ)=CNETs1(NZ)-RCO2T
    RECOs1=RECOs1-RCO2T
    TRAUs1=TRAUs1-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
    DO 6470 M=1,4
      CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(1,M,NZ)*(RDNDLC+RDNSNC)
      ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(1,M,NZ)*(RDNDLN+RDNSNN)
      PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(1,M,NZ)*(RDNDLP+RDNSNP)
6470  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNFB=branch N2 fixation
!
    CPOLNBs1(NB,NZ)=CPOLNBs1(NB,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLC
    ZPOLNBs1(NB,NZ)=ZPOLNBs1(NB,NZ)-ZADDN+RCNDLN+RCNSNN+RUPNFB
    PPOLNBs1(NB,NZ)=PPOLNBs1(NB,NZ)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
    WTNDBs1(NB,NZ)=WTNDBs1(NB,NZ)+GRNDG-RXNDLC-RXNSNC
    WTNDBNs1(NB,NZ)=WTNDBNs1(NB,NZ)+ZADDN-RXNDLN-RXNSNN
    WTNDBPs1(NB,NZ)=WTNDBPs1(NB,NZ)+PADDN-RXNDLP-RXNSNP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2121)'NODGR',I,J,NZ,NB
!    2,RCNDL,RMNDL,RGNDL,RGN2P,RCO2T,RXNDLC,SPNDLI
!    2,RGN2P,RGN2F,CGNDL,RSNDL,GRNDG,RGNDG,RCNSNC
!    3,ZADDN,PADDN,RCCC,RCCN
!    8,RCCP,RDNDLC,RDNDLN,RDNDLP,RDNDLX,WTLSBs1(NB,NZ)
!    3,WTNDBs1(NB,NZ),WTNDBNs1(NB,NZ),WTNDBPs1(NB,NZ)
!    4,CPOLNBs1(NB,NZ),ZPOLNBs1(NB,NZ),PPOLNBs1(NB,NZ)
!    5,CCPOLN,CZPOLN,CPPOLN
!    6,TFN3s1(NZ),FCNPF,WFNG,CCNDLB,RDNDBX,CPOOLNX,VMXOX
!2121  FORMAT(A8,4I4,60F16.8)
!     ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     WTLSB=leaf+petiole C mass
!     WTNDB=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGB=parameter to calculate nonstructural C,N,P exchange
!     CCNDLB=bacteria:leaf+petiole ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!
    IF(CPOOLs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.WTLSBs1(NB,NZ).GT.ZEROLs1(NZ))THEN
      CCNDLB=WTNDBs1(NB,NZ)/WTLSBs1(NB,NZ)
      WTLSB1=WTLSBs1(NB,NZ)
      WTNDB1=AMIN1(WTLSBs1(NB,NZ),AMAX1(WTNDI*AREA3s1(NUs1),WTNDBs1(NB,NZ)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROPs1(NZ))THEN
        FXRNX=FXRN(INTYPs1(NZ))/(1.0+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(INTYPs1(NZ))))
        CPOOLD=(CPOOLs1(NB,NZ)*WTNDB1-CPOLNBs1(NB,NZ)*WTLSB1)/WTLSBT
        XFRC=FXRNX*CPOOLD
        CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)-XFRC
        CPOLNBs1(NB,NZ)=CPOLNBs1(NB,NZ)+XFRC
        CPOOLT=CPOOLs1(NB,NZ)+CPOLNBs1(NB,NZ)
        IF(CPOOLT.GT.ZEROPs1(NZ))THEN
          ZPOOLD=(ZPOOLs1(NB,NZ)*CPOLNBs1(NB,NZ)-ZPOLNBs1(NB,NZ)*CPOOLs1(NB,NZ))/CPOOLT
          XFRN=FXRNX*ZPOOLD
          PPOOLD=(PPOOLs1(NB,NZ)*CPOLNBs1(NB,NZ) &
            -PPOLNBs1(NB,NZ)*CPOOLs1(NB,NZ))/CPOOLT
          XFRP=FXRNX*PPOOLD
          ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-XFRN
          PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-XFRP
          ZPOLNBs1(NB,NZ)=ZPOLNBs1(NB,NZ)+XFRN
          PPOLNBs1(NB,NZ)=PPOLNBs1(NB,NZ)+XFRP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2120)'NODEX',I,J,NZ,NB,IFLGAs1(NB,NZ)
!    2,XFRC,XFRN,XFRP
!    3,WTLSBs1(NB,NZ),WTNDBs1(NB,NZ),CPOOLT,CCNDLB,FXRNX
!    4,CPOLNBs1(NB,NZ),ZPOLNBs1(NB,NZ),PPOLNBs1(NB,NZ)
!    4,CPOOLs1(NB,NZ),ZPOOLs1(NB,NZ),PPOOLs1(NB,NZ)
!    5,WTLSB1,WTNDB1,WTLSBT
!2120  FORMAT(A8,5I4,40E12.4)
!     ENDIF
!     WRITE(*,2121)'NODBAL',I,J,NZ,NB,CPOLNBs1(NB,NZ)
!    2,WTNDBs1(NB,NZ),CPOLNBs1(NB,NZ)+WTNDBs1(NB,NZ)
!    3,RMNDL,RCNDL,RGN2F,CGNDL,RCNDLC,GRNDG,RXNDLC,RXNSNC,RCO2T
!    4,RGNDG,RGNDL,RCNSNC
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  end subroutine CanopyNoduleBiochemistry
!------------------------------------------------------------------------------------------

  subroutine C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: CH2O3(25),CH2O4(25)
  integer :: K
  real(r8) :: CPL4M,CCBS,CPL3K,CO2LK

!     begin_execution

  DO 170 K=1,JNODS1
    IF(WGLFs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     MESOPHYLL TO BUNDLE SHEATH TRANSFER
!
!     WGLF=node leaf C mass
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CH2O3,CH2O4=total CO2 fixation in bundle sheath,mesophyll
!     CPL4M=mesophyll to bundle sheath transfer of nonstructural C4
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!
      CPOOL3s1(K,NB,NZ)=CPOOL3s1(K,NB,NZ)-CH2O3(K)
      CPOOL4s1(K,NB,NZ)=CPOOL4s1(K,NB,NZ)+CH2O4(K)
      CPL4M=1.0_r8*(CPOOL4s1(K,NB,NZ)*WGLFs1(K,NB,NZ)*FBS &
        -CPOOL3s1(K,NB,NZ)*WGLFs1(K,NB,NZ)*FMP) &
        /(WGLFs1(K,NB,NZ)*(FBS+FMP))
      CPOOL4s1(K,NB,NZ)=CPOOL4s1(K,NB,NZ)-CPL4M
      CPOOL3s1(K,NB,NZ)=CPOOL3s1(K,NB,NZ)+CPL4M
!
!     BUNDLE SHEATH CO2 DECARBOXYLATION
!
!     CCBS=CO2 concn in bundle sheath (uM)
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WGLF=node leaf C mass
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!     CPL3K=bundle sheath CO2 decarboxylation
!     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll in C4 (uM)
!     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
!     CPOOL3=C4 nonstructural C mass in bundle sheath
!
      CCBS=AMAX1(0.0,0.083E+09*CO2Bs1(K,NB,NZ)/(WGLFs1(K,NB,NZ)*FBS))
      CPL3K=2.5E-02*CPOOL3s1(K,NB,NZ)/(1.0+CCBS/CO2KI)
      CPOOL3s1(K,NB,NZ)=CPOOL3s1(K,NB,NZ)-CPL3K
      CO2Bs1(K,NB,NZ)=CO2Bs1(K,NB,NZ)+FCO2B*CPL3K
      HCOBs1(K,NB,NZ)=HCOBs1(K,NB,NZ)+FHCOB*CPL3K
!
!     BUNDLE SHEATH LEAKAGE
!
!     CO2LK=bundle sheath CO2 leakage
!     CCBS=CO2 concn in bundle sheath (uM)
!     CO2L=intercellular CO2 concentration (uM)
!     WGLF=node leaf C mass
!     FBS=leaf water content in bundle sheath
!     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
      CO2LK=5.0E-07*(CCBS-CO2Ls1(NZ))*WGLFs1(K,NB,NZ)*FBS
      CO2Bs1(K,NB,NZ)=CO2Bs1(K,NB,NZ)-FCO2B*CO2LK
      HCOBs1(K,NB,NZ)=HCOBs1(K,NB,NZ)-FHCOB*CO2LK
!     IF(NB.EQ.1.AND.K.EQ.14)THEN
!     WRITE(*,6667)'CO2K',I,J,NB,K,CO2LK,CO2Bs1(K,NB,NZ)
!    2,HCOBs1(K,NB,NZ),CPOOL3s1(K,NB,NZ),CH2O3(K),CH2O4(K)
!    3,CCBS,CO2Ls1(NZ),WGLFs1(K,NB,NZ),CPL3Z
!    4,CPL3K,CPL4M,CPOOL4s1(K,NB,NZ),FBS,FMP
!6667  FORMAT(A8,4I4,30E14.6)
!     ENDIF
!
!     TOTAL C EXCHANGE
!
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!     CO2LK=bundle sheath CO2 leakage
!
      TCO2Ts1(NZ)=TCO2Ts1(NZ)-CO2LK
      TCO2As1(NZ)=TCO2As1(NZ)-CO2LK
      CNETs1(NZ)=CNETs1(NZ)-CO2LK
      RECOs1=RECOs1-CO2LK
      TRAUs1=TRAUs1-CO2LK
    ENDIF
170 CONTINUE
  end subroutine C4PhotoProductTransfer
!------------------------------------------------------------------------------------------

  subroutine ComputeGPP(I,J,NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
    WTSHXN,TFN5,WFNG,WFNC,WFNSG,CH2O3,CH2O4,CNPG,rco2c,RMNCS,SNCR,CGROS,CNRDM,CNRDA)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFN6(JZ1)
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,TFN5,WFNG
  real(r8), intent(in) :: WFNC,WFNSG
  real(r8), intent(out) :: rco2c,RMNCS
  real(r8), intent(out) :: CH2O3(25),CH2O4(25)
  REAL(R8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: SNCR
  real(r8), intent(out) :: CGROS
  real(r8), intent(out) :: CNRDM,CNRDA
  real(r8) :: CO2F,ZADDB,PADDB,CH2O
  integer :: K

! begin_execution

! FDBK=N,P feedback inhibition on C3 CO2 fixation
! SSIN=sine of solar angle
! RADP=total PAR absorbed by canopy
! CO2Q=canopy air CO2 concentration
!
  IF(IDAYs1(1,NB,NZ).NE.0)THEN
!   IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!     WRITE(*,5651)'CHECK1',I,J,NZ,NB,IDAYs1(1,NB,NZ)
!    2,FDBKs1(NB,NZ),RADPs1(NZ),CO2Qs1(NZ)
!    3,ARLF1s1(1,NB,NZ)
!5651  FORMAT(A8,5I4,12E12.4)
!   ENDIF

    IF(abs(FDBKs1(NB,NZ)).GT.0)THEN
      IF(SSINs1.GT.0.0.AND.RADPs1(NZ).GT.0.0 &
        .AND.CO2Qs1(NZ).GT.0.0)THEN
        CO2F=0._r8
        CH2O=0._r8
        IF(IGTYPs1(NZ).NE.0.OR.WFNC.GT.0.0)THEN
!
!         FOR EACH NODE
!
          DO 100 K=1,JNODS1
            CH2O3(K)=0._r8
            CH2O4(K)=0._r8
            IF(ARLF1s1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!             C4 PHOTOSYNTHESIS
!
!             ARLF,ARLFL=leaf area
!             ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!             VCGR4=PEP carboxylation rate unlimited by CO2
!
              IF(ICTYPs1(NZ).EQ.4.AND.VCGR4s1(K,NB,NZ).GT.0.0)THEN
!
                CALL ComputeGPP_C4(K,NB,NZ,WFNG,WFNC,CH2O3,CH2O4,CO2F,CH2O)
!
!               C3 PHOTOSYNTHESIS
!
              ELSEIF(ICTYPs1(NZ).NE.4.AND.VCGROs1(K,NB,NZ).GT.0.0)THEN
                call ComputeGPP_C3(K,NB,NZ,WFNG,WFNC,CH2O3,CO2F,CH2O)

              ENDIF
            ENDIF
100       CONTINUE
!
!         CO2F,CH2O=total CO2 fixation,CH2O production
!
          CO2F=CO2F*0.0432
          CH2O=CH2O*0.0432
!
!         CONVERT UMOL M-2 S-1 TO G C M-2 H-1
!
          DO 150 K=1,JNODS1
            CH2O3(K)=CH2O3(K)*0.0432
            CH2O4(K)=CH2O4(K)*0.0432
150       CONTINUE
        ELSE
          CO2F=0._r8
          CH2O=0._r8
          IF(ICTYPs1(NZ).EQ.4)THEN
            DO 155 K=1,JNODS1
              CH2O3(K)=0._r8
              CH2O4(K)=0._r8
155         CONTINUE
          ENDIF
        ENDIF
      ELSE
        CO2F=0._r8
        CH2O=0._r8
        IF(ICTYPs1(NZ).EQ.4)THEN
          DO 160 K=1,JNODS1
            CH2O3(K)=0._r8
            CH2O4(K)=0._r8
160       CONTINUE
        ENDIF
      ENDIF
    ELSE
      CO2F=0._r8
      CH2O=0._r8
      IF(ICTYPs1(NZ).EQ.4)THEN
        DO 165 K=1,JNODS1
          CH2O3(K)=0._r8
          CH2O4(K)=0._r8
165     CONTINUE
      ENDIF
    ENDIF
!
!   SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
!
    call ComputeRAutoAfEmergence(NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
      CO2F,CH2O,TFN5,WFNG,WFNSG,WTSHXN,ZADDB,CNPG,PADDB,RCO2C,RMNCS,SNCR,CGROS,CNRDA)

!   SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
!
  ELSE
    call ComputeRAutoBfEmergence(NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,&
      WFNG,WFNSG,ZADDB,CNPG,PADDB,RCO2C,RMNCS,SNCR,CGROS,CNRDM,CNRDA,CH2O)
  ENDIF

!   REMOVE C,N,P USED IN MAINTENANCE + GROWTH REPIRATION AND GROWTH
!   FROM NON-STRUCTURAL POOLS
!
!   CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!   CH2O=total CH2O production
!   RMNCS=maintenance respiration
!   RCO2C=respiration from non-structural C
!   CGROS=total non-structural C used in growth and respiration
!   CNRDA=respiration for N assimilation
!   ZADDB,PADDB=nonstructural N,P used in growth
!   RNH3B=NH3 flux between atmosphere and branch from uptake.f
!   XFRC,XFRN,XFRP=branch-root layer C,N,P transfer
!
    CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+CH2O-AMIN1(RMNCS,RCO2C)-CGROS-CNRDA
    ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-ZADDB+RNH3Bs1(NB,NZ)
    PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-PADDB
  end subroutine ComputeGPP
!------------------------------------------------------------------------------------------

  subroutine CalcPartitionCoeff(I,J,NB,NZ,part,PTRT,IFLGZ,IFLGY)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  integer, intent(out) :: IFLGZ,IFLGY
  real(r8), intent(out):: PART(7),PTRT
  REAL(R8) :: ARLFI

  integer :: N
  real(r8) :: FPARTL
  real(r8) :: PARTS
  real(r8) :: PARTX
  real(r8) :: TOTAL
  real(r8) :: PSILY(0:3)
  real(r8), parameter :: FPART1=1.00_r8
  real(r8), parameter :: FPART2=0.40_r8

  PSILY=real((/-200.0,-2.0,-2.0,-2.0/),r8)

!     begin_execution
!
!     PARTITION GROWTH WITHIN EACH BRANCH FROM GROWTH STAGE
!     1=LEAF,2=SHEATH OR PETIOLE,3=STALK,4=RESERVE,
!     5,6=REPRODUCTIVE ORGANS,7=GRAIN
!
!     PART=organ partitioning fraction
!

  TOTAL=0._r8
  DO 10 N=1,7
    PART(N)=0._r8
10    CONTINUE
!
!     IF BEFORE FLORAL INDUCTION
!
!     IDAYs1(2,=floral initiation date
!
  IF(IDAYs1(2,NB,NZ).EQ.0)THEN
    PART(1)=0.725
    PART(2)=0.275
!
!     IF BEFORE ANTHESIS
!
!     IDAYs1(6,=start of anthesis and setting final seed number
!     TGSTGI=total change in vegv node number normalized for maturity group
!
  ELSEIF(IDAYs1(6,NB,NZ).EQ.0)THEN
    PART(1)=AMAX1(PART1X,0.725-FPART1*TGSTGIs1(NB,NZ))
    PART(2)=AMAX1(PART2X,0.275-FPART2*TGSTGIs1(NB,NZ))
    PARTS=1.0_r8-PART(1)-PART(2)
    PART(3)=0.60*PARTS
    PART(4)=0.30*PARTS
    PARTX=PARTS-PART(3)-PART(4)
    PART(5)=0.5*PARTX
    PART(6)=0.5*PARTX
!
!     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     IDAYs1(7,=start of grain filling and setting max seed size
!     IDTYP=growth habit:0=determinate,1=indetermimate from PFT file
!     TGSTGF=total change in reprv node number normalized for maturity group
!
  ELSEIF(IDAYs1(7,NB,NZ).EQ.0)THEN
    IF(IDTYPs1(NZ).EQ.0)THEN
      PART(1)=0._r8
      PART(2)=0._r8
    ELSE
      PART(1)=AMAX1(PART1X,(0.725-FPART1)*(1.0-TGSTGFs1(NB,NZ)))
      PART(2)=AMAX1(PART2X,(0.275-FPART2)*(1.0-TGSTGFs1(NB,NZ)))
    ENDIF
    PARTS=1.0_r8-PART(1)-PART(2)
    PART(3)=AMAX1(0.0,0.60*PARTS*(1.0-TGSTGFs1(NB,NZ)))
    PART(4)=AMAX1(0.0,0.30*PARTS*(1.0-TGSTGFs1(NB,NZ)))
    PARTX=PARTS-PART(3)-PART(4)
    PART(5)=0.5*PARTX
    PART(6)=0.5*PARTX
!
!     DURING GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDTYP=growth habit:0=determinate,1=indetermimate
!
  ELSE
    IF(IDTYPs1(NZ).EQ.0)THEN
      PART(7)=1.0_r8
    ELSE
      PART(1)=PART1X
      PART(2)=PART2X
      PARTS=1.0_r8-PART(1)-PART(2)
      IF(ISTYPs1(NZ).EQ.0)THEN
        PART(3)=0.125*PARTS
        PART(5)=0.125*PARTS
        PART(6)=0.125*PARTS
        PART(7)=0.625*PARTS
      ELSE
        PART(3)=0.75*PARTS
        PART(7)=0.25*PARTS
      ENDIF
    ENDIF
  ENDIF
!
!     IF AFTER GRAIN FILLING
!
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IDAYs1(10,=physiological maturity date
!
  IF(IBTYPs1(NZ).EQ.0.AND.IDAYs1(10,NB,NZ).NE.0)THEN
    IF(ISTYPs1(NZ).EQ.0)THEN
      PART(4)=0._r8
      PART(3)=0._r8
      PART(7)=0._r8
    ELSE
      PART(4)=PART(4)+PART(3)
      PART(3)=0._r8
      PART(7)=0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM STALK TO STALK RESERVES IF RESERVES BECOME LOW
!
!     WTRSVB,WVSTKB=stalk reserve,sapwood mass
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!
  IF(IDAYs1(2,NB,NZ).NE.0)THEN
    IF(WTRSVBs1(NB,NZ).LT.XFRX*WVSTKBs1(NB,NZ))THEN
      DO 1020 N=1,7
        IF(N.NE.4)THEN
          PART(4)=PART(4)+0.10*PART(N)
          PART(N)=PART(N)-0.10*PART(N)
        ENDIF
1020  CONTINUE
!
!     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
!
    ELSEIF(WTRSVBs1(NB,NZ).GT.1.0*WVSTKBs1(NB,NZ))THEN
      PART(3)=PART(3)+PART(4)+PART(7)
      PART(4)=0._r8
      PART(7)=0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
!
!     ARLFP=PFT leaf area
!
  ARLFI=ARLFPs1(NZ)/AREA3s1(NUs1)
  IF(ARLFI.GT.5.0)THEN
    FPARTL=AMAX1(0.0,(10.0-ARLFI)/5.0)
    PART(3)=PART(3)+(1.0-FPARTL)*(PART(1)+PART(2))
    PART(1)=FPARTL*PART(1)
    PART(2)=FPARTL*PART(2)
  ENDIF
  IF(NB.EQ.NB1s1(NZ))THEN
    PTRT=PART(1)+PART(2)
  ENDIF
!
!     DECIDUOUS LEAF FALL AFTER GRAIN FILL IN DETERMINATES,
!     AFTER AUTUMNIZATION IN INDETERMINATES, OR AFTER SUSTAINED
!     WATER STRESS
!
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     FVRN=fraction of hours required for leafoff to initiate remobilization
!     IDAYs1(8,=end date for setting final seed number
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IFLGY,IFLGZ=remobilization flags
!     FLGZ=control rate of remobilization
!
  IF((ISTYPs1(NZ).NE.0.AND.VRNFs1(NB,NZ) &
    .GE.FVRNs1(IWTYPs1(NZ))*VRNXs1(NB,NZ)) &
    .OR.(ISTYPs1(NZ).EQ.0 &
    .AND.IDAYs1(8,NB,NZ).NE.0))THEN
    IFLGZ=1
    IF(ISTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.0)THEN
      IFLGY=1
      FLGZs1(NB,NZ)=FLGZs1(NB,NZ)+1.0
    ELSEIF((IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3) &
      .AND.TCCs1(NZ).LT.CTCs1(NZ))THEN
      IFLGY=1
      FLGZs1(NB,NZ)=FLGZs1(NB,NZ)+1.0
    ELSEIF(IWTYPs1(NZ).GE.2 &
      .AND.PSILTs1(NZ).LT.PSILY(IGTYPs1(NZ)))THEN
      IFLGY=1
      FLGZs1(NB,NZ)=FLGZs1(NB,NZ)+1.0
    ENDIF
    IF(ISTYPs1(NZ).NE.0.AND.IWTYPs1(NZ).NE.0)THEN
      PART(3)=PART(3)+0.5*(PART(1)+PART(2))
      PART(4)=PART(4)+0.5*(PART(1)+PART(2))
      PART(1)=0._r8
      PART(2)=0._r8
    ENDIF
  ELSE
    IFLGZ=0
    IFLGY=0
    FLGZs1(NB,NZ)=0._r8
  ENDIF
!
!     CHECK PARTITIONING COEFFICIENTS
!
  DO 1000 N=1,7
    IF(N.EQ.3.AND.test_aeqb(SNL1s1(NZ),0._r8))THEN
      PART(N)=0._r8
    ELSE
      PART(N)=AMAX1(0.0,PART(N))
    ENDIF
    TOTAL=TOTAL+PART(N)
1000  CONTINUE
  IF(TOTAL.GT.ZEROs1)THEN
    DO 1010 N=1,7
      PART(N)=PART(N)/TOTAL
1010  CONTINUE
  ELSE
    DO 1015 N=1,7
      PART(N)=0._r8
1015  CONTINUE
  ENDIF
  end subroutine CalcPartitionCoeff
!------------------------------------------------------------------------------------------

  subroutine GrowOneBranch(I,J,NB,NZ,TFN6,ZCX,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WFNG,WFNC,WFNS,WFNSG,PTRT,UPNFC,IFLGZ)
  implicit none
  integer, intent(in)  :: I,J,NB,NZ
  REAL(R8), INTENT(IN) :: TFN6(JZ1)
  real(r8), intent(in) :: ZCX(JP1)
  real(r8), intent(in) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WFNG,WFNC,WFNS,WFNSG
  real(r8), intent(inout) :: UPNFC(JP1)
  real(r8), intent(out) :: PTRT
  integer, intent(out) :: IFLGZ
  real(r8) :: DMSHD
  integer  :: K,KNOD,KK,kx,K1,K2,KVSTGX,KSNC
  integer  :: KN,MNNOD,NN,MXNOD,M,N,NNOD1
  integer  :: NBK,NBY,NBL,NBX
  real(r8) :: ZPOOLD,XFRN1,XFRP1
  real(r8) :: XKSNC
  REAL(R8) :: ch2o3(25),ch2o4(25)
  integer  :: NBZ(10),IFLGY
  REAL(R8) :: PART(7)
  real(r8) :: ALLOCL,ALLOCS
  REAL(R8) :: ALLOCN
  REAL(R8) :: CNPG
  real(r8) :: CCE
  real(r8) :: CCC,CNC,CPC
  REAL(R8) :: cpoolt
  real(r8) :: DMLFB
  real(r8) :: DMSHB
  real(r8) :: DMSHT
  real(r8) :: CNLFB
  real(r8) :: CPLFB
  real(r8) :: ETOL
  real(r8) :: FSNC
  real(r8) :: FSNCL
  real(r8) :: FSNCS
  real(r8) :: GROGR
  real(r8) :: GROLF
  real(r8) :: GRORSV,GROHSK,GROEAR
  real(r8) :: GROSHT,GROSHE,GROLFN,GROSHN
  real(r8) :: GROSTN,GRORSN,GROHSN
  real(r8) :: GROEAN,GROSTK
  real(r8) :: GROGRN,GROLFP
  real(r8) :: GROSHP,GROSTP
  real(r8) :: GRORSP,GROHSP
  real(r8) :: GROEAP,GROGRP
  real(r8) :: GNOD
  REAL(R8) :: GRO,GRON,GROP
  REAL(R8) :: GSLA,GROA
  real(r8) :: GSSL,GROS
  real(r8) :: GROH
  real(r8) :: PPOOLD,RCO2C
  REAL(R8) :: RMNCS,RCO2V
  real(r8) :: SNCR
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SLA,SSL,SNL
  real(r8) :: SNCZ,SNCX
  real(r8) :: SNCF
  real(r8) :: CNSHB,CPSHB
  real(r8) :: CNLFM,CPLFM
  real(r8) :: CNLFX,CPLFX,CNSHX,CPSHX
  real(r8) :: CGROS
  real(r8) :: CNRDM,CNRDA
  real(r8) :: WTSHXN
  real(r8) :: XFRC,XFRN,XFRP
! begin_execution

  WTLSBs1(NB,NZ)=AMAX1(0.0_r8,WTLFBs1(NB,NZ)+WTSHEBs1(NB,NZ))

  IF(IDTHBs1(NB,NZ).EQ.0)THEN
    call CalcPartitionCoeff(I,J,NB,NZ,PART,PTRT,IFLGY,IFLGZ)
!
!   SHOOT COEFFICIENTS FOR GROWTH RESPIRATION AND N,P CONTENTS
!   FROM GROWTH YIELDS ENTERED IN 'READQ', AND FROM PARTITIONING
!   COEFFICIENTS ABOVE
!
!   DM*B=C production vs nonstructural C consumption
!   CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!   *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
!   *EAR=ear,*GR=grain from PFT file,*SH=shoot
!   DMSHT=branch C production vs nonstructural C consumption
!   DMSHD=branch C respiration vs nonstructural C consumption
!   CN*M,CP*M=min N,P production vs nonstructural C consumption
!   CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
!   CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!
    IF(IDAYs1(1,NB,NZ).NE.0)THEN
      DMLFB=DMLFs1(NZ)
      DMSHB=DMSHEs1(NZ)
      CNLFB=CNLFW
      CNSHB=CNSHW
      CPLFB=CPLFW
      CPSHB=CPSHW
    ELSE
      DMLFB=DMRTs1(NZ)
      DMSHB=DMRTs1(NZ)
      CNLFB=CNRTW
      CNSHB=CNRTW
      CPLFB=CPRTW
      CPSHB=CPRTW
    ENDIF
    DMSHT=PART(1)*DMLFB+PART(2)*DMSHB+PART(3)*DMSTKs1(NZ) &
      +PART(4)*DMRSVs1(NZ)+PART(5)*DMHSKs1(NZ) &
      +PART(6)*DMEARs1(NZ)+PART(7)*DMGRs1(NZ)
    DMSHD=1.0_r8-DMSHT
    CNLFM=PART(1)*DMLFB*ZPLFM*CNLFB
    CPLFM=PART(1)*DMLFB*ZPLFM*CPLFB
    CNLFX=PART(1)*DMLFB*ZPLFD*CNLFB
    CPLFX=PART(1)*DMLFB*ZPLFD*CPLFB
    CNSHX=PART(2)*DMSHB*CNSHB &
      +PART(3)*DMSTKs1(NZ)*CNSTKs1(NZ) &
      +PART(4)*DMRSVs1(NZ)*CNRSVs1(NZ) &
      +PART(5)*DMHSKs1(NZ)*CNHSKs1(NZ) &
      +PART(6)*DMEARs1(NZ)*CNEARs1(NZ) &
      +PART(7)*DMGRs1(NZ)*CNRSVs1(NZ)
    CPSHX=PART(2)*DMSHB*CPSHB &
      +PART(3)*DMSTKs1(NZ)*CPSTKs1(NZ) &
      +PART(4)*DMRSVs1(NZ)*CPRSVs1(NZ) &
      +PART(5)*DMHSKs1(NZ)*CPHSKs1(NZ) &
      +PART(6)*DMEARs1(NZ)*CPEARs1(NZ) &
      +PART(7)*DMGRs1(NZ)*CPRSVs1(NZ)
!
!   TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
!
!   WTSHXN=shoot structural N mass
!   WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,petiole,husk,ear,grain N mass
!   CNSTK,WVSTKB=stalk N:C,sapwood mass
!   IDAYs1(10=date of physiological maturity
!
    WTSHXN=AMAX1(0.0,WTLFBNs1(NB,NZ)+WTSHBNs1(NB,NZ) &
      +CNSTKs1(NZ)*WVSTKBs1(NB,NZ))
    IF(IDAYs1(10,NB,NZ).EQ.0)THEN
      WTSHXN=WTSHXN+AMAX1(0.0,WTHSBNs1(NB,NZ) &
        +WTEABNs1(NB,NZ)+WTGRBNs1(NB,NZ))
    ENDIF
!
!   GROSS PRIMARY PRODUCTIVITY
!
    call ComputeGPP(I,J,NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,TFN5,&
      WFNG,WFNC,WFNSG,CH2O3,CH2O4,CNPG,rco2c,RMNCS,SNCR,CGROS,CNRDM,CNRDA)
!

!
!   TRANSFER OF C4 FIXATION PRODUCTS FROM NON-STRUCTURAL POOLS
!   IN MESOPHYLL TO THOSE IN BUNDLE SHEATH, DECARBOXYLATION
!   OF C4 FIXATION PRODUCTS IN BUNDLE SHEATH, LEAKAGE OF DECARBOXYLATION
!   PRODUCTS BACK TO MESOPHYLL IN C4 PLANTS
!
!   ICTYP=photosynthesis type:3=C3,4=C4
!
    IF(ICTYPs1(NZ).EQ.4)THEN
      call C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
    ENDIF
!
!   C,N,P GROWTH OF LEAF, SHEATH OR PETIOLE, STALK,
!   STALK RESERVES, REPRODUCTIVE ORGANS, GRAIN
!
!   GRO*,GRO*N,GRO*P=organ C,N,P growth rate
!   DM*=C production vs nonstructural C consumption
!   organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
!   HSK=husk,EAR=ear,GR=grain,SHT=shoot
!   PART=organ partitioning fraction
!   CGROS=total non-structural C used in growth and growth respiration
!   CN*,CP*=N:C,P:C ratios in plant organs
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!   CNPG=N,P constraint on growth respiration
!   WT*,WT*N,WT*P=organ C,N,P mass
!
    GROLF=PART(1)*CGROS*DMLFB
    GROSHE=PART(2)*CGROS*DMSHB
    GROSTK=PART(3)*CGROS*DMSTKs1(NZ)
    GRORSV=PART(4)*CGROS*DMRSVs1(NZ)
    GROHSK=PART(5)*CGROS*DMHSKs1(NZ)
    GROEAR=PART(6)*CGROS*DMEARs1(NZ)
    GROGR=PART(7)*CGROS*DMGRs1(NZ)
    GROSHT=CGROS*DMSHT
    GROLFN=GROLF*CNLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHN=GROSHE*CNSHB
    GROSTN=GROSTK*CNSTKs1(NZ)
    GRORSN=GRORSV*CNRSVs1(NZ)
    GROHSN=GROHSK*CNHSKs1(NZ)
    GROEAN=GROEAR*CNEARs1(NZ)
    GROGRN=GROGR*CNRSVs1(NZ)
    GROLFP=GROLF*CPLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHP=GROSHE*CPSHB
    GROSTP=GROSTK*CPSTKs1(NZ)
    GRORSP=GRORSV*CPRSVs1(NZ)
    GROHSP=GROHSK*CPHSKs1(NZ)
    GROEAP=GROEAR*CPEARs1(NZ)
    GROGRP=GROGR*CPRSVs1(NZ)
    WTLFBs1(NB,NZ)=WTLFBs1(NB,NZ)+GROLF
    WTSHEBs1(NB,NZ)=WTSHEBs1(NB,NZ)+GROSHE
    WTSTKBs1(NB,NZ)=WTSTKBs1(NB,NZ)+GROSTK
    WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+GRORSV
    WTHSKBs1(NB,NZ)=WTHSKBs1(NB,NZ)+GROHSK
    WTEARBs1(NB,NZ)=WTEARBs1(NB,NZ)+GROEAR
    WTLFBNs1(NB,NZ)=WTLFBNs1(NB,NZ)+GROLFN
    WTSHBNs1(NB,NZ)=WTSHBNs1(NB,NZ)+GROSHN
    WTSTBNs1(NB,NZ)=WTSTBNs1(NB,NZ)+GROSTN
    WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+GRORSN
    WTHSBNs1(NB,NZ)=WTHSBNs1(NB,NZ)+GROHSN
    WTEABNs1(NB,NZ)=WTEABNs1(NB,NZ)+GROEAN
    WTLFBPs1(NB,NZ)=WTLFBPs1(NB,NZ)+GROLFP
    WTSHBPs1(NB,NZ)=WTSHBPs1(NB,NZ)+GROSHP
    WTSTBPs1(NB,NZ)=WTSTBPs1(NB,NZ)+GROSTP
    WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+GRORSP
    WTHSBPs1(NB,NZ)=WTHSBPs1(NB,NZ)+GROHSP
    WTEABPs1(NB,NZ)=WTEABPs1(NB,NZ)+GROEAP
!
!   ETOLIATION
!
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
!   ETOL=coefficient for etoliation effects on expansion,extension
!
    CCE=AMIN1(safe_adb(CZPOLBs1(NB,NZ),CZPOLBs1(NB,NZ)+CCPOLBs1(NB,NZ)*CNKI) &
      ,safe_adb(CPPOLBs1(NB,NZ),CPPOLBs1(NB,NZ)+CCPOLBs1(NB,NZ)*CPKI))

    ETOL=1.0_r8+CCE
!
!   DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
!
!   MXNOD,MNNOD=max,min node number currently growing
!   KVSTG=integer of most recent leaf number
!   KNOD,GNOD=number of currently growing nodes
!   ALLOCL=fraction of leaf growth allocated to each node
!   GRO,GRON,GROP=leaf C,N,P growth at each node
!   GSLA=allocation of leaf area growth to each node
!   FNOD=scales node number for perennial vegetation (e.g. trees)
!   NNOD=number of concurrently growing nodes
!
    IF(NB.EQ.NB1s1(NZ).AND.HTCTLs1(NZ).LE.SDPTHs1(NZ))THEN
      NNOD1=0
    ELSE
      NNOD1=1
    ENDIF
    IF(GROLF.GT.0.0)THEN
      MXNOD=KVSTGs1(NB,NZ)
      MNNOD=MAX(NNOD1,MXNOD-NNODs1(NZ)+1)
      MXNOD=MAX(MXNOD,MNNOD)
      KNOD=MXNOD-MNNOD+1
      GNOD=KNOD
      ALLOCL=1.0_r8/GNOD
      GRO=ALLOCL*GROLF
      GRON=ALLOCL*GROLFN
      GROP=ALLOCL*GROLFP
      GSLA=ALLOCL*FNODs1(NZ)*NNODs1(NZ)
!
!     GROWTH AT EACH CURRENT NODE
!
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     GRO,GRON,GROP=leaf C,N,P growth at each node
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
      DO 490 KK=MNNOD,MXNOD
        K=MOD(KK,JNODS1)
        IF(K.EQ.0.AND.KK.NE.0)K=25
          WGLFs1(K,NB,NZ)=WGLFs1(K,NB,NZ)+GRO
          WGLFNs1(K,NB,NZ)=WGLFNs1(K,NB,NZ)+GRON
          WGLFPs1(K,NB,NZ)=WGLFPs1(K,NB,NZ)+GROP
          WSLFs1(K,NB,NZ)=WSLFs1(K,NB,NZ) &
            +AMIN1(GRON*CNWSs1(NZ),GROP*CPWSs1(NZ))
!
!         SPECIFIC LEAF AREA FUNCTION OF CURRENT LEAF MASS
!         AT EACH NODE
!
!         SLA=specific area of leaf growth
!         ETOL=coefficient for etoliation effects on expansion,extension
!         SLA1=growth in leaf area vs mass from PFT file
!         SLA2=parameter for calculating leaf area expansion
!         WGLF=leaf C mass
!         PP=PFT population
!         GSLA=allocation of leaf area growth to each node
!         WFNS=turgor expansion,extension function
!         GROA,GRO=leaf area,mass growth
!         ARLFB,ARLF=branch,node leaf area
!
          SLA=ETOL*SLA1s1(NZ)*(AMAX1(ZEROLs1(NZ) &
            ,WGLFs1(K,NB,NZ))/(PPs1(NZ)*GSLA))**SLA2*WFNS
          GROA=GRO*SLA
          ARLFBs1(NB,NZ)=ARLFBs1(NB,NZ)+GROA
          ARLF1s1(K,NB,NZ)=ARLF1s1(K,NB,NZ)+GROA
490     CONTINUE
      ENDIF
!
!     DISTRIBUTE SHEATH OR PETIOLE GROWTH AMONG CURRENTLY GROWING NODES
!
!     MXNOD,MNNOD=max,min node number currently growing
!     KVSTG=integer of most recent leaf number
!     GNOD=number of currently growing nodes
!     ALLOCS=fraction of petiole growth allocated to each node
!     GRO,GRON,GROP=petiole C,N,P growth at each node
!     GSSL=allocation of petiole length growth to each node
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
      IF(GROSHE.GT.0.0)THEN
        MXNOD=KVSTGs1(NB,NZ)
        MNNOD=MAX(NNOD1,MXNOD-NNODs1(NZ)+1)
        MXNOD=MAX(MXNOD,MNNOD)
        GNOD=MXNOD-MNNOD+1
        ALLOCS=1.0_r8/GNOD
        GRO=ALLOCS*GROSHE
        GRON=ALLOCS*GROSHN
        GROP=ALLOCS*GROSHP
        GSSL=ALLOCL*FNODs1(NZ)*NNODs1(NZ)
!
!       GROWTH AT EACH CURRENT NODE
!
!       WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!       GRO,GRON,GROP=petiole C,N,P growth at each node
!       CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
        DO 505 KK=MNNOD,MXNOD
          K=MOD(KK,JNODS1)
          IF(K.EQ.0.AND.KK.NE.0)K=25
            WGSHEs1(K,NB,NZ)=WGSHEs1(K,NB,NZ)+GRO
            WGSHNs1(K,NB,NZ)=WGSHNs1(K,NB,NZ)+GRON
            WGSHPs1(K,NB,NZ)=WGSHPs1(K,NB,NZ)+GROP
            WSSHEs1(K,NB,NZ)=WSSHEs1(K,NB,NZ) &
              +AMIN1(GRON*CNWSs1(NZ),GROP*CPWSs1(NZ))
!
!           SPECIFIC SHEATH OR PETIOLE LENGTH FUNCTION OF CURRENT MASS
!           AT EACH NODE
        !
        !   SSL=specific length of petiole growth
        !   ETOL=coefficient for etoliation effects on expansion,extension
        !   SSL1=growth in petiole length vs mass from PFT file
        !   SSL2=parameter for calculating petiole extension
        !   WGSHE=petiole C mass
        !   PP=PFT population
        !   GSSL=allocation of petiole length growth to each node
        !   WFNS=turgor expansion,extension function
        !   GROS,GRO=petiole length,mass growth
        !   HTSHE=petiole length
!
            IF(WGLFs1(K,NB,NZ).GT.0.0)THEN
              SSL=ETOL*SSL1s1(NZ)*(AMAX1(ZEROLs1(NZ) &
                ,WGSHEs1(K,NB,NZ))/(PPs1(NZ)*GSSL))**SSL2*WFNS
              GROS=GRO/PPs1(NZ)*SSL
              HTSHEs1(K,NB,NZ)=HTSHEs1(K,NB,NZ)+GROS*ANGSHs1(NZ)
        !     IF(I.EQ.120.AND.J.EQ.24)THEN
        !     WRITE(*,2526)'HTSHE',I,J,NZ,NB,K,SSL,WGSHEs1(K,NB,NZ)
        !    2,HTSHEs1(K,NB,NZ),PPs1(NZ),SSL1s1(NZ)
        !    3,GSLA,SSL3,WFNS,GROS,GRO,ANGSHs1(NZ),ZEROLs1(NZ)
        !    4,CCPOLBs1(NB,NZ),ETOL
        !2526  FORMAT(A8,5I4,20E12.4)
       !     ENDIF
            ENDIF
505       CONTINUE
        ENDIF
!
    !   DISTRIBUTE STALK GROWTH AMONG CURRENTLY GROWING NODES
    !
    !   MXNOD,MNNOD=max,min node number currently growing
    !   KVSTG=integer of most recent leaf number
    !   GNOD=number of currently growing nodes
    !   ALLOCN=fraction of stalk growth allocated to each node
    !   GRO,GRON,GROP=stalk C,N,P growth at each node
!
        IF(IDAYs1(1,NB,NZ).EQ.0)THEN
          NN=0
        ELSE
          NN=1
        ENDIF
        MXNOD=KVSTGs1(NB,NZ)
        MNNOD=MAX(MIN(NN,MAX(NN,MXNOD-NNODs1(NZ))),KVSTGs1(NB,NZ)-23)
        MXNOD=MAX(MXNOD,MNNOD)
        IF(GROSTK.GT.0.0)THEN
          GNOD=MXNOD-MNNOD+1
          ALLOCN=1.0_r8/GNOD
          GRO=ALLOCN*GROSTK
          GRON=ALLOCN*GROSTN
          GROP=ALLOCN*GROSTP
    !
    !     SPECIFIC INTERNODE LENGTH FUNCTION OF CURRENT STALK MASS
    !     AT EACH NODE
    !
    !     SNL=specific length of stalk growth
    !     ETOL=coefficient for etoliation effects on expansion,extension
    !     SNL1=growth in stalk length vs mass from PFT file
    !     SNL2=parameter for calculating stalk extension
    !     WTSKB=stalk C mass
    !     PP=PFT population
    !     GROH,GRO=stalk length,mass growth
!
          SNL=ETOL*SNL1s1(NZ)*(WTSTKBs1(NB,NZ)/PPs1(NZ))**SNL2
          GROH=GRO/PPs1(NZ)*SNL
          KX=MOD(MNNOD-1,JNODS1)
          IF(KX.EQ.0.AND.MNNOD-1.NE.0)KX=25
!
    !     GROWTH AT EACH CURRENT NODE
    !
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !     GRO,GRON,GROP=stalk C,N,P growth at each node
    !     HTNODX,HTNODE=stalk height,stalk internode length
    !     ANGBR=sine of stalk angle from horizontal from PFT file
!
          DO 510 KK=MNNOD,MXNOD
            K1=MOD(KK,JNODS1)
            IF(K1.EQ.0.AND.KK.NE.0)K1=25
            K2=MOD(KK-1,JNODS1)
            IF(K2.EQ.0.AND.KK-1.NE.0)K2=25
            WGNODEs1(K1,NB,NZ)=WGNODEs1(K1,NB,NZ)+GRO
            WGNODNs1(K1,NB,NZ)=WGNODNs1(K1,NB,NZ)+GRON
            WGNODPs1(K1,NB,NZ)=WGNODPs1(K1,NB,NZ)+GROP
            HTNODXs1(K1,NB,NZ)=HTNODXs1(K1,NB,NZ)+GROH*ANGBRs1(NZ)
            IF(K1.NE.0)THEN
              HTNODEs1(K1,NB,NZ)=HTNODXs1(K1,NB,NZ) &
                +HTNODEs1(K2,NB,NZ)
            ELSE
              HTNODEs1(K1,NB,NZ)=HTNODXs1(K1,NB,NZ)
            ENDIF

510       CONTINUE
        ENDIF
!
    !   RECOVERY OF REMOBILIZABLE N,P DURING REMOBILIZATION DEPENDS
    !   ON SHOOT NON-STRUCTURAL C:N:P
    !
    !   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
    !   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
    !   RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !   RCCZ,RCCY=min,max fractions for shoot C recycling
    !   RCCX,RCCQ=max fractions for shoot N,P recycling
    !   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
        IF(IDAYs1(1,NB,NZ).NE.0 &
         .AND.CCPOLBs1(NB,NZ).GT.ZEROs1)THEN
          CCC=AMAX1(0.0,AMIN1(1.0 &
            ,safe_adb(CZPOLBs1(NB,NZ),CZPOLBs1(NB,NZ) &
            +CCPOLBs1(NB,NZ)*CNKI) &
            ,safe_adb(CPPOLBs1(NB,NZ),CPPOLBs1(NB,NZ) &
            +CCPOLBs1(NB,NZ)*CPKI)))
          CNC=AMAX1(0.0,AMIN1(1.0 &
            ,safe_adb(CCPOLBs1(NB,NZ),CCPOLBs1(NB,NZ) &
            +CZPOLBs1(NB,NZ)/CNKI)))
          CPC=AMAX1(0.0,AMIN1(1.0 &
            ,safe_adb(CCPOLBs1(NB,NZ),CCPOLBs1(NB,NZ) &
            +CPPOLBs1(NB,NZ)/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        RCCC=RCCZ(IGTYPs1(NZ))+CCC*RCCY(IGTYPs1(NZ))
        RCCN=CNC*RCCX(IGTYPs1(NZ))
        RCCP=CPC*RCCQ(IGTYPs1(NZ))
!
!       WITHDRAW REMOBILIZABLE C,N,P FROM LOWEST NODE AFTER
!       MAXIMUM NODE NUMBER OF 25 IS REACHED
!
!       IFLGG=PFT senescence flag
!       KVSTG=integer of most recent leaf number
!       TFN3=temperature function for canopy growth
!       XRLA=rate of leaf appearance at 25 oC (h-1)
!       FSNC=fraction of lowest leaf to be remobilized
!
        IF(IFLGGs1(NB,NZ).EQ.1)THEN
          KVSTGX=MAX(0,KVSTGs1(NB,NZ)-JNODS1+1)
          IF(KVSTGX.GT.0)THEN
            K=MOD(KVSTGX,JNODS1)
            IF(K.EQ.0.AND.KVSTGX.GT.0)K=JNODS1
            KX=MOD(KVSTGs1(NB,NZ),JNODS1)
            IF(KX.EQ.0.AND.KVSTGs1(NB,NZ).NE.0)KX=JNODS1
            FSNC=TFN3s1(NZ)*XRLAs1(NZ)
!
        !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
        !
        !   IFLGP=flag for remobilization
        !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !   ARLF=node leaf area
        !   RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
!
            IF(IFLGPs1(NB,NZ).EQ.1)THEN
              WGLFXs1(NB,NZ)=AMAX1(0.0,WGLFs1(K,NB,NZ))
              WGLFNXs1(NB,NZ)=AMAX1(0.0,WGLFNs1(K,NB,NZ))
              WGLFPXs1(NB,NZ)=AMAX1(0.0,WGLFPs1(K,NB,NZ))
              ARLFZs1(NB,NZ)=AMAX1(0.0,ARLF1s1(K,NB,NZ))
              IF(WGLFXs1(NB,NZ).GT.ZEROPs1(NZ))THEN
                RCCLXs1(NB,NZ)=RCCC*WGLFXs1(NB,NZ)
                RCZLXs1(NB,NZ)=WGLFNXs1(NB,NZ)*(RCCN+(1.0-RCCN)*RCCC)
                RCPLXs1(NB,NZ)=WGLFPXs1(NB,NZ)*(RCCP+(1.0-RCCP)*RCCC)
              ELSE
                RCCLXs1(NB,NZ)=0._r8
                RCZLXs1(NB,NZ)=0._r8
                RCPLXs1(NB,NZ)=0._r8
              ENDIF
            ENDIF
!
    !       FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !       FSNC,FSNCL=fraction of lowest leaf to be remobilized
    !
            IF(FSNC*WGLFXs1(NB,NZ).GT.WGLFs1(K,NB,NZ) &
              .AND.WGLFXs1(NB,NZ).GT.ZEROPs1(NZ))THEN
              FSNCL=AMAX1(0.0,WGLFs1(K,NB,NZ)/WGLFXs1(NB,NZ))
            ELSE
              FSNCL=FSNC
            ENDIF
!
    !       NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !       TO FRACTIONS SET IN 'STARTQ'
    !
    !       CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !       CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !       FSNCL=fraction of lowest leaf to be remobilized
    !       RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
    !       WGLFX,WGLFNX,WGLFPX=senescing leaf C,N,P mass
    !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !       FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
            DO 6300 M=1,4
              CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(5,M,NZ) &
                *FSNCL*(WGLFXs1(NB,NZ)-RCCLXs1(NB,NZ))*FWODBs1(0)
              ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(5,M,NZ) &
                *FSNCL*(WGLFNXs1(NB,NZ)-RCZLXs1(NB,NZ))*FWODLNs1(0)
              PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(5,M,NZ) &
                *FSNCL*(WGLFPXs1(NB,NZ)-RCPLXs1(NB,NZ))*FWODLPs1(0)
              CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(1,M,NZ) &
                *FSNCL*(WGLFXs1(NB,NZ)-RCCLXs1(NB,NZ))*FWODBs1(1)
              ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(1,M,NZ) &
                *FSNCL*(WGLFNXs1(NB,NZ)-RCZLXs1(NB,NZ))*FWODLNs1(1)
              PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(1,M,NZ) &
                *FSNCL*(WGLFPXs1(NB,NZ)-RCPLXs1(NB,NZ))*FWODLPs1(1)
6300        CONTINUE
!
    !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !       FSNCL=fraction of lowest leaf to be remobilized
    !       ARLFB,ARLFZ=branch living,senescing leaf area
    !       WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in living,senescing leaf
    !       WSLF=leaf protein mass
    !       CNWS,CPWS=protein:N,protein:P ratios from startq.f
    !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
    !       RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
!
            ARLFBs1(NB,NZ)=ARLFBs1(NB,NZ) &
              -FSNCL*ARLFZs1(NB,NZ)
            WTLFBs1(NB,NZ)=WTLFBs1(NB,NZ) &
              -FSNCL*WGLFXs1(NB,NZ)
            WTLFBNs1(NB,NZ)=WTLFBNs1(NB,NZ) &
              -FSNCL*WGLFNXs1(NB,NZ)
            WTLFBPs1(NB,NZ)=WTLFBPs1(NB,NZ) &
              -FSNCL*WGLFPXs1(NB,NZ)
            ARLF1s1(K,NB,NZ)=ARLF1s1(K,NB,NZ) &
              -FSNCL*ARLFZs1(NB,NZ)
            WGLFs1(K,NB,NZ)=WGLFs1(K,NB,NZ) &
              -FSNCL*WGLFXs1(NB,NZ)
            WGLFNs1(K,NB,NZ)=WGLFNs1(K,NB,NZ) &
              -FSNCL*WGLFNXs1(NB,NZ)
            WGLFPs1(K,NB,NZ)=WGLFPs1(K,NB,NZ) &
              -FSNCL*WGLFPXs1(NB,NZ)
            WSLFs1(K,NB,NZ)=AMAX1(0.0,WSLFs1(K,NB,NZ) &
              -FSNCL*AMAX1(WGLFNXs1(NB,NZ)*CNWSs1(NZ) &
              ,WGLFPXs1(NB,NZ)*CPWSs1(NZ)))
            CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+FSNCL*RCCLXs1(NB,NZ)
            ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+FSNCL*RCZLXs1(NB,NZ)
            PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+FSNCL*RCPLXs1(NB,NZ)
!
    !       REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P ALSO DEPENDS ON
    !       STRUCTURAL C:N:P
    !
    !       IFLGP=flag for remobilization
    !       WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !       HTSHEX=petiole length
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
!
            IF(IFLGPs1(NB,NZ).EQ.1)THEN
              WGSHEXs1(NB,NZ)=AMAX1(0.0,WGSHEs1(K,NB,NZ))
              WGSHNXs1(NB,NZ)=AMAX1(0.0,WGSHNs1(K,NB,NZ))
              WGSHPXs1(NB,NZ)=AMAX1(0.0,WGSHPs1(K,NB,NZ))
              HTSHEXs1(NB,NZ)=AMAX1(0.0,HTSHEs1(K,NB,NZ))
              IF(WGSHEXs1(NB,NZ).GT.ZEROPs1(NZ))THEN
                RCCSXs1(NB,NZ)=RCCC*WGSHEXs1(NB,NZ)
                RCZSXs1(NB,NZ)=WGSHNXs1(NB,NZ) &
                  *(RCCN+(1.0-RCCN)*RCCSXs1(NB,NZ)/WGSHEXs1(NB,NZ))
                RCPSXs1(NB,NZ)=WGSHPXs1(NB,NZ) &
                  *(RCCP+(1.0-RCCP)*RCCSXs1(NB,NZ)/WGSHEXs1(NB,NZ))
              ELSE
                RCCSXs1(NB,NZ)=0._r8
                RCZSXs1(NB,NZ)=0._r8
                RCPSXs1(NB,NZ)=0._r8
              ENDIF
              WTSTXBs1(NB,NZ)=WTSTXBs1(NB,NZ)+WGNODEs1(K,NB,NZ)
              WTSTXNs1(NB,NZ)=WTSTXNs1(NB,NZ)+WGNODNs1(K,NB,NZ)
              WTSTXPs1(NB,NZ)=WTSTXPs1(NB,NZ)+WGNODPs1(K,NB,NZ)
              WGNODEs1(K,NB,NZ)=0._r8
              WGNODNs1(K,NB,NZ)=0._r8
              WGNODPs1(K,NB,NZ)=0._r8
              HTNODXs1(K,NB,NZ)=0._r8
            ENDIF
!
    !       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !
            IF(FSNC*WGSHEXs1(NB,NZ).GT.WGSHEs1(K,NB,NZ) &
              .AND.WGSHEXs1(NB,NZ).GT.ZEROPs1(NZ))THEN
              FSNCS=AMAX1(0.0,WGSHEs1(K,NB,NZ)/WGSHEXs1(NB,NZ))
            ELSE
              FSNCS=FSNC
            ENDIF
!
    !       NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !       TO FRACTIONS SET IN 'STARTQ'
    !
    !       CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !       CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !       FSNCS=fraction of lowest petiole to be remobilized
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
    !       WGSHX,WGSHNX,WGSHPX=senescing petiole C,N,P mass
    !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !       FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!
            DO 6305 M=1,4
              CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(5,M,NZ) &
                *FSNCS*(WGSHEXs1(NB,NZ)-RCCSXs1(NB,NZ))*FWODBs1(0)
              ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(5,M,NZ) &
                *FSNCS*(WGSHNXs1(NB,NZ)-RCZSXs1(NB,NZ))*FWODSNs1(0)
              PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(5,M,NZ) &
                *FSNCS*(WGSHPXs1(NB,NZ)-RCPSXs1(NB,NZ))*FWODSPs1(0)
              CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(2,M,NZ) &
                *FSNCS*(WGSHEXs1(NB,NZ)-RCCSXs1(NB,NZ))*FWODBs1(1)
              ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(2,M,NZ) &
                *FSNCS*(WGSHNXs1(NB,NZ)-RCZSXs1(NB,NZ))*FWODSNs1(1)
              PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(2,M,NZ) &
                *FSNCS*(WGSHPXs1(NB,NZ)-RCPSXs1(NB,NZ))*FWODSPs1(1)
6305        CONTINUE
!
    !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !       HTSHE,HTSHEX=living,senescing petiole length
    !       WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in living,senescing petiole
    !       WSSHE=petiole protein mass
    !       CNWS,CPWS=protein:N,protein:P ratios from startq.f
    !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
!
            WTSHEBs1(NB,NZ)=WTSHEBs1(NB,NZ) &
              -FSNCS*WGSHEXs1(NB,NZ)
            WTSHBNs1(NB,NZ)=WTSHBNs1(NB,NZ) &
              -FSNCS*WGSHNXs1(NB,NZ)
            WTSHBPs1(NB,NZ)=WTSHBPs1(NB,NZ) &
              -FSNCS*WGSHPXs1(NB,NZ)
            HTSHEs1(K,NB,NZ)=HTSHEs1(K,NB,NZ) &
              -FSNCS*HTSHEXs1(NB,NZ)
            WGSHEs1(K,NB,NZ)=WGSHEs1(K,NB,NZ) &
              -FSNCS*WGSHEXs1(NB,NZ)
            WGSHNs1(K,NB,NZ)=WGSHNs1(K,NB,NZ) &
              -FSNCS*WGSHNXs1(NB,NZ)
            WGSHPs1(K,NB,NZ)=WGSHPs1(K,NB,NZ) &
              -FSNCS*WGSHPXs1(NB,NZ)
            WSSHEs1(K,NB,NZ)=AMAX1(0.0,WSSHEs1(K,NB,NZ) &
              -FSNCS*AMAX1(WGSHNXs1(NB,NZ)*CNWSs1(NZ) &
              ,WGSHPXs1(NB,NZ)*CPWSs1(NZ)))
            CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+FSNCS*RCCSXs1(NB,NZ)
            ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+FSNCS*RCZSXs1(NB,NZ)
            PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+FSNCS*RCPSXs1(NB,NZ)
          ENDIF
        ENDIF
!
  !     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0
  !
  !     SNCR=excess maintenance respiration
  !     WTRSVB=stalk reserve C mass
  !     RCO2V=remobilization of stalk reserve C
  !     VMXC=rate constant for nonstructural C oxidation in respiration
  !     TFN3=temperature function for canopy growth
  !
        IF(IFLGZ.EQ.0)THEN
          IF(SNCR.GT.0.0.AND.WTRSVBs1(NB,NZ).GT.0.0)THEN
            RCO2V=AMIN1(SNCR,VMXC*WTRSVBs1(NB,NZ)*TFN3s1(NZ))
            WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)-RCO2V
            SNCR=SNCR-RCO2V
          ENDIF
        ENDIF
!
!       TOTAL REMOBILIZATION = GROWTH RESPIRATION < 0 + DECIDUOUS LEAF
!       FALL DURING AUTUMN + REMOBILZATION DURING GRAIN FILL IN ANNUALS
!
!       ISTYP=growth habit:0=annual,1=perennial from PFT file
!       IFLGY,IFLGZ=remobilization flags
!       SNCZ=phenologically-driven respiration senescence during late-season
!       FXFB=rate constant for plant-storage nonstructural C,N,P exchange
!       IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!       WTLSB=leaf+petiole mass
!       FLGZ=control rate of remobilization
!       FLGZX=number of hours until full senescence after physl maturity
!       SNCX=total senescence respiration
!       KVSTG,KVSTGN=integer of highest,lowest leaf number currently growing
!       KSNC=number of nodes undergoing remobilization
!       SNCF=ratio of phenologically-driven vs total senescence respiration
!
        IF(IFLGZ.EQ.1.AND.IFLGY.EQ.1.AND.ISTYPs1(NZ).NE.0)THEN
          SNCZ=FXFB(IBTYPs1(NZ)) &
            *WTLSBs1(NB,NZ)*AMIN1(1.0,FLGZs1(NB,NZ)/FLGZX)
        ELSE
          SNCZ=0._r8
        ENDIF
        SNCX=SNCR+SNCZ
        IF(SNCX.GT.ZEROPs1(NZ))THEN
          SNCF=SNCZ/SNCX
          KSNC=INT(0.5*(KVSTGs1(NB,NZ)-KVSTGNs1(NB,NZ)))+1
          XKSNC=KSNC
          KN=MAX(0,KVSTGNs1(NB,NZ)-1)
    !
    !     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
    !     IF MAIN STEM POOLS ARE DEPLETED
    !
    !     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
    !     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
    !
          IF(IBTYPs1(NZ).NE.0.AND.IGTYPs1(NZ).GT.1 &
            .AND.NB.EQ.NB1s1(NZ).AND.test_aeqb(SNCF,0._r8))THEN
            NBY=0
          DO 584 NBL=1,NBRs1(NZ)
            NBZ(NBL)=0
584       CONTINUE
          DO 586 NBL=1,NBRs1(NZ)
            NBX=KVSTGs1(NB,NZ)
            DO 585 NBK=1,NBRs1(NZ)
              IF(IDTHBs1(NBK,NZ).EQ.0.AND.NBK.NE.NB1s1(NZ) &
                .AND.NBTBs1(NBK,NZ).LT.NBX &
                .AND.NBTBs1(NBK,NZ).GT.NBY)THEN
                NBZ(NBL)=NBK
                NBX=NBTBs1(NBK,NZ)
              ENDIF
585         CONTINUE
            IF(NBZ(NBL).NE.0)THEN
              NBY=NBTBs1(NBZ(NBL),NZ)
            ENDIF
586       CONTINUE
          DO 580 NBL=1,NBRs1(NZ)
            IF(NBZ(NBL).NE.0)THEN
              IF(NBTBs1(NBZ(NBL),NZ).LT.KK)THEN
                IF(CPOOLs1(NBZ(NBL),NZ).GT.0)THEN
                  XFRC=1.0E-02_r8*AMIN1(SNCX,CPOOLs1(NBZ(NBL),NZ))
                  XFRN=XFRC*ZPOOLs1(NBZ(NBL),NZ)/CPOOLs1(NBZ(NBL),NZ)
                  XFRP=XFRC*PPOOLs1(NBZ(NBL),NZ)/CPOOLs1(NBZ(NBL),NZ)
                ELSE
                  XFRC=0._r8
                  XFRN=1.0E-02_r8*ZPOOLs1(NBZ(NBL),NZ)
                  XFRP=1.0E-02_r8*PPOOLs1(NBZ(NBL),NZ)
                ENDIF
                CPOOLs1(NBZ(NBL),NZ)=CPOOLs1(NBZ(NBL),NZ)-XFRC
                ZPOOLs1(NBZ(NBL),NZ)=ZPOOLs1(NBZ(NBL),NZ)-XFRN
                PPOOLs1(NBZ(NBL),NZ)=PPOOLs1(NBZ(NBL),NZ)-XFRP
                CPOOLs1(NB1s1(NZ),NZ)=CPOOLs1(NB1s1(NZ),NZ) &
                  +XFRC*SNCF
                ZPOOLs1(NB1s1(NZ),NZ)=ZPOOLs1(NB1s1(NZ),NZ) &
                  +XFRN
                PPOOLs1(NB1s1(NZ),NZ)=PPOOLs1(NB1s1(NZ),NZ) &
                  +XFRP
                SNCX=SNCX-XFRC
                IF(SNCX.LE.0.0)exit
              ENDIF
            ENDIF
580       CONTINUE
        ENDIF
!
        IF(SNCX.GT.0.0)call RemobilizeLeafLayers(KN,KSNC,NB,nz,XKSNC,SNCX,RCCC,RCCN,RCCP,SNCF)
    ENDIF

!
!   DEATH IF MAIN STALK OF TREE DIES
!
!   IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   IDTHB=branch living flag: 0=alive,1=dead
!   KVSTGX,KVSTG=integer of lowest,highest leaf number currently growing
!   WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   CNLFB,CPLFB=N:C,P:C ratios in leaf
!
    IF(IBTYPs1(NZ).NE.0.AND.IGTYPs1(NZ).GT.1 &
      .AND.IDTHBs1(NB1s1(NZ),NZ).EQ.1)IDTHBs1(NB,NZ)=1
!
!     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
!
    KVSTGX=MAX(0,KVSTGs1(NB,NZ)-24)
    DO 495 KK=KVSTGX,KVSTGs1(NB,NZ)
      K=MOD(KK,JNODS1)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      IF(WGLFs1(K,NB,NZ).GT.0.0)THEN
        CPOOLT=WGLFs1(K,NB,NZ)+CPOOLs1(NB,NZ)
        IF(CPOOLT.GT.ZEROPs1(NZ))THEN
          ZPOOLD=WGLFNs1(K,NB,NZ)*CPOOLs1(NB,NZ) &
            -ZPOOLs1(NB,NZ)*WGLFs1(K,NB,NZ)
          XFRN1=AMAX1(0.0,AMIN1(1.0E-03*ZPOOLD/CPOOLT,WGLFNs1(K,NB,NZ) &
            -ZPLFM*CNLFB*WGLFs1(K,NB,NZ)))
          PPOOLD=WGLFPs1(K,NB,NZ)*CPOOLs1(NB,NZ) &
            -PPOOLs1(NB,NZ)*WGLFs1(K,NB,NZ)
          XFRP1=AMAX1(0.0,AMIN1(1.0E-03*PPOOLD/CPOOLT,WGLFPs1(K,NB,NZ) &
            -ZPLFM*CPLFB*WGLFs1(K,NB,NZ)))
          XFRN=AMAX1(XFRN1,10.0*XFRP1)
          XFRP=AMAX1(XFRP1,0.10*XFRN1)
          WGLFNs1(K,NB,NZ)=WGLFNs1(K,NB,NZ)-XFRN
          WTLFBNs1(NB,NZ)=WTLFBNs1(NB,NZ)-XFRN
          ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+XFRN
          WGLFPs1(K,NB,NZ)=WGLFPs1(K,NB,NZ)-XFRP
          WTLFBPs1(NB,NZ)=WTLFBPs1(NB,NZ)-XFRP
          PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+XFRP
          WSLFs1(K,NB,NZ)=AMAX1(0.0,WSLFs1(K,NB,NZ) &
            -AMAX1(XFRN*CNWSs1(NZ),XFRP*CPWSs1(NZ)))
        ENDIF
      ENDIF
495 CONTINUE
!   KK inherits value from loop 495, is it right?
    call AllocateLeafToCanopyLayers(NB,NZ,ZCX)

!
    !     ALLOCATE LEAF AREA TO INCLINATION CLASSES ACCORDING TO
    !     DISTRIBUTION ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
    !
    !     SSIN=sine of solar angle
    !     SURF=leaf node surface area in canopy layer
    !     ARLF,ARLFL=leaf node surface area in canopy layer
    !     ZC,DPTHS=canopy,snowpack height
    !     CLASS=leaf inclination class
    !
    IF(SSINNs1.GT.0.0)THEN
      call LeafClassAllocation(NB,NZ)
    ENDIF

    call GrainFilling(I,NB,NZ,GROGR,GROSTK,GROGRN,GROGRP)
!
    call PhenologyReset(I,NB,NZ)
!
    call CarbNutInBranchTransfer(I,J,NB,NZ,IFLGZ,WFNG,WFNSG)

!   CANOPY N2 FIXATION (CYANOBACTERIA)
!
    call CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WFNG,UPNFC)
  ENDIF
  end subroutine GrowOneBranch
!------------------------------------------------------------------------------------------
  subroutine RemobilizeLeafLayers(KN,KSNC,NB,nz,XKSNC,SNCX,RCCC,RCCN,RCCP,SNCF)
  implicit none
  integer, intent(inout) :: KN
  INTEGER, intent(in)    :: nb,nz,KSNC
  real(r8), intent(in)   :: XKSNC,SNCX
  REAL(R8), INTENT(IN) :: RCCC,RCCN,RCCP
  real(r8), intent(inout):: SNCF
  integer :: N,M,K,KK,MXNOD,MNNOD
  real(r8) :: FSNCL,FSNCS
  real(r8) :: FNCLF,FSNAL,FSNAS,FRCC
  real(r8) :: FSNCK
  real(r8) :: FSNCR
  real(r8) :: HTNODZ
  real(r8) :: RCCL,RCZL,RCPL
  real(r8) :: RCCS,RCZS,RCPS
  real(r8) :: RCSC,RCSN,RCSP
  real(r8) :: RCCK,RCZK,RCPK
  real(r8) :: SNCR
  real(r8) :: SNCZ
  real(r8) :: SNCT
  real(r8) :: SNCLF
  real(r8) :: SNCSH
! begin_execution
!     REMOBILIZATION AND LITTERFALL WHEN GROWTH RESPIRATION < 0
!     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
!
!     SNCX,SNCT=branch,node senescence respiration
!     KSNC=number of nodes undergoing remobilization
!
!     IF(NZ.EQ.2.OR.NZ.EQ.3)THEN
!     WRITE(*,1266)'SNCX1',I,J,NX,NY,NZ,NB,SNCY,SNCR,SNCX,SNCF
!    2,CPOOLs1(NB,NZ),WTLSBs1(NB,NZ),RCCC
!     ENDIF

  DO 575 N=1,KSNC
    SNCT=SNCX/XKSNC
    DO 650 KK=KN,KVSTGs1(NB,NZ)
      SNCLF=0._r8
      SNCSH=0._r8
      K=MOD(KK,JNODS1)
      IF(K.EQ.0.AND.KK.NE.0)K=25
!
!       REMOBILIZATION OF LEAF C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
!
!       WGLF,WGSHE=node leaf,petiole C mass
    !       SCNF,SCNSH=leaf,petiole senescence respiration
    !       RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
    !       RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !       RCCZ,RCCY=min,max fractions for shoot C recycling
    !
      IF(WGLFs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
        FNCLF=WGLFs1(K,NB,NZ)/(WGLFs1(K,NB,NZ)+WGSHEs1(K,NB,NZ))
        SNCLF=FNCLF*SNCT
        SNCSH=SNCT-SNCLF
        RCCL=RCCC*WGLFs1(K,NB,NZ)
        RCZL=WGLFNs1(K,NB,NZ)*(RCCN+(1.0-RCCN)*RCCC)
        RCPL=WGLFPs1(K,NB,NZ)*(RCCP+(1.0-RCCP)*RCCC)
!
    !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !         FSNCL,FSNAL=fraction of current leaf C,area to be remobilized
    !
        IF(RCCL.GT.ZEROPs1(NZ))THEN
          FSNCL=AMAX1(0.0,AMIN1(1.0,SNCLF/RCCL))
        ELSE
          FSNCL=1.0_r8
        ENDIF
        FSNAL=FSNCL
!
        !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
        !     TO FRACTIONS SET IN 'STARTQ'
        !
        !     CSNC,ZSNC,PSNC=literfall C,N,P
        !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
        !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
        !     FSNCL=fraction of current leaf to be remobilized
        !     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
        !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
        !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
        !     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
        !
        DO 6310 M=1,4
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(5,M,NZ) &
            *FSNCL*(WGLFs1(K,NB,NZ)-RCCL)*FWODBs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(5,M,NZ) &
            *FSNCL*(WGLFNs1(K,NB,NZ)-RCZL)*FWODLNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(5,M,NZ) &
            *FSNCL*(WGLFPs1(K,NB,NZ)-RCPL)*FWODLPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(1,M,NZ) &
            *FSNCL*(WGLFs1(K,NB,NZ)-RCCL)*FWODBs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(1,M,NZ) &
            *FSNCL*(WGLFNs1(K,NB,NZ)-RCZL)*FWODLNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(1,M,NZ) &
            *FSNCL*(WGLFPs1(K,NB,NZ)-RCPL)*FWODLPs1(1)
6310    CONTINUE
        IF(K.NE.0)THEN
          CSNCs1(2,1,0,NZ)=CSNCs1(2,1,0,NZ) &
            +FSNCL*(CPOOL3s1(K,NB,NZ)+CPOOL4s1(K,NB,NZ))
          CPOOL3s1(K,NB,NZ)=CPOOL3s1(K,NB,NZ) &
            -FSNCL*CPOOL3s1(K,NB,NZ)
          CPOOL4s1(K,NB,NZ)=CPOOL4s1(K,NB,NZ) &
            -FSNCL*CPOOL4s1(K,NB,NZ)
        ENDIF
!
!     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
!
!     ARLFB=leaf area
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     FSNCL=fraction of current leaf to be remobilized
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
        ARLFBs1(NB,NZ)=AMAX1(0.0,ARLFBs1(NB,NZ)-FSNAL*ARLF1s1(K,NB,NZ))
        WTLFBs1(NB,NZ)=AMAX1(0.0,WTLFBs1(NB,NZ)-FSNCL*WGLFs1(K,NB,NZ))
        WTLFBNs1(NB,NZ)=AMAX1(0.0,WTLFBNs1(NB,NZ)-FSNCL*WGLFNs1(K,NB,NZ))
        WTLFBPs1(NB,NZ)=AMAX1(0.0,WTLFBPs1(NB,NZ)-FSNCL*WGLFPs1(K,NB,NZ))
        ARLF1s1(K,NB,NZ)=ARLF1s1(K,NB,NZ)-FSNAL*ARLF1s1(K,NB,NZ)
        WGLFs1(K,NB,NZ)=WGLFs1(K,NB,NZ)-FSNCL*WGLFs1(K,NB,NZ)
        WGLFNs1(K,NB,NZ)=WGLFNs1(K,NB,NZ)-FSNCL*WGLFNs1(K,NB,NZ)
        WGLFPs1(K,NB,NZ)=WGLFPs1(K,NB,NZ)-FSNCL*WGLFPs1(K,NB,NZ)
        WSLFs1(K,NB,NZ)=AMAX1(0.0,WSLFs1(K,NB,NZ) &
          -FSNCL*AMAX1(WGLFNs1(K,NB,NZ)*CNWSs1(NZ) &
          ,WGLFPs1(K,NB,NZ)*CPWSs1(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
!     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     FSNCL=fraction of current leaf C to be remobilized
!     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
!     SNCLF,SNCT=remaining senescence respiration carried to next node
!
        CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+FSNCL*RCCL*SNCF
        ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+FSNCL*RCZL
        PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+FSNCL*RCPL
        SNCLF=SNCLF-FSNCL*RCCL
        SNCT=SNCT-FSNCL*RCCL
        IF(WTLFBs1(NB,NZ).LE.ZEROLs1(NZ))THEN
          WTLFBs1(NB,NZ)=0._r8
          ARLFBs1(NB,NZ)=0._r8
        ENDIF
      !
        !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
        !
!        IF(SNCLF.LE.ZEROPs1(NZ))GO TO 564
        !
        !     OTHERWISE REMAINING C,N,P IN LEAF GOES TO LITTERFALL
        !
        !     CSNC,ZSNC,PSNC=literfall C,N,P
        !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
        !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
        !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
        !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
        !     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
        !     ARLFB=leaf area
        !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
        !     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
        !
      ELSE
        DO 6315 M=1,4
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(5,M,NZ) &
            *WGLFs1(K,NB,NZ)*FWODBs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(5,M,NZ) &
            *WGLFNs1(K,NB,NZ)*FWODLNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(5,M,NZ) &
            *WGLFPs1(K,NB,NZ)*FWODLPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(1,M,NZ) &
            *WGLFs1(K,NB,NZ)*FWODBs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(1,M,NZ) &
            *WGLFNs1(K,NB,NZ)*FWODLNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(1,M,NZ) &
            *WGLFPs1(K,NB,NZ)*FWODLPs1(1)
6315    CONTINUE
        IF(K.NE.0)THEN
          CSNCs1(2,1,0,NZ)=CSNCs1(2,1,0,NZ)+CPOOL3s1(K,NB,NZ)+CPOOL4s1(K,NB,NZ)
          CPOOL3s1(K,NB,NZ)=0._r8
          CPOOL4s1(K,NB,NZ)=0._r8
        ENDIF
        ARLFBs1(NB,NZ)=AMAX1(0.0,ARLFBs1(NB,NZ)-ARLF1s1(K,NB,NZ))
        WTLFBs1(NB,NZ)=AMAX1(0.0,WTLFBs1(NB,NZ)-WGLFs1(K,NB,NZ))
        WTLFBNs1(NB,NZ)=AMAX1(0.0,WTLFBNs1(NB,NZ)-WGLFNs1(K,NB,NZ))
        WTLFBPs1(NB,NZ)=AMAX1(0.0,WTLFBPs1(NB,NZ)-WGLFPs1(K,NB,NZ))
        ARLF1s1(K,NB,NZ)=0._r8
        WGLFs1(K,NB,NZ)=0._r8
        WGLFNs1(K,NB,NZ)=0._r8
        WGLFPs1(K,NB,NZ)=0._r8
        WSLFs1(K,NB,NZ)=0._r8
        IF(WTLFBs1(NB,NZ).LE.ZEROLs1(NZ))THEN
          WTLFBs1(NB,NZ)=0._r8
          ARLFBs1(NB,NZ)=0._r8
        ENDIF
      ENDIF
!564   CONTINUE
!
    !     REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P DEPENDS ON
    !     NON-STRUCTURAL C:N:P
    !
    !     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !
      IF(WGSHEs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
        RCCS=RCCC*WGSHEs1(K,NB,NZ)
        RCZS=WGSHNs1(K,NB,NZ)*(RCCN+(1.0-RCCN)*RCCC)
        RCPS=WGSHPs1(K,NB,NZ)*(RCCP+(1.0-RCCP)*RCCC)
!
      !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
      !     OR PETIOLE
      !
      !     FSNCS,FSNAS=fraction of current petiole C,length to be remobilized
      !
        IF(RCCS.GT.ZEROPs1(NZ))THEN
          FSNCS=AMAX1(0.0,AMIN1(1.0,SNCSH/RCCS))
        ELSE
          FSNCS=1.0_r8
        ENDIF
        FSNAS=FSNCS
    !
      !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FSNCS=fraction of current petiole to be remobilized
      !     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
      !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
      !
        DO 6320 M=1,4
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(5,M,NZ) &
            *FSNCS*(WGSHEs1(K,NB,NZ)-RCCS)*FWODBs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(5,M,NZ) &
            *FSNCS*(WGSHNs1(K,NB,NZ)-RCZS)*FWODSNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(5,M,NZ) &
            *FSNCS*(WGSHPs1(K,NB,NZ)-RCPS)*FWODSPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(2,M,NZ) &
            *FSNCS*(WGSHEs1(K,NB,NZ)-RCCS)*FWODBs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(2,M,NZ) &
            *FSNCS*(WGSHNs1(K,NB,NZ)-RCZS)*FWODSNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(2,M,NZ) &
            *FSNCS*(WGSHPs1(K,NB,NZ)-RCPS)*FWODSPs1(1)
6320    CONTINUE
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     HTSHE=petiole length
      !     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !     FSNCS=fraction of current petiole to be remobilized
      !     CNWS,CPWS=protein:N,protein:P ratios from startq.f
      !
        WTSHEBs1(NB,NZ)=AMAX1(0.0,WTSHEBs1(NB,NZ) &
          -FSNCS*WGSHEs1(K,NB,NZ))
        WTSHBNs1(NB,NZ)=AMAX1(0.0,WTSHBNs1(NB,NZ) &
          -FSNCS*WGSHNs1(K,NB,NZ))
        WTSHBPs1(NB,NZ)=AMAX1(0.0,WTSHBPs1(NB,NZ) &
          -FSNCS*WGSHPs1(K,NB,NZ))
        HTSHEs1(K,NB,NZ)=HTSHEs1(K,NB,NZ) &
          -FSNAS*HTSHEs1(K,NB,NZ)
        WGSHEs1(K,NB,NZ)=WGSHEs1(K,NB,NZ) &
          -FSNCS*WGSHEs1(K,NB,NZ)
        WGSHNs1(K,NB,NZ)=WGSHNs1(K,NB,NZ) &
          -FSNCS*WGSHNs1(K,NB,NZ)
        WGSHPs1(K,NB,NZ)=WGSHPs1(K,NB,NZ) &
          -FSNCS*WGSHPs1(K,NB,NZ)
        WSSHEs1(K,NB,NZ)=AMAX1(0.0,WSSHEs1(K,NB,NZ) &
          -FSNCS*AMAX1(WGSHNs1(K,NB,NZ)*CNWSs1(NZ) &
          ,WGSHPs1(K,NB,NZ)*CPWSs1(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
  !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
  !
  !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  !     FSNCS=fraction of current petiole C to be remobilized
  !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
  !     SNCSH,SNCT=remaining senescence respiration carried to next node
  !
        CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+FSNCS*RCCS*SNCF
        ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+FSNCS*RCZS
        PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+FSNCS*RCPS
        SNCSH=SNCSH-FSNCS*RCCS
        SNCT=SNCT-FSNCS*RCCS
        IF(WTSHEBs1(NB,NZ).LE.ZEROLs1(NZ))THEN
          WTSHEBs1(NB,NZ)=0._r8
        ENDIF
!
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
        IF(SNCSH.LE.ZEROPs1(NZ))GO TO 565
  !
      !     OTHERWISE REMAINING C,N,P IN SHEATH OR PETIOLE GOES TO LITTERFALL
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FWODB=C woody fraction in branch:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
      !     HTSHE=petiole length
      !     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !
      ELSE
        DO 6325 M=1,4
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(5,M,NZ) &
            *WGSHEs1(K,NB,NZ)*FWODBs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(5,M,NZ) &
            *WGSHNs1(K,NB,NZ)*FWODSNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(5,M,NZ) &
            *WGSHPs1(K,NB,NZ)*FWODSPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(2,M,NZ) &
            *WGSHEs1(K,NB,NZ)*FWODBs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(2,M,NZ) &
            *WGSHNs1(K,NB,NZ)*FWODSNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(2,M,NZ) &
            *WGSHPs1(K,NB,NZ)*FWODSPs1(1)
6325    CONTINUE
        WTSHEBs1(NB,NZ)=AMAX1(0.0,WTSHEBs1(NB,NZ)-WGSHEs1(K,NB,NZ))
        WTSHBNs1(NB,NZ)=AMAX1(0.0,WTSHBNs1(NB,NZ)-WGSHNs1(K,NB,NZ))
        WTSHBPs1(NB,NZ)=AMAX1(0.0,WTSHBPs1(NB,NZ)-WGSHPs1(K,NB,NZ))
        HTSHEs1(K,NB,NZ)=0._r8
        WGSHEs1(K,NB,NZ)=0._r8
        WGSHNs1(K,NB,NZ)=0._r8
        WGSHPs1(K,NB,NZ)=0._r8
        WSSHEs1(K,NB,NZ)=0._r8
        IF(WTSHEBs1(NB,NZ).LE.ZEROLs1(NZ))THEN
          WTSHEBs1(NB,NZ)=0._r8
        ENDIF
      ENDIF
650 CONTINUE
    KN=KN+1
    SNCR=SNCT*(1.0-SNCF)
!
!     REMOBILIZATION OF RESERVE C
!
!     WTRSVB=stalk reserve C mass
!     SNCR=excess maintenance respiration
!
    IF(WTRSVBs1(NB,NZ).GT.SNCR)THEN
      WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)-SNCR
      SNCR=0._r8
      cycle
    ENDIF
!
!     REMOBILIZATION OF STALK C,N,P
!
!     FXFS=rate constant for remobilization of stalk C,N,P (h-1)
!     SNCZ=phenologically-driven respiration senescence during late-season
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     WTSTKB,WVSTKB=stalk,sapwood C mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     MXNOD,MNNOD=max,min node number currently growing
!     NNOD=number of concurrently growing nodes
!     KVSTG=integer of most recent leaf number
!
    SNCZ=FXFS*SNCR
    SNCT=SNCR+SNCZ
    IF(ISTYPs1(NZ).NE.0.AND.SNCT.GT.ZEROPs1(NZ) &
      .AND.WTSTKBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      SNCF=SNCZ/SNCT
      FRCC=WVSTKBs1(NB,NZ)/WTSTKBs1(NB,NZ)
      RCSC=RCCC*FRCC
      RCSN=RCCN*FRCC
      RCSP=RCCP*FRCC
      MXNOD=KVSTGs1(NB,NZ)
      MNNOD=MAX(MIN(0,MAX(0,MXNOD-NNODs1(NZ))),KVSTGs1(NB,NZ)-23)
      MXNOD=MAX(MXNOD,MNNOD)
      DO 1650 KK=MXNOD,MNNOD,-1
        K=MOD(KK,JNODS1)
        IF(K.EQ.0.AND.KK.NE.0)K=25
!     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
    !     WRITE(*,2356)'WGNODE1',I,J,NX,NY,NZ,NB,K,KK,MXNOD,MNNOD
    !    2,KSNC,RCCC,FRCC,RCSC,SNCT,WGNODEs1(K,NB,NZ)
    !    3,HTNODXs1(K,NB,NZ),WTSTKBs1(NB,NZ)
    !    4,CPOOLs1(NB,NZ)
    !     ENDIF
    !
    !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !
        IF(WGNODEs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
          RCCK=RCSC*WGNODEs1(K,NB,NZ)
          RCZK=WGNODNs1(K,NB,NZ)*(RCSN+(1.0-RCSN)*RCSC)
          RCPK=WGNODPs1(K,NB,NZ)*(RCSP+(1.0-RCSP)*RCSC)
    !
      !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
      !
      !     FSNCS=fraction of lowest internode to be remobilized
!
          IF(RCCK.GT.ZEROPs1(NZ))THEN
            FSNCK=AMAX1(0.0,AMIN1(1.0,SNCT/RCCK))
          ELSE
            FSNCK=1.0_r8
          ENDIF
    !
      !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
      !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
      !     FSNCK=fraction of lowest internode to be remobilized
      !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
      !     WGNODE,WGNODN,WGNODP=senescing internode C,N,P mass
!
          DO 7310 M=1,4
            CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
              *FSNCK*(WGNODEs1(K,NB,NZ)-RCCK)*FWOODs1(0)
            ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
              *FSNCK*(WGNODNs1(K,NB,NZ)-RCZK)*FWOODNs1(0)
            PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
              *FSNCK*(WGNODPs1(K,NB,NZ)-RCPK)*FWOODPs1(0)
            CSNCs1(M,1,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
              *FSNCK*(WGNODEs1(K,NB,NZ)-RCCK)*FWOODs1(1)
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
              *FSNCK*(WGNODNs1(K,NB,NZ)-RCZK)*FWOODNs1(1)
            PSNCs1(M,1,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
              *FSNCK*(WGNODPs1(K,NB,NZ)-RCPK)*FWOODPs1(1)
7310      CONTINUE
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     FSNCK=fraction of lowest internode to be remobilized
      !     HTNODE,HTNODX=living,senescing internode length
      !     WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
      !
          WTSTKBs1(NB,NZ)=AMAX1(0.0,WTSTKBs1(NB,NZ) &
            -FSNCK*WGNODEs1(K,NB,NZ))
          WTSTBNs1(NB,NZ)=AMAX1(0.0,WTSTBNs1(NB,NZ) &
            -FSNCK*WGNODNs1(K,NB,NZ))
          WTSTBPs1(NB,NZ)=AMAX1(0.0,WTSTBPs1(NB,NZ) &
            -FSNCK*WGNODPs1(K,NB,NZ))
          HTNODEs1(K,NB,NZ)=HTNODEs1(K,NB,NZ) &
            -FSNCK*HTNODXs1(K,NB,NZ)
          WGNODEs1(K,NB,NZ)=WGNODEs1(K,NB,NZ) &
            -FSNCK*WGNODEs1(K,NB,NZ)
          WGNODNs1(K,NB,NZ)=WGNODNs1(K,NB,NZ) &
            -FSNCK*WGNODNs1(K,NB,NZ)
          WGNODPs1(K,NB,NZ)=WGNODPs1(K,NB,NZ) &
            -FSNCK*WGNODPs1(K,NB,NZ)
          HTNODXs1(K,NB,NZ)=HTNODXs1(K,NB,NZ) &
            -FSNCK*HTNODXs1(K,NB,NZ)
!
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
      !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
      !     FSNCK=fraction of lowest internode to be remobilized
      !     SNCT=remaining node senescence respiration
      !
          WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+FSNCK*RCCK*SNCF
          WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+FSNCK*RCZK
          WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+FSNCK*RCPK
          SNCT=SNCT-FSNCK*RCCK
      !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
      !     WRITE(*,2356)'WGNODE2',I,J,NX,NY,NZ,NB,K,KK,MXNOD,MNNOD
      !    2,KSNC,RCCC,FRCC,RCSC,SNCT,WGNODEs1(K,NB,NZ)
      !    3,HTNODXs1(K,NB,NZ),WTSTKBs1(NB,NZ)
      !    4,CPOOLs1(NB,NZ)
      !2356  FORMAT(A8,11I4,12E16.8)
      !     ENDIF
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
          IF(SNCT.LE.ZEROPs1(NZ))GO TO 565
    !
    !       OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
    !
    !       CSNC,ZSNC,PSNC=literfall C,N,P
    !       CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !       WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
    !       HTNODE,HTNODX=living,senescing internode length
    !
        ELSE
          DO 7315 M=1,4
            CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
              *WGNODEs1(K,NB,NZ)*FWOODs1(0)
            ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
              *WGNODNs1(K,NB,NZ)*FWOODNs1(0)
            PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
              *WGNODPs1(K,NB,NZ)*FWOODPs1(0)
            CSNCs1(M,1,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
              *WGNODEs1(K,NB,NZ)*FWOODs1(1)
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
              *WGNODNs1(K,NB,NZ)*FWOODNs1(1)
            PSNCs1(M,1,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
              *WGNODPs1(K,NB,NZ)*FWOODPs1(1)
7315      CONTINUE
          WTSTKBs1(NB,NZ)=AMAX1(0.0,WTSTKBs1(NB,NZ) &
            -WGNODEs1(K,NB,NZ))
          WTSTBNs1(NB,NZ)=AMAX1(0.0,WTSTBNs1(NB,NZ) &
            -WGNODNs1(K,NB,NZ))
          WTSTBPs1(NB,NZ)=AMAX1(0.0,WTSTBPs1(NB,NZ) &
            -WGNODPs1(K,NB,NZ))
          HTNODEs1(K,NB,NZ)=HTNODEs1(K,NB,NZ) &
            -HTNODXs1(K,NB,NZ)
          WGNODEs1(K,NB,NZ)=0._r8
          WGNODNs1(K,NB,NZ)=0._r8
          WGNODPs1(K,NB,NZ)=0._r8
          HTNODXs1(K,NB,NZ)=0._r8
        ENDIF
1650  CONTINUE
!
    !   RESIDUAL STALK
    !
    !   RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
      IF(WTSTXBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
        RCCK=RCSC*WTSTXBs1(NB,NZ)
        RCZK=WTSTXNs1(NB,NZ)*(RCSN+(1.0-RCSN)*RCSC)
        RCPK=WTSTXPs1(NB,NZ)*(RCSP+(1.0-RCSP)*RCSC)
    !
    !     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
        IF(RCCK.GT.ZEROPs1(NZ))THEN
          FSNCR=AMAX1(0.0,AMIN1(1.0,SNCT/RCCK))
        ELSE
          FSNCR=1.0_r8
        ENDIF
    !
    !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !     TO FRACTIONS SET IN 'STARTQ'
    !
        DO 8310 M=1,4
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
            *FSNCR*(WTSTXBs1(NB,NZ)-RCCK)*FWOODs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
            *FSNCR*(WTSTXNs1(NB,NZ)-RCZK)*FWOODNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
            *FSNCR*(WTSTXPs1(NB,NZ)-RCPK)*FWOODPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
            *FSNCR*(WTSTXBs1(NB,NZ)-RCCK)*FWOODs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
            *FSNCR*(WTSTXNs1(NB,NZ)-RCZK)*FWOODNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
            *FSNCR*(WTSTXPs1(NB,NZ)-RCPK)*FWOODPs1(1)
8310    CONTINUE
    !
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !     HTNODE,HTNODX=living,senescing internode length
    !
        WTSTKBs1(NB,NZ)=AMAX1(0.0,WTSTKBs1(NB,NZ) &
          -FSNCR*WTSTXBs1(NB,NZ))
        WTSTBNs1(NB,NZ)=AMAX1(0.0,WTSTBNs1(NB,NZ) &
          -FSNCR*WTSTXNs1(NB,NZ))
        WTSTBPs1(NB,NZ)=AMAX1(0.0,WTSTBPs1(NB,NZ) &
          -FSNCR*WTSTXPs1(NB,NZ))
        WTSTXBs1(NB,NZ)=AMAX1(0.0,WTSTXBs1(NB,NZ) &
          -FSNCR*WTSTXBs1(NB,NZ))
        WTSTXNs1(NB,NZ)=AMAX1(0.0,WTSTXNs1(NB,NZ) &
          -FSNCR*WTSTXNs1(NB,NZ))
        WTSTXPs1(NB,NZ)=AMAX1(0.0,WTSTXPs1(NB,NZ) &
          -FSNCR*WTSTXPs1(NB,NZ))
        HTNODZ=0._r8
        DO 8320 K=0,JNODS1
          HTNODZ=AMAX1(HTNODZ,HTNODEs1(K,NB,NZ))
8320    CONTINUE
        HTNODZ=AMAX1(0.0,HTNODZ-FSNCR*HTNODZ)
        DO 8325 K=0,JNODS1
          HTNODEs1(K,NB,NZ)=AMIN1(HTNODZ,HTNODEs1(K,NB,NZ))
8325    CONTINUE
!
    !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
    !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
    !
    !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
    !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !     FSNCR=fraction of residual stalk to be remobilized
    !     SNCT=remaining node senescence respiration
    !
        WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+FSNCR*RCCK*SNCF
        WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+FSNCR*RCZK
        WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+FSNCR*RCPK
        SNCT=SNCT-FSNCR*RCCK
      ENDIF
  !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
  !     WRITE(*,2357)'WTSTXB1',I,J,NZ,NB,K,FSNCR,SNCT
  !    3,WTSTKBs1(NB,NZ),WTSTXBs1(NB,NZ)
  !    4,(HTNODEs1(K,NB,NZ),K=0,JNODS1)
  !2357  FORMAT(A8,5I4,40E12.4)
  !     ENDIF
  !
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
      IF(SNCT.LE.ZEROPs1(NZ))cycle
  !
  !     OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
  !
  !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
  !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
  !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
  !
    ELSE
      DO 8315 M=1,4
        CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
          *WTSTXBs1(NB,NZ)*FWOODs1(0)
        ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
          *WTSTXNs1(NB,NZ)*FWOODNs1(0)
        PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
          *WTSTXPs1(NB,NZ)*FWOODPs1(0)
        CSNCs1(M,1,0,NZ)=CSNCs1(M,0,0,NZ)+CFOPCs1(3,M,NZ) &
          *WTSTXBs1(NB,NZ)*FWOODs1(1)
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ) &
          *WTSTXNs1(NB,NZ)*FWOODNs1(1)
        PSNCs1(M,1,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ) &
          *WTSTXPs1(NB,NZ)*FWOODPs1(1)
8315  CONTINUE
      WTSTKBs1(NB,NZ)=AMAX1(0.0,WTSTKBs1(NB,NZ)-WTSTXBs1(NB,NZ))
      WTSTBNs1(NB,NZ)=AMAX1(0.0,WTSTBNs1(NB,NZ)-WTSTXNs1(NB,NZ))
      WTSTBPs1(NB,NZ)=AMAX1(0.0,WTSTBPs1(NB,NZ)-WTSTXPs1(NB,NZ))
      WTSTXBs1(NB,NZ)=0._r8
      WTSTXNs1(NB,NZ)=0._r8
      WTSTXPs1(NB,NZ)=0._r8
    ENDIF
!
!     REMOBILIZATION OF STORAGE C,N,P
!
!     WTRVC=storage C
!     IDTHB=branch living flag: 0=alive,1=dead
!     SNCR=remaining excess maintenance respiration
!
    SNCR=SNCT/(1.0+FXFS)
    IF(WTRVCs1(NZ).GT.SNCR)THEN
      WTRVCs1(NZ)=WTRVCs1(NZ)-SNCR
      SNCR=0._r8
    ELSEIF(ISTYPs1(NZ).NE.0)THEN
      IDTHBs1(NB,NZ)=1
    ENDIF
565 CONTINUE
575 CONTINUE

  end subroutine RemobilizeLeafLayers

!------------------------------------------------------------------------------------------

  subroutine PhenologyReset(I,NB,NZ)

  implicit none
  integer, intent(in) :: I,nb,nz
  integer :: K,M
! begin_execution
!   RESET PHENOLOGY AT EMERGENCE ('VRNS' > 'VRNL')
!   AND END OF SEASON ('VRNF' > 'VRNX')
!
!   ISTYP=growth habit:0=annual,1=perennial from PFT file
!   IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!   IFLGE=flag for enabling leafout:0=enable,1=disable
!   VRNS,VRNL=leafout hours,hours required for leafout
!   IFLGF=flag for enabling leafoff:0=enable,1=disable
!   VRNF,VRNX=leafoff hours,hours required for leafoff
!
  IF(ISTYPs1(NZ).NE.0.OR.(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0))THEN
    IF((IFLGEs1(NB,NZ).EQ.0.AND.VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ)) &
      .OR.(IFLGFs1(NB,NZ).EQ.0.AND.VRNFs1(NB,NZ).GE.VRNXs1(NB,NZ)))THEN
!
  !    SPRING PHENOLOGY RESET
  !
  !    GROUP,GROUPI=node number required for floral initiation
  !    PSTGI=node number at floral initiation
  !    PSTGF=node number at flowering
  !    VSTGX=leaf number on date of floral initiation
  !    TGSTGI=total change in vegve node number normalized for maturity group
  !    TGSTGF=total change in reprve node number normalized for maturity group
  !    IDAYs1(1,=emergence date
!
      IF((IFLGEs1(NB,NZ).EQ.0.AND.ISTYPs1(NZ).NE.0) &
        .AND.(VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ)))THEN
        IF(ISTYPs1(NZ).EQ.0)THEN
          GROUPs1(NB,NZ)=AMAX1(0.0,GROUPIs1(NZ)-NBTBs1(NB,NZ))
        ELSE
          GROUPs1(NB,NZ)=GROUPIs1(NZ)
        ENDIF
        PSTGIs1(NB,NZ)=PSTGs1(NB,NZ)
        PSTGFs1(NB,NZ)=0._r8
        VSTGXs1(NB,NZ)=0._r8
        TGSTGIs1(NB,NZ)=0._r8
        TGSTGFs1(NB,NZ)=0._r8
        IDAYs1(1,NB,NZ)=I
        DO 2005 M=2,10
          IDAYs1(M,NB,NZ)=0
2005    CONTINUE
        IF(NB.EQ.NB1s1(NZ))THEN
          WSTRs1(NZ)=0._r8
        ENDIF
    !
    !   SPRING LEAF AND SHEATH RESET
    !
    !   IFLGA,IFLGE=flag for initializing,enabling leafout
    !   VRNS,VRNL=leafout hours,hours required for leafout
    !   PSTG=node number
    !   VSTG=number of leaves appeared
    !     KVSTG=integer of most recent leaf number currently growing
    !     FLG4=number of hours with no grain fill
    !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !     WT*,WT*N,WT*P=branch organ C,N,P mass
    !     WG,WG*N,WG*P=node organ C,N,P mass
    !     organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
    !     HSK=husk,EAR=ear,GR=grain,SHT=shoot
    !     ARLFB,ARLF=branch,node leaf area
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !     HTNODX,HTNODE=stalk height,stalk internode length
    !     GRNOB=seed set number
    !     GRNXB=potential number of seed set sites
    !     GRWTB=individual seed size
!
        IF(IFLGEs1(NB,NZ).EQ.0.AND.ISTYPs1(NZ).NE.0 &
          .AND.VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ))THEN
          IF(IBTYPs1(NZ).EQ.0)THEN
            PSTGs1(NB,NZ)=XTLIs1(NZ)
            VSTGs1(NB,NZ)=0._r8
            KLEAFs1(NB,NZ)=1
            KVSTGs1(NB,NZ)=1
            FLG4s1(NB,NZ)=0._r8
            DO 5330 M=1,4
              CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ) &
                +CFOPCs1(5,M,NZ)*WTLFBs1(NB,NZ)*FWODBs1(0) &
                +CFOPCs1(5,M,NZ)*WTSHEBs1(NB,NZ)*FWODBs1(0)
              ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ) &
                +CFOPNs1(5,M,NZ)*WTLFBNs1(NB,NZ)*FWODLNs1(0) &
                +CFOPNs1(5,M,NZ)*WTSHBNs1(NB,NZ)*FWODSNs1(0)
              PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ) &
                +CFOPPs1(5,M,NZ)*WTLFBPs1(NB,NZ)*FWODLPs1(0) &
                +CFOPPs1(5,M,NZ)*WTSHBPs1(NB,NZ)*FWODSPs1(0)
              CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
                +CFOPCs1(1,M,NZ)*WTLFBs1(NB,NZ)*FWODBs1(1) &
                +CFOPCs1(2,M,NZ)*WTSHEBs1(NB,NZ)*FWODBs1(1)
              ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
                +CFOPNs1(1,M,NZ)*WTLFBNs1(NB,NZ)*FWODLNs1(1) &
                +CFOPNs1(2,M,NZ)*WTSHBNs1(NB,NZ)*FWODSNs1(1)
              PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
                +CFOPPs1(1,M,NZ)*WTLFBPs1(NB,NZ)*FWODLPs1(1) &
                +CFOPPs1(2,M,NZ)*WTSHBPs1(NB,NZ)*FWODSPs1(1)
5330        CONTINUE
            ARLFBs1(NB,NZ)=0._r8
            WTLFBs1(NB,NZ)=0._r8
            WTLFBNs1(NB,NZ)=0._r8
            WTLFBPs1(NB,NZ)=0._r8
            WTSHEBs1(NB,NZ)=0._r8
            WTSHBNs1(NB,NZ)=0._r8
            WTSHBPs1(NB,NZ)=0._r8
            DO 5335 K=0,JNODS1
              ARLF1s1(K,NB,NZ)=0._r8
              HTSHEs1(K,NB,NZ)=0._r8
              WGLFs1(K,NB,NZ)=0._r8
              WSLFs1(K,NB,NZ)=0._r8
              WGLFNs1(K,NB,NZ)=0._r8
              WGLFPs1(K,NB,NZ)=0._r8
              WGSHEs1(K,NB,NZ)=0._r8
              WSSHEs1(K,NB,NZ)=0._r8
              WGSHNs1(K,NB,NZ)=0._r8
              WGSHPs1(K,NB,NZ)=0._r8
5335        CONTINUE
          ENDIF
        ENDIF
    !
    !     RESIDUAL STALKS BECOME LITTERFALL IN GRASSES, SHRUBS AT
    !     START OF SEASON
    !
        IF((IFLGEs1(NB,NZ).EQ.0.AND.ISTYPs1(NZ).NE.0) &
          .AND.VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ))THEN
          DO 6245 M=1,4
            CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(2,M,NZ) &
              *(WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ)+WTGRBs1(NB,NZ))
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(2,M,NZ) &
              *(WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ)+WTGRBNs1(NB,NZ))
            PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(2,M,NZ) &
              *(WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ)+WTGRBPs1(NB,NZ))
6245      CONTINUE
          WTHSKBs1(NB,NZ)=0._r8
          WTEARBs1(NB,NZ)=0._r8
          WTGRBs1(NB,NZ)=0._r8
          WTHSBNs1(NB,NZ)=0._r8
          WTEABNs1(NB,NZ)=0._r8
          WTGRBNs1(NB,NZ)=0._r8
          WTHSBPs1(NB,NZ)=0._r8
          WTEABPs1(NB,NZ)=0._r8
          WTGRBPs1(NB,NZ)=0._r8
          GRNXBs1(NB,NZ)=0._r8
          GRNOBs1(NB,NZ)=0._r8
          GRWTBs1(NB,NZ)=0._r8
          IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
            DO 6345 M=1,4
              CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
              ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
              PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
6345        CONTINUE
            WTSTKBs1(NB,NZ)=0._r8
            WTSTBNs1(NB,NZ)=0._r8
            WTSTBPs1(NB,NZ)=0._r8
            WTSTXBs1(NB,NZ)=0._r8
            WTSTXNs1(NB,NZ)=0._r8
            WTSTXPs1(NB,NZ)=0._r8
            DO 6340 K=0,JNODS1
              HTNODEs1(K,NB,NZ)=0._r8
              HTNODXs1(K,NB,NZ)=0._r8
              WGNODEs1(K,NB,NZ)=0._r8
              WGNODNs1(K,NB,NZ)=0._r8
              WGNODPs1(K,NB,NZ)=0._r8
6340        CONTINUE
          ENDIF
        ENDIF
      ENDIF
!
  !   SPRING OR FALL FLAG RESET
  !
      IF(IFLGEs1(NB,NZ).EQ.0 &
        .AND.VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ))THEN
        IFLGEs1(NB,NZ)=1
        IFLGFs1(NB,NZ)=0
        IFLGRs1(NB,NZ)=0
        IFLGQs1(NB,NZ)=0
      ELSE
        IFLGEs1(NB,NZ)=0
        IFLGFs1(NB,NZ)=1
        IFLGRs1(NB,NZ)=1
        IFLGQs1(NB,NZ)=0
        IFLGAs1(NB,NZ)=0
      ENDIF
    ENDIF
  ENDIF
!
!   REPRODUCTIVE MATERIAL BECOMES LITTERFALL AT END OF SEASON
!
  IF(IFLGRs1(NB,NZ).EQ.1)THEN
    IFLGQs1(NB,NZ)=IFLGQs1(NB,NZ)+1
    IF(IFLGQs1(NB,NZ).EQ.IFLGQX)THEN
      IFLGRs1(NB,NZ)=0
      IFLGQs1(NB,NZ)=0
    ENDIF
    DO 6330 M=1,4
      CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+FSNR*CFOPCs1(2,M,NZ) &
        *(WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ))
      ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+FSNR*CFOPNs1(2,M,NZ) &
        *(WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ))
      PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+FSNR*CFOPPs1(2,M,NZ) &
        *(WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ))
      IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
        WTRVCs1(NZ)=WTRVCs1(NZ) &
          +FSNR*CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
        WTRVNs1(NZ)=WTRVNs1(NZ) &
          +FSNR*CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
        WTRVPs1(NZ)=WTRVPs1(NZ) &
          +FSNR*CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
      ELSE
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
          +FSNR*CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
          +FSNR*CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
          +FSNR*CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
      ENDIF
6330  CONTINUE
    WTHSKBs1(NB,NZ)=(1.0-FSNR)*WTHSKBs1(NB,NZ)
    WTEARBs1(NB,NZ)=(1.0-FSNR)*WTEARBs1(NB,NZ)
    WTGRBs1(NB,NZ)=(1.0-FSNR)*WTGRBs1(NB,NZ)
    WTHSBNs1(NB,NZ)=(1.0-FSNR)*WTHSBNs1(NB,NZ)
    WTEABNs1(NB,NZ)=(1.0-FSNR)*WTEABNs1(NB,NZ)
    WTGRBNs1(NB,NZ)=(1.0-FSNR)*WTGRBNs1(NB,NZ)
    WTHSBPs1(NB,NZ)=(1.0-FSNR)*WTHSBPs1(NB,NZ)
    WTEABPs1(NB,NZ)=(1.0-FSNR)*WTEABPs1(NB,NZ)
    WTGRBPs1(NB,NZ)=(1.0-FSNR)*WTGRBPs1(NB,NZ)
    GRNXBs1(NB,NZ)=(1.0-FSNR)*GRNXBs1(NB,NZ)
    GRNOBs1(NB,NZ)=(1.0-FSNR)*GRNOBs1(NB,NZ)
    GRWTBs1(NB,NZ)=(1.0-FSNR)*GRWTBs1(NB,NZ)
!
!     STALKS BECOME LITTERFALL IN GRASSES AT END OF SEASON
!
    IF((IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1).AND.ISTYPs1(NZ).NE.0)THEN
      DO 6335 M=1,4
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+FSNR*CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+FSNR*CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+FSNR*CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
6335    CONTINUE
      WTSTKBs1(NB,NZ)=(1.0-FSNR)*WTSTKBs1(NB,NZ)
      WTSTBNs1(NB,NZ)=(1.0-FSNR)*WTSTBNs1(NB,NZ)
      WTSTBPs1(NB,NZ)=(1.0-FSNR)*WTSTBPs1(NB,NZ)
      WTSTXBs1(NB,NZ)=(1.0-FSNR)*WTSTXBs1(NB,NZ)
      WTSTXNs1(NB,NZ)=(1.0-FSNR)*WTSTXNs1(NB,NZ)
      WTSTXPs1(NB,NZ)=(1.0-FSNR)*WTSTXPs1(NB,NZ)
      DO 2010 K=0,JNODS1
    !     HTNODEs1(K,NB,NZ)=(1.0-FSNR)*HTNODEs1(K,NB,NZ)
        HTNODXs1(K,NB,NZ)=(1.0-FSNR)*HTNODXs1(K,NB,NZ)
        WGNODEs1(K,NB,NZ)=(1.0-FSNR)*WGNODEs1(K,NB,NZ)
        WGNODNs1(K,NB,NZ)=(1.0-FSNR)*WGNODNs1(K,NB,NZ)
        WGNODPs1(K,NB,NZ)=(1.0-FSNR)*WGNODPs1(K,NB,NZ)
2010  CONTINUE
    ENDIF
!
!     SELF-SEEDING ANNUALS IF COLD OR DROUGHT DECIDUOUS
!
!     ISTYP=growth habit:0=annual,1=perennial
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IDAYH,IYRH=day,year of harvesting
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     IDAY0,IYR0=day,year of planting
!     IFLGI=PFT initialization flag:0=no,1=yes
!
!     IF(J.EQ.INT(ZNOONs1))THEN
    IF(NB.EQ.NB1s1(NZ))THEN
      IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
        IDAYHs1(NZ)=I
        IYRHs1(NZ)=IYRCs1
        IHVSTs1(NZ)=1
        JHVSTs1(NZ)=2
        HVSTs1(NZ)=0._r8
        THINs1(NZ)=0._r8
        EHVSTs1(1,1,NZ)=1.0_r8
        EHVSTs1(1,2,NZ)=1.0_r8
        EHVSTs1(1,3,NZ)=1.0_r8
        EHVSTs1(1,4,NZ)=1.0_r8
        EHVSTs1(2,1,NZ)=0._r8
        EHVSTs1(2,2,NZ)=1.0_r8
        EHVSTs1(2,3,NZ)=0._r8
        EHVSTs1(2,4,NZ)=0._r8
        IDAY0s1(NZ)=-1E+06
        IYR0s1(NZ)=-1E+06
        IFLGIs1(NZ)=1
!     WRITE(*,3366)'HVST',I,J,IYRCs1,IDAYHs1(NZ),IYRHs1(NZ)
!    2,IHVSTs1(NZ),JHVSTs1(NZ),IFLGIs1(NZ)
!3366  FORMAT(A8,8I8)
      ENDIF
    ENDIF
!     ENDIF
  ENDIF
  end subroutine PhenologyReset

!------------------------------------------------------------------------------------------

  subroutine AllocateLeafToCanopyLayers(NB,NZ,ZCX)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: ZCX(JP1)
  integer  :: LL,LU,L,K,k1,k2,KK
  integer  :: KVSTGX,KVSTG1,LHTLFU,LHTLFL
  integer  :: LHTBRU,LHTBRL,N
  real(r8) :: ZSTK
  real(r8) :: YWGLFN,YWGLFP
  real(r8) :: YARLF,YWGLF,YLGLF,XLGLF
  real(r8) :: ARSTKB,ASTV
  real(r8) :: FRACL
  real(r8) :: HTBR
  real(r8) :: HTSTK
  real(r8) :: HTLF,HTLFL,HTLFU
  real(r8) :: HTLFB
  real(r8) :: RSTK
  real(r8) :: TLGLF
! begin_execution
!   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
!
!   HTCTL=hypocotyledon height
!   SDPTH=seeding depth
!   ARLF=node leaf area
!   HTSHE=petiole length
!
  KVSTGNs1(NB,NZ)=0
  IF(HTCTLs1(NZ).LE.SDPTHs1(NZ) &
    .AND.ARLF1s1(0,NB1s1(NZ),NZ).GT.0.0)THEN
    XLGLF=SQRT(1.0E+02*ARLF1s1(0,NB1s1(NZ),NZ)/PPs1(NZ))
    HTCTLs1(NZ)=XLGLF+HTSHEs1(0,NB1s1(NZ),NZ) &
      +HTNODEs1(0,NB1s1(NZ),NZ)
  ENDIF
!
! IF CANOPY HAS EMERGED
!
  IF(HTCTLs1(NZ).GT.SDPTHs1(NZ))THEN
    DO 540 K=0,JNODS1
      DO  L=1,JC1
        ARLFLs1(L,K,NB,NZ)=0._r8
        WGLFLs1(L,K,NB,NZ)=0._r8
        WGLFLNs1(L,K,NB,NZ)=0._r8
        WGLFLPs1(L,K,NB,NZ)=0._r8
      enddo
540   CONTINUE
    DO 535 L=1,JC1
      ARSTKs1(L,NB,NZ)=0._r8
535 CONTINUE
!
!   BRANCH HEIGHT
!
!   IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   KVSTG1,KVSTGN=integer of highest,lowest leaf number currently growing
!   HTNODE=internode length
!   HTBR=branch base height
!
    IF(IBTYPs1(NZ).NE.0.AND.IGTYPs1(NZ).GT.1)THEN
      IF(NB.NE.NB1s1(NZ))THEN
        KVSTG1=MAX(1,KVSTGs1(NB1s1(NZ),NZ)-24)
        IF(NBTBs1(NB,NZ).GE.KVSTG1)THEN
          K=MOD(NBTBs1(NB,NZ),JNODS1)
!        IF(K.EQ.0.AND.KK.NE.0)K=JNODS1
          IF(K.EQ.0)K=JNODS1
          HTBR=HTNODEs1(K,NB1s1(NZ),NZ)
        ELSE
          HTBR=0._r8
        ENDIF
      ELSE
        HTBR=0._r8
      ENDIF
    ELSE
      HTBR=0._r8
    ENDIF
    KVSTGX=MAX(0,KVSTGs1(NB,NZ)-24)
!
!   FOR ALL LEAFED NODES
!
    DO 560 KK=KVSTGX,KVSTGs1(NB,NZ)
      K=MOD(KK,JNODS1)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      !
      !     HEIGHT OF STALK INTERNODE + SHEATH OR PETIOLE
      !     AND LENGTH OF LEAF
      !
      !     HTSTK=stalk height
      !     HTNODE=internode length
      !     HTLF=leaf node height
      !     ARLF=leaf node area
      !     PP=plant population
      !     FNOD=scales node number for perennial vegetation (e.g. trees)
      !     XLGLF=leaf length
!
      HTSTK=HTBR+HTNODEs1(K,NB,NZ)
      HTLF=HTSTK+HTSHEs1(K,NB,NZ)
      XLGLF=AMAX1(0.0,SQRT(WDLFs1(NZ)*AMAX1(0.0 &
        ,ARLF1s1(K,NB,NZ))/(PPs1(NZ)*FNODs1(NZ))))
      TLGLF=0._r8
      !
      !   ALLOCATE FRACTIONS OF LEAF IN EACH INCLINATION CLASS
      !     FROM HIGHEST TO LOWEST TO CANOPY LAYER
      !
      !     YLGLF=leaf elevation
      !     CLASS=leaf inclination class
      !     XLGLF=leaf length
      !     ZC,ZCX=canopy height
      !     HTLFL,HTLFU=height of leaf base,tip
      !     ZL=height to bottom of each canopy layer
      !     LHTLFL,LHTLFU=layer number of leaf base,tip
      !     FRACL=leaf fraction in each layer
!
      DO 555 N=JLI1,1,-1
        YLGLF=ZSINs1(N)*CLASSs1(N,NZ)*XLGLF
        HTLFL=AMIN1(ZCX(NZ)+0.01-YLGLF,HTLF+TLGLF)
        HTLFU=AMIN1(ZCX(NZ)+0.01,HTLFL+YLGLF)
        LU=0
        LL=0
        DO 550 L=JC1,1,-1
          IF(LU.EQ.1.AND.LL.EQ.1)exit
          IF((HTLFU.GT.ZLs1(L-1).OR.ZLs1(L-1).LE.ZEROs1).AND.LU.EQ.0)THEN
            LHTLFU=MAX(1,L)
            LU=1
          ENDIF
          IF((HTLFL.GT.ZLs1(L-1).OR.ZLs1(L-1).LE.ZEROs1).AND.LL.EQ.0)THEN
            LHTLFL=MAX(1,L)
            LL=1
          ENDIF
550     CONTINUE
        DO 570 L=LHTLFL,LHTLFU
          IF(LHTLFU.EQ.LHTLFL)THEN
            FRACL=CLASSs1(N,NZ)
          ELSEIF(HTLFU.GT.HTLFL.AND.ZLs1(L).GT.HTLFL)THEN
            FRACL=CLASSs1(N,NZ)*(AMIN1(HTLFU,ZLs1(L)) &
              -AMAX1(HTLFL,ZLs1(L-1)))/(HTLFU-HTLFL)
          ELSE
            FRACL=CLASSs1(N,NZ)
          ENDIF
          YARLF=FRACL*ARLF1s1(K,NB,NZ)
          YWGLF=FRACL*WGLFs1(K,NB,NZ)
          YWGLFN=FRACL*WGLFNs1(K,NB,NZ)
          YWGLFP=FRACL*WGLFPs1(K,NB,NZ)
!
    !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
    !
    !     ARLFL=leaf node area in canopy layer
    !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
    !     ARLFV,WGLFV=total leaf area,C in canopy layer
    !     HTNODE=internode length
    !
          ARLFLs1(L,K,NB,NZ)=ARLFLs1(L,K,NB,NZ)+YARLF
          WGLFLs1(L,K,NB,NZ)=WGLFLs1(L,K,NB,NZ)+YWGLF
          WGLFLNs1(L,K,NB,NZ)=WGLFLNs1(L,K,NB,NZ)+YWGLFN
          WGLFLPs1(L,K,NB,NZ)=WGLFLPs1(L,K,NB,NZ)+YWGLFP
          ARLFVs1(L,NZ)=ARLFVs1(L,NZ)+YARLF
          WGLFVs1(L,NZ)=WGLFVs1(L,NZ)+YWGLF
          !     IF(NZ.EQ.2)THEN
          !     WRITE(*,4813)'GRO',I,J,NX,NY,NZ,NB,K,KK,L,LHTLFL,LHTLFU
          !    2,FRACL,HTLFU,HTLFL,ZLs1(L-1),ARLFBs1(NB,NZ)
          !    3,ARLF1s1(K,NB,NZ),WTLFBs1(NB,NZ),WGLFs1(K,NB,NZ)
          !    4,ARLFPs1(NZ),ZLs1(L),HTLF,TLGLF,HTSTK,HTBR
          !    4,HTNODEs1(K,NB,NZ),HTSHEs1(K,NB,NZ),YLGLF
          !    5,ZSINs1(N),CLASSs1(N,NZ),XLGLF,ZCs1(NZ)
          !    6,ZCX(NZ)
          !4813  FORMAT(A8,11I4,30E12.4)
          !     ENDIF
570     CONTINUE
        TLGLF=TLGLF+YLGLF
        ZCs1(NZ)=AMAX1(ZCs1(NZ),HTLFU)
555   CONTINUE
      IF(WSSHEs1(K,NB,NZ).GT.0.0)THEN
        IF(KVSTGNs1(NB,NZ).EQ.0)KVSTGNs1(NB,NZ)=MIN(KK,KVSTGs1(NB,NZ))
      ENDIF
560 CONTINUE
    IF(KVSTGNs1(NB,NZ).EQ.0)KVSTGNs1(NB,NZ)=KVSTGs1(NB,NZ)
    K1=MOD(KVSTGs1(NB,NZ),JNODS1)
    IF(K1.EQ.0.AND.KVSTGs1(NB,NZ).NE.0)K1=JNODS1
    K2=MOD(KVSTGs1(NB,NZ)-1,JNODS1)
    IF(K2.EQ.0.AND.KVSTGs1(NB,NZ)-1.NE.0)K2=JNODS1
    IF(test_aeqb(HTNODEs1(K1,NB,NZ),0._r8))THEN
      HTNODEs1(K1,NB,NZ)=HTNODEs1(K2,NB,NZ)
    ENDIF
    HTLFB=HTBR+AMAX1(0.0,HTNODEs1(K1,NB,NZ))
!
  !     ALLOCATE STALK SURFACE AREA TO CANOPY LAYERS
  !
  !     HTNODE=internode length
  !     HTLFB=leaf base height
  !     ZL=height to bottom of each canopy layer
  !     LHTBRL,LHTBRU=layer number of branch base,tip
  !     WTSTKB,ARSTKB=branch stalk mass,surface area
  !     FSTK=fraction of stalk area contributing to water,heat flow
  !     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
  !     WVSTKB=stalk sapwood mass
  !     FRACL=stalk fraction in each layer
  !     ARSTK=total branch stalk surface area in each layer
  !
  !     IF(NZ.EQ.1)THEN
  !     WRITE(*,6679)'K1',I,J,NZ,NB,K1,KVSTGs1(NB,NZ)
  !    2,HTNODEs1(K1,NB,NZ)
  !6679  FORMAT(A8,6I4,12E12.4)
  !     ENDIF
    IF(HTNODEs1(K1,NB,NZ).GT.0.0)THEN
      LU=0
      LL=0
      DO 545 L=JC1,1,-1
        IF(LU.EQ.1.AND.LL.EQ.1)exit
        IF((HTLFB.GT.ZLs1(L-1).OR.ZLs1(L-1).LE.ZEROs1).AND.LU.EQ.0)THEN
          LHTBRU=MAX(1,L)
          LU=1
        ENDIF
        IF((HTBR.GT.ZLs1(L-1).OR.ZLs1(L-1).LE.ZEROs1).AND.LL.EQ.0)THEN
          LHTBRL=MAX(1,L)
          LL=1
        ENDIF
545   CONTINUE
      RSTK=SQRT(VSTK*(AMAX1(0.0_r8,WTSTKBs1(NB,NZ))/PPs1(NZ)) &
        /(PICON*HTNODEs1(K1,NB,NZ)))
      ARSTKB=PICON*HTNODEs1(K1,NB,NZ)*PPs1(NZ)*RSTK
      IF(ISTYPs1(NZ).EQ.0)THEN
        WVSTKBs1(NB,NZ)=WTSTKBs1(NB,NZ)
      ELSE
        ZSTK=AMIN1(ZSTX,FSTK*RSTK)
        ASTV=PICON*(2.0_r8*RSTK*ZSTK-ZSTK**2)
        WVSTKBs1(NB,NZ)=ASTV/VSTK*HTNODEs1(K1,NB,NZ)*PPs1(NZ)
      ENDIF
    !     IF(NZ.EQ.1)THEN
        !     WRITE(*,6677)'WVSTK',I,J,NX,NY,NZ,NB,K1,WVSTKBs1(NB,NZ)
        !    2,ASTV,VSTK,HTNODEs1(K1,NB,NZ),PPs1(NZ)
        !6677  FORMAT(A8,7I4,12E12.4)
        !     ENDIF
      DO 445 L=LHTBRL,LHTBRU
        IF(HTLFB.GT.HTBR)THEN
          IF(HTLFB.GT.ZLs1(L-1))THEN
            FRACL=(AMIN1(HTLFB,ZLs1(L))-AMAX1(HTBR &
              ,ZLs1(L-1)))/(HTLFB-HTBR)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        ARSTKs1(L,NB,NZ)=FRACL*ARSTKB
445   CONTINUE
    ELSE
      WVSTKBs1(NB,NZ)=0._r8
      DO 450 L=1,JC1
        ARSTKs1(L,NB,NZ)=0._r8
450   CONTINUE
    ENDIF
  ELSE
    WVSTKBs1(NB,NZ)=0._r8
    DO 455 L=1,JC1
      ARSTKs1(L,NB,NZ)=0._r8
455 CONTINUE
  ENDIF
end subroutine AllocateLeafToCanopyLayers

!------------------------------------------------------------------------------------------

  subroutine GrainFilling(I,NB,NZ,GROGR,GROSTK,GROGRN,GROGRP)
  implicit none
  integer, intent(in) :: I,NB,NZ
  real(r8), intent(in) :: GROGR,GROSTK,GROGRN,GROGRP
  real(r8) :: ZPGRP,ZPGRN,ZNPGP,ZNPGN
  real(r8) :: XLOCM,XLOCC,XLOCN,XLOCP
  real(r8) :: FGRNX
  real(r8) :: GRMXB
  real(r8) :: GROLM
  real(r8) :: GROLC
  real(r8) :: SET
! begin_execution

!
!   SET MAXIMUM GRAIN NUMBER FROM SHOOT MASS BEFORE ANTHESIS
!
!   IDAYs1(3,=start of stem elongation and setting max seed number
!   IDAYs1(6,=start of anthesis and setting final seed number
!   GRNXB=potential number of seed set sites
!   STMX=maximum potential seed number from PFT file
!     GROSTK=stalk growth rate
!
  IF(IDAYs1(3,NB,NZ).NE.0.AND.IDAYs1(6,NB,NZ).EQ.0)THEN
    GRNXBs1(NB,NZ)=GRNXBs1(NB,NZ)+STMXs1(NZ)*AMAX1(0.0,GROSTK)
!     WRITE(*,4246)'GRNX',I,J,NZ,NB,IDAYs1(3,NB,NZ)
!    2,GRNXBs1(NB,NZ),STMXs1(NZ),CGROS,GROSTK
  ENDIF
!
!   SET FINAL GRAIN NUMBER FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!   IDAYs1(6,=start of anthesis and setting final seed number
!   IDAYs1(7,=start of grain filling and setting max seed size
!   IDAYs1(8,=end date setting for final seed number
!   IDAYs1(9,=end of setting max seed size
!   SET=seed set limited by nonstructural C,N,P
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   TCC=canopy temperature
!   CTC=chilling temperature for CO2 fixation, seed loss (oC)
!   HTC=high temperature threshold for grain number loss
!   FGRNX=loss of seed set
!   SSTX=sensitivity to TCC > HTC,TCC < CTC from startq.f (seeds oC-1)
!   GRNOB=seed set number
!   GRNXB=potential number of seed set sites
!   SDMX=maximum seed number per STMX from PFT file
!   DGSTGF=change in reproductive node number normalized for maturity group
!
  IF(IDAYs1(6,NB,NZ).NE.0.AND.IDAYs1(9,NB,NZ).EQ.0)THEN
    SET=AMIN1(CCPOLBs1(NB,NZ)/(CCPOLBs1(NB,NZ)+SETC) &
      ,CZPOLBs1(NB,NZ)/(CZPOLBs1(NB,NZ)+SETN) &
      ,CPPOLBs1(NB,NZ)/(CPPOLBs1(NB,NZ)+SETP))
    IF(TCCs1(NZ).LT.CTCs1(NZ))THEN
      IF(IDAYs1(7,NB,NZ).EQ.0)THEN
        FGRNX=SSTXs1(NZ)*(CTCs1(NZ)-TCCs1(NZ))
      ELSEIF(IDAYs1(8,NB,NZ).EQ.0)THEN
        FGRNX=SSTXs1(NZ)*(CTCs1(NZ)-TCCs1(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSEIF(TCCs1(NZ).GT.HTCs1(NZ))THEN
      IF(IDAYs1(7,NB,NZ).EQ.0)THEN
        FGRNX=SSTXs1(NZ)*(TCCs1(NZ)-HTCs1(NZ))
      ELSEIF(IDAYs1(8,NB,NZ).EQ.0)THEN
        FGRNX=SSTXs1(NZ)*(TCCs1(NZ)-HTCs1(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSE
      FGRNX=0._r8
    ENDIF
    IF(IDAYs1(6,NB,NZ).NE.0.AND.IDAYs1(8,NB,NZ).EQ.0)THEN
!     GRNXBs1(NB,NZ)=GRNXBs1(NB,NZ)*FGRNX
      GRNOBs1(NB,NZ)=AMIN1(SDMXs1(NZ)*GRNXBs1(NB,NZ) &
        ,GRNOBs1(NB,NZ)+SDMXs1(NZ)*GRNXBs1(NB,NZ) &
        *SET*DGSTGFs1(NB,NZ)-FGRNX*GRNOBs1(NB,NZ))
!     IF(FGRNX.LT.1.0)THEN
!     WRITE(*,4246)'GRNO',I,J,NZ,NB,IDAYs1(7,NB,NZ),TCCs1(NZ)
!    2,HTCs1(NZ),FGRNX,GRNXBs1(NB,NZ),GRNOBs1(NB,NZ)
!    3,SET,CCPOLBs1(NB,NZ),CZPOLBs1(NB,NZ)
!    4,CPPOLBs1(NB,NZ)
!4246  FORMAT(A8,5I4,20E12.4)
!     ENDIF
    ENDIF
!

!     SET MAXIMUM GRAIN SIZE FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!     GRMX=maximum individual seed size from PFT file (g)
!     DGSTGF=change in reproductive node number normalized for maturity group
!     GRWTB=individual seed size
!     SET=seed set limited by nonstructural C,N,P
!
    IF(IDAYs1(7,NB,NZ).NE.0.AND.IDAYs1(9,NB,NZ).EQ.0)THEN
      GRMXB=GRMXs1(NZ)
      GRWTBs1(NB,NZ)=AMIN1(GRMXs1(NZ),GRWTBs1(NB,NZ) &
        +GRMXB*AMAX1(0.50,SET**0.25)*DGSTGFs1(NB,NZ))
!       IF(FGRNX.LT.1.0)THEN
!       WRITE(*,4246)'GRWT',I,J,NZ,NB,IDAYs1(8,NB,NZ),TCCs1(NZ)
!      2,HTCs1(NZ),FGRNX,GRMXs1(NZ),GRWTBs1(NB,NZ)
!       ENDIF
    ENDIF
  ENDIF
!
!   GRAIN FILL BY TRANSLOCATION FROM STALK RESERVES
!   UNTIL GRAIN SINK (=FINAL GRAIN NUMBER X MAXIMUM
!   GRAIN SIZE) IS FILLED OR RESERVES ARE EXHAUSTED
!
!   IDAYs1(7,=start of grain filling and setting max seed size
!   WTGRB=total seed C mass
!   GRWTB=individual seed size
!   GRNOB=seed set number
!   GROLM=maximum grain fill rate
!   GFILL=grain filling rate at 25 oC from PFT file
!   TFN3=temperature function for canopy growth
!   TFN4=temperature function for root growth
!
  IF(IDAYs1(7,NB,NZ).NE.0)THEN
    IF(WTGRBs1(NB,NZ).GE.GRWTBs1(NB,NZ)*GRNOBs1(NB,NZ))THEN
      GROLM=0._r8
    ELSEIF(IRTYPs1(NZ).EQ.0)THEN
      GROLM=AMAX1(0.0,GFILLs1(NZ)*GRNOBs1(NB,NZ)*SQRT(TFN3s1(NZ)))
    ELSE
      GROLM=AMAX1(0.0,GFILLs1(NZ)*GRNOBs1(NB,NZ) &
        *SQRT(TFN4s1(NGs1(NZ),NZ)))
    ENDIF
!
!     GRAIN FILL RATE MAY BE CONSTRAINED BY HIGH GRAIN C:N OR C:P
!
!     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     GROLM,GROLC=maximum,actual grain fill rate
!     XLOCM,XLOCC=maximum,actual C translocation rate from reserve to grain
!
    IF(WTGRBNs1(NB,NZ).LT.ZPGRM*CNGRs1(NZ) &
      *WTGRBs1(NB,NZ).OR.WTGRBPs1(NB,NZ).LT.ZPGRM &
      *CPGRs1(NZ)*WTGRBs1(NB,NZ))THEN
      GROLC=0._r8
    ELSE
      GROLC=GROLM
    ENDIF
    XLOCM=AMIN1(GROLM,WTRSVBs1(NB,NZ))
    XLOCC=AMIN1(GROLC,WTRSVBs1(NB,NZ))
!
!     GRAIN N OR P FILL RATE MAY BE LIMITED BY C:N OR C:P RATIOS
!     OF STALK RESERVES
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     ZNPGN,ZNPGP=effect of reserve N:C,P:C on grain fill N:C,P:C
!     SETN,SETP=Km for nonstructural N,P concn on seed set (g g-1)
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     ZPGRD=1.0_r8-ZPGRM
!     ZPGRN,ZPGRP=N:C,P:C ratios during grain fill
!     XLOCM,XLOCC=maximum,actual C translocation rate from reserve to grain
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     XLOCN,XLOCP=N,P translocation rate from reserve to grain
!
    IF(WTRSVBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      ZNPGN=WTRSBNs1(NB,NZ)/(WTRSBNs1(NB,NZ) &
        +SETN*WTRSVBs1(NB,NZ))
      ZNPGP=WTRSBPs1(NB,NZ)/(WTRSBPs1(NB,NZ) &
        +SETP*WTRSVBs1(NB,NZ))
      ZPGRN=ZPGRM+ZPGRD*AMAX1(0.0,AMIN1(1.0,ZNPGN))
      ZPGRP=ZPGRM+ZPGRD*AMAX1(0.0,AMIN1(1.0,ZNPGP))
      XLOCN=AMIN1(XLOCM*CNGRs1(NZ) &
        ,AMAX1(0.0,WTRSBNs1(NB,NZ)*ZPGRN) &
        ,(WTGRBs1(NB,NZ)+XLOCC)*CNGRs1(NZ)-WTGRBNs1(NB,NZ))
      XLOCP=AMIN1(XLOCM*CPGRs1(NZ) &
        ,AMAX1(0.0,WTRSBPs1(NB,NZ)*ZPGRP) &
        ,(WTGRBs1(NB,NZ)+XLOCC)*CPGRs1(NZ)-WTGRBPs1(NB,NZ))
    ELSE
      XLOCN=0._r8
      XLOCP=0._r8
    ENDIF
!     IF(NX.EQ.1.AND.NY.EQ.6.AND.NZ.EQ.3)THEN
!     WRITE(*,85)'XLOC',I,J,NZ,NB,WTGRBs1(NB,NZ),WTGRBNs1(NB,NZ)
!    2,WTRSVBs1(NB,NZ),WTRSBNs1(NB,NZ),XLOCC,XLOCN,XLOCP,XLOCM
!    3,CNGRs1(NZ),ZPGRX,ZNPG,GROLC,GROLM,GROGR,GROGRN
!    3,XLOCM*CNGRs1(NZ),AMAX1(0.0,WTRSBNs1(NB,NZ)*ZPGRX)
!    4,(WTGRBs1(NB,NZ)+XLOCC)*CNGRs1(NZ)-WTGRBNs1(NB,NZ)
!    4,GRNOBs1(NB,NZ),GRWTBs1(NB,NZ),GFILLs1(NZ)
!    5,SQRT(TFN3s1(NZ)),FLG4s1(NB,NZ)
!85    FORMAT(A8,4I4,20E12.4)
!     ENDIF
!
!     TRANSLOCATE C,N,P FROM STALK RESERVES TO GRAIN
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     GROGR=grain growth rate
!     XLOCC,XLOCN,XLOCP=C,N,P translocation rate from reserve to grain
!
    WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+GROGR-XLOCC
    WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+GROGRN-XLOCN
    WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+GROGRP-XLOCP
    WTGRBs1(NB,NZ)=WTGRBs1(NB,NZ)+XLOCC
    WTGRBNs1(NB,NZ)=WTGRBNs1(NB,NZ)+XLOCN
    WTGRBPs1(NB,NZ)=WTGRBPs1(NB,NZ)+XLOCP
  ELSE
    XLOCC=0._r8
    XLOCN=0._r8
    XLOCP=0._r8
  ENDIF
!
!   SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
!   HAS STOPPED FOR SET PERIOD OF TIME
!
!   IDAYs1(8,=end date setting for final seed number
!   XLOCC=C translocation rate from reserve to grain
!   PP=PFT population
!   FLG4=number of hours with no grain fill
!   FLG4X=number of hours with no grain filling until physl maturity
!   IDAYs1(10,=date of physiological maturity
!
  IF(IDAYs1(8,NB,NZ).NE.0)THEN
    IF(XLOCC.LE.1.0E-09*PPs1(NZ))THEN
      FLG4s1(NB,NZ)=FLG4s1(NB,NZ)+1.0
    ELSE
      FLG4s1(NB,NZ)=0._r8
    ENDIF
    IF(FLG4s1(NB,NZ).GE.FLG4X)THEN
      IF(IDAYs1(10,NB,NZ).EQ.0)THEN
        IDAYs1(10,NB,NZ)=I
      ENDIF
    ENDIF
!
!     TERMINATE ANNUALS AFTER GRAIN FILL
!
!     ISTYP=growth habit:0=annual,1=perennial
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     FLG4=number of hours with no grain fill
!     FLG4X=number of hours with no grain filling until physl maturity
!     FLG4Y=number of hours after physiol maturity required for senescence
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!
    IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
      IF(FLG4s1(NB,NZ).GT.FLG4X+FLG4Y(IWTYPs1(NZ)))THEN
        VRNFs1(NB,NZ)=VRNXs1(NB,NZ)+0.5
      ENDIF
    ENDIF
  ENDIF
  end subroutine GrainFilling

!------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8) :: dangle
  integer :: L,K,N
  ! begin_execution
  DO 900 K=1,JNODS1
    DO  L=1,JC1
      DO  N=1,JLI1
        SURFs1(N,L,K,NB,NZ)=0._r8
      enddo
    enddo
900 CONTINUE
! ARLFXB=0._r8
! ARLFXL=0._r8
! SURFXX=0._r8
  DO 500 K=1,JNODS1
!     ARLFXB=ARLFXB+ARLF1s1(K,NB,NZ)
    IF(ARLF1s1(K,NB,NZ).GT.0.0)THEN
      DO 700 L=JC1,1,-1
!       ARLFXL=ARLFXL+ARLFLs1(L,K,NB,NZ)
        DO 800 N=1,JLI1
          SURFs1(N,L,K,NB,NZ)=AMAX1(0.0_r8,CLASSs1(N,NZ) &
            *0.25_r8*ARLFLs1(L,K,NB,NZ))
  !       SURFXX=SURFXX+SURFs1(N,L,K,NB,NZ)

800     CONTINUE
700   CONTINUE
    ENDIF
500 CONTINUE
!
! ALLOCATE STALK AREA TO INCLINATION CLASSES ACCORDING TO
! BRANCH ANGLE ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
!
! SURFB=stalk surface area in canopy layer
! ANGBR=stem angle from horizontal
! ARSTK=total branch stalk surface area in each layer
!
  DO 910 L=1,JC1
    DO  N=1,JLI1
      SURFBs1(N,L,NB,NZ)=0._r8
    enddo
910  CONTINUE

  IF(NB.EQ.NB1s1(NZ))THEN
    N=JLI1
  ELSE
    dangle=PICON2/real(JLI1,r8)
    N=MIN(JLI1,INT(ASIN(ANGBRs1(NZ))/dangle)+1)
  ENDIF
  DO 710 L=JC1,1,-1
    SURFBs1(N,L,NB,NZ)=ARSTKs1(L,NB,NZ)/real(JLI1,r8)
710   CONTINUE
  end subroutine LeafClassAllocation
!------------------------------------------------------------------------------------------
  subroutine CarbNutInBranchTransfer(I,J,NB,NZ,IFLGZ,WFNG,WFNSG)
  implicit none
  integer, intent(in) :: I,J,NB,NZ,IFLGZ
  real(r8), intent(in) :: WFNG,WFNSG
  integer :: L
  real(r8) :: ZPOOLM,ZPOOLD
  real(r8) :: XFRPX,XFRCX,XFRNX
  real(r8) :: ATRPPD
  REAL(R8) :: cpoolt
  real(r8) :: CPOOLM,CH2OH
  real(r8) :: CWTRSV
  REAL(R8) :: CWTRSN,CWTRSP
  real(r8) :: CNR,CPR
  real(r8) :: CNL,CPL
  real(r8) :: CPOOLD
  real(r8) :: DATRP
  real(r8) :: FXFC,FXFN
  real(r8) :: FWTBR
  real(r8) :: GFNX
  real(r8) :: PPOOLD
  real(r8) :: PPDX
  real(r8) :: PPOOLM
  real(r8) :: UPNH4B,UPPO4B
  real(r8) :: UPNH4R,UPPO4R
  real(r8) :: WTRTM,WFNSP
  real(r8) :: WTPLTT,WTRTRX
  real(r8) :: WTPLTX,WVSTBX
  real(r8) :: WTRTTX,WTRSBX
  real(r8) :: WTRVCX
  real(r8) :: XFRC,XFRN,XFRP
  ! begin_execution
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!
!   IF(NZ.EQ.2)THEN
!   WRITE(*,2322)'VRNS',I,J,NX,NY,NZ,NB,NB1s1(NZ),IFLGZ
!    2,ISTYPs1(NZ),IFLGIs1(NZ),IFLGEs1(NB,NZ)
!    3,IFLGFs1(NB,NZ),IDAY0s1(NZ),IYR0s1(NZ)
!    3,VRNSs1(NB1s1(NZ),NZ),VRNLs1(NB,NZ)
!    3,VRNFs1(NB,NZ),VRNXs1(NB,NZ)
!2322  FORMAT(A8,14I4,20E12.4)
!   ENDIF
  IF((ISTYPs1(NZ).EQ.0.AND.IFLGIs1(NZ).EQ.0) &
    .OR.(I.GE.IDAY0s1(NZ).AND.IYRCs1.EQ.IYR0s1(NZ) &
    .AND.VRNFs1(NB,NZ) &
    .LT.FVRNs1(IWTYPs1(NZ))*VRNXs1(NB,NZ)) &
    .OR.(VRNSs1(NB1s1(NZ),NZ).GE.VRNLs1(NB,NZ) &
    .AND.VRNFs1(NB,NZ) &
    .LT.FVRNs1(IWTYPs1(NZ))*VRNXs1(NB,NZ)))THEN
    WTRTM=0._r8
    CPOOLM=0._r8
    DO 4 L=NUs1,NIs1(NZ)
      WTRTM=WTRTM+AMAX1(0.0,WTRTDs1(1,L,NZ))
      CPOOLM=CPOOLM+AMAX1(0.0,CPOOLRs1(1,L,NZ))
4   CONTINUE
!
  ! RESET TIME COUNTER
  !
  ! ATRP=hourly leafout counter
  ! IFLGA=flag for initializing leafout
  !
    IF(IFLGAs1(NB,NZ).EQ.0)THEN
      ATRPs1(NB,NZ)=0._r8
      IFLGAs1(NB,NZ)=1
    ENDIF
  !
  ! INCREMENT TIME COUNTER
  !
  ! IPTYP=photoperiod type:0=day neutral,1=short day,2=long day
  ! IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  ! XDL=critical photoperiod (h):<0=maximum daylength from site file
  ! XPPD=photoperiod sensitivity (node h-1)
  ! DYLN=daylength
  ! WFNSG=expansion,extension function of canopy water potential
  ! TFN3=temperature function for canopy growth
  ! ATRPX=number of hours required to initiate remobilization of storage C for leafout
  !
    IF(NB.EQ.NB1s1(NZ))THEN
      IF(IPTYPs1(NZ).EQ.2.AND.(IWTYPs1(NZ).EQ.1.OR.IWTYPs1(NZ).EQ.3))THEN
        PPDX=AMAX1(0.0,XDLs1(NZ)-XPPDs1(NZ)-DYLNs1)
        ATRPPD=EXP(-0.0*PPDX)
      ELSE
        ATRPPD=1.0_r8
      ENDIF
      IF(IGTYPs1(NZ).NE.0)THEN
        WFNSP=WFNSG
      ELSE
        WFNSP=1.0_r8
      ENDIF
      DATRP=ATRPPD*TFN3s1(NZ)*WFNSP
      ATRPs1(NB,NZ)=ATRPs1(NB,NZ)+DATRP
!     IF(NZ.EQ.2)THEN
!     WRITE(*,2323)'ATRP',I,J,NX,NY,NZ,NB,ATRPs1(NB,NZ),DATRP
!    2,ATRPPD,TFN3s1(NZ),WFNSG,PPDX,XDLs1(NZ),XPPDs1(NZ)
!    3,DYLNs1,WTLFBs1(NB,NZ),ARLFBs1(NB,NZ)
!    4,HTCTLs1(NZ)
!2323  FORMAT(A8,6I4,20E12.4)
!     ENDIF
      IF(ATRPs1(NB,NZ).LE.ATRPX(ISTYPs1(NZ)) &
        .OR.(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).EQ.0))THEN
        IF(WTRVCs1(NZ).GT.ZEROPs1(NZ))THEN
          CPOOLT=CPOOLM+CPOOLs1(NB,NZ)
  !
  !       REMOBILIZE C FROM SEASONAL STORAGE AT FIRST-ORDER RATE
  !       MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
  !
    !     GVMX=specific oxidation rate of storage C during leafout at 25 C
    !     WTRVC=storage C
    !     CH2OH=storage C oxidation rate during leafout
    !     CPOOL,CPOOLR=non-structural C mass in branch,root
    !     FXSH,FXRT=shoot-root partitioning of storage C during leafout
    !     WTRTD=root C mass
!
          GFNX=GVMX(ISTYPs1(NZ))*DATRP
          CH2OH=AMAX1(0.0,GFNX*WTRVCs1(NZ))
      !   IF(NZ.EQ.2)THEN
      !   WRITE(*,2123)'GERM0',I,J,NX,NY,NZ,NB
      !    2,GFNX,CH2OH,WTRVCs1(NZ)
      !    2,CPOOLs1(NB,NZ),CPOOLRs1(1,NGs1(NZ),NZ)
      !    3,FXSH(ISTYPs1(NZ)),FXRT(ISTYPs1(NZ))
      !2123  FORMAT(A8,6I4,20E12.4)
      !   ENDIF
          WTRVCs1(NZ)=WTRVCs1(NZ)-CH2OH
          CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+CH2OH*FXSH(ISTYPs1(NZ))
          IF(WTRTM.GT.ZEROPs1(NZ).AND.CPOOLM.GT.ZEROPs1(NZ))THEN
            DO 50 L=NUs1,NIs1(NZ)
              FXFC=AMAX1(0.0,WTRTDs1(1,L,NZ))/WTRTM
              CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ) &
                +FXFC*CH2OH*FXRT(ISTYPs1(NZ))
50          CONTINUE
          ELSE
            CPOOLRs1(1,NGs1(NZ),NZ)=CPOOLRs1(1,NGs1(NZ),NZ)+CH2OH*FXRT(ISTYPs1(NZ))
          ENDIF
        ELSE
          CH2OH=0._r8
        ENDIF
      ELSE
        CH2OH=0._r8
      ENDIF
      !
      !     REMOBILIZE N,P FROM SEASONAL STORAGE AT FIRST-ORDER RATE
      !     MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
      !
      !     WTRVC,WTRVN,WTRVP=storage C,N,P
      !     ISTYP=growth habit:0=annual,1=perennial from PFT file
      !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
      !     UPNH4B,UPPO4B=N,P transfer from storage to shoot
      !     CH2OH=storage C oxidation rate during leafout
      !     FRSV=rate constant for remobiln of storage C,N,P during leafout C
      !     FXSH=shoot partitioning of storage C during leafout
      !
      IF(WTRVCs1(NZ).GT.ZEROPs1(NZ))THEN
        IF(ISTYPs1(NZ).NE.0)THEN
          CPOOLT=AMAX1(0.0,WTRVCs1(NZ)+CPOOLs1(NB,NZ))
          ZPOOLD=(WTRVNs1(NZ)*CPOOLs1(NB,NZ)-ZPOOLs1(NB,NZ)*WTRVCs1(NZ))/CPOOLT
          PPOOLD=(WTRVPs1(NZ)*CPOOLs1(NB,NZ)-PPOOLs1(NB,NZ)*WTRVCs1(NZ))/CPOOLT
          UPNH4B=AMAX1(0.0,FRSV(IBTYPs1(NZ))*ZPOOLD)
          UPPO4B=AMAX1(0.0,FRSV(IBTYPs1(NZ))*PPOOLD)
        ELSE
          UPNH4B=AMAX1(0.0,FXSH(ISTYPs1(NZ))*CH2OH*WTRVNs1(NZ)/WTRVCs1(NZ))
          UPPO4B=AMAX1(0.0,FXSH(ISTYPs1(NZ))*CH2OH*WTRVPs1(NZ)/WTRVCs1(NZ))
        ENDIF
      ELSE
        UPNH4B=AMAX1(0.0,FXSH(ISTYPs1(NZ))*WTRVNs1(NZ))
        UPPO4B=AMAX1(0.0,FXSH(ISTYPs1(NZ))*WTRVPs1(NZ))
      ENDIF
    !
    ! ADD TO NON-STRUCTURAL POOLS IN ROOT
    !
    ! CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
    ! WTRVC,WTRVN,WTRVP=storage C,N,P
    ! ISTYP=growth habit:0=annual,1=perennial from PFT file
    ! UPNH4R,UPPO4R=N,P transfer from storage to root
    ! FRSV=rate constant for remobiln of storage C,N,P during leafout
    ! FXRT=root partitioning of storage C during leafout
    !
      CPOOLM=0._r8
      ZPOOLM=0._r8
      PPOOLM=0._r8
      DO 3 L=NUs1,NIs1(NZ)
        CPOOLM=CPOOLM+AMAX1(0.0,CPOOLRs1(1,L,NZ))
        ZPOOLM=ZPOOLM+AMAX1(0.0,ZPOOLRs1(1,L,NZ))
        PPOOLM=PPOOLM+AMAX1(0.0,PPOOLRs1(1,L,NZ))
3     CONTINUE
      IF(WTRVCs1(NZ).GT.ZEROPs1(NZ))THEN
        IF(ISTYPs1(NZ).NE.0)THEN
          CPOOLT=AMAX1(ZEROPs1(NZ),WTRVCs1(NZ)+CPOOLM)
          ZPOOLD=(WTRVNs1(NZ)*CPOOLM-ZPOOLM*WTRVCs1(NZ))/CPOOLT
          PPOOLD=(WTRVPs1(NZ)*CPOOLM-PPOOLM*WTRVCs1(NZ))/CPOOLT
          UPNH4R=AMAX1(0.0,FRSV(IBTYPs1(NZ))*ZPOOLD)
          UPPO4R=AMAX1(0.0,FRSV(IBTYPs1(NZ))*PPOOLD)
!         IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!         WRITE(*,9878)'GERM1',I,J,NZ,UPNH4R,FRSV(IBTYPs1(NZ))
!        2,ZPOOLD,WTRVNs1(NZ),CPOOLM,ZPOOLM,WTRVCs1(NZ)
!        3,CPOOLT
!9878     FORMAT(A8,3I4,12E24.16)
!         ENDIF
        ELSE
          UPNH4R=AMAX1(0.0,FXRT(ISTYPs1(NZ))*CH2OH*WTRVNs1(NZ)/WTRVCs1(NZ))
          UPPO4R=AMAX1(0.0,FXRT(ISTYPs1(NZ))*CH2OH*WTRVPs1(NZ)/WTRVCs1(NZ))
        ENDIF
      ELSE
        UPNH4R=AMAX1(0.0,FXRT(ISTYPs1(NZ))*WTRVNs1(NZ))
        UPPO4R=AMAX1(0.0,FXRT(ISTYPs1(NZ))*WTRVPs1(NZ))
      ENDIF
!
!     TRANSFER STORAGE FLUXES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     UPNH4B,UPPO4B=N,P transfer from storage to shoot
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     UPNH4R,UPPO4R=N,P transfer from storage to root
!     FXFN=root layer allocation
!
      WTRVNs1(NZ)=WTRVNs1(NZ)-UPNH4B-UPNH4R
      WTRVPs1(NZ)=WTRVPs1(NZ)-UPPO4B-UPPO4R
      ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+UPNH4B
      PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+UPPO4B
      IF(WTRTM.GT.ZEROPs1(NZ).AND.CPOOLM.GT.ZEROPs1(NZ))THEN
        DO 51 L=NUs1,NIs1(NZ)
          FXFN=AMAX1(0.0,CPOOLRs1(1,L,NZ))/CPOOLM
!         IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!         WRITE(*,9879)'GERM2',I,J,NZ,L,UPNH4R,FXFN
!        2,ZPOOLRs1(1,L,NZ),CPOOLRs1(1,L,NZ),CPOOLM
!9879      FORMAT(A8,4I4,12E24.16)
!         ENDIF
          ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)+FXFN*UPNH4R
          PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)+FXFN*UPPO4R
51      CONTINUE
      ELSE
  !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
  !     WRITE(*,9879)'GERM3',I,J,NZ,L,UPNH4R,FXFN
  !    2,ZPOOLRs1(1,L,NZ),CPOOLRs1(1,L,NZ),CPOOLM
  !     ENDIF
        ZPOOLRs1(1,NGs1(NZ),NZ)=ZPOOLRs1(1,NGs1(NZ),NZ)+UPNH4R
        PPOOLRs1(1,NGs1(NZ),NZ)=PPOOLRs1(1,NGs1(NZ),NZ)+UPPO4R
      ENDIF
    ENDIF
  !
  ! REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
  !
  ! ATRP=hourly leafout counter
  ! TFN3=temperature function for canopy growth
  ! ATRPX=number of hours required for remobilization of storage C during leafout
  ! WFNG=growth function of canopy water potential
  ! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  ! XFRC,XFRN,XFRC=nonstructural C,N,P transfer
  !
    IF(NB.NE.NB1s1(NZ).AND.ATRPs1(NB,NZ) &
      .LE.ATRPX(ISTYPs1(NZ)))THEN
      ATRPs1(NB,NZ)=ATRPs1(NB,NZ)+TFN3s1(NZ)*WFNG
      XFRC=AMAX1(0.0,0.05*TFN3s1(NZ) &
        *(0.5*(CPOOLs1(NB1s1(NZ),NZ)+CPOOLs1(NB,NZ))-CPOOLs1(NB,NZ)))
      XFRN=AMAX1(0.0,0.05*TFN3s1(NZ) &
        *(0.5*(ZPOOLs1(NB1s1(NZ),NZ)+ZPOOLs1(NB,NZ))-ZPOOLs1(NB,NZ)))
      XFRP=AMAX1(0.0,0.05*TFN3s1(NZ) &
      *(0.5*(PPOOLs1(NB1s1(NZ),NZ)+PPOOLs1(NB,NZ))-PPOOLs1(NB,NZ)))
      CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+XFRC
      ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+XFRN
      PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+XFRP
      CPOOLs1(NB1s1(NZ),NZ)=CPOOLs1(NB1s1(NZ),NZ)-XFRC
      ZPOOLs1(NB1s1(NZ),NZ)=ZPOOLs1(NB1s1(NZ),NZ)-XFRN
      PPOOLs1(NB1s1(NZ),NZ)=PPOOLs1(NB1s1(NZ),NZ)-XFRP
    ENDIF
  ENDIF

!
! TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
! IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
! IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
!
! ISTYP=growth habit:0=annual,1=perennial from PFT file
! IFLGZ=remobilization flag
! WVSTKB=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! FXFB=rate constant for plant-storage nonstructural C,N,P exchange
! XFRC,XFRN,XFRC=nonstructural C,N,P transfer
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! ISTYP=growth habit:0=annual,1=perennial from PFT file
! IFLGZ=remobilization flag
! WVSTKB=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! XFRC,XFRN,XFRC=nonstructural C,N,P transfer
! WTRVC,WTRVN,WTRVP=storage C,N,P
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!
  IF(IFLGZ.EQ.1.AND.ISTYPs1(NZ).NE.0)THEN
    IF(WVSTKBs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.WTRSVBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      CWTRSV=AMAX1(0.0,WTRSVBs1(NB,NZ)/WVSTKBs1(NB,NZ))
      CWTRSN=AMAX1(0.0,WTRSBNs1(NB,NZ)/WVSTKBs1(NB,NZ))
      CWTRSP=AMAX1(0.0,WTRSBPs1(NB,NZ)/WVSTKBs1(NB,NZ))
      CNR=CWTRSV/(CWTRSV+CWTRSN/CNKI)
      CPR=CWTRSV/(CWTRSV+CWTRSP/CPKI)
    ELSE
      CNR=0._r8
      CPR=0._r8
    ENDIF
    XFRCX=FXFB(IBTYPs1(NZ))*AMAX1(0.0,WTRSVBs1(NB,NZ))
    XFRNX=FXFB(IBTYPs1(NZ))*AMAX1(0.0,WTRSBNs1(NB,NZ))*(1.0+CNR)
    XFRPX=FXFB(IBTYPs1(NZ))*AMAX1(0.0,WTRSBPs1(NB,NZ))*(1.0+CPR)
    XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
    WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)-XFRC
    WTRVCs1(NZ)=WTRVCs1(NZ)+XFRC
    WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)-XFRN
    WTRVNs1(NZ)=WTRVNs1(NZ)+XFRN
    WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)-XFRP
    WTRVPs1(NZ)=WTRVPs1(NZ)+XFRP
    IF(CCPOLBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      CNL=CCPOLBs1(NB,NZ)/(CCPOLBs1(NB,NZ)+CZPOLBs1(NB,NZ)/CNKI)
      CPL=CCPOLBs1(NB,NZ)/(CCPOLBs1(NB,NZ)+CPPOLBs1(NB,NZ)/CPKI)
    ELSE
      CNL=0._r8
      CPL=0._r8
    ENDIF
    XFRCX=FXFB(IBTYPs1(NZ))*AMAX1(0.0,CPOOLs1(NB,NZ))
    XFRNX=FXFB(IBTYPs1(NZ))*AMAX1(0.0,ZPOOLs1(NB,NZ))*(1.0+CNL)
    XFRPX=FXFB(IBTYPs1(NZ))*AMAX1(0.0,PPOOLs1(NB,NZ))*(1.0+CPL)
    XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
    CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)-XFRC
    WTRVCs1(NZ)=WTRVCs1(NZ)+XFRC
    ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-XFRN
    WTRVNs1(NZ)=WTRVNs1(NZ)+XFRN
    PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-XFRP
    WTRVPs1(NZ)=WTRVPs1(NZ)+XFRP
!   IF(NZ.EQ.1)THEN
!     WRITE(*,4490)'RSV',I,J,NZ,NB,XFRC,XFRN,WTRSVBs1(NB,NZ)
!    2,WTRSBNs1(NB,NZ),WTRVCs1(NZ),WTRVNs1(NZ)
!    3,CNR,CNL,CPOOLs1(NB,NZ),ZPOOLs1(NB,NZ)
!    4,FXFB(IBTYPs1(NZ))
!4490  FORMAT(A8,4I4,20E12.4)
!     ENDIF
  ENDIF
!
!   TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES AND ROOTS TO RESERVES
!   IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN STALK RESERVES
!   AND LEAVES IN PERENNIALS ACCORDING TO CONCENTRATION DIFFERENCES
!
!   ISTYP=growth habit:0=annual,1=perennial from PFT file
!   IDAYs1(3,=start of stem elongation and setting max seed number
!   IDAYs1(8,=end date setting for final seed number
!   WTLSB=leaf+petiole mass
!   WVSTKB=stalk sapwood mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   FXFY,FXFZ=rate constant for plant-reserve nonstructural C,N,P exchange
!   XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!   CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  IF((ISTYPs1(NZ).EQ.0.AND.IDAYs1(8,NB,NZ).NE.0) &
    .OR.(ISTYPs1(NZ).EQ.1.AND.IDAYs1(3,NB,NZ).NE.0))THEN
    WTPLTT=WTLSBs1(NB,NZ)+WVSTKBs1(NB,NZ)
    CPOOLT=CPOOLs1(NB,NZ)+WTRSVBs1(NB,NZ)
    IF(WTPLTT.GT.ZEROPs1(NZ))THEN
      CPOOLD=(CPOOLs1(NB,NZ)*WVSTKBs1(NB,NZ) &
        -WTRSVBs1(NB,NZ)*WTLSBs1(NB,NZ))/WTPLTT
      XFRC=FXFY(ISTYPs1(NZ))*CPOOLD
      CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)-XFRC
      WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+XFRC
    ENDIF
    IF(CPOOLT.GT.ZEROPs1(NZ))THEN
      ZPOOLD=(ZPOOLs1(NB,NZ)*WTRSVBs1(NB,NZ) &
        -WTRSBNs1(NB,NZ)*CPOOLs1(NB,NZ))/CPOOLT
      PPOOLD=(PPOOLs1(NB,NZ)*WTRSVBs1(NB,NZ) &
        -WTRSBPs1(NB,NZ)*CPOOLs1(NB,NZ))/CPOOLT
      XFRN=FXFZ(ISTYPs1(NZ))*ZPOOLD
      XFRP=FXFZ(ISTYPs1(NZ))*PPOOLD
      ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-XFRN
      WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+XFRN
      PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-XFRP
      WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+XFRP
    ENDIF
!     IF(NZ.EQ.1)THEN
!     WRITE(*,4488)'EXCHC',I,J,NX,NY,NZ,NB,NS,XFRC,XFRN
!    2,FXFZ(ISTYPs1(NZ)),WTRSVBs1(NB,NZ),CPOOLs1(NB,NZ)
!    3,WVSTKBs1(NB,NZ),WTLSBs1(NB,NZ)
!    4,CPOOLT,CPOOLD,ZPOOLs1(NB,NZ),WTRSBNs1(NB,NZ)
!4488  FORMAT(A8,7I4,12E12.4)
!     ENDIF
    IF(ISTYPs1(NZ).EQ.0.AND.IDAYs1(8,NB,NZ).NE.0)THEN
      DO 2050 L=NUs1,NIs1(NZ)
        IF(VOLXs1(L).GT.ZEROS2s1)THEN
          WTRTRX=AMAX1(ZEROPs1(NZ),WTRTLs1(1,L,NZ)*FWODRs1(1))
          WTPLTX=WTRTRX+WVSTKBs1(NB,NZ)
          IF(WTPLTX.GT.ZEROPs1(NZ))THEN
            CPOOLD=(CPOOLRs1(1,L,NZ)*WVSTKBs1(NB,NZ)-WTRSVBs1(NB,NZ)*WTRTRX)/WTPLTX
            XFRC=AMAX1(0.0,FXFY(ISTYPs1(NZ))*CPOOLD)
            CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)-XFRC
            WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+XFRC
            CPOOLT=CPOOLRs1(1,L,NZ)+WTRSVBs1(NB,NZ)
            IF(CPOOLT.GT.ZEROPs1(NZ))THEN
              ZPOOLD=(ZPOOLRs1(1,L,NZ)*WTRSVBs1(NB,NZ) &
                -WTRSBNs1(NB,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              PPOOLD=(PPOOLRs1(1,L,NZ)*WTRSVBs1(NB,NZ) &
                -WTRSBPs1(NB,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRN=AMAX1(0.0,FXFZ(ISTYPs1(NZ))*ZPOOLD)
              XFRP=AMAX1(0.0,FXFZ(ISTYPs1(NZ))*PPOOLD)
              ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)-XFRN
              WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+XFRN
              PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)-XFRP
              WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+XFRP
          !     IF(NZ.EQ.1)THEN
          !     WRITE(*,4489)'EXCHC',I,J,NZ,NB,L,WTRSVBs1(NB,NZ)
          !    2,WVSTKBs1(NB,NZ),CPOOLRs1(1,L,NZ)
          !    3,WTRTLs1(1,L,NZ),FWOODs1(1),WTRTRX,WTPLTX
          !    4,CPOOLT,CPOOLD,XFRC,FXFZ(ISTYPs1(NZ))
          !4489  FORMAT(A8,5I4,12E16.8)
          !     ENDIF
          !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
          !     WRITE(*,4489)'EXCHN',I,J,NZ,NB,L,WTRSBNs1(NB,NZ)
          !    2,WTRSVBs1(NB,NZ),ZPOOLRs1(1,L,NZ)
          !    3,CPOOLRs1(1,L,NZ),FWOODs1(1),ZPOOLD,XFRN
          !     ENDIF
            ENDIF
          ENDIF
        ENDIF
2050  CONTINUE
    ENDIF
  ENDIF
!
!   REPLENISH BRANCH NON-STRUCTURAL POOL FROM
!   SEASONAL STORAGE POOL
!
!   WVSTKB,WVSTK=stalk,total stalk sapwood mass
!   WTRT=total root mass
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   XFRX=maximum storage C content for remobiln from stalk,root reserves
!   XFRC=C transfer
!
  IF(WVSTKBs1(NB,NZ).GT.ZEROPs1(NZ) &
    .AND.WVSTKs1(NZ).GT.ZEROPs1(NZ) &
    .AND.WTRTs1(NZ).GT.ZEROPs1(NZ) &
    .AND.WTRSVBs1(NB,NZ).LE.XFRX*WVSTKBs1(NB,NZ))THEN
    FWTBR=WVSTKBs1(NB,NZ)/WVSTKs1(NZ)
    WVSTBX=WVSTKBs1(NB,NZ)
    WTRTTX=WTRTs1(NZ)*FWTBR
    WTPLTT=WVSTBX+WTRTTX
    WTRSBX=AMAX1(0.0,WTRSVBs1(NB,NZ))
    WTRVCX=AMAX1(0.0,WTRVCs1(NZ)*FWTBR)
    CPOOLD=(WTRVCX*WVSTBX-WTRSBX*WTRTTX)/WTPLTT
    XFRC=AMAX1(0.0,XFRY*CPOOLD)
    WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+XFRC
    WTRVCs1(NZ)=WTRVCs1(NZ)-XFRC
  ENDIF
end subroutine CarbNutInBranchTransfer
!------------------------------------------------------------------------------------------
  subroutine StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,CNSHW,&
    CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,WFNC,WFNS,WFNSG)
  integer, intent(in) :: I,J,NZ
  integer, intent(out):: ICHK1(2,JZ1)
  integer, intent(out):: NRX(2,JZ1)
  REAL(R8), INTENT(OUT):: TFN6(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,WFNC
  real(r8), intent(out) :: WFNS,WFNSG
  integer :: L,NR,N
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM
!     begin_execution

  DO 2 L=1,JC1
    ARLFVs1(L,NZ)=0._r8
    WGLFVs1(L,NZ)=0._r8
    ARSTVs1(L,NZ)=0._r8
2 CONTINUE
  DO 5 NR=1,NRTs1(NZ)
    DO  N=1,MYs1(NZ)
      NRX(N,NR)=0
      ICHK1(N,NR)=0
    enddo
5 CONTINUE
  DO 9 N=1,MYs1(NZ)
    DO 6 L=NUs1,NJs1
      WSRTLs1(N,L,NZ)=0._r8
      RTN1s1(N,L,NZ)=0._r8
      RTNLs1(N,L,NZ)=0._r8
      RCO2Ms1(N,L,NZ)=0._r8
      RCO2Ns1(N,L,NZ)=0._r8
      RCO2As1(N,L,NZ)=0._r8
6   CONTINUE
9 CONTINUE
!
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     WTSTK,WVSTK=stalk,sapwood mass
!     FWOOD,FWODB=C woody fraction in stalk,other organs:0=woody,1=non-woody
!     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
!     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!     *LF=leaf,*SHE=petiole,*STK=stalk,*RT=root
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
!
  IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1 &
    .OR.WTSTKs1(NZ).LE.ZEROPs1(NZ))THEN
    FWODBs1(1)=1.0_r8
    FWOODs1(1)=1.0_r8
    FWODRs1(1)=1.0_r8
  ELSE
    FWODBs1(1)=1.0_r8
    FWOODs1(1)=SQRT(WVSTKs1(NZ)/WTSTKs1(NZ))
    FWODRs1(1)=SQRT(FRTX*WVSTKs1(NZ)/WTSTKs1(NZ))
  ENDIF
  FWODBs1(0)=1.0_r8-FWODBs1(1)
  FWOODs1(0)=1.0_r8-FWOODs1(1)
  FWODRs1(0)=1.0_r8-FWODRs1(1)
  CNLFW=FWODBs1(0)*CNSTKs1(NZ)+FWODBs1(1)*CNLFs1(NZ)
  CPLFW=FWODBs1(0)*CPSTKs1(NZ)+FWODBs1(1)*CPLFs1(NZ)
  CNSHW=FWODBs1(0)*CNSTKs1(NZ)+FWODBs1(1)*CNSHEs1(NZ)
  CPSHW=FWODBs1(0)*CPSTKs1(NZ)+FWODBs1(1)*CPSHEs1(NZ)
  CNRTW=FWODRs1(0)*CNSTKs1(NZ)+FWODRs1(1)*CNRTs1(NZ)
  CPRTW=FWODRs1(0)*CPSTKs1(NZ)+FWODRs1(1)*CPRTs1(NZ)
  FWODLNs1(0)=FWODBs1(0)*CNSTKs1(NZ)/CNLFW
  FWODLPs1(0)=FWODBs1(0)*CPSTKs1(NZ)/CPLFW
  FWODSNs1(0)=FWODBs1(0)*CNSTKs1(NZ)/CNSHW
  FWODSPs1(0)=FWODBs1(0)*CPSTKs1(NZ)/CPSHW
  FWOODNs1(0)=FWOODs1(0)*CNSTKs1(NZ)/CNRTW
  FWOODPs1(0)=FWOODs1(0)*CPSTKs1(NZ)/CPRTW
  FWODRNs1(0)=FWODRs1(0)*CNRTs1(NZ)/CNRTW
  FWODRPs1(0)=FWODRs1(0)*CPRTs1(NZ)/CPRTW
  FWODLNs1(1)=1.0_r8-FWODLNs1(0)
  FWODLPs1(1)=1.0_r8-FWODLPs1(0)
  FWODSNs1(1)=1.0_r8-FWODSNs1(0)
  FWODSPs1(1)=1.0_r8-FWODSPs1(0)
  FWOODNs1(1)=1.0_r8-FWOODNs1(0)
  FWOODPs1(1)=1.0_r8-FWOODPs1(0)
  FWODRNs1(1)=1.0_r8-FWODRNs1(0)
  FWODRPs1(1)=1.0_r8-FWODRPs1(0)
!
!     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
!     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
!
!     TKC,TKCM=canopy temperature,canopy temp used in Arrhenius eqn
!     TKS,TKSM=soil temperature,soil temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN5,TFN6=temperature function for canopy,root mntc respn (25 oC =1)
!     8.3143,710.0=gas constant,enthalpy
!     62500,195000,232500=energy of activn,high,low temp inactivn(KJ mol-1)
!
  TKCM=TKCs1(NZ)+OFFSTs1(NZ)
  RTK=8.3143*TKCM
  STK=710.0*TKCM
  ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
  TFN5=EXP(25.214-62500/RTK)/ACTVM
  DO 7 L=NUs1,NJs1
    TKSM=TKSs1(L)+OFFSTs1(NZ)
    RTK=8.3143*TKSM
    STK=710.0*TKSM
    ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
    TFN6(L)=EXP(25.214-62500/RTK)/ACTVM
7 CONTINUE
!
!     PRIMARY ROOT NUMBER
!
!     WTRTA=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     XRTN1=multiplier for number of primary root axes
!
  WTRTAs1(NZ)=AMAX1(0.999992087*WTRTAs1(NZ),WTRTs1(NZ)/PPs1(NZ))
  XRTN1=AMAX1(1.0,WTRTAs1(NZ)**0.667)*PPs1(NZ)
!
!     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
!     FROM CANOPY TURGOR
!
!     WFNS=turgor expansion,extension function
!     PSILG,PSILM=current,minimum canopy turgor potl for expansion,extension
!     WFNC=stomatal resistance function of canopy turgor
!     PSILT=canopy water potential
!     WFNG=growth function of canopy water potential
!     WFNSG=expansion,extension function of canopy water potential
!
  WFNS=AMIN1(1.0,AMAX1(0.0,PSILGs1(NZ)-PSILM))
  IF(IGTYPs1(NZ).EQ.0)THEN
    WFNC=1.0_r8
    WFNG=EXP(0.05*PSILTs1(NZ))
    WFNSG=WFNS**0.10
  ELSE
    WFNC=EXP(RCSs1(NZ)*PSILGs1(NZ))
    WFNG=EXP(0.10*PSILTs1(NZ))
    WFNSG=WFNS**0.25
  ENDIF
  end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

  subroutine NonstructlBiomTransfer(I,J,NZ,PTRT,RLNT,RTSK1,RTSK2,RTNT,IFLGZ)
  implicit none
  integer, intent(in) :: I,J,NZ,IFLGZ
  real(r8), intent(in):: PTRT
  real(r8), INTENT(IN) :: RLNT(2,JZ1)
  real(r8),intent(in) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  real(r8),intent(in) :: RTNT(2)
  integer :: L,NB,N,NR
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD
  real(r8) :: XFRPX,XFRCX,XFRNX
  real(r8) :: WTLSBZ(JC1)
  real(r8) :: CPOOLZ(JC1),ZPOOLZ(JC1),PPOOLZ(JC1)
  REAL(R8) :: FWTR(JZ1),FWTB(JP1)
  real(r8) :: CPOOLT
  real(r8) :: CNL,CPL,CPOOLD
  real(r8) :: CPOOLB,CPOOLS
  real(r8) :: FWTC
  real(r8) :: FWTS
  real(r8) :: PPOOLB
  real(r8) :: PPOOLD
  real(r8) :: PPOOLT
  real(r8) :: PTSHTR
  real(r8) :: PPOOLS
  real(r8) :: WTPLTT
  real(r8) :: WTRTLX
  real(r8) :: WTRTD1
  real(r8) :: WTSTKT
  real(r8) :: WTRSVT,WTRSNT,WTRSPT
  real(r8) :: WTRSVD,WTRSND,WTRSPD
  real(r8) :: WTRTD2,WTLSBX,WTLSBB
  real(r8) :: WTRTLR
  real(r8) :: XFRC,XFRN,XFRP
!     begin_execution
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ATRP=hourly leafout counter
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     WTLSB=leaf+petiole mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!
  IF(NBRs1(NZ).GT.1)THEN
    WTPLTT=0._r8
    CPOOLT=0._r8
    ZPOOLT=0._r8
    PPOOLT=0._r8
    DO 300 NB=1,NBRs1(NZ)
      IF(IDTHBs1(NB,NZ).EQ.0)THEN
        IF(ATRPs1(NB,NZ).GT.ATRPX(ISTYPs1(NZ)))THEN
          WTLSBZ(NB)=AMAX1(0.0,WTLSBs1(NB,NZ))
          CPOOLZ(NB)=AMAX1(0.0,CPOOLs1(NB,NZ))
          ZPOOLZ(NB)=AMAX1(0.0,ZPOOLs1(NB,NZ))
          PPOOLZ(NB)=AMAX1(0.0,PPOOLs1(NB,NZ))
          WTPLTT=WTPLTT+WTLSBZ(NB)
          CPOOLT=CPOOLT+CPOOLZ(NB)
          ZPOOLT=ZPOOLT+ZPOOLZ(NB)
          PPOOLT=PPOOLT+PPOOLZ(NB)
        ENDIF
      ENDIF
300 CONTINUE
    DO 305 NB=1,NBRs1(NZ)
      IF(IDTHBs1(NB,NZ).EQ.0)THEN
        IF(ATRPs1(NB,NZ).GT.ATRPX(ISTYPs1(NZ)))THEN
          IF(WTPLTT.GT.ZEROPs1(NZ) &
            .AND.CPOOLT.GT.ZEROPs1(NZ))THEN
            CPOOLD=CPOOLT*WTLSBZ(NB)-CPOOLZ(NB)*WTPLTT
            ZPOOLD=ZPOOLT*CPOOLZ(NB)-ZPOOLZ(NB)*CPOOLT
            PPOOLD=PPOOLT*CPOOLZ(NB)-PPOOLZ(NB)*CPOOLT
            XFRC=0.01*CPOOLD/WTPLTT
            XFRN=0.01*ZPOOLD/CPOOLT
            XFRP=0.01*PPOOLD/CPOOLT
            CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+XFRC
            ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+XFRN
            PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+XFRP
          ENDIF
        ENDIF
      ENDIF
305 CONTINUE
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     WVSTKB=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IDAYs1(7,=start of grain filling and setting max seed size
!
  IF(NBRs1(NZ).GT.1)THEN
    WTSTKT=0._r8
    WTRSVT=0._r8
    WTRSNT=0._r8
    WTRSPT=0._r8
    DO 330 NB=1,NBRs1(NZ)
      IF(IDTHBs1(NB,NZ).EQ.0)THEN
        IF(IDAYs1(7,NB,NZ).NE.0)THEN
          WTSTKT=WTSTKT+WVSTKBs1(NB,NZ)
          WTRSVT=WTRSVT+WTRSVBs1(NB,NZ)
          WTRSNT=WTRSNT+WTRSBNs1(NB,NZ)
          WTRSPT=WTRSPT+WTRSBPs1(NB,NZ)
        ENDIF
      ENDIF
330 CONTINUE
    IF(WTSTKT.GT.ZEROPs1(NZ) &
      .AND.WTRSVT.GT.ZEROPs1(NZ))THEN
      DO 335 NB=1,NBRs1(NZ)
        IF(IDTHBs1(NB,NZ).EQ.0)THEN
          IF(IDAYs1(7,NB,NZ).NE.0)THEN
            WTRSVD=WTRSVT*WVSTKBs1(NB,NZ)-WTRSVBs1(NB,NZ)*WTSTKT
            XFRC=0.1*WTRSVD/WTSTKT
            WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+XFRC
            WTRSND=WTRSNT*WTRSVBs1(NB,NZ)-WTRSBNs1(NB,NZ)*WTRSVT
            XFRN=0.1*WTRSND/WTRSVT
            WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+XFRN
            WTRSPD=WTRSPT*WTRSVBs1(NB,NZ)-WTRSBPs1(NB,NZ)*WTRSVT
            XFRP=0.1*WTRSPD/WTRSVT
            WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+XFRP
          ENDIF
        ENDIF
335   CONTINUE
    ENDIF
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BWTWEEN ROOT AND MYCORRHIZAE
!     IN EACH ROOTED SOIL LAYER FROM NON-STRUCTURAL C,N,P
!     CONCENTRATION DIFFERENCES
!
!     MY=mycorrhizal:1=no,2=yes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in 1:root,2:mycorrhizae
!     WTRTD=1:root,2:mycorrhizal C mass
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
!
  IF(MYs1(NZ).EQ.2)THEN
    DO 425 L=NUs1,NIXs1(NZ)
      IF(CPOOLRs1(1,L,NZ).GT.ZEROPs1(NZ) &
        .AND.WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
        WTRTD1=WTRTDs1(1,L,NZ)
        WTRTD2=AMIN1(WTRTDs1(1,L,NZ),AMAX1(FSNK &
          *WTRTDs1(1,L,NZ),WTRTDs1(2,L,NZ)))
        WTPLTT=WTRTD1+WTRTD2
        IF(WTPLTT.GT.ZEROPs1(NZ))THEN
          CPOOLD=(CPOOLRs1(1,L,NZ)*WTRTD2 &
            -CPOOLRs1(2,L,NZ)*WTRTD1)/WTPLTT
          XFRC=FMYC*CPOOLD
          CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)-XFRC
          CPOOLRs1(2,L,NZ)=CPOOLRs1(2,L,NZ)+XFRC
          CPOOLT=CPOOLRs1(1,L,NZ)+CPOOLRs1(2,L,NZ)
          IF(CPOOLT.GT.ZEROPs1(NZ))THEN
            ZPOOLD=(ZPOOLRs1(1,L,NZ)*CPOOLRs1(2,L,NZ) &
              -ZPOOLRs1(2,L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
            XFRN=FMYC*ZPOOLD
            PPOOLD=(PPOOLRs1(1,L,NZ)*CPOOLRs1(2,L,NZ) &
              -PPOOLRs1(2,L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
            XFRP=FMYC*PPOOLD
            ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)-XFRN
            ZPOOLRs1(2,L,NZ)=ZPOOLRs1(2,L,NZ)+XFRN
            PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)-XFRP
            PPOOLRs1(2,L,NZ)=PPOOLRs1(2,L,NZ)+XFRP
    !     IF(L.EQ.NIXs1(NZ))THEN
    !     WRITE(*,9873)'MYCO',I,J,NZ,L,XFRC,XFRN,XFRP
!    2,CPOOLRs1(1,L,NZ),WTRTDs1(1,L,NZ)
!    3,CPOOLRs1(2,L,NZ),WTRTD2
!    3,WTPLTT,ZPOOLRs1(1,L,NZ),ZPOOLRs1(2,L,NZ)
!    4,PPOOLRs1(1,L,NZ),PPOOLRs1(2,L,NZ),CPOOLT
!9873  FORMAT(A8,4I4,20E24.16)
!     ENDIF
          ENDIF
        ENDIF
      ENDIF
425 CONTINUE
  ENDIF
!
!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
  IF(IFLGZ.EQ.1.AND.ISTYPs1(NZ).NE.0)THEN
    DO 5545 N=1,MYs1(NZ)
      DO 5550 L=NUs1,NIs1(NZ)
        IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
          CNL=CCPOLRs1(N,L,NZ)/(CCPOLRs1(N,L,NZ) &
            +CZPOLRs1(N,L,NZ)/CNKI)
          CPL=CCPOLRs1(N,L,NZ)/(CCPOLRs1(N,L,NZ) &
            +CPPOLRs1(N,L,NZ)/CPKI)
        ELSE
          CNL=0._r8
          CPL=0._r8
        ENDIF
        XFRCX=FXFR(IGTYPs1(NZ)) &
          *AMAX1(0.0,CPOOLRs1(N,L,NZ))
        XFRNX=FXFR(IGTYPs1(NZ)) &
          *AMAX1(0.0,ZPOOLRs1(N,L,NZ))*(1.0+CNL)
        XFRPX=FXFR(IGTYPs1(NZ)) &
          *AMAX1(0.0,PPOOLRs1(N,L,NZ))*(1.0+CPL)
        XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
        XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
        XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
        CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-XFRC
        WTRVCs1(NZ)=WTRVCs1(NZ)+XFRC
        ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-XFRN
        WTRVNs1(NZ)=WTRVNs1(NZ)+XFRN
        PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-XFRP
        WTRVPs1(NZ)=WTRVPs1(NZ)+XFRP
5550  CONTINUE
5545  CONTINUE
  ENDIF
!
!     ROOT AND NODULE TOTALS
!
!     WTRTL,WTRTD=active,actual root C mass
!     WTRT1,WTRT2=primary,secondary root C mass in soil layer
!     TCO2T=total PFT respiration
!     RCO2A=total root respiration
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
  DO 5445 N=1,MYs1(NZ)
    DO 5450 L=NUs1,NIs1(NZ)
      WTRTLs1(N,L,NZ)=0._r8
      WTRTDs1(N,L,NZ)=0._r8
      DO 5460 NR=1,NRTs1(NZ)
        WTRTLs1(N,L,NZ)=WTRTLs1(N,L,NZ)+WTRT2s1(N,L,NR,NZ)
        WTRTDs1(N,L,NZ)=WTRTDs1(N,L,NZ)+WTRT2s1(N,L,NR,NZ) &
          +WTRT1s1(N,L,NR,NZ)
5460  CONTINUE
      TCO2Ts1(NZ)=TCO2Ts1(NZ)+RCO2As1(N,L,NZ)
      RECOs1=RECOs1+RCO2As1(N,L,NZ)
      TRAUs1=TRAUs1+RCO2As1(N,L,NZ)
5450  CONTINUE
    DO 5470 NR=1,NRTs1(NZ)
      WTRTLs1(N,NINRs1(NR,NZ),NZ)=WTRTLs1(N,NINRs1(NR,NZ),NZ) &
        +RTWT1s1(N,NR,NZ)
5470  CONTINUE
5445  CONTINUE
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
!
!     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
!     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     WTLS,WTRT=total PFT leaf+petiole,root C mass
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     RLNT,RTNT=root layer,root system sink strength
!
!     IF(ISTYPs1(NZ).EQ.1)THEN
  IF(WTLSs1(NZ).GT.ZEROPs1(NZ))THEN
    FWTC=AMIN1(1.0,0.667*WTRTs1(NZ)/WTLSs1(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF
  IF(WTRTs1(NZ).GT.ZEROPs1(NZ))THEN
    FWTS=AMIN1(1.0,WTLSs1(NZ)/(0.667*WTRTs1(NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
  DO 290 L=NUs1,NIs1(NZ)
    IF(RTNT(1).GT.ZEROPs1(NZ))THEN
      FWTR(L)=AMAX1(0.0,RLNT(1,L)/RTNT(1))
    ELSE
      FWTR(L)=1.0_r8
    ENDIF
290 CONTINUE
!     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
!     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
!
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!
  WTLSs1(NZ)=0._r8
  DO 309 NB=1,NBRs1(NZ)
    WTLSs1(NZ)=WTLSs1(NZ)+WTLSBs1(NB,NZ)
309 CONTINUE
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAYs1(8,=end date for setting final seed number
!     FWTB=branch sink weighting factor
!     PTSHT=rate constant for equilibrating shoot-root nonstructural C concn from PFT file
!     PTRT=allocation to leaf+petiole used to modify PTSHT in annuals
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     FWODB=C woody fraction in branch:0=woody,1=non-woody
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  DO 310 NB=1,NBRs1(NZ)
    IF(IDTHBs1(NB,NZ).EQ.0)THEN
      IF(WTLSs1(NZ).GT.ZEROPs1(NZ))THEN
        FWTB(NB)=AMAX1(0.0,WTLSBs1(NB,NZ)/WTLSs1(NZ))
      ELSE
        FWTB(NB)=1.0_r8
      ENDIF
      IF(ISTYPs1(NZ).EQ.0)THEN
        PTSHTR=PTSHTs1(NZ)*PTRT**0.167
      ELSE
        PTSHTR=PTSHTs1(NZ)
      ENDIF
      DO 415 L=NUs1,NIs1(NZ)
        WTLSBX=WTLSBs1(NB,NZ)*FWODBs1(1)*FWTR(L)*FWTC
        WTRTLX=WTRTLs1(1,L,NZ)*FWODRs1(1)*FWTB(NB)*FWTS
        WTLSBB=AMAX1(0.0,WTLSBX,FSNK*WTRTLX)
        WTRTLR=AMAX1(0.0,WTRTLX,FSNK*WTLSBX)
        WTPLTT=WTLSBB+WTRTLR
        IF(WTPLTT.GT.ZEROPs1(NZ))THEN
          CPOOLB=AMAX1(0.0,CPOOLs1(NB,NZ)*FWTR(L))
          CPOOLS=AMAX1(0.0,CPOOLRs1(1,L,NZ)*FWTB(NB))
          CPOOLD=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
          XFRC=PTSHTR*CPOOLD
          CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)-XFRC
          CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)+XFRC
          CPOOLT=CPOOLS+CPOOLB
          IF(CPOOLT.GT.ZEROPs1(NZ))THEN
            ZPOOLB=AMAX1(0.0,ZPOOLs1(NB,NZ)*FWTR(L))
            ZPOOLS=AMAX1(0.0,ZPOOLRs1(1,L,NZ)*FWTB(NB))
            ZPOOLD=(ZPOOLB*CPOOLS-ZPOOLS*CPOOLB)/CPOOLT
            XFRN=PTSHTR*ZPOOLD
            PPOOLB=AMAX1(0.0,PPOOLs1(NB,NZ)*FWTR(L))
            PPOOLS=AMAX1(0.0,PPOOLRs1(1,L,NZ)*FWTB(NB))
            PPOOLD=(PPOOLB*CPOOLS-PPOOLS*CPOOLB)/CPOOLT
            XFRP=PTSHTR*PPOOLD
          ELSE
            XFRN=0._r8
            XFRP=0._r8
          ENDIF
          ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-XFRN
          ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)+XFRN
          PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-XFRP
          PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)+XFRP
!     IF(NZ.EQ.2.AND.NB.EQ.1)THEN
!     WRITE(*,3344)'ROOT',I,J,NX,NY,NZ,NB,L
!    2,XFRC,XFRN,XFRP,CPOOLs1(NB,NZ)
!    3,CPOOLRs1(1,L,NZ),ZPOOLs1(NB,NZ)
!    3,ZPOOLRs1(1,L,NZ),FWTB(NB),FWTR(L)
!    3,FWTC,FWTS,WTLSBX,WTRTLX,FSNK,FDBKs1(NB,NZ)
!    4,CPOOLD,CPOOLB,WTLSBB,CPOOLS,WTRTLR
!    5,FWOODs1(1),FWODBs1(1),WTRTLs1(1,L,NZ)
!    6,WTLSBs1(NB,NZ),RLNT(1,L),RTNT(1)
!3344  FORMAT(A8,7I4,30E12.4)
!     ENDIF
        ENDIF
415   CONTINUE
    ENDIF
310   CONTINUE
  end subroutine NonstructlBiomTransfer
!------------------------------------------------------------------------------------------

  subroutine ComputeTotalBiom(NZ,CPOOLK)

  integer, intent(in) :: NZ
  real(r8), intent(out) :: CPOOLK(JC1,JP1)
  integer :: L,K,N,NB
!     begin_execution
!     TOTAL C,N,P IN EACH BRANCH
!
!     CPOOLK=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     CPOOL,ZPOOL,PPOOL=C3 non-structural C,N,P mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
!
  DO 320 NB=1,NBRs1(NZ)
    CPOOLK(NB,NZ)=0._r8
    DO 325 K=1,JNODS1
      CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
        +CPOOL3s1(K,NB,NZ)+CPOOL4s1(K,NB,NZ) &
        +CO2Bs1(K,NB,NZ)+HCOBs1(K,NB,NZ)
325   CONTINUE
    WTSHTBs1(NB,NZ)=WTLFBs1(NB,NZ) &
      +WTSHEBs1(NB,NZ)+WTSTKBs1(NB,NZ)+WTRSVBs1(NB,NZ) &
      +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ)+WTGRBs1(NB,NZ) &
      +CPOOLs1(NB,NZ)+CPOOLK(NB,NZ)
    WTSHTNs1(NB,NZ)=WTLFBNs1(NB,NZ) &
      +WTSHBNs1(NB,NZ)+WTSTBNs1(NB,NZ)+WTRSBNs1(NB,NZ) &
      +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ)+WTGRBNs1(NB,NZ) &
      +ZPOOLs1(NB,NZ)
    WTSHTPs1(NB,NZ)=WTLFBPs1(NB,NZ) &
      +WTSHBPs1(NB,NZ)+WTSTBPs1(NB,NZ)+WTRSBPs1(NB,NZ) &
      +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ)+WTGRBPs1(NB,NZ) &
      +PPOOLs1(NB,NZ)
320   CONTINUE
!
!     TOTAL C,N,P IN ROOTS AND MYCORRHIZAE IN EACH SOIL LAYER
!
!     WTRTD=root C mass
!     CPOOLR=non-structural C mass in root
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
!     UPNF=PFT N2 fixation
!
  DO 345 N=1,MYs1(NZ)
    DO  L=NUs1,NIs1(NZ)
      WTRTDs1(N,L,NZ)=WTRTDs1(N,L,NZ)+CPOOLRs1(N,L,NZ)
    enddo
345   CONTINUE
  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,UPNFC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: UPNFC(JP1)
  integer :: L,NR,N,NB
!     begin_execution

!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     ARLFB=branch leaf area
!     ARSTK=total branch stalk surface area in each layer
!     GRNOB=seed set number
!
  CPOOLPs1(NZ)=sum(CPOOLs1(1:NBRs1(NZ),NZ))
  ZPOOLPs1(NZ)=sum(ZPOOLs1(1:NBRs1(NZ),NZ))
  PPOOLPs1(NZ)=sum(PPOOLs1(1:NBRs1(NZ),NZ))
  WTSHTs1(NZ)=sum(WTSHTBs1(1:NBRs1(NZ),NZ))
  WTSHNs1(NZ)=sum(WTSHTNs1(1:NBRs1(NZ),NZ))
  WTSHPs1(NZ)=sum(WTSHTPs1(1:NBRs1(NZ),NZ))
  WTLFs1(NZ)=sum(WTLFBs1(1:NBRs1(NZ),NZ))
  WTSHEs1(NZ)=sum(WTSHEBs1(1:NBRs1(NZ),NZ))
  WTSTKs1(NZ)=sum(WTSTKBs1(1:NBRs1(NZ),NZ))
  WVSTKs1(NZ)=sum(WVSTKBs1(1:NBRs1(NZ),NZ))
  WTRSVs1(NZ)=sum(WTRSVBs1(1:NBRs1(NZ),NZ))
  WTHSKs1(NZ)=sum(WTHSKBs1(1:NBRs1(NZ),NZ))
  WTEARs1(NZ)=sum(WTEARBs1(1:NBRs1(NZ),NZ))
  WTGRs1(NZ)=sum(WTGRBs1(1:NBRs1(NZ),NZ))
  WTLSs1(NZ)=sum(WTLSBs1(1:NBRs1(NZ),NZ))
  WTLFNs1(NZ)=sum(WTLFBNs1(1:NBRs1(NZ),NZ))
  WTSHENs1(NZ)=sum(WTSHBNs1(1:NBRs1(NZ),NZ))
  WTSTKNs1(NZ)=sum(WTSTBNs1(1:NBRs1(NZ),NZ))
  WTRSVNs1(NZ)=sum(WTRSBNs1(1:NBRs1(NZ),NZ))
  WTHSKNs1(NZ)=sum(WTHSBNs1(1:NBRs1(NZ),NZ))
  WTEARNs1(NZ)=sum(WTEABNs1(1:NBRs1(NZ),NZ))
  WTGRNNs1(NZ)=sum(WTGRBNs1(1:NBRs1(NZ),NZ))
  WTLFPs1(NZ)=sum(WTLFBPs1(1:NBRs1(NZ),NZ))
  WTSHEPs1(NZ)=sum(WTSHBPs1(1:NBRs1(NZ),NZ))
  WTSTKPs1(NZ)=sum(WTSTBPs1(1:NBRs1(NZ),NZ))
  WTRSVPs1(NZ)=sum(WTRSBPs1(1:NBRs1(NZ),NZ))
  WTHSKPs1(NZ)=sum(WTHSBPs1(1:NBRs1(NZ),NZ))
  WTEARPs1(NZ)=sum(WTEABPs1(1:NBRs1(NZ),NZ))
  WTGRNPs1(NZ)=sum(WTGRBPs1(1:NBRs1(NZ),NZ))
  GRNOs1(NZ)=sum(GRNOBs1(1:NBRs1(NZ),NZ))
  ARLFPs1(NZ)=sum(ARLFBs1(1:NBRs1(NZ),NZ))
  ARSTPs1(NZ)=sum(ARSTKs1(1:JC1,1:NBRs1(NZ),NZ))
  ARSTVs1(L,NZ)=sum(ARSTKs1(1:JC1,1:NBRs1(NZ),NZ))
!
!     ACCUMULATE ROOT STATE VARIABLES FROM ROOT LAYER STATE VARIABLES
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!

  WTRTs1(NZ)=sum(CPOOLRs1(1:MYs1(NZ),NUs1:NJs1,NZ))
  WTRTNs1(NZ)=sum(ZPOOLRs1(1:MYs1(NZ),NUs1:NJs1,NZ))
  WTRTPs1(NZ)=sum(PPOOLRs1(1:MYs1(NZ),NUs1:NJs1,NZ))
  WTRTSs1(NZ)=sum(WTRT1s1(1:MYs1(NZ),NUs1:NJs1,1:NRTs1(NZ),NZ)) &
    +sum(WTRT2s1(1:MYs1(NZ),NUs1:NJs1,1:NRTs1(NZ),NZ))
  WTRTSNs1(NZ)=sum(WTRT1Ns1(1:MYs1(NZ),NUs1:NJs1,1:NRTs1(NZ),NZ)) &
    +sum(WTRT2Ns1(1:MYs1(NZ),NUs1:NJs1,1:NRTs1(NZ),NZ))
  WTRTSPs1(NZ)=sum(WTRT1Ps1(1:MYs1(NZ),NUs1:NJs1,1:NRTs1(NZ),NZ)) &
    +sum(WTRT2Ps1(1:MYs1(NZ),NUs1:NJs1,1:NRTs1(NZ),NZ))
  WTRTs1(NZ)=WTRTs1(NZ)+WTRTSs1(NZ)
  WTRTNs1(NZ)=WTRTNs1(NZ)+WTRTSNs1(NZ)
  WTRTPs1(NZ)=WTRTPs1(NZ)+WTRTSPs1(NZ)
!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  IF(INTYPs1(NZ).NE.0)THEN
    IF(INTYPs1(NZ).GE.4)THEN
      DO 7950 NB=1,NBRs1(NZ)
        CPOLNPs1(NZ)=CPOLNPs1(NZ)+CPOLNBs1(NB,NZ)
        ZPOLNPs1(NZ)=ZPOLNPs1(NZ)+ZPOLNBs1(NB,NZ)
        PPOLNPs1(NZ)=PPOLNPs1(NZ)+PPOLNBs1(NB,NZ)
7950  CONTINUE
      WTNDs1(NZ)=sum(WTNDBs1(1:NBRs1(NZ),NZ))+&
        sum(CPOLNBs1(1:NBRs1(NZ),NZ))
      WTNDNs1(NZ)=sum(WTNDBNs1(1:NBRs1(NZ),NZ))+&
        sum(ZPOLNBs1(1:NBRs1(NZ),NZ))
      WTNDPs1(NZ)=sum(WTNDBPs1(1:NBRs1(NZ),NZ))+&
        sum(PPOLNBs1(1:NBRs1(NZ),NZ))
    ELSEIF(INTYPs1(NZ).GE.1.AND.INTYPs1(NZ).LE.3)THEN
      WTNDs1(NZ)=sum(WTNDLs1(NUs1:NIs1(NZ),NZ))+&
        sum(CPOOLNs1(NUs1:NIs1(NZ),NZ))
      WTNDNs1(NZ)=sum(WTNDLNs1(NUs1:NIs1(NZ),NZ))+&
        sum(ZPOOLNs1(NUs1:NIs1(NZ),NZ))
      WTNDPs1(NZ)=sum(WTNDLPs1(NUs1:NIs1(NZ),NZ))+&
        sum(PPOOLNs1(NUs1:NIs1(NZ),NZ))
    ENDIF
  ENDIF
!
!     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
!
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
!     UPNF=PFT N2 fixation
!     TCUPTK,TZUPTK,TPUPTK=cumulative PFT root-soil C,N,P exchange
!     TZUPFX=cumulative PFT N2 fixation
!
  HCUPTKs1(NZ)=UPOMCs1(NZ)
  HZUPTKs1(NZ)=UPOMNs1(NZ)+UPNH4s1(NZ)+UPNO3s1(NZ)+UPNFs1(NZ)
  HPUPTKs1(NZ)=UPOMPs1(NZ)+UPH2Ps1(NZ)+UPH1Ps1(NZ)
  TCUPTKs1(NZ)=TCUPTKs1(NZ)+UPOMCs1(NZ)
  TZUPTKs1(NZ)=TZUPTKs1(NZ)+UPOMNs1(NZ)+UPNH4s1(NZ)+UPNO3s1(NZ)
  TPUPTKs1(NZ)=TPUPTKs1(NZ)+UPOMPs1(NZ)+UPH2Ps1(NZ)+UPH1Ps1(NZ)
  TZUPFXs1(NZ)=TZUPFXs1(NZ)+UPNFs1(NZ)+UPNFC(NZ)
  end subroutine AccumulateStates

!------------------------------------------------------------------------------------------

  subroutine ResetBranchRootStates(NZ,CPOOLK)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: CPOOLK(JC1,JP1)
  integer :: L,NR,N,NB
!     begin_execution
!     RESET BRANCH STATE VARIABLES
!
  DO 8835 NB=1,NBRs1(NZ)
    CPOOLs1(NB,NZ)=0._r8
    CPOOLK(NB,NZ)=0._r8
    ZPOOLs1(NB,NZ)=0._r8
    PPOOLs1(NB,NZ)=0._r8
    CPOLNBs1(NB,NZ)=0._r8
    ZPOLNBs1(NB,NZ)=0._r8
    PPOLNBs1(NB,NZ)=0._r8
    WTSHTBs1(NB,NZ)=0._r8
    WTLFBs1(NB,NZ)=0._r8
    WTNDBs1(NB,NZ)=0._r8
    WTSHEBs1(NB,NZ)=0._r8
    WTSTKBs1(NB,NZ)=0._r8
    WVSTKBs1(NB,NZ)=0._r8
    WTRSVBs1(NB,NZ)=0._r8
    WTHSKBs1(NB,NZ)=0._r8
    WTEARBs1(NB,NZ)=0._r8
    WTGRBs1(NB,NZ)=0._r8
    WTLSBs1(NB,NZ)=0._r8
    WTSHTNs1(NB,NZ)=0._r8
    WTLFBNs1(NB,NZ)=0._r8
    WTNDBNs1(NB,NZ)=0._r8
    WTSHBNs1(NB,NZ)=0._r8
    WTSTBNs1(NB,NZ)=0._r8
    WTRSBNs1(NB,NZ)=0._r8
    WTHSBNs1(NB,NZ)=0._r8
    WTEABNs1(NB,NZ)=0._r8
    WTGRBNs1(NB,NZ)=0._r8
    WTSHTPs1(NB,NZ)=0._r8
    WTLFBPs1(NB,NZ)=0._r8
    WTNDBPs1(NB,NZ)=0._r8
    WTSHBPs1(NB,NZ)=0._r8
    WTSTBPs1(NB,NZ)=0._r8
    WTRSBPs1(NB,NZ)=0._r8
    WTHSBPs1(NB,NZ)=0._r8
    WTEABPs1(NB,NZ)=0._r8
    WTGRBPs1(NB,NZ)=0._r8
    WTSTXBs1(NB,NZ)=0._r8
    WTSTXNs1(NB,NZ)=0._r8
    WTSTXPs1(NB,NZ)=0._r8
8835  CONTINUE
!
!     RESET ROOT STATE VARIABLES
!
  DO 6416 L=NUs1,NJs1
    DO  N=1,MYs1(NZ)
      CPOOLRs1(N,L,NZ)=0._r8
      ZPOOLRs1(N,L,NZ)=0._r8
      PPOOLRs1(N,L,NZ)=0._r8
      DO  NR=1,NRTs1(NZ)
        WTRT1s1(N,L,NR,NZ)=0._r8
        WTRT1Ns1(N,L,NR,NZ)=0._r8
        WTRT1Ps1(N,L,NR,NZ)=0._r8
        WTRT2s1(N,L,NR,NZ)=0._r8
        WTRT2Ns1(N,L,NR,NZ)=0._r8
        WTRT2Ps1(N,L,NR,NZ)=0._r8
        RTWT1s1(N,NR,NZ)=0._r8
        RTWT1Ns1(N,NR,NZ)=0._r8
        RTWT1Ps1(N,NR,NZ)=0._r8
        RTLG1s1(N,L,NR,NZ)=0._r8
        RTLG2s1(N,L,NR,NZ)=0._r8
        RTN2s1(N,L,NR,NZ)=0._r8
      enddo
    enddo
6416  CONTINUE
  WTRVCs1(NZ)=0._r8
  WTRVNs1(NZ)=0._r8
  WTRVPs1(NZ)=0._r8
  IDTHs1(NZ)=1
  end subroutine ResetBranchRootStates
!------------------------------------------------------------------------------------------

  subroutine ResetDeadRootStates(NB,NZ,CPOOLK)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: CPOOLK(JC1,JP1)
  integer :: L,K,N
!     begin_execution
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     GRNOB=seed set number
!     GRNXB=potential number of seed set sites
!     GRWTB=individual seed size
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WSLF=leaf protein mass
!     ARLFB=branch leaf area
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
  CPOOLs1(NB,NZ)=0._r8
  CPOOLK(NB,NZ)=0._r8
  ZPOOLs1(NB,NZ)=0._r8
  PPOOLs1(NB,NZ)=0._r8
  CPOLNBs1(NB,NZ)=0._r8
  ZPOLNBs1(NB,NZ)=0._r8
  PPOLNBs1(NB,NZ)=0._r8
  WTSHTBs1(NB,NZ)=0._r8
  WTLFBs1(NB,NZ)=0._r8
  WTNDBs1(NB,NZ)=0._r8
  WTSHEBs1(NB,NZ)=0._r8
  WTSTKBs1(NB,NZ)=0._r8
  WVSTKBs1(NB,NZ)=0._r8
  WTRSVBs1(NB,NZ)=0._r8
  WTHSKBs1(NB,NZ)=0._r8
  WTEARBs1(NB,NZ)=0._r8
  WTGRBs1(NB,NZ)=0._r8
  WTLSBs1(NB,NZ)=0._r8
  WTSHTNs1(NB,NZ)=0._r8
  WTLFBNs1(NB,NZ)=0._r8
  WTNDBNs1(NB,NZ)=0._r8
  WTSHBNs1(NB,NZ)=0._r8
  WTSTBNs1(NB,NZ)=0._r8
  WTRSBNs1(NB,NZ)=0._r8
  WTHSBNs1(NB,NZ)=0._r8
  WTEABNs1(NB,NZ)=0._r8
  WTGRBNs1(NB,NZ)=0._r8
  WTSHTPs1(NB,NZ)=0._r8
  WTLFBPs1(NB,NZ)=0._r8
  WTNDBPs1(NB,NZ)=0._r8
  WTSHBPs1(NB,NZ)=0._r8
  WTSTBPs1(NB,NZ)=0._r8
  WTRSBPs1(NB,NZ)=0._r8
  WTHSBPs1(NB,NZ)=0._r8
  WTEABPs1(NB,NZ)=0._r8
  WTGRBPs1(NB,NZ)=0._r8
  GRNXBs1(NB,NZ)=0._r8
  GRNOBs1(NB,NZ)=0._r8
  GRWTBs1(NB,NZ)=0._r8
  ARLFBs1(NB,NZ)=0._r8
  WTSTXBs1(NB,NZ)=0._r8
  WTSTXNs1(NB,NZ)=0._r8
  WTSTXPs1(NB,NZ)=0._r8
  DO 8855 K=0,JNODS1
    IF(K.NE.0)THEN
      CPOOL3s1(K,NB,NZ)=0._r8
      CPOOL4s1(K,NB,NZ)=0._r8
      CO2Bs1(K,NB,NZ)=0._r8
      HCOBs1(K,NB,NZ)=0._r8
    ENDIF
    ARLF1s1(K,NB,NZ)=0._r8
    HTNODEs1(K,NB,NZ)=0._r8
    HTNODXs1(K,NB,NZ)=0._r8
    HTSHEs1(K,NB,NZ)=0._r8
    WGLFs1(K,NB,NZ)=0._r8
    WSLFs1(K,NB,NZ)=0._r8
    WGLFNs1(K,NB,NZ)=0._r8
    WGLFPs1(K,NB,NZ)=0._r8
    WGSHEs1(K,NB,NZ)=0._r8
    WSSHEs1(K,NB,NZ)=0._r8
    WGSHNs1(K,NB,NZ)=0._r8
    WGSHPs1(K,NB,NZ)=0._r8
    WGNODEs1(K,NB,NZ)=0._r8
    WGNODNs1(K,NB,NZ)=0._r8
    WGNODPs1(K,NB,NZ)=0._r8
    DO 8865 L=1,JC1
      ARLFVs1(L,NZ)=ARLFVs1(L,NZ)-ARLFLs1(L,K,NB,NZ)
      WGLFVs1(L,NZ)=WGLFVs1(L,NZ)-WGLFLs1(L,K,NB,NZ)
      ARLFLs1(L,K,NB,NZ)=0._r8
      WGLFLs1(L,K,NB,NZ)=0._r8
      WGLFLNs1(L,K,NB,NZ)=0._r8
      WGLFLPs1(L,K,NB,NZ)=0._r8
      IF(K.NE.0)THEN
        DO 8860 N=1,JLI1
          SURFs1(N,L,K,NB,NZ)=0._r8
8860    CONTINUE
      ENDIF
8865  CONTINUE
8855  CONTINUE
  DO 8875 L=1,JC1
    ARSTKs1(L,NB,NZ)=0._r8
    DO  N=1,JLI1
      SURFBs1(N,L,NB,NZ)=0._r8
    enddo
8875  CONTINUE

  end subroutine ResetDeadRootStates
!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C4(K,NB,NZ,WFNG,WFNC,CH2O3,CH2O4,CO2F,CH2O)
  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(in):: WFNG,WFNC
  real(r8), intent(inout) :: CH2O3(25),CH2O4(25),CO2F,CH2O
  integer :: L,NN,M,N
  real(r8) :: WFN4
  real(r8) :: WFNB
  real(r8) :: CO2X,CO2C,CO2Y
  real(r8) :: CBXNX
  real(r8) :: DIFF
  real(r8) :: ETLF4
  real(r8) :: EGRO4
  real(r8) :: EGROX
  real(r8) :: ETLF
  real(r8) :: EGRO
  real(r8) :: GSL
  real(r8) :: PARX,PARJ
  real(r8) :: RS,RSL,VL
  real(r8) :: VGROX
  real(r8) :: VA,VG
! begin_execution

! FOR EACH CANOPY LAYER
!
  DO 110 L=JC1,1,-1
    IF(ARLFLs1(L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      DO 115 N =1,JLI1
        DO 120 M =1,JSA1
!
!         CO2 FIXATION IN MESOPHYLL BY SUNLIT LEAVES
!
!         SURFX=unself-shaded leaf surface area
!
          IF(SURFXs1(N,L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
            IF(PARs1(N,M,L,NZ).GT.0.0)THEN
!
!             C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             ETGR4=light saturated e- transport rate from stomate.f
!             ETLF4=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO4=light-limited PEP carboxylation rate
!             CBXN4=PEP caboxylation efficiency
!             VL=PEP carboxylation rate limited by light,CO2,N,P
!             VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
!             FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARs1(N,M,L,NZ)
              PARJ=PARX+ETGR4s1(K,NB,NZ)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4s1(K,NB,NZ)))/CURV2
              EGRO4=ETLF4*CBXN4s1(K,NB,NZ)
              VL=AMIN1(VGRO4s1(K,NB,NZ),EGRO4)*FDBK4s1(K,NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFN4=RS/RSL
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP carboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               VCGR4,VGROX=PEP carboxylation rate unlimited,limited by CO2
!               XKCO24=Km for VCMX4 from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 125 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
                  VGROX=VCGR4s1(K,NB,NZ)*CO2Y/(CO2C+XKCO24s1(NZ))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4s1(K,NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
125             CONTINUE
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O4(K)=CH2O4(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAUSs1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1))*0.0432
!               IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.3.AND.K.EQ.KLEAFs1(NB,NZ)
!              2.AND.(I/10)*10.EQ.I.AND.J.EQ.12)THEN
!               WRITE(*,4444)'VLD4',IYRCs1,I,J,NZ,L,M,N,K,VL,PARs1(N,M,L,NZ)
!              2,PARs1(N,M,L,NZ)*TAUSs1(L+1)+PARDIFs1(N,M,L,NZ)
!              3*TAU0s1(L+1)
!              2,RAPS,TKCs1(NZ),CO2Qs1(NZ),ETGR4s1(K,NB,NZ)
!              3,CBXN4s1(K,NB,NZ),VGRO4s1(K,NB,NZ),EGRO
!              3,FDBK4s1(K,NB,NZ),CH2O4(K),WFN4,VGROX,EGROX
!              4,VCGR4s1(K,NB,NZ),CO2X,CO2C,CBXNX
!              5,RS,RSL,SURFXs1(N,L,K,NB,NZ)
!4444            FORMAT(A8,8I4,40E12.4)
!               ENDIF
!
!               C3 CARBOXYLATION REACTIONS IN BUNDLE SHEATH OF C4 PLANTS
!
!               ETGRO=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               CBXN=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!               FDBK=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+ETGROs1(K,NB,NZ)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
                EGRO=ETLF*CBXNs1(K,NB,NZ)
                VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*WFNB*FDBKs1(NB,NZ)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAUSs1(L+1)
!               IF(L.EQ.NC-1.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
!               WRITE(*,4445)'VLD3',IYRCs1,I,J,NZ,L,M,N,K,VL,PARs1(N,M,L,NZ)
!              2,RAPS,TKCs1(NZ),CO2Qs1(NZ),ETGROs1(K,NB,NZ)
!              3,CBXNs1(K,NB,NZ),VGROs1(K,NB,NZ),EGRO
!              3,FDBKs1(NB,NZ),WFNB
!4445            FORMAT(A8,8I4,20E12.4)
!               ENDIF
              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDIFs1(N,M,L,NZ).GT.0.0)THEN
!
!           C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!           QNTM=quantum efficiency
!           PARDIFs1=diffuse PAR flux
!           ETGR4=light saturated e- transport rate from stomate.f
!           ETLF4=light-limited e- transport rate
!           CURV=shape parameter for e- transport response to PAR
!           EGRO4=light-limited PEP carboxylation rate
!           CBXN4=PEP caboxylation efficiency
!           VL=PEP carboxylation rate limited by light,CO2,N,P
!           VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
!           FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARDIFs1(N,M,L,NZ)
              PARJ=PARX+ETGR4s1(K,NB,NZ)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4s1(K,NB,NZ)))/CURV2
              EGRO4=ETLF4*CBXN4s1(K,NB,NZ)
              VL=AMIN1(VGRO4s1(K,NB,NZ),EGRO4)*FDBK4s1(K,NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFN4,WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFN4=(RS/RSL)**1.00
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP caboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
!               XKCO24=Km for VCMX4 from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 135 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
                  VGROX=VCGR4s1(K,NB,NZ)*CO2Y/(CO2C+XKCO24s1(NZ))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4s1(K,NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
135             CONTINUE
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0s1=fraction of diffuse radiation transmitted from layer above
!
                CH2O4(K)=CH2O4(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAU0s1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAU0s1(L+1))*0.0432
!               WRITE(*,4455)'VLB4',IYRCs1,I,J,NZ,L,M,N,K,VL,PARs1(N,M,L,NZ)
!              2,RAPS,TKCs1(NZ),CO2Qs1(NZ),ETGR4s1(K,NB,NZ)
!              3,CBXN4s1(K,NB,NZ),VGRO4s1(K,NB,NZ),EGRO
!              3,FDBK4s1(K,NB,NZ),CH2O4(K),WFN4,VGROX,EGROX
!              4,VCGR4s1(K,NB,NZ),CO2X,CO2C,CBXNX
!              5,RS,RSL,SURFXs1(N,L,K,NB,NZ)
!4455            FORMAT(A8,8I4,40E12.4)
!
!               C3 CARBOXYLATION REACTIONS IN IN BUNDLE SHEATH OF C4 PLANTS
!
!               ETGRO=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               CBXN=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!               FDBK=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+ETGROs1(K,NB,NZ)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
                EGRO=ETLF*CBXNs1(K,NB,NZ)
                VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*WFNB*FDBKs1(NB,NZ)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0s1=fraction of diffuse radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAU0s1(L+1)
!               IF(J.EQ.13.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
!               WRITE(*,4444)'VLB4',IYRCs1,I,J,NZ,L,K,VL,PARDIFs1(N,M,L,NZ)
!              2,RAPY,TKCs1(NZ),CO2Qs1(NZ),CO2X,FMOLs1(NZ)/GSL
!              3,VCGROs1(K,NB,NZ),ETLF,FDBKs1(NB,NZ),WFNB
!               ENDIF
              ENDIF
            ENDIF
          ENDIF
120     CONTINUE
115   CONTINUE
    ENDIF
110   CONTINUE
  CO2F=CO2F+CH2O4(K)
  CH2O=CH2O+CH2O3(K)
  end subroutine ComputeGPP_C4

!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C3(K,NB,NZ,WFNG,WFNC,CH2O3,CO2F,CH2O)

  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(in) :: WFNG,WFNC
  real(r8), intent(inout) :: CH2O3(25),CO2F,CH2O
  integer :: L,NN,M,N
  real(r8) :: WFNB
  real(r8) :: CO2X,CO2C,CO2Y
  real(r8) :: CBXNX
  real(r8) :: DIFF
  real(r8) :: EGROX
  real(r8) :: ETLF
  real(r8) :: EGRO
  real(r8) :: GSL
  real(r8) :: PARX,PARJ
  real(r8) :: RS,RSL
  real(r8) :: VL,VGROX
  real(r8) :: VA,VG
!
! FOR EACH CANOPY LAYER
!
  DO 210 L=JC1,1,-1
    IF(ARLFLs1(L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      DO 215 N=1,JLI1
        DO 220 M=1,JSA1
!
!         CO2 FIXATION BY SUNLIT LEAVES
!
!         SURFX=unself-shaded leaf surface area
!
          IF(SURFXs1(N,L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
            IF(PARs1(N,M,L,NZ).GT.0.0)THEN
!
!             C3 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             ETGRO=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             CBXN=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!             FDBK=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARs1(N,M,L,NZ)
              PARJ=PARX+ETGROs1(K,NB,NZ)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
              EGRO=ETLF*CBXNs1(K,NB,NZ)
              VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*FDBKs1(NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMPL=C3 CO2 compensation point (uM)
!               CBXNX=rubisco carboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco
!               VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               XKCO2O=Km for rubisco carboxylation
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 225 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMPLs1(K,NB,NZ))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPLs1(K,NB,NZ))
                  VGROX=VCGROs1(K,NB,NZ)*CO2Y/(CO2C+XKCO2Os1(NZ))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*FDBKs1(NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
225             CONTINUE
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAUSs1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1))*0.0432
!               IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1.AND.K.EQ.KLEAFs1(NB,NZ)-1
!              2.AND.J.EQ.12)THEN
!               WRITE(20,3335)'VLD',IYRCs1,I,J,NZ,L,M,N,K,VL,PARs1(N,M,L,NZ)
!              2,RAPS,TKCs1(NZ),TKA,CO2Qs1(NZ),CO2X,CO2C,FMOLs1(NZ)
!              3/GSL,VGROX,EGROX,ETLF,CBXNX,FDBKs1(NB,NZ),WFNB,PSILGs1(NZ)
!              4,VCGROs1(K,NB,NZ),ETGROs1(K,NB,NZ),COMPLs1(K,NB,NZ)
!              5,SURFXs1(N,L,K,NB,NZ),TAUSs1(L+1),CH2O3(K)
!3335            FORMAT(A8,8I4,30E12.4)
!               ENDIF
              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDIFs1(N,M,L,NZ).GT.0.0)THEN
!
!             C3 CARBOXYLATION REACTIONS USING VARIABLES FROM 'STOMATE'
!
!             QNTM=quantum efficiency
!             PARDIFs1=diffuse PAR flux
!             ETGRO=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             CBXN=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!             FDBK=N,P feedback inhibition on C3 CO2 fixation
!
              PARX=QNTM*PARDIFs1(N,M,L,NZ)
              PARJ=PARX+ETGROs1(K,NB,NZ)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
              EGRO=ETLF*CBXNs1(K,NB,NZ)
              VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*FDBKs1(NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C3 CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMPL=C3 CO2 compensation point (uM)
!               CBXNX=rubisco caboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco carboxylase
!               VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               XKCO2O=Km for rubisco carboxylation from stomate.f (uM)
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 235 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMPLs1(K,NB,NZ))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPLs1(K,NB,NZ))
                  VGROX=VCGROs1(K,NB,NZ)*CO2Y/(CO2C+XKCO2Os1(NZ))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*FDBKs1(NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
235             CONTINUE
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0s1=fraction of diffuse radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAU0s1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAU0s1(L+1))*0.0432
!               IF(J.EQ.13.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
!               WRITE(*,3335)'VLB',IYRCs1,I,J,NZ,L,K,VL,PARDIFs1(N,M,L,NZ)
!              2,RAPY,TKCs1(NZ),CO2Qs1(NZ),CO2X,FMOLs1(NZ)/GSL
!              3,VCGROs1(K,NB,NZ),ETLF,FDBKs1(NB,NZ),WFNB
!               ENDIF
              ENDIF
            ENDIF
          ENDIF
220     CONTINUE
215   CONTINUE
    ENDIF
210   CONTINUE
  CO2F=CO2F+CH2O3(K)
  CH2O=CH2O+CH2O3(K)
  end subroutine ComputeGPP_C3

!------------------------------------------------------------------------------------------

  subroutine ComputeRAutoAfEmergence(NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,CO2F,&
    CH2O,TFN5,WFNG,WFNSG,WTSHXN,ZADDB,CNPG,PADDB,RCO2C,RMNCS,SNCR,CGROS,CNRDA)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(out) :: ZADDB
  real(r8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: PADDB,RCO2C,RMNCS,SNCR,CGROS,CNRDA
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,CO2F,CH2O,TFN5
  real(r8), intent(in) :: WFNG,WFNSG
  real(r8) :: ZPOOLB
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y,RCO2G
  real(r8) :: RCO2T,RCO2CM
! begin_execution
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
!
  IF(CCPOLBs1(NB,NZ).GT.ZEROs1)THEN
    CNPG=AMIN1(CZPOLBs1(NB,NZ)/(CZPOLBs1(NB,NZ) &
      +CCPOLBs1(NB,NZ)*CNKI),CPPOLBs1(NB,NZ)/(CPPOLBs1(NB,NZ) &
      +CCPOLBs1(NB,NZ)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P
!
! RCO2C=respiration from non-structural C
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! TFN3=temperature function for canopy growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! FDBKX=termination feedback inhibition on C3 CO2
!
  RCO2C=AMAX1(0.0,VMXC*CPOOLs1(NB,NZ) &
    *TFN3s1(NZ))*CNPG*FDBKXs1(NB,NZ)*WFNG
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RMNCS=maintenance respiration
! TFN5=temperature function for canopy maintenance respiration
! WTSHXN=shoot structural N mass
! IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WFNG=growth function of canopy water potential
!
  RMNCS=AMAX1(0.0,RMPLT*TFN5*WTSHXN)
  IF(IGTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.2)THEN
    RMNCS=RMNCS*WFNG
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X=difference between non-structural C respn and mntc respn
! RCO2Y=growth respiration unlimited by N,P
! WFNSG=expansion,extension function of canopy water potential
! SNCR=excess maintenance respiration
!
  RCO2X=RCO2C-RMNCS
  RCO2Y=AMAX1(0.0,RCO2X)*WFNSG
  SNCR=AMAX1(0.0,-RCO2X)
!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2Y,RCO2G=growth respiration unlimited,limited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! DMSHD=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
!
  IF(RCO2Y.GT.0.0.AND.(CNSHX.GT.0.0.OR.CNLFX.GT.0.0))THEN
    ZPOOLB=AMAX1(0.0,ZPOOLs1(NB,NZ))
    PPOOLB=AMAX1(0.0,PPOOLs1(NB,NZ))
    RCO2G=AMIN1(RCO2Y,ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
  ELSE
    RCO2G=0._r8
  ENDIF
!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! CGROS=total non-structural C used in growth and growth respiration
! RCO2G=growth respiration limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! ZADDB,PADDB=nonstructural N,P used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! CNRDA=respiration for N assimilation
! CH2O=total CH2O production
!
  CGROS=RCO2G/DMSHD
  ZADDB=AMAX1(0.0,AMIN1(ZPOOLs1(NB,NZ),CGROS*(CNSHX+CNLFM+CNLFX*CNPG)))
  PADDB=AMAX1(0.0,AMIN1(PPOOLs1(NB,NZ) &
    ,CGROS*(CPSHX+CPLFM+CPLFX*CNPG)))
  CNRDA=AMAX1(0.0,1.70*ZADDB-0.025*CH2O)
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! RCO2T=total C respiration
! RMNCS=maintenance respiration
! RCO2C=respiration from non-structural C
! RCO2G=growth respiration limited by N,P
! SNCR=excess maintenance respiration
! CNRDA=respiration for N assimilation
! CARBN=total PFT CO2 fixation
! CO2F=total CO2 fixation
! TCO2T,TCO2A=total,above-ground PFT respiration
! CNET=PFT net CO2 fixation
! TGPP=ecosystem GPP
! RECO=ecosystem respiration
! TRAU=total autotrophic respiration
!
  RCO2T=AMIN1(RMNCS,RCO2C)+RCO2G+SNCR+CNRDA
  CARBNs1(NZ)=CARBNs1(NZ)+CO2F
  TCO2Ts1(NZ)=TCO2Ts1(NZ)-RCO2T
  TCO2As1(NZ)=TCO2As1(NZ)-RCO2T
  CNETs1(NZ)=CNETs1(NZ)+CO2F-RCO2T
  TGPPs1=TGPPs1+CO2F
  RECOs1=RECOs1-RCO2T
  TRAUs1=TRAUs1-RCO2T
! IF(NZ.EQ.1)THEN
!   WRITE(*,4477)'RCO2',I,J,NX,NY,NZ,NB,IFLGZ,CPOOLs1(NB,NZ)
!    2,CH2O,RMNCS,RCO2C,CGROS,CNRDA,CNPG,RCO2T,RCO2X,SNCR
!    3,RCO2G,DMSHD,ZADDB,PART(1),PART(2),DMLFB,DMSHB
!    4,WTRSVBs1(NB,NZ),WVSTKBs1(NB,NZ),WTSHXN
!    5,ZPOOLs1(NB,NZ),PPOOLs1(NB,NZ),PSILTs1(NZ)
!    6,ZADDB,RNH3Bs1(NB,NZ),WFRs1(1,NGs1(NZ),NZ)
!    7,WFNG,TFN3s1(NZ),TFN5,FDBKXs1(NB,NZ),VMXC
!4477  FORMAT(A8,7I4,40E12.4)
! ENDIF
!
  end subroutine ComputeRAutoAfEmergence


!------------------------------------------------------------------------------------------

  subroutine ComputeRAutoBfEmergence(NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,&
    CNLFX,CPLFX,WTSHXN,WFNG,WFNSG,ZADDB,CNPG,PADDB,RCO2C,RMNCS,SNCR,CGROS,CNRDM,CNRDA,CH2O)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(in) :: TFN6(JZ1)
  real(r8), intent(out) :: ZADDB,rco2c
  real(r8), INTENT(OUT) :: CNPG,PADDB,RMNCS,SNCR,CGROS,CNRDM,CNRDA,CH2O
  real(r8), intent(in) :: DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,WFNG
  real(r8), intent(in) :: WFNSG
  real(r8) :: ZPOOLB,ZADDBM,CGROSM
  real(r8) :: FNP
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y,RCO2G
  real(r8) :: RCO2T,RCO2CM
  real(r8) :: RCO2XM,RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: SNCRM
! begin_execution
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  IF(CCPOLBs1(NB,NZ).GT.ZEROs1)THEN
    CNPG=AMIN1(CZPOLBs1(NB,NZ)/(CZPOLBs1(NB,NZ) &
      +CCPOLBs1(NB,NZ)*CNKI) &
      ,CPPOLBs1(NB,NZ)/(CPPOLBs1(NB,NZ) &
      +CCPOLBs1(NB,NZ)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P, O2 UPTAKE
!
! RCO2CM,RCO2C=respiration from non-structural C unlimited,limited by O2
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! TFN4=temperature function for root growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! FDBKX=termination feedback inhibition on C3 CO2
! WFR=constraint by O2 consumption on all root processes
!
  RCO2CM=AMAX1(0.0,VMXC*CPOOLs1(NB,NZ) &
    *TFN4s1(NGs1(NZ),NZ))*CNPG*WFNG*FDBKXs1(NB,NZ)
  RCO2C=RCO2CM*WFRs1(1,NGs1(NZ),NZ)
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RMNCS=maintenance respiration
! TFN6=temperature function for root maintenance respiration
! WTSHXN=shoot structural N mass
! IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WFNG=growth function of canopy water potential
!
  RMNCS=AMAX1(0.0,RMPLT*TFN6(NGs1(NZ))*WTSHXN)
  IF(IGTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.2)THEN
    RMNCS=RMNCS*WFNG
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! WFNSG=expansion,extension function of canopy water potential
! SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!
  RCO2XM=RCO2CM-RMNCS
  RCO2X=RCO2C-RMNCS
  RCO2YM=AMAX1(0.0,RCO2XM)*WFNSG
  RCO2Y=AMAX1(0.0,RCO2X)*WFNSG
  SNCRM=AMAX1(0.0,-RCO2XM)
  SNCR=AMAX1(0.0,-RCO2X)
!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! FNP=growth respiration limited by O2 and N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! RCO2GM,RCO2G=growth respiration unltd,ltd by O2 and limited by N,P
! WFR=constraint by O2 consumption on all root processes
!
  IF(CNSHX.GT.0.0.OR.CNLFX.GT.0.0)THEN
    ZPOOLB=AMAX1(0.0,ZPOOLs1(NB,NZ))
    PPOOLB=AMAX1(0.0,PPOOLs1(NB,NZ))
    FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
    IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
    ELSE
      RCO2GM=0._r8
    ENDIF
    IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFRs1(1,NGs1(NZ),NZ))
    ELSE
      RCO2G=0._r8
    ENDIF
  ELSE
    RCO2GM=0._r8
    RCO2G=0._r8
  ENDIF
!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! CGROSM,CGROS=total non-structural C used in growth and respn unltd,ltd by O2
! RCO2GM,RCO2G=growth respiration unltd,ltd by O2 and limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! ZADDBM,ZADDB,PADDB=nonstructural N,P unltd,ltd by O2 used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
  CGROSM=RCO2GM/DMSHD
  CGROS=RCO2G/DMSHD
  ZADDBM=AMAX1(0.0,CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
  ZADDB=AMAX1(0.0,CGROS*(CNSHX+CNLFM+CNLFX*CNPG))
  PADDB=AMAX1(0.0,CGROS*(CPSHX+CPLFM+CPLFX*CNPG))
  CNRDM=AMAX1(0.0,1.70*ZADDBM)
  CNRDA=AMAX1(0.0,1.70*ZADDB)
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! RCO2TM,RCO2T=total C respiration unltd,ltd by O2
! RMNCS=maintenance respiration
! RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
! SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
! CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
! RCO2A=total root respiration
! RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
  RCO2TM=RMNCS+RCO2GM+SNCRM+CNRDM
  RCO2T=RMNCS+RCO2G+SNCR+CNRDA
  RCO2Ms1(1,NGs1(NZ),NZ)=RCO2Ms1(1,NGs1(NZ),NZ)+RCO2TM
  RCO2Ns1(1,NGs1(NZ),NZ)=RCO2Ns1(1,NGs1(NZ),NZ)+RCO2T
  RCO2As1(1,NGs1(NZ),NZ)=RCO2As1(1,NGs1(NZ),NZ)-RCO2T
  CH2O=0._r8
  end subroutine ComputeRAutoBfEmergence

end module grosubsMod
