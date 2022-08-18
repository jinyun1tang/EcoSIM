module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : test_aeqb,safe_adb
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  use PhotoSynsMod
  use RootMod, only : RootBGCModel
  use NoduleBGCMod
  use LitterFallMod
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
  associate(                            &
    ZCs1       => plt_morph%ZCs1        &
  )
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
! TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS
  call LiveDeadTransformation(I,J)
  end associate
  END subroutine grosubs

!------------------------------------------------------------------------------------------

  subroutine LiveDeadTransformation(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: L,K,NZ,M,NB
  real(r8) :: XFRC,XFRN,XFRP
!     begin_execution

  associate(                           &
    IGTYPs1   => plt_pheno%IGTYPs1   , &
    IFLGIs1   => plt_pheno%IFLGIs1   , &
    IFLGEs1   => plt_pheno%IFLGEs1   , &
    IBTYPs1   => plt_pheno%IBTYPs1   , &
    VRNLs1    => plt_pheno%VRNLs1    , &
    VRNSs1    => plt_pheno%VRNSs1    , &
    SDPTHIs1  => plt_morph%SDPTHIs1  , &
    NBRs1     => plt_morph%NBRs1       &
  )
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
    ENDIF
9975  CONTINUE
  end associate
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
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT
  real(r8) :: XRTN1
  real(r8) :: TFN5
  real(r8) :: WFNG
  real(r8) :: WFNC
  real(r8) :: WFNS,WFNSG
! begin_execution
  associate(                              &
    IGTYPs1  => plt_pheno%IGTYPs1       , &
    IDTHRs1  => plt_pheno%IDTHRs1       , &
    IDTHPs1  => plt_pheno%IDTHPs1       , &
    SDARs1   => plt_morph%SDARs1        , &
    SDVLs1   => plt_morph%SDVLs1        , &
    NBRs1    => plt_morph%NBRs1         , &
    NRTs1    => plt_morph%NRTs1           &
  )
  IF(IDTHPs1(NZ).EQ.0.OR.IDTHRs1(NZ).EQ.0)THEN
    UPNFC(NZ)=0._r8
    IFLGZ = 0
    call StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,WFNC,WFNS,WFNSG)
!
!     CALCULATE GROWTH OF EACH BRANCH
!
!     WTLFB,WTSHEB,WTLSB=leaf,petiole,leaf+petiole mass
!     IDTHB=branch living flag: 0=alive,1=dead
!
    DO 105 NB=1,NBRs1(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6,ZCX,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WFNG,WFNC,WFNS,WFNSG,PTRT,UPNFC,IFLGZ)
105 CONTINUE
!
    call RootBGCModel(I,J,NZ,IFLGZ,ICHK1,IDTHRN,NRX,PTRT,TFN6,CNRTW,CPRTW,XRTN1)

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
  end associate
  end subroutine GrowPlant


!------------------------------------------------------------------------------------------

  subroutine C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: CH2O3(25),CH2O4(25)
  integer :: K
  real(r8) :: CPL4M,CCBS,CPL3K,CO2LK

!     begin_execution
  associate(                              &
    CO2Ls1    =>   plt_photo%CO2Ls1       &
  )
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
  end associate
  end subroutine C4PhotoProductTransfer

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

  associate(                              &
    VRNXs1    =>  plt_pheno%VRNXs1      , &
    FLGZs1    =>  plt_pheno%FLGZs1      , &
    IDAYs1    =>  plt_pheno%IDAYs1      , &
    IGTYPs1   =>  plt_pheno%IGTYPs1     , &
    TGSTGFs1  =>  plt_pheno%TGSTGFs1    , &
    IDTYPs1   =>  plt_pheno%IDTYPs1     , &
    IBTYPs1   =>  plt_pheno%IBTYPs1     , &
    VRNFs1    =>  plt_pheno%VRNFs1      , &
    TGSTGIs1  =>  plt_pheno%TGSTGIs1    , &
    ISTYPs1   =>   plt_pheno%ISTYPs1    , &
    IWTYPs1   =>  plt_pheno%IWTYPs1     , &
    SNL1s1    =>  plt_morph%SNL1s1      , &
    ARLFPs1   =>  plt_morph%ARLFPs1     , &
    NB1s1     =>  plt_morph%NB1s1         &
  )

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
      PART(1)=AMAX1(PART1X,(0.725-FPART1)*(1.0_r8-TGSTGFs1(NB,NZ)))
      PART(2)=AMAX1(PART2X,(0.275-FPART2)*(1.0_r8-TGSTGFs1(NB,NZ)))
    ENDIF
    PARTS=1.0_r8-PART(1)-PART(2)
    PART(3)=AMAX1(0.0,0.60*PARTS*(1.0_r8-TGSTGFs1(NB,NZ)))
    PART(4)=AMAX1(0.0,0.30*PARTS*(1.0_r8-TGSTGFs1(NB,NZ)))
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
    PART(3)=PART(3)+(1.0_r8-FPARTL)*(PART(1)+PART(2))
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
  end associate
  end subroutine CalcPartitionCoeff

!------------------------------------------------------------------------------------------

  subroutine UpdatePhotosynthates(I,J,NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
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
  associate(                            &
    IDAYs1     =>  plt_pheno%IDAYs1     &
  )
! begin_execution
! FDBK=N,P feedback inhibition on C3 CO2 fixation
! SSIN=sine of solar angle
! RADP=total PAR absorbed by canopy
! CO2Q=canopy air CO2 concentration
!
  IF(IDAYs1(1,NB,NZ).NE.0)THEN

    call ComputeGPP(NB,NZ,WFNG,WFNC,CH2O3,CH2O4,CH2O)
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
  end associate
  end subroutine UpdatePhotosynthates

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
  associate(                                &
    IDTHBs1      =>  plt_pheno%IDTHBs1    , &
    KVSTGNs1     =>  plt_pheno%KVSTGNs1   , &
    XRLAs1       =>  plt_pheno%XRLAs1     , &
    KVSTGs1      =>  plt_pheno%KVSTGs1    , &
    IGTYPs1      =>  plt_pheno%IGTYPs1    , &
    IFLGPs1      =>  plt_pheno%IFLGPs1    , &
    FLGZs1       =>  plt_pheno%FLGZs1     , &
    IDAYs1       =>  plt_pheno%IDAYs1     , &
    IBTYPs1      =>  plt_pheno%IBTYPs1    , &
    IFLGGs1      =>  plt_pheno%IFLGGs1    , &
    ISTYPs1      =>  plt_pheno%ISTYPs1    , &
    SSINNs1      =>   plt_rad%SSINNs1     , &
    NB1s1        =>   plt_morph%NB1s1     , &
    NBTBs1       =>   plt_morph%NBTBs1    , &
    NBRs1        =>   plt_morph%NBRs1     , &
    SLA1s1       =>   plt_morph%SLA1s1    , &
    SNL1s1       =>   plt_morph%SNL1s1    , &
    ANGSHs1      =>   plt_morph%ANGSHs1   , &
    ARLFZs1      =>   plt_morph%ARLFZs1   , &
    ARLFBs1      =>   plt_morph%ARLFBs1   , &
    SDPTHs1      =>   plt_morph%SDPTHs1   , &
    HTSHEXs1     =>   plt_morph%HTSHEXs1  , &
    ARLF1s1      =>   plt_morph%ARLF1s1   , &
    HTSHEs1      =>   plt_morph%HTSHEs1   , &
    ANGBRs1      =>   plt_morph%ANGBRs1   , &
    SSL1s1       =>   plt_morph%SSL1s1    , &
    HTNODEs1     =>   plt_morph%HTNODEs1  , &
    HTNODXs1     =>   plt_morph%HTNODXs1  , &
    NNODs1       =>   plt_morph%NNODs1      &
  )
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
    call UpdatePhotosynthates(I,J,NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,&
      CNLFX,CPLFX,WTSHXN,TFN5,WFNG,WFNC,WFNSG,CH2O3,CH2O4,CNPG,rco2c,RMNCS,&
      SNCR,CGROS,CNRDM,CNRDA)
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
              HTNODEs1(K1,NB,NZ)=HTNODXs1(K1,NB,NZ)+HTNODEs1(K2,NB,NZ)
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
                RCZLXs1(NB,NZ)=WGLFNXs1(NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
                RCPLXs1(NB,NZ)=WGLFPXs1(NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
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
            ARLFBs1(NB,NZ)=ARLFBs1(NB,NZ)-FSNCL*ARLFZs1(NB,NZ)
            WTLFBs1(NB,NZ)=WTLFBs1(NB,NZ)-FSNCL*WGLFXs1(NB,NZ)
            WTLFBNs1(NB,NZ)=WTLFBNs1(NB,NZ)-FSNCL*WGLFNXs1(NB,NZ)
            WTLFBPs1(NB,NZ)=WTLFBPs1(NB,NZ)-FSNCL*WGLFPXs1(NB,NZ)
            ARLF1s1(K,NB,NZ)=ARLF1s1(K,NB,NZ)-FSNCL*ARLFZs1(NB,NZ)
            WGLFs1(K,NB,NZ)=WGLFs1(K,NB,NZ)-FSNCL*WGLFXs1(NB,NZ)
            WGLFNs1(K,NB,NZ)=WGLFNs1(K,NB,NZ)-FSNCL*WGLFNXs1(NB,NZ)
            WGLFPs1(K,NB,NZ)=WGLFPs1(K,NB,NZ)-FSNCL*WGLFPXs1(NB,NZ)
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
                  *(RCCN+(1.0_r8-RCCN)*RCCSXs1(NB,NZ)/WGSHEXs1(NB,NZ))
                RCPSXs1(NB,NZ)=WGSHPXs1(NB,NZ) &
                  *(RCCP+(1.0_r8-RCCP)*RCCSXs1(NB,NZ)/WGSHEXs1(NB,NZ))
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
  end associate
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
  associate(                                 &
    KVSTGs1      =>   plt_pheno%KVSTGs1    , &
    IDTHBs1      =>   plt_pheno%IDTHBs1    , &
    ISTYPs1      =>   plt_pheno%ISTYPs1    , &
    ARLFBs1      =>   plt_morph%ARLFBs1    , &
    NNODs1       =>   plt_morph%NNODs1     , &
    HTNODEs1     =>   plt_morph%HTNODEs1   , &
    HTSHEs1      =>   plt_morph%HTSHEs1    , &
    ARLF1s1      =>   plt_morph%ARLF1s1    , &
    HTNODXs1     =>   plt_morph%HTNODXs1     &
  )
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
        RCZL=WGLFNs1(K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCPL=WGLFPs1(K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
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
        RCZS=WGSHNs1(K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCPS=WGSHPs1(K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
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
        WTSHEBs1(NB,NZ)=AMAX1(0.0,WTSHEBs1(NB,NZ)-FSNCS*WGSHEs1(K,NB,NZ))
        WTSHBNs1(NB,NZ)=AMAX1(0.0,WTSHBNs1(NB,NZ)-FSNCS*WGSHNs1(K,NB,NZ))
        WTSHBPs1(NB,NZ)=AMAX1(0.0,WTSHBPs1(NB,NZ)-FSNCS*WGSHPs1(K,NB,NZ))
        HTSHEs1(K,NB,NZ)=HTSHEs1(K,NB,NZ)-FSNAS*HTSHEs1(K,NB,NZ)
        WGSHEs1(K,NB,NZ)=WGSHEs1(K,NB,NZ)-FSNCS*WGSHEs1(K,NB,NZ)
        WGSHNs1(K,NB,NZ)=WGSHNs1(K,NB,NZ)-FSNCS*WGSHNs1(K,NB,NZ)
        WGSHPs1(K,NB,NZ)=WGSHPs1(K,NB,NZ)-FSNCS*WGSHPs1(K,NB,NZ)
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
    SNCR=SNCT*(1.0_r8-SNCF)
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
    !
    !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !
        IF(WGNODEs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
          RCCK=RCSC*WGNODEs1(K,NB,NZ)
          RCZK=WGNODNs1(K,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
          RCPK=WGNODPs1(K,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
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
          WTSTKBs1(NB,NZ)=AMAX1(0.0,WTSTKBs1(NB,NZ)-FSNCK*WGNODEs1(K,NB,NZ))
          WTSTBNs1(NB,NZ)=AMAX1(0.0,WTSTBNs1(NB,NZ)-FSNCK*WGNODNs1(K,NB,NZ))
          WTSTBPs1(NB,NZ)=AMAX1(0.0,WTSTBPs1(NB,NZ)-FSNCK*WGNODPs1(K,NB,NZ))
          HTNODEs1(K,NB,NZ)=HTNODEs1(K,NB,NZ)-FSNCK*HTNODXs1(K,NB,NZ)
          WGNODEs1(K,NB,NZ)=WGNODEs1(K,NB,NZ)-FSNCK*WGNODEs1(K,NB,NZ)
          WGNODNs1(K,NB,NZ)=WGNODNs1(K,NB,NZ)-FSNCK*WGNODNs1(K,NB,NZ)
          WGNODPs1(K,NB,NZ)=WGNODPs1(K,NB,NZ)-FSNCK*WGNODPs1(K,NB,NZ)
          HTNODXs1(K,NB,NZ)=HTNODXs1(K,NB,NZ)-FSNCK*HTNODXs1(K,NB,NZ)
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
          WTSTKBs1(NB,NZ)=AMAX1(0.0,WTSTKBs1(NB,NZ)-WGNODEs1(K,NB,NZ))
          WTSTBNs1(NB,NZ)=AMAX1(0.0,WTSTBNs1(NB,NZ)-WGNODNs1(K,NB,NZ))
          WTSTBPs1(NB,NZ)=AMAX1(0.0,WTSTBPs1(NB,NZ)-WGNODPs1(K,NB,NZ))
          HTNODEs1(K,NB,NZ)=HTNODEs1(K,NB,NZ)-HTNODXs1(K,NB,NZ)
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
        RCZK=WTSTXNs1(NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
        RCPK=WTSTXPs1(NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
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
  end associate
  end subroutine RemobilizeLeafLayers

!------------------------------------------------------------------------------------------

  subroutine PhenologyReset(I,NB,NZ)

  implicit none
  integer, intent(in) :: I,nb,nz
  integer :: K,M
! begin_execution
  associate(                               &
    VRNXs1     =>  plt_pheno%VRNXs1      , &
    IWTYPs1    =>  plt_pheno%IWTYPs1     , &
    IGTYPs1    =>  plt_pheno%IGTYPs1     , &
    VRNFs1     =>  plt_pheno%VRNFs1      , &
    VRNLs1     =>  plt_pheno%VRNLs1      , &
    ISTYPs1    =>  plt_pheno%ISTYPs1     , &
    IFLGFs1    =>  plt_pheno%IFLGFs1     , &
    IFLGIs1    =>  plt_pheno%IFLGIs1     , &
    IDAYs1     =>  plt_pheno%IDAYs1      , &
    IFLGEs1    =>  plt_pheno%IFLGEs1     , &
    IBTYPs1    =>  plt_pheno%IBTYPs1     , &
    TGSTGIs1   =>  plt_pheno%TGSTGIs1    , &
    TGSTGFs1   =>  plt_pheno%TGSTGFs1    , &
    VRNSs1     =>  plt_pheno%VRNSs1      , &
    GROUPs1    =>  plt_pheno%GROUPs1     , &
    VSTGXs1    =>  plt_pheno%VSTGXs1     , &
    FLG4s1     =>  plt_pheno%FLG4s1      , &
    IFLGRs1    =>  plt_pheno%IFLGRs1     , &
    IFLGAs1    =>  plt_pheno%IFLGAs1     , &
    IFLGQs1    =>  plt_pheno%IFLGQs1     , &
    KVSTGs1    =>  plt_pheno%KVSTGs1     , &
    VSTGs1     =>  plt_morph%VSTGs1      , &
    XTLIs1     =>  plt_morph%XTLIs1      , &
    HTSHEs1    =>  plt_morph%HTSHEs1     , &
    KLEAFs1    =>  plt_morph%KLEAFs1     , &
    HTNODXs1   =>  plt_morph%HTNODXs1    , &
    HTNODEs1   =>  plt_morph%HTNODEs1    , &
    GRNXBs1    =>  plt_morph%GRNXBs1     , &
    PSTGFs1    =>  plt_morph%PSTGFs1     , &
    NB1s1      =>  plt_morph%NB1s1       , &
    NBTBs1     =>  plt_morph%NBTBs1      , &
    ARLFBs1    =>  plt_morph%ARLFBs1     , &
    ARLF1s1    =>  plt_morph%ARLF1s1     , &
    PSTGs1     =>  plt_morph%PSTGs1      , &
    PSTGIs1    =>  plt_morph%PSTGIs1     , &
    GRNOBs1    =>  plt_morph%GRNOBs1       &
  )
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
          DO 6245 M=1,JP1
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
    WTHSKBs1(NB,NZ)=(1.0_r8-FSNR)*WTHSKBs1(NB,NZ)
    WTEARBs1(NB,NZ)=(1.0_r8-FSNR)*WTEARBs1(NB,NZ)
    WTGRBs1(NB,NZ)=(1.0_r8-FSNR)*WTGRBs1(NB,NZ)
    WTHSBNs1(NB,NZ)=(1.0_r8-FSNR)*WTHSBNs1(NB,NZ)
    WTEABNs1(NB,NZ)=(1.0_r8-FSNR)*WTEABNs1(NB,NZ)
    WTGRBNs1(NB,NZ)=(1.0_r8-FSNR)*WTGRBNs1(NB,NZ)
    WTHSBPs1(NB,NZ)=(1.0_r8-FSNR)*WTHSBPs1(NB,NZ)
    WTEABPs1(NB,NZ)=(1.0_r8-FSNR)*WTEABPs1(NB,NZ)
    WTGRBPs1(NB,NZ)=(1.0_r8-FSNR)*WTGRBPs1(NB,NZ)
    GRNXBs1(NB,NZ)=(1.0_r8-FSNR)*GRNXBs1(NB,NZ)
    GRNOBs1(NB,NZ)=(1.0_r8-FSNR)*GRNOBs1(NB,NZ)
    GRWTBs1(NB,NZ)=(1.0_r8-FSNR)*GRWTBs1(NB,NZ)
!
!     STALKS BECOME LITTERFALL IN GRASSES AT END OF SEASON
!
    IF((IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1).AND.ISTYPs1(NZ).NE.0)THEN
      DO 6335 M=1,4
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+FSNR*CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+FSNR*CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+FSNR*CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
6335    CONTINUE
      WTSTKBs1(NB,NZ)=(1.0_r8-FSNR)*WTSTKBs1(NB,NZ)
      WTSTBNs1(NB,NZ)=(1.0_r8-FSNR)*WTSTBNs1(NB,NZ)
      WTSTBPs1(NB,NZ)=(1.0_r8-FSNR)*WTSTBPs1(NB,NZ)
      WTSTXBs1(NB,NZ)=(1.0_r8-FSNR)*WTSTXBs1(NB,NZ)
      WTSTXNs1(NB,NZ)=(1.0_r8-FSNR)*WTSTXNs1(NB,NZ)
      WTSTXPs1(NB,NZ)=(1.0_r8-FSNR)*WTSTXPs1(NB,NZ)
      DO 2010 K=0,JNODS1
    !     HTNODEs1(K,NB,NZ)=(1.0_r8-FSNR)*HTNODEs1(K,NB,NZ)
        HTNODXs1(K,NB,NZ)=(1.0_r8-FSNR)*HTNODXs1(K,NB,NZ)
        WGNODEs1(K,NB,NZ)=(1.0_r8-FSNR)*WGNODEs1(K,NB,NZ)
        WGNODNs1(K,NB,NZ)=(1.0_r8-FSNR)*WGNODNs1(K,NB,NZ)
        WGNODPs1(K,NB,NZ)=(1.0_r8-FSNR)*WGNODPs1(K,NB,NZ)
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
      ENDIF
    ENDIF
!     ENDIF
  ENDIF
  end associate
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
  associate(                            &
    IGTYPs1    =>  plt_pheno%IGTYPs1  , &
    ISTYPs1    =>  plt_pheno%ISTYPs1  , &
    KVSTGs1    =>  plt_pheno%KVSTGs1  , &
    IBTYPs1    =>  plt_pheno%IBTYPs1  , &
    KVSTGNs1   =>  plt_pheno%KVSTGNs1 , &
    SDPTHs1    => plt_morph%SDPTHs1   , &
    NB1s1      => plt_morph%NB1s1     , &
    ARLFLs1    => plt_morph%ARLFLs1   , &
    HTSHEs1    => plt_morph%HTSHEs1   , &
    ARSTKs1    => plt_morph%ARSTKs1   , &
    NBTBs1     => plt_morph%NBTBs1    , &
    HTNODEs1   => plt_morph%HTNODEs1  , &
    ZLs1       => plt_morph%ZLs1      , &
    ZCs1       => plt_morph%ZCs1      , &
    CLASSs1    => plt_morph%CLASSs1   , &
    ARLF1s1    => plt_morph%ARLF1s1   , &
    ARLFVs1    => plt_morph%ARLFVs1   , &
    ZSINs1     => plt_rad%ZSINs1        &
  )
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
    HTCTLs1(NZ)=XLGLF+HTSHEs1(0,NB1s1(NZ),NZ)+HTNODEs1(0,NB1s1(NZ),NZ)
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
  end associate
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
  associate(                              &
    IWTYPs1    =>  plt_pheno%IWTYPs1    , &
    VRNFs1     =>  plt_pheno%VRNFs1     , &
    FLG4s1     =>  plt_pheno%FLG4s1     , &
    VRNXs1     =>  plt_pheno%VRNXs1     , &
    ISTYPs1    =>  plt_pheno%ISTYPs1    , &
    DGSTGFs1   =>  plt_pheno%DGSTGFs1   , &
    IDAYs1     =>  plt_pheno%IDAYs1     , &
    NGs1       =>  plt_morph%NGs1       , &
    GRNXBs1    =>  plt_morph%GRNXBs1    , &
    GRNOBs1    =>  plt_morph%GRNOBs1      &
  )
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
  end associate
  end subroutine GrainFilling

!------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8) :: dangle
  integer :: L,K,N
  ! begin_execution
  associate(                             &
    ANGBRs1   =>  plt_morph%ANGBRs1    , &
    ARLF1s1   =>  plt_morph%ARLF1s1    , &
    ARSTKs1   =>  plt_morph%ARSTKs1    , &
    CLASSs1   =>  plt_morph%CLASSs1    , &
    ARLFLs1   =>  plt_morph%ARLFLs1    , &
    SURFBs1   =>  plt_morph%SURFBs1    , &
    NB1s1     =>  plt_morph%NB1s1      , &
    SURFs1    =>  plt_morph%SURFs1       &
  )
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
  end associate
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
  associate(                          &
    IDAYs1   =>  plt_pheno%IDAYs1   , &
    VRNXs1   =>  plt_pheno%VRNXs1   , &
    VRNSs1   =>  plt_pheno%VRNSs1   , &
    VRNLs1   =>  plt_pheno%VRNLs1   , &
    IFLGIs1  =>  plt_pheno%IFLGIs1  , &
    XPPDs1   =>  plt_pheno%XPPDs1   , &
    ISTYPs1  =>  plt_pheno%ISTYPs1  , &
    IPTYPs1  =>  plt_pheno%IPTYPs1  , &
    XDLs1    =>  plt_pheno%XDLs1    , &
    VRNFs1   =>  plt_pheno%VRNFs1   , &
    IBTYPs1  =>  plt_pheno%IBTYPs1  , &
    IGTYPs1  =>  plt_pheno%IGTYPs1  , &
    IWTYPs1  =>  plt_pheno%IWTYPs1  , &
    IFLGAs1  =>  plt_pheno%IFLGAs1  , &
    ATRPs1   =>   plt_pheno%ATRPs1  , &
    NGs1     =>   plt_morph%NGs1    , &
    NB1s1    =>  plt_morph%NB1s1    , &
    NIs1     =>  plt_morph%NIs1       &
  )
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!

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
  end associate
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

  associate(                            &
    IBTYPs1    =>  plt_pheno%IBTYPs1  , &
    IGTYPs1    =>  plt_pheno%IGTYPs1  , &
    RCSs1    =>  plt_photo%RCSs1      , &
    ARLFVs1  =>  plt_morph%ARLFVs1    , &
    ARSTVs1  =>  plt_morph%ARSTVs1    , &
    NRTs1    =>  plt_morph%NRTs1        &
  )
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
  end associate
  end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

  subroutine ComputeTotalBiom(NZ,CPOOLK)

  integer, intent(in) :: NZ
  real(r8), intent(out) :: CPOOLK(JC1,JP1)
  integer :: L,K,N,NB
!     begin_execution
  associate(                                 &
    NIs1         =>   plt_morph%NIs1       , &
    NBRs1        =>  plt_morph%NBRs1         &
  )
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
  end associate
  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,UPNFC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: UPNFC(JP1)
  integer :: L,NR,N,NB
!     begin_execution
  associate(                            &
    NBRs1      =>  plt_morph%NBRs1    , &
    NIs1       =>  plt_morph%NIs1     , &
    NRTs1      =>  plt_morph%NRTs1    , &
    ARLFBs1    =>  plt_morph%ARLFBs1  , &
    ARSTPs1    =>  plt_morph%ARSTPs1  , &
    ARSTKs1    =>  plt_morph%ARSTKs1  , &
    GRNOBs1    =>  plt_morph%GRNOBs1  , &
    ARLFPs1    =>  plt_morph%ARLFPs1  , &
    ARSTVs1    =>  plt_morph%ARSTVs1  , &
    GRNOs1     =>  plt_morph%GRNOs1     &
  )
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
  GRNOs1(NZ)  =sum(GRNOBs1(1:NBRs1(NZ),NZ))
  ARLFPs1(NZ)=sum(ARLFBs1(1:NBRs1(NZ),NZ))
  ARSTPs1(NZ)=sum(ARSTKs1(1:JC1,1:NBRs1(NZ),NZ))
  ARSTVs1(1:JC1,1:NBRs1(NZ))=0._r8
  DO NB=1,NBRs1(NZ)
    DO L=1,JC1
      ARSTVs1(L,NZ)=ARSTVs1(L,NZ)+ARSTKs1(L,NB,NZ)
    ENDDO
  ENDDO
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
  end associate
  end subroutine AccumulateStates



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
  associate(                            &
    IGTYPs1    =>  plt_pheno%IGTYPs1  , &
    IWTYPs1     =>  plt_pheno%IWTYPs1 , &
    FDBKXs1     => plt_photo%FDBKXs1    &
  )
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

  end associate
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
  associate(                              &
    IGTYPs1     =>  plt_pheno%IGTYPs1   , &
    IWTYPs1     =>  plt_pheno%IWTYPs1   , &
    NGs1        =>  plt_morph%NGs1      , &
    FDBKXs1     =>  plt_photo%FDBKXs1     &
  )
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  IF(CCPOLBs1(NB,NZ).GT.ZEROs1)THEN
    CNPG=AMIN1(CZPOLBs1(NB,NZ)/(CZPOLBs1(NB,NZ)+CCPOLBs1(NB,NZ)*CNKI), &
      CPPOLBs1(NB,NZ)/(CPPOLBs1(NB,NZ)+CCPOLBs1(NB,NZ)*CPKI))
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
  RCO2CM=AMAX1(0.0_r8,VMXC*CPOOLs1(NB,NZ) &
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
  RMNCS=AMAX1(0.0_r8,RMPLT*TFN6(NGs1(NZ))*WTSHXN)
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
  RCO2YM=AMAX1(0.0_r8,RCO2XM)*WFNSG
  RCO2Y=AMAX1(0.0_r8,RCO2X)*WFNSG
  SNCRM=AMAX1(0.0_r8,-RCO2XM)
  SNCR=AMAX1(0.0_r8,-RCO2X)
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
  IF(CNSHX.GT.0.0_r8.OR.CNLFX.GT.0.0_r8)THEN
    ZPOOLB=AMAX1(0.0_r8,ZPOOLs1(NB,NZ))
    PPOOLB=AMAX1(0.0_r8,PPOOLs1(NB,NZ))
    FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
    IF(RCO2YM.GT.0.0_r8)THEN
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
  ZADDBM=AMAX1(0.0_r8,CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
  ZADDB=AMAX1(0.0_r8,CGROS*(CNSHX+CNLFM+CNLFX*CNPG))
  PADDB=AMAX1(0.0_r8,CGROS*(CPSHX+CPLFM+CPLFX*CNPG))
  CNRDM=AMAX1(0.0_r8,1.70_r8*ZADDBM)
  CNRDA=AMAX1(0.0_r8,1.70_r8*ZADDB)
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
  end associate
  end subroutine ComputeRAutoBfEmergence

end module grosubsMod
