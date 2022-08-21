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
  use LitterFallMod
  use PlantBranchMod
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
    IFLGCs1    => plt_pheno%IFLGCs1   , &
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
    WTSTGPs1  => plt_biom%WTSTGPs1   , &
    WTSHPs1   => plt_biom%WTSHPs1    , &
    WTSHNs1   => plt_biom%WTSHNs1    , &
    WTRTs1    => plt_biom%WTRTs1     , &
    WTSHTs1   => plt_biom%WTSHTs1    , &
    WTNDs1    => plt_biom%WTNDs1     , &
    WTRTNs1   => plt_biom%WTRTNs1    , &
    WTNDNs1   => plt_biom%WTNDNs1    , &
    WTNDPs1   => plt_biom%WTNDPs1    , &
    WTRTPs1   => plt_biom%WTRTPs1    , &
    WTRVNs1   => plt_biom%WTRVNs1    , &
    WTRVCs1   => plt_biom%WTRVCs1    , &
    WTRVPs1   => plt_biom%WTRVPs1    , &
    WTSTGNs1  => plt_biom%WTSTGNs1   , &
    WTSTGs1   => plt_biom%WTSTGs1    , &
    RSETCs1   => plt_pheno%RSETCs1   , &
    RSETNs1   => plt_pheno%RSETNs1   , &
    RSETPs1   => plt_pheno%RSETPs1   , &
    IFLGCs1   => plt_pheno%IFLGCs1   , &
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
    WVSTKs1  =>  plt_biom%WVSTKs1     , &
    WSRTLs1  =>  plt_biom%WSRTLs1     , &
    WTRTAs1  =>  plt_biom%WTRTAs1     , &
    WTSTKs1  =>  plt_biom%WTSTKs1     , &
    WTRTs1   =>  plt_biom%WTRTs1      , &
    IBTYPs1  =>  plt_pheno%IBTYPs1    , &
    IGTYPs1  =>  plt_pheno%IGTYPs1    , &
    CNRTs1   =>  plt_allom%CNRTs1     , &
    CPRTs1   =>  plt_allom%CPRTs1     , &
    FVRNs1   =>  plt_allom%FVRNs1     , &
    FWODRNs1 =>  plt_allom%FWODRNs1   , &
    FWODLNs1 =>  plt_allom%FWODLNs1   , &
    FWODLPs1 =>  plt_allom%FWODLPs1   , &
    FWOODNs1 =>  plt_allom%FWOODNs1   , &
    FWODSPs1 =>  plt_allom%FWODSPs1   , &
    FWODRPs1 =>  plt_allom%FWODRPs1   , &
    FWOODPs1 =>  plt_allom%FWOODPs1   , &
    FWODSNs1 =>  plt_allom%FWODSNs1   , &
    FWODBs1  =>  plt_allom%FWODBs1    , &
    FWOODs1  =>  plt_allom%FWOODs1    , &
    FWODRs1  =>  plt_allom%FWODRs1    , &
    CNLFs1   =>  plt_allom%CNLFs1     , &
    CPLFs1   =>  plt_allom%CPLFs1     , &
    CNSHEs1  =>  plt_allom%CNSHEs1    , &
    CPSHEs1  =>  plt_allom%CPSHEs1    , &
    CNSTKs1  =>  plt_allom%CNSTKs1    , &
    CPSTKs1  =>  plt_allom%CPSTKs1    , &
    RCSs1    =>  plt_photo%RCSs1      , &
    MYs1     =>  plt_morph%MYs1       , &
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
  WFNS=AMIN1(1.0,AMAX1(0.0_r8,PSILGs1(NZ)-PSILM))
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
    WTLFBs1      =>  plt_biom%WTLFBs1      , &
    WTSHBNs1     =>  plt_biom%WTSHBNs1     , &
    WTLFBNs1     =>  plt_biom%WTLFBNs1     , &
    WTLFBPs1     =>  plt_biom%WTLFBPs1     , &
    WTGRBs1      =>  plt_biom%WTGRBs1      , &
    WTSHTNs1     =>  plt_biom%WTSHTNs1     , &
    WTGRBNs1     =>  plt_biom%WTGRBNs1     , &
    WTGRBPs1     =>  plt_biom%WTGRBPs1     , &
    CPOOLRs1     =>  plt_biom%CPOOLRs1     , &
    WTRTDs1      =>  plt_biom%WTRTDs1      , &
    WTSHTBs1     =>  plt_biom%WTSHTBs1     , &
    WTSHTPs1     =>  plt_biom%WTSHTPs1     , &
    WTEARBs1     =>  plt_biom%WTEARBs1     , &
    WTSTKBs1     =>  plt_biom%WTSTKBs1     , &
    WTRSBNs1     =>  plt_biom%WTRSBNs1     , &
    WTSHEBs1     =>  plt_biom%WTSHEBs1     , &
    WTSHBPs1     =>  plt_biom%WTSHBPs1     , &
    WTSTBNs1     =>  plt_biom%WTSTBNs1     , &
    WTSTBPs1     =>  plt_biom%WTSTBPs1     , &
    WTHSKBs1     =>  plt_biom%WTHSKBs1     , &
    WTEABNs1     =>  plt_biom%WTEABNs1     , &
    WTHSBPs1     =>  plt_biom%WTHSBPs1     , &
    WTRSVBs1     =>  plt_biom%WTRSVBs1     , &
    WTHSBNs1     =>  plt_biom%WTHSBNs1     , &
    WTRSBPs1     =>  plt_biom%WTRSBPs1     , &
    WTEABPs1     =>  plt_biom%WTEABPs1     , &
    CPOOLs1      =>  plt_biom%CPOOLs1      , &
    ZPOOLs1      =>  plt_biom%ZPOOLs1      , &
    PPOOLs1      =>  plt_biom%PPOOLs1      , &
    MYs1         =>  plt_morph%MYs1        , &
    NIs1         =>  plt_morph%NIs1        , &
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
    CPOOLs1    =>  plt_biom%CPOOLs1   , &
    ZPOOLs1    =>  plt_biom%ZPOOLs1   , &
    PPOOLs1    =>  plt_biom%PPOOLs1   , &
    CPOOLRs1   =>  plt_biom%CPOOLRs1  , &
    ZPOOLRs1   =>  plt_biom%ZPOOLRs1  , &
    ZPOLNPs1   =>  plt_biom%ZPOLNPs1  , &
    PPOLNBs1   =>  plt_biom%PPOLNBs1  , &
    CPOLNBs1   =>  plt_biom%CPOLNBs1  , &
    ZPOOLPs1   =>  plt_biom%ZPOOLPs1  , &
    ZPOLNBs1   =>  plt_biom%ZPOLNBs1  , &
    PPOOLRs1   =>  plt_biom%PPOOLRs1  , &
    WTRTSNs1   =>  plt_biom%WTRTSNs1  , &
    WTRTSPs1   =>  plt_biom%WTRTSPs1  , &
    WTNDBPs1   =>  plt_biom%WTNDBPs1  , &
    WTNDPs1    =>  plt_biom%WTNDPs1   , &
    WTNDs1     =>  plt_biom%WTNDs1    , &
    WTRTSs1    =>  plt_biom%WTRTSs1   , &
    WTRTNs1    =>  plt_biom%WTRTNs1   , &
    WTRTs1     =>  plt_biom%WTRTs1    , &
    WTNDBNs1   =>  plt_biom%WTNDBNs1  , &
    WTNDBs1    =>  plt_biom%WTNDBs1   , &
    WTNDNs1    =>  plt_biom%WTNDNs1   , &
    WTSHTBs1   =>  plt_biom%WTSHTBs1  , &
    WTSHTNs1   =>  plt_biom%WTSHTNs1  , &
    WTSHTPs1   =>  plt_biom%WTSHTPs1  , &
    WTSTKBs1   =>  plt_biom%WTSTKBs1  , &
    WTHSKBs1   =>  plt_biom%WTHSKBs1  , &
    WTRSVBs1   =>  plt_biom%WTRSVBs1  , &
    WTEARBs1   =>  plt_biom%WTEARBs1  , &
    WTLSBs1    =>  plt_biom%WTLSBs1   , &
    WTLFBNs1   =>  plt_biom%WTLFBNs1  , &
    WTSTBNs1   =>  plt_biom%WTSTBNs1  , &
    WTEABNs1   =>  plt_biom%WTEABNs1  , &
    WTGRBNs1   =>  plt_biom%WTGRBNs1  , &
    WTHSBPs1   =>  plt_biom%WTHSBPs1  , &
    WTEABPs1   =>  plt_biom%WTEABPs1  , &
    WTGRBPs1   =>  plt_biom%WTGRBPs1  , &
    WTSTBPs1   =>  plt_biom%WTSTBPs1  , &
    WTSHBPs1   =>  plt_biom%WTSHBPs1  , &
    WTRSBPs1   =>  plt_biom%WTRSBPs1  , &
    WTLFBPs1   =>  plt_biom%WTLFBPs1  , &
    WTHSBNs1   =>  plt_biom%WTHSBNs1  , &
    WTRSBNs1   =>  plt_biom%WTRSBNs1  , &
    WTSHBNs1   =>  plt_biom%WTSHBNs1  , &
    WTGRBs1    =>  plt_biom%WTGRBs1   , &
    WTLFBs1    =>  plt_biom%WTLFBs1   , &
    WVSTKBs1   =>  plt_biom%WVSTKBs1  , &
    WTSHEBs1   =>  plt_biom%WTSHEBs1  , &
    WTSHTs1    =>  plt_biom%WTSHTs1   , &
    WTSHNs1    =>  plt_biom%WTSHNs1   , &
    WTSHPs1    =>  plt_biom%WTSHPs1   , &
    WTLFs1     =>  plt_biom%WTLFs1    , &
    WTSHEs1    =>  plt_biom%WTSHEs1   , &
    WTSTKs1    =>  plt_biom%WTSTKs1   , &
    WTRTPs1    =>  plt_biom%WTRTPs1   , &
    WVSTKs1    =>  plt_biom%WVSTKs1   , &
    WTRSVs1    =>  plt_biom%WTRSVs1   , &
    WTHSKs1    =>  plt_biom%WTHSKs1   , &
    WTEARs1    =>  plt_biom%WTEARs1   , &
    WTGRs1     =>  plt_biom%WTGRs1    , &
    WTLSs1     =>  plt_biom%WTLSs1    , &
    WTLFNs1    =>  plt_biom%WTLFNs1   , &
    WTSHENs1   =>  plt_biom%WTSHENs1  , &
    WTSTKNs1   =>  plt_biom%WTSTKNs1  , &
    PPOOLPs1   =>  plt_biom%PPOOLPs1  , &
    PPOLNPs1   =>  plt_biom%PPOLNPs1  , &
    WTRSVNs1   =>  plt_biom%WTRSVNs1  , &
    WTHSKNs1   =>  plt_biom%WTHSKNs1  , &
    WTEARNs1   =>  plt_biom%WTEARNs1  , &
    WTGRNNs1   =>  plt_biom%WTGRNNs1  , &
    WTLFPs1    =>  plt_biom%WTLFPs1   , &
    WTSHEPs1   =>  plt_biom%WTSHEPs1  , &
    WTSTKPs1   =>  plt_biom%WTSTKPs1  , &
    WTRSVPs1   =>  plt_biom%WTRSVPs1  , &
    WTHSKPs1   =>  plt_biom%WTHSKPs1  , &
    WTEARPs1   =>  plt_biom%WTEARPs1  , &
    WTGRNPs1   =>  plt_biom%WTGRNPs1  , &
    CPOOLPs1   =>  plt_biom%CPOOLPs1  , &
    CPOLNPs1   =>  plt_biom%CPOLNPs1  , &
    NBRs1      =>  plt_morph%NBRs1    , &
    MYs1       =>  plt_morph%MYs1     , &
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

end module grosubsMod
