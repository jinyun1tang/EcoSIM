module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : test_aeqb,safe_adb,AZMAX1
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
    IFLGC    => plt_pheno%IFLGC   , &
    NP       => plt_site%NP       , &
    NP0      => plt_site%NP0      , &
    NJ       => plt_site%NJ       , &
    CNET     => plt_bgcr%CNET     , &
    HESNC    => plt_bgcr%HESNC    , &
    ESNC     => plt_bgcr%ESNC     , &
    ZC       => plt_morph%ZC        &
  )
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
!

      D9980: DO NZ=1,NP0
        D1: DO L=0,NJ
          DO K=0,1
            DO M=1,jsken
              ESNC(M,K,L,NZ,1:npelms)=0._r8
            ENDDO
          ENDDO
        ENDDO D1
        HESNC(NZ,1:npelms)=0._r8
        CNET(NZ)=0._r8
        ZCX(NZ)=ZC(NZ)
        ZC(NZ)=0._r8
      ENDDO D9980
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
      DO 9985 NZ=1,NP

! IFLGC= flag for living pft
        IF(IFLGC(NZ).EQ.1)THEN
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

  integer :: L,K,NZ,M,NB,NE
  real(r8) :: XFRC,XFRN,XFRP
!     begin_execution

  associate(                           &
    IDAY0   => plt_distb%IDAY0   , &
    IYR0    => plt_distb%IYR0    , &
    THVSTE  => plt_distb%THVSTE  , &
    HVSTE   => plt_distb%HVSTE   , &
    VPO4F   => plt_distb%VPO4F   , &
    VN2OF   => plt_distb%VN2OF   , &
    VCH4F   => plt_distb%VCH4F   , &
    VCO2F   => plt_distb%VCO2F   , &
    VNH3F   => plt_distb%VNH3F   , &
    WTSTDE  => plt_biom%WTSTDE   , &
    WTRTE   => plt_biom%WTRTE    , &
    WTSHTE  => plt_biom%WTSHTE  , &
    WTNDE   => plt_biom%WTNDE    , &
    WTRVE   => plt_biom%WTRVE    , &
    WTSTGE  => plt_biom%WTSTGE   , &
    TFN3    => plt_pheno%TFN3    , &
    RSETC   => plt_pheno%RSETC   , &
    RSETN   => plt_pheno%RSETN   , &
    RSETP   => plt_pheno%RSETP   , &
    IFLGC   => plt_pheno%IFLGC   , &
    IGTYP   => plt_pheno%IGTYP   , &
    IFLGI   => plt_pheno%IFLGI   , &
    IFLGE   => plt_pheno%IFLGE   , &
    IBTYP   => plt_pheno%IBTYP   , &
    VRNL    => plt_pheno%VRNL    , &
    VRNS    => plt_pheno%VRNS    , &
    BALE    => plt_site%BALE     , &
    NP0     => plt_site%NP0      , &
    NJ      => plt_site%NJ       , &
    IYRC    => plt_site%IYRC     , &
    ESNC    => plt_bgcr%ESNC     , &
    ZNPP    => plt_bgcr%ZNPP     , &
    TZUPFX  => plt_bgcr%TZUPFX   , &
    HESNC   => plt_bgcr%HESNC    , &
    TNH3C   => plt_bgcr%TNH3C    , &
    TESNC   => plt_bgcr%TESNC    , &
    TESN0   => plt_bgcr%TESN0    , &
    TCO2T   => plt_bgcr%TCO2T    , &
    CARBN   => plt_bgcr%CARBN    , &
    TEUPTK  => plt_rbgc%TEUPTK   , &
    CDPTHZ  => plt_site%CDPTHZ   , &
    SDPTHI  => plt_morph%SDPTHI  , &
    NBR     => plt_morph%NBR       &
  )
  D9975: DO NZ=1,NP0
!
!     ACTIVATE DORMANT SEEDS
!
    D205: DO NB=1,NBR(NZ)
      IF(IFLGI(NZ).EQ.1)THEN
        IF(IFLGE(NB,NZ).EQ.0 &
          .AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
          IDAY0(NZ)=I
          IYR0(NZ)=IYRC
          SDPTHI(NZ)=0.005_r8+CDPTHZ(0)
          IFLGI(NZ)=0
        ENDIF
      ENDIF
    ENDDO D205
!
!     LITTERFALL FROM STANDING DEAD
!
!     XFRC,XFRN,XFRP=litterfall from standing dead
!     TFN3=temperature function for canopy growth
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall
!
    D6235: DO M=1,jsken
      XFRC=1.5814E-05*TFN3(NZ)*WTSTDE(M,NZ,ielmc)
      XFRN=1.5814E-05*TFN3(NZ)*WTSTDE(M,NZ,ielmn)
      XFRP=1.5814E-05*TFN3(NZ)*WTSTDE(M,NZ,ielmp)
      IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
        ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc)+XFRC
        ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn)+XFRN
        ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp)+XFRP
      ELSE
        ESNC(M,0,0,NZ,ielmc)=ESNC(M,0,0,NZ,ielmc)+XFRC
        ESNC(M,0,0,NZ,ielmn)=ESNC(M,0,0,NZ,ielmn)+XFRN
        ESNC(M,0,0,NZ,ielmp)=ESNC(M,0,0,NZ,ielmp)+XFRP
      ENDIF
      WTSTDE(M,NZ,ielmc)=WTSTDE(M,NZ,ielmc)-XFRC
      WTSTDE(M,NZ,ielmn)=WTSTDE(M,NZ,ielmn)-XFRN
      WTSTDE(M,NZ,ielmp)=WTSTDE(M,NZ,ielmp)-XFRP
    ENDDO D6235
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LITTERFALL
!
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
!
    DO NE=1,npelms
      D6430: DO M=1,jsken
        DO K=0,1
          TESN0(NZ,NE)=TESN0(NZ,ielmc)+ESNC(M,K,0,NZ,NE)
          D8955: DO L=0,NJ
            HESNC(NZ,NE)=HESNC(NZ,NE)+ESNC(M,K,L,NZ,NE)
            TESNC(NZ,NE)=TESNC(NZ,NE)+ESNC(M,K,L,NZ,NE)
          ENDDO D8955
        enddo
      ENDDO D6430
    ENDDO
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    WTSTGE(NZ,ielmc)=sum(WTSTDE(1:4,NZ,ielmc))
    WTSTGE(NZ,ielmn)=sum(WTSTDE(1:4,NZ,ielmn))
    WTSTGE(NZ,ielmp)=sum(WTSTDE(1:4,NZ,ielmp))
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
    ZNPP(NZ)=CARBN(NZ)+TCO2T(NZ)
    IF(IFLGC(NZ).EQ.1)THEN
      BALE(NZ,ielmc)=WTSHTE(NZ,ielmc)+WTRTE(NZ,ielmc)+WTNDE(NZ,ielmc) &
        +WTRVE(NZ,ielmc)-ZNPP(NZ)+TESNC(NZ,ielmc)-TEUPTK(NZ,ielmc) &
        -RSETC(NZ)+WTSTGE(NZ,ielmc)+THVSTE(NZ,ielmc) &
        +HVSTE(NZ,ielmc)-VCO2F(NZ)-VCH4F(NZ)
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
      BALE(NZ,ielmn)=WTSHTE(NZ,ielmn)+WTRTE(NZ,ielmn)+WTNDE(NZ,ielmn) &
        +WTRVE(NZ,ielmn)+TESNC(NZ,ielmn)-TEUPTK(NZ,ielmn)-TNH3C(NZ) &
        -RSETN(NZ)+WTSTGE(NZ,ielmn)+HVSTE(NZ,ielmn)+THVSTE(NZ,ielmn) &
        -VNH3F(NZ)-VN2OF(NZ)-TZUPFX(NZ)
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
      BALE(NZ,ielmp)=WTSHTE(NZ,ielmp)+WTRTE(NZ,ielmp)+WTNDE(NZ,ielmp) &
        +WTRVE(NZ,ielmp)+TESNC(NZ,ielmp)-TEUPTK(NZ,ielmp) &
        -RSETP(NZ)+WTSTDE(1,NZ,ielmp)+WTSTGE(NZ,ielmp) &
        +HVSTE(NZ,ielmp)+THVSTE(NZ,ielmp)-VPO4F(NZ)
    ENDIF
  ENDDO D9975
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
    IGTYP  => plt_pheno%IGTYP       , &
    IDTHR  => plt_pheno%IDTHR       , &
    IDTHP  => plt_pheno%IDTHP       , &
    UPNF   => plt_rbgc%UPNF         , &
    UPH2P  => plt_rbgc%UPH2P        , &
    UPNH4  => plt_rbgc%UPNH4        , &
    UPH1P  => plt_rbgc%UPH1P        , &
    UPNO3  => plt_rbgc%UPNO3        , &
    HPUPTK => plt_rbgc%HPUPTK       , &
    HZUPTK => plt_rbgc%HZUPTK       , &
    HCUPTK => plt_rbgc%HCUPTK       , &
    UPOMC  => plt_rbgc%UPOMC        , &
    UPOMN  => plt_rbgc%UPOMN        , &
    UPOMP  => plt_rbgc%UPOMP        , &
    SDAR   => plt_morph%SDAR        , &
    SDVL   => plt_morph%SDVL        , &
    NBR    => plt_morph%NBR         , &
    NRT    => plt_morph%NRT           &
  )
  IF(IDTHP(NZ).EQ.0.OR.IDTHR(NZ).EQ.0)THEN
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
    DO  NB=1,NBR(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6,ZCX,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WFNG,WFNC,WFNS,WFNSG,PTRT,UPNFC,IFLGZ)
    ENDDO
!
    call RootBGCModel(I,J,NZ,IFLGZ,ICHK1,IDTHRN,NRX,PTRT,TFN6,CNRTW,CPRTW,XRTN1)

!
    call ComputeTotalBiom(NZ,CPOOLK)
  ELSE
    HCUPTK(NZ)=UPOMC(NZ)
    HZUPTK(NZ)=UPOMN(NZ)+UPNH4(NZ)+UPNO3(NZ)+UPNF(NZ)
    HPUPTK(NZ)=UPOMP(NZ)+UPH2P(NZ)+UPH1P(NZ)
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
    TKC    =>  plt_ew%TKC         , &
    TKS    =>  plt_ew%TKS         , &
    PSILT  =>  plt_ew%PSILT       , &
    PSILG  =>  plt_ew%PSILG       , &
    PP     =>  plt_site%PP        , &
    NU     =>  plt_site%NU        , &
    NJ     =>  plt_site%NJ        , &
    OFFST  =>  plt_pheno%OFFST    , &
    ZEROP  =>  plt_biom%ZEROP     , &
    WVSTK  =>  plt_biom%WVSTK     , &
    WSRTL  =>  plt_biom%WSRTL     , &
    WTRTA  =>  plt_biom%WTRTA     , &
    WTSTKE =>  plt_biom%WTSTKE    , &
    WTRTE  =>  plt_biom%WTRTE     , &
    WGLFV  =>  plt_biom%WGLFV     , &
    IBTYP  =>  plt_pheno%IBTYP    , &
    IGTYP  =>  plt_pheno%IGTYP    , &
    RCO2A  =>  plt_rbgc%RCO2A     , &
    RCO2M  =>  plt_rbgc%RCO2M     , &
    RCO2N  =>  plt_rbgc%RCO2N     , &
    CNRT   =>  plt_allom%CNRT     , &
    CPRT   =>  plt_allom%CPRT     , &
    FVRN   =>  plt_allom%FVRN     , &
    FWODRN =>  plt_allom%FWODRN   , &
    FWODLN =>  plt_allom%FWODLN   , &
    FWODLP =>  plt_allom%FWODLP   , &
    FWOODN =>  plt_allom%FWOODN   , &
    FWODSP =>  plt_allom%FWODSP   , &
    FWODRP =>  plt_allom%FWODRP   , &
    FWOODP =>  plt_allom%FWOODP   , &
    FWODSN =>  plt_allom%FWODSN   , &
    FWODB  =>  plt_allom%FWODB    , &
    FWOOD  =>  plt_allom%FWOOD    , &
    FWODR  =>  plt_allom%FWODR    , &
    CNLF   =>  plt_allom%CNLF     , &
    CPLF   =>  plt_allom%CPLF     , &
    CNSHE  =>  plt_allom%CNSHE    , &
    CPSHE  =>  plt_allom%CPSHE    , &
    CNSTK  =>  plt_allom%CNSTK    , &
    CPSTK  =>  plt_allom%CPSTK    , &
    RCS    =>  plt_photo%RCS      , &
    RTN1   =>  plt_morph%RTN1     , &
    RTNL   =>  plt_morph%RTNL     , &
    MY     =>  plt_morph%MY       , &
    ARLFV  =>  plt_morph%ARLFV    , &
    ARSTV  =>  plt_morph%ARSTV    , &
    NRT    =>  plt_morph%NRT        &
  )
  DO 2 L=1,JC1
    ARLFV(L,NZ)=0._r8
    WGLFV(L,NZ)=0._r8
    ARSTV(L,NZ)=0._r8
2 CONTINUE
  DO 5 NR=1,NRT(NZ)
    DO  N=1,MY(NZ)
      NRX(N,NR)=0
      ICHK1(N,NR)=0
    enddo
5 CONTINUE
  DO 9 N=1,MY(NZ)
    DO 6 L=NU,NJ
      WSRTL(N,L,NZ)=0._r8
      RTN1(N,L,NZ)=0._r8
      RTNL(N,L,NZ)=0._r8
      RCO2M(N,L,NZ)=0._r8
      RCO2N(N,L,NZ)=0._r8
      RCO2A(N,L,NZ)=0._r8
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
  IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1 &
    .OR.WTSTKE(NZ,ielmc).LE.ZEROP(NZ))THEN
    FWODB(1)=1.0_r8
    FWOOD(1)=1.0_r8
    FWODR(1)=1.0_r8
  ELSE
    FWODB(1)=1.0_r8
    FWOOD(1)=SQRT(WVSTK(NZ)/WTSTKE(NZ,ielmc))
    FWODR(1)=SQRT(FRTX*WVSTK(NZ)/WTSTKE(NZ,ielmc))
  ENDIF
  FWODB(0)=1.0_r8-FWODB(1)
  FWOOD(0)=1.0_r8-FWOOD(1)
  FWODR(0)=1.0_r8-FWODR(1)
  CNLFW=FWODB(0)*CNSTK(NZ)+FWODB(1)*CNLF(NZ)
  CPLFW=FWODB(0)*CPSTK(NZ)+FWODB(1)*CPLF(NZ)
  CNSHW=FWODB(0)*CNSTK(NZ)+FWODB(1)*CNSHE(NZ)
  CPSHW=FWODB(0)*CPSTK(NZ)+FWODB(1)*CPSHE(NZ)
  CNRTW=FWODR(0)*CNSTK(NZ)+FWODR(1)*CNRT(NZ)
  CPRTW=FWODR(0)*CPSTK(NZ)+FWODR(1)*CPRT(NZ)
  FWODLN(0)=FWODB(0)*CNSTK(NZ)/CNLFW
  FWODLP(0)=FWODB(0)*CPSTK(NZ)/CPLFW
  FWODSN(0)=FWODB(0)*CNSTK(NZ)/CNSHW
  FWODSP(0)=FWODB(0)*CPSTK(NZ)/CPSHW
  FWOODN(0)=FWOOD(0)*CNSTK(NZ)/CNRTW
  FWOODP(0)=FWOOD(0)*CPSTK(NZ)/CPRTW
  FWODRN(0)=FWODR(0)*CNRT(NZ)/CNRTW
  FWODRP(0)=FWODR(0)*CPRT(NZ)/CPRTW
  FWODLN(1)=1.0_r8-FWODLN(0)
  FWODLP(1)=1.0_r8-FWODLP(0)
  FWODSN(1)=1.0_r8-FWODSN(0)
  FWODSP(1)=1.0_r8-FWODSP(0)
  FWOODN(1)=1.0_r8-FWOODN(0)
  FWOODP(1)=1.0_r8-FWOODP(0)
  FWODRN(1)=1.0_r8-FWODRN(0)
  FWODRP(1)=1.0_r8-FWODRP(0)
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
  TKCM=TKC(NZ)+OFFST(NZ)
  RTK=8.3143*TKCM
  STK=710.0*TKCM
  ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
  TFN5=EXP(25.214-62500/RTK)/ACTVM
  DO 7 L=NU,NJ
    TKSM=TKS(L)+OFFST(NZ)
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
  WTRTA(NZ)=AMAX1(0.999992087*WTRTA(NZ),WTRTE(NZ,ielmc)/PP(NZ))
  XRTN1=AMAX1(1.0,WTRTA(NZ)**0.667)*PP(NZ)
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
  WFNS=AMIN1(1.0,AZMAX1(PSILG(NZ)-PSILM))
  IF(IGTYP(NZ).EQ.0)THEN
    WFNC=1.0_r8
    WFNG=EXP(0.05*PSILT(NZ))
    WFNSG=WFNS**0.10
  ELSE
    WFNC=EXP(RCS(NZ)*PSILG(NZ))
    WFNG=EXP(0.10*PSILT(NZ))
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
    WTLFB      =>  plt_biom%WTLFB      , &
    WTSHBN     =>  plt_biom%WTSHBN     , &
    WTLFBN     =>  plt_biom%WTLFBN     , &
    WTLFBP     =>  plt_biom%WTLFBP     , &
    WTGRB      =>  plt_biom%WTGRB      , &
    WTSHTN     =>  plt_biom%WTSHTN     , &
    WTGRBN     =>  plt_biom%WTGRBN     , &
    WTGRBP     =>  plt_biom%WTGRBP     , &
    EPOOLR     =>  plt_biom%EPOOLR     , &
    WTRTD      =>  plt_biom%WTRTD      , &
    WTSHTB     =>  plt_biom%WTSHTB     , &
    WTSHTP     =>  plt_biom%WTSHTP     , &
    WTEARB     =>  plt_biom%WTEARB     , &
    WTSTKB     =>  plt_biom%WTSTKB     , &
    WTRSBN     =>  plt_biom%WTRSBN     , &
    WTSHEB     =>  plt_biom%WTSHEB     , &
    WTSHBP     =>  plt_biom%WTSHBP     , &
    WTSTBN     =>  plt_biom%WTSTBN     , &
    WTSTBP     =>  plt_biom%WTSTBP     , &
    WTHSKB     =>  plt_biom%WTHSKB     , &
    WTEABN     =>  plt_biom%WTEABN     , &
    WTHSBP     =>  plt_biom%WTHSBP     , &
    WTRSVB     =>  plt_biom%WTRSVB     , &
    WTHSBN     =>  plt_biom%WTHSBN     , &
    WTRSBP     =>  plt_biom%WTRSBP     , &
    WTEABP     =>  plt_biom%WTEABP     , &
    CPOOL      =>  plt_biom%CPOOL      , &
    ZPOOL      =>  plt_biom%ZPOOL      , &
    PPOOL      =>  plt_biom%PPOOL      , &
    NU         =>  plt_site%NU         , &
    CPOOL3     =>  plt_photo%CPOOL3    , &
    CPOOL4     =>  plt_photo%CPOOL4    , &
    HCOB       =>  plt_photo%HCOB      , &
    CO2B       =>  plt_photo%CO2B      , &
    MY         =>  plt_morph%MY        , &
    NI         =>  plt_morph%NI        , &
    NBR        =>  plt_morph%NBR         &
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
  DO 320 NB=1,NBR(NZ)
    CPOOLK(NB,NZ)=0._r8
    DO 325 K=1,JNODS1
      CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
        +CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
        +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
325   CONTINUE
    WTSHTB(NB,NZ)=WTLFB(NB,NZ) &
      +WTSHEB(NB,NZ)+WTSTKB(NB,NZ)+WTRSVB(NB,NZ) &
      +WTHSKB(NB,NZ)+WTEARB(NB,NZ)+WTGRB(NB,NZ) &
      +CPOOL(NB,NZ)+CPOOLK(NB,NZ)
    WTSHTN(NB,NZ)=WTLFBN(NB,NZ) &
      +WTSHBN(NB,NZ)+WTSTBN(NB,NZ)+WTRSBN(NB,NZ) &
      +WTHSBN(NB,NZ)+WTEABN(NB,NZ)+WTGRBN(NB,NZ) &
      +ZPOOL(NB,NZ)
    WTSHTP(NB,NZ)=WTLFBP(NB,NZ) &
      +WTSHBP(NB,NZ)+WTSTBP(NB,NZ)+WTRSBP(NB,NZ) &
      +WTHSBP(NB,NZ)+WTEABP(NB,NZ)+WTGRBP(NB,NZ) &
      +PPOOL(NB,NZ)
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
  D345: DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      WTRTD(N,L,NZ)=WTRTD(N,L,NZ)+EPOOLR(N,L,NZ,ielmc)
    enddo
  ENDDO D345
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
    CPOOL    =>  plt_biom%CPOOL   , &
    ZPOOL    =>  plt_biom%ZPOOL   , &
    PPOOL    =>  plt_biom%PPOOL   , &
    EPOOLR   =>  plt_biom%EPOOLR  , &
    ZPOLNP   =>  plt_biom%ZPOLNP  , &
    PPOLNB   =>  plt_biom%PPOLNB  , &
    CPOLNB   =>  plt_biom%CPOLNB  , &
    ZPOLNB   =>  plt_biom%ZPOLNB  , &
    WTNDBP   =>  plt_biom%WTNDBP  , &
    WTNDE    =>  plt_biom%WTNDE   , &
    WTRTSE   =>  plt_biom%WTRTSE  , &
    WTRTE    =>  plt_biom%WTRTE   , &
    WTNDBN   =>  plt_biom%WTNDBN  , &
    WTNDB    =>  plt_biom%WTNDB   , &
    WTSHTB   =>  plt_biom%WTSHTB  , &
    WTSHTN   =>  plt_biom%WTSHTN  , &
    WTSHTP   =>  plt_biom%WTSHTP  , &
    WTSTKB   =>  plt_biom%WTSTKB  , &
    WTHSKB   =>  plt_biom%WTHSKB  , &
    WTRSVB   =>  plt_biom%WTRSVB  , &
    WTEARB   =>  plt_biom%WTEARB  , &
    WTLSB    =>  plt_biom%WTLSB   , &
    WTLFBN   =>  plt_biom%WTLFBN  , &
    WTSTBN   =>  plt_biom%WTSTBN  , &
    WTEABN   =>  plt_biom%WTEABN  , &
    WTGRBN   =>  plt_biom%WTGRBN  , &
    WTHSBP   =>  plt_biom%WTHSBP  , &
    WTEABP   =>  plt_biom%WTEABP  , &
    WTGRBP   =>  plt_biom%WTGRBP  , &
    WTSTBP   =>  plt_biom%WTSTBP  , &
    WTSHBP   =>  plt_biom%WTSHBP  , &
    WTRSBP   =>  plt_biom%WTRSBP  , &
    WTLFBP   =>  plt_biom%WTLFBP  , &
    WTHSBN   =>  plt_biom%WTHSBN  , &
    WTRSBN   =>  plt_biom%WTRSBN  , &
    WTSHBN   =>  plt_biom%WTSHBN  , &
    WTGRB    =>  plt_biom%WTGRB   , &
    WTLFB    =>  plt_biom%WTLFB   , &
    WVSTKB   =>  plt_biom%WVSTKB  , &
    WTSHEB   =>  plt_biom%WTSHEB  , &
    WTSHTE   =>  plt_biom%WTSHTE  , &
    WTLFE    =>  plt_biom%WTLFE   , &
    WTSHEE   =>  plt_biom%WTSHEE  , &
    WTSTKE   =>  plt_biom%WTSTKE  , &
    WVSTK    =>  plt_biom%WVSTK   , &
    WTRSVE   =>  plt_biom%WTRSVE  , &
    WTHSKE   =>  plt_biom%WTHSKE  , &
    WTEARE   =>  plt_biom%WTEARE  , &
    WTGRE    =>  plt_biom%WTGRE   , &
    WTLS     =>  plt_biom%WTLS    , &
    PPOLNP   =>  plt_biom%PPOLNP  , &
    WTNDL    =>  plt_biom%WTNDL   , &
    WTNDLN   =>  plt_biom%WTNDLN  , &
    EPOOLP   =>  plt_biom%EPOOLP  , &
    CPOLNP   =>  plt_biom%CPOLNP  , &
    WTRT1    =>  plt_biom%WTRT1   , &
    WTRT1N   =>  plt_biom%WTRT1N  , &
    WTRT1P   =>  plt_biom%WTRT1P  , &
    WTNDLP   =>  plt_biom%WTNDLP  , &
    WTRT2    =>  plt_biom%WTRT2   , &
    WTRT2N   =>  plt_biom%WTRT2N  , &
    WTRT2P   =>  plt_biom%WTRT2P  , &
    CPOOLN   =>  plt_biom%CPOOLN  , &
    ZPOOLN   =>  plt_biom%ZPOOLN  , &
    PPOOLN   =>  plt_biom%PPOOLN  , &
    TZUPFX   =>  plt_bgcr%TZUPFX  , &
    TEUPTK   =>  plt_rbgc%TEUPTK  , &
    UPH1P    =>  plt_rbgc%UPH1P   , &
    HCUPTK   =>  plt_rbgc%HCUPTK  , &
    HZUPTK   =>  plt_rbgc%HZUPTK  , &
    HPUPTK   =>  plt_rbgc%HPUPTK  , &
    UPNF     =>  plt_rbgc%UPNF    , &
    UPH2P    =>  plt_rbgc%UPH2P   , &
    UPNO3    =>  plt_rbgc%UPNO3   , &
    UPNH4    =>  plt_rbgc%UPNH4   , &
    UPOMC    =>  plt_rbgc%UPOMC   , &
    UPOMN    =>  plt_rbgc%UPOMN   , &
    UPOMP    =>  plt_rbgc%UPOMP   , &
    NJ       =>  plt_site%NJ      , &
    NU       =>  plt_site%NU      , &
    NBR      =>  plt_morph%NBR    , &
    MY       =>  plt_morph%MY     , &
    NI       =>  plt_morph%NI     , &
    NRT      =>  plt_morph%NRT    , &
    ARLFB    =>  plt_morph%ARLFB  , &
    ARSTP    =>  plt_morph%ARSTP  , &
    ARSTK    =>  plt_morph%ARSTK  , &
    GRNOB    =>  plt_morph%GRNOB  , &
    ARLFP    =>  plt_morph%ARLFP  , &
    ARSTV    =>  plt_morph%ARSTV  , &
    INTYP    =>  plt_morph%INTYP  , &
    GRNO     =>  plt_morph%GRNO     &
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
  EPOOLP(NZ,ielmc)=sum(CPOOL(1:NBR(NZ),NZ))
  EPOOLP(NZ,ielmn)=sum(ZPOOL(1:NBR(NZ),NZ))
  EPOOLP(NZ,ielmp)=sum(PPOOL(1:NBR(NZ),NZ))
  WTSHTE(NZ,ielmc)=sum(WTSHTB(1:NBR(NZ),NZ))
  WTSHTE(NZ,ielmn)=sum(WTSHTN(1:NBR(NZ),NZ))
  WTSHTE(NZ,ielmp)=sum(WTSHTP(1:NBR(NZ),NZ))
  WTLFE(NZ,ielmc)=sum(WTLFB(1:NBR(NZ),NZ))
  WTSHEE(NZ,ielmc)=sum(WTSHEB(1:NBR(NZ),NZ))
  WTSTKE(NZ,ielmc)=sum(WTSTKB(1:NBR(NZ),NZ))
  WVSTK(NZ)=sum(WVSTKB(1:NBR(NZ),NZ))
  WTRSVE(NZ,ielmc)=sum(WTRSVB(1:NBR(NZ),NZ))
  WTHSKE(NZ,ielmc)=sum(WTHSKB(1:NBR(NZ),NZ))
  WTEARE(NZ,ielmc)=sum(WTEARB(1:NBR(NZ),NZ))
  WTGRE(NZ,ielmc)=sum(WTGRB(1:NBR(NZ),NZ))
  WTLS(NZ)=sum(WTLSB(1:NBR(NZ),NZ))
  WTLFE(NZ,ielmn)=sum(WTLFBN(1:NBR(NZ),NZ))
  WTSHEE(NZ,ielmn)=sum(WTSHBN(1:NBR(NZ),NZ))
  WTSTKE(NZ,ielmn)=sum(WTSTBN(1:NBR(NZ),NZ))
  WTRSVE(NZ,ielmn)=sum(WTRSBN(1:NBR(NZ),NZ))
  WTHSKE(NZ,ielmn)=sum(WTHSBN(1:NBR(NZ),NZ))
  WTEARE(NZ,ielmn)=sum(WTEABN(1:NBR(NZ),NZ))
  WTGRE(NZ,ielmn)=sum(WTGRBN(1:NBR(NZ),NZ))
  WTLFE(NZ,ielmp)=sum(WTLFBP(1:NBR(NZ),NZ))
  WTSHEE(NZ,ielmp)=sum(WTSHBP(1:NBR(NZ),NZ))
  WTSTKE(NZ,ielmp)=sum(WTSTBP(1:NBR(NZ),NZ))
  WTRSVE(NZ,ielmp)=sum(WTRSBP(1:NBR(NZ),NZ))
  WTHSKE(NZ,ielmp)=sum(WTHSBP(1:NBR(NZ),NZ))
  WTEARE(NZ,ielmp)=sum(WTEABP(1:NBR(NZ),NZ))
  WTGRE(NZ,ielmp)=sum(WTGRBP(1:NBR(NZ),NZ))
  GRNO(NZ)  =sum(GRNOB(1:NBR(NZ),NZ))
  ARLFP(NZ)=sum(ARLFB(1:NBR(NZ),NZ))
  ARSTP(NZ)=sum(ARSTK(1:JC1,1:NBR(NZ),NZ))
  ARSTV(1:JC1,1:NBR(NZ))=0._r8
  DO NB=1,NBR(NZ)
    DO L=1,JC1
      ARSTV(L,NZ)=ARSTV(L,NZ)+ARSTK(L,NB,NZ)
    ENDDO
  ENDDO
!
!     ACCUMULATE ROOT STATE VARIABLES FROM ROOT LAYER STATE VARIABLES
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!

  WTRTE(NZ,ielmc)=sum(EPOOLR(1:MY(NZ),NU:NJ,NZ,ielmc))
  WTRTE(NZ,ielmn)=sum(EPOOLR(1:MY(NZ),NU:NJ,NZ,ielmn))
  WTRTE(NZ,ielmp)=sum(EPOOLR(1:MY(NZ),NU:NJ,NZ,ielmp))
  WTRTSE(NZ,ielmc)=sum(WTRT1(1:MY(NZ),NU:NJ,1:NRT(NZ),NZ)) &
    +sum(WTRT2(1:MY(NZ),NU:NJ,1:NRT(NZ),NZ))
  WTRTSE(NZ,ielmn)=sum(WTRT1N(1:MY(NZ),NU:NJ,1:NRT(NZ),NZ)) &
    +sum(WTRT2N(1:MY(NZ),NU:NJ,1:NRT(NZ),NZ))
  WTRTSE(NZ,ielmp)=sum(WTRT1P(1:MY(NZ),NU:NJ,1:NRT(NZ),NZ)) &
    +sum(WTRT2P(1:MY(NZ),NU:NJ,1:NRT(NZ),NZ))
  WTRTE(NZ,ielmc)=WTRTE(NZ,ielmc)+WTRTSE(NZ,ielmc)
  WTRTE(NZ,ielmn)=WTRTE(NZ,ielmn)+WTRTSE(NZ,ielmn)
  WTRTE(NZ,ielmp)=WTRTE(NZ,ielmp)+WTRTSE(NZ,ielmp)
!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  IF(INTYP(NZ).NE.0)THEN
    IF(INTYP(NZ).GE.4)THEN
      DO 7950 NB=1,NBR(NZ)
        CPOLNP(NZ)=CPOLNP(NZ)+CPOLNB(NB,NZ)
        ZPOLNP(NZ)=ZPOLNP(NZ)+ZPOLNB(NB,NZ)
        PPOLNP(NZ)=PPOLNP(NZ)+PPOLNB(NB,NZ)
7950  CONTINUE
      WTNDE(NZ,ielmc)=sum(WTNDB(1:NBR(NZ),NZ))+&
        sum(CPOLNB(1:NBR(NZ),NZ))
      WTNDE(NZ,ielmn)=sum(WTNDBN(1:NBR(NZ),NZ))+&
        sum(ZPOLNB(1:NBR(NZ),NZ))
      WTNDE(NZ,ielmp)=sum(WTNDBP(1:NBR(NZ),NZ))+&
        sum(PPOLNB(1:NBR(NZ),NZ))
    ELSEIF(INTYP(NZ).GE.1.AND.INTYP(NZ).LE.3)THEN
      WTNDE(NZ,ielmc)=sum(WTNDL(NU:NI(NZ),NZ))+&
        sum(CPOOLN(NU:NI(NZ),NZ))
      WTNDE(NZ,ielmn)=sum(WTNDLN(NU:NI(NZ),NZ))+&
        sum(ZPOOLN(NU:NI(NZ),NZ))
      WTNDE(NZ,ielmp)=sum(WTNDLP(NU:NI(NZ),NZ))+&
        sum(PPOOLN(NU:NI(NZ),NZ))
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
  HCUPTK(NZ)=UPOMC(NZ)
  HZUPTK(NZ)=UPOMN(NZ)+UPNH4(NZ)+UPNO3(NZ)+UPNF(NZ)
  HPUPTK(NZ)=UPOMP(NZ)+UPH2P(NZ)+UPH1P(NZ)
  TEUPTK(NZ,ielmc)=TEUPTK(NZ,ielmc)+UPOMC(NZ)
  TEUPTK(NZ,ielmn)=TEUPTK(NZ,ielmn)+UPOMN(NZ)+UPNH4(NZ)+UPNO3(NZ)
  TEUPTK(NZ,ielmp)=TEUPTK(NZ,ielmp)+UPOMP(NZ)+UPH2P(NZ)+UPH1P(NZ)
  TZUPFX(NZ)=TZUPFX(NZ)+UPNF(NZ)+UPNFC(NZ)
  end associate
  end subroutine AccumulateStates

end module grosubsMod
