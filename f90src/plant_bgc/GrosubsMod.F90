module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : test_aeqb,safe_adb,AZMAX1
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use EcoSiMParDataMod, only : pltpar
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

  subroutine InitGrosub(jpstgs,JRS)

  implicit none
  integer, intent(out) :: jpstgs,JRS

  call InitVegPars(pltpar)
  jpstgs = pltpar%jpstgs
  jrs = pltpar%JRS


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
      DO K=1,pltpar%n_pltlitrk
        DO M=1,jsken
          ESNC(M,1:npelms,K,L,NZ)=0._r8
        ENDDO
      ENDDO
    ENDDO D1
    HESNC(1:npelms,NZ)=0._r8
    CNET(NZ)=0._r8
    ZCX(NZ)=ZC(NZ)
    ZC(NZ)=0._r8
  ENDDO D9980
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
  D9985: DO NZ=1,NP

! IFLGC= flag for living pft
    IF(IFLGC(NZ).EQ.1)THEN
      call GrowPlant(I,J,NZ,ZCX,CPOOLK)
    ENDIF

!     HARVEST STANDING DEAD

    call RemoveBiomassByDisturbance(I,J,NZ,CPOOLK)
  ENDDO D9985
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
  real(r8) :: XFRC,XFRN,XFRP,XFRE
!     begin_execution

  associate(                       &
    k_fine_litr=> pltpar%k_fine_litr ,&
    k_woody_litr=> pltpar%k_woody_litr,&
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
    RSETE   => plt_pheno%RSETE   , &
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
        IF(IFLGE(NB,NZ).EQ.0.AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
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
    DO NE=1,npelms
      D6235: DO M=1,jsken
        XFRE=1.5814E-05_r8*TFN3(NZ)*WTSTDE(M,NE,NZ)
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          ESNC(M,NE,k_fine_litr,0,NZ)=ESNC(M,NE,k_fine_litr,0,NZ)+XFRE
        ELSE
          ESNC(M,NE,k_woody_litr,0,NZ)=ESNC(M,NE,k_woody_litr,0,NZ)+XFRE
        ENDIF
        WTSTDE(M,NE,NZ)=WTSTDE(M,NE,NZ)-XFRE
      ENDDO D6235
    ENDDO
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LITTERFALL
!
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
!
    DO K=1,pltpar%n_pltlitrk
      DO NE=1,npelms
        D6430: DO M=1,jsken
          TESN0(NE,NZ)=TESN0(NE,NZ)+ESNC(M,NE,K,0,NZ)
          D8955: DO L=0,NJ
            HESNC(NE,NZ)=HESNC(NE,NZ)+ESNC(M,NE,K,L,NZ)
            TESNC(NE,NZ)=TESNC(NE,NZ)+ESNC(M,NE,K,L,NZ)
          ENDDO D8955
        ENDDO D6430
      enddo
    ENDDO
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    DO NE=1,npelms
      WTSTGE(NE,NZ)=sum(WTSTDE(1:jsken,NE,NZ))
    ENDDO
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
      DO NE=1,npelms
        BALE(NE,NZ)=WTSHTE(NE,NZ)+WTRTE(NE,NZ)+WTNDE(NE,NZ) &
          +WTRVE(NE,NZ)+TESNC(NE,NZ)-TEUPTK(NE,NZ) &
          -RSETE(NE,NZ)+WTSTGE(NE,NZ)+HVSTE(NE,NZ)+THVSTE(NE,NZ)
      ENDDO
      BALE(ielmc,NZ)=BALE(ielmc,NZ)-ZNPP(NZ)-VCO2F(NZ)-VCH4F(NZ)
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
      BALE(ielmn,NZ)=BALE(ielmn,NZ)-TNH3C(NZ)-VNH3F(NZ)-VN2OF(NZ)-TZUPFX(NZ)
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
      BALE(ielmp,NZ)=BALE(ielmp,NZ)-VPO4F(NZ)
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
    HEUPTK => plt_rbgc%HEUPTK       , &
    UPOME  => plt_rbgc%UPOME        , &
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
    HEUPTK(1:npelms,NZ)=UPOME(1:npelms,NZ)
    HEUPTK(ielmn,NZ)=HEUPTK(ielmn,NZ)+UPNH4(NZ)+UPNO3(NZ)+UPNF(NZ)
    HEUPTK(ielmp,NZ)=HEUPTK(ielmp,NZ)+UPH2P(NZ)+UPH1P(NZ)
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
    FWODLE =>  plt_allom%FWODLE   , &
    FWODBE =>  plt_allom%FWODBE   , &
    FWOODE =>  plt_allom%FWOODE   , &
    FWODRE =>  plt_allom%FWODRE   , &
    CNLF   =>  plt_allom%CNLF     , &
    CPLF   =>  plt_allom%CPLF     , &
    CNSHE  =>  plt_allom%CNSHE    , &
    CPSHE  =>  plt_allom%CPSHE    , &
    CNSTK  =>  plt_allom%CNSTK    , &
    CPSTK  =>  plt_allom%CPSTK    , &
    k_fine_litr=> pltpar%k_fine_litr,&
    k_woody_litr=> pltpar%k_woody_litr,&
    RCS    =>  plt_photo%RCS      , &
    RTN1   =>  plt_morph%RTN1     , &
    RTNL   =>  plt_morph%RTNL     , &
    MY     =>  plt_morph%MY       , &
    ARLFV  =>  plt_morph%ARLFV    , &
    ARSTV  =>  plt_morph%ARSTV    , &
    NRT    =>  plt_morph%NRT        &
  )
  D2: DO L=1,JC1
    ARLFV(L,NZ)=0._r8
    WGLFV(L,NZ)=0._r8
    ARSTV(L,NZ)=0._r8
  ENDDO D2
  D5: DO NR=1,NRT(NZ)
    DO  N=1,MY(NZ)
      NRX(N,NR)=0
      ICHK1(N,NR)=0
    enddo
  ENDDO D5
  D9: DO N=1,MY(NZ)
    D6: DO L=NU,NJ
      WSRTL(N,L,NZ)=0._r8
      RTN1(N,L,NZ)=0._r8
      RTNL(N,L,NZ)=0._r8
      RCO2M(N,L,NZ)=0._r8
      RCO2N(N,L,NZ)=0._r8
      RCO2A(N,L,NZ)=0._r8
    ENDDO D6
  ENDDO D9
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
  IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1.OR.WTSTKE(ielmc,NZ).LE.ZEROP(NZ))THEN
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=1.0_r8
    FWODRE(ielmc,k_fine_litr)=1.0_r8
  ELSE
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=SQRT(WVSTK(NZ)/WTSTKE(ielmc,NZ))
    FWODRE(ielmc,k_fine_litr)=SQRT(FRTX*WVSTK(NZ)/WTSTKE(ielmc,NZ))
  ENDIF
  FWODBE(ielmc,k_woody_litr)=1.0_r8-FWODBE(ielmc,k_fine_litr)
  FWOODE(ielmc,k_woody_litr)=1.0_r8-FWOODE(ielmc,k_fine_litr)
  FWODRE(ielmc,k_woody_litr)=1.0_r8-FWODRE(ielmc,k_fine_litr)
  CNLFW=FWODBE(ielmc,k_woody_litr)*CNSTK(NZ)+FWODBE(ielmc,k_fine_litr)*CNLF(NZ)
  CPLFW=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)+FWODBE(ielmc,k_fine_litr)*CPLF(NZ)
  CNSHW=FWODBE(ielmc,k_woody_litr)*CNSTK(NZ)+FWODBE(ielmc,k_fine_litr)*CNSHE(NZ)
  CPSHW=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)+FWODBE(ielmc,k_fine_litr)*CPSHE(NZ)
  CNRTW=FWODRE(ielmc,k_woody_litr)*CNSTK(NZ)+FWODRE(ielmc,k_fine_litr)*CNRT(NZ)
  CPRTW=FWODRE(ielmc,k_woody_litr)*CPSTK(NZ)+FWODRE(ielmc,k_fine_litr)*CPRT(NZ)
  FWODLE(ielmc,1:n_pltlitrk)=FWODBE(ielmc,1:n_pltlitrk)
  FWODLE(ielmn,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*CNSTK(NZ)/CNLFW
  FWODLE(ielmp,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)/CPLFW
  FWODBE(ielmn,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*CNSTK(NZ)/CNSHW
  FWODBE(ielmp,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)/CPSHW
  FWOODE(ielmn,k_woody_litr)=FWOODE(ielmc,k_woody_litr)*CNSTK(NZ)/CNRTW
  FWOODE(ielmp,k_woody_litr)=FWOODE(ielmc,k_woody_litr)*CPSTK(NZ)/CPRTW
  FWODRE(ielmn,k_woody_litr)=FWODRE(ielmc,k_woody_litr)*CNRT(NZ)/CNRTW
  FWODRE(ielmp,k_woody_litr)=FWODRE(ielmc,k_woody_litr)*CPRT(NZ)/CPRTW
  FWODLE(ielmn,k_fine_litr)=1.0_r8-FWODLE(ielmn,k_woody_litr)
  FWODLE(ielmp,k_fine_litr)=1.0_r8-FWODLE(ielmp,k_woody_litr)
  FWODBE(ielmn,k_fine_litr)=1.0_r8-FWODBE(ielmn,k_woody_litr)
  FWODBE(ielmp,k_fine_litr)=1.0_r8-FWODBE(ielmp,k_woody_litr)
  FWOODE(ielmn,k_fine_litr)=1.0_r8-FWOODE(ielmn,k_woody_litr)
  FWOODE(ielmp,k_fine_litr)=1.0_r8-FWOODE(ielmp,k_woody_litr)
  FWODRE(ielmn,k_fine_litr)=1.0_r8-FWODRE(ielmn,k_woody_litr)
  FWODRE(ielmp,k_fine_litr)=1.0_r8-FWODRE(ielmp,k_woody_litr)
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
  RTK=RGAS*TKCM
  STK=710.0_r8*TKCM
  ACTVM=1._r8+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TFN5=EXP(25.214_r8-62500._r8/RTK)/ACTVM
  D7: DO L=NU,NJ
    TKSM=TKS(L)+OFFST(NZ)
    RTK=RGAS*TKSM
    STK=710.0_r8*TKSM
    ACTVM=1+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
    TFN6(L)=EXP(25.214_r8-62500._r8/RTK)/ACTVM
  ENDDO D7
!
!     PRIMARY ROOT NUMBER
!
!     WTRTA=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     XRTN1=multiplier for number of primary root axes
!
  WTRTA(NZ)=AMAX1(0.999992087_r8*WTRTA(NZ),WTRTE(ielmc,NZ)/PP(NZ))
  XRTN1=AMAX1(1.0_r8,WTRTA(NZ)**0.667_r8)*PP(NZ)
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
  WFNS=AMIN1(1.0_r8,AZMAX1(PSILG(NZ)-PSILM))
  IF(IGTYP(NZ).EQ.0)THEN
    WFNC=1.0_r8
    WFNG=EXP(0.05_r8*PSILT(NZ))
    WFNSG=WFNS**0.10_r8
  ELSE
    WFNC=EXP(RCS(NZ)*PSILG(NZ))
    WFNG=EXP(0.10_r8*PSILT(NZ))
    WFNSG=WFNS**0.25_r8
  ENDIF
  end associate
  end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

  subroutine ComputeTotalBiom(NZ,CPOOLK)

  integer, intent(in) :: NZ
  real(r8), intent(out) :: CPOOLK(JC1,JP1)
  integer :: L,K,N,NB,NE
!     begin_execution
  associate(                                 &
    WTLFBE     =>  plt_biom%WTLFBE     , &
    WTGRBE     =>  plt_biom%WTGRBE     , &
    EPOOLR     =>  plt_biom%EPOOLR     , &
    WTRTD      =>  plt_biom%WTRTD      , &
    WTSHTBE    =>  plt_biom%WTSHTBE    , &
    WTEARBE    =>  plt_biom%WTEARBE    , &
    WTSTKBE    =>  plt_biom%WTSTKBE    , &
    WTSHEBE    =>  plt_biom%WTSHEBE    , &
    WTHSKBE    =>  plt_biom%WTHSKBE    , &
    WTRSVBE    =>  plt_biom%WTRSVBE    , &
    EPOOL      =>  plt_biom%EPOOL      , &
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
!     WTRSVBE,WTRSBN,WTRSBP=stalk reserve C,N,P mass
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
  DO NE=1,npelms
    DO NB=1,NBR(NZ)
      WTSHTBE(NB,NE,NZ)=WTLFBE(NB,NE,NZ) &
        +WTSHEBE(NB,NE,NZ)+WTSTKBE(NB,NE,NZ)+WTRSVBE(NB,NE,NZ) &
        +WTHSKBE(NB,NE,NZ)+WTEARBE(NB,NE,NZ)+WTGRBE(NB,NE,NZ) &
        +EPOOL(NB,NE,NZ)
    ENDDO
  ENDDO

  D320: DO NB=1,NBR(NZ)
    CPOOLK(NB,NZ)=0._r8
    D325: DO K=1,JNODS1
      CPOOLK(NB,NZ)=CPOOLK(NB,NZ)+CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
        +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
    ENDDO D325
    WTSHTBE(NB,ielmc,NZ)=WTSHTBE(NB,ielmc,NZ)+CPOOLK(NB,NZ)
  ENDDO D320

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
      WTRTD(N,L,NZ)=WTRTD(N,L,NZ)+EPOOLR(ielmc,N,L,NZ)
    enddo
  ENDDO D345
  end associate
  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,UPNFC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: UPNFC(JP1)
  integer :: L,NR,N,NB,NE
!     begin_execution
  associate(                            &
    EPOOL    =>  plt_biom%EPOOL   , &
    EPOOLR   =>  plt_biom%EPOOLR  , &
    EPOLNB   =>  plt_biom%EPOLNB  , &
    WTNDE    =>  plt_biom%WTNDE   , &
    WTRTSE   =>  plt_biom%WTRTSE  , &
    WTRTE    =>  plt_biom%WTRTE   , &
    WTNDBE   =>  plt_biom%WTNDBE  , &
    WTSHTBE  =>  plt_biom%WTSHTBE , &
    WTSTKBE  =>  plt_biom%WTSTKBE , &
    WTHSKBE  =>  plt_biom%WTHSKBE , &
    WTRSVBE   =>  plt_biom%WTRSVBE  , &
    WTEARBE  =>  plt_biom%WTEARBE , &
    WTLSB    =>  plt_biom%WTLSB   , &
    WTGRBE   =>  plt_biom%WTGRBE  , &
    WTLFBE   =>  plt_biom%WTLFBE  , &
    WVSTKB   =>  plt_biom%WVSTKB  , &
    WTSHEBE  =>  plt_biom%WTSHEBE , &
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
    WTNDLE   =>  plt_biom%WTNDLE  , &
    EPOOLP   =>  plt_biom%EPOOLP  , &
    EPOLNP   =>  plt_biom%EPOLNP  , &
    WTRT1E   =>  plt_biom%WTRT1E  , &
    WTRT2E   =>  plt_biom%WTRT2E  , &
    EPOOLN   =>  plt_biom%EPOOLN  , &
    TZUPFX   =>  plt_bgcr%TZUPFX  , &
    TEUPTK   =>  plt_rbgc%TEUPTK  , &
    UPH1P    =>  plt_rbgc%UPH1P   , &
    HEUPTK   =>  plt_rbgc%HEUPTK  , &
    UPNF     =>  plt_rbgc%UPNF    , &
    UPH2P    =>  plt_rbgc%UPH2P   , &
    UPNO3    =>  plt_rbgc%UPNO3   , &
    UPNH4    =>  plt_rbgc%UPNH4   , &
    UPOME    =>  plt_rbgc%UPOME   , &
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
!     WTRSVBE,WTRSBN,WTRSBP=stalk reserve C,N,P mass
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
  DO NE=1,npelms
    EPOOLP(NE,NZ)=sum(EPOOL(1:NBR(NZ),NE,NZ))
    WTSHTE(NE,NZ)=sum(WTSHTBE(1:NBR(NZ),NE,NZ))
    WTSHEE(NE,NZ)=sum(WTSHEBE(1:NBR(NZ),NE,NZ))
    WTSTKE(NE,NZ)=sum(WTSTKBE(1:NBR(NZ),NE,NZ))
    WTLFE(NE,NZ)=sum(WTLFBE(1:NBR(NZ),NE,NZ))
    WTRSVE(NE,NZ)=sum(WTRSVBE(1:NBR(NZ),NE,NZ))
    WTHSKE(NE,NZ)=sum(WTHSKBE(1:NBR(NZ),NE,NZ))
    WTGRE(NE,NZ)=sum(WTGRBE(1:NBR(NZ),NE,NZ))
    WTEARE(NE,NZ)=sum(WTEARBE(1:NBR(NZ),NE,NZ))
!root state variables
    WTRTE(NE,NZ)=sum(EPOOLR(NE,1:MY(NZ),NU:NJ,NZ))
    WTRTSE(NE,NZ)=sum(WTRT1E(NE,1:MY(NZ),NU:NJ,1:NRT(NZ),NZ)) &
      +sum(WTRT2E(NE,1:MY(NZ),NU:NJ,1:NRT(NZ),NZ))
    WTRTE(NE,NZ)=WTRTE(NE,NZ)+WTRTSE(NE,NZ)
  ENDDO

  WVSTK(NZ)=sum(WVSTKB(1:NBR(NZ),NZ))
  WTLS(NZ) =sum(WTLSB(1:NBR(NZ),NZ))
  GRNO(NZ) =sum(GRNOB(1:NBR(NZ),NZ))
  ARLFP(NZ)=sum(ARLFB(1:NBR(NZ),NZ))
  ARSTP(NZ)=sum(ARSTK(1:JC1,1:NBR(NZ),NZ))
  ARSTV(1:JC1,1:NBR(NZ))=0._r8
  DO NB=1,NBR(NZ)
    DO L=1,JC1
      ARSTV(L,NZ)=ARSTV(L,NZ)+ARSTK(L,NB,NZ)
    ENDDO
  ENDDO

!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  IF(INTYP(NZ).NE.0)THEN
    IF(INTYP(NZ).GE.4)THEN
      DO NE=1,npelms
        D7950: DO NB=1,NBR(NZ)
          EPOLNP(NE,NZ)=EPOLNP(NE,NZ)+EPOLNB(NB,NE,NZ)
        ENDDO D7950
        WTNDE(NE,NZ)=sum(WTNDBE(1:NBR(NZ),NE,NZ))+sum(EPOLNB(1:NBR(NZ),NE,NZ))
      ENDDO
    ELSEIF(INTYP(NZ).GE.1.AND.INTYP(NZ).LE.3)THEN
      DO NE=1,npelms
        WTNDE(NE,NZ)=sum(WTNDLE(NU:NI(NZ),NE,NZ))+sum(EPOOLN(NU:NI(NZ),NE,NZ))
      ENDDO
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
  HEUPTK(1:npelms,NZ)=UPOME(1:npelms,NZ)
  HEUPTK(ielmn,NZ)=HEUPTK(ielmn,NZ)+UPNH4(NZ)+UPNO3(NZ)+UPNF(NZ)
  HEUPTK(ielmp,NZ)=HEUPTK(ielmp,NZ)+UPH2P(NZ)+UPH1P(NZ)

  TEUPTK(1:npelms,NZ)=TEUPTK(1:npelms,NZ)+UPOME(1:npelms,NZ)
  TEUPTK(ielmn,NZ)=TEUPTK(ielmn,NZ)+UPNH4(NZ)+UPNO3(NZ)
  TEUPTK(ielmp,NZ)=TEUPTK(ielmp,NZ)+UPH2P(NZ)+UPH1P(NZ)
  TZUPFX(NZ)=TZUPFX(NZ)+UPNF(NZ)+UPNFC(NZ)
  end associate
  end subroutine AccumulateStates

end module grosubsMod
