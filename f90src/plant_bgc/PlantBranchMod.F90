module PlantBranchMod

! Description:
! module for plant biological transformations
  use minimathmod, only : test_aeqb,safe_adb,AZMAX1
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
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: GrowOneBranch
  contains
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
  associate(                            &
    DMLF       =>  plt_allom%DMLF     , &
    DMSHE      =>  plt_allom%DMSHE    , &
    DMEAR      =>  plt_allom%DMEAR    , &
    CNEAR      =>  plt_allom%CNEAR    , &
    CPRSV      =>  plt_allom%CPRSV    , &
    CPHSK      =>  plt_allom%CPHSK    , &
    CPSTK      =>  plt_allom%CPSTK    , &
    CNSTK      =>  plt_allom%CNSTK    , &
    DMSTK      =>  plt_allom%DMSTK    , &
    FWODSN     =>  plt_allom%FWODSN   , &
    FWODSP     =>  plt_allom%FWODSP   , &
    DMHSK      =>  plt_allom%DMHSK    , &
    FWODB      =>  plt_allom%FWODB    , &
    FWODLN     =>  plt_allom%FWODLN   , &
    FWODLP     =>  plt_allom%FWODLP   , &
    CNHSK      =>  plt_allom%CNHSK    , &
    CPEAR      =>  plt_allom%CPEAR    , &
    DMRSV      =>  plt_allom%DMRSV    , &
    DMGR       =>  plt_allom%DMGR     , &
    CNWS       =>  plt_allom%CNWS     , &
    CPWS       =>  plt_allom%CPWS     , &
    CNRSV      =>  plt_allom%CNRSV    , &
    DMRT       =>  plt_allom%DMRT     , &
    FNOD       =>  plt_allom%FNOD     , &
    ESNC       =>  plt_bgcr%ESNC      , &
    CFOPC      =>  plt_soilchem%CFOPC , &
    CFOPN      =>  plt_soilchem%CFOPN , &
    CFOPP      =>  plt_soilchem%CFOPP , &
    WTSTKBE    =>  plt_biom%WTSTKBE   , &
    WTRSVBE    =>  plt_biom%WTRSVBE   , &
    WTHSKBE    =>  plt_biom%WTHSKBE   , &
    WTEARBE    =>  plt_biom%WTEARBE   , &
    WTSTXB     =>  plt_biom%WTSTXB    , &
    WTSTXN     =>  plt_biom%WTSTXN    , &
    WTSTXP     =>  plt_biom%WTSTXP    , &
    WGSHEXE    =>  plt_biom%WGSHEXE   , &
    WVSTKB     =>  plt_biom%WVSTKB    , &
    EPOOL      =>  plt_biom%EPOOL     , &
    WTLSB      =>  plt_biom%WTLSB     , &
    WGNODN     =>  plt_biom%WGNODN    , &
    WGLFX      =>  plt_biom%WGLFX     , &
    WGLFPX     =>  plt_biom%WGLFPX    , &
    WGLFNX     =>  plt_biom%WGLFNX    , &
    WTSHEBE    =>  plt_biom%WTSHEBE   , &
    WTLFBE     =>  plt_biom%WTLFBE    , &
    WSLF       =>  plt_biom%WSLF      , &
    WGLFE      =>  plt_biom%WGLFE     , &
    WGSHE      =>  plt_biom%WGSHE     , &
    CEPOLB     =>  plt_biom%CEPOLB    , &
    WGNODE     =>  plt_biom%WGNODE    , &
    WGNODP     =>  plt_biom%WGNODP    , &
    WSSHE      =>  plt_biom%WSSHE     , &
    WTGRBE     =>  plt_biom%WTGRBE    , &
    ZEROP      =>  plt_biom%ZEROP     , &
    ZEROL      =>  plt_biom%ZEROL     , &
    ICTYP      =>  plt_photo%ICTYP    , &
    IDTHB      =>  plt_pheno%IDTHB    , &
    KVSTGN     =>  plt_pheno%KVSTGN   , &
    XRLA       =>  plt_pheno%XRLA     , &
    RCCLX      =>  plt_pheno%RCCLX    , &
    RCZLX      =>  plt_pheno%RCZLX    , &
    RCPLX      =>  plt_pheno%RCPLX    , &
    KVSTG      =>  plt_pheno%KVSTG    , &
    IGTYP      =>  plt_pheno%IGTYP    , &
    TFN3       =>  plt_pheno%TFN3     , &
    IFLGP      =>  plt_pheno%IFLGP    , &
    FLGZ       =>  plt_pheno%FLGZ     , &
    IDAY       =>  plt_pheno%IDAY     , &
    RCCSX      =>  plt_pheno%RCCSX    , &
    RCZSX      =>  plt_pheno%RCZSX    , &
    RCPSX      =>  plt_pheno%RCPSX    , &
    IBTYP      =>  plt_pheno%IBTYP    , &
    IFLGG      =>  plt_pheno%IFLGG    , &
    ISTYP      =>  plt_pheno%ISTYP    , &
    CTC        =>  plt_pheno%CTC      , &
    SSINN      =>  plt_rad%SSINN      , &
    PP         =>  plt_site%PP        , &
    ZERO       =>  plt_site%ZERO      , &
    PSILT      =>  plt_ew%PSILT       , &
    NB1        =>  plt_morph%NB1      , &
    HTCTL      =>   plt_morph%HTCTL   , &
    NBTB       =>   plt_morph%NBTB    , &
    NBR        =>   plt_morph%NBR     , &
    SLA1       =>   plt_morph%SLA1    , &
    SNL1       =>   plt_morph%SNL1    , &
    ANGSH      =>   plt_morph%ANGSH   , &
    ARLFZ      =>   plt_morph%ARLFZ   , &
    ARLFB      =>   plt_morph%ARLFB   , &
    SDPTH      =>   plt_morph%SDPTH   , &
    HTSHEX     =>   plt_morph%HTSHEX  , &
    ARLF1      =>   plt_morph%ARLF1   , &
    HTSHE      =>   plt_morph%HTSHE   , &
    ANGBR      =>   plt_morph%ANGBR   , &
    SSL1       =>   plt_morph%SSL1    , &
    HTNODE     =>   plt_morph%HTNODE  , &
    HTNODX     =>   plt_morph%HTNODX  , &
    NNOD       =>   plt_morph%NNOD      &
  )
  WTLSB(NB,NZ)=AZMAX1(WTLFBE(NB,ielmc,NZ)+WTSHEBE(NB,ielmc,NZ))

  IF(IDTHB(NB,NZ).EQ.0)THEN
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
    IF(IDAY(1,NB,NZ).NE.0)THEN
      DMLFB=DMLF(NZ)
      DMSHB=DMSHE(NZ)
      CNLFB=CNLFW
      CNSHB=CNSHW
      CPLFB=CPLFW
      CPSHB=CPSHW
    ELSE
      DMLFB=DMRT(NZ)
      DMSHB=DMRT(NZ)
      CNLFB=CNRTW
      CNSHB=CNRTW
      CPLFB=CPRTW
      CPSHB=CPRTW
    ENDIF
    DMSHT=PART(1)*DMLFB+PART(2)*DMSHB+PART(3)*DMSTK(NZ) &
      +PART(4)*DMRSV(NZ)+PART(5)*DMHSK(NZ) &
      +PART(6)*DMEAR(NZ)+PART(7)*DMGR(NZ)
    DMSHD=1.0_r8-DMSHT
    CNLFM=PART(1)*DMLFB*ZPLFM*CNLFB
    CPLFM=PART(1)*DMLFB*ZPLFM*CPLFB
    CNLFX=PART(1)*DMLFB*ZPLFD*CNLFB
    CPLFX=PART(1)*DMLFB*ZPLFD*CPLFB
    CNSHX=PART(2)*DMSHB*CNSHB &
      +PART(3)*DMSTK(NZ)*CNSTK(NZ) &
      +PART(4)*DMRSV(NZ)*CNRSV(NZ) &
      +PART(5)*DMHSK(NZ)*CNHSK(NZ) &
      +PART(6)*DMEAR(NZ)*CNEAR(NZ) &
      +PART(7)*DMGR(NZ)*CNRSV(NZ)
    CPSHX=PART(2)*DMSHB*CPSHB &
      +PART(3)*DMSTK(NZ)*CPSTK(NZ) &
      +PART(4)*DMRSV(NZ)*CPRSV(NZ) &
      +PART(5)*DMHSK(NZ)*CPHSK(NZ) &
      +PART(6)*DMEAR(NZ)*CPEAR(NZ) &
      +PART(7)*DMGR(NZ)*CPRSV(NZ)
!
!   TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
!
!   WTSHXN=shoot structural N mass
!   WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,petiole,husk,ear,grain N mass
!   CNSTK,WVSTKB=stalk N:C,sapwood mass
!   IDAY(10=date of physiological maturity
!
    WTSHXN=AZMAX1(WTLFBE(NB,ielmn,NZ)+WTSHEBE(NB,ielmn,NZ) &
      +CNSTK(NZ)*WVSTKB(NB,NZ))
    IF(IDAY(10,NB,NZ).EQ.0)THEN
      WTSHXN=WTSHXN+AZMAX1(WTHSKBE(NB,ielmn,NZ) &
        +WTEARBE(NB,ielmn,NZ)+WTGRBE(NB,ielmn,NZ))
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
    IF(ICTYP(NZ).EQ.4)THEN
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
    GROSTK=PART(3)*CGROS*DMSTK(NZ)
    GRORSV=PART(4)*CGROS*DMRSV(NZ)
    GROHSK=PART(5)*CGROS*DMHSK(NZ)
    GROEAR=PART(6)*CGROS*DMEAR(NZ)
    GROGR=PART(7)*CGROS*DMGR(NZ)
    GROSHT=CGROS*DMSHT
    GROLFN=GROLF*CNLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHN=GROSHE*CNSHB
    GROSTN=GROSTK*CNSTK(NZ)
    GRORSN=GRORSV*CNRSV(NZ)
    GROHSN=GROHSK*CNHSK(NZ)
    GROEAN=GROEAR*CNEAR(NZ)
    GROGRN=GROGR*CNRSV(NZ)
    GROLFP=GROLF*CPLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHP=GROSHE*CPSHB
    GROSTP=GROSTK*CPSTK(NZ)
    GRORSP=GRORSV*CPRSV(NZ)
    GROHSP=GROHSK*CPHSK(NZ)
    GROEAP=GROEAR*CPEAR(NZ)
    GROGRP=GROGR*CPRSV(NZ)
    WTLFBE(NB,ielmc,NZ)=WTLFBE(NB,ielmc,NZ)+GROLF
    WTSHEBE(NB,ielmc,NZ)=WTSHEBE(NB,ielmc,NZ)+GROSHE
    WTSTKBE(NB,ielmc,NZ)=WTSTKBE(NB,ielmc,NZ)+GROSTK
    WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+GRORSV
    WTHSKBE(NB,ielmc,NZ)=WTHSKBE(NB,ielmc,NZ)+GROHSK
    WTEARBE(NB,ielmc,NZ)=WTEARBE(NB,ielmc,NZ)+GROEAR
    WTLFBE(NB,ielmn,NZ)=WTLFBE(NB,ielmn,NZ)+GROLFN
    WTSHEBE(NB,ielmn,NZ)=WTSHEBE(NB,ielmn,NZ)+GROSHN
    WTSTKBE(NB,ielmn,NZ)=WTSTKBE(NB,ielmn,NZ)+GROSTN
    WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)+GRORSN
    WTHSKBE(NB,ielmn,NZ)=WTHSKBE(NB,ielmn,NZ)+GROHSN
    WTEARBE(NB,ielmn,NZ)=WTEARBE(NB,ielmn,NZ)+GROEAN
    WTLFBE(NB,ielmp,NZ)=WTLFBE(NB,ielmp,NZ)+GROLFP
    WTSHEBE(NB,ielmp,NZ)=WTSHEBE(NB,ielmp,NZ)+GROSHP
    WTSTKBE(NB,ielmp,NZ)=WTSTKBE(NB,ielmp,NZ)+GROSTP
    WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)+GRORSP
    WTHSKBE(NB,ielmp,NZ)=WTHSKBE(NB,ielmp,NZ)+GROHSP
    WTEARBE(NB,ielmp,NZ)=WTEARBE(NB,ielmp,NZ)+GROEAP
!
!   ETOLIATION
!
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
!   ETOL=coefficient for etoliation effects on expansion,extension
!
    CCE=AMIN1(safe_adb(CEPOLB(NB,ielmn,NZ),CEPOLB(NB,ielmn,NZ)+CEPOLB(NB,ielmc,NZ)*CNKI) &
      ,safe_adb(CEPOLB(NB,ielmp,NZ),CEPOLB(NB,ielmp,NZ)+CEPOLB(NB,ielmc,NZ)*CPKI))

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
    IF(NB.EQ.NB1(NZ).AND.HTCTL(NZ).LE.SDPTH(NZ))THEN
      NNOD1=0
    ELSE
      NNOD1=1
    ENDIF
    IF(GROLF.GT.0.0)THEN
      MXNOD=KVSTG(NB,NZ)
      MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ)+1)
      MXNOD=MAX(MXNOD,MNNOD)
      KNOD=MXNOD-MNNOD+1
      GNOD=KNOD
      ALLOCL=1.0_r8/GNOD
      GRO=ALLOCL*GROLF
      GRON=ALLOCL*GROLFN
      GROP=ALLOCL*GROLFP
      GSLA=ALLOCL*FNOD(NZ)*NNOD(NZ)
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
          WGLFE(K,NB,ielmc,NZ)=WGLFE(K,NB,ielmc,NZ)+GRO
          WGLFE(K,NB,ielmn,NZ)=WGLFE(K,NB,ielmn,NZ)+GRON
          WGLFE(K,NB,ielmp,NZ)=WGLFE(K,NB,ielmp,NZ)+GROP
          WSLF(K,NB,NZ)=WSLF(K,NB,NZ) &
            +AMIN1(GRON*CNWS(NZ),GROP*CPWS(NZ))
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
          SLA=ETOL*SLA1(NZ)*(AMAX1(ZEROL(NZ) &
            ,WGLFE(K,NB,ielmc,NZ))/(PP(NZ)*GSLA))**SLA2*WFNS
          GROA=GRO*SLA
          ARLFB(NB,NZ)=ARLFB(NB,NZ)+GROA
          ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)+GROA
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
        MXNOD=KVSTG(NB,NZ)
        MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ)+1)
        MXNOD=MAX(MXNOD,MNNOD)
        GNOD=MXNOD-MNNOD+1
        ALLOCS=1.0_r8/GNOD
        GRO=ALLOCS*GROSHE
        GRON=ALLOCS*GROSHN
        GROP=ALLOCS*GROSHP
        GSSL=ALLOCL*FNOD(NZ)*NNOD(NZ)
!
!       GROWTH AT EACH CURRENT NODE
!
!       WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!       GRO,GRON,GROP=petiole C,N,P growth at each node
!       CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
        D505: DO KK=MNNOD,MXNOD
          K=MOD(KK,JNODS1)
          IF(K.EQ.0.AND.KK.NE.0)K=25
            WGSHE(K,NB,ielmc,NZ)=WGSHE(K,NB,ielmc,NZ)+GRO
            WGSHE(K,NB,ielmn,NZ)=WGSHE(K,NB,ielmn,NZ)+GRON
            WGSHE(K,NB,ielmp,NZ)=WGSHE(K,NB,ielmp,NZ)+GROP
            WSSHE(K,NB,NZ)=WSSHE(K,NB,NZ) &
              +AMIN1(GRON*CNWS(NZ),GROP*CPWS(NZ))
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
            IF(WGLFE(K,NB,ielmc,NZ).GT.0.0_r8)THEN
              SSL=ETOL*SSL1(NZ)*(AMAX1(ZEROL(NZ) &
                ,WGSHE(K,NB,ielmc,NZ))/(PP(NZ)*GSSL))**SSL2*WFNS
              GROS=GRO/PP(NZ)*SSL
              HTSHE(K,NB,NZ)=HTSHE(K,NB,NZ)+GROS*ANGSH(NZ)
            ENDIF
          ENDDO D505
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
        IF(IDAY(1,NB,NZ).EQ.0)THEN
          NN=0
        ELSE
          NN=1
        ENDIF
        MXNOD=KVSTG(NB,NZ)
        MNNOD=MAX(MIN(NN,MAX(NN,MXNOD-NNOD(NZ))),KVSTG(NB,NZ)-23)
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
          SNL=ETOL*SNL1(NZ)*(WTSTKBE(NB,ielmc,NZ)/PP(NZ))**SNL2
          GROH=GRO/PP(NZ)*SNL
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
            WGNODE(K1,NB,NZ)=WGNODE(K1,NB,NZ)+GRO
            WGNODN(K1,NB,NZ)=WGNODN(K1,NB,NZ)+GRON
            WGNODP(K1,NB,NZ)=WGNODP(K1,NB,NZ)+GROP
            HTNODX(K1,NB,NZ)=HTNODX(K1,NB,NZ)+GROH*ANGBR(NZ)
            IF(K1.NE.0)THEN
              HTNODE(K1,NB,NZ)=HTNODX(K1,NB,NZ)+HTNODE(K2,NB,NZ)
            ELSE
              HTNODE(K1,NB,NZ)=HTNODX(K1,NB,NZ)
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
        IF(IDAY(1,NB,NZ).NE.0 &
         .AND.CEPOLB(NB,ielmc,NZ).GT.ZERO)THEN
          CCC=AZMAX1(AMIN1(1.0 &
            ,safe_adb(CEPOLB(NB,ielmn,NZ),CEPOLB(NB,ielmn,NZ) &
            +CEPOLB(NB,ielmc,NZ)*CNKI) &
            ,safe_adb(CEPOLB(NB,ielmp,NZ),CEPOLB(NB,ielmp,NZ) &
            +CEPOLB(NB,ielmc,NZ)*CPKI)))
          CNC=AZMAX1(AMIN1(1.0 &
            ,safe_adb(CEPOLB(NB,ielmc,NZ),CEPOLB(NB,ielmc,NZ) &
            +CEPOLB(NB,ielmn,NZ)/CNKI)))
          CPC=AZMAX1(AMIN1(1.0 &
            ,safe_adb(CEPOLB(NB,ielmc,NZ),CEPOLB(NB,ielmc,NZ) &
            +CEPOLB(NB,ielmp,NZ)/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        RCCC=RCCZ(IGTYP(NZ))+CCC*RCCY(IGTYP(NZ))
        RCCN=CNC*RCCX(IGTYP(NZ))
        RCCP=CPC*RCCQ(IGTYP(NZ))
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
        IF(IFLGG(NB,NZ).EQ.1)THEN
          KVSTGX=MAX(0,KVSTG(NB,NZ)-JNODS1+1)
          IF(KVSTGX.GT.0)THEN
            K=MOD(KVSTGX,JNODS1)
            IF(K.EQ.0.AND.KVSTGX.GT.0)K=JNODS1
            KX=MOD(KVSTG(NB,NZ),JNODS1)
            IF(KX.EQ.0.AND.KVSTG(NB,NZ).NE.0)KX=JNODS1
            FSNC=TFN3(NZ)*XRLA(NZ)
!
        !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
        !
        !   IFLGP=flag for remobilization
        !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !   ARLF=node leaf area
        !   RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
!
            IF(IFLGP(NB,NZ).EQ.1)THEN
              WGLFX(NB,NZ)=AZMAX1(WGLFE(K,NB,ielmc,NZ))
              WGLFNX(NB,NZ)=AZMAX1(WGLFE(K,NB,ielmn,NZ))
              WGLFPX(NB,NZ)=AZMAX1(WGLFE(K,NB,ielmp,NZ))
              ARLFZ(NB,NZ)=AZMAX1(ARLF1(K,NB,NZ))
              IF(WGLFX(NB,NZ).GT.ZEROP(NZ))THEN
                RCCLX(NB,NZ)=RCCC*WGLFX(NB,NZ)
                RCZLX(NB,NZ)=WGLFNX(NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
                RCPLX(NB,NZ)=WGLFPX(NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
              ELSE
                RCCLX(NB,NZ)=0._r8
                RCZLX(NB,NZ)=0._r8
                RCPLX(NB,NZ)=0._r8
              ENDIF
            ENDIF
!
    !       FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !       FSNC,FSNCL=fraction of lowest leaf to be remobilized
    !
            IF(FSNC*WGLFX(NB,NZ).GT.WGLFE(K,NB,ielmc,NZ) &
              .AND.WGLFX(NB,NZ).GT.ZEROP(NZ))THEN
              FSNCL=AZMAX1(WGLFE(K,NB,ielmc,NZ)/WGLFX(NB,NZ))
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
            D6300: DO M=1,jsken
              ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(5,M,NZ) &
                *FSNCL*(WGLFX(NB,NZ)-RCCLX(NB,NZ))*FWODB(0)
              ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(5,M,NZ) &
                *FSNCL*(WGLFNX(NB,NZ)-RCZLX(NB,NZ))*FWODLN(0)
              ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(5,M,NZ) &
                *FSNCL*(WGLFPX(NB,NZ)-RCPLX(NB,NZ))*FWODLP(0)
              ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(1,M,NZ) &
                *FSNCL*(WGLFX(NB,NZ)-RCCLX(NB,NZ))*FWODB(1)
              ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(1,M,NZ) &
                *FSNCL*(WGLFNX(NB,NZ)-RCZLX(NB,NZ))*FWODLN(1)
              ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(1,M,NZ) &
                *FSNCL*(WGLFPX(NB,NZ)-RCPLX(NB,NZ))*FWODLP(1)
            ENDDO D6300
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
            ARLFB(NB,NZ)=ARLFB(NB,NZ)-FSNCL*ARLFZ(NB,NZ)
            WTLFBE(NB,ielmc,NZ)=WTLFBE(NB,ielmc,NZ)-FSNCL*WGLFX(NB,NZ)
            WTLFBE(NB,ielmn,NZ)=WTLFBE(NB,ielmn,NZ)-FSNCL*WGLFNX(NB,NZ)
            WTLFBE(NB,ielmp,NZ)=WTLFBE(NB,ielmp,NZ)-FSNCL*WGLFPX(NB,NZ)
            ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)-FSNCL*ARLFZ(NB,NZ)
            WGLFE(K,NB,ielmc,NZ)=WGLFE(K,NB,ielmc,NZ)-FSNCL*WGLFX(NB,NZ)
            WGLFE(K,NB,ielmn,NZ)=WGLFE(K,NB,ielmn,NZ)-FSNCL*WGLFNX(NB,NZ)
            WGLFE(K,NB,ielmp,NZ)=WGLFE(K,NB,ielmp,NZ)-FSNCL*WGLFPX(NB,NZ)
            WSLF(K,NB,NZ)=AZMAX1(WSLF(K,NB,NZ) &
              -FSNCL*AMAX1(WGLFNX(NB,NZ)*CNWS(NZ) &
              ,WGLFPX(NB,NZ)*CPWS(NZ)))
            EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+FSNCL*RCCLX(NB,NZ)
            EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+FSNCL*RCZLX(NB,NZ)
            EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+FSNCL*RCPLX(NB,NZ)
!
    !       REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P ALSO DEPENDS ON
    !       STRUCTURAL C:N:P
    !
    !       IFLGP=flag for remobilization
    !       WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !       HTSHEX=petiole length
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
!
            IF(IFLGP(NB,NZ).EQ.1)THEN
              WGSHEXE(NB,ielmc,NZ)=AZMAX1(WGSHE(K,NB,ielmc,NZ))
              WGSHEXE(NB,ielmn,NZ)=AZMAX1(WGSHE(K,NB,ielmn,NZ))
              WGSHEXE(NB,ielmp,NZ)=AZMAX1(WGSHE(K,NB,ielmp,NZ))
              HTSHEX(NB,NZ)=AZMAX1(HTSHE(K,NB,NZ))
              IF(WGSHEXE(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
                RCCSX(NB,NZ)=RCCC*WGSHEXE(NB,ielmc,NZ)
                RCZSX(NB,NZ)=WGSHEXE(NB,ielmn,NZ) &
                  *(RCCN+(1.0_r8-RCCN)*RCCSX(NB,NZ)/WGSHEXE(NB,ielmc,NZ))
                RCPSX(NB,NZ)=WGSHEXE(NB,ielmp,NZ) &
                  *(RCCP+(1.0_r8-RCCP)*RCCSX(NB,NZ)/WGSHEXE(NB,ielmc,NZ))
              ELSE
                RCCSX(NB,NZ)=0._r8
                RCZSX(NB,NZ)=0._r8
                RCPSX(NB,NZ)=0._r8
              ENDIF
              WTSTXB(NB,NZ)=WTSTXB(NB,NZ)+WGNODE(K,NB,NZ)
              WTSTXN(NB,NZ)=WTSTXN(NB,NZ)+WGNODN(K,NB,NZ)
              WTSTXP(NB,NZ)=WTSTXP(NB,NZ)+WGNODP(K,NB,NZ)
              WGNODE(K,NB,NZ)=0._r8
              WGNODN(K,NB,NZ)=0._r8
              WGNODP(K,NB,NZ)=0._r8
              HTNODX(K,NB,NZ)=0._r8
            ENDIF
!
    !       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !
            IF(FSNC*WGSHEXE(NB,ielmc,NZ).GT.WGSHE(K,NB,ielmc,NZ) &
              .AND.WGSHEXE(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
              FSNCS=AZMAX1(WGSHE(K,NB,ielmc,NZ)/WGSHEXE(NB,ielmc,NZ))
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
            D6305: DO M=1,jsken
              ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(5,M,NZ) &
                *FSNCS*(WGSHEXE(NB,ielmc,NZ)-RCCSX(NB,NZ))*FWODB(0)
              ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(5,M,NZ) &
                *FSNCS*(WGSHEXE(NB,ielmn,NZ)-RCZSX(NB,NZ))*FWODSN(0)
              ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(5,M,NZ) &
                *FSNCS*(WGSHEXE(NB,ielmp,NZ)-RCPSX(NB,NZ))*FWODSP(0)
              ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(2,M,NZ) &
                *FSNCS*(WGSHEXE(NB,ielmc,NZ)-RCCSX(NB,NZ))*FWODB(1)
              ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(2,M,NZ) &
                *FSNCS*(WGSHEXE(NB,ielmn,NZ)-RCZSX(NB,NZ))*FWODSN(1)
              ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(2,M,NZ) &
                *FSNCS*(WGSHEXE(NB,ielmp,NZ)-RCPSX(NB,NZ))*FWODSP(1)
            ENDDO D6305
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
            WTSHEBE(NB,ielmc,NZ)=WTSHEBE(NB,ielmc,NZ)-FSNCS*WGSHEXE(NB,ielmc,NZ)
            WTSHEBE(NB,ielmn,NZ)=WTSHEBE(NB,ielmn,NZ)-FSNCS*WGSHEXE(NB,ielmn,NZ)
            WTSHEBE(NB,ielmp,NZ)=WTSHEBE(NB,ielmp,NZ)-FSNCS*WGSHEXE(NB,ielmp,NZ)
            HTSHE(K,NB,NZ)=HTSHE(K,NB,NZ)-FSNCS*HTSHEX(NB,NZ)
            WGSHE(K,NB,ielmc,NZ)=WGSHE(K,NB,ielmc,NZ)-FSNCS*WGSHEXE(NB,ielmc,NZ)
            WGSHE(K,NB,ielmn,NZ)=WGSHE(K,NB,ielmn,NZ)-FSNCS*WGSHEXE(NB,ielmn,NZ)
            WGSHE(K,NB,ielmp,NZ)=WGSHE(K,NB,ielmp,NZ)-FSNCS*WGSHEXE(NB,ielmp,NZ)
            WSSHE(K,NB,NZ)=AZMAX1(WSSHE(K,NB,NZ) &
              -FSNCS*AMAX1(WGSHEXE(NB,ielmn,NZ)*CNWS(NZ) &
              ,WGSHEXE(NB,ielmp,NZ)*CPWS(NZ)))
            EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+FSNCS*RCCSX(NB,NZ)
            EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+FSNCS*RCZSX(NB,NZ)
            EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+FSNCS*RCPSX(NB,NZ)
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
          IF(SNCR.GT.0.0.AND.WTRSVBE(NB,ielmc,NZ).GT.0.0)THEN
            RCO2V=AMIN1(SNCR,VMXC*WTRSVBE(NB,ielmc,NZ)*TFN3(NZ))
            WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)-RCO2V
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
        IF(IFLGZ.EQ.1.AND.IFLGY.EQ.1.AND.ISTYP(NZ).NE.0)THEN
          SNCZ=FXFB(IBTYP(NZ)) &
            *WTLSB(NB,NZ)*AMIN1(1.0,FLGZ(NB,NZ)/FLGZX)
        ELSE
          SNCZ=0._r8
        ENDIF
        SNCX=SNCR+SNCZ
        IF(SNCX.GT.ZEROP(NZ))THEN
          SNCF=SNCZ/SNCX
          KSNC=INT(0.5*(KVSTG(NB,NZ)-KVSTGN(NB,NZ)))+1
          XKSNC=KSNC
          KN=MAX(0,KVSTGN(NB,NZ)-1)
    !
    !     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
    !     IF MAIN STEM POOLS ARE DEPLETED
    !
    !     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
    !     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
    !
          IF(IBTYP(NZ).NE.0.AND.IGTYP(NZ).GT.1 &
            .AND.NB.EQ.NB1(NZ).AND.test_aeqb(SNCF,0._r8))THEN
            NBY=0
          DO 584 NBL=1,NBR(NZ)
            NBZ(NBL)=0
584       CONTINUE
          DO 586 NBL=1,NBR(NZ)
            NBX=KVSTG(NB,NZ)
            DO 585 NBK=1,NBR(NZ)
              IF(IDTHB(NBK,NZ).EQ.0.AND.NBK.NE.NB1(NZ) &
                .AND.NBTB(NBK,NZ).LT.NBX &
                .AND.NBTB(NBK,NZ).GT.NBY)THEN
                NBZ(NBL)=NBK
                NBX=NBTB(NBK,NZ)
              ENDIF
585         CONTINUE
            IF(NBZ(NBL).NE.0)THEN
              NBY=NBTB(NBZ(NBL),NZ)
            ENDIF
586       CONTINUE
          D580: DO NBL=1,NBR(NZ)
            IF(NBZ(NBL).NE.0)THEN
              IF(NBTB(NBZ(NBL),NZ).LT.KK)THEN
                IF(EPOOL(NBZ(NBL),ielmc,NZ).GT.0)THEN
                  XFRC=1.0E-02_r8*AMIN1(SNCX,EPOOL(NBZ(NBL),ielmc,NZ))
                  XFRN=XFRC*EPOOL(NBZ(NBL),ielmn,NZ)/EPOOL(NBZ(NBL),ielmc,NZ)
                  XFRP=XFRC*EPOOL(NBZ(NBL),ielmp,NZ)/EPOOL(NBZ(NBL),ielmc,NZ)
                ELSE
                  XFRC=0._r8
                  XFRN=1.0E-02_r8*EPOOL(NBZ(NBL),ielmn,NZ)
                  XFRP=1.0E-02_r8*EPOOL(NBZ(NBL),ielmp,NZ)
                ENDIF
                EPOOL(NBZ(NBL),ielmc,NZ)=EPOOL(NBZ(NBL),ielmc,NZ)-XFRC
                EPOOL(NBZ(NBL),ielmn,NZ)=EPOOL(NBZ(NBL),ielmn,NZ)-XFRN
                EPOOL(NBZ(NBL),ielmp,NZ)=EPOOL(NBZ(NBL),ielmp,NZ)-XFRP
                EPOOL(NB1(NZ),ielmc,NZ)=EPOOL(NB1(NZ),ielmc,NZ)+XFRC*SNCF
                EPOOL(NB1(NZ),ielmn,NZ)=EPOOL(NB1(NZ),ielmn,NZ)+XFRN
                EPOOL(NB1(NZ),ielmp,NZ)=EPOOL(NB1(NZ),ielmp,NZ)+XFRP
                SNCX=SNCX-XFRC
                IF(SNCX.LE.0.0_r8)exit
              ENDIF
            ENDIF
          ENDDO D580
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
    IF(IBTYP(NZ).NE.0.AND.IGTYP(NZ).GT.1 &
      .AND.IDTHB(NB1(NZ),NZ).EQ.1)IDTHB(NB,NZ)=1
!
!     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
!
    KVSTGX=MAX(0,KVSTG(NB,NZ)-24)
    D495: DO KK=KVSTGX,KVSTG(NB,NZ)
      K=MOD(KK,JNODS1)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      IF(WGLFE(K,NB,ielmc,NZ).GT.0.0)THEN
        CPOOLT=WGLFE(K,NB,ielmc,NZ)+EPOOL(NB,ielmc,NZ)
        IF(CPOOLT.GT.ZEROP(NZ))THEN
          ZPOOLD=WGLFE(K,NB,ielmn,NZ)*EPOOL(NB,ielmc,NZ) &
            -EPOOL(NB,ielmn,NZ)*WGLFE(K,NB,ielmc,NZ)
          XFRN1=AZMAX1(AMIN1(1.0E-03*ZPOOLD/CPOOLT,WGLFE(K,NB,ielmn,NZ) &
            -ZPLFM*CNLFB*WGLFE(K,NB,ielmc,NZ)))
          PPOOLD=WGLFE(K,NB,ielmp,NZ)*EPOOL(NB,ielmc,NZ) &
            -EPOOL(NB,ielmp,NZ)*WGLFE(K,NB,ielmc,NZ)
          XFRP1=AZMAX1(AMIN1(1.0E-03*PPOOLD/CPOOLT,WGLFE(K,NB,ielmp,NZ) &
            -ZPLFM*CPLFB*WGLFE(K,NB,ielmc,NZ)))
          XFRN=AMAX1(XFRN1,10.0*XFRP1)
          XFRP=AMAX1(XFRP1,0.10*XFRN1)
          WGLFE(K,NB,ielmn,NZ)=WGLFE(K,NB,ielmn,NZ)-XFRN
          WTLFBE(NB,ielmn,NZ)=WTLFBE(NB,ielmn,NZ)-XFRN
          EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+XFRN
          WGLFE(K,NB,ielmp,NZ)=WGLFE(K,NB,ielmp,NZ)-XFRP
          WTLFBE(NB,ielmp,NZ)=WTLFBE(NB,ielmp,NZ)-XFRP
          EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+XFRP
          WSLF(K,NB,NZ)=AZMAX1(WSLF(K,NB,NZ) &
            -AMAX1(XFRN*CNWS(NZ),XFRP*CPWS(NZ)))
        ENDIF
      ENDIF
    ENDDO D495
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
    IF(SSINN.GT.0.0)THEN
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
    TCC     =>  plt_ew%TCC          , &
    PSILT   =>  plt_ew%PSILT        , &
    WVSTKB  =>  plt_biom%WVSTKB     , &
    WTRSVBE =>  plt_biom%WTRSVBE    , &
    ZERO    =>  plt_site%ZERO       , &
    NU      =>  plt_site%NU         , &
    AREA3   =>  plt_site%AREA3      , &
    FVRN    =>  plt_allom%FVRN      , &
    VRNX    =>  plt_pheno%VRNX      , &
    IDAY    =>  plt_pheno%IDAY      , &
    IGTYP   =>  plt_pheno%IGTYP     , &
    TGSTGF  =>  plt_pheno%TGSTGF    , &
    IDTYP   =>  plt_pheno%IDTYP     , &
    FLGZ    =>  plt_pheno%FLGZ      , &
    IBTYP   =>  plt_pheno%IBTYP     , &
    VRNF    =>  plt_pheno%VRNF      , &
    TGSTGI  =>  plt_pheno%TGSTGI    , &
    ISTYP   =>  plt_pheno%ISTYP     , &
    IWTYP   =>  plt_pheno%IWTYP     , &
    CTC     =>  plt_pheno%CTC       , &
    RNH3B   =>  plt_rbgc%RNH3B      , &
    CNET    =>  plt_bgcr%CNET       , &
    SNL1    =>  plt_morph%SNL1      , &
    ARLFP   =>  plt_morph%ARLFP     , &
    NB1     =>  plt_morph%NB1         &
  )

  PSILY=real((/-200.0_r8,-2.0,-2.0,-2.0/),r8)

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
!     IDAY(2,=floral initiation date
!
  IF(IDAY(2,NB,NZ).EQ.0)THEN
    PART(1)=0.725
    PART(2)=0.275
!
!     IF BEFORE ANTHESIS
!
!     IDAY(6,=start of anthesis and setting final seed number
!     TGSTGI=total change in vegv node number normalized for maturity group
!
  ELSEIF(IDAY(6,NB,NZ).EQ.0)THEN
    PART(1)=AMAX1(PART1X,0.725-FPART1*TGSTGI(NB,NZ))
    PART(2)=AMAX1(PART2X,0.275-FPART2*TGSTGI(NB,NZ))
    PARTS=1.0_r8-PART(1)-PART(2)
    PART(3)=0.60_r8*PARTS
    PART(4)=0.30_r8*PARTS
    PARTX=PARTS-PART(3)-PART(4)
    PART(5)=0.5_r8*PARTX
    PART(6)=0.5_r8*PARTX
!
!     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     IDAY(7,=start of grain filling and setting max seed size
!     IDTYP=growth habit:0=determinate,1=indetermimate from PFT file
!     TGSTGF=total change in reprv node number normalized for maturity group
!
  ELSEIF(IDAY(7,NB,NZ).EQ.0)THEN
    IF(IDTYP(NZ).EQ.0)THEN
      PART(1)=0._r8
      PART(2)=0._r8
    ELSE
      PART(1)=AMAX1(PART1X,(0.725-FPART1)*(1.0_r8-TGSTGF(NB,NZ)))
      PART(2)=AMAX1(PART2X,(0.275-FPART2)*(1.0_r8-TGSTGF(NB,NZ)))
    ENDIF
    PARTS=1.0_r8-PART(1)-PART(2)
    PART(3)=AZMAX1(0.60*PARTS*(1.0_r8-TGSTGF(NB,NZ)))
    PART(4)=AZMAX1(0.30*PARTS*(1.0_r8-TGSTGF(NB,NZ)))
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
    IF(IDTYP(NZ).EQ.0)THEN
      PART(7)=1.0_r8
    ELSE
      PART(1)=PART1X
      PART(2)=PART2X
      PARTS=1.0_r8-PART(1)-PART(2)
      IF(ISTYP(NZ).EQ.0)THEN
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
!     IDAY(10,=physiological maturity date
!
  IF(IBTYP(NZ).EQ.0.AND.IDAY(10,NB,NZ).NE.0)THEN
    IF(ISTYP(NZ).EQ.0)THEN
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
  IF(IDAY(2,NB,NZ).NE.0)THEN
    IF(WTRSVBE(NB,ielmc,NZ).LT.XFRX*WVSTKB(NB,NZ))THEN
      DO 1020 N=1,7
        IF(N.NE.4)THEN
          PART(4)=PART(4)+0.10*PART(N)
          PART(N)=PART(N)-0.10*PART(N)
        ENDIF
1020  CONTINUE
!
!     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
!
    ELSEIF(WTRSVBE(NB,ielmc,NZ).GT.1.0*WVSTKB(NB,NZ))THEN
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
  ARLFI=ARLFP(NZ)/AREA3(NU)
  IF(ARLFI.GT.5.0)THEN
    FPARTL=AZMAX1((10.0-ARLFI)/5.0)
    PART(3)=PART(3)+(1.0_r8-FPARTL)*(PART(1)+PART(2))
    PART(1)=FPARTL*PART(1)
    PART(2)=FPARTL*PART(2)
  ENDIF
  IF(NB.EQ.NB1(NZ))THEN
    PTRT=PART(1)+PART(2)
  ENDIF
!
!     DECIDUOUS LEAF FALL AFTER GRAIN FILL IN DETERMINATES,
!     AFTER AUTUMNIZATION IN INDETERMINATES, OR AFTER SUSTAINED
!     WATER STRESS
!
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     FVRN=fraction of hours required for leafoff to initiate remobilization
!     IDAY(8,=end date for setting final seed number
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IFLGY,IFLGZ=remobilization flags
!     FLGZ=control rate of remobilization
!
  IF((ISTYP(NZ).NE.0.AND.VRNF(NB,NZ) &
    .GE.FVRN(IWTYP(NZ))*VRNX(NB,NZ)) &
    .OR.(ISTYP(NZ).EQ.0 &
    .AND.IDAY(8,NB,NZ).NE.0))THEN
    IFLGZ=1
    IF(ISTYP(NZ).EQ.0.OR.IWTYP(NZ).EQ.0)THEN
      IFLGY=1
      FLGZ(NB,NZ)=FLGZ(NB,NZ)+1.0
    ELSEIF((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.TCC(NZ).LT.CTC(NZ))THEN
      IFLGY=1
      FLGZ(NB,NZ)=FLGZ(NB,NZ)+1.0
    ELSEIF(IWTYP(NZ).GE.2 &
      .AND.PSILT(NZ).LT.PSILY(IGTYP(NZ)))THEN
      IFLGY=1
      FLGZ(NB,NZ)=FLGZ(NB,NZ)+1.0
    ENDIF
    IF(ISTYP(NZ).NE.0.AND.IWTYP(NZ).NE.0)THEN
      PART(3)=PART(3)+0.5*(PART(1)+PART(2))
      PART(4)=PART(4)+0.5*(PART(1)+PART(2))
      PART(1)=0._r8
      PART(2)=0._r8
    ENDIF
  ELSE
    IFLGZ=0
    IFLGY=0
    FLGZ(NB,NZ)=0._r8
  ENDIF
!
!     CHECK PARTITIONING COEFFICIENTS
!
  DO 1000 N=1,7
    IF(N.EQ.3.AND.test_aeqb(SNL1(NZ),0._r8))THEN
      PART(N)=0._r8
    ELSE
      PART(N)=AZMAX1(PART(N))
    ENDIF
    TOTAL=TOTAL+PART(N)
1000  CONTINUE
  IF(TOTAL.GT.ZERO)THEN
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
  associate(                        &
    RNH3B    =>  plt_rbgc%RNH3B   , &
    EPOOL    =>  plt_biom%EPOOL   , &
    IDAY     =>  plt_pheno%IDAY     &
  )
! begin_execution
! FDBK=N,P feedback inhibition on C3 CO2 fixation
! SSIN=sine of solar angle
! RADP=total PAR absorbed by canopy
! CO2Q=canopy air CO2 concentration
!
  IF(IDAY(1,NB,NZ).NE.0)THEN

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
  EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+CH2O-AMIN1(RMNCS,RCO2C)-CGROS-CNRDA
  EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)-ZADDB+RNH3B(NB,NZ)
  EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)-PADDB
  end associate
  end subroutine UpdatePhotosynthates
!------------------------------------------------------------------------------------------

  subroutine C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: CH2O3(25),CH2O4(25)
  integer :: K
  real(r8) :: CPL4M,CCBS,CPL3K,CO2LK

!     begin_execution
  associate(                              &
    CNET    =>  plt_bgcr%CNET       , &
    TCO2T   =>  plt_bgcr%TCO2T      , &
    TCO2A   =>  plt_bgcr%TCO2A      , &
    TRAU    =>  plt_bgcr%TRAU       , &
    RECO    =>  plt_bgcr%RECO       , &
    WGLFE   =>  plt_biom%WGLFE      , &
    ZEROP   =>  plt_biom%ZEROP      , &
    HCOB    =>  plt_photo%HCOB      , &
    CO2B    =>  plt_photo%CO2B      , &
    CPOOL3  =>  plt_photo%CPOOL3    , &
    CPOOL4  =>  plt_photo%CPOOL4    , &
    CO2L    =>  plt_photo%CO2L        &
  )
  DO 170 K=1,JNODS1
    IF(WGLFE(K,NB,ielmc,NZ).GT.ZEROP(NZ))THEN
!
!     MESOPHYLL TO BUNDLE SHEATH TRANSFER
!
!     WGLF=node leaf C mass
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CH2O3,CH2O4=total CO2 fixation in bundle sheath,mesophyll
!     CPL4M=mesophyll to bundle sheath transfer of nonstructural C4
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!
      CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)-CH2O3(K)
      CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)+CH2O4(K)
      CPL4M=1.0_r8*(CPOOL4(K,NB,NZ)*WGLFE(K,NB,ielmc,NZ)*FBS &
        -CPOOL3(K,NB,NZ)*WGLFE(K,NB,ielmc,NZ)*FMP) &
        /(WGLFE(K,NB,ielmc,NZ)*(FBS+FMP))
      CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)-CPL4M
      CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)+CPL4M
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
      CCBS=AZMAX1(0.083E+09*CO2B(K,NB,NZ)/(WGLFE(K,NB,ielmc,NZ)*FBS))
      CPL3K=2.5E-02*CPOOL3(K,NB,NZ)/(1.0+CCBS/CO2KI)
      CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)-CPL3K
      CO2B(K,NB,NZ)=CO2B(K,NB,NZ)+FCO2B*CPL3K
      HCOB(K,NB,NZ)=HCOB(K,NB,NZ)+FHCOB*CPL3K
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
      CO2LK=5.0E-07*(CCBS-CO2L(NZ))*WGLFE(K,NB,ielmc,NZ)*FBS
      CO2B(K,NB,NZ)=CO2B(K,NB,NZ)-FCO2B*CO2LK
      HCOB(K,NB,NZ)=HCOB(K,NB,NZ)-FHCOB*CO2LK

!
!     TOTAL C EXCHANGE
!
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!     CO2LK=bundle sheath CO2 leakage
!
      TCO2T(NZ)=TCO2T(NZ)-CO2LK
      TCO2A(NZ)=TCO2A(NZ)-CO2LK
      CNET(NZ)=CNET(NZ)-CO2LK
      RECO=RECO-CO2LK
      TRAU=TRAU-CO2LK
    ENDIF
170 CONTINUE
  end associate
  end subroutine C4PhotoProductTransfer
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
  real(r8) :: RCEL(1:npelms)
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
    WTRSVBE    =>  plt_biom%WTRSVBE    , &
    WTSTKBE    =>  plt_biom%WTSTKBE    , &
    WGNODE     =>  plt_biom%WGNODE     , &
    WGNODP     =>  plt_biom%WGNODP     , &
    WVSTKB     =>  plt_biom%WVSTKB     , &
    WGNODN     =>  plt_biom%WGNODN     , &
    WGSHE      =>  plt_biom%WGSHE      , &
    WGLFE      =>  plt_biom%WGLFE      , &
    WSLF       =>  plt_biom%WSLF       , &
    WTRVE      =>  plt_biom%WTRVE      , &
    WTLFBE     =>  plt_biom%WTLFBE     , &
    WTSHEBE    =>  plt_biom%WTSHEBE    , &
    WTSTXN     =>  plt_biom%WTSTXN     , &
    WTSTXP     =>  plt_biom%WTSTXP     , &
    WTSTXB     =>  plt_biom%WTSTXB     , &
    WSSHE      =>  plt_biom%WSSHE      , &
    ZEROP      =>  plt_biom%ZEROP      , &
    ZEROL      =>  plt_biom%ZEROL      , &
    EPOOL      =>  plt_biom%EPOOL      , &
    PP         =>   plt_site%PP        , &
    CPOOL3     =>  plt_photo%CPOOL3    , &
    CPOOL4     =>  plt_photo%CPOOL4    , &
    FWOOD      =>  plt_allom%FWOOD     , &
    FWOODN     =>  plt_allom%FWOODN    , &
    FWODSN     =>  plt_allom%FWODSN    , &
    FWODSP     =>  plt_allom%FWODSP    , &
    FWODB      =>  plt_allom%FWODB     , &
    FWODLN     =>  plt_allom%FWODLN    , &
    FWODLP     =>  plt_allom%FWODLP    , &
    FWODRN     =>  plt_allom%FWODRN    , &
    FWOODP     =>  plt_allom%FWOODP    , &
    CNWS       =>  plt_allom%CNWS      , &
    CPWS       =>  plt_allom%CPWS      , &
    ESNC       =>  plt_bgcr%ESNC       , &
    CFOPC      =>  plt_soilchem%CFOPC  , &
    CFOPN      =>  plt_soilchem%CFOPN  , &
    CFOPP      =>  plt_soilchem%CFOPP  , &
    KVSTG      =>   plt_pheno%KVSTG    , &
    IDTHB      =>   plt_pheno%IDTHB    , &
    ISTYP      =>   plt_pheno%ISTYP    , &
    ARLFB      =>   plt_morph%ARLFB    , &
    NNOD       =>   plt_morph%NNOD     , &
    HTNODE     =>   plt_morph%HTNODE   , &
    HTSHE      =>   plt_morph%HTSHE    , &
    ARLF1      =>   plt_morph%ARLF1    , &
    HTNODX     =>   plt_morph%HTNODX     &
  )
!     REMOBILIZATION AND LITTERFALL WHEN GROWTH RESPIRATION < 0
!     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
!
!     SNCX,SNCT=branch,node senescence respiration
!     KSNC=number of nodes undergoing remobilization
!

  DO 575 N=1,KSNC
    SNCT=SNCX/XKSNC
    D650: DO KK=KN,KVSTG(NB,NZ)
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
      IF(WGLFE(K,NB,ielmc,NZ).GT.ZEROP(NZ))THEN
        FNCLF=WGLFE(K,NB,ielmc,NZ)/(WGLFE(K,NB,ielmc,NZ)+WGSHE(K,NB,ielmc,NZ))
        SNCLF=FNCLF*SNCT
        SNCSH=SNCT-SNCLF
        RCCL=RCCC*WGLFE(K,NB,ielmc,NZ)
        RCZL=WGLFE(K,NB,ielmn,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCPL=WGLFE(K,NB,ielmp,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
    !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !         FSNCL,FSNAL=fraction of current leaf C,area to be remobilized
    !
        IF(RCCL.GT.ZEROP(NZ))THEN
          FSNCL=AZMAX1(AMIN1(1.0,SNCLF/RCCL))
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
        D6310: DO M=1,jsken
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(5,M,NZ) &
            *FSNCL*(WGLFE(K,NB,ielmc,NZ)-RCCL)*FWODB(0)
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(5,M,NZ) &
            *FSNCL*(WGLFE(K,NB,ielmn,NZ)-RCZL)*FWODLN(0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(5,M,NZ) &
            *FSNCL*(WGLFE(K,NB,ielmp,NZ)-RCPL)*FWODLP(0)
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(1,M,NZ) &
            *FSNCL*(WGLFE(K,NB,ielmc,NZ)-RCCL)*FWODB(1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(1,M,NZ) &
            *FSNCL*(WGLFE(K,NB,ielmn,NZ)-RCZL)*FWODLN(1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(1,M,NZ) &
            *FSNCL*(WGLFE(K,NB,ielmp,NZ)-RCPL)*FWODLP(1)
        ENDDO D6310
        IF(K.NE.0)THEN
          ESNC(2,ielmc,1,0,NZ)=ESNC(2,ielmc,1,0,NZ) &
            +FSNCL*(CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ))
          CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ) &
            -FSNCL*CPOOL3(K,NB,NZ)
          CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ) &
            -FSNCL*CPOOL4(K,NB,NZ)
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
        ARLFB(NB,NZ)=AZMAX1(ARLFB(NB,NZ)-FSNAL*ARLF1(K,NB,NZ))
        WTLFBE(NB,ielmc,NZ)=AZMAX1(WTLFBE(NB,ielmc,NZ)-FSNCL*WGLFE(K,NB,ielmc,NZ))
        WTLFBE(NB,ielmn,NZ)=AZMAX1(WTLFBE(NB,ielmn,NZ)-FSNCL*WGLFE(K,NB,ielmn,NZ))
        WTLFBE(NB,ielmp,NZ)=AZMAX1(WTLFBE(NB,ielmp,NZ)-FSNCL*WGLFE(K,NB,ielmp,NZ))
        ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)-FSNAL*ARLF1(K,NB,NZ)
        WGLFE(K,NB,ielmc,NZ)=WGLFE(K,NB,ielmc,NZ)-FSNCL*WGLFE(K,NB,ielmc,NZ)
        WGLFE(K,NB,ielmn,NZ)=WGLFE(K,NB,ielmn,NZ)-FSNCL*WGLFE(K,NB,ielmn,NZ)
        WGLFE(K,NB,ielmp,NZ)=WGLFE(K,NB,ielmp,NZ)-FSNCL*WGLFE(K,NB,ielmp,NZ)
        WSLF(K,NB,NZ)=AZMAX1(WSLF(K,NB,NZ) &
          -FSNCL*AMAX1(WGLFE(K,NB,ielmn,NZ)*CNWS(NZ) &
          ,WGLFE(K,NB,ielmp,NZ)*CPWS(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
!     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     FSNCL=fraction of current leaf C to be remobilized
!     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
!     SNCLF,SNCT=remaining senescence respiration carried to next node
!
        EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+FSNCL*RCCL*SNCF
        EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+FSNCL*RCZL
        EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+FSNCL*RCPL
        SNCLF=SNCLF-FSNCL*RCCL
        SNCT=SNCT-FSNCL*RCCL
        IF(WTLFBE(NB,ielmc,NZ).LE.ZEROL(NZ))THEN
          WTLFBE(NB,ielmc,NZ)=0._r8
          ARLFB(NB,NZ)=0._r8
        ENDIF
      !
        !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
        !
!        IF(SNCLF.LE.ZEROP(NZ))GO TO 564
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
        D6315: DO M=1,jsken
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(5,M,NZ) &
            *WGLFE(K,NB,ielmc,NZ)*FWODB(0)
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(5,M,NZ) &
            *WGLFE(K,NB,ielmn,NZ)*FWODLN(0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(5,M,NZ) &
            *WGLFE(K,NB,ielmp,NZ)*FWODLP(0)
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(1,M,NZ) &
            *WGLFE(K,NB,ielmc,NZ)*FWODB(1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(1,M,NZ) &
            *WGLFE(K,NB,ielmn,NZ)*FWODLN(1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(1,M,NZ) &
            *WGLFE(K,NB,ielmp,NZ)*FWODLP(1)
        ENDDO D6315
        IF(K.NE.0)THEN
          ESNC(2,ielmc,1,0,NZ)=ESNC(2,ielmc,1,0,NZ)+CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ)
          CPOOL3(K,NB,NZ)=0._r8
          CPOOL4(K,NB,NZ)=0._r8
        ENDIF
        ARLFB(NB,NZ)=AZMAX1(ARLFB(NB,NZ)-ARLF1(K,NB,NZ))
        WTLFBE(NB,ielmc,NZ)=AZMAX1(WTLFBE(NB,ielmc,NZ)-WGLFE(K,NB,ielmc,NZ))
        WTLFBE(NB,ielmn,NZ)=AZMAX1(WTLFBE(NB,ielmn,NZ)-WGLFE(K,NB,ielmn,NZ))
        WTLFBE(NB,ielmp,NZ)=AZMAX1(WTLFBE(NB,ielmp,NZ)-WGLFE(K,NB,ielmp,NZ))
        ARLF1(K,NB,NZ)=0._r8
        WGLFE(K,NB,ielmc,NZ)=0._r8
        WGLFE(K,NB,ielmn,NZ)=0._r8
        WGLFE(K,NB,ielmp,NZ)=0._r8
        WSLF(K,NB,NZ)=0._r8
        IF(WTLFBE(NB,ielmc,NZ).LE.ZEROL(NZ))THEN
          WTLFBE(NB,ielmc,NZ)=0._r8
          ARLFB(NB,NZ)=0._r8
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
      IF(WGSHE(K,NB,ielmc,NZ).GT.ZEROP(NZ))THEN
        RCCS=RCCC*WGSHE(K,NB,ielmc,NZ)
        RCZS=WGSHE(K,NB,ielmn,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCPS=WGSHE(K,NB,ielmp,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
      !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
      !     OR PETIOLE
      !
      !     FSNCS,FSNAS=fraction of current petiole C,length to be remobilized
      !
        IF(RCCS.GT.ZEROP(NZ))THEN
          FSNCS=AZMAX1(AMIN1(1.0,SNCSH/RCCS))
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
        D6320: DO M=1,jsken
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(5,M,NZ) &
            *FSNCS*(WGSHE(K,NB,ielmc,NZ)-RCCS)*FWODB(0)
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(5,M,NZ) &
            *FSNCS*(WGSHE(K,NB,ielmn,NZ)-RCZS)*FWODSN(0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(5,M,NZ) &
            *FSNCS*(WGSHE(K,NB,ielmp,NZ)-RCPS)*FWODSP(0)
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(2,M,NZ) &
            *FSNCS*(WGSHE(K,NB,ielmc,NZ)-RCCS)*FWODB(1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(2,M,NZ) &
            *FSNCS*(WGSHE(K,NB,ielmn,NZ)-RCZS)*FWODSN(1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(2,M,NZ) &
            *FSNCS*(WGSHE(K,NB,ielmp,NZ)-RCPS)*FWODSP(1)
        ENDDO D6320
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     HTSHE=petiole length
      !     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !     FSNCS=fraction of current petiole to be remobilized
      !     CNWS,CPWS=protein:N,protein:P ratios from startq.f
      !
        WTSHEBE(NB,ielmc,NZ)=AZMAX1(WTSHEBE(NB,ielmc,NZ)-FSNCS*WGSHE(K,NB,ielmc,NZ))
        WTSHEBE(NB,ielmn,NZ)=AZMAX1(WTSHEBE(NB,ielmn,NZ)-FSNCS*WGSHE(K,NB,ielmn,NZ))
        WTSHEBE(NB,ielmp,NZ)=AZMAX1(WTSHEBE(NB,ielmp,NZ)-FSNCS*WGSHE(K,NB,ielmp,NZ))
        HTSHE(K,NB,NZ)=HTSHE(K,NB,NZ)-FSNAS*HTSHE(K,NB,NZ)
        WGSHE(K,NB,ielmc,NZ)=WGSHE(K,NB,ielmc,NZ)-FSNCS*WGSHE(K,NB,ielmc,NZ)
        WGSHE(K,NB,ielmn,NZ)=WGSHE(K,NB,ielmn,NZ)-FSNCS*WGSHE(K,NB,ielmn,NZ)
        WGSHE(K,NB,ielmp,NZ)=WGSHE(K,NB,ielmp,NZ)-FSNCS*WGSHE(K,NB,ielmp,NZ)
        WSSHE(K,NB,NZ)=AZMAX1(WSSHE(K,NB,NZ) &
          -FSNCS*AMAX1(WGSHE(K,NB,ielmn,NZ)*CNWS(NZ) &
          ,WGSHE(K,NB,ielmp,NZ)*CPWS(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
  !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
  !
  !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  !     FSNCS=fraction of current petiole C to be remobilized
  !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
  !     SNCSH,SNCT=remaining senescence respiration carried to next node
  !
        EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+FSNCS*RCCS*SNCF
        EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+FSNCS*RCZS
        EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+FSNCS*RCPS
        SNCSH=SNCSH-FSNCS*RCCS
        SNCT=SNCT-FSNCS*RCCS
        IF(WTSHEBE(NB,ielmc,NZ).LE.ZEROL(NZ))THEN
          WTSHEBE(NB,ielmc,NZ)=0._r8
        ENDIF
!
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
        IF(SNCSH.LE.ZEROP(NZ))GO TO 565
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
        D6325: DO M=1,jsken
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(5,M,NZ) &
            *WGSHE(K,NB,ielmc,NZ)*FWODB(0)
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(5,M,NZ) &
            *WGSHE(K,NB,ielmn,NZ)*FWODSN(0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(5,M,NZ) &
            *WGSHE(K,NB,ielmp,NZ)*FWODSP(0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPC(2,M,NZ) &
            *WGSHE(K,NB,ielmc,NZ)*FWODB(1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(2,M,NZ) &
            *WGSHE(K,NB,ielmn,NZ)*FWODSN(1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(2,M,NZ) &
            *WGSHE(K,NB,ielmp,NZ)*FWODSP(1)
        ENDDO D6325
        WTSHEBE(NB,ielmc,NZ)=AZMAX1(WTSHEBE(NB,ielmc,NZ)-WGSHE(K,NB,ielmc,NZ))
        WTSHEBE(NB,ielmn,NZ)=AZMAX1(WTSHEBE(NB,ielmn,NZ)-WGSHE(K,NB,ielmn,NZ))
        WTSHEBE(NB,ielmp,NZ)=AZMAX1(WTSHEBE(NB,ielmp,NZ)-WGSHE(K,NB,ielmp,NZ))
        HTSHE(K,NB,NZ)=0._r8
        WGSHE(K,NB,ielmc,NZ)=0._r8
        WGSHE(K,NB,ielmn,NZ)=0._r8
        WGSHE(K,NB,ielmp,NZ)=0._r8
        WSSHE(K,NB,NZ)=0._r8
        IF(WTSHEBE(NB,ielmc,NZ).LE.ZEROL(NZ))THEN
          WTSHEBE(NB,ielmc,NZ)=0._r8
        ENDIF
      ENDIF
    ENDDO D650
    KN=KN+1
    SNCR=SNCT*(1.0_r8-SNCF)
!
!     REMOBILIZATION OF RESERVE C
!
!     WTRSVB=stalk reserve C mass
!     SNCR=excess maintenance respiration
!
    IF(WTRSVBE(NB,ielmc,NZ).GT.SNCR)THEN
      WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)-SNCR
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
    IF(ISTYP(NZ).NE.0.AND.SNCT.GT.ZEROP(NZ) &
      .AND.WTSTKBE(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
      SNCF=SNCZ/SNCT
      FRCC=WVSTKB(NB,NZ)/WTSTKBE(NB,ielmc,NZ)
      RCSC=RCCC*FRCC
      RCSN=RCCN*FRCC
      RCSP=RCCP*FRCC
      MXNOD=KVSTG(NB,NZ)
      MNNOD=MAX(MIN(0,MAX(0,MXNOD-NNOD(NZ))),KVSTG(NB,NZ)-23)
      MXNOD=MAX(MXNOD,MNNOD)
      D1650: DO KK=MXNOD,MNNOD,-1
        K=MOD(KK,JNODS1)
        IF(K.EQ.0.AND.KK.NE.0)K=25
    !
    !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !
        IF(WGNODE(K,NB,NZ).GT.ZEROP(NZ))THEN
          RCCK=RCSC*WGNODE(K,NB,NZ)
          RCZK=WGNODN(K,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
          RCPK=WGNODP(K,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
      !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
      !
      !     FSNCS=fraction of lowest internode to be remobilized
!
          IF(RCCK.GT.ZEROP(NZ))THEN
            FSNCK=AZMAX1(AMIN1(1.0,SNCT/RCCK))
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
          D7310: DO M=1,jsken
            ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
              *FSNCK*(WGNODE(K,NB,NZ)-RCCK)*FWOOD(0)
            ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
              *FSNCK*(WGNODN(K,NB,NZ)-RCZK)*FWOODN(0)
            ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
              *FSNCK*(WGNODP(K,NB,NZ)-RCPK)*FWOODP(0)
            ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
              *FSNCK*(WGNODE(K,NB,NZ)-RCCK)*FWOOD(1)
            ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
              *FSNCK*(WGNODN(K,NB,NZ)-RCZK)*FWOODN(1)
            ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
              *FSNCK*(WGNODP(K,NB,NZ)-RCPK)*FWOODP(1)
          ENDDO D7310
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     FSNCK=fraction of lowest internode to be remobilized
      !     HTNODE,HTNODX=living,senescing internode length
      !     WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
      !
          WTSTKBE(NB,ielmc,NZ)=AZMAX1(WTSTKBE(NB,ielmc,NZ)-FSNCK*WGNODE(K,NB,NZ))
          WTSTKBE(NB,ielmn,NZ)=AZMAX1(WTSTKBE(NB,ielmn,NZ)-FSNCK*WGNODN(K,NB,NZ))
          WTSTKBE(NB,ielmp,NZ)=AZMAX1(WTSTKBE(NB,ielmp,NZ)-FSNCK*WGNODP(K,NB,NZ))
          HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)-FSNCK*HTNODX(K,NB,NZ)
          WGNODE(K,NB,NZ)=WGNODE(K,NB,NZ)-FSNCK*WGNODE(K,NB,NZ)
          WGNODN(K,NB,NZ)=WGNODN(K,NB,NZ)-FSNCK*WGNODN(K,NB,NZ)
          WGNODP(K,NB,NZ)=WGNODP(K,NB,NZ)-FSNCK*WGNODP(K,NB,NZ)
          HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ)-FSNCK*HTNODX(K,NB,NZ)
!
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
      !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
      !     FSNCK=fraction of lowest internode to be remobilized
      !     SNCT=remaining node senescence respiration
      !
          WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+FSNCK*RCCK*SNCF
          WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)+FSNCK*RCZK
          WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)+FSNCK*RCPK
          SNCT=SNCT-FSNCK*RCCK
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
          IF(SNCT.LE.ZEROP(NZ))GO TO 565
    !
    !       OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
    !
    !       CSNC,ZSNC,PSNC=literfall C,N,P
    !       CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !       WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
    !       HTNODE,HTNODX=living,senescing internode length
    !
        ELSE
          D7315: DO M=1,jsken
            ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
              *WGNODE(K,NB,NZ)*FWOOD(0)
            ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
              *WGNODN(K,NB,NZ)*FWOODN(0)
            ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
              *WGNODP(K,NB,NZ)*FWOODP(0)
            ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
              *WGNODE(K,NB,NZ)*FWOOD(1)
            ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
              *WGNODN(K,NB,NZ)*FWOODN(1)
            ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
              *WGNODP(K,NB,NZ)*FWOODP(1)
          ENDDO D7315
          WTSTKBE(NB,ielmc,NZ)=AZMAX1(WTSTKBE(NB,ielmc,NZ)-WGNODE(K,NB,NZ))
          WTSTKBE(NB,ielmn,NZ)=AZMAX1(WTSTKBE(NB,ielmn,NZ)-WGNODN(K,NB,NZ))
          WTSTKBE(NB,ielmp,NZ)=AZMAX1(WTSTKBE(NB,ielmp,NZ)-WGNODP(K,NB,NZ))
          HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)-HTNODX(K,NB,NZ)
          WGNODE(K,NB,NZ)=0._r8
          WGNODN(K,NB,NZ)=0._r8
          WGNODP(K,NB,NZ)=0._r8
          HTNODX(K,NB,NZ)=0._r8
        ENDIF
      ENDDO D1650
!
    !   RESIDUAL STALK
    !
    !   RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
      IF(WTSTXB(NB,NZ).GT.ZEROP(NZ))THEN
        RCCK=RCSC*WTSTXB(NB,NZ)
        RCZK=WTSTXN(NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
        RCPK=WTSTXP(NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
    !     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
        IF(RCCK.GT.ZEROP(NZ))THEN
          FSNCR=AZMAX1(AMIN1(1.0,SNCT/RCCK))
        ELSE
          FSNCR=1.0_r8
        ENDIF
    !
    !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !     TO FRACTIONS SET IN 'STARTQ'
    !
        D8310: DO M=1,jsken
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
            *FSNCR*(WTSTXB(NB,NZ)-RCCK)*FWOOD(0)
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
            *FSNCR*(WTSTXN(NB,NZ)-RCZK)*FWOODN(0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
            *FSNCR*(WTSTXP(NB,NZ)-RCPK)*FWOODP(0)
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
            *FSNCR*(WTSTXB(NB,NZ)-RCCK)*FWOOD(1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
            *FSNCR*(WTSTXN(NB,NZ)-RCZK)*FWOODN(1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
            *FSNCR*(WTSTXP(NB,NZ)-RCPK)*FWOODP(1)
        ENDDO D8310
    !
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !     HTNODE,HTNODX=living,senescing internode length
    !
        WTSTKBE(NB,ielmc,NZ)=AZMAX1(WTSTKBE(NB,ielmc,NZ) &
          -FSNCR*WTSTXB(NB,NZ))
        WTSTKBE(NB,ielmn,NZ)=AZMAX1(WTSTKBE(NB,ielmn,NZ) &
          -FSNCR*WTSTXN(NB,NZ))
        WTSTKBE(NB,ielmp,NZ)=AZMAX1(WTSTKBE(NB,ielmp,NZ) &
          -FSNCR*WTSTXP(NB,NZ))
        WTSTXB(NB,NZ)=AZMAX1(WTSTXB(NB,NZ) &
          -FSNCR*WTSTXB(NB,NZ))
        WTSTXN(NB,NZ)=AZMAX1(WTSTXN(NB,NZ) &
          -FSNCR*WTSTXN(NB,NZ))
        WTSTXP(NB,NZ)=AZMAX1(WTSTXP(NB,NZ) &
          -FSNCR*WTSTXP(NB,NZ))
        HTNODZ=0._r8
        DO 8320 K=0,JNODS1
          HTNODZ=AMAX1(HTNODZ,HTNODE(K,NB,NZ))
8320    CONTINUE
        HTNODZ=AZMAX1(HTNODZ-FSNCR*HTNODZ)
        DO 8325 K=0,JNODS1
          HTNODE(K,NB,NZ)=AMIN1(HTNODZ,HTNODE(K,NB,NZ))
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
        WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+FSNCR*RCCK*SNCF
        WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)+FSNCR*RCZK
        WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)+FSNCR*RCPK
        SNCT=SNCT-FSNCR*RCCK
      ENDIF
  !
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
      IF(SNCT.LE.ZEROP(NZ))cycle
  !
  !     OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
  !
  !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
  !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
  !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
  !
    ELSE
      D8315: DO M=1,jsken
        ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
          *WTSTXB(NB,NZ)*FWOOD(0)
        ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
          *WTSTXN(NB,NZ)*FWOODN(0)
        ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
          *WTSTXP(NB,NZ)*FWOODP(0)
        ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,0,0,NZ)+CFOPC(3,M,NZ) &
          *WTSTXB(NB,NZ)*FWOOD(1)
        ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPN(3,M,NZ) &
          *WTSTXN(NB,NZ)*FWOODN(1)
        ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPP(3,M,NZ) &
          *WTSTXP(NB,NZ)*FWOODP(1)
      ENDDO D8315
      WTSTKBE(NB,ielmc,NZ)=AZMAX1(WTSTKBE(NB,ielmc,NZ)-WTSTXB(NB,NZ))
      WTSTKBE(NB,ielmn,NZ)=AZMAX1(WTSTKBE(NB,ielmn,NZ)-WTSTXN(NB,NZ))
      WTSTKBE(NB,ielmp,NZ)=AZMAX1(WTSTKBE(NB,ielmp,NZ)-WTSTXP(NB,NZ))
      WTSTXB(NB,NZ)=0._r8
      WTSTXN(NB,NZ)=0._r8
      WTSTXP(NB,NZ)=0._r8
    ENDIF
!
!     REMOBILIZATION OF STORAGE C,N,P
!
!     WTRVC=storage C
!     IDTHB=branch living flag: 0=alive,1=dead
!     SNCR=remaining excess maintenance respiration
!
    SNCR=SNCT/(1.0+FXFS)
    IF(WTRVE(ielmc,NZ).GT.SNCR)THEN
      WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)-SNCR
      SNCR=0._r8
    ELSEIF(ISTYP(NZ).NE.0)THEN
      IDTHB(NB,NZ)=1
    ENDIF
565 CONTINUE
575 CONTINUE
  end associate
  end subroutine RemobilizeLeafLayers


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
    WGLFL    =>  plt_biom%WGLFL   , &
    WGLFE    =>  plt_biom%WGLFE   , &
    WGLFV    =>  plt_biom%WGLFV   , &
    WGLFLP   =>  plt_biom%WGLFLP  , &
    WVSTKB   =>  plt_biom%WVSTKB  , &
    WTSTKBE  =>  plt_biom%WTSTKBE , &
    WGLFLN   =>  plt_biom%WGLFLN  , &
    WSSHE    =>  plt_biom%WSSHE   , &
    FNOD     =>  plt_allom%FNOD   , &
    IGTYP    =>  plt_pheno%IGTYP  , &
    ISTYP    =>  plt_pheno%ISTYP  , &
    KVSTG    =>  plt_pheno%KVSTG  , &
    IBTYP    =>  plt_pheno%IBTYP  , &
    KVSTGN   =>  plt_pheno%KVSTGN , &
    WDLF     => plt_morph%WDLF    , &
    SDPTH    => plt_morph%SDPTH   , &
    NB1      => plt_morph%NB1     , &
    ARLFL    => plt_morph%ARLFL   , &
    HTSHE    => plt_morph%HTSHE   , &
    ARSTK    => plt_morph%ARSTK   , &
    NBTB     => plt_morph%NBTB    , &
    HTNODE   => plt_morph%HTNODE  , &
    ZL       => plt_morph%ZL      , &
    ZC       => plt_morph%ZC      , &
    HTCTL    => plt_morph%HTCTL   , &
    CLASS    => plt_morph%CLASS   , &
    ARLF1    => plt_morph%ARLF1   , &
    ARLFV    => plt_morph%ARLFV   , &
    PP       => plt_site%PP       , &
    ZERO     => plt_site%ZERO     , &
    ZSIN     => plt_rad%ZSIN        &
  )
!   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
!
!   HTCTL=hypocotyledon height
!   SDPTH=seeding depth
!   ARLF=node leaf area
!   HTSHE=petiole length
!
  KVSTGN(NB,NZ)=0
  IF(HTCTL(NZ).LE.SDPTH(NZ) &
    .AND.ARLF1(0,NB1(NZ),NZ).GT.0.0)THEN
    XLGLF=SQRT(1.0E+02*ARLF1(0,NB1(NZ),NZ)/PP(NZ))
    HTCTL(NZ)=XLGLF+HTSHE(0,NB1(NZ),NZ)+HTNODE(0,NB1(NZ),NZ)
  ENDIF
!
! IF CANOPY HAS EMERGED
!
  IF(HTCTL(NZ).GT.SDPTH(NZ))THEN
    D540: DO K=0,JNODS1
      DO  L=1,JC1
        ARLFL(L,K,NB,NZ)=0._r8
        WGLFL(L,K,NB,NZ)=0._r8
        WGLFLN(L,K,NB,NZ)=0._r8
        WGLFLP(L,K,NB,NZ)=0._r8
      enddo
    ENDDO D540
    D535: DO L=1,JC1
      ARSTK(L,NB,NZ)=0._r8
    ENDDO D535
!
!   BRANCH HEIGHT
!
!   IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   KVSTG1,KVSTGN=integer of highest,lowest leaf number currently growing
!   HTNODE=internode length
!   HTBR=branch base height
!
    IF(IBTYP(NZ).NE.0.AND.IGTYP(NZ).GT.1)THEN
      IF(NB.NE.NB1(NZ))THEN
        KVSTG1=MAX(1,KVSTG(NB1(NZ),NZ)-24)
        IF(NBTB(NB,NZ).GE.KVSTG1)THEN
          K=MOD(NBTB(NB,NZ),JNODS1)
!        IF(K.EQ.0.AND.KK.NE.0)K=JNODS1
          IF(K.EQ.0)K=JNODS1
          HTBR=HTNODE(K,NB1(NZ),NZ)
        ELSE
          HTBR=0._r8
        ENDIF
      ELSE
        HTBR=0._r8
      ENDIF
    ELSE
      HTBR=0._r8
    ENDIF
    KVSTGX=MAX(0,KVSTG(NB,NZ)-24)
!
!   FOR ALL LEAFED NODES
!
    D560: DO KK=KVSTGX,KVSTG(NB,NZ)
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
      HTSTK=HTBR+HTNODE(K,NB,NZ)
      HTLF=HTSTK+HTSHE(K,NB,NZ)
      XLGLF=AZMAX1(SQRT(WDLF(NZ)*AMAX1(0.0 &
        ,ARLF1(K,NB,NZ))/(PP(NZ)*FNOD(NZ))))
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
      D555: DO N=JLI1,1,-1
        YLGLF=ZSIN(N)*CLASS(N,NZ)*XLGLF
        HTLFL=AMIN1(ZCX(NZ)+0.01-YLGLF,HTLF+TLGLF)
        HTLFU=AMIN1(ZCX(NZ)+0.01,HTLFL+YLGLF)
        LU=0
        LL=0
        DO 550 L=JC1,1,-1
          IF(LU.EQ.1.AND.LL.EQ.1)exit
          IF((HTLFU.GT.ZL(L-1).OR.ZL(L-1).LE.ZERO).AND.LU.EQ.0)THEN
            LHTLFU=MAX(1,L)
            LU=1
          ENDIF
          IF((HTLFL.GT.ZL(L-1).OR.ZL(L-1).LE.ZERO).AND.LL.EQ.0)THEN
            LHTLFL=MAX(1,L)
            LL=1
          ENDIF
550     CONTINUE
        D570: DO L=LHTLFL,LHTLFU
          IF(LHTLFU.EQ.LHTLFL)THEN
            FRACL=CLASS(N,NZ)
          ELSEIF(HTLFU.GT.HTLFL.AND.ZL(L).GT.HTLFL)THEN
            FRACL=CLASS(N,NZ)*(AMIN1(HTLFU,ZL(L)) &
              -AMAX1(HTLFL,ZL(L-1)))/(HTLFU-HTLFL)
          ELSE
            FRACL=CLASS(N,NZ)
          ENDIF
          YARLF=FRACL*ARLF1(K,NB,NZ)
          YWGLF=FRACL*WGLFE(K,NB,ielmc,NZ)
          YWGLFN=FRACL*WGLFE(K,NB,ielmn,NZ)
          YWGLFP=FRACL*WGLFE(K,NB,ielmp,NZ)
!
    !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
    !
    !     ARLFL=leaf node area in canopy layer
    !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
    !     ARLFV,WGLFV=total leaf area,C in canopy layer
    !     HTNODE=internode length
    !
          ARLFL(L,K,NB,NZ)=ARLFL(L,K,NB,NZ)+YARLF
          WGLFL(L,K,NB,NZ)=WGLFL(L,K,NB,NZ)+YWGLF
          WGLFLN(L,K,NB,NZ)=WGLFLN(L,K,NB,NZ)+YWGLFN
          WGLFLP(L,K,NB,NZ)=WGLFLP(L,K,NB,NZ)+YWGLFP
          ARLFV(L,NZ)=ARLFV(L,NZ)+YARLF
          WGLFV(L,NZ)=WGLFV(L,NZ)+YWGLF
        ENDDO D570
        TLGLF=TLGLF+YLGLF
        ZC(NZ)=AMAX1(ZC(NZ),HTLFU)
      ENDDO D555
      IF(WSSHE(K,NB,NZ).GT.0.0)THEN
        IF(KVSTGN(NB,NZ).EQ.0)KVSTGN(NB,NZ)=MIN(KK,KVSTG(NB,NZ))
      ENDIF
    ENDDO D560
    IF(KVSTGN(NB,NZ).EQ.0)KVSTGN(NB,NZ)=KVSTG(NB,NZ)
    K1=MOD(KVSTG(NB,NZ),JNODS1)
    IF(K1.EQ.0.AND.KVSTG(NB,NZ).NE.0)K1=JNODS1
    K2=MOD(KVSTG(NB,NZ)-1,JNODS1)
    IF(K2.EQ.0.AND.KVSTG(NB,NZ)-1.NE.0)K2=JNODS1
    IF(test_aeqb(HTNODE(K1,NB,NZ),0._r8))THEN
      HTNODE(K1,NB,NZ)=HTNODE(K2,NB,NZ)
    ENDIF
    HTLFB=HTBR+AZMAX1(HTNODE(K1,NB,NZ))
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
  !     WRITE(*,6679)'K1',I,J,NZ,NB,K1,KVSTG(NB,NZ)
  !    2,HTNODE(K1,NB,NZ)
  !6679  FORMAT(A8,6I4,12E12.4)
  !     ENDIF
    IF(HTNODE(K1,NB,NZ).GT.0.0)THEN
      LU=0
      LL=0
      DO 545 L=JC1,1,-1
        IF(LU.EQ.1.AND.LL.EQ.1)exit
        IF((HTLFB.GT.ZL(L-1).OR.ZL(L-1).LE.ZERO).AND.LU.EQ.0)THEN
          LHTBRU=MAX(1,L)
          LU=1
        ENDIF
        IF((HTBR.GT.ZL(L-1).OR.ZL(L-1).LE.ZERO).AND.LL.EQ.0)THEN
          LHTBRL=MAX(1,L)
          LL=1
        ENDIF
545   CONTINUE
      RSTK=SQRT(VSTK*(AZMAX1(WTSTKBE(NB,ielmc,NZ))/PP(NZ)) &
        /(PICON*HTNODE(K1,NB,NZ)))
      ARSTKB=PICON*HTNODE(K1,NB,NZ)*PP(NZ)*RSTK
      IF(ISTYP(NZ).EQ.0)THEN
        WVSTKB(NB,NZ)=WTSTKBE(NB,ielmc,NZ)
      ELSE
        ZSTK=AMIN1(ZSTX,FSTK*RSTK)
        ASTV=PICON*(2.0_r8*RSTK*ZSTK-ZSTK**2)
        WVSTKB(NB,NZ)=ASTV/VSTK*HTNODE(K1,NB,NZ)*PP(NZ)
      ENDIF
    !     IF(NZ.EQ.1)THEN
        !     WRITE(*,6677)'WVSTK',I,J,NX,NY,NZ,NB,K1,WVSTKB(NB,NZ)
        !    2,ASTV,VSTK,HTNODE(K1,NB,NZ),PP(NZ)
        !6677  FORMAT(A8,7I4,12E12.4)
        !     ENDIF
      DO 445 L=LHTBRL,LHTBRU
        IF(HTLFB.GT.HTBR)THEN
          IF(HTLFB.GT.ZL(L-1))THEN
            FRACL=(AMIN1(HTLFB,ZL(L))-AMAX1(HTBR &
              ,ZL(L-1)))/(HTLFB-HTBR)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        ARSTK(L,NB,NZ)=FRACL*ARSTKB
445   CONTINUE
    ELSE
      WVSTKB(NB,NZ)=0._r8
      DO 450 L=1,JC1
        ARSTK(L,NB,NZ)=0._r8
450   CONTINUE
    ENDIF
  ELSE
    WVSTKB(NB,NZ)=0._r8
    DO 455 L=1,JC1
      ARSTK(L,NB,NZ)=0._r8
455 CONTINUE
  ENDIF
  end associate
  end subroutine AllocateLeafToCanopyLayers
!------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8) :: dangle
  integer :: L,K,N
  ! begin_execution
  associate(                             &
    ANGBR   =>  plt_morph%ANGBR    , &
    ARLF1   =>  plt_morph%ARLF1    , &
    ARSTK   =>  plt_morph%ARSTK    , &
    CLASS   =>  plt_morph%CLASS    , &
    ARLFL   =>  plt_morph%ARLFL    , &
    SURFB   =>  plt_morph%SURFB    , &
    NB1     =>  plt_morph%NB1      , &
    SURF    =>  plt_morph%SURF       &
  )
  DO 900 K=1,JNODS1
    DO  L=1,JC1
      DO  N=1,JLI1
        SURF(N,L,K,NB,NZ)=0._r8
      enddo
    enddo
900 CONTINUE
! ARLFXB=0._r8
! ARLFXL=0._r8
! SURFXX=0._r8
  DO 500 K=1,JNODS1
!     ARLFXB=ARLFXB+ARLF1(K,NB,NZ)
    IF(ARLF1(K,NB,NZ).GT.0.0)THEN
      DO 700 L=JC1,1,-1
!       ARLFXL=ARLFXL+ARLFL(L,K,NB,NZ)
        DO 800 N=1,JLI1
          SURF(N,L,K,NB,NZ)=AZMAX1(CLASS(N,NZ) &
            *0.25_r8*ARLFL(L,K,NB,NZ))
  !       SURFXX=SURFXX+SURF(N,L,K,NB,NZ)

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
      SURFB(N,L,NB,NZ)=0._r8
    enddo
910  CONTINUE

  IF(NB.EQ.NB1(NZ))THEN
    N=JLI1
  ELSE
    dangle=PICON2/real(JLI1,r8)
    N=MIN(JLI1,INT(ASIN(ANGBR(NZ))/dangle)+1)
  ENDIF
  DO 710 L=JC1,1,-1
    SURFB(N,L,NB,NZ)=ARSTK(L,NB,NZ)/real(JLI1,r8)
710   CONTINUE
  end associate
  end subroutine LeafClassAllocation

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
    TCC      =>  plt_ew%TCC         , &
    CEPOLB   =>  plt_biom%CEPOLB    , &
    WTGRBE   =>  plt_biom%WTGRBE    , &
    WTRSVBE  =>  plt_biom%WTRSVBE   , &
    ZEROP    =>  plt_biom%ZEROP     , &
    GRWTB    =>  plt_allom%GRWTB    , &
    CNGR     =>  plt_allom%CNGR     , &
    CPGR     =>  plt_allom%CPGR     , &
    TFN4     =>  plt_pheno%TFN4     , &
    TFN3     =>  plt_pheno%TFN3     , &
    IWTYP    =>  plt_pheno%IWTYP    , &
    VRNF     =>  plt_pheno%VRNF     , &
    FLG4     =>  plt_pheno%FLG4     , &
    VRNX     =>  plt_pheno%VRNX     , &
    GFILL    =>  plt_pheno%GFILL    , &
    SSTX     =>  plt_pheno%SSTX     , &
    ISTYP    =>  plt_pheno%ISTYP    , &
    HTC      =>  plt_pheno%HTC      , &
    DGSTGF   =>  plt_pheno%DGSTGF   , &
    IDAY     =>  plt_pheno%IDAY     , &
    CTC      =>  plt_pheno%CTC      , &
    PP       =>  plt_site%PP        , &
    SDMX     =>  plt_morph%SDMX     , &
    GRMX     =>  plt_morph%GRMX     , &
    STMX     =>  plt_morph%STMX     , &
    NG       =>  plt_morph%NG       , &
    GRNXB    =>  plt_morph%GRNXB    , &
    IRTYP    =>  plt_morph%IRTYP    , &
    GRNOB    =>  plt_morph%GRNOB      &
  )
!
!   SET MAXIMUM GRAIN NUMBER FROM SHOOT MASS BEFORE ANTHESIS
!
!   IDAY(3,=start of stem elongation and setting max seed number
!   IDAY(6,=start of anthesis and setting final seed number
!   GRNXB=potential number of seed set sites
!   STMX=maximum potential seed number from PFT file
!     GROSTK=stalk growth rate
!
  IF(IDAY(3,NB,NZ).NE.0.AND.IDAY(6,NB,NZ).EQ.0)THEN
    GRNXB(NB,NZ)=GRNXB(NB,NZ)+STMX(NZ)*AZMAX1(GROSTK)
  ENDIF
!
!   SET FINAL GRAIN NUMBER FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!   IDAY(6,=start of anthesis and setting final seed number
!   IDAY(7,=start of grain filling and setting max seed size
!   IDAY(8,=end date setting for final seed number
!   IDAY(9,=end of setting max seed size
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
  IF(IDAY(6,NB,NZ).NE.0.AND.IDAY(9,NB,NZ).EQ.0)THEN
    SET=AMIN1(CEPOLB(NB,ielmc,NZ)/(CEPOLB(NB,ielmc,NZ)+SETC) &
      ,CEPOLB(NB,ielmn,NZ)/(CEPOLB(NB,ielmn,NZ)+SETN) &
      ,CEPOLB(NB,ielmp,NZ)/(CEPOLB(NB,ielmp,NZ)+SETP))
    IF(TCC(NZ).LT.CTC(NZ))THEN
      IF(IDAY(7,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(CTC(NZ)-TCC(NZ))
      ELSEIF(IDAY(8,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(CTC(NZ)-TCC(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSEIF(TCC(NZ).GT.HTC(NZ))THEN
      IF(IDAY(7,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(TCC(NZ)-HTC(NZ))
      ELSEIF(IDAY(8,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(TCC(NZ)-HTC(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSE
      FGRNX=0._r8
    ENDIF
    IF(IDAY(6,NB,NZ).NE.0.AND.IDAY(8,NB,NZ).EQ.0)THEN
!     GRNXB(NB,NZ)=GRNXB(NB,NZ)*FGRNX
      GRNOB(NB,NZ)=AMIN1(SDMX(NZ)*GRNXB(NB,NZ) &
        ,GRNOB(NB,NZ)+SDMX(NZ)*GRNXB(NB,NZ) &
        *SET*DGSTGF(NB,NZ)-FGRNX*GRNOB(NB,NZ))
    ENDIF
!

!     SET MAXIMUM GRAIN SIZE FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!     GRMX=maximum individual seed size from PFT file (g)
!     DGSTGF=change in reproductive node number normalized for maturity group
!     GRWTB=individual seed size
!     SET=seed set limited by nonstructural C,N,P
!
    IF(IDAY(7,NB,NZ).NE.0.AND.IDAY(9,NB,NZ).EQ.0)THEN
      GRMXB=GRMX(NZ)
      GRWTB(NB,NZ)=AMIN1(GRMX(NZ),GRWTB(NB,NZ) &
        +GRMXB*AMAX1(0.50,SET**0.25)*DGSTGF(NB,NZ))

    ENDIF
  ENDIF
!
!   GRAIN FILL BY TRANSLOCATION FROM STALK RESERVES
!   UNTIL GRAIN SINK (=FINAL GRAIN NUMBER X MAXIMUM
!   GRAIN SIZE) IS FILLED OR RESERVES ARE EXHAUSTED
!
!   IDAY(7,=start of grain filling and setting max seed size
!   WTGRB=total seed C mass
!   GRWTB=individual seed size
!   GRNOB=seed set number
!   GROLM=maximum grain fill rate
!   GFILL=grain filling rate at 25 oC from PFT file
!   TFN3=temperature function for canopy growth
!   TFN4=temperature function for root growth
!
  IF(IDAY(7,NB,NZ).NE.0)THEN
    IF(WTGRBE(NB,ielmc,NZ).GE.GRWTB(NB,NZ)*GRNOB(NB,NZ))THEN
      GROLM=0._r8
    ELSEIF(IRTYP(NZ).EQ.0)THEN
      GROLM=AZMAX1(GFILL(NZ)*GRNOB(NB,NZ)*SQRT(TFN3(NZ)))
    ELSE
      GROLM=AZMAX1(GFILL(NZ)*GRNOB(NB,NZ) &
        *SQRT(TFN4(NG(NZ),NZ)))
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
    IF(WTGRBE(NB,ielmn,NZ).LT.ZPGRM*CNGR(NZ) &
      *WTGRBE(NB,ielmc,NZ).OR.WTGRBE(NB,ielmp,NZ).LT.ZPGRM &
      *CPGR(NZ)*WTGRBE(NB,ielmc,NZ))THEN
      GROLC=0._r8
    ELSE
      GROLC=GROLM
    ENDIF
    XLOCM=AMIN1(GROLM,WTRSVBE(NB,ielmc,NZ))
    XLOCC=AMIN1(GROLC,WTRSVBE(NB,ielmc,NZ))
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
    IF(WTRSVBE(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
      ZNPGN=WTRSVBE(NB,ielmn,NZ)/(WTRSVBE(NB,ielmn,NZ) &
        +SETN*WTRSVBE(NB,ielmc,NZ))
      ZNPGP=WTRSVBE(NB,ielmp,NZ)/(WTRSVBE(NB,ielmp,NZ) &
        +SETP*WTRSVBE(NB,ielmc,NZ))
      ZPGRN=ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0,ZNPGN))
      ZPGRP=ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0,ZNPGP))
      XLOCN=AMIN1(XLOCM*CNGR(NZ) &
        ,AZMAX1(WTRSVBE(NB,ielmn,NZ)*ZPGRN) &
        ,(WTGRBE(NB,ielmc,NZ)+XLOCC)*CNGR(NZ)-WTGRBE(NB,ielmn,NZ))
      XLOCP=AMIN1(XLOCM*CPGR(NZ) &
        ,AZMAX1(WTRSVBE(NB,ielmp,NZ)*ZPGRP) &
        ,(WTGRBE(NB,ielmc,NZ)+XLOCC)*CPGR(NZ)-WTGRBE(NB,ielmp,NZ))
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
    WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+GROGR-XLOCC
    WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)+GROGRN-XLOCN
    WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)+GROGRP-XLOCP
    WTGRBE(NB,ielmc,NZ)=WTGRBE(NB,ielmc,NZ)+XLOCC
    WTGRBE(NB,ielmn,NZ)=WTGRBE(NB,ielmn,NZ)+XLOCN
    WTGRBE(NB,ielmp,NZ)=WTGRBE(NB,ielmp,NZ)+XLOCP
  ELSE
    XLOCC=0._r8
    XLOCN=0._r8
    XLOCP=0._r8
  ENDIF
!
!   SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
!   HAS STOPPED FOR SET PERIOD OF TIME
!
!   IDAY(8,=end date setting for final seed number
!   XLOCC=C translocation rate from reserve to grain
!   PP=PFT population
!   FLG4=number of hours with no grain fill
!   FLG4X=number of hours with no grain filling until physl maturity
!   IDAY(10,=date of physiological maturity
!
  IF(IDAY(8,NB,NZ).NE.0)THEN
    IF(XLOCC.LE.1.0E-09*PP(NZ))THEN
      FLG4(NB,NZ)=FLG4(NB,NZ)+1.0
    ELSE
      FLG4(NB,NZ)=0._r8
    ENDIF
    IF(FLG4(NB,NZ).GE.FLG4X)THEN
      IF(IDAY(10,NB,NZ).EQ.0)THEN
        IDAY(10,NB,NZ)=I
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
    IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
      IF(FLG4(NB,NZ).GT.FLG4X+FLG4Y(IWTYP(NZ)))THEN
        VRNF(NB,NZ)=VRNX(NB,NZ)+0.5
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine GrainFilling
!------------------------------------------------------------------------------------------

  subroutine PhenologyReset(I,NB,NZ)

  implicit none
  integer, intent(in) :: I,nb,nz
  integer :: K,M
! begin_execution
  associate(                           &
    EHVST    =>  plt_distb%EHVST     , &
    IYRH     =>  plt_distb%IYRH      , &
    THIN     =>  plt_distb%THIN      , &
    HVST     =>  plt_distb%HVST      , &
    IHVST    =>  plt_distb%IHVST     , &
    JHVST    =>  plt_distb%JHVST     , &
    IDAY0    =>  plt_distb%IDAY0     , &
    IDAYH    =>  plt_distb%IDAYH     , &
    IYR0     =>  plt_distb%IYR0      , &
    WTGRBE   =>  plt_biom%WTGRBE     , &
    WGLFE    =>  plt_biom%WGLFE      , &
    WTRTD    =>  plt_biom%WTRTD      , &
    WTSTXB   =>  plt_biom%WTSTXB     , &
    WTSTXP   =>  plt_biom%WTSTXP     , &
    WTSTKBE  =>  plt_biom%WTSTKBE    , &
    EPOOLR   =>  plt_biom%EPOOLR     , &
    WTSHEBE  =>  plt_biom%WTSHEBE    , &
    WTEARBE  =>  plt_biom%WTEARBE    , &
    WSLF     =>  plt_biom%WSLF       , &
    WSSHE    =>  plt_biom%WSSHE      , &
    WTHSKBE  =>  plt_biom%WTHSKBE    , &
    WGSHE    =>  plt_biom%WGSHE      , &
    WTLFBE   =>  plt_biom%WTLFBE     , &
    WTRVE    =>  plt_biom%WTRVE      , &
    WTSTXN   =>  plt_biom%WTSTXN     , &
    WGNODE   =>  plt_biom%WGNODE     , &
    WGNODN   =>  plt_biom%WGNODN     , &
    WGNODP   =>  plt_biom%WGNODP     , &
    FWODLP   =>  plt_allom%FWODLP    , &
    FWODLN   =>  plt_allom%FWODLN    , &
    FWODSN   =>  plt_allom%FWODSN    , &
    FWODB    =>  plt_allom%FWODB     , &
    FWODSP   =>  plt_allom%FWODSP    , &
    GRWTB    => plt_allom%GRWTB      , &
    VRNX     =>  plt_pheno%VRNX      , &
    IWTYP    =>  plt_pheno%IWTYP     , &
    IGTYP    =>  plt_pheno%IGTYP     , &
    VRNF     =>  plt_pheno%VRNF      , &
    VRNL     =>  plt_pheno%VRNL      , &
    ISTYP    =>  plt_pheno%ISTYP     , &
    IFLGF    =>  plt_pheno%IFLGF     , &
    IFLGI    =>  plt_pheno%IFLGI     , &
    IDAY     =>  plt_pheno%IDAY      , &
    IFLGE    =>  plt_pheno%IFLGE     , &
    IBTYP    =>  plt_pheno%IBTYP     , &
    TGSTGI   =>  plt_pheno%TGSTGI    , &
    TGSTGF   =>  plt_pheno%TGSTGF    , &
    VRNS     =>  plt_pheno%VRNS      , &
    GROUP    =>  plt_pheno%GROUP     , &
    VSTGX    =>  plt_pheno%VSTGX     , &
    FLG4     =>  plt_pheno%FLG4      , &
    IFLGR    =>  plt_pheno%IFLGR     , &
    IFLGA    =>  plt_pheno%IFLGA     , &
    GROUPI   =>  plt_pheno%GROUPI    , &
    WSTR     =>   plt_pheno%WSTR     , &
    IFLGQ    =>  plt_pheno%IFLGQ     , &
    KVSTG    =>  plt_pheno%KVSTG     , &
    CFOPC    =>  plt_soilchem%CFOPC  , &
    CFOPN    =>  plt_soilchem%CFOPN  , &
    CFOPP    =>  plt_soilchem%CFOPP  , &
    IYRC     =>  plt_site%IYRC       , &
    ESNC     =>  plt_bgcr%ESNC       , &
    MY       =>  plt_morph%MY        , &
    VSTG     =>  plt_morph%VSTG      , &
    XTLI     =>  plt_morph%XTLI      , &
    HTSHE    =>  plt_morph%HTSHE     , &
    KLEAF    =>  plt_morph%KLEAF     , &
    HTNODX   =>  plt_morph%HTNODX    , &
    HTNODE   =>  plt_morph%HTNODE    , &
    GRNXB    =>  plt_morph%GRNXB     , &
    PSTGF    =>  plt_morph%PSTGF     , &
    NB1      =>  plt_morph%NB1       , &
    NBTB     =>  plt_morph%NBTB      , &
    ARLFB    =>  plt_morph%ARLFB     , &
    ARLF1    =>  plt_morph%ARLF1     , &
    PSTG     =>  plt_morph%PSTG      , &
    PSTGI    =>  plt_morph%PSTGI     , &
    GRNOB    =>  plt_morph%GRNOB       &
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
  IF(ISTYP(NZ).NE.0.OR.(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).GT.1))THEN
    IF((IFLGE(NB,NZ).EQ.0.AND.VRNS(NB,NZ).GE.VRNL(NB,NZ)) &
      .OR.(IFLGF(NB,NZ).EQ.0.AND.VRNF(NB,NZ).GE.VRNX(NB,NZ)))THEN
!
  !    SPRING PHENOLOGY RESET
  !
  !    GROUP,GROUPI=node number required for floral initiation
  !    PSTGI=node number at floral initiation
  !    PSTGF=node number at flowering
  !    VSTGX=leaf number on date of floral initiation
  !    TGSTGI=total change in vegve node number normalized for maturity group
  !    TGSTGF=total change in reprve node number normalized for maturity group
  !    IDAY(1,=emergence date
!
      IF((IFLGE(NB,NZ).EQ.0.AND.ISTYP(NZ).NE.0) &
        .AND.(VRNS(NB,NZ).GE.VRNL(NB,NZ)))THEN
        IF(ISTYP(NZ).EQ.0)THEN
          GROUP(NB,NZ)=AZMAX1(GROUPI(NZ)-NBTB(NB,NZ))
        ELSE
          GROUP(NB,NZ)=GROUPI(NZ)
        ENDIF
        PSTGI(NB,NZ)=PSTG(NB,NZ)
        PSTGF(NB,NZ)=0._r8
        VSTGX(NB,NZ)=0._r8
        TGSTGI(NB,NZ)=0._r8
        TGSTGF(NB,NZ)=0._r8
        IDAY(1,NB,NZ)=I
        DO 2005 M=2,10
          IDAY(M,NB,NZ)=0
2005    CONTINUE
        IF(NB.EQ.NB1(NZ))THEN
          WSTR(NZ)=0._r8
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
        IF(IFLGE(NB,NZ).EQ.0.AND.ISTYP(NZ).NE.0 &
          .AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
          IF(IBTYP(NZ).EQ.0)THEN
            PSTG(NB,NZ)=XTLI(NZ)
            VSTG(NB,NZ)=0._r8
            KLEAF(NB,NZ)=1
            KVSTG(NB,NZ)=1
            FLG4(NB,NZ)=0._r8
            D5330: DO M=1,jsken
              ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ) &
                +CFOPC(5,M,NZ)*WTLFBE(NB,ielmc,NZ)*FWODB(0) &
                +CFOPC(5,M,NZ)*WTSHEBE(NB,ielmc,NZ)*FWODB(0)
              ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ) &
                +CFOPN(5,M,NZ)*WTLFBE(NB,ielmn,NZ)*FWODLN(0) &
                +CFOPN(5,M,NZ)*WTSHEBE(NB,ielmn,NZ)*FWODSN(0)
              ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ) &
                +CFOPP(5,M,NZ)*WTLFBE(NB,ielmp,NZ)*FWODLP(0) &
                +CFOPP(5,M,NZ)*WTSHEBE(NB,ielmp,NZ)*FWODSP(0)
              ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
                +CFOPC(1,M,NZ)*WTLFBE(NB,ielmc,NZ)*FWODB(1) &
                +CFOPC(2,M,NZ)*WTSHEBE(NB,ielmc,NZ)*FWODB(1)
              ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
                +CFOPN(1,M,NZ)*WTLFBE(NB,ielmn,NZ)*FWODLN(1) &
                +CFOPN(2,M,NZ)*WTSHEBE(NB,ielmn,NZ)*FWODSN(1)
              ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
                +CFOPP(1,M,NZ)*WTLFBE(NB,ielmp,NZ)*FWODLP(1) &
                +CFOPP(2,M,NZ)*WTSHEBE(NB,ielmp,NZ)*FWODSP(1)
            ENDDO D5330
            ARLFB(NB,NZ)=0._r8
            WTLFBE(NB,1:npelms,NZ)=0._r8
            WTSHEBE(NB,1:npelms,NZ)=0._r8
            D5335: DO K=0,JNODS1
              ARLF1(K,NB,NZ)=0._r8
              HTSHE(K,NB,NZ)=0._r8
              WGLFE(K,NB,ielmc,NZ)=0._r8
              WSLF(K,NB,NZ)=0._r8
              WGLFE(K,NB,ielmn,NZ)=0._r8
              WGLFE(K,NB,ielmp,NZ)=0._r8
              WGSHE(K,NB,1:npelms,NZ)=0._r8
              WSSHE(K,NB,NZ)=0._r8
            ENDDO D5335
          ENDIF
        ENDIF
    !
    !     RESIDUAL STALKS BECOME LITTERFALL IN GRASSES, SHRUBS AT
    !     START OF SEASON
    !
        IF((IFLGE(NB,NZ).EQ.0.AND.ISTYP(NZ).NE.0).AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
          D6245: DO M=1,jsken
            ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(2,M,NZ) &
              *(WTHSKBE(NB,ielmc,NZ)+WTEARBE(NB,ielmc,NZ)+WTGRBE(NB,ielmc,NZ))
            ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(2,M,NZ) &
              *(WTHSKBE(NB,ielmn,NZ)+WTEARBE(NB,ielmn,NZ)+WTGRBE(NB,ielmn,NZ))
            ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(2,M,NZ) &
              *(WTHSKBE(NB,ielmp,NZ)+WTEARBE(NB,ielmp,NZ)+WTGRBE(NB,ielmp,NZ))
          ENDDO D6245
          WTHSKBE(NB,ielmc,NZ)=0._r8
          WTEARBE(NB,ielmc,NZ)=0._r8
          WTGRBE(NB,ielmc,NZ)=0._r8
          WTHSKBE(NB,ielmn,NZ)=0._r8
          WTEARBE(NB,ielmn,NZ)=0._r8
          WTGRBE(NB,ielmn,NZ)=0._r8
          WTHSKBE(NB,ielmp,NZ)=0._r8
          WTEARBE(NB,ielmp,NZ)=0._r8
          WTGRBE(NB,ielmp,NZ)=0._r8
          GRNXB(NB,NZ)=0._r8
          GRNOB(NB,NZ)=0._r8
          GRWTB(NB,NZ)=0._r8
          IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
            D6345: DO M=1,jsken
              ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(3,M,NZ)*WTSTKBE(NB,ielmc,NZ)
              ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(3,M,NZ)*WTSTKBE(NB,ielmn,NZ)
              ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(3,M,NZ)*WTSTKBE(NB,ielmp,NZ)
            ENDDO D6345
            WTSTKBE(NB,ielmc,NZ)=0._r8
            WTSTKBE(NB,ielmn,NZ)=0._r8
            WTSTKBE(NB,ielmp,NZ)=0._r8
            WTSTXB(NB,NZ)=0._r8
            WTSTXN(NB,NZ)=0._r8
            WTSTXP(NB,NZ)=0._r8
            D6340: DO K=0,JNODS1
              HTNODE(K,NB,NZ)=0._r8
              HTNODX(K,NB,NZ)=0._r8
              WGNODE(K,NB,NZ)=0._r8
              WGNODN(K,NB,NZ)=0._r8
              WGNODP(K,NB,NZ)=0._r8
            ENDDO D6340
          ENDIF
        ENDIF
      ENDIF
!
  !   SPRING OR FALL FLAG RESET
  !
      IF(IFLGE(NB,NZ).EQ.0 &
        .AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
        IFLGE(NB,NZ)=1
        IFLGF(NB,NZ)=0
        IFLGR(NB,NZ)=0
        IFLGQ(NB,NZ)=0
      ELSE
        IFLGE(NB,NZ)=0
        IFLGF(NB,NZ)=1
        IFLGR(NB,NZ)=1
        IFLGQ(NB,NZ)=0
        IFLGA(NB,NZ)=0
      ENDIF
    ENDIF
  ENDIF
!
!   REPRODUCTIVE MATERIAL BECOMES LITTERFALL AT END OF SEASON
!
  IF(IFLGR(NB,NZ).EQ.1)THEN
    IFLGQ(NB,NZ)=IFLGQ(NB,NZ)+1
    IF(IFLGQ(NB,NZ).EQ.IFLGQX)THEN
      IFLGR(NB,NZ)=0
      IFLGQ(NB,NZ)=0
    ENDIF
    D6330: DO M=1,jsken
      ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+FSNR*CFOPC(2,M,NZ) &
        *(WTHSKBE(NB,ielmc,NZ)+WTEARBE(NB,ielmc,NZ))
      ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+FSNR*CFOPN(2,M,NZ) &
        *(WTHSKBE(NB,ielmn,NZ)+WTEARBE(NB,ielmn,NZ))
      ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+FSNR*CFOPP(2,M,NZ) &
        *(WTHSKBE(NB,ielmp,NZ)+WTEARBE(NB,ielmp,NZ))
      IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
        WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+FSNR*CFOPC(2,M,NZ)*WTGRBE(NB,ielmc,NZ)
        WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+FSNR*CFOPN(2,M,NZ)*WTGRBE(NB,ielmn,NZ)
        WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+FSNR*CFOPP(2,M,NZ)*WTGRBE(NB,ielmp,NZ)
      ELSE
        ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+FSNR*CFOPC(2,M,NZ)*WTGRBE(NB,ielmc,NZ)
        ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+FSNR*CFOPN(2,M,NZ)*WTGRBE(NB,ielmn,NZ)
        ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+FSNR*CFOPP(2,M,NZ)*WTGRBE(NB,ielmp,NZ)
      ENDIF
    ENDDO D6330
    WTHSKBE(NB,ielmc,NZ)=(1.0_r8-FSNR)*WTHSKBE(NB,ielmc,NZ)
    WTEARBE(NB,ielmc,NZ)=(1.0_r8-FSNR)*WTEARBE(NB,ielmc,NZ)
    WTGRBE(NB,ielmc,NZ)=(1.0_r8-FSNR)*WTGRBE(NB,ielmc,NZ)
    WTHSKBE(NB,ielmn,NZ)=(1.0_r8-FSNR)*WTHSKBE(NB,ielmn,NZ)
    WTEARBE(NB,ielmn,NZ)=(1.0_r8-FSNR)*WTEARBE(NB,ielmn,NZ)
    WTGRBE(NB,ielmn,NZ)=(1.0_r8-FSNR)*WTGRBE(NB,ielmn,NZ)
    WTHSKBE(NB,ielmp,NZ)=(1.0_r8-FSNR)*WTHSKBE(NB,ielmp,NZ)
    WTEARBE(NB,ielmp,NZ)=(1.0_r8-FSNR)*WTEARBE(NB,ielmp,NZ)
    WTGRBE(NB,ielmp,NZ)=(1.0_r8-FSNR)*WTGRBE(NB,ielmp,NZ)
    GRNXB(NB,NZ)=(1.0_r8-FSNR)*GRNXB(NB,NZ)
    GRNOB(NB,NZ)=(1.0_r8-FSNR)*GRNOB(NB,NZ)
    GRWTB(NB,NZ)=(1.0_r8-FSNR)*GRWTB(NB,NZ)
!
!     STALKS BECOME LITTERFALL IN GRASSES AT END OF SEASON
!
    IF((IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1).AND.ISTYP(NZ).NE.0)THEN
      D6335: DO M=1,jsken
        ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+FSNR*CFOPC(3,M,NZ)*WTSTKBE(NB,ielmc,NZ)
        ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+FSNR*CFOPN(3,M,NZ)*WTSTKBE(NB,ielmn,NZ)
        ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+FSNR*CFOPP(3,M,NZ)*WTSTKBE(NB,ielmp,NZ)
      ENDDO D6335
      WTSTKBE(NB,ielmc,NZ)=(1.0_r8-FSNR)*WTSTKBE(NB,ielmc,NZ)
      WTSTKBE(NB,ielmn,NZ)=(1.0_r8-FSNR)*WTSTKBE(NB,ielmn,NZ)
      WTSTKBE(NB,ielmp,NZ)=(1.0_r8-FSNR)*WTSTKBE(NB,ielmp,NZ)
      WTSTXB(NB,NZ)=(1.0_r8-FSNR)*WTSTXB(NB,NZ)
      WTSTXN(NB,NZ)=(1.0_r8-FSNR)*WTSTXN(NB,NZ)
      WTSTXP(NB,NZ)=(1.0_r8-FSNR)*WTSTXP(NB,NZ)
      D2010: DO K=0,JNODS1
    !     HTNODE(K,NB,NZ)=(1.0_r8-FSNR)*HTNODE(K,NB,NZ)
        HTNODX(K,NB,NZ)=(1.0_r8-FSNR)*HTNODX(K,NB,NZ)
        WGNODE(K,NB,NZ)=(1.0_r8-FSNR)*WGNODE(K,NB,NZ)
        WGNODN(K,NB,NZ)=(1.0_r8-FSNR)*WGNODN(K,NB,NZ)
        WGNODP(K,NB,NZ)=(1.0_r8-FSNR)*WGNODP(K,NB,NZ)
      ENDDO D2010
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
!     IF(J.EQ.INT(ZNOON))THEN
    IF(NB.EQ.NB1(NZ))THEN
      IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
        IDAYH(NZ)=I
        IYRH(NZ)=IYRC
        IHVST(NZ)=1
        JHVST(NZ)=2
        HVST(NZ)=0._r8
        THIN(NZ)=0._r8
        EHVST(1,1,NZ)=1.0_r8
        EHVST(1,2,NZ)=1.0_r8
        EHVST(1,3,NZ)=1.0_r8
        EHVST(1,4,NZ)=1.0_r8
        EHVST(2,1,NZ)=0._r8
        EHVST(2,2,NZ)=1.0_r8
        EHVST(2,3,NZ)=0._r8
        EHVST(2,4,NZ)=0._r8
        IDAY0(NZ)=-1E+06
        IYR0(NZ)=-1E+06
        IFLGI(NZ)=1
      ENDIF
    ENDIF
!     ENDIF
  ENDIF
  end associate
  end subroutine PhenologyReset
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
    IDAY0  =>  plt_distb%IDAY0  , &
    IYR0   =>  plt_distb%IYR0   , &
    IYRC   =>  plt_site%IYRC    , &
    ZEROS2 => plt_site%ZEROS2   , &
    DYLN   =>  plt_site%DYLN    , &
    NU     =>  plt_site%NU      , &
    FVRN   =>  plt_allom%FVRN   , &
    FWODR  =>  plt_allom%FWODR  , &
    WTRTE  =>  plt_biom%WTRTE   , &
    WTRTD  =>  plt_biom%WTRTD   , &
    WTRTL  =>  plt_biom%WTRTL   , &
    WVSTK  =>  plt_biom%WVSTK   , &
    EPOOL  =>  plt_biom%EPOOL   , &
    EPOOLR =>  plt_biom%EPOOLR  , &
    WTRVE  =>  plt_biom%WTRVE   , &
    CEPOLB =>  plt_biom%CEPOLB  , &
    WTLSB  =>  plt_biom%WTLSB   , &
    ZEROP  =>  plt_biom%ZEROP   , &
    WTRSVBE=>  plt_biom%WTRSVBE , &
    WVSTKB =>  plt_biom%WVSTKB  , &
    VOLX   =>  plt_soilchem%VOLX, &
    IDAY   =>  plt_pheno%IDAY   , &
    TFN3   =>  plt_pheno%TFN3   , &
    VRNX   =>  plt_pheno%VRNX   , &
    VRNS   =>  plt_pheno%VRNS   , &
    VRNL   =>  plt_pheno%VRNL   , &
    IFLGI  =>  plt_pheno%IFLGI  , &
    XPPD   =>  plt_pheno%XPPD   , &
    ISTYP  =>  plt_pheno%ISTYP  , &
    IPTYP  =>  plt_pheno%IPTYP  , &
    XDL    =>  plt_pheno%XDL    , &
    VRNF   =>  plt_pheno%VRNF   , &
    IBTYP  =>  plt_pheno%IBTYP  , &
    IGTYP  =>  plt_pheno%IGTYP  , &
    IWTYP  =>  plt_pheno%IWTYP  , &
    IFLGA  =>  plt_pheno%IFLGA  , &
    ATRP   =>   plt_pheno%ATRP  , &
    NG     =>   plt_morph%NG    , &
    NB1    =>  plt_morph%NB1    , &
    NI     =>  plt_morph%NI       &
  )
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!

  IF((ISTYP(NZ).EQ.0.AND.IFLGI(NZ).EQ.0) &
    .OR.(I.GE.IDAY0(NZ).AND.IYRC.EQ.IYR0(NZ) &
    .AND.VRNF(NB,NZ).LT.FVRN(IWTYP(NZ))*VRNX(NB,NZ)) &
    .OR.(VRNS(NB1(NZ),NZ).GE.VRNL(NB,NZ) &
    .AND.VRNF(NB,NZ) &
    .LT.FVRN(IWTYP(NZ))*VRNX(NB,NZ)))THEN
    WTRTM=0._r8
    CPOOLM=0._r8
    D4: DO L=NU,NI(NZ)
      WTRTM=WTRTM+AZMAX1(WTRTD(1,L,NZ))
      CPOOLM=CPOOLM+AZMAX1(EPOOLR(ielmc,1,L,NZ))
    ENDDO D4
!
  ! RESET TIME COUNTER
  !
  ! ATRP=hourly leafout counter
  ! IFLGA=flag for initializing leafout
  !
    IF(IFLGA(NB,NZ).EQ.0)THEN
      ATRP(NB,NZ)=0._r8
      IFLGA(NB,NZ)=1
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
    IF(NB.EQ.NB1(NZ))THEN
      IF(IPTYP(NZ).EQ.2.AND.(IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3))THEN
        PPDX=AZMAX1(XDL(NZ)-XPPD(NZ)-DYLN)
        ATRPPD=EXP(-0.0*PPDX)
      ELSE
        ATRPPD=1.0_r8
      ENDIF
      IF(IGTYP(NZ).NE.0)THEN
        WFNSP=WFNSG
      ELSE
        WFNSP=1.0_r8
      ENDIF
      DATRP=ATRPPD*TFN3(NZ)*WFNSP
      ATRP(NB,NZ)=ATRP(NB,NZ)+DATRP
      IF(ATRP(NB,NZ).LE.ATRPX(ISTYP(NZ)) &
        .OR.(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).EQ.0))THEN
        IF(WTRVE(ielmc,NZ).GT.ZEROP(NZ))THEN
          CPOOLT=CPOOLM+EPOOL(NB,ielmc,NZ)
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
          GFNX=GVMX(ISTYP(NZ))*DATRP
          CH2OH=AZMAX1(GFNX*WTRVE(ielmc,NZ))
          WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)-CH2OH
          EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+CH2OH*FXSH(ISTYP(NZ))
          IF(WTRTM.GT.ZEROP(NZ).AND.CPOOLM.GT.ZEROP(NZ))THEN
            D50: DO L=NU,NI(NZ)
              FXFC=AZMAX1(WTRTD(1,L,NZ))/WTRTM
              EPOOLR(ielmc,1,L,NZ)=EPOOLR(ielmc,1,L,NZ)+FXFC*CH2OH*FXRT(ISTYP(NZ))
            ENDDO D50
          ELSE
            EPOOLR(ielmc,1,NG(NZ),NZ)=EPOOLR(ielmc,1,NG(NZ),NZ)+CH2OH*FXRT(ISTYP(NZ))
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
      IF(WTRVE(ielmc,NZ).GT.ZEROP(NZ))THEN
        IF(ISTYP(NZ).NE.0)THEN
          CPOOLT=AZMAX1(WTRVE(ielmc,NZ)+EPOOL(NB,ielmc,NZ))
          ZPOOLD=(WTRVE(ielmn,NZ)*EPOOL(NB,ielmc,NZ)-EPOOL(NB,ielmn,NZ)*WTRVE(ielmc,NZ))/CPOOLT
          PPOOLD=(WTRVE(ielmp,NZ)*EPOOL(NB,ielmc,NZ)-EPOOL(NB,ielmp,NZ)*WTRVE(ielmc,NZ))/CPOOLT
          UPNH4B=AZMAX1(FRSV(IBTYP(NZ))*ZPOOLD)
          UPPO4B=AZMAX1(FRSV(IBTYP(NZ))*PPOOLD)
        ELSE
          UPNH4B=AZMAX1(FXSH(ISTYP(NZ))*CH2OH*WTRVE(ielmn,NZ)/WTRVE(ielmc,NZ))
          UPPO4B=AZMAX1(FXSH(ISTYP(NZ))*CH2OH*WTRVE(ielmp,NZ)/WTRVE(ielmc,NZ))
        ENDIF
      ELSE
        UPNH4B=AZMAX1(FXSH(ISTYP(NZ))*WTRVE(ielmn,NZ))
        UPPO4B=AZMAX1(FXSH(ISTYP(NZ))*WTRVE(ielmp,NZ))
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
      D3: DO L=NU,NI(NZ)
        CPOOLM=CPOOLM+AZMAX1(EPOOLR(ielmc,1,L,NZ))
        ZPOOLM=ZPOOLM+AZMAX1(EPOOLR(ielmn,1,L,NZ))
        PPOOLM=PPOOLM+AZMAX1(EPOOLR(ielmp,1,L,NZ))
      ENDDO D3
      IF(WTRVE(ielmc,NZ).GT.ZEROP(NZ))THEN
        IF(ISTYP(NZ).NE.0)THEN
          CPOOLT=AMAX1(ZEROP(NZ),WTRVE(ielmc,NZ)+CPOOLM)
          ZPOOLD=(WTRVE(ielmn,NZ)*CPOOLM-ZPOOLM*WTRVE(ielmc,NZ))/CPOOLT
          PPOOLD=(WTRVE(ielmp,NZ)*CPOOLM-PPOOLM*WTRVE(ielmc,NZ))/CPOOLT
          UPNH4R=AZMAX1(FRSV(IBTYP(NZ))*ZPOOLD)
          UPPO4R=AZMAX1(FRSV(IBTYP(NZ))*PPOOLD)
        ELSE
          UPNH4R=AZMAX1(FXRT(ISTYP(NZ))*CH2OH*WTRVE(ielmn,NZ)/WTRVE(ielmc,NZ))
          UPPO4R=AZMAX1(FXRT(ISTYP(NZ))*CH2OH*WTRVE(ielmp,NZ)/WTRVE(ielmc,NZ))
        ENDIF
      ELSE
        UPNH4R=AZMAX1(FXRT(ISTYP(NZ))*WTRVE(ielmn,NZ))
        UPPO4R=AZMAX1(FXRT(ISTYP(NZ))*WTRVE(ielmp,NZ))
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
      WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)-UPNH4B-UPNH4R
      WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)-UPPO4B-UPPO4R
      EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+UPNH4B
      EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+UPPO4B
      IF(WTRTM.GT.ZEROP(NZ).AND.CPOOLM.GT.ZEROP(NZ))THEN
        D51: DO L=NU,NI(NZ)
          FXFN=AZMAX1(EPOOLR(ielmc,1,L,NZ))/CPOOLM

          EPOOLR(ielmn,1,L,NZ)=EPOOLR(ielmn,1,L,NZ)+FXFN*UPNH4R
          EPOOLR(ielmp,1,L,NZ)=EPOOLR(ielmp,1,L,NZ)+FXFN*UPPO4R
        ENDDO D51
      ELSE

        EPOOLR(ielmn,1,NG(NZ),NZ)=EPOOLR(ielmn,1,NG(NZ),NZ)+UPNH4R
        EPOOLR(ielmp,1,NG(NZ),NZ)=EPOOLR(ielmp,1,NG(NZ),NZ)+UPPO4R
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
    IF(NB.NE.NB1(NZ).AND.ATRP(NB,NZ) &
      .LE.ATRPX(ISTYP(NZ)))THEN
      ATRP(NB,NZ)=ATRP(NB,NZ)+TFN3(NZ)*WFNG
      XFRC=AZMAX1(0.05*TFN3(NZ) &
        *(0.5*(EPOOL(NB1(NZ),ielmc,NZ)+EPOOL(NB,ielmc,NZ))-EPOOL(NB,ielmc,NZ)))
      XFRN=AZMAX1(0.05*TFN3(NZ) &
        *(0.5*(EPOOL(NB1(NZ),ielmn,NZ)+EPOOL(NB,ielmn,NZ))-EPOOL(NB,ielmn,NZ)))
      XFRP=AZMAX1(0.05*TFN3(NZ) &
      *(0.5*(EPOOL(NB1(NZ),ielmp,NZ)+EPOOL(NB,ielmp,NZ))-EPOOL(NB,ielmp,NZ)))
      EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)+XFRC
      EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)+XFRN
      EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)+XFRP
      EPOOL(NB1(NZ),ielmc,NZ)=EPOOL(NB1(NZ),ielmc,NZ)-XFRC
      EPOOL(NB1(NZ),ielmn,NZ)=EPOOL(NB1(NZ),ielmn,NZ)-XFRN
      EPOOL(NB1(NZ),ielmp,NZ)=EPOOL(NB1(NZ),ielmp,NZ)-XFRP
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
  IF(IFLGZ.EQ.1.AND.ISTYP(NZ).NE.0)THEN
    IF(WVSTKB(NB,NZ).GT.ZEROP(NZ) &
      .AND.WTRSVBE(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
      CWTRSV=AZMAX1(WTRSVBE(NB,ielmc,NZ)/WVSTKB(NB,NZ))
      CWTRSN=AZMAX1(WTRSVBE(NB,ielmn,NZ)/WVSTKB(NB,NZ))
      CWTRSP=AZMAX1(WTRSVBE(NB,ielmp,NZ)/WVSTKB(NB,NZ))
      CNR=CWTRSV/(CWTRSV+CWTRSN/CNKI)
      CPR=CWTRSV/(CWTRSV+CWTRSP/CPKI)
    ELSE
      CNR=0._r8
      CPR=0._r8
    ENDIF
    XFRCX=FXFB(IBTYP(NZ))*AZMAX1(WTRSVBE(NB,ielmc,NZ))
    XFRNX=FXFB(IBTYP(NZ))*AZMAX1(WTRSVBE(NB,ielmn,NZ))*(1.0+CNR)
    XFRPX=FXFB(IBTYP(NZ))*AZMAX1(WTRSVBE(NB,ielmp,NZ))*(1.0+CPR)
    XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
    WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)-XFRC
    WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+XFRC
    WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)-XFRN
    WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+XFRN
    WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)-XFRP
    WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+XFRP
    IF(CEPOLB(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
      CNL=CEPOLB(NB,ielmc,NZ)/(CEPOLB(NB,ielmc,NZ)+CEPOLB(NB,ielmn,NZ)/CNKI)
      CPL=CEPOLB(NB,ielmc,NZ)/(CEPOLB(NB,ielmc,NZ)+CEPOLB(NB,ielmp,NZ)/CPKI)
    ELSE
      CNL=0._r8
      CPL=0._r8
    ENDIF
    XFRCX=FXFB(IBTYP(NZ))*AZMAX1(EPOOL(NB,ielmc,NZ))
    XFRNX=FXFB(IBTYP(NZ))*AZMAX1(EPOOL(NB,ielmn,NZ))*(1.0+CNL)
    XFRPX=FXFB(IBTYP(NZ))*AZMAX1(EPOOL(NB,ielmp,NZ))*(1.0+CPL)
    XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
    EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)-XFRC
    WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+XFRC
    EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)-XFRN
    WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+XFRN
    EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)-XFRP
    WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+XFRP
  ENDIF
!
!   TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES AND ROOTS TO RESERVES
!   IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN STALK RESERVES
!   AND LEAVES IN PERENNIALS ACCORDING TO CONCENTRATION DIFFERENCES
!
!   ISTYP=growth habit:0=annual,1=perennial from PFT file
!   IDAY(3,=start of stem elongation and setting max seed number
!   IDAY(8,=end date setting for final seed number
!   WTLSB=leaf+petiole mass
!   WVSTKB=stalk sapwood mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   FXFY,FXFZ=rate constant for plant-reserve nonstructural C,N,P exchange
!   XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!   CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  IF((ISTYP(NZ).EQ.0.AND.IDAY(8,NB,NZ).NE.0) &
    .OR.(ISTYP(NZ).EQ.1.AND.IDAY(3,NB,NZ).NE.0))THEN
    WTPLTT=WTLSB(NB,NZ)+WVSTKB(NB,NZ)
    CPOOLT=EPOOL(NB,ielmc,NZ)+WTRSVBE(NB,ielmc,NZ)
    IF(WTPLTT.GT.ZEROP(NZ))THEN
      CPOOLD=(EPOOL(NB,ielmc,NZ)*WVSTKB(NB,NZ) &
        -WTRSVBE(NB,ielmc,NZ)*WTLSB(NB,NZ))/WTPLTT
      XFRC=FXFY(ISTYP(NZ))*CPOOLD
      EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)-XFRC
      WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+XFRC
    ENDIF
    IF(CPOOLT.GT.ZEROP(NZ))THEN
      ZPOOLD=(EPOOL(NB,ielmn,NZ)*WTRSVBE(NB,ielmc,NZ) &
        -WTRSVBE(NB,ielmn,NZ)*EPOOL(NB,ielmc,NZ))/CPOOLT
      PPOOLD=(EPOOL(NB,ielmp,NZ)*WTRSVBE(NB,ielmc,NZ) &
        -WTRSVBE(NB,ielmp,NZ)*EPOOL(NB,ielmc,NZ))/CPOOLT
      XFRN=FXFZ(ISTYP(NZ))*ZPOOLD
      XFRP=FXFZ(ISTYP(NZ))*PPOOLD
      EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)-XFRN
      WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)+XFRN
      EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)-XFRP
      WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)+XFRP
    ENDIF

    IF(ISTYP(NZ).EQ.0.AND.IDAY(8,NB,NZ).NE.0)THEN
      D2050: DO L=NU,NI(NZ)
        IF(VOLX(L).GT.ZEROS2)THEN
          WTRTRX=AMAX1(ZEROP(NZ),WTRTL(1,L,NZ)*FWODR(1))
          WTPLTX=WTRTRX+WVSTKB(NB,NZ)
          IF(WTPLTX.GT.ZEROP(NZ))THEN
            CPOOLD=(EPOOLR(ielmc,1,L,NZ)*WVSTKB(NB,NZ)-WTRSVBE(NB,ielmc,NZ)*WTRTRX)/WTPLTX
            XFRC=AZMAX1(FXFY(ISTYP(NZ))*CPOOLD)
            EPOOLR(ielmc,1,L,NZ)=EPOOLR(ielmc,1,L,NZ)-XFRC
            WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+XFRC
            CPOOLT=EPOOLR(ielmc,1,L,NZ)+WTRSVBE(NB,ielmc,NZ)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              ZPOOLD=(EPOOLR(ielmn,1,L,NZ)*WTRSVBE(NB,ielmc,NZ)-WTRSVBE(NB,ielmn,NZ)*EPOOLR(ielmc,1,L,NZ))/CPOOLT
              PPOOLD=(EPOOLR(ielmp,1,L,NZ)*WTRSVBE(NB,ielmc,NZ)-WTRSVBE(NB,ielmp,NZ)*EPOOLR(ielmc,1,L,NZ))/CPOOLT
              XFRN=AZMAX1(FXFZ(ISTYP(NZ))*ZPOOLD)
              XFRP=AZMAX1(FXFZ(ISTYP(NZ))*PPOOLD)
              EPOOLR(ielmn,1,L,NZ)=EPOOLR(ielmn,1,L,NZ)-XFRN
              WTRSVBE(NB,ielmn,NZ)=WTRSVBE(NB,ielmn,NZ)+XFRN
              EPOOLR(ielmp,1,L,NZ)=EPOOLR(ielmp,1,L,NZ)-XFRP
              WTRSVBE(NB,ielmp,NZ)=WTRSVBE(NB,ielmp,NZ)+XFRP

            ENDIF
          ENDIF
        ENDIF
      ENDDO D2050
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
  IF(WVSTKB(NB,NZ).GT.ZEROP(NZ) &
    .AND.WVSTK(NZ).GT.ZEROP(NZ) &
    .AND.WTRTE(ielmc,NZ).GT.ZEROP(NZ) &
    .AND.WTRSVBE(NB,ielmc,NZ).LE.XFRX*WVSTKB(NB,NZ))THEN
    FWTBR=WVSTKB(NB,NZ)/WVSTK(NZ)
    WVSTBX=WVSTKB(NB,NZ)
    WTRTTX=WTRTE(ielmc,NZ)*FWTBR
    WTPLTT=WVSTBX+WTRTTX
    WTRSBX=AZMAX1(WTRSVBE(NB,ielmc,NZ))
    WTRVCX=AZMAX1(WTRVE(ielmc,NZ)*FWTBR)
    CPOOLD=(WTRVCX*WVSTBX-WTRSBX*WTRTTX)/WTPLTT
    XFRC=AZMAX1(XFRY*CPOOLD)
    WTRSVBE(NB,ielmc,NZ)=WTRSVBE(NB,ielmc,NZ)+XFRC
    WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)-XFRC
  ENDIF
  end associate
  end subroutine CarbNutInBranchTransfer
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
  associate(                             &
    CNET      =>  plt_bgcr%CNET    , &
    TCO2T     =>  plt_bgcr%TCO2T   , &
    TRAU      =>  plt_bgcr%TRAU    , &
    RECO      =>  plt_bgcr%RECO    , &
    TGPP      =>  plt_bgcr%TGPP    , &
    CARBN     =>  plt_bgcr%CARBN   , &
    TCO2A     =>  plt_bgcr%TCO2A   , &
    EPOOL     =>  plt_biom%EPOOL   , &
    CEPOLB    =>  plt_biom%CEPOLB  , &
    ZERO      => plt_site%ZERO     , &
    IGTYP     =>  plt_pheno%IGTYP  , &
    TFN3      =>  plt_pheno%TFN3   , &
    IWTYP     =>  plt_pheno%IWTYP  , &
    FDBKX     =>  plt_photo%FDBKX    &
  )
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
!
  IF(CEPOLB(NB,ielmc,NZ).GT.ZERO)THEN
    CNPG=AMIN1(CEPOLB(NB,ielmn,NZ)/(CEPOLB(NB,ielmn,NZ) &
      +CEPOLB(NB,ielmc,NZ)*CNKI),CEPOLB(NB,ielmp,NZ)/(CEPOLB(NB,ielmp,NZ) &
      +CEPOLB(NB,ielmc,NZ)*CPKI))
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
  RCO2C=AZMAX1(VMXC*EPOOL(NB,ielmc,NZ) &
    *TFN3(NZ))*CNPG*FDBKX(NB,NZ)*WFNG
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
  RMNCS=AZMAX1(RMPLT*TFN5*WTSHXN)
  IF(IGTYP(NZ).EQ.0.OR.IWTYP(NZ).EQ.2)THEN
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
  RCO2Y=AZMAX1(RCO2X)*WFNSG
  SNCR=AZMAX1(-RCO2X)
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
    ZPOOLB=AZMAX1(EPOOL(NB,ielmn,NZ))
    PPOOLB=AZMAX1(EPOOL(NB,ielmp,NZ))
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
  ZADDB=AZMAX1(AMIN1(EPOOL(NB,ielmn,NZ),CGROS*(CNSHX+CNLFM+CNLFX*CNPG)))
  PADDB=AZMAX1(AMIN1(EPOOL(NB,ielmp,NZ) &
    ,CGROS*(CPSHX+CPLFM+CPLFX*CNPG)))
  CNRDA=AZMAX1(1.70*ZADDB-0.025*CH2O)
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
  CARBN(NZ)=CARBN(NZ)+CO2F
  TCO2T(NZ)=TCO2T(NZ)-RCO2T
  TCO2A(NZ)=TCO2A(NZ)-RCO2T
  CNET(NZ)=CNET(NZ)+CO2F-RCO2T
  TGPP=TGPP+CO2F
  RECO=RECO-RCO2T
  TRAU=TRAU-RCO2T

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
  associate(                          &
    CEPOLB    =>  plt_biom%CEPOLB   , &
    EPOOL     =>  plt_biom%EPOOL    , &
    IGTYP     =>  plt_pheno%IGTYP   , &
    TFN4      =>  plt_pheno%TFN4    , &
    WFR       =>  plt_rbgc%WFR      , &
    IWTYP     =>  plt_pheno%IWTYP   , &
    CNET      =>  plt_bgcr%CNET     , &
    RCO2M     =>  plt_rbgc%RCO2M    , &
    RCO2N     =>  plt_rbgc%RCO2N    , &
    RCO2A     =>  plt_rbgc%RCO2A    , &
    ZERO      => plt_site%ZERO      , &
    NG        =>  plt_morph%NG      , &
    FDBKX     =>  plt_photo%FDBKX     &
  )
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  IF(CEPOLB(NB,ielmc,NZ).GT.ZERO)THEN
    CNPG=AMIN1(CEPOLB(NB,ielmn,NZ)/(CEPOLB(NB,ielmn,NZ)+CEPOLB(NB,ielmc,NZ)*CNKI), &
      CEPOLB(NB,ielmp,NZ)/(CEPOLB(NB,ielmp,NZ)+CEPOLB(NB,ielmc,NZ)*CPKI))
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
  RCO2CM=AZMAX1(VMXC*EPOOL(NB,ielmc,NZ) &
    *TFN4(NG(NZ),NZ))*CNPG*WFNG*FDBKX(NB,NZ)
  RCO2C=RCO2CM*WFR(1,NG(NZ),NZ)
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
  RMNCS=AZMAX1(RMPLT*TFN6(NG(NZ))*WTSHXN)
  IF(IGTYP(NZ).EQ.0.OR.IWTYP(NZ).EQ.2)THEN
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
  RCO2YM=AZMAX1(RCO2XM)*WFNSG
  RCO2Y=AZMAX1(RCO2X)*WFNSG
  SNCRM=AZMAX1(-RCO2XM)
  SNCR=AZMAX1(-RCO2X)
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
    ZPOOLB=AZMAX1(EPOOL(NB,ielmn,NZ))
    PPOOLB=AZMAX1(EPOOL(NB,ielmp,NZ))
    FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
    IF(RCO2YM.GT.0.0_r8)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
    ELSE
      RCO2GM=0._r8
    ENDIF
    IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(1,NG(NZ),NZ))
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
  ZADDBM=AZMAX1(CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
  ZADDB=AZMAX1(CGROS*(CNSHX+CNLFM+CNLFX*CNPG))
  PADDB=AZMAX1(CGROS*(CPSHX+CPLFM+CPLFX*CNPG))
  CNRDM=AZMAX1(1.70_r8*ZADDBM)
  CNRDA=AZMAX1(1.70_r8*ZADDB)
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
  RCO2M(1,NG(NZ),NZ)=RCO2M(1,NG(NZ),NZ)+RCO2TM
  RCO2N(1,NG(NZ),NZ)=RCO2N(1,NG(NZ),NZ)+RCO2T
  RCO2A(1,NG(NZ),NZ)=RCO2A(1,NG(NZ),NZ)-RCO2T
  CH2O=0._r8
  end associate
  end subroutine ComputeRAutoBfEmergence

end module PlantBranchMod
