module PlantBranchMod

! Description:
! module for plant biological transformations
  use minimathmod, only : isclose,safe_adb,AZMAX1
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  use PhotoSynsMod
  use RootMod, only : RootBGCModel
  use NoduleBGCMod
  use LitterFallMod
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: GrowOneBranch
  contains
!------------------------------------------------------------------------------------------

  subroutine GrowOneBranch(I,J,NB,NZ,TFN6,ZCX,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WFNG,&
    Stomata_Activity,WFNS,WFNSG,PTRT,UPNFC,IFLGZ)
  implicit none
  integer, intent(in)  :: I,J,NB,NZ
  REAL(R8), INTENT(IN) :: TFN6(JZ1)
  real(r8), intent(in) :: ZCX(JP1)
  real(r8), intent(in) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG
  real(r8), intent(inout) :: UPNFC(JP1)
  real(r8), intent(out) :: PTRT
  integer, intent(out) :: IFLGZ
  real(r8) :: DMSHD
  integer  :: K,KNOD,KK,kx,K1,K2,KVSTGX,KSNC
  integer  :: KN,MNNOD,NN,MXNOD,M,N,NNOD1,NE
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
  real(r8) :: GROGRE(NumOfPlantChemElements)
  real(r8) :: GROLFE(NumOfPlantChemElements)
  real(r8) :: GRORSVE(NumOfPlantChemElements),GROHSKE(NumOfPlantChemElements)
  real(r8) :: GROEARE(NumOfPlantChemElements)
  real(r8) :: GROSHT,GROSHE(NumOfPlantChemElements)
  real(r8) :: GROSTKE(NumOfPlantChemElements)
  real(r8) :: GNOD
  REAL(R8) :: GROE(NumOfPlantChemElements)
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
  real(r8) :: XFRE(1:NumOfPlantChemElements)
! begin_execution
  associate(                            &
    instruct =>  pltpar%instruct  , &
    ifoliar  =>  pltpar%ifoliar   , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    infoliar =>  pltpar%infoliar  , &
    icwood   =>  pltpar%icwood    , &
    k_woody_litr=> pltpar%k_woody_litr, &
    k_fine_litr=> pltpar%k_fine_litr  , &
    DMLF       =>  plt_allom%DMLF     , &
    DMSHE      =>  plt_allom%DMSHE    , &
    DMEAR      =>  plt_allom%DMEAR    , &
    CNEAR      =>  plt_allom%CNEAR    , &
    CPRSV      =>  plt_allom%CPRSV    , &
    CPHSK      =>  plt_allom%CPHSK    , &
    CPSTK      =>  plt_allom%CPSTK    , &
    CNSTK      =>  plt_allom%CNSTK    , &
    DMSTK      =>  plt_allom%DMSTK    , &
    DMHSK      =>  plt_allom%DMHSK    , &
    FWODBE     =>  plt_allom%FWODBE   , &
    FWODLE     =>  plt_allom%FWODLE   , &
    CNHSK      =>  plt_allom%CNHSK    , &
    CPEAR      =>  plt_allom%CPEAR    , &
    DMRSV      =>  plt_allom%DMRSV    , &
    DMGR       =>  plt_allom%DMGR     , &
    CNWS       =>  plt_allom%CNWS     , &
    CPWS       =>  plt_allom%CPWS     , &
    CNRSV      =>  plt_allom%CNRSV    , &
    BiomGrowthYieldRoot       =>  plt_allom%BiomGrowthYieldRoot     , &
    FNOD       =>  plt_allom%FNOD     , &
    ESNC       =>  plt_bgcr%ESNC      , &
    CFOPE      =>  plt_soilchem%CFOPE , &
    WTSTKBE    =>  plt_biom%WTSTKBE   , &
    WTRSVBE    =>  plt_biom%WTRSVBE   , &
    WTHSKBE    =>  plt_biom%WTHSKBE   , &
    WTEARBE    =>  plt_biom%WTEARBE   , &
    WTSTXBE    =>  plt_biom%WTSTXBE   , &
    WGSHEXE    =>  plt_biom%WGSHEXE   , &
    CanPBStalkC     =>  plt_biom%CanPBStalkC    , &
    EPOOL      =>  plt_biom%EPOOL     , &
    CanPBLeafShethC      =>  plt_biom%CanPBLeafShethC     , &
    WGLFEX     =>  plt_biom%WGLFEX    , &
    WTSHEBE    =>  plt_biom%WTSHEBE   , &
    WTLFBE     =>  plt_biom%WTLFBE    , &
    WSLF       =>  plt_biom%WSLF      , &
    WGLFE      =>  plt_biom%WGLFE     , &
    WGSHE      =>  plt_biom%WGSHE     , &
    CEPOLB     =>  plt_biom%CEPOLB    , &
    WGNODE     =>  plt_biom%WGNODE    , &
    WSSHE      =>  plt_biom%WSSHE     , &
    WTGRBE     =>  plt_biom%WTGRBE    , &
    ZEROP      =>  plt_biom%ZEROP     , &
    ZEROL      =>  plt_biom%ZEROL     , &
    ICTYP      =>  plt_photo%ICTYP    , &
    IDTHB      =>  plt_pheno%IDTHB    , &
    KVSTGN     =>  plt_pheno%KVSTGN   , &
    XRLA       =>  plt_pheno%XRLA     , &
    RCELX      =>  plt_pheno%RCELX    , &
    KVSTG      =>  plt_pheno%KVSTG    , &
    IGTYP      =>  plt_pheno%IGTYP    , &
    fTgrowCanP       =>  plt_pheno%fTgrowCanP     , &
    IFLGP      =>  plt_pheno%IFLGP    , &
    FLGZ       =>  plt_pheno%FLGZ     , &
    IDAY       =>  plt_pheno%IDAY     , &
    RCESX      =>  plt_pheno%RCESX    , &
    IBTYP      =>  plt_pheno%IBTYP    , &
    IFLGG      =>  plt_pheno%IFLGG    , &
    ISTYP      =>  plt_pheno%ISTYP    , &
    CTC        =>  plt_pheno%CTC      , &
    SSINN      =>  plt_rad%SSINN      , &
    pftPlantPopulation         =>  plt_site%pftPlantPopulation        , &
    ZERO       =>  plt_site%ZERO      , &
    PSICanP      =>  plt_ew%PSICanP       , &
    NB1        =>  plt_morph%NB1      , &
    HypoctoylHeight      =>   plt_morph%HypoctoylHeight   , &
    BranchNumber_brchpft       =>   plt_morph%BranchNumber_brchpft    , &
    NumOfBranches_pft        =>   plt_morph%NumOfBranches_pft     , &
    SLA1       =>   plt_morph%SLA1    , &
    SNL1       =>   plt_morph%SNL1    , &
    ANGSH      =>   plt_morph%ANGSH   , &
    ARLFZ      =>   plt_morph%ARLFZ   , &
    CanopyBranchLeafA_pft      =>   plt_morph%CanopyBranchLeafA_pft   , &
    SeedinDepth      =>   plt_morph%SeedinDepth   , &
    CanPBranchHeight     =>   plt_morph%CanPBranchHeight  , &
    ARLF1      =>   plt_morph%ARLF1   , &
    CanPSheathHeight      =>   plt_morph%CanPSheathHeight   , &
    ANGBR      =>   plt_morph%ANGBR   , &
    SSL1       =>   plt_morph%SSL1    , &
    HTNODE     =>   plt_morph%HTNODE  , &
    HTNODX     =>   plt_morph%HTNODX  , &
    NNOD       =>   plt_morph%NNOD      &
  )
  CanPBLeafShethC(NB,NZ)=AZMAX1(WTLFBE(ielmc,NB,NZ)+WTSHEBE(ielmc,NB,NZ))

  IF(IDTHB(NB,NZ).EQ.ibralive)THEN
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
      DMLFB=BiomGrowthYieldRoot(NZ)
      DMSHB=BiomGrowthYieldRoot(NZ)
      CNLFB=CNRTW
      CNSHB=CNRTW
      CPLFB=CPRTW
      CPSHB=CPRTW
    ENDIF
!part 1(leaf), (2) sheath, (3) stalk, (4) reserve, (5) husk, (6) ear, (7) grain
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
!   CNSTK,CanPBStalkC=stalk N:C,sapwood mass
!   IDAY(10=date of physiological maturity
!
    WTSHXN=AZMAX1(WTLFBE(ielmn,NB,NZ)+WTSHEBE(ielmn,NB,NZ) &
      +CNSTK(NZ)*CanPBStalkC(NB,NZ))
    IF(IDAY(10,NB,NZ).EQ.0)THEN
      WTSHXN=WTSHXN+AZMAX1(WTHSKBE(ielmn,NB,NZ) &
        +WTEARBE(ielmn,NB,NZ)+WTGRBE(ielmn,NB,NZ))
    ENDIF
!
!   GROSS PRIMARY PRODUCTIVITY
!
    call UpdatePhotosynthates(I,J,NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,&
      CNLFX,CPLFX,WTSHXN,TFN5,WFNG,Stomata_Activity,WFNSG,CH2O3,CH2O4,CNPG,rco2c,RMNCS,&
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
    IF(ICTYP(NZ).EQ.ic4_photo)THEN
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
    GROLFE(ielmc)=PART(1)*CGROS*DMLFB
    GROSHE(ielmc)=PART(2)*CGROS*DMSHB
    GROSTKE(ielmc)=PART(3)*CGROS*DMSTK(NZ)
    GRORSVE(ielmc)=PART(4)*CGROS*DMRSV(NZ)
    GROHSKE(ielmc)=PART(5)*CGROS*DMHSK(NZ)
    GROEARE(ielmc)=PART(6)*CGROS*DMEAR(NZ)
    GROGRE(ielmc)=PART(7)*CGROS*DMGR(NZ)
    GROSHT=CGROS*DMSHT

    GROLFE(ielmn)=GROLFE(ielmc)*CNLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHE(ielmn)=GROSHE(ielmc)*CNSHB
    GROSTKE(ielmn)=GROSTKE(ielmc)*CNSTK(NZ)
    GRORSVE(ielmn)=GRORSVE(ielmc)*CNRSV(NZ)
    GROHSKE(ielmn)=GROHSKE(ielmc)*CNHSK(NZ)
    GROEARE(ielmn)=GROEARE(ielmc)*CNEAR(NZ)
    GROGRE(ielmn)=GROGRE(ielmc)*CNRSV(NZ)

    GROLFE(ielmp)=GROLFE(ielmc)*CPLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHE(ielmp)=GROSHE(ielmc)*CPSHB
    GROSTKE(ielmp)=GROSTKE(ielmc)*CPSTK(NZ)
    GRORSVE(ielmp)=GRORSVE(ielmc)*CPRSV(NZ)
    GROHSKE(ielmp)=GROHSKE(ielmc)*CPHSK(NZ)
    GROEARE(ielmp)=GROEARE(ielmc)*CPEAR(NZ)
    GROGRE(ielmp)=GROGRE(ielmc)*CPRSV(NZ)

    DO NE=1,NumOfPlantChemElements
      WTLFBE(NE,NB,NZ)=WTLFBE(NE,NB,NZ)+GROLFE(NE)
      WTSHEBE(NE,NB,NZ)=WTSHEBE(NE,NB,NZ)+GROSHE(NE)
      WTSTKBE(NE,NB,NZ)=WTSTKBE(NE,NB,NZ)+GROSTKE(NE)
      WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)+GRORSVE(NE)
      WTHSKBE(NE,NB,NZ)=WTHSKBE(NE,NB,NZ)+GROHSKE(NE)
      WTEARBE(NE,NB,NZ)=WTEARBE(NE,NB,NZ)+GROEARE(NE)
    ENDDO

!
!   ETOLIATION
!
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
!   ETOL=coefficient for etoliation effects on expansion,extension
!
    CCE=AMIN1(safe_adb(CEPOLB(ielmn,NB,NZ),CEPOLB(ielmn,NB,NZ)+CEPOLB(ielmc,NB,NZ)*CNKI) &
      ,safe_adb(CEPOLB(ielmp,NB,NZ),CEPOLB(ielmp,NB,NZ)+CEPOLB(ielmc,NB,NZ)*CPKI))

    ETOL=1.0_r8+CCE
!
!   DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
!
!   MXNOD,MNNOD=max,min node number currently growing
!   KVSTG=integer of most recent leaf number
!   KNOD,GNOD=number of currently growing nodes
!   ALLOCL=fraction of leaf growth allocated to each node
!   GRO,GROE(ielmn),GROE(ielmp)=leaf C,N,P growth at each node
!   GSLA=allocation of leaf area growth to each node
!   FNOD=scales node number for perennial vegetation (e.g. trees)
!   NNOD=number of concurrently growing nodes
!
    IF(NB.EQ.NB1(NZ).AND.HypoctoylHeight(NZ).LE.SeedinDepth(NZ))THEN
      NNOD1=0
    ELSE
      NNOD1=1
    ENDIF
    IF(GROLFE(ielmc).GT.0.0_r8)THEN
      MXNOD=KVSTG(NB,NZ)
      MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ)+1)
      MXNOD=MAX(MXNOD,MNNOD)
      KNOD=MXNOD-MNNOD+1
      GNOD=KNOD
      ALLOCL=1.0_r8/GNOD
      DO NE=1,NumOfPlantChemElements
        GROE(NE)=ALLOCL*GROLFE(NE)
      ENDDO
      GSLA=ALLOCL*FNOD(NZ)*NNOD(NZ)
!
!     GROWTH AT EACH CURRENT NODE
!
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     GRO,GROE(ielmn),GROE(ielmp)=leaf C,N,P growth at each node
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
      D490: DO KK=MNNOD,MXNOD
        K=MOD(KK,JNODS1)
        IF(K.EQ.0.AND.KK.NE.0)K=25
          DO NE=1,NumOfPlantChemElements
            WGLFE(NE,K,NB,NZ)=WGLFE(NE,K,NB,NZ)+GROE(NE)
          ENDDO
          WSLF(K,NB,NZ)=WSLF(K,NB,NZ)+AMIN1(GROE(ielmn)*CNWS(NZ),GROE(ielmp)*CPWS(NZ))
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
!         CanopyBranchLeafA_pft,ARLF=branch,node leaf area
!
          SLA=ETOL*SLA1(NZ)*(AMAX1(ZEROL(NZ) &
            ,WGLFE(ielmc,K,NB,NZ))/(pftPlantPopulation(NZ)*GSLA))**SLA2*WFNS
          GROA=GROE(ielmc)*SLA
          CanopyBranchLeafA_pft(NB,NZ)=CanopyBranchLeafA_pft(NB,NZ)+GROA
          ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)+GROA
        ENDDO D490
      ENDIF
!
!     DISTRIBUTE SHEATH OR PETIOLE GROWTH AMONG CURRENTLY GROWING NODES
!
!     MXNOD,MNNOD=max,min node number currently growing
!     KVSTG=integer of most recent leaf number
!     GNOD=number of currently growing nodes
!     ALLOCS=fraction of petiole growth allocated to each node
!     GRO,GROE(ielmn),GROE(ielmp)=petiole C,N,P growth at each node
!     GSSL=allocation of petiole length growth to each node
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
      IF(GROSHE(ielmc).GT.0.0)THEN
        MXNOD=KVSTG(NB,NZ)
        MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ)+1)
        MXNOD=MAX(MXNOD,MNNOD)
        GNOD=MXNOD-MNNOD+1
        ALLOCS=1.0_r8/GNOD
        DO NE=1,NumOfPlantChemElements
          GROE(NE)=ALLOCS*GROSHE(NE)
        ENDDO
        GSSL=ALLOCL*FNOD(NZ)*NNOD(NZ)
!
!       GROWTH AT EACH CURRENT NODE
!
!       WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!       GRO,GROE(ielmn),GROE(ielmp)=petiole C,N,P growth at each node
!       CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
        D505: DO KK=MNNOD,MXNOD
          K=MOD(KK,JNODS1)
          IF(K.EQ.0.AND.KK.NE.0)K=25
            DO NE=1,NumOfPlantChemElements
              WGSHE(NE,K,NB,NZ)=WGSHE(NE,K,NB,NZ)+GROE(NE)
            ENDDO
            WSSHE(K,NB,NZ)=WSSHE(K,NB,NZ) &
              +AMIN1(GROE(ielmn)*CNWS(NZ),GROE(ielmp)*CPWS(NZ))
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
        !   CanPSheathHeight=petiole length
!
            IF(WGLFE(ielmc,K,NB,NZ).GT.0.0_r8)THEN
              SSL=ETOL*SSL1(NZ)*(AMAX1(ZEROL(NZ) &
                ,WGSHE(ielmc,K,NB,NZ))/(pftPlantPopulation(NZ)*GSSL))**SSL2*WFNS
              GROS=GROE(ielmc)/pftPlantPopulation(NZ)*SSL
              CanPSheathHeight(K,NB,NZ)=CanPSheathHeight(K,NB,NZ)+GROS*ANGSH(NZ)
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
    !   GRO,GROE(ielmn),GROE(ielmp)=stalk C,N,P growth at each node
!
        IF(IDAY(1,NB,NZ).EQ.0)THEN
          NN=0
        ELSE
          NN=1
        ENDIF
        MXNOD=KVSTG(NB,NZ)
        MNNOD=MAX(MIN(NN,MAX(NN,MXNOD-NNOD(NZ))),KVSTG(NB,NZ)-23)
        MXNOD=MAX(MXNOD,MNNOD)
        IF(GROSTKE(ielmc).GT.0.0_r8)THEN
          GNOD=MXNOD-MNNOD+1
          ALLOCN=1.0_r8/GNOD
          DO NE=1,NumOfPlantChemElements
            GROE(NE)=ALLOCN*GROSTKE(NE)
          ENDDO
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
          SNL=ETOL*SNL1(NZ)*(WTSTKBE(ielmc,NB,NZ)/pftPlantPopulation(NZ))**SNL2
          GROH=GROE(ielmc)/pftPlantPopulation(NZ)*SNL
          KX=MOD(MNNOD-1,JNODS1)
          IF(KX.EQ.0.AND.MNNOD-1.NE.0)KX=25
!
    !     GROWTH AT EACH CURRENT NODE
    !
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !     GRO,GROE(ielmn),GROE(ielmp)=stalk C,N,P growth at each node
    !     HTNODX,HTNODE=stalk height,stalk internode length
    !     ANGBR=sine of stalk angle from horizontal from PFT file
!
          D510: DO KK=MNNOD,MXNOD
            K1=MOD(KK,JNODS1)
            IF(K1.EQ.0.AND.KK.NE.0)K1=25
            K2=MOD(KK-1,JNODS1)
            IF(K2.EQ.0.AND.KK-1.NE.0)K2=25
            DO NE=1,NumOfPlantChemElements
              WGNODE(NE,K1,NB,NZ)=WGNODE(NE,K1,NB,NZ)+GROE(NE)
            ENDDO
            HTNODX(K1,NB,NZ)=HTNODX(K1,NB,NZ)+GROH*ANGBR(NZ)
            IF(K1.NE.0)THEN
              HTNODE(K1,NB,NZ)=HTNODX(K1,NB,NZ)+HTNODE(K2,NB,NZ)
            ELSE
              HTNODE(K1,NB,NZ)=HTNODX(K1,NB,NZ)
            ENDIF
          ENDDO D510
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
        IF(IDAY(1,NB,NZ).NE.0.AND.CEPOLB(ielmc,NB,NZ).GT.ZERO)THEN
          CCC=AZMAX1(AMIN1(1.0_r8 &
            ,safe_adb(CEPOLB(ielmn,NB,NZ),CEPOLB(ielmn,NB,NZ) &
            +CEPOLB(ielmc,NB,NZ)*CNKI) &
            ,safe_adb(CEPOLB(ielmp,NB,NZ),CEPOLB(ielmp,NB,NZ) &
            +CEPOLB(ielmc,NB,NZ)*CPKI)))
          CNC=AZMAX1(AMIN1(1.0_r8 &
            ,safe_adb(CEPOLB(ielmc,NB,NZ),CEPOLB(ielmc,NB,NZ) &
            +CEPOLB(ielmn,NB,NZ)/CNKI)))
          CPC=AZMAX1(AMIN1(1.0_r8 &
            ,safe_adb(CEPOLB(ielmc,NB,NZ),CEPOLB(ielmc,NB,NZ) &
            +CEPOLB(ielmp,NB,NZ)/CPKI)))
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
!       fTgrowCanP=temperature function for canopy growth
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
            FSNC=fTgrowCanP(NZ)*XRLA(NZ)
!
        !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
        !
        !   IFLGP=flag for remobilization
        !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !   ARLF=node leaf area
        !   RCEL(ielmc)X,RCEL(ielmn)X,RCEL(ielmp)X=remobilization of C,N,P from senescing leaf
!
            IF(IFLGP(NB,NZ).EQ.1)THEN
              DO NE=1,NumOfPlantChemElements
                WGLFEX(NE,NB,NZ)=AZMAX1(WGLFE(NE,K,NB,NZ))
              ENDDO
              ARLFZ(NB,NZ)=AZMAX1(ARLF1(K,NB,NZ))
              IF(WGLFEX(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
                RCELX(ielmc,NB,NZ)=WGLFEX(ielmc,NB,NZ)*RCCC
                RCELX(ielmn,NB,NZ)=WGLFEX(ielmn,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
                RCELX(ielmp,NB,NZ)=WGLFEX(ielmp,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
              ELSE
                RCELX(1:NumOfPlantChemElements,NB,NZ)=0._r8
              ENDIF
            ENDIF
!
    !       FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !       FSNC,FSNCL=fraction of lowest leaf to be remobilized
    !
            IF(FSNC*WGLFEX(ielmc,NB,NZ).GT.WGLFE(ielmc,K,NB,NZ) &
              .AND.WGLFEX(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
              FSNCL=AZMAX1(WGLFE(ielmc,K,NB,NZ)/WGLFEX(ielmc,NB,NZ))
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
    !       RCEL(ielmc)X,RCEL(ielmn)X,RCEL(ielmp)X=remobilization of C,N,P from senescing leaf
    !       WGLFX,WGLFNX,WGLFPX=senescing leaf C,N,P mass
    !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !       FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
            DO NE=1,NumOfPlantChemElements
              D6300: DO M=1,jsken
                ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
                  *FSNCL*(WGLFEX(NE,NB,NZ)-RCELX(NE,NB,NZ))*FWODLE(NE,k_woody_litr)
                ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,ifoliar,M,NZ) &
                  *FSNCL*(WGLFEX(NE,NB,NZ)-RCELX(NE,NB,NZ))*FWODBE(NE,k_fine_litr)
              ENDDO D6300
            ENDDO
!
    !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !       FSNCL=fraction of lowest leaf to be remobilized
    !       CanopyBranchLeafA_pft,ARLFZ=branch living,senescing leaf area
    !       WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in living,senescing leaf
    !       WSLF=leaf protein mass
    !       CNWS,CPWS=protein:N,protein:P ratios from startq.f
    !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
    !       RCEL(ielmc)X,RCEL(ielmn)X,RCEL(ielmp)X=remobilization of C,N,P from senescing leaf
!
            CanopyBranchLeafA_pft(NB,NZ)=CanopyBranchLeafA_pft(NB,NZ)-FSNCL*ARLFZ(NB,NZ)
            DO NE=1,NumOfPlantChemElements
              WTLFBE(NE,NB,NZ)=WTLFBE(NE,NB,NZ)-FSNCL*WGLFEX(NE,NB,NZ)
              WGLFE(NE,K,NB,NZ)=WGLFE(NE,K,NB,NZ)-FSNCL*WGLFEX(NE,NB,NZ)
              EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+FSNCL*RCELX(NE,NB,NZ)
            ENDDO
            ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)-FSNCL*ARLFZ(NB,NZ)
            WSLF(K,NB,NZ)=AZMAX1(WSLF(K,NB,NZ) &
              -FSNCL*AMAX1(WGLFEX(ielmn,NB,NZ)*CNWS(NZ) &
              ,WGLFEX(ielmp,NB,NZ)*CPWS(NZ)))
!
    !       REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P ALSO DEPENDS ON
    !       STRUCTURAL C:N:P
    !
    !       IFLGP=flag for remobilization
    !       WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !       CanPBranchHeight=petiole length
    !       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
!
            IF(IFLGP(NB,NZ).EQ.1)THEN
              DO NE=1,NumOfPlantChemElements
                WGSHEXE(NE,NB,NZ)=AZMAX1(WGSHE(NE,K,NB,NZ))
              ENDDO
              CanPBranchHeight(NB,NZ)=AZMAX1(CanPSheathHeight(K,NB,NZ))
              IF(WGSHEXE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
                RCESX(ielmc,NB,NZ)=RCCC*WGSHEXE(ielmc,NB,NZ)
                RCESX(ielmn,NB,NZ)=WGSHEXE(ielmn,NB,NZ) &
                  *(RCCN+(1.0_r8-RCCN)*RCESX(ielmc,NB,NZ)/WGSHEXE(ielmc,NB,NZ))
                RCESX(ielmp,NB,NZ)=WGSHEXE(ielmp,NB,NZ) &
                  *(RCCP+(1.0_r8-RCCP)*RCESX(ielmc,NB,NZ)/WGSHEXE(ielmc,NB,NZ))
              ELSE
                RCESX(1:NumOfPlantChemElements,NB,NZ)=0._r8
              ENDIF
              DO NE=1,NumOfPlantChemElements
                WTSTXBE(NE,NB,NZ)=WTSTXBE(NE,NB,NZ)+WGNODE(NE,K,NB,NZ)
              ENDDO
              WGNODE(1:NumOfPlantChemElements,K,NB,NZ)=0._r8
              HTNODX(K,NB,NZ)=0._r8
            ENDIF
!
    !       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !
            IF(FSNC*WGSHEXE(ielmc,NB,NZ).GT.WGSHE(ielmc,K,NB,NZ) &
              .AND.WGSHEXE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
              FSNCS=AZMAX1(WGSHE(ielmc,K,NB,NZ)/WGSHEXE(ielmc,NB,NZ))
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
    !       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
    !       WGSHX,WGSHNX,WGSHPX=senescing petiole C,N,P mass
    !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !       FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!
            DO NE=1,NumOfPlantChemElements
              D6305: DO M=1,jsken
                ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
                  *FSNCS*(WGSHEXE(NE,NB,NZ)-RCESX(NE,NB,NZ))*FWODBE(NE,k_woody_litr)

                ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,infoliar,M,NZ) &
                  *FSNCS*(WGSHEXE(NE,NB,NZ)-RCESX(NE,NB,NZ))*FWODBE(NE,k_fine_litr)
              ENDDO D6305
            ENDDO
!
    !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !       CanPSheathHeight,CanPBranchHeight=living,senescing petiole length
    !       WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in living,senescing petiole
    !       WSSHE=petiole protein mass
    !       CNWS,CPWS=protein:N,protein:P ratios from startq.f
    !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
!
            DO NE=1,NumOfPlantChemElements
              WTSHEBE(NE,NB,NZ)=WTSHEBE(NE,NB,NZ)-FSNCS*WGSHEXE(NE,NB,NZ)
              WGSHE(NE,K,NB,NZ)=WGSHE(NE,K,NB,NZ)-FSNCS*WGSHEXE(NE,NB,NZ)
              EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+FSNCS*RCESX(NE,NB,NZ)
            ENDDO
            CanPSheathHeight(K,NB,NZ)=CanPSheathHeight(K,NB,NZ)-FSNCS*CanPBranchHeight(NB,NZ)

            WSSHE(K,NB,NZ)=AZMAX1(WSSHE(K,NB,NZ) &
              -FSNCS*AMAX1(WGSHEXE(ielmn,NB,NZ)*CNWS(NZ) &
              ,WGSHEXE(ielmp,NB,NZ)*CPWS(NZ)))

          ENDIF
        ENDIF
!
  !     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0
  !
  !     SNCR=excess maintenance respiration
  !     WTRSVB=stalk reserve C mass
  !     RCO2V=remobilization of stalk reserve C
  !     VMXC=rate constant for nonstructural C oxidation in respiration
  !     fTgrowCanP=temperature function for canopy growth
  !
        IF(IFLGZ.EQ.0)THEN
          IF(SNCR.GT.0.0.AND.WTRSVBE(ielmc,NB,NZ).GT.0.0)THEN
            RCO2V=AMIN1(SNCR,VMXC*WTRSVBE(ielmc,NB,NZ)*fTgrowCanP(NZ))
            WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)-RCO2V
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
!       CanPBLeafShethC=leaf+petiole mass
!       FLGZ=control rate of remobilization
!       FLGZX=number of hours until full senescence after physl maturity
!       SNCX=total senescence respiration
!       KVSTG,KVSTGN=integer of highest,lowest leaf number currently growing
!       KSNC=number of nodes undergoing remobilization
!       SNCF=ratio of phenologically-driven vs total senescence respiration
!
        IF(IFLGZ.EQ.1.AND.IFLGY.EQ.1.AND.ISTYP(NZ).NE.iplt_annual)THEN
          SNCZ=FXFB(IBTYP(NZ))*CanPBLeafShethC(NB,NZ)*AMIN1(1.0_r8,FLGZ(NB,NZ)/FLGZX)
        ELSE
          SNCZ=0._r8
        ENDIF
        SNCX=SNCR+SNCZ
        IF(SNCX.GT.ZEROP(NZ))THEN
          SNCF=SNCZ/SNCX
          KSNC=INT(0.5_r8*(KVSTG(NB,NZ)-KVSTGN(NB,NZ)))+1
          XKSNC=KSNC
          KN=MAX(0,KVSTGN(NB,NZ)-1)
    !
    !     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
    !     IF MAIN STEM POOLS ARE DEPLETED
    !
    !     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
    !     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !     XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
    !
          IF(IBTYP(NZ).NE.0.AND.IGTYP(NZ).GT.1 &
            .AND.NB.EQ.NB1(NZ).AND.isclose(SNCF,0._r8))THEN
            NBY=0
          D584: DO NBL=1,NumOfBranches_pft(NZ)
            NBZ(NBL)=0
          ENDDO D584
          D586: DO NBL=1,NumOfBranches_pft(NZ)
            NBX=KVSTG(NB,NZ)
            D585: DO NBK=1,NumOfBranches_pft(NZ)
              IF(IDTHB(NBK,NZ).EQ.0.AND.NBK.NE.NB1(NZ) &
                .AND.BranchNumber_brchpft(NBK,NZ).LT.NBX.AND.BranchNumber_brchpft(NBK,NZ).GT.NBY)THEN
                NBZ(NBL)=NBK
                NBX=BranchNumber_brchpft(NBK,NZ)
              ENDIF
            ENDDO D585
            IF(NBZ(NBL).NE.0)THEN
              NBY=BranchNumber_brchpft(NBZ(NBL),NZ)
            ENDIF
          ENDDO D586
          D580: DO NBL=1,NumOfBranches_pft(NZ)
            IF(NBZ(NBL).NE.0)THEN
              IF(BranchNumber_brchpft(NBZ(NBL),NZ).LT.KK)THEN
                IF(EPOOL(ielmc,NBZ(NBL),NZ).GT.0)THEN
                  XFRE(ielmc)=1.0E-02_r8*AMIN1(SNCX,EPOOL(ielmc,NBZ(NBL),NZ))
                  DO NE=2,NumOfPlantChemElements
                    XFRE(NE)=XFRE(NE)*EPOOL(NE,NBZ(NBL),NZ)/EPOOL(ielmc,NBZ(NBL),NZ)
                  ENDDO
                ELSE
                  XFRE(ielmc)=0._r8
                  DO NE=2,NumOfPlantChemElements
                    XFRE(NE)=1.0E-02_r8*EPOOL(NE,NBZ(NBL),NZ)
                  ENDDO
                ENDIF
                DO NE=1,NumOfPlantChemElements
                  EPOOL(NE,NBZ(NBL),NZ)=EPOOL(NE,NBZ(NBL),NZ)-XFRE(NE)
                ENDDO
                EPOOL(ielmc,NB1(NZ),NZ)=EPOOL(ielmc,NB1(NZ),NZ)+XFRE(ielmc)*SNCF
                DO NE=2,NumOfPlantChemElements
                  EPOOL(NE,NB1(NZ),NZ)=EPOOL(NE,NB1(NZ),NZ)+XFRE(NE)
                ENDDO
                SNCX=SNCX-XFRE(ielmc)
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
    IF(IBTYP(NZ).NE.0.AND.IGTYP(NZ).GT.1.AND.IDTHB(NB1(NZ),NZ).EQ.1)then
      IDTHB(NB,NZ)=ibrdead
    endif
!
!     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
!
    KVSTGX=MAX(0,KVSTG(NB,NZ)-24)
    D495: DO KK=KVSTGX,KVSTG(NB,NZ)
      K=MOD(KK,JNODS1)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      IF(WGLFE(ielmc,K,NB,NZ).GT.0.0)THEN
        CPOOLT=WGLFE(ielmc,K,NB,NZ)+EPOOL(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZEROP(NZ))THEN
          ZPOOLD=WGLFE(ielmn,K,NB,NZ)*EPOOL(ielmc,NB,NZ) &
            -EPOOL(ielmn,NB,NZ)*WGLFE(ielmc,K,NB,NZ)
          XFRN1=AZMAX1(AMIN1(1.0E-03_r8*ZPOOLD/CPOOLT,WGLFE(ielmn,K,NB,NZ) &
            -ZPLFM*CNLFB*WGLFE(ielmc,K,NB,NZ)))
          PPOOLD=WGLFE(ielmp,K,NB,NZ)*EPOOL(ielmc,NB,NZ) &
            -EPOOL(ielmp,NB,NZ)*WGLFE(ielmc,K,NB,NZ)
          XFRP1=AZMAX1(AMIN1(1.0E-03*PPOOLD/CPOOLT,WGLFE(ielmp,K,NB,NZ) &
            -ZPLFM*CPLFB*WGLFE(ielmc,K,NB,NZ)))
          XFRE(ielmn)=AMAX1(XFRN1,10.0_r8*XFRP1)
          XFRE(ielmp)=AMAX1(XFRP1,0.10_r8*XFRN1)
          DO NE=2,NumOfPlantChemElements
            WGLFE(NE,K,NB,NZ)=WGLFE(NE,K,NB,NZ)-XFRE(NE)
            WTLFBE(NE,NB,NZ)=WTLFBE(NE,NB,NZ)-XFRE(NE)
            EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+XFRE(NE)
          ENDDO
          WSLF(K,NB,NZ)=AZMAX1(WSLF(K,NB,NZ)-AMAX1(XFRE(ielmn)*CNWS(NZ),XFRE(ielmp)*CPWS(NZ)))
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
    !     LeafA_lyrnodbrchpft=leaf node surface area in canopy layer
    !     ARLF,CanPLNBLA=leaf node surface area in canopy layer
    !     ZC,DPTHS=canopy,snowpack height
    !     CLASS=leaf inclination class
    !
    IF(SSINN.GT.0.0)THEN
      call LeafClassAllocation(NB,NZ)
    ENDIF

    call GrainFilling(I,NB,NZ,GROGRE,GROSTKE(ielmc))
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
    PSICanP   =>  plt_ew%PSICanP        , &
    CanPBStalkC  =>  plt_biom%CanPBStalkC     , &
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
    CO2NetFix_pft    =>  plt_bgcr%CO2NetFix_pft       , &
    SNL1    =>  plt_morph%SNL1      , &
    CanopyLeafA_pft   =>  plt_morph%CanopyLeafA_pft     , &
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
  D10: DO N=1,JPRT
    PART(N)=0._r8
  ENDDO D10
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
      IF(ISTYP(NZ).EQ.iplt_annual)THEN
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
    IF(ISTYP(NZ).EQ.iplt_annual)THEN
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
!     WTRSVB,CanPBStalkC=stalk reserve,sapwood mass
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!
  IF(IDAY(2,NB,NZ).NE.0)THEN
    IF(WTRSVBE(ielmc,NB,NZ).LT.XFRX*CanPBStalkC(NB,NZ))THEN
      D1020: DO N=1,JPRT
        IF(N.NE.4)THEN
          PART(4)=PART(4)+0.10*PART(N)
          PART(N)=PART(N)-0.10*PART(N)
        ENDIF
      ENDDO D1020
!
!     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
!
    ELSEIF(WTRSVBE(ielmc,NB,NZ).GT.1.0*CanPBStalkC(NB,NZ))THEN
      PART(3)=PART(3)+PART(4)+PART(7)
      PART(4)=0._r8
      PART(7)=0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
!
!     CanopyLeafA_pft=PFT leaf area
!
  ARLFI=CanopyLeafA_pft(NZ)/AREA3(NU)
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
  IF((ISTYP(NZ).NE.iplt_annual.AND.VRNF(NB,NZ).GE.FVRN(IWTYP(NZ))*VRNX(NB,NZ)) &
    .OR.(ISTYP(NZ).EQ.iplt_annual.AND.IDAY(8,NB,NZ).NE.0))THEN
    IFLGZ=1
    IF(ISTYP(NZ).EQ.iplt_annual.OR.IWTYP(NZ).EQ.0)THEN
      IFLGY=1
      FLGZ(NB,NZ)=FLGZ(NB,NZ)+1.0_r8
    ELSEIF((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3).AND.TCC(NZ).LT.CTC(NZ))THEN
      IFLGY=1
      FLGZ(NB,NZ)=FLGZ(NB,NZ)+1.0_r8
    ELSEIF(IWTYP(NZ).GE.2.AND.PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
      IFLGY=1
      FLGZ(NB,NZ)=FLGZ(NB,NZ)+1.0_r8
    ENDIF
    IF(ISTYP(NZ).NE.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
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
  D1000: DO N=1,JPRT
    IF(N.EQ.3.AND.isclose(SNL1(NZ),0._r8))THEN
      PART(N)=0._r8
    ELSE
      PART(N)=AZMAX1(PART(N))
    ENDIF
    TOTAL=TOTAL+PART(N)
  ENDDO D1000
  IF(TOTAL.GT.ZERO)THEN
    D1010: DO N=1,JPRT
      PART(N)=PART(N)/TOTAL
    ENDDO D1010
  ELSE
    D1015: DO N=1,JPRT
      PART(N)=0._r8
    ENDDO D1015
  ENDIF
  end associate
  end subroutine CalcPartitionCoeff

!------------------------------------------------------------------------------------------

  subroutine UpdatePhotosynthates(I,J,NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
    WTSHXN,TFN5,WFNG,Stomata_Activity,WFNSG,CH2O3,CH2O4,CNPG,rco2c,RMNCS,SNCR,CGROS,CNRDM,CNRDA)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFN6(JZ1)
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,TFN5,WFNG
  real(r8), intent(in) :: Stomata_Activity,WFNSG
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

    call ComputeGPP(NB,NZ,WFNG,Stomata_Activity,CH2O3,CH2O4,CH2O,CO2F)
!
!   SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
!
    call ComputeRAutoAfEmergence(NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
      CO2F,CH2O,TFN5,WFNG,WFNSG,WTSHXN,ZADDB,CNPG,PADDB,RCO2C,RMNCS,SNCR,CGROS,CNRDA)

!   SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
!
  ELSE
    call ComputeRAutoB4Emergence(NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,WTSHXN,&
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
!   XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=branch-root layer C,N,P transfer
!
  EPOOL(ielmc,NB,NZ)=EPOOL(ielmc,NB,NZ)+CH2O-AMIN1(RMNCS,RCO2C)-CGROS-CNRDA
  EPOOL(ielmn,NB,NZ)=EPOOL(ielmn,NB,NZ)-ZADDB+RNH3B(NB,NZ)
  EPOOL(ielmp,NB,NZ)=EPOOL(ielmp,NB,NZ)-PADDB
  end associate
  end subroutine UpdatePhotosynthates
!------------------------------------------------------------------------------------------

  subroutine C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: CH2O3(25),CH2O4(25)
  integer :: K
  real(r8) :: CPL4M,CCBS,CPL3K,CO2LeakFromBundsheth

!     begin_execution
  associate(                              &
    CO2NetFix_pft    =>  plt_bgcr%CO2NetFix_pft       , &
    TCO2T   =>  plt_bgcr%TCO2T      , &
    TCO2A   =>  plt_bgcr%TCO2A      , &
    Eco_AutoR_col    =>  plt_bgcr%Eco_AutoR_col       , &
    ECO_ER_col    =>  plt_bgcr%ECO_ER_col       , &
    WGLFE   =>  plt_biom%WGLFE      , &
    ZEROP   =>  plt_biom%ZEROP      , &
    HCOB    =>  plt_photo%HCOB      , &
    CO2B    =>  plt_photo%CO2B      , &
    CPOOL3  =>  plt_photo%CPOOL3    , &
    CPOOL4  =>  plt_photo%CPOOL4    , &
    CO2L    =>  plt_photo%CO2L        &
  )
  D170: DO K=1,JNODS1
    IF(WGLFE(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
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
      CPL4M=1.0_r8*(CPOOL4(K,NB,NZ)*WGLFE(ielmc,K,NB,NZ)*FBS &
        -CPOOL3(K,NB,NZ)*WGLFE(ielmc,K,NB,NZ)*FMP) &
        /(WGLFE(ielmc,K,NB,NZ)*(FBS+FMP))
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
      CCBS=AZMAX1(0.083E+09*CO2B(K,NB,NZ)/(WGLFE(ielmc,K,NB,NZ)*FBS))
      CPL3K=2.5E-02*CPOOL3(K,NB,NZ)/(1.0+CCBS/CO2KI)
      CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)-CPL3K
      CO2B(K,NB,NZ)=CO2B(K,NB,NZ)+FCO2B*CPL3K
      HCOB(K,NB,NZ)=HCOB(K,NB,NZ)+FHCOB*CPL3K
!
!     BUNDLE SHEATH LEAKAGE
!
!     CO2LeakFromBundsheth=bundle sheath CO2 leakage
!     CCBS=CO2 concn in bundle sheath (uM)
!     CO2L=intercellular CO2 concentration (uM)
!     WGLF=node leaf C mass
!     FBS=leaf water content in bundle sheath
!     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
      CO2LeakFromBundsheth=5.0E-07_r8*(CCBS-CO2L(NZ))*WGLFE(ielmc,K,NB,NZ)*FBS
      CO2B(K,NB,NZ)=CO2B(K,NB,NZ)-FCO2B*CO2LeakFromBundsheth
      HCOB(K,NB,NZ)=HCOB(K,NB,NZ)-FHCOB*CO2LeakFromBundsheth

!
!     TOTAL C EXCHANGE
!
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CO2NetFix_pft=PFT net CO2 fixation
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!     CO2LeakFromBundsheth=bundle sheath CO2 leakage
!
      TCO2T(NZ)=TCO2T(NZ)-CO2LeakFromBundsheth
      TCO2A(NZ)=TCO2A(NZ)-CO2LeakFromBundsheth
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-CO2LeakFromBundsheth
      ECO_ER_col=ECO_ER_col-CO2LeakFromBundsheth
      Eco_AutoR_col=Eco_AutoR_col-CO2LeakFromBundsheth
    ENDIF
  ENDDO D170
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
  integer :: N,M,K,KK,MXNOD,MNNOD,NE
  real(r8) :: FSNCL,FSNCS
  real(r8) :: FNCLF,FSNAL,FSNAS,FRCC
  real(r8) :: FSNCK
  real(r8) :: FSNCR
  real(r8) :: HTNODZ
  real(r8) :: RCEL(1:NumOfPlantChemElements)
  real(r8) :: RCES(1:NumOfPlantChemElements)

  real(r8) :: RCSC,RCSN,RCSP
  real(r8) :: RCEK(NumOfPlantChemElements)
  real(r8) :: SNCR
  real(r8) :: SNCZ
  real(r8) :: SNCT
  real(r8) :: SNCLF
  real(r8) :: SNCSH
! begin_execution
  associate(                             &
    instruct   =>  pltpar%instruct     , &
    ifoliar    =>  pltpar%ifoliar      , &
    istalk     =>  pltpar%istalk       , &
    iroot      =>  pltpar%iroot        , &
    infoliar   =>  pltpar%infoliar     , &
    icwood     =>  pltpar%icwood       , &
    k_fine_litr=>  pltpar%k_fine_litr  , &
    k_woody_litr=> pltpar%k_woody_litr , &
    WTRSVBE    =>  plt_biom%WTRSVBE    , &
    WTSTKBE    =>  plt_biom%WTSTKBE    , &
    WGNODE     =>  plt_biom%WGNODE     , &
    CanPBStalkC     =>  plt_biom%CanPBStalkC     , &
    WGSHE      =>  plt_biom%WGSHE      , &
    WGLFE      =>  plt_biom%WGLFE      , &
    WSLF       =>  plt_biom%WSLF       , &
    WTRVE      =>  plt_biom%WTRVE      , &
    WTLFBE     =>  plt_biom%WTLFBE     , &
    WTSHEBE    =>  plt_biom%WTSHEBE    , &
    WTSTXBE    =>  plt_biom%WTSTXBE    , &
    WSSHE      =>  plt_biom%WSSHE      , &
    ZEROP      =>  plt_biom%ZEROP      , &
    ZEROL      =>  plt_biom%ZEROL      , &
    EPOOL      =>  plt_biom%EPOOL      , &
    pftPlantPopulation         =>   plt_site%pftPlantPopulation        , &
    CPOOL3     =>  plt_photo%CPOOL3    , &
    CPOOL4     =>  plt_photo%CPOOL4    , &
    FWOODE     =>  plt_allom%FWOODE    , &
    FWODBE     =>  plt_allom%FWODBE    , &
    FWODLE     =>  plt_allom%FWODLE    , &
    CNWS       =>  plt_allom%CNWS      , &
    CPWS       =>  plt_allom%CPWS      , &
    ESNC       =>  plt_bgcr%ESNC       , &
    CFOPE      =>  plt_soilchem%CFOPE  , &
    icarbhyro  =>  pltpar%icarbhyro    , &
    KVSTG      =>   plt_pheno%KVSTG    , &
    IDTHB      =>   plt_pheno%IDTHB    , &
    ISTYP      =>   plt_pheno%ISTYP    , &
    CanopyBranchLeafA_pft      =>   plt_morph%CanopyBranchLeafA_pft    , &
    NNOD       =>   plt_morph%NNOD     , &
    HTNODE     =>   plt_morph%HTNODE   , &
    CanPSheathHeight      =>   plt_morph%CanPSheathHeight    , &
    ARLF1      =>   plt_morph%ARLF1    , &
    HTNODX     =>   plt_morph%HTNODX     &
  )
!     REMOBILIZATION AND LITTERFALL WHEN GROWTH RESPIRATION < 0
!     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
!
!     SNCX,SNCT=branch,node senescence respiration
!     KSNC=number of nodes undergoing remobilization
!

  D575: DO N=1,KSNC
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
    !       RCEL(ielmc),RCEL(ielmn),RCEL(ielmp)=remobilization of C,N,P from senescing leaf
    !       RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !       RCCZ,RCCY=min,max fractions for shoot C recycling
    !
      IF(WGLFE(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
        FNCLF=WGLFE(ielmc,K,NB,NZ)/(WGLFE(ielmc,K,NB,NZ)+WGSHE(ielmc,K,NB,NZ))
        SNCLF=FNCLF*SNCT
        SNCSH=SNCT-SNCLF
        RCEL(ielmc)=RCCC*WGLFE(ielmc,K,NB,NZ)
        RCEL(ielmn)=WGLFE(ielmn,K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCEL(ielmp)=WGLFE(ielmp,K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
    !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !         FSNCL,FSNAL=fraction of current leaf C,area to be remobilized
    !
        IF(RCEL(ielmc).GT.ZEROP(NZ))THEN
          FSNCL=AZMAX1(AMIN1(1.0,SNCLF/RCEL(ielmc)))
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
        !     RCEL(ielmc),RCEL(ielmn),RCEL(ielmp)=remobilization of C,N,P from senescing leaf
        !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
        !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
        !     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
        !
        DO NE=1,NumOfPlantChemElements
          D6310: DO M=1,jsken
            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *FSNCL*(WGLFE(NE,K,NB,NZ)-RCEL(NE))*FWODLE(NE,k_woody_litr)

            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,ifoliar,M,NZ) &
              *FSNCL*(WGLFE(NE,K,NB,NZ)-RCEL(NE))*FWODLE(NE,k_fine_litr)
          ENDDO D6310
        ENDDO
        IF(K.NE.0)THEN
          ESNC(ielmc,icarbhyro,k_fine_litr,0,NZ)=ESNC(ielmc,icarbhyro,k_fine_litr,0,NZ) &
            +FSNCL*(CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ))
          CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)-FSNCL*CPOOL3(K,NB,NZ)
          CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)-FSNCL*CPOOL4(K,NB,NZ)
        ENDIF
!
!     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
!
!     CanopyBranchLeafA_pft=leaf area
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     FSNCL=fraction of current leaf to be remobilized
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
        CanopyBranchLeafA_pft(NB,NZ)=AZMAX1(CanopyBranchLeafA_pft(NB,NZ)-FSNAL*ARLF1(K,NB,NZ))
        ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)-FSNAL*ARLF1(K,NB,NZ)

        DO NE=1,NumOfPlantChemElements
          WTLFBE(NE,NB,NZ)=AZMAX1(WTLFBE(NE,NB,NZ)-FSNCL*WGLFE(NE,K,NB,NZ))
          WGLFE(NE,K,NB,NZ)=WGLFE(NE,K,NB,NZ)-FSNCL*WGLFE(NE,K,NB,NZ)
        ENDDO
        WSLF(K,NB,NZ)=AZMAX1(WSLF(K,NB,NZ) &
          -FSNCL*AMAX1(WGLFE(ielmn,K,NB,NZ)*CNWS(NZ) &
          ,WGLFE(ielmp,K,NB,NZ)*CPWS(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
!     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     FSNCL=fraction of current leaf C to be remobilized
!     RCEL(ielmc),RCEL(ielmn),RCEL(ielmp)=remobilization of C,N,P from senescing leaf
!     SNCLF,SNCT=remaining senescence respiration carried to next node
!
        EPOOL(ielmc,NB,NZ)=EPOOL(ielmc,NB,NZ)+FSNCL*RCEL(ielmc)*SNCF
        DO NE=2,NumOfPlantChemElements
          EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+FSNCL*RCEL(NE)
        ENDDO
        SNCLF=SNCLF-FSNCL*RCEL(ielmc)
        SNCT=SNCT-FSNCL*RCEL(ielmc)
        IF(WTLFBE(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          WTLFBE(ielmc,NB,NZ)=0._r8
          CanopyBranchLeafA_pft(NB,NZ)=0._r8
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
        !     CanopyBranchLeafA_pft=leaf area
        !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
        !     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
        !
      ELSE
        DO NE=1,NumOfPlantChemElements
          D6315: DO M=1,jsken
            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *WGLFE(NE,K,NB,NZ)*FWODLE(NE,k_woody_litr)

            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,ifoliar,M,NZ) &
              *WGLFE(NE,K,NB,NZ)*FWODLE(NE,k_fine_litr)
          ENDDO D6315
        ENDDO
        IF(K.NE.0)THEN
          ESNC(ielmc,icarbhyro,k_fine_litr,0,NZ)=ESNC(ielmc,icarbhyro,k_fine_litr,0,NZ)+CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ)
          CPOOL3(K,NB,NZ)=0._r8
          CPOOL4(K,NB,NZ)=0._r8
        ENDIF
        CanopyBranchLeafA_pft(NB,NZ)=AZMAX1(CanopyBranchLeafA_pft(NB,NZ)-ARLF1(K,NB,NZ))
        DO NE=1,NumOfPlantChemElements
          WTLFBE(NE,NB,NZ)=AZMAX1(WTLFBE(NE,NB,NZ)-WGLFE(NE,K,NB,NZ))
          WGLFE(NE,K,NB,NZ)=0._r8
        ENDDO
        ARLF1(K,NB,NZ)=0._r8
        WSLF(K,NB,NZ)=0._r8
        IF(WTLFBE(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          WTLFBE(ielmc,NB,NZ)=0._r8
          CanopyBranchLeafA_pft(NB,NZ)=0._r8
        ENDIF
      ENDIF
!564   CONTINUE
!
    !     REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P DEPENDS ON
    !     NON-STRUCTURAL C:N:P
    !
    !     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !
      IF(WGSHE(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
        RCES(ielmc)=WGSHE(ielmc,K,NB,NZ)*RCCC
        RCES(ielmn)=WGSHE(ielmn,K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCES(ielmp)=WGSHE(ielmp,K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
      !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
      !     OR PETIOLE
      !
      !     FSNCS,FSNAS=fraction of current petiole C,length to be remobilized
      !
        IF(RCES(ielmc).GT.ZEROP(NZ))THEN
          FSNCS=AZMAX1(AMIN1(1.0,SNCSH/RCES(ielmc)))
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
      !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
      !
        DO NE=1,NumOfPlantChemElements
          D6320: DO M=1,jsken
            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *FSNCS*(WGSHE(NE,K,NB,NZ)-RCES(NE))*FWODBE(NE,k_woody_litr)
            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,infoliar,M,NZ) &
              *FSNCS*(WGSHE(NE,K,NB,NZ)-RCES(NE))*FWODBE(NE,k_fine_litr)
          ENDDO D6320
          WTSHEBE(NE,NB,NZ)=AZMAX1(WTSHEBE(NE,NB,NZ)-FSNCS*WGSHE(NE,K,NB,NZ))
          WGSHE(NE,K,NB,NZ)=WGSHE(NE,K,NB,NZ)-FSNCS*WGSHE(NE,K,NB,NZ)
        ENDDO
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     CanPSheathHeight=petiole length
      !     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !     FSNCS=fraction of current petiole to be remobilized
      !     CNWS,CPWS=protein:N,protein:P ratios from startq.f
      !

        CanPSheathHeight(K,NB,NZ)=CanPSheathHeight(K,NB,NZ)-FSNAS*CanPSheathHeight(K,NB,NZ)
        WSSHE(K,NB,NZ)=AZMAX1(WSSHE(K,NB,NZ) &
          -FSNCS*AMAX1(WGSHE(ielmn,K,NB,NZ)*CNWS(NZ) &
          ,WGSHE(ielmp,K,NB,NZ)*CPWS(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
  !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
  !
  !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  !     FSNCS=fraction of current petiole C to be remobilized
  !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
  !     SNCSH,SNCT=remaining senescence respiration carried to next node
  !
        EPOOL(ielmc,NB,NZ)=EPOOL(ielmc,NB,NZ)+FSNCS*RCES(ielmc)*SNCF
        DO NE=2,NumOfPlantChemElements
          EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+FSNCS*RCES(NE)
        ENDDO
        SNCSH=SNCSH-FSNCS*RCES(ielmc)
        SNCT=SNCT-FSNCS*RCES(ielmc)
        IF(WTSHEBE(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          WTSHEBE(ielmc,NB,NZ)=0._r8
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
      !     CanPSheathHeight=petiole length
      !     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !
      ELSE
        DO NE=1,NumOfPlantChemElements
          D6325: DO M=1,jsken
            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *WGSHE(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,infoliar,M,NZ) &
              *WGSHE(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)
          ENDDO D6325
          WTSHEBE(NE,NB,NZ)=AZMAX1(WTSHEBE(NE,NB,NZ)-WGSHE(NE,K,NB,NZ))
          WGSHE(NE,K,NB,NZ)=0._r8
        ENDDO
        CanPSheathHeight(K,NB,NZ)=0._r8
        WSSHE(K,NB,NZ)=0._r8
        IF(WTSHEBE(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          WTSHEBE(ielmc,NB,NZ)=0._r8
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
    IF(WTRSVBE(ielmc,NB,NZ).GT.SNCR)THEN
      WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)-SNCR
      SNCR=0._r8
      cycle
    ENDIF
!
!     REMOBILIZATION OF STALK C,N,P
!
!     FXFS=rate constant for remobilization of stalk C,N,P (h-1)
!     SNCZ=phenologically-driven respiration senescence during late-season
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     WTSTKB,CanPBStalkC=stalk,sapwood C mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     MXNOD,MNNOD=max,min node number currently growing
!     NNOD=number of concurrently growing nodes
!     KVSTG=integer of most recent leaf number
!
    SNCZ=FXFS*SNCR
    SNCT=SNCR+SNCZ
    IF(ISTYP(NZ).NE.iplt_annual.AND.SNCT.GT.ZEROP(NZ) &
      .AND.WTSTKBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      SNCF=SNCZ/SNCT
      FRCC=CanPBStalkC(NB,NZ)/WTSTKBE(ielmc,NB,NZ)
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
    !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !
        IF(WGNODE(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
          RCEK(ielmc)=WGNODE(ielmc,K,NB,NZ)*RCSC
          RCEK(ielmn)=WGNODE(ielmn,K,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
          RCEK(ielmp)=WGNODE(ielmp,K,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
      !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
      !
      !     FSNCS=fraction of lowest internode to be remobilized
!
          IF(RCEK(ielmc).GT.ZEROP(NZ))THEN
            FSNCK=AZMAX1(AMIN1(1.0,SNCT/RCEK(ielmc)))
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
      !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
      !     WGNODE,WGNODN,WGNODP=senescing internode C,N,P mass
!
          DO NE=1,NumOfPlantChemElements
            D7310: DO M=1,jsken
              ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
                *FSNCK*(WGNODE(NE,K,NB,NZ)-RCEK(NE))*FWOODE(NE,k_woody_litr)

              ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
                *FSNCK*(WGNODE(NE,K,NB,NZ)-RCEK(NE))*FWOODE(NE,k_fine_litr)
            ENDDO D7310

!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     FSNCK=fraction of lowest internode to be remobilized
      !     HTNODE,HTNODX=living,senescing internode length
      !     WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
      !
            WTSTKBE(NE,NB,NZ)=AZMAX1(WTSTKBE(NE,NB,NZ)-FSNCK*WGNODE(NE,K,NB,NZ))
            WGNODE(NE,K,NB,NZ)=WGNODE(NE,K,NB,NZ)-FSNCK*WGNODE(NE,K,NB,NZ)
          ENDDO
          HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)-FSNCK*HTNODX(K,NB,NZ)
          HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ)-FSNCK*HTNODX(K,NB,NZ)
!
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
      !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
      !     FSNCK=fraction of lowest internode to be remobilized
      !     SNCT=remaining node senescence respiration
      !
          WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)+FSNCK*RCEK(ielmc)*SNCF
          DO NE=2,NumOfPlantChemElements
            WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)+FSNCK*RCEK(NE)
          ENDDO
          SNCT=SNCT-FSNCK*RCEK(ielmc)
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
          DO NE=1,NumOfPlantChemElements
            D7315: DO M=1,jsken
              ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
                *WGNODE(NE,K,NB,NZ)*FWOODE(NE,k_woody_litr)

              ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
                *WGNODE(NE,K,NB,NZ)*FWOODE(NE,k_fine_litr)
            ENDDO D7315
            WTSTKBE(NE,NB,NZ)=AZMAX1(WTSTKBE(NE,NB,NZ)-WGNODE(NE,K,NB,NZ))
            WGNODE(NE,K,NB,NZ)=0._r8
          ENDDO
          HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)-HTNODX(K,NB,NZ)
          HTNODX(K,NB,NZ)=0._r8
        ENDIF
      ENDDO D1650
!
    !   RESIDUAL STALK
    !
    !   RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
      IF(WTSTXBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
        RCEK(ielmc)=WTSTXBE(ielmc,NB,NZ)*RCSC
        RCEK(ielmn)=WTSTXBE(ielmn,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
        RCEK(ielmp)=WTSTXBE(ielmp,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
    !     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
        IF(RCEK(ielmc).GT.ZEROP(NZ))THEN
          FSNCR=AZMAX1(AMIN1(1.0,SNCT/RCEK(ielmc)))
        ELSE
          FSNCR=1.0_r8
        ENDIF
    !
    !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !     TO FRACTIONS SET IN 'STARTQ'
    !
        DO NE=1,NumOfPlantChemElements
          D8310: DO M=1,jsken
            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
              *FSNCR*(WTSTXBE(NE,NB,NZ)-RCEK(NE))*FWOODE(NE,k_woody_litr)

            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
              *FSNCR*(WTSTXBE(NE,NB,NZ)-RCEK(NE))*FWOODE(NE,k_fine_litr)
          ENDDO D8310

    !
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !     HTNODE,HTNODX=living,senescing internode length
    !
          WTSTKBE(NE,NB,NZ)=AZMAX1(WTSTKBE(NE,NB,NZ)-FSNCR*WTSTXBE(NE,NB,NZ))

          WTSTXBE(NE,NB,NZ)=AZMAX1(WTSTXBE(NE,NB,NZ)-FSNCR*WTSTXBE(NE,NB,NZ))
        ENDDO
        HTNODZ=0._r8
        D8320: DO K=0,JNODS1
          HTNODZ=AMAX1(HTNODZ,HTNODE(K,NB,NZ))
        ENDDO D8320
        HTNODZ=AZMAX1(HTNODZ-FSNCR*HTNODZ)
        D8325: DO K=0,JNODS1
          HTNODE(K,NB,NZ)=AMIN1(HTNODZ,HTNODE(K,NB,NZ))
        ENDDO D8325
!
    !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
    !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
    !
    !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
    !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !     FSNCR=fraction of residual stalk to be remobilized
    !     SNCT=remaining node senescence respiration
    !
        WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)+FSNCR*RCEK(ielmc)*SNCF
        DO NE=2,NumOfPlantChemElements
          WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)+FSNCR*RCEK(NE)
        ENDDO
        SNCT=SNCT-FSNCR*RCEK(ielmc)
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
      DO NE=1,NumOfPlantChemElements
        D8315: DO M=1,jsken
          ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
            *WTSTXBE(NE,NB,NZ)*FWOODE(NE,k_woody_litr)

          ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
            *WTSTXBE(NE,NB,NZ)*FWOODE(NE,k_fine_litr)
        ENDDO D8315
        WTSTKBE(NE,NB,NZ)=AZMAX1(WTSTKBE(NE,NB,NZ)-WTSTXBE(NE,NB,NZ))
        WTSTXBE(NE,NB,NZ)=0._r8
      ENDDO
    ENDIF
!
!     REMOBILIZATION OF STORAGE C,N,P
!
!     WTRVC=storage C
!     IDTHB=branch living flag: 0=alive,1=dead
!     SNCR=remaining excess maintenance respiration
!
    SNCR=SNCT/(1.0_r8+FXFS)
    IF(WTRVE(ielmc,NZ).GT.SNCR)THEN
      WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)-SNCR
      SNCR=0._r8
    ELSEIF(ISTYP(NZ).NE.iplt_annual)THEN
      IDTHB(NB,NZ)=ibrdead
    ENDIF
565 CONTINUE
  ENDDO D575
  end associate
  end subroutine RemobilizeLeafLayers


!------------------------------------------------------------------------------------------

  subroutine AllocateLeafToCanopyLayers(NB,NZ,ZCX)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: ZCX(JP1)
  integer  :: LL,LU,L,K,k1,k2,KK,NE
  integer  :: KVSTGX,KVSTG1,LHTLFU,LHTLFL
  integer  :: LHTBRU,LHTBRL,N
  real(r8) :: ZSTK
  real(r8) :: YWGLFE(NumOfPlantChemElements)
  real(r8) :: YARLF,YLGLF,XLGLF
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
    WGLFLE   =>  plt_biom%WGLFLE  , &
    WGLFE    =>  plt_biom%WGLFE   , &
    CanopyLeafCpft_lyr    =>  plt_biom%CanopyLeafCpft_lyr   , &
    CanPBStalkC   =>  plt_biom%CanPBStalkC  , &
    WTSTKBE  =>  plt_biom%WTSTKBE , &
    WSSHE    =>  plt_biom%WSSHE   , &
    FNOD     =>  plt_allom%FNOD   , &
    IGTYP    =>  plt_pheno%IGTYP  , &
    ISTYP    =>  plt_pheno%ISTYP  , &
    KVSTG    =>  plt_pheno%KVSTG  , &
    IBTYP    =>  plt_pheno%IBTYP  , &
    KVSTGN   =>  plt_pheno%KVSTGN , &
    WDLF     => plt_morph%WDLF    , &
    SeedinDepth    => plt_morph%SeedinDepth   , &
    NB1      => plt_morph%NB1     , &
    CanPLNBLA    => plt_morph%CanPLNBLA   , &
    CanPSheathHeight    => plt_morph%CanPSheathHeight   , &
    CanopyBranchStemApft_lyr    => plt_morph%CanopyBranchStemApft_lyr   , &
    BranchNumber_brchpft     => plt_morph%BranchNumber_brchpft    , &
    HTNODE   => plt_morph%HTNODE  , &
    CanopyHeightz       => plt_morph%CanopyHeightz      , &
    CanopyHeight       => plt_morph%CanopyHeight      , &
    HypoctoylHeight    => plt_morph%HypoctoylHeight   , &
    CLASS    => plt_morph%CLASS   , &
    ARLF1    => plt_morph%ARLF1   , &
    CanopyLeafApft_lyr    => plt_morph%CanopyLeafApft_lyr   , &
    pftPlantPopulation       => plt_site%pftPlantPopulation       , &
    ZERO     => plt_site%ZERO     , &
    ZSIN     => plt_rad%ZSIN        &
  )
!   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
!
!   HypoctoylHeight=hypocotyledon height
!   SeedinDepth=seeding depth
!   ARLF=node leaf area
!   CanPSheathHeight=petiole length
!
  KVSTGN(NB,NZ)=0
  IF(HypoctoylHeight(NZ).LE.SeedinDepth(NZ).AND.ARLF1(0,NB1(NZ),NZ).GT.0.0)THEN
    XLGLF=SQRT(1.0E+02*ARLF1(0,NB1(NZ),NZ)/pftPlantPopulation(NZ))
    HypoctoylHeight(NZ)=XLGLF+CanPSheathHeight(0,NB1(NZ),NZ)+HTNODE(0,NB1(NZ),NZ)
  ENDIF
!
! IF CANOPY HAS EMERGED
!
  IF(HypoctoylHeight(NZ).GT.SeedinDepth(NZ))THEN
    D540: DO K=0,JNODS1
      DO  L=1,NumOfCanopyLayers1
        CanPLNBLA(L,K,NB,NZ)=0._r8
        WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ)=0._r8
      enddo
    ENDDO D540
    D535: DO L=1,NumOfCanopyLayers1
      CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
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
        IF(BranchNumber_brchpft(NB,NZ).GE.KVSTG1)THEN
          K=MOD(BranchNumber_brchpft(NB,NZ),JNODS1)
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
      HTLF=HTSTK+CanPSheathHeight(K,NB,NZ)
      XLGLF=AZMAX1(SQRT(WDLF(NZ)*AMAX1(0.0 &
        ,ARLF1(K,NB,NZ))/(pftPlantPopulation(NZ)*FNOD(NZ))))
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
        D550: DO L=NumOfCanopyLayers1,1,-1
          IF(LU.EQ.1.AND.LL.EQ.1)exit
          IF((HTLFU.GT.CanopyHeightz(L-1).OR.CanopyHeightz(L-1).LE.ZERO).AND.LU.EQ.0)THEN
            LHTLFU=MAX(1,L)
            LU=1
          ENDIF
          IF((HTLFL.GT.CanopyHeightz(L-1).OR.CanopyHeightz(L-1).LE.ZERO).AND.LL.EQ.0)THEN
            LHTLFL=MAX(1,L)
            LL=1
          ENDIF
        ENDDO D550
        D570: DO L=LHTLFL,LHTLFU
          IF(LHTLFU.EQ.LHTLFL)THEN
            FRACL=CLASS(N,NZ)
          ELSEIF(HTLFU.GT.HTLFL.AND.CanopyHeightz(L).GT.HTLFL)THEN
            FRACL=CLASS(N,NZ)*(AMIN1(HTLFU,CanopyHeightz(L)) &
              -AMAX1(HTLFL,CanopyHeightz(L-1)))/(HTLFU-HTLFL)
          ELSE
            FRACL=CLASS(N,NZ)
          ENDIF
          YARLF=FRACL*ARLF1(K,NB,NZ)
          DO NE=1,NumOfPlantChemElements
            YWGLFE(NE)=FRACL*WGLFE(NE,K,NB,NZ)
          ENDDO
!
    !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
    !
    !     CanPLNBLA=leaf node area in canopy layer
    !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
    !     CanopyLeafApft_lyr,CanopyLeafCpft_lyr=total leaf area,C in canopy layer
    !     HTNODE=internode length
    !
          CanPLNBLA(L,K,NB,NZ)=CanPLNBLA(L,K,NB,NZ)+YARLF
          DO NE=1,NumOfPlantChemElements
            WGLFLE(NE,L,K,NB,NZ)=WGLFLE(NE,L,K,NB,NZ)+YWGLFE(NE)
          ENDDO
          CanopyLeafApft_lyr(L,NZ)=CanopyLeafApft_lyr(L,NZ)+YARLF
          CanopyLeafCpft_lyr(L,NZ)=CanopyLeafCpft_lyr(L,NZ)+YWGLFE(ielmc)
        ENDDO D570
        TLGLF=TLGLF+YLGLF
        CanopyHeight(NZ)=AMAX1(CanopyHeight(NZ),HTLFU)
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
    IF(isclose(HTNODE(K1,NB,NZ),0._r8))THEN
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
  !     CanPBStalkC=stalk sapwood mass
  !     FRACL=stalk fraction in each layer
  !     CanopyBranchStemApft_lyr=total branch stalk surface area in each layer
  !
  !     IF(NZ.EQ.1)THEN
  !     WRITE(*,6679)'K1',I,J,NZ,NB,K1,KVSTG(NB,NZ)
  !    2,HTNODE(K1,NB,NZ)
  !6679  FORMAT(A8,6I4,12E12.4)
  !     ENDIF
    IF(HTNODE(K1,NB,NZ).GT.0.0)THEN
      LU=0
      LL=0
      D545: DO L=NumOfCanopyLayers1,1,-1
        IF(LU.EQ.1.AND.LL.EQ.1)exit
        IF((HTLFB.GT.CanopyHeightz(L-1).OR.CanopyHeightz(L-1).LE.ZERO).AND.LU.EQ.0)THEN
          LHTBRU=MAX(1,L)
          LU=1
        ENDIF
        IF((HTBR.GT.CanopyHeightz(L-1).OR.CanopyHeightz(L-1).LE.ZERO).AND.LL.EQ.0)THEN
          LHTBRL=MAX(1,L)
          LL=1
        ENDIF
      ENDDO D545
      RSTK=SQRT(VSTK*(AZMAX1(WTSTKBE(ielmc,NB,NZ))/pftPlantPopulation(NZ))/(PICON*HTNODE(K1,NB,NZ)))
      ARSTKB=PICON*HTNODE(K1,NB,NZ)*pftPlantPopulation(NZ)*RSTK
      IF(ISTYP(NZ).EQ.iplt_annual)THEN
        CanPBStalkC(NB,NZ)=WTSTKBE(ielmc,NB,NZ)
      ELSE
        ZSTK=AMIN1(ZSTX,FSTK*RSTK)
        ASTV=PICON*(2.0_r8*RSTK*ZSTK-ZSTK**2)
        CanPBStalkC(NB,NZ)=ASTV/VSTK*HTNODE(K1,NB,NZ)*pftPlantPopulation(NZ)
      ENDIF
    !     IF(NZ.EQ.1)THEN
        !     WRITE(*,6677)'WVSTK',I,J,NX,NY,NZ,NB,K1,CanPBStalkC(NB,NZ)
        !    2,ASTV,VSTK,HTNODE(K1,NB,NZ),pftPlantPopulation(NZ)
        !6677  FORMAT(A8,7I4,12E12.4)
        !     ENDIF
      D445: DO L=LHTBRL,LHTBRU
        IF(HTLFB.GT.HTBR)THEN
          IF(HTLFB.GT.CanopyHeightz(L-1))THEN
            FRACL=(AMIN1(HTLFB,CanopyHeightz(L))-AMAX1(HTBR &
              ,CanopyHeightz(L-1)))/(HTLFB-HTBR)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        CanopyBranchStemApft_lyr(L,NB,NZ)=FRACL*ARSTKB
      ENDDO D445
    ELSE
      CanPBStalkC(NB,NZ)=0._r8
      D450: DO L=1,NumOfCanopyLayers1
        CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
      ENDDO D450
    ENDIF
  ELSE
    CanPBStalkC(NB,NZ)=0._r8
    D455: DO L=1,NumOfCanopyLayers1
      CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
    ENDDO D455
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
    CanopyBranchStemApft_lyr   =>  plt_morph%CanopyBranchStemApft_lyr    , &
    CLASS   =>  plt_morph%CLASS    , &
    CanPLNBLA   =>  plt_morph%CanPLNBLA    , &
    StemA_lyrnodbrchpft   =>  plt_morph%StemA_lyrnodbrchpft    , &
    NB1     =>  plt_morph%NB1      , &
    LeafA_lyrnodbrchpft    =>  plt_morph%LeafA_lyrnodbrchpft       &
  )
  D900: DO K=1,JNODS1
    DO  L=1,NumOfCanopyLayers1
      DO  N=1,JLI1
        LeafA_lyrnodbrchpft(N,L,K,NB,NZ)=0._r8
      enddo
    enddo
  ENDDO D900
! ARLFXB=0._r8
! ARLFXL=0._r8
! SURFXX=0._r8
  D500: DO K=1,JNODS1
!     ARLFXB=ARLFXB+ARLF1(K,NB,NZ)
    IF(ARLF1(K,NB,NZ).GT.0.0)THEN
      D700: DO L=NumOfCanopyLayers1,1,-1
!       ARLFXL=ARLFXL+CanPLNBLA(L,K,NB,NZ)
        D800: DO N=1,JLI1
          LeafA_lyrnodbrchpft(N,L,K,NB,NZ)=AZMAX1(CLASS(N,NZ)*0.25_r8*CanPLNBLA(L,K,NB,NZ))
  !       SURFXX=SURFXX+LeafA_lyrnodbrchpft(N,L,K,NB,NZ)
        ENDDO D800
      ENDDO D700
    ENDIF
  ENDDO D500
!
! ALLOCATE STALK AREA TO INCLINATION CLASSES ACCORDING TO
! BRANCH ANGLE ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
!
! StemA_lyrnodbrchpft=stalk surface area in canopy layer
! ANGBR=stem angle from horizontal
! CanopyBranchStemApft_lyr=total branch stalk surface area in each layer
!
  D910: DO L=1,NumOfCanopyLayers1
    DO  N=1,JLI1
      StemA_lyrnodbrchpft(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D910

  IF(NB.EQ.NB1(NZ))THEN
    N=JLI1
  ELSE
    dangle=PICON2h/real(JLI1,r8)
    N=MIN(JLI1,INT(ASIN(ANGBR(NZ))/dangle)+1)
  ENDIF
  D710: DO L=NumOfCanopyLayers1,1,-1
    StemA_lyrnodbrchpft(N,L,NB,NZ)=CanopyBranchStemApft_lyr(L,NB,NZ)/real(JLI1,r8)
  ENDDO D710
  end associate
  end subroutine LeafClassAllocation

!------------------------------------------------------------------------------------------

  subroutine GrainFilling(I,NB,NZ,GROGRE,GROSTKC)
  implicit none
  integer, intent(in) :: I,NB,NZ
  real(r8), intent(in) :: GROGRE(NumOfPlantChemElements),GROSTKC
  real(r8) :: ZPGRP,ZPGRN,ZNPGP,ZNPGN
  real(r8) :: XLOCM,XLOCE(NumOfPlantChemElements)
  real(r8) :: FGRNX
  real(r8) :: GRMXB
  real(r8) :: GROLM
  real(r8) :: GROLC
  real(r8) :: SeedSET
  integer :: NE
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
    fTgrowRootP     =>  plt_pheno%fTgrowRootP     , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP     , &
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
    pftPlantPopulation       =>  plt_site%pftPlantPopulation        , &
    SDMX     =>  plt_morph%SDMX     , &
    MaxSeedCMass    =>  plt_morph%MaxSeedCMass    , &
    STMX     =>  plt_morph%STMX     , &
    NGTopRootLayer      =>  plt_morph%NGTopRootLayer      , &
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
!     GROSTKC=stalk growth rate
!
  IF(IDAY(3,NB,NZ).NE.0.AND.IDAY(6,NB,NZ).EQ.0)THEN
    GRNXB(NB,NZ)=GRNXB(NB,NZ)+STMX(NZ)*AZMAX1(GROSTKC)
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
    SeedSET=AMIN1(CEPOLB(ielmc,NB,NZ)/(CEPOLB(ielmc,NB,NZ)+SETC) &
      ,CEPOLB(ielmn,NB,NZ)/(CEPOLB(ielmn,NB,NZ)+SETN) &
      ,CEPOLB(ielmp,NB,NZ)/(CEPOLB(ielmp,NB,NZ)+SETP))

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
        *SeedSET*DGSTGF(NB,NZ)-FGRNX*GRNOB(NB,NZ))
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
      GRMXB=MaxSeedCMass(NZ)
      GRWTB(NB,NZ)=AMIN1(MaxSeedCMass(NZ),GRWTB(NB,NZ) &
        +GRMXB*AMAX1(0.50_r8,SeedSET**0.25_r8)*DGSTGF(NB,NZ))
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
!   fTgrowCanP=temperature function for canopy growth
!   fTgrowRootP=temperature function for root growth
!
  IF(IDAY(7,NB,NZ).NE.0)THEN
    IF(WTGRBE(ielmc,NB,NZ).GE.GRWTB(NB,NZ)*GRNOB(NB,NZ))THEN
      GROLM=0._r8
    ELSEIF(IRTYP(NZ).EQ.0)THEN
      GROLM=AZMAX1(GFILL(NZ)*GRNOB(NB,NZ)*SQRT(fTgrowCanP(NZ)))
    ELSE
      GROLM=AZMAX1(GFILL(NZ)*GRNOB(NB,NZ)*SQRT(fTgrowRootP(NGTopRootLayer(NZ),NZ)))
    ENDIF
!
!     GRAIN FILL RATE MAY BE CONSTRAINED BY HIGH GRAIN C:N OR C:P
!
!     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     GROLM,GROLC=maximum,actual grain fill rate
!     XLOCM,XLOCE(ielmc)=maximum,actual C translocation rate from reserve to grain
!
    IF(WTGRBE(ielmn,NB,NZ).LT.ZPGRM*CNGR(NZ) &
      *WTGRBE(ielmc,NB,NZ).OR.WTGRBE(ielmp,NB,NZ).LT.ZPGRM &
      *CPGR(NZ)*WTGRBE(ielmc,NB,NZ))THEN
      GROLC=0._r8
    ELSE
      GROLC=GROLM
    ENDIF
    XLOCM=AMIN1(GROLM,WTRSVBE(ielmc,NB,NZ))
    XLOCE(ielmc)=AMIN1(GROLC,WTRSVBE(ielmc,NB,NZ))
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
!     XLOCM,XLOCE(ielmc)=maximum,actual C translocation rate from reserve to grain
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     XLOCE(ielmn),XLOCE(ielmp)=N,P translocation rate from reserve to grain
!
    IF(WTRSVBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      ZNPGN=WTRSVBE(ielmn,NB,NZ)/(WTRSVBE(ielmn,NB,NZ) &
        +SETN*WTRSVBE(ielmc,NB,NZ))
      ZNPGP=WTRSVBE(ielmp,NB,NZ)/(WTRSVBE(ielmp,NB,NZ) &
        +SETP*WTRSVBE(ielmc,NB,NZ))
      ZPGRN=ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGN))
      ZPGRP=ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGP))
      XLOCE(ielmn)=AMIN1(XLOCM*CNGR(NZ) &
        ,AZMAX1(WTRSVBE(ielmn,NB,NZ)*ZPGRN) &
        ,(WTGRBE(ielmc,NB,NZ)+XLOCE(ielmc))*CNGR(NZ)-WTGRBE(ielmn,NB,NZ))
      XLOCE(ielmp)=AMIN1(XLOCM*CPGR(NZ) &
        ,AZMAX1(WTRSVBE(ielmp,NB,NZ)*ZPGRP) &
        ,(WTGRBE(ielmc,NB,NZ)+XLOCE(ielmc))*CPGR(NZ)-WTGRBE(ielmp,NB,NZ))
    ELSE
      XLOCE(ielmn)=0._r8
      XLOCE(ielmp)=0._r8
    ENDIF
!
!     TRANSLOCATE C,N,P FROM STALK RESERVES TO GRAIN
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     GROGR=grain growth rate
!     XLOCE(ielmc),XLOCE(ielmn),XLOCE(ielmp)=C,N,P translocation rate from reserve to grain
!
    DO NE=1,NumOfPlantChemElements
      WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)+GROGRE(NE)-XLOCE(NE)
      WTGRBE(NE,NB,NZ)=WTGRBE(NE,NB,NZ)+XLOCE(NE)
    ENDDO
  ELSE
    XLOCE(1:NumOfPlantChemElements)=0._r8
  ENDIF
!
!   SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
!   HAS STOPPED FOR SET PERIOD OF TIME
!
!   IDAY(8,=end date setting for final seed number
!   XLOCE(ielmc)=C translocation rate from reserve to grain
!   PP=PFT population
!   FLG4=number of hours with no grain fill
!   Hours4PhyslMature=number of hours with no grain filling until physl maturity
!   IDAY(10,=date of physiological maturity
!
  IF(IDAY(8,NB,NZ).NE.0)THEN
    IF(XLOCE(ielmc).LE.1.0E-09*pftPlantPopulation(NZ))THEN
      FLG4(NB,NZ)=FLG4(NB,NZ)+1.0
    ELSE
      FLG4(NB,NZ)=0._r8
    ENDIF
    IF(FLG4(NB,NZ).GE.Hours4PhyslMature)THEN
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
!     Hours4PhyslMature=number of hours with no grain filling until physiological maturity
!     Hours4SenesAftMature=number of hours after physiol maturity required for senescence
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!
    IF(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
      IF(FLG4(NB,NZ).GT.Hours4PhyslMature+Hours4SenesAftMature(IWTYP(NZ)))THEN
        VRNF(NB,NZ)=VRNX(NB,NZ)+0.5_r8
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine GrainFilling
!------------------------------------------------------------------------------------------

  subroutine PhenologyReset(I,NB,NZ)

  implicit none
  integer, intent(in) :: I,nb,nz
  integer :: K,M,NE
  real(r8) :: FSNR1
! begin_execution
  associate(                        &
    instruct =>  pltpar%instruct  , &
    ifoliar  =>  pltpar%ifoliar   , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    infoliar =>  pltpar%infoliar  , &
    icwood   =>  pltpar%icwood    , &
    k_fine_litr => pltpar%k_fine_litr, &
    k_woody_litr=> pltpar%k_woody_litr, &
    EHVST    =>  plt_distb%EHVST     , &
    IYRH     =>  plt_distb%IYRH      , &
    THIN_pft     =>  plt_distb%THIN_pft      , &
    HVST     =>  plt_distb%HVST      , &
    IHVST    =>  plt_distb%IHVST     , &
    JHVST    =>  plt_distb%JHVST     , &
    IDAY0    =>  plt_distb%IDAY0     , &
    IDAYH    =>  plt_distb%IDAYH     , &
    IYR0     =>  plt_distb%IYR0      , &
    WTGRBE   =>  plt_biom%WTGRBE     , &
    WGLFE    =>  plt_biom%WGLFE      , &
    PopPlantRootC_vr   =>  plt_biom%PopPlantRootC_vr     , &
    WTSTXBE  =>  plt_biom%WTSTXBE    , &
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
    WGNODE   =>  plt_biom%WGNODE     , &
    FWODLE   =>  plt_allom%FWODLE    , &
    FWODBE   =>  plt_allom%FWODBE    , &
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
    CFOPE    =>  plt_soilchem%CFOPE  , &
    IYRC     =>  plt_site%IYRC       , &
    ESNC     =>  plt_bgcr%ESNC       , &
    MY       =>  plt_morph%MY        , &
    VSTG     =>  plt_morph%VSTG      , &
    XTLI     =>  plt_morph%XTLI      , &
    CanPSheathHeight    =>  plt_morph%CanPSheathHeight     , &
    KLEAF    =>  plt_morph%KLEAF     , &
    HTNODX   =>  plt_morph%HTNODX    , &
    HTNODE   =>  plt_morph%HTNODE    , &
    GRNXB    =>  plt_morph%GRNXB     , &
    PSTGF    =>  plt_morph%PSTGF     , &
    NB1      =>  plt_morph%NB1       , &
    BranchNumber_brchpft     =>  plt_morph%BranchNumber_brchpft      , &
    CanopyBranchLeafA_pft    =>  plt_morph%CanopyBranchLeafA_pft     , &
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
  IF(ISTYP(NZ).NE.iplt_annual.OR.(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).GT.1))THEN
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
      IF((IFLGE(NB,NZ).EQ.0.AND.ISTYP(NZ).NE.iplt_annual).AND.(VRNS(NB,NZ).GE.VRNL(NB,NZ)))THEN
        IF(ISTYP(NZ).EQ.iplt_annual)THEN
          GROUP(NB,NZ)=AZMAX1(GROUPI(NZ)-BranchNumber_brchpft(NB,NZ))
        ELSE
          GROUP(NB,NZ)=GROUPI(NZ)
        ENDIF
        PSTGI(NB,NZ)=PSTG(NB,NZ)
        PSTGF(NB,NZ)=0._r8
        VSTGX(NB,NZ)=0._r8
        TGSTGI(NB,NZ)=0._r8
        TGSTGF(NB,NZ)=0._r8
        IDAY(1,NB,NZ)=I
        D2005: DO M=2,NumGrothStages
          IDAY(M,NB,NZ)=0
        ENDDO D2005
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
    !     CanopyBranchLeafA_pft,ARLF=branch,node leaf area
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !     HTNODX,HTNODE=stalk height,stalk internode length
    !     GRNOB=seed set number
    !     GRNXB=potential number of seed set sites
    !     GRWTB=individual seed size
!
        IF(IFLGE(NB,NZ).EQ.0.AND.ISTYP(NZ).NE.iplt_annual.AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
          IF(IBTYP(NZ).EQ.0)THEN
            PSTG(NB,NZ)=XTLI(NZ)
            VSTG(NB,NZ)=0._r8
            KLEAF(NB,NZ)=1
            KVSTG(NB,NZ)=1
            FLG4(NB,NZ)=0._r8
            D5330: DO M=1,jsken
              ESNC(ielmc,M,k_woody_litr,0,NZ)=ESNC(ielmc,M,k_woody_litr,0,NZ) &
                +CFOPE(ielmc,icwood,M,NZ)*WTLFBE(ielmc,NB,NZ)*FWODBE(ielmc,k_woody_litr) &
                +CFOPE(ielmc,icwood,M,NZ)*WTSHEBE(ielmc,NB,NZ)*FWODBE(ielmc,k_woody_litr)
              ESNC(ielmn,M,k_woody_litr,0,NZ)=ESNC(ielmn,M,k_woody_litr,0,NZ) &
                +CFOPE(ielmn,icwood,M,NZ)*WTLFBE(ielmn,NB,NZ)*FWODLE(ielmn,k_woody_litr) &
                +CFOPE(ielmn,icwood,M,NZ)*WTSHEBE(ielmn,NB,NZ)*FWODBE(ielmn,k_woody_litr)
              ESNC(ielmp,M,k_woody_litr,0,NZ)=ESNC(ielmp,M,k_woody_litr,0,NZ) &
                +CFOPE(ielmp,icwood,M,NZ)*WTLFBE(ielmp,NB,NZ)*FWODLE(ielmp,k_woody_litr) &
                +CFOPE(ielmp,icwood,M,NZ)*WTSHEBE(ielmp,NB,NZ)*FWODBE(ielmp,k_woody_litr)
              ESNC(ielmc,M,k_fine_litr,0,NZ)=ESNC(ielmc,M,k_fine_litr,0,NZ) &
                +CFOPE(ielmc,ifoliar,M,NZ)*WTLFBE(ielmc,NB,NZ)*FWODBE(ielmc,k_fine_litr) &
                +CFOPE(ielmc,infoliar,M,NZ)*WTSHEBE(ielmc,NB,NZ)*FWODBE(ielmc,k_fine_litr)
              ESNC(ielmn,M,k_fine_litr,0,NZ)=ESNC(ielmn,M,k_fine_litr,0,NZ) &
                +CFOPE(ielmn,ifoliar,M,NZ)*WTLFBE(ielmn,NB,NZ)*FWODLE(ielmn,k_fine_litr) &
                +CFOPE(ielmn,infoliar,M,NZ)*WTSHEBE(ielmn,NB,NZ)*FWODBE(ielmn,k_fine_litr)
              ESNC(ielmp,M,k_fine_litr,0,NZ)=ESNC(ielmp,M,k_fine_litr,0,NZ) &
                +CFOPE(ielmp,ifoliar,M,NZ)*WTLFBE(ielmp,NB,NZ)*FWODLE(ielmp,k_fine_litr) &
                +CFOPE(ielmp,infoliar,M,NZ)*WTSHEBE(ielmp,NB,NZ)*FWODBE(ielmp,1)
            ENDDO D5330
            CanopyBranchLeafA_pft(NB,NZ)=0._r8
            WTLFBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
            WTSHEBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
            D5335: DO K=0,JNODS1
              ARLF1(K,NB,NZ)=0._r8
              CanPSheathHeight(K,NB,NZ)=0._r8
              WSLF(K,NB,NZ)=0._r8
              WGLFE(1:NumOfPlantChemElements,K,NB,NZ)=0._r8
              WGSHE(1:NumOfPlantChemElements,K,NB,NZ)=0._r8
              WSSHE(K,NB,NZ)=0._r8
            ENDDO D5335
          ENDIF
        ENDIF
    !
    !     RESIDUAL STALKS BECOME LITTERFALL IN GRASSES, SHRUBS AT
    !     START OF SEASON
    !
        IF((IFLGE(NB,NZ).EQ.0.AND.ISTYP(NZ).NE.iplt_annual).AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN

          D6245: DO M=1,jsken
            DO NE=1,NumOfPlantChemElements
              ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,infoliar,M,NZ) &
                *(WTHSKBE(NE,NB,NZ)+WTEARBE(NE,NB,NZ)+WTGRBE(NE,NB,NZ))
            ENDDO
          ENDDO D6245
            
          WTHSKBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
          WTEARBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
          WTGRBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
          
          GRNXB(NB,NZ)=0._r8
          GRNOB(NB,NZ)=0._r8
          GRWTB(NB,NZ)=0._r8
          IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
            D6345: DO M=1,jsken
              DO NE=1,NumOfPlantChemElements
                ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,istalk,M,NZ)*WTSTKBE(NE,NB,NZ)
              ENDDO
            ENDDO D6345
            WTSTKBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
            WTSTXBE(1:NumOfPlantChemElements,NB,NZ)=0._r8
            DO K=0,JNODS1
              DO NE=1,NumOfPlantChemElements
                WGNODE(NE,K,NB,NZ)=0._r8
              ENDDO
            ENDDO
            D6340: DO K=0,JNODS1
              HTNODE(K,NB,NZ)=0._r8
              HTNODX(K,NB,NZ)=0._r8
            ENDDO D6340
          ENDIF
        ENDIF
      ENDIF
!
  !   SPRING OR FALL FLAG RESET
  !
      IF(IFLGE(NB,NZ).EQ.0.AND.VRNS(NB,NZ).GE.VRNL(NB,NZ))THEN
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
    FSNR1=1.0_r8-FSNR

    D6330: DO M=1,jsken
      DO NE=1,NumOfPlantChemElements
        ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+FSNR*CFOPE(NE,infoliar,M,NZ) &
          *(WTHSKBE(NE,NB,NZ)+WTEARBE(NE,NB,NZ))
        IF(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
          WTRVE(NE,NZ)=WTRVE(NE,NZ)+FSNR*CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
        ELSE
          ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+FSNR*CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
        ENDIF
      ENDDO
    ENDDO D6330

    DO NE=1,NumOfPlantChemElements  
      WTHSKBE(NE,NB,NZ)=FSNR1*WTHSKBE(NE,NB,NZ)
      WTEARBE(NE,NB,NZ)=FSNR1*WTEARBE(NE,NB,NZ)
      WTGRBE(NE,NB,NZ)=FSNR1*WTGRBE(NE,NB,NZ)
    ENDDO
    GRNXB(NB,NZ)=FSNR1*GRNXB(NB,NZ)
    GRNOB(NB,NZ)=FSNR1*GRNOB(NB,NZ)
    GRWTB(NB,NZ)=FSNR1*GRWTB(NB,NZ)
!
!     STALKS BECOME LITTERFALL IN GRASSES AT END OF SEASON
!
    IF((IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1).AND.ISTYP(NZ).NE.iplt_annual)THEN

      D6335: DO M=1,jsken
        DO NE=1,NumOfPlantChemElements
          ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+FSNR*CFOPE(NE,istalk,M,NZ)*WTSTKBE(NE,NB,NZ)
        ENDDO
      ENDDO D6335
      DO NE=1,NumOfPlantChemElements  
        WTSTKBE(NE,NB,NZ)=FSNR1*WTSTKBE(NE,NB,NZ)
        WTSTXBE(NE,NB,NZ)=FSNR1*WTSTXBE(NE,NB,NZ)
      ENDDO
      DO K=0,JNODS1
        DO NE=1,NumOfPlantChemElements
          WGNODE(NE,K,NB,NZ)=FSNR1*WGNODE(NE,K,NB,NZ)
        ENDDO
      ENDDO
      D2010: DO K=0,JNODS1
    !     HTNODE(K,NB,NZ)=FSNR1*HTNODE(K,NB,NZ)
        HTNODX(K,NB,NZ)=FSNR1*HTNODX(K,NB,NZ)
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
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
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
      !deciduous annual plant
      IF(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
        IDAYH(NZ)=I
        IYRH(NZ)=IYRC
        IHVST(NZ)=1
        JHVST(NZ)=ihv_tmareseed
        HVST(NZ)=0._r8
        THIN_pft(NZ)=0._r8
        EHVST(1,ipld_leaf,NZ)=1.0_r8
        EHVST(1,ipld_nofoliar,NZ)=1.0_r8
        EHVST(1,ipld_woody,NZ)=1.0_r8
        EHVST(1,ipld_stdead,NZ)=1.0_r8
        EHVST(2,ipld_leaf,NZ)=0._r8
        EHVST(2,ipld_nofoliar,NZ)=1.0_r8
        EHVST(2,ipld_woody,NZ)=0._r8
        EHVST(2,ipld_stdead,NZ)=0._r8
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
  integer :: L,NE
  real(r8) :: ZPOOLD
  real(r8) :: XFRPX,XFRCX,XFRNX
  real(r8) :: ATRPPD
  REAL(R8) :: cpoolt
  real(r8) :: CH2OH
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
  real(r8) :: EPOOLM(1:NumOfPlantChemElements)
  real(r8) :: UPNH4B,UPPO4B
  real(r8) :: UPNH4R,UPPO4R
  real(r8) :: WTRTM,WFNSP
  real(r8) :: WTPLTT,WTRTRX
  real(r8) :: WTPLTX,WVSTBX
  real(r8) :: WTRTTX,WTRSBX
  real(r8) :: WTRVCX
  real(r8) :: XFRE(1:NumOfPlantChemElements)
  ! begin_execution
  associate(                          &
    IDAY0  =>  plt_distb%IDAY0  , &
    IYR0   =>  plt_distb%IYR0   , &
    IYRC   =>  plt_site%IYRC    , &
    ZEROS2 => plt_site%ZEROS2   , &
    DYLN   =>  plt_site%DYLN    , &
    NU     =>  plt_site%NU      , &
    FVRN   =>  plt_allom%FVRN   , &
    FWODRE =>  plt_allom%FWODRE , &
    WTRTE  =>  plt_biom%WTRTE   , &
    PopPlantRootC_vr =>  plt_biom%PopPlantRootC_vr  , &
    WTRTL  =>  plt_biom%WTRTL   , &
    CanPStalkC  =>  plt_biom%CanPStalkC   , &
    EPOOL  =>  plt_biom%EPOOL   , &
    EPOOLR =>  plt_biom%EPOOLR  , &
    WTRVE  =>  plt_biom%WTRVE   , &
    CEPOLB =>  plt_biom%CEPOLB  , &
    CanPBLeafShethC  =>  plt_biom%CanPBLeafShethC   , &
    ZEROP  =>  plt_biom%ZEROP   , &
    WTRSVBE=>  plt_biom%WTRSVBE , &
    CanPBStalkC =>  plt_biom%CanPBStalkC  , &
    VLSoilPoreMicP   =>  plt_soilchem%VLSoilPoreMicP, &
    IDAY   =>  plt_pheno%IDAY   , &
    fTgrowCanP   =>  plt_pheno%fTgrowCanP   , &
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
    HourCounter4LeafOut_brch   =>   plt_pheno%HourCounter4LeafOut_brch  , &
    NGTopRootLayer    =>   plt_morph%NGTopRootLayer   , &
    NB1    =>  plt_morph%NB1    , &
    NI     =>  plt_morph%NI       &
  )
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!

  IF((ISTYP(NZ).EQ.iplt_annual.AND.IFLGI(NZ).EQ.0) &
    .OR.(I.GE.IDAY0(NZ).AND.IYRC.EQ.IYR0(NZ) &
    .AND.VRNF(NB,NZ).LT.FVRN(IWTYP(NZ))*VRNX(NB,NZ)) &
    .OR.(VRNS(NB1(NZ),NZ).GE.VRNL(NB,NZ) &
    .AND.VRNF(NB,NZ).LT.FVRN(IWTYP(NZ))*VRNX(NB,NZ)))THEN
    WTRTM=0._r8
    EPOOLM(ielmc)=0._r8
    D4: DO L=NU,NI(NZ)
      WTRTM=WTRTM+AZMAX1(PopPlantRootC_vr(ipltroot,L,NZ))
      EPOOLM(ielmc)=EPOOLM(ielmc)+AZMAX1(EPOOLR(ielmc,ipltroot,L,NZ))
    ENDDO D4
!
  ! RESET TIME COUNTER
  !
  ! HourCounter4LeafOut_brch=hourly leafout counter
  ! IFLGA=flag for initializing leafout
  !
    IF(IFLGA(NB,NZ).EQ.0)THEN
      HourCounter4LeafOut_brch(NB,NZ)=0._r8
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
  ! fTgrowCanP=temperature function for canopy growth
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
      DATRP=ATRPPD*fTgrowCanP(NZ)*WFNSP
      HourCounter4LeafOut_brch(NB,NZ)=HourCounter4LeafOut_brch(NB,NZ)+DATRP
      IF(HourCounter4LeafOut_brch(NB,NZ).LE.ATRPX(ISTYP(NZ)) &
        .OR.(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).EQ.0))THEN
        IF(WTRVE(ielmc,NZ).GT.ZEROP(NZ))THEN
          CPOOLT=EPOOLM(ielmc)+EPOOL(ielmc,NB,NZ)
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
          EPOOL(ielmc,NB,NZ)=EPOOL(ielmc,NB,NZ)+CH2OH*FXSH(ISTYP(NZ))
          IF(WTRTM.GT.ZEROP(NZ).AND.EPOOLM(ielmc).GT.ZEROP(NZ))THEN
            D50: DO L=NU,NI(NZ)
              FXFC=AZMAX1(PopPlantRootC_vr(ipltroot,L,NZ))/WTRTM
              EPOOLR(ielmc,ipltroot,L,NZ)=EPOOLR(ielmc,ipltroot,L,NZ)+FXFC*CH2OH*FXRT(ISTYP(NZ))
            ENDDO D50
          ELSE
            EPOOLR(ielmc,ipltroot,NGTopRootLayer(NZ),NZ)=EPOOLR(ielmc,ipltroot,NGTopRootLayer(NZ),NZ)+CH2OH*FXRT(ISTYP(NZ))
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
        IF(ISTYP(NZ).NE.iplt_annual)THEN
          CPOOLT=AZMAX1(WTRVE(ielmc,NZ)+EPOOL(ielmc,NB,NZ))
          ZPOOLD=(WTRVE(ielmn,NZ)*EPOOL(ielmc,NB,NZ)-EPOOL(ielmn,NB,NZ)*WTRVE(ielmc,NZ))/CPOOLT
          PPOOLD=(WTRVE(ielmp,NZ)*EPOOL(ielmc,NB,NZ)-EPOOL(ielmp,NB,NZ)*WTRVE(ielmc,NZ))/CPOOLT
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

      DO NE=1,NumOfPlantChemElements
        EPOOLM(NE)=0._r8
        D3: DO L=NU,NI(NZ)
          EPOOLM(NE)=EPOOLM(NE)+AZMAX1(EPOOLR(NE,1,L,NZ))
        ENDDO D3
      ENDDO
      IF(WTRVE(ielmc,NZ).GT.ZEROP(NZ))THEN
        IF(ISTYP(NZ).NE.iplt_annual)THEN
          CPOOLT=AMAX1(ZEROP(NZ),WTRVE(ielmc,NZ)+EPOOLM(ielmc))
          ZPOOLD=(WTRVE(ielmn,NZ)*EPOOLM(ielmc)-EPOOLM(ielmn)*WTRVE(ielmc,NZ))/CPOOLT
          PPOOLD=(WTRVE(ielmp,NZ)*EPOOLM(ielmc)-EPOOLM(ielmp)*WTRVE(ielmc,NZ))/CPOOLT
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
      EPOOL(ielmn,NB,NZ)=EPOOL(ielmn,NB,NZ)+UPNH4B
      EPOOL(ielmp,NB,NZ)=EPOOL(ielmp,NB,NZ)+UPPO4B
      IF(WTRTM.GT.ZEROP(NZ).AND.EPOOLM(ielmc).GT.ZEROP(NZ))THEN
        D51: DO L=NU,NI(NZ)
          FXFN=AZMAX1(EPOOLR(ielmc,ipltroot,L,NZ))/EPOOLM(ielmc)

          EPOOLR(ielmn,ipltroot,L,NZ)=EPOOLR(ielmn,ipltroot,L,NZ)+FXFN*UPNH4R
          EPOOLR(ielmp,ipltroot,L,NZ)=EPOOLR(ielmp,ipltroot,L,NZ)+FXFN*UPPO4R
        ENDDO D51
      ELSE

        EPOOLR(ielmn,ipltroot,NGTopRootLayer(NZ),NZ)=EPOOLR(ielmn,ipltroot,NGTopRootLayer(NZ),NZ)+UPNH4R
        EPOOLR(ielmp,ipltroot,NGTopRootLayer(NZ),NZ)=EPOOLR(ielmp,ipltroot,NGTopRootLayer(NZ),NZ)+UPPO4R
      ENDIF
    ENDIF
  !
  ! REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
  !
  ! ATRP=hourly leafout counter
  ! fTgrowCanP=temperature function for canopy growth
  ! ATRPX=number of hours required for remobilization of storage C during leafout
  ! WFNG=growth function of canopy water potential
  ! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  ! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
  !
    IF(NB.NE.NB1(NZ).AND.HourCounter4LeafOut_brch(NB,NZ).LE.ATRPX(ISTYP(NZ)))THEN
      HourCounter4LeafOut_brch(NB,NZ)=HourCounter4LeafOut_brch(NB,NZ)+fTgrowCanP(NZ)*WFNG
      DO NE=1,NumOfPlantChemElements
        XFRE(NE)=AZMAX1(0.05_r8*fTgrowCanP(NZ) &
          *(0.5_r8*(EPOOL(NE,NB1(NZ),NZ)+EPOOL(NE,NB,NZ))-EPOOL(NE,NB,NZ)))
        EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+XFRE(NE)
        EPOOL(NE,NB1(NZ),NZ)=EPOOL(NE,NB1(NZ),NZ)-XFRE(NE)
      ENDDO
    ENDIF
  ENDIF

!
! TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
! IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
! IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
!
! ISTYP=growth habit:0=annual,1=perennial from PFT file
! IFLGZ=remobilization flag
! CanPBStalkC=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! FXFB=rate constant for plant-storage nonstructural C,N,P exchange
! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! ISTYP=growth habit:0=annual,1=perennial from PFT file
! IFLGZ=remobilization flag
! CanPBStalkC=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
! WTRVC,WTRVN,WTRVP=storage C,N,P
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!
  IF(IFLGZ.EQ.1.AND.ISTYP(NZ).NE.iplt_annual)THEN
    IF(CanPBStalkC(NB,NZ).GT.ZEROP(NZ).AND.WTRSVBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      CWTRSV=AZMAX1(WTRSVBE(ielmc,NB,NZ)/CanPBStalkC(NB,NZ))
      CWTRSN=AZMAX1(WTRSVBE(ielmn,NB,NZ)/CanPBStalkC(NB,NZ))
      CWTRSP=AZMAX1(WTRSVBE(ielmp,NB,NZ)/CanPBStalkC(NB,NZ))
      CNR=CWTRSV/(CWTRSV+CWTRSN/CNKI)
      CPR=CWTRSV/(CWTRSV+CWTRSP/CPKI)
    ELSE
      CNR=0._r8
      CPR=0._r8
    ENDIF
    XFRCX=FXFB(IBTYP(NZ))*AZMAX1(WTRSVBE(ielmc,NB,NZ))
    XFRNX=FXFB(IBTYP(NZ))*AZMAX1(WTRSVBE(ielmn,NB,NZ))*(1.0+CNR)
    XFRPX=FXFB(IBTYP(NZ))*AZMAX1(WTRSVBE(ielmp,NB,NZ))*(1.0+CPR)
    XFRE(ielmc)=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRE(ielmn)=AMIN1(XFRNX,XFRE(ielmc)*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRE(ielmp)=AMIN1(XFRPX,XFRE(ielmc)*CPMX,XFRNX*CPMX/CNMN*0.5)
    DO NE=1,NumOfPlantChemElements
      WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)-XFRE(NE)
      WTRVE(NE,NZ)=WTRVE(NE,NZ)+XFRE(NE)
    ENDDO
    IF(CEPOLB(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      CNL=CEPOLB(ielmc,NB,NZ)/(CEPOLB(ielmc,NB,NZ)+CEPOLB(ielmn,NB,NZ)/CNKI)
      CPL=CEPOLB(ielmc,NB,NZ)/(CEPOLB(ielmc,NB,NZ)+CEPOLB(ielmp,NB,NZ)/CPKI)
    ELSE
      CNL=0._r8
      CPL=0._r8
    ENDIF
    XFRCX=FXFB(IBTYP(NZ))*AZMAX1(EPOOL(ielmc,NB,NZ))
    XFRNX=FXFB(IBTYP(NZ))*AZMAX1(EPOOL(ielmn,NB,NZ))*(1.0+CNL)
    XFRPX=FXFB(IBTYP(NZ))*AZMAX1(EPOOL(ielmp,NB,NZ))*(1.0+CPL)
    XFRE(ielmc)=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRE(ielmn)=AMIN1(XFRNX,XFRE(ielmc)*CNMX,XFRPX*CNMX/CPMN*0.5_r8)
    XFRE(ielmp)=AMIN1(XFRPX,XFRE(ielmc)*CPMX,XFRNX*CPMX/CNMN*0.5_r8)
    DO NE=1,NumOfPlantChemElements
      EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)-XFRE(NE)
      WTRVE(NE,NZ)=WTRVE(NE,NZ)+XFRE(NE)
    ENDDO
  ENDIF
!
!   TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES AND ROOTS TO RESERVES
!   IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN STALK RESERVES
!   AND LEAVES IN PERENNIALS ACCORDING TO CONCENTRATION DIFFERENCES
!
!   ISTYP=growth habit:0=annual,1=perennial from PFT file
!   IDAY(3,=start of stem elongation and setting max seed number
!   IDAY(8,=end date setting for final seed number
!   CanPBLeafShethC=leaf+petiole mass
!   CanPBStalkC=stalk sapwood mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   FXFY,FXFZ=rate constant for plant-reserve nonstructural C,N,P exchange
!   XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
!   CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  IF((ISTYP(NZ).EQ.iplt_annual.AND.IDAY(8,NB,NZ).NE.0) &
    .OR.(ISTYP(NZ).EQ.iplt_preanu.AND.IDAY(3,NB,NZ).NE.0))THEN
    WTPLTT=CanPBLeafShethC(NB,NZ)+CanPBStalkC(NB,NZ)
    CPOOLT=EPOOL(ielmc,NB,NZ)+WTRSVBE(ielmc,NB,NZ)
    IF(WTPLTT.GT.ZEROP(NZ))THEN
      CPOOLD=(EPOOL(ielmc,NB,NZ)*CanPBStalkC(NB,NZ) &
        -WTRSVBE(ielmc,NB,NZ)*CanPBLeafShethC(NB,NZ))/WTPLTT
      XFRE(ielmc)=FXFY(ISTYP(NZ))*CPOOLD
      EPOOL(ielmc,NB,NZ)=EPOOL(ielmc,NB,NZ)-XFRE(ielmc)
      WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)+XFRE(ielmc)
    ENDIF
    IF(CPOOLT.GT.ZEROP(NZ))THEN
      ZPOOLD=(EPOOL(ielmn,NB,NZ)*WTRSVBE(ielmc,NB,NZ) &
        -WTRSVBE(ielmn,NB,NZ)*EPOOL(ielmc,NB,NZ))/CPOOLT
      PPOOLD=(EPOOL(ielmp,NB,NZ)*WTRSVBE(ielmc,NB,NZ) &
        -WTRSVBE(ielmp,NB,NZ)*EPOOL(ielmc,NB,NZ))/CPOOLT
      XFRE(ielmn)=FXFZ(ISTYP(NZ))*ZPOOLD
      XFRE(ielmp)=FXFZ(ISTYP(NZ))*PPOOLD
      DO NE=2,NumOfPlantChemElements
        EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)-XFRE(NE)
        WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)+XFRE(NE)
      ENDDO
    ENDIF

    IF(ISTYP(NZ).EQ.iplt_annual.AND.IDAY(8,NB,NZ).NE.0)THEN
      D2050: DO L=NU,NI(NZ)
        IF(VLSoilPoreMicP(L).GT.ZEROS2)THEN
          WTRTRX=AMAX1(ZEROP(NZ),WTRTL(ipltroot,L,NZ)*FWODRE(ielmc,1))
          WTPLTX=WTRTRX+CanPBStalkC(NB,NZ)
          IF(WTPLTX.GT.ZEROP(NZ))THEN
            CPOOLD=(EPOOLR(ielmc,ipltroot,L,NZ)*CanPBStalkC(NB,NZ)-WTRSVBE(ielmc,NB,NZ)*WTRTRX)/WTPLTX
            XFRE(ielmc)=AZMAX1(FXFY(ISTYP(NZ))*CPOOLD)
            EPOOLR(ielmc,ipltroot,L,NZ)=EPOOLR(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
            WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)+XFRE(ielmc)
            CPOOLT=EPOOLR(ielmc,ipltroot,L,NZ)+WTRSVBE(ielmc,NB,NZ)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              ZPOOLD=(EPOOLR(ielmn,ipltroot,L,NZ)*WTRSVBE(ielmc,NB,NZ)-WTRSVBE(ielmn,NB,NZ)*EPOOLR(ielmc,ipltroot,L,NZ))/CPOOLT
              PPOOLD=(EPOOLR(ielmp,ipltroot,L,NZ)*WTRSVBE(ielmc,NB,NZ)-WTRSVBE(ielmp,NB,NZ)*EPOOLR(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(ielmn)=AZMAX1(FXFZ(ISTYP(NZ))*ZPOOLD)
              XFRE(ielmp)=AZMAX1(FXFZ(ISTYP(NZ))*PPOOLD)
              DO NE=2,NumOfPlantChemElements
                EPOOLR(NE,1,L,NZ)=EPOOLR(NE,1,L,NZ)-XFRE(NE)
                WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)+XFRE(NE)
              ENDDO
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
!   CanPBStalkC,WVSTK=stalk,total stalk sapwood mass
!   WTRT=total root mass
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   XFRX=maximum storage C content for remobiln from stalk,root reserves
!   XFRE(ielmc)=C transfer
!
  IF(CanPBStalkC(NB,NZ).GT.ZEROP(NZ).AND.CanPStalkC(NZ).GT.ZEROP(NZ) &
    .AND.WTRTE(ielmc,NZ).GT.ZEROP(NZ) &
    .AND.WTRSVBE(ielmc,NB,NZ).LE.XFRX*CanPBStalkC(NB,NZ))THEN
    FWTBR=CanPBStalkC(NB,NZ)/CanPStalkC(NZ)
    WVSTBX=CanPBStalkC(NB,NZ)
    WTRTTX=WTRTE(ielmc,NZ)*FWTBR
    WTPLTT=WVSTBX+WTRTTX
    WTRSBX=AZMAX1(WTRSVBE(ielmc,NB,NZ))
    WTRVCX=AZMAX1(WTRVE(ielmc,NZ)*FWTBR)
    CPOOLD=(WTRVCX*WVSTBX-WTRSBX*WTRTTX)/WTPLTT
    XFRE(ielmc)=AZMAX1(XFRY*CPOOLD)
    WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)+XFRE(ielmc)
    WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)-XFRE(ielmc)
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
  real(r8) :: Rauto_pft,RCO2CM
! begin_execution
  associate(                             &
    CO2NetFix_pft      =>  plt_bgcr%CO2NetFix_pft    , &
    TCO2T     =>  plt_bgcr%TCO2T   , &
    Eco_AutoR_col      =>  plt_bgcr%Eco_AutoR_col    , &
    ECO_ER_col      =>  plt_bgcr%ECO_ER_col    , &
    Eco_GPP_col      =>  plt_bgcr%Eco_GPP_col    , &
    CARBN     =>  plt_bgcr%CARBN   , &
    TCO2A     =>  plt_bgcr%TCO2A   , &
    EPOOL     =>  plt_biom%EPOOL   , &
    CEPOLB    =>  plt_biom%CEPOLB  , &
    ZERO      => plt_site%ZERO     , &
    IGTYP     =>  plt_pheno%IGTYP  , &
    fTgrowCanP      =>  plt_pheno%fTgrowCanP   , &
    IWTYP     =>  plt_pheno%IWTYP  , &
    FDBKX     =>  plt_photo%FDBKX    &
  )
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
!
  IF(CEPOLB(ielmc,NB,NZ).GT.ZERO)THEN
    CNPG=AMIN1(CEPOLB(ielmn,NB,NZ)/(CEPOLB(ielmn,NB,NZ) &
      +CEPOLB(ielmc,NB,NZ)*CNKI),CEPOLB(ielmp,NB,NZ)/(CEPOLB(ielmp,NB,NZ) &
      +CEPOLB(ielmc,NB,NZ)*CPKI))
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
! fTgrowCanP=temperature function for canopy growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! FDBKX=termination feedback inhibition on C3 CO2
!
  RCO2C=AZMAX1(VMXC*EPOOL(ielmc,NB,NZ) &
    *fTgrowCanP(NZ))*CNPG*FDBKX(NB,NZ)*WFNG
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
  IF(RCO2Y.GT.0.0_r8.AND.(CNSHX.GT.0.0_r8.OR.CNLFX.GT.0.0_r8))THEN
    ZPOOLB=AZMAX1(EPOOL(ielmn,NB,NZ))
    PPOOLB=AZMAX1(EPOOL(ielmp,NB,NZ))
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
  ZADDB=AZMAX1(AMIN1(EPOOL(ielmn,NB,NZ),CGROS*(CNSHX+CNLFM+CNLFX*CNPG)))
  PADDB=AZMAX1(AMIN1(EPOOL(ielmp,NB,NZ) &
    ,CGROS*(CPSHX+CPLFM+CPLFX*CNPG)))
  CNRDA=AZMAX1(1.70*ZADDB-0.025*CH2O)
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! Rauto_pft=total C respiration
! RMNCS=maintenance respiration
! RCO2C=respiration from non-structural C
! RCO2G=growth respiration limited by N,P
! SNCR=excess maintenance respiration
! CNRDA=respiration for N assimilation
! CARBN=total PFT CO2 fixation
! CO2F=total CO2 fixation
! TCO2T,TCO2A=total,above-ground PFT respiration
! CO2NetFix_pft=PFT net CO2 fixation
! Eco_GPP_col=ecosystem GPP
! ECO_ER_col=ecosystem respiration
! Eco_AutoR_col=total autotrophic respiration
!
  Rauto_pft=AMIN1(RMNCS,RCO2C)+RCO2G+SNCR+CNRDA
  CARBN(NZ)=CARBN(NZ)+CO2F
  TCO2T(NZ)=TCO2T(NZ)-Rauto_pft
  TCO2A(NZ)=TCO2A(NZ)-Rauto_pft
  CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)+CO2F-Rauto_pft
  Eco_GPP_col=Eco_GPP_col+CO2F
  ECO_ER_col=ECO_ER_col-Rauto_pft
  Eco_AutoR_col=Eco_AutoR_col-Rauto_pft

  end associate
  end subroutine ComputeRAutoAfEmergence

!------------------------------------------------------------------------------------------

  subroutine ComputeRAutoB4Emergence(NB,NZ,TFN6,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,&
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
  real(r8) :: Rauto_pft,RCO2CM
  real(r8) :: RCO2XM,RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: SNCRM
! begin_execution
  associate(                          &
    CEPOLB    =>  plt_biom%CEPOLB   , &
    EPOOL     =>  plt_biom%EPOOL    , &
    IGTYP     =>  plt_pheno%IGTYP   , &
    fTgrowRootP      =>  plt_pheno%fTgrowRootP    , &
    WFR       =>  plt_rbgc%WFR      , &
    IWTYP     =>  plt_pheno%IWTYP   , &
    CO2NetFix_pft      =>  plt_bgcr%CO2NetFix_pft     , &
    RCO2M     =>  plt_rbgc%RCO2M    , &
    RCO2N     =>  plt_rbgc%RCO2N    , &
    RCO2A     =>  plt_rbgc%RCO2A    , &
    ZERO      => plt_site%ZERO      , &
    NGTopRootLayer       =>  plt_morph%NGTopRootLayer     , &
    FDBKX     =>  plt_photo%FDBKX     &
  )
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  IF(CEPOLB(ielmc,NB,NZ).GT.ZERO)THEN
    CNPG=AMIN1(CEPOLB(ielmn,NB,NZ)/(CEPOLB(ielmn,NB,NZ)+CEPOLB(ielmc,NB,NZ)*CNKI), &
      CEPOLB(ielmp,NB,NZ)/(CEPOLB(ielmp,NB,NZ)+CEPOLB(ielmc,NB,NZ)*CPKI))
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
! fTgrowRootP=temperature function for root growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! FDBKX=termination feedback inhibition on C3 CO2
! WFR=constraint by O2 consumption on all root processes
!
  RCO2CM=AZMAX1(VMXC*EPOOL(ielmc,NB,NZ) &
    *fTgrowRootP(NGTopRootLayer(NZ),NZ))*CNPG*WFNG*FDBKX(NB,NZ)
  RCO2C=RCO2CM*WFR(1,NGTopRootLayer(NZ),NZ)
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
  RMNCS=AZMAX1(RMPLT*TFN6(NGTopRootLayer(NZ))*WTSHXN)
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
    ZPOOLB=AZMAX1(EPOOL(ielmn,NB,NZ))
    PPOOLB=AZMAX1(EPOOL(ielmp,NB,NZ))
    FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG),PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))

    IF(RCO2YM.GT.0.0_r8)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
    ELSE
      RCO2GM=0._r8
    ENDIF

    IF(RCO2Y.GT.0.0_r8)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(1,NGTopRootLayer(NZ),NZ))
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
! RCO2TM,Rauto_pft=total C respiration unltd,ltd by O2
! RMNCS=maintenance respiration
! RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
! SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
! CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
! RCO2A=total root respiration
! RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
  RCO2TM=RMNCS+RCO2GM+SNCRM+CNRDM
  Rauto_pft=RMNCS+RCO2G+SNCR+CNRDA
  RCO2M(1,NGTopRootLayer(NZ),NZ)=RCO2M(1,NGTopRootLayer(NZ),NZ)+RCO2TM
  RCO2N(1,NGTopRootLayer(NZ),NZ)=RCO2N(1,NGTopRootLayer(NZ),NZ)+Rauto_pft
  RCO2A(1,NGTopRootLayer(NZ),NZ)=RCO2A(1,NGTopRootLayer(NZ),NZ)-Rauto_pft
  CH2O=0._r8
  end associate
  end subroutine ComputeRAutoB4Emergence

end module PlantBranchMod
