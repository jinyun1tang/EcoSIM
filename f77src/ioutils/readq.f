
      SUBROUTINE readq(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE READS INPUT DATA FROM PLANT SPECIES
C     AND MANAGEMENT FILES IDENTIFIED IN 'ROUTQ'
C
      use data_kind_mod, only : r8 => SHR_KIND_R8
      implicit none
      integer, intent(in) :: NEX
      integer, intent(in) :: NA(1:NEX),ND(1:NEX)
      integer, intent(in) :: NT,NE,NAX,NDX,NTX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS

      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk17.h"

      real(r8), PARAMETER :: TWILGT=0.06976
      integer :: IDX,IMO,IYR,IDY,ICUT,IDYE,IDYG,IDYS,JCUT,LPY,M
      integer :: NX,NY,NZ,NB,N,NN
      real(r8) :: DY,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23
      real(r8) :: ECUT24,HCUT,PCUT,VRNLI
      real(r8) :: VRNXI

C     execution begins here

      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      DO 9985 NZ=1,NP(NY,NX)
C
C     OPEN PFT(11), PFT MANAGEMENT(12) FILES FROM
C     FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
C
C     PREFIX=path for files in current or higher level directory
C     DATAP=PFT file name
C
      IF(DATAP(NZ,NY,NX).NE.'NO')THEN
C     WRITE(*,2233)'READQ',NX,NY,NZ,IETYP(NY,NX),DATAP(NZ,NY,NX)
2233  FORMAT(A8,4I4,2A16)
      OPEN(11,FILE=TRIM(PREFIX)//DATAP(NZ,NY,NX),STATUS='OLD')
      if(lverb)then
      write(*,*)'read parameters for pft ',NZ,' from ',DATAP(NZ,NY,NX)
      endif
      ENDIF
      IF(DATAM(NZ,NY,NX).NE.'NO')THEN
      OPEN(12,FILE=TRIM(PREFIX)//DATAM(NZ,NY,NX),STATUS='OLD')
      ENDIF
C
C     READ INPUTS FOR EACH PLANT SPECIES
C
      IF(DATAP(NZ,NY,NX).NE.'NO')THEN
C
C     PLANT FUNCTIONAL TYPE
C
C     ICTYP=photosynthesis type:3=C3,4=C4
C     IGTYP=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
C     ISTYP=growth habit:0=annual,1=perennial
C     IDTYP=growth habit:0=determinate,1=indetermimate
C     INTYP=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
C     4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
C     IPTYP=photoperiod type:0=day neutral,1=short day,2=long day
C     IBTYP=turnover:if IGTYP=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
C                   :if IGTYP=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
C     IRTYP=storage organ:0=above ground,1=below ground
C     MY=mycorrhizal:1=no,2=yes
C     ZTYPI=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
C     3=warm temperate,4=subtropical,5=tropical
C
      READ(11,*)ICTYP(NZ,NY,NX),IGTYP(NZ,NY,NX),ISTYP(NZ,NY,NX)
     2,IDTYP(NZ,NY,NX),INTYP(NZ,NY,NX),IWTYP(NZ,NY,NX)
     3,IPTYP(NZ,NY,NX),IBTYP(NZ,NY,NX),IRTYP(NZ,NY,NX),MY(NZ,NY,NX)
     4,ZTYPI(NZ,NY,NX)
      if(lverb)then
      write(*,*)'PLANT FUNCTIONAL TYPE'
      write(*,*)'photosynthesis type, 3=C3, 4=C4: ICTYP ',
     2ICTYP(NZ,NY,NX)
      write(*,*)'root profile, 0=shallow (eg bryophytes),'//
     2' 1=intermediate(eg herbs),2=deep (eg trees) : IGTYP ',
     3IGTYP(NZ,NY,NX)
      write(*,*)'ISTYP=growth habit, 0=annual,1=perennial: ISTYP',
     2ISTYP(NZ,NY,NX)
      write(*,*)'growth habit, 0=determinate,1=indetermimate: IDTYP ',
     2IDTYP(NZ,NY,NX)
      write(*,*)'N2 fixation, 1,2,3=rapid to slow root symbiosis '//
     2'(e.g.legumes), 4,5,6=rapid to slow canopy symbiosis (e.g.'//
     3'cyanobacteria): INTYP ',INTYP(NZ,NY,NX)
      write(*,*)'phenology type, 0=evergreen,1=cold deciduous,'//
     2'2=drought deciduous,3=1+2: IWTYP',IWTYP(NZ,NY,NX)
      write(*,*)'photoperiod type, 0=day neutral,1=short day'//
     2',2=long day: IPTYP ',IPTYP(NZ,NY,NX)
      write(*,*)'plant biom turnover, if IGTYP=0 or 1:all '//
     2'above-ground:0,1=rapid(deciduous),2=very slow(evergreen)'//
     3',3=slow(semi-deciduous): IBTYP',IBTYP(NZ,NY,NX)
      write(*,*)'storage organ:0=above ground,1=below '//
     2'ground: IRTYP',IRTYP(NZ,NY,NX)
      write(*,*)'mycorrhizal, 1=no,2=yes: MY',MY(NZ,NY,NX)
      write(*,*)'thermal adaptation zone:1=arctic,boreal,'//
     2'2=cool temperate: ZTYPI',ZTYPI(NZ,NY,NX)
      endif
C
C     PHOTOSYNTHETIC PROPERTIES
C
C     VCMX,VOMX=specific rubisco carboxylase,oxygenase activity (umol C,O g-1 s-1)
C     VCMX4=specific PEP carboxylase activity (umol g-1 s-1)
C     XKCO2,XKO2,XKCO24=Km for VCMX,VOMX,VCMX4 (uM)
C     RUBP,PEPC=fraction of leaf protein in rubisco, PEP carboxylase
C     ETMX=specific chlorophyll activity (umol e- g-1 s-1)
C     CHL=fraction of leaf protein in mesophyll(C3), bundle sheath(C4) chlorophyll
C     CHL4=fraction of leaf protein in mesophyll chlorophyll(C4)
C     FCO2=intercellular:atmospheric CO2 concentration ratio
C
      READ(11,*)VCMX(NZ,NY,NX),VOMX(NZ,NY,NX),VCMX4(NZ,NY,NX)
     2,XKCO2(NZ,NY,NX),XKO2(NZ,NY,NX),XKCO24(NZ,NY,NX)
     3,RUBP(NZ,NY,NX),PEPC(NZ,NY,NX),ETMX(NZ,NY,NX),CHL(NZ,NY,NX)
     4,CHL4(NZ,NY,NX),FCO2(NZ,NY,NX)

      if(lverb)then
      write(*,*)'PHOTOSYNTHETIC PROPERTIES'
      write(*,*)'specific rubisco carboxylase (umol C g-1 s-1): VCMX'
     2,VCMX(NZ,NY,NX)
      write(*,*)'specific rubisco oxygenase (umol O g-1 s-1): VOMX'
     2,VOMX(NZ,NY,NX)
      write(*,*)'specific PEP carboxylase activity (umol g-1 s-1): '//
     2' VCMX4',VCMX4(NZ,NY,NX)
      write(*,*)'Km for VCMX (uM): XKCO2',XKCO2(NZ,NY,NX)
      write(*,*)'Km for VOMX (uM): XKO2',XKO2(NZ,NY,NX)
      write(*,*)'KM for VCMX4 (uM): XKCO24',XKCO24(NZ,NY,NX)
      write(*,*)'fraction of leaf protein in rubisco: RUBP'
     2,RUBP(NZ,NY,NX)
      write(*,*)'fraction of leaf protein in PEP carboxylase: PEPC'
     2,PEPC(NZ,NY,NX)
      write(*,*)'specific chlorophyll activity (umol e- g-1 s-1): ETMX'
     2,ETMX(NZ,NY,NX)
      write(*,*)'fraction of leaf protein in mesophyll(C3),'
     2//' bundle sheath(C4) chlorophyll: CHL',CHL(NZ,NY,NX)
      write(*,*)'fraction of leaf protein in mesophyll chlorophyll(C4)'
     2,CHL4(NZ,NY,NX)
      write(*,*)'intercellular:atmospheric CO2 concentration '//
     2'ratio:FCO2 ',FCO2(NZ,NY,NX)
      endif
C
C     OPTICAL PROPERTIES
C
C     ALBP,ALBP,TAUR,TAUP=leaf SW,PAR albedo,transmission
C
      READ(11,*)ALBR(NZ,NY,NX),ALBP(NZ,NY,NX),TAUR(NZ,NY,NX)
     2,TAUP(NZ,NY,NX)
      if(lverb)then
      write(*,*)'OPTICAL PROPERTIES'
      write(*,*)'leaf SW albedo: ALBR',ALBR(NZ,NY,NX)
      write(*,*)'leaf PAR albedo: ALBP',ALBP(NZ,NY,NX)
      write(*,*)'leaf SW transmission: TAUR',TAUR(NZ,NY,NX)
      write(*,*)'leaf PAR transmission: TAUP',TAUP(NZ,NY,NX)
      endif
C
C     PHENOLOGICAL PROPERTIES
C
C     XRNI,XRLA=rate of node initiation,leaf appearance at 25oC (h-1)
C     CTC=chilling temperature for CO2 fixation, seed loss (oC)
C     VRNLI,VRNXI=hour requirement for spring leafout,autumn leafoff
C     WDLF=leaf length:width ratio
C     PB=nonstructural C concentration needed for branching
C     GROUPX,XTLI=node number required for floral initiation,at planting
C     XDL=critical photoperiod (h):<0=maximum daylength from site file
C     XPPD=photoperiod sensitivity (node h-1)
C
      READ(11,*)XRNI(NZ,NY,NX),XRLA(NZ,NY,NX),CTC(NZ,NY,NX)
     2,VRNLI,VRNXI,WDLF(NZ,NY,NX),PB(NZ,NY,NX)
      READ(11,*)GROUPX(NZ,NY,NX),XTLI(NZ,NY,NX),XDL(NZ,NY,NX)
     2,XPPD(NZ,NY,NX)
      if(lverb)then
      write(*,*)'PHENOLOGICAL PROPERTIES'
      write(*,*)'rate of node initiation at 25oC (h-1): XRNI'
     2,XRNI(NZ,NY,NX)
      write(*,*)'rate of leaf appearance at 25oC (h-1): XRLA'
     2,XRLA(NZ,NY,NX)
      write(*,*)'chilling temperature for CO2 fixation, '//
     2'seed loss (oC): CTC',CTC(NZ,NY,NX)
      write(*,*)'hour requirement for spring leafout: VRNLI',VRNLI
      write(*,*)'hour requirement for autumn leafoff: VRNXI',VRNXI
      write(*,*)'leaf length:width ratio: WDLF',WDLF(NZ,NY,NX)
      write(*,*)'nonstructural C concentration needed for branching'//
     2':PB',PB(NZ,NY,NX)
      endif
C
C     MORPHOLOGICAL PROPERTIES
C
C     SLA1,SSL1,SNL1=growth in leaf area,petiole length,internode length vs mass
C     CLASS=fraction of leaf area in 0-22.5,45,67.5,90o inclination classes
C     CFI=initial clumping factor
C     ANGBR,ANGSH=stem,petiole angle from horizontal
C     STMX=maximum potential seed mumber from pre-anthesis stalk growth
C     SDMX=maximum seed number per STMX
C     GRMX=maximum seed size per SDMX (g)
C     GRDM=seed size at planting (g)
C     GFILL=grain filling rate at 25 oC (g seed-1 h-1)
C     WTSTDI=mass of dead standing biomass at planting
C
      READ(11,*)SLA1(NZ,NY,NX),SSL1(NZ,NY,NX),SNL1(NZ,NY,NX)
      READ(11,*)(CLASS(N,NZ,NY,NX),N=1,4),CFI(NZ,NY,NX),ANGBR(NZ,NY,NX)
     2,ANGSH(NZ,NY,NX)
      READ(11,*)STMX(NZ,NY,NX),SDMX(NZ,NY,NX),GRMX(NZ,NY,NX)
     2,GRDM(NZ,NY,NX),GFILL(NZ,NY,NX),WTSTDI(NZ,NY,NX)
      if(lverb)then
      write(*,*)'MORPHOLOGICAL PROPERTIES'
      write(*,*)'growth in leaf area vs mass: SLA1 ',SLA1(NZ,NY,NX)
      write(*,*)'growth in petiole length vs mass: SSL1 ',SSL1(NZ,NY,NX)
      write(*,*)'growth in internode length vs mass: SNL1'
     2,SNL1(NZ,NY,NX)
      write(*,*)'fraction of leaf area in 0-22.5,45,67.5,90o '//
     2'inclination classes: CLASS',(CLASS(N,NZ,NY,NX),N=1,4)
      write(*,*)'initial clumping factor: CFI',CFI(NZ,NY,NX)
      write(*,*)'stem angle from horizontal: ANGBR',ANGBR(NZ,NY,NX)
      write(*,*)'petiole angle from horizontal: ANGSH',ANGSH(NZ,NY,NX)
      write(*,*)'maximum potential seed mumber from '//
     2'pre-anthesis stalk growth: STMX',STMX(NZ,NY,NX)
      write(*,*)'maximum seed number per STMX: SDMX',SDMX(NZ,NY,NX)
      write(*,*)'maximum seed size per SDMX (g): GRMX',GRMX(NZ,NY,NX)
      write(*,*)'seed size at planting (g): GRDM',GRDM(NZ,NY,NX)
      write(*,*)'grain filling rate at 25 oC (g seed-1 h-1): GFILL'
     2,GFILL(NZ,NY,NX)
      write(*,*)'mass of dead standing biomass at planting: WTSTDI'
     2,WTSTDI(NZ,NY,NX)
      endif
C
C     ROOT CHARACTERISTICS
C
C     RRAD1M,RRAD2M=radius of primary,secondary roots
C     PORT=root porosity
C     PR=nonstructural C concentration needed for root branching
C     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
C     PTSHT=rate constant for equilibrating shoot-root nonstructural C concn
C     RTFQ=root branching frequency (m-1)
C
      READ(11,*)RRAD1M(1,NZ,NY,NX),RRAD2M(1,NZ,NY,NX),PORT(1,NZ,NY,NX)
     2,PR(NZ,NY,NX),RSRR(1,NZ,NY,NX),RSRA(1,NZ,NY,NX)
     3,PTSHT(NZ,NY,NX),RTFQ(NZ,NY,NX)
      if(lverb)then
      write(*,*)'ROOT CHARACTERISTICS'
      write(*,*)'radius of primary roots: RRAD1M',RRAD1M(1,NZ,NY,NX)
      write(*,*)'radius of secondary roots: RRAD2M',RRAD2M(1,NZ,NY,NX)
      write(*,*)'root porosity: PORT',PORT(1,NZ,NY,NX)
      write(*,*)'nonstructural C concentration needed for root'//
     2' branching: PR',PR(NZ,NY,NX)
      write(*,*)'radial root resistivity (m2 MPa-1 h-1): RSRR'
     2,RSRR(1,NZ,NY,NX)
      write(*,*)'axial root resistivity (m2 MPa-1 h-1): RSRA'
     2,RSRA(1,NZ,NY,NX)
      write(*,*)'rate constant for equilibrating shoot-root '//
     2'nonstructural C concn: PTSHT',PTSHT(NZ,NY,NX)
      write(*,*)'root branching frequency (m-1): RTFQ',RTFQ(NZ,NY,NX)
      endif
C
C     ROOT UPTAKE PARAMETERS
C
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake (g m-2 h-1),Km (uM), min concn (uM)
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake (g m-2 h-1),Km (uM), min concn (uM)
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake (g m-2 h-1),Km (uM), min concn (uM)
C
      READ(11,*)UPMXZH(1,NZ,NY,NX),UPKMZH(1,NZ,NY,NX),UPMNZH(1,NZ,NY,NX)
      READ(11,*)UPMXZO(1,NZ,NY,NX),UPKMZO(1,NZ,NY,NX),UPMNZO(1,NZ,NY,NX)
      READ(11,*)UPMXPO(1,NZ,NY,NX),UPKMPO(1,NZ,NY,NX),UPMNPO(1,NZ,NY,NX)
      if(lverb)then
      write(*,*)'ROOT UPTAKE PARAMETERS'
      write(*,*)'NH4 max uptake (g m-2 h-1): UPMXZH',UPMXZH(1,NZ,NY,NX)
      write(*,*)'NH4 uptake Km (uM): UPKMZH',UPKMZH(1,NZ,NY,NX)
      write(*,*)'NH4 uptake min concn (uM): UPMNZH',UPMNZH(1,NZ,NY,NX)
      write(*,*)'NO3 max uptake (g m-2 h-1): UPMXZO',UPMXZO(1,NZ,NY,NX)
      write(*,*)'NO3 uptake Km (uM): UPKMZO',UPKMZO(1,NZ,NY,NX)
      write(*,*)'NO3 uptake min concn (uM): UPMNZO',UPMNZO(1,NZ,NY,NX)
      write(*,*)'H2PO4 max uptake (g m-2 h-1): UPMXPO'
     2,UPMXPO(1,NZ,NY,NX)
      write(*,*)'H2PO4 uptake Km (uM): UPKMPO',UPKMPO(1,NZ,NY,NX)
      write(*,*)'H2PO4 uptake min conc (uM): UPMNPO',UPMNPO(1,NZ,NY,NX)
      endif
C
C     WATER RELATIONS
C
C     OSMO=leaf osmotic potential at zero leaf water potential (MPa)
C     RCS=shape parameter for stomatal resistance vs leaf turgor potential
C     RSMX=cuticular resistance (s m-1)
C
      READ(11,*)OSMO(NZ,NY,NX),RCS(NZ,NY,NX),RSMX(NZ,NY,NX)
      if(lverb)then
      write(*,*)'WATER RELATIONS'
      write(*,*)'leaf osmotic potential at zero leaf water '//
     2'potential (MPa): OSMO',OSMO(NZ,NY,NX)
      write(*,*)'shape parameter for stomatal resistance vs '//
     2'leaf turgor potential: RCS',RCS(NZ,NY,NX)
      write(*,*)'cuticular resistance (s m-1): RSMX',RSMX(NZ,NY,NX)
      endif
C
C     ORGAN GROWTH YIELDS
C
C     DM*=DM-C production vs nonstructural C consumption (g g-1)
C     *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
C     *EAR=ear,*GR=grain,*RT=root,*ND=bacteria in root nodule,canopy
C
      READ(11,*)DMLF(NZ,NY,NX),DMSHE(NZ,NY,NX),DMSTK(NZ,NY,NX)
     2,DMRSV(NZ,NY,NX),DMHSK(NZ,NY,NX),DMEAR(NZ,NY,NX)
     3,DMGR(NZ,NY,NX),DMRT(NZ,NY,NX),DMND(NZ,NY,NX)
C
C     ORGAN N AND P CONCENTRATIONS
C
C     CN*,CP*=N:C,P:C ratios in plant organs
C
      READ(11,*)CNLF(NZ,NY,NX),CNSHE(NZ,NY,NX),CNSTK(NZ,NY,NX)
     2,CNRSV(NZ,NY,NX),CNHSK(NZ,NY,NX),CNEAR(NZ,NY,NX)
     3,CNGR(NZ,NY,NX),CNRT(NZ,NY,NX),CNND(NZ,NY,NX)
      READ(11,*)CPLF(NZ,NY,NX),CPSHE(NZ,NY,NX),CPSTK(NZ,NY,NX)
     2,CPRSV(NZ,NY,NX),CPHSK(NZ,NY,NX),CPEAR(NZ,NY,NX)
     3,CPGR(NZ,NY,NX),CPRT(NZ,NY,NX),CPND(NZ,NY,NX)

      if(lverb)then
      write(*,*)'ORGAN GROWTH YIELDS'
      write(*,*)'leaf dry matter C production vs '//
     2'nonstructural C consumption (g g-1): DMLF',DMLF(NZ,NY,NX)
      write(*,*)'petiole dry matter C production vs '//
     2'nonstructural C consumption (g g-1): DMSHE',DMSHE(NZ,NY,NX)
      write(*,*)'stalk dry matter C production vs '//
     2'nonstructural C consumption (g g-1): DMSTK',DMSTK(NZ,NY,NX)
      write(*,*)'stalk reserve C production vs '//
     2'nonstructural C consumption (g g-1): DMRSV',DMRSV(NZ,NY,NX)
      write(*,*)'husk dry matter C production vs '//
     2'nonstructural Cconsumption (g g-1): DMHSK',DMHSK(NZ,NY,NX)
      write(*,*)'ear dry matter C production vs '//
     2'nonstructural Cconsumption (g g-1): DMEAR',DMEAR(NZ,NY,NX)
      write(*,*)'grain C production vs nonstructural C'//
     2' consumption (g g-1): DMGR',DMGR(NZ,NY,NX)
      write(*,*)'root dry matter C production vs nonstructural C'//
     2' consumption (g g-1): DMRT',DMRT(NZ,NY,NX)
      write(*,*)'nodule bacteria in root nodule,canopy dry matter'//
     2'C production vs nonstructural C consumption (g g-1): DMND'
     2,DMND(NZ,NY,NX)
      write(*,*)'ORGAN N AND P CONCENTRATIONS'
      write(*,*)'NC ratio in plant leaves: CNLF',CNLF(NZ,NY,NX)
      write(*,*)'NC ratio in plant petiole: CNSHE',CNSHE(NZ,NY,NX)
      write(*,*)'NC ratio in plant stalk: CNSTK',CNSTK(NZ,NY,NX)
      write(*,*)'NC ratio in plant stalk reserve: CNRSV',CNRSV(NZ,NY,NX)
      write(*,*)'NC ratio in plant husk: CNHSK',CNHSK(NZ,NY,NX)
      write(*,*)'NC ratio in plant ear: CNEAR',CNEAR(NZ,NY,NX)
      write(*,*)'NC ratio in plant grain: CNGR',CNGR(NZ,NY,NX)
      write(*,*)'NC ratio in plant root: CNRT',CNRT(NZ,NY,NX)
      write(*,*)'NC ratio in plant module: CNND',CNND(NZ,NY,NX)
      write(*,*)'PC ratio in plant leaves: CPLF',CPLF(NZ,NY,NX)
      write(*,*)'PC ratio in plant petiole: CPSHE',CPSHE(NZ,NY,NX)
      write(*,*)'PC ratio in plant stalk: CPSTK',CPSTK(NZ,NY,NX)
      write(*,*)'PC ratio in plant stalk reserve: CPRSV',CPRSV(NZ,NY,NX)
      write(*,*)'PC ratio in plant husk: CPHSK',CPHSK(NZ,NY,NX)
      write(*,*)'PC ratio in plant ear: CPEAR',CPEAR(NZ,NY,NX)
      write(*,*)'PC ratio in plant grain: CPGR',CPGR(NZ,NY,NX)
      write(*,*)'PC ratio in plant root: CPRT',CPRT(NZ,NY,NX)
      write(*,*)'PC ratio in plant module: CPND',CPND(NZ,NY,NX)
      endif
C
C     RE-CALCULATE PLANT INPUTS IN MODEL UNITS
C
C     ABSR,ABSP=leaf SW,PAR absorbtivity
C
      VCMX(NZ,NY,NX)=2.5*VCMX(NZ,NY,NX)
      VOMX(NZ,NY,NX)=2.5*VOMX(NZ,NY,NX)
      VCMX4(NZ,NY,NX)=2.5*VCMX4(NZ,NY,NX)
      ETMX(NZ,NY,NX)=2.5*ETMX(NZ,NY,NX)
      ABSR(NZ,NY,NX)=1.0-ALBR(NZ,NY,NX)-TAUR(NZ,NY,NX)
      ABSP(NZ,NY,NX)=1.0-ALBP(NZ,NY,NX)-TAUP(NZ,NY,NX)
      ALBR(NZ,NY,NX)=ALBR(NZ,NY,NX)/ABSR(NZ,NY,NX)
      ALBP(NZ,NY,NX)=ALBP(NZ,NY,NX)/ABSP(NZ,NY,NX)
      TAUR(NZ,NY,NX)=TAUR(NZ,NY,NX)/ABSR(NZ,NY,NX)
      TAUP(NZ,NY,NX)=TAUP(NZ,NY,NX)/ABSP(NZ,NY,NX)
      ANGBR(NZ,NY,NX)=SIN(ANGBR(NZ,NY,NX)/57.29578)
      ANGSH(NZ,NY,NX)=SIN(ANGSH(NZ,NY,NX)/57.29578)
      GROUPI(NZ,NY,NX)=GROUPX(NZ,NY,NX)
      IF(IBTYP(NZ,NY,NX).NE.0)THEN
      XRNI(NZ,NY,NX)=XRNI(NZ,NY,NX)/25.0
      XRLA(NZ,NY,NX)=XRLA(NZ,NY,NX)/25.0
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)/25.0
      XTLI(NZ,NY,NX)=XTLI(NZ,NY,NX)/25.0
      ENDIF
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)-XTLI(NZ,NY,NX)
      IF(XDL(NZ,NY,NX).LT.0.0)THEN
      XDL(NZ,NY,NX)=DYLM(NY,NX)
      ENDIF
      DO 5 NB=1,10
      IF(IWTYP(NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      VRNL(NB,NZ,NY,NX)=AMIN1(4380.0,VRNLI+144.0*ZTYPI(NZ,NY,NX)*(NB-1))
      VRNX(NB,NZ,NY,NX)=AMIN1(4380.0,VRNXI+144.0*ZTYPI(NZ,NY,NX)*(NB-1))
      ELSE
      VRNL(NB,NZ,NY,NX)=VRNLI
      VRNX(NB,NZ,NY,NX)=VRNXI
      ENDIF
5     CONTINUE
      CLOSE(11)
      ENDIF
C     WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,XDL(NZ,NY,NX)
1111  FORMAT(A20,2I8,E12.4)
C
C     READ INPUTS FOR PLANT MANAGEMENT
C
      DO 15 M=1,366
      IHVST(NZ,M,NY,NX)=-1
      JHVST(NZ,M,NY,NX)=0
      HVST(NZ,M,NY,NX)=1.0E+06
      THIN(NZ,M,NY,NX)=-1.0
      EHVST(1,1,NZ,M,NY,NX)=1.0
      EHVST(1,2,NZ,M,NY,NX)=1.0
      EHVST(1,3,NZ,M,NY,NX)=1.0
      EHVST(1,4,NZ,M,NY,NX)=1.0
      EHVST(2,1,NZ,M,NY,NX)=1.0
      EHVST(2,2,NZ,M,NY,NX)=1.0
      EHVST(2,3,NZ,M,NY,NX)=1.0
      EHVST(2,4,NZ,M,NY,NX)=1.0
15    CONTINUE
      if(lverb)then
      write(*,*)'read pft management file ',DATAM(NZ,NY,NX)
      endif
      IF(DATAM(NZ,NY,NX).NE.'NO')THEN
      N=0
      NN=0
      IDAY0(NZ,NY,NX)=-1E+06
      IYR0(NZ,NY,NX)=-1E+06
      IDAYH(NZ,NY,NX)=1E+06
      IYRH(NZ,NY,NX)=1E+06
10    IF(N.EQ.0)THEN
C
C     PLANTING
C
C     DY=planting date DDMMYYYY
C     PPI=initial planting density (m-2)
C     SDPTHI=seeding depth (m)
C
      READ(12,*,END=540)DY,PPI(NZ,NY,NX),SDPTHI(NZ,NY,NX)
      SDPTHI(NZ,NY,NX)=AMAX1(1.0E-06,SDPTHI(NZ,NY,NX))
      if(lverb)then
      write(*,*)'planting date DDMMYYYY: DY',DY
      write(*,*)'initial planting density (m-2): PPI',PPI(NZ,NY,NX)
      write(*,*)'seeding depth (m): SDPTHI',SDPTHI(NZ,NY,NX)
      endif
      ELSE
C
C     HARVESTING
C
C     DY=if harvesting: harvesting date DDMMYYYY
C       =if grazing: first grazing date DDMMYYYY, followed by last grazing date in next line
C     ICUT=harvest type:0=none,1=grain,2=all above-ground,3=pruning
C     4=grazing,5=fire,6=insect
C     JCUT=terminate PFT:0=no,1=yes,2=yes,but reseed
C     HCUT=if harvesting: >0=cutting height,<0=fraction of LAI removed
C         =if grazing: grazer biomass (g FW m-2)
C     PCUT=if harvesting: thinning,fraction of population removed
C         =if grazing: grazer consumption rate (g DW g FW-1 d-1)
C     ECUT11,ECUT12,ECUT13,ECUT14=fraction of leaf,non-foliar,woody,standing dead
C     removed from PFT
C     ECUT21,ECUT22,ECUT23,ECUT124=fraction of leaf,non-foliar,woody,standing dead
C     removed from ecosystem
C
      READ(12,*,END=540)DY,ICUT,JCUT,HCUT,PCUT
     2,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23,ECUT24

      if(lverb)then
      write(*,*)'HARVESTING'
      write(*,*)'if harvesting: harvesting date DDMMYYYY '//
     2'if grazing: first grazing date DDMMYYYY, followed by '//
     3'last grazing date in next line: DY',DY
      write(*,*)'harvest type:0=none,1=grain,2=all '//
     2'above-ground,3=pruning: ICUT',ICUT
      write(*,*)'terminate PFT:0=no,1=yes,2=yes,but reseed: JCUT',JCUT
      write(*,*)'if harvesting: >0=cutting height,<0=fraction'//
     2' of LAI removed; if grazing: grazer consumption rate'//
     3' (g DW g FW-1 d-1)): PCUT',PCUT
      write(*,*)'fraction of leaf removed from PFT: ECUT11',ECUT11
      write(*,*)'fraction of non-foliar removed from PFT: ECUT12'
     2,ECUT12
      write(*,*)'fraction of wood removed from PFT: ECUT13',ECUT13
      write(*,*)'fraction of standing dead removed from PFT: ECUT14'
     2,ECUT14
      write(*,*)'fraction of leaf removed from ecosystem: ECUT21',ECUT21
      write(*,*)'fraction of non-foliar removed from ecosystem: ECUT22'
     2,ECUT22
      write(*,*)'fraction of woody removed from ecosystem: ECUT23'
     2,ECUT23
      write(*,*)'fraction of standing dead removed from ecosystem:'//
     2' ECUT24',ECUT24
      endif
      ENDIF
C
C     DERIVE PLANTING,HARVESTING DATES,YEARS
C
C     IDAY0,IYR0,IDAYH,IYRH=day,year of planting,harvesting
C
      IDX=INT(DY/1.0E+06)
      IMO=INT(DY/1.0E+04-IDX*1.0E+02)
      IYR=INT(DY-(IDX*1.0E+06+IMO*1.0E+04))
      LPY=0
      IF(MOD(IYR,4))520,510,520
510   IF(IMO.GT.2)LPY=1
520   IF(IMO.EQ.1)GO TO 525
      IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
      GO TO 527
525   IDY=IDX
527   IF(N.EQ.0)THEN
      IF(IDY.GT.0.AND.IYR.GT.0)THEN
C     IDY=IDY-0.5*(NTX-1)
      IDAY0(NZ,NY,NX)=IDY
C     IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      IYR=IYR+(NT-1)*NF+(NTX-1)*NFX-NTZX
C     ENDIF
      IYR0(NZ,NY,NX)=MIN(IYR,IYRC)
      IDAYX(NZ,NY,NX)=IDAY0(NZ,NY,NX)
      IYRX(NZ,NY,NX)=IYR0(NZ,NY,NX)
      PPZ(NZ,NY,NX)=PPI(NZ,NY,NX)
      ENDIF
      ELSE
      IF(IDY.GT.0.AND.JCUT.EQ.1)THEN
C     IDY=IDY+0.5*(NTX-1)
      IDAYH(NZ,NY,NX)=IDY
      IYR=IYR+(NT-1)*NF+(NTX-1)*NFX-NTZX
      IYRH(NZ,NY,NX)=MIN(IYR,IYRC)
      ENDIF
C
C     LOAD MODEL ARRAYS WITH PFT MANAGEMENT FOR SPECIFIED DATES
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground,3=pruning,4=grazing,5=fire,6=herbivory
C     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m2),IHVST=5:fire
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from PFT
C     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from ecosystem
C
      IHVST(NZ,IDY,NY,NX)=ICUT
      JHVST(NZ,IDY,NY,NX)=JCUT
      HVST(NZ,IDY,NY,NX)=HCUT
      THIN(NZ,IDY,NY,NX)=PCUT
      EHVST(1,1,NZ,IDY,NY,NX)=ECUT11
      EHVST(1,2,NZ,IDY,NY,NX)=ECUT12
      EHVST(1,3,NZ,IDY,NY,NX)=ECUT13
      EHVST(1,4,NZ,IDY,NY,NX)=ECUT14
      EHVST(2,1,NZ,IDY,NY,NX)=ECUT21
      EHVST(2,2,NZ,IDY,NY,NX)=ECUT22
      EHVST(2,3,NZ,IDY,NY,NX)=ECUT23
      EHVST(2,4,NZ,IDY,NY,NX)=ECUT24
C
C     IF ANIMAL OR INSECT GRAZING LOAD ALL DATES BETWEEN FIRST AND LAST
C
      IF(IHVST(NZ,IDY,NY,NX).EQ.4.OR.IHVST(NZ,IDY,NY,NX).EQ.6)THEN
      NN=NN+1
      IF(MOD(NN,2))570,560,570
560   IDYE=IDY
      DO 580 IDYG=IDYS+1,IDYE-1
      IHVST(NZ,IDYG,NY,NX)=ICUT
      JHVST(NZ,IDYG,NY,NX)=JCUT
      HVST(NZ,IDYG,NY,NX)=HCUT
      THIN(NZ,IDYG,NY,NX)=PCUT
      EHVST(1,1,NZ,IDYG,NY,NX)=ECUT11
      EHVST(1,2,NZ,IDYG,NY,NX)=ECUT12
      EHVST(1,3,NZ,IDYG,NY,NX)=ECUT13
      EHVST(1,4,NZ,IDYG,NY,NX)=ECUT14
      EHVST(2,1,NZ,IDYG,NY,NX)=ECUT21
      EHVST(2,2,NZ,IDYG,NY,NX)=ECUT22
      EHVST(2,3,NZ,IDYG,NY,NX)=ECUT23
      EHVST(2,4,NZ,IDYG,NY,NX)=ECUT24
580   CONTINUE
570   IDYS=IDY
      ENDIF
      ENDIF
      N=N+1
      GO TO 10
540   CLOSE(12)
      ELSE
      SDPTHI(NZ,NY,NX)=1.0E-06
      ENDIF
      IDAYY(NZ,NY,NX)=IDAYH(NZ,NY,NX)
      IYRY(NZ,NY,NX)=IYRH(NZ,NY,NX)
      IFLGC(NZ,NY,NX)=0
      IDTH(NZ,NY,NX)=0
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END
