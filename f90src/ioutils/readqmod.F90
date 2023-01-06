module readqmod
!!
! Description:
! code to read plant relevant files
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use fileUtil, only : open_safe, check_read
  use minimathmod, only : isLeap
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use CanopyDataType
  use PlantTraitDataType
  use PlantMngmtDataType
  use EcosimConst
  use RootDataType
  use EcoSIMHistMod
  use CanopyRadDataType
  use GridDataType
  implicit none
  private


  character(len=*), parameter :: mod_filename = __FILE__

  integer :: IDX,IMO,IYR,IDY,ICUT,IDYE,IDYG,IDYS,JCUT,LPY
  PUBLIC :: readq
  CONTAINS

  SUBROUTINE readq(NA,ND,NT,NE,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)
!!
! Description
! THIS SUBROUTINE READS INPUT DATA FROM PLANT SPECIES
!     AND MANAGEMENT FILES IDENTIFIED IN 'ROUTQ'
!
  implicit none
  integer, intent(in) :: NEX
  integer, intent(in) :: NA(1:NEX),ND(1:NEX)
  integer, intent(in) :: NT,NE,NTX,NF,NFX,NTZ &
    ,NTZX,NHW,NHE,NVN,NVS


  integer :: NX,NY,NZ

! begin_execution

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      D9985: DO NZ=1,NP(NY,NX)
!
!       OPEN PFT(11), PFT MANAGEMENT(12) FILES FROM
!       FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!
!       PREFIX=path for files in current or higher level directory
!       DATAP=PFT file name
        call ReadPlantProperties(NZ,NY,NX)

!
!  READ INPUTS FOR PLANT MANAGEMENT
        call ReadPlantManagement(NZ,NY,NX,NF,NFX,NT,NTX,NTZX)

        IFLGC(NZ,NY,NX)=0
        IDTH(NZ,NY,NX)=0
      ENDDO D9985
    ENDDO D9990
  ENDDO D9995
  RETURN
  END SUBROUTINE readq


!------------------------------------------------------------------------------------------

  subroutine ReadPlantManagement(NZ,NY,NX,NF,NFX,NT,NTX,NTZX)
  implicit none
  integer, intent(in) :: NZ,NY,NX,NF,NFX,NT,NTX,NTZX

  integer :: M,NN,N
  real(r8) :: DY,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23
  real(r8) :: ECUT24,HCUT,PCUT

! begin_execution
  IF(DATAM(NZ,NY,NX).NE.'NO')THEN
    call OPEN_safe(12,PREFIX,DATAM(NZ,NY,NX),'OLD',mod_filename,__LINE__)
  ENDIF

!
  D15: DO M=1,366
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
  ENDDO D15
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
    do while(.TRUE.)
    IF(N.EQ.0)THEN
!
!     PLANTING
!
!     DY=planting date DDMMYYYY
!     PPI=initial planting density (m-2)
!     SDPTHI=seeding depth (m)
!
      READ(12,*,END=540)DY,PPI(NZ,NY,NX),SDPTHI(NZ,NY,NX)
      SDPTHI(NZ,NY,NX)=AMAX1(ppmc,SDPTHI(NZ,NY,NX))
      if(lverb)then
        write(*,*)'planting date DDMMYYYY: DY',DY
        write(*,*)'initial planting density (m-2): PPI',PPI(NZ,NY,NX)
        write(*,*)'seeding depth (m): SDPTHI',SDPTHI(NZ,NY,NX)
      endif
    ELSE
!
!     HARVESTING
!
!     DY=if harvesting: harvesting date DDMMYYYY
!     =if grazing: first grazing date DDMMYYYY, followed by last grazing date in next line
!     ICUT=harvest type:0=none,1=grain,2=all above-ground,3=pruning
!     4=grazing,5=fire,6=insect
!     JCUT=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HCUT=if harvesting: >0=cutting height,<0=fraction of LAI removed
!         =if grazing: grazer biomass (g FW m-2)
!     PCUT=if harvesting: thinning,fraction of population removed
!         =if grazing: grazer consumption rate (g DW g FW-1 d-1)
!     ECUT11,ECUT12,ECUT13,ECUT14=fraction of leaf,non-foliar,woody,standing dead
!     removed from PFT
!     ECUT21,ECUT22,ECUT23,ECUT124=fraction of leaf,non-foliar,woody,standing dead
!     removed from ecosystem
!
      READ(12,*,END=540)DY,ICUT,JCUT,HCUT,PCUT &
        ,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23,ECUT24

      if(lverb)then
        write(*,*)'HARVESTING'
        write(*,*)'if harvesting: harvesting date DDMMYYYY '// &
          'if grazing: first grazing date DDMMYYYY, followed by '// &
          'last grazing date in next line: DY',DY
        write(*,*)'harvest type:0=none,1=grain,2=all '// &
          'above-ground,3=pruning: ICUT',ICUT
        write(*,*)'terminate PFT:0=no,1=yes,2=yes,but reseed: JCUT',JCUT
        write(*,*)'if harvesting: >0=cutting height,<0=fraction'// &
          ' of LAI removed; if grazing: grazer consumption rate'// &
          ' (g DW g FW-1 d-1)): PCUT',PCUT
        write(*,*)'fraction of leaf removed from PFT: ECUT11',ECUT11
        write(*,*)'fraction of non-foliar removed from PFT: ECUT12' &
          ,ECUT12
        write(*,*)'fraction of wood removed from PFT: ECUT13',ECUT13
        write(*,*)'fraction of standing dead removed from PFT: ECUT14' &
          ,ECUT14
        write(*,*)'fraction of leaf removed from ecosystem: ECUT21',ECUT21
        write(*,*)'fraction of non-foliar removed from ecosystem: ECUT22' &
          ,ECUT22
        write(*,*)'fraction of woody removed from ecosystem: ECUT23' &
          ,ECUT23
        write(*,*)'fraction of standing dead removed from ecosystem:'// &
          ' ECUT24',ECUT24
      endif
    ENDIF
!
!   DERIVE PLANTING,HARVESTING DATES,YEARS
!
!   IDAY0,IYR0,IDAYH,IYRH=day,year of planting,harvesting
!
    IDX=INT(DY/1.0E+06)
    IMO=INT(DY/1.0E+04-IDX*1.0E+02)
    IYR=INT(DY-(IDX*1.0E+06+IMO*1.0E+04))
    LPY=0

    if(isLeap(iyr) .and. IMO.GT.2)LPY=1
    IF(IMO.EQ.1)then
      IDY=IDX
    else
      IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
    endif
    IF(N.EQ.0)THEN
      IF(IDY.GT.0.AND.IYR.GT.0)THEN
!       IDY=IDY-0.5*(NTX-1)
        IDAY0(NZ,NY,NX)=IDY
!       IF(ISTYP(NZ,NY,NX).EQ.0)THEN
        IYR=IYR+(NT-1)*NF+(NTX-1)*NFX-NTZX
!       ENDIF
        IYR0(NZ,NY,NX)=MIN(IYR,IYRC)
        IDAYX(NZ,NY,NX)=IDAY0(NZ,NY,NX)
        IYRX(NZ,NY,NX)=IYR0(NZ,NY,NX)
        PPZ(NZ,NY,NX)=PPI(NZ,NY,NX)
      ENDIF
    ELSE
      IF(IDY.GT.0.AND.JCUT.EQ.1)THEN
!       IDY=IDY+0.5*(NTX-1)
        IDAYH(NZ,NY,NX)=IDY
        IYR=IYR+(NT-1)*NF+(NTX-1)*NFX-NTZX
        IYRH(NZ,NY,NX)=MIN(IYR,IYRC)
      ENDIF
!
!     LOAD MODEL ARRAYS WITH PFT MANAGEMENT FOR SPECIFIED DATES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground,3=pruning,4=grazing,5=fire,6=herbivory
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!
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
!
!     IF ANIMAL OR INSECT GRAZING LOAD ALL DATES BETWEEN FIRST AND LAST
!
      IF(IHVST(NZ,IDY,NY,NX).EQ.4.OR.IHVST(NZ,IDY,NY,NX).EQ.6)THEN
        NN=NN+1

        if(mod(nn,2)==0)then
          IDYE=IDY

          D580: DO IDYG=IDYS+1,IDYE-1
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
          ENDDO D580
        endif
!570     IDYS=IDY
        IDYS=IDY
      ENDIF
    ENDIF
    N=N+1
    enddo
540 CLOSE(12)
  ELSE
    SDPTHI(NZ,NY,NX)=ppmc
  ENDIF
  IDAYY(NZ,NY,NX)=IDAYH(NZ,NY,NX)
  IYRY(NZ,NY,NX)=IYRH(NZ,NY,NX)
  end subroutine ReadPlantManagement


!------------------------------------------------------------------------------------------

  subroutine ReadPlantProperties(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX

  integer :: N, NB, ierr
  real(r8) :: VRNXI,VRNLI
  character(len=200) :: tline

! begin_execution

  IF(DATAP(NZ,NY,NX).NE.'NO')THEN
!   WRITE(*,2233)'READQ',NX,NY,NZ,IETYP(NY,NX),DATAP(NZ,NY,NX)
!2233     FORMAT(A8,4I4,2A16)
    call OPEN_safe(11,PREFIX,DATAP(NZ,NY,NX),'OLD',mod_filename,__LINE__)
    if(lverb)then
      write(*,*)'read parameters for pft ',NZ,' from ',DATAP(NZ,NY,NX)
    endif
  ENDIF
!
! READ INPUTS FOR EACH PLANT SPECIES
!
  IF(DATAP(NZ,NY,NX).NE.'NO')THEN
!
!   PLANT FUNCTIONAL TYPE
!
!   ICTYP=photosynthesis type:3=C3,4=C4
!   IGTYP=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
!   ISTYP=growth habit:0=annual,1=perennial
!   IDTYP=growth habit:0=determinate,1=indetermimate
!   INTYP=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
!   4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
!   IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
!   IPTYP=photoperiod type:0=day neutral,1=short day,2=long day
!   IBTYP=turnover:if IGTYP=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
!                   :if IGTYP=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
!   IRTYP=storage organ:0=above ground,1=below ground
!   MY=mycorrhizal:1=no,2=yes
!   ZTYPI=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
!   3=warm temperate,4=subtropical,5=tropical
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)ICTYP(NZ,NY,NX),IGTYP(NZ,NY,NX),ISTYP(NZ,NY,NX) &
      ,IDTYP(NZ,NY,NX),INTYP(NZ,NY,NX),IWTYP(NZ,NY,NX) &
      ,IPTYP(NZ,NY,NX),IBTYP(NZ,NY,NX),IRTYP(NZ,NY,NX),MY(NZ,NY,NX) &
      ,ZTYPI(NZ,NY,NX)
    call check_read(ierr,11,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      call pft_display(NZ,NY,NX)
    endif
!
!   PHOTOSYNTHETIC PROPERTIES
!
!   VCMX,VOMX=specific rubisco carboxylase,oxygenase activity (umol C,O g-1 s-1)
!   VCMX4=specific PEP carboxylase activity (umol g-1 s-1)
!   XKCO2,XKO2,XKCO24=Km for VCMX,VOMX,VCMX4 (uM)
!   RUBP,PEPC=fraction of leaf protein in rubisco, PEP carboxylase
!   ETMX=specific chlorophyll activity (umol e- g-1 s-1)
!   CHL=fraction of leaf protein in mesophyll(C3), bundle sheath(C4) chlorophyll
!   CHL4=fraction of leaf protein in mesophyll chlorophyll(C4)
!   FCO2=intercellular:atmospheric CO2 concentration ratio
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)VCMX(NZ,NY,NX),VOMX(NZ,NY,NX),VCMX4(NZ,NY,NX) &
      ,XKCO2(NZ,NY,NX),XKO2(NZ,NY,NX),XKCO24(NZ,NY,NX) &
      ,RUBP(NZ,NY,NX),PEPC(NZ,NY,NX),ETMX(NZ,NY,NX),CHL(NZ,NY,NX) &
      ,CHL4(NZ,NY,NX),FCO2(NZ,NY,NX)
    call check_read(ierr,12,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      call photosyns_trait_disp(NZ,NY,NX)
    endif
!
!   OPTICAL PROPERTIES
!
!   ALBP,ALBP,TAUR,TAUP=leaf SW,PAR albedo,transmission
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)ALBR(NZ,NY,NX),ALBP(NZ,NY,NX),TAUR(NZ,NY,NX),TAUP(NZ,NY,NX)
    call check_read(ierr,4,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      write(*,*)'OPTICAL PROPERTIES'
      write(*,*)'leaf SW albedo: ALBR',ALBR(NZ,NY,NX)
      write(*,*)'leaf PAR albedo: ALBP',ALBP(NZ,NY,NX)
      write(*,*)'leaf SW transmission: TAUR',TAUR(NZ,NY,NX)
      write(*,*)'leaf PAR transmission: TAUP',TAUP(NZ,NY,NX)
    endif
!
!   PHENOLOGICAL PROPERTIES
!
!   XRNI,XRLA=rate of node initiation,leaf appearance at 25oC (h-1)
!   CTC=chilling temperature for CO2 fixation, seed loss (oC)
!   VRNLI,VRNXI=hour requirement for spring leafout,autumn leafoff
!   WDLF=leaf length:width ratio
!   PB=nonstructural C concentration needed for branching
!   GROUPX,XTLI=node number required for floral initiation,at planting
!   XDL=critical photoperiod (h):<0=maximum daylength from site file
!   XPPD=photoperiod sensitivity (node h-1)
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)XRNI(NZ,NY,NX),XRLA(NZ,NY,NX),CTC(NZ,NY,NX)  &
      ,VRNLI,VRNXI,WDLF(NZ,NY,NX),PB(NZ,NY,NX)
    call check_read(ierr,7,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)GROUPX(NZ,NY,NX),XTLI(NZ,NY,NX),XDL(NZ,NY,NX) &
      ,XPPD(NZ,NY,NX)
    call check_read(ierr,4,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      call Phenology_trait_disp(NZ,NY,NX,VRNLI,VRNXI)
    endif
!
!   MORPHOLOGICAL PROPERTIES
!
!   SLA1,SSL1,SNL1=growth in leaf area,petiole length,internode length vs mass
!   CLASS=fraction of leaf area in 0-22.5,45,67.5,90o inclination classes
!   CFI=initial clumping factor
!   ANGBR,ANGSH=stem,petiole angle from horizontal
!   STMX=maximum potential seed mumber from pre-anthesis stalk growth
!   SDMX=maximum seed number per STMX
!   GRMX=maximum seed size per SDMX (g)
!   GRDM=seed size at planting (g)
!   GFILL=grain filling rate at 25 oC (g seed-1 h-1)
!   WTSTDI=mass of dead standing biomass at planting
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)SLA1(NZ,NY,NX),SSL1(NZ,NY,NX),SNL1(NZ,NY,NX)
    call check_read(ierr,3,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)(CLASS(N,NZ,NY,NX),N=1,JLI),CFI(NZ,NY,NX),ANGBR(NZ,NY,NX) &
      ,ANGSH(NZ,NY,NX)
    call check_read(ierr,JLI+3,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)STMX(NZ,NY,NX),SDMX(NZ,NY,NX),GRMX(NZ,NY,NX) &
      ,GRDM(NZ,NY,NX),GFILL(NZ,NY,NX),WTSTDI(NZ,NY,NX)
    call check_read(ierr,6,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      call morphology_trait_disp(NZ,NY,NX)
    endif
!
!   ROOT CHARACTERISTICS
!
!   RRAD1M,RRAD2M=radius of primary,secondary roots
!   PORT=root porosity
!   PR=nonstructural C concentration needed for root branching
!   RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
!   PTSHT=rate constant for equilibrating shoot-root nonstructural C concn
!   RTFQ=root branching frequency (m-1)
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)RRAD1M(1,NZ,NY,NX),RRAD2M(1,NZ,NY,NX),PORT(1,NZ,NY,NX) &
      ,PR(NZ,NY,NX),RSRR(1,NZ,NY,NX),RSRA(1,NZ,NY,NX),PTSHT(NZ,NY,NX),RTFQ(NZ,NY,NX)
    call check_read(ierr,8,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      call Root_trait_disp(NZ,NY,NX)
    endif
!
!   ROOT UPTAKE PARAMETERS
!
!   UPMXZH,UPKMZH,UPMNZH=NH4 max uptake (g m-2 h-1),Km (uM), min concn (uM)
!   UPMXZO,UPKMZO,UPMNZO=NO3 max uptake (g m-2 h-1),Km (uM), min concn (uM)
!   UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake (g m-2 h-1),Km (uM), min concn (uM)
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)UPMXZH(1,NZ,NY,NX),UPKMZH(1,NZ,NY,NX),UPMNZH(1,NZ,NY,NX)
    call check_read(ierr,3,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)UPMXZO(1,NZ,NY,NX),UPKMZO(1,NZ,NY,NX),UPMNZO(1,NZ,NY,NX)
    call check_read(ierr,3,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)UPMXPO(1,NZ,NY,NX),UPKMPO(1,NZ,NY,NX),UPMNPO(1,NZ,NY,NX)
    call check_read(ierr,3,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      call Root_nutrient_trait_disp(NZ,NY,NX)
    endif
!
!   WATER RELATIONS
!
!   OSMO=leaf osmotic potential at zero leaf water potential (MPa)
!   RCS=shape parameter for stomatal resistance vs leaf turgor potential
!   RSMX=cuticular resistance (s m-1)
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)OSMO(NZ,NY,NX),RCS(NZ,NY,NX),RSMX(NZ,NY,NX)
    call check_read(ierr,3,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      write(*,*)'WATER RELATIONS'
      write(*,*)'leaf osmotic potential at zero leaf water '// &
         'potential (MPa): OSMO',OSMO(NZ,NY,NX)
      write(*,*)'shape parameter for stomatal resistance vs '// &
         'leaf turgor potential: RCS',RCS(NZ,NY,NX)
      write(*,*)'cuticular resistance (s m-1): RSMX',RSMX(NZ,NY,NX)
    endif
!
!   ORGAN GROWTH YIELDS
!
!   DM*=DM-C production vs nonstructural C consumption (g g-1)
!     *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
!     *EAR=ear,*GR=grain,*RT=root,*ND=bacteria in root nodule,canopy
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)DMLF(NZ,NY,NX),DMSHE(NZ,NY,NX),DMSTK(NZ,NY,NX) &
      ,DMRSV(NZ,NY,NX),DMHSK(NZ,NY,NX),DMEAR(NZ,NY,NX) &
      ,DMGR(NZ,NY,NX),DMRT(NZ,NY,NX),DMND(NZ,NY,NX)
    call check_read(ierr,9,DATAP(NZ,NY,NX),__LINE__,mod_filename)

!
!   ORGAN N AND P CONCENTRATIONS
!
!   CN*,CP*=N:C,P:C ratios in plant organs
!
    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)CNLF(NZ,NY,NX),CNSHE(NZ,NY,NX),CNSTK(NZ,NY,NX) &
      ,CNRSV(NZ,NY,NX),CNHSK(NZ,NY,NX),CNEAR(NZ,NY,NX) &
      ,CNGR(NZ,NY,NX),CNRT(NZ,NY,NX),CNND(NZ,NY,NX)
    call check_read(ierr,9,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    read(11,'(A)')tline
    READ(tline,*,iostat=ierr)CPLF(NZ,NY,NX),CPSHE(NZ,NY,NX),CPSTK(NZ,NY,NX) &
      ,CPRSV(NZ,NY,NX),CPHSK(NZ,NY,NX),CPEAR(NZ,NY,NX) &
      ,CPGR(NZ,NY,NX),CPRT(NZ,NY,NX),CPND(NZ,NY,NX)
    call check_read(ierr,9,DATAP(NZ,NY,NX),__LINE__,mod_filename)

    if(lverb)then
      write(*,*)'ORGAN GROWTH YIELDS'
      write(*,*)'leaf dry matter C production vs '// &
        'nonstructural C consumption (g g-1): DMLF',DMLF(NZ,NY,NX)
      write(*,*)'petiole dry matter C production vs '// &
        'nonstructural C consumption (g g-1): DMSHE',DMSHE(NZ,NY,NX)
      write(*,*)'stalk dry matter C production vs '// &
        'nonstructural C consumption (g g-1): DMSTK',DMSTK(NZ,NY,NX)
      write(*,*)'stalk reserve C production vs '// &
        'nonstructural C consumption (g g-1): DMRSV',DMRSV(NZ,NY,NX)
      write(*,*)'husk dry matter C production vs '// &
        'nonstructural Cconsumption (g g-1): DMHSK',DMHSK(NZ,NY,NX)
      write(*,*)'ear dry matter C production vs '// &
        'nonstructural Cconsumption (g g-1): DMEAR',DMEAR(NZ,NY,NX)
      write(*,*)'grain C production vs nonstructural C'// &
        ' consumption (g g-1): DMGR',DMGR(NZ,NY,NX)
      write(*,*)'root dry matter C production vs nonstructural C'// &
        ' consumption (g g-1): DMRT',DMRT(NZ,NY,NX)
      write(*,*)'nodule bacteria in root nodule,canopy dry matter'// &
        'C production vs nonstructural C consumption (g g-1): DMND' &
        ,DMND(NZ,NY,NX)
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
!
!   RE-CALCULATE PLANT INPUTS IN MODEL UNITS
!
!   ABSR,ABSP=leaf SW,PAR absorbtivity
!
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
    ANGBR(NZ,NY,NX)=SIN(ANGBR(NZ,NY,NX)/57.29578_r8)
    ANGSH(NZ,NY,NX)=SIN(ANGSH(NZ,NY,NX)/57.29578_r8)
    GROUPI(NZ,NY,NX)=GROUPX(NZ,NY,NX)
    IF(IBTYP(NZ,NY,NX).NE.0)THEN
      XRNI(NZ,NY,NX)=XRNI(NZ,NY,NX)/25.0_r8
      XRLA(NZ,NY,NX)=XRLA(NZ,NY,NX)/25.0_r8
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)/25.0_r8
      XTLI(NZ,NY,NX)=XTLI(NZ,NY,NX)/25.0_r8
    ENDIF
    GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)-XTLI(NZ,NY,NX)
    IF(XDL(NZ,NY,NX).LT.0.0_r8)THEN
      XDL(NZ,NY,NX)=DYLM(NY,NX)
    ENDIF
    D5: DO NB=1,JC
      IF(IWTYP(NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.iplt_annual)THEN
        VRNL(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNLI+144.0_r8*ZTYPI(NZ,NY,NX)*(NB-1))
        VRNX(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNXI+144.0_r8*ZTYPI(NZ,NY,NX)*(NB-1))
      ELSE
        VRNL(NB,NZ,NY,NX)=VRNLI
        VRNX(NB,NZ,NY,NX)=VRNXI
      ENDIF
    ENDDO D5
    CLOSE(11)
  ENDIF
! WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,XDL(NZ,NY,NX)
!1111    FORMAT(A20,2I8,E12.4)
  end subroutine ReadPlantProperties
!------------------------------------------------------------------------------------------
  subroutine pft_display(NZ,NY,NX)
  use abortutils , only : endrun
  implicit none
  integer, intent(in) :: NZ,NY,NX
!   ICTYP=photosynthesis type:3=C3,4=C4
!   IGTYP=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
!   ISTYP=growth habit:0=annual,1=perennial
!   IDTYP=growth habit:0=determinate,1=indetermimate
!   INTYP=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
!   4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
!   IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
!   IPTYP=photoperiod type:0=day neutral,1=short day,2=long day
!   IBTYP=turnover:if IGTYP=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
!                   :if IGTYP=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
!   IRTYP=storage organ:0=above ground,1=below ground
!   MY=mycorrhizal:1=no,2=yes
!   ZTYPI=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
!   3=warm temperate,4=subtropical,5=tropical

  write(*,*)'PLANT FUNCTIONAL TYPE for (NZ,NY,NX)=',NZ,NY,NX
  select CASE (ICTYP(NZ,NY,NX))
  case (3)
    write(*,*)'C3 photosynthesis'
  case (4)
    write(*,*)'C4 photosynthesis'
  case default
    write(*,*)'photosynthesis type not defined'
  end select

  select case(IGTYP(NZ,NY,NX))
  case (0)
    write(*,*)'shallow root profile, like bryophytes'
  case (1)
    write(8,*)'intermediate root profile, like herbs'
  case (2)
    write(*,*)'deep root profile, like trees'
  case default
    write(*,*)'root profile not defined'
  end select

  select case (ISTYP(NZ,NY,NX))
  case (0)
    write(*,*)'Annual plant'
  case (1)
    write(*,*)'perennial plant'
  case default
    write(*,*)'growth habit not defined'
  end select

  select case (IDTYP(NZ,NY,NX))
  case (0)
    write(*,*)'determinate growth habit'
  case (1)
    write(*,*)'indetermimate growth habit'
  case default
    write(*,*)'growth habit not defined'
  end select

  select case (INTYP(NZ,NY,NX))
! 1,2, 3, e.g. legumes
  case (1)
    write(*,*)'Rapid root N2 fixation symbiosis'
  case (2)
    write(*,*)'Intermediate root N2 fixation symbiosis'
  case (3)
    write(*,*)'Slow root N2 fixation symbiosis'
!4,5,6, e.g. cyanobacteria
  case (4)
    write(*,*)'Rapid canopy N2 fixation symbiosis'
  case (5)
    write(*,*)'Intermediate canopy N2 fixation symbiosis'
  case (6)
    write(*,*)'Slow canopy N2 fixation symbiosis'
  case default
    write(*,*)'No N2 fixation symbiosis defined'
  end select

  select case(IWTYP(NZ,NY,NX))
  case (0)
    write(*,*)'phenology type: evergreen'
  case (1)
    write(*,*)'phenology type: cold deciduous'
  case (2)
    write(*,*)'phenology type: drought deciduous'
  case (3)
    write(*,*)'phenology type: cold+drought deciduous'
  case default
    write(*,*)'phenology type not defined'
  end select

  select case(IPTYP(NZ,NY,NX))
  case (0)
    write(*,*)'day neutral photoperiod'
  case (1)
    write(*,*)'short day photoperiod'
  case (2)
    write(*,*)'long day photoperiod'
  case default
    write(*,*)'photoperiod not defined'
  end select

  select case(IBTYP(NZ,NY,NX))
  case (0, 1)
    write(*,*)'plant biome turnover rapid (deciduous)'
  case (2)
    write(*,*)'Plant biome turnover very slow (needleleaf evergreen)'
  case (3)
    write(*,*)'Plant biome turnover slow broadleaf evergreen'
  case (4)
    write(*,*)'Plant biome turnover semi-deciduous'
  case (5)
    write(*,*)'Plant biome turnover semi-evergreen'
  case default
    write(*,*)'Plant biome turnover not defined'
  end select

  select case(IRTYP(NZ,NY,NX))
  case (0)
    write(*,*)'Above ground storage organ'
  case (1)
    write(*,*)'Belowground storage organ'
  case default
    write(*,*)'Storage organ not defined'
  end select

  select case(MY(NZ,NY,NX))
  case (1)
    write(*,*)'No mycorrhizal'
  case (2)
    write(*,*)'Mycorrhizal'
  case default
    write(*,*)'Wrong option for mycorrhizae'
  end select

  select case(INT(ZTYPI(NZ,NY,NX)))
  case (1)
    write(*,*)'thermal adaptation zone: arctic, boreal'
  case (2)
    write(*,*)'thermal adaptation zone: cool temperate'
  case (3)
    write(*,*)'thermal adaptation zone: warm temperate'
  case (4)
    write(*,*)'thermal adaptation zone: subtropical'
  case (5)
    write(*,*)'thermal adaptation zone: tropical'
  case default
    write(*,*)'Not thermal adaptation zone defined'
  end select

  end subroutine pft_display

!------------------------------------------------------------------------------------------
  subroutine photosyns_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'PHOTOSYNTHETIC PROPERTIES'
  write(*,*)'Specific rubisco carboxylase (umol C g-1 s-1): ',VCMX(NZ,NY,NX)
  write(*,*)'Specific rubisco oxygenase (umol O g-1 s-1): ',VOMX(NZ,NY,NX)
  write(*,*)'Specific PEP carboxylase activity (umol g-1 s-1): '&
    ,VCMX4(NZ,NY,NX)
  write(*,*)'Km for VCMX (uM): ',XKCO2(NZ,NY,NX)
  write(*,*)'Km for VOMX (uM): ',XKO2(NZ,NY,NX)
  write(*,*)'KM for VCMX4 (uM): ',XKCO24(NZ,NY,NX)
  write(*,*)'Fraction of leaf protein in rubisco: ',RUBP(NZ,NY,NX)
  write(*,*)'Fraction of leaf protein in PEP carboxylase: ',PEPC(NZ,NY,NX)
  write(*,*)'Specific chlorophyll activity (umol e- g-1 s-1): ',ETMX(NZ,NY,NX)
  if(ICTYP(NZ,NY,NX).eq.ic3_photo)then
    write(*,*)'Fraction of leaf protein as chlorophyll in mesophyll (C3) ',CHL(NZ,NY,NX)
  elseif(ICTYP(NZ,NY,NX).eq.ic4_photo)then
    write(*,*)'Fraction of leaf protein as chlorophyll in bundle sheath(C4)',CHL(NZ,NY,NX)
    write(*,*)'fraction of leaf protein in mesophyll chlorophyll(C4)',CHL4(NZ,NY,NX)
  endif
  write(*,*)'intercellular:atmospheric CO2 concentration ratio:',FCO2(NZ,NY,NX)
  end subroutine photosyns_trait_disp

!------------------------------------------------------------------------------------------

  subroutine Phenology_trait_disp(NZ,NY,NX,VRNLI,VRNXI)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: VRNLI,VRNXI
  write(*,*)'PHENOLOGICAL PROPERTIES'
  write(*,*)'rate of node initiation at 25oC (h-1): XRNI',XRNI(NZ,NY,NX)
  write(*,*)'rate of leaf appearance at 25oC (h-1): XRLA',XRLA(NZ,NY,NX)
  write(*,*)'chilling temperature for CO2 fixation, '// &
    'seed loss (oC): CTC',CTC(NZ,NY,NX)
  write(*,*)'hour requirement for spring leafout: VRNLI',VRNLI
  write(*,*)'hour requirement for autumn leafoff: VRNXI',VRNXI
  write(*,*)'leaf length:width ratio: WDLF',WDLF(NZ,NY,NX)
  write(*,*)'nonstructural C concentration needed for branching'// &
    ':PB',PB(NZ,NY,NX)
  end subroutine Phenology_trait_disp
!------------------------------------------------------------------------------------------
  subroutine morphology_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer :: N

  write(*,*)'MORPHOLOGICAL PROPERTIES'
  write(*,*)'growth in leaf area vs mass: SLA1 ',SLA1(NZ,NY,NX)
  write(*,*)'growth in petiole length vs mass: SSL1 ',SSL1(NZ,NY,NX)
  write(*,*)'growth in internode length vs mass: SNL1',SNL1(NZ,NY,NX)
  write(*,*)'fraction of leaf area in 0-22.5,45,67.5,90o '// &
    'inclination classes: CLASS',(CLASS(N,NZ,NY,NX),N=1,JLI)
  write(*,*)'initial clumping factor: CFI',CFI(NZ,NY,NX)
  write(*,*)'stem angle from horizontal: ANGBR',ANGBR(NZ,NY,NX)
  write(*,*)'petiole angle from horizontal: ANGSH',ANGSH(NZ,NY,NX)
  write(*,*)'maximum potential seed mumber from '// &
    'pre-anthesis stalk growth: STMX',STMX(NZ,NY,NX)
  write(*,*)'maximum seed number per STMX: SDMX',SDMX(NZ,NY,NX)
  write(*,*)'maximum seed size per SDMX (g): GRMX',GRMX(NZ,NY,NX)
  write(*,*)'seed size at planting (g): GRDM',GRDM(NZ,NY,NX)
  write(*,*)'grain filling rate at 25 oC (g seed-1 h-1): GFILL' &
    ,GFILL(NZ,NY,NX)
  write(*,*)'mass of dead standing biomass at planting: WTSTDI' &
    ,WTSTDI(NZ,NY,NX)
  end subroutine morphology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'ROOT CHARACTERISTICS'
  write(*,*)'radius of primary roots: RRAD1M',RRAD1M(1,NZ,NY,NX)
  write(*,*)'radius of secondary roots: RRAD2M',RRAD2M(1,NZ,NY,NX)
  write(*,*)'primary/fine root porosity: PORT',PORT(1,NZ,NY,NX)
  write(*,*)'nonstructural C concentration needed for root'// &
    ' branching: PR',PR(NZ,NY,NX)
  write(*,*)'radial root resistivity for water uptake (m2 MPa-1 h-1): RSRR',RSRR(1,NZ,NY,NX)
  write(*,*)'axial root resistivity for water uptake (m2 MPa-1 h-1): RSRA',RSRA(1,NZ,NY,NX)
  write(*,*)'rate constant for equilibrating shoot-root '// &
    'nonstructural C concn: PTSHT',PTSHT(NZ,NY,NX)
  write(*,*)'root branching frequency (m-1): RTFQ',RTFQ(NZ,NY,NX)
  end subroutine Root_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_nutrient_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'ROOT UPTAKE PARAMETERS'
  write(*,*)'NH4 max uptake (g m-2 h-1): UPMXZH',UPMXZH(1,NZ,NY,NX)
  write(*,*)'NH4 uptake Km (uM): UPKMZH',UPKMZH(1,NZ,NY,NX)
  write(*,*)'NH4 uptake min concn (uM): UPMNZH',UPMNZH(1,NZ,NY,NX)
  write(*,*)'NO3 max uptake (g m-2 h-1): UPMXZO',UPMXZO(1,NZ,NY,NX)
  write(*,*)'NO3 uptake Km (uM): UPKMZO',UPKMZO(1,NZ,NY,NX)
  write(*,*)'NO3 uptake min concn (uM): UPMNZO',UPMNZO(1,NZ,NY,NX)
  write(*,*)'H2PO4 max uptake (g m-2 h-1): UPMXPO',UPMXPO(1,NZ,NY,NX)
  write(*,*)'H2PO4 uptake Km (uM): UPKMPO',UPKMPO(1,NZ,NY,NX)
  write(*,*)'H2PO4 uptake min conc (uM): UPMNPO',UPMNPO(1,NZ,NY,NX)
  end subroutine Root_nutrient_trait_disp

end module readqmod
