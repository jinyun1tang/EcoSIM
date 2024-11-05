module readiMod
!!
! code to read site, topographic data
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils   , only : endrun
  use ncdio_pio
  use EcoSIMCtrlMod
  use fileUtil     , only : open_safe, check_read
  use minimathmod  , only : isclose, AZMAX1
  use MiniFuncMod  , only : GetDayLength
  use EcoSIMConfig , only : column_mode
  use EcoSiMParDataMod, only : micpar
  use SOMDataType
  use CanopyRadDataType
  use EcosimConst
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use ChemTranspDataType
  use LandSurfDataType
  use ClimForcDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SurfLitterDataType
  use SnowDataType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use GridDataType
  use SoilHydroParaMod, only : ComputeSoilHydroPars
  use SoilPhysParaMod, only : SetDeepSoil
  implicit none
  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8) :: datav(40)
  CHARACTER(len=16) :: OUTW,OUTI,OUTT,OUTN,OUTF
  CHARACTER(len=4) :: CHARY
  CHARACTER(len=1) :: TTYPE,CTYPE,IVAR(20),VAR(50),TYP(50)
  character(len=3), parameter :: model_status(0:1)=(/'off','on '/)

  integer :: ll
  real(r8) :: DAT(50),DATK(50)
  real(r8) :: ALATG,ATCAG,AZI,ASPX,CO2EIG,CH4EG,WTBLDepz_nat,WTBLDepz_tile
  real(r8) :: DTBLGG,DEC,initSnowDepth,OXYEG,RCHQNG,RCHQEG
  real(r8) :: RCHQSG,RCHQWG,RCHGNUG,RCHGEUG,RCHGSUG,RCHGWUG
  real(r8) :: RCHGNTG,RCHGETG,RCHGSTG,RCHGWTG,RCHGDG
  real(r8) :: SL0,Z2GEG,Z2OEG,ZNH3EG,SLX,SL1,SL2

  integer :: L,NH1,NH2,NV1,NV2,NL1
  integer :: NL2
  type(file_desc_t) :: grid_nfid

  public :: readi,erosion_model_status
  public :: GridConectionMode
  contains

  SUBROUTINE readi(NHW,NHE,NVN,NVS)
!!
! Description:
! THIS SUBROUTINE READS ALL SOIL AND TOPOGRAPHIC INPUT FILES
!
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: jj,NX,NY
  integer :: ierr
  character(len=200) :: tline
!
!
! OPEN OUTPUT LOGFILES,AND SITE,TOPOGRAPHY FILES FROM
! FILE NAMES IN DATA ARRAYS LOADED IN 'MAIN'
!
  OPEN(18,FILE=trim(outdir)//'logfile1',STATUS='UNKNOWN')
  OPEN(19,FILE=trim(outdir)//'logfile2',STATUS='UNKNOWN')
  OPEN(20,FILE=trim(outdir)//'logfile3',STATUS='UNKNOWN')

  WRITE(18,5000)' 08 JUL 2021'
5000  FORMAT(A16)

  call ncd_pio_openfile(grid_nfid, grid_file_in, ncd_nowrite)

  call readsiteNC(NHW,NHE,NVN,NVS)

  call readTopoNC

  if(.not.column_mode)then
    DO  NX=NHW,NHE
      NL(NVS+1,NX)=NL(NVS,NX)
    ENDDO
    DO  NY=NVN,NVS
      NL(NY,NHE+1)=NL(NY,NHE)
    ENDDO
    NL(NVS+1,NHE+1)=NL(NVS,NHE)
  endif
  IOLD=0
  call ncd_pio_closefile(grid_nfid)


  END SUBROUTINE readi

!------------------------------------------------------------------------------------------

  function erosion_model_status(flag)result(status)

  implicit none
  integer, intent(in) :: flag
  character(len=40) :: status
  !0 means allowing  freeze- thaw to change elevation
  !1 means allowing freeze-thaw plus erosion to change elevation
  !2 means allowing freeze-thaw plus SOC accumulation to change elevation
  !3 means allowing freeze-thaw plus SOC accumulation, plus erosion to change elevation
  !-1 means no change in elevation.
  select case (flag)
  case (-1)
    status='no elv change'
  case (0)
    status='freeze-thaw elv change'
  case (1)
    status='freeze-thaw+erosion elv change'
  case (2)
    status='freeze-thaw+SOC accum elv change'
  case (3)
    status='freeze-thaw+SOC+erosion elv change'
  case default
    status=''
    call endrun('wrong erosion model option in '//trim(mod_filename)//' at line',__LINE__)
  end select
  end function erosion_model_status

!------------------------------------------------------------------------------------------

  function GridConectionMode(NCNG)result(status)
  implicit none
  integer, intent(in) :: NCNG
  character(len=40) :: status

  select case(NCNG)
  case (1)
    status='3D grid'   !3D
  case (2)
    status='2D Along the north-south direction' !2d, no x direction
  case (3)
    status='1D vertical column'           !1d, vertical only
  case default
    status=''
    call endrun('wrong option for NCNG in '//trim(mod_filename)//' at line',__LINE__)
  end select

  end function GridConectionMode

!------------------------------------------------------------------------------------------
  function WaterTableStatus(iWaterTabelMode)result(status)

  implicit none
  integer, intent(in) :: iWaterTabelMode

  character(len=64) :: status

  select case(iWaterTabelMode)
  case (0)
    status='No external water table'
  case (1)
    status='Natural stationary external water table'
  case (2)
    status='Natural mobile external water table'
  case (3)
    status='Artificial stationary extneral water table'
  case (4)
    status='Artificial mobile external water table'
  case default
    call endrun('wrong option for iWaterTabelMode in '//trim(mod_filename)//' at line',__LINE__)
  end select
  end function WaterTableStatus

!------------------------------------------------------------------------------------------

  subroutine readsiteNC(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: DHI(JX),DVI(JY)
  character(len=200) :: tline
  integer :: XI
  integer :: NY,NX
  integer :: IETYPG
  integer :: ierr,jj,loc
  integer :: iWaterTabelMode
!
! READ SITE DATA
!
! ALATG,ALTIG,ATCAG=latitude,altitude,MAT(oC)
! iWaterTabelMode=water table flag
! :0=none
! :1,2=natural stationary,mobile
! :3,4=artificial stationary,mobile
! OXYEG,Z2GEG,CO2EIG,CH4EG,Z2OEG,ZNH3EG=atm O2,N2,CO2,CH4,N2O,NH3 (ppm)
! IETYPG,iErosionMode=Koppen climate zone,erosion options
! NCNG=1:lateral connections between grid cells,3:no connections
! WTBLDepz_nat,WTBLDepz_tile=depth of natural,artificial (tile) water table (iWaterTabelMode)
! DTBLGG=slope of natural water table relative to landscape surface
! RCHQNG,RCHQEG,RCHQSG,RCHQWG=boundary condns for N,E,S,W surface runoff
! RCHGNUG,RCHGEUG,RCHGSUG,RCHGWUG=bound condns for N,E,S,W subsurf flow
! RCHGNTG,RCHGETG,RCHGSTG,RCHGWTG=N,E,S,W distance to water table (m)
! RCHGDG=lower boundary conditions for water flow
! DHI=width of each W-E landscape column
! DVI=width of each N-S landscape row

! assuming only one grid, at the moment
  loc=1
  call ncd_getvar(grid_nfid,'ALATG',loc,ALATG)
  call ncd_getvar(grid_nfid,'ALTIG',loc,ALTIG)
  call ncd_getvar(grid_nfid,'ATCAG',loc,ATCAG)
  call ncd_getvar(grid_nfid,'IDTBLG',loc,iWaterTabelMode)

  call ncd_getvar(grid_nfid,'IETYPG',loc,IETYPG)
  call ncd_getvar(grid_nfid,'DTBLIG',loc,WTBLDepz_nat)
  call ncd_getvar(grid_nfid,'DTBLDIG',loc,WTBLDepz_tile)
  call ncd_getvar(grid_nfid,'DTBLGG',loc,DTBLGG)

  call ncd_getvar(grid_nfid,'RCHQNG',loc,RCHQNG)
  call ncd_getvar(grid_nfid,'RCHQEG',loc,RCHQEG)
  call ncd_getvar(grid_nfid,'RCHQSG',loc,RCHQSG)
  call ncd_getvar(grid_nfid,'RCHQWG',loc,RCHQWG)
  call ncd_getvar(grid_nfid,'RCHGNUG',loc,RCHGNUG)
  call ncd_getvar(grid_nfid,'RCHGEUG',loc,RCHGEUG)
  call ncd_getvar(grid_nfid,'RCHGSUG',loc,RCHGSUG)
  call ncd_getvar(grid_nfid,'RCHGWUG',loc,RCHGWUG)
  call ncd_getvar(grid_nfid,'RCHGNTG',loc,RCHGNTG)
  call ncd_getvar(grid_nfid,'RCHGETG',loc,RCHGETG)
  call ncd_getvar(grid_nfid,'RCHGSTG',loc,RCHGSTG)
  call ncd_getvar(grid_nfid,'RCHGWTG',loc,RCHGWTG)
  call ncd_getvar(grid_nfid,'RCHGDG',loc,RCHGDG)

  call ncd_getvar(grid_nfid,'DHI',loc,DHI(1:NHE))
  call ncd_getvar(grid_nfid,'DVI',loc,DVI(1:NVS))

  if(lverb)then
    write(*,*)'read site data file: ',DATA1(1)
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'Latitude (o): ALATG',ALATG
    write(*,*)'Altitude (m): ALTIG',ALTIG
    write(*,*)'Mean annual temperaure (oC): ATCAG',ATCAG
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'atmospheric NH3 (ppm): ZNH3EG',ZNH3EG
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'Koppen climate zone: IETYPG',IETYPG

    write(*,*)'depth of natural water table: DTBLIG',WTBLDepz_nat
    write(*,*)'depth of artificial water table: DTBLDIG',WTBLDepz_tile
    write(*,*)'slope of natural water table relative to landscape '// &
      'surface: DTBLGG',DTBLGG
    write(*,*)'boundary condns for N surface runoff: RCHQNG',RCHQNG
    write(*,*)'boundary condns for E surface runoff: RCHQEG',RCHQEG
    write(*,*)'boundary condns for S surface runoff: RCHQSG',RCHQSG
    write(*,*)'boundary condns for W surface runoff: RCHQWG',RCHQWG
    write(*,*)'bound condns for N subsurf flow: RCHGNUG',RCHGNUG
    write(*,*)'bound condns for E subsurf flow: RCHGEUG',RCHGEUG
    write(*,*)'bound condns for S subsurf flow: RCHGSUG',RCHGSUG
    write(*,*)'bound condns for W subsurf flow: RCHGWUG',RCHGWUG
    write(*,*)'N distance to water table (m): RCHGNTG',RCHGNTG
    write(*,*)'E distance to water table (m): RCHGETG',RCHGETG
    write(*,*)'S distance to water table (m): RCHGSTG',RCHGSTG
    write(*,*)'W distance to water table (m): RCHGWTG',RCHGWTG
    write(*,*)'lower boundary conditions for water flow:RCHGDG', RCHGDG
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'width of each W-E landscape column: DHI'
    write(*,*)(DHI(NX),NX=1,NHE)
    write(*,*)'width of each N-S landscape row: DVI'
    write(*,*)(DVI(NY),NY=1,NVS)
    write(*,'(100A)')('=',ll=1,100)
  endif

  D9895: DO NX=NHW,NHE
    D9890: DO NY=NVN,NVS
      ALAT(NY,NX)=ALATG
      PBOT_col(NY,NX)=PBOT_col(NY,NX)*exp(-ALT(NY,NX)/hpresc)
      ALTI(NY,NX)=ALTIG
      ATCAI(NY,NX)=ATCAG
      IDWaterTable(NY,NX)=iWaterTabelMode
      OXYE(NY,NX)=ao2_ppm
      Z2GE(NY,NX)=an2_ppm
      CO2EI(NY,NX)=aco2_ppm
      CH4E_col(NY,NX)=ach4_ppm
      Z2OE(NY,NX)=an2o_ppm
      ZNH3E_col(NY,NX)=ZNH3EG
      KoppenClimZone_col(NY,NX)   = IETYPG
      FlowDirIndicator(NY,NX)     = grid_mode
      NatWtblDepz_col(NY,NX)      = WTBLDepz_nat
      WtblDepzTile_col(NY,NX)     = WTBLDepz_tile
      WaterTBLSlope(NY,NX)        = DTBLGG

      RechargNorthSurf(NY,NX)     = RCHQNG
      RechargEastSurf(NY,NX)      = RCHQEG
      RechargSouthSurf(NY,NX)     = RCHQSG
      RechargWestSurf(NY,NX)      = RCHQWG

      RechargNorthSubSurf(NY,NX)  = RCHGNUG
      RechargEastSubSurf(NY,NX)   = RCHGEUG
      RechargSouthSubSurf(NY,NX)  = RCHGSUG
      RechargWestSubSurf(NY,NX)   = RCHGWUG

      RechargRateNorthWTBL(NY,NX) = RCHGNTG
      RechargRateEastWTBL(NY,NX)  = RCHGETG
      RechargRateSouthWTBL(NY,NX) = RCHGSTG
      RechargRateWestWTBL(NY,NX)  = RCHGWTG
      RCHGD(NY,NX)                = RCHGDG
      DH(NY,NX)=DHI(NX)
      DV(NY,NX)=DVI(NY)
      CO2E_col(NY,NX)=CO2EI(NY,NX)
      H2GE(NY,NX)=1.0E-03_r8
!
!     CALCULATE MAXIMUM DAYLENTH FOR PLANT PHENOLOGY
!
!     DayLenthMax=maximum daylength (h)
!
      IF(ALAT(NY,NX).GT.0.0_r8)THEN
        XI=173
      ELSE
        XI=356
      ENDIF
      DayLenthMax(NY,NX)=GetDayLength(ALAT(NY,NX),XI)

    ENDDO D9890
  ENDDO D9895

  end subroutine readsiteNC

!------------------------------------------------------------------------------------------

  subroutine readTopoNC()
!
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer :: jj,NY,NX,ll
  character(len=200) :: tline
  integer :: NM(JY,JX),ntp,ntopus
  real(r8) :: dat1(1:JZ)
  real(r8) :: corrector
! begin_execution
  associate(                            &
  k_woody_litr => micpar%k_woody_litr , &
  k_fine_litr  => micpar%k_fine_litr  , &
  k_manure     => micpar%k_manure       &
  )

!
! READ TOPOGRAPHY DATA AND SOIL FILE NAME FOR EACH GRID CELL
!
! for each unit within the landscape:
! NH1,NV1,NH2,NV2=NW,SE column,row
! ASPX=N,E,S,W aspect (o)
! SL0= slope (o)
! DPTHSX=initial snowpack depth
! SLX: not used

  ntopus=get_dim_len(grid_nfid, 'ntopou')
  if(first_topou)ntopus=1
  DO ntp=1,ntopus
    call ncd_getvar(grid_nfid, 'NH1', ntp, NH1)
    call ncd_getvar(grid_nfid, 'NV1', ntp, NV1)
    call ncd_getvar(grid_nfid, 'NH2', ntp, NH2)
    call ncd_getvar(grid_nfid, 'NV2', ntp, NV2)
    call ncd_getvar(grid_nfid, 'ASPX', ntp,ASPX)
    call ncd_getvar(grid_nfid, 'SL0', ntp,SL0)
    call ncd_getvar(grid_nfid, 'DPTHSX', ntp,initSnowDepth)

!
! OPEN AND READ SOIL FILE
!

    DO  NX=NH1,NH2
      DO  NY=NV1,NV2
!
!     SURFACE SLOPES AND ASPECTS
!
        ASP_col(NY,NX)=ASPX
        SL(NY,NX)=SL0
        SnowDepth_col(NY,NX)=initSnowDepth
!
!     CONVERT ASPECT from geographic format TO GEOMETRIC FORMAT
!
!     what is geometric format mean? geographic format 0 is north, 90 east, 180 south,
!     geometric format 0/360 is east, counting as counterclock wise
        ASP_col(NY,NX)=450.0_r8-ASP_col(NY,NX)
        IF(ASP_col(NY,NX).GE.360.0_r8)ASP_col(NY,NX)=ASP_col(NY,NX)-360.0_r8
      ENDDO
    ENDDO

    call ncd_getvar(grid_nfid, 'PSIFC', ntp,PSIAtFldCapacity(NV1,NH1))
    call ncd_getvar(grid_nfid, 'PSIWP', ntp,PSIAtWiltPoint(NV1,NH1))
    call ncd_getvar(grid_nfid, 'ALBS',  ntp,SoilAlbedo_col(NV1,NH1))
    call ncd_getvar(grid_nfid, 'PH0',   ntp,PH(0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSCf',  ntp,RSC(k_fine_litr,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSNf',  ntp,RSN(k_fine_litr,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSPf',  ntp,RSP(k_fine_litr,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSCw',  ntp,RSC(k_woody_litr,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSNw',  ntp,RSN(k_woody_litr,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSPw',  ntp,RSP(k_woody_litr,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSCm',  ntp,RSC(k_manure,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSNm',  ntp,RSN(k_manure,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSPm',  ntp,RSP(k_manure,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'IXTYP1',ntp,IXTYP(1,NV1,NH1))
    call ncd_getvar(grid_nfid, 'IXTYP2',ntp,IXTYP(2,NV1,NH1))
    call ncd_getvar(grid_nfid, 'NUI'   ,ntp,NUI(NV1,NH1))
    call ncd_getvar(grid_nfid, 'NJ'    ,ntp,MaxNumRootLays(NV1,NH1))
    call ncd_getvar(grid_nfid, 'NL1'   ,ntp,NL1)
    call ncd_getvar(grid_nfid, 'NL2'   ,ntp,NL2)
    call ncd_getvar(grid_nfid, 'ISOILR',ntp,ISOILR(NV1,NH1))

    NU(NV1,NH1)=NUI(NV1,NH1)
    NK(NV1,NH1)=MaxNumRootLays(NV1,NH1)+1
    NM(NV1,NH1)=MaxNumRootLays(NV1,NH1)+NL1
!  the extra soil layer below root zone cannot be greater than what is allowed
    NL2=min0(JZ-NM(NV1,NH1),NL2)
    NLI(NV1,NH1)=NM(NV1,NH1)+NL2
    NL(NV1,NH1)=NLI(NV1,NH1)

    call ncd_getvar(grid_nfid, 'CDPTH',ntp,CumDepz2LayerBot_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'BKDSI',ntp,SoiBulkDensityt0_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'FC', ntp,FieldCapacity_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'WP', ntp,WiltPoint_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'SCNV', ntp,SatHydroCondVert_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'SCNH', ntp,SatHydroCondHrzn_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CSAND',ntp,CSAND(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CSILT',ntp,CSILT(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'FHOL',ntp,SoilFracAsMacP_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'ROCK',ntp,ROCK_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'PH',ntp,PH(1:JZ,NV1,NH1))
    !meq/100g means 1.e-3 mol/100 g =1 cmol/kg , cmol ~ 1.e-2 mol
    call ncd_getvar(grid_nfid, 'CEC',ntp,CEC_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'AEC',ntp,AEC_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CORGC',ntp,CSoilOrgM_vr(ielmc,1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CORGR',ntp,COMLitrC_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CORGN',ntp,CSoilOrgM_vr(ielmn,1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CORGP',ntp,CSoilOrgM_vr(ielmp,1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CNH4',ntp,CNH4(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CNO3',ntp,CNO3(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CPO4',ntp,CPO4(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CAL' ,ntp,CAL(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CFE' ,ntp,CFE(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCA' ,ntp,CCA(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CMG' ,ntp,CMG(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CNA' ,ntp,CNA(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CKA' ,ntp,CKA(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CSO4',ntp,CSO4(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCL',ntp,CCL(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CALPO',ntp,CALPO(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CFEPO',ntp,CFEPO(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCAPD',ntp,CCAPD(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCAPH',ntp,CCAPH(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CALOH',ntp,CALOH(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CFEOH',ntp,CFEOH(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCACO',ntp,CCACO(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCASO',ntp,CCASO(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'GKC4',ntp,GKC4(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'GKCH',ntp,GKCH(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'GKCA',ntp,GKCA(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'GKCM',ntp,GKCM(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'GKCN',ntp,GKCN(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'GKCK',ntp,GKCK(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'THW',ntp,THW(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'THI',ntp,THI(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'RSCfL',ntp,dat1(1:JZ));RSC(k_fine_litr,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSNfL',ntp,dat1(1:JZ));RSN(k_fine_litr,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSPfL',ntp,dat1(1:JZ));RSP(k_fine_litr,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSCwL',ntp,dat1(1:JZ));RSC(k_woody_litr,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSNwL',ntp,dat1(1:JZ));RSN(k_woody_litr,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSPwL',ntp,dat1(1:JZ));RSP(k_woody_litr,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSCmL',ntp,dat1(1:JZ));RSC(k_manure,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSNmL',ntp,dat1(1:JZ));RSN(k_manure,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSPmL',ntp,dat1(1:JZ));RSP(k_manure,1:JZ,NV1,NH1)=dat1(1:JZ)

    DO  NX=NH1,NH2
      DO  NY=NV1,NV2
        IF (NX/=NH1 .OR. NY/=NV1) THEN
!
!     SURFACE PROPERTIES
!
!     PSIAtFldCapacity,PSIAtWiltPoint=water potentials at field capacity,wilting point (MPa)
!     SoilAlbedo=wet soil albedo
!     PH=litter pH
!     RSC,RSC,RSP=C,N,P in fine(1,0),woody(0,0),manure(2,0) surface litter (g m-2)
!     IXTYP=surface litter type:1=plant,2=manure
!     NUI,MaxNumRootLays=number of soil surface layer,maximum rooting layer
!     NL1,NL2=number of additional layers below NJ with,without data in file
!     ISOILR=natural(0),reconstructed(1) soil profile
!
          PSIAtFldCapacity(NY,NX)=PSIAtFldCapacity(NV1,NH1)
          PSIAtWiltPoint(NY,NX)=PSIAtWiltPoint(NV1,NH1)
          SoilAlbedo_col(NY,NX) =SoilAlbedo_col(NV1,NH1)
          PH(0,NY,NX) =PH(0,NV1,NH1)
          RSC(k_fine_litr,0,NY,NX) =RSC(k_fine_litr,0,NV1,NH1)
          RSN(k_fine_litr,0,NY,NX) =RSN(k_fine_litr,0,NV1,NH1)
          RSP(k_fine_litr,0,NY,NX) =RSP(k_fine_litr,0,NV1,NH1)
          RSC(k_woody_litr,0,NY,NX) =RSC(k_woody_litr,0,NV1,NH1)
          RSN(k_woody_litr,0,NY,NX) =RSN(k_woody_litr,0,NV1,NH1)
          RSP(k_woody_litr,0,NY,NX) =RSP(k_woody_litr,0,NV1,NH1)
          RSC(k_manure,0,NY,NX) =RSC(k_manure,0,NV1,NH1)
          RSN(k_manure,0,NY,NX) =RSN(k_manure,0,NV1,NH1)
          RSP(k_manure,0,NY,NX) =RSP(k_manure,0,NV1,NH1)
          IXTYP(1,NY,NX) =IXTYP(1,NV1,NH1)
          IXTYP(2,NY,NX) =IXTYP(2,NV1,NH1)
          NUI(NY,NX) = NUI(NV1,NH1)
          MaxNumRootLays(NY,NX)  = MaxNumRootLays(NV1,NH1)
          ISOILR(NY,NX)=ISOILR(NV1,NH1)
          NU(NY,NX)=NU(NV1,NH1)
          NK(NY,NX)=NK(NV1,NH1)
          NM(NY,NX)=NM(NV1,NH1)
!  the extra soil layer below root zone cannot be greater than what is allowed
          NLI(NY,NX)=NLI(NV1,NH1)
          NL(NY,NX)=NLI(NV1,NH1)
        ENDIF

!
!     PHYSICAL PROPERTIES
!
!     CDPTH=depth to bottom (m) > 0
!     SoiBulkDensityt0_vr=initial bulk density (Mg m-3,0=water), it refers to solid matter
!
!
        IF (NX/=NH1 .OR. NY/=NV1) THEN
          DO L=NU(NY,NX),NM(NY,NX)
            CumDepz2LayerBot_vr(L,NY,NX) = CumDepz2LayerBot_vr(L,NV1,NH1)
            SoiBulkDensityt0_vr(L,NY,NX) = SoiBulkDensityt0_vr(L,NV1,NH1)
            FieldCapacity_vr(L,NY,NX)    = FieldCapacity_vr(L,NV1,NH1)
            WiltPoint_vr(L,NY,NX)        = WiltPoint_vr(L,NV1,NH1)
            SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L,NV1,NH1)
            SatHydroCondHrzn_vr(L,NY,NX) = SatHydroCondHrzn_vr(L,NV1,NH1)
            CSAND(L,NY,NX)               = CSAND(L,NV1,NH1)
            CSILT(L,NY,NX)               = CSILT(L,NV1,NH1)
            SoilFracAsMacP_vr(L,NY,NX)   = SoilFracAsMacP_vr(L,NV1,NH1)
            ROCK_vr(L,NY,NX)             = ROCK_vr(L,NV1,NH1)
            PH(L,NY,NX)                  = PH(L,NV1,NH1)
            CEC_vr(L,NY,NX)                 = CEC_vr(L,NV1,NH1)
            AEC_vr(L,NY,NX)                 = AEC_vr(L,NV1,NH1)
            CSoilOrgM_vr(ielmc,L,NY,NX)  = CSoilOrgM_vr(ielmc,L,NV1,NH1)
            COMLitrC_vr(L,NY,NX)         = COMLitrC_vr(L,NV1,NH1)
            CSoilOrgM_vr(ielmn,L,NY,NX)  = CSoilOrgM_vr(ielmn,L,NV1,NH1)
            CSoilOrgM_vr(ielmp,L,NY,NX)  = CSoilOrgM_vr(ielmp,L,NV1,NH1)

            CNH4(L,NY,NX) = CNH4(L,NV1,NH1)
            CNO3(L,NY,NX) = CNO3(L,NV1,NH1)
            CPO4(L,NY,NX) = CPO4(L,NV1,NH1)

            CAL(L,NY,NX)  = CAL(L,NV1,NH1)
            CFE(L,NY,NX)  = CFE(L,NV1,NH1)
            CCA(L,NY,NX)  = CCA(L,NV1,NH1)
            CMG(L,NY,NX)  = CMG(L,NV1,NH1)
            CNA(L,NY,NX)  = CNA(L,NV1,NH1)
            CKA(L,NY,NX)  = CKA(L,NV1,NH1)
            CSO4(L,NY,NX) = CSO4(L,NV1,NH1)
            CCL(L,NY,NX)  = CCL(L,NV1,NH1)

            CALPO(L,NY,NX) = CALPO(L,NV1,NH1)
            CFEPO(L,NY,NX) = CFEPO(L,NV1,NH1)
            CCAPD(L,NY,NX) = CCAPD(L,NV1,NH1)
            CCAPH(L,NY,NX) = CCAPH(L,NV1,NH1)
            CALOH(L,NY,NX) = CALOH(L,NV1,NH1)
            CFEOH(L,NY,NX) = CFEOH(L,NV1,NH1)
            CCACO(L,NY,NX) = CCACO(L,NV1,NH1)
            CCASO(L,NY,NX) = CCASO(L,NV1,NH1)

            GKC4(L,NY,NX) = GKC4(L,NV1,NH1)
            GKCH(L,NY,NX) = GKCH(L,NV1,NH1)
            GKCA(L,NY,NX) = GKCA(L,NV1,NH1)
            GKCM(L,NY,NX) = GKCM(L,NV1,NH1)
            GKCN(L,NY,NX) = GKCN(L,NV1,NH1)
            GKCK(L,NY,NX) = GKCK(L,NV1,NH1)

            THW(L,NY,NX) = THW(L,NV1,NH1)
            THI(L,NY,NX) = THI(L,NV1,NH1)

            RSC(k_fine_litr,L,NY,NX)  = RSC(k_fine_litr,L,NV1,NH1)
            RSN(k_fine_litr,L,NY,NX)  = RSN(k_fine_litr,L,NV1,NH1)
            RSP(k_fine_litr,L,NY,NX)  = RSP(k_fine_litr,L,NV1,NH1)
            RSC(k_woody_litr,L,NY,NX) = RSC(k_woody_litr,L,NV1,NH1)
            RSN(k_woody_litr,L,NY,NX) = RSN(k_woody_litr,L,NV1,NH1)
            RSP(k_woody_litr,L,NY,NX) = RSP(k_woody_litr,L,NV1,NH1)
            RSC(k_manure,L,NY,NX)     = RSC(k_manure,L,NV1,NH1)
            RSN(k_manure,L,NY,NX)     = RSN(k_manure,L,NV1,NH1)
            RSP(k_manure,L,NY,NX)     = RSP(k_manure,L,NV1,NH1)

          ENDDO
        ENDIF

        if(lverb)then
          CALL Disp_topo_charc(NY,NX,NU(NY,NX),NM(NY,NX))
        endif
        RSC(k_fine_litr,0,NY,NX)     = AMAX1(ppmc,RSC(k_fine_litr,0,NY,NX))
        RSN(k_fine_litr,0,NY,NX)     = AMAX1(0.04E-06_r8,RSN(k_fine_litr,0,NY,NX))
        RSP(k_fine_litr,0,NY,NX)     = AMAX1(0.004E-06_r8,RSP(k_fine_litr,0,NY,NX))
        SatHydroCondVert_vr(0,NY,NX) = 10.0_r8*0.098_r8
!
!     SET FLAGS FOR ESTIMATING FC,WP,SCNV,SCNH IF UNKNOWN
!
!     ISOIL=flag for calculating FC(1),WiltPoint_vr(2),SatHydroCondVert_vr(3),SatHydroCondHrzn_vr(4)
!
        call ComputeSoilHydroPars(NY,NX,NU(NY,NX),NM(NY,NX))

!
!     FILL OUT SOIL BOUNDARY LAYERS ABOVE ROOTING ZONE (NOT USED)
!     below is for soil repacking, whence NU>1
!     root zone from NU to NM, why use 0.025?
!
        IF(NU(NY,NX).GT.1)THEN
          DO  L=NU(NY,NX)-1,0,-1
            IF(SoiBulkDensityt0_vr(L+1,NY,NX).GT.0.025_r8)THEN
              CumDepz2LayerBot_vr(L,NY,NX)=CumDepz2LayerBot_vr(L+1,NY,NX)-0.01_r8
            ELSE
              CumDepz2LayerBot_vr(L,NY,NX)=CumDepz2LayerBot_vr(L+1,NY,NX)-0.02_r8
            ENDIF
            IF(L.GT.0)THEN
              SoiBulkDensityt0_vr(L,NY,NX) = SoiBulkDensityt0_vr(L+1,NY,NX)
              FieldCapacity_vr(L,NY,NX)    = FieldCapacity_vr(L+1,NY,NX)
              WiltPoint_vr(L,NY,NX)        = WiltPoint_vr(L+1,NY,NX)
              SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L+1,NY,NX)
              SatHydroCondHrzn_vr(L,NY,NX) = SatHydroCondHrzn_vr(L+1,NY,NX)
              CSAND(L,NY,NX)               = CSAND(L+1,NY,NX)
              CSILT(L,NY,NX)               = CSILT(L+1,NY,NX)
              CCLAY(L,NY,NX)               = CCLAY(L+1,NY,NX)
              SoilFracAsMacP_vr(L,NY,NX)   = SoilFracAsMacP_vr(L+1,NY,NX)
              ROCK_vr(L,NY,NX)             = ROCK_vr(L+1,NY,NX)
              PH(L,NY,NX)                  = PH(L+1,NY,NX)
              CEC_vr(L,NY,NX)                 = CEC_vr(L+1,NY,NX)
              AEC_vr(L,NY,NX)                 = AEC_vr(L+1,NY,NX)
              CSoilOrgM_vr(ielmc,L,NY,NX)  = 1.0_r8*CSoilOrgM_vr(ielmc,L+1,NY,NX)
              COMLitrC_vr(L,NY,NX)         = 1.0_r8*COMLitrC_vr(L+1,NY,NX)
              CSoilOrgM_vr(ielmn,L,NY,NX)  = 1.0_r8*CSoilOrgM_vr(ielmn,L+1,NY,NX)
              CSoilOrgM_vr(ielmp,L,NY,NX)  = 1.0_r8*CSoilOrgM_vr(ielmp,L+1,NY,NX)
              CNH4(L,NY,NX)                = CNH4(L+1,NY,NX)
              CNO3(L,NY,NX)                = CNO3(L+1,NY,NX)
              CPO4(L,NY,NX)                = CPO4(L+1,NY,NX)
              CAL(L,NY,NX)                 = CAL(L+1,NY,NX)
              CFE(L,NY,NX)                 = CFE(L+1,NY,NX)
              CCA(L,NY,NX)                 = CCA(L+1,NY,NX)
              CMG(L,NY,NX)                 = CMG(L+1,NY,NX)
              CNA(L,NY,NX)                 = CNA(L+1,NY,NX)
              CKA(L,NY,NX)                 = CKA(L+1,NY,NX)
              CSO4(L,NY,NX)                = CSO4(L+1,NY,NX)
              CCL(L,NY,NX)                 = CCL(L+1,NY,NX)
              CALOH(L,NY,NX)               = CALOH(L+1,NY,NX)
              CFEOH(L,NY,NX)               = CFEOH(L+1,NY,NX)
              CCACO(L,NY,NX)               = CCACO(L+1,NY,NX)
              CCASO(L,NY,NX)               = CCASO(L+1,NY,NX)
              CALPO(L,NY,NX)               = CALPO(L+1,NY,NX)
              CFEPO(L,NY,NX)               = CFEPO(L+1,NY,NX)
              CCAPD(L,NY,NX)               = CCAPD(L+1,NY,NX)
              CCAPH(L,NY,NX)               = CCAPH(L+1,NY,NX)
              GKC4(L,NY,NX)                = GKC4(L+1,NY,NX)
              GKCH(L,NY,NX)                = GKCH(L+1,NY,NX)
              GKCA(L,NY,NX)                = GKCA(L+1,NY,NX)
              GKCM(L,NY,NX)                = GKCM(L+1,NY,NX)
              GKCN(L,NY,NX)                = GKCN(L+1,NY,NX)
              GKCK(L,NY,NX)                = GKCK(L+1,NY,NX)
              THW(L,NY,NX)                 = THW(L+1,NY,NX)
              THI(L,NY,NX)                 = THI(L+1,NY,NX)
              ISOIL(1:4,L,NY,NX)           = ISOIL(1:4,L+1,NY,NX)
              RSC(k_fine_litr,L,NY,NX)     = 0.0_r8
              RSN(k_fine_litr,L,NY,NX)     = 0.0_r8
              RSP(k_fine_litr,L,NY,NX)     = 0.0_r8
              RSC(k_woody_litr,L,NY,NX)    = 0.0_r8
              RSN(k_woody_litr,L,NY,NX)    = 0.0_r8
              RSP(k_woody_litr,L,NY,NX)    = 0.0_r8
              RSC(k_manure,L,NY,NX)        = 0.0_r8
              RSN(k_manure,L,NY,NX)        = 0.0_r8
              RSP(k_manure,L,NY,NX)        = 0.0_r8
            ENDIF
          ENDDO
        ENDIF
!
!     ADD SOIL BOUNDARY LAYERS BELOW SOIL ZONE
!     depth of layer (L-1) is at the middle between that of layer L-2 and L
        if(lverb)write(*,*)'SetDeepSoil'
        call SetDeepSoil(NY,NX,NM(NY,NX),JZ)

!
!   CALCULATE DERIVED SOIL PROPERTIES FROM INPUT SOIL PROPERTIES
!
!   FracSoiAsMicP_vr=micropore fraction excluding macropore,rock
!   SCNV,SCNH=vertical,lateral Ksat converted to m2 MPa-1 h-1
!   CSAND,CSILT,CCLAY=sand,silt,clay content converted to g Mg-1
!   CORGC,CORGR=SOC,POC converted to g Mg-1
!   CEC,AEC=cation,anion exchange capacity converted to mol Mg-1
!   CNH4...=solute concentrations converted to mol Mg-1
!   SoiBulkDensityt0_vr: initial bulk density

        DO  L=1,NL(NY,NX)
  !   SoilFracAsMacP: macropore fraction
  !     SoiBulkDensityt0_vr(L,NY,NX)=SoiBulkDensityt0_vr(L,NY,NX)/(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
          SoiBulkDensity_vr(L,NY,NX)=SoiBulkDensityt0_vr(L,NY,NX)
          IF(isclose(SoiBulkDensity_vr(L,NY,NX),0.0_r8))SoilFracAsMacP_vr(L,NY,NX)=0.0_r8
  !     fraction of soil as micropore
          FracSoiAsMicP_vr(L,NY,NX)=(1.0_r8-ROCK_vr(L,NY,NX))*(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
  !  Macropore correction is off, when reporting from measurements, FieldCapacity includes contribution from
  !  both macropores and micropores    
  !     FieldCapacity_vr(L,NY,NX)=FieldCapacity_vr(L,NY,NX)/(1.0-SoilFracAsMacP_vr(L,NY,NX))
  !     WiltPoint_vr(L,NY,NX)=WiltPoint_vr(L,NY,NX)/(1.0-SoilFracAsMacP_vr(L,NY,NX))
  !
          SatHydroCondVert_vr(L,NY,NX)=0.098_r8*SatHydroCondVert_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
          SatHydroCondHrzn_vr(L,NY,NX)=0.098_r8*SatHydroCondHrzn_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
          CCLAY(L,NY,NX)=AZMAX1(1.0E+03_r8-(CSAND(L,NY,NX)+CSILT(L,NY,NX)))
          CSoilOrgM_vr(ielmc,L,NY,NX)=CSoilOrgM_vr(ielmc,L,NY,NX)*1.0E+03_r8   !convert from Kg to g C
          COMLitrC_vr(L,NY,NX)=COMLitrC_vr(L,NY,NX)*1.0E+03_r8   !convert from Kg to g C
          CORGCI(L,NY,NX)=CSoilOrgM_vr(ielmc,L,NY,NX)
          SoilFracAsMacPt0_vr(L,NY,NX)=SoilFracAsMacP_vr(L,NY,NX)
  ! soil texture is reported based on mass basis soley for mineral component of the soil
          corrector=1.0E-03_r8*AZMAX1((1.0_r8-CSoilOrgM_vr(ielmc,L,NY,NX)/orgcden))
          CSAND(L,NY,NX)=CSAND(L,NY,NX)*corrector
          CSILT(L,NY,NX)=CSILT(L,NY,NX)*corrector
          CCLAY(L,NY,NX)=CCLAY(L,NY,NX)*corrector
          CEC_vr(L,NY,NX)=CEC_vr(L,NY,NX)*10.0_r8   !convert from meq/100g to cmol/kg
          AEC_vr(L,NY,NX)=AEC_vr(L,NY,NX)*10.0_r8   !convert from meq/100g to cmol/kg
          CNH4(L,NY,NX)=CNH4(L,NY,NX)/natomw
          CNO3(L,NY,NX)=CNO3(L,NY,NX)/natomw
          CPO4(L,NY,NX)=CPO4(L,NY,NX)/patomw
          CAL(L,NY,NX)=CAL(L,NY,NX)/27.0_r8
          CFE(L,NY,NX)=CFE(L,NY,NX)/56.0_r8
          CCA(L,NY,NX)=CCA(L,NY,NX)/40.0_r8
          CMG(L,NY,NX)=CMG(L,NY,NX)/24.3_r8
          CNA(L,NY,NX)=CNA(L,NY,NX)/23.0_r8
          CKA(L,NY,NX)=CKA(L,NY,NX)/39.1_r8
          CSO4(L,NY,NX)=CSO4(L,NY,NX)/32.0_r8
          CCL(L,NY,NX)=CCL(L,NY,NX)/35.5_r8
          CALPO(L,NY,NX)=CALPO(L,NY,NX)/patomw
          CFEPO(L,NY,NX)=CFEPO(L,NY,NX)/patomw
          CCAPD(L,NY,NX)=CCAPD(L,NY,NX)/patomw
          CCAPH(L,NY,NX)=CCAPH(L,NY,NX)/(patomw*3.0_r8)
          CALOH(L,NY,NX)=CALOH(L,NY,NX)/27.0_r8
          CFEOH(L,NY,NX)=CFEOH(L,NY,NX)/56.0_r8
          CCACO(L,NY,NX)=CCACO(L,NY,NX)/40.0_r8
          CCASO(L,NY,NX)=CCASO(L,NY,NX)/40.0_r8
  !
  !     ESTIMATE SON,SOP,CEC IF UNKNOWN
  !     BIOCHEMISTRY 130:117-131
  !
          IF(CSoilOrgM_vr(ielmn,L,NY,NX).LT.0.0_r8)THEN
  !  default ORGN parameterization
            CSoilOrgM_vr(ielmn,L,NY,NX)=AMIN1(0.125_r8*CSoilOrgM_vr(ielmc,L,NY,NX),&
              8.9E+02_r8*(CSoilOrgM_vr(ielmc,L,NY,NX)/1.0E+04_r8)**0.80_r8)
          ENDIF
          IF(CSoilOrgM_vr(ielmp,L,NY,NX).LT.0.0_r8)THEN
            CSoilOrgM_vr(ielmp,L,NY,NX)=AMIN1(0.0125_r8*CSoilOrgM_vr(ielmc,L,NY,NX),&
              1.2E+02_r8*(CSoilOrgM_vr(ielmc,L,NY,NX)/1.0E+04_r8)**0.52_r8)
          ENDIF
          IF(CEC_vr(L,NY,NX).LT.0.0_r8)THEN
            !estimate from input data
            CEC_vr(L,NY,NX)=10.0_r8*(200.0_r8*2.0_r8*CSoilOrgM_vr(ielmc,L,NY,NX)/1.0E+06_r8 &
              +80.0_r8*CCLAY(L,NY,NX)+20.0_r8*CSILT(L,NY,NX) &
              +5.0_r8*CSAND(L,NY,NX))
          ENDIF
        ENDDO
        CSoilOrgM_vr(ielmc,0,NY,NX)=orgcden
        FracSoiAsMicP_vr(0,NY,NX)=1.0_r8
      ENDDO
    ENDDO
  ENDDO
  end associate
  end subroutine readTopoNC

!------------------------------------------------------------------------------------------

  subroutine Disp_topo_charc(NY,NX,NU,NM)

  implicit none
  integer, intent(in) :: NY,NX,NU,NM
  integer :: LL,L

  associate(                            &
  k_woody_litr => micpar%k_woody_litr , &
  k_fine_litr  => micpar%k_fine_litr  , &
  k_manure     => micpar%k_manure       &
  )

  write(*,'(40A)')('-',ll=1,40)
  write(*,*)'Topographic characterization'
  write(*,'(40A)')('-',ll=1,40)
  write(*,*)'NY, NX =',NY,NX
  write(*,*)'Aspect (o): ASPX',ASP_col(NY,NX)
  write(*,*)'Slope (o): SL0',SL(NY,NX)
  write(*,*)'Initial snowpack depth: initSnowDepth',SnowDepth_col(NY,NX)
  write(*,'(100A)')('=',ll=1,100)

  write(*,*)''
  write(*,*)'NY,NX=',NY,NX
  write(*,*)'Water potential at field capacity (MPa)',PSIAtFldCapacity(NY,NX)
  write(*,*)'Water potential at wilting point (MPa)',PSIAtWiltPoint(NY,NX)
  write(*,*)'Wet soil albedo',SoilAlbedo_col(NY,NX)

  write(*,*)'Litter pH',PH(0,NY,NX)
  write(*,*)'C in surface fine litter (g m-2)',RSC(k_fine_litr,0,NY,NX)
  write(*,*)'N in surface fine litter (g m-2)',RSN(k_fine_litr,0,NY,NX)
  write(*,*)'P in surface fine litter (g m-2)',RSP(k_fine_litr,0,NY,NX)
  write(*,*)'C in surface woody litter (g m-2)',RSC(k_woody_litr,0,NY,NX)
  write(*,*)'N in surface woody litter (g m-2)',RSN(k_woody_litr,0,NY,NX)
  write(*,*)'P in surface woody litter (g m-2)',RSP(k_woody_litr,0,NY,NX)
  write(*,*)'C in surface manure litter (g m-2)',RSC(k_manure,0,NY,NX)
  write(*,*)'N in surface manure litter (g m-2)',RSN(k_manure,0,NY,NX)
  write(*,*)'P in surface manure litter (g m-2)',RSP(k_manure,0,NY,NX)
  write(*,*)'surface litter type:1=plant,2=manure',IXTYP(1,NY,NX),IXTYP(2,NY,NX)
  write(*,*)'layer number of soil surface layer NUI',NUI(NY,NX)
  write(*,*)'layer number of maximum rooting layer NJ',MaxNumRootLays(NY,NX)
  write(*,*)'number of layers involved in model calculation',NL(NY,NX)
  write(*,*)'Flag for natural(0),reconstructed(1) soil profile', ISOILR(NY,NX)
  write(*,*)
  write(*,'(A,I2,A,I2)')'read data for layers from layer NU ',NU,' to layer NM ',NM

  write(*,*)'Depth to bottom of soil layer (m): CDPTH'
  write(*,*)(CumDepz2LayerBot_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Initial bulk density (Mg m-3, 0=water): SoiBulkDensityt0_vr'
  write(*,*)(SoiBulkDensityt0_vr(L,NY,NX),L=NU,NM)
!
!     HYDROLOGIC PROPERTIES
!
!     FC,WP=field capacity,wilting point:<0=unknown (m3 m-3)
!     SCNV,SCNH=vertical,lateral Ksat:<0=unknown (mm h-1)
!

  write(*,*)''
  write(*,*)'NY,NX=',NY,NX
  write(*,*)'L=',NU,NM
  write(*,*)'Field capacity (m3 m-3): FC'
  write(*,*)(FieldCapacity_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Wilting point (m3 m-3): WP'
  write(*,*)(WiltPoint_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Vertical Ksat (mm h-1): SCNV'
  write(*,*)(SatHydroCondVert_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Lateral Ksat (mm h-1): SCNH'
  write(*,*)(SatHydroCondHrzn_vr(L,NY,NX),L=NU,NM)
!
!     PHYSICAL PROPERTIES
!
!     CSAND,CSILT=sand,silt contents (kg Mg-1)
!     SoilFracAsMacP_vr,ROCK=macropore,rock fraction
!

  write(*,*)''
  write(*,*)'Sand (kg Mg-1): CSAND'
  write(*,*)(CSAND(L,NY,NX),L=NU,NM)
  write(*,*)'Silt (kg Mg-1): CSILT'
  write(*,*)(CSILT(L,NY,NX),L=NU,NM)
  write(*,*)'Macropore fraction (0-1): FOHL'
  write(*,*)(SoilFracAsMacP_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Rock fraction (0-1): ROCK'
  write(*,*)(ROCK_vr(L,NY,NX),L=NU,NM)
!
!     CHEMICAL PROPERTIES
!
!     PH=pH
!     CEC,AEC=cation,anion exchange capacity:CEC<0=unknown (cmol Kg-1)
!

  write(*,*)''
  write(*,*)'pH'
  write(*,*)(PH(L,NY,NX),L=NU,NM)
  write(*,*)'Cation exchange capacity (cmol Kg-1): CEC'
  write(*,*)(CEC_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Anion exchange capacity (cmol Kg-1): AEC'
  write(*,*)(AEC_vr(L,NY,NX),L=NU,NM)
!
!     ORGANIC C, N AND P CONCENTRATIONS
!
!     CORGC,COMLitrC_vr=total SOC,POC(part of SOC) (kg Mg-1)
!     CORGN,CORGP=SON,SOP:<0=unknown (g Mg-1)
!
  write(*,*)''
  write(*,*)'Total SOC (kg C/Mg soil): CORGC'
  write(*,*)(CSoilOrgM_vr(ielmc,L,NY,NX),L=NU,NM)
  write(*,*)'POC (part of SOC) (kg C/Mg soil): CORGR '
  write(*,*)(COMLitrC_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Total SON (g N/Mg soil): CORGN '
  write(*,*)(CSoilOrgM_vr(ielmn,L,NY,NX),L=NU,NM)
  write(*,*)'Total SOP (g P/Mg soil): CORGP'
  write(*,*)(CSoilOrgM_vr(ielmp,L,NY,NX),L=NU,NM)
!
!     INORGANIC N AND P CONCENTRATIONS
!
!     CNH4,CNO3,CPO4=soluble+exchangeable NH4,NO3,H2PO4 (g Mg-1)
!
  write(*,*)''
  write(*,*)'Total soil NH4 concentration (gN Mg-1): CNH4 '
  write(*,*)(CNH4(L,NY,NX),L=NU,NM)
  write(*,*)'Total soil NO3 concentration (gN Mg-1): CNO3'
  write(*,*)(CNO3(L,NY,NX),L=NU,NM)
  write(*,*)'Total soil H2PO4 concentration (gP Mg-1): CPO4 '
  write(*,*)(CPO4(L,NY,NX),L=NU,NM)
!
!     CATION AND ANION CONCENTRATIONS
!
!     C*=soluble concentration from sat. paste extract (g Mg-1)
!     AL,FE,CA,MG,NA,KA,SO4,CL=Al,Fe,Ca,Mg,Na,K,SO4-S,Cl
!
  write(*,*)''
  write(*,*)'Soluble soil Al content (g Mg-1): CAL'
  write(*,*)(CAL(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Fe content (g Mg-1): CFE'
  write(*,*)(CFE(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Ca content (g Mg-1): CCA'
  write(*,*)(CCA(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Na content (g Mg-1): CNA'
  write(*,*)(CNA(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil K content (g Mg-1): CKA'
  write(*,*)(CKA(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil SO4 content (gS Mg-1): CSO4'
  write(*,*)(CSO4(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Cl content (g Mg-1): CCL'
  write(*,*)(CCL(L,NY,NX),L=NU,NM)
!
!     PRECIPITATED MINERAL CONCENTRATIONS
!
!     CALPO,CFEPO,CCAPD,CCAPH=AlPO4,FePO4,CaHPO4,apatite (g Mg-1)
!     CALOH,CFEOH,CCACO,CCASO=AlOH3,FeOH3,CaSO4,CaCO3 (g Mg-1)
!
  write(*,*)''
  write(*,*)'Soil AlPO4 content (g Mg-1): CALPO'
  write(*,*)(CALPO(L,NY,NX),L=NU,NM)
  write(*,*)'Soil FePO4 content (g Mg-1): CFEPO'
  write(*,*)(CFEPO(L,NY,NX),L=NU,NM)
  write(*,*)'Soil CaHPO4 content (g Mg-1): CCAPD'
  write(*,*)(CCAPD(L,NY,NX),L=NU,NM)
  write(*,*)'Soil apatite content (g Mg-1): CCAPH '
  write(*,*)(CCAPH(L,NY,NX),L=NU,NM)
  write(*,*)'Soil Al(OH)3 content (g Mg-1): CALOH '
  write(*,*)(CALOH(L,NY,NX),L=NU,NM)
  write(*,*)'Soil Fe(OH)3 content (g Mg-1): CFEOH '
  write(*,*)(CFEOH(L,NY,NX),L=NU,NM)
  write(*,*)'Soil CaCO3 content (g Mg-1): CCACO'
  write(*,*)(CCACO(L,NY,NX),L=NU,NM)
  write(*,*)'Soil CaSO4 content (g Mg-1): CCASO'
  write(*,*)(CCASO(L,NY,NX),L=NU,NM)
!
!     GAPON SELECTIVITY CO-EFFICIENTS
!
!     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
!     Ca-NH4,Ca-H,Ca-Al,Ca-Mg,Ca-Na,Ca-K
!
  write(*,*)''
  write(*,*)'Ca-NH4 Gapon selectivity coefficient: GKC4'
  write(*,*)(GKC4(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-H Gapon selectivity coefficient: GKCH'
  write(*,*)(GKCH(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-Al Gapon selectivity coefficient: GKCA'
  write(*,*)(GKCA(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-Mg Gapon selectivity coefficient: GKCM'
  write(*,*)(GKCM(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-Na Gapon selectivity coefficient: GKCN'
  write(*,*)(GKCN(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-K Gapon selectivity coefficient: GKCK'
  write(*,*)(GKCK(L,NY,NX),L=NU,NM)

!
!     INITIAL WATER, ICE CONTENTS
!
  write(*,*)'Initial soil water content (m3/m3): THW'
  write(*,*)(THW(L,NY,NX),L=NU,NM)
  write(*,*)'Initial soil ice content (m3/m3): THI'
  write(*,*)(THI(L,NY,NX),L=NU,NM)
!
!     THW,THI=initial water,ice:>1=satd,1=FC,0=WP,<0=0,0-1=m3 m-3
!
!     INITIAL PLANT AND ANIMAL RESIDUE C, N AND P
!
!     RSC,RSC,RSP=C,N,P in fine(1),woody(0),manure(2) litter (g m-2)
!
  write(*,*)''
  write(*,*)'Initial fine litter C (gC m-2): RSC(k_fine_litr)'
  write(*,*)(RSC(k_fine_litr,L,NY,NX),L=NU,NM)
  write(*,*)'Initial fine litter N (gN m-2): RSN(k_fine_litr)'
  write(*,*)(RSN(k_fine_litr,L,NY,NX),L=NU,NM)
  write(*,*)'Initial fine litter P (gP m-2): RSP(k_fine_litr)'
  write(*,*)(RSP(k_fine_litr,L,NY,NX),L=NU,NM)
  write(*,*)'Initial woody liter C (gC m-2): RSC(k_woody_litr)'
  write(*,*)(RSC(k_woody_litr,L,NY,NX),L=NU,NM)
  write(*,*)'Initial woody litter N (gN m-2): RSN(k_woody_litr)'
  write(*,*)(RSN(k_woody_litr,L,NY,NX),L=NU,NM)
  write(*,*)'Initial woody litter P (gP m-2): RSP(k_woody_litr)'
  write(*,*)(RSP(k_woody_litr,L,NY,NX),L=NU,NM)
  write(*,*)'Initial manure liter C (gC m-2): RSC(k_manure)'
  write(*,*)(RSC(k_manure,L,NY,NX),L=NU,NM)
  write(*,*)'Initial manure litter N (gN m-2): RSN(k_manure)'
  write(*,*)(RSN(k_manure,L,NY,NX),L=NU,NM)
  write(*,*)'Initial manure litter P (gP m-2): RSP(k_manure)'
  write(*,*)(RSP(k_manure,L,NY,NX),L=NU,NM)
  write(*,'(100A)')('=',ll=1,100)

  end associate
  end subroutine Disp_topo_charc

end module readiMod
