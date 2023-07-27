module EcoSIMAPI
  USE EcoSIMCtrlDataType
  use timings      , only : start_timer, end_timer
  use ErosionMod   , only : erosion
  use Hour1Mod     , only : hour1
  use RedistMod    , only : redist
  use GeochemAPI   , only : soluteModel
  use PlantAPI     , only : PlantModel
  use MicBGCAPI    , only : MicrobeModel, MicAPI_Init, MicAPI_cleanup
  use TrnsfrMod    , only : trnsfr
  use TrnsfrsMod   , only : trnsfrs
  use EcoSIMCtrlMod, only : lverb,plant_model,soichem_model,micb_model,salt_model
  use WatsubMod    , only : watsub
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__
  public :: soil
  public :: readnamelist
  public :: regressiontest
contains

  subroutine Run_EcoSIM_one_step(I,J,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8) :: t1

  if(lverb)WRITE(*,334)'HOUR1'
  call start_timer(t1)
  CALL HOUR1(I,J,NHW,NHE,NVN,NVS)
  call end_timer('HOUR1',t1)
  !
  !   CALCULATE SOIL ENERGY BALANCE, WATER AND HEAT FLUXES IN 'WATSUB'
  !
  if(lverb)WRITE(*,334)'WAT'
  call start_timer(t1)
  CALL WATSUB(I,J,NHW,NHE,NVN,NVS)
  call end_timer('WAT',t1)
  !
  !   CALCULATE SOIL BIOLOGICAL TRANSFORMATIONS IN 'NITRO'
  !
  if(micb_model)then
    if(lverb)WRITE(*,334)'NIT'
    call start_timer(t1)
    CALL MicrobeModel(I,J,NHW,NHE,NVN,NVS)
    call end_timer('NIT',t1)
  endif
  !
  !   UPDATE PLANT biogeochemistry
  !
  if(lverb)WRITE(*,334)'PlantModel'
  if(plant_model)then
    call PlantModel(I,J,NHW,NHE,NVN,NVS)
  endif
  !
  !
  !   CALCULATE SOLUTE EQUILIBRIA IN 'SOLUTE'
  !
  
  if(soichem_model)then
    if(lverb)WRITE(*,334)'SOL'  
    call start_timer(t1)
    CALL soluteModel(I,J,NHW,NHE,NVN,NVS)
    call end_timer('SOL',t1)
  endif
  !
  !   CALCULATE GAS AND SOLUTE FLUXES IN 'TRNSFR'
  !
  if(lverb)WRITE(*,334)'TRN'
  call start_timer(t1)
  CALL TRNSFR(I,J,NHW,NHE,NVN,NVS)
  call end_timer('TRN',t1)
  !
  !   CALCULATE ADDITIONAL SOLUTE FLUXES IN 'TRNSFRS' IF SALT OPTION SELECTED
  !
  if(lverb)WRITE(*,334)'TRNS'
  !    if(I>=170)print*,TKS(0,NVN,NHW)'
  if(salt_model)then
    call start_timer(t1)
    CALL TRNSFRS(I,J,NHW,NHE,NVN,NVS)
    call end_timer('TRNSFRS',t1)
  endif
  !
  !   CALCULATE SOIL SEDIMENT TRANSPORT IN 'EROSION'
  !
  if(lverb)WRITE(*,334)'EROSION'
  !    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL EROSION(I,J,NHW,NHE,NVN,NVS)
  call end_timer('EROSION',t1)
  !
  !   UPDATE ALL SOIL STATE VARIABLES FOR WATER, HEAT, GAS, SOLUTE
  !   AND SEDIMENT FLUXES IN 'REDIST'
  !
  if(lverb)WRITE(*,334)'RED'
  !    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL REDIST(I,J,NHW,NHE,NVN,NVS)
  call end_timer('RED',t1)

334   FORMAT(A8)

  end subroutine Run_EcoSIM_one_step
! ----------------------------------------------------------------------
 subroutine readnamelist(nml_buffer, case_name, prefix,LYRG,nmicbguilds)
!!
! Description:
! read control namelist
  use abortutils     , only : endrun
  use EcoSIMConfig   , only : transport_on,column_mode, do_instequil,ref_date
  use EcoSIMConfig   , only : finidat,restartFileFullPath,brnch_retain_casename,start_date
  use ForcWriterMod  , only : bgc_forc_conf,do_bgcforc_write
  use fileUtil       , only : iulog,ecosim_namelist_buffer_size
  use EcoSIMHistMod  , only : DATAC
  use EcoSIMCtrlMod
  use HistFileMod
  implicit none
  character(len=*), parameter :: mod_filename = __FILE__
  character(len=ecosim_namelist_buffer_size), intent(in) :: nml_buffer
  character(len=36)    , intent(out) :: case_name
  character(len=80)    , intent(out) :: prefix
  integer              , intent(out) :: LYRG
  integer              , intent(out) :: nmicbguilds

  logical :: do_regression_test
  integer :: num_of_simdays
  logical :: lverbose
  integer :: num_microbial_guilds
  integer :: do_doy,do_year,do_layer
  character(len=64) :: bgc_fname

  namelist /ecosim/case_name, prefix, do_regression_test, &
    num_of_simdays,lverbose,num_microbial_guilds,transport_on,column_mode,&
    do_instequil,salt_model, pft_file_in,grid_file_in,pft_mgmt_in, clm_factor_in,&
    clm_file_in,soil_mgmt_in,forc_periods,&
    NPXS,NPYS,JOUTS,continue_run,visual_out,restart_out,&
    finidat,restartFileFullPath,brnch_retain_casename,plant_model,micb_model,&
    soichem_model,atm_ghg_in,aco2_ppm,ao2_ppm,an2_ppm,an2_ppm,ach4_ppm,anh3_ppm
  
  namelist /ecosim/hist_nhtfrq,hist_mfilt,hist_fincl1,hist_fincl2,hist_yrclose, &
    do_budgets,ref_date,start_date

  logical :: laddband
  namelist /bbgcforc/do_bgcforc_write,do_year,do_doy,laddband,do_layer,&
    bgc_fname

  !local variables
  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
  integer :: year0

  continue_run=.false.
  NPXS=30   !number of cycles per hour for water,heat,solute flux calcns
  NPYS=20   !number of cycles per NPX for gas flux calcns
  JOUTS=1   !frequency on hourly scale

  visual_out =.false.
  restart_out=.false.  
  do_budgets =.false.
  plant_model=.true.
  soichem_model=.true.
  micb_model=.true.
  ref_date  = '18000101000000'   !place holder for future
  start_date= '18000101000000'   !start date of the simulation, differ from the forcing date
  finidat=' '
  restartFileFullPath = ' '                   !contains the full pathname of the restart file

  brnch_retain_casename=.false.
  hist_yrclose=.false.
  
  forc_periods=(/1980,1980,1,1981,1988,2,1989,2008,1/)

  num_of_simdays=-1
  do_year=-1
  do_doy=0
  do_layer=1
  salt_model=.false.
  laddband=.false.
  do_regression_test=.false.
  lverbose=.false.
  num_microbial_guilds=1
  do_bgcforc_write=.false.
  bgc_fname='bbforc.nc'
  do_instequil=.false.

  clm_factor_in=''
  pft_file_in=''
  grid_file_in=''
  pft_mgmt_in=''
  clm_file_in=''
  soil_mgmt_in=''
  atm_ghg_in=''
  aco2_ppm  = 280._r8
  ach4_ppm  = 1.144_r8
  an2o_ppm  = 0.270_r8
  ao2_ppm   = 0.209e6_r8
  an2_ppm   = 0.78e6_r8
  anh3_ppm  = 5.e-3_r8  

  read(nml_buffer, nml=ecosim, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading ecosim namelist ",nml_error,ioerror_msg
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  read(nml_buffer, nml=bbgcforc, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading bbgcforc namelist "
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,ecosim)
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,bbgcforc)
    write(iulog, *)
    write(iulog, *) '--------------------'

  endif
  read(start_date,'(I4)')year0
  call etimer%Init(nml_buffer,year0=year0)
  if(do_bgcforc_write)then
    bgc_forc_conf%doy =do_doy
    bgc_forc_conf%year=do_year
    bgc_forc_conf%laddband=laddband
    bgc_forc_conf%Layer=do_layer
    bgc_forc_conf%bgc_fname=bgc_fname
  endif
  do_rgres=do_regression_test
  LYRG=num_of_simdays
  lverb=lverbose
  nmicbguilds=num_microbial_guilds
  if(.not. soichem_model)then
    salt_model=.false.
  endif
end subroutine readnamelist
! ----------------------------------------------------------------------

subroutine soil(NE,NEX,NHW,NHE,NVN,NVS,nlend)
!!
! Description:
! THIS IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use DayMod       , only : day
  use ExecMod      , only : exec
  use StarteMod    , only : starte
  use StartqMod    , only : startq
  use StartsMod    , only : starts
  use VisualMod    , only : visual
  use WthrMod      , only : wthr
  use RestartMod   , only : restFile
  use PlantInfoMod , only : ReadPlantInfo
  use readsmod     , only : reads
  use timings      , only : init_timer, start_timer, end_timer,end_timer_loop
  use InitEcoSIM   , only : InitModules2
  use EcoSIMCtrlMod
  use PlantAPI     , only : PlantModel
  use MicBGCAPI    , only : MicrobeModel, MicAPI_Init, MicAPI_cleanup
  use ForcWriterMod, only : do_bgcforc_write,WriteBBGCForc
  use GridConsts
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use HistDataType , only : hist_ecosim
  use HistFileMod  , only : hist_update_hbuf,hist_htapes_wrapup

  implicit none
  integer :: yearc, yeari
  integer, intent(in) :: NE,NEX,NHW,NHE,NVN,NVS
  logical, intent(out) :: nlend
  character(len=*), parameter :: mod_filename = __FILE__
  real(r8) :: t1
  integer :: I,J
  integer :: idaz
  character(len=14) :: ymdhs
  logical :: rstwr, lnyr
! begin_execution
!
! READ INPUT DATA FOR SITE, SOILS AND MANAGEMENT IN 'READS'
! AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTS'
!
333   FORMAT(A8)

  is_first_year=frectyp%yearacc.EQ.0

  call init_timer(outdir)

  if(lverb)WRITE(*,333)'READS: read climate forcing'
  CALL READS(frectyp%yearcur,frectyp%yearclm,NE,NEX,NHW,NHE,NVN,NVS)

  !temporary set up for setting mass balance check
  IBEGIN=1;ISTART=1;ILAST=0
  
  call etimer%get_ymdhs(ymdhs)
  
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4))THEN
!   the first simulation year
    if(lverb)WRITE(*,333)'STARTS'
    CALL STARTS(NHW,NHE,NVN,NVS)
!
!   RECOVER VALUES OF ALL SOIL STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTS' IF NEEDED

  ENDIF
!
  if(plant_model)then
    !plant information is read in every year, but the active flags
    !are set using the checkpoint file.
    if(lverb)WRITE(*,333)'ReadPlantInfo'  
    WRITE(*,333)'ReadPlantInfo'  
    call ReadPlantInfo(frectyp%yearcur,frectyp%yearclm,NE,NEX,NHW,NHE,NVN,NVS)
  endif

! INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
!
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4) .and. plant_model)THEN
    !initialize by year
    if(lverb)WRITE(*,333)'STARTQ'
    CALL STARTQ(NHW,NHE,NVN,NVS,1,JP)

  ENDIF

  if(soichem_model)then
! INITIALIZE ALL SOIL CHEMISTRY VARIABLES IN 'STARTE' 
! This is done done every year, because tracer concentrations
! in rainfall vary every year. In a more reasonable way, e.g., 
! when coupled to atmospheric chemistry code, it should be done by 
! hour
    if(lverb)WRITE(*,333)'STARTE'
    CALL STARTE(NHW,NHE,NVN,NVS)
  endif

  iyear_cur=frectyp%yearcur
  LYRC=etimer%get_days_cur_year()
  
  DO I=1,LYRC
    IF(do_rgres .and. I.eq.LYRG)RETURN
    !   UPDATE DAILY VARIABLES SUCH AS MANAGEMENT INPUTS
    !
    if(lverb)WRITE(*,333)'DAY'
    CALL DAY(I,NHW,NHE,NVN,NVS)

    DO J=1,24
      call etimer%get_ymdhs(ymdhs)
      
      if(ymdhs==frectyp%ymdhs0)then
        frectyp%lskip_loop=.false. 
        if(is_restart())then          
          call restFile(flag='read')
        endif
      endif        
      if(frectyp%lskip_loop)then
        call etimer%update_time_stamp()      
        cycle
      endif

    !
    !   UPDATE HOURLY VARIABLES IN 'HOUR1'
    !   set up climate forcing for the new hour
      if(lverb)WRITE(*,333)'WTHR'
      call start_timer(t1)
      CALL WTHR(I,J,NHW,NHE,NVN,NVS)
      call end_timer('WTHR',t1)
      
      if(lverb)WRITE(*,333)'Run_EcoSIM_one_step'
      call Run_EcoSIM_one_step(I,J,NHW,NHE,NVN,NVS)
  !
  !   WRITE OUTPUT FOR DYNAMIC VISUALIZATION
  !
      IF(visual_out)THEN
        IF((J/JOUT)*JOUT.EQ.J)THEN
          if(lverb)WRITE(*,333)'VIS'
          CALL VISUAL(I,J,NHW,NHE,NVN,NVS)
        ENDIF
      ENDIF
          
      call end_timer_loop()
      
      call hist_ecosim%hist_update(bounds)
      
      call hist_update_hbuf(bounds)      
      call etimer%update_time_stamp()      

      nlend=etimer%its_time_to_exit()
      rstwr=etimer%its_time_to_write_restart()
      lnyr=etimer%its_a_new_year().and.hist_yrclose
      
      call hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )      
      if(rstwr)then
        call restFile(flag='write')
      endif
      if(nlend)exit      
    END DO
    
!
! PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
!
    if(lverb)WRITE(*,333)'EXEC'
    CALL EXEC(I)
!

    if(do_bgcforc_write)then
      call WriteBBGCForc(I,IYRR)
    endif
    if(nlend)exit
  END DO

  RETURN
END subroutine soil
! ----------------------------------------------------------------------

subroutine regressiontest(nmfile,case_name, NX, NY)
!!
! Description:
! subroutine to conduct regression tests

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use TestMod      , only : regression
  use minimathmod  , only : safe_adb
  use GridConsts
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use PlantDataRateType
  use SoilBGCDataType
  use GridDataType
  implicit none

  character(len=*), parameter :: mod_filename = __FILE__
  character(len=*), intent(in) :: nmfile
  character(len=*), intent(in) :: case_name
  integer, intent(in) :: NX
  integer, intent(in) :: NY
  !local variables
  character(len=128) :: category
  character(len=128) :: name
  integer :: nz, ll
  real(r8) :: datv(12)

  call regression%Init(trim(nmfile),case_name)

  if(regression%write_regression_output)then
    write(*,*)'write regression file'
    call regression%OpenOutput()

    do NZ=1,NP(NY,NX)
      IF(IFLGC(NZ,NY,NX).EQ.PlantIsActive)THEN

        category = 'flux'
        name = 'NH4_UPTK (g m^-3 h^-1)'
        datv=0._r8
        do ll=1,12
          if(AREA(3,ll,NY,NX)>0._r8)then
            datv(ll)=safe_adb(RUPNH4(1,ll,NZ,NY,NX)+RUPNH4(2,ll,NZ,NY,NX) &
              +RUPNHB(1,ll,NZ,NY,NX)+RUPNHB(2,ll,NZ,NY,NX),AREA(3,ll,NY,NX))
          endif
        enddo
        call regression%writedata(category,name, datv)
        exit
      endif

    enddo

    category = 'state'
    name = 'aqueous soil O2 (g m^3)'
    datv=trc_solcl(idg_O2,1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'liquid soil water (m^3 m^-3)'
    datv=THETWZ(1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'soil temperature (oC)'
    datv=TCS(1:12,NY,NX)
    call regression%writedata(category,name,datv)
    write(*,*)'Finish regression file writing'
    call regression%CloseOutput()
  endif
end subroutine regressiontest

end module EcoSIMAPI
