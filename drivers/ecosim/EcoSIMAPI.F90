module EcoSIMAPI
  USE EcoSIMCtrlDataType
  use timings,         only: start_timer,        end_timer
  use ErosionMod,      only: erosion
  use Hour1Mod,        only: hour1
  use RedistMod,       only: redist
  use GeochemAPI,      only: soluteModel
  use PlantMod,        only: PlantModel
  use MicBGCAPI,       only: MicrobeModel,       MicAPI_Init,      MicAPI_cleanup
  use TranspNoSaltMod, only: TranspNoSalt
  use TranspSaltMod,   only: TranspSalt
  use SoilBGCDataType, only: trcs_soHml_vr,       trcs_solml_vr
  use TracerIDMod,     only: ids_NO2B,           ids_NO2,          idg_O2
  use PerturbationMod, only: check_Soil_Warming, set_soil_warming, config_soil_warming
  use WatsubMod,       only: watsub
  use PlantMgmtDataType, only: NP  
  use SoilWaterDataType
  use EcoSIMCtrlMod  
  use BalancesMod  
  use SoilDiagsMod  
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: soil
  public :: readnamelist
  public :: regressiontest,write_modelconfig
  logical :: do_timing
contains

  subroutine Run_EcoSIM_one_step(I,J,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8) :: t1

  if(lverb)WRITE(*,334)'HOUR1'
  if(do_timing)call start_timer(t1)
  !print*,'hour1'
  CALL HOUR1(I,J,NHW,NHE,NVN,NVS)

  if(do_timing)call end_timer('HOUR1',t1)
  !
  !   CALCULATE SOIL ENERGY BALANCE, WATER AND HEAT FLUXES IN 'WATSUB'
  !
  !print*,'watsub'
  if(lverb)WRITE(*,334)'WAT'
  if(do_timing)call start_timer(t1)
  CALL WATSUB(I,J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('WAT',t1)
  !
  !   CALCULATE SOIL BIOLOGICAL TRANSFORMATIONS IN 'NITRO'
  !     
  if(microbial_model)then
    if(lverb)WRITE(*,334)'NIT'
    if(do_timing)call start_timer(t1)
    CALL MicrobeModel(I,J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('NIT',t1)
  endif
!  print*,'plant model'
  !
  !   UPDATE PLANT biogeochemistry
  !
  if(lverb)WRITE(*,334)'PlantModel'
  if(plant_model)then
    if(do_timing)call start_timer(t1)  
    call PlantModel(I,J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('PlantModel',t1)    
  endif

  !
  !   CALCULATE SOLUTE EQUILIBRIA IN 'SOLUTE'
  !
  if(soichem_model)then
    if(lverb)WRITE(*,334)'SOL'
    if(do_timing)call start_timer(t1)
    CALL soluteModel(I,J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('soluteModel',t1)
  endif

  !
  !   CALCULATE GAS AND SOLUTE FLUXES IN 'TranspNoSalt'
  !
  if(lverb)WRITE(*,334)'TRN'
  if(do_timing)call start_timer(t1)
  CALL TranspNoSalt(I,J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('TranspNoSalt',t1)

  !
  !   CALCULATE ADDITIONAL SOLUTE FLUXES IN 'TranspSalt' IF SALT OPTION SELECTED
  !
  if(lverb)WRITE(*,334)'TRNS'

  if(salt_model)then
    if(do_timing)call start_timer(t1)
    CALL TranspSalt(I,J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('TranspSalt',t1)
  endif

  !
  !   CALCULATE SOIL SEDIMENT TRANSPORT IN 'EROSION'
  !
  if(lverb)WRITE(*,334)'EROSION'

  if(do_timing)call start_timer(t1)
  CALL EROSION(I,J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('EROSION',t1)
  !
  !   UPDATE ALL SOIL STATE VARIABLES FOR WATER, HEAT, GAS, SOLUTE
  !   AND SEDIMENT FLUXES IN 'REDIST'
  !
  if(lverb)WRITE(*,334)'REDIST'
  if(do_timing)call start_timer(t1)    
  CALL REDIST(I,J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('REDIST',t1)
334   FORMAT(A8)

  call DiagSoilGasPressure(I,J,NHW,NHE,NVN,NVS)    

  call EndCheckBalances(I,J,NHW,NHE,NVN,NVS)

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
  character(len=*), parameter :: mod_filename = &
  __FILE__
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
    clm_hour_file_in,clm_day_file_in,soil_mgmt_in,forc_periods,NCYC_LITR,NCYC_SNOW,&
    NPXS,NPYS,continue_run,visual_out,restart_out,&
    finidat,restartFileFullPath,brnch_retain_casename,plant_model,microbial_model,&
    soichem_model,atm_ghg_in,aco2_ppm,ao2_ppm,an2_ppm,an2_ppm,ach4_ppm,anh3_ppm,&
    snowRedist_model,disp_planttrait,iErosionMode,grid_mode,atm_ch4_fix,atm_n2o_fix,&
    atm_co2_fix,first_topou,first_pft,fixWaterLevel,arg_ppm,ldebug_day

  namelist /ecosim/hist_nhtfrq,hist_mfilt,hist_fincl1,hist_fincl2,hist_yrclose, &
    do_budgets,ref_date,start_date,do_timing,warming_exp

  logical :: laddband
  namelist /bbgcforc/do_bgcforc_write,do_year,do_doy,laddband,do_layer,&
    bgc_fname

  !local variables
  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
  integer :: year0

  do_timing    = .false.
  continue_run = .false.
  NPXS         = 30   !number of cycles per hour for water, heat, solute flux calcns
  NPYS         = 20   !number of cycles per NPX for gas flux calculations
  
  ldebug_day  =-1
  NCYC_LITR             = 20
  NCYC_SNOW             = 20
  grid_mode             = 3
  iErosionMode          = -1
  visual_out            = .false.
  restart_out           = .false.
  do_budgets            = .false.
  plant_model           = .true.
  soichem_model         = .true.
  microbial_model       = .true.
  disp_planttrait       = .false.
  ref_date              = '18000101000000'   !place holder for future
  start_date            = '18000101000000'   !start date of the simulation, differ from the forcing date
  finidat               = ''
  restartFileFullPath   = ''                 !contains the full pathname of the restart file
  warming_exp           = ''
  brnch_retain_casename = .false.
  hist_yrclose          = .false.
  fixWaterLevel         = .false.
  forc_periods      = 0
  forc_periods(1:9) = (/1980,1980,1,1981,1988,2,1989,2008,1/)

  num_of_simdays       = -1
  do_year              = -1
  do_doy               = 0
  do_layer             = 1
  salt_model           = .false.
  laddband             = .false.
  do_regression_test   = .false.
  lverbose             = .false.
  num_microbial_guilds = 1
  do_bgcforc_write     = .false.
  bgc_fname            = 'bbforc.nc'
  do_instequil         = .false.

  clm_factor_in    = 'NO'
  pft_file_in      = ''
  grid_file_in     = ''
  pft_mgmt_in      = 'NO'
  clm_hour_file_in = ''
  clm_day_file_in  = ''
  soil_mgmt_in     = 'NO'
  atm_ghg_in       = ''
  aco2_ppm         = 280._r8
  ach4_ppm         = 1.144_r8
  an2o_ppm         = 0.270_r8
  ao2_ppm          = 0.209e6_r8
  arg_ppm          = 0.00934e6_r8
  an2_ppm          = 0.78e6_r8
  anh3_ppm         = 5.e-3_r8
  atm_co2_fix      = -100._r8
  atm_n2o_fix      = -100._r8
  atm_ch4_fix      = -100._r8
  first_topou      = .false.
  
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
  erosion_model=iErosionMode>=0
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
  if(.not.soichem_model)then
    salt_model=.false.
  endif

  call config_soil_warming(warming_exp)
end subroutine readnamelist
! ----------------------------------------------------------------------

subroutine soil(NHW,NHE,NVN,NVS,nlend)
!!
! Description:
! THIS IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod,   only: r8 => DAT_KIND_R8
  use DayMod,          only: day
  use ExecMod,         only: exec
  use StarteMod,       only: starte
  use StartqMod,       only: startq
  use StartsMod,       only: starts
  use VisualMod,       only: visual
  use WthrMod,         only: wthr
  use RestartMod,      only: restFile
  use PlantInfoMod,    only: ReadPlantInfo
  use readsmod,        only: ReadClimSoilForcing
  use YearMod,         only: SetAnnualAccumlators
  use timings,         only: init_timer,         start_timer, end_timer, end_timer_loop
  use ForcWriterMod,   only: do_bgcforc_write,   WriteBBGCForc
  use HistDataType,    only: hist_ecosim
  use HistFileMod,     only: hist_update_hbuf,   hist_htapes_wrapup
  use ClimReadMod,     only: read_soil_warming_Tref
  use EcoSIMCtrlMod
  use GridConsts
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMConfig

  implicit none
  integer :: yearc, yeari
  integer, intent(in) :: NHW,NHE,NVN,NVS
  logical, intent(out) :: nlend
  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8) :: t1
  integer :: I,J
  integer :: idaz
  character(len=14) :: ymdhs
  logical :: rstwr, lnyr
  logical :: lverb0
! begin_execution
!
! READ INPUT DATA FOR SITE, SOILS AND MANAGEMENT IN 'READS'
! AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTS'
!
333   FORMAT(A8)

  is_first_year=frectyp%yearacc.EQ.0

  if(do_timing)call init_timer(outdir)

  if(lverb)WRITE(*,333)'READS: read climate forcing'
  CALL ReadClimSoilForcing(frectyp%yearcur,frectyp%yearclm,NHW,NHE,NVN,NVS)

  !temporary set up for setting mass balance check
  IBEGIN=1;ISTART=1;ILAST=0

  call etimer%get_ymdhs(ymdhs)

  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4))THEN
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
    call ReadPlantInfo(frectyp%yearcur,frectyp%yearclm,NHW,NHE,NVN,NVS)
  endif

! INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
!
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4) .and. plant_model)THEN
    !initialize by year
    if(lverb)WRITE(*,333)'STARTQ'
    CALL STARTQ(NHW,NHE,NVN,NVS,1,JP)
  ENDIF
  if(ymdhs(1:4)==frectyp%ymdhs0(1:4) .and. soichem_model)then
! INITIALIZE ALL SOIL CHEMISTRY VARIABLES IN 'STARTE'
! This is done done every year, because tracer concentrations
! in rainfall vary every year. In a more reasonable way, e.g.,
! when coupled to atmospheric chemistry code, it should be done by
! hour
    if(lverb)WRITE(*,333)'STARTE'
    CALL STARTE(NHW,NHE,NVN,NVS)
  endif

  iYearCurrent=frectyp%yearcur

  if(check_Soil_Warming(iYearCurrent,1))then
    call read_soil_warming_Tref(iYearCurrent,NHW,NHE,NVN,NVS)    

    call set_soil_warming(iYearCurrent,NHW,NHE,NVN,NVS)
  endif
  lverb0 = lverb
  DazCurrYear=etimer%get_days_cur_year()
  DO I=1,DazCurrYear    
    if(ldebug_day==I)then
      lverb=.true.      
    else
      lverb=lverb0
    endif  
    IF(do_rgres .and. I.eq.LYRG)RETURN
    !   UPDATE DAILY VARIABLES SUCH AS MANAGEMENT INPUTS
    !
    if(lverb)WRITE(*,333)'DAY'
    CALL DAY(I,NHW,NHE,NVN,NVS)
    
    call SetAnnualAccumlators(I, NHW, NHE, NVN, NVS)
    
    DO J=1,24
      call etimer%get_ymdhs(ymdhs)
      
      if(ymdhs==frectyp%ymdhs0)then
        frectyp%lskip_loop=.false.        
        if(is_restart() .or. is_branch())then
          call restFile(flag='read')
          if (j==1)call SetAnnualAccumlators(I, NHW, NHE, NVN, NVS)

          call SummarizeTracerMass(I,J,NHW,NHE,NVN,NVS)
          
        endif
      endif
      if(frectyp%lskip_loop)then
        call etimer%update_time_stamp()
        cycle
      endif
    !
    !   UPDATE HOURLY VARIABLES IN 'HOUR1'
    !   set up climate forcing for the new hour

      if(do_timing)call start_timer(t1)
      CALL WTHR(I,J,NHW,NHE,NVN,NVS)
      if(do_timing)call end_timer('WTHR',t1)

      if(lverb)WRITE(*,333)'Run_EcoSIM_one_step'
      call Run_EcoSIM_one_step(I,J,NHW,NHE,NVN,NVS)

      if(do_timing)call end_timer_loop()
      
      if(lverb)write(*,*)'hist_update'
      call hist_ecosim%hist_update(I,J,bounds)

      call hist_update_hbuf(bounds)

      call etimer%update_time_stamp()

      nlend=etimer%its_time_to_exit()
      rstwr=etimer%its_time_to_write_restart(nlend)
      lnyr=etimer%its_a_new_year()

      call hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )

      if(rstwr)then
        call restFile(flag='write')
      endif      
      if(nlend)exit
    END DO

!
! PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
!
    if(frectyp%lskip_loop)cycle

    if(do_bgcforc_write)then
      call WriteBBGCFORC(I,IYRR)
    endif
    if(nlend)exit
  END DO
  if(lverb)write(*,333)'exit soil'
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

  character(len=*), parameter :: mod_filename = &
  __FILE__
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
      IF(IsPlantActive_pft(NZ,NY,NX).EQ.iActive)THEN

        category = 'flux'
        name = 'NH4_UPTK (g m^-3 h^-1)'
        datv=0._r8
        do ll=1,12
          if(AREA(3,ll,NY,NX)>0._r8)then
            datv(ll)=safe_adb(RootNutUptake_pvr(ids_NH4,1,ll,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4,2,ll,NZ,NY,NX) &
              +RootNutUptake_pvr(ids_NH4B,1,ll,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4B,2,ll,NZ,NY,NX),AREA(3,ll,NY,NX))
          endif
        enddo
        call regression%writedata(category,name, datv)
        exit
      endif

    enddo

    category = 'state'
    name = 'aqueous soil O2 (g m^3)'
    datv=trc_solcl_vr(idg_O2,1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'liquid soil water (m^3 m^-3)'
    datv=ThetaH2OZ_vr(1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'soil temperature (oC)'
    datv=TCS_vr(1:12,NY,NX)
    call regression%writedata(category,name,datv)
    write(*,*)'Finish regression file writing'
    call regression%CloseOutput()
  endif
end subroutine regressiontest

!------------------------------------------------------------------------------------------
  subroutine write_modelconfig()
  use fileUtil, only : getavu, relavu, opnfil
  use readiMod
  implicit none
  character(len=64) :: fnm_loc
  integer :: nu_cfg

  if(disp_modelconfig)then
    nu_cfg=getavu()
    write(fnm_loc,'(A,I4,A)')'ecosim.setup'
    call opnfil(fnm_loc,nu_cfg,'f')    

    write(nu_cfg,*)'microbial_model=',microbial_model
    write(nu_cfg,*)'plant_model=',plant_model
    write(nu_cfg,*)'salt_model=',salt_model
    write(nu_cfg,*)'erosion_model=',erosion_model_status(iErosionMode)
    write(nu_cfg,*)'grid_mode=',GridConectionMode(grid_mode)
    call relavu(nu_cfg)
  endif
  end subroutine write_modelconfig

end module EcoSIMAPI
