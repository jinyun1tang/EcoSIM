module EcoSIMCtrlMod
  use ncdio_pio,       only: file_desc_t
  use data_kind_mod,   only: r8 => DAT_KIND_R8
  use ecosim_Time_Mod, only: ecosim_time_type
  use fileUtil,        only: datestrlen, iulog
  use abortutils,      only: endrun
implicit none
  save
  character(len=*),private, parameter :: mod_filename =&
   __FILE__  
  integer :: ldebug_day  = -1 
  logical :: salt_model   =.false.    !toggle for salt model
  logical :: erosion_model=.false.
  logical :: first_topou   = .false.  !only simulate first topo unit
  logical :: first_pft     = .false.  !only simulate first pft
  integer :: iErosionMode =-1         !erosion option  
  logical :: plant_model       =.true.
  logical :: microbial_model   =.true.
  logical :: soichem_model     =.true.
  logical :: snowRedist_model  =.true.
  logical :: ATS_cpl_mode      =.false.
  logical :: plantOM4Heat      =.false.
  logical :: fixWaterLevel     =.false.
  integer :: yearf1       !first year of daily climate forcing
  integer :: yearf2       !first year of hourly climate forcing
  integer :: nyeardal1    !number of daily climate forcing
  real(r8) :: aco2_ppm  = 280._r8
  real(r8) :: ach4_ppm  = 1.144_r8
  real(r8) :: an2o_ppm  = 0.270_r8
  real(r8) :: ao2_ppm   = 0.209e6_r8
  real(r8) :: arg_ppm   = 0.00934e6_r8
  real(r8) :: an2_ppm   = 0.78e6_r8
  real(r8) :: anh3_ppm  = 5.e-3_r8
  real(r8) :: atm_co2_fix=-100._r8
  real(r8) :: atm_ch4_fix=-100._r8
  real(r8) :: atm_n2o_fix=-100._r8
  character(len=256) :: warming_exp=''
  character(len=300) :: pft_file_in
  character(len=300) :: pft_mgmt_in
  character(len=300) :: grid_file_in
  character(len=300) :: clm_hour_file_in =''     !file for hourly climate forcing
  character(len=300) :: clm_day_file_in  =''    !file for daily climate forcing
  character(len=300) :: soil_mgmt_in     !file for soil management information
  character(len=300) :: clm_factor_in    !file for climate change factors
  character(len=300) :: atm_ghg_in       !file for atmospheric GHG concentrations
  logical :: do_budgets = .false.
  type(file_desc_t)  :: pft_nfid 
  type(ecosim_time_type) :: etimer
  logical :: Lirri_auto=.false.

  character(len=8) :: diag_opt ='nsteps'
  integer :: diag_frq=-999999999
  logical :: continue_run
  logical :: restart_out
  logical :: visual_out
  logical :: hist_yrclose
  character(len=16) :: hist_config(10)
  character(len=8)  :: sim_yyyymmdd
  integer :: forc_periods(15)
  integer :: NPXS(5),NPYS(5)
  integer :: NCYC_LITR  !number of subcycles for litr
  integer :: NCYC_SNOW  !number of subcycles for snow
  logical :: lverb           !logical switch for verbose output
  logical :: disp_planttrait =.true.
  logical :: disp_modelconfig=.true.
  logical :: do_rgres        !logical switch for regression tests
  integer :: grid_mode = 3  !vertical only

  type, public :: forc_data_rec_type
  integer :: pft_rec
  integer :: yearclm
  integer :: yearacc      !year passed
  integer :: yearcur      !current year 
  integer :: yearpre      !previous year
  integer :: yearrst      !restart year
  logical :: lskip_loop   !logical switch to skip a loop
  character(len=datestrlen) :: ymdhs0  !the beginning yyyymmddhhmmss
  contains
  procedure, public :: Init  => init_frectyp
  end type forc_data_rec_type

  type(forc_data_rec_type), public :: frectyp

  contains

  subroutine Init_frectyp(this)
  implicit none
  class(forc_data_rec_type)  :: this

  this%pft_rec=0    !rec index for pft data
  this%yearclm=0    !year index for climate data
  this%yearacc=0    !accumulated years for the whole simulation
  this%yearcur=0    !current year
  this%yearpre=0    !previous year
  this%yearrst=0    !restarting year
  this%lskip_loop=.true.
  this%ymdhs0='00000000000000'
  end subroutine Init_frectyp

  !-----------------------------------------------------------------------
  integer function get_sim_len(forc_periods,nperiods)
  implicit none
  integer, dimension(:), intent(in) :: forc_periods
  integer :: nperiods
  integer :: nn1,id,nelms

  nelms = size(forc_periods)
  nperiods   = nelms/3

  get_sim_len=0
  DO nn1=0,nperiods-1
    id=nn1*3+1
    get_sim_len=get_sim_len+(abs(forc_periods(id+1)-forc_periods(id))+1)*forc_periods(id+2)
  enddo
  if(get_sim_len<0)then
  call endrun('Negative simulation length, check forc_periods set up in '//mod_filename,__LINE__)
  endif
  end function get_sim_len  


end module EcoSIMCtrlMod
