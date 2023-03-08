module EcoSIMCtrlMod
  use ncdio_pio, only : file_desc_t
  use ecosim_Time_Mod, only : ecosim_time_type
implicit none
  save
  logical :: salt_model=.false.
  logical :: erosion_model=.false.
  character(len=300) :: pft_file_in
  character(len=300) :: pft_mgmt_in
  character(len=300) :: grid_file_in
  character(len=300) :: clm_file_in
  character(len=300) :: soil_mgmt_in
  character(len=300) :: clm_factor_in          !file for climate change factors
  type(file_desc_t)  :: pft_nfid 
  type(ecosim_time_type) :: etimer
  logical :: Lirri_auto=.false.

  logical :: continue_run
  logical :: restart_out
  logical :: visual_out
  logical :: cold_run
  character(len=16) :: hist_config(10)
  character(len=8)  :: sim_yyyymmdd
  integer :: forc_periods(9)
  integer :: NPXS(3),NPYS(3),JOUTS(3),IOUTS(3),KOUTS(3)
  logical :: lverb           !logical switch for verbose output
  logical :: do_rgres        !logical switch for regression tests

  type, public :: forc_data_rec_type
  integer :: pft_rec
  integer :: yearclm
  integer :: yearacc      !year passed
  integer :: yearcur      !current year 
  integer :: yearpre      !previous year
  logical :: lskip_loop   !logical switch to skip a loop
  character(len=14) :: ymdhs0  !the beginning yyyymmddhhmmss
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
  this%lskip_loop=.true.
  this%ymdhs0='00000000000000'
  end subroutine Init_frectyp
end module EcoSIMCtrlMod
