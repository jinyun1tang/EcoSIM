module EcoSIMCtrlMod
  use ncdio_pio, only : file_desc_t
  use ecosim_Time_Mod, only : ecosim_time_type
implicit none
  save
  logical :: salt_model=.false.
  logical :: erosion_model=.false.
  character(len=300) :: pft_file_in
  character(len=300) :: pft_mgmt_in
  type(file_desc_t) :: pft_nfid
  character(len=300) :: grid_file_in
  character(len=300) :: clm_file_in
  type(ecosim_time_type) :: etimer
! define the simulation controler
! sim_periods(:,1)=(/year1,year2,1/)     cold start
! sim_periods(:,2)=(/year1,year2,ncycs/) spin up
! sim_periods(:,3)=(/year1,year2,1/)     regular sims
  integer :: sim_periods(3,3)


end module EcoSIMCtrlMod
