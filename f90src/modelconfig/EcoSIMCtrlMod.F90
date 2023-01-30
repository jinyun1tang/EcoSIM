module EcoSIMCtrlMod
  use ncdio_pio, only : file_desc_t
implicit none
  save
  logical :: salt_model=.false.
  logical :: erosion_model=.false.
  character(len=300) :: pft_file_in
  character(len=300) :: pft_mgmt_in
  type(file_desc_t) :: pft_nfid
  character(len=300) :: grid_file_in

end module EcoSIMCtrlMod
