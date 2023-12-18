program GridReader
  use ncdio_pio
  use netcdf
  use abortutils, only : endrun
implicit none
  type(file_desc_t) :: grid_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: NHW,NHE,NVN,NVS,ngrid,ntopou
  character(len=120) :: ncd_file
  character(len=*), parameter :: mod_filename = &
  __FILE__

  CALL GETARG(1,ncd_file)

  call ncd_pio_openfile(grid_nfid, ncd_file, ncd_nowrite)
  write(*,*)'open file ok',ncd_file
  ngrid=get_dim_len(grid_nfid, 'ngrid')

  ntopou=get_dim_len(grid_nfid, 'ntopou')

  write(*,*)ngrid
  call check_var(grid_nfid, 'NHW', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NHW in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NHW), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NHE', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NHE in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NHE), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NVN', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NVN in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NVN), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NVS', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NVS in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NVS), 'in '//trim(mod_filename))

  call ncd_pio_closefile(grid_nfid)




end program GridReader