program main
  use abortutils, only : endrun
  use fileUtil, only : ecosim_filename_length
implicit none
  character(len=*), parameter :: mod_filename = __FILE__

  integer :: arg_count
  character(len=ecosim_filename_length)      :: namelist_filename
  character(len=ecosim_namelist_buffer_size) :: namelist_buffer
  character(len=ecosim_filename_length)      :: base_filename

  arg_count = command_argument_count()
  if (arg_count /= 1) then
     write(stdout, '(a, i3)') 'ERROR: must pass exactly one arguement one the command line, received ', arg_count
     call usage()
     call endrun()
  end if

  call get_command_argument(1, namelist_filename)

  write(stdout, '(a, a)') 'Reading namelist filename : ', trim(namelist_filename)
  namelist_buffer = ''

  call namelist_to_buffer(namelist_filename, namelist_buffer)

!  call readNML(namelist_filename)

!  call RunModel(namelist_buffer)

end program main
