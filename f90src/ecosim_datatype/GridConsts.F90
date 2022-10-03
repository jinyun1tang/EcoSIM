module GridConsts
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  integer :: JX=6
  integer :: JY=4
  integer :: JZ=20
  integer :: JH
  integer :: JV
  integer :: JD

  integer, PARAMETER :: JP=5
  integer, PARAMETER :: JC=10   !# of canopy layers
  integer, PARAMETER :: JS=5
  integer, PARAMETER :: JLI=4   !# of sectors for the leaf zenith [0,pi/2]
  integer, PARAMETER :: JLA=4   !# of sectors for the leaf azimuth, [0,pi]
  integer, PARAMETER :: JSA=4   !# of sectors for the sky azimuth  [0,2*pi]
  integer, parameter :: JNODS=25!# of nodes for plant canopy
  integer  :: JG
end module GridConsts
