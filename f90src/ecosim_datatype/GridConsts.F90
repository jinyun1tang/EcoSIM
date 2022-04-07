module GridConsts
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  integer, PARAMETER :: JX=4
  integer, PARAMETER :: JY=4
  integer, PARAMETER :: JZ=20
  integer, PARAMETER :: JH=JX+1
  integer, PARAMETER :: JV=JY+1
  integer, PARAMETER :: JD=JZ+1
  integer, PARAMETER :: JP=5
  integer, PARAMETER :: JC=10
  integer, PARAMETER :: JS=5
  integer, PARAMETER :: JLI=4   !# of sectors for the leaf zenith [0,pi/2]
  integer, PARAMETER :: JLA=4   !# of sectors for the leaf azimuth, [0,pi]
  integer, PARAMETER :: JSA=4   !# of sectors for the sky azimuth  [0,2*pi]
  integer, parameter :: JNODS=25
  integer, parameter :: jcplx=5 !# of microbe-substrate complexes
  integer, parameter :: jcplx1=jcplx-1
  integer, parameter :: jsken=4 !# of kinetic components of the substrates
end module GridConsts
