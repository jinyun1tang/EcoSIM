module PhenologyDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8

  implicit none
  public
  save

  real(r8) :: FVRN(0:5)
  real(r8) :: FWOOD(0:1)
  real(r8) :: FWOODN(0:1)
  real(r8) :: FWOODP(0:1)
  REAL(R8) :: FWODB(0:1)
  real(r8) :: FWODLN(0:1)
  real(r8) :: FWODLP(0:1)
  REAL(R8) :: FWODSN(0:1)
  real(r8) :: FWODSP(0:1)
  real(r8) :: FWODR(0:1)
  real(r8) :: FWODRN(0:1)
  real(r8) :: FWODRP(0:1)
  
  contains
  subroutine InitPhenologyData

  implicit none


  FVRN =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)
  end subroutine InitPhenologyData
end module PhenologyDataType
