module MiscUtilMod

! DESCRIPTION:
!  subroutines/functions to do miscellaneous things
implicit none

  public :: addone
  contains

! ----------------------------------------------------------------------

  function addone(itemp)result(ans)
!
!  DESCRIPTION
! increase itemp by one
  implicit none
  integer, intent(inout) :: itemp

  integer :: ans

  itemp=itemp+1
  ans=itemp
  end function addone


end module MiscUtilMod
