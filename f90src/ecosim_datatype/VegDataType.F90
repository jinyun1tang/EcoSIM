module VegDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public

  real(r8) :: ZSIN(4)
  real(r8) :: ZCOS(4)

  contains

  subroutine  InitVegData

  implicit none

  ZSIN=real((/0.195,0.556,0.831,0.981/),r8)
  ZCOS=real((/0.981,0.831,0.556,0.195/),r8)

  end subroutine InitVegData
end module VegDataType
