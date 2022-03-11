module GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  public
  save

  integer, PARAMETER :: JX=4
  integer, PARAMETER :: JY=4
  integer, PARAMETER :: JZ=20
  integer, PARAMETER :: JH=JX+1
  integer, PARAMETER :: JV=JY+1
  integer, PARAMETER :: JD=JZ+1
  integer, PARAMETER :: JP=5
  integer, PARAMETER :: JC=10
  integer, PARAMETER :: JS=5
  integer, PARAMETER :: JLI=4
  integer, PARAMETER :: JLA=4
  integer, PARAMETER :: JSA=4

  integer :: NU(JY,JX)
  integer :: NUI(JY,JX)
  integer :: NJ(JY,JX)
  integer :: NK(JY,JX)
  integer :: NLI(JV,JH)
  integer :: NL(JV,JH)
  integer :: NUM(JY,JX)


end module GridDataType
