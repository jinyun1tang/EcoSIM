module EcoSIMHistMod

!!
! data types of plant characteristics
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  public
  CHARACTER(len=16) :: DATA(30)
  CHARACTER(len=16) :: DATAC(30,250,250)
  CHARACTER(len=16) :: DATAP(JP,JY,JX)
  CHARACTER(len=16) :: DATAM(JP,JY,JX)
  CHARACTER(len=16) :: DATAX(JP),DATAY(JP)
  CHARACTER(len=16) :: DATAZ(JP,JY,JX)
  CHARACTER(len=16) :: OUTS(10)
  CHARACTER(len=16) :: OUTP(10)
  CHARACTER(len=16) :: OUTFILS(10,JY,JX)
  CHARACTER(len=16) :: OUTFILP(10,JP,JY,JX)
  CHARACTER(len=3) :: CHOICE(102,20)
  CHARACTER(len=8) :: CDATE
  CHARACTER(len=640):: PREFIX
  character(len=1280):: outdir
  integer :: IDATA(60),NOUTS(10),NOUTP(10)



end module EcoSIMHistMod
