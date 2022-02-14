subroutine regressiontest(nmfile,case_name, NX, NY)
!!
! Description:
! subroutine to conduct regression tests

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use TestMod, only : regression
  use minimathmod, only : safe_adb
  use SoilChemDataType
  implicit none

  include "parameters.h"
  include "files.h"
  include "blkc.h"
  include "blk1cp.h"
  include "blk1cr.h"
  include "blk1g.h"
  include "blk1n.h"
  include "blk1p.h"
  include "blk3.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk9a.h"
  include "blk9b.h"
  include "blk9c.h"
  include "blk11a.h"
  include "blk12a.h"
  include "blk12b.h"
  include "blk13b.h"
  include "blk14.h"

  character(len=*), parameter :: mod_filename = __FILE__
  character(len=*), intent(in) :: nmfile
  character(len=*), intent(in) :: case_name
  integer, intent(in) :: NX
  integer, intent(in) :: NY
  !local variables
  character(len=128) :: category
  character(len=128) :: name
  integer :: nz, ll
  real(r8) :: datv(12)

  call regression%Init(trim(nmfile),case_name)

  if(regression%write_regression_output)then
    write(*,*)'write regression file'
    call regression%OpenOutput()

    do NZ=1,NP(NY,NX)
      IF(IFLGC(NZ,NY,NX).EQ.1)THEN

        category = 'flux'
        name = 'NH4_UPTK (g m^-3 h^-1)'
        datv=0._r8
        do ll=1,12
          if(AREA(3,ll,NY,NX)>0._r8)then
            datv(ll)=safe_adb(RUPNH4(1,ll,NZ,NY,NX)+RUPNH4(2,ll,NZ,NY,NX) &
              +RUPNHB(1,ll,NZ,NY,NX)+RUPNHB(2,ll,NZ,NY,NX),AREA(3,ll,NY,NX))
          endif
        enddo
        call regression%writedata(category,name, datv)
        exit
      endif

    enddo

    category = 'state'
    name = 'aqueous soil O2 (g m^3)'
    datv=COXYS(1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'liquid soil water (m^3 m^-3)'
    datv=THETWZ(1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'soil temperature (oC)'
    datv=TCS(1:12,NY,NX)
    call regression%writedata(category,name,datv)
    write(*,*)'Finish regression file writing'
    call regression%CloseOutput()
  endif
end subroutine regressiontest
