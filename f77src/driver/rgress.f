      subroutine rgress(nmfile,case_name)

      use TestMod, only : regression              
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk3.h"

      character(len=*) :: nmfile
      character(len=*) :: case_name

      !local variables
      character(len=128) :: category
      character(len=128) :: name

      call regression%Init(trim(nmfile),case_name)

      if(regression%write_regression_output)then
        write(*,*)'write regression file'
        call regression%OpenOutput()


        category = 'Flux (g m^-3 h^-1)'
        name = 'NH4_UPTK'
        call regression%writedata(category,name,(/1.,1.,1.,1./))

        category = 'concentration (g m^3)'
        name = 'O2'
        call regression%writedata(category,name,(/1.,1.,1.,1./))

        category = 'Soil water (m^3 m^-3)'
        name = 'THETWZ'
        call regression%writedata(category,name,(/1.,1.,1.,1./))

        category = 'Soil Temperature (oC)'
        name = 'THETWZ'
        call regression%writedata(category,name,(/1.,1.,1.,1./))        

        call regression%CloseOutput()
      endif        
      end subroutine rgress    
