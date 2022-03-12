
      SUBROUTINE woutq(I,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES OUT THE NAMES OF ALL PLANT SPECIES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!
      use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  use FlagDataType
      implicit none
      integer, intent(in) :: I,NHW,NHE,NVN,NVS


      include "filec.h"
      include "files.h"
      include "blkc.h"
      integer :: NX,NY,NZ
!     execution begins here

      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      WRITE(30,90)I,IDATA(3),NP(NY,NX) &
      ,(DATAP(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NP(NY,NX))
90    FORMAT(2I4,1I3,5(A16,I4))
9990  CONTINUE
9995  CONTINUE
!1000  RETURN
      END
