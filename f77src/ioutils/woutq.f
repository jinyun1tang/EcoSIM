
      SUBROUTINE woutq(I,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE WRITES OUT THE NAMES OF ALL PLANT SPECIES TO
C     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
C     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
C     TO INITIALIZE LATER MODEL RUNS
C
      use data_kind_mod, only : r8 => SHR_KIND_R8
      integer, intent(in) :: I,NHW,NHE,NVN,NVS

      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"

C     execution begins here

      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      WRITE(30,90)I,IDATA(3),NP(NY,NX)
     2,(DATAP(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NP(NY,NX))
90    FORMAT(2I4,1I3,5(A16,I4))
9990  CONTINUE
9995  CONTINUE
1000  RETURN
      END
