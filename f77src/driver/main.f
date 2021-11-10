      PROGRAM main
C     THIS SUBROUTINE READS THE RUNSCRIPT AND ENTERS FILENAMES INTO DATA ARRAYS
C     FOR USE IN 'READS' AND 'READQ'. WHEN FINISHED THIS SUBROUTINE CALLS
C     'SOIL' WHICH IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
C
      use data_kind_mod, only : r8 => SHR_KIND_R8
      use TestMod, only : regression

      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"

      DIMENSION NA(250),ND(250)
C      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
C     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
C     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
C      CHARACTER*3 CHOICE(102,20)
C      CHARACTER*8 CDATE
C      CHARACTER*80 PREFIX
      CHARACTER(len=80):: BUF
      character(len=80):: runfile
      character(len=36):: case_name
      character(len=36):: nmlfile
      logical :: is_dos

      is_dos=.false.

      print*,'obtain working directory'
      CALL GETCWD(BUF)
C
C     IDENTIFY OPERATING SYSTEM: DOS OR UNIX
C
      CALL GETARG(1,nmlfile)
      print*,'read namelist'
      call readnamelist(trim(nmlfile),runfile, case_name, prefix,
     2do_rgres,LYRG,lverb)
      print*,'read runfile'
      OPEN(5,FILE=runfile,STATUS='OLD')
      IF((.NOT.(BUF(1:1).EQ.'/'.OR.BUF(1:1).EQ.'~'))
     2.AND.BUF(2:2).EQ.':')THEN
      print*,'dos system'
      is_dos=.true.
C      CALL GETARG(1,BUF)
C      OPEN(5,FILE=BUF,STATUS='OLD')
      CALL GETARG(2,BUF)
      CALL CHDIR(BUF)
      PREFIX='.\\'  ! location where input files are put
C     make output directory
      ELSE
      print*,'unix system'
C     make output directory
      ENDIF
C
C     READ INPUT FILES
C
10    FORMAT(A16)

      if(is_dos)then
      outdir=trim(buf)//'\\'//trim(case_name)//'_outputs\\'
      else
      outdir=trim(buf)//'/'//trim(case_name)//'_outputs/'
      endif
      call system('mkdir -p '//trim(outdir))
      print*,'input files at:',trim(prefix)
      IGO=0

C
C     NUMBER OF COLUMNS AND ROWS
C
      READ(5,*)NHW,NVN,NHE,NVS
C
C     SITE FILE
C
      READ(5,10)DATA(1)
C
C     TOPOGRAPHY FILE
C
      READ(5,10)DATA(2)
C
C     READ THE NUMBER OF TIMES THE SCENARIOS IN THE MODEL RUN
C     ARE TO BE EXECUTED
C
100   READ(5,*)NAX,NDX
      IF(NAX.EQ.0.AND.NDX.EQ.0)GO TO 1000
C
C     NUMBER OF SCENES IN THE NEXT SCENARIO OF THE MODEL RUN
C     AND THE NUMBER OF TIMES THIS SCENARIO IS TO BE EXECUTED
C
      DO 105 NEX=1,NAX
      READ(5,*)NAY,NDY
      NA(NEX)=NAY
      ND(NEX)=NDY
C
C     FOR EACH SCENE IN THIS SCENARIO:
C
      DO 110 NE=1,NA(NEX)
C
C     WEATHER FILE
C
      READ(5,10,END=1000)DATAC(3,NE,NEX)
C
C     WEATHER OPTIONS
C
      READ(5,10)DATAC(4,NE,NEX)
C
C     LAND MANAGEMENT FILE
C
      READ(5,10)DATAC(9,NE,NEX)
C
C     PLANT MANAGEMENT FILE
C
      READ(5,10)DATAC(10,NE,NEX)
C
C     OUTPUT DATA CONTROL
C
      DO 115 N=21,30
      READ(5,10)DATAC(N,NE,NEX)
115   CONTINUE
110   CONTINUE
105   CONTINUE
C
C     RUN THIS SCRIPT
C
      DO 120 NTX=1,NDX
      DO 120 NEX=1,NAX
      DO 120 NT=1,ND(NEX)
      DO 120 NE=1,NA(NEX)
      CALL SOIL(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      IGO=IGO+1
120   CONTINUE
C
C     SCRIPT COMPLETED, START NEXT SCRIPT
C
      GO TO 100
      close(5)
1000  if(do_rgres)then
        call regressiontest(trim(nmlfile),trim(case_name),NHW,NVN)
      endif
      STOP
      END
