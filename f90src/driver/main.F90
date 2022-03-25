PROGRAM main
!!
! Description:
! THIS SUBROUTINE READS THE RUNSCRIPT AND ENTERS FILENAMES INTO DATA ARRAYS
! FOR USE IN 'READS' AND 'READQ'. WHEN FINISHED THIS SUBROUTINE CALLS
! 'SOIL' WHICH IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod     , only : r8 => SHR_KIND_R8
  use TestMod           , only : regression
  use InitEcoSIM        , only :  InitModules
  use EcoSIMDesctruct   , only : DestructEcoSIM
  use EcoSIMCtrlDataType
  use GridConsts
  use EcoSIMHistMod
  implicit none

  character(len=*), parameter :: mod_filename = __FILE__
  integer :: NA(250),ND(250)
  integer :: NAX,NDX,NEX,NAY,NDY,NE,N,NTX,NT
  integer :: NHW,NVN,NHE,NVS

  CHARACTER(len=640):: BUF
  character(len=80):: runfile
  character(len=36):: case_name
  character(len=36):: nmlfile
  logical :: is_dos
  integer :: nmicbguilds
!!
! begin_execution

  is_dos=.false.

  write(*,*)'obtain working directory'
  CALL GETCWD(BUF)
!
! IDENTIFY OPERATING SYSTEM: DOS OR UNIX
!
  CALL GETARG(1,nmlfile)
  write(*,*)'read namelist'
  call readnamelist(trim(nmlfile),runfile, case_name, prefix, &
    do_rgres,LYRG,lverb, nmicbguilds)

  call  InitModules(nmicbguilds)

  write(*,*)'read runfile'
  OPEN(5,FILE=runfile,STATUS='OLD')
  IF((.NOT.(BUF(1:1).EQ.'/'.OR.BUF(1:1).EQ.'~')).AND.BUF(2:2).EQ.':')THEN
    write(*,*)'dos system'
    is_dos=.true.
!   CALL GETARG(1,BUF)
!   OPEN(5,FILE=BUF,STATUS='OLD')
    CALL GETARG(2,BUF)
    CALL CHDIR(BUF)
    PREFIX='.\\'  ! location where input files are put
!   make output directory
  ELSE
    write(*,*)'unix system'
!   make output directory
  ENDIF
!
! READ INPUT FILES
!
10    FORMAT(A16)

  if(is_dos)then
    outdir=trim(buf)//'\\'//trim(case_name)//'_outputs\\'
  else
    outdir=trim(buf)//'/'//trim(case_name)//'_outputs/'
  endif
  call system('mkdir -p '//trim(outdir))
  write(*,*)'input files at:',trim(prefix)
  IGO=0
!
! NUMBER OF COLUMNS AND ROWS
!
  READ(5,*)NHW,NVN,NHE,NVS
!
! SITE FILE
!
  READ(5,10)DATA(1)
!
! TOPOGRAPHY FILE
!
  READ(5,10)DATA(2)
!
! READ THE NUMBER OF TIMES THE SCENARIOS IN THE MODEL RUN ARE TO BE EXECUTED
!
100   READ(5,*)NAX,NDX
  IF(NAX.EQ.0.AND.NDX.EQ.0)GO TO 1000
!
!   NUMBER OF SCENES IN THE NEXT SCENARIO OF THE MODEL RUN
!   AND THE NUMBER OF TIMES THIS SCENARIO IS TO BE EXECUTED
!
  DO 105 NEX=1,NAX
    READ(5,*)NAY,NDY
    NA(NEX)=NAY
    ND(NEX)=NDY
!
!   FOR EACH SCENE IN THIS SCENARIO:
!
    DO 110 NE=1,NA(NEX)
!
!     WEATHER FILE
!
      READ(5,10,END=1000)DATAC(3,NE,NEX)
!
!     WEATHER OPTIONS
!
      READ(5,10)DATAC(4,NE,NEX)
!
!     LAND MANAGEMENT FILE
!
      READ(5,10)DATAC(9,NE,NEX)
!
!     PLANT MANAGEMENT FILE
!
      READ(5,10)DATAC(10,NE,NEX)
!
!     OUTPUT DATA CONTROL
!
      DO 115 N=21,30
        READ(5,10)DATAC(N,NE,NEX)
115   CONTINUE
110   CONTINUE
105   CONTINUE
!
!  RUN THIS SCRIPT
!
  DO 120 NTX=1,NDX
    DO  NEX=1,NAX
      DO  NT=1,ND(NEX)
        DO  NE=1,NA(NEX)
          CALL SOIL(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
          IGO=IGO+1
        enddo
      enddo
    enddo
120   CONTINUE
!
!   SCRIPT COMPLETED, START NEXT SCRIPT
!
  GO TO 100
  close(5)
1000  continue
  if(do_rgres)then
    call regressiontest(trim(nmlfile),trim(case_name),NHW,NVN)
  endif
  call DestructEcoSIM
END program main
