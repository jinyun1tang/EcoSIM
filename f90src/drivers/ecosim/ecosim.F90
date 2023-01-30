PROGRAM main
!!
! Description:
! THIS SUBROUTINE READS THE RUNSCRIPT AND ENTERS FILENAMES INTO DATA ARRAYS
! FOR USE IN 'READS' AND 'READQ'. WHEN FINISHED THIS SUBROUTINE CALLS
! 'SOIL' WHICH IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod     , only : r8 => SHR_KIND_R8
  use TestMod           , only : regression
  use InitEcoSIM        , only : InitModules
  use EcoSIMDesctruct   , only : DestructEcoSIM
  use EcoSIMAPI         , only : SetMesh
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcosimConst
  USE fileUtil          , ONLY : iulog
  use EcoSIMConfig      , only : case_name
  implicit none

  character(len=*), parameter :: mod_filename = __FILE__
  integer :: NA(250),ND(250)
  integer :: NAX,NDX,NEX,NAY,NDY,NE,N,NTX,NT
  integer :: NHW,NVN,NHE,NVS

  CHARACTER(len=640):: BUF
  character(len=80):: runfile
  character(len=36):: nmlfile
  logical :: is_dos
  integer :: nmicbguilds
!!
! begin_execution

  is_dos=.false.

  write(iulog,*)'obtain working directory'
  CALL GETCWD(BUF)
!
! IDENTIFY OPERATING SYSTEM: DOS OR UNIX

  IF((.NOT.(BUF(1:1).EQ.'/'.OR.BUF(1:1).EQ.'~')).AND.BUF(2:2).EQ.':')THEN
    write(iulog,*)'Running ecosim on dos system'
    is_dos=.true.
!   CALL GETARG(1,BUF)
!   OPEN(5,FILE=BUF,STATUS='OLD')
    CALL GETARG(2,BUF)
    CALL CHDIR(BUF)
    PREFIX='.\\'  ! location where input files are put
!   make output directory
  ELSE
    write(iulog,*)'Running ecosim on unix system'
!   make output directory
  ENDIF
!
  CALL GETARG(1,nmlfile)

  write(iulog,*)'read namelist'
  call readnamelist(trim(nmlfile),runfile, case_name, prefix, &
    do_rgres,LYRG,lverb, nmicbguilds)

  write(iulog,*)'read runfile',trim(runfile)
  OPEN(5,FILE=runfile,STATUS='OLD')
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
  write(iulog,*)'input files at:',trim(prefix)
  IGO=0
!
! NUMBER OF COLUMNS AND ROWS for the whole land scape
!
  READ(5,*)NHW,NVN,NHE,NVS

  call SetMesh(NHW,NVN,NHE,NVS)

  call  InitModules(nmicbguilds)

!
! SITE FILE
!
  READ(5,10)DATA1(1)
!
! TOPOGRAPHY FILE
!
  READ(5,10)DATA1(2)
!
! READ THE NUMBER OF TIMES THE SCENARIOS IN THE MODEL RUN ARE TO BE EXECUTED
!
100   READ(5,*)NAX,NDX
  IF(NAX.EQ.0.AND.NDX.EQ.0)GO TO 1000
!
!   NUMBER OF SCENES IN THE NEXT SCENARIO OF THE MODEL RUN
!   AND THE NUMBER OF TIMES THIS SCENARIO IS TO BE EXECUTED
!
  D105: DO NEX=1,NAX
    READ(5,*)NAY,NDY
    NA(NEX)=NAY
    ND(NEX)=NDY
!
!   FOR EACH SCENE IN THIS SCENARIO:
!
    D110: DO NE=1,NA(NEX)
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
      D115: DO N=21,30
        READ(5,10)DATAC(N,NE,NEX)
      ENDDO D115
    ENDDO D110
  ENDDO D105
!
!  RUN THIS SCRIPT
!
  D120: DO NTX=1,NDX      !repeat scenario NDX times, usually set to 1, though multiple scenes is possible
    DO  NEX=1,NAX         !run nax scenes for scenario ntx, each scene could have periods, where each period has multiple years
      DO  NT=1,ND(NEX)    !ND(NEX)=NDY, period NT
        DO  NE=1,NA(NEX)  !NA(NEX)=NAY, year NE in period NT
! run simulation for one year
! each year has its own climate forcing, land/pft management and io setup
          CALL SOIL(NA,ND,NT,NE,NAX,NTX,NEX,NHW,NHE,NVN,NVS)
          IGO=IGO+1
        enddo
      enddo
    enddo
  ENDDO D120
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
