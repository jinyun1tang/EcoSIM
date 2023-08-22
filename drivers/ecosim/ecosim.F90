PROGRAM main
!!
! Description:
! THIS SUBROUTINE READS THE RUNSCRIPT AND ENTERS FILENAMES INTO DATA ARRAYS
! FOR USE IN 'READS' AND 'READQ'. WHEN FINISHED THIS SUBROUTINE CALLS
! 'SOIL' WHICH IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use TestMod           , only : regression
  use InitEcoSIM        , only : InitModules
  use EcoSIMDesctruct   , only : DestructEcoSIM
  use GridMod           , only : SetMesh
  use readiMod          , only : readi
  USE fileUtil          , ONLY : iulog,ecosim_namelist_buffer_size,namelist_to_buffer
  use HistFileMod       , only : hist_htapes_build
  use EcoSIMConfig      , only : case_name,set_sim_type,start_date,is_restart,datestrlen
  use StartsMod         , only : set_ecosim_solver
  use RestartMod        , only : get_restart_date
  use MicBGCAPI         , only : MicAPI_Init, MicAPI_cleanup
  use EcoSIMCtrlMod  
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMAPI         , only : soil,readnamelist,regressiontest     
  use EcosimConst
  implicit none

  character(len=*), parameter :: mod_filename = __FILE__

  integer :: NA(250),ND(250)
  integer :: NAX,NDX,NEX,NAY,NDY,NE,N,NTX,NT
  integer :: NHW,NVN,NHE,NVS
  integer :: nn1,nn2,nn3
  integer :: year_ini,nyr1,yeari,nstopyr
  CHARACTER(len=640):: BUF
  character(len=36):: nmlfile
  character(len=14) :: ymdhs
  character(len=14) :: ymdhs0 
  logical :: is_dos,nlend
  character(len=datestrlen) :: curr_date
  integer :: nmicbguilds
  character(len=ecosim_namelist_buffer_size) :: nml_buffer

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

  call namelist_to_buffer(nmlfile,nml_buffer)

  write(iulog,*)'read namelist'
  call readnamelist(nml_buffer, case_name, prefix, LYRG, nmicbguilds)

  if(is_dos)then
    outdir=trim(buf)//'\\'//trim(case_name)//'_outputs\\'
  else
    outdir=trim(buf)//'/'//trim(case_name)//'_outputs/'
  endif
  call system('mkdir -p '//trim(outdir))

  read(start_date,'(I4)')year_ini  

! NUMBER OF COLUMNS AND ROWS for the whole land scape
!
  call SetMesh(NHW,NVN,NHE,NVS)

  call  InitModules(nmicbguilds)

  if(lverb)WRITE(*,*)'read initialization information READI'
  CALL readi(NE,NEX,NHW,NHE,NVN,NVS)

  NE=1;NEX=1

  call set_sim_type()

  nstopyr=get_sim_len(forc_periods)

  call frectyp%Init()

  if(is_restart())then
    print*,'read restart/checkpoint info file: ecosim_rst'
    call get_restart_date(curr_date)
    frectyp%ymdhs0=curr_date    !the actual simulation beginning year
  else  
    frectyp%ymdhs0=start_date   !the actual simulation beginning year
  endif
  
  call hist_htapes_build()

  IGO=0
  yeari=year_ini
  print*,frectyp%ymdhs0,yeari
  
  DO nn1=1,3
    call set_ecosim_solver(NPXS(NN1),NPYS(NN1))

   !set up output frequency
    JOUT=JOUTS(NN1)   !frequency on hourly scale

    call MicAPI_Init
    
    do nn2=1,forc_periods(nn1*3)
      nn3=(nn1-1)*3
      do nyr1=forc_periods(nn3+1),forc_periods(nn3+2)
        
        frectyp%yearclm=nyr1
        frectyp%yearcur=etimer%get_curr_yearAD()
        nlend=.false.
        if(frectyp%yearcur==yeari)then
          call soil(NE,NEX,NHW,NHE,NVN,NVS,nlend)
        endif
        if(nlend)exit
        frectyp%yearacc=frectyp%yearacc+1
        call etimer%get_ymdhs(ymdhs)        
        if(.not.frectyp%lskip_loop)print*,frectyp%yearcur,nyr1,ymdhs
        yeari=yeari+1
      end do  
      if(nlend)exit
    end do
    call MicAPI_cleanup
    if(nlend)exit
  end do  

  if(do_rgres)then
    call regressiontest(trim(nmlfile),trim(case_name),NHW,NVN)
  endif
  call DestructEcoSIM
END program main

