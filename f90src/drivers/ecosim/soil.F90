SUBROUTINE soil(NE,NEX,NHW,NHE,NVN,NVS)
!!
! Description:
! THIS IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use DayMod       , only : day
  use ExecMod      , only : exec
  use StarteMod    , only : starte
  use StartqMod    , only : startq
  use StartsMod    , only : starts
  use VisualMod    , only : visual
  use WthrMod      , only : wthr
  use RestartMod   , only : restart
  use PlantInfoMod , only : ReadPlantInfo
  use readsmod     , only : reads
  use Hist1Mod     , only : fouts,foutp,outpd,outph,outsd,outsh
  use timings      , only : init_timer, start_timer, end_timer,end_timer_loop
  use InitEcoSIM   , only : InitModules2
  use EcoSIMCtrlMod
  use PlantAPI     , only : PlantModel
  use MicBGCAPI    , only : MicrobeModel, MicAPI_Init, MicAPI_cleanup
  use ForcWriterMod, only : do_bgcforc_write,WriteBBGCForc
  use EcoSIMAPI    , only : Run_EcoSIM_one_step
  use GridConsts
  use EcoSIMCtrlMod, only : frectyp
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use HistDataType , only : hist_ecosim
  use HistFileMod  , only : hist_update_hbuf,hist_htapes_wrapup
  use GridConsts   , only : bounds

  implicit none
  integer :: yearc, yeari
  integer, intent(in) :: NE,NEX,NHW,NHE,NVN,NVS

  character(len=*), parameter :: mod_filename = __FILE__
  real(r8) :: t1
  integer :: I,J
  integer :: idaz
  character(len=14) :: ymdhs
  logical :: nlend, rstwr, lnyr
! begin_execution
!
! READ INPUT DATA FOR SITE, SOILS AND MANAGEMENT IN 'READS'
! AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTS'
!
333   FORMAT(A8)

  is_first_year=frectyp%yearacc.EQ.0

  call init_timer(outdir)

  if(lverb)WRITE(*,333)'READS: read climate forcing'
  CALL READS(frectyp%yearcur,frectyp%yearclm,NE,NEX,NHW,NHE,NVN,NVS)

  !temporary set up for setting mass balance check
  IBEGIN=1;ISTART=1;ILAST=0
  
!  if(lverb)WRITE(*,333)'FOUTS'
!  CALL FOUTS(NE,NEX,NHW,NHE,NVN,NVS)

  call etimer%get_ymdhs(ymdhs)
  
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4))THEN
!   the first simulation year
    if(lverb)WRITE(*,333)'STARTS'
    CALL STARTS(NHW,NHE,NVN,NVS)
!
!   RECOVER VALUES OF ALL SOIL STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTS' IF NEEDED
!
    IF(continue_run)THEN
!IRUN: start date of current scenario
      if(lverb)WRITE(*,333)'ROUTS'
      CALL ROUTS(NHW,NHE,NVN,NVS)
    ENDIF
  ENDIF
!
  if(lverb)WRITE(*,333)'ReadPlantInfo'
  call ReadPlantInfo(frectyp%yearcur,frectyp%yearclm,NE,NEX,NHW,NHE,NVN,NVS)

  if(lverb)WRITE(*,333)'FOUTP'
  CALL FOUTP(NE,NEX,NHW,NHE,NVN,NVS)

! INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
!
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4))THEN
    if(lverb)WRITE(*,333)'STARTQ'
    CALL STARTQ(NHW,NHE,NVN,NVS,1,JP)
!
!   RECOVER VALUES OF ALL PLANT STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTP' IF NEEDED
!
    IF(continue_run)THEN
      if(lverb)WRITE(*,333)'ROUTP'
      CALL ROUTP(NHW,NHE,NVN,NVS)
    ENDIF
  ENDIF

! INITIALIZE ALL SOIL CHEMISTRY VARIABLES IN 'STARTE'
!
  if(lverb)WRITE(*,333)'STARTE'
  CALL STARTE(NHW,NHE,NVN,NVS)

  iyear_cur=frectyp%yearcur
  LYRC=etimer%get_days_cur_year()
  
  DO I=1,LYRC
    IF(do_rgres .and. I.eq.LYRG)RETURN
    !   UPDATE DAILY VARIABLES SUCH AS MANAGEMENT INPUTS
    !
    if(lverb)WRITE(*,333)'DAY'
    CALL DAY(I,NHW,NHE,NVN,NVS)

    DO J=1,24
      call etimer%get_ymdhs(ymdhs)
      if(ymdhs==frectyp%ymdhs0)frectyp%lskip_loop=.false.       
      if(frectyp%lskip_loop)cycle      

    !
    !   UPDATE HOURLY VARIABLES IN 'HOUR1'
    !   set up climate forcing for the new hour
      if(lverb)WRITE(*,333)'WTHR'
      call start_timer(t1)
      CALL WTHR(I,J,NHW,NHE,NVN,NVS)
      call end_timer('WTHR',t1)

      if(lverb)WRITE(*,333)'Run_EcoSIM_one_step'
      call Run_EcoSIM_one_step(I,J,NHW,NHE,NVN,NVS)
  !
  !   WRITE HOURLY SOIL AND PLANT OUTPUT IN 'OUTSH' AND 'OUTPH'
  !
!      IF((J/JOUT)*JOUT.EQ.J)THEN
!        if(lverb)WRITE(*,333)'OUTSH'
!        CALL OUTSH(I,J,NE,NEX,NHW,NHE,NVN,NVS)

!        if(lverb)WRITE(*,333)'OUTPH'
!        CALL OUTPH(I,J,NE,NEX,NHW,NHE,NVN,NVS)
!      ENDIF

  !   WRITE OUTPUT FOR DYNAMIC VISUALIZATION
  !
      IF(visual_out)THEN
        IF((J/JOUT)*JOUT.EQ.J)THEN
          if(lverb)WRITE(*,333)'VIS'
          CALL VISUAL(I,J,NHW,NHE,NVN,NVS)
        ENDIF
      ENDIF
    
      call end_timer_loop()
      call hist_ecosim%hist_update(bounds)
      call hist_update_hbuf(bounds)      
      call etimer%update_time_stamp()      

      nlend=etimer%its_time_to_exit()
      rstwr=.false.
      lnyr=etimer%its_a_new_year().and.hist_yrclose
      call hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )      
    END DO
    
    IF(restart_out.AND.KOUT.GT.0)THEN
!
!   WRITE ALL SOIL AND PLANT STATE VARIABLES AND OTHER INFORMATION
!   NEEDED TO RE-INITIALIZE THE MODEL TO CHECKPOINT FILES
!   IN 'WOUTS', 'WOUTP' AND 'WOUTQ'
!
      IF((I/KOUT)*KOUT.EQ.I.OR.I.EQ.IFIN)THEN
        call restart(I,NHW,NHE,NVN,NVS)
      ENDIF
    ENDIF
!
!   WRITE DAILY SOIL AND PLANT OUTPUT IN 'OUTSD' AND 'OUTPD'
!
!    IF((I/IOUT)*IOUT.EQ.I)THEN
!      if(lverb)WRITE(*,333)'OUTSD'
!      CALL OUTSD(I,NE,NEX,NHW,NHE,NVN,NVS)
!      if(lverb)WRITE(*,333)'OUTPD'
!      CALL OUTPD(I,NE,NEX,NHW,NHE,NVN,NVS)
!    ENDIF
!
! PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
!
    if(lverb)WRITE(*,333)'EXEC'
    CALL EXEC(I)
!
! RE-INITIALIZE MODEL FROM CHECKPOINT FILES IF NEEDED
!
!    IF(NYR.NE.1.AND.IDAYR.NE.IOLD)THEN
!   reinitialize the model from checkfiles
!      if(lverb)WRITE(*,333)'ROUTS'
!      CALL ROUTS(NHW,NHE,NVN,NVS)
!      if(lverb)WRITE(*,333)'ROUTP'
!      CALL ROUTP(NHW,NHE,NVN,NVS)
!    ENDIF

    if(do_bgcforc_write)then
      call WriteBBGCForc(I,IYRR)
    endif

  END DO
!
! WRITE OUTPUT FILES FOR EACH GRID CELL IN 'SPLIT'
!
!ifdef _WIN_
!  CALL SPLIT(NE,NEX,NHW,NHE,NVN,NVS)
!else
! CALL SPLITC(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!endif
! WRITE(*,333)'SPLIT'

  RETURN
END subroutine soil
