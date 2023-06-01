subroutine soil(NE,NEX,NHW,NHE,NVN,NVS,nlend)
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
  use RestartMod   , only : restFile
  use PlantInfoMod , only : ReadPlantInfo
  use readsmod     , only : reads
  use timings      , only : init_timer, start_timer, end_timer,end_timer_loop
  use InitEcoSIM   , only : InitModules2
  use EcoSIMCtrlMod
  use PlantAPI     , only : PlantModel
  use MicBGCAPI    , only : MicrobeModel, MicAPI_Init, MicAPI_cleanup
  use ForcWriterMod, only : do_bgcforc_write,WriteBBGCForc
  use EcoSIMAPI    , only : Run_EcoSIM_one_step
  use GridConsts
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use HistDataType , only : hist_ecosim
  use HistFileMod  , only : hist_update_hbuf,hist_htapes_wrapup

  implicit none
  integer :: yearc, yeari
  integer, intent(in) :: NE,NEX,NHW,NHE,NVN,NVS
  logical, intent(out) :: nlend
  character(len=*), parameter :: mod_filename = __FILE__
  real(r8) :: t1
  integer :: I,J
  integer :: idaz
  character(len=14) :: ymdhs
  logical :: rstwr, lnyr
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
  
  call etimer%get_ymdhs(ymdhs)
  
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4))THEN
!   the first simulation year
    if(lverb)WRITE(*,333)'STARTS'
    CALL STARTS(NHW,NHE,NVN,NVS)
!
!   RECOVER VALUES OF ALL SOIL STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTS' IF NEEDED

  ENDIF
!
  if(plant_model)then
    !plant information is read in every year, but the active flags
    !are set using the checkpoint file.
    if(lverb)WRITE(*,333)'ReadPlantInfo'  
    WRITE(*,333)'ReadPlantInfo'  
    call ReadPlantInfo(frectyp%yearcur,frectyp%yearclm,NE,NEX,NHW,NHE,NVN,NVS)
  endif

! INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
!
  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4) .and. plant_model)THEN
    !initialize by year
    if(lverb)WRITE(*,333)'STARTQ'
    CALL STARTQ(NHW,NHE,NVN,NVS,1,JP)

  ENDIF

  if(soichem_model)then
! INITIALIZE ALL SOIL CHEMISTRY VARIABLES IN 'STARTE' 
! This is done done every year, because tracer concentrations
! in rainfall vary every year. In a more reasonable way, e.g., 
! when coupled to atmospheric chemistry code, it should be done by 
! hour
    if(lverb)WRITE(*,333)'STARTE'
    CALL STARTE(NHW,NHE,NVN,NVS)
  endif

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
      
      if(ymdhs==frectyp%ymdhs0)then
        frectyp%lskip_loop=.false. 
        if(is_restart())then          
          call restFile(flag='read')
        endif
      endif        
      if(frectyp%lskip_loop)then
        call etimer%update_time_stamp()      
        cycle
      endif

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
      rstwr=etimer%its_time_to_write_restart()
      lnyr=etimer%its_a_new_year().and.hist_yrclose
      
      call hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )      
      if(rstwr)then
        call restFile(flag='write')
      endif
      if(nlend)exit      
    END DO
    
!
! PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
!
    if(lverb)WRITE(*,333)'EXEC'
    CALL EXEC(I)
!

    if(do_bgcforc_write)then
      call WriteBBGCForc(I,IYRR)
    endif
    if(nlend)exit
  END DO

  RETURN
END subroutine soil
