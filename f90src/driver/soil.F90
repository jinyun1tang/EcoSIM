SUBROUTINE soil(NA,ND,NT,NE,NAX,NTX,NEX,NHW,NHE,NVN,NVS)
!!
! Description:
! THIS IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use DayMod       , only : day
  use ErosionMod   , only : erosion
  use ExecMod      , only : exec
  use Hour1Mod     , only : hour1
  use RedistMod    , only : redist
  use soluteMod    , only : solute
  use StarteMod    , only : starte
  use StartqMod    , only : startq
  use StartsMod    , only : starts
  use TrnsfrMod    , only : trnsfr
  use TrnsfrsMod   , only : trnsfrs
  use VisualMod    , only : visual
  use WatsubMod    , only : watsub
  use WthrMod      , only : wthr
  use readiMod     , only : readi
  use readqmod     , only : readq
  use readsmod     , only : reads
  use Hist1Mod     , only : fouts,foutp,outpd,outph,outsd,outsh
  use timings      , only : init_timer, start_timer, end_timer,end_timer_loop
  use GridConsts
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use PlantAPI     , only : PlantModel
  use MicBGCAPI    , only : MicrobeModel
  implicit none

  integer, intent(in) :: NT,NE,NAX,NTX,NEX,NHW,NHE,NVN,NVS
  integer, intent(in) :: NA(1:NEX),ND(1:NEX)

  character(len=*), parameter :: mod_filename = __FILE__

  integer :: I,J
  integer, SAVE :: NF,NX,NTZ,NTZX,NFX
  DATA NF,NX,NTZ,NTZX,NFX/0,0,0,0,0/
  real(r8) :: t1

! begin_execution
!
! READ INPUT DATA FOR SITE, SOILS AND MANAGEMENT IN 'READS'
! AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTS'
!
333   FORMAT(A8)

  is_first_year=IGO.EQ.0
  call init_timer(outdir)

  IF(is_first_year)THEN
    if(lverb)WRITE(*,333)'READI'
    CALL READI(NA,ND,NT,NE,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)
  ENDIF

  if(lverb)WRITE(*,333)'READS'
  CALL READS(NA,ND,NT,NE,NAX,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)

  if(lverb)WRITE(*,333)'FOUTS'
  CALL FOUTS(NT,NE,NTX,NEX,NF,NFX,NHW,NHE,NVN,NVS)
!
! INITIALIZE ALL SOIL VARIABLES IN 'STARTS'
!
  IF((is_restart_run.AND.is_first_year).OR.IDAYR.NE.IOLD)THEN
    if(lverb)WRITE(*,333)'STARTS'
    CALL STARTS(NHW,NHE,NVN,NVS)
!
!   RECOVER VALUES OF ALL SOIL STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTS' IF NEEDED
!
    IF(is_restart_run)THEN
      IF((IDAYR.GE.IRUN.AND.IYRR.EQ.IDATA(9)).OR.IYRR.GT.IDATA(9))THEN
        if(lverb)WRITE(*,333)'ROUTS'
        CALL ROUTS(NHW,NHE,NVN,NVS)
      ENDIF
    ENDIF
  ENDIF
!
! RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
!
  if(lverb)WRITE(*,333)'ROUTQ'
  CALL ROUTQ(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!
!   READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
!   AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
!
  if(lverb)WRITE(*,333)'READQ'
  CALL READQ(NA,ND,NT,NE,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)

  if(lverb)WRITE(*,333)'FOUTP'
  CALL FOUTP(NT,NE,NTX,NEX,NF,NFX,NHW,NHE,NVN,NVS)
!
! INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
!
  IF((is_restart_run.AND.is_first_year).OR.IDAYR.NE.IOLD)THEN
    if(lverb)WRITE(*,333)'STARTQ'
    CALL STARTQ(NHW,NHE,NVN,NVS,1,5)
!
!   RECOVER VALUES OF ALL PLANT STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTP' IF NEEDED
!
    IF(is_restart_run)THEN
      IF((IDAYR.GT.IRUN.AND.IYRR.EQ.IDATA(9)).OR.IYRR.GT.IDATA(9))THEN
        if(lverb)WRITE(*,333)'ROUTP'
        CALL ROUTP(NHW,NHE,NVN,NVS)
      ENDIF
    ENDIF
  ENDIF
!
! INITIALIZE ALL SOIL CHEMISTRY VARIABLES IN 'STARTE'
!
  if(lverb)WRITE(*,333)'STARTE'
  CALL STARTE(NHW,NHE,NVN,NVS)
!
!   BEGIN DAILY TIME STEP
!
9000  I=IDAYR+1
  IF(I.NE.LYRC.OR.NYR.NE.0)THEN
    IF(I.GT.LYRX.AND.NYR.NE.0)THEN
      I=I-LYRX
      IYRR=IYRR+1
    ENDIF
  ENDIF
  IF(do_rgres .and. I.eq.LYRG)RETURN
  if(lverb)write(*,'(2(A,I4))')'I=',I,' IFIN=',IFIN
  IF(I.GT.IFIN)GO TO 9999
  if(lverb)WRITE(*,333)'FINISH'
!
!   UPDATE DAILY VARIABLES SUCH AS MANAGEMENT INPUTS
!
  if(lverb)WRITE(*,333)'DAY'
  CALL DAY(I,NHW,NHE,NVN,NVS)
  DO 9995 J=1,24
!
!   UPDATE HOURLY VARIABLES IN 'HOUR1'
!
    if(lverb)WRITE(*,333)'WTHR'
!    if(I>=170)print*,I,J,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL WTHR(I,J,NHW,NHE,NVN,NVS)
    call end_timer('WTHR',t1)

    if(lverb)WRITE(*,333)'HOUR1'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL HOUR1(I,J,NHW,NHE,NVN,NVS)
    call end_timer('HOUR1',t1)
!
!   CALCULATE SOIL ENERGY BALANCE, WATER AND HEAT FLUXES IN 'WATSUB'
!
    if(lverb)WRITE(*,333)'WAT'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL WATSUB(I,J,NHW,NHE,NVN,NVS)
    call end_timer('WAT',t1)
!
!   CALCULATE SOIL BIOLOGICAL TRANSFORMATIONS IN 'NITRO'
!
    if(lverb)WRITE(*,333)'NIT'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL MicrobeModel(I,J,NHW,NHE,NVN,NVS)
    call end_timer('NIT',t1)
!
!   UPDATE PLANT biogeochemistry
!
    if(lverb)WRITE(*,333)'PlantModel'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call PlantModel(I,J,NHW,NHE,NVN,NVS)
!
!
!   CALCULATE SOLUTE EQUILIBRIA IN 'SOLUTE'
!
    if(lverb)WRITE(*,333)'SOL'
    call start_timer(t1)
    CALL SOLUTE(I,J,NHW,NHE,NVN,NVS)
    call end_timer('SOL',t1)
!
!   CALCULATE GAS AND SOLUTE FLUXES IN 'TRNSFR'
!
    if(lverb)WRITE(*,333)'TRN'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL TRNSFR(I,J,NHW,NHE,NVN,NVS)
    call end_timer('TRN',t1)
!
!   CALCULATE ADDITIONAL SOLUTE FLUXES IN 'TRNSFRS' IF SALT OPTION SELECTED
!
    if(lverb)WRITE(*,333)'TRNS'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL TRNSFRS(I,J,NHW,NHE,NVN,NVS)
    call end_timer('TRNSFRS',t1)
!
!   CALCULATE SOIL SEDIMENT TRANSPORT IN 'EROSION'
!
    if(lverb)WRITE(*,333)'EROSION'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL EROSION(I,J,NHW,NHE,NVN,NVS)
    call end_timer('EROSION',t1)
!
!   UPDATE ALL SOIL STATE VARIABLES FOR WATER, HEAT, GAS, SOLUTE
!   AND SEDIMENT FLUXES IN 'REDIST'
!
    if(lverb)WRITE(*,333)'RED'
!    if(I>=170)print*,TKS(0,NVN,NHW)
    call start_timer(t1)
    CALL REDIST(I,J,NHW,NHE,NVN,NVS)
    call end_timer('RED',t1)
!    if(I>=170)print*,'afred',TKS(0,NVN,NHW)
!   WRITE(*,333)'END'
!
!   WRITE HOURLY SOIL AND PLANT OUTPUT IN 'OUTSH' AND 'OUTPH'
!
    IF((J/JOUT)*JOUT.EQ.J)THEN
      if(lverb)WRITE(*,333)'OUTSH'
      CALL OUTSH(I,J,NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)

      if(lverb)WRITE(*,333)'OUTPH'
      CALL OUTPH(I,J,NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
    ENDIF
!
!   WRITE OUTPUT FOR DYNAMIC VISUALIZATION
!
    IF(DATA1(18).EQ.'YES')THEN
      IF((J/JOUT)*JOUT.EQ.J)THEN
        if(lverb)WRITE(*,333)'VIS'
        CALL VISUAL(I,J,NHW,NHE,NVN,NVS)
      ENDIF
    ENDIF
    call end_timer_loop()

9995  CONTINUE
  IF(DATA1(19).EQ.'YES'.AND.KOUT.GT.0)THEN
!
!   WRITE ALL SOIL AND PLANT STATE VARIABLES AND OTHER INFORMATION
!   NEEDED TO RE-INITIALIZE THE MODEL TO CHECKPOINT FILES
!   IN 'WOUTS', 'WOUTP' AND 'WOUTQ'
!
    IF((I/KOUT)*KOUT.EQ.I.OR.I.EQ.IFIN)THEN
      if(lverb)WRITE(*,333)'WOUTS'
      CALL WOUTS(I,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'WOUTP'
      CALL WOUTP(I,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'WOUTQ'
      CALL WOUTQ(I,NHW,NHE,NVN,NVS)
    ENDIF
  ENDIF
!
!   WRITE DAILY SOIL AND PLANT OUTPUT IN 'OUTSD' AND 'OUTPD'
!
  IF((I/IOUT)*IOUT.EQ.I)THEN
    if(lverb)WRITE(*,333)'OUTSD'
    CALL OUTSD(I,NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
    if(lverb)WRITE(*,333)'OUTPD'
    CALL OUTPD(I,NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
  ENDIF
!
! PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
!
  if(lverb)WRITE(*,333)'EXEC'
  CALL EXEC(I)
!
! RE-INITIALIZE MODEL FROM CHECKPOINT FILES IF NEEDED
!
  IF(NYR.NE.1.AND.IDAYR.NE.IOLD)THEN
!   reinitialize the model from checkfiles
    if(lverb)WRITE(*,333)'ROUTS'
    CALL ROUTS(NHW,NHE,NVN,NVS)
    if(lverb)WRITE(*,333)'ROUTP'
    CALL ROUTP(NHW,NHE,NVN,NVS)
  ENDIF
  GO TO 9000
9999  CONTINUE
! WRITE(*,333)'LOOP'
!
! WRITE OUTPUT FILES FOR EACH GRID CELL IN 'SPLIT'
!
!ifdef _WIN_
  CALL SPLIT(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!else
! CALL SPLITC(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!endif
! WRITE(*,333)'SPLIT'
  RETURN
END subroutine soil
