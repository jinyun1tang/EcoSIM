      SUBROUTINE soil(NA,ND,NT,NE,NAX,NDX,NTX,NEX
     2,NHW,NHE,NVN,NVS)
C
C     THIS IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
C
      use data_kind_mod, only : r8 => SHR_KIND_R8
      use DayMod       , only : day
      use ErosionMod   , only : erosion
      use ExecMod      , only : exec
      use ExtractMod   , only : extract
      use grosubMod    , only : grosub
      use HfuncMod     , only : hfunc
      use Hour1Mod     , only : hour1
      use nitroMod     , only : nitro
      use RedistMod    , only : redist
      use soluteMod    , only : solute
      use StarteMod    , only : starte
      use StartqMod    , only : startq
      use StartsMod    , only : starts
      use StomateMod   , only : stomate
      use TrnsfrMod    , only : trnsfr
      use TrnsfrsMod   , only : trnsfrs
      use UptakeMod    , only : uptake
      use VisualMod    , only : visual
      use WatsubMod    , only : watsub
      use WthrMod      , only : wthr
      use :: timings

      implicit none

      integer, intent(in) :: NA(1:NEX),ND(1:NEX)
      integer, intent(in) :: NT,NE,NAX,NDX,NTX,NEX
     2,NHW,NHE,NVN,NVS

      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"

      integer :: I,J
      integer, SAVE :: NF,NX,NTZ,NTZX, NFX
      DATA NF,NX,NTZ,NTZX/0,0,0,0/
      real*8 :: t1, t2

C     execution begins here

C
C     READ INPUT DATA FOR SITE, SOILS AND MANAGEMENT IN 'READS'
C     AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTS'
C
      call init_timer(outdir)

      IF(IGO.EQ.0)THEN
      if(lverb)WRITE(*,333)'READI'
      CALL READI(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
      ENDIF
      if(lverb)WRITE(*,333)'READS'
      CALL READS(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'FOUTS'
      CALL FOUTS(NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NHW,NHE,NVN,NVS)
C
C     INITIALIZE ALL SOIL VARIABLES IN 'STARTS'
C
      IF((DATA(20).EQ.'YES'.AND.IGO.EQ.0).OR.IDAYR.NE.IOLD)THEN
      if(lverb)WRITE(*,333)'STARTS'
      CALL STARTS(NHW,NHE,NVN,NVS)
C
C     RECOVER VALUES OF ALL SOIL STATE VARIABLES FROM EARLIER RUN
C     IN 'ROUTS' IF NEEDED
C
      IF(DATA(20).EQ.'YES')THEN
      IF((IDAYR.GE.IRUN.AND.IYRR.EQ.IDATA(9))
     2.OR.IYRR.GT.IDATA(9))THEN
      if(lverb)WRITE(*,333)'ROUTS'
      CALL ROUTS(NHW,NHE,NVN,NVS)
      ENDIF
      ENDIF
      ENDIF
C
C     RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
C
      if(lverb)WRITE(*,333)'ROUTQ'
      CALL ROUTQ(NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
C
C     READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
C     AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
C
      if(lverb)WRITE(*,333)'READQ'
      CALL READQ(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'FOUTP'
      CALL FOUTP(NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NHW,NHE,NVN,NVS)
C
C     INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
C
      IF((DATA(20).EQ.'YES'.AND.IGO.EQ.0).OR.IDAYR.NE.IOLD)THEN
      if(lverb)WRITE(*,333)'STARTQ'
      CALL STARTQ(NHW,NHE,NVN,NVS,1,5)
C
C     RECOVER VALUES OF ALL PLANT STATE VARIABLES FROM EARLIER RUN
C     IN 'ROUTP' IF NEEDED
C
      IF(DATA(20).EQ.'YES')THEN
      IF((IDAYR.GT.IRUN.AND.IYRR.EQ.IDATA(9))
     2.OR.IYRR.GT.IDATA(9))THEN
      if(lverb)WRITE(*,333)'ROUTP'
      CALL ROUTP(NHW,NHE,NVN,NVS)
      ENDIF
      ENDIF
      ENDIF
C
C     INITIALIZE ALL SOIL CHEMISTRY VARIABLES IN 'STARTE'
C
      if(lverb)WRITE(*,333)'STARTE'
      CALL STARTE(NHW,NHE,NVN,NVS)
C
C     BEGIN DAILY TIME STEP
C
9000  I=IDAYR+1
      IF(I.NE.LYRC.OR.NYR.NE.0)THEN
      IF(I.GT.LYRX.AND.NYR.NE.0)THEN
      I=I-LYRX
      IYRR=IYRR+1
      ENDIF
      ENDIF
      IF(do_rgres .and. I.eq.LYRG)RETURN
      if(lverb)write(*,'(2(A,X,I4))')'I=',I,' IFIN=',IFIN
      IF(I.GT.IFIN)GO TO 9999
      if(lverb)WRITE(*,333)'FINISH'
C
C     UPDATE DAILY VARIABLES SUCH AS MANAGEMENT INPUTS
C
      if(lverb)WRITE(*,333)'DAY'
      CALL DAY(I,NHW,NHE,NVN,NVS)
      DO 9995 J=1,24
C
C     UPDATE HOURLY VARIABLES IN 'HOUR1'
C
      if(lverb)WRITE(*,333)'WTHR'
      call start_timer(t1)
      CALL WTHR(I,J,NHW,NHE,NVN,NVS)
      call end_timer('WTHR',t1)
      if(lverb)WRITE(*,333)'HOUR1'
333   FORMAT(A8)
      call start_timer(t1)
      CALL HOUR1(I,J,NHW,NHE,NVN,NVS)
      call end_timer('HOUR1',t1)
C
C     CALCULATE SOIL ENERGY BALANCE, WATER AND HEAT FLUXES IN 'WATSUB'
C
      if(lverb)WRITE(*,333)'WAT'
      call start_timer(t1)
      CALL WATSUB(I,J,NHW,NHE,NVN,NVS)
      call end_timer('WAT',t1)
C
C     CALCULATE SOIL BIOLOGICAL TRANSFORMATIONS IN 'NITRO'
C
      if(lverb)WRITE(*,333)'NIT'
      call start_timer(t1)
      CALL NITRO(I,J,NHW,NHE,NVN,NVS)
      call end_timer('NIT',t1)
C
C     UPDATE PLANT PHENOLOGY IN 'HFUNC'
C
      if(lverb)WRITE(*,333)'HFUNC'
      call start_timer(t1)
      CALL HFUNC(I,J,NHW,NHE,NVN,NVS)
      call end_timer('HFUNC',t1)
C
C     CALCULATE CANOPY CO2 UPTAKE AT FULL TURGOR, CANOPY WATER POTENTIAL,
C     HYDRAULIC AND STOMATAL RESISTANCES,AND CANOPY ENERGY BALANCE IN 'UPTAKE'
C     CALCULATE ROOT UPTAKE OF WATER, OXYGEN, NH4, NO3 AND PO4 IN 'UPTAKE'
C
      if(lverb)WRITE(*,333)'UPTK'
      call start_timer(t1)
      CALL UPTAKE(I,J,NHW,NHE,NVN,NVS)
      call end_timer('UPTK',t1)
C
C     CALCULATE CANOPY CO2 UPTAKE AT AMBIENT TURGOR, AUTOTROPHIC AND GROWTH
C     RESPIRATION, PLANT C ALLOCATION, CANOPY AND ROOT GROWTH IN 'GROSUB'
C
      if(lverb)WRITE(*,333)'GRO'
      call start_timer(t1)
      CALL GROSUB(I,J,NHW,NHE,NVN,NVS)
      call end_timer('GRO',t1)
C
C     CALCULATE ROOT-SOIL C AND NUTRIENT EXCHANGE FOR ALL PLANT SPECIES
C     IN 'EXTRACT'
C
      if(lverb)WRITE(*,333)'EXTR'
      call start_timer(t1)
      CALL EXTRACT(I,J,NHW,NHE,NVN,NVS)
      call end_timer('EXTR',t1)
C
C     CALCULATE SOLUTE EQUILIBRIA IN 'SOLUTE'
C
      if(lverb)WRITE(*,333)'SOL'
      call start_timer(t1)
      CALL SOLUTE(I,J,NHW,NHE,NVN,NVS)
      call end_timer('SOL',t1)
C
C     CALCULATE GAS AND SOLUTE FLUXES IN 'TRNSFR'
C
      if(lverb)WRITE(*,333)'TRN'
      call start_timer(t1)
      CALL TRNSFR(I,J,NHW,NHE,NVN,NVS)
      call end_timer('TRN',t1)
C
C     CALCULATE ADDITIONAL SOLUTE FLUXES IN 'TRNSFRS' IF SALT OPTION SELECTED
C
      if(lverb)WRITE(*,333)'TRNS'
      call start_timer(t1)
      CALL TRNSFRS(I,J,NHW,NHE,NVN,NVS)
      call end_timer('TRNSFRS',t1)
C
C     CALCULATE SOIL SEDIMENT TRANSPORT IN 'EROSION'
C
      if(lverb)WRITE(*,333)'EROSION'
      call start_timer(t1)
      CALL EROSION(I,J,NHW,NHE,NVN,NVS)
      call end_timer('EROSION',t1)
C
C     UPDATE ALL SOIL STATE VARIABLES FOR WATER, HEAT, GAS, SOLUTE
C     AND SEDIMENT FLUXES IN 'REDIST'
C
      if(lverb)WRITE(*,333)'RED'
      call start_timer(t1)
      CALL REDIST(I,J,NHW,NHE,NVN,NVS)
      call end_timer('RED',t1)
C     WRITE(*,333)'END'
C
C     WRITE HOURLY SOIL AND PLANT OUTPUT IN 'OUTSH' AND 'OUTPH'
C
      IF((J/JOUT)*JOUT.EQ.J)THEN
      if(lverb)WRITE(*,333)'OUTSH'
      CALL OUTSH(I,J,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'OUTPH'
      CALL OUTPH(I,J,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      ENDIF
C
C     WRITE OUTPUT FOR DYNAMIC VISUALIZATION
C
      IF(DATA(18).EQ.'YES')THEN
      IF((J/JOUT)*JOUT.EQ.J)THEN
      if(lverb)WRITE(*,333)'VIS'
      CALL VISUAL(I,J,NHW,NHE,NVN,NVS)
      ENDIF
      ENDIF
      call end_timer_loop()

9995  CONTINUE
      IF(DATA(19).EQ.'YES'.AND.KOUT.GT.0)THEN
C
C     WRITE ALL SOIL AND PLANT STATE VARIABLES AND OTHER INFORMATION
C     NEEDED TO RE-INITIALIZE THE MODEL TO CHECKPOINT FILES
C     IN 'WOUTS', 'WOUTP' AND 'WOUTQ'
C
      IF((I/KOUT)*KOUT.EQ.I.OR.I.EQ.IFIN)THEN
      if(lverb)WRITE(*,333)'WOUTS'
      CALL WOUTS(I,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'WOUTP'
      CALL WOUTP(I,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'WOUTQ'
      CALL WOUTQ(I,NHW,NHE,NVN,NVS)
      ENDIF
      ENDIF
C
C     WRITE DAILY SOIL AND PLANT OUTPUT IN 'OUTSD' AND 'OUTPD'
C
      IF((I/IOUT)*IOUT.EQ.I)THEN
      if(lverb)WRITE(*,333)'OUTSD'
      CALL OUTSD(I,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'OUTPD'
      CALL OUTPD(I,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      ENDIF
C
C     PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
C
      if(lverb)WRITE(*,333)'EXEC'
      CALL EXEC(I)
C
C     RE-INITIALIZE MODEL FROM CHECKPOINT FILES IF NEEDED
C
      IF(NYR.NE.1.AND.IDAYR.NE.IOLD)THEN
C     reinitialize the model from checkfiles
      if(lverb)WRITE(*,333)'ROUTS'
      CALL ROUTS(NHW,NHE,NVN,NVS)
      if(lverb)WRITE(*,333)'ROUTP'
      CALL ROUTP(NHW,NHE,NVN,NVS)
      ENDIF
      GO TO 9000
9999  CONTINUE
C     WRITE(*,333)'LOOP'
C
C     WRITE OUTPUT FILES FOR EACH GRID CELL IN 'SPLIT'
C
Cifdef _WIN_
      CALL SPLIT(NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
Celse
C     CALL SPLITC(NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
Cendif
C     WRITE(*,333)'SPLIT'
      RETURN
      END
