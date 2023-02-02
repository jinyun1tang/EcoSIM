module readsmod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils , only : endrun
  use fileUtil   , only : open_safe
  use minimathmod, only : isLeap,AZMAX1
  use GridConsts
  use FlagDataType
  use FertilizerDataType
  use ClimForcDataType
  use SoilWaterDataType
  use LandSurfDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use IrrigationDataType
  use GridDataType
  use EcoSIMConfig
  implicit none
  private

  character(len=*), parameter :: mod_filename = __FILE__
  integer, SAVE :: N1,N2,N1X,N2X,IFLGY,IYRX,IYRD

  real(r8) :: CN4RIG,CNORIG,CN4RG,CNORG,CPORG,CALRG
  real(r8) :: CFERG,CCARG,CMGRG,CNARG,CKARG,CSORG,CCLRG
  real(r8) :: Z0G,ZNOONG
  real(r8) :: PHRG

  CHARACTER(len=1) :: IVAR(20),VAR(50),TYP(50),CTYPE
  CHARACTER(len=16) :: OUTW,OUTI,OUTT,OUTN,OUTF
  CHARACTER(len=4) :: CHARY
  integer :: IDAT(20)
  real(r8) :: DAT(50),DATK(50),OUT(50)
  real(r8) :: datav(40)
  DATA IFLGY,IYRX,IYRD/0,0,0/

  integer :: IDATE,IDY,IFLG3,IX,I,ICHECK



  public :: reads
  contains

  SUBROUTINE reads(NA,ND,NT,NE,NAX,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)
!
! THIS SUBROUTINE READS ALL SOIL AND PLANT MANAGEMENT INPUT FILES
!
  use ReadManagementMod, only : ReadManagementFiles
  implicit none
  integer, intent(in) :: NEX
  integer, intent(in) :: NA(1:NEX),ND(1:NEX)
  integer, intent(in) :: NT,NE,NAX,NTX,NHW,NHE,NVN,NVS
  integer, intent(inout) :: NF, NFX, NTZ
  integer, intent(out) :: NTZX

  CHARACTER(len=1) :: TTYPE
  integer :: kk,N,L,NY,NX,NZ,K
  integer :: LL,J
  integer :: LPY


! begin_execution

!
! OPTIONS(4, AND LAND MANAGEMENT(9, FILES FROM
! FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!
! PREFIX=path for files in current or higher level directory
!
  call OPEN_safe(4,PREFIX,DATAC(4,NE,NEX),'OLD',mod_filename,__LINE__)

!
! ARTIFICIAL SOIL WARMING
!
! soiltemp=file with hourly soil temperatures from baseline run
! OUT=hourly soil temperatures from baseline run (oC)
! TKSZ=temperature used to calculate additional heat flux
! for warming in watsub.f
!
! OPEN(6,FILE='soiltemp',STATUS='OLD')
!23   READ(6,'(F8.3,4X,A8,I8,50E16.7E3)',END=27)DOY,CDATE,J,2,(OUT(L),L=1,33)
! IF(J.EQ.24)THEN
!   I=INT(DOY)
! ELSE
!   I=INT(DOY)+1
! ENDIF
! DO 24 L=1,20
!   TKSZ(I,J,L)=OUT(L+13)+4.0+TC2K
!24    CONTINUE
! GO TO 23
!27    CONTINUE
!
! READ START AND END DATES, WEATHER OPTIONS
!
! IDATA(1),IDATA(2),IDATA(3)=start date of scenario DDMMYYYY
! IDATA(4),IDATA(5),IDATA(6)=end date of scenario DDMMYYYY
! IDATA(7),IDATA(8),IDATA(9)=start date of run DDMMYYYY
! DATA1(18),DATA1(19),DATA1(20)=options for visualization in visual.f
! generating checkpoint files,resuming from earlier checkpt files
! DRAD,DTMPX,DTMPN,DHUM,DPREC,DIRRI,DWIND,DCO2E,DCNR4,DCNOR
!   =annual changes in radiation,max+min temperature,humidity,
! precip,irrign,windspeed,atm CO2 concn,NH4,NO3 concn in precip
! NPX=number of cycles per hour for water,heat,solute flux calcns
! NPY=number of cycles per NPX for gas flux calcns
! JOUT,IOUT,KOUT=output frequency for hourly,daily,checkpoint data
! ICLM=changes to weather data (0=none,1=step,2=transient)

! this reads the option file that sets up frequency of model output and # of
! iterations used by the ecosim solvers

  READ(4,'(2I2,I4)')IDATA(1),IDATA(2),IDATA(3)       !beginning of the current year
  READ(4,'(2I2,I4)')IDATA(4),IDATA(5),IDATA(6)       !end of current year
  READ(4,'(2I2,I4)')IDATA(7),IDATA(8),IDATA(9)       !beginning of the past year
  READ(4,'(A3)')DATA1(18)
  READ(4,'(A3)')DATA1(19)
  READ(4,'(A3)')DATA1(20)
  if(lverb)then
    write(*,'(100A)')('-',ll=1,100)
    write(*,*)'read option file ',DATAC(4,NE,NEX)
    write(*,*)'start date of scenario DDMMYYYY: IDATA(1:3) ' &
      ,IDATA(1),IDATA(2),IDATA(3)
    write(*,*)'end date of scenario DDMMYYYY: IDATA(4:6) ' &
      ,IDATA(4),IDATA(5),IDATA(6)
    write(*,*)'start date of run DDMMYYYY: IDATA(7:9) ' &
      ,IDATA(7),IDATA(8),IDATA(9)
    write(*,*)'ouput hourly visualization in visual.f?: '// &
      'DATA1(18) ',DATA1(18)
    write(*,*)'write checkpoint file?: DATA1(19) ',DATA1(19)
    write(*,*)'read chkpt file?: DATA1(20) ',DATA1(20)
  endif
! determine whether to read checkpoing file (i.e. an actual restart run)
  is_restart_run=(data1(20)=='YES')
! read changing factor for climate variables
  D25: DO N=1,4
    READ(4,*)DRAD(N),DTMPX(N),DTMPN(N),DHUM(N),DPREC(N) &
      ,DIRRI(N),DWIND(N),DCO2E(N),DCN4R(N),DCNOR(N)
  ENDDO D25

  if(lverb)then
    write(*,*)'annual changes in radiation: DRAD(1:4)',DRAD(1:4)
    write(*,*)'annual changes in max temperature: DTMPX(1:4)',DTMPX(1:4)
    write(*,*)'annual changes in min temperature: DTMPN(1:4)',DTMPN(1:4)
    write(*,*)'annual changes in humidity: DHUM(1:4)',DHUM(1:4)
    write(*,*)'annual changes in precip: DPREC(1:4)',DPREC(1:4)
    write(*,*)'annual changes in irrigation: DIRRI(1:4)',DIRRI(1:4)
    write(*,*)'annual changes in wind speed: DWIND(1:4)',DWIND(1:4)
    write(*,*)'annual changes in atm CO2 conc: DCO2E(1:4)',DCO2E(1:4)
    write(*,*)'annual changes in atm NH4 conc: DCN4R(1:4)',DCN4R(1:4)
    write(*,*)'annual changes in atm NO3 conc: DCNOR(1:4)',DCNOR(1:4)
  endif

  D26: DO N=5,12
    DRAD(N)=DRAD(N-1)
    DTMPX(N)=DTMPX(N-1)
    DTMPN(N)=DTMPN(N-1)
    DHUM(N)=DHUM(N-1)
    DPREC(N)=DPREC(N-1)
    DIRRI(N)=DIRRI(N-1)
    DWIND(N)=DWIND(N-1)
    DCO2E(N)=DCO2E(N-1)
    DCN4R(N)=DCN4R(N-1)
    DCNOR(N)=DCNOR(N-1)
  ENDDO D26
  READ(4,*)NPX,NPY,JOUT,IOUT,KOUT,ICLM
  if(lverb)then
    write(*,*)'number of cycles per hour for water, heat, and '// &
      'solute flux calcns: NPX ',NPX
    write(*,*)'number of cycles per NPX for gas flux calcns: NPY',NPY
    write(*,*)'output frequency for hourly data: JOUT ',JOUT
    write(*,*)'output frequency for daily data: IOUT ',IOUT
    write(*,*)'output frequency for checkpoint data: KOUT ',KOUT
    write(*,*)'changes to weather data (0=none,1=step,'// &
      '2=transient): ICLM ',ICLM
  endif
  CLOSE(4)

!
! INCREMENTS IN START AND END DATES FOR SUCCESSIVE SCENARIOS
! FROM LOOPS FOR SCENES, SCENARIOS IN RUNSCRIPT SET IN MAIN.F
!
! IDATA(3),IDATA(6)=start,end year of current scene
!
  NTZX=NTZ
  IF(is_first_year.OR.IDATA(3).NE.0)THEN
    IDATA(3)=IDATA(3)+(NT-1)*NF+(NTX-1)*NFX-NTZX
    IDATA(6)=IDATA(6)+(NT-1)*NF+(NTX-1)*NFX-NTZX
    IYRC=IDATA(3)
  ELSE
    !not the first year
    IF(IDATA(1).EQ.1.AND.IDATA(2).EQ.1)THEN
      !jan 1st
      IDATA(3)=IYRC+1
    ELSE
      IDATA(3)=IYRC
    ENDIF
    IDATA(6)=IDATA(3)
    IYRC=IDATA(3)
  ENDIF

  IF(NE.EQ.1)THEN
    N1=IDATA(3)
  ENDIF
  IF(NE.EQ.NA(NEX))THEN
    N2=IDATA(6)
    NF=N2-N1+1
    IF(IDATA(4).NE.31.OR.IDATA(5).NE.12)THEN
      NTZ=NTZ+1
    ENDIF
  ENDIF
  IF(NE.EQ.1.AND.NT.EQ.1.AND.NEX.EQ.1)THEN
    N1X=IDATA(3)
  ENDIF
  IF(NE.EQ.NA(NEX).AND.NT.EQ.ND(NEX).AND.NEX.EQ.NAX)THEN
    N2X=IDATA(6)
    NFX=N2X-N1X+1
    IF(NE.NE.NA(NEX))THEN
      IF(IDATA(4).NE.31.OR.IDATA(5).NE.12)THEN
        NTZ=NTZ+1
      ENDIF
    ENDIF
  ENDIF
! WRITE(*,7766)'IDATA3',IGO,IDATA(3),IDATA(6),IYRR,IYRC
!    2,NE,NT,NEX,NF,NTX,NFX,NTZ,NTZX,N1,N2,N1X,N2X
!    3,NA(NEX),ND(NEX),NAX
!7766  FORMAT(A8,30I8)
!
! OPEN CHECKPOINT FILES FOR SOIL VARIABLES
!
! IDATE=year label for checkpoint files
! DATA(1)=site file name
! W,N=water+heat,nutrient checkpoint files
!

  IF(is_first_year)THEN
    IF(is_restart_run)THEN
      IDATE=IDATA(9)
    ELSE
      IDATE=IDATA(3)
    ENDIF
    WRITE(CHARY,'(I4)')IDATE
    OUTW='W'//DATA1(1)(1:2)//CHARY(1:4)
    OUTN='N'//DATA1(1)(1:2)//CHARY(1:4)
    OPEN(21,FILE=trim(outdir)//OUTW,STATUS='UNKNOWN')
    OPEN(22,FILE=trim(outdir)//OUTN,STATUS='UNKNOWN')
  ENDIF
!
! CALCULATE START AND FINISH DATES IN JULIAN DAYS
! FROM DATE INPUTS IN OPTION FILE
!
! ISTART,IBEGIN=start dates for current scene
! IFIN,ILAST=end dates for current scene
! LYRC=number of days in current year
!
  LPY=0
  LYRC=365   ! # days for current year
  LYRX=365   ! # days for last year
!
  D575: DO N=1,7,3
    IF(isLeap(IDATA(N+2)))then
      IF(IDATA(N+1).GT.2)LPY=1
      IF(N.EQ.1)LYRC=366
    endif

    IF(IDATA(N+1).EQ.1)then
! Jan
      IDY=IDATA(N)
    else
! total ordinal days counted till month IDATA(N+1)
      IDY=30*(IDATA(N+1)-1)+ICOR(IDATA(N+1)-1)+IDATA(N)+LPY
    endif

    IF(N.EQ.1)ISTART=IDY      !
    IF(N.EQ.4)IFIN=IDY        !
    IF(N.EQ.7)IRUN=IDY        !starting date of the restart run

    IF(isLeap(IDATA(N+2)-1))then
      IF(N.EQ.1)LYRX=366
    endif
  ENDDO D575

  IF(is_first_year)THEN
    IF(.NOT.is_restart_run)IRUN=ISTART
    L=1
    ILAST=ISTART-1
    ITERM=IFIN
  ELSE
    L=2
    ILAST=MIN(ISTART-1,ITERM,IEND)
    ITERM=IFIN
  ENDIF

!
! READ WEATHER DATA
!
! DATAC(3=weather file
! TTYPE,CTYPE=time step,calendar format
! NI,NN=number of time,weather data variables
! IVAR,VAR=time,weather variable type
! TYP=weather variable units
! Z0G,IFLGW=windspeed meast height,flag for raising Z0G with vegn
! ZNOONG=time of solar noon
! PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG,CKARG,
! CSORG,CCLRG=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
! concentration in precipitation
! IDAT,DAT=time,weather variable
!
  IF(DATAC(3,NE,NEX).NE.'NO')THEN

    call ReadClim(NE,NEX,NTX,NFX,L,I,IX,TTYPE)

!
! ACCOUNT FOR LEAP YEAR
!
    IF(I.EQ.365)THEN
      IF(TTYPE.EQ.'D')THEN
! daily
        TMPX(I+1)=TMPX(I)
        TMPN(I+1)=TMPN(I)
        SRAD(I+1)=SRAD(I)
        WIND(I+1)=WIND(I)
        DWPT(1,I+1)=DWPT(1,I)
        DWPT(2,I+1)=DWPT(2,I)
        RAIN(I+1)=RAIN(I)
      ELSE
        D130: DO J=1,24
          TMPH(J,I+1)=TMPH(J,I)
          SRADH(J,I+1)=SRADH(J,I)
          WINDH(J,I+1)=WINDH(J,I)
          DWPTH(J,I+1)=DWPTH(J,I)
          RAINH(J,I+1)=RAINH(J,I)
          XRADH(J,I+1)=XRADH(J,I)
        ENDDO D130
      ENDIF
      IX=I+1
    ENDIF

  ELSE
    IFLGW=1
    Z0G=2.0_r8
    ZNOONG=12.0_r8
    PHRG=7.0_r8
    CN4RIG=0.0_r8
    CNORIG=0.0_r8
    CN4RG=CN4RIG
    CNORG=CNORIG
    CPORG=0.0_r8
    CALRG=0.0_r8
    CFERG=0.0_r8
    CCARG=0.0_r8
    CMGRG=0.0_r8
    CNARG=0.0_r8
    CKARG=0.0_r8
    CSORG=0.0_r8
    CCLRG=0.0_r8
    IX=365
  ENDIF
!
! CALCULATE PRECIPITATION CONCENTRATIONS IN MOLE UNITS
!
  CN4RIG=CN4RIG/14.0_r8
  CNORIG=CNORIG/14.0_r8
  CN4RG=CN4RIG
  CNORG=CNORIG
  CPORG=CPORG/31.0_r8
  CALRG=CALRG/27.0_r8
  CFERG=CFERG/55.8_r8
  CCARG=CCARG/40.0_r8
  CMGRG=CMGRG/24.3_r8
  CNARG=CNARG/23.0_r8
  CKARG=CKARG/39.1_r8
  CSORG=CSORG/32.0_r8
  CCLRG=CCLRG/35.5_r8

  D8970: DO NX=NHW,NHE
    D8975: DO NY=NVN,NVS
      Z0(NY,NX)=Z0G         !windspeed meast height
      ZNOON(NY,NX)=ZNOONG
      PHR(NY,NX)=PHRG
      CN4RI(NY,NX)=CN4RIG
      CNORI(NY,NX)=CNORIG
      CN4R(NY,NX)=CN4RIG
      CNOR(NY,NX)=CNORIG
      CPOR(NY,NX)=CPORG
      CALR(NY,NX)=CALRG
      CFER(NY,NX)=CFERG
      CCAR(NY,NX)=CCARG
      CMGR(NY,NX)=CMGRG
      CNAR(NY,NX)=CNARG
      CKAR(NY,NX)=CKARG
      CSOR(NY,NX)=CSORG
      CCLR(NY,NX)=CCLRG
    ENDDO D8975
  ENDDO D8970
!
! DERIVE END DATES FROM TIME VARIABLES
!

  ICHECK=0
  IF(TTYPE.EQ.'H'.AND.J.NE.24)ICHECK=1
  IEND=IX-ICHECK
  IFIN=MIN(IFIN,IEND)
  IDAYR=MIN(ISTART-1,ILAST) !day of recovery from earlier run
  IYRR=IDATA(3)
  NYR=0
  IF(IDAYR.EQ.0)THEN
    IDAYR=LYRX
    IYRR=IDATA(3)-1
    NYR=1
  ENDIF
  IFLGY=0
!
!     READ LAND MANAGEMENT FILE NAMES FOR EACH GRID CELL
!
  D9980: DO NX=NHW,NHE
    D9985: DO NY=NVN,NVS
      ROWN(NY,NX)=0.0
      ROWO(NY,NX)=0.0
      ROWP(NY,NX)=0.0
      D325: DO I=1,366
        ITILL(I,NY,NX)=0
        DCORP(I,NY,NX)=0.0
      ENDDO D325
      D40: DO I=1,366
        D45: DO N=1,20
          FERT(N,I,NY,NX)=0.0
        ENDDO D45
        D35: DO N=0,2
          IYTYP(N,I,NY,NX)=0
        ENDDO D35
        FDPTH(I,NY,NX)=0.0
      ENDDO D40
      D125: DO I=1,366
        D120: DO J=1,24
          RRIG(J,I,NY,NX)=0.0_r8
        ENDDO D120
        PHQ(I,NY,NX)=7.0_r8
        CN4Q(I,NY,NX)=0.0_r8
        CNOQ(I,NY,NX)=0.0_r8
        CPOQ(I,NY,NX)=0.0_r8
        CALQ(I,NY,NX)=0.0_r8
        CFEQ(I,NY,NX)=0.0_r8
        CCAQ(I,NY,NX)=0.0_r8
        CMGQ(I,NY,NX)=0.0_r8
        CNAQ(I,NY,NX)=0.0_r8
        CKAQ(I,NY,NX)=0.0_r8
        CSOQ(I,NY,NX)=0.0_r8
        CCLQ(I,NY,NX)=0.0_r8
        WDPTH(I,NY,NX)=0.0_r8
        ROWI(I,NY,NX)=0.0_r8
      ENDDO D125
    ENDDO D9985
  ENDDO D9980
!
! READ LAND MANAGEMENT FILE DATAC(9 LOADED IN 'MAIN'.
! THIS FILE CONTAINS NAMES OF TILLAGE, IRRIGATION
! AND FERTILIZER FILES
!
  IF(DATAC(9,NE,NEX).NE.'NO')THEN
!
    call ReadManagementFiles(DATAC(9,NE,NEX))

  ENDIF
  IMNG=1
  RETURN
  END subroutine reads
!------------------------------------------------------------------------------------------

  subroutine interp3hourweather(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: II,JJ
  JJ=J-3
  II=I
  IF(JJ.EQ.0)THEN
    JJ=24
    II=I-1
  ENDIF
!
! INFILL 3-HOURLY WEATHER DATA
!
  IF(II.LT.1)THEN
    TMPH(J-2,I)=TMPH(J,I)
    TMPH(J-1,I)=TMPH(J,I)
    SRADH(J-2,I)=SRADH(J,I)
    SRADH(J-1,I)=SRADH(J,I)
    WINDH(J-2,I)=WINDH(J,I)
    WINDH(J-1,I)=WINDH(J,I)
    DWPTH(J-2,I)=DWPTH(J,I)
    DWPTH(J-1,I)=DWPTH(J,I)
    RAINH(J,I)=RAINH(J,I)/3.0
    RAINH(J-2,I)=RAINH(J,I)
    RAINH(J-1,I)=RAINH(J,I)
    XRADH(J-2,I)=XRADH(J,I)
    XRADH(J-1,I)=XRADH(J,I)
  ELSE
    TMPH(J-2,I)=0.667*TMPH(JJ,II)+0.333*TMPH(J,I)
    TMPH(J-1,I)=0.333*TMPH(JJ,II)+0.667*TMPH(J,I)
    SRADH(J-2,I)=0.667*SRADH(JJ,II)+0.333*SRADH(J,I)
    SRADH(J-1,I)=0.333*SRADH(JJ,II)+0.667*SRADH(J,I)
    WINDH(J-2,I)=0.667*WINDH(JJ,II)+0.333*WINDH(J,I)
    WINDH(J-1,I)=0.333*WINDH(JJ,II)+0.667*WINDH(J,I)
    DWPTH(J-2,I)=0.667*DWPTH(JJ,II)+0.333*DWPTH(J,I)
    DWPTH(J-1,I)=0.333*DWPTH(JJ,II)+0.667*DWPTH(J,I)
    RAINH(J,I)=RAINH(J,I)/3.0_r8
    RAINH(J-2,I)=RAINH(J,I)
    RAINH(J-1,I)=RAINH(J,I)
    XRADH(J-2,I)=0.667*XRADH(JJ,II)+0.333*XRADH(J,I)
    XRADH(J-1,I)=0.333*XRADH(JJ,II)+0.667*XRADH(J,I)
  ENDIF
  end subroutine interp3hourweather

!------------------------------------------------------------------------------------------

  subroutine readhourweather(I,J,IH,go60,L,NN,NI,TTYPE)

!     read hourly weather data
!
!     DERIVE DAY I AND HOUR J FROM TIME VARIABLES IVAR
!
  implicit none
  integer, intent(inout) :: I      !julian day
  integer, intent(inout) :: J      !hour
  integer, intent(inout) :: IH

  integer, intent(in) :: L,NN,NI
  character(len=1),intent(in) :: TTYPE
  logical, intent(out) :: go60

  integer :: LPY
  integer :: K
  integer :: M
  integer :: N

  go60=.false.

  IWTHR(L)=2
  D190: DO K=1,NI
    IF(IVAR(K).EQ.'M')THEN
     !  month
      M=IDAT(K)
    ELSEIF(IVAR(K).EQ.'D')THEN
     ! julian day
      N=IDAT(K)
    ELSEIF(IVAR(K).EQ.'H')THEN
     ! hour
      J=IDAT(K)
    ENDIF
    IF(IVAR(K).EQ.'Y')THEN
      IFLGY=1
    ENDIF
  ENDDO D190

!  write(*,'(4(A,X,I4))')'IFLGY=',IFLGY,' IYRX=',IYRX,&
!    ' IYRC=',IYRC,' J=',J
   IF(IFLGY.EQ.1.AND.IYRX.LT.IYRC)then
     GO60=.true.
     return
   endif
   IF(CTYPE.EQ.'J')THEN
     I=N
   ELSE
     LPY=0
     IF(isLeap(IDATA(3)).and.M.GT.2)LPY=1
     IF(M.EQ.1)THEN
       I=N
     ELSE
       I=30*(M-1)+ICOR(M-1)+N+LPY
     ENDIF
   ENDIF

   IF(J.GT.24.AND.(J/100)*100.NE.J)THEN
     D80: DO K=1,NN
       DATK(K)=DATK(K)+DAT(K)
     ENDDO D80
     IH=IH+1
     go60=.true.
     return
   ENDIF

   IF(J.GT.24)J=INT(J/100)
     IF(J.EQ.0)THEN
       J=24
       I=I-1
       IF(I.LT.1)then
         GO60=.true.
         return
        endif
      ENDIF
!
!     DERIVE START DATE FROM TIME VARIABLES
!
      IF(IFLG3.EQ.0)THEN
        IBEGIN=N
        ISTART=MAX(ISTART,IBEGIN)
        IFLG3=1
      ENDIF

      IF(L.NE.1)THEN
        IF(I.LE.ILAST)then
          GO60=.true.
          return
        endif
      ENDIF
      XRADH(J,I)=0.0
!
!     CONVERT HOURLY WEATHER VARIABLES TO MODEL UNITS
!     AND ENTER INTO MODEL ARRAYS
!
!     TMPH=temperature (oC)
!     SRADH=solar radiation (MJ m-2 h-1)
!     WINDH=windspeed (m h-1)
!     DWPTH=vapor pressure (kPa)
!     RAINH=precipitation (m h-1)
!     XRADH=longwave radiation (MJ m-2 h-1)
!
      D95: DO K=1,NN
!
!       TEMPERATURE
!
        IF(VAR(K).EQ.'T')THEN
          IF(TYP(K).EQ.'F')THEN
! Fahrenheit to celcius
            TMPH(J,I)=((DAT(K)+DATK(K))/IH-32.0_r8)*0.556_r8
          ELSEIF(TYP(K).EQ.'K')THEN
! Temperature given as kelvin
            TMPH(J,I)=(DAT(K)+DATK(K))/IH-TC2K
          ELSE
            TMPH(J,I)=(DAT(K)+DATK(K))/IH
          ENDIF
!
!         SOLAR RADIATION
!
        ELSEIF(VAR(K).EQ.'R')THEN
! shortwave radiation
          IF(TYP(K).EQ.'W')THEN
! W m-2 to MJ m-2 per hour (3600 seconds per hour * 1.e-6)
            SRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.0036_r8)
          ELSEIF(TYP(K).EQ.'J')THEN
! 1.e4 J m-2
            SRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.01_r8)
          ELSEIF(TYP(K).EQ.'K')THEN
! kJ m-2 per hour
            SRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.001_r8)
          ELSEIF(TYP(K).EQ.'P')THEN
            SRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.0036_r8*0.457_r8)
!     ELSEIF(TYP(K).EQ.'M')THEN
!     SRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*3.6*0.457)
          ELSE
            SRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)
          ENDIF
!
!     WIND SPEED
!
        ELSEIF(VAR(K).EQ.'W')THEN
          IF(TYP(K).EQ.'S')THEN
! given as m s-1, into m per hour
            WINDH(J,I)=(DAT(K)+DATK(K))/IH*3600.0_r8
          ELSEIF(TYP(K).EQ.'H')THEN
! km per hour into m per hour
            WINDH(J,I)=(DAT(K)+DATK(K))/IH*1000.0_r8
          ELSEIF(TYP(K).EQ.'M')THEN
! mile per hour into m per hour, it should be 1609.34 though
            WINDH(J,I)=(DAT(K)+DATK(K))/IH*1600.0_r8
          ELSE
            WINDH(J,I)=(DAT(K)+DATK(K))/IH
          ENDIF
        ELSEIF(VAR(K).EQ.'H')THEN
!
!       VAPOR PRESSURE
!
          IF(TYP(K).EQ.'D')THEN
! given as celcius degree
            DWPTH(J,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+(DAT(K)+DATK(K))/IH)))
          ELSEIF(TYP(K).EQ.'F')THEN
! given as Fahrenheit
            DAT(K)=(DAT(K)-32.0_r8)*0.556_r8
            DWPTH(J,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+(DAT(K)+DATK(K))/IH)))
          ELSEIF(TYP(K).EQ.'H')THEN
! given as relative humidity, [0, 1]
            DWPTH(J,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+TMPH(J,I)))) &
              *AZMAX1(AMIN1(1.0_r8,(DAT(K)+DATK(K))/IH))
          ELSEIF(TYP(K).EQ.'R')THEN
! given as relative humidity [0, 100]
            DWPTH(J,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+TMPH(J,I)))) &
              *AZMAX1(AMIN1(100.0_r8,(DAT(K)+DATK(K))/IH))*0.01_r8
          ELSEIF(TYP(K).EQ.'S')THEN
! what is the unit? mass mixing ratio? [0,100]
            DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)*0.0289_r8/18.0_r8*101.325_r8 &
              *EXP(-ALTIG/7272.0_r8)*288.15_r8/(TC2K+TMPH(J,I))
          ELSEIF(TYP(K).EQ.'G')THEN
! what is the unit? mass mixing ratio? [0, 1] ALTIG is doing altitude correction,
            DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)*28.9_r8/18.0_r8*101.325_r8 &
              *EXP(-ALTIG/7272.0_r8)*288.15_r8/(TC2K+TMPH(J,I))
          ELSEIF(TYP(K).EQ.'M')THEN
! given as hPa
            DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.1_r8)
          ELSE
            DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)
          ENDIF
!
!       PRECIPITATION, [m h-1]
!
      ELSEIF(VAR(K).EQ.'P')THEN
        IF(TYP(K).EQ.'M')THEN
!   [mm m-2] (mm to m, 1.e-3)
          RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))/1.0E+03
        ELSEIF(TYP(K).EQ.'C')THEN
!   [cm m-2]
          RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))/1.0E+02
        ELSEIF(TYP(K).EQ.'I')THEN
!  [inch m-2] -> [m m-2]
          RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))*0.0254
        ELSEIF(TYP(K).EQ.'S')THEN
!  [second m-2]
          IF(TTYPE.EQ.'H')THEN
! given in hour
            RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))*3.6
          ELSE
! given in half hour
            RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))*1.8
          ENDIF
        ELSE
          RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))
        ENDIF
!
!       LONGWAVE RADIATION (OPTIONAL)
!
      ELSEIF(VAR(K).EQ.'L')THEN
        IF(TYP(K).EQ.'W')THEN
          ! watss m-2, into MJ per hour, *3600 seconds * 1.0-6 = 0.0036
          XRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.0036)
        ELSEIF(TYP(K).EQ.'J')THEN
          XRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.01)
        ELSEIF(TYP(K).EQ.'K')THEN
          ! kJ m-2, into MJ per hour, * 0.001
          XRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.001)
        ELSE
          XRADH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)
        ENDIF
      ENDIF
      DATK(K)=0.0
  ENDDO D95
  end subroutine readhourweather
!------------------------------------------------------------------------------------------

  subroutine readdayweather(I,L,go110,go60,NTX,NFX,NN,NI)
!     read daily weather data
!
!     DERIVE DAY I FROM TIME VARIABLES IVAR
!
!     IWTHR=weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene
!
  implicit none

  integer, intent(in) :: L,NTX,NFX,NN,NI

  integer, intent(out) :: I
  logical, intent(out) :: go110, go60
  integer :: LPY
  integer :: K,M,N

  go110=.false.
  go60=.false.
  IWTHR(L)=1
  D160: DO K=1,NI
    IF(IVAR(K).EQ.'M')THEN
      M=IDAT(K)
    ELSEIF(IVAR(K).EQ.'D')THEN
      N=IDAT(K)
    ENDIF
    IF(IVAR(K).EQ.'Y')THEN
      IFLGY=1
      IYRX=IDAT(K)+(NTX-1)*NFX
      IF(isLeap(IDAT(K)))then
        IYRD=366
      else
        IYRD=365
      endif
    ENDIF
  ENDDO D160
  IF(IFLGY.EQ.1.AND.IYRX.LT.IYRC)THEN
    GO60=.TRUE.
    RETURN
  ENDIF
! CTYPE for calendar format
  IF(CTYPE.EQ.'J')THEN
    I=N
  ELSE
    LPY=0
    IF(isLeap(IDATA(3)).and.M.GT.2)LPY=1
    IF(M.EQ.1)THEN
      I=N
    ELSE
      I=30*(M-1)+ICOR(M-1)+N+LPY
    ENDIF
  ENDIF
!
! DERIVE START DATE FROM TIME VARIABLES
!
  IF(IFLG3.EQ.0)THEN
    IBEGIN=I
    ISTART=MAX(ISTART,IBEGIN)
    IFLG3=1
  ENDIF
  IF(L.NE.1)THEN
    IF(I.LE.ILAST)THEN
      GO60=.TRUE.
    ENDIF
  ENDIF
!
! CONVERT DAILY WEATHER VARIABLES TO MODEL UNITS
! AND ENTER INTO MODEL ARRAYS
!
! TMPX,TMPN=maximum,minimum temperature (OC)
! SRAD=solar radiation (MJ m-2 d-1)
! WIND=windspeed (m h-1)
! DWPT=vapor pressure (kPa)
! RAIN=precipitation (mm d-1)
!
  D65: DO K=1,NN
!
!   MAX,MIN REMPERATURE
!
    IF(VAR(K).EQ.'M')THEN
      IF(TYP(K).EQ.'F')THEN
        TMPX(I)=(DAT(K)-32.0)*0.556
      ELSEIF(TYP(K).EQ.'K')THEN
        TMPX(I)=DAT(K)-273.16
      ELSE
        TMPX(I)=DAT(K)
      ENDIF
    ELSEIF(VAR(K).EQ.'N')THEN
      IF(TYP(K).EQ.'F')THEN
        TMPN(I)=(DAT(K)-32.0)*0.556
      ELSEIF(TYP(K).EQ.'K')THEN
        TMPN(I)=DAT(K)-273.16
      ELSE
        TMPN(I)=DAT(K)
      ENDIF
!
!     SOLAR RADIATION
!
    ELSEIF(VAR(K).EQ.'R')THEN
      IF(TYP(K).EQ.'L')THEN
        SRAD(I)=AZMAX1(DAT(K)/23.87)
      ELSEIF(TYP(K).EQ.'J')THEN
        SRAD(I)=AZMAX1(DAT(K)*0.01)
      ELSE
        SRAD(I)=AZMAX1(DAT(K))
      ENDIF
!
!     WIND SPEED
!
    ELSEIF(VAR(K).EQ.'W')THEN
      IF(TYP(K).EQ.'S')THEN
        WIND(I)=ABS(DAT(K))*3600.0
      ELSEIF(TYP(K).EQ.'H')THEN
        WIND(I)=ABS(DAT(K))*1000.0
      ELSEIF(TYP(K).EQ.'D')THEN
        WIND(I)=ABS(DAT(K))*1000.0/24.0
      ELSEIF(TYP(K).EQ.'M')THEN
        WIND(I)=ABS(DAT(K))*1600.0
      ELSE
        WIND(I)=ABS(DAT(K))
      ENDIF
!
!     VAPOR PRESSURE
!
    ELSEIF(VAR(K).EQ.'H')THEN
      IF(TYP(K).EQ.'D')THEN
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+DAT(K))))
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+DAT(K))))
      ELSEIF(TYP(K).EQ.'F')THEN
        DAT(K)=(DAT(K)-32.0_r8)*0.556_r8
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+DAT(K))))
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+DAT(K))))
      ELSEIF(TYP(K).EQ.'H')THEN
        DAT(K)=AZMAX1(AMIN1(1.0_r8,DAT(K)))
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+(TMPN(I)+TMPX(I))/2)))*DAT(K)
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+TMPN(I))))
      ELSEIF(TYP(K).EQ.'R')THEN
        DAT(K)=AZMAX1(AMIN1(100.0_r8,DAT(K)))
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+(TMPN(I)+TMPX(I))/2)))*DAT(K)*0.01_r8
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/(TC2K+TMPN(I))))
      ELSEIF(TYP(K).EQ.'S')THEN
        DWPT(1,I)=AZMAX1(DAT(K))*0.0289_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/(TC2K+(TMPN(I)+TMPX(I))/2)
        DWPT(2,I)=AZMAX1(DAT(K))*0.0289_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/(TC2K+TMPN(I))
      ELSEIF(TYP(K).EQ.'G')THEN
        DWPT(1,I)=AZMAX1(DAT(K))*28.9_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/(TC2K+(TMPN(I)+TMPX(I))/2._r8)
        DWPT(2,I)=AZMAX1(DAT(K))*28.9_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/(TC2K+TMPN(I))
      ELSEIF(TYP(K).EQ.'M')THEN
        DWPT(1,I)=AZMAX1(DAT(K)*0.1_r8)
        DWPT(2,I)=AZMAX1(DAT(K)*0.1_r8)
      ELSE
        DWPT(1,I)=AZMAX1(DAT(K))
        DWPT(2,I)=AZMAX1(DAT(K))
      ENDIF
!
!     PRECIPITATION
!
    ELSEIF(VAR(K).EQ.'P')THEN
      IF(TYP(K).EQ.'M')THEN
        RAIN(I)=AZMAX1(DAT(K))/1.0E+03_r8
      ELSEIF(TYP(K).EQ.'C')THEN
        RAIN(I)=AZMAX1(DAT(K))/1.0E+02_r8
      ELSEIF(TYP(K).EQ.'I')THEN
        RAIN(I)=AZMAX1(DAT(K))*0.0254_r8
      ELSE
        RAIN(I)=AZMAX1(DAT(K))
      ENDIF
    ENDIF
  ENDDO D65
  IX=I
  IF(IFLGY.EQ.1.AND.I.EQ.IYRD)THEN
    GO110=.true.
    return
  ENDIF
  GO60=.TRUE.
  end subroutine readdayweather

!------------------------------------------------------------------------------------------

  function getclimttype(ttype)result(ans)
!  return climate data time step
  implicit none
  character(len=1), intent(in):: ttype
  character(len=16) :: ans

  select case (ttype)
  case ('D')
    ans='daily'
  case ('H')
    ans='hourly'
  case ('3')
    ans='3-hourly'
  end select
  end function getclimttype

!------------------------------------------------------------------------------------------

  function getclimctype(ctype)result(ans)
! return climate calendar type
  implicit none
  character(len=1), intent(in) :: ctype
  character(len=16) :: ans

  select case(ctype)
  case ('J')
    ans='Julian day'
  case default
    ans='not julian day'
  end select
  end function getclimctype

!------------------------------------------------------------------------------------------

  subroutine ReadClim(NE,NEX,NTX,NFX,L,I,IX,TTYPE)
  implicit none
  integer, intent(in) :: NE,NEX,NTX,NFX,L
  integer, intent(out) :: I,IX

  CHARACTER(len=1),intent(out) :: TTYPE
  integer :: K, KK, IH, NI, NN, J
  integer :: LL
  LOGICAL :: GO110,GO60

! OPEN WEATHER(3,
  call OPEN_safe(3,PREFIX,DATAC(3,NE,NEX),'OLD',mod_filename,__LINE__)
  IFLG3=0
  READ(3,'(2A1,2I2,50A1)')TTYPE,CTYPE,NI,NN,(IVAR(K),K=1,NI),(VAR(K),K=1,NN)
  READ(3,'(50A1)')(TYP(K),K=1,NN)
  read(3,*)(datav(kk),kk=1,3)
  Z0G=datav(1)
  IFLGW=int(datav(2))
  ZNOONG=datav(3)

!  fourth line in the weather file
  READ(3,*)PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG,CKARG,CSORG,CCLRG
  if(lverb)then
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'read weather file head from ',DATAC(3,NE,NEX)
    write(*,*)'time step format: TTYPE ',TTYPE,trim(getclimttype(TTYPE))
    write(*,*)'calendar format: CTYPE ',CTYPE,trim(getclimctype(CTYPE))
    write(*,*)'number of time variables: NI ',NI
    write(*,*)'time var type: IVAR ',(IVAR(K),K=1,NI)
    write(*,*)'number of weather data variables: NN ',NN
    write(*,*)'weather var type: VAR ',(VAR(K),K=1,NN)
    write(*,*)'weather var units: TYP ',(TYP(K),K=1,NN)
    write(*,*)'windspeed measurement height: Z0G',Z0G
    write(*,*)'flag for raising Z0G with vegn: IFLGW ',IFLGW
    write(*,*)'time of solar noon: ZNOONG ',ZNOONG
    write(*,*)'pH in precipitation: PHRG ',PHRG
    write(*,*)'NH4 conc in precip: CN4RIG ',CN4RIG
    write(*,*)'NO3 conc in precip: CNORIG ',CNORIG
    write(*,*)'H2PO4 conc in precip: CPORG ',CPORG
    write(*,*)'Al conc in precip: CALRG ',CALRG
    write(*,*)'Fe conc in precip: CFERG ',CFERG
    write(*,*)'Ca conc in precip: CCARG ',CCARG
    write(*,*)'Mg conc in precip: CMGRG ',CMGRG
    write(*,*)'Na conc in precip: CNARG ',CNARG
    write(*,*)'K conc in precip: CKARG ',CKARG
    write(*,*)'SO4 conc in precip: CSORG ',CSORG
    write(*,*)'Cl conc in precip: CCLRG ',CCLRG
    write(*,*)'weather data are in the format time,weather variable'
  endif
  D55: DO K=1,NN
    DATK(K)=0.0_r8
  ENDDO D55
  IH=1

! the file reading loop
  DO while(.TRUE.)
    read(3,*,END=110)(datav(k),k=1,NI),(DAT(K),K=1,NN)
    do k = 1, ni
      idat(k)=int(datav(k))
    enddo
!
    IF(TTYPE.EQ.'D')THEN
!   READ DAILY WEATHER DATA AND CONVERT TO MODEL UNITS
      call readdayweather(I,L,GO110,GO60,NTX,NFX,NN,NI)
      IF(GO60)cycle
      IF(GO110)EXIT
!!
    ELSE
!     READ HOURLY WEATHER DATA AND CONVERT TO MODEL UNITS
      call readhourweather(I,J,IH,go60,L,NN,NI,TTYPE)
!     write(*,*)'goto60=',go60
      IF(GO60)cycle
      IH=1
      IX=I
!     write(*,*)'344IX=',IX
      IF(TTYPE.EQ.'3')THEN
!       weather data is every 3 hrs, do interpolation
        call interp3hourweather(I,J)
      ENDIF

      IF(IFLGY.EQ.1.AND.I.EQ.IYRD.AND.J.EQ.24)THEN
        EXIT
      ENDIF
      cycle
    ENDIF
  enddo
110 CONTINUE
  CLOSE(3)
  end subroutine ReadClim
end module readsmod
