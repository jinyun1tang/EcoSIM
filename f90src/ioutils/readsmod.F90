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
  use ClimReadMod
  implicit none
  private

  character(len=*), parameter :: mod_filename = __FILE__
  integer, SAVE :: N1,N2,N1X,N2X

  real(r8) :: CN4RIG,CNORIG,CN4RG,CNORG,CPORG,CALRG
  real(r8) :: CFERG,CCARG,CMGRG,CNARG,CKARG,CSORG,CCLRG
  real(r8) :: Z0G,ZNOONG
  real(r8) :: PHRG

  CHARACTER(len=16) :: OUTW,OUTI,OUTT,OUTN,OUTF
  CHARACTER(len=4) :: CHARY

  integer :: IDATE,IDY,IFLG3,I,ICHECK



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
  integer :: LPY,IX
  type(atm_forc_type) :: atmf

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
    N1X=IDATA(3)   !beginning of current year
  ENDIF

  IF(NE.EQ.NA(NEX).AND.NT.EQ.ND(NEX).AND.NEX.EQ.NAX)THEN
    N2X=IDATA(6)   !end of current year
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
! Z0G,IFLGW=windspeed mean height,flag for raising Z0G with vegn
! ZNOONG=time of solar noon
! PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG,CKARG,
! CSORG,CCLRG=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
! concentration in precipitation
! IDAT,DAT=time,weather variable
!
  IF(DATAC(3,NE,NEX).NE.'NO')THEN

    call ReadClim(IDATA(3),DATAC(3,NE,NEX),NTX,NFX,L,I,IX,TTYPE,atmf)

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
    atmf%Z0G=2.0_r8
    atmf%ZNOONG=12.0_r8
    atmf%PHRG=7.0_r8
    atmf%CN4RIG=0.0_r8
    atmf%CNORIG=0.0_r8
    CN4RG=atmf%CN4RIG
    CNORG=atmf%CNORIG
    atmf%CPORG=0.0_r8
    atmf%CALRG=0.0_r8
    atmf%CFERG=0.0_r8
    atmf%CCARG=0.0_r8
    atmf%CMGRG=0.0_r8
    atmf%CNARG=0.0_r8
    atmf%CKARG=0.0_r8
    atmf%CSORG=0.0_r8
    atmf%CCLRG=0.0_r8
    IX=365
  ENDIF
!
! CALCULATE PRECIPITATION CONCENTRATIONS IN MOLE UNITS
!
  CN4RIG=atmf%CN4RIG/14.0_r8
  CNORIG=atmf%CNORIG/14.0_r8
  CN4RG=CN4RIG
  CNORG=CNORIG
  CPORG=atmf%CPORG/31.0_r8
  CALRG=atmf%CALRG/27.0_r8
  CFERG=atmf%CFERG/55.8_r8
  CCARG=atmf%CCARG/40.0_r8
  CMGRG=atmf%CMGRG/24.3_r8
  CNARG=atmf%CNARG/23.0_r8
  CKARG=atmf%CKARG/39.1_r8
  CSORG=atmf%CSORG/32.0_r8
  CCLRG=atmf%CCLRG/35.5_r8
  PHRG=atmf%PHRG
  Z0G=atmf%Z0G
  ZNOONG=atmf%ZNOONG
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

end module readsmod
