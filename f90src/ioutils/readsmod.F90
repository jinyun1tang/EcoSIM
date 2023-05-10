module readsmod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils , only : endrun
  use fileUtil   , only : open_safe
  use minimathmod, only : isLeap,AZMAX1
  use EcoSIMCtrlMod, only : lverb
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

  SUBROUTINE reads(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)
!
! THIS SUBROUTINE READS ALL SOIL AND PLANT MANAGEMENT INPUT FILES
!
  use ReadManagementMod, only : ReadManagementFiles
  use EcoSIMCtrlMod, only : soil_mgmt_in
  implicit none
  integer, intent(in) :: yearc   !current model year
  integer, intent(in) :: yeari   !current data year
  integer, intent(in) :: NE,NEX,NHW,NHE,NVN,NVS

  CHARACTER(len=1) :: TTYPE
  integer :: kk,N,L,NY,NX,NZ,K
  integer :: LL,J
  integer :: LPY,IX
  type(atm_forc_type) :: atmf

  call readCLMfactors(yeari)

! begin_execution

!
! OPTIONS(4, AND LAND MANAGEMENT(9, FILES FROM
! FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!
! PREFIX=path for files in current or higher level directory
!
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

! generating checkpoint files,resuming from earlier checkpt files
! DRAD,DTMPX,DTMPN,DHUM,DPREC,DIRRI,DWIND,DCO2E,DCNR4,DCNOR
! =annual changes in radiation,max+min temperature,humidity,
! precip,irrign,windspeed,atm CO2 concn,NH4,NO3 concn in precip
! NPX=number of cycles per hour for water,heat,solute flux calcns
! NPY=number of cycles per NPX for gas flux calcns
! JOUT,IOUT,KOUT=output frequency for hourly,daily,checkpoint data
! ICLM=changes to weather data (0=none,1=step,2=transient)

! this reads the option file that sets up frequency of model output and # of
! iterations used by the ecosim solvers

! determine whether to read checkpoing file (i.e. an actual restart run)
! read changing factor for climate variables
! read from netcdf file

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
    write(*,*)'number of cycles per hour for water, heat, and '// &
      'solute flux calcns: NPX ',NPX
    write(*,*)'number of cycles per NPX for gas flux calcns: NPY',NPY
    write(*,*)'output frequency for hourly data: JOUT ',JOUT
    write(*,*)'output frequency for daily data: IOUT ',IOUT
    write(*,*)'output frequency for checkpoint data: KOUT ',KOUT
    write(*,*)'changes to weather data (0=none,1=step,'// &
      '2=transient): ICLM ',ICLM
  endif

!
! INCREMENTS IN START AND END DATES FOR SUCCESSIVE SCENARIOS
! FROM LOOPS FOR SCENES, SCENARIOS IN RUNSCRIPT SET IN MAIN.F
!
! IDATA(3),IDATA(6)=start,end year of current scene
!
  IYRC=yearc

!
! OPEN CHECKPOINT FILES FOR SOIL VARIABLES
!
! IDATE=year label for checkpoint files
! DATA(1)=site file name
! W,N=water+heat,nutrient checkpoint files
!
  IF(is_first_year)THEN
    idate=yearc

    WRITE(CHARY,'(I4)')IDATE
    OUTW='W'//DATA1(1)(1:2)//CHARY(1:4)
    OUTN='N'//DATA1(1)(1:2)//CHARY(1:4)
    OPEN(21,FILE=trim(outdir)//OUTW,STATUS='UNKNOWN')
    OPEN(22,FILE=trim(outdir)//OUTN,STATUS='UNKNOWN')
  ENDIF

  IF(is_first_year)THEN
    L=1
  ELSE
    L=2
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

  call ReadClimNC(yearc,yeari,L,atmf)

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
      ROWN(NY,NX)=0.0_r8
      ROWO(NY,NX)=0.0_r8
      ROWP(NY,NX)=0.0_r8
      D325: DO I=1,366
        ITILL(I,NY,NX)=0
        DCORP(I,NY,NX)=0.0_r8
      ENDDO D325
      D40: DO I=1,366
        D45: DO N=1,20
          FERT(N,I,NY,NX)=0.0_r8
        ENDDO D45
        D35: DO N=0,2
          IYTYP(N,I,NY,NX)=0
        ENDDO D35
        FDPTH(I,NY,NX)=0.0_r8
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
  IF(trim(soil_mgmt_in).NE.'NO')THEN
    call ReadManagementFiles(yeari)
  ENDIF
  IMNG=1
  RETURN
  END subroutine reads
!------------------------------------------------------------------------------------------
  subroutine readCLMfactors(yeari)
  !
  !DESCRIPTION
  !read climate change factors
  use EcoSIMCtrlMod, only : clm_factor_in
  use netcdf
  use ncdio_pio
  use fileUtil, only : file_exists
  implicit none
  integer, intent(in) :: yeari
  integer :: iyear,year
  type(file_desc_t) :: clm_factor_nfid
!  type(Var_desc_t) :: vardesc
!  logical :: readvar
  INTEGER :: N

  if(.not.file_exists(clm_factor_in))then
    call endrun("Do not find clm_factor_in file "//trim(clm_factor_in)//" in "//mod_filename, __LINE__)
  endif
  call ncd_pio_openfile(clm_factor_nfid, clm_factor_in, ncd_nowrite)
  
  iyear=1
  DO while(.true.)
    call ncd_getvar(clm_factor_nfid,'year',iyear,year)  
    if(year==yeari)exit
    iyear=iyear+1
  ENDDO
  call ncd_getvar(clm_factor_nfid,'DRAD',iyear,DRAD(1))
  call ncd_getvar(clm_factor_nfid,'DTMPX',iyear,DTMPX(1))
  call ncd_getvar(clm_factor_nfid,'DTMPN',iyear,DTMPN(1))
  call ncd_getvar(clm_factor_nfid,'DHUM',iyear,DHUM(1))
  call ncd_getvar(clm_factor_nfid,'DPREC',iyear,DPREC(1))
  call ncd_getvar(clm_factor_nfid,'DIRRI',iyear,DIRRI(1))
  call ncd_getvar(clm_factor_nfid,'DWIND',iyear,DWIND(1))
  call ncd_getvar(clm_factor_nfid,'DCO2E',iyear,DCO2E(1))
  call ncd_getvar(clm_factor_nfid,'DCN4R',iyear,DCN4R(1))
  call ncd_getvar(clm_factor_nfid,'DCNOR',iyear,DCNOR(1))
  call ncd_getvar(clm_factor_nfid,'ICLM',iyear,ICLM)

  call ncd_pio_closefile(clm_factor_nfid)

  DO N=2,12
    DRAD(N)=DRAD(1)
    DTMPX(N)=DTMPX(1)
    DTMPN(N)=DTMPN(1)
    DHUM(N)=DHUM(1)
    DPREC(N)=DPREC(1)
    DIRRI(N)=DIRRI(1)
    DWIND(N)=DWIND(1)
    DCO2E(N)=DCO2E(1)
    DCN4R(N)=DCN4R(1)
    DCNOR(N)=DCNOR(1)
  ENDDO

  end subroutine readCLMfactors
end module readsmod
