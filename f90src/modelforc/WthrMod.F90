module WthrMod
  !
  ! Description:
  ! code to process the weather forcing
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use MiniMathMod, only : safe_adb,vapsat0,test_aneb,test_aeqb
  use EcosimConst
  use CanopyRadDataType
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use LandSurfDataType
  use FertilizerDataType
  use PlantTraitDataType
  use AqueChemDatatype
  use IrrigationDataType
  use GridDataType
  use EcoSIMConfig
  use EcosimConst, only : TWILGT
  use MiniMathMod, only : AZMAX1
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__

  real(r8) :: AMP,CLD,DTA,DHR,DTS,EMM,RADX,RADZ,VPX,XJ

  !
  !     CDIR,CDIF=fraction of solar SW,sky diffuse radiation in visible
  !     PDIR,PDIF=PAR:SW ratio (umol s-1/(MJ h-1))
  !     TSNOW=temperature below which precipitation is snow (oC)
  !
  real(r8), PARAMETER :: CDIR=0.42,CDIF=0.58,PDIR=1269.4,PDIF=1269.4
  real(r8), PARAMETER :: TSNOW=-0.25

  public :: wthr
  contains

  SUBROUTINE wthr(I,J,NHW,NHE,NVN,NVS)
  !
  !     Description:
  !
  !     REINITIALIZES WEATHER VARIABLES USED IN OTHER
  !     SUBROUTINES
  !
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: ITYPE,NX,NY,N,NZ
  real(r8) :: PRECRI(JY,JX)
  real(r8) :: PRECWI(JY,JX)
  real(r8) :: PRECII(JY,JX)
  real(r8) :: PRECUI(JY,JX)
  real(r8) :: RADN(JY,JX)
  real(r8) :: VPS(JY,JX)

  !     execution begins here

  XJ=J
  DOY=I-1+XJ/24
  !
  !     SWITCH OUT ECOSYS WEATHER HERE IF CESM WEATHER IS READ IN
  !
  !     IF
  !
  !     IWTHR=weather data type in first(1) or second(2) scene
  !     ITYPE=weather data type:1=daily,2=hourly
  !
  IF(is_first_year.OR.I.LE.ILAST)THEN
    ITYPE=IWTHR(1)
  ELSE
    ITYPE=IWTHR(2)
  ENDIF
  !
  !     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
  !     AND PRECIPITATION FROM DAILY WEATHER ARRAYS LOADED IN 'READS'
  !
  IF(ITYPE.EQ.1)THEN
!
    call DailyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PRECRI,PRECWI,VPS)
    !     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
    !     AND PRECIPITATION FROM HOURLY WEATHER ARRAYS LOADED IN 'READS'
    !
  ELSE
    call HourlyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PRECRI,PRECWI,VPS)
  ENDIF
!
  call CalcRadiation(I,J,NHW,NHE,NVN,NVS,RADN,PRECUI,PRECII)
!

  IF(ICLM.EQ.1.OR.ICLM.EQ.2)THEN
    call ApplyClimateCorrection(I,J,NHW,NHE,NVN,NVS,PRECUI,&
      PRECRI,PRECII,PRECWI,VPS)
  ENDIF
!

  call SummaryForOutput(NHW,NHE,NVN,NVS,PRECUI,PRECRI,PRECII,PRECWI)

  END subroutine wthr
!------------------------------------------------------------------------------------------

  subroutine DailyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PRECRI,PRECWI,VPS)
!
!     Description:
!
  implicit none
  integer,  intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RADN(JY,JX),PRECRI(JY,JX),PRECWI(JY,JX),VPS(JY,JX)
  integer :: NY,NX
  !     begin_execution

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !     IETYP=Koppen climate zone:-2=phytotron
      !     RADN=hourky SW radiation
      !     RMAX=maximum hourly radiation
      !     ZNOON=time of solar noon
      !     DYLN=daylength
!
      IF(IETYP(NY,NX).NE.-2)THEN
        IF(DYLN(NY,NX).GT.ZERO)THEN
          RADN(NY,NX)=AZMAX1(RMAX*SIN((J-(ZNOON(NY,NX) &
            -DYLN(NY,NX)/2.0_r8))*PICON/DYLN(NY,NX)))
        ELSE
          RADN(NY,NX)=0.0_r8
        ENDIF
      ELSE
        RADN(NY,NX)=RMAX/24.0_r8
      ENDIF
      !
      !     TCA,TKA=air temperature (oC,K)
      !     TAVG*,AMP*=daily averages, amplitudes from day.f
      !
      IF(J.LT.(ZNOON(NY,NX)-DYLN(NY,NX)/2))THEN
        TCA(NY,NX)=TAVG1+AMP1*SIN(((J+ZNOON(NY,NX)-3.0)*PICON &
          /(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+PICON2)
      ELSEIF(J.GT.ZNOON(NY,NX)+3)THEN
        TCA(NY,NX)=TAVG3+AMP3*SIN(((J-ZNOON(NY,NX)-3.0)*PICON &
          /(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+PICON2)
      ELSE
        TCA(NY,NX)=TAVG2+AMP2*SIN(((J-(ZNOON(NY,NX) &
          -DYLN(NY,NX)/2.0))*PICON/(3.0+DYLN(NY,NX)/2.0))-PICON2)
      ENDIF
      TKA(NY,NX)=TCA(NY,NX)+TC2K
      !
      !     VPK,VPS=ambient,saturated vapor pressure
      !     VAVG*,VAMP*=daily averages, amplitudes from day.f
      !      ALTI=altitude
      !
      IF(J.LT.(ZNOON(NY,NX)-DYLN(NY,NX)/2))THEN
        VPK(NY,NX)=VAVG1+VMP1*SIN(((J+ZNOON(NY,NX)-3.0)*PICON &
          /(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+PICON2)
      ELSEIF(J.GT.ZNOON(NY,NX)+3)THEN
        VPK(NY,NX)=VAVG3+VMP3*SIN(((J-ZNOON(NY,NX)-3.0)*PICON &
          /(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+PICON2)
      ELSE
        VPK(NY,NX)=VAVG2+VMP2*SIN(((J-(ZNOON(NY,NX) &
          -DYLN(NY,NX)/2.0))*PICON /(3.0+DYLN(NY,NX)/2.0))-PICON2)
      ENDIF
      !VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TKA(NY,NX))) &
      VPS(NY,NX)=vapsat0(tka(ny,nx))*EXP(-ALTI(NY,NX)/7272.0)
      VPK(NY,NX)=AMIN1(VPS(NY,NX),VPK(NY,NX))
!
      !     UA=wind speed
!
      UA(NY,NX)=AMAX1(3600.0,WIND(I))
!
!     TSNOW=temperature below which precipitation is snow (oC)
!     PRECRI,PRECWI=rainfall,snowfall
!
      IF(J.GE.13.AND.J.LE.24)THEN
        IF(TCA(NY,NX).GT.TSNOW)THEN
          PRECRI(NY,NX)=RAIN(I)/12.0
          PRECWI(NY,NX)=0.0_r8
        ELSE
          PRECRI(NY,NX)=0.0_r8
          PRECWI(NY,NX)=RAIN(I)/12.0
    !     IF(PRECWI(NY,NX).LT.0.25E-03)PRECWI(NY,NX)=0.0_r8
        ENDIF
      ELSE
        PRECRI(NY,NX)=0.0_r8
        PRECWI(NY,NX)=0.0_r8
      ENDIF
    ENDDO
  enddo
!
  end subroutine DailyWeather
!------------------------------------------------------------------------------------------

  subroutine HourlyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PRECRI,PRECWI,VPS)
!
!     Description:
!
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RADN(JY,JX),PRECRI(JY,JX),PRECWI(JY,JX)
  real(r8), intent(out) :: VPS(JY,JX)
  integer :: NY,NX
  !     begin_execution

  DO  NX=NHW,NHE
    DO NY=NVN,NVS
 !
      !     RADN=SW radiation at horizontal surface
      !     TCA,TKA=air temperature (oC,K)
      !     VPK,VPS=ambient,saturated vapor pressure
      !     UA=wind speed
      !     TSNOW=temperature below which precipitation is snow (oC)
      !     PRECRI,PRECWI=rainfall,snowfall
!
      RADN(NY,NX)=SRADH(J,I)
      TCA(NY,NX)=TMPH(J,I)
      TKA(NY,NX)=TCA(NY,NX)+TC2K
      vps(ny,nx)=vapsat0(tka(ny,nx))*EXP(-ALTI(NY,NX)/7272.0_r8)
      VPK(NY,NX)=AMIN1(DWPTH(J,I),VPS(NY,NX))
      UA(NY,NX)=AMAX1(3600.0,WINDH(J,I))
      IF(TCA(NY,NX).GT.TSNOW)THEN
        PRECRI(NY,NX)=RAINH(J,I)
        PRECWI(NY,NX)=0.0_r8
      ELSE
        PRECRI(NY,NX)=0.0_r8
        PRECWI(NY,NX)=RAINH(J,I)
      ENDIF
    enddo
  enddo
  end subroutine HourlyWeather
!------------------------------------------------------------------------------------------

  subroutine CalcRadiation(I,J,NHW,NHE,NVN,NVS,RADN,PRECUI,PRECII)
!
  !     Description:
!
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: RADN(JY,JX)
  real(r8), intent(out) :: PRECUI(JY,JX),PRECII(JY,JX)
  integer :: NY,NX
  real(r8) :: AZI  !solar azimuth
  REAL(R8) :: DEC  !solar declination
  !     begin_execution
  !     CALCULATE DIRECT, DIFFUSE AND LONGWAVE RADIATION FROM
  !     INCOMING RADIATION READ IN 'READS', SOLAR ANGLE, HUMIDITY,
  !     TEMPERATURE AND CLOUDINESS
  !
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     IF OUTDOORS
!
!     SSIN,SSINN=sine solar angle of current,next hour
!     RADX=solar constant at horizontal surface
!     RADN=SW radiation at horizontal surface
!
      IF(IETYP(NY,NX).GE.-1)THEN
        AZI=SIN(ALAT(NY,NX)*1.7453E-02_r8)*SIN(DECLIN*1.7453E-02_r8)
        DEC=COS(ALAT(NY,NX)*1.7453E-02_r8)*COS(DECLIN*1.7453E-02_r8)

        SSIN(NY,NX)=AZMAX1(AZI+DEC*COS(.2618_r8*(ZNOON(NY,NX)-(J-0.5_r8))))
        SSINN(NY,NX)=AZMAX1(AZI+DEC*COS(.2618_r8*(ZNOON(NY,NX)-(J+0.5_r8))))

        !     IF(SSIN(NY,NX).GT.0.0.AND.SSIN(NY,NX).LT.TWILGT)SSIN(NY,NX)=TWILGT
        IF(RADN(NY,NX).LE.0.0)SSIN(NY,NX)=0.0_r8
        IF(SSIN(NY,NX).LE.-TWILGT)RADN(NY,NX)=0.0_r8
        RADX=4.896_r8*AZMAX1(SSIN(NY,NX))
        RADN(NY,NX)=AMIN1(RADX,RADN(NY,NX))
!
        !     DIRECT VS DIFFUSE RADIATION IN SOLAR OR SKY BEAMS
        !
        !     RADZ=diffuse radiation at horizontal surface
        !     RADS,RADY,RAPS,RAPY=direct,diffuse SW,PAR in solar beam
!

        RADZ=AMIN1(RADN(NY,NX),0.5_r8*(RADX-RADN(NY,NX)))
        RADS(NY,NX)=safe_adb(RADN(NY,NX)-RADZ,SSIN(NY,NX))
        RADS(NY,NX)=AMIN1(4.167_r8,RADS(NY,NX))
        RADY(NY,NX)=RADZ/TYSIN
        RAPS(NY,NX)=RADS(NY,NX)*CDIR*PDIR
        RAPY(NY,NX)=RADY(NY,NX)*CDIF*PDIF
        !
        !     ATMOSPHERIC RADIATIVE PROPERTIES AFM 139:171
        !
        !     CLD=cloudiness factor for EMM
        !     EMM=sky emissivity
        !     VPK,TKA=vapor pressure,temperature
        !
        IF(RADX.GT.ZERO)THEN
          CLD=AMIN1(1.0_r8,AMAX1(0.2_r8,2.33_r8-3.33_r8*RADN(NY,NX)/RADX))
        ELSE
          CLD=0.2_r8
        ENDIF
        EMM=0.625_r8*AMAX1(1.0_r8,(1.0E+03_r8*VPK(NY,NX)/TKA(NY,NX))**0.131_r8)
        EMM=EMM*(1.0_r8+0.242_r8*CLD**0.583_r8)
        !
        !     IF PHYTOTRON
        !
      ELSE
        IF(RADN(NY,NX).LE.0.0_r8)THEN
          SSIN(NY,NX)=0.0_r8
        ELSE
          SSIN(NY,NX)=1.0_r8
        ENDIF
        SSINN(NY,NX)=1.0_r8
        CLD=0.0_r8
        EMM=0.96
      ENDIF
!
      !     LONGWAVE RADIATION
      !
      !     XRADH=longwave radiation
      !     THSX=longwave radiation from weather file or calculated from
      !     atmospheric properties
!
      IF(XRADH(J,I).GT.0.0)THEN
        !     THSX(NY,NX)=EMM*(2.04E-10*TKA(NY,NX)**4)
        !     THSX(NY,NX)=THSX(NY,NX)+XRADH(J,I)
        THSX(NY,NX)=XRADH(J,I)
      ELSE
        THSX(NY,NX)=EMM*(2.04E-10_r8*TKA(NY,NX)**4._r8)
      ENDIF
!
      !     INSERT CESM WEATHER HERE
      !
      !     ELSE
      !     RADS=DIRECT SW RADIATION (MJ M-2 H-1)
      !     RADY=INDIRECT SW RADIATION (MJ M-2 H-1)
      !     RAPS=DIRECT PAR (UMOL M-2 S-1)
      !     RAPY=INDIRECT PAR (UMOL M-2 S-1)
      !     THSX=LW RADIATION (MJ M-2 H-1)
      !     TCA=AIR TEMPERATURE (C)
      !     TKA=AIR TEMPERATURE (K)
      !     VPK=VAPOR PRESSURE (KPA)
      !     UA=WINDSPEED (M H-1)
      !     PRECRI(NY,NX)=RAIN (M H-1)
      !     PRECWI(NY,NX)=SNOW (M H-1)
      !     SSIN=SOLAR ANGLE CURRENT HOUR (SINE)
      !     SSINN=SOLAR ANGLE NEXT HOUR (SINE)
      !     ENDIF
      !
      !     ADD IRRIGATION
      !
      !     PRECII,PRECUI=surface,subsurface irrigation
      !     RRIG=irrigation from soil management file in reads.f
      !
      WDPTHD=WDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
!     IF(WDPTHD.LE.CDPTH(NU(NY,NX),NY,NX))THEN
      PRECII(NY,NX)=RRIG(J,I,NY,NX)
      PRECUI(NY,NX)=0.0_r8
!     ELSE
!     PRECII(NY,NX)=0.0_r8
!     PRECUI(NY,NX)=RRIG(J,I,NY,NX)
!     ENDIF
    ENDDO
  ENDDO
  end subroutine CalcRadiation
!------------------------------------------------------------------------------------------

  subroutine ApplyClimateCorrection(I,J,NHW,NHE,NVN,NVS,PRECUI,&
    PRECRI,PRECII,PRECWI,VPS)
  !
  !     DESCRIPTION:
  !     IMPLEMENT CLIMATE CHANGES READ IN 'READS' TO HOURLY TEMPERATURE,
  !     RADIATION, WINDSPEED,VAPOR PRESSURE, PRECIPITATION, IRRIGATION
  !     AND CO2
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: PRECUI(JY,JX),PRECRI(JY,JX),PRECII(JY,JX)
  real(r8), intent(inout) :: PRECWI(JY,JX),VPS(JY,JX)
  integer :: NY,NX,NZ,N
  !     begin_execution
  !
  !     ICLM=changes to weather data (0=none,1=step,2=transient)
  !     N=season or month
  !
  !
  !     SEASONAL CHANGES
  !
  IF(I.GT.334.OR.I.LE.59)THEN
    N=1
  ELSEIF(I.GT.59.AND.I.LE.151)THEN
    N=2
  ELSEIF(I.GT.151.AND.I.LE.243)THEN
    N=3
  ELSE
    N=4
  ENDIF
  !
  !     MONTHLY CHANGES
  !
  !     IF(I.GT.0.AND.I.LE.31)THEN
  !     N=1
  !     ELSEIF(I.GT.31.AND.I.LE.59)THEN
  !     N=2
  !     ELSEIF(I.GT.59.AND.I.LE.90)THEN
  !     N=3
  !     ELSEIF(I.GT.90.AND.I.LE.120)THEN
  !     N=4
  !     ELSEIF(I.GT.120.AND.I.LE.151)THEN
  !     N=5
  !     ELSEIF(I.GT.151.AND.I.LE.181)THEN
  !     N=6
  !     ELSEIF(I.GT.181.AND.I.LE.212)THEN
  !     N=7
  !     ELSEIF(I.GT.212.AND.I.LE.243)THEN
  !     N=8
  !     ELSEIF(I.GT.243.AND.I.LE.273)THEN
  !     N=9
  !     ELSEIF(I.GT.273.AND.I.LE.304)THEN
  !     N=10
  !     ELSEIF(I.GT.304.AND.I.LE.334)THEN
  !     N=11
  !     ELSE
  !     N=12
  !     ENDIF
  DO 9925 NX=NHW,NHE
    DO 9920 NY=NVN,NVS
!
      !     TEMPERATURE CHANGE
      !
      !     TDTPX,TDTPN=change in max,min temperature
      !     DTA,AMP,DHR=change in daily average,amplitude of air temperature
      !     DHR=diurnal effect on AMP
      !
      IF(test_aneb(TDTPX(NY,NX,N),0.0_r8).OR.test_aneb(TDTPN(NY,NX,N),0.0_r8))THEN
        DTA=0.5_r8*(TDTPX(NY,NX,N)+TDTPN(NY,NX,N))
        AMP=0.5_r8*(TDTPX(NY,NX,N)-TDTPN(NY,NX,N))
        DHR=SIN(0.2618_r8*(J-(ZNOON(NY,NX)+3.0_r8))+PICON2)
        TCA(NY,NX)=TCA(NY,NX)+DTA+AMP*DHR
        TKA(NY,NX)=TCA(NY,NX)+TC2K
!
        !     ACCLIMATION TO GRADUAL CLIMATE CHANGE
        !
        !     DTS=change in daily average soil temperature
        !     ATCA,ATCS=mean annual air,soil temperature
        !     OFFSET=shift in Arrhenius curve for MFT activity in nitro.f
        !     OFFST=shift in Arrhenius curve for PFT activity in uptake.f
        !     ZTYP=PFT thermal adaptation zone
        !     HTC=high temperature threshold for grain number loss (oC)
        !     GROUPI,XTLI=node number at floral initiation,planting (maturity group)
!
        IF(ICLM.EQ.2.AND.J.EQ.1)THEN
          DTS=0.5_r8*DTA
          ATCA(NY,NX)=ATCAI(NY,NX)+DTA
          ATCS(NY,NX)=ATCAI(NY,NX)+DTS
          OFFSET(NY,NX)=0.33*(12.5-AZMAX1(AMIN1(25.0,ATCS(NY,NX))))
          DO NZ=1,NP(NY,NX)
            ZTYP(NZ,NY,NX)=ZTYPI(NZ,NY,NX)+0.30/2.667*DTA
            OFFST(NZ,NY,NX)=2.667*(2.5-ZTYP(NZ,NY,NX))
            !     TCZ(NZ,NY,NX)=TCZD-OFFST(NZ,NY,NX)
            !     TCX(NZ,NY,NX)=AMIN1(15.0,TCZ(NZ,NY,NX)+TCXD)
            IF(ICTYP(NZ,NY,NX).EQ.3)THEN
              HTC(NZ,NY,NX)=27.0+3.0*ZTYP(NZ,NY,NX)
            ELSE
              HTC(NZ,NY,NX)=30.0+3.0*ZTYP(NZ,NY,NX)
            ENDIF
            GROUPI(NZ,NY,NX)=GROUPX(NZ,NY,NX)+0.30*DTA
            IF(IBTYP(NZ,NY,NX).NE.0)THEN
              GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)/25.0
            ENDIF
            GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)-XTLI(NZ,NY,NX)

          ENDDO
        ENDIF
!
        !     ADJUST VAPOR PRESSURE FOR TEMPERATURE CHANGE
!
        IF(test_aeqb(DHUM(N),1.0_r8))THEN
          VPX=VPS(NY,NX)
          !VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TKA(NY,NX))) &
          vps(ny,ny)=vapsat0(tka(ny,nx))*EXP(-ALTI(NY,NX)/7272.0)
          VPK(NY,NX)=VPK(NY,NX)*VPS(NY,NX)/VPX
        ENDIF
      ENDIF
!
      !     CHANGES IN VAPOR PRESSURE,RADIATION,WIND SPEED,PRECIPITATION
      !     IRRIGATION ,CO2,NH4,NO3
      !
      !     TDRAD,TDWND,TDHUM=change in radiation,windspeed,vapor pressure
      !     TDPRC,TDIRRI=change in precipitation,irrigation
      !     TDCO2,TDCN4,TDCNO=change in atm CO2,NH4,NO3 concn in precipitation
!
      RADS(NY,NX)=RADS(NY,NX)*TDRAD(NY,NX,N)
      RADY(NY,NX)=RADY(NY,NX)*TDRAD(NY,NX,N)
      RAPS(NY,NX)=RAPS(NY,NX)*TDRAD(NY,NX,N)
      RAPY(NY,NX)=RAPY(NY,NX)*TDRAD(NY,NX,N)
      UA(NY,NX)=UA(NY,NX)*TDWND(NY,NX,N)
      VPK(NY,NX)=AMIN1(VPS(NY,NX),VPK(NY,NX)*TDHUM(NY,NX,N))
      PRECRI(NY,NX)=PRECRI(NY,NX)*TDPRC(NY,NX,N)
      PRECWI(NY,NX)=PRECWI(NY,NX)*TDPRC(NY,NX,N)
      PRECII(NY,NX)=PRECII(NY,NX)*TDIRI(NY,NX,N)
      PRECUI(NY,NX)=PRECUI(NY,NX)*TDIRI(NY,NX,N)
      CO2E(NY,NX)=CO2EI(NY,NX)*TDCO2(NY,NX,N)
      CN4R(NY,NX)=CN4RI(NY,NX)*TDCN4(NY,NX,N)
      CNOR(NY,NX)=CNORI(NY,NX)*TDCNO(NY,NX,N)
9920  CONTINUE
9925  CONTINUE
  end subroutine ApplyClimateCorrection
!------------------------------------------------------------------------------------------

  subroutine SummaryForOutput(NHW,NHE,NVN,NVS,PRECUI,PRECRI,PRECII,PRECWI)
  !
  !     DESCRIPTION:
  !     DAILY WEATHER TOTALS, MAXIMA AND MINIMA FOR DAILY OUTPUT
  !     CHECK AGAINST INPUTS FROM WEATHER FILE
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(in):: PRECUI(JY,JX),PRECRI(JY,JX)
  real(r8), intent(in) :: PRECII(JY,JX),PRECWI(JY,JX)
  integer :: NY,NX
  !
  !     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(SSIN(NY,NX).GT.0.0_r8)TRAD(NY,NX)=TRAD(NY,NX)+RADS(NY,NX) &
        *SSIN(NY,NX)+RADY(NY,NX)*TYSIN
      TAMX(NY,NX)=AMAX1(TAMX(NY,NX),TCA(NY,NX))
      TAMN(NY,NX)=AMIN1(TAMN(NY,NX),TCA(NY,NX))
      HUDX(NY,NX)=AMAX1(HUDX(NY,NX),VPK(NY,NX))
      HUDN(NY,NX)=AMIN1(HUDN(NY,NX),VPK(NY,NX))
      TWIND(NY,NX)=TWIND(NY,NX)+UA(NY,NX)
      VPA(NY,NX)=VPK(NY,NX)*2.173E-03_r8/TKA(NY,NX)
      TRAI(NY,NX)=TRAI(NY,NX)+(PRECRI(NY,NX)+PRECWI(NY,NX) &
        +PRECII(NY,NX)+PRECUI(NY,NX))*1000.0_r8
!
      !     WATER AND HEAT INPUTS TO GRID CELLS
      !
      !     AREA=area of grid cell (m2)
      !
      !     PRECR,PRECW=rainfall,snowfall
      !     PRECI,PRECU=surface,subsurface irrigation
      !     PRECA,PRECQ=rain+irrigation,rain+snow
      !     THS=sky LW radiation
!
      PRECR(NY,NX)=PRECRI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECW(NY,NX)=PRECWI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECI(NY,NX)=PRECII(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECU(NY,NX)=PRECUI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECA(NY,NX)=PRECR(NY,NX)+PRECI(NY,NX)
      PRECQ(NY,NX)=PRECR(NY,NX)+PRECW(NY,NX)
      THS(NY,NX)=THSX(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
    ENDDO
  ENDDO
  end subroutine SummaryForOutput

end module WthrMod
