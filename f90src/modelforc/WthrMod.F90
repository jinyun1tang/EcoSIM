module WthrMod
  !
  ! Description:
  ! code to process the weather forcing
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use MiniMathMod  , only : safe_adb,vapsat0,isclose
  use MiniFuncMod  , only : get_sun_declin
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
  use MiniMathMod, only : AZMAX1
  use UnitMod    , only : units
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: AMP,CLD,DTA,DHR,DTS,EMM,RADX,RADZ,VPX,XJ

  !
  !     CDIR,CDIF=fraction of solar SW,sky diffuse radiation in visible
  !     PDIR,PDIF=PAR:SW ratio (umol s-1/(MJ h-1))
  !     TSNOW=temperature below which precipitation is snow (oC)
  !
  real(r8), parameter :: SolConst=4.896_r8 !MJ/m2/hour, solar constant
  real(r8), PARAMETER :: CDIR=0.42_r8   !
  real(r8), parameter :: CDIF=0.58_r8
  real(r8), parameter :: PDIR=1269.4_r8
  real(r8), parameter :: PDIF=1269.4_r8
  real(r8), PARAMETER :: TSNOW=-0.25_r8  !oC, threshold temperature for snowfall

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
  real(r8) :: PrecAsRain(JY,JX)
  real(r8) :: PrecAsSnow(JY,JX)
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
    call DailyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PrecAsRain,PrecAsSnow,VPS)
    !     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
    !     AND PRECIPITATION FROM HOURLY WEATHER ARRAYS LOADED IN 'READS'
    !
  ELSE
    call HourlyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PrecAsRain,PrecAsSnow,VPS)
  ENDIF
!
  call CalcRadiation(I,J,NHW,NHE,NVN,NVS,RADN,PRECUI,PRECII)
!

  IF(ICLM.EQ.1.OR.ICLM.EQ.2)THEN
    !climate data will be perturbed
    call CorrectClimate(I,J,NHW,NHE,NVN,NVS,PRECUI,PrecAsRain,PRECII,PrecAsSnow,VPS)
  ENDIF
!

  call SummaryForOutput(NHW,NHE,NVN,NVS,PRECUI,PrecAsRain,PRECII,PrecAsSnow)

  END subroutine wthr
!------------------------------------------------------------------------------------------

  subroutine DailyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PrecAsRain,PrecAsSnow,VPS)
!
!     Description:
!
  implicit none
  integer,  intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RADN(JY,JX),PrecAsRain(JY,JX)
  real(r8), intent(out) :: PrecAsSnow(JY,JX),VPS(JY,JX)
  integer :: NY,NX
  !     begin_execution

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !     KoppenClimZone=Koppen climate zone:-2=phytotron
      !     RADN=hourky SW radiation
      !     RMAX=maximum hourly radiation
      !     SolarNoonHour_col=time of solar noon
      !     DayLenthCurrent=daylength
!
      IF(KoppenClimZone(NY,NX).NE.-2)THEN
!       not phytotron
        IF(DayLenthCurrent(NY,NX).GT.ZERO)THEN
          RADN(NY,NX)=AZMAX1(RMAX*SIN((J-(SolarNoonHour_col(NY,NX) &
            -DayLenthCurrent(NY,NX)/2.0_r8))*PICON/DayLenthCurrent(NY,NX)))
        ELSE
          RADN(NY,NX)=0.0_r8
        ENDIF
      ELSE
!       phytotron
        RADN(NY,NX)=RMAX/24.0_r8
      ENDIF
      !
      !     TCA,TairK=air temperature (oC,K)
      !     TAVG*,AMP*=daily averages, amplitudes from day.f
      !
      IF(J.LT.(SolarNoonHour_col(NY,NX)-DayLenthCurrent(NY,NX)/2))THEN
        TCA(NY,NX)=TAVG1+AMP1*SIN(((J+SolarNoonHour_col(NY,NX)-3.0)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0-DayLenthCurrent(NY,NX)/2.0))+PICON2h)
      ELSEIF(J.GT.SolarNoonHour_col(NY,NX)+3)THEN
        TCA(NY,NX)=TAVG3+AMP3*SIN(((J-SolarNoonHour_col(NY,NX)-3.0)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0-DayLenthCurrent(NY,NX)/2.0))+PICON2h)
      ELSE
        TCA(NY,NX)=TAVG2+AMP2*SIN(((J-(SolarNoonHour_col(NY,NX) &
          -DayLenthCurrent(NY,NX)/2.0_r8))*PICON/(3.0_r8+DayLenthCurrent(NY,NX)/2.0_r8))-PICON2h)
      ENDIF
      TairK(NY,NX)=units%Celcius2Kelvin(TCA(NY,NX))
      if(abs(TairK(NY,NX))>400._r8)then
        print*,'air temperature problematic',TairK(NY,NX),TCA(NY,NX)
      endif
      !
      !     VPK,VPS=ambient,saturated vapor pressure
      !     VAVG*,VAMP*=daily averages, amplitudes from day.f
      !      ALTI=altitude
      !
      IF(J.LT.(SolarNoonHour_col(NY,NX)-DayLenthCurrent(NY,NX)/2))THEN
        VPK(NY,NX)=VAVG1+VMP1*SIN(((J+SolarNoonHour_col(NY,NX)-3.0_r8)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0_r8-DayLenthCurrent(NY,NX)/2.0_r8))+PICON2h)
      ELSEIF(J.GT.SolarNoonHour_col(NY,NX)+3)THEN
        VPK(NY,NX)=VAVG3+VMP3*SIN(((J-SolarNoonHour_col(NY,NX)-3.0_r8)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0_r8-DayLenthCurrent(NY,NX)/2.0_r8))+PICON2h)
      ELSE
        VPK(NY,NX)=VAVG2+VMP2*SIN(((J-(SolarNoonHour_col(NY,NX) &
          -DayLenthCurrent(NY,NX)/2.0_r8))*PICON /(3.0_r8+DayLenthCurrent(NY,NX)/2.0_r8))-PICON2h)
      ENDIF
      !VPS(NY,NX)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TairK(NY,NX))) &
      VPS(NY,NX)=vapsat0(TairK(ny,nx))*EXP(-ALTI(NY,NX)/7272.0_r8)
      VPK(NY,NX)=AMIN1(VPS(NY,NX),VPK(NY,NX))
!
      !     UA=wind speed
!
      WindSpeedAtm(NY,NX)=AMAX1(3600.0_r8,WIND(I))
!
!     TSNOW=temperature below which precipitation is snow (oC)
!     PrecAsRain,PrecAsSnow=rainfall,snowfall
!
      IF(J.GE.13.AND.J.LE.24)THEN
        IF(TCA(NY,NX).GT.TSNOW)THEN
          PrecAsRain(NY,NX)=RAIN(I)/12.0_r8   !rainfall
          PrecAsSnow(NY,NX)=0.0_r8   !snow fall
        ELSE
          PrecAsRain(NY,NX)=0.0_r8
          PrecAsSnow(NY,NX)=RAIN(I)/12.0_r8   !precip as snowfall
    !     IF(PrecAsSnow(NY,NX).LT.0.25E-03_r8)PrecAsSnow(NY,NX)=0.0_r8
        ENDIF
      ELSE
        PrecAsRain(NY,NX)=0.0_r8
        PrecAsSnow(NY,NX)=0.0_r8
      ENDIF
    ENDDO
  enddo
!
  end subroutine DailyWeather
!------------------------------------------------------------------------------------------

  subroutine HourlyWeather(I,J,NHW,NHE,NVN,NVS,RADN,PrecAsRain,PrecAsSnow,VPS)
!
!     Description:
!
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RADN(JY,JX),PrecAsRain(JY,JX),PrecAsSnow(JY,JX)
  real(r8), intent(out) :: VPS(JY,JX)
  integer :: NY,NX
  !     begin_execution

  DO  NX=NHW,NHE
    DO NY=NVN,NVS
 !
      !     RADN=SW radiation at horizontal surface,MJ/m2/hour
      !     TCA,TairK=air temperature (oC,K)
      !     VPK,VPS=ambient,saturated vapor pressure
      !     UA=wind speed
      !     TSNOW=temperature below which precipitation is snow (oC)
      !     PrecAsRain,PrecAsSnow=rainfall,snowfall
!
      RADN(NY,NX)=SWRad_hrly(J,I)  
      TCA(NY,NX)=TMPH(J,I)

      TairK(NY,NX)=units%Celcius2Kelvin(TCA(NY,NX))
      !elevation corrected saturated air vapor pressure
      VPS(NY,NX)=vapsat0(TairK(ny,nx))*EXP(-ALTI(NY,NX)/7272.0_r8)
      VPK(NY,NX)=AMIN1(DWPTH(J,I),VPS(NY,NX))
      WindSpeedAtm(NY,NX)=AMAX1(3600.0_r8,WINDH(J,I))

      !snowfall is determined by air tempeature
      IF(TCA(NY,NX).GT.TSNOW)THEN
        PrecAsRain(NY,NX)=RAINH(J,I)
        PrecAsSnow(NY,NX)=0.0_r8
      ELSE
        PrecAsRain(NY,NX)=0.0_r8
        PrecAsSnow(NY,NX)=RAINH(J,I)
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
  real(r8) :: AZI  !solar azimuth contribution
  REAL(R8) :: DEC  !solar declination contribution
  real(r8) :: DECLIN
  real(r8), PARAMETER :: PICON12=PICON/12._r8
  !     begin_execution
  !     CALCULATE DIRECT, DIFFUSE AND LONGWAVE RADIATION FROM
  !     INCOMING RADIATION READ IN 'READS', SOLAR ANGLE, HUMIDITY,
  !     TEMPERATURE AND CLOUDINESS
  !

  DECLIN=get_sun_declin(I)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     IF OUTDOORS
!
!     SineSolarIncliAngle,SineSolarIncliAngleNextHour=sine solar angle of current,next hour
!     RADX=solar constant at horizontal surface
!     RADN=SW radiation at horizontal surface, MJ/
!     IETYP: koppen climate zone
      IF(KoppenClimZone(NY,NX).GE.-1)THEN
        AZI=SIN(ALAT(NY,NX)*RadianPerDegree)*SIN(DECLIN*RadianPerDegree)
        DEC=COS(ALAT(NY,NX)*RadianPerDegree)*COS(DECLIN*RadianPerDegree)
        !check eq.(11.1) in Campbell and Norman, 1998, p168.
        SineSolarIncliAngle(NY,NX)=AZMAX1(AZI+DEC*COS(PICON12*(SolarNoonHour_col(NY,NX)-(J-0.5_r8))))
        SineSolarIncliAngleNextHour(NY,NX)=AZMAX1(AZI+DEC*COS(PICON12*(SolarNoonHour_col(NY,NX)-(J+0.5_r8))))

        !IF(SineSolarIncliAngle(NY,NX).GT.0.0_r8.AND.SineSolarIncliAngle(NY,NX).LT.TWILGT)SineSolarIncliAngle(NY,NX)=TWILGT
        IF(RADN(NY,NX).LE.0.0_r8)SineSolarIncliAngle(NY,NX)=0.0_r8
        IF(SineSolarIncliAngle(NY,NX).LE.-TWILGT)RADN(NY,NX)=0.0_r8
        RADX=SolConst*AZMAX1(SineSolarIncliAngle(NY,NX))
        RADN(NY,NX)=AMIN1(RADX,RADN(NY,NX))
!
        !     DIRECT VS DIFFUSE RADIATION IN SOLAR OR SKY BEAMS
        !
        !     RADZ=diffuse radiation at horizontal surface
        !     RADS,RADY,RAPS,PARDiffus_col=direct,diffuse SW,PAR in solar beam
!       !

        RADZ=AMIN1(RADN(NY,NX),0.5_r8*(RADX-RADN(NY,NX)))
        RadSWDirect_col(NY,NX)=safe_adb(RADN(NY,NX)-RADZ,SineSolarIncliAngle(NY,NX))
        RadSWDirect_col(NY,NX)=AMIN1(4.167_r8,RadSWDirect_col(NY,NX))
        RadSWDiffus_col(NY,NX)=RADZ/TotSineSkyAngles_grd
        PARDirect_col(NY,NX)=RadSWDirect_col(NY,NX)*CDIR*PDIR  !MJ/m2/hr
        PARDiffus_col(NY,NX)=RadSWDiffus_col(NY,NX)*CDIF*PDIF  !MJ/m2/hr
        !
        !     ATMOSPHERIC RADIATIVE PROPERTIES 
        !  Duarte et al., 2006, AFM 139:171
        !     CLD=cloudiness factor for EMM
        !     EMM=sky emissivity
        !     VPK,TairK=vapor pressure,temperature
        !
        IF(RADX.GT.ZERO)THEN
          CLD=AMIN1(1.0_r8,AMAX1(0.2_r8,2.33_r8-3.33_r8*RADN(NY,NX)/RADX))
        ELSE
          CLD=0.2_r8
        ENDIF
        EMM=0.625_r8*AMAX1(1.0_r8,(1.0E+03_r8*VPK(NY,NX)/TairK(NY,NX))**0.131_r8)
        EMM=EMM*(1.0_r8+0.242_r8*CLD**0.583_r8)
        !
        !     IF PHYTOTRON
        !
      ELSE
        IF(RADN(NY,NX).LE.0.0_r8)THEN
          SineSolarIncliAngle(NY,NX)=0.0_r8
        ELSE
          SineSolarIncliAngle(NY,NX)=1.0_r8
        ENDIF
        SineSolarIncliAngleNextHour(NY,NX)=1.0_r8
        CLD=0.0_r8
        EMM=0.96
      ENDIF
!
      !     LONGWAVE RADIATION
      !
      !     RadLWClm=longwave radiation
      !     THSX=longwave radiation from weather file or calculated from
      !     atmospheric properties
!
      IF(RadLWClm(J,I).GT.0.0_r8)THEN
        !     THSX(NY,NX)=EMM*(stefboltz_const*TairK(NY,NX)**4)
        !     THSX(NY,NX)=THSX(NY,NX)+RadLWClm(J,I)
        THSX(NY,NX)=RadLWClm(J,I)
      ELSE
        THSX(NY,NX)=EMM*(stefboltz_const*TairK(NY,NX)**4._r8)
      ENDIF
!
      !     INSERT CESM WEATHER HERE
      !
      !     ELSE
      !     RADS=DIRECT SW RADIATION (MJ M-2 H-1)
      !     RadSWDiffus_col=INDIRECT SW RADIATION (MJ M-2 H-1)
      !     PARDirect_col=DIRECT PAR (UMOL M-2 S-1)
      !     PARDiffus_col=INDIRECT PAR (UMOL M-2 S-1)
      !     THSX=LW RADIATION (MJ M-2 H-1)
      !     TCA=AIR TEMPERATURE (C)
      !     TairK=AIR TEMPERATURE (K)
      !     VPK=VAPOR PRESSURE (KPA)
      !     UA=WINDSPEED (M H-1)
      !     PrecAsRain(NY,NX)=RAIN (M H-1)
      !     PrecAsSnow(NY,NX)=SNOW (M H-1)
      !     SineSolarIncliAngle=SOLAR ANGLE CURRENT HOUR (SINE)
      !     SineSolarIncliAngleNextHour=SOLAR ANGLE NEXT HOUR (SINE)
      !     ENDIF
      !
      !     ADD IRRIGATION
      !
      !     PRECII,PRECUI=surface,subsurface irrigation
      !     RRIG=irrigation from soil management file in reads.f
      !
      WDPTHD=WDPTH(I,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
!     IF(WDPTHD.LE.CumDepth2LayerBottom(NU(NY,NX),NY,NX))THEN
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

  subroutine CorrectClimate(I,J,NHW,NHE,NVN,NVS,PRECUI,PrecAsRain,PRECII,PrecAsSnow,VPS)
  !
  !     DESCRIPTION:
  !     IMPLEMENT CLIMATE CHANGES READ IN 'READS' TO HOURLY TEMPERATURE,
  !     RADIATION, WINDSPEED,VAPOR PRESSURE, PRECIPITATION, IRRIGATION
  !     AND CO2
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: PRECUI(JY,JX),PrecAsRain(JY,JX),PRECII(JY,JX)
  real(r8), intent(inout) :: PrecAsSnow(JY,JX),VPS(JY,JX)
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
  D9925: DO NX=NHW,NHE
    D9920: DO NY=NVN,NVS
!
      !     TEMPERATURE CHANGE
      !
      !     TDTPX,TDTPN=change in max,min temperature
      !     DTA,AMP,DHR=change in daily average,amplitude of air temperature
      !     DHR=diurnal effect on AMP
      !
      IF(.not.isclose(TDTPX(N,NY,NX),0.0_r8).OR.(.not.isclose(TDTPN(N,NY,NX),0.0_r8)))THEN
        DTA=0.5_r8*(TDTPX(N,NY,NX)+TDTPN(N,NY,NX))
        AMP=0.5_r8*(TDTPX(N,NY,NX)-TDTPN(N,NY,NX))
        DHR=SIN(0.2618_r8*(J-(SolarNoonHour_col(NY,NX)+3.0_r8))+PICON2h)
        TCA(NY,NX)=TCA(NY,NX)+DTA+AMP*DHR
        TairK(NY,NX)=units%Celcius2Kelvin(TCA(NY,NX))
!
        !     ACCLIMATION TO GRADUAL CLIMATE CHANGE
        !
        !     DTS=change in daily average soil temperature
        !     ATCA,ATCS=mean annual air,soil temperature
        !     OFFSET=shift in Arrhenius curve for MFT activity in nitro.f
        !     OFFST=shift in Arrhenius curve for PFT activity in uptake.f
        !     iPlantThermoAdaptZone=PFT thermal adaptation zone
        !     HTC=high temperature threshold for grain number loss (oC)
        !     GROUPI,XTLI=node number at floral initiation,planting (maturity group)
!
        IF(ICLM.EQ.2.AND.J.EQ.1)THEN
          DTS=0.5_r8*DTA
          ATCA(NY,NX)=ATCAI(NY,NX)+DTA
          ATCS(NY,NX)=ATCAI(NY,NX)+DTS
          OFFSET(NY,NX)=0.33*(12.5-AZMAX1(AMIN1(25.0,ATCS(NY,NX))))
          DO NZ=1,NP(NY,NX)
            iPlantThermoAdaptZone(NZ,NY,NX)=iPlantInitThermoAdaptZone(NZ,NY,NX)+0.30/2.667*DTA
            OFFST(NZ,NY,NX)=2.667*(2.5-iPlantThermoAdaptZone(NZ,NY,NX))
            !     TCelsChill4Leaf_pft(NZ,NY,NX)=TCZD-OFFST(NZ,NY,NX)
            !     TCelcius4LeafOffHarden_pft(NZ,NY,NX)=AMIN1(15.0,TCelsChill4Leaf_pft(NZ,NY,NX)+TCXD)
            IF(iPlantPhotosynthesisType(NZ,NY,NX).EQ.3)THEN
              HighTCLimtSeed_pft(NZ,NY,NX)=27.0+3.0*iPlantThermoAdaptZone(NZ,NY,NX)
            ELSE
              HighTCLimtSeed_pft(NZ,NY,NX)=30.0+3.0*iPlantThermoAdaptZone(NZ,NY,NX)
            ENDIF
            MatureGroup_pft(NZ,NY,NX)=GROUPX(NZ,NY,NX)+0.30*DTA
            IF(iPlantTurnoverPattern_pft(NZ,NY,NX).NE.0)THEN
              MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)/25.0
            ENDIF
            MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)-XTLI(NZ,NY,NX)

          ENDDO
        ENDIF
!
        !     ADJUST VAPOR PRESSURE FOR TEMPERATURE CHANGE
!
        IF(isclose(DHUM(N),1.0_r8))THEN
          VPX=VPS(NY,NX)
          !VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TairK(NY,NX))) &
          vps(ny,ny)=vapsat0(TairK(ny,nx))*EXP(-ALTI(NY,NX)/7272.0)
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
      RadSWDirect_col(NY,NX)=RadSWDirect_col(NY,NX)*TDRAD(N,NY,NX)
      RadSWDiffus_col(NY,NX)=RadSWDiffus_col(NY,NX)*TDRAD(N,NY,NX)
      PARDirect_col(NY,NX)=PARDirect_col(NY,NX)*TDRAD(N,NY,NX)
      PARDiffus_col(NY,NX)=PARDiffus_col(NY,NX)*TDRAD(N,NY,NX)
      WindSpeedAtm(NY,NX)=WindSpeedAtm(NY,NX)*TDWND(N,NY,NX)
      VPK(NY,NX)=AMIN1(VPS(NY,NX),VPK(NY,NX)*TDHUM(N,NY,NX))
      PrecAsRain(NY,NX)=PrecAsRain(NY,NX)*TDPRC(N,NY,NX)
      PrecAsSnow(NY,NX)=PrecAsSnow(NY,NX)*TDPRC(N,NY,NX)
      PRECII(NY,NX)=PRECII(NY,NX)*TDIRI(N,NY,NX)
      PRECUI(NY,NX)=PRECUI(NY,NX)*TDIRI(N,NY,NX)
      CO2E(NY,NX)=CO2EI(NY,NX)   !used in photosynthesis, soil CO2 transport
      NH4_rain_conc(NY,NX)=CN4RI(NY,NX)*TDCN4(N,NY,NX)
      NO3_rain_conc(NY,NX)=CNORI(NY,NX)*TDCNO(N,NY,NX)
    ENDDO D9920
  ENDDO D9925
  end subroutine CorrectClimate
!------------------------------------------------------------------------------------------

  subroutine SummaryForOutput(NHW,NHE,NVN,NVS,PRECUI,PrecAsRain,PRECII,PrecAsSnow)
  !
  !     DESCRIPTION:
  !     DAILY WEATHER TOTALS, MAXIMA AND MINIMA FOR DAILY OUTPUT
  !     CHECK AGAINST INPUTS FROM WEATHER FILE
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(in):: PRECUI(JY,JX),PrecAsRain(JY,JX)
  real(r8), intent(in) :: PRECII(JY,JX),PrecAsSnow(JY,JX)
  integer :: NY,NX
  !
  !     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(SineSolarIncliAngle(NY,NX).GT.0.0_r8)TRAD(NY,NX)=TRAD(NY,NX)+RadSWDirect_col(NY,NX) &
        *SineSolarIncliAngle(NY,NX)+RadSWDiffus_col(NY,NX)*TotSineSkyAngles_grd
      TAMX(NY,NX)=AMAX1(TAMX(NY,NX),TCA(NY,NX))
      TAMN(NY,NX)=AMIN1(TAMN(NY,NX),TCA(NY,NX))
      HUDX(NY,NX)=AMAX1(HUDX(NY,NX),VPK(NY,NX))
      HUDN(NY,NX)=AMIN1(HUDN(NY,NX),VPK(NY,NX))
      TWIND(NY,NX)=TWIND(NY,NX)+WindSpeedAtm(NY,NX)
      VPA(NY,NX)=VPK(NY,NX)*2.173E-03_r8/TairK(NY,NX)
      TRAI(NY,NX)=TRAI(NY,NX)+(PrecAsRain(NY,NX)+PrecAsSnow(NY,NX) &
        +PRECII(NY,NX)+PRECUI(NY,NX))*1000.0_r8
!
      !     WATER AND HEAT INPUTS TO GRID CELLS
      !
      !     AREA=area of grid cell (m2)
      !
      !     PRECR,PRECW=rainfall,snowfall
      !     PRECI,PRECU=surface,subsurface irrigation
      !     PRECA,PrecAtm=rain+irrigation,rain+snow
      !     THS=sky LW radiation
!
      RainFalPrec(NY,NX)=PrecAsRain(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      SnoFalPrec(NY,NX)=PrecAsSnow(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      IrrigSurface(NY,NX)=PRECII(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      IrrigSubsurf(NY,NX)=PRECUI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PrecRainAndSurfirrig(NY,NX)=RainFalPrec(NY,NX)+IrrigSurface(NY,NX)
      PrecAtm(NY,NX)=RainFalPrec(NY,NX)+SnoFalPrec(NY,NX)
      LWRadSky(NY,NX)=THSX(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
    ENDDO
  ENDDO
  end subroutine SummaryForOutput

end module WthrMod
