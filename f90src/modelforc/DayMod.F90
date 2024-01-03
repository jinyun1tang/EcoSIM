
  module DayMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use minimathmod  , only : isLeap,AZMAX1
  use MiniFuncMod  , only : GetDayLength
  use EcoSIMCtrlMod
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use PlantTraitDataType
  use SurfLitterDataType, only : XCORP
  use PlantDataRateType
  use CanopyDataType
  use RootDataType
  use PlantMngmtDataType
  use SOMDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use GridDataType
  use EcoSIMConfig
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  CHARACTER(len=3) :: CHARN1,CHARN2
  CHARACTER(len=4) :: CHARN3

  real(r8) :: CORP,DIRRA1,DIRRA2,FW,FZ
  real(r8) :: RR,TFZ,TWP,TVW,XI

  integer :: ITYPE,I2,I3,J,L,M,N,NN,N1,N2,N3,NX,NY,NZ

  public :: day
  contains

  SUBROUTINE day(I,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE REINITIALIZES DAILY VARIABLES USED IN OTHER
!     SUBROUTINES E.G. LAND MANAGEMENT
!
  implicit none

  integer, intent(in) :: I
  integer, intent(in) :: NHW,NHE,NVN,NVS
!     execution begins here
!     begin_execution
!
!     WRITE DATE
!
!     CDATE=DDMMYYYY
!
  N=0
  NN=0
  D500: DO M=1,12
    N=30*M+ICOR(M)

!  leap year February.
    IF(isLeap(iYearCurrent) .and. M.GE.2)N=N+1
      IF(I.LE.N)THEN
        N1=I-NN
        N2=M
        N3=iYearCurrent
        WRITE(CHARN1,'(I3)')N1+100
        WRITE(CHARN2,'(I3)')N2+100
        WRITE(CHARN3,'(I4)')N3
        WRITE(CDATE,'(2A2,A4)')CHARN1(2:3),CHARN2(2:3),CHARN3(1:4)
        call UpdateDailyAccumulators(I, NHW, NHE, NVN, NVS)
        exit
      ENDIF
      NN=N
  ENDDO D500

  call TillageandIrrigationEvents(I, NHW, NHE, NVN, NVS)
  RETURN

  END subroutine day
!------------------------------------------------------------------------------------------

  subroutine SetAnnualAccumlators(NY,NX)
  
  !set accumulators at the beginning of year
  implicit none
  integer, intent(in) :: NY,NX

  integer :: NZ,NE

  AmendCFlx_col(NY,NX)=0._r8
  LiterfalOrgC_col(NY,NX)=0._r8
  UCOP(NY,NX)=0._r8
  Eco_NBP_col(NY,NX)=0._r8
  URAIN(NY,NX)=0._r8
  UEVAP(NY,NX)=0._r8
  URUN(NY,NX)=0._r8
  USEDOU(NY,NX)=0._r8
  UVOLO(NY,NX)=0._r8
  HydroIonFlx_col(NY,NX)=0._r8
  FertNFlx_col(NY,NX)=0._r8
  LiterfalOrgN_col(NY,NX)=0._r8
  HydroSufDONFlx_col(NY,NX)=0._r8
  HydroSufDINFlx_col(NY,NX)=0._r8
  FerPFlx_col(NY,NX)=0._r8
  LiterfalOrgP_col(NY,NX)=0._r8
  HydroSufDOPFlx_col(NY,NX)=0._r8
  HydroSufDIPFlx_col(NY,NX)=0._r8
  CO2byFire_col(NY,NX)=0._r8
  CH4byFire_col(NY,NX)=0._r8
  O2byFire_col(NY,NX)=0._r8
  N2ObyFire_col(NY,NX)=0._r8
  NH3byFire_col(NY,NX)=0._r8
  PO4byFire_col(NY,NX)=0._r8
  Eco_HR_col(NY,NX)=0._r8
  Eco_GPP_col(NY,NX)=0._r8
  Eco_NPP_col(NY,NX)=0._r8
  Eco_AutoR_col(NY,NX)=0._r8
  EcoHavstElmnt_col(:,NY,NX)=0._r8
  NetNH4Mineralize_col(NY,NX)=0._r8
  NetPO4Mineralize_col(NY,NX)=0._r8

  D960: DO NZ=1,NP0(NY,NX)
  !NetCumElmntFlx2Plant_pft: effect of canopy element status on seed set
    DO NE=1,NumPlantChemElmnts
      NetCumElmntFlx2Plant_pft(NE,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(NE,NZ,NY,NX)+PlantExudChemElmntCum_pft(NE,NZ,NY,NX)-LitrfallChemElmnts_pft(NE,NZ,NY,NX)
      EcoHavstElmntCum_pft(NE,NZ,NY,NX)=EcoHavstElmntCum_pft(NE,NZ,NY,NX)+EcoHavstElmnt_pft(NE,NZ,NY,NX)
    enddo
    NetCumElmntFlx2Plant_pft(ielmc,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmc,NZ,NY,NX)+GrossCO2Fix_pft(NZ,NY,NX)+GrossResp_pft(NZ,NY,NX) &
      -CO2ByFire_pft(NZ,NY,NX)-CH4ByFire_pft(NZ,NY,NX)
    NetCumElmntFlx2Plant_pft(ielmn,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmn,NZ,NY,NX)+NH3EmiCum_pft(NZ,NY,NX)+PlantN2FixCum_pft(NZ,NY,NX) &
      -NH3byFire_pft(NZ,NY,NX)-N2ObyFire_pft(NZ,NY,NX)
    NetCumElmntFlx2Plant_pft(ielmp,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmp,NZ,NY,NX)-PO4byFire_pft(NZ,NY,NX)

! the following variables are accumulated daily
    GrossCO2Fix_pft(NZ,NY,NX)=0._r8
    PlantExudChemElmntCum_pft(:,NZ,NY,NX)=0._r8
    GrossResp_pft(NZ,NY,NX)=0._r8
    CanopyPlusNoduRespC_pft(NZ,NY,NX)=0._r8
    ETCanopy_pft(NZ,NY,NX)=0._r8
    PlantN2FixCum_pft(NZ,NY,NX)=0._r8
    NH3EmiCum_pft(NZ,NY,NX)=0._r8
    CO2ByFire_pft(NZ,NY,NX)=0._r8
    CH4ByFire_pft(NZ,NY,NX)=0._r8
    VOXYF(NZ,NY,NX)=0._r8
    NH3byFire_pft(NZ,NY,NX)=0._r8
    N2ObyFire_pft(NZ,NY,NX)=0._r8
    PO4byFire_pft(NZ,NY,NX)=0._r8

    EcoHavstElmnt_pft(:,NZ,NY,NX)=0._r8
    SurfLitrfallChemElmnts_pft(:,NZ,NY,NX)=0._r8
    LitrfallChemElmnts_pft(:,NZ,NY,NX)=0._r8
  ENDDO D960
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    TSED(NY,NX)=0._r8
  ENDIF
  end subroutine SetAnnualAccumlators        
!------------------------------------------------------------------------------------------     
  subroutine UpdateDailyAccumulators(I, NHW, NHE, NVN, NVS)
!     WRITE DAILY MAX MIN ACCUMULATORS FOR WEATHER VARIABLES

  implicit none
  integer, intent(in) :: I, NHW, NHE, NVN, NVS

!  real(r8) :: AZI
!  REAL(R8) :: DEC
  integer :: NE

  D955: DO NX=NHW,NHE
    D950: DO NY=NVN,NVS
      TRAD(NY,NX)=0._r8
      TAMX(NY,NX)=-100.0_r8
      TAMN(NY,NX)=100.0_r8
      HUDX(NY,NX)=0._r8
      HUDN(NY,NX)=100.0_r8
      TWIND(NY,NX)=0._r8
      TRAI(NY,NX)=0._r8
!
!     RESET ANNUAL FLUX ACCUMULATORS AT START OF ANNUAL CYCLE
!     ALAT=latitude +ve=N,-ve=S
!
      IF((ALAT(NY,NX).GE.0.0.AND.I.EQ.1).OR.(ALAT(NY,NX).LT.0.0.AND.I.EQ.1))THEN        
        call SetAnnualAccumlators(NY,NX)
      ENDIF
!
!     CALCULATE DAYLENGTH FROM SOLAR ANGLE
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     ALAT=latitude
!
      DayLenthPrev(NY,NX)=DayLenthCurrent(NY,NX)

      DayLenthCurrent(NY,NX)=GetDayLength(ALAT(NY,NX),I)
!
!     TIME STEP OF WEARHER DATA
!     ITYPE 1=daily,2=hourly
!
      IF(is_first_year)THEN
      !weather input is in daily format
        ITYPE=IWTHR(1)
      ELSE
      !weather input is in hourly format
        ITYPE=IWTHR(2)
      ENDIF

!
!     PARAMETERS FOR CALCULATING HOURLY RADIATION, TEMPERATURE
!     AND VAPOR PRESSURE FROM DAILY VALUES
!
!     DLYN=daylength
!     SRAD=daily radiation from weather file
!     RMAX=maximum hourly radiation
!     I2,I,I3=previous,current,next day
!     TMPX,TMPN=maximum,minimum daily temperature from weather file
!     DWPT=daily vapor pressure from weather file
!     TAVG*,AMP*,VAVG*,VMP*=daily avgs, amps to calc hourly values in wthr.f
!     KoppenClimZone=Koppen climate zone

      IF(ITYPE.EQ.1)THEN
        IF(KoppenClimZone(NY,NX).GE.-1)THEN
          IF(DayLenthCurrent(NY,NX).GT.ZERO)THEN
            RMAX=SRAD(I)/(DayLenthCurrent(NY,NX)*0.658_r8)
          ELSE
            RMAX=0._r8
          ENDIF
        ELSE
          RMAX=SRAD(I)
        ENDIF
        I2=I-1
        I3=I+1
        IF(I.EQ.1)I2=LYRX
        IF(I.EQ.IBEGIN)I2=I
        IF(I.EQ.LYRC)I3=I
        TAVG1=(TMPX(I2)+TMPN(I))/2._r8
        TAVG2=(TMPX(I)+TMPN(I))/2._r8
        TAVG3=(TMPX(I)+TMPN(I3))/2._r8
        AMP1=TAVG1-TMPN(I)
        AMP2=TAVG2-TMPN(I)
        AMP3=TAVG3-TMPN(I3)
        VAVG1=(DWPT(1,I2)+DWPT(2,I))/2._r8
        VAVG2=(DWPT(1,I)+DWPT(2,I))/2._r8
        VAVG3=(DWPT(1,I)+DWPT(2,I3))/2._r8
        VMP1=VAVG1-DWPT(2,I)
        VMP2=VAVG2-DWPT(2,I)
        VMP3=VAVG3-DWPT(2,I3)
      ENDIF
!
!     MODIFIERS TO TEMPERATURE, RADIATION, WIND, HUMIDITY, PRECIPITATION,
!     IRRIGATION AND CO2 INPUTS FROM CLIMATE CHANGES ENTERED IN OPTION
!     FILE IN 'READS'
!
!     ICLM=type of climate change 1=step,2=incremental
!     TDTPX,TDTPN=change in max,min temperature
!     TDRAD,TDWND,TDHUM=change in radiation,windspeed,vapor pressure
!     TDPRC,TDIRRI=change in precipitation,irrigation
!     TDCO2,TDCN4,TDCNO=change in atm CO2,NH4,NO3 concn in precipitation
!
      D600: DO N=1,12
!
!     STEP CHANGES
!
        IF(ICLM.EQ.1)THEN
          TDTPX(N,NY,NX)=DTMPX(N)
          TDTPN(N,NY,NX)=DTMPN(N)
          TDRAD(N,NY,NX)=DRAD(N)
          TDWND(N,NY,NX)=DWIND(N)
          TDHUM(N,NY,NX)=DHUM(N)
          TDPRC(N,NY,NX)=DPREC(N)
          TDIRI(N,NY,NX)=DIRRI(N)
          TDCO2(N,NY,NX)=DCO2E(N)
          TDCN4(N,NY,NX)=DCN4R(N)
          TDCNO(N,NY,NX)=DCNOR(N)
!
!     INCRENENTAL CHANGES
!
        ELSEIF(ICLM.EQ.2)THEN
! LYRC: number of days in current year
          TDTPX(N,NY,NX)=TDTPX(N,NY,NX)+DTMPX(N)/LYRC
          TDTPN(N,NY,NX)=TDTPN(N,NY,NX)+DTMPN(N)/LYRC
          TDRAD(N,NY,NX)=TDRAD(N,NY,NX)+(DRAD(N)-1.0_r8)/LYRC
          TDWND(N,NY,NX)=TDWND(N,NY,NX)+(DWIND(N)-1.0_r8)/LYRC
          TDHUM(N,NY,NX)=TDHUM(N,NY,NX)+(DHUM(N)-1.0_r8)/LYRC
          TDPRC(N,NY,NX)=TDPRC(N,NY,NX)+(DPREC(N)-1.0_r8)/LYRC
          TDIRI(N,NY,NX)=TDIRI(N,NY,NX)+(DIRRI(N)-1.0_r8)/LYRC
          TDCO2(N,NY,NX)=TDCO2(N,NY,NX)*EXP(LOG(DCO2E(N))/LYRC)
          TDCN4(N,NY,NX)=TDCN4(N,NY,NX)+(DCN4R(N)-1.0_r8)/LYRC
          TDCNO(N,NY,NX)=TDCNO(N,NY,NX)+(DCNOR(N)-1.0_r8)/LYRC
        ENDIF
      ENDDO D600
    ENDDO D950
  ENDDO D955

  END subroutine UpdateDailyAccumulators
!-----------------------------------------------------------------------------------------

  subroutine TillageandIrrigationEvents(I, NHW, NHE, NVN, NVS)
!
  use EcoSIMCtrlMod, only : Lirri_auto
  implicit none

  integer, intent(in) :: I, NHW, NHE, NVN, NVS

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!     ATTRIBUTE MIXING COEFFICIENTS AND SURFACE ROUGHNESS PARAMETERS
!     TO TILLAGE EVENTS FROM SOIL MANAGEMENT FILE IN 'READS'
!
!     ITILL=soil disturbance type, 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     CORP=soil mixing fraction used in redist.f
!
      IF(ITILL(I,NY,NX).LE.10)THEN
      ! type-1 tillage
        CORP=AMIN1(1.0_r8,AZMAX1(ITILL(I,NY,NX)/10.0_r8))
      ELSEIF(ITILL(I,NY,NX).LE.20)THEN
      ! type-2 tillage
        CORP=AMIN1(1.0_r8,AZMAX1((ITILL(I,NY,NX)-10.0_r8)/10.0_r8))
      ENDIF

      XCORP(NY,NX)=1.0_r8-CORP
!     WRITE(*,2227)'TILL',I,ITILL(I,NY,NX),CORP,XCORP(NY,NX)
!2227  FORMAT(A8,2I4,12E12.4)
!
!     AUTOMATIC IRRIGATION IF SELECTED
!
!     DATA1(6)=irrigation file name
!     IIRRA=start,finish dates(1,2),hours(3,4) of automated irrigation
!     DIRRX=depth to which water depletion and rewatering is calculated(1)
!     DIRRA=depth to,at which irrigation is applied(1,2)
!     POROS,FC,WP=water content at saturation,field capacity,wilting point
!     CIRRA= fraction of FC to which irrigation will raise SWC
!     FW=fraction of soil layer in irrigation zone
!     FZ=SWC at which irrigation is triggered
!     VLSoilPoreMicP,VOLW,VOLI=total,water,ice volume
!     IFLGV=flag for irrigation criterion,0=SWC,1=canopy water potential
!     FIRRA=depletion of SWC from CIRRA to WP(IFLGV=0),or minimum canopy
!     water potential(IFLGV=1), to trigger irrigation
!     RR=total irrigation requirement
!     RRIG=hourly irrigation amount applied in wthr.f
!
      IF(Lirri_auto)THEN
      !automated irrigation
        IF(I.GE.IIRRA(1,NY,NX).AND.I.LE.IIRRA(2,NY,NX))THEN
          TFZ=0._r8
          TWP=0._r8
          TVW=0._r8
          DIRRA1=DIRRA(1,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
          DIRRA2=DIRRA(2,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)

          D165: DO L=NU(NY,NX),NL(NY,NX)
            IF(CumDepth2LayerBottom(L-1,NY,NX).LT.DIRRA1)THEN
              FW=AMIN1(1.0_r8,(DIRRA1-CumDepth2LayerBottom(L-1,NY,NX))/(CumDepth2LayerBottom(L,NY,NX)-CumDepth2LayerBottom(L-1,NY,NX)))
              FZ=AMIN1(POROS(L,NY,NX),WiltPoint(L,NY,NX)+CIRRA(NY,NX)*(FieldCapacity(L,NY,NX)-WiltPoint(L,NY,NX)))
              TFZ=TFZ+FW*FZ*VLSoilPoreMicP(L,NY,NX)
              TWP=TWP+FW*WiltPoint(L,NY,NX)*VLSoilPoreMicP(L,NY,NX)
              TVW=TVW+FW*(VLWatMicP(L,NY,NX)+VLiceMicP(L,NY,NX))
            ENDIF
          ENDDO D165

          IF((IFLGV(NY,NX).EQ.0 .AND. TVW.LT.TWP+FIRRA(NY,NX)*(TFZ-TWP)) &
            .OR.(IFLGV(NY,NX).EQ.1.AND.PSICanPDailyMin(1,NY,NX).LT.FIRRA(NY,NX)))THEN
            RR=AZMAX1(TFZ-TVW)
            IF(RR.GT.0.0_r8)THEN
              D170: DO J=IIRRA(3,NY,NX),IIRRA(4,NY,NX)
                RRIG(J,I,NY,NX)=RR/(IIRRA(4,NY,NX)-IIRRA(3,NY,NX)+1)
              ENDDO D170
              WDPTH(I,NY,NX)=DIRRA(2,NY,NX)
              WRITE(*,2222)'auto',iYearCurrent,I,IIRRA(3,NY,NX),IIRRA(4,NY,NX) &
                ,IFLGV(NY,NX),RR,TFZ,TVW,TWP,FIRRA(NY,NX),PSICanPDailyMin(1,NY,NX) &
                ,CIRRA(NY,NX),DIRRA1,WDPTH(I,NY,NX)

2222  FORMAT(A8,5I6,40E12.4)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D9990
  ENDDO D9995
  end subroutine TillageandIrrigationEvents
END module DayMod
