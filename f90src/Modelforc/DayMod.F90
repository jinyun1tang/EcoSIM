
  module DayMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod  , only : isLeap,AZMAX1
  use MiniFuncMod  , only : GetDayLength
  use CanopyRadDataType
  use EcosimConst  
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
  use SurfLitterDataType, only : XTillCorp_col
  use PlantDataRateType
  use CanopyDataType
  use RootDataType
  use PlantMgmtDataType
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
      call UpdateDailyAccumulators(I, NHW, NHE, NVN, NVS)
      exit
    ENDIF
    NN=N
  ENDDO D500

  call TillageandIrrigationEvents(I, NHW, NHE, NVN, NVS)

  
  if(ldo_sp_mode)call PrescribePhenologyInterp(I, NHW, NHE, NVN, NVS)

  RETURN

  END subroutine day

!------------------------------------------------------------------------------------------     
  subroutine PrescribePhenologyInterp(I, NHW, NHE, NVN, NVS)

  implicit none
  integer, intent(in) :: I
  Integer, intent(in) ::  NHW, NHE, NVN, NVS
  real(r8) :: t
  integer :: it(2)
  integer :: months(2)
  integer :: kmo,dofmon,ndaysmon,NY,NX,NZ,L
  real(r8) :: timwt(2)
  real(r8) :: ZL1(0:NumCanopyLayers)
  real(r8) :: AreaInterval,AreaL
  real(r8) :: ARX  !interval canopy area: leaf+stem
  real(r8) :: DZL  !canopy interval height 
  real(r8) :: lai(12)=(/1.1852, 1.1821, 1.1554, 1.2433, 1.2922, 1.3341, 1.2296, 1.4118, 1.4343, 1.3941, 1.2721, 1.2218/)
  real(r8) :: sai(12)=(/0.3190, 0.3058, 0.3058, 0.3032, 0.3058, 0.3117, 0.3433, 0.3032, 0.3084, 0.3292, 0.3656, 0.3249/)

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        !FOR test only
        tlai_mon_pft(:,NZ,NY,NX)       = LAI
        tsai_mon_pft(:,NZ,NY,NX)       = SAI
        height_top_mon_pft(:,NZ,NY,NX) = 17._r8
        LeafAngleClass_pft(:,NZ,NY,NX) = 1./real(NumLeafZenithSectors,kind=r8)
      ENDDO
    ENDDO
  ENDDO    

  ndaysmon = etimer%get_curr_mon_days()
  dofmon   = etimer%get_curr_dom()
  kmo      = etimer%get_curr_mon()
  t = (dofmon-0.5_r8) / ndaysmon
  it(1) = t + 0.5_r8
  it(2) = it(1) + 1
  months(1) = kmo + it(1) - 1
  months(2) = kmo + it(2) - 1
  if (months(1) <  1) months(1) = 12
  if (months(2) > 12) months(2) = 1

  timwt(1) = (it(1)+0.5_r8) - t
  timwt(2) = 1._r8-timwt(1)
  
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      CanopyLeafArea_col(NY,NX) = 0._r8
      StemArea_col(NY,NX)       = 0._r8
      CanopyHeight_col(NY,NX)   = 0.0_r8
      NP_col(NY,NX)=1
      DO NZ=1,NP_col(NY,NX)
        tlai_day_pft(NZ,NY,NX)     = timwt(1)*tlai_mon_pft(months(1),NZ,NY,NX)+timwt(2)*tlai_mon_pft(months(2),NZ,NY,NX)
        tsai_day_pft(NZ,NY,NX)     = timwt(1)*tsai_mon_pft(months(1),NZ,NY,NX)+timwt(2)*tsai_mon_pft(months(2),NZ,NY,NX)
        CanopyHeight_pft(NZ,NY,NX) = timwt(1)*height_top_mon_pft(months(1),NZ,NY,NX)+timwt(2)*height_top_mon_pft(months(2),NZ,NY,NX)
        CanopyLeafArea_col(NY,NX)  = CanopyLeafArea_col(NY,NX)+tlai_day_pft(NZ,NY,NX)
        StemArea_col(NY,NX)        = StemArea_col(NY,NX)+tsai_day_pft(NZ,NY,NX)
        CanopyHeight_col(NY,NX)    = AMAX1(CanopyHeight_col(NY,NX),CanopyHeight_pft(NZ,NY,NX))
        CanopyLeafAreaZ_pft(1:NumCanopyLayers,NZ,NY,NX)=tlai_day_pft(NZ,NY,NX)/real(NumCanopyLayers,kind=r8)
        CanopyStemAreaZ_pft(1:NumCanopyLayers,NZ,NY,NX)=tsai_day_pft(NZ,NY,NX)/real(NumCanopyLayers,kind=r8)

      ENDDO

      !set vertical desitribution of LAI and             
      CanopyLeafAareZ_col(1:NumCanopyLayers,NY,NX)=CanopyLeafArea_col(NY,NX)/real(NumCanopyLayers,kind=r8)
      CanopyStemAareZ_col(1:NumCanopyLayers,NY,NX)=StemArea_col(NY,NX)/real(NumCanopyLayers,kind=r8)

      !divide canopy height      
      CanopyHeightZ_col(NumCanopyLayers,NY,NX) = CanopyHeight_col(NY,NX)+0.01_r8
      ZL1(NumCanopyLayers)               = CanopyHeightZ_col(NumCanopyLayers,NY,NX)
      ZL1(0)                                = 0.0_r8
      
      !for simplicity, right now unifrom division is used
      !divide total are into NumCanopyLayers1, from top to bottom
      AreaInterval=(CanopyLeafArea_col(NY,NX)+StemArea_col(NY,NX))/NumCanopyLayers  
      DZL=CanopyHeightZ_col(NumCanopyLayers,NY,NX)/real(NumCanopyLayers,kind=r8)

      IF(AreaInterval.GT.ZEROS(NY,NX))THEN
         DO L=NumCanopyLayers,1,-1
           CanopyHeightZ_col(L-1,NY,NX)=CanopyHeightZ_col(L,NY,NX)-DZL
         ENDDO
!        DO L=NumCanopyLayers,2,-1
!          AreaL=CanopyLeafAareZ_col(L,NY,NX)+CanopyStemAareZ_col(L,NY,NX)
!          print*,'AreaInterval',AreaInterval,AreaL
!          !greater than the mean leaf area or area-interval
!          IF(AreaL.GT.1.01_r8*AreaInterval)THEN
!            DZL      = CanopyHeightZ_col(L,NY,NX)-CanopyHeightZ_col(L-1,NY,NX)
!            ZL1(L-1) = CanopyHeightZ_col(L-1,NY,NX)+0.5_r8*AMIN1(1.0_r8,(AreaL-AreaInterval)/AreaL)*DZL
!          ELSEIF(AreaL.LT.0.99_r8*AreaInterval)THEN
!            ARX = CanopyLeafAareZ_col(L-1,NY,NX)+CanopyStemAareZ_col(L-1,NY,NX)
!            DZL = CanopyHeightZ_col(L-1,NY,NX)-CanopyHeightZ_col(L-2,NY,NX)
            
!            !layer L-1 has significant leaf+stem (canopy) area
!            IF(ARX.GT.ZEROS(NY,NX))THEN
!              ZL1(L-1)=CanopyHeightZ_col(L-1,NY,NX)-0.5_r8*AMIN1(1.0_r8,(AreaInterval-AreaL)/ARX)*DZL
!            ELSE
!              ZL1(L-1)=CanopyHeightZ_col(L-1,NY,NX)
!            ENDIF
!          ELSE
!            ZL1(L-1)=CanopyHeightZ_col(L-1,NY,NX)
!          ENDIF
!        ENDDO

!        DO L=NumCanopyLayers,2,-1
!          CanopyHeightZ_col(L-1,NY,NX)=ZL1(L-1)
!        ENDDO
      ENDIF

    ENDDO
  ENDDO  

  end subroutine PrescribePhenologyInterp
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
!     RESET ANNUAL FLUX ACCUMULATORS AT START OF ANNUAL CYCLE
!     ALAT=latitude +ve=N,-ve=S
!
      HUDX_col(NY,NX)=0._r8
      HUDN_col(NY,NX)=100.0_r8
      TWIND_col(NY,NX)=0._r8
!      PrecDaily_col(NY,NX)=0._r8
!
!
!     CALCULATE DAYLENGTH FROM SOLAR ANGLE
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     ALAT=latitude
!
      DayLenthPrev_col(NY,NX)=DayLensCurr_col(NY,NX)

      DayLensCurr_col(NY,NX)=GetDayLength(ALAT_col(NY,NX),I)
!
!     TIME STEP OF WEARHER DATA
!     ITYPE 1=daily,2=hourly
!
      ITYPE=IWTHR


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
        IF(KoppenClimZone_col(NY,NX).GE.-1)THEN
          IF(DayLensCurr_col(NY,NX).GT.ZERO)THEN
            RMAX=SRAD(I)/(DayLensCurr_col(NY,NX)*0.658_r8)
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
        IF(I.EQ.DazCurrYear)I3=I
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
!     TDCN4,TDCNO=change in atm CO2,NH4,NO3 concn in precipitation
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
          TDCN4(N,NY,NX)=DCN4R(N)
          TDCNO(N,NY,NX)=DCNOR(N)
!
!     INCRENENTAL CHANGES
!
        ELSEIF(ICLM.EQ.2)THEN
! DazCurrYear: number of days in current year
          TDTPX(N,NY,NX)=TDTPX(N,NY,NX)+DTMPX(N)/DazCurrYear
          TDTPN(N,NY,NX)=TDTPN(N,NY,NX)+DTMPN(N)/DazCurrYear
          TDRAD(N,NY,NX)=TDRAD(N,NY,NX)+(DRAD(N)-1.0_r8)/DazCurrYear
          TDWND(N,NY,NX)=TDWND(N,NY,NX)+(DWIND(N)-1.0_r8)/DazCurrYear
          TDHUM(N,NY,NX)=TDHUM(N,NY,NX)+(DHUM(N)-1.0_r8)/DazCurrYear
          TDPRC(N,NY,NX)=TDPRC(N,NY,NX)+(DPREC(N)-1.0_r8)/DazCurrYear
          TDIRI(N,NY,NX)=TDIRI(N,NY,NX)+(DIRRI(N)-1.0_r8)/DazCurrYear
          TDCN4(N,NY,NX)=TDCN4(N,NY,NX)+(DCN4R(N)-1.0_r8)/DazCurrYear
          TDCNO(N,NY,NX)=TDCNO(N,NY,NX)+(DCNOR(N)-1.0_r8)/DazCurrYear
        ENDIF
      ENDDO D600
    ENDDO D950
  ENDDO D955

  END subroutine UpdateDailyAccumulators
!-----------------------------------------------------------------------------------------

  subroutine TillageandIrrigationEvents(I, NHW, NHE, NVN, NVS)
  !
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
      IF(iSoilDisturbType_col(I,NY,NX).LE.10)THEN
        ! type-1 tillage
        CORP=AMIN1(1.0_r8,AZMAX1(iSoilDisturbType_col(I,NY,NX)/10.0_r8))
      ELSEIF(iSoilDisturbType_col(I,NY,NX).LE.20)THEN
        ! type-2 tillage
        CORP=AMIN1(1.0_r8,AZMAX1((iSoilDisturbType_col(I,NY,NX)-10.0_r8)/10.0_r8))
      ENDIF
      !fraction of material not mixed
      XTillCorp_col(NY,NX)=1.0_r8-CORP
!     WRITE(*,2227)'TILL',I,iSoilDisturbType_col(I,NY,NX),CORP,XTillCorp_col(NY,NX)
!2227  FORMAT(A8,2I4,12E12.4)
!
!     AUTOMATIC IRRIGATION IF SELECTED
!
!     DATA1(6)=irrigation file name
!     IIRRA=start,finish dates(1,2),hours(3,4) of automated irrigation
!     DIRRX=depth to which water depletion and rewatering is calculated(1)
!     DIRRA=depth to,at which irrigation is applied(1,2)
!     POROS,FC,WP=water content at saturation,field capacity,wilting point
!     CIRRA_col= fraction of FC to which irrigation will raise SWC
!     FW=fraction of soil layer in irrigation zone
!     FZ=SWC at which irrigation is triggered
!     VLSoilPoreMicP_vr,VOLW,VOLI=total,water,ice volume
!     IFLGV_col=flag for irrigation criterion,0=SWC,1=canopy water potential
!     FIRRA_col=depletion of SWC from CIRRA_col to WP(IFLGV_col=0),or minimum canopy
!     water potential(IFLGV_col=1), to trigger irrigation
!     RR=total irrigation requirement
!     RRIG=hourly irrigation amount applied in wthr.f
!
      IF(Lirri_auto)THEN
      !automated irrigation
        IF(I.GE.IIRRA(1,NY,NX).AND.I.LE.IIRRA(2,NY,NX))THEN
          TFZ=0._r8
          TWP=0._r8
          TVW=0._r8
          DIRRA1=DIRRA(1,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
          DIRRA2=DIRRA(2,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)

          D165: DO L=NU_col(NY,NX),NL_col(NY,NX)
            IF(CumDepz2LayBottom_vr(L-1,NY,NX).LT.DIRRA1)THEN
              FW=AMIN1(1.0_r8,(DIRRA1-CumDepz2LayBottom_vr(L-1,NY,NX)) &
                /(CumDepz2LayBottom_vr(L,NY,NX)-CumDepz2LayBottom_vr(L-1,NY,NX)))
              FZ=AMIN1(POROS_vr(L,NY,NX),WiltPoint_vr(L,NY,NX)+CIRRA_col(NY,NX)*(FieldCapacity_vr(L,NY,NX)-WiltPoint_vr(L,NY,NX)))
              TFZ=TFZ+FW*FZ*VLSoilPoreMicP_vr(L,NY,NX)
              TWP=TWP+FW*WiltPoint_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
              TVW=TVW+FW*(VLWatMicP_vr(L,NY,NX)+VLiceMicP_vr(L,NY,NX))
            ENDIF
          ENDDO D165

          IF((IFLGV_col(NY,NX).EQ.0 .AND. TVW.LT.TWP+FIRRA_col(NY,NX)*(TFZ-TWP)) &
            .OR.(IFLGV_col(NY,NX).EQ.1.AND.PSICanPDailyMin_pft(1,NY,NX).LT.FIRRA_col(NY,NX)))THEN
            RR=AZMAX1(TFZ-TVW)
            IF(RR.GT.0.0_r8)THEN
              D170: DO J=IIRRA(3,NY,NX),IIRRA(4,NY,NX)
                RRIG(J,I,NY,NX)=RR/(IIRRA(4,NY,NX)-IIRRA(3,NY,NX)+1)
              ENDDO D170
              WDPTH(I,NY,NX)=DIRRA(2,NY,NX)
              WRITE(*,2222)'auto',iYearCurrent,I,IIRRA(3,NY,NX),IIRRA(4,NY,NX) &
                ,IFLGV_col(NY,NX),RR,TFZ,TVW,TWP,FIRRA_col(NY,NX),PSICanPDailyMin_pft(1,NY,NX) &
                ,CIRRA_col(NY,NX),DIRRA1,WDPTH(I,NY,NX)

2222  FORMAT(A8,5I6,40E12.4)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D9990
  ENDDO D9995
  end subroutine TillageandIrrigationEvents
END module DayMod
