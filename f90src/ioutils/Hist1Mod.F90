module Hist1Mod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use fileUtil, only : open_safe
  use ElmIDMod
  use GridConsts
  use GridDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use FlagDataType
  use PlantTraitDataType
  use PlantDataRateType
  use ClimForcDataType
  use CanopyDataType
  use RootDataType
  use SOMDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use SnowDataType
  use ChemTranspDataType
  use PlantMngmtDataType
  use EcosimBGCFluxType
  use SoilPropertyDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SurfSoilDataType
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: fouts,foutp, outpd, outph, outsd, outsh
  contains

  SUBROUTINE foutp(NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE OPENS AND LABELS OUTPUT FILES
!     FOR PLANT DATA
! write file head for plant variable output
  implicit none
  integer, intent(in) :: NE,NEX
  integer, intent(in) :: NHW,NHE,NVN,NVS


  integer, SAVE :: LUN

  CHARACTER(len=8):: CDOY,DATE,HOUR
  CHARACTER(len=1):: CHARM,CHARN,CHARZ,CHARO
  CHARACTER(len=2):: CHARX,CHARY
  CHARACTER(len=4):: CHARR
  CHARACTER(len=16):: HEAD(50)

  integer :: L,M,NX,NY,NZ,N,K
!     execution begins here

  CDOY='DOY     '
  DATE='DATE    '
  HOUR='HOUR    '
  CHARM='0'
  CHARO='1'
  DO NX=NHW,NHE
    DO   NY=NVN,NVS
      DO   NZ=1,NP0(NY,NX)
        DO   N=1,10
          OUTFILP(N,NZ,NY,NX)= '                '
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!
!     OPEN AND NAME OUTPUT FILES
!
  D1010: DO N=21,30
    IF(DATAC(N,NE,NEX).NE.'NO')THEN
      WRITE(CHARR,'(I4)')IYRC
      OUTP(N-20)=CHARO//CHARR//DATAC(N,NE,NEX)
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP0(NY,NX)
            IF(NX.LT.10)THEN
              WRITE(CHARN,'(I1)')NX
              CHARX=CHARM//CHARN
            ELSE
              WRITE(CHARX,'(I2)')NX
            ENDIF
            IF(NY.LT.10)THEN
              WRITE(CHARN,'(I1)')NY
              CHARY=CHARM//CHARN
            ELSE
              WRITE(CHARY,'(I2)')NY
            ENDIF
            WRITE(CHARZ,'(I1)')NZ
            OUTFILP(N-20,NZ,NY,NX)=CHARX//CHARY//CHARZ//CHARR//DATAC(N,NE,NEX)
          ENDDO
        ENDDO
      ENDDO
      LUN=N+20
      OPEN(LUN,FILE=trim(outdir)//OUTP(N-20),STATUS='UNKNOWN')
!
!     WRITE HEADINGS TO OUTPUT FILES
!
! hourly plant carbon
      M=0
      IF(N.EQ.21)THEN
        DO  L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='CAN_CO2_FLUX'
            IF(L.EQ.52)HEAD(M)='CAN_GPP'
            IF(L.EQ.53)HEAD(M)='CAN_RA'
            IF(L.EQ.54)HEAD(M)='[TNC]'
            IF(L.EQ.55)HEAD(M)='STOML_RSC'
            IF(L.EQ.56)HEAD(M)='BLYR_RSC'
            IF(L.EQ.57)HEAD(M)='[CAN_CO2]'
            IF(L.EQ.58)HEAD(M)='LAI'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
! hourly plant water
      IF(N.EQ.22)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='PSI_CAN'
            IF(L.EQ.52)HEAD(M)='TURG_CAN'
            IF(L.EQ.53)HEAD(M)='STOM_RSC'
            IF(L.EQ.54)HEAD(M)='BLYR_RSC'
            IF(L.EQ.55)HEAD(M)='TRANSPN'
            IF(L.EQ.56)HEAD(M)='O2_STRESS'
            IF(L.EQ.57)HEAD(M)='PSI_RT_1'
            IF(L.EQ.58)HEAD(M)='PSI_RT_2'
            IF(L.EQ.59)HEAD(M)='PSI_RT_3'
            IF(L.EQ.60)HEAD(M)='PSI_RT_4'
            IF(L.EQ.61)HEAD(M)='PSI_RT_5'
            IF(L.EQ.62)HEAD(M)='PSI_RT_6'
            IF(L.EQ.63)HEAD(M)='PSI_RT_7'
            IF(L.EQ.64)HEAD(M)='PSI_RT_8'
            IF(L.EQ.65)HEAD(M)='PSI_RT_9'
            IF(L.EQ.66)HEAD(M)='PSI_RT_10'
            IF(L.EQ.67)HEAD(M)='PSI_RT_11'
            IF(L.EQ.68)HEAD(M)='PSI_RT_12'
            IF(L.EQ.69)HEAD(M)='PSI_RT_13'
            IF(L.EQ.70)HEAD(M)='PSI_RT_14'
            IF(L.EQ.71)HEAD(M)='PSI_RT_15'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!hourly plant nitrogen
      IF(N.EQ.23)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='NH4_UPTK'
            IF(L.EQ.52)HEAD(M)='NO3_UPTK'
            IF(L.EQ.53)HEAD(M)='N2_FIXN'
            IF(L.EQ.54)HEAD(M)='[TNN]'
            IF(L.EQ.55)HEAD(M)='NH3_FLUX'
            IF(L.EQ.56)HEAD(M)='UP_NH4_1'
            IF(L.EQ.57)HEAD(M)='UP_NH4_2'
            IF(L.EQ.58)HEAD(M)='UP_NH4_3'
            IF(L.EQ.59)HEAD(M)='UP_NH4_4'
            IF(L.EQ.60)HEAD(M)='UP_NH4_5'
            IF(L.EQ.61)HEAD(M)='UP_NH4_6'
            IF(L.EQ.62)HEAD(M)='UP_NH4_7'
            IF(L.EQ.63)HEAD(M)='UP_NH4_8'
            IF(L.EQ.64)HEAD(M)='UP_NH4_9'
            IF(L.EQ.65)HEAD(M)='UP_NH4_10'
            IF(L.EQ.66)HEAD(M)='UP_NH4_11'
            IF(L.EQ.67)HEAD(M)='UP_NH4_12'
            IF(L.EQ.68)HEAD(M)='UP_NH4_13'
            IF(L.EQ.69)HEAD(M)='UP_NH4_14'
            IF(L.EQ.70)HEAD(M)='UP_NH4_15'
            IF(L.EQ.71)HEAD(M)='UP_NO3_1'
            IF(L.EQ.72)HEAD(M)='UP_NO3_2'
            IF(L.EQ.73)HEAD(M)='UP_NO3_3'
            IF(L.EQ.74)HEAD(M)='UP_NO3_4'
            IF(L.EQ.75)HEAD(M)='UP_NO3_5'
            IF(L.EQ.76)HEAD(M)='UP_NO3_6'
            IF(L.EQ.77)HEAD(M)='UP_NO3_7'
            IF(L.EQ.78)HEAD(M)='UP_NO3_8'
            IF(L.EQ.79)HEAD(M)='UP_NO3_9'
            IF(L.EQ.80)HEAD(M)='UP_NO3_10'
            IF(L.EQ.81)HEAD(M)='UP_NO3_11'
            IF(L.EQ.82)HEAD(M)='UP_NO3_12'
            IF(L.EQ.83)HEAD(M)='UP_NO3_13'
            IF(L.EQ.84)HEAD(M)='UP_NO3_14'
            IF(L.EQ.85)HEAD(M)='UP_NO3_15'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!hourly plant Phosphorus
      IF(N.EQ.24)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='PO4_UPTK'
            IF(L.EQ.52)HEAD(M)='[TNP]'
            IF(L.EQ.53)HEAD(M)='UP_PO4_1'
            IF(L.EQ.54)HEAD(M)='UP_PO4_2'
            IF(L.EQ.55)HEAD(M)='UP_PO4_3'
            IF(L.EQ.56)HEAD(M)='UP_PO4_4'
            IF(L.EQ.57)HEAD(M)='UP_PO4_5'
            IF(L.EQ.58)HEAD(M)='UP_PO4_6'
            IF(L.EQ.59)HEAD(M)='UP_PO4_7'
            IF(L.EQ.60)HEAD(M)='UP_PO4_8'
            IF(L.EQ.61)HEAD(M)='UP_PO4_9'
            IF(L.EQ.62)HEAD(M)='UP_PO4_10'
            IF(L.EQ.63)HEAD(M)='UP_PO4_11'
            IF(L.EQ.64)HEAD(M)='UP_PO4_12'
            IF(L.EQ.65)HEAD(M)='UP_PO4_13'
            IF(L.EQ.66)HEAD(M)='UP_PO4_14'
            IF(L.EQ.67)HEAD(M)='UP_PO4_15'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!hourly plant heat
      IF(N.EQ.25)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='CAN_RN'
            IF(L.EQ.52)HEAD(M)='CAN_LE'
            IF(L.EQ.53)HEAD(M)='CAN_H'
            IF(L.EQ.54)HEAD(M)='CAN_G'
            IF(L.EQ.55)HEAD(M)='CAN_TEMP'
            IF(L.EQ.56)HEAD(M)='TEMP_FN'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!daily plant carbon
      IF(N.EQ.26)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='SHOOT_C'
            IF(L.EQ.52)HEAD(M)='LEAF_C'
            IF(L.EQ.53)HEAD(M)='SHTH_C'
            IF(L.EQ.54)HEAD(M)='STALK_C'
            IF(L.EQ.55)HEAD(M)='RESERVE_C'
            IF(L.EQ.56)HEAD(M)='REPRO_C'
            IF(L.EQ.57)HEAD(M)='GRAIN_C'
            IF(L.EQ.58)HEAD(M)='ROOT_C'
            IF(L.EQ.59)HEAD(M)='NDL_C'
            IF(L.EQ.60)HEAD(M)='STORED_C'
            IF(L.EQ.61)HEAD(M)='GRAIN_NO.'
            IF(L.EQ.62)HEAD(M)='LAI'
            IF(L.EQ.63)HEAD(M)='GPP'
            IF(L.EQ.64)HEAD(M)='EXUD_C'
            IF(L.EQ.65)HEAD(M)='LITTER_C'
            IF(L.EQ.66)HEAD(M)='SF_LIT_C'
            IF(L.EQ.67)HEAD(M)='AUTO_RESP'
            IF(L.EQ.68)HEAD(M)='ABV_GRD_RESP'
            IF(L.EQ.69)HEAD(M)='[SOL._C]'
            IF(L.EQ.70)HEAD(M)='HVST_C'
            IF(L.EQ.71)HEAD(M)='DNS_RT_1'
            IF(L.EQ.72)HEAD(M)='DNS_RT_2'
            IF(L.EQ.73)HEAD(M)='DNS_RT_3'
            IF(L.EQ.74)HEAD(M)='DNS_RT_4'
            IF(L.EQ.75)HEAD(M)='DNS_RT_5'
            IF(L.EQ.76)HEAD(M)='DNS_RT_6'
            IF(L.EQ.77)HEAD(M)='DNS_RT_7'
            IF(L.EQ.78)HEAD(M)='DNS_RT_8'
            IF(L.EQ.79)HEAD(M)='DNS_RT_9'
            IF(L.EQ.80)HEAD(M)='DNS_RT_10'
            IF(L.EQ.81)HEAD(M)='DNS_RT_11'
            IF(L.EQ.82)HEAD(M)='DNS_RT_12'
            IF(L.EQ.83)HEAD(M)='DNS_RT_13'
            IF(L.EQ.84)HEAD(M)='DNS_RT_14'
            IF(L.EQ.85)HEAD(M)='DNS_RT_15'
            IF(L.EQ.86)HEAD(M)='C_BALANCE'
            IF(L.EQ.87)HEAD(M)='STG_DEAD_C'
            IF(L.EQ.88)HEAD(M)='FIRE_CO2'
            IF(L.EQ.89)HEAD(M)='FIRE_CH4'
            IF(L.EQ.90)HEAD(M)='NPP'
            IF(L.EQ.91)HEAD(M)='CAN_HT'
            IF(L.EQ.92)HEAD(M)='POPN'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!daily plant water
      IF(N.EQ.27)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='TRANSPN'
            IF(L.EQ.52)HEAD(M)='WTR_STRESS'
            IF(L.EQ.53)HEAD(M)='OXY_STRESS'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!daily plant nitrogen
      IF(N.EQ.28)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='SHOOT_N'
            IF(L.EQ.52)HEAD(M)='LEAF_N'
            IF(L.EQ.53)HEAD(M)='SHTH_N'
            IF(L.EQ.54)HEAD(M)='STALK_N'
            IF(L.EQ.55)HEAD(M)='RESERVE_N'
            IF(L.EQ.56)HEAD(M)='HUSK_N'
            IF(L.EQ.57)HEAD(M)='GRAIN_N'
            IF(L.EQ.58)HEAD(M)='ROOT_N'
            IF(L.EQ.59)HEAD(M)='NDL_N'
            IF(L.EQ.60)HEAD(M)='STORED_N'
            IF(L.EQ.61)HEAD(M)='UPTK_N'
            IF(L.EQ.62)HEAD(M)='LITTER_N'
            IF(L.EQ.63)HEAD(M)='TL_N_FIXED'
            IF(L.EQ.64)HEAD(M)='[SOL_N]'
            IF(L.EQ.65)HEAD(M)='N_STRESS'
            IF(L.EQ.66)HEAD(M)='[SOL_P]'
            IF(L.EQ.67)HEAD(M)='P_STRESS'
            IF(L.EQ.68)HEAD(M)='NH3_FLUX'
            IF(L.EQ.69)HEAD(M)='HVST_N'
            IF(L.EQ.70)HEAD(M)='N_BALANCE'
            IF(L.EQ.71)HEAD(M)='STG_DEAD_N'
            IF(L.EQ.72)HEAD(M)='FIRE_N'
            IF(L.EQ.73)HEAD(M)='SF_LIT_N'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!daily plant Phosphorus
      IF(N.EQ.29)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='SHOOT_P'
            IF(L.EQ.52)HEAD(M)='LEAF_P'
            IF(L.EQ.53)HEAD(M)='SHTH_P'
            IF(L.EQ.54)HEAD(M)='STALK_P'
            IF(L.EQ.55)HEAD(M)='RESERVE_P'
            IF(L.EQ.56)HEAD(M)='HUSK_P'
            IF(L.EQ.57)HEAD(M)='GRAIN_P'
            IF(L.EQ.58)HEAD(M)='ROOT_P'
            IF(L.EQ.59)HEAD(M)='NDL_P'
            IF(L.EQ.60)HEAD(M)='STORED_P'
            IF(L.EQ.61)HEAD(M)='UPTK_P'
            IF(L.EQ.62)HEAD(M)='LITTER_P'
            IF(L.EQ.63)HEAD(M)='[SOL_P]'
            IF(L.EQ.64)HEAD(M)='P_STRESS'
            IF(L.EQ.65)HEAD(M)='HVST_P'
            IF(L.EQ.66)HEAD(M)='P_BALANCE'
            IF(L.EQ.67)HEAD(M)='STG_DEAD_P'
            IF(L.EQ.68)HEAD(M)='FIRE_P'
            IF(L.EQ.69)HEAD(M)='SF_LIT_P'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
!daily plant heat
      IF(N.EQ.30)THEN
        DO L=51,100
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.51)HEAD(M)='GROWTH_STG'
            IF(L.EQ.52)HEAD(M)='BRANCH_NO.'
            IF(L.EQ.53)HEAD(M)='NODE_NO.'
            IF(L.EQ.54)HEAD(M)='RUB_ACTVN'
            IF(L.EQ.55)HEAD(M)='LEAF_N/C'
            IF(L.EQ.56)HEAD(M)='LEAF_P/C'
            IF(L.EQ.57)HEAD(M)='MIN_LWP'
            IF(L.EQ.58)HEAD(M)='O2_STRESS'
            IF(L.EQ.59)HEAD(M)='TEMP_STRESS'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTP(N-20)=M
      ENDIF
      IF(N.LE.25)WRITE(LUN,'(A12,2A8,50A16)')CDOY,DATE,HOUR,(HEAD(K),K=1,M)
      IF(N.GE.26)WRITE(LUN,'(A12,A8,50A16)')CDOY,DATE,(HEAD(K),K=1,M)
    ENDIF
  ENDDO D1010
  RETURN
  END SUBROUTINE foutp

!------------------------------------------------------------------------


  SUBROUTINE fouts(NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE OPENS AND LABELS OUTPUT FILES FOR SOIL DATA
!

  implicit none
  integer, intent(in) :: NE,NEX
  integer, intent(in) :: NHW,NHE,NVN,NVS


  character(len=*), parameter :: mod_filename = __FILE__
  CHARACTER(len=8) :: CDOY,DATE,HOUR
  CHARACTER(len=1) :: CHARM,CHARN,CHARZ,CHARO
  CHARACTER(len=2) :: CHARX,CHARY
  CHARACTER(len=4) :: CHARR
  CHARACTER(len=16):: HEAD(50)
  integer :: IDY1,IDY2,IDY,LPY,LUN,L,M,NX,NY,N,k
  real(r8):: dy


!     begin_execution


  CDOY='DOY     '
  DATE='DATE    '
  HOUR='HOUR    '
  CHARM='0'
  CHARO='0'
  CHARZ='0'
  DO NX=NHW,NHE
    DO   NY=NVN,NVS
      DO   N=1,10
        OUTFILS(N,NY,NX)= '                '
      ENDDO
    ENDDO
  ENDDO
!
!     OPEN AND NAME OUTPUT FILES
!
  DO N=21,30
    IF(DATAC(N,NE,NEX).NE.'NO')THEN
      call OPEN_safe(15,PREFIX,DATAC(N,NE,NEX),'OLD',mod_filename,__LINE__)
      WRITE(CHARR,'(I4)')IYRC
      OUTS(N-20)=CHARO//CHARR//DATAC(N,NE,NEX)
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NX.LT.10)THEN
            WRITE(CHARN,'(I1)')NX
            CHARX=CHARM//CHARN
          ELSE
            WRITE(CHARX,'(I2)')NX
          ENDIF
          IF(NY.LT.10)THEN
            WRITE(CHARN,'(I1)')NY
            CHARY=CHARM//CHARN
          ELSE
            WRITE(CHARY,'(I2)')NY
          ENDIF
          OUTFILS(N-20,NY,NX)=CHARX//CHARY//CHARZ//CHARR//DATAC(N,NE,NEX)
        ENDDO
      ENDDO
      DO M=1,2
        READ(15,*)DY
        IDY1=INT(DY/1.0E+02)   !get day, the date is given as 0101 for Jan 1st, or 3112, Dec 31th
        IDY2=INT(DY/1.0E+00-IDY1*1.0E+02) !get month
        IF(IDY2.GT.2)LPY=1
        IF(IDY2.EQ.1)THEN
          IDY=IDY1
        ELSE
          IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
        ENDIF
        IF(M.EQ.1)IDATA(N)=IDY
        IF(M.EQ.2)IDATA(N+20)=IDY
      ENDDO
!read the variable options in the form of YES/NO
      M=0
      DO while(.True.)
        M=M+1
        READ(15,'(A3)',END=4030)CHOICE(M,N-20)
      ENDDO
4030  CONTINUE
      CLOSE(15)
      LUN=N+10
      OPEN(LUN,FILE=trim(outdir)//OUTS(N-20),STATUS='UNKNOWN')
!
!     WRITE HEADINGS TO OUTPUT FILES
!
!hourly soil carbon
      M=0
      IF(N.EQ.21)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='SOIL_CO2_FLUX'
            IF(L.EQ.2)HEAD(M)='ECO_CO2_FLUX'
            IF(L.EQ.3)HEAD(M)='CH4_FLUX'
            IF(L.EQ.4)HEAD(M)='O2_FLUX'
            IF(L.EQ.5)HEAD(M)='CO2_1'
            IF(L.EQ.6)HEAD(M)='CO2_2'
            IF(L.EQ.7)HEAD(M)='CO2_3'
            IF(L.EQ.8)HEAD(M)='CO2_4'
            IF(L.EQ.9)HEAD(M)='CO2_5'
            IF(L.EQ.10)HEAD(M)='CO2_6'
            IF(L.EQ.11)HEAD(M)='CO2_7'
            IF(L.EQ.12)HEAD(M)='CO2_8'
            IF(L.EQ.13)HEAD(M)='CO2_9'
            IF(L.EQ.14)HEAD(M)='CO2_10'
            IF(L.EQ.15)HEAD(M)='CO2_11'
            IF(L.EQ.16)HEAD(M)='CO2_12'
            IF(L.EQ.17)HEAD(M)='CO2_13'
            IF(L.EQ.18)HEAD(M)='CO2_14'
            IF(L.EQ.19)HEAD(M)='CO2_LIT'
            IF(L.EQ.20)HEAD(M)='CH4_1'
            IF(L.EQ.21)HEAD(M)='CH4_2'
            IF(L.EQ.22)HEAD(M)='CH4_3'
            IF(L.EQ.23)HEAD(M)='CH4_4'
            IF(L.EQ.24)HEAD(M)='CH4_5'
            IF(L.EQ.25)HEAD(M)='CH4_6'
            IF(L.EQ.26)HEAD(M)='CH4_7'
            IF(L.EQ.27)HEAD(M)='CH4_8'
            IF(L.EQ.28)HEAD(M)='CH4_9'
            IF(L.EQ.29)HEAD(M)='CH4_10'
            IF(L.EQ.30)HEAD(M)='CH4_11'
            IF(L.EQ.31)HEAD(M)='CH4_12'
            IF(L.EQ.32)HEAD(M)='CH4_13'
            IF(L.EQ.33)HEAD(M)='CH4_14'
            IF(L.EQ.34)HEAD(M)='CH4_15'
            IF(L.EQ.35)HEAD(M)='O2_1'
            IF(L.EQ.36)HEAD(M)='O2_2'
            IF(L.EQ.37)HEAD(M)='O2_3'
            IF(L.EQ.38)HEAD(M)='O2_4'
            IF(L.EQ.39)HEAD(M)='O2_5'
            IF(L.EQ.40)HEAD(M)='O2_6'
            IF(L.EQ.41)HEAD(M)='O2_7'
            IF(L.EQ.42)HEAD(M)='O2_8'
            IF(L.EQ.43)HEAD(M)='O2_9'
            IF(L.EQ.44)HEAD(M)='O2_10'
            IF(L.EQ.45)HEAD(M)='O2_11'
            IF(L.EQ.46)HEAD(M)='O2_12'
            IF(L.EQ.47)HEAD(M)='O2_13'
            IF(L.EQ.48)HEAD(M)='O2_14'
            IF(L.EQ.49)HEAD(M)='O2_15'
            IF(L.EQ.50)HEAD(M)='O2_LIT'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!hourly soil water
      IF(N.EQ.22)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='EVAPN'
            IF(L.EQ.2)HEAD(M)='RUNOFF'
            IF(L.EQ.3)HEAD(M)='SEDIMENT'
            IF(L.EQ.4)HEAD(M)='TTL_SWC'
            IF(L.EQ.5)HEAD(M)='DISCHG'
            IF(L.EQ.6)HEAD(M)='SNOWPACK'
            IF(L.EQ.7)HEAD(M)='WTR_1'
            IF(L.EQ.8)HEAD(M)='WTR_2'
            IF(L.EQ.9)HEAD(M)='WTR_3'
            IF(L.EQ.10)HEAD(M)='WTR_4'
            IF(L.EQ.11)HEAD(M)='WTR_5'
            IF(L.EQ.12)HEAD(M)='WTR_6'
            IF(L.EQ.13)HEAD(M)='WTR_7'
            IF(L.EQ.14)HEAD(M)='WTR_8'
            IF(L.EQ.15)HEAD(M)='WTR_9'
            IF(L.EQ.16)HEAD(M)='WTR_10'
            IF(L.EQ.17)HEAD(M)='WTR_11'
            IF(L.EQ.18)HEAD(M)='WTR_12'
            IF(L.EQ.19)HEAD(M)='WTR_13'
            IF(L.EQ.20)HEAD(M)='WTR_14'
            IF(L.EQ.21)HEAD(M)='WTR_15'
            IF(L.EQ.22)HEAD(M)='WTR_16'
            IF(L.EQ.23)HEAD(M)='WTR_17'
            IF(L.EQ.24)HEAD(M)='WTR_18'
            IF(L.EQ.25)HEAD(M)='WTR_19'
            IF(L.EQ.26)HEAD(M)='WTR_20'
            IF(L.EQ.27)HEAD(M)='SURF_WTR'
            IF(L.EQ.28)HEAD(M)='ICE_1'
            IF(L.EQ.29)HEAD(M)='ICE_2'
            IF(L.EQ.30)HEAD(M)='ICE_3'
            IF(L.EQ.31)HEAD(M)='ICE_4'
            IF(L.EQ.32)HEAD(M)='ICE_5'
            IF(L.EQ.33)HEAD(M)='ICE_6'
            IF(L.EQ.34)HEAD(M)='ICE_7'
            IF(L.EQ.35)HEAD(M)='ICE_8'
            IF(L.EQ.36)HEAD(M)='ICE_9'
            IF(L.EQ.37)HEAD(M)='ICE_10'
            IF(L.EQ.38)HEAD(M)='ICE_11'
            IF(L.EQ.39)HEAD(M)='ICE_12'
            IF(L.EQ.40)HEAD(M)='ICE_13'
            IF(L.EQ.41)HEAD(M)='ICE_14'
            IF(L.EQ.42)HEAD(M)='ICE_15'
            IF(L.EQ.43)HEAD(M)='ICE_16'
            IF(L.EQ.44)HEAD(M)='ICE_17'
            IF(L.EQ.45)HEAD(M)='ICE_18'
            IF(L.EQ.46)HEAD(M)='ICE_19'
            IF(L.EQ.47)HEAD(M)='ICE_20'
            IF(L.EQ.48)HEAD(M)='SURF_ICE'
            IF(L.EQ.49)HEAD(M)='ACTV_LYR'
            IF(L.EQ.50)HEAD(M)='WTR_TBL'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!hourly soil nitrogen
      IF(N.EQ.23)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='N2O_FLUX'
            IF(L.EQ.2)HEAD(M)='N2G_FLUX'
            IF(L.EQ.3)HEAD(M)='NH3_FLUX'
            IF(L.EQ.4)HEAD(M)='SURF_N_FLUX'
            IF(L.EQ.5)HEAD(M)='SUBS_N_FLUX'
            IF(L.EQ.6)HEAD(M)='N2O_1'
            IF(L.EQ.7)HEAD(M)='N2O_2'
            IF(L.EQ.8)HEAD(M)='N2O_3'
            IF(L.EQ.9)HEAD(M)='N2O_4'
            IF(L.EQ.10)HEAD(M)='N2O_5'
            IF(L.EQ.11)HEAD(M)='N2O_6'
            IF(L.EQ.12)HEAD(M)='N2O_7'
            IF(L.EQ.13)HEAD(M)='N2O_8'
            IF(L.EQ.14)HEAD(M)='N2O_9'
            IF(L.EQ.15)HEAD(M)='N2O_10'
            IF(L.EQ.16)HEAD(M)='N2O_11'
            IF(L.EQ.17)HEAD(M)='N2O_12'
            IF(L.EQ.18)HEAD(M)='N2O_13'
            IF(L.EQ.19)HEAD(M)='N2O_14'
            IF(L.EQ.20)HEAD(M)='N2O_15'
            IF(L.EQ.21)HEAD(M)='NH3_1'
            IF(L.EQ.22)HEAD(M)='NH3_2'
            IF(L.EQ.23)HEAD(M)='NH3_3'
            IF(L.EQ.24)HEAD(M)='NH3_4'
            IF(L.EQ.25)HEAD(M)='NH3_5'
            IF(L.EQ.26)HEAD(M)='NH3_6'
            IF(L.EQ.27)HEAD(M)='NH3_7'
            IF(L.EQ.28)HEAD(M)='NH3_8'
            IF(L.EQ.29)HEAD(M)='NH3_9'
            IF(L.EQ.30)HEAD(M)='NH3_10'
            IF(L.EQ.31)HEAD(M)='NH3_11'
            IF(L.EQ.32)HEAD(M)='NH3_12'
            IF(L.EQ.33)HEAD(M)='NH3_13'
            IF(L.EQ.34)HEAD(M)='NH3_14'
            IF(L.EQ.35)HEAD(M)='NH3_15'
            IF(L.EQ.36)HEAD(M)='N2O_LIT'
            IF(L.EQ.37)HEAD(M)='NH3_LIT'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!hourly soil Phosphorus
      IF(N.EQ.24)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='SURF_P_FLUX'
            IF(L.EQ.2)HEAD(M)='SUBS_P_FLUX'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!hourly soil heat
      IF(N.EQ.25)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='SOL_RADN'
            IF(L.EQ.2)HEAD(M)='AIR_TEMP'
            IF(L.EQ.3)HEAD(M)='HUM'
            IF(L.EQ.4)HEAD(M)='WIND'
            IF(L.EQ.5)HEAD(M)='PREC'
            IF(L.EQ.6)HEAD(M)='SOIL_RN'
            IF(L.EQ.7)HEAD(M)='SOIL_LE'
            IF(L.EQ.8)HEAD(M)='SOIL_H'
            IF(L.EQ.9)HEAD(M)='SOIL_G'
            IF(L.EQ.10)HEAD(M)='ECO_RN'
            IF(L.EQ.11)HEAD(M)='ECO_LE'
            IF(L.EQ.12)HEAD(M)='ECO_H'
            IF(L.EQ.13)HEAD(M)='ECO_G'
            IF(L.EQ.14)HEAD(M)='TEMP_1'
            IF(L.EQ.15)HEAD(M)='TEMP_2'
            IF(L.EQ.16)HEAD(M)='TEMP_3'
            IF(L.EQ.17)HEAD(M)='TEMP_4'
            IF(L.EQ.18)HEAD(M)='TEMP_5'
            IF(L.EQ.19)HEAD(M)='TEMP_6'
            IF(L.EQ.20)HEAD(M)='TEMP_7'
            IF(L.EQ.21)HEAD(M)='TEMP_8'
            IF(L.EQ.22)HEAD(M)='TEMP_9'
            IF(L.EQ.23)HEAD(M)='TEMP_10'
            IF(L.EQ.24)HEAD(M)='TEMP_11'
            IF(L.EQ.25)HEAD(M)='TEMP_12'
            IF(L.EQ.26)HEAD(M)='TEMP_13'
            IF(L.EQ.27)HEAD(M)='TEMP_14'
            IF(L.EQ.28)HEAD(M)='TEMP_15'
            IF(L.EQ.29)HEAD(M)='TEMP_16'
            IF(L.EQ.30)HEAD(M)='TEMP_17'
            IF(L.EQ.31)HEAD(M)='TEMP_18'
            IF(L.EQ.32)HEAD(M)='TEMP_19'
            IF(L.EQ.33)HEAD(M)='TEMP_20'
            IF(L.EQ.34)HEAD(M)='TEMP_LITTER'
            IF(L.EQ.35)HEAD(M)='TEMP_SNOW'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!daily soil carbon
      IF(N.EQ.26)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='RESIDUE_C'
            IF(L.EQ.2)HEAD(M)='HUMUS_C'
            IF(L.EQ.3)HEAD(M)='AMENDED_C'
            IF(L.EQ.4)HEAD(M)='LITTER_C'
            IF(L.EQ.5)HEAD(M)='CO2_FLUX'
            IF(L.EQ.6)HEAD(M)='O2_FLUX'
            IF(L.EQ.7)HEAD(M)='AUTO_RESP'
            IF(L.EQ.8)HEAD(M)='MICRO_C'
            IF(L.EQ.9)HEAD(M)='SURF_RES'
            IF(L.EQ.10)HEAD(M)='CH4_FLUX'
            IF(L.EQ.11)HEAD(M)='SUR_DOC+SED_FLX'
            IF(L.EQ.12)HEAD(M)='SUB_DOC_FLX'
            IF(L.EQ.13)HEAD(M)='SUR_DIC_FLX'
            IF(L.EQ.14)HEAD(M)='SUB_DIC_FLX'
            IF(L.EQ.15)HEAD(M)='ATM_CO2'
            IF(L.EQ.16)HEAD(M)='NBP'
            IF(L.EQ.17)HEAD(M)='FIRE_CO2'
            IF(L.EQ.18)HEAD(M)='SOC_1'
            IF(L.EQ.19)HEAD(M)='SOC_2'
            IF(L.EQ.20)HEAD(M)='SOC_3'
            IF(L.EQ.21)HEAD(M)='SOC_4'
            IF(L.EQ.22)HEAD(M)='SOC_5'
            IF(L.EQ.23)HEAD(M)='SOC_6'
            IF(L.EQ.24)HEAD(M)='SOC_7'
            IF(L.EQ.25)HEAD(M)='SOC_8'
            IF(L.EQ.26)HEAD(M)='SOC_9'
            IF(L.EQ.27)HEAD(M)='SOC_10'
            IF(L.EQ.28)HEAD(M)='SOC_11'
            IF(L.EQ.29)HEAD(M)='SOC_12'
            IF(L.EQ.30)HEAD(M)='SOC_13'
            IF(L.EQ.31)HEAD(M)='SOC_14'
            IF(L.EQ.32)HEAD(M)='SOC_15'
            IF(L.EQ.41)HEAD(M)='H2_FLUX'
            IF(L.EQ.42)HEAD(M)='ECO_HVST_C'
            IF(L.EQ.43)HEAD(M)='ECO_LAI'
            IF(L.EQ.44)HEAD(M)='ECO_GPP'
            IF(L.EQ.45)HEAD(M)='ECO_RA'
            IF(L.EQ.46)HEAD(M)='ECO_NPP'
            IF(L.EQ.47)HEAD(M)='ECO_RH'
            IF(L.EQ.48)HEAD(M)='FIRE_CH4'
            IF(L.EQ.49)HEAD(M)='TTL_DIC'
            IF(L.EQ.50)HEAD(M)='STG_DEAD'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!daily soil water
      IF(N.EQ.27)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='PRECN'
            IF(L.EQ.2)HEAD(M)='ET'
            IF(L.EQ.3)HEAD(M)='RUNOFF'
            IF(L.EQ.4)HEAD(M)='WATER'
            IF(L.EQ.5)HEAD(M)='DISCHG'
            IF(L.EQ.6)HEAD(M)='SNOWPACK'
            IF(L.EQ.7)HEAD(M)='WTR_1'
            IF(L.EQ.8)HEAD(M)='WTR_2'
            IF(L.EQ.9)HEAD(M)='WTR_3'
            IF(L.EQ.10)HEAD(M)='WTR_4'
            IF(L.EQ.11)HEAD(M)='WTR_5'
            IF(L.EQ.12)HEAD(M)='WTR_6'
            IF(L.EQ.13)HEAD(M)='WTR_7'
            IF(L.EQ.14)HEAD(M)='WTR_8'
            IF(L.EQ.15)HEAD(M)='WTR_9'
            IF(L.EQ.16)HEAD(M)='WTR_10'
            IF(L.EQ.17)HEAD(M)='WTR_11'
            IF(L.EQ.18)HEAD(M)='WTR_12'
            IF(L.EQ.19)HEAD(M)='WTR_13'
            IF(L.EQ.20)HEAD(M)='SURF_WTR'
            IF(L.EQ.21)HEAD(M)='ICE_1'
            IF(L.EQ.22)HEAD(M)='ICE_2'
            IF(L.EQ.23)HEAD(M)='ICE_3'
            IF(L.EQ.24)HEAD(M)='ICE_4'
            IF(L.EQ.25)HEAD(M)='ICE_5'
            IF(L.EQ.26)HEAD(M)='ICE_6'
            IF(L.EQ.27)HEAD(M)='ICE_7'
            IF(L.EQ.28)HEAD(M)='ICE_8'
            IF(L.EQ.29)HEAD(M)='ICE_9'
            IF(L.EQ.30)HEAD(M)='ICE_10'
            IF(L.EQ.31)HEAD(M)='ICE_11'
            IF(L.EQ.32)HEAD(M)='ICE_12'
            IF(L.EQ.33)HEAD(M)='ICE_13'
            IF(L.EQ.34)HEAD(M)='SURF_ICE'
            IF(L.EQ.35)HEAD(M)='PSI_1'
            IF(L.EQ.36)HEAD(M)='PSI_2'
            IF(L.EQ.37)HEAD(M)='PSI_3'
            IF(L.EQ.38)HEAD(M)='PSI_4'
            IF(L.EQ.39)HEAD(M)='PSI_5'
            IF(L.EQ.40)HEAD(M)='PSI_6'
            IF(L.EQ.41)HEAD(M)='PSI_7'
            IF(L.EQ.42)HEAD(M)='PSI_8'
            IF(L.EQ.43)HEAD(M)='PSI_9'
            IF(L.EQ.44)HEAD(M)='PSI_10'
            IF(L.EQ.45)HEAD(M)='PSI_11'
            IF(L.EQ.46)HEAD(M)='SEDIMENT'
            IF(L.EQ.47)HEAD(M)='PSI_SURF'
            IF(L.EQ.48)HEAD(M)='SURF_ELEV'
            IF(L.EQ.49)HEAD(M)='ACTV_LYR'
            IF(L.EQ.50)HEAD(M)='WTR_TBL'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!daily soil nitrogen
      IF(N.EQ.28)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='RESIDUE_N'
            IF(L.EQ.2)HEAD(M)='HUMUS_N'
            IF(L.EQ.3)HEAD(M)='FERTZR_N'
            IF(L.EQ.4)HEAD(M)='NET_PL_EXCH_N'
            IF(L.EQ.5)HEAD(M)='NH4'
            IF(L.EQ.6)HEAD(M)='NO3'
            IF(L.EQ.7)HEAD(M)='SUR_DON+SED_FLX'
            IF(L.EQ.8)HEAD(M)='SUB_DON_FLX'
            IF(L.EQ.9)HEAD(M)='SUR_DIN_FLX'
            IF(L.EQ.10)HEAD(M)='SUB_DIN_FLX'
            IF(L.EQ.11)HEAD(M)='N2O_FLUX'
            IF(L.EQ.12)HEAD(M)='NH3_FLUX'
            IF(L.EQ.13)HEAD(M)='N2_FIXN'
            IF(L.EQ.14)HEAD(M)='MICRO_N'
            IF(L.EQ.15)HEAD(M)='NH4_1'
            IF(L.EQ.16)HEAD(M)='NH4_2'
            IF(L.EQ.17)HEAD(M)='NH4_3'
            IF(L.EQ.18)HEAD(M)='NH4_4'
            IF(L.EQ.19)HEAD(M)='NH4_5'
            IF(L.EQ.20)HEAD(M)='NH4_6'
            IF(L.EQ.21)HEAD(M)='NH4_7'
            IF(L.EQ.22)HEAD(M)='NH4_8'
            IF(L.EQ.23)HEAD(M)='NH4_9'
            IF(L.EQ.24)HEAD(M)='NH4_10'
            IF(L.EQ.25)HEAD(M)='NH4_11'
            IF(L.EQ.26)HEAD(M)='NH4_12'
            IF(L.EQ.27)HEAD(M)='NH4_13'
            IF(L.EQ.28)HEAD(M)='NH4_14'
            IF(L.EQ.29)HEAD(M)='NH4_15'
            IF(L.EQ.30)HEAD(M)='NO3_1'
            IF(L.EQ.31)HEAD(M)='NO3_2'
            IF(L.EQ.32)HEAD(M)='NO3_3'
            IF(L.EQ.33)HEAD(M)='NO3_4'
            IF(L.EQ.34)HEAD(M)='NO3_5'
            IF(L.EQ.35)HEAD(M)='NO3_6'
            IF(L.EQ.36)HEAD(M)='NO3_7'
            IF(L.EQ.37)HEAD(M)='NO3_8'
            IF(L.EQ.38)HEAD(M)='NO3_9'
            IF(L.EQ.39)HEAD(M)='NO3_10'
            IF(L.EQ.40)HEAD(M)='NO3_11'
            IF(L.EQ.41)HEAD(M)='NO3_12'
            IF(L.EQ.42)HEAD(M)='NO3_13'
            IF(L.EQ.43)HEAD(M)='NO3_14'
            IF(L.EQ.44)HEAD(M)='NO3_15'
            IF(L.EQ.45)HEAD(M)='NH4_RES'
            IF(L.EQ.46)HEAD(M)='NO3_RES'
            IF(L.EQ.47)HEAD(M)='ECO_HVST_N'
            IF(L.EQ.48)HEAD(M)='NET_N_MIN'
            IF(L.EQ.49)HEAD(M)='FIRE_N'
            IF(L.EQ.50)HEAD(M)='N2_FLUX'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!daily soil phosphorus
      IF(N.EQ.29)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='RESIDUE_P'
            IF(L.EQ.2)HEAD(M)='HUMUS_P'
            IF(L.EQ.3)HEAD(M)='FERTZR_P'
            IF(L.EQ.4)HEAD(M)='NET_PL_EXCH_P'
            IF(L.EQ.5)HEAD(M)='EXCH_PO4'
            IF(L.EQ.6)HEAD(M)='SUR_DOP+SED_FLX'
            IF(L.EQ.7)HEAD(M)='SUB_DOP_FLX'
            IF(L.EQ.8)HEAD(M)='SUR_DIP_FLX'
            IF(L.EQ.9)HEAD(M)='SUB_DIP_FLX'
            IF(L.EQ.10)HEAD(M)='PRECIP_P'
            IF(L.EQ.11)HEAD(M)='MICRO_P'
            IF(L.EQ.12)HEAD(M)='FIRE_P'
            IF(L.EQ.13)HEAD(M)='PO4_1'
            IF(L.EQ.14)HEAD(M)='PO4_2'
            IF(L.EQ.15)HEAD(M)='PO4_3'
            IF(L.EQ.16)HEAD(M)='PO4_4'
            IF(L.EQ.17)HEAD(M)='PO4_5'
            IF(L.EQ.18)HEAD(M)='PO4_6'
            IF(L.EQ.19)HEAD(M)='PO4_7'
            IF(L.EQ.20)HEAD(M)='PO4_8'
            IF(L.EQ.21)HEAD(M)='PO4_9'
            IF(L.EQ.22)HEAD(M)='PO4_10'
            IF(L.EQ.23)HEAD(M)='PO4_11'
            IF(L.EQ.24)HEAD(M)='PO4_12'
            IF(L.EQ.25)HEAD(M)='PO4_13'
            IF(L.EQ.26)HEAD(M)='PO4_14'
            IF(L.EQ.27)HEAD(M)='PO4_15'
            IF(L.EQ.28)HEAD(M)='EXCH_P_1'
            IF(L.EQ.29)HEAD(M)='EXCH_P_2'
            IF(L.EQ.30)HEAD(M)='EXCH_P_3'
            IF(L.EQ.31)HEAD(M)='EXCH_P_4'
            IF(L.EQ.32)HEAD(M)='EXCH_P_5'
            IF(L.EQ.33)HEAD(M)='EXCH_P_6'
            IF(L.EQ.34)HEAD(M)='EXCH_P_7'
            IF(L.EQ.35)HEAD(M)='EXCH_P_8'
            IF(L.EQ.36)HEAD(M)='EXCH_P_9'
            IF(L.EQ.37)HEAD(M)='EXCH_P_10'
            IF(L.EQ.38)HEAD(M)='EXCH_P_11'
            IF(L.EQ.39)HEAD(M)='EXCH_P_12'
            IF(L.EQ.40)HEAD(M)='EXCH_P_13'
            IF(L.EQ.41)HEAD(M)='EXCH_P_14'
            IF(L.EQ.42)HEAD(M)='EXCH_P_15'
            IF(L.EQ.43)HEAD(M)='PO4_RES'
            IF(L.EQ.44)HEAD(M)='EXCH_P_RES'
            IF(L.EQ.47)HEAD(M)='ECO_HVST_P'
            IF(L.EQ.48)HEAD(M)='NET_P_MIN'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
!daily soil heat
      IF(N.EQ.30)THEN
        DO L=1,50
          IF(CHOICE(L,N-20).EQ.'YES')THEN
            M=M+1
            IF(L.EQ.1)HEAD(M)='RADN'
            IF(L.EQ.2)HEAD(M)='TMAX_AIR'
            IF(L.EQ.3)HEAD(M)='TMIN_AIR'
            IF(L.EQ.4)HEAD(M)='HMAX_AIR'
            IF(L.EQ.5)HEAD(M)='HMIN_AIR'
            IF(L.EQ.6)HEAD(M)='WIND'
            IF(L.EQ.7)HEAD(M)='PRECN'
            IF(L.EQ.8)HEAD(M)='TMAX_SOIL_1'
            IF(L.EQ.9)HEAD(M)='TMIN_SOIL_1'
            IF(L.EQ.10)HEAD(M)='TMAX_SOIL_2'
            IF(L.EQ.11)HEAD(M)='TMIN_SOIL_2'
            IF(L.EQ.12)HEAD(M)='TMAX_SOIL_3'
            IF(L.EQ.13)HEAD(M)='TMIN_SOIL_3'
            IF(L.EQ.14)HEAD(M)='TMAX_SOIL_4'
            IF(L.EQ.15)HEAD(M)='TMIN_SOIL_4'
            IF(L.EQ.16)HEAD(M)='TMAX_SOIL_5'
            IF(L.EQ.17)HEAD(M)='TMIN_SOIL_5'
            IF(L.EQ.18)HEAD(M)='TMAX_SOIL_6'
            IF(L.EQ.19)HEAD(M)='TMIN_SOIL_6'
            IF(L.EQ.20)HEAD(M)='TMAX_SOIL_7'
            IF(L.EQ.21)HEAD(M)='TMIN_SOIL_7'
            IF(L.EQ.22)HEAD(M)='TMAX_SOIL_8'
            IF(L.EQ.23)HEAD(M)='TMIN_SOIL_8'
            IF(L.EQ.24)HEAD(M)='TMAX_SOIL_9'
            IF(L.EQ.25)HEAD(M)='TMIN_SOIL_9'
            IF(L.EQ.26)HEAD(M)='TMAX_SOIL_10'
            IF(L.EQ.27)HEAD(M)='TMIN_SOIL_10'
            IF(L.EQ.28)HEAD(M)='TMAX_SOIL_11'
            IF(L.EQ.29)HEAD(M)='TMIN_SOIL_11'
            IF(L.EQ.30)HEAD(M)='TMAX_SOIL_12'
            IF(L.EQ.31)HEAD(M)='TMIN_SOIL_12'
            IF(L.EQ.32)HEAD(M)='TMAX_SOIL_13'
            IF(L.EQ.33)HEAD(M)='TMIN_SOIL_13'
            IF(L.EQ.34)HEAD(M)='TMAX_SOIL_14'
            IF(L.EQ.35)HEAD(M)='TMIN_SOIL_14'
            IF(L.EQ.36)HEAD(M)='TMAX_LITTER'
            IF(L.EQ.37)HEAD(M)='TMIN_LITTER'
            IF(L.EQ.38)HEAD(M)='ECND_1'
            IF(L.EQ.39)HEAD(M)='ECND_2'
            IF(L.EQ.40)HEAD(M)='ECND_3'
            IF(L.EQ.41)HEAD(M)='ECND_4'
            IF(L.EQ.42)HEAD(M)='ECND_5'
            IF(L.EQ.43)HEAD(M)='ECND_6'
            IF(L.EQ.44)HEAD(M)='ECND_7'
            IF(L.EQ.45)HEAD(M)='ECND_8'
            IF(L.EQ.46)HEAD(M)='ECND_9'
            IF(L.EQ.47)HEAD(M)='ECND_10'
            IF(L.EQ.48)HEAD(M)='ECND_11'
            IF(L.EQ.49)HEAD(M)='ECND_12'
            IF(L.EQ.50)HEAD(M)='TTL_SALT_DISCHG'
!            print*,M,HEAD(M)
          ENDIF
        ENDDO
        NOUTS(N-20)=M
      ENDIF
      IF(N.LE.25)WRITE(LUN,'(A12,2A8,50A16)')CDOY,DATE,HOUR,(HEAD(K),K=1,M)
      IF(N.GE.26)WRITE(LUN,'(A12,A8,50A16)')CDOY,DATE,(HEAD(K),K=1,M)
    ENDIF
  ENDDO
  RETURN
  END SUBROUTINE fouts


!------------------------------------------------------------------------

  SUBROUTINE outpd(I,NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES DAILY OUTPUT FOR PLANT
!     C, N, P, WATER AND HEAT TO OUTPUT FILES DEPENDING
!     ON SELECTIONS MADE IN OUTPUT CONTROL FILES IN
!     THE RUN SCRIPT
!

  implicit none
  integer, intent(in) :: I,NE
  integer, intent(in) :: NEX,NHW,NHE,NVN,NVS


  CHARACTER(len=16) :: CHEAD
  real(r8) :: HEAD(50)
  integer :: IHEAD,K,KN,LUN,M,N,NX,NY,NZ
!     execution begins here
! second half of output
  DO N=26,30
    IF(DATAC(N,NE,NEX).NE.'NO')THEN
      DO M=1,50
        HEAD(M)=0.0
      ENDDO
      LUN=N+20
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP0(NY,NX)
            IF(IFLGC(NZ,NY,NX).EQ.1)THEN
!
!     WRITE DAILY CROP DATA TO OUTPUT FILES
!
              M=0
!
!     DAILY PLANT C OUTPUT
!
              IF(N.EQ.26)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=WTSHTE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=WTLFE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.53)HEAD(M)=WTSHEE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.54)HEAD(M)=WTSTKE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.55)HEAD(M)=WTRSVE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.56)HEAD(M)=(WTHSKE(ielmc,NZ,NY,NX)+WTEARE(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.57)HEAD(M)=WTGRE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.58)HEAD(M)=WTRTE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.59)HEAD(M)=WTNDE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.60)HEAD(M)=WTRVE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.61)HEAD(M)=GRNO(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.62)HEAD(M)=ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.63)HEAD(M)=CARBN(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.64)HEAD(M)=TEUPTK(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.65)HEAD(M)=TESNC(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.66)HEAD(M)=TESN0(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.67)HEAD(M)=TCO2T(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.68)HEAD(M)=TCO2A(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.69)HEAD(M)=CEPOLP(ielmc,NZ,NY,NX)
                      IF(K.EQ.70)HEAD(M)=HVSTE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.71)HEAD(M)=RTDNP(1,1,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.72)HEAD(M)=RTDNP(1,2,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.73)HEAD(M)=RTDNP(1,3,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.74)HEAD(M)=RTDNP(1,4,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.75)HEAD(M)=RTDNP(1,5,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.76)HEAD(M)=RTDNP(1,6,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.77)HEAD(M)=RTDNP(1,7,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.78)HEAD(M)=RTDNP(1,8,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.79)HEAD(M)=RTDNP(1,9,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.80)HEAD(M)=RTDNP(1,10,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.81)HEAD(M)=RTDNP(1,11,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.82)HEAD(M)=RTDNP(1,12,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.83)HEAD(M)=RTDNP(1,13,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.84)HEAD(M)=RTDNP(1,14,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.85)HEAD(M)=RTDNP(1,15,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.86)HEAD(M)=BALE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.87)HEAD(M)=WTSTGE(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.88)HEAD(M)=VCO2F(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.89)HEAD(M)=VCH4F(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.90)HEAD(M)=ZNPP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.91)HEAD(M)=ZC(NZ,NY,NX)
                      IF(K.EQ.92)HEAD(M)=PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     DAILY PLANT WATER OUTPUT
!
              IF(N.EQ.27)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=-CTRAN(NZ,NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=WSTR(NZ,NY,NX)
                      IF(K.EQ.53)HEAD(M)=OSTR(NZ,NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     DAILY PLANT N OUTPUT
!
              IF(N.EQ.28)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=WTSHTE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=WTLFE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.53)HEAD(M)=WTSHEE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.54)HEAD(M)=WTSTKE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.55)HEAD(M)=WTRSVE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.56)HEAD(M)=(WTHSKE(ielmn,NZ,NY,NX)+WTEARE(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.57)HEAD(M)=WTGRE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.58)HEAD(M)=WTRTE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.59)HEAD(M)=WTNDE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.60)HEAD(M)=WTRVE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.61)HEAD(M)=TEUPTK(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.62)HEAD(M)=TESNC(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.63)HEAD(M)=TZUPFX(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.64)HEAD(M)=CEPOLP(ielmn,NZ,NY,NX)
                      IF(K.EQ.65)THEN
                        IF(WTLFE(ielmc,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                          HEAD(M)=(WTLFE(ielmn,NZ,NY,NX)+EPOOLP(ielmn,NZ,NY,NX))/(WTLFE(ielmc,NZ,NY,NX)+EPOOLP(ielmc,NZ,NY,NX))
                        ELSE
                          HEAD(M)=0.0
                        ENDIF
                      ENDIF
                      IF(K.EQ.66)HEAD(M)=CEPOLP(ielmp,NZ,NY,NX)
                      IF(K.EQ.67)THEN
                        IF(WTLFE(ielmc,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                          HEAD(M)=(WTLFE(ielmn,NZ,NY,NX)+EPOOLP(ielmp,NZ,NY,NX))/(WTLFE(ielmc,NZ,NY,NX)+EPOOLP(ielmc,NZ,NY,NX))
                        ELSE
                          HEAD(M)=0.0
                        ENDIF
                      ENDIF
                      IF(K.EQ.68)HEAD(M)=TNH3C(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.69)HEAD(M)=HVSTE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.70)HEAD(M)=BALE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.71)HEAD(M)=WTSTGE(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.72)HEAD(M)=VNH3F(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.73)HEAD(M)=TESN0(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     DAILY PLANT P OUTPUT
!
              IF(N.EQ.29)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=WTSHTE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=WTLFE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.53)HEAD(M)=WTSHEE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.54)HEAD(M)=WTSTKE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.55)HEAD(M)=WTRSVE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.56)HEAD(M)=(WTHSKE(ielmp,NZ,NY,NX)+WTEARE(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.57)HEAD(M)=WTGRE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.58)HEAD(M)=WTRTE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.59)HEAD(M)=WTNDE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.60)HEAD(M)=WTRVE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.61)HEAD(M)=TEUPTK(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.62)HEAD(M)=TESNC(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.63)HEAD(M)=CEPOLP(ielmp,NZ,NY,NX)
                      IF(K.EQ.64)THEN
                        IF(WTLFE(ielmc,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                          HEAD(M)=(WTLFE(ielmp,NZ,NY,NX)+EPOOLP(ielmp,NZ,NY,NX))/(WTLFE(ielmc,NZ,NY,NX)+EPOOLP(ielmc,NZ,NY,NX))
                        ELSE
                          HEAD(M)=0.0
                        ENDIF
                      ENDIF
                      IF(K.EQ.65)HEAD(M)=HVSTE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.66)HEAD(M)=BALE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.67)HEAD(M)=WTSTGE(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.68)HEAD(M)=VPO4F(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.69)HEAD(M)=TESN0(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     DAILY PLANT HEAT OUTPUT
!
              IF(N.EQ.30)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  CHEAD=' '
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)THEN
                        DO KN=10,0,-1
                          IF(KN.EQ.0)THEN
                            CHEAD='  PLANTING'
                            exit
                          ELSE
                            IF(IDAY(KN,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
                              CYCLE
                            ELSE
                              IF(KN.EQ.1)CHEAD='  EMERGENCE'
                              IF(KN.EQ.2)CHEAD='  FLORAL_INIT.'
                              IF(KN.EQ.3)CHEAD='  JOINTING'
                              IF(KN.EQ.4)CHEAD='  ELONGATION'
                              IF(KN.EQ.5)CHEAD='  HEADING'
                              IF(KN.EQ.6)CHEAD='  ANTHESIS'
                              IF(KN.EQ.7)CHEAD='  SEED_FILL'
                              IF(KN.EQ.8)CHEAD='  SEED_NO._SET'
                              IF(KN.EQ.9)CHEAD='  SEED_MASS_SET'
                              IF(KN.EQ.10)CHEAD='  END_SEED_FILL'
                              exit
                            ENDIF
                          ENDIF
                        ENDDO
                      ENDIF
                      IF(K.EQ.52)IHEAD=NBR(NZ,NY,NX)
                      IF(K.EQ.53)HEAD(M)=VSTG(NB1(NZ,NY,NX),NZ,NY,NX)
                      IF(K.EQ.54)HEAD(M)=FDBK(NB1(NZ,NY,NX),NZ,NY,NX)
                      IF(K.EQ.55)THEN
                        IF(WTLFE(ielmc,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                          HEAD(M)=(WTLFE(ielmn,NZ,NY,NX)+EPOOLP(ielmn,NZ,NY,NX))/(WTLFE(ielmc,NZ,NY,NX)+EPOOLP(ielmc,NZ,NY,NX))
                        ELSE
                          HEAD(M)=0.0_r8
                        ENDIF
                      ENDIF
                      IF(K.EQ.56)THEN
                        IF(WTLFE(ielmc,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                          HEAD(M)=(WTLFE(ielmp,NZ,NY,NX)+EPOOLP(ielmp,NZ,NY,NX))/(WTLFE(ielmc,NZ,NY,NX)+EPOOLP(ielmc,NZ,NY,NX))
                        ELSE
                          HEAD(M)=0.0_r8
                        ENDIF
                      ENDIF
                      IF(K.EQ.57)HEAD(M)=PSILZ(NZ,NY,NX)
                      IF(K.EQ.58)HEAD(M)=OSTR(NZ,NY,NX)
                      IF(K.EQ.59)HEAD(M)=TFN3(NZ,NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,A16,I16,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,CHEAD,IHEAD,(HEAD(K),K=3,M)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  RETURN
  END SUBROUTINE outpd

!------------------------------------------------------------------------


  SUBROUTINE outph(I,J,NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES HOURLY OUTPUT FOR PLANT
!     C, N, P, WATER AND HEAT TO OUTPUT FILES DEPENDING
!     ON SELECTIONS MADE IN OUTPUT CONTROL FILES IN
!     THE RUN SCRIPT
!

  implicit none
  integer, intent(in) :: I,J,NE
  integer, intent(in) :: NEX,NHW,NHE,NVN,NVS


  real(r8) :: HEAD(50)
  integer :: K,LUN,M,N,NX,NY,NZ
!     execution begins here
! first half of output
  DO N=21,25
    IF(DATAC(N,NE,NEX).NE.'NO')THEN
      DO M=1,50
        HEAD(M)=0.0
      ENDDO
      LUN=N+20
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            IF(IFLGC(NZ,NY,NX).EQ.1)THEN
!
!     WRITE HOURLY CROP DATA TO OUTPUT FILES
!
              M=0
!
!     HOURLY PLANT C OUTPUT
!
              IF(N.EQ.21)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO  K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=CNET(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.148
                      IF(K.EQ.52)HEAD(M)=CARBN(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.53)HEAD(M)=TCO2A(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.54)HEAD(M)=CEPOLP(ielmc,NZ,NY,NX)
                      IF(K.EQ.55)HEAD(M)=RC(NZ,NY,NX)*1.56*3600.0
                      IF(K.EQ.56)HEAD(M)=RA(NZ,NY,NX)*1.34*3600.0
                      IF(K.EQ.57)HEAD(M)=CO2Q(NZ,NY,NX)
                      IF(K.EQ.58)HEAD(M)=ARLFS(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     HOURLY PLANT WATER OUTPUT
!
              IF(N.EQ.22)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=PSILT(NZ,NY,NX)
                      IF(K.EQ.52)HEAD(M)=PSILG(NZ,NY,NX)
                      IF(K.EQ.53)HEAD(M)=RC(NZ,NY,NX)*3600.0_r8
                      IF(K.EQ.54)HEAD(M)=RA(NZ,NY,NX)*3600.0_r8
                      IF(K.EQ.55)HEAD(M)=EP(NZ,NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.56)HEAD(M)=OSTR(NZ,NY,NX)
                      IF(K.EQ.57)HEAD(M)=PSIRT(1,1,NZ,NY,NX)
                      IF(K.EQ.58)HEAD(M)=PSIRT(1,2,NZ,NY,NX)
                      IF(K.EQ.59)HEAD(M)=PSIRT(1,3,NZ,NY,NX)
                      IF(K.EQ.60)HEAD(M)=PSIRT(1,4,NZ,NY,NX)
                      IF(K.EQ.61)HEAD(M)=PSIRT(1,5,NZ,NY,NX)
                      IF(K.EQ.62)HEAD(M)=PSIRT(1,6,NZ,NY,NX)
                      IF(K.EQ.63)HEAD(M)=PSIRT(1,7,NZ,NY,NX)
                      IF(K.EQ.64)HEAD(M)=PSIRT(1,8,NZ,NY,NX)
                      IF(K.EQ.65)HEAD(M)=PSIRT(1,9,NZ,NY,NX)
                      IF(K.EQ.66.AND.JZ>=10)HEAD(M)=PSIRT(1,10,NZ,NY,NX)
                      IF(K.EQ.67.AND.JZ>=11)HEAD(M)=PSIRT(1,11,NZ,NY,NX)
                      IF(K.EQ.68.AND.JZ>=12)HEAD(M)=PSIRT(1,12,NZ,NY,NX)
                      IF(K.EQ.69.AND.JZ>=13)HEAD(M)=PSIRT(1,13,NZ,NY,NX)
                      IF(K.EQ.70.AND.JZ>=14)HEAD(M)=PSIRT(1,14,NZ,NY,NX)
                      IF(K.EQ.71.AND.JZ>=15)HEAD(M)=PSIRT(1,15,NZ,NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     HOURLY PLANT N OUTPUT
!
              IF(N.EQ.23)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=UPNH4(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=UPNO3(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.53)HEAD(M)=UPNF(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.54)HEAD(M)=CEPOLP(ielmn,NZ,NY,NX)
                      IF(K.EQ.55)HEAD(M)=RNH3C(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.56)HEAD(M)=(RUPNH4(1,1,NZ,NY,NX)+RUPNH4(2,1,NZ,NY,NX) &
                        +RUPNHB(1,1,NZ,NY,NX)+RUPNHB(2,1,NZ,NY,NX))/AREA(3,1,NY,NX)
                      IF(K.EQ.57)HEAD(M)=(RUPNH4(1,2,NZ,NY,NX)+RUPNH4(2,2,NZ,NY,NX) &
                        +RUPNHB(1,2,NZ,NY,NX)+RUPNHB(2,2,NZ,NY,NX))/AREA(3,2,NY,NX)
                      IF(K.EQ.58)HEAD(M)=(RUPNH4(1,3,NZ,NY,NX)+RUPNH4(2,3,NZ,NY,NX) &
                        +RUPNHB(1,3,NZ,NY,NX)+RUPNHB(2,3,NZ,NY,NX))/AREA(3,3,NY,NX)
                      IF(K.EQ.59)HEAD(M)=(RUPNH4(1,4,NZ,NY,NX)+RUPNH4(2,4,NZ,NY,NX) &
                        +RUPNHB(1,4,NZ,NY,NX)+RUPNHB(2,4,NZ,NY,NX))/AREA(3,4,NY,NX)
                      IF(K.EQ.60)HEAD(M)=(RUPNH4(1,5,NZ,NY,NX)+RUPNH4(2,5,NZ,NY,NX) &
                        +RUPNHB(1,5,NZ,NY,NX)+RUPNHB(2,5,NZ,NY,NX))/AREA(3,5,NY,NX)
                      IF(K.EQ.61)HEAD(M)=(RUPNH4(1,6,NZ,NY,NX)+RUPNH4(2,6,NZ,NY,NX) &
                        +RUPNHB(1,6,NZ,NY,NX)+RUPNHB(2,6,NZ,NY,NX))/AREA(3,6,NY,NX)
                      IF(K.EQ.62)HEAD(M)=(RUPNH4(1,7,NZ,NY,NX)+RUPNH4(2,7,NZ,NY,NX) &
                        +RUPNHB(1,7,NZ,NY,NX)+RUPNHB(2,7,NZ,NY,NX))/AREA(3,7,NY,NX)
                      IF(K.EQ.63)HEAD(M)=(RUPNH4(1,8,NZ,NY,NX)+RUPNH4(2,8,NZ,NY,NX) &
                        +RUPNHB(1,8,NZ,NY,NX)+RUPNHB(2,8,NZ,NY,NX))/AREA(3,8,NY,NX)
                      IF(K.EQ.64)HEAD(M)=(RUPNH4(1,9,NZ,NY,NX)+RUPNH4(2,9,NZ,NY,NX) &
                        +RUPNHB(1,9,NZ,NY,NX)+RUPNHB(2,9,NZ,NY,NX))/AREA(3,9,NY,NX)
                      IF(K.EQ.65.AND.JZ>=10)HEAD(M)=(RUPNH4(1,10,NZ,NY,NX)+RUPNH4(2,10,NZ,NY,NX) &
                        +RUPNHB(1,10,NZ,NY,NX)+RUPNHB(2,10,NZ,NY,NX))/AREA(3,10,NY,NX)
                      IF(K.EQ.66.AND.JZ>=11)HEAD(M)=(RUPNH4(1,11,NZ,NY,NX)+RUPNH4(2,11,NZ,NY,NX) &
                        +RUPNHB(1,11,NZ,NY,NX)+RUPNHB(2,11,NZ,NY,NX))/AREA(3,11,NY,NX)
                      IF(K.EQ.67.AND.JZ>=12)HEAD(M)=(RUPNH4(1,12,NZ,NY,NX)+RUPNH4(2,12,NZ,NY,NX) &
                        +RUPNHB(1,12,NZ,NY,NX)+RUPNHB(2,12,NZ,NY,NX))/AREA(3,12,NY,NX)
                      IF(K.EQ.68.AND.JZ>=13)HEAD(M)=(RUPNH4(1,13,NZ,NY,NX)+RUPNH4(2,13,NZ,NY,NX) &
                        +RUPNHB(1,13,NZ,NY,NX)+RUPNHB(2,13,NZ,NY,NX))/AREA(3,13,NY,NX)
                      IF(K.EQ.69.AND.JZ>=14)HEAD(M)=(RUPNH4(1,14,NZ,NY,NX)+RUPNH4(2,14,NZ,NY,NX) &
                        +RUPNHB(1,14,NZ,NY,NX)+RUPNHB(2,14,NZ,NY,NX))/AREA(3,14,NY,NX)
                      IF(K.EQ.70.AND.JZ>=15)HEAD(M)=(RUPNH4(1,15,NZ,NY,NX)+RUPNH4(2,15,NZ,NY,NX) &
                        +RUPNHB(1,15,NZ,NY,NX)+RUPNHB(2,15,NZ,NY,NX))/AREA(3,15,NY,NX)

                      IF(K.EQ.71)HEAD(M)=(RUPNO3(1,1,NZ,NY,NX)+RUPNO3(2,1,NZ,NY,NX) &
                        +RUPNOB(1,1,NZ,NY,NX)+RUPNOB(2,1,NZ,NY,NX))/AREA(3,1,NY,NX)
                      IF(K.EQ.72)HEAD(M)=(RUPNO3(1,2,NZ,NY,NX)+RUPNO3(2,2,NZ,NY,NX) &
                        +RUPNOB(1,2,NZ,NY,NX)+RUPNOB(2,2,NZ,NY,NX))/AREA(3,2,NY,NX)
                      IF(K.EQ.73)HEAD(M)=(RUPNO3(1,3,NZ,NY,NX)+RUPNO3(2,3,NZ,NY,NX) &
                        +RUPNOB(1,3,NZ,NY,NX)+RUPNOB(2,3,NZ,NY,NX))/AREA(3,3,NY,NX)
                      IF(K.EQ.74)HEAD(M)=(RUPNO3(1,4,NZ,NY,NX)+RUPNO3(2,4,NZ,NY,NX) &
                        +RUPNOB(1,4,NZ,NY,NX)+RUPNOB(2,4,NZ,NY,NX))/AREA(3,4,NY,NX)
                      IF(K.EQ.75)HEAD(M)=(RUPNO3(1,5,NZ,NY,NX)+RUPNO3(2,5,NZ,NY,NX) &
                        +RUPNOB(1,5,NZ,NY,NX)+RUPNOB(2,5,NZ,NY,NX))/AREA(3,5,NY,NX)
                      IF(K.EQ.76)HEAD(M)=(RUPNO3(1,6,NZ,NY,NX)+RUPNO3(2,6,NZ,NY,NX) &
                        +RUPNOB(1,6,NZ,NY,NX)+RUPNOB(2,6,NZ,NY,NX))/AREA(3,6,NY,NX)
                      IF(K.EQ.77)HEAD(M)=(RUPNO3(1,7,NZ,NY,NX)+RUPNO3(2,7,NZ,NY,NX) &
                        +RUPNOB(1,7,NZ,NY,NX)+RUPNOB(2,7,NZ,NY,NX))/AREA(3,7,NY,NX)
                      IF(K.EQ.78)HEAD(M)=(RUPNO3(1,8,NZ,NY,NX)+RUPNO3(2,8,NZ,NY,NX) &
                        +RUPNOB(1,8,NZ,NY,NX)+RUPNOB(2,8,NZ,NY,NX))/AREA(3,8,NY,NX)
                      IF(K.EQ.79)HEAD(M)=(RUPNO3(1,9,NZ,NY,NX)+RUPNO3(2,9,NZ,NY,NX) &
                        +RUPNOB(1,9,NZ,NY,NX)+RUPNOB(2,9,NZ,NY,NX))/AREA(3,9,NY,NX)
                      IF(K.EQ.80.AND.JZ>=10)HEAD(M)=(RUPNO3(1,10,NZ,NY,NX)+RUPNO3(2,10,NZ,NY,NX) &
                        +RUPNOB(1,10,NZ,NY,NX)+RUPNOB(2,10,NZ,NY,NX))/AREA(3,10,NY,NX)
                      IF(K.EQ.81.AND.JZ>=11)HEAD(M)=(RUPNO3(1,11,NZ,NY,NX)+RUPNO3(2,11,NZ,NY,NX) &
                        +RUPNOB(1,11,NZ,NY,NX)+RUPNOB(2,11,NZ,NY,NX))/AREA(3,11,NY,NX)
                      IF(K.EQ.82.AND.JZ>=12)HEAD(M)=(RUPNO3(1,12,NZ,NY,NX)+RUPNO3(2,12,NZ,NY,NX) &
                        +RUPNOB(1,12,NZ,NY,NX)+RUPNOB(2,12,NZ,NY,NX))/AREA(3,12,NY,NX)
                      IF(K.EQ.83.AND.JZ>=13)HEAD(M)=(RUPNO3(1,13,NZ,NY,NX)+RUPNO3(2,13,NZ,NY,NX) &
                        +RUPNOB(1,13,NZ,NY,NX)+RUPNOB(2,13,NZ,NY,NX))/AREA(3,13,NY,NX)
                      IF(K.EQ.84.AND.JZ>=14)HEAD(M)=(RUPNO3(1,14,NZ,NY,NX)+RUPNO3(2,14,NZ,NY,NX) &
                        +RUPNOB(1,14,NZ,NY,NX)+RUPNOB(2,14,NZ,NY,NX))/AREA(3,14,NY,NX)
                      IF(K.EQ.85.AND.JZ>=15)HEAD(M)=(RUPNO3(1,15,NZ,NY,NX)+RUPNO3(2,15,NZ,NY,NX) &
                        +RUPNOB(1,15,NZ,NY,NX)+RUPNOB(2,15,NZ,NY,NX))/AREA(3,15,NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     HOURLY PLANT P OUTPUT
!
              IF(N.EQ.24)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=UPH2P(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=CEPOLP(ielmp,NZ,NY,NX)
                      IF(K.EQ.53)HEAD(M)=(RUPH2P(1,1,NZ,NY,NX)+RUPH2P(2,1,NZ,NY,NX) &
                        +RUPH2B(1,1,NZ,NY,NX)+RUPH2B(2,1,NZ,NY,NX))/AREA(3,1,NY,NX)
                      IF(K.EQ.54)HEAD(M)=(RUPH2P(1,2,NZ,NY,NX)+RUPH2P(2,2,NZ,NY,NX) &
                        +RUPH2B(1,2,NZ,NY,NX)+RUPH2B(2,2,NZ,NY,NX))/AREA(3,2,NY,NX)
                      IF(K.EQ.55)HEAD(M)=(RUPH2P(1,3,NZ,NY,NX)+RUPH2P(2,3,NZ,NY,NX) &
                        +RUPH2B(1,3,NZ,NY,NX)+RUPH2B(2,3,NZ,NY,NX))/AREA(3,3,NY,NX)
                      IF(K.EQ.56)HEAD(M)=(RUPH2P(1,4,NZ,NY,NX)+RUPH2P(2,4,NZ,NY,NX) &
                        +RUPH2B(1,4,NZ,NY,NX)+RUPH2B(2,4,NZ,NY,NX))/AREA(3,4,NY,NX)
                      IF(K.EQ.57)HEAD(M)=(RUPH2P(1,5,NZ,NY,NX)+RUPH2P(2,5,NZ,NY,NX) &
                        +RUPH2B(1,5,NZ,NY,NX)+RUPH2B(2,5,NZ,NY,NX))/AREA(3,5,NY,NX)
                      IF(K.EQ.58)HEAD(M)=(RUPH2P(1,6,NZ,NY,NX)+RUPH2P(2,6,NZ,NY,NX) &
                        +RUPH2B(1,6,NZ,NY,NX)+RUPH2B(2,6,NZ,NY,NX))/AREA(3,6,NY,NX)
                      IF(K.EQ.59)HEAD(M)=(RUPH2P(1,7,NZ,NY,NX)+RUPH2P(2,7,NZ,NY,NX) &
                        +RUPH2B(1,7,NZ,NY,NX)+RUPH2B(2,7,NZ,NY,NX))/AREA(3,7,NY,NX)
                      IF(K.EQ.60)HEAD(M)=(RUPH2P(1,8,NZ,NY,NX)+RUPH2P(2,8,NZ,NY,NX) &
                        +RUPH2B(1,8,NZ,NY,NX)+RUPH2B(2,8,NZ,NY,NX))/AREA(3,8,NY,NX)
                      IF(K.EQ.61)HEAD(M)=(RUPH2P(1,9,NZ,NY,NX)+RUPH2P(2,9,NZ,NY,NX) &
                        +RUPH2B(1,9,NZ,NY,NX)+RUPH2B(2,9,NZ,NY,NX))/AREA(3,9,NY,NX)
                      IF(K.EQ.62)HEAD(M)=(RUPH2P(1,10,NZ,NY,NX)+RUPH2P(2,10,NZ,NY,NX) &
                        +RUPH2B(1,10,NZ,NY,NX)+RUPH2B(2,10,NZ,NY,NX))/AREA(3,10,NY,NX)
                      IF(K.EQ.63)HEAD(M)=(RUPH2P(1,11,NZ,NY,NX)+RUPH2P(2,11,NZ,NY,NX) &
                        +RUPH2B(1,11,NZ,NY,NX)+RUPH2B(2,11,NZ,NY,NX))/AREA(3,11,NY,NX)
                      IF(K.EQ.64)HEAD(M)=(RUPH2P(1,12,NZ,NY,NX)+RUPH2P(2,12,NZ,NY,NX) &
                        +RUPH2B(1,12,NZ,NY,NX)+RUPH2B(2,12,NZ,NY,NX))/AREA(3,12,NY,NX)
                      IF(K.EQ.65)HEAD(M)=(RUPH2P(1,13,NZ,NY,NX)+RUPH2P(2,13,NZ,NY,NX) &
                        +RUPH2B(1,13,NZ,NY,NX)+RUPH2B(2,13,NZ,NY,NX))/AREA(3,13,NY,NX)
                      IF(K.EQ.66)HEAD(M)=(RUPH2P(1,14,NZ,NY,NX)+RUPH2P(2,14,NZ,NY,NX) &
                        +RUPH2B(1,14,NZ,NY,NX)+RUPH2B(2,14,NZ,NY,NX))/AREA(3,14,NY,NX)
                      IF(K.EQ.67)HEAD(M)=(RUPH2P(1,15,NZ,NY,NX)+RUPH2P(2,15,NZ,NY,NX) &
                        +RUPH2B(1,15,NZ,NY,NX)+RUPH2B(2,15,NZ,NY,NX))/AREA(3,15,NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
!
!     HOURLY PLANT HEAT OUTPUT
!
              IF(N.EQ.25)THEN
                IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
                  DO K=51,100
                    IF(CHOICE(K,N-20).EQ.'YES')THEN
                      M=M+1
                      IF(K.EQ.51)HEAD(M)=277.8*RAD1(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.52)HEAD(M)=277.8*EFLXC(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.53)HEAD(M)=277.8*SFLXC(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.54)HEAD(M)=277.8*HFLXC(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                      IF(K.EQ.55)HEAD(M)=TCC(NZ,NY,NX)
                      IF(K.EQ.56)HEAD(M)=TFN3(NZ,NY,NX)
                    ENDIF
                  ENDDO
                  WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILP(N-20,NZ,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  RETURN
  END SUBROUTINE outph
!------------------------------------------------------------------------


  SUBROUTINE outsd(I,NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES DAILY OUTPUT FOR SOIL
!     C, N, P, WATER AND HEAT TO OUTPUT FILES DEPENDING
!     ON SELECTIONS MADE IN OUTPUT CONTROL FILES IN
!     THE RUN SCRIPT
!

  implicit none
  integer, intent(in) :: I,NE,NEX,NHW,NHE,NVN,NVS


  real(r8) :: HEAD(50)
  integer :: K,LUN,M,N,NX,NY
!     execution begins here
! second half of soil vars output
  DO N=26,30
    IF(DATAC(N,NE,NEX).NE.'NO')THEN
      LUN=N+10
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          !
          !     WRITE DAILY SOIL DATA TO OUTPUT FILES
          !
          DO M=1,50
            HEAD(M)=0.0
          ENDDO

          M=0
          !
          !     DAILY SOIL C OUTPUT
!
          IF(N.EQ.26)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=URSDC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.2)HEAD(M)=UORGC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.3)HEAD(M)=UORGF(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.4)HEAD(M)=UXCSN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.5)HEAD(M)=UCO2G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.6)HEAD(M)=UOXYG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.7)HEAD(M)=UCOP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.8)HEAD(M)=TOMT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.9)HEAD(M)=ORGC(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.10)HEAD(M)=UCH4G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.11)HEAD(M)=UDOCQ(NY,NX)/TAREA
                  IF(K.EQ.12)HEAD(M)=UDOCD(NY,NX)/TAREA
                  IF(K.EQ.13)HEAD(M)=UDICQ(NY,NX)/TAREA
                  IF(K.EQ.14)HEAD(M)=UDICD(NY,NX)/TAREA
                  IF(K.EQ.15)HEAD(M)=CO2E(NY,NX)
                  IF(K.EQ.16)HEAD(M)=TNBP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.17)HEAD(M)=UCO2F(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.18)HEAD(M)=ORGC(1,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.19)HEAD(M)=ORGC(2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.20)HEAD(M)=ORGC(3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.21)HEAD(M)=ORGC(4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.22)HEAD(M)=ORGC(5,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.23)HEAD(M)=ORGC(6,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.24)HEAD(M)=ORGC(7,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.25)HEAD(M)=ORGC(8,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.26)HEAD(M)=ORGC(9,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.27.AND.JZ>=10)HEAD(M)=ORGC(10,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.28.AND.JZ>=10)HEAD(M)=ORGC(11,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.29.AND.JZ>=10)HEAD(M)=ORGC(12,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.30.AND.JZ>=10)HEAD(M)=ORGC(13,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.31.AND.JZ>=10)HEAD(M)=ORGC(14,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.32.AND.JZ>=10)HEAD(M)=ORGC(15,NY,NX)/AREA(3,NU(NY,NX),NY,NX)

                  IF(K.EQ.41)HEAD(M)=UH2GG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.42)HEAD(M)=XHVSTE(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.43)HEAD(M)=ARLFC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.44)HEAD(M)=TGPP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.45)HEAD(M)=TRAU(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.46)HEAD(M)=TNPP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.47)HEAD(M)=THRE(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.48)HEAD(M)=UCH4F(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.49)HEAD(M)=UCO2S(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.50)HEAD(M)=WTSTGT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
    !     DAILY SOIL WATER OUTPUT
    !
          IF(N.EQ.27)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=1000.0_r8*URAIN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.2)HEAD(M)=1000.0_r8*UEVAP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.3)HEAD(M)=1000.0_r8*URUN(NY,NX)/TAREA
                  IF(K.EQ.4)HEAD(M)=1000.0_r8*UVOLW(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.5)HEAD(M)=1000.0_r8*UVOLO(NY,NX)/TAREA
                  IF(K.EQ.6)HEAD(M)=1000.0_r8*DPTHS(NY,NX)
                  IF(K.EQ.7)HEAD(M)=THETWZ(1,NY,NX)
                  IF(K.EQ.8)HEAD(M)=THETWZ(2,NY,NX)
                  IF(K.EQ.9)HEAD(M)=THETWZ(3,NY,NX)
                  IF(K.EQ.10)HEAD(M)=THETWZ(4,NY,NX)
                  IF(K.EQ.11)HEAD(M)=THETWZ(5,NY,NX)
                  IF(K.EQ.12)HEAD(M)=THETWZ(6,NY,NX)
                  IF(K.EQ.13)HEAD(M)=THETWZ(7,NY,NX)
                  IF(K.EQ.14)HEAD(M)=THETWZ(8,NY,NX)
                  IF(K.EQ.15)HEAD(M)=THETWZ(9,NY,NX)
                  IF(K.EQ.16)HEAD(M)=THETWZ(10,NY,NX)
                  IF(K.EQ.17)HEAD(M)=THETWZ(11,NY,NX)
                  IF(K.EQ.18)HEAD(M)=THETWZ(12,NY,NX)
                  IF(K.EQ.19)HEAD(M)=THETWZ(13,NY,NX)
                  IF(K.EQ.20)HEAD(M)=THETWZ(0,NY,NX)
                  IF(K.EQ.21)HEAD(M)=THETIZ(1,NY,NX)
                  IF(K.EQ.22)HEAD(M)=THETIZ(2,NY,NX)
                  IF(K.EQ.23)HEAD(M)=THETIZ(3,NY,NX)
                  IF(K.EQ.24)HEAD(M)=THETIZ(4,NY,NX)
                  IF(K.EQ.25)HEAD(M)=THETIZ(5,NY,NX)
                  IF(K.EQ.26)HEAD(M)=THETIZ(6,NY,NX)
                  IF(K.EQ.27)HEAD(M)=THETIZ(7,NY,NX)
                  IF(K.EQ.28)HEAD(M)=THETIZ(8,NY,NX)
                  IF(K.EQ.29)HEAD(M)=THETIZ(9,NY,NX)
                  IF(K.EQ.30)HEAD(M)=THETIZ(10,NY,NX)
                  IF(K.EQ.31)HEAD(M)=THETIZ(11,NY,NX)
                  IF(K.EQ.32)HEAD(M)=THETIZ(12,NY,NX)
                  IF(K.EQ.33)HEAD(M)=THETIZ(13,NY,NX)
                  IF(K.EQ.34)HEAD(M)=THETIZ(0,NY,NX)
                  IF(K.EQ.35)HEAD(M)=PSISM(1,NY,NX)+PSISO(1,NY,NX)
                  IF(K.EQ.36)HEAD(M)=PSISM(2,NY,NX)+PSISO(2,NY,NX)
                  IF(K.EQ.37)HEAD(M)=PSISM(3,NY,NX)+PSISO(3,NY,NX)
                  IF(K.EQ.38)HEAD(M)=PSISM(4,NY,NX)+PSISO(4,NY,NX)
                  IF(K.EQ.39)HEAD(M)=PSISM(5,NY,NX)+PSISO(5,NY,NX)
                  IF(K.EQ.40)HEAD(M)=PSISM(6,NY,NX)+PSISO(6,NY,NX)
                  IF(K.EQ.41)HEAD(M)=PSISM(7,NY,NX)+PSISO(7,NY,NX)
                  IF(K.EQ.42)HEAD(M)=PSISM(8,NY,NX)+PSISO(8,NY,NX)
                  IF(K.EQ.43)HEAD(M)=PSISM(9,NY,NX)+PSISO(9,NY,NX)
                  IF(K.EQ.44.AND.JZ>=10)HEAD(M)=PSISM(10,NY,NX)+PSISO(10,NY,NX)
                  IF(K.EQ.45.AND.JZ>=11)HEAD(M)=PSISM(11,NY,NX)+PSISO(11,NY,NX)
                  IF(K.EQ.46)HEAD(M)=USEDOU(NY,NX)*1000.0/TAREA
                  IF(K.EQ.47)HEAD(M)=PSISM(0,NY,NX)
                  IF(K.EQ.48)HEAD(M)=-CDPTH(NU(NY,NX)-1,NY,NX)+DLYR(3,0,NY,NX)
                  IF(K.EQ.49)HEAD(M)=-(DPTHA(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
                  IF(K.EQ.50)HEAD(M)=-(DPTHT(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
    !     DAILY SOIL N OUTPUT
    !
          IF(N.EQ.28)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=URSDN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.2)HEAD(M)=UORGN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.3)HEAD(M)=UFERTN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.4)HEAD(M)=UXZSN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.5)HEAD(M)=UNH4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.6)HEAD(M)=UNO3(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.7)HEAD(M)=UDONQ(NY,NX)/TAREA
                  IF(K.EQ.8)HEAD(M)=UDOND(NY,NX)/TAREA
                  IF(K.EQ.9)HEAD(M)=UDINQ(NY,NX)/TAREA
                  IF(K.EQ.10)HEAD(M)=UDIND(NY,NX)/TAREA
                  IF(K.EQ.11)HEAD(M)=UN2OG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.12)HEAD(M)=UNH3G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.13)HEAD(M)=UN2GS(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.14)HEAD(M)=TONT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.15)THEN
                    IF(BKVL(1,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,1,NY,NX)+trc_solml(ids_NH4B,1,NY,NX)+14.0*(trcx_solml(idx_NH4,1,NY,NX)+trcx_solml(idx_NH4B,1,NY,NX)))/BKVL(1,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.16)THEN
                    IF(BKVL(2,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,2,NY,NX)+trc_solml(ids_NH4B,2,NY,NX)+14.0*(trcx_solml(idx_NH4,2,NY,NX)+trcx_solml(idx_NH4B,2,NY,NX)))/BKVL(2,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.17)THEN
                    IF(BKVL(3,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,3,NY,NX)+trc_solml(ids_NH4B,3,NY,NX)+14.0*(trcx_solml(idx_NH4,3,NY,NX)+trcx_solml(idx_NH4B,3,NY,NX)))/BKVL(3,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.18)THEN
                    IF(BKVL(4,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,4,NY,NX)+trc_solml(ids_NH4B,4,NY,NX)+14.0*(trcx_solml(idx_NH4,4,NY,NX)+trcx_solml(idx_NH4B,4,NY,NX)))/BKVL(4,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.19)THEN
                    IF(BKVL(5,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,5,NY,NX)+trc_solml(ids_NH4B,5,NY,NX)+14.0*(trcx_solml(idx_NH4,5,NY,NX)+trcx_solml(idx_NH4B,5,NY,NX)))/BKVL(5,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.20)THEN
                    IF(BKVL(6,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,6,NY,NX)+trc_solml(ids_NH4B,6,NY,NX)+14.0*(trcx_solml(idx_NH4,6,NY,NX)+trcx_solml(idx_NH4B,6,NY,NX)))/BKVL(6,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.21)THEN
                    IF(BKVL(7,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,7,NY,NX)+trc_solml(ids_NH4B,7,NY,NX)+14.0*(trcx_solml(idx_NH4,7,NY,NX)+trcx_solml(idx_NH4B,7,NY,NX)))/BKVL(7,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.22)THEN
                    IF(BKVL(8,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,8,NY,NX)+trc_solml(ids_NH4B,8,NY,NX)+14.0*(trcx_solml(idx_NH4,8,NY,NX)+trcx_solml(idx_NH4B,8,NY,NX)))/BKVL(8,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.23)THEN
                    IF(BKVL(9,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,9,NY,NX)+trc_solml(ids_NH4B,9,NY,NX)+14.0*(trcx_solml(idx_NH4,9,NY,NX)+trcx_solml(idx_NH4B,9,NY,NX)))/BKVL(9,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.24.AND.JZ>=10)THEN
                    IF(BKVL(10,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,10,NY,NX)+trc_solml(ids_NH4B,10,NY,NX)+14.0*(trcx_solml(idx_NH4,10,NY,NX)+trcx_solml(idx_NH4B,10,NY,NX)))/BKVL(10,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.25.AND.JZ>=11)THEN
                    IF(BKVL(11,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,11,NY,NX)+trc_solml(ids_NH4B,11,NY,NX)+14.0*(trcx_solml(idx_NH4,11,NY,NX)+trcx_solml(idx_NH4B,11,NY,NX)))/BKVL(11,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.26.AND.JZ>=12)THEN
                    IF(BKVL(12,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,12,NY,NX)+trc_solml(ids_NH4B,12,NY,NX)+14.0*(trcx_solml(idx_NH4,12,NY,NX)+trcx_solml(idx_NH4B,12,NY,NX)))/BKVL(12,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.27.AND.JZ>=13)THEN
                    IF(BKVL(13,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,13,NY,NX)+trc_solml(ids_NH4B,13,NY,NX)+14.0*(trcx_solml(idx_NH4,13,NY,NX)+trcx_solml(idx_NH4B,13,NY,NX)))/BKVL(13,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.28.AND.JZ>=14)THEN
                    IF(BKVL(14,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,14,NY,NX)+trc_solml(ids_NH4B,14,NY,NX)+14.0*(trcx_solml(idx_NH4,14,NY,NX)+trcx_solml(idx_NH4B,14,NY,NX)))/BKVL(14,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.29.AND.JZ>=15)THEN
                    IF(BKVL(15,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,15,NY,NX)+trc_solml(ids_NH4B,15,NY,NX)+14.0*(trcx_solml(idx_NH4,15,NY,NX)+trcx_solml(idx_NH4B,15,NY,NX)))/BKVL(15,NY,NX)
                    ENDIF
                  ENDIF

                  IF(K.EQ.30)THEN
                    IF(BKVL(1,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,1,NY,NX)+trc_solml(ids_NO3B,1,NY,NX)+trc_solml(ids_NO2,1,NY,NX)+trc_solml(ids_NO2B,1,NY,NX))/BKVL(1,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.31)THEN
                    IF(BKVL(2,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,2,NY,NX)+trc_solml(ids_NO3B,2,NY,NX)+trc_solml(ids_NO2,2,NY,NX)+trc_solml(ids_NO2B,2,NY,NX))/BKVL(2,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.32)THEN
                    IF(BKVL(3,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,3,NY,NX)+trc_solml(ids_NO3B,3,NY,NX)+trc_solml(ids_NO2,3,NY,NX)+trc_solml(ids_NO2B,3,NY,NX))/BKVL(3,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.33)THEN
                    IF(BKVL(4,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,4,NY,NX)+trc_solml(ids_NO3B,4,NY,NX)+trc_solml(ids_NO2,4,NY,NX)+trc_solml(ids_NO2B,4,NY,NX))/BKVL(4,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.34)THEN
                    IF(BKVL(5,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,5,NY,NX)+trc_solml(ids_NO3B,5,NY,NX)+trc_solml(ids_NO2,5,NY,NX)+trc_solml(ids_NO2B,5,NY,NX))/BKVL(5,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.35)THEN
                    IF(BKVL(6,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,6,NY,NX)+trc_solml(ids_NO3B,6,NY,NX)+trc_solml(ids_NO2,6,NY,NX)+trc_solml(ids_NO2B,6,NY,NX))/BKVL(6,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.36)THEN
                    IF(BKVL(7,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,7,NY,NX)+trc_solml(ids_NO3B,7,NY,NX)+trc_solml(ids_NO2,7,NY,NX)+trc_solml(ids_NO2B,7,NY,NX))/BKVL(7,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.37)THEN
                    IF(BKVL(8,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,8,NY,NX)+trc_solml(ids_NO3B,8,NY,NX)+trc_solml(ids_NO2,8,NY,NX)+trc_solml(ids_NO2B,8,NY,NX))/BKVL(8,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.38)THEN
                    IF(BKVL(9,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,9,NY,NX)+trc_solml(ids_NO3B,9,NY,NX)+trc_solml(ids_NO2,9,NY,NX)+trc_solml(ids_NO2B,9,NY,NX))/BKVL(9,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.39.AND.JZ>=10)THEN
                    IF(BKVL(10,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,10,NY,NX)+trc_solml(ids_NO3B,10,NY,NX)+trc_solml(ids_NO2,10,NY,NX)+trc_solml(ids_NO2B,10,NY,NX))/BKVL(10,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.40.AND.JZ>=11)THEN
                    IF(BKVL(11,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,11,NY,NX)+trc_solml(ids_NO3B,11,NY,NX)+trc_solml(ids_NO2,11,NY,NX)+trc_solml(ids_NO2B,11,NY,NX))/BKVL(11,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.41.AND.JZ>=12)THEN
                    IF(BKVL(12,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,12,NY,NX)+trc_solml(ids_NO3B,12,NY,NX)+trc_solml(ids_NO2,12,NY,NX)+trc_solml(ids_NO2B,12,NY,NX))/BKVL(12,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.42.AND.JZ>=13)THEN
                    IF(BKVL(13,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,13,NY,NX)+trc_solml(ids_NO3B,13,NY,NX)+trc_solml(ids_NO2,13,NY,NX)+trc_solml(ids_NO2B,13,NY,NX))/BKVL(13,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.43.AND.JZ>=14)THEN
                    IF(BKVL(14,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,14,NY,NX)+trc_solml(ids_NO3B,14,NY,NX)+trc_solml(ids_NO2,14,NY,NX)+trc_solml(ids_NO2B,14,NY,NX))/BKVL(14,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.44.AND.JZ>=15)THEN
                    IF(BKVL(15,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,15,NY,NX)+trc_solml(ids_NO3B,15,NY,NX)+trc_solml(ids_NO2,15,NY,NX)+trc_solml(ids_NO2B,15,NY,NX))/BKVL(15,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.45)THEN
                    IF(BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NH4,0,NY,NX)+14.0*trcx_solml(idx_NH4,0,NY,NX))/BKVL(0,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.46)THEN
                    IF(BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_NO3,0,NY,NX)+trc_solml(ids_NO2,0,NY,NX))/BKVL(0,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.47)HEAD(M)=XHVSTE(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.48)HEAD(M)=-TRINH4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.49)HEAD(M)=UNH3F(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.50)HEAD(M)=UN2GG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
!     DAILY SOIL P OUTPUT
!
          IF(N.EQ.29)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=URSDP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.2)HEAD(M)=UORGP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.3)HEAD(M)=UFERTP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.4)HEAD(M)=UXPSN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.5)HEAD(M)=UPO4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.6)HEAD(M)=UDOPQ(NY,NX)/TAREA
                  IF(K.EQ.7)HEAD(M)=UDOPD(NY,NX)/TAREA
                  IF(K.EQ.8)HEAD(M)=UDIPQ(NY,NX)/TAREA
                  IF(K.EQ.9)HEAD(M)=UDIPD(NY,NX)/TAREA
                  IF(K.EQ.10)HEAD(M)=UPP4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.11)HEAD(M)=TOPT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.12)HEAD(M)=UPO4F(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.13)THEN
                    IF(BKVL(1,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,1,NY,NX)+trc_solml(ids_H1PO4B,1,NY,NX)+trc_solml(ids_H2PO4,1,NY,NX)+trc_solml(ids_H2PO4B,1,NY,NX))/VOLW(1,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.14)THEN
                    IF(BKVL(2,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,2,NY,NX)+trc_solml(ids_H2PO4B,1,NY,NX)+trc_solml(ids_H2PO4,2,NY,NX)+trc_solml(ids_H2PO4B,2,NY,NX))/VOLW(2,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.15)THEN
                    IF(BKVL(3,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,3,NY,NX)+trc_solml(ids_H1PO4B,3,NY,NX)+trc_solml(ids_H2PO4,3,NY,NX)+trc_solml(ids_H2PO4B,3,NY,NX))/VOLW(3,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.16)THEN
                    IF(BKVL(4,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,4,NY,NX)+trc_solml(ids_H1PO4B,4,NY,NX)+trc_solml(ids_H2PO4,4,NY,NX)+trc_solml(ids_H2PO4B,4,NY,NX))/VOLW(4,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.17)THEN
                    IF(BKVL(5,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,5,NY,NX)+trc_solml(ids_H1PO4B,5,NY,NX)+trc_solml(ids_H2PO4,5,NY,NX)+trc_solml(ids_H2PO4B,5,NY,NX))/VOLW(5,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.18)THEN
                    IF(BKVL(6,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,6,NY,NX)+trc_solml(ids_H1PO4B,6,NY,NX)+trc_solml(ids_H2PO4,6,NY,NX)+trc_solml(ids_H2PO4B,6,NY,NX))/VOLW(6,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.19)THEN
                    IF(BKVL(7,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,7,NY,NX)+trc_solml(ids_H1PO4B,7,NY,NX)+trc_solml(ids_H2PO4,7,NY,NX)+trc_solml(ids_H2PO4B,7,NY,NX))/VOLW(7,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.20)THEN
                    IF(BKVL(8,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,8,NY,NX)+trc_solml(ids_H1PO4B,8,NY,NX)+trc_solml(ids_H2PO4,8,NY,NX)+trc_solml(ids_H2PO4B,8,NY,NX))/VOLW(8,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.21)THEN
                    IF(BKVL(9,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,9,NY,NX)+trc_solml(ids_H1PO4B,9,NY,NX)+trc_solml(ids_H2PO4,9,NY,NX)+trc_solml(ids_H2PO4B,9,NY,NX))/VOLW(9,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.22.AND.JZ>=10)THEN
                    IF(BKVL(10,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,10,NY,NX)+trc_solml(ids_H1PO4B,10,NY,NX)+trc_solml(ids_H2PO4,10,NY,NX)+trc_solml(ids_H2PO4B,10,NY,NX))/VOLW(10,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.23.AND.JZ>=11)THEN
                    IF(BKVL(11,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,11,NY,NX)+trc_solml(ids_H1PO4B,11,NY,NX)+trc_solml(ids_H2PO4,11,NY,NX)+trc_solml(ids_H2PO4B,11,NY,NX))/VOLW(11,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.24.AND.JZ>=12)THEN
                    IF(BKVL(12,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,12,NY,NX)+trc_solml(ids_H1PO4B,1,NY,NX)+trc_solml(ids_H2PO4,12,NY,NX)+trc_solml(ids_H2PO4B,12,NY,NX))/VOLW(12,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.25.AND.JZ>=13)THEN
                    IF(BKVL(13,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,13,NY,NX)+trc_solml(ids_H1PO4B,13,NY,NX)+trc_solml(ids_H2PO4,13,NY,NX)+trc_solml(ids_H2PO4B,13,NY,NX))/VOLW(13,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.26.AND.JZ>=14)THEN
                    IF(BKVL(14,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,14,NY,NX)+trc_solml(ids_H1PO4B,14,NY,NX)+trc_solml(ids_H2PO4,14,NY,NX)+trc_solml(ids_H2PO4B,14,NY,NX))/VOLW(14,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.27.AND.JZ>=15)THEN
                    IF(BKVL(15,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=(trc_solml(ids_H1PO4,15,NY,NX)+trc_solml(ids_H1PO4B,15,NY,NX)+trc_solml(ids_H2PO4,15,NY,NX)+trc_solml(ids_H2PO4B,15,NY,NX))/VOLW(15,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.28)THEN
                    IF(BKVL(1,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,1,NY,NX)+trcx_solml(idx_H2PO4,1,NY,NX)+trcx_solml(idx_HPO4B,1,NY,NX)+trcx_solml(idx_H2PO4B,1,NY,NX))/BKVL(1,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.29)THEN
                    IF(BKVL(2,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,2,NY,NX)+trcx_solml(idx_H2PO4,2,NY,NX)+trcx_solml(idx_HPO4B,2,NY,NX)+trcx_solml(idx_H2PO4B,2,NY,NX))/BKVL(2,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.30)THEN
                    IF(BKVL(3,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,3,NY,NX)+trcx_solml(idx_H2PO4,3,NY,NX)+trcx_solml(idx_HPO4B,3,NY,NX)+trcx_solml(idx_H2PO4B,3,NY,NX))/BKVL(3,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.31)THEN
                    IF(BKVL(4,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,4,NY,NX)+trcx_solml(idx_H2PO4,4,NY,NX)+trcx_solml(idx_HPO4B,4,NY,NX)+trcx_solml(idx_H2PO4B,4,NY,NX))/BKVL(4,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.32)THEN
                    IF(BKVL(5,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,5,NY,NX)+trcx_solml(idx_H2PO4,5,NY,NX)+trcx_solml(idx_HPO4B,5,NY,NX)+trcx_solml(idx_H2PO4B,5,NY,NX))/BKVL(5,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.33)THEN
                    IF(BKVL(6,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,6,NY,NX)+trcx_solml(idx_H2PO4,6,NY,NX)+trcx_solml(idx_HPO4B,6,NY,NX)+trcx_solml(idx_H2PO4B,6,NY,NX))/BKVL(6,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.34)THEN
                    IF(BKVL(7,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,7,NY,NX)+trcx_solml(idx_H2PO4,7,NY,NX)+trcx_solml(idx_HPO4B,7,NY,NX)+trcx_solml(idx_H2PO4B,7,NY,NX))/BKVL(7,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.35)THEN
                    IF(BKVL(8,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,8,NY,NX)+trcx_solml(idx_H2PO4,8,NY,NX)+trcx_solml(idx_HPO4B,8,NY,NX)+trcx_solml(idx_H2PO4B,8,NY,NX))/BKVL(8,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.36)THEN
                    IF(BKVL(9,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,9,NY,NX)+trcx_solml(idx_H2PO4,9,NY,NX)+trcx_solml(idx_HPO4B,9,NY,NX)+trcx_solml(idx_H2PO4B,9,NY,NX))/BKVL(9,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.37.AND.JZ>=10)THEN
                    IF(BKVL(10,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,10,NY,NX)+trcx_solml(idx_H2PO4,10,NY,NX)+trcx_solml(idx_HPO4B,10,NY,NX)+trcx_solml(idx_H2PO4B,10,NY,NX))/BKVL(10,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.38.AND.JZ>=11)THEN
                    IF(BKVL(11,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,11,NY,NX)+trcx_solml(idx_H2PO4,11,NY,NX)+trcx_solml(idx_HPO4B,11,NY,NX)+trcx_solml(idx_H2PO4B,11,NY,NX))/BKVL(11,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.39.AND.JZ>=12)THEN
                    IF(BKVL(12,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,12,NY,NX)+trcx_solml(idx_H2PO4,12,NY,NX)+trcx_solml(idx_HPO4B,12,NY,NX)+trcx_solml(idx_H2PO4B,12,NY,NX))/BKVL(12,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.40.AND.JZ>=13)THEN
                    IF(BKVL(13,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,13,NY,NX)+trcx_solml(idx_H2PO4,13,NY,NX)+trcx_solml(idx_HPO4B,13,NY,NX)+trcx_solml(idx_H2PO4B,13,NY,NX))/BKVL(13,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.41.AND.JZ>=14)THEN
                    IF(BKVL(14,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,14,NY,NX)+trcx_solml(idx_H2PO4,14,NY,NX)+trcx_solml(idx_HPO4B,14,NY,NX)+trcx_solml(idx_H2PO4B,14,NY,NX))/BKVL(14,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.42.AND.JZ>=15)THEN
                    IF(BKVL(15,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,15,NY,NX)+trcx_solml(idx_H2PO4,15,NY,NX)+trcx_solml(idx_HPO4B,15,NY,NX)+trcx_solml(idx_H2PO4B,15,NY,NX))/BKVL(15,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.43)THEN
                    IF(BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=trc_solml(ids_H2PO4,0,NY,NX)/BKVL(0,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.44)THEN
                    IF(BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
                      HEAD(M)=31.0*(trcx_solml(idx_HPO4,0,NY,NX)+trcx_solml(idx_H2PO4,0,NY,NX))/BKVL(0,NY,NX)
                    ENDIF
                  ENDIF
                  IF(K.EQ.47)HEAD(M)=XHVSTE(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.48)HEAD(M)=-TRIPO4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
!     DAILY SOIL HEAT OUTPUT
!
          IF(N.EQ.30)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=TRAD(NY,NX)
                  IF(K.EQ.2)HEAD(M)=TAMX(NY,NX)
                  IF(K.EQ.3)HEAD(M)=TAMN(NY,NX)
                  IF(K.EQ.4)HEAD(M)=HUDX(NY,NX)
                  IF(K.EQ.5)HEAD(M)=HUDN(NY,NX)
                  IF(K.EQ.6)HEAD(M)=TWIND(NY,NX)*0.001
                  IF(K.EQ.7)HEAD(M)=TRAI(NY,NX)
                  IF(K.EQ.8)HEAD(M)=TSMX(1,NY,NX)
                  IF(K.EQ.9)HEAD(M)=TSMN(1,NY,NX)
                  IF(K.EQ.10)HEAD(M)=TSMX(2,NY,NX)
                  IF(K.EQ.11)HEAD(M)=TSMN(2,NY,NX)
                  IF(K.EQ.12)HEAD(M)=TSMX(3,NY,NX)
                  IF(K.EQ.13)HEAD(M)=TSMN(3,NY,NX)
                  IF(K.EQ.14)HEAD(M)=TSMX(4,NY,NX)
                  IF(K.EQ.15)HEAD(M)=TSMN(4,NY,NX)
                  IF(K.EQ.16)HEAD(M)=TSMX(5,NY,NX)
                  IF(K.EQ.17)HEAD(M)=TSMN(5,NY,NX)
                  IF(K.EQ.18)HEAD(M)=TSMX(6,NY,NX)
                  IF(K.EQ.19)HEAD(M)=TSMN(6,NY,NX)
                  IF(K.EQ.20)HEAD(M)=TSMX(7,NY,NX)
                  IF(K.EQ.21)HEAD(M)=TSMN(7,NY,NX)
                  IF(K.EQ.22)HEAD(M)=TSMX(8,NY,NX)
                  IF(K.EQ.23)HEAD(M)=TSMN(8,NY,NX)
                  IF(K.EQ.24)HEAD(M)=TSMX(9,NY,NX)
                  IF(K.EQ.25)HEAD(M)=TSMN(9,NY,NX)
                  IF(K.EQ.26.AND.JZ>=10)HEAD(M)=TSMX(10,NY,NX)
                  IF(K.EQ.27.AND.JZ>=10)HEAD(M)=TSMN(10,NY,NX)
                  IF(K.EQ.28.AND.JZ>=11)HEAD(M)=TSMX(11,NY,NX)
                  IF(K.EQ.29.AND.JZ>=11)HEAD(M)=TSMN(11,NY,NX)
                  IF(K.EQ.30.AND.JZ>=12)HEAD(M)=TSMX(12,NY,NX)
                  IF(K.EQ.31.AND.JZ>=12)HEAD(M)=TSMN(12,NY,NX)
                  IF(K.EQ.32.AND.JZ>=13)HEAD(M)=TSMX(13,NY,NX)
                  IF(K.EQ.33.AND.JZ>=13)HEAD(M)=TSMN(13,NY,NX)
                  IF(K.EQ.34.AND.JZ>=14)HEAD(M)=TSMX(14,NY,NX)
                  IF(K.EQ.35.AND.JZ>=14)HEAD(M)=TSMN(14,NY,NX)
                  IF(K.EQ.36)HEAD(M)=TSMX(0,NY,NX)
                  IF(K.EQ.37)HEAD(M)=TSMN(0,NY,NX)
                  IF(K.EQ.38)HEAD(M)=ECND(1,NY,NX)
                  IF(K.EQ.39)HEAD(M)=ECND(2,NY,NX)
                  IF(K.EQ.40)HEAD(M)=ECND(3,NY,NX)
                  IF(K.EQ.41)HEAD(M)=ECND(4,NY,NX)
                  IF(K.EQ.42)HEAD(M)=ECND(5,NY,NX)
                  IF(K.EQ.43)HEAD(M)=ECND(6,NY,NX)
                  IF(K.EQ.44)HEAD(M)=ECND(7,NY,NX)
                  IF(K.EQ.45)HEAD(M)=ECND(8,NY,NX)
                  IF(K.EQ.46)HEAD(M)=ECND(9,NY,NX)
                  IF(K.EQ.47.AND.JZ>=10)HEAD(M)=ECND(10,NY,NX)
                  IF(K.EQ.48.AND.JZ>=11)HEAD(M)=ECND(11,NY,NX)
                  IF(K.EQ.49.AND.JZ>=12)HEAD(M)=ECND(12,NY,NX)
            !     IF(K.EQ.50.AND.JZ>=13)HEAD(M)=ECND(13,NY,NX)
                  IF(K.EQ.50)HEAD(M)=UIONOU(NY,NX)/TAREA
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  RETURN
  END subroutine outsd


!------------------------------------------------------------------------

  SUBROUTINE outsh(I,J,NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES HOURLY OUTPUT FOR SOIL
!     C, N, P, WATER AND HEAT TO OUTPUT FILES DEPENDING
!     ON SELECTIONS MADE IN OUTPUT CONTROL FILES IN
!     THE RUN SCRIPT
!

  implicit none
  integer, intent(in) :: I,J,NE
  integer, intent(in) :: NEX,NHW,NHE,NVN,NVS


  real(r8) :: HEAD(50)
  integer :: K,LUN,M,N,NX,NY
!     execution begins here
! first half of soil vars output
  DO N=21,25
    IF(DATAC(N,NE,NEX).NE.'NO')THEN
      DO  M=1,50
        HEAD(M)=0.0
      ENDDO
      LUN=N+10
      DO NX=NHW,NHE
        DO NY=NVN,NVS
!
!     WRITE HOURLY SOIL DATA TO OUTPUT FILES
!
          M=0
!
!     HOURLY SOIL C OUTPUT
!
          IF(N.EQ.21)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=HCO2G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
                  IF(K.EQ.2)HEAD(M)=TCNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
                  IF(K.EQ.3)HEAD(M)=HCH4G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
                  IF(K.EQ.4)HEAD(M)=HOXYG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*8.68056
                  IF(K.EQ.5)HEAD(M)=trc_solcl(idg_CO2,1,NY,NX)
                  IF(K.EQ.6)HEAD(M)=trc_solcl(idg_CO2,2,NY,NX)
                  IF(K.EQ.7)HEAD(M)=trc_solcl(idg_CO2,3,NY,NX)
                  IF(K.EQ.8)HEAD(M)=trc_solcl(idg_CO2,4,NY,NX)
                  IF(K.EQ.9)HEAD(M)=trc_solcl(idg_CO2,5,NY,NX)
                  IF(K.EQ.10)HEAD(M)=trc_solcl(idg_CO2,6,NY,NX)
                  IF(K.EQ.11)HEAD(M)=trc_solcl(idg_CO2,7,NY,NX)
                  IF(K.EQ.12)HEAD(M)=trc_solcl(idg_CO2,8,NY,NX)
                  IF(K.EQ.13)HEAD(M)=trc_solcl(idg_CO2,9,NY,NX)
                  IF(K.EQ.14.AND.JZ>=10)HEAD(M)=trc_solcl(idg_CO2,10,NY,NX)
                  IF(K.EQ.15.AND.JZ>=11)HEAD(M)=trc_solcl(idg_CO2,11,NY,NX)
                  IF(K.EQ.16.AND.JZ>=12)HEAD(M)=trc_solcl(idg_CO2,12,NY,NX)
                  IF(K.EQ.17.AND.JZ>=13)HEAD(M)=trc_solcl(idg_CO2,13,NY,NX)
                  IF(K.EQ.18.AND.JZ>=14)HEAD(M)=trc_solcl(idg_CO2,14,NY,NX)

                  IF(K.EQ.19)HEAD(M)=trc_solcl(idg_CO2,0,NY,NX)
                  IF(K.EQ.20)HEAD(M)=trc_solcl(idg_CH4,1,NY,NX)
                  IF(K.EQ.21)HEAD(M)=trc_solcl(idg_CH4,2,NY,NX)
                  IF(K.EQ.22)HEAD(M)=trc_solcl(idg_CH4,3,NY,NX)
                  IF(K.EQ.23)HEAD(M)=trc_solcl(idg_CH4,4,NY,NX)
                  IF(K.EQ.24)HEAD(M)=trc_solcl(idg_CH4,5,NY,NX)
                  IF(K.EQ.25)HEAD(M)=trc_solcl(idg_CH4,6,NY,NX)
                  IF(K.EQ.26)HEAD(M)=trc_solcl(idg_CH4,7,NY,NX)
                  IF(K.EQ.27)HEAD(M)=trc_solcl(idg_CH4,8,NY,NX)
                  IF(K.EQ.28)HEAD(M)=trc_solcl(idg_CH4,9,NY,NX)
                  IF(K.EQ.29.AND.JZ>=10)HEAD(M)=trc_solcl(idg_CH4,10,NY,NX)
                  IF(K.EQ.30.AND.JZ>=11)HEAD(M)=trc_solcl(idg_CH4,11,NY,NX)
                  IF(K.EQ.31.AND.JZ>=12)HEAD(M)=trc_solcl(idg_CH4,12,NY,NX)
                  IF(K.EQ.32.AND.JZ>=13)HEAD(M)=trc_solcl(idg_CH4,13,NY,NX)
                  IF(K.EQ.33.AND.JZ>=14)HEAD(M)=trc_solcl(idg_CH4,14,NY,NX)
                  IF(K.EQ.34.AND.JZ>=15)HEAD(M)=trc_solcl(idg_CH4,15,NY,NX)

                  IF(K.EQ.35)HEAD(M)=trc_solcl(idg_O2,1,NY,NX)
                  IF(K.EQ.36)HEAD(M)=trc_solcl(idg_O2,2,NY,NX)
                  IF(K.EQ.37)HEAD(M)=trc_solcl(idg_O2,3,NY,NX)
                  IF(K.EQ.38)HEAD(M)=trc_solcl(idg_O2,4,NY,NX)
                  IF(K.EQ.39)HEAD(M)=trc_solcl(idg_O2,5,NY,NX)
                  IF(K.EQ.40)HEAD(M)=trc_solcl(idg_O2,6,NY,NX)
                  IF(K.EQ.41)HEAD(M)=trc_solcl(idg_O2,7,NY,NX)
                  IF(K.EQ.42)HEAD(M)=trc_solcl(idg_O2,8,NY,NX)
                  IF(K.EQ.43)HEAD(M)=trc_solcl(idg_O2,9,NY,NX)
                  IF(K.EQ.44.AND.JZ>=10)HEAD(M)=trc_solcl(idg_O2,10,NY,NX)
                  IF(K.EQ.45.AND.JZ>=11)HEAD(M)=trc_solcl(idg_O2,11,NY,NX)
                  IF(K.EQ.46.AND.JZ>=12)HEAD(M)=trc_solcl(idg_O2,12,NY,NX)
                  IF(K.EQ.47.AND.JZ>=13)HEAD(M)=trc_solcl(idg_O2,13,NY,NX)
                  IF(K.EQ.48.AND.JZ>=14)HEAD(M)=trc_solcl(idg_O2,14,NY,NX)
                  IF(K.EQ.49.AND.JZ>=15)HEAD(M)=trc_solcl(idg_O2,15,NY,NX)
                  IF(K.EQ.50)HEAD(M)=trc_solcl(idg_O2,0,NY,NX)
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
!     HOURLY SOIL WATER OUTPUT
!
          IF(N.EQ.22)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=TEVAPG(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.2)HEAD(M)=-WQRH(NY,NX)*1000.0/TAREA
                  IF(K.EQ.3)HEAD(M)=USEDOU(NY,NX)*1000.0/TAREA
                  IF(K.EQ.4)HEAD(M)=UVOLW(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.5)HEAD(M)=HVOLO(NY,NX)*1000.0/TAREA
                  IF(K.EQ.6)HEAD(M)=AMAX1(1.0E-64,(VOLSS(NY,NX)+VOLIS(NY,NX) &
                    *DENSI+VOLWS(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX))

                  IF(K.EQ.7)HEAD(M)=THETWZ(1,NY,NX)
                  IF(K.EQ.8)HEAD(M)=THETWZ(2,NY,NX)
                  IF(K.EQ.9)HEAD(M)=THETWZ(3,NY,NX)
                  IF(K.EQ.10)HEAD(M)=THETWZ(4,NY,NX)
                  IF(K.EQ.11)HEAD(M)=THETWZ(5,NY,NX)
                  IF(K.EQ.12)HEAD(M)=THETWZ(6,NY,NX)
                  IF(K.EQ.13)HEAD(M)=THETWZ(7,NY,NX)
                  IF(K.EQ.14)HEAD(M)=THETWZ(8,NY,NX)
                  IF(K.EQ.15)HEAD(M)=THETWZ(9,NY,NX)
                  IF(K.EQ.16.AND.JZ>=10)HEAD(M)=THETWZ(10,NY,NX)
                  IF(K.EQ.17.AND.JZ>=11)HEAD(M)=THETWZ(11,NY,NX)
                  IF(K.EQ.18.AND.JZ>=12)HEAD(M)=THETWZ(12,NY,NX)
                  IF(K.EQ.19.AND.JZ>=13)HEAD(M)=THETWZ(13,NY,NX)
                  IF(K.EQ.20.AND.JZ>=14)HEAD(M)=THETWZ(14,NY,NX)
                  IF(K.EQ.21.AND.JZ>=15)HEAD(M)=THETWZ(15,NY,NX)
                  IF(K.EQ.22.AND.JZ>=16)HEAD(M)=THETWZ(16,NY,NX)
                  IF(K.EQ.23.AND.JZ>=17)HEAD(M)=THETWZ(17,NY,NX)
                  IF(K.EQ.24.AND.JZ>=18)HEAD(M)=THETWZ(18,NY,NX)
                  IF(K.EQ.25.AND.JZ>=19)HEAD(M)=THETWZ(19,NY,NX)
                  IF(K.EQ.26.AND.JZ>=20)HEAD(M)=THETWZ(20,NY,NX)

                  IF(K.EQ.27)HEAD(M)=THETWZ(0,NY,NX)
                  IF(K.EQ.28)HEAD(M)=THETIZ(1,NY,NX)
                  IF(K.EQ.29)HEAD(M)=THETIZ(2,NY,NX)
                  IF(K.EQ.30)HEAD(M)=THETIZ(3,NY,NX)
                  IF(K.EQ.31)HEAD(M)=THETIZ(4,NY,NX)
                  IF(K.EQ.32)HEAD(M)=THETIZ(5,NY,NX)
                  IF(K.EQ.33)HEAD(M)=THETIZ(6,NY,NX)
                  IF(K.EQ.34)HEAD(M)=THETIZ(7,NY,NX)
                  IF(K.EQ.35)HEAD(M)=THETIZ(8,NY,NX)
                  IF(K.EQ.36)HEAD(M)=THETIZ(9,NY,NX)
                  IF(K.EQ.37.AND.JZ>=10)HEAD(M)=THETIZ(10,NY,NX)
                  IF(K.EQ.38.AND.JZ>=11)HEAD(M)=THETIZ(11,NY,NX)
                  IF(K.EQ.39.AND.JZ>=12)HEAD(M)=THETIZ(12,NY,NX)
                  IF(K.EQ.40.AND.JZ>=13)HEAD(M)=THETIZ(13,NY,NX)
                  IF(K.EQ.41.AND.JZ>=14)HEAD(M)=THETIZ(14,NY,NX)
                  IF(K.EQ.42.AND.JZ>=15)HEAD(M)=THETIZ(15,NY,NX)
                  IF(K.EQ.43.AND.JZ>=16)HEAD(M)=THETIZ(16,NY,NX)
                  IF(K.EQ.44.AND.JZ>=17)HEAD(M)=THETIZ(17,NY,NX)
                  IF(K.EQ.45.AND.JZ>=18)HEAD(M)=THETIZ(18,NY,NX)
                  IF(K.EQ.46.AND.JZ>=19)HEAD(M)=THETIZ(19,NY,NX)
                  IF(K.EQ.47.AND.JZ>=20)HEAD(M)=THETIZ(20,NY,NX)

                  IF(K.EQ.48)HEAD(M)=THETIZ(0,NY,NX)
                  IF(K.EQ.49)HEAD(M)=-(DPTHA(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
                  IF(K.EQ.50)HEAD(M)=-(DPTHT(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
!     HOURLY SOIL N OUTPUT
!
          IF(N.EQ.23)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=HN2OG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.2)HEAD(M)=HN2GG(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.3)HEAD(M)=HNH3G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.4)HEAD(M)=UDINQ(NY,NX)/TAREA
                  IF(K.EQ.5)HEAD(M)=UDIND(NY,NX)/TAREA
                  IF(K.EQ.6)HEAD(M)=trc_solcl(idg_N2O,1,NY,NX)
                  IF(K.EQ.7)HEAD(M)=trc_solcl(idg_N2O,2,NY,NX)
                  IF(K.EQ.8)HEAD(M)=trc_solcl(idg_N2O,3,NY,NX)
                  IF(K.EQ.9)HEAD(M)=trc_solcl(idg_N2O,4,NY,NX)
                  IF(K.EQ.10)HEAD(M)=trc_solcl(idg_N2O,5,NY,NX)
                  IF(K.EQ.11)HEAD(M)=trc_solcl(idg_N2O,6,NY,NX)
                  IF(K.EQ.12)HEAD(M)=trc_solcl(idg_N2O,7,NY,NX)
                  IF(K.EQ.13)HEAD(M)=trc_solcl(idg_N2O,8,NY,NX)
                  IF(K.EQ.14)HEAD(M)=trc_solcl(idg_N2O,9,NY,NX)
                  IF(K.EQ.15.AND.JZ>=10)HEAD(M)=trc_solcl(idg_N2O,10,NY,NX)
                  IF(K.EQ.16.AND.JZ>=11)HEAD(M)=trc_solcl(idg_N2O,11,NY,NX)
                  IF(K.EQ.17.AND.JZ>=12)HEAD(M)=trc_solcl(idg_N2O,12,NY,NX)
                  IF(K.EQ.18.AND.JZ>=13)HEAD(M)=trc_solcl(idg_N2O,13,NY,NX)
                  IF(K.EQ.19.AND.JZ>=14)HEAD(M)=trc_solcl(idg_N2O,14,NY,NX)
                  IF(K.EQ.20.AND.JZ>=15)HEAD(M)=trc_solcl(idg_N2O,15,NY,NX)

                  IF(K.EQ.21)HEAD(M)=trc_solcl(idg_NH3,1,NY,NX)
                  IF(K.EQ.22)HEAD(M)=trc_solcl(idg_NH3,2,NY,NX)
                  IF(K.EQ.23)HEAD(M)=trc_solcl(idg_NH3,3,NY,NX)
                  IF(K.EQ.24)HEAD(M)=trc_solcl(idg_NH3,4,NY,NX)
                  IF(K.EQ.25)HEAD(M)=trc_solcl(idg_NH3,5,NY,NX)
                  IF(K.EQ.26)HEAD(M)=trc_solcl(idg_NH3,6,NY,NX)
                  IF(K.EQ.27)HEAD(M)=trc_solcl(idg_NH3,7,NY,NX)
                  IF(K.EQ.28)HEAD(M)=trc_solcl(idg_NH3,8,NY,NX)
                  IF(K.EQ.29)HEAD(M)=trc_solcl(idg_NH3,9,NY,NX)
                  IF(K.EQ.30.AND.JZ>=10)HEAD(M)=trc_solcl(idg_NH3,10,NY,NX)
                  IF(K.EQ.31.AND.JZ>=11)HEAD(M)=trc_solcl(idg_NH3,11,NY,NX)
                  IF(K.EQ.32.AND.JZ>=12)HEAD(M)=trc_solcl(idg_NH3,12,NY,NX)
                  IF(K.EQ.33.AND.JZ>=13)HEAD(M)=trc_solcl(idg_NH3,13,NY,NX)
                  IF(K.EQ.34.AND.JZ>=14)HEAD(M)=trc_solcl(idg_NH3,14,NY,NX)
                  IF(K.EQ.35.AND.JZ>=15)HEAD(M)=trc_solcl(idg_NH3,15,NY,NX)
                  IF(K.EQ.36)HEAD(M)=trc_solcl(idg_N2O,0,NY,NX)
                  IF(K.EQ.37)HEAD(M)=trc_solcl(idg_NH3,0,NY,NX)
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
!     HOURLY SOIL P OUTPUT
!
          IF(N.EQ.24)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=UDIPQ(NY,NX)/TAREA
                  IF(K.EQ.2)HEAD(M)=UDIPD(NY,NX)/TAREA
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
!
!     HOURLY SOIL HEAT OUTPUT
!
          IF(N.EQ.25)THEN
            IF(I.GE.IDATA(N).AND.I.LE.IDATA(N+20))THEN
              DO K=1,50
                IF(CHOICE(K,N-20).EQ.'YES')THEN
                  M=M+1
                  IF(K.EQ.1)HEAD(M)=RAD(NY,NX)*277.8
                  IF(K.EQ.2)HEAD(M)=TCA(NY,NX)
                  IF(K.EQ.3)HEAD(M)=VPK(NY,NX)
                  IF(K.EQ.4)HEAD(M)=UA(NY,NX)/3600.0
                  IF(K.EQ.5)HEAD(M)=(PRECR(NY,NX)+PRECW(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.6)HEAD(M)=HEATI(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.7)HEAD(M)=HEATE(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.8)HEAD(M)=HEATS(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.9)HEAD(M)=-(HEATH(NY,NX)-HEATV(NY,NX))*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.10)HEAD(M)=TRN(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.11)HEAD(M)=TLE(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.12)HEAD(M)=TSH(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.13)HEAD(M)=TGH(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
                  IF(K.EQ.14)HEAD(M)=TCS(1,NY,NX)
                  IF(K.EQ.15)HEAD(M)=TCS(2,NY,NX)
                  IF(K.EQ.16)HEAD(M)=TCS(3,NY,NX)
                  IF(K.EQ.17)HEAD(M)=TCS(4,NY,NX)
                  IF(K.EQ.18)HEAD(M)=TCS(5,NY,NX)
                  IF(K.EQ.19)HEAD(M)=TCS(6,NY,NX)
                  IF(K.EQ.20)HEAD(M)=TCS(7,NY,NX)
                  IF(K.EQ.21)HEAD(M)=TCS(8,NY,NX)
                  IF(K.EQ.22)HEAD(M)=TCS(9,NY,NX)
                  IF(K.EQ.23.AND.JZ>=10)HEAD(M)=TCS(10,NY,NX)
                  IF(K.EQ.24.AND.JZ>=11)HEAD(M)=TCS(11,NY,NX)
                  IF(K.EQ.25.AND.JZ>=12)HEAD(M)=TCS(12,NY,NX)
                  IF(K.EQ.26.AND.JZ>=13)HEAD(M)=TCS(13,NY,NX)
                  IF(K.EQ.27.AND.JZ>=14)HEAD(M)=TCS(14,NY,NX)
                  IF(K.EQ.28.AND.JZ>=15)HEAD(M)=TCS(15,NY,NX)
                  IF(K.EQ.29.AND.JZ>=16)HEAD(M)=TCS(16,NY,NX)
                  IF(K.EQ.30.AND.JZ>=17)HEAD(M)=TCS(17,NY,NX)
                  IF(K.EQ.31.AND.JZ>=18)HEAD(M)=TCS(18,NY,NX)
                  IF(K.EQ.32.AND.JZ>=19)HEAD(M)=TCS(19,NY,NX)
                  IF(K.EQ.33.AND.JZ>=20)HEAD(M)=TCS(20,NY,NX)
                  IF(K.EQ.34)HEAD(M)=TCS(0,NY,NX)
                  IF(K.EQ.35)HEAD(M)=TCW(1,NY,NX)
                ENDIF
              ENDDO
              WRITE(LUN,'(A16,F8.3,4X,A8,I8,50E16.7E3)')OUTFILS(N-20,NY,NX),DOY,CDATE,J,(HEAD(K),K=1,M)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  RETURN
  END SUBROUTINE outsh

end module Hist1Mod
