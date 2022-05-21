module RedistMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : padr, print_info,endrun,destroy
  use minimathmod, only : safe_adb
  use MicBGCPars, only : micpar
  use EcosimConst
  use FlagDataType
  use GridDataType
  use ErosionBalMod
  use RunoffBalMod
  use TFlxTypeMod
  use SoilLayerDynMod
  use SOMDataType
  USE SoilPropertyDataType
  use SoilWaterDataType
  USE SurfLitterDataType
  USE SoilHeatDataType
  USE EcoSIMCtrlDataType
  USE SoilBGCDataType
  USE EcosysBGCFluxType
  use PlantDataRateType
  USE CanopyDataType
  USE ClimForcDataType
  USE SurfSoilDataType
  USE SnowDataType
  USE LandSurfDataType
  USE FertilizerDataType
  USE AqueChemDatatype
  USE ChemTranspDataType
  USE IrrigationDataType
  USE RootDataType
  USE SoilPhysDataType
  use MicrobialDataType
  USE SedimentDataType
  USE EcoSimSumDataType
  USE LateralTranspMod
  use SnowBalMod
  implicit none

  private

  character(len=*), parameter :: mod_filename = __FILE__
  real(r8) :: CI,CH,CO,CX
  real(r8) :: CS,CIB,CHB,COB,CORP,CORP0,DCORPW
  real(r8) :: DORGP,DC,DN,DP,DDLYXS,DDLYXX,DDLYRS
  real(r8) :: DCORPZ
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA,ECCO,ECHC
  real(r8) :: ECSO,ECCL
  real(r8) :: ENGYZ,ENGYR,ENGYM
  real(r8) :: ENGYV,ENGYL
  real(r8) :: FI,HRAINR,HFLXO,HI,HO
  real(r8) :: HS,HGB,HOB,HFLXD,OI,OO,OS,OIB,OOB
  real(r8) :: OC,ON,OP,PSS
  real(r8) :: PI,PXB,POS,POX,POP,RAINR
  real(r8) :: SIN,SGN,SIP,SNB
  real(r8) :: SPB,SNM0,SPM0,SIR,SII,SBU,SSW,SSS,SNM,SPM,SSB
  real(r8) :: SD,SSH,SSF,SSX,SSP,SST,TFLWT,TVHCP
  real(r8) :: TVHCM,TVOLW,TVOLWH,TVOLI,TVOLIH,TENGY
  real(r8) :: TBKDX,TFC,TWP,TSCNV,TSCNH,TSAND
  real(r8) :: TSILT,TCLAY,TXCEC,TXAEC,TGKC4,TGKCA,TGKCM,TGKCN
  real(r8) :: TGKCK,TNFNIH,TNH4FA,TNH3FA,TNHUFA,TNO3FA,TNH4FB
  real(r8) :: TNH3FB,TNHUFB,TNO3FB,TNH4S,TNH4B,TNH3S,TNH3B,TNO3S
  real(r8) :: TNO3B,TNO2S,TNO2B,TZAL,TZFE,TZHY,TZCA,TZMG,TZNA
  real(r8) :: TZKA,TZOH,TZSO4,TZCL,TZCO3,TZHCO3,TZALO1,TZALO2
  real(r8) :: TZALO3,TZALO4,TZALS,TZFEO1,TZFEO2,TZFEO3,TZFEO4
  real(r8) :: TZFES,TZCAO,TZCAC,TZCAH,TZCAS,TZMGO,TZMGC,TZMGH
  real(r8) :: TZMGS,TZNAC,TZNAS,TZKAS,TH0PO4,TH1PO4,TH2PO4,TH3PO4
  real(r8) :: TZFE1P,TZFE2P,TZCA0P,TZCA1P,TZCA2P,TZMG1P,TH0POB
  real(r8) :: TH1POB,TH2POB,TH3POB,TFE1PB,TFE2PB,TCA0PB,TCA1PB
  real(r8) :: TCA2PB,TMG1PB,TXNH4,TXNHB,TXHY,TXAL,TXFE,TXCA
  real(r8) :: TXMG,TXNA,TXKA,TXHC,TXAL2,TXFE2,TXOH0,TXOH1,TXOH2
  real(r8) :: TXH1P,TXH2P,TXOH0B,TXOH1B,TXOH2B,TXH1PB,TXH2PB
  real(r8) :: TPALOH,TPFEOH,TPCACO,TPCASO,TPALPO,TPFEPO,TPCAPD
  real(r8) :: TPCAPH,TPCAPM,TPALPB,TPFEPB,TPCPDB,TPCPHB,TPCPMB
  real(r8) :: TCO2G,TCH4G,TCOZS,TCHFS,TOXYG,TOXYS,TZ2GG,TZ2GS
  real(r8) :: TZ2OG,TZ2OS,TZNH3G,TH2GG,TH2GS,TZNFN2,TZNFNI,TCO2GS
  real(r8) :: TCH4GS,TOXYGS,TZ2GSG,TZ2OGS,TH2GGS,TNH4GS,TNH3GS
  real(r8) :: TNO3GS,TNO2GS,TP14GS,TPO4GS,TXN4G,TXOH0G,TXOH1G
  real(r8) :: TXOH2G,TXH1PG,TXH2PG,TALPOG,TFEPOG,TCAPDG,TCAPHG
  real(r8) :: TCAPMG,TNH4FG,TNH3FG,TNHUFG,TNO3FG,TZNFNG,TVOLWR
  real(r8) :: TENGYR,TL,TI,TX
  real(r8) :: VHCPZ,VHCPY,VHCPO,VHCPX
  real(r8) :: VOLWXX,VOLIXX,VOLSLX,VOLXX,VOLTX,WX,WO
  real(r8) :: WI,WS
  real(r8) :: XCORP0
  real(r8) :: ZSI
  real(r8) :: ZXB,ZGI,ZNGGIN,ZN2OIN,ZNH3IN,ZG,Z4S,Z4X,Z4F,ZOS
  real(r8) :: ZOF,ZGB,Z2B,ZHB,ZNHUX0,ZNHUXI,ZNFNX0

  integer :: IFLGLS,LS,LS2

  real(r8) :: TOSGC(4,0:2),TOSGA(4,0:2),TOSGN(4,0:2),TOSGP(4,0:2) &
    ,TORXC(2,0:2),TORXN(2,0:2),TORXP(2,0:2),TOQGC(0:2),TOQGN(0:2) &
    ,TOQGP(0:2),TOQHC(0:2),TOQHN(0:2),TOQHP(0:2),TOHGC(0:2) &
    ,TOHGN(0:2),TOHGP(0:2),TOHGA(0:2),TOQGA(0:2),TOQHA(0:2) &
    ,TXCO2(JY,JX),DORGE(JY,JX),OMCL(0:JZ,JY,JX),OMNL(0:JZ,JY,JX) &
    ,THETCX(0:1),DVOLW(JZ,JY,JX)

  real(r8),allocatable :: TOMC(:,:,:,:)
  real(r8),allocatable :: TOMN(:,:,:,:)
  real(r8),allocatable :: TOMP(:,:,:,:)
  real(r8),allocatable :: TOMGC(:,:,:,:)
  real(r8),allocatable :: TOMGN(:,:,:,:)
  real(r8),allocatable :: TOMGP(:,:,:,:)

  real(r8),allocatable :: TOMCERff(:,:,:,:,:)
  real(r8),allocatable :: TOMNERff(:,:,:,:,:)
  real(r8),allocatable :: TOMPERff(:,:,:,:,:)
  real(r8),allocatable :: TOMCff(:,:,:)
  real(r8),allocatable :: TOMNff(:,:,:)
  real(r8),allocatable :: TOMPff(:,:,:)
  real(r8),allocatable :: TOMGCff(:,:,:)
  real(r8),allocatable :: TOMGNff(:,:,:)
  real(r8),allocatable :: TOMGPff(:,:,:)

  real(r8),allocatable ::  TORC(:,:)
  real(r8),allocatable ::  TORN(:,:)
  real(r8),allocatable ::  TORP(:,:)
  real(r8),allocatable ::  TOQC(:)
  real(r8),allocatable ::  TOQN(:)
  real(r8),allocatable ::  TOQP(:)
  real(r8),allocatable ::  TOQA(:)
  real(r8),allocatable ::  TOHC(:)
  real(r8),allocatable ::  TOHN(:)
  real(r8),allocatable ::  TOHP(:)
  real(r8),allocatable ::  TOHA(:)
  real(r8),allocatable ::  TOSC(:,:)
  real(r8),allocatable ::  TOSA(:,:)
  real(r8),allocatable ::  TOSN(:,:)
  real(r8),allocatable ::  TOSP(:,:)
  DATA THETCX/8.0E-06,8.0E-06/

  integer :: curday, curhour
  public :: redist, InitRedist, DestructRedist
  contains

  subroutine InitRedist
  implicit none

  call InitTflxType()

  allocate(TOMC(3,JG,7,0:jcplx1))
  allocate(TOMN(3,JG,7,0:jcplx1))
  allocate(TOMP(3,JG,7,0:jcplx1))
  allocate(TOMGC(3,JG,7,0:jcplx1))
  allocate(TOMGN(3,JG,7,0:jcplx1))
  allocate(TOMGP(3,JG,7,0:jcplx1))

  allocate(TOMCERff(3,JG,7,JY,JX))
  allocate(TOMNERff(3,JG,7,JY,JX))
  allocate(TOMPERff(3,JG,7,JY,JX))
  allocate(TOMCff(3,JG,7))
  allocate(TOMNff(3,JG,7))
  allocate(TOMPff(3,JG,7))
  allocate(TOMGCff(3,JG,7))
  allocate(TOMGNff(3,JG,7))
  allocate(TOMGPff(3,JG,7))

  allocate(TORC(2,0:jcplx1));   TORC=0._r8
  allocate(TORN(2,0:jcplx1));   TORN=0._r8
  allocate(TORP(2,0:jcplx1));   TORP=0._r8
  allocate(TOQC(0:jcplx1));     TOQC=0._r8
  allocate(TOQN(0:jcplx1));     TOQN=0._r8
  allocate(TOQP(0:jcplx1));     TOQP=0._r8
  allocate(TOQA(0:jcplx1));     TOQA=0._r8
  allocate(TOHC(0:jcplx1));     TOHC=0._r8
  allocate(TOHN(0:jcplx1));     TOHN=0._r8
  allocate(TOHP(0:jcplx1));     TOHP=0._r8
  allocate(TOHA(0:jcplx1));     TOHA=0._r8
  allocate(TOSC(jsken,0:jcplx1));TOSC=0._r8
  allocate(TOSA(jsken,0:jcplx1));TOSA=0._r8
  allocate(TOSN(jsken,0:jcplx1));TOSN=0._r8
  allocate(TOSP(jsken,0:jcplx1));TOSP=0._r8
  end subroutine InitRedist
!------------------------------------------------------------------------------------------

  subroutine DestructRedist
  implicit none
  call destroy(TOMC)
  call destroy(TOMN)
  call destroy(TOMP)
  call destroy(TOMGC)
  call destroy(TOMGN)
  call destroy(TOMGP)
  call destroy(TOMCERff)
  call destroy(TOMNERff)
  call destroy(TOMPERff)
  call destroy(TOMCff)
  call destroy(TOMNff)
  call destroy(TOMPff)
  call destroy(TOMGCff)
  call destroy(TOMGNff)
  call destroy(TOMGPff)

  call destroy(TORC)
  call destroy(TORN)
  call destroy(TORP)
  call destroy(TOQC)
  call destroy(TOQN)
  call destroy(TOQP)
  call destroy(TOQA)
  call destroy(TOHC)
  call destroy(TOHN)
  call destroy(TOHP)
  call destroy(TOHA)
  call destroy(TOSC)
  call destroy(TOSA)
  call destroy(TOSN)
  call destroy(TOSP)

  end subroutine DestructRedist
!------------------------------------------------------------------------------------------

  SUBROUTINE redist(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE UPDATES SOIL STATE VARIABLES WITH WATER, HEAT,
!     C, N, P, SOLUTE FLUXES CALCULATED IN EARLIER SUBROUTINES
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L,LL,K,M,N,N1,N2,N3,N4,N5,N6,NN,NO
  integer :: NZ,NUX,NR,L0,L1,LY,N4B,N5B
  integer :: LG
  real(r8) :: DORGC(JZ,JY,JX),DVOLI(JZ,JY,JX)
  real(r8) :: UDVOLI,UDLYXF
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: ENGYW
!     execution begins here
  curday=I
  curhour=J
  VOLISO=0.0_r8
  UDVOLI=0.0_r8
  UDLYXF=0.0_r8
  TFLWT=0.0_r8
  VOLPT=0.0_r8
  VOLTT=0.0_r8

  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      TXCO2(NY,NX)=0.0_r8
      DORGE(NY,NX)=0.0_r8
      WQRH(NY,NX)=0.0_r8
!
      call AddFluxToSurfaceResidue(NY,NX)
!
      call SinkChemicals(NY,NX)
!
!     RUNOFF AND SUBSURFACE BOUNDARY FLUXES
!
      call RunoffBal(I,J,NY,NX,NHW,NHE,NVN,NVS)
!
!     CHANGE EXTERNAL WATER TABLE DEPTH THROUGH DISTURBANCE
!
      call ModifyExWTBLByDisturbance(I,J,NY,NX)

      call LateralTranspt(I,J,NY,NX,LG)

      call SnowDynUpdate(NY,NX)

      call HandleSurfaceBoundary(I,NY,NX)
!
!     OVERLAND FLOW
      call OverlandFlow(NY,NX)

      call ChemicalOverlandFlow(NY,NX)
!
      call ChemicalBySnowRedistribution(NY,NX)
      !
      !     UPDATE STATE VARIABLES WITH TOTAL FLUXES CALCULATED ABOVE
      !
      !     IF(J.EQ.24)THEN
      !
      !     SNOWPACK VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
      !     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
      !     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS
      !
      DO 9785 L=1,JS
        WS=VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)*DENSI

        VOLWSO=VOLWSO+WS
        UVOLW(NY,NX)=UVOLW(NY,NX)+WS
        ENGYW=VHCPW(L,NY,NX)*TKW(L,NY,NX)
        HEATSO=HEATSO+ENGYW
        TLCO2G=TLCO2G+CO2W(L,NY,NX)+CH4W(L,NY,NX)
        UCO2S(NY,NX)=UCO2S(NY,NX)+CO2W(L,NY,NX)+CH4W(L,NY,NX)
        OXYGSO=OXYGSO+OXYW(L,NY,NX)
        TLN2G=TLN2G+ZNGW(L,NY,NX)+ZN2W(L,NY,NX)
        TLNH4=TLNH4+ZN4W(L,NY,NX)+ZN3W(L,NY,NX)
        TLNO3=TLNO3+ZNOW(L,NY,NX)
        TLPO4=TLPO4+Z1PW(L,NY,NX)+ZHPW(L,NY,NX)
        IF(ISALTG.NE.0)THEN
          SSW=ZALW(L,NY,NX)+ZFEW(L,NY,NX)+ZHYW(L,NY,NX)+ZCAW(L,NY,NX) &
            +ZMGW(L,NY,NX)+ZNAW(L,NY,NX)+ZKAW(L,NY,NX)+ZOHW(L,NY,NX) &
            +ZSO4W(L,NY,NX)+ZCLW(L,NY,NX)+ZCO3W(L,NY,NX)+H0PO4W(L,NY,NX) &
            +2.0*(ZHCO3W(L,NY,NX)+ZALH1W(L,NY,NX) &
            +ZALSW(L,NY,NX)+ZFEH1W(L,NY,NX)+ZFESW(L,NY,NX)+ZCAOW(L,NY,NX) &
            +ZCACW(L,NY,NX)+ZCASW(L,NY,NX)+ZMGOW(L,NY,NX)+ZMGCW(L,NY,NX) &
            +ZMGSW(L,NY,NX)+ZNACW(L,NY,NX)+ZNASW(L,NY,NX)+ZKASW(L,NY,NX) &
            +ZCA0PW(L,NY,NX)) &
            +3.0*(ZALH2W(L,NY,NX)+ZFEH2W(L,NY,NX)+ZCAHW(L,NY,NX) &
            +ZMGHW(L,NY,NX)+ZFE1PW(L,NY,NX)+ZCA1PW(L,NY,NX) &
            +ZMG1PW(L,NY,NX)) &
            +4.0*(ZALH3W(L,NY,NX)+ZFEH3W(L,NY,NX)+H3PO4W(L,NY,NX) &
            +ZFE2PW(L,NY,NX)+ZCA2PW(L,NY,NX)) &
            +5.0*(ZALH4W(L,NY,NX)+ZFEH4W(L,NY,NX))
          TION=TION+SSW
          !     WRITE(*,3335)'SSW',I,J,L,SSW
          !    2,ZALW(L,NY,NX),ZFEW(L,NY,NX),ZHYW(L,NY,NX),ZCAW(L,NY,NX)
          !    2,ZMGW(L,NY,NX),ZNAW(L,NY,NX),ZKAW(L,NY,NX),ZOHW(L,NY,NX)
          !    3,ZSO4W(L,NY,NX),ZCLW(L,NY,NX),ZCO3W(L,NY,NX),H0PO4W(L,NY,NX)
          !    4,ZHCO3W(L,NY,NX),ZALH1W(L,NY,NX)
          !    5,ZALSW(L,NY,NX),ZFEH1W(L,NY,NX),ZFESW(L,NY,NX),ZCAOW(L,NY,NX)
          !    6,ZCACW(L,NY,NX),ZCASW(L,NY,NX),ZMGOW(L,NY,NX),ZMGCW(L,NY,NX)
          !    7,ZMGSW(L,NY,NX),ZNACW(L,NY,NX),ZNASW(L,NY,NX),ZKASW(L,NY,NX)
          !    8,ZCA0PW(L,NY,NX)
          !    9,ZALH2W(L,NY,NX),ZFEH2W(L,NY,NX),ZCAHW(L,NY,NX)
          !    1,ZMGHW(L,NY,NX),ZFE1PW(L,NY,NX),ZCA1PW(L,NY,NX)
          !    2,ZMG1PW(L,NY,NX)
          !    2,ZALH3W(L,NY,NX),ZFEH3W(L,NY,NX),H3PO4W(L,NY,NX)
          !    3,ZFE2PW(L,NY,NX),ZCA2PW(L,NY,NX)
          !    5,ZALH4W(L,NY,NX),ZFEH4W(L,NY,NX)
        ENDIF
9785  CONTINUE
!
      call SumupSurfaceChemicalMass(NY,NX)
!
      call UpdateChemInLayers(NY,NX,LG,VOLISO,DORGC,DVOLI)
!
!     SNOWPACK LAYERING
      call SnowpackLayering(NY,NX)

!     watch out 'L' in the following call
      call RelayerSoilProfile(L,NY,NX,DORGC,DVOLI,UDVOLI,UDLYXF)

      TRN(NY,NX)=TRN(NY,NX)+HEATI(NY,NX)
      TLE(NY,NX)=TLE(NY,NX)+HEATE(NY,NX)
      TSH(NY,NX)=TSH(NY,NX)+HEATS(NY,NX)
      TGH(NY,NX)=TGH(NY,NX)-(HEATH(NY,NX)-HEATV(NY,NX))
      TLEC(NY,NX)=TLEC(NY,NX)+HEATE(NY,NX)*RAC(NY,NX)
      TSHC(NY,NX)=TSHC(NY,NX)+HEATS(NY,NX)*RAC(NY,NX)
      TCNET(NY,NX)=TCCAN(NY,NX)+HCO2G(NY,NX)
      RECO(NY,NX)=RECO(NY,NX)+HCO2G(NY,NX)
      TCAN(NY,NX)=TCAN(NY,NX)+TCCAN(NY,NX)
      TNPP(NY,NX)=TGPP(NY,NX)+TRAU(NY,NX)
      TNBP(NY,NX)=TCAN(NY,NX)+UCO2G(NY,NX)+UCH4G(NY,NX) &
        -UDOCQ(NY,NX)-UDICQ(NY,NX)-UDOCD(NY,NX)-UDICD(NY,NX) &
        +TXCO2(NY,NX)
      !     WRITE(*,6646)'TNBP',I,J,NX,NY,TNBP(NY,NX)
      !    2,TCAN(NY,NX),TCCAN(NY,NX),TCNET(NY,NX)
      !    2,TNPP(NY,NX),THRE(NY,NX),UCO2G(NY,NX),UCH4G(NY,NX)
      !    2,UDOCQ(NY,NX),UDICQ(NY,NX),UDOCD(NY,NX),UDICD(NY,NX)
      !    3,TXCO2(NY,NX)
      IF(NU(NY,NX).GT.NUI(NY,NX))THEN
        DO 235 L=NUI(NY,NX),NU(NY,NX)-1
          IF(VOLX(L,NY,NX).LE.ZEROS2(NY,NX))THEN
            TKS(L,NY,NX)=TKS(NU(NY,NX),NY,NX)
            TCS(L,NY,NX)=TKS(L,NY,NX)-TC2K
          ENDIF
235     CONTINUE
      ENDIF
      !
      !     MIX ALL SOIL STATE VARIABLES AND INCORPORATE ALL SURFACE
      !     RESIDUE STATE VARIABLES WITHIN THE TILLAGE ZONE TO THE EXTENT
      !     ASSOCIATED IN 'DAY' WITH EACH TILLAGE EVENT ENTERED IN THE
      !     TILLAGE FILE
      call ApplyMixing(I,J,NY,NX)

      !
      !     OUTPUT FOR SOIL WATER, ICE CONTENTS
      !
      THETWZ(0,NY,NX)=AMAX1(0.0,(VOLW(0,NY,NX)-VOLWRX(NY,NX))/AREA(3,0,NY,NX))
      THETIZ(0,NY,NX)=AMAX1(0.0,(VOLI(0,NY,NX)-VOLWRX(NY,NX))/AREA(3,0,NY,NX))
      !     THETWZ(0,NY,NX)=AMAX1(0.0,AMIN1(1.0
      !     2,VOLW(0,NY,NX)/VOLR(NY,NX)))
      !     THETIZ(0,NY,NX)=AMAX1(0.0,AMIN1(1.0
      !    2,VOLI(0,NY,NX)/VOLR(NY,NX)))
      DO 9945 L=NUI(NY,NX),NL(NY,NX)
        VOLXX=AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)*FMPR(L,NY,NX)
        VOLTX=VOLXX+VOLAH(L,NY,NX)
        THETWZ(L,NY,NX)=safe_adb(VOLW(L,NY,NX)+AMIN1(VOLAH(L,NY,NX) &
          ,VOLWH(L,NY,NX)),VOLTX)
        THETIZ(L,NY,NX)=safe_adb(VOLI(L,NY,NX)+AMIN1(VOLAH(L,NY,NX) &
          ,VOLIH(L,NY,NX)),VOLTX)
        !     IF(NX.EQ.1)THEN
        !     WRITE(*,1189)'THETWZ',I,J,NX,NY,L
        !    2,THETWZ(L,NY,NX),THETIZ(L,NY,NX),VOLW(L,NY,NX)
        !    2,VOLWH(L,NY,NX),VOLI(L,NY,NX),VOLIH(L,NY,NX)
        !    2,VOLTX,VOLXX,DLYR(3,L,NY,NX)
!1189  FORMAT(A8,5I4,20E12.4)
!       ENDIF
9945  CONTINUE
!
!     CHECK MATERIAL BALANCES
!
      IF(I.EQ.365.AND.J.EQ.24)THEN
        WRITE(19,2221)'ORGC',I,J,IYRC,NX,NY &
          ,(ORGC(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
        WRITE(20,2221)'ORGN',I,J,IYRC,NX,NY &
          ,(ORGN(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
2221    FORMAT(A8,5I6,21E14.6)
      ENDIF
      !     IF(J.LE.2)THEN
      !     WRITE(20,2221)'OMCL',I,J,IYRC,NX,NY
      !    2,(OMCL(L,NY,NX),L=0,NL(NY,NX))
      !     WRITE(20,2221)'OMNL',I,J,IYRC,NX,NY
      !    2,(OMNL(L,NY,NX),L=0,NL(NY,NX))
      !     WRITE(20,2222)'TLC',I,J,IYRC,NX,NY,TLRSDC+TLORGC+TLCO2G-CO2GIN
      !    2+TCOU-TORGF-XCSN,TLRSDC,TLORGC,TLCO2G,CO2GIN,TCOU,TORGF,XCSN
      !    3,(ORGR(L,NY,NX),L=0,NL(NY,NX))
      !    5,XCODFS(NY,NX),XCOFLG(3,NU(NY,NX),NY,NX),TCO2Z(NY,NX)
      !    2,FLQGQ(NY,NX)*CCOR(NY,NX),FLQGI(NY,NX)*CCOQ(NY,NX),XCODFG(0,NY,NX)
      !     3,XCODFR(NY,NX),XCHDFS(NY,NX),XCHFLG(3,NU(NY,NX),NY,NX)
      !    2,FLQGQ(NY,NX)*CCHR(NY,NX),FLQGI(NY,NX)*CCHQ(NY,NX),XCHDFG(0,NY,NX)
      !    3,XCHDFR(NY,NX),PRECU(NY,NX)*CCOQ(NY,NX),PRECU(NY,NX)*CCHQ(NY,NX)
      !    6,TCOQRS(NY,NX),TCHQRS(NY,NX),XCOFLS(1,0,NY,NX+1)
      !    7,XCOFLS(2,0,NY+1,NX)
      !    3,UCOP(NY,NX),UDOCQ(NY,NX),UDICQ(NY,NX),UDOCD(NY,NX),UDICD(NY,NX)
      !    2,(((CSNT(M,K,L,NY,NX),M=1,4),K=0,1),L=0,NJ(NY,NX))
      !    3,(TCO2P(L,NY,NX),L=1,NJ(NY,NX)),(TCO2S(L,NY,NX),L=1,NJ(NY,NX))
      !    4,CQ,ZCSNC(NY,NX)
      !     WRITE(20,2222)'TLW',I,J,IYRC,NX,NY,VOLWSO-CRAIN+CRUN+CEVAP+VOLWOU
      !    2,VOLWSO,CRAIN,CRUN,CEVAP,VOLWOU
      !    3,TVOLWC(NY,NX),TVOLWP(NY,NX),VOLW(0,NY,NX),VOLI(0,NY,NX)*DENSI
      !    4,TFLWC(NY,NX),TEVAPC(NY,NX),TEVAPG(NY,NX),TEVAPP(NY,NX)
      !    5,VOLSS(NY,NX),VOLWS(NY,NX),VOLIS(NY,NX)*DENSI
      !    6,(VOLW(L,NY,NX),L=0,NL(NY,NX))
      !    6,(VOLI(L,NY,NX)*DENSI,L=0,NL(NY,NX))
      !    6,TQS(NY,NX),TQW(NY,NX),TQI(NY,NX),TVOLWC(NY,NX)
      !    7,TVOLWP(NY,NX),(TUPWTR(L,NY,NX),L=1,JZ)
      !     WRITE(20,2222)'TLH',I,J,IYRC,NX,NY,HEATSO-HEATIN+HEATOU
      !    2,HEATSO,HEATIN,HEATOU,HEATH(NY,NX),HRAINR
      !    I,HTHAWR(NY,NX),HFLXD,HFLXO,cpw*TKA(NY,NX)*PRECA(NY,NX)
      !    I,cps*TKA(NY,NX)*PRECW(NY,NX),HEATH(NY,NX),THFLXC(NY,NX)
      !    I,(XTHAWW(L,NY,NX),L=1,JS)
      !    I,(THTHAW(L,NY,NX),L=NU(NY,NX),NL(NY,NX))
      !    I,(TUPHT(L,NY,NX),L=NU(NY,NX),NL(NY,NX))
      !    O,cpw*TKA(NY,NX)*PRECU(NY,NX),THQS(NY,NX)
      !    S,((VHCP(L,NY,NX)*TKS(L,NY,NX)),L=0,NL(NY,NX))
      !    S,((VHCPW(L,NY,NX)*TKW(L,NY,NX)),L=1,JS)
      !    S,VHCP(0,NY,NX),VHCPRX(NY,NX)
      !     WRITE(19,2224)'TLO',I,J,IYRC,NX,NY,OXYGSO-OXYGIN+OXYGOU
      !    2,OXYGSO,OXYGIN,OXYGOU
      !    3,XOXDFS(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX)
      !    3,XOXDFG(0,NY,NX),TOXYZ(NY,NX),FLQGQ(NY,NX)*COXR(NY,NX)
      !    4,FLQGI(NY,NX)*COXQ,PRECU(NY,NX)*COXQ
      !    3,(OXYG(L,NY,NX),L=1,JZ)
      !    4,(OXYS(L,NY,NX),L=0,JZ)
      !    4,(OXYSH(L,NY,NX),L=1,JZ)
      !    2,(RUPOXO(L,NY,NX),L=1,NJ(NY,NX))
      !    3,(TUPOXP(L,NY,NX),L=1,NJ(NY,NX)),(TOXFLA(L,NY,NX),L=1,NJ(NY,NX))
      !     WRITE(20,2222)'TLN',I,J,IYRC,NX,NY,TLRSDN+TLORGN+TLN2G+TLNH4
      !    2+TLNO3-ZN2GIN-TZIN+TZOU-TORGN-XZSN,TLRSDN,TLORGN,TLN2G,TLNH4
      !    3,TLNO3,ZN2GIN,TZIN,TZOU,TORGN,XZSN
      !    5,(XN4(L,NY,NX),L=0,NL(NY,NX))
      !    4,(((ZSNT(M,K,L,NY,NX),M=1,4),K=0,1),L=0,JZ)
      !    5,(TUPNH4(L,NY,NX),L=1,JZ)
      !    6,(TUPNO3(L,NY,NX),L=1,JZ),(TNHFLA(L,NY,NX),L=1,JZ)
      !    7,XN3DFS(NY,NX),XNBDFS(NY,NX)
      !    8,XN3FLG(3,NU(NY,NX),NY,NX),TNH3Z(NY,NX),UN2GS(NY,NX)
      !    9,(XN2GS(L,NY,NX),L=0,JZ)
      !    4,PRECQ(NY,NX),PRECR(NY,NX)
      !    4,PRECW(NY,NX),PRECI(NY,NX),FLQGM(NY,NX),FLQRM(NY,NX)
      !     WRITE(20,2222)'TLP',I,J,IYRC,NX,NY,TLRSDP+TLORGP
      !    2+TLPO4-TPIN+TPOU-TORGP-XPSN,TLRSDP,TLORGP
      !    2,TLPO4,TPIN,TPOU,TORGP,XPSN
      !    2,(Z1PW(L,NY,NX),L=1,JS),(ZHPW(L,NY,NX),L=1,JS)
      !    3,(H1PO4(L,NY,NX),L=0,JZ),(H2PO4(L,NY,NX),L=0,JZ)
      !    4,(H1PO4H(L,NY,NX),L=1,JZ),(H2PO4H(L,NY,NX),L=1,JZ)
      !    6,FLQGQ(NY,NX),FLQRQ(NY,NX),CPOR(NY,NX),CH1PR(NY,NX)
      !    7,FLQGI(NY,NX),FLQRI(NY,NX),CPOQ(I,NY,NX),CH1PQ(I,NY,NX)
      !     WRITE(20,2222)'TLI',I,J,IYRC,NX,NY,TION-TIONIN+TIONOU
      !    2,TION,TIONIN,TIONOU
      !    3,PRECQ(NY,NX),XHGDFS(NY,NX),XHGFLG(3,NU(NY,NX),NY,NX)
      !    4,TH2GZ(NY,NX)
      !    4,(XHGQRS(N,NY,NX),N=1,2),(RH2GO(L,NY,NX),L=1,JZ)
      !    5,(THGFLA(L,NY,NX),L=1,JZ)
      !    6,(H2GG(L,NY,NX),L=1,JZ),(TLH2GP(L,NY,NX),L=1,JZ)
      !     WRITE(20,2224)'TLG',I,J,IYRC,NX,NY,TLH2G-H2GIN+H2GOU,TLH2G
      !    2,H2GIN,H2GOU,(H2GG(L,NY,NX),L=0,NJ(NY,NX))
      !    3,(H2GS(L,NY,NX),L=0,NJ(NY,NX))
      !    3,(H2GSH(L,NY,NX),L=1,NJ(NY,NX)),(TLH2GP(L,NY,NX),L=1,NJ(NY,NX))
      !    4,XHGDFS(NY,NX),TH2GZ(NY,NX),(THGFLA(L,NY,NX),L=1,NJ(NY,NX))
      !    2,XHGDFG(0,NY,NX),XHGDFR(NY,NX),(XHGBBL(L,NY,NX),L=1,NJ(NY,NX))
      !    3,(RH2GO(L,NY,NX),L=0,NJ(NY,NX)),(TUPHGS(L,NY,NX),L=1,NJ(NY,NX))
      !    4,(XHGQRS(N,NY,NX),N=1,2),((XHGFLS(N,L,NY,NX),N=1,3),L=0,NJ(NY,NX))
      !    5,((XHGFHS(N,L,NY,NX),N=1,3),L=0,NJ(NY,NX))
      !    6,((XHGFLG(N,L,NY,NX),N=1,3),L=0,NJ(NY,NX))
      !     WRITE(*,2223)'TLS',I,J,IYRC,NX,NY,NU(NY,NX),TSEDSO+TSEDOU
      !    2,TSEDSO,TSEDOU,USEDOU(NY,NX),DLYR(3,NU(NY,NX),NY,NX)
      !    3,BKVL(NU(NY,NX),NY,NX),SAND(NU(NY,NX),NY,NX)
      !    4,SILT(NU(NY,NX),NY,NX),CLAY(NU(NY,NX),NY,NX)
      !    5,ORGC(NU(NY,NX),NY,NX)
      !2222  FORMAT(A8,5I6,240F18.6)
      !2223  FORMAT(A8,6I6,160F16.6)
      !2224  FORMAT(A8,5I6,160F16.6)
      !     ENDIF
9990  CONTINUE
9995  CONTINUE
  RETURN

  END subroutine redist

!------------------------------------------------------------------------------------------

  subroutine ModifyExWTBLByDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
!     begin_execution
!     ZNOON=hour of solar noon from weather file
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     DCORP=intensity (fire) or depth (tillage,drainage) of disturbance
!     CDPTH(NU=soil surface elevation
!     DTBLI,DTBLDI=depth of natural,artificial water table from readi.f
!     DTBLX,DTBLZ=current,initial natural water table depth
!     DTBLY,DTBLD=current,initial artificial water table depth
!     IDTBL=water table flag from readi.f
!        :0=none
!        :1,2=natural stationary,mobile
!        :3,4=artificial stationary,mobile
!     HVOLO=hourly water loss through lateral and lower boundaries
!
  IF(J.EQ.INT(ZNOON(NY,NX)).AND.ITILL(I,NY,NX).EQ.23)THEN
    DCORPW=DCORP(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
    DTBLI(NY,NX)=DCORPW
    DTBLZ(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX)) &
      *(1.0-DTBLG(NY,NX))
    DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
  ENDIF
  IF(J.EQ.INT(ZNOON(NY,NX)).AND.ITILL(I,NY,NX).EQ.24)THEN
    DCORPW=DCORP(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
    IF(IDTBL(NY,NX).EQ.1)THEN
      IDTBL(NY,NX)=3
    ELSEIF(IDTBL(NY,NX).EQ.2)THEN
      IDTBL(NY,NX)=4
    ENDIF
    DTBLDI(NY,NX)=DCORPW
    DTBLD(NY,NX)=AMAX1(0.0,DTBLDI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX)) &
      *(1.0-DTBLG(NY,NX)))
    DTBLY(NY,NX)=DTBLD(NY,NX)
  ENDIF
!
!     SET DEPTH OF MOBILE EXTERNAL WATER TABLE
!
  IF(IDTBL(NY,NX).EQ.2.OR.IDTBL(NY,NX).EQ.4)THEN
!     DTBLX(NY,NX)=DTBLX(NY,NX)-0.0*HVOLO(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
!    2-0.00167*(DTBLX(NY,NX)-DTBLZ(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
    DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
  ENDIF
  IF(IDTBL(NY,NX).EQ.4)THEN
    DTBLY(NY,NX)=DTBLY(NY,NX)-0.0*HVOLO(NY,NX)/AREA(3,NU(NY,NX),NY,NX) &
      -0.00167*(DTBLY(NY,NX)-DTBLD(NY,NX))
  ENDIF
  end subroutine ModifyExWTBLByDisturbance
!------------------------------------------------------------------------------------------

  subroutine HandleSurfaceBoundary(I,NY,NX)
  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: L,K
  real(r8):: tkspre,vhcp1s
  ! begin_execution
  !
  ! CALCULATE SURFACE RESIDUE TEMPERATURE FROM ITS CHANGE
  ! IN HEAT STORAGE
  !
  VHCPZ=VHCP(0,NY,NX)             !old heat capacity
  VHCPY=cpw*VOLW(0,NY,NX)+cpi*VOLI(0,NY,NX)+cpo*ORGC(0,NY,NX) !new heat capacity
  VHCPO=VHCPY-VHCPZ               !change in heat capacity
  HFLXO=VHCPO*TKA(NY,NX)          !tka: air temperature in kelvin, hflxo represents incoming heat
  !update water and ice content in residue
  VOLW(0,NY,NX)=max(VOLW(0,NY,NX)+FLWR(NY,NX)+THAWR(NY,NX)+TQR(NY,NX),0._r8)
  VOLI(0,NY,NX)=max(VOLI(0,NY,NX)-THAWR(NY,NX)/DENSI,0._r8)
  ENGYZ=VHCPZ*TKS(0,NY,NX)
  !update heat caapcity
  VHCP(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX)+cpi*VOLI(0,NY,NX)
  IF(VHCP(0,NY,NX).GT.VHCPRX(NY,NX))THEN
    !when there are still significant heat capacity of the residual layer
    tkspre=TKS(0,NY,NX)
    TKS(0,NY,NX)=(ENGYZ+HFLWR(NY,NX)+HTHAWR(NY,NX)+HFLXO &
      +THQR(NY,NX))/VHCP(0,NY,NX)
    HEATIN=HEATIN+HFLXO
    Ls=NUM(NY,NX)
    !if(curday>=175)write(*,*)'at line',__LINE__,TKS(0,NY,NX),tks(Ls,ny,nx),tkspre
    if(abs(TKS(0,NY,NX)-tks(Ls,ny,nx))>10._r8)then
      vhcp1s=VHCM(Ls,NY,NX)+cpw*(VOLW(Ls,NY,NX)+VOLWH(Ls,NY,NX))+&
        cpi*(VOLI(Ls,NY,NX)+VOLIH(Ls,NY,NX))
      !if(curday>=175)then
      !  write(*,*)'LS=',LS
      !  write(*,*)'VHCM(Ls,NY,NX),BKDS(L,NY,NX)=',VHCM(Ls,NY,NX),BKDS(Ls,NY,NX)
      !  write(*,*)'VOLW(Ls,NY,NX),VOLWH(Ls,NY,NX)=',VOLW(Ls,NY,NX),VOLWH(Ls,NY,NX)
      !  write(*,*)'VOLI(Ls,NY,NX),VOLIH(Ls,NY,NX)=',VOLI(Ls,NY,NX),VOLIH(Ls,NY,NX)
      !  write(*,*)'at line',__LINE__,vhcp1s
      !endif
      TKS(0,NY,NX)=(tks(0,ny,nx)*VHCP(0,NY,NX)+TKS(Ls,NY,NX)*vhcp1s)/(VHCP(0,NY,NX)+vhcp1s)
      TKS(Ls,NY,NX)=TKS(0,NY,NX)
      !if(curday>=175)write(*,*)'at line',__LINE__,TKS(0,NY,NX),tks(Ls,ny,nx),tkspre
!      write(*,*)'\n'
!      write(*,*)'curday       =',curday,curhour
!      write(*,*)'TKS(0,NY,NX) =',tkspre,TKS(0,NY,NX),TKS(Ls,NY,NX),tka(ny,nx)
!      write(*,*)'ENGYZ        =',ENGYZ
!      write(*,*)'HFLWR(NY,NX) =',HFLWR(NY,NX)
!      write(*,*)'HTHAWR(NY,NX)=',HTHAWR(NY,NX)
!      write(*,*)'HFLXO        =',HFLXO
!      write(*,*)'THQR(NY,NX)  =',THQR(NY,NX)
!      write(*,*)'forc         =',(ENGYZ+HFLWR(NY,NX)+HTHAWR(NY,NX)+HFLXO+THQR(NY,NX))
!      write(*,*)'VHCP(0,NY,NX)=',VHCP(0,NY,NX),vhcp1s, &

!      write(*,*)'VHCPZ        =',VHCPZ
!      write(*,*)'ORGC(0,NY,NX)=',ORGC(0,NY,NX)
!      write(*,*)'VOLW(0,NY,NX)=',VOLW(0,NY,NX)
!      write(*,*)'VOLI(0,NY,NX)=',VOLI(0,NY,NX)
!      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ELSE
    HEATIN=HEATIN+HFLXO+(TKS(NUM(NY,NX),NY,NX)-TKS(0,NY,NX))*VHCP(0,NY,NX)
    TKS(0,NY,NX)=TKS(NUM(NY,NX),NY,NX)
    if(TKS(0,NY,NX)>350._r8 .or. tks(0,ny,nx)<200._r8)then
       write(*,*)'TKS(0,NY,NX)=',TKS(0,NY,NX),NUM(NY,NX)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDIF
  ENGYR=VHCP(0,NY,NX)*TKS(0,NY,NX)
  HEATSO=HEATSO+ENGYR
  HEATIN=HEATIN+HTHAWR(NY,NX)
  TCS(0,NY,NX)=TKS(0,NY,NX)-TC2K
  TSMX(0,NY,NX)=AMAX1(TSMX(0,NY,NX),TCS(0,NY,NX))
  TSMN(0,NY,NX)=AMIN1(TSMN(0,NY,NX),TCS(0,NY,NX))
  !     UVOLW(NY,NX)=UVOLW(NY,NX)-VOLW(0,NY,NX)-VOLI(0,NY,NX)*DENSI
  !     IF(I.GT.350.AND.NX.EQ.1)THEN
  !     WRITE(*,6634)'VOLWR',I,J,NX,NY,VOLW(0,NY,NX),VOLI(0,NY,NX)
  !    2,ORGC(0,NY,NX),FLWR(NY,NX),TQR(NY,NX),THAWR(NY,NX)
  !    3,TKS(0,NY,NX),HFLWR(NY,NX),HFLXO,THQR(NY,NX)
  !    4,HTHAWR(NY,NX),VHCP(0,NY,NX),TKS(NUM(NY,NX),NY,NX)
  !    5,ENGYR,ENGYZ,VHCPZ,VHCPRX(NY,NX)
!6634  FORMAT(A8,4I4,30E14.6)
  !     ENDIF
  !
  !     SURFACE BOUNDARY WATER FLUXES
  !
  WI=PRECQ(NY,NX)+PRECI(NY,NX)   !total incoming water flux=rain/snowfall + irrigation
  CRAIN=CRAIN+WI
  URAIN(NY,NX)=URAIN(NY,NX)+WI
  WO=TEVAPG(NY,NX)+TEVAPP(NY,NX) !total outgoing water flux
  !  if(WO/=WO)then
  !    write(*,*)'TEVAPG=',TEVAPG(NY,NX)
  !    write(*,*)'TEVAPP=',TEVAPP(NY,NX)
  !    call endrun(msg='NaN encounterd in '//mod_filename,line=__LINE__)
  !  endif
  CEVAP=CEVAP-WO
  UEVAP(NY,NX)=UEVAP(NY,NX)-WO
  VOLWOU=VOLWOU-PRECU(NY,NX)
  HVOLO(NY,NX)=HVOLO(NY,NX)-PRECU(NY,NX)
  UVOLO(NY,NX)=UVOLO(NY,NX)-PRECU(NY,NX)
  UDRAIN(NY,NX)=UDRAIN(NY,NX)+FLW(3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY HEAT FLUXES
  !
  HEATIN=HEATIN+cpw*TKA(NY,NX)*PRECA(NY,NX) &
      +cps*TKA(NY,NX)*PRECW(NY,NX)
  HEATIN=HEATIN+HEATH(NY,NX)+THFLXC(NY,NX)
  DO 5150 L=1,JS
    HEATIN=HEATIN+XTHAWW(L,NY,NX)
5150  CONTINUE
  HEATOU=HEATOU-cpw*TKA(NY,NX)*PRECU(NY,NX)
!
! SURFACE BOUNDARY CO2, CH4 AND DOC FLUXES
!
  CI=XCODFS(NY,NX)+XCOFLG(3,NU(NY,NX),NY,NX)+TCO2Z(NY,NX) &
      +(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX) &
      +XCODFG(0,NY,NX)+XCODFR(NY,NX)
  CH=XCHDFS(NY,NX)+XCHFLG(3,NU(NY,NX),NY,NX)+TCH4Z(NY,NX) &
      +(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCHR(NY,NX) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX))*CCHQ(NY,NX) &
      +XCHDFG(0,NY,NX)+XCHDFR(NY,NX)
  CO=-PRECU(NY,NX)*CCOQ(NY,NX)
  CX=-PRECU(NY,NX)*CCHQ(NY,NX)
  UCO2G(NY,NX)=UCO2G(NY,NX)+CI
  HCO2G(NY,NX)=HCO2G(NY,NX)+CI
  UCH4G(NY,NX)=UCH4G(NY,NX)+CH
  HCH4G(NY,NX)=HCH4G(NY,NX)+CH
  CO2GIN=CO2GIN+CI+CH
  TCOU=TCOU+CO+CX
  !
  !     SURFACE BOUNDARY O2 FLUXES
  !
  OI=XOXDFS(NY,NX)+XOXFLG(3,NU(NY,NX),NY,NX)+TOXYZ(NY,NX) &
      +(FLQGQ(NY,NX)+FLQRQ(NY,NX))*COXR(NY,NX) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX))*COXQ(NY,NX) &
      +XOXDFG(0,NY,NX)+XOXDFR(NY,NX)
  OXYGIN=OXYGIN+OI
  OO=RUPOXO(0,NY,NX)-PRECU(NY,NX)*COXQ(NY,NX)
  OXYGOU=OXYGOU+OO
  UOXYG(NY,NX)=UOXYG(NY,NX)+OI
  HOXYG(NY,NX)=HOXYG(NY,NX)+OI
  HI=XHGDFS(NY,NX)+XHGFLG(3,NU(NY,NX),NY,NX)+TH2GZ(NY,NX) &
      +XHGDFG(0,NY,NX)+XHGDFR(NY,NX)
  H2GIN=H2GIN+HI
  HO=RH2GO(0,NY,NX)
  H2GOU=H2GOU+HO
  !     IF(I.EQ.256)THEN
  !     WRITE(*,6646)'UCO2G',I,J,NX,NY,UCO2G(NY,NX),UOXYG(NY,NX),CI,OI
  !    2,XCODFS(NY,NX),XCOFLG(3,NU(NY,NX),NY,NX),TCO2Z(NY,NX)
  !    2,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX)
  !    3,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX)
  !    4,XCODFG(0,NY,NX),XCODFR(NY,NX)
  !    5,XOXDFS(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX),TOXYZ(NY,NX)
  !    2,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*COXR(NY,NX)
  !    3,(FLQGI(NY,NX)+FLQRI(NY,NX))*COXQ(NY,NX)
  !    4,XOXDFG(0,NY,NX),XOXDFR(NY,NX)
  !    5,(TLCO2P(L,NY,NX),L=1,10)
  !    6,(TLOXYP(L,NY,NX),L=1,10)
  !     WRITE(*,6646)'UCH4G',I,J,NX,NY,UCH4G(NY,NX),CH
  !    2,XCHDFS(NY,NX),XCHFLG(3,NU(NY,NX),NY,NX),TCH4Z(NY,NX)
  !    2,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCHR(NY,NX)
  !    3,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCHQ(NY,NX)
  !    4,XCHDFG(0,NY,NX),XCHDFR(NY,NX)
  !    5,(TLCH4P(L,NY,NX),L=1,10)
!6646  FORMAT(A8,4I4,60E12.4)
  !     ENDIF
  !
  !     SURFACE BOUNDARY N2, N2O, NH3, NH4, NO3, AND DON FLUXES
  !
  ZSI=((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(CN4R(NY,NX)+CN3R(NY,NX)+CNOR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))*14.0
  ZXB=-PRECU(NY,NX)*(CNNQ(NY,NX)+CN2Q(NY,NX))-PRECU(NY,NX) &
      *(CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX))*14.0
  TZIN=TZIN+ZSI
  TZOU=TZOU+ZXB
  ZGI=(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX)) &
      +XNGDFS(NY,NX)+XN2DFS(NY,NX)+XN3DFS(NY,NX) &
      +XNBDFS(NY,NX)+XNGFLG(3,NU(NY,NX),NY,NX) &
      +XN2FLG(3,NU(NY,NX),NY,NX)+XN3FLG(3,NU(NY,NX),NY,NX) &
      +TN2OZ(NY,NX)+TNH3Z(NY,NX) &
      +XN2DFG(0,NY,NX)+XNGDFG(0,NY,NX)+XN3DFG(0,NY,NX) &
      +XNGDFR(NY,NX)+XN2DFR(NY,NX)+XN3DFR(NY,NX)
  ZN2GIN=ZN2GIN+ZGI
  ZDRAIN(NY,NX)=ZDRAIN(NY,NX)+XN4FLW(3,NK(NY,NX),NY,NX) &
      +XN3FLW(3,NK(NY,NX),NY,NX)+XNOFLW(3,NK(NY,NX),NY,NX) &
      +XNXFLS(3,NK(NY,NX),NY,NX)+XN4FLB(3,NK(NY,NX),NY,NX) &
      +XN3FLB(3,NK(NY,NX),NY,NX)+XNOFLB(3,NK(NY,NX),NY,NX) &
      +XNXFLB(3,NK(NY,NX),NY,NX)
  ZNGGIN=XNGDFS(NY,NX)+XNGFLG(3,NU(NY,NX),NY,NX)+XNGDFG(0,NY,NX)
  ZN2OIN=XN2DFS(NY,NX)+XN2FLG(3,NU(NY,NX),NY,NX)+XN2DFG(0,NY,NX)
  ZNH3IN=XN3DFS(NY,NX)+XNBDFS(NY,NX)+XN3FLG(3,NU(NY,NX),NY,NX) &
      +XN3DFG(0,NY,NX)
  ! UN2GG(NY,NX)=UN2GG(NY,NX)+ZNGGIN
  !     HN2GG(NY,NX)=HN2GG(NY,NX)+ZNGGIN
  UN2OG(NY,NX)=UN2OG(NY,NX)+ZN2OIN
  HN2OG(NY,NX)=HN2OG(NY,NX)+ZN2OIN
  UNH3G(NY,NX)=UNH3G(NY,NX)+ZNH3IN
  HNH3G(NY,NX)=HNH3G(NY,NX)+ZNH3IN
  UN2GS(NY,NX)=UN2GS(NY,NX)+XN2GS(0,NY,NX)
  UH2GG(NY,NX)=UH2GG(NY,NX)+HI
  !
  !     SURFACE BOUNDARY PO4 AND DOP FLUXES
  !
  PI=31.0_r8*((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(CPOR(NY,NX)+CH1PR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(CPOQ(I,NY,NX)+CH1PQ(I,NY,NX)))
  PXB=-31.0_r8*PRECU(NY,NX)*(CPOQ(I,NY,NX)+CH1PQ(I,NY,NX))
  TPIN=TPIN+PI
  TPOU=TPOU+PXB
  PDRAIN(NY,NX)=PDRAIN(NY,NX)+XH2PFS(3,NK(NY,NX),NY,NX) &
      +XH2BFB(3,NK(NY,NX),NY,NX)+XH1PFS(3,NK(NY,NX),NY,NX) &
      +XH1BFB(3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY ION FLUXES
  !
  SIN=((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(2.0_r8*CN4R(NY,NX)+CN3R(NY,NX)+CNOR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(2.0*CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))
  SGN=(2.0_r8*(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX)) &
      +2.0_r8*(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX)) &
      +2.0_r8*(XNGDFS(NY,NX)+XN2DFS(NY,NX))+XN3DFS(NY,NX) &
      +XNBDFS(NY,NX)+2.0_r8*(XNGFLG(3,NU(NY,NX),NY,NX) &
      +XN2FLG(3,NU(NY,NX),NY,NX))+XN3FLG(3,NU(NY,NX),NY,NX) &
      +2.0_r8*TN2OZ(NY,NX)+TNH3Z(NY,NX) &
      +2.0_r8*(XN2DFG(0,NY,NX)+XNGDFG(0,NY,NX))+XN3DFG(0,NY,NX) &
      +2.0_r8*(XNGDFR(NY,NX)+XN2DFR(NY,NX))+XN3DFR(NY,NX))/14.0_r8
  SIP=((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(3.0_r8*CPOR(NY,NX)+2.0_r8*CH1PR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(3.0_r8*CPOQ(I,NY,NX)+2.0_r8*CH1PQ(I,NY,NX)))
  SNB=-PRECU(NY,NX)*(CNNQ(NY,NX)+CN2Q(NY,NX))-PRECU(NY,NX) &
      *(2.0_r8*CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX))
      SPB=-PRECU(NY,NX)*(3.0*CPOQ(I,NY,NX)+2.0*CH1PQ(I,NY,NX))
  SNM0=(2.0_r8*XNH4S(0,NY,NX)+XNO3S(0,NY,NX)+XNO2S(0,NY,NX) &
      -2.0_r8*XN2GS(0,NY,NX))/14.0_r8
  SPM0=(2.0_r8*XH1PS(0,NY,NX)+3.0_r8*XH2PS(0,NY,NX))/31.0_r8
  !
  !     ACCUMULATE PLANT LITTERFALL FLUXES
  !
  XCSN=XCSN+ZCSNC(NY,NX)
  XZSN=XZSN+ZZSNC(NY,NX)
  XPSN=XPSN+ZPSNC(NY,NX)
  UXCSN(NY,NX)=UXCSN(NY,NX)+ZCSNC(NY,NX)
  UXZSN(NY,NX)=UXZSN(NY,NX)+ZZSNC(NY,NX)
  UXPSN(NY,NX)=UXPSN(NY,NX)+ZPSNC(NY,NX)
  !
  !     SURFACE BOUNDARY SALT FLUXES FROM RAINFALL AND SURFACE IRRIGATION
  !
  IF(ISALTG.NE.0)THEN
    SIR=PRECQ(NY,NX)*(CALR(NY,NX)+CFER(NY,NX)+CHYR(NY,NX)+CCAR(NY,NX) &
      +CMGR(NY,NX)+CNAR(NY,NX)+CKAR(NY,NX)+COHR(NY,NX)+CSOR(NY,NX) &
      +CCLR(NY,NX)+CC3R(NY,NX)+CH0PR(NY,NX) &
      +2.0_r8*(CHCR(NY,NX)+CAL1R(NY,NX)+CALSR(NY,NX)+CFE1R(NY,NX) &
      +CFESR(NY,NX)+CCAOR(NY,NX)+CCACR(NY,NX)+CCASR(NY,NX)+CMGOR(NY,NX) &
      +CMGCR(NY,NX)+CMGSR(NY,NX)+CNACR(NY,NX)+CNASR(NY,NX) &
      +CKASR(NY,NX)+CC0PR(NY,NX)) &
      +3.0_r8*(CAL2R(NY,NX)+CFE2R(NY,NX)+CCAHR(NY,NX)+CMGHR(NY,NX) &
      +CF1PR(NY,NX)+CC1PR(NY,NX)+CM1PR(NY,NX)) &
      +4.0_r8*(CAL3R(NY,NX)+CFE3R(NY,NX)+CH3PR(NY,NX)+CF2PR(NY,NX) &
      +CC2PR(NY,NX)) &
      +5.0_r8*(CAL4R(NY,NX)+CFE4R(NY,NX)))
    SII=PRECI(NY,NX)*(CALQ(I,NY,NX)+CFEQ(I,NY,NX)+CHYQ(I,NY,NX) &
      +CCAQ(I,NY,NX)+CMGQ(I,NY,NX)+CNAQ(I,NY,NX)+CKAQ(I,NY,NX) &
      +COHQ(I,NY,NX)+CSOQ(I,NY,NX)+CCLQ(I,NY,NX)+CC3Q(I,NY,NX) &
      +CH0PQ(I,NY,NX) &
      +2.0_r8*(CHCQ(I,NY,NX)+CAL1Q(I,NY,NX)+CALSQ(I,NY,NX) &
      +CFE1Q(I,NY,NX)+CFESQ(I,NY,NX)+CCAOQ(I,NY,NX)+CCACQ(I,NY,NX) &
      +CCASQ(I,NY,NX)+CMGOQ(I,NY,NX)+CMGCQ(I,NY,NX)+CMGSQ(I,NY,NX) &
      +CNACQ(I,NY,NX)+CNASQ(I,NY,NX)+CKASQ(I,NY,NX)+CC0PQ(I,NY,NX)) &
      +3.0_r8*(CAL2Q(I,NY,NX)+CFE2Q(I,NY,NX)+CCAHQ(I,NY,NX) &
      +CMGHQ(I,NY,NX)+CF1PQ(I,NY,NX)+CC1PQ(I,NY,NX)+CM1PQ(I,NY,NX)) &
      +4.0_r8*(CAL3Q(I,NY,NX)+CFE3Q(I,NY,NX) &
      +CH3PQ(I,NY,NX)+CF2PQ(I,NY,NX)+CC2PQ(I,NY,NX)) &
      +5.0_r8*(CAL4Q(I,NY,NX)+CFE4Q(I,NY,NX)))
    TIONIN=TIONIN+SIR+SII
    !     WRITE(*,3338)'SIR',I,J,SIR,PRECQ(NY,NX)
    !    2,CALR(NY,NX),CFER(NY,NX),CHYR(NY,NX),CCAR(NY,NX)
    !    2,CMGR(NY,NX),CNAR(NY,NX),CKAR(NY,NX),COHR(NY,NX),CSOR(NY,NX)
    !    3,CCLR(NY,NX),CC3R(NY,NX),CH0PR(NY,NX)
    !    4,CHCR(NY,NX),CAL1R(NY,NX),CALSR(NY,NX),CFE1R(NY,NX)
    !    5,CFESR(NY,NX),CCAOR(NY,NX),CCACR(NY,NX),CCASR(NY,NX),CMGOR(NY,NX)
    !    6,CMGCR(NY,NX),CMGSR(NY,NX),CNACR(NY,NX),CNASR(NY,NX)
    !    7,CKASR(NY,NX),CC0PR(NY,NX)
    !    8,CAL2R(NY,NX),CFE2R(NY,NX),CCAHR(NY,NX),CMGHR(NY,NX)
    !    9,CF1PR(NY,NX),CC1PR(NY,NX),CM1PR(NY,NX)
    !    1,CAL3R(NY,NX),CFE3R(NY,NX),CH3PR(NY,NX),CF2PR(NY,NX)
    !    2,CC2PR(NY,NX),CAL4R(NY,NX),CFE4R(NY,NX)
    !
    !     SUBSURFACE BOUNDARY SALT FLUXES FROM SUBSURFACE IRRIGATION
    !
    SBU=-PRECU(NY,NX)*(CALQ(I,NY,NX)+CFEQ(I,NY,NX)+CHYQ(I,NY,NX) &
      +CCAQ(I,NY,NX)+CMGQ(I,NY,NX)+CNAQ(I,NY,NX)+CKAQ(I,NY,NX) &
      +COHQ(I,NY,NX)+CSOQ(I,NY,NX)+CCLQ(I,NY,NX)+CC3Q(I,NY,NX) &
      +CH0PQ(I,NY,NX) &
      +2.0_r8*(CHCQ(I,NY,NX)+CAL1Q(I,NY,NX)+CALSQ(I,NY,NX) &
      +CFE1Q(I,NY,NX)+CFESQ(I,NY,NX)+CCAOQ(I,NY,NX)+CCACQ(I,NY,NX) &
      +CCASQ(I,NY,NX)+CMGOQ(I,NY,NX)+CMGCQ(I,NY,NX)+CMGSQ(I,NY,NX) &
      +CNACQ(I,NY,NX)+CNASQ(I,NY,NX)+CKASQ(I,NY,NX)+CC0PQ(I,NY,NX)) &
      +3.0_r8*(CAL2Q(I,NY,NX)+CFE2Q(I,NY,NX)+CCAHQ(I,NY,NX)+CMGHQ(I,NY,NX) &
      +CF1PQ(I,NY,NX)+CC1PQ(I,NY,NX)+CM1PQ(I,NY,NX)) &
      +4.0_r8*(CAL3Q(I,NY,NX)+CFE3Q(I,NY,NX) &
      +CH3PQ(I,NY,NX)+CF2PQ(I,NY,NX)+CC2PQ(I,NY,NX)) &
      +5.0_r8*(CAL4Q(I,NY,NX)+CFE4Q(I,NY,NX)))
    TIONOU=TIONOU+SBU
  ENDIF
  !
  ! GAS EXCHANGE FROM SURFACE VOLATILIZATION-DISSOLUTION
  !
  DO 9680 K=0,2
    OQC(K,0,NY,NX)=OQC(K,0,NY,NX)+XOCFLS(K,3,0,NY,NX)
    OQN(K,0,NY,NX)=OQN(K,0,NY,NX)+XONFLS(K,3,0,NY,NX)
    OQP(K,0,NY,NX)=OQP(K,0,NY,NX)+XOPFLS(K,3,0,NY,NX)
    OQA(K,0,NY,NX)=OQA(K,0,NY,NX)+XOAFLS(K,3,0,NY,NX)
9680  CONTINUE
  CO2S(0,NY,NX)=CO2S(0,NY,NX)+XCODFR(NY,NX)+XCOFLS(3,0,NY,NX) &
    +XCODFG(0,NY,NX)-RCO2O(0,NY,NX)
  CH4S(0,NY,NX)=CH4S(0,NY,NX)+XCHDFR(NY,NX)+XCHFLS(3,0,NY,NX) &
    +XCHDFG(0,NY,NX)-RCH4O(0,NY,NX)
  OXYS(0,NY,NX)=OXYS(0,NY,NX)+XOXDFR(NY,NX)+XOXFLS(3,0,NY,NX) &
    +XOXDFG(0,NY,NX)-RUPOXO(0,NY,NX)
  Z2GS(0,NY,NX)=Z2GS(0,NY,NX)+XNGDFR(NY,NX)+XNGFLS(3,0,NY,NX) &
    +XNGDFG(0,NY,NX)-RN2G(0,NY,NX)-XN2GS(0,NY,NX)
  Z2OS(0,NY,NX)=Z2OS(0,NY,NX)+XN2DFR(NY,NX)+XN2FLS(3,0,NY,NX) &
    +XN2DFG(0,NY,NX)-RN2O(0,NY,NX)
  H2GS(0,NY,NX)=H2GS(0,NY,NX)+XHGDFR(NY,NX)+XHGFLS(3,0,NY,NX) &
    +XHGDFG(0,NY,NX)-RH2GO(0,NY,NX)
  ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+XN4FLW(3,0,NY,NX) &
    +XNH4S(0,NY,NX)+TRN4S(0,NY,NX)
  ZNH3S(0,NY,NX)=ZNH3S(0,NY,NX)+XN3DFR(NY,NX)+XN3FLW(3,0,NY,NX) &
    +XN3DFG(0,NY,NX)+TRN3S(0,NY,NX)
  ZNO3S(0,NY,NX)=ZNO3S(0,NY,NX)+XNOFLW(3,0,NY,NX) &
    +XNO3S(0,NY,NX)+TRNO3(0,NY,NX)
  ZNO2S(0,NY,NX)=ZNO2S(0,NY,NX)+XNXFLS(3,0,NY,NX)+XNO2S(0,NY,NX)
  H1PO4(0,NY,NX)=H1PO4(0,NY,NX)+TRH1P(0,NY,NX)+XH1PFS(3,0,NY,NX) &
    +XH1PS(0,NY,NX)
  H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+TRH2P(0,NY,NX)+XH2PFS(3,0,NY,NX) &
    +XH2PS(0,NY,NX)
  CO2S(NU(NY,NX),NY,NX)=CO2S(NU(NY,NX),NY,NX)+XCODFS(NY,NX)
  CH4S(NU(NY,NX),NY,NX)=CH4S(NU(NY,NX),NY,NX)+XCHDFS(NY,NX)
  OXYS(NU(NY,NX),NY,NX)=OXYS(NU(NY,NX),NY,NX)+XOXDFS(NY,NX)
  Z2GS(NU(NY,NX),NY,NX)=Z2GS(NU(NY,NX),NY,NX)+XNGDFS(NY,NX)
  Z2OS(NU(NY,NX),NY,NX)=Z2OS(NU(NY,NX),NY,NX)+XN2DFS(NY,NX)
  ZNH3S(NU(NY,NX),NY,NX)=ZNH3S(NU(NY,NX),NY,NX)+XN3DFS(NY,NX)
  ZNH3B(NU(NY,NX),NY,NX)=ZNH3B(NU(NY,NX),NY,NX)+XNBDFS(NY,NX)
  H2GS(NU(NY,NX),NY,NX)=H2GS(NU(NY,NX),NY,NX)+XHGDFS(NY,NX)
  THRE(NY,NX)=THRE(NY,NX)+RCO2O(0,NY,NX)+RCH4O(0,NY,NX)
  UN2GG(NY,NX)=UN2GG(NY,NX)+RN2G(0,NY,NX)
  HN2GG(NY,NX)=HN2GG(NY,NX)+RN2G(0,NY,NX)
  ROXYF(0,NY,NX)=XOXDFG(0,NY,NX)
  RCO2F(0,NY,NX)=XCODFG(0,NY,NX)
  RCH4F(0,NY,NX)=XCHDFG(0,NY,NX)
  ROXYL(0,NY,NX)=XOXDFR(NY,NX)+XOXFLS(3,0,NY,NX) &
    -(FLQRQ(NY,NX)*COXR(NY,NX)+FLQRI(NY,NX)*COXQ(NY,NX))
  RCH4L(0,NY,NX)=XCHDFR(NY,NX)+XCHFLS(3,0,NY,NX) &
    -(FLQRQ(NY,NX)*CCHR(NY,NX)+FLQRI(NY,NX)*CCHQ(NY,NX))
  ROXYL(NU(NY,NX),NY,NX)=ROXYL(NU(NY,NX),NY,NX)+XOXDFS(NY,NX)
  RCH4L(NU(NY,NX),NY,NX)=RCH4L(NU(NY,NX),NY,NX)+XCHDFS(NY,NX)
  !     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
  !     WRITE(*,6644)'OXYS0',I,J,NX,NY,NU(NY,NX)
  !    2,OXYS(0,NY,NX),XOXDFR(NY,NX),XOXFLS(3,0,NY,NX)
  !    2,XOXDFG(0,NY,NX),RUPOXO(0,NY,NX)
  !     WRITE(*,6644)'CO2',I,J,NX,NY,NU(NY,NX)
  !    2,CO2S(NU(NY,NX),NY,NX),HCO2G(NY,NX)
  !    2,CI,XCODFS(NY,NX),XCOFLG(3,NU(NY,NX),NY,NX),TCO2Z(NY,NX)
  !    3,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX)
  !    4,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX)
  !    5,CO2S(0,NY,NX),XCODFR(NY,NX),XCOFLS(3,0,NY,NX)
  !    2,XCODFG(0,NY,NX),RCO2O(0,NY,NX)
  !    5,XCODFG(0,NY,NX),XCODFR(NY,NX),VOLP(0,NY,NX)
  !    6,VOLP(NU(NY,NX),NY,NX)
  !     WRITE(*,6644)'OXYS1',I,J,NX,NY,NU(NY,NX)
  !    2,OXYS(NU(NY,NX),NY,NX),HOXYG(NY,NX)
  !    2,XOXDFS(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX),TOXYZ(NY,NX)
  !    5,XOXDFG(NU(NY,NX),NY,NX),VOLP(NU(NY,NX),NY,NX)
  !    3,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX)
  !    4,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX)
  !     WRITE(*,6644)'CH4',I,J,NX,NY,NU(NY,NX),CH,XCHDFS(NY,NX)
  !    2,XCHFLG(3,NU(NY,NX),NY,NX),TCH4Z(NY,NX),FLQGQ(NY,NX)
  !    3,FLQRQ(NY,NX),FLQGI(NY,NX),FLQRI(NY,NX),CCHR(NY,NX),CCHQ(NY,NX)
  !    4,XCHDFG(0,NY,NX),XCHDFR(NY,NX),CH4S(NU(NY,NX),NY,NX)
!6644  FORMAT(A8,5I4,30E12.4)
  !     ENDIF
  !
  !     SURFACE LITTER ION EXCHANGE AND PRECIPITATION
  !
  XN4(0,NY,NX)=XN4(0,NY,NX)+TRXN4(0,NY,NX)
  XOH0(0,NY,NX)=XOH0(0,NY,NX)+TRXH0(0,NY,NX)
  XOH1(0,NY,NX)=XOH1(0,NY,NX)+TRXH1(0,NY,NX)
  XOH2(0,NY,NX)=XOH2(0,NY,NX)+TRXH2(0,NY,NX)
  XH1P(0,NY,NX)=XH1P(0,NY,NX)+TRX1P(0,NY,NX)
  XH2P(0,NY,NX)=XH2P(0,NY,NX)+TRX2P(0,NY,NX)
  PALPO(0,NY,NX)=PALPO(0,NY,NX)+TRALPO(0,NY,NX)
  PFEPO(0,NY,NX)=PFEPO(0,NY,NX)+TRFEPO(0,NY,NX)
  PCAPD(0,NY,NX)=PCAPD(0,NY,NX)+TRCAPD(0,NY,NX)
  PCAPH(0,NY,NX)=PCAPH(0,NY,NX)+TRCAPH(0,NY,NX)
  PCAPM(0,NY,NX)=PCAPM(0,NY,NX)+TRCAPM(0,NY,NX)
  !
  !     SURFACE LITTER OUTPUTS
  !
  !     IF(I.GE.350)THEN
  !     WRITE(*,1119)'CO2S0',I,J,NX,NY,CO2S(0,NY,NX),XCODFS(NY,NX)
  !    2,XCODFR(NY,NX),XCOFLS(3,0,NY,NX),XCODFG(0,NY,NX),RCO2O(0,NY,NX)
  !    3,VOLT(0,NY,NX),CVRD(NY,NX)
  !     WRITE(*,1119)'CH4S0',I,J,NX,NY,CH4S(0,NY,NX),XCHDFS(NY,NX)
  !    2,XCHDFR(NY,NX),XCHFLS(3,0,NY,NX),XCHDFG(0,NY,NX),RCH4O(0,NY,NX)
  !    3,RCH4L(0,NY,NX)
  !     WRITE(*,1119)'OXYS0',I,J,NX,NY,OXYS(0,NY,NX),XOXDFR(NY,NX)
  !    2,XOXFLS(3,0,NY,NX),XOXDFG(0,NY,NX),RUPOXO(0,NY,NX)
  !    3,ROXYL(0,NY,NX),TOXQRS(NY,NX),COXYS(0,NY,NX)
!1119  FORMAT(A8,4I4,12E12.4)
  !     ENDIF
  !     IF(NX.EQ.5)THEN
  !     WRITE(*,5533)'NH30',I,J,NX,NY,ZNH4S(0,NY,NX),XN4FLW(3,0,NY,NX)
  !    2,XNH4S(0,NY,NX),XN3FLW(3,0,NY,NX),TRN4S(0,NY,NX),TRXN4(0,NY,NX)
  !    3,ZNH3S(0,NY,NX),TRN3S(0,NY,NX),XN3DFG(0,NY,NX),XN3DFR(NY,NX)
  !     4,ZNHUFA(0,NY,NX),XNO2S(0,NY,NX),XN4(0,NY,NX)*14.0
  !     WRITE(*,5533)'ZNO3S0',I,J,NX,NY,ZNO3S(0,NY,NX),XNOFLW(3,0,NY,NX)
  !    2,XNO3S(0,NY,NX),TRNO3(0,NY,NX),ZNO2S(0,NY,NX),XNXFLS(3,0,NY,NX)
  !    3,XNO2S(0,NY,NX)
  !     WRITE(*,5533)'H2PO40',I,J,NX,NY,H2PO4(0,NY,NX)
  !    2,XH2PFS(3,0,NY,NX),XH2PS(0,NY,NX),TRH2P(0,NY,NX)
!5533  FORMAT(A8,4I4,20F12.4)
  !     ENDIF
  !     WRITE(*,5544)'HP140',I,J,NX,NY,H1PO4(0,NY,NX)
  !    2,XH1P(0,NY,NX),TRH1P(0,NY,NX),XH1PFS(3,0,NY,NX)
  !    2,XH1PS(0,NY,NX),TP1QRS(NY,NX)
  !     WRITE(*,5544)'HP240',I,J,NX,NY,H2PO4(0,NY,NX)
  !    2,XH2P(0,NY,NX),TRH2P(0,NY,NX),XH2PFS(3,0,NY,NX)
  !    2,XH2PS(0,NY,NX),TPOQRS(NY,NX)
!5544  FORMAT(A8,4I4,40E12.4)
  end subroutine HandleSurfaceBoundary
!------------------------------------------------------------------------------------------

  subroutine OverlandFlow(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K
  !     begin_execution
  !
  IF(ABS(TQR(NY,NX)).GT.ZEROS(NY,NX))THEN
    !
    !     DOC, DON, DOP
    !
    DO 8570 K=0,2
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)+TOCQRS(K,NY,NX)
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)+TONQRS(K,NY,NX)
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)+TOPQRS(K,NY,NX)
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)+TOAQRS(K,NY,NX)
!     IF(NX.EQ.1.AND.NY.EQ.6)THEN
!     WRITE(*,2626)'OQCR',I,J,NX,NY,K,OQC(K,0,NY,NX)
!    2,TOCQRS(K,NY,NX),OQN(K,0,NY,NX),TONQRS(K,NY,NX)
!2626  FORMAT(A8,5I4,20E12.4)
!     ENDIF
8570  CONTINUE
!
    !    SOLUTES
!
    CO2S(0,NY,NX)=CO2S(0,NY,NX)+TCOQRS(NY,NX)
    CH4S(0,NY,NX)=CH4S(0,NY,NX)+TCHQRS(NY,NX)
    OXYS(0,NY,NX)=OXYS(0,NY,NX)+TOXQRS(NY,NX)
    Z2GS(0,NY,NX)=Z2GS(0,NY,NX)+TNGQRS(NY,NX)
    Z2OS(0,NY,NX)=Z2OS(0,NY,NX)+TN2QRS(NY,NX)
    H2GS(0,NY,NX)=H2GS(0,NY,NX)+THGQRS(NY,NX)
    ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+TN4QRS(NY,NX)
    ZNH3S(0,NY,NX)=ZNH3S(0,NY,NX)+TN3QRS(NY,NX)
    ZNO3S(0,NY,NX)=ZNO3S(0,NY,NX)+TNOQRS(NY,NX)
    ZNO2S(0,NY,NX)=ZNO2S(0,NY,NX)+TNXQRS(NY,NX)
    H1PO4(0,NY,NX)=H1PO4(0,NY,NX)+TP1QRS(NY,NX)
    H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+TPOQRS(NY,NX)
    IF(ISALTG.NE.0)THEN
      ZAL(0,NY,NX)=ZAL(0,NY,NX)+TQRAL(NY,NX)
      ZFE(0,NY,NX)=ZFE(0,NY,NX)+TQRFE(NY,NX)
      ZHY(0,NY,NX)=ZHY(0,NY,NX)+TQRHY(NY,NX)
      ZCA(0,NY,NX)=ZCA(0,NY,NX)+TQRCA(NY,NX)
      ZMG(0,NY,NX)=ZMG(0,NY,NX)+TQRMG(NY,NX)
      ZNA(0,NY,NX)=ZNA(0,NY,NX)+TQRNA(NY,NX)
      ZKA(0,NY,NX)=ZKA(0,NY,NX)+TQRKA(NY,NX)
      ZOH(0,NY,NX)=ZOH(0,NY,NX)+TQROH(NY,NX)
      ZSO4(0,NY,NX)=ZSO4(0,NY,NX)+TQRSO(NY,NX)
      ZCL(0,NY,NX)=ZCL(0,NY,NX)+TQRCL(NY,NX)
      ZCO3(0,NY,NX)=ZCO3(0,NY,NX)+TQRC3(NY,NX)
      ZHCO3(0,NY,NX)=ZHCO3(0,NY,NX)+TQRHC(NY,NX)
      ZALOH1(0,NY,NX)=ZALOH1(0,NY,NX)+TQRAL1(NY,NX)
      ZALOH2(0,NY,NX)=ZALOH2(0,NY,NX)+TQRAL2(NY,NX)
      ZALOH3(0,NY,NX)=ZALOH3(0,NY,NX)+TQRAL3(NY,NX)
      ZALOH4(0,NY,NX)=ZALOH4(0,NY,NX)+TQRAL4(NY,NX)
      ZALS(0,NY,NX)=ZALS(0,NY,NX)+TQRALS(NY,NX)
      ZFEOH1(0,NY,NX)=ZFEOH1(0,NY,NX)+TQRFE1(NY,NX)
      ZFEOH2(0,NY,NX)=ZFEOH2(0,NY,NX)+TQRFE2(NY,NX)
      ZFEOH3(0,NY,NX)=ZFEOH3(0,NY,NX)+TQRFE3(NY,NX)
      ZFEOH4(0,NY,NX)=ZFEOH4(0,NY,NX)+TQRFE4(NY,NX)
      ZFES(0,NY,NX)=ZFES(0,NY,NX)+TQRFES(NY,NX)
      ZCAO(0,NY,NX)=ZCAO(0,NY,NX)+TQRCAO(NY,NX)
      ZCAC(0,NY,NX)=ZCAC(0,NY,NX)+TQRCAC(NY,NX)
      ZCAH(0,NY,NX)=ZCAH(0,NY,NX)+TQRCAH(NY,NX)
      ZCAS(0,NY,NX)=ZCAS(0,NY,NX)+TQRCAS(NY,NX)
      ZMGO(0,NY,NX)=ZMGO(0,NY,NX)+TQRMGO(NY,NX)
      ZMGC(0,NY,NX)=ZMGC(0,NY,NX)+TQRMGC(NY,NX)
      ZMGH(0,NY,NX)=ZMGH(0,NY,NX)+TQRMGH(NY,NX)
      ZMGS(0,NY,NX)=ZMGS(0,NY,NX)+TQRMGS(NY,NX)
      ZNAC(0,NY,NX)=ZNAC(0,NY,NX)+TQRNAC(NY,NX)
      ZNAS(0,NY,NX)=ZNAS(0,NY,NX)+TQRNAS(NY,NX)
      ZKAS(0,NY,NX)=ZKAS(0,NY,NX)+TQRKAS(NY,NX)
      H0PO4(0,NY,NX)=H0PO4(0,NY,NX)+TQRH0P(NY,NX)
      H3PO4(0,NY,NX)=H3PO4(0,NY,NX)+TQRH3P(NY,NX)
      ZFE1P(0,NY,NX)=ZFE1P(0,NY,NX)+TQRF1P(NY,NX)
      ZFE2P(0,NY,NX)=ZFE2P(0,NY,NX)+TQRF2P(NY,NX)
      ZCA0P(0,NY,NX)=ZCA0P(0,NY,NX)+TQRC0P(NY,NX)
      ZCA1P(0,NY,NX)=ZCA1P(0,NY,NX)+TQRC1P(NY,NX)
      ZCA2P(0,NY,NX)=ZCA2P(0,NY,NX)+TQRC2P(NY,NX)
      ZMG1P(0,NY,NX)=ZMG1P(0,NY,NX)+TQRM1P(NY,NX)
    ENDIF
  ENDIF
  end subroutine OverlandFlow
!------------------------------------------------------------------------------------------

  subroutine ChemicalOverlandFlow(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,NO,M,NGL
  ! begin_execution
  !
  ! INTERNAL SURFACE SEDIMENT TRANSPORT
  !
  IF((IERSNG.EQ.1.OR.IERSNG.EQ.3) &
    .AND.ABS(TSEDER(NY,NX)).GT.ZEROS(NY,NX))THEN
    TSED(NY,NX)=TSED(NY,NX)+TSEDER(NY,NX)
    !
    !     SOIL MINERAL FRACTIONS
    !
    SAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)+TSANER(NY,NX)
    SILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)+TSILER(NY,NX)
    CLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)+TCLAER(NY,NX)
    XCEC(NU(NY,NX),NY,NX)=XCEC(NU(NY,NX),NY,NX)+TCECER(NY,NX)
    XAEC(NU(NY,NX),NY,NX)=XAEC(NU(NY,NX),NY,NX)+TAECER(NY,NX)
    !
    !     FERTILIZER POOLS
!
    ZNH4FA(NU(NY,NX),NY,NX)=ZNH4FA(NU(NY,NX),NY,NX)+TNH4ER(NY,NX)
    ZNH3FA(NU(NY,NX),NY,NX)=ZNH3FA(NU(NY,NX),NY,NX)+TNH3ER(NY,NX)
    ZNHUFA(NU(NY,NX),NY,NX)=ZNHUFA(NU(NY,NX),NY,NX)+TNHUER(NY,NX)
    ZNO3FA(NU(NY,NX),NY,NX)=ZNO3FA(NU(NY,NX),NY,NX)+TNO3ER(NY,NX)
    ZNH4FB(NU(NY,NX),NY,NX)=ZNH4FB(NU(NY,NX),NY,NX)+TNH4EB(NY,NX)
    ZNH3FB(NU(NY,NX),NY,NX)=ZNH3FB(NU(NY,NX),NY,NX)+TNH3EB(NY,NX)
    ZNHUFB(NU(NY,NX),NY,NX)=ZNHUFB(NU(NY,NX),NY,NX)+TNHUEB(NY,NX)
    ZNO3FB(NU(NY,NX),NY,NX)=ZNO3FB(NU(NY,NX),NY,NX)+TNO3EB(NY,NX)
!
    !   EXCHANGEABLE CATIONS AND ANIONS
!
    XN4(NU(NY,NX),NY,NX)=XN4(NU(NY,NX),NY,NX)+TN4ER(NY,NX)
    XNB(NU(NY,NX),NY,NX)=XNB(NU(NY,NX),NY,NX)+TNBER(NY,NX)
    XHY(NU(NY,NX),NY,NX)=XHY(NU(NY,NX),NY,NX)+THYER(NY,NX)
    XAL(NU(NY,NX),NY,NX)=XAL(NU(NY,NX),NY,NX)+TALER(NY,NX)
    XFE(NU(NY,NX),NY,NX)=XFE(NU(NY,NX),NY,NX)+TFEER(NY,NX)
    XCA(NU(NY,NX),NY,NX)=XCA(NU(NY,NX),NY,NX)+TCAER(NY,NX)
    XMG(NU(NY,NX),NY,NX)=XMG(NU(NY,NX),NY,NX)+TMGER(NY,NX)
    XNA(NU(NY,NX),NY,NX)=XNA(NU(NY,NX),NY,NX)+TNAER(NY,NX)
    XKA(NU(NY,NX),NY,NX)=XKA(NU(NY,NX),NY,NX)+TKAER(NY,NX)
    XHC(NU(NY,NX),NY,NX)=XHC(NU(NY,NX),NY,NX)+THCER(NY,NX)
    XALO2(NU(NY,NX),NY,NX)=XALO2(NU(NY,NX),NY,NX)+TAL2ER(NY,NX)
    XFEO2(NU(NY,NX),NY,NX)=XFEO2(NU(NY,NX),NY,NX)+TFE2ER(NY,NX)
    XOH0(NU(NY,NX),NY,NX)=XOH0(NU(NY,NX),NY,NX)+TOH0ER(NY,NX)
    XOH1(NU(NY,NX),NY,NX)=XOH1(NU(NY,NX),NY,NX)+TOH1ER(NY,NX)
    XOH2(NU(NY,NX),NY,NX)=XOH2(NU(NY,NX),NY,NX)+TOH2ER(NY,NX)
    XH1P(NU(NY,NX),NY,NX)=XH1P(NU(NY,NX),NY,NX)+TH1PER(NY,NX)
    XH2P(NU(NY,NX),NY,NX)=XH2P(NU(NY,NX),NY,NX)+TH2PER(NY,NX)
    XOH0B(NU(NY,NX),NY,NX)=XOH0B(NU(NY,NX),NY,NX)+TOH0EB(NY,NX)
    XOH1B(NU(NY,NX),NY,NX)=XOH1B(NU(NY,NX),NY,NX)+TOH1EB(NY,NX)
    XOH2B(NU(NY,NX),NY,NX)=XOH2B(NU(NY,NX),NY,NX)+TOH2EB(NY,NX)
    XH1PB(NU(NY,NX),NY,NX)=XH1PB(NU(NY,NX),NY,NX)+TH1PEB(NY,NX)
    XH2PB(NU(NY,NX),NY,NX)=XH2PB(NU(NY,NX),NY,NX)+TH2PEB(NY,NX)
!
    !     PRECIPITATES
!
    PALOH(NU(NY,NX),NY,NX)=PALOH(NU(NY,NX),NY,NX)+TALOER(NY,NX)
    PFEOH(NU(NY,NX),NY,NX)=PFEOH(NU(NY,NX),NY,NX)+TFEOER(NY,NX)
    PCACO(NU(NY,NX),NY,NX)=PCACO(NU(NY,NX),NY,NX)+TCACER(NY,NX)
    PCASO(NU(NY,NX),NY,NX)=PCASO(NU(NY,NX),NY,NX)+TCASER(NY,NX)
    PALPO(NU(NY,NX),NY,NX)=PALPO(NU(NY,NX),NY,NX)+TALPER(NY,NX)
    PFEPO(NU(NY,NX),NY,NX)=PFEPO(NU(NY,NX),NY,NX)+TFEPER(NY,NX)
    PCAPD(NU(NY,NX),NY,NX)=PCAPD(NU(NY,NX),NY,NX)+TCPDER(NY,NX)
    PCAPH(NU(NY,NX),NY,NX)=PCAPH(NU(NY,NX),NY,NX)+TCPHER(NY,NX)
    PCAPM(NU(NY,NX),NY,NX)=PCAPM(NU(NY,NX),NY,NX)+TCPMER(NY,NX)
    PALPB(NU(NY,NX),NY,NX)=PALPB(NU(NY,NX),NY,NX)+TALPEB(NY,NX)
    PFEPB(NU(NY,NX),NY,NX)=PFEPB(NU(NY,NX),NY,NX)+TFEPEB(NY,NX)
    PCPDB(NU(NY,NX),NY,NX)=PCPDB(NU(NY,NX),NY,NX)+TCPDEB(NY,NX)
    PCPHB(NU(NY,NX),NY,NX)=PCPHB(NU(NY,NX),NY,NX)+TCPHEB(NY,NX)
    PCPMB(NU(NY,NX),NY,NX)=PCPMB(NU(NY,NX),NY,NX)+TCPMEB(NY,NX)
!
    !   ORGANIC CONSTITUENTS
!
    DORGP=0.0_r8
    DO 9280 K=0,jcplx1
      DO  NO=1,7
        DO  M=1,3
          DO NGL=1,JG
            OMC(M,NGL,NO,K,NU(NY,NX),NY,NX)=OMC(M,NGL,NO,K,NU(NY,NX),NY,NX)+TOMCER(M,NGL,NO,K,NY,NX)
            OMN(M,NGL,NO,K,NU(NY,NX),NY,NX)=OMN(M,NGL,NO,K,NU(NY,NX),NY,NX)+TOMNER(M,NGL,NO,K,NY,NX)
            OMP(M,NGL,NO,K,NU(NY,NX),NY,NX)=OMP(M,NGL,NO,K,NU(NY,NX),NY,NX)+TOMPER(M,NGL,NO,K,NY,NX)
            DORGE(NY,NX)=DORGE(NY,NX)+TOMCER(M,NGL,NO,K,NY,NX)
            DORGP=DORGP+TOMPER(M,NGL,NO,K,NY,NX)
          enddo
        enddo
      enddo
9280  CONTINUE
      DO  NO=1,7
        DO  M=1,3
          DO NGL=1,JG
            OMCff(M,NGL,NO,NU(NY,NX),NY,NX)=OMCff(M,NGL,NO,NU(NY,NX),NY,NX)+TOMCERff(M,NGL,NO,NY,NX)
            OMNff(M,NGL,NO,NU(NY,NX),NY,NX)=OMNff(M,NGL,NO,NU(NY,NX),NY,NX)+TOMNERff(M,NGL,NO,NY,NX)
            OMPff(M,NGL,NO,NU(NY,NX),NY,NX)=OMPff(M,NGL,NO,NU(NY,NX),NY,NX)+TOMPERff(M,NGL,NO,NY,NX)
            DORGE(NY,NX)=DORGE(NY,NX)+TOMCERff(M,NGL,NO,NY,NX)
            DORGP=DORGP+TOMPER(M,NGL,NO,K,NY,NX)
          enddo
        enddo
      enddo

    DO 9275 K=0,jcplx1
      DO 9270 M=1,2
        ORC(M,K,NU(NY,NX),NY,NX)=ORC(M,K,NU(NY,NX),NY,NX)+TORCER(M,K,NY,NX)
        ORN(M,K,NU(NY,NX),NY,NX)=ORN(M,K,NU(NY,NX),NY,NX)+TORNER(M,K,NY,NX)
        ORP(M,K,NU(NY,NX),NY,NX)=ORP(M,K,NU(NY,NX),NY,NX)+TORPER(M,K,NY,NX)
        DORGE(NY,NX)=DORGE(NY,NX)+TORCER(M,K,NY,NX)
        DORGP=DORGP+TORPER(M,K,NY,NX)
9270  CONTINUE
      OHC(K,NU(NY,NX),NY,NX)=OHC(K,NU(NY,NX),NY,NX)+TOHCER(K,NY,NX)
      OHN(K,NU(NY,NX),NY,NX)=OHN(K,NU(NY,NX),NY,NX)+TOHNER(K,NY,NX)
      OHP(K,NU(NY,NX),NY,NX)=OHP(K,NU(NY,NX),NY,NX)+TOHPER(K,NY,NX)
      OHA(K,NU(NY,NX),NY,NX)=OHA(K,NU(NY,NX),NY,NX)+TOHAER(K,NY,NX)
      DORGE(NY,NX)=DORGE(NY,NX)+TOHCER(K,NY,NX)+TOHAER(K,NY,NX)
      DORGP=DORGP+TOHPER(K,NY,NX)
      DO 9265 M=1,jsken
        OSC(M,K,NU(NY,NX),NY,NX)=OSC(M,K,NU(NY,NX),NY,NX)+TOSCER(M,K,NY,NX)
        OSA(M,K,NU(NY,NX),NY,NX)=OSA(M,K,NU(NY,NX),NY,NX)+TOSAER(M,K,NY,NX)
        OSN(M,K,NU(NY,NX),NY,NX)=OSN(M,K,NU(NY,NX),NY,NX)+TOSNER(M,K,NY,NX)
        OSP(M,K,NU(NY,NX),NY,NX)=OSP(M,K,NU(NY,NX),NY,NX)+TOSPER(M,K,NY,NX)
        DORGE(NY,NX)=DORGE(NY,NX)+TOSCER(M,K,NY,NX)
        DORGP=DORGP+TOSPER(M,K,NY,NX)
9265  CONTINUE
9275  CONTINUE
  ENDIF
  end subroutine ChemicalOverlandFlow
!------------------------------------------------------------------------------------------

  subroutine ChemicalBySnowRedistribution(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TQS(NY,NX))>0._r8)THEN
    CO2W(1,NY,NX)=CO2W(1,NY,NX)+TCOQSS(NY,NX)
    CH4W(1,NY,NX)=CH4W(1,NY,NX)+TCHQSS(NY,NX)
    OXYW(1,NY,NX)=OXYW(1,NY,NX)+TOXQSS(NY,NX)
    ZNGW(1,NY,NX)=ZNGW(1,NY,NX)+TNGQSS(NY,NX)
    ZN2W(1,NY,NX)=ZN2W(1,NY,NX)+TN2QSS(NY,NX)
    ZN4W(1,NY,NX)=ZN4W(1,NY,NX)+TN4QSS(NY,NX)
    ZN3W(1,NY,NX)=ZN3W(1,NY,NX)+TN3QSS(NY,NX)
    ZNOW(1,NY,NX)=ZNOW(1,NY,NX)+TNOQSS(NY,NX)
    Z1PW(1,NY,NX)=Z1PW(1,NY,NX)+TP1QSS(NY,NX)
    ZHPW(1,NY,NX)=ZHPW(1,NY,NX)+TPOQSS(NY,NX)
    IF(ISALTG.NE.0)THEN
      ZALW(1,NY,NX)=ZALW(1,NY,NX)+TQSAL(NY,NX)
      ZFEW(1,NY,NX)=ZFEW(1,NY,NX)+TQSFE(NY,NX)
      ZHYW(1,NY,NX)=ZHYW(1,NY,NX)+TQSHY(NY,NX)
      ZCAW(1,NY,NX)=ZCAW(1,NY,NX)+TQSCA(NY,NX)
      ZMGW(1,NY,NX)=ZMGW(1,NY,NX)+TQSMG(NY,NX)
      ZNAW(1,NY,NX)=ZNAW(1,NY,NX)+TQSNA(NY,NX)
      ZKAW(1,NY,NX)=ZKAW(1,NY,NX)+TQSKA(NY,NX)
      ZOHW(1,NY,NX)=ZOHW(1,NY,NX)+TQSOH(NY,NX)
      ZSO4W(1,NY,NX)=ZSO4W(1,NY,NX)+TQSSO(NY,NX)
      ZCLW(1,NY,NX)=ZCLW(1,NY,NX)+TQSCL(NY,NX)
      ZCO3W(1,NY,NX)=ZCO3W(1,NY,NX)+TQSC3(NY,NX)
      ZHCO3W(1,NY,NX)=ZHCO3W(1,NY,NX)+TQSHC(NY,NX)
      ZALH1W(1,NY,NX)=ZALH1W(1,NY,NX)+TQSAL1(NY,NX)
      ZALH2W(1,NY,NX)=ZALH2W(1,NY,NX)+TQSAL2(NY,NX)
      ZALH3W(1,NY,NX)=ZALH3W(1,NY,NX)+TQSAL3(NY,NX)
      ZALH4W(1,NY,NX)=ZALH4W(1,NY,NX)+TQSAL4(NY,NX)
      ZALSW(1,NY,NX)=ZALSW(1,NY,NX)+TQSALS(NY,NX)
      ZFEH1W(1,NY,NX)=ZFEH1W(1,NY,NX)+TQSFE1(NY,NX)
      ZFEH2W(1,NY,NX)=ZFEH2W(1,NY,NX)+TQSFE2(NY,NX)
      ZFEH3W(1,NY,NX)=ZFEH3W(1,NY,NX)+TQSFE3(NY,NX)
      ZFEH4W(1,NY,NX)=ZFEH4W(1,NY,NX)+TQSFE4(NY,NX)
      ZFESW(1,NY,NX)=ZFESW(1,NY,NX)+TQSFES(NY,NX)
      ZCAOW(1,NY,NX)=ZCAOW(1,NY,NX)+TQSCAO(NY,NX)
      ZCACW(1,NY,NX)=ZCACW(1,NY,NX)+TQSCAC(NY,NX)
      ZCAHW(1,NY,NX)=ZCAHW(1,NY,NX)+TQSCAH(NY,NX)
      ZCASW(1,NY,NX)=ZCASW(1,NY,NX)+TQSCAS(NY,NX)
      ZMGOW(1,NY,NX)=ZMGOW(1,NY,NX)+TQSMGO(NY,NX)
      ZMGCW(1,NY,NX)=ZMGCW(1,NY,NX)+TQSMGC(NY,NX)
      ZMGHW(1,NY,NX)=ZMGHW(1,NY,NX)+TQSMGH(NY,NX)
      ZMGSW(1,NY,NX)=ZMGSW(1,NY,NX)+TQSMGS(NY,NX)
      ZNACW(1,NY,NX)=ZNACW(1,NY,NX)+TQSNAC(NY,NX)
      ZNASW(1,NY,NX)=ZNASW(1,NY,NX)+TQSNAS(NY,NX)
      ZKASW(1,NY,NX)=ZKASW(1,NY,NX)+TQSKAS(NY,NX)
      H0PO4W(1,NY,NX)=H0PO4W(1,NY,NX)+TQSH0P(NY,NX)
      H3PO4W(1,NY,NX)=H3PO4W(1,NY,NX)+TQSH3P(NY,NX)
      ZFE1PW(1,NY,NX)=ZFE1PW(1,NY,NX)+TQSF1P(NY,NX)
      ZFE2PW(1,NY,NX)=ZFE2PW(1,NY,NX)+TQSF2P(NY,NX)
      ZCA0PW(1,NY,NX)=ZCA0PW(1,NY,NX)+TQSC0P(NY,NX)
      ZCA1PW(1,NY,NX)=ZCA1PW(1,NY,NX)+TQSC1P(NY,NX)
      ZCA2PW(1,NY,NX)=ZCA2PW(1,NY,NX)+TQSC2P(NY,NX)
      ZMG1PW(1,NY,NX)=ZMG1PW(1,NY,NX)+TQSM1P(NY,NX)
    ENDIF
  ENDIF
  end subroutine ChemicalBySnowRedistribution
!------------------------------------------------------------------------------------------

  subroutine SumupSurfaceChemicalMass(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,N,M,NGL
!     begin_execution
!     TOTAL C,N,P, SALTS IN SURFACE RESIDUE
!
  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  DO 6975 K=0,jcplx1
    RC0(K,NY,NX)=0.0_r8
6975  CONTINUE
  RC0ff(NY,NX)=0.0_r8

  OMCL(0,NY,NX)=0.0_r8
  OMNL(0,NY,NX)=0.0_r8
  DO 6970 K=0,jcplx1
    IF(K.NE.3.AND.K.NE.4)THEN
      !
      ! TOTAL MICROBIAL C,N,P
      !
      DO 6960 N=1,7
        DO  M=1,3
          DO NGL=1,JG
            DC=DC+OMC(M,NGL,N,K,0,NY,NX)
            DN=DN+OMN(M,NGL,N,K,0,NY,NX)
            DP=DP+OMP(M,NGL,N,K,0,NY,NX)
            RC0(K,NY,NX)=RC0(K,NY,NX)+OMC(M,NGL,N,K,0,NY,NX)
            TOMT(NY,NX)=TOMT(NY,NX)+OMC(M,NGL,N,K,0,NY,NX)
            TONT(NY,NX)=TONT(NY,NX)+OMN(M,NGL,N,K,0,NY,NX)
            TOPT(NY,NX)=TOPT(NY,NX)+OMP(M,NGL,N,K,0,NY,NX)
            OMCL(0,NY,NX)=OMCL(0,NY,NX)+OMC(M,NGL,N,K,0,NY,NX)
            OMNL(0,NY,NX)=OMNL(0,NY,NX)+OMN(M,NGL,N,K,0,NY,NX)
          enddo
        enddo
6960  CONTINUE
    ENDIF
6970  CONTINUE

      !
      ! TOTAL MICROBIAL C,N,P
      !
      DO N=1,7
        DO  M=1,3
          DO NGL=1,JG
            DC=DC+OMCff(M,NGL,N,0,NY,NX)
            DN=DN+OMNff(M,NGL,N,0,NY,NX)
            DP=DP+OMPff(M,NGL,N,0,NY,NX)
            RC0ff(NY,NX)=RC0ff(NY,NX)+OMCff(M,NGL,N,0,NY,NX)
            TOMT(NY,NX)=TOMT(NY,NX)+OMCff(M,NGL,N,0,NY,NX)
            TONT(NY,NX)=TONT(NY,NX)+OMNff(M,NGL,N,0,NY,NX)
            TOPT(NY,NX)=TOPT(NY,NX)+OMPff(M,NGL,N,0,NY,NX)
            OMCL(0,NY,NX)=OMCL(0,NY,NX)+OMCff(M,NGL,N,0,NY,NX)
            OMNL(0,NY,NX)=OMNL(0,NY,NX)+OMNff(M,NGL,N,0,NY,NX)
          enddo
        enddo
    ENDDO

  !
  !     TOTAL MICROBIAL RESIDUE C,N,P
  !
  DO 6900 K=0,2
    DO 6940 M=1,2
      DC=DC+ORC(M,K,0,NY,NX)
      DN=DN+ORN(M,K,0,NY,NX)
      DP=DP+ORP(M,K,0,NY,NX)
      RC0(K,NY,NX)=RC0(K,NY,NX)+ORC(M,K,0,NY,NX)
6940  CONTINUE
!
    !     TOTAL DOC, DON, DOP
!
    DC=DC+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX)+OHC(K,0,NY,NX) &
      +OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
    DN=DN+OQN(K,0,NY,NX)+OQNH(K,0,NY,NX)+OHN(K,0,NY,NX)
    DP=DP+OQP(K,0,NY,NX)+OQPH(K,0,NY,NX)+OHP(K,0,NY,NX)
    RC0(K,NY,NX)=RC0(K,NY,NX)+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX) &
      +OHC(K,0,NY,NX)+OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
!
    !     TOTAL PLANT RESIDUE C,N,P
!
    DO 6930 M=1,4
      DC=DC+OSC(M,K,0,NY,NX)
      DN=DN+OSN(M,K,0,NY,NX)
      DP=DP+OSP(M,K,0,NY,NX)
      RC0(K,NY,NX)=RC0(K,NY,NX)+OSC(M,K,0,NY,NX)
6930  CONTINUE
6900  CONTINUE
  ORGC(0,NY,NX)=DC
  ORGN(0,NY,NX)=DN
  ORGR(0,NY,NX)=DC
  TLRSDC=TLRSDC+DC
  URSDC(NY,NX)=URSDC(NY,NX)+DC
  TLRSDN=TLRSDN+DN
  URSDN(NY,NX)=URSDN(NY,NX)+DN
  TLRSDP=TLRSDP+DP
  URSDP(NY,NX)=URSDP(NY,NX)+DP
  WS=TVOLWC(NY,NX)+TVOLWP(NY,NX)+VOLW(0,NY,NX)+VOLI(0,NY,NX)*DENSI
  VOLWSO=VOLWSO+WS
  UVOLW(NY,NX)=UVOLW(NY,NX)+WS
  HEATSO=HEATSO+TENGYC(NY,NX)
  CS=CO2S(0,NY,NX)+CH4S(0,NY,NX)
  TLCO2G=TLCO2G+CS
  UCO2S(NY,NX)=UCO2S(NY,NX)+CS
  HS=H2GS(0,NY,NX)
  TLH2G=TLH2G+HS
  OS=OXYS(0,NY,NX)
  OXYGSO=OXYGSO+OS
  ZG=Z2GS(0,NY,NX)+Z2OS(0,NY,NX)
  TLN2G=TLN2G+ZG
  Z4S=ZNH4S(0,NY,NX)+ZNH3S(0,NY,NX)
  Z4X=14.0*XN4(0,NY,NX)
  Z4F=14.0*(ZNH4FA(0,NY,NX)+ZNHUFA(0,NY,NX)+ZNH3FA(0,NY,NX))
  TLNH4=TLNH4+Z4S+Z4X+Z4F
  UNH4(NY,NX)=UNH4(NY,NX)+Z4S+Z4X

  ZOS=ZNO3S(0,NY,NX)+ZNO2S(0,NY,NX)
  ZOF=14.0*ZNO3FA(0,NY,NX)
  TLNO3=TLNO3+ZOS+ZOF
  UNO3(NY,NX)=UNO3(NY,NX)+ZOS
  POS=H1PO4(0,NY,NX)+H2PO4(0,NY,NX)
  POX=31.0*(XH1P(0,NY,NX)+XH2P(0,NY,NX))
  POP=31.0*(PALPO(0,NY,NX)+PFEPO(0,NY,NX) &
    +PCAPD(0,NY,NX)) &
    +62.0*PCAPM(0,NY,NX) &
    +93.0*PCAPH(0,NY,NX)
  TLPO4=TLPO4+POS+POX+POP
  UPO4(NY,NX)=UPO4(NY,NX)+POX
  UPP4(NY,NX)=UPP4(NY,NX)+POP
  IF(ISALTG.NE.0)THEN
    ZAL(0,NY,NX)=ZAL(0,NY,NX)+XALFLS(3,0,NY,NX)
    ZFE(0,NY,NX)=ZFE(0,NY,NX)+XFEFLS(3,0,NY,NX)
    ZHY(0,NY,NX)=ZHY(0,NY,NX)+XHYFLS(3,0,NY,NX)
    ZCA(0,NY,NX)=ZCA(0,NY,NX)+XCAFLS(3,0,NY,NX)
    ZMG(0,NY,NX)=ZMG(0,NY,NX)+XMGFLS(3,0,NY,NX)
    ZNA(0,NY,NX)=ZNA(0,NY,NX)+XNAFLS(3,0,NY,NX)
    ZKA(0,NY,NX)=ZKA(0,NY,NX)+XKAFLS(3,0,NY,NX)
    ZOH(0,NY,NX)=ZOH(0,NY,NX)+XOHFLS(3,0,NY,NX)

    ZSO4(0,NY,NX)=ZSO4(0,NY,NX)+XSOFLS(3,0,NY,NX)
    ZCL(0,NY,NX)=ZCL(0,NY,NX)+XCLFLS(3,0,NY,NX)
    ZCO3(0,NY,NX)=ZCO3(0,NY,NX)+XC3FLS(3,0,NY,NX)
    ZHCO3(0,NY,NX)=ZHCO3(0,NY,NX)+XHCFLS(3,0,NY,NX)
    ZALOH1(0,NY,NX)=ZALOH1(0,NY,NX)+XAL1FS(3,0,NY,NX)
    ZALOH2(0,NY,NX)=ZALOH2(0,NY,NX)+XAL2FS(3,0,NY,NX)
    ZALOH3(0,NY,NX)=ZALOH3(0,NY,NX)+XAL3FS(3,0,NY,NX)
    ZALOH4(0,NY,NX)=ZALOH4(0,NY,NX)+XAL4FS(3,0,NY,NX)
    ZALS(0,NY,NX)=ZALS(0,NY,NX)+XALSFS(3,0,NY,NX)
    ZFEOH1(0,NY,NX)=ZFEOH1(0,NY,NX)+XFE1FS(3,0,NY,NX)
    ZFEOH2(0,NY,NX)=ZFEOH2(0,NY,NX)+XFE2FS(3,0,NY,NX)
    ZFEOH3(0,NY,NX)=ZFEOH3(0,NY,NX)+XFE3FS(3,0,NY,NX)
    ZFEOH4(0,NY,NX)=ZFEOH4(0,NY,NX)+XFE4FS(3,0,NY,NX)
    ZFES(0,NY,NX)=ZFES(0,NY,NX)+XFESFS(3,0,NY,NX)
    ZCAO(0,NY,NX)=ZCAO(0,NY,NX)+XCAOFS(3,0,NY,NX)
    ZCAC(0,NY,NX)=ZCAC(0,NY,NX)+XCACFS(3,0,NY,NX)
    ZCAH(0,NY,NX)=ZCAH(0,NY,NX)+XCAHFS(3,0,NY,NX)
    ZCAS(0,NY,NX)=ZCAS(0,NY,NX)+XCASFS(3,0,NY,NX)
    ZMGO(0,NY,NX)=ZMGO(0,NY,NX)+XMGOFS(3,0,NY,NX)
    ZMGC(0,NY,NX)=ZMGC(0,NY,NX)+XMGCFS(3,0,NY,NX)
    ZMGH(0,NY,NX)=ZMGH(0,NY,NX)+XMGHFS(3,0,NY,NX)
    ZMGS(0,NY,NX)=ZMGS(0,NY,NX)+XMGSFS(3,0,NY,NX)
    ZNAC(0,NY,NX)=ZNAC(0,NY,NX)+XNACFS(3,0,NY,NX)
    ZNAS(0,NY,NX)=ZNAS(0,NY,NX)+XNASFS(3,0,NY,NX)
    ZKAS(0,NY,NX)=ZKAS(0,NY,NX)+XKASFS(3,0,NY,NX)
    H0PO4(0,NY,NX)=H0PO4(0,NY,NX)+XH0PFS(3,0,NY,NX)
    H3PO4(0,NY,NX)=H3PO4(0,NY,NX)+XH3PFS(3,0,NY,NX)
    ZFE1P(0,NY,NX)=ZFE1P(0,NY,NX)+XF1PFS(3,0,NY,NX)
    ZFE2P(0,NY,NX)=ZFE2P(0,NY,NX)+XF2PFS(3,0,NY,NX)
    ZCA0P(0,NY,NX)=ZCA0P(0,NY,NX)+XC0PFS(3,0,NY,NX)
    ZCA1P(0,NY,NX)=ZCA1P(0,NY,NX)+XC1PFS(3,0,NY,NX)
    ZCA2P(0,NY,NX)=ZCA2P(0,NY,NX)+XC2PFS(3,0,NY,NX)
    ZMG1P(0,NY,NX)=ZMG1P(0,NY,NX)+XM1PFS(3,0,NY,NX)
    PSS=31.0*(H0PO4(0,NY,NX)+H3PO4(0,NY,NX)+ZFE1P(0,NY,NX) &
      +ZFE2P(0,NY,NX)+ZCA0P(0,NY,NX)+ZCA1P(0,NY,NX) &
      +ZCA2P(0,NY,NX)+ZMG1P(0,NY,NX))
    TLPO4=TLPO4+PSS
    SSS=ZAL(0,NY,NX)+ZFE(0,NY,NX)+ZHY(0,NY,NX)+ZCA(0,NY,NX) &
      +ZMG(0,NY,NX)+ZNA(0,NY,NX)+ZKA(0,NY,NX)+ZOH(0,NY,NX) &
      +ZSO4(0,NY,NX)+ZCL(0,NY,NX)+ZCO3(0,NY,NX)+H0PO4(0,NY,NX) &
      +2.0*(ZHCO3(0,NY,NX)+ZALOH1(0,NY,NX)+ZALS(0,NY,NX) &
      +ZFEOH1(0,NY,NX)+ZFES(0,NY,NX)+ZCAO(0,NY,NX)+ZCAC(0,NY,NX) &
      +ZCAS(0,NY,NX)+ZMGO(0,NY,NX)+ZMGC(0,NY,NX)+ZMGS(0,NY,NX) &
      +ZNAC(0,NY,NX)+ZNAS(0,NY,NX)+ZKAS(0,NY,NX)+ZCA0P(0,NY,NX)) &
      +3.0*(ZALOH2(0,NY,NX)+ZFEOH2(0,NY,NX)+ZCAH(0,NY,NX) &
      +ZMGH(0,NY,NX)+ZFE1P(0,NY,NX)+ZCA1P(0,NY,NX)+ZMG1P(0,NY,NX)) &
      +4.0*(ZALOH3(0,NY,NX)+ZFEOH3(0,NY,NX)+H3PO4(0,NY,NX) &
      +ZFE2P(0,NY,NX)+ZCA2P(0,NY,NX)) &
      +5.0*(ZALOH4(0,NY,NX)+ZFEOH4(0,NY,NX))
    TION=TION+SSS
    UION(NY,NX)=UION(NY,NX)+SSS
  ENDIF

  end subroutine SumupSurfaceChemicalMass
!------------------------------------------------------------------------------------------

  subroutine UpdateChemInLayers(NY,NX,LG,VOLISO,DORGC,DVOLI)
  implicit none
  integer, intent(in) :: NY,NX,LG
  real(r8), intent(inout) :: VOLISO
  real(r8),intent(out) :: DORGC(JZ,JY,JX)
  REAL(R8),INTENT(OUT) :: DVOLI(JZ,JY,JX)
  integer  :: L,K,M,N,LL,NGL
  real(r8) :: TKS00,TKSX
  real(r8) :: ECNO,ENGY
  !     begin_execution
  !     UPDATE SOIL LAYER VARIABLES WITH TOTAL FLUXES
  !
  TVHCP=0.0_r8
  TVHCM=0.0_r8
  TVOLW=0.0_r8
  TVOLWH=0.0_r8
  TVOLI=0.0_r8
  TVOLIH=0.0_r8
  TENGY=0.0_r8
  DO 125 L=NU(NY,NX),NL(NY,NX)
    !
    !     WATER, ICE, HEAT, TEMPERATURE
    !
    TKSX=TKS(L,NY,NX)
    VHCPX=VHCP(L,NY,NX)
    VOLWXX=VOLW(L,NY,NX)
    VOLIXX=VOLI(L,NY,NX)
    VOLW(L,NY,NX)=VOLW(L,NY,NX)+TFLW(L,NY,NX)+FINH(L,NY,NX) &
      +TTHAW(L,NY,NX)+TUPWTR(L,NY,NX)+FLU(L,NY,NX)
    VOLWX(L,NY,NX)=VOLWX(L,NY,NX)+TFLWX(L,NY,NX)+FINH(L,NY,NX) &
      +TTHAW(L,NY,NX)+TUPWTR(L,NY,NX)+FLU(L,NY,NX)
    VOLWX(L,NY,NX)=AMIN1(VOLW(L,NY,NX) &
      ,VOLWX(L,NY,NX)+0.01*(VOLW(L,NY,NX)-VOLWX(L,NY,NX)))
    VOLI(L,NY,NX)=VOLI(L,NY,NX)-TTHAW(L,NY,NX)/DENSI
    VOLWH(L,NY,NX)=VOLWH(L,NY,NX)+TFLWH(L,NY,NX)-FINH(L,NY,NX) &
      +TTHAWH(L,NY,NX)
    VOLIH(L,NY,NX)=VOLIH(L,NY,NX)-TTHAWH(L,NY,NX)/DENSI
    DVOLW(L,NY,NX)=VOLW1(L,NY,NX)+VOLWH1(L,NY,NX) &
      -VOLW(L,NY,NX)-VOLWH(L,NY,NX)
    DVOLI(L,NY,NX)=VOLI1(L,NY,NX)+VOLIH1(L,NY,NX) &
      -VOLI(L,NY,NX)-VOLIH(L,NY,NX)
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX)-VOLI(L,NY,NX) &
        +VOLAH(L,NY,NX)-VOLWH(L,NY,NX)-VOLIH(L,NY,NX))
    ELSE
      VOLP(L,NY,NX)=0.0_r8
!     VOLA(L,NY,NX)=VOLW(L,NY,NX)+VOLI(L,NY,NX)
!    2+DVOLW(L,NY,NX)+DVOLI(L,NY,NX)
!     VOLX(L,NY,NX)=VOLA(L,NY,NX)
!     VOLT(L,NY,NX)=VOLA(L,NY,NX)
    ENDIF
    ENGY=VHCPX*TKSX
    VHCP(L,NY,NX)=VHCM(L,NY,NX) &
      +cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX)) &
      +cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
    TVHCP=TVHCP+VHCP(L,NY,NX)
    TVHCM=TVHCM+VHCM(L,NY,NX)
    TVOLW=TVOLW+VOLW(L,NY,NX)
    TVOLWH=TVOLWH+VOLWH(L,NY,NX)
    TVOLI=TVOLI+VOLI(L,NY,NX)
    TVOLIH=TVOLIH+VOLIH(L,NY,NX)
    TENGY=TENGY+ENGY
    !
    !     ARTIFICIAL SOIL WARMING
    !
    !     IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NU(NY,NX)
    !    3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
    !     THFLW(L,NY,NX)=THFLW(L,NY,NX)
    !    2+(TKSZ(I,J,L)-TKS(L,NY,NX))*VHCP(L,NY,NX)
    !     WRITE(*,3379)'TKSZ',I,J,NX,NY,L,TKSZ(I,J,L)
    !    2,TKS(L,NY,NX),VHCP(L,NY,NX),THFLW(L,NY,NX)
    !3379  FORMAT(A8,6I4,12E12.4)
    !     ENDIF
    !
    !     END ARTIFICIAL SOIL WARMING
    !
    IF(VHCP(L,NY,NX).GT.ZEROS(NY,NX))THEN
      TKS00=(ENGY+THFLW(L,NY,NX)+THTHAW(L,NY,NX) &
        +TUPHT(L,NY,NX)+HWFLU(L,NY,NX))/VHCP(L,NY,NX)

      if(abs(TKS00-TKSX)>100._r8)then
        TKS(L,NY,NX)=TKS(NUM(NY,NX),NY,NX)
        write(*,*)'line',__LINE__,L,TKS(L,NY,NX),TKS(NUM(NY,NX),NY,NX)
        write(*,*)'ENGY+THFLW(L,NY,NX)+THTHAW(L,NY,NX)=',ENGY,THFLW(L,NY,NX),THTHAW(L,NY,NX)
        write(*,*)'TUPHT(L,NY,NX)+HWFLU(L,NY,NX)=',TUPHT(L,NY,NX),HWFLU(L,NY,NX)
        write(*,*)'VHCP(L,NY,NX)=',VHCP(L,NY,NX)
        write(*,*)'NUM',NUM(NY,NX),'L=',L
        write(*,*)'water',VOLW(L,NY,NX),VOLWH(L,NY,NX)
        write(*,*)'ice',VOLI(L,NY,NX),VOLIH(L,NY,NX)
!        call endrun(trim(mod_filename)//' at line',__LINE__)
      else
        TKS(L,NY,NX)=TKS00
      endif
    ELSE
      TKS(L,NY,NX)=TKS(NUM(NY,NX),NY,NX)
    ENDIF
    TCS(L,NY,NX)=TKS(L,NY,NX)-TC2K
    TSMX(L,NY,NX)=AMAX1(TSMX(L,NY,NX),TCS(L,NY,NX))
    TSMN(L,NY,NX)=AMIN1(TSMN(L,NY,NX),TCS(L,NY,NX))
    UN2GS(NY,NX)=UN2GS(NY,NX)+XN2GS(L,NY,NX)
    !     IF(NX.EQ.1)THEN
    !     WRITE(*,6547)'VOLW',I,J,NX,NY,L,VOLW(L,NY,NX),VOLI(L,NY,NX)
    !    2,VOLW(L,NY,NX)+VOLI(L,NY,NX),VOLP(L,NY,NX),VOLA(L,NY,NX)
    !    3,VOLWH(L,NY,NX),VOLIH(L,NY,NX),VOLAH(L,NY,NX)
    !    4,DVOLI(L,NY,NX),TTHAW(L,NY,NX),TTHAWH(L,NY,NX)
    !    5,TFLW(L,NY,NX),FLW(3,L,NY,NX),FLW(3,L+1,NY,NX)
    !    2,FLW(2,L,NY+1,NX),FLW(1,L,NY,NX+1)
    !    3,FINH(L,NY,NX),TUPWTR(L,NY,NX),FLU(L,NY,NX),TQR(NY,NX)
    !    4,PSISM(L,NY,NX),PSISM(L+1,NY,NX),VOLX(L,NY,NX),VOLT(L,NY,NX)
    !     WRITE(*,6547)'VOLWH',I,J,NX,NY,L,VOLWH(L,NY,NX)
    !    5,FLWH(3,L,NY,NX),FLWH(3,L+1,NY,NX)
    !    2,TFLWH(L,NY,NX),FINH(L,NY,NX),TTHAWH(L,NY,NX)
    !    3,VOLIH(L,NY,NX),VOLAH(L,NY,NX)
    !6547  FORMAT(A8,5I4,40E16.8)
    !     WRITE(*,6633)'TKS',I,J,NX,NY,L,TKS(L,NY,NX),ENGY,THFLW(L,NY,NX)
    !    2,THTHAW(L,NY,NX),TUPHT(L,NY,NX),HWFLU(L,NY,NX),VHCP(L,NY,NX)
    !    3,TVHCP,TVHCM,TVOLW,TVOLWH,TVOLI,TVOLIH,TENGY,TKSX,VHCPX
    !    3,VOLWXX,VOLIXX,VOLW(L,NY,NX),VOLWH(L,NY,NX),VOLI(L,NY,NX)
    !    4,VOLIH(L,NY,NX),TFLW(L,NY,NX),FINH(L,NY,NX),TTHAW(L,NY,NX)
    !    5,TUPWTR(L,NY,NX),FLU(L,NY,NX),TQR(NY,NX)
    !    6,HFLW(3,L,NY,NX),HFLW(3,L+1,NY,NX)
    !    7,ENGY+THFLW(L,NY,NX)+THTHAW(L,NY,NX)+TUPHT(L,NY,NX)+HWFLU(L,NY,NX)
    !6633  FORMAT(A8,5I4,40E16.8)
    !     ENDIF
    !
    !     RESIDUE FROM PLANT LITTERFALL
!
    DO 8565 K=0,1
      DO  M=1,4
        OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)+CSNT(M,K,L,NY,NX)
        DO NGL=1,JG
          OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)+CSNT(M,K,L,NY,NX)*micpar%OMCI(1+(NGL-1)*3,K)
        ENDDO
        OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+ZSNT(M,K,L,NY,NX)
        OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+PSNT(M,K,L,NY,NX)
!     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
!     WRITE(*,8484)'OSC',I,J,L,K,M,OSC(M,K,L,NY,NX)
!    2,OSN(M,K,L,NY,NX),OSP(M,K,L,NY,NX),CSNT(M,K,L,NY,NX)
!    3,ZSNT(M,K,L,NY,NX),PSNT(M,K,L,NY,NX)
!8484  FORMAT(A8,5I4,12E12.4)
!     ENDIF
      enddo
8565  CONTINUE
!
    !     DOC, DON, DOP FROM AQUEOUS TRANSPORT
!
    DO 8560 K=0,jcplx1
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+TOCFLS(K,L,NY,NX)+XOCFXS(K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+TONFLS(K,L,NY,NX)+XONFXS(K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+TOPFLS(K,L,NY,NX)+XOPFXS(K,L,NY,NX)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+TOAFLS(K,L,NY,NX)+XOAFXS(K,L,NY,NX)
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)+TOCFHS(K,L,NY,NX)-XOCFXS(K,L,NY,NX)
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)+TONFHS(K,L,NY,NX)-XONFXS(K,L,NY,NX)
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)+TOPFHS(K,L,NY,NX)-XOPFXS(K,L,NY,NX)
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)+TOAFHS(K,L,NY,NX)-XOAFXS(K,L,NY,NX)
      !     IF(L.EQ.1)THEN
      !      WRITE(*,2627)'OQCL',I,J,NX,NY,L,K
      !    2,OQC(K,L,NY,NX),TOCFLS(K,L,NY,NX),XOCFXS(K,L,NY,NX)
      !    4,OQCH(K,L,NY,NX),TOCFHS(K,L,NY,NX),XOCFXS(K,L,NY,NX)
      !    3,OQN(K,L,NY,NX),TONFLS(K,L,NY,NX),XONFXS(K,L,NY,NX)
      !    4,OQNH(K,L,NY,NX),TONFHS(K,L,NY,NX),XONFXS(K,L,NY,NX)
      !    5,OQA(K,L,NY,NX),TOAFLS(K,L,NY,NX),XOAFXS(K,L,NY,NX)
      !    5,OQAH(K,L,NY,NX),TOAFHS(K,L,NY,NX),XOAFXS(K,L,NY,NX)
      !2627  FORMAT(A8,6I4,20E12.4)
      !     ENDIF
8560  CONTINUE
    !
    !     DOC, DON, DOP FROM PLANT EXUDATION
    !
    DO 195 K=0,jcplx1
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+TDFOMC(K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+TDFOMN(K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+TDFOMP(K,L,NY,NX)
195   CONTINUE
    !
    !     SOIL SOLUTES FROM AQUEOUS TRANSPORT, MICROBIAL AND ROOT
    !     EXCHANGE, EQUILIBRIUM REACTIONS, GAS EXCHANGE,
    !     MICROPORE-MACROPORE EXCHANGE,
    !
    CO2S(L,NY,NX)=CO2S(L,NY,NX)+TCOFLS(L,NY,NX)+XCODFG(L,NY,NX) &
      -RCO2O(L,NY,NX)-TCO2S(L,NY,NX)+RCOFLU(L,NY,NX)+XCOFXS(L,NY,NX) &
      +TRCO2(L,NY,NX)+XCOBBL(L,NY,NX)
    CH4S(L,NY,NX)=CH4S(L,NY,NX)+TCHFLS(L,NY,NX)+XCHDFG(L,NY,NX) &
      -RCH4O(L,NY,NX)-TUPCHS(L,NY,NX)+RCHFLU(L,NY,NX) &
      +XCHFXS(L,NY,NX)+XCHBBL(L,NY,NX)
    OXYS(L,NY,NX)=OXYS(L,NY,NX)+TOXFLS(L,NY,NX)+XOXDFG(L,NY,NX) &
      -RUPOXO(L,NY,NX)-TUPOXS(L,NY,NX)+ROXFLU(L,NY,NX) &
      +XOXFXS(L,NY,NX)+XOXBBL(L,NY,NX)
    !     IF(I.EQ.287)THEN
    !     WRITE(*,5432)'CO2SL',I,J,NX,NY,L,CO2S(L,NY,NX),TCOFLS(L,NY,NX)
    !    2,XCODFG(L,NY,NX),RCO2O(L,NY,NX),TCO2S(L,NY,NX)
    !    3,RCOFLU(L,NY,NX),XCOFXS(L,NY,NX),TRCO2(L,NY,NX)
    !    4,XCOBBL(L,NY,NX),CO2G(L,NY,NX),ORGC(L,NY,NX)
    !     WRITE(*,5432)'CH4SL',I,J,NX,NY,L,CH4S(L,NY,NX),CH4G(L,NY,NX)
    !    2,TCHFLS(L,NY,NX),XCHDFG(L,NY,NX),RCH4O(L,NY,NX),TUPCHS(L,NY,NX)
    !    3,RCHFLU(L,NY,NX),XCHFXS(L,NY,NX),XCHBBL(L,NY,NX)
    !    4,XCOBBL(L,NY,NX),XCHFLS(3,L,NY,NX),XCHFLS(3,L+1,NY,NX)
    !     WRITE(*,5432)'OXYSL',I,J,NX,NY,L,OXYS(L,NY,NX),OXYG(L,NY,NX)
    !    4,XOXFLS(3,L,NY,NX),XOXFLS(3,L+1,NY,NX),TOXFLS(L,NY,NX)
    !    2,XOXDFG(L,NY,NX),RUPOXO(L,NY,NX),TUPOXS(L,NY,NX)
    !    3,ROXFLU(L,NY,NX),XOXFXS(L,NY,NX),XOXBBL(L,NY,NX)
    !    5,XOXDFS(NY,NX),COXYS(L,NY,NX)
    !5432  FORMAT(A8,5I4,20E12.4)
    !     ENDIF
    Z2GS(L,NY,NX)=Z2GS(L,NY,NX)+TNGFLS(L,NY,NX)+XNGDFG(L,NY,NX) &
      -RN2G(L,NY,NX)-TUPNF(L,NY,NX)+RNGFLU(L,NY,NX)+XNGFXS(L,NY,NX) &
      -XN2GS(L,NY,NX)+XNGBBL(L,NY,NX)
    Z2OS(L,NY,NX)=Z2OS(L,NY,NX)+TN2FLS(L,NY,NX)+XN2DFG(L,NY,NX) &
      -RN2O(L,NY,NX)-TUPN2S(L,NY,NX)+RN2FLU(L,NY,NX)+XN2FXS(L,NY,NX) &
      +XN2BBL(L,NY,NX)
    !     IF(NY.EQ.5.AND.L.EQ.1)THEN
    !     WRITE(*,4444)'Z2OS',I,J,NX,NY,L,Z2OS(L,NY,NX),TN2FLS(L,NY,NX)
    !    2,XN2DFG(L,NY,NX),RN2O(L,NY,NX),TUPN2S(L,NY,NX),RN2FLU(L,NY,NX)
    !    3,XN2FXS(L,NY,NX),XN2BBL(L,NY,NX),XN2FLS(3,L,NY,NX)
    !    4,XN2FLS(3,L+1,NY,NX),XN2DFS(NY,NX)
    !    3,Z2GS(L,NY,NX),TNGFLS(L,NY,NX),XNGDFG(L,NY,NX)
    !    4,RN2G(L,NY,NX),TUPNF(L,NY,NX),RNGFLU(L,NY,NX),XNGFXS(L,NY,NX)
    !    5,XN2GS(L,NY,NX),XNGBBL(L,NY,NX),XNGDFS(NY,NX)
    !     ENDIF
    H2GS(L,NY,NX)=H2GS(L,NY,NX)+THGFLS(L,NY,NX)+XHGDFG(L,NY,NX) &
      -RH2GO(L,NY,NX)-TUPHGS(L,NY,NX)+RHGFLU(L,NY,NX) &
      +XHGFXS(L,NY,NX)+XHGBBL(L,NY,NX)
    ZNH3S(L,NY,NX)=ZNH3S(L,NY,NX)+TN3FLS(L,NY,NX)+XN3DFG(L,NY,NX) &
      +TRN3S(L,NY,NX)-TUPN3S(L,NY,NX)+RN3FLU(L,NY,NX) &
      +XN3FXW(L,NY,NX)+XN3BBL(L,NY,NX)
    ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+TN4FLS(L,NY,NX)+XNH4S(L,NY,NX) &
      +TRN4S(L,NY,NX)-TUPNH4(L,NY,NX)+RN4FLU(L,NY,NX) &
      +XN4FXW(L,NY,NX)
    !     IF(IYRC.EQ.2006.AND.I.EQ.361.AND.NX.EQ.1)THEN
    !     WRITE(*,4443)'H2GS',I,J,NX,NY,L,H2GS(L,NY,NX),THGFLS(L,NY,NX)
    !    2,XHGDFG(L,NY,NX),RH2GO(L,NY,NX),TUPHGS(L,NY,NX),RHGFLU(L,NY,NX)
    !    3,XHGFXS(L,NY,NX),XHGBBL(L,NY,NX),XHGDFS(NY,NX)
    !     WRITE(*,4444)'NH3',I,J,NX,NY,L,ZNH3S(L,NY,NX),TN3FLS(L,NY,NX)
    !    2,XN3DFG(L,NY,NX),TRN3S(L,NY,NX),TUPN3S(L,NY,NX)
    !    3,RN3FLU(L,NY,NX),XN3FXW(L,NY,NX),XN3BBL(L,NY,NX),XN3DFS(NY,NX)
    !    4,ZNH4S(L,NY,NX)
    !    4,TN4FLS(L,NY,NX),XNH4S(L,NY,NX),TRN4S(L,NY,NX),TUPNH4(L,NY,NX)
    !    5,RN4FLU(L,NY,NX),XN4FXW(L,NY,NX),TN4QRS(NY,NX),TN3QRS(NY,NX)
    !    6,ZNH3SH(L,NY,NX),ZNH4SH(L,NY,NX),14.0*XN4(L,NY,NX)
    !4443  FORMAT(A8,5I4,30F16.8)
    !4444  FORMAT(A8,5I4,30E12.4)
    !     ENDIF
    ZNO3S(L,NY,NX)=ZNO3S(L,NY,NX)+TNOFLS(L,NY,NX)+XNO3S(L,NY,NX) &
      +TRNO3(L,NY,NX)-TUPNO3(L,NY,NX)+RNOFLU(L,NY,NX) &
      +XNOFXW(L,NY,NX)
    ZNO2S(L,NY,NX)=ZNO2S(L,NY,NX)+TNXFLS(L,NY,NX)+XNO2S(L,NY,NX) &
      +TRNO2(L,NY,NX)+XNXFXS(L,NY,NX)
    !     IF(L.EQ.NU(NY,NX))THEN
    !     WRITE(*,5545)'NO3',I,J,NX,NY,L,ZNO3S(L,NY,NX),TNOFLS(L,NY,NX)
    !    2,XNO3S(L,NY,NX),TRNO3(L,NY,NX),TUPNO3(L,NY,NX),RNOFLU(L,NY,NX)
    !    3,XNOFXW(L,NY,NX),ZNO2S(L,NY,NX),TNXFLS(L,NY,NX)
    !    4,XNO2S(L,NY,NX),TRNO2(L,NY,NX),XNXFXS(L,NY,NX),TNXQRS(NY,NX)
    !5545  FORMAT(A8,5I4,40E12.4)
    !     ENDIF
    H1PO4(L,NY,NX)=H1PO4(L,NY,NX)+TP1FLS(L,NY,NX)+XH1PS(L,NY,NX) &
      +TRH1P(L,NY,NX)-TUPH1P(L,NY,NX)+RH1PFU(L,NY,NX)+XH1PXS(L,NY,NX)
    H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+TPOFLS(L,NY,NX)+XH2PS(L,NY,NX) &
      +TRH2P(L,NY,NX)-TUPH2P(L,NY,NX)+RH2PFU(L,NY,NX)+XH2PXS(L,NY,NX)
    ZNH3B(L,NY,NX)=ZNH3B(L,NY,NX)+TN3FLB(L,NY,NX)+XNBDFG(L,NY,NX) &
      +TRN3B(L,NY,NX)-TUPN3B(L,NY,NX)+RN3FBU(L,NY,NX) &
      +XN3FXB(L,NY,NX)+XNBBBL(L,NY,NX)
    ZNH4B(L,NY,NX)=ZNH4B(L,NY,NX)+TN4FLB(L,NY,NX)+XNH4B(L,NY,NX) &
      +TRN4B(L,NY,NX)-TUPNHB(L,NY,NX)+RN4FBU(L,NY,NX) &
      +XN4FXB(L,NY,NX)
    ZNO3B(L,NY,NX)=ZNO3B(L,NY,NX)+TNOFLB(L,NY,NX)+XNO3B(L,NY,NX) &
      +TRNOB(L,NY,NX)-TUPNOB(L,NY,NX)+RNOFBU(L,NY,NX) &
      +XNOFXB(L,NY,NX)
    ZNO2B(L,NY,NX)=ZNO2B(L,NY,NX)+TNXFLB(L,NY,NX)+XNO2B(L,NY,NX) &
      +TRN2B(L,NY,NX)+XNXFXB(L,NY,NX)
    H1POB(L,NY,NX)=H1POB(L,NY,NX)+TH1BFB(L,NY,NX)+XH1BS(L,NY,NX) &
      +TRH1B(L,NY,NX)-TUPH1B(L,NY,NX)+RH1BBU(L,NY,NX) &
      +XH1BXB(L,NY,NX)
    H2POB(L,NY,NX)=H2POB(L,NY,NX)+TH2BFB(L,NY,NX)+XH2BS(L,NY,NX) &
      +TRH2B(L,NY,NX) -TUPH2B(L,NY,NX)+RH2BBU(L,NY,NX) &
      +XH2BXB(L,NY,NX)
    THRE(NY,NX)=THRE(NY,NX)+RCO2O(L,NY,NX)+RCH4O(L,NY,NX)
!    print*,'curday',RCO2O(L,NY,NX),RCH4O(L,NY,NX)
    UN2GG(NY,NX)=UN2GG(NY,NX)+RN2G(L,NY,NX)
    HN2GG(NY,NX)=HN2GG(NY,NX)+RN2G(L,NY,NX)
    !
    !     EXCHANGEABLE CATIONS AND ANIONS FROM EXCHANGE REACTIONS
    !
    XN4(L,NY,NX)=XN4(L,NY,NX)+TRXN4(L,NY,NX)
    XNB(L,NY,NX)=XNB(L,NY,NX)+TRXNB(L,NY,NX)
    XOH0(L,NY,NX)=XOH0(L,NY,NX)+TRXH0(L,NY,NX)
    XOH1(L,NY,NX)=XOH1(L,NY,NX)+TRXH1(L,NY,NX)
    XOH2(L,NY,NX)=XOH2(L,NY,NX)+TRXH2(L,NY,NX)
    XH1P(L,NY,NX)=XH1P(L,NY,NX)+TRX1P(L,NY,NX)
    XH2P(L,NY,NX)=XH2P(L,NY,NX)+TRX2P(L,NY,NX)
    XOH0B(L,NY,NX)=XOH0B(L,NY,NX)+TRBH0(L,NY,NX)
    XOH1B(L,NY,NX)=XOH1B(L,NY,NX)+TRBH1(L,NY,NX)
    XOH2B(L,NY,NX)=XOH2B(L,NY,NX)+TRBH2(L,NY,NX)
    XH1PB(L,NY,NX)=XH1PB(L,NY,NX)+TRB1P(L,NY,NX)
    XH2PB(L,NY,NX)=XH2PB(L,NY,NX)+TRB2P(L,NY,NX)
    !     IF(J.EQ.12.AND.L.LE.4)THEN
    !     WRITE(*,4445)'NHB',I,J,NX,NY,L,ZNH3B(L,NY,NX),TN3FLB(L,NY,NX)
    !    2,XNBDFG(L,NY,NX),TRN3B(L,NY,NX),TUPN3B(L,NY,NX)
    !    3,RN3FBU(L,NY,NX),XN3FXB(L,NY,NX),XNBBBL(L,NY,NX),TUPNHB(L,NY,NX)
    !    4,ZNH4B(L,NY,NX),TN4FLB(L,NY,NX),XNH4B(L,NY,NX)
    !    5,TRN4B(L,NY,NX),TUPNHB(L,NY,NX),RN4FBU(L,NY,NX),XNB(L,NY,NX)*14.0
    !     WRITE(*,4445)'NOB',I,J,NX,NY,L,ZNO2B(L,NY,NX),TNXFLB(L,NY,NX)
    !    2,XNO2B(L,NY,NX),TRN2B(L,NY,NX),XNXFXB(L,NY,NX)
    !4445  FORMAT(A8,5I4,20E12.4)
    !     ENDIF
    !
    !     PRECIPITATES FROM PRECIPITATION-DISSOLUTION REACTIONS
    !
    PALPO(L,NY,NX)=PALPO(L,NY,NX)+TRALPO(L,NY,NX)
    PFEPO(L,NY,NX)=PFEPO(L,NY,NX)+TRFEPO(L,NY,NX)
    PCAPD(L,NY,NX)=PCAPD(L,NY,NX)+TRCAPD(L,NY,NX)
    PCAPH(L,NY,NX)=PCAPH(L,NY,NX)+TRCAPH(L,NY,NX)
    PCAPM(L,NY,NX)=PCAPM(L,NY,NX)+TRCAPM(L,NY,NX)
    PALPB(L,NY,NX)=PALPB(L,NY,NX)+TRALPB(L,NY,NX)
    PFEPB(L,NY,NX)=PFEPB(L,NY,NX)+TRFEPB(L,NY,NX)
    PCPDB(L,NY,NX)=PCPDB(L,NY,NX)+TRCPDB(L,NY,NX)
    PCPHB(L,NY,NX)=PCPHB(L,NY,NX)+TRCPHB(L,NY,NX)
    PCPMB(L,NY,NX)=PCPMB(L,NY,NX)+TRCPMB(L,NY,NX)
    !
    !     MACROPORE SOLUTES FROM MACROPORE-MICROPORE EXCHANGE
    !
    CO2SH(L,NY,NX)=CO2SH(L,NY,NX)+TCOFHS(L,NY,NX)-XCOFXS(L,NY,NX)
    CH4SH(L,NY,NX)=CH4SH(L,NY,NX)+TCHFHS(L,NY,NX)-XCHFXS(L,NY,NX)
    OXYSH(L,NY,NX)=OXYSH(L,NY,NX)+TOXFHS(L,NY,NX)-XOXFXS(L,NY,NX)
    Z2GSH(L,NY,NX)=Z2GSH(L,NY,NX)+TNGFHS(L,NY,NX)-XNGFXS(L,NY,NX)
    Z2OSH(L,NY,NX)=Z2OSH(L,NY,NX)+TN2FHS(L,NY,NX)-XN2FXS(L,NY,NX)
    H2GSH(L,NY,NX)=H2GSH(L,NY,NX)+THGFHS(L,NY,NX)-XHGFXS(L,NY,NX)
    ZNH4SH(L,NY,NX)=ZNH4SH(L,NY,NX)+TN4FHS(L,NY,NX)-XN4FXW(L,NY,NX)
    ZNH3SH(L,NY,NX)=ZNH3SH(L,NY,NX)+TN3FHS(L,NY,NX)-XN3FXW(L,NY,NX)
    ZNO3SH(L,NY,NX)=ZNO3SH(L,NY,NX)+TNOFHS(L,NY,NX)-XNOFXW(L,NY,NX)
    ZNO2SH(L,NY,NX)=ZNO2SH(L,NY,NX)+TNXFHS(L,NY,NX)-XNXFXS(L,NY,NX)
    H1PO4H(L,NY,NX)=H1PO4H(L,NY,NX)+TP1FHS(L,NY,NX)-XH1PXS(L,NY,NX)
    H2PO4H(L,NY,NX)=H2PO4H(L,NY,NX)+TPOFHS(L,NY,NX)-XH2PXS(L,NY,NX)
    ZNH4BH(L,NY,NX)=ZNH4BH(L,NY,NX)+TN4FHB(L,NY,NX)-XN4FXB(L,NY,NX)
    ZNH3BH(L,NY,NX)=ZNH3BH(L,NY,NX)+TN3FHB(L,NY,NX)-XN3FXB(L,NY,NX)
    ZNO3BH(L,NY,NX)=ZNO3BH(L,NY,NX)+TNOFHB(L,NY,NX)-XNOFXB(L,NY,NX)
    ZNO2BH(L,NY,NX)=ZNO2BH(L,NY,NX)+TNXFHB(L,NY,NX)-XNXFXB(L,NY,NX)
    H1POBH(L,NY,NX)=H1POBH(L,NY,NX)+TH1BHB(L,NY,NX)-XH1BXB(L,NY,NX)
    H2POBH(L,NY,NX)=H2POBH(L,NY,NX)+TH2BHB(L,NY,NX)-XH2BXB(L,NY,NX)
    !     IF(L.EQ.4)THEN
    !     WRITE(*,4747)'ZNH4SH',I,J,NX,NY,L
    !    2,ZNH4SH(L,NY,NX),TN4FHS(L,NY,NX),XN4FXW(L,NY,NX)
    !    2,ZNH3SH(L,NY,NX),TN3FHS(L,NY,NX),XN3FXW(L,NY,NX)
    !     WRITE(*,4747)'ZNO3SH',I,J,NX,NY,L,ZNO3SH(L,NY,NX)
    !    2,TNOFHS(L,NY,NX),XNOFXW(L,NY,NX)
    !    3,ZNO2SH(L,NY,NX),TNXFHS(L,NY,NX),XNXFXS(L,NY,NX)
    !4747  FORMAT(A8,5I4,12E12.4)
    !     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
    !     WRITE(*,5545)'HP14',I,J,NX,NY,L,H1PO4(L,NY,NX),TP1FLS(L,NY,NX)
    !    2,XH1PS(L,NY,NX),TRH1P(L,NY,NX),TUPH1P(L,NY,NX),RH1PFU(L,NY,NX)
    !    3,XH1PXS(L,NY,NX),XH1P(L,NY,NX),H1POB(L,NY,NX)
    !    4,TH1BFB(L,NY,NX),XH1BS(L,NY,NX),TRH1B(L,NY,NX),TUPH1B(L,NY,NX)
    !    2,RH1BBU(L,NY,NX),XH1BXB(L,NY,NX),XH1PB(L,NY,NX)
    !    2,H1PO4H(L,NY,NX),TP1FHS(L,NY,NX),XH1PXS(L,NY,NX)
    !    2,H1POBH(L,NY,NX),TH1BHB(L,NY,NX),XH1BXB(L,NY,NX)
    !     WRITE(*,5545)'HP24',I,J,NX,NY,L,H2PO4(L,NY,NX),TPOFLS(L,NY,NX)
    !    2,XH2PS(L,NY,NX),TRH2P(L,NY,NX),TUPH2P(L,NY,NX),RH2PFU(L,NY,NX)
    !    3,XH2PXS(L,NY,NX),XH2P(L,NY,NX),H2POB(L,NY,NX)
    !    4,TH2BFB(L,NY,NX),XH2BS(L,NY,NX),TRH2B(L,NY,NX),TUPH2B(L,NY,NX)
    !    5,RH2BBU(L,NY,NX),XH2BXB(L,NY,NX),XH2PB(L,NY,NX)
    !    2,H2PO4H(L,NY,NX),TPOFHS(L,NY,NX),XH2PXS(L,NY,NX)
    !    2,H2POBH(L,NY,NX),TH2BHB(L,NY,NX),XH2BXB(L,NY,NX)
    !     ENDIF
    !     ENDIF
    !
    !     GASES FROM VOLATILIZATION-DISSOLUTION AND GAS TRANSFER
!
    CO2G(L,NY,NX)=CO2G(L,NY,NX)+TCOFLG(L,NY,NX)-XCODFG(L,NY,NX)
    CH4G(L,NY,NX)=CH4G(L,NY,NX)+TCHFLG(L,NY,NX)-XCHDFG(L,NY,NX)
    OXYG(L,NY,NX)=OXYG(L,NY,NX)+TOXFLG(L,NY,NX)-XOXDFG(L,NY,NX)
    Z2GG(L,NY,NX)=Z2GG(L,NY,NX)+TNGFLG(L,NY,NX)-XNGDFG(L,NY,NX)
    Z2OG(L,NY,NX)=Z2OG(L,NY,NX)+TN2FLG(L,NY,NX)-XN2DFG(L,NY,NX)
    ZNH3G(L,NY,NX)=ZNH3G(L,NY,NX)+TNHFLG(L,NY,NX)-XN3DFG(L,NY,NX) &
      -XNBDFG(L,NY,NX)+TRN3G(L,NY,NX)
    H2GG(L,NY,NX)=H2GG(L,NY,NX)+THGFLG(L,NY,NX)-XHGDFG(L,NY,NX)
    ROXYF(L,NY,NX)=TOXFLG(L,NY,NX)
    RCO2F(L,NY,NX)=TCOFLG(L,NY,NX)
    RCH4F(L,NY,NX)=TCHFLG(L,NY,NX)
    ROXYL(L,NY,NX)=TOXFLS(L,NY,NX)+ROXFLU(L,NY,NX)+XOXFXS(L,NY,NX) &
      +XOXBBL(L,NY,NX)
    RCH4L(L,NY,NX)=TCHFLS(L,NY,NX)+RCHFLU(L,NY,NX)+XCHFXS(L,NY,NX) &
      +XCHBBL(L,NY,NX)
    !     IF(L.EQ.1)THEN
    !     WRITE(*,5432)'CO2GL',I,J,NX,NY,L,CO2G(L,NY,NX),TCOFLG(L,NY,NX)
    !    2,XCODFG(L,NY,NX),THETP(L,NY,NX)
    !     WRITE(*,5432)'OXYGL',I,J,NX,NY,L,OXYG(L,NY,NX),TOXFLG(L,NY,NX)
    !    2,XOXDFG(L,NY,NX),COXYG(L,NY,NX),XOXFLG(3,L,NY,NX)
    !    3,XOXFLG(3,L+1,NY,NX),XOXFLG(1,L,NY,NX+1)
    !     WRITE(*,5432)'CH4GL',I,J,NX,NY,L,CH4G(L,NY,NX),TCHFLG(L,NY,NX)
    !    2,XCHDFG(L,NY,NX),CCH4G(L,NY,NX),XCHFLG(3,L,NY,NX)
    !    3,XCHFLG(3,L+1,NY,NX),XCHDFS(NY,NX),RCH4F(L,NY,NX)
    !    4,RCH4L(L,NY,NX),TCHFLS(L,NY,NX),RCHFLU(L,NY,NX)
    !    5,XCHFXS(L,NY,NX),XCHBBL(L,NY,NX)
    !     ENDIF
    !
    !     GRID CELL BOUNDARY FLUXES FROM ROOT GAS TRANSFER
!
    HEATIN=HEATIN+THTHAW(L,NY,NX)+TUPHT(L,NY,NX)
    CIB=TCOFLA(L,NY,NX)
    CHB=TCHFLA(L,NY,NX)
    OIB=TOXFLA(L,NY,NX)
!    HGB=THGFLA(L,NY,NX)
    HGB=0.0_r8
    ZGB=0.0_r8
    Z2B=TN2FLA(L,NY,NX)
    ZHB=TNHFLA(L,NY,NX)
!
!     GRID CELL BOUNDARY FLUXES BUBBLING
!
    IF(LG.EQ.0)THEN
      LL=0
      CIB=CIB+XCOBBL(L,NY,NX)
      CHB=CHB+XCHBBL(L,NY,NX)
      OIB=OIB+XOXBBL(L,NY,NX)
      ZGB=ZGB+XNGBBL(L,NY,NX)
      Z2B=Z2B+XN2BBL(L,NY,NX)
      ZHB=ZHB+XN3BBL(L,NY,NX)+XNBBBL(L,NY,NX)
      HGB=HGB+XHGBBL(L,NY,NX)
    ELSE
      LL=MIN(L,LG)
      CO2G(LL,NY,NX)=CO2G(LL,NY,NX)-XCOBBL(L,NY,NX)
      CH4G(LL,NY,NX)=CH4G(LL,NY,NX)-XCHBBL(L,NY,NX)
      OXYG(LL,NY,NX)=OXYG(LL,NY,NX)-XOXBBL(L,NY,NX)
      Z2GG(LL,NY,NX)=Z2GG(LL,NY,NX)-XNGBBL(L,NY,NX)
      Z2OG(LL,NY,NX)=Z2OG(LL,NY,NX)-XN2BBL(L,NY,NX)
      ZNH3G(LL,NY,NX)=ZNH3G(LL,NY,NX)-XN3BBL(L,NY,NX)-XNBBBL(L,NY,NX)
      H2GG(LL,NY,NX)=H2GG(LL,NY,NX)-XHGBBL(L,NY,NX)
      IF(LG.LT.L)THEN
        TLCO2G=TLCO2G-XCOBBL(L,NY,NX)-XCHBBL(L,NY,NX)
        UCO2S(NY,NX)=UCO2S(NY,NX)-XCOBBL(L,NY,NX)-XCHBBL(L,NY,NX)
        OXYGSO=OXYGSO-XOXBBL(L,NY,NX)
        TLN2G=TLN2G-XNGBBL(L,NY,NX)-XN2BBL(L,NY,NX) &
          -XN3BBL(L,NY,NX)-XNBBBL(L,NY,NX)
        TLH2G=TLH2G-XHGBBL(L,NY,NX)
      ENDIF
    ENDIF
    CO2GIN=CO2GIN+CIB+CHB
    COB=TCO2P(L,NY,NX)+TCO2S(L,NY,NX)-TRCO2(L,NY,NX)
    TCOU=TCOU+COB
    HCO2G(NY,NX)=HCO2G(NY,NX)+CIB
    UCO2G(NY,NX)=UCO2G(NY,NX)+CIB
    HCH4G(NY,NX)=HCH4G(NY,NX)+CHB
    UCH4G(NY,NX)=UCH4G(NY,NX)+CHB
    UCOP(NY,NX)=UCOP(NY,NX)+TCO2P(L,NY,NX)+TCO2S(L,NY,NX)
    UDICD(NY,NX)=UDICD(NY,NX)-12.0*TBCO2(L,NY,NX)
    TXCO2(NY,NX)=TXCO2(NY,NX)+12.0*TBCO2(L,NY,NX)
    OXYGIN=OXYGIN+OIB
    OOB=RUPOXO(L,NY,NX)+TUPOXP(L,NY,NX)+TUPOXS(L,NY,NX)
    OXYGOU=OXYGOU+OOB
    UOXYG(NY,NX)=UOXYG(NY,NX)+OIB
    HOXYG(NY,NX)=HOXYG(NY,NX)+OIB
    H2GIN=H2GIN+HGB
    HOB=RH2GO(L,NY,NX)+TUPHGS(L,NY,NX)
    H2GOU=H2GOU+HOB
    ZN2GIN=ZN2GIN+ZGB+Z2B+ZHB
!     UN2GG(NY,NX)=UN2GG(NY,NX)+ZGB
!     HN2GG(NY,NX)=HN2GG(NY,NX)+ZGB
    UN2OG(NY,NX)=UN2OG(NY,NX)+Z2B
    HN2OG(NY,NX)=HN2OG(NY,NX)+Z2B
    UNH3G(NY,NX)=UNH3G(NY,NX)+ZHB
    HNH3G(NY,NX)=HNH3G(NY,NX)+ZHB
    UH2GG(NY,NX)=UH2GG(NY,NX)+HGB
    !     IF(NY.EQ.5)THEN
    !     WRITE(*,6645)'PLT',I,J,NX,NY,L,LG,LL
    !    2,TXCO2NY,NX),TBCO2(L,NY,NX)
    !    2,HCH4G(NY,NX),CHB,TCHFLA(L,NY,NX),XCHBBL(L,NY,NX)
    !    2,HOXYG(NY,NX),OIB
    !    3,XOXBBL(L,NY,NX),TUPOXP(L,NY,NX),TUPOXS(L,NY,NX)
    !    4,TOXFLA(L,NY,NX),OXYG(L,NY,NX),SOXYL(L,NY,NX)
    !    4,HCO2G(NY,NX),CIB,TCOFLA(L,NY,NX),XCOBBL(L,NY,NX)
    !    4,TRCO2(L,NY,NX)
    !    2,UN2OG(NY,NX),ZGI,XN2BBL(L,NY,NX)
    !    5,TN2FLA(L,NY,NX),TNHFLA(L,NY,NX),THGFLA(L,NY,NX)
    !    2,UN2GG(NY,NX),ZGI,XNGBBL(L,NY,NX)
    !    5,TN2FLA(L,NY,NX),TNHFLA(L,NY,NX),THGFLA(L,NY,NX)
    !    6,CH4G(LL,NY,NX)
    !6645  FORMAT(A8,7I4,30E12.4)
    !      ENDIF
    !
    !     GRID CELL BOUNDARY FLUXES FROM EQUILIBRIUM REACTIONS
!
    SNM=(2.0*(XNH4S(L,NY,NX)+XNH4B(L,NY,NX)-TUPNH4(L,NY,NX) &
      -TUPNHB(L,NY,NX)-XN2GS(L,NY,NX))-TUPN3S(L,NY,NX)-TUPN3B(L,NY,NX) &
      +XNO3S(L,NY,NX)+XNO3B(L,NY,NX)-TUPNO3(L,NY,NX)-TUPNOB(L,NY,NX) &
      +XNO2S(L,NY,NX)+XNO2B(L,NY,NX))/14.0
    SPM=(2.0*(XH1PS(L,NY,NX)+XH1BS(L,NY,NX)-TUPH1P(L,NY,NX) &
      -TUPH1B(L,NY,NX))+3.0*(XH2PS(L,NY,NX)+XH2BS(L,NY,NX) &
      -TUPH2P(L,NY,NX)-TUPH2B(L,NY,NX)))/31.0
    SSB=TRH2O(L,NY,NX)+TBCO2(L,NY,NX)+XZHYS(L,NY,NX) &
      +TBION(L,NY,NX)
      TIONOU=TIONOU-SSB
    !     UIONOU(NY,NX)=UIONOU(NY,NX)-SSB
    !     WRITE(20,3339)'SSB',I,J,L,SSB,TRH2O(L,NY,NX)
    !    2,TBCO2(L,NY,NX),XZHYS(L,NY,NX),TBION(L,NY,NX)
    !
    !     GAS AND SOLUTE EXCHANGE WITHIN GRID CELL ADDED TO ECOSYSTEM
    !     TOTALS FOR CALCULATING COMPETITION CONSTRAINTS ON MICROBIAL
    !     AND ROOT POPULATIONS
!
    DO 7990 K=0,jcplx1
      DO 7980 N=1,7
        DO NGL=1,JG
          ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+ROXYS(NGL,N,K,L,NY,NX)
          RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+RVMX4(NGL,N,K,L,NY,NX)+RINHO(NGL,N,K,L,NY,NX)
          RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+RVMX3(NGL,N,K,L,NY,NX)+RINOO(NGL,N,K,L,NY,NX)
          RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+RVMX2(NGL,N,K,L,NY,NX)
          RN2OX(L,NY,NX)=RN2OX(L,NY,NX)+RVMX1(NGL,N,K,L,NY,NX)
          RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+RIPOO(NGL,N,K,L,NY,NX)
          RP14X(L,NY,NX)=RP14X(L,NY,NX)+RIPO1(NGL,N,K,L,NY,NX)
          RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+RVMB4(NGL,N,K,L,NY,NX)+RINHB(NGL,N,K,L,NY,NX)
          RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+RVMB3(NGL,N,K,L,NY,NX)+RINOB(NGL,N,K,L,NY,NX)
          RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+RVMB2(NGL,N,K,L,NY,NX)
          RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+RIPBO(NGL,N,K,L,NY,NX)
          RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+RIPB1(NGL,N,K,L,NY,NX)
          ROQCX(K,L,NY,NX)=ROQCX(K,L,NY,NX)+ROQCS(NGL,N,K,L,NY,NX)
          ROQAX(K,L,NY,NX)=ROQAX(K,L,NY,NX)+ROQAS(NGL,N,K,L,NY,NX)
        enddo
7980  CONTINUE
7990  CONTINUE
      DO  N=1,7
        DO NGL=1,JG
          ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+ROXYSff(NGL,N,L,NY,NX)
          RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+RVMX4ff(NGL,N,L,NY,NX)+RINHOff(NGL,N,L,NY,NX)
          RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+RVMX3ff(NGL,N,L,NY,NX)+RINOOff(NGL,N,L,NY,NX)
          RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+RVMX2ff(NGL,N,L,NY,NX)
          RN2OX(L,NY,NX)=RN2OX(L,NY,NX)+RVMX1ff(NGL,N,L,NY,NX)
          RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+RIPOOff(NGL,N,L,NY,NX)
          RP14X(L,NY,NX)=RP14X(L,NY,NX)+RIPO1ff(NGL,N,L,NY,NX)
          RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+RVMB4ff(NGL,N,L,NY,NX)+RINHBff(NGL,N,L,NY,NX)
          RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+RVMB3ff(NGL,N,L,NY,NX)+RINOBff(NGL,N,L,NY,NX)
          RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+RVMB2ff(NGL,N,L,NY,NX)
          RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+RIPBOff(NGL,N,L,NY,NX)
          RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+RIPB1ff(NGL,N,L,NY,NX)
        enddo
    ENDDO

    RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+RVMXC(L,NY,NX)
    RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+RVMBC(L,NY,NX)
    !
    !     GRID CELL VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
    !     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
    !     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS SINCE START OF RUN
    !
    !     IF(J.EQ.24)THEN
    WS=VOLW(L,NY,NX)+VOLWH(L,NY,NX) &
      +(VOLI(L,NY,NX)+VOLIH(L,NY,NX))*DENSI
    VOLWSO=VOLWSO+WS
    VOLISO=VOLISO+VOLI(L,NY,NX)+VOLIH(L,NY,NX)
    UVOLW(NY,NX)=UVOLW(NY,NX)+WS
!    2-WP(L,NY,NX)*VOLX(L,NY,NX)
    HEATSO=HEATSO+VHCP(L,NY,NX)*TKS(L,NY,NX)
    SD=SAND(L,NY,NX)+SILT(L,NY,NX)+CLAY(L,NY,NX)
    TSEDSO=TSEDSO+SD
    CS=CO2G(L,NY,NX)+CO2S(L,NY,NX)+CO2SH(L,NY,NX)+TLCO2P(L,NY,NX) &
      +CH4G(L,NY,NX)+CH4S(L,NY,NX)+CH4SH(L,NY,NX)+TLCH4P(L,NY,NX)
    TLCO2G=TLCO2G+CS
    UCO2S(NY,NX)=UCO2S(NY,NX)+CS
    HS=H2GG(L,NY,NX)+H2GS(L,NY,NX)+H2GSH(L,NY,NX)+TLH2GP(L,NY,NX)
    TLH2G=TLH2G+HS
    !     IF(NX.EQ.1.AND.NY.EQ.1)THEN
    !     WRITE(*,8642)'TLCO2G',I,J,L,TLCO2G,CS,CO2G(L,NY,NX),CO2S(L,NY,NX)
    !    2,CO2SH(L,NY,NX),TLCO2P(L,NY,NX),CH4G(L,NY,NX),CH4S(L,NY,NX)
    !    3,CH4SH(L,NY,NX),TLCH4P(L,NY,NX)
    !8642  FORMAT(A8,3I4,20F16.6)
    !     ENDIF
    OS=OXYG(L,NY,NX)+OXYS(L,NY,NX)+OXYSH(L,NY,NX)+TLOXYP(L,NY,NX)
    OXYGSO=OXYGSO+OS
    ZG=Z2GG(L,NY,NX)+Z2GS(L,NY,NX)+Z2GSH(L,NY,NX)+TLN2OP(L,NY,NX) &
      +Z2OG(L,NY,NX)+Z2OS(L,NY,NX)+Z2OSH(L,NY,NX)+TLNH3P(L,NY,NX) &
      +ZNH3G(L,NY,NX)
    TLN2G=TLN2G+ZG
    Z4S=ZNH4S(L,NY,NX)+ZNH4SH(L,NY,NX)+ZNH4B(L,NY,NX) &
      +ZNH4BH(L,NY,NX)+ZNH3S(L,NY,NX)+ZNH3SH(L,NY,NX) &
      +ZNH3B(L,NY,NX)+ZNH3BH(L,NY,NX)
    Z4X=14.0*(XN4(L,NY,NX)+XNB(L,NY,NX))
    Z4F=14.0*(ZNH4FA(L,NY,NX)+ZNHUFA(L,NY,NX)+ZNH3FA(L,NY,NX) &
      +ZNH4FB(L,NY,NX)+ZNHUFB(L,NY,NX)+ZNH3FB(L,NY,NX))
    TLNH4=TLNH4+Z4S+Z4X+Z4F
    UNH4(NY,NX)=UNH4(NY,NX)+Z4S+Z4X
    !     IF(I.EQ.168)THEN
    !     WRITE(*,5455)'TLN2GL',I,J,NX,NY,L,TLN2G
    !    2,ZG,Z2GG(L,NY,NX),Z2GS(L,NY,NX),Z2GSH(L,NY,NX),TLN2OP(L,NY,NX)
    !    2,Z2OG(L,NY,NX),Z2OS(L,NY,NX),Z2OSH(L,NY,NX),TLNH3P(L,NY,NX)
    !    3,ZNH3G(L,NY,NX)
    !     WRITE(*,5455)'TLNH4L',I,J,NX,NY,L,TLNH4,UNH4(NY,NX)
    !    2,Z4S,Z4X,Z4F,XN4(L,NY,NX)
    !    2,XNB(L,NY,NX),ZNH4S(L,NY,NX),ZNH4SH(L,NY,NX)
    !    3,ZNH4B(L,NY,NX),ZNH4BH(L,NY,NX),ZNH3S(L,NY,NX),ZNH3SH(L,NY,NX)
    !    4,ZNH3B(L,NY,NX),ZNH3BH(L,NY,NX),TN4FHB(L,NY,NX),XN4FXB(L,NY,NX)
    !5455  FORMAT(A8,5I4,30E12.4)
    !     ENDIF
    ZOS=ZNO3S(L,NY,NX)+ZNO3SH(L,NY,NX)+ZNO3B(L,NY,NX) &
      +ZNO3BH(L,NY,NX)+ZNO2S(L,NY,NX)+ZNO2SH(L,NY,NX) &
      +ZNO2B(L,NY,NX)+ZNO2BH(L,NY,NX)
    ZOF=14.0*(ZNO3FA(L,NY,NX)+ZNO3FA(L,NY,NX))
    TLNO3=TLNO3+ZOS+ZOF
    UNO3(NY,NX)=UNO3(NY,NX)+ZOS
    POS=H2PO4(L,NY,NX)+H2PO4H(L,NY,NX)+H2POB(L,NY,NX) &
      +H2POBH(L,NY,NX)+H1PO4(L,NY,NX)+H1PO4H(L,NY,NX) &
      +H1POB(L,NY,NX)+H1POBH(L,NY,NX)
    POX=31.0*(XH1P(L,NY,NX)+XH2P(L,NY,NX) &
      +XH1PB(L,NY,NX)+XH2PB(L,NY,NX))
    POP=31.0*(PALPO(L,NY,NX)+PFEPO(L,NY,NX)+PCAPD(L,NY,NX) &
      +PALPB(L,NY,NX)+PFEPB(L,NY,NX)+PCPDB(L,NY,NX)) &
      +62.0*(PCAPM(L,NY,NX)+PCPMB(L,NY,NX)) &
      +93.0*(PCAPH(L,NY,NX)+PCPHB(L,NY,NX))
    TLPO4=TLPO4+POS+POX+POP
    UPO4(NY,NX)=UPO4(NY,NX)+POX
    UPP4(NY,NX)=UPP4(NY,NX)+POP
    !     IF(L.EQ.NU(NY,NX))THEN
    !     WRITE(*,2233)'TLPO4',I,J,NX,NY,L
    !    2,POS,POX,POP,TLPO4,TSEDER(NY,NX)
    !    3,TH1PEB(NY,NX),TH2PEB(NY,NX),XH1PB(L,NY,NX),XH2PB(L,NY,NX)
    !2233  FORMAT(A8,5I4,30F17.8)
    !     ENDIF
    !
    !     TOTAL SOC,SON,SOP
    !
    !     OMC=microbial biomass, ORC=microbial residue
    !     OQC,OQCH=DOC in micropores,macropores
    !     OQA,OQAH=acetate in micropores,macropores
    !     OHC,OHA=adsorbed SOC,acetate
    !     OSC=SOC(K=0:woody litter, K=1:non-woody litter,
    !     K=2:manure, K=3:POC, K=4:humus)
!
    DC=0.0_r8
    DN=0.0_r8
    DP=0.0_r8
    OC=0.0_r8
    ON=0.0_r8
    OP=0.0_r8
    OMCL(L,NY,NX)=0.0_r8
    OMNL(L,NY,NX)=0.0_r8

    DO 7970 K=0,jcplx1
      IF(K.LE.2)THEN
        DO 7960 N=1,7
          DO  M=1,3
            DO NGL=1,JG
              DC=DC+OMC(M,NGL,N,K,L,NY,NX)
              DN=DN+OMN(M,NGL,N,K,L,NY,NX)
              DP=DP+OMP(M,NGL,N,K,L,NY,NX)
              TOMT(NY,NX)=TOMT(NY,NX)+OMC(M,NGL,N,K,L,NY,NX)
              TONT(NY,NX)=TONT(NY,NX)+OMN(M,NGL,N,K,L,NY,NX)
              TOPT(NY,NX)=TOPT(NY,NX)+OMP(M,NGL,N,K,L,NY,NX)
              OMCL(L,NY,NX)=OMCL(L,NY,NX)+OMC(M,NGL,N,K,L,NY,NX)
              OMNL(L,NY,NX)=OMNL(L,NY,NX)+OMN(M,NGL,N,K,L,NY,NX)
            ENDDO
          enddo
7960    CONTINUE
      ELSE
        DO 7950 N=1,7
          DO  M=1,3
            DO NGL=1,JG
              OC=OC+OMC(M,NGL,N,K,L,NY,NX)
              ON=ON+OMN(M,NGL,N,K,L,NY,NX)
              OP=OP+OMP(M,NGL,N,K,L,NY,NX)
              TOMT(NY,NX)=TOMT(NY,NX)+OMC(M,NGL,N,K,L,NY,NX)
              TONT(NY,NX)=TONT(NY,NX)+OMN(M,NGL,N,K,L,NY,NX)
              TOPT(NY,NX)=TOPT(NY,NX)+OMP(M,NGL,N,K,L,NY,NX)
              OMCL(L,NY,NX)=OMCL(L,NY,NX)+OMC(M,NGL,N,K,L,NY,NX)
              OMNL(L,NY,NX)=OMNL(L,NY,NX)+OMN(M,NGL,N,K,L,NY,NX)
            enddo
          enddo
7950    CONTINUE
      ENDIF
7970  CONTINUE

        DO  N=1,7
          DO  M=1,3
            DO NGL=1,JG
              OC=OC+OMCff(M,NGL,N,L,NY,NX)
              ON=ON+OMNff(M,NGL,N,L,NY,NX)
              OP=OP+OMPff(M,NGL,N,L,NY,NX)
              TOMT(NY,NX)=TOMT(NY,NX)+OMCff(M,NGL,N,L,NY,NX)
              TONT(NY,NX)=TONT(NY,NX)+OMNff(M,NGL,N,L,NY,NX)
              TOPT(NY,NX)=TOPT(NY,NX)+OMPff(M,NGL,N,L,NY,NX)
              OMCL(L,NY,NX)=OMCL(L,NY,NX)+OMCff(M,NGL,N,L,NY,NX)
              OMNL(L,NY,NX)=OMNL(L,NY,NX)+OMNff(M,NGL,N,L,NY,NX)
            enddo
          enddo
        ENDDO

    DO 7900 K=0,jcplx1
      IF(K.LE.2)THEN
        DO 7940 M=1,2
          DC=DC+ORC(M,K,L,NY,NX)
          DN=DN+ORN(M,K,L,NY,NX)
          DP=DP+ORP(M,K,L,NY,NX)
7940    CONTINUE
        DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
          +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
        DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
        DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
        DO 7930 M=1,4
          DC=DC+OSC(M,K,L,NY,NX)
          DN=DN+OSN(M,K,L,NY,NX)
          DP=DP+OSP(M,K,L,NY,NX)
7930    CONTINUE
      ELSE
        DO 7920 M=1,2
          OC=OC+ORC(M,K,L,NY,NX)
          ON=ON+ORN(M,K,L,NY,NX)
          OP=OP+ORP(M,K,L,NY,NX)
7920    CONTINUE
        OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
          +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
        ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
        OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
        DO 7910 M=1,4
          OC=OC+OSC(M,K,L,NY,NX)
          ON=ON+OSN(M,K,L,NY,NX)
          OP=OP+OSP(M,K,L,NY,NX)
7910    CONTINUE
      ENDIF
7900  CONTINUE

    ORGC(L,NY,NX)=DC+OC
    ORGN(L,NY,NX)=DN+ON
    ORGR(L,NY,NX)=DC
    IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
      DORGC(L,NY,NX)=ORGCX(L,NY,NX)-ORGC(L,NY,NX)
      IF(L.EQ.NU(NY,NX))THEN
        DORGC(L,NY,NX)=DORGC(L,NY,NX)+DORGE(NY,NX)
      ENDIF
    ELSE
      DORGC(L,NY,NX)=0.0_r8
    ENDIF

    TLRSDC=TLRSDC+DC
    URSDC(NY,NX)=URSDC(NY,NX)+DC
    TLRSDN=TLRSDN+DN
    URSDN(NY,NX)=URSDN(NY,NX)+DN
    TLRSDP=TLRSDP+DP
    URSDP(NY,NX)=URSDP(NY,NX)+DP
    TLORGC=TLORGC+OC
    UORGC(NY,NX)=UORGC(NY,NX)+OC
    TLORGN=TLORGN+ON
    UORGN(NY,NX)=UORGN(NY,NX)+ON
    TLORGP=TLORGP+OP
    UORGP(NY,NX)=UORGP(NY,NX)+OP
    TSEDSO=TSEDSO+(DC+OC)*1.0E-06

!
!     TOTAL SALT IONS
!
    IF(ISALTG.NE.0)THEN
      ZAL(L,NY,NX)=ZAL(L,NY,NX)+TRAL(L,NY,NX)+TALFLS(L,NY,NX) &
        +RALFLU(L,NY,NX)+XALFXS(L,NY,NX)
      ZFE(L,NY,NX)=ZFE(L,NY,NX)+TRFE(L,NY,NX)+TFEFLS(L,NY,NX) &
        +RFEFLU(L,NY,NX)+XFEFXS(L,NY,NX)
      ZHY(L,NY,NX)=ZHY(L,NY,NX)+TRHY(L,NY,NX)+THYFLS(L,NY,NX) &
        +RHYFLU(L,NY,NX)+XHYFXS(L,NY,NX)+XZHYS(L,NY,NX)
      ZCA(L,NY,NX)=ZCA(L,NY,NX)+TRCA(L,NY,NX)+TCAFLS(L,NY,NX) &
        +RCAFLU(L,NY,NX)+XCAFXS(L,NY,NX)
      ZMG(L,NY,NX)=ZMG(L,NY,NX)+TRMG(L,NY,NX)+TMGFLS(L,NY,NX) &
        +RMGFLU(L,NY,NX)+XMGFXS(L,NY,NX)
      ZNA(L,NY,NX)=ZNA(L,NY,NX)+TRNA(L,NY,NX)+TNAFLS(L,NY,NX) &
        +RNAFLU(L,NY,NX)+XNAFXS(L,NY,NX)
      ZKA(L,NY,NX)=ZKA(L,NY,NX)+TRKA(L,NY,NX)+TKAFLS(L,NY,NX) &
        +RKAFLU(L,NY,NX)+XKAFXS(L,NY,NX)
      ZOH(L,NY,NX)=ZOH(L,NY,NX)+TROH(L,NY,NX)+TOHFLS(L,NY,NX) &
        +ROHFLU(L,NY,NX)+XOHFXS(L,NY,NX)

      ZSO4(L,NY,NX)=ZSO4(L,NY,NX)+TRSO4(L,NY,NX)+TSOFLS(L,NY,NX) &
        +RSOFLU(L,NY,NX)+XSOFXS(L,NY,NX)
      ZCL(L,NY,NX)=ZCL(L,NY,NX)+TCLFLS(L,NY,NX)+RCLFLU(L,NY,NX) &
        +XCLFXS(L,NY,NX)
      ZCO3(L,NY,NX)=ZCO3(L,NY,NX)+TRCO3(L,NY,NX)+TC3FLS(L,NY,NX) &
        +XC3FXS(L,NY,NX)
      ZHCO3(L,NY,NX)=ZHCO3(L,NY,NX)+TRHCO(L,NY,NX)+THCFLS(L,NY,NX) &
        +XHCFXS(L,NY,NX)
      ZALOH1(L,NY,NX)=ZALOH1(L,NY,NX)+TRAL1(L,NY,NX)+TAL1FS(L,NY,NX) &
        +XAL1XS(L,NY,NX)
      ZALOH2(L,NY,NX)=ZALOH2(L,NY,NX)+TRAL2(L,NY,NX)+TAL2FS(L,NY,NX) &
        +XAL2XS(L,NY,NX)-TRXAL2(L,NY,NX)
      ZALOH3(L,NY,NX)=ZALOH3(L,NY,NX)+TRAL3(L,NY,NX)+TAL3FS(L,NY,NX) &
        +XAL3XS(L,NY,NX)
      ZALOH4(L,NY,NX)=ZALOH4(L,NY,NX)+TRAL4(L,NY,NX)+TAL4FS(L,NY,NX) &
        +XAL4XS(L,NY,NX)
      ZALS(L,NY,NX)=ZALS(L,NY,NX)+TRALS(L,NY,NX)+TALSFS(L,NY,NX) &
        +XALSXS(L,NY,NX)
      ZFEOH1(L,NY,NX)=ZFEOH1(L,NY,NX)+TRFE1(L,NY,NX)+TFE1FS(L,NY,NX) &
        +XFE1XS(L,NY,NX)
      ZFEOH2(L,NY,NX)=ZFEOH2(L,NY,NX)+TRFE2(L,NY,NX)+TFE2FS(L,NY,NX) &
        +XFE2XS(L,NY,NX)-TRXFE2(L,NY,NX)
      ZFEOH3(L,NY,NX)=ZFEOH3(L,NY,NX)+TRFE3(L,NY,NX)+TFE3FS(L,NY,NX) &
        +XFE3XS(L,NY,NX)
      ZFEOH4(L,NY,NX)=ZFEOH4(L,NY,NX)+TRFE4(L,NY,NX)+TFE4FS(L,NY,NX) &
        +XFE4XS(L,NY,NX)
      ZFES(L,NY,NX)=ZFES(L,NY,NX)+TRFES(L,NY,NX)+TFESFS(L,NY,NX) &
        +XFESXS(L,NY,NX)
      ZCAO(L,NY,NX)=ZCAO(L,NY,NX)+TRCAO(L,NY,NX)+TCAOFS(L,NY,NX) &
        +XCAOXS(L,NY,NX)
      ZCAC(L,NY,NX)=ZCAC(L,NY,NX)+TRCAC(L,NY,NX)+TCACFS(L,NY,NX) &
        +XCACXS(L,NY,NX)
      ZCAH(L,NY,NX)=ZCAH(L,NY,NX)+TRCAH(L,NY,NX)+TCAHFS(L,NY,NX) &
        +XCAHXS(L,NY,NX)
      ZCAS(L,NY,NX)=ZCAS(L,NY,NX)+TRCAS(L,NY,NX)+TCASFS(L,NY,NX) &
        +XCASXS(L,NY,NX)
      ZMGO(L,NY,NX)=ZMGO(L,NY,NX)+TRMGO(L,NY,NX)+TMGOFS(L,NY,NX) &
        +XMGOXS(L,NY,NX)
      ZMGC(L,NY,NX)=ZMGC(L,NY,NX)+TRMGC(L,NY,NX)+TMGCFS(L,NY,NX) &
        +XMGCXS(L,NY,NX)
      ZMGH(L,NY,NX)=ZMGH(L,NY,NX)+TRMGH(L,NY,NX)+TMGHFS(L,NY,NX) &
        +XMGHXS(L,NY,NX)
      ZMGS(L,NY,NX)=ZMGS(L,NY,NX)+TRMGS(L,NY,NX)+TMGSFS(L,NY,NX) &
        +XMGSXS(L,NY,NX)
      ZNAC(L,NY,NX)=ZNAC(L,NY,NX)+TRNAC(L,NY,NX)+TNACFS(L,NY,NX) &
        +XNACXS(L,NY,NX)
      ZNAS(L,NY,NX)=ZNAS(L,NY,NX)+TRNAS(L,NY,NX)+TNASFS(L,NY,NX) &
        +XNASXS(L,NY,NX)
      ZKAS(L,NY,NX)=ZKAS(L,NY,NX)+TRKAS(L,NY,NX)+TKASFS(L,NY,NX) &
        +XKASXS(L,NY,NX)
      H0PO4(L,NY,NX)=H0PO4(L,NY,NX)+TRH0P(L,NY,NX)+TH0PFS(L,NY,NX) &
        +XH0PXS(L,NY,NX)
      H3PO4(L,NY,NX)=H3PO4(L,NY,NX)+TRH3P(L,NY,NX)+TH3PFS(L,NY,NX) &
        +XH3PXS(L,NY,NX)
      ZFE1P(L,NY,NX)=ZFE1P(L,NY,NX)+TRF1P(L,NY,NX)+TF1PFS(L,NY,NX) &
        +XF1PXS(L,NY,NX)
      ZFE2P(L,NY,NX)=ZFE2P(L,NY,NX)+TRF2P(L,NY,NX)+TF2PFS(L,NY,NX) &
        +XF2PXS(L,NY,NX)
      ZCA0P(L,NY,NX)=ZCA0P(L,NY,NX)+TRC0P(L,NY,NX)+TC0PFS(L,NY,NX) &
        +XC0PXS(L,NY,NX)
      ZCA1P(L,NY,NX)=ZCA1P(L,NY,NX)+TRC1P(L,NY,NX)+TC1PFS(L,NY,NX) &
        +XC1PXS(L,NY,NX)
      ZCA2P(L,NY,NX)=ZCA2P(L,NY,NX)+TRC2P(L,NY,NX)+TC2PFS(L,NY,NX) &
        +XC2PXS(L,NY,NX)
      ZMG1P(L,NY,NX)=ZMG1P(L,NY,NX)+TRM1P(L,NY,NX)+TM1PFS(L,NY,NX) &
        +XM1PXS(L,NY,NX)
      H0POB(L,NY,NX)=H0POB(L,NY,NX)+TRH0B(L,NY,NX)+TH0BFB(L,NY,NX) &
        +XH0BXB(L,NY,NX)
      H3POB(L,NY,NX)=H3POB(L,NY,NX)+TRH3B(L,NY,NX)+TH3BFB(L,NY,NX) &
        +XH3BXB(L,NY,NX)
      ZFE1PB(L,NY,NX)=ZFE1PB(L,NY,NX)+TRF1B(L,NY,NX)+TF1BFB(L,NY,NX) &
        +XF1BXB(L,NY,NX)
      ZFE2PB(L,NY,NX)=ZFE2PB(L,NY,NX)+TRF2B(L,NY,NX)+TF2BFB(L,NY,NX) &
        +XF2BXB(L,NY,NX)
      ZCA0PB(L,NY,NX)=ZCA0PB(L,NY,NX)+TRC0B(L,NY,NX)+TC0BFB(L,NY,NX) &
        +XC0BXB(L,NY,NX)
      ZCA1PB(L,NY,NX)=ZCA1PB(L,NY,NX)+TRC1B(L,NY,NX)+TC1BFB(L,NY,NX) &
        +XC1BXB(L,NY,NX)
      ZCA2PB(L,NY,NX)=ZCA2PB(L,NY,NX)+TRC2B(L,NY,NX)+TC2BFB(L,NY,NX) &
        +XC2BXB(L,NY,NX)
      ZMG1PB(L,NY,NX)=ZMG1PB(L,NY,NX)+TRM1B(L,NY,NX)+TM1BFB(L,NY,NX) &
        +XM1BXB(L,NY,NX)
      ZALH(L,NY,NX)=ZALH(L,NY,NX)+TALFHS(L,NY,NX)-XALFXS(L,NY,NX)
      ZFEH(L,NY,NX)=ZFEH(L,NY,NX)+TFEFHS(L,NY,NX)-XFEFXS(L,NY,NX)
      ZHYH(L,NY,NX)=ZHYH(L,NY,NX)+THYFHS(L,NY,NX)-XHYFXS(L,NY,NX)
      ZCCH(L,NY,NX)=ZCCH(L,NY,NX)+TCAFHS(L,NY,NX)-XCAFXS(L,NY,NX)
      ZMAH(L,NY,NX)=ZMAH(L,NY,NX)+TMGFHS(L,NY,NX)-XMGFXS(L,NY,NX)
      ZNAH(L,NY,NX)=ZNAH(L,NY,NX)+TNAFHS(L,NY,NX)-XNAFXS(L,NY,NX)
      ZKAH(L,NY,NX)=ZKAH(L,NY,NX)+TKAFHS(L,NY,NX)-XKAFXS(L,NY,NX)
      ZOHH(L,NY,NX)=ZOHH(L,NY,NX)+TOHFHS(L,NY,NX)-XOHFXS(L,NY,NX)
      ZSO4H(L,NY,NX)=ZSO4H(L,NY,NX)+TSOFHS(L,NY,NX)-XSOFXS(L,NY,NX)
      ZCLH(L,NY,NX)=ZCLH(L,NY,NX)+TCLFHS(L,NY,NX)-XCLFXS(L,NY,NX)
      ZCO3H(L,NY,NX)=ZCO3H(L,NY,NX)+TC3FHS(L,NY,NX)-XC3FXS(L,NY,NX)
      ZHCO3H(L,NY,NX)=ZHCO3H(L,NY,NX)+THCFHS(L,NY,NX)-XHCFXS(L,NY,NX)
      ZALO1H(L,NY,NX)=ZALO1H(L,NY,NX)+TAL1HS(L,NY,NX)-XAL1XS(L,NY,NX)
      ZALO2H(L,NY,NX)=ZALO2H(L,NY,NX)+TAL2HS(L,NY,NX)-XAL2XS(L,NY,NX)
      ZALO3H(L,NY,NX)=ZALO3H(L,NY,NX)+TAL3HS(L,NY,NX)-XAL3XS(L,NY,NX)
      ZALO4H(L,NY,NX)=ZALO4H(L,NY,NX)+TAL4HS(L,NY,NX)-XAL4XS(L,NY,NX)
      ZALSH(L,NY,NX)=ZALSH(L,NY,NX)+TALSHS(L,NY,NX)-XALSXS(L,NY,NX)
      ZFEO1H(L,NY,NX)=ZFEO1H(L,NY,NX)+TFE1HS(L,NY,NX)-XFE1XS(L,NY,NX)
      ZFEO2H(L,NY,NX)=ZFEO2H(L,NY,NX)+TFE2HS(L,NY,NX)-XFE2XS(L,NY,NX)
      ZFEO3H(L,NY,NX)=ZFEO3H(L,NY,NX)+TFE3HS(L,NY,NX)-XFE3XS(L,NY,NX)
      ZFEO4H(L,NY,NX)=ZFEO4H(L,NY,NX)+TFE4HS(L,NY,NX)-XFE4XS(L,NY,NX)
      ZFESH(L,NY,NX)=ZFESH(L,NY,NX)+TFESHS(L,NY,NX)-XFESXS(L,NY,NX)
      ZCAOH(L,NY,NX)=ZCAOH(L,NY,NX)+TCAOHS(L,NY,NX)-XCAOXS(L,NY,NX)
      ZCACH(L,NY,NX)=ZCACH(L,NY,NX)+TCACHS(L,NY,NX)-XCACXS(L,NY,NX)
      ZCAHH(L,NY,NX)=ZCAHH(L,NY,NX)+TCAHHS(L,NY,NX)-XCAHXS(L,NY,NX)
      ZCASH(L,NY,NX)=ZCASH(L,NY,NX)+TCASHS(L,NY,NX)-XCASXS(L,NY,NX)
      ZMGOH(L,NY,NX)=ZMGOH(L,NY,NX)+TMGOHS(L,NY,NX)-XMGOXS(L,NY,NX)
      ZMGCH(L,NY,NX)=ZMGCH(L,NY,NX)+TMGCHS(L,NY,NX)-XMGCXS(L,NY,NX)
      ZMGHH(L,NY,NX)=ZMGHH(L,NY,NX)+TMGHHS(L,NY,NX)-XMGHXS(L,NY,NX)
      ZMGSH(L,NY,NX)=ZMGSH(L,NY,NX)+TMGSHS(L,NY,NX)-XMGSXS(L,NY,NX)
      ZNACH(L,NY,NX)=ZNACH(L,NY,NX)+TNACHS(L,NY,NX)-XNACXS(L,NY,NX)
      ZNASH(L,NY,NX)=ZNASH(L,NY,NX)+TNASHS(L,NY,NX)-XNASXS(L,NY,NX)
      ZKASH(L,NY,NX)=ZKASH(L,NY,NX)+TKASHS(L,NY,NX)-XKASXS(L,NY,NX)
      H0PO4H(L,NY,NX)=H0PO4H(L,NY,NX)+TH0PHS(L,NY,NX)-XH0PXS(L,NY,NX)
      H3PO4H(L,NY,NX)=H3PO4H(L,NY,NX)+TH3PHS(L,NY,NX)-XH3PXS(L,NY,NX)
      ZFE1PH(L,NY,NX)=ZFE1PH(L,NY,NX)+TF1PHS(L,NY,NX)-XF1PXS(L,NY,NX)
      ZFE2PH(L,NY,NX)=ZFE2PH(L,NY,NX)+TF2PHS(L,NY,NX)-XF2PXS(L,NY,NX)
      ZCA0PH(L,NY,NX)=ZCA0PH(L,NY,NX)+TC0PHS(L,NY,NX)-XC0PXS(L,NY,NX)
      ZCA1PH(L,NY,NX)=ZCA1PH(L,NY,NX)+TC1PHS(L,NY,NX)-XC1PXS(L,NY,NX)
      ZCA2PH(L,NY,NX)=ZCA2PH(L,NY,NX)+TC2PHS(L,NY,NX)-XC2PXS(L,NY,NX)
      ZMG1PH(L,NY,NX)=ZMG1PH(L,NY,NX)+TM1PHS(L,NY,NX)-XM1PXS(L,NY,NX)
      H0POBH(L,NY,NX)=H0POBH(L,NY,NX)+TH0BHB(L,NY,NX)-XH0BXB(L,NY,NX)
      H3POBH(L,NY,NX)=H3POBH(L,NY,NX)+TH3BHB(L,NY,NX)-XH3BXB(L,NY,NX)
      ZFE1BH(L,NY,NX)=ZFE1BH(L,NY,NX)+TF1BHB(L,NY,NX)-XF1BXB(L,NY,NX)
      ZFE2BH(L,NY,NX)=ZFE2BH(L,NY,NX)+TF2BHB(L,NY,NX)-XF2BXB(L,NY,NX)
      ZCA0BH(L,NY,NX)=ZCA0BH(L,NY,NX)+TC0BHB(L,NY,NX)-XC0BXB(L,NY,NX)
      ZCA1BH(L,NY,NX)=ZCA1BH(L,NY,NX)+TC1BHB(L,NY,NX)-XC1BXB(L,NY,NX)
      ZCA2BH(L,NY,NX)=ZCA2BH(L,NY,NX)+TC2BHB(L,NY,NX)-XC2BXB(L,NY,NX)
      ZMG1BH(L,NY,NX)=ZMG1BH(L,NY,NX)+TM1BHB(L,NY,NX)-XM1BXB(L,NY,NX)
      XHY(L,NY,NX)=XHY(L,NY,NX)+TRXHY(L,NY,NX)
      XAL(L,NY,NX)=XAL(L,NY,NX)+TRXAL(L,NY,NX)
      XFE(L,NY,NX)=XFE(L,NY,NX)+TRXFE(L,NY,NX)
      XCA(L,NY,NX)=XCA(L,NY,NX)+TRXCA(L,NY,NX)
      XMG(L,NY,NX)=XMG(L,NY,NX)+TRXMG(L,NY,NX)
      XNA(L,NY,NX)=XNA(L,NY,NX)+TRXNA(L,NY,NX)
      XKA(L,NY,NX)=XKA(L,NY,NX)+TRXKA(L,NY,NX)
      XHC(L,NY,NX)=XHC(L,NY,NX)+TRXHC(L,NY,NX)
      XALO2(L,NY,NX)=XALO2(L,NY,NX)+TRXAL2(L,NY,NX)
      XFEO2(L,NY,NX)=XFEO2(L,NY,NX)+TRXFE2(L,NY,NX)
      PALOH(L,NY,NX)=PALOH(L,NY,NX)+TRALOH(L,NY,NX)
      PFEOH(L,NY,NX)=PFEOH(L,NY,NX)+TRFEOH(L,NY,NX)
      PCACO(L,NY,NX)=PCACO(L,NY,NX)+TRCACO(L,NY,NX)
      PCASO(L,NY,NX)=PCASO(L,NY,NX)+TRCASO(L,NY,NX)

      PSS=31.0*(H0PO4(L,NY,NX)+H3PO4(L,NY,NX)+ZFE1P(L,NY,NX) &
        +ZFE2P(L,NY,NX)+ZCA0P(L,NY,NX)+ZCA1P(L,NY,NX) &
        +ZCA2P(L,NY,NX)+ZMG1P(L,NY,NX)+H0POB(L,NY,NX) &
        +H3POB(L,NY,NX)+ZFE1PB(L,NY,NX)+ZFE2PB(L,NY,NX) &
        +ZCA0PB(L,NY,NX)+ZCA1PB(L,NY,NX)+ZCA2PB(L,NY,NX) &
        +ZMG1PB(L,NY,NX)+H0PO4H(L,NY,NX)+H3PO4H(L,NY,NX) &
        +ZFE1PH(L,NY,NX)+ZFE2PH(L,NY,NX)+ZCA0PH(L,NY,NX) &
        +ZCA1PH(L,NY,NX)+ZCA2PH(L,NY,NX)+ZMG1PH(L,NY,NX) &
        +H0POBH(L,NY,NX)+H3POBH(L,NY,NX)+ZFE1BH(L,NY,NX) &
        +ZFE2BH(L,NY,NX)+ZCA0BH(L,NY,NX)+ZCA1BH(L,NY,NX) &
        +ZCA2BH(L,NY,NX)+ZMG1BH(L,NY,NX))
      TLPO4=TLPO4+PSS

      SSS=ZAL(L,NY,NX)+ZFE(L,NY,NX)+ZHY(L,NY,NX)+ZCA(L,NY,NX) &
        +ZMG(L,NY,NX)+ZNA(L,NY,NX)+ZKA(L,NY,NX)+ZOH(L,NY,NX) &
        +ZSO4(L,NY,NX)+ZCL(L,NY,NX)+ZCO3(L,NY,NX)+H0PO4(L,NY,NX) &
        +H0POB(L,NY,NX) &
        +2.0_r8*(ZHCO3(L,NY,NX)+ZALOH1(L,NY,NX) &
        +ZALS(L,NY,NX)+ZFEOH1(L,NY,NX)+ZFES(L,NY,NX)+ZCAO(L,NY,NX) &
        +ZCAC(L,NY,NX)+ZCAS(L,NY,NX)+ZMGO(L,NY,NX)+ZMGC(L,NY,NX) &
        +ZMGS(L,NY,NX)+ZNAC(L,NY,NX)+ZNAS(L,NY,NX)+ZKAS(L,NY,NX) &
        +ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX)) &
        +3.0_r8*(ZALOH2(L,NY,NX)+ZFEOH2(L,NY,NX)+ZCAH(L,NY,NX) &
        +ZMGH(L,NY,NX)+ZFE1P(L,NY,NX)+ZCA1P(L,NY,NX)+ZMG1P(L,NY,NX) &
        +ZFE1PB(L,NY,NX)+ZCA1PB(L,NY,NX)+ZMG1PB(L,NY,NX)) &
        +4.0*(ZALOH3(L,NY,NX)+ZFEOH3(L,NY,NX)+H3PO4(L,NY,NX) &
        +ZFE2P(L,NY,NX)+ZCA2P(L,NY,NX)+H3POB(L,NY,NX)+ZFE2PB(L,NY,NX) &
        +ZCA2PB(L,NY,NX)) &
        +5.0*(ZALOH4(L,NY,NX)+ZFEOH4(L,NY,NX))
      SSH=ZALH(L,NY,NX)+ZFEH(L,NY,NX)+ZHYH(L,NY,NX)+ZCCH(L,NY,NX) &
        +ZMAH(L,NY,NX)+ZNAH(L,NY,NX)+ZKAH(L,NY,NX)+ZOHH(L,NY,NX) &
        +ZSO4H(L,NY,NX)+ZCLH(L,NY,NX)+ZCO3H(L,NY,NX) +H0PO4H(L,NY,NX) &
        +H0POBH(L,NY,NX) &
        +2.0*(ZHCO3H(L,NY,NX)+ZALO1H(L,NY,NX) &
        +ZALSH(L,NY,NX)+ZFEO1H(L,NY,NX)+ZFESH(L,NY,NX)+ZCAOH(L,NY,NX) &
        +ZCACH(L,NY,NX)+ZCASH(L,NY,NX)+ZMGOH(L,NY,NX)+ZMGCH(L,NY,NX) &
        +ZMGSH(L,NY,NX)+ZNACH(L,NY,NX)+ZNASH(L,NY,NX)+ZKASH(L,NY,NX) &
        +ZCA0PH(L,NY,NX)+ZCA0BH(L,NY,NX)) &
        +3.0*(ZALO2H(L,NY,NX)+ZFEO2H(L,NY,NX)+ZCAHH(L,NY,NX) &
        +ZMGHH(L,NY,NX)+ZFE1PH(L,NY,NX)+ZCA1PH(L,NY,NX)+ZMG1PH(L,NY,NX) &
        +ZFE1BH(L,NY,NX)+ZCA1BH(L,NY,NX)+ZMG1BH(L,NY,NX)) &
        +4.0*(ZALO3H(L,NY,NX)+ZFEO3H(L,NY,NX)+H3PO4H(L,NY,NX) &
        +ZFE2PH(L,NY,NX)+ZCA2PH(L,NY,NX)+H3POBH(L,NY,NX) &
        +ZFE2BH(L,NY,NX)+ZCA2BH(L,NY,NX)) &
        +5.0*(ZALO4H(L,NY,NX)+ZFEO4H(L,NY,NX))
!
!     TOTAL FERILIZER,EXCHANGEABLE CATIONS AND ANIONS, PRECIPITATES
!
      SSF=ZNH3FA(L,NY,NX)+ZNHUFA(L,NY,NX)+ZNO3FA(L,NY,NX) &
        +ZNH3FB(L,NY,NX)+ZNHUFB(L,NY,NX)+ZNO3FB(L,NY,NX) &
        +2.0*(ZNH4FA(L,NY,NX)+ZNH4FB(L,NY,NX))
      SSX=XHY(L,NY,NX)+XAL(L,NY,NX) &
        +XFE(L,NY,NX)+XCA(L,NY,NX)+XMG(L,NY,NX) &
        +XNA(L,NY,NX)+XKA(L,NY,NX)+XHC(L,NY,NX) &
        +XOH0(L,NY,NX)+XOH0B(L,NY,NX) &
        +2.0*(XN4(L,NY,NX)+XNB(L,NY,NX) &
        +XOH1(L,NY,NX)+XOH1B(L,NY,NX)) &
        +3.0*(XALO2(L,NY,NX)+XFEO2(L,NY,NX) &
        +XOH2(L,NY,NX)+XOH2B(L,NY,NX) &
        +XH1P(L,NY,NX)+XH1PB(L,NY,NX)) &
        +4.0*(XH2P(L,NY,NX)+XH2PB(L,NY,NX))
      SSP=2.0*(PCACO(L,NY,NX)+PCASO(L,NY,NX) &
        +PALPO(L,NY,NX)+PFEPO(L,NY,NX) &
        +PALPB(L,NY,NX)+PFEPB(L,NY,NX)) &
        +3.0*(PCAPD(L,NY,NX)+PCPDB(L,NY,NX)) &
        +4.0*(PALOH(L,NY,NX)+PFEOH(L,NY,NX)) &
        +7.0*(PCAPM(L,NY,NX)+PCPMB(L,NY,NX)) &
        +9.0*(PCAPH(L,NY,NX)+PCPHB(L,NY,NX))
      SST=SSS+SSH+SSF+SSX+SSP
      TION=TION+SST
      UION(NY,NX)=UION(NY,NX)+SST

!
!     SOIL ELECTRICAL CONDUCTIVITY
!
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      ECHY=0.337*AMAX1(0.0,ZHY(L,NY,NX)/VOLW(L,NY,NX))
      ECOH=0.192*AMAX1(0.0,ZOH(L,NY,NX)/VOLW(L,NY,NX))
      ECAL=0.056*AMAX1(0.0,ZAL(L,NY,NX)*3.0/VOLW(L,NY,NX))
      ECFE=0.051*AMAX1(0.0,ZFE(L,NY,NX)*3.0/VOLW(L,NY,NX))
      ECCA=0.060*AMAX1(0.0,ZCA(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECMG=0.053*AMAX1(0.0,ZMG(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECNA=0.050*AMAX1(0.0,ZNA(L,NY,NX)/VOLW(L,NY,NX))
      ECKA=0.070*AMAX1(0.0,ZKA(L,NY,NX)/VOLW(L,NY,NX))
      ECCO=0.072*AMAX1(0.0,ZCO3(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECHC=0.044*AMAX1(0.0,ZHCO3(L,NY,NX)/VOLW(L,NY,NX))
      ECSO=0.080*AMAX1(0.0,ZSO4(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECCL=0.076*AMAX1(0.0,ZCL(L,NY,NX)/VOLW(L,NY,NX))
      ECNO=0.071*AMAX1(0.0,ZNO3S(L,NY,NX)/(VOLW(L,NY,NX)*14.0))
      ECND(L,NY,NX)=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA &
      +ECCO+ECHC+ECSO+ECCL+ECNO

      ELSE
      ECND(L,NY,NX)=0.0_r8
      ENDIF
      ENDIF

125   CONTINUE
  end subroutine UpdateChemInLayers
!------------------------------------------------------------------------------------------

  subroutine SnowpackLayering(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: FX,FY
  integer :: L,L1,L0
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
!     begin_execution
!
  IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
    DO 325 L=1,JS-1
      VOLSLX=VOLSL(L,NY,NX)
      IF(VOLSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        DDLYXS=(VOLSI(L,NY,NX)-VOLSSL(L,NY,NX)/DENSS(L,NY,NX) &
          -VOLWSL(L,NY,NX)-VOLISL(L,NY,NX))/AREA(3,L,NY,NX)
        DDLYXX=DDLYXS
        IF(DDLYXS.LT.-ZERO.OR.DLYRS(L+1,NY,NX).GT.ZERO)THEN
          DDLYRS=AMIN1(DDLYXS,DLYRS(L+1,NY,NX))
          IFLGLS=1
        ELSE
          DDLYXS=(VOLSL(L,NY,NX)-VOLSSL(L,NY,NX)/DENSS(L,NY,NX) &
            -VOLWSL(L,NY,NX)-VOLISL(L,NY,NX))/AREA(3,L,NY,NX)
          DDLYRS=DDLYXS
          IFLGLS=2
        ENDIF
      ELSE
        DDLYRS=0.0_r8
        IFLGLS=0
      ENDIF
      !
      !     RESET SNOW LAYER DEPTHS
      !
      CDPTHS(L,NY,NX)=CDPTHS(L,NY,NX)+DDLYRS
      DLYRS(L,NY,NX)=CDPTHS(L,NY,NX)-CDPTHS(L-1,NY,NX)
!
      !     TRANSFER STATE VARIABLES BETWEEN LAYERS
      !
      IF(ABS(DDLYRS).GT.ZERO)THEN
        !     WRITE(*,1113)'DDLYRS',I,J,NX,NY,L,IFLGLS,DDLYRS,DDLYXS
        !    2,CDPTHS(L,NY,NX),DLYRS(L,NY,NX),VOLSI(L,NY,NX),VOLSL(L,NY,NX)
        !    3,VOLSLX,DDLYXX,CDPTHS(L+1,NY,NX),DLYRS(L+1,NY,NX)
!1113  FORMAT(A8,6I4,30E14.6)
        IF(DDLYRS.GT.0.0)THEN
          L1=L
          L0=L+1
          IF(DDLYRS.LT.DDLYXS)THEN
            FX=1.0
          ELSE
            FX=AMIN1(1.0,DDLYRS*AREA(3,L0,NY,NX)/VOLSL(L0,NY,NX))
          ENDIF
        ELSE
          L1=L+1
          L0=L
          IF(VOLSL(L0,NY,NX).LT.VOLSI(L0,NY,NX))THEN
            FX=0.0_r8
          ELSE
            FX=AMIN1(1.0,-DDLYRS*AREA(3,L0,NY,NX)/VOLSL(L0,NY,NX))
          ENDIF
        ENDIF
      IF(FX.GT.0.0)THEN
      FY=1.0-FX
!     IF(IYRC.EQ.2006.AND.I.EQ.361.AND.NX.EQ.1)THEN
!     WRITE(*,5596)'SNOW1',I,J,NX,NY,L,NU(NY,NX),L1,L0,FX,FY
!    3,DDLYRS,VOLSI(L0,NY,NX),VOLSL(L0,NY,NX),VOLSSL(L0,NY,NX)
!    3,VOLWSL(L0,NY,NX),VOLISL(L0,NY,NX),VOLSI(L1,NY,NX)
!    4,VOLSL(L1,NY,NX),VOLSSL(L1,NY,NX),VOLWSL(L1,NY,NX)
!    4,VOLISL(L1,NY,NX),CDPTHS(L0,NY,NX),CDPTHS(L1,NY,NX)
!    5,DENSS(L1,NY,NX),DENSS(L0,NY,NX)
!    5,TKW(L0,NY,NX),TKW(L1,NY,NX)
!    5,VHCPW(L0,NY,NX),VHCPW(L1,NY,NX)
!    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
!    5,VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
!    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
!    5+VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
!5596  FORMAT(A8,8I4,100E14.6)
!     ENDIF
!
!     TARGET SNOW LAYER
!
      VOLSSL(L1,NY,NX)=VOLSSL(L1,NY,NX)+FX*VOLSSL(L0,NY,NX)
      VOLWSL(L1,NY,NX)=VOLWSL(L1,NY,NX)+FX*VOLWSL(L0,NY,NX)
      VOLISL(L1,NY,NX)=VOLISL(L1,NY,NX)+FX*VOLISL(L0,NY,NX)
      VOLSL(L1,NY,NX)=VOLSSL(L1,NY,NX)/DENSS(L1,NY,NX) &
      +VOLWSL(L1,NY,NX)+VOLISL(L1,NY,NX)
      ENGY1X=VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
      ENGY0X=VHCPW(L0,NY,NX)*TKW(L0,NY,NX)
      ENGY1=ENGY1X+FX*ENGY0X
      VHCPW(L1,NY,NX)=cps*VOLSSL(L1,NY,NX)+cpw*VOLWSL(L1,NY,NX) &
      +cpi*VOLISL(L1,NY,NX)
      IF(VHCPW(L1,NY,NX).GT.ZEROS(NY,NX))THEN
      TKW(L1,NY,NX)=ENGY1/VHCPW(L1,NY,NX)
      ELSE
      TKW(L1,NY,NX)=TKW(L0,NY,NX)
      ENDIF
      TCW(L1,NY,NX)=TKW(L1,NY,NX)-TC2K
      CO2W(L1,NY,NX)=CO2W(L1,NY,NX)+FX*CO2W(L0,NY,NX)
      CH4W(L1,NY,NX)=CH4W(L1,NY,NX)+FX*CH4W(L0,NY,NX)
      OXYW(L1,NY,NX)=OXYW(L1,NY,NX)+FX*OXYW(L0,NY,NX)
      ZNGW(L1,NY,NX)=ZNGW(L1,NY,NX)+FX*ZNGW(L0,NY,NX)
      ZN2W(L1,NY,NX)=ZN2W(L1,NY,NX)+FX*ZN2W(L0,NY,NX)
      ZN4W(L1,NY,NX)=ZN4W(L1,NY,NX)+FX*ZN4W(L0,NY,NX)
      ZN3W(L1,NY,NX)=ZN3W(L1,NY,NX)+FX*ZN3W(L0,NY,NX)
      ZNOW(L1,NY,NX)=ZNOW(L1,NY,NX)+FX*ZNOW(L0,NY,NX)
      Z1PW(L1,NY,NX)=Z1PW(L1,NY,NX)+FX*Z1PW(L0,NY,NX)
      ZHPW(L1,NY,NX)=ZHPW(L1,NY,NX)+FX*ZHPW(L0,NY,NX)
      IF(ISALTG.NE.0)THEN
      ZALW(L1,NY,NX)=ZALW(L1,NY,NX)+FX*ZALW(L0,NY,NX)
      ZFEW(L1,NY,NX)=ZFEW(L1,NY,NX)+FX*ZFEW(L0,NY,NX)
      ZHYW(L1,NY,NX)=ZHYW(L1,NY,NX)+FX*ZHYW(L0,NY,NX)
      ZCAW(L1,NY,NX)=ZCAW(L1,NY,NX)+FX*ZCAW(L0,NY,NX)
      ZMGW(L1,NY,NX)=ZMGW(L1,NY,NX)+FX*ZMGW(L0,NY,NX)
      ZNAW(L1,NY,NX)=ZNAW(L1,NY,NX)+FX*ZNAW(L0,NY,NX)
      ZKAW(L1,NY,NX)=ZKAW(L1,NY,NX)+FX*ZKAW(L0,NY,NX)
      ZOHW(L1,NY,NX)=ZOHW(L1,NY,NX)+FX*ZOHW(L0,NY,NX)
      ZSO4W(L1,NY,NX)=ZSO4W(L1,NY,NX)+FX*ZSO4W(L0,NY,NX)
      ZCLW(L1,NY,NX)=ZCLW(L1,NY,NX)+FX*ZCLW(L0,NY,NX)
      ZCO3W(L1,NY,NX)=ZCO3W(L1,NY,NX)+FX*ZCO3W(L0,NY,NX)
      ZHCO3W(L1,NY,NX)=ZHCO3W(L1,NY,NX)+FX*ZHCO3W(L0,NY,NX)
      ZALH1W(L1,NY,NX)=ZALH1W(L1,NY,NX)+FX*ZALH1W(L0,NY,NX)
      ZALH2W(L1,NY,NX)=ZALH2W(L1,NY,NX)+FX*ZALH2W(L0,NY,NX)
      ZALH3W(L1,NY,NX)=ZALH3W(L1,NY,NX)+FX*ZALH3W(L0,NY,NX)
      ZALH4W(L1,NY,NX)=ZALH4W(L1,NY,NX)+FX*ZALH4W(L0,NY,NX)
      ZALSW(L1,NY,NX)=ZALSW(L1,NY,NX)+FX*ZALSW(L0,NY,NX)
      ZFEH1W(L1,NY,NX)=ZFEH1W(L1,NY,NX)+FX*ZFEH1W(L0,NY,NX)
      ZFEH2W(L1,NY,NX)=ZFEH2W(L1,NY,NX)+FX*ZFEH2W(L0,NY,NX)
      ZFEH3W(L1,NY,NX)=ZFEH3W(L1,NY,NX)+FX*ZFEH3W(L0,NY,NX)
      ZFEH4W(L1,NY,NX)=ZFEH4W(L1,NY,NX)+FX*ZFEH4W(L0,NY,NX)
      ZFESW(L1,NY,NX)=ZFESW(L1,NY,NX)+FX*ZFESW(L0,NY,NX)
      ZCAOW(L1,NY,NX)=ZCAOW(L1,NY,NX)+FX*ZCAOW(L0,NY,NX)
      ZCACW(L1,NY,NX)=ZCACW(L1,NY,NX)+FX*ZCACW(L0,NY,NX)
      ZCAHW(L1,NY,NX)=ZCAHW(L1,NY,NX)+FX*ZCAHW(L0,NY,NX)
      ZCASW(L1,NY,NX)=ZCASW(L1,NY,NX)+FX*ZCASW(L0,NY,NX)
      ZMGOW(L1,NY,NX)=ZMGOW(L1,NY,NX)+FX*ZMGOW(L0,NY,NX)
      ZMGCW(L1,NY,NX)=ZMGCW(L1,NY,NX)+FX*ZMGCW(L0,NY,NX)
      ZMGHW(L1,NY,NX)=ZMGHW(L1,NY,NX)+FX*ZMGHW(L0,NY,NX)
      ZMGSW(L1,NY,NX)=ZMGSW(L1,NY,NX)+FX*ZMGSW(L0,NY,NX)
      ZNACW(L1,NY,NX)=ZNACW(L1,NY,NX)+FX*ZNACW(L0,NY,NX)
      ZNASW(L1,NY,NX)=ZNASW(L1,NY,NX)+FX*ZNASW(L0,NY,NX)
      ZKASW(L1,NY,NX)=ZKASW(L1,NY,NX)+FX*ZKASW(L0,NY,NX)
      H0PO4W(L1,NY,NX)=H0PO4W(L1,NY,NX)+FX*H0PO4W(L0,NY,NX)
      H3PO4W(L1,NY,NX)=H3PO4W(L1,NY,NX)+FX*H3PO4W(L0,NY,NX)
      ZFE1PW(L1,NY,NX)=ZFE1PW(L1,NY,NX)+FX*ZFE1PW(L0,NY,NX)
      ZFE2PW(L1,NY,NX)=ZFE2PW(L1,NY,NX)+FX*ZFE2PW(L0,NY,NX)
      ZCA0PW(L1,NY,NX)=ZCA0PW(L1,NY,NX)+FX*ZCA0PW(L0,NY,NX)
      ZCA1PW(L1,NY,NX)=ZCA1PW(L1,NY,NX)+FX*ZCA1PW(L0,NY,NX)
      ZCA2PW(L1,NY,NX)=ZCA2PW(L1,NY,NX)+FX*ZCA2PW(L0,NY,NX)
      ZMG1PW(L1,NY,NX)=ZMG1PW(L1,NY,NX)+FX*ZMG1PW(L0,NY,NX)
      ENDIF
!
!     SOURCE SNOW LAYER
!
      VOLSSL(L0,NY,NX)=FY*VOLSSL(L0,NY,NX)
      VOLWSL(L0,NY,NX)=FY*VOLWSL(L0,NY,NX)
      VOLISL(L0,NY,NX)=FY*VOLISL(L0,NY,NX)
      VOLSL(L0,NY,NX)=VOLSSL(L0,NY,NX)/DENSS(L0,NY,NX) &
      +VOLWSL(L0,NY,NX)+VOLISL(L0,NY,NX)
      ENGY0=FY*ENGY0X
      VHCPW(L0,NY,NX)=cps*VOLSSL(L0,NY,NX)+cpw*VOLWSL(L0,NY,NX) &
      +cpi*VOLISL(L0,NY,NX)
      IF(VHCPW(L0,NY,NX).GT.ZEROS(NY,NX))THEN
      TKW(L0,NY,NX)=ENGY0/VHCPW(L0,NY,NX)
      ELSE
      TKW(L0,NY,NX)=TKW(L1,NY,NX)
      ENDIF
      TCW(L0,NY,NX)=TKW(L0,NY,NX)-TC2K
      CO2W(L0,NY,NX)=FY*CO2W(L0,NY,NX)
      CH4W(L0,NY,NX)=FY*CH4W(L0,NY,NX)
      OXYW(L0,NY,NX)=FY*OXYW(L0,NY,NX)
      ZNGW(L0,NY,NX)=FY*ZNGW(L0,NY,NX)
      ZN2W(L0,NY,NX)=FY*ZN2W(L0,NY,NX)
      ZN4W(L0,NY,NX)=FY*ZN4W(L0,NY,NX)
      ZN3W(L0,NY,NX)=FY*ZN3W(L0,NY,NX)
      ZNOW(L0,NY,NX)=FY*ZNOW(L0,NY,NX)
      Z1PW(L0,NY,NX)=FY*Z1PW(L0,NY,NX)
      ZHPW(L0,NY,NX)=FY*ZHPW(L0,NY,NX)
      IF(ISALTG.NE.0)THEN
      ZALW(L0,NY,NX)=FY*ZALW(L0,NY,NX)
      ZFEW(L0,NY,NX)=FY*ZFEW(L0,NY,NX)
      ZHYW(L0,NY,NX)=FY*ZHYW(L0,NY,NX)
      ZCAW(L0,NY,NX)=FY*ZCAW(L0,NY,NX)
      ZMGW(L0,NY,NX)=FY*ZMGW(L0,NY,NX)
      ZNAW(L0,NY,NX)=FY*ZNAW(L0,NY,NX)
      ZKAW(L0,NY,NX)=FY*ZKAW(L0,NY,NX)
      ZOHW(L0,NY,NX)=FY*ZOHW(L0,NY,NX)
      ZSO4W(L0,NY,NX)=FY*ZSO4W(L0,NY,NX)
      ZCLW(L0,NY,NX)=FY*ZCLW(L0,NY,NX)
      ZCO3W(L0,NY,NX)=FY*ZCO3W(L0,NY,NX)
      ZHCO3W(L0,NY,NX)=FY*ZHCO3W(L0,NY,NX)
      ZALH1W(L0,NY,NX)=FY*ZALH1W(L0,NY,NX)
      ZALH2W(L0,NY,NX)=FY*ZALH2W(L0,NY,NX)
      ZALH3W(L0,NY,NX)=FY*ZALH3W(L0,NY,NX)
      ZALH4W(L0,NY,NX)=FY*ZALH4W(L0,NY,NX)
      ZALSW(L0,NY,NX)=FY*ZALSW(L0,NY,NX)
      ZFEH1W(L0,NY,NX)=FY*ZFEH1W(L0,NY,NX)
      ZFEH2W(L0,NY,NX)=FY*ZFEH2W(L0,NY,NX)
      ZFEH3W(L0,NY,NX)=FY*ZFEH3W(L0,NY,NX)
      ZFEH4W(L0,NY,NX)=FY*ZFEH4W(L0,NY,NX)
      ZFESW(L0,NY,NX)=FY*ZFESW(L0,NY,NX)
      ZCAOW(L0,NY,NX)=FY*ZCAOW(L0,NY,NX)
      ZCACW(L0,NY,NX)=FY*ZCACW(L0,NY,NX)
      ZCAHW(L0,NY,NX)=FY*ZCAHW(L0,NY,NX)
      ZCASW(L0,NY,NX)=FY*ZCASW(L0,NY,NX)
      ZMGOW(L0,NY,NX)=FY*ZMGOW(L0,NY,NX)
      ZMGCW(L0,NY,NX)=FY*ZMGCW(L0,NY,NX)
      ZMGHW(L0,NY,NX)=FY*ZMGHW(L0,NY,NX)
      ZMGSW(L0,NY,NX)=FY*ZMGSW(L0,NY,NX)
      ZNACW(L0,NY,NX)=FY*ZNACW(L0,NY,NX)
      ZNASW(L0,NY,NX)=FY*ZNASW(L0,NY,NX)
      ZKASW(L0,NY,NX)=FY*ZKASW(L0,NY,NX)
      H0PO4W(L0,NY,NX)=FY*H0PO4W(L0,NY,NX)
      H3PO4W(L0,NY,NX)=FY*H3PO4W(L0,NY,NX)
      ZFE1PW(L0,NY,NX)=FY*ZFE1PW(L0,NY,NX)
      ZFE2PW(L0,NY,NX)=FY*ZFE2PW(L0,NY,NX)
      ZCA0PW(L0,NY,NX)=FY*ZCA0PW(L0,NY,NX)
      ZCA1PW(L0,NY,NX)=FY*ZCA1PW(L0,NY,NX)
      ZCA2PW(L0,NY,NX)=FY*ZCA2PW(L0,NY,NX)
      ZMG1PW(L0,NY,NX)=FY*ZMG1PW(L0,NY,NX)
      ENDIF
!     IF(VOLWSL(L0,NY,NX)+VOLISL(L0,NY,NX)
!    2+VOLSSL(L0,NY,NX).LE.ZEROS(NY,NX))THEN
!     CDPTHS(L1,NY,NX)=CDPTHS(L0,NY,NX)
!     ENDIF
!     IF(IYRC.EQ.2006.AND.I.EQ.361.AND.NX.EQ.1)THEN
!     WRITE(*,5596)'SNOW2',I,J,NX,NY,L,NU(NY,NX),L1,L0,FX,FY
!    3,DDLYRS,VOLSI(L0,NY,NX),VOLSL(L0,NY,NX),VOLSSL(L0,NY,NX)
!    3,VOLWSL(L0,NY,NX),VOLISL(L0,NY,NX),VOLSI(L1,NY,NX)
!    4,VOLSL(L1,NY,NX),VOLSSL(L1,NY,NX),VOLWSL(L1,NY,NX)
!    4,VOLISL(L1,NY,NX),CDPTHS(L0,NY,NX),CDPTHS(L1,NY,NX)
!    5,DENSS(L1,NY,NX),DENSS(L0,NY,NX)
!    5,TKW(L0,NY,NX),TKW(L1,NY,NX)
!    5,VHCPW(L0,NY,NX),VHCPW(L1,NY,NX)
!    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
!    5,VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
!    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
!    5+VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
!     ENDIF
      ENDIF
      ENDIF
325   CONTINUE
      ENDIF
      end subroutine SnowpackLayering
!------------------------------------------------------------------------------------------

      subroutine ApplyMixing(I,J,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NY,NX

      integer :: K,N,M,L,LL,NGL
!     begin_execution
!
      IF(J.EQ.INT(ZNOON(NY,NX)).AND.XCORP(NY,NX).LT.1.0 &
      .AND.DCORP(I,NY,NX).GT.0.0)THEN
!
!     EXTENT OF MIXING
!
      IFLGS(NY,NX)=1
      CORP=1.0-XCORP(NY,NX)
      ENGYP(NY,NX)=0.0_r8
!
!     TEMPORARY ACCUMULATORS
!
      TBKDX=0.0_r8
      TFC=0.0_r8
      TWP=0.0_r8
      TSCNV=0.0_r8
      TSCNH=0.0_r8
      TSAND=0.0_r8
      TSILT=0.0_r8
      TCLAY=0.0_r8
      TXCEC=0.0_r8
      TXAEC=0.0_r8
      TGKC4=0.0_r8
      TGKCA=0.0_r8
      TGKCM=0.0_r8
      TGKCN=0.0_r8
      TGKCK=0.0_r8
      TVOLW=0.0_r8
      TVOLI=0.0_r8
!     TVOLP=0.0_r8
!     TVOLA=0.0_r8
      TENGY=0.0_r8
      TVHCM=0.0_r8
      TNFNIH=0.0_r8
      TNH4FA=0.0_r8
      TNH3FA=0.0_r8
      TNHUFA=0.0_r8
      TNO3FA=0.0_r8
      TNH4FB=0.0_r8
      TNH3FB=0.0_r8
      TNHUFB=0.0_r8
      TNO3FB=0.0_r8
      TNH4S=0.0_r8
      TNH4B=0.0_r8
      TNH3S=0.0_r8
      TNH3B=0.0_r8
      TNO3S=0.0_r8
      TNO3B=0.0_r8
      TNO2S=0.0_r8
      TNO2B=0.0_r8
      TZAL=0.0_r8
      TZFE=0.0_r8
      TZHY=0.0_r8
      TZCA=0.0_r8
      TZMG=0.0_r8
      TZNA=0.0_r8
      TZKA=0.0_r8
      TZOH=0.0_r8
      TZSO4=0.0_r8
      TZCL=0.0_r8
      TZCO3=0.0_r8
      TZHCO3=0.0_r8
      TZALO1=0.0_r8
      TZALO2=0.0_r8
      TZALO3=0.0_r8
      TZALO4=0.0_r8
      TZALS=0.0_r8
      TZFEO1=0.0_r8
      TZFEO2=0.0_r8
      TZFEO3=0.0_r8
      TZFEO4=0.0_r8
      TZFES=0.0_r8
      TZCAO=0.0_r8
      TZCAC=0.0_r8
      TZCAH=0.0_r8
      TZCAS=0.0_r8
      TZMGO=0.0_r8
      TZMGC=0.0_r8
      TZMGH=0.0_r8
      TZMGS=0.0_r8
      TZNAC=0.0_r8
      TZNAS=0.0_r8
      TZKAS=0.0_r8
      TH0PO4=0.0_r8
      TH1PO4=0.0_r8
      TH2PO4=0.0_r8
      TH3PO4=0.0_r8
      TZFE1P=0.0_r8
      TZFE2P=0.0_r8
      TZCA0P=0.0_r8
      TZCA1P=0.0_r8
      TZCA2P=0.0_r8
      TZMG1P=0.0_r8
      TH0POB=0.0_r8
      TH1POB=0.0_r8
      TH2POB=0.0_r8
      TH3POB=0.0_r8
      TFE1PB=0.0_r8
      TFE2PB=0.0_r8
      TCA0PB=0.0_r8
      TCA1PB=0.0_r8
      TCA2PB=0.0_r8
      TMG1PB=0.0_r8
      TXNH4=0.0_r8
      TXNHB=0.0_r8
      TXHY=0.0_r8
      TXAL=0.0_r8
      TXFE=0.0_r8
      TXCA=0.0_r8
      TXMG=0.0_r8
      TXNA=0.0_r8
      TXKA=0.0_r8
      TXHC=0.0_r8
      TXAL2=0.0_r8
      TXFE2=0.0_r8
      TXOH0=0.0_r8
      TXOH1=0.0_r8
      TXOH2=0.0_r8
      TXH1P=0.0_r8
      TXH2P=0.0_r8
      TXOH0B=0.0_r8
      TXOH1B=0.0_r8
      TXOH2B=0.0_r8
      TXH1PB=0.0_r8
      TXH2PB=0.0_r8
      TPALOH=0.0_r8
      TPFEOH=0.0_r8
      TPCACO=0.0_r8
      TPCASO=0.0_r8
      TPALPO=0.0_r8
      TPFEPO=0.0_r8
      TPCAPD=0.0_r8
      TPCAPH=0.0_r8
      TPCAPM=0.0_r8
      TPALPB=0.0_r8
      TPFEPB=0.0_r8
      TPCPDB=0.0_r8
      TPCPHB=0.0_r8
      TPCPMB=0.0_r8
      TCO2G=0.0_r8
      TCH4G=0.0_r8
      TCOZS=0.0_r8
      TCHFS=0.0_r8
      TOXYG=0.0_r8
      TOXYS=0.0_r8
      TZ2GG=0.0_r8
      TZ2GS=0.0_r8
      TZ2OG=0.0_r8
      TZ2OS=0.0_r8
      TZNH3G=0.0_r8
      TH2GG=0.0_r8
      TH2GS=0.0_r8

      TOMC(:,:,:,:)=0.0_r8
      TOMN(:,:,:,:)=0.0_r8
      TOMP(:,:,:,:)=0.0_r8

      TOMCff(:,:,:)=0.0_r8
      TOMNff(:,:,:)=0.0_r8
      TOMPff(:,:,:)=0.0_r8

      TORC(:,:)=0.0_r8
      TORN(:,:)=0.0_r8
      TORP(:,:)=0.0_r8

      TOQC(:)=0.0_r8
      TOQN(:)=0.0_r8
      TOQP(:)=0.0_r8
      TOQA(:)=0.0_r8
      TOHC(:)=0.0_r8
      TOHN(:)=0.0_r8
      TOHP(:)=0.0_r8
      TOHA(:)=0.0_r8
      TOSC(:,:)=0.0_r8
      TOSA(:,:)=0.0_r8
      TOSN(:,:)=0.0_r8
      TOSP(:,:)=0.0_r8
      TZNFN2=0.0_r8
      TZNFNI=0.0_r8
      ZNHUX0=0.0_r8
      ZNHUXI=0.0_r8
      ZNFNX0=0.0_r8
!
!     ACCUMULATE STATE VARIABLES IN SURFACE RESIDUE FOR ADDITION
!     TO SOIL IN TILLAGE MIXING ZONE
!
!     IF(ORGC(0,NY,NX).GT.ZEROS(NY,NX))THEN
!     XCORP0=AMAX1(XCORP(NY,NX),AMIN1(1.0,
!    2(VHCPRX(NY,NX)/cpo)/ORGC(0,NY,NX)))
      XCORP0=AMAX1(0.001,XCORP(NY,NX))
!     ELSE
!     XCORP0=1.0
!     ENDIF
      CORP0=1.0-XCORP0
      DC=0.0_r8
      DN=0.0_r8
      DP=0.0_r8

      DO 3950 K=0,jcplx1
      IF(K.NE.3.AND.K.NE.4)THEN
      DO 3945 N=1,7
      DO  M=1,3
      DO NGL=1,JG
      TOMGC(M,NGL,N,K)=OMC(M,NGL,N,K,0,NY,NX)*CORP0
      TOMGN(M,NGL,N,K)=OMN(M,NGL,N,K,0,NY,NX)*CORP0
      TOMGP(M,NGL,N,K)=OMP(M,NGL,N,K,0,NY,NX)*CORP0
      OMC(M,NGL,N,K,0,NY,NX)=OMC(M,NGL,N,K,0,NY,NX)*XCORP0
      OMN(M,NGL,N,K,0,NY,NX)=OMN(M,NGL,N,K,0,NY,NX)*XCORP0
      OMP(M,NGL,N,K,0,NY,NX)=OMP(M,NGL,N,K,0,NY,NX)*XCORP0
      DC=DC+OMC(M,NGL,N,K,0,NY,NX)
      DN=DN+OMN(M,NGL,N,K,0,NY,NX)
      DP=DP+OMP(M,NGL,N,K,0,NY,NX)
      enddo
      enddo
3945  CONTINUE
      ENDIF
3950  CONTINUE

      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      TOMGCff(M,NGL,N)=OMCff(M,NGL,N,0,NY,NX)*CORP0
      TOMGNff(M,NGL,N)=OMNff(M,NGL,N,0,NY,NX)*CORP0
      TOMGPff(M,NGL,N)=OMPff(M,NGL,N,0,NY,NX)*CORP0
      OMCff(M,NGL,N,0,NY,NX)=OMCff(M,NGL,N,0,NY,NX)*XCORP0
      OMNff(M,NGL,N,0,NY,NX)=OMNff(M,NGL,N,0,NY,NX)*XCORP0
      OMPff(M,NGL,N,0,NY,NX)=OMPff(M,NGL,N,0,NY,NX)*XCORP0
      DC=DC+OMCff(M,NGL,N,0,NY,NX)
      DN=DN+OMNff(M,NGL,N,0,NY,NX)
      DP=DP+OMPff(M,NGL,N,0,NY,NX)
      enddo
      enddo
    ENDDO


      DO 3940 K=0,2
      DO 3935 M=1,2
      TORXC(M,K)=ORC(M,K,0,NY,NX)*CORP0
      TORXN(M,K)=ORN(M,K,0,NY,NX)*CORP0
      TORXP(M,K)=ORP(M,K,0,NY,NX)*CORP0
      ORC(M,K,0,NY,NX)=ORC(M,K,0,NY,NX)*XCORP0
      ORN(M,K,0,NY,NX)=ORN(M,K,0,NY,NX)*XCORP0
      ORP(M,K,0,NY,NX)=ORP(M,K,0,NY,NX)*XCORP0
      DC=DC+ORC(M,K,0,NY,NX)
      DN=DN+ORN(M,K,0,NY,NX)
      DP=DP+ORP(M,K,0,NY,NX)
3935  CONTINUE
      TOQGC(K)=OQC(K,0,NY,NX)*CORP0
      TOQGN(K)=OQN(K,0,NY,NX)*CORP0
      TOQGP(K)=OQP(K,0,NY,NX)*CORP0
      TOQGA(K)=OQA(K,0,NY,NX)*CORP0
      TOQHC(K)=OQCH(K,0,NY,NX)*CORP0
      TOQHN(K)=OQNH(K,0,NY,NX)*CORP0
      TOQHP(K)=OQPH(K,0,NY,NX)*CORP0
      TOQHA(K)=OQAH(K,0,NY,NX)*CORP0
      TOHGC(K)=OHC(K,0,NY,NX)*CORP0
      TOHGN(K)=OHN(K,0,NY,NX)*CORP0
      TOHGP(K)=OHP(K,0,NY,NX)*CORP0
      TOHGA(K)=OHA(K,0,NY,NX)*CORP0
!
!     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
!
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)*XCORP0
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)*XCORP0
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)*XCORP0
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)*XCORP0
      OQCH(K,0,NY,NX)=OQCH(K,0,NY,NX)*XCORP0
      OQNH(K,0,NY,NX)=OQNH(K,0,NY,NX)*XCORP0
      OQPH(K,0,NY,NX)=OQPH(K,0,NY,NX)*XCORP0
      OQAH(K,0,NY,NX)=OQAH(K,0,NY,NX)*XCORP0
      OHC(K,0,NY,NX)=OHC(K,0,NY,NX)*XCORP0
      OHN(K,0,NY,NX)=OHN(K,0,NY,NX)*XCORP0
      OHP(K,0,NY,NX)=OHP(K,0,NY,NX)*XCORP0
      OHA(K,0,NY,NX)=OHA(K,0,NY,NX)*XCORP0
      DC=DC+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX)+OHC(K,0,NY,NX) &
      +OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
      DN=DN+OQN(K,0,NY,NX)+OQNH(K,0,NY,NX)+OHN(K,0,NY,NX)
      DP=DP+OQP(K,0,NY,NX)+OQPH(K,0,NY,NX)+OHP(K,0,NY,NX)
      DO 3965 M=1,4
      TOSGC(M,K)=OSC(M,K,0,NY,NX)*CORP0
      TOSGA(M,K)=OSA(M,K,0,NY,NX)*CORP0
      TOSGN(M,K)=OSN(M,K,0,NY,NX)*CORP0
      TOSGP(M,K)=OSP(M,K,0,NY,NX)*CORP0
      OSC(M,K,0,NY,NX)=OSC(M,K,0,NY,NX)*XCORP0
      OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)*XCORP0
      OSN(M,K,0,NY,NX)=OSN(M,K,0,NY,NX)*XCORP0
      OSP(M,K,0,NY,NX)=OSP(M,K,0,NY,NX)*XCORP0
      DC=DC+OSC(M,K,0,NY,NX)
      DN=DN+OSN(M,K,0,NY,NX)
      DP=DP+OSP(M,K,0,NY,NX)
3965  CONTINUE
3940  CONTINUE
      TCO2GS=CO2S(0,NY,NX)*CORP0
      TCH4GS=CH4S(0,NY,NX)*CORP0
      TOXYGS=OXYS(0,NY,NX)*CORP0
      TZ2GSG=Z2GS(0,NY,NX)*CORP0
      TZ2OGS=Z2OS(0,NY,NX)*CORP0
      TH2GGS=H2GS(0,NY,NX)*CORP0
      TNH4GS=ZNH4S(0,NY,NX)*CORP0
      TNH3GS=ZNH3S(0,NY,NX)*CORP0
      TNO3GS=ZNO3S(0,NY,NX)*CORP0
      TNO2GS=ZNO2S(0,NY,NX)*CORP0
      TP14GS=H1PO4(0,NY,NX)*CORP0
      TPO4GS=H2PO4(0,NY,NX)*CORP0
      TXN4G=XN4(0,NY,NX)*CORP0
      TXOH0G=XOH0(0,NY,NX)*CORP0
      TXOH1G=XOH1(0,NY,NX)*CORP0
      TXOH2G=XOH2(0,NY,NX)*CORP0
      TXH1PG=XH1P(0,NY,NX)*CORP0
      TXH2PG=XH2P(0,NY,NX)*CORP0
      TALPOG=PALPO(0,NY,NX)*CORP0
      TFEPOG=PFEPO(0,NY,NX)*CORP0
      TCAPDG=PCAPD(0,NY,NX)*CORP0
      TCAPHG=PCAPH(0,NY,NX)*CORP0
      TCAPMG=PCAPM(0,NY,NX)*CORP0
      TNH4FG=ZNH4FA(0,NY,NX)*CORP0
      TNH3FG=ZNH3FA(0,NY,NX)*CORP0
      TNHUFG=ZNHUFA(0,NY,NX)*CORP0
      TNO3FG=ZNO3FA(0,NY,NX)*CORP0
      TZNFNG=ZNFNI(0,NY,NX)*CORP0
      TVOLWR=VOLW(0,NY,NX)*CORP0
      HFLXD=cpo*ORGC(0,NY,NX)*CORP0*TKS(0,NY,NX)
      HEATIN=HEATIN-HFLXD
      HEATSO=HEATSO-HFLXD
      TENGYR=cpw*TVOLWR*TKS(0,NY,NX)
      ORGC(0,NY,NX)=DC
      ORGN(0,NY,NX)=DN
      ORGR(0,NY,NX)=DC
      CO2S(0,NY,NX)=CO2S(0,NY,NX)*XCORP0
      CH4S(0,NY,NX)=CH4S(0,NY,NX)*XCORP0
      OXYS(0,NY,NX)=OXYS(0,NY,NX)*XCORP0
      Z2GS(0,NY,NX)=Z2GS(0,NY,NX)*XCORP0
      Z2OS(0,NY,NX)=Z2OS(0,NY,NX)*XCORP0
      H2GS(0,NY,NX)=H2GS(0,NY,NX)*XCORP0
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)*XCORP0
      ZNH3S(0,NY,NX)=ZNH3S(0,NY,NX)*XCORP0
      ZNO3S(0,NY,NX)=ZNO3S(0,NY,NX)*XCORP0
      ZNO2S(0,NY,NX)=ZNO2S(0,NY,NX)*XCORP0
      H1PO4(0,NY,NX)=H1PO4(0,NY,NX)*XCORP0
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)*XCORP0
      XN4(0,NY,NX)=XN4(0,NY,NX)*XCORP0
      XOH0(0,NY,NX)=XOH0(0,NY,NX)*XCORP0
      XOH1(0,NY,NX)=XOH1(0,NY,NX)*XCORP0
      XOH2(0,NY,NX)=XOH2(0,NY,NX)*XCORP0
      XH1P(0,NY,NX)=XH1P(0,NY,NX)*XCORP0
      XH2P(0,NY,NX)=XH2P(0,NY,NX)*XCORP0
      PALPO(0,NY,NX)=PALPO(0,NY,NX)*XCORP0
      PFEPO(0,NY,NX)=PFEPO(0,NY,NX)*XCORP0
      PCAPD(0,NY,NX)=PCAPD(0,NY,NX)*XCORP0
      PCAPH(0,NY,NX)=PCAPH(0,NY,NX)*XCORP0
      PCAPM(0,NY,NX)=PCAPM(0,NY,NX)*XCORP0
      ZNH4FA(0,NY,NX)=ZNH4FA(0,NY,NX)*XCORP0
      ZNH3FA(0,NY,NX)=ZNH3FA(0,NY,NX)*XCORP0
      ZNHUFA(0,NY,NX)=ZNHUFA(0,NY,NX)*XCORP0
      ZNO3FA(0,NY,NX)=ZNO3FA(0,NY,NX)*XCORP0
      VOLW(0,NY,NX)=VOLW(0,NY,NX)*XCORP0
      VHCP(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX) &
      +cpi*VOLI(0,NY,NX)
      VOLR(NY,NX)=VOLR(NY,NX)*XCORP0
      VOLT(0,NY,NX)=VOLT(0,NY,NX)*XCORP0
      ZNHUX0=AMAX1(ZNHUX0,ZNHU0(0,NY,NX))
      ZNHUXI=AMAX1(ZNHUXI,ZNHUI(0,NY,NX))
      ZNFNX0=AMAX1(ZNFNX0,ZNFN0(0,NY,NX))
      LL=NU(NY,NX)
!
!     REDISTRIBUTE SOIL STATE VARIABLES DURING TILLAGE
!
      DCORPZ=AMIN1(DCORP(I,NY,NX),CDPTHZ(NL(NY,NX),NY,NX))
!
!     ACCUMULATE SOIL STATE VARIABLES WITHIN TILLAGE MIXING ZONE
!
      DO 1000 L=NU(NY,NX),NL(NY,NX)
      IF(CDPTHZ(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ &
      .AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
      TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CDPTHZ(L,NY,NX) &
      -DLYR(3,L,NY,NX)))
      FI=TL/DCORPZ
      TI=TL/DLYR(3,L,NY,NX)
      TBKDX=TBKDX+FI*BKDSI(L,NY,NX)
      TFC=TFC+FI*FC(L,NY,NX)
      TWP=TWP+FI*WP(L,NY,NX)
      TSCNV=TSCNV+FI*SCNV(L,NY,NX)
      TSCNH=TSCNH+FI*SCNH(L,NY,NX)
      TSAND=TSAND+TI*SAND(L,NY,NX)
      TSILT=TSILT+TI*SILT(L,NY,NX)
      TCLAY=TCLAY+TI*CLAY(L,NY,NX)
      TXCEC=TXCEC+TI*XCEC(L,NY,NX)
      TXAEC=TXAEC+TI*XAEC(L,NY,NX)
      TGKC4=TGKC4+FI*GKC4(L,NY,NX)
      TGKCA=TGKCA+FI*GKCA(L,NY,NX)
      TGKCM=TGKCM+FI*GKCM(L,NY,NX)
      TGKCN=TGKCN+FI*GKCN(L,NY,NX)
      TGKCK=TGKCK+FI*GKCK(L,NY,NX)
      TVOLW=TVOLW+TI*VOLW(L,NY,NX)
      TVOLI=TVOLI+TI*VOLI(L,NY,NX)
!     TVOLP=TVOLP+TI*VOLP(L,NY,NX)
!     TVOLA=TVOLA+TI*VOLA(L,NY,NX)
      TENGY=TENGY+TI*(cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX)) &
      +cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX)))*TKS(L,NY,NX)
      TNH4FA=TNH4FA+TI*ZNH4FA(L,NY,NX)
      TNH3FA=TNH3FA+TI*ZNH3FA(L,NY,NX)
      TNHUFA=TNHUFA+TI*ZNHUFA(L,NY,NX)
      TNO3FA=TNO3FA+TI*ZNO3FA(L,NY,NX)
      TNH4FB=TNH4FB+TI*ZNH4FB(L,NY,NX)
      TNH3FB=TNH3FB+TI*ZNH3FB(L,NY,NX)
      TNHUFB=TNHUFB+TI*ZNHUFB(L,NY,NX)
      TNO3FB=TNO3FB+TI*ZNO3FB(L,NY,NX)
      TNH4S=TNH4S+TI*ZNH4S(L,NY,NX)
      TNH4B=TNH4B+TI*ZNH4B(L,NY,NX)
      TNH3S=TNH3S+TI*ZNH3S(L,NY,NX)
      TNH3B=TNH3B+TI*ZNH3B(L,NY,NX)
      TNO3S=TNO3S+TI*ZNO3S(L,NY,NX)
      TNO3B=TNO3B+TI*ZNO3B(L,NY,NX)
      TNO2S=TNO2S+TI*ZNO2S(L,NY,NX)
      TNO2B=TNO2B+TI*ZNO2B(L,NY,NX)
      TZAL=TZAL+TI*ZAL(L,NY,NX)
      TZFE=TZFE+TI*ZFE(L,NY,NX)
      TZHY=TZHY+TI*ZHY(L,NY,NX)
      TZCA=TZCA+TI*ZCA(L,NY,NX)
      TZMG=TZMG+TI*ZMG(L,NY,NX)
      TZNA=TZNA+TI*ZNA(L,NY,NX)
      TZKA=TZKA+TI*ZKA(L,NY,NX)
      TZOH=TZOH+TI*ZOH(L,NY,NX)
      TZSO4=TZSO4+TI*ZSO4(L,NY,NX)
      TZCL=TZCL+TI*ZCL(L,NY,NX)
      TZCO3=TZCO3+TI*ZCO3(L,NY,NX)
      TZHCO3=TZHCO3+TI*ZHCO3(L,NY,NX)
      TZALO1=TZALO1+TI*ZALOH1(L,NY,NX)
      TZALO2=TZALO2+TI*ZALOH2(L,NY,NX)
      TZALO3=TZALO3+TI*ZALOH3(L,NY,NX)
      TZALO4=TZALO4+TI*ZALOH4(L,NY,NX)
      TZALS=TZALS+TI*ZALS(L,NY,NX)
      TZFEO1=TZFEO1+TI*ZFEOH1(L,NY,NX)
      TZFEO2=TZFEO2+TI*ZFEOH2(L,NY,NX)
      TZFEO3=TZFEO3+TI*ZFEOH3(L,NY,NX)
      TZFEO4=TZFEO4+TI*ZFEOH4(L,NY,NX)
      TZFES=TZFES+TI*ZFES(L,NY,NX)
      TZCAO=TZCAO+TI*ZCAO(L,NY,NX)
      TZCAC=TZCAC+TI*ZCAC(L,NY,NX)
      TZCAH=TZCAH+TI*ZCAH(L,NY,NX)
      TZCAS=TZCAS+TI*ZCAS(L,NY,NX)
      TZMGO=TZMGO+TI*ZMGO(L,NY,NX)
      TZMGC=TZMGC+TI*ZMGC(L,NY,NX)
      TZMGH=TZMGH+TI*ZMGH(L,NY,NX)
      TZMGS=TZMGS+TI*ZMGS(L,NY,NX)
      TZNAC=TZNAC+TI*ZNAC(L,NY,NX)
      TZNAS=TZNAS+TI*ZNAS(L,NY,NX)
      TZKAS=TZKAS+TI*ZKAS(L,NY,NX)
      TH0PO4=TH0PO4+TI*H0PO4(L,NY,NX)
      TH1PO4=TH1PO4+TI*H1PO4(L,NY,NX)
      TH2PO4=TH2PO4+TI*H2PO4(L,NY,NX)
      TH3PO4=TH3PO4+TI*H3PO4(L,NY,NX)
      TZFE1P=TZFE1P+TI*ZFE1P(L,NY,NX)
      TZFE2P=TZFE2P+TI*ZFE2P(L,NY,NX)
      TZCA0P=TZCA0P+TI*ZCA0P(L,NY,NX)
      TZCA1P=TZCA1P+TI*ZCA1P(L,NY,NX)
      TZCA2P=TZCA2P+TI*ZCA2P(L,NY,NX)
      TZMG1P=TZMG1P+TI*ZMG1P(L,NY,NX)
      TH0POB=TH0POB+TI*H0POB(L,NY,NX)
      TH1POB=TH1POB+TI*H1POB(L,NY,NX)
      TH2POB=TH2POB+TI*H2POB(L,NY,NX)
      TH3POB=TH3POB+TI*H3POB(L,NY,NX)
      TFE1PB=TFE1PB+TI*ZFE1PB(L,NY,NX)
      TFE2PB=TFE2PB+TI*ZFE2PB(L,NY,NX)
      TCA0PB=TCA0PB+TI*ZCA0PB(L,NY,NX)
      TCA1PB=TCA1PB+TI*ZCA1PB(L,NY,NX)
      TCA2PB=TCA2PB+TI*ZCA2PB(L,NY,NX)
      TMG1PB=TMG1PB+TI*ZMG1PB(L,NY,NX)
      TXNH4=TXNH4+TI*XN4(L,NY,NX)
      TXNHB=TXNHB+TI*XNB(L,NY,NX)
      TXHY=TXHY+TI*XHY(L,NY,NX)
      TXAL=TXAL+TI*XAL(L,NY,NX)
      TXFE=TXFE+TI*XFE(L,NY,NX)
      TXCA=TXCA+TI*XCA(L,NY,NX)
      TXMG=TXMG+TI*XMG(L,NY,NX)
      TXNA=TXNA+TI*XNA(L,NY,NX)
      TXKA=TXKA+TI*XKA(L,NY,NX)
      TXHC=TXHC+TI*XHC(L,NY,NX)
      TXAL2=TXAL2+TI*XALO2(L,NY,NX)
      TXFE2=TXFE2+TI*XFEO2(L,NY,NX)
      TXOH0=TXOH0+TI*XOH0(L,NY,NX)
      TXOH1=TXOH1+TI*XOH1(L,NY,NX)
      TXOH2=TXOH2+TI*XOH2(L,NY,NX)
      TXH1P=TXH1P+TI*XH1P(L,NY,NX)
      TXH2P=TXH2P+TI*XH2P(L,NY,NX)
      TXOH0B=TXOH0B+TI*XOH0B(L,NY,NX)
      TXOH1B=TXOH1B+TI*XOH1B(L,NY,NX)
      TXOH2B=TXOH2B+TI*XOH2B(L,NY,NX)
      TXH1PB=TXH1PB+TI*XH1PB(L,NY,NX)
      TXH2PB=TXH2PB+TI*XH2PB(L,NY,NX)
      TPALOH=TPALOH+TI*PALOH(L,NY,NX)
      TPFEOH=TPFEOH+TI*PFEOH(L,NY,NX)
      TPCACO=TPCACO+TI*PCACO(L,NY,NX)
      TPCASO=TPCASO+TI*PCASO(L,NY,NX)
      TPALPO=TPALPO+TI*PALPO(L,NY,NX)
      TPFEPO=TPFEPO+TI*PFEPO(L,NY,NX)
      TPCAPD=TPCAPD+TI*PCAPD(L,NY,NX)
      TPCAPH=TPCAPH+TI*PCAPH(L,NY,NX)
      TPCAPM=TPCAPM+TI*PCAPM(L,NY,NX)
      TPALPB=TPALPB+TI*PALPB(L,NY,NX)
      TPFEPB=TPFEPB+TI*PFEPB(L,NY,NX)
      TPCPDB=TPCPDB+TI*PCPDB(L,NY,NX)
      TPCPHB=TPCPHB+TI*PCPHB(L,NY,NX)
      TPCPMB=TPCPMB+TI*PCPMB(L,NY,NX)
      TCO2G=TCO2G+TI*CO2G(L,NY,NX)
      TCH4G=TCH4G+TI*CH4G(L,NY,NX)
      TCOZS=TCOZS+TI*CO2S(L,NY,NX)
      TCHFS=TCHFS+TI*CH4S(L,NY,NX)
      TOXYG=TOXYG+TI*OXYG(L,NY,NX)
      TOXYS=TOXYS+TI*OXYS(L,NY,NX)
      TZ2GG=TZ2GG+TI*Z2GG(L,NY,NX)
      TZ2GS=TZ2GS+TI*Z2GS(L,NY,NX)
      TZ2OG=TZ2OG+TI*Z2OG(L,NY,NX)
      TZ2OS=TZ2OS+TI*Z2OS(L,NY,NX)
      TZNH3G=TZNH3G+TI*ZNH3G(L,NY,NX)
      TH2GG=TH2GG+TI*H2GG(L,NY,NX)
      TH2GS=TH2GS+TI*H2GS(L,NY,NX)
      DO 4985 K=0,jcplx1
      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      TOMC(M,NGL,N,K)=TOMC(M,NGL,N,K)+TI*OMC(M,NGL,N,K,L,NY,NX)
      TOMN(M,NGL,N,K)=TOMN(M,NGL,N,K)+TI*OMN(M,NGL,N,K,L,NY,NX)
      TOMP(M,NGL,N,K)=TOMP(M,NGL,N,K)+TI*OMP(M,NGL,N,K,L,NY,NX)
      enddo
      enddo
      enddo
4985  CONTINUE

      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      TOMCff(M,NGL,N)=TOMCff(M,NGL,N)+TI*OMCff(M,NGL,N,L,NY,NX)
      TOMNff(M,NGL,N)=TOMNff(M,NGL,N)+TI*OMNff(M,NGL,N,L,NY,NX)
      TOMPff(M,NGL,N)=TOMPff(M,NGL,N)+TI*OMPff(M,NGL,N,L,NY,NX)
      enddo
      enddo
      enddo

      DO 4980 K=0,jcplx1
      DO 4975 M=1,2
      TORC(M,K)=TORC(M,K)+TI*ORC(M,K,L,NY,NX)
      TORN(M,K)=TORN(M,K)+TI*ORN(M,K,L,NY,NX)
      TORP(M,K)=TORP(M,K)+TI*ORP(M,K,L,NY,NX)
4975  CONTINUE
      TOQC(K)=TOQC(K)+TI*OQC(K,L,NY,NX)
      TOQN(K)=TOQN(K)+TI*OQN(K,L,NY,NX)
      TOQP(K)=TOQP(K)+TI*OQP(K,L,NY,NX)
      TOQA(K)=TOQA(K)+TI*OQA(K,L,NY,NX)
      TOHC(K)=TOHC(K)+TI*OHC(K,L,NY,NX)
      TOHN(K)=TOHN(K)+TI*OHN(K,L,NY,NX)
      TOHP(K)=TOHP(K)+TI*OHP(K,L,NY,NX)
      TOHA(K)=TOHA(K)+TI*OHA(K,L,NY,NX)
      DO 4970 M=1,4
      TOSC(M,K)=TOSC(M,K)+TI*OSC(M,K,L,NY,NX)
      TOSA(M,K)=TOSA(M,K)+TI*OSA(M,K,L,NY,NX)
      TOSN(M,K)=TOSN(M,K)+TI*OSN(M,K,L,NY,NX)
      TOSP(M,K)=TOSP(M,K)+TI*OSP(M,K,L,NY,NX)
4970  CONTINUE
4980  CONTINUE
      ZNHUX0=AMAX1(ZNHUX0,ZNHU0(L,NY,NX))
      ZNHUXI=AMAX1(ZNHUXI,ZNHUI(L,NY,NX))
      ZNFNX0=AMAX1(ZNFNX0,ZNFN0(L,NY,NX))
      TZNFNI=TZNFNI+ZNFNI(L,NY,NX)
      LL=L
      ENDIF
1000  CONTINUE
!
!     CHANGE SOIL STATE VARIABLES IN TILLAGE MIXING ZONE
!     TO ACCOUNT FOR REDISTRIBUTION FROM MIXING
!
      DO 2000 L=NU(NY,NX),LL
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
      TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CDPTHZ(L,NY,NX) &
      -DLYR(3,L,NY,NX)))
      FI=TL/DCORPZ
      TI=TL/DLYR(3,L,NY,NX)
      TX=1.0-TI
      BKDSI(L,NY,NX)=TI*(BKDSI(L,NY,NX)+CORP*(TBKDX-BKDSI(L,NY,NX))) &
      +TX*BKDSI(L,NY,NX)
      FC(L,NY,NX)=TI*(FC(L,NY,NX)+CORP*(TFC-FC(L,NY,NX))) &
      +TX*FC(L,NY,NX)
      WP(L,NY,NX)=TI*(WP(L,NY,NX)+CORP*(TWP-WP(L,NY,NX))) &
      +TX*WP(L,NY,NX)
      SCNV(L,NY,NX)=TI*(SCNV(L,NY,NX)+CORP*(TSCNV-SCNV(L,NY,NX))) &
      +TX*SCNV(L,NY,NX)
      SCNH(L,NY,NX)=TI*(SCNH(L,NY,NX)+CORP*(TSCNH-SCNH(L,NY,NX))) &
      +TX*SCNH(L,NY,NX)
      SAND(L,NY,NX)=TI*SAND(L,NY,NX)+CORP*(FI*TSAND-TI*SAND(L,NY,NX)) &
      +TX*SAND(L,NY,NX)
      SILT(L,NY,NX)=TI*SILT(L,NY,NX)+CORP*(FI*TSILT-TI*SILT(L,NY,NX)) &
      +TX*SILT(L,NY,NX)
      CLAY(L,NY,NX)=TI*CLAY(L,NY,NX)+CORP*(FI*TCLAY-TI*CLAY(L,NY,NX)) &
      +TX*CLAY(L,NY,NX)
      XCEC(L,NY,NX)=TI*XCEC(L,NY,NX)+CORP*(FI*TXCEC-TI*XCEC(L,NY,NX)) &
      +TX*XCEC(L,NY,NX)
      XAEC(L,NY,NX)=TI*XAEC(L,NY,NX)+CORP*(FI*TXAEC-TI*XAEC(L,NY,NX)) &
      +TX*XAEC(L,NY,NX)
      GKC4(L,NY,NX)=TI*(GKC4(L,NY,NX)+CORP*(TGKC4-GKC4(L,NY,NX))) &
      +TX*GKC4(L,NY,NX)
      GKCA(L,NY,NX)=TI*(GKCA(L,NY,NX)+CORP*(TGKCA-GKCA(L,NY,NX))) &
      +TX*GKCA(L,NY,NX)
      GKCM(L,NY,NX)=TI*(GKCM(L,NY,NX)+CORP*(TGKCM-GKCM(L,NY,NX))) &
      +TX*GKCM(L,NY,NX)
      GKCN(L,NY,NX)=TI*(GKCN(L,NY,NX)+CORP*(TGKCN-GKCN(L,NY,NX))) &
      +TX*GKCN(L,NY,NX)
      GKCK(L,NY,NX)=TI*(GKCK(L,NY,NX)+CORP*(TGKCK-GKCK(L,NY,NX))) &
      +TX*GKCK(L,NY,NX)
      ENGYM=VHCM(L,NY,NX)*TKS(L,NY,NX)
      ENGYV=(cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX)) &
      +cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX)))*TKS(L,NY,NX)
      VOLW(L,NY,NX)=TI*VOLW(L,NY,NX)+CORP*(FI*TVOLW-TI*VOLW(L,NY,NX)) &
      +TX*VOLW(L,NY,NX)+FI*TVOLWR
      VOLI(L,NY,NX)=TI*VOLI(L,NY,NX)+CORP*(FI*TVOLI-TI*VOLI(L,NY,NX)) &
      +TX*VOLI(L,NY,NX)
      VOLWX(L,NY,NX)=VOLW(L,NY,NX)
!     VOLW(L,NY,NX)=VOLW(L,NY,NX)+CORP*VOLWH(L,NY,NX)
!     VOLI(L,NY,NX)=VOLI(L,NY,NX)+CORP*VOLIH(L,NY,NX)
!     VOLA(L,NY,NX)=VOLA(L,NY,NX)+CORP*VOLAH(L,NY,NX)
!     VOLWH(L,NY,NX)=XCORP(NY,NX)*VOLWH(L,NY,NX)
!     VOLIH(L,NY,NX)=XCORP(NY,NX)*VOLIH(L,NY,NX)
!     VOLAH(L,NY,NX)=XCORP(NY,NX)*VOLAH(L,NY,NX)
!     FHOL(L,NY,NX)=XCORP(NY,NX)*FHOL(L,NY,NX)
      ENGYL=TI*ENGYV+CORP*(FI*TENGY-TI*ENGYV)+TX*ENGYV+FI*TENGYR
      VHCP(L,NY,NX)=VHCM(L,NY,NX)+cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX)) &
      +cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
      TKS(L,NY,NX)=(ENGYM+ENGYL)/VHCP(L,NY,NX)
      TCS(L,NY,NX)=TKS(L,NY,NX)-TC2K
      ZNH4FA(L,NY,NX)=TI*ZNH4FA(L,NY,NX)+CORP*(FI*TNH4FA &
      -TI*ZNH4FA(L,NY,NX))+TX*ZNH4FA(L,NY,NX)
      ZNH3FA(L,NY,NX)=TI*ZNH3FA(L,NY,NX)+CORP*(FI*TNH3FA &
      -TI*ZNH3FA(L,NY,NX))+TX*ZNH3FA(L,NY,NX)
      ZNHUFA(L,NY,NX)=TI*ZNHUFA(L,NY,NX)+CORP*(FI*TNHUFA &
      -TI*ZNHUFA(L,NY,NX))+TX*ZNHUFA(L,NY,NX)
      ZNO3FA(L,NY,NX)=TI*ZNO3FA(L,NY,NX)+CORP*(FI*TNO3FA &
      -TI*ZNO3FA(L,NY,NX))+TX*ZNO3FA(L,NY,NX)
      ZNH4FB(L,NY,NX)=TI*ZNH4FB(L,NY,NX)+CORP*(FI*TNH4FB &
      -TI*ZNH4FB(L,NY,NX))+TX*ZNH4FB(L,NY,NX)
      ZNH3FB(L,NY,NX)=TI*ZNH3FB(L,NY,NX)+CORP*(FI*TNH3FB &
      -TI*ZNH3FB(L,NY,NX))+TX*ZNH3FB(L,NY,NX)
      ZNHUFB(L,NY,NX)=TI*ZNHUFB(L,NY,NX)+CORP*(FI*TNHUFB &
      -TI*ZNHUFB(L,NY,NX))+TX*ZNHUFB(L,NY,NX)
      ZNO3FB(L,NY,NX)=TI*ZNO3FB(L,NY,NX)+CORP*(FI*TNO3FB &
      -TI*ZNO3FB(L,NY,NX))+TX*ZNO3FB(L,NY,NX)
      ZNH4S(L,NY,NX)=TI*ZNH4S(L,NY,NX)+CORP*(FI*TNH4S-TI*ZNH4S(L,NY,NX)) &
      +TX*ZNH4S(L,NY,NX)+CORP*ZNH4SH(L,NY,NX)
      ZNH4B(L,NY,NX)=TI*ZNH4B(L,NY,NX)+CORP*(FI*TNH4B-TI*ZNH4B(L,NY,NX)) &
      +TX*ZNH4B(L,NY,NX)+CORP*ZNH4BH(L,NY,NX)
      ZNH3S(L,NY,NX)=TI*ZNH3S(L,NY,NX)+CORP*(FI*TNH3S-TI*ZNH3S(L,NY,NX)) &
      +TX*ZNH3S(L,NY,NX)+CORP*ZNH3SH(L,NY,NX)
      ZNH3B(L,NY,NX)=TI*ZNH3B(L,NY,NX)+CORP*(FI*TNH3B-TI*ZNH3B(L,NY,NX)) &
      +TX*ZNH3B(L,NY,NX)+CORP*ZNH3BH(L,NY,NX)
      ZNO3S(L,NY,NX)=TI*ZNO3S(L,NY,NX)+CORP*(FI*TNO3S-TI*ZNO3S(L,NY,NX)) &
      +TX*ZNO3S(L,NY,NX)+CORP*ZNO3SH(L,NY,NX)
      ZNO3B(L,NY,NX)=TI*ZNO3B(L,NY,NX)+CORP*(FI*TNO3B-TI*ZNO3B(L,NY,NX)) &
      +TX*ZNO3B(L,NY,NX)+CORP*ZNO3BH(L,NY,NX)
      ZNO2S(L,NY,NX)=TI*ZNO2S(L,NY,NX)+CORP*(FI*TNO2S-TI*ZNO2S(L,NY,NX)) &
      +TX*ZNO2S(L,NY,NX)+CORP*ZNO2SH(L,NY,NX)
      ZNO2B(L,NY,NX)=TI*ZNO2B(L,NY,NX)+CORP*(FI*TNO2B-TI*ZNO2B(L,NY,NX)) &
      +TX*ZNO2B(L,NY,NX)+CORP*ZNO2BH(L,NY,NX)
      ZAL(L,NY,NX)=TI*ZAL(L,NY,NX)+CORP*(FI*TZAL-TI*ZAL(L,NY,NX)) &
      +TX*ZAL(L,NY,NX)+CORP*ZALH(L,NY,NX)
      ZFE(L,NY,NX)=TI*ZFE(L,NY,NX)+CORP*(FI*TZFE-TI*ZFE(L,NY,NX)) &
      +TX*ZFE(L,NY,NX)+CORP*ZFEH(L,NY,NX)
      ZHY(L,NY,NX)=TI*ZHY(L,NY,NX)+CORP*(FI*TZHY-TI*ZHY(L,NY,NX)) &
      +TX*ZHY(L,NY,NX)+CORP*ZHYH(L,NY,NX)
      ZCA(L,NY,NX)=TI*ZCA(L,NY,NX)+CORP*(FI*TZCA-TI*ZCA(L,NY,NX)) &
      +TX*ZCA(L,NY,NX)+CORP*ZCCH(L,NY,NX)
      ZMG(L,NY,NX)=TI*ZMG(L,NY,NX)+CORP*(FI*TZMG-TI*ZMG(L,NY,NX)) &
      +TX*ZMG(L,NY,NX)+CORP*ZMAH(L,NY,NX)
      ZNA(L,NY,NX)=TI*ZNA(L,NY,NX)+CORP*(FI*TZNA-TI*ZNA(L,NY,NX)) &
      +TX*ZNA(L,NY,NX)+CORP*ZNAH(L,NY,NX)
      ZKA(L,NY,NX)=TI*ZKA(L,NY,NX)+CORP*(FI*TZKA-TI*ZKA(L,NY,NX)) &
      +TX*ZKA(L,NY,NX)+CORP*ZKAH(L,NY,NX)
      ZOH(L,NY,NX)=TI*ZOH(L,NY,NX)+CORP*(FI*TZOH-TI*ZOH(L,NY,NX)) &
      +TX*ZOH(L,NY,NX)+CORP*ZOHH(L,NY,NX)
      ZSO4(L,NY,NX)=TI*ZSO4(L,NY,NX)+CORP*(FI*TZSO4-TI*ZSO4(L,NY,NX)) &
      +TX*ZSO4(L,NY,NX)+CORP*ZSO4H(L,NY,NX)
      ZCL(L,NY,NX)=TI*ZCL(L,NY,NX)+CORP*(FI*TZCL-TI*ZCL(L,NY,NX)) &
      +TX*ZCL(L,NY,NX)+CORP*ZCLH(L,NY,NX)
      ZCO3(L,NY,NX)=TI*ZCO3(L,NY,NX)+CORP*(FI*TZCO3-TI*ZCO3(L,NY,NX)) &
      +TX*ZCO3(L,NY,NX)+CORP*ZCO3H(L,NY,NX)
      ZHCO3(L,NY,NX)=TI*ZHCO3(L,NY,NX)+CORP*(FI*TZHCO3 &
      -TI*ZHCO3(L,NY,NX))+TX*ZHCO3(L,NY,NX)+CORP*ZHCO3H(L,NY,NX)
      ZALOH1(L,NY,NX)=TI*ZALOH1(L,NY,NX)+CORP*(FI*TZALO1 &
      -TI*ZALOH1(L,NY,NX))+TX*ZALOH1(L,NY,NX)+CORP*ZALO1H(L,NY,NX)
      ZALOH2(L,NY,NX)=TI*ZALOH2(L,NY,NX)+CORP*(FI*TZALO2 &
      -TI*ZALOH2(L,NY,NX))+TX*ZALOH2(L,NY,NX)+CORP*ZALO2H(L,NY,NX)
      ZALOH3(L,NY,NX)=TI*ZALOH3(L,NY,NX)+CORP*(FI*TZALO3 &
      -TI*ZALOH3(L,NY,NX))+TX*ZALOH3(L,NY,NX)+CORP*ZALO3H(L,NY,NX)
      ZALOH4(L,NY,NX)=TI*ZALOH4(L,NY,NX)+CORP*(FI*TZALO4 &
      -TI*ZALOH4(L,NY,NX))+TX*ZALOH4(L,NY,NX)+CORP*ZALO4H(L,NY,NX)
      ZALS(L,NY,NX)=TI*ZALS(L,NY,NX)+CORP*(FI*TZALS-TI*ZALS(L,NY,NX)) &
      +TX*ZALS(L,NY,NX)+CORP*ZALSH(L,NY,NX)
      ZFEOH1(L,NY,NX)=TI*ZFEOH1(L,NY,NX)+CORP*(FI*TZFEO1 &
      -TI*ZFEOH1(L,NY,NX))+TX*ZFEOH1(L,NY,NX)+CORP*ZFEO1H(L,NY,NX)
      ZFEOH2(L,NY,NX)=TI*ZFEOH2(L,NY,NX)+CORP*(FI*TZFEO2 &
      -TI*ZFEOH2(L,NY,NX))+TX*ZFEOH2(L,NY,NX)+CORP*ZFEO2H(L,NY,NX)
      ZFEOH3(L,NY,NX)=TI*ZFEOH3(L,NY,NX)+CORP*(FI*TZFEO3 &
      -TI*ZFEOH3(L,NY,NX))+TX*ZFEOH3(L,NY,NX)+CORP*ZFEO3H(L,NY,NX)
      ZFEOH4(L,NY,NX)=TI*ZFEOH4(L,NY,NX)+CORP*(FI*TZFEO4 &
      -TI*ZFEOH4(L,NY,NX))+TX*ZFEOH4(L,NY,NX)+CORP*ZFEO4H(L,NY,NX)
      ZFES(L,NY,NX)=TI*ZFES(L,NY,NX)+CORP*(FI*TZFES-TI*ZFES(L,NY,NX)) &
      +TX*ZFES(L,NY,NX)+CORP*ZFESH(L,NY,NX)
      ZCAO(L,NY,NX)=TI*ZCAO(L,NY,NX)+CORP*(FI*TZCAO-TI*ZCAO(L,NY,NX)) &
      +TX*ZCAO(L,NY,NX)+CORP*ZCAOH(L,NY,NX)
      ZCAC(L,NY,NX)=TI*ZCAC(L,NY,NX)+CORP*(FI*TZCAC-TI*ZCAC(L,NY,NX)) &
      +TX*ZCAC(L,NY,NX)+CORP*ZCACH(L,NY,NX)
      ZCAH(L,NY,NX)=TI*ZCAH(L,NY,NX)+CORP*(FI*TZCAH-TI*ZCAH(L,NY,NX)) &
      +TX*ZCAH(L,NY,NX)+CORP*ZCAHH(L,NY,NX)
      ZCAS(L,NY,NX)=TI*ZCAS(L,NY,NX)+CORP*(FI*TZCAS-TI*ZCAS(L,NY,NX)) &
      +TX*ZCAS(L,NY,NX)+CORP*ZCASH(L,NY,NX)
      ZMGO(L,NY,NX)=TI*ZMGO(L,NY,NX)+CORP*(FI*TZMGO-TI*ZMGO(L,NY,NX)) &
      +TX*ZMGO(L,NY,NX)+CORP*ZMGOH(L,NY,NX)
      ZMGC(L,NY,NX)=TI*ZMGC(L,NY,NX)+CORP*(FI*TZMGC-TI*ZMGC(L,NY,NX)) &
      +TX*ZMGC(L,NY,NX)+CORP*ZMGCH(L,NY,NX)
      ZMGH(L,NY,NX)=TI*ZMGH(L,NY,NX)+CORP*(FI*TZMGH-TI*ZMGH(L,NY,NX)) &
      +TX*ZMGH(L,NY,NX)+CORP*ZMGHH(L,NY,NX)
      ZMGS(L,NY,NX)=TI*ZMGS(L,NY,NX)+CORP*(FI*TZMGS-TI*ZMGS(L,NY,NX)) &
      +TX*ZMGS(L,NY,NX)+CORP*ZMGSH(L,NY,NX)
      ZNAC(L,NY,NX)=TI*ZNAC(L,NY,NX)+CORP*(FI*TZNAC-TI*ZNAC(L,NY,NX)) &
      +TX*ZNAC(L,NY,NX)+CORP*ZNACH(L,NY,NX)
      ZNAS(L,NY,NX)=TI*ZNAS(L,NY,NX)+CORP*(FI*TZNAS-TI*ZNAS(L,NY,NX)) &
      +TX*ZNAS(L,NY,NX)+CORP*ZNASH(L,NY,NX)
      ZKAS(L,NY,NX)=TI*ZKAS(L,NY,NX)+CORP*(FI*TZKAS-TI*ZKAS(L,NY,NX)) &
      +TX*ZKAS(L,NY,NX)+CORP*ZKASH(L,NY,NX)
      H0PO4(L,NY,NX)=TI*H0PO4(L,NY,NX)+CORP*(FI*TH0PO4 &
      -TI*H0PO4(L,NY,NX))+TX*H0PO4(L,NY,NX)+CORP*H0PO4H(L,NY,NX)
      H1PO4(L,NY,NX)=TI*H1PO4(L,NY,NX)+CORP*(FI*TH1PO4 &
      -TI*H1PO4(L,NY,NX))+TX*H1PO4(L,NY,NX)+CORP*H1PO4H(L,NY,NX)
      H2PO4(L,NY,NX)=TI*H2PO4(L,NY,NX)+CORP*(FI*TH2PO4 &
      -TI*H2PO4(L,NY,NX))+TX*H2PO4(L,NY,NX)+CORP*H2PO4H(L,NY,NX)
      H3PO4(L,NY,NX)=TI*H3PO4(L,NY,NX)+CORP*(FI*TH3PO4 &
      -TI*H3PO4(L,NY,NX))+TX*H3PO4(L,NY,NX)+CORP*H3PO4H(L,NY,NX)
      ZFE1P(L,NY,NX)=TI*ZFE1P(L,NY,NX)+CORP*(FI*TZFE1P &
      -TI*ZFE1P(L,NY,NX))+TX*ZFE1P(L,NY,NX)+CORP*ZFE1PH(L,NY,NX)
      ZFE2P(L,NY,NX)=TI*ZFE2P(L,NY,NX)+CORP*(FI*TZFE2P &
      -TI*ZFE2P(L,NY,NX))+TX*ZFE2P(L,NY,NX)+CORP*ZFE2PH(L,NY,NX)
      ZCA0P(L,NY,NX)=TI*ZCA0P(L,NY,NX)+CORP*(FI*TZCA0P &
      -TI*ZCA0P(L,NY,NX))+TX*ZCA0P(L,NY,NX)+CORP*ZCA0PH(L,NY,NX)
      ZCA1P(L,NY,NX)=TI*ZCA1P(L,NY,NX)+CORP*(FI*TZCA1P &
      -TI*ZCA1P(L,NY,NX))+TX*ZCA1P(L,NY,NX)+CORP*ZCA1PH(L,NY,NX)
      ZCA2P(L,NY,NX)=TI*ZCA2P(L,NY,NX)+CORP*(FI*TZCA2P &
      -TI*ZCA2P(L,NY,NX))+TX*ZCA2P(L,NY,NX)+CORP*ZCA2PH(L,NY,NX)
      ZMG1P(L,NY,NX)=TI*ZMG1P(L,NY,NX)+CORP*(FI*TZMG1P &
      -TI*ZMG1P(L,NY,NX))+TX*ZMG1P(L,NY,NX)+CORP*ZMG1PH(L,NY,NX)
      H0POB(L,NY,NX)=TI*H0POB(L,NY,NX)+CORP*(FI*TH0POB &
      -TI*H0POB(L,NY,NX))+TX*H0POB(L,NY,NX)+CORP*H0POBH(L,NY,NX)
      H1POB(L,NY,NX)=TI*H1POB(L,NY,NX)+CORP*(FI*TH1POB &
      -TI*H1POB(L,NY,NX))+TX*H1POB(L,NY,NX)+CORP*H1POBH(L,NY,NX)
      H2POB(L,NY,NX)=TI*H2POB(L,NY,NX)+CORP*(FI*TH2POB &
      -TI*H2POB(L,NY,NX))+TX*H2POB(L,NY,NX)+CORP*H2POBH(L,NY,NX)
      H3POB(L,NY,NX)=TI*H3POB(L,NY,NX)+CORP*(FI*TH3POB &
      -TI*H3POB(L,NY,NX))+TX*H3POB(L,NY,NX)+CORP*H3POBH(L,NY,NX)
      ZFE1PB(L,NY,NX)=TI*ZFE1PB(L,NY,NX)+CORP*(FI*TFE1PB &
      -TI*ZFE1PB(L,NY,NX))+TX*ZFE1PB(L,NY,NX)+CORP*ZFE1BH(L,NY,NX)
      ZFE2PB(L,NY,NX)=TI*ZFE2PB(L,NY,NX)+CORP*(FI*TFE2PB &
      -TI*ZFE2PB(L,NY,NX))+TX*ZFE2PB(L,NY,NX)+CORP*ZFE2BH(L,NY,NX)
      ZCA0PB(L,NY,NX)=TI*ZCA0PB(L,NY,NX)+CORP*(FI*TCA0PB &
      -TI*ZCA0PB(L,NY,NX))+TX*ZCA0PB(L,NY,NX)+CORP*ZCA0BH(L,NY,NX)
      ZCA1PB(L,NY,NX)=TI*ZCA1PB(L,NY,NX)+CORP*(FI*TCA1PB &
      -TI*ZCA1PB(L,NY,NX))+TX*ZCA1PB(L,NY,NX)+CORP*ZCA1BH(L,NY,NX)
      ZCA2PB(L,NY,NX)=TI*ZCA2PB(L,NY,NX)+CORP*(FI*TCA2PB &
      -TI*ZCA2PB(L,NY,NX))+TX*ZCA2PB(L,NY,NX)+CORP*ZCA2BH(L,NY,NX)
      ZMG1PB(L,NY,NX)=TI*ZMG1PB(L,NY,NX)+CORP*(FI*TMG1PB &
      -TI*ZMG1PB(L,NY,NX))+TX*ZMG1PB(L,NY,NX)+CORP*ZMG1BH(L,NY,NX)
      XN4(L,NY,NX)=TI*XN4(L,NY,NX)+CORP*(FI*TXNH4-TI*XN4(L,NY,NX)) &
      +TX*XN4(L,NY,NX)
      XNB(L,NY,NX)=TI*XNB(L,NY,NX)+CORP*(FI*TXNHB-TI*XNB(L,NY,NX)) &
      +TX*XNB(L,NY,NX)
      XHY(L,NY,NX)=TI*XHY(L,NY,NX)+CORP*(FI*TXHY-TI*XHY(L,NY,NX)) &
      +TX*XHY(L,NY,NX)
      XAL(L,NY,NX)=TI*XAL(L,NY,NX)+CORP*(FI*TXAL-TI*XAL(L,NY,NX)) &
      +TX*XAL(L,NY,NX)
      XFE(L,NY,NX)=TI*XFE(L,NY,NX)+CORP*(FI*TXFE-TI*XFE(L,NY,NX)) &
      +TX*XFE(L,NY,NX)
      XCA(L,NY,NX)=TI*XCA(L,NY,NX)+CORP*(FI*TXCA-TI*XCA(L,NY,NX)) &
      +TX*XCA(L,NY,NX)
      XMG(L,NY,NX)=TI*XMG(L,NY,NX)+CORP*(FI*TXMG-TI*XMG(L,NY,NX)) &
      +TX*XMG(L,NY,NX)
      XNA(L,NY,NX)=TI*XNA(L,NY,NX)+CORP*(FI*TXNA-TI*XNA(L,NY,NX)) &
      +TX*XNA(L,NY,NX)
      XKA(L,NY,NX)=TI*XKA(L,NY,NX)+CORP*(FI*TXKA-TI*XKA(L,NY,NX)) &
      +TX*XKA(L,NY,NX)
      XHC(L,NY,NX)=TI*XHC(L,NY,NX)+CORP*(FI*TXHC-TI*XHC(L,NY,NX)) &
      +TX*XHC(L,NY,NX)
      XALO2(L,NY,NX)=TI*XALO2(L,NY,NX)+CORP*(FI*TXAL2-TI*XALO2(L,NY,NX)) &
      +TX*XALO2(L,NY,NX)
      XFEO2(L,NY,NX)=TI*XFEO2(L,NY,NX)+CORP*(FI*TXFE2-TI*XFEO2(L,NY,NX)) &
      +TX*XFEO2(L,NY,NX)
      XOH0(L,NY,NX)=TI*XOH0(L,NY,NX)+CORP*(FI*TXOH0-TI*XOH0(L,NY,NX)) &
      +TX*XOH0(L,NY,NX)
      XOH1(L,NY,NX)=TI*XOH1(L,NY,NX)+CORP*(FI*TXOH1-TI*XOH1(L,NY,NX)) &
      +TX*XOH1(L,NY,NX)
      XOH2(L,NY,NX)=TI*XOH2(L,NY,NX)+CORP*(FI*TXOH2-TI*XOH2(L,NY,NX)) &
      +TX*XOH2(L,NY,NX)
      XH1P(L,NY,NX)=TI*XH1P(L,NY,NX)+CORP*(FI*TXH1P-TI*XH1P(L,NY,NX)) &
      +TX*XH1P(L,NY,NX)
      XH2P(L,NY,NX)=TI*XH2P(L,NY,NX)+CORP*(FI*TXH2P-TI*XH2P(L,NY,NX)) &
      +TX*XH2P(L,NY,NX)
      XOH0B(L,NY,NX)=TI*XOH0B(L,NY,NX)+CORP*(FI*TXOH0B &
      -TI*XOH0B(L,NY,NX))+TX*XOH0B(L,NY,NX)
      XOH1B(L,NY,NX)=TI*XOH1B(L,NY,NX)+CORP*(FI*TXOH1B &
      -TI*XOH1B(L,NY,NX))+TX*XOH1B(L,NY,NX)
      XOH2B(L,NY,NX)=TI*XOH2B(L,NY,NX)+CORP*(FI*TXOH2B &
      -TI*XOH2B(L,NY,NX))+TX*XOH2B(L,NY,NX)
      XH1PB(L,NY,NX)=TI*XH1PB(L,NY,NX)+CORP*(FI*TXH1PB &
      -TI*XH1PB(L,NY,NX))+TX*XH1PB(L,NY,NX)
      XH2PB(L,NY,NX)=TI*XH2PB(L,NY,NX)+CORP*(FI*TXH2PB &
      -TI*XH2PB(L,NY,NX))+TX*XH2PB(L,NY,NX)
      PALOH(L,NY,NX)=TI*PALOH(L,NY,NX)+CORP*(FI*TPALOH &
      -TI*PALOH(L,NY,NX))+TX*PALOH(L,NY,NX)
      PFEOH(L,NY,NX)=TI*PFEOH(L,NY,NX)+CORP*(FI*TPFEOH &
      -TI*PFEOH(L,NY,NX))+TX*PFEOH(L,NY,NX)
      PCACO(L,NY,NX)=TI*PCACO(L,NY,NX)+CORP*(FI*TPCACO &
      -TI*PCACO(L,NY,NX))+TX*PCACO(L,NY,NX)
      PCASO(L,NY,NX)=TI*PCASO(L,NY,NX)+CORP*(FI*TPCASO &
      -TI*PCASO(L,NY,NX))+TX*PCASO(L,NY,NX)
      PALPO(L,NY,NX)=TI*PALPO(L,NY,NX)+CORP*(FI*TPALPO &
      -TI*PALPO(L,NY,NX))+TX*PALPO(L,NY,NX)
      PFEPO(L,NY,NX)=TI*PFEPO(L,NY,NX)+CORP*(FI*TPFEPO &
      -TI*PFEPO(L,NY,NX))+TX*PFEPO(L,NY,NX)
      PCAPD(L,NY,NX)=TI*PCAPD(L,NY,NX)+CORP*(FI*TPCAPD &
      -TI*PCAPD(L,NY,NX))+TX*PCAPD(L,NY,NX)
      PCAPH(L,NY,NX)=TI*PCAPH(L,NY,NX)+CORP*(FI*TPCAPH &
      -TI*PCAPH(L,NY,NX))+TX*PCAPH(L,NY,NX)
      PCAPM(L,NY,NX)=TI*PCAPM(L,NY,NX)+CORP*(FI*TPCAPM &
      -TI*PCAPM(L,NY,NX))+TX*PCAPM(L,NY,NX)
      PALPB(L,NY,NX)=TI*PALPB(L,NY,NX)+CORP*(FI*TPALPB &
      -TI*PALPB(L,NY,NX))+TX*PALPB(L,NY,NX)
      PFEPB(L,NY,NX)=TI*PFEPB(L,NY,NX)+CORP*(FI*TPFEPB &
      -TI*PFEPB(L,NY,NX))+TX*PFEPB(L,NY,NX)
      PCPDB(L,NY,NX)=TI*PCPDB(L,NY,NX)+CORP*(FI*TPCPDB &
      -TI*PCPDB(L,NY,NX))+TX*PCPDB(L,NY,NX)
      PCPHB(L,NY,NX)=TI*PCPHB(L,NY,NX)+CORP*(FI*TPCPHB &
      -TI*PCPHB(L,NY,NX))+TX*PCPHB(L,NY,NX)
      PCPMB(L,NY,NX)=TI*PCPMB(L,NY,NX)+CORP*(FI*TPCPMB &
      -TI*PCPMB(L,NY,NX))+TX*PCPMB(L,NY,NX)
      CO2G(L,NY,NX)=TI*CO2G(L,NY,NX)+CORP*(FI*TCO2G-TI*CO2G(L,NY,NX)) &
      +TX*CO2G(L,NY,NX)
      CH4G(L,NY,NX)=TI*CH4G(L,NY,NX)+CORP*(FI*TCH4G-TI*CH4G(L,NY,NX)) &
      +TX*CH4G(L,NY,NX)
      CO2S(L,NY,NX)=TI*CO2S(L,NY,NX)+CORP*(FI*TCOZS-TI*CO2S(L,NY,NX)) &
      +TX*CO2S(L,NY,NX)+CORP*CO2SH(L,NY,NX)
      CH4S(L,NY,NX)=TI*CH4S(L,NY,NX)+CORP*(FI*TCHFS-TI*CH4S(L,NY,NX)) &
      +TX*CH4S(L,NY,NX)+CORP*CH4SH(L,NY,NX)
      OXYG(L,NY,NX)=TI*OXYG(L,NY,NX)+CORP*(FI*TOXYG-TI*OXYG(L,NY,NX)) &
      +TX*OXYG(L,NY,NX)
      OXYS(L,NY,NX)=TI*OXYS(L,NY,NX)+CORP*(FI*TOXYS-TI*OXYS(L,NY,NX)) &
      +TX*OXYS(L,NY,NX)+CORP*OXYSH(L,NY,NX)
      Z2GG(L,NY,NX)=TI*Z2GG(L,NY,NX)+CORP*(FI*TZ2GG-TI*Z2GG(L,NY,NX)) &
      +TX*Z2GG(L,NY,NX)
      Z2GS(L,NY,NX)=TI*Z2GS(L,NY,NX)+CORP*(FI*TZ2GS-TI*Z2GS(L,NY,NX)) &
      +TX*Z2GS(L,NY,NX)+CORP*Z2GSH(L,NY,NX)
      Z2OG(L,NY,NX)=TI*Z2OG(L,NY,NX)+CORP*(FI*TZ2OG-TI*Z2OG(L,NY,NX)) &
      +TX*Z2OG(L,NY,NX)
      Z2OS(L,NY,NX)=TI*Z2OS(L,NY,NX)+CORP*(FI*TZ2OS-TI*Z2OS(L,NY,NX)) &
      +TX*Z2OS(L,NY,NX)+CORP*Z2OSH(L,NY,NX)
      ZNH3G(L,NY,NX)=TI*ZNH3G(L,NY,NX)+CORP*(FI*TZNH3G &
      -TI*ZNH3G(L,NY,NX))+TX*ZNH3G(L,NY,NX)
      H2GG(L,NY,NX)=TI*H2GG(L,NY,NX)+CORP*(FI*TH2GG-TI*H2GG(L,NY,NX)) &
      +TX*H2GG(L,NY,NX)
      H2GS(L,NY,NX)=TI*H2GS(L,NY,NX)+CORP*(FI*TH2GS-TI*H2GS(L,NY,NX)) &
      +TX*H2GS(L,NY,NX)+CORP*H2GSH(L,NY,NX)
      ZNH4SH(L,NY,NX)=XCORP(NY,NX)*ZNH4SH(L,NY,NX)
      ZNH3SH(L,NY,NX)=XCORP(NY,NX)*ZNH3SH(L,NY,NX)
      ZNO3SH(L,NY,NX)=XCORP(NY,NX)*ZNO3SH(L,NY,NX)
      ZNO2SH(L,NY,NX)=XCORP(NY,NX)*ZNO2SH(L,NY,NX)
      H1PO4H(L,NY,NX)=XCORP(NY,NX)*H1PO4H(L,NY,NX)
      H2PO4H(L,NY,NX)=XCORP(NY,NX)*H2PO4H(L,NY,NX)
      ZNH4BH(L,NY,NX)=XCORP(NY,NX)*ZNH4BH(L,NY,NX)
      ZNH3BH(L,NY,NX)=XCORP(NY,NX)*ZNH3BH(L,NY,NX)
      ZNO3BH(L,NY,NX)=XCORP(NY,NX)*ZNO3BH(L,NY,NX)
      ZNO2BH(L,NY,NX)=XCORP(NY,NX)*ZNO2BH(L,NY,NX)
      H1POBH(L,NY,NX)=XCORP(NY,NX)*H1POBH(L,NY,NX)
      H2POBH(L,NY,NX)=XCORP(NY,NX)*H2POBH(L,NY,NX)
      ZALH(L,NY,NX)=XCORP(NY,NX)*ZALH(L,NY,NX)
      ZFEH(L,NY,NX)=XCORP(NY,NX)*ZFEH(L,NY,NX)
      ZHYH(L,NY,NX)=XCORP(NY,NX)*ZHYH(L,NY,NX)
      ZCCH(L,NY,NX)=XCORP(NY,NX)*ZCCH(L,NY,NX)
      ZMAH(L,NY,NX)=XCORP(NY,NX)*ZMAH(L,NY,NX)
      ZNAH(L,NY,NX)=XCORP(NY,NX)*ZNAH(L,NY,NX)
      ZKAH(L,NY,NX)=XCORP(NY,NX)*ZKAH(L,NY,NX)
      ZOHH(L,NY,NX)=XCORP(NY,NX)*ZOHH(L,NY,NX)
      ZSO4H(L,NY,NX)=XCORP(NY,NX)*ZSO4H(L,NY,NX)
      ZCLH(L,NY,NX)=XCORP(NY,NX)*ZCLH(L,NY,NX)
      ZCO3H(L,NY,NX)=XCORP(NY,NX)*ZCO3H(L,NY,NX)
      ZHCO3H(L,NY,NX)=XCORP(NY,NX)*ZHCO3H(L,NY,NX)
      ZALO1H(L,NY,NX)=XCORP(NY,NX)*ZALO1H(L,NY,NX)
      ZALO2H(L,NY,NX)=XCORP(NY,NX)*ZALO2H(L,NY,NX)
      ZALO3H(L,NY,NX)=XCORP(NY,NX)*ZALO3H(L,NY,NX)
      ZALO4H(L,NY,NX)=XCORP(NY,NX)*ZALO4H(L,NY,NX)
      ZALSH(L,NY,NX)=XCORP(NY,NX)*ZALSH(L,NY,NX)
      ZFEO1H(L,NY,NX)=XCORP(NY,NX)*ZFEO1H(L,NY,NX)
      ZFEO2H(L,NY,NX)=XCORP(NY,NX)*ZFEO2H(L,NY,NX)
      ZFEO3H(L,NY,NX)=XCORP(NY,NX)*ZFEO3H(L,NY,NX)
      ZFEO4H(L,NY,NX)=XCORP(NY,NX)*ZFEO4H(L,NY,NX)
      ZFESH(L,NY,NX)=XCORP(NY,NX)*ZFESH(L,NY,NX)
      ZCAOH(L,NY,NX)=XCORP(NY,NX)*ZCAOH(L,NY,NX)
      ZCACH(L,NY,NX)=XCORP(NY,NX)*ZCACH(L,NY,NX)
      ZCAHH(L,NY,NX)=XCORP(NY,NX)*ZCAHH(L,NY,NX)
      ZCASH(L,NY,NX)=XCORP(NY,NX)*ZCASH(L,NY,NX)
      ZMGOH(L,NY,NX)=XCORP(NY,NX)*ZMGOH(L,NY,NX)
      ZMGCH(L,NY,NX)=XCORP(NY,NX)*ZMGCH(L,NY,NX)
      ZMGHH(L,NY,NX)=XCORP(NY,NX)*ZMGHH(L,NY,NX)
      ZMGSH(L,NY,NX)=XCORP(NY,NX)*ZMGSH(L,NY,NX)
      ZNACH(L,NY,NX)=XCORP(NY,NX)*ZNACH(L,NY,NX)
      ZNASH(L,NY,NX)=XCORP(NY,NX)*ZNASH(L,NY,NX)
      ZKASH(L,NY,NX)=XCORP(NY,NX)*ZKASH(L,NY,NX)
      H0PO4H(L,NY,NX)=XCORP(NY,NX)*H0PO4H(L,NY,NX)
      H3PO4H(L,NY,NX)=XCORP(NY,NX)*H3PO4H(L,NY,NX)
      ZFE1PH(L,NY,NX)=XCORP(NY,NX)*ZFE1PH(L,NY,NX)
      ZFE2PH(L,NY,NX)=XCORP(NY,NX)*ZFE2PH(L,NY,NX)
      ZCA0PH(L,NY,NX)=XCORP(NY,NX)*ZCA0PH(L,NY,NX)
      ZCA1PH(L,NY,NX)=XCORP(NY,NX)*ZCA1PH(L,NY,NX)
      ZCA2PH(L,NY,NX)=XCORP(NY,NX)*ZCA2PH(L,NY,NX)
      ZMG1PH(L,NY,NX)=XCORP(NY,NX)*ZMG1PH(L,NY,NX)
      H0POBH(L,NY,NX)=XCORP(NY,NX)*H0POBH(L,NY,NX)
      H1POBH(L,NY,NX)=XCORP(NY,NX)*H1POBH(L,NY,NX)
      H3POBH(L,NY,NX)=XCORP(NY,NX)*H3POBH(L,NY,NX)
      ZFE1BH(L,NY,NX)=XCORP(NY,NX)*ZFE1BH(L,NY,NX)
      ZFE2BH(L,NY,NX)=XCORP(NY,NX)*ZFE2BH(L,NY,NX)
      ZCA0BH(L,NY,NX)=XCORP(NY,NX)*ZCA0BH(L,NY,NX)
      ZCA1BH(L,NY,NX)=XCORP(NY,NX)*ZCA1BH(L,NY,NX)
      ZCA2BH(L,NY,NX)=XCORP(NY,NX)*ZCA2BH(L,NY,NX)
      ZMG1BH(L,NY,NX)=XCORP(NY,NX)*ZMG1BH(L,NY,NX)
      CO2SH(L,NY,NX)=XCORP(NY,NX)*CO2SH(L,NY,NX)
      CH4SH(L,NY,NX)=XCORP(NY,NX)*CH4SH(L,NY,NX)
      OXYSH(L,NY,NX)=XCORP(NY,NX)*OXYSH(L,NY,NX)
      Z2GSH(L,NY,NX)=XCORP(NY,NX)*Z2GSH(L,NY,NX)
      Z2OSH(L,NY,NX)=XCORP(NY,NX)*Z2OSH(L,NY,NX)
      H2GSH(L,NY,NX)=XCORP(NY,NX)*H2GSH(L,NY,NX)

      DO 5965 K=0,jcplx1
      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      OMC(M,NGL,N,K,L,NY,NX)=TI*OMC(M,NGL,N,K,L,NY,NX)+CORP*(FI*TOMC(M,NGL,N,K) &
        -TI*OMC(M,NGL,N,K,L,NY,NX))+TX*OMC(M,NGL,N,K,L,NY,NX)
      OMN(M,NGL,N,K,L,NY,NX)=TI*OMN(M,NGL,N,K,L,NY,NX)+CORP*(FI*TOMN(M,NGL,N,K) &
        -TI*OMN(M,NGL,N,K,L,NY,NX))+TX*OMN(M,NGL,N,K,L,NY,NX)
      OMP(M,NGL,N,K,L,NY,NX)=TI*OMP(M,NGL,N,K,L,NY,NX)+CORP*(FI*TOMP(M,NGL,N,K) &
        -TI*OMP(M,NGL,N,K,L,NY,NX))+TX*OMP(M,NGL,N,K,L,NY,NX)
      enddo
      enddo
      enddo
5965  CONTINUE
      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      OMCff(M,NGL,N,L,NY,NX)=TI*OMCff(M,NGL,N,L,NY,NX)+CORP*(FI*TOMCff(M,NGL,N) &
        -TI*OMCff(M,NGL,N,L,NY,NX))+TX*OMCff(M,NGL,N,L,NY,NX)
      OMNff(M,NGL,N,L,NY,NX)=TI*OMNff(M,NGL,N,L,NY,NX)+CORP*(FI*TOMNff(M,NGL,N) &
        -TI*OMNff(M,NGL,N,L,NY,NX))+TX*OMNff(M,NGL,N,L,NY,NX)
      OMPff(M,NGL,N,L,NY,NX)=TI*OMPff(M,NGL,N,L,NY,NX)+CORP*(FI*TOMPff(M,NGL,N) &
        -TI*OMPff(M,NGL,N,L,NY,NX))+TX*OMPff(M,NGL,N,L,NY,NX)
      enddo
      enddo
      enddo

      DO 5980 K=0,jcplx1
      DO 5975 M=1,2
      ORC(M,K,L,NY,NX)=TI*ORC(M,K,L,NY,NX)+CORP*(FI*TORC(M,K) &
      -TI*ORC(M,K,L,NY,NX))+TX*ORC(M,K,L,NY,NX)
      ORN(M,K,L,NY,NX)=TI*ORN(M,K,L,NY,NX)+CORP*(FI*TORN(M,K) &
      -TI*ORN(M,K,L,NY,NX))+TX*ORN(M,K,L,NY,NX)
      ORP(M,K,L,NY,NX)=TI*ORP(M,K,L,NY,NX)+CORP*(FI*TORP(M,K) &
      -TI*ORP(M,K,L,NY,NX))+TX*ORP(M,K,L,NY,NX)
5975  CONTINUE
      OQC(K,L,NY,NX)=TI*OQC(K,L,NY,NX)+CORP*(FI*TOQC(K) &
      -TI*OQC(K,L,NY,NX))+TX*OQC(K,L,NY,NX)+CORP*OQCH(K,L,NY,NX)
      OQN(K,L,NY,NX)=TI*OQN(K,L,NY,NX)+CORP*(FI*TOQN(K) &
      -TI*OQN(K,L,NY,NX))+TX*OQN(K,L,NY,NX)+CORP*OQNH(K,L,NY,NX)
      OQP(K,L,NY,NX)=TI*OQP(K,L,NY,NX)+CORP*(FI*TOQP(K) &
      -TI*OQP(K,L,NY,NX))+TX*OQP(K,L,NY,NX)+CORP*OQPH(K,L,NY,NX)
      OQA(K,L,NY,NX)=TI*OQA(K,L,NY,NX)+CORP*(FI*TOQA(K) &
      -TI*OQA(K,L,NY,NX))+TX*OQA(K,L,NY,NX)+CORP*OQAH(K,L,NY,NX)
      OQCH(K,L,NY,NX)=XCORP(NY,NX)*OQCH(K,L,NY,NX)
      OQNH(K,L,NY,NX)=XCORP(NY,NX)*OQNH(K,L,NY,NX)
      OQPH(K,L,NY,NX)=XCORP(NY,NX)*OQPH(K,L,NY,NX)
      OQAH(K,L,NY,NX)=XCORP(NY,NX)*OQAH(K,L,NY,NX)
      OHC(K,L,NY,NX)=TI*OHC(K,L,NY,NX)+CORP*(FI*TOHC(K) &
      -TI*OHC(K,L,NY,NX))+TX*OHC(K,L,NY,NX)
      OHN(K,L,NY,NX)=TI*OHN(K,L,NY,NX)+CORP*(FI*TOHN(K) &
      -TI*OHN(K,L,NY,NX))+TX*OHN(K,L,NY,NX)
      OHP(K,L,NY,NX)=TI*OHP(K,L,NY,NX)+CORP*(FI*TOHP(K) &
      -TI*OHP(K,L,NY,NX))+TX*OHP(K,L,NY,NX)
      OHA(K,L,NY,NX)=TI*OHA(K,L,NY,NX)+CORP*(FI*TOHA(K) &
      -TI*OHA(K,L,NY,NX))+TX*OHA(K,L,NY,NX)
      DO 5970 M=1,4
      OSC(M,K,L,NY,NX)=TI*OSC(M,K,L,NY,NX)+CORP*(FI*TOSC(M,K) &
      -TI*OSC(M,K,L,NY,NX))+TX*OSC(M,K,L,NY,NX)
      OSA(M,K,L,NY,NX)=TI*OSA(M,K,L,NY,NX)+CORP*(FI*TOSA(M,K) &
      -TI*OSA(M,K,L,NY,NX))+TX*OSA(M,K,L,NY,NX)
      OSN(M,K,L,NY,NX)=TI*OSN(M,K,L,NY,NX)+CORP*(FI*TOSN(M,K) &
      -TI*OSN(M,K,L,NY,NX))+TX*OSN(M,K,L,NY,NX)
      OSP(M,K,L,NY,NX)=TI*OSP(M,K,L,NY,NX)+CORP*(FI*TOSP(M,K) &
      -TI*OSP(M,K,L,NY,NX))+TX*OSP(M,K,L,NY,NX)
5970  CONTINUE
5980  CONTINUE
!
!     ADD STATE VARIABLES IN SURFACE RESIDUE INCORPORATED
!     WITHIN TILLAGE MIXING ZONE
!
      DO 5910 K=0,jcplx1
      IF(K.NE.3.AND.K.NE.4)THEN
      DO 5915 N=1,7
      DO M=1,3
      DO NGL=1,JG
      OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)+FI*TOMGC(M,NGL,N,K)
      OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)+FI*TOMGN(M,NGL,N,K)
      OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)+FI*TOMGP(M,NGL,N,K)
      enddo
      enddo
5915  CONTINUE
      ENDIF
5910  CONTINUE
      DO  N=1,7
      DO M=1,3
      DO NGL=1,JG
      OMCff(M,NGL,N,L,NY,NX)=OMCff(M,NGL,N,L,NY,NX)+FI*TOMGCff(M,NGL,N)
      OMNff(M,NGL,N,L,NY,NX)=OMNff(M,NGL,N,L,NY,NX)+FI*TOMGNff(M,NGL,N)
      OMPff(M,NGL,N,L,NY,NX)=OMPff(M,NGL,N,L,NY,NX)+FI*TOMGPff(M,NGL,N)
      enddo
      enddo
      ENDDO

      DO 5920 K=0,2
      DO 5925 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+FI*TORXC(M,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+FI*TORXN(M,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+FI*TORXP(M,K)
5925  CONTINUE
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+FI*TOQGC(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+FI*TOQGN(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+FI*TOQGP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+FI*TOQGA(K)
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)+FI*TOQHC(K)
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)+FI*TOQHN(K)
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)+FI*TOQHP(K)
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)+FI*TOQHA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+FI*TOHGC(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+FI*TOHGN(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+FI*TOHGP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+FI*TOHGA(K)
      DO 5930 M=1,4
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)+FI*TOSGC(M,K)
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)+FI*TOSGA(M,K)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+FI*TOSGN(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+FI*TOSGP(M,K)
5930  CONTINUE
5920  CONTINUE
      OC=0.0_r8
      ON=0.0_r8
      OP=0.0_r8
      DC=0.0_r8
      DN=0.0_r8
      DP=0.0_r8

      DO 5985 K=0,jcplx1
      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      OC=OC+OMC(M,NGL,N,K,L,NY,NX)
      ON=ON+OMN(M,NGL,N,K,L,NY,NX)
      OP=OP+OMP(M,NGL,N,K,L,NY,NX)
      enddo
      enddo
      enddo
5985  CONTINUE

      DO  K=0,2
      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      DC=DC+OMC(M,NGL,N,K,L,NY,NX)
      DN=DN+OMN(M,NGL,N,K,L,NY,NX)
      DP=DP+OMP(M,NGL,N,K,L,NY,NX)
      enddo
      enddo
      enddo
   ENDDO
      DO  N=1,7
      DO  M=1,3
      DO NGL=1,JG
      OC=OC+OMCff(M,NGL,N,L,NY,NX)
      ON=ON+OMNff(M,NGL,N,L,NY,NX)
      OP=OP+OMPff(M,NGL,N,L,NY,NX)
      enddo
      enddo
      enddo

      DO 6995 K=0,jcplx1
      DO 6985 M=1,2
      OC=OC+ORC(M,K,L,NY,NX)
      ON=ON+ORN(M,K,L,NY,NX)
      OP=OP+ORP(M,K,L,NY,NX)
      IF(K.LE.2)THEN
      DC=DC+ORC(M,K,L,NY,NX)
      DN=DN+ORN(M,K,L,NY,NX)
      DP=DP+ORP(M,K,L,NY,NX)
      ENDIF
6985  CONTINUE

      OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
      +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      IF(K.LE.2)THEN
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
      +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
      ENDIF
      DO 6980 M=1,4
      OC=OC+OSC(M,K,L,NY,NX)
      ON=ON+OSN(M,K,L,NY,NX)
      OP=OP+OSP(M,K,L,NY,NX)
      IF(K.LE.2)THEN
      DC=DC+OSC(M,K,L,NY,NX)
      DN=DN+OSN(M,K,L,NY,NX)
      DP=DP+OSP(M,K,L,NY,NX)
      ENDIF
6980  CONTINUE
6995  CONTINUE
      ORGC(L,NY,NX)=OC
      ORGN(L,NY,NX)=ON
      ORGR(L,NY,NX)=DC
      CO2S(L,NY,NX)=CO2S(L,NY,NX)+FI*TCO2GS
      CH4S(L,NY,NX)=CH4S(L,NY,NX)+FI*TCH4GS
      OXYS(L,NY,NX)=OXYS(L,NY,NX)+FI*TOXYGS
      Z2GS(L,NY,NX)=Z2GS(L,NY,NX)+FI*TZ2GSG
      Z2OS(L,NY,NX)=Z2OS(L,NY,NX)+FI*TZ2OGS
      H2GS(L,NY,NX)=H2GS(L,NY,NX)+FI*TH2GGS
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+FI*TNH4GS
      ZNH3S(L,NY,NX)=ZNH3S(L,NY,NX)+FI*TNH3GS
      ZNO3S(L,NY,NX)=ZNO3S(L,NY,NX)+FI*TNO3GS
      ZNO2S(L,NY,NX)=ZNO2S(L,NY,NX)+FI*TNO2GS
      H1PO4(L,NY,NX)=H1PO4(L,NY,NX)+FI*TP14GS
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+FI*TPO4GS
      XN4(L,NY,NX)=XN4(L,NY,NX)+FI*TXN4G
      XOH0(L,NY,NX)=XOH0(L,NY,NX)+FI*TXOH0G
      XOH1(L,NY,NX)=XOH1(L,NY,NX)+FI*TXOH1G
      XOH2(L,NY,NX)=XOH2(L,NY,NX)+FI*TXOH2G
      XH1P(L,NY,NX)=XH1P(L,NY,NX)+FI*TXH1PG
      XH2P(L,NY,NX)=XH2P(L,NY,NX)+FI*TXH2PG
      PALPO(L,NY,NX)=PALPO(L,NY,NX)+FI*TALPOG
      PFEPO(L,NY,NX)=PFEPO(L,NY,NX)+FI*TFEPOG
      PCAPD(L,NY,NX)=PCAPD(L,NY,NX)+FI*TCAPDG
      PCAPH(L,NY,NX)=PCAPH(L,NY,NX)+FI*TCAPHG
      PCAPM(L,NY,NX)=PCAPM(L,NY,NX)+FI*TCAPMG
      ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)+FI*TNH4FG
      ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)+FI*TNH3FG
      ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)+FI*TNHUFG
      ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)+FI*TNO3FG
      ZNHU0(L,NY,NX)=ZNHUX0
      ZNHUI(L,NY,NX)=ZNHUXI
      ZNFN0(L,NY,NX)=ZNFNX0
      ZNFNI(L,NY,NX)=(TI*ZNFNI(L,NY,NX)+CORP*(FI*TZNFNI &
      -TI*ZNFNI(L,NY,NX))+TX*ZNFNI(L,NY,NX)+FI*TZNFNG)/FI
      TZNFN2=TZNFN2+ZNFNI(L,NY,NX)
      ENDIF
2000  CONTINUE
      ZNFN0(0,NY,NX)=ZNFNX0
      ZNFNI(0,NY,NX)=ZNFNI(0,NY,NX)*XCORP0
      TZNFN2=TZNFN2+TZNFNG
      TZNFNI=TZNFNI+TZNFNG
      DO 2001 L=NU(NY,NX),LL
      IF(TZNFN2.GT.ZERO)THEN
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*TZNFNI/TZNFN2
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX) &
      +0.5*(ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX))
      ENDIF
2001  CONTINUE
      ENDIF
      end subroutine ApplyMixing
!------------------------------------------------------------------------------------------

  subroutine AddFluxToSurfaceResidue(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: M,N,K,NGL
!     begin_execution
!     ADD ABOVE-GROUND LITTERFALL FROM EXTRACT.F TO SURFACE RESIDUE
!
!     OSC,OSN,OSP=SOC,SON,SOP
!     CSNT,ZSNT,PSNT=total C,N,P litterfall
!     ORGC=total SOC
!     RAINR,HRAINR=water,heat in litterfall
!     FLWR,HFLWR=water,heat flux into litter
!     HEATIN=cumulative net surface heat transfer
!
  DO 6965 K=0,1
    DO  M=1,jsken
      OSC(M,K,0,NY,NX)=OSC(M,K,0,NY,NX)+CSNT(M,K,0,NY,NX)
      DO NGL=1,JG
        OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)+CSNT(M,K,0,NY,NX)*micpar%OMCI(1+(NGL-1)*3,K)
      ENDDO
      OSN(M,K,0,NY,NX)=OSN(M,K,0,NY,NX)+ZSNT(M,K,0,NY,NX)
      OSP(M,K,0,NY,NX)=OSP(M,K,0,NY,NX)+PSNT(M,K,0,NY,NX)
      ORGC(0,NY,NX)=ORGC(0,NY,NX)+CSNT(M,K,0,NY,NX)
      RAINR=CSNT(M,K,0,NY,NX)*THETCX(K)
      HRAINR=RAINR*cpw*TKA(NY,NX)
      FLWR(NY,NX)=FLWR(NY,NX)+RAINR
      HFLWR(NY,NX)=HFLWR(NY,NX)+HRAINR
      !if(HFLWR(NY,NX)<-1.8_r8)then
      !  write(*,*)'HFLWR(NY,NX)+HRAINR',HFLWR(NY,NX),HRAINR
      !  write(*,*)'RAINR*cpw*TKA(NY,NX)',RAINR,TKA(NY,NX)
      !  call print_info(trim(mod_filename)//' at line',__LINE__)
      !endif
      CRAIN=CRAIN+RAINR
      HEATIN=HEATIN+HRAINR
    enddo
6965  CONTINUE
  !
  !     GAS AND SOLUTE EXCHANGE WITHIN SURFACE LITTER ADDED TO ECOSYSTEM
  !     TOTALS FOR CALCULATING COMPETITION CONSTRAINTS ON MICROBIAL
  !     AND ROOT POPULATIONS IN NITRO.F AND UPTAKE.F
  !
  !     ROXYX=O2 demand by all microbial,root,myco populations
  !     RNH4X=NH4 demand in non-band by all microbial,root,myco populations
  !     RNO3X=NO3 demand in non-band by all microbial,root,myco populations
  !     RPO4X=H2PO4 demand in non-band by all microbial,root,myco populations
  !     RP14X=HPO4 demand in non-band by all microbial,root,myco populations
  !     RNHBX=NH4 demand in band by all microbial,root,myco populations
  !     RN3BX=NO3 demand in band by all microbial,root,myco populations
  !     RPOBX=H2PO4 demand in band by all microbial,root,myco populations
  !     RP1BX=HPO4 demand in band by all microbial,root,myco populations
  !     ROXYS=O2 demand from DOC,DOA oxidation
  !     RVMX4=demand for NH4 oxidation
  !     RVMX3=demand for NO3 reduction
  !     RVMX2=demand for NO2 oxidation
  !     RVMX1=demand for N2O reduction
  !     RINHO,RINHOR=substrate-unlimited NH4 immobilization
  !     RINOO,RINOOR=substrate-unlimited NO3 immobilization
  !     RIPOO,RIPOOR=substrate-unlimited H2PO4 immobilization
  !     RIPO1,RIPO1R=substrate-unlimited HPO4 immobilization
  !     ROQCX,ROQAX=total DOC,DOA demand from DOC,DOA oxidation
  !     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
  !     RNO2X=total demand for NO2 reduction
  !     RVMXC=demand for NO2 reduction
!
  DO 8990 K=0,jcplx1
    IF(K.NE.3.AND.K.NE.4)THEN
      DO 8980 N=1,NFGs
        DO NGL=1,JG
          ROXYX(0,NY,NX)=ROXYX(0,NY,NX)+ROXYS(NGL,N,K,0,NY,NX)
          RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RVMX4(NGL,N,K,0,NY,NX)
          RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RVMX3(NGL,N,K,0,NY,NX)
          RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMX2(NGL,N,K,0,NY,NX)
          RN2OX(0,NY,NX)=RN2OX(0,NY,NX)+RVMX1(NGL,N,K,0,NY,NX)
          RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RINHO(NGL,N,K,0,NY,NX)
          RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RINOO(NGL,N,K,0,NY,NX)
          RPO4X(0,NY,NX)=RPO4X(0,NY,NX)+RIPOO(NGL,N,K,0,NY,NX)
          RP14X(0,NY,NX)=RP14X(0,NY,NX)+RIPO1(NGL,N,K,0,NY,NX)
          RNH4X(NU(NY,NX),NY,NX)=RNH4X(NU(NY,NX),NY,NX)+RINHOR(NGL,N,K,NY,NX)
          RNO3X(NU(NY,NX),NY,NX)=RNO3X(NU(NY,NX),NY,NX)+RINOOR(NGL,N,K,NY,NX)
          RPO4X(NU(NY,NX),NY,NX)=RPO4X(NU(NY,NX),NY,NX)+RIPOOR(NGL,N,K,NY,NX)
          RP14X(NU(NY,NX),NY,NX)=RP14X(NU(NY,NX),NY,NX)+RIPO1R(NGL,N,K,NY,NX)
          ROQCX(K,0,NY,NX)=ROQCX(K,0,NY,NX)+ROQCS(NGL,N,K,0,NY,NX)
          ROQAX(K,0,NY,NX)=ROQAX(K,0,NY,NX)+ROQAS(NGL,N,K,0,NY,NX)
        ENDDO
8980  CONTINUE
    ENDIF
8990  CONTINUE

      DO  N=1,NFGs
        DO NGL=1,JG
          ROXYX(0,NY,NX)=ROXYX(0,NY,NX)+ROXYSff(NGL,N,0,NY,NX)
          RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RVMX4ff(NGL,N,0,NY,NX)
          RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RVMX3ff(NGL,N,0,NY,NX)
          RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMX2ff(NGL,N,0,NY,NX)
          RN2OX(0,NY,NX)=RN2OX(0,NY,NX)+RVMX1ff(NGL,N,0,NY,NX)
          RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RINHOff(NGL,N,0,NY,NX)
          RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RINOOff(NGL,N,0,NY,NX)
          RPO4X(0,NY,NX)=RPO4X(0,NY,NX)+RIPOOff(NGL,N,0,NY,NX)
          RP14X(0,NY,NX)=RP14X(0,NY,NX)+RIPO1ff(NGL,N,0,NY,NX)
          RNH4X(NU(NY,NX),NY,NX)=RNH4X(NU(NY,NX),NY,NX)+RINHORff(NGL,N,NY,NX)
          RNO3X(NU(NY,NX),NY,NX)=RNO3X(NU(NY,NX),NY,NX)+RINOORff(NGL,N,NY,NX)
          RPO4X(NU(NY,NX),NY,NX)=RPO4X(NU(NY,NX),NY,NX)+RIPOORff(NGL,N,NY,NX)
          RP14X(NU(NY,NX),NY,NX)=RP14X(NU(NY,NX),NY,NX)+RIPO1Rff(NGL,N,NY,NX)

        ENDDO
      ENDDO

  RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMXC(0,NY,NX)
  end subroutine AddFluxToSurfaceResidue

end module RedistMod
