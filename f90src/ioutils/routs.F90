  SUBROUTINE routs(NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE READS SOIL CHECKPOINT FILES TO
!     RE-INITILIAZE THE MODEL FROM A SELECTED DATE IN AN EARLIER RUN
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StartsMod    , only : starts
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use LandSurfDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use PlantTraitDataType
  use PlantDataRateType
  use SnowDataType
  USE CanopyDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use SoilBGCDataType
  use EcoSimSumDataType
  use PlantMngmtDataType
  use RootDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use AqueChemDatatype
  use SedimentDataType
  use GridDataType
  use EcoSIMConfig
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: IDATE,IYR,K,M,NX,NY,N,L,NGL
!     execution begins here

  REWIND(21)
  REWIND(22)
  REWIND(23)
  REWIND(24)
  REWIND(25)
!
  !     INITIALIZE MANAGEMENT ARRAYS
!
  IF(IMNG.EQ.0)THEN
    DO 7980 NX=NHW,NHE
      DO 7975 NY=NVN,NVS
        IFLGS(NY,NX)=1
        DO 7990 M=1,366
          ITILL(M,NY,NX)=0
          DCORP(M,NY,NX)=0.0
          DO 7985 N=1,24
            RRIG(N,M,NY,NX)=0.0
7985      CONTINUE
          DO 7995 N=1,20
            FERT(N,M,NY,NX)=0.0
7995      CONTINUE
7990    CONTINUE
7975  CONTINUE
7980  CONTINUE
  ENDIF
!
! READ CHECKPOINT FILES UP TO DATE OF RE-INITIALIZATION
!
8000  READ(21,90,END=1001)IDATE,IYR,CRAIN,TSEDOU &
      ,HEATIN,OXYGIN,TORGF,TORGN,TORGP,CO2GIN,ZN2GIN,VOLWOU,CEVAP &
      ,CRUN,HEATOU,OXYGOU,TCOU,TZOU,TPOU,TZIN,TPIN,XCSN,XZSN,XPSN
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      READ(21,95)IDATE,IYR,(TDTPX(NY,NX,N),N=1,12) &
        ,(TDTPN(NY,NX,N),N=1,12),(TDRAD(NY,NX,N),N=1,12) &
        ,(TDWND(NY,NX,N),N=1,12),(TDHUM(NY,NX,N),N=1,12) &
        ,(TDPRC(NY,NX,N),N=1,12),(TDIRI(NY,NX,N),N=1,12) &
        ,(TDCO2(NY,NX,N),N=1,12),(TDCN4(NY,NX,N),N=1,12) &
        ,(TDCNO(NY,NX,N),N=1,12)
      READ(21,93)IDATE,IYR,IFLGT(NY,NX),IFNHB(NY,NX),IDTBL(NY,NX) &
        ,IFNOB(NY,NX),IFPOB(NY,NX),IUTYP(NY,NX) &
        ,ZT(NY,NX),TFLWC(NY,NX),TSED(NY,NX) &
        ,ZS(NY,NX),THRMC(NY,NX),THRMG(NY,NX),TCNET(NY,NX),TVOLWC(NY,NX) &
        ,VOLWG(NY,NX),URAIN(NY,NX),ARLFC(NY,NX),ARSTC(NY,NX),PPT(NY,NX) &
        ,VOLWD(NY,NX),ZM(NY,NX),UCO2G(NY,NX),UCH4G(NY,NX),UOXYG(NY,NX) &
        ,UN2GG(NY,NX),UN2OG(NY,NX),UNH3G(NY,NX),UN2GS(NY,NX) &
        ,UCO2F(NY,NX),UCH4F(NY,NX),UOXYF(NY,NX),UN2OF(NY,NX) &
        ,UNH3F(NY,NX),UPO4F(NY,NX),UORGF(NY,NX),UFERTN(NY,NX) &
        ,UFERTP(NY,NX),UVOLO(NY,NX),UEVAP(NY,NX),URUN(NY,NX) &
        ,UDOCQ(NY,NX),UDOCD(NY,NX),UDONQ(NY,NX),UDOND(NY,NX) &
        ,UDOPQ(NY,NX),UDOPD(NY,NX),UDICQ(NY,NX),UDICD(NY,NX) &
        ,UDINQ(NY,NX),UDIND(NY,NX),UDIPQ(NY,NX),UDIPD(NY,NX) &
        ,UXCSN(NY,NX),UXZSN(NY,NX),UXPSN(NY,NX),UCOP(NY,NX) &
        ,UDRAIN(NY,NX),ZDRAIN(NY,NX),PDRAIN(NY,NX),DPNH4(NY,NX) &
        ,DPNO3(NY,NX),DPPO4(NY,NX),DETS(NY,NX),DETE(NY,NX) &
        ,CER(NY,NX),XER(NY,NX),USEDOU(NY,NX),ROWN(NY,NX),ROWO(NY,NX) &
        ,ROWP(NY,NX),DTBLZ(NY,NX),DTBLD(NY,NX),TNBP(NY,NX),VOLR(NY,NX) &
        ,SED(NY,NX),TGPP(NY,NX),TRAU(NY,NX),TNPP(NY,NX),THRE(NY,NX) &
        ,TCAN(NY,NX),TLEC(NY,NX),TSHC(NY,NX),DYLN(NY,NX),DYLX(NY,NX) &
        ,DTBLX(NY,NX),DTBLY(NY,NX),(RC0(K,NY,NX),K=1,jcplx),RC0ff(NY,NX) &
        ,VOLSS(NY,NX),VOLWS(NY,NX),VOLIS(NY,NX),VOLS(NY,NX) &
        ,DPTHS(NY,NX),ENGYP(NY,NX)
      READ(21,91)IDATE,IYR,(VOLSSL(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(VOLISL(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(VOLWSL(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(VOLSL(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(DENSS(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(DLYRS(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(VHCPW(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(TKW(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(TCW(L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcg_solsml(idg_CO2,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcg_solsml(idg_CH4,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcg_solsml(idg_O2,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcg_solsml(idg_N2,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcg_solsml(idg_N2O,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcn_solsml(ids_NH4,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcg_solsml(idg_NH3,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcn_solsml(ids_NO3,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcn_solsml(ids_H1PO4,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(trcn_solsml(ids_H2PO4,L,NY,NX),L=1,JS)
      READ(21,91)IDATE,IYR,(FHOL(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(DLYR(3,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(CDPTH(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(CDPTHZ(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(BKDSI(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(BKDS(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(CORGC(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(POROS(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(FC(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(WP(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(SCNV(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(SCNH(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(SAND(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(SILT(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(CLAY(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLW(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLWX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLI(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLP(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLA(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLY(L,NY,NX),L=0,NLI(NY,NX))
!     READ(21,91)IDATE,IYR,(GKC4(L,NY,NX),L=1,NLI(NY,NX))
!     READ(21,91)IDATE,IYR,(GKCH(L,NY,NX),L=1,NLI(NY,NX))
!     READ(21,91)IDATE,IYR,(GKCA(L,NY,NX),L=1,NLI(NY,NX))
!     READ(21,91)IDATE,IYR,(GKCM(L,NY,NX),L=1,NLI(NY,NX))
!     READ(21,91)IDATE,IYR,(GKCN(L,NY,NX),L=1,NLI(NY,NX))
!     READ(21,91)IDATE,IYR,(GKCK(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XCEC(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XAEC(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(TCS(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(TKS(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VHCP(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VHCM(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_gasml(idg_CO2,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_solml(idg_CO2,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_soHml(idg_CO2,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_gasml(idg_CH4,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_solml(idg_CH4,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_soHml(idg_CH4,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(ROXYF(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RCO2F(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(ROXYL(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RCH4F(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RCH4L(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_gasml(idg_H2,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_solml(idg_H2,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_gasml(idg_O2,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_solml(idg_O2,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trc_soHml(idg_O2,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(ROXYX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RNH4X(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RNO3X(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RNO2X(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RN2OX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RPO4X(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RP14X(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RNHBX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RN3BX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RN2BX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RPOBX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RP1BX(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,95)IDATE,IYR,((ROQCX(K,L,NY,NX),L=0,NLI(NY,NX)),K=1,jcplx)
      READ(21,95)IDATE,IYR,((ROQAX(K,L,NY,NX),L=0,NLI(NY,NX)),K=1,jcplx)
      READ(21,91)IDATE,IYR,(VOLWH(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLIH(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(VOLAH(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(RTDNT(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(ZL(L,NY,NX),L=0,JC)
      READ(21,91)IDATE,IYR,(ARLFT(L,NY,NX),L=1,JC)
      READ(21,91)IDATE,IYR,(ARSTT(L,NY,NX),L=1,JC)
      READ(21,91)IDATE,IYR,(WGLFT(L,NY,NX),L=1,JC)
      READ(21,91)IDATE,IYR,(trcs_VLN(ids_NH4,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcs_VLN(ids_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcs_VLN(ids_NO3,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcs_VLN(ids_NO3B,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcs_VLN(ids_H1PO4,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcs_VLN(ids_H1PO4B,L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XN4(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XNB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XOH0(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XOH1(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XOH2(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XH1P(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XH2P(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XOH0B(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XOH1B(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XOH2B(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XH1PB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(XH2PB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PALPO(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PFEPO(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PCAPD(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PCAPH(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PCAPM(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PALPB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PFEPB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PCPDB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PCPHB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(PCPMB(L,NY,NX),L=0,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(CION(L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Al,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Fe,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Hp,L,NY,NX),L=1,NLI(NY,NX))
      READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Ca,L,NY,NX),L=1,NLI(NY,NX))
      IF(ISALTG.NE.0)THEN
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Al,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Fe,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Hp,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Ca,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Mg,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Na,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_K,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_OH,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_SO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_Cl,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_HCO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_AlOH,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_AlOH2,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_AlOH3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_AlOH4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_AlSO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeOH,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeOH2,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeOH3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeOH4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeSO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaOH2,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaCO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaHCO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaSO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_MgOH2,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_MgCO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_MgHCO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_MgSO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_NaCO3,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_NaSO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_KSO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_H0PO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_H3PO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeHPO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_FeH2PO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaPO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaHPO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_CaH2PO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcs_solsml(idsa_MgHPO4,L,NY,NX),L=1,JS)
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Mg,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Na,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_K,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_OH,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_SO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_Cl,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_HCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_AlOH,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_AlOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_AlOH4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_AlSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeOH,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeOH4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaHCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_NaCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_NaSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_KSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_H0PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_H3PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeHPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaHPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgHPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_H0PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_H3PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_FeH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_CaH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_solml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Al,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Fe,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Hp,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Ca,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Mg,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Na,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_K,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_OH,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_SO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_Cl,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_HCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_AlOH,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_AlOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_AlOH4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_AlSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeOH,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeOH4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaHCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_MgOH2,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_MgCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_MgSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_NaCO3,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_NaSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_KSO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_H0PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_H3PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeHPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaHPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_MgHPO4,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_H0PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_H3PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_FeH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_CaH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(trcsa_soHml(idsa_MgHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XHY(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XAL(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XFE(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XCA(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XMG(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XNA(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(XKA(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(PALOH(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(PFEOH(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(PCACO(L,NY,NX),L=1,NLI(NY,NX))
        READ(21,91)IDATE,IYR,(PCASO(L,NY,NX),L=1,NLI(NY,NX))
      ENDIF
      D9985: DO K=1,jcplx
        DO  N=1,NFGs
          DO NGL=JGnio(N),JGnfo(N)
            READ(22,91)IDATE,IYR,(ROXYS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMX4(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMX3(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMX2(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMX1(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMB4(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMB3(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RVMB2(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RINHO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RINOO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RIPOO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RINHB(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RINOB(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(RIPBO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(ROQCS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(ROQAS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,RINHOR(NGL,K,NY,NX)
            READ(22,91)IDATE,IYR,RINOOR(NGL,K,NY,NX)
            READ(22,91)IDATE,IYR,RIPOOR(NGL,K,NY,NX)

            DO  M=1,nlbiomcp
              READ(22,91)IDATE,IYR,(OMC(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
              READ(22,91)IDATE,IYR,(OMN(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
              READ(22,91)IDATE,IYR,(OMP(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            enddo
          enddo
        enddo
      ENDDO D9985

      DO  N=1,NFGs
        DO NGL=JGniA(N),JGnfA(N)
          READ(22,91)IDATE,IYR,(ROXYSff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMX4ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMX3ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMX2ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMX1ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMB4ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMB3ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RVMB2ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RINHOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RINOOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RIPOOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RINHBff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RINOBff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(RIPBOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,RINHORff(NGL,NY,NX)
          READ(22,91)IDATE,IYR,RINOORff(NGL,NY,NX)
          READ(22,91)IDATE,IYR,RIPOORff(NGL,NY,NX)

          DO  M=1,nlbiomcp
            READ(22,91)IDATE,IYR,(OMCff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(OMNff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
            READ(22,91)IDATE,IYR,(OMPff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
          enddo
        ENDDO
      ENDDO

      D9980: DO K=1,jcplx
        D9975: DO M=1,ndbiomcp
          READ(22,91)IDATE,IYR,(ORC(M,K,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(ORN(M,K,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(ORP(M,K,L,NY,NX),L=0,NLI(NY,NX))
        ENDDO D9975
        READ(22,91)IDATE,IYR,(OQC(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQN(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQP(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQA(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQCH(K,L,NY,NX),L=1,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQNH(K,L,NY,NX),L=1,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQPH(K,L,NY,NX),L=1,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OQAH(K,L,NY,NX),L=1,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OHC(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OHN(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OHP(K,L,NY,NX),L=0,NLI(NY,NX))
        READ(22,91)IDATE,IYR,(OHA(K,L,NY,NX),L=0,NLI(NY,NX))

        D9970: DO M=1,jsken
          READ(22,91)IDATE,IYR,(OSC(M,K,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(OSA(M,K,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(OSN(M,K,L,NY,NX),L=0,NLI(NY,NX))
          READ(22,91)IDATE,IYR,(OSP(M,K,L,NY,NX),L=0,NLI(NY,NX))
        ENDDO D9970
      ENDDO D9980

      READ(22,91)IDATE,IYR,(RVMXC(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ORGC(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ORGR(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNH4FA(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNH3FA(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNHUFA(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNO3FA(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNH4FB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNH3FB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNHUFB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNO3FB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_gasml(idg_N2,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(idg_N2,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(idg_N2,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_gasml(idg_N2O,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(idg_N2O,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(idg_N2O,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_gasml(idg_NH3,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_NH4,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_NH4,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(idg_NH3,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(idg_NH3,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_NO3,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_NO3,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_NO2,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_NO2,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_H1PO4,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_H2PO4,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_H1PO4,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_H2PO4,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_NH4B,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(idg_NH3B,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(idg_NH3B,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_NO3B,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_NO3B,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_NO2B,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_NO2B,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_H1PO4B,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_solml(ids_H2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_H1PO4B,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(trc_soHml(ids_H2PO4B,L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(WDNHB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(DPNHB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(WDNOB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(DPNOB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(WDPOB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(DPPOB(L,NY,NX),L=1,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNHUI(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNHU0(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNFNI(L,NY,NX),L=0,NLI(NY,NX))
      READ(22,91)IDATE,IYR,(ZNFN0(L,NY,NX),L=0,NLI(NY,NX))
    ENDDO D9990
  ENDDO D9995
  IF(IDATE.LT.IDAYR.OR.IYR.LT.IYRR)THEN
    GO TO 8000
  ELSEIF(IDATE.GE.IDAYR.AND.IYR.EQ.IYRR)THEN
    IDAYR=IDATE
    ISTART=IDATE+1
    GO TO 1000
  ENDIF
1001  CONTINUE
  IDAYR=IDATE
  ISTART=IDATE+1
  IYR=IYRR
  DATA1(20)='NO'
  is_restart_run=(DATA1(20)=='YES')
  CALL STARTS(NHW,NHE,NVN,NVS)
  GO TO 1000
1000  RETURN
90    FORMAT(2I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
91    FORMAT(2I4,21E17.8E3)
92    FORMAT(2I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
93    FORMAT(8I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
95    FORMAT(2I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
  RETURN
  END SUBROUTINE routs
