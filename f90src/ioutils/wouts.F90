  SUBROUTINE wouts(I,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES OUT ALL SOIL VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcoSIMConfig, only : jcplx=> jcplxc, NFGs=> NFGsc,nlbiomcp=>nlbiomcpc
  use EcoSIMConfig, only : ndbiomcp=>ndbiomcpc,jsken=>jskenc
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use ClimForcDataType
  use LandSurfDataType
  use PlantTraitDataType
  use PlantDataRateType
  use SnowDataType
  use SurfLitterDataType
  use SurfSoilDataType
  use CanopyDataType
  use SoilBGCDataType
  use PlantMngmtDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use AqueChemDatatype
  use SedimentDataType
  use GridDataType

  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: K,M,NX,NY,N,L,NGL
!     execution begins here

  WRITE(21,90)I,IDATA(3),CRAIN,TSEDOU &
    ,HEATIN,OXYGIN,TORGF,TORGN,TORGP,CO2GIN,ZN2GIN,VOLWOU,CEVAP &
    ,CRUN,HEATOU,OXYGOU,TCOU,TZOU,TPOU,TZIN,TPIN,XCSN,XZSN,XPSN
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      WRITE(21,95)I,IDATA(3),(TDTPX(NY,NX,N),N=1,12) &
        ,(TDTPN(NY,NX,N),N=1,12),(TDRAD(NY,NX,N),N=1,12) &
        ,(TDWND(NY,NX,N),N=1,12),(TDHUM(NY,NX,N),N=1,12) &
        ,(TDPRC(NY,NX,N),N=1,12),(TDIRI(NY,NX,N),N=1,12) &
        ,(TDCO2(NY,NX,N),N=1,12),(TDCN4(NY,NX,N),N=1,12) &
        ,(TDCNO(NY,NX,N),N=1,12)
      WRITE(21,93)I,IDATA(3),IFLGT(NY,NX),IFNHB(NY,NX),IDTBL(NY,NX) &
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
      WRITE(21,91)I,IDATA(3),(VOLSSL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(VOLISL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(VOLWSL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(VOLSL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(DENSS(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(DLYRS(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(VHCPW(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(TKW(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(TCW(L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcg_solsml(idg_CO2,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcg_solsml(idg_CH4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcg_solsml(idg_O2,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcg_solsml(idg_N2,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcg_solsml(idg_N2O,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcn_solsml(ids_NH4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcg_solsml(idg_NH3,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcn_solsml(ids_NO3,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcn_solsml(ids_H1PO4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(trcn_solsml(ids_H2PO4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,IDATA(3),(FHOL(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(DLYR(3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(CDPTH(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(CDPTHZ(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(BKDSI(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(BKDS(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(CORGC(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(POROS(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(FC(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(WP(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(SCNV(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(SCNH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(SAND(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(SILT(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(CLAY(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLW(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLWX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLI(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLP(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLA(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLY(L,NY,NX),L=0,NLI(NY,NX))
!     WRITE(21,91)I,IDATA(3),(GKC4(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,IDATA(3),(GKCH(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,IDATA(3),(GKCA(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,IDATA(3),(GKCM(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,IDATA(3),(GKCN(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,IDATA(3),(GKCK(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_CEC,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_AEC,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(TCS(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(TKS(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VHCP(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VHCM(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_gasml(idg_CO2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_solml(idg_CO2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_soHml(idg_CO2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_gasml(idg_CH4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_solml(idg_CH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_soHml(idg_CH4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(ROXYF(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RCO2F(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(ROXYL(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RCH4F(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RCH4L(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_gasml(idg_H2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_solml(idg_H2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_gasml(idg_O2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_solml(idg_O2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trc_soHml(idg_O2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(ROXYX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RNH4X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RNO3X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RNO2X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RN2OX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RPO4X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RP14X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RNHBX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RN3BX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RN2BX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RPOBX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RP1BX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,95)I,IDATA(3),((ROQCX(K,L,NY,NX),L=0,NLI(NY,NX)),K=1,jcplx)
      WRITE(21,95)I,IDATA(3),((ROQAX(K,L,NY,NX),L=0,NLI(NY,NX)),K=1,jcplx)
      WRITE(21,91)I,IDATA(3),(VOLWH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLIH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(VOLAH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(RTDNT(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(ZL(L,NY,NX),L=0,JC)
      WRITE(21,91)I,IDATA(3),(ARLFT(L,NY,NX),L=1,JC)
      WRITE(21,91)I,IDATA(3),(ARSTT(L,NY,NX),L=1,JC)
      WRITE(21,91)I,IDATA(3),(WGLFT(L,NY,NX),L=1,JC)
      WRITE(21,91)I,IDATA(3),(trcs_VLN(ids_NH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcs_VLN(ids_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcs_VLN(ids_NO3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcs_VLN(ids_NO3B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcs_VLN(ids_H1PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcs_VLN(ids_H1PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_NH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_OHe,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_OH,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_OHp,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_HPO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_H2PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_OHeB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_OHB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_OHpB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_HPO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcx_solml(idx_H2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_AlPO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_FePO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_CaHPO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_HA,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_CaH2PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_AlPO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_FePO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_CaHPO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_HAB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_CaH2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(CION(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Al,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Fe,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Hp,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Ca,L,NY,NX),L=1,NLI(NY,NX))
      IF(ISALTG.NE.0)THEN
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Al,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Fe,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Hp,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Ca,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Mg,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Na,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_K,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_OH,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_SO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_Cl,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_HCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_AlOH,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_AlOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_AlOH3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_AlOH4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_AlSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeOH,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeOH3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeOH4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaHCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_MgOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_MgCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_MgHCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_MgSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_NaCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_NaSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_KSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_H0PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_H3PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeHPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_FeH2PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaHPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_CaH2PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcs_solsml(idsa_MgHPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Mg,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Na,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_OH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_SO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_Cl,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_HCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_AlOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_AlOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_AlOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_AlSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_NaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_NaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_KSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_H0PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_H3PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_H0PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_H3PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_FeH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_CaH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_solml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Al,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Fe,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Hp,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Ca,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Mg,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Na,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_OH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_SO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_Cl,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_HCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_AlOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_AlOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_AlOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_AlSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_MgOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_MgCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_MgSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_NaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_NaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_KSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_H0PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_H3PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_MgHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_H0PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_H3PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_FeH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_CaH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcsa_soHml(idsa_MgHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_Hp,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_Al,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_Fe,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_Ca,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_Mg,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_Na,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcx_solml(idx_K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,IDATA(3),(trcp_salml(idsp_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
      ENDIF
      D9985: DO K=1,jcplx
        DO  N=1,NFGs
          DO NGL=JGnio(N),JGnfo(N)
            WRITE(22,91)I,IDATA(3),(ROXYS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMX4(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMX3(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMX2(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMX1(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMB4(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMB3(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RVMB2(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RINHO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RINOO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RIPOO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RINHB(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RINOB(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(RIPBO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(ROQCS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(ROQAS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),RINHOR(NGL,K,NY,NX)
            WRITE(22,91)I,IDATA(3),RINOOR(NGL,K,NY,NX)
            WRITE(22,91)I,IDATA(3),RIPOOR(NGL,K,NY,NX)
            DO M=1,nlbiomcp
              WRITE(22,91)I,IDATA(3),(OMC(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
              WRITE(22,91)I,IDATA(3),(OMN(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
              WRITE(22,91)I,IDATA(3),(OMP(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            enddo
          enddo
        enddo
      ENDDO D9985

      DO  N=1,NFGs
        DO NGL=JGniA(N),JGnfA(N)
          WRITE(22,91)I,IDATA(3),(ROXYSff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMX4ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMX3ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMX2ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMX1ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMB4ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMB3ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RVMB2ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RINHOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RINOOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RIPOOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RINHBff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RINOBff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(RIPBOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),RINHORff(NGL,NY,NX)
          WRITE(22,91)I,IDATA(3),RINOORff(NGL,NY,NX)
          WRITE(22,91)I,IDATA(3),RIPOORff(NGL,NY,NX)
          DO M=1,nlbiomcp
            WRITE(22,91)I,IDATA(3),(OMCff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(OMNff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,IDATA(3),(OMPff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
          enddo
        enddo
      enddo


      D9980: DO K=1,jcplx
        D9975: DO M=1,ndbiomcp
          WRITE(22,91)I,IDATA(3),(ORC(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(ORN(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(ORP(M,K,L,NY,NX),L=0,NLI(NY,NX))
        ENDDO D9975
        WRITE(22,91)I,IDATA(3),(OQC(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQN(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQP(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQA(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQCH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQNH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQPH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OQAH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OHC(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OHN(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OHP(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,IDATA(3),(OHA(K,L,NY,NX),L=0,NLI(NY,NX))
        D9970: DO M=1,jsken
          WRITE(22,91)I,IDATA(3),(OSC(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(OSA(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(OSN(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,IDATA(3),(OSP(M,K,L,NY,NX),L=0,NLI(NY,NX))
        ENDDO D9970
      ENDDO D9980
      WRITE(22,91)I,IDATA(3),(RVMXC(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ORGC(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ORGR(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNH4FA(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNH3FA(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNHUFA(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNO3FA(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNH4FB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNH3FB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNHUFB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNO3FB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_gasml(idg_N2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(idg_N2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(idg_N2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_gasml(idg_N2O,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(idg_N2O,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(idg_N2O,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_gasml(idg_NH3,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_NH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_NH4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(idg_NH3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(idg_NH3,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_NO3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_NO3,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_NO2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_NO2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_H1PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_H2PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_H1PO4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_H2PO4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_NH4B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(idg_NH3B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(idg_NH3B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_NO3B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_NO3B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_NO2B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_NO2B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_H1PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_solml(ids_H2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_H1PO4B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(trc_soHml(ids_H2PO4B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(WDNHB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(DPNHB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(WDNOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(DPNOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(WDPOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(DPPOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNHUI(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNHU0(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNFNI(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(ZNFN0(L,NY,NX),L=0,NLI(NY,NX))
9990  CONTINUE
9995  CONTINUE
90    FORMAT(2I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
91    FORMAT(2I4,21E17.8E3)
!92    FORMAT(2I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
!      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
93    FORMAT(8I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
95    FORMAT(2I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3 &
      ,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
  RETURN
  END SUBROUTINE wouts
