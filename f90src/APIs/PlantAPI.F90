module PlantAPI

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ExtractsMod   , only : extracts
  use grosubsMod    , only : grosubs
  use HfuncsMod    , only : hfuncs
  use EcoSIMSolverPar
  use UptakesMod    , only : uptakes
  use EcoSiMParDataMod, only : micpar, pltpar
  use PlantDisturbMod, only : PrepLandscapeGrazing
  use timings      , only : start_timer, end_timer
  use EcoSIMHistMod
  use SnowDataType
  use TracerIDMod
  use SoilPhysDataType, only : ALBX
  use SurfLitterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use EcoSimSumDataType
  use SoilHeatDataType
  use SOMDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use GridDataType
  use RootDataType
  use SoilWaterDataType
  use CanopyDataType
  use PlantDataRateType
  use PlantTraitDataType
  use CanopyRadDataType
  use FlagDataType
  use EcosimBGCFluxType
  use FertilizerDataType
  use SoilBGCDataType
  use PlantMngmtDataType
  use PlantAPIData
implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: PlantModel
  public :: PlantCanopyRadsModel

  contains


  subroutine PlantModel(I,J,NHW,NHE,NVN,NVS)


  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1
  integer :: NY,NX

333   FORMAT(A8)

  call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)

  plt_site%TBALE(:)=TBALE(:)
  DO NX=NHW,NHE
    DO NY=NVN,NVS

    call  PlantAPISend(I,J,NY,NX)

!   UPDATE PLANT PHENOLOGY IN 'HFUNC'
!
    !if(lverb)WRITE(*,333)'HFUNC'
    !call start_timer(t1)

      CALL HFUNCs(I,J)
    !call end_timer('HFUNC',t1)

      CALL UPTAKES(I,J)

      CALL GROSUBs(I,J)

      CALL EXTRACTs(I,J)

      call PlantAPIRecv(I,J,NY,NX)
    ENDDO
  ENDDO
  TBALE(:)=plt_site%TBALE(:)


  end subroutine PlantModel
!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DPTH0)
  use CanopyCondsMod
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DPTH0

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DPTH0)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel
!------------------------------------------------------------------------------------------

  subroutine PlantAPIRecv(I,J,NY,NX)
  !
  !DESCRIPTION
  !
  use EcoSIMConfig, only : jsken=>jskenc,jcplx => jcplxc
  use PlantAPIData, only : plt_rad
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: NB,NR,NZ,K,L,M,N,I1,NE

  I1=I+1;if(I1>LYRC)I1=1
  IFLGT(NY,NX)=plt_site%IFLGT
  PPT(NY,NX) =plt_site%PPT
  RECO(NY,NX)=plt_bgcr%RECO
  TNBP(NY,NX)=plt_bgcr%TNBP
  TGPP(NY,NX)=plt_bgcr%TGPP
  TLEC(NY,NX)=plt_ew%TLEC
  TSHC(NY,NX)=plt_ew%TSHC
  TRAU(NY,NX)=plt_bgcr%TRAU
  ZESNC(1:npelms,NY,NX) =plt_bgcr%ZESNC(1:npelms)
  XHVSTE(1:npelms,NY,NX)=plt_distb%XHVSTE(1:npelms)
  TVOLWC(NY,NX)=plt_ew%TVOLWC
  TSH(NY,NX)   =plt_ew%TSH
  WTSTGT(NY,NX)=plt_biom%WTSTGT
  UVOLO(NY,NX) =plt_ew%UVOLO
  ARSTC(NY,NX) =plt_morph%ARSTC
  ARLFC(NY,NX) =plt_morph%ARLFC
  TRN(NY,NX)   =plt_rad%TRN
  TLE(NY,NX)   =plt_ew%TLE
  TGH(NY,NX)   =plt_ew%TGH
  TEVAPP(NY,NX)=plt_ew%TEVAPP
  THRMC(NY,NX) =plt_ew%THRMC
  TEVAPC(NY,NX)=plt_ew%TEVAPC
  THFLXC(NY,NX)=plt_ew%THFLXC
  TVOLWP(NY,NX)=plt_ew%TVOLWP
  TENGYC(NY,NX)=plt_ew%TENGYC
  TH2GZ(NY,NX) =plt_bgcr%TH2GZ
  TN2OZ(NY,NX) =plt_rbgc%TN2OZ
  TCO2Z(NY,NX) =plt_rbgc%TCO2Z

  TNH3Z(NY,NX)=plt_rbgc%TNH3Z
  TOXYZ(NY,NX)=plt_rbgc%TOXYZ

  TCCAN(NY,NX)=plt_bgcr%TCCAN
  TCH4Z(NY,NX)=plt_rbgc%TCH4Z

  FERT(17:19,I1,NY,NX)=plt_distb%FERT(17:19)
  FERT(3,I1,NY,NX) =plt_distb%FERT(3)
  IYTYP(2,I1,NY,NX)=plt_distb%IYTYP
  FWOODE(1:npelms,1:n_pltlitrk) =plt_allom%FWOODE(1:npelms,1:n_pltlitrk)
  FWODRE(1:npelms,1:n_pltlitrk) =plt_allom%FWODRE(1:npelms,1:n_pltlitrk)
  FWODBE(1:npelms,1:n_pltlitrk) =plt_allom%FWODBE(1:npelms,1:n_pltlitrk)
  FWODLE(1:npelms,1:n_pltlitrk)=plt_allom%FWODLE(1:npelms,1:n_pltlitrk)
  VOLWOU   =plt_site%VOLWOU
  DO L=1,JC
    WGLFT(L,NY,NX)=plt_biom%WGLFT(L)
    ARSTT(L,NY,NX)=plt_morph%ARSTT(L)
    ARLFT(L,NY,NX)=plt_morph%ARLFT(L)
  ENDDO
  DO L=NU(NY,NX),NL(NY,NX)
    DO K=1,jcplx
      XOQCS(K,L,NY,NX)=plt_bgcr%XOQCS(K,L)
      XOQNS(K,L,NY,NX)=plt_bgcr%XOQNS(K,L)
      XOQPS(K,L,NY,NX)=plt_bgcr%XOQPS(K,L)
    ENDDO
    DO M=1,NPH
      ROXSK(M,L,NY,NX)=plt_rbgc%ROXSK(M,L)
    ENDDO
  ENDDO

  DO L=0,JZ
    RPOBX(L,NY,NX) =plt_bgcr%RPOBX(L)
    RP1BX(L,NY,NX) =plt_bgcr%RP1BX(L)
    RN3BX(L,NY,NX) =plt_bgcr%RN3BX(L)
    RNHBX(L,NY,NX) =plt_bgcr%RNHBX(L)
    RP14X(L,NY,NX) =plt_bgcr%RP14X(L)
    RPO4X(L,NY,NX) =plt_bgcr%RPO4X(L)
    RNO3X(L,NY,NX) =plt_bgcr%RNO3X(L)
    RNH4X(L,NY,NX) =plt_bgcr%RNH4X(L)
    ROXYX(L,NY,NX) =plt_bgcr%ROXYX(L)
    TUPHT(L,NY,NX) =plt_ew%TUPHT(L)
    TUPWTR(L,NY,NX)=plt_ew%TUPWTR(L)
    DO  K=1,micpar%n_pltlitrk
      DO NE=1,npelms
        DO  M=1,jsken
          ESNT(M,NE,K,L,NY,NX)=plt_bgcr%ESNT(M,NE,K,L)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,JZ
    TUPNF(L,NY,NX) =plt_rbgc%TUPNF(L)
    RTDNT(L,NY,NX) =plt_morph%RTDNT(L)
    TLOXYP(L,NY,NX)=plt_rbgc%TLOXYP(L)
    TLCO2P(L,NY,NX)=plt_rbgc%TLCO2P(L)
    TLCH4P(L,NY,NX)=plt_rbgc%TLCH4P(L)
    TLN2OP(L,NY,NX)=plt_rbgc%TLN2OP(L)
    TLNH3P(L,NY,NX)=plt_rbgc%TLNH3P(L)
    TLH2GP(L,NY,NX)=plt_bgcr%TLH2GP(L)
    TCOFLA(L,NY,NX)=plt_rbgc%TCOFLA(L)
    TOXFLA(L,NY,NX)=plt_rbgc%TOXFLA(L)
    TCHFLA(L,NY,NX)=plt_rbgc%TCHFLA(L)
    TN2FLA(L,NY,NX)=plt_rbgc%TN2FLA(L)
    TNHFLA(L,NY,NX)=plt_rbgc%TNHFLA(L)
    THGFLA(L,NY,NX)=plt_bgcr%THGFLA(L)
    TCO2P(L,NY,NX) =plt_bgcr%TCO2P(L)
    TUPOXP(L,NY,NX)=plt_bgcr%TUPOXP(L)
    TCO2S(L,NY,NX) =plt_bgcr%TCO2S(L)
    TUPOXS(L,NY,NX)=plt_bgcr%TUPOXS(L)
    TUPCHS(L,NY,NX)=plt_bgcr%TUPCHS(L)
    TUPN2S(L,NY,NX)=plt_bgcr%TUPN2S(L)
    TUPN3S(L,NY,NX)=plt_bgcr%TUPN3S(L)
    TUPN3B(L,NY,NX)=plt_bgcr%TUPN3B(L)
    TUPHGS(L,NY,NX)=plt_bgcr%TUPHGS(L)
    TUPNH4(L,NY,NX)=plt_bgcr%TUPNH4(L)
    TUPNO3(L,NY,NX)=plt_bgcr%TUPNO3(L)
    TUPH1B(L,NY,NX)=plt_bgcr%TUPH1B(L)
    TUPH2B(L,NY,NX)=plt_bgcr%TUPH2B(L)
    TUPNHB(L,NY,NX)=plt_bgcr%TUPNHB(L)
    TUPNOB(L,NY,NX)=plt_bgcr%TUPNOB(L)
    TUPH1P(L,NY,NX)=plt_bgcr%TUPH1P(L)
    TUPH2P(L,NY,NX)=plt_bgcr%TUPH2P(L)
    DO  K=1,jcplx
      TDFOMP(K,L,NY,NX)=plt_bgcr%TDFOMP(K,L)
      TDFOMN(K,L,NY,NX)=plt_bgcr%TDFOMN(K,L)
      TDFOMC(K,L,NY,NX)=plt_bgcr%TDFOMC(K,L)
    ENDDO
  ENDDO
  DO NZ=1,NP0(NY,NX)
    WTRTE(1:npelms,NZ,NY,NX) =plt_biom%WTRTE(1:npelms,NZ)
    BALE(1:npelms,NZ,NY,NX)  =plt_site%BALE(1:npelms,NZ)
    EPOOLP(1:npelms,NZ,NY,NX)=plt_biom%EPOOLP(1:npelms,NZ)
    EPOLNP(1:npelms,NZ,NY,NX)=plt_biom%EPOLNP(1:npelms,NZ)
    CEPOLP(1:npelms,NZ,NY,NX)=plt_biom%CEPOLP(1:npelms,NZ)
    HESNC(1:npelms,NZ,NY,NX) =plt_bgcr%HESNC(1:npelms,NZ)
    HVSTE(1:npelms,NZ,NY,NX) =plt_distb%HVSTE(1:npelms,NZ)
    RSETE(1:npelms,NZ,NY,NX) =plt_pheno%RSETE(1:npelms,NZ)
    TESN0(1:npelms,NZ,NY,NX) =plt_bgcr%TESN0(1:npelms,NZ)
    TESNC(1:npelms,NZ,NY,NX) =plt_bgcr%TESNC(1:npelms,NZ)
    THVSTE(1:npelms,NZ,NY,NX) =plt_distb%THVSTE(1:npelms,NZ)
    TEUPTK(1:npelms,NZ,NY,NX) =plt_rbgc%TEUPTK(1:npelms,NZ)
    UPOME(1:npelms,NZ,NY,NX)  =plt_rbgc%UPOME(1:npelms,NZ)
    WTSTGE(1:npelms,NZ,NY,NX)  =plt_biom%WTSTGE(1:npelms,NZ)
    WTRVE(1:npelms,NZ,NY,NX)  =plt_biom%WTRVE(1:npelms,NZ)
    WTSHTE(1:npelms,NZ,NY,NX)  =plt_biom%WTSHTE(1:npelms,NZ)
    WTLFE(1:npelms,NZ,NY,NX)   =plt_biom%WTLFE(1:npelms,NZ)
    WTSHEE(1:npelms,NZ,NY,NX)  =plt_biom%WTSHEE(1:npelms,NZ)
    WTSTKE(1:npelms,NZ,NY,NX)  =plt_biom%WTSTKE(1:npelms,NZ)
    WTRSVE(1:npelms,NZ,NY,NX)  =plt_biom%WTRSVE(1:npelms,NZ)
    WTHSKE(1:npelms,NZ,NY,NX)  =plt_biom%WTHSKE(1:npelms,NZ)
    WTEARE(1:npelms,NZ,NY,NX)  =plt_biom%WTEARE(1:npelms,NZ)
    WTGRE(1:npelms,NZ,NY,NX)   =plt_biom%WTGRE(1:npelms,NZ)
    WTRTSE(1:npelms,NZ,NY,NX)  =plt_biom%WTRTSE(1:npelms,NZ)
    WTNDE(1:npelms,NZ,NY,NX)   =plt_biom%WTNDE(1:npelms,NZ)
    HEUPTK(1:npelms,NZ,NY,NX)=plt_rbgc%HEUPTK(1:npelms,NZ)
    ARLFP(NZ,NY,NX) =plt_morph%ARLFP(NZ)
    ARSTP(NZ,NY,NX) =plt_morph%ARSTP(NZ)
    CCPLNP(NZ,NY,NX)=plt_biom%CCPLNP(NZ)
    CNET(NZ,NY,NX)  =plt_bgcr%CNET(NZ)
    CO2Q(NZ,NY,NX)  =plt_photo%CO2Q(NZ)
    CO2I(NZ,NY,NX)  =plt_photo%CO2I(NZ)
    CO2L(NZ,NY,NX)  =plt_photo%CO2L(NZ)
    CF(NZ,NY,NX)    =plt_morph%CF(NZ)
    CNWS(NZ,NY,NX)  =plt_allom%CNWS(NZ)
    CPWS(NZ,NY,NX)  =plt_allom%CPWS(NZ)
    CWSRT(NZ,NY,NX) =plt_allom%CWSRT(NZ)
    CNRTS(NZ,NY,NX) =plt_allom%CNRTS(NZ)
    CPRTS(NZ,NY,NX) =plt_allom%CPRTS(NZ)
    CTRAN(NZ,NY,NX) =plt_ew%CTRAN(NZ)
    CARBN(NZ,NY,NX) =plt_bgcr%CARBN(NZ)
    CHILL(NZ,NY,NX) =plt_photo%CHILL(NZ)
    DCO2(NZ,NY,NX)  =plt_photo%DCO2(NZ)
    DTKC(NZ,NY,NX)  =plt_ew%DTKC(NZ)
    ENGYX(NZ,NY,NX) =plt_ew%ENGYX(NZ)
    EP(NZ,NY,NX)    =plt_ew%EP(NZ)
    EVAPC(NZ,NY,NX) =plt_ew%EVAPC(NZ)
    EFLXC(NZ,NY,NX) =plt_ew%EFLXC(NZ)
    FMOL(NZ,NY,NX)  =plt_photo%FMOL(NZ)
    FRADP(NZ,NY,NX) =plt_rad%FRADP(NZ)
    FNOD(NZ,NY,NX)  =plt_allom%FNOD(NZ)
    GRNO(NZ,NY,NX)  =plt_morph%GRNO(NZ)
    HTCTL(NZ,NY,NX) =plt_morph%HTCTL(NZ)
    HTC(NZ,NY,NX)   =plt_pheno%HTC(NZ)
    HFLXC(NZ,NY,NX) =plt_ew%HFLXC(NZ)
    HTSTZ(NZ,NY,NX) =plt_morph%HTSTZ(NZ)
    IFLGC(NZ,NY,NX) =plt_pheno%IFLGC(NZ)
    IDTH(NZ,NY,NX)  =plt_pheno%IDTH(NZ)
    IDTHP(NZ,NY,NX) =plt_pheno%IDTHP(NZ)
    IDTHR(NZ,NY,NX) =plt_pheno%IDTHR(NZ)
    IDAY0(NZ,NY,NX) =plt_distb%IDAY0(NZ)
    IYR0(NZ,NY,NX)  =plt_distb%IYR0(NZ)
    IFLGI(NZ,NY,NX) =plt_pheno%IFLGI(NZ)
    IDAYH(NZ,NY,NX) =plt_distb%IDAYH(NZ)
    IYRH(NZ,NY,NX)  =plt_distb%IYRH(NZ)
    NIX(NZ,NY,NX)   =plt_morph%NIX(NZ)
    NBT(NZ,NY,NX)   =plt_morph%NBT(NZ)
    NBR(NZ,NY,NX)   =plt_morph%NBR(NZ)
    NRT(NZ,NY,NX)   =plt_morph%NRT(NZ)
    NI(NZ,NY,NX)    =plt_morph%NI(NZ)
    NG(NZ,NY,NX)    =plt_morph%NG(NZ)
    NB1(NZ,NY,NX)   =plt_morph%NB1(NZ)
    NNOD(NZ,NY,NX)  =plt_morph%NNOD(NZ)
    O2L(NZ,NY,NX)   =plt_photo%O2L(NZ)
    O2I(NZ,NY,NX)   =plt_photo%O2I(NZ)
    OFFST(NZ,NY,NX) =plt_pheno%OFFST(NZ)
    OSTR(NZ,NY,NX)  =plt_pheno%OSTR(NZ)
    PPX(NZ,NY,NX)   =plt_site%PPX(NZ)
    PP(NZ,NY,NX)    =plt_site%PP(NZ)
    PPI(NZ,NY,NX)   =plt_site%PPI(NZ)
    PSILT(NZ,NY,NX) =plt_ew%PSILT(NZ)
    PSILO(NZ,NY,NX) =plt_ew%PSILO(NZ)
    PSILG(NZ,NY,NX) =plt_ew%PSILG(NZ)
    PSILZ(NZ,NY,NX) =plt_ew%PSILZ(NZ)
    RCO2Z(NZ,NY,NX) =plt_bgcr%RCO2Z(NZ)
    ROXYZ(NZ,NY,NX) =plt_bgcr%ROXYZ(NZ)
    RCH4Z(NZ,NY,NX) =plt_bgcr%RCH4Z(NZ)
    RN2OZ(NZ,NY,NX) =plt_bgcr%RN2OZ(NZ)
    RNH3Z(NZ,NY,NX) =plt_bgcr%RNH3Z(NZ)
    RH2GZ(NZ,NY,NX) =plt_bgcr%RH2GZ(NZ)
    RSMN(NZ,NY,NX)  =plt_photo%RSMN(NZ)
    RSMH(NZ,NY,NX)  =plt_photo%RSMH(NZ)
    RCMX(NZ,NY,NX)  =plt_photo%RCMX(NZ)
    RNH3C(NZ,NY,NX) =plt_bgcr%RNH3C(NZ)
    RAD1(NZ,NY,NX)  =plt_rad%RAD1(NZ)
    RAZ(NZ,NY,NX)   =plt_ew%RAZ(NZ)
    RC(NZ,NY,NX)    =plt_photo%RC(NZ)
    RA(NZ,NY,NX)    =plt_photo%RA(NZ)
    SDPTHI(NZ,NY,NX)=plt_morph%SDPTHI(NZ)
    SCO2(NZ,NY,NX)  =plt_photo%SCO2(NZ)
    SO2(NZ,NY,NX)   =plt_photo%SO2(NZ)
    SSTX(NZ,NY,NX)  =plt_pheno%SSTX(NZ)
    SDVL(NZ,NY,NX)  =plt_morph%SDVL(NZ)
    SDLG(NZ,NY,NX)  =plt_morph%SDLG(NZ)
    SDAR(NZ,NY,NX)  =plt_morph%SDAR(NZ)
    SDPTH(NZ,NY,NX) =plt_morph%SDPTH(NZ)
    SFLXC(NZ,NY,NX) =plt_ew%SFLXC(NZ)
    TCO2T(NZ,NY,NX) =plt_bgcr%TCO2T(NZ)
    TCO2A(NZ,NY,NX) =plt_bgcr%TCO2A(NZ)
    TCZ(NZ,NY,NX)   =plt_pheno%TCZ(NZ)
    TCX(NZ,NY,NX)   =plt_pheno%TCX(NZ)
    TNH3C(NZ,NY,NX)  =plt_bgcr%TNH3C(NZ)
    TZUPFX(NZ,NY,NX) =plt_bgcr%TZUPFX(NZ)
    TKC(NZ,NY,NX)    =plt_ew%TKC(NZ)
    TCC(NZ,NY,NX)    =plt_ew%TCC(NZ)
    THRM1(NZ,NY,NX)  =plt_rad%THRM1(NZ)
    TKCZ(NZ,NY,NX)   =plt_ew%TKCZ(NZ)
    TKG(NZ,NY,NX)    =plt_pheno%TKG(NZ)
    TCG(NZ,NY,NX)    =plt_pheno%TCG(NZ)
    TFN3(NZ,NY,NX)   =plt_pheno%TFN3(NZ)
    UPNF(NZ,NY,NX)   =plt_rbgc%UPNF(NZ)
    UPNH4(NZ,NY,NX)  =plt_rbgc%UPNH4(NZ)
    UPNO3(NZ,NY,NX)  =plt_rbgc%UPNO3(NZ)
    UPH2P(NZ,NY,NX)  =plt_rbgc%UPH2P(NZ)
    UPH1P(NZ,NY,NX)  =plt_rbgc%UPH1P(NZ)
    VHCPC(NZ,NY,NX)  =plt_ew%VHCPC(NZ)
    VOLWP(NZ,NY,NX)  =plt_ew%VOLWP(NZ)
    VCO2F(NZ,NY,NX)  =plt_distb%VCO2F(NZ)
    VCH4F(NZ,NY,NX)  =plt_distb%VCH4F(NZ)
    VOXYF(NZ,NY,NX)  =plt_distb%VOXYF(NZ)
    VNH3F(NZ,NY,NX)  =plt_distb%VNH3F(NZ)
    VN2OF(NZ,NY,NX)  =plt_distb%VN2OF(NZ)
    VPO4F(NZ,NY,NX)  =plt_distb%VPO4F(NZ)
    VOLWC(NZ,NY,NX)  =plt_ew%VOLWC(NZ)
    WSTR(NZ,NY,NX)   =plt_pheno%WSTR(NZ)
    WTRVX(NZ,NY,NX)  =plt_biom%WTRVX(NZ)
    WVSTK(NZ,NY,NX)  =plt_biom%WVSTK(NZ)
    WTLS(NZ,NY,NX)   =plt_biom%WTLS(NZ)
    WTRTA(NZ,NY,NX)  =plt_biom%WTRTA(NZ)
    XKCO2L(NZ,NY,NX) =plt_photo%XKCO2L(NZ)
    XKCO2O(NZ,NY,NX) =plt_photo%XKCO2O(NZ)
    ZC(NZ,NY,NX)     =plt_morph%ZC(NZ)
    ZNPP(NZ,NY,NX)   =plt_bgcr%ZNPP(NZ)
    ZEROP(NZ,NY,NX)  =plt_biom%ZEROP(NZ)
    ZEROQ(NZ,NY,NX)  =plt_rbgc%ZEROQ(NZ)
    ZEROL(NZ,NY,NX)  =plt_biom%ZEROL(NZ)
    ZTYP(NZ,NY,NX)   =plt_pheno%ZTYP(NZ)
    HVST(NZ,I,NY,NX) =plt_distb%HVST(NZ)
    IHVST(NZ,I,NY,NX)=plt_distb%IHVST(NZ)
    JHVST(NZ,I,NY,NX)=plt_distb%JHVST(NZ)
    THIN(NZ,I,NY,NX) =plt_distb%THIN(NZ)
    DO L=1,JZ
      WTNDLE(L,1:npelms,NZ,NY,NX) =plt_biom%WTNDLE(L,1:npelms,NZ)
      EPOOLN(L,1:npelms,NZ,NY,NX)=plt_biom%EPOOLN(L,1:npelms,NZ)
      RUPNF(L,NZ,NY,NX) =plt_bgcr%RUPNF(L,NZ)
      TFN4(L,NZ,NY,NX)  =plt_pheno%TFN4(L,NZ)
    ENDDO
    DO L=1,JC
      ARLFV(L,NZ,NY,NX)=plt_morph%ARLFV(L,NZ)
      WGLFV(L,NZ,NY,NX)=plt_biom%WGLFV(L,NZ)
      ARSTV(L,NZ,NY,NX)=plt_morph%ARSTV(L,NZ)
    ENDDO

    DO L=0,JZ
      DO K=1,micpar%n_pltlitrk
        DO M=1,jsken
          ESNC(M,1:npelms,K,L,NZ,NY,NX)=plt_bgcr%ESNC(M,1:npelms,K,L,NZ)
        ENDDO
      ENDDO
    ENDDO

    DO NE=1,npelms
      DO NB=1,NBR(NZ,NY,NX)
        EPOOL(NB,NE,NZ,NY,NX) =plt_biom%EPOOL(NB,NE,NZ)
        WTEARBE(NB,NE,NZ,NY,NX)=plt_biom%WTEARBE(NB,NE,NZ)
      ENDDO
    ENDDO


    DO NE=1,npelms
      DO NB=1,NBR(NZ,NY,NX)
        EPOLNB(NB,NE,NZ,NY,NX)=plt_biom%EPOLNB(NB,NE,NZ)
        WTSHTBE(NB,NE,NZ,NY,NX)=plt_biom%WTSHTBE(NB,NE,NZ)
        WTSHEBE(NB,NE,NZ,NY,NX)=plt_biom%WTSHEBE(NB,NE,NZ)
        WTSTKBE(NB,NE,NZ,NY,NX)=plt_biom%WTSTKBE(NB,NE,NZ)
        CEPOLB(NB,NE,NZ,NY,NX)=plt_biom%CEPOLB(NB,NE,NZ)
        WTLFBE(NB,NE,NZ,NY,NX) =plt_biom%WTLFBE(NB,NE,NZ)
        WTRSVBE(NB,NE,NZ,NY,NX)=plt_biom%WTRSVBE(NB,NE,NZ)
        WTHSKBE(NB,NE,NZ,NY,NX)=plt_biom%WTHSKBE(NB,NE,NZ)
        WTGRBE(NB,NE,NZ,NY,NX) =plt_biom%WTGRBE(NB,NE,NZ)
        WTNDBE(NB,NE,NZ,NY,NX) =plt_biom%WTNDBE(NB,NE,NZ)
        WGSHEXE(NB,NE,NZ,NY,NX)=plt_biom%WGSHEXE(NB,NE,NZ)
      ENDDO
    ENDDO
    DO NB=1,NBR(NZ,NY,NX)
      DO L=1,JC
        ARSTK(L,NB,NZ,NY,NX)=plt_morph%ARSTK(L,NB,NZ)
      ENDDO
      ATRP(NB,NZ,NY,NX)  =plt_pheno%ATRP(NB,NZ)
      ARLFB(NB,NZ,NY,NX) =plt_morph%ARLFB(NB,NZ)
      ARLFZ(NB,NZ,NY,NX) =plt_morph%ARLFZ(NB,NZ)

      DGSTGI(NB,NZ,NY,NX)=plt_pheno%DGSTGI(NB,NZ)
      DGSTGF(NB,NZ,NY,NX)=plt_pheno%DGSTGF(NB,NZ)
      FLG4(NB,NZ,NY,NX)  =plt_pheno%FLG4(NB,NZ)
      FDBK(NB,NZ,NY,NX)  =plt_photo%FDBK(NB,NZ)
      FDBKX(NB,NZ,NY,NX) =plt_photo%FDBKX(NB,NZ)
      FLGZ(NB,NZ,NY,NX)  =plt_pheno%FLGZ(NB,NZ)
      GROUP(NB,NZ,NY,NX) =plt_pheno%GROUP(NB,NZ)
      GSTGI(NB,NZ,NY,NX) =plt_pheno%GSTGI(NB,NZ)
      GSTGF(NB,NZ,NY,NX) =plt_pheno%GSTGF(NB,NZ)
      GRNXB(NB,NZ,NY,NX) =plt_morph%GRNXB(NB,NZ)
      GRNOB(NB,NZ,NY,NX) =plt_morph%GRNOB(NB,NZ)
      GRWTB(NB,NZ,NY,NX) =plt_allom%GRWTB(NB,NZ)
      HTSHEX(NB,NZ,NY,NX)=plt_morph%HTSHEX(NB,NZ)
      IFLGP(NB,NZ,NY,NX) =plt_pheno%IFLGP(NB,NZ)
      IDTHB(NB,NZ,NY,NX) =plt_pheno%IDTHB(NB,NZ)
      IFLGF(NB,NZ,NY,NX) =plt_pheno%IFLGF(NB,NZ)
      IFLGE(NB,NZ,NY,NX) =plt_pheno%IFLGE(NB,NZ)
      IFLGA(NB,NZ,NY,NX) =plt_pheno%IFLGA(NB,NZ)
      IFLGG(NB,NZ,NY,NX) =plt_pheno%IFLGG(NB,NZ)
      IFLGR(NB,NZ,NY,NX) =plt_pheno%IFLGR(NB,NZ)
      IFLGQ(NB,NZ,NY,NX) =plt_pheno%IFLGQ(NB,NZ)
      KVSTG(NB,NZ,NY,NX) =plt_pheno%KVSTG(NB,NZ)
      KLEAF(NB,NZ,NY,NX) =plt_morph%KLEAF(NB,NZ)
      KLEAFX(NB,NZ,NY,NX)=plt_morph%KLEAFX(NB,NZ)
      KVSTGN(NB,NZ,NY,NX)=plt_pheno%KVSTGN(NB,NZ)
      NBTB(NB,NZ,NY,NX)  =plt_morph%NBTB(NB,NZ)
      PSTG(NB,NZ,NY,NX)  =plt_morph%PSTG(NB,NZ)
      PSTGI(NB,NZ,NY,NX) =plt_morph%PSTGI(NB,NZ)
      PSTGF(NB,NZ,NY,NX) =plt_morph%PSTGF(NB,NZ)
      RCELX(NB,1:npelms,NZ,NY,NX) =plt_pheno%RCELX(NB,1:npelms,NZ)
      RCESX(NB,1:npelms,NZ,NY,NX) =plt_pheno%RCESX(NB,1:npelms,NZ)
      RNH3B(NB,NZ,NY,NX) =plt_rbgc%RNH3B(NB,NZ)
      TGSTGI(NB,NZ,NY,NX)=plt_pheno%TGSTGI(NB,NZ)
      TGSTGF(NB,NZ,NY,NX)=plt_pheno%TGSTGF(NB,NZ)
      VRNY(NB,NZ,NY,NX)  =plt_pheno%VRNY(NB,NZ)
      VRNZ(NB,NZ,NY,NX)  =plt_pheno%VRNZ(NB,NZ)
      VRNS(NB,NZ,NY,NX)  =plt_pheno%VRNS(NB,NZ)
      VRNF(NB,NZ,NY,NX)  =plt_pheno%VRNF(NB,NZ)
      VSTG(NB,NZ,NY,NX)  =plt_morph%VSTG(NB,NZ)
      VSTGX(NB,NZ,NY,NX) =plt_pheno%VSTGX(NB,NZ)
      WTLSB(NB,NZ,NY,NX) =plt_biom%WTLSB(NB,NZ)

      WGLFEX(NB,1:npelms,NZ,NY,NX) =plt_biom%WGLFEX(NB,1:npelms,NZ)
      WTSTXBE(NB,1:npelms,NZ,NY,NX)=plt_biom%WTSTXBE(NB,1:npelms,NZ)
      WVSTKB(NB,NZ,NY,NX)=plt_biom%WVSTKB(NB,NZ)

      DO K=0,JNODS
        ARLF(K,NB,NZ,NY,NX)=plt_morph%ARLF1(K,NB,NZ)
        HTNODX(K,NB,NZ,NY,NX)=plt_morph%HTNODX(K,NB,NZ)
        HTNODE(K,NB,NZ,NY,NX)=plt_morph%HTNODE(K,NB,NZ)
        HTSHE(K,NB,NZ,NY,NX) =plt_morph%HTSHE(K,NB,NZ)
        WGNODE(K,NB,1:npelms,NZ,NY,NX)=plt_biom%WGNODE(K,NB,1:npelms,NZ)
        WGLFE(K,NB,1:npelms,NZ,NY,NX)  =plt_biom%WGLFE(K,NB,1:npelms,NZ)
        WSLF(K,NB,NZ,NY,NX)  =plt_biom%WSLF(K,NB,NZ)
        WGSHE(K,NB,1:npelms,NZ,NY,NX) =plt_biom%WGSHE(K,NB,1:npelms,NZ)
        WSSHE(K,NB,NZ,NY,NX) =plt_biom%WSSHE(K,NB,NZ)
      ENDDO
      DO  L=1,JC
        DO N=1,JLI
          SURFB(N,L,NB,NZ,NY,NX)=plt_morph%SURFB(N,L,NB,NZ)
        ENDDO
      ENDDO
      DO K=0,JNODS
        DO  L=1,JC
          ARLFL(L,K,NB,NZ,NY,NX) =plt_morph%ARLFL(L,K,NB,NZ)
          WGLFLE(L,K,NB,1:npelms,NZ,NY,NX) =plt_biom%WGLFLE(L,K,NB,1:npelms,NZ)
        ENDDO
      ENDDO
      DO M=1,pltpar%jpstgs
        IDAY(M,NB,NZ,NY,NX)=plt_pheno%IDAY(M,NB,NZ)
      ENDDO
      DO K=1,JNODS
        DO  L=1,JC
          DO N=1,JLI
            SURF(N,L,K,NB,NZ,NY,NX) =plt_morph%SURF(N,L,K,NB,NZ)
            SURFX(N,L,K,NB,NZ,NY,NX)=plt_photo%SURFX(N,L,K,NB,NZ)
          ENDDO
        ENDDO
        CPOOL3(K,NB,NZ,NY,NX)=plt_photo%CPOOL3(K,NB,NZ)
        CPOOL4(K,NB,NZ,NY,NX)=plt_photo%CPOOL4(K,NB,NZ)
        CO2B(K,NB,NZ,NY,NX)  =plt_photo%CO2B(K,NB,NZ)
        COMPL(K,NB,NZ,NY,NX) =plt_photo%COMPL(K,NB,NZ)
        CBXN(K,NB,NZ,NY,NX)  =plt_photo%CBXN(K,NB,NZ)
        CBXN4(K,NB,NZ,NY,NX) =plt_photo%CBXN4(K,NB,NZ)
        ETGRO(K,NB,NZ,NY,NX) =plt_photo%ETGRO(K,NB,NZ)
        ETGR4(K,NB,NZ,NY,NX) =plt_photo%ETGR4(K,NB,NZ)
        FDBK4(K,NB,NZ,NY,NX) =plt_photo%FDBK4(K,NB,NZ)
        HCOB(K,NB,NZ,NY,NX)  =plt_photo%HCOB(K,NB,NZ)
        VCGRO(K,NB,NZ,NY,NX) =plt_photo%VCGRO(K,NB,NZ)
        VGRO(K,NB,NZ,NY,NX)  =plt_photo%VGRO(K,NB,NZ)
        VCGR4(K,NB,NZ,NY,NX) =plt_photo%VCGR4(K,NB,NZ)
        VGRO4(K,NB,NZ,NY,NX) =plt_photo%VGRO4(K,NB,NZ)
      ENDDO
    ENDDO
    DO M=1,jsken
      WTSTDE(M,1:npelms,NZ,NY,NX)=plt_biom%WTSTDE(M,1:npelms,NZ)
    ENDDO

    DO  L=1,JZ
      DO N=1,pltpar%jroots
        EPOOLR(1:npelms,N,L,NZ,NY,NX)=plt_biom%EPOOLR(1:npelms,N,L,NZ)
        CCPOLR(N,L,NZ,NY,NX)=plt_biom%CCPOLR(N,L,NZ)
        CZPOLR(N,L,NZ,NY,NX)=plt_biom%CZPOLR(N,L,NZ)
        CPPOLR(N,L,NZ,NY,NX)=plt_biom%CPPOLR(N,L,NZ)
        CWSRTL(N,L,NZ,NY,NX)=plt_biom%CWSRTL(N,L,NZ)
        trcg_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)  =plt_rbgc%trcg_rootml(idg_beg:idg_end-1,N,L,NZ)
        trcs_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)  =plt_rbgc%trcs_rootml(idg_beg:idg_end-1,N,L,NZ)
        PSIRT(N,L,NZ,NY,NX) =plt_ew%PSIRT(N,L,NZ)
        PSIRO(N,L,NZ,NY,NX) =plt_ew%PSIRO(N,L,NZ)
        PSIRG(N,L,NZ,NY,NX) =plt_ew%PSIRG(N,L,NZ)
        RTN1(N,L,NZ,NY,NX)  =plt_morph%RTN1(N,L,NZ)
        RTNL(N,L,NZ,NY,NX)  =plt_morph%RTNL(N,L,NZ)
        RTLGP(N,L,NZ,NY,NX) =plt_morph%RTLGP(N,L,NZ)
        RTDNP(N,L,NZ,NY,NX) =plt_morph%RTDNP(N,L,NZ)
        RTVLP(N,L,NZ,NY,NX) =plt_morph%RTVLP(N,L,NZ)
        RTVLW(N,L,NZ,NY,NX) =plt_morph%RTVLW(N,L,NZ)
        RRAD1(N,L,NZ,NY,NX) =plt_morph%RRAD1(N,L,NZ)
        RRAD2(N,L,NZ,NY,NX) =plt_morph%RRAD2(N,L,NZ)
        RTARP(N,L,NZ,NY,NX) =plt_morph%RTARP(N,L,NZ)
        RTLGA(N,L,NZ,NY,NX) =plt_morph%RTLGA(N,L,NZ)
        RCO2M(N,L,NZ,NY,NX) =plt_rbgc%RCO2M(N,L,NZ)
        RCO2N(N,L,NZ,NY,NX) =plt_rbgc%RCO2N(N,L,NZ)
        RCO2A(N,L,NZ,NY,NX) =plt_rbgc%RCO2A(N,L,NZ)
        RCO2P(N,L,NZ,NY,NX) =plt_rbgc%RCO2P(N,L,NZ)
        RUPOXP(N,L,NZ,NY,NX)=plt_rbgc%RUPOXP(N,L,NZ)
        RCO2S(N,L,NZ,NY,NX) =plt_rbgc%RCO2S(N,L,NZ)
        RUPOXS(N,L,NZ,NY,NX)=plt_rbgc%RUPOXS(N,L,NZ)
        RUPCHS(N,L,NZ,NY,NX)=plt_rbgc%RUPCHS(N,L,NZ)
        RUPN2S(N,L,NZ,NY,NX)=plt_rbgc%RUPN2S(N,L,NZ)
        RUPN3S(N,L,NZ,NY,NX)=plt_rbgc%RUPN3S(N,L,NZ)
        RUPN3B(N,L,NZ,NY,NX)=plt_rbgc%RUPN3B(N,L,NZ)
        RUPHGS(N,L,NZ,NY,NX)=plt_rbgc%RUPHGS(N,L,NZ)
        RCOFLA(N,L,NZ,NY,NX)=plt_rbgc%RCOFLA(N,L,NZ)
        ROXFLA(N,L,NZ,NY,NX)=plt_rbgc%ROXFLA(N,L,NZ)
        RCHFLA(N,L,NZ,NY,NX)=plt_rbgc%RCHFLA(N,L,NZ)
        RN2FLA(N,L,NZ,NY,NX)=plt_rbgc%RN2FLA(N,L,NZ)
        RNHFLA(N,L,NZ,NY,NX)=plt_rbgc%RNHFLA(N,L,NZ)
        RHGFLA(N,L,NZ,NY,NX)=plt_rbgc%RHGFLA(N,L,NZ)
        RCODFA(N,L,NZ,NY,NX)=plt_rbgc%RCODFA(N,L,NZ)
        ROXDFA(N,L,NZ,NY,NX)=plt_rbgc%ROXDFA(N,L,NZ)
        RCHDFA(N,L,NZ,NY,NX)=plt_rbgc%RCHDFA(N,L,NZ)
        RN2DFA(N,L,NZ,NY,NX)=plt_rbgc%RN2DFA(N,L,NZ)
        RNHDFA(N,L,NZ,NY,NX)=plt_rbgc%RNHDFA(N,L,NZ)
        RHGDFA(N,L,NZ,NY,NX)=plt_rbgc%RHGDFA(N,L,NZ)
        RUNNHP(N,L,NZ,NY,NX)=plt_rbgc%RUNNHP(N,L,NZ)
        RUPNH4(N,L,NZ,NY,NX)=plt_rbgc%RUPNH4(N,L,NZ)
        RUONH4(N,L,NZ,NY,NX)=plt_rbgc%RUONH4(N,L,NZ)
        RUCNH4(N,L,NZ,NY,NX)=plt_rbgc%RUCNH4(N,L,NZ)
        RUNNBP(N,L,NZ,NY,NX)=plt_rbgc%RUNNBP(N,L,NZ)
        RUPNHB(N,L,NZ,NY,NX)=plt_rbgc%RUPNHB(N,L,NZ)
        RUONHB(N,L,NZ,NY,NX)=plt_rbgc%RUONHB(N,L,NZ)
        RUCNHB(N,L,NZ,NY,NX)=plt_rbgc%RUCNHB(N,L,NZ)
        RUNNOP(N,L,NZ,NY,NX)=plt_rbgc%RUNNOP(N,L,NZ)
        RUPNO3(N,L,NZ,NY,NX)=plt_rbgc%RUPNO3(N,L,NZ)
        RUONO3(N,L,NZ,NY,NX)=plt_rbgc%RUONO3(N,L,NZ)
        RUCNO3(N,L,NZ,NY,NX)=plt_rbgc%RUCNO3(N,L,NZ)
        RUNNXP(N,L,NZ,NY,NX)=plt_rbgc%RUNNXP(N,L,NZ)
        RUPNOB(N,L,NZ,NY,NX)=plt_rbgc%RUPNOB(N,L,NZ)
        RUONOB(N,L,NZ,NY,NX)=plt_rbgc%RUONOB(N,L,NZ)
        RUCNOB(N,L,NZ,NY,NX)=plt_rbgc%RUCNOB(N,L,NZ)
        RUPP2P(N,L,NZ,NY,NX)=plt_rbgc%RUPP2P(N,L,NZ)
        RUPH2P(N,L,NZ,NY,NX)=plt_rbgc%RUPH2P(N,L,NZ)
        RUOH2P(N,L,NZ,NY,NX)=plt_rbgc%RUOH2P(N,L,NZ)
        RUCH2P(N,L,NZ,NY,NX)=plt_rbgc%RUCH2P(N,L,NZ)
        RUPP2B(N,L,NZ,NY,NX)=plt_rbgc%RUPP2B(N,L,NZ)
        RUPH2B(N,L,NZ,NY,NX)=plt_rbgc%RUPH2B(N,L,NZ)
        RUOH2B(N,L,NZ,NY,NX)=plt_rbgc%RUOH2B(N,L,NZ)
        RUCH2B(N,L,NZ,NY,NX)=plt_rbgc%RUCH2B(N,L,NZ)
        RUPP1P(N,L,NZ,NY,NX)=plt_rbgc%RUPP1P(N,L,NZ)
        RUPH1P(N,L,NZ,NY,NX)=plt_rbgc%RUPH1P(N,L,NZ)
        RUOH1P(N,L,NZ,NY,NX)=plt_rbgc%RUOH1P(N,L,NZ)
        RUCH1P(N,L,NZ,NY,NX)=plt_rbgc%RUCH1P(N,L,NZ)
        RUPP1B(N,L,NZ,NY,NX)=plt_rbgc%RUPP1B(N,L,NZ)
        RUPH1B(N,L,NZ,NY,NX)=plt_rbgc%RUPH1B(N,L,NZ)
        RUOH1B(N,L,NZ,NY,NX)=plt_rbgc%RUOH1B(N,L,NZ)
        RUCH1B(N,L,NZ,NY,NX)=plt_rbgc%RUCH1B(N,L,NZ)
        ROXYP(N,L,NZ,NY,NX) =plt_rbgc%ROXYP(N,L,NZ)
        UPWTR(N,L,NZ,NY,NX) =plt_ew%UPWTR(N,L,NZ)
        WTRTL(N,L,NZ,NY,NX) =plt_biom%WTRTL(N,L,NZ)
        WTRTD(N,L,NZ,NY,NX) =plt_biom%WTRTD(N,L,NZ)
        WSRTL(N,L,NZ,NY,NX) =plt_biom%WSRTL(N,L,NZ)
        WFR(N,L,NZ,NY,NX)   =plt_rbgc%WFR(N,L,NZ)
      ENDDO
    ENDDO

    DO NR=1,NRT(NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=plt_morph%NINR(NR,NZ)
      DO N=1,2
        RTWT1E(N,NR,1:npelms,NZ,NY,NX) =plt_biom%RTWT1E(N,NR,1:npelms,NZ)
        RTDP1(N,NR,NZ,NY,NX) =plt_morph%RTDP1(N,NR,NZ)
      ENDDO
      DO L=1,JZ
        DO N=1,2
          WTRT1E(1:npelms,N,L,NR,NZ,NY,NX) =plt_biom%WTRT1E(1:npelms,N,L,NR,NZ)
          WTRT2E(1:npelms,N,L,NR,NZ,NY,NX) =plt_biom%WTRT2E(1:npelms,N,L,NR,NZ)
          RTLG1(N,L,NR,NZ,NY,NX) =plt_morph%RTLG1(N,L,NR,NZ)
          RTLG2(N,L,NR,NZ,NY,NX) =plt_morph%RTLG2(N,L,NR,NZ)
          RTN2(N,L,NR,NZ,NY,NX)  =plt_morph%RTN2(N,L,NR,NZ)
        ENDDO
      ENDDO
    ENDDO

    EHVST(1:2,1:4,NZ,I,NY,NX)=plt_distb%EHVST(1:2,1:4,NZ)


    DO L=NU(NY,NX),NI(NZ,NY,NX)

      DO K=1,jcplx
        DO N=1,2
          DO NE=1,npelms
            RDFOME(NE,N,K,L,NZ,NY,NX)=plt_rbgc%RDFOME(NE,N,K,L,NZ)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO NE=1,npelms
      DO M=1,jsken
        DO N=0,JP
          CFOPE(N,M,NE,NZ,NY,NX)=plt_soilchem%CFOPE(N,M,NE,NZ)
        enddo
      enddo
    ENDDO
    RRAD1M(2,NZ,NY,NX)=plt_morph%RRAD1M(2,NZ)
    RRAD2M(2,NZ,NY,NX)=plt_morph%RRAD2M(2,NZ)
    PORT(2,NZ,NY,NX)  =plt_morph%PORT(2,NZ)
    UPMXZH(2,NZ,NY,NX)=plt_rbgc%UPMXZH(2,NZ)
    UPKMZH(2,NZ,NY,NX)=plt_rbgc%UPKMZH(2,NZ)
    UPMNZH(2,NZ,NY,NX)=plt_rbgc%UPMNZH(2,NZ)
    UPMXZO(2,NZ,NY,NX)=plt_rbgc%UPMXZO(2,NZ)
    UPKMZO(2,NZ,NY,NX)=plt_rbgc%UPKMZO(2,NZ)
    UPMNZO(2,NZ,NY,NX)=plt_rbgc%UPMNZO(2,NZ)
    UPMXPO(2,NZ,NY,NX)=plt_rbgc%UPMXPO(2,NZ)
    UPKMPO(2,NZ,NY,NX)=plt_rbgc%UPKMPO(2,NZ)
    UPMNPO(2,NZ,NY,NX)=plt_rbgc%UPMNPO(2,NZ)
    RSRR(2,NZ,NY,NX)  =plt_morph%RSRR(2,NZ)
    RSRA(2,NZ,NY,NX)  =plt_morph%RSRA(2,NZ)
    DO N=1,2
      PORTX(N,NZ,NY,NX)=plt_morph%PORTX(N,NZ)
      RRADP(N,NZ,NY,NX)=plt_morph%RRADP(N,NZ)
      DMVL(N,NZ,NY,NX) =plt_morph%DMVL(N,NZ)
      RTLG1X(N,NZ,NY,NX)=plt_morph%RTLG1X(N,NZ)
      RTLG2X(N,NZ,NY,NX)=plt_morph%RTLG2X(N,NZ)
      RRAD1X(N,NZ,NY,NX)=plt_morph%RRAD1X(N,NZ)
      RRAD2X(N,NZ,NY,NX)=plt_morph%RRAD2X(N,NZ)
      RTAR1X(N,NZ,NY,NX)=plt_morph%RTAR1X(N,NZ)
      RTAR2X(N,NZ,NY,NX)=plt_morph%RTAR2X(N,NZ)
    enDDO
  ENDDO

  end subroutine PlantAPIRecv


!------------------------------------------------------------------------------------------

  subroutine PlantAPISend(I,J,NY,NX)
  !
  !DESCRIPTION
  !Send data to plant model
  use PlantAPIData, only : plt_rad
  use EcoSIMConfig, only : jsken=>jskenc,jcplx => jcplxc
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer :: K,L,M,N,NB,NZ,NR,I1,NE

  plt_site%LYRC=LYRC
  I1=I+1;if(I1>LYRC)I1=1
  plt_site%ZERO=ZERO
  plt_site%ZERO2=ZERO2
  plt_site%ALAT=ALAT(NY,NX)
  plt_site%ATCA=ATCA(NY,NX)
  plt_morph%ARLSS=ARLSS(NY,NX)
  plt_morph%ARLFC=ARLFC(NY,NX)
  plt_site%ALT=ALT(NY,NX)
  plt_site%CCO2EI=CCO2EI(NY,NX)
  plt_site%CO2EI=CO2EI(NY,NX)
  plt_site%COXYE=AtmGgms(idg_O2,NY,NX)
  plt_bgcr%CNETX=CNETX(NY,NX)
  plt_site%CO2E=CO2E(NY,NX)
  plt_site%CCO2E=AtmGgms(idg_CO2,NY,NX)
  plt_site%CCH4E=AtmGgms(idg_CH4,NY,NX)
  plt_site%CZ2OE=AtmGgms(idg_N2O,NY,NX)
  plt_site%CH2GE=AtmGgms(idg_H2,NY,NX)
  plt_site%CNH3E=AtmGgms(idg_NH3,NY,NX)
  plt_site%DYLX=DYLX(NY,NX)
  plt_site%DYLN=DYLN(NY,NX)
  plt_ew%DPTHS=DPTHS(NY,NX)
  plt_site%DYLM=DYLM(NY,NX)
  plt_site%IETYP=IETYP(NY,NX)
  plt_site%IYRC=IYRC
  plt_site%NL=NL(NY,NX)
  plt_site%NP0=NP0(NY,NX)
  plt_site%NJ=NJ(NY,NX)
  plt_site%NP=NP(NY,NX)
  plt_site%NU=NU(NY,NX)
  plt_site%OXYE=OXYE(NY,NX)
  plt_ew%RAB=RAB(NY,NX)
  plt_ew%RIB=RIB(NY,NX)
  plt_rad%SSINN=SSINN(NY,NX)
  plt_rad%SSIN=SSIN(NY,NX)
  plt_ew%TKW=TKW(1,NY,NX)
  plt_ew%TKA=TKA(NY,NX)
  plt_rad%THRMGX=THRMGX(NY,NX)
  plt_rad%THS=THS(NY,NX)
  plt_ew%VPA=VPA(NY,NX)
  plt_distb%XCORP=XCORP(NY,NX)
  plt_site%ZNOON =ZNOON(NY,NX)
  plt_site%ZEROS2=ZEROS2(NY,NX)
  plt_site%ZEROS =ZEROS(NY,NX)
  plt_ew%ZR=ZR(NY,NX)
  plt_morph%ZT=ZT(NY,NX)
  plt_ew%ZD=ZD(NY,NX)
  plt_site%IDATA(:)=IDATA(:)
  plt_distb%DCORP=DCORP(I,NY,NX)
  plt_distb%ITILL=ITILL(I,NY,NX)

  plt_morph%ZL(0)=ZL(0,NY,NX)
  DO  L=1,JC
    plt_morph%ZL(L)=ZL(L,NY,NX)
    plt_rad%TAUS(L)=TAUS(L,NY,NX)
    plt_rad%TAU0(L)=TAU0(L,NY,NX)
  ENDDO
  plt_rad%TAUS(JC+1)=TAUS(JC+1,NY,NX)
  plt_rad%TAU0(JC+1)=TAU0(JC+1,NY,NX)

  DO N=1,JLI
    plt_rad%ZSIN(N)=ZSIN(N)
  ENDDO

  DO L=1,NL(NY,NX)
    plt_soilchem%CNDU(L)=CNDU(L,NY,NX)
    plt_soilchem%GasDifc(idg_CO2,L)=GasDifc(idg_CO2,L,NY,NX)
    plt_soilchem%GasDifc(idg_CH4,L)=GasDifc(idg_CH4,L,NY,NX)
    plt_soilchem%GasDifc(idg_H2,L) =GasDifc(idg_H2,L,NY,NX)
    plt_soilchem%GasDifc(idg_O2,L) =GasDifc(idg_O2,L,NY,NX)
    plt_soilchem%GasDifc(idg_N2O,L)=GasDifc(idg_N2O,L,NY,NX)
    plt_soilchem%GasDifc(idg_NH3,L)=GasDifc(idg_NH3,L,NY,NX)
    plt_soilchem%RSCS(L)=RSCS(L,NY,NX)
    plt_site%DPTHZ(L)=DPTHZ(L,NY,NX)
  ENDDO
  plt_allom%FVRN(:) = FVRN(:)
  DO L=1,NL(NY,NX)
    plt_soilchem%trc_gasml(idg_CO2,L)=trc_gasml(idg_CO2,L,NY,NX)
    plt_soilchem%trc_gasml(idg_O2,L)=trc_gasml(idg_O2,L,NY,NX)
  ENDDO

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)     =AREA(3,L,NY,NX)
    plt_soilchem%BKDS(L)  =BKDS(L,NY,NX)

    plt_soilchem%trc_solcl(ids_H2PO4,L) =trc_solcl(ids_H2PO4,L,NY,NX)
    plt_soilchem%trc_solcl(ids_H1PO4B,L)=trc_solcl(ids_H1PO4B,L,NY,NX)
    plt_soilchem%trc_solcl(ids_H2PO4B,L)=trc_solcl(ids_H2PO4B,L,NY,NX)
    plt_soilchem%trc_solcl(ids_NO3,L) =trc_solcl(ids_NO3,L,NY,NX)
    plt_soilchem%trc_solcl(ids_H1PO4,L) =trc_solcl(ids_H1PO4,L,NY,NX)
    plt_soilchem%trc_solcl(ids_NO3B,L) =trc_solcl(ids_NO3B,L,NY,NX)
    plt_soilchem%trc_solcl(ids_NH4,L) =trc_solcl(ids_NH4,L,NY,NX)
    plt_soilchem%trc_solcl(ids_NH4B,L) =trc_solcl(ids_NH4B,L,NY,NX)

    plt_soilchem%trc_solcl(idg_CH4,L) =trc_solcl(idg_CH4,L,NY,NX)
    plt_soilchem%trc_solcl(idg_N2O,L) =trc_solcl(idg_N2O,L,NY,NX)
    plt_soilchem%trc_solcl(idg_NH3,L) =trc_solcl(idg_NH3,L,NY,NX)
    plt_soilchem%trc_solcl(idg_NH3B,L) =trc_solcl(idg_NH3B,L,NY,NX)
    plt_soilchem%trc_solcl(idg_H2,L) =trc_solcl(idg_H2,L,NY,NX)

    plt_soilchem%SolDifc(idg_CO2,L) =SolDifc(idg_CO2,L,NY,NX)
    plt_soilchem%SolDifc(idg_CH4,L) =SolDifc(idg_CH4,L,NY,NX)
    plt_soilchem%SolDifc(idg_H2,L) =SolDifc(idg_H2,L,NY,NX)
    plt_soilchem%SolDifc(idg_O2,L) =SolDifc(idg_O2,L,NY,NX)
    plt_soilchem%SolDifc(ids_H1PO4,L) =SolDifc(ids_H1PO4,L,NY,NX)
    plt_soilchem%SolDifc(ids_NO3,L) =SolDifc(ids_NO3,L,NY,NX)
    plt_soilchem%SolDifc(idg_NH3,L) =SolDifc(idg_NH3,L,NY,NX)
    plt_soilchem%SolDifc(idg_N2O,L) =SolDifc(idg_N2O,L,NY,NX)

    plt_soilchem%trc_gascl(idg_CH4,L) =trc_gascl(idg_CH4,L,NY,NX)
    plt_soilchem%trc_gascl(idg_N2O,L) =trc_gascl(idg_N2O,L,NY,NX)
    plt_soilchem%trc_gascl(idg_NH3,L) =trc_gascl(idg_NH3,L,NY,NX)
    plt_soilchem%trc_gascl(idg_H2,L) =trc_gascl(idg_H2,L,NY,NX)
    plt_soilchem%CORGC(L) =CORGC(L,NY,NX)
    plt_site%CDPTHZ(L)    =CDPTHZ(L,NY,NX)
    plt_site%FMPR(L)      =FMPR(L,NY,NX)

    plt_soilchem%trc_solml(idg_CO2,L)  =trc_solml(idg_CO2,L,NY,NX)
    plt_soilchem%trc_solml(idg_CH4,L)  =trc_solml(idg_CH4,L,NY,NX)
    plt_soilchem%trc_solml(idg_H2,L)  =trc_solml(idg_H2,L,NY,NX)
    plt_soilchem%trc_solml(idg_N2O,L)  =trc_solml(idg_N2O,L,NY,NX)
    plt_soilchem%trc_solml(idg_O2,L)  =trc_solml(idg_O2,L,NY,NX)
    plt_soilchem%trc_solml(idg_NH3,L) =trc_solml(idg_NH3,L,NY,NX)
    plt_soilchem%trc_solml(idg_NH3B,L) =trc_solml(idg_NH3B,L,NY,NX)
    plt_soilchem%trc_solml(ids_H1PO4,L) =trc_solml(ids_H1PO4,L,NY,NX)
    plt_soilchem%trc_solml(ids_H2PO4,L) =trc_solml(ids_H2PO4,L,NY,NX)
    plt_soilchem%trc_solml(ids_H1PO4B,L) =trc_solml(ids_H1PO4B,L,NY,NX)
    plt_soilchem%trc_solml(ids_H2PO4B,L) =trc_solml(ids_H2PO4B,L,NY,NX)
    plt_soilchem%trc_solml(ids_NO3,L) =trc_solml(ids_NO3,L,NY,NX)
    plt_soilchem%trc_solml(ids_NO3B,L) =trc_solml(ids_NO3B,L,NY,NX)
    plt_soilchem%trc_solml(ids_NH4,L) =trc_solml(ids_NH4,L,NY,NX)
    plt_soilchem%trc_solml(ids_NH4B,L) =trc_solml(ids_NH4B,L,NY,NX)

    plt_soilchem%GSolbility(idg_CO2,L) =GSolbility(idg_CO2,L,NY,NX)
    plt_soilchem%GSolbility(idg_O2,L) =GSolbility(idg_O2,L,NY,NX)
    plt_soilchem%GSolbility(idg_CH4,L) =GSolbility(idg_CH4,L,NY,NX)
    plt_soilchem%GSolbility(idg_N2O,L) =GSolbility(idg_N2O,L,NY,NX)
    plt_soilchem%GSolbility(idg_NH3,L) =GSolbility(idg_NH3,L,NY,NX)
    plt_soilchem%GSolbility(idg_H2,L) =GSolbility(idg_H2,L,NY,NX)

    plt_ew%PSIST(L)       =PSIST(L,NY,NX)
    plt_bgcr%RPO4Y(L)     =RPO4Y(L,NY,NX)
    plt_bgcr%RPOBY(L)     =RPOBY(L,NY,NX)
    plt_bgcr%RP14Y(L)     =RP14Y(L,NY,NX)
    plt_bgcr%RP1BY(L)     =RP1BY(L,NY,NX)
    plt_bgcr%RNO3Y(L)     =RNO3Y(L,NY,NX)
    plt_bgcr%RNH4Y(L)     =RNH4Y(L,NY,NX)
    plt_bgcr%RNHBY(L)     =RNHBY(L,NY,NX)
    plt_bgcr%RN3BY(L)     =RN3BY(L,NY,NX)
    plt_bgcr%ROXYF(L)     =ROXYF(L,NY,NX)
    plt_bgcr%RCO2F(L)     =RCO2F(L,NY,NX)
    plt_bgcr%ROXYL(L)     =ROXYL(L,NY,NX)
    plt_bgcr%ROXYY(L)     =ROXYY(L,NY,NX)
    plt_ew%TKS(L)         =TKS(L,NY,NX)


    plt_soilchem%THETW(L) =THETW(L,NY,NX)
    plt_soilchem%THETY(L) =THETY(L,NY,NX)
    plt_soilchem%TFND(L)  =TFND(L,NY,NX)
    plt_soilchem%VOLX(L)  =VOLX(L,NY,NX)
    plt_soilchem%trcs_VLN(ids_H1PO4B,L) =trcs_VLN(ids_H1PO4B,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NO3,L) =trcs_VLN(ids_NO3,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_H1PO4,L) =trcs_VLN(ids_H1PO4,L,NY,NX)
    plt_soilchem%VOLY(L)  =VOLY(L,NY,NX)
    plt_soilchem%VOLI(L)  =VOLI(L,NY,NX)
    plt_soilchem%VOLW(L)  =VOLW(L,NY,NX)
    plt_soilchem%VOLA(L)  =VOLA(L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NO3B,L) =trcs_VLN(ids_NO3B,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NH4,L) =trcs_VLN(ids_NH4,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NH4B,L) =trcs_VLN(ids_NH4B,L,NY,NX)

    plt_site%DLYR3(L)     =DLYR(3,L,NY,NX)
    DO K=1,jcplx
      plt_soilchem%FOSRH(K,L)=FOSRH(K,L,NY,NX)
      plt_soilchem%OQC(K,L)=OQC(K,L,NY,NX)
      plt_soilchem%OQN(K,L)=OQN(K,L,NY,NX)
      plt_soilchem%OQP(K,L)=OQP(K,L,NY,NX)
    ENDDO
  ENDDO

  DO NZ=1,NP0(NY,NX)
!plant properties begin
    plt_photo%ICTYP(NZ)=ICTYP(NZ,NY,NX)
    plt_pheno%IGTYP(NZ)=IGTYP(NZ,NY,NX)
    plt_pheno%ISTYP(NZ)=ISTYP(NZ,NY,NX)
    plt_pheno%IDTYP(NZ)=IDTYP(NZ,NY,NX)
    plt_pheno%IWTYP(NZ)=IWTYP(NZ,NY,NX)
    plt_pheno%IPTYP(NZ)=IPTYP(NZ,NY,NX)
    plt_pheno%IBTYP(NZ)=IBTYP(NZ,NY,NX)
    plt_pheno%ZTYPI(NZ)=ZTYPI(NZ,NY,NX)
    plt_morph%IRTYP(NZ)=IRTYP(NZ,NY,NX)
    plt_morph%INTYP(NZ)=INTYP(NZ,NY,NX)
    plt_morph%MY(NZ)=MY(NZ,NY,NX)
    plt_photo%VCMX(NZ)=VCMX(NZ,NY,NX)
    plt_photo%VOMX(NZ)=VOMX(NZ,NY,NX)
    plt_photo%VCMX4(NZ)=VCMX4(NZ,NY,NX)
    plt_photo%XKCO2(NZ)=XKCO2(NZ,NY,NX)
    plt_photo%XKO2(NZ)=XKO2(NZ,NY,NX)
    plt_photo%XKCO24(NZ)=XKCO24(NZ,NY,NX)
    plt_photo%RUBP(NZ)=RUBP(NZ,NY,NX)
    plt_photo%PEPC(NZ)=PEPC(NZ,NY,NX)
    plt_photo%ETMX(NZ)=ETMX(NZ,NY,NX)
    plt_photo%CHL(NZ)=CHL(NZ,NY,NX)
    plt_photo%CHL4(NZ)=CHL4(NZ,NY,NX)
    plt_photo%FCO2(NZ)=FCO2(NZ,NY,NX)

    plt_pheno%XRNI(NZ)=XRNI(NZ,NY,NX)
    plt_pheno%XRLA(NZ)=XRLA(NZ,NY,NX)
    plt_pheno%CTC(NZ)=CTC(NZ,NY,NX)
    plt_morph%WDLF(NZ)=WDLF(NZ,NY,NX)
    plt_pheno%PB(NZ)=PB(NZ,NY,NX)
    plt_morph%XTLI(NZ)=XTLI(NZ,NY,NX)
    plt_pheno%XDL(NZ)=XDL(NZ,NY,NX)
    plt_pheno%XPPD(NZ)=XPPD(NZ,NY,NX)
    plt_morph%SLA1(NZ)=SLA1(NZ,NY,NX)
    plt_morph%SSL1(NZ)=SSL1(NZ,NY,NX)
    plt_morph%SNL1(NZ)=SNL1(NZ,NY,NX)
    DO  N=1,JLI
      plt_morph%CLASS(N,NZ)=CLASS(N,NZ,NY,NX)
    ENDDO
    plt_morph%CFI(NZ)=CFI(NZ,NY,NX)
    plt_morph%ANGBR(NZ)=ANGBR(NZ,NY,NX)
    plt_morph%ANGSH(NZ)=ANGSH(NZ,NY,NX)
    plt_morph%STMX(NZ)=STMX(NZ,NY,NX)
    plt_morph%SDMX(NZ)=SDMX(NZ,NY,NX)
    plt_morph%GRMX(NZ)=GRMX(NZ,NY,NX)
    plt_morph%GRDM(NZ)=GRDM(NZ,NY,NX)
    plt_pheno%GFILL(NZ)=GFILL(NZ,NY,NX)
    plt_biom%WTSTDI(NZ)=WTSTDI(NZ,NY,NX)

!initial root values
    DO N=1,MY(NZ,NY,NX)
      plt_morph%RRAD1M(N,NZ)=RRAD1M(N,NZ,NY,NX)
      plt_morph%RRAD2M(N,NZ)=RRAD2M(N,NZ,NY,NX)
      plt_morph%PORT(N,NZ)  =PORT(N,NZ,NY,NX)
      plt_morph%RSRR(N,NZ)  =RSRR(N,NZ,NY,NX)
      plt_morph%RSRA(N,NZ)=RSRA(N,NZ,NY,NX)
      plt_rbgc%UPMXZH(N,NZ)=UPMXZH(N,NZ,NY,NX)
      plt_rbgc%UPKMZH(N,NZ)=UPKMZH(N,NZ,NY,NX)
      plt_rbgc%UPMNZH(N,NZ)=UPMNZH(N,NZ,NY,NX)
      plt_rbgc%UPMXZO(N,NZ)=UPMXZO(N,NZ,NY,NX)
      plt_rbgc%UPKMZO(N,NZ)=UPKMZO(N,NZ,NY,NX)
      plt_rbgc%UPMNZO(N,NZ)=UPMNZO(N,NZ,NY,NX)
      plt_rbgc%UPMXPO(N,NZ)=UPMXPO(N,NZ,NY,NX)
      plt_rbgc%UPKMPO(N,NZ)=UPKMPO(N,NZ,NY,NX)
      plt_rbgc%UPMNPO(N,NZ)=UPMNPO(N,NZ,NY,NX)
    ENDDO
    plt_pheno%PR(NZ)=PR(NZ,NY,NX)
    plt_pheno%PTSHT(NZ)=PTSHT(NZ,NY,NX)
    plt_morph%RTFQ(NZ)=RTFQ(NZ,NY,NX)
    plt_ew%OSMO(NZ)=OSMO(NZ,NY,NX)
    plt_photo%RCS(NZ)=RCS(NZ,NY,NX)
    plt_photo%RSMX(NZ)=RSMX(NZ,NY,NX)

    plt_allom%DMLF(NZ) =DMLF(NZ,NY,NX)
    plt_allom%DMSHE(NZ)=DMSHE(NZ,NY,NX)
    plt_allom%DMSTK(NZ)=DMSTK(NZ,NY,NX)
    plt_allom%DMRSV(NZ)=DMRSV(NZ,NY,NX)
    plt_allom%DMHSK(NZ)=DMHSK(NZ,NY,NX)
    plt_allom%DMEAR(NZ)=DMEAR(NZ,NY,NX)
    plt_allom%DMGR(NZ) =DMGR(NZ,NY,NX)
    plt_allom%DMRT(NZ) =DMRT(NZ,NY,NX)
    plt_allom%DMND(NZ) =DMND(NZ,NY,NX)
    plt_allom%CNLF(NZ) =CNLF(NZ,NY,NX)
    plt_allom%CNSHE(NZ)=CNSHE(NZ,NY,NX)
    plt_allom%CNSTK(NZ)=CNSTK(NZ,NY,NX)
    plt_allom%CNRSV(NZ)=CNRSV(NZ,NY,NX)
    plt_allom%CNHSK(NZ)=CNHSK(NZ,NY,NX)
    plt_allom%CNEAR(NZ)=CNEAR(NZ,NY,NX)
    plt_allom%CNGR(NZ)=CNGR(NZ,NY,NX)
    plt_allom%CNRT(NZ)=CNRT(NZ,NY,NX)
    plt_allom%CNND(NZ)=CNND(NZ,NY,NX)
    plt_allom%CPLF(NZ)=CPLF(NZ,NY,NX)
    plt_allom%CPSHE(NZ)=CPSHE(NZ,NY,NX)
    plt_allom%CPSTK(NZ)=CPSTK(NZ,NY,NX)
    plt_allom%CPRSV(NZ)=CPRSV(NZ,NY,NX)
    plt_allom%CPHSK(NZ)=CPHSK(NZ,NY,NX)
    plt_allom%CPEAR(NZ)=CPEAR(NZ,NY,NX)
    plt_allom%CPGR(NZ)=CPGR(NZ,NY,NX)
    plt_allom%CPRT(NZ)=CPRT(NZ,NY,NX)
    plt_allom%CPND(NZ)=CPND(NZ,NY,NX)

!plant properties end

    plt_morph%ARLFS(NZ)=ARLFS(NZ,NY,NX)
    plt_distb%IYRX(NZ)=IYRX(NZ,NY,NX)
    plt_distb%IDAYX(NZ)=IDAYX(NZ,NY,NX)
    plt_distb%IYRY(NZ)=IYRY(NZ,NY,NX)
    plt_rad%RADP(NZ)=RADP(NZ,NY,NX)
    plt_rad%RADC(NZ)=RADC(NZ,NY,NX)
    plt_ew%FLWC(NZ)=FLWC(NZ,NY,NX)

    plt_site%PPZ(NZ)=PPZ(NZ,NY,NX)
    plt_distb%IDAYY(NZ)=IDAYY(NZ,NY,NX)
    plt_morph%CFX(NZ)=CFX(NZ,NY,NX)
    plt_site%DATAP(NZ)=DATAP(NZ,NY,NX)
    plt_pheno%GROUPI(NZ)=GROUPI(NZ,NY,NX)
    plt_biom%WTSHTA(NZ)=WTSHTA(NZ,NY,NX)


    DO NB=1,NBR(NZ,NY,NX)
      plt_pheno%VRNL(NB,NZ)=VRNL(NB,NZ,NY,NX)
      plt_pheno%VRNX(NB,NZ)=VRNX(NB,NZ,NY,NX)
    ENDDO

    DO L=1,JC
      DO  M=1,JSA
        DO  N=1,JLI
          plt_rad%PAR(N,M,L,NZ)=PAR(N,M,L,NZ,NY,NX)
          plt_rad%PARDIF(N,M,L,NZ)=PARDIF(N,M,L,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NL(NY,NX)
    DO M=1,NPH
      plt_site%VOLWM(M,L)=VOLWM(M,L,NY,NX)
      plt_site%VOLPM(M,L)=VOLPM(M,L,NY,NX)
      plt_site%TORT(M,L)=TORT(M,L,NY,NX)
      plt_site%FILM(M,L)=FILM(M,L,NY,NX)
      plt_soilchem%DFGS(M,L)=DFGS(M,L,NY,NX)
    ENDDO
  ENDDO

! sent variables also modified
  plt_site%IFLGT=IFLGT(NY,NX)
  plt_site%VOLWOU=VOLWOU
  plt_site%PPT=PPT(NY,NX)
  plt_bgcr%RECO=RECO(NY,NX)
  plt_biom%WTSTGT=WTSTGT(NY,NX)
  plt_bgcr%ZESNC(1:npelms)=ZESNC(1:npelms,NY,NX)
  plt_morph%ARSTC=ARSTC(NY,NX)
  plt_ew%TSH=TSH(NY,NX)
  plt_ew%TVOLWC=TVOLWC(NY,NX)
  plt_bgcr%TNBP=TNBP(NY,NX)
  plt_bgcr%TGPP=TGPP(NY,NX)
  plt_ew%TLEC=TLEC(NY,NX)
  plt_ew%TSHC=TSHC(NY,NX)
  plt_bgcr%TRAU=TRAU(NY,NX)
  plt_ew%UVOLO    =UVOLO(NY,NX)
  plt_distb%XHVSTE(1:npelms)=XHVSTE(1:npelms,NY,NX)
  plt_rad%TRN     =TRN(NY,NX)
  plt_rbgc%TN2OZ=TN2OZ(NY,NX)
  plt_ew%TEVAPC=TEVAPC(NY,NX)
  plt_ew%TLE=TLE(NY,NX)
  plt_rbgc%TCO2Z=TCO2Z(NY,NX)
  plt_ew%TGH    =TGH(NY,NX)
  plt_ew%TEVAPP =TEVAPP(NY,NX)
  plt_bgcr%TH2GZ  =TH2GZ(NY,NX)
  plt_ew%THFLXC =THFLXC(NY,NX)
  plt_ew%THRMC  =THRMC(NY,NX)
  plt_rbgc%TOXYZ=TOXYZ(NY,NX)
  plt_rbgc%TNH3Z=TNH3Z(NY,NX)
  plt_ew%TVOLWP =TVOLWP(NY,NX)
  plt_ew%TENGYC =TENGYC(NY,NX)
  plt_bgcr%TCCAN=TCCAN(NY,NX)
  plt_rbgc%TCH4Z=TCH4Z(NY,NX)
  plt_distb%FERT(1:20)=FERT(1:20,I1,NY,NX)

  DO  L=1,JC
    plt_morph%ARSTT(L)=ARSTT(L,NY,NX)
    plt_biom%WGLFT(L)=WGLFT(L,NY,NX)
    plt_morph%ARLFT(L)=ARLFT(L,NY,NX)
  ENDDO

  plt_allom%FWOODE(1:npelms,1:n_pltlitrk)=FWOODE(1:npelms,1:n_pltlitrk)
  plt_allom%FWODRE(1:npelms,1:n_pltlitrk)=FWODRE(1:npelms,1:n_pltlitrk)
  plt_allom%FWODBE(1:npelms,1:n_pltlitrk)=FWODBE(1:npelms,1:n_pltlitrk)
  plt_allom%FWODLE(1:npelms,1:n_pltlitrk)=FWODLE(1:npelms,1:n_pltlitrk)

  DO L=0,NL(NY,NX)
    DO K=1,jcplx
      plt_bgcr%XOQCS(K,L)=XOQCS(K,L,NY,NX)
      plt_bgcr%XOQNS(K,L)=XOQNS(K,L,NY,NX)
      plt_bgcr%XOQPS(K,L)=XOQPS(K,L,NY,NX)
    ENDDO
  ENDDO

  DO L=0,JZ
    plt_bgcr%RPOBX(L)=RPOBX(L,NY,NX)
    plt_bgcr%RP1BX(L)=RP1BX(L,NY,NX)
    plt_bgcr%RN3BX(L)=RN3BX(L,NY,NX)
    plt_bgcr%RNHBX(L)=RNHBX(L,NY,NX)
    plt_bgcr%RP14X(L)=RP14X(L,NY,NX)
    plt_bgcr%RPO4X(L)=RPO4X(L,NY,NX)
    plt_bgcr%RNO3X(L)=RNO3X(L,NY,NX)
    plt_bgcr%RNH4X(L)=RNH4X(L,NY,NX)
    plt_bgcr%ROXYX(L)=ROXYX(L,NY,NX)
    plt_ew%TUPHT(L)  =TUPHT(L,NY,NX)
    plt_ew%TUPWTR(L) =TUPWTR(L,NY,NX)
    DO  K=1,micpar%n_pltlitrk
      DO NE=1,npelms
        DO  M=1,jsken
          plt_bgcr%ESNT(M,NE,K,L)=ESNT(M,NE,K,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO L=1,JZ
    plt_rbgc%TUPNF(L)=TUPNF(L,NY,NX)
    plt_morph%RTDNT(L)=RTDNT(L,NY,NX)
    plt_rbgc%TLOXYP(L)=TLOXYP(L,NY,NX)
    plt_rbgc%TLCO2P(L)=TLCO2P(L,NY,NX)
    plt_rbgc%TLCH4P(L)=TLCH4P(L,NY,NX)
    plt_rbgc%TLN2OP(L)=TLN2OP(L,NY,NX)
    plt_rbgc%TLNH3P(L)=TLNH3P(L,NY,NX)
    plt_bgcr%TLH2GP(L)=TLH2GP(L,NY,NX)
    plt_rbgc%TCOFLA(L)=TCOFLA(L,NY,NX)
    plt_rbgc%TOXFLA(L)=TOXFLA(L,NY,NX)
    plt_rbgc%TCHFLA(L)=TCHFLA(L,NY,NX)
    plt_rbgc%TN2FLA(L)=TN2FLA(L,NY,NX)
    plt_rbgc%TNHFLA(L)=TNHFLA(L,NY,NX)
    plt_bgcr%THGFLA(L)=THGFLA(L,NY,NX)
    plt_bgcr%TCO2P(L) =TCO2P(L,NY,NX)
    plt_bgcr%TUPOXP(L)=TUPOXP(L,NY,NX)
    plt_bgcr%TCO2S(L) =TCO2S(L,NY,NX)
    plt_bgcr%TUPOXS(L)=TUPOXS(L,NY,NX)
    plt_bgcr%TUPCHS(L)=TUPCHS(L,NY,NX)
    plt_bgcr%TUPN2S(L)=TUPN2S(L,NY,NX)
    plt_bgcr%TUPN3S(L)=TUPN3S(L,NY,NX)
    plt_bgcr%TUPN3B(L)=TUPN3B(L,NY,NX)
    plt_bgcr%TUPHGS(L)=TUPHGS(L,NY,NX)
    plt_bgcr%TUPNH4(L)=TUPNH4(L,NY,NX)
    plt_bgcr%TUPNO3(L)=TUPNO3(L,NY,NX)
    plt_bgcr%TUPH1B(L)=TUPH1B(L,NY,NX)
    plt_bgcr%TUPH2B(L)=TUPH2B(L,NY,NX)
    plt_bgcr%TUPNHB(L)=TUPNHB(L,NY,NX)
    plt_bgcr%TUPNOB(L)=TUPNOB(L,NY,NX)
    plt_bgcr%TUPH1P(L)=TUPH1P(L,NY,NX)
    plt_bgcr%TUPH2P(L)=TUPH2P(L,NY,NX)
    DO  K=1,jcplx
      plt_bgcr%TDFOMP(K,L)=TDFOMP(K,L,NY,NX)
      plt_bgcr%TDFOMN(K,L)=TDFOMN(K,L,NY,NX)
      plt_bgcr%TDFOMC(K,L)=TDFOMC(K,L,NY,NX)
    ENDDO
  ENDDO


  DO NZ=1,NP0(NY,NX)
    plt_biom%WTRTE(1:npelms,NZ)=WTRTE(1:npelms,NZ,NY,NX)
    plt_biom%WTLFE(1:npelms,NZ)  =WTLFE(1:npelms,NZ,NY,NX)
    plt_biom%WTEARE(1:npelms,NZ) =WTEARE(1:npelms,NZ,NY,NX)
    plt_biom%WTRTSE(1:npelms,NZ) =WTRTSE(1:npelms,NZ,NY,NX)
    plt_biom%WTHSKE(1:npelms,NZ) =WTHSKE(1:npelms,NZ,NY,NX)
    plt_biom%WTRSVE(1:npelms,NZ) =WTRSVE(1:npelms,NZ,NY,NX)
    plt_biom%WTSTKE(1:npelms,NZ) =WTSTKE(1:npelms,NZ,NY,NX)
    plt_biom%WTSHEE(1:npelms,NZ) =WTSHEE(1:npelms,NZ,NY,NX)
    plt_biom%WTGRE(1:npelms,NZ)=WTGRE(1:npelms,NZ,NY,NX)
    plt_site%BALE(1:npelms,NZ)  =BALE(1:npelms,NZ,NY,NX)
    plt_biom%CEPOLP(1:npelms,NZ)=CEPOLP(1:npelms,NZ,NY,NX)
    plt_biom%EPOOLP(1:npelms,NZ)=EPOOLP(1:npelms,NZ,NY,NX)
    plt_biom%EPOLNP(1:npelms,NZ)=EPOLNP(1:npelms,NZ,NY,NX)
    plt_bgcr%HESNC(1:npelms,NZ) =HESNC(1:npelms,NZ,NY,NX)
    plt_distb%HVSTE(1:npelms,NZ)=HVSTE(1:npelms,NZ,NY,NX)
    plt_pheno%RSETE(1:npelms,NZ)=RSETE(1:npelms,NZ,NY,NX)
    plt_bgcr%TESN0(1:npelms,NZ)  =TESN0(1:npelms,NZ,NY,NX)
    plt_bgcr%TESNC(1:npelms,NZ)  =TESNC(1:npelms,NZ,NY,NX)
    plt_rbgc%TEUPTK(1:npelms,NZ) =TEUPTK(1:npelms,NZ,NY,NX)
    plt_distb%THVSTE(1:npelms,NZ)=THVSTE(1:npelms,NZ,NY,NX)
    plt_rbgc%UPOME(1:npelms,NZ)=UPOME(1:npelms,NZ,NY,NX)
    plt_biom%WTRVE(1:npelms,NZ) =WTRVE(1:npelms,NZ,NY,NX)
    plt_biom%WTSHTE(1:npelms,NZ) =WTSHTE(1:npelms,NZ,NY,NX)
    plt_biom%WTSTGE(1:npelms,NZ) =WTSTGE(1:npelms,NZ,NY,NX)
    plt_biom%WTNDE(1:npelms,NZ)  =WTNDE(1:npelms,NZ,NY,NX)

    plt_ew%TKCZ(NZ)=TKCZ(NZ,NY,NX)
    plt_photo%SO2(NZ)=SO2(NZ,NY,NX)
    plt_ew%PSILZ(NZ)=PSILZ(NZ,NY,NX)
    plt_ew%HFLXC(NZ)=HFLXC(NZ,NY,NX)
    plt_photo%CO2Q(NZ)=CO2Q(NZ,NY,NX)
    plt_photo%O2L(NZ)=O2L(NZ,NY,NX)

    plt_photo%CO2L(NZ)=CO2L(NZ,NY,NX)
    plt_distb%EHVST(1:2,1:4,NZ)=EHVST(1:2,1:4,NZ,I,NY,NX)

    plt_pheno%IDTH(NZ)=IDTH(NZ,NY,NX)
    plt_distb%IYR0(NZ)=IYR0(NZ,NY,NX)
    plt_morph%NNOD(NZ)=NNOD(NZ,NY,NX)

    plt_pheno%IDTHP(NZ)=IDTHP(NZ,NY,NX)
    plt_pheno%IDTHR(NZ)=IDTHR(NZ,NY,NX)
    plt_allom%CNRTS(NZ)=CNRTS(NZ,NY,NX)
    plt_allom%CPRTS(NZ)=CPRTS(NZ,NY,NX)

    plt_rbgc%ZEROQ(NZ)=ZEROQ(NZ,NY,NX)
    plt_pheno%SSTX(NZ)=SSTX(NZ,NY,NX)
    plt_rad%FRADP(NZ) =FRADP(NZ,NY,NX)
    plt_bgcr%RNH3C(NZ)=RNH3C(NZ,NY,NX)
    plt_ew%DTKC(NZ)   =DTKC(NZ,NY,NX)
    plt_pheno%ZTYP(NZ)=ZTYP(NZ,NY,NX)
    plt_pheno%HTC(NZ) =HTC(NZ,NY,NX)
    plt_ew%EP(NZ)      =EP(NZ,NY,NX)

    plt_biom%WVSTK(NZ) =WVSTK(NZ,NY,NX)

    plt_photo%CHILL(NZ)=CHILL(NZ,NY,NX)
    plt_ew%PSILO(NZ)=PSILO(NZ,NY,NX)

    plt_ew%TCC(NZ)=TCC(NZ,NY,NX)
    plt_allom%CWSRT(NZ)=CWSRT(NZ,NY,NX)
    plt_photo%RSMH(NZ)=RSMH(NZ,NY,NX)
    plt_biom%WTRTA(NZ)=WTRTA(NZ,NY,NX)
    plt_morph%CF(NZ)=CF(NZ,NY,NX)

    plt_site%PPI(NZ)=PPI(NZ,NY,NX)
    plt_site%PPX(NZ)=PPX(NZ,NY,NX)
    plt_distb%VOXYF(NZ)=VOXYF(NZ,NY,NX)
    plt_ew%ENGYX(NZ)=ENGYX(NZ,NY,NX)
    plt_bgcr%TCO2A(NZ)=TCO2A(NZ,NY,NX)
    plt_ew%CTRAN(NZ)=CTRAN(NZ,NY,NX)

    plt_morph%NBT(NZ)=NBT(NZ,NY,NX)
    plt_morph%NG(NZ)=NG(NZ,NY,NX)
    plt_pheno%IFLGI(NZ)=IFLGI(NZ,NY,NX)
    plt_morph%NIX(NZ)=NIX(NZ,NY,NX)
    plt_morph%NRT(NZ)=NRT(NZ,NY,NX)
    plt_morph%NB1(NZ)=NB1(NZ,NY,NX)
    plt_morph%NBR(NZ)=NBR(NZ,NY,NX)

    plt_pheno%IFLGC(NZ)=IFLGC(NZ,NY,NX)
    plt_distb%IDAY0(NZ)=IDAY0(NZ,NY,NX)
    plt_distb%IDAYH(NZ)=IDAYH(NZ,NY,NX)
    plt_distb%IYRH(NZ)=IYRH(NZ,NY,NX)


    plt_distb%HVST(NZ) =HVST(NZ,I,NY,NX)
    plt_distb%IHVST(NZ)=IHVST(NZ,I,NY,NX)
    plt_distb%JHVST(NZ)=JHVST(NZ,I,NY,NX)
    plt_distb%THIN(NZ) =THIN(NZ,I,NY,NX)
    plt_morph%ARSTP(NZ)=ARSTP(NZ,NY,NX)
    plt_morph%ARLFP(NZ)=ARLFP(NZ,NY,NX)

    plt_photo%O2I(NZ)  =O2I(NZ,NY,NX)
    plt_photo%CO2I(NZ) =CO2I(NZ,NY,NX)

    plt_biom%CCPLNP(NZ)=CCPLNP(NZ,NY,NX)
    plt_bgcr%CNET(NZ)=CNET(NZ,NY,NX)
    plt_bgcr%CARBN(NZ)=CARBN(NZ,NY,NX)

    plt_allom%CNWS(NZ)=CNWS(NZ,NY,NX)
    plt_allom%CPWS(NZ)=CPWS(NZ,NY,NX)

    plt_photo%DCO2(NZ)=DCO2(NZ,NY,NX)
    plt_photo%FMOL(NZ)=FMOL(NZ,NY,NX)
    plt_allom%FNOD(NZ)=FNOD(NZ,NY,NX)

    plt_morph%HTCTL(NZ)=HTCTL(NZ,NY,NX)
    plt_rbgc%HEUPTK(1:npelms,NZ)=HEUPTK(1:npelms,NZ,NY,NX)
    plt_morph%HTSTZ(NZ)=HTSTZ(NZ,NY,NX)
    plt_morph%NI(NZ)   =NI(NZ,NY,NX)
    plt_photo%RA(NZ)   =RA(NZ,NY,NX)
    plt_photo%RC(NZ)   =RC(NZ,NY,NX)
    plt_ew%TKC(NZ)     =TKC(NZ,NY,NX)
    plt_ew%SFLXC(NZ)   =SFLXC(NZ,NY,NX)
    plt_rad%RAD1(NZ)   =RAD1(NZ,NY,NX)
    plt_rad%THRM1(NZ)  =THRM1(NZ,NY,NX)
    plt_ew%EFLXC(NZ)   =EFLXC(NZ,NY,NX)
    plt_ew%EVAPC(NZ)   =EVAPC(NZ,NY,NX)
    plt_photo%RSMN(NZ) =RSMN(NZ,NY,NX)
    plt_pheno%OFFST(NZ)=OFFST(NZ,NY,NX)
    plt_pheno%OSTR(NZ)=OSTR(NZ,NY,NX)
    plt_site%PP(NZ)=PP(NZ,NY,NX)
    plt_ew%PSILT(NZ)=PSILT(NZ,NY,NX)
    plt_ew%PSILG(NZ)=PSILG(NZ,NY,NX)
    plt_photo%RCMX(NZ)=RCMX(NZ,NY,NX)
    plt_bgcr%RCO2Z(NZ)=RCO2Z(NZ,NY,NX)
    plt_bgcr%ROXYZ(NZ)=ROXYZ(NZ,NY,NX)
    plt_bgcr%RCH4Z(NZ)=RCH4Z(NZ,NY,NX)
    plt_bgcr%RN2OZ(NZ)=RN2OZ(NZ,NY,NX)
    plt_bgcr%RNH3Z(NZ)=RNH3Z(NZ,NY,NX)
    plt_bgcr%RH2GZ(NZ)=RH2GZ(NZ,NY,NX)
    plt_ew%RAZ(NZ)=RAZ(NZ,NY,NX)
    plt_photo%SCO2(NZ)=SCO2(NZ,NY,NX)
    plt_morph%SDPTH(NZ)=SDPTH(NZ,NY,NX)
    plt_morph%SDPTHI(NZ)=SDPTHI(NZ,NY,NX)
    plt_morph%SDLG(NZ)=SDLG(NZ,NY,NX)
    plt_morph%SDVL(NZ)=SDVL(NZ,NY,NX)
    plt_morph%SDAR(NZ)=SDAR(NZ,NY,NX)
    plt_pheno%TCZ(NZ)    =TCZ(NZ,NY,NX)
    plt_pheno%TCG(NZ)    =TCG(NZ,NY,NX)
    plt_pheno%TCX(NZ)    =TCX(NZ,NY,NX)
    plt_pheno%TKG(NZ)      =TKG(NZ,NY,NX)
    plt_pheno%TFN3(NZ)  =TFN3(NZ,NY,NX)

    plt_bgcr%TCO2T(NZ)  =TCO2T(NZ,NY,NX)
    plt_photo%XKCO2L(NZ)=XKCO2L(NZ,NY,NX)
    plt_photo%XKCO2O(NZ)=XKCO2O(NZ,NY,NX)
    plt_bgcr%TNH3C(NZ)  =TNH3C(NZ,NY,NX)
    plt_bgcr%TZUPFX(NZ) =TZUPFX(NZ,NY,NX)
    plt_rbgc%UPNF(NZ)=UPNF(NZ,NY,NX)
    plt_rbgc%UPNO3(NZ)=UPNO3(NZ,NY,NX)
    plt_rbgc%UPNH4(NZ)=UPNH4(NZ,NY,NX)
    plt_rbgc%UPH1P(NZ)=UPH1P(NZ,NY,NX)
    plt_rbgc%UPH2P(NZ)=UPH2P(NZ,NY,NX)
    plt_ew%VOLWC(NZ)=VOLWC(NZ,NY,NX)
    plt_ew%VHCPC(NZ)=VHCPC(NZ,NY,NX)
    plt_distb%VCH4F(NZ)=VCH4F(NZ,NY,NX)
    plt_distb%VCO2F(NZ)=VCO2F(NZ,NY,NX)
    plt_distb%VN2OF(NZ)=VN2OF(NZ,NY,NX)
    plt_distb%VNH3F(NZ)=VNH3F(NZ,NY,NX)
    plt_distb%VPO4F(NZ)=VPO4F(NZ,NY,NX)
    plt_ew%VOLWP(NZ) =VOLWP(NZ,NY,NX)
    plt_pheno%WSTR(NZ) =WSTR(NZ,NY,NX)
    plt_biom%WTRVX(NZ) =WTRVX(NZ,NY,NX)
    plt_biom%WTLS(NZ)  =WTLS(NZ,NY,NX)

    plt_biom%ZEROL(NZ) =ZEROL(NZ,NY,NX)
    plt_biom%ZEROP(NZ) =ZEROP(NZ,NY,NX)
    plt_morph%ZC(NZ)   =ZC(NZ,NY,NX)
    plt_bgcr%ZNPP(NZ)  =ZNPP(NZ,NY,NX)


    DO L=1,NL(NY,NX)
      DO K=1,jcplx
        DO N=1,2
          DO NE=1,npelms
            plt_rbgc%RDFOME(NE,N,K,L,NZ)=RDFOME(NE,N,K,L,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DO NE=1,npelms
      DO NB=1,NBR(NZ,NY,NX)
        plt_biom%EPOOL(NB,NE,NZ)=EPOOL(NB,NE,NZ,NY,NX)
        plt_biom%EPOLNB(NB,NE,NZ)=EPOLNB(NB,NE,NZ,NY,NX)
        plt_biom%WTSHTBE(NB,NE,NZ)=WTSHTBE(NB,NE,NZ,NY,NX)
        plt_biom%CEPOLB(NB,NE,NZ)=CEPOLB(NB,NE,NZ,NY,NX)
        plt_biom%WTSHEBE(NB,NE,NZ)=WTSHEBE(NB,NE,NZ,NY,NX)
        plt_biom%WTSTKBE(NB,NE,NZ)=WTSTKBE(NB,NE,NZ,NY,NX)
        plt_biom%WTLFBE(NB,NE,NZ)=WTLFBE(NB,NE,NZ,NY,NX)
        plt_biom%WTRSVBE(NB,NE,NZ)=WTRSVBE(NB,NE,NZ,NY,NX)
        plt_biom%WTHSKBE(NB,NE,NZ)=WTHSKBE(NB,NE,NZ,NY,NX)
        plt_biom%WTGRBE(NB,NE,NZ)=WTGRBE(NB,NE,NZ,NY,NX)
        plt_biom%WTEARBE(NB,NE,NZ)=WTEARBE(NB,NE,NZ,NY,NX)
        plt_biom%WTNDBE(NB,NE,NZ)=WTNDBE(NB,NE,NZ,NY,NX)
        plt_biom%WGSHEXE(NB,NE,NZ)=WGSHEXE(NB,NE,NZ,NY,NX)
      ENDDO
    ENDDO

    DO NB=1,NBR(NZ,NY,NX)
      plt_photo%FDBK(NB,NZ)=FDBK(NB,NZ,NY,NX)
      plt_photo%FDBKX(NB,NZ)=FDBKX(NB,NZ,NY,NX)
      plt_pheno%ATRP(NB,NZ)=ATRP(NB,NZ,NY,NX)
      plt_morph%ARLFZ(NB,NZ)=ARLFZ(NB,NZ,NY,NX)
      plt_morph%ARLFB(NB,NZ)=ARLFB(NB,NZ,NY,NX)

      plt_pheno%DGSTGI(NB,NZ)=DGSTGI(NB,NZ,NY,NX)
      plt_pheno%DGSTGF(NB,NZ)=DGSTGF(NB,NZ,NY,NX)
      plt_pheno%FLG4(NB,NZ)=FLG4(NB,NZ,NY,NX)
      plt_pheno%FLGZ(NB,NZ)=FLGZ(NB,NZ,NY,NX)
      plt_pheno%GROUP(NB,NZ)=GROUP(NB,NZ,NY,NX)
      plt_pheno%GSTGI(NB,NZ)=GSTGI(NB,NZ,NY,NX)
      plt_pheno%GSTGF(NB,NZ)=GSTGF(NB,NZ,NY,NX)
      plt_morph%GRNXB(NB,NZ)=GRNXB(NB,NZ,NY,NX)
      plt_morph%GRNOB(NB,NZ)=GRNOB(NB,NZ,NY,NX)
      plt_allom%GRWTB(NB,NZ)=GRWTB(NB,NZ,NY,NX)
      plt_morph%HTSHEX(NB,NZ)=HTSHEX(NB,NZ,NY,NX)
      plt_pheno%IDTHB(NB,NZ)=IDTHB(NB,NZ,NY,NX)
      plt_pheno%IFLGP(NB,NZ)=IFLGP(NB,NZ,NY,NX)
      plt_pheno%IFLGF(NB,NZ)=IFLGF(NB,NZ,NY,NX)
      plt_pheno%IFLGE(NB,NZ)=IFLGE(NB,NZ,NY,NX)
      plt_pheno%IFLGA(NB,NZ)=IFLGA(NB,NZ,NY,NX)
      plt_pheno%IFLGG(NB,NZ)=IFLGG(NB,NZ,NY,NX)
      plt_pheno%IFLGR(NB,NZ)=IFLGR(NB,NZ,NY,NX)
      plt_pheno%IFLGQ(NB,NZ)=IFLGQ(NB,NZ,NY,NX)
      plt_pheno%KVSTG(NB,NZ)=KVSTG(NB,NZ,NY,NX)
      plt_pheno%KVSTGN(NB,NZ)=KVSTGN(NB,NZ,NY,NX)
      plt_morph%NBTB(NB,NZ)=NBTB(NB,NZ,NY,NX)
      plt_morph%PSTG(NB,NZ)=PSTG(NB,NZ,NY,NX)
      plt_morph%PSTGI(NB,NZ)=PSTGI(NB,NZ,NY,NX)
      plt_morph%PSTGF(NB,NZ)=PSTGF(NB,NZ,NY,NX)
      plt_pheno%RCELX(NB,1:npelms,NZ)=RCELX(NB,1:npelms,NZ,NY,NX)
      plt_pheno%RCESX(NB,1:npelms,NZ)=RCESX(NB,1:npelms,NZ,NY,NX)
      plt_rbgc%RNH3B(NB,NZ)=RNH3B(NB,NZ,NY,NX)
      plt_pheno%TGSTGI(NB,NZ)=TGSTGI(NB,NZ,NY,NX)
      plt_pheno%TGSTGF(NB,NZ)=TGSTGF(NB,NZ,NY,NX)
      plt_pheno%VSTGX(NB,NZ)=VSTGX(NB,NZ,NY,NX)
      plt_morph%VSTG(NB,NZ)=VSTG(NB,NZ,NY,NX)
      plt_pheno%VRNY(NB,NZ)=VRNY(NB,NZ,NY,NX)
      plt_pheno%VRNZ(NB,NZ)=VRNZ(NB,NZ,NY,NX)
      plt_pheno%VRNS(NB,NZ)=VRNS(NB,NZ,NY,NX)
      plt_pheno%VRNF(NB,NZ)=VRNF(NB,NZ,NY,NX)
      plt_biom%WTLSB(NB,NZ)=WTLSB(NB,NZ,NY,NX)
      plt_biom%WGLFEX(NB,1:npelms,NZ) =WGLFEX(NB,1:npelms,NZ,NY,NX)
      plt_biom%WTSTXBE(NB,1:npelms,NZ)=WTSTXBE(NB,1:npelms,NZ,NY,NX)
      plt_biom%WVSTKB(NB,NZ)=WVSTKB(NB,NZ,NY,NX)
      DO M=1,pltpar%jpstgs
        plt_pheno%IDAY(M,NB,NZ)=IDAY(M,NB,NZ,NY,NX)
      ENDDO
      DO K=1,JNODS
        DO  L=1,JC
          DO N=1,JLI
            plt_morph%SURF(N,L,K,NB,NZ) =SURF(N,L,K,NB,NZ,NY,NX)
            plt_photo%SURFX(N,L,K,NB,NZ)=SURFX(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
        plt_photo%CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ,NY,NX)
        plt_photo%CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ,NY,NX)
        plt_photo%CO2B(K,NB,NZ)  =CO2B(K,NB,NZ,NY,NX)
        plt_photo%COMPL(K,NB,NZ) =COMPL(K,NB,NZ,NY,NX)
        plt_photo%CBXN(K,NB,NZ)  =CBXN(K,NB,NZ,NY,NX)
        plt_photo%CBXN4(K,NB,NZ) =CBXN4(K,NB,NZ,NY,NX)
        plt_photo%ETGRO(K,NB,NZ) =ETGRO(K,NB,NZ,NY,NX)
        plt_photo%ETGR4(K,NB,NZ) =ETGR4(K,NB,NZ,NY,NX)
        plt_photo%FDBK4(K,NB,NZ) =FDBK4(K,NB,NZ,NY,NX)
        plt_photo%HCOB(K,NB,NZ)  =HCOB(K,NB,NZ,NY,NX)
        plt_photo%VCGRO(K,NB,NZ) =VCGRO(K,NB,NZ,NY,NX)
        plt_photo%VGRO(K,NB,NZ)  =VGRO(K,NB,NZ,NY,NX)
        plt_photo%VCGR4(K,NB,NZ) =VCGR4(K,NB,NZ,NY,NX)
        plt_photo%VGRO4(K,NB,NZ) =VGRO4(K,NB,NZ,NY,NX)
      ENDDO
      DO K=0,JNODS
        plt_morph%ARLF1(K,NB,NZ) =ARLF(K,NB,NZ,NY,NX)
        plt_morph%HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ,NY,NX)
        plt_morph%HTSHE(K,NB,NZ) =HTSHE(K,NB,NZ,NY,NX)
        plt_morph%HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ,NY,NX)
        plt_biom%WGNODE(K,NB,1:npelms,NZ) =WGNODE(K,NB,1:npelms,NZ,NY,NX)
        plt_biom%WGLFE(K,NB,1:npelms,NZ)   =WGLFE(K,NB,1:npelms,NZ,NY,NX)
        plt_biom%WSLF(K,NB,NZ)   =WSLF(K,NB,NZ,NY,NX)
        plt_biom%WGSHE(K,NB,1:npelms,NZ)  =WGSHE(K,NB,1:npelms,NZ,NY,NX)
        plt_biom%WSSHE(K,NB,NZ)  =WSSHE(K,NB,NZ,NY,NX)
      ENDDO

      DO K=0,JNODS
        DO  L=1,JC
          plt_morph%ARLFL(L,K,NB,NZ)=ARLFL(L,K,NB,NZ,NY,NX)
          plt_biom%WGLFLE(L,K,NB,1:npelms,NZ) =WGLFLE(L,K,NB,1:npelms,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,JC
        plt_morph%ARSTK(L,NB,NZ)=ARSTK(L,NB,NZ,NY,NX)
      ENDDO
    enddo
    DO NE=1,npelms
      DO L=1,NL(NY,NX)
        plt_biom%WTNDLE(L,NE,NZ) =WTNDLE(L,NE,NZ,NY,NX)
      ENDDO
    ENDDO
    DO L=1,NL(NY,NX)
      DO N=1,MY(NZ,NY,NX)
        plt_biom%EPOOLR(1:npelms,N,L,NZ)=EPOOLR(1:npelms,N,L,NZ,NY,NX)
        plt_rbgc%trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=trcs_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)
        plt_rbgc%trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=trcg_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)

        plt_biom%CCPOLR(N,L,NZ)=CCPOLR(N,L,NZ,NY,NX)
        plt_biom%CZPOLR(N,L,NZ)=CZPOLR(N,L,NZ,NY,NX)
        plt_biom%CPPOLR(N,L,NZ)=CPPOLR(N,L,NZ,NY,NX)
        plt_biom%CWSRTL(N,L,NZ)=CWSRTL(N,L,NZ,NY,NX)

        plt_ew%PSIRT(N,L,NZ)=PSIRT(N,L,NZ,NY,NX)
        plt_ew%PSIRO(N,L,NZ)=PSIRO(N,L,NZ,NY,NX)
        plt_ew%PSIRG(N,L,NZ)=PSIRG(N,L,NZ,NY,NX)
        plt_rbgc%RUPNHB(N,L,NZ)=RUPNHB(N,L,NZ,NY,NX)
        plt_rbgc%RUPNH4(N,L,NZ)=RUPNH4(N,L,NZ,NY,NX)
        plt_rbgc%RUPH2P(N,L,NZ)=RUPH2P(N,L,NZ,NY,NX)
        plt_rbgc%RUPNOB(N,L,NZ)=RUPNOB(N,L,NZ,NY,NX)
        plt_rbgc%RUPNO3(N,L,NZ)=RUPNO3(N,L,NZ,NY,NX)
        plt_rbgc%RUPH1B(N,L,NZ)=RUPH1B(N,L,NZ,NY,NX)
        plt_rbgc%RUPH1P(N,L,NZ)=RUPH1P(N,L,NZ,NY,NX)
        plt_rbgc%RUPH2B(N,L,NZ)=RUPH2B(N,L,NZ,NY,NX)
        plt_rbgc%RUONHB(N,L,NZ)=RUONHB(N,L,NZ,NY,NX)
        plt_rbgc%RUONH4(N,L,NZ)=RUONH4(N,L,NZ,NY,NX)
        plt_rbgc%RUOH2P(N,L,NZ)=RUOH2P(N,L,NZ,NY,NX)
        plt_rbgc%RUONOB(N,L,NZ)=RUONOB(N,L,NZ,NY,NX)
        plt_rbgc%RUONO3(N,L,NZ)=RUONO3(N,L,NZ,NY,NX)
        plt_rbgc%RUOH1B(N,L,NZ)=RUOH1B(N,L,NZ,NY,NX)
        plt_rbgc%RUOH1P(N,L,NZ)=RUOH1P(N,L,NZ,NY,NX)
        plt_rbgc%RUOH2B(N,L,NZ)=RUOH2B(N,L,NZ,NY,NX)
        plt_rbgc%RUCNHB(N,L,NZ)=RUCNHB(N,L,NZ,NY,NX)
        plt_rbgc%RUCNH4(N,L,NZ)=RUCNH4(N,L,NZ,NY,NX)
        plt_rbgc%RUCH2P(N,L,NZ)=RUCH2P(N,L,NZ,NY,NX)
        plt_rbgc%RUCNOB(N,L,NZ)=RUCNOB(N,L,NZ,NY,NX)
        plt_rbgc%RUCNO3(N,L,NZ)=RUCNO3(N,L,NZ,NY,NX)
        plt_rbgc%RUCH1B(N,L,NZ)=RUCH1B(N,L,NZ,NY,NX)
        plt_rbgc%RUCH1P(N,L,NZ)=RUCH1P(N,L,NZ,NY,NX)
        plt_rbgc%RUCH2B(N,L,NZ)=RUCH2B(N,L,NZ,NY,NX)
        plt_rbgc%RCO2M(N,L,NZ)=RCO2M(N,L,NZ,NY,NX)
        plt_rbgc%RCO2N(N,L,NZ)=RCO2N(N,L,NZ,NY,NX)
        plt_rbgc%RCO2A(N,L,NZ)=RCO2A(N,L,NZ,NY,NX)
        plt_morph%RTN1(N,L,NZ)=RTN1(N,L,NZ,NY,NX)
        plt_morph%RTNL(N,L,NZ)=RTNL(N,L,NZ,NY,NX)
        plt_morph%RTLGP(N,L,NZ)=RTLGP(N,L,NZ,NY,NX)
        plt_morph%RTDNP(N,L,NZ)=RTDNP(N,L,NZ,NY,NX)
        plt_morph%RTVLP(N,L,NZ)=RTVLP(N,L,NZ,NY,NX)
        plt_morph%RTVLW(N,L,NZ)=RTVLW(N,L,NZ,NY,NX)
        plt_morph%RRAD1(N,L,NZ)=RRAD1(N,L,NZ,NY,NX)
        plt_morph%RRAD2(N,L,NZ)=RRAD2(N,L,NZ,NY,NX)
        plt_morph%RTARP(N,L,NZ)=RTARP(N,L,NZ,NY,NX)
        plt_morph%RTLGA(N,L,NZ)=RTLGA(N,L,NZ,NY,NX)
        plt_rbgc%RCO2P(N,L,NZ)=RCO2P(N,L,NZ,NY,NX)
        plt_rbgc%RUPOXP(N,L,NZ)=RUPOXP(N,L,NZ,NY,NX)
        plt_rbgc%RCO2S(N,L,NZ)=RCO2S(N,L,NZ,NY,NX)
        plt_rbgc%RUPOXS(N,L,NZ)=RUPOXS(N,L,NZ,NY,NX)
        plt_rbgc%RUPCHS(N,L,NZ)=RUPCHS(N,L,NZ,NY,NX)
        plt_rbgc%RUPN2S(N,L,NZ)=RUPN2S(N,L,NZ,NY,NX)
        plt_rbgc%RUPN3S(N,L,NZ)=RUPN3S(N,L,NZ,NY,NX)
        plt_rbgc%RUPN3B(N,L,NZ)=RUPN3B(N,L,NZ,NY,NX)
        plt_rbgc%RUPHGS(N,L,NZ)=RUPHGS(N,L,NZ,NY,NX)
        plt_rbgc%RCOFLA(N,L,NZ)=RCOFLA(N,L,NZ,NY,NX)
        plt_rbgc%ROXFLA(N,L,NZ)=ROXFLA(N,L,NZ,NY,NX)
        plt_rbgc%RCHFLA(N,L,NZ)=RCHFLA(N,L,NZ,NY,NX)
        plt_rbgc%RN2FLA(N,L,NZ)=RN2FLA(N,L,NZ,NY,NX)
        plt_rbgc%RNHFLA(N,L,NZ)=RNHFLA(N,L,NZ,NY,NX)
        plt_rbgc%RHGFLA(N,L,NZ)=RHGFLA(N,L,NZ,NY,NX)
        plt_rbgc%RCODFA(N,L,NZ)=RCODFA(N,L,NZ,NY,NX)
        plt_rbgc%ROXDFA(N,L,NZ)=ROXDFA(N,L,NZ,NY,NX)
        plt_rbgc%RCHDFA(N,L,NZ)=RCHDFA(N,L,NZ,NY,NX)
        plt_rbgc%RN2DFA(N,L,NZ)=RN2DFA(N,L,NZ,NY,NX)
        plt_rbgc%RNHDFA(N,L,NZ)=RNHDFA(N,L,NZ,NY,NX)
        plt_rbgc%RHGDFA(N,L,NZ)=RHGDFA(N,L,NZ,NY,NX)
        plt_rbgc%ROXYP(N,L,NZ) =ROXYP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNHP(N,L,NZ)=RUNNHP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNBP(N,L,NZ)=RUNNBP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNOP(N,L,NZ)=RUNNOP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNXP(N,L,NZ)=RUNNXP(N,L,NZ,NY,NX)
        plt_rbgc%RUPP2P(N,L,NZ)=RUPP2P(N,L,NZ,NY,NX)
        plt_rbgc%RUPP2B(N,L,NZ)=RUPP2B(N,L,NZ,NY,NX)
        plt_rbgc%RUPP1P(N,L,NZ)=RUPP1P(N,L,NZ,NY,NX)
        plt_rbgc%RUPP1B(N,L,NZ)=RUPP1B(N,L,NZ,NY,NX)
        plt_ew%UPWTR(N,L,NZ)   =UPWTR(N,L,NZ,NY,NX)
        plt_rbgc%WFR(N,L,NZ)   =WFR(N,L,NZ,NY,NX)
        plt_biom%WTRTL(N,L,NZ) =WTRTL(N,L,NZ,NY,NX)
        plt_biom%WTRTD(N,L,NZ) =WTRTD(N,L,NZ,NY,NX)
        plt_biom%WSRTL(N,L,NZ) =WSRTL(N,L,NZ,NY,NX)

      enddo
      plt_bgcr%RUPNF(L,NZ) =RUPNF(L,NZ,NY,NX)
      plt_pheno%TFN4(L,NZ) =TFN4(L,NZ,NY,NX)
      plt_biom%EPOOLN(L,1:npelms,NZ)=EPOOLN(L,1:npelms,NZ,NY,NX)
    ENDDO
    DO L=1,JC
      plt_morph%ARSTV(L,NZ)=ARSTV(L,NZ,NY,NX)
      plt_morph%ARLFV(L,NZ)=ARLFV(L,NZ,NY,NX)
      plt_biom%WGLFV(L,NZ) =WGLFV(L,NZ,NY,NX)
    ENDDO
    DO N=1,MY(NZ,NY,NX)
      plt_morph%DMVL(N,NZ)=DMVL(N,NZ,NY,NX)
      plt_morph%PORTX(N,NZ)=PORTX(N,NZ,NY,NX)
      plt_morph%RTAR2X(N,NZ)=RTAR2X(N,NZ,NY,NX)
      plt_morph%RTAR1X(N,NZ)=RTAR1X(N,NZ,NY,NX)
      plt_morph%RRAD1X(N,NZ)=RRAD1X(N,NZ,NY,NX)
      plt_morph%RRAD2X(N,NZ)=RRAD2X(N,NZ,NY,NX)
      plt_morph%RTLG1X(N,NZ)=RTLG1X(N,NZ,NY,NX)
      plt_morph%RRADP(N,NZ) =RRADP(N,NZ,NY,NX)
      plt_morph%RTLG2X(N,NZ)=RTLG2X(N,NZ,NY,NX)
    ENDDO

    DO NR=1,NRT(NZ,NY,NX)
      plt_morph%NINR(NR,NZ)=NINR(NR,NZ,NY,NX)
      DO L=1,NJ(NY,NX)
        DO N=1,MY(NZ,NY,NX)
          plt_morph%RTLG1(N,L,NR,NZ)=RTLG1(N,L,NR,NZ,NY,NX)
          plt_morph%RTLG2(N,L,NR,NZ)=RTLG2(N,L,NR,NZ,NY,NX)
          plt_morph%RTN2(N,L,NR,NZ) =RTN2(N,L,NR,NZ,NY,NX)
          plt_biom%WTRT2E(1:npelms,N,L,NR,NZ) =WTRT2E(1:npelms,N,L,NR,NZ,NY,NX)
          plt_biom%WTRT1E(1:npelms,N,L,NR,NZ) =WTRT1E(1:npelms,N,L,NR,NZ,NY,NX)
        enddo
      enddo
      DO N=1,MY(NZ,NY,NX)
        plt_morph%RTDP1(N,NR,NZ)=RTDP1(N,NR,NZ,NY,NX)
        plt_biom%RTWT1E(N,NR,1:npelms,NZ)=RTWT1E(N,NR,1:npelms,NZ,NY,NX)
      enddo
    enddo
    DO NE=1,npelms
      DO M=1,jsken
        plt_biom%WTSTDE(M,NE,NZ)=WTSTDE(M,NE,NZ,NY,NX)
      ENDDO
    ENDDO
    DO L=0,NJ(NY,NX)
      DO K=1,micpar%n_pltlitrk
        DO M=1,jsken
          plt_bgcr%ESNC(M,1:npelms,K,L,NZ)=ESNC(M,1:npelms,K,L,NZ,NY,NX)
        enddo
      enddo
    ENDDO
    DO NE=1,npelms
      DO M=1,jsken
        DO N=0,JP
          plt_soilchem%CFOPE(N,M,NE,NZ)=CFOPE(N,M,NE,NZ,NY,NX)
        enddo
      enddo
    ENDDO
  ENDDO

  DO L=1,NL(NY,NX)
    DO M=1,NPH
      plt_rbgc%ROXSK(M,L)=ROXSK(M,L,NY,NX)
      plt_soilchem%THETPM(M,L)=THETPM(M,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine PlantAPISend


!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMSend(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,N,M,NN,NZ,K,NB
!  Integers
  plt_site%IETYP=IETYP(NY,NX)
  plt_site%NP=NP(NY,NX)
  plt_site%NU=NU(NY,NX)
  plt_morph%ARSTC=ARSTC(NY,NX)
  plt_morph%ARLFC=ARLFC(NY,NX)
  plt_site%ZEROS =ZEROS(NY,NX)
  plt_site%ZERO  =ZERO
  plt_ew%DPTHS   =DPTHS(NY,NX)
  plt_ew%TKA     =TKA(NY,NX)
  plt_morph%ZT   =ZT(NY,NX)
  plt_site%Z0=Z0(NY,NX)
  plt_ew%ZD=ZD(NY,NX)
  plt_site%UA=UA(NY,NX)
  plt_ew%VHCPWX=VHCPWX(NY,NX)
  plt_ew%VHCPW1=VHCPW(1,NY,NX)
  plt_morph%ZL(0)=ZL(0,NY,NX)
  DO L=1,JC
    plt_morph%ARSTT(L)=ARSTT(L,NY,NX)
    plt_morph%ARLFT(L)=ARLFT(L,NY,NX)
    plt_morph%ZL(L)=ZL(L,NY,NX)
    plt_rad%TAUS(L)=TAUS(L,NY,NX)
  ENDDO
  plt_rad%TAUS(JC+1)=TAUS(JC+1,NY,NX)

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)=AREA(3,L,NY,NX)
    plt_soilchem%VOLX(L)=VOLX(L,NY,NX)
    plt_soilchem%VOLY(L)=VOLY(L,NY,NX)
    plt_soilchem%VOLW(L)=VOLW(L,NY,NX)
  ENDDO
  plt_rad%SSIN=SSIN(NY,NX)
  plt_site%ZNOON=ZNOON(NY,NX)
  plt_morph%ARLSS=ARLSS(NY,NX)
  plt_rad%GAZI=GAZI(NY,NX)
  plt_rad%GCOS=GCOS(NY,NX)
  plt_rad%GSIN=GSIN(NY,NX)
  plt_rad%RADY=RADY(NY,NX)
  plt_rad%RAPY=RAPY(NY,NX)
  plt_rad%RADS=RADS(NY,NX)
  plt_rad%RAPS=RAPS(NY,NX)
  plt_site%ZS=ZS(NY,NX)
  plt_ew%VOLWS=VOLWS(NY,NX)
  plt_ew%VOLIS=VOLIS(NY,NX)
  plt_ew%VOLSS=VOLSS(NY,NX)
  plt_rad%TYSIN  =TYSIN
  plt_rad%ALBS   =ALBS(NY,NX)
  plt_rad%ALBX   =ALBX(NY,NX)
  plt_site%ZEROS2=ZEROS2(NY,NX)
  plt_site%POROS1=POROS(NU(NY,NX),NY,NX)
  DO NZ=1,NP(NY,NX)
    plt_morph%ARLFP(NZ)=ARLFP(NZ,NY,NX)
    plt_morph%ZC(NZ)=ZC(NZ,NY,NX)
    plt_morph%CFX(NZ)=CFX(NZ,NY,NX)
    plt_rad%ABSR(NZ)=ABSR(NZ,NY,NX)
    plt_rad%ABSP(NZ)=ABSP(NZ,NY,NX)
    plt_rad%TAUR(NZ)=TAUR(NZ,NY,NX)
    plt_rad%ALBR(NZ)=ALBR(NZ,NY,NX)
    plt_rad%TAUP(NZ)=TAUP(NZ,NY,NX)
    plt_rad%ALBP(NZ)=ALBP(NZ,NY,NX)
    plt_morph%NBR(NZ)=NBR(NZ,NY,NX)
    plt_morph%CF(NZ)=CF(NZ,NY,NX)
    DO NB=1,NBR(NZ,NY,NX)
      DO K=0,JNODS
        DO  L=1,JC
          plt_morph%ARLFL(L,K,NB,NZ)=ARLFL(L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,JC
        plt_morph%ARSTK(L,NB,NZ)=ARSTK(L,NB,NZ,NY,NX)
      ENDDO
      DO K=1,JNODS
        DO  L=1,JC
          DO N=1,JLI
            plt_morph%SURF(N,L,K,NB,NZ)=SURF(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      DO  L=1,JC
        DO N=1,JLI
          plt_morph%SURFB(N,L,NB,NZ)=SURFB(N,L,NB,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO N=1,JSA
    plt_rad%OMEGAG(N)=OMEGAG(N,NY,NX)
  ENDDO
  DO N=1,JLI
    plt_rad%ZCOS(N)=ZCOS(N)
    plt_rad%ZSIN(N)=ZSIN(N)
  ENDDO
  DO NN=1,JLA
    DO M=1,JLI
      DO N=1,JSA
        plt_rad%OMEGA(N,M,NN)=OMEGA(N,M,NN)
        plt_rad%OMEGX(N,M,NN)=OMEGX(N,M,NN)
        plt_rad%IALBY(N,M,NN)=IALBY(N,M,NN)
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMSend

!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMRecv(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer :: N,M,NN,L,NZ,K,NB

  ZD(NY,NX)=plt_ew%ZD
  ZR(NY,NX)=plt_ew%ZR
  RAB(NY,NX)=plt_ew%RAB
  RIB(NY,NX)=plt_ew%RIB
  ZT(NY,NX)=plt_morph%ZT
  RADS(NY,NX)=plt_rad%RADS
  RADY(NY,NX)=plt_rad%RADY
  RAPS(NY,NX)=plt_rad%RAPS
  RAPY(NY,NX)=plt_rad%RAPY
  RADG(NY,NX)=plt_rad%RADG
  FRADG(NY,NX)=plt_rad%FRADG
  RAD(NY,NX)=plt_rad%RAD0
  RAP(NY,NX)=plt_rad%RAP0
  DO L=0,JC
    ZL(L,NY,NX)=plt_morph%ZL(L)
  ENDDO
  DO L=1,JC
    TAUS(L,NY,NX)=plt_rad%TAUS(L)
    TAU0(L,NY,NX)=plt_rad%TAU0(L)
  ENDDO
  ARLSS(NY,NX)=plt_morph%ARLSS
  DO NZ=1,NP(NY,NX)
    ARLFS(NZ,NY,NX)=plt_morph%ARLFS(NZ)
    RADC(NZ,NY,NX) =plt_rad%RADC(NZ)
    RADP(NZ,NY,NX) =plt_rad%RADP(NZ)
    CFX(NZ,NY,NX)  =plt_morph%CFX(NZ)
    FRADP(NZ,NY,NX)=plt_rad%FRADP(NZ)

    DO L=1,JC
      DO M=1,JSA
        DO  N=1,JLI
          PARDIF(N,M,L,NZ,NY,NX)=plt_rad%PARDIF(N,M,L,NZ)
          PAR(N,M,L,NZ,NY,NX)   =plt_rad%PAR(N,M,L,NZ)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMRecv

end module PlantAPI
