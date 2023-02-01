module RestartMod
!
! DESCRIPTION
! code to read/write restart files
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcoSIMConfig, only : jcplx=> jcplxc, NFGs=> NFGsc,nlbiomcp=>nlbiomcpc
  use EcoSIMConfig, only : ndbiomcp=>ndbiomcpc,jsken=>jskenc
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use CanopyRadDataType
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
  private
  save
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: restart
  contains

  subroutine restart(I,NHW,NHE,NVN,NVS)
  use EcoSIMCtrlDataType, only : lverb
  implicit none
  integer, intent (in) :: I,NHW,NHE,NVN,NVS

  if(lverb)WRITE(*,334)'WOUTS'
  CALL WOUTS(I,NHW,NHE,NVN,NVS)
  if(lverb)WRITE(*,334)'WOUTP'
  CALL WOUTP(I,NHW,NHE,NVN,NVS)
  if(lverb)WRITE(*,334)'WOUTQ'
  CALL WOUTQ(I,NHW,NHE,NVN,NVS)

334   FORMAT(A8)

  end subroutine restart

!------------------------------------------------------------------------------------------
  SUBROUTINE woutp(I,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES OUT ALL PLANT VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: K,L,NX,NY,NZ,N,NB,NR,M
!     execution begins here

  D9990: DO NX=NHW,NHE
    D9985: DO NY=NVN,NVS
      WRITE(26,90)I,IDATA(3),NP(NY,NX),IFLGT(NY,NX) &
        ,(IFLGC(NZ,NY,NX),NZ=1,NP(NY,NX)) &
        ,(IDTH(NZ,NY,NX),NZ=1,NP(NY,NX))
      WRITE(26,97)I,IDATA(3),NP(NY,NX) &
        ,((WTSTDE(M,ielmc,NZ,NY,NX),M=1,jsken),NZ=1,NP(NY,NX)) &
        ,((WTSTDE(M,ielmn,NZ,NY,NX),M=1,jsken),NZ=1,NP(NY,NX)) &
        ,((WTSTDE(M,ielmp,NZ,NY,NX),M=1,jsken),NZ=1,NP(NY,NX))
      IF(IFLGT(NY,NX).GT.0)THEN
        D9980: DO NZ=1,NP(NY,NX)
          IF(IFLGC(NZ,NY,NX).NE.0)THEN
            WRITE(26,91)I,IDATA(3),NZ,IYR0(NZ,NY,NX),IDAY0(NZ,NY,NX) &
              ,IYRH(NZ,NY,NX),IDAYH(NZ,NY,NX),NG(NZ,NY,NX),IDTHP(NZ,NY,NX) &
              ,IDTHR(NZ,NY,NX),NBR(NZ,NY,NX),NBT(NZ,NY,NX),NB1(NZ,NY,NX) &
              ,IFLGI(NZ,NY,NX),NRT(NZ,NY,NX),NIX(NZ,NY,NX),MY(NZ,NY,NX) &
              ,(NINR(NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
            WRITE(26,93)I,IDATA(3),NZ,TEUPTK(ielmc,NZ,NY,NX),TESNC(ielmc,NZ,NY,NX) &
              ,TEUPTK(ielmn,NZ,NY,NX),TESNC(ielmn,NZ,NY,NX),TEUPTK(ielmp,NZ,NY,NX),TESNC(ielmp,NZ,NY,NX) &
              ,TZUPFX(NZ,NY,NX),TCO2T(NZ,NY,NX),CTRAN(NZ,NY,NX),CARBN(NZ,NY,NX) &
              ,TCC(NZ,NY,NX),TKC(NZ,NY,NX),TCG(NZ,NY,NX),TKG(NZ,NY,NX) &
              ,TFN3(NZ,NY,NX),WVSTK(NZ,NY,NX),VOLWP(NZ,NY,NX) &
              ,PSILT(NZ,NY,NX),PSILO(NZ,NY,NX),PSILG(NZ,NY,NX),WTRTA(NZ,NY,NX) &
              ,ARLFP(NZ,NY,NX),ARSTP(NZ,NY,NX),PP(NZ,NY,NX),HVSTE(ielmc,NZ,NY,NX) &
              ,HVSTE(ielmn,NZ,NY,NX),HVSTE(ielmp,NZ,NY,NX),THVSTE(ielmc,NZ,NY,NX),THVSTE(ielmn,NZ,NY,NX) &
              ,THVSTE(ielmp,NZ,NY,NX),TCO2A(NZ,NY,NX),RSETE(ielmc,NZ,NY,NX),RSETE(ielmn,NZ,NY,NX) &
              ,RSETE(ielmp,NZ,NY,NX),TNH3C(NZ,NY,NX),TESN0(ielmc,NZ,NY,NX),PPI(NZ,NY,NX) &
              ,PPX(NZ,NY,NX),WTSTGE(ielmc,NZ,NY,NX),WTSTGE(ielmn,NZ,NY,NX),WTSTGE(ielmp,NZ,NY,NX) &
              ,ZC(NZ,NY,NX),VOLWC(NZ,NY,NX),CF(NZ,NY,NX),VCO2F(NZ,NY,NX) &
              ,VCH4F(NZ,NY,NX),VOXYF(NZ,NY,NX),VNH3F(NZ,NY,NX),VN2OF(NZ,NY,NX) &
              ,VPO4F(NZ,NY,NX)
            WRITE(26,92)I,IDATA(3),NZ,(GROUP(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(PSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(PSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(PSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(VSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(VSTGX(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(GSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(GSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(TGSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(TGSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(VRNS(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(VRNF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(VRNY(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(VRNZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(ATRP(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(FLG4(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,IDATA(3),NZ,(FLGZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(IFLGA(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(IFLGE(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(IFLGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(IFLGR(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(IFLGQ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(IDTHB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(NBTB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(KLEAF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(KVSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,IDATA(3),NZ,(KVSTGN(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            DO N=1,jpstgs
              WRITE(26,95)I,IDATA(3),NZ,(IDAY(N,NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            ENDDO
            WRITE(27,96)I,IDATA(3),NZ,EPOOLP(ielmc,NZ,NY,NX),EPOLNP(ielmc,NZ,NY,NX) &
              ,VHCPC(NZ,NY,NX),DTKC(NZ,NY,NX),WTSHTE(ielmc,NZ,NY,NX) &
              ,WTLS(NZ,NY,NX),WTRTE(ielmc,NZ,NY,NX),WTSHTE(ielmn,NZ,NY,NX),WTSHTE(ielmp,NZ,NY,NX) &
              ,WTLFE(ielmc,NZ,NY,NX),WTSHEE(ielmc,NZ,NY,NX),WTSTKE(ielmc,NZ,NY,NX),WTRSVE(ielmc,NZ,NY,NX) &
              ,WTHSKE(ielmc,NZ,NY,NX),WTEARE(ielmc,NZ,NY,NX),WTGRE(ielmc,NZ,NY,NX),WTNDE(ielmc,NZ,NY,NX) &
              ,WTRVE(ielmc,NZ,NY,NX),WTRVE(ielmn,NZ,NY,NX),WTRVE(ielmp,NZ,NY,NX),HTCTL(NZ,NY,NX) &
              ,SDPTH(NZ,NY,NX),WSTR(NZ,NY,NX) &
              ,CHILL(NZ,NY,NX),WTRTSE(ielmc,NZ,NY,NX),FRADP(NZ,NY,NX)
            WRITE(27,92)I,IDATA(3),NZ,(EPOOL(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(EPOOL(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(EPOOL(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(EPOLNB(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(EPOLNB(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(EPOLNB(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSHTBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTLFBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTNDBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSHEBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSTKBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WVSTKB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTRSVBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTHSKBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTEARBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTGRBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTLSB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSHTBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTLFBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTNDBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSHEBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSTKBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTRSVBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTHSKBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTEARBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTGRBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSHTBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTLFBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTNDBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSHEBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSTKBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTRSVBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTHSKBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTEARBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTGRBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(GRNXB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(GRNOB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(GRWTB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(ARLFB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WGLFEX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WGLFEX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WGLFEX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(ARLFZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(RCELX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(RCELX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(RCELX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WGSHEXE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WGSHEXE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WGSHEXE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(HTSHEX(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(RCESX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(RCESX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(RCESX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSTXBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSTXBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,IDATA(3),NZ,(WTSTXBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,94)I,IDATA(3),NZ,(ARLFV(L,NZ,NY,NX),L=1,JC)
            WRITE(27,94)I,IDATA(3),NZ,(ARSTV(L,NZ,NY,NX),L=1,JC)

            D9945: DO NB=1,NBR(NZ,NY,NX)
              WRITE(28,93)I,IDATA(3),NZ,(CPOOL4(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(CPOOL3(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(CO2B(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(HCOB(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(ARLF(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGLFE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WSLF(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(HTSHE(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGSHE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WSSHE(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(HTNODE(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(HTNODX(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGNODE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGLFE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGSHE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGNODE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGLFE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGSHE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,IDATA(3),NZ,(WGNODE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              D9950: DO K=0,JNODS
                WRITE(28,94)I,IDATA(3),NZ,(ARLFL(L,K,NB,NZ,NY,NX),L=1,JC)
                WRITE(28,94)I,IDATA(3),NZ,(WGLFLE(L,K,NB,ielmc,NZ,NY,NX),L=1,JC)
                WRITE(28,94)I,IDATA(3),NZ,(WGLFLE(L,K,NB,ielmn,NZ,NY,NX),L=1,JC)
                WRITE(28,94)I,IDATA(3),NZ,(WGLFLE(L,K,NB,ielmp,NZ,NY,NX),L=1,JC)
                IF(K.NE.0)THEN
                  D9940: DO N=1,JLI
                    WRITE(28,94)I,IDATA(3),NZ,(SURF(N,L,K,NB,NZ,NY,NX),L=1,JC)
                  ENDDO D9940
                ENDIF
              ENDDO D9950
              WRITE(28,94)I,IDATA(3),NZ,(ARSTK(L,NB,NZ,NY,NX),L=1,JC)
              DO  L=1,JC
                WRITE(28,92)I,IDATA(3),NZ,(SURFB(N,L,NB,NZ,NY,NX),N=1,JLI)
              ENDDO
            ENDDO D9945
            D9970: DO N=1,MY(NZ,NY,NX)
              WRITE(29,94)I,IDATA(3),NZ,(PSIRT(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(PSIRO(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(PSIRG(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTN1(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTNL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTLGP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTDNP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTVLP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTVLW(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RRAD1(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RRAD2(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTARP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTLGA(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RCO2M(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RCO2A(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RCO2N(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcg_rootml(idg_CO2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcg_rootml(idg_O2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcg_rootml(idg_CH4,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcg_rootml(idg_N2O,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcg_rootml(idg_NH3,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcg_rootml(idg_H2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcs_rootml(idg_CO2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcs_rootml(idg_O2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcs_rootml(idg_CH4,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcs_rootml(idg_N2O,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcs_rootml(idg_H2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(trcs_rootml(idg_NH3,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(WTRTL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(WTRTD(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(WSRTL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(ROXYP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUNNHP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUNNOP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUPP2P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUPP1P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUNNBP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUNNXP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUPP2B(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RUPP1B(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(WFR(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(EPOOLR(ielmc,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(EPOOLR(ielmn,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(EPOOLR(ielmp,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTDP1(N,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTWT1E(N,NR,ielmc,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTWT1E(N,NR,ielmn,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              WRITE(29,94)I,IDATA(3),NZ,(RTWT1E(N,NR,ielmp,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              D9965: DO NR=1,NRT(NZ,NY,NX)
                WRITE(29,94)I,IDATA(3),NZ,(RTN2(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(RTLG1(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(WTRT1E(ielmc,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(WTRT1E(ielmn,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(WTRT1E(ielmp,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(WTRT2E(ielmc,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(WTRT2E(ielmn,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(WTRT2E(ielmp,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,IDATA(3),NZ,(RTLG2(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
              ENDDO D9965
            ENDDO D9970
            WRITE(29,94)I,IDATA(3),NZ,(EPOOLN(L,ielmc,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,IDATA(3),NZ,(EPOOLN(L,ielmn,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,IDATA(3),NZ,(EPOOLN(L,ielmp,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,IDATA(3),NZ,(WTNDLE(L,ielmc,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,IDATA(3),NZ,(WTNDLE(L,ielmn,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,IDATA(3),NZ,(WTNDLE(L,ielmp,NZ,NY,NX),L=1,NJ(NY,NX))
          ENDIF
        ENDDO D9980
      ENDIF
    ENDDO D9985
  ENDDO D9990
90    FORMAT(30I4)
91    FORMAT(3I4,4I8,30I4)
92    FORMAT(3I4,10E17.8E3)
93    FORMAT(3I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
94    FORMAT(3I4,21E17.8E3)
95    FORMAT(3I4,10I6)
96    FORMAT(3I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
97    FORMAT(3I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
  RETURN
  END SUBROUTINE woutp

!------------------------------------------------------------------------------------------


  SUBROUTINE woutq(I,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES OUT THE NAMES OF ALL PLANT SPECIES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!
  implicit none
  integer, intent(in) :: I   !day information
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NX,NY,NZ
!     execution begins here
! IFLGC: flag for living pft
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      WRITE(30,90)I,IDATA(3),NP(NY,NX),(DATAP(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NP(NY,NX))
90    FORMAT(2I4,1I3,5(A16,I4))
    ENDDO D9990
  ENDDO D9995
!1000  RETURN
  END SUBROUTINE woutq
!------------------------------------------------------------------------------------------


  SUBROUTINE wouts(I,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE WRITES OUT ALL SOIL VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!

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
      IF(salt_model)THEN
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
      WRITE(22,91)I,IDATA(3),(FertN_soil(ifert_nh4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_soil(ifert_nh3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_soil(ifert_urea,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_soil(ifert_no3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_band(ifert_nh4_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_band(ifert_nh3_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_band(ifert_urea_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,IDATA(3),(FertN_band(ifert_no3_band,L,NY,NX),L=1,NLI(NY,NX))
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

end module restartMod
