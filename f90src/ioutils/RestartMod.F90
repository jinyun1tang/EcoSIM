module RestartMod
!
! DESCRIPTION
! code to read/write restart files
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcoSIMConfig, only : jcplx=> jcplxc, NFGs=> NFGsc,nlbiomcp=>nlbiomcpc
  use EcoSIMConfig, only : ndbiomcp=>ndbiomcpc,jsken=>jskenc
  use EcoSIMCtrlDataType, only : iyear_cur
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
  use EcoSIMCtrlMod, only : lverb
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
      WRITE(26,90)I,iyear_cur,NP(NY,NX),IFLGT(NY,NX) &
        ,(IFLGC(NZ,NY,NX),NZ=1,NP(NY,NX)) &
        ,(IDTH(NZ,NY,NX),NZ=1,NP(NY,NX))
      WRITE(26,97)I,iyear_cur,NP(NY,NX) &
        ,((WTSTDE(M,ielmc,NZ,NY,NX),M=1,jsken),NZ=1,NP(NY,NX)) &
        ,((WTSTDE(M,ielmn,NZ,NY,NX),M=1,jsken),NZ=1,NP(NY,NX)) &
        ,((WTSTDE(M,ielmp,NZ,NY,NX),M=1,jsken),NZ=1,NP(NY,NX))
      IF(IFLGT(NY,NX).GT.0)THEN
        D9980: DO NZ=1,NP(NY,NX)
          IF(IFLGC(NZ,NY,NX).NE.0)THEN
            WRITE(26,91)I,iyear_cur,NZ,IYR0(NZ,NY,NX),IDAY0(NZ,NY,NX) &
              ,IYRH(NZ,NY,NX),IDAYH(NZ,NY,NX),NG(NZ,NY,NX),IDTHP(NZ,NY,NX) &
              ,IDTHR(NZ,NY,NX),NBR(NZ,NY,NX),NBT(NZ,NY,NX),NB1(NZ,NY,NX) &
              ,IFLGI(NZ,NY,NX),NRT(NZ,NY,NX),NIX(NZ,NY,NX),MY(NZ,NY,NX) &
              ,(NINR(NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
            WRITE(26,93)I,iyear_cur,NZ,TEUPTK(ielmc,NZ,NY,NX),TESNC(ielmc,NZ,NY,NX) &
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
            WRITE(26,92)I,iyear_cur,NZ,(GROUP(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(PSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(PSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(PSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(VSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(VSTGX(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(GSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(GSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(TGSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(TGSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(VRNS(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(VRNF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(VRNY(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(VRNZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(ATRP(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(FLG4(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,92)I,iyear_cur,NZ,(FLGZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(IFLGA(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(IFLGE(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(IFLGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(IFLGR(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(IFLGQ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(IDTHB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(NBTB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(KLEAF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(KVSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(26,95)I,iyear_cur,NZ,(KVSTGN(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            DO N=1,jpstgs
              WRITE(26,95)I,iyear_cur,NZ,(IDAY(N,NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            ENDDO
            WRITE(27,96)I,iyear_cur,NZ,EPOOLP(ielmc,NZ,NY,NX),EPOLNP(ielmc,NZ,NY,NX) &
              ,VHCPC(NZ,NY,NX),DTKC(NZ,NY,NX),WTSHTE(ielmc,NZ,NY,NX) &
              ,WTLS(NZ,NY,NX),WTRTE(ielmc,NZ,NY,NX),WTSHTE(ielmn,NZ,NY,NX),WTSHTE(ielmp,NZ,NY,NX) &
              ,WTLFE(ielmc,NZ,NY,NX),WTSHEE(ielmc,NZ,NY,NX),WTSTKE(ielmc,NZ,NY,NX),WTRSVE(ielmc,NZ,NY,NX) &
              ,WTHSKE(ielmc,NZ,NY,NX),WTEARE(ielmc,NZ,NY,NX),WTGRE(ielmc,NZ,NY,NX),WTNDE(ielmc,NZ,NY,NX) &
              ,WTRVE(ielmc,NZ,NY,NX),WTRVE(ielmn,NZ,NY,NX),WTRVE(ielmp,NZ,NY,NX),HTCTL(NZ,NY,NX) &
              ,SDPTH(NZ,NY,NX),WSTR(NZ,NY,NX) &
              ,CHILL(NZ,NY,NX),WTRTSE(ielmc,NZ,NY,NX),FRADP(NZ,NY,NX)
            WRITE(27,92)I,iyear_cur,NZ,(EPOOL(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(EPOOL(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(EPOOL(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(EPOLNB(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(EPOLNB(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(EPOLNB(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSHTBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTLFBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTNDBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSHEBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSTKBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WVSTKB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTRSVBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTHSKBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTEARBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTGRBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTLSB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSHTBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTLFBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTNDBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSHEBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSTKBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTRSVBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTHSKBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTEARBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTGRBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSHTBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTLFBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTNDBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSHEBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSTKBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTRSVBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTHSKBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTEARBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTGRBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(GRNXB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(GRNOB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(GRWTB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(ARLFB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WGLFEX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WGLFEX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WGLFEX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(ARLFZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(RCELX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(RCELX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(RCELX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WGSHEXE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WGSHEXE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WGSHEXE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(HTSHEX(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(RCESX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(RCESX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(RCESX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSTXBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSTXBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,92)I,iyear_cur,NZ,(WTSTXBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            WRITE(27,94)I,iyear_cur,NZ,(ARLFV(L,NZ,NY,NX),L=1,JC)
            WRITE(27,94)I,iyear_cur,NZ,(ARSTV(L,NZ,NY,NX),L=1,JC)

            D9945: DO NB=1,NBR(NZ,NY,NX)
              WRITE(28,93)I,iyear_cur,NZ,(CPOOL4(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(CPOOL3(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(CO2B(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(HCOB(K,NB,NZ,NY,NX),K=1,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(ARLF(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGLFE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WSLF(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(HTSHE(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGSHE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WSSHE(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(HTNODE(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(HTNODX(K,NB,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGNODE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGLFE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGSHE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGNODE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGLFE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGSHE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              WRITE(28,93)I,iyear_cur,NZ,(WGNODE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              D9950: DO K=0,JNODS
                WRITE(28,94)I,iyear_cur,NZ,(ARLFL(L,K,NB,NZ,NY,NX),L=1,JC)
                WRITE(28,94)I,iyear_cur,NZ,(WGLFLE(L,K,NB,ielmc,NZ,NY,NX),L=1,JC)
                WRITE(28,94)I,iyear_cur,NZ,(WGLFLE(L,K,NB,ielmn,NZ,NY,NX),L=1,JC)
                WRITE(28,94)I,iyear_cur,NZ,(WGLFLE(L,K,NB,ielmp,NZ,NY,NX),L=1,JC)
                IF(K.NE.0)THEN
                  D9940: DO N=1,JLI
                    WRITE(28,94)I,iyear_cur,NZ,(SURF(N,L,K,NB,NZ,NY,NX),L=1,JC)
                  ENDDO D9940
                ENDIF
              ENDDO D9950
              WRITE(28,94)I,iyear_cur,NZ,(ARSTK(L,NB,NZ,NY,NX),L=1,JC)
              DO  L=1,JC
                WRITE(28,92)I,iyear_cur,NZ,(SURFB(N,L,NB,NZ,NY,NX),N=1,JLI)
              ENDDO
            ENDDO D9945
            D9970: DO N=1,MY(NZ,NY,NX)
              WRITE(29,94)I,iyear_cur,NZ,(PSIRT(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(PSIRO(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(PSIRG(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTN1(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTNL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTLGP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTDNP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTVLP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTVLW(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RRAD1(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RRAD2(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTARP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTLGA(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RCO2M(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RCO2A(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RCO2N(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcg_rootml(idg_CO2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcg_rootml(idg_O2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcg_rootml(idg_CH4,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcg_rootml(idg_N2O,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcg_rootml(idg_NH3,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcg_rootml(idg_H2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcs_rootml(idg_CO2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcs_rootml(idg_O2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcs_rootml(idg_CH4,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcs_rootml(idg_N2O,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcs_rootml(idg_H2,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(trcs_rootml(idg_NH3,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(WTRTL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(WTRTD(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(WSRTL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(ROXYP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUNNHP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUNNOP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUPP2P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUPP1P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUNNBP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUNNXP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUPP2B(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RUPP1B(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(WFR(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(EPOOLR(ielmc,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(EPOOLR(ielmn,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(EPOOLR(ielmp,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTDP1(N,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTWT1E(N,NR,ielmc,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTWT1E(N,NR,ielmn,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              WRITE(29,94)I,iyear_cur,NZ,(RTWT1E(N,NR,ielmp,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              D9965: DO NR=1,NRT(NZ,NY,NX)
                WRITE(29,94)I,iyear_cur,NZ,(RTN2(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(RTLG1(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(WTRT1E(ielmc,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(WTRT1E(ielmn,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(WTRT1E(ielmp,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(WTRT2E(ielmc,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(WTRT2E(ielmn,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(WTRT2E(ielmp,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                WRITE(29,94)I,iyear_cur,NZ,(RTLG2(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
              ENDDO D9965
            ENDDO D9970
            WRITE(29,94)I,iyear_cur,NZ,(EPOOLN(L,ielmc,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,iyear_cur,NZ,(EPOOLN(L,ielmn,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,iyear_cur,NZ,(EPOOLN(L,ielmp,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,iyear_cur,NZ,(WTNDLE(L,ielmc,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,iyear_cur,NZ,(WTNDLE(L,ielmn,NZ,NY,NX),L=1,NJ(NY,NX))
            WRITE(29,94)I,iyear_cur,NZ,(WTNDLE(L,ielmp,NZ,NY,NX),L=1,NJ(NY,NX))
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
      WRITE(30,90)I,iyear_cur,NP(NY,NX),(DATAP(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NP(NY,NX))
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

  WRITE(21,90)I,iyear_cur,CRAIN,TSEDOU &
    ,HEATIN,OXYGIN,TORGF,TORGN,TORGP,CO2GIN,ZN2GIN,VOLWOU,CEVAP &
    ,CRUN,HEATOU,OXYGOU,TCOU,TZOU,TPOU,TZIN,TPIN,XCSN,XZSN,XPSN
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      WRITE(21,95)I,iyear_cur,(TDTPX(NY,NX,N),N=1,12) &
        ,(TDTPN(NY,NX,N),N=1,12),(TDRAD(NY,NX,N),N=1,12) &
        ,(TDWND(NY,NX,N),N=1,12),(TDHUM(NY,NX,N),N=1,12) &
        ,(TDPRC(NY,NX,N),N=1,12),(TDIRI(NY,NX,N),N=1,12) &
        ,(TDCO2(NY,NX,N),N=1,12),(TDCN4(NY,NX,N),N=1,12) &
        ,(TDCNO(NY,NX,N),N=1,12)
      WRITE(21,93)I,iyear_cur,IFLGT(NY,NX),IFNHB(NY,NX),IDTBL(NY,NX) &
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
      WRITE(21,91)I,iyear_cur,(VOLSSL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(VOLISL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(VOLWSL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(VOLSL(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(DENSS(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(DLYRS(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(VHCPW(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(TKW(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(TCW(L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcg_solsml(idg_CO2,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcg_solsml(idg_CH4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcg_solsml(idg_O2,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcg_solsml(idg_N2,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcg_solsml(idg_N2O,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcn_solsml(ids_NH4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcg_solsml(idg_NH3,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcn_solsml(ids_NO3,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcn_solsml(ids_H1PO4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(trcn_solsml(ids_H2PO4,L,NY,NX),L=1,JS)
      WRITE(21,91)I,iyear_cur,(FHOL(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(DLYR(3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(CDPTH(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(CDPTHZ(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(BKDSI(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(BKDS(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(CORGC(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(POROS(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(FC(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(WP(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(SCNV(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(SCNH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(SAND(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(SILT(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(CLAY(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLW(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLWX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLI(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLP(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLA(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLY(L,NY,NX),L=0,NLI(NY,NX))
!     WRITE(21,91)I,iyear_cur,(GKC4(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,iyear_cur,(GKCH(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,iyear_cur,(GKCA(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,iyear_cur,(GKCM(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,iyear_cur,(GKCN(L,NY,NX),L=1,NLI(NY,NX))
!     WRITE(21,91)I,iyear_cur,(GKCK(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_CEC,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_AEC,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(TCS(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(TKS(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VHCP(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VHCM(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_gasml(idg_CO2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_solml(idg_CO2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_soHml(idg_CO2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_gasml(idg_CH4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_solml(idg_CH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_soHml(idg_CH4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(ROXYF(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RCO2F(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(ROXYL(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RCH4F(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RCH4L(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_gasml(idg_H2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_solml(idg_H2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_gasml(idg_O2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_solml(idg_O2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trc_soHml(idg_O2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(ROXYX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RNH4X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RNO3X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RNO2X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RN2OX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RPO4X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RP14X(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RNHBX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RN3BX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RN2BX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RPOBX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RP1BX(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,95)I,iyear_cur,((ROQCX(K,L,NY,NX),L=0,NLI(NY,NX)),K=1,jcplx)
      WRITE(21,95)I,iyear_cur,((ROQAX(K,L,NY,NX),L=0,NLI(NY,NX)),K=1,jcplx)
      WRITE(21,91)I,iyear_cur,(VOLWH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLIH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(VOLAH(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(RTDNT(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(ZL(L,NY,NX),L=0,JC)
      WRITE(21,91)I,iyear_cur,(ARLFT(L,NY,NX),L=1,JC)
      WRITE(21,91)I,iyear_cur,(ARSTT(L,NY,NX),L=1,JC)
      WRITE(21,91)I,iyear_cur,(WGLFT(L,NY,NX),L=1,JC)
      WRITE(21,91)I,iyear_cur,(trcs_VLN(ids_NH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcs_VLN(ids_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcs_VLN(ids_NO3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcs_VLN(ids_NO3B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcs_VLN(ids_H1PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcs_VLN(ids_H1PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_NH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_OHe,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_OH,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_OHp,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_HPO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_H2PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_OHeB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_OHB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_OHpB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_HPO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcx_solml(idx_H2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_AlPO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_FePO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_CaHPO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_HA,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_CaH2PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_AlPO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_FePO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_CaHPO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_HAB,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_CaH2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(CION(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Al,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Fe,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Hp,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Ca,L,NY,NX),L=1,NLI(NY,NX))
      IF(salt_model)THEN
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Al,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Fe,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Hp,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Ca,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Mg,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Na,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_K,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_OH,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_SO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_Cl,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_HCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_AlOH,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_AlOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_AlOH3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_AlOH4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_AlSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeOH,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeOH3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeOH4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaHCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_MgOH2,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_MgCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_MgHCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_MgSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_NaCO3,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_NaSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_KSO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_H0PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_H3PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeHPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_FeH2PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaHPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_CaH2PO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcs_solsml(idsa_MgHPO4,L,NY,NX),L=1,JS)
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Mg,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Na,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_OH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_SO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_Cl,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_HCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_AlOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_AlOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_AlOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_AlSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_NaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_NaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_KSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_H0PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_H3PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_H0PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_H3PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_FeH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_CaH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_solml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Al,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Fe,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Hp,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Ca,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Mg,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Na,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_OH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_SO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_Cl,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_HCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_AlOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_AlOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_AlOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_AlSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeOH,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeOH4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_MgOH2,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_MgCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_MgHCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_MgSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_NaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_NaSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_KSO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_H0PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_H3PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaH2PO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_MgHPO4,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_H0PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_H3PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_FeH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_CaH2PO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcsa_soHml(idsa_MgHPO4B,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_Hp,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_Al,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_Fe,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_Ca,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_Mg,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_Na,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcx_solml(idx_K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_AlOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_FeOH3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_CaCO3,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(21,91)I,iyear_cur,(trcp_salml(idsp_CaSO4,L,NY,NX),L=1,NLI(NY,NX))
      ENDIF
      D9985: DO K=1,jcplx
        DO  N=1,NFGs
          DO NGL=JGnio(N),JGnfo(N)
            WRITE(22,91)I,iyear_cur,(ROXYS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMX4(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMX3(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMX2(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMX1(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMB4(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMB3(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RVMB2(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RINHO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RINOO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RIPOO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RINHB(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RINOB(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(RIPBO(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(ROQCS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(ROQAS(NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,RINHOR(NGL,K,NY,NX)
            WRITE(22,91)I,iyear_cur,RINOOR(NGL,K,NY,NX)
            WRITE(22,91)I,iyear_cur,RIPOOR(NGL,K,NY,NX)
            DO M=1,nlbiomcp
              WRITE(22,91)I,iyear_cur,(OMC(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
              WRITE(22,91)I,iyear_cur,(OMN(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
              WRITE(22,91)I,iyear_cur,(OMP(M,NGL,K,L,NY,NX),L=0,NLI(NY,NX))
            enddo
          enddo
        enddo
      ENDDO D9985

      DO  N=1,NFGs
        DO NGL=JGniA(N),JGnfA(N)
          WRITE(22,91)I,iyear_cur,(ROXYSff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMX4ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMX3ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMX2ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMX1ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMB4ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMB3ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RVMB2ff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RINHOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RINOOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RIPOOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RINHBff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RINOBff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(RIPBOff(NGL,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,RINHORff(NGL,NY,NX)
          WRITE(22,91)I,iyear_cur,RINOORff(NGL,NY,NX)
          WRITE(22,91)I,iyear_cur,RIPOORff(NGL,NY,NX)
          DO M=1,nlbiomcp
            WRITE(22,91)I,iyear_cur,(OMCff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(OMNff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
            WRITE(22,91)I,iyear_cur,(OMPff(M,NGL,L,NY,NX),L=0,NLI(NY,NX))
          enddo
        enddo
      enddo


      D9980: DO K=1,jcplx
        D9975: DO M=1,ndbiomcp
          WRITE(22,91)I,iyear_cur,(ORC(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(ORN(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(ORP(M,K,L,NY,NX),L=0,NLI(NY,NX))
        ENDDO D9975
        WRITE(22,91)I,iyear_cur,(OQC(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQN(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQP(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQA(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQCH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQNH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQPH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OQAH(K,L,NY,NX),L=1,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OHC(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OHN(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OHP(K,L,NY,NX),L=0,NLI(NY,NX))
        WRITE(22,91)I,iyear_cur,(OHA(K,L,NY,NX),L=0,NLI(NY,NX))
        D9970: DO M=1,jsken
          WRITE(22,91)I,iyear_cur,(OSC(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(OSA(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(OSN(M,K,L,NY,NX),L=0,NLI(NY,NX))
          WRITE(22,91)I,iyear_cur,(OSP(M,K,L,NY,NX),L=0,NLI(NY,NX))
        ENDDO D9970
      ENDDO D9980
      WRITE(22,91)I,iyear_cur,(RVMXC(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(ORGC(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(ORGR(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_soil(ifert_nh4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_soil(ifert_nh3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_soil(ifert_urea,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_soil(ifert_no3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_band(ifert_nh4_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_band(ifert_nh3_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_band(ifert_urea_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(FertN_band(ifert_no3_band,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_gasml(idg_N2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(idg_N2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(idg_N2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_gasml(idg_N2O,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(idg_N2O,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(idg_N2O,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_gasml(idg_NH3,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_NH4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_NH4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(idg_NH3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(idg_NH3,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_NO3,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_NO3,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_NO2,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_NO2,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_H1PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_H2PO4,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_H1PO4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_H2PO4,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_NH4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_NH4B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(idg_NH3B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(idg_NH3B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_NO3B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_NO3B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_NO2B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_NO2B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_H1PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_solml(ids_H2PO4B,L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_H1PO4B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(trc_soHml(ids_H2PO4B,L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(WDNHB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(DPNHB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(WDNOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(DPNOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(WDPOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(DPPOB(L,NY,NX),L=1,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(ZNHUI(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(ZNHU0(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(ZNFNI(L,NY,NX),L=0,NLI(NY,NX))
      WRITE(22,91)I,iyear_cur,(ZNFN0(L,NY,NX),L=0,NLI(NY,NX))
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
