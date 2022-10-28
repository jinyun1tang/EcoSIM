
  SUBROUTINE routp(NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE READS PLANT CHECKPOINT FILES TO
!     RE-INITIALIZE THE MODEL FROM A SELECTED DATE IN AN EARLIER RUN
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StartqMod    , only : startq
  use GridConsts
  use FlagDataType
  use PlantDataRateType
  use EcoSIMCtrlDataType
  use PlantTraitDataType
  use PlantMngmtDataType
  use CanopyDataType
  use RootDataType
  use EcoSIMHistMod
  use CanopyRadDataType
  use GridDataType
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS



  integer :: NPP(JY,JX),IYR02(JP,JY,JX),IDAY02(JP,JY,JX) &
  ,IYRH2(JP,JY,JX),IDAYH2(JP,JY,JX)
  integer :: IDATE,IYR,K,L,M,NX,NY,NQ,N,NB,NR,NZ
!     execution begins here

  REWIND(26)
  REWIND(27)
  REWIND(28)
  REWIND(29)
!
!     READ CHECKPOINT FILES UP TO DATE OF RE-INITIALIZATION
!
8000  CONTINUE
  DO 9990 NX=NHW,NHE
    DO 9985 NY=NVN,NVS
      READ(26,90,END=1001)IDATE,IYR,NPP(NY,NX),IFLGT(NY,NX) &
        ,(IFLGC(NZ,NY,NX),NZ=1,NPP(NY,NX)) &
        ,(IDTH(NZ,NY,NX),NZ=1,NPP(NY,NX))
      READ(26,97)IDATE,IYR,NPP(NY,NX) &
        ,((WTSTDE(M,ielmc,NZ,NY,NX),M=1,jsken),NZ=1,NPP(NY,NX)) &
        ,((WTSTDE(M,ielmn,NZ,NY,NX),M=1,jsken),NZ=1,NPP(NY,NX)) &
        ,((WTSTDE(M,ielmp,NZ,NY,NX),M=1,jsken),NZ=1,NPP(NY,NX))
      IF(IFLGT(NY,NX).GT.0)THEN
        DO 9980 NQ=1,NPP(NY,NX)
          IF(IFLGC(NQ,NY,NX).EQ.1)THEN
            READ(26,91)IDATE,IYR,NZ,IYR02(NZ,NY,NX),IDAY02(NZ,NY,NX) &
              ,IYRH2(NZ,NY,NX),IDAYH2(NZ,NY,NX),NG(NZ,NY,NX),IDTHP(NZ,NY,NX) &
              ,IDTHR(NZ,NY,NX),NBR(NZ,NY,NX),NBT(NZ,NY,NX),NB1(NZ,NY,NX) &
              ,IFLGI(NZ,NY,NX),NRT(NZ,NY,NX),NIX(NZ,NY,NX),MY(NZ,NY,NX) &
              ,(NINR(NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
            READ(26,93)IDATE,IYR,NZ,TEUPTK(ielmc,NZ,NY,NX),TESNC(ielmc,NZ,NY,NX) &
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
            READ(26,92)IDATE,IYR,NZ,(GROUP(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(PSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(PSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(PSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(VSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(VSTGX(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(GSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(GSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(TGSTGI(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(TGSTGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(VRNS(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(VRNF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(VRNY(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(VRNZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(ATRP(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(FLG4(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,92)IDATE,IYR,NZ,(FLGZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(IFLGA(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(IFLGE(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(IFLGF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(IFLGR(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(IFLGQ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(IDTHB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(NBTB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(KLEAF(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(KVSTG(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(26,95)IDATE,IYR,NZ,(KVSTGN(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            DO  N=1,10
              READ(26,95)IDATE,IYR,NZ,(IDAY(N,NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            ENDDO
            READ(27,96)IDATE,IYR,NZ,EPOOLP(ielmc,NZ,NY,NX),EPOLNP(ielmc,NZ,NY,NX) &
              ,VHCPC(NZ,NY,NX),DTKC(NZ,NY,NX),WTSHTE(ielmc,NZ,NY,NX) &
              ,WTLS(NZ,NY,NX),WTRTE(ielmc,NZ,NY,NX),WTSHTE(ielmn,NZ,NY,NX),WTSHTE(ielmp,NZ,NY,NX) &
              ,WTLFE(ielmc,NZ,NY,NX),WTSHEE(ielmc,NZ,NY,NX),WTSTKE(ielmc,NZ,NY,NX),WTRSVE(ielmc,NZ,NY,NX) &
              ,WTHSKE(ielmc,NZ,NY,NX),WTEARE(ielmc,NZ,NY,NX),WTGRE(ielmc,NZ,NY,NX),WTNDE(ielmc,NZ,NY,NX) &
              ,WTRVE(ielmc,NZ,NY,NX),WTRVE(ielmn,NZ,NY,NX),WTRVE(ielmp,NZ,NY,NX),HTCTL(NZ,NY,NX) &
              ,SDPTH(NZ,NY,NX),WSTR(NZ,NY,NX) &
              ,CHILL(NZ,NY,NX),WTRTSE(ielmc,NZ,NY,NX),FRADP(NZ,NY,NX)
            READ(27,92)IDATE,IYR,NZ,(EPOOL(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(EPOOL(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(EPOOL(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(EPOLNB(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(EPOLNB(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(EPOLNB(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSHTBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTLFBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTNDBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSHEBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSTKBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WVSTKB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTRSVBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTHSKBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTEARBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTGRBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTLSB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSHTBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTLFBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTNDBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSHEBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSTKBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTRSVBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTHSKBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTEARBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTGRBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSHTBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTLFBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTNDBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSHEBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSTKBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTRSVBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTHSKBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTEARBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTGRBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(GRNXB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(GRNOB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(GRWTB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(ARLFB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WGLFEX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WGLFEX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WGLFEX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(ARLFZ(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(RCELX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(RCELX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(RCELX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WGSHEXE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WGSHEXE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WGSHEXE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(HTSHEX(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(RCESX(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(RCESX(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(RCESX(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSTXBE(NB,ielmc,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSTXBE(NB,ielmn,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,92)IDATE,IYR,NZ,(WTSTXBE(NB,ielmp,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
            READ(27,94)IDATE,IYR,NZ,(ARLFV(L,NZ,NY,NX),L=1,JC)
            READ(27,94)IDATE,IYR,NZ,(ARSTV(L,NZ,NY,NX),L=1,JC)
            D9945: DO NB=1,NBR(NQ,NY,NX)
              READ(28,93)IDATE,IYR,NZ,(CPOOL4(K,NB,NZ,NY,NX),K=1,JNODS)
              READ(28,93)IDATE,IYR,NZ,(CPOOL3(K,NB,NZ,NY,NX),K=1,JNODS)
              READ(28,93)IDATE,IYR,NZ,(CO2B(K,NB,NZ,NY,NX),K=1,JNODS)
              READ(28,93)IDATE,IYR,NZ,(HCOB(K,NB,NZ,NY,NX),K=1,JNODS)
              READ(28,93)IDATE,IYR,NZ,(ARLF(K,NB,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGLFE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WSLF(K,NB,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(HTSHE(K,NB,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGSHE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WSSHE(K,NB,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(HTNODE(K,NB,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(HTNODX(K,NB,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGNODE(K,NB,ielmc,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGLFE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGSHE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGNODE(K,NB,ielmn,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGLFE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGSHE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              READ(28,93)IDATE,IYR,NZ,(WGNODE(K,NB,ielmp,NZ,NY,NX),K=0,JNODS)
              D9950: DO K=0,JNODS
                READ(28,94)IDATE,IYR,NZ,(ARLFL(L,K,NB,NZ,NY,NX),L=1,JC)
                READ(28,94)IDATE,IYR,NZ,(WGLFLE(L,K,NB,ielmc,NZ,NY,NX),L=1,JC)
                READ(28,94)IDATE,IYR,NZ,(WGLFLE(L,K,NB,ielmn,NZ,NY,NX),L=1,JC)
                READ(28,94)IDATE,IYR,NZ,(WGLFLE(L,K,NB,ielmp,NZ,NY,NX),L=1,JC)
                IF(K.NE.0)THEN
                  D9940: DO N=1,4
                    READ(28,94)IDATE,IYR,NZ,(SURF(N,L,K,NB,NZ,NY,NX),L=1,JC)
                  ENDDO D9940
                ENDIF
              ENDDO D9950
              READ(28,94)IDATE,IYR,NZ,(ARSTK(L,NB,NZ,NY,NX),L=1,JC)
              DO L=1,JC
                READ(28,92)IDATE,IYR,NZ,(SURFB(N,L,NB,NZ,NY,NX),N=1,JLI)
              ENDDO
            ENDDO D9945
            DO 9970 N=1,MY(NQ,NY,NX)
              READ(29,94)IDATE,IYR,NZ,(PSIRT(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(PSIRO(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(PSIRG(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTN1(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTNL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTLGP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTDNP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTVLP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTVLW(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RRAD1(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RRAD2(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTARP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTLGA(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RCO2M(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RCO2A(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RCO2N(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(CO2A(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(OXYA(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(CH4A(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(Z2OA(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(ZH3A(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(CO2P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(OXYP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(CH4P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(Z2OP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(ZH3P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(WTRTL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(WTRTD(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(WSRTL(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(ROXYP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUNNHP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUNNOP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUPP2P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUPP1P(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUNNBP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUNNXP(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUPP2B(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RUPP1B(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(WFR(N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(EPOOLR(ielmc,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(EPOOLR(ielmn,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(EPOOLR(ielmp,N,L,NZ,NY,NX),L=1,NJ(NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTDP1(N,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTWT1E(N,NR,ielmc,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTWT1E(N,NR,ielmn,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              READ(29,94)IDATE,IYR,NZ,(RTWT1E(N,NR,ielmp,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
              D9965: DO NR=1,NRT(NZ,NY,NX)
                READ(29,94)IDATE,IYR,NZ,(RTN2(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(RTLG1(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(WTRT1E(ielmc,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(WTRT1E(ielmn,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(WTRT1E(ielmp,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(WTRT2E(ielmc,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(WTRT2E(ielmn,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(WTRT2E(ielmp,N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
                READ(29,94)IDATE,IYR,NZ,(RTLG2(N,L,NR,NZ,NY,NX),L=1,NJ(NY,NX))
              ENDDO D9965
9970        CONTINUE
            READ(29,94)IDATE,IYR,NZ,(EPOOLN(L,ielmc,NZ,NY,NX),L=1,NJ(NY,NX))
            READ(29,94)IDATE,IYR,NZ,(EPOOLN(L,ielmn,NZ,NY,NX),L=1,NJ(NY,NX))
            READ(29,94)IDATE,IYR,NZ,(EPOOLN(L,ielmp,NZ,NY,NX),L=1,NJ(NY,NX))
            READ(29,94)IDATE,IYR,NZ,(WTNDLE(L,ielmc,NZ,NY,NX),L=1,NJ(NY,NX))
            READ(29,94)IDATE,IYR,NZ,(WTNDLE(L,ielmn,NZ,NY,NX),L=1,NJ(NY,NX))
            READ(29,94)IDATE,IYR,NZ,(WTNDLE(L,ielmp,NZ,NY,NX),L=1,NJ(NY,NX))
          ENDIF
9980    CONTINUE
      ENDIF
9985  CONTINUE
9990  CONTINUE
  IF(IDATE.LT.IDAYR.OR.IYR.LT.IYRR)THEN
    GO TO 8000
  ELSEIF(IDATE.GE.IDAYR.AND.IYR.EQ.IYRR)THEN
    DO 9890 NX=NHW,NHE
      DO 9885 NY=NVN,NVS
        IF(IYR0(NZ,NY,NX).EQ.-1E+06)THEN
          IDAY0(NZ,NY,NX)=IDAY02(NZ,NY,NX)
          IYR0(NZ,NY,NX)=IYR02(NZ,NY,NX)
        ENDIF
        IF(IYRH(NZ,NY,NX).EQ.1E+06)THEN
          IDAYH(NZ,NY,NX)=IDAYH2(NZ,NY,NX)
          IYRH(NZ,NY,NX)=IYRH2(NZ,NY,NX)
        ENDIF
        IF(IFLGT(NY,NX).GT.0)THEN
          D9880: DO NZ=1,NP(NY,NX)
            IF(IFLGC(NZ,NY,NX).EQ.1)THEN
              D9600: DO NB=NBR(NZ,NY,NX)+1,10
                IFLGA(NB,NZ,NY,NX)=0
                IFLGE(NB,NZ,NY,NX)=1
                IFLGF(NB,NZ,NY,NX)=0
                IFLGR(NB,NZ,NY,NX)=0
                IFLGQ(NB,NZ,NY,NX)=0
                GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
                PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
                PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
                PSTGF(NB,NZ,NY,NX)=0.0
                VSTG(NB,NZ,NY,NX)=0.0
                VSTGX(NB,NZ,NY,NX)=0.0
                KLEAF(NB,NZ,NY,NX)=1
                KLEAFX(NB,NZ,NY,NX)=1
                KVSTG(NB,NZ,NY,NX)=1
                KVSTGN(NB,NZ,NY,NX)=0
                GSTGI(NB,NZ,NY,NX)=0.0
                GSTGF(NB,NZ,NY,NX)=0.0
                TGSTGI(NB,NZ,NY,NX)=0.0
                TGSTGF(NB,NZ,NY,NX)=0.0
                VRNY(NB,NZ,NY,NX)=VRNL(NB,NZ,NY,NX)+0.5
                VRNZ(NB,NZ,NY,NX)=0.0
                VRNS(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)
                VRNF(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)
                ATRP(NB,NZ,NY,NX)=0.0
                FDBK(NB,NZ,NY,NX)=1.0
                FDBKX(NB,NZ,NY,NX)=1.0
                FLG4(NB,NZ,NY,NX)=0
                NBTB(NB,NZ,NY,NX)=0
                IDTHB(NB,NZ,NY,NX)=1
                D15: DO M=1,10
                  IDAY(M,NB,NZ,NY,NX)=0
                ENDDO D15
                EPOOL(NB,1:npelms,NZ,NY,NX)=0.0
                WTSHTBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTSHEBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTSTKBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTLFBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTRSVBE(NB,1:npelms,NZ,NY,NX)=0.0
                WVSTKB(NB,NZ,NY,NX)=0.0
                WTHSKBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTEARBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTGRBE(NB,1:npelms,NZ,NY,NX)=0.0
                WTLSB(NB,NZ,NY,NX)=0.0
                GRNXB(NB,NZ,NY,NX)=0.0
                GRNOB(NB,NZ,NY,NX)=0.0
                GRWTB(NB,NZ,NY,NX)=0.0
                ARLFB(NB,NZ,NY,NX)=0.0
                RNH3B(NB,NZ,NY,NX)=0.0
                RCELX(NB,1:npelms,NZ,NY,NX)=0.0
                WGLFEX(NB,1:npelms,NZ,NY,NX)=0.0
                ARLFZ(NB,NZ,NY,NX)=0.0
                RCESX(NB,1:npelms,NZ,NY,NX)=0.0
                WTSTXBE(NB,1:npelms,NZ,NY,NX)=0.0
                WGSHEXE(NB,1:npelms,NZ,NY,NX)=0.0
                HTSHEX(NB,NZ,NY,NX)=0.0
                D5: DO L=1,JZ
                  ARSTK(L,NB,NZ,NY,NX)=0.0
                  DO N=1,JLI
                    SURFB(N,L,NB,NZ,NY,NX)=0.0
                  ENDDO
                ENDDO D5
              ENDDO D9600
              D9610: DO NB=NBR(NZ,NY,NX)+1,JBR
                D9620: DO K=0,JNODS
                  ARLF(K,NB,NZ,NY,NX)=0.0
                  HTNODE(K,NB,NZ,NY,NX)=0.0
                  HTNODX(K,NB,NZ,NY,NX)=0.0
                  HTSHE(K,NB,NZ,NY,NX)=0.0
                  WGLFE(K,NB,1:npelms,NZ,NY,NX)=0.0
                  WSLF(K,NB,NZ,NY,NX)=0.0
                  WGSHE(K,NB,1:npelms,NZ,NY,NX)=0.0
                  WSSHE(K,NB,NZ,NY,NX)=0.0
                  WGNODE(K,NB,1:npelms,NZ,NY,NX)=0.0
                  D55: DO L=1,JZ
                    ARLFL(L,K,NB,NZ,NY,NX)=0.0
                    WGLFLE(L,K,NB,1:npelms,NZ,NY,NX)=0.0
                  ENDDO D55
                  IF(K.NE.0)THEN
                    CPOOL3(K,NB,NZ,NY,NX)=0.0
                    CO2B(K,NB,NZ,NY,NX)=0.0
                    HCOB(K,NB,NZ,NY,NX)=0.0
                    CPOOL4(K,NB,NZ,NY,NX)=0.0
                    D45: DO L=1,JC
                      DO N=1,4
                        SURF(N,L,K,NB,NZ,NY,NX)=0.0
                      ENDDO
                    ENDDO D45
                  ENDIF
                ENDDO D9620
              ENDDO D9610
              D9700: DO NR=NRT(NZ,NY,NX)+1,10
                D9705: DO N=1,MY(NZ,NY,NX)
                  RTDP1(N,NR,NZ,NY,NX)=AMAX1(0.0,SDPTH(NZ,NY,NX) &
                    +CDPTH(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX))
                  D9710: DO L=1,JZ
                    WTRT1E(1:npelms,N,L,NR,NZ,NY,NX)=0.0
                    RTLG1(N,L,NR,NZ,NY,NX)=0.0
                    WTRT2E(1:npelms,N,L,NR,NZ,NY,NX)=0.0
                    RTLG2(N,L,NR,NZ,NY,NX)=0.0
                  ENDDO D9710
                ENDDO D9705
              ENDDO D9700
            ENDIF
          ENDDO D9880
        ENDIF
9885  CONTINUE
9890  CONTINUE
    !for pft 1 to 5, for each grid, only 5 pfts are allowed at maximum
    CALL STARTQ(NHW,NHE,NVN,NVS,1,5)
    GO TO 1000
  ELSE
    CALL STARTQ(NHW,NHE,NVN,NVS,1,5)
    GO TO 1000
  ENDIF
90    FORMAT(30I4)
91    FORMAT(3I4,4I8,30I4)
92    FORMAT(3I4,10E17.8E3)
93    FORMAT(3I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
94    FORMAT(3I4,21E17.8E3)
95    FORMAT(3I4,10I6)
96    FORMAT(3I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
97    FORMAT(3I4,15E17.8E3,/,15E17.8E3,/,15E17.8E3,/,15E17.8E3)
1001  CONTINUE
1000  RETURN
  END SUBROUTINE routp
