module Hour1Mod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod  , only : test_aeqb
  use abortutils   , only : endrun, print_info
  use EcosimConst
  use MicrobialDataType
  use SOMDataType
  use SoilChemDataType
  use FertilizerDataType
  use VegDataType
  implicit none

  private
  include "parameters.h"
  include "files.h"
  include "filec.h"
  include "blkc.h"
  include "blk1cp.h"
  include "blk1g.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk3.h"
  include "blk5.h"
  include "blk6.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk9a.h"
  include "blk9b.h"
  include "blk9c.h"
  include "blk10.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk15a.h"
  include "blk15b.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"
  include "blk19a.h"
  include "blk19b.h"
  include "blk19c.h"
  include "blk19d.h"
  include "blk20a.h"
  include "blk20b.h"
  include "blk20c.h"
  include "blk20d.h"
  include "blk20e.h"
  include "blk20f.h"
  include "blk21a.h"
  include "blk21b.h"
  include "blk22b.h"
  include "blk22c.h"

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__
  real(r8) :: OFC(2),OFN(2),OFP(2),CNOF(4),CPOF(4)
  real(r8) :: PSISK(0:100),THETK(100),TAUY(0:JC+1),RABSL(0:JC+1) &
    ,RABPL(0:JC+1),RAFSL(0:JC+1),RAFPL(0:JC+1),RADSL(JP,JY,JX) &
    ,RADPL(JP,JY,JX),RADS1(JP,JY,JX),RADP1(JP,JY,JX),RADS2(JP,JY,JX) &
    ,RADP2(JP,JY,JX),RAYSL(JP,JY,JX),RAYPL(JP,JY,JX),RAYS1(JP,JY,JX) &
    ,RAYP1(JP,JY,JX),RAYS2(JP,JY,JX),RAYP2(JP,JY,JX),RADSW(JP,JY,JX) &
    ,RADPW(JP,JY,JX),RADW1(JP,JY,JX),RADQ1(JP,JY,JX),RADW2(JP,JY,JX) &
    ,RADQ2(JP,JY,JX),RAYSW(JP,JY,JX),RAYPW(JP,JY,JX),RAYW1(JP,JY,JX) &
    ,RAYQ1(JP,JY,JX),RAYW2(JP,JY,JX),RAYQ2(JP,JY,JX),RADSA(JP,JY,JX) &
    ,RAPSA(JP,JY,JX),RDNDIR(4,4,JP,JY,JX),PARDIR(4,4,JP,JY,JX) &
    ,RADWA(JP,JY,JX),RAPWA(JP,JY,JX),RDNDIW(4,4,JP,JY,JX) &
    ,PARDIW(4,4,JP,JY,JX),TSURF(4,JZ,JP,JY,JX),TSURFB(4,JZ,JP,JY,JX) &
    ,XVOLWC(0:3),ZL1(0:JZ,JY,JX),THETPZ(JZ,JY,JX),TRADC(JY,JX) &
    ,TRAPC(JY,JX),TRADG(JY,JX),TRAPG(JY,JX),THETRX(0:2),DPTH0(JY,JX)
  integer :: IALBS(4,4)
!
!     *SG=diffusivity (m2 h-1):CG=CO2g,CL=CO2s,CH=CH4g,CQ=CH4s,OG=O2g
!     OL=O2s,ZG=N2g,ZL=N2s,Z2=N2Og,ZV=N2Os,ZH=NH3g,ZN=NH3s,ZO=NO3
!     PO=H2PO4,OC=DOC,ON=DON,OP=DOP,OA=acetate,WG=H2Og,AL=Al,FE=Fe
!     HY=H,CA=Ca,GM=Mg,AN=Na,AK=K,OH=OH,C3=CO3,HC=HCO3,SO=SO4,CL=Cl
!     HG=H2g,HL=H2s
!
  real(r8), PARAMETER :: CGSG=4.68E-02,CLSG=4.25E-06 &
    ,CHSG=7.80E-02,CQSG=7.08E-06,OGSG=6.43E-02,OLSG=8.57E-06 &
    ,ZGSG=5.57E-02,ZLSG=7.34E-06,Z2SG=5.57E-02,ZVSG=5.72E-06 &
    ,ZHSG=6.67E-02,ZNSG=4.00E-06,ZOSG=6.00E-06,POSG=3.00E-06 &
    ,OCSG=1.0E-08,ONSG=1.0E-08,OPSG=1.0E-08,OASG=3.64E-06 &
    ,WGSG=7.70E-02,ALSG=5.0E-06,FESG=5.0E-06,HYSG=5.0E-06 &
    ,CASG=5.0E-06,GMSG=5.0E-06,ANSG=5.0E-06,AKSG=5.0E-06 &
    ,OHSG=5.0E-06,C3SG=5.0E-06,HCSG=5.0E-06,SOSG=5.0E-06 &
    ,CLSX=5.0E-06,HGSG=5.57E-02,HLSG=7.34E-06
!
!     SC*X=solubility (g m-3/g m-3):CO2=CO2,CH4=CH4,OXY=O2,N2G=N2g
!     N2O=N2O,NH3=NH3,H2G=H2
!     AC*X=activity (g m-3):CO2=CO2,CH4=CH4,OXY=O2,N2G=N2g
!     N2O=N2O,NH3=NH3,H2G=H2
!
  real(r8), PARAMETER :: SCO2X=7.391E-01,SCH4X=3.156E-02 &
    ,SOXYX=2.925E-02,SN2GX=1.510E-02,SN2OX=5.241E-01 &
    ,SNH3X=2.852E+02,SH2GX=3.156E-02 &
    ,ACO2X=0.14,ACH4X=0.14,AOXYX=0.31,AN2GX=0.23,AN2OX=0.23 &
    ,ANH3X=0.07,AH2GX=0.14
!
!     ALBRW,ALBPW=stalk albedo for shortwave,PAR
!     VISCW=water viscosity (Mg m-1 s)
!     BKDSX=maximm soil bulk density, ZW=snowpack surface roughness (m)
!     CFW=stalk clumping factor,FORGW=minimum SOC or organic soil (g Mg-1)
!     THETPW=minimum air-filled porosity for saturation (m3 m-3)
!     RAM=minimum boundary layer resistance (h m-1)
!
  real(r8), PARAMETER :: ALBRW=0.1,ALBPW=0.1,ABSRW=1.0-ALBRW &
    ,ABSPW=1.0-ALBPW
  real(r8), PARAMETER :: VISCW=1.0E-06,BKDSX=1.89,ZW=0.01,CFW=0.5 &
    ,FORGW=0.25E+06,DTHETW=1.0E-06,THETPW=0.01,THETWP=1.0-THETPW &
    ,RAM=2.78E-03
!
!     XVOLWC=foliar water retention capacity (m3 m-2)
!     THETRX=litter water retention capacity (m3 g C-1)
!
  DATA XVOLWC/5.0E-04,2.5E-04,2.5E-04,2.5E-04/
  DATA THETRX/4.0E-06,8.0E-06,8.0E-06/

  REAL(r8) :: TFACL,TFACG,TFACW,TFACR,TFACA
  real(r8) :: ALBW,ALBG,ART,ARL,ARX,ARLSC,ARLSG,BKVLNX,BETAG
  real(r8) :: BETY,BETZ,BAREF,CORGCM,CORGM,COHS,CAC,CAS,CVRDF
  real(r8) :: CACX,CASX,CORGCX,CNOFT,CPOFT,D50,DC,DGAZI,DAZI
  real(r8) :: DZL,FVOLR,FH2O,FCX,FCLX,FCDX,FSNOW,FRADPT,FDPTHF
  real(r8) :: FDPTHM,FRNT,FRPT,H0PO4T,H1PO4T,H2PO4T,H3PO4T,OC
  real(r8) :: OSCI,OSNI,OSPI,OSCX,OSNX,OSPX,OMC1,OMN1,OMP1,OQC1
  real(r8) :: OQN1,OQP1,OSC1,OSN1,OSP1,PTDS,PSDX,PSIS1,PMA,PMB
  real(r8) :: PHA,PALPOT,PFEPOT,PCAPDT,PCAPHT,PCAPMT,PMAX,PMBX
  real(r8) :: PHAX,RADYL,RAPYL,RADYN,RADYW,RAPYN,RAPYW,RADST
  real(r8) :: RADWT,RADPT,RADQT,RA1ST,RA1WT,RA1PT,RA1QT,RA2ST
  real(r8) :: RA2WT,RA2PT,RA2QT,RADSG,RADYG,RAPSG,RAPYG,RASG
  real(r8) :: RAPG,RNT,RPT,SUM2,SUM1,SAGL,THETF,THETW1,THETWM
  real(r8) :: THETPX,TVOLG0,TVOLWI,THETWR,TSURFY,TSURFZ,TSURFS
  real(r8) :: TSURFX,TSURWY,TSURWZ,TSURWS,TSURWX,VISCWL,VORGC
  real(r8) :: VMINL,VSAND,VOLIT,VOLAT,VOLITL,VOLWTL,VOLATL,VOLWRZ
  real(r8) :: VOLIRZ,VOLWCX,WPX,WPLX,XJ,XK,XVOLW0,XVOLI0,XAREA
  real(r8) :: XTAUS,XTAUY,XN4T,XOH0T,XOH1T,XOH2T,XH1PT,XH2PT
  real(r8) :: YK,YAREA,ZD50,ZC3,ZA3,ZC2,ZA2,ZC1,ZA1,ZN,ZION1
  real(r8) :: ZAZI,ZAGL,ZX,ZY,ZE,ZZ,Z4A,Z3A,ZUA,ZOA,Z4B,Z3B
  real(r8) :: ZUB,ZOB,ZNH4T,ZNH3T,ZNO3T,ZNO2T,ZFE1PT,ZFE2PT
  real(r8) :: ZCA0PT,ZCA1PT,ZCA2PT,ZMG1PT,Z4AX,Z3AX,ZUAX,ZOAX
  real(r8) :: Z4BX,Z3BX,ZUBX,ZOBX
!     adding reals missed by the analyzer by hand
  real(r8) :: STOPS, STOPX, STOPY, STOPZ, STOPSZ, STOPYZ

  integer :: IFLGL,IFLGY,ICHKA,LFDPTH

  public :: hour1
  contains

  SUBROUTINE hour1(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE REINITIALIZES HOURLY VARIABLES USED IN OTHER
!     SUBROUTINES
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: L,NX,NY

  integer :: NZ,NR
!     execution begins here

  XJ=J
  DOY=I-1+XJ/24
!
!     RESET HOURLY SOIL ACCUMULATORS FOR WATER, HEAT, GASES, SOLUTES
!
  VOLWSO=0.0_r8
  HEATSO=0.0_r8
  OXYGSO=0.0_r8
  TLH2G=0.0_r8
  TSEDSO=0.0_r8
  TLRSDC=0.0_r8
  TLORGC=0.0_r8
  TLCO2G=0.0_r8
  TLRSDN=0.0_r8
  TLORGN=0.0_r8
  TLN2G=0.0_r8
  TLRSDP=0.0_r8
  TLORGP=0.0_r8
  TLNH4=0.0_r8
  TLNO3=0.0_r8
  TLPO4=0.0_r8
  TION=0.0_r8
  TBALC=0.0_r8
  TBALN=0.0_r8
  TBALP=0.0_r8
  call SetAtmosTracerConcentration(I,NHW,NHE,NVN,NVS)
!
!     RESET FLUX ARRAYS USED IN OTHER SUBROUTINES
!
  call ResetFluxArrays(I,NHW,NHE,NVN,NVS)
!
!     IF SALT FLAG SET
!
  call ResetSaltModelArrays(NHW,NHE,NVN,NVS)
!
!     RESET SOIL PROPERTIES AND PEDOTRANSFER FUNCTIONS
!     FOLLOWING ANY SOIL DISTURBANCE
!
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      IF(J.EQ.1)THEN
        IFLGT(NY,NX)=0
        DO 9905 NZ=1,NP(NY,NX)
          PSILZ(NZ,NY,NX)=0.0_r8
9905    CONTINUE
      ENDIF
!
!     HYDROLOGICAL PRPOERTIES OR SURFACE LITTER
!
!     VOLWRX=liter water holding capacity
!     VOLR=dry litter volume
!     POROS0,FC,WP=litter porosity,field capacity,wilting point
!
      VOLWRX(NY,NX)=AMAX1(0.0,THETRX(0)*RC0(0,NY,NX) &
      +THETRX(1)*RC0(1,NY,NX)+THETRX(2)*RC0(2,NY,NX))
      VOLR(NY,NX)=AMAX1(0.0,RC0(0,NY,NX)*1.0E-06/BKRS(0) &
      +RC0(1,NY,NX)*1.0E-06/BKRS(1)+RC0(2,NY,NX)*1.0E-06/BKRS(2))
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
        FVOLR=VOLWRX(NY,NX)/VOLR(NY,NX)
      ELSE
        FVOLR=THETRX(1)/BKRS(1)
      ENDIF
      POROS0(NY,NX)=FVOLR
      FC(0,NY,NX)=0.500*FVOLR
      WP(0,NY,NX)=0.125*FVOLR
      PSL(0,NY,NX)=LOG(POROS0(NY,NX))
      FCL(0,NY,NX)=LOG(FC(0,NY,NX))
      WPL(0,NY,NX)=LOG(WP(0,NY,NX))
      PSD(0,NY,NX)=PSL(0,NY,NX)-FCL(0,NY,NX)
      FCD(0,NY,NX)=FCL(0,NY,NX)-WPL(0,NY,NX)
      SRP(0,NY,NX)=1.00
!
!     RESET SURFACE LITTER PHYSICAL PROPERTIES (DENSITY, TEXTURE)
!     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
!
      call SetLiterSoilPropAftDisturbance(I,J,NY,NX)
!
!     END OF RESET AFTER DISTURBANCE
!
!     PARAMETERS FOR COHESION, EROSIVITY, AND ROUGHNESS OF SURFACE SOIL USED
!     FOR SURFACE WATER AND SEDIMENT TRANSPORT IN 'EROSION'

      call SetSurfaceProperty4SedErosion(NY,NX)
!
!     RESET HOURLY ACCUMULATORS
      call SetHourlyAccumulators(NY,NX)
!
!     RESET ARRAYS TO TRANSFER MATERIALS WITHIN SOILS
!     AND BETWEEN SOILS AND PLANTS
!
      DO 9875 L=0,NL(NY,NX)

        call SetArrays4PlantSoilTransfer(L,NY,NX)
!
!     IF SOC FLAG IS SET
!
        IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
          call UpdateTotalSOC(L,NY,NX)
        ENDIF

9875  CONTINUE
      IFLGL=0
      IFLGY=0
      ICHKA=0
      DO 9985 L=NUI(NY,NX),NLI(NY,NX)

        call ZeroHourlyArrays(L,NY,NX)
!
        call GetChemicalConcsInSoil(L,NY,NX)
!     WRITE(*,1113)'CORGC',I,J,NX,NY,L,CORGC(L,NY,NX)
!    2,ORGC(L,NY,NX),BKVL(L,NY,NX),BKDS(L,NY,NX),VOLX(L,NY,NX)
!
        call GetSoluteConcentrations(L,NY,NX)

        call PrepVars4PlantMicrobeUptake(L,NY,NX)
!
        call GetSoilHydraulicVars(L,NY,NX)
!
!     CALCULATE ACTIVE LAYER DEPTH

        call DiagActiveLayerDepth(L,NY,NX)
!
!     OUTPUT FOR WATER TABLE DEPTH
!
        call GetOutput4WaterTableDepth(L,NY,NX)
9985  CONTINUE
!

      call GetSurfResidualProperties(NY,NX)

      call SetTracerPropertyInLiterAir(NY,NX)
!
      call MultiLayerSurfaceRadiation(I,J,NY,NX)
!     IF(NX.EQ.4.AND.NY.EQ.5)THEN
!     DO 140 NZ=1,NP(NY,NX)
!     WRITE(*,1926)'CANOPY',IYRC,I,J,NX,NY,NZ,FRADP(NZ,NY,NX)
!    2,RADP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX),RAP(NY,NX)
!    2,FRADG(NY,NX),ARLFS(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
!    3,ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
!    4,ARSTP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
!    4,ARLSS(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
!    5,SSIN(NY,NX),DPTHS(NY,NX),FRADPT
!1926  FORMAT(A10,6I6,30E12.4)
!140   CONTINUE
!     ENDIF

      call DivideCanopyLayerByLAI(NY,NX)

      call CalcBoundaryLayerProperties(NY,NX)
!
!     RESET HOURLY INDICATORS
!
      THRMCX(NY,NX)=THRMC(NY,NX)
      THRMGX(NY,NX)=THRMG(NY,NX)
      CNETX(NY,NX)=TCNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      THRMC(NY,NX)=0.0_r8
      THRMG(NY,NX)=0.0_r8
      TLEX(NY,NX)=TLEC(NY,NX)
      TSHX(NY,NX)=TSHC(NY,NX)
      TLEC(NY,NX)=0.0_r8
      TSHC(NY,NX)=0.0_r8
      TRN(NY,NX)=0.0_r8
      TLE(NY,NX)=0.0_r8
      TSH(NY,NX)=0.0_r8
      TGH(NY,NX)=0.0_r8
      TCCAN(NY,NX)=0.0_r8
      TCNET(NY,NX)=0.0_r8
      RECO(NY,NX)=0.0_r8
!
!     CANOPY RETENTION OF PRECIPITATION
!
!     XVOLWC=foliar surface water retention capacity
!     ARLFP,ARSTP=leaf,stalk area of PFT
!     FLWC,TFLWC=water retention of PFT,combined canopy
!     PRECA=precipitation+irrigation
!     FRADP=fraction of radiation received by each PFT canopy
!     VOLWC=canopy surface water retention
!
      DO 1930 NZ=1,NP(NY,NX)
        VOLWCX=XVOLWC(IGTYP(NZ,NY,NX))*(ARLFP(NZ,NY,NX)+ARSTP(NZ,NY,NX))
        FLWC(NZ,NY,NX)=AMAX1(0.0,AMIN1(PRECA(NY,NX)*FRADP(NZ,NY,NX) &
          ,VOLWCX-VOLWC(NZ,NY,NX)))
        TFLWCI(NY,NX)=TFLWCI(NY,NX)+PRECA(NY,NX)*FRADP(NZ,NY,NX)
        TFLWC(NY,NX)=TFLWC(NY,NX)+FLWC(NZ,NY,NX)
!
!     NUMBERS OF TOP AND BOTTOM ROOTED SOIL LAYERS
!
!     NG=number of uppermost rooted layer
!     NINR=number of lowest rooted layer
!
        NG(NZ,NY,NX)=MAX(NG(NZ,NY,NX),NU(NY,NX))
        NIX(NZ,NY,NX)=MAX(NIX(NZ,NY,NX),NU(NY,NX))
        DO 9790 NR=1,10
          NINR(NR,NZ,NY,NX)=MAX(NINR(NR,NZ,NY,NX),NU(NY,NX))
9790    CONTINUE
1930  CONTINUE
!
!     WRITE SW AND PAR ALBEDO
!
!     IF(ABS(J-ZNOON(NY,NX)).LT.1)THEN
!     IF(RAD(NY,NX).GT.0.0.AND.RAP(NY,NX).GT.0.0)THEN
!     WRITE(19,1927)'ALBEDO',IYRC,I,J,NX,NY
!    2,RAD(NY,NX),RAP(NY,NX),TRADC(NY,NX),TRAPC(NY,NX)
!    3,TRADG(NY,NX),TRAPG(NY,NX)
!    4,(RAD(NY,NX)-TRADC(NY,NX)-TRADG(NY,NX))/RAD(NY,NX)
!    5,(RAP(NY,NX)-TRAPC(NY,NX)-TRAPG(NY,NX))/RAP(NY,NX)
!1927  FORMAT(A10,5I6,30E12.4)
!     ENDIF
!     ENDIF
9990  CONTINUE
9995  CONTINUE
!
!     FERTILIZER APPLICATIONS OCCUR AT SOLAR NOON
  call ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)
  RETURN

  END subroutine hour1
!------------------------------------------------------------------------------------------

  subroutine SetAtmosTracerConcentration(I,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: NX, NY
!     begin_execution
!
!     CONCENTRATIONS OF CO2, CH4, O2, N2, N2O, NH3, H2 IN ATMOSPHERE,
!     PRECIPITATION AND IRRIGATION FROM MIXING RATIOS READ IN 'READS'
!
!     C*E,C*R,C*Q=atmospheric,precipitation,irrigation solute concentrations
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
  DO 9145 NX=NHW,NHE
    DO 9140 NY=NVN,NVS
      CCO2E(NY,NX)=CO2E(NY,NX)*5.36E-04*TREF/TKA(NY,NX)
      CCH4E(NY,NX)=CH4E(NY,NX)*5.36E-04*TREF/TKA(NY,NX)
      COXYE(NY,NX)=OXYE(NY,NX)*1.43E-03*TREF/TKA(NY,NX)
      CZ2GE(NY,NX)=Z2GE(NY,NX)*1.25E-03*TREF/TKA(NY,NX)
      CZ2OE(NY,NX)=Z2OE(NY,NX)*1.25E-03*TREF/TKA(NY,NX)
      CNH3E(NY,NX)=ZNH3E(NY,NX)*6.25E-04*TREF/TKA(NY,NX)
      CH2GE(NY,NX)=H2GE(NY,NX)*8.92E-05*TREF/TKA(NY,NX)
      CCOR(NY,NX)=CCO2E(NY,NX)*SCO2X/(EXP(ACO2X*CSTRR(NY,NX))) &
        *EXP(0.843-0.0281*TCA(NY,NX))
      CCHR(NY,NX)=CCH4E(NY,NX)*SCH4X/(EXP(ACH4X*CSTRR(NY,NX))) &
        *EXP(0.597-0.0199*TCA(NY,NX))
      COXR(NY,NX)=COXYE(NY,NX)*SOXYX/(EXP(AOXYX*CSTRR(NY,NX))) &
        *EXP(0.516-0.0172*TCA(NY,NX))
      CNNR(NY,NX)=CZ2GE(NY,NX)*SN2GX/(EXP(AN2GX*CSTRR(NY,NX))) &
        *EXP(0.456-0.0152*TCA(NY,NX))
      CN2R(NY,NX)=CZ2OE(NY,NX)*SN2OX/(EXP(AN2OX*CSTRR(NY,NX))) &
        *EXP(0.897-0.0299*TCA(NY,NX))
      CCOQ(NY,NX)=CCO2E(NY,NX)*SCO2X/(EXP(ACO2X*CSTRQ(I,NY,NX))) &
        *EXP(0.843-0.0281*TCA(NY,NX))
      CCHQ(NY,NX)=CCH4E(NY,NX)*SCH4X/(EXP(ACH4X*CSTRQ(I,NY,NX))) &
        *EXP(0.597-0.0199*TCA(NY,NX))
      COXQ(NY,NX)=COXYE(NY,NX)*SOXYX/(EXP(AOXYX*CSTRQ(I,NY,NX))) &
        *EXP(0.516-0.0172*TCA(NY,NX))
      CNNQ(NY,NX)=CZ2GE(NY,NX)*SN2GX/(EXP(AN2GX*CSTRQ(I,NY,NX))) &
        *EXP(0.456-0.0152*TCA(NY,NX))
      CN2Q(NY,NX)=CZ2OE(NY,NX)*SN2OX/(EXP(AN2OX*CSTRQ(I,NY,NX))) &
        *EXP(0.897-0.0299*TCA(NY,NX))
9140  CONTINUE
9145  CONTINUE
  end subroutine SetAtmosTracerConcentration
!------------------------------------------------------------------------------------------

  subroutine ResetFluxArrays(I,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: L,N,NX,NY,K,NN,NO,M,NGL

!     begin_execution
  DO 9895 NX=NHW,NHE+1
    DO 9890 NY=NVN,NVS+1
      DO 9885 L=0,NL(NY,NX)+1
        DO 9880 N=1,3
!
!     WATER,SNOW,SOLUTE RUNOFF
!
        IF(L.EQ.0.AND.N.NE.3)THEN
          DO 9835 NN=1,2
            QR(N,NN,NY,NX)=0.0_r8
            HQR(N,NN,NY,NX)=0.0_r8
            DO 9870 K=0,4
              XOCQRS(K,N,NN,NY,NX)=0.0_r8
              XONQRS(K,N,NN,NY,NX)=0.0_r8
              XOPQRS(K,N,NN,NY,NX)=0.0_r8
              XOAQRS(K,N,NN,NY,NX)=0.0_r8
9870        CONTINUE
            XCOQRS(N,NN,NY,NX)=0.0_r8
            XCHQRS(N,NN,NY,NX)=0.0_r8
            XOXQRS(N,NN,NY,NX)=0.0_r8
            XNGQRS(N,NN,NY,NX)=0.0_r8
            XN2QRS(N,NN,NY,NX)=0.0_r8
            XHGQRS(N,NN,NY,NX)=0.0_r8
            XN4QRW(N,NN,NY,NX)=0.0_r8
            XN3QRW(N,NN,NY,NX)=0.0_r8
            XNOQRW(N,NN,NY,NX)=0.0_r8
            XNXQRS(N,NN,NY,NX)=0.0_r8
            XP1QRW(N,NN,NY,NX)=0.0_r8
            XP4QRW(N,NN,NY,NX)=0.0_r8
9835      CONTINUE
          QS(N,NY,NX)=0.0_r8
          QW(N,NY,NX)=0.0_r8
          QI(N,NY,NX)=0.0_r8
          HQS(N,NY,NX)=0.0_r8
          XCOQSS(N,NY,NX)=0.0_r8
          XCHQSS(N,NY,NX)=0.0_r8
          XOXQSS(N,NY,NX)=0.0_r8
          XNGQSS(N,NY,NX)=0.0_r8
          XN2QSS(N,NY,NX)=0.0_r8
          XN4QSS(N,NY,NX)=0.0_r8
          XN3QSS(N,NY,NX)=0.0_r8
          XNOQSS(N,NY,NX)=0.0_r8
          XP1QSS(N,NY,NX)=0.0_r8
          XP4QSS(N,NY,NX)=0.0_r8
!
!     IF EROSION FLAG SET
!
          IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
            DO 9855 NN=1,2
              XSEDER(N,NN,NY,NX)=0.0_r8
              XSANER(N,NN,NY,NX)=0.0_r8
              XSILER(N,NN,NY,NX)=0.0_r8
              XCLAER(N,NN,NY,NX)=0.0_r8
              XCECER(N,NN,NY,NX)=0.0_r8
              XAECER(N,NN,NY,NX)=0.0_r8
              XNH4ER(N,NN,NY,NX)=0.0_r8
              XNH3ER(N,NN,NY,NX)=0.0_r8
              XNHUER(N,NN,NY,NX)=0.0_r8
              XNO3ER(N,NN,NY,NX)=0.0_r8
              XNH4EB(N,NN,NY,NX)=0.0_r8
              XNH3EB(N,NN,NY,NX)=0.0_r8
              XNHUEB(N,NN,NY,NX)=0.0_r8
              XNO3EB(N,NN,NY,NX)=0.0_r8
              XN4ER(N,NN,NY,NX)=0.0_r8
              XNBER(N,NN,NY,NX)=0.0_r8
              XHYER(N,NN,NY,NX)=0.0_r8
              XALER(N,NN,NY,NX)=0.0_r8
              XCAER(N,NN,NY,NX)=0.0_r8
              XMGER(N,NN,NY,NX)=0.0_r8
              XNAER(N,NN,NY,NX)=0.0_r8
              XKAER(N,NN,NY,NX)=0.0_r8
              XHCER(N,NN,NY,NX)=0.0_r8
              XAL2ER(N,NN,NY,NX)=0.0_r8
              XOH0ER(N,NN,NY,NX)=0.0_r8
              XOH1ER(N,NN,NY,NX)=0.0_r8
              XOH2ER(N,NN,NY,NX)=0.0_r8
              XH1PER(N,NN,NY,NX)=0.0_r8
              XH2PER(N,NN,NY,NX)=0.0_r8
              XOH0EB(N,NN,NY,NX)=0.0_r8
              XOH1EB(N,NN,NY,NX)=0.0_r8
              XOH2EB(N,NN,NY,NX)=0.0_r8
              XH1PEB(N,NN,NY,NX)=0.0_r8
              XH2PEB(N,NN,NY,NX)=0.0_r8
              PALOER(N,NN,NY,NX)=0.0_r8
              PFEOER(N,NN,NY,NX)=0.0_r8
              PCACER(N,NN,NY,NX)=0.0_r8
              PCASER(N,NN,NY,NX)=0.0_r8
              PALPER(N,NN,NY,NX)=0.0_r8
              PFEPER(N,NN,NY,NX)=0.0_r8
              PCPDER(N,NN,NY,NX)=0.0_r8
              PCPHER(N,NN,NY,NX)=0.0_r8
              PCPMER(N,NN,NY,NX)=0.0_r8
              PALPEB(N,NN,NY,NX)=0.0_r8
              PFEPEB(N,NN,NY,NX)=0.0_r8
              PCPDEB(N,NN,NY,NX)=0.0_r8
              PCPHEB(N,NN,NY,NX)=0.0_r8
              PCPMEB(N,NN,NY,NX)=0.0_r8
              DO 9480 K=0,5
                DO  NO=1,7
                  DO NGL=1,JG
                  OMCER(3+(NGL-1)*3,NO,K,N,NN,NY,NX)=0.0_r8
                  DO  M=1,2
                    OMCER(M+(NGL-1)*3,NO,K,N,NN,NY,NX)=0.0_r8
                    OMNER(M+(NGL-1)*3,NO,K,N,NN,NY,NX)=0.0_r8
                    OMPER(M+(NGL-1)*3,NO,K,N,NN,NY,NX)=0.0_r8
                  enddo
                enddo
                ENDDO
9480          CONTINUE
              DO 9475 K=0,4
                DO 9470 M=1,2
                  ORCER(M,K,N,NN,NY,NX)=0.0_r8
                  ORNER(M,K,N,NN,NY,NX)=0.0_r8
                  ORPER(M,K,N,NN,NY,NX)=0.0_r8
9470            CONTINUE
                OHCER(K,N,NN,NY,NX)=0.0_r8
                OHNER(K,N,NN,NY,NX)=0.0_r8
                OHPER(K,N,NN,NY,NX)=0.0_r8
                DO 9465 M=1,4
                  OSCER(M,K,N,NN,NY,NX)=0.0_r8
                  OSAER(M,K,N,NN,NY,NX)=0.0_r8
                  OSNER(M,K,N,NN,NY,NX)=0.0_r8
                  OSPER(M,K,N,NN,NY,NX)=0.0_r8
9465            CONTINUE
9475          CONTINUE
9855        CONTINUE
          ENDIF
        ENDIF
!
!     GAS AND SOLUTE FLUXES
!
        XCOFLS(N,L,NY,NX)=0.0_r8
        XCHFLS(N,L,NY,NX)=0.0_r8
        XOXFLS(N,L,NY,NX)=0.0_r8
        XNGFLS(N,L,NY,NX)=0.0_r8
        XN2FLS(N,L,NY,NX)=0.0_r8
        XHGFLS(N,L,NY,NX)=0.0_r8
        XN4FLW(N,L,NY,NX)=0.0_r8
        XN3FLW(N,L,NY,NX)=0.0_r8
        XNOFLW(N,L,NY,NX)=0.0_r8
        XNXFLS(N,L,NY,NX)=0.0_r8
        XH1PFS(N,L,NY,NX)=0.0_r8
        XH2PFS(N,L,NY,NX)=0.0_r8
        DO 9860 K=0,4
          XOCFLS(K,N,L,NY,NX)=0.0_r8
          XONFLS(K,N,L,NY,NX)=0.0_r8
          XOPFLS(K,N,L,NY,NX)=0.0_r8
          XOAFLS(K,N,L,NY,NX)=0.0_r8
9860    CONTINUE
9880  CONTINUE
!
!     BAND AND MACROPORE FLUXES
!
      IF(L.NE.0)THEN
        DO 9840 N=1,3
          FLW(N,L,NY,NX)=0.0_r8
          FLWX(N,L,NY,NX)=0.0_r8
          FLWH(N,L,NY,NX)=0.0_r8
          HFLW(N,L,NY,NX)=0.0_r8
          XN4FLB(N,L,NY,NX)=0.0_r8
          XN3FLB(N,L,NY,NX)=0.0_r8
          XNOFLB(N,L,NY,NX)=0.0_r8
          XNXFLB(N,L,NY,NX)=0.0_r8
          XH1BFB(N,L,NY,NX)=0.0_r8
          XH2BFB(N,L,NY,NX)=0.0_r8
          XCOFHS(N,L,NY,NX)=0.0_r8
          XCHFHS(N,L,NY,NX)=0.0_r8
          XOXFHS(N,L,NY,NX)=0.0_r8
          XNGFHS(N,L,NY,NX)=0.0_r8
          XN2FHS(N,L,NY,NX)=0.0_r8
          XHGFHS(N,L,NY,NX)=0.0_r8
          XN4FHW(N,L,NY,NX)=0.0_r8
          XN3FHW(N,L,NY,NX)=0.0_r8
          XNOFHW(N,L,NY,NX)=0.0_r8
          XNXFHS(N,L,NY,NX)=0.0_r8
          XH1PHS(N,L,NY,NX)=0.0_r8
          XH2PHS(N,L,NY,NX)=0.0_r8
          XN4FHB(N,L,NY,NX)=0.0_r8
          XN3FHB(N,L,NY,NX)=0.0_r8
          XNOFHB(N,L,NY,NX)=0.0_r8
          XNXFHB(N,L,NY,NX)=0.0_r8
          XH1BHB(N,L,NY,NX)=0.0_r8
          XH2BHB(N,L,NY,NX)=0.0_r8
          XCOFLG(N,L,NY,NX)=0.0_r8
          XCHFLG(N,L,NY,NX)=0.0_r8
          XOXFLG(N,L,NY,NX)=0.0_r8
          XNGFLG(N,L,NY,NX)=0.0_r8
          XN2FLG(N,L,NY,NX)=0.0_r8
          XN3FLG(N,L,NY,NX)=0.0_r8
          XHGFLG(N,L,NY,NX)=0.0_r8
          DO 9820 K=0,4
            XOCFHS(K,N,L,NY,NX)=0.0_r8
            XONFHS(K,N,L,NY,NX)=0.0_r8
            XOPFHS(K,N,L,NY,NX)=0.0_r8
            XOAFHS(K,N,L,NY,NX)=0.0_r8
9820      CONTINUE
9840    CONTINUE
      ENDIF
9885  CONTINUE
9890  CONTINUE
9895  CONTINUE
  end subroutine ResetFluxArrays
!------------------------------------------------------------------------------------------

  subroutine ResetSaltModelArrays(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: N,NX,NY,L,NN
!     begin_execution

  DO 8895 NX=NHW,NHE+1
    DO 8890 NY=NVN,NVS+1
      IF(ISALTG.NE.0)THEN
        DO 8885 L=1,NL(NY,NX)+1
          DO 8880 N=1,3
            IF(L.EQ.1.AND.N.NE.3)THEN
              DO 9836 NN=1,2
                XQRAL(N,NN,NY,NX)=0.0_r8
                XQRFE(N,NN,NY,NX)=0.0_r8
                XQRHY(N,NN,NY,NX)=0.0_r8
                XQRCA(N,NN,NY,NX)=0.0_r8
                XQRMG(N,NN,NY,NX)=0.0_r8
                XQRNA(N,NN,NY,NX)=0.0_r8
                XQRKA(N,NN,NY,NX)=0.0_r8
                XQROH(N,NN,NY,NX)=0.0_r8
                XQRSO(N,NN,NY,NX)=0.0_r8
                XQRCL(N,NN,NY,NX)=0.0_r8
                XQRC3(N,NN,NY,NX)=0.0_r8
                XQRHC(N,NN,NY,NX)=0.0_r8
                XQRAL1(N,NN,NY,NX)=0.0_r8
                XQRAL2(N,NN,NY,NX)=0.0_r8
                XQRAL3(N,NN,NY,NX)=0.0_r8
                XQRAL4(N,NN,NY,NX)=0.0_r8
                XQRALS(N,NN,NY,NX)=0.0_r8
                XQRFE1(N,NN,NY,NX)=0.0_r8
                XQRFE2(N,NN,NY,NX)=0.0_r8
                XQRFE3(N,NN,NY,NX)=0.0_r8
                XQRFE4(N,NN,NY,NX)=0.0_r8
                XQRFES(N,NN,NY,NX)=0.0_r8
                XQRCAO(N,NN,NY,NX)=0.0_r8
                XQRCAC(N,NN,NY,NX)=0.0_r8
                XQRCAH(N,NN,NY,NX)=0.0_r8
                XQRCAS(N,NN,NY,NX)=0.0_r8
                XQRMGO(N,NN,NY,NX)=0.0_r8
                XQRMGC(N,NN,NY,NX)=0.0_r8
                XQRMGH(N,NN,NY,NX)=0.0_r8
                XQRMGS(N,NN,NY,NX)=0.0_r8
                XQRNAC(N,NN,NY,NX)=0.0_r8
                XQRNAS(N,NN,NY,NX)=0.0_r8
                XQRKAS(N,NN,NY,NX)=0.0_r8
                XQRH0P(N,NN,NY,NX)=0.0_r8
                XQRH3P(N,NN,NY,NX)=0.0_r8
                XQRF1P(N,NN,NY,NX)=0.0_r8
                XQRF2P(N,NN,NY,NX)=0.0_r8
                XQRC0P(N,NN,NY,NX)=0.0_r8
                XQRC1P(N,NN,NY,NX)=0.0_r8
                XQRC2P(N,NN,NY,NX)=0.0_r8
                XQRM1P(N,NN,NY,NX)=0.0_r8
9836          CONTINUE
              XQSAL(N,NY,NX)=0.0_r8
              XQSFE(N,NY,NX)=0.0_r8
              XQSHY(N,NY,NX)=0.0_r8
              XQSCA(N,NY,NX)=0.0_r8
              XQSMG(N,NY,NX)=0.0_r8
              XQSNA(N,NY,NX)=0.0_r8
              XQSKA(N,NY,NX)=0.0_r8
              XQSOH(N,NY,NX)=0.0_r8
              XQSSO(N,NY,NX)=0.0_r8
              XQSCL(N,NY,NX)=0.0_r8
              XQSC3(N,NY,NX)=0.0_r8
              XQSHC(N,NY,NX)=0.0_r8
              XQSAL1(N,NY,NX)=0.0_r8
              XQSAL2(N,NY,NX)=0.0_r8
              XQSAL3(N,NY,NX)=0.0_r8
              XQSAL4(N,NY,NX)=0.0_r8
              XQSALS(N,NY,NX)=0.0_r8
              XQSFE1(N,NY,NX)=0.0_r8
              XQSFE2(N,NY,NX)=0.0_r8
              XQSFE3(N,NY,NX)=0.0_r8
              XQSFE4(N,NY,NX)=0.0_r8
              XQSFES(N,NY,NX)=0.0_r8
              XQSCAO(N,NY,NX)=0.0_r8
              XQSCAC(N,NY,NX)=0.0_r8
              XQSCAH(N,NY,NX)=0.0_r8
              XQSCAS(N,NY,NX)=0.0_r8
              XQSMGO(N,NY,NX)=0.0_r8
              XQSMGC(N,NY,NX)=0.0_r8
              XQSMGH(N,NY,NX)=0.0_r8
              XQSMGS(N,NY,NX)=0.0_r8
              XQSNAC(N,NY,NX)=0.0_r8
              XQSNAS(N,NY,NX)=0.0_r8
              XQSKAS(N,NY,NX)=0.0_r8
              XQSH0P(N,NY,NX)=0.0_r8
              XQSH1P(N,NY,NX)=0.0_r8
              XQSH3P(N,NY,NX)=0.0_r8
              XQSF1P(N,NY,NX)=0.0_r8
              XQSF2P(N,NY,NX)=0.0_r8
              XQSC0P(N,NY,NX)=0.0_r8
              XQSC1P(N,NY,NX)=0.0_r8
              XQSC2P(N,NY,NX)=0.0_r8
              XQSM1P(N,NY,NX)=0.0_r8
            ENDIF
            XALFLS(N,L,NY,NX)=0.0_r8
            XFEFLS(N,L,NY,NX)=0.0_r8
            XHYFLS(N,L,NY,NX)=0.0_r8
            XCAFLS(N,L,NY,NX)=0.0_r8
            XMGFLS(N,L,NY,NX)=0.0_r8
            XNAFLS(N,L,NY,NX)=0.0_r8
            XKAFLS(N,L,NY,NX)=0.0_r8
            XOHFLS(N,L,NY,NX)=0.0_r8
            XSOFLS(N,L,NY,NX)=0.0_r8
            XCLFLS(N,L,NY,NX)=0.0_r8
            XC3FLS(N,L,NY,NX)=0.0_r8
            XHCFLS(N,L,NY,NX)=0.0_r8
            XAL1FS(N,L,NY,NX)=0.0_r8
            XAL2FS(N,L,NY,NX)=0.0_r8
            XAL3FS(N,L,NY,NX)=0.0_r8
            XAL4FS(N,L,NY,NX)=0.0_r8
            XALSFS(N,L,NY,NX)=0.0_r8
            XFE1FS(N,L,NY,NX)=0.0_r8
            XFE2FS(N,L,NY,NX)=0.0_r8
            XFE3FS(N,L,NY,NX)=0.0_r8
            XFE4FS(N,L,NY,NX)=0.0_r8
            XFESFS(N,L,NY,NX)=0.0_r8
            XCAOFS(N,L,NY,NX)=0.0_r8
            XCACFS(N,L,NY,NX)=0.0_r8
            XCAHFS(N,L,NY,NX)=0.0_r8
            XCASFS(N,L,NY,NX)=0.0_r8
            XMGOFS(N,L,NY,NX)=0.0_r8
            XMGCFS(N,L,NY,NX)=0.0_r8
            XMGHFS(N,L,NY,NX)=0.0_r8
            XMGSFS(N,L,NY,NX)=0.0_r8
            XNACFS(N,L,NY,NX)=0.0_r8
            XNASFS(N,L,NY,NX)=0.0_r8
            XKASFS(N,L,NY,NX)=0.0_r8
            XH0PFS(N,L,NY,NX)=0.0_r8
            XH3PFS(N,L,NY,NX)=0.0_r8
            XF1PFS(N,L,NY,NX)=0.0_r8
            XF2PFS(N,L,NY,NX)=0.0_r8
            XC0PFS(N,L,NY,NX)=0.0_r8
            XC1PFS(N,L,NY,NX)=0.0_r8
            XC2PFS(N,L,NY,NX)=0.0_r8
            XM1PFS(N,L,NY,NX)=0.0_r8
            XH0BFB(N,L,NY,NX)=0.0_r8
            XH3BFB(N,L,NY,NX)=0.0_r8
            XF1BFB(N,L,NY,NX)=0.0_r8
            XF2BFB(N,L,NY,NX)=0.0_r8
            XC0BFB(N,L,NY,NX)=0.0_r8
            XC1BFB(N,L,NY,NX)=0.0_r8
            XC2BFB(N,L,NY,NX)=0.0_r8
            XM1BFB(N,L,NY,NX)=0.0_r8
            XALFHS(N,L,NY,NX)=0.0_r8
            XFEFHS(N,L,NY,NX)=0.0_r8
            XHYFHS(N,L,NY,NX)=0.0_r8
            XCAFHS(N,L,NY,NX)=0.0_r8
            XMGFHS(N,L,NY,NX)=0.0_r8
            XNAFHS(N,L,NY,NX)=0.0_r8
            XKAFHS(N,L,NY,NX)=0.0_r8
            XOHFHS(N,L,NY,NX)=0.0_r8
            XSOFHS(N,L,NY,NX)=0.0_r8
            XCLFHS(N,L,NY,NX)=0.0_r8
            XC3FHS(N,L,NY,NX)=0.0_r8
            XHCFHS(N,L,NY,NX)=0.0_r8
            XAL1HS(N,L,NY,NX)=0.0_r8
            XAL2HS(N,L,NY,NX)=0.0_r8
            XAL3HS(N,L,NY,NX)=0.0_r8
            XAL4HS(N,L,NY,NX)=0.0_r8
            XALSHS(N,L,NY,NX)=0.0_r8
            XFE1HS(N,L,NY,NX)=0.0_r8
            XFE2HS(N,L,NY,NX)=0.0_r8
            XFE3HS(N,L,NY,NX)=0.0_r8
            XFE4HS(N,L,NY,NX)=0.0_r8
            XFESHS(N,L,NY,NX)=0.0_r8
            XCAOHS(N,L,NY,NX)=0.0_r8
            XCACHS(N,L,NY,NX)=0.0_r8
            XCAHHS(N,L,NY,NX)=0.0_r8
            XCASHS(N,L,NY,NX)=0.0_r8
            XMGOHS(N,L,NY,NX)=0.0_r8
            XMGCHS(N,L,NY,NX)=0.0_r8
            XMGHHS(N,L,NY,NX)=0.0_r8
            XMGSHS(N,L,NY,NX)=0.0_r8
            XNACHS(N,L,NY,NX)=0.0_r8
            XNASHS(N,L,NY,NX)=0.0_r8
            XKASHS(N,L,NY,NX)=0.0_r8
            XH0PHS(N,L,NY,NX)=0.0_r8
            XH3PHS(N,L,NY,NX)=0.0_r8
            XF1PHS(N,L,NY,NX)=0.0_r8
            XF2PHS(N,L,NY,NX)=0.0_r8
            XC0PHS(N,L,NY,NX)=0.0_r8
            XC1PHS(N,L,NY,NX)=0.0_r8
            XC2PHS(N,L,NY,NX)=0.0_r8
            XM1PHS(N,L,NY,NX)=0.0_r8
            XH0BHB(N,L,NY,NX)=0.0_r8
            XH3BHB(N,L,NY,NX)=0.0_r8
            XF1BHB(N,L,NY,NX)=0.0_r8
            XF2BHB(N,L,NY,NX)=0.0_r8
            XC0BHB(N,L,NY,NX)=0.0_r8
            XC1BHB(N,L,NY,NX)=0.0_r8
            XC2BHB(N,L,NY,NX)=0.0_r8
            XM1BHB(N,L,NY,NX)=0.0_r8
8880      CONTINUE
8885    CONTINUE
      ENDIF
8890  CONTINUE
8895  CONTINUE
  end subroutine ResetSaltModelArrays
!------------------------------------------------------------------------------------------

  subroutine SetSoilPropertyAftDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L,K,N,M
  !     begin_execution
  DO 9975 L=NUI(NY,NX),NLI(NY,NX)
    !
    !     AREA,DLYR=lateral(1,2), vertical(3) area,thickness of soil layer
    !     VOLT,VOLX,VOLY=layer volume including,excluding rock,macropores
    !
    IF(BKDS(L,NY,NX).LE.ZERO.AND.DLYR(3,L,NY,NX).LE.ZERO2)THEN
      VOLW(L,NY,NX)=0.0_r8
      VOLI(L,NY,NX)=0.0_r8
    ENDIF
    AREA(1,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
    AREA(2,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
    VOLT(L,NY,NX)=AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)
    VOLX(L,NY,NX)=VOLT(L,NY,NX)*FMPR(L,NY,NX)
    IF(BKDS(L,NY,NX).LE.ZERO)THEN
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
    ENDIF
    !     IF(NX.EQ.1)THEN
    !     WRITE(*,443)'VOLT',I,J,NX,NY,L,VOLT(L,NY,NX)
    !    2,VOLX(L,NY,NX),DLYR(3,L,NY,NX),AREA(3,L,NY,NX)
    !443   FORMAT(A8,5I4,20E14.6)
    !     ENDIF
    !
    !     BKVL=soil mass
    !     C*=concentration,ORGC=SOC,SAND=sand,SILT=silt,CLAY=clay
    !     PTDS=particle density
    !     PTDSNU=particle density of surface layer for use in erosion.f
    !     POROS=porosity used in diffusivity
    !     VOLA,VOLW,VOLI,VOLP=total,water-,ice-,air-filled micropore volume
    !     VOLAH,VOLWH,VOLIH,VOLPH=total,water-,ice-,air-filled macropore volume
    !     EHUM=fraction of microbial decomposition product allocated to humus
    !     EPOC=fraction of SOC decomposition product allocated to POC
    !     SRP=parameter for deviation from linear log-log water retention
    !     PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
    !     PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
    !     FC,WP=water contents at field capacity,wilting point
    !     FCL,WPL=log FC,WP
    !     FCD,PSD=FCL-WPL,log(POROS)-FCL
    !
      BKVL(L,NY,NX)=BKDS(L,NY,NX)*VOLX(L,NY,NX)
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
        CORGC(L,NY,NX)=AMIN1(0.55E+06,ORGC(L,NY,NX)/BKVL(L,NY,NX))
        CSAND(L,NY,NX)=SAND(L,NY,NX)/BKVL(L,NY,NX)
        CSILT(L,NY,NX)=SILT(L,NY,NX)/BKVL(L,NY,NX)
        CCLAY(L,NY,NX)=CLAY(L,NY,NX)/BKVL(L,NY,NX)
      ELSE
        CORGC(L,NY,NX)=0.0_r8
        CSAND(L,NY,NX)=0.0_r8
        CSILT(L,NY,NX)=0.0_r8
        CCLAY(L,NY,NX)=0.0_r8
      ENDIF
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
        CORGCM=AMAX1(0.0,AMIN1(1.0,1.82E-06*CORGC(L,NY,NX)))
        PTDS=1.30*CORGCM+2.66*(1.0-CORGCM)
        IF(L.EQ.NU(NY,NX))THEN
          POROS(L,NY,NX)=AMAX1(POROS(L,NY,NX),1.0-(BKDS(L,NY,NX)/PTDS))
        ELSE
          POROS(L,NY,NX)=1.0-(BKDS(L,NY,NX)/PTDS)
        ENDIF
      ELSE
        PTDS=0.0_r8
        POROS(L,NY,NX)=1.0
      ENDIF
      !     VOLA(L,NY,NX)=AMAX1(POROS(L,NY,NX)*VOLY(L,NY,NX)
      !    2,VOLW(L,NY,NX)+VOLI(L,NY,NX))
      !     VOLAH(L,NY,NX)=AMAX1(FHOL(L,NY,NX)*VOLT(L,NY,NX)
      !    2,VOLWH(L,NY,NX)+VOLIH(L,NY,NX))
      VOLA(L,NY,NX)=POROS(L,NY,NX)*VOLY(L,NY,NX)
      VOLAH(L,NY,NX)=FHOL(L,NY,NX)*VOLT(L,NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
        VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX) &
          -VOLI(L,NY,NX))+AMAX1(0.0,VOLAH(L,NY,NX)-VOLWH(L,NY,NX) &
          -VOLIH(L,NY,NX))
      ELSE
        VOLP(L,NY,NX)=0.0_r8
      ENDIF
      EHUM(L,NY,NX)=0.200+0.333*AMIN1(0.5,CCLAY(L,NY,NX))
      !    2+0.167E-06*CORGC(L,NY,NX)
      !     WRITE(*,3331)'EHUM',I,J,L,EHUM(L,NY,NX),CCLAY(L,NY,NX)
      !    3,CSILT(L,NY,NX),CSAND(L,NY,NX),CORGC(L,NY,NX),ORGC(L,NY,NX)
      !    4,VOLA(L,NY,NX),POROS(L,NY,NX),VOLX(L,NY,NX),BKVL(L,NY,NX)
!3331  FORMAT(A8,3I4,12E12.4)
      EPOC(L,NY,NX)=1.0
      IF(CORGC(L,NY,NX).GT.FORGC)THEN
        SRP(L,NY,NX)=0.25
      ELSEIF(CORGC(L,NY,NX).GT.0.5*FORGC)THEN
        SRP(L,NY,NX)=0.33
      ELSE
        SRP(L,NY,NX)=1.00
      ENDIF
      PSL(L,NY,NX)=LOG(POROS(L,NY,NX))
      IF((ISOIL(1,L,NY,NX).EQ.0.AND.ISOIL(2,L,NY,NX).EQ.0) &
        .OR.DATA(20).EQ.'YES')THEN
        FCL(L,NY,NX)=LOG(FC(L,NY,NX))
        WPL(L,NY,NX)=LOG(WP(L,NY,NX))
        PSD(L,NY,NX)=PSL(L,NY,NX)-FCL(L,NY,NX)
        FCD(L,NY,NX)=FCL(L,NY,NX)-WPL(L,NY,NX)
      ELSE
        !
        !     DEFAULT SOIL HYDROLOGIC PPTYS (FIELD CAPACITY, WILTING POINT)
        !     IF ACTUAL VALUES WERE NOT INPUT TO THE SOIL FILE
        !
        !     THW,THI=initial soil water,ice content from soil file
        !
        IF(DATA(20).EQ.'NO')THEN
          IF(ISOIL(1,L,NY,NX).EQ.1.OR.ISOIL(2,L,NY,NX).EQ.1)THEN
            IF(CORGC(L,NY,NX).LT.FORGW)THEN
              FC(L,NY,NX)=0.2576-0.20*CSAND(L,NY,NX) &
                +0.36*CCLAY(L,NY,NX)+0.60E-06*CORGC(L,NY,NX)
            ELSE
              IF(BKDS(L,NY,NX).LT.0.075)THEN
                FC(L,NY,NX)=0.27
              ELSEIF(BKDS(L,NY,NX).LT.0.195)THEN
                FC(L,NY,NX)=0.62
              ELSE
                FC(L,NY,NX)=0.71
              ENDIF
            ENDIF
            FC(L,NY,NX)=FC(L,NY,NX)/(1.0-FHOL(L,NY,NX))
            FC(L,NY,NX)=AMIN1(0.75*POROS(L,NY,NX),FC(L,NY,NX))
            !     WRITE(*,3332)'FC',IYRC,I,J,L,FC(L,NY,NX),CCLAY(L,NY,NX)
            !    2,CORGC(L,NY,NX),CSAND(L,NY,NX)
!3332  FORMAT(A8,4I6,20E12.4)
            IF(CORGC(L,NY,NX).LT.FORGW)THEN
              WP(L,NY,NX)=0.0260+0.50*CCLAY(L,NY,NX)+0.32E-06*CORGC(L,NY,NX)
            ELSE
              IF(BKDS(L,NY,NX).LT.0.075)THEN
                WP(L,NY,NX)=0.04
              ELSEIF(BKDS(L,NY,NX).LT.0.195)THEN
                WP(L,NY,NX)=0.15
              ELSE
                WP(L,NY,NX)=0.22
              ENDIF
            ENDIF
            WP(L,NY,NX)=WP(L,NY,NX)/(1.0-FHOL(L,NY,NX))
            WP(L,NY,NX)=AMIN1(0.75*FC(L,NY,NX),WP(L,NY,NX))
            !     WRITE(*,3332)'WP',IYRC,I,J,L,WP(L,NY,NX),CCLAY(L,NY,NX)
            !    2,CORGC(L,NY,NX),FC(L,NY,NX)
          ENDIF
          FCL(L,NY,NX)=LOG(FC(L,NY,NX))
          WPL(L,NY,NX)=LOG(WP(L,NY,NX))
          PSD(L,NY,NX)=PSL(L,NY,NX)-FCL(L,NY,NX)
          FCD(L,NY,NX)=FCL(L,NY,NX)-WPL(L,NY,NX)
        ENDIF
        IF(I.EQ.IBEGIN.AND.J.EQ.1.AND.IYRC.EQ.IDATA(9))THEN
          IF(THW(L,NY,NX).GT.1.0.OR.DPTH(L,NY,NX).GE.DTBLZ(NY,NX))THEN
            THETW(L,NY,NX)=POROS(L,NY,NX)
          ELSEIF(test_aeqb(THW(L,NY,NX),1._r8))THEN
            THETW(L,NY,NX)=FC(L,NY,NX)
          ELSEIF(test_aeqb(THW(L,NY,NX),0._r8))THEN
            THETW(L,NY,NX)=WP(L,NY,NX)
          ELSEIF(THW(L,NY,NX).LT.0.0)THEN
            THETW(L,NY,NX)=0.0_r8
          ENDIF
          IF(THI(L,NY,NX).GT.1.0.OR.DPTH(L,NY,NX).GE.DTBLZ(NY,NX))THEN
            THETI(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX) &
              ,POROS(L,NY,NX)-THW(L,NY,NX)))
          ELSEIF(test_aeqb(THI(L,NY,NX),1._r8))THEN
            THETI(L,NY,NX)=AMAX1(0.0,AMIN1(FC(L,NY,NX) &
              ,POROS(L,NY,NX)-THW(L,NY,NX)))
          ELSEIF(test_aeqb(THI(L,NY,NX),0._r8))THEN
            THETI(L,NY,NX)=AMAX1(0.0,AMIN1(WP(L,NY,NX) &
              ,POROS(L,NY,NX)-THW(L,NY,NX)))
          ELSEIF(THI(L,NY,NX).LT.0.0)THEN
            THETI(L,NY,NX)=0.0_r8
          ENDIF
          IF(DATA(20).EQ.'NO')THEN
            VOLW(L,NY,NX)=THETW(L,NY,NX)*VOLX(L,NY,NX)
            VOLWX(L,NY,NX)=VOLW(L,NY,NX)
            VOLWH(L,NY,NX)=THETW(L,NY,NX)*VOLAH(L,NY,NX)
            VOLI(L,NY,NX)=THETI(L,NY,NX)*VOLX(L,NY,NX)
            VOLIH(L,NY,NX)=THETI(L,NY,NX)*VOLAH(L,NY,NX)
            VHCP(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW(L,NY,NX) &
              +VOLWH(L,NY,NX))+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
            THETWZ(L,NY,NX)=THETW(L,NY,NX)
            THETIZ(L,NY,NX)=THETI(L,NY,NX)
          ENDIF
        ENDIF
      ENDIF
      VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX) &
        -VOLI(L,NY,NX))+AMAX1(0.0,VOLAH(L,NY,NX)-VOLWH(L,NY,NX) &
        -VOLIH(L,NY,NX))
      IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        THETP(L,NY,NX)=VOLP(L,NY,NX)/VOLY(L,NY,NX)
      ELSE
        THETP(L,NY,NX)=0.0_r8
      ENDIF
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
        THETY(L,NY,NX)=EXP((PSIMX(NY,NX)-LOG(-PSIHY)) &
          *FCD(L,NY,NX)/PSIMD(NY,NX)+FCL(L,NY,NX))
      ELSE
        THETY(L,NY,NX)=ZERO2
      ENDIF
      !
      !     SATURATED HYDRAULIC CONDUCTIVITY FROM SWC AT SATURATION VS.
      !     -0.033 MPA (MINERAL SOILS) IF NOT ENTERED IN SOIL FILE IN 'READS'
      !
      !     SCNV,SCNH=vertical,lateral saturated hydraulic conductivity
!
      IF(ISOIL(3,L,NY,NX).EQ.1)THEN
        IF(CORGC(L,NY,NX).LT.FORGW)THEN
          THETF=AMIN1(POROS(L,NY,NX),EXP((PSIMS(NY,NX)-LOG(0.033)) &
            *(PSL(L,NY,NX)-FCL(L,NY,NX))/PSISD(NY,NX)+PSL(L,NY,NX)))
          SCNV(L,NY,NX)=1.54*((POROS(L,NY,NX)-THETF)/THETF)**2
        ELSE
          SCNV(L,NY,NX)=0.10+75.0*1.0E-15**BKDS(L,NY,NX)
          SCNV(L,NY,NX)=SCNV(L,NY,NX)*FMPR(L,NY,NX)
        ENDIF
        !     WRITE(*,3332)'SCNV',IYRC,I,J,L,SCNV(L,NY,NX),POROS(L,NY,NX)
        !    2,THETF,FMPR(L,NY,NX),PSIMS(NY,NX),LOG(0.033)
        !    3,PSL(L,NY,NX),FCL(L,NY,NX),PSISD(NY,NX)
      ENDIF
      IF(ISOIL(4,L,NY,NX).EQ.1)THEN
        IF(CORGC(L,NY,NX).LT.FORGW)THEN
          THETF=AMIN1(POROS(L,NY,NX),EXP((PSIMS(NY,NX)-LOG(0.033)) &
            *(PSL(L,NY,NX)-FCL(L,NY,NX))/PSISD(NY,NX)+PSL(L,NY,NX)))
          SCNH(L,NY,NX)=1.54*((POROS(L,NY,NX)-THETF)/THETF)**2
        ELSE
          SCNH(L,NY,NX)=0.10+75.0*1.0E-15**BKDS(L,NY,NX)
          SCNH(L,NY,NX)=SCNH(L,NY,NX)*FMPR(L,NY,NX)
        ENDIF
        !     WRITE(*,3332)'SCNH',IYRC,I,J,L,SCNH(L,NY,NX),POROS(L,NY,NX)
        !    2,THETF,FMPR(L,NY,NX)
      ENDIF
      !     WRITE(*,3333)'PPTYS',I,J,NX,NY,L,IBEGIN,IYRC,IDATA(9)
      !    2,ISOIL(1,L,NY,NX),ISOIL(2,L,NY,NX)
      !    3,ISOIL(3,L,NY,NX),ISOIL(4,L,NY,NX)
      !    3,SCNV(L,NY,NX),SCNH(L,NY,NX),POROS(L,NY,NX)
      !    2,FC(L,NY,NX),WP(L,NY,NX),BKDS(L,NY,NX),THW(L,NY,NX)
      !    3,THETW(L,NY,NX),THI(L,NY,NX),THETI(L,NY,NX),VOLI(L,NY,NX)
!3333  FORMAT(A8,12I6,20E12.4)
      !
      !     HYDRAULIC CONDUCTIVITY FUNCTION FROM KSAT AND SOIL WATER RELEASE CURVE
      !
      !     THETK,PSISK=micropore class water content,potential
      !     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
      !
      !     IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      SUM2=0.0_r8
      DO 1320 K=1,100
        XK=K-1
        THETK(K)=POROS(L,NY,NX)-(XK/100.0*POROS(L,NY,NX))
        IF(THETK(K).LT.FC(L,NY,NX))THEN
          PSISK(K)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
            +((FCL(L,NY,NX)-LOG(THETK(K))) &
            /FCD(L,NY,NX)*PSIMD(NY,NX))))
        ELSEIF(THETK(K).LT.POROS(L,NY,NX)-DTHETW)THEN
          PSISK(K)=-EXP(PSIMS(NY,NX) &
            +(((PSL(L,NY,NX)-LOG(THETK(K))) &
            /PSD(L,NY,NX))**SRP(L,NY,NX)*PSISD(NY,NX)))
        ELSE
          PSISK(K)=PSISE(L,NY,NX)
        ENDIF
        SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
1320  CONTINUE
      DO 1335 K=1,100
      SUM1=0.0_r8
      XK=K-1
      YK=((100.0-XK)/100.0)**1.33
      DO 1330 M=K,100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2)
1330  CONTINUE
      DO 1340 N=1,3
      IF(N.EQ.3)THEN
      HCND(N,K,L,NY,NX)=SCNV(L,NY,NX)*YK*SUM1/SUM2
      IF(K.GT.1.AND.PSISK(K).LT.PSISA(L,NY,NX) &
      .AND.PSISK(K-1).GE.PSISA(L,NY,NX))THEN
      THETS(L,NY,NX)=THETK(K)
      ENDIF
!     WRITE(*,3536)'PSI',L,K,THETK(K),PSISK(K),HCND(N,K,L,NY,NX)
!    2,PSL(L,NY,NX),LOG(THETK(K)),POROS(L,NY,NX),FC(L,NY,NX),WP(L,NY,NX)
!    3,SRP(L,NY,NX),THETS(L,NY,NX),PSISA(L,NY,NX)
!3536  FORMAT(A8,2I4,12E12.4)
      ELSE
      HCND(N,K,L,NY,NX)=SCNH(L,NY,NX)*YK*SUM1/SUM2
      ENDIF
1340  CONTINUE
1335  CONTINUE
!2340  CONTINUE
!2335  CONTINUE
!     ENDIF
!
!     SOIL MACROPORE DIMENSIONS AND CONDUCTIVITY FROM MACROPORE FRACTION
!     ENTERED IN 'READS'
!
!     PHOL,NHOL,HRAD=path length between, number,radius of macropores
!     CNDH=macropore hydraulic conductivity
!
      HRAD(L,NY,NX)=0.5E-03
      NHOL(L,NY,NX)=INT(VOLAH(L,NY,NX)/(PICON*HRAD(L,NY,NX)**2 &
      *VOLTI(L,NY,NX)))
      IF(NHOL(L,NY,NX).GT.0.0)THEN
      PHOL(L,NY,NX)=1.0/(SQRT(PICON*NHOL(L,NY,NX)))
      ELSE
      PHOL(L,NY,NX)=1.0
      ENDIF
      VISCWL=VISCW*EXP(0.533-0.0267*TCS(L,NY,NX))
      CNDH(L,NY,NX)=3.6E+03*PICON*NHOL(L,NY,NX)*HRAD(L,NY,NX)**4 &
      /(8.0*VISCWL)
!
!     SOIL HEAT CAPACITY AND THERMAL CONDUCTIVITY OF SOLID PHASE
!     FROM SOC AND TEXTURE
!
!     VORGC,VMINL,VSAND=volume fractions of SOC,mineral,sand
!     STC,DTC=weighted thermal conductivity of soil solid component
!
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VORGC=CORGCM*BKDS(L,NY,NX)/PTDS
      VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*BKDS(L,NY,NX)/PTDS
      VSAND=CSAND(L,NY,NX)*BKDS(L,NY,NX)/PTDS
      STC(L,NY,NX)=(1.253*VORGC*9.050E-04+0.514*VMINL*1.056E-02 &
      +0.386*VSAND*2.112E-02)*FMPR(L,NY,NX) &
      +0.514*ROCK(L,NY,NX)*1.056E-02
      DTC(L,NY,NX)=(1.253*VORGC+0.514*VMINL+0.386*VSAND) &
      *FMPR(L,NY,NX)+0.514*ROCK(L,NY,NX)
      ELSE
      STC(L,NY,NX)=0.0_r8
      DTC(L,NY,NX)=0.0_r8
      ENDIF
9975  CONTINUE
      end subroutine SetSoilPropertyAftDisturbance
!------------------------------------------------------------------------------------------

      subroutine ResetSurfResidualProperty(NY,NX)
      implicit none
      integer, intent(in) :: NY,NX
!     begin_execution
!
!     FCR=litter water content at -0.01 MPa
!     THETY=litter hygroscopic water content
!
      CORGC(0,NY,NX)=0.55E+06
!
!     SOIL SURFACE WATER STORAGE CAPACITY
!
!     IDTBL=water table flag from site file
!     DTBLX,DTBLZ=current,initial natural water table depth
!     DTBLY,DTBLD=current,initial artificial water table depth
!     ZS,ZW=soil,water surface roughness
!     VOLWD=soil surface water retention capacity
!     VOLWG=VOLWD accounting for above-ground water table
!     EHUM=fraction of microbial decompn product allocated to surface humus
!     EPOC=fraction of SOC decomposition product allocated to surface POC
!
      IF(IDTBL(NY,NX).LE.1.OR.IDTBL(NY,NX).EQ.3)THEN
      DTBLX(NY,NX)=DTBLZ(NY,NX)
      ELSEIF(IDTBL(NY,NX).EQ.2.OR.IDTBL(NY,NX).EQ.4)THEN
      DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      ENDIF
      IF(IDTBL(NY,NX).EQ.3.OR.IDTBL(NY,NX).EQ.4)THEN
      DTBLY(NY,NX)=DTBLD(NY,NX)
      ENDIF
!     IF(J.EQ.24)THEN
!     WRITE(*,1114)'DTBLX',I,J,IYRC,NX,NY,IDTBL(NY,NX),DTBLX(NY,NX)
!    2,DTBLZ(NY,NX),CDPTH(NU(NY,NX)-1,NY,NX)
!    3,(DPTH(L,NY,NX),L=1,NL(NY,NX))
!    3,BKDS(NU(NY,NX),NY,NX),BKDSI(NU(NY,NX),NY,NX)
!    4,DLYR(3,NU(NY,NX),NY,NX),DLYRI(3,NU(NY,NX),NY,NX)
!1114  FORMAT(A8,6I4,30E12.4)
!     ENDIF
      IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
      ZS(NY,NX)=0.020
      ELSE
      ZS(NY,NX)=ZW
      ENDIF
      VOLWD(NY,NX)=AMAX1(0.001,0.112*ZS(NY,NX)+3.10*ZS(NY,NX)**2 &
      -0.012*ZS(NY,NX)*SLOPE(0,NY,NX))*AREA(3,NU(NY,NX),NY,NX)
      VOLWG(NY,NX)=AMAX1(VOLWD(NY,NX) &
      ,-(DTBLX(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)) &
      *AREA(3,NU(NY,NX),NY,NX))
!     WRITE(*,6631)'VOLWG',I,J,NY,NX,NU(NY,NX),VOLWG(NY,NX)
!    2,DTBLX(NY,NX),CDPTH(NU(NY,NX)-1,NY,NX)
!    3,AREA(3,NU(NY,NX),NY,NX),VOLWD(NY,NX),ZS(NY,NX),SLOPE(0,NY,NX)
!6631  FORMAT(A8,5I4,12E12.4)
      DPTH(NU(NY,NX),NY,NX)=CDPTH(NU(NY,NX),NY,NX) &
      -0.5*DLYR(3,NU(NY,NX),NY,NX)
      IF(BKVL(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX) &
      /BKVL(NU(NY,NX),NY,NX)
      CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX) &
      /BKVL(NU(NY,NX),NY,NX)
      CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX) &
      /BKVL(NU(NY,NX),NY,NX)
      ELSE
      CCLAY(NU(NY,NX),NY,NX)=0.0_r8
      CSILT(NU(NY,NX),NY,NX)=0.0_r8
      CSAND(NU(NY,NX),NY,NX)=0.0_r8
      ENDIF
      EHUM(0,NY,NX)=0.200+0.333*AMIN1(0.5,CCLAY(NU(NY,NX),NY,NX))
!    2+0.167E-06*CORGC(NU(NY,NX),NY,NX)
      EPOC(0,NY,NX)=0.150
      end subroutine ResetSurfResidualProperty
!------------------------------------------------------------------------------------------

      subroutine SetLiterSoilPropAftDisturbance(I,J,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NY,NX

      integer :: K,M
!     begin_execution
!     IFLGS=disturbance flag
!     BKDS,BKDSI=current,initial bulk density
!
!     WRITE(*,1116)'IFLGS',IYRC,I,J,IFLGS(NY,NX)
!1116  FORMAT(A8,4I6)
      IF(IFLGS(NY,NX).NE.0)THEN
      IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      BKDS(0,NY,NX)=BKVL(0,NY,NX)/VOLT(0,NY,NX)
      ELSE
      BKDS(0,NY,NX)=BKRS(1)
      ENDIF
      THETY(0,NY,NX)=EXP((PSIMX(NY,NX)-LOG(-PSIHY)) &
      *FCD(0,NY,NX)/PSIMD(NY,NX)+FCL(0,NY,NX))
      SUM2=0.0_r8
      DO 1220 K=1,100
      XK=K-1
      THETK(K)=POROS0(NY,NX)-(XK/100.0*POROS0(NY,NX))
      IF(THETK(K).LT.FC(0,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
      +((FCL(0,NY,NX)-LOG(THETK(K))) &
      /FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETK(K).LT.POROS0(NY,NX))THEN
      PSISK(K)=-EXP(PSIMS(NY,NX) &
      +(((PSL(0,NY,NX)-LOG(THETK(K))) &
      /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISK(K)=PSISE(0,NY,NX)
      ENDIF
      SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
1220  CONTINUE
      DO 1235 K=1,100
      SUM1=0.0_r8
      XK=K-1
      YK=((100.0-XK)/100.0)**1.33
      DO 1230 M=K,100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2)
1230  CONTINUE
      HCND(3,K,0,NY,NX)=SCNV(0,NY,NX)*YK*SUM1/SUM2
      HCND(1,K,0,NY,NX)=0.0_r8
      HCND(2,K,0,NY,NX)=0.0_r8
      IF(K.GT.1.AND.PSISK(K).LT.PSISA(0,NY,NX) &
      .AND.PSISK(K-1).GE.PSISA(0,NY,NX))THEN
      THETS(0,NY,NX)=THETK(K)
      ENDIF
!     WRITE(*,3534)'PSI0',K,THETK(K),PSISK(K),HCND(3,K,0,NY,NX)
!    2,PSL(0,NY,NX),LOG(THETK(K)),POROS0(NY,NX),FC(0,NY,NX),WP(0,NY,NX)
!    3,SRP(0,NY,NX),VOLWRX(NY,NX),VOLR(NY,NX),THETS(0,NY,NX)
!    4,PSISA(0,NY,NX)
!3534  FORMAT(A8,1I4,20E12.4)
1235  CONTINUE
!
!     RESET SOIL PHYSICAL PROPERTIES (DENSITY, TEXTURE)
!     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
!
      call SetSoilPropertyAftDisturbance(I,J,NY,NX)
!
!     SURFACE RESIDUE PROPERTIES
      call ResetSurfResidualProperty(NY,NX)
!
!     IFLGS=reset disturbance flag
!
      IFLGS(NY,NX)=0
      ENDIF
      end subroutine SetLiterSoilPropAftDisturbance
!------------------------------------------------------------------------------------------

      subroutine SetHourlyAccumulators(NY,NX)
!     implicit none
      integer, intent(in) :: NX,NY

      integer :: L
!     begin_execution

      UCO2S(NY,NX)=0.0_r8
      TOMT(NY,NX)=0.0_r8
      TONT(NY,NX)=0.0_r8
      TOPT(NY,NX)=0.0_r8
      UVOLW(NY,NX)=0.0_r8
      URSDC(NY,NX)=0.0_r8
      UORGC(NY,NX)=0.0_r8
      URSDN(NY,NX)=0.0_r8
      UORGN(NY,NX)=0.0_r8
      URSDP(NY,NX)=0.0_r8
      UORGP(NY,NX)=0.0_r8
      UNH4(NY,NX)=0.0_r8
      UNO3(NY,NX)=0.0_r8
      UPO4(NY,NX)=0.0_r8
      UPP4(NY,NX)=0.0_r8
      UION(NY,NX)=0.0_r8
      HVOLO(NY,NX)=0.0_r8
      HCO2G(NY,NX)=0.0_r8
      HCH4G(NY,NX)=0.0_r8
      HOXYG(NY,NX)=0.0_r8
      HN2GG(NY,NX)=0.0_r8
      HN2OG(NY,NX)=0.0_r8
      HNH3G(NY,NX)=0.0_r8
      FLWR(NY,NX)=0.0_r8
      HFLWR(NY,NX)=0.0_r8
      THAWR(NY,NX)=0.0_r8
      HTHAWR(NY,NX)=0.0_r8
      HEATI(NY,NX)=0.0_r8
      HEATS(NY,NX)=0.0_r8
      HEATE(NY,NX)=0.0_r8
      HEATV(NY,NX)=0.0_r8
      HEATH(NY,NX)=0.0_r8
      TEVAPG(NY,NX)=0.0_r8
      XCODFS(NY,NX)=0.0_r8
      XCHDFS(NY,NX)=0.0_r8
      XOXDFS(NY,NX)=0.0_r8
      XNGDFS(NY,NX)=0.0_r8
      XN2DFS(NY,NX)=0.0_r8
      XN3DFS(NY,NX)=0.0_r8
      XNBDFS(NY,NX)=0.0_r8
      XHGDFS(NY,NX)=0.0_r8
      XCODFR(NY,NX)=0.0_r8
      XCHDFR(NY,NX)=0.0_r8
      XOXDFR(NY,NX)=0.0_r8
      XNGDFR(NY,NX)=0.0_r8
      XN2DFR(NY,NX)=0.0_r8
      XN3DFR(NY,NX)=0.0_r8
      XHGDFR(NY,NX)=0.0_r8
      TVOLWP(NY,NX)=0.0_r8
      TVOLWC(NY,NX)=0.0_r8
      TFLWCI(NY,NX)=0.0_r8
      TFLWC(NY,NX)=0.0_r8
      TEVAPP(NY,NX)=0.0_r8
      TEVAPC(NY,NX)=0.0_r8
      THFLXC(NY,NX)=0.0_r8
      TENGYC(NY,NX)=0.0_r8
      TCO2Z(NY,NX)=0.0_r8
      TOXYZ(NY,NX)=0.0_r8
      TCH4Z(NY,NX)=0.0_r8
      TN2OZ(NY,NX)=0.0_r8
      TNH3Z(NY,NX)=0.0_r8
      TH2GZ(NY,NX)=0.0_r8
      ZCSNC(NY,NX)=0.0_r8
      ZZSNC(NY,NX)=0.0_r8
      ZPSNC(NY,NX)=0.0_r8
      WTSTGT(NY,NX)=0.0_r8
      PPT(NY,NX)=0.0_r8
      DO 9865 L=1,JS
      FLSW(L,NY,NX)=0.0_r8
      FLSWH(L,NY,NX)=0.0_r8
      HFLSW(L,NY,NX)=0.0_r8
      FLSWR(L,NY,NX)=0.0_r8
      HFLSWR(L,NY,NX)=0.0_r8
      XFLWS(L,NY,NX)=0.0_r8
      XFLWW(L,NY,NX)=0.0_r8
      XFLWI(L,NY,NX)=0.0_r8
      XHFLWW(L,NY,NX)=0.0_r8
      XWFLXS(L,NY,NX)=0.0_r8
      XWFLXI(L,NY,NX)=0.0_r8
      XTHAWW(L,NY,NX)=0.0_r8
      XCOBLS(L,NY,NX)=0.0_r8
      XCHBLS(L,NY,NX)=0.0_r8
      XOXBLS(L,NY,NX)=0.0_r8
      XNGBLS(L,NY,NX)=0.0_r8
      XN2BLS(L,NY,NX)=0.0_r8
      XN4BLW(L,NY,NX)=0.0_r8
      XN3BLW(L,NY,NX)=0.0_r8
      XNOBLW(L,NY,NX)=0.0_r8
      XH1PBS(L,NY,NX)=0.0_r8
      XH2PBS(L,NY,NX)=0.0_r8
      IF(ISALTG.NE.0)THEN
      XALBLS(L,NY,NX)=0.0_r8
      XFEBLS(L,NY,NX)=0.0_r8
      XHYBLS(L,NY,NX)=0.0_r8
      XCABLS(L,NY,NX)=0.0_r8
      XMGBLS(L,NY,NX)=0.0_r8
      XNABLS(L,NY,NX)=0.0_r8
      XKABLS(L,NY,NX)=0.0_r8
      XOHBLS(L,NY,NX)=0.0_r8
      XSOBLS(L,NY,NX)=0.0_r8
      XCLBLS(L,NY,NX)=0.0_r8
      XC3BLS(L,NY,NX)=0.0_r8
      XHCBLS(L,NY,NX)=0.0_r8
      XAL1BS(L,NY,NX)=0.0_r8
      XAL2BS(L,NY,NX)=0.0_r8
      XAL3BS(L,NY,NX)=0.0_r8
      XAL4BS(L,NY,NX)=0.0_r8
      XALSBS(L,NY,NX)=0.0_r8
      XFE1BS(L,NY,NX)=0.0_r8
      XFE2BS(L,NY,NX)=0.0_r8
      XFE3BS(L,NY,NX)=0.0_r8
      XFE4BS(L,NY,NX)=0.0_r8
      XFESBS(L,NY,NX)=0.0_r8
      XCAOBS(L,NY,NX)=0.0_r8
      XCACBS(L,NY,NX)=0.0_r8
      XCAHBS(L,NY,NX)=0.0_r8
      XCASBS(L,NY,NX)=0.0_r8
      XMGOBS(L,NY,NX)=0.0_r8
      XMGCBS(L,NY,NX)=0.0_r8
      XMGHBS(L,NY,NX)=0.0_r8
      XMGSBS(L,NY,NX)=0.0_r8
      XNACBS(L,NY,NX)=0.0_r8
      XNASBS(L,NY,NX)=0.0_r8
      XKASBS(L,NY,NX)=0.0_r8
      XH0PBS(L,NY,NX)=0.0_r8
      XH3PBS(L,NY,NX)=0.0_r8
      XF1PBS(L,NY,NX)=0.0_r8
      XF2PBS(L,NY,NX)=0.0_r8
      XC0PBS(L,NY,NX)=0.0_r8
      XC1PBS(L,NY,NX)=0.0_r8
      XC2PBS(L,NY,NX)=0.0_r8
      XM1PBS(L,NY,NX)=0.0_r8
      ENDIF
9865  CONTINUE
      end subroutine SetHourlyAccumulators
!------------------------------------------------------------------------------------------

      subroutine SetArrays4PlantSoilTransfer(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX

      integer :: K,M
!     begin_execution

      DO 9950 K=0,1
      DO M=1,4
      CSNT(M,K,L,NY,NX)=0.0_r8
      ZSNT(M,K,L,NY,NX)=0.0_r8
      PSNT(M,K,L,NY,NX)=0.0_r8
      enddo
9950  CONTINUE
      DO 7775 K=0,4
      XOQCS(K,L,NY,NX)=0.0_r8
      XOQNS(K,L,NY,NX)=0.0_r8
      XOQPS(K,L,NY,NX)=0.0_r8
      XOQAS(K,L,NY,NX)=0.0_r8
7775  CONTINUE
      XZHYS(L,NY,NX)=0.0_r8
      TRN4S(L,NY,NX)=0.0_r8
      TRN3S(L,NY,NX)=0.0_r8
      TRN3G(L,NY,NX)=0.0_r8
      TRNO3(L,NY,NX)=0.0_r8
      TRNO2(L,NY,NX)=0.0_r8
      TRH1P(L,NY,NX)=0.0_r8
      TRH2P(L,NY,NX)=0.0_r8
      TRXN4(L,NY,NX)=0.0_r8
      TRXH0(L,NY,NX)=0.0_r8
      TRXH1(L,NY,NX)=0.0_r8
      TRXH2(L,NY,NX)=0.0_r8
      TRX1P(L,NY,NX)=0.0_r8
      TRX2P(L,NY,NX)=0.0_r8
      TRALPO(L,NY,NX)=0.0_r8
      TRFEPO(L,NY,NX)=0.0_r8
      TRCAPD(L,NY,NX)=0.0_r8
      TRCAPH(L,NY,NX)=0.0_r8
      TRCAPM(L,NY,NX)=0.0_r8
      TUPWTR(L,NY,NX)=0.0_r8
      TUPHT(L,NY,NX)=0.0_r8
      XCODFG(L,NY,NX)=0.0_r8
      XCHDFG(L,NY,NX)=0.0_r8
      XOXDFG(L,NY,NX)=0.0_r8
      XNGDFG(L,NY,NX)=0.0_r8
      XN2DFG(L,NY,NX)=0.0_r8
      XN3DFG(L,NY,NX)=0.0_r8
      XNBDFG(L,NY,NX)=0.0_r8
      XHGDFG(L,NY,NX)=0.0_r8
      IF(L.GE.NU(NY,NX))THEN
      DO 195 K=0,4
      TDFOMC(K,L,NY,NX)=0.0_r8
      TDFOMN(K,L,NY,NX)=0.0_r8
      TDFOMP(K,L,NY,NX)=0.0_r8
195   CONTINUE
      ENDIF
      DO 9795 M=1,NPH
      ROXSK(M,L,NY,NX)=0.0_r8
9795  CONTINUE
      end subroutine SetArrays4PlantSoilTransfer
!------------------------------------------------------------------------------------------

      subroutine GetOutput4WaterTableDepth(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX

      integer :: LL
!     begin_execution
!     IDTBL=water table flag from site file
!     THETPZ,THETPW=current,minimum air-filled, porosity for water table
!     DPTH,DTBLX=depth of soil layer midpoint, water table
!     PSIS1=water potential in hydraulic equilibrium with layer below
!     THETW1,THETWP=water content at PSIS1,minimum SWC for water table
!     DPTHT=water table depth
!
      IF(IDTBL(NY,NX).NE.0)THEN
      IF(IFLGY.EQ.0)THEN
      IF(THETPZ(L,NY,NX).LT.THETPW.OR.L.EQ.NL(NY,NX))THEN
      IFLGY=1
      IF(DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
      DO 5705 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
      IF(THETPZ(LL,NY,NX).GE.THETPW.AND.LL.NE.NL(NY,NX))THEN
      IFLGY=0
      GO TO 5706
      ELSEIF(DPTH(LL,NY,NX).GE.DTBLX(NY,NX))THEN
      GO TO 5706
      ENDIF
5705  CONTINUE
      ENDIF
5706  CONTINUE
      IF(IFLGY.EQ.1)THEN
      IF(THETPZ(L,NY,NX).GE.THETPW.AND.L.NE.NL(NY,NX))THEN
      PSIS1=PSISM(L+1,NY,NX)-0.0098*(DPTH(L+1,NY,NX)-DPTH(L,NY,NX))
      THETWM=THETWP*POROS(L,NY,NX)
      THETW1=AMIN1(THETWM,EXP((PSIMS(NY,NX)-LOG(-PSIS1)) &
      *PSD(L,NY,NX)/PSISD(NY,NX)+PSL(L,NY,NX)))
      IF(THETWM.GT.THETW1)THEN
      THETPX=AMIN1(1.0,AMAX1(0.0,(THETWM-THETW(L,NY,NX)) &
      /(THETWM-THETW1)))
      DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)*(1.0-THETPX)
      ELSE
      DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
      ENDIF
      ELSEIF(L.GT.NU(NY,NX))THEN
      PSIS1=PSISM(L,NY,NX)-0.0098*(DPTH(L,NY,NX)-DPTH(L-1,NY,NX))
      THETWM=THETWP*POROS(L-1,NY,NX)
      THETW1=AMIN1(THETWM,EXP((PSIMS(NY,NX)-LOG(-PSIS1)) &
      *PSD(L-1,NY,NX)/PSISD(NY,NX)+PSL(L-1,NY,NX)))
      IF(THETWM.GT.THETW1)THEN
      THETPX=AMIN1(1.0,AMAX1(0.0,(THETWM-THETW(L-1,NY,NX)) &
      /(THETWM-THETW1)))
      DPTHT(NY,NX)=CDPTH(L-1,NY,NX)-DLYR(3,L-1,NY,NX)*(1.0-THETPX)
      ELSE
      DPTHT(NY,NX)=CDPTH(L-1,NY,NX)-DLYR(3,L-1,NY,NX)
      ENDIF
      ELSE
      DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
!     IF(NX.EQ.3.AND.NY.EQ.3)THEN
!     WRITE(*,5353)'DPTHT',I,J,NX,NY,L,LL,IFLGY,DPTHT(NY,NX)
!    2,CDPTH(L,NY,NX),DLYR(3,L,NY,NX),DPTH(L,NY,NX),DPTH(L-1,NY,NX)
!    3,PSIS1,PSISM(L,NY,NX),THETWM,THETW1,THETW(L-1,NY,NX)
!    4,THETPZ(L,NY,NX),THETPW
!5353  FORMAT(A8,7I4,30E12.4)
!     ENDIF
      ENDIF
      end subroutine GetOutput4WaterTableDepth
!------------------------------------------------------------------------------------------

  subroutine SetSurfaceProperty4SedErosion(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  !     begin_execution
  !
  !     DETS=soil detachability from rainfall impact
  !     D50=average particle size
  !     CER,XER=parameters for runoff transport capacity
  !     ZD50=particle size effect on surface roughness
  !     VLS=hourly sinking rate
  !     COHS=soil cohesion
  !     DETE=soil detachability
  !     ZM=surface roughness used in runoff velocity calculation in watsub.f
  !
  BKVLNU(NY,NX)=AMAX1(0.0_r8,BKVLNM(NY,NX)+1.82E-06_r8*ORGC(NU(NY,NX),NY,NX))
  !     WRITE(*,2423)'BKVLH',I,J,NX,NY,BKVL(NU(N2,N1),N2,N1)
  !    2,BKVLNM(N2,N1),BKVLNU(N2,N1),ORGR(NU(N2,N1),N2,N1)
  !    3,ORGC(NU(N2,N1),N2,N1),SAND(NU(N2,N1),N2,N1)
  !    3,SILT(NU(N2,N1),N2,N1),CLAY(NU(N2,N1),N2,N1)
!2423  FORMAT(A8,4I4,20E12.4)
  BKVLNX=SAND(NU(NY,NX),NY,NX)+SILT(NU(NY,NX),NY,NX) &
    +CLAY(NU(NY,NX),NY,NX)+1.82E-06*ORGC(NU(NY,NX),NY,NX)
  IF(BKVLNX.GT.ZEROS(NY,NX))THEN
    CORGM=1.82E-06_r8*ORGC(NU(NY,NX),NY,NX)/BKVLNX
    CORGC(NU(NY,NX),NY,NX)=0.55E+06_r8*CORGM
    CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/BKVLNX
    CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/BKVLNX
    CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/BKVLNX
  ELSE
    CORGM=0.0_r8
    CORGC(NU(NY,NX),NY,NX)=0.0_r8
    CSAND(NU(NY,NX),NY,NX)=0.0_r8
    CSILT(NU(NY,NX),NY,NX)=1.0
    CCLAY(NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
    D50=1.0_r8*CCLAY(NU(NY,NX),NY,NX)+10.0_r8*CSILT(NU(NY,NX),NY,NX) &
      +100.0_r8*CSAND(NU(NY,NX),NY,NX)+100.0_r8*CORGM
    ZD50=0.041*(1.0E-06_r8*D50)**0.167_r8
    ZM(NY,NX)=ZS(NY,NX)+ZD50+1.0_r8*VOLR(NY,NX)/AREA(3,0,NY,NX)
    CER(NY,NX)=((D50+5.0_r8)/0.32_r8)**(-0.6_r8)
    XER(NY,NX)=((D50+5.0_r8)/300.0_r8)**0.25_r8
    DETS(NY,NX)=1.0E-06_r8*(1.0_r8+2.0_r8*(1.0_r8-CSILT(NU(NY,NX),NY,NX)-CORGM))
    COHS=2.0_r8+10.0_r8*(CCLAY(NU(NY,NX),NY,NX)+CORGM) &
      +5.0_r8*(1.0_r8-EXP(-2.0E-06_r8*RTDNT(NU(NY,NX),NY,NX)))
    DETE(NY,NX)=0.79_r8*EXP(-0.85_r8*AMAX1(1.0_r8,COHS))
    PTDSNU(NY,NX)=1.30_r8*CORGM+2.66_r8*(1.0_r8-CORGM)
    VISCWL=VISCW*EXP(0.533_r8-0.0267_r8*TCS(0,NY,NX))
    VLS(NY,NX)=3.6E+03_r8*9.8_r8*(PTDSNU(NY,NX)-1.0_r8) &
      *(1.0E-06_r8*D50)**2/(18.0_r8*VISCWL)
!     WRITE(*,1118)'COHS',I,J,NX,NY,NU(NY,NX),COHS,DETE(NY,NX)
!    2,ZM(NY,NX),VLS(NY,NX),D50,ZD50,PTDSNU(NY,NX)
!    3,RTDNT(NU(NY,NX),NY,NX),VOLR(NY,NX)/AREA(3,0,NY,NX)
!    3,ORGC(0,NY,NX)
!    3,CCLAY(NU(NY,NX),NY,NX),CSILT(NU(NY,NX),NY,NX)
!    4,CSAND(NU(NY,NX),NY,NX),CORGM,CORGC(NU(NY,NX),NY,NX)
!    5,BKVL(NU(NY,NX),NY,NX),BKVLNX
!    3,VISCWL,TCS(0,NY,NX)
!1118  FORMAT(A8,5I4,20E12.4)
  ENDIF
  end subroutine SetSurfaceProperty4SedErosion
!------------------------------------------------------------------------------------------

  subroutine UpdateTotalSOC(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: K,M,N,NGL
  !     begin_execution
  !
  !     TOTAL SOC FOR CALCULATING CHANGES IN SOC CALCULATED IN NITRO.F
  !
  !     OMC=microbial biomass, ORC=microbial residue
  !     OQC,OQCH=DOC in micropores,macropores
  !     OQA,OQAH=acetate in micropores,macropores
  !     OHC,OHA=adsorbed SOC,acetate
  !     OSC=SOC(K=0:woody litter, K=1:non-woody litter,
  !     K=2:manure, K=3:POC, K=4:humus)
  !
  DC=0.0_r8
  OC=0.0_r8
  DO 7970 K=0,5
    DO 7950 N=1,7
      DO NGL=1,JG
        DO  M=1,3
          OC=OC+OMC(M,NGL,N,K,L,NY,NX)
        enddo
      ENDDO
7950  CONTINUE
7970  CONTINUE
  DO 7900 K=0,4
    DO 7920 M=1,2
      OC=OC+ORC(M,K,L,NY,NX)
7920  CONTINUE
    OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
      +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
    DO 7910 M=1,4
      OC=OC+OSC(M,K,L,NY,NX)
7910  CONTINUE
7900  CONTINUE
  ORGCX(L,NY,NX)=OC
  end subroutine UpdateTotalSOC
!------------------------------------------------------------------------------------------

  subroutine GetSoilHydraulicVars(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: K
  ! begin_execution
  ! WATER POTENTIALS
  !
  ! FC,WP=water contents at field capacity,wilting point,saturation
  ! PSISM,PSISE=matric,saturation water potential
  ! SRP=parameter for deviation from linear log-log water retention
  ! FC,WP=water contents at field capacity,wilting point from soil file
  ! FCL,WPL=log FC,WP
  ! FCD,PSD=FCL-WPL,log(POROS)-FCL
  ! FCI,WPI=FC,WP of ice
  ! THETIX=ice concentration
!
  IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX) &
    .AND.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
    THETW1=AMAX1(0.0,AMIN1(POROS(L,NY,NX) &
      ,VOLW(L,NY,NX)/VOLY(L,NY,NX)))
    IF(THETW1.LT.FC(L,NY,NX))THEN
      PSISM(L,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
        +((FCL(L,NY,NX)-LOG(THETW1)) &
        /FCD(L,NY,NX)*PSIMD(NY,NX))))
    ELSEIF(THETW1.LT.POROS(L,NY,NX)-DTHETW)THEN
      PSISM(L,NY,NX)=-EXP(PSIMS(NY,NX) &
        +(((PSL(L,NY,NX)-LOG(THETW1)) &
        /PSD(L,NY,NX))**SRP(L,NY,NX)*PSISD(NY,NX)))
    ELSE
      PSISM(L,NY,NX)=PSISE(L,NY,NX)
    ENDIF
  ELSEIF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX).and.THETI(L,NY,NX)>ZEROS2(NY,NX))THEN
    FCX=FCI*THETI(L,NY,NX)
    WPX=WPI*THETI(L,NY,NX)
    FCLX=LOG(FCX)
    WPLX=LOG(WPX)
    PSDX=PSL(L,NY,NX)-FCLX
    FCDX=FCLX-WPLX
    IF(THETW(L,NY,NX).LT.FCX)THEN
      PSISM(L,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
        +((FCLX-LOG(THETW(L,NY,NX))) &
        /FCDX*PSIMD(NY,NX))))
    ELSEIF(THETW(L,NY,NX).LT.POROS(L,NY,NX)-DTHETW)THEN
      PSISM(L,NY,NX)=-EXP(PSIMS(NY,NX) &
        +(((PSL(L,NY,NX)-LOG(THETW(L,NY,NX))) &
        /PSDX)*PSISD(NY,NX)))
    ELSE
      PSISM(L,NY,NX)=PSISE(L,NY,NX)
    ENDIF
  ELSE
    PSISM(L,NY,NX)=PSISE(L,NY,NX)
  ENDIF
!     WRITE(*,443)'PSISM',I,J,NX,NY,L,PSISM(L,NY,NX)
!    2,THETW(L,NY,NX),THETI(L,NY,NX),FCX,WPX,POROS(L,NY,NX)
!
!     SOIL OSMOTIC, GRAVIMETRIC AND MATRIC WATER POTENTIALS
!
!     PSISM,PSISO,PSISH,PSIST=matric,osmotic,gravimetric,total water potential
!
  PSISO(L,NY,NX)=-8.3143E-06*TKS(L,NY,NX)*CION(L,NY,NX)
  PSISH(L,NY,NX)=0.0098*(ALT(NY,NX)-DPTH(L,NY,NX))
  PSIST(L,NY,NX)=AMIN1(0.0,PSISM(L,NY,NX)+PSISO(L,NY,NX)+PSISH(L,NY,NX))
!     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
!     WRITE(*,1113)'PSISM1',I,J,NX,NY,L,PSISM(L,NY,NX)
!    2,THETW(L,NY,NX),THETW1,VOLX(L,NY,NX)
!    2,FC(L,NY,NX),WP(L,NY,NX),POROS(L,NY,NX),VOLP(L,NY,NX)
!    3,VOLW(L,NY,NX),VOLI(L,NY,NX),VOLA(L,NY,NX),VOLT(L,NY,NX)
!    4,CDPTH(L,NY,NX),DPTH(L,NY,NX),CDPTHZ(L,NY,NX),DPTHZ(L,NY,NX)
!    5,DLYR(3,L,NY,NX),PSIST(L,NY,NX),PSISO(L,NY,NX)
!    2,PSISH(L,NY,NX),TKS(L,NY,NX),CION(L,NY,NX)
!1113  FORMAT(A8,5I4,50E12.4)
!     ENDIF
!
!     SOIL RESISTANCE TO ROOT PENETRATION
!
!     RSCS=soil resistance to root penetration (MPa)
!
!     IF(BKDS(L,NY,NX).GT.ZERO)THEN
!     CCLAYT=CCLAY(L,NY,NX)*1.0E+02
!     CORGCT=CORGC(L,NY,NX)*1.0E-04
!     CC=EXP(-3.6733-0.1447*CCLAYT+0.7653*CORGCT)
!     DD=-0.4805-0.1239*CCLAYT+0.2080*CORGCT
!     EE=3.8521+0.0963*CCLAYT
!     RSCS(L,NY,NX)=CC*THETW(L,NY,NX)**DD*BKDS(L,NY,NX)**EE
!     ELSE
      RSCS(L,NY,NX)=0.0_r8
!     ENDIF
!     WRITE(*,2442)'RSCS',I,J,NX,NY,L,RSCS(L,NY,NX),THETW(L,NY,NX)
!    2,BKDS(L,NY,NX),CCLAY(L,NY,NX),CORGC(L,NY,NX)
!2442  FORMAT(A8,5I4,12E12.4)
!
!     SOIL HYDRAULIC CONDUCTIVITIES FROM AMBIENT SOIL WATER CONTENTS
!
!     CNDU=soil hydraulic conductivity for root uptake
!
  K=MAX(1,MIN(100,INT(100.0*(POROS(L,NY,NX)-THETW(L,NY,NX))/POROS(L,NY,NX))+1))
  CNDU(L,NY,NX)=0.5*(HCND(1,K,L,NY,NX)+HCND(3,K,L,NY,NX))
  end subroutine GetSoilHydraulicVars
!------------------------------------------------------------------------------------------

      subroutine DiagActiveLayerDepth(L,NY,NX)

      implicit none
      integer, intent(in) :: L,NY,NX

      integer :: LL
!     begin_execution
!
!     VOLI,VOLIH=ice volume in micropores,macropores
!     VOLW,VOLWH=water volume in micropores,macropores
!     VOLA,VOLAH=total volume in micropores,macropores
!     DPTHA=active layer depth
!     CDPTH,DLYR=depth to bottom,thickness of soil layer
!
      IF(ICHKA.EQ.0)THEN
      VOLIT=VOLI(L,NY,NX)+VOLIH(L,NY,NX)
      VOLAT=VOLA(L,NY,NX)+VOLAH(L,NY,NX)
      IF(VOLAT.GT.ZEROS2(NY,NX).AND.VOLIT.GT.0.01*VOLAT)THEN
      DO 5700 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
      VOLITL=VOLI(LL,NY,NX)+VOLIH(LL,NY,NX)
      VOLWTL=VOLW(LL,NY,NX)+VOLWH(LL,NY,NX)
      VOLATL=VOLA(LL,NY,NX)+VOLAH(LL,NY,NX)
      IF(VOLATL.GT.ZEROS2(NY,NX).AND.VOLITL.LT.0.01*VOLATL)THEN
      GO TO 5701
      ENDIF
5700  CONTINUE
      IF(VOLAT.GT.ZEROS2(NY,NX))THEN
      DPTHA(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX) &
      *AMIN1(1.0,VOLIT/VOLAT)
      ELSE
      DPTHA(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
      ENDIF
      ICHKA=1
      GO TO 5702
5701  DPTHA(NY,NX)=9999.0
5702  CONTINUE
      ENDIF
      ENDIF
      end subroutine DiagActiveLayerDepth
!------------------------------------------------------------------------------------------

  subroutine SetTracerPropertyInLiterAir(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,L
!     begin_execution
!
!     LITTER GAS CONCENTRATIOS
!
!     C*G=soil gas gaseous concentration
!     *E=atmospheric concentration
!     TKS,TCS=litter temperature (K,C)
!     S*L=gas solubility
!     C*S=soil gas aqueous concentration
!
  CCO2G(0,NY,NX)=CO2E(NY,NX)*5.36E-04*TREF/TKS(0,NY,NX)
  CCH4G(0,NY,NX)=CH4E(NY,NX)*5.36E-04*TREF/TKS(0,NY,NX)
  COXYG(0,NY,NX)=OXYE(NY,NX)*1.43E-03*TREF/TKS(0,NY,NX)
  CZ2GG(0,NY,NX)=Z2GE(NY,NX)*1.25E-03*TREF/TKS(0,NY,NX)
  CZ2OG(0,NY,NX)=Z2OE(NY,NX)*1.25E-03*TREF/TKS(0,NY,NX)
  CNH3G(0,NY,NX)=ZNH3E(NY,NX)*6.25E-04*TREF/TKS(0,NY,NX)
  CH2GG(0,NY,NX)=H2GE(NY,NX)*8.92E-05*TREF/TKS(0,NY,NX)
  CNH4B(0,NY,NX)=0.0_r8
  CNH3B(0,NY,NX)=0.0_r8
  CNO3B(0,NY,NX)=0.0_r8
  CNO2B(0,NY,NX)=0.0_r8
  CH2P4B(0,NY,NX)=0.0_r8
  SCO2L(0,NY,NX)=SCO2X*EXP(0.843-0.0281*TCS(0,NY,NX))
  SCH4L(0,NY,NX)=SCH4X*EXP(0.597-0.0199*TCS(0,NY,NX))
  SOXYL(0,NY,NX)=SOXYX*EXP(0.516-0.0172*TCS(0,NY,NX))
  SN2GL(0,NY,NX)=SN2GX*EXP(0.456-0.0152*TCS(0,NY,NX))
  SN2OL(0,NY,NX)=SN2OX*EXP(0.897-0.0299*TCS(0,NY,NX))
  SNH3L(0,NY,NX)=SNH3X*EXP(0.513-0.0171*TCS(0,NY,NX))
  SH2GL(0,NY,NX)=SH2GX*EXP(0.597-0.0199*TCS(0,NY,NX))
  IF(VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    CCO2S(0,NY,NX)=AMAX1(0.0,CO2S(0,NY,NX)/VOLW(0,NY,NX))
    CCH4S(0,NY,NX)=AMAX1(0.0,CH4S(0,NY,NX)/VOLW(0,NY,NX))
    COXYS(0,NY,NX)=AMAX1(0.0,OXYS(0,NY,NX)/VOLW(0,NY,NX))
    CZ2GS(0,NY,NX)=AMAX1(0.0,Z2GS(0,NY,NX)/VOLW(0,NY,NX))
    CZ2OS(0,NY,NX)=AMAX1(0.0,Z2OS(0,NY,NX)/VOLW(0,NY,NX))
    CH2GS(0,NY,NX)=AMAX1(0.0,H2GS(0,NY,NX)/VOLW(0,NY,NX))
  ELSE
    CCO2S(0,NY,NX)=0.0_r8
    CCH4S(0,NY,NX)=0.0_r8
    COXYS(0,NY,NX)=0.0_r8
    CZ2GS(0,NY,NX)=0.0_r8
    CZ2OS(0,NY,NX)=0.0_r8
    CH2GS(0,NY,NX)=0.0_r8
  ENDIF
!
!     TFACL=temperature effect on diffusivity
!     *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
!     *SG PARAMETER statement above
!
  TFACL=(TKS(0,NY,NX)/298.15)**6
  TFND(0,NY,NX)=TFACL
  CQSGL(0,NY,NX)=CQSG*TFACL
  OLSGL(0,NY,NX)=OLSG*TFACL
  ZLSGL(0,NY,NX)=ZLSG*TFACL
  ZNSGL(0,NY,NX)=ZNSG*TFACL
  HLSGL(0,NY,NX)=HLSG*TFACL
  ZVSGL(0,NY,NX)=ZVSG*TFACL
  ZOSGL(0,NY,NX)=ZOSG*TFACL
  POSGL(0,NY,NX)=POSG*TFACL
  OCSGL(0,NY,NX)=OCSG*TFACL
  ONSGL(0,NY,NX)=ONSG*TFACL
  OPSGL(0,NY,NX)=OPSG*TFACL
  OASGL(0,NY,NX)=OASG*TFACL
!
!     R*Y,R*X=total substrate uptake from previous,current hour
!     used in nitro.f, uptake.f
!
  ROXYY(0,NY,NX)=ROXYX(0,NY,NX)
  RNH4Y(0,NY,NX)=RNH4X(0,NY,NX)
  RNO3Y(0,NY,NX)=RNO3X(0,NY,NX)
  RNO2Y(0,NY,NX)=RNO2X(0,NY,NX)
  RN2OY(0,NY,NX)=RN2OX(0,NY,NX)
  RP14Y(0,NY,NX)=RP14X(0,NY,NX)
  RPO4Y(0,NY,NX)=RPO4X(0,NY,NX)
  ROXYX(0,NY,NX)=0.0_r8
  RNH4X(0,NY,NX)=0.0_r8
  RNO3X(0,NY,NX)=0.0_r8
  RNO2X(0,NY,NX)=0.0_r8
  RN2OX(0,NY,NX)=0.0_r8
  RP14X(0,NY,NX)=0.0_r8
  RPO4X(0,NY,NX)=0.0_r8
  DO 5055 K=0,4
    ROQCY(K,0,NY,NX)=ROQCX(K,0,NY,NX)
    ROQAY(K,0,NY,NX)=ROQAX(K,0,NY,NX)
    ROQCX(K,0,NY,NX)=0.0_r8
    ROQAX(K,0,NY,NX)=0.0_r8
5055  CONTINUE
!
!     WGSGA,WGSGR,WGSGW=vapor diffusivity in air,litter,snowpack
!
  TFACA=(TKA(NY,NX)/298.15)**1.75
  WGSGA(NY,NX)=WGSG*TFACA
  if(TKS(0,NY,NX)<0._r8)then
    write(*,*)'TKS(0,NY,NX)=',TKS(0,NY,NX)
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif
  TFACR=(TKS(0,NY,NX)/298.15)**1.75
  WGSGR(NY,NX)=WGSG*TFACR
  DO 5060 L=1,JS
    TFACW=(TKW(L,NY,NX)/298.15)**1.75
    WGSGW(L,NY,NX)=WGSG*TFACW
5060  CONTINUE
  end subroutine SetTracerPropertyInLiterAir
!------------------------------------------------------------------------------------------

  subroutine GetSoluteConcentrations(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
!     begin_execution
!     CALCULATE SOIL CONCENTRATIONS OF NH4, NH3, NO3, PO4
!     IN BAND AND NON-BAND ZONES
!
!     C*S=solute concentration in non-band
!     CH1P4,CH2P4=HPO4,H2PO4 concentration in non-band
!     Z*P=P ion pair amounts in non-band (see solute.f)
!     VLNH4,VLNO3,VLPO4=fraction of soil volume in NH4,NO3,PO4 non-band
!
  IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    IF(VLNH4(L,NY,NX).GT.ZERO)THEN
      CNH4S(L,NY,NX)=AMAX1(0.0,ZNH4S(L,NY,NX) &
        /(VOLW(L,NY,NX)*VLNH4(L,NY,NX)))
      CNH3S(L,NY,NX)=AMAX1(0.0,ZNH3S(L,NY,NX) &
        /(VOLW(L,NY,NX)*VLNH4(L,NY,NX)))
    ELSE
      CNH4S(L,NY,NX)=0.0_r8
      CNH3S(L,NY,NX)=0.0_r8
    ENDIF
    IF(VLNO3(L,NY,NX).GT.ZERO)THEN
      CNO3S(L,NY,NX)=AMAX1(0.0,ZNO3S(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLNO3(L,NY,NX)))
      CNO2S(L,NY,NX)=AMAX1(0.0,ZNO2S(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLNO3(L,NY,NX)))
      ELSE
      CNO3S(L,NY,NX)=0.0_r8
      CNO2S(L,NY,NX)=0.0_r8
      ENDIF
      IF(VLPO4(L,NY,NX).GT.ZERO)THEN
      CH1P4(L,NY,NX)=AMAX1(0.0,H1PO4(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLPO4(L,NY,NX)))
      CH2P4(L,NY,NX)=AMAX1(0.0,H2PO4(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLPO4(L,NY,NX)))
      CPO4S(L,NY,NX)=AMAX1(0.0,((H0PO4(L,NY,NX)+H3PO4(L,NY,NX) &
      +ZFE1P(L,NY,NX)+ZFE2P(L,NY,NX)+ZCA0P(L,NY,NX) &
      +ZCA1P(L,NY,NX)+ZCA2P(L,NY,NX)+ZMG1P(L,NY,NX))*31.0 &
      +H1PO4(L,NY,NX)+H2PO4(L,NY,NX))/(VOLW(L,NY,NX)*VLPO4(L,NY,NX)))
      ELSE
      CH1P4(L,NY,NX)=0.0_r8
      CH2P4(L,NY,NX)=0.0_r8
      CPO4S(L,NY,NX)=0.0_r8
      ENDIF
!
!     C*B=solute concentration in band
!     CH1PB,CH2PB=HPO4,H2PO4 concentration in band
!     Z*B=P ion pair amounts in band (see solute.f)
!     VLNHB,VLNOB,VLPOB=fraction of soil volume in NH4,NO3,PO4 band
!
      IF(VLNHB(L,NY,NX).GT.ZERO)THEN
      CNH4B(L,NY,NX)=AMAX1(0.0,ZNH4B(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLNHB(L,NY,NX)))
      CNH3B(L,NY,NX)=AMAX1(0.0,ZNH3B(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLNHB(L,NY,NX)))
      ELSE
      CNH4B(L,NY,NX)=0.0_r8
      CNH3B(L,NY,NX)=0.0_r8
      ENDIF
      IF(VLNOB(L,NY,NX).GT.ZERO)THEN
      CNO3B(L,NY,NX)=AMAX1(0.0,ZNO3B(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLNOB(L,NY,NX)))
      CNO2B(L,NY,NX)=AMAX1(0.0,ZNO2B(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLNOB(L,NY,NX)))
      ELSE
      CNO3B(L,NY,NX)=0.0_r8
      CNO2B(L,NY,NX)=0.0_r8
      ENDIF
      IF(VLPOB(L,NY,NX).GT.ZERO)THEN
      CH1P4B(L,NY,NX)=AMAX1(0.0,H1POB(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLPOB(L,NY,NX)))
      CH2P4B(L,NY,NX)=AMAX1(0.0,H2POB(L,NY,NX) &
      /(VOLW(L,NY,NX)*VLPOB(L,NY,NX)))
      CPO4B(L,NY,NX)=AMAX1(0.0,((H0POB(L,NY,NX)+H3POB(L,NY,NX) &
      +ZFE1PB(L,NY,NX)+ZFE2PB(L,NY,NX)+ZCA0PB(L,NY,NX) &
      +ZCA1PB(L,NY,NX)+ZCA2PB(L,NY,NX)+ZMG1PB(L,NY,NX))*31.0 &
      +H1POB(L,NY,NX)+H2POB(L,NY,NX))/(VOLW(L,NY,NX)*VLPOB(L,NY,NX)))
      ELSE
      CH1P4B(L,NY,NX)=0.0_r8
      CH2P4B(L,NY,NX)=0.0_r8
      CPO4B(L,NY,NX)=0.0_r8
      ENDIF
      ELSE
      CNH4S(L,NY,NX)=0.0_r8
      CNH3S(L,NY,NX)=0.0_r8
      CNO3S(L,NY,NX)=0.0_r8
      CNO2S(L,NY,NX)=0.0_r8
      CH1P4(L,NY,NX)=0.0_r8
      CH2P4(L,NY,NX)=0.0_r8
      CPO4S(L,NY,NX)=0.0_r8
      CNH4B(L,NY,NX)=0.0_r8
      CNH3B(L,NY,NX)=0.0_r8
      CNO3B(L,NY,NX)=0.0_r8
      CNO2B(L,NY,NX)=0.0_r8
      CH1P4B(L,NY,NX)=0.0_r8
      CH2P4B(L,NY,NX)=0.0_r8
      CPO4B(L,NY,NX)=0.0_r8
      ENDIF
      end subroutine GetSoluteConcentrations
!------------------------------------------------------------------------------------------

  subroutine PrepVars4PlantMicrobeUptake(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: K
! begin_execution
!
! PREPARE ARRAYS FOR TOTAL O2 UPTAKE AND NH4,NO3.NO2,N2O,HPO4,H2PO4
! UPTAKE IN NON-BAND,BAND AND DOC,DON,DOP,ACETATE UPTAKE
!
! R*Y,R*X=total substrate uptake from previous,current hour
! used in nitro.f, uptake.f
!
  ROXYY(L,NY,NX)=ROXYX(L,NY,NX)
  RNH4Y(L,NY,NX)=RNH4X(L,NY,NX)
  RNO3Y(L,NY,NX)=RNO3X(L,NY,NX)
  RNO2Y(L,NY,NX)=RNO2X(L,NY,NX)
  RN2OY(L,NY,NX)=RN2OX(L,NY,NX)
  RP14Y(L,NY,NX)=RP14X(L,NY,NX)
  RPO4Y(L,NY,NX)=RPO4X(L,NY,NX)
  RNHBY(L,NY,NX)=RNHBX(L,NY,NX)
  RN3BY(L,NY,NX)=RN3BX(L,NY,NX)
  RN2BY(L,NY,NX)=RN2BX(L,NY,NX)
  RP1BY(L,NY,NX)=RP1BX(L,NY,NX)
  RPOBY(L,NY,NX)=RPOBX(L,NY,NX)
  ROXYX(L,NY,NX)=0.0_r8
  RNH4X(L,NY,NX)=0.0_r8
  RNO3X(L,NY,NX)=0.0_r8
  RNO2X(L,NY,NX)=0.0_r8
  RN2OX(L,NY,NX)=0.0_r8
  RP14X(L,NY,NX)=0.0_r8
  RPO4X(L,NY,NX)=0.0_r8
  RNHBX(L,NY,NX)=0.0_r8
  RN3BX(L,NY,NX)=0.0_r8
  RN2BX(L,NY,NX)=0.0_r8
  RP1BX(L,NY,NX)=0.0_r8
  RPOBX(L,NY,NX)=0.0_r8
  DO 5050 K=0,4
    ROQCY(K,L,NY,NX)=ROQCX(K,L,NY,NX)
    ROQAY(K,L,NY,NX)=ROQAX(K,L,NY,NX)
    ROQCX(K,L,NY,NX)=0.0_r8
    ROQAX(K,L,NY,NX)=0.0_r8
5050  CONTINUE
!
! DIFFUSIVITY
!
! TFACG,TFACL=temperature effects on gaseous,aqueous diffusivity
!
! *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
! *SG PARAMETER statement above
!
  if(TKS(L,NY,NX)<0._r8)THEN
    WRITE(*,*)'TKS=',L,TKS(L,NY,NX)
    CALL ENDRUN(TRIM(MOD_FILENAME)//' at line',__LINE__)
  ENDIF
  TFACG=(TKS(L,NY,NX)/298.15)**1.75
  TFACL=(TKS(L,NY,NX)/298.15)**6
  TFND(L,NY,NX)=TFACL
  CGSGL(L,NY,NX)=CGSG*TFACG
  CHSGL(L,NY,NX)=CHSG*TFACG
  OGSGL(L,NY,NX)=OGSG*TFACG
  ZGSGL(L,NY,NX)=ZGSG*TFACG
  Z2SGL(L,NY,NX)=Z2SG*TFACG
  ZHSGL(L,NY,NX)=ZHSG*TFACG
  HGSGL(L,NY,NX)=HGSG*TFACG
  CLSGL(L,NY,NX)=CLSG*TFACL
  CQSGL(L,NY,NX)=CQSG*TFACL
  OLSGL(L,NY,NX)=OLSG*TFACL
  ZLSGL(L,NY,NX)=ZLSG*TFACL
  ZNSGL(L,NY,NX)=ZNSG*TFACL
  HLSGL(L,NY,NX)=HLSG*TFACL
  ZVSGL(L,NY,NX)=ZVSG*TFACL
  ZOSGL(L,NY,NX)=ZOSG*TFACL
  POSGL(L,NY,NX)=POSG*TFACL
  OCSGL(L,NY,NX)=OCSG*TFACL
  ONSGL(L,NY,NX)=ONSG*TFACL
  OPSGL(L,NY,NX)=OPSG*TFACL
  OASGL(L,NY,NX)=OASG*TFACL
  WGSGL(L,NY,NX)=WGSG*TFACG
  IF(ISALTG.NE.0)THEN
    ALSGL(L,NY,NX)=ALSG*TFACL
    FESGL(L,NY,NX)=FESG*TFACL
    HYSGL(L,NY,NX)=HYSG*TFACL
    CASGL(L,NY,NX)=CASG*TFACL
    GMSGL(L,NY,NX)=GMSG*TFACL
    ANSGL(L,NY,NX)=ANSG*TFACL
    AKSGL(L,NY,NX)=AKSG*TFACL
    OHSGL(L,NY,NX)=OHSG*TFACL
    C3SGL(L,NY,NX)=C3SG*TFACL
    HCSGL(L,NY,NX)=HCSG*TFACL
    SOSGL(L,NY,NX)=SOSG*TFACL
    CLSXL(L,NY,NX)=CLSX*TFACL
!
!   TOTAL ION CONCENTRATION
!
!   ZC3,ZA3,ZC2,ZA2,ZC1,ZA1=total tri-,di-,univalent cations C,anions A
!   CSTR,CION=ion strength, total ion concentration
!
    ZC3=ZAL(L,NY,NX)+ZFE(L,NY,NX)
    ZA3=H0PO4(L,NY,NX)+H0POB(L,NY,NX)
    ZC2=ZCA(L,NY,NX)+ZMG(L,NY,NX)+ZALOH1(L,NY,NX)+ZFEOH1(L,NY,NX) &
      +ZFE2P(L,NY,NX)+ZFE2PB(L,NY,NX)
    ZA2=ZSO4(L,NY,NX)+ZCO3(L,NY,NX)+H1PO4(L,NY,NX)+H1POB(L,NY,NX)
    ZC1=(ZNH4S(L,NY,NX)+ZNH4B(L,NY,NX))/14.0+ZHY(L,NY,NX) &
      +ZNA(L,NY,NX)+ZKA(L,NY,NX)+ZALOH2(L,NY,NX)+ZFEOH2(L,NY,NX) &
      +ZALS(L,NY,NX)+ZFES(L,NY,NX)+ZCAO(L,NY,NX)+ZCAH(L,NY,NX) &
      +ZMGO(L,NY,NX)+ZMGH(L,NY,NX)+ZFE1P(L,NY,NX)+ZFE1PB(L,NY,NX) &
      +ZCA2P(L,NY,NX)+ZCA2PB(L,NY,NX)
    ZA1=(ZNO3S(L,NY,NX)+ZNO3B(L,NY,NX))/14.0+ZOH(L,NY,NX) &
      +ZHCO3(L,NY,NX)+ZCL(L,NY,NX)+ZALOH4(L,NY,NX)+ZFEOH4(L,NY,NX) &
      +ZNAC(L,NY,NX)+ZNAS(L,NY,NX)+ZKAS(L,NY,NX)+(H2PO4(L,NY,NX) &
      +H2POB(L,NY,NX))/31.0+ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX)
    ZN=CO2S(L,NY,NX)/12.0+CH4S(L,NY,NX)/12.0+OXYS(L,NY,NX)/32.0 &
      +(Z2GS(L,NY,NX)+Z2OS(L,NY,NX)+ZNH3S(L,NY,NX)+ZNH3B(L,NY,NX))/14.0 &
      +ZALOH3(L,NY,NX)+ZFEOH3(L,NY,NX)+ZCAC(L,NY,NX)+ZCAS(L,NY,NX) &
      +ZMGC(L,NY,NX)+ZMGS(L,NY,NX)+H3PO4(L,NY,NX)+ZCA1P(L,NY,NX) &
      +ZMG1P(L,NY,NX)+H3POB(L,NY,NX)+ZCA1PB(L,NY,NX)+ZMG1PB(L,NY,NX)
    ZION1=ABS(3.0*(ZC3-ZA3)+2.0*(ZC2-ZA2)+ZC1-ZA1)
    IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CSTR(L,NY,NX)=AMAX1(0.0_r8,0.5E-03_r8*(9.0_r8*(ZC3+ZA3)+4.0_r8*(ZC2+ZA2) &
        +ZC1+ZA1+ZION1)/VOLW(L,NY,NX))
      CION(L,NY,NX)=AMAX1(0.0_r8,(ZC3+ZA3+ZC2+ZA2+ZC1+ZA1+ZN)/VOLW(L,NY,NX))
    ELSE
      CSTR(L,NY,NX)=0.0_r8
      CION(L,NY,NX)=0.0_r8
    ENDIF
  ENDIF
! IF(L.EQ.1)THEN
!     WRITE(*,1113)'CION',I,J,NX,NY,L,CION(L,NY,NX)
!    2,CSTR(L,NY,NX),ZC3,ZA3,ZC2,ZA2,ZC1,ZA1,ZN,VOLW(L,NY,NX)
!    3,ZAL(L,NY,NX),ZFE(L,NY,NX),ZNO3S(L,NY,NX),ZNO3B(L,NY,NX)
!    4,ZOH(L,NY,NX),ZHCO3(L,NY,NX),ZCL(L,NY,NX),ZALOH4(L,NY,NX)
!    5,ZFEOH4(L,NY,NX),ZNAC(L,NY,NX),ZNAS(L,NY,NX)
!    6,ZKAS(L,NY,NX),H2PO4(L,NY,NX),H2POB(L,NY,NX)
!    7,ZCA0P(L,NY,NX),ZCA0PB(L,NY,NX),ZCO3(L,NY,NX)
! ENDIF
!
! OSTWALD COEFFICIENTS FOR CO2, CH4, O2, N2, N2O, NH3 AND H2
! SOLUBILITY IN WATER
!
! S*L=solubility of gas in water
! TCS=soil temperature (oC)
!
  FH2O=5.56E+04_r8/(5.56E+04_r8+CION(L,NY,NX))
  SCO2L(L,NY,NX)=SCO2X/(EXP(ACO2X*CSTR(L,NY,NX))) &
      *EXP(0.843_r8-0.0281_r8*TCS(L,NY,NX))*FH2O
  SCH4L(L,NY,NX)=SCH4X/(EXP(ACH4X*CSTR(L,NY,NX))) &
      *EXP(0.597_r8-0.0199_r8*TCS(L,NY,NX))*FH2O
  SOXYL(L,NY,NX)=SOXYX/(EXP(AOXYX*CSTR(L,NY,NX))) &
      *EXP(0.516_r8-0.0172_r8*TCS(L,NY,NX))*FH2O
  SN2GL(L,NY,NX)=SN2GX/(EXP(AN2GX*CSTR(L,NY,NX))) &
      *EXP(0.456_r8-0.0152_r8*TCS(L,NY,NX))*FH2O
  SN2OL(L,NY,NX)=SN2OX/(EXP(AN2OX*CSTR(L,NY,NX))) &
      *EXP(0.897_r8-0.0299_r8*TCS(L,NY,NX))*FH2O
  SNH3L(L,NY,NX)=SNH3X/(EXP(ANH3X*CSTR(L,NY,NX))) &
      *EXP(0.513_r8-0.0171_r8*TCS(L,NY,NX))*FH2O
  SH2GL(L,NY,NX)=SH2GX/(EXP(AH2GX*CSTR(L,NY,NX))) &
      *EXP(0.597_r8-0.0199_r8*TCS(L,NY,NX))*FH2O
  end subroutine PrepVars4PlantMicrobeUptake
!------------------------------------------------------------------------------------------

  subroutine GetSurfResidualProperties(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
! begin_execution
! PHYSICAL PROPERTIES, AND WATER, GAS, AND MINERAL CONTENTS
! OF SURFACE RESIDUE
!
! VOLWRX=liter water holding capacity
! THETRX,RC0=specific WHC,mass of woody(0),fine(1),manure(2) litter
! VOLR=dry litter volume
! BKRS=dry bulk density of woody(0),fine(1),manure(2) litter
! TVOLG0=excess litter water+ice
! VOLT,VOLX=wet litter volume
! BKVL=litter mass
! VOLW,VOLI,VOLA,VOLP=litter water,ice,porosity,air volume
! THETW,THETI,THETA,THETP=litter water,ice,porosity,air concentration
! POROS=litter porosity
! THETW0,THETI0,DPTH0=litter excess water,ice,water+ice depth
! DLYR=litter thickness
! PSISM,PSISE=litter matric,saturation water potential
!
  TVOLG0=AMAX1(0.0,VOLW(0,NY,NX)+VOLI(0,NY,NX)-VOLWRX(NY,NX))
  VOLT(0,NY,NX)=TVOLG0+VOLR(NY,NX)
  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLX(0,NY,NX)=VOLT(0,NY,NX)
    BKVL(0,NY,NX)=1.82E-06_r8*ORGC(0,NY,NX)
    VOLA(0,NY,NX)=AMAX1(0.0_r8,VOLR(NY,NX)-BKVL(0,NY,NX)/1.30_r8)
    VOLP(0,NY,NX)=AMAX1(0.0_r8,VOLA(0,NY,NX)-VOLW(0,NY,NX)-VOLI(0,NY,NX))
    IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
      POROS(0,NY,NX)=VOLA(0,NY,NX)/VOLR(NY,NX)
      THETW(0,NY,NX)=AMAX1(0.0_r8,AMIN1(1.0_r8,VOLW(0,NY,NX)/VOLR(NY,NX)))
      THETI(0,NY,NX)=AMAX1(0.0_r8,AMIN1(1.0_r8,VOLI(0,NY,NX)/VOLR(NY,NX)))
      THETP(0,NY,NX)=AMAX1(0.0_r8,AMIN1(1.0_r8,VOLP(0,NY,NX)/VOLR(NY,NX)))
    ELSE
      POROS(0,NY,NX)=1.0_r8
      THETW(0,NY,NX)=0.0_r8
      THETI(0,NY,NX)=0.0_r8
      THETP(0,NY,NX)=0.0_r8
    ENDIF
    TVOLWI=VOLW(0,NY,NX)+VOLI(0,NY,NX)
    IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VOLWRZ=VOLW(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      VOLIRZ=VOLI(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      XVOLW0=AMAX1(0.0_r8,VOLW(0,NY,NX)-VOLWRZ)/AREA(3,NU(NY,NX),NY,NX)
      XVOLI0=AMAX1(0.0_r8,VOLI(0,NY,NX)-VOLIRZ)/AREA(3,NU(NY,NX),NY,NX)
    ELSE
      XVOLW0=0.0_r8
      XVOLI0=0.0_r8
    ENDIF
    DPTH0(NY,NX)=XVOLW0+XVOLI0
    DLYR(3,0,NY,NX)=VOLX(0,NY,NX)/AREA(3,0,NY,NX)
    IF(VOLR(NY,NX).GT.ZEROS(NY,NX).AND.VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
        PSISM(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)+((FCL(0,NY,NX)-LOG(THETWR)) &
          /FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS(0,NY,NX))THEN
        PSISM(0,NY,NX)=-EXP(PSIMS(NY,NX) &
          +(((PSL(0,NY,NX)-LOG(THETWR)) &
          /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
        PSISM(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      PSISO(0,NY,NX)=0.0
      PSISH(0,NY,NX)=0.0098*(ALT(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX) &
        +0.5*DLYR(3,0,NY,NX))
      PSIST(0,NY,NX)=AMIN1(0.0,PSISM(0,NY,NX)+PSISO(0,NY,NX) &
        +PSISH(0,NY,NX))
!
!     LITTER NH4,NH3,NO3,NO2,HPO4,H2PO4 CONCENTRATIONS
!
!     C*=litter solute concentrations
!
      CNH4S(0,NY,NX)=AMAX1(0.0,ZNH4S(0,NY,NX)/VOLW(0,NY,NX))
      CNH3S(0,NY,NX)=AMAX1(0.0,ZNH3S(0,NY,NX)/VOLW(0,NY,NX))
      CNO3S(0,NY,NX)=AMAX1(0.0,ZNO3S(0,NY,NX)/VOLW(0,NY,NX))
      CNO2S(0,NY,NX)=AMAX1(0.0,ZNO2S(0,NY,NX)/VOLW(0,NY,NX))
      CH1P4(0,NY,NX)=AMAX1(0.0,H1PO4(0,NY,NX)/VOLW(0,NY,NX))
      CH2P4(0,NY,NX)=AMAX1(0.0,H2PO4(0,NY,NX)/VOLW(0,NY,NX))
    ELSE
      PSISM(0,NY,NX)=PSISM(NU(NY,NX),NY,NX)
      CNH4S(0,NY,NX)=0.0_r8
      CNH3S(0,NY,NX)=0.0_r8
      CNO3S(0,NY,NX)=0.0_r8
      CNO2S(0,NY,NX)=0.0_r8
      CH1P4(0,NY,NX)=0.0_r8
      CH2P4(0,NY,NX)=0.0_r8
      ENDIF
      ELSE
      VOLX(0,NY,NX)=0.0
      BKVL(0,NY,NX)=0.0
      VOLA(0,NY,NX)=0.0
      VOLP(0,NY,NX)=0.0
      POROS(0,NY,NX)=1.0
      DLYR(3,0,NY,NX)=0.0
      THETW(0,NY,NX)=0.0
      THETI(0,NY,NX)=0.0
      THETP(0,NY,NX)=1.0
      VOLWRX(NY,NX)=0.0
      PSISM(0,NY,NX)=PSISM(NU(NY,NX),NY,NX)
      CNH4S(0,NY,NX)=0.0
      CNH3S(0,NY,NX)=0.0
      CNO3S(0,NY,NX)=0.0
      CNO2S(0,NY,NX)=0.0
      CH1P4(0,NY,NX)=0.0
      CH2P4(0,NY,NX)=0.0
      CCO2S(0,NY,NX)=0.0
      CCH4S(0,NY,NX)=0.0
      COXYS(0,NY,NX)=0.0
      CZ2GS(0,NY,NX)=0.0
      CZ2OS(0,NY,NX)=0.0
      CH2GS(0,NY,NX)=0.0
      ENDIF
      end subroutine GetSurfResidualProperties
!------------------------------------------------------------------------------------------

  subroutine MultiLayerSurfaceRadiation(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: NB,NZ,L,K,M,N,NN
  !     begin_execution
  !     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
  !     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
  !     NODE K, BRANCH NB, PFT NZ
  !
  !     ARLFS,ARLSS=leaf+stalk area of combined,each PFT canopy
  !     ZL=height to bottom of canopy layer
  !     DPTHS,DPTH0=snowpack,surface water depths
  !     ARLFL,ARSTK=leaf,stalk areas of PFT
  !     RAD,RAP=vertical direct+diffuse SW,PAR
  !     RADS,RADY,RAPS,RAPY=solar beam direct,diffuse SW,PAR
  !     SSIN,TYSIN=sine of solar,sky angles
  !     RADC,RADP=total SW,PAR absorbed by canopy
  !     CFX=clumping factor for self-shading
  !
  ARLSS(NY,NX)=0.0
  DO 1135 NZ=1,NP(NY,NX)
    ARLFS(NZ,NY,NX)=0.0
    DO  NB=1,NBR(NZ,NY,NX)
      DO  L=1,JC
        IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO &
          .AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
          DO 1130 K=1,25
            ARLFS(NZ,NY,NX)=ARLFS(NZ,NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
            ARLSS(NY,NX)=ARLSS(NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
1130      CONTINUE
          ARLFS(NZ,NY,NX)=ARLFS(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
          ARLSS(NY,NX)=ARLSS(NY,NX)+ARSTK(L,NB,NZ,NY,NX)
        ENDIF
      enddo
    enddo
1135  CONTINUE
  IF(SSIN(NY,NX).GT.ZERO)THEN
    RAD(NY,NX)=RADS(NY,NX)*SSIN(NY,NX)+RADY(NY,NX)*TYSIN
    RAP(NY,NX)=RAPS(NY,NX)*SSIN(NY,NX)+RAPY(NY,NX)*TYSIN
  ELSE
    RADS(NY,NX)=0.0
    RADY(NY,NX)=0.0
    RAPS(NY,NX)=0.0
    RAPY(NY,NX)=0.0
    RAD(NY,NX)=0.0
    RAP(NY,NX)=0.0
  ENDIF
  TRADC(NY,NX)=0.0
  TRAPC(NY,NX)=0.0
  DO 1025 NZ=1,NP(NY,NX)
    RADC(NZ,NY,NX)=0.0
    RADP(NZ,NY,NX)=0.0
    CFX(NZ,NY,NX)=CF(NZ,NY,NX)*(1.0-0.025 &
      *ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX))
1025  CONTINUE
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     SAZI,SCOS=solar azimuth,cosine of solar angle
  !     BETAG=incident solar angle at ground surface
  !     GCOS,GSIN=cos,sin of ground surface
  !     ZNOON=hour of solar noon from weather file
  !
  IF(SSIN(NY,NX).GT.ZERO)THEN
    SAZI=0.2618*(ZNOON(NY,NX)-J)+4.7124
    SCOS=SQRT(1.0-SSIN(NY,NX)**2)
    DGAZI=COS(GAZI(NY,NX)-SAZI)
    BETAG=AMAX1(0.0,AMIN1(1.0,GCOS(NY,NX)*SSIN(NY,NX) &
      +GSIN(NY,NX)*SCOS*DGAZI))
    IF(ARLSS(NY,NX).GT.0.0)THEN
      SAGL=ASIN(SSIN(NY,NX))
      !
      !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
      !
      !     RADSA,RADWA,RAPSA,RAPWA=SW,PAR absorbed at leaf,stalk surface
      !     perpendicular to incoming radiation
      !
      DO 1050 NZ=1,NP(NY,NX)
        RADSA(NZ,NY,NX)=RADS(NY,NX)*ABSR(NZ,NY,NX)
        RADWA(NZ,NY,NX)=RADS(NY,NX)*ABSRW
        RAPSA(NZ,NY,NX)=RAPS(NY,NX)*ABSP(NZ,NY,NX)
        RAPWA(NZ,NY,NX)=RAPS(NY,NX)*ABSPW
1050  CONTINUE
      !
      !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
      !
      !     ZAZI=leaf azimuth
      !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
      !     ZAGL=determines forward vs backscattering
      !     IALBS=flag for forward vs backscattering
      !
      DO 1100 M=1,4
        ZAZI=SAZI+(M-0.5)*PICON/4
        DAZI=COS(ZAZI-SAZI)
        DO  N=1,4
          BETY=ZCOS(N)*SSIN(NY,NX)+ZSIN(N)*SCOS*DAZI
          BETA(N,M)=ABS(BETY)
          BETX(N,M)=BETA(N,M)/SSIN(NY,NX)
          IF(ZCOS(N).GT.SSIN(NY,NX))THEN
            BETZ=ACOS(BETY)
          ELSE
            BETZ=-ACOS(BETY)
          ENDIF
          IF(BETZ.GT.-PICON2)THEN
            ZAGL=SAGL+2.0*BETZ
          ELSE
            ZAGL=SAGL-2.0*(PICON+BETZ)
          ENDIF
          IF(ZAGL.GT.0.0.AND.ZAGL.LT.PICON)THEN
            IALBS(N,M)=1
          ELSE
            IALBS(N,M)=2
          ENDIF
!
          !     INTENSITY OF ABSORBED DIRECT RADIATION AT LEAF SURFACES
          !
          !     RDNDIR,RDNDIW,PARDIR,PARDIW=SW,PAR flux absorbed by leaf,stalk surfaces
          !     PAR,PARDIF=direct,diffuse PAR flux
          !     RADYL,RAPYL=solar beam diffuse SW,PAR flux
          !     RAFYL,RAFPL=forward scattered diffuse SW,PAR flux
          !     TAUS,TAUY=fraction of direct,diffuse radiation transmitted
!
          DO  NZ=1,NP(NY,NX)
            RDNDIR(N,M,NZ,NY,NX)=RADSA(NZ,NY,NX)*ABS(BETA(N,M))
            RDNDIW(N,M,NZ,NY,NX)=RADWA(NZ,NY,NX)*ABS(BETA(N,M))
            PARDIR(N,M,NZ,NY,NX)=RAPSA(NZ,NY,NX)*ABS(BETA(N,M))
            PARDIW(N,M,NZ,NY,NX)=RAPWA(NZ,NY,NX)*ABS(BETA(N,M))
            DO L=1,JC
              PARDIF(N,M,L,NZ,NY,NX)=0.0
              PAR(N,M,L,NZ,NY,NX)=PARDIR(N,M,NZ,NY,NX)
            enddo
          enddo
        enddo
1100  CONTINUE
      XAREA=1.00/AREA(3,NU(NY,NX),NY,NX)
      YAREA=0.25/AREA(3,NU(NY,NX),NY,NX)
      RADYL=RADY(NY,NX)
      RAPYL=RAPY(NY,NX)
      TAUS(JC+1,NY,NX)=1.0
      TAUY(JC+1)=1.0
      RAFSL(JC+1)=0.0
      RAFPL(JC+1)=0.0
      STOPS=0.0
      !
      !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
      !     LAYERS AND ANGLE CLASSES
      !
      !     TSURF,TSURFB,SURF,SURFB=leaf,stalk total,PFT surface area
!
      DO 1150 NZ=1,NP(NY,NX)
        DO  L=1,JC
          DO  N=1,4
            TSURF(N,L,NZ,NY,NX)=0.0
            TSURFB(N,L,NZ,NY,NX)=0.0
          enddo
        enddo
1150  CONTINUE
      DO 1200 NZ=1,NP(NY,NX)
        DO  NB=1,NBR(NZ,NY,NX)
          DO  L=1,JC
            IF(ZL(L-1,NY,NX).GT.DPTHS(NY,NX)-ZERO &
             .AND.ZL(L-1,NY,NX).GT.DPTH0(NY,NX)-ZERO)THEN
              DO 1205 N=1,4
                DO 1210 K=1,25
                  TSURF(N,L,NZ,NY,NX)=TSURF(N,L,NZ,NY,NX)+SURF(N,L,K,NB,NZ,NY,NX)
1210            CONTINUE
                TSURFB(N,L,NZ,NY,NX)=TSURFB(N,L,NZ,NY,NX) &
                  +SURFB(N,L,NB,NZ,NY,NX)
1205          CONTINUE
            ENDIF
          enddo
        enddo
1200  CONTINUE
      !
      !     CALCULATE ABSORPTION, REFLECTION AND TRANSMISSION OF DIRECT AND
      !     DIFFUSE DOWNWARD TOTAL AND VISIBLE RADIATION BY EACH SPECIES
      !     NZ IN EACH LAYER L
      !
      !     RAFYL,RAFPL=forward scattered diffuse SW,PAR
      !     RABYL,RABPL=backscattered diffuse SW,PAR
      !     RADYL,RAPYL=solar beam diffuse SW,PAR
      !     STOPY,STOPSZ,STOPYZ=fraction of direct,diffuse radiation intercepted
      !
      DO 1800 L=JC,1,-1
        IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO &
          .AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
          RADYL=RADYL*TAUY(L+1)+RAFSL(L+1)
          RAPYL=RAPYL*TAUY(L+1)+RAFPL(L+1)
          RAFSL(L)=0.0
          RAFPL(L)=0.0
          RABSL(L)=0.0
          RABPL(L)=0.0
          STOPY=0.0
          STOPSZ=0.0
          STOPYZ=0.0
    !
          !     RESET ACCUMULATORS OB ABSORBED, REFLECTED AND TRANSMITTED RADIATION
          !
          !     RADSL,RADSW,RADPL,RADPW=direct atmosph SW,PAR absbd by leaf,stalk surf
          !     RAYSL,RAYSW,RAYPL,RAYPW=diffuse atmosph SW,PAR absbd by leaf,stalk surf
          !     RADS1,RADW1,RADP1,RADQ1=backscattered direct SW,PAR absbd by leaf,stalk surf
          !     RAYS1,RAYW1,RAYP1,RAYQ1=backscattered diffuse SW,PAR absbd by leaf,stalk surf
          !     RADS2,RADW2,RADP2,RADQ2=fwd scattered direct SW,PAR absbd by leaf,stalk surf
          !     RAYS2,RAYW2,RAYP2,RAYQ2=fwd scattered diffuse SW,PAR absbd by leaf,stalk surf
    !
          DO 1500 NZ=1,NP(NY,NX)
            RADSL(NZ,NY,NX)=0.0
            RADSW(NZ,NY,NX)=0.0
            RADPL(NZ,NY,NX)=0.0
            RADPW(NZ,NY,NX)=0.0
            RAYSL(NZ,NY,NX)=0.0
            RAYSW(NZ,NY,NX)=0.0
            RAYPL(NZ,NY,NX)=0.0
            RAYPW(NZ,NY,NX)=0.0
            RADS1(NZ,NY,NX)=0.0
            RADW1(NZ,NY,NX)=0.0
            RADP1(NZ,NY,NX)=0.0
            RADQ1(NZ,NY,NX)=0.0
            RAYS1(NZ,NY,NX)=0.0
            RAYW1(NZ,NY,NX)=0.0
            RAYP1(NZ,NY,NX)=0.0
            RAYQ1(NZ,NY,NX)=0.0
            RADS2(NZ,NY,NX)=0.0
            RADW2(NZ,NY,NX)=0.0
            RADP2(NZ,NY,NX)=0.0
            RADQ2(NZ,NY,NX)=0.0
            RAYS2(NZ,NY,NX)=0.0
            RAYW2(NZ,NY,NX)=0.0
            RAYP2(NZ,NY,NX)=0.0
            RAYQ2(NZ,NY,NX)=0.0
    !
      !     LEAF SURFACE AREA IN EACH INCLINATION CLASS N, AZIMUTH CLASS M,
      !     LAYER L AND SPECIES NZ
      !
      !     TSURFY=unself-shaded leaf area
      !     TSURFZ=unself-shaded leaf area m-2 in each azimuth class
      !     TSURFS=TSURFY with shading from canopy layers above
      !     TSURFX=TSURFS m-2
      !     TSURWY=unself-shaded stalk area
      !     TSURWZ=unself-shaded stalk area m-2 in each azimuth class
      !     TSURWS=TSURWY with shading from canopy layers above
      !     TSURWX=TSURWS m-2
      !
            DO 1600 N=1,4
              TSURFY=TSURF(N,L,NZ,NY,NX)*CFX(NZ,NY,NX)
              TSURFZ=TSURFY*YAREA
              TSURFS=TSURFY*TAUS(L+1,NY,NX)
              TSURFX=TSURFS*XAREA
              TSURWY=TSURFB(N,L,NZ,NY,NX)*CFW
              TSURWZ=TSURWY*YAREA
              TSURWS=TSURWY*TAUS(L+1,NY,NX)
              TSURWX=TSURWS*XAREA
              !
              !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
              !
              !     STOPZ=accumulated horizontal area of intercepted direct radiation
              !
              DO 1700 M=1,4
                RADSL(NZ,NY,NX)=RADSL(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
                RADSW(NZ,NY,NX)=RADSW(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
                RADPL(NZ,NY,NX)=RADPL(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
                RADPW(NZ,NY,NX)=RADPW(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
                STOPSZ=STOPSZ+(TSURFX+TSURWX)*BETX(N,M)
!
          !     BACKSCATTERING OF REFLECTED DIRECT RADIATION
          !
                IF(IALBS(N,M).EQ.1)THEN
                  RADS1(NZ,NY,NX)=RADS1(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
                  RADW1(NZ,NY,NX)=RADW1(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
                  RADP1(NZ,NY,NX)=RADP1(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
                  RADQ1(NZ,NY,NX)=RADQ1(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
                  !
                  ! FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
                  !
                ELSE
                  RADS2(NZ,NY,NX)=RADS2(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
                  RADW2(NZ,NY,NX)=RADW2(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
                  RADP2(NZ,NY,NX)=RADP2(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
                  RADQ2(NZ,NY,NX)=RADQ2(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
                ENDIF
    !
                !     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
                !
                !     RADYN,RADYW,RAPYN,RAPYW=diffuse SW,PAR flux absorbed by leaf,stalk surf
                !     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface
!
                DO 1750 NN=1,4
                  RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ,NY,NX)
                  RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ,NY,NX)
                  RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
                  PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)+RAPYN
                  PAR(N,M,L,NZ,NY,NX)=PAR(N,M,L,NZ,NY,NX)+RAPYN
!
                  !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
                  !
                  !     STOPYZ=accumulated horizontal area of intercepted diffuse radiation
                  !
                  RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)+TSURFY*RADYN
                  RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)+TSURWY*RADYW
                  RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)+TSURFY*RAPYN
                  RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)+TSURWY*RAPYW
                  STOPYZ=STOPYZ+(TSURFZ+TSURWZ)*OMEGX(M,N,NN)
    !
                  !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
                  !
                  IF(IALBY(M,N,NN).EQ.1)THEN
                    RAYS1(NZ,NY,NX)=RAYS1(NZ,NY,NX)+TSURFY*RADYN
                    RAYW1(NZ,NY,NX)=RAYW1(NZ,NY,NX)+TSURWY*RADYW
                    RAYP1(NZ,NY,NX)=RAYP1(NZ,NY,NX)+TSURFY*RAPYN
                    RAYQ1(NZ,NY,NX)=RAYQ1(NZ,NY,NX)+TSURWY*RAPYW
                    !
                    !     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
                    !
                  ELSE
                    RAYS2(NZ,NY,NX)=RAYS2(NZ,NY,NX)+TSURFY*RADYN
                    RAYW2(NZ,NY,NX)=RAYW2(NZ,NY,NX)+TSURWY*RADYW
                    RAYP2(NZ,NY,NX)=RAYP2(NZ,NY,NX)+TSURFY*RAPYN
                    RAYQ2(NZ,NY,NX)=RAYQ2(NZ,NY,NX)+TSURWY*RAPYW
                  ENDIF
1750            CONTINUE
1700          CONTINUE
1600        CONTINUE
1500      CONTINUE
          !
          !     ACCUMULATED INTERCEPTION BY CANOPY LAYER
          !
          !     XTAUS=interception of direct radiation in current layer
          !     STOPZ=accumulated interception of direct radiation from topmost layer
          !     TAUS=transmission of direct radiation to next lower layer
!
          IF(STOPS+STOPSZ.GT.1.0)THEN
            IF(STOPSZ.GT.ZERO)THEN
              XTAUS=(1.0-STOPS)/((1.0-STOPS)-(1.0-STOPS-STOPSZ))
            ELSE
              XTAUS=0.0
            ENDIF
            TAUS(L+1,NY,NX)=TAUS(L+1,NY,NX)*XTAUS
            STOPSZ=STOPSZ*XTAUS
            DO 1510 NZ=1,NP(NY,NX)
              RADSL(NZ,NY,NX)=RADSL(NZ,NY,NX)*XTAUS
              RADSW(NZ,NY,NX)=RADSW(NZ,NY,NX)*XTAUS
              RADPL(NZ,NY,NX)=RADPL(NZ,NY,NX)*XTAUS
              RADPW(NZ,NY,NX)=RADPW(NZ,NY,NX)*XTAUS
              RADS1(NZ,NY,NX)=RADS1(NZ,NY,NX)*XTAUS
              RADW1(NZ,NY,NX)=RADW1(NZ,NY,NX)*XTAUS
              RADP1(NZ,NY,NX)=RADP1(NZ,NY,NX)*XTAUS
              RADQ1(NZ,NY,NX)=RADQ1(NZ,NY,NX)*XTAUS
              RADS2(NZ,NY,NX)=RADS2(NZ,NY,NX)*XTAUS
              RADW2(NZ,NY,NX)=RADW2(NZ,NY,NX)*XTAUS
              RADP2(NZ,NY,NX)=RADP2(NZ,NY,NX)*XTAUS
              RADQ2(NZ,NY,NX)=RADQ2(NZ,NY,NX)*XTAUS
1510        CONTINUE
          ENDIF
!
          !     XTAUY=interception of diffuse radiation in current layer
          !     STOPYZ=accumulated interception of diffuse radiation from topmost layer
          !     TAUY=transmission of diffuse radiation to next lower layer
!
          IF(STOPY+STOPYZ.GT.1.0)THEN
            XTAUY=(1.0-STOPY)/((1.0-STOPY)-(1.0-STOPY-STOPYZ))
            TAUY(L+1)=TAUY(L+1)*XTAUY
            STOPYZ=STOPYZ*XTAUY
            DO 1520 NZ=1,NP(NY,NX)
              RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)*XTAUY
              RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)*XTAUY
              RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)*XTAUY
              RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)*XTAUY
              RAYS1(NZ,NY,NX)=RAYS1(NZ,NY,NX)*XTAUY
              RAYW1(NZ,NY,NX)=RAYW1(NZ,NY,NX)*XTAUY
              RAYP1(NZ,NY,NX)=RAYP1(NZ,NY,NX)*XTAUY
              RAYQ1(NZ,NY,NX)=RAYQ1(NZ,NY,NX)*XTAUY
              RAYS2(NZ,NY,NX)=RAYS2(NZ,NY,NX)*XTAUY
              RAYW2(NZ,NY,NX)=RAYW2(NZ,NY,NX)*XTAUY
              RAYP2(NZ,NY,NX)=RAYP2(NZ,NY,NX)*XTAUY
              RAYQ2(NZ,NY,NX)=RAYQ2(NZ,NY,NX)*XTAUY
              DO 1730 N=1,4
                DO  M=1,4
                  PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)*XTAUY
                  PAR(N,M,L,NZ,NY,NX)=PARDIR(N,M,NZ,NY,NX)+PARDIF(N,M,L,NZ,NY,NX)
                enddo
1730          CONTINUE
1520        CONTINUE
          ENDIF
          !
          !     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
          !
          !     RADST,RADWT,RADPT,RADQT=total atmospheric SW,PAR absorbed by leaf,stalk
          !     RA1ST,RA1WT,RA1PT,RA1QT=total backscattered SW,PAR absd by leaf,stalk
          !     RA2ST,RA2WT,RA2PT,RA2QT=total fwd scattered SW,PAR absd by leaf,stalk
          !     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
          !     RABSL,RABPL=total back scattered SW,PAR to next layer
          !     RADC,TRADC,RADP,TRADP=total atmospheric SW,PAR absbd by each,all PFT
          !     STOPS,STOPY=accumulated interception of direct,diffuse radiation
          !     TAUS,TAUY=transmission of direct,diffuse radiation to next lower layer
          !
          DO 1530 NZ=1,NP(NY,NX)
            RADST=RADSL(NZ,NY,NX)+RAYSL(NZ,NY,NX)
            RADWT=RADSW(NZ,NY,NX)+RAYSW(NZ,NY,NX)
            RADPT=RADPL(NZ,NY,NX)+RAYPL(NZ,NY,NX)
            RADQT=RADPW(NZ,NY,NX)+RAYPW(NZ,NY,NX)
            RA1ST=RADS1(NZ,NY,NX)+RAYS1(NZ,NY,NX)
            RA1WT=RADW1(NZ,NY,NX)+RAYW1(NZ,NY,NX)
            RA1PT=RADP1(NZ,NY,NX)+RAYP1(NZ,NY,NX)
            RA1QT=RADP1(NZ,NY,NX)+RAYQ1(NZ,NY,NX)
            RA2ST=RADS2(NZ,NY,NX)+RAYS2(NZ,NY,NX)
            RA2WT=RADW2(NZ,NY,NX)+RAYW2(NZ,NY,NX)
            RA2PT=RADP2(NZ,NY,NX)+RAYP2(NZ,NY,NX)
            RA2QT=RADQ2(NZ,NY,NX)+RAYQ2(NZ,NY,NX)
            RAFSL(L)=RAFSL(L)+(RADST*TAUR(NZ,NY,NX) &
              +RA2ST*ALBR(NZ,NY,NX)+RA2WT*ALBRW)*YAREA
            RAFPL(L)=RAFPL(L)+(RADPT*TAUP(NZ,NY,NX) &
              +RA2PT*ALBP(NZ,NY,NX)+RA2QT*ALBPW)*YAREA
            RABSL(L)=RABSL(L)+(RA1ST*ALBR(NZ,NY,NX)+RA1WT*ALBRW)*YAREA
            RABPL(L)=RABPL(L)+(RA1PT*ALBP(NZ,NY,NX)+RA1QT*ALBPW)*YAREA
            RADC(NZ,NY,NX)=RADC(NZ,NY,NX)+RADST+RADWT
            RADP(NZ,NY,NX)=RADP(NZ,NY,NX)+RADPT+RADQT
            TRADC(NY,NX)=TRADC(NY,NX)+RADST+RADWT
            TRAPC(NY,NX)=TRAPC(NY,NX)+RADPT+RADQT
1530      CONTINUE
          STOPS=STOPS+STOPSZ
          STOPY=STOPY+STOPYZ
          TAUS(L,NY,NX)=1.0-STOPS
          TAU0(L,NY,NX)=1.0-TAUS(L,NY,NX)
          TAUY(L)=1.0-STOPY
        ELSE
          RAFSL(L)=RAFSL(L+1)
          RAFPL(L)=RAFPL(L+1)
          TAUS(L,NY,NX)=TAUS(L+1,NY,NX)
          TAU0(L,NY,NX)=1.0-TAUS(L,NY,NX)
          TAUY(L)=TAUY(L+1)
        ENDIF
1800  CONTINUE
      !
      !     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
      !
      !     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horizl ground surface
      !     RADS,RAPS =solar beam direct SW,PAR flux
      !     TAUS,TAUY=transmission of direct,diffuse radiation below canopy
      !     RADYL,RAPYL=solar beam diffuse SW,PAR flux
      !     RASG,RAPG=total SW,PAR at ground surface
      !     BETAG,OMEGAG=incident solar,sky angle at ground surface
!
      RADSG=RADS(NY,NX)*TAUS(1,NY,NX)
      RADYG=RADYL*TAUY(1)+RAFSL(1)
      RAPSG=RAPS(NY,NX)*TAUS(1,NY,NX)
      RAPYG=RAPYL*TAUY(1)+RAFPL(1)
      RASG=ABS(BETAG)*RADSG
      RAPG=ABS(BETAG)*RAPSG
      DO 20 N=1,4
        RASG=RASG+ABS(OMEGAG(N,NY,NX))*RADYG
        RAPG=RAPG+ABS(OMEGAG(N,NY,NX))*RAPYG
20    CONTINUE
      RADG(NY,NX)=RASG*AREA(3,NU(NY,NX),NY,NX)
!
      !     RADIATION REFLECTED FROM GROUND SURFACE
      !
      !     VHCPW,VHCPWX=current,minimum snowpack heat capacity
      !     ALBW,VOLSS,VOLWS,VOLIS=snowpack albedo,snow,water,ice volume
      !     ALBG,ALBS,FSNOW=ground,soil albedo,snow cover fraction
      !     THETW1=soil surface water content
      !     RABSL,RADPL=SW,PAR backscatter from ground surface
      !     TRADG,TRAPG=SW,PAR absorbed by ground surface
!
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
        ALBW=(0.80*VOLSS(NY,NX)+0.30*VOLIS(NY,NX)+0.06*VOLWS(NY,NX)) &
          /(VOLSS(NY,NX)+VOLIS(NY,NX)+VOLWS(NY,NX))
        FSNOW=AMIN1((DPTHS(NY,NX)/0.07)**2,1.0)
        ALBG=FSNOW*ALBW+(1.0-FSNOW)*ALBS(NY,NX)
      ELSE
        IF(VOLX(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
          THETW1=AMIN1(POROS(NU(NY,NX),NY,NX) &
            ,VOLW(NU(NY,NX),NY,NX)/VOLY(NU(NY,NX),NY,NX))
        ELSE
          THETW1=0.0
        ENDIF
        ALBG=AMIN1(ALBX(NY,NX),ALBS(NY,NX) &
          +AMAX1(0.0,ALBX(NY,NX)-THETW1))
      ENDIF
      RABSL(0)=RASG*ALBG*0.25
      RABPL(0)=RAPG*ALBG*0.25
      TRADG(NY,NX)=(1.0-ALBG)*RASG*AREA(3,NU(NY,NX),NY,NX)
      TRAPG(NY,NX)=(1.0-ALBG)*RAPG*AREA(3,NU(NY,NX),NY,NX)
!
      !     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
      !
      !     RABSL,RABPL=total backscattered SW,PAR to next layer
      !     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
      !     RADYN,RADYW,RAPYN,RAPYW=leaf,stalk SW,PAR absbd fwd+back flux
      !     RAYSL,RAYSW,RAYPL,RAYPW=total leaf,stalk SW,PAR absbd fwd+back
      !     RADC,TRADC,RADP,TRADP=total SW,PAR absbd by each,all PFT
!
      RADYL=0.0
      RAPYL=0.0
      TAUY(0)=1.0
      RAFSL(0)=0.0
      RAFPL(0)=0.0
      DO 2800 L=1,JC
        IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO &
         .AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
          RADYL=RADYL*TAUY(L-1)+RAFSL(L-1)+RABSL(L-1)
          RAPYL=RAPYL*TAUY(L-1)+RAFPL(L-1)+RABPL(L-1)
          RAFSL(L)=0.0
          RAFPL(L)=0.0
          DO 2500 NZ=1,NP(NY,NX)
            RAYSL(NZ,NY,NX)=0.0
            RAYSW(NZ,NY,NX)=0.0
            RAYPL(NZ,NY,NX)=0.0
            RAYPW(NZ,NY,NX)=0.0
            DO 2600 N=1,4
              TSURFY=TSURF(N,L,NZ,NY,NX)*CFX(NZ,NY,NX)
              TSURWY=TSURFB(N,L,NZ,NY,NX)*CFW
              DO 2700 M=1,4
                DO 2750 NN=1,4
                  RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ,NY,NX)
                  RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ,NY,NX)
                  RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
                  PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)+RAPYN
                  PAR(N,M,L,NZ,NY,NX)=PAR(N,M,L,NZ,NY,NX)+RAPYN
                  RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)+TSURFY*RADYN
                  RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)+TSURWY*RADYW
                  RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)+TSURFY*RAPYN
                  RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)+TSURWY*RAPYW
2750            CONTINUE
2700          CONTINUE
2600        CONTINUE
            RAFSL(L)=RAFSL(L)+RAYSL(NZ,NY,NX)*TAUR(NZ,NY,NX)*YAREA
            RAFPL(L)=RAFPL(L)+RAYPL(NZ,NY,NX)*TAUP(NZ,NY,NX)*YAREA
            RADC(NZ,NY,NX)=RADC(NZ,NY,NX)+RAYSL(NZ,NY,NX)+RAYSW(NZ,NY,NX)
            RADP(NZ,NY,NX)=RADP(NZ,NY,NX)+RAYPL(NZ,NY,NX)+RAYPW(NZ,NY,NX)
            TRADC(NY,NX)=TRADC(NY,NX)+RAYSL(NZ,NY,NX)+RAYSW(NZ,NY,NX)
            TRAPC(NY,NX)=TRAPC(NY,NX)+RAYPL(NZ,NY,NX)+RAYPW(NZ,NY,NX)
2500      CONTINUE
        ELSE
          RAFSL(L)=RAFSL(L-1)
          RAFPL(L)=RAFPL(L-1)
          RABSL(L)=RABSL(L-1)
          RABPL(L)=RABPL(L-1)
        ENDIF
2800  CONTINUE
!
      !     RADIATION AT GROUND SURFACE IF NO CANOPY
!
    ELSE
      RASG=ABS(BETAG)*RADS(NY,NX)
      DO 120 N=1,4
        RASG=RASG+ABS(OMEGAG(N,NY,NX))*RADY(NY,NX)
120   CONTINUE
      RADG(NY,NX)=RASG*AREA(3,NU(NY,NX),NY,NX)
      DO 135 NZ=1,NP(NY,NX)
        RADC(NZ,NY,NX)=0.0
        RADP(NZ,NY,NX)=0.0
135   CONTINUE
    ENDIF
!
    !     IF NO RADIATION
!
  ELSE
    RADG(NY,NX)=0.0
    DO 125 NZ=1,NP(NY,NX)
      RADC(NZ,NY,NX)=0.0
      RADP(NZ,NY,NX)=0.0
125 CONTINUE
  ENDIF
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FRADG=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     ARLSS,ARLFS=leaf+stalk area of all PFTs,each PFT
  !
  FRADG(NY,NX)=1.0
  IF(ARLSS(NY,NX).GT.ZEROS(NY,NX))THEN
    FRADPT=1.0-EXP(-0.65*ARLSS(NY,NX)/AREA(3,NU(NY,NX),NY,NX))
    DO 145 NZ=1,NP(NY,NX)
      FRADP(NZ,NY,NX)=FRADPT*ARLFS(NZ,NY,NX)/ARLSS(NY,NX)
      FRADG(NY,NX)=FRADG(NY,NX)-FRADP(NZ,NY,NX)
145 CONTINUE
  ELSE
    FRADG(NY,NX)=1.0
    DO 146 NZ=1,NP(NY,NX)
      FRADP(NZ,NY,NX)=0.0
146 CONTINUE
  ENDIF
  end subroutine MultiLayerSurfaceRadiation
!------------------------------------------------------------------------------------------

  subroutine DivideCanopyLayerByLAI(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: NZ,L
  !     begin_execution
  !
  !     DIVISION OF CANOPY INTO JC LAYERS WITH EQUAL LAI
  !
  !     ZT,ZC=heights of combined canopy,PFT canopy
  !     ZL=height to bottom of each canopy layer
  !     ARLFC,ARSTC=leaf,stalk area of combined canopy
  !     ARLFT,ARSTT=leaf,stalk area of combined canopy layer
  !
  ZT(NY,NX)=0.0
  DO 9685 NZ=1,NP(NY,NX)
    ZT(NY,NX)=AMAX1(ZT(NY,NX),ZC(NZ,NY,NX))
9685  CONTINUE
  ZL(JC,NY,NX)=ZT(NY,NX)+0.01
  ZL1(JC,NY,NX)=ZL(JC,NY,NX)
  ZL1(0,NY,NX)=0.0
  ART=(ARLFC(NY,NX)+ARSTC(NY,NX))/JC
  IF(ART.GT.ZEROS(NY,NX))THEN
    DO 2765 L=JC,2,-1
      ARL=ARLFT(L,NY,NX)+ARSTT(L,NY,NX)
      IF(ARL.GT.1.01*ART)THEN
        DZL=ZL(L,NY,NX)-ZL(L-1,NY,NX)
        ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)+0.5*AMIN1(1.0,(ARL-ART)/ARL)*DZL
      ELSEIF(ARL.LT.0.99*ART)THEN
        ARX=ARLFT(L-1,NY,NX)+ARSTT(L-1,NY,NX)
        DZL=ZL(L-1,NY,NX)-ZL(L-2,NY,NX)
        IF(ARX.GT.ZEROS(NY,NX))THEN
          ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)-0.5*AMIN1(1.0,(ART-ARL)/ARX)*DZL
        ENDIF
      ELSE
        ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)
      ENDIF
      !     IF(J.EQ.12)THEN
      !     WRITE(*,3233)'ZL',I,J,NX,NY,L,ZL1(L,NY,NX),ZL1(L-1,NY,NX)
      !    3,ZL(L,NY,NX),ZL(L-1,NY,NX),ART,ARL,ARX
      !    2,DZL,ARLFC(NY,NX),ARSTC(NY,NX),ARLFT(L,NY,NX),ARSTT(L,NY,NX)
      !    3,ARLFT(L-1,NY,NX),ARSTT(L-1,NY,NX)
      !    3,ZL(JC,NY,NX),ZT(NY,NX),ZC(1,NY,NX)
      !3233  FORMAT(A8,5I4,30E12.4)
      !     ENDIF
2765  CONTINUE
    DO 2770 L=JC,2,-1
      ZL(L-1,NY,NX)=ZL1(L-1,NY,NX)
!     ZL(L-1,NY,NX)=AMAX1(0.0,AMIN1(ZL(L,NY,NX)-1.0E-06
!    2,ZL(L-1,NY,NX)))
2770  CONTINUE
  ENDIF
  end subroutine DivideCanopyLayerByLAI
!------------------------------------------------------------------------------------------

  subroutine CalcBoundaryLayerProperties(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
!     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
!
!     ARLFC,ARSTC=leaf,stalk area of combined canopy
!     DPTHS,DPTH0=snowpack,surface water depths
!     ZT,ZD,ZR=canopy,zero plane displacement,roughness height
!     ZZ=reference height for wind speed
!
  ARLSC=ARLFC(NY,NX)+ARSTC(NY,NX)
  IF(ARLSC.GT.ZEROS(NY,NX) &
    .AND.ZT(NY,NX).GE.DPTHS(NY,NX)-ZERO &
    .AND.ZT(NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
    ARLSG=ARLSC/AREA(3,NU(NY,NX),NY,NX)
    ZX=EXP(-0.5*ARLSG)
    ZY=1.0-ZX
    ZD(NY,NX)=ZT(NY,NX)*AMAX1(0.0,1.0-2.0/ARLSG*ZY)
    ZE=ZT(NY,NX)*AMAX1(0.05,ZX*ZY)
  ELSE
    ZD(NY,NX)=0.0
    ZE=0.0
  ENDIF
  IF(IFLGW.EQ.1)THEN
    ZZ=Z0(NY,NX)+ZT(NY,NX)
  ELSE
    ZZ=AMAX1(Z0(NY,NX),ZD(NY,NX)+2.0)
  ENDIF
  IF(IETYP(NY,NX).GE.0)THEN
    IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      ZR(NY,NX)=AMAX1(0.001,ZE,ZW)
    ELSE
      ZR(NY,NX)=AMAX1(0.001,ZE,ZS(NY,NX))
    ENDIF
!
!     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
!
!     RAB,RAM=biome canopy,minimum isothermal boundary layer resistance
!     UA=wind speed
!     RIB=canopy isothermal Richardson number
!
    RAB(NY,NX)=AMAX1(RAM,(LOG((ZZ-ZD(NY,NX))/ZR(NY,NX)))**2/(0.168*UA(NY,NX)))
    RIB(NY,NX)=1.27E+08*(ZZ-ZR(NY,NX))/(UA(NY,NX)**2*TKA(NY,NX))
  ELSE
    RAB(NY,NX)=RAM
    RIB(NY,NX)=0.0
  ENDIF
  end subroutine CalcBoundaryLayerProperties
!------------------------------------------------------------------------------------------

  subroutine ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)
!
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  integer :: NX,NY
!     begin_execution

  DO 8990 NX=NHW,NHE
    DO 8995 NY=NVN,NVS
      IF(J.EQ.INT(ZNOON(NY,NX)))THEN

        call ApplyMineralFertilizer(I,J,NY,NX)
!
!     SOIL LAYER NUMBER IN WHICH PLANT OR ANIMAL RESIDUES ARE APPLIED
!
        call ApplyPlantAnimalResidue(I,J,NY,NX)
!
!     FERTILIZER UREA, NITRIFICATION INHIBITORS
        call ApplyUreaNitrifierInhibitor(I,J,NY,NX)

      ENDIF
8995  CONTINUE
8990  CONTINUE
  end subroutine ApplyFertilizerAtNoon
!------------------------------------------------------------------------------------------

  subroutine ApplyUreaNitrifierInhibitor(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L
!     begin_execution
!
!     IYTYP=fertilizer release type from fertilizer input file
!     FERT=fertilizer type from fertilizer input file
!     IUTYP=urea hydrolysis inhibitor type (1=no,2=yes)
!     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity
!     ZNFN0,ZNFNI=initial,current nitrification inhibition activity
!
  IF(FERT(3,I,NY,NX).GT.0.0.OR.FERT(7,I,NY,NX).GT.0.0)THEN
    IF(IYTYP(0,I,NY,NX).EQ.0)THEN
      IUTYP(NY,NX)=0
    ELSEIF(IYTYP(0,I,NY,NX).EQ.1.OR.IYTYP(0,I,NY,NX).EQ.3)THEN
      IUTYP(NY,NX)=1
    ELSE
      IUTYP(NY,NX)=2
    ENDIF
    DO 9964 L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNHU0(L,NY,NX)=1.0
        ZNHUI(L,NY,NX)=1.0
      ELSE
        ZNHU0(L,NY,NX)=0.0
        ZNHUI(L,NY,NX)=0.0
      ENDIF
9964  CONTINUE
  ENDIF
  IF(IYTYP(0,I,NY,NX).EQ.3.OR.IYTYP(0,I,NY,NX).EQ.4)THEN
    DO 9965 L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNFN0(L,NY,NX)=1.0
        ZNFNI(L,NY,NX)=1.0
      ELSE
        ZNFN0(L,NY,NX)=0.0
        ZNFNI(L,NY,NX)=0.0
      ENDIF
9965  CONTINUE
  ENDIF
  end subroutine ApplyUreaNitrifierInhibitor
!------------------------------------------------------------------------------------------

  subroutine ApplyPlantAnimalResidue(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L,K,M,N,NN,NGL
!     begin_execution
!     LFDPTH=layer number
!
  IF(OFC(1)+OFC(2).GT.0.0)THEN
    DO 2985 L=0,JZ
      FDPTHM=FDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      IF(FDPTHM.LE.0.0)THEN
        LFDPTH=0
        exit
      ELSEIF(CDPTH(L,NY,NX).GE.FDPTHM)THEN
        LFDPTH=L
        exit
      ENDIF
2985  CONTINUE
!
!     ALLOCATION OF PLANT RESIDUE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     CFOSC=fraction of litter allocated to protein(1)
!     soluble CH2O(2), cellulose(3) and lignin(4)
!     ITYPE=litter type entered in fertilizer input file
!
!     MAIZE
!
    IF(IYTYP(1,I,NY,NX).EQ.1)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.080
      CFOSC(2,1,LFDPTH,NY,NX)=0.245
      CFOSC(3,1,LFDPTH,NY,NX)=0.613
      CFOSC(4,1,LFDPTH,NY,NX)=0.062
!
!     WHEAT
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.2)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.125
      CFOSC(2,1,LFDPTH,NY,NX)=0.171
      CFOSC(3,1,LFDPTH,NY,NX)=0.560
      CFOSC(4,1,LFDPTH,NY,NX)=0.144
!
!     SOYBEAN
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.3)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.138
      CFOSC(2,1,LFDPTH,NY,NX)=0.426
      CFOSC(3,1,LFDPTH,NY,NX)=0.316
      CFOSC(4,1,LFDPTH,NY,NX)=0.120
!
!     OLD STRAW
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.4)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.075
      CFOSC(2,1,LFDPTH,NY,NX)=0.125
      CFOSC(3,1,LFDPTH,NY,NX)=0.550
      CFOSC(4,1,LFDPTH,NY,NX)=0.250
!
!     STRAW
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.5)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.036
      CFOSC(2,1,LFDPTH,NY,NX)=0.044
      CFOSC(3,1,LFDPTH,NY,NX)=0.767
      CFOSC(4,1,LFDPTH,NY,NX)=0.153
!
!     COMPOST
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.6)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.143
      CFOSC(2,1,LFDPTH,NY,NX)=0.015
      CFOSC(3,1,LFDPTH,NY,NX)=0.640
      CFOSC(4,1,LFDPTH,NY,NX)=0.202
!
!     GREEN MANURE
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.7)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.202
      CFOSC(2,1,LFDPTH,NY,NX)=0.013
      CFOSC(3,1,LFDPTH,NY,NX)=0.560
      CFOSC(4,1,LFDPTH,NY,NX)=0.225
!
!     SIMPLE SUBSTRATE
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.10)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.000
      CFOSC(2,1,LFDPTH,NY,NX)=1.000
      CFOSC(3,1,LFDPTH,NY,NX)=0.000
      CFOSC(4,1,LFDPTH,NY,NX)=0.000
    ELSE
      CFOSC(1,1,LFDPTH,NY,NX)=0.075
      CFOSC(2,1,LFDPTH,NY,NX)=0.125
      CFOSC(3,1,LFDPTH,NY,NX)=0.550
      CFOSC(4,1,LFDPTH,NY,NX)=0.250
    ENDIF
!
!     ALLOCATION OF ANIMAL MANURE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     RUMINANT
!
    IF(IYTYP(2,I,NY,NX).EQ.1)THEN
      CFOSC(1,2,LFDPTH,NY,NX)=0.036
      CFOSC(2,2,LFDPTH,NY,NX)=0.044
      CFOSC(3,2,LFDPTH,NY,NX)=0.630
      CFOSC(4,2,LFDPTH,NY,NX)=0.290
!
!     NON-RUMINANT
!
    ELSEIF(IYTYP(2,I,NY,NX).EQ.2)THEN
      CFOSC(1,2,LFDPTH,NY,NX)=0.138
      CFOSC(2,2,LFDPTH,NY,NX)=0.401
      CFOSC(3,2,LFDPTH,NY,NX)=0.316
      CFOSC(4,2,LFDPTH,NY,NX)=0.145
!
!     GRAZING
!
    ELSEIF(IYTYP(2,I,NY,NX).EQ.3)THEN
      CFOSC(1,2,LFDPTH,NY,NX)=0.036
      CFOSC(2,2,LFDPTH,NY,NX)=0.044
      CFOSC(3,2,LFDPTH,NY,NX)=0.630
      CFOSC(4,2,LFDPTH,NY,NX)=0.290
!
!     OTHER
!
    ELSE
      CFOSC(1,2,LFDPTH,NY,NX)=0.138
      CFOSC(2,2,LFDPTH,NY,NX)=0.401
      CFOSC(3,2,LFDPTH,NY,NX)=0.316
      CFOSC(4,2,LFDPTH,NY,NX)=0.145
    ENDIF
!
!     DISTRIBUTE RESIDUE APPLICATION AMONG COMPONENTS OF RESIDUE COMPLEX
!
!     OFC,OFN,OFP=litter C,N,P application from fertilizer file
!
    DO 2965 K=1,2
      OSCI=OFC(K)*AREA(3,LFDPTH,NY,NX)
      OSNI=OFN(K)*AREA(3,LFDPTH,NY,NX)
      OSPI=OFP(K)*AREA(3,LFDPTH,NY,NX)
      IF(BKVL(LFDPTH,NY,NX).GT.ZEROS(NY,NX))THEN
        CORGCX=OSCI/BKVL(LFDPTH,NY,NX)
      ELSE
        CORGCX=0.55E+06
      ENDIF
      OSCX=0.0
      OSNX=0.0
      OSPX=0.0
!
!     BIOMASSES OF MICROBIAL POPULATIONS IN RESIDUE
!
!     OMC,OMN,OMP=microbial biomass in litter application
!     OMCI=microbial biomass content in litter
!     OMCF,OMCA=hetero,autotrophic biomass composition in litter
!
      DO 2960 N=1,7
        DO NGL=1,JG
          DO 2961 M=1,3
            OMC1=AMAX1(0.0,AMIN1(OSCI*OMCI(M+(NGL-1)*3,K)*OMCF(N),OSCI-OSCX))
            OMN1=AMAX1(0.0,AMIN1(OMC1*CNOMC(M,NGL,N,K),OSNI-OSNX))
            OMP1=AMAX1(0.0,AMIN1(OMC1*CPOMC(M,NGL,N,K),OSPI-OSPX))
            OMC(M,NGL,N,K,LFDPTH,NY,NX)=OMC(M,NGL,N,K,LFDPTH,NY,NX)+OMC1
            OMN(M,NGL,N,K,LFDPTH,NY,NX)=OMN(M,NGL,N,K,LFDPTH,NY,NX)+OMN1
            OMP(M,NGL,N,K,LFDPTH,NY,NX)=OMP(M,NGL,N,K,LFDPTH,NY,NX)+OMP1
!     WRITE(*,2345)'OMCI',I,J,LFDPTH,K,N,M
!    2,OMC1,OMN1,OMP1,OSCI,OMCI(M,K)
!    2,OMCF(N),OSCX,CNOMC(M,NGL,N,K),CPOMC(M,NGL,N,K),OSNI,OSPI
!    2,OMC(M,NGL,N,K,LFDPTH,NY,NX),OMN(M,NGL,N,K,LFDPTH,NY,NX)
!2345  FORMAT(A8,6I4,20E12.4)
            OSCX=OSCX+OMC1
            OSNX=OSNX+OMN1
            OSPX=OSPX+OMP1
            DO 2962 NN=1,7
              OMC(M,NGL,NN,5,LFDPTH,NY,NX)=OMC(M,NGL,NN,5,LFDPTH,NY,NX)+OMC1*OMCA(NN)
              OMN(M,NGL,NN,5,LFDPTH,NY,NX)=OMN(M,NGL,NN,5,LFDPTH,NY,NX)+OMN1*OMCA(NN)
              OMP(M,NGL,NN,5,LFDPTH,NY,NX)=OMP(M,NGL,NN,5,LFDPTH,NY,NX)+OMP1*OMCA(NN)
!     WRITE(*,2346)'OMCA',I,J,LFDPTH,K,N,NN,M
!    2,OMC1,OMCA(NN),OMC(M,NN,5,LFDPTH,NY,NX)
!2346  FORMAT(A8,7I4,20E12.4)
              OSCX=OSCX+OMC1*OMCA(NN)
              OSNX=OSNX+OMN1*OMCA(NN)
              OSPX=OSPX+OMP1*OMCA(NN)
2962        CONTINUE
2961      CONTINUE
        ENDDO
2960  CONTINUE
!
!     DOC, DON AND DOP IN RESIDUE
!
!     OQC,OQN,OQP=DOC,DON,DOP in litter
!
      OQC1=AMIN1(0.1*OSCX,OSCI-OSCX)
      OQN1=AMIN1(0.1*OSNX,OSNI-OSNX)
      OQP1=AMIN1(0.1*OSPX,OSPI-OSPX)
      OQC(K,LFDPTH,NY,NX)=OQC(K,LFDPTH,NY,NX)+OQC1
      OQN(K,LFDPTH,NY,NX)=OQN(K,LFDPTH,NY,NX)+OQN1
      OQP(K,LFDPTH,NY,NX)=OQP(K,LFDPTH,NY,NX)+OQP1
!
!     REMAINDER DISTRIBUTED TO RESIDUE FRACTIONS
!
!     OSC,OSN,OSP,OSA=SOC,SON,SOP,colonized SOC in litter
!     VOLT=litter volume
!     UORGF,UFERTN,UFERTP=accumulated litter C,N,P application
!     TNBP=accumulated net biome productivity
!
      OSCX=OSCX+OQC1
      OSNX=OSNX+OQN1
      OSPX=OSPX+OQP1
      CNOFT=0.0
      CPOFT=0.0
      IF(OSCI-OSCX.GT.ZEROS(NY,NX))THEN
        RNT=0.0
        RPT=0.0
        DO 965 M=1,4
          RNT=RNT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*CNOFC(M,K)
          RPT=RPT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*CPOFC(M,K)
965     CONTINUE
        FRNT=(OSNI-OSNX)/RNT
        FRPT=(OSPI-OSPX)/RPT
        DO 970 M=1,4
          CNOF(M)=CNOFC(M,K)*FRNT
          CPOF(M)=CPOFC(M,K)*FRPT
          CNOFT=CNOFT+CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)
          CPOFT=CPOFT+CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)
970     CONTINUE
      ELSE
        DO 975 M=1,4
          CNOF(M)=0.0
          CPOF(M)=0.0
975     CONTINUE
      ENDIF
      DO 2970 M=1,4
        OSC1=CFOSC(M,K,LFDPTH,NY,NX)*(OSCI-OSCX)
        IF(CNOFT.GT.ZERO)THEN
          OSN1=CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)/CNOFT*(OSNI-OSNX)
        ELSE
          OSN1=0.0
        ENDIF
        IF(CPOFT.GT.ZERO)THEN
          OSP1=CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)/CPOFT*(OSPI-OSPX)
        ELSE
          OSP1=0.0
        ENDIF
        OSC(M,K,LFDPTH,NY,NX)=OSC(M,K,LFDPTH,NY,NX)+OSC1
        DO NGL=1,JG
          OSA(M,K,LFDPTH,NY,NX)=OSA(M,K,LFDPTH,NY,NX)+OSC1*OMCI(1+(NGL-1)*3,K)
        ENDDO
        OSN(M,K,LFDPTH,NY,NX)=OSN(M,K,LFDPTH,NY,NX)+OSN1
        OSP(M,K,LFDPTH,NY,NX)=OSP(M,K,LFDPTH,NY,NX)+OSP1
        IF(LFDPTH.EQ.0)THEN
          VOLT(LFDPTH,NY,NX)=VOLT(LFDPTH,NY,NX)+OSC1*1.0E-06/BKRS(1)
        ENDIF
2970  CONTINUE
      TORGF=TORGF+OSCI
      TORGN=TORGN+OSNI
      TORGP=TORGP+OSPI
      UORGF(NY,NX)=UORGF(NY,NX)+OSCI
      UFERTN(NY,NX)=UFERTN(NY,NX)+OSNI
      UFERTP(NY,NX)=UFERTP(NY,NX)+OSPI
      IF(IYTYP(2,I,NY,NX).LT.3)THEN
        TNBP(NY,NX)=TNBP(NY,NX)+OSCI
      ENDIF
2965  CONTINUE
  ENDIF
  end subroutine ApplyPlantAnimalResidue
!------------------------------------------------------------------------------------------

  subroutine ApplyMineralFertilizer(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L
!     begin_execution
!
!     NH4,NH3,UREA,NO3 FERTILIZER APPLICATION
!
!     *A,*B=broadcast,banded
!     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
!
  Z4A=FERT(1,I,NY,NX)
  Z3A=FERT(2,I,NY,NX)
  ZUA=FERT(3,I,NY,NX)
  ZOA=FERT(4,I,NY,NX)
  Z4B=FERT(5,I,NY,NX)
  Z3B=FERT(6,I,NY,NX)
  ZUB=FERT(7,I,NY,NX)
  ZOB=FERT(8,I,NY,NX)
!
!     MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE
!
!     PM*,PH*=Ca(H2PO4)2,apatite
!
  PMA=FERT(9,I,NY,NX)
  PMB=FERT(10,I,NY,NX)
  PHA=FERT(11,I,NY,NX)
!
!     LIME AND GYPSUM
!
!     CAC,CAS=CaCO3,CaSO4
!
  CAC=FERT(12,I,NY,NX)
  CAS=FERT(13,I,NY,NX)
!
!     PLANT(1) AND ANIMAL(2) RESIDUE C, N AND P
!
  OFC(1)=FERT(14,I,NY,NX)
  OFN(1)=FERT(15,I,NY,NX)
  OFP(1)=FERT(16,I,NY,NX)
  OFC(2)=FERT(17,I,NY,NX)
  OFN(2)=FERT(18,I,NY,NX)
  OFP(2)=FERT(19,I,NY,NX)
!
!     SOIL LAYER NUMBER AT DEPTH OF FERTILIZER APPLICATION
!
!     LFDPTH=layer number
!     CVRDF=fraction of fertilizer applied to surface litter
!
  IF(Z4A+Z3A+ZUA+ZOA+Z4B+Z3B+ZUB+ZOB+PMA+PMB+PHA+CAC+CAS.GT.0.0)THEN
    FDPTHF=FDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
    IF(FDPTHF.LE.0.0.AND.test_aeqb(Z4B+Z3B+ZUB+ZOB+PMB,0._r8))THEN
      LFDPTH=0
      CVRDF=1.0-EXP(-0.8E-02*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))
    ELSE
      DO 65 L=NUI(NY,NX),JZ
        IF(CDPTH(L,NY,NX).GE.FDPTHF)THEN
          LFDPTH=L
          CVRDF=1.0
          exit
        ENDIF
65    CONTINUE
    ENDIF
    BAREF=1.0-CVRDF
!
!     RESET WIDTH AND DEPTH OF NH4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWN=width of NH4 band row
!     DPNHB,WDNHB=depth,width of NH4 band
!     VLNHB,VLNH4=soil volume in NH4 band,non-band
!
    IF((Z4B+Z3B+ZUB.GT.0.0).OR.((ZNH4B(LFDPTH,NY,NX).GT.0.0 &
      .OR.ZNH3B(LFDPTH,NY,NX).GT.0.0).AND.IFNHB(NY,NX).EQ.0))THEN
      IFNHB(NY,NX)=1
      ROWN(NY,NX)=ROWI(I,NY,NX)
      DO 50 L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          DPNHB(L,NY,NX)=DLYR(3,L,NY,NX)
          WDNHB(L,NY,NX)=0.0
        ELSEIF(L.EQ.LFDPTH)THEN
          DPNHB(L,NY,NX)=AMAX1(0.025,FDPTHF-CDPTH(L-1,NY,NX))
          WDNHB(L,NY,NX)=AMIN1(0.025,ROWN(NY,NX))
        ELSE
          DPNHB(L,NY,NX)=0.0
          WDNHB(L,NY,NX)=0.0
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          VLNHB(L,NY,NX)=AMIN1(0.999,WDNHB(L,NY,NX)/ROWN(NY,NX) &
            *DPNHB(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          VLNHB(L,NY,NX)=0.0
        ENDIF
        VLNH4(L,NY,NX)=1.0-VLNHB(L,NY,NX)
        ZNH4T=ZNH4S(L,NY,NX)+ZNH4B(L,NY,NX)
        ZNH3T=ZNH3S(L,NY,NX)+ZNH3B(L,NY,NX)
        XN4T=XN4(L,NY,NX)+XNB(L,NY,NX)
        ZNH4S(L,NY,NX)=ZNH4T*VLNH4(L,NY,NX)
        ZNH3S(L,NY,NX)=ZNH3T*VLNH4(L,NY,NX)
        ZNH4B(L,NY,NX)=ZNH4T*VLNHB(L,NY,NX)
        ZNH3B(L,NY,NX)=ZNH3T*VLNHB(L,NY,NX)
        XN4(L,NY,NX)=XN4T*VLNH4(L,NY,NX)
        XNB(L,NY,NX)=XN4T*VLNHB(L,NY,NX)
50    CONTINUE
      DPNH4(NY,NX)=DPNHB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF NO3 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWO=width of NO3 band row
!     DPNOB,WDNOB=depth,width of NO3 band
!     VLNOB,VLNO3=soil volume in NO3 band,non-band
!
    IF((Z4B+Z3B+ZUB+ZOB.GT.0.0).OR.((ZNO3B(LFDPTH,NY,NX).GT.0.0 &
      .OR.ZNO2B(LFDPTH,NY,NX).GT.0.0).AND.IFNOB(NY,NX).EQ.0))THEN
      IFNOB(NY,NX)=1
      ROWO(NY,NX)=ROWI(I,NY,NX)
      DO 45 L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          DPNOB(L,NY,NX)=DLYR(3,L,NY,NX)
          WDNOB(L,NY,NX)=0.0
        ELSEIF(L.EQ.LFDPTH)THEN
          DPNOB(L,NY,NX)=AMAX1(0.01,FDPTHF-CDPTH(L-1,NY,NX))
          WDNOB(L,NY,NX)=AMIN1(0.01,ROWO(NY,NX))
        ELSE
          DPNOB(L,NY,NX)=0.0
          WDNOB(L,NY,NX)=0.0
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          VLNOB(L,NY,NX)=AMIN1(0.999,WDNOB(L,NY,NX)/ROWO(NY,NX) &
            *DPNOB(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          VLNOB(L,NY,NX)=0.0
        ENDIF
        VLNO3(L,NY,NX)=1.0-VLNOB(L,NY,NX)
        ZNO3T=ZNO3S(L,NY,NX)+ZNO3B(L,NY,NX)
        ZNO2T=ZNO2S(L,NY,NX)+ZNO2B(L,NY,NX)
        ZNO3S(L,NY,NX)=ZNO3T*VLNO3(L,NY,NX)
        ZNO2S(L,NY,NX)=ZNO2T*VLNO3(L,NY,NX)
        ZNO3B(L,NY,NX)=ZNO3T*VLNOB(L,NY,NX)
        ZNO2B(L,NY,NX)=ZNO2T*VLNOB(L,NY,NX)
45    CONTINUE
      DPNO3(NY,NX)=DPNOB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF PO4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWP=width of H2PO4 band row
!     DPPOB,WDPOB=depth,width of H2PO4 band
!     VLPOB,VLPO4=soil volume in H2PO4 band,non-band
!
    IF((PMB.GT.0.0).OR.(H2POB(LFDPTH,NY,NX).GT.0.0 &
      .AND.IFPOB(NY,NX).EQ.0))THEN
      IFPOB(NY,NX)=1
      ROWP(NY,NX)=ROWI(I,NY,NX)
      DO 40 L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          DPPOB(L,NY,NX)=DLYR(3,L,NY,NX)
          WDPOB(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
        ELSEIF(L.EQ.LFDPTH)THEN
          DPPOB(L,NY,NX)=AMAX1(0.01,FDPTHF-CDPTH(L-1,NY,NX))
          WDPOB(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
        ELSE
          DPPOB(L,NY,NX)=0.0
          WDPOB(L,NY,NX)=0.0
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          VLPOB(L,NY,NX)=AMIN1(0.999,WDPOB(L,NY,NX)/ROWP(NY,NX) &
          *DPPOB(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          VLPOB(L,NY,NX)=0.0
        ENDIF
        VLPO4(L,NY,NX)=1.0-VLPOB(L,NY,NX)
        H0PO4T=H0PO4(L,NY,NX)+H0POB(L,NY,NX)
        H1PO4T=H1PO4(L,NY,NX)+H1POB(L,NY,NX)
        H2PO4T=H2PO4(L,NY,NX)+H2POB(L,NY,NX)
        H3PO4T=H3PO4(L,NY,NX)+H3POB(L,NY,NX)
        ZFE1PT=ZFE1P(L,NY,NX)+ZFE1PB(L,NY,NX)
        ZFE2PT=ZFE2P(L,NY,NX)+ZFE2PB(L,NY,NX)
        ZCA0PT=ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX)
      ZCA1PT=ZCA1P(L,NY,NX)+ZCA1PB(L,NY,NX)
      ZCA2PT=ZCA2P(L,NY,NX)+ZCA2PB(L,NY,NX)
      ZMG1PT=ZMG1P(L,NY,NX)+ZMG1PB(L,NY,NX)
      XOH0T=XOH0(L,NY,NX)+XOH0B(L,NY,NX)
      XOH1T=XOH1(L,NY,NX)+XOH1B(L,NY,NX)
      XOH2T=XOH2(L,NY,NX)+XOH2B(L,NY,NX)
      XH1PT=XH1P(L,NY,NX)+XH1PB(L,NY,NX)
      XH2PT=XH2P(L,NY,NX)+XH2PB(L,NY,NX)
      PALPOT=PALPO(L,NY,NX)+PALPB(L,NY,NX)
      PFEPOT=PFEPO(L,NY,NX)+PFEPB(L,NY,NX)
      PCAPDT=PCAPD(L,NY,NX)+PCPDB(L,NY,NX)
      PCAPHT=PCAPH(L,NY,NX)+PCPHB(L,NY,NX)
      PCAPMT=PCAPM(L,NY,NX)+PCPMB(L,NY,NX)
      H0PO4(L,NY,NX)=H0PO4T*VLPO4(L,NY,NX)
      H1PO4(L,NY,NX)=H1PO4T*VLPO4(L,NY,NX)
      H2PO4(L,NY,NX)=H2PO4T*VLPO4(L,NY,NX)
      H3PO4(L,NY,NX)=H3PO4T*VLPO4(L,NY,NX)
      ZFE1P(L,NY,NX)=ZFE1PT*VLPO4(L,NY,NX)
      ZFE2P(L,NY,NX)=ZFE2PT*VLPO4(L,NY,NX)
      ZCA0P(L,NY,NX)=ZCA0PT*VLPO4(L,NY,NX)
      ZCA1P(L,NY,NX)=ZCA1PT*VLPO4(L,NY,NX)
      ZCA2P(L,NY,NX)=ZCA2PT*VLPO4(L,NY,NX)
      ZMG1P(L,NY,NX)=ZMG1PT*VLPO4(L,NY,NX)
      H0POB(L,NY,NX)=H0PO4T*VLPOB(L,NY,NX)
      H1POB(L,NY,NX)=H1PO4T*VLPOB(L,NY,NX)
      H2POB(L,NY,NX)=H2PO4T*VLPOB(L,NY,NX)
      H3POB(L,NY,NX)=H3PO4T*VLPOB(L,NY,NX)
      ZFE1PB(L,NY,NX)=ZFE1PT*VLPOB(L,NY,NX)
      ZFE2PB(L,NY,NX)=ZFE2PT*VLPOB(L,NY,NX)
      ZCA0PB(L,NY,NX)=ZCA0PT*VLPOB(L,NY,NX)
      ZCA1PB(L,NY,NX)=ZCA1PT*VLPOB(L,NY,NX)
      ZCA2PB(L,NY,NX)=ZCA2PT*VLPOB(L,NY,NX)
      ZMG1PB(L,NY,NX)=ZMG1PT*VLPOB(L,NY,NX)
      XOH0(L,NY,NX)=XOH0T*VLPO4(L,NY,NX)
      XOH1(L,NY,NX)=XOH1T*VLPO4(L,NY,NX)
      XOH2(L,NY,NX)=XOH2T*VLPO4(L,NY,NX)
      XH1P(L,NY,NX)=XH1PT*VLPO4(L,NY,NX)
      XH2P(L,NY,NX)=XH2PT*VLPO4(L,NY,NX)
      XOH0B(L,NY,NX)=XOH0T*VLPOB(L,NY,NX)
      XOH1B(L,NY,NX)=XOH1T*VLPOB(L,NY,NX)
      XOH2B(L,NY,NX)=XOH2T*VLPOB(L,NY,NX)
      XH1PB(L,NY,NX)=XH1PT*VLPOB(L,NY,NX)
      XH2PB(L,NY,NX)=XH2PT*VLPOB(L,NY,NX)
      PALPO(L,NY,NX)=PALPOT*VLPO4(L,NY,NX)
      PFEPO(L,NY,NX)=PFEPOT*VLPO4(L,NY,NX)
      PCAPD(L,NY,NX)=PCAPDT*VLPO4(L,NY,NX)
      PCAPH(L,NY,NX)=PCAPHT*VLPO4(L,NY,NX)
      PCAPM(L,NY,NX)=PCAPMT*VLPO4(L,NY,NX)
      PALPB(L,NY,NX)=PALPOT*VLPOB(L,NY,NX)
      PFEPB(L,NY,NX)=PFEPOT*VLPOB(L,NY,NX)
      PCPDB(L,NY,NX)=PCAPDT*VLPOB(L,NY,NX)
      PCPHB(L,NY,NX)=PCAPHT*VLPOB(L,NY,NX)
      PCPMB(L,NY,NX)=PCAPMT*VLPOB(L,NY,NX)
40    CONTINUE
      DPPO4(NY,NX)=DPPOB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
      ENDIF
!
!     UPDATE STATE VARIABLES FOR BROADCAST AND BANDED FERTILIZER
!     NH4, NH3, UREA, NO3, PO4, LIME AND GYPSUM IN SOIL
!     AND CONVERT FROM G TO MOLE
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=bdcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH2PO4,CaHPO4,apatite in non-band
!     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH2PO4,CaHPO4,apatite in band
!     PCACO,PCASO=precipitated CaCO3,CaSO4
!
      Z4AX=Z4A*AREA(3,LFDPTH,NY,NX)/14.0
      Z3AX=Z3A*AREA(3,LFDPTH,NY,NX)/14.0
      ZUAX=ZUA*AREA(3,LFDPTH,NY,NX)/14.0
      ZOAX=ZOA*AREA(3,LFDPTH,NY,NX)/14.0
      Z4BX=Z4B*AREA(3,LFDPTH,NY,NX)/14.0
      Z3BX=Z3B*AREA(3,LFDPTH,NY,NX)/14.0
      ZUBX=ZUB*AREA(3,LFDPTH,NY,NX)/14.0
      ZOBX=ZOB*AREA(3,LFDPTH,NY,NX)/14.0
      PMAX=PMA*AREA(3,LFDPTH,NY,NX)/62.0
      PMBX=PMB*AREA(3,LFDPTH,NY,NX)/62.0
      PHAX=PHA*AREA(3,LFDPTH,NY,NX)/93.0
      CACX=CAC*AREA(3,LFDPTH,NY,NX)/40.0
      CASX=CAS*AREA(3,LFDPTH,NY,NX)/40.0
      ZNH4FA(LFDPTH,NY,NX)=ZNH4FA(LFDPTH,NY,NX)+Z4AX*CVRDF
      ZNHUFA(LFDPTH,NY,NX)=ZNHUFA(LFDPTH,NY,NX)+ZUAX*CVRDF
      ZNO3FA(LFDPTH,NY,NX)=ZNO3FA(LFDPTH,NY,NX)+ZOAX*CVRDF
      ZNH4FB(LFDPTH,NY,NX)=ZNH4FB(LFDPTH,NY,NX)+Z4BX*CVRDF
      ZNHUFB(LFDPTH,NY,NX)=ZNHUFB(LFDPTH,NY,NX)+ZUBX*CVRDF
      ZNO3FB(LFDPTH,NY,NX)=ZNO3FB(LFDPTH,NY,NX)+ZOBX*CVRDF
      PCAPM(LFDPTH,NY,NX)=PCAPM(LFDPTH,NY,NX) &
      +PMAX*VLPO4(LFDPTH,NY,NX)*CVRDF
      PCPMB(LFDPTH,NY,NX)=PCPMB(LFDPTH,NY,NX) &
      +PMAX*VLPOB(LFDPTH,NY,NX)*CVRDF+PMBX*CVRDF
      PCAPH(LFDPTH,NY,NX)=PCAPH(LFDPTH,NY,NX) &
      +PHAX*VLPO4(LFDPTH,NY,NX)*CVRDF
      PCPHB(LFDPTH,NY,NX)=PCPHB(LFDPTH,NY,NX) &
      +PHAX*VLPOB(LFDPTH,NY,NX)*CVRDF
      IF(LFDPTH.EQ.0)THEN
      ZNH4FA(NU(NY,NX),NY,NX)=ZNH4FA(NU(NY,NX),NY,NX)+Z4AX*BAREF
      ZNH3FA(NU(NY,NX),NY,NX)=ZNH3FA(NU(NY,NX),NY,NX)+Z3AX
      ZNHUFA(NU(NY,NX),NY,NX)=ZNHUFA(NU(NY,NX),NY,NX)+ZUAX*BAREF
      ZNO3FA(NU(NY,NX),NY,NX)=ZNO3FA(NU(NY,NX),NY,NX)+ZOAX*BAREF
      ZNH4FB(NU(NY,NX),NY,NX)=ZNH4FB(NU(NY,NX),NY,NX)+Z4BX*BAREF
      ZNH3FB(NU(NY,NX),NY,NX)=ZNH3FB(NU(NY,NX),NY,NX)+Z3BX
      ZNHUFB(NU(NY,NX),NY,NX)=ZNHUFB(NU(NY,NX),NY,NX)+ZUBX*BAREF
      ZNO3FB(NU(NY,NX),NY,NX)=ZNO3FB(NU(NY,NX),NY,NX)+ZOBX*BAREF
      PCAPM(NU(NY,NX),NY,NX)=PCAPM(NU(NY,NX),NY,NX) &
      +PMAX*VLPO4(NU(NY,NX),NY,NX)*BAREF
      PCPMB(NU(NY,NX),NY,NX)=PCPMB(NU(NY,NX),NY,NX) &
      +PMAX*VLPOB(NU(NY,NX),NY,NX)*BAREF+PMBX*BAREF
      PCAPH(NU(NY,NX),NY,NX)=PCAPH(NU(NY,NX),NY,NX) &
      +PHAX*VLPO4(NU(NY,NX),NY,NX)*BAREF
      PCPHB(NU(NY,NX),NY,NX)=PCPHB(NU(NY,NX),NY,NX) &
      +PHAX*VLPOB(NU(NY,NX),NY,NX)*BAREF
      ELSE
      ZNH3FA(LFDPTH,NY,NX)=ZNH3FA(LFDPTH,NY,NX)+Z3AX*CVRDF
      ZNH3FB(LFDPTH,NY,NX)=ZNH3FB(LFDPTH,NY,NX)+Z3BX*CVRDF
      ENDIF
      PCACO(NU(NY,NX),NY,NX)=PCACO(NU(NY,NX),NY,NX)+CACX
      PCASO(NU(NY,NX),NY,NX)=PCASO(NU(NY,NX),NY,NX)+CASX
      TZIN=TZIN+14.0*(Z4AX+Z3AX+ZUAX+ZOAX+Z4BX+Z3BX+ZUBX+ZOBX)
      TPIN=TPIN+62.0*(PMAX+PMBX)+93.0*PHAX
      TIONIN=TIONIN+2.0*(CACX+CASX)
      UFERTN(NY,NX)=UFERTN(NY,NX)+14.0*(Z4AX+Z4BX+Z3AX+Z3BX &
      +ZUAX+ZUBX+ZOAX+ZOBX)
      UFERTP(NY,NX)=UFERTP(NY,NX)+62.0*(PMAX+PMBX)+93.0*PHAX
      ENDIF
      end subroutine ApplyMineralFertilizer
!------------------------------------------------------------------------------------------

      subroutine GetChemicalConcsInSoil(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX
!     begin_execution
!     CALCULATE SOIL CONCENTRATIONS OF SOLUTES, GASES
!
!     THETW,THETI,THETP=soil micropore water,ice,air concentration
!     THETPZ=soil micropore+macropore air concn for output
!
      IF(VOLX(L,NY,NX).LE.ZEROS(NY,NX))THEN
      THETW(L,NY,NX)=POROS(L,NY,NX)
      THETI(L,NY,NX)=0.0
      THETP(L,NY,NX)=0.0
      ELSE
      THETW(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX) &
      ,VOLW(L,NY,NX)/VOLY(L,NY,NX)))
      THETI(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX) &
      ,VOLI(L,NY,NX)/VOLY(L,NY,NX)))
      THETP(L,NY,NX)=AMAX1(0.0,VOLP(L,NY,NX)/VOLY(L,NY,NX))
      ENDIF
      THETPZ(L,NY,NX)=AMAX1(0.0,POROS(L,NY,NX)-THETW(L,NY,NX) &
      -THETI(L,NY,NX))
!     IF(L.EQ.7)THEN
!     WRITE(*,1117)'BKDS',I,J,L
!    3,BKDS(L,NY,NX),BKDSI(L,NY,NX),BKVL(L,NY,NX)
!    4,DLYR(3,L,NY,NX),DLYRI(3,L,NY,NX)
!    2,CORGC(L,NY,NX),ORGC(L,NY,NX)
!    6,VOLT(L,NY,NX),VOLX(L,NY,NX),VOLA(L,NY,NX)
!    7,VOLY(L,NY,NX),VOLW(L,NY,NX),VOLI(L,NY,NX)
!    7,VOLWH(L,NY,NX),VOLIH(L,NY,NX)
!    8,THETWZ(L,NY,NX),THETIZ(L,NY,NX)
!    9,THETWZ(L,NY,NX)+THETIZ(L,NY,NX)
!1117  FORMAT(A8,3I4,30E14.6)
!     ENDIF
!
!     GAS CONCENTRATIONS
!
!     C*G=soil gas gaseous concentration
!     C*S=soil gas aqueous concentration
!
      IF(THETP(L,NY,NX).GT.THETX)THEN
      CCO2G(L,NY,NX)=AMAX1(0.0,CO2G(L,NY,NX)/VOLP(L,NY,NX))
      CCH4G(L,NY,NX)=AMAX1(0.0,CH4G(L,NY,NX)/VOLP(L,NY,NX))
      COXYG(L,NY,NX)=AMAX1(0.0,OXYG(L,NY,NX)/VOLP(L,NY,NX))
      CZ2GG(L,NY,NX)=AMAX1(0.0,Z2GG(L,NY,NX)/VOLP(L,NY,NX))
      CZ2OG(L,NY,NX)=AMAX1(0.0,Z2OG(L,NY,NX)/VOLP(L,NY,NX))
      CNH3G(L,NY,NX)=AMAX1(0.0,ZNH3G(L,NY,NX)/VOLP(L,NY,NX))
      CH2GG(L,NY,NX)=AMAX1(0.0,H2GG(L,NY,NX)/VOLP(L,NY,NX))
      ELSE
      CCO2G(L,NY,NX)=0.0
      CCH4G(L,NY,NX)=0.0
      COXYG(L,NY,NX)=0.0
      CZ2GG(L,NY,NX)=0.0
      CZ2OG(L,NY,NX)=0.0
      CNH3G(L,NY,NX)=0.0
      CH2GG(L,NY,NX)=0.0
      ENDIF
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CCO2S(L,NY,NX)=AMAX1(0.0,CO2S(L,NY,NX)/VOLW(L,NY,NX))
      CCH4S(L,NY,NX)=AMAX1(0.0,CH4S(L,NY,NX)/VOLW(L,NY,NX))
      COXYS(L,NY,NX)=AMAX1(0.0,OXYS(L,NY,NX)/VOLW(L,NY,NX))
      CZ2GS(L,NY,NX)=AMAX1(0.0,Z2GS(L,NY,NX)/VOLW(L,NY,NX))
      CZ2OS(L,NY,NX)=AMAX1(0.0,Z2OS(L,NY,NX)/VOLW(L,NY,NX))
      CH2GS(L,NY,NX)=AMAX1(0.0,H2GS(L,NY,NX)/VOLW(L,NY,NX))
      ELSE
      CCO2S(L,NY,NX)=0.0
      CCH4S(L,NY,NX)=0.0
      COXYS(L,NY,NX)=0.0
      CZ2GS(L,NY,NX)=0.0
      CZ2OS(L,NY,NX)=0.0
      CH2GS(L,NY,NX)=0.0
      ENDIF
!
!     CORGC=SOC concentration
!
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CORGC(L,NY,NX)=AMIN1(0.55E+06,ORGC(L,NY,NX)/BKVL(L,NY,NX))
      ELSE
      CORGC(L,NY,NX)=0.0
      ENDIF
      end subroutine GetChemicalConcsInSoil
!------------------------------------------------------------------------------------------

      subroutine ZeroHourlyArrays(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX

      integer :: K
!     begin_execution
      FINH(L,NY,NX)=0.0
      TCO2S(L,NY,NX)=0.0
      TCO2P(L,NY,NX)=0.0
      TCOFLA(L,NY,NX)=0.0
      TCHFLA(L,NY,NX)=0.0
      TLCO2P(L,NY,NX)=0.0
      TUPOXP(L,NY,NX)=0.0
      TUPOXS(L,NY,NX)=0.0
      TUPCHS(L,NY,NX)=0.0
      TUPN2S(L,NY,NX)=0.0
      TUPN3S(L,NY,NX)=0.0
      TUPN3B(L,NY,NX)=0.0
      TUPHGS(L,NY,NX)=0.0
      TOXFLA(L,NY,NX)=0.0
      TCHFLA(L,NY,NX)=0.0
      TN2FLA(L,NY,NX)=0.0
      TNHFLA(L,NY,NX)=0.0
      THGFLA(L,NY,NX)=0.0
      TLOXYP(L,NY,NX)=0.0
      TLCH4P(L,NY,NX)=0.0
      TLN2OP(L,NY,NX)=0.0
      TLNH3P(L,NY,NX)=0.0
      TLH2GP(L,NY,NX)=0.0
      TUPNH4(L,NY,NX)=0.0
      TUPNO3(L,NY,NX)=0.0
      TUPH2P(L,NY,NX)=0.0
      TUPH1P(L,NY,NX)=0.0
      TUPNHB(L,NY,NX)=0.0
      TUPNOB(L,NY,NX)=0.0
      TUPH2B(L,NY,NX)=0.0
      TUPH1B(L,NY,NX)=0.0
      TUPNF(L,NY,NX)=0.0
      TRN4B(L,NY,NX)=0.0
      TRN3B(L,NY,NX)=0.0
      TRNOB(L,NY,NX)=0.0
      TRN2B(L,NY,NX)=0.0
      TRH1B(L,NY,NX)=0.0
      TRH2B(L,NY,NX)=0.0
      TRAL(L,NY,NX)=0.0
      TRFE(L,NY,NX)=0.0
      TRHY(L,NY,NX)=0.0
      TRCA(L,NY,NX)=0.0
      TRMG(L,NY,NX)=0.0
      TRNA(L,NY,NX)=0.0
      TRKA(L,NY,NX)=0.0
      TROH(L,NY,NX)=0.0
      TRSO4(L,NY,NX)=0.0
      TRCO3(L,NY,NX)=0.0
      TRHCO(L,NY,NX)=0.0
      TRCO2(L,NY,NX)=0.0
      TBCO2(L,NY,NX)=0.0
      TRAL1(L,NY,NX)=0.0
      TRAL2(L,NY,NX)=0.0
      TRAL3(L,NY,NX)=0.0
      TRAL4(L,NY,NX)=0.0
      TRALS(L,NY,NX)=0.0
      TRFE1(L,NY,NX)=0.0
      TRFE2(L,NY,NX)=0.0
      TRFE3(L,NY,NX)=0.0
      TRFE4(L,NY,NX)=0.0
      TRFES(L,NY,NX)=0.0
      TRCAO(L,NY,NX)=0.0
      TRCAC(L,NY,NX)=0.0
      TRCAH(L,NY,NX)=0.0
      TRCAS(L,NY,NX)=0.0
      TRMGO(L,NY,NX)=0.0
      TRMGC(L,NY,NX)=0.0
      TRMGH(L,NY,NX)=0.0
      TRMGS(L,NY,NX)=0.0
      TRNAC(L,NY,NX)=0.0
      TRNAS(L,NY,NX)=0.0
      TRKAS(L,NY,NX)=0.0
      TRH0P(L,NY,NX)=0.0
      TRH3P(L,NY,NX)=0.0
      TRC0P(L,NY,NX)=0.0
      TRF1P(L,NY,NX)=0.0
      TRF2P(L,NY,NX)=0.0
      TRC1P(L,NY,NX)=0.0
      TRC2P(L,NY,NX)=0.0
      TRM1P(L,NY,NX)=0.0
      TRH0B(L,NY,NX)=0.0
      TRH3B(L,NY,NX)=0.0
      TRF1B(L,NY,NX)=0.0
      TRF2B(L,NY,NX)=0.0
      TRC0B(L,NY,NX)=0.0
      TRC1B(L,NY,NX)=0.0
      TRC2B(L,NY,NX)=0.0
      TRM1B(L,NY,NX)=0.0
      TRXNB(L,NY,NX)=0.0
      TRXHY(L,NY,NX)=0.0
      TRXAL(L,NY,NX)=0.0
      TRXFE(L,NY,NX)=0.0
      TRXCA(L,NY,NX)=0.0
      TRXMG(L,NY,NX)=0.0
      TRXNA(L,NY,NX)=0.0
      TRXKA(L,NY,NX)=0.0
      TRXHC(L,NY,NX)=0.0
      TRXAL2(L,NY,NX)=0.0
      TRXFE2(L,NY,NX)=0.0
      TRBH0(L,NY,NX)=0.0
      TRBH1(L,NY,NX)=0.0
      TRBH2(L,NY,NX)=0.0
      TRB1P(L,NY,NX)=0.0
      TRB2P(L,NY,NX)=0.0
      TRALOH(L,NY,NX)=0.0
      TRFEOH(L,NY,NX)=0.0
      TRCACO(L,NY,NX)=0.0
      TRCASO(L,NY,NX)=0.0
      TRALPB(L,NY,NX)=0.0
      TRFEPB(L,NY,NX)=0.0
      TRCPDB(L,NY,NX)=0.0
      TRCPHB(L,NY,NX)=0.0
      TRCPMB(L,NY,NX)=0.0
      XCOFXS(L,NY,NX)=0.0
      XCHFXS(L,NY,NX)=0.0
      XOXFXS(L,NY,NX)=0.0
      XNGFXS(L,NY,NX)=0.0
      XN2FXS(L,NY,NX)=0.0
      XHGFXS(L,NY,NX)=0.0
      XN4FXW(L,NY,NX)=0.0
      XN3FXW(L,NY,NX)=0.0
      XNOFXW(L,NY,NX)=0.0
      XNXFXS(L,NY,NX)=0.0
      XH1PXS(L,NY,NX)=0.0
      XH2PXS(L,NY,NX)=0.0
      XN4FXB(L,NY,NX)=0.0
      XN3FXB(L,NY,NX)=0.0
      XNOFXB(L,NY,NX)=0.0
      XNXFXB(L,NY,NX)=0.0
      XH1BXB(L,NY,NX)=0.0
      XH2BXB(L,NY,NX)=0.0
      XALFXS(L,NY,NX)=0.0
      XFEFXS(L,NY,NX)=0.0
      XHYFXS(L,NY,NX)=0.0
      XCAFXS(L,NY,NX)=0.0
      XMGFXS(L,NY,NX)=0.0
      XNAFXS(L,NY,NX)=0.0
      XKAFXS(L,NY,NX)=0.0
      XOHFXS(L,NY,NX)=0.0
      XSOFXS(L,NY,NX)=0.0
      XCLFXS(L,NY,NX)=0.0
      XC3FXS(L,NY,NX)=0.0
      XHCFXS(L,NY,NX)=0.0
      XAL1XS(L,NY,NX)=0.0
      XAL2XS(L,NY,NX)=0.0
      XAL3XS(L,NY,NX)=0.0
      XAL4XS(L,NY,NX)=0.0
      XALSXS(L,NY,NX)=0.0
      XFE1XS(L,NY,NX)=0.0
      XFE2XS(L,NY,NX)=0.0
      XFE3XS(L,NY,NX)=0.0
      XFE4XS(L,NY,NX)=0.0
      XFESXS(L,NY,NX)=0.0
      XCAOXS(L,NY,NX)=0.0
      XCACXS(L,NY,NX)=0.0
      XCAHXS(L,NY,NX)=0.0
      XCASXS(L,NY,NX)=0.0
      XMGOXS(L,NY,NX)=0.0
      XMGCXS(L,NY,NX)=0.0
      XMGHXS(L,NY,NX)=0.0
      XMGSXS(L,NY,NX)=0.0
      XNACXS(L,NY,NX)=0.0
      XNASXS(L,NY,NX)=0.0
      XKASXS(L,NY,NX)=0.0
      XH0PXS(L,NY,NX)=0.0
      XH3PXS(L,NY,NX)=0.0
      XF1PXS(L,NY,NX)=0.0
      XF2PXS(L,NY,NX)=0.0
      XC0PXS(L,NY,NX)=0.0
      XC1PXS(L,NY,NX)=0.0
      XC2PXS(L,NY,NX)=0.0
      XM1PXS(L,NY,NX)=0.0
      XH0BXB(L,NY,NX)=0.0
      XH3BXB(L,NY,NX)=0.0
      XF1BXB(L,NY,NX)=0.0
      XF2BXB(L,NY,NX)=0.0
      XC0BXB(L,NY,NX)=0.0
      XC1BXB(L,NY,NX)=0.0
      XC2BXB(L,NY,NX)=0.0
      XM1BXB(L,NY,NX)=0.0
      DO 9955 K=0,4
      XOCFXS(K,L,NY,NX)=0.0
      XONFXS(K,L,NY,NX)=0.0
      XOPFXS(K,L,NY,NX)=0.0
      XOAFXS(K,L,NY,NX)=0.0
9955  CONTINUE
      THAW(L,NY,NX)=0.0
      THAWH(L,NY,NX)=0.0
      HTHAW(L,NY,NX)=0.0
      XCOBBL(L,NY,NX)=0.0
      XCHBBL(L,NY,NX)=0.0
      XOXBBL(L,NY,NX)=0.0
      XNGBBL(L,NY,NX)=0.0
      XN2BBL(L,NY,NX)=0.0
      XN3BBL(L,NY,NX)=0.0
      XNBBBL(L,NY,NX)=0.0
      XHGBBL(L,NY,NX)=0.0
      RTDNT(L,NY,NX)=0.0
      end subroutine ZeroHourlyArrays

end module Hour1Mod
