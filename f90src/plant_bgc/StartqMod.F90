module StartqMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  implicit none

  private
  include "parameters.h"
  include "filec.h"
  include "files.h"
  include "blkc.h"
  include "blk1cp.h"
  include "blk1cr.h"
  include "blk1g.h"
  include "blk1n.h"
  include "blk1p.h"
  include "blk1s.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk3.h"
  include "blk5.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk9a.h"
  include "blk9b.h"
  include "blk9c.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk12a.h"
  include "blk12b.h"
  include "blk14.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"

  real(r8) :: CNOPC(4),CPOPC(4)
  real(r8) :: CNOPCT,CPOPCT,CCO2A,CCO2P,COXYA,COXYP,FDM,WTSTDX

  integer :: K,L,M,NX,NY,NZ2X,NZ,N,NR,NB

  public :: startq
  contains

  SUBROUTINE startq(NHWQ,NHEQ,NVNQ,NVSQ,NZ1Q,NZ2Q)
!
!     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
!
  implicit none
  integer, intent(in) :: NHWQ,NHEQ,NVNQ,NVSQ,NZ1Q,NZ2Q
!     begin_execution
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IFLGC=PFT flag:0=not active,1=active
!     IYR0,IDAY0,IYRH,IDAYH=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,CFI=current,initial clumping factor
!     RSMH=cuticular resistance to water (h m-1)
!     RCMX=cuticular resistance to CO2 (s m-1)
!     CNWS,CPWS=protein:N,protein:P ratios
!     CWSRT=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

  DO 9995 NX=NHWQ,NHEQ
    DO 9990 NY=NVNQ,NVSQ
      NZ2X=MIN(NZ2Q,NP(NY,NX))
      DO 9985 NZ=NZ1Q,NZ2X
        IF(IFLGC(NZ,NY,NX).EQ.0)THEN

          call InitShootGrowth(NZ,NY,NX)

          call PlantLitterFractions(NZ,NY,NX)

          call PFTThermalAcclimation(NZ,NY,NX)

          call InitDimensionsandUptake(NZ,NY,NX)

          call InitPlantPhenoMorphoBio(NZ,NY,NX)

          call InitMassBalance(NZ,NY,NX)

          call InitPlantHeatandWater(NZ,NY,NX)

          call InitRootMychorMorphoBio(NZ,NY,NX)

          call InitSeedMorphoBio(NZ,NY,NX)
        !     ENDIF
        ENDIF
        ZEROP(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)
        ZEROQ(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        ZEROL(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)*1.0E+06
9985  CONTINUE
!
!     FILL OUT UNUSED ARRAYS
!
      DO 9986 NZ=NP(NY,NX)+1,5
        TCSN0(NZ,NY,NX)=0._r8
        TZSN0(NZ,NY,NX)=0._r8
        TPSN0(NZ,NY,NX)=0._r8
        TCSNC(NZ,NY,NX)=0._r8
        TZSNC(NZ,NY,NX)=0._r8
        TPSNC(NZ,NY,NX)=0._r8
        WTSTG(NZ,NY,NX)=0._r8
        WTSTGN(NZ,NY,NX)=0._r8
        WTSTGP(NZ,NY,NX)=0._r8
        DO 6401 L=1,NL(NY,NX)
          DO  K=0,1
            DO  M=1,4
              CSNC(M,K,L,NZ,NY,NX)=0._r8
              ZSNC(M,K,L,NZ,NY,NX)=0._r8
              PSNC(M,K,L,NZ,NY,NX)=0._r8
            enddo
          enddo
6401    CONTINUE
9986  CONTINUE
9990  CONTINUE
9995  CONTINUE
  RETURN
  END subroutine startq
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX

  IYR0(NZ,NY,NX)=IYRX(NZ,NY,NX)
  IDAY0(NZ,NY,NX)=IDAYX(NZ,NY,NX)
  IYRH(NZ,NY,NX)=IYRY(NZ,NY,NX)
  IDAYH(NZ,NY,NX)=IDAYY(NZ,NY,NX)
  PPI(NZ,NY,NX)=PPZ(NZ,NY,NX)
  PPX(NZ,NY,NX)=PPI(NZ,NY,NX)
  CF(NZ,NY,NX)=CFI(NZ,NY,NX)
!     WRITE(*,3232)'STARTQ',IYRC,NX,NY,NZ
!    2,IDAY0(NZ,NY,NX),IYR0(NZ,NY,NX)
!    3,IDAYH(NZ,NY,NX),IYRH(NZ,NY,NX)
!    4,IYRC,IDAYX(NZ,NY,NX),IDAYY(NZ,NY,NX)
!    5,IYRX(NZ,NY,NX),IYRY(NZ,NY,NX),IFLGC(NZ,NY,NX)
!    5,PPI(NZ,NY,NX),PPX(NZ,NY,NX),CFI(NZ,NY,NX),CF(NZ,NY,NX)
!3232  FORMAT(A8,15I8,20E12.4)
!     IF(DATAP(NZ,NY,NX).NE.'NO')THEN
  RSMH(NZ,NY,NX)=RSMX(NZ,NY,NX)/3600.0
  RCMX(NZ,NY,NX)=RSMX(NZ,NY,NX)*1.56
  CNWS(NZ,NY,NX)=2.5
  CPWS(NZ,NY,NX)=25.0
  CWSRT(NZ,NY,NX)=AMIN1(CNRT(NZ,NY,NX)*CNWS(NZ,NY,NX),CPRT(NZ,NY,NX)*CPWS(NZ,NY,NX))
  IF(ICTYP(NZ,NY,NX).EQ.3)THEN
    O2I(NZ,NY,NX)=2.10E+05
  ELSE
    O2I(NZ,NY,NX)=3.96E+05
  ENDIF
  end subroutine InitShootGrowth
!------------------------------------------------------------------------------------------

  subroutine PlantLitterFractions(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX
!
!     FRACTIONS OF PLANT LITTER ALLOCATED TO KINETIC COMPONENTS
!     PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
!
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!
!     NONSTRUCTURAL
!
  CFOPC(0,1,NZ,NY,NX)=0.0_r8
  CFOPC(0,2,NZ,NY,NX)=0.67_r8
  CFOPC(0,3,NZ,NY,NX)=0.33_r8
  CFOPC(0,4,NZ,NY,NX)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(IGTYP(NZ,NY,NX).EQ.0)THEN
    CFOPC(1,1,NZ,NY,NX)=0.07_r8
    CFOPC(1,2,NZ,NY,NX)=0.25_r8
    CFOPC(1,3,NZ,NY,NX)=0.30_r8
    CFOPC(1,4,NZ,NY,NX)=0.38_r8
    CFOPC(2,1,NZ,NY,NX)=0.07_r8
    CFOPC(2,2,NZ,NY,NX)=0.25_r8
    CFOPC(2,3,NZ,NY,NX)=0.30_r8
    CFOPC(2,4,NZ,NY,NX)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(INTYP(NZ,NY,NX).NE.0)THEN
    CFOPC(1,1,NZ,NY,NX)=0.16_r8
    CFOPC(1,2,NZ,NY,NX)=0.38_r8
    CFOPC(1,3,NZ,NY,NX)=0.34_r8
    CFOPC(1,4,NZ,NY,NX)=0.12_r8
    CFOPC(2,1,NZ,NY,NX)=0.07_r8
    CFOPC(2,2,NZ,NY,NX)=0.41_r8
    CFOPC(2,3,NZ,NY,NX)=0.37_r8
    CFOPC(2,4,NZ,NY,NX)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    CFOPC(1,1,NZ,NY,NX)=0.08_r8
    CFOPC(1,2,NZ,NY,NX)=0.41_r8
    CFOPC(1,3,NZ,NY,NX)=0.36_r8
    CFOPC(1,4,NZ,NY,NX)=0.15_r8
    CFOPC(2,1,NZ,NY,NX)=0.07_r8
    CFOPC(2,2,NZ,NY,NX)=0.41_r8
    CFOPC(2,3,NZ,NY,NX)=0.36_r8
    CFOPC(2,4,NZ,NY,NX)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.1.OR.IBTYP(NZ,NY,NX).EQ.3)THEN
    CFOPC(1,1,NZ,NY,NX)=0.07_r8
    CFOPC(1,2,NZ,NY,NX)=0.34_r8
    CFOPC(1,3,NZ,NY,NX)=0.36_r8
    CFOPC(1,4,NZ,NY,NX)=0.23_r8
    CFOPC(2,1,NZ,NY,NX)=0.0_r8
    CFOPC(2,2,NZ,NY,NX)=0.045_r8
    CFOPC(2,3,NZ,NY,NX)=0.660_r8
    CFOPC(2,4,NZ,NY,NX)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPC(1,1,NZ,NY,NX)=0.07_r8
    CFOPC(1,2,NZ,NY,NX)=0.25_r8
    CFOPC(1,3,NZ,NY,NX)=0.38_r8
    CFOPC(1,4,NZ,NY,NX)=0.30_r8
    CFOPC(2,1,NZ,NY,NX)=0.0_r8
    CFOPC(2,2,NZ,NY,NX)=0.045_r8
    CFOPC(2,3,NZ,NY,NX)=0.660_r8
    CFOPC(2,4,NZ,NY,NX)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ,NY,NX).EQ.0)THEN
    CFOPC(3,1,NZ,NY,NX)=0.07_r8
    CFOPC(3,2,NZ,NY,NX)=0.25_r8
    CFOPC(3,3,NZ,NY,NX)=0.30_r8
    CFOPC(3,4,NZ,NY,NX)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    CFOPC(3,1,NZ,NY,NX)=0.03_r8
    CFOPC(3,2,NZ,NY,NX)=0.25_r8
    CFOPC(3,3,NZ,NY,NX)=0.57_r8
    CFOPC(3,4,NZ,NY,NX)=0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    CFOPC(3,1,NZ,NY,NX)=0.0_r8
    CFOPC(3,2,NZ,NY,NX)=0.045_r8
    CFOPC(3,3,NZ,NY,NX)=0.660_r8
    CFOPC(3,4,NZ,NY,NX)=0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ,NY,NX).EQ.0)THEN
    CFOPC(4,1,NZ,NY,NX)=0.07_r8
    CFOPC(4,2,NZ,NY,NX)=0.25_r8
    CFOPC(4,3,NZ,NY,NX)=0.30_r8
    CFOPC(4,4,NZ,NY,NX)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    CFOPC(4,1,NZ,NY,NX)=0.057_r8
    CFOPC(4,2,NZ,NY,NX)=0.263_r8
    CFOPC(4,3,NZ,NY,NX)=0.542_r8
    CFOPC(4,4,NZ,NY,NX)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.1.OR.IBTYP(NZ,NY,NX).EQ.3)THEN
    CFOPC(4,1,NZ,NY,NX)=0.059_r8
    CFOPC(4,2,NZ,NY,NX)=0.308_r8
    CFOPC(4,3,NZ,NY,NX)=0.464_r8
    CFOPC(4,4,NZ,NY,NX)=0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPC(4,1,NZ,NY,NX)=0.059_r8
    CFOPC(4,2,NZ,NY,NX)=0.308_r8
    CFOPC(4,3,NZ,NY,NX)=0.464_r8
    CFOPC(4,4,NZ,NY,NX)=0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  CFOPC(5,1,NZ,NY,NX)=0.00_r8
  CFOPC(5,2,NZ,NY,NX)=0.045_r8
  CFOPC(5,3,NZ,NY,NX)=0.660_r8
  CFOPC(5,4,NZ,NY,NX)=0.295_r8
!
!     INITIALIZE C-N AND C-P RATIOS IN PLANT LITTER
!
!     CNOPC,CPOPC=fractions to allocate N,P to kinetic components
!     CFOPN,CFOPP=distribution of litter N,P to kinetic components
!
  CNOPC(1)=0.020_r8
  CNOPC(2)=0.010_r8
  CNOPC(3)=0.010_r8
  CNOPC(4)=0.020_r8
  CPOPC(1)=0.0020_r8
  CPOPC(2)=0.0010_r8
  CPOPC(3)=0.0010_r8
  CPOPC(4)=0.0020_r8
  DO 110 N=0,5
    CNOPCT=0.0_r8
    CPOPCT=0.0_r8
    DO 100 M=1,4
      CNOPCT=CNOPCT+CFOPC(N,M,NZ,NY,NX)*CNOPC(M)
      CPOPCT=CPOPCT+CFOPC(N,M,NZ,NY,NX)*CPOPC(M)
100 CONTINUE
    DO 105 M=1,4
      CFOPN(N,M,NZ,NY,NX)=CFOPC(N,M,NZ,NY,NX)*CNOPC(M)/CNOPCT
      CFOPP(N,M,NZ,NY,NX)=CFOPC(N,M,NZ,NY,NX)*CPOPC(M)/CPOPCT
105 CONTINUE
110 CONTINUE
!
!     CONCURRENT NODE GROWTH
!
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
  IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    FNOD(NZ,NY,NX)=1.0
    IF(GROUPI(NZ,NY,NX).LE.10)THEN
      NNOD(NZ,NY,NX)=3
    ELSEIF(GROUPI(NZ,NY,NX).LE.15)THEN
      NNOD(NZ,NY,NX)=4
    ELSE
      NNOD(NZ,NY,NX)=5
    ENDIF
  ELSE
    FNOD(NZ,NY,NX)=AMAX1(1.0,0.04/XRLA(NZ,NY,NX))
    NNOD(NZ,NY,NX)=24
  ENDIF
  end subroutine PlantLitterFractions
!------------------------------------------------------------------------------------------

      subroutine PFTThermalAcclimation(NZ,NY,NX)

      implicit none
      integer, intent(in) :: NZ, NY, NX
!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,ZTYPI=dynamic,initial thermal adaptation zone from PFT file
!     OFFST=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TCX=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SSTX=sensitivity to HTC (seeds oC-1 above HTC)
!
      TCZD=5.00
      TCXD=12.00
      ZTYP(NZ,NY,NX)=ZTYPI(NZ,NY,NX)
      OFFST(NZ,NY,NX)=2.667*(2.5-ZTYP(NZ,NY,NX))
      TCZ(NZ,NY,NX)=TCZD-OFFST(NZ,NY,NX)
      TCX(NZ,NY,NX)=AMIN1(15.0,TCXD-OFFST(NZ,NY,NX))
      IF(ICTYP(NZ,NY,NX).EQ.3)THEN
      IF(DATAP(NZ,NY,NX)(1:4).EQ.'soyb')THEN
      HTC(NZ,NY,NX)=30.0_r8+3.0_r8*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.002_r8
      ELSE
      HTC(NZ,NY,NX)=27.0_r8+3.0_r8*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.002_r8
      ENDIF
      ELSE
      HTC(NZ,NY,NX)=27.0_r8+3.0_r8*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.005_r8
      ENDIF
      end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

      subroutine InitDimensionsandUptake(NZ,NY,NX)

      implicit none
      integer, intent(in) :: NZ, NY, NX
!
!     SEED CHARACTERISTICS
!
!     SDVL,SDLG,SDAR=seed volume(m3),length(m),area(m2)
!     GRDM=seed C mass (g) from PFT file
!
      SDVL(NZ,NY,NX)=GRDM(NZ,NY,NX)*5.0E-06
      SDLG(NZ,NY,NX)=2.0*(0.75*SDVL(NZ,NY,NX)/PICON)**0.33
      SDAR(NZ,NY,NX)=4.0*PICON*(SDLG(NZ,NY,NX)/2.0)**2
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SDPTH=seeding depth(m) from PFT management file
!     CDPTHZ=depth to soil layer bottom from surface(m)
!     NG,NIX,NINR=seeding,upper,lower rooting layer
!     CNRTS,CPRTS=N,P root growth yield
!     RRAD1M,RRAD2M=maximum primary,secondary mycorrhizal radius (m)
!     PORT=mycorrhizal porosity
!     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake(g m-2 h-1),Km(uM), min concn (uM)
!     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
!
      SDPTH(NZ,NY,NX)=SDPTHI(NZ,NY,NX)
      DO 9795 L=NU(NY,NX),NL(NY,NX)
      IF(SDPTH(NZ,NY,NX).GE.CDPTHZ(L-1,NY,NX) &
      .AND.SDPTH(NZ,NY,NX).LT.CDPTHZ(L,NY,NX))THEN
      NG(NZ,NY,NX)=L
      NIX(NZ,NY,NX)=L
      DO 9790 NR=1,10
      NINR(NR,NZ,NY,NX)=L
9790  CONTINUE
      ENDIF
9795  CONTINUE
      CNRTS(NZ,NY,NX)=CNRT(NZ,NY,NX)*DMRT(NZ,NY,NX)
      CPRTS(NZ,NY,NX)=CPRT(NZ,NY,NX)*DMRT(NZ,NY,NX)
      RRAD1M(2,NZ,NY,NX)=5.0E-06
      RRAD2M(2,NZ,NY,NX)=5.0E-06
      PORT(2,NZ,NY,NX)=PORT(1,NZ,NY,NX)
      UPMXZH(2,NZ,NY,NX)=UPMXZH(1,NZ,NY,NX)
      UPKMZH(2,NZ,NY,NX)=UPKMZH(1,NZ,NY,NX)
      UPMNZH(2,NZ,NY,NX)=UPMNZH(1,NZ,NY,NX)
      UPMXZO(2,NZ,NY,NX)=UPMXZO(1,NZ,NY,NX)
      UPKMZO(2,NZ,NY,NX)=UPKMZO(1,NZ,NY,NX)
      UPMNZO(2,NZ,NY,NX)=UPMNZO(1,NZ,NY,NX)
      UPMXPO(2,NZ,NY,NX)=UPMXPO(1,NZ,NY,NX)
      UPKMPO(2,NZ,NY,NX)=UPKMPO(1,NZ,NY,NX)
      UPMNPO(2,NZ,NY,NX)=UPMNPO(1,NZ,NY,NX)
      RSRR(2,NZ,NY,NX)=1.0E+04
      RSRA(2,NZ,NY,NX)=1.0E+12
!
!     PORTX=tortuosity for gas transport
!     RRADP=path length for radial diffusion within root (m)
!     DMVL=volume:C ratio (m3 g-1)
!     RTLG1X,RTLG2X=specific primary,secondary root length (m g-1)
!     RTAR1X,RTAR2X=specific primary,secondary root area (m2 g-1)
!
      DO 500 N=1,2
      PORTX(N,NZ,NY,NX)=PORT(N,NZ,NY,NX)**1.33
      RRADP(N,NZ,NY,NX)=LOG(1.0/SQRT(AMAX1(0.01,PORT(N,NZ,NY,NX))))
      DMVL(N,NZ,NY,NX)=1.0E-06/(0.05*(1.0-PORT(N,NZ,NY,NX)))
      RTLG1X(N,NZ,NY,NX)=DMVL(N,NZ,NY,NX)/(3.142*RRAD1M(N,NZ,NY,NX)**2)
      RTLG2X(N,NZ,NY,NX)=DMVL(N,NZ,NY,NX)/(3.142*RRAD2M(N,NZ,NY,NX)**2)
      RRAD1X(N,NZ,NY,NX)=RRAD1M(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-PORT(N,NZ,NY,NX)))
      RRAD2X(N,NZ,NY,NX)=RRAD2M(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-PORT(N,NZ,NY,NX)))
      RTAR1X(N,NZ,NY,NX)=3.142*RRAD1X(N,NZ,NY,NX)**2
      RTAR2X(N,NZ,NY,NX)=3.142*RRAD2X(N,NZ,NY,NX)**2
500   CONTINUE
      end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

      subroutine InitPlantPhenoMorphoBio(NZ,NY,NX)

      implicit none
      integer, intent(in) :: NZ, NY, NX
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
      PP(NZ,NY,NX)=PPX(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      IFLGI(NZ,NY,NX)=0
      IDTHP(NZ,NY,NX)=0
      IDTHR(NZ,NY,NX)=0
      NBT(NZ,NY,NX)=0
      NBR(NZ,NY,NX)=0
      HTCTL(NZ,NY,NX)=0._r8
      ZC(NZ,NY,NX)=0._r8
      DO 10 NB=1,10
      IFLGA(NB,NZ,NY,NX)=0
      IFLGE(NB,NZ,NY,NX)=0
      IFLGF(NB,NZ,NY,NX)=0
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0._r8
      VSTG(NB,NZ,NY,NX)=0._r8
      VSTGX(NB,NZ,NY,NX)=0._r8
      KLEAF(NB,NZ,NY,NX)=1
      KLEAFX(NB,NZ,NY,NX)=1
      KVSTG(NB,NZ,NY,NX)=1
      KVSTGN(NB,NZ,NY,NX)=0
      GSTGI(NB,NZ,NY,NX)=0._r8
      GSTGF(NB,NZ,NY,NX)=0._r8
      TGSTGI(NB,NZ,NY,NX)=0._r8
      TGSTGF(NB,NZ,NY,NX)=0._r8
      VRNY(NB,NZ,NY,NX)=0._r8
      VRNZ(NB,NZ,NY,NX)=0._r8
      VRNS(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)
      VRNF(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)
      ATRP(NB,NZ,NY,NX)=0._r8
      FDBK(NB,NZ,NY,NX)=1.0
      FDBKX(NB,NZ,NY,NX)=1.0
      FLG4(NB,NZ,NY,NX)=0
      FLGZ(NB,NZ,NY,NX)=0
      NBTB(NB,NZ,NY,NX)=0
      IDTHB(NB,NZ,NY,NX)=1
      DO 15 M=1,10
      IDAY(M,NB,NZ,NY,NX)=0
15    CONTINUE
10    CONTINUE
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
      WSTR(NZ,NY,NX)=0._r8
      CHILL(NZ,NY,NX)=0._r8
      DO 25 NB=1,10
      CPOOL(NB,NZ,NY,NX)=0._r8
      ZPOOL(NB,NZ,NY,NX)=0._r8
      PPOOL(NB,NZ,NY,NX)=0._r8
      CPOLNB(NB,NZ,NY,NX)=0._r8
      ZPOLNB(NB,NZ,NY,NX)=0._r8
      PPOLNB(NB,NZ,NY,NX)=0._r8
      WTSHTB(NB,NZ,NY,NX)=0._r8
      WTLFB(NB,NZ,NY,NX)=0._r8
      WTNDB(NB,NZ,NY,NX)=0._r8
      WTSHEB(NB,NZ,NY,NX)=0._r8
      WTSTKB(NB,NZ,NY,NX)=0._r8
      WVSTKB(NB,NZ,NY,NX)=0._r8
      WTRSVB(NB,NZ,NY,NX)=0._r8
      WTHSKB(NB,NZ,NY,NX)=0._r8
      WTEARB(NB,NZ,NY,NX)=0._r8
      WTGRB(NB,NZ,NY,NX)=0._r8
      WTLSB(NB,NZ,NY,NX)=0._r8
      WTSHTN(NB,NZ,NY,NX)=0._r8
      WTLFBN(NB,NZ,NY,NX)=0._r8
      WTNDBN(NB,NZ,NY,NX)=0._r8
      WTSHBN(NB,NZ,NY,NX)=0._r8
      WTSTBN(NB,NZ,NY,NX)=0._r8
      WTRSBN(NB,NZ,NY,NX)=0._r8
      WTHSBN(NB,NZ,NY,NX)=0._r8
      WTEABN(NB,NZ,NY,NX)=0._r8
      WTGRBN(NB,NZ,NY,NX)=0._r8
      WTSHTP(NB,NZ,NY,NX)=0._r8
      WTLFBP(NB,NZ,NY,NX)=0._r8
      WTNDBP(NB,NZ,NY,NX)=0._r8
      WTSHBP(NB,NZ,NY,NX)=0._r8
      WTSTBP(NB,NZ,NY,NX)=0._r8
      WTRSBP(NB,NZ,NY,NX)=0._r8
      WTHSBP(NB,NZ,NY,NX)=0._r8
      WTEABP(NB,NZ,NY,NX)=0._r8
      WTGRBP(NB,NZ,NY,NX)=0._r8
      GRNXB(NB,NZ,NY,NX)=0._r8
      GRNOB(NB,NZ,NY,NX)=0._r8
      GRWTB(NB,NZ,NY,NX)=0._r8
      ARLFB(NB,NZ,NY,NX)=0._r8
      RNH3B(NB,NZ,NY,NX)=0._r8
      RCZLX(NB,NZ,NY,NX)=0._r8
      RCPLX(NB,NZ,NY,NX)=0._r8
      RCCLX(NB,NZ,NY,NX)=0._r8
      WGLFX(NB,NZ,NY,NX)=0._r8
      WGLFNX(NB,NZ,NY,NX)=0._r8
      WGLFPX(NB,NZ,NY,NX)=0._r8
      ARLFZ(NB,NZ,NY,NX)=0._r8
      RCZSX(NB,NZ,NY,NX)=0._r8
      RCPSX(NB,NZ,NY,NX)=0._r8
      RCCSX(NB,NZ,NY,NX)=0._r8
      WTSTXB(NB,NZ,NY,NX)=0._r8
      WTSTXN(NB,NZ,NY,NX)=0._r8
      WTSTXP(NB,NZ,NY,NX)=0._r8
      WGSHEX(NB,NZ,NY,NX)=0._r8
      WGSHNX(NB,NZ,NY,NX)=0._r8
      WGSHPX(NB,NZ,NY,NX)=0._r8
      HTSHEX(NB,NZ,NY,NX)=0._r8
      DO 5 L=1,NL(NY,NX)
      ARSTK(L,NB,NZ,NY,NX)=0._r8
      DO N=1,4
      SURFB(N,L,NB,NZ,NY,NX)=0._r8
      enddo
5     CONTINUE
      DO K=0,25
      ARLF(K,NB,NZ,NY,NX)=0._r8
      HTNODE(K,NB,NZ,NY,NX)=0._r8
      HTNODX(K,NB,NZ,NY,NX)=0._r8
      HTSHE(K,NB,NZ,NY,NX)=0._r8
      WGLF(K,NB,NZ,NY,NX)=0._r8
      WSLF(K,NB,NZ,NY,NX)=0._r8
      WGLFN(K,NB,NZ,NY,NX)=0._r8
      WGLFP(K,NB,NZ,NY,NX)=0._r8
      WGSHE(K,NB,NZ,NY,NX)=0._r8
      WSSHE(K,NB,NZ,NY,NX)=0._r8
      WGSHN(K,NB,NZ,NY,NX)=0._r8
      WGSHP(K,NB,NZ,NY,NX)=0._r8
      WGNODE(K,NB,NZ,NY,NX)=0._r8
      WGNODN(K,NB,NZ,NY,NX)=0._r8
      WGNODP(K,NB,NZ,NY,NX)=0._r8
      DO 55 L=1,NL(NY,NX)
      ARLFL(L,K,NB,NZ,NY,NX)=0._r8
      WGLFL(L,K,NB,NZ,NY,NX)=0._r8
      WGLFLN(L,K,NB,NZ,NY,NX)=0._r8
      WGLFLP(L,K,NB,NZ,NY,NX)=0._r8
55    CONTINUE
      IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ,NY,NX)=0._r8
      CO2B(K,NB,NZ,NY,NX)=0._r8
      HCOB(K,NB,NZ,NY,NX)=0._r8
      CPOOL4(K,NB,NZ,NY,NX)=0._r8
      DO 45 L=1,JC
      DO N=1,4
      SURF(N,L,K,NB,NZ,NY,NX)=0._r8
      enddo
45    CONTINUE
      ENDIF
      enddo
25    CONTINUE
      DO 35 L=1,NL(NY,NX)
      ARLFV(L,NZ,NY,NX)=0._r8
      WGLFV(L,NZ,NY,NX)=0._r8
      ARSTV(L,NZ,NY,NX)=0._r8
35    CONTINUE
      CPOOLP(NZ,NY,NX)=0._r8
      ZPOOLP(NZ,NY,NX)=0._r8
      PPOOLP(NZ,NY,NX)=0._r8
      CCPOLP(NZ,NY,NX)=0._r8
      CCPLNP(NZ,NY,NX)=0._r8
      CZPOLP(NZ,NY,NX)=0._r8
      CPPOLP(NZ,NY,NX)=0._r8
      WTSHT(NZ,NY,NX)=0._r8
      WTLF(NZ,NY,NX)=0._r8
      WTSHE(NZ,NY,NX)=0._r8
      WTSTK(NZ,NY,NX)=0._r8
      WVSTK(NZ,NY,NX)=0._r8
      WTRSV(NZ,NY,NX)=0._r8
      WTHSK(NZ,NY,NX)=0._r8
      WTEAR(NZ,NY,NX)=0._r8
      WTGR(NZ,NY,NX)=0._r8
      WTRT(NZ,NY,NX)=0._r8
      WTRTS(NZ,NY,NX)=0._r8
      WTND(NZ,NY,NX)=0._r8
      WTLS(NZ,NY,NX)=0._r8
      WTSHN(NZ,NY,NX)=0._r8
      WTLFN(NZ,NY,NX)=0._r8
      WTSHEN(NZ,NY,NX)=0._r8
      WTSTKN(NZ,NY,NX)=0._r8
      WTRSVN(NZ,NY,NX)=0._r8
      WTHSKN(NZ,NY,NX)=0._r8
      WTEARN(NZ,NY,NX)=0._r8
      WTGRNN(NZ,NY,NX)=0._r8
      WTNDN(NZ,NY,NX)=0._r8
      WTSHP(NZ,NY,NX)=0._r8
      WTLFP(NZ,NY,NX)=0._r8
      WTSHEP(NZ,NY,NX)=0._r8
      WTSTKP(NZ,NY,NX)=0._r8
      WTRSVP(NZ,NY,NX)=0._r8
      WTHSKP(NZ,NY,NX)=0._r8
      WTEARP(NZ,NY,NX)=0._r8
      WTGRNP(NZ,NY,NX)=0._r8
      WTNDP(NZ,NY,NX)=0._r8
      ARLFP(NZ,NY,NX)=0._r8
      WTRTA(NZ,NY,NX)=0._r8
      ARSTP(NZ,NY,NX)=0._r8
      end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

      subroutine InitMassBalance(NZ,NY,NX)

      implicit none
      integer, intent(in) :: NZ, NY, NX
!
!     INITIALIZE MASS BALANCE CHECKS
!
      IF(DATA(20).EQ.'NO'.AND.IGO.EQ.0)THEN
      CARBN(NZ,NY,NX)=0._r8
      TCSN0(NZ,NY,NX)=0._r8
      TZSN0(NZ,NY,NX)=0._r8
      TPSN0(NZ,NY,NX)=0._r8
      TCO2T(NZ,NY,NX)=0._r8
      TCO2A(NZ,NY,NX)=0._r8
      TCUPTK(NZ,NY,NX)=0._r8
      TCSNC(NZ,NY,NX)=0._r8
      TZUPTK(NZ,NY,NX)=0._r8
      TZSNC(NZ,NY,NX)=0._r8
      TPUPTK(NZ,NY,NX)=0._r8
      TPSNC(NZ,NY,NX)=0._r8
      TZUPFX(NZ,NY,NX)=0._r8
      RNH3C(NZ,NY,NX)=0._r8
      TNH3C(NZ,NY,NX)=0._r8
      VCO2F(NZ,NY,NX)=0._r8
      VCH4F(NZ,NY,NX)=0._r8
      VOXYF(NZ,NY,NX)=0._r8
      VNH3F(NZ,NY,NX)=0._r8
      VN2OF(NZ,NY,NX)=0._r8
      VPO4F(NZ,NY,NX)=0._r8
      THVSTC(NZ,NY,NX)=0._r8
      THVSTN(NZ,NY,NX)=0._r8
      THVSTP(NZ,NY,NX)=0._r8
      HVSTC(NZ,NY,NX)=0._r8
      HVSTN(NZ,NY,NX)=0._r8
      HVSTP(NZ,NY,NX)=0._r8
      RSETC(NZ,NY,NX)=0._r8
      RSETN(NZ,NY,NX)=0._r8
      RSETP(NZ,NY,NX)=0._r8
      CTRAN(NZ,NY,NX)=0._r8
      WTSTG(NZ,NY,NX)=0._r8
      WTSTGN(NZ,NY,NX)=0._r8
      WTSTGP(NZ,NY,NX)=0._r8
      WTSTDX=WTSTDI(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      DO 155 M=1,4
      WTSTDG(M,NZ,NY,NX)=WTSTDX*CFOPC(5,M,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=WTSTDX*CNSTK(NZ,NY,NX) &
      *CFOPN(5,M,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=WTSTDX*CPSTK(NZ,NY,NX) &
      *CFOPP(5,M,NZ,NY,NX)
      WTSTG(NZ,NY,NX)=WTSTG(NZ,NY,NX)+WTSTDG(M,NZ,NY,NX)
      WTSTGN(NZ,NY,NX)=WTSTGN(NZ,NY,NX)+WTSTDN(M,NZ,NY,NX)
      WTSTGP(NZ,NY,NX)=WTSTGP(NZ,NY,NX)+WTSTDP(M,NZ,NY,NX)
155   CONTINUE
      ENDIF
      end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

      subroutine InitPlantHeatandWater(NZ,NY,NX)

      implicit none
      integer, intent(in) :: NZ, NY, NX
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHCPC=canopy heat capacity (MJ m-3 K-1)
!     TCC,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSILT,PSILO,PSILG=canopy total,osmotic,turgor water potl(MPa)
!
      VHCPC(NZ,NY,NX)=cpw*WTSHT(NZ,NY,NX)*10.0E-06
      ENGYX(NZ,NY,NX)=0._r8
      DTKC(NZ,NY,NX)=0._r8
      TCC(NZ,NY,NX)=ATCA(NY,NX)
      TKC(NZ,NY,NX)=TCC(NZ,NY,NX)+TC2K
      TCG(NZ,NY,NX)=TCC(NZ,NY,NX)
      TKG(NZ,NY,NX)=TCG(NZ,NY,NX)+TC2K
      TFN3(NZ,NY,NX)=1.0
      PSILT(NZ,NY,NX)=-1.0E-03
      PSILO(NZ,NY,NX)=OSMO(NZ,NY,NX)+PSILT(NZ,NY,NX)
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      EP(NZ,NY,NX)=0._r8
      FRADP(NZ,NY,NX)=0._r8
      end subroutine InitPlantHeatandWater
!------------------------------------------------------------------------------------------

  subroutine InitRootMychorMorphoBio(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ, NY, NX
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRT,PSIRO,PSIRG=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NRT(NZ,NY,NX)=0
  UPNH4(NZ,NY,NX)=0._r8
  UPNO3(NZ,NY,NX)=0._r8
  UPH2P(NZ,NY,NX)=0._r8
  UPH1P(NZ,NY,NX)=0._r8
  UPNF(NZ,NY,NX)=0._r8
  DO 40 N=1,2
    DO 20 L=1,NL(NY,NX)
      UPWTR(N,L,NZ,NY,NX)=0._r8
      PSIRT(N,L,NZ,NY,NX)=-0.01
      PSIRO(N,L,NZ,NY,NX)=OSMO(NZ,NY,NX)+PSIRT(N,L,NZ,NY,NX)
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)-PSIRO(N,L,NZ,NY,NX))
      CPOOLR(N,L,NZ,NY,NX)=0._r8
      ZPOOLR(N,L,NZ,NY,NX)=0._r8
      PPOOLR(N,L,NZ,NY,NX)=0._r8
      CCPOLR(N,L,NZ,NY,NX)=0._r8
      CZPOLR(N,L,NZ,NY,NX)=0._r8
      CPPOLR(N,L,NZ,NY,NX)=0._r8
      CWSRTL(N,L,NZ,NY,NX)=CWSRT(NZ,NY,NX)
      WTRTL(N,L,NZ,NY,NX)=0._r8
      WTRTD(N,L,NZ,NY,NX)=0._r8
      WSRTL(N,L,NZ,NY,NX)=0._r8
      RTN1(N,L,NZ,NY,NX)=0._r8
      RTNL(N,L,NZ,NY,NX)=0._r8
      RTLGP(N,L,NZ,NY,NX)=0._r8
      RTDNP(N,L,NZ,NY,NX)=0._r8
      RTVLP(N,L,NZ,NY,NX)=0._r8
      RTVLW(N,L,NZ,NY,NX)=0._r8
      RRAD1(N,L,NZ,NY,NX)=RRAD1M(N,NZ,NY,NX)
      RRAD2(N,L,NZ,NY,NX)=RRAD2M(N,NZ,NY,NX)
      RTARP(N,L,NZ,NY,NX)=0._r8
      RTLGA(N,L,NZ,NY,NX)=1.0E-03
      RUPNH4(N,L,NZ,NY,NX)=0._r8
      RUPNO3(N,L,NZ,NY,NX)=0._r8
      RUPH2P(N,L,NZ,NY,NX)=0._r8
      RUPH1P(N,L,NZ,NY,NX)=0._r8
      RUPNHB(N,L,NZ,NY,NX)=0._r8
      RUPNOB(N,L,NZ,NY,NX)=0._r8
      RUPH2B(N,L,NZ,NY,NX)=0._r8
      RUPH1B(N,L,NZ,NY,NX)=0._r8
      ROXYP(N,L,NZ,NY,NX)=0._r8
      RUNNHP(N,L,NZ,NY,NX)=0._r8
      RUNNBP(N,L,NZ,NY,NX)=0._r8
      RUNNOP(N,L,NZ,NY,NX)=0._r8
      RUNNXP(N,L,NZ,NY,NX)=0._r8
      RUPP2P(N,L,NZ,NY,NX)=0._r8
      RUPP1P(N,L,NZ,NY,NX)=0._r8
      RUPP2B(N,L,NZ,NY,NX)=0._r8
      RUPP1B(N,L,NZ,NY,NX)=0._r8
      CCO2A=CCO2EI(NY,NX)
      CCO2P=0.030*EXP(-2.621-0.0317*ATCA(NY,NX))*CO2EI(NY,NX)
      CO2A(N,L,NZ,NY,NX)=CCO2A*RTVLP(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=CCO2P*RTVLW(N,L,NZ,NY,NX)
      RCOFLA(N,L,NZ,NY,NX)=0._r8
      RCODFA(N,L,NZ,NY,NX)=0._r8
      RCO2S(N,L,NZ,NY,NX)=0._r8
      RCO2P(N,L,NZ,NY,NX)=0._r8
      COXYA=COXYE(NY,NX)
      COXYP=0.032*EXP(-6.175-0.0211*ATCA(NY,NX))*OXYE(NY,NX)
      OXYA(N,L,NZ,NY,NX)=COXYA*RTVLP(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=COXYP*RTVLW(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=0._r8
      CH4P(N,L,NZ,NY,NX)=0._r8
      Z2OA(N,L,NZ,NY,NX)=0._r8
      Z2OP(N,L,NZ,NY,NX)=0._r8
      ZH3A(N,L,NZ,NY,NX)=0._r8
      ZH3P(N,L,NZ,NY,NX)=0._r8
      H2GA(N,L,NZ,NY,NX)=0._r8
      H2GP(N,L,NZ,NY,NX)=0._r8
      WFR(N,L,NZ,NY,NX)=1.0
      DO 30 NR=1,10
        RTN2(N,L,NR,NZ,NY,NX)=0._r8
        RTLG1(N,L,NR,NZ,NY,NX)=0._r8
        WTRT1(N,L,NR,NZ,NY,NX)=0._r8
        WTRT1N(N,L,NR,NZ,NY,NX)=0._r8
        WTRT1P(N,L,NR,NZ,NY,NX)=0._r8
        RTLG2(N,L,NR,NZ,NY,NX)=0._r8
        WTRT2(N,L,NR,NZ,NY,NX)=0._r8
        WTRT2N(N,L,NR,NZ,NY,NX)=0._r8
        WTRT2P(N,L,NR,NZ,NY,NX)=0._r8
        RTDP1(N,NR,NZ,NY,NX)=SDPTH(NZ,NY,NX)
        RTWT1(N,NR,NZ,NY,NX)=0._r8
        RTWT1N(N,NR,NZ,NY,NX)=0._r8
        RTWT1P(N,NR,NZ,NY,NX)=0._r8
30    CONTINUE
      IF(N.EQ.1)THEN
        DO 6400 K=0,1
          DO  M=1,4
            CSNC(M,K,L,NZ,NY,NX)=0._r8
            ZSNC(M,K,L,NZ,NY,NX)=0._r8
            PSNC(M,K,L,NZ,NY,NX)=0._r8
          enddo
6400    CONTINUE
        CPOOLN(L,NZ,NY,NX)=0._r8
        ZPOOLN(L,NZ,NY,NX)=0._r8
        PPOOLN(L,NZ,NY,NX)=0._r8
        WTNDL(L,NZ,NY,NX)=0._r8
        WTNDLN(L,NZ,NY,NX)=0._r8
        WTNDLP(L,NZ,NY,NX)=0._r8
        RUPNF(L,NZ,NY,NX)=0._r8
      ENDIF
20  CONTINUE
40  CONTINUE

  RUPNH4(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RUPNHB(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RUPH2P(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RUPH2B(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RTDNP(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     WTLSB=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     VOLWP,VOLWC=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     WTRTL,WTRTD=total root C mass (g)
!     WSRTL=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  WTRVX(NZ,NY,NX)=GRDM(NZ,NY,NX)*PP(NZ,NY,NX)
  WTRVC(NZ,NY,NX)=WTRVX(NZ,NY,NX)
  WTRVN(NZ,NY,NX)=CNGR(NZ,NY,NX)*WTRVC(NZ,NY,NX)
  WTRVP(NZ,NY,NX)=CPGR(NZ,NY,NX)*WTRVC(NZ,NY,NX)
  WTLFBN(1,NZ,NY,NX)=CNGR(NZ,NY,NX)*WTLFB(1,NZ,NY,NX)
  WTLFBP(1,NZ,NY,NX)=CPGR(NZ,NY,NX)*WTLFB(1,NZ,NY,NX)
  WTLSB(1,NZ,NY,NX)=WTLFB(1,NZ,NY,NX)+WTSHEB(1,NZ,NY,NX)
  WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(1,NZ,NY,NX)
  FDM=AMIN1(1.0,0.16-0.045*PSILT(NZ,NY,NX))
  VOLWP(NZ,NY,NX)=1.0E-06*WTLS(NZ,NY,NX)/FDM
  VOLWC(NZ,NY,NX)=0._r8
  ZPOOL(1,NZ,NY,NX)=CNGR(NZ,NY,NX)*CPOOL(1,NZ,NY,NX)
  PPOOL(1,NZ,NY,NX)=CPGR(NZ,NY,NX)*CPOOL(1,NZ,NY,NX)
  WTRT1N(1,NG(NZ,NY,NX),1,NZ,NY,NX)=CNGR(NZ,NY,NX)*WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
  WTRT1P(1,NG(NZ,NY,NX),1,NZ,NY,NX)=CPGR(NZ,NY,NX)*WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
  RTWT1N(1,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*RTWT1(1,1,NZ,NY,NX)
  RTWT1P(1,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*RTWT1(1,1,NZ,NY,NX)
  WTRTL(1,NG(NZ,NY,NX),NZ,NY,NX)=WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
  WTRTD(1,NG(NZ,NY,NX),NZ,NY,NX)=WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
  WSRTL(1,NG(NZ,NY,NX),NZ,NY,NX)=WTRTL(1,NG(NZ,NY,NX),NZ,NY,NX)*CWSRT(NZ,NY,NX)
  ZPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=CNGR(NZ,NY,NX)*CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
  PPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=CPGR(NZ,NY,NX)*CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
  end subroutine InitSeedMorphoBio

  end module StartqMod
