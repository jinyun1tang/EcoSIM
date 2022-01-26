
      module StomateMod
      use data_kind_mod, only : r8 => SHR_KIND_R8
      implicit none

      private
      include "parameters.h"
      include "blkc.h"
      include "blk1cp.h"
      include "blk1g.h"
      include "blk1n.h"
      include "blk1p.h"
      include "blk2a.h"
      include "blk3.h"
      include "blk5.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk1u.h"

      real(r8) :: ACTV,CH2O,CC4M,CCBS,ETDN4,ETLF4,EGRO4,ETDN,ETLF
      real(r8) :: EGRO,PARX,PARJ,RI,RAC,RTK,RSX,STK,TCCZ,TKCO,TFN1
      real(r8) :: TFN2,TFNE,VCDN4,VL,VCDN,VOGRO,WSDN,XKO2L

      real(r8) :: FLG4Y(0:5)
!
!     QNTM=quantum efficiency (umol e- umol-1 PAR)
!     CURV=shape parameter for e- transport response to PAR
!     ELEC3,ELEC4=e- requirement for CO2 fixn by rubisco,PEP carboxylase
!     (umol e- umol CO2)
!     CNKI,CPKI=nonstruct N,P inhibition constant on rubisco (g N,P g-1 C)
!     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
!     ATRPZ=hours to full dehardening of conifers in spring (h)
!     COMP4=C4 CO2 compensation point (uM)
!     FDML=leaf water content (g H2O g-1 C)
!     FBS,FMP=leaf water content in bundle sheath, mesophyll in C4 CO2 fixn
!     C4KI=nonstructural C inhibition constant on PEP carboxylase (uM)
!     FLG4Y=number of hours with no grain fill to terminate annuals
!
      real(r8), PARAMETER :: QNTM=0.45,CURV=0.70,CURV2=2.0*CURV &
      ,CURV4=4.0*CURV,ELEC3=4.5,ELEC4=3.0
      real(r8), PARAMETER :: CNKI=1.0E+02,CPKI=1.0E+03
      real(r8), PARAMETER :: RSMY=2.78E-03,ATRPZ=276.9
      real(r8), PARAMETER :: COMP4=0.5,FDML=6.0,FBS=0.2*FDML &
      ,FMP=0.8*FDML,C4KI=5.0E+06
      DATA FLG4Y/336.0,672.0,672.0,672.0,672.0,672.0/

      public :: stomate
      contains

      Subroutine stomate(I,J,NZ,NY,NX)
!
!     THIS SUBROUTINE CALCULATES CANOPY STOMATAL RESISTANCE AT MAXIMUM
!     CANOPY TURGOR FOR USE IN ENERGY BALANCE EQUATIONS IN 'UPTAKE'
!
      implicit none
      integer, intent(in) :: I, J
      integer, intent(in) :: NZ,NY,NX

      integer :: K,L,M,NB,N
!     begin_execution
!
!     CANOPY TEMPERATURE + OFFSET FOR THERMAL ADAPTATION FROM 'READQ'
!
!     CANOPY BOUNDARY LAYER RESISTANCE
!
!     RI=Richardson's number
!     RIB=canopy isothermal Richardson�s number
!     TKA,TKCZ=air,canopy temperature
!     RAZ=canopy isothermal boundary later resistance
!     RAC=canopy boundary layer resistance to CO2
!     FMOL=number of moles of air per m3
!
      RI=AMAX1(-0.3,AMIN1(0.075,RIB(NY,NX)*(TKA(NY,NX)-TKCZ(NZ,NY,NX))))
      RAC=1.34*AMAX1(5.56E-03,RAZ(NZ,NY,NX)/(1.0-10.0*RI))
      FMOL(NZ,NY,NX)=1.2194E+04/TKCZ(NZ,NY,NX)
!
!     CANOPY CO2 CONCENTRATION FROM CO2 INFLUXES AND EFFLUXES
!
!     CO2Q,CO2E=CO2 concentrations in canopy air,atmosphere, umol mol-1 (ppmv)
!     CNETX=net CO2 flux in canopy air from soil,plants, g d-2 h-1
!
      CO2Q(NZ,NY,NX)=CO2E(NY,NX)-8.33E+04*CNETX(NY,NX) &
      *RAC/FMOL(NZ,NY,NX)
      CO2Q(NZ,NY,NX)=AMIN1(CO2E(NY,NX)+200.0 &
      ,AMAX1(0.0,CO2E(NY,NX)-200.0,CO2Q(NZ,NY,NX)))
!
!     MESOPHYLL CO2 CONCENTRATION FROM CI:CA RATIO ENTERED IN 'READQ'
!
!     CO2I=intercellular CO2 concentration
!     FCO2=intercellular:atmospheric CO2 concn ratio from PFT file, parameter
!     SSIN=sine of solar angle
!     ARLFP=PFT leaf area
!
      CO2I(NZ,NY,NX)=FCO2(NZ,NY,NX)*CO2Q(NZ,NY,NX)

      IF(SSIN(NY,NX).GT.0.0.AND.ARLFP(NZ,NY,NX) &
      .GT.ZEROP(NZ,NY,NX))THEN
!
      call PhotoActivePFT(NZ,NY,NX)
      ELSE
!
      RSMN(NZ,NY,NX)=RSMH(NZ,NY,NX)
      ENDIF
!     IF(ICTYP(NZ,NY,NX).EQ.3)THEN
!     WRITE(19,3010)'CH2O',I,J,CH2O
!     ELSEIF(ICTYP(NZ,NY,NX).EQ.4)THEN
!     WRITE(20,3010)'CH2O',I,J,CH2O
!     ENDIF
3010  FORMAT(A8,2I4,1E12.4)
      RETURN
      END subroutine stomate
!------------------------------------------------------------------------------------------

      subroutine C3ShadedLeaves(K,N,M,L,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: K,N,M,L,NB,NZ,NY,NX
!     begin_execution
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PARDIF=diffuse PAR flux
!     ETGR=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
      PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     VGRO=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     CH2O=total rubisco carboxylation rate
!     FDBK=N,P feedback inhibition on C3 CO2 fixation
!     SURFX=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAU0(L+1,NY,NX)
      end subroutine C3ShadedLeaves
!------------------------------------------------------------------------------------------

      subroutine C3SunlitLeaves(K,N,M,L,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: K,N,M,L,NB,NZ,NY,NX
!     begin_execution
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     ETGRO=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
      PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     VGRO=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     FDBK=N,P feedback inhibition on C3 CO2 fixation
!     CH2O=total rubisco carboxylation rate
!     SURFX=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAUS(L+1,NY,NX)
!
!     IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1.AND.K.EQ.KLEAF(NB,NZ,NY,NX)-1
!    2.AND.J.EQ.14)THEN
!     WRITE(20,6798)'STD',I,J,L,M,N,K,NB,VL,PAR(N,M,L,NZ,NY,NX),RAPS
!    2,TKCZ(NZ,NY,NX),CO2Q(NZ,NY,NX),ETGRO(K,NB,NZ,NY,NX)
!    3,CBXN(K,NB,NZ,NY,NX),VGRO(K,NB,NZ,NY,NX),EGRO
!    3,FDBK(NB,NZ,NY,NX),CH2O,TFN1,TFN2,TFNE,WSDN
!    3,VCGRO(K,NB,NZ,NY,NX),VCDN,CO2I(NZ,NY,NX),CO2L(NZ,NY,NX)
6798  FORMAT(A8,7I4,40E12.4)
!     ENDIF
      end subroutine C3SunlitLeaves
!------------------------------------------------------------------------------------------

      subroutine C3PhotosynsCanopyLayerL(L,K,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: L,K,NB,NZ,NY,NX

      integer :: N,M
!     begin_execution
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
      DO 3600 N=1,4
      DO 3500 M=1,4
      IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!     SUNLIT LEAVES
!
      IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
      call C3SunlitLeaves(K,N,M,L,NB,NZ,NY,NX)
      ENDIF
!
!     SHADED LEAVES
      IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
      call C3ShadedLeaves(K,N,M,L,NB,NZ,NY,NX)
      ENDIF
!
      ENDIF
3500  CONTINUE
3600  CONTINUE
      end subroutine C3PhotosynsCanopyLayerL
!------------------------------------------------------------------------------------------

      subroutine C3Photosynthesis(K,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: K,NB,NZ,NY,NX

      integer :: L
!     begin_execution
!
!     SURFICIAL DENSITY OF RUBISCO AND ITS CHLOROPHYLL
!
!     VCDN=surficial density of rubisco in mesophyll
!     ETDN=surficial density of chlorophyll in esophyll
!     RUBP=fraction of leaf protein in rubisco
!     CHL=fraction of leaf protein in mesophyll chlorophyll
!     WSDN=leaf protein surficial density
!
      VCDN=RUBP(NZ,NY,NX)*WSDN
      ETDN=CHL(NZ,NY,NX)*WSDN
!
!     CO2-LIMITED C3 CARBOXYLATION RATES
!
!     VCGRO=rubisco carboxylation rate unlimited by CO2
!     VCMX=specific rubisco carboxylation activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN=surficial density of rubisco in mesophyll
!     VOGRO=rubisco oxygenation rate
!     TFN2=temperature function for oxygenation
!     COMPL=C3 CO2 compensation point (uM)
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!     VGRO=rubisco carboxylation rate limited by CO2
!
      VCGRO(K,NB,NZ,NY,NX)=VCMX(NZ,NY,NX)*TFN1*VCDN
      VOGRO=VOMX(NZ,NY,NX)*TFN2*VCDN
      COMPL(K,NB,NZ,NY,NX)=0.5*O2L(NZ,NY,NX)*VOGRO*XKCO2L(NZ,NY,NX) &
      /(VCGRO(K,NB,NZ,NY,NX)*XKO2L)
      VGRO(K,NB,NZ,NY,NX)=AMAX1(0.0,VCGRO(K,NB,NZ,NY,NX) &
      *(CO2L(NZ,NY,NX)-COMPL(K,NB,NZ,NY,NX)) &
      /(CO2L(NZ,NY,NX)+XKCO2O(NZ,NY,NX)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     ETGRO=light-limited rubisco carboxylation rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN=surficial density of chlorophyll in mesophyll
!     CBXN=rubisco caboxylation efficiency
!     CO2L=intercellular CO2 concentrations (uM)
!     COMPL=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
      ETGRO(K,NB,NZ,NY,NX)=ETMX(NZ,NY,NX)*TFNE*ETDN
      CBXN(K,NB,NZ,NY,NX)=AMAX1(0.0,(CO2L(NZ,NY,NX) &
      -COMPL(K,NB,NZ,NY,NX))/(ELEC3*CO2L(NZ,NY,NX) &
      +10.5*COMPL(K,NB,NZ,NY,NX)))
!
!     FOR EACH CANOPY LAYER
!
!     ARLFL=leaf area
!     SURFX=unself-shaded leaf surface area
!
      DO 3700 L=JC,1,-1
      IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
      call C3PhotosynsCanopyLayerL(L,K,NB,NZ,NY,NX)
!
      ENDIF
3700  CONTINUE
      end subroutine C3Photosynthesis
!------------------------------------------------------------------------------------------

      subroutine C4Photosynthesis(K,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: K,NB,NZ,NY,NX

      integer :: L
!     begin_execution
!
!     FEEDBACK ON C4 CARBOXYLATION FROM C4 NON-STRUCTURAL C
!
!     CC4M,CCBS=C4 nonstruct C concn in mesophyll,bundle sheath (uM)
!     CPOOL4,CO2B=C4 nonstructural C mass in mesophyll,bundle sheath
!     WGLF=leaf C mass
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
      CC4M=AMAX1(0.0,0.021E+09*CPOOL4(K,NB,NZ,NY,NX) &
      /(WGLF(K,NB,NZ,NY,NX)*FMP))
      CCBS=AMAX1(0.0,0.083E+09*CO2B(K,NB,NZ,NY,NX) &
      /(WGLF(K,NB,NZ,NY,NX)*FBS))
      FDBK4(K,NB,NZ,NY,NX)=1.0/(1.0+CC4M/C4KI)
      FDBK4(K,NB,NZ,NY,NX)=FDBK4(K,NB,NZ,NY,NX)*FDBKX(NB,NZ,NY,NX)
!
!     SURFICIAL DENSITY OF PEPC AND ITS CHLOROPHYLL
!
!     VCDN4=surficial density of PEP carboxylase in mesophyll
!     ETDN4=surficial density of chlorophyll in mesophyll
!     PEPC=fraction of leaf protein in PEP carboxylase
!     CHL4=fraction of leaf protein in mesophyll chlorophyll
!     WSDN=leaf protein surficial density
!
      VCDN4=PEPC(NZ,NY,NX)*WSDN
      ETDN4=CHL4(NZ,NY,NX)*WSDN
!
!     CO2-LIMITED C4 CARBOXYLATION RATES
!
!     VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
!     VCMX4=specific PEP carboxylase activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN4=surficial density of PEP carboxylase in mesophyll
!     CO2L=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     XKCO24=Km for VCMX4 from PFT file (uM)
!
      VCGR4(K,NB,NZ,NY,NX)=VCMX4(NZ,NY,NX)*TFN1*VCDN4
      VGRO4(K,NB,NZ,NY,NX)=AMAX1(0.0,VCGR4(K,NB,NZ,NY,NX) &
      *(CO2L(NZ,NY,NX)-COMP4)/(CO2L(NZ,NY,NX)+XKCO24(NZ,NY,NX)))
!
!     C4 ELECTRON TRANSFER RATES
!
!     ETGR4=light saturated e- transport rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN4=surficial density of chlorophyll in mesophyll
!     CBXN4=PEP caboxylation efficiency
!     CO2L=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!
      ETGR4(K,NB,NZ,NY,NX)=ETMX(NZ,NY,NX)*TFNE*ETDN4
      CBXN4(K,NB,NZ,NY,NX)=AMAX1(0.0,(CO2L(NZ,NY,NX)-COMP4) &
      /(ELEC4*CO2L(NZ,NY,NX)+10.5*COMP4))
!
!     FOR EACH CANOPY LAYER
!
!     ARLFL=leaf area
!     SURFX=unself-shaded leaf surface area
!
      DO 2700 L=JC,1,-1
      IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
      call C4PhotosynsCanopyLayerL(L,K,NB,NZ,NY,NX)
      ENDIF
2700  CONTINUE
!
!     VARIABLES FOR C3 PHOTOSYNTHESIS DRIVEN BY C4
!
!     VCDN=surficial density of rubisco in bundle sheath
!     ETDN=surficial density of chlorophyll in bundle sheath
!     RUBP=fraction of leaf protein in rubisco
!     CHL=fraction of leaf protein in bundle sheath chlorophyll
!     WSDN=leaf protein surficial density
!
      VCDN=RUBP(NZ,NY,NX)*WSDN
      ETDN=CHL(NZ,NY,NX)*WSDN
!
!     CO2-LIMITED C3 CARBOXYLATION RATES
!
!     VCGRO=rubisco carboxylation rate unlimited by CO2
!     VCMX=specific rubisco carboxylation activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN=surficial density of rubisco in bundle sheath
!     VOGRO=rubisco oxygenation rate
!     TFN2=temperature function for oxygenation
!     COMPL=C3 CO2 compensation point (uM)
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!     VGRO=rubisco carboxylation rate limited by CO2
!     CCBS=C4 nonstruct C concn in bundle sheath (uM)
!
      VCGRO(K,NB,NZ,NY,NX)=VCMX(NZ,NY,NX)*TFN1*VCDN
      VOGRO=VOMX(NZ,NY,NX)*TFN2*VCDN
      COMPL(K,NB,NZ,NY,NX)=0.5*O2L(NZ,NY,NX)*VOGRO*XKCO2L(NZ,NY,NX) &
      /(VCGRO(K,NB,NZ,NY,NX)*XKO2L)
      VGRO(K,NB,NZ,NY,NX)=AMAX1(0.0,VCGRO(K,NB,NZ,NY,NX) &
      *(CCBS-COMPL(K,NB,NZ,NY,NX))/(CCBS+XKCO2O(NZ,NY,NX)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     ETGRO=light-limited rubisco carboxylation rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN=surficial density of chlorophyll in bundle sheath
!     CBXN=rubisco caboxylation efficiency
!     CO2L=intercellular CO2 concentrations (uM)
!     COMPL=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
      ETGRO(K,NB,NZ,NY,NX)=ETMX(NZ,NY,NX)*TFNE*ETDN
      CBXN(K,NB,NZ,NY,NX)=AMAX1(0.0,(CCBS-COMPL(K,NB,NZ,NY,NX)) &
      /(ELEC3*CCBS+10.5*COMPL(K,NB,NZ,NY,NX)))

      end subroutine C4Photosynthesis
!------------------------------------------------------------------------------------------

      subroutine C4PhotosynsCanopyLayerL(L,K,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: L,K,NB,NZ,NY,NX

      integer :: M,N
!     begin_execution
!
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
      DO 2600 N=1,4
      DO 2500 M=1,4
      IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!     SUNLIT LEAVES
      IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
      call C4SunlitLeaves(K,N,M,L,NB,NZ,NY,NX)
      ENDIF
!
!     SHADED LEAVES
      IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
      call C4ShadedLeaves(K,N,M,L,NB,NZ,NY,NX)
      ENDIF
      ENDIF
2500  CONTINUE
2600  CONTINUE

      end subroutine C4PhotosynsCanopyLayerL
!------------------------------------------------------------------------------------------

      subroutine C4ShadedLeaves(K,N,M,L,NB,NZ,NY,NX)

      implicit none
      integer, intent(in) :: K,N,M,L,NB,NZ,NY,NX
!     begin_execution
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PARDIF=diffuse PAR flux
!     ETGR4=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
      PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
      ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
      EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     VGRO4=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     CH2O=total PEP carboxylation rate
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!     SURFX=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
      VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAU0(L+1,NY,NX)
!     WRITE(*,6799)'STB',I,J,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX),RAPS
!    2,TKCZ(NZ,NY,NX),CO2Q(NZ,NY,NX),ETGR4(K,NB,NZ,NY,NX)
!    3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO4
!    3,FDBK4(K,NB,NZ,NY,NX),CH2O,VGRO4(K,NB,NZ,NY,NX),EGRO4
!    3,VCGR4(K,NB,NZ,NY,NX),CO2I(NZ,NY,NX),CO2L(NZ,NY,NX)
6799  FORMAT(A8,6I4,40E12.4)
      end subroutine C4ShadedLeaves
!------------------------------------------------------------------------------------------

      subroutine C4SunlitLeaves(K,N,M,L,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: K,N,M,L,NB,NZ,NY,NX
!     begin_execution
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     ETGR4=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
      PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
      ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
      EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     VGRO4=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!     CH2O=total PEP carboxylation rate
!     SURFX=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
      VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAUS(L+1,NY,NX)
!     IF(L.GT.NC-4.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.3)THEN
!     WRITE(*,6789)'STO',I,J,L,M,N,K,L,VL
!    2,PAR(N,M,L,NZ,NY,NX),RAPS
!    2,TKCZ(NZ,NY,NX),CO2Q(NZ,NY,NX),ETGR4(K,NB,NZ,NY,NX)
!    3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO4
!    3,FDBK4(K,NB,NZ,NY,NX),CH2O,VGRO4(K,NB,NZ,NY,NX),EGRO4
!    3,SURFX(N,L,K,NB,NZ,NY,NX)
!    3,VCGR4(K,NB,NZ,NY,NX),CO2I(NZ,NY,NX),CO2L(NZ,NY,NX),TFN1,TFN2
!    4,TFNE,WSDN,VCDN4
6789  FORMAT(A8,7I4,40E12.4)
!     ENDIF
      end subroutine C4SunlitLeaves
!------------------------------------------------------------------------------------------

      subroutine LivingBranch(NB,NZ,NY,NX)
      implicit none
      integer, intent(in):: NB,NZ,NY,NX

      integer :: K
!     begin_execution

      DO 2800 K=1,25
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WSDN=WSLF(K,NB,NZ,NY,NX)/ARLF(K,NB,NZ,NY,NX)
      ELSE
      WSDN=0.0
      ENDIF
!     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
!     WRITE(*,2125)'WSDN',I,J,NX,NY,NZ,NB,K,WSDN
!    2,WGLF(K,NB,NZ,NY,NX),WSLF(K,NB,NZ,NY,NX)
!    3,ARLF(K,NB,NZ,NY,NX)
2125  FORMAT(A8,7I4,12E12.4)
!     ENDIF

      IF(WSDN.GT.ZERO)THEN
!
!     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
!     C4 PHOTOSYNTHESIS
      call C4Photosynthesis(K,NB,NZ,NY,NX)
      ELSE
!     C3 PHOTOSYNTHESIS
      call C3Photosynthesis(K,NB,NZ,NY,NX)
      ENDIF
!
      ELSE
      VCGR4(K,NB,NZ,NY,NX)=0.0
      VCGRO(K,NB,NZ,NY,NX)=0.0
      ENDIF
2800  CONTINUE
      end subroutine LivingBranch
!------------------------------------------------------------------------------------------

      subroutine PhenoActiveBranch(NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: NB,NZ,NY,NX
!     begin_execution
!
!     FEEDBACK ON C3 CARBOXYLATION FROM NON-STRUCTURAL C:N:P
!
!     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
!     FDBK=N,P feedback inhibition on C3 CO2 fixation
!     CNKI,CPKI=nonstructural N,P inhibition constant on rubisco
!
      IF(CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
      FDBK(NB,NZ,NY,NX)=AMIN1(CZPOLB(NB,NZ,NY,NX) &
      /(CZPOLB(NB,NZ,NY,NX)+CCPOLB(NB,NZ,NY,NX)/CNKI) &
      ,CPPOLB(NB,NZ,NY,NX) &
      /(CPPOLB(NB,NZ,NY,NX)+CCPOLB(NB,NZ,NY,NX)/CPKI))
      ELSE
      FDBK(NB,NZ,NY,NX)=1.0
      ENDIF
!
!     CHILLING
!
!     CHILL=accumulated chilling hours used to limit CO2 fixn
!
!     FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)/(1.0+0.25*CHILL(NZ,NY,NX))
!
!     DEHARDENING OF EVERGREENS IN SPRING
!
!     ATRP=hours above threshold temperature for dehardening since leafout
!     ATRPZ=hours to full dehardening of conifers in spring
!
      IF(IWTYP(NZ,NY,NX).NE.0.AND.IBTYP(NZ,NY,NX).GE.2)THEN
      FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)*AMAX1(0.0,AMIN1(1.0 &
      ,ATRP(NB,NZ,NY,NX)/(0.9*ATRPZ)))
      ENDIF
!
!     TERMINATION OF ANNUALS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     FLG4=number of hours with no grain fill after start of grain fill
!     FLG4Y=number of hours with no grain fill to terminate annuals
!
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.FLG4(NB,NZ,NY,NX).GT.0.0)THEN
      FDBKX(NB,NZ,NY,NX)=AMAX1(0.0 &
      ,1.0-FLG4(NB,NZ,NY,NX)/FLG4Y(IWTYP(NZ,NY,NX)))
      ELSE
      FDBKX(NB,NZ,NY,NX)=1.0
      ENDIF
      FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)*FDBKX(NB,NZ,NY,NX)
!     IF(NZ.EQ.4)THEN
!     WRITE(*,4242)'FDBK',I,J,NZ,NB,IDTHB(NB,NZ,NY,NX)
!    2,FDBK(NB,NZ,NY,NX),VRNS(NB,NZ,NY,NX),VRNF(NB,NZ,NY,NX)
!    3,CCPOLB(NB,NZ,NY,NX),CZPOLB(NB,NZ,NY,NX),CPPOLB(NB,NZ,NY,NX)
!    3,FDBKX(NB,NZ,NY,NX),ATRP(NB,NZ,NY,NX)
4242  FORMAT(A8,5I4,12E20.4)
!     ENDIF
!
!     FOR EACH NODE
!
!     IDTHB=branch life flag:0=living,1=dead
!     ARLF,WGLF,WSLF=leaf area,C mass,protein mass
!     WSDN=leaf protein surficial density
!
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      call LivingBranch(NB,NZ,NY,NX)
      ENDIF
      end subroutine PhenoActiveBranch
!------------------------------------------------------------------------------------------

      subroutine PrepPhotosynthesis(NZ,NY,NX)
      implicit none
      integer, intent(in) :: NZ,NY,NX
!     begin_execution
!
!     CO2 AND O2 AQUEOUS SOLUBILITY
!
!     TCCZ=canopy temperature
!     SCO2,SO2=solubility of CO2,O2 (uM/(umol mol-1))
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!
      TCCZ=TKCZ(NZ,NY,NX)-273.15
      SCO2(NZ,NY,NX)=EXP(-2.621-0.0317*TCCZ)
      SO2(NZ,NY,NX)=EXP(-6.175-0.0211*TCCZ)
      CO2L(NZ,NY,NX)=CO2I(NZ,NY,NX)*SCO2(NZ,NY,NX)
      O2L(NZ,NY,NX)=O2I(NZ,NY,NX)*SO2(NZ,NY,NX)
      DCO2(NZ,NY,NX)=FMOL(NZ,NY,NX)*(CO2Q(NZ,NY,NX)-CO2I(NZ,NY,NX))
!
!     ARRHENIUS FUNCTIONS FOR CARBOXYLATION AND OXYGENATION
!
!     TKCZ,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN1,TFN2,TFNE=temperature function for carboxylation,
!     oxygenation,e- transport (25 oC =1)
!     8.313,710.0=gas constant,enthalpy
!     197500,222500=energy of high,low temp inactivn(KJ mol-1)
!     65000,60000,43000=activation energy for carboxylation,
!     oxygenation,e- transport
!
      CH2O=0.0
      TKCO=TKCZ(NZ,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKCO
      STK=710.0*TKCO
      ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFN1=EXP(26.237-65000/RTK)/ACTV
      TFN2=EXP(24.220-60000/RTK)/ACTV
      TFNE=EXP(17.362-43000/RTK)/ACTV
!
!     M-M CONSTANT FOR CARBOXYLATION FROM 'READQ' ADJUSTED FOR TEMPERATURE
!
!     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!
      XKCO2L(NZ,NY,NX)=XKCO2(NZ,NY,NX)*EXP(16.136-40000/RTK)
      XKO2L=XKO2(NZ,NY,NX)*EXP(8.067-20000/RTK)
      XKCO2O(NZ,NY,NX)=XKCO2L(NZ,NY,NX)*(1.0+O2L(NZ,NY,NX)/XKO2L)
      end subroutine PrepPhotosynthesis
!------------------------------------------------------------------------------------------

      subroutine PhotoActivePFT(NZ,NY,NX)
      implicit none
      integer, intent(in) :: NZ,NY,NX

      integer :: NB,K
!     begin_execution

      call PrepPhotosynthesis(NZ,NY,NX)
!
!     FOR EACH BRANCH
!
      DO 2900 NB=1,NBR(NZ,NY,NX)
!
!     FEEDBACK ON CO2 FIXATION
!
!     IWTYP=phenology type from PFT file
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!
      IF(IWTYP(NZ,NY,NX).EQ.0 &
      .OR.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX) &
      .OR.VRNF(NB,NZ,NY,NX).LT.VRNX(NB,NZ,NY,NX))THEN

      call PhenoActiveBranch(NB,NZ,NY,NX)
      ELSE
      FDBK(NB,NZ,NY,NX)=0.0
      FDBKX(NB,NZ,NY,NX)=1.0
      DO 2805 K=1,25
      VCGR4(K,NB,NZ,NY,NX)=0.0
      VCGRO(K,NB,NZ,NY,NX)=0.0
2805  CONTINUE
      ENDIF
2900  CONTINUE
!
!     MINIMUM CANOPY STOMATAL RESISTANCE FROM CO2 CONCENTRATION
!     DIFFERENCE DIVIDED BY TOTAL CO2 FIXATION
!
!     RSX,RSMN=minimum canopy stomatal resistance to CO2,H2O (h m-1)
!     CH2O=total PEP(C4) or rubisco(C3) carboxylation rate
!     FRADP=fraction of radiation received by each PFT canopy
!     DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!     AREA=area of grid cell
!     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
!
      IF(CH2O.GT.ZEROP(NZ,NY,NX))THEN
      RSX=FRADP(NZ,NY,NX)*DCO2(NZ,NY,NX) &
      *AREA(3,NU(NY,NX),NY,NX)/(CH2O*3600.0)
      ELSE
      RSX=RSMH(NZ,NY,NX)*1.56
      ENDIF
      RSMN(NZ,NY,NX)=AMIN1(RSMH(NZ,NY,NX),AMAX1(RSMY,RSX*0.641))
      end subroutine PhotoActivePFT

      end module StomateMod
