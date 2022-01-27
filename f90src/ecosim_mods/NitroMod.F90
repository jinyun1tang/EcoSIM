module nitroMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : safe_adb
  implicit none

  private
  include "parameters.h"
  include "blkc.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk10.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk13d.h"
  include "blk15a.h"
  include "blk15b.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"
  include "blk19a.h"
  include "blk19b.h"
  include "blk21b.h"

  character(len=*), parameter :: mod_filename = __FILE__

! SUBSTRATE DECOMPOSITION BY MICROBIAL POPULATIONS
!
! ORAD=microbial radius (m), BIOS=microbial density (n m-3)
! BIOA=microbial surface area (m2 m-3), DCKI=inhibition of
! decomposition by microbial concentration (g C m-3)
! RCCX=maximum remobilization of microbial N (-)
! RCCY=maximum remobilization of microbial P (-)
!     RCCZ, RCCY = minimum, maximum remobilization of microbial C (-)
!     FPRIM, FPRIMM=fraction of nonstructural, microbial C,N,P
!     transferred with priming (-), OMGR=rate constant for
!     transferring nonstructural to structural microbial C (h-1)
!     OQKI=DOC product inhibition constant for decomposition (g C m-3)
!     H2KI=H2 product inhibition for methanogenesis (g H m-3)
!     COMKI, COMKM= Km to slow microbial decomposition, maintenance
!     respiration with low microbial C (g micr C g-1 subs C)
!     CKC=controls C remobilization of microbial C (g C g-1 C)
!     FOSCZ0, FOSCZL=rate constants for mixing surface (0) and
!     subsurface (L) litter (h-1),FMN=minimum ratio of total
!     biological demand for any substrate by any microbial population
!     DCKM0, DCKML=Km for SOC decomposition (g C g-1 soil)
!
      real(r8), PARAMETER :: ORAD=1.0E-06,BIOS=1.0E-06/(4.19*ORAD**3) &
      ,BIOA=BIOS*12.57*ORAD**2,DCKI=2.5,RCCX=0.833 &
      ,RCCQ=0.833,RCCZ=0.167,RCCY=0.833,FPRIM=5.0E-02,FPRIMM=1.0E-06 &
      ,OMGR=2.5E-01,OQKI=1.2E+03,H2KI=1.0,OAKI=12.0,COMKI=1.0E-03 &
      ,COMKM=1.0E-04,CKC=1.0E-03,FOSCZ0=2.0E-02,FOSCZL=2.0E-06 &
      ,FMN=1.0E-03,DCKM0=1.0E+03,DCKML=1.0E+03
!
!     SPECIFIC RESPIRATION RATES, M-M UPTAKE CONSTANTS,
!     STOICHIOMETRIC CONSTANTS FOR MICROBIAL REDOX REACTIONS
!
!     VMX*=specific oxidation rates (g C g-1C h-1_
!        O=all bacteria, F=fungi, M=acetotrophic methanogens
!        H=ammonia oxidizers, N=nitrite oxidizers, 4=methanotrophs
!        C=hydrogenotrophic methanogens
!     OQK*=Km for DOC uptake by heterotrophs (g C m-3)
!        M=all bacteria and fungi, A=acetate by fermenters
!        AM=acetate by acetotrophic methanogens
!     CCKM=Km for CO2 uptake, CCK4=Km for CH4 uptake (g C m-3)
!     Z*KM=Km for N uptake (g N m-3)
!        H=NH4 by nitrifiers, N=NO2 by nitrifiers
!        3=NO3 by denitrifiers, 2=NO2 by denitrifiers
!        1=N2O uptake by denitrifiers
!     Z4*=NH4 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
!       MX=maximum uptake rate, KU=Km, MN= minimum concentration
!     ZO*=NO3 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
!       MX=maximum uptake rate, KU=Km, MN= minimum concentration
!     HP*=H2PO4 uptake kinetics by all MFTs(g P m-2 h-1, g P m-3)
!       MX=maximum uptake rate, KU=Km, MN= minimum concentration
!     ZFKM=Km for N2 uptake by diazotrophs (g N m-3)
!     H2KM=Km for H2 uptake by hydrogenotrophic methanogens (g H m-3)
!     ECNH=efficiency CO2 conversion to biomass by ammonia oxidizers
!     ECNO=efficiency CO2 conversion to biomass by nitrite oxidizers
!     ECHO=efficiency CO2 conversion to biomass by methane oxidizers
!     ECN3,ECN2,ECN1=N2:O2 ratios for e- transfers to NO3, NO2 and N2O
!     by denitrifiers, RNFNI=parameter for nitrification inhibition
!     ZHKI=inhibition of nitrification inhibition by NH3 (g N m-3)
!     VMKI=product inhibn for NOx reduction by denitrifiers(g N m-3)
!     VHKI=product inhibn for NH3 oxidation by nitrifiers (g N m-3)
!     OXKA=Km for O2 uptake by nitrifiers(g O m-3)
!
      real(r8), PARAMETER :: VMXO=0.125,VMXF=0.125,VMXM=0.125 &
      ,VMXH=0.375,VMXN=0.25,VMX4=0.375,VMXC=0.125,OQKM=1.2E+01 &
      ,OQKA=1.2E+01,OQKAM=1.2E+01,CCKM=0.15,CCK4=1.2E-03,ZHKM=1.4 &
      ,ZNKM=1.4,Z3KM=1.4,Z2KM=1.4,Z1KM=0.014,Z4MX=5.0E-03 &
      ,Z4KU=0.40,Z4MN=0.0125,ZOMX=5.0E-03,ZOKU=0.35,ZOMN=0.03 &
      ,HPMX=1.0E-03,HPKU=0.075,HPMN=0.002,ZFKM=0.14,H2KM=0.01 &
      ,ECNH=0.30,ECNO=0.10,ECHO=0.75,ECN3=0.857,ECN2=0.857,ECN1=0.429 &
      ,RNFNI=2.0E-04,ZHKI=7.0E+03,VMKI=0.25,VHKI=15.0,OXKA=0.16
!
!     ENERGY REQUIREMENTS FOR MICROBIAL GROWTH AND
!     ENERGY YIELDS FROM REDUCTION OF O2, OC, CH4, NO3, N2
!
!     EOM*=energy requirements for microbial growth (kJ g-1 C)
!        C=aerobic bacteria, D=denitrifiers, G=fungi, F=fermenters
!        H=methanogens, N=diazotrophs
!     G*=free energy yields of redox reactions (kJ g-1 C or N)
!        O2X=DOC-CO2, H4X=CO2-CH4, CHX=DOC-acetate, O2A=acetate-CO2
!        C4X=acetate-CH4, COX=CO2-CH4, NOX=NO3-NO2,NO2-N2O,N2O-N2
!        N2X=N2-NH3
!     E*=growth respiration efficiency (-)(growth yield=1.0-E*)
!        N2X=aerobic N2 fixation, N2Y=anaerobic N2 fixation
!        O2X=aerobic bacteria (DOC), H4X=fermenters, O2G=fungi
!        O2D=denitrifiers (aerobic), NFX=diazotrophs
!        NOX= denitrifiers (anaerobic),O2A=aerobic bacteria (acetate)
!
      real(r8), PARAMETER :: EOMC=25.0,EOMD=37.5,EOMG=37.5,EOMF=75.0 &
      ,EOMH=25.0,EOMN=75.0,GO2X=37.5,GH4X=66.5,GCHX=4.50 &
      ,GO2A=GO2X-GCHX,GC4X=3.00,GCOX=11.00,GNOX=10.0 &
      ,GN2X=187.5,EN2X=GO2X/GN2X,EN2Y=GCHX/GN2X &
      ,EO2X=1.0/(1.0+GO2X/EOMC),EH4X=1.0/(1.0+GH4X/EOMC) &
      ,EO2G=1.0/(1.0+GO2X/EOMG),EO2D=1.0/(1.0+GO2X/EOMD) &
      ,ENFX=1.0/(1.0+GO2X/EOMN),ENOX=1.0/(1.0+GNOX/EOMC) &
      ,EO2A=1.0/(1.0+GO2A/EOMC)
!
!     SORPTION COEFFICIENTS
!
      real(r8), PARAMETER :: TSORP=0.5,HSORP=1.0
!
!     SPECIFIC DECOMPOSITION RATES
!
!     DOSA=rate constant for litter colonization by heterotrophs (g C g-1 C)
!     SP*= specific decomposition rate constant (g subs. C g-1 micr. C)
!       OHC=adsorbed SOC, OHA=adsorbed acetate, OSC=SOC
!       (K=0,M=1,4 woody litter, K=1,M=1,4 non-woody litter,
!       K=2,M=1,4 manure, K=3,M=1,1 POC, K=4,M=1,2 humus)
!       ORC (M=1,2) microbial residue, OMC (M=1,2) microbial biomass
!     RMOM=specific maintenance respiration (g C g-1 N h-1)
!
  real(r8), PARAMETER :: SPOHC=0.25,SPOHA=0.25,RMOM=0.010

  real(r8) :: ACTV,ACTVM,AECX,B,CHY1,CHNO2,CHNOB,COMC,CH4G1
  real(r8) :: CH4S1,CCH4S1,COXYS1,C,CNH4X,CNH4Y,CNO3X,CNO3Y
  real(r8) :: CH2PX,CH2PY,CH1PX,CH1PY,CGOXC
  real(r8) :: C3C,CNC,CPC,CGOMZ,CNOMX,CPOMX,COQCK,COSC,CNR
  real(r8) :: CPR,CGROMC,DIFOX,DCKD,DFNS,DOSAK,DCORPC,DC,DN
  real(r8) :: DP,EO2Q,ECHZ,FNH4S,FNHBS,FNO3S,FNO3B,FNO2S,FNO2B
      real(r8) :: FH1PS,FH1PB,FH2PS,FH2PB,FOXYX,FNH4X,FNB4X,FNO3X
      real(r8) :: FNB3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA,FSBSTC
      real(r8) :: FSBSTA,FSBST,FGOCP,FGOAP,FNH4,FNB4,FCN4S,FCN4B
      real(r8) :: FNO2,FNB2,FCN2S,FCN2B,FNO3,FNB3,FVMXDX,FN2O,FPH
      real(r8) :: FOSCXS,FOSCXD,FORGCX,GH2X,GH2F,GOAX,GOAF,GHAX
      real(r8) :: GOMX,GOMM,GH2H,H2GSX,H2POM,H2PBM,H1POM,H1PBM,H2P4M
      real(r8) :: H1P4M,HFLXD,OXYI,ORGCL,OHCQ,OXKX,OLSGL1,OXYG1
      real(r8) :: OXYS1,OQCZ3,OQCD3,OQCD3S,OQCD3B,OQCZ2,OQCD2,OQCD2S
      real(r8) :: OQCD2B,OQCZ1,OQCD1,OSRT,OQCI,OQCX,OQNX,OQPX,OQAX
      real(r8) :: OHCX,OHNX,OHPX,OHAX,ORGRL,ORGRLL,OSCXD,OMCXS,OMNXS
      real(r8) :: OMPXS,ORCXS,ORNXS,ORPXS,OQCXS,OQCHXS,OHCXS,OQAXS
      real(r8) :: OQAHXS,OHAXS,OQNXS,OQNHXS,OHNXS,OQPXS,OQPHXS,OHPXS
      real(r8) :: OSCXS,OSAXS,OSNXS,OSPXS,OC,ON,OP,OCH,ONH,OPH,ONX
      real(r8) :: OPX,OCA,OAH,RTK,RGOCY,RGOCZ,RGOAZ,RGOCX,RGOAX
      real(r8) :: RGOCP,RGOAP,RGOMP,ROXYSX,ROQCSX,ROQASX,RGOFY,RGOFZ
      real(r8) :: RGOFX,RGOGY,RGOGZ,RGOGX,RNNH4,RNNHB,RVOXP,RVOXPA
      real(r8) :: RVOXPB,RNNO2,RNNOB,RCH4L1,RCH4F1,RCH4S1,RVOXP1
      real(r8) :: RGOMP1,RCHDF,RUPMX,ROXYFX,ROXYLX,RRADO,RMPOX,ROXDFQ
      real(r8) :: RH2GZ,ROXYD,RDNO3X,RDNOBX,RDNOX,RDNOT,RGOM3X,RGOMD3
      real(r8) :: RDNO2X,RDN2X,RDN2T,RGOM2X,RGOMD2,RDN2OX,RGOM1X
      real(r8) :: RGOMD1,RINHP,RINHX,RINOP,RINOX,RIPOP,RIPOX,RIP1P
      real(r8) :: RIP1X,RINHPR,RINOPR,RIPOPR,RIP1PR,RMOMX,RMOMT
      real(r8) :: RGOMT,RXOMT,RGN2P,RCNO2,RCNOB,RCN2O
      real(r8) :: RCN2B,RCNO3,RCN3B,RCOQN,RHOSCM,STK,SPOMX,THETR
      real(r8) :: THETZ,TKSO,TOSC,TOSA,TORC,TOHC,TSRH,TOMA,TOMN
      real(r8) :: TFOXYX,TFNH4X,TFNO3X,TFNO2X,TFN2OX,TFP14X,TFPO4X
      real(r8) :: TFNH4B,TFNO3B,TFNO2B,TFP14B,TFPO4B,TCH4H,TCH4A
      real(r8) :: TFOQC,TFOQA,TRH2G,THETW1,TRINH,TRINO,TRIPO,TRIP1
      real(r8) :: TRINB,TRIOB,TRIPB,TRIB1,TRGOM,TRGOC,TRGOD,TRGOA
      real(r8) :: TRGOH,TUPOX,TRDN3,TRDNB,TRDN2,TRD2B,TRDNO,TRN2F
      real(r8) :: VOLWZ,VMXX,VMXA,VMX4S,VMX4B,VMX2S,VMX2B,VMXA1
      real(r8) :: VOLWCH,VOLWPM,VOLWOX,VOLPOX,VMXD3,VMXDXS,VMXDXB
      real(r8) :: VMXDXT,VMXD3S,VMXD3B,VMXD2,VMXD2S,VMXD2B,VMXD1
      real(r8) :: VMXD1S,VMXD4,VMXD4S,VMXD4B,VMXC4S,VMXC4B,VOLXX
      real(r8) :: VOLXW,VOLCX,VOLCW,VOLAX,VOLAW,XCO2,X,XFRK,XFRC
      real(r8) :: XFRN,XFRP,XFRA,XFMC,XFMN,XFMP,ZNFN4S,ZNFN4B,ZNO3SX
      real(r8) :: ZNO3BX,ZNO2SX,ZNO2BX,Z2OSX,ZNH4M,ZNHBM,ZNO3M,ZNOBM

      real(r8) :: CNOMA(7,0:5),CPOMA(7,0:5),OMA(7,0:5),FOMA(7,0:5) &
      ,FOMN(7,0:5),RDOSC(4,0:4),RDOSN(4,0:4),RDOSP(4,0:4),RHOSC(4,0:4) &
      ,RHOSN(4,0:4),RHOSP(4,0:4),RCOSC(4,0:4),RCOSN(4,0:4),RCOSP(4,0:4) &
      ,SPOSC(4,0:4),RDORC(2,0:4),RDORN(2,0:4),RDORP(2,0:4),SPORC(2) &
      ,RDOHC(0:4),RDOHN(0:4),RDOHP(0:4),RDOHA(0:4),CSORP(0:4),ZSORP(0:4) &
      ,PSORP(0:4),CSORPA(0:4),OSRH(0:4),RUPOX(7,0:5),RGN2F(7,0:5) &
      ,RGOMO(7,0:5),ROXYM(7,0:5),ROXYP(7,0:5),ROXYO(7,0:5) &
      ,RDNO3(7,0:5),RDNOB(7,0:5),RDNO2(7,0:5),RDN2B(7,0:5),RDN2O(7,0:5) &
      ,RGOMD(7,0:5),RMOMC(2,7,0:5),RINH4(7,0:5),RINO3(7,0:5) &
      ,RIPO4(7,0:5),RINB4(7,0:5),RINB3(7,0:5),RIPOB(7,0:5),FOMK(7,0:5) &
      ,RDOMC(2,7,0:5),RDOMN(2,7,0:5),RDOMP(2,7,0:5),RHOMC(2,7,0:5) &
      ,RHOMN(2,7,0:5),RHOMP(2,7,0:5),RCOMC(2,7,0:5),RCOMN(2,7,0:5) &
      ,RCOMP(2,7,0:5),CGOMC(7,0:5),CGOMN(7,0:5),RH2GX(7,0:5) &
      ,CGOMP(7,0:5),RDMMC(2,7,0:5),RHMMC(2,7,0:5),RCMMC(2,7,0:5) &
      ,RDMMN(2,7,0:5),RHMMN(2,7,0:5),RCMMN(2,7,0:5),RDMMP(2,7,0:5) &
      ,RHMMP(2,7,0:5),RCMMP(2,7,0:5),RCCMC(2,7,0:4) &
      ,RCCMN(2,7,0:4),RCCMP(2,7,0:4),RN2FX(7,0:5),TOMK(0:5) &
      ,TONK(0:5),TOPK(0:5),SPOMC(2),OMC2(7,0:5),TFNG(7,0:5),TFNR(7,0:5) &
      ,OMN2(7,0:5),FOM2(7,0:5),FOCA(0:4),FOAA(0:4),RXOMC(2,7,0:5) &
      ,RXOMN(2,7,0:5),RXOMP(2,7,0:5),R3OMC(2,7,0:5),R3OMN(2,7,0:5) &
      ,R3OMP(2,7,0:5),RXMMC(2,7,0:5),RXMMN(2,7,0:5),RXMMP(2,7,0:5) &
      ,R3MMC(2,7,0:5),R3MMN(2,7,0:5),R3MMP(2,7,0:5),WFN(7,0:5)

      real(r8) :: CGOQC(7,0:5),CGOAC(7,0:5),ROQCK(0:4),XOQCK(0:4) &
      ,EN2F(7),ORCT(0:4),OSCT(0:4),OSAT(0:4),ZNH4T(0:JZ),ZNO3T(0:JZ) &
      ,ZNO2T(0:JZ),H2P4T(0:JZ),RINH4R(7,0:5),RINO3R(7,0:5) &
      ,RIPO4R(7,0:5),FNH4XR(7,0:5),FNO3XR(7,0:5),FPO4XR(7,0:5) &
      ,RGOMY(7,0:5),CNQ(0:4),CPQ(0:4),CNH(0:4),CPH(0:4) &
      ,CNS(4,0:4),CPS(4,0:4),ROQCD(7,0:4),FORC(0:5),SPOMK(2),RMOMK(2) &
      ,CGOMS(2,7,0:5),CGONS(2,7,0:5),CGOPS(2,7,0:5),H1P4T(0:JZ) &
      ,TONX(0:5),TOPX(0:5),FCNK(0:4),FCPK(0:4),FP14XR(7,0:5) &
      ,RCO2X(7,0:5),RCH3X(7,0:5),RCH4X(7,0:5),RVOXA(7),RVOXB(7) &
      ,XOQCZ(0:4),XOQNZ(0:4),XOQPZ(0:4),XOQAZ(0:4) &
      ,XOMCZ(3,7,0:4),XOMNZ(3,7,0:4),XOMPZ(3,7,0:4) &
      ,FCN(7,0:5),FCP(7,0:5),FCNP(7,0:5),RIP14(7,0:5),RIP1B(7,0:5) &
      ,TCGOQC(0:5),TCGOAC(0:5),TCGOMN(0:5),TCGOMP(0:5) &
      ,TRN2ON(JY,JX),TRN2OD(JY,JX),TRN2GD(JY,JX),RIP14R(7,0:5)

  real(r8) :: ONL(4,0:4),OPL(4,0:4),EFIRE(2,21:22),DOSA(0:4)
!

      DATA DOSA/0.25E-03,0.25,0.25,0.25,0.25/
      DATA SPOSC/7.5,7.5,1.5,0.5,7.5,7.5,1.5,0.5 &
      ,7.5,7.5,1.5,0.5,0.05,0.00,0.00,0.00 &
      ,0.05,0.0167,0.00,0.00/
      DATA SPORC/7.5,1.5/
      DATA SPOMC/1.0E-02,0.1E-02/
      DATA EN2F/0.0,0.0,0.0,0.0,0.0,EN2X,EN2Y/
      DATA EFIRE/1.0,1.0,0.917,0.167/
      REAL(r8) :: WFNG,TFNX,TFNY,CNSHZ,CPSHZ,FRM

      public :: nitro
  contains

    SUBROUTINE nitro(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOIL BIOLOGICAL TRANSFORMATIONS
!
    implicit none
    integer, intent(in) :: I, J
    integer, intent(in) :: NHW,NHE,NVN,NVS

    integer :: L,NX,NY
!   begin_execution

    DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
!       IF(I.EQ.1.AND.J.EQ.1)THEN
!         TRN2ON(NY,NX)=0.0
!         TRN2OD(NY,NX)=0.0
!         TRN2GD(NY,NX)=0.0
!       ENDIF
!
!       VOLWZ=water volume used to calculate aqueous microbial
!       concentrations that drive microbial density effects on
!       decomposition
!
        DO 998 L=0,NL(NY,NX)
          call SoilBGCOneLayer(I,J,L,NY,NX)
998     CONTINUE
!       WRITE(20,3434)'RN2O',IYRC,I,J,(RN2O(L,NY,NX),L=0,NL(NY,NX))
3434    FORMAT(A8,3I4,20E12.4)
!
!       SOC LOSS IF FIRE OR REMOVAL EVENT IS ENTERED IN DISTURBANCE FILE
!
        call SOMRemovalByDisturbance(I,J,NY,NX)
9990  CONTINUE
9995  CONTINUE
    RETURN
    END subroutine nitro
!------------------------------------------------------------------------------------------

    subroutine SoilBGCOneLayer(I,J,L,NY,NX)
    implicit none
    integer, intent(in) :: I,J,L,NY,NX

      integer :: LL,K,KL
      integer :: M,N
!     begin_execution

      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
!     get KL, the number of mic-om complexes
      IF(L.EQ.0)THEN
!     surface litter layer
      KL=2
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETR=VOLW(0,NY,NX)/VOLR(NY,NX)
      THETZ=AMAX1(0.0,THETR-THETY(L,NY,NX))
      VOLWZ=THETZ*VOLR(NY,NX)
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,8825)'THETZ',I,J,L,THETR,THETZ,VOLWZ,VOLWRX(NY,NX)
!    2,VOLW(0,NY,NX),POROS(L,NY,NX),FC(0,NY,NX),WP(0,NY,NX)
!    3,THETY(L,NY,NX),PSISM(0,NY,NX),ORGC(0,NY,NX),VOLR(NY,NX)
8825  FORMAT(A8,3I4,20E12.4)
!     ENDIF
      ELSE
      VOLWZ=0.0
      ENDIF
      ELSE
!     non-surface layer
      KL=4
      THETZ=AMAX1(0.0,(AMIN1(AMAX1(0.5*POROS(L,NY,NX),FC(L,NY,NX)) &
      ,THETW(L,NY,NX))-THETY(L,NY,NX)))
      VOLWZ=THETZ*VOLY(L,NY,NX)
!     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.6)THEN
!     WRITE(*,8824)'THETZ',I,J,NX,NY,L,THETZ,THETW(L,NY,NX),VOLWZ
!    2,POROS(L,NY,NX),FC(L,NY,NX),WP(L,NY,NX),THETY(L,NY,NX)
!    3,VOLW(L,NY,NX),VOLWH(L,NY,NX),VOLY(L,NY,NX),VOLT(L,NY,NX)
!    4,DTBLX(NY,NX)
8824  FORMAT(A8,5I4,20E12.4)
!     ENDIF
      ENDIF
!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH OFFSET FOR THERMAL ADAPTATION
!

      call StageBGCEnvironCondition(KL,L,NY,NX)
!
!     RESPIRATION BY MICROBIAL POPULATIONS
!
      TFOXYX=0.0
      TFNH4X=0.0
      TFNO3X=0.0
      TFNO2X=0.0
      TFN2OX=0.0
      TFP14X=0.0
      TFPO4X=0.0
      TFNH4B=0.0
      TFNO3B=0.0
      TFNO2B=0.0
      TFP14B=0.0
      TFPO4B=0.0
      TCH4H=0.0
      TCH4A=0.0
      TFOQC=0.0
      TFOQA=0.0
      TRH2G=0.0
      IF(L.NE.0)THEN
      LL=L
      ELSE
      LL=NU(NY,NX)
      ENDIF

      call MicrobialCatabolism(L,LL,NY,NX)
!
!     ChemoDenitrification
      call ChemoDenitrification(L,NY,NX)
!
!     DECOMPOSITION
!
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
      DO 1870 K=0,KL
      ROQCK(K)=0.0
      DO 1875 N=1,7
      ROQCK(K)=ROQCK(K)+ROQCD(N,K)
1875  CONTINUE
      XOQCK(K)=0.0
      XOQCZ(K)=0.0
      XOQNZ(K)=0.0
      XOQPZ(K)=0.0
      XOQAZ(K)=0.0
      DO 845 N=1,7
      DO 845 M=1,3
      XOMCZ(M,N,K)=0.0
      XOMNZ(M,N,K)=0.0
      XOMPZ(M,N,K)=0.0
845   CONTINUE
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.LE.1)THEN
!     WRITE(*,4443)'PRIM1',I,J,NX,NY,L,K,ROQCK(K)
!    2,XOQCK(K),OQC(K,L,NY,NX),XOQCZ(K),OQN(K,L,NY,NX),XOQNZ(K)
!    3,OQP(K,L,NY,NX),XOQPZ(K),OQA(K,L,NY,NX),XOQAZ(K)
!     ENDIF
1870  CONTINUE
!
!     PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C
      call OMTransferForPriming(KL,L,NY,NX)
!
!     TRANSFER ALL PRIMING AMONG ALL K
!
!     TOQCK=total respiration of DOC+DOA in soil layer
!     ROQCK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
      DO 1790 K=0,KL
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
!     WRITE(*,4443)'PRIM2',I,J,NX,NY,L,K,ROQCK(K)
!    2,XOQCK(K),OQC(K,L,NY,NX),XOQCZ(K),OQN(K,L,NY,NX),XOQNZ(K)
!    3,OQP(K,L,NY,NX),XOQPZ(K),OQA(K,L,NY,NX),XOQAZ(K),TOMK(K)
!    3,TONK(K),TOPK(K),TONX(K),TOPX(K),CNOMX,CPOMX,FCNK(K),FCPK(K)
!    4,TOQCK(L,NY,NX)
4443  FORMAT(A8,6I4,20E12.4)
!     ENDIF
!
!     DECOMPOSITION OF ORGANIC SUBSTRATES
!
      call SolidOMDecomposition(K,L,NY,NX)
!
!     DOC ADSORPTION - DESORPTION
!
      call DOMSorption(K,L,NY,NX)
!     IF(I.EQ.116)THEN
!     WRITE(*,591)'CSORP',I,J,NX,NY,L,K,CSORP(K),CSORPA(K)
!    1,OQC(K,L,NY,NX),OHC(K,L,NY,NX),OQA(K,L,NY,NX),OHA(K,L,NY,NX)
!    2,OQC(K,L,NY,NX)/VOLWM(NPH,L,NY,NX)
!    2,OQA(K,L,NY,NX)/VOLWM(NPH,L,NY,NX)
!    3,OHC(K,L,NY,NX)/BKVL(L,NY,NX),OHA(K,L,NY,NX)/BKVL(L,NY,NX)
!    4,BKVL(L,NY,NX),VOLWM(NPH,L,NY,NX),FOCA(K),FOAA(K)
!    5,FOSRH(K,L,NY,NX),TCGOQC(K),OQCX
591   FORMAT(A8,6I4,40E12.4)
!     ENDIF
1790  CONTINUE

      call RedistDecompositionProduct(KL,L,NY,NX)
!
!     MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION

      call MicrobialAnabolicUpdate(L,NY,NX)
!
!     MICROBIAL COLONIZATION OF NEW LITTER
!
      call MicrobialLitterColonization(KL,L,NY,NX)
!
!     AGGREGATE ALL TRANSFORMATIONS CALCULATED ABOVE FOR EACH N,K
!
      call AggregateTransformations(L,NY,NX)
!     IF(ISALTG.NE.0)THEN
!     XZHYS(L,NY,NX)=XZHYS(L,NY,NX)+0.1429*(RVOXA(1)+RVOXB(1)
!    2-TRDN3-TRDNB)-0.0714*(TRDN2+TRD2B+TRDNO)
!     ENDIF
!     IF(L.EQ.0)THEN
!     WRITE(*,2323)'XNH4S',I,J,L,XNH4S(L,NY,NX)
!    2,TRINH,RVOXA(1),VLNH4(L,NY,NX),TRDN2
!     WRITE(*,2323)'XNO3S',I,J,L,XNO3S(L,NY,NX)
!    2,TRINO,RVOXA(2),VLNO3(L,NY,NX),TRDN3,RCNO3
!     WRITE(*,2323)'XH2PS',I,J,L,XH2PS(L,NY,NX)
!    2,RIPOT,TRIPO,VLPO4(L,NY,NX)
!     WRITE(*,2323)'XNO2B',I,J,L,XNO2B(L,NY,NX),RVOXB(1)
!    2,VLNHB(L,NY,NX),RVOXB(2),VLNOB(L,NY,NX),TRDNB,TRD2B,RCNOB
!     ENDIF
!     WRITE(*,2324)'XOQCS',I,J,NX,NY,L,(XOQCS(K,L,NY,NX),K=0,4)
2324  FORMAT(A8,5I4,12E12.4)
      ELSE
      RCO2O(L,NY,NX)=0.0
      RCH4O(L,NY,NX)=0.0
      RH2GO(L,NY,NX)=0.0
      RUPOXO(L,NY,NX)=0.0
      RN2G(L,NY,NX)=0.0
      RN2O(L,NY,NX)=0.0
      XNH4S(L,NY,NX)=0.0
      XNO3S(L,NY,NX)=0.0
      XNO2S(L,NY,NX)=0.0
      XH2PS(L,NY,NX)=0.0
      XH1PS(L,NY,NX)=0.0
      XNH4B(L,NY,NX)=0.0
      XNO3B(L,NY,NX)=0.0
      XNO2B(L,NY,NX)=0.0
      XH2BS(L,NY,NX)=0.0
      XH1BS(L,NY,NX)=0.0
      XN2GS(L,NY,NX)=0.0
      ENDIF
!
!     MIX LITTER C BETWEEN ADJACENT SOIL LAYERS L AND LL
!
      call VerticalLitterMixLvsLL(I,J,L,NY,NX)
      ELSE
      RCO2O(L,NY,NX)=0.0
      RCH4O(L,NY,NX)=0.0
      RH2GO(L,NY,NX)=0.0
      RUPOXO(L,NY,NX)=0.0
      RN2G(L,NY,NX)=0.0
      RN2O(L,NY,NX)=0.0
      XNH4S(L,NY,NX)=0.0
      XNO3S(L,NY,NX)=0.0
      XNO2S(L,NY,NX)=0.0
      XH2PS(L,NY,NX)=0.0
      XH1PS(L,NY,NX)=0.0
      XNH4B(L,NY,NX)=0.0
      XNO3B(L,NY,NX)=0.0
      XNO2B(L,NY,NX)=0.0
      XH2BS(L,NY,NX)=0.0
      XH1BS(L,NY,NX)=0.0
      XN2GS(L,NY,NX)=0.0
      ENDIF
      end subroutine SoilBGCOneLayer
!------------------------------------------------------------------------------------------

      subroutine SOMRemovalByDisturbance(I,J,NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: I,J,NY,NX

      integer :: L,K,M,N,IFLGJ,NLL
!     begin_execution

      IF(J.EQ.INT(ZNOON(NY,NX)).AND.(ITILL(I,NY,NX).EQ.21 &
      .OR.ITILL(I,NY,NX).EQ.22))THEN
      IF(ITILL(I,NY,NX).EQ.22)THEN
      IFLGS(NY,NX)=1
      IFLGJ=0
      NLL=-1
      DO 2945 L=0,NL(NY,NX)
!     WRITE(*,9494)'FIRE',I,J,L,IFLGJ,NLL,ORGC(L,NY,NX),THETW(L,NY,NX)
!    2,FVLWB,CORGC(L,NY,NX),FORGC,DPTH(L,NY,NX),BKDS(L,NY,NX)
!    3,VOLY(L,NY,NX),DTBLX(NY,NX),DCORP(I,NY,NX)
9494  FORMAT(A8,5I6,12E12.4)
      IF(L.EQ.0.OR.L.GE.NUM(NY,NX))THEN
      IF(IFLGJ.EQ.1)THEN
      GO TO 2946
      ELSEIF(THETW(L,NY,NX).GT.FVLWB.OR.CORGC(L,NY,NX).LE.FORGC)THEN
      IFLGJ=1
      ELSE
      NLL=L
      ENDIF
      ENDIF
2945  CONTINUE
      ELSE
      NLL=0
      ENDIF
2946  CONTINUE
      DO 2950 L=0,NLL
      IF(NLL.GE.0)THEN
      IF(ITILL(I,NY,NX).EQ.22)THEN
      IF(L.EQ.0)THEN
      FORGCX=0.0
      ELSE
      FORGCX=FORGC
      ENDIF
      DCORPC=AMIN1(0.999,DCORP(I,NY,NX))*(CORGC(L,NY,NX)-FORGCX) &
      /(AMAX1(CORGC(L,NY,NX),0.55E+06)-FORGCX)
      ELSE
      DCORPC=AMIN1(0.999,DCORP(I,NY,NX))
      ENDIF
!     VOLWOU=VOLWOU+DCORPC*VOLW(L,NY,NX)
!     HEATOU=HEATOU+DCORPC*4.19*TKS(L,NY,NX)*VOLW(L,NY,NX)
!     VOLW(L,NY,NX)=VOLW(L,NY,NX)-DCORPC*VOLW(L,NY,NX)
!     WRITE(*,9696)'BURN',I,J,L,NLL,ITILL(I,NY,NX)
!    2,CORGC(L,NY,NX),ORGC(L,NY,NX)
!    2,FORGCX,DCORPC,DCORP(I,NY,NX),VOLW(L,NY,NX),BKDS(L,NY,NX)
9696  FORMAT(A8,5I6,12E12.4)
      OC=0.0
      ON=0.0
      OP=0.0
      DC=0.0
      DN=0.0
      DP=0.0
      DO 2955 K=0,4
      DO 2955 M=1,4
      ONL(M,K)=0.0
      OPL(M,K)=0.0
2955  CONTINUE
      DO 2970 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
      DO 2960 N=1,7
      DO 2960 M=1,3
      OCH=DCORPC*OMC(M,N,K,L,NY,NX)
      ONH=DCORPC*OMN(M,N,K,L,NY,NX)
      OPH=DCORPC*OMP(M,N,K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSEIF(K.LE.4)THEN
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ELSEIF(K.EQ.5)THEN
      ONL(4,1)=ONL(4,1)+ONH-ONX
      OPL(4,1)=OPL(4,1)+OPH-OPX
      ENDIF
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)-OCH
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)-ONH
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)-OPH
      DC=DC+OMC(M,N,K,L,NY,NX)
      DN=DN+OMN(M,N,K,L,NY,NX)
      DP=DP+OMP(M,N,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
2960  CONTINUE
      ENDIF
2970  CONTINUE
!
!     REMOVE MICROBIAL RESIDUE
!
      DO 2900 K=0,4
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 2940 M=1,2
      OCH=DCORPC*ORC(M,K,L,NY,NX)
      ONH=DCORPC*ORN(M,K,L,NY,NX)
      OPH=DCORPC*ORP(M,K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-OCH
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ONH
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-OPH
      DC=DC+ORC(M,K,L,NY,NX)
      DN=DN+ORN(M,K,L,NY,NX)
      DP=DP+ORP(M,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
2940  CONTINUE
!
!     REMOVE DOC, DON, DOP
!
      OCH=DCORPC*OQC(K,L,NY,NX)
      OCA=DCORPC*OQA(K,L,NY,NX)
      ONH=DCORPC*OQN(K,L,NY,NX)
      OPH=DCORPC*OQP(K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OCH
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OCA
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ONH
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OPH
      OC=OC+OCH+OCA
      ON=ON+ONX
      OP=OP+OPX
      OCH=DCORPC*OQCH(K,L,NY,NX)
      ONH=DCORPC*OQNH(K,L,NY,NX)
      OPH=DCORPC*OQPH(K,L,NY,NX)
      OAH=DCORPC*OQAH(K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OCH
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-ONH
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OPH
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OAH
      OC=OC+OCH+OAH
      ON=ON+ONX
      OP=OP+OPX
!
!     REMOVE ADSORBED OM
!
      OCH=DCORPC*OHC(K,L,NY,NX)
      ONH=DCORPC*OHN(K,L,NY,NX)
      OPH=DCORPC*OHP(K,L,NY,NX)
      OAH=DCORPC*OHA(K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OCH
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-ONH
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OPH
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OAH
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
      +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
!
!     REMOVE RESIDUE
!
      DO 2930 M=1,4
      OCH=DCORPC*OSC(M,K,L,NY,NX)
      OCA=DCORPC*OSA(M,K,L,NY,NX)
      ONH=DCORPC*OSN(M,K,L,NY,NX)
      OPH=DCORPC*OSP(M,K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      ONL(M,K)=ONL(M,K)+ONH-ONX
      OPL(M,K)=OPL(M,K)+OPH-OPX
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OCH
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OCA
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-ONH
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OPH
      DC=DC+OSC(M,K,L,NY,NX)
      DN=DN+OSN(M,K,L,NY,NX)
      DP=DP+OSP(M,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
2930  CONTINUE
      ENDIF
2900  CONTINUE
!
!     ADD UNBURNED N,P TO ORG N, ORG P
!
      DO 2905 K=0,4
      DO 2905 M=1,4
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+ONL(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+OPL(M,K)
      DN=DN+ONL(M,K)
      DP=DP+OPL(M,K)
2905  CONTINUE
!
!     REMOVE FERTILIZER IN RESIDUE
!
      IF(ITILL(I,NY,NX).EQ.21)THEN
      ON=ON+DCORPC*(ZNH4S(L,NY,NX)+ZNH3S(L,NY,NX) &
      +ZNO3S(L,NY,NX)+ZNO2S(L,NY,NX))
      OP=OP+DCORPC*(H1PO4(L,NY,NX)+H2PO4(L,NY,NX))
      ZNH4S(L,NY,NX)=(1.0-DCORPC)*ZNH4S(L,NY,NX)
      ZNH3S(L,NY,NX)=(1.0-DCORPC)*ZNH3S(L,NY,NX)
      ZNO3S(L,NY,NX)=(1.0-DCORPC)*ZNO3S(L,NY,NX)
      ZNO2S(L,NY,NX)=(1.0-DCORPC)*ZNO2S(L,NY,NX)
      H1PO4(L,NY,NX)=(1.0-DCORPC)*H1PO4(L,NY,NX)
      H2PO4(L,NY,NX)=(1.0-DCORPC)*H2PO4(L,NY,NX)
      XN4(L,NY,NX)=(1.0-DCORPC)*XN4(L,NY,NX)
      PALPO(L,NY,NX)=(1.0-DCORPC)*PALPO(L,NY,NX)
      PFEPO(L,NY,NX)=(1.0-DCORPC)*PFEPO(L,NY,NX)
      PCAPD(L,NY,NX)=(1.0-DCORPC)*PCAPD(L,NY,NX)
      PCAPH(L,NY,NX)=(1.0-DCORPC)*PCAPH(L,NY,NX)
      PCAPM(L,NY,NX)=(1.0-DCORPC)*PCAPM(L,NY,NX)
      ZNH4FA(L,NY,NX)=(1.0-DCORPC)*ZNH4FA(L,NY,NX)
      ZNH3FA(L,NY,NX)=(1.0-DCORPC)*ZNH3FA(L,NY,NX)
      ZNHUFA(L,NY,NX)=(1.0-DCORPC)*ZNHUFA(L,NY,NX)
      ZNO3FA(L,NY,NX)=(1.0-DCORPC)*ZNO3FA(L,NY,NX)
      ENDIF
      ORGC(L,NY,NX)=DC
      ORGN(L,NY,NX)=DN
      IF(L.EQ.0)THEN
      HFLXD=4.19E-06*(ORGCX(L,NY,NX)-ORGC(L,NY,NX))*TKS(L,NY,NX)
      HEATOU=HEATOU+HFLXD
      ENDIF
!     IF(L.EQ.0)THEN
!     VHCP(0,NY,NX)=2.496E-06*ORGC(0,NY,NX)+4.19*VOLW(0,NY,NX)
!    2+1.9274*VOLI(0,NY,NX)
!     ELSE
!     VHCP(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW(L,NY,NX)+VOLWH(L,NY,NX))
!    2+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
!     ENDIF
      IF(ITILL(I,NY,NX).EQ.21)THEN
      TCOU=TCOU+OC
      TZOU=TZOU+ON
      TPOU=TPOU+OP
      UDOCQ(NY,NX)=UDOCQ(NY,NX)+OC
      UDONQ(NY,NX)=UDONQ(NY,NX)+ON
      UDOPQ(NY,NX)=UDOPQ(NY,NX)+OP
      TNBP(NY,NX)=TNBP(NY,NX)-OC
      ELSEIF(ITILL(I,NY,NX).EQ.22)THEN
      CO2GIN=CO2GIN-OC
      OXYGIN=OXYGIN+2.667*OC
      OXYGOU=OXYGOU+2.667*OC
      TZOU=TZOU+ON
      TPOU=TPOU+OP
      UCO2F(NY,NX)=UCO2F(NY,NX)-(1.0-FCH4F)*OC
      UCH4F(NY,NX)=UCH4F(NY,NX)-FCH4F*OC
      UOXYF(NY,NX)=UOXYF(NY,NX)+(1.0-FCH4F)*2.667*OC
      UNH3F(NY,NX)=UNH3F(NY,NX)-ON
      UN2OF(NY,NX)=UN2OF(NY,NX)-0.0
      UPO4F(NY,NX)=UPO4F(NY,NX)-OP
      TNBP(NY,NX)=TNBP(NY,NX)-OC
      ENDIF
      ENDIF
2950  CONTINUE
      ENDIF
      end subroutine SOMRemovalByDisturbance
!------------------------------------------------------------------------------------------

      subroutine VerticalLitterMixLvsLL(I,J,L,NY,NX)

      implicit none
      integer, intent(in) :: I,J,L,NY,NX

      integer :: LL,LN
!     begin_execution

      IF(FOSCZ0.GT.ZERO)THEN
!     ORGR=total litter C
!     FOSCZ0=rate constant for mixing surface litter
!     FOSCXS=mixing fraction for surface litter
!     TOQCK=total active biomass respiration activity
!     TFNX=temperature function
!     VOLX=soil layer volume
!     OSCXD=mixing required for equilibrating litter concentration
!     FOSCXD=mixing fraction for equilibrating subsurface litter
!     FOSCXS=mixing fraction for subsurface litter
!
!     IF(I.EQ.116)THEN
!     WRITE(*,336)'LAYER',I,J,L,TOQCK(L,NY,NX),TOMA,TFNX,TOMA*TFNX
336   FORMAT(A8,3I4,20E12.4)
!     ENDIF
      IF(L.LT.NL(NY,NX))THEN
!     get mixing rate
      IF(L.EQ.0)THEN
      LL=NU(NY,NX)
      IF(ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOSCXS=AMIN1(1.0,FOSCZ0/ORGR(L,NY,NX)*TOQCK(L,NY,NX))
      ELSE
      FOSCXS=0.0
      ENDIF
      ELSE
      DO 1100 LN=L+1,NL(NY,NX)
      IF(VOLX(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
      LL=LN
      GO TO 1101
      ENDIF
1100  CONTINUE
1101  CONTINUE
      ORGRL=AMAX1(0.0,ORGR(L,NY,NX))
      ORGRLL=AMAX1(0.0,ORGR(LL,NY,NX))
      OSCXD=(ORGRL*VOLT(LL,NY,NX)-ORGRLL*VOLT(L,NY,NX)) &
      /(VOLT(L,NY,NX)+VOLT(LL,NY,NX))
      IF(OSCXD.GT.0.0.AND.ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOSCXD=OSCXD/ORGR(L,NY,NX)
      ELSEIF(OSCXD.LT.0.0.AND.ORGR(LL,NY,NX).GT.ZEROS(NY,NX))THEN
      FOSCXD=OSCXD/ORGR(LL,NY,NX)
      ELSE
      FOSCXD=0.0
      ENDIF
      IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      FOSCXS=FOSCZL*FOSCXD*TOQCK(L,NY,NX)/VOLT(L,NY,NX)
      ELSE
      FOSCXS=0.0
      ENDIF
      ENDIF
!     IF(L.EQ.3)THEN
!     WRITE(*,1115)'MIX',I,J,NX,NY,L,LL,FOSCXS,FOSCZ0
!    2,FOSCZL,FOSCXD,TOQCK(L,NY,NX),VOLT(L,NY,NX)
!    2,OSCXD,ORGR(L,NY,NX),ORGR(LL,NY,NX)
!    3,TKS(L,NY,NX)
1115  FORMAT(A8,6I4,30E12.4)
!     ENDIF
!     apply mixing
      call ApplyVerticalMix(FOSCXS,L,LL,NY,NX)
      ENDIF
!     IF((I/1)*1.EQ.I.AND.J.EQ.19.AND.L.LE.5)THEN
!     WRITE(*,2123)'TOTALL',I,J,NX,NY,L,TFOXYX,TFNH4X
!    2,TFNO3X,TFPO4X,TFNH4B,TFNO3B,TFPO4B,TFNO2X,TFNO2B
!    3,TFOQC,TFOQA
2123  FORMAT(A8,5I4,12E15.4)
!     ENDIF
      ENDIF
      end subroutine VerticalLitterMixLvsLL
!------------------------------------------------------------------------------------------

      subroutine ApplyVerticalMix(FOSCXS,L,LL,NY,NX)

      implicit none
      real(r8), intent(in) :: FOSCXS
      integer, intent(in) :: L,LL,NY,NX

      integer :: K,M,N
!     begin_execution
      IF(FOSCXS.GT.ZERO)THEN
      DO 7971 K=1,2
      DO 7961 N=1,7
      DO 7962 M=1,3
      IF(FOSCXS.GT.0.0)THEN
      OMCXS=FOSCXS*AMAX1(0.0,OMC(M,N,K,L,NY,NX))
      OMNXS=FOSCXS*AMAX1(0.0,OMN(M,N,K,L,NY,NX))
      OMPXS=FOSCXS*AMAX1(0.0,OMP(M,N,K,L,NY,NX))
      ELSE
      OMCXS=FOSCXS*AMAX1(0.0,OMC(M,N,K,LL,NY,NX))
      OMNXS=FOSCXS*AMAX1(0.0,OMN(M,N,K,LL,NY,NX))
      OMPXS=FOSCXS*AMAX1(0.0,OMP(M,N,K,LL,NY,NX))
      ENDIF
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)-OMCXS
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)-OMNXS
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)-OMPXS
      OMC(M,N,K,LL,NY,NX)=OMC(M,N,K,LL,NY,NX)+OMCXS
      OMN(M,N,K,LL,NY,NX)=OMN(M,N,K,LL,NY,NX)+OMNXS
      OMP(M,N,K,LL,NY,NX)=OMP(M,N,K,LL,NY,NX)+OMPXS
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,5558)'OMX',I,J,NX,NY,L,LL,K,N,M,OMC(M,N,K,L,NY,NX)
!    2,OMN(M,N,K,L,NY,NX),OMP(M,N,K,L,NY,NX),OMC(M,N,K,LL,NY,NX)
!    2,OMN(M,N,K,LL,NY,NX),OMP(M,N,K,LL,NY,NX)
!    3,OMCXS,OMNXS,OMPXS,FOSCXS
5558  FORMAT(A8,9I4,12E12.4)
!     ENDIF
7962  CONTINUE
7961  CONTINUE
7971  CONTINUE
      DO 7901 K=1,2
      DO 7941 M=1,2
      IF(FOSCXS.GT.0.0)THEN
      ORCXS=FOSCXS*AMAX1(0.0,ORC(M,K,L,NY,NX))
      ORNXS=FOSCXS*AMAX1(0.0,ORN(M,K,L,NY,NX))
      ORPXS=FOSCXS*AMAX1(0.0,ORP(M,K,L,NY,NX))
      ELSE
      ORCXS=FOSCXS*AMAX1(0.0,ORC(M,K,LL,NY,NX))
      ORNXS=FOSCXS*AMAX1(0.0,ORN(M,K,LL,NY,NX))
      ORPXS=FOSCXS*AMAX1(0.0,ORP(M,K,LL,NY,NX))
      ENDIF
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-ORCXS
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ORNXS
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-ORPXS
      ORC(M,K,LL,NY,NX)=ORC(M,K,LL,NY,NX)+ORCXS
      ORN(M,K,LL,NY,NX)=ORN(M,K,LL,NY,NX)+ORNXS
      ORP(M,K,LL,NY,NX)=ORP(M,K,LL,NY,NX)+ORPXS
!     IF(L.EQ.3.AND.K.EQ.2)THEN
!     WRITE(*,7942)'ORC',I,J,L,LL,K,M,ORC(M,K,L,NY,NX)
!    2,ORC(M,K,LL,NY,NX),ORCXS,FOSCXS
7942  FORMAT(A8,6I4,20E12.4)
!     ENDIF
7941  CONTINUE
      IF(FOSCXS.GT.0.0)THEN
      OQCXS=FOSCXS*AMAX1(0.0,OQC(K,L,NY,NX))
      OQCHXS=FOSCXS*AMAX1(0.0,OQCH(K,L,NY,NX))
      OHCXS=FOSCXS*AMAX1(0.0,OHC(K,L,NY,NX))
      OQAXS=FOSCXS*AMAX1(0.0,OQA(K,L,NY,NX))
      OQAHXS=FOSCXS*AMAX1(0.0,OQAH(K,L,NY,NX))
      OHAXS=FOSCXS*AMAX1(0.0,OHA(K,L,NY,NX))
      OQNXS=FOSCXS*AMAX1(0.0,OQN(K,L,NY,NX))
      OQNHXS=FOSCXS*AMAX1(0.0,OQNH(K,L,NY,NX))
      OHNXS=FOSCXS*AMAX1(0.0,OHN(K,L,NY,NX))
      OQPXS=FOSCXS*AMAX1(0.0,OQP(K,L,NY,NX))
      OQPHXS=FOSCXS*AMAX1(0.0,OQPH(K,L,NY,NX))
      OHPXS=FOSCXS*AMAX1(0.0,OHP(K,L,NY,NX))
      ELSE
      OQCXS=FOSCXS*AMAX1(0.0,OQC(K,LL,NY,NX))
      OQCHXS=FOSCXS*AMAX1(0.0,OQCH(K,LL,NY,NX))
      OHCXS=FOSCXS*AMAX1(0.0,OHC(K,LL,NY,NX))
      OQAXS=FOSCXS*AMAX1(0.0,OQA(K,LL,NY,NX))
      OQAHXS=FOSCXS*AMAX1(0.0,OQAH(K,LL,NY,NX))
      OHAXS=FOSCXS*AMAX1(0.0,OHA(K,LL,NY,NX))
      OQNXS=FOSCXS*AMAX1(0.0,OQN(K,LL,NY,NX))
      OQNHXS=FOSCXS*AMAX1(0.0,OQNH(K,LL,NY,NX))
      OHNXS=FOSCXS*AMAX1(0.0,OHN(K,LL,NY,NX))
      OQPXS=FOSCXS*AMAX1(0.0,OQP(K,LL,NY,NX))
      OQPHXS=FOSCXS*AMAX1(0.0,OQPH(K,LL,NY,NX))
      OHPXS=FOSCXS*AMAX1(0.0,OHP(K,LL,NY,NX))
      ENDIF
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OQCXS
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OQCHXS
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OHCXS
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OQAXS
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OQAHXS
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OHAXS
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-OQNXS
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-OQNHXS
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-OHNXS
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OQPXS
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OQPHXS
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OHPXS
      OQC(K,LL,NY,NX)=OQC(K,LL,NY,NX)+OQCXS
      OQCH(K,LL,NY,NX)=OQCH(K,LL,NY,NX)+OQCHXS
      OHC(K,LL,NY,NX)=OHC(K,LL,NY,NX)+OHCXS
      OQA(K,LL,NY,NX)=OQA(K,LL,NY,NX)+OQAXS
      OQAH(K,LL,NY,NX)=OQAH(K,LL,NY,NX)+OQAHXS
      OHA(K,LL,NY,NX)=OHA(K,LL,NY,NX)+OHAXS
      OQN(K,LL,NY,NX)=OQN(K,LL,NY,NX)+OQNXS
      OQNH(K,LL,NY,NX)=OQNH(K,LL,NY,NX)+OQNHXS
      OHN(K,LL,NY,NX)=OHN(K,LL,NY,NX)+OHNXS
      OQP(K,LL,NY,NX)=OQP(K,LL,NY,NX)+OQPXS
      OQPH(K,LL,NY,NX)=OQPH(K,LL,NY,NX)+OQPHXS
      OHP(K,LL,NY,NX)=OHP(K,LL,NY,NX)+OHPXS
      DO 7931 M=1,4
      IF(FOSCXS.GT.0.0)THEN
      OSCXS=FOSCXS*AMAX1(0.0,OSC(M,K,L,NY,NX))
      OSAXS=FOSCXS*AMAX1(0.0,OSA(M,K,L,NY,NX))
      OSNXS=FOSCXS*AMAX1(0.0,OSN(M,K,L,NY,NX))
      OSPXS=FOSCXS*AMAX1(0.0,OSP(M,K,L,NY,NX))
      ELSE
      OSCXS=FOSCXS*AMAX1(0.0,OSC(M,K,LL,NY,NX))
      OSAXS=FOSCXS*AMAX1(0.0,OSA(M,K,LL,NY,NX))
      OSNXS=FOSCXS*AMAX1(0.0,OSN(M,K,LL,NY,NX))
      OSPXS=FOSCXS*AMAX1(0.0,OSP(M,K,LL,NY,NX))
      ENDIF
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OSCXS
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OSAXS
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-OSNXS
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OSPXS
      OSC(M,K,LL,NY,NX)=OSC(M,K,LL,NY,NX)+OSCXS
      OSA(M,K,LL,NY,NX)=OSA(M,K,LL,NY,NX)+OSAXS
      OSN(M,K,LL,NY,NX)=OSN(M,K,LL,NY,NX)+OSNXS
      OSP(M,K,LL,NY,NX)=OSP(M,K,LL,NY,NX)+OSPXS
7931  CONTINUE
7901  CONTINUE
      ENDIF
      end subroutine ApplyVerticalMix
!------------------------------------------------------------------------------------------

      subroutine StageBGCEnvironCondition(KL,L,NY,NX)
      implicit none
      integer, intent(in) :: KL,L,NY,NX

      integer :: K
      integer :: M,N
!     begin_execution
!     TKS=soil temperature
!     OFFSET=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TFNX,TFNY=temperature function for growth,maintenance respiration
!
      TKSO=TKS(L,NY,NX)+OFFSET(NY,NX)
      RTK=8.3143*TKSO
      STK=710.0*TKSO
      ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFNX=EXP(25.229-62500/RTK)/ACTV
      ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
      TFNY=EXP(25.214-62500/RTK)/ACTVM
!
!     OXYI=inhibition of fermenters by O2
!     ORGCL=SOC used to calculate microbial concentration
!
      OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS(L,NY,NX)+2.5)))
      ORGCL=AMIN1(1.0E+05*BKVL(L,NY,NX),ORGC(L,NY,NX))
!
!     TOTAL MINERAL NH4, NO3 AND PO4
!
!     allocate NH4, NO3, HPO4, H2PO4 to non-band and band fractions
!
      ZNH4T(L)=AMAX1(0.0,ZNH4S(L,NY,NX))+AMAX1(0.0,ZNH4B(L,NY,NX))
!     IF(ZNH4T(L).GT.ZEROS(NY,NX))THEN
!     FNH4S=AMAX1(0.0,ZNH4S(L,NY,NX))/ZNH4T(L)
!     FNHBS=AMAX1(0.0,ZNH4B(L,NY,NX))/ZNH4T(L)
!     ELSE
      FNH4S=VLNH4(L,NY,NX)
      FNHBS=VLNHB(L,NY,NX)
!     ENDIF
      ZNO3T(L)=AMAX1(0.0,ZNO3S(L,NY,NX))+AMAX1(0.0,ZNO3B(L,NY,NX))
!     IF(ZNO3T(L).GT.ZEROS(NY,NX))THEN
!     FNO3S=AMAX1(0.0,ZNO3S(L,NY,NX))/ZNO3T(L)
!     FNO3B=AMAX1(0.0,ZNO3B(L,NY,NX))/ZNO3T(L)
!     ELSE
      FNO3S=VLNO3(L,NY,NX)
      FNO3B=VLNOB(L,NY,NX)
!     ENDIF
      ZNO2T(L)=AMAX1(0.0,ZNO2S(L,NY,NX))+AMAX1(0.0,ZNO2B(L,NY,NX))
!     IF(ZNO2T(L).GT.ZEROS(NY,NX))THEN
!     FNO2S=AMAX1(0.0,ZNO2S(L,NY,NX))/ZNO2T(L)
!     FNO2B=AMAX1(0.0,ZNO2B(L,NY,NX))/ZNO2T(L)
!     ELSE
      FNO2S=VLNO3(L,NY,NX)
      FNO2B=VLNOB(L,NY,NX)
!     ENDIF
      H1P4T(L)=AMAX1(0.0,H1PO4(L,NY,NX))+AMAX1(0.0,H1POB(L,NY,NX))
!     IF(H1P4T(L).GT.ZEROS(NY,NX))THEN
!     FH1PS=AMAX1(0.0,H1PO4(L,NY,NX))/H1P4T(L)
!     FH1PB=AMAX1(0.0,H1POB(L,NY,NX))/H1P4T(L)
!     ELSE
      FH1PS=VLPO4(L,NY,NX)
      FH1PB=VLPOB(L,NY,NX)
!     ENDIF
      H2P4T(L)=AMAX1(0.0,H2PO4(L,NY,NX))+AMAX1(0.0,H2POB(L,NY,NX))
!     IF(H2P4T(L).GT.ZEROS(NY,NX))THEN
!     FH2PS=AMAX1(0.0,H2PO4(L,NY,NX))/H2P4T(L)
!     FH2PB=AMAX1(0.0,H2POB(L,NY,NX))/H2P4T(L)
!     ELSE
      FH2PS=VLPO4(L,NY,NX)
      FH2PB=VLPOB(L,NY,NX)
!     ENDIF
!
!     CCO2S=aqueous CO2 concentration
!
      XCO2=CCO2S(L,NY,NX)/(CCO2S(L,NY,NX)+CCKM)
!
!     TOTAL SUBSTRATE
!
!     TOSC=total SOC, TOSA=total colonized SOC
!     TORC=total microbial residue, TOHC=total adsorbed C
!     in each K:
!     OSCT=total SOC n each K, OSAT=total colonized SOC
!     ORCT=total microbial residue, OHCT=total adsorbed C
!
      TOSC=0.0
      TOSA=0.0
      TORC=0.0
      TOHC=0.0
!
!     TOTAL SOLID SUBSTRATE
!
      DO 870 K=0,KL
      OSCT(K)=0.0
      OSAT(K)=0.0
      DO 865 M=1,4
      OSCT(K)=OSCT(K)+OSC(M,K,L,NY,NX)
      OSAT(K)=OSAT(K)+OSA(M,K,L,NY,NX)
865   CONTINUE
      TOSC=TOSC+OSCT(K)
      TOSA=TOSA+OSAT(K)
870   CONTINUE
!
!     TOTAL BIORESIDUE
!
      DO 880 K=0,KL
      ORCT(K)=0.0
      DO 875 M=1,2
      ORCT(K)=ORCT(K)+ORC(M,K,L,NY,NX)
!     IF(L.EQ.4.AND.K.EQ.2)THEN
!     WRITE(*,876)'ORCT',I,J,NX,NY,L,K,M,ORCT(K)
!    2,ORC(M,K,L,NY,NX)
876   FORMAT(A8,7I4,60E12.4)
!     ENDIF
875   CONTINUE
      TORC=TORC+ORCT(K)
!
!     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
!
!     OSRH=total SOC
!
      TOHC=TOHC+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
880   CONTINUE
      DO 860 K=0,KL
      OSRH(K)=OSAT(K)+ORCT(K)+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,861)'OSRH',I,J,NX,NY,L,K,OSRH(K),OSCT(K)
!    2,OSAT(K),ORCT(K),OHC(K,L,NY,NX),OHA(K,L,NY,NX)
861   FORMAT(A8,6I4,20E12.4)
!     ENDIF
860   CONTINUE
      TSRH=TOSA+TORC+TOHC
!
!     C:N AND C:P RATIOS OF TOTAL BIOMASS
!     CNOMA,CPOMA=N,P contents of active biomass OMA
!     FCN,FCP=effects of N,P limitations on biomass activity
!
      TOMA=0.0
      TOMN=0.0
      DO 890 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 895 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      IF(OMC(1,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNOMA(N,K)=AMAX1(0.0,OMN(1,N,K,L,NY,NX)/OMC(1,N,K,L,NY,NX))
      CPOMA(N,K)=AMAX1(0.0,OMP(1,N,K,L,NY,NX)/OMC(1,N,K,L,NY,NX))
      ELSE
      CNOMA(N,K)=CNOMC(1,N,K)
      CPOMA(N,K)=CPOMC(1,N,K)
      ENDIF
      OMA(N,K)=AMAX1(0.0,OMC(1,N,K,L,NY,NX)/FL(1))
      FCN(N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CNOMA(N,K)/CNOMC(1,N,K))))
      FCP(N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CPOMA(N,K)/CPOMC(1,N,K))))
      FCNP(N,K)=AMIN1(FCN(N,K),FCP(N,K))
!
!     TOTAL BIOMASS
!     OMC2=active biomass in recalcitrant fraction
!
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      TOMA=TOMA+OMA(N,K)
      ENDIF
      IF((K.LE.4.AND.N.EQ.2).OR.(K.EQ.5.AND.N.EQ.1))THEN
      TOMN=TOMN+OMA(N,K)
      ENDIF
      OMC2(N,K)=AMAX1(0.0,AMIN1(OMA(N,K)*FL(2),OMC(2,N,K,L,NY,NX)))
      IF(OMC(2,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOM2(N,K)=AMAX1(0.0,OMC2(N,K)/OMC(2,N,K,L,NY,NX))
      OMN2(N,K)=AMAX1(0.0,FOM2(N,K)*OMN(2,N,K,L,NY,NX))
      ELSE
      FOM2(N,K)=0.0
      OMN2(N,K)=0.0
      ENDIF
      ENDIF
895   CONTINUE
      ENDIF
890   CONTINUE
      DO 690 K=0,KL
      TOMK(K)=0.0
      TONK(K)=0.0
      TOPK(K)=0.0
      TONX(K)=0.0
      TOPX(K)=0.0
      DO 685 N=1,7
      TOMK(K)=TOMK(K)+OMA(N,K)
      TONK(K)=TONK(K)+OMA(N,K)*CNOMA(N,K)
      TOPK(K)=TOPK(K)+OMA(N,K)*CPOMA(N,K)
      TONX(K)=TONX(K)+OMA(N,K)*CNOMC(1,N,K)
      TOPX(K)=TOPX(K)+OMA(N,K)*CPOMC(1,N,K)
685   CONTINUE
690   CONTINUE
!
!     FOSRH=fraction of total SOC in each substrate complex K
!
      DO 790 K=0,KL
      IF(TSRH.GT.ZEROS(NY,NX))THEN
      FOSRH(K,L,NY,NX)=OSRH(K)/TSRH
      ELSE
      FOSRH(K,L,NY,NX)=1.0
      ENDIF
!
!     DOC CONCENTRATIONS
!
!     COQC,COQA=aqueous DOC,acetate concentrations
!     VOLWM=soil water content, FOSRH=fraction of total SOC
!     occupied by each substrate complex K
!
      IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(FOSRH(K,L,NY,NX).GT.ZERO)THEN
      COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX) &
      /(VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)))
      COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX) &
      /(VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)))
      ELSE
      COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX)/VOLWM(NPH,L,NY,NX))
      COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX)/VOLWM(NPH,L,NY,NX))
      ENDIF
      ELSE
      COQC(K,L,NY,NX)=0.0
      COQA(K,L,NY,NX)=0.0
      OHCQ=0.0
      ENDIF
!
!     CNQ,CPQ=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
!
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNQ(K)=AMAX1(0.0,OQN(K,L,NY,NX)/OQC(K,L,NY,NX))
      CPQ(K)=AMAX1(0.0,OQP(K,L,NY,NX)/OQC(K,L,NY,NX))
      ELSE
      CNQ(K)=0.0
      CPQ(K)=0.0
      ENDIF
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX).AND.OQA(K,L,NY,NX) &
      .GT.ZEROS(NY,NX))THEN
      FOCA(K)=OQC(K,L,NY,NX)/(OQC(K,L,NY,NX)+OQA(K,L,NY,NX))
      FOAA(K)=1.0-FOCA(K)
      ELSEIF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOCA(K)=1.0
      FOAA(K)=0.0
      ELSE
      FOCA(K)=0.0
      FOAA(K)=1.0
      ENDIF
790   CONTINUE
!
!     nitrous acid concn CHNO2 and energy yield of hydrogenotrophic
!     methanogenesis GH2X at ambient H2 concentration CH2GS
!
      CHY1=AMAX1(ZERO,10.0**(-(PH(L,NY,NX)-3.0)))
      CHNO2=CNO2S(L,NY,NX)*CHY1/0.5
      CHNOB=CNO2B(L,NY,NX)*CHY1/0.5
      GH2X=8.3143E-03*TKS(L,NY,NX) &
      *LOG((AMAX1(1.0E-03,CH2GS(L,NY,NX))/H2KI)**4)
      end subroutine StageBGCEnvironCondition
!------------------------------------------------------------------------------------------

      subroutine MicrobialCatabolism(L,LL,NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: L,LL,NY,NX

      integer :: K,M,N
!     begin_execution

      DO 760 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      TCGOQC(K)=0.0
      TCGOAC(K)=0.0
      TCGOMN(K)=0.0
      TCGOMP(K)=0.0
      DO 750 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      IF(K.LE.4)THEN
      IF(N.EQ.3)THEN
!
!     WFNG=water potential (PSISM) effect on microbial respiration
!     OXKX=Km for O2 uptake
!     OXKM=Km for heterotrophic O2 uptake set in starts.f
!     TFNG=combined temp and water stress effect on growth respiration
!     TFNR=temperature effect on maintenance respiration
!
      WFNG=EXP(0.1*PSISM(L,NY,NX))
      ELSE
      WFNG=EXP(0.2*PSISM(L,NY,NX))
      ENDIF
      OXKX=OXKM
      ELSE
      WFNG=EXP(0.2*PSISM(L,NY,NX))
      OXKX=OXKA
      ENDIF
      TFNG(N,K)=TFNX*WFNG
      TFNR(N,K)=TFNY
!
!     FOMA,FOMN=fraction of total active biomass C,N in each N and K
!
      IF(OMA(N,K).GT.0.0)THEN
      IF(TOMA.GT.ZEROS(NY,NX))THEN
      FOMA(N,K)=OMA(N,K)/TOMA
      ELSE
      FOMA(N,K)=1.0
      ENDIF
      IF(TOMN.GT.ZEROS(NY,NX))THEN
      FOMN(N,K)=OMA(N,K)/TOMN
      ELSE
      FOMN(N,K)=1.0
      ENDIF
      IF(TOMK(K).GT.ZEROS(NY,NX))THEN
      FOMK(N,K)=OMA(N,K)/TOMK(K)
      ELSE
      FOMK(N,K)=1.0
      ENDIF
!
!     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
!
!     COMC=microbial C concentration relative to substrate
!     SPOMK=effect of microbial C concentration on microbial decay
!     RMOMK=effect of microbial C concentration on maintenance respn
!
      IF(ORGCL.GT.ZEROS(NY,NX))THEN
      DO 765 M=1,2
      COMC=OMC(M,N,K,L,NY,NX)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
765   CONTINUE
      ELSE
      DO 770 M=1,2
      SPOMK(M)=1.0
      RMOMK(M)=1.0
770   CONTINUE
      ENDIF
!
!     FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
!     AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
!
      call SubstrateCompetitionFactors(N,K,L,NY,NX)
!
!     HETEROTROPHIC BIOMASS RESPIRATION
!
      IF(K.LE.4)THEN
!
!     RESPIRATION BY HETEROTROPHIC AEROBES:
!     N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!    (6)N2 FIXERS
!
      IF(N.LE.3.OR.N.EQ.6)THEN
      call AerobicHeterotrophCatabolism(N,K,L,NY,NX)
!     RESPIRATION BY HETEROTROPHIC ANAEROBES:
!     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
!
!     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND
!     ACETATE CONCENTRATION
!
!     GH2F=energy yield of acetotrophic methanogenesis per g C
!     GHAX=H2 effect on energy yield of fermentation
!     GOAX=acetate effect on energy yield of fermentation
!     ECHZ=growth respiration efficiency of fermentation
!
      ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
      call AnaerobCatabolism(N,K,L,NY,NX)
!     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
!
!     GOMX=acetate effect on energy yield
!     ECHZ=growth respiration efficiency of aceto. methanogenesis
!
      ELSEIF(N.EQ.5)THEN
      call AcetoMethanogenCatabolism(N,K,L,NY,NX)
      ENDIF
!
!     RESPIRATION RATES BY AUTOTROPHS 'RGOMP' FROM SPECIFIC
!     OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR. N=(1) NH4 OXIDIZERS (2) NO2
!     OXIDIZERS,(3) CH4 OXIDIZERS, (5) H2TROPHIC METHANOGENS
!
      ELSEIF(K.EQ.5)THEN
!
!     NH3 OXIDIZERS
!
      IF(N.EQ.1)THEN
      call NH3OxidizerCatabolism(N,K,L,NY,NX)
!     NO2 OXIDIZERS
!
      ELSEIF(N.EQ.2)THEN
!
      call NO2OxidizerCatabolism(N,K,L,NY,NX)
!     H2TROPHIC METHANOGENS
!
      ELSEIF(N.EQ.5)THEN

      call H2MethanogensCatabolism(N,K,L,NY,NX)
!     METHANOTROPHS
!
      ELSEIF(N.EQ.3)THEN
      call MethanotrophCatabolism(N,K,L,NY,NX)
      ELSE
      RGOMP=0.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ENDIF
      ELSE
      RGOMP=0.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ENDIF
!
!     O2 UPTAKE BY AEROBES
!
!     RUPOX, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
!     RUPMX=O2-unlimited rate of O2 uptake
!     FOXYX=fraction of O2 uptake by N,K relative to total
!     XNPG=1/(NPH*NPT)
!     ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
!     OLSGL=aqueous O2 diffusivity
!     OXYG,OXYS=gaseous, aqueous O2 amounts
!     FLQRQ,FLQRI=surface water flux from precipitation, irrigation
!     COXR,COXQ=O2 concentration in FLQRQ,FLQRI
!
      RUPOX(N,K)=0.0
      IF(N.LE.3.OR.N.EQ.6)THEN
!     N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!    (6)N2 FIXERS
      call AerobsO2Uptake(N,K,L,NY,NX)
      ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
      RGOMO(N,K)=RGOMP
      RCO2X(N,K)=0.333*RGOMO(N,K)
      RCH3X(N,K)=0.667*RGOMO(N,K)
      RCH4X(N,K)=0.0
      ROXYO(N,K)=ROXYM(N,K)
      IF(K.LE.4)THEN
      RH2GX(N,K)=0.111*RGOMO(N,K)
      ELSE
      RH2GX(N,K)=0.0
      ENDIF
      ELSEIF(N.EQ.5)THEN
      RGOMO(N,K)=RGOMP
      IF(K.LE.4)THEN
      RCO2X(N,K)=0.50*RGOMO(N,K)
      RCH3X(N,K)=0.00
      RCH4X(N,K)=0.50*RGOMO(N,K)
      ROXYO(N,K)=ROXYM(N,K)
      RH2GX(N,K)=0.0
      ELSEIF(K.EQ.5)THEN
      RCO2X(N,K)=0.00
      RCH3X(N,K)=0.00
      RCH4X(N,K)=RGOMO(N,K)
      ROXYO(N,K)=ROXYM(N,K)
      RH2GX(N,K)=0.0
      RH2GZ=0.667*RGOMO(N,K)
      ENDIF
      ENDIF
!
!     HETEROTROPHIC DENITRIFICATION
!
      IF(K.LE.4.AND.N.EQ.2.AND.ROXYM(N,K).GT.0.0 &
      .AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN
      call HeteroDenitrificCatabolism(N,K,L,NY,NX)
!     AUTOTROPHIC DENITRIFICATION
!
      ELSEIF(K.EQ.5.AND.N.EQ.1.AND.ROXYM(N,K).GT.0.0 &
      .AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN

      call AutotrophDenitrificCatabolism(N,K,L,NY,NX)
      ELSE
      RDNO3(N,K)=0.0
      RDNOB(N,K)=0.0
      RDNO2(N,K)=0.0
      RDN2B(N,K)=0.0
      RDN2O(N,K)=0.0
      RGOMY(N,K)=0.0
      RGOMD(N,K)=0.0
      ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
      call BiomassMineralization(N,K,L,NY,NX)
!
      call GatherMicrobialRespiration(N,K,L,NY,NX)
!
      call GetMicrobialAnabolismFlux(N,K,L,NY,NX)

      ELSE
      call DealNoActiveMicrobes(N,K,L)
      ENDIF
      ENDIF
750   CONTINUE
      ENDIF
760   CONTINUE
      end subroutine MicrobialCatabolism
!------------------------------------------------------------------------------------------

      subroutine ChemoDenitrification(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX
!     begin_execution
!
!     FNO2,FNB2=fraction of total NO2 demand in non-band,band
!     VMXC4S,VMXC4B=substrate-unlimited NO2 reduction in non-band,band
!     CHNO2,CHNOB=nitrous acid concentration in non-band,band
!     VOLWM=soil water content
!     FNO3S,FNO3B=fractions of NO2 in non-band,band
!     TFNX=temperature stress function
!     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
!     RCN2O,RCN2B=N2O production from nitrous acid reduction in non-band,band
!     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
!     RCOQN=DON production from nitrous acid reduction
!     RVMXC,RVMBC=demand for NO2 reduction in non-band,band
!
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMXC(L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=FMN*VLNO3(L,NY,NX)
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMBC(L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=FMN*VLNOB(L,NY,NX)
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
      VMXC4S=7.5E-02*CHNO2*VOLWM(NPH,L,NY,NX)*FNO3S*TFNX
      VMXC4B=7.5E-02*CHNOB*VOLWM(NPH,L,NY,NX)*FNO3B*TFNX
      RCNO2=AMAX1(0.0,AMIN1(ZNO2S(L,NY,NX)*FNO2,VMXC4S))
      RCNOB=AMAX1(0.0,AMIN1(ZNO2B(L,NY,NX)*FNB2,VMXC4B))
      RCN2O=0.10*RCNO2
      RCN2B=0.10*RCNOB
      RCNO3=0.80*RCNO2
      RCN3B=0.80*RCNOB
      RCOQN=0.10*(RCNO2+RCNOB)
      RVMXC(L,NY,NX)=VMXC4S
      RVMBC(L,NY,NX)=VMXC4B
!     IF((I/1)*1.EQ.I.AND.L.LE.5)THEN
!     WRITE(*,7779)'CHEMO',I,J,L,RCNO2,RCNOB,CHY1,CHNO2,CHNOB
!    2,CNO2S(L,NY,NX),CNO2B(L,NY,NX),VOLWM(NPH,L,NY,NX),FNO2
!    3,VMXC4S,VMXC4B,RVMXC(L,NY,NX),RNO2Y(L,NY,NX),RCN2O,RCN2B
!    4,RCNO3,RCNOB,RCOQN,VLNO3(L,NY,NX),VLNOB(L,NY,NX)
7779  FORMAT(A8,3I4,30E12.4)
!     ENDIF
      end subroutine ChemoDenitrification
!------------------------------------------------------------------------------------------

      subroutine OMTransferForPriming(KL,L,NY,NX)

      implicit none
      integer, intent(in) :: KL,L,NY,NX

      integer :: K,M,N,KK
!     begin_execution
!
!     OSRH=total SOC in each K
!     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
!     between each K and KK, FPRIM=priming transfer rate constant
!     TFND=temperature effect on priming transfers
!     ROQCK,OQC,OQN,OQP=respiration,DOC,DON,DOP
!     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
!
      DO 795 K=0,KL
      IF(K.LE.KL-1)THEN
      DO 800 KK=K+1,KL
      OSRT=OSRH(K)+OSRH(KK)
      IF(OSRH(K).GT.ZEROS(NY,NX).AND.OSRH(KK).GT.ZEROS(NY,NX))THEN
      XFRK=FPRIM*TFND(L,NY,NX)*(ROQCK(K)*OSRH(KK) &
      -ROQCK(KK)*OSRH(K))/OSRT
      XFRC=FPRIM*TFND(L,NY,NX)*(OQC(K,L,NY,NX)*OSRH(KK) &
      -OQC(KK,L,NY,NX)*OSRH(K))/OSRT
      XFRN=FPRIM*TFND(L,NY,NX)*(OQN(K,L,NY,NX)*OSRH(KK) &
      -OQN(KK,L,NY,NX)*OSRH(K))/OSRT
      XFRP=FPRIM*TFND(L,NY,NX)*(OQP(K,L,NY,NX)*OSRH(KK) &
      -OQP(KK,L,NY,NX)*OSRH(K))/OSRT
      XFRA=FPRIM*TFND(L,NY,NX)*(OQA(K,L,NY,NX)*OSRH(KK) &
      -OQA(KK,L,NY,NX)*OSRH(K))/OSRT
      IF(ROQCK(K)+XOQCK(K)-XFRK.GT.0.0 &
      .AND.ROQCK(KK)+XOQCK(KK)+XFRK.GT.0.0)THEN
      XOQCK(K)=XOQCK(K)-XFRK
      XOQCK(KK)=XOQCK(KK)+XFRK
!     IF(I.EQ.116)THEN
!     WRITE(*,4442)'XOQCK',I,J,NX,NY,L,K,KK,XFRC,ROQCK(K)
!    2,OSRH(K),ROQCK(KK),OSRH(KK),XOQCK(K),XOQCK(KK)
!    3,OQC(K,L,NY,NX),OQC(KK,L,NY,NX)
4442  FORMAT(A8,7I4,12E12.4)
!     ENDIF
      ENDIF
      IF(OQC(K,L,NY,NX)+XOQCZ(K)-XFRC.GT.0.0 &
      .AND.OQC(KK,L,NY,NX)+XOQCZ(KK)+XFRC.GT.0.0)THEN
      XOQCZ(K)=XOQCZ(K)-XFRC
      XOQCZ(KK)=XOQCZ(KK)+XFRC
!     IF(I.EQ.116)THEN
!     WRITE(*,4442)'XOQCZ',I,J,NX,NY,L,K,KK,XFRC,OQC(K,L,NY,NX)
!    2,OSRH(K),OQC(KK,L,NY,NX),OSRH(KK),XOQCZ(K),XOQCZ(KK)
!    3,OQC(K,L,NY,NX),OQC(KK,8,NY,NX)
!     ENDIF
      ENDIF
      IF(OQN(K,L,NY,NX)+XOQNZ(K)-XFRN.GT.0.0 &
      .AND.OQN(KK,L,NY,NX)+XOQNZ(KK)+XFRN.GT.0.0)THEN
      XOQNZ(K)=XOQNZ(K)-XFRN
      XOQNZ(KK)=XOQNZ(KK)+XFRN
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
!     WRITE(*,4442)'XOQNZ',I,J,NX,NY,L,K,KK,XFRN,OQN(K,L,NY,NX)
!    2,OSRH(K),OQN(KK,L,NY,NX),OSRH(KK),XOQNZ(K),XOQNZ(KK)
!     ENDIF
      ENDIF
      IF(OQP(K,L,NY,NX)+XOQPZ(K)-XFRP.GT.0.0 &
      .AND.OQP(KK,L,NY,NX)+XOQPZ(KK)+XFRP.GT.0.0)THEN
      XOQPZ(K)=XOQPZ(K)-XFRP
      XOQPZ(KK)=XOQPZ(KK)+XFRP
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
!     WRITE(*,4442)'XOQPZ',I,J,NX,NY,L,K,KK,XFRP,OQP(K,L,NY,NX)
!    2,OSRH(K),OQP(KK,L,NY,NX),OSRH(KK),XOQPZ(K),XOQPZ(KK)
!     ENDIF
      ENDIF
      IF(OQA(K,L,NY,NX)+XOQAZ(K)-XFRA.GT.0.0 &
      .AND.OQA(KK,L,NY,NX)+XOQAZ(KK)+XFRA.GT.0.0)THEN
      XOQAZ(K)=XOQAZ(K)-XFRA
      XOQAZ(KK)=XOQAZ(KK)+XFRA
!     IF((I/1)*1.EQ.I.AND.L.EQ.3.AND.K.EQ.1)THEN
!     WRITE(*,4442)'XOQAZ',I,J,NX,NY,L,K,KK,XFRA,OQA(K,L,NY,NX)
!    2,OSRH(K),OQA(KK,L,NY,NX),OSRH(KK),XOQAZ(K),XOQAZ(KK)
!     ENDIF
      ENDIF
!
!     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
!
!     XFMC,XFMN,XFMP=transfer of microbial C,N,P
!     between each K and KK, FPRIMM=priming transfer rate constant
!     TFNG=temperature+water effect
!     OMC,OMN,OMP=microbial C,N,P
!     OSRH=total SOC in each K
!     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
!
      DO 850 N=1,7
      DO 850 M=1,3
      XFMC=FPRIMM*TFNG(N,K)*(OMC(M,N,K,L,NY,NX)*OSRH(KK) &
      -OMC(M,N,KK,L,NY,NX)*OSRH(K))/OSRT
      XFMN=FPRIMM*TFNG(N,K)*(OMN(M,N,K,L,NY,NX)*OSRH(KK) &
      -OMN(M,N,KK,L,NY,NX)*OSRH(K))/OSRT
      XFMP=FPRIMM*TFNG(N,K)*(OMP(M,N,K,L,NY,NX)*OSRH(KK) &
      -OMP(M,N,KK,L,NY,NX)*OSRH(K))/OSRT
      IF(OMC(M,N,K,L,NY,NX)+XOMCZ(M,N,K)-XFMC.GT.0.0 &
      .AND.OMC(M,N,KK,L,NY,NX)+XOMCZ(M,N,KK)+XFMC.GT.0.0)THEN
      XOMCZ(M,N,K)=XOMCZ(M,N,K)-XFMC
      XOMCZ(M,N,KK)=XOMCZ(M,N,KK)+XFMC
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,4447)'XOMCZ',I,J,NX,NY,L,K,KK,N,M,XFMC
!    2,OMC(M,N,K,L,NY,NX)
!    2,OQC(K,L,NY,NX),OMC(M,N,KK,L,NY,NX),OQC(KK,8,NY,NX)
!    3,XOMCZ(M,N,K),XOMCZ(M,N,KK)
4447  FORMAT(A8,9I4,20E12.4)
!     ENDIF
      ENDIF
      IF(OMN(M,N,K,L,NY,NX)+XOMNZ(M,N,K)-XFMN.GT.0.0 &
      .AND.OMN(M,N,KK,L,NY,NX)+XOMNZ(M,N,KK)+XFMN.GT.0.0)THEN
      XOMNZ(M,N,K)=XOMNZ(M,N,K)-XFMN
      XOMNZ(M,N,KK)=XOMNZ(M,N,KK)+XFMN
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,4447)'XOMNZ',I,J,NX,NY,L,K,KK,N,M,XFMN
!    2,OMN(M,N,K,L,NY,NX)
!    2,OSRH(K),OMN(M,N,KK,L,NY,NX),OSRH(KK),XOMNZ(M,N,K),XOMNZ(M,N,KK)
!     ENDIF
      ENDIF
      IF(OMP(M,N,K,L,NY,NX)+XOMPZ(M,N,K)-XFMP.GT.0.0 &
      .AND.OMP(M,N,KK,L,NY,NX)+XOMPZ(M,N,KK)+XFMP.GT.0.0)THEN
      XOMPZ(M,N,K)=XOMPZ(M,N,K)-XFMP
      XOMPZ(M,N,KK)=XOMPZ(M,N,KK)+XFMP
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,4447)'XOMPZ',I,J,NX,NY,L,K,KK,N,M,XFMP
!    2,OMP(M,N,K,L,NY,NX),OSRH(K),OMP(M,N,KK,L,NY,NX),OSRH(KK)
!    2,XOMPZ(M,N,K),XOMPZ(M,N,KK)
!     ENDIF
      ENDIF
850   CONTINUE
      ENDIF
800   CONTINUE
      ENDIF
795   CONTINUE
!
!     TRANSFER ALL PRIMING AMONG ALL K
!
!     TOQCK=total respiration of DOC+DOA in soil layer
!     ROQCK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
      TOQCK(L,NY,NX)=0.0
      DO 840 K=0,KL
      ROQCK(K)=ROQCK(K)+XOQCK(K)
      TOQCK(L,NY,NX)=TOQCK(L,NY,NX)+ROQCK(K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+XOQCZ(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+XOQNZ(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+XOQPZ(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+XOQAZ(K)
      DO 840 N=1,7
      DO 840 M=1,3
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)+XOMCZ(M,N,K)
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)+XOMNZ(M,N,K)
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)+XOMPZ(M,N,K)
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,5559)'XOM',I,J,NX,NY,L,K,N,M,OMC(M,N,K,L,NY,NX)
!    2,OMN(M,N,K,L,NY,NX),OMP(M,N,K,L,NY,NX)
!    3,XOMCZ(M,N,K),XOMNZ(M,N,K),XOMPZ(M,N,K)
5559  FORMAT(A8,8I4,12E12.4)
!     ENDIF
840   CONTINUE

      end subroutine OMTransferForPriming
!------------------------------------------------------------------------------------------

      subroutine DOMSorption(K,L,NY,NX)
      implicit none
      integer, intent(in) :: K,L,NY,NX
!     begin_execution
!     VOLWM=soil water content, FOSRH=fraction of total SOC
!     AEC,AECX=anion exchange capacity
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     TCGOQC,TCGOMN,TCGOMP,TCGOAC=total uptake of DOC,DON,DOP,acetate
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     TSORP,HSORP=sorption rate constant and coefficient for OHC
!     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
      IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX) &
      .AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      IF(L.EQ.0)THEN
      AECX=0.5E+03
      ELSE
      AECX=AEC(L,NY,NX)
      ENDIF
      OQCX=AMAX1(ZEROS(NY,NX),OQC(K,L,NY,NX)-TCGOQC(K))
      OQNX=AMAX1(ZEROS(NY,NX),OQN(K,L,NY,NX)-TCGOMN(K))
      OQPX=AMAX1(ZEROS(NY,NX),OQP(K,L,NY,NX)-TCGOMP(K))
      OQAX=AMAX1(ZEROS(NY,NX),OQA(K,L,NY,NX)-TCGOAC(K))
      OHCX=AMAX1(ZEROS(NY,NX),OHC(K,L,NY,NX))
      OHNX=AMAX1(ZEROS(NY,NX),OHN(K,L,NY,NX))
      OHPX=AMAX1(ZEROS(NY,NX),OHP(K,L,NY,NX))
      OHAX=AMAX1(ZEROS(NY,NX),OHA(K,L,NY,NX))
      VOLXX=BKVL(L,NY,NX)*AECX*HSORP*FOSRH(K,L,NY,NX)
      VOLXW=VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)
      IF(FOCA(K).GT.ZERO)THEN
      VOLCX=FOCA(K)*VOLXX
      VOLCW=FOCA(K)*VOLXW
      CSORP(K)=TSORP*(OQCX*VOLCX-OHCX*VOLCW)/(VOLCX+VOLCW)
      ELSE
      CSORP(K)=TSORP*(OQCX*VOLXX-OHCX*VOLXW)/(VOLXX+VOLXW)
      ENDIF
      IF(FOAA(K).GT.ZERO)THEN
      VOLAX=FOAA(K)*VOLXX
      VOLAW=FOAA(K)*VOLXW
      CSORPA(K)=TSORP*(OQAX*VOLAX-OHAX*VOLAW)/(VOLAX+VOLAW)
      ELSE
      CSORPA(K)=TSORP*(OQAX*VOLXX-OHAX*VOLXW)/(VOLXX+VOLXW)
      ENDIF
      ZSORP(K)=TSORP*(OQNX*VOLXX-OHNX*VOLXW)/(VOLXX+VOLXW)
      PSORP(K)=TSORP*(OQPX*VOLXX-OHPX*VOLXW)/(VOLXX+VOLXW)
      ELSE
      CSORP(K)=0.0
      CSORPA(K)=0.0
      ZSORP(K)=0.0
      PSORP(K)=0.0
      ENDIF
      end subroutine DOMSorption
!------------------------------------------------------------------------------------------

      subroutine SolidOMDecomposition(K,L,NY,NX)

      implicit none
      integer, intent(in) :: K,L,NY,NX

      integer :: M
!     begin_execution
!     FCPK=N,P limitation to microbial activity in each K
!     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
!     COQCK=aqueous concentration of microbial activity
!     DCKD=Km for decomposition of SOC at current COQCK
!     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
!     DCKI=inhibition of decomposition by microbial concentration
!     OSRH=total SOC
!     COSC=concentration of total SOC
!     BKVL,VOLX=mass, volume of soil layer
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     OQKI=DOC product inhibition constant for decomposition
!
      IF(TOMK(K).GT.ZEROS(NY,NX))THEN
      CNOMX=TONK(K)/TONX(K)
      CPOMX=TOPK(K)/TOPX(K)
      FCNK(K)=AMIN1(1.0,AMAX1(0.50,CNOMX))
      FCPK(K)=AMIN1(1.0,AMAX1(0.50,CPOMX))
      ELSE
      FCNK(K)=1.0
      FCPK(K)=1.0
      ENDIF
!
!     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
!     CONSTANT FOR DECOMPOSITION
!
      IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
      COQCK=AMIN1(0.1E+06,ROQCK(K)/VOLWZ)
      ELSE
      COQCK=0.1E+06
      ENDIF
      IF(L.EQ.0)THEN
      DCKD=DCKM0*(1.0+COQCK/DCKI)
      ELSE
      DCKD=DCKML*(1.0+COQCK/DCKI)
      ENDIF
      IF(OSRH(K).GT.ZEROS(NY,NX))THEN
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      COSC=OSRH(K)/BKVL(L,NY,NX)
      ELSE
      COSC=OSRH(K)/VOLY(L,NY,NX)
      ENDIF
      DFNS=COSC/(COSC+DCKD)
      OQCI=1.0/(1.0+COQC(K,L,NY,NX)/OQKI)
!     IF(L.EQ.0.AND.J.EQ.15)THEN
!     WRITE(*,4242)'COSC',I,J,L,K,DFNS,COSC,COQCK,DCKD,OSRH(K)
!    2,OSAT(K),OSCT(K),ORCT(K),OHC(K,L,NY,NX),BKVL(L,NY,NX),ROQCK(K)
!    3,VOLWZ,VOLWRX(NY,NX),VOLW(0,NY,NX),FCR(NY,NX)
!    4,THETY(L,NY,NX)
4242  FORMAT(A8,4I4,30E12.4)
!     ENDIF
!
!     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES 'RDOS*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     CNS,CPS=N:C,P:C ratios of SOC
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OSA,OSN,OSP=active biomass C,N,P
!     SPOSC=specific decomposition rate constant
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
      DO 785 M=1,4
      IF(OSC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNS(M,K)=AMAX1(0.0,OSN(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
      CPS(M,K)=AMAX1(0.0,OSP(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
      RDOSC(M,K)=AMAX1(0.0,AMIN1(0.5*OSA(M,K,L,NY,NX) &
      ,SPOSC(M,K)*ROQCK(K)*DFNS*OQCI*TFNX*OSA(M,K,L,NY,NX)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      RDOSN(M,K)=AMAX1(0.0,AMIN1(OSN(M,K,L,NY,NX) &
      ,CNS(M,K)*RDOSC(M,K)))/FCNK(K)
      RDOSP(M,K)=AMAX1(0.0,AMIN1(OSP(M,K,L,NY,NX) &
      ,CPS(M,K)*RDOSC(M,K)))/FCPK(K)
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,4444)'RDOSC',I,J,NX,NY,L,K,M,RDOSC(M,K),RDOSN(M,K)
!    2,RDOSP(M,K),CNS(M,K),CPS(M,K),SPOSC(M,K),ROQCK(K),DFNS,TFNX
!    3,OQCI,OSA(M,K,L,NY,NX),OSRH(K),COSC,COQCK,DCKD,VOLWZ
!    4,COXYS(L,NY,NX),TKS(L,NY,NX),PSISM(L,NY,NX),THETW(L,NY,NX)
!    5,WFN(1,K),WFN(3,K),OXYI,VOLW(0,NY,NX),VOLWRX(NY,NX)
!    4,FOSRH(K,L,NY,NX),VOLY(L,NY,NX),ORGC(L,NY,NX),OSC(M,K,L,NY,NX)
!    2,OSN(M,K,L,NY,NX),OSP(M,K,L,NY,NX),TONK(K),TONX(K),FCNK(K)
!    6,FCPK(K),WFN(1,K),WFN(3,K),WFN(4,K),COQC(K,L,NY,NX)
!    7,OQC(K,L,NY,NX),VOLWM(NPH,L,NY,NX)
4444  FORMAT(A8,7I4,50E12.4)
!     ENDIF
      ELSE
      CNS(M,K)=CNOSC(M,K,L,NY,NX)
      CPS(M,K)=CPOSC(M,K,L,NY,NX)
      RDOSC(M,K)=0.0
      RDOSN(M,K)=0.0
      RDOSP(M,K)=0.0
      ENDIF
785   CONTINUE
!
!     HUMIFICATION OF DECOMPOSED RESIDUE LIGNIN WITH PROTEIN,
!     CH2O AND CELLULOSE 'RHOS*' WITH REMAINDER 'RCOS*' TO DOC,DON,DOP
!
!     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP
!     CNRH,CPRH=N:C,P:C in POC
!     EPOC=fraction of RDOSC allocated to POC from hour1.f
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
      IF(K.LE.2)THEN
      RHOSC(4,K)=AMAX1(0.0,AMIN1(RDOSN(4,K)/CNRH(3) &
      ,RDOSP(4,K)/CPRH(3),EPOC(L,NY,NX)*RDOSC(4,K)))
      RHOSCM=0.10*RHOSC(4,K)
      RHOSC(1,K)=AMAX1(0.0,AMIN1(RDOSC(1,K),RDOSN(1,K)/CNRH(3) &
      ,RDOSP(1,K)/CPRH(3),RHOSCM))
      RHOSC(2,K)=AMAX1(0.0,AMIN1(RDOSC(2,K),RDOSN(2,K)/CNRH(3) &
      ,RDOSP(2,K)/CPRH(3),RHOSCM))
      RHOSC(3,K)=AMAX1(0.0,AMIN1(RDOSC(3,K),RDOSN(3,K)/CNRH(3) &
      ,RDOSP(3,K)/CPRH(3),RHOSCM-RHOSC(2,K)))
      DO 805 M=1,4
      RHOSN(M,K)=AMIN1(RDOSN(M,K),RHOSC(M,K)*CNRH(3))
      RHOSP(M,K)=AMIN1(RDOSP(M,K),RHOSC(M,K)*CPRH(3))
      RCOSC(M,K)=RDOSC(M,K)-RHOSC(M,K)
      RCOSN(M,K)=RDOSN(M,K)-RHOSN(M,K)
      RCOSP(M,K)=RDOSP(M,K)-RHOSP(M,K)
805   CONTINUE
      ELSE
      DO 810 M=1,4
      RHOSC(M,K)=0.0
      RHOSN(M,K)=0.0
      RHOSP(M,K)=0.0
      RCOSC(M,K)=RDOSC(M,K)
      RCOSN(M,K)=RDOSN(M,K)
      RCOSP(M,K)=RDOSP(M,K)
810   CONTINUE
      ENDIF
      ELSE
      DO 780 M=1,4
      RDOSC(M,K)=0.0
      RDOSN(M,K)=0.0
      RDOSP(M,K)=0.0
      RHOSC(M,K)=0.0
      RHOSN(M,K)=0.0
      RHOSP(M,K)=0.0
      RCOSC(M,K)=0.0
      RCOSN(M,K)=0.0
      RCOSP(M,K)=0.0
780   CONTINUE
      ENDIF
!
!     C, N, P DECOMPOSITION RATE OF BIORESIDUE 'RDOR*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     CNR,CPR=N:C,P:C ratios of microbial residue
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     SPORC=specific decomposition rate constant for microbial residue
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
      IF(OSRH(K).GT.ZEROS(NY,NX))THEN
      DO 775 M=1,2
      IF(ORC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNR=AMAX1(0.0,ORN(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
      CPR=AMAX1(0.0,ORP(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
      RDORC(M,K)=AMAX1(0.0,AMIN1(ORC(M,K,L,NY,NX) &
      ,SPORC(M)*ROQCK(K)*DFNS*OQCI*TFNX*ORC(M,K,L,NY,NX)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      RDORN(M,K)=AMAX1(0.0,AMIN1(ORN(M,K,L,NY,NX),CNR*RDORC(M,K))) &
      /FCNK(K)
      RDORP(M,K)=AMAX1(0.0,AMIN1(ORP(M,K,L,NY,NX),CPR*RDORC(M,K))) &
      /FCPK(K)
      ELSE
      RDORC(M,K)=0.0
      RDORN(M,K)=0.0
      RDORP(M,K)=0.0
      ENDIF
775   CONTINUE
      ELSE
      DO 776 M=1,2
      RDORC(M,K)=0.0
      RDORN(M,K)=0.0
      RDORP(M,K)=0.0
776   CONTINUE
      ENDIF
!
!     C, N, P DECOMPOSITION RATE OF SORBED SUBSTRATES 'RDOH*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     CNH,CPH=N:C,P:C ratios of adsorbed C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     SPOHC=specific decomposition rate constant for adsorbed C
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
      IF(OSRH(K).GT.ZEROS(NY,NX))THEN
      IF(OHC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNH(K)=AMAX1(0.0,OHN(K,L,NY,NX)/OHC(K,L,NY,NX))
      CPH(K)=AMAX1(0.0,OHP(K,L,NY,NX)/OHC(K,L,NY,NX))
      RDOHC(K)=AMAX1(0.0,AMIN1(OHC(K,L,NY,NX) &
      ,SPOHC*ROQCK(K)*DFNS*OQCI*TFNX*OHC(K,L,NY,NX)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      RDOHN(K)=AMAX1(0.0,AMIN1(OHN(K,L,NY,NX),CNH(K)*RDOHC(K))) &
      /FCNK(K)
      RDOHP(K)=AMAX1(0.0,AMIN1(OHP(K,L,NY,NX),CPH(K)*RDOHC(K))) &
      /FCPK(K)
      RDOHA(K)=AMAX1(0.0,AMIN1(OHA(K,L,NY,NX) &
      ,SPOHA*ROQCK(K)*DFNS*TFNX*OHA(K,L,NY,NX)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      ELSE
      CNH(K)=0.0
      CPH(K)=0.0
      RDOHC(K)=0.0
      RDOHN(K)=0.0
      RDOHP(K)=0.0
      RDOHA(K)=0.0
      ENDIF
      ELSE
      CNH(K)=0.0
      CPH(K)=0.0
      RDOHC(K)=0.0
      RDOHN(K)=0.0
      RDOHP(K)=0.0
      RDOHA(K)=0.0
      ENDIF
      end subroutine SolidOMDecomposition
!------------------------------------------------------------------------------------------

      subroutine RedistDecompositionProduct(KL,L,NY,NX)

      implicit none
      integer, intent(in) :: KL,L,NY,NX

      integer :: K,M,N
!     begin_execution
!
!     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
!     HETEROTROPHIC SUBSTRATE-MICROBE COMPLEXES
!
!     FORC=fraction of total microbial residue
!     ORCT=microbial residue
!     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
      DO 1690 K=0,KL
      IF(TORC.GT.ZEROS(NY,NX))THEN
      FORC(K)=ORCT(K)/TORC
      ELSE
      IF(K.EQ.3)THEN
      FORC(K)=1.0
      ELSE
      FORC(K)=0.0
      ENDIF
      ENDIF
      DO 1685 N=1,7
      DO 1680 M=1,2
      RCCMC(M,N,K)=(RCOMC(M,N,5)+RCMMC(M,N,5))*FORC(K)
      RCCMN(M,N,K)=(RCOMN(M,N,5)+RCMMN(M,N,5))*FORC(K)
      RCCMP(M,N,K)=(RCOMP(M,N,5)+RCMMP(M,N,5))*FORC(K)
!     IF(L.EQ.0)THEN
!     WRITE(*,8821)'RCCMC',I,J,L,K,N,M,RCCMC(M,N,K)
!    2,RCOMC(M,N,5),RCMMC(M,N,5),FORC(K)
!     ENDIF
1680  CONTINUE
1685  CONTINUE
1690  CONTINUE
!
!     REDISTRIBUTE C,N AND P TRANSFORMATIONS AMONG STATE
!     VARIABLES IN SUBSTRATE-MICROBE COMPLEXES
!

      DO 590 K=0,KL
      DO 580 M=1,4
!
!     SUBSTRATE DECOMPOSITION PRODUCTS
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OQC,OQN,OQP,OQA=DOC,DON,DOP
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-RDOSC(M,K)
!     OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-RDOSC(M,K)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-RDOSN(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-RDOSP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RCOSC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RCOSN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RCOSP(M,K)
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.K.EQ.4)THEN
!     WRITE(*,4444)'RDOSC',I,J,NX,NY,L,K,M,OSC(M,K,L,NY,NX)
!    2,RDOSC(M,K)
!     ENDIF
!
!     LIGNIFICATION PRODUCTS
!
!     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!
      IF(L.NE.0)THEN
      OSC(1,3,L,NY,NX)=OSC(1,3,L,NY,NX)+RHOSC(M,K)
!     OSA(1,3,L,NY,NX)=OSA(1,3,L,NY,NX)+RHOSC(M,K)
      OSN(1,3,L,NY,NX)=OSN(1,3,L,NY,NX)+RHOSN(M,K)
      OSP(1,3,L,NY,NX)=OSP(1,3,L,NY,NX)+RHOSP(M,K)
      ELSE
      OSC(1,3,NU(NY,NX),NY,NX)=OSC(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
!     OSA(1,3,NU(NY,NX),NY,NX)=OSA(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
      OSN(1,3,NU(NY,NX),NY,NX)=OSN(1,3,NU(NY,NX),NY,NX)+RHOSN(M,K)
      OSP(1,3,NU(NY,NX),NY,NX)=OSP(1,3,NU(NY,NX),NY,NX)+RHOSP(M,K)
      ENDIF
580   CONTINUE
!
!     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     RCOQN=DON production from nitrous acid reduction
!
      DO 575 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-RDORC(M,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-RDORN(M,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-RDORP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RDORC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RDORN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RDORP(M,K)
575   CONTINUE
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RDOHC(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RDOHN(K)+RCOQN*FORC(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RDOHP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+RDOHA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-RDOHC(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-RDOHN(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-RDOHP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-RDOHA(K)
!
!     MICROBIAL UPTAKE OF DISSOLVED C, N, P
!
!     CGOQC,CGOAC,CGOMN,CGOMP=DOC,acetate,DON,DOP uptake
!     RCH3X=acetate production from fermentation
!
      DO 570 N=1,7
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CGOQC(N,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-CGOMN(N,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-CGOMP(N,K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CGOAC(N,K)+RCH3X(N,K)
!
!     MICROBIAL DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
      DO 565 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+RCOMC(M,N,K)+RCCMC(M,N,K) &
      +RCMMC(M,N,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+RCOMN(M,N,K)+RCCMN(M,N,K) &
      +RCMMN(M,N,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+RCOMP(M,N,K)+RCCMP(M,N,K) &
      +RCMMP(M,N,K)
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4.AND.K.EQ.2)THEN
!     WRITE(*,8821)'ORC',I,J,L,K,N,M,ORC(M,K,L,NY,NX)
!    2,RCOMC(M,N,K),RCCMC(M,N,K),RCMMC(M,N,K),RDORC(M,K)
!     WRITE(*,8821)'ORP',I,J,L,K,N,M,ORP(M,K,L,NY,NX)
!    2,RCOMP(M,N,K),RCCMP(M,N,K),RCMMP(M,N,K),RDORP(M,K)
8821  FORMAT(A8,6I4,40E12.4)
!     ENDIF
565   CONTINUE
570   CONTINUE
!
!     SORPTION PRODUCTS
!
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CSORP(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ZSORP(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-PSORP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CSORPA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+CSORP(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+ZSORP(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+PSORP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+CSORPA(K)
!     IF(L.EQ.1)THEN
!     WRITE(*,592)'OQC',I,J,NX,NY,L,K,OQC(K,L,NY,NX)
!    2,OQA(K,L,NY,NX),(RCOSC(M,K),M=1,4),(RDORC(M,K),M=1,2)
!    3,RDOHC(K),(CGOQC(N,K),N=1,7),CSORP(K),OHC(K,L,NY,NX)
!    4,(WFN(N,K),N=1,7),RDOHA(K),(RCH3X(N,K),N=1,7)
!    3,(CGOAC(N,K),N=1,7),CSORPA(K),OHA(K,L,NY,NX)
!     WRITE(*,592)'OQN',I,J,NX,NY,L,K,OQN(K,L,NY,NX)
!    2,(RCOSN(M,K),M=1,4),(RDORN(M,K),M=1,2),RDOHN(K)
!    2,RCOQN*FORC(K),(CGOMN(N,K),N=1,7),ZSORP(K),OHN(K,L,NY,NX)
592   FORMAT(A8,6I4,80E12.4)
!     ENDIF
590   CONTINUE
      end subroutine RedistDecompositionProduct
!------------------------------------------------------------------------------------------

      subroutine MicrobialAnabolicUpdate(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX

      integer :: K,M,N
!     begin_execution
!
!     OMC,OMN,OMP=microbial C,N,P
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
!     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
!

      DO 550 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 545 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      DO 540 M=1,2
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)+CGOMS(M,N,K) &
      -RXOMC(M,N,K)-RXMMC(M,N,K)
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)+CGONS(M,N,K) &
      -RXOMN(M,N,K)-RXMMN(M,N,K)
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)+CGOPS(M,N,K) &
      -RXOMP(M,N,K)-RXMMP(M,N,K)
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,4488)'OMN2',I,J,NX,NY,L,K,N,M,OMC(M,N,K,L,NY,NX)
!    2,OMN(M,N,K,L,NY,NX),OMP(M,N,K,L,NY,NX)
!     WRITE(*,4488)'RDOMC',I,J,NX,NY,L,K,N,M,CGOMS(M,N,K),CGOQC(N,K)
!    4,CGOAC(N,K),RGOMO(N,K),RGOMD(N,K),RXOMC(M,N,K),RXMMC(M,N,K)
!    3,RMOMC(M,N,K),TFNX,OMGR,OMC(3,N,K,L,NY,NX),WFN(N,K)
!    3,OMC(M,N,K,L,NY,NX),OMA(N,K),TSRH
!    4,RCH3X(N,K),RH2GZ,RH2GX(4,K),FOCA(K),FOAA(K)
!    6,OQA(K,L,NY,NX),OHA(K,L,NY,NX),OQC(K,L,NY,NX),OHC(K,L,NY,NX)
!    7,OMP(M,N,K,L,NY,NX),CGOPS(M,N,K),RDOMP(M,N,K),RDMMP(M,N,K)
!    8,OMP(3,N,K,L,NY,NX),CGOMP(N,K),RIPO4(N,K)
4488  FORMAT(A8,8I4,40E12.4)
!     ENDIF
!
!     HUMIFICATION PRODUCTS
!
!     CFOMC=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!
      IF(L.NE.0)THEN
      OSC(1,4,L,NY,NX)=OSC(1,4,L,NY,NX)+CFOMC(1,L,NY,NX) &
      *(RHOMC(M,N,K)+RHMMC(M,N,K))
!     OSA(1,4,L,NY,NX)=OSA(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
!    2*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(1,4,L,NY,NX)=OSN(1,4,L,NY,NX)+CFOMC(1,L,NY,NX) &
      *(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(1,4,L,NY,NX)=OSP(1,4,L,NY,NX)+CFOMC(1,L,NY,NX) &
      *(RHOMP(M,N,K)+RHMMP(M,N,K))
      OSC(2,4,L,NY,NX)=OSC(2,4,L,NY,NX)+CFOMC(2,L,NY,NX) &
      *(RHOMC(M,N,K)+RHMMC(M,N,K))
!     OSA(2,4,L,NY,NX)=OSA(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
!    2*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(2,4,L,NY,NX)=OSN(2,4,L,NY,NX)+CFOMC(2,L,NY,NX) &
      *(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(2,4,L,NY,NX)=OSP(2,4,L,NY,NX)+CFOMC(2,L,NY,NX) &
      *(RHOMP(M,N,K)+RHMMP(M,N,K))
!     IF((I/10)*10.EQ.I.AND.J.EQ.24)THEN
!     WRITE(*,4445)'RHOMC',I,J,NX,NY,L,K,M,N,OSC(1,4,L,NY,NX)
!    2,OSC(2,4,L,NY,NX),CFOMC(1,L,NY,NX),CFOMC(2,L,NY,NX)
!    3,RHOMC(M,N,K),RHMMC(M,N,K)
4445  FORMAT(A8,8I4,40E12.4)
!     ENDIF
      ELSE
      OSC(1,4,NU(NY,NX),NY,NX)=OSC(1,4,NU(NY,NX),NY,NX) &
      +CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
!     OSA(1,4,NU(NY,NX),NY,NX)=OSA(1,4,NU(NY,NX),NY,NX)
!    2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(1,4,NU(NY,NX),NY,NX)=OSN(1,4,NU(NY,NX),NY,NX) &
      +CFOMC(1,NU(NY,NX),NY,NX)*(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(1,4,NU(NY,NX),NY,NX)=OSP(1,4,NU(NY,NX),NY,NX) &
      +CFOMC(1,NU(NY,NX),NY,NX)*(RHOMP(M,N,K)+RHMMP(M,N,K))
      OSC(2,4,NU(NY,NX),NY,NX)=OSC(2,4,NU(NY,NX),NY,NX) &
      +CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
!     OSA(2,4,NU(NY,NX),NY,NX)=OSA(2,4,NU(NY,NX),NY,NX)
!    2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(2,4,NU(NY,NX),NY,NX)=OSN(2,4,NU(NY,NX),NY,NX) &
      +CFOMC(2,NU(NY,NX),NY,NX)*(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(2,4,NU(NY,NX),NY,NX)=OSP(2,4,NU(NY,NX),NY,NX) &
      +CFOMC(2,NU(NY,NX),NY,NX)*(RHOMP(M,N,K)+RHMMP(M,N,K))
      ENDIF
540   CONTINUE
!
!     INPUTS TO NONSTRUCTURAL POOLS
!
!     CGOMC=total DOC+acetate uptake
!     RGOMO=total respiration
!     RGOMD=respiration for denitrifcation
!     RGN2F=respiration for N2 fixation
!     RCO2X=total CO2 emission
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
!     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!     CGOMN,CGOMP=DON, DOP uptake
!     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
!     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
!     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
!     RINH4R,RINO3R =substrate-limited NH4,NO3 mineraln-immobiln
!     RIPO4R,RIP14R=substrate-limited H2PO4,HPO4 mineraln-immobiln
!
      CGROMC=CGOMC(N,K)-RGOMO(N,K)-RGOMD(N,K)-RGN2F(N,K)
      RCO2X(N,K)=RCO2X(N,K)+RGN2F(N,K)
      DO 555 M=1,2
      OMC(3,N,K,L,NY,NX)=OMC(3,N,K,L,NY,NX)-CGOMS(M,N,K) &
      +R3OMC(M,N,K)
      OMN(3,N,K,L,NY,NX)=OMN(3,N,K,L,NY,NX)-CGONS(M,N,K) &
      +R3OMN(M,N,K)+R3MMN(M,N,K)
      OMP(3,N,K,L,NY,NX)=OMP(3,N,K,L,NY,NX)-CGOPS(M,N,K) &
      +R3OMP(M,N,K)+R3MMP(M,N,K)
      RCO2X(N,K)=RCO2X(N,K)+R3MMC(M,N,K)
555   CONTINUE
      OMC(3,N,K,L,NY,NX)=OMC(3,N,K,L,NY,NX)+CGROMC
      OMN(3,N,K,L,NY,NX)=OMN(3,N,K,L,NY,NX)+CGOMN(N,K) &
      +RINH4(N,K)+RINB4(N,K)+RINO3(N,K)+RINB3(N,K)+RN2FX(N,K)
      OMP(3,N,K,L,NY,NX)=OMP(3,N,K,L,NY,NX)+CGOMP(N,K) &
      +RIPO4(N,K)+RIPOB(N,K)+RIP14(N,K)+RIP1B(N,K)
      IF(L.EQ.0)THEN
      OMN(3,N,K,L,NY,NX)=OMN(3,N,K,L,NY,NX)+RINH4R(N,K)+RINO3R(N,K)
      OMP(3,N,K,L,NY,NX)=OMP(3,N,K,L,NY,NX)+RIPO4R(N,K)+RIP14R(N,K)
      ENDIF
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,5556)'OMC3',I,J,NX,NY,L,K,N,OMC(3,N,K,L,NY,NX)
!    2,OMN(3,N,K,L,NY,NX),OMP(3,N,K,L,NY,NX),OMC(1,N,K,L,NY,NX)
!    2,OMN(1,N,K,L,NY,NX),OMP(1,N,K,L,NY,NX),WFN(N,K),OXYI
!    2,COXYS(L,NY,NX)
!    2,CGOMS(1,N,K),CGOMS(2,N,K),CGROMC
!    3,CGOPS(1,N,K),CGOPS(2,N,K),CGOMP(N,K),RIPO4(N,K)
!    4,CGOMC(N,K),RGOMO(N,K),RGOMD(N,K),RMOMT,WFN(N,K)
!    5,(CGONS(M,N,K),M=1,2),(R3OMN(M,N,K),M=1,2),(R3MMN(M,N,K)
!    6,M=1,2),(XOMCZ(M,N,K),M=1,2)
!    6,CGOMN(N,K),RINH4(N,K),RINB4(N,K),RINO3(N,K),RINB3(N,K)
!    7,RN2FX(N,K)
5556  FORMAT(A8,7I4,60E12.4)
!     ENDIF
      ENDIF
545   CONTINUE
      ENDIF
550   CONTINUE
      end subroutine MicrobialAnabolicUpdate
!------------------------------------------------------------------------------------------

      subroutine MicrobialLitterColonization(KL,L,NY,NX)
      implicit none
      integer, intent(in) :: KL,L,NY,NX

      integer :: K,M
!     begin_execution
!     OSCT,OSAT,OSCX=total,colonized,uncolonized SOC
!     OSA,OSC=colonized,total litter
!     DOSA=rate constant for litter colonization
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
      DO 475 K=0,KL
      OSCT(K)=0.0
      OSAT(K)=0.0
      DO 475 M=1,4
      OSCT(K)=OSCT(K)+OSC(M,K,L,NY,NX)
      OSAT(K)=OSAT(K)+OSA(M,K,L,NY,NX)
475   CONTINUE
      DO 480 K=0,KL
      IF(OSCT(K).GT.ZEROS(NY,NX))THEN
      DOSAK=DOSA(K)*AMAX1(0.0,ROQCK(K))
      DO 485 M=1,4
      OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX) &
      ,OSA(M,K,L,NY,NX)+DOSAK*OSC(M,K,L,NY,NX)/OSCT(K))
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,8822)'OSA',I,J,NX,NY,L,K,M,OSA(M,K,L,NY,NX)
!    2,OSC(M,K,L,NY,NX),DOSA(K),ROQCK(K),DOSAK,OSAT(K),OSCT(K)
8822  FORMAT(A8,7I4,30E12.4)
!     ENDIF
485   CONTINUE
      ELSE
      DO 490 M=1,4
      OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX),OSA(M,K,L,NY,NX))
490   CONTINUE
      ENDIF
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,8823)'OSC',I,J,L,K,((OMC(M,N,K,L,NY,NX),N=1,7),M=1,3)
!    2,(ORC(M,K,L,NY,NX),M=1,2),OQC(K,L,NY,NX),OQCH(K,L,NY,NX)
!    3,OHC(K,L,NY,NX),OQA(K,L,NY,NX),OQAH(K,L,NY,NX),OHA(K,L,NY,NX)
!    4,(OSC(M,K,L,NY,NX),M=1,4)
8823  FORMAT(A8,4I4,100E12.4)
!     ENDIF
480   CONTINUE
      end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

      subroutine AggregateTransformations(L,NY,NX)
      implicit none
      integer, intent(in) :: L,NY,NX

      integer :: K,M,N
!     begin_execution

      TRINH=0.0
      TRINO=0.0
      TRIPO=0.0
      TRIP1=0.0
      TRINB=0.0
      TRIOB=0.0
      TRIPB=0.0
      TRIB1=0.0
      TRGOM=0.0
      TRGOC=0.0
      TRGOD=0.0
      TRGOA=0.0
      TRGOH=0.0
      TUPOX=0.0
      TRDN3=0.0
      TRDNB=0.0
      TRDN2=0.0
      TRD2B=0.0
      TRDNO=0.0
      TRN2F=0.0
      DO 650 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 640 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      TRINH=TRINH+RINH4(N,K)
      TRINO=TRINO+RINO3(N,K)
      TRIPO=TRIPO+RIPO4(N,K)
      TRIP1=TRIP1+RIP14(N,K)
      TRINB=TRINB+RINB4(N,K)
      TRIOB=TRIOB+RINB3(N,K)
      TRIPB=TRIPB+RIPOB(N,K)
      TRIB1=TRIB1+RIP1B(N,K)
      TRN2F=TRN2F+RN2FX(N,K)
      IF(L.EQ.NU(NY,NX))THEN
      TRINH=TRINH+RINH4R(N,K)
      TRINO=TRINO+RINO3R(N,K)
      TRIPO=TRIPO+RIPO4R(N,K)
      TRIP1=TRIP1+RIP14R(N,K)
      ENDIF
!     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
!     WRITE(*,4469)'TRINH',I,J,NX,NY,L,K,N,TRINH,RINH4(N,K),RINH4R(N,K)
!     WRITE(*,4469)'TRIPO',I,J,NX,NY,L,K,N,TRIPO,RIPO4(N,K),RIPO4R(N,K)
!    2,CGOMP(N,K)
4469  FORMAT(A8,7I4,20E12.4)
!     ENDIF
      TRGOM=TRGOM+RCO2X(N,K)
      TRGOC=TRGOC+RCH4X(N,K)
      TRGOD=TRGOD+RGOMD(N,K)
      TUPOX=TUPOX+RUPOX(N,K)
      TRDN3=TRDN3+RDNO3(N,K)
      TRDNB=TRDNB+RDNOB(N,K)
      TRDN2=TRDN2+RDNO2(N,K)
      TRD2B=TRD2B+RDN2B(N,K)
      TRDNO=TRDNO+RDN2O(N,K)
      TRGOH=TRGOH+RH2GX(N,K)
!     IF(IYRC.EQ.2012.AND.I.EQ.151.AND.NX.EQ.1)THEN
!     WRITE(*,3333)'TRGOM',I,J,NX,NY,L,K,N,TRGOM
!    2,RCO2X(N,K),TRGOA,RGOMO(N,K),WFN(N,K),RGOMP
!     WRITE(*,3333)'TUPOX',I,J,NX,NY,L,K,N,TUPOX,RUPOX(N,K)
!     ENDIF
!     IF(J.EQ.12.AND.L.LE.4)THEN
!     WRITE(*,3333)'N2O',I,J,NX,NY,L,K,N,TRDN2,TRD2B,TRDNO
!    2,RDNO2(N,K),RDN2B(N,K),RDN2O(N,K),COXYS(L,NY,NX)
!    3,COXYG(L,NY,NX)
!     WRITE(*,3333)'TRGOH',I,J,NX,NY,L,K,N,TRGOH,RH2GX(N,K)
!    2,RGOMO(N,K)
3333  FORMAT(A8,7I4,20E12.4)
!     ENDIF
      ENDIF
640   CONTINUE
      ENDIF
650   CONTINUE
      DO 645 N=1,7
      IF(N.LE.3.OR.N.EQ.5)THEN
      IF(N.NE.3)THEN
      TRGOA=TRGOA+CGOMC(N,5)
      ENDIF
      ENDIF
645   CONTINUE
!
!     ALLOCATE AGGREGATED TRANSFORMATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2O=net CO2 uptake
!     TRGOA=total CO2 uptake by autotrophs
!     TRGOM total CO2 emission by heterotrophs reducing O2
!     TRGOD=total CO2 emission by denitrifiers reducing NOx
!     RVOXA(3)=CH4 oxidation
!     RCH4O=net CH4 uptake
!     CGOMC=total CH4 uptake by autotrophs
!     TRGOC=total CH4 emission
!     RH2GO=net H2 uptake
!     RH2GZ,TRGOH=total H2 uptake, emission
!     RUPOXO,TUPOX=total O2 uptake
!     RN2G=total N2 production
!     TRDNO=total N2O reduction
!     RN2O=total N2O uptake
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RCN2O,RCN2B=nitrous acid reduction in non-band,band
!
      RCO2O(L,NY,NX)=TRGOA-TRGOM-TRGOD-RVOXA(3)
      RCH4O(L,NY,NX)=RVOXA(3)+CGOMC(3,5)-TRGOC
      RH2GO(L,NY,NX)=RH2GZ-TRGOH
      RUPOXO(L,NY,NX)=TUPOX
      RN2G(L,NY,NX)=-TRDNO
      RN2O(L,NY,NX)=-TRDN2-TRD2B-RCN2O-RCN2B+TRDNO
!     IF(L.EQ.10)THEN
!     WRITE(*,2468)'RCO2O',I,J,NX,NY,L,RCO2O(L,NY,NX)
!    2,TRGOA,TRGOM,TRGOD,RVOXA(3),RCH4O(L,NY,NX)
!    3,CGOMC(3,5),TRGOC
!     WRITE(*,2468)'RN2O',I,J,NX,NY,L
!    2,RN2O(L,NY,NX),TRDN2,TRD2B,RCN2O,RCN2B,TRDNO
!    2,RCH4O(L,NY,NX),RVOXA(3)
!    2,CGOMC(3,5),TRGOC,(OMA(N,1),N=1,7)
2468  FORMAT(A8,5I4,20E12.4)
!     ENDIF
!
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate
!
      DO 655 K=0,4
      DO 660 M=1,4
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RCOSC(M,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RCOSN(M,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RCOSP(M,K)
660   CONTINUE
      DO 665 M=1,2
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RDORC(M,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RDORN(M,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RDORP(M,K)
665   CONTINUE
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RDOHC(K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RDOHN(K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RDOHP(K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)+RDOHA(K)
      DO 670 N=1,7
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CGOQC(N,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-CGOMN(N,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-CGOMP(N,K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CGOAC(N,K)+RCH3X(N,K)
670   CONTINUE
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CSORP(K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-ZSORP(K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-PSORP(K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CSORPA(K)
655   CONTINUE
!
!     XNH4S,XNH4B=net change in NH4 in band,non-band
!     TRINH,TRINB=total NH4 mineraln-immobn in non-band,band
!     RVOXA(1),RVOXB(1)=total NH4 oxidation in non-band,band
!     XNO3S,XNO3B=net change in NO3 in band,non-band
!     TRINO,TRIOB=total NO3 immobn in non-band,band
!     RVOXA(2),RVOXB(2)=total NO2 oxidation in non-band,band
!     TRDN3,TRDNB=total NO3 reduction in non-band,band
!     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
!     XNO2S,XNO2B=net change in NO3 in band,non-band
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
!     XH2PS,XH2BS=net change in H2PO4 in band,non-band
!     TRIPO,TRIPB=total H2PO4 mineraln-immobn in non-band,band
!     XH1PS,XH1BS=net change in HPO4 in band,non-band
!     TRIP1,TRIB1=total HPO4 mineraln-immobn in non-band,band
!     XN2GS=total N2 fixation
!     XZHYS=total H+ production
!     TRN2F=total N2 fixation
!
      XNH4S(L,NY,NX)=-TRINH-RVOXA(1)
      XNO3S(L,NY,NX)=-TRINO+RVOXA(2)-TRDN3+RCNO3
      XNO2S(L,NY,NX)=RVOXA(1)-RVOXA(2)+TRDN3-TRDN2-RCNO2
      XH2PS(L,NY,NX)=-TRIPO
      XH1PS(L,NY,NX)=-TRIP1
      XNH4B(L,NY,NX)=-TRINB-RVOXB(1)
      XNO3B(L,NY,NX)=-TRIOB+RVOXB(2)-TRDNB+RCN3B
      XNO2B(L,NY,NX)=RVOXB(1)-RVOXB(2)+TRDNB-TRD2B-RCNOB
      XH2BS(L,NY,NX)=-TRIPB
      XH1BS(L,NY,NX)=-TRIB1
      XN2GS(L,NY,NX)=TRN2F
      TFNQ(L,NY,NX)=TFNX
      VOLQ(L,NY,NX)=VOLWZ
      end subroutine AggregateTransformations
!------------------------------------------------------------------------------------------

      subroutine SubstrateCompetitionFactors(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!     F*=fraction of substrate uptake relative to total uptake from
!     previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
!     NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
!     POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
!     oxidation, OQA=acetate oxidation
!
      IF(ROXYY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOXYX=AMAX1(FMN,ROXYS(N,K,L,NY,NX)/ROXYY(L,NY,NX))
      ELSE
      FOXYX=AMAX1(FMN,FOMA(N,K))
      ENDIF
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4X=AMAX1(FMN,RINHO(N,K,L,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4X=AMAX1(FMN,FOMA(N,K)*VLNH4(L,NY,NX))
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB4X=AMAX1(FMN,RINHB(N,K,L,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNB4X=AMAX1(FMN,FOMA(N,K)*VLNHB(L,NY,NX))
      ENDIF
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3X=AMAX1(FMN,RINOO(N,K,L,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3X=AMAX1(FMN,FOMA(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB3X=AMAX1(FMN,RINOB(N,K,L,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNB3X=AMAX1(FMN,FOMA(N,K)*VLNOB(L,NY,NX))
      ENDIF
      IF(RPO4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4X=AMAX1(FMN,RIPOO(N,K,L,NY,NX)/RPO4Y(L,NY,NX))
      ELSE
      FPO4X=AMAX1(FMN,FOMA(N,K)*VLPO4(L,NY,NX))
      ENDIF
      IF(RPOBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPOBX=AMAX1(FMN,RIPBO(N,K,L,NY,NX)/RPOBY(L,NY,NX))
      ELSE
      FPOBX=AMAX1(FMN,FOMA(N,K)*VLPOB(L,NY,NX))
      ENDIF
      IF(RP14Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP14X=AMAX1(FMN,RIPO1(N,K,L,NY,NX)/RP14Y(L,NY,NX))
      ELSE
      FP14X=AMAX1(FMN,FOMA(N,K)*VLPO4(L,NY,NX))
      ENDIF
      IF(RP1BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP1BX=AMAX1(FMN,RIPB1(N,K,L,NY,NX)/RP1BY(L,NY,NX))
      ELSE
      FP1BX=AMAX1(FMN,FOMA(N,K)*VLPOB(L,NY,NX))
      ENDIF
      IF(K.LE.4)THEN
      IF(ROQCY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQC=AMAX1(FMN,ROQCS(N,K,L,NY,NX)/ROQCY(K,L,NY,NX))
      ELSE
      FOQC=AMAX1(FMN,FOMK(N,K))
      ENDIF
      TFOQC=TFOQC+FOQC
      IF(ROQAY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQA=AMAX1(FMN,ROQAS(N,K,L,NY,NX)/ROQAY(K,L,NY,NX))
      ELSE
      FOQA=AMAX1(FMN,FOMK(N,K))
      ENDIF
      TFOQA=TFOQA+FOQA
      ENDIF
      TFOXYX=TFOXYX+FOXYX
      TFNH4X=TFNH4X+FNH4X
      TFNO3X=TFNO3X+FNO3X
      TFPO4X=TFPO4X+FPO4X
      TFP14X=TFP14X+FP14X
      TFNH4B=TFNH4B+FNB4X
      TFNO3B=TFNO3B+FNB3X
      TFPO4B=TFPO4B+FPOBX
      TFP14B=TFP14B+FP1BX
!
!     FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
!     MICROBIAL POPULATIONS IN SURFACE RESIDUE
!     F*=fraction of substrate uptake relative to total uptake from
!     previous hour in surface litter, labels as for soil layers above
!
      IF(L.EQ.0)THEN
      IF(RNH4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4XR(N,K)=AMAX1(FMN,RINHOR(N,K,NY,NX) &
      /RNH4Y(NU(NY,NX),NY,NX))
      ELSE
      FNH4XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      IF(RNO3Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3XR(N,K)=AMAX1(FMN,RINOOR(N,K,NY,NX) &
      /RNO3Y(NU(NY,NX),NY,NX))
      ELSE
      FNO3XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      IF(RPO4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4XR(N,K)=AMAX1(FMN,RIPOOR(N,K,NY,NX) &
      /RPO4Y(NU(NY,NX),NY,NX))
      ELSE
      FPO4XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      IF(RP14Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FP14XR(N,K)=AMAX1(FMN,RIPO1R(N,K,NY,NX) &
      /RP14Y(NU(NY,NX),NY,NX))
      ELSE
      FP14XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      ENDIF
      IF(L.EQ.NU(NY,NX).AND.K.NE.3.AND.K.NE.4 &
      .AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
      TFNH4X=TFNH4X+FNH4XR(N,K)
      TFNO3X=TFNO3X+FNO3XR(N,K)
      TFPO4X=TFPO4X+FPO4XR(N,K)
      TFP14X=TFP14X+FP14XR(N,K)
      ENDIF
      end subroutine SubstrateCompetitionFactors
!------------------------------------------------------------------------------------------

      subroutine NH3OxidizerCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!
!     FACTOR TO REGULATE COMPETITION FOR NH4 AMONG DIFFERENT
!     MICROBIAL AND ROOT POPULATIONS FNH4
!
!     FNH4,FNB4=frac of total biol demand for NH4 in non-band, band
!
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4=AMAX1(FMN,RVMX4(N,K,L,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4=AMAX1(FMN,VLNH4(L,NY,NX)*FOMA(N,K))
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB4=AMAX1(FMN,RVMB4(N,K,L,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNB4=AMAX1(FMN,VLNHB(L,NY,NX)*FOMA(N,K))
      ENDIF
      TFNH4X=TFNH4X+FNH4
      TFNH4B=TFNH4B+FNB4
!
!     NITRIFICATION INHIBITION
!
!     ZNFN0=inhibition when fertilizer added
!     ZNFNI=reduction in inhibition since fertilizer added
!     CNH4S,CNH4B=NH4 concentrations in non-band, band
!     TFNX=temperature effect
!     RNFNI=rate constant for inhibition decline
!     ZHKI=inhibition from high CNH4
!     ZNFN4S,ZNFN4B=inhibition in non-band, band
!
      IF(ZNFN0(L,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*(1.0-RNFNI*TFNX)
      ZNFN4S=ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX)/(1.0+CNH4S(L,NY,NX)/ZHKI)
      ZNFN4B=ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX)/(1.0+CNH4B(L,NY,NX)/ZHKI)
      ELSE
      ZNFN4S=1.0
      ZNFN4B=1.0
      ENDIF
!
!     NH3 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     NH3 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXX=potential NH3 oxidation, VMXH=specific oxidation
!     TFNG=temperature+water limitation, FCNP=N,P limitation
!     XCO2=aqueous CO2 limitation, OMA=active biomass
!     VMXA= non-substrate limited NH3 oxidation
!     VHKI=nonlinear increase in VMXA with VMXH
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     CNH4S,CNH4B=NH4 concentration in non-band, band
!     ZHKM=Km for NH4 uptake
!     FNH4,FNB4=fractions of total NH4 demand in non-band, band
!     ZNH4S,ZNH4B=NH4 amount in non-band, band
!     RNNH4,RNNHB=NH3 oxidation in non-band, band
!     RGOMP=O2-unlimited respiration
!     ECNH=efficiency CO2 conversion to biomass
!     RVMX4,RVMXB=nitrifier demand for NH4 in non-band, band
!
      ECHZ=EO2X
      VMXX=VMXH*TFNG(N,K)*FCNP(N,K)*XCO2*OMA(N,K)
      IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
      VMXA=VMXX/(1.0+VMXX/(VHKI*VOLWZ))
      ELSE
      VMXA=0.0
      ENDIF
      FCN4S=FNH4S*CNH4S(L,NY,NX)/(CNH4S(L,NY,NX)+ZHKM)
      FCN4B=FNHBS*CNH4B(L,NY,NX)/(CNH4B(L,NY,NX)+ZHKM)
      FSBST=FCN4S+FCN4B
      VMX4S=VMXA*FCN4S
      VMX4B=VMXA*FCN4B
      RNNH4=AMAX1(0.0,AMIN1(VMX4S,FNH4*ZNH4S(L,NY,NX)))*ZNFN4S
      RNNHB=AMAX1(0.0,AMIN1(VMX4B,FNB4*ZNH4B(L,NY,NX)))*ZNFN4B
      RVOXP=RNNH4+RNNHB
      RVOXPA=RNNH4
      RVOXPB=RNNHB
      RGOMP=AMAX1(0.0,RVOXP*ECNH*ECHZ)
      RVMX4(N,K,L,NY,NX)=VMX4S
      RVMB4(N,K,L,NY,NX)=VMX4B
!
!     O2 DEMAND FROM NH3 OXIDATION
!
!     ROXYM=O2 demand from respiration by nitrifiers
!     ROXYP,ROXYM=O2 demand from respiration + NH3 oxidation
!
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)+3.429*RVOXP
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
!     IF(IYRC.EQ.2012.AND.I.EQ.151.AND.NX.EQ.1)THEN
!     WRITE(*,6666)'NITRI',I,J,L,K,N,RNNH4,RNNHB,VMXX,VMXA,VOLWZ
!    2,CNH4S(L,NY,NX),CNH4B(L,NY,NX)
!    2,14.0*XN4(L,NY,NX),14.0*XNB(L,NY,NX)
!    3,ZNH4S(L,NY,NX),ZNH4B(L,NY,NX),COXYS(L,NY,NX),RGOMP
!    4,PH(L,NY,NX),TFNX,FCNP(N,K),XCO2,ROXYM(N,K)
!    5,VMX4S,VMX4B,FCN4S,FCN4B,FNH4S,FNHBS,OMA(N,K)
!    6,FNH4,FNB4,ZNFN4S,ZNFN4B,ZNFNI(L,NY,NX),ZNFN0(L,NY,NX)
6666  FORMAT(A8,5I4,40E12.4)
!     ENDIF
!
      end subroutine NH3OxidizerCatabolism
!------------------------------------------------------------------------------------------

      subroutine NO2OxidizerCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
!     MICROBIAL POPULATIONS
!
!     FNO2=fraction of total biological demand for NO2 in non-band, band
!
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMN(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMN(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
!
!     NO2 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     NO2 CONCENTRATIONS
!
!     ECHZ=growth respiration efficiency
!     VMXA= non-substrate limited NH3 oxidation
!     VMXN=specific oxidation
!     TFNG=temperature+water limitation, FCNP=N,P limitation
!     XCO2=aqueous CO2 limitation, OMA=active biomass
!     OMA=active biomass
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     CNO2S,CNO2B=NO2 concentration in non-band, band
!     ZNKM=Km for NO2 uptake
!     FNO2,FNB2=fractions of total NO2 demand in non-band, band
!     ZNO2S,ZNO2B=NO2 amount in non-band, band
!     RNNO2,RNNOB=NO2 oxidation in non-band, band
!     RGOMP=O2-unlimited respiration
!     ECNO=efficiency CO2 conversion to biomass
!     RVMX2,RVMB2=nitrifier demand for NO2 in non-band, band
!
      ECHZ=EO2X
      VMXA=TFNG(N,K)*FCNP(N,K)*XCO2*OMA(N,K)*VMXN
      FCN2S=FNH4S*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+ZNKM)
      FCN2B=FNHBS*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+ZNKM)
      FSBST=FCN2S+FCN2B
      VMX2S=VMXA*FCN2S
      VMX2B=VMXA*FCN2B
      RNNO2=AMAX1(0.0,AMIN1(VMX2S,FNO2*ZNO2S(L,NY,NX)))
      RNNOB=AMAX1(0.0,AMIN1(VMX2B,FNB2*ZNO2B(L,NY,NX)))
      RVOXP=RNNO2+RNNOB
      RVOXPA=RNNO2
      RVOXPB=RNNOB
      RGOMP=AMAX1(0.0,RVOXP*ECNO*ECHZ)
      RVMX2(N,K,L,NY,NX)=VMX2S
      RVMB2(N,K,L,NY,NX)=VMX2B
!
!     O2 DEMAND FROM NO2 OXIDATION
!
!     ROXYM=O2 demand from respiration by nitrifiers
!     ROXYP,ROXYM=O2 demand from respiration + NO2 oxidation
!
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)+1.143*RVOXP
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.LE.6)THEN
!     WRITE(*,6667)'NO2OX',I,J,L,K,N,RNNO2,RNNOB,ZNO2S(L,NY,NX)
!    2,ZNO2B(L,NY,NX),CNO2S(L,NY,NX),CNO2B(L,NY,NX),CNH3S(L,NY,NX)
!    3,CNH3B(L,NY,NX),CNH4S(L,NY,NX),CNH4B(L,NY,NX),CNO3S(L,NY,NX)
!    3,CNO3B(L,NY,NX),CHNO2,CHNOB,VMXA,TFNG(N,K),FCNP(N,K),VMXN,ZNKM
!    4,FCN2S,FCN2B,OMA(N,K),FOMN(N,K),TOMN,RVMX2(N,K,L,NY,NX)
!    5,RNO2Y(L,NY,NX),FNO2,FNB2,ROXYM(N,K),ROXYP(N,K)
!    6,ROXYS(N,K,L,NY,NX),VLNHB(L,NY,NX),VLNOB(L,NY,NX)
!    7,SPOMK(1),RMOMK(1)
6667  FORMAT(A8,5I4,50E12.4)
!     ENDIF
!
      end subroutine NO2OxidizerCatabolism
!------------------------------------------------------------------------------------------

      subroutine H2MethanogensCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!
!     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
!
!     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
!     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
!     VMXA=substrate-unlimited H2 oxidation rate
!     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (TRH2G)
!     CH2GS=H2 concentration, H2KM=Km for H2 uptake
!     RGOMP=H2 oxidation, ROXY*=O2 demand
!
      GH2H=GH2X/12.0
      ECHZ=AMAX1(EO2X,AMIN1(1.0 &
      ,1.0/(1.0+AMAX1(0.0,(GCOX+GH2H))/EOMH)))
      VMXA=TFNG(N,K)*FCNP(N,K)*XCO2*OMA(N,K)*VMXC
      H2GSX=H2GS(L,NY,NX)+0.111*TRH2G
      FSBST=CH2GS(L,NY,NX)/(CH2GS(L,NY,NX)+H2KM)
      RGOMP=AMAX1(0.0,AMIN1(1.5*H2GSX,VMXA*FSBST))
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      TCH4A=TCH4A+RGOMP
!     IF((I/30)*30.EQ.I.AND.NX.EQ.3.AND.NY.EQ.1.AND.J.EQ.24)THEN
!     WRITE(*,5553)'H2METH',I,J,NX,NY,L,K,N,RGOMP,H2GS(L,NY,NX)
!    2,H2GSX,CH2GS(L,NY,NX),VMXA,TFNG(N,K),FCNP(N,K),XCO2
!    3,OMA(N,K),VMXC,ECHZ,GCOX,GH2H,TKS(L,NY,NX),FSBST
!    4,SPOMK(1),RMOMK(1)
5553  FORMAT(A8,7I4,20E12.4)
!     ENDIF
!
      end subroutine H2MethanogensCatabolism
!------------------------------------------------------------------------------------------

      subroutine MethanotrophCatabolism(N,K,L,NY,NX)

      implicit none
      integer, intent(in) :: N,K,L,NY,NX

      integer :: M,MM
!     begin_execution
!
!     CH4 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     CH4 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXA=potential oxidation
!     TFNG=temperature+water effect,FCNP=N,P limitation
!     OMA=active biomass,VMX4=specific respiration rate
!     RCH4L=total aqueous CH4 exchange from previous hour
!     RCH4F=total gaseous CH4 exchange from previous hour
!     TCH4H+TCH4A=total CH4 generated from methanogenesis
!     XNPG=1.0/(NPH*NPT)
!     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
!     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
!     VOLPM,VOLWM=air,water-filled porosity
!     SCH4L=CH4 aqueous solubility
!     CCK4=Km for CH4 uptake
!     ECHO=efficiency CO2 conversion to biomass
!     RGOMP1=substrate-limited CH4 oxidation
!     RCHDF=gaseous-aqueous CH4 exchange
!     DFGS=rate constant for gaseous-aqueous exchange
!
      ECHZ=EH4X
      VMXA=TFNG(N,K)*FCNP(N,K)*OMA(N,K)*VMX4
      RCH4L1=RCH4L(L,NY,NX)*XNPG
      RCH4F1=RCH4F(L,NY,NX)*XNPG
      RCH4S1=(TCH4H+TCH4A)*XNPG
      IF(L.EQ.0)THEN
      CH4G1=CCH4E(NY,NX)*VOLPM(1,L,NY,NX)
      ELSE
      CH4G1=CCH4G(L,NY,NX)*VOLPM(1,L,NY,NX)
      ENDIF
      CH4S1=CH4S(L,NY,NX)
      VMXA1=VMXA*XNPG
      RVOXP=0.0
      RGOMP=0.0
!
!     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
!     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
!
      DO 320 M=1,NPH
      IF(VOLWM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLWCH=VOLWM(M,L,NY,NX)*SCH4L(L,NY,NX)
      VOLWPM=VOLWCH+VOLPM(M,L,NY,NX)
      DO 325 MM=1,NPT
      CH4G1=CH4G1+RCH4F1
      CH4S1=CH4S1+RCH4L1+RCH4S1
      CCH4S1=AMAX1(0.0,safe_adb(CH4S1,VOLWM(M,L,NY,NX)))
      FSBST=CCH4S1/(CCH4S1+CCK4)
      RVOXP1=AMIN1(AMAX1(0.0,CH4S1)/(1.0+ECHO*ECHZ) &
      ,VMXA1*FSBST)
      RGOMP1=RVOXP1*ECHO*ECHZ
      CH4S1=CH4S1-RVOXP1-RGOMP1
      IF(THETPM(M,L,NY,NX).GT.THETX)THEN
      RCHDF=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),CH4G1)*VOLWCH &
      -CH4S1*VOLPM(M,L,NY,NX))/VOLWPM
      ELSE
      RCHDF=0.0
      ENDIF
      CH4G1=CH4G1-RCHDF
      CH4S1=CH4S1+RCHDF
      RVOXP=RVOXP+RVOXP1
      RGOMP=RGOMP+RGOMP1
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.0
!    2.AND.MM.EQ.NPT)THEN
!     WRITE(*,5547)'CH4OX',I,J,NX,NY,L,K,N,M,MM,RVOXP1,RGOMP1,CH4G1
!     2,CH4S1,VMXA1,RVOXP,RGOMP,RCHDF,RCH4L1,RCH4F1,RCH4S1,CCH4S1
!    3,ECHO,ECHZ,OMA(N,K),VOLWM(M,L,NY,NX),VOLPM(M,L,NY,NX),VOLWCH
!    4,THETPM(M,L,NY,NX),SCH4L(L,NY,NX),DFGS(M,L,NY,NX)
!    5,COXYS(L,NY,NX),CCH4E(NY,NX),FSBST,SPOMK(1),RMOMK(1)
!    6,CH4G1/VOLPM(M,L,NY,NX)
5547  FORMAT(A8,9I4,30E12.4)
!     ENDIF
325   CONTINUE
      ENDIF
320   CONTINUE
      RVOXPA=RVOXP
      RVOXPB=0.0
!
!     O2 DEMAND FROM CH4 OXIDATION
!
!     ROXYM=O2 demand from respiration
!     ROXYP=O2 demand from respiration + CH4 oxidation
!
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)+4.00*RVOXP
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
      end subroutine MethanotrophCatabolism
!------------------------------------------------------------------------------------------

      subroutine AcetoMethanogenCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
      GOMX=8.3143E-03*TKS(L,NY,NX) &
      *LOG((AMAX1(ZERO,COQA(K,L,NY,NX))/OAKI))
      GOMM=GOMX/24.0
      ECHZ=AMAX1(EO2X,AMIN1(1.0 &
      ,1.0/(1.0+AMAX1(0.0,(GC4X+GOMM))/EOMH)))
!
!     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     COQA=DOA concentration
!     OQKAM=Km for acetate uptake,FCNP=N,P limitation
!     VMXM=specific respiration rate
!     WFNG=water stress effect, OMA=active biomass
!     TFNX=temp stress effect, FOQA= acetate limitation
!     RGOGX=substrate-limited respiration of acetate
!     RGOGX=competition-limited respiration of acetate
!     OQA=acetate, FOQA=fraction of biological demand for acetate
!     RGOMP=O2-unlimited respiration of acetate
!     ROXY*=O2 demand, ROQCS,ROQCA=DOC, acetate demand
!     ROQCD=microbial respiration used to represent microbial activity
!
      FSBST=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKAM)
      RGOGY=AMAX1(0.0,FCNP(N,K)*VMXM*WFNG*OMA(N,K))
      RGOGZ=RGOGY*FSBST*TFNX
      RGOGX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*ECHZ)
      RGOMP=AMIN1(RGOGX,RGOGZ)
      FGOCP=0.0
      FGOAP=1.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ROQCS(N,K,L,NY,NX)=0.0
      ROQAS(N,K,L,NY,NX)=RGOGZ
      ROQCD(N,K)=0.0
      TCH4H=TCH4H+0.5*RGOMP
!     IF((I/30)*30.EQ.I.AND.NX.EQ.3.AND.NY.EQ.1.AND.J.EQ.24)THEN
!     WRITE(*,5552)'ACMETH',I,J,NX,NY,L,K,N,RGOMP,RGOGZ,RGOGX,GOMM
!    2,ECHZ,FCNP(N,K),TFNG(N,K),OMA(N,K),FOQA,COQA(K,L,NY,NX)
!    2,OQA(K,L,NY,NX)
!    3,OMC(1,N,K,L,NY,NX),OMC(2,N,K,L,NY,NX),OMC(3,N,K,L,NY,NX)
!    3,OMN(1,N,K,L,NY,NX),OMN(2,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX)
!    5,VOLWM(NPH,L,NY,NX),PSISM(L,NY,NX),WFNG,COXYS(L,NY,NX)
!    6,OHA(K,L,NY,NX),FSBST,SPOMK(1),RMOMK(1)
5552  FORMAT(A8,7I4,40E12.4)
!     ENDIF
      end subroutine AcetoMethanogenCatabolism
!------------------------------------------------------------------------------------------

      subroutine AerobicHeterotrophCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!     ENERGY YIELDS OF O2 REDOX REACTIONS
!     E* = growth respiration efficiency calculated in PARAMETERS
!
      IF(N.EQ.1)THEN
      EO2Q=EO2X
      ELSEIF(N.EQ.2)THEN
      EO2Q=EO2D
      ELSEIF(N.EQ.3)THEN
      EO2Q=EO2G
      ELSEIF(N.EQ.6)THEN
      EO2Q=ENFX
      ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     'RGO*Z'FROM SPECIFIC RESPIRATION RATE, ACTIVE BIOMASS, DOC OR
!     ACETATE CONCENTRATION,MICROBIAL C:N:P FACTOR, AND TEMPERATURE
!     FOLLOWED BY POTENTIAL RESPIRATION RATES 'RGO*P' WITH UNLIMITED
!     SUBSTRATE USED FOR MICROBIAL COMPETITION FACTOR
!
!     COQC,COQA=DOC,DOA concentration, FOCA,FOAA=DOC,DOA vs DOC+DOA
!     FCNP=N,P limitation,VMXO=specific respiration rate
!     WFNG=water stress effect, OMA=active biomass
!     TFNX=temp stress effect,FOQC,FOQA=OQC,OQA limitation
!     RGOMP=O2-unlimited respiration of DOC+DOA
!     RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
!
      FSBSTC=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)
      FSBSTA=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKA)
      FSBST=FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
      RGOCY=AMAX1(0.0,FCNP(N,K)*VMXO*WFNG*OMA(N,K))
      RGOCZ=RGOCY*FSBSTC*FOCA(K)*TFNX
      RGOAZ=RGOCY*FSBSTA*FOAA(K)*TFNX
      RGOCX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*EO2Q)
      RGOAX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*EO2A)
      RGOCP=AMIN1(RGOCX,RGOCZ)
      RGOAP=AMIN1(RGOAX,RGOAZ)
      RGOMP=RGOCP+RGOAP
      IF(RGOMP.GT.ZEROS(NY,NX))THEN
      FGOCP=RGOCP/RGOMP
      FGOAP=RGOAP/RGOMP
      ELSE
      FGOCP=1.0
      FGOAP=0.0
      ENDIF
!
!     ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
!     BY HETEROTROPHIC AEROBES
!
!     ECHZ=growth respiration yield
!     ROXYM,ROXYP,ROXYS=O2 demand from DOC,DOA oxidation
!     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
!     ROQCD=microbial respiration used to represent microbial activity
!
      ECHZ=EO2Q*FGOCP+EO2A*FGOAP
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)
      ROXYSX=ROXYS(N,K,L,NY,NX)
      ROQCSX=ROQCS(N,K,L,NY,NX)
      ROQASX=ROQAS(N,K,L,NY,NX)
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
      ROQCS(N,K,L,NY,NX)=RGOCZ
      ROQAS(N,K,L,NY,NX)=RGOAZ
      ROQCD(N,K)=RGOCY
!     IF((I/1)*1.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,5555)'RGOMP',I,J,NX,NY,L,K,N,RGOMP,RGOCX,RGOAX,RGOCZ
!    2,RGOAZ,RGOCX,RGOAX,FCNP(N,K),TFNG(N,K),VMXO,OMA(N,K),OSRH(K)
!    2,FOQC,FOQA,COQC(K,L,NY,NX),OQC(K,L,NY,NX),EO2Q,TKS(L,NY,NX)
!    3,COXYS(L,NY,NX),OQKM,OMC(1,N,K,L,NY,NX),OMC(2,N,K,L,NY,NX)
!    4,OMC(3,N,K,L,NY,NX),VOLWM(NPH,L,NY,NX),FOSRH(K,L,NY,NX)
!    5,FSBST,SPOMK(1),RMOMK(1),ROQCD(N,K),ROXYSX,ROXYS(N,K,L,NY,NX)
!    6,ROQCSX,ROQCS(N,K,L,NY,NX),ROQASX,ROQAS(N,K,L,NY,NX)
!    7,TFNX,WFNG,PSISM(L,NY,NX),TKS(L,NY,NX),TKSO
5555  FORMAT(A8,7I4,60E12.4)
!     ENDIF
!
      end subroutine AerobicHeterotrophCatabolism
!------------------------------------------------------------------------------------------

      subroutine AnaerobCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
      GH2F=GH2X/72.0
      GOAX=8.3143E-03*TKS(L,NY,NX) &
      *LOG((AMAX1(ZERO,COQA(K,L,NY,NX))/OAKI)**2)
      GOAF=GOAX/72.0
      GHAX=GH2F+GOAF
      IF(N.EQ.4)THEN
      ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0 &
      /(1.0+AMAX1(0.0,(GCHX-GHAX))/EOMF)))
      ELSE
      ECHZ=AMAX1(ENFX,AMIN1(1.0,1.0 &
      /(1.0+AMAX1(0.0,(GCHX-GHAX))/EOMN)))
      ENDIF
!
!     RESPIRATION RATES BY HETEROTROPHIC ANAEROBES 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     OXYI=O2 inhibition of fermentation
!     FCNP=N,P limitation on respiration
!     VMXF=maximum respiration rate by fermenters
!     WFNG=water stress effect on respiration
!     OMA=active fermenter biomass
!     TFNX=temp stress effect, FOQC=OQC limitation
!     RFOMP=O2-unlimited respiration of DOC
!     ROQCD=microbial respiration used to represent microbial activity
!
      FSBST=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)*OXYI
      RGOFY=AMAX1(0.0,FCNP(N,K)*VMXF*WFNG*OMA(N,K))
      RGOFZ=RGOFY*FSBST*TFNX
      RGOFX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*ECHZ)
      RGOMP=AMIN1(RGOFX,RGOFZ)
      FGOCP=1.0
      FGOAP=0.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ROQCS(N,K,L,NY,NX)=RGOFZ
      ROQAS(N,K,L,NY,NX)=0.0
      ROQCD(N,K)=RGOFY
      TRH2G=TRH2G+RGOMP
!     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.6)THEN
!     WRITE(*,5554)'FERM',I,J,NX,NY,L,K,N,RGOMP,RGOFX,RGOFZ,GHAX,GOAF
!    2,ECHZ,FCNP(N,K),TFNG(N,K),OMA(N,K),OSRH(K),FOQC,COQC(K,L,NY,NX)
!    3,OQKM,OMC(1,N,K,L,NY,NX),OMC(2,N,K,L,NY,NX),OMC(3,N,K,L,NY,NX)
!    3,OMN(1,N,K,L,NY,NX),OMN(2,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX)
!    5,VOLWM(NPH,L,NY,NX),PSISM(L,NY,NX),WFNG,COXYS(L,NY,NX),OXYI
!    6,FSBST,FOSRH(K,L,NY,NX),SPOMK(1),RMOMK(1),ROQCD(N,K)
5554  FORMAT(A8,7I4,60E12.4)
!     ENDIF
!
      end subroutine AnaerobCatabolism
!------------------------------------------------------------------------------------------

      subroutine HeteroDenitrificCatabolism(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!
!     FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
!     AND ROOT POPULATIONS
!
!     FNO3,FNB3=fraction of total biological demand for NO3
!
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3=AMAX1(FMN,RVMX3(N,K,L,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3=AMAX1(FMN,FOMA(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB3=AMAX1(FMN,RVMB3(N,K,L,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNB3=AMAX1(FMN,FOMA(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO3X=TFNO3X+FNO3
      TFNO3B=TFNO3B+FNB3
!
!     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
!
!     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
!     VMXD3=demand for NO3-N reduction
!     VMXDXS,VMXDXB=maximum NO3 reduction in non-band, band
!     FNO3S,FNO3B=fractions of total NO3 in non-band, band
!     CNO3S,CNO3B=NO3 concentrations in non-band, band
!     Z3KM,Z2KM=Km for NO3, NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD3S,VMXD3B=substrate-unlimited NO3 reduction in non-band,band
!     OQCD3S,OQCD3B=DOC limitation to NO3 reduction in non-band, band
!     RDNO3,RDNOB=substrate-limited NO3 reduction in non-band,band
!     RGOM3X,RGOMD3=substrate-unltd,-ltd respn from NO3 reduction
!     RVMX3,RVMB3=demand for NO3 reduction in non-band,band
!
      ROXYD=AMAX1(0.0,ROXYM(N,K)-ROXYO(N,K))
      VMXD3=0.875*ROXYD
      IF(CNO3S(L,NY,NX).GT.ZERO)THEN
      VMXDXS=FNO3S*VMXD3*CNO3S(L,NY,NX)/(CNO3S(L,NY,NX)+Z3KM) &
      /(1.0+(CNO2S(L,NY,NX)*Z3KM)/(CNO3S(L,NY,NX)*Z2KM))
      ELSE
      VMXDXS=0.0
      ENDIF
      IF(CNO3B(L,NY,NX).GT.ZERO)THEN
      VMXDXB=FNO3B*VMXD3*CNO3B(L,NY,NX)/(CNO3B(L,NY,NX)+Z3KM) &
      /(1.0+(CNO2B(L,NY,NX)*Z3KM)/(CNO3B(L,NY,NX)*Z2KM))
      ELSE
      VMXDXB=0.0
      ENDIF
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD3S=VMXDXS*FVMXDX
      VMXD3B=VMXDXB*FVMXDX
      OQCZ3=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC-RGOCP*WFN(N,K))
      OQCD3=OQCZ3/ECN3
      OQCD3S=OQCD3*FNO3S
      OQCD3B=OQCD3*FNO3B
      ZNO3SX=ZNO3S(L,NY,NX)*FNO3
      ZNO3BX=ZNO3B(L,NY,NX)*FNB3
      RDNO3X=AMAX1(0.0,AMIN1(ZNO3SX,VMXD3S))
      RDNOBX=AMAX1(0.0,AMIN1(ZNO3BX,VMXD3B))
      RDNO3(N,K)=AMAX1(0.0,AMIN1(VMXD3S,OQCD3S,ZNO3SX))
      RDNOB(N,K)=AMAX1(0.0,AMIN1(VMXD3B,OQCD3B,ZNO3BX))
      RDNOX=RDNO3X+RDNOBX
      RDNOT=RDNO3(N,K)+RDNOB(N,K)
      RGOM3X=ECN3*RDNOX
      RGOMD3=ECN3*RDNOT
      RVMX3(N,K,L,NY,NX)=VMXD3S
      RVMB3(N,K,L,NY,NX)=VMXD3B
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMA(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMA(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 AND NO3 IN BAND AND NON-BAND SOIL ZONES
!
!     VMXD2=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM,Z1KM=Km for NO2, N2O uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD2S,VMXD2B=substrate-unlimited NO2 reduction in non-band,band
!     OQCD2S,OQCD2B=DOC limitation to NO2 reduction in non-band, band
!     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
!     RGOM2X,RGOMD2=substrate-unltd,-ltd respn from NO2 reduction
!
      VMXD2=VMXD3-RDNOT
      IF(CNO2S(L,NY,NX).GT.ZERO)THEN
      VMXDXS=FNO2S*VMXD2*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM) &
      /(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2S(L,NY,NX)*Z1KM))
      ELSE
      VMXDXS=0.0
      ENDIF
      IF(CNO2B(L,NY,NX).GT.ZERO)THEN
      VMXDXB=FNO2B*VMXD2*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM) &
      /(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2B(L,NY,NX)*Z1KM))
      ELSE
      VMXDXB=0.0
      ENDIF
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD2S=VMXDXS*FVMXDX
      VMXD2B=VMXDXB*FVMXDX
      OQCZ2=AMAX1(0.0,OQCZ3-RGOMD3)
      OQCD2=OQCZ2/ECN2
      OQCD2S=OQCD2*FNO3S
      OQCD2B=OQCD2*FNO3B
      ZNO2SX=(ZNO2S(L,NY,NX)+RDNO3(N,K))*FNO2
      ZNO2BX=(ZNO2B(L,NY,NX)+RDNOB(N,K))*FNB2
      RDNO2X=AMAX1(0.0,AMIN1(ZNO2SX,VMXD2S))
      RDNOBX=AMAX1(0.0,AMIN1(ZNO2BX,VMXD2B))
      RDNO2(N,K)=AMAX1(0.0,AMIN1(VMXD2S,OQCD2S,ZNO2SX))
      RDN2B(N,K)=AMAX1(0.0,AMIN1(VMXD2B,OQCD2B,ZNO2BX))
      RDN2X=RDNO2X+RDNOBX
      RDN2T=RDNO2(N,K)+RDN2B(N,K)
      RGOM2X=ECN2*RDN2X
      RGOMD2=ECN2*RDN2T
      RVMX2(N,K,L,NY,NX)=VMXD2S
      RVMB2(N,K,L,NY,NX)=VMXD2B
!
!     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
!     AND ROOT POPULATIONS
!
!     FN2O=fraction of total biological demand for N2O
!
      IF(RN2OY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FN2O=AMAX1(FMN,RVMX1(N,K,L,NY,NX)/RN2OY(L,NY,NX))
      ELSE
      FN2O=AMAX1(FMN,FOMA(N,K))
      ENDIF
      TFN2OX=TFN2OX+FN2O
!
!     N2O REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS N2O
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2, NO3 AND NO2 IN BAND AND NON-BAND SOIL ZONES
!
!     VMXD1=demand for N2O-N reduction
!     VMXDXS=maximum N2O reduction
!     CZ2OS=N2O concentrations
!     Z1KM=Km for N2O uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD1S=substrate-unlimited N2O reduction
!     OQCD1=DOC limitation to N2O reduction
!     RDN2O=substrate-limited N2O reduction
!     RGOM1X,RGOMD1=substrate-unltd,-ltd  respn from N2O reduction
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NOx reduction
!     RVMX1=demand for N2O reduction
!
      VMXD1=(VMXD2-RDN2T)*2.0
      VMXDXS=VMXD1*CZ2OS(L,NY,NX)/(CZ2OS(L,NY,NX)+Z1KM)
      IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXS/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD1S=VMXDXS*FVMXDX
      OQCZ1=AMAX1(0.0,OQCZ2-RGOMD2)
      OQCD1=OQCZ1/ECN1
      Z2OSX=(Z2OS(L,NY,NX)+RDN2T)*FN2O
      RDN2OX=AMAX1(0.0,AMIN1(Z2OSX,VMXD1S))
      RDN2O(N,K)=AMAX1(0.0,AMIN1(VMXD1S,OQCD1,Z2OSX))
      RGOM1X=ECN1*RDN2OX
      RGOMD1=ECN1*RDN2O(N,K)
      RGOMY(N,K)=RGOM3X+RGOM2X+RGOM1X
      RGOMD(N,K)=RGOMD3+RGOMD2+RGOMD1
      RVMX1(N,K,L,NY,NX)=VMXD1S
!     TRN2OD(NY,NX)=TRN2OD(NY,NX)+RDNO2(N,K)+RDN2B(N,K)
!     TRN2GD(NY,NX)=TRN2GD(NY,NX)+RDN2O(N,K)
!     IF((I/1)*1.EQ.I.AND.L.LE.5)THEN
!     WRITE(*,2222)'DENIT',I,J,L,K,N,RDNO3(N,K),RDNOB(N,K),RDNO2(N,K)
!    2,RDN2B(N,K),RDN2O(N,K),TRN2OD(NY,NX),TRN2GD(NY,NX)
!    3,COXYS(L,NY,NX),COXYG(L,NY,NX),ROXYM(N,K)
!    3,ROXYO(N,K),OMA(N,K),VMXD,CNO3S(L,NY,NX),CNO3B(L,NY,NX)
!    4,CNO2S(L,NY,NX),CNO2B(L,NY,NX),CZ2OS(L,NY,NX),VLNO3(L,NY,NX)
!    5,VLNOB(L,NY,NX),THETW(L,NY,NX),THETI(L,NY,NX),FOMA(N,K)
!    5,ZNO3S(L,NY,NX),ZNO3B(L,NY,NX),ZNO2S(L,NY,NX),ZNO2B(L,NY,NX)
!    6,Z2OS(L,NY,NX),RGOMY(N,K),RGOMD(N,K),TOMA,FOXYX,FNO23S,FNO23B
!    7,OQC(K,L,NY,NX),FOQC,RGOCP,WFN(N,K),VOLWZ,FOSRH(K,L,NY,NX),ZERO
!    9,RGOM3X,RGOM2X,RGOM1X,FNO3,FNO2,FN2O,ZNO3SX,ZNO2SX,Z2OSX
!    3,OQCD3S,OQCD2S,OQCD1,VMXD3S,VMXD2S,VMXD1S,VMXD3,VMXD2,VMXD1
!    4,ROXYD,VMXDX,TFNX,WFNG,TFNG(N,K),PSISM(L,NY,NX)
!    2,(1.0+(CNO2S(L,NY,NX)*Z3KM)/(CNO3S(L,NY,NX)*Z2KM))
!    2,(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2S(L,NY,NX)*Z1KM))
2222  FORMAT(A8,5I4,70E12.4)
!     ENDIF
!
      end subroutine HeteroDenitrificCatabolism
!------------------------------------------------------------------------------------------

      subroutine AutotrophDenitrificCatabolism(N,K,L,NY,NX)

      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMN(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMN(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2
!
!     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
!     VMXD4=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM=Km for NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD4S,VMXD4B=substrate-unlimited NO2 reduction in non-band,band
!     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NO2 reduction
!     ECNO=efficiency CO2 conversion to biomass
!     ECHZ=growth respiration efficiency
!     RVOXA,RVOXB=total O2-limited (1)NH4,(2)NO2,(3)CH4 oxidation
!
      ROXYD=AMAX1(0.0,ROXYM(N,K)-ROXYO(N,K))
      VMXD4=0.875*ROXYD*XCO2
      VMXDXS=FNO2S*VMXD4*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM)
      VMXDXB=FNO2B*VMXD4*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM)
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD4S=VMXDXS*FVMXDX
      VMXD4B=VMXDXB*FVMXDX
      ZNO2SX=ZNO2S(L,NY,NX)+RVOXA(1)
      ZNO2BX=ZNO2B(L,NY,NX)+RVOXB(1)
      RDNO2(N,K)=AMAX1(0.0,AMIN1(VMXD4S,ZNO2SX))
      RDN2B(N,K)=AMAX1(0.0,AMIN1(VMXD4B,ZNO2BX))
      RDNOT=RDNO2(N,K)+RDN2B(N,K)
      RGOMY(N,K)=0.0
      RGOMD(N,K)=RDNOT*ECNO*ENOX
      RDNO3(N,K)=0.0
      RDNOB(N,K)=0.0
      RDN2O(N,K)=0.0
      RVMX2(N,K,L,NY,NX)=VMXD4S
      RVMB2(N,K,L,NY,NX)=VMXD4B
      RVOXA(N)=RVOXA(N)+0.333*RDNO2(N,K)
      RVOXB(N)=RVOXB(N)+0.333*RDN2B(N,K)
!     TRN2ON(NY,NX)=TRN2ON(NY,NX)+RDNO2(N,K)+RDN2B(N,K)
!     IF((I/1)*1.EQ.I.AND.J.EQ.19.AND.L.LE.5)THEN
!     WRITE(*,7777)'AUTO',I,J,L,K,N,RDNO2(N,K)
!    2,RDN2B(N,K),TRN2ON(NY,NX)
!    2,CNO2S(L,NY,NX),CNO2B(L,NY,NX),CNO3S(L,NY,NX),CNO3B(L,NY,NX)
!    3,Z2OS(L,NY,NX),VLNOB(L,NY,NX),ZNO2S(L,NY,NX),ZNO2B(L,NY,NX)
!    3,XCO2,FNO2,FNB2,TFNG(N,K),OMA(N,K),ROXYP(N,K)
!    2,ROXYM(N,K),ROXYO(N,K),WFN(N,K),FOXYX
!    3,THETW(L,NY,NX),COXYS(L,NY,NX),COXYG(L,NY,NX)
!    4,ROXYD,VMXD4,VMXDXS,VMXDXB,VMXD4S,VMXD4B,FNO2S,FNO2B
!    5,ZNFN4S,ZNFN4B
7777  FORMAT(A8,5I4,50E12.4)
!     ENDIF
      end subroutine AutotrophDenitrificCatabolism
!------------------------------------------------------------------------------------------

      subroutine AerobsO2Uptake(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX

      integer :: M,MX
!     begin_execution
      IF(ROXYP(N,K).GT.ZEROS(NY,NX).AND.FOXYX.GT.ZERO)THEN
      IF(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
!
!     MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC
!     POPULATION
!
      RUPMX=ROXYP(N,K)*XNPG
      ROXYFX=ROXYF(L,NY,NX)*XNPG*FOXYX
      OLSGL1=OLSGL(L,NY,NX)*XNPG
      IF(L.NE.0)THEN
      OXYG1=OXYG(L,NY,NX)*FOXYX
      ROXYLX=ROXYL(L,NY,NX)*XNPG*FOXYX
      ELSE
      OXYG1=COXYG(L,NY,NX)*VOLPM(1,L,NY,NX)*FOXYX
      ROXYLX=(ROXYL(L,NY,NX)+FLQRQ(NY,NX)*COXR(NY,NX) &
      +FLQRI(NY,NX)*COXQ(NY,NX))*XNPG*FOXYX
      ENDIF
      OXYS1=OXYS(L,NY,NX)*FOXYX
!
!     O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
!     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
!
      DO 420 M=1,NPH
!
!     ACTUAL REDUCTION OF AQUEOUS BY AEROBES CALCULATED
!     FROM MASS FLOW PLUS DIFFUSION = ACTIVE UPTAKE
!     COUPLED WITH DISSOLUTION OF GASEOUS O2 DURING REDUCTION
!     OF AQUEOUS O2 FROM DISSOLUTION RATE CONSTANT 'DFGS'
!     CALCULATED IN 'WATSUB'
!
!     VOLWM,VOLPM,VOLX=water, air and total volumes
!     ORAD=microbial radius,FILM=water film thickness
!     DIFOX=aqueous O2 diffusion, TORT=tortuosity
!     BIOS=microbial number, OMA=active biomass
!     SOXYL=O2 solubility, OXKX=Km for O2 uptake
!     OXYS,COXYS=aqueous O2 amount, concentration
!     OXYG,COXYG=gaseous O2 amount, concentration
!     RMPOX,ROXSK=O2 uptake
!
      THETW1=AMAX1(0.0,VOLWM(M,L,NY,NX)/VOLY(L,NY,NX))
      RRADO=ORAD*(FILM(M,L,NY,NX)+ORAD)/FILM(M,L,NY,NX)
      DIFOX=TORT(M,L,NY,NX)*OLSGL1*12.57*BIOS*OMA(N,K)*RRADO
      VOLWOX=VOLWM(M,L,NY,NX)*SOXYL(L,NY,NX)
      VOLPOX=VOLPM(M,L,NY,NX)
      VOLWPM=VOLWOX+VOLPOX
      DO 425 MX=1,NPT
      OXYG1=OXYG1+ROXYFX
      OXYS1=OXYS1+ROXYLX
      COXYS1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX) &
      ,AMAX1(0.0,safe_adb(OXYS1,(VOLWM(M,L,NY,NX)*FOXYX))))
      X=DIFOX*COXYS1
      IF(X.GT.ZEROS(NY,NX).AND.OXYS1.GT.ZEROS(NY,NX))THEN
      B=-RUPMX-DIFOX*OXKX-X
      C=X*RUPMX
      RMPOX=(-B-SQRT(B*B-4.0*C))/2.0
      ELSE
      RMPOX=0.0
      ENDIF
      OXYS1=OXYS1-RMPOX
      IF(THETPM(M,L,NY,NX).GT.THETX.AND.VOLPOX.GT.ZEROS(NY,NX))THEN
      ROXDFQ=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),OXYG1)*VOLWOX &
      -OXYS1*VOLPOX)/VOLWPM
      ELSE
      ROXDFQ=0.0
      ENDIF
      OXYG1=OXYG1-ROXDFQ
      OXYS1=OXYS1+ROXDFQ
      RUPOX(N,K)=RUPOX(N,K)+RMPOX
      ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RMPOX
!     IF(I.EQ.151.AND.J.EQ.24.AND.L.LE.5.AND.M.EQ.NPH.AND.MX.EQ.NPT)THEN
!     WRITE(*,5545)'RMPOX',I,J,L,K,N,M,MX,OXYS1,ROXDFQ,ROXYLX,RMPOX
!    2,DFGS(M,L,NY,NX),OXYG1,VOLWOX,VOLPOX,VOLWPM,X,B,C
!    3,RUPMX,DIFOX,OXKX,COXYS1,FOXYX,ROXYL(L,NY,NX)
!    4,ROXSK(M,L,NY,NX),VOLWM(M,L,NY,NX)/VOLY(L,NY,NX)
!    5,OXYS(L,NY,NX)
5545  FORMAT(A8,7I4,30E16.6)
!     ENDIF
!     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.3
!    2.AND.K.GE.3.AND.N.EQ.3)THEN
!     WRITE(*,5544)'OXY',I,J,L,K,N,M,MX,RUPOX(N,K),ROXYP(N,K)
!    2,ROXSK(M,L,NY,NX),RUPMX,RMPOX,DIFOX,OLSGL1,BIOS,OMA(N,K),X
!    2,ROXDFQ,ROXYLX,ROXYFX,FOXYX,COXYS1,OXYS1,OXYG1,OXYS1
!    4/(VOLWM(M,L,NY,NX)*FOXYX),OXYG1/(VOLPM(M,L,NY,NX)*FOXYX)
!    5,THETW1,THETPM(M,L,NY,NX),DFGS(M,L,NY,NX),ROXSK(M,L,NY,NX)
!    6,VOLPM(M,L,NY,NX),VOLWM(M,L,NY,NX),VOLA(L,NY,NX)
!    7,COXYS(L,NY,NX),COXYG(L,NY,NX),ROXYY(L,NY,NX)
5544  FORMAT(A8,7I4,50E12.4)
!     ENDIF
425   CONTINUE
420   CONTINUE
!
!     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
!
!     WFN=ratio of O2-limited to O2-unlimited uptake
!     RVMX4,RVNHB,RVMX2,RVMB2=NH3,NO2 oxidation in non-band, band
!
      WFN(N,K)=AMIN1(1.0,AMAX1(0.0,RUPOX(N,K)/ROXYP(N,K)))
!     IF(K.LE.4)THEN
!     ROQCS(N,K,L,NY,NX)=ROQCS(N,K,L,NY,NX)*WFN(N,K)
!     ROQAS(N,K,L,NY,NX)=ROQAS(N,K,L,NY,NX)*WFN(N,K)
!     ROQCD(N,K)=ROQCD(N,K)*WFN(N,K)
!     ENDIF
      IF(K.EQ.5)THEN
      IF(N.EQ.1)THEN
      RVMX4(N,K,L,NY,NX)=RVMX4(N,K,L,NY,NX)*WFN(N,K)
      RVMB4(N,K,L,NY,NX)=RVMB4(N,K,L,NY,NX)*WFN(N,K)
      ELSEIF(N.EQ.2)THEN
      RVMX2(N,K,L,NY,NX)=RVMX2(N,K,L,NY,NX)*WFN(N,K)
      RVMB2(N,K,L,NY,NX)=RVMB2(N,K,L,NY,NX)*WFN(N,K)
      ENDIF
      ENDIF
      ELSE
      RUPOX(N,K)=ROXYP(N,K)
      WFN(N,K)=1.0
      ENDIF
      ELSE
      RUPOX(N,K)=0.0
      WFN(N,K)=1.0
      ENDIF
!
!     RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2
!
!     RGOMO,RGOMP=O2-limited, O2-unlimited respiration
!     RCO2X,RCH3X,RCH4X,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
!     ROXYO=O2-limited O2 uptake
!     RVOXA,RVOXB=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
!
      RGOMO(N,K)=RGOMP*WFN(N,K)
      RCO2X(N,K)=RGOMO(N,K)
      RCH3X(N,K)=0.0
      RCH4X(N,K)=0.0
      ROXYO(N,K)=ROXYM(N,K)*WFN(N,K)
      RH2GX(N,K)=0.0
      IF(K.EQ.5)THEN
      RVOXA(N)=RVOXPA*WFN(N,K)
      RVOXB(N)=RVOXPB*WFN(N,K)
      ENDIF
      end subroutine AerobsO2Uptake
!------------------------------------------------------------------------------------------

      subroutine BiomassMineralization(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
!     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMN=microbial nonstructural C,N
!     CNOMC=maximum microbial N:C ratio
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     RINHX=microbially limited NH4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHO,RINHB=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
!     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
!

      RINHP=(OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K)-OMN(3,N,K,L,NY,NX))
      IF(RINHP.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(L,NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(L,NY,NX)-Z4MN)
      RINHX=AMIN1(RINHP,BIOA*OMA(N,K)*TFNG(N,K)*Z4MX)
      RINHO(N,K,L,NY,NX)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
      RINHB(N,K,L,NY,NX)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
      ZNH4M=Z4MN*VOLW(L,NY,NX)*FNH4S
      ZNHBM=Z4MN*VOLW(L,NY,NX)*FNHBS
      RINH4(N,K)=AMIN1(FNH4X*AMAX1(0.0,(ZNH4S(L,NY,NX)-ZNH4M)) &
      ,RINHO(N,K,L,NY,NX))
      RINB4(N,K)=AMIN1(FNB4X*AMAX1(0.0,(ZNH4B(L,NY,NX)-ZNHBM)) &
      ,RINHB(N,K,L,NY,NX))
      ELSE
      RINHO(N,K,L,NY,NX)=0.0
      RINHB(N,K,L,NY,NX)=0.0
      RINH4(N,K)=RINHP*FNH4S
      RINB4(N,K)=RINHP*FNHBS
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+(RINH4(N,K)+RINB4(N,K))
!    2/AREA(3,L,NY,NX)
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,7776)'RINH4',I,J,NX,NY,L,K,N,RINH4(N,K),RINHP
!    1,BIOA*OMA(N,K)*Z4MX*TFNG(N,K),BIOA,OMA(N,K),Z4MX,TFNG(N,K)
!    2,OMC(3,N,K,L,NY,NX),CNOMC(3,N,K),OMN(3,N,K,L,NY,NX)
!    3,RINHO(N,K,L,NY,NX),CNH4S(L,NY,NX),FNH4X,ZNH4S(L,NY,NX)
!    4,ZNH4B(L,NY,NX),ZNH4T(L),OQN(K,L,NY,NX),TRINH4(NY,NX)
7776  FORMAT(A8,7I4,30E12.4)
!     ENDIF
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
!     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINOP=NO3 immobilization (+ve) demand
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     min NO3 concentration and Km for NO3 uptake
!     RINOX=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RINOO,RINOB=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
!     TRINH4=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
      RINOP=AMAX1(0.0,RINHP-RINH4(N,K)-RINB4(N,K))
      IF(RINOP.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(L,NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(L,NY,NX)-ZOMN)
      RINOX=AMIN1(RINOP,BIOA*OMA(N,K)*TFNG(N,K)*ZOMX)
      RINOO(N,K,L,NY,NX)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
      RINOB(N,K,L,NY,NX)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
      ZNO3M=ZOMN*VOLW(L,NY,NX)*FNO3S
      ZNOBM=ZOMN*VOLW(L,NY,NX)*FNO3B
      RINO3(N,K)=AMIN1(FNO3X*AMAX1(0.0,(ZNO3S(L,NY,NX)-ZNO3M)) &
      ,RINOO(N,K,L,NY,NX))
      RINB3(N,K)=AMIN1(FNB3X*AMAX1(0.0,(ZNO3B(L,NY,NX)-ZNOBM)) &
      ,RINOB(N,K,L,NY,NX))
      ELSE
      RINOO(N,K,L,NY,NX)=0.0
      RINOB(N,K,L,NY,NX)=0.0
      RINO3(N,K)=RINOP*FNO3S
      RINB3(N,K)=RINOP*FNO3B
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+(RINO3(N,K)+RINB3(N,K))
!     IF(RINO3(N,K).LT.0.0.OR.RINB3(N,K).LT.0.0)THEN
!     WRITE(*,4321)'RINO3',I,J,NX,NY,L,K,N
!    2,RINOO(N,K,L,NY,NX),RINO3(N,K)
!    2,RINOP,BIOA,OMA(N,K),TFNG(N,K),ZOMX,WFN(N,K),FNO3X,FNO3B
!    2,VLNO3(L,NY,NX),VLNOB(L,NY,NX),CNO3S(L,NY,NX),CNO3B(L,NY,NX)
!    3,RINOB(N,K,L,NY,NX),RINB3,ZNO3S(L,NY,NX),ZNO3B(L,NY,NX)
!    3,OMC(3,N,K,L,NY,NX),CPOMC(3,N,K),OMP(3,N,K,L,NY,NX),WFNG
4321  FORMAT(A8,7I4,60E12.4)
!     ENDIF
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMP=microbial nonstructural C,P
!     CPOMC=maximum microbial P:C ratio
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     min H2PO4 concentration and Km for H2PO4 uptake
!     RIPOX=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RIPOO,RIPBO=substrate-unlimited H2PO4 mineraln-immobiln
!     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
!     VOLW=water content
!     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
!     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
      RIPOP=(OMC(3,N,K,L,NY,NX)*CPOMC(3,N,K)-OMP(3,N,K,L,NY,NX))
      IF(RIPOP.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(L,NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(L,NY,NX)-HPMN)
      RIPOX=AMIN1(RIPOP,BIOA*OMA(N,K)*TFNG(N,K)*HPMX)
      RIPOO(N,K,L,NY,NX)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
      RIPBO(N,K,L,NY,NX)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
      H2POM=HPMN*VOLW(L,NY,NX)*FH2PS
      H2PBM=HPMN*VOLW(L,NY,NX)*FH2PB
      RIPO4(N,K)=AMIN1(FPO4X*AMAX1(0.0,(H2PO4(L,NY,NX)-H2POM)) &
      ,RIPOO(N,K,L,NY,NX))
      RIPOB(N,K)=AMIN1(FPOBX*AMAX1(0.0,(H2POB(L,NY,NX)-H2PBM)) &
      ,RIPBO(N,K,L,NY,NX))
      ELSE
      RIPOO(N,K,L,NY,NX)=0.0
      RIPBO(N,K,L,NY,NX)=0.0
      RIPO4(N,K)=RIPOP*FH2PS
      RIPOB(N,K)=RIPOP*FH2PB
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIPO4(N,K)+RIPOB(N,K))
!     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
!     WRITE(*,4322)'RIPO4',I,J,NX,NY,L,K,N,RIPO4(N,K),FPO4X,H2P4T(L)
!    2,RIPOO(N,K,L,NY,NX),RIPOP,BIOA,OMA(N,K),TFNG(N,K),HPMX,WFN(N,K)
!    2,VLPO4(L,NY,NX),VLPOB(L,NY,NX),CH2P4(L,NY,NX),CH2P4B(L,NY,NX)
!    3,OMC(3,N,K,L,NY,NX),CPOMC(3,N,K),OMP(3,N,K,L,NY,NX),WFNG
4322  FORMAT(A8,7I4,30E12.4)
!     ENDIF
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIP1P=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     min HPO4 concentration and Km for HPO4 uptake
!     RIP1X=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RIPO1,RIPB1=substrate-unlimited HPO4 mineraln-immobiln
!     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
!     VOLW=water content
!     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
!     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
!     TRIPO4=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
      RIP1P=0.1*AMAX1(0.0,RIPOP-RIPO4(N,K)-RIPOB(N,K))
      IF(RIP1P.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(L,NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(L,NY,NX)-HPMN)
      RIP1X=AMIN1(RIP1P,BIOA*OMA(N,K)*TFNG(N,K)*HPMX)
      RIPO1(N,K,L,NY,NX)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
      RIPB1(N,K,L,NY,NX)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
      H1POM=HPMN*VOLW(L,NY,NX)*FH1PS
      H1PBM=HPMN*VOLW(L,NY,NX)*FH1PB
      RIP14(N,K)=AMIN1(FP14X*AMAX1(0.0,(H1PO4(L,NY,NX)-H1POM)) &
      ,RIPO1(N,K,L,NY,NX))
      RIP1B(N,K)=AMIN1(FP1BX*AMAX1(0.0,(H1POB(L,NY,NX)-H1PBM)) &
      ,RIPB1(N,K,L,NY,NX))
      ELSE
      RIPO1(N,K,L,NY,NX)=0.0
      RIPB1(N,K,L,NY,NX)=0.0
      RIP14(N,K)=RIP1P*FH1PS
      RIP1B(N,K)=RIP1P*FH1PB
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIP14(N,K)+RIP1B(N,K))
!     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
!     WRITE(*,4323)'RIP14',I,J,NX,NY,L,K,N,RIP14(N,K),FP14X,H1P4T(L)
!    2,RIPO1(N,K,L,NY,NX),RIP1P,BIOA,OMA(N,K),TFNG(N,K),HPMX,WFN(N,K)
!    2,VLPO4(L,NY,NX),VLPOB(L,NY,NX),CH1P4(L,NY,NX),CH1P4B(L,NY,NX)
!    3,OMC(3,N,K,L,NY,NX),CPOMC(3,N,K),OMP(3,N,K,L,NY,NX),WFNG
4323  FORMAT(A8,7I4,30E12.4)
!     ENDIF
!
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RINHPR=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHOR=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M=NH4 not available for uptake
!     FNH4XR=fractions of biological NH4 demand
!     RINH4R=substrate-limited NH4 mineraln-immobiln
!     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
!
      IF(L.EQ.0)THEN
      RINHPR=RINHP-RINH4(N,K)-RINO3(N,K)
      IF(RINHPR.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(NU(NY,NX),NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(NU(NY,NX),NY,NX)-Z4MN)
      RINHOR(N,K,NY,NX)=AMIN1(RINHPR,BIOA*OMA(N,K)*TFNG(N,K)*Z4MX) &
      *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y &
      /(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VOLW(NU(NY,NX),NY,NX)
      RINH4R(N,K)=AMIN1(FNH4XR(N,K)*AMAX1(0.0 &
      ,(ZNH4T(NU(NY,NX))-ZNH4M)),RINHOR(N,K,NY,NX))
      ELSE
      RINHOR(N,K,NY,NX)=0.0
      RINH4R(N,K)=RINHPR
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+RINH4R(N,K)
!    2/AREA(3,L,NY,NX)
!     IF(K.EQ.2.AND.N.EQ.1)THEN
!     WRITE(*,7778)'RINH4R',I,J,NX,NY,L,K,N,RINH4R(N,K),RINHPR
!    2,BIOA*OMA(N,K)*Z4MX,RINHP,RINH4(N,K),RINO3(N,K)
!    3,RINHOR(N,K,NY,NX),CNH4S(NU(NY,NX),NY,NX),FNH4XR(N,K)
!    4,ZNH4T(NU(NY,NX))
7778  FORMAT(A8,7I4,20E12.4)
!     ENDIF
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RINOPR=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     minimum NO3 concentration and Km for NO3 uptake
!     RINOOR=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RINO3R=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M=NO3 not available for uptake
!     FNO3XR=fraction of biological NO3 demand
!     RINO3R=substrate-limited NO3 immobiln
!     TRINH4=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
      RINOPR=AMAX1(0.0,RINHPR-RINH4R(N,K))
      IF(RINOPR.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(NU(NY,NX),NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(NU(NY,NX),NY,NX)-ZOMN)
      RINOOR(N,K,NY,NX)=AMAX1(RINOPR,BIOA*OMA(N,K)*TFNG(N,K)*ZOMX) &
      *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y &
      /(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VOLW(NU(NY,NX),NY,NX)
      RINO3R(N,K)=AMIN1(FNO3XR(N,K)*AMAX1(0.0 &
      ,(ZNO3T(NU(NY,NX))-ZNO3M)),RINOOR(N,K,NY,NX))
      ELSE
      RINOOR(N,K,NY,NX)=0.0
      RINO3R(N,K)=RINOPR
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+RINO3R(N,K)
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RIPOPR=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     minimum H2PO4 concentration and Km for H2PO4 uptake
!     RIPOOR=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RIPOOR=substrate-unlimited H2PO4 mineraln-immobiln
!     VOLW=water content
!     H2P4M=H2PO4 not available for uptake
!     FPO4XR=fractions of biological H2PO4 demand
!     RIPO4R=substrate-limited H2PO4 mineraln-immobiln
!     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
      RIPOPR=RIPOP-RIPO4(N,K)
      IF(RIPOPR.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(NU(NY,NX),NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPOOR(N,K,NY,NX)=AMIN1(RIPOPR,BIOA*OMA(N,K)*TFNG(N,K)*HPMX) &
      *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY &
      /(CH2PY+HPKU))
      H2P4M=HPMN*VOLW(NU(NY,NX),NY,NX)
      RIPO4R(N,K)=AMIN1(FPO4XR(N,K)*AMAX1(0.0 &
      ,(H2P4T(NU(NY,NX))-H2P4M)),RIPOOR(N,K,NY,NX))
      ELSE
      RIPOOR(N,K,NY,NX)=0.0
      RIPO4R(N,K)=RIPOPR
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIPO4R(N,K)
!     WRITE(*,7778)'RIPO4R',I,J,NX,NY,L,K,N,RIPO4R(N,K),FPO4XR(N,K)
!    2,H2P4T(NU(NY,NX)),RIPOOR(N,K,NY,NX),RIPOPR
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RIP1PR=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     minimum HPO4 concentration and Km for HPO4 uptake
!     RIPO1R=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RIPO1R=substrate-unlimited HPO4 mineraln-immobiln
!     VOLW=water content
!     H1P4M=HPO4 not available for uptake
!     FP14XR=fraction of biological HPO4 demand
!     RIP14R=substrate-limited HPO4 minereraln-immobiln
!     TRIPO4=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
      RIP1PR=0.1*AMAX1(0.0,RIPOPR-RIPO4R(N,K))
      IF(RIP1PR.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(NU(NY,NX),NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPO1R(N,K,NY,NX)=AMIN1(RIP1PR,BIOA*OMA(N,K)*TFNG(N,K)*HPMX) &
      *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY &
      /(CH1PY+HPKU))
      H1P4M=HPMN*VOLW(NU(NY,NX),NY,NX)
      RIP14R(N,K)=AMIN1(FP14XR(N,K)*AMAX1(0.0 &
      ,(H1P4T(NU(NY,NX))-H1P4M)),RIPO1R(N,K,NY,NX))
      ELSE
      RIPO1R(N,K,NY,NX)=0.0
      RIP14R(N,K)=RIP1PR
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIP14R(N,K)
!     WRITE(*,7778)'RIP14R',I,J,NX,NY,L,K,N,RIP14R(N,K),FP14XR(N,K)
!    2,H1P4T(NU(NY,NX)),RIPO1R(N,K,NY,NX),RIP1PR
      ENDIF
      end subroutine BiomassMineralization
!------------------------------------------------------------------------------------------

      subroutine GatherMicrobialRespiration(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX
!     begin_execution
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TFNR=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
      FPH=1.0+AMAX1(0.0,0.25*(6.5-PH(L,NY,NX)))
      RMOMX=RMOM*TFNR(N,K)*FPH
      RMOMC(1,N,K)=OMN(1,N,K,L,NY,NX)*RMOMX*RMOMK(1)
      RMOMC(2,N,K)=OMN2(N,K)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMOMT=total maintenance respiration
!     RGOMT=growth respiration
!     RXOMT=senescence respiration
!
      RMOMT=RMOMC(1,N,K)+RMOMC(2,N,K)
      RGOMT=AMAX1(0.0,RGOMO(N,K)-RMOMT)
      RXOMT=AMAX1(0.0,RMOMT-RGOMO(N,K))
!
!     N2 FIXATION: N=(6) AEROBIC, (7) ANAEROBIC
!     FROM GROWTH RESPIRATION, FIXATION ENERGY REQUIREMENT,
!     MICROBIAL N REQUIREMENT IN LABILE (1) AND
!     RESISTANT (2) FRACTIONS
!
!     RGN2P=respiration to meet N2 fixation demand
!     OMC,OMN=microbial nonstructural C,N
!     CNOMC=maximum microbial N:C ratio
!     EN2F=N2 fixation yield per unit nonstructural C
!     RGOMT=growth respiration
!     RGN2F=respiration for N2 fixation
!     CZ2GS=aqueous N2 concentration
!     ZFKM=Km for N2 uptake
!     OMGR*OMC(3,N,K,L,NY,NX)=nonstructural C limitation to RGN2F
!     RN2FX=N2 fixation rate
!
      IF(K.LE.4.AND.(N.EQ.6.OR.N.EQ.7))THEN
      RGN2P=AMAX1(0.0,OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K) &
      -OMN(3,N,K,L,NY,NX))/EN2F(N)
      IF(RGOMT.GT.ZEROS(NY,NX))THEN
      RGN2F(N,K)=AMIN1(RGOMT*RGN2P/(RGOMT+RGN2P) &
      *CZ2GS(L,NY,NX)/(CZ2GS(L,NY,NX)+ZFKM) &
      ,OMGR*OMC(3,N,K,L,NY,NX))
      ELSE
      RGN2F(N,K)=0.0
      ENDIF
      RN2FX(N,K)=RGN2F(N,K)*EN2F(N)
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,5566)'N2 FIX',I,J,NX,NY,L,K,N,RN2FX(N,K),EN2F(N)
!    2,OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K),OMN(3,N,K,L,NY,NX)
!    3,RINH4(N,K),RINO3(N,K),RGN2P,RGN2F(N,K),RGOMT
!    4,CZ2GS(L,NY,NX)
5566  FORMAT(A8,7I4,30E12.4)
!     ENDIF
      ELSE
      RN2FX(N,K)=0.0
      RGN2F(N,K)=0.0
      ENDIF
      end subroutine GatherMicrobialRespiration
!------------------------------------------------------------------------------------------

      subroutine GetMicrobialAnabolismFlux(N,K,L,NY,NX)
      implicit none
      integer, intent(in) :: N,K,L,NY,NX

      integer :: M
      real(r8) :: RCCC,RCCN,RCCP
      real(r8) :: CCC,CGOMX,CGOMD
      real(r8) :: CXC
!     begin_execution
!     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
!     FROM O2, NOX AND C REDUCTION
!
!     CGOMX=DOC+acetate uptake from aerobic growth respiration
!     CGOMD=DOC+acetate uptake from denitrifier growth respiration
!     RMOMT=maintenance respiration
!     RGOMO=total respiration
!     RGOMD=respiration for denitrifcation
!     RGN2F=respiration for N2 fixation
!     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
!     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake(heterotrophs
!     CGOMC=total CO2,CH4 uptake (autotrophs)
!     CGOMN,CGOMP=DON, DOP uptake
!     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
!     OQN,OPQ=DON,DOP
!     FOMK=faction of OMA in total OMA
!     CNQ,CPQ=DON/DOC, DOP/DOC
!     FCN,FCP=limitation from N,P
!
      CGOMX=AMIN1(RMOMT,RGOMO(N,K))+RGN2F(N,K) &
      +(RGOMT-RGN2F(N,K))/ECHZ
      CGOMD=RGOMD(N,K)/ENOX
      CGOMC(N,K)=CGOMX+CGOMD
      IF(K.LE.4)THEN
      CGOQC(N,K)=CGOMX*FGOCP+CGOMD
      CGOAC(N,K)=CGOMX*FGOAP
      CGOXC=CGOQC(N,K)+CGOAC(N,K)
      CGOMN(N,K)=AMAX1(0.0,AMIN1(OQN(K,L,NY,NX)*FOMK(N,K) &
      ,CGOXC*CNQ(K)/FCN(N,K)))
      CGOMP(N,K)=AMAX1(0.0,AMIN1(OQP(K,L,NY,NX)*FOMK(N,K) &
      ,CGOXC*CPQ(K)/FCP(N,K)))
      ELSE
      CGOQC(N,K)=CGOMX+CGOMD
      CGOAC(N,K)=0.0
      CGOMN(N,K)=0.0
      CGOMP(N,K)=0.0
      ENDIF
      TCGOQC(K)=TCGOQC(K)+CGOQC(N,K)
      TCGOAC(K)=TCGOAC(K)+CGOAC(N,K)
      TCGOMN(K)=TCGOMN(K)+CGOMN(N,K)
      TCGOMP(K)=TCGOMP(K)+CGOMP(N,K)
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.3)THEN
!     WRITE(*,5557)'CGOQC',I,J,NX,NY,L,K,N,CGOQC(N,K),CGOMX
!    2,FGOCP,FGOAP,CGOMD,RMOMT,RGN2F(N,K),ECHZ
!    3,RGOMD(N,K),ENOX,RGOMO(N,K),WFN(N,K),FOXYX
!     WRITE(*,5557)'CGOMP',I,J,NX,NY,L,K,N,CGOMP(N,K),OQP(K,L,NY,NX)
!    2,FOMK(N,K),CGOXC,CPQ(K),FCP(N,K),CGOQC(N,K),CGOAC(N,K)
5557  FORMAT(A8,7I4,30E12.4)
!     ENDIF
!
!     TRANSFER UPTAKEN C,N,P FROM STORAGE TO ACTIVE BIOMASS
!
!     OMC,OMN,OMP=nonstructural C,N,P
!     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
!     CNOMC,CPOMC=maximum microbial N:C, P:C ratios
!     RCCC,RCCN,RCCP=C,N,P recycling fractions
!     RCCZ,RCCY=min, max C recycling fractions
!     RCCX,RCCQ=max N,P recycling fractions
!
      IF(OMC(3,N,K,L,NY,NX).GT.ZEROS(NY,NX) &
      .AND.OMC(1,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CCC=AMAX1(0.0,AMIN1(1.0 &
      ,OMN(3,N,K,L,NY,NX)/(OMN(3,N,K,L,NY,NX) &
      +OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K)) &
      ,OMP(3,N,K,L,NY,NX)/(OMP(3,N,K,L,NY,NX) &
      +OMC(3,N,K,L,NY,NX)*CPOMC(3,N,K))))
      CXC=OMC(3,N,K,L,NY,NX)/OMC(1,N,K,L,NY,NX)
      C3C=1.0/(1.0+CXC/CKC)
      CNC=AMAX1(0.0,AMIN1(1.0 &
      ,OMC(3,N,K,L,NY,NX)/(OMC(3,N,K,L,NY,NX) &
      +OMN(3,N,K,L,NY,NX)/CNOMC(3,N,K))))
      CPC=AMAX1(0.0,AMIN1(1.0 &
      ,OMC(3,N,K,L,NY,NX)/(OMC(3,N,K,L,NY,NX) &
      +OMP(3,N,K,L,NY,NX)/CPOMC(3,N,K))))
      RCCC=RCCZ+AMAX1(CCC,C3C)*RCCY
      RCCN=CNC*RCCX
      RCCP=CPC*RCCQ
      ELSE
      RCCC=RCCZ
      RCCN=0.0
      RCCP=0.0
      ENDIF
!     IF((I/120)*120.EQ.I.AND.J.EQ.24)THEN
!     WRITE(*,5555)'RCCC',I,J,NX,NY,L,K,N,RCCC,RCCN,RCCP
!    2,OMC(3,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX),OMP(3,N,K,L,NY,NX)
!    3,CCC,C3C,CNC,CPC
!     ENDIF
!
!     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
!
!     CGOMZ=transfer from nonstructural to structural microbial C
!     TFNG=temperature+water stress function
!     OMGR=rate constant for transferring nonstructural to structural C
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     FL=partitioning between labile and resistant microbial components
!     OMC,OMN,OMP=nonstructural microbial C,N,P
!
      CGOMZ=TFNG(N,K)*OMGR*AMAX1(0.0,OMC(3,N,K,L,NY,NX))
      DO 745 M=1,2
      CGOMS(M,N,K)=FL(M)*CGOMZ
      IF(OMC(3,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CGONS(M,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMN(3,N,K,L,NY,NX)) &
      ,CGOMS(M,N,K)*OMN(3,N,K,L,NY,NX)/OMC(3,N,K,L,NY,NX))
      CGOPS(M,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMP(3,N,K,L,NY,NX)) &
      ,CGOMS(M,N,K)*OMP(3,N,K,L,NY,NX)/OMC(3,N,K,L,NY,NX))
      ELSE
      CGONS(M,N,K)=0.0
      CGOPS(M,N,K)=0.0
      ENDIF
!
!     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
!     RATE, TEMPERATURE
!
!     SPOMX=rate constant for microbial decomposition
!     SPOMC=basal decomposition rate
!     SPOMK=effect of low microbial C concentration on microbial decay
!     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
!     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall
!     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
!
      SPOMX=SQRT(TFNG(N,K))*SPOMC(M)*SPOMK(M)
      RXOMC(M,N,K)=AMAX1(0.0,OMC(M,N,K,L,NY,NX)*SPOMX)
      RXOMN(M,N,K)=AMAX1(0.0,OMN(M,N,K,L,NY,NX)*SPOMX)
      RXOMP(M,N,K)=AMAX1(0.0,OMP(M,N,K,L,NY,NX)*SPOMX)
      RDOMC(M,N,K)=RXOMC(M,N,K)*(1.0-RCCC)
      RDOMN(M,N,K)=RXOMN(M,N,K)*(1.0-RCCC)*(1.0-RCCN)
      RDOMP(M,N,K)=RXOMP(M,N,K)*(1.0-RCCC)*(1.0-RCCP)
      R3OMC(M,N,K)=RXOMC(M,N,K)-RDOMC(M,N,K)
      R3OMN(M,N,K)=RXOMN(M,N,K)-RDOMN(M,N,K)
      R3OMP(M,N,K)=RXOMP(M,N,K)-RDOMP(M,N,K)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!
      RHOMC(M,N,K)=AMAX1(0.0,RDOMC(M,N,K)*EHUM(L,NY,NX))
      RHOMN(M,N,K)=AMAX1(0.0,RDOMN(M,N,K)*EHUM(L,NY,NX))
      RHOMP(M,N,K)=AMAX1(0.0,RDOMP(M,N,K)*EHUM(L,NY,NX))
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,8821)'RHOMC',I,J,L,K,N,M
!    2,RXOMC(M,N,K),RXOMN(M,N,K),RXOMP(M,N,K)
!    2,RDOMC(M,N,K),RDOMN(M,N,K),RDOMP(M,N,K)
!    2,R3OMC(M,N,K),R3OMN(M,N,K),R3OMP(M,N,K)
!    2,RHOMC(M,N,K),RHOMN(M,N,K),RHOMP(M,N,K)
!    4,OMC(M,N,K,L,NY,NX),OMN(M,N,K,L,NY,NX)
!    5,OMP(M,N,K,L,NY,NX)
!    4,OMC(3,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX)
!    5,OMP(3,N,K,L,NY,NX)
!    6,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
!    2,SPOMX,RCCC,RCCN,RCCP
!     ENDIF
!
!     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
!
      RCOMC(M,N,K)=RDOMC(M,N,K)-RHOMC(M,N,K)
      RCOMN(M,N,K)=RDOMN(M,N,K)-RHOMN(M,N,K)
      RCOMP(M,N,K)=RDOMP(M,N,K)-RHOMP(M,N,K)
745   CONTINUE
!
!     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
!     EXCEEDS UPTAKE
!
!     OMC,OMN,OMP=microbial C,N,P
!     RMOMT=total maintenance respiration
!     RXOMT=senescence respiration
!     RCCC=C recycling fraction
!     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
!     RMOMC=maintenance respiration
!     CNOMA,CPOMA=N:C,P:C ratios of active biomass
!     RDMMC,RDMMN,RDMMP=microbial C,N,P litterfall from senescence
!     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!
      IF(RXOMT.GT.ZEROS(NY,NX).AND.RMOMT.GT.ZEROS(NY,NX) &
      .AND.RCCC.GT.ZERO)THEN
      FRM=RXOMT/RMOMT
      DO 730 M=1,2
      RXMMC(M,N,K)=AMIN1(OMC(M,N,K,L,NY,NX) &
      ,AMAX1(0.0,FRM*RMOMC(M,N,K)/RCCC))
      RXMMN(M,N,K)=AMIN1(OMN(M,N,K,L,NY,NX) &
      ,AMAX1(0.0,RXMMC(M,N,K)*CNOMA(N,K)))
      RXMMP(M,N,K)=AMIN1(OMP(M,N,K,L,NY,NX) &
      ,AMAX1(0.0,RXMMC(M,N,K)*CPOMA(N,K)))
      RDMMC(M,N,K)=RXMMC(M,N,K)*(1.0-RCCC)
      RDMMN(M,N,K)=RXMMN(M,N,K)*(1.0-RCCN)*(1.0-RCCC)
      RDMMP(M,N,K)=RXMMP(M,N,K)*(1.0-RCCP)*(1.0-RCCC)
      R3MMC(M,N,K)=RXMMC(M,N,K)-RDMMC(M,N,K)
      R3MMN(M,N,K)=RXMMN(M,N,K)-RDMMN(M,N,K)
      R3MMP(M,N,K)=RXMMP(M,N,K)-RDMMP(M,N,K)
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!     EHUM=humus transfer fraction
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
      RHMMC(M,N,K)=AMAX1(0.0,RDMMC(M,N,K)*EHUM(L,NY,NX))
      RHMMN(M,N,K)=AMAX1(0.0,RDMMN(M,N,K)*EHUM(L,NY,NX))
      RHMMP(M,N,K)=AMAX1(0.0,RDMMP(M,N,K)*EHUM(L,NY,NX))
      RCMMC(M,N,K)=RDMMC(M,N,K)-RHMMC(M,N,K)
      RCMMN(M,N,K)=RDMMN(M,N,K)-RHMMN(M,N,K)
      RCMMP(M,N,K)=RDMMP(M,N,K)-RHMMP(M,N,K)
!     IF(L.EQ.11.AND.K.EQ.1)THEN
!     WRITE(*,8821)'RCMMC',I,J,L,K,N,M,RCMMC(M,N,K)
!    2,RDMMC(M,N,K),RHMMC(M,N,K),OMC(M,N,K,L,NY,NX)
!    3,FRM,RMOMC(M,N,K),OMN(1,N,K,L,NY,NX),OMN2(N,K)
!    4,RMOM,TFNR(N,K),FPH,RDMMN(M,N,K),CNSHZ,RDMMP(M,N,K)
!    5,CPSHZ,EHUM(L,NY,NX),RXOMT,RMOMT,RMOMT,RGOMO(N,K)
!    6,RGOMP,WFN(N,K)
!     WRITE(*,8821)'RCMMP',I,J,L,K,N,M,RCMMP(M,N,K)
!    2,RDMMP(M,N,K),RHMMP(M,N,K),EHUM(L,NY,NX)
!    3,RCCC,RCCN,RCCP,RXMMP(M,N,K)
!     ENDIF
730   CONTINUE
      ELSE
      DO 720 M=1,2
      RXMMC(M,N,K)=0.0
      RXMMN(M,N,K)=0.0
      RXMMP(M,N,K)=0.0
      RDMMC(M,N,K)=0.0
      RDMMN(M,N,K)=0.0
      RDMMP(M,N,K)=0.0
      R3MMC(M,N,K)=0.0
      R3MMN(M,N,K)=0.0
      R3MMP(M,N,K)=0.0
      RHMMC(M,N,K)=0.0
      RHMMN(M,N,K)=0.0
      RHMMP(M,N,K)=0.0
      RCMMC(M,N,K)=0.0
      RCMMN(M,N,K)=0.0
      RCMMP(M,N,K)=0.0
720   CONTINUE
      ENDIF

      end subroutine GetMicrobialAnabolismFlux
!------------------------------------------------------------------------------------------

      subroutine DealNoActiveMicrobes(N,K,L)
!
!     zero out within modulue flux arrays
      implicit none
      integer, intent(in) :: N,K,L

      integer :: M
!     begin_execution

      RUPOX(N,K)=0.0
      RGOMO(N,K)=0.0
      RCO2X(N,K)=0.0
      RCH3X(N,K)=0.0
      RCH4X(N,K)=0.0
      RGOMY(N,K)=0.0
      RGOMD(N,K)=0.0
      CGOMC(N,K)=0.0
      CGOMN(N,K)=0.0
      CGOMP(N,K)=0.0
      CGOQC(N,K)=0.0
      CGOAC(N,K)=0.0
      RDNO3(N,K)=0.0
      RDNOB(N,K)=0.0
      RDNO2(N,K)=0.0
      RDN2B(N,K)=0.0
      RDN2O(N,K)=0.0
      RN2FX(N,K)=0.0
      RINH4(N,K)=0.0
      RINO3(N,K)=0.0
      RIPO4(N,K)=0.0
      RIP14(N,K)=0.0
      RINB4(N,K)=0.0
      RINB3(N,K)=0.0
      RIPOB(N,K)=0.0
      RIP1B(N,K)=0.0
      IF(L.EQ.0)THEN
      RINH4R(N,K)=0.0
      RINO3R(N,K)=0.0
      RIPO4R(N,K)=0.0
      RIP14R(N,K)=0.0
      FNH4XR(N,K)=0.0
      FNO3XR(N,K)=0.0
      FPO4XR(N,K)=0.0
      FP14XR(N,K)=0.0
      ENDIF
      DO 725 M=1,2
      CGOMS(M,N,K)=0.0
      CGONS(M,N,K)=0.0
      CGOPS(M,N,K)=0.0
      RMOMC(M,N,K)=0.0
      RXMMC(M,N,K)=0.0
      RXMMN(M,N,K)=0.0
      RXMMP(M,N,K)=0.0
      RDMMC(M,N,K)=0.0
      RDMMN(M,N,K)=0.0
      RDMMP(M,N,K)=0.0
      R3MMC(M,N,K)=0.0
      R3MMN(M,N,K)=0.0
      R3MMP(M,N,K)=0.0
      RHMMC(M,N,K)=0.0
      RHMMN(M,N,K)=0.0
      RHMMP(M,N,K)=0.0
      RCMMC(M,N,K)=0.0
      RCMMN(M,N,K)=0.0
      RCMMP(M,N,K)=0.0
      RXOMC(M,N,K)=0.0
      RXOMN(M,N,K)=0.0
      RXOMP(M,N,K)=0.0
      RDOMC(M,N,K)=0.0
      RDOMN(M,N,K)=0.0
      RDOMP(M,N,K)=0.0
      R3OMC(M,N,K)=0.0
      R3OMN(M,N,K)=0.0
      R3OMP(M,N,K)=0.0
      RHOMC(M,N,K)=0.0
      RHOMN(M,N,K)=0.0
      RHOMP(M,N,K)=0.0
      RCOMC(M,N,K)=0.0
      RCOMN(M,N,K)=0.0
      RCOMP(M,N,K)=0.0
725   CONTINUE
      RH2GX(N,K)=0.0
      IF(K.EQ.5)THEN
      RVOXA(N)=0.0
      RVOXB(N)=0.0
      IF(N.EQ.5)THEN
      RH2GZ=0.0
      ENDIF
      ENDIF
      end subroutine DealNoActiveMicrobes

      end module NitroMod
