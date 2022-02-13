module grosubMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : test_aeqb,safe_adb
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use SOMDataType
  implicit none

  private
  include "parameters.h"
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
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk14.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"
! end_include_section

  character(len=*), private, parameter :: mod_filename = __FILE__
! DIMENSION VCO2(400,366,05)
!
! PART1X,PART2X=modifiers to organ partitioning coefficients
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! FSNR=rate constant for litterfall at end of growing season (h-1)
! FLG4X=number of hours with no grain filling required for physl maturity
! FLGZX=number of hours until full senescence after physl maturity
! XFRX=maximum storage C content for remobiln from stalk,root reserves
! XFRY=rate const for remobiln to storage from stalk,root reserves (h-1)
! IFLGQ,IFLGQX=current,required hours after physl maturity until start of litterfall
! FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
! FXFS=rate constant for remobilization of stalk C,N,P (h-1)
! FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
!
  real(r8), PARAMETER :: PART1X=0.05_r8,PART2X=0.02_r8 &
      ,VMXC=0.015_r8,FSNR=2.884E-03_r8,FLG4X=168.0_r8 &
      ,FLGZX=240.0_r8,XFRX=2.5E-02_r8,XFRY=2.5E-03_r8 &
      ,FSNK=0.05_r8,FXFS=1.0_r8,FMYC=0.1_r8
  integer, parameter :: IFLGQX=960
!
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth (g N,P g-1 C)
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     PSILM=minimum water potential for organ expansion,extension (MPa)
!     RCMN=minimum stomatal resistance to CO2 (s m-1)
!     RTDPX=distance behind growing point for secondary roots (m)
!     RTLGAX=minimum average secondary root length (m)
!     EMODR=root modulus of elasticity (MPa)
!
  real(r8), PARAMETER :: CNKI=1.0E-01_r8,CPKI=1.0E-02_r8
  real(r8), PARAMETER :: RMPLT=0.010_r8,PSILM=0.1_r8,RCMN=1.560E+01_r8 &
      ,RTDPX=0.00_r8,RTLGAX=1.0E-03_r8,EMODR=5.0_r8
!
!     QNTM=quantum efficiency (umol e- umol-1 PAR)
!     CURV=shape parameter for e- transport response to PAR
!     ELEC3,ELEC4=e- requirement for CO2 fixn by rubisco,PEP carboxylase
!     (umol e- umol CO2)
!     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll in C4 (uM)
!     FCO2B,FHCOB=partition decarboxylation,leakage to CO2,HCO3 in C4
!     COMP4=C4 CO2 compensation point (uM)
!     FDML=leaf water content (g H2O g-1 C)
!     FBS,FMP=leaf water content in bundle sheath, mesophyll in C4 CO2 fixn
!
      real(r8), PARAMETER :: QNTM=0.45_r8,CURV=0.70_r8,CURV2=2.0_r8*CURV &
      ,CURV4=4.0_r8*CURV,ELEC3=4.5_r8,ELEC4=3.0_r8,CO2KI=1.0E+03_r8,FCO2B=0.02_r8 &
      ,FHCOB=1.0_r8-FCO2B
      real(r8), PARAMETER :: COMP4=0.5_r8,FDML=6.0_r8,FBS=0.2_r8*FDML &
      ,FMP=0.8_r8*FDML
!
!     ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     GY,GZ=partitioning of grazed material to removal,respiration
!     FSTK=fraction of stalk area contributing to water,heat flow
!     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
!     FRTX=fraction used to calculate woody faction of stalk,root
!
      real(r8),PARAMETER :: ZPLFM=0.33_r8,ZPLFD=1.0_r8-ZPLFM,ZPGRM=0.75_r8 &
      ,ZPGRD=1.0_r8-ZPGRM,GY=1.0_r8,GZ=1.0_r8-GY
      real(r8), PARAMETER :: FSTK=0.05_r8,ZSTX=1.0E-03_r8,DSTK=0.225_r8 &
      ,VSTK=1.0E-06_r8/DSTK,FRTX=1.0_r8/(1.0_r8-(1.0_r8-FSTK)**2)
!
!     SETC,SETN,SETP=Km for nonstructural C,N,P concn on seed set (g g-1)
!     SLA2,SSL2,SNL2=parameter for calculating leaf area expansion, petiole
!     and internode extension vs leaf, petiole, internode growth
!     CNMX,CPMX,CNMN,CPMN=max,min N:C,P:C for nonstructural C,N,P transfers
!
      real(r8), PARAMETER :: SETC=1.0E-02_r8,SETN=1.0E-03_r8,SETP=1.0E-04_r8
      real(r8), PARAMETER :: SLA2=-0.33_r8,SSL2=-0.50_r8,SNL2=-0.67_r8
      real(r8), PARAMETER :: CNMX=0.20_r8,CPMX=0.020_r8,CNMN=0.050_r8,CPMN=0.005_r8
!
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     VMXO=specific respiration rate by bacterial N2 fixers (g g-1 h-1)
!     SPNDL=specific decomposition rate by canopy,root bacterial N2 fixers (g g-1 h-1)
!     CCNGB,CCNGR=parameters to calculate nonstructural C,N,P exchange between bacteria and branch,root
!     WTNDI=initial bacterial mass at infection (g C m-2)
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria (g N,P g-1 C)
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     RCCZ,RCCY=min,max fractions for shoot,bacteria C recycling
!     RCCX,RCCQ=max fractions for shoot,bacteria N,P recycling
!
      real(r8), PARAMETER :: EN2F=0.20_r8,VMXO=0.125_r8,SPNDLK=0.01_r8 &
      ,SPNDL=5.0E-04_r8,CCNGR=2.5E-01_r8,CCNGB=6.0E-04_r8 &
      ,WTNDI=1.0E-03_r8,CZKM=2.5E-03_r8,CPKM=2.5E-04_r8 &
      ,RCCZR=0.056_r8,RCCYR=0.167_r8,RCCXR=0.833_r8,RCCQR=0.833_r8 &
      ,RCCZN=0.167_r8,RCCYN=0.833_r8,RCCXN=0.833_r8,RCCQN=0.833_r8

      integer :: ICHK1(2,JZ),NRX(2,JZ),NBZ(10)
      real(r8) :: PART(7),TFN6(JZ),ARSTKB(10) &
      ,FXFB(0:3),FXFR(0:3),FXRT(0:1),FXSH(0:1),FXRN(6) &
      ,WTLSBZ(10),CPOOLZ(10),ZPOOLZ(10),PPOOLZ(10) &
      ,ZCX(JP,JY,JX),UPNFC(JP,JY,JX),FRSV(0:3),FXFY(0:1),FXFZ(0:1)

      real(r8) :: RTNT(2),RLNT(2,JZ),RTSK1(2,JZ,10),RTSK2(2,JZ,10) &
      ,RTDPL(10,JZ),FWTR(JZ),FWTB(JP),FRTDP(0:3),RCCX(0:3),RCCQ(0:3) &
      ,RCCZ(0:3),RCCY(0:3),EFIRE(2,5:5),WGLFBL(JZ,10,JP,JY,JX) &
      ,WTSHTA(JZ,JY,JX),FLG4Y(0:3),ATRPX(0:1),GVMX(0:1),RTSK(0:3)

      real(r8) :: CH2O3(25),CH2O4(25),CPOOLK(10,JP,JY,JX),FHVSTK(0:25) &
      ,FHVSHK(0:25),WFNGR(2,JZ),PSILY(0:3),FVRN(0:5)
      real(r8) :: FWOOD(0:1),FWOODN(0:1),FWOODP(0:1) &
      ,FWODB(0:1),FWODLN(0:1),FWODLP(0:1),FWODSN(0:1),FWODSP(0:1) &
      ,FWODR(0:1),FWODRN(0:1),FWODRP(0:1)

      real(r8) :: ACTVM,ARLFI,ALLOCL,ALLOCS,ALLOCN,ASTV,ATRPPD,ARLFY
      real(r8) :: ARLFR,ARLFG,APSILT,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW
      real(r8) :: CPRTW,CNLFB,CNSHB,CPLFB,CPSHB,CNLFM,CPLFM,CNLFX
      real(r8) :: CPLFX,CNSHX,CPSHX,CO2F,CH2O,CO2X,CO2C,CO2Y,CBXNX
      real(r8) :: CNPG,CGROS,CNRDA,CGROSM,CNRDM,CPL4M,CCBS,CPL3K
      real(r8) :: CO2LK,CCE,CCC,CNC,CPC,CPOOLT,CPOOLM,CH2OH,CWTRSV
      real(r8) :: CWTRSN,CWTRSP,CNR,CPR,CNL,CPL,CPOOLD,CCPOLN,CZPOLN
      real(r8) :: CPPOLN,CGNDL,CCNDLB,CUPRL,CUPRO,CUPRC,CGRORM,CGROR
      real(r8) :: CPOOLX,CPOOLNX,CCNDLR,CPOOLB,CPOOLS,CCPOLX,CCPLNX
      real(r8) :: CPOOLG,CPOLNG,CPOLNX,DMLFB,DMSHB,DMSHT,DMSHD,DIFF
      real(r8) :: DATRP,DMRTD,DMRTR,ETLF4,EGRO4,EGROX,ETLF,EGRO
      real(r8) :: ETOL,FPART1,FPART2,FPARTL,FNP,FSNC,FSNCL,FSNCS
      real(r8) :: FNCLF,FSNAL,FSNAS,FRCC,FSNCK,FSNCR,FRACL,FGRNX
      real(r8) :: FXFC,FXFN,FWTBR,FCNPF,FXRNX,FRTN,FSNC1,FRCO2,FSNCM
      real(r8) :: FSNCP,FGROL,FGROZ,FWTRT,FWTC,FWTS,FHGT,FHVST,FHVSH
      real(r8) :: FHVST4,FHGTK,FHVSTS,FHVSTG,FHVSHG,FHVSTH,FHVSTE
      real(r8) :: FHVSHH,FHVSHE,FDM,FFIRE,FFIRN,FFIRP,FHVSN,FHVSP
      real(r8) :: GROGR,GSL,GROLF,GROSHE,GROSTK,GRORSV,GROHSK,GROEAR
      real(r8) :: GROSHT,GROLFN,GROSHN,GROSTN,GRORSN,GROHSN,GROEAN
      real(r8) :: GROGRN,GROLFP,GROSHP,GROSTP,GRORSP,GROHSP,GROEAP
      real(r8) :: GROGRP,GNOD,GRO,GRON,GROP,GSLA,GROA,GSSL,GROS
      real(r8) :: GROH,GRMXB,GROLM,GROLC,GFNX,GRNDG,GRTWGM,GRTWTG
      real(r8) :: GRTLGL,GRTWTL,GRTWTN,GRTWTP,GRTWTM,HTNODZ,HTBR
      real(r8) :: HTSTK,HTLF,HTLFL,HTLFU,HTLFB,HTSTKX,PARTS,PARTX
      real(r8) :: PTRT,PARX,PARJ,PPOOLB,PADDB,PPOOLD,PPDX,PPOOLM
      real(r8) :: PADDN,PADD2,PADD1,PPOOLT,PTSHTR,PPOOLS,PPOOLG
      real(r8) :: PPOLNG,PPOOLX,PPOLNX,RTK,RS,RSL,RCO2C,RMNCS,RCO2X
      real(r8) :: RCO2Y,RCO2G,RCO2T,RCO2CM,RCO2XM,RCO2YM,RCO2GM
      real(r8) :: RCO2TM,RCCC,RCCN,RCCP,RCO2V,RCCL,RCZL,RCPL,RCCS
      real(r8) :: RCZS,RCPS,RCSC,RCSN,RCSP,RCCK,RCZK,RCPK,RSTK,RCNDL
      real(r8) :: RMNDL,RXNDL,RGNDL,RSNDL,RGN2P,RGN2F,RUPNFB,RXNDLC
      real(r8) :: RXNDLN,RXNDLP,RDNDLC,RDNDLN,RDNDLP,RCNDLC,RCNDLN
      real(r8) :: RCNDLP,RGNDG,RXNSNC,RXNSNN,RXNSNP,RDNSNC,RDNSNN
      real(r8) :: RDNSNP,RCNSNC,RCNSNN,RCNSNP,RTDPP,RTDPS,RTSKP
      real(r8) :: RTSKS,RSCS2,RTLGL,RTLGZ,RMNCR,RCO2RM,RCO2R,RCCR
      real(r8) :: RCZR,RCPR,RTN2X,RTN2Y,RTDP1X,RSCS1,RTLGX,RTLGT
      real(r8) :: RTVL,RTAR,RCNDLM,RXNDLM,RGNDLM,RSNDLM,STK,SNCR
      real(r8) :: SNCRM,SLA,SSL,SNL,SNCZ,SNCX,SNCF,SNCT,SNCLF,SNCSH
      real(r8) :: SET,SPNDLI,SPNDX,TKCM,TKSM,TOTAL,TLGLF,TFRCO2
      real(r8) :: UPNH4B,UPPO4B,UPNH4R,UPPO4R,VL,VGROX,VG,VA,VOLWPX
      real(r8) :: WTSHTZ,WHVSBL,WTHTH0,WTHNH0,WTHPH0,WTHTH1,WTHNH1
      real(r8) :: WTHPH1,WTHTH2,WTHNH2,WTHPH2,WTHTH3,WTHNH3,WTHPH3
      real(r8) :: WTHTH4,WTHNH4,WTHPH4,WTHTR1,WTHNR1,WTHPR1,WTHTR2
      real(r8) :: WTHNR2,WTHPR2,WTHTR3,WTHNR3,WTHPR3,WTHTR4,WTHNR4
      real(r8) :: WTHPR4,WTHTX0,WTHNX0,WTHPX0,WTHTX1,WTHNX1,WTHPX1
      real(r8) :: WTHTX2,WTHNX2,WTHPX2,WTHTX3,WTHNX3,WTHPX3,WTHTX4
      real(r8) :: WTHNX4,WTHPX4,WTHTG,WTHNG,WTHPG,WTSHXN,WTRTM,WFNSP
      real(r8) :: WTPLTT,WTRTRX,WTPLTX,WVSTBX,WTRTTX,WTRSBX,WTRVCX
      real(r8) :: WTLSB1,WTNDB1,WTLSBT,WTRTX,WTRTZ,WTRTLX,WTRTTT
      real(r8) :: WTRTT,WTRTD1,WTNDL1,WTRTDT,WTSTKT,WTRSVT,WTRSNT
      real(r8) :: WTRSPT,WTRSVD,WTRSND,WTRSPD,WTRTD2,WTLSBX,WTLSBB
      real(r8) :: WTRTLR,WHVSTT,WHVSLF,WHVHSH,WHVEAH,WHVGRH,WHVSCP
      real(r8) :: WHVSTH,WHVRVH,WHVSLX,WHVSLY,WHVSCL,WHVSNL,WHVXXX
      real(r8) :: WHVSSX,WTSHTT,WHVSHX,WHVSHY,WHVSHH,WHVSCS,WHVSNS
      real(r8) :: WHVHSX,WHVHSY,WHVEAX,WHVEAY,WHVGRX,WHVGRY,WHVSNP
      real(r8) :: WHVSKX,WHVSTX,WHVSTY,WHVRVX,WHVRVY,WTNDG,WTNDNG
      real(r8) :: WTNDPG,WGLFGX,WGSHGX,WGLFGY,WGSHGY,WGLFG,WGLFNG
      real(r8) :: WGLFPG,WHVSBS,WHVSCX,WHVSNX,WVPLT,WHVSTD,WTHTR0
      real(r8) :: WTHNR0,WTHPR0,WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2
      real(r8) :: WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4,WTHTHT,WTHTRT
      real(r8) :: WTHNHT,WTHNRT,WTHPHT,WTHPRT,WTHTXT,WTHNXT,WTHPXT
      real(r8) :: XHVST,XRTN1,XKSNC,XFRC,XFRN,XFRP,XFRN1,XFRP1,XLGLF
      real(r8) :: XLOCM,XLOCC,XLOCN,XLOCP,XFRCX,XFRNX,XFRPX,XFRW
      real(r8) :: XFRD,XHVSN,XHVSP,YLGLF,YARLF,YWGLF,YWGLFN,YWGLFP
      real(r8) :: ZPOOLB,ZADDB,ZADDBM,ZPOOLD,ZSTK,ZNPGN,ZNPGP,ZPGRN
      real(r8) :: ZPGRP,ZPOOLM,ZADDN,ZADD2M,ZADD2,ZADD1M,ZADD1,ZPOOLT
      real(r8) :: ZPOOLS,ZPOOLG,ZPOLNG,ZPOOLX,ZPOLNX

      integer :: IFLGZ,IFLGY,IDTHRN,IDTHY,K,KNOD,KK,KX,K1,K2
      integer :: KVSTGX,KSNC,KN,KVSTG1,L,LU,LL,LHTLFU,LHTLFL,LHTBRU
      integer :: LHTBRL,LZ,L1,LX,MNNOD,NX,NY,NZ,NN,NX1
      integer :: M,MXNOD,NY1,NR,N,NB,NNOD1,NBY,NBL,NBX,NBK
!     DATA TC4,TLK/0.0,0.0/
      REAL(r8) :: TFN5,WFNG,WFNC,WFNS,WFNSG,WFN4,WFNB &
      ,WFNR,WFNRG,FSNC2

      DATA RCCZ/0.167_r8,0.167_r8,0.167_r8,0.056_r8/
      DATA RCCY/0.333_r8,0.333_r8,0.167,0.333_r8/
      DATA RCCX/0.417_r8,0.833_r8,0.833,0.833_r8/
      DATA RCCQ/0.417_r8,0.833_r8,0.833,0.833_r8/
!
!     RTSK=relative primary root sink strength 0.25=shallow,4.0=deep root profile
!     FXRN=rate constant for plant-bacteria nonstructl C,N,P exchange (h-1)
!     FXFB=rate constant for leaf-storage nonstructl C,N,P exchange (h-1)
!     FXFR=rate constant for root-storage nonstructl C,N,P exchange (h-1)
!     FPART=parameter for allocating nonstructural C to shoot organs
!     FXSH,FXRT=shoot-root partitioning of storage C during leafout
!     FRSV=rate constant for remobiln of storage C,N,P during leafout (h-1)
!     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P exchange (h-1)
!     EFIRE=combustion  of N,P relative to C
!     PSILY=canopy water potential below which leafoff is induced (MPa)
!     FLG4Y=number of hours after physiol maturity required for senescence
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!     FVRN=fraction of hours required for leafoff to initiate remobilization
!
  DATA RTSK/0.25_r8,1.0_r8,4.0_r8,4.0_r8/
  DATA FXRN/0.25_r8,0.125_r8,0.0625_r8,0.225_r8,0.075_r8,0.025_r8/
  DATA FXFB/1.0E-02_r8,1.0E-02_r8,1.0E-05_r8,5.0E-05_r8/
  DATA FXFR/1.0E-02_r8,1.0E-02_r8,1.0E-05_r8,5.0E-05_r8/
  DATA FPART1/1.00_r8/,FPART2/0.40_r8/
  DATA FXSH/0.50_r8,0.75_r8/,FXRT/0.50_r8,0.25_r8/
  DATA FRSV/0.025_r8,0.025_r8,0.001_r8,0.001_r8/
  DATA FXFY/0.025_r8,0.005_r8/,FXFZ/0.25_r8,0.05_r8/
  DATA EFIRE/0.917_r8,0.167_r8/
  DATA PSILY/-200.0_r8,-2.0_r8,-2.0_r8,-2.0_r8/
  DATA FLG4Y/360.0_r8,1440.0_r8,720.0_r8,720.0_r8/
  DATA ATRPX/68.96_r8,276.9_r8/,GVMX/0.010_r8,0.0025_r8/
  DATA FVRN/0.75_r8,0.5_r8,0.5_r8,0.5_r8,0.5_r8,0.5_r8/

  public :: grosub
  contains

  SUBROUTINE grosub(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

! begin_execution
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
  call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)
!
!     INITIALIZE SENESCENCE ARRAYS
!
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      DO 9980 NZ=1,NP0(NY,NX)
        DO 1 L=0,NJ(NY,NX)
          DO K=0,1
            DO M=1,4
              CSNC(M,K,L,NZ,NY,NX)=0._r8
              ZSNC(M,K,L,NZ,NY,NX)=0._r8
              PSNC(M,K,L,NZ,NY,NX)=0._r8
            ENDDO
          ENDDO
1       CONTINUE
        HCSNC(NZ,NY,NX)=0._r8
        HZSNC(NZ,NY,NX)=0._r8
        HPSNC(NZ,NY,NX)=0._r8
        CNET(NZ,NY,NX)=0._r8
        UPNFC(NZ,NY,NX)=0._r8
        ZCX(NZ,NY,NX)=ZC(NZ,NY,NX)
        ZC(NZ,NY,NX)=0._r8
9980  CONTINUE
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
      DO 9985 NZ=1,NP(NY,NX)
!     IF(J.EQ.INT(ZNOON(NY,NX)))THEN
        XHVST=1.0_r8
        WHVSBL=0._r8
        WTHTH0=0._r8
        WTHNH0=0._r8
        WTHPH0=0._r8
        WTHTH1=0._r8
        WTHNH1=0._r8
        WTHPH1=0._r8
        WTHTH2=0._r8
        WTHNH2=0._r8
        WTHPH2=0._r8
        WTHTH3=0._r8
        WTHNH3=0._r8
        WTHPH3=0._r8
        WTHTH4=0._r8
        WTHNH4=0._r8
        WTHPH4=0._r8
        WTHTR1=0._r8
        WTHNR1=0._r8
        WTHPR1=0._r8
        WTHTR2=0._r8
        WTHNR2=0._r8
        WTHPR2=0._r8
        WTHTR3=0._r8
        WTHNR3=0._r8
        WTHPR3=0._r8
        WTHTR4=0._r8
        WTHNR4=0._r8
        WTHPR4=0._r8
        WTHTX0=0._r8
        WTHNX0=0._r8
        WTHPX0=0._r8
        WTHTX1=0._r8
        WTHNX1=0._r8
        WTHPX1=0._r8
        WTHTX2=0._r8
        WTHNX2=0._r8
        WTHPX2=0._r8
        WTHTX3=0._r8
        WTHNX3=0._r8
        WTHPX3=0._r8
        WTHTX4=0._r8
        WTHNX4=0._r8
        WTHPX4=0._r8
        WTHTG=0._r8
        WTHNG=0._r8
        WTHPG=0._r8
!     ENDIF
!     IF(NX.EQ.4.AND.NY.EQ.4.AND.NZ.EQ.2)THEN
!     WRITE(*,2328)'IFLGC1',I,J,NZ,IFLGC(NZ,NY,NX)
!    2,IDTHP(NZ,NY,NX),IDTHR(NZ,NY,NX)
!2328  FORMAT(A8,10I4)
!     ENDIF
        call GrowPlant(I,J,NZ,NY,NX)
!
!     HARVEST STANDING DEAD
!
        call RemoveBiomassByDisturbance(I,J,NZ,NY,NX)
9985  CONTINUE
!
!     TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS
      call LiveDeadTransformation(I,J,NY,NX)
9990  CONTINUE
9995  CONTINUE
      RETURN
      END subroutine grosub
!------------------------------------------------------------------------------------------

      subroutine PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)
      implicit none
      integer, intent(in) :: I, J
      integer, intent(in) :: NHW,NHE,NVN,NVS
!     begin_execution

      DO 2995 NX=NHW,NHE
      DO 2990 NY=NVN,NVS
      DO 2985 NZ=1,NP(NY,NX)
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     LSG=landscape grazing section number
!     WTSHTZ,WTSHTA=total,average biomass in landscape grazing section
!
      IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      WTSHTZ=0
      NN=0
      DO 1995 NX1=NHW,NHE
      DO 1990 NY1=NVN,NVS
      IF(LSG(NZ,NY1,NX1).EQ.LSG(NZ,NY,NX))THEN
      IF(IFLGC(NZ,NY1,NX1).EQ.1)THEN
      WTSHTZ=WTSHTZ+WTSHT(NZ,NY1,NX1)
      NN=NN+1
      ENDIF
      ENDIF
1990  CONTINUE
1995  CONTINUE
      IF(NN.GT.0)THEN
      WTSHTA(NZ,NY,NX)=WTSHTZ/NN
      ELSE
      WTSHTA(NZ,NY,NX)=WTSHT(NZ,NY,NX)
      ENDIF
      ENDIF
2985  CONTINUE
2990  CONTINUE
2995  CONTINUE
      end subroutine PrepLandscapeGrazing
!------------------------------------------------------------------------------------------

      subroutine RemoveBiomassByDisturbance(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     THIN=thinning:fraction of population removed
!     FHVST=fraction of standing dead mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     WTHTH4,WTHNH4,WTHPH4=harvested standing dead C,N,P
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
!
      IF(IHVST(NZ,I,NY,NX).GE.0)THEN
      IF(J.EQ.INT(ZNOON(NY,NX)).AND.IHVST(NZ,I,NY,NX).NE.4 &
      .AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      FHVST=AMAX1(0.0,1.0-EHVST(1,4,NZ,I,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSH=AMAX1(0.0,1.0-EHVST(1,4,NZ,I,NY,NX)*THIN(NZ,I,NY,NX))
      ELSE
      FHVSH=FHVST
      ENDIF
      ENDIF
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WTSTG(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSTD=HVST(NZ,I,NY,NX)*THIN(NZ,I,NY,NX)*0.45/24.0 &
      *AREA(3,NU(NY,NX),NY,NX)*EHVST(1,4,NZ,I,NY,NX)
      FHVST=AMAX1(0.0,1.0-WHVSTD/WTSTG(NZ,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
      ENDIF
      ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
      ENDIF
      DO 6475 M=1,4
      WTHTH4=WTHTH4+(1.0-FHVSH)*WTSTDG(M,NZ,NY,NX)
      WTHNH4=WTHNH4+(1.0-FHVSH)*WTSTDN(M,NZ,NY,NX)
      WTHPH4=WTHPH4+(1.0-FHVSH)*WTSTDP(M,NZ,NY,NX)
      WTHTX4=WTHTX4+(FHVSH-FHVST)*WTSTDG(M,NZ,NY,NX)
      WTHNX4=WTHNX4+(FHVSH-FHVST)*WTSTDN(M,NZ,NY,NX)
      WTHPX4=WTHPX4+(FHVSH-FHVST)*WTSTDP(M,NZ,NY,NX)
      WTSTDG(M,NZ,NY,NX)=FHVST*WTSTDG(M,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=FHVST*WTSTDN(M,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=FHVST*WTSTDP(M,NZ,NY,NX)
6475  CONTINUE
!
      call ApplyDisturbanceBiomRemoval(I,J,NZ,NY,NX)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
      call TotalBiomRemovalByDisturbance(I,J,NZ,NY,NX)
!
!     ABOVE-GROUND LITTERFALL FROM HARVESTING
!
      call LiterfallByDisturbance(I,J,NZ,NY,NX)

      ZEROP(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)
      ZEROQ(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      ZEROL(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)*1.0E+06
      ENDIF
      end subroutine RemoveBiomassByDisturbance
!------------------------------------------------------------------------------------------

      subroutine LiveDeadTransformation(I,J,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NY,NX
!     begin_execution

      DO 9975 NZ=1,NP0(NY,NX)
!
!     ACTIVATE DORMANT SEEDS
!
      DO 205 NB=1,NBR(NZ,NY,NX)
      IF(IFLGI(NZ,NY,NX).EQ.1)THEN
      IF(IFLGE(NB,NZ,NY,NX).EQ.0 &
      .AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
      IDAY0(NZ,NY,NX)=I
      IYR0(NZ,NY,NX)=IYRC
      SDPTHI(NZ,NY,NX)=0.005_r8+CDPTHZ(0,NY,NX)
      IFLGI(NZ,NY,NX)=0
      ENDIF
      ENDIF
205   CONTINUE
!
!     LITTERFALL FROM STANDING DEAD
!
!     XFRC,XFRN,XFRP=litterfall from standing dead
!     TFN3=temperature function for canopy growth
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall
!
      DO 6235 M=1,4
      XFRC=1.5814E-05*TFN3(NZ,NY,NX)*WTSTDG(M,NZ,NY,NX)
      XFRN=1.5814E-05*TFN3(NZ,NY,NX)*WTSTDN(M,NZ,NY,NX)
      XFRP=1.5814E-05*TFN3(NZ,NY,NX)*WTSTDP(M,NZ,NY,NX)
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+XFRC
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+XFRN
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+XFRP
      ELSE
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+XFRC
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+XFRN
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+XFRP
      ENDIF
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)-XFRC
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)-XFRN
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)-XFRP
6235  CONTINUE
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LITTERFALL
!
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
!
      DO 6430 M=1,4
      DO K=0,1
      TCSN0(NZ,NY,NX)=TCSN0(NZ,NY,NX)+CSNC(M,K,0,NZ,NY,NX)
      TZSN0(NZ,NY,NX)=TZSN0(NZ,NY,NX)+ZSNC(M,K,0,NZ,NY,NX)
      TPSN0(NZ,NY,NX)=TPSN0(NZ,NY,NX)+PSNC(M,K,0,NZ,NY,NX)
      DO 8955 L=0,NJ(NY,NX)
      HCSNC(NZ,NY,NX)=HCSNC(NZ,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
      HZSNC(NZ,NY,NX)=HZSNC(NZ,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
      HPSNC(NZ,NY,NX)=HPSNC(NZ,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
      TCSNC(NZ,NY,NX)=TCSNC(NZ,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
      TZSNC(NZ,NY,NX)=TZSNC(NZ,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
      TPSNC(NZ,NY,NX)=TPSNC(NZ,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
8955  CONTINUE
      enddo
6430  CONTINUE
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
      WTSTG(NZ,NY,NX)=WTSTDG(1,NZ,NY,NX)+WTSTDG(2,NZ,NY,NX) &
      +WTSTDG(3,NZ,NY,NX)+WTSTDG(4,NZ,NY,NX)
      WTSTGN(NZ,NY,NX)=WTSTDN(1,NZ,NY,NX)+WTSTDN(2,NZ,NY,NX) &
      +WTSTDN(3,NZ,NY,NX)+WTSTDN(4,NZ,NY,NX)
      WTSTGP(NZ,NY,NX)=WTSTDP(1,NZ,NY,NX)+WTSTDP(2,NZ,NY,NX) &
      +WTSTDP(3,NZ,NY,NX)+WTSTDP(4,NZ,NY,NX)
!
!     PLANT C BALANCE = TOTAL C STATE VARIABLES + TOTAL
!     AUTOTROPHIC RESPIRATION + TOTAL LITTERFALL - TOTAL EXUDATION
!     - TOTAL CO2 FIXATION
!
!     BALC=PFT C balance
!     WTSHT,WTRT,WTND,WTRVC,WTSTG=PFT shoot,root,bacteria,storage,standing dead C
!     ZNPP=cumulative PFT NPP
!     TCSNC=cumulative PFT C litterfall
!     TCUPTK=cumulative PFT root-soil C exchange
!     RSETC=cumulative C balance from previous year
!     THVSTC=cumulative PFT C removed from ecosystem from previous year
!     HVSTC=total PFT C removed from ecosystem in current year
!     VCO2F,VCH4F=CO2,CH4 emission from disturbance
!
      ZNPP(NZ,NY,NX)=CARBN(NZ,NY,NX)+TCO2T(NZ,NY,NX)
      IF(IFLGC(NZ,NY,NX).EQ.1)THEN
      BALC(NZ,NY,NX)=WTSHT(NZ,NY,NX)+WTRT(NZ,NY,NX)+WTND(NZ,NY,NX) &
      +WTRVC(NZ,NY,NX)-ZNPP(NZ,NY,NX)+TCSNC(NZ,NY,NX)-TCUPTK(NZ,NY,NX) &
      -RSETC(NZ,NY,NX)+WTSTG(NZ,NY,NX)+THVSTC(NZ,NY,NX) &
      +HVSTC(NZ,NY,NX)-VCO2F(NZ,NY,NX)-VCH4F(NZ,NY,NX)
!     IF(NZ.EQ.2)THEN
!     WRITE(*,1111)'BALC',I,J,NX,NY,NZ,NRT(NZ,NY,NX)
!    2,BALC(NZ,NY,NX),WTSHT(NZ,NY,NX)
!    2,WTRT(NZ,NY,NX),WTND(NZ,NY,NX),WTRVC(NZ,NY,NX),ZNPP(NZ,NY,NX)
!    3,TCSNC(NZ,NY,NX),TCUPTK(NZ,NY,NX),RSETC(NZ,NY,NX)
!    3,WTSTG(NZ,NY,NX)
!    2,THVSTC(NZ,NY,NX),HVSTC(NZ,NY,NX),VCO2F(NZ,NY,NX)
!    3,VCH4F(NZ,NY,NX),CARBN(NZ,NY,NX),TCO2T(NZ,NY,NX)
!    3,((CSNC(M,1,L,NZ,NY,NX),M=1,4),L=0,NJ(NY,NX))
!    3,WTLF(NZ,NY,NX),WTSHE(NZ,NY,NX),WTSTK(NZ,NY,NX),WTRSV(NZ,NY,NX)
!    3,WTHSK(NZ,NY,NX),WTEAR(NZ,NY,NX),WTGR(NZ,NY,NX)
!    4,((CPOOLR(N,L,NZ,NY,NX),L=NU(NY,NX),NJ(NY,NX)),N=1,2)
!    5,(((WTRT1(N,L,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
!    2,L=NU(NY,NX),NJ(NY,NX)),N=1,2)
!    5,(((WTRT2(N,L,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
!    2,L=NU(NY,NX),NJ(NY,NX)),N=1,2)
!1111  FORMAT(A8,6I4,360F16.8)
!     ENDIF
!
!     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LITTERFALL
!     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
!
!     BALN=PFT N balance
!     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT shoot,root,bacteria,storage,standing dead N
!     TZSNC=cumulative PFT N litterfall
!     TZUPTK=cumulative PFT root-soil N exchange
!     TNH3C=cumulative NH3 exchange
!     RSETN=cumulative N balance from previous year
!     THVSTN=cumulative PFT N removed from ecosystem from previous year
!     HVSTN=total PFT N removed from ecosystem in current year
!     VNH3F,VN2OF=NH3,N2O emission from disturbance
!     TZUPFX=cumulative PFT N2 fixation
!
      BALN(NZ,NY,NX)=WTSHN(NZ,NY,NX)+WTRTN(NZ,NY,NX)+WTNDN(NZ,NY,NX) &
      +WTRVN(NZ,NY,NX)+TZSNC(NZ,NY,NX)-TZUPTK(NZ,NY,NX)-TNH3C(NZ,NY,NX) &
      -RSETN(NZ,NY,NX)+WTSTGN(NZ,NY,NX)+HVSTN(NZ,NY,NX)+THVSTN(NZ,NY,NX) &
      -VNH3F(NZ,NY,NX)-VN2OF(NZ,NY,NX)-TZUPFX(NZ,NY,NX)
!     IF(NZ.EQ.1)THEN
!     WRITE(*,1112)'BALN',I,J,NX,NY,NZ,BALN(NZ,NY,NX),WTSHN(NZ,NY,NX)
!    2,WTRTN(NZ,NY,NX),WTNDN(NZ,NY,NX),WTRVN(NZ,NY,NX),TZSNC(NZ,NY,NX)
!    3,TZUPTK(NZ,NY,NX),TNH3C(NZ,NY,NX),RSETN(NZ,NY,NX),HVSTN(NZ,NY,NX)
!    4,WTSTGN(NZ,NY,NX),WTLFN(NZ,NY,NX),WTSHEN(NZ,NY,NX)
!    5,WTSTKN(NZ,NY,NX),WTRSVN(NZ,NY,NX),WTHSKN(NZ,NY,NX)
!    3,WTEARN(NZ,NY,NX),WTGRNN(NZ,NY,NX),UPOMN(NZ,NY,NX),UPNH4(NZ,NY,NX)
!    2,UPNO3(NZ,NY,NX),VNH3F(NZ,NY,NX),VN2OF(NZ,NY,NX)
!    4,((RDFOMN(N,L,NZ,NY,NX),N=1,2),L=NU(NY,NX),NI(NZ,NY,NX))
!    4,((ZPOOLR(N,L,NZ,NY,NX),N=1,2),L=NU(NY,NX),NI(NZ,NY,NX))
!1112  FORMAT(A8,5I4,200F18.6)
!     ENDIF
!
!     PLANT P BALANCE = TOTAL P STATE VARIABLES + TOTAL P LITTERFALL
!     - TOTAL P UPTAKE FROM SOIL
!
!     BALP=PFT N balance
!     WTSHP,WTRTP,WTNDP,WTRVP,WTSTGP=PFT shoot,root,bacteria,storage,standing dead P
!     TPSNC=cumulative PFT P litterfall
!     TPUPTK=cumulative PFT root-soil P exchange
!     RSETP=cumulative P balance from previous year
!     THVSTP=cumulative PFT P removed from ecosystem from previous year
!     HVSTP=total PFT P removed from ecosystem in current year
!     VPO4F=PO4 emission from disturbance
!
      BALP(NZ,NY,NX)=WTSHP(NZ,NY,NX)+WTRTP(NZ,NY,NX)+WTNDP(NZ,NY,NX) &
      +WTRVP(NZ,NY,NX)+TPSNC(NZ,NY,NX)-TPUPTK(NZ,NY,NX) &
      -RSETP(NZ,NY,NX)+WTSTDP(1,NZ,NY,NX)+WTSTGP(NZ,NY,NX) &
      +HVSTP(NZ,NY,NX)+THVSTP(NZ,NY,NX)-VPO4F(NZ,NY,NX)
!     IF(NZ.EQ.4)THEN
!     WRITE(*,1112)'BALP',I,J,NX,NY,NZ,BALP(NZ,NY,NX),WTSHP(NZ,NY,NX)
!    2,WTRTP(NZ,NY,NX),WTNDP(NZ,NY,NX),WTRVP(NZ,NY,NX),TPSNC(NZ,NY,NX)
!    3,TPUPTK(NZ,NY,NX),RSETP(NZ,NY,NX)
!    4,WTSTDP(1,NZ,NY,NX),WTSTGP(NZ,NY,NX),HVSTP(NZ,NY,NX)
!    5,THVSTP(NZ,NY,NX),VPO4F(NZ,NY,NX)
!     ENDIF
      ENDIF
9975  CONTINUE
      end subroutine LiveDeadTransformation
!------------------------------------------------------------------------------------------

  subroutine GrowPlant(I,J,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX

! begin_execution
  IF(IFLGC(NZ,NY,NX).EQ.1)THEN

   IF(IDTHP(NZ,NY,NX).EQ.0.OR.IDTHR(NZ,NY,NX).EQ.0)THEN

      call StagePlantForGrowth(I,J,NZ,NY,NX)
!
!     CALCULATE GROWTH OF EACH BRANCH
!
!     WTLFB,WTSHEB,WTLSB=leaf,petiole,leaf+petiole mass
!     IDTHB=branch living flag: 0=alive,1=dead
!
      DO 105 NB=1,NBR(NZ,NY,NX)
      call GrowOneBranch(I,J,NB,NZ,NY,NX)
105   CONTINUE
!
!     ROOT GROWTH
!
      call RootBiochemistry(I,J,NZ,NY,NX)
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
      RTLGP(1,NG(NZ,NY,NX),NZ,NY,NX)=RTLGP(1,NG(NZ,NY,NX),NZ,NY,NX) &
      +SDLG(NZ,NY,NX)
      IF(DLYR(3,NG(NZ,NY,NX),NY,NX).GT.ZERO)THEN
      RTDNP(1,NG(NZ,NY,NX),NZ,NY,NX)=RTLGP(1,NG(NZ,NY,NX),NZ,NY,NX) &
      /DLYR(3,NG(NZ,NY,NX),NY,NX)
      ELSE
      RTDNP(1,NG(NZ,NY,NX),NZ,NY,NX)=0._r8
      ENDIF
      RTVL=RTVLP(1,NG(NZ,NY,NX),NZ,NY,NX) &
      +RTVLW(1,NG(NZ,NY,NX),NZ,NY,NX)+SDVL(NZ,NY,NX)*PP(NZ,NY,NX)
      RTVLP(1,NG(NZ,NY,NX),NZ,NY,NX)=PORT(1,NZ,NY,NX)*RTVL
      RTVLW(1,NG(NZ,NY,NX),NZ,NY,NX)=(1.0-PORT(1,NZ,NY,NX))*RTVL
      RTARP(1,NG(NZ,NY,NX),NZ,NY,NX)=RTARP(1,NG(NZ,NY,NX),NZ,NY,NX) &
      +SDAR(NZ,NY,NX)
      IF(IDTHRN.EQ.NRT(NZ,NY,NX).OR.(WTRVC(NZ,NY,NX) &
      .LE.ZEROL(NZ,NY,NX).AND.ISTYP(NZ,NY,NX).NE.0))THEN
      IDTHR(NZ,NY,NX)=1
      IDTHP(NZ,NY,NX)=1
      ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
!
      call RootNoduleBiomchemistry(I,J,NZ,NY,NX)

      call NonstructlBiomTransfer(I,J,NZ,NY,NX)
!
!
      call ComputeTotalBiom(NZ,NY,NX)
      ELSE
      HCUPTK(NZ,NY,NX)=UPOMC(NZ,NY,NX)
      HZUPTK(NZ,NY,NX)=UPOMN(NZ,NY,NX)+UPNH4(NZ,NY,NX)+UPNO3(NZ,NY,NX) &
      +UPNF(NZ,NY,NX)
      HPUPTK(NZ,NY,NX)=UPOMP(NZ,NY,NX)+UPH2P(NZ,NY,NX)+UPH1P(NZ,NY,NX)
      ENDIF
!
!     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
!
      call RemoveBiomByHarvest(I,J,NZ,NY,NX)
!
!     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
!
      call RemoveBiomByTillage(I,J,NZ,NY,NX)
!
!     RESET DEAD BRANCHES
      call ResetDeadBranch(I,J,NZ,NY,NX)
!
      call AccumulateStates(I,J,NZ,NY,NX)

      ENDIF
      end subroutine GrowPlant
!------------------------------------------------------------------------------------------

      subroutine ResetDeadBranch(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!
!     ZNOON=hour of solar noon
!     IDAY(1,=emergence date
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAYH,IYRH=day,year of harvesting
!     IYRC=current year
!     IDTHB=branch living flag: 0=alive,1=dead
!     GROUP=node number required for floral initiation
!     PSTGI=node number at floral initiation
!     PSTGF=node number at flowering
!     VSTG=number of leaves appeared
!     KVSTG=integer of most recent leaf number currently growing
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     IFLGA=flag for initializing leafout
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     ATRP=hourly leafout counter
!     FDBK,FDBKX=N,P feedback inhibition on C3 CO2 fixation
!     IFLGA,IFLGE=flag for initializing,enabling leafout
!     IFLGF=flag for enabling leafoff:0=enable,1=disable
!     IFLGQ=current hours after physl maturity until start of litterfall
!
      IF(J.EQ.INT(ZNOON(NY,NX)) &
      .AND.IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0 &
      .AND.(ISTYP(NZ,NY,NX).NE.0.OR.(I.GE.IDAYH(NZ,NY,NX) &
      .AND.IYRC.GE.IYRH(NZ,NY,NX))))THEN
      IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
      call LiterfallFromDeadBranches(I,J,NZ,NY,NX)

      IF(IDTHY.EQ.NBR(NZ,NY,NX))THEN
      IDTHP(NZ,NY,NX)=1
      NBT(NZ,NY,NX)=0
      WSTR(NZ,NY,NX)=0._r8
      IF(IFLGI(NZ,NY,NX).EQ.1)THEN
      NBR(NZ,NY,NX)=1
      ELSE
      NBR(NZ,NY,NX)=0
      ENDIF
      HTCTL(NZ,NY,NX)=0._r8
      VOLWOU=VOLWOU+VOLWP(NZ,NY,NX)
      UVOLO(NY,NX)=UVOLO(NY,NX)+VOLWP(NZ,NY,NX)
      VOLWP(NZ,NY,NX)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     ISTYP=growth habit:0=annual,1=perennial
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(WTRVC(NZ,NY,NX).LT.1.0E-04*WTRT(NZ,NY,NX) &
      .AND.ISTYP(NZ,NY,NX).NE.0)IDTHR(NZ,NY,NX)=1
      IF(ISTYP(NZ,NY,NX).EQ.0)IDTHR(NZ,NY,NX)=1
      IF(JHVST(NZ,I,NY,NX).NE.0)IDTHR(NZ,NY,NX)=1
      IF(PP(NZ,NY,NX).LE.0.0)IDTHR(NZ,NY,NX)=1
      IF(IDTHR(NZ,NY,NX).EQ.1)IDTHP(NZ,NY,NX)=1
      ENDIF
!
!     DEAD ROOTS
!
!
!     LITTERFALL FROM DEAD ROOTS
!
      call LiterfallFromDeadRoots(I,J,NZ,NY,NX)
!
      call LiterfallFromRootShootStorage(I,J,NZ,NY,NX)
      ENDIF
      end subroutine ResetDeadBranch
!------------------------------------------------------------------------------------------

      subroutine LiterfallFromRootShootStorage(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM SHOOT AT DEATH
!
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!     IFLGI=PFT initialization flag:0=no,1=yes
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
      IF(IDTHP(NZ,NY,NX).EQ.1.AND.IDTHR(NZ,NY,NX).EQ.1)THEN
      IF(IFLGI(NZ,NY,NX).EQ.0)THEN
      DO 6425 M=1,4
      DO 8825 NB=1,NBR(NZ,NY,NX)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*(CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX) &
      +CPOOLK(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)) &
      +CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1) &
      +WTNDB(NB,NZ,NY,NX)) &
      +CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1) &
      +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX))
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0) &
      +WTSHEB(NB,NZ,NY,NX)*FWODB(0))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*(ZPOOL(NB,NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX) &
      +WTRSBN(NB,NZ,NY,NX)) &
      +CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1) &
      +WTNDBN(NB,NZ,NY,NX)) &
      +CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1) &
      +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX))
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0) &
      +WTSHBN(NB,NZ,NY,NX)*FWODSN(0))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*(PPOOL(NB,NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX) &
      +WTRSBP(NB,NZ,NY,NX)) &
      +CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1) &
      +WTNDBP(NB,NZ,NY,NX)) &
      +CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1) &
      +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX))
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0) &
      +WTSHBP(NB,NZ,NY,NX)*FWODSP(0))
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX) &
      +CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX) &
      +CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX) &
      +CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)
      ELSE
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)
      ENDIF
8825  CONTINUE
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM ROOT AND STORGE AT DEATH
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
      DO 6415 L=NU(NY,NX),NJ(NY,NX)
      DO N=1,MY(NZ,NY,NX)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX) &
      *CPOOLR(N,L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX) &
      *ZPOOLR(N,L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX) &
      *PPOOLR(N,L,NZ,NY,NX)
      DO NR=1,NRT(NZ,NY,NX)
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
      *(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
      *(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
      *(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1)
      ENDDO
      ENDDO
6415  CONTINUE
      CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX)*FWOOD(0)
      ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX)*FWOODN(0)
      PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX)*FWOODP(0)
      CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX)*FWOOD(1)
      ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX)*FWOODN(1)
      PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX)*FWOODP(1)
6425  CONTINUE
!
      call ResetBranchRootStates(NZ,NY,NX)
      ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     LYRC=number of days in current year
!     IDAY0,IYR0=day,year of planting
!
      IF(ISTYP(NZ,NY,NX).NE.0.AND.JHVST(NZ,I,NY,NX).EQ.0)THEN
      IF(I.LT.LYRC)THEN
      IDAY0(NZ,NY,NX)=I+1
      IYR0(NZ,NY,NX)=IDATA(3)
      ELSE
      IDAY0(NZ,NY,NX)=1
      IYR0(NZ,NY,NX)=IDATA(3)+1
      ENDIF
      ENDIF
      ENDIF
      end subroutine LiterfallFromRootShootStorage
!------------------------------------------------------------------------------------------

      subroutine LiterfallFromDeadRoots(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     IDTHR=PFT root living flag: 0=alive,1=dead
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
      IF(IDTHR(NZ,NY,NX).EQ.1)THEN
      DO 8900 N=1,MY(NZ,NY,NX)
      DO 8895 L=NU(NY,NX),NJ(NY,NX)
      DO 6410 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX) &
      *CPOOLR(N,L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX) &
      *ZPOOLR(N,L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX) &
      *PPOOLR(N,L,NZ,NY,NX)
      DO  NR=1,NRT(NZ,NY,NX)
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
      *(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
      *(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
      *(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1)
      enddo
6410  CONTINUE
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-CO2A(N,L,NZ,NY,NX) &
      -CO2P(N,L,NZ,NY,NX)
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-OXYA(N,L,NZ,NY,NX) &
      -OXYP(N,L,NZ,NY,NX)
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-CH4A(N,L,NZ,NY,NX) &
      -CH4P(N,L,NZ,NY,NX)
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-Z2OA(N,L,NZ,NY,NX) &
      -Z2OP(N,L,NZ,NY,NX)
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-ZH3A(N,L,NZ,NY,NX) &
      -ZH3P(N,L,NZ,NY,NX)
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-H2GA(N,L,NZ,NY,NX) &
      -H2GP(N,L,NZ,NY,NX)
      CO2A(N,L,NZ,NY,NX)=0._r8
      OXYA(N,L,NZ,NY,NX)=0._r8
      CH4A(N,L,NZ,NY,NX)=0._r8
      Z2OA(N,L,NZ,NY,NX)=0._r8
      ZH3A(N,L,NZ,NY,NX)=0._r8
      H2GA(N,L,NZ,NY,NX)=0._r8
      CO2P(N,L,NZ,NY,NX)=0._r8
      OXYP(N,L,NZ,NY,NX)=0._r8
      CH4P(N,L,NZ,NY,NX)=0._r8
      Z2OP(N,L,NZ,NY,NX)=0._r8
      ZH3P(N,L,NZ,NY,NX)=0._r8
      H2GP(N,L,NZ,NY,NX)=0._r8
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RTLGA=average secondary root length
!
      DO 8870 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=0._r8
      WTRT1N(N,L,NR,NZ,NY,NX)=0._r8
      WTRT1P(N,L,NR,NZ,NY,NX)=0._r8
      WTRT2(N,L,NR,NZ,NY,NX)=0._r8
      WTRT2N(N,L,NR,NZ,NY,NX)=0._r8
      WTRT2P(N,L,NR,NZ,NY,NX)=0._r8
      RTWT1(N,NR,NZ,NY,NX)=0._r8
      RTWT1N(N,NR,NZ,NY,NX)=0._r8
      RTWT1P(N,NR,NZ,NY,NX)=0._r8
      RTLG1(N,L,NR,NZ,NY,NX)=0._r8
      RTLG2(N,L,NR,NZ,NY,NX)=0._r8
      RTN2(N,L,NR,NZ,NY,NX)=0._r8
8870  CONTINUE
      CPOOLR(N,L,NZ,NY,NX)=0._r8
      ZPOOLR(N,L,NZ,NY,NX)=0._r8
      PPOOLR(N,L,NZ,NY,NX)=0._r8
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
      RTLGA(N,L,NZ,NY,NX)=RTLGAX
!
!     LITTERFALL AND STATE VARIABLES FROM DEAD NODULES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
      IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
      DO 6420 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *WTNDL(L,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *WTNDLN(L,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *WTNDLP(L,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX)
6420  CONTINUE
      WTNDL(L,NZ,NY,NX)=0._r8
      WTNDLN(L,NZ,NY,NX)=0._r8
      WTNDLP(L,NZ,NY,NX)=0._r8
      CPOOLN(L,NZ,NY,NX)=0._r8
      ZPOOLN(L,NZ,NY,NX)=0._r8
      PPOOLN(L,NZ,NY,NX)=0._r8
      ENDIF
8895  CONTINUE
8900  CONTINUE
!
!     RESET DEPTH VARIABLES OF DEAD ROOTS
!
!     NINR=deepest root layer
!     RTDP1=primary root depth from soil surface
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!
      DO 8795 NR=1,NRT(NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=NG(NZ,NY,NX)
      DO 8790 N=1,MY(NZ,NY,NX)
      RTDP1(N,NR,NZ,NY,NX)=SDPTH(NZ,NY,NX)
      RTWT1(N,NR,NZ,NY,NX)=0._r8
      RTWT1N(N,NR,NZ,NY,NX)=0._r8
      RTWT1P(N,NR,NZ,NY,NX)=0._r8
8790  CONTINUE
8795  CONTINUE
      NIX(NZ,NY,NX)=NG(NZ,NY,NX)
      NRT(NZ,NY,NX)=0
      ENDIF
      end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

      subroutine LiterfallFromDeadBranches(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution

      DO 8845 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.1)THEN
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0._r8
      VSTG(NB,NZ,NY,NX)=0._r8
      VSTGX(NB,NZ,NY,NX)=0._r8
      KLEAF(NB,NZ,NY,NX)=1
      KVSTG(NB,NZ,NY,NX)=1
      TGSTGI(NB,NZ,NY,NX)=0._r8
      TGSTGF(NB,NZ,NY,NX)=0._r8
      VRNS(NB,NZ,NY,NX)=0._r8
      VRNF(NB,NZ,NY,NX)=0._r8
      VRNY(NB,NZ,NY,NX)=0._r8
      VRNZ(NB,NZ,NY,NX)=0._r8
      ATRP(NB,NZ,NY,NX)=0._r8
      FLG4(NB,NZ,NY,NX)=0._r8
      FDBK(NB,NZ,NY,NX)=1.0_r8
      FDBKX(NB,NZ,NY,NX)=1.0_r8
      IFLGA(NB,NZ,NY,NX)=0
      IFLGE(NB,NZ,NY,NX)=1
      IFLGF(NB,NZ,NY,NX)=0
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
      NBTB(NB,NZ,NY,NX)=0
      DO 8850 M=1,10
      IDAY(M,NB,NZ,NY,NX)=0
8850  CONTINUE
!
!     LITTERFALL FROM DEAD BRANCHES
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      DO 6405 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*CPOLNB(NB,NZ,NY,NX) &
      +CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1) &
      +WTNDB(NB,NZ,NY,NX)) &
      +CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1) &
      +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX))
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0) &
      +WTSHEB(NB,NZ,NY,NX)*FWODB(0))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*ZPOLNB(NB,NZ,NY,NX) &
      +CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1) &
      +WTNDBN(NB,NZ,NY,NX)) &
      +CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1) &
      +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX))
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0) &
      +WTSHBN(NB,NZ,NY,NX)*FWODSN(0))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*PPOLNB(NB,NZ,NY,NX) &
      +CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1) &
      +WTNDBP(NB,NZ,NY,NX)) &
      +CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1) &
      +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX))
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0) &
      +WTSHBP(NB,NZ,NY,NX)*FWODSP(0))
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX) &
      +CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX) &
      +CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX) &
      +CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)
      ELSE
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)
      ENDIF
6405  CONTINUE
!
!     RECOVER NON-STRUCTURAL C,N,P FROM BRANCH TO
!     SEASONAL STORAGE RESERVES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX) &
      +CPOOLK(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+PPOOL(NB,NZ,NY,NX)
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      DO 6406 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*WTRSBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*WTRSBP(NB,NZ,NY,NX)
6406  CONTINUE
      ELSE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
      ENDIF
!
      call ResetDeadRootStates(NB,NZ,NY,NX)

      IDTHY=IDTHY+1
      ENDIF
8845  CONTINUE
      end subroutine LiterfallFromDeadBranches
!------------------------------------------------------------------------------------------

      subroutine RemoveBiomByTillage(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     ZNOON=hour of solar noon
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IDAY0,IYR0=day,year of planting
!     IYRC=current year
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!     PPX,PP=PFT population per m2,grid cell
!     FRADP=fraction of radiation received by each PFT canopy
!     VHCPC=canopy heat capacity
!
      IF(J.EQ.INT(ZNOON(NY,NX)).AND.(IBTYP(NZ,NY,NX).EQ.0 &
      .OR.IGTYP(NZ,NY,NX).LE.1).AND.(I.NE.IDAY0(NZ,NY,NX) &
      .OR.IYRC.NE.IYR0(NZ,NY,NX)))THEN
      IF(ITILL(I,NY,NX).LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.IDAY0(NZ,NY,NX).OR.IYRC.GT.IYR0(NZ,NY,NX))THEN
      XHVST=XCORP(NY,NX)
      PPX(NZ,NY,NX)=PPX(NZ,NY,NX)*XHVST
      PP(NZ,NY,NX)=PP(NZ,NY,NX)*XHVST
      FRADP(NZ,NY,NX)=FRADP(NZ,NY,NX)*XHVST
      VHCPC(NZ,NY,NX)=VHCPC(NZ,NY,NX)*XHVST
      WTLS(NZ,NY,NX)=0._r8
      WVSTK(NZ,NY,NX)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
      DO 8975 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(PP(NZ,NY,NX).LE.0.0)IDTHB(NB,NZ,NY,NX)=1
!
!     LITTERFALL FROM BRANCHES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      DO 6380 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPC(0,M,NZ,NY,NX)*(CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX) &
      +CPOOLK(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)) &
      +CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1) &
      +WTNDB(NB,NZ,NY,NX)) &
      +CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1) &
      +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)))
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0) &
      +WTSHEB(NB,NZ,NY,NX)*FWODB(0))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPN(0,M,NZ,NY,NX)*(ZPOOL(NB,NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX) &
      +WTRSBN(NB,NZ,NY,NX)) &
      +CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1) &
      +WTNDBN(NB,NZ,NY,NX)) &
      +CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1) &
      +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)))
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0) &
      +WTSHBN(NB,NZ,NY,NX)*FWODSN(0))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPP(0,M,NZ,NY,NX)*(PPOOL(NB,NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX) &
      +WTRSBP(NB,NZ,NY,NX)) &
      +CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1) &
      +WTNDBP(NB,NZ,NY,NX)) &
      +CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1) &
      +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)))
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0) &
      +WTSHBP(NB,NZ,NY,NX)*FWODSP(0))
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(5,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)*FWOOD(0)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(5,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)*FWOODN(0)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(5,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)*FWOODP(0)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)*FWOOD(1)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)*FWOODN(1)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)*FWOODP(1)
6380  CONTINUE
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)*XHVST
      CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX)*XHVST
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)*XHVST
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)*XHVST
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)*XHVST
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)*XHVST
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)*XHVST
      WTSHTB(NB,NZ,NY,NX)=WTSHTB(NB,NZ,NY,NX)*XHVST
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)*XHVST
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)*XHVST
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)*XHVST
      WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)*XHVST
      WVSTKB(NB,NZ,NY,NX)=WVSTKB(NB,NZ,NY,NX)*XHVST
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)*XHVST
      WTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)*XHVST
      WTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)*XHVST
      WTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)*XHVST
      WTSHTN(NB,NZ,NY,NX)=WTSHTN(NB,NZ,NY,NX)*XHVST
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)*XHVST
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)*XHVST
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)*XHVST
      WTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)*XHVST
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)*XHVST
      WTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)*XHVST
      WTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)*XHVST
      WTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)*XHVST
      WTSHTP(NB,NZ,NY,NX)=WTSHTP(NB,NZ,NY,NX)*XHVST
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)*XHVST
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)*XHVST
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)*XHVST
      WTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)*XHVST
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)*XHVST
      WTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)*XHVST
      WTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)*XHVST
      WTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)*XHVST
      GRNXB(NB,NZ,NY,NX)=GRNXB(NB,NZ,NY,NX)*XHVST
      GRNOB(NB,NZ,NY,NX)=GRNOB(NB,NZ,NY,NX)*XHVST
      GRWTB(NB,NZ,NY,NX)=GRWTB(NB,NZ,NY,NX)*XHVST
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)*XHVST
      WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX) &
      +WTSHEB(NB,NZ,NY,NX))
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
      WTSTXB(NB,NZ,NY,NX)=WTSTXB(NB,NZ,NY,NX)*XHVST
      WTSTXN(NB,NZ,NY,NX)=WTSTXN(NB,NZ,NY,NX)*XHVST
      WTSTXP(NB,NZ,NY,NX)=WTSTXP(NB,NZ,NY,NX)*XHVST
      WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      DO 8970 K=0,25
      IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)*XHVST
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)*XHVST
      CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)*XHVST
      HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)*XHVST
      ENDIF
      ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)*XHVST
      WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)*XHVST
      WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)*XHVST
!     HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)*XHVST
      WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)*XHVST
      WSSHE(K,NB,NZ,NY,NX)=WSSHE(K,NB,NZ,NY,NX)*XHVST
!     HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX)*XHVST
!     HTNODX(K,NB,NZ,NY,NX)=HTNODX(K,NB,NZ,NY,NX)*XHVST
      WGNODE(K,NB,NZ,NY,NX)=WGNODE(K,NB,NZ,NY,NX)*XHVST
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)*XHVST
      WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)*XHVST
      WGNODN(K,NB,NZ,NY,NX)=WGNODN(K,NB,NZ,NY,NX)*XHVST
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)*XHVST
      WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)*XHVST
      WGNODP(K,NB,NZ,NY,NX)=WGNODP(K,NB,NZ,NY,NX)*XHVST
      DO 8965 L=1,JC
      ARLFL(L,K,NB,NZ,NY,NX)=ARLFL(L,K,NB,NZ,NY,NX)*XHVST
      WGLFL(L,K,NB,NZ,NY,NX)=WGLFL(L,K,NB,NZ,NY,NX)*XHVST
      WGLFLN(L,K,NB,NZ,NY,NX)=WGLFLN(L,K,NB,NZ,NY,NX)*XHVST
      WGLFLP(L,K,NB,NZ,NY,NX)=WGLFLP(L,K,NB,NZ,NY,NX)*XHVST
8965  CONTINUE
8970  CONTINUE
      ENDIF
8975  CONTINUE
!
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
      VOLWPX=VOLWP(NZ,NY,NX)
      WVPLT=AMAX1(0.0_r8,WTLS(NZ,NY,NX)+WVSTK(NZ,NY,NX))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
      VOLWP(NZ,NY,NX)=1.0E-06_r8*WVPLT/FDM
      VOLWOU=VOLWOU+VOLWPX-VOLWP(NZ,NY,NX)
      UVOLO(NY,NX)=UVOLO(NY,NX)+VOLWPX-VOLWP(NZ,NY,NX)
!
!     TERMINATE ROOTS IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     PP=PFT population
!     IDTHR,IDTHP=PFT root,shoot living flag: 0=alive,1=dead
!     IDTH=PFT living flag: 0=alive,1=dead
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     IDAYH,IYRH=day,year of harvesting
!     IYRC=current year
!
      IF(PP(NZ,NY,NX).LE.0.0)THEN
      IDTHR(NZ,NY,NX)=1
      IDTHP(NZ,NY,NX)=1
      IDTH(NZ,NY,NX)=1
      JHVST(NZ,I,NY,NX)=1
      IDAYH(NZ,NY,NX)=I
      IYRH(NZ,NY,NX)=IYRC
      ENDIF
!
!     LITTERFALL FROM ROOTS DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!
      DO 8985 N=1,MY(NZ,NY,NX)
      DO 8980 L=NU(NY,NX),NJ(NY,NX)
      DO 6385 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(0,M,NZ,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(0,M,NZ,NY,NX)*ZPOOLR(N,L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(0,M,NZ,NY,NX)*PPOOLR(N,L,NZ,NY,NX)
      DO NR=1,NRT(NZ,NY,NX)
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(5,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
      +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(5,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
      +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(5,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
      +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPC(4,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
      +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPN(4,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
      +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *CFOPP(4,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
      +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1)
      ENDDO
6385  CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(1.0-XHVST) &
      *(CO2A(N,L,NZ,NY,NX)+CO2P(N,L,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(1.0-XHVST) &
      *(OXYA(N,L,NZ,NY,NX)+OXYP(N,L,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(1.0-XHVST) &
      *(CH4A(N,L,NZ,NY,NX)+CH4P(N,L,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(1.0-XHVST) &
      *(Z2OA(N,L,NZ,NY,NX)+Z2OP(N,L,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(1.0-XHVST) &
      *(ZH3A(N,L,NZ,NY,NX)+ZH3P(N,L,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(1.0-XHVST) &
      *(H2GA(N,L,NZ,NY,NX)+H2GP(N,L,NZ,NY,NX))
      CO2A(N,L,NZ,NY,NX)=XHVST*CO2A(N,L,NZ,NY,NX)
      OXYA(N,L,NZ,NY,NX)=XHVST*OXYA(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=XHVST*CH4A(N,L,NZ,NY,NX)
      Z2OA(N,L,NZ,NY,NX)=XHVST*Z2OA(N,L,NZ,NY,NX)
      ZH3A(N,L,NZ,NY,NX)=XHVST*ZH3A(N,L,NZ,NY,NX)
      H2GA(N,L,NZ,NY,NX)=XHVST*H2GA(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=XHVST*CO2P(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=XHVST*OXYP(N,L,NZ,NY,NX)
      CH4P(N,L,NZ,NY,NX)=XHVST*CH4P(N,L,NZ,NY,NX)
      Z2OP(N,L,NZ,NY,NX)=XHVST*Z2OP(N,L,NZ,NY,NX)
      ZH3P(N,L,NZ,NY,NX)=XHVST*ZH3P(N,L,NZ,NY,NX)
      H2GP(N,L,NZ,NY,NX)=XHVST*H2GP(N,L,NZ,NY,NX)
!
!     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
      DO 8960 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*XHVST
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*XHVST
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*XHVST
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*XHVST
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*XHVST
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*XHVST
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*XHVST
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*XHVST
      RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*XHVST
8960  CONTINUE
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*XHVST
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*XHVST
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*XHVST
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*XHVST
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*XHVST
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*XHVST
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*XHVST
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*XHVST
      RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*XHVST
      RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*XHVST
      RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*XHVST
      RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*XHVST
      RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*XHVST
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)*XHVST
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)*XHVST
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)*XHVST
!
!     LITTERFALL AND STATE VARIABLES FOR NODULES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
      IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
      DO 6395 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPC(4,M,NZ,NY,NX)*WTNDL(L,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX))
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPN(4,M,NZ,NY,NX)*WTNDLN(L,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX))
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPP(4,M,NZ,NY,NX)*WTNDLP(L,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX))
6395  CONTINUE
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*XHVST
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*XHVST
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*XHVST
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*XHVST
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*XHVST
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*XHVST
      ENDIF
8980  CONTINUE
8985  CONTINUE
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
!     DURING TILLAGE
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
      DO 6400 M=1,4
      CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(0)
      ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(0)
      PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(0)
      CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(1)
      ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(1)
      PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(1)
6400  CONTINUE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*XHVST
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*XHVST
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*XHVST
      ENDIF
      ENDIF
      ENDIF
      end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------

      subroutine RemoveBiomByHarvest(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
      IF((IHVST(NZ,I,NY,NX).GE.0.AND.J.EQ.INT(ZNOON(NY,NX)) &
      .AND.IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6) &
      .OR.(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6))THEN
!
!     ACCUMULATE ALL HARVESTED MATERIAL ABOVE CUTTING HEIGHT
!     ACCOUNTING FOR HARVEST EFFICIENCY ENTERED IN 'READQ'
!
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     PPX,PP=PFT population per m2,grid cell
!     THIN=thinning:fraction of population removed
!     CF=clumping factor
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     ARLFC,ARLFT=leaf area of combined canopy, canopy layer
!     ARLFR,ARLFY=leaf area harvested,remaining
!     ZL=height to bottom of each canopy layer
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(JHVST(NZ,I,NY,NX).NE.2)THEN
      PPX(NZ,NY,NX)=PPX(NZ,NY,NX)*(1.0-THIN(NZ,I,NY,NX))
      PP(NZ,NY,NX)=PP(NZ,NY,NX)*(1.0-THIN(NZ,I,NY,NX))
      ELSE
!     PPI(NZ,NY,NX)=AMAX1(1.0,0.5*(PPI(NZ,NY,NX)+GRNO(NZ,NY,NX)
!    2/AREA(3,NU(NY,NX),NY,NX)))
      PPX(NZ,NY,NX)=PPI(NZ,NY,NX)
      PP(NZ,NY,NX)=PPX(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      ENDIF
      IF(IHVST(NZ,I,NY,NX).EQ.3)THEN
      CF(NZ,NY,NX)=CF(NZ,NY,NX)*HVST(NZ,I,NY,NX)
      ENDIF
      IF(IHVST(NZ,I,NY,NX).LE.2.AND.HVST(NZ,I,NY,NX).LT.0.0)THEN
      ARLFY=(1.0-ABS(HVST(NZ,I,NY,NX)))*ARLFC(NY,NX)
      ARLFR=0._r8
      DO 9875 L=1,JC
      IF(ZL(L,NY,NX).GT.ZL(L-1,NY,NX) &
      .AND.ARLFT(L,NY,NX).GT.ZEROS(NY,NX) &
      .AND.ARLFR.LT.ARLFY)THEN
      IF(ARLFR+ARLFT(L,NY,NX).GT.ARLFY)THEN
      HVST(NZ,I,NY,NX)=ZL(L-1,NY,NX)+((ARLFY-ARLFR) &
      /ARLFT(L,NY,NX))*(ZL(L,NY,NX)-ZL(L-1,NY,NX))
      ENDIF
      ELSE
      HVST(NZ,I,NY,NX)=0._r8
      ENDIF
      ARLFR=ARLFR+ARLFT(L,NY,NX)
!     WRITE(*,6544)'HVST',I,J,L,NZ,IHVST(NZ,I,NY,NX),ARLFC(NY,NX)
!    2,ARLFT(L,NY,NX),ARLFY,ARLFR,ZL(L,NY,NX),ZL(L-1,NY,NX)
!    3,ARLFV(L,NZ,NY,NX),HVST(NZ,I,NY,NX)
!6544  FORMAT(A8,5I4,20E12.4)
9875  CONTINUE
      ENDIF
      WHVSTT=0._r8
      WHVSLF=0._r8
      WHVHSH=0._r8
      WHVEAH=0._r8
      WHVGRH=0._r8
      WHVSCP=0._r8
      WHVSTH=0._r8
      WHVRVH=0._r8
      ELSE
!
!     GRAZING REMOVAL
!
!     WTSHTA=average biomass in landscape grazing section
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WHVSTT=total phytomass grazed, removed
!     TFN3=temperature function for canopy growth
!     CCPOLP=nonstructural C concentration in canopy
!     CCPLNP=nonstructural C concentration in canopy nodules
!
      IF(WTSHTA(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSTT=HVST(NZ,I,NY,NX)*THIN(NZ,I,NY,NX)*0.45/24.0 &
      *AREA(3,NU(NY,NX),NY,NX)*WTSHT(NZ,NY,NX)/WTSHTA(NZ,NY,NX)
      ELSE
      WHVSTT=0._r8
      ENDIF
      IF(IHVST(NZ,I,NY,NX).EQ.6)THEN
      WHVSTT=WHVSTT*TFN3(NZ,NY,NX)
      ENDIF
      CCPOLX=CCPOLP(NZ,NY,NX)/(1.0+CCPOLP(NZ,NY,NX))
      CCPLNX=CCPLNP(NZ,NY,NX)/(1.0+CCPLNP(NZ,NY,NX))
!
!     LEAF,BACTERIA GRAZED,REMOVED
!
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosyst
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WTLF=PFT leaf C mass
!     WHVXXX=grazing requirement unmet by leaf
!
      WHVSLX=WHVSTT*EHVST(1,1,NZ,I,NY,NX)
      WHVSLY=AMIN1(WTLF(NZ,NY,NX),WHVSLX)
      WHVSLF=WHVSLY*(1.0-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVSLX-WHVSLY)
      WHVSSX=WHVSTT*EHVST(1,2,NZ,I,NY,NX)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     WHVXXX=grazing requirement unmet by non-foliar removal
!
      WTSHTT=WTSHE(NZ,NY,NX)+WTHSK(NZ,NY,NX)+WTEAR(NZ,NY,NX) &
      +WTGR(NZ,NY,NX)
      IF(WTSHTT.GT.ZEROP(NZ,NY,NX))THEN
      WHVSHX=WHVSSX*WTSHE(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVSHY=AMIN1(WTSHE(NZ,NY,NX),WHVSHX)
      WHVSHH=WHVSHY*(1.0-CCPOLX)
      WHVSCS=WHVSHY*CCPOLX
      WHVSNS=WHVSHY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVSHX-WHVSHY)
      WHVHSX=WHVSSX*WTHSK(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVHSY=AMIN1(WTHSK(NZ,NY,NX),WHVHSX)
      WHVHSH=WHVHSY
      WHVXXX=AMAX1(0.0,WHVHSX-WHVHSY)
      WHVEAX=WHVSSX*WTEAR(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVEAY=AMIN1(WTEAR(NZ,NY,NX),WHVEAX)
      WHVEAH=WHVEAY
      WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
      WHVGRX=WHVSSX*WTGR(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVGRY=AMIN1(WTGR(NZ,NY,NX),WHVGRX)
      WHVGRH=WHVGRY
      WHVXXX=AMAX1(0.0,WHVGRX-WHVGRY)
      ELSE
      WHVSHH=0._r8
      WHVSCS=0._r8
      WHVSNS=0._r8
      WHVHSH=0._r8
      WHVEAH=0._r8
      WHVGRH=0._r8
      WHVXXX=WHVXXX+WHVSSX
      ENDIF
      WHVSCP=WHVSCL+WHVSCS
      WHVSNP=WHVSNL+WHVSNS
      WHVSKX=WHVSTT*EHVST(1,3,NZ,I,NY,NX)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     WHVXXX=grazing requirement unmet by stalk,reserve
!
      WTSTKT=WTSTK(NZ,NY,NX)+WTRSV(NZ,NY,NX)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
      WHVSTX=WHVSKX*WTSTK(NZ,NY,NX)/WTSTKT+WHVXXX
      WHVSTY=AMIN1(WTSTK(NZ,NY,NX),WHVSTX)
      WHVSTH=WHVSTY
      WHVXXX=AMAX1(0.0,WHVSTX-WHVSTY)
      WHVRVX=WHVSKX*WTRSV(NZ,NY,NX)/WTSTKT+WHVXXX
      WHVRVY=AMIN1(WTRSV(NZ,NY,NX),WHVRVX)
      WHVRVH=WHVRVY
      WHVXXX=AMAX1(0.0,WHVRVX-WHVRVY)
      ELSE
      WHVSTH=0._r8
      WHVRVH=0._r8
      WHVXXX=AMAX1(0.0,WHVSKX)
!
!     ALLOCATE UNMET DEMAND FOR GRAZING TO LEAF,PETIOLE,HUSK
!     EAR,GRAIN
!
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WHVSH*,WHVHS,WHVEA,WHVGR,WHVSC=
!            petiole,husk,ear,grain,nonstructural C removed
!
      IF(WHVXXX.GT.0.0)THEN
      WHVSLY=AMIN1(WTLF(NZ,NY,NX)-WHVSLF-WHVSCL,WHVXXX)
      WHVSLF=WHVSLF+WHVSLY*(1.0-CCPOLX)
      WHVSCL=WHVSCL+WHVSLY*CCPOLX
      WHVSNL=WHVSNL+WHVSLY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVXXX-WHVSLY)
      IF(WTSHTT.GT.ZEROP(NZ,NY,NX))THEN
      WHVSHX=WHVXXX*WTSHE(NZ,NY,NX)/WTSHTT
      WHVSHY=AMIN1(WTSHE(NZ,NY,NX),WHVSHX)
      WHVSHH=WHVSHH+WHVSHY*(1.0-CCPOLX)
      WHVSCS=WHVSCS+WHVSHY*CCPOLX
      WHVSNS=WHVSNS+WHVSHY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVXXX-WHVSHY)
      WHVHSX=WHVXXX*WTHSK(NZ,NY,NX)/WTSHTT
      WHVHSY=AMIN1(WTHSK(NZ,NY,NX),WHVHSX)
      WHVHSH=WHVHSH+WHVHSY
      WHVXXX=AMAX1(0.0,WHVXXX-WHVHSY)
      WHVEAX=WHVXXX*WTEAR(NZ,NY,NX)/WTSHTT
      WHVEAY=AMIN1(WTEAR(NZ,NY,NX),WHVEAX)
      WHVEAH=WHVEAH+WHVEAY
      WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
      WHVGRX=WHVXXX*WTGR(NZ,NY,NX)/WTSHTT
      WHVGRY=AMIN1(WTGR(NZ,NY,NX),WHVGRX)
      WHVGRH=WHVGRH+WHVGRY
      WHVXXX=AMAX1(0.0,WHVGRX-WHVGRY)
      ENDIF
      ENDIF
      ENDIF
!
!     ALL HARVEST REMOVALS
!
!     WGLFBL=branch leaf C mass in canopy layer
!
      DO 9860 NB=1,NBR(NZ,NY,NX)
      DO  L=1,JC
      DO  K=0,25
      WGLFBL(L,NB,NZ,NY,NX)=0._r8
      enddo
      enddo
9860  CONTINUE
      DO 9870 NB=1,NBR(NZ,NY,NX)
      DO  L=1,JC
      DO  K=0,25
      WGLFBL(L,NB,NZ,NY,NX)=WGLFBL(L,NB,NZ,NY,NX) &
      +WGLFL(L,K,NB,NZ,NY,NX)
      enddo
      enddo
9870  CONTINUE
      ENDIF
!
!     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     ZL=height to bottom of each canopy layer
!     FHGT=fraction of canopy layer height not harvested
!     FHVST=fraction of canopy layer mass not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
      DO 9865 L=JC,1,-1
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(IHVST(NZ,I,NY,NX).NE.3)THEN
      IF(ZL(L,NY,NX).GT.ZL(L-1,NY,NX))THEN
      FHGT=AMAX1(0.0,AMIN1(1.0,1.0-((ZL(L,NY,NX)) &
      -HVST(NZ,I,NY,NX))/(ZL(L,NY,NX)-ZL(L-1,NY,NX))))
      ELSE
      FHGT=1.0_r8
      ENDIF
      ELSE
      FHGT=0._r8
      ENDIF
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      FHVST=AMAX1(0.0,1.0-(1.0-FHGT)*EHVST(1,1,NZ,I,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSH=1.0_r8-(1.0-FHGT)*EHVST(1,1,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ELSE
      FHVSH=FHVST
      ENDIF
      ENDIF
      ELSE
      FHVST=0._r8
      FHVSH=0._r8
      ENDIF
!
!     CUT LEAVES AT HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTLF=PFT leaf C mass
!     WGLFBL=branch leaf C mass in canopy layer
!     WHVBSL,WHVSLF=layer,total leaf C mass removed
!     WGLFL=leaf node C in canopy layer
!     FHVST=fraction of leaf node mass not harvested
!
      DO 9855 NB=1,NBR(NZ,NY,NX)
      IF((IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6) &
      .AND.WTLF(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
      WHVSBL=WHVSLF*AMAX1(0.0,WGLFBL(L,NB,NZ,NY,NX))/WTLF(NZ,NY,NX)
      ELSE
      WHVSBL=0._r8
      ENDIF
      DO 9845 K=25,0,-1
      IF((IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6) &
      .OR.WHVSBL.GT.0.0)THEN
      IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WGLFL(L,K,NB,NZ,NY,NX).GT.WHVSBL)THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFL(L,K,NB,NZ,NY,NX)-WHVSBL) &
      /WGLFL(L,K,NB,NZ,NY,NX)))
      FHVSH=FHVST
      ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
      ENDIF
      ENDIF
!
!     HARVESTED LEAF AREA, C, N, P
!
!     FHVST=fraction of leaf node mass not harvested
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     ARLFL,ARSTK=leaf,stalk node area in canopy layer
!     WTHTH1,WTHNH1,WTHPH1=harvested leaf C,N,P
!     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
!     WTHTH3,WTHNH3,WTHPH3=harvested woody C,N,P
!     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
      WHVSBL=WHVSBL-(1.0-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)
      WTHTH1=WTHTH1+(1.0-FHVSH)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(1)
      WTHNH1=WTHNH1+(1.0-FHVSH)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(1)
      WTHPH1=WTHPH1+(1.0-FHVSH)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(1)
      WTHTX1=WTHTX1+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(1)
      WTHNX1=WTHNX1+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(1)
      WTHPX1=WTHPX1+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(1)
      WTHTH3=WTHTH3+(1.0-FHVSH)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(0)
      WTHNH3=WTHNH3+(1.0-FHVSH)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(0)
      WTHPH3=WTHPH3+(1.0-FHVSH)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(0)
      WTHTX3=WTHTX3+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(0)
      WTHNX3=WTHNX3+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(0)
      WTHPX3=WTHPX3+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(0)
!
!     REMAINING LEAF C,N,P AND AREA
!
      WGLFL(L,K,NB,NZ,NY,NX)=FHVST*WGLFL(L,K,NB,NZ,NY,NX)
      WGLFLN(L,K,NB,NZ,NY,NX)=FHVST*WGLFLN(L,K,NB,NZ,NY,NX)
      WGLFLP(L,K,NB,NZ,NY,NX)=FHVST*WGLFLP(L,K,NB,NZ,NY,NX)
      ARLFL(L,K,NB,NZ,NY,NX)=FHVST*ARLFL(L,K,NB,NZ,NY,NX)
      IF(K.EQ.1)THEN
      ARSTK(L,NB,NZ,NY,NX)=FHVST*ARSTK(L,NB,NZ,NY,NX)
      ENDIF
      ENDIF
!     IF(I.EQ.262.AND.K.EQ.5)THEN
!     WRITE(*,6543)'GRAZ',I,J,NZ,NB,K,L,IHVST(NZ,I,NY,NX)
!    2,ZL(L,NY,NX),ZL(L-1,NY,NX),HVST(NZ,I,NY,NX),FHVST,FHVSH
!    5,WGLFBL(L,NB,NZ,NY,NX),WTLF(NZ,NY,NX),CPOOLP(NZ,NY,NX)
!    6,ARLFL(L,K,NB,NZ,NY,NX),WGLF(K,NB,NZ,NY,NX),ARLF(K,NB,NZ,NY,NX)
!    7,HTNODE(K,NB,NZ,NY,NX)
!    7,WTSHTA(NZ,NY,NX),WHVSBL,WHVSTT,WHVSLF,WHVSHH
!    3,WHVHSH,WHVEAH,WHVGRH,WHVSCP,WHVSTH,WHVRVH,WHVXXX
!    4,WTSHTT,WHVSSX,CCPOLX
!6543  FORMAT(A8,7I4,30E12.4)
!     ENDIF
9845  CONTINUE
9855  CONTINUE
      ARLFV(L,NZ,NY,NX)=0._r8
      WGLFV(L,NZ,NY,NX)=0._r8
      ARSTV(L,NZ,NY,NX)=ARSTV(L,NZ,NY,NX)*FHVST
9865  CONTINUE
      DO 9835 NB=1,NBR(NZ,NY,NX)
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
      WGLFGX=0._r8
      WGSHGX=0._r8
      WGLFGY=0._r8
      WGSHGY=0._r8
      DO 9825 K=0,25
      ARLFG=0._r8
      WGLFG=0._r8
      WGLFNG=0._r8
      WGLFPG=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     ARLFL,ARLFV=leaf node,total area in canopy layer
!
      DO 9815 L=1,JC
      ARLFG=ARLFG+ARLFL(L,K,NB,NZ,NY,NX)
      WGLFG=WGLFG+WGLFL(L,K,NB,NZ,NY,NX)
      WGLFNG=WGLFNG+WGLFLN(L,K,NB,NZ,NY,NX)
      WGLFPG=WGLFPG+WGLFLP(L,K,NB,NZ,NY,NX)
      ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
      WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)+WGLFL(L,K,NB,NZ,NY,NX)
9815  CONTINUE
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FHVSTK=fraction of internode layer mass not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.EHVST(1,1,NZ,I,NY,NX).GT.0.0)THEN
      FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(1.0-(1.0-AMAX1(0.0,WGLFG) &
      /WGLF(K,NB,NZ,NY,NX))*EHVST(1,2,NZ,I,NY,NX) &
      /EHVST(1,1,NZ,I,NY,NX))))
      FHVSHK(K)=FHVSTK(K)
      ELSE
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      FHVSTK(K)=1.0_r8-EHVST(1,2,NZ,I,NY,NX)
      FHVSHK(K)=FHVSTK(K)
      ELSE
      FHVSTK(K)=1.0_r8-THIN(NZ,I,NY,NX)
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSHK(K)=1.0_r8-EHVST(1,2,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ELSE
      FHVSHK(K)=FHVSTK(K)
      ENDIF
      ENDIF
      ENDIF
      ELSE
      FHVSTK(K)=0._r8
      FHVSHK(K)=0._r8
      ENDIF
!
!     ACCUMULATE REMAINING BRANCH LEAF AREA, C, N, P
!
!     WGLF=leaf node C mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     ARLFB,ARLF=branch,node leaf area
!     WSLF=leaf protein mass
!
      WGLFGY=WGLFGY+WGLF(K,NB,NZ,NY,NX)
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX) &
      -WGLF(K,NB,NZ,NY,NX)+WGLFG
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX) &
      -WGLFN(K,NB,NZ,NY,NX)+WGLFNG
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX) &
      -WGLFP(K,NB,NZ,NY,NX)+WGLFPG
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)-ARLF(K,NB,NZ,NY,NX)+ARLFG
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX) &
      *ARLFG/ARLF(K,NB,NZ,NY,NX)
      ELSE
      WSLF(K,NB,NZ,NY,NX)=0._r8
      ENDIF
      ARLF(K,NB,NZ,NY,NX)=ARLFG
      WGLF(K,NB,NZ,NY,NX)=WGLFG
      WGLFN(K,NB,NZ,NY,NX)=WGLFNG
      WGLFP(K,NB,NZ,NY,NX)=WGLFPG
      WGLFGX=WGLFGX+WGLF(K,NB,NZ,NY,NX)
9825  CONTINUE
!
!     CUT SHEATHS OR PETIOLES AND STALKS HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSHE,WTSHEB=PFT,branch petiole C mass
!     WHVSBS,WHVSHH=branch, PFT petiole C mass removed
!     HTNODE=internode length
!     HTSTKX=internode length removed
!
      HTSTKX=0._r8
      IF((IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6) &
      .AND.WTSHE(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSBS=WHVSHH*WTSHEB(NB,NZ,NY,NX)/WTSHE(NZ,NY,NX)
      ELSE
      WHVSBS=0._r8
      ENDIF
      DO 9805 K=25,0,-1
!112   FORMAT(A8,8I4,12E12.4)
      IF(HTNODE(K,NB,NZ,NY,NX).GT.0.0) &
      HTSTKX=AMAX1(HTSTKX,HTNODE(K,NB,NZ,NY,NX))
!     WRITE(*,112)'VSTG',I,J,NX,NY,NZ,NB,K,IDTHB(NB,NZ,NY,NX)
!    2,VSTG(NB,NZ,NY,NX),FHVSTK(K),HTSTKX,HTNODE(K,NB,NZ,NY,NX)
!    3,HVST(NZ,I,NY,NX)
!
!     HARVESTED SHEATH OR PETIOLE C,N,P
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WHVSBS=branch petiole C mass removed
!     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!     FHVSTK=fraction of internode layer mass not harvested
!     WTHTH2,WTHNH2,WTHPH2=harvested petiole C,N,P
!     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     HTSHE,HTNODE=petiole,internode length
!
      IF((IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6) &
      .OR.WHVSBS.GT.0.0)THEN
      IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WGSHE(K,NB,NZ,NY,NX).GT.WHVSBS)THEN
      FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(WGSHE(K,NB,NZ,NY,NX)-WHVSBS) &
      /WGSHE(K,NB,NZ,NY,NX)))
      FHVSHK(K)=FHVSTK(K)
      ELSE
      FHVSTK(K)=0._r8
      FHVSHK(K)=0._r8
      ENDIF
      ENDIF
      WHVSBS=WHVSBS-(1.0-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
      WTHTH2=WTHTH2+(1.0-FHVSHK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(1)
      WTHNH2=WTHNH2+(1.0-FHVSHK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(1)
      WTHPH2=WTHPH2+(1.0-FHVSHK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(1)
      WTHTX2=WTHTX2+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX) &
      *FWODB(1)
      WTHNX2=WTHNX2+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX) &
      *FWODSN(1)
      WTHPX2=WTHPX2+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX) &
      *FWODSP(1)
      WTHTH3=WTHTH3+(1.0-FHVSHK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(0)
      WTHNH3=WTHNH3+(1.0-FHVSHK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(0)
      WTHPH3=WTHPH3+(1.0-FHVSHK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(0)
      WTHTX3=WTHTX3+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX) &
      *FWODB(0)
      WTHNX3=WTHNX3+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX) &
      *FWODSN(0)
      WTHPX3=WTHPX3+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX) &
      *FWODSP(0)
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     WGSHE=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     HTSHE=node petiole height
!     WSSHE=petiole protein mass
!
      WGSHGY=WGSHGY+WGSHE(K,NB,NZ,NY,NX)
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX) &
      -(1.0-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX) &
      -(1.0-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX) &
      -(1.0-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)
      WGSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHE(K,NB,NZ,NY,NX)
      WSSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WSSHE(K,NB,NZ,NY,NX)
      WGSHN(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHN(K,NB,NZ,NY,NX)
      WGSHP(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHP(K,NB,NZ,NY,NX)
      WSSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WSSHE(K,NB,NZ,NY,NX)
      IF(IHVST(NZ,I,NY,NX).LE.2 &
      .AND.HTSHE(K,NB,NZ,NY,NX).GT.0.0)THEN
      FHGT=AMAX1(0.0,AMIN1(1.0,(HTNODE(K,NB,NZ,NY,NX) &
      +HTSHE(K,NB,NZ,NY,NX)-HVST(NZ,I,NY,NX))/HTSHE(K,NB,NZ,NY,NX)))
      HTSHE(K,NB,NZ,NY,NX)=(1.0-FHGT)*HTSHE(K,NB,NZ,NY,NX)
      ELSE
      HTSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*HTSHE(K,NB,NZ,NY,NX)
      ENDIF
      WGSHGX=WGSHGX+WGSHE(K,NB,NZ,NY,NX)
!     IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
!     IF(HTNODE(K,NB,NZ,NY,NX).GT.HVST(NZ,I,NY,NX)
!    2.OR.IHVST(NZ,I,NY,NX).EQ.3)THEN
!     IF(test_aeqb(FHVSTK(K),0._r8).AND.K.GT.0)THEN
!     IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
!     VSTG(NB,NZ,NY,NX)=AMAX1(0.0,VSTG(NB,NZ,NY,NX)-1.0)
!     ELSE
!     VSTG(NB,NZ,NY,NX)=AMAX1(0.0,VSTG(NB,NZ,NY,NX)-0.04)
!     ENDIF
!     ENDIF
!     ENDIF
!     ENDIF
      ENDIF
9805  CONTINUE
!
!     CUT NON-STRUCTURAL C,N,P IN HARVESTED BRANCHES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     FHVST=fraction of leaf+petiole node mass not harvested
!     CPOOLG,ZPOOLG,PPOOLG=branch non-structural C,N,P mass after harvest
!     CPOLNG,ZPOLNG,PPOLNG=nonstructural C,N,P in bacteria after harvest
!     WTNDG,WTNDNG,WTNDPG=bacterial C,N,P mass after harvest
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
      CPOOLX=AMAX1(0.0,CPOOL(NB,NZ,NY,NX))
      ZPOOLX=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      PPOOLX=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
      CPOLNX=AMAX1(0.0,CPOLNB(NB,NZ,NY,NX))
      ZPOLNX=AMAX1(0.0,ZPOLNB(NB,NZ,NY,NX))
      PPOLNX=AMAX1(0.0,PPOLNB(NB,NZ,NY,NX))
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(WGLFGY+WGSHGY.GT.ZEROP(NZ,NY,NX))THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFGX+WGSHGX) &
      /(WGLFGY+WGSHGY)))
      CPOOLG=CPOOLX*FHVST
      ZPOOLG=ZPOOLX*FHVST
      PPOOLG=PPOOLX*FHVST
      CPOLNG=CPOLNX*FHVST
      ZPOLNG=ZPOLNX*FHVST
      PPOLNG=PPOLNX*FHVST
      WTNDG=WTNDB(NB,NZ,NY,NX)*FHVST
      WTNDNG=WTNDBN(NB,NZ,NY,NX)*FHVST
      WTNDPG=WTNDBP(NB,NZ,NY,NX)*FHVST
      ELSE
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
      ENDIF
      ELSE
      IF(WTLS(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
      WTLSBX=AMAX1(0.0,WTLSB(NB,NZ,NY,NX))
      IF(CPOOL(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSCX=AMAX1(0.0,WHVSCP)*WTLSBX/WTLS(NZ,NY,NX)
      CPOOLG=AMAX1(0.0,CPOOLX-WHVSCX)
      ZPOOLG=AMAX1(0.0,ZPOOLX-WHVSCX*ZPOOLX/CPOOL(NB,NZ,NY,NX))
      PPOOLG=AMAX1(0.0,PPOOLX-WHVSCX*PPOOLX/CPOOL(NB,NZ,NY,NX))
      ELSE
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      ENDIF
      IF(CPOLNB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSNX=AMAX1(0.0,WHVSNP)*WTLSBX/WTLS(NZ,NY,NX)
      CPOLNG=AMAX1(0.0,CPOLNX-WHVSNX)
      ZPOLNG=AMAX1(0.0,ZPOLNX-WHVSNX*ZPOLNX/CPOLNB(NB,NZ,NY,NX))
      PPOLNG=AMAX1(0.0,PPOLNX-WHVSNX*PPOLNX/CPOLNB(NB,NZ,NY,NX))
      WTNDG=WTNDB(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
      WTNDNG=WTNDBN(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
      WTNDPG=WTNDBP(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
      ELSE
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
      ENDIF
      ELSE
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
      ENDIF
      ENDIF
!
!     HARVESTED NON-STRUCTURAL C, N, P
!
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!
      WTHTH0=WTHTH0+CPOOLX-CPOOLG+CPOLNX-CPOLNG
      WTHNH0=WTHNH0+ZPOOLX-ZPOOLG+ZPOLNX-ZPOLNG
      WTHPH0=WTHPH0+PPOOLX-PPOOLG+PPOLNX-PPOLNG
      WTHTH0=WTHTH0+WTNDB(NB,NZ,NY,NX)-WTNDG
      WTHNH0=WTHNH0+WTNDBN(NB,NZ,NY,NX)-WTNDNG
      WTHPH0=WTHPH0+WTNDBP(NB,NZ,NY,NX)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
      CPOOL(NB,NZ,NY,NX)=CPOOLG
      ZPOOL(NB,NZ,NY,NX)=ZPOOLG
      PPOOL(NB,NZ,NY,NX)=PPOOLG
      CPOLNB(NB,NZ,NY,NX)=CPOLNG
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNG
      PPOLNB(NB,NZ,NY,NX)=PPOLNG
      WTNDB(NB,NZ,NY,NX)=WTNDG
      WTNDBN(NB,NZ,NY,NX)=WTNDNG
      WTNDBP(NB,NZ,NY,NX)=WTNDPG
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CPOOLG=branch non-structural C mass after harvest
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
      IF(ICTYP(NZ,NY,NX).EQ.4.AND.CPOOLX.GT.ZEROP(NZ,NY,NX))THEN
      FHVST4=CPOOLG/CPOOLX
      DO 9810 K=1,25
      WTHTH0=WTHTH0+(1.0-FHVST4)*CPOOL3(K,NB,NZ,NY,NX)
      WTHTH0=WTHTH0+(1.0-FHVST4)*CPOOL4(K,NB,NZ,NY,NX)
      WTHTH0=WTHTH0+(1.0-FHVST4)*CO2B(K,NB,NZ,NY,NX)
      WTHTH0=WTHTH0+(1.0-FHVST4)*HCOB(K,NB,NZ,NY,NX)
      CPOOL3(K,NB,NZ,NY,NX)=FHVST4*CPOOL3(K,NB,NZ,NY,NX)
      CPOOL4(K,NB,NZ,NY,NX)=FHVST4*CPOOL4(K,NB,NZ,NY,NX)
      CO2B(K,NB,NZ,NY,NX)=FHVST4*CO2B(K,NB,NZ,NY,NX)
      HCOB(K,NB,NZ,NY,NX)=FHVST4*HCOB(K,NB,NZ,NY,NX)
9810  CONTINUE
      ENDIF
!
!     CUT STALKS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTSTKX=internode length removed
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     FHGT=fraction of canopy layer height not harvested
!     FHVST=fraction of canopy layer mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WTSTK=stalk C mass
!
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(HTSTKX.GT.ZERO)THEN
      IF(IHVST(NZ,I,NY,NX).NE.3)THEN
      FHGT=AMAX1(0.0,AMIN1(1.0,HVST(NZ,I,NY,NX)/HTSTKX))
      ELSE
      FHGT=0._r8
      ENDIF
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      FHVST=AMAX1(0.0,1.0-(1.0-FHGT)*EHVST(1,3,NZ,I,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSH=1.0_r8-(1.0-FHGT)*EHVST(1,3,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ELSE
      FHVSH=FHVST
      ENDIF
      ENDIF
      ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
      ENDIF
      ELSE
      IF(WTSTK(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,1.0-WHVSTH/WTSTK(NZ,NY,NX)))
      FHVSH=FHVST
      ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
      ENDIF
      ENDIF
!
!     HARVESTED STALK C,N,P
!
!     WTHTH3,WTHNH3,WTHPH3=harvested stalk C,N,P
!     WTHTX3,WTHNX3,WTHPX3=harvested stalk C,N,P to litter
!     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in harvested stalk
!
      WTHTH3=WTHTH3+(1.0-FHVSH)*WTSTKB(NB,NZ,NY,NX)
      WTHNH3=WTHNH3+(1.0-FHVSH)*WTSTBN(NB,NZ,NY,NX)
      WTHPH3=WTHPH3+(1.0-FHVSH)*WTSTBP(NB,NZ,NY,NX)
      WTHTX3=WTHTX3+(FHVSH-FHVST)*WTSTKB(NB,NZ,NY,NX)
      WTHNX3=WTHNX3+(FHVSH-FHVST)*WTSTBN(NB,NZ,NY,NX)
      WTHPX3=WTHPX3+(FHVSH-FHVST)*WTSTBP(NB,NZ,NY,NX)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
      WTSTKB(NB,NZ,NY,NX)=FHVST*WTSTKB(NB,NZ,NY,NX)
      WTSTBN(NB,NZ,NY,NX)=FHVST*WTSTBN(NB,NZ,NY,NX)
      WTSTBP(NB,NZ,NY,NX)=FHVST*WTSTBP(NB,NZ,NY,NX)
      WVSTKB(NB,NZ,NY,NX)=FHVST*WVSTKB(NB,NZ,NY,NX)
      WTSTXB(NB,NZ,NY,NX)=FHVST*WTSTXB(NB,NZ,NY,NX)
      WTSTXN(NB,NZ,NY,NX)=FHVST*WTSTXN(NB,NZ,NY,NX)
      WTSTXP(NB,NZ,NY,NX)=FHVST*WTSTXP(NB,NZ,NY,NX)
!
!     CUT STALK NODES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTNODX,HTNODE=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTSTK=stalk C mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
      DO 9820 K=25,0,-1
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(HTNODX(K,NB,NZ,NY,NX).GT.ZERO)THEN
      IF(IHVST(NZ,I,NY,NX).NE.3)THEN
      FHGTK=AMAX1(0.0,AMIN1(1.0,(HTNODE(K,NB,NZ,NY,NX) &
      -HVST(NZ,I,NY,NX))/HTNODX(K,NB,NZ,NY,NX)))
      ELSE
      FHGTK=0._r8
      ENDIF
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      FHVSTS=AMAX1(0.0,1.0-FHGTK*EHVST(1,3,NZ,I,NY,NX))
      ELSE
      FHVSTS=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      ENDIF
      ELSE
      FHVSTS=1.0_r8
      ENDIF
      ELSE
      IF(WTSTK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTS=AMAX1(0.0,AMIN1(1.0,1.0-WHVSTH/WTSTK(NZ,NY,NX)))
      ELSE
      FHVSTS=1.0_r8
      ENDIF
      ENDIF
      WGNODE(K,NB,NZ,NY,NX)=FHVSTS*WGNODE(K,NB,NZ,NY,NX)
      WGNODN(K,NB,NZ,NY,NX)=FHVSTS*WGNODN(K,NB,NZ,NY,NX)
      WGNODP(K,NB,NZ,NY,NX)=FHVSTS*WGNODP(K,NB,NZ,NY,NX)
      IF(IHVST(NZ,I,NY,NX).LE.2.AND.test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      HTNODX(K,NB,NZ,NY,NX)=FHVSTS*HTNODX(K,NB,NZ,NY,NX)
      HTNODE(K,NB,NZ,NY,NX)=AMIN1(HTNODE(K,NB,NZ,NY,NX) &
      ,HVST(NZ,I,NY,NX))
      ENDIF
!     IF(NZ.EQ.2)THEN
!     WRITE(*,4811)'STK2',I,J,NX,NY,NZ,NB,K,IHVST(NZ,I,NY,NX)
!    2,HTNODX(K,NB,NZ,NY,NX),HTNODE(K,NB,NZ,NY,NX)
!    3,HVST(NZ,I,NY,NX),FHGTK,FHVSTS,ARLF(K,NB,NZ,NY,NX)
!    4,EHVST(1,3,NZ,I,NY,NX),THIN(NZ,I,NY,NX)
!4811  FORMAT(A8,8I4,12E12.4)
!     ENDIF
9820  CONTINUE
!
!     CUT STALK RESERVES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSTKB=C mass remaining in harvested stalk
!     WTRSV=stalk reserve C mass
!     WHVRVH=remaining stalk reserve C mass
!     FHVST=fraction of reserve mass not harvested
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(WTSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVST=FHVST
      FHVSH=FHVSH
      ELSE
      FHVST=0._r8
      FHVSH=0._r8
      ENDIF
      ELSE
      IF(WTRSV(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,1.0-WHVRVH/WTRSV(NZ,NY,NX)))
      FHVSH=FHVST
      ELSE
      FHVST=0._r8
      FHVSH=0._r8
      ENDIF
      ENDIF
!
!     HARVESTED STALK RESERVE C,N,P
!
!     WTHTH3,WTHNH3,WTHPH3=harvested stalk C,N,P
!     WTHTX3,WTHNX3,WTHPX3=harvested stalk C,N,P to litter
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!
      WTHTH3=WTHTH3+(1.0-FHVSH)*WTRSVB(NB,NZ,NY,NX)
      WTHNH3=WTHNH3+(1.0-FHVSH)*WTRSBN(NB,NZ,NY,NX)
      WTHPH3=WTHPH3+(1.0-FHVSH)*WTRSBP(NB,NZ,NY,NX)
      WTHTX3=WTHTX3+(FHVSH-FHVST)*WTRSVB(NB,NZ,NY,NX)
      WTHNX3=WTHNX3+(FHVSH-FHVST)*WTRSBN(NB,NZ,NY,NX)
      WTHPX3=WTHPX3+(FHVSH-FHVST)*WTRSBP(NB,NZ,NY,NX)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
      WTRSVB(NB,NZ,NY,NX)=FHVST*WTRSVB(NB,NZ,NY,NX)
      WTRSBN(NB,NZ,NY,NX)=FHVST*WTRSBN(NB,NZ,NY,NX)
      WTRSBP(NB,NZ,NY,NX)=FHVST*WTRSBP(NB,NZ,NY,NX)
!
!     CUT REPRODUCTIVE ORGANS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FHVSTG,FHVSTH,FHVSTE=fraction of grain,husk,ear mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(HVST(NZ,I,NY,NX).LT.HTSTKX &
      .OR.IHVST(NZ,I,NY,NX).EQ.1 &
      .OR.IHVST(NZ,I,NY,NX).EQ.3)THEN
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
      FHVSTG=1.0_r8-EHVST(1,2,NZ,I,NY,NX)
      FHVSHG=FHVSTG
      ELSE
      FHVSTG=1.0_r8-THIN(NZ,I,NY,NX)
      FHVSHG=1.0_r8-EHVST(1,2,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ENDIF
      ELSE
      FHVSTG=1.0_r8-THIN(NZ,I,NY,NX)
      FHVSHG=FHVSTG
      ENDIF
      FHVSTH=FHVSTG
      FHVSTE=FHVSTG
      FHVSHH=FHVSHG
      FHVSHE=FHVSHG
      ELSE
      IF(WTHSK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTH=AMAX1(0.0,AMIN1(1.0,1.0-WHVHSH/WTHSK(NZ,NY,NX)))
      FHVSHH=FHVSTH
      ELSE
      FHVSTH=1.0_r8
      FHVSHH=1.0_r8
      ENDIF
      IF(WTEAR(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTE=AMAX1(0.0,AMIN1(1.0,1.0-WHVEAH/WTEAR(NZ,NY,NX)))
      FHVSHE=FHVSTE
      ELSE
      FHVSTE=1.0_r8
      FHVSHE=1.0_r8
      ENDIF
      IF(WTGR(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTG=AMAX1(0.0,AMIN1(1.0,1.0-WHVGRH/WTGR(NZ,NY,NX)))
      FHVSHG=FHVSTG
      ELSE
      FHVSTG=1.0_r8
      FHVSHG=1.0_r8
      ENDIF
      ENDIF
!
!     HARVESTED REPRODUCTIVE C,N,P
!
!     WTHTH2,WTHNH2,WTHPH2=reproductive C,N,P removed
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTHTG,WTHNG,WTHPG=grain harvested
!
      WTHTH2=WTHTH2+(1.0-FHVSHH)*WTHSKB(NB,NZ,NY,NX)+(1.0-FHVSHE) &
      *WTEARB(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRB(NB,NZ,NY,NX)
      WTHNH2=WTHNH2+(1.0-FHVSHH)*WTHSBN(NB,NZ,NY,NX)+(1.0-FHVSHE) &
      *WTEABN(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRBN(NB,NZ,NY,NX)
      WTHPH2=WTHPH2+(1.0-FHVSHH)*WTHSBP(NB,NZ,NY,NX)+(1.0-FHVSHE) &
      *WTEABP(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRBP(NB,NZ,NY,NX)
      WTHTX2=WTHTX2+(FHVSHH-FHVSTH)*WTHSKB(NB,NZ,NY,NX)+(FHVSHE-FHVSTE) &
      *WTEARB(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRB(NB,NZ,NY,NX)
      WTHNX2=WTHNX2+(FHVSHH-FHVSTH)*WTHSBN(NB,NZ,NY,NX)+(FHVSHE-FHVSTE) &
      *WTEABN(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRBN(NB,NZ,NY,NX)
      WTHPX2=WTHPX2+(FHVSHH-FHVSTH)*WTHSBP(NB,NZ,NY,NX)+(FHVSHE-FHVSTE) &
      *WTEABP(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRBP(NB,NZ,NY,NX)
      WTHTG=WTHTG+(1.0-FHVSTG)*WTGRB(NB,NZ,NY,NX)
      WTHNG=WTHNG+(1.0-FHVSTG)*WTGRBN(NB,NZ,NY,NX)
      WTHPG=WTHPG+(1.0-FHVSTG)*WTGRBP(NB,NZ,NY,NX)
!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
      WTHSKB(NB,NZ,NY,NX)=FHVSTH*WTHSKB(NB,NZ,NY,NX)
      WTEARB(NB,NZ,NY,NX)=FHVSTE*WTEARB(NB,NZ,NY,NX)
      WTGRB(NB,NZ,NY,NX)=FHVSTG*WTGRB(NB,NZ,NY,NX)
      WTHSBN(NB,NZ,NY,NX)=FHVSTH*WTHSBN(NB,NZ,NY,NX)
      WTEABN(NB,NZ,NY,NX)=FHVSTE*WTEABN(NB,NZ,NY,NX)
      WTGRBN(NB,NZ,NY,NX)=FHVSTG*WTGRBN(NB,NZ,NY,NX)
      WTHSBP(NB,NZ,NY,NX)=FHVSTH*WTHSBP(NB,NZ,NY,NX)
      WTEABP(NB,NZ,NY,NX)=FHVSTE*WTEABP(NB,NZ,NY,NX)
      WTGRBP(NB,NZ,NY,NX)=FHVSTG*WTGRBP(NB,NZ,NY,NX)
      GRNXB(NB,NZ,NY,NX)=FHVSTG*GRNXB(NB,NZ,NY,NX)
      GRNOB(NB,NZ,NY,NX)=FHVSTG*GRNOB(NB,NZ,NY,NX)
      GRWTB(NB,NZ,NY,NX)=FHVSTG*GRWTB(NB,NZ,NY,NX)
!
!     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
!
!     CPOOLK=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WTLSB=leaf+petiole mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WVSTKB=stalk sapwood mass
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
      CPOOLK(NB,NZ,NY,NX)=0._r8
      DO 1325 K=1,25
      CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX) &
      +CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX) &
      +CO2B(K,NB,NZ,NY,NX)+HCOB(K,NB,NZ,NY,NX)
1325  CONTINUE
      WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX) &
      +WTSHEB(NB,NZ,NY,NX))
      WTSHTB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX) &
      +WTSHEB(NB,NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX) &
      +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX) &
      +CPOOL(NB,NZ,NY,NX)+CPOOLK(NB,NZ,NY,NX))
      WTSHTN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX) &
      +WTSHBN(NB,NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX) &
      +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX) &
      +ZPOOL(NB,NZ,NY,NX))
      WTSHTP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX) &
      +WTSHBP(NB,NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX) &
      +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX) &
      +PPOOL(NB,NZ,NY,NX))
      VOLWPX=VOLWP(NZ,NY,NX)
      WVPLT=AMAX1(0.0_r8,WTLS(NZ,NY,NX)+WVSTK(NZ,NY,NX))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
      VOLWP(NZ,NY,NX)=1.0E-06_r8*WVPLT/FDM
      VOLWOU=VOLWOU+VOLWPX-VOLWP(NZ,NY,NX)
      UVOLO(NY,NX)=UVOLO(NY,NX)+VOLWPX-VOLWP(NZ,NY,NX)
!
!     RESET PHENOLOGY, GROWTH STAGE IF STALKS ARE CUT
!
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     ZC=canopy height
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IDAY(1,=emergence date
!     GROUP=node number required for floral initiation
!     PSTGI=node number at floral initiation
!     PSTGF=node number at flowering
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     IFLGA=flag for initializing leafout
!
      IF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1) &
      .AND.(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6) &
      .AND.ZC(NZ,NY,NX).GT.HVST(NZ,I,NY,NX))THEN
      IF((IWTYP(NZ,NY,NX).NE.0.AND.VRNF(NB,NZ,NY,NX) &
      .LE.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX)) &
      .OR.(IWTYP(NZ,NY,NX).EQ.0 &
      .AND.IDAY(1,NB,NZ,NY,NX).NE.0))THEN
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0._r8
      VSTGX(NB,NZ,NY,NX)=0._r8
      TGSTGI(NB,NZ,NY,NX)=0._r8
      TGSTGF(NB,NZ,NY,NX)=0._r8
      FLG4(NB,NZ,NY,NX)=0._r8
      IDAY(1,NB,NZ,NY,NX)=I
      DO 3005 M=2,10
      IDAY(M,NB,NZ,NY,NX)=0
3005  CONTINUE
      IFLGA(NB,NZ,NY,NX)=0
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      DO 3010 NBX=1,NBR(NZ,NY,NX)
      IF(NBX.NE.NB1(NZ,NY,NX))THEN
      GROUP(NBX,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTGI(NBX,NZ,NY,NX)=PSTG(NBX,NZ,NY,NX)
      PSTGF(NBX,NZ,NY,NX)=0._r8
      VSTGX(NBX,NZ,NY,NX)=0._r8
      TGSTGI(NBX,NZ,NY,NX)=0._r8
      TGSTGF(NBX,NZ,NY,NX)=0._r8
      FLG4(NBX,NZ,NY,NX)=0._r8
      IDAY(1,NBX,NZ,NY,NX)=I
      DO 3015 M=2,10
      IDAY(M,NBX,NZ,NY,NX)=0
3015  CONTINUE
      IFLGA(NBX,NZ,NY,NX)=0
      ENDIF
3010  CONTINUE
      ENDIF
      ENDIF
      ENDIF
!
!     DEATH OF BRANCH IF KILLING HARVEST ENTERED IN 'READQ'
!
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!     WTLS=total PFT leaf+petiole C mass
!     WTSTK=total PFT stalk C mass
!     WVSTK=total PFT sapwood C mass
!     ARSTK=total PFT stalk surface area
!
      IF(JHVST(NZ,I,NY,NX).NE.0)IDTHB(NB,NZ,NY,NX)=1
      IF(PP(NZ,NY,NX).LE.0.0)IDTHB(NB,NZ,NY,NX)=1
9835  CONTINUE
      WTLS(NZ,NY,NX)=0._r8
      WTSTK(NZ,NY,NX)=0._r8
      WVSTK(NZ,NY,NX)=0._r8
      ARSTP(NZ,NY,NX)=0._r8
      DO 9840 NB=1,NBR(NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
      WTSTK(NZ,NY,NX)=WTSTK(NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)
      WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      DO 9830 L=1,JC
      ARSTP(NZ,NY,NX)=ARSTP(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
9830  CONTINUE
9840  CONTINUE
!
!     ROOT LITTERFALL FROM HARVESTING OR FIRE
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     THETW=soil water concentration
!     CORGC=SOC concentration
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     EFIRE=combustion  of N,P relative to C
!     FHVST,FHVSN,FHVSP=fraction of root layer C,N,P not removed by disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CNET=PFT net CO2 fixation
!     TNBP=total net biome productivity
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      XHVST=1.0_r8-THIN(NZ,I,NY,NX)
      DO 3985 N=1,MY(NZ,NY,NX)
      DO 3980 L=NU(NY,NX),NJ(NY,NX)
      IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      XHVST=1.0_r8-THIN(NZ,I,NY,NX)
      XHVSN=XHVST
      XHVSP=XHVST
      FFIRE=0._r8
      FFIRN=0._r8
      FFIRP=0._r8
      ELSE
      IF(THETW(L,NY,NX).GT.FVLWB.OR.CORGC(L,NY,NX).LE.FORGC &
      .OR.ITILL(I,NY,NX).NE.22)THEN
      XHVST=1.0_r8
      XHVSN=XHVST
      XHVSP=XHVST
      FFIRE=0._r8
      FFIRN=0._r8
      FFIRP=0._r8
      ELSE
      XHVST=1.0_r8-DCORP(I,NY,NX)*EHVST(1,3,NZ,I,NY,NX) &
      *AMIN1(1.0,(CORGC(L,NY,NX)-FORGC)/(0.55E+06-FORGC))
      XHVSN=XHVST
      XHVSP=XHVST
      FFIRE=EHVST(2,3,NZ,I,NY,NX)
      FFIRN=FFIRE*EFIRE(1,IHVST(NZ,I,NY,NX))
      FFIRP=FFIRE*EFIRE(2,IHVST(NZ,I,NY,NX))
      ENDIF
      ENDIF
      DO 3385 M=1,4
      FHVST=(1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
      FHVSN=(1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*ZPOOLR(N,L,NZ,NY,NX)
      FHVSP=(1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*PPOOLR(N,L,NZ,NY,NX)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRE)*FHVST
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRN)*FHVSN
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRP)*FHVSP
      VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
      VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*FFIRE*FHVST
      VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST*2.667
      VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-FFIRN*FHVSN
      VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
      VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-FFIRP*FHVSP
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
      TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*FFIRE*FHVST
      DO NR=1,NRT(NZ,NY,NX)
      FHVST=(1.0-XHVST)*CFOPC(5,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
      +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0)
      FHVSN=(1.0-XHVSN)*CFOPN(5,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
      +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0)
      FHVSP=(1.0-XHVSP)*CFOPP(5,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
      +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRE)*FHVST
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRN)*FHVSN
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRP)*FHVSP
      VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
      VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*FFIRE*FHVST
      VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST*2.667
      VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-FFIRN*FHVSN
      VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
      VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-FFIRP*FHVSP
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
      TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*FFIRE*FHVST
      FHVST=(1.0-XHVST)*CFOPC(4,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
      +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1)
      FHVSN=(1.0-XHVSN)*CFOPN(4,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
      +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1)
      FHVSP=(1.0-XHVSP)*CFOPP(4,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
      +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRE)*FHVST
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRN)*FHVSN
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRP)*FHVSP
      VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
      VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*FFIRE*FHVST
      VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST*2.667
      VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-FFIRN*FHVSN
      VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
      VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-FFIRP*FHVSP
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
      TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*FFIRE*FHVST
      enddo
3385  CONTINUE
!     WRITE(*,6161)'FIRE',I,J,NZ,L,N,M,VCO2F(NZ,NY,NX),FFIRE
!    2,FHVST,CFOPC(4,M,NZ,NY,NX),CPOOLR(N,L,NZ,NY,NX),THETW(L,NY,NX)
!    3,CORGC(L,NY,NX)
!6161  FORMAT(A8,6I4,20E12.4)
!
!     RELEASE ROOT GAS CONTENTS DURING HARVESTING
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(1.0-XHVST) &
      *(CO2A(N,L,NZ,NY,NX)+CO2P(N,L,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(1.0-XHVST) &
      *(OXYA(N,L,NZ,NY,NX)+OXYP(N,L,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(1.0-XHVST) &
      *(CH4A(N,L,NZ,NY,NX)+CH4P(N,L,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(1.0-XHVST) &
      *(Z2OA(N,L,NZ,NY,NX)+Z2OP(N,L,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(1.0-XHVST) &
      *(ZH3A(N,L,NZ,NY,NX)+ZH3P(N,L,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(1.0-XHVST) &
      *(H2GA(N,L,NZ,NY,NX)+H2GP(N,L,NZ,NY,NX))
      CO2A(N,L,NZ,NY,NX)=XHVST*CO2A(N,L,NZ,NY,NX)
      OXYA(N,L,NZ,NY,NX)=XHVST*OXYA(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=XHVST*CH4A(N,L,NZ,NY,NX)
      Z2OA(N,L,NZ,NY,NX)=XHVST*Z2OA(N,L,NZ,NY,NX)
      ZH3A(N,L,NZ,NY,NX)=XHVST*ZH3A(N,L,NZ,NY,NX)
      H2GA(N,L,NZ,NY,NX)=XHVST*H2GA(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=XHVST*CO2P(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=XHVST*OXYP(N,L,NZ,NY,NX)
      CH4P(N,L,NZ,NY,NX)=XHVST*CH4P(N,L,NZ,NY,NX)
      Z2OP(N,L,NZ,NY,NX)=XHVST*Z2OP(N,L,NZ,NY,NX)
      ZH3P(N,L,NZ,NY,NX)=XHVST*ZH3P(N,L,NZ,NY,NX)
      H2GP(N,L,NZ,NY,NX)=XHVST*H2GP(N,L,NZ,NY,NX)
!
!     REDUCE ROOT STATE VARIABLES DURING HARVESTING
!
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
      DO 3960 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*XHVST
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*XHVSN
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*XHVSN
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*XHVSP
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*XHVSP
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*XHVST
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*XHVST
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*XHVST
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*XHVST
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*XHVST
      RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*XHVST
3960  CONTINUE
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*XHVST
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*XHVSN
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*XHVSP
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*XHVST
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*XHVST
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*XHVST
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*XHVST
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*XHVST
      RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*XHVST
      RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*XHVST
      RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*XHVST
      RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*XHVST
      RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*XHVST
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)*XHVST
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)*XHVST
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)*XHVST
!
!     NODULE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
      IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
      DO 3395 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
      *(CFOPC(4,M,NZ,NY,NX)*WTNDL(L,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX))
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVSN) &
      *(CFOPN(4,M,NZ,NY,NX)*WTNDLN(L,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX))
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVSP) &
      *(CFOPP(4,M,NZ,NY,NX)*WTNDLP(L,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX))
3395  CONTINUE
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*XHVST
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*XHVSN
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*XHVSP
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*XHVST
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*XHVSN
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*XHVSP
      ENDIF
3980  CONTINUE
3985  CONTINUE
!
!     STORAGE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
      IF(ISTYP(NZ,NY,NX).NE.0)THEN
      DO 3400 M=1,4
      CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(0)
      ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(0)
      PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(0)
      CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(1)
      ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(1)
      PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
      +((1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(1)
3400  CONTINUE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*XHVST
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*XHVSN
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*XHVSP
      ENDIF
      ENDIF
      ENDIF
      end subroutine RemoveBiomByHarvest
!------------------------------------------------------------------------------------------

      subroutine RootNoduleBiomchemistry(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
      IF(INTYP(NZ,NY,NX).GE.1.AND.INTYP(NZ,NY,NX).LE.3)THEN
      DO 5400 L=NU(NY,NX),NIX(NZ,NY,NX)
      IF(WTRTD(1,L,NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
!
!     INITIAL INFECTION
!
      IF(WTNDL(L,NZ,NY,NX).LE.0.0)THEN
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX) &
      +WTNDI*AREA(3,NU(NY,NX),NY,NX)
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX) &
      +WTNDI*AREA(3,NU(NY,NX),NY,NX)*CNND(NZ,NY,NX)
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX) &
      +WTNDI*AREA(3,NU(NY,NX),NY,NX)*CPND(NZ,NY,NX)
      ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
      IF(WTNDL(L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CCPOLN=AMAX1(0.0,CPOOLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX))
      CZPOLN=AMAX1(0.0,ZPOOLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX))
      CPPOLN=AMAX1(0.0,PPOOLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX))
      ELSE
      CCPOLN=1.0_r8
      CZPOLN=1.0_r8
      CPPOLN=1.0_r8
      ENDIF
      IF(CCPOLN.GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0 &
      ,CZPOLN/(CZPOLN+CCPOLN*CNKI) &
      ,CPPOLN/(CPPOLN+CCPOLN*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLN/(CCPOLN+CZPOLN/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLN/(CCPOLN+CPPOLN/CPKI)))
      ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
      ENDIF
      IF(WTNDL(L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FCNPF=AMIN1(1.0 &
      ,SQRT(WTNDLN(L,NZ,NY,NX)/(WTNDL(L,NZ,NY,NX)*CNND(NZ,NY,NX))) &
      ,SQRT(WTNDLP(L,NZ,NY,NX)/(WTNDL(L,NZ,NY,NX)*CPND(NZ,NY,NX))))
      ELSE
      FCNPF=1.0_r8
      ENDIF
      SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDLM=respiration from non-structural C unltd by O2
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDL=bacterial C mass
!     TFN4=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNGR=growth function of root water potential
!
      RCNDLM=AMAX1(0.0,AMIN1(CPOOLN(L,NZ,NY,NX) &
      ,VMXO*WTNDL(L,NZ,NY,NX))*FCNPF*TFN4(L,NZ,NY,NX)*WFNGR(1,L))
      CPOOLNX=CPOOLN(L,NZ,NY,NX)
!
!     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCNDL=respiration from non-structural C ltd by O2
!     WFR=constraint by O2 consumption on all root processes
!
      RCNDL=RCNDLM*WFR(1,L,NZ,NY,NX)
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN6=temperature function for root maintenance respiration
!     WTNDLN=bacterial N mass
!
      RMNDL=AMAX1(0.0,RMPLT*TFN6(L)*WTNDLN(L,NZ,NY,NX))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDLM,RXNDL=difference between non-structural C respn and mntc respn unltd,ltd by O2
!     RGNDLM,RGNDL=growth respiration unlimited by N,P and unltd,ltd by O2
!     RSNDLM,RSNDL=excess maintenance respiration unltd,ltd by O2
!
      RXNDLM=RCNDLM-RMNDL
      RXNDL=RCNDL-RMNDL
      RGNDLM=AMAX1(0.0,RXNDLM)
      RGNDL=AMAX1(0.0,RXNDL)
      RSNDLM=AMAX1(0.0,-RXNDLM)
      RSNDL=AMAX1(0.0,-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDL,WTNDLN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNF,UPNF=layer,total root N2 fixation
!
      RGN2P=AMAX1(0.0,WTNDL(L,NZ,NY,NX)*CNND(NZ,NY,NX) &
      -WTNDLN(L,NZ,NY,NX))/EN2F
      IF(RGNDL.GT.ZEROP(NZ,NY,NX))THEN
      RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
      ELSE
      RGN2F=0._r8
      ENDIF
      RUPNF(L,NZ,NY,NX)=RGN2F*EN2F
      UPNF(NZ,NY,NX)=UPNF(NZ,NY,NX)+RUPNF(L,NZ,NY,NX)
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTRTD=root C mass
!     CCNDLR=bacteria:root ratio
!     RDNDLX=effect of CCNDLR on bacteria decomposition rate
!     CCNKR=Km for bacterial vs root mass in decomposition
!     SPNDX=specific bacterial decomposition rate at current CCNDLR
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
      RCCC=RCCZN+CCC*RCCYN
      RCCN=CNC*RCCXN
      RCCP=CPC*RCCQN
      SPNDX=SPNDL*SQRT(TFN4(L,NZ,NY,NX)*WFNGR(1,L))
      RXNDLC=SPNDX*WTNDL(L,NZ,NY,NX)
      RXNDLN=SPNDX*WTNDLN(L,NZ,NY,NX)
      RXNDLP=SPNDX*WTNDLP(L,NZ,NY,NX)
      RDNDLC=RXNDLC*(1.0-RCCC)
      RDNDLN=RXNDLN*(1.0-RCCC)*(1.0-RCCN)
      RDNDLP=RXNDLP*(1.0-RCCC)*(1.0-RCCP)
      RCNDLC=RXNDLC-RDNDLC
      RCNDLN=RXNDLN-RDNDLN
      RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
      CGNDL=AMIN1(CPOOLN(L,NZ,NY,NX)-AMIN1(RMNDL,RCNDL) &
      -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0-DMND(NZ,NY,NX)))
      GRNDG=CGNDL*DMND(NZ,NY,NX)
      RGNDG=RGN2F+CGNDL*(1.0-DMND(NZ,NY,NX))
      ZADDN=AMAX1(0.0,AMIN1(ZPOOLN(L,NZ,NY,NX) &
      ,GRNDG*CNND(NZ,NY,NX)))*CZPOLN/(CZPOLN+CZKM)
      PADDN=AMAX1(0.0,AMIN1(PPOOLN(L,NZ,NY,NX) &
      ,GRNDG*CPND(NZ,NY,NX)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
      IF(RSNDL.GT.0.0.AND.WTNDL(L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.RCCC.GT.ZERO)THEN
      RXNSNC=RSNDL/RCCC
      RXNSNN=RXNSNC*WTNDLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX)
      RXNSNP=RXNSNC*WTNDLP(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX)
      RDNSNC=RXNSNC*(1.0-RCCC)
      RDNSNN=RXNSNN*(1.0-RCCC)*(1.0-RCCN)
      RDNSNP=RXNSNP*(1.0-RCCC)*(1.0-RCCP)
      RCNSNC=RXNSNC-RDNSNC
      RCNSNN=RXNSNN-RDNSNN
      RCNSNP=RXNSNP-RDNSNP
      ELSE
      RXNSNC=0._r8
      RXNSNN=0._r8
      RXNSNP=0._r8
      RDNSNC=0._r8
      RDNSNN=0._r8
      RDNSNP=0._r8
      RCNSNC=0._r8
      RCNSNN=0._r8
      RCNSNP=0._r8
      ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     RCO2A=total root respiration
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
      RCO2TM=AMIN1(RMNDL,RCNDLM)+RGNDLM+RCNSNC
      RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
      RCO2M(1,L,NZ,NY,NX)=RCO2M(1,L,NZ,NY,NX)+RCO2TM
      RCO2N(1,L,NZ,NY,NX)=RCO2N(1,L,NZ,NY,NX)+RCO2T
      RCO2A(1,L,NZ,NY,NX)=RCO2A(1,L,NZ,NY,NX)-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
      DO 6370 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *(RDNDLC+RDNSNC)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *(RDNDLN+RDNSNN)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *(RDNDLP+RDNSNP)
6370  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNF=root N2 fixation
!
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)-AMIN1(RMNDL,RCNDL) &
      -RGN2F-CGNDL+RCNDLC
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)-ZADDN+RCNDLN+RCNSNN &
      +RUPNF(L,NZ,NY,NX)
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)+GRNDG-RXNDLC-RXNSNC
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)+ZADDN-RXNDLN-RXNSNN
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)+PADDN-RXNDLP-RXNSNP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2122)'NODGR',I,J,NZ,L,RCNDLM,RCNDL,RMNDL,RGNDL,RGN2P
!    2,RGN2F,CGNDL,RSNDL,GRNDG,ZADDN,PADDN,RCCC,RCCN,RCCP
!    8,RDNDLC,RDNDLN,RDNDLP,RCNDLC,RDNDLX,WFR(1,L,NZ,NY,NX)
!    3,WTNDL(L,NZ,NY,NX),WTNDLN(L,NZ,NY,NX),WTNDLP(L,NZ,NY,NX)
!    2,CPOOLN(L,NZ,NY,NX),ZPOOLN(L,NZ,NY,NX),PPOOLN(L,NZ,NY,NX)
!    5,FCNPF,TFN4(L,NZ,NY,NX),WFNGR(1,L),PSIRT(1,L,NZ,NY,NX)
!    5,CCPOLN,CZPOLN,CPPOLN,CPOOLNX
!    6,VMXO*WTNDL(L,NZ,NY,NX)*TFN4(L,NZ,NY,NX)*FCNPF*WFNGR(1,L)
!2122  FORMAT(A8,4I4,60E14.6)
!     ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
!     WTRTD=root C mass
!     WTNDL=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGR=parameter to calculate nonstructural C,N,P exchange
!     CCNDLR=bacteria:root ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
      IF(CPOOLR(1,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WTRTD(1,L,NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
      CCNDLR=WTNDL(L,NZ,NY,NX)/WTRTD(1,L,NZ,NY,NX)
      WTRTD1=WTRTD(1,L,NZ,NY,NX)
      WTNDL1=AMIN1(WTRTD(1,L,NZ,NY,NX) &
      ,AMAX1(WTNDI*AREA(3,NU(NY,NX),NY,NX),WTNDL(L,NZ,NY,NX)))
      WTRTDT=WTRTD1+WTNDL1
      IF(WTRTDT.GT.ZEROP(NZ,NY,NX))THEN
      FXRNX=FXRN(INTYP(NZ,NY,NX))/(1.0+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(INTYP(NZ,NY,NX))))
      CPOOLD=(CPOOLR(1,L,NZ,NY,NX)*WTNDL1 &
      -CPOOLN(L,NZ,NY,NX)*WTRTD1)/WTRTDT
      XFRC=FXRNX*CPOOLD
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLR(1,L,NZ,NY,NX)+CPOOLN(L,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOLR(1,L,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX) &
      -ZPOOLN(L,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRN=FXRNX*ZPOOLD
      PPOOLD=(PPOOLR(1,L,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX) &
      -PPOOLN(L,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRP=FXRNX*PPOOLD
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)+XFRN
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)+XFRP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2122)'NODEX',I,J,NZ,L,XFRC,XFRN,XFRP
!    3,WTRTD(1,L,NZ,NY,NX),WTNDL1,CPOOLT,CCNDLR,FXRNX
!    4,WTNDL(L,NZ,NY,NX),WTNDLN(L,NZ,NY,NX),WTNDLP(L,NZ,NY,NX)
!    2,CPOOLN(L,NZ,NY,NX),ZPOOLN(L,NZ,NY,NX),PPOOLN(L,NZ,NY,NX)
!    3,CPOOLR(1,L,NZ,NY,NX),ZPOOLR(1,L,NZ,NY,NX),PPOOLR(1,L,NZ,NY,NX)
!     ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
5400  CONTINUE
      ENDIF
      end subroutine RootNoduleBiomchemistry
!------------------------------------------------------------------------------------------

      subroutine RootBiochemistry(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution

      NIX(NZ,NY,NX)=NG(NZ,NY,NX)
      IDTHRN=0
!
!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER
!
      DO 4995 N=1,MY(NZ,NY,NX)
      DO 4990 L=NU(NY,NX),NI(NZ,NY,NX)
!
!     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
!     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
!
!     VOLX=soil layer volume excluding macropore, rocks
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!     RUONH4,RUONHB,RUON03,RUONOB=uptake from non-band,band of NH4,NO3 unlimited by O2
!     RUOH2P,RUOH2B,RUOH1P,RUOH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by O2
!     RUCNH4,RUCNHB,RUCN03,RUCNOB=uptake from non-band,band of NH4,NO3 unlimited by nonstructural C
!     RUCH2P,RUCH2B,RUCH1P,RUCH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by nonstructural C
!
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CUPRL=0.86*(RUPNH4(N,L,NZ,NY,NX)+RUPNHB(N,L,NZ,NY,NX) &
      +RUPNO3(N,L,NZ,NY,NX)+RUPNOB(N,L,NZ,NY,NX)+RUPH2P(N,L,NZ,NY,NX) &
      +RUPH2B(N,L,NZ,NY,NX)+RUPH1P(N,L,NZ,NY,NX)+RUPH1B(N,L,NZ,NY,NX))
      CUPRO=0.86*(RUONH4(N,L,NZ,NY,NX)+RUONHB(N,L,NZ,NY,NX) &
      +RUONO3(N,L,NZ,NY,NX)+RUONOB(N,L,NZ,NY,NX)+RUOH2P(N,L,NZ,NY,NX) &
      +RUOH2B(N,L,NZ,NY,NX)+RUOH1P(N,L,NZ,NY,NX)+RUOH1B(N,L,NZ,NY,NX))
      CUPRC=0.86*(RUCNH4(N,L,NZ,NY,NX)+RUCNHB(N,L,NZ,NY,NX) &
      +RUCNO3(N,L,NZ,NY,NX)+RUCNOB(N,L,NZ,NY,NX)+RUCH2P(N,L,NZ,NY,NX) &
      +RUCH2B(N,L,NZ,NY,NX)+RUCH1P(N,L,NZ,NY,NX)+RUCH1B(N,L,NZ,NY,NX))
!
!     ACCUMULATE RESPIRATION IN FLUX ARRAYS
!
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)+CUPRO
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)+CUPRC
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)-CUPRL
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-CUPRL
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!
      DO 195 K=0,4
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+RDFOMC(N,K,L,NZ,NY,NX)
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)+RDFOMN(N,K,L,NZ,NY,NX)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)+RDFOMP(N,K,L,NZ,NY,NX)
195   CONTINUE
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX) &
      +RUPNH4(N,L,NZ,NY,NX)+RUPNHB(N,L,NZ,NY,NX) &
      +RUPNO3(N,L,NZ,NY,NX)+RUPNOB(N,L,NZ,NY,NX)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX) &
      +RUPH2P(N,L,NZ,NY,NX)+RUPH2B(N,L,NZ,NY,NX) &
      +RUPH1P(N,L,NZ,NY,NX)+RUPH1B(N,L,NZ,NY,NX)
!     IF(L.EQ.1)THEN
!     WRITE(*,9881)'CUPNH4',I,J,NZ,L,N,CPOOLR(N,L,NZ,NY,NX)
!    2,ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX),CUPRL
!    2,RDFOMC(N,L,NZ,NY,NX),RDFOMN(N,L,NZ,NY,NX),RDFOMP(N,L,NZ,NY,NX)
!    2,RUPNH4(N,L,NZ,NY,NX),RUPNHB(N,L,NZ,NY,NX),RUPNO3(N,L,NZ,NY,NX)
!    2,RUPNOB(N,L,NZ,NY,NX),RUPH2P(N,L,NZ,NY,NX),RUPH2B(N,L,NZ,NY,NX)
!    3,RUPH12P(N,L,NZ,NY,NX),RUPH1B(N,L,NZ,NY,NX),WFR(N,L,NZ,NY,NX)
!9881  FORMAT(A8,5I4,30E24.16)
!     ENDIF
!
!     GROWTH OF EACH ROOT AXIS
!
      DO 4985 NR=1,NRT(NZ,NY,NX)
!
!     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
!
!     RTDP1=primary root depth from soil surface
!     RTDPP=primary root depth from canopy
!     CDPTHZ=depth from soil surface to layer bottom
!     HTSTZ=canopy height for water uptake
!     RTSK=relative primary root sink strength
!     RTSK1=primary root sink strength
!     XRTN1=number of primary root axes
!     RRAD1,RRAD2=primary, secondary root radius
!     RTNT,RLNT=total root sink strength
!
      IF(N.EQ.1)THEN
      IF(RTDP1(N,NR,NZ,NY,NX).GT.CDPTHZ(L-1,NY,NX))THEN
      IF(RTDP1(N,NR,NZ,NY,NX).LE.CDPTHZ(L,NY,NX))THEN
      RTDPP=RTDP1(N,NR,NZ,NY,NX)+HTSTZ(NZ,NY,NX)
      RTSK1(N,L,NR)=RTSK(IGTYP(NZ,NY,NX))*XRTN1 &
      *RRAD1(N,L,NZ,NY,NX)**2/RTDPP
      RTNT(N)=RTNT(N)+RTSK1(N,L,NR)
      RLNT(N,L)=RLNT(N,L)+RTSK1(N,L,NR)
      ENDIF
      ENDIF
      ENDIF
!
!     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
!     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
!
!     RTDPL=depth of primary root axis in layer
!     RTDP1=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     RTDPX=distance behind growing point for secondary roots
!     DLYR=layer thickness
!     SDPTH=seeding depth
!     HTCTL=hypocotyledon height
!     HTSTZ=canopy height for water uptake
!     RTDPS=secondary root depth from canopy
!     RTSKP,RTSKS=primary,secondary root sink strength
!     RTN2=number of secondary root axes
!     RTSK2=total secondary root sink strength
!     RTLGA=average secondary root length
!     RTNT,RLNT=total root sink strength
!
      IF(N.EQ.1)THEN
      RTDPL(NR,L)=AMAX1(0.0,RTDP1(1,NR,NZ,NY,NX)-CDPTHZ(L-1,NY,NX) &
      -RTDPX)
      RTDPL(NR,L)=AMAX1(0.0,AMIN1(DLYR(3,L,NY,NX),RTDPL(NR,L)) &
      -AMAX1(0.0,SDPTH(NZ,NY,NX)-CDPTHZ(L-1,NY,NX)-HTCTL(NZ,NY,NX)))
      RTDPS=AMAX1(SDPTH(NZ,NY,NX),CDPTHZ(L-1,NY,NX)) &
      +0.5*RTDPL(NR,L)+HTSTZ(NZ,NY,NX)
      IF(RTDPS.GT.ZERO)THEN
      RTSKP=XRTN1*RRAD1(N,L,NZ,NY,NX)**2/RTDPS
      RTSKS=RTN2(N,L,NR,NZ,NY,NX)*RRAD2(N,L,NZ,NY,NX)**2 &
      /RTLGA(N,L,NZ,NY,NX)
      IF(RTSKP+RTSKS.GT.ZEROP(NZ,NY,NX))THEN
      RTSK2(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
      ELSE
      RTSK2(N,L,NR)=0._r8
      ENDIF
      ELSE
      RTSK2(N,L,NR)=0._r8
      ENDIF
      ELSE
      RTSK2(N,L,NR)=RTN2(N,L,NR,NZ,NY,NX)*RRAD2(N,L,NZ,NY,NX)**2 &
      /RTLGA(N,L,NZ,NY,NX)
      ENDIF
      RTNT(N)=RTNT(N)+RTSK2(N,L,NR)
      RLNT(N,L)=RLNT(N,L)+RTSK2(N,L,NR)
!     IF(IYRC.EQ.2000.AND.I.LE.160)THEN
!     WRITE(*,3341)'SINK',I,J,NX,NY,NZ,L,NR,N,RTDP1(N,NR,NZ,NY,NX)
!    2,HTCTL(NZ,NY,NX),RTSK1(N,L,NR),RTSK2(N,L,NR),RLNT(N,L),RTNT(N)
!    3,XRTN1,PP(NZ,NY,NX),RRAD1(N,L,NZ,NY,NX),RTDPS,RTDPP
!    4,RTDPL(NR,L),RTN2(N,L,NR,NZ,NY,NX),RRAD2(N,L,NZ,NY,NX)
!    2,RTLGA(N,L,NZ,NY,NX),CDPTHZ(L-1,NY,NX),CDPTHZ(L,NY,NX)
!3341  FORMAT(A8,8I4,30E12.4)
!     ENDIF
4985  CONTINUE
      ENDIF
4990  CONTINUE
4995  CONTINUE
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
      DO 5010 N=1,MY(NZ,NY,NX)
      DO 5000 L=NU(NY,NX),NI(NZ,NY,NX)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VOLX=soil layer volume excluding macropore, rocks
!
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
!     WRITE(*,4994)'5004',I,J,NZ,N,L,NI(NZ,NY,NX)
!    2,NL(NY,NX),VOLX(L,NY,NX),CDPTHZ(L-1,NY,NX)
      DO 5003 LZ=L+1,NL(NY,NX)
!     WRITE(*,4994)'5003',I,J,NZ,N,L,LZ
!    2,LZ,VOLX(L,NY,NX),CDPTHZ(LZ,NY,NX)
      IF(VOLX(LZ,NY,NX).GT.ZEROS2(NY,NX) &
      .OR.LZ.EQ.NL(NY,NX))THEN
      L1=LZ
      GO TO 5004
      ENDIF
5003  CONTINUE
5004  CONTINUE
!     WRITE(*,4994)'5005',I,J,NZ,N,L,LZ
!    2,L1,VOLX(L,NY,NX),CDPTHZ(L1,NY,NX)
!4994  FORMAT(A8,7I4,12E12.4)
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     RSCS,RSCS2=soil resistance to secondary root penetration (MPa)
!     RRAD2=secondary root radius
!     WFNR=water function for root extension
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WFNGR,WFNRG=growth,respiration function of root water potential
!     PSIRT,PSIRG=root total,turgor water potential
!     DMRT=root growth yield
!
      RSCS2=RSCS(L,NY,NX)*RRAD2(N,L,NZ,NY,NX)/1.0E-03
      WFNR=AMIN1(1.0,AMAX1(0.0,PSIRG(N,L,NZ,NY,NX)-PSILM-RSCS2))
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      WFNGR(N,L)=EXP(0.05*PSIRT(N,L,NZ,NY,NX))
      WFNRG=WFNR**0.10
      ELSE
      WFNGR(N,L)=EXP(0.10*PSIRT(N,L,NZ,NY,NX))
      WFNRG=WFNR**0.25
      ENDIF
      DMRTD=1.0_r8-DMRT(NZ,NY,NX)
      RTLGL=0._r8
      RTLGZ=0._r8
      WTRTX=0._r8
      WTRTZ=0._r8
!
!     FOR EACH ROOT AXIS
!
      DO 5050 NR=1,NRT(NZ,NY,NX)
!
!     SECONDARY ROOT EXTENSION
!
      IF(L.LE.NINR(NR,NZ,NY,NX).AND.NRX(N,NR).EQ.0)THEN
!
!     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
!     TO CURRENT AXIS
!
!     RTSK2=total secondary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of secondary root sink strength in axis
!
      IF(RLNT(N,L).GT.ZEROP(NZ,NY,NX))THEN
      FRTN=RTSK2(N,L,NR)/RLNT(N,L)
      ELSE
      FRTN=1.0_r8
      ENDIF
!
!     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
      IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CNPG=AMIN1(CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX) &
      +CCPOLR(N,L,NZ,NY,NX)*CNKI),CPPOLR(N,L,NZ,NY,NX) &
      /(CPPOLR(N,L,NZ,NY,NX)+CCPOLR(N,L,NZ,NY,NX)*CPKI))
      ELSE
      CNPG=1.0_r8
      ENDIF
!
!     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT2N=secondary root N mass
!     TFN6=temperature function for root maintenance respiration
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WFNGR=growth function of root water potential
!
      RMNCR=AMAX1(0.0,RMPLT*WTRT2N(N,L,NR,NZ,NY,NX))*TFN6(L)
      IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
      RMNCR=RMNCR*WFNGR(N,L)
      ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     TFN4=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     FDBKX=termination feedback inhibition on C3 CO2
!     WFNGR=growth function of root water potential
!
      RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLR(N,L,NZ,NY,NX) &
      *TFN4(L,NZ,NY,NX))*CNPG*FDBKX(NB1(NZ,NY,NX),NZ,NY,NX) &
      *WFNGR(N,L)
!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RCO2R=RCO2RM*WFR(N,L,NZ,NY,NX)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AMAX1(0.0,RCO2XM)*WFNRG
      RCO2Y=AMAX1(0.0,RCO2X)*WFNRG
!
!     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AMAX1(0.0,ZPOOLR(N,L,NZ,NY,NX))
      PPOOLB=AMAX1(0.0,PPOOLR(N,L,NZ,NY,NX))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ,NY,NX) &
      ,PPOOLB*DMRTR/CPRTS(NZ,NY,NX))
      IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
      RCO2GM=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(N,L,NZ,NY,NX))
      ELSE
      RCO2G=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,GRTWTG=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD2M,ZADD2,PADD2=nonstructural N,P unlimited,limited by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*DMRT(NZ,NY,NX)
      GRTWTG=CGROR*DMRT(NZ,NY,NX)
      ZADD2M=AMAX1(0.0,GRTWGM*CNRTW)
      ZADD2=AMAX1(0.0,AMIN1(FRTN*ZPOOLR(N,L,NZ,NY,NX),GRTWTG*CNRTW))
      PADD2=AMAX1(0.0,AMIN1(FRTN*PPOOLR(N,L,NZ,NY,NX),GRTWTG*CPRTW))
      CNRDM=AMAX1(0.0,1.70*ZADD2M)
      CNRDA=AMAX1(0.0,1.70*ZADD2)
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAY(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0 &
      .AND.CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0 &
      ,CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX) &
      +CCPOLR(N,L,NZ,NY,NX)*CNKI) &
      ,CPPOLR(N,L,NZ,NY,NX)/(CPPOLR(N,L,NZ,NY,NX) &
      +CCPOLR(N,L,NZ,NY,NX)*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX) &
      +CZPOLR(N,L,NZ,NY,NX)/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX) &
      +CPPOLR(N,L,NZ,NY,NX)/CPKI)))
      ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
      ENDIF
      RCCC=RCCZR+CCC*RCCYR
      RCCN=CNC*RCCXR
      RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
!     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     WFR=constraint by O2 consumption on all root processes
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC2=fraction of secondary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0)THEN
      IF(-RCO2XM.LT.WTRT2(N,L,NR,NZ,NY,NX)*RCCC)THEN
      SNCRM=-RCO2XM
      ELSE
      SNCRM=AMAX1(0.0,WTRT2(N,L,NR,NZ,NY,NX)*RCCC)
      ENDIF
      ELSE
      SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0)THEN
      IF(-RCO2X.LT.WTRT2(N,L,NR,NZ,NY,NX)*RCCC)THEN
      SNCR=-RCO2X
      ELSE
      SNCR=AMAX1(0.0,WTRT2(N,L,NR,NZ,NY,NX)*RCCC) &
      *WFR(N,L,NZ,NY,NX)
      ENDIF
      ELSE
      SNCR=0._r8
      ENDIF
      IF(SNCR.GT.0.0.AND.WTRT2(N,L,NR,NZ,NY,NX) &
      .GT.ZEROP(NZ,NY,NX))THEN
      RCCR=RCCC*WTRT2(N,L,NR,NZ,NY,NX)
      RCZR=WTRT2N(N,L,NR,NZ,NY,NX)*(RCCN+(1.0-RCCN) &
      *RCCR/WTRT2(N,L,NR,NZ,NY,NX))
      RCPR=WTRT2P(N,L,NR,NZ,NY,NX)*(RCCP+(1.0-RCCP) &
      *RCCR/WTRT2(N,L,NR,NZ,NY,NX))
      IF(RCCR.GT.ZEROP(NZ,NY,NX))THEN
      FSNC2=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
      ELSE
      FSNC2=1.0_r8
      ENDIF
      ELSE
      RCCR=0._r8
      RCZR=0._r8
      RCPR=0._r8
      FSNC2=0._r8
      ENDIF
!
!     SECONDARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC2=fraction of secondary root C to be remobilized
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      DO 6350 M=1,4
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
      *FSNC2*(WTRT2(N,L,NR,NZ,NY,NX)-RCCR)*FWODR(0)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
      *FSNC2*(WTRT2N(N,L,NR,NZ,NY,NX)-RCZR)*FWODRN(0)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
      *FSNC2*(WTRT2P(N,L,NR,NZ,NY,NX)-RCPR)*FWODRP(0)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *FSNC2*(WTRT2(N,L,NR,NZ,NY,NX)-RCCR)*FWODR(1)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *FSNC2*(WTRT2N(N,L,NR,NZ,NY,NX)-RCZR)*FWODRN(1)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *FSNC2*(WTRT2P(N,L,NR,NZ,NY,NX)-RCPR)*FWODRP(1)
6350  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC2=fraction of secondary root C to be remobilized
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-AMIN1(RMNCR,RCO2R) &
      -CGROR-CNRDA-SNCR+FSNC2*RCCR
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)-ZADD2+FSNC2*RCZR
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)-PADD2+FSNC2*RCPR
!
!     TOTAL SECONDARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)+RCO2TM
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)+RCO2T
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)-RCO2T
!
!     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=secondary root length extension
!     GRTWTG=secondary root C growth ltd by O2
!     RTLG2X=specific secondary root length from startq.f
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     FSNC2=fraction of secondary root C to be remobilized
!     RTLG2=secondary root length
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      GRTLGL=GRTWTG*RTLG2X(N,NZ,NY,NX)*WFNR*FWODR(1) &
      -FSNC2*RTLG2(N,L,NR,NZ,NY,NX)
      GRTWTL=GRTWTG-FSNC2*WTRT2(N,L,NR,NZ,NY,NX)
      GRTWTN=ZADD2-FSNC2*WTRT2N(N,L,NR,NZ,NY,NX)
      GRTWTP=PADD2-FSNC2*WTRT2P(N,L,NR,NZ,NY,NX)
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     RTLG2=secondary root length
!     GRTLGL=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTFQ=root branching frequency from PFT file
!     RTN2,RTNL=number of secondary root axes
!
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)+GRTLGL
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)+GRTWTL
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)+GRTWTN
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)+GRTWTP
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX) &
      +AMIN1(CNWS(NZ,NY,NX)*WTRT2N(N,L,NR,NZ,NY,NX) &
      ,CPWS(NZ,NY,NX)*WTRT2P(N,L,NR,NZ,NY,NX))
      RTLGL=RTLGL+RTLG2(N,L,NR,NZ,NY,NX)
      WTRTX=WTRTX+WTRT2(N,L,NR,NZ,NY,NX)
      RTN2X=RTFQ(NZ,NY,NX)*XRTN1
      RTN2Y=RTFQ(NZ,NY,NX)*RTN2X
      RTN2(N,L,NR,NZ,NY,NX)=(RTN2X+RTN2Y)*DLYR(3,L,NY,NX)
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)+RTN2(N,L,NR,NZ,NY,NX)
!     IF((I/10)*10.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,9876)'RCO22',I,J,NZ,NR,L,N,NINR(NR,NZ,NY,NX)
!    2,RCO2TM,RCO2T,RMNCR,RCO2RM,RCO2R,RCO2GM,RCO2G
!    3,RCO2XM,RCO2X,CGROR,SNCRM,SNCR,CNRDA,CPOOLR(N,L,NZ,NY,NX),FRTN
!    4,TFN4(L,NZ,NY,NX),CNPG,FDBKX(NB1(NZ,NY,NX),NZ,NY,NX),WFNGR(N,L)
!    5,TFN6(L),GRTWTG,GRTWTL,GRTLGL,RTLG2(N,L,NR,NZ,NY,NX)
!     5,WTRT2(N,L,NR,NZ,NY,NX)
!    4,RCO2M(N,L,NZ,NY,NX),RCO2A(N,L,NZ,NY,NX),WFR(N,L,NZ,NY,NX)
!    8,ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
!    9,FSNC2,RLNT(N,L),RTSK1(N,L,NR),RTSK2(N,L,NR)
!    4,RTN2X,RTN2Y,XRTN1
!    5,RTDPL(NR,L),RTDNP(N,L,NZ,NY,NX)
!    5,RTDP1(1,NR,NZ,NY,NX),CDPTHZ(L-1,NY,NX),DLYR(3,L,NY,NX)
!    6,SDPTH(NZ,NY,NX),HTCTL(NZ,NY,NX)
!    5,WFNRG,RSCS2,PSILM,PSIRG(N,L,NZ,NY,NX),PSIRT(N,L,NZ,NY,NX)
!    6,FNP,RTLGP(N,L,NZ,NY,NX),ZADD2,PADD2,CUPRO,CUPRL
!    7,RUPNH4(N,L,NZ,NY,NX),RUPNHB(N,L,NZ,NY,NX)
!    8,RUPNO3(N,L,NZ,NY,NX),RUPNOB(N,L,NZ,NY,NX)
!    9,RUPH2P(N,L,NZ,NY,NX),RUPH2B(N,L,NZ,NY,NX)
!    9,RUPH1P(N,L,NZ,NY,NX),RUPH1B(N,L,NZ,NY,NX)
!    6,RDFOMN(N,L,NZ,NY,NX),RDFOMP(N,L,NZ,NY,NX)
!    2,RTN1(N,L,NZ,NY,NX),RTN2(N,L,NR,NZ,NY,NX)
!    3,RTNL(N,L,NZ,NY,NX),DLYR(3,L,NY,NX)
!9876  FORMAT(A8,7I4,100F16.8)
!     ENDIF
!
!     PRIMARY ROOT EXTENSION
!
!     BKDS=soil bulk density
!     RTDP1,RTDP1X=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     ICHKL=flag for identifying layer with primary root tip
!     RTN1=number of primary root axes
!     XRTN1=multiplier for number of primary root axes
!
      IF(N.EQ.1)THEN
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      RTDP1X=RTDP1(N,NR,NZ,NY,NX)-CDPTHZ(0,NY,NX)
      ELSE
      RTDP1X=RTDP1(N,NR,NZ,NY,NX)
      ENDIF
      IF(RTDP1X.GT.CDPTHZ(L-1,NY,NX).AND.ICHK1(N,NR).EQ.0)THEN
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)+XRTN1
      IF(RTDP1X.LE.CDPTHZ(L,NY,NX).OR.L.EQ.NJ(NY,NX))THEN
      ICHK1(N,NR)=1
!     IF(J.EQ.24.AND.NZ.EQ.2)THEN
!     WRITE(*,9874)'RTDP1',I,J,NZ,NR,L,L-1,L1,N,NINR(NR,NZ,NY,NX)
!    2,ICHK1(N,NR),RTDP1(N,NR,NZ,NY,NX),RTDP1X,RTN1(N,L,NZ,NY,NX)
!    3,CDPTHZ(L-1,NY,NX),CDPTHZ(L,NY,NX)
!9874  FORMAT(A8,10I4,12E12.4)
!     ENDIF
!
!     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
!     ATTRIBUTED TO CURRENT AXIS
!
!     RTSK1=primary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of primary root sink strength in axis
!
      IF(RLNT(N,L).GT.ZEROP(NZ,NY,NX))THEN
      FRTN=RTSK1(N,L,NR)/RLNT(N,L)
      ELSE
      FRTN=1.0_r8
      ENDIF
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     RSCS,RSCS1=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
!     WFNR=water function for root extension
!     WFNRG=respiration function of root water potential
!
      RSCS1=RSCS(L,NY,NX)*RRAD1(N,L,NZ,NY,NX)/1.0E-03
      WFNR=AMIN1(1.0,AMAX1(0.0,PSIRG(N,L,NZ,NY,NX)-PSILM-RSCS1))
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      WFNRG=WFNR**0.10
      ELSE
      WFNRG=WFNR**0.25
      ENDIF
!
!     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
      IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CNPG=AMIN1(CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX) &
      +CCPOLR(N,L,NZ,NY,NX)*CNKI),CPPOLR(N,L,NZ,NY,NX) &
      /(CPPOLR(N,L,NZ,NY,NX)+CCPOLR(N,L,NZ,NY,NX)*CPKI))
      ELSE
      CNPG=1.0_r8
      ENDIF
!
!     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT1N=primary root N mass
!     TFN6=temperature function for root maintenance respiration
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WFNGR=growth function of root water potential
!
      RMNCR=AMAX1(0.0,RMPLT*RTWT1N(N,NR,NZ,NY,NX))*TFN6(L)
      IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
      RMNCR=RMNCR*WFNGR(N,L)
      ENDIF
!
!     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of primary root sink strength in axis
!     CPOOL=non-structural C mass
!     TFN4=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     FDBKX=termination feedback inhibition on C3 CO2
!     WFNGR=growth function of root water potential
!
      RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLR(N,L,NZ,NY,NX) &
      *TFN4(L,NZ,NY,NX))*CNPG*FDBKX(NB1(NZ,NY,NX),NZ,NY,NX) &
      *WFNGR(N,L)
      IF(RTDP1X.GE.CDPTHZ(NJ(NY,NX),NY,NX))THEN
      RCO2RM=AMIN1(RMNCR,RCO2RM)
      ENDIF
!
!     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RCO2R=RCO2RM*WFR(N,L,NZ,NY,NX)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AMAX1(0.0,RCO2XM)*WFNRG
      RCO2Y=AMAX1(0.0,RCO2X)*WFNRG
!
!     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AMAX1(0.0,ZPOOLR(N,L,NZ,NY,NX))
      PPOOLB=AMAX1(0.0,PPOOLR(N,L,NZ,NY,NX))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ,NY,NX) &
      ,PPOOLB*DMRTR/CPRTS(NZ,NY,NX))
      IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
      RCO2GM=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(N,L,NZ,NY,NX))
      ELSE
      RCO2G=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,GRTWTG=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD1M,ZADD1,PADD1=nonstructural N,P unltd,ltd by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*DMRT(NZ,NY,NX)
      GRTWTG=CGROR*DMRT(NZ,NY,NX)
      ZADD1M=AMAX1(0.0,GRTWGM*CNRTW)
      ZADD1=AMAX1(0.0,AMIN1(FRTN*ZPOOLR(N,L,NZ,NY,NX),GRTWTG*CNRTW))
      PADD1=AMAX1(0.0,AMIN1(FRTN*PPOOLR(N,L,NZ,NY,NX),GRTWTG*CPRTW))
      CNRDM=AMAX1(0.0,1.70*ZADD1M)
      CNRDA=AMAX1(0.0,1.70*ZADD1)
!
!     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAY(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0 &
      .AND.CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0 &
      ,CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX) &
      +CCPOLR(N,L,NZ,NY,NX)*CNKI) &
      ,CPPOLR(N,L,NZ,NY,NX)/(CPPOLR(N,L,NZ,NY,NX) &
      +CCPOLR(N,L,NZ,NY,NX)*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX) &
      +CZPOLR(N,L,NZ,NY,NX)/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX) &
      +CPPOLR(N,L,NZ,NY,NX)/CPKI)))
      ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
      ENDIF
      RCCC=RCCZR+CCC*RCCYR
      RCCN=CNC*RCCXR
      RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
!     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     WFR=constraint by O2 consumption on all root processes
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC1=fraction of primary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0)THEN
      IF(-RCO2XM.LT.RTWT1(N,NR,NZ,NY,NX)*RCCC)THEN
      SNCRM=-RCO2XM
      ELSE
      SNCRM=AMAX1(0.0,RTWT1(N,NR,NZ,NY,NX)*RCCC)
      ENDIF
      ELSE
      SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0)THEN
      IF(-RCO2X.LT.RTWT1(N,NR,NZ,NY,NX)*RCCC)THEN
      SNCR=-RCO2X
      ELSE
      SNCR=AMAX1(0.0,RTWT1(N,NR,NZ,NY,NX)*RCCC) &
      *WFR(N,L,NZ,NY,NX)
      ENDIF
      ELSE
      SNCR=0._r8
      ENDIF
      IF(SNCR.GT.0.0.AND.RTWT1(N,NR,NZ,NY,NX) &
      .GT.ZEROP(NZ,NY,NX))THEN
      RCCR=RCCC*RTWT1(N,NR,NZ,NY,NX)
      RCZR=RTWT1N(N,NR,NZ,NY,NX)*(RCCN+(1.0-RCCN) &
      *RCCR/RTWT1(N,NR,NZ,NY,NX))
      RCPR=RTWT1P(N,NR,NZ,NY,NX)*(RCCP+(1.0-RCCP) &
      *RCCR/RTWT1(N,NR,NZ,NY,NX))
      IF(RCCR.GT.ZEROP(NZ,NY,NX))THEN
      FSNC1=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
      ELSE
      FSNC1=1.0_r8
      ENDIF
      ELSE
      RCCR=0._r8
      RCZR=0._r8
      RCPR=0._r8
      FSNC1=0._r8
      ENDIF
!
!     PRIMARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      DO 6355 M=1,4
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
      *FSNC1*(RTWT1(N,NR,NZ,NY,NX)-RCCR)*FWODR(0)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
      *FSNC1*(RTWT1N(N,NR,NZ,NY,NX)-RCZR)*FWODRN(0)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
      *FSNC1*(RTWT1P(N,NR,NZ,NY,NX)-RCPR)*FWODRP(0)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *FSNC1*(RTWT1(N,NR,NZ,NY,NX)-RCCR)*FWODR(1)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *FSNC1*(RTWT1N(N,NR,NZ,NY,NX)-RCZR)*FWODRN(1)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *FSNC1*(RTWT1P(N,NR,NZ,NY,NX)-RCPR)*FWODRP(1)
6355  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC1=fraction of primary root C to be remobilized
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-AMIN1(RMNCR,RCO2R) &
      -CGROR-CNRDA-SNCR+FSNC1*RCCR
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)-ZADD1+FSNC1*RCZR
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)-PADD1+FSNC1*RCPR
!
!     TOTAL PRIMARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
!
!     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
!     THROUGH WHICH PRIMARY ROOTS GROW
!
!     RTDP1=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     RTLG1=primary root length
!     SDPTH=seeding depth
!     FRCO2=fraction of primary root respiration attributed to layer
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!
      IF(RTDP1(N,NR,NZ,NY,NX).GT.CDPTHZ(NG(NZ,NY,NX),NY,NX))THEN
      TFRCO2=0._r8
      DO 5100 LL=NG(NZ,NY,NX),NINR(NR,NZ,NY,NX)
      IF(LL.LT.NINR(NR,NZ,NY,NX))THEN
      FRCO2=AMIN1(1.0,RTLG1(N,LL,NR,NZ,NY,NX) &
      /(RTDP1(N,NR,NZ,NY,NX)-SDPTH(NZ,NY,NX)))
      ELSE
      FRCO2=1.0_r8-TFRCO2
      ENDIF
      TFRCO2=TFRCO2+FRCO2
      RCO2M(N,LL,NZ,NY,NX)=RCO2M(N,LL,NZ,NY,NX)+RCO2TM*FRCO2
      RCO2N(N,LL,NZ,NY,NX)=RCO2N(N,LL,NZ,NY,NX)+RCO2T*FRCO2
      RCO2A(N,LL,NZ,NY,NX)=RCO2A(N,LL,NZ,NY,NX)-RCO2T*FRCO2
!     IF(NZ.EQ.2)THEN
!     WRITE(*,9877)'RCO2A',I,J,NZ,NR,L,LL,N,NG(NZ,NY,NX)
!    2,NINR(NR,NZ,NY,NX),RCO2T,FRCO2,TFRCO2,RCO2A(N,LL,NZ,NY,NX)
!    3,RTLG1(N,LL,NR,NZ,NY,NX),RTDP1(N,NR,NZ,NY,NX)
!    4,SDPTH(NZ,NY,NX)
!     ENDIF
5100  CONTINUE
      ELSE
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)+RCO2TM
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)+RCO2T
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)-RCO2T
      ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     GRTWTG=primary root C growth ltd by O2
!     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RTLG2=secondary root length
!
      GRTWTL=GRTWTG-FSNC1*RTWT1(N,NR,NZ,NY,NX)
      GRTWTN=ZADD1-FSNC1*RTWT1N(N,NR,NZ,NY,NX)
      GRTWTP=PADD1-FSNC1*RTWT1P(N,NR,NZ,NY,NX)
      IF(GRTWTL.LT.0.0)THEN
      LX=MAX(1,L-1)
      DO 5105 LL=L,LX,-1
      GRTWTM=GRTWTL
      IF(GRTWTL.LT.0.0)THEN
      IF(GRTWTL.GT.-WTRT2(N,LL,NR,NZ,NY,NX))THEN
      RTLG2(N,LL,NR,NZ,NY,NX)=RTLG2(N,LL,NR,NZ,NY,NX)+GRTWTL &
      *RTLG2(N,LL,NR,NZ,NY,NX)/WTRT2(N,LL,NR,NZ,NY,NX)
      WTRT2(N,LL,NR,NZ,NY,NX)=WTRT2(N,LL,NR,NZ,NY,NX)+GRTWTL
      GRTWTL=0._r8
      ELSE
      GRTWTL=GRTWTL+WTRT2(N,LL,NR,NZ,NY,NX)
      RTLG2(N,LL,NR,NZ,NY,NX)=0._r8
      WTRT2(N,LL,NR,NZ,NY,NX)=0._r8
      ENDIF
      ENDIF
      IF(GRTWTN.LT.0.0)THEN
      IF(GRTWTN.GT.-WTRT2N(N,LL,NR,NZ,NY,NX))THEN
      WTRT2N(N,LL,NR,NZ,NY,NX)=WTRT2N(N,LL,NR,NZ,NY,NX)+GRTWTN
      GRTWTN=0._r8
      ELSE
      GRTWTN=GRTWTN+WTRT2N(N,LL,NR,NZ,NY,NX)
      WTRT2N(N,LL,NR,NZ,NY,NX)=0._r8
      ENDIF
      ENDIF
      IF(GRTWTP.LT.0.0)THEN
      IF(GRTWTP.GT.-WTRT2P(N,LL,NR,NZ,NY,NX))THEN
      WTRT2P(N,LL,NR,NZ,NY,NX)=WTRT2P(N,LL,NR,NZ,NY,NX)+GRTWTP
      GRTWTP=0._r8
      ELSE
      GRTWTP=GRTWTP+WTRT2P(N,LL,NR,NZ,NY,NX)
      WTRT2P(N,LL,NR,NZ,NY,NX)=0._r8
      ENDIF
      ENDIF
!     WRITE(*,9876)'WTRT2',I,J,NZ,NR,LL,N
!    2,GRTWTL,GRTWTM,GRTWTG,FSNC1,SNCR,RCCR,RTWT1(N,NR,NZ,NY,NX)
!    3,WTRT2(1,LL,NR,NZ,NY,NX),WTRTL(1,LL,NZ,NY,NX)
!    3,WTRT2(2,LL,NR,NZ,NY,NX),WTRTL(2,LL,NZ,NY,NX)
!    4,RTLG2(1,LL,NR,NZ,NY,NX),RTLG1(1,LL,NR,NZ,NY,NX)
!    4,RTLG2(2,LL,NR,NZ,NY,NX),RTLG1(2,LL,NR,NZ,NY,NX)
!
!     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
!     OF SECONDARY ROOTS
!
!     GRTWTM=negative primary root C growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     FSNCM,FSNCP=fraction of mycorrhizal structural,nonstructural C to be remobilized
!     WTRTL=active root C mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     RTLG2=mycorrhizal length
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
!
      IF(GRTWTM.LT.0.0)THEN
      IF(WTRT2(1,LL,NR,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FSNCM=AMIN1(1.0,ABS(GRTWTM)/WTRT2(1,LL,NR,NZ,NY,NX))
      ELSE
      FSNCM=1.0_r8
      ENDIF
      IF(WTRTL(1,LL,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FSNCP=AMIN1(1.0,ABS(GRTWTM)/WTRTL(1,LL,NZ,NY,NX))
      ELSE
      FSNCP=1.0_r8
      ENDIF
      DO 6450 M=1,4
      CSNC(M,0,LL,NZ,NY,NX)=CSNC(M,0,LL,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
      *FSNCM*AMAX1(0.0,WTRT2(2,LL,NR,NZ,NY,NX))*FWODR(0)
      ZSNC(M,0,LL,NZ,NY,NX)=ZSNC(M,0,LL,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
      *FSNCM*AMAX1(0.0,WTRT2N(2,LL,NR,NZ,NY,NX))*FWODRN(0)
      PSNC(M,0,LL,NZ,NY,NX)=PSNC(M,0,LL,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
      *FSNCM*AMAX1(0.0,WTRT2P(2,LL,NR,NZ,NY,NX))*FWODRP(0)
      CSNC(M,1,LL,NZ,NY,NX)=CSNC(M,1,LL,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX) &
      *FSNCM*AMAX1(0.0,WTRT2(2,LL,NR,NZ,NY,NX))*FWODR(1)
      ZSNC(M,1,LL,NZ,NY,NX)=ZSNC(M,1,LL,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX) &
      *FSNCM*AMAX1(0.0,WTRT2N(2,LL,NR,NZ,NY,NX))*FWODRN(1)
      PSNC(M,1,LL,NZ,NY,NX)=PSNC(M,1,LL,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX) &
      *FSNCM*AMAX1(0.0,WTRT2P(2,LL,NR,NZ,NY,NX))*FWODRP(1)
      CSNC(M,1,LL,NZ,NY,NX)=CSNC(M,1,LL,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX) &
      *FSNCP*AMAX1(0.0,CPOOLR(2,LL,NZ,NY,NX))
      ZSNC(M,1,LL,NZ,NY,NX)=ZSNC(M,1,LL,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX) &
      *FSNCP*AMAX1(0.0,ZPOOLR(2,LL,NZ,NY,NX))
      PSNC(M,1,LL,NZ,NY,NX)=PSNC(M,1,LL,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX) &
      *FSNCP*AMAX1(0.0,PPOOLR(2,LL,NZ,NY,NX))
6450  CONTINUE
      RTLG2(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,RTLG2(2,LL,NR,NZ,NY,NX)) &
      *(1.0-FSNCM)
      WTRT2(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,WTRT2(2,LL,NR,NZ,NY,NX)) &
      *(1.0-FSNCM)
      WTRT2N(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,WTRT2N(2,LL,NR,NZ,NY,NX)) &
      *(1.0-FSNCM)
      WTRT2P(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,WTRT2P(2,LL,NR,NZ,NY,NX)) &
      *(1.0-FSNCM)
      CPOOLR(2,LL,NZ,NY,NX)=AMAX1(0.0,CPOOLR(2,LL,NZ,NY,NX)) &
      *(1.0-FSNCP)
      ZPOOLR(2,LL,NZ,NY,NX)=AMAX1(0.0,ZPOOLR(2,LL,NZ,NY,NX)) &
      *(1.0-FSNCP)
      PPOOLR(2,LL,NZ,NY,NX)=AMAX1(0.0,PPOOLR(2,LL,NZ,NY,NX)) &
      *(1.0-FSNCP)
      ENDIF
5105  CONTINUE
      ENDIF
!
!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=primary root length extension
!     GRTWTG=primary root C growth ltd by O2
!     RTLG1X=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SDPTH=seeding depth
!     FSNC1=fraction of primary root C to be remobilized
!     RTLG1=primary root length
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
      IF(GRTWTL.LT.0.0.AND.RTWT1(N,NR,NZ,NY,NX) &
      .GT.ZEROP(NZ,NY,NX))THEN
      GRTLGL=GRTWTG*RTLG1X(N,NZ,NY,NX)/PP(NZ,NY,NX)*WFNR*FWODR(1) &
      +GRTWTL*(RTDP1(N,NR,NZ,NY,NX)-SDPTH(NZ,NY,NX)) &
      /RTWT1(N,NR,NZ,NY,NX)
      ELSE
      GRTLGL=GRTWTG*RTLG1X(N,NZ,NY,NX)/PP(NZ,NY,NX)*WFNR*FWODR(1)
      ENDIF
      IF(L.LT.NJ(NY,NX))THEN
      GRTLGL=AMIN1(DLYR(3,L1,NY,NX),GRTLGL)
      ENDIF
!
!     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
!     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
!     BOUNDARY OF CURRENT LAYER
!
!     GRTLGL=primary root length extension
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!
      IF(GRTLGL.GT.ZEROP(NZ,NY,NX).AND.L.LT.NJ(NY,NX))THEN
      FGROL=AMAX1(0.0,AMIN1(1.0,(CDPTHZ(L,NY,NX) &
      -RTDP1(N,NR,NZ,NY,NX))/GRTLGL))
      IF(FGROL.LT.1.0)FGROL=0._r8
      FGROZ=AMAX1(0.0,1.0-FGROL)
      ELSE
      FGROL=1.0_r8
      FGROZ=0._r8
      ENDIF
!
!     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
!     AND AXIS NUMBER
!
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     GRTLGL=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTLG1=primary root length
!
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)+GRTWTL
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)+GRTWTN
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)+GRTWTP
      RTDP1(N,NR,NZ,NY,NX)=RTDP1(N,NR,NZ,NY,NX)+GRTLGL
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)+GRTWTL*FGROL
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)+GRTWTN*FGROL
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)+GRTWTP*FGROL
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX) &
      +AMIN1(CNWS(NZ,NY,NX)*WTRT1N(N,L,NR,NZ,NY,NX) &
      ,CPWS(NZ,NY,NX)*WTRT1P(N,L,NR,NZ,NY,NX))
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)+GRTLGL*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of GRTLGL in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     WTRTD=root C mass
!     RTLG1=primary root length
!     GRTLGL=primary root length extension
!     FRTN=fraction of primary root sink strength in axis
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     PSIRT,PSIRG,PSIRO=root total,turgor,osmotic water potential
!     NINR=deepest root layer
!
      IF(FGROZ.GT.0.0)THEN
      WTRT1(N,L1,NR,NZ,NY,NX)=WTRT1(N,L1,NR,NZ,NY,NX) &
      +GRTWTL*FGROZ
      WTRT1N(N,L1,NR,NZ,NY,NX)=WTRT1N(N,L1,NR,NZ,NY,NX) &
      +GRTWTN*FGROZ
      WTRT1P(N,L1,NR,NZ,NY,NX)=WTRT1P(N,L1,NR,NZ,NY,NX) &
      +GRTWTP*FGROZ
      WSRTL(N,L1,NZ,NY,NX)=WSRTL(N,L1,NZ,NY,NX) &
      +AMIN1(CNWS(NZ,NY,NX)*WTRT1N(N,L1,NR,NZ,NY,NX) &
      ,CPWS(NZ,NY,NX)*WTRT1P(N,L1,NR,NZ,NY,NX))
      WTRTD(N,L1,NZ,NY,NX)=WTRTD(N,L1,NZ,NY,NX) &
      +WTRT1(N,L1,NR,NZ,NY,NX)
      RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)+GRTLGL*FGROZ
      RRAD1(N,L1,NZ,NY,NX)=RRAD1(N,L,NZ,NY,NX)
      RTLGZ=RTLGZ+RTLG1(N,L1,NR,NZ,NY,NX)
      WTRTZ=WTRTZ+WTRT1(N,L1,NR,NZ,NY,NX)
      XFRC=FRTN*CPOOLR(N,L,NZ,NY,NX)
      XFRN=FRTN*ZPOOLR(N,L,NZ,NY,NX)
      XFRP=FRTN*PPOOLR(N,L,NZ,NY,NX)
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-XFRC
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)-XFRN
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)-XFRP
      CPOOLR(N,L1,NZ,NY,NX)=CPOOLR(N,L1,NZ,NY,NX)+XFRC
      ZPOOLR(N,L1,NZ,NY,NX)=ZPOOLR(N,L1,NZ,NY,NX)+XFRN
      PPOOLR(N,L1,NZ,NY,NX)=PPOOLR(N,L1,NZ,NY,NX)+XFRP
      PSIRT(N,L1,NZ,NY,NX)=PSIRT(N,L,NZ,NY,NX)
      PSIRO(N,L1,NZ,NY,NX)=PSIRO(N,L,NZ,NY,NX)
      PSIRG(N,L1,NZ,NY,NX)=PSIRG(N,L,NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=MAX(NG(NZ,NY,NX),L+1)
!     IF(NZ.EQ.2)THEN
!     WRITE(*,9877)'INFIL',I,J,NZ,NR,L,N,NINR(NR,NZ,NY,NX)
!    2,FRTN,WTRTD(N,L1,NZ,NY,NX),CPOOLR(N,L1,NZ,NY,NX)
!    2,FGROZ,RTDP1(N,NR,NZ,NY,NX),GRTLGL,CDPTHZ(L,NY,NX)
!     ENDIF
      ENDIF
!     IF(NZ.EQ.2.AND.NR.EQ.5)THEN
!     WRITE(*,9877)'RCO21',I,J,NZ,NR,L,L-1,L1,N,NINR(NR,NZ,NY,NX)
!    2,CDPTHZ(L,NY,NX),CDPTHZ(L-1,NY,NX),CDPTHZ(L1,NY,NX)
!    2,RCO2TM,RCO2T,RMNCR,RCO2RM,RCO2R,RCO2GM,RCO2G
!    3,RCO2XM,RCO2X,CGROR,SNCRM,SNCR,CNRDA,CPOOLR(N,L,NZ,NY,NX),FRTN
!    4,TFN4(L,NZ,NY,NX),CNPG,FDBKX(NB1(NZ,NY,NX),NZ,NY,NX),WFNGR(N,L)
!    5,TFN6(L),GRTWTG,GRTWTL,GRTLGL,FGROL,RTLG1(N,L,NR,NZ,NY,NX)
!    6,WTRT1(N,L,NR,NZ,NY,NX),RTDP1(N,NR,NZ,NY,NX),RTDP1X
!    3,RCO2M(N,L,NZ,NY,NX),RCO2A(N,L,NZ,NY,NX),WFR(N,L,NZ,NY,NX)
!    4,RTSK1(N,L,NR),RRAD1(N,L,NZ,NY,NX),RTDPP
!    5,PSIRG(N,L,NZ,NY,NX),WFNR,WFNRG,FWOOD(1)
!    6,FGROZ,RTWT1(N,NR,NZ,NY,NX),FSNC1
!    9,ZADD1,PADD1,ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
!    1,RUPNH4(N,L,NZ,NY,NX),RUPNO3(N,L,NZ,NY,NX)
!9877  FORMAT(A8,9I4,100E12.4)
!     ENDIF
      ENDIF
!
!     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
!     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
!     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
!     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
!
!     NINR=deepest root layer
!     VOLX=soil layer volume excluding macropore, rocks
!     RTDP1X=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!     FRTN=fraction of primary root sink strength in axis
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTLG1=primary root length
!     WSRTL=root protein C mass
!     WTRTD=root C mass
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
      IF(L.EQ.NINR(NR,NZ,NY,NX))THEN
      DO 5115 LL=L,NG(NZ,NY,NX)+1,-1
      IF(VOLX(LL-1,NY,NX).GT.ZEROS2(NY,NX) &
      .AND.(RTDP1X.LT.CDPTHZ(LL-1,NY,NX) &
      .OR.RTDP1X.LT.SDPTH(NZ,NY,NX)))THEN
      IF(RLNT(N,LL).GT.ZEROP(NZ,NY,NX))THEN
      FRTN=(RTSK1(N,LL,NR)+RTSK2(N,LL,NR))/RLNT(N,LL)
      ELSE
      FRTN=1.0_r8
      ENDIF
      DO 5110 NN=1,MY(NZ,NY,NX)
      WTRT1(NN,LL-1,NR,NZ,NY,NX)=WTRT1(NN,LL-1,NR,NZ,NY,NX) &
      +WTRT1(NN,LL,NR,NZ,NY,NX)
      WTRT1N(NN,LL-1,NR,NZ,NY,NX)=WTRT1N(NN,LL-1,NR,NZ,NY,NX) &
      +WTRT1N(NN,LL,NR,NZ,NY,NX)
      WTRT1P(NN,LL-1,NR,NZ,NY,NX)=WTRT1P(NN,LL-1,NR,NZ,NY,NX) &
      +WTRT1P(NN,LL,NR,NZ,NY,NX)
      WTRT2(NN,LL-1,NR,NZ,NY,NX)=WTRT2(NN,LL-1,NR,NZ,NY,NX) &
      +WTRT2(NN,LL,NR,NZ,NY,NX)
      WTRT2N(NN,LL-1,NR,NZ,NY,NX)=WTRT2N(NN,LL-1,NR,NZ,NY,NX) &
      +WTRT2N(NN,LL,NR,NZ,NY,NX)
      WTRT2P(NN,LL-1,NR,NZ,NY,NX)=WTRT2P(NN,LL-1,NR,NZ,NY,NX) &
      +WTRT2P(NN,LL,NR,NZ,NY,NX)
      RTLG1(NN,LL-1,NR,NZ,NY,NX)=RTLG1(NN,LL-1,NR,NZ,NY,NX) &
      +RTLG1(NN,LL,NR,NZ,NY,NX)
      WTRT1(NN,LL,NR,NZ,NY,NX)=0._r8
      WTRT1N(NN,LL,NR,NZ,NY,NX)=0._r8
      WTRT1P(NN,LL,NR,NZ,NY,NX)=0._r8
      WTRT2(NN,LL,NR,NZ,NY,NX)=0._r8
      WTRT2N(NN,LL,NR,NZ,NY,NX)=0._r8
      WTRT2P(NN,LL,NR,NZ,NY,NX)=0._r8
      RTLG1(NN,LL,NR,NZ,NY,NX)=0._r8
      XFRC=FRTN*CPOOLR(NN,LL,NZ,NY,NX)
      XFRN=FRTN*ZPOOLR(NN,LL,NZ,NY,NX)
      XFRP=FRTN*PPOOLR(NN,LL,NZ,NY,NX)
      XFRW=FRTN*WSRTL(NN,L,NZ,NY,NX)
      XFRD=FRTN*WTRTD(NN,LL,NZ,NY,NX)
      CPOOLR(NN,LL,NZ,NY,NX)=CPOOLR(NN,LL,NZ,NY,NX)-XFRC
      ZPOOLR(NN,LL,NZ,NY,NX)=ZPOOLR(NN,LL,NZ,NY,NX)-XFRN
      PPOOLR(NN,LL,NZ,NY,NX)=PPOOLR(NN,LL,NZ,NY,NX)-XFRP
      WSRTL(NN,LL,NZ,NY,NX)=WSRTL(NN,LL,NZ,NY,NX)-XFRW
      WTRTD(NN,LL,NZ,NY,NX)=WTRTD(NN,LL,NZ,NY,NX)-XFRD
      CPOOLR(NN,LL-1,NZ,NY,NX)=CPOOLR(NN,LL-1,NZ,NY,NX)+XFRC
      ZPOOLR(NN,LL-1,NZ,NY,NX)=ZPOOLR(NN,LL-1,NZ,NY,NX)+XFRN
      PPOOLR(NN,LL-1,NZ,NY,NX)=PPOOLR(NN,LL-1,NZ,NY,NX)+XFRP
      WSRTL(NN,LL-1,NZ,NY,NX)=WSRTL(NN,LL-1,NZ,NY,NX)+XFRW
      WTRTD(NN,LL-1,NZ,NY,NX)=WTRTD(NN,LL-1,NZ,NY,NX)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-FRTN*(CO2A(NN,LL,NZ,NY,NX) &
      +CO2P(NN,LL,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-FRTN*(OXYA(NN,LL,NZ,NY,NX) &
      +OXYP(NN,LL,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-FRTN*(CH4A(NN,LL,NZ,NY,NX) &
      +CH4P(NN,LL,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-FRTN*(Z2OA(NN,LL,NZ,NY,NX) &
      +Z2OP(NN,LL,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-FRTN*(ZH3A(NN,LL,NZ,NY,NX) &
      +ZH3P(NN,LL,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-FRTN*(H2GA(NN,LL,NZ,NY,NX) &
      +H2GP(NN,LL,NZ,NY,NX))
      CO2A(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CO2A(NN,LL,NZ,NY,NX)
      OXYA(NN,LL,NZ,NY,NX)=(1.0-FRTN)*OXYA(NN,LL,NZ,NY,NX)
      CH4A(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CH4A(NN,LL,NZ,NY,NX)
      Z2OA(NN,LL,NZ,NY,NX)=(1.0-FRTN)*Z2OA(NN,LL,NZ,NY,NX)
      ZH3A(NN,LL,NZ,NY,NX)=(1.0-FRTN)*ZH3A(NN,LL,NZ,NY,NX)
      H2GA(NN,LL,NZ,NY,NX)=(1.0-FRTN)*H2GA(NN,LL,NZ,NY,NX)
      CO2P(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CO2P(NN,LL,NZ,NY,NX)
      OXYP(NN,LL,NZ,NY,NX)=(1.0-FRTN)*OXYP(NN,LL,NZ,NY,NX)
      CH4P(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CH4P(NN,LL,NZ,NY,NX)
      Z2OP(NN,LL,NZ,NY,NX)=(1.0-FRTN)*Z2OP(NN,LL,NZ,NY,NX)
      ZH3P(NN,LL,NZ,NY,NX)=(1.0-FRTN)*ZH3P(NN,LL,NZ,NY,NX)
      H2GP(NN,LL,NZ,NY,NX)=(1.0-FRTN)*H2GP(NN,LL,NZ,NY,NX)
!     IF(NZ.EQ.2)THEN
!     WRITE(*,9868)'WITHDR',I,J,NZ,NR,LL,NN,NINR(NR,NZ,NY,NX)
!    2,FRTN,RTSK1(N,LL,NR),RTSK2(N,LL,NR),RLNT(N,LL)
!    2,WTRTD(NN,LL-1,NZ,NY,NX),WTRTD(NN,LL,NZ,NY,NX)
!    2,RTLG1(NN,LL-1,NR,NZ,NY,NX),RTLG1(NN,LL,NR,NZ,NY,NX)
!    2,RTLG2(NN,LL-1,NR,NZ,NY,NX),RTLG2(NN,LL,NR,NZ,NY,NX)
!    3,RTDP1(N,NR,NZ,NY,NX),RTDP1(NN,NR,NZ,NY,NX)
!    4,CPOOLR(NN,LL-1,NZ,NY,NX),CPOOLR(NN,LL,NZ,NY,NX)
!    4,WTRT1(NN,LL-1,NR,NZ,NY,NX),WTRT1(NN,LL,NR,NZ,NY,NX)
!    4,WTRT2(NN,LL-1,NR,NZ,NY,NX),WTRT2(NN,LL,NR,NZ,NY,NX)
!9868  FORMAT(A8,7I4,100E24.16)
!      ENDIF
5110  CONTINUE
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     RTN2,RTNL=number of secondary root axes
!     RTN1=number of primary root axes
!     RTLG1=primary root length
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!
      RTNL(N,LL,NZ,NY,NX)=RTNL(N,LL,NZ,NY,NX) &
      -RTN2(N,LL,NR,NZ,NY,NX)
      RTNL(N,LL-1,NZ,NY,NX)=RTNL(N,LL-1,NZ,NY,NX) &
      +RTN2(N,LL,NR,NZ,NY,NX)
      RTN2(N,LL,NR,NZ,NY,NX)=0._r8
      RTN1(N,LL,NZ,NY,NX)=RTN1(N,LL,NZ,NY,NX)-XRTN1
      IF(LL-1.GT.NG(NZ,NY,NX))THEN
      RTLG1(N,LL-1,NR,NZ,NY,NX)=DLYR(3,LL-1,NY,NX) &
      -(CDPTHZ(LL-1,NY,NX)-RTDP1(N,NR,NZ,NY,NX))
      ELSE
      RTLG1(N,LL-1,NR,NZ,NY,NX)=DLYR(3,LL-1,NY,NX) &
      -(CDPTHZ(LL-1,NY,NX)-RTDP1(N,NR,NZ,NY,NX)) &
      -(SDPTH(NZ,NY,NX)-CDPTHZ(LL-2,NY,NX))
      ENDIF
!
!     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     FRTN=fraction of primary root sink strength in axis
!     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
!
      IF(INTYP(NZ,NY,NX).GE.1.AND.INTYP(NZ,NY,NX).LE.3)THEN
      XFRC=FRTN*WTNDL(LL,NZ,NY,NX)
      XFRN=FRTN*WTNDLN(LL,NZ,NY,NX)
      XFRP=FRTN*WTNDLP(LL,NZ,NY,NX)
      WTNDL(LL,NZ,NY,NX)=WTNDL(LL,NZ,NY,NX)-XFRC
      WTNDLN(LL,NZ,NY,NX)=WTNDLN(LL,NZ,NY,NX)-XFRN
      WTNDLP(LL,NZ,NY,NX)=WTNDLP(LL,NZ,NY,NX)-XFRP
      WTNDL(LL-1,NZ,NY,NX)=WTNDL(LL-1,NZ,NY,NX)+XFRC
      WTNDLN(LL-1,NZ,NY,NX)=WTNDLN(LL-1,NZ,NY,NX)+XFRN
      WTNDLP(LL-1,NZ,NY,NX)=WTNDLP(LL-1,NZ,NY,NX)+XFRP
      XFRC=FRTN*CPOOLN(LL,NZ,NY,NX)
      XFRN=FRTN*ZPOOLN(LL,NZ,NY,NX)
      XFRP=FRTN*PPOOLN(LL,NZ,NY,NX)
      CPOOLN(LL,NZ,NY,NX)=CPOOLN(LL,NZ,NY,NX)-XFRC
      ZPOOLN(LL,NZ,NY,NX)=ZPOOLN(LL,NZ,NY,NX)-XFRN
      PPOOLN(LL,NZ,NY,NX)=PPOOLN(LL,NZ,NY,NX)-XFRP
      CPOOLN(LL-1,NZ,NY,NX)=CPOOLN(LL-1,NZ,NY,NX)+XFRC
      ZPOOLN(LL-1,NZ,NY,NX)=ZPOOLN(LL-1,NZ,NY,NX)+XFRN
      PPOOLN(LL-1,NZ,NY,NX)=PPOOLN(LL-1,NZ,NY,NX)+XFRP
!     WRITE(*,9868)'WITHDRN',I,J,NZ,NR,LL,NN,NINR(NR,NZ,NY,NX)
!    2,WTNDL(LL,NZ,NY,NX),CPOOLN(LL,NZ,NY,NX),RTDP1(N,NR,NZ,NY,NX)
      ENDIF
      NINR(NR,NZ,NY,NX)=MAX(NG(NZ,NY,NX),LL-1)
      ELSE
      GO TO 5120
      ENDIF
5115  CONTINUE
      ENDIF
5120  CONTINUE
!
!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!
      IF(WTRT1(N,L,NR,NZ,NY,NX).LT.0.0)THEN
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=0._r8
      ENDIF
      IF(WTRT2(N,L,NR,NZ,NY,NX).LT.0.0)THEN
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX)
      WTRT2(N,L,NR,NZ,NY,NX)=0._r8
      ENDIF
!
!     TOTAL PRIMARY ROOT LENGTH AND MASS
!
!     RTLGZ=total primary root length
!     WTRTZ=total primary root C mass
!     RTLG1=primary root length in soil layer
!     WTRT1=primary root C mass in soil layer
!     NINR=deepest root layer
!
      RTLGZ=RTLGZ+RTLG1(N,L,NR,NZ,NY,NX)
      WTRTZ=WTRTZ+WTRT1(N,L,NR,NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=MIN(NINR(NR,NZ,NY,NX),NJ(NY,NX))
      IF(L.EQ.NINR(NR,NZ,NY,NX))NRX(N,NR)=1
      ENDIF
      ENDIF
      RTLGZ=RTLGZ+RTLG1(N,L,NR,NZ,NY,NX)
      WTRTZ=WTRTZ+WTRT1(N,L,NR,NZ,NY,NX)
!     ENDIF
      ENDIF
      NIX(NZ,NY,NX)=MAX(NIX(NZ,NY,NX),NINR(NR,NZ,NY,NX))
5050  CONTINUE
!
!     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
!     SEASONAL STORAGE POOL IS DEPLETED
!
!     WTRTL,WTRT=total root C mass
!     WTRVC=storage C
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!     CPOOLR=non-structural C mass in root
!
      IF(L.LE.NIX(NZ,NY,NX))THEN
      IF(WTRTL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WTRT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WTRVC(NZ,NY,NX).LT.XFRX*WTRT(NZ,NY,NX))THEN
      FWTRT=WTRTL(N,L,NZ,NY,NX)/WTRT(NZ,NY,NX)
      WTRTLX=WTRTL(N,L,NZ,NY,NX)
      WTRTTX=WTRT(NZ,NY,NX)*FWTRT
      WTRTTT=WTRTLX+WTRTTX
      CPOOLX=AMAX1(0.0,CPOOLR(N,L,NZ,NY,NX))
      WTRVCX=AMAX1(0.0,WTRVC(NZ,NY,NX)*FWTRT)
      CPOOLD=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
      XFRC=AMIN1(0.0,XFRY*CPOOLD)
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-XFRC
!     WRITE(*,3471)'RVC',I,J,NX,NY,NZ,L
!    2,XFRC,CPOOLR(N,L,NZ,NY,NX),WTRTD(N,L,NZ,NY,NX)
!    3,WTRVC(NZ,NY,NX),WTRT(NZ,NY,NX),FWTRT
!3471  FORMAT(A8,6I4,12E12.4)
      ENDIF
      ENDIF
!
!     ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
!     TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
!
!     RTLGZ=total primary root length
!     WTRTZ=total primary root C mass
!     RTLGL=total secondary root length
!     WTRTX=total secondary root C mass
!     RTLGT=total root length
!     WTRTT=total root C mass
!     FWOOD=C woody fraction in root:0=woody,1=non-woody
!     PP=PFT population
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVL,RTVLW,RTVLP=root or myco total,aqueous,gaseous volume
!     RRAD1,RRAD2=primary,secondary root radius
!     RTARP=root surface area per plant
!     RTLGA=average secondary root length
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!
      IF(N.EQ.1)THEN
      RTLGZ=RTLGZ*FWODR(1)
      RTLGL=RTLGL*FWODR(1)
      ENDIF
      RTLGX=RTLGZ*PP(NZ,NY,NX)
      RTLGT=RTLGL+RTLGX
      WTRTT=WTRTX+WTRTZ
      IF(RTLGT.GT.ZEROP(NZ,NY,NX).AND.WTRTT.GT.ZEROP(NZ,NY,NX) &
      .AND.PP(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RTLGP(N,L,NZ,NY,NX)=RTLGT/PP(NZ,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
      RTDNP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)/DLYR(3,L,NY,NX)
      ELSE
      RTDNP(N,L,NZ,NY,NX)=0._r8
      ENDIF
      RTVL=AMAX1(RTAR1X(N,NZ,NY,NX)*RTLGX+RTAR2X(N,NZ,NY,NX)*RTLGL &
      ,WTRTT*DMVL(N,NZ,NY,NX)*PSIRG(N,L,NZ,NY,NX))
      RTVLP(N,L,NZ,NY,NX)=PORT(N,NZ,NY,NX)*RTVL
      RTVLW(N,L,NZ,NY,NX)=(1.0-PORT(N,NZ,NY,NX))*RTVL
      RRAD1(N,L,NZ,NY,NX)=AMAX1(RRAD1X(N,NZ,NY,NX) &
      ,(1.0+PSIRT(N,L,NZ,NY,NX)/EMODR)*RRAD1M(N,NZ,NY,NX))
      RRAD2(N,L,NZ,NY,NX)=AMAX1(RRAD2X(N,NZ,NY,NX) &
      ,(1.0+PSIRT(N,L,NZ,NY,NX)/EMODR)*RRAD2M(N,NZ,NY,NX))
      RTAR=6.283*RRAD1(N,L,NZ,NY,NX)*RTLGX &
      +6.283*RRAD2(N,L,NZ,NY,NX)*RTLGL
      IF(RTNL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RTLGA(N,L,NZ,NY,NX)=AMAX1(RTLGAX,RTLGL/RTNL(N,L,NZ,NY,NX))
      ELSE
      RTLGA(N,L,NZ,NY,NX)=RTLGAX
      ENDIF
      RTARP(N,L,NZ,NY,NX)=RTAR/PP(NZ,NY,NX)
!     IF(N.EQ.1)THEN
!     RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*RTLGAX/RTLGA(N,L,NZ,NY,NX)
!     ENDIF
!     IF((I/10)*10.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2124)'RTLGA',I,J,NZ,L,N
!    2,RTLGAX,RTLGA(N,L,NZ,NY,NX),RTLGP(N,L,NZ,NY,NX),RTLGT,PP(NZ,NY,NX)
!    3,RTLGL,CPOOLR(N,L,NZ,NY,NX),WTRTD(N,L,NZ,NY,NX)
!2124  FORMAT(A8,5I4,12E12.4)
!     ENDIF
      ELSE
      RTLGP(N,L,NZ,NY,NX)=0._r8
      RTDNP(N,L,NZ,NY,NX)=0._r8
      RTVLP(N,L,NZ,NY,NX)=0._r8
      RTVLW(N,L,NZ,NY,NX)=0._r8
      RRAD1(N,L,NZ,NY,NX)=RRAD1M(N,NZ,NY,NX)
      RRAD2(N,L,NZ,NY,NX)=RRAD2M(N,NZ,NY,NX)
      RTARP(N,L,NZ,NY,NX)=0._r8
      RTLGA(N,L,NZ,NY,NX)=RTLGAX
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(CO2A(N,L,NZ,NY,NX) &
      +CO2P(N,L,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(OXYA(N,L,NZ,NY,NX) &
      +OXYP(N,L,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(CH4A(N,L,NZ,NY,NX) &
      +CH4P(N,L,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(Z2OA(N,L,NZ,NY,NX) &
      +Z2OP(N,L,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(ZH3A(N,L,NZ,NY,NX) &
      +ZH3P(N,L,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(H2GA(N,L,NZ,NY,NX) &
      +H2GP(N,L,NZ,NY,NX))
      CO2A(N,L,NZ,NY,NX)=0._r8
      OXYA(N,L,NZ,NY,NX)=0._r8
      CH4A(N,L,NZ,NY,NX)=0._r8
      Z2OA(N,L,NZ,NY,NX)=0._r8
      ZH3A(N,L,NZ,NY,NX)=0._r8
      H2GA(N,L,NZ,NY,NX)=0._r8
      CO2P(N,L,NZ,NY,NX)=0._r8
      OXYP(N,L,NZ,NY,NX)=0._r8
      CH4P(N,L,NZ,NY,NX)=0._r8
      Z2OP(N,L,NZ,NY,NX)=0._r8
      ZH3P(N,L,NZ,NY,NX)=0._r8
      H2GP(N,L,NZ,NY,NX)=0._r8
      ENDIF
      ENDIF
5000  CONTINUE
5010  CONTINUE
      end subroutine RootBiochemistry
!------------------------------------------------------------------------------------------

      subroutine CanopyNoduleBiochemistry(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     INTYP=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
      IF(INTYP(NZ,NY,NX).GE.4)THEN
!
!     INITIAL INFECTION
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
      IF(WTNDB(NB,NZ,NY,NX).LE.0.0)THEN
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX) &
      +WTNDI*AREA(3,NU(NY,NX),NY,NX)
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX) &
      +WTNDI*AREA(3,NU(NY,NX),NY,NX)*CNND(NZ,NY,NX)
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX) &
      +WTNDI*AREA(3,NU(NY,NX),NY,NX)*CPND(NZ,NY,NX)
      ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
      IF(WTNDB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CCPOLN=AMAX1(0.0,CPOLNB(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX))
      CZPOLN=AMAX1(0.0,ZPOLNB(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX))
      CPPOLN=AMAX1(0.0,PPOLNB(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX))
      ELSE
      CCPOLN=1.0_r8
      CZPOLN=1.0_r8
      CPPOLN=1.0_r8
      ENDIF
      IF(CCPOLN.GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0 &
      ,CZPOLN/(CZPOLN+CCPOLN*CNKI) &
      ,CPPOLN/(CPPOLN+CCPOLN*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLN/(CCPOLN+CZPOLN/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0 &
      ,CCPOLN/(CCPOLN+CPPOLN/CPKI)))
      ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
      ENDIF
      IF(WTNDB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FCNPF=AMIN1(1.0 &
      ,SQRT(WTNDBN(NB,NZ,NY,NX)/(WTNDB(NB,NZ,NY,NX)*CNND(NZ,NY,NX))) &
      ,SQRT(WTNDBP(NB,NZ,NY,NX)/(WTNDB(NB,NZ,NY,NX)*CPND(NZ,NY,NX))))
      ELSE
      FCNPF=1.0_r8
      ENDIF
      SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDL=respiration from non-structural C
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDB=bacterial C mass
!     TFN3=temperature function for canopy growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNG=growth function of canopy water potential
!
      RCNDL=AMAX1(0.0,AMIN1(CPOLNB(NB,NZ,NY,NX) &
      ,VMXO*WTNDB(NB,NZ,NY,NX))*FCNPF*TFN3(NZ,NY,NX)*WFNG)
!     CPOOLNX=CPOLNB(NB,NZ,NY,NX)
!     VMXOX=VMXO*WTNDB(NB,NZ,NY,NX)*FCNPF*TFN3(NZ,NY,NX)*WFNG
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN5=temperature function for canopy maintenance respiration
!     WTNDBN=bacterial N mass
!
      RMNDL=AMAX1(0.0,RMPLT*TFN5*WTNDBN(NB,NZ,NY,NX))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDL=difference between non-structural C respn and mntc respn
!     RGNDL=growth respiration unlimited by N,P
!     RSNDL=excess maintenance respiration
!
      RXNDL=RCNDL-RMNDL
      RGNDL=AMAX1(0.0,RXNDL)
      RSNDL=AMAX1(0.0,-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDB,WTNDBN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNFB,UPNFC=branch,total N2 fixation
!
      RGN2P=AMAX1(0.0,WTNDB(NB,NZ,NY,NX)*CNND(NZ,NY,NX) &
      -WTNDBN(NB,NZ,NY,NX))/EN2F
      IF(RGNDL.GT.ZEROP(NZ,NY,NX))THEN
      RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
      ELSE
      RGN2F=0._r8
      ENDIF
      RUPNFB=RGN2F*EN2F
      UPNFC(NZ,NY,NX)=UPNFC(NZ,NY,NX)+RUPNFB
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION AND LEAKAGE
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTLSB=leaf+petiole mass
!     CCNDLB=bacteria:leaf+petiole ratio
!     RDNDBX=effect of CCNDLB on bacteria decomposition rate
!     SPNDX=specific bacterial decomposition rate at current CCNDLB
!     SPNDL=specific decomposition rate by bacterial N2 fixers
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
      RCCC=RCCZN+CCC*RCCYN
      RCCN=CNC*RCCXN
      RCCP=CPC*RCCQN
      SPNDX=SPNDL*SQRT(TFN3(NZ,NY,NX)*WFNG)
      RXNDLC=SPNDX*WTNDB(NB,NZ,NY,NX)
      RXNDLN=SPNDX*WTNDBN(NB,NZ,NY,NX)
      RXNDLP=SPNDX*WTNDBP(NB,NZ,NY,NX)
      RDNDLC=RXNDLC*(1.0-RCCC)
      RDNDLN=RXNDLN*(1.0-RCCC)*(1.0-RCCN)
      RDNDLP=RXNDLP*(1.0-RCCC)*(1.0-RCCP)
      RCNDLC=RXNDLC-RDNDLC
      RCNDLN=RXNDLN-RDNDLN
      RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
      CGNDL=AMIN1(CPOLNB(NB,NZ,NY,NX)-AMIN1(RMNDL,RCNDL) &
      -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0-DMND(NZ,NY,NX)))
      GRNDG=CGNDL*DMND(NZ,NY,NX)
      RGNDG=RGN2F+CGNDL*(1.0-DMND(NZ,NY,NX))
      ZADDN=AMAX1(0.0,AMIN1(ZPOLNB(NB,NZ,NY,NX) &
      ,GRNDG*CNND(NZ,NY,NX)))*CZPOLN/(CZPOLN+CZKM)
      PADDN=AMAX1(0.0,AMIN1(PPOLNB(NB,NZ,NY,NX) &
      ,GRNDG*CPND(NZ,NY,NX)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
      IF(RSNDL.GT.0.0.AND.WTNDB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.RCCC.GT.ZERO)THEN
      RXNSNC=RSNDL/RCCC
      RXNSNN=RXNSNC*WTNDBN(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX)
      RXNSNP=RXNSNC*WTNDBP(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX)
      RDNSNC=RXNSNC*(1.0-RCCC)
      RDNSNN=RXNSNN*(1.0-RCCC)*(1.0-RCCN)
      RDNSNP=RXNSNP*(1.0-RCCC)*(1.0-RCCP)
      RCNSNC=RXNSNC-RDNSNC
      RCNSNN=RXNSNN-RDNSNN
      RCNSNP=RXNSNP-RDNSNP
      ELSE
      RXNSNC=0._r8
      RXNSNN=0._r8
      RXNSNP=0._r8
      RDNSNC=0._r8
      RDNSNN=0._r8
      RDNSNP=0._r8
      RCNSNC=0._r8
      RCNSNN=0._r8
      RCNSNP=0._r8
      ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2T=total C respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
      RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-RCO2T
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-RCO2T
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-RCO2T
      RECO(NY,NX)=RECO(NY,NX)-RCO2T
      TRAU(NY,NX)=TRAU(NY,NX)-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
      DO 6470 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX) &
      *(RDNDLC+RDNSNC)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX) &
      *(RDNDLN+RDNSNN)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX) &
      *(RDNDLP+RDNSNP)
6470  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNFB=branch N2 fixation
!
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)-AMIN1(RMNDL,RCNDL) &
      -RGN2F-CGNDL+RCNDLC
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)-ZADDN+RCNDLN+RCNSNN &
      +RUPNFB
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)+GRNDG-RXNDLC-RXNSNC
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)+ZADDN-RXNDLN-RXNSNN
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)+PADDN-RXNDLP-RXNSNP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2121)'NODGR',I,J,NZ,NB
!    2,RCNDL,RMNDL,RGNDL,RGN2P,RCO2T,RXNDLC,SPNDLI
!    2,RGN2P,RGN2F,CGNDL,RSNDL,GRNDG,RGNDG,RCNSNC
!    3,ZADDN,PADDN,RCCC,RCCN
!    8,RCCP,RDNDLC,RDNDLN,RDNDLP,RDNDLX,WTLSB(NB,NZ,NY,NX)
!    3,WTNDB(NB,NZ,NY,NX),WTNDBN(NB,NZ,NY,NX),WTNDBP(NB,NZ,NY,NX)
!    4,CPOLNB(NB,NZ,NY,NX),ZPOLNB(NB,NZ,NY,NX),PPOLNB(NB,NZ,NY,NX)
!    5,CCPOLN,CZPOLN,CPPOLN
!    6,TFN3(NZ,NY,NX),FCNPF,WFNG,CCNDLB,RDNDBX,CPOOLNX,VMXOX
!2121  FORMAT(A8,4I4,60F16.8)
!     ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     WTLSB=leaf+petiole C mass
!     WTNDB=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGB=parameter to calculate nonstructural C,N,P exchange
!     CCNDLB=bacteria:leaf+petiole ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!
      IF(CPOOL(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WTLSB(NB,NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
      CCNDLB=WTNDB(NB,NZ,NY,NX)/WTLSB(NB,NZ,NY,NX)
      WTLSB1=WTLSB(NB,NZ,NY,NX)
      WTNDB1=AMIN1(WTLSB(NB,NZ,NY,NX) &
      ,AMAX1(WTNDI*AREA(3,NU(NY,NX),NY,NX),WTNDB(NB,NZ,NY,NX)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROP(NZ,NY,NX))THEN
      FXRNX=FXRN(INTYP(NZ,NY,NX))/(1.0+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(INTYP(NZ,NY,NX))))
      CPOOLD=(CPOOL(NB,NZ,NY,NX)*WTNDB1 &
      -CPOLNB(NB,NZ,NY,NX)*WTLSB1)/WTLSBT
      XFRC=FXRNX*CPOOLD
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)+XFRC
      CPOOLT=CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOL(NB,NZ,NY,NX)*CPOLNB(NB,NZ,NY,NX) &
      -ZPOLNB(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      XFRN=FXRNX*ZPOOLD
      PPOOLD=(PPOOL(NB,NZ,NY,NX)*CPOLNB(NB,NZ,NY,NX) &
      -PPOLNB(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      XFRP=FXRNX*PPOOLD
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)+XFRN
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)+XFRP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2120)'NODEX',I,J,NZ,NB,IFLGA(NB,NZ,NY,NX)
!    2,XFRC,XFRN,XFRP
!    3,WTLSB(NB,NZ,NY,NX),WTNDB(NB,NZ,NY,NX),CPOOLT,CCNDLB,FXRNX
!    4,CPOLNB(NB,NZ,NY,NX),ZPOLNB(NB,NZ,NY,NX),PPOLNB(NB,NZ,NY,NX)
!    4,CPOOL(NB,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX),PPOOL(NB,NZ,NY,NX)
!    5,WTLSB1,WTNDB1,WTLSBT
!2120  FORMAT(A8,5I4,40E12.4)
!     ENDIF
!     WRITE(*,2121)'NODBAL',I,J,NZ,NB,CPOLNB(NB,NZ,NY,NX)
!    2,WTNDB(NB,NZ,NY,NX),CPOLNB(NB,NZ,NY,NX)+WTNDB(NB,NZ,NY,NX)
!    3,RMNDL,RCNDL,RGN2F,CGNDL,RCNDLC,GRNDG,RXNDLC,RXNSNC,RCO2T
!    4,RGNDG,RGNDL,RCNSNC
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      end subroutine CanopyNoduleBiochemistry
!------------------------------------------------------------------------------------------

      subroutine C4PhotoProductTransfer(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution

      DO 170 K=1,25
      IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!     MESOPHYLL TO BUNDLE SHEATH TRANSFER
!
!     WGLF=node leaf C mass
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CH2O3,CH2O4=total CO2 fixation in bundle sheath,mesophyll
!     CPL4M=mesophyll to bundle sheath transfer of nonstructural C4
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)-CH2O3(K)
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)+CH2O4(K)
      CPL4M=1.0_r8*(CPOOL4(K,NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)*FBS &
      -CPOOL3(K,NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)*FMP) &
      /(WGLF(K,NB,NZ,NY,NX)*(FBS+FMP))
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)-CPL4M
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)+CPL4M
!
!     BUNDLE SHEATH CO2 DECARBOXYLATION
!
!     CCBS=CO2 concn in bundle sheath (uM)
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WGLF=node leaf C mass
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!     CPL3K=bundle sheath CO2 decarboxylation
!     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll in C4 (uM)
!     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
!     CPOOL3=C4 nonstructural C mass in bundle sheath
!
      CCBS=AMAX1(0.0,0.083E+09*CO2B(K,NB,NZ,NY,NX) &
      /(WGLF(K,NB,NZ,NY,NX)*FBS))
      CPL3K=2.5E-02*CPOOL3(K,NB,NZ,NY,NX)/(1.0+CCBS/CO2KI)
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)-CPL3K
      CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)+FCO2B*CPL3K
      HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)+FHCOB*CPL3K
!
!     BUNDLE SHEATH LEAKAGE
!
!     CO2LK=bundle sheath CO2 leakage
!     CCBS=CO2 concn in bundle sheath (uM)
!     CO2L=intercellular CO2 concentration (uM)
!     WGLF=node leaf C mass
!     FBS=leaf water content in bundle sheath
!     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
      CO2LK=5.0E-07*(CCBS-CO2L(NZ,NY,NX))*WGLF(K,NB,NZ,NY,NX)*FBS
      CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)-FCO2B*CO2LK
      HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)-FHCOB*CO2LK
!     IF(NB.EQ.1.AND.K.EQ.14)THEN
!     WRITE(*,6667)'CO2K',I,J,NB,K,CO2LK,CO2B(K,NB,NZ,NY,NX)
!    2,HCOB(K,NB,NZ,NY,NX),CPOOL3(K,NB,NZ,NY,NX),CH2O3(K),CH2O4(K)
!    3,CCBS,CO2L(NZ,NY,NX),WGLF(K,NB,NZ,NY,NX),CPL3Z
!    4,CPL3K,CPL4M,CPOOL4(K,NB,NZ,NY,NX),FBS,FMP
!6667  FORMAT(A8,4I4,30E14.6)
!     ENDIF
!
!     TOTAL C EXCHANGE
!
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!     CO2LK=bundle sheath CO2 leakage
!
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-CO2LK
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-CO2LK
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-CO2LK
      RECO(NY,NX)=RECO(NY,NX)-CO2LK
      TRAU(NY,NX)=TRAU(NY,NX)-CO2LK
      ENDIF
170   CONTINUE
      end subroutine C4PhotoProductTransfer
!------------------------------------------------------------------------------------------

  subroutine ComputeGPP(I,J,NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NB,NZ,NY,NX

! begin_execution

! FDBK=N,P feedback inhibition on C3 CO2 fixation
! SSIN=sine of solar angle
! RADP=total PAR absorbed by canopy
! CO2Q=canopy air CO2 concentration
!
  IF(IDAY(1,NB,NZ,NY,NX).NE.0)THEN
!   IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!     WRITE(*,5651)'CHECK1',I,J,NZ,NB,IDAY(1,NB,NZ,NY,NX)
!    2,FDBK(NB,NZ,NY,NX),RADP(NZ,NY,NX),CO2Q(NZ,NY,NX)
!    3,ARLF(1,NB,NZ,NY,NX)
!5651  FORMAT(A8,5I4,12E12.4)
!   ENDIF

    IF(abs(FDBK(NB,NZ,NY,NX)).GT.0)THEN
      IF(SSIN(NY,NX).GT.0.0.AND.RADP(NZ,NY,NX).GT.0.0 &
        .AND.CO2Q(NZ,NY,NX).GT.0.0)THEN
        CO2F=0._r8
        CH2O=0._r8
        IF(IGTYP(NZ,NY,NX).NE.0.OR.WFNC.GT.0.0)THEN
!
!         FOR EACH NODE
!
          DO 100 K=1,25
            CH2O3(K)=0._r8
            CH2O4(K)=0._r8
            IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!             C4 PHOTOSYNTHESIS
!
!             ARLF,ARLFL=leaf area
!             ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!             VCGR4=PEP carboxylation rate unlimited by CO2
!
              IF(ICTYP(NZ,NY,NX).EQ.4.AND.VCGR4(K,NB,NZ,NY,NX).GT.0.0)THEN
!
                CALL ComputeGPP_C4(K,NB,NZ,NY,NX)
!
!               C3 PHOTOSYNTHESIS
!
              ELSEIF(ICTYP(NZ,NY,NX).NE.4.AND.VCGRO(K,NB,NZ,NY,NX).GT.0.0)THEN
                call ComputeGPP_C3(K,NB,NZ,NY,NX)

              ENDIF
            ENDIF
100       CONTINUE
!
!         CO2F,CH2O=total CO2 fixation,CH2O production
!
          CO2F=CO2F*0.0432
          CH2O=CH2O*0.0432
!
!         CONVERT UMOL M-2 S-1 TO G C M-2 H-1
!
          DO 150 K=1,25
            CH2O3(K)=CH2O3(K)*0.0432
            CH2O4(K)=CH2O4(K)*0.0432
150       CONTINUE
        ELSE
          CO2F=0._r8
          CH2O=0._r8
          IF(ICTYP(NZ,NY,NX).EQ.4)THEN
            DO 155 K=1,25
              CH2O3(K)=0._r8
              CH2O4(K)=0._r8
155         CONTINUE
          ENDIF
        ENDIF
      ELSE
        CO2F=0._r8
        CH2O=0._r8
        IF(ICTYP(NZ,NY,NX).EQ.4)THEN
          DO 160 K=1,25
            CH2O3(K)=0._r8
            CH2O4(K)=0._r8
160       CONTINUE
        ENDIF
      ENDIF
    ELSE
      CO2F=0._r8
      CH2O=0._r8
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
        DO 165 K=1,25
          CH2O3(K)=0._r8
          CH2O4(K)=0._r8
165     CONTINUE
      ENDIF
    ENDIF
!
!   SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
!
    call ComputeRAutoAfEmergence(NB,NZ,NY,NX)


!   SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
!
  ELSE
    call ComputeRAutoBfEmergence(NB,NZ,NY,NX)
  ENDIF
  end subroutine ComputeGPP
!------------------------------------------------------------------------------------------

      subroutine CalcPartitionCoeff(I,J,NB,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NB,NZ,NY,NX
!     begin_execution
!
!     PARTITION GROWTH WITHIN EACH BRANCH FROM GROWTH STAGE
!     1=LEAF,2=SHEATH OR PETIOLE,3=STALK,4=RESERVE,
!     5,6=REPRODUCTIVE ORGANS,7=GRAIN
!
!     PART=organ partitioning fraction
!
      ARSTKB(NB)=0._r8
      TOTAL=0._r8
      DO 10 N=1,7
      PART(N)=0._r8
10    CONTINUE
!
!     IF BEFORE FLORAL INDUCTION
!
!     IDAY(2,=floral initiation date
!
      IF(IDAY(2,NB,NZ,NY,NX).EQ.0)THEN
      PART(1)=0.725
      PART(2)=0.275
!
!     IF BEFORE ANTHESIS
!
!     IDAY(6,=start of anthesis and setting final seed number
!     TGSTGI=total change in vegv node number normalized for maturity group
!
      ELSEIF(IDAY(6,NB,NZ,NY,NX).EQ.0)THEN
      PART(1)=AMAX1(PART1X,0.725-FPART1*TGSTGI(NB,NZ,NY,NX))
      PART(2)=AMAX1(PART2X,0.275-FPART2*TGSTGI(NB,NZ,NY,NX))
      PARTS=1.0_r8-PART(1)-PART(2)
      PART(3)=0.60*PARTS
      PART(4)=0.30*PARTS
      PARTX=PARTS-PART(3)-PART(4)
      PART(5)=0.5*PARTX
      PART(6)=0.5*PARTX
!
!     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     IDAY(7,=start of grain filling and setting max seed size
!     IDTYP=growth habit:0=determinate,1=indetermimate from PFT file
!     TGSTGF=total change in reprv node number normalized for maturity group
!
      ELSEIF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
      IF(IDTYP(NZ,NY,NX).EQ.0)THEN
      PART(1)=0._r8
      PART(2)=0._r8
      ELSE
      PART(1)=AMAX1(PART1X,(0.725-FPART1)*(1.0-TGSTGF(NB,NZ,NY,NX)))
      PART(2)=AMAX1(PART2X,(0.275-FPART2)*(1.0-TGSTGF(NB,NZ,NY,NX)))
      ENDIF
      PARTS=1.0_r8-PART(1)-PART(2)
      PART(3)=AMAX1(0.0,0.60*PARTS*(1.0-TGSTGF(NB,NZ,NY,NX)))
      PART(4)=AMAX1(0.0,0.30*PARTS*(1.0-TGSTGF(NB,NZ,NY,NX)))
      PARTX=PARTS-PART(3)-PART(4)
      PART(5)=0.5*PARTX
      PART(6)=0.5*PARTX
!
!     DURING GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDTYP=growth habit:0=determinate,1=indetermimate
!
      ELSE
      IF(IDTYP(NZ,NY,NX).EQ.0)THEN
      PART(7)=1.0_r8
      ELSE
      PART(1)=PART1X
      PART(2)=PART2X
      PARTS=1.0_r8-PART(1)-PART(2)
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      PART(3)=0.125*PARTS
      PART(5)=0.125*PARTS
      PART(6)=0.125*PARTS
      PART(7)=0.625*PARTS
      ELSE
      PART(3)=0.75*PARTS
      PART(7)=0.25*PARTS
      ENDIF
      ENDIF
      ENDIF
!
!     IF AFTER GRAIN FILLING
!
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IDAY(10,=physiological maturity date
!
      IF(IBTYP(NZ,NY,NX).EQ.0.AND.IDAY(10,NB,NZ,NY,NX).NE.0)THEN
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      PART(4)=0._r8
      PART(3)=0._r8
      PART(7)=0._r8
      ELSE
      PART(4)=PART(4)+PART(3)
      PART(3)=0._r8
      PART(7)=0._r8
      ENDIF
      ENDIF
!
!     REDIRECT FROM STALK TO STALK RESERVES IF RESERVES BECOME LOW
!
!     WTRSVB,WVSTKB=stalk reserve,sapwood mass
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!
      IF(IDAY(2,NB,NZ,NY,NX).NE.0)THEN
      IF(WTRSVB(NB,NZ,NY,NX).LT.XFRX*WVSTKB(NB,NZ,NY,NX))THEN
      DO 1020 N=1,7
      IF(N.NE.4)THEN
      PART(4)=PART(4)+0.10*PART(N)
      PART(N)=PART(N)-0.10*PART(N)
      ENDIF
1020  CONTINUE
!
!     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
!
      ELSEIF(WTRSVB(NB,NZ,NY,NX).GT.1.0*WVSTKB(NB,NZ,NY,NX))THEN
      PART(3)=PART(3)+PART(4)+PART(7)
      PART(4)=0._r8
      PART(7)=0._r8
      ENDIF
      ENDIF
!
!     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
!
!     ARLFP=PFT leaf area
!
      ARLFI=ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      IF(ARLFI.GT.5.0)THEN
      FPARTL=AMAX1(0.0,(10.0-ARLFI)/5.0)
      PART(3)=PART(3)+(1.0-FPARTL)*(PART(1)+PART(2))
      PART(1)=FPARTL*PART(1)
      PART(2)=FPARTL*PART(2)
      ENDIF
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      PTRT=PART(1)+PART(2)
      ENDIF
!
!     DECIDUOUS LEAF FALL AFTER GRAIN FILL IN DETERMINATES,
!     AFTER AUTUMNIZATION IN INDETERMINATES, OR AFTER SUSTAINED
!     WATER STRESS
!
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     FVRN=fraction of hours required for leafoff to initiate remobilization
!     IDAY(8,=end date for setting final seed number
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IFLGY,IFLGZ=remobilization flags
!     FLGZ=control rate of remobilization
!
      IF((ISTYP(NZ,NY,NX).NE.0.AND.VRNF(NB,NZ,NY,NX) &
      .GE.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX)) &
      .OR.(ISTYP(NZ,NY,NX).EQ.0 &
      .AND.IDAY(8,NB,NZ,NY,NX).NE.0))THEN
      IFLGZ=1
      IF(ISTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.0)THEN
      IFLGY=1
      FLGZ(NB,NZ,NY,NX)=FLGZ(NB,NZ,NY,NX)+1.0
      ELSEIF((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3) &
      .AND.TCC(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
      IFLGY=1
      FLGZ(NB,NZ,NY,NX)=FLGZ(NB,NZ,NY,NX)+1.0
      ELSEIF(IWTYP(NZ,NY,NX).GE.2 &
      .AND.PSILT(NZ,NY,NX).LT.PSILY(IGTYP(NZ,NY,NX)))THEN
      IFLGY=1
      FLGZ(NB,NZ,NY,NX)=FLGZ(NB,NZ,NY,NX)+1.0
      ENDIF
      IF(ISTYP(NZ,NY,NX).NE.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      PART(3)=PART(3)+0.5*(PART(1)+PART(2))
      PART(4)=PART(4)+0.5*(PART(1)+PART(2))
      PART(1)=0._r8
      PART(2)=0._r8
      ENDIF
      ELSE
      IFLGZ=0
      IFLGY=0
      FLGZ(NB,NZ,NY,NX)=0._r8
      ENDIF
!
!     CHECK PARTITIONING COEFFICIENTS
!
      DO 1000 N=1,7
      IF(N.EQ.3.AND.test_aeqb(SNL1(NZ,NY,NX),0._r8))THEN
      PART(N)=0._r8
      ELSE
      PART(N)=AMAX1(0.0,PART(N))
      ENDIF
      TOTAL=TOTAL+PART(N)
1000  CONTINUE
      IF(TOTAL.GT.ZERO)THEN
      DO 1010 N=1,7
      PART(N)=PART(N)/TOTAL
1010  CONTINUE
      ELSE
      DO 1015 N=1,7
      PART(N)=0._r8
1015  CONTINUE
      ENDIF
      end subroutine CalcPartitionCoeff
!------------------------------------------------------------------------------------------

  subroutine GrowOneBranch(I,J,NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NB,NZ,NY,NX
! begin_execution

  WTLSB(NB,NZ,NY,NX)=AMAX1(0.0_r8,WTLFB(NB,NZ,NY,NX)+WTSHEB(NB,NZ,NY,NX))

  IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
    call CalcPartitionCoeff(I,J,NB,NZ,NY,NX)
!
!   SHOOT COEFFICIENTS FOR GROWTH RESPIRATION AND N,P CONTENTS
!   FROM GROWTH YIELDS ENTERED IN 'READQ', AND FROM PARTITIONING
!   COEFFICIENTS ABOVE
!
!   DM*B=C production vs nonstructural C consumption
!   CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!   *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
!   *EAR=ear,*GR=grain from PFT file,*SH=shoot
!   DMSHT=branch C production vs nonstructural C consumption
!   DMSHD=branch C respiration vs nonstructural C consumption
!   CN*M,CP*M=min N,P production vs nonstructural C consumption
!   CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
!   CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!
    IF(IDAY(1,NB,NZ,NY,NX).NE.0)THEN
      DMLFB=DMLF(NZ,NY,NX)
      DMSHB=DMSHE(NZ,NY,NX)
      CNLFB=CNLFW
      CNSHB=CNSHW
      CPLFB=CPLFW
      CPSHB=CPSHW
    ELSE
      DMLFB=DMRT(NZ,NY,NX)
      DMSHB=DMRT(NZ,NY,NX)
      CNLFB=CNRTW
      CNSHB=CNRTW
      CPLFB=CPRTW
      CPSHB=CPRTW
    ENDIF
    DMSHT=PART(1)*DMLFB+PART(2)*DMSHB+PART(3)*DMSTK(NZ,NY,NX) &
      +PART(4)*DMRSV(NZ,NY,NX)+PART(5)*DMHSK(NZ,NY,NX) &
      +PART(6)*DMEAR(NZ,NY,NX)+PART(7)*DMGR(NZ,NY,NX)
    DMSHD=1.0_r8-DMSHT
    CNLFM=PART(1)*DMLFB*ZPLFM*CNLFB
    CPLFM=PART(1)*DMLFB*ZPLFM*CPLFB
    CNLFX=PART(1)*DMLFB*ZPLFD*CNLFB
    CPLFX=PART(1)*DMLFB*ZPLFD*CPLFB
    CNSHX=PART(2)*DMSHB*CNSHB &
      +PART(3)*DMSTK(NZ,NY,NX)*CNSTK(NZ,NY,NX) &
      +PART(4)*DMRSV(NZ,NY,NX)*CNRSV(NZ,NY,NX) &
      +PART(5)*DMHSK(NZ,NY,NX)*CNHSK(NZ,NY,NX) &
      +PART(6)*DMEAR(NZ,NY,NX)*CNEAR(NZ,NY,NX) &
      +PART(7)*DMGR(NZ,NY,NX)*CNRSV(NZ,NY,NX)
    CPSHX=PART(2)*DMSHB*CPSHB &
      +PART(3)*DMSTK(NZ,NY,NX)*CPSTK(NZ,NY,NX) &
      +PART(4)*DMRSV(NZ,NY,NX)*CPRSV(NZ,NY,NX) &
      +PART(5)*DMHSK(NZ,NY,NX)*CPHSK(NZ,NY,NX) &
      +PART(6)*DMEAR(NZ,NY,NX)*CPEAR(NZ,NY,NX) &
      +PART(7)*DMGR(NZ,NY,NX)*CPRSV(NZ,NY,NX)
!
!   TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
!
!   WTSHXN=shoot structural N mass
!   WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,petiole,husk,ear,grain N mass
!   CNSTK,WVSTKB=stalk N:C,sapwood mass
!   IDAY(10=date of physiological maturity
!
    WTSHXN=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)+WTSHBN(NB,NZ,NY,NX) &
      +CNSTK(NZ,NY,NX)*WVSTKB(NB,NZ,NY,NX))
    IF(IDAY(10,NB,NZ,NY,NX).EQ.0)THEN
      WTSHXN=WTSHXN+AMAX1(0.0,WTHSBN(NB,NZ,NY,NX) &
        +WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX))
    ENDIF
!
!   GROSS PRIMARY PRODUCTIVITY
!
    call ComputeGPP(I,J,NB,NZ,NY,NX)
!
!   REMOVE C,N,P USED IN MAINTENANCE + GROWTH REPIRATION AND GROWTH
!   FROM NON-STRUCTURAL POOLS
!
!   CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!   CH2O=total CH2O production
!   RMNCS=maintenance respiration
!   RCO2C=respiration from non-structural C
!   CGROS=total non-structural C used in growth and respiration
!   CNRDA=respiration for N assimilation
!   ZADDB,PADDB=nonstructural N,P used in growth
!   RNH3B=NH3 flux between atmosphere and branch from uptake.f
!   XFRC,XFRN,XFRP=branch-root layer C,N,P transfer
!
    CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+CH2O-AMIN1(RMNCS,RCO2C)-CGROS-CNRDA
    ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-ZADDB+RNH3B(NB,NZ,NY,NX)
    PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-PADDB
!
!   TRANSFER OF C4 FIXATION PRODUCTS FROM NON-STRUCTURAL POOLS
!   IN MESOPHYLL TO THOSE IN BUNDLE SHEATH, DECARBOXYLATION
!   OF C4 FIXATION PRODUCTS IN BUNDLE SHEATH, LEAKAGE OF DECARBOXYLATION
!   PRODUCTS BACK TO MESOPHYLL IN C4 PLANTS
!
!   ICTYP=photosynthesis type:3=C3,4=C4
!
    IF(ICTYP(NZ,NY,NX).EQ.4)THEN
      call C4PhotoProductTransfer(I,J,NZ,NY,NX)
    ENDIF
!
!   C,N,P GROWTH OF LEAF, SHEATH OR PETIOLE, STALK,
!   STALK RESERVES, REPRODUCTIVE ORGANS, GRAIN
!
!   GRO*,GRO*N,GRO*P=organ C,N,P growth rate
!   DM*=C production vs nonstructural C consumption
!   organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
!   HSK=husk,EAR=ear,GR=grain,SHT=shoot
!   PART=organ partitioning fraction
!   CGROS=total non-structural C used in growth and growth respiration
!   CN*,CP*=N:C,P:C ratios in plant organs
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!   CNPG=N,P constraint on growth respiration
!   WT*,WT*N,WT*P=organ C,N,P mass
!
    GROLF=PART(1)*CGROS*DMLFB
    GROSHE=PART(2)*CGROS*DMSHB
    GROSTK=PART(3)*CGROS*DMSTK(NZ,NY,NX)
    GRORSV=PART(4)*CGROS*DMRSV(NZ,NY,NX)
    GROHSK=PART(5)*CGROS*DMHSK(NZ,NY,NX)
    GROEAR=PART(6)*CGROS*DMEAR(NZ,NY,NX)
    GROGR=PART(7)*CGROS*DMGR(NZ,NY,NX)
    GROSHT=CGROS*DMSHT
    GROLFN=GROLF*CNLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHN=GROSHE*CNSHB
    GROSTN=GROSTK*CNSTK(NZ,NY,NX)
    GRORSN=GRORSV*CNRSV(NZ,NY,NX)
    GROHSN=GROHSK*CNHSK(NZ,NY,NX)
    GROEAN=GROEAR*CNEAR(NZ,NY,NX)
    GROGRN=GROGR*CNRSV(NZ,NY,NX)
    GROLFP=GROLF*CPLFB*(ZPLFM+ZPLFD*CNPG)
    GROSHP=GROSHE*CPSHB
    GROSTP=GROSTK*CPSTK(NZ,NY,NX)
    GRORSP=GRORSV*CPRSV(NZ,NY,NX)
    GROHSP=GROHSK*CPHSK(NZ,NY,NX)
    GROEAP=GROEAR*CPEAR(NZ,NY,NX)
    GROGRP=GROGR*CPRSV(NZ,NY,NX)
    WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)+GROLF
    WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)+GROSHE
    WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)+GROSTK
    WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+GRORSV
    WTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)+GROHSK
    WTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)+GROEAR
    WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)+GROLFN
    WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)+GROSHN
    WTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)+GROSTN
    WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+GRORSN
    WTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)+GROHSN
    WTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)+GROEAN
    WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)+GROLFP
    WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)+GROSHP
    WTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)+GROSTP
    WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+GRORSP
    WTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)+GROHSP
    WTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)+GROEAP
!
!   ETOLIATION
!
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
!   ETOL=coefficient for etoliation effects on expansion,extension
!
    CCE=AMIN1(safe_adb(CZPOLB(NB,NZ,NY,NX),CZPOLB(NB,NZ,NY,NX)+CCPOLB(NB,NZ,NY,NX)*CNKI) &
      ,safe_adb(CPPOLB(NB,NZ,NY,NX),CPPOLB(NB,NZ,NY,NX)+CCPOLB(NB,NZ,NY,NX)*CPKI))

    ETOL=1.0_r8+CCE
!
!   DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
!
!   MXNOD,MNNOD=max,min node number currently growing
!   KVSTG=integer of most recent leaf number
!   KNOD,GNOD=number of currently growing nodes
!   ALLOCL=fraction of leaf growth allocated to each node
!   GRO,GRON,GROP=leaf C,N,P growth at each node
!   GSLA=allocation of leaf area growth to each node
!   FNOD=scales node number for perennial vegetation (e.g. trees)
!   NNOD=number of concurrently growing nodes
!
    IF(NB.EQ.NB1(NZ,NY,NX).AND.HTCTL(NZ,NY,NX).LE.SDPTH(NZ,NY,NX))THEN
      NNOD1=0
    ELSE
      NNOD1=1
    ENDIF
    IF(GROLF.GT.0.0)THEN
      MXNOD=KVSTG(NB,NZ,NY,NX)
      MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ,NY,NX)+1)
      MXNOD=MAX(MXNOD,MNNOD)
      KNOD=MXNOD-MNNOD+1
      GNOD=KNOD
      ALLOCL=1.0_r8/GNOD
      GRO=ALLOCL*GROLF
      GRON=ALLOCL*GROLFN
      GROP=ALLOCL*GROLFP
      GSLA=ALLOCL*FNOD(NZ,NY,NX)*NNOD(NZ,NY,NX)
!
!     GROWTH AT EACH CURRENT NODE
!
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     GRO,GRON,GROP=leaf C,N,P growth at each node
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
      DO 490 KK=MNNOD,MXNOD
        K=MOD(KK,25)
        IF(K.EQ.0.AND.KK.NE.0)K=25
          WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)+GRO
          WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)+GRON
          WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)+GROP
          WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX) &
            +AMIN1(GRON*CNWS(NZ,NY,NX),GROP*CPWS(NZ,NY,NX))
!
!         SPECIFIC LEAF AREA FUNCTION OF CURRENT LEAF MASS
!         AT EACH NODE
!
!         SLA=specific area of leaf growth
!         ETOL=coefficient for etoliation effects on expansion,extension
!         SLA1=growth in leaf area vs mass from PFT file
!         SLA2=parameter for calculating leaf area expansion
!         WGLF=leaf C mass
!         PP=PFT population
!         GSLA=allocation of leaf area growth to each node
!         WFNS=turgor expansion,extension function
!         GROA,GRO=leaf area,mass growth
!         ARLFB,ARLF=branch,node leaf area
!
          SLA=ETOL*SLA1(NZ,NY,NX)*(AMAX1(ZEROL(NZ,NY,NX) &
            ,WGLF(K,NB,NZ,NY,NX))/(PP(NZ,NY,NX)*GSLA))**SLA2*WFNS
          GROA=GRO*SLA
          ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)+GROA
          ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)+GROA
490     CONTINUE
      ENDIF
!
!     DISTRIBUTE SHEATH OR PETIOLE GROWTH AMONG CURRENTLY GROWING NODES
!
!     MXNOD,MNNOD=max,min node number currently growing
!     KVSTG=integer of most recent leaf number
!     GNOD=number of currently growing nodes
!     ALLOCS=fraction of petiole growth allocated to each node
!     GRO,GRON,GROP=petiole C,N,P growth at each node
!     GSSL=allocation of petiole length growth to each node
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
      IF(GROSHE.GT.0.0)THEN
        MXNOD=KVSTG(NB,NZ,NY,NX)
        MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ,NY,NX)+1)
        MXNOD=MAX(MXNOD,MNNOD)
        GNOD=MXNOD-MNNOD+1
        ALLOCS=1.0_r8/GNOD
        GRO=ALLOCS*GROSHE
        GRON=ALLOCS*GROSHN
        GROP=ALLOCS*GROSHP
        GSSL=ALLOCL*FNOD(NZ,NY,NX)*NNOD(NZ,NY,NX)
!
!       GROWTH AT EACH CURRENT NODE
!
!       WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!       GRO,GRON,GROP=petiole C,N,P growth at each node
!       CNWS,CPWS=protein:N,protein:P ratios from startq.f
!
        DO 505 KK=MNNOD,MXNOD
          K=MOD(KK,25)
          IF(K.EQ.0.AND.KK.NE.0)K=25
            WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)+GRO
            WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)+GRON
            WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)+GROP
            WSSHE(K,NB,NZ,NY,NX)=WSSHE(K,NB,NZ,NY,NX) &
              +AMIN1(GRON*CNWS(NZ,NY,NX),GROP*CPWS(NZ,NY,NX))
!
!           SPECIFIC SHEATH OR PETIOLE LENGTH FUNCTION OF CURRENT MASS
!           AT EACH NODE
        !
        !   SSL=specific length of petiole growth
        !   ETOL=coefficient for etoliation effects on expansion,extension
        !   SSL1=growth in petiole length vs mass from PFT file
        !   SSL2=parameter for calculating petiole extension
        !   WGSHE=petiole C mass
        !   PP=PFT population
        !   GSSL=allocation of petiole length growth to each node
        !   WFNS=turgor expansion,extension function
        !   GROS,GRO=petiole length,mass growth
        !   HTSHE=petiole length
!
            IF(WGLF(K,NB,NZ,NY,NX).GT.0.0)THEN
              SSL=ETOL*SSL1(NZ,NY,NX)*(AMAX1(ZEROL(NZ,NY,NX) &
                ,WGSHE(K,NB,NZ,NY,NX))/(PP(NZ,NY,NX)*GSSL))**SSL2*WFNS
              GROS=GRO/PP(NZ,NY,NX)*SSL
              HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)+GROS*ANGSH(NZ,NY,NX)
        !     IF(I.EQ.120.AND.J.EQ.24)THEN
        !     WRITE(*,2526)'HTSHE',I,J,NZ,NB,K,SSL,WGSHE(K,NB,NZ,NY,NX)
        !    2,HTSHE(K,NB,NZ,NY,NX),PP(NZ,NY,NX),SSL1(NZ,NY,NX)
        !    3,GSLA,SSL3,WFNS,GROS,GRO,ANGSH(NZ,NY,NX),ZEROL(NZ,NY,NX)
        !    4,CCPOLB(NB,NZ,NY,NX),ETOL
        !2526  FORMAT(A8,5I4,20E12.4)
       !     ENDIF
            ENDIF
505       CONTINUE
        ENDIF
!
    !   DISTRIBUTE STALK GROWTH AMONG CURRENTLY GROWING NODES
    !
    !   MXNOD,MNNOD=max,min node number currently growing
    !   KVSTG=integer of most recent leaf number
    !   GNOD=number of currently growing nodes
    !   ALLOCN=fraction of stalk growth allocated to each node
    !   GRO,GRON,GROP=stalk C,N,P growth at each node
!
        IF(IDAY(1,NB,NZ,NY,NX).EQ.0)THEN
          NN=0
        ELSE
          NN=1
        ENDIF
        MXNOD=KVSTG(NB,NZ,NY,NX)
        MNNOD=MAX(MIN(NN,MAX(NN,MXNOD-NNOD(NZ,NY,NX))) &
          ,KVSTG(NB,NZ,NY,NX)-23)
        MXNOD=MAX(MXNOD,MNNOD)
        IF(GROSTK.GT.0.0)THEN
          GNOD=MXNOD-MNNOD+1
          ALLOCN=1.0_r8/GNOD
          GRO=ALLOCN*GROSTK
          GRON=ALLOCN*GROSTN
          GROP=ALLOCN*GROSTP
    !
    !     SPECIFIC INTERNODE LENGTH FUNCTION OF CURRENT STALK MASS
    !     AT EACH NODE
    !
    !     SNL=specific length of stalk growth
    !     ETOL=coefficient for etoliation effects on expansion,extension
    !     SNL1=growth in stalk length vs mass from PFT file
    !     SNL2=parameter for calculating stalk extension
    !     WTSKB=stalk C mass
    !     PP=PFT population
    !     GROH,GRO=stalk length,mass growth
!
          SNL=ETOL*SNL1(NZ,NY,NX)*(WTSTKB(NB,NZ,NY,NX)/PP(NZ,NY,NX))**SNL2
          GROH=GRO/PP(NZ,NY,NX)*SNL
          KX=MOD(MNNOD-1,25)
          IF(KX.EQ.0.AND.MNNOD-1.NE.0)KX=25
!
    !     GROWTH AT EACH CURRENT NODE
    !
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !     GRO,GRON,GROP=stalk C,N,P growth at each node
    !     HTNODX,HTNODE=stalk height,stalk internode length
    !     ANGBR=sine of stalk angle from horizontal from PFT file
!
          DO 510 KK=MNNOD,MXNOD
            K1=MOD(KK,25)
            IF(K1.EQ.0.AND.KK.NE.0)K1=25
            K2=MOD(KK-1,25)
            IF(K2.EQ.0.AND.KK-1.NE.0)K2=25
            WGNODE(K1,NB,NZ,NY,NX)=WGNODE(K1,NB,NZ,NY,NX)+GRO
            WGNODN(K1,NB,NZ,NY,NX)=WGNODN(K1,NB,NZ,NY,NX)+GRON
            WGNODP(K1,NB,NZ,NY,NX)=WGNODP(K1,NB,NZ,NY,NX)+GROP
            HTNODX(K1,NB,NZ,NY,NX)=HTNODX(K1,NB,NZ,NY,NX)+GROH*ANGBR(NZ,NY,NX)
            IF(K1.NE.0)THEN
              HTNODE(K1,NB,NZ,NY,NX)=HTNODX(K1,NB,NZ,NY,NX) &
                +HTNODE(K2,NB,NZ,NY,NX)
            ELSE
              HTNODE(K1,NB,NZ,NY,NX)=HTNODX(K1,NB,NZ,NY,NX)
            ENDIF
        !   IF(NZ.EQ.1)THEN
        !     WRITE(*,515)'HTNODE',I,J,NX,NY,NZ,NB,KK,K1,K2,MNNOD,MXNOD
        !    1,NNOD(NZ,NY,NX),ARLF(K1,NB,NZ,NY,NX)
        !    2,HTNODE(K1,NB,NZ,NY,NX),HTNODE(K2,NB,NZ,NY,NX),SNL,GRO
        !    3,ALLOCN,WTSTKB(NB,NZ,NY,NX),WGNODE(K1,NB,NZ,NY,NX)
        !    4,HTNODX(K1,NB,NZ,NY,NX),PP(NZ,NY,NX),GROSTK
        !515   FORMAT(A8,12I4,20E12.4)
        !   ENDIF
510       CONTINUE
        ENDIF
!
    !   RECOVERY OF REMOBILIZABLE N,P DURING REMOBILIZATION DEPENDS
    !   ON SHOOT NON-STRUCTURAL C:N:P
    !
    !   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
    !   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
    !   RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !   RCCZ,RCCY=min,max fractions for shoot C recycling
    !   RCCX,RCCQ=max fractions for shoot N,P recycling
    !   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
        IF(IDAY(1,NB,NZ,NY,NX).NE.0 &
         .AND.CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
          CCC=AMAX1(0.0,AMIN1(1.0 &
            ,CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX) &
            +CCPOLB(NB,NZ,NY,NX)*CNKI) &
            ,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX) &
            +CCPOLB(NB,NZ,NY,NX)*CPKI)))
          CNC=AMAX1(0.0,AMIN1(1.0 &
            ,CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX) &
            +CZPOLB(NB,NZ,NY,NX)/CNKI)))
          CPC=AMAX1(0.0,AMIN1(1.0 &
            ,CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX) &
            +CPPOLB(NB,NZ,NY,NX)/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        RCCC=RCCZ(IGTYP(NZ,NY,NX))+CCC*RCCY(IGTYP(NZ,NY,NX))
        RCCN=CNC*RCCX(IGTYP(NZ,NY,NX))
        RCCP=CPC*RCCQ(IGTYP(NZ,NY,NX))
!
!       WITHDRAW REMOBILIZABLE C,N,P FROM LOWEST NODE AFTER
!       MAXIMUM NODE NUMBER OF 25 IS REACHED
!
!       IFLGG=PFT senescence flag
!       KVSTG=integer of most recent leaf number
!       TFN3=temperature function for canopy growth
!       XRLA=rate of leaf appearance at 25 oC (h-1)
!       FSNC=fraction of lowest leaf to be remobilized
!
        IF(IFLGG(NB,NZ,NY,NX).EQ.1)THEN
          KVSTGX=MAX(0,KVSTG(NB,NZ,NY,NX)-24)
          IF(KVSTGX.GT.0)THEN
            K=MOD(KVSTGX,25)
            IF(K.EQ.0.AND.KVSTGX.GT.0)K=25
            KX=MOD(KVSTG(NB,NZ,NY,NX),25)
            IF(KX.EQ.0.AND.KVSTG(NB,NZ,NY,NX).NE.0)KX=25
            FSNC=TFN3(NZ,NY,NX)*XRLA(NZ,NY,NX)
!
        !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
        !
        !   IFLGP=flag for remobilization
        !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !   ARLF=node leaf area
        !   RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
!
            IF(IFLGP(NB,NZ,NY,NX).EQ.1)THEN
              WGLFX(NB,NZ,NY,NX)=AMAX1(0.0,WGLF(K,NB,NZ,NY,NX))
              WGLFNX(NB,NZ,NY,NX)=AMAX1(0.0,WGLFN(K,NB,NZ,NY,NX))
              WGLFPX(NB,NZ,NY,NX)=AMAX1(0.0,WGLFP(K,NB,NZ,NY,NX))
              ARLFZ(NB,NZ,NY,NX)=AMAX1(0.0,ARLF(K,NB,NZ,NY,NX))
              IF(WGLFX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                RCCLX(NB,NZ,NY,NX)=RCCC*WGLFX(NB,NZ,NY,NX)
                RCZLX(NB,NZ,NY,NX)=WGLFNX(NB,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
                RCPLX(NB,NZ,NY,NX)=WGLFPX(NB,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
              ELSE
                RCCLX(NB,NZ,NY,NX)=0._r8
                RCZLX(NB,NZ,NY,NX)=0._r8
                RCPLX(NB,NZ,NY,NX)=0._r8
              ENDIF
            ENDIF
!
    !       FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !       FSNC,FSNCL=fraction of lowest leaf to be remobilized
    !
            IF(FSNC*WGLFX(NB,NZ,NY,NX).GT.WGLF(K,NB,NZ,NY,NX) &
              .AND.WGLFX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
              FSNCL=AMAX1(0.0,WGLF(K,NB,NZ,NY,NX)/WGLFX(NB,NZ,NY,NX))
            ELSE
              FSNCL=FSNC
            ENDIF
!
    !       NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !       TO FRACTIONS SET IN 'STARTQ'
    !
    !       CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !       CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !       FSNCL=fraction of lowest leaf to be remobilized
    !       RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
    !       WGLFX,WGLFNX,WGLFPX=senescing leaf C,N,P mass
    !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !       FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
            DO 6300 M=1,4
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
                *FSNCL*(WGLFX(NB,NZ,NY,NX)-RCCLX(NB,NZ,NY,NX))*FWODB(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
                *FSNCL*(WGLFNX(NB,NZ,NY,NX)-RCZLX(NB,NZ,NY,NX))*FWODLN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
                *FSNCL*(WGLFPX(NB,NZ,NY,NX)-RCPLX(NB,NZ,NY,NX))*FWODLP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX) &
                *FSNCL*(WGLFX(NB,NZ,NY,NX)-RCCLX(NB,NZ,NY,NX))*FWODB(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX) &
                *FSNCL*(WGLFNX(NB,NZ,NY,NX)-RCZLX(NB,NZ,NY,NX))*FWODLN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX) &
                *FSNCL*(WGLFPX(NB,NZ,NY,NX)-RCPLX(NB,NZ,NY,NX))*FWODLP(1)
6300        CONTINUE
!
    !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !       FSNCL=fraction of lowest leaf to be remobilized
    !       ARLFB,ARLFZ=branch living,senescing leaf area
    !       WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in living,senescing leaf
    !       WSLF=leaf protein mass
    !       CNWS,CPWS=protein:N,protein:P ratios from startq.f
    !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
    !       RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
!
            ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX) &
              -FSNCL*ARLFZ(NB,NZ,NY,NX)
            WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX) &
              -FSNCL*WGLFX(NB,NZ,NY,NX)
            WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX) &
              -FSNCL*WGLFNX(NB,NZ,NY,NX)
            WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX) &
              -FSNCL*WGLFPX(NB,NZ,NY,NX)
            ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX) &
              -FSNCL*ARLFZ(NB,NZ,NY,NX)
            WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX) &
              -FSNCL*WGLFX(NB,NZ,NY,NX)
            WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX) &
              -FSNCL*WGLFNX(NB,NZ,NY,NX)
            WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX) &
              -FSNCL*WGLFPX(NB,NZ,NY,NX)
            WSLF(K,NB,NZ,NY,NX)=AMAX1(0.0,WSLF(K,NB,NZ,NY,NX) &
              -FSNCL*AMAX1(WGLFNX(NB,NZ,NY,NX)*CNWS(NZ,NY,NX) &
              ,WGLFPX(NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
            CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+FSNCL*RCCLX(NB,NZ,NY,NX)
            ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+FSNCL*RCZLX(NB,NZ,NY,NX)
            PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+FSNCL*RCPLX(NB,NZ,NY,NX)
!
    !       REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P ALSO DEPENDS ON
    !       STRUCTURAL C:N:P
    !
    !       IFLGP=flag for remobilization
    !       WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !       HTSHEX=petiole length
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
!
            IF(IFLGP(NB,NZ,NY,NX).EQ.1)THEN
              WGSHEX(NB,NZ,NY,NX)=AMAX1(0.0,WGSHE(K,NB,NZ,NY,NX))
              WGSHNX(NB,NZ,NY,NX)=AMAX1(0.0,WGSHN(K,NB,NZ,NY,NX))
              WGSHPX(NB,NZ,NY,NX)=AMAX1(0.0,WGSHP(K,NB,NZ,NY,NX))
              HTSHEX(NB,NZ,NY,NX)=AMAX1(0.0,HTSHE(K,NB,NZ,NY,NX))
              IF(WGSHEX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
                RCCSX(NB,NZ,NY,NX)=RCCC*WGSHEX(NB,NZ,NY,NX)
                RCZSX(NB,NZ,NY,NX)=WGSHNX(NB,NZ,NY,NX) &
                  *(RCCN+(1.0-RCCN)*RCCSX(NB,NZ,NY,NX)/WGSHEX(NB,NZ,NY,NX))
                RCPSX(NB,NZ,NY,NX)=WGSHPX(NB,NZ,NY,NX) &
                  *(RCCP+(1.0-RCCP)*RCCSX(NB,NZ,NY,NX)/WGSHEX(NB,NZ,NY,NX))
              ELSE
                RCCSX(NB,NZ,NY,NX)=0._r8
                RCZSX(NB,NZ,NY,NX)=0._r8
                RCPSX(NB,NZ,NY,NX)=0._r8
              ENDIF
              WTSTXB(NB,NZ,NY,NX)=WTSTXB(NB,NZ,NY,NX)+WGNODE(K,NB,NZ,NY,NX)
              WTSTXN(NB,NZ,NY,NX)=WTSTXN(NB,NZ,NY,NX)+WGNODN(K,NB,NZ,NY,NX)
              WTSTXP(NB,NZ,NY,NX)=WTSTXP(NB,NZ,NY,NX)+WGNODP(K,NB,NZ,NY,NX)
        !     IF(NZ.EQ.2)THEN
        !     WRITE(*,2358)'WTSTXB',I,J,NZ,NB,K,WTSTXB(NB,NZ,NY,NX)
        !    2,WTSTKB(NB,NZ,NY,NX),WGNODE(K,NB,NZ,NY,NX)
        !2358  FORMAT(A8,5I4,12E12.4)
        !     ENDIF
              WGNODE(K,NB,NZ,NY,NX)=0._r8
              WGNODN(K,NB,NZ,NY,NX)=0._r8
              WGNODP(K,NB,NZ,NY,NX)=0._r8
              HTNODX(K,NB,NZ,NY,NX)=0._r8
            ENDIF
!
    !       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !
            IF(FSNC*WGSHEX(NB,NZ,NY,NX).GT.WGSHE(K,NB,NZ,NY,NX) &
              .AND.WGSHEX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
              FSNCS=AMAX1(0.0,WGSHE(K,NB,NZ,NY,NX)/WGSHEX(NB,NZ,NY,NX))
            ELSE
              FSNCS=FSNC
            ENDIF
!
    !       NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !       TO FRACTIONS SET IN 'STARTQ'
    !
    !       CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !       CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !       FSNCS=fraction of lowest petiole to be remobilized
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
    !       WGSHX,WGSHNX,WGSHPX=senescing petiole C,N,P mass
    !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !       FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!
            DO 6305 M=1,4
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
                *FSNCS*(WGSHEX(NB,NZ,NY,NX)-RCCSX(NB,NZ,NY,NX))*FWODB(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
                *FSNCS*(WGSHNX(NB,NZ,NY,NX)-RCZSX(NB,NZ,NY,NX))*FWODSN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
                *FSNCS*(WGSHPX(NB,NZ,NY,NX)-RCPSX(NB,NZ,NY,NX))*FWODSP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX) &
                *FSNCS*(WGSHEX(NB,NZ,NY,NX)-RCCSX(NB,NZ,NY,NX))*FWODB(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX) &
                *FSNCS*(WGSHNX(NB,NZ,NY,NX)-RCZSX(NB,NZ,NY,NX))*FWODSN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX) &
                *FSNCS*(WGSHPX(NB,NZ,NY,NX)-RCPSX(NB,NZ,NY,NX))*FWODSP(1)
6305        CONTINUE
!
    !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !       FSNCS=fraction of lowest petiole to be remobilized
    !       HTSHE,HTSHEX=living,senescing petiole length
    !       WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in living,senescing petiole
    !       WSSHE=petiole protein mass
    !       CNWS,CPWS=protein:N,protein:P ratios from startq.f
    !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !       RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
!
            WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX) &
              -FSNCS*WGSHEX(NB,NZ,NY,NX)
            WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX) &
              -FSNCS*WGSHNX(NB,NZ,NY,NX)
            WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX) &
              -FSNCS*WGSHPX(NB,NZ,NY,NX)
            HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX) &
              -FSNCS*HTSHEX(NB,NZ,NY,NX)
            WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX) &
              -FSNCS*WGSHEX(NB,NZ,NY,NX)
            WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX) &
              -FSNCS*WGSHNX(NB,NZ,NY,NX)
            WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX) &
              -FSNCS*WGSHPX(NB,NZ,NY,NX)
            WSSHE(K,NB,NZ,NY,NX)=AMAX1(0.0,WSSHE(K,NB,NZ,NY,NX) &
              -FSNCS*AMAX1(WGSHNX(NB,NZ,NY,NX)*CNWS(NZ,NY,NX) &
              ,WGSHPX(NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
            CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+FSNCS*RCCSX(NB,NZ,NY,NX)
            ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+FSNCS*RCZSX(NB,NZ,NY,NX)
            PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+FSNCS*RCPSX(NB,NZ,NY,NX)
          ENDIF
        ENDIF
!
  !     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0
  !
  !     SNCR=excess maintenance respiration
  !     WTRSVB=stalk reserve C mass
  !     RCO2V=remobilization of stalk reserve C
  !     VMXC=rate constant for nonstructural C oxidation in respiration
  !     TFN3=temperature function for canopy growth
  !
        IF(IFLGZ.EQ.0)THEN
          IF(SNCR.GT.0.0.AND.WTRSVB(NB,NZ,NY,NX).GT.0.0)THEN
            RCO2V=AMIN1(SNCR,VMXC*WTRSVB(NB,NZ,NY,NX)*TFN3(NZ,NY,NX))
            WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-RCO2V
            SNCR=SNCR-RCO2V
          ENDIF
        ENDIF
!
!       TOTAL REMOBILIZATION = GROWTH RESPIRATION < 0 + DECIDUOUS LEAF
!       FALL DURING AUTUMN + REMOBILZATION DURING GRAIN FILL IN ANNUALS
!
!       ISTYP=growth habit:0=annual,1=perennial from PFT file
!       IFLGY,IFLGZ=remobilization flags
!       SNCZ=phenologically-driven respiration senescence during late-season
!       FXFB=rate constant for plant-storage nonstructural C,N,P exchange
!       IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!       WTLSB=leaf+petiole mass
!       FLGZ=control rate of remobilization
!       FLGZX=number of hours until full senescence after physl maturity
!       SNCX=total senescence respiration
!       KVSTG,KVSTGN=integer of highest,lowest leaf number currently growing
!       KSNC=number of nodes undergoing remobilization
!       SNCF=ratio of phenologically-driven vs total senescence respiration
!
        IF(IFLGZ.EQ.1.AND.IFLGY.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
          SNCZ=FXFB(IBTYP(NZ,NY,NX)) &
            *WTLSB(NB,NZ,NY,NX)*AMIN1(1.0,FLGZ(NB,NZ,NY,NX)/FLGZX)
        ELSE
          SNCZ=0._r8
        ENDIF
        SNCX=SNCR+SNCZ
        IF(SNCX.GT.ZEROP(NZ,NY,NX))THEN
          SNCF=SNCZ/SNCX
          KSNC=INT(0.5*(KVSTG(NB,NZ,NY,NX)-KVSTGN(NB,NZ,NY,NX)))+1
          XKSNC=KSNC
          KN=MAX(0,KVSTGN(NB,NZ,NY,NX)-1)
    !     IF(NZ.EQ.2.OR.NZ.EQ.3)THEN
    !     WRITE(*,1266)'SNCX0',I,J,NX,NY,NZ,NB,SNCY,SNCR,SNCX,SNCF
    !    2,CPOOL(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX),RCCC
    !1266  FORMAT(A8,6I4,12E16.8)
    !     ENDIF
    !
    !     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
    !     IF MAIN STEM POOLS ARE DEPLETED
    !
    !     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
    !     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
    !     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
    !
          IF(IBTYP(NZ,NY,NX).NE.0.AND.IGTYP(NZ,NY,NX).GT.1 &
            .AND.NB.EQ.NB1(NZ,NY,NX).AND.test_aeqb(SNCF,0._r8))THEN
            NBY=0
          DO 584 NBL=1,NBR(NZ,NY,NX)
            NBZ(NBL)=0
584       CONTINUE
          DO 586 NBL=1,NBR(NZ,NY,NX)
            NBX=KVSTG(NB,NZ,NY,NX)
            DO 585 NBK=1,NBR(NZ,NY,NX)
              IF(IDTHB(NBK,NZ,NY,NX).EQ.0.AND.NBK.NE.NB1(NZ,NY,NX) &
                .AND.NBTB(NBK,NZ,NY,NX).LT.NBX &
                .AND.NBTB(NBK,NZ,NY,NX).GT.NBY)THEN
                NBZ(NBL)=NBK
                NBX=NBTB(NBK,NZ,NY,NX)
              ENDIF
585         CONTINUE
            IF(NBZ(NBL).NE.0)THEN
              NBY=NBTB(NBZ(NBL),NZ,NY,NX)
            ENDIF
586       CONTINUE
          DO 580 NBL=1,NBR(NZ,NY,NX)
            IF(NBZ(NBL).NE.0)THEN
              IF(NBTB(NBZ(NBL),NZ,NY,NX).LT.KK)THEN
                IF(CPOOL(NBZ(NBL),NZ,NY,NX).GT.0)THEN
                  XFRC=1.0E-02_r8*AMIN1(SNCX,CPOOL(NBZ(NBL),NZ,NY,NX))
                  XFRN=XFRC*ZPOOL(NBZ(NBL),NZ,NY,NX)/CPOOL(NBZ(NBL),NZ,NY,NX)
                  XFRP=XFRC*PPOOL(NBZ(NBL),NZ,NY,NX)/CPOOL(NBZ(NBL),NZ,NY,NX)
                ELSE
                  XFRC=0._r8
                  XFRN=1.0E-02_r8*ZPOOL(NBZ(NBL),NZ,NY,NX)
                  XFRP=1.0E-02_r8*PPOOL(NBZ(NBL),NZ,NY,NX)
                ENDIF
                CPOOL(NBZ(NBL),NZ,NY,NX)=CPOOL(NBZ(NBL),NZ,NY,NX)-XFRC
                ZPOOL(NBZ(NBL),NZ,NY,NX)=ZPOOL(NBZ(NBL),NZ,NY,NX)-XFRN
                PPOOL(NBZ(NBL),NZ,NY,NX)=PPOOL(NBZ(NBL),NZ,NY,NX)-XFRP
                CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=CPOOL(NB1(NZ,NY,NX),NZ,NY,NX) &
                  +XFRC*SNCF
                ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX) &
                  +XFRN
                PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=PPOOL(NB1(NZ,NY,NX),NZ,NY,NX) &
                  +XFRP
                SNCX=SNCX-XFRC
                IF(SNCX.LE.0.0)GO TO 595
              ENDIF
            ENDIF
580       CONTINUE
        ENDIF
!
  !     REMOBILIZATION AND LITTERFALL WHEN GROWTH RESPIRATION < 0
  !     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
  !
  !     SNCX,SNCT=branch,node senescence respiration
  !     KSNC=number of nodes undergoing remobilization
  !
  !     IF(NZ.EQ.2.OR.NZ.EQ.3)THEN
  !     WRITE(*,1266)'SNCX1',I,J,NX,NY,NZ,NB,SNCY,SNCR,SNCX,SNCF
  !    2,CPOOL(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX),RCCC
  !     ENDIF
        DO 575 N=1,KSNC
          SNCT=SNCX/XKSNC
          DO 650 KK=KN,KVSTG(NB,NZ,NY,NX)
            SNCLF=0._r8
            SNCSH=0._r8
            K=MOD(KK,25)
            IF(K.EQ.0.AND.KK.NE.0)K=25
    !
    !       REMOBILIZATION OF LEAF C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !       WGLF,WGSHE=node leaf,petiole C mass
    !       SCNF,SCNSH=leaf,petiole senescence respiration
    !       RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
    !       RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !       RCCZ,RCCY=min,max fractions for shoot C recycling
    !
            IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
              FNCLF=WGLF(K,NB,NZ,NY,NX)/(WGLF(K,NB,NZ,NY,NX) &
                +WGSHE(K,NB,NZ,NY,NX))
              SNCLF=FNCLF*SNCT
              SNCSH=SNCT-SNCLF
              RCCL=RCCC*WGLF(K,NB,NZ,NY,NX)
              RCZL=WGLFN(K,NB,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
              RCPL=WGLFP(K,NB,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
!
    !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !         FSNCL,FSNAL=fraction of current leaf C,area to be remobilized
    !
              IF(RCCL.GT.ZEROP(NZ,NY,NX))THEN
                FSNCL=AMAX1(0.0,AMIN1(1.0,SNCLF/RCCL))
              ELSE
                FSNCL=1.0_r8
              ENDIF
              FSNAL=FSNCL
!
        !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
        !     TO FRACTIONS SET IN 'STARTQ'
        !
        !     CSNC,ZSNC,PSNC=literfall C,N,P
        !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
        !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
        !     FSNCL=fraction of current leaf to be remobilized
        !     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
        !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
        !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
        !     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
        !
              DO 6310 M=1,4
                CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
                  *FSNCL*(WGLF(K,NB,NZ,NY,NX)-RCCL)*FWODB(0)
                ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
                  *FSNCL*(WGLFN(K,NB,NZ,NY,NX)-RCZL)*FWODLN(0)
                PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
                  *FSNCL*(WGLFP(K,NB,NZ,NY,NX)-RCPL)*FWODLP(0)
                CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX) &
                  *FSNCL*(WGLF(K,NB,NZ,NY,NX)-RCCL)*FWODB(1)
                ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX) &
                  *FSNCL*(WGLFN(K,NB,NZ,NY,NX)-RCZL)*FWODLN(1)
                PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX) &
                  *FSNCL*(WGLFP(K,NB,NZ,NY,NX)-RCPL)*FWODLP(1)
6310          CONTINUE
              IF(K.NE.0)THEN
                CSNC(2,1,0,NZ,NY,NX)=CSNC(2,1,0,NZ,NY,NX) &
                  +FSNCL*(CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX))
                CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX) &
                  -FSNCL*CPOOL3(K,NB,NZ,NY,NX)
                CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX) &
                  -FSNCL*CPOOL4(K,NB,NZ,NY,NX)
              ENDIF
!
        !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
        !
        !     ARLFB=leaf area
        !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
        !     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
        !     FSNCL=fraction of current leaf to be remobilized
        !     CNWS,CPWS=protein:N,protein:P ratios from startq.f
        !
              ARLFB(NB,NZ,NY,NX)=AMAX1(0.0,ARLFB(NB,NZ,NY,NX)-FSNAL*ARLF(K,NB,NZ,NY,NX))
              WTLFB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)-FSNCL*WGLF(K,NB,NZ,NY,NX))
              WTLFBN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)-FSNCL*WGLFN(K,NB,NZ,NY,NX))
              WTLFBP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX)-FSNCL*WGLFP(K,NB,NZ,NY,NX))
              ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)-FSNAL*ARLF(K,NB,NZ,NY,NX)
              WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)-FSNCL*WGLF(K,NB,NZ,NY,NX)
              WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)-FSNCL*WGLFN(K,NB,NZ,NY,NX)
              WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)-FSNCL*WGLFP(K,NB,NZ,NY,NX)
              WSLF(K,NB,NZ,NY,NX)=AMAX1(0.0,WSLF(K,NB,NZ,NY,NX) &
                -FSNCL*AMAX1(WGLFN(K,NB,NZ,NY,NX)*CNWS(NZ,NY,NX) &
                ,WGLFP(K,NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
!
        !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
        !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
        !
        !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
        !     FSNCL=fraction of current leaf C to be remobilized
        !     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
        !     SNCLF,SNCT=remaining senescence respiration carried to next node
        !
              CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+FSNCL*RCCL*SNCF
              ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+FSNCL*RCZL
              PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+FSNCL*RCPL
              SNCLF=SNCLF-FSNCL*RCCL
              SNCT=SNCT-FSNCL*RCCL
              IF(WTLFB(NB,NZ,NY,NX).LE.ZEROL(NZ,NY,NX))THEN
                WTLFB(NB,NZ,NY,NX)=0._r8
                ARLFB(NB,NZ,NY,NX)=0._r8
              ENDIF
      !
        !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
        !
              IF(SNCLF.LE.ZEROP(NZ,NY,NX))GO TO 564
        !
        !     OTHERWISE REMAINING C,N,P IN LEAF GOES TO LITTERFALL
        !
        !     CSNC,ZSNC,PSNC=literfall C,N,P
        !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
        !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
        !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
        !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
        !     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
        !     ARLFB=leaf area
        !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
        !     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
        !
            ELSE
              DO 6315 M=1,4
                CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
                  *WGLF(K,NB,NZ,NY,NX)*FWODB(0)
                ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
                  *WGLFN(K,NB,NZ,NY,NX)*FWODLN(0)
                PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
                  *WGLFP(K,NB,NZ,NY,NX)*FWODLP(0)
                CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX) &
                  *WGLF(K,NB,NZ,NY,NX)*FWODB(1)
                ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX) &
                  *WGLFN(K,NB,NZ,NY,NX)*FWODLN(1)
                PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX) &
                  *WGLFP(K,NB,NZ,NY,NX)*FWODLP(1)
6315          CONTINUE
              IF(K.NE.0)THEN
                CSNC(2,1,0,NZ,NY,NX)=CSNC(2,1,0,NZ,NY,NX) &
                  +CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX)
                CPOOL3(K,NB,NZ,NY,NX)=0._r8
                CPOOL4(K,NB,NZ,NY,NX)=0._r8
              ENDIF
              ARLFB(NB,NZ,NY,NX)=AMAX1(0.0,ARLFB(NB,NZ,NY,NX)-ARLF(K,NB,NZ,NY,NX))
              WTLFB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)-WGLF(K,NB,NZ,NY,NX))
              WTLFBN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)-WGLFN(K,NB,NZ,NY,NX))
              WTLFBP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX)-WGLFP(K,NB,NZ,NY,NX))
              ARLF(K,NB,NZ,NY,NX)=0._r8
              WGLF(K,NB,NZ,NY,NX)=0._r8
              WGLFN(K,NB,NZ,NY,NX)=0._r8
              WGLFP(K,NB,NZ,NY,NX)=0._r8
              WSLF(K,NB,NZ,NY,NX)=0._r8
              IF(WTLFB(NB,NZ,NY,NX).LE.ZEROL(NZ,NY,NX))THEN
                WTLFB(NB,NZ,NY,NX)=0._r8
                ARLFB(NB,NZ,NY,NX)=0._r8
              ENDIF
            ENDIF
564       CONTINUE
    !
    !     REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P DEPENDS ON
    !     NON-STRUCTURAL C:N:P
    !
    !     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
    !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !
          IF(WGSHE(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            RCCS=RCCC*WGSHE(K,NB,NZ,NY,NX)
            RCZS=WGSHN(K,NB,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
            RCPS=WGSHP(K,NB,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
!
      !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
      !     OR PETIOLE
      !
      !     FSNCS,FSNAS=fraction of current petiole C,length to be remobilized
      !
            IF(RCCS.GT.ZEROP(NZ,NY,NX))THEN
              FSNCS=AMAX1(0.0,AMIN1(1.0,SNCSH/RCCS))
            ELSE
              FSNCS=1.0_r8
            ENDIF
            FSNAS=FSNCS
    !
      !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FSNCS=fraction of current petiole to be remobilized
      !     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
      !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
      !
            DO 6320 M=1,4
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
                *FSNCS*(WGSHE(K,NB,NZ,NY,NX)-RCCS)*FWODB(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
                *FSNCS*(WGSHN(K,NB,NZ,NY,NX)-RCZS)*FWODSN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
                *FSNCS*(WGSHP(K,NB,NZ,NY,NX)-RCPS)*FWODSP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX) &
                *FSNCS*(WGSHE(K,NB,NZ,NY,NX)-RCCS)*FWODB(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX) &
                *FSNCS*(WGSHN(K,NB,NZ,NY,NX)-RCZS)*FWODSN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX) &
                *FSNCS*(WGSHP(K,NB,NZ,NY,NX)-RCPS)*FWODSP(1)
6320        CONTINUE
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     HTSHE=petiole length
      !     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !     FSNCS=fraction of current petiole to be remobilized
      !     CNWS,CPWS=protein:N,protein:P ratios from startq.f
      !
            WTSHEB(NB,NZ,NY,NX)=AMAX1(0.0,WTSHEB(NB,NZ,NY,NX) &
             -FSNCS*WGSHE(K,NB,NZ,NY,NX))
            WTSHBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBN(NB,NZ,NY,NX) &
              -FSNCS*WGSHN(K,NB,NZ,NY,NX))
            WTSHBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBP(NB,NZ,NY,NX) &
              -FSNCS*WGSHP(K,NB,NZ,NY,NX))
            HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX) &
              -FSNAS*HTSHE(K,NB,NZ,NY,NX)
            WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX) &
              -FSNCS*WGSHE(K,NB,NZ,NY,NX)
            WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX) &
              -FSNCS*WGSHN(K,NB,NZ,NY,NX)
            WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX) &
              -FSNCS*WGSHP(K,NB,NZ,NY,NX)
            WSSHE(K,NB,NZ,NY,NX)=AMAX1(0.0,WSSHE(K,NB,NZ,NY,NX) &
              -FSNCS*AMAX1(WGSHN(K,NB,NZ,NY,NX)*CNWS(NZ,NY,NX) &
              ,WGSHP(K,NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
    !
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
      !     FSNCS=fraction of current petiole C to be remobilized
      !     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
      !     SNCSH,SNCT=remaining senescence respiration carried to next node
      !
            CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+FSNCS*RCCS*SNCF
            ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+FSNCS*RCZS
            PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+FSNCS*RCPS
            SNCSH=SNCSH-FSNCS*RCCS
            SNCT=SNCT-FSNCS*RCCS
            IF(WTSHEB(NB,NZ,NY,NX).LE.ZEROL(NZ,NY,NX))THEN
              WTSHEB(NB,NZ,NY,NX)=0._r8
            ENDIF
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
            IF(SNCSH.LE.ZEROP(NZ,NY,NX))GO TO 565
      !
      !     OTHERWISE REMAINING C,N,P IN SHEATH OR PETIOLE GOES TO LITTERFALL
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FWODB=C woody fraction in branch:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
      !     HTSHE=petiole length
      !     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
      !     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
      !
          ELSE
          DO 6325 M=1,4
            CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX) &
              *WGSHE(K,NB,NZ,NY,NX)*FWODB(0)
            ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX) &
              *WGSHN(K,NB,NZ,NY,NX)*FWODSN(0)
            PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX) &
              *WGSHP(K,NB,NZ,NY,NX)*FWODSP(0)
            CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX) &
              *WGSHE(K,NB,NZ,NY,NX)*FWODB(1)
            ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX) &
              *WGSHN(K,NB,NZ,NY,NX)*FWODSN(1)
            PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX) &
              *WGSHP(K,NB,NZ,NY,NX)*FWODSP(1)
6325      CONTINUE
          WTSHEB(NB,NZ,NY,NX)=AMAX1(0.0,WTSHEB(NB,NZ,NY,NX) &
            -WGSHE(K,NB,NZ,NY,NX))
          WTSHBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBN(NB,NZ,NY,NX) &
            -WGSHN(K,NB,NZ,NY,NX))
          WTSHBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBP(NB,NZ,NY,NX) &
            -WGSHP(K,NB,NZ,NY,NX))
          HTSHE(K,NB,NZ,NY,NX)=0._r8
          WGSHE(K,NB,NZ,NY,NX)=0._r8
          WGSHN(K,NB,NZ,NY,NX)=0._r8
          WGSHP(K,NB,NZ,NY,NX)=0._r8
          WSSHE(K,NB,NZ,NY,NX)=0._r8
          IF(WTSHEB(NB,NZ,NY,NX).LE.ZEROL(NZ,NY,NX))THEN
            WTSHEB(NB,NZ,NY,NX)=0._r8
          ENDIF
        ENDIF
650   CONTINUE
      KN=KN+1
      SNCR=SNCT*(1.0-SNCF)
!
!     REMOBILIZATION OF RESERVE C
!
!     WTRSVB=stalk reserve C mass
!     SNCR=excess maintenance respiration
!
      IF(WTRSVB(NB,NZ,NY,NX).GT.SNCR)THEN
        WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-SNCR
        SNCR=0._r8
        GO TO 565
      ENDIF
!
!     REMOBILIZATION OF STALK C,N,P
!
!     FXFS=rate constant for remobilization of stalk C,N,P (h-1)
!     SNCZ=phenologically-driven respiration senescence during late-season
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     WTSTKB,WVSTKB=stalk,sapwood C mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     MXNOD,MNNOD=max,min node number currently growing
!     NNOD=number of concurrently growing nodes
!     KVSTG=integer of most recent leaf number
!
      SNCZ=FXFS*SNCR
      SNCT=SNCR+SNCZ
      IF(ISTYP(NZ,NY,NX).NE.0.AND.SNCT.GT.ZEROP(NZ,NY,NX) &
        .AND.WTSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        SNCF=SNCZ/SNCT
        FRCC=WVSTKB(NB,NZ,NY,NX)/WTSTKB(NB,NZ,NY,NX)
        RCSC=RCCC*FRCC
        RCSN=RCCN*FRCC
        RCSP=RCCP*FRCC
        MXNOD=KVSTG(NB,NZ,NY,NX)
        MNNOD=MAX(MIN(0,MAX(0,MXNOD-NNOD(NZ,NY,NX))) &
          ,KVSTG(NB,NZ,NY,NX)-23)
        MXNOD=MAX(MXNOD,MNNOD)
        DO 1650 KK=MXNOD,MNNOD,-1
          K=MOD(KK,25)
          IF(K.EQ.0.AND.KK.NE.0)K=25
    !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
    !     WRITE(*,2356)'WGNODE1',I,J,NX,NY,NZ,NB,K,KK,MXNOD,MNNOD
    !    2,KSNC,RCCC,FRCC,RCSC,SNCT,WGNODE(K,NB,NZ,NY,NX)
    !    3,HTNODX(K,NB,NZ,NY,NX),WTSTKB(NB,NZ,NY,NX)
    !    4,CPOOL(NB,NZ,NY,NX)
    !     ENDIF
    !
    !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !
          IF(WGNODE(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            RCCK=RCSC*WGNODE(K,NB,NZ,NY,NX)
            RCZK=WGNODN(K,NB,NZ,NY,NX)*(RCSN+(1.0-RCSN)*RCSC)
            RCPK=WGNODP(K,NB,NZ,NY,NX)*(RCSP+(1.0-RCSP)*RCSC)
    !
      !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
      !
      !     FSNCS=fraction of lowest internode to be remobilized
!
            IF(RCCK.GT.ZEROP(NZ,NY,NX))THEN
              FSNCK=AMAX1(0.0,AMIN1(1.0,SNCT/RCCK))
            ELSE
              FSNCK=1.0_r8
            ENDIF
    !
      !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
      !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
      !     FSNCK=fraction of lowest internode to be remobilized
      !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
      !     WGNODE,WGNODN,WGNODP=senescing internode C,N,P mass
!
            DO 7310 M=1,4
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
                *FSNCK*(WGNODE(K,NB,NZ,NY,NX)-RCCK)*FWOOD(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
                *FSNCK*(WGNODN(K,NB,NZ,NY,NX)-RCZK)*FWOODN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
                *FSNCK*(WGNODP(K,NB,NZ,NY,NX)-RCPK)*FWOODP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
                *FSNCK*(WGNODE(K,NB,NZ,NY,NX)-RCCK)*FWOOD(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
                *FSNCK*(WGNODN(K,NB,NZ,NY,NX)-RCZK)*FWOODN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
                *FSNCK*(WGNODP(K,NB,NZ,NY,NX)-RCPK)*FWOODP(1)
7310        CONTINUE
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
      !
      !     FSNCK=fraction of lowest internode to be remobilized
      !     HTNODE,HTNODX=living,senescing internode length
      !     WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
      !
            WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX) &
              -FSNCK*WGNODE(K,NB,NZ,NY,NX))
            WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX) &
              -FSNCK*WGNODN(K,NB,NZ,NY,NX))
            WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX) &
              -FSNCK*WGNODP(K,NB,NZ,NY,NX))
            HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX) &
              -FSNCK*HTNODX(K,NB,NZ,NY,NX)
            WGNODE(K,NB,NZ,NY,NX)=WGNODE(K,NB,NZ,NY,NX) &
              -FSNCK*WGNODE(K,NB,NZ,NY,NX)
            WGNODN(K,NB,NZ,NY,NX)=WGNODN(K,NB,NZ,NY,NX) &
              -FSNCK*WGNODN(K,NB,NZ,NY,NX)
            WGNODP(K,NB,NZ,NY,NX)=WGNODP(K,NB,NZ,NY,NX) &
              -FSNCK*WGNODP(K,NB,NZ,NY,NX)
            HTNODX(K,NB,NZ,NY,NX)=HTNODX(K,NB,NZ,NY,NX) &
              -FSNCK*HTNODX(K,NB,NZ,NY,NX)
!
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
      !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
      !     FSNCK=fraction of lowest internode to be remobilized
      !     SNCT=remaining node senescence respiration
      !
            WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+FSNCK*RCCK*SNCF
            WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+FSNCK*RCZK
            WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+FSNCK*RCPK
            SNCT=SNCT-FSNCK*RCCK
      !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
      !     WRITE(*,2356)'WGNODE2',I,J,NX,NY,NZ,NB,K,KK,MXNOD,MNNOD
      !    2,KSNC,RCCC,FRCC,RCSC,SNCT,WGNODE(K,NB,NZ,NY,NX)
      !    3,HTNODX(K,NB,NZ,NY,NX),WTSTKB(NB,NZ,NY,NX)
      !    4,CPOOL(NB,NZ,NY,NX)
      !2356  FORMAT(A8,11I4,12E16.8)
      !     ENDIF
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
            IF(SNCT.LE.ZEROP(NZ,NY,NX))GO TO 565
    !
    !       OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
    !
    !       CSNC,ZSNC,PSNC=literfall C,N,P
    !       CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !       WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in senescing internode
    !       HTNODE,HTNODX=living,senescing internode length
    !
          ELSE
            DO 7315 M=1,4
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
                *WGNODE(K,NB,NZ,NY,NX)*FWOOD(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
                *WGNODN(K,NB,NZ,NY,NX)*FWOODN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
                *WGNODP(K,NB,NZ,NY,NX)*FWOODP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
                *WGNODE(K,NB,NZ,NY,NX)*FWOOD(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
                *WGNODN(K,NB,NZ,NY,NX)*FWOODN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
                *WGNODP(K,NB,NZ,NY,NX)*FWOODP(1)
7315        CONTINUE
            WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX) &
              -WGNODE(K,NB,NZ,NY,NX))
            WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX) &
              -WGNODN(K,NB,NZ,NY,NX))
            WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX) &
              -WGNODP(K,NB,NZ,NY,NX))
            HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX) &
              -HTNODX(K,NB,NZ,NY,NX)
            WGNODE(K,NB,NZ,NY,NX)=0._r8
            WGNODN(K,NB,NZ,NY,NX)=0._r8
            WGNODP(K,NB,NZ,NY,NX)=0._r8
            HTNODX(K,NB,NZ,NY,NX)=0._r8
          ENDIF
1650    CONTINUE
!
    !   RESIDUAL STALK
    !
    !   RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
        IF(WTSTXB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
          RCCK=RCSC*WTSTXB(NB,NZ,NY,NX)
          RCZK=WTSTXN(NB,NZ,NY,NX)*(RCSN+(1.0-RCSN)*RCSC)
          RCPK=WTSTXP(NB,NZ,NY,NX)*(RCSP+(1.0-RCSP)*RCSC)
    !
    !     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
          IF(RCCK.GT.ZEROP(NZ,NY,NX))THEN
            FSNCR=AMAX1(0.0,AMIN1(1.0,SNCT/RCCK))
          ELSE
            FSNCR=1.0_r8
          ENDIF
    !
    !     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
    !     TO FRACTIONS SET IN 'STARTQ'
    !
          DO 8310 M=1,4
            CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
              *FSNCR*(WTSTXB(NB,NZ,NY,NX)-RCCK)*FWOOD(0)
            ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
              *FSNCR*(WTSTXN(NB,NZ,NY,NX)-RCZK)*FWOODN(0)
            PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
              *FSNCR*(WTSTXP(NB,NZ,NY,NX)-RCPK)*FWOODP(0)
            CSNC(M,1,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
              *FSNCR*(WTSTXB(NB,NZ,NY,NX)-RCCK)*FWOOD(1)
            ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
              *FSNCR*(WTSTXN(NB,NZ,NY,NX)-RCZK)*FWOODN(1)
            PSNC(M,1,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
              *FSNCR*(WTSTXP(NB,NZ,NY,NX)-RCPK)*FWOODP(1)
8310      CONTINUE
    !
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !     HTNODE,HTNODX=living,senescing internode length
    !
          WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX) &
            -FSNCR*WTSTXB(NB,NZ,NY,NX))
          WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX) &
            -FSNCR*WTSTXN(NB,NZ,NY,NX))
          WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX) &
            -FSNCR*WTSTXP(NB,NZ,NY,NX))
          WTSTXB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTXB(NB,NZ,NY,NX) &
            -FSNCR*WTSTXB(NB,NZ,NY,NX))
          WTSTXN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTXN(NB,NZ,NY,NX) &
            -FSNCR*WTSTXN(NB,NZ,NY,NX))
          WTSTXP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTXP(NB,NZ,NY,NX) &
            -FSNCR*WTSTXP(NB,NZ,NY,NX))
          HTNODZ=0._r8
          DO 8320 K=0,25
            HTNODZ=AMAX1(HTNODZ,HTNODE(K,NB,NZ,NY,NX))
8320      CONTINUE
          HTNODZ=AMAX1(0.0,HTNODZ-FSNCR*HTNODZ)
          DO 8325 K=0,25
            HTNODE(K,NB,NZ,NY,NX)=AMIN1(HTNODZ,HTNODE(K,NB,NZ,NY,NX))
8325      CONTINUE
!
    !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
    !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
    !
    !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
    !     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
    !     FSNCR=fraction of residual stalk to be remobilized
    !     SNCT=remaining node senescence respiration
    !
          WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+FSNCR*RCCK*SNCF
          WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+FSNCR*RCZK
          WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+FSNCR*RCPK
          SNCT=SNCT-FSNCR*RCCK
        ENDIF
  !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
  !     WRITE(*,2357)'WTSTXB1',I,J,NZ,NB,K,FSNCR,SNCT
  !    3,WTSTKB(NB,NZ,NY,NX),WTSTXB(NB,NZ,NY,NX)
  !    4,(HTNODE(K,NB,NZ,NY,NX),K=0,25)
  !2357  FORMAT(A8,5I4,40E12.4)
  !     ENDIF
  !
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
        IF(SNCT.LE.ZEROP(NZ,NY,NX))GO TO 565
  !
  !     OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
  !
  !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
  !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
  !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
  !
      ELSE
        DO 8315 M=1,4
          CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
            *WTSTXB(NB,NZ,NY,NX)*FWOOD(0)
          ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
            *WTSTXN(NB,NZ,NY,NX)*FWOODN(0)
          PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
            *WTSTXP(NB,NZ,NY,NX)*FWOODP(0)
          CSNC(M,1,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX) &
            *WTSTXB(NB,NZ,NY,NX)*FWOOD(1)
          ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX) &
            *WTSTXN(NB,NZ,NY,NX)*FWOODN(1)
          PSNC(M,1,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX) &
            *WTSTXP(NB,NZ,NY,NX)*FWOODP(1)
8315    CONTINUE
        WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX)-WTSTXB(NB,NZ,NY,NX))
        WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX)-WTSTXN(NB,NZ,NY,NX))
        WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX)-WTSTXP(NB,NZ,NY,NX))
        WTSTXB(NB,NZ,NY,NX)=0._r8
        WTSTXN(NB,NZ,NY,NX)=0._r8
        WTSTXP(NB,NZ,NY,NX)=0._r8
      ENDIF
!
!     REMOBILIZATION OF STORAGE C,N,P
!
!     WTRVC=storage C
!     IDTHB=branch living flag: 0=alive,1=dead
!     SNCR=remaining excess maintenance respiration
!
      SNCR=SNCT/(1.0+FXFS)
      IF(WTRVC(NZ,NY,NX).GT.SNCR)THEN
        WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-SNCR
        SNCR=0._r8
        GO TO 565
      ELSEIF(ISTYP(NZ,NY,NX).NE.0)THEN
        IDTHB(NB,NZ,NY,NX)=1
      ENDIF
565   CONTINUE
575   CONTINUE
    ENDIF
595 CONTINUE
!
!   DEATH IF MAIN STALK OF TREE DIES
!
!   IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   IDTHB=branch living flag: 0=alive,1=dead
!   KVSTGX,KVSTG=integer of lowest,highest leaf number currently growing
!   WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   CNLFB,CPLFB=N:C,P:C ratios in leaf
!
    IF(IBTYP(NZ,NY,NX).NE.0.AND.IGTYP(NZ,NY,NX).GT.1 &
      .AND.IDTHB(NB1(NZ,NY,NX),NZ,NY,NX).EQ.1)IDTHB(NB,NZ,NY,NX)=1
!
!     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
!
    KVSTGX=MAX(0,KVSTG(NB,NZ,NY,NX)-24)
    DO 495 KK=KVSTGX,KVSTG(NB,NZ,NY,NX)
      K=MOD(KK,25)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      IF(WGLF(K,NB,NZ,NY,NX).GT.0.0)THEN
        CPOOLT=WGLF(K,NB,NZ,NY,NX)+CPOOL(NB,NZ,NY,NX)
        IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
          ZPOOLD=WGLFN(K,NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX) &
            -ZPOOL(NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)
          XFRN1=AMAX1(0.0,AMIN1(1.0E-03*ZPOOLD/CPOOLT,WGLFN(K,NB,NZ,NY,NX) &
            -ZPLFM*CNLFB*WGLF(K,NB,NZ,NY,NX)))
          PPOOLD=WGLFP(K,NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX) &
            -PPOOL(NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)
          XFRP1=AMAX1(0.0,AMIN1(1.0E-03*PPOOLD/CPOOLT,WGLFP(K,NB,NZ,NY,NX) &
            -ZPLFM*CPLFB*WGLF(K,NB,NZ,NY,NX)))
          XFRN=AMAX1(XFRN1,10.0*XFRP1)
          XFRP=AMAX1(XFRP1,0.10*XFRN1)
          WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)-XFRN
          WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)-XFRN
          ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+XFRN
          WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)-XFRP
          WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)-XFRP
          PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+XFRP
          WSLF(K,NB,NZ,NY,NX)=AMAX1(0.0,WSLF(K,NB,NZ,NY,NX) &
            -AMAX1(XFRN*CNWS(NZ,NY,NX),XFRP*CPWS(NZ,NY,NX)))
        ENDIF
      ENDIF
495 CONTINUE
!
    call AllocateLeafToCanopyLayers(NB,NZ,NY,NX)

!
    !     ALLOCATE LEAF AREA TO INCLINATION CLASSES ACCORDING TO
    !     DISTRIBUTION ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
    !
    !     SSIN=sine of solar angle
    !     SURF=leaf node surface area in canopy layer
    !     ARLF,ARLFL=leaf node surface area in canopy layer
    !     ZC,DPTHS=canopy,snowpack height
    !     CLASS=leaf inclination class
    !
    IF(SSINN(NY,NX).GT.0.0)THEN
      call LeafClassAllocation(NB,NZ,NY,NX)
    ENDIF

    call GrainFilling(I,NB,NZ,NY,NX)
!
    call PhenologyReset(I,NB,NZ,NY,NX)
!
    call CarbNutInBranchTransfer(I,J,NB,NZ,NY,NX)

!   CANOPY N2 FIXATION (CYANOBACTERIA)
!
    call CanopyNoduleBiochemistry(I,J,NZ,NY,NX)
  ENDIF
  end subroutine GrowOneBranch
!------------------------------------------------------------------------------------------

  subroutine PhenologyReset(I,NB,NZ,NY,NX)

  implicit none
  integer, intent(in) :: I,nb,nz,ny,nx

!   RESET PHENOLOGY AT EMERGENCE ('VRNS' > 'VRNL')
!   AND END OF SEASON ('VRNF' > 'VRNX')
!
!   ISTYP=growth habit:0=annual,1=perennial from PFT file
!   IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!   IFLGE=flag for enabling leafout:0=enable,1=disable
!   VRNS,VRNL=leafout hours,hours required for leafout
!   IFLGF=flag for enabling leafoff:0=enable,1=disable
!   VRNF,VRNX=leafoff hours,hours required for leafoff
!
  IF(ISTYP(NZ,NY,NX).NE.0.OR.(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0))THEN
    IF((IFLGE(NB,NZ,NY,NX).EQ.0.AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX)) &
      .OR.(IFLGF(NB,NZ,NY,NX).EQ.0.AND.VRNF(NB,NZ,NY,NX).GE.VRNX(NB,NZ,NY,NX)))THEN
!
  !    SPRING PHENOLOGY RESET
  !
  !    GROUP,GROUPI=node number required for floral initiation
  !    PSTGI=node number at floral initiation
  !    PSTGF=node number at flowering
  !    VSTGX=leaf number on date of floral initiation
  !    TGSTGI=total change in vegve node number normalized for maturity group
  !    TGSTGF=total change in reprve node number normalized for maturity group
  !    IDAY(1,=emergence date
!
      IF((IFLGE(NB,NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0) &
        .AND.(VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX)))THEN
        IF(ISTYP(NZ,NY,NX).EQ.0)THEN
          GROUP(NB,NZ,NY,NX)=AMAX1(0.0,GROUPI(NZ,NY,NX)-NBTB(NB,NZ,NY,NX))
        ELSE
          GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
        ENDIF
        PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
        PSTGF(NB,NZ,NY,NX)=0._r8
        VSTGX(NB,NZ,NY,NX)=0._r8
        TGSTGI(NB,NZ,NY,NX)=0._r8
        TGSTGF(NB,NZ,NY,NX)=0._r8
        IDAY(1,NB,NZ,NY,NX)=I
        DO 2005 M=2,10
          IDAY(M,NB,NZ,NY,NX)=0
2005    CONTINUE
        IF(NB.EQ.NB1(NZ,NY,NX))THEN
          WSTR(NZ,NY,NX)=0._r8
        ENDIF
    !
    !   SPRING LEAF AND SHEATH RESET
    !
    !   IFLGA,IFLGE=flag for initializing,enabling leafout
    !   VRNS,VRNL=leafout hours,hours required for leafout
    !   PSTG=node number
    !   VSTG=number of leaves appeared
    !     KVSTG=integer of most recent leaf number currently growing
    !     FLG4=number of hours with no grain fill
    !     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
    !     WT*,WT*N,WT*P=branch organ C,N,P mass
    !     WG,WG*N,WG*P=node organ C,N,P mass
    !     organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
    !     HSK=husk,EAR=ear,GR=grain,SHT=shoot
    !     ARLFB,ARLF=branch,node leaf area
    !     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
    !     HTNODX,HTNODE=stalk height,stalk internode length
    !     GRNOB=seed set number
    !     GRNXB=potential number of seed set sites
    !     GRWTB=individual seed size
!
        IF(IFLGE(NB,NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0 &
          .AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
          IF(IBTYP(NZ,NY,NX).EQ.0)THEN
            PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
            VSTG(NB,NZ,NY,NX)=0._r8
            KLEAF(NB,NZ,NY,NX)=1
            KVSTG(NB,NZ,NY,NX)=1
            FLG4(NB,NZ,NY,NX)=0._r8
            DO 5330 M=1,4
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX) &
                +CFOPC(5,M,NZ,NY,NX)*WTLFB(NB,NZ,NY,NX)*FWODB(0) &
                +CFOPC(5,M,NZ,NY,NX)*WTSHEB(NB,NZ,NY,NX)*FWODB(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX) &
                +CFOPN(5,M,NZ,NY,NX)*WTLFBN(NB,NZ,NY,NX)*FWODLN(0) &
                +CFOPN(5,M,NZ,NY,NX)*WTSHBN(NB,NZ,NY,NX)*FWODSN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX) &
                +CFOPP(5,M,NZ,NY,NX)*WTLFBP(NB,NZ,NY,NX)*FWODLP(0) &
                +CFOPP(5,M,NZ,NY,NX)*WTSHBP(NB,NZ,NY,NX)*FWODSP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
                +CFOPC(1,M,NZ,NY,NX)*WTLFB(NB,NZ,NY,NX)*FWODB(1) &
                +CFOPC(2,M,NZ,NY,NX)*WTSHEB(NB,NZ,NY,NX)*FWODB(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
                +CFOPN(1,M,NZ,NY,NX)*WTLFBN(NB,NZ,NY,NX)*FWODLN(1) &
                +CFOPN(2,M,NZ,NY,NX)*WTSHBN(NB,NZ,NY,NX)*FWODSN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
                +CFOPP(1,M,NZ,NY,NX)*WTLFBP(NB,NZ,NY,NX)*FWODLP(1) &
                +CFOPP(2,M,NZ,NY,NX)*WTSHBP(NB,NZ,NY,NX)*FWODSP(1)
5330        CONTINUE
            ARLFB(NB,NZ,NY,NX)=0._r8
            WTLFB(NB,NZ,NY,NX)=0._r8
            WTLFBN(NB,NZ,NY,NX)=0._r8
            WTLFBP(NB,NZ,NY,NX)=0._r8
            WTSHEB(NB,NZ,NY,NX)=0._r8
            WTSHBN(NB,NZ,NY,NX)=0._r8
            WTSHBP(NB,NZ,NY,NX)=0._r8
            DO 5335 K=0,25
              ARLF(K,NB,NZ,NY,NX)=0._r8
              HTSHE(K,NB,NZ,NY,NX)=0._r8
              WGLF(K,NB,NZ,NY,NX)=0._r8
              WSLF(K,NB,NZ,NY,NX)=0._r8
              WGLFN(K,NB,NZ,NY,NX)=0._r8
              WGLFP(K,NB,NZ,NY,NX)=0._r8
              WGSHE(K,NB,NZ,NY,NX)=0._r8
              WSSHE(K,NB,NZ,NY,NX)=0._r8
              WGSHN(K,NB,NZ,NY,NX)=0._r8
              WGSHP(K,NB,NZ,NY,NX)=0._r8
5335        CONTINUE
          ENDIF
        ENDIF
    !
    !     RESIDUAL STALKS BECOME LITTERFALL IN GRASSES, SHRUBS AT
    !     START OF SEASON
    !
        IF((IFLGE(NB,NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0) &
          .AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
          DO 6245 M=1,4
            CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX) &
              *(WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX))
            ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX) &
              *(WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX))
            PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX) &
              *(WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX))
6245      CONTINUE
          WTHSKB(NB,NZ,NY,NX)=0._r8
          WTEARB(NB,NZ,NY,NX)=0._r8
          WTGRB(NB,NZ,NY,NX)=0._r8
          WTHSBN(NB,NZ,NY,NX)=0._r8
          WTEABN(NB,NZ,NY,NX)=0._r8
          WTGRBN(NB,NZ,NY,NX)=0._r8
          WTHSBP(NB,NZ,NY,NX)=0._r8
          WTEABP(NB,NZ,NY,NX)=0._r8
          WTGRBP(NB,NZ,NY,NX)=0._r8
          GRNXB(NB,NZ,NY,NX)=0._r8
          GRNOB(NB,NZ,NY,NX)=0._r8
          GRWTB(NB,NZ,NY,NX)=0._r8
          IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
            DO 6345 M=1,4
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)
6345        CONTINUE
            WTSTKB(NB,NZ,NY,NX)=0._r8
            WTSTBN(NB,NZ,NY,NX)=0._r8
            WTSTBP(NB,NZ,NY,NX)=0._r8
            WTSTXB(NB,NZ,NY,NX)=0._r8
            WTSTXN(NB,NZ,NY,NX)=0._r8
            WTSTXP(NB,NZ,NY,NX)=0._r8
            DO 6340 K=0,25
              HTNODE(K,NB,NZ,NY,NX)=0._r8
              HTNODX(K,NB,NZ,NY,NX)=0._r8
              WGNODE(K,NB,NZ,NY,NX)=0._r8
              WGNODN(K,NB,NZ,NY,NX)=0._r8
              WGNODP(K,NB,NZ,NY,NX)=0._r8
6340        CONTINUE
          ENDIF
        ENDIF
      ENDIF
!
  !   SPRING OR FALL FLAG RESET
  !
      IF(IFLGE(NB,NZ,NY,NX).EQ.0 &
        .AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
        IFLGE(NB,NZ,NY,NX)=1
        IFLGF(NB,NZ,NY,NX)=0
        IFLGR(NB,NZ,NY,NX)=0
        IFLGQ(NB,NZ,NY,NX)=0
      ELSE
        IFLGE(NB,NZ,NY,NX)=0
        IFLGF(NB,NZ,NY,NX)=1
        IFLGR(NB,NZ,NY,NX)=1
        IFLGQ(NB,NZ,NY,NX)=0
        IFLGA(NB,NZ,NY,NX)=0
      ENDIF
    ENDIF
  ENDIF
!
!   REPRODUCTIVE MATERIAL BECOMES LITTERFALL AT END OF SEASON
!
  IF(IFLGR(NB,NZ,NY,NX).EQ.1)THEN
    IFLGQ(NB,NZ,NY,NX)=IFLGQ(NB,NZ,NY,NX)+1
    IF(IFLGQ(NB,NZ,NY,NX).EQ.IFLGQX)THEN
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
    ENDIF
    DO 6330 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+FSNR*CFOPC(2,M,NZ,NY,NX) &
        *(WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+FSNR*CFOPN(2,M,NZ,NY,NX) &
        *(WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+FSNR*CFOPP(2,M,NZ,NY,NX) &
        *(WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX))
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
        WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX) &
          +FSNR*CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
        WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX) &
          +FSNR*CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
        WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX) &
          +FSNR*CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
        CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
          +FSNR*CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
        ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
          +FSNR*CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
        PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
          +FSNR*CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
6330  CONTINUE
    WTHSKB(NB,NZ,NY,NX)=(1.0-FSNR)*WTHSKB(NB,NZ,NY,NX)
    WTEARB(NB,NZ,NY,NX)=(1.0-FSNR)*WTEARB(NB,NZ,NY,NX)
    WTGRB(NB,NZ,NY,NX)=(1.0-FSNR)*WTGRB(NB,NZ,NY,NX)
    WTHSBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTHSBN(NB,NZ,NY,NX)
    WTEABN(NB,NZ,NY,NX)=(1.0-FSNR)*WTEABN(NB,NZ,NY,NX)
    WTGRBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTGRBN(NB,NZ,NY,NX)
    WTHSBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTHSBP(NB,NZ,NY,NX)
    WTEABP(NB,NZ,NY,NX)=(1.0-FSNR)*WTEABP(NB,NZ,NY,NX)
    WTGRBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTGRBP(NB,NZ,NY,NX)
    GRNXB(NB,NZ,NY,NX)=(1.0-FSNR)*GRNXB(NB,NZ,NY,NX)
    GRNOB(NB,NZ,NY,NX)=(1.0-FSNR)*GRNOB(NB,NZ,NY,NX)
    GRWTB(NB,NZ,NY,NX)=(1.0-FSNR)*GRWTB(NB,NZ,NY,NX)
!
!     STALKS BECOME LITTERFALL IN GRASSES AT END OF SEASON
!
    IF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1).AND.ISTYP(NZ,NY,NX).NE.0)THEN
      DO 6335 M=1,4
        CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+FSNR*CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)
        ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+FSNR*CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)
        PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+FSNR*CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)
6335    CONTINUE
      WTSTKB(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTKB(NB,NZ,NY,NX)
      WTSTBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTBN(NB,NZ,NY,NX)
      WTSTBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTBP(NB,NZ,NY,NX)
      WTSTXB(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTXB(NB,NZ,NY,NX)
      WTSTXN(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTXN(NB,NZ,NY,NX)
      WTSTXP(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTXP(NB,NZ,NY,NX)
      DO 2010 K=0,25
    !     HTNODE(K,NB,NZ,NY,NX)=(1.0-FSNR)*HTNODE(K,NB,NZ,NY,NX)
        HTNODX(K,NB,NZ,NY,NX)=(1.0-FSNR)*HTNODX(K,NB,NZ,NY,NX)
        WGNODE(K,NB,NZ,NY,NX)=(1.0-FSNR)*WGNODE(K,NB,NZ,NY,NX)
        WGNODN(K,NB,NZ,NY,NX)=(1.0-FSNR)*WGNODN(K,NB,NZ,NY,NX)
        WGNODP(K,NB,NZ,NY,NX)=(1.0-FSNR)*WGNODP(K,NB,NZ,NY,NX)
2010  CONTINUE
    ENDIF
!
!     SELF-SEEDING ANNUALS IF COLD OR DROUGHT DECIDUOUS
!
!     ISTYP=growth habit:0=annual,1=perennial
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IDAYH,IYRH=day,year of harvesting
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     IDAY0,IYR0=day,year of planting
!     IFLGI=PFT initialization flag:0=no,1=yes
!
!     IF(J.EQ.INT(ZNOON(NY,NX)))THEN
    IF(NB.EQ.NB1(NZ,NY,NX))THEN
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
        IDAYH(NZ,NY,NX)=I
        IYRH(NZ,NY,NX)=IYRC
        IHVST(NZ,I,NY,NX)=1
        JHVST(NZ,I,NY,NX)=2
        HVST(NZ,I,NY,NX)=0._r8
        THIN(NZ,I,NY,NX)=0._r8
        EHVST(1,1,NZ,I,NY,NX)=1.0_r8
        EHVST(1,2,NZ,I,NY,NX)=1.0_r8
        EHVST(1,3,NZ,I,NY,NX)=1.0_r8
        EHVST(1,4,NZ,I,NY,NX)=1.0_r8
        EHVST(2,1,NZ,I,NY,NX)=0._r8
        EHVST(2,2,NZ,I,NY,NX)=1.0_r8
        EHVST(2,3,NZ,I,NY,NX)=0._r8
        EHVST(2,4,NZ,I,NY,NX)=0._r8
        IDAY0(NZ,NY,NX)=-1E+06
        IYR0(NZ,NY,NX)=-1E+06
        IFLGI(NZ,NY,NX)=1
!     WRITE(*,3366)'HVST',I,J,IYRC,IDAYH(NZ,NY,NX),IYRH(NZ,NY,NX)
!    2,IHVST(NZ,I,NY,NX),JHVST(NZ,I,NY,NX),IFLGI(NZ,NY,NX)
!3366  FORMAT(A8,8I8)
      ENDIF
    ENDIF
!     ENDIF
  ENDIF
  end subroutine PhenologyReset

!------------------------------------------------------------------------------------------

  subroutine AllocateLeafToCanopyLayers(NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NB,NZ,NY,NX


! begin_execution
!   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
!
!   HTCTL=hypocotyledon height
!   SDPTH=seeding depth
!   ARLF=node leaf area
!   HTSHE=petiole length
!
  KVSTGN(NB,NZ,NY,NX)=0
  IF(HTCTL(NZ,NY,NX).LE.SDPTH(NZ,NY,NX) &
    .AND.ARLF(0,NB1(NZ,NY,NX),NZ,NY,NX).GT.0.0)THEN
    XLGLF=SQRT(1.0E+02*ARLF(0,NB1(NZ,NY,NX),NZ,NY,NX)/PP(NZ,NY,NX))
    HTCTL(NZ,NY,NX)=XLGLF+HTSHE(0,NB1(NZ,NY,NX),NZ,NY,NX) &
      +HTNODE(0,NB1(NZ,NY,NX),NZ,NY,NX)
  ENDIF
!
! IF CANOPY HAS EMERGED
!
  IF(HTCTL(NZ,NY,NX).GT.SDPTH(NZ,NY,NX))THEN
    DO 540 K=0,25
      DO  L=1,JC
        ARLFL(L,K,NB,NZ,NY,NX)=0._r8
        WGLFL(L,K,NB,NZ,NY,NX)=0._r8
        WGLFLN(L,K,NB,NZ,NY,NX)=0._r8
        WGLFLP(L,K,NB,NZ,NY,NX)=0._r8
      enddo
540   CONTINUE
    DO 535 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0._r8
535 CONTINUE
!
!   BRANCH HEIGHT
!
!   IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!   IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   KVSTG1,KVSTGN=integer of highest,lowest leaf number currently growing
!   HTNODE=internode length
!   HTBR=branch base height
!
    IF(IBTYP(NZ,NY,NX).NE.0.AND.IGTYP(NZ,NY,NX).GT.1)THEN
      IF(NB.NE.NB1(NZ,NY,NX))THEN
        KVSTG1=MAX(1,KVSTG(NB1(NZ,NY,NX),NZ,NY,NX)-24)
        IF(NBTB(NB,NZ,NY,NX).GE.KVSTG1)THEN
        K=MOD(NBTB(NB,NZ,NY,NX),25)
        IF(K.EQ.0.AND.KK.NE.0)K=25
          HTBR=HTNODE(K,NB1(NZ,NY,NX),NZ,NY,NX)
        ELSE
          HTBR=0._r8
        ENDIF
      ELSE
        HTBR=0._r8
      ENDIF
    ELSE
      HTBR=0._r8
    ENDIF
    KVSTGX=MAX(0,KVSTG(NB,NZ,NY,NX)-24)
!
!   FOR ALL LEAFED NODES
!
    DO 560 KK=KVSTGX,KVSTG(NB,NZ,NY,NX)
      K=MOD(KK,25)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      !
      !     HEIGHT OF STALK INTERNODE + SHEATH OR PETIOLE
      !     AND LENGTH OF LEAF
      !
      !     HTSTK=stalk height
      !     HTNODE=internode length
      !     HTLF=leaf node height
      !     ARLF=leaf node area
      !     PP=plant population
      !     FNOD=scales node number for perennial vegetation (e.g. trees)
      !     XLGLF=leaf length
!
      HTSTK=HTBR+HTNODE(K,NB,NZ,NY,NX)
      HTLF=HTSTK+HTSHE(K,NB,NZ,NY,NX)
      XLGLF=AMAX1(0.0,SQRT(WDLF(NZ,NY,NX)*AMAX1(0.0 &
        ,ARLF(K,NB,NZ,NY,NX))/(PP(NZ,NY,NX)*FNOD(NZ,NY,NX))))
      TLGLF=0._r8
      !
      !   ALLOCATE FRACTIONS OF LEAF IN EACH INCLINATION CLASS
      !     FROM HIGHEST TO LOWEST TO CANOPY LAYER
      !
      !     YLGLF=leaf elevation
      !     CLASS=leaf inclination class
      !     XLGLF=leaf length
      !     ZC,ZCX=canopy height
      !     HTLFL,HTLFU=height of leaf base,tip
      !     ZL=height to bottom of each canopy layer
      !     LHTLFL,LHTLFU=layer number of leaf base,tip
      !     FRACL=leaf fraction in each layer
!
      DO 555 N=4,1,-1
        YLGLF=ZSIN(N)*CLASS(N,NZ,NY,NX)*XLGLF
        HTLFL=AMIN1(ZCX(NZ,NY,NX)+0.01-YLGLF,HTLF+TLGLF)
        HTLFU=AMIN1(ZCX(NZ,NY,NX)+0.01,HTLFL+YLGLF)
        LU=0
        LL=0
        DO 550 L=JC,1,-1
          IF(LU.EQ.1.AND.LL.EQ.1)GO TO 551
          IF((HTLFU.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO).AND.LU.EQ.0)THEN
            LHTLFU=MAX(1,L)
            LU=1
          ENDIF
          IF((HTLFL.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO).AND.LL.EQ.0)THEN
            LHTLFL=MAX(1,L)
            LL=1
          ENDIF
550     CONTINUE
551     CONTINUE
        DO 570 L=LHTLFL,LHTLFU
          IF(LHTLFU.EQ.LHTLFL)THEN
            FRACL=CLASS(N,NZ,NY,NX)
          ELSEIF(HTLFU.GT.HTLFL.AND.ZL(L,NY,NX).GT.HTLFL)THEN
            FRACL=CLASS(N,NZ,NY,NX)*(AMIN1(HTLFU,ZL(L,NY,NX)) &
              -AMAX1(HTLFL,ZL(L-1,NY,NX)))/(HTLFU-HTLFL)
          ELSE
            FRACL=CLASS(N,NZ,NY,NX)
          ENDIF
          YARLF=FRACL*ARLF(K,NB,NZ,NY,NX)
          YWGLF=FRACL*WGLF(K,NB,NZ,NY,NX)
          YWGLFN=FRACL*WGLFN(K,NB,NZ,NY,NX)
          YWGLFP=FRACL*WGLFP(K,NB,NZ,NY,NX)
!
    !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
    !
    !     ARLFL=leaf node area in canopy layer
    !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
    !     ARLFV,WGLFV=total leaf area,C in canopy layer
    !     HTNODE=internode length
    !
          ARLFL(L,K,NB,NZ,NY,NX)=ARLFL(L,K,NB,NZ,NY,NX)+YARLF
          WGLFL(L,K,NB,NZ,NY,NX)=WGLFL(L,K,NB,NZ,NY,NX)+YWGLF
          WGLFLN(L,K,NB,NZ,NY,NX)=WGLFLN(L,K,NB,NZ,NY,NX)+YWGLFN
          WGLFLP(L,K,NB,NZ,NY,NX)=WGLFLP(L,K,NB,NZ,NY,NX)+YWGLFP
          ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)+YARLF
          WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)+YWGLF
          !     IF(NZ.EQ.2)THEN
          !     WRITE(*,4813)'GRO',I,J,NX,NY,NZ,NB,K,KK,L,LHTLFL,LHTLFU
          !    2,FRACL,HTLFU,HTLFL,ZL(L-1,NY,NX),ARLFB(NB,NZ,NY,NX)
          !    3,ARLF(K,NB,NZ,NY,NX),WTLFB(NB,NZ,NY,NX),WGLF(K,NB,NZ,NY,NX)
          !    4,ARLFP(NZ,NY,NX),ZL(L,NY,NX),HTLF,TLGLF,HTSTK,HTBR
          !    4,HTNODE(K,NB,NZ,NY,NX),HTSHE(K,NB,NZ,NY,NX),YLGLF
          !    5,ZSIN(N),CLASS(N,NZ,NY,NX),XLGLF,ZC(NZ,NY,NX)
          !    6,ZCX(NZ,NY,NX)
          !4813  FORMAT(A8,11I4,30E12.4)
          !     ENDIF
570     CONTINUE
        TLGLF=TLGLF+YLGLF
        ZC(NZ,NY,NX)=AMAX1(ZC(NZ,NY,NX),HTLFU)
555   CONTINUE
      IF(WSSHE(K,NB,NZ,NY,NX).GT.0.0)THEN
        IF(KVSTGN(NB,NZ,NY,NX).EQ.0)KVSTGN(NB,NZ,NY,NX)=MIN(KK,KVSTG(NB,NZ,NY,NX))
      ENDIF
560 CONTINUE
    IF(KVSTGN(NB,NZ,NY,NX).EQ.0)KVSTGN(NB,NZ,NY,NX)=KVSTG(NB,NZ,NY,NX)
    K1=MOD(KVSTG(NB,NZ,NY,NX),25)
    IF(K1.EQ.0.AND.KVSTG(NB,NZ,NY,NX).NE.0)K1=25
    K2=MOD(KVSTG(NB,NZ,NY,NX)-1,25)
    IF(K2.EQ.0.AND.KVSTG(NB,NZ,NY,NX)-1.NE.0)K2=25
    IF(test_aeqb(HTNODE(K1,NB,NZ,NY,NX),0._r8))THEN
      HTNODE(K1,NB,NZ,NY,NX)=HTNODE(K2,NB,NZ,NY,NX)
    ENDIF
    HTLFB=HTBR+AMAX1(0.0,HTNODE(K1,NB,NZ,NY,NX))
!
  !     ALLOCATE STALK SURFACE AREA TO CANOPY LAYERS
  !
  !     HTNODE=internode length
  !     HTLFB=leaf base height
  !     ZL=height to bottom of each canopy layer
  !     LHTBRL,LHTBRU=layer number of branch base,tip
  !     WTSTKB,ARSTKB=branch stalk mass,surface area
  !     FSTK=fraction of stalk area contributing to water,heat flow
  !     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
  !     WVSTKB=stalk sapwood mass
  !     FRACL=stalk fraction in each layer
  !     ARSTK=total branch stalk surface area in each layer
  !
  !     IF(NZ.EQ.1)THEN
  !     WRITE(*,6679)'K1',I,J,NZ,NB,K1,KVSTG(NB,NZ,NY,NX)
  !    2,HTNODE(K1,NB,NZ,NY,NX)
  !6679  FORMAT(A8,6I4,12E12.4)
  !     ENDIF
    IF(HTNODE(K1,NB,NZ,NY,NX).GT.0.0)THEN
      LU=0
      LL=0
      DO 545 L=JC,1,-1
        IF(LU.EQ.1.AND.LL.EQ.1)GO TO 546
        IF((HTLFB.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO).AND.LU.EQ.0)THEN
          LHTBRU=MAX(1,L)
          LU=1
        ENDIF
        IF((HTBR.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO).AND.LL.EQ.0)THEN
          LHTBRL=MAX(1,L)
          LL=1
        ENDIF
545   CONTINUE
546   CONTINUE
      RSTK=SQRT(VSTK*(AMAX1(0.0_r8,WTSTKB(NB,NZ,NY,NX))/PP(NZ,NY,NX)) &
        /(3.1416_r8*HTNODE(K1,NB,NZ,NY,NX)))
      ARSTKB(NB)=3.1416_r8*HTNODE(K1,NB,NZ,NY,NX)*PP(NZ,NY,NX)*RSTK
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
        WVSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)
      ELSE
        ZSTK=AMIN1(ZSTX,FSTK*RSTK)
        ASTV=3.1416_r8*(2.0_r8*RSTK*ZSTK-ZSTK**2)
        WVSTKB(NB,NZ,NY,NX)=ASTV/VSTK*HTNODE(K1,NB,NZ,NY,NX)*PP(NZ,NY,NX)
      ENDIF
    !     IF(NZ.EQ.1)THEN
        !     WRITE(*,6677)'WVSTK',I,J,NX,NY,NZ,NB,K1,WVSTKB(NB,NZ,NY,NX)
        !    2,ASTV,VSTK,HTNODE(K1,NB,NZ,NY,NX),PP(NZ,NY,NX)
        !6677  FORMAT(A8,7I4,12E12.4)
        !     ENDIF
      DO 445 L=LHTBRL,LHTBRU
        IF(HTLFB.GT.HTBR)THEN
          IF(HTLFB.GT.ZL(L-1,NY,NX))THEN
            FRACL=(AMIN1(HTLFB,ZL(L,NY,NX))-AMAX1(HTBR &
              ,ZL(L-1,NY,NX)))/(HTLFB-HTBR)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        ARSTK(L,NB,NZ,NY,NX)=FRACL*ARSTKB(NB)
445   CONTINUE
    ELSE
      WVSTKB(NB,NZ,NY,NX)=0._r8
      DO 450 L=1,JC
        ARSTK(L,NB,NZ,NY,NX)=0._r8
450   CONTINUE
    ENDIF
  ELSE
    WVSTKB(NB,NZ,NY,NX)=0._r8
    DO 455 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0._r8
455 CONTINUE
  ENDIF
end subroutine AllocateLeafToCanopyLayers

!------------------------------------------------------------------------------------------

  subroutine GrainFilling(I,NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,NB,NZ,NY,NX

! begin_execution

!
!   SET MAXIMUM GRAIN NUMBER FROM SHOOT MASS BEFORE ANTHESIS
!
!   IDAY(3,=start of stem elongation and setting max seed number
!   IDAY(6,=start of anthesis and setting final seed number
!   GRNXB=potential number of seed set sites
!   STMX=maximum potential seed number from PFT file
!     GROSTK=stalk growth rate
!
  IF(IDAY(3,NB,NZ,NY,NX).NE.0.AND.IDAY(6,NB,NZ,NY,NX).EQ.0)THEN
    GRNXB(NB,NZ,NY,NX)=GRNXB(NB,NZ,NY,NX)+STMX(NZ,NY,NX)*AMAX1(0.0,GROSTK)
!     WRITE(*,4246)'GRNX',I,J,NZ,NB,IDAY(3,NB,NZ,NY,NX)
!    2,GRNXB(NB,NZ,NY,NX),STMX(NZ,NY,NX),CGROS,GROSTK
  ENDIF
!
!   SET FINAL GRAIN NUMBER FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!   IDAY(6,=start of anthesis and setting final seed number
!   IDAY(7,=start of grain filling and setting max seed size
!   IDAY(8,=end date setting for final seed number
!   IDAY(9,=end of setting max seed size
!   SET=seed set limited by nonstructural C,N,P
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   TCC=canopy temperature
!   CTC=chilling temperature for CO2 fixation, seed loss (oC)
!   HTC=high temperature threshold for grain number loss
!   FGRNX=loss of seed set
!   SSTX=sensitivity to TCC > HTC,TCC < CTC from startq.f (seeds oC-1)
!   GRNOB=seed set number
!   GRNXB=potential number of seed set sites
!   SDMX=maximum seed number per STMX from PFT file
!   DGSTGF=change in reproductive node number normalized for maturity group
!
  IF(IDAY(6,NB,NZ,NY,NX).NE.0.AND.IDAY(9,NB,NZ,NY,NX).EQ.0)THEN
    SET=AMIN1(CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX)+SETC) &
      ,CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX)+SETN) &
      ,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX)+SETP))
    IF(TCC(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
      IF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
        FGRNX=SSTX(NZ,NY,NX)*(CTC(NZ,NY,NX)-TCC(NZ,NY,NX))
      ELSEIF(IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
        FGRNX=SSTX(NZ,NY,NX)*(CTC(NZ,NY,NX)-TCC(NZ,NY,NX))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSEIF(TCC(NZ,NY,NX).GT.HTC(NZ,NY,NX))THEN
      IF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
        FGRNX=SSTX(NZ,NY,NX)*(TCC(NZ,NY,NX)-HTC(NZ,NY,NX))
      ELSEIF(IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
        FGRNX=SSTX(NZ,NY,NX)*(TCC(NZ,NY,NX)-HTC(NZ,NY,NX))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSE
      FGRNX=0._r8
    ENDIF
    IF(IDAY(6,NB,NZ,NY,NX).NE.0.AND.IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
!     GRNXB(NB,NZ,NY,NX)=GRNXB(NB,NZ,NY,NX)*FGRNX
      GRNOB(NB,NZ,NY,NX)=AMIN1(SDMX(NZ,NY,NX)*GRNXB(NB,NZ,NY,NX) &
        ,GRNOB(NB,NZ,NY,NX)+SDMX(NZ,NY,NX)*GRNXB(NB,NZ,NY,NX) &
        *SET*DGSTGF(NB,NZ,NY,NX)-FGRNX*GRNOB(NB,NZ,NY,NX))
!     IF(FGRNX.LT.1.0)THEN
!     WRITE(*,4246)'GRNO',I,J,NZ,NB,IDAY(7,NB,NZ,NY,NX),TCC(NZ,NY,NX)
!    2,HTC(NZ,NY,NX),FGRNX,GRNXB(NB,NZ,NY,NX),GRNOB(NB,NZ,NY,NX)
!    3,SET,CCPOLB(NB,NZ,NY,NX),CZPOLB(NB,NZ,NY,NX)
!    4,CPPOLB(NB,NZ,NY,NX)
!4246  FORMAT(A8,5I4,20E12.4)
!     ENDIF
    ENDIF
!

!     SET MAXIMUM GRAIN SIZE FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!     GRMX=maximum individual seed size from PFT file (g)
!     DGSTGF=change in reproductive node number normalized for maturity group
!     GRWTB=individual seed size
!     SET=seed set limited by nonstructural C,N,P
!
    IF(IDAY(7,NB,NZ,NY,NX).NE.0.AND.IDAY(9,NB,NZ,NY,NX).EQ.0)THEN
      GRMXB=GRMX(NZ,NY,NX)
      GRWTB(NB,NZ,NY,NX)=AMIN1(GRMX(NZ,NY,NX),GRWTB(NB,NZ,NY,NX) &
        +GRMXB*AMAX1(0.50,SET**0.25)*DGSTGF(NB,NZ,NY,NX))
!       IF(FGRNX.LT.1.0)THEN
!       WRITE(*,4246)'GRWT',I,J,NZ,NB,IDAY(8,NB,NZ,NY,NX),TCC(NZ,NY,NX)
!      2,HTC(NZ,NY,NX),FGRNX,GRMX(NZ,NY,NX),GRWTB(NB,NZ,NY,NX)
!       ENDIF
    ENDIF
  ENDIF
!
!   GRAIN FILL BY TRANSLOCATION FROM STALK RESERVES
!   UNTIL GRAIN SINK (=FINAL GRAIN NUMBER X MAXIMUM
!   GRAIN SIZE) IS FILLED OR RESERVES ARE EXHAUSTED
!
!   IDAY(7,=start of grain filling and setting max seed size
!   WTGRB=total seed C mass
!   GRWTB=individual seed size
!   GRNOB=seed set number
!   GROLM=maximum grain fill rate
!   GFILL=grain filling rate at 25 oC from PFT file
!   TFN3=temperature function for canopy growth
!   TFN4=temperature function for root growth
!
  IF(IDAY(7,NB,NZ,NY,NX).NE.0)THEN
    IF(WTGRB(NB,NZ,NY,NX).GE.GRWTB(NB,NZ,NY,NX)*GRNOB(NB,NZ,NY,NX))THEN
      GROLM=0._r8
    ELSEIF(IRTYP(NZ,NY,NX).EQ.0)THEN
      GROLM=AMAX1(0.0,GFILL(NZ,NY,NX)*GRNOB(NB,NZ,NY,NX)*SQRT(TFN3(NZ,NY,NX)))
    ELSE
      GROLM=AMAX1(0.0,GFILL(NZ,NY,NX)*GRNOB(NB,NZ,NY,NX) &
        *SQRT(TFN4(NG(NZ,NY,NX),NZ,NY,NX)))
    ENDIF
!
!     GRAIN FILL RATE MAY BE CONSTRAINED BY HIGH GRAIN C:N OR C:P
!
!     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     GROLM,GROLC=maximum,actual grain fill rate
!     XLOCM,XLOCC=maximum,actual C translocation rate from reserve to grain
!
    IF(WTGRBN(NB,NZ,NY,NX).LT.ZPGRM*CNGR(NZ,NY,NX) &
      *WTGRB(NB,NZ,NY,NX).OR.WTGRBP(NB,NZ,NY,NX).LT.ZPGRM &
      *CPGR(NZ,NY,NX)*WTGRB(NB,NZ,NY,NX))THEN
      GROLC=0._r8
    ELSE
      GROLC=GROLM
    ENDIF
    XLOCM=AMIN1(GROLM,WTRSVB(NB,NZ,NY,NX))
    XLOCC=AMIN1(GROLC,WTRSVB(NB,NZ,NY,NX))
!
!     GRAIN N OR P FILL RATE MAY BE LIMITED BY C:N OR C:P RATIOS
!     OF STALK RESERVES
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     ZNPGN,ZNPGP=effect of reserve N:C,P:C on grain fill N:C,P:C
!     SETN,SETP=Km for nonstructural N,P concn on seed set (g g-1)
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     ZPGRD=1.0_r8-ZPGRM
!     ZPGRN,ZPGRP=N:C,P:C ratios during grain fill
!     XLOCM,XLOCC=maximum,actual C translocation rate from reserve to grain
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     XLOCN,XLOCP=N,P translocation rate from reserve to grain
!
    IF(WTRSVB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      ZNPGN=WTRSBN(NB,NZ,NY,NX)/(WTRSBN(NB,NZ,NY,NX) &
        +SETN*WTRSVB(NB,NZ,NY,NX))
      ZNPGP=WTRSBP(NB,NZ,NY,NX)/(WTRSBP(NB,NZ,NY,NX) &
        +SETP*WTRSVB(NB,NZ,NY,NX))
      ZPGRN=ZPGRM+ZPGRD*AMAX1(0.0,AMIN1(1.0,ZNPGN))
      ZPGRP=ZPGRM+ZPGRD*AMAX1(0.0,AMIN1(1.0,ZNPGP))
      XLOCN=AMIN1(XLOCM*CNGR(NZ,NY,NX) &
        ,AMAX1(0.0,WTRSBN(NB,NZ,NY,NX)*ZPGRN) &
        ,(WTGRB(NB,NZ,NY,NX)+XLOCC)*CNGR(NZ,NY,NX)-WTGRBN(NB,NZ,NY,NX))
      XLOCP=AMIN1(XLOCM*CPGR(NZ,NY,NX) &
        ,AMAX1(0.0,WTRSBP(NB,NZ,NY,NX)*ZPGRP) &
        ,(WTGRB(NB,NZ,NY,NX)+XLOCC)*CPGR(NZ,NY,NX)-WTGRBP(NB,NZ,NY,NX))
    ELSE
      XLOCN=0._r8
      XLOCP=0._r8
    ENDIF
!     IF(NX.EQ.1.AND.NY.EQ.6.AND.NZ.EQ.3)THEN
!     WRITE(*,85)'XLOC',I,J,NZ,NB,WTGRB(NB,NZ,NY,NX),WTGRBN(NB,NZ,NY,NX)
!    2,WTRSVB(NB,NZ,NY,NX),WTRSBN(NB,NZ,NY,NX),XLOCC,XLOCN,XLOCP,XLOCM
!    3,CNGR(NZ,NY,NX),ZPGRX,ZNPG,GROLC,GROLM,GROGR,GROGRN
!    3,XLOCM*CNGR(NZ,NY,NX),AMAX1(0.0,WTRSBN(NB,NZ,NY,NX)*ZPGRX)
!    4,(WTGRB(NB,NZ,NY,NX)+XLOCC)*CNGR(NZ,NY,NX)-WTGRBN(NB,NZ,NY,NX)
!    4,GRNOB(NB,NZ,NY,NX),GRWTB(NB,NZ,NY,NX),GFILL(NZ,NY,NX)
!    5,SQRT(TFN3(NZ,NY,NX)),FLG4(NB,NZ,NY,NX)
!85    FORMAT(A8,4I4,20E12.4)
!     ENDIF
!
!     TRANSLOCATE C,N,P FROM STALK RESERVES TO GRAIN
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     GROGR=grain growth rate
!     XLOCC,XLOCN,XLOCP=C,N,P translocation rate from reserve to grain
!
    WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+GROGR-XLOCC
    WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+GROGRN-XLOCN
    WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+GROGRP-XLOCP
    WTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)+XLOCC
    WTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)+XLOCN
    WTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)+XLOCP
  ELSE
    XLOCC=0._r8
    XLOCN=0._r8
    XLOCP=0._r8
  ENDIF
!
!   SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
!   HAS STOPPED FOR SET PERIOD OF TIME
!
!   IDAY(8,=end date setting for final seed number
!   XLOCC=C translocation rate from reserve to grain
!   PP=PFT population
!   FLG4=number of hours with no grain fill
!   FLG4X=number of hours with no grain filling until physl maturity
!   IDAY(10,=date of physiological maturity
!
  IF(IDAY(8,NB,NZ,NY,NX).NE.0)THEN
    IF(XLOCC.LE.1.0E-09*PP(NZ,NY,NX))THEN
      FLG4(NB,NZ,NY,NX)=FLG4(NB,NZ,NY,NX)+1.0
    ELSE
      FLG4(NB,NZ,NY,NX)=0._r8
    ENDIF
    IF(FLG4(NB,NZ,NY,NX).GE.FLG4X)THEN
      IF(IDAY(10,NB,NZ,NY,NX).EQ.0)THEN
        IDAY(10,NB,NZ,NY,NX)=I
      ENDIF
    ENDIF
!
!     TERMINATE ANNUALS AFTER GRAIN FILL
!
!     ISTYP=growth habit:0=annual,1=perennial
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     FLG4=number of hours with no grain fill
!     FLG4X=number of hours with no grain filling until physl maturity
!     FLG4Y=number of hours after physiol maturity required for senescence
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!
    IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      IF(FLG4(NB,NZ,NY,NX).GT.FLG4X+FLG4Y(IWTYP(NZ,NY,NX)))THEN
        VRNF(NB,NZ,NY,NX)=VRNX(NB,NZ,NY,NX)+0.5
      ENDIF
    ENDIF
  ENDIF
  end subroutine GrainFilling

!------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(nb,nz,ny,nx)
  implicit none
  integer, intent(in) :: NB,NZ,NY,NX

  ! begin_execution
  DO 900 K=1,25
    DO  L=1,JC
      DO  N=1,4
        SURF(N,L,K,NB,NZ,NY,NX)=0._r8
      enddo
    enddo
900 CONTINUE
! ARLFXB=0._r8
! ARLFXL=0._r8
! SURFXX=0._r8
  DO 500 K=1,25
!     ARLFXB=ARLFXB+ARLF(K,NB,NZ,NY,NX)
    IF(ARLF(K,NB,NZ,NY,NX).GT.0.0)THEN
      DO 700 L=JC,1,-1
!       ARLFXL=ARLFXL+ARLFL(L,K,NB,NZ,NY,NX)
        DO 800 N=1,4
          SURF(N,L,K,NB,NZ,NY,NX)=AMAX1(0.0_r8,CLASS(N,NZ,NY,NX) &
            *0.25_r8*ARLFL(L,K,NB,NZ,NY,NX))
  !       SURFXX=SURFXX+SURF(N,L,K,NB,NZ,NY,NX)
  !       IF(NZ.EQ.2)THEN
  !       WRITE(*,6363)'SURF',I,J,NX,NY,NZ,NB,K,L,N
  !      2,ARLFL(L,K,NB,NZ,NY,NX)
  !      2,SURF(N,L,K,NB,NZ,NY,NX),CLASS(N,NZ,NY,NX),ARLF(K,NB,NZ,NY,NX)
  !      3,ARLFXB,ARLFXL,SURFXX,ARLF(0,NB,NZ,NY,NX)
  !      4,ARLFB(NB,NZ,NY,NX),ZC(NZ,NY,NX)
  !6363    FORMAT(A8,9I4,20E12.4)
  !       ENDIF
800     CONTINUE
700   CONTINUE
    ENDIF
500 CONTINUE
!
! ALLOCATE STALK AREA TO INCLINATION CLASSES ACCORDING TO
! BRANCH ANGLE ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
!
! SURFB=stalk surface area in canopy layer
! ANGBR=stem angle from horizontal
! ARSTK=total branch stalk surface area in each layer
!
  DO 910 L=1,JC
    DO  N=1,4
      SURFB(N,L,NB,NZ,NY,NX)=0._r8
    enddo
910   CONTINUE
  IF(NB.EQ.NB1(NZ,NY,NX))THEN
    N=4
  ELSE
    N=MIN(4,INT(ASIN(ANGBR(NZ,NY,NX))/0.3927)+1)
  ENDIF
  DO 710 L=JC,1,-1
    SURFB(N,L,NB,NZ,NY,NX)=0.25*ARSTK(L,NB,NZ,NY,NX)
710   CONTINUE
  end subroutine LeafClassAllocation
!------------------------------------------------------------------------------------------
  subroutine CarbNutInBranchTransfer(I,J,NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NB,NZ,NY,NX

  ! begin_execution
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!
!   IF(NZ.EQ.2)THEN
!   WRITE(*,2322)'VRNS',I,J,NX,NY,NZ,NB,NB1(NZ,NY,NX),IFLGZ
!    2,ISTYP(NZ,NY,NX),IFLGI(NZ,NY,NX),IFLGE(NB,NZ,NY,NX)
!    3,IFLGF(NB,NZ,NY,NX),IDAY0(NZ,NY,NX),IYR0(NZ,NY,NX)
!    3,VRNS(NB1(NZ,NY,NX),NZ,NY,NX),VRNL(NB,NZ,NY,NX)
!    3,VRNF(NB,NZ,NY,NX),VRNX(NB,NZ,NY,NX)
!2322  FORMAT(A8,14I4,20E12.4)
!   ENDIF
  IF((ISTYP(NZ,NY,NX).EQ.0.AND.IFLGI(NZ,NY,NX).EQ.0) &
    .OR.(I.GE.IDAY0(NZ,NY,NX).AND.IYRC.EQ.IYR0(NZ,NY,NX) &
    .AND.VRNF(NB,NZ,NY,NX) &
    .LT.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX)) &
    .OR.(VRNS(NB1(NZ,NY,NX),NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX) &
    .AND.VRNF(NB,NZ,NY,NX) &
    .LT.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX)))THEN
    WTRTM=0._r8
    CPOOLM=0._r8
    DO 4 L=NU(NY,NX),NI(NZ,NY,NX)
      WTRTM=WTRTM+AMAX1(0.0,WTRTD(1,L,NZ,NY,NX))
      CPOOLM=CPOOLM+AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))
4   CONTINUE
!
  ! RESET TIME COUNTER
  !
  ! ATRP=hourly leafout counter
  ! IFLGA=flag for initializing leafout
  !
    IF(IFLGA(NB,NZ,NY,NX).EQ.0)THEN
      ATRP(NB,NZ,NY,NX)=0._r8
      IFLGA(NB,NZ,NY,NX)=1
    ENDIF
  !
  ! INCREMENT TIME COUNTER
  !
  ! IPTYP=photoperiod type:0=day neutral,1=short day,2=long day
  ! IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  ! XDL=critical photoperiod (h):<0=maximum daylength from site file
  ! XPPD=photoperiod sensitivity (node h-1)
  ! DYLN=daylength
  ! WFNSG=expansion,extension function of canopy water potential
  ! TFN3=temperature function for canopy growth
  ! ATRPX=number of hours required to initiate remobilization of storage C for leafout
  !
    IF(NB.EQ.NB1(NZ,NY,NX))THEN
      IF(IPTYP(NZ,NY,NX).EQ.2.AND.(IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3))THEN
        PPDX=AMAX1(0.0,XDL(NZ,NY,NX)-XPPD(NZ,NY,NX)-DYLN(NY,NX))
        ATRPPD=EXP(-0.0*PPDX)
      ELSE
        ATRPPD=1.0_r8
      ENDIF
      IF(IGTYP(NZ,NY,NX).NE.0)THEN
        WFNSP=WFNSG
      ELSE
        WFNSP=1.0_r8
      ENDIF
      DATRP=ATRPPD*TFN3(NZ,NY,NX)*WFNSP
      ATRP(NB,NZ,NY,NX)=ATRP(NB,NZ,NY,NX)+DATRP
!     IF(NZ.EQ.2)THEN
!     WRITE(*,2323)'ATRP',I,J,NX,NY,NZ,NB,ATRP(NB,NZ,NY,NX),DATRP
!    2,ATRPPD,TFN3(NZ,NY,NX),WFNSG,PPDX,XDL(NZ,NY,NX),XPPD(NZ,NY,NX)
!    3,DYLN(NY,NX),WTLFB(NB,NZ,NY,NX),ARLFB(NB,NZ,NY,NX)
!    4,HTCTL(NZ,NY,NX)
!2323  FORMAT(A8,6I4,20E12.4)
!     ENDIF
      IF(ATRP(NB,NZ,NY,NX).LE.ATRPX(ISTYP(NZ,NY,NX)) &
        .OR.(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).EQ.0))THEN
        IF(WTRVC(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
          CPOOLT=CPOOLM+CPOOL(NB,NZ,NY,NX)
  !
  !       REMOBILIZE C FROM SEASONAL STORAGE AT FIRST-ORDER RATE
  !       MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
  !
    !     GVMX=specific oxidation rate of storage C during leafout at 25 C
    !     WTRVC=storage C
    !     CH2OH=storage C oxidation rate during leafout
    !     CPOOL,CPOOLR=non-structural C mass in branch,root
    !     FXSH,FXRT=shoot-root partitioning of storage C during leafout
    !     WTRTD=root C mass
!
          GFNX=GVMX(ISTYP(NZ,NY,NX))*DATRP
          CH2OH=AMAX1(0.0,GFNX*WTRVC(NZ,NY,NX))
      !   IF(NZ.EQ.2)THEN
      !   WRITE(*,2123)'GERM0',I,J,NX,NY,NZ,NB
      !    2,GFNX,CH2OH,WTRVC(NZ,NY,NX)
      !    2,CPOOL(NB,NZ,NY,NX),CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
      !    3,FXSH(ISTYP(NZ,NY,NX)),FXRT(ISTYP(NZ,NY,NX))
      !2123  FORMAT(A8,6I4,20E12.4)
      !   ENDIF
          WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-CH2OH
          CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+CH2OH*FXSH(ISTYP(NZ,NY,NX))
          IF(WTRTM.GT.ZEROP(NZ,NY,NX).AND.CPOOLM.GT.ZEROP(NZ,NY,NX))THEN
            DO 50 L=NU(NY,NX),NI(NZ,NY,NX)
              FXFC=AMAX1(0.0,WTRTD(1,L,NZ,NY,NX))/WTRTM
              CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX) &
                +FXFC*CH2OH*FXRT(ISTYP(NZ,NY,NX))
50          CONTINUE
          ELSE
            CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)+CH2OH*FXRT(ISTYP(NZ,NY,NX))
          ENDIF
        ELSE
          CH2OH=0._r8
        ENDIF
      ELSE
        CH2OH=0._r8
      ENDIF
      !
      !     REMOBILIZE N,P FROM SEASONAL STORAGE AT FIRST-ORDER RATE
      !     MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
      !
      !     WTRVC,WTRVN,WTRVP=storage C,N,P
      !     ISTYP=growth habit:0=annual,1=perennial from PFT file
      !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
      !     UPNH4B,UPPO4B=N,P transfer from storage to shoot
      !     CH2OH=storage C oxidation rate during leafout
      !     FRSV=rate constant for remobiln of storage C,N,P during leafout C
      !     FXSH=shoot partitioning of storage C during leafout
      !
      IF(WTRVC(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        IF(ISTYP(NZ,NY,NX).NE.0)THEN
          CPOOLT=AMAX1(0.0,WTRVC(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX))
          ZPOOLD=(WTRVN(NZ,NY,NX)*CPOOL(NB,NZ,NY,NX)-ZPOOL(NB,NZ,NY,NX)*WTRVC(NZ,NY,NX))/CPOOLT
          PPOOLD=(WTRVP(NZ,NY,NX)*CPOOL(NB,NZ,NY,NX)-PPOOL(NB,NZ,NY,NX)*WTRVC(NZ,NY,NX))/CPOOLT
          UPNH4B=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*ZPOOLD)
          UPPO4B=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*PPOOLD)
        ELSE
          UPNH4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))*CH2OH*WTRVN(NZ,NY,NX)/WTRVC(NZ,NY,NX))
          UPPO4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))*CH2OH*WTRVP(NZ,NY,NX)/WTRVC(NZ,NY,NX))
        ENDIF
      ELSE
        UPNH4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))*WTRVN(NZ,NY,NX))
        UPPO4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))*WTRVP(NZ,NY,NX))
      ENDIF
    !
    ! ADD TO NON-STRUCTURAL POOLS IN ROOT
    !
    ! CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
    ! WTRVC,WTRVN,WTRVP=storage C,N,P
    ! ISTYP=growth habit:0=annual,1=perennial from PFT file
    ! UPNH4R,UPPO4R=N,P transfer from storage to root
    ! FRSV=rate constant for remobiln of storage C,N,P during leafout
    ! FXRT=root partitioning of storage C during leafout
    !
      CPOOLM=0._r8
      ZPOOLM=0._r8
      PPOOLM=0._r8
      DO 3 L=NU(NY,NX),NI(NZ,NY,NX)
        CPOOLM=CPOOLM+AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))
        ZPOOLM=ZPOOLM+AMAX1(0.0,ZPOOLR(1,L,NZ,NY,NX))
        PPOOLM=PPOOLM+AMAX1(0.0,PPOOLR(1,L,NZ,NY,NX))
3     CONTINUE
      IF(WTRVC(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        IF(ISTYP(NZ,NY,NX).NE.0)THEN
          CPOOLT=AMAX1(ZEROP(NZ,NY,NX),WTRVC(NZ,NY,NX)+CPOOLM)
          ZPOOLD=(WTRVN(NZ,NY,NX)*CPOOLM-ZPOOLM*WTRVC(NZ,NY,NX))/CPOOLT
          PPOOLD=(WTRVP(NZ,NY,NX)*CPOOLM-PPOOLM*WTRVC(NZ,NY,NX))/CPOOLT
          UPNH4R=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*ZPOOLD)
          UPPO4R=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*PPOOLD)
!         IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!         WRITE(*,9878)'GERM1',I,J,NZ,UPNH4R,FRSV(IBTYP(NZ,NY,NX))
!        2,ZPOOLD,WTRVN(NZ,NY,NX),CPOOLM,ZPOOLM,WTRVC(NZ,NY,NX)
!        3,CPOOLT
!9878     FORMAT(A8,3I4,12E24.16)
!         ENDIF
        ELSE
          UPNH4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))*CH2OH*WTRVN(NZ,NY,NX)/WTRVC(NZ,NY,NX))
          UPPO4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))*CH2OH*WTRVP(NZ,NY,NX)/WTRVC(NZ,NY,NX))
        ENDIF
      ELSE
        UPNH4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))*WTRVN(NZ,NY,NX))
        UPPO4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))*WTRVP(NZ,NY,NX))
      ENDIF
!
!     TRANSFER STORAGE FLUXES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     UPNH4B,UPPO4B=N,P transfer from storage to shoot
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     UPNH4R,UPPO4R=N,P transfer from storage to root
!     FXFN=root layer allocation
!
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)-UPNH4B-UPNH4R
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)-UPPO4B-UPPO4R
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+UPNH4B
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+UPPO4B
      IF(WTRTM.GT.ZEROP(NZ,NY,NX).AND.CPOOLM.GT.ZEROP(NZ,NY,NX))THEN
        DO 51 L=NU(NY,NX),NI(NZ,NY,NX)
          FXFN=AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))/CPOOLM
!         IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!         WRITE(*,9879)'GERM2',I,J,NZ,L,UPNH4R,FXFN
!        2,ZPOOLR(1,L,NZ,NY,NX),CPOOLR(1,L,NZ,NY,NX),CPOOLM
!9879      FORMAT(A8,4I4,12E24.16)
!         ENDIF
          ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)+FXFN*UPNH4R
          PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)+FXFN*UPPO4R
51      CONTINUE
      ELSE
  !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
  !     WRITE(*,9879)'GERM3',I,J,NZ,L,UPNH4R,FXFN
  !    2,ZPOOLR(1,L,NZ,NY,NX),CPOOLR(1,L,NZ,NY,NX),CPOOLM
  !     ENDIF
        ZPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=ZPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)+UPNH4R
        PPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=PPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)+UPPO4R
      ENDIF
    ENDIF
  !
  ! REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
  !
  ! ATRP=hourly leafout counter
  ! TFN3=temperature function for canopy growth
  ! ATRPX=number of hours required for remobilization of storage C during leafout
  ! WFNG=growth function of canopy water potential
  ! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  ! XFRC,XFRN,XFRC=nonstructural C,N,P transfer
  !
    IF(NB.NE.NB1(NZ,NY,NX).AND.ATRP(NB,NZ,NY,NX) &
      .LE.ATRPX(ISTYP(NZ,NY,NX)))THEN
      ATRP(NB,NZ,NY,NX)=ATRP(NB,NZ,NY,NX)+TFN3(NZ,NY,NX)*WFNG
      XFRC=AMAX1(0.0,0.05*TFN3(NZ,NY,NX) &
        *(0.5*(CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)+CPOOL(NB,NZ,NY,NX))-CPOOL(NB,NZ,NY,NX)))
      XFRN=AMAX1(0.0,0.05*TFN3(NZ,NY,NX) &
        *(0.5*(ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX))-ZPOOL(NB,NZ,NY,NX)))
      XFRP=AMAX1(0.0,0.05*TFN3(NZ,NY,NX) &
      *(0.5*(PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)+PPOOL(NB,NZ,NY,NX))-PPOOL(NB,NZ,NY,NX)))
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+XFRC
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+XFRP
      CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)-XFRC
      ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)-XFRN
      PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)-XFRP
    ENDIF
  ENDIF

!
! TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
! IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
! IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
!
! ISTYP=growth habit:0=annual,1=perennial from PFT file
! IFLGZ=remobilization flag
! WVSTKB=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! FXFB=rate constant for plant-storage nonstructural C,N,P exchange
! XFRC,XFRN,XFRC=nonstructural C,N,P transfer
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! ISTYP=growth habit:0=annual,1=perennial from PFT file
! IFLGZ=remobilization flag
! WVSTKB=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! XFRC,XFRN,XFRC=nonstructural C,N,P transfer
! WTRVC,WTRVN,WTRVP=storage C,N,P
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!
  IF(IFLGZ.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
    IF(WVSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WTRSVB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CWTRSV=AMAX1(0.0,WTRSVB(NB,NZ,NY,NX)/WVSTKB(NB,NZ,NY,NX))
      CWTRSN=AMAX1(0.0,WTRSBN(NB,NZ,NY,NX)/WVSTKB(NB,NZ,NY,NX))
      CWTRSP=AMAX1(0.0,WTRSBP(NB,NZ,NY,NX)/WVSTKB(NB,NZ,NY,NX))
      CNR=CWTRSV/(CWTRSV+CWTRSN/CNKI)
      CPR=CWTRSV/(CWTRSV+CWTRSP/CPKI)
    ELSE
      CNR=0._r8
      CPR=0._r8
    ENDIF
    XFRCX=FXFB(IBTYP(NZ,NY,NX))*AMAX1(0.0,WTRSVB(NB,NZ,NY,NX))
    XFRNX=FXFB(IBTYP(NZ,NY,NX))*AMAX1(0.0,WTRSBN(NB,NZ,NY,NX))*(1.0+CNR)
    XFRPX=FXFB(IBTYP(NZ,NY,NX))*AMAX1(0.0,WTRSBP(NB,NZ,NY,NX))*(1.0+CPR)
    XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
    WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-XFRC
    WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+XFRC
    WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)-XFRN
    WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+XFRN
    WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)-XFRP
    WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+XFRP
    IF(CCPOLB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CNL=CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX)+CZPOLB(NB,NZ,NY,NX)/CNKI)
      CPL=CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX)+CPPOLB(NB,NZ,NY,NX)/CPKI)
    ELSE
      CNL=0._r8
      CPL=0._r8
    ENDIF
    XFRCX=FXFB(IBTYP(NZ,NY,NX))*AMAX1(0.0,CPOOL(NB,NZ,NY,NX))
    XFRNX=FXFB(IBTYP(NZ,NY,NX))*AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))*(1.0+CNL)
    XFRPX=FXFB(IBTYP(NZ,NY,NX))*AMAX1(0.0,PPOOL(NB,NZ,NY,NX))*(1.0+CPL)
    XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
    XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
    CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
    WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+XFRC
    ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
    WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+XFRN
    PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
    WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+XFRP
!   IF(NZ.EQ.1)THEN
!     WRITE(*,4490)'RSV',I,J,NZ,NB,XFRC,XFRN,WTRSVB(NB,NZ,NY,NX)
!    2,WTRSBN(NB,NZ,NY,NX),WTRVC(NZ,NY,NX),WTRVN(NZ,NY,NX)
!    3,CNR,CNL,CPOOL(NB,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX)
!    4,FXFB(IBTYP(NZ,NY,NX))
!4490  FORMAT(A8,4I4,20E12.4)
!     ENDIF
  ENDIF
!
!   TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES AND ROOTS TO RESERVES
!   IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN STALK RESERVES
!   AND LEAVES IN PERENNIALS ACCORDING TO CONCENTRATION DIFFERENCES
!
!   ISTYP=growth habit:0=annual,1=perennial from PFT file
!   IDAY(3,=start of stem elongation and setting max seed number
!   IDAY(8,=end date setting for final seed number
!   WTLSB=leaf+petiole mass
!   WVSTKB=stalk sapwood mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   FXFY,FXFZ=rate constant for plant-reserve nonstructural C,N,P exchange
!   XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!   CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  IF((ISTYP(NZ,NY,NX).EQ.0.AND.IDAY(8,NB,NZ,NY,NX).NE.0) &
    .OR.(ISTYP(NZ,NY,NX).EQ.1.AND.IDAY(3,NB,NZ,NY,NX).NE.0))THEN
    WTPLTT=WTLSB(NB,NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
    CPOOLT=CPOOL(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
    IF(WTPLTT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=(CPOOL(NB,NZ,NY,NX)*WVSTKB(NB,NZ,NY,NX) &
        -WTRSVB(NB,NZ,NY,NX)*WTLSB(NB,NZ,NY,NX))/WTPLTT
      XFRC=FXFY(ISTYP(NZ,NY,NX))*CPOOLD
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
    ENDIF
    IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOL(NB,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX) &
        -WTRSBN(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      PPOOLD=(PPOOL(NB,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX) &
        -WTRSBP(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      XFRN=FXFZ(ISTYP(NZ,NY,NX))*ZPOOLD
      XFRP=FXFZ(ISTYP(NZ,NY,NX))*PPOOLD
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+XFRP
    ENDIF
!     IF(NZ.EQ.1)THEN
!     WRITE(*,4488)'EXCHC',I,J,NX,NY,NZ,NB,NS,XFRC,XFRN
!    2,FXFZ(ISTYP(NZ,NY,NX)),WTRSVB(NB,NZ,NY,NX),CPOOL(NB,NZ,NY,NX)
!    3,WVSTKB(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX)
!    4,CPOOLT,CPOOLD,ZPOOL(NB,NZ,NY,NX),WTRSBN(NB,NZ,NY,NX)
!4488  FORMAT(A8,7I4,12E12.4)
!     ENDIF
    IF(ISTYP(NZ,NY,NX).EQ.0.AND.IDAY(8,NB,NZ,NY,NX).NE.0)THEN
      DO 2050 L=NU(NY,NX),NI(NZ,NY,NX)
        IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          WTRTRX=AMAX1(ZEROP(NZ,NY,NX),WTRTL(1,L,NZ,NY,NX)*FWODR(1))
          WTPLTX=WTRTRX+WVSTKB(NB,NZ,NY,NX)
          IF(WTPLTX.GT.ZEROP(NZ,NY,NX))THEN
            CPOOLD=(CPOOLR(1,L,NZ,NY,NX)*WVSTKB(NB,NZ,NY,NX)-WTRSVB(NB,NZ,NY,NX)*WTRTRX)/WTPLTX
            XFRC=AMAX1(0.0,FXFY(ISTYP(NZ,NY,NX))*CPOOLD)
            CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
            WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
            CPOOLT=CPOOLR(1,L,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
            IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
              ZPOOLD=(ZPOOLR(1,L,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX) &
                -WTRSBN(NB,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
              PPOOLD=(PPOOLR(1,L,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX) &
                -WTRSBP(NB,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
              XFRN=AMAX1(0.0,FXFZ(ISTYP(NZ,NY,NX))*ZPOOLD)
              XFRP=AMAX1(0.0,FXFZ(ISTYP(NZ,NY,NX))*PPOOLD)
              ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
              WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+XFRN
              PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
              WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+XFRP
          !     IF(NZ.EQ.1)THEN
          !     WRITE(*,4489)'EXCHC',I,J,NZ,NB,L,WTRSVB(NB,NZ,NY,NX)
          !    2,WVSTKB(NB,NZ,NY,NX),CPOOLR(1,L,NZ,NY,NX)
          !    3,WTRTL(1,L,NZ,NY,NX),FWOOD(1),WTRTRX,WTPLTX
          !    4,CPOOLT,CPOOLD,XFRC,FXFZ(ISTYP(NZ,NY,NX))
          !4489  FORMAT(A8,5I4,12E16.8)
          !     ENDIF
          !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
          !     WRITE(*,4489)'EXCHN',I,J,NZ,NB,L,WTRSBN(NB,NZ,NY,NX)
          !    2,WTRSVB(NB,NZ,NY,NX),ZPOOLR(1,L,NZ,NY,NX)
          !    3,CPOOLR(1,L,NZ,NY,NX),FWOOD(1),ZPOOLD,XFRN
          !     ENDIF
            ENDIF
          ENDIF
        ENDIF
2050  CONTINUE
    ENDIF
  ENDIF
!
!   REPLENISH BRANCH NON-STRUCTURAL POOL FROM
!   SEASONAL STORAGE POOL
!
!   WVSTKB,WVSTK=stalk,total stalk sapwood mass
!   WTRT=total root mass
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   XFRX=maximum storage C content for remobiln from stalk,root reserves
!   XFRC=C transfer
!
  IF(WVSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
    .AND.WVSTK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
    .AND.WTRT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
    .AND.WTRSVB(NB,NZ,NY,NX).LE.XFRX*WVSTKB(NB,NZ,NY,NX))THEN
    FWTBR=WVSTKB(NB,NZ,NY,NX)/WVSTK(NZ,NY,NX)
    WVSTBX=WVSTKB(NB,NZ,NY,NX)
    WTRTTX=WTRT(NZ,NY,NX)*FWTBR
    WTPLTT=WVSTBX+WTRTTX
    WTRSBX=AMAX1(0.0,WTRSVB(NB,NZ,NY,NX))
    WTRVCX=AMAX1(0.0,WTRVC(NZ,NY,NX)*FWTBR)
    CPOOLD=(WTRVCX*WVSTBX-WTRSBX*WTRTTX)/WTPLTT
    XFRC=AMAX1(0.0,XFRY*CPOOLD)
    WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
    WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-XFRC
  ENDIF
end subroutine CarbNutInBranchTransfer
!------------------------------------------------------------------------------------------
      subroutine StagePlantForGrowth(I,J,NZ,NY,NX)
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     IF(I.EQ.1.AND.J.EQ.1)THEN
!     DO 87 II=1,366
!     DO 87 N=1,400
!     VCO2(N,II,NZ)=0._r8
!     enddo
!87    CONTINUE
!     ENDIF
!     IF(IYRC.GE.2099)THEN
!     IF(I.EQ.365.AND.J.EQ.24)THEN
!     DO 88 N=1,400
!     WRITE(19,12)IYRC,NZ,N,(VCO2(N,II,NZ),II=1,181)
!     WRITE(20,12)IYRC,NZ,N,(VCO2(N,II,NZ),II=182,365)
!12    FORMAT(3I8,365E12.4)
!88    CONTINUE
!     ENDIF
!     ENDIF
      IFLGZ=0
      IFLGY=0
      DO 2 L=1,JC
      ARLFV(L,NZ,NY,NX)=0._r8
      WGLFV(L,NZ,NY,NX)=0._r8
      ARSTV(L,NZ,NY,NX)=0._r8
2     CONTINUE
      DO 5 NR=1,NRT(NZ,NY,NX)
      DO  N=1,MY(NZ,NY,NX)
      NRX(N,NR)=0
      ICHK1(N,NR)=0
      enddo
5     CONTINUE
      DO 9 N=1,MY(NZ,NY,NX)
      RTNT(N)=0._r8
      DO 6 L=NU(NY,NX),NJ(NY,NX)
      WSRTL(N,L,NZ,NY,NX)=0._r8
      RTN1(N,L,NZ,NY,NX)=0._r8
      RTNL(N,L,NZ,NY,NX)=0._r8
      RCO2M(N,L,NZ,NY,NX)=0._r8
      RCO2N(N,L,NZ,NY,NX)=0._r8
      RCO2A(N,L,NZ,NY,NX)=0._r8
      RLNT(N,L)=0._r8
      DO  NR=1,NRT(NZ,NY,NX)
      RTSK1(N,L,NR)=0._r8
      RTSK2(N,L,NR)=0._r8
      enddo
6     CONTINUE
9     CONTINUE
!
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     WTSTK,WVSTK=stalk,sapwood mass
!     FWOOD,FWODB=C woody fraction in stalk,other organs:0=woody,1=non-woody
!     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
!     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!     *LF=leaf,*SHE=petiole,*STK=stalk,*RT=root
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
!
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1 &
      .OR.WTSTK(NZ,NY,NX).LE.ZEROP(NZ,NY,NX))THEN
      FWODB(1)=1.0_r8
      FWOOD(1)=1.0_r8
      FWODR(1)=1.0_r8
      ELSE
      FWODB(1)=1.0_r8
      FWOOD(1)=SQRT(WVSTK(NZ,NY,NX)/WTSTK(NZ,NY,NX))
      FWODR(1)=SQRT(FRTX*WVSTK(NZ,NY,NX)/WTSTK(NZ,NY,NX))
      ENDIF
      FWODB(0)=1.0_r8-FWODB(1)
      FWOOD(0)=1.0_r8-FWOOD(1)
      FWODR(0)=1.0_r8-FWODR(1)
!     WRITE(*,8822)'FWOOD',I,J,NX,NY,NZ,IBTYP(NZ,NY,NX),IGTYP(NZ,NY,NX)
!    2,FWOOD(0),FWOOD(1),FWODB(0),FWODB(1),FWODR(0),FWODR(1)
!    3,WVSTK(NZ,NY,NX),WTSTK(NZ,NY,NX)
!8822  FORMAT(A8,7I4,12E12.4)
      CNLFW=FWODB(0)*CNSTK(NZ,NY,NX)+FWODB(1)*CNLF(NZ,NY,NX)
      CPLFW=FWODB(0)*CPSTK(NZ,NY,NX)+FWODB(1)*CPLF(NZ,NY,NX)
      CNSHW=FWODB(0)*CNSTK(NZ,NY,NX)+FWODB(1)*CNSHE(NZ,NY,NX)
      CPSHW=FWODB(0)*CPSTK(NZ,NY,NX)+FWODB(1)*CPSHE(NZ,NY,NX)
      CNRTW=FWODR(0)*CNSTK(NZ,NY,NX)+FWODR(1)*CNRT(NZ,NY,NX)
      CPRTW=FWODR(0)*CPSTK(NZ,NY,NX)+FWODR(1)*CPRT(NZ,NY,NX)
      FWODLN(0)=FWODB(0)*CNSTK(NZ,NY,NX)/CNLFW
      FWODLP(0)=FWODB(0)*CPSTK(NZ,NY,NX)/CPLFW
      FWODSN(0)=FWODB(0)*CNSTK(NZ,NY,NX)/CNSHW
      FWODSP(0)=FWODB(0)*CPSTK(NZ,NY,NX)/CPSHW
      FWOODN(0)=FWOOD(0)*CNSTK(NZ,NY,NX)/CNRTW
      FWOODP(0)=FWOOD(0)*CPSTK(NZ,NY,NX)/CPRTW
      FWODRN(0)=FWODR(0)*CNRT(NZ,NY,NX)/CNRTW
      FWODRP(0)=FWODR(0)*CPRT(NZ,NY,NX)/CPRTW
      FWODLN(1)=1.0_r8-FWODLN(0)
      FWODLP(1)=1.0_r8-FWODLP(0)
      FWODSN(1)=1.0_r8-FWODSN(0)
      FWODSP(1)=1.0_r8-FWODSP(0)
      FWOODN(1)=1.0_r8-FWOODN(0)
      FWOODP(1)=1.0_r8-FWOODP(0)
      FWODRN(1)=1.0_r8-FWODRN(0)
      FWODRP(1)=1.0_r8-FWODRP(0)
!
!     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
!     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
!
!     TKC,TKCM=canopy temperature,canopy temp used in Arrhenius eqn
!     TKS,TKSM=soil temperature,soil temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN5,TFN6=temperature function for canopy,root mntc respn (25 oC =1)
!     8.3143,710.0=gas constant,enthalpy
!     62500,195000,232500=energy of activn,high,low temp inactivn(KJ mol-1)
!
      TKCM=TKC(NZ,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKCM
      STK=710.0*TKCM
      ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
      TFN5=EXP(25.214-62500/RTK)/ACTVM
      DO 7 L=NU(NY,NX),NJ(NY,NX)
      TKSM=TKS(L,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKSM
      STK=710.0*TKSM
      ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
      TFN6(L)=EXP(25.214-62500/RTK)/ACTVM
7     CONTINUE
      GROGR=0._r8
!
!     PRIMARY ROOT NUMBER
!
!     WTRTA=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     XRTN1=multiplier for number of primary root axes
!
      WTRTA(NZ,NY,NX)=AMAX1(0.999992087*WTRTA(NZ,NY,NX) &
      ,WTRT(NZ,NY,NX)/PP(NZ,NY,NX))
      XRTN1=AMAX1(1.0,WTRTA(NZ,NY,NX)**0.667)*PP(NZ,NY,NX)
!
!     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
!     FROM CANOPY TURGOR
!
!     WFNS=turgor expansion,extension function
!     PSILG,PSILM=current,minimum canopy turgor potl for expansion,extension
!     WFNC=stomatal resistance function of canopy turgor
!     PSILT=canopy water potential
!     WFNG=growth function of canopy water potential
!     WFNSG=expansion,extension function of canopy water potential
!
      WFNS=AMIN1(1.0,AMAX1(0.0,PSILG(NZ,NY,NX)-PSILM))
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      WFNC=1.0_r8
      WFNG=EXP(0.05*PSILT(NZ,NY,NX))
      WFNSG=WFNS**0.10
      ELSE
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      WFNG=EXP(0.10*PSILT(NZ,NY,NX))
      WFNSG=WFNS**0.25
      ENDIF
      end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

      subroutine NonstructlBiomTransfer(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ATRP=hourly leafout counter
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     WTLSB=leaf+petiole mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!
      IF(NBR(NZ,NY,NX).GT.1)THEN
      WTPLTT=0._r8
      CPOOLT=0._r8
      ZPOOLT=0._r8
      PPOOLT=0._r8
      DO 300 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(ATRP(NB,NZ,NY,NX).GT.ATRPX(ISTYP(NZ,NY,NX)))THEN
      WTLSBZ(NB)=AMAX1(0.0,WTLSB(NB,NZ,NY,NX))
      CPOOLZ(NB)=AMAX1(0.0,CPOOL(NB,NZ,NY,NX))
      ZPOOLZ(NB)=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      PPOOLZ(NB)=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
      WTPLTT=WTPLTT+WTLSBZ(NB)
      CPOOLT=CPOOLT+CPOOLZ(NB)
      ZPOOLT=ZPOOLT+ZPOOLZ(NB)
      PPOOLT=PPOOLT+PPOOLZ(NB)
      ENDIF
      ENDIF
300   CONTINUE
      DO 305 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(ATRP(NB,NZ,NY,NX).GT.ATRPX(ISTYP(NZ,NY,NX)))THEN
      IF(WTPLTT.GT.ZEROP(NZ,NY,NX) &
      .AND.CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=CPOOLT*WTLSBZ(NB)-CPOOLZ(NB)*WTPLTT
      ZPOOLD=ZPOOLT*CPOOLZ(NB)-ZPOOLZ(NB)*CPOOLT
      PPOOLD=PPOOLT*CPOOLZ(NB)-PPOOLZ(NB)*CPOOLT
      XFRC=0.01*CPOOLD/WTPLTT
      XFRN=0.01*ZPOOLD/CPOOLT
      XFRP=0.01*PPOOLD/CPOOLT
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+XFRC
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+XFRP
      ENDIF
      ENDIF
      ENDIF
305   CONTINUE
      ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     WVSTKB=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IDAY(7,=start of grain filling and setting max seed size
!
      IF(NBR(NZ,NY,NX).GT.1)THEN
      WTSTKT=0._r8
      WTRSVT=0._r8
      WTRSNT=0._r8
      WTRSPT=0._r8
      DO 330 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(IDAY(7,NB,NZ,NY,NX).NE.0)THEN
      WTSTKT=WTSTKT+WVSTKB(NB,NZ,NY,NX)
      WTRSVT=WTRSVT+WTRSVB(NB,NZ,NY,NX)
      WTRSNT=WTRSNT+WTRSBN(NB,NZ,NY,NX)
      WTRSPT=WTRSPT+WTRSBP(NB,NZ,NY,NX)
      ENDIF
      ENDIF
330   CONTINUE
      IF(WTSTKT.GT.ZEROP(NZ,NY,NX) &
      .AND.WTRSVT.GT.ZEROP(NZ,NY,NX))THEN
      DO 335 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(IDAY(7,NB,NZ,NY,NX).NE.0)THEN
      WTRSVD=WTRSVT*WVSTKB(NB,NZ,NY,NX) &
      -WTRSVB(NB,NZ,NY,NX)*WTSTKT
      XFRC=0.1*WTRSVD/WTSTKT
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
      WTRSND=WTRSNT*WTRSVB(NB,NZ,NY,NX) &
      -WTRSBN(NB,NZ,NY,NX)*WTRSVT
      XFRN=0.1*WTRSND/WTRSVT
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+XFRN
      WTRSPD=WTRSPT*WTRSVB(NB,NZ,NY,NX) &
      -WTRSBP(NB,NZ,NY,NX)*WTRSVT
      XFRP=0.1*WTRSPD/WTRSVT
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+XFRP
      ENDIF
      ENDIF
335   CONTINUE
      ENDIF
      ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BWTWEEN ROOT AND MYCORRHIZAE
!     IN EACH ROOTED SOIL LAYER FROM NON-STRUCTURAL C,N,P
!     CONCENTRATION DIFFERENCES
!
!     MY=mycorrhizal:1=no,2=yes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in 1:root,2:mycorrhizae
!     WTRTD=1:root,2:mycorrhizal C mass
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
!
      IF(MY(NZ,NY,NX).EQ.2)THEN
      DO 425 L=NU(NY,NX),NIX(NZ,NY,NX)
      IF(CPOOLR(1,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.WTRTD(1,L,NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
      WTRTD1=WTRTD(1,L,NZ,NY,NX)
      WTRTD2=AMIN1(WTRTD(1,L,NZ,NY,NX),AMAX1(FSNK &
      *WTRTD(1,L,NZ,NY,NX),WTRTD(2,L,NZ,NY,NX)))
      WTPLTT=WTRTD1+WTRTD2
      IF(WTPLTT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=(CPOOLR(1,L,NZ,NY,NX)*WTRTD2 &
      -CPOOLR(2,L,NZ,NY,NX)*WTRTD1)/WTPLTT
      XFRC=FMYC*CPOOLD
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
      CPOOLR(2,L,NZ,NY,NX)=CPOOLR(2,L,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLR(1,L,NZ,NY,NX)+CPOOLR(2,L,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOLR(1,L,NZ,NY,NX)*CPOOLR(2,L,NZ,NY,NX) &
      -ZPOOLR(2,L,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRN=FMYC*ZPOOLD
      PPOOLD=(PPOOLR(1,L,NZ,NY,NX)*CPOOLR(2,L,NZ,NY,NX) &
      -PPOOLR(2,L,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRP=FMYC*PPOOLD
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
      ZPOOLR(2,L,NZ,NY,NX)=ZPOOLR(2,L,NZ,NY,NX)+XFRN
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
      PPOOLR(2,L,NZ,NY,NX)=PPOOLR(2,L,NZ,NY,NX)+XFRP
!     IF(L.EQ.NIX(NZ,NY,NX))THEN
!     WRITE(*,9873)'MYCO',I,J,NZ,L,XFRC,XFRN,XFRP
!    2,CPOOLR(1,L,NZ,NY,NX),WTRTD(1,L,NZ,NY,NX)
!    3,CPOOLR(2,L,NZ,NY,NX),WTRTD2
!    3,WTPLTT,ZPOOLR(1,L,NZ,NY,NX),ZPOOLR(2,L,NZ,NY,NX)
!    4,PPOOLR(1,L,NZ,NY,NX),PPOOLR(2,L,NZ,NY,NX),CPOOLT
!9873  FORMAT(A8,4I4,20E24.16)
!     ENDIF
      ENDIF
      ENDIF
      ENDIF
425   CONTINUE
      ENDIF
!
!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
      IF(IFLGZ.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      DO 5545 N=1,MY(NZ,NY,NX)
      DO 5550 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CNL=CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX) &
      +CZPOLR(N,L,NZ,NY,NX)/CNKI)
      CPL=CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX) &
      +CPPOLR(N,L,NZ,NY,NX)/CPKI)
      ELSE
      CNL=0._r8
      CPL=0._r8
      ENDIF
      XFRCX=FXFR(IGTYP(NZ,NY,NX)) &
      *AMAX1(0.0,CPOOLR(N,L,NZ,NY,NX))
      XFRNX=FXFR(IGTYP(NZ,NY,NX)) &
      *AMAX1(0.0,ZPOOLR(N,L,NZ,NY,NX))*(1.0+CNL)
      XFRPX=FXFR(IGTYP(NZ,NY,NX)) &
      *AMAX1(0.0,PPOOLR(N,L,NZ,NY,NX))*(1.0+CPL)
      XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
      XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
      XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+XFRC
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)-XFRN
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+XFRN
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)-XFRP
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+XFRP
5550  CONTINUE
5545  CONTINUE
      ENDIF
!
!     ROOT AND NODULE TOTALS
!
!     WTRTL,WTRTD=active,actual root C mass
!     WTRT1,WTRT2=primary,secondary root C mass in soil layer
!     TCO2T=total PFT respiration
!     RCO2A=total root respiration
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
      DO 5445 N=1,MY(NZ,NY,NX)
      DO 5450 L=NU(NY,NX),NI(NZ,NY,NX)
      WTRTL(N,L,NZ,NY,NX)=0._r8
      WTRTD(N,L,NZ,NY,NX)=0._r8
      DO 5460 NR=1,NRT(NZ,NY,NX)
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX)
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX) &
      +WTRT1(N,L,NR,NZ,NY,NX)
5460  CONTINUE
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)+RCO2A(N,L,NZ,NY,NX)
      RECO(NY,NX)=RECO(NY,NX)+RCO2A(N,L,NZ,NY,NX)
      TRAU(NY,NX)=TRAU(NY,NX)+RCO2A(N,L,NZ,NY,NX)
5450  CONTINUE
      DO 5470 NR=1,NRT(NZ,NY,NX)
      WTRTL(N,NINR(NR,NZ,NY,NX),NZ,NY,NX) &
      =WTRTL(N,NINR(NR,NZ,NY,NX),NZ,NY,NX) &
      +RTWT1(N,NR,NZ,NY,NX)
5470  CONTINUE
5445  CONTINUE
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
!
!     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
!     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     WTLS,WTRT=total PFT leaf+petiole,root C mass
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     RLNT,RTNT=root layer,root system sink strength
!
!     IF(ISTYP(NZ,NY,NX).EQ.1)THEN
      IF(WTLS(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FWTC=AMIN1(1.0,0.667*WTRT(NZ,NY,NX)/WTLS(NZ,NY,NX))
      ELSE
      FWTC=1.0_r8
      ENDIF
      IF(WTRT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FWTS=AMIN1(1.0,WTLS(NZ,NY,NX)/(0.667*WTRT(NZ,NY,NX)))
      ELSE
      FWTS=1.0_r8
      ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
      DO 290 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(RTNT(1).GT.ZEROP(NZ,NY,NX))THEN
      FWTR(L)=AMAX1(0.0,RLNT(1,L)/RTNT(1))
      ELSE
      FWTR(L)=1.0_r8
      ENDIF
290   CONTINUE
!     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
!     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
!
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!
      WTLS(NZ,NY,NX)=0._r8
      DO 309 NB=1,NBR(NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
309   CONTINUE
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAY(8,=end date for setting final seed number
!     FWTB=branch sink weighting factor
!     PTSHT=rate constant for equilibrating shoot-root nonstructural C concn from PFT file
!     PTRT=allocation to leaf+petiole used to modify PTSHT in annuals
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     FWODB=C woody fraction in branch:0=woody,1=non-woody
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
      DO 310 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(WTLS(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FWTB(NB)=AMAX1(0.0,WTLSB(NB,NZ,NY,NX)/WTLS(NZ,NY,NX))
      ELSE
      FWTB(NB)=1.0_r8
      ENDIF
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      PTSHTR=PTSHT(NZ,NY,NX)*PTRT**0.167
      ELSE
      PTSHTR=PTSHT(NZ,NY,NX)
      ENDIF
      DO 415 L=NU(NY,NX),NI(NZ,NY,NX)
      WTLSBX=WTLSB(NB,NZ,NY,NX)*FWODB(1)*FWTR(L)*FWTC
      WTRTLX=WTRTL(1,L,NZ,NY,NX)*FWODR(1)*FWTB(NB)*FWTS
      WTLSBB=AMAX1(0.0,WTLSBX,FSNK*WTRTLX)
      WTRTLR=AMAX1(0.0,WTRTLX,FSNK*WTLSBX)
      WTPLTT=WTLSBB+WTRTLR
      IF(WTPLTT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLB=AMAX1(0.0,CPOOL(NB,NZ,NY,NX)*FWTR(L))
      CPOOLS=AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX)*FWTB(NB))
      CPOOLD=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
      XFRC=PTSHTR*CPOOLD
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLS+CPOOLB
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX)*FWTR(L))
      ZPOOLS=AMAX1(0.0,ZPOOLR(1,L,NZ,NY,NX)*FWTB(NB))
      ZPOOLD=(ZPOOLB*CPOOLS-ZPOOLS*CPOOLB)/CPOOLT
      XFRN=PTSHTR*ZPOOLD
      PPOOLB=AMAX1(0.0,PPOOL(NB,NZ,NY,NX)*FWTR(L))
      PPOOLS=AMAX1(0.0,PPOOLR(1,L,NZ,NY,NX)*FWTB(NB))
      PPOOLD=(PPOOLB*CPOOLS-PPOOLS*CPOOLB)/CPOOLT
      XFRP=PTSHTR*PPOOLD
      ELSE
      XFRN=0._r8
      XFRP=0._r8
      ENDIF
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)+XFRP
!     IF(NZ.EQ.2.AND.NB.EQ.1)THEN
!     WRITE(*,3344)'ROOT',I,J,NX,NY,NZ,NB,L
!    2,XFRC,XFRN,XFRP,CPOOL(NB,NZ,NY,NX)
!    3,CPOOLR(1,L,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX)
!    3,ZPOOLR(1,L,NZ,NY,NX),FWTB(NB),FWTR(L)
!    3,FWTC,FWTS,WTLSBX,WTRTLX,FSNK,FDBK(NB,NZ,NY,NX)
!    4,CPOOLD,CPOOLB,WTLSBB,CPOOLS,WTRTLR
!    5,FWOOD(1),FWODB(1),WTRTL(1,L,NZ,NY,NX)
!    6,WTLSB(NB,NZ,NY,NX),RLNT(1,L),RTNT(1)
!3344  FORMAT(A8,7I4,30E12.4)
!     ENDIF
      ENDIF
415   CONTINUE
      ENDIF
310   CONTINUE
      end subroutine NonstructlBiomTransfer
!------------------------------------------------------------------------------------------

      subroutine ComputeTotalBiom(NZ,NY,NX)

      integer, intent(in) :: NZ,NY,NX
!     begin_execution
!     TOTAL C,N,P IN EACH BRANCH
!
!     CPOOLK=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     CPOOL,ZPOOL,PPOOL=C3 non-structural C,N,P mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
!
      DO 320 NB=1,NBR(NZ,NY,NX)
      CPOOLK(NB,NZ,NY,NX)=0._r8
      DO 325 K=1,25
      CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX) &
      +CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX) &
      +CO2B(K,NB,NZ,NY,NX)+HCOB(K,NB,NZ,NY,NX)
325   CONTINUE
      WTSHTB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX) &
      +WTSHEB(NB,NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX) &
      +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX) &
      +CPOOL(NB,NZ,NY,NX)+CPOOLK(NB,NZ,NY,NX)
      WTSHTN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX) &
      +WTSHBN(NB,NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX) &
      +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX) &
      +ZPOOL(NB,NZ,NY,NX)
      WTSHTP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX) &
      +WTSHBP(NB,NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX) &
      +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX) &
      +PPOOL(NB,NZ,NY,NX)
320   CONTINUE
!
!     TOTAL C,N,P IN ROOTS AND MYCORRHIZAE IN EACH SOIL LAYER
!
!     WTRTD=root C mass
!     CPOOLR=non-structural C mass in root
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
!     UPNF=PFT N2 fixation
!
      DO 345 N=1,MY(NZ,NY,NX)
      DO  L=NU(NY,NX),NI(NZ,NY,NX)
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
      enddo
345   CONTINUE
      end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

      subroutine AccumulateStates(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     INITIALIZE PFT STATE VARIABLES
!
      CPOOLP(NZ,NY,NX)=0._r8
      ZPOOLP(NZ,NY,NX)=0._r8
      PPOOLP(NZ,NY,NX)=0._r8
      WTSHT(NZ,NY,NX)=0._r8
      WTSHN(NZ,NY,NX)=0._r8
      WTSHP(NZ,NY,NX)=0._r8
      WTLF(NZ,NY,NX)=0._r8
      WTSHE(NZ,NY,NX)=0._r8
      WTSTK(NZ,NY,NX)=0._r8
      WVSTK(NZ,NY,NX)=0._r8
      WTRSV(NZ,NY,NX)=0._r8
      WTHSK(NZ,NY,NX)=0._r8
      WTEAR(NZ,NY,NX)=0._r8
      WTGR(NZ,NY,NX)=0._r8
      WTLS(NZ,NY,NX)=0._r8
      WTRT(NZ,NY,NX)=0._r8
      WTRTS(NZ,NY,NX)=0._r8
      WTRTN(NZ,NY,NX)=0._r8
      WTRTP(NZ,NY,NX)=0._r8
      WTLFN(NZ,NY,NX)=0._r8
      WTSHEN(NZ,NY,NX)=0._r8
      WTSTKN(NZ,NY,NX)=0._r8
      WTRSVN(NZ,NY,NX)=0._r8
      WTHSKN(NZ,NY,NX)=0._r8
      WTEARN(NZ,NY,NX)=0._r8
      WTGRNN(NZ,NY,NX)=0._r8
      WTLFP(NZ,NY,NX)=0._r8
      WTSHEP(NZ,NY,NX)=0._r8
      WTSTKP(NZ,NY,NX)=0._r8
      WTRSVP(NZ,NY,NX)=0._r8
      WTHSKP(NZ,NY,NX)=0._r8
      WTEARP(NZ,NY,NX)=0._r8
      WTGRNP(NZ,NY,NX)=0._r8
      GRNO(NZ,NY,NX)=0._r8
      ARLFP(NZ,NY,NX)=0._r8
      ARSTP(NZ,NY,NX)=0._r8
      DO 8940 L=1,JC
      ARSTV(L,NZ,NY,NX)=0._r8
8940  CONTINUE
!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     ARLFB=branch leaf area
!     ARSTK=total branch stalk surface area in each layer
!     GRNOB=seed set number
!
      DO 8950 NB=1,NBR(NZ,NY,NX)
      CPOOLP(NZ,NY,NX)=CPOOLP(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX)
      ZPOOLP(NZ,NY,NX)=ZPOOLP(NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX)
      PPOOLP(NZ,NY,NX)=PPOOLP(NZ,NY,NX)+PPOOL(NB,NZ,NY,NX)
      WTSHT(NZ,NY,NX)=WTSHT(NZ,NY,NX)+WTSHTB(NB,NZ,NY,NX)
      WTLF(NZ,NY,NX)=WTLF(NZ,NY,NX)+WTLFB(NB,NZ,NY,NX)
      WTSHE(NZ,NY,NX)=WTSHE(NZ,NY,NX)+WTSHEB(NB,NZ,NY,NX)
      WTSTK(NZ,NY,NX)=WTSTK(NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)
      WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      WTRSV(NZ,NY,NX)=WTRSV(NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
      WTHSK(NZ,NY,NX)=WTHSK(NZ,NY,NX)+WTHSKB(NB,NZ,NY,NX)
      WTEAR(NZ,NY,NX)=WTEAR(NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)
      WTGR(NZ,NY,NX)=WTGR(NZ,NY,NX)+WTGRB(NB,NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
      WTSHN(NZ,NY,NX)=WTSHN(NZ,NY,NX)+WTSHTN(NB,NZ,NY,NX)
      WTLFN(NZ,NY,NX)=WTLFN(NZ,NY,NX)+WTLFBN(NB,NZ,NY,NX)
      WTSHEN(NZ,NY,NX)=WTSHEN(NZ,NY,NX)+WTSHBN(NB,NZ,NY,NX)
      WTSTKN(NZ,NY,NX)=WTSTKN(NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)
      WTRSVN(NZ,NY,NX)=WTRSVN(NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
      WTHSKN(NZ,NY,NX)=WTHSKN(NZ,NY,NX)+WTHSBN(NB,NZ,NY,NX)
      WTEARN(NZ,NY,NX)=WTEARN(NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)
      WTGRNN(NZ,NY,NX)=WTGRNN(NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX)
      WTSHP(NZ,NY,NX)=WTSHP(NZ,NY,NX)+WTSHTP(NB,NZ,NY,NX)
      WTLFP(NZ,NY,NX)=WTLFP(NZ,NY,NX)+WTLFBP(NB,NZ,NY,NX)
      WTSHEP(NZ,NY,NX)=WTSHEP(NZ,NY,NX)+WTSHBP(NB,NZ,NY,NX)
      WTSTKP(NZ,NY,NX)=WTSTKP(NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)
      WTRSVP(NZ,NY,NX)=WTRSVP(NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
      WTHSKP(NZ,NY,NX)=WTHSKP(NZ,NY,NX)+WTHSBP(NB,NZ,NY,NX)
      WTEARP(NZ,NY,NX)=WTEARP(NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)
      WTGRNP(NZ,NY,NX)=WTGRNP(NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX)
      ARLFP(NZ,NY,NX)=ARLFP(NZ,NY,NX)+ARLFB(NB,NZ,NY,NX)
      GRNO(NZ,NY,NX)=GRNO(NZ,NY,NX)+GRNOB(NB,NZ,NY,NX)
      DO 8945 L=1,JC
      ARSTP(NZ,NY,NX)=ARSTP(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
      ARSTV(L,NZ,NY,NX)=ARSTV(L,NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
8945  CONTINUE
8950  CONTINUE
!
!     ACCUMULATE ROOT STATE VARIABLES FROM ROOT LAYER STATE VARIABLES
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!
      DO 8925 N=1,MY(NZ,NY,NX)
      DO 8930 L=NU(NY,NX),NJ(NY,NX)
      WTRT(NZ,NY,NX)=WTRT(NZ,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
      WTRTN(NZ,NY,NX)=WTRTN(NZ,NY,NX)+ZPOOLR(N,L,NZ,NY,NX)
      WTRTP(NZ,NY,NX)=WTRTP(NZ,NY,NX)+PPOOLR(N,L,NZ,NY,NX)
      DO 8935 NR=1,NRT(NZ,NY,NX)
      WTRT(NZ,NY,NX)=WTRT(NZ,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX) &
      +WTRT2(N,L,NR,NZ,NY,NX)
      WTRTS(NZ,NY,NX)=WTRTS(NZ,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX) &
      +WTRT2(N,L,NR,NZ,NY,NX)
      WTRTN(NZ,NY,NX)=WTRTN(NZ,NY,NX)+WTRT1N(N,L,NR,NZ,NY,NX) &
      +WTRT2N(N,L,NR,NZ,NY,NX)
      WTRTP(NZ,NY,NX)=WTRTP(NZ,NY,NX)+WTRT1P(N,L,NR,NZ,NY,NX) &
      +WTRT2P(N,L,NR,NZ,NY,NX)
8935  CONTINUE
8930  CONTINUE
8925  CONTINUE
!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
      IF(INTYP(NZ,NY,NX).NE.0)THEN
      WTND(NZ,NY,NX)=0._r8
      WTNDN(NZ,NY,NX)=0._r8
      WTNDP(NZ,NY,NX)=0._r8
      IF(INTYP(NZ,NY,NX).GE.4)THEN
      DO 7950 NB=1,NBR(NZ,NY,NX)
      CPOLNP(NZ,NY,NX)=CPOLNP(NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
      ZPOLNP(NZ,NY,NX)=ZPOLNP(NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX)
      PPOLNP(NZ,NY,NX)=PPOLNP(NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX)
      WTND(NZ,NY,NX)=WTND(NZ,NY,NX)+WTNDB(NB,NZ,NY,NX) &
      +CPOLNB(NB,NZ,NY,NX)
      WTNDN(NZ,NY,NX)=WTNDN(NZ,NY,NX)+WTNDBN(NB,NZ,NY,NX) &
      +ZPOLNB(NB,NZ,NY,NX)
      WTNDP(NZ,NY,NX)=WTNDP(NZ,NY,NX)+WTNDBP(NB,NZ,NY,NX) &
      +PPOLNB(NB,NZ,NY,NX)
7950  CONTINUE
      ELSEIF(INTYP(NZ,NY,NX).GE.1.AND.INTYP(NZ,NY,NX).LE.3)THEN
      DO 8920 L=NU(NY,NX),NI(NZ,NY,NX)
      WTND(NZ,NY,NX)=WTND(NZ,NY,NX)+WTNDL(L,NZ,NY,NX) &
      +CPOOLN(L,NZ,NY,NX)
      WTNDN(NZ,NY,NX)=WTNDN(NZ,NY,NX)+WTNDLN(L,NZ,NY,NX) &
      +ZPOOLN(L,NZ,NY,NX)
      WTNDP(NZ,NY,NX)=WTNDP(NZ,NY,NX)+WTNDLP(L,NZ,NY,NX) &
      +PPOOLN(L,NZ,NY,NX)
8920  CONTINUE
      ENDIF
      ENDIF
!
!     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
!
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
!     UPNF=PFT N2 fixation
!     TCUPTK,TZUPTK,TPUPTK=cumulative PFT root-soil C,N,P exchange
!     TZUPFX=cumulative PFT N2 fixation
!
      HCUPTK(NZ,NY,NX)=UPOMC(NZ,NY,NX)
      HZUPTK(NZ,NY,NX)=UPOMN(NZ,NY,NX)+UPNH4(NZ,NY,NX)+UPNO3(NZ,NY,NX) &
      +UPNF(NZ,NY,NX)
      HPUPTK(NZ,NY,NX)=UPOMP(NZ,NY,NX)+UPH2P(NZ,NY,NX)+UPH1P(NZ,NY,NX)
      TCUPTK(NZ,NY,NX)=TCUPTK(NZ,NY,NX)+UPOMC(NZ,NY,NX)
      TZUPTK(NZ,NY,NX)=TZUPTK(NZ,NY,NX)+UPOMN(NZ,NY,NX)+UPNH4(NZ,NY,NX) &
      +UPNO3(NZ,NY,NX)
      TPUPTK(NZ,NY,NX)=TPUPTK(NZ,NY,NX)+UPOMP(NZ,NY,NX)+UPH2P(NZ,NY,NX) &
      +UPH1P(NZ,NY,NX)
      TZUPFX(NZ,NY,NX)=TZUPFX(NZ,NY,NX)+UPNF(NZ,NY,NX)+UPNFC(NZ,NY,NX)
      end subroutine AccumulateStates
!------------------------------------------------------------------------------------------

      subroutine LiterfallByDisturbance(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
!     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
!     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      DO 6375 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*(WTHTR0+WTHTX0) &
      +CFOPC(1,M,NZ,NY,NX)*(WTHTR1+WTHTX1) &
      +CFOPC(2,M,NZ,NY,NX)*(WTHTR2+WTHTX2)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*(WTHNR0+WTHNX0) &
      +CFOPN(1,M,NZ,NY,NX)*(WTHNR1+WTHNX1) &
      +CFOPN(2,M,NZ,NY,NX)*(WTHNR2+WTHNX2)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*(WTHPR0+WTHPX0) &
      +CFOPP(1,M,NZ,NY,NX)*(WTHPR1+WTHPX1) &
      +CFOPP(2,M,NZ,NY,NX)*(WTHPR2+WTHPX2)
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(3,M,NZ,NY,NX)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*(WTHNR3+WTHNX3+WTHNR4+WTHNX4)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*(WTHPR3+WTHPX3+WTHPR4+WTHPX4)
      ELSE
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*(WTHTX3+WTHTX4)
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTHNX3+WTHNX4)
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTHPX3+WTHPX4)
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTR4)*FWOOD(0)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNR4)*FWOODN(0)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPR4)*FWOODP(0)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTR4)*FWOOD(1)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNR4)*FWOODN(1)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPR4)*FWOODP(0)
      ENDIF
6375  CONTINUE
!
!     ABOVE-GROUND LITTERFALL FROM FIRE
!
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
!     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
!     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
!     WTHNL0,WTHPL0=nonstructural N,P to litter
!     WTHNL1,WTHPL1=leaf N,P to litter
!     WTHNL2,WTHPL2=fine,non-leaf N,P to litter
!     WTHNL3,WTHPL3=woody N,P to litter
!     WTHNL4,WTHPL4=standing dead N,P to litter
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
      ELSE
      DO 6485 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(0,M,NZ,NY,NX)*(WTHTR0+WTHTX0) &
      +CFOPC(1,M,NZ,NY,NX)*(WTHTR1+WTHTX1) &
      +CFOPC(2,M,NZ,NY,NX)*(WTHTR2+WTHTX2)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*WTHNL0 &
      +CFOPN(1,M,NZ,NY,NX)*WTHNL1 &
      +CFOPN(2,M,NZ,NY,NX)*WTHNL2
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*WTHPL0 &
      +CFOPP(1,M,NZ,NY,NX)*WTHPL1 &
      +CFOPP(2,M,NZ,NY,NX)*WTHPL2
      ZSNC(4,1,0,NZ,NY,NX)=ZSNC(4,1,0,NZ,NY,NX) &
      +CFOPN(0,M,NZ,NY,NX)*(WTHNR0+WTHNX0-WTHNL0) &
      +CFOPN(1,M,NZ,NY,NX)*(WTHNR1+WTHNX1-WTHNL1) &
      +CFOPN(2,M,NZ,NY,NX)*(WTHNR2+WTHNX2-WTHNL2)
      PSNC(4,1,0,NZ,NY,NX)=PSNC(4,1,0,NZ,NY,NX) &
      +CFOPP(0,M,NZ,NY,NX)*(WTHPR0+WTHPX0-WTHPL0) &
      +CFOPP(1,M,NZ,NY,NX)*(WTHPR1+WTHPX1-WTHPL1) &
      +CFOPP(2,M,NZ,NY,NX)*(WTHPR2+WTHPX2-WTHPL2)
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(3,M,NZ,NY,NX)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*(WTHNL3+WTHNL4)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*(WTHPL3+WTHPL4)
      ZSNC(4,1,0,NZ,NY,NX)=ZSNC(4,1,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*(WTHNR3+WTHNX3-WTHNL3+WTHNR4+WTHNX4-WTHNL4)
      PSNC(4,1,0,NZ,NY,NX)=PSNC(4,1,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*(WTHPR3+WTHPX3-WTHPL3+WTHPR4+WTHPX4-WTHPL4)
      ELSE
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX) &
      +CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTX3)
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*WTHNL3
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*WTHPL3
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX) &
      *CFOPC(3,M,NZ,NY,NX)*(WTHTR4+WTHTX4)*FWOOD(0)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*WTHNL4*FWOODN(0)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*WTHPL4*FWOODP(0)
      ZSNC(4,0,0,NZ,NY,NX)=ZSNC(4,0,0,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNX3-WTHNL3 &
      +WTHNR4+WTHNX4-WTHNL4)*FWOODN(0)
      PSNC(4,0,0,NZ,NY,NX)=PSNC(4,0,0,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPX3-WTHPL3 &
      +WTHPR4+WTHPX4-WTHPL4)*FWOODP(0)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
      +CFOPC(3,M,NZ,NY,NX)*(WTHTR4+WTHTX4)*FWOOD(1)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
      +CFOPN(3,M,NZ,NY,NX)*WTHNL4*FWOODN(1)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
      +CFOPP(3,M,NZ,NY,NX)*WTHPL4*FWOODP(1)
      ZSNC(4,1,0,NZ,NY,NX)=ZSNC(4,1,0,NZ,NY,NX) &
      +CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNX3-WTHNL3 &
      +WTHNR4+WTHNX4-WTHNL4)*FWOODN(1)
      PSNC(4,1,0,NZ,NY,NX)=PSNC(4,1,0,NZ,NY,NX) &
      +CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPX3-WTHPL3 &
      +WTHPR4+WTHPX4-WTHPL4)*FWOODP(1)
      ENDIF
6485  CONTINUE
      ENDIF
      ELSE
!
!     ABOVE-GROUND LITTERFALL FROM GRAZING
!
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!
      TCSNC(NZ,NY,NX)=TCSNC(NZ,NY,NX)+WTHTRT+WTHTXT
      TZSNC(NZ,NY,NX)=TZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
      TPSNC(NZ,NY,NX)=TPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
      TCSN0(NZ,NY,NX)=TCSN0(NZ,NY,NX)+WTHTRT+WTHTXT
      TZSN0(NZ,NY,NX)=TZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
      TPSN0(NZ,NY,NX)=TPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
      ENDIF
      end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

      subroutine TotalBiomRemovalByDisturbance(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!
!     WTHTHT,WTHNHT,WTHPHT=total C,N,P removed
!     WTHTRT,WTHNRT,WTHPRT=total C,N,P to litter
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
      WTHTHT=WTHTH0+WTHTH1+WTHTH2+WTHTH3+WTHTH4
      WTHTRT=WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4
      WTHNHT=WTHNH0+WTHNH1+WTHNH2+WTHNH3+WTHNH4
      WTHNRT=WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4
      WTHPHT=WTHPH0+WTHPH1+WTHPH2+WTHPH3+WTHPH4
      WTHPRT=WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4
      WTHTXT=WTHTX0+WTHTX1+WTHTX2+WTHTX3+WTHTX4
      WTHNXT=WTHNX0+WTHNX1+WTHNX2+WTHNX3+WTHNX4
      WTHPXT=WTHPX0+WTHPX1+WTHPX2+WTHPX3+WTHPX4

      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      IF(JHVST(NZ,I,NY,NX).NE.2)THEN
      HVSTC(NZ,NY,NX)=HVSTC(NZ,NY,NX)+WTHTHT-WTHTRT
      HVSTN(NZ,NY,NX)=HVSTN(NZ,NY,NX)+WTHNHT-WTHNRT
      HVSTP(NZ,NY,NX)=HVSTP(NZ,NY,NX)+WTHPHT-WTHPRT
      TNBP(NY,NX)=TNBP(NY,NX)+WTHTRT-WTHTHT
      XHVSTC(NY,NX)=XHVSTC(NY,NX)+WTHTHT-WTHTRT
      XHVSTN(NY,NX)=XHVSTN(NY,NX)+WTHNHT-WTHNRT
      XHVSTP(NY,NX)=XHVSTP(NY,NX)+WTHPHT-WTHPRT
      ELSE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+WTHTHT-WTHTRT
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+WTHNHT-WTHNRT
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+WTHPHT-WTHPRT
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CNET=PFT net CO2 fixation
!     TNBP=total net biome productivity
!
      ELSE
      VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*(WTHTHT-WTHTRT)
      VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*(WTHTHT-WTHTRT)
      VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*(WTHTHT-WTHTRT)*2.667
      VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-WTHNHT+WTHNRT
      VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
      VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-WTHPHT+WTHPRT
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*(WTHTHT-WTHTRT)
      TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*(WTHTHT-WTHTRT)
!     WRITE(*,5679)'FIRE2',I,J,NZ,VCO2F(NZ,NY,NX),FCH4F,WTHNH0
!    2,WTHNH1,WTHNH2,WTHNH3,WTHNH4,WTHNR0,WTHNR1,WTHNR2
!    3,WTHNR3,WTHNR4,WTHNHT,WTHNRT
!5679  FORMAT(A8,3I4,20E12.4)
      ENDIF
!
!     C,N,P REMOVED FROM GRAZING
!
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     GY=growth yield of grazers
!     WTHTHT,WTHNHT,WTHPHT=total C,N,P removed
!     WTHTRT,WTHNRT,WTHPRT=total C,N,P to litter
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
      ELSE
      HVSTC(NZ,NY,NX)=HVSTC(NZ,NY,NX)+GY*(WTHTHT-WTHTRT)
      HVSTN(NZ,NY,NX)=HVSTN(NZ,NY,NX)+WTHNHT-WTHNRT
      HVSTP(NZ,NY,NX)=HVSTP(NZ,NY,NX)+WTHPHT-WTHPRT
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-GZ*(WTHTHT-WTHTRT)
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-GZ*(WTHTHT-WTHTRT)
!     TNBP(NY,NX)=TNBP(NY,NX)+GY*(WTHTRT-WTHTHT)
!     CNET(NZ,NY,NX)=CNET(NZ,NY,NX)+GZ*(WTHTRT-WTHTHT)
      XHVSTC(NY,NX)=XHVSTC(NY,NX)+GY*(WTHTHT-WTHTRT)
      XHVSTN(NY,NX)=XHVSTN(NY,NX)+WTHNHT-WTHNRT
      XHVSTP(NY,NX)=XHVSTP(NY,NX)+WTHPHT-WTHPRT
      RECO(NY,NX)=RECO(NY,NX)-GZ*(WTHTHT-WTHTRT)
      TRAU(NY,NX)=TRAU(NY,NX)-GZ*(WTHTHT-WTHTRT)
!     WRITE(*,6542)'GRAZ',I,J,NX,NY,NZ,HVSTC(NZ,NY,NX)
!    2,GY,GZ,WTHTHT,WTHTRT
      ENDIF
      end subroutine TotalBiomRemovalByDisturbance
!------------------------------------------------------------------------------------------

      subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,NY,NX)
      implicit none
      integer, intent(in) :: I,J,NZ,NY,NX
!     begin_execution
!     IF NO PLANT C,N,P REMOVED AT HARVEST (ALL RESIDUE RETURNED)
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!     WTHTH1,WTHNH1,WTHPH1=leaf C,N,P removed
!     WTHTH2,WTHNH2,WTHPH2=fine,non-leaf C,N,P removed
!     WTHTH3,WTHNH3,WTHPH3=woody C,N,P removed
!     WTHTH4,WTHNH4,WTHPH4=standing dead C,N,P removed
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.1)THEN
      WTHTR0=WTHTH0
      WTHNR0=WTHNH0
      WTHPR0=WTHPH0
      WTHTR1=WTHTH1
      WTHNR1=WTHNH1
      WTHPR1=WTHPH1
      WTHTR2=WTHTH2-WTHTG*EHVST(2,2,NZ,I,NY,NX)
      WTHNR2=WTHNH2-WTHNG*EHVST(2,2,NZ,I,NY,NX)
      WTHPR2=WTHPH2-WTHPG*EHVST(2,2,NZ,I,NY,NX)
      WTHTR3=WTHTH3
      WTHNR3=WTHNH3
      WTHPR3=WTHPH3
      WTHTR4=WTHTH4
      WTHNR4=WTHNH4
      WTHPR4=WTHPH4
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.2)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.3)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX)*0.5)
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX)*0.5)
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX)*0.5)
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX)*0.5)
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX)*0.5)
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX)*0.5)
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
      FERT(17,I+1,NY,NX)=FERT(17,I+1,NY,NX) &
      +(WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4)/AREA(3,NU(NY,NX),NY,NX)
      FERT(18,I+1,NY,NX)=FERT(18,I+1,NY,NX) &
      +(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA(3,NU(NY,NX),NY,NX)*0.5
      FERT(3,I+1,NY,NX)=FERT(3,I+1,NY,NX) &
      +(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA(3,NU(NY,NX),NY,NX)*0.5
      FERT(19,I+1,NY,NX)=FERT(19,I+1,NY,NX) &
      +(WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4)/AREA(3,NU(NY,NX),NY,NX)
      IYTYP(2,I+1,NY,NX)=3
!     IF(NX.EQ.2)THEN
!     WRITE(*,6542)'MANURE',I,J,NX,NY,NZ,FERT(17,I+1,NY,NX)
!    2,WTHTR1,WTHTR2,WTHTR3,WTHTR4,WTHNR1,WTHNR2,WTHNR3,WTHNR4
!6542  FORMAT(A8,5I4,20E12.4)
!     ENDIF
!
!     REMOVALS BY FIRE
!
!     EFIRE=combustion  of N,P relative to C
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!     WTHTH1,WTHNH1,WTHPH1=leaf C,N,P removed
!     WTHTH2,WTHNH2,WTHPH2=fine,non-leaf C,N,P removed
!     WTHTH3,WTHNH3,WTHPH3=woody C,N,P removed
!     WTHTH4,WTHNH4,WTHPH4=standing dead C,N,P removed
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     WTHTL0,WTHNL0,WTHPL0=nonstructural C,N,P removed from ecosystem
!     WTHTL1,WTHNL1,WTHPL1=leaf C,N,P removed from ecosystem
!     WTHTL2,WTHNL2,WTHPL2=fine,non-leaf C,N,P removed from ecosystem
!     WTHTL3,WTHNL3,WTHPL3=woody C,N,P removed from ecosystem
!     WTHTL4,WTHNL4,WTHPL4=standing dead C,N,P removed from ecosystem
!
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.5)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,1,NZ,I,NY,NX))
      WTHNL0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPL0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,1,NZ,I,NY,NX))
      WTHNL1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPL1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,2,NZ,I,NY,NX))
      WTHNL2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPL2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,3,NZ,I,NY,NX))
      WTHNL3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPL3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX)) &
      *EHVST(2,4,NZ,I,NY,NX))
      WTHNL4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPL4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      ENDIF
      end subroutine ApplyDisturbanceBiomRemoval
!------------------------------------------------------------------------------------------

      subroutine ResetBranchRootStates(NZ,NY,NX)
      implicit none
      integer, intent(in) :: NZ,NY,NX
!     begin_execution
!     RESET BRANCH STATE VARIABLES
!
      DO 8835 NB=1,NBR(NZ,NY,NX)
      CPOOL(NB,NZ,NY,NX)=0._r8
      CPOOLK(NB,NZ,NY,NX)=0._r8
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
      WTSTXB(NB,NZ,NY,NX)=0._r8
      WTSTXN(NB,NZ,NY,NX)=0._r8
      WTSTXP(NB,NZ,NY,NX)=0._r8
8835  CONTINUE
!
!     RESET ROOT STATE VARIABLES
!
      DO 6416 L=NU(NY,NX),NJ(NY,NX)
      DO  N=1,MY(NZ,NY,NX)
      CPOOLR(N,L,NZ,NY,NX)=0._r8
      ZPOOLR(N,L,NZ,NY,NX)=0._r8
      PPOOLR(N,L,NZ,NY,NX)=0._r8
      DO  NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=0._r8
      WTRT1N(N,L,NR,NZ,NY,NX)=0._r8
      WTRT1P(N,L,NR,NZ,NY,NX)=0._r8
      WTRT2(N,L,NR,NZ,NY,NX)=0._r8
      WTRT2N(N,L,NR,NZ,NY,NX)=0._r8
      WTRT2P(N,L,NR,NZ,NY,NX)=0._r8
      RTWT1(N,NR,NZ,NY,NX)=0._r8
      RTWT1N(N,NR,NZ,NY,NX)=0._r8
      RTWT1P(N,NR,NZ,NY,NX)=0._r8
      RTLG1(N,L,NR,NZ,NY,NX)=0._r8
      RTLG2(N,L,NR,NZ,NY,NX)=0._r8
      RTN2(N,L,NR,NZ,NY,NX)=0._r8
      enddo
      enddo
6416  CONTINUE
      WTRVC(NZ,NY,NX)=0._r8
      WTRVN(NZ,NY,NX)=0._r8
      WTRVP(NZ,NY,NX)=0._r8
      IDTH(NZ,NY,NX)=1
      end subroutine ResetBranchRootStates
!------------------------------------------------------------------------------------------

      subroutine ResetDeadRootStates(NB,NZ,NY,NX)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
      implicit none
      integer, intent(in) :: NB,NZ,NY,NX
!     begin_execution
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     GRNOB=seed set number
!     GRNXB=potential number of seed set sites
!     GRWTB=individual seed size
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WSLF=leaf protein mass
!     ARLFB=branch leaf area
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
      CPOOL(NB,NZ,NY,NX)=0._r8
      CPOOLK(NB,NZ,NY,NX)=0._r8
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
      WTSTXB(NB,NZ,NY,NX)=0._r8
      WTSTXN(NB,NZ,NY,NX)=0._r8
      WTSTXP(NB,NZ,NY,NX)=0._r8
      DO 8855 K=0,25
      IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ,NY,NX)=0._r8
      CPOOL4(K,NB,NZ,NY,NX)=0._r8
      CO2B(K,NB,NZ,NY,NX)=0._r8
      HCOB(K,NB,NZ,NY,NX)=0._r8
      ENDIF
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
      DO 8865 L=1,JC
      ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)-ARLFL(L,K,NB,NZ,NY,NX)
      WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)-WGLFL(L,K,NB,NZ,NY,NX)
      ARLFL(L,K,NB,NZ,NY,NX)=0._r8
      WGLFL(L,K,NB,NZ,NY,NX)=0._r8
      WGLFLN(L,K,NB,NZ,NY,NX)=0._r8
      WGLFLP(L,K,NB,NZ,NY,NX)=0._r8
      IF(K.NE.0)THEN
      DO 8860 N=1,4
      SURF(N,L,K,NB,NZ,NY,NX)=0._r8
8860  CONTINUE
      ENDIF
8865  CONTINUE
8855  CONTINUE
      DO 8875 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0._r8
      DO  N=1,4
      SURFB(N,L,NB,NZ,NY,NX)=0._r8
      enddo
8875  CONTINUE

      end subroutine ResetDeadRootStates
!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C4(K,NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: K,NB,NZ,NY,NX

! begin_execution

! FOR EACH CANOPY LAYER
!
  DO 110 L=JC,1,-1
    IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      DO 115 N = 1,4
        DO 120 M = 1,4
!
!         CO2 FIXATION IN MESOPHYLL BY SUNLIT LEAVES
!
!         SURFX=unself-shaded leaf surface area
!
          IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
!
!             C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             ETGR4=light saturated e- transport rate from stomate.f
!             ETLF4=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO4=light-limited PEP carboxylation rate
!             CBXN4=PEP caboxylation efficiency
!             VL=PEP carboxylation rate limited by light,CO2,N,P
!             VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
!             FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
              PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
              EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
              VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
                RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOL(NZ,NY,NX)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(IGTYP(NZ,NY,NX).NE.0)THEN
                  WFN4=RS/RSL
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP carboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               VCGR4,VGROX=PEP carboxylation rate unlimited,limited by CO2
!               XKCO24=Km for VCMX4 from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2I(NZ,NY,NX)
                DO 125 NN=1,100
                  CO2C=CO2X*SCO2(NZ,NY,NX)
                  CO2Y=AMAX1(0.0,CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
                  VGROX=VCGR4(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO24(NZ,NY,NX))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4(K,NB,NZ,NY,NX)
                  VG=(CO2Q(NZ,NY,NX)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)GO TO 130
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Q(NZ,NY,NX)-VA/GSL
                  ELSE
                    VL=0._r8
                    GO TO 130
                  ENDIF
125             CONTINUE
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
130             CH2O4(K)=CH2O4(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX) &
                  *TAUS(L+1,NY,NX)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAUS(L+1,NY,NX))*0.0432
!               IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.3.AND.K.EQ.KLEAF(NB,NZ,NY,NX)
!              2.AND.(I/10)*10.EQ.I.AND.J.EQ.12)THEN
!               WRITE(*,4444)'VLD4',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
!              2,PAR(N,M,L,NZ,NY,NX)*TAUS(L+1,NY,NX)+PARDIF(N,M,L,NZ,NY,NX)
!              3*TAU0(L+1,NY,NX)
!              2,RAPS,TKC(NZ,NY,NX),CO2Q(NZ,NY,NX),ETGR4(K,NB,NZ,NY,NX)
!              3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO
!              3,FDBK4(K,NB,NZ,NY,NX),CH2O4(K),WFN4,VGROX,EGROX
!              4,VCGR4(K,NB,NZ,NY,NX),CO2X,CO2C,CBXNX
!              5,RS,RSL,SURFX(N,L,K,NB,NZ,NY,NX)
!4444            FORMAT(A8,8I4,40E12.4)
!               ENDIF
!
!               C3 CARBOXYLATION REACTIONS IN BUNDLE SHEATH OF C4 PLANTS
!
!               ETGRO=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               CBXN=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!               FDBK=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
                EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
                VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*WFNB*FDBK(NB,NZ,NY,NX)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX) &
                  *TAUS(L+1,NY,NX)
!               IF(L.EQ.NC-1.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
!               WRITE(*,4445)'VLD3',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
!              2,RAPS,TKC(NZ,NY,NX),CO2Q(NZ,NY,NX),ETGRO(K,NB,NZ,NY,NX)
!              3,CBXN(K,NB,NZ,NY,NX),VGRO(K,NB,NZ,NY,NX),EGRO
!              3,FDBK(NB,NZ,NY,NX),WFNB
!4445            FORMAT(A8,8I4,20E12.4)
!               ENDIF
              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
!
!           C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!           QNTM=quantum efficiency
!           PARDIF=diffuse PAR flux
!           ETGR4=light saturated e- transport rate from stomate.f
!           ETLF4=light-limited e- transport rate
!           CURV=shape parameter for e- transport response to PAR
!           EGRO4=light-limited PEP carboxylation rate
!           CBXN4=PEP caboxylation efficiency
!           VL=PEP carboxylation rate limited by light,CO2,N,P
!           VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
!           FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
              PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
              EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
              VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
                RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOL(NZ,NY,NX)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFN4,WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(IGTYP(NZ,NY,NX).NE.0)THEN
                  WFN4=(RS/RSL)**1.00
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP caboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
!               XKCO24=Km for VCMX4 from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2I(NZ,NY,NX)
                DO 135 NN=1,100
                  CO2C=CO2X*SCO2(NZ,NY,NX)
                  CO2Y=AMAX1(0.0,CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
                  VGROX=VCGR4(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO24(NZ,NY,NX))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4(K,NB,NZ,NY,NX)
                  VG=(CO2Q(NZ,NY,NX)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)GO TO 140
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Q(NZ,NY,NX)-VA/GSL
                  ELSE
                    VL=0._r8
                    GO TO 140
                  ENDIF
135             CONTINUE
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0=fraction of diffuse radiation transmitted from layer above
!
140             CH2O4(K)=CH2O4(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX) &
                  *TAU0(L+1,NY,NX)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAU0(L+1,NY,NX))*0.0432
!               WRITE(*,4455)'VLB4',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
!              2,RAPS,TKC(NZ,NY,NX),CO2Q(NZ,NY,NX),ETGR4(K,NB,NZ,NY,NX)
!              3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO
!              3,FDBK4(K,NB,NZ,NY,NX),CH2O4(K),WFN4,VGROX,EGROX
!              4,VCGR4(K,NB,NZ,NY,NX),CO2X,CO2C,CBXNX
!              5,RS,RSL,SURFX(N,L,K,NB,NZ,NY,NX)
!4455            FORMAT(A8,8I4,40E12.4)
!
!               C3 CARBOXYLATION REACTIONS IN IN BUNDLE SHEATH OF C4 PLANTS
!
!               ETGRO=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               CBXN=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!               FDBK=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
                EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
                VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*WFNB*FDBK(NB,NZ,NY,NX)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0=fraction of diffuse radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX) &
                  *TAU0(L+1,NY,NX)
!               IF(J.EQ.13.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
!               WRITE(*,4444)'VLB4',IYRC,I,J,NZ,L,K,VL,PARDIF(N,M,L,NZ,NY,NX)
!              2,RAPY,TKC(NZ,NY,NX),CO2Q(NZ,NY,NX),CO2X,FMOL(NZ,NY,NX)/GSL
!              3,VCGRO(K,NB,NZ,NY,NX),ETLF,FDBK(NB,NZ,NY,NX),WFNB
!               ENDIF
              ENDIF
            ENDIF
          ENDIF
120     CONTINUE
115   CONTINUE
    ENDIF
110   CONTINUE
  CO2F=CO2F+CH2O4(K)
  CH2O=CH2O+CH2O3(K)
  end subroutine ComputeGPP_C4

!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C3(K,NB,NZ,NY,NX)

  implicit none
  integer, intent(in) :: K,NB,NZ,NY,NX

!
! FOR EACH CANOPY LAYER
!
  DO 210 L=JC,1,-1
    IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      DO 215 N=1,4
        DO 220 M=1,4
!
!         CO2 FIXATION BY SUNLIT LEAVES
!
!         SURFX=unself-shaded leaf surface area
!
          IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
!
!             C3 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             ETGRO=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             CBXN=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!             FDBK=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
              PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
              EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
              VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
                RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOL(NZ,NY,NX)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on CO2 fixation
!
                IF(IGTYP(NZ,NY,NX).NE.0)THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMPL=C3 CO2 compensation point (uM)
!               CBXNX=rubisco carboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco
!               VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               XKCO2O=Km for rubisco carboxylation
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2I(NZ,NY,NX)
                DO 225 NN=1,100
                  CO2C=CO2X*SCO2(NZ,NY,NX)
                  CO2Y=AMAX1(0.0,CO2C-COMPL(K,NB,NZ,NY,NX))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPL(K,NB,NZ,NY,NX))
                  VGROX=VCGRO(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO2O(NZ,NY,NX))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*FDBK(NB,NZ,NY,NX)
                  VG=(CO2Q(NZ,NY,NX)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)GO TO 230
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Q(NZ,NY,NX)-VA/GSL
                  ELSE
                    VL=0._r8
                    GO TO 230
                  ENDIF
225             CONTINUE
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
230             CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX) &
                  *TAUS(L+1,NY,NX)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAUS(L+1,NY,NX))*0.0432
!               IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1.AND.K.EQ.KLEAF(NB,NZ,NY,NX)-1
!              2.AND.J.EQ.12)THEN
!               WRITE(20,3335)'VLD',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
!              2,RAPS,TKC(NZ,NY,NX),TKA,CO2Q(NZ,NY,NX),CO2X,CO2C,FMOL(NZ,NY,NX)
!              3/GSL,VGROX,EGROX,ETLF,CBXNX,FDBK(NB,NZ,NY,NX),WFNB,PSILG(NZ,NY,NX)
!              4,VCGRO(K,NB,NZ,NY,NX),ETGRO(K,NB,NZ,NY,NX),COMPL(K,NB,NZ,NY,NX)
!              5,SURFX(N,L,K,NB,NZ,NY,NX),TAUS(L+1,NY,NX),CH2O3(K)
!3335            FORMAT(A8,8I4,30E12.4)
!               ENDIF
              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
!
!             C3 CARBOXYLATION REACTIONS USING VARIABLES FROM 'STOMATE'
!
!             QNTM=quantum efficiency
!             PARDIF=diffuse PAR flux
!             ETGRO=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             CBXN=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!             FDBK=N,P feedback inhibition on C3 CO2 fixation
!
              PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
              PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
              EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
              VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
                RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOL(NZ,NY,NX)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C3 CO2 fixation
!
                IF(IGTYP(NZ,NY,NX).NE.0)THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMPL=C3 CO2 compensation point (uM)
!               CBXNX=rubisco caboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco carboxylase
!               VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               XKCO2O=Km for rubisco carboxylation from stomate.f (uM)
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2I(NZ,NY,NX)
                DO 235 NN=1,100
                  CO2C=CO2X*SCO2(NZ,NY,NX)
                  CO2Y=AMAX1(0.0,CO2C-COMPL(K,NB,NZ,NY,NX))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPL(K,NB,NZ,NY,NX))
                  VGROX=VCGRO(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO2O(NZ,NY,NX))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*FDBK(NB,NZ,NY,NX)
                  VG=(CO2Q(NZ,NY,NX)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)GO TO 240
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Q(NZ,NY,NX)-VA/GSL
                  ELSE
                    VL=0._r8
                    GO TO 240
                  ENDIF
235             CONTINUE
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0=fraction of diffuse radiation transmitted from layer above
!
240             CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX) &
                  *TAU0(L+1,NY,NX)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAU0(L+1,NY,NX))*0.0432
!               IF(J.EQ.13.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
!               WRITE(*,3335)'VLB',IYRC,I,J,NZ,L,K,VL,PARDIF(N,M,L,NZ,NY,NX)
!              2,RAPY,TKC(NZ,NY,NX),CO2Q(NZ,NY,NX),CO2X,FMOL(NZ,NY,NX)/GSL
!              3,VCGRO(K,NB,NZ,NY,NX),ETLF,FDBK(NB,NZ,NY,NX),WFNB
!               ENDIF
              ENDIF
            ENDIF
          ENDIF
220     CONTINUE
215   CONTINUE
    ENDIF
210   CONTINUE
  CO2F=CO2F+CH2O3(K)
  CH2O=CH2O+CH2O3(K)
  end subroutine ComputeGPP_C3

!------------------------------------------------------------------------------------------

  subroutine ComputeRAutoAfEmergence(NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NB,NZ,NY,NX

! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
!
  IF(CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
    CNPG=AMIN1(CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX) &
      +CCPOLB(NB,NZ,NY,NX)*CNKI),CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX) &
      +CCPOLB(NB,NZ,NY,NX)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P
!
! RCO2C=respiration from non-structural C
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! TFN3=temperature function for canopy growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! FDBKX=termination feedback inhibition on C3 CO2
!
  RCO2C=AMAX1(0.0,VMXC*CPOOL(NB,NZ,NY,NX) &
    *TFN3(NZ,NY,NX))*CNPG*FDBKX(NB,NZ,NY,NX)*WFNG
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RMNCS=maintenance respiration
! TFN5=temperature function for canopy maintenance respiration
! WTSHXN=shoot structural N mass
! IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WFNG=growth function of canopy water potential
!
  RMNCS=AMAX1(0.0,RMPLT*TFN5*WTSHXN)
  IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
    RMNCS=RMNCS*WFNG
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X=difference between non-structural C respn and mntc respn
! RCO2Y=growth respiration unlimited by N,P
! WFNSG=expansion,extension function of canopy water potential
! SNCR=excess maintenance respiration
!
  RCO2X=RCO2C-RMNCS
  RCO2Y=AMAX1(0.0,RCO2X)*WFNSG
  SNCR=AMAX1(0.0,-RCO2X)
!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2Y,RCO2G=growth respiration unlimited,limited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! DMSHD=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
!
  IF(RCO2Y.GT.0.0.AND.(CNSHX.GT.0.0.OR.CNLFX.GT.0.0))THEN
    ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
    PPOOLB=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
    RCO2G=AMIN1(RCO2Y,ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
  ELSE
    RCO2G=0._r8
  ENDIF
!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! CGROS=total non-structural C used in growth and growth respiration
! RCO2G=growth respiration limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! ZADDB,PADDB=nonstructural N,P used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! CNRDA=respiration for N assimilation
! CH2O=total CH2O production
!
  CGROS=RCO2G/DMSHD
  ZADDB=AMAX1(0.0,AMIN1(ZPOOL(NB,NZ,NY,NX) &
    ,CGROS*(CNSHX+CNLFM+CNLFX*CNPG)))
  PADDB=AMAX1(0.0,AMIN1(PPOOL(NB,NZ,NY,NX) &
    ,CGROS*(CPSHX+CPLFM+CPLFX*CNPG)))
  CNRDA=AMAX1(0.0,1.70*ZADDB-0.025*CH2O)
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! RCO2T=total C respiration
! RMNCS=maintenance respiration
! RCO2C=respiration from non-structural C
! RCO2G=growth respiration limited by N,P
! SNCR=excess maintenance respiration
! CNRDA=respiration for N assimilation
! CARBN=total PFT CO2 fixation
! CO2F=total CO2 fixation
! TCO2T,TCO2A=total,above-ground PFT respiration
! CNET=PFT net CO2 fixation
! TGPP=ecosystem GPP
! RECO=ecosystem respiration
! TRAU=total autotrophic respiration
!
  RCO2T=AMIN1(RMNCS,RCO2C)+RCO2G+SNCR+CNRDA
  CARBN(NZ,NY,NX)=CARBN(NZ,NY,NX)+CO2F
  TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-RCO2T
  TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-RCO2T
  CNET(NZ,NY,NX)=CNET(NZ,NY,NX)+CO2F-RCO2T
  TGPP(NY,NX)=TGPP(NY,NX)+CO2F
  RECO(NY,NX)=RECO(NY,NX)-RCO2T
  TRAU(NY,NX)=TRAU(NY,NX)-RCO2T
! IF(NZ.EQ.1)THEN
!   WRITE(*,4477)'RCO2',I,J,NX,NY,NZ,NB,IFLGZ,CPOOL(NB,NZ,NY,NX)
!    2,CH2O,RMNCS,RCO2C,CGROS,CNRDA,CNPG,RCO2T,RCO2X,SNCR
!    3,RCO2G,DMSHD,ZADDB,PART(1),PART(2),DMLFB,DMSHB
!    4,WTRSVB(NB,NZ,NY,NX),WVSTKB(NB,NZ,NY,NX),WTSHXN
!    5,ZPOOL(NB,NZ,NY,NX),PPOOL(NB,NZ,NY,NX),PSILT(NZ,NY,NX)
!    6,ZADDB,RNH3B(NB,NZ,NY,NX),WFR(1,NG(NZ,NY,NX),NZ,NY,NX)
!    7,WFNG,TFN3(NZ,NY,NX),TFN5,FDBKX(NB,NZ,NY,NX),VMXC
!4477  FORMAT(A8,7I4,40E12.4)
! ENDIF
!
  end subroutine ComputeRAutoAfEmergence


!------------------------------------------------------------------------------------------

  subroutine ComputeRAutoBfEmergence(NB,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NB,NZ,NY,NX
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  IF(CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
    CNPG=AMIN1(CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX) &
      +CCPOLB(NB,NZ,NY,NX)*CNKI) &
      ,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX) &
      +CCPOLB(NB,NZ,NY,NX)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P, O2 UPTAKE
!
! RCO2CM,RCO2C=respiration from non-structural C unlimited,limited by O2
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! TFN4=temperature function for root growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! FDBKX=termination feedback inhibition on C3 CO2
! WFR=constraint by O2 consumption on all root processes
!
  RCO2CM=AMAX1(0.0,VMXC*CPOOL(NB,NZ,NY,NX) &
    *TFN4(NG(NZ,NY,NX),NZ,NY,NX))*CNPG*WFNG*FDBKX(NB,NZ,NY,NX)
  RCO2C=RCO2CM*WFR(1,NG(NZ,NY,NX),NZ,NY,NX)
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RMNCS=maintenance respiration
! TFN6=temperature function for root maintenance respiration
! WTSHXN=shoot structural N mass
! IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WFNG=growth function of canopy water potential
!
  RMNCS=AMAX1(0.0,RMPLT*TFN6(NG(NZ,NY,NX))*WTSHXN)
  IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
    RMNCS=RMNCS*WFNG
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! WFNSG=expansion,extension function of canopy water potential
! SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!
  RCO2XM=RCO2CM-RMNCS
  RCO2X=RCO2C-RMNCS
  RCO2YM=AMAX1(0.0,RCO2XM)*WFNSG
  RCO2Y=AMAX1(0.0,RCO2X)*WFNSG
  SNCRM=AMAX1(0.0,-RCO2XM)
  SNCR=AMAX1(0.0,-RCO2X)
!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! FNP=growth respiration limited by O2 and N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! RCO2GM,RCO2G=growth respiration unltd,ltd by O2 and limited by N,P
! WFR=constraint by O2 consumption on all root processes
!
  IF(CNSHX.GT.0.0.OR.CNLFX.GT.0.0)THEN
    ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
    PPOOLB=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
    FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
    IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
    ELSE
      RCO2GM=0._r8
    ENDIF
    IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(1,NG(NZ,NY,NX),NZ,NY,NX))
    ELSE
      RCO2G=0._r8
    ENDIF
  ELSE
    RCO2GM=0._r8
    RCO2G=0._r8
  ENDIF
!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! CGROSM,CGROS=total non-structural C used in growth and respn unltd,ltd by O2
! RCO2GM,RCO2G=growth respiration unltd,ltd by O2 and limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! ZADDBM,ZADDB,PADDB=nonstructural N,P unltd,ltd by O2 used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
  CGROSM=RCO2GM/DMSHD
  CGROS=RCO2G/DMSHD
  ZADDBM=AMAX1(0.0,CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
  ZADDB=AMAX1(0.0,CGROS*(CNSHX+CNLFM+CNLFX*CNPG))
  PADDB=AMAX1(0.0,CGROS*(CPSHX+CPLFM+CPLFX*CNPG))
  CNRDM=AMAX1(0.0,1.70*ZADDBM)
  CNRDA=AMAX1(0.0,1.70*ZADDB)
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! RCO2TM,RCO2T=total C respiration unltd,ltd by O2
! RMNCS=maintenance respiration
! RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
! SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
! CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
! RCO2A=total root respiration
! RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
  RCO2TM=RMNCS+RCO2GM+SNCRM+CNRDM
  RCO2T=RMNCS+RCO2G+SNCR+CNRDA
  RCO2M(1,NG(NZ,NY,NX),NZ,NY,NX)=RCO2M(1,NG(NZ,NY,NX),NZ,NY,NX)+RCO2TM
  RCO2N(1,NG(NZ,NY,NX),NZ,NY,NX)=RCO2N(1,NG(NZ,NY,NX),NZ,NY,NX)+RCO2T
  RCO2A(1,NG(NZ,NY,NX),NZ,NY,NX)=RCO2A(1,NG(NZ,NY,NX),NZ,NY,NX)-RCO2T
  CH2O=0._r8
  end subroutine ComputeRAutoBfEmergence

end module grosubMod
