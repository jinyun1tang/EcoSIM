module GrosubPars

! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = __FILE__
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
!
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth (g N,P g-1 C)
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     PSILM=minimum water potential for organ expansion,extension (MPa)
!     RCMN=minimum stomatal resistance to CO2 (s m-1)
!     RTDPX=distance behind growing point for secondary roots (m)
!     RTLGAX=minimum average secondary root length (m)
!     EMODR=root modulus of elasticity (MPa)
!
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
!
!     ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     GY,GZ=partitioning of grazed material to removal,respiration
!     FSTK=fraction of stalk area contributing to water,heat flow
!     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
!     FRTX=fraction used to calculate woody faction of stalk,root
!
!
!     SETC,SETN,SETP=Km for nonstructural C,N,P concn on seed set (g g-1)
!     SLA2,SSL2,SNL2=parameter for calculating leaf area expansion, petiole
!     and internode extension vs leaf, petiole, internode growth
!     CNMX,CPMX,CNMN,CPMN=max,min N:C,P:C for nonstructural C,N,P transfers
!
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
  real(r8) :: PART1X
  real(r8) :: PART2X
  real(r8) :: VMXC
  real(r8) :: FSNR
  real(r8) :: FLG4X
  real(r8) :: FLGZX
  real(r8) :: XFRX
  real(r8) :: XFRY
  real(r8) :: FSNK
  real(r8) :: FXFS
  real(r8) :: FMYC
  real(r8) :: CNKI
  real(r8) :: CPKI
  real(r8) :: RMPLT
  real(r8) :: PSILM
  real(r8) :: RCMN
  real(r8) :: RTDPX
  real(r8) :: RTLGAX
  real(r8) :: EMODR
  real(r8) :: QNTM
  real(r8) :: CURV
  real(r8) :: CURV2
  real(r8) :: CURV4
  real(r8) :: ELEC3
  real(r8) :: ELEC4
  real(r8) :: CO2KI
  real(r8) :: FCO2B
  real(r8) :: FHCOB
  real(r8) :: COMP4
  real(r8) :: FDML
  real(r8) :: FBS
  real(r8) :: FMP
  real(r8) :: ZPLFM
  real(r8) :: ZPLFD
  real(r8) :: ZPGRM
  real(r8) :: ZPGRD
  real(r8) :: GY
  real(r8) :: GZ
  real(r8) :: FSTK
  real(r8) :: ZSTX
  real(r8) :: DSTK
  real(r8) :: VSTK
  real(r8) :: FRTX
  real(r8) :: SETC
  real(r8) :: SETN
  real(r8) :: SETP
  real(r8) :: SLA2
  real(r8) :: SSL2
  real(r8) :: SNL2
  real(r8) :: CNMX
  real(r8) :: CPMX
  real(r8) :: CNMN
  real(r8) :: CPMN
  real(r8) :: EN2F
  real(r8) :: VMXO
  real(r8) :: SPNDLK
  real(r8) :: SPNDL
  real(r8) :: CCNGR
  real(r8) :: CCNGB
  real(r8) :: WTNDI
  real(r8) :: CZKM
  real(r8) :: CPKM
  real(r8) :: RCCZR
  real(r8) :: RCCYR
  real(r8) :: RCCXR
  real(r8) :: RCCQR
  real(r8) :: RCCZN
  real(r8) :: RCCYN
  real(r8) :: RCCXN
  real(r8) :: RCCQN

  integer :: IFLGQX
  REAL(R8) :: FRSV(0:3),FXFY(0:1),FXFZ(0:1)
  real(r8) :: FXFB(0:3),FXFR(0:3),FXRT(0:1),FXSH(0:1),FXRN(6)
  REAL(R8) :: RCCX(0:3),RCCQ(0:3)
  REAL(R8) :: RCCZ(0:3),RCCY(0:3)
  real(r8) :: FLG4Y(0:3)
  real(r8) :: ATRPX(0:1)
  real(r8) :: GVMX(0:1)
  real(r8) :: RTSK(0:3)

  type, public :: plant_bgc_par_type
   !nonstructural(0,*),
   !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
  integer :: instruct
  integer :: ifoliar
  integer :: infoliar
  integer :: istalk
  integer :: iroot
  integer :: icwood
  integer :: jpstgs       !number of growth stages
  integer :: JRS          !maximum number of root layers
  integer  :: JP1         !number of plants
  integer  :: JBR         !maximum number of branches
  integer  :: JSA1        !number of sectors for the sky azimuth  [0,2*pi]
  integer  :: jcplx       !number of organo-microbial complexes
  integer  :: JLA1        !number of sectors for the leaf azimuth, [0,pi]
  integer  :: JC1         !number of canopy layers
  integer  :: JZ1         !number of soil layers
  integer  :: JLI1        !number of sectors for the leaf zenith [0,pi/2]
  integer  :: JNODS1      !number of canopy nodes
  integer  :: jsken       !number of kinetic components in litter,  PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer  :: Jlitgrp     !number of litter groups nonstructural(0,*),
                          !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
  integer  :: JPRT
  integer  :: n_pltlitrk  !number of plant litter microbial-om complexes
  integer :: iprotein
  integer :: icarbhyro
  integer :: icellulos
  integer :: ilignin
  integer :: k_woody_litr
  integer :: k_fine_litr
  integer :: jroots
  end type plant_bgc_par_type


  contains

  subroutine InitVegPars(pltpar)
  implicit none
  type(plant_bgc_par_type)  :: pltpar

  pltpar%instruct=0
  pltpar%ifoliar=1
  pltpar%infoliar=2
  pltpar%istalk=3
  pltpar%iroot=4
  pltpar%icwood=5
  pltpar%jpstgs=10
  pltpar%JRS=10
  pltpar%JPRT=7
  pltpar%jroots=2

  PART1X=0.05_r8
  PART2X=0.02_r8
  VMXC=0.015_r8
  FSNR=2.884E-03_r8
  FLG4X=168.0_r8
  FLGZX=240.0_r8
  XFRX=2.5E-02_r8
  XFRY=2.5E-03_r8
  FSNK=0.05_r8
  FXFS=1.0_r8
  FMYC=0.1_r8
  CNKI=1.0E-01_r8
  CPKI=1.0E-02_r8
  RMPLT=0.010_r8
  PSILM=0.1_r8
  RCMN=1.560E+01_r8
  RTDPX=0.00_r8
  RTLGAX=1.0E-03_r8
  EMODR=5.0_r8
  QNTM=0.45_r8
  CURV=0.70_r8
  CURV2=2.0_r8*CURV
  CURV4=4.0_r8*CURV
  ELEC3=4.5_r8
  ELEC4=3.0_r8
  CO2KI=1.0E+03_r8
  FCO2B=0.02_r8
  FHCOB=1.0_r8-FCO2B
  COMP4=0.5_r8
  FDML=6.0_r8
  FBS=0.2_r8*FDML
  FMP=0.8_r8*FDML
  ZPLFM=0.33_r8
  ZPLFD=1.0_r8-ZPLFM
  ZPGRM=0.75_r8
  ZPGRD=1.0_r8-ZPGRM
  GY=1.0_r8
  GZ=1.0_r8-GY
  FSTK=0.05_r8
  ZSTX=1.0E-03_r8
  DSTK=0.225_r8
  VSTK=1.0E-06_r8/DSTK
  FRTX=1.0_r8/(1.0_r8-(1.0_r8-FSTK)**2)
  SETC=1.0E-02_r8
  SETN=1.0E-03_r8
  SETP=1.0E-04_r8
  SLA2=-0.33_r8
  SSL2=-0.50_r8
  SNL2=-0.67_r8
  CNMX=0.20_r8
  CPMX=0.020_r8
  CNMN=0.050_r8
  CPMN=0.005_r8
  EN2F=0.20_r8
  VMXO=0.125_r8
  SPNDLK=0.01_r8
  SPNDL=5.0E-04_r8
  CCNGR=2.5E-01_r8
  CCNGB=6.0E-04_r8
  WTNDI=1.0E-03_r8
  CZKM=2.5E-03_r8
  CPKM=2.5E-04_r8
  RCCZR=0.056_r8
  RCCYR=0.167_r8
  RCCXR=0.833_r8
  RCCQR=0.833_r8
  RCCZN=0.167_r8
  RCCYN=0.833_r8
  RCCXN=0.833_r8
  RCCQN=0.833_r8

  IFLGQX=960


  RCCZ=real((/0.167,0.167,0.167,0.056/),r8)
  RCCY=real((/0.333,0.333,0.167,0.333/),r8)
  RCCX=real((/0.417,0.833,0.833,0.833/),r8)
  RCCQ=real((/0.417,0.833,0.833,0.833/),r8)

  RTSK=real((/0.25,1.0,4.0,4.0/),r8)
  FXRN=real((/0.25,0.125,0.0625,0.225,0.075,0.025/),r8)
  FXFB=real((/1.0E-02,1.0E-02,1.0E-05,5.0E-05/),r8)
  FXFR=real((/1.0E-02,1.0E-02,1.0E-05,5.0E-05/),r8)
  FXSH=real((/0.50,0.75/),r8);FXRT=real((/0.50,0.25/),r8)
  FRSV=real((/0.025,0.025,0.001,0.001/),r8)
  FXFY=real((/0.025,0.005/),r8);FXFZ=real((/0.25,0.05/),r8)
  FLG4Y=real((/360.0,1440.0,720.0,720.0/),r8)
  ATRPX=real((/68.96,276.9/),r8);GVMX=real((/0.010,0.0025/),r8)
  end subroutine InitVegPars


end module GrosubPars
