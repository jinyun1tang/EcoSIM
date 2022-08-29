module UptakesMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,test_aneb
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use NutUptakeMod
  use PlantAPIData
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__

  public :: uptakes
  public :: InitUptake
  contains

  subroutine InitUptake

  implicit none

  call InitUptakePars

  end subroutine InitUptake
!------------------------------------------------------------------------------------------

  subroutine uptakes(I,J)
!
!     THIS subroutine CALCULATES EXCHANGES OF ENERGY, C, N AND P
!     BETWEEN THE CANOPY AND THE ATMOSPHERE AND BETWEEN ROOTS AND THE SOIL
!
  implicit none
  integer, intent(in) :: I, J

  integer :: NN,N,NZ,K,L
  real(r8) :: FPC
  real(r8) :: OSTRN,OSTRD,PARHC
  real(r8) :: WVPLT
  real(r8) :: PSIST1(JZ1),PATH(2,JZ1)
  real(r8) :: RRADL(2,JZ1),RSRT(2,JZ1)
  real(r8) :: RSRG(2,JZ1),RSR1(2,JZ1)
  real(r8) :: RSR2(2,JZ1),RSSX(2,JZ1)
  real(r8) :: RSRS(2,JZ1),WTRTG(JZ1)
  real(r8) :: FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8) :: FRTDPX(JZ1,JP1),RTARR(2,JZ1)
  real(r8) :: VOLPU(JZ1),VOLWU(JZ1)
  real(r8) :: TKCX,CNDT,HFLWC1,VHCPX
  real(r8) :: DIFF,PSILH,UPRT,VFLXC
  real(r8) :: FDMP
  integer :: ILYR(2,JZ1)
!     begin_execution
  associate(                         &
    RA     => plt_photo%RA     , &
    VOLWP  => plt_ew%VOLWP     , &
    TCC    => plt_ew%TCC       , &
    TKCZ   => plt_ew%TKCZ      , &
    FLWC   => plt_ew%FLWC      , &
    EFLXC  => plt_ew%EFLXC     , &
    TKC    => plt_ew%TKC       , &
    PSILT  => plt_ew%PSILT     , &
    SFLXC  => plt_ew%SFLXC     , &
    TSHC   => plt_ew%TSHC      , &
    DTKC   => plt_ew%DTKC      , &
    EVAPC  => plt_ew%EVAPC     , &
    TKA    => plt_ew%TKA       , &
    TLEC   => plt_ew%TLEC      , &
    EP     => plt_ew%EP        , &
    VOLWC  => plt_ew%VOLWC     , &
    CDPTHZ => plt_site%CDPTHZ  , &
    PPT    => plt_site%PPT     , &
    NP     => plt_site%NP      , &
    PP     => plt_site%PP      , &
    NU     => plt_site%NU      , &
    AREA3  => plt_site%AREA3   , &
    ZEROS  => plt_site%ZEROS   , &
    WTRTD  => plt_biom%WTRTD   , &
    ZEROL  => plt_biom%ZEROL   , &
    ZEROP  => plt_biom%ZEROP   , &
    WTLS   => plt_biom%WTLS    , &
    WVSTK  => plt_biom%WVSTK   , &
    IDAY   => plt_pheno%IDAY   , &
    OSTR   => plt_pheno%OSTR   , &
    IFLGC  => plt_pheno%IFLGC  , &
    ARLFC  => plt_morph%ARLFC  , &
    RTDP1  => plt_morph%RTDP1  , &
    ARLSS  => plt_morph%ARLSS  , &
    ARLFS  => plt_morph%ARLFS  , &
    SDPTH  => plt_morph%SDPTH  , &
    NB1    => plt_morph%NB1    , &
    FRADP  => plt_rad%FRADP      &
  )

  call PrepUptake(PSIST1,WTRTG,VOLPU,VOLWU)
!
!     IF PLANT SPECIES EXISTS
!
  DO NZ=1,NP
    OSTRN=0.0_r8
    OSTRD=0.0_r8
    IF(IFLGC(NZ).EQ.1.AND.PP(NZ).GT.0.0)THEN

      call UpdateCanopyProperty(NZ)

!     STOMATE=solve for minimum canopy stomatal resistance
      CALL STOMATEs(I,J,NZ)
!
!     CALCULATE VARIABLES USED IN ROOT UPTAKE OF WATER AND NUTRIENTS
      call UpdateRootProperty(NZ,PATH,RRADL,WTRTG,FPQ,FPP,FRTDPX,RTARR)
!
!     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
!     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER CONTENT
!
!     (AG: - originally this line had a N0B1 here )
      IF((IDAY(1,NB1(NZ),NZ).NE.0).AND.(ARLFS(NZ).GT.ZEROL(NZ) &
        .AND.FRADP(NZ).GT.0.0_r8).AND.(RTDP1(1,1,NZ).GT.SDPTH(NZ)+CDPTHZ(0)))THEN
!
        call CalcResistance(NZ,PATH,RRADL,RTARR,RSRT,RSRG,RSR1,RSR2,RSSX,RSRS,CNDT,PSILH,ILYR)
!
!     INITIALIZE CANOPY WATER POTENTIAL, OTHER VARIABLES USED IN ENERGY
!     BALANCE THAT DON'T NEED TO BE RECALCULATED DURING CONVERGENCE
!
!     PSILT=initial estimate of total canopy water potential
!     FLWC,HFLWC1=convective water,heat flux from precip to canopy
!     FTHRM,FDTHS=LW emitted,absorbed by canopy
!     FPC=fraction of ground area AREA occupied by PFT
!     CCPOLT=total nonstructural canopy C,N,P concentration
!     CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy
!     OSWT=molar mass of CCPOLT
!     TKCX=intermediate estimate of TKC used in convergence solution
!     WVPLT=leaf+petiole+stalk mass
!     VSTK=specific stalk volume
!     VHCPX=canopy heat capacity
!     VOLWP,VOLWC=water volume in canopy,on canopy surfaces
!
        PSILT(NZ)=AMIN1(-1.0E-06,0.667*PSILT(NZ))
        EP(NZ)=0.0_r8
        EVAPC(NZ)=0.0_r8
        HFLWC1=FLWC(NZ)*cpw*TKA

        IF(ARLSS.GT.ZEROS)THEN
          FPC=ARLFS(NZ)/ARLSS*AMIN1(1.0,0.5*ARLFC/AREA3(NU))
        ELSEIF(PPT.GT.ZEROS)THEN
          FPC=PP(NZ)/PPT
        ELSE
          FPC=1.0/NP
        ENDIF

        TKCX=TKC(NZ)
        WVPLT=AMAX1(0.0,WTLS(NZ)+WVSTK(NZ))
        VHCPX=cpw*(WVPLT*VSTK+VOLWC(NZ)+VOLWP(NZ))
!
!     CONVERGENCE SOLUTION
!
        NN=CanopyEnergyH2OIteration(I,J,NZ,FPC,WVPLT,&
          PSIST1,PARHC,DIFF,UPRT,VFLXC,FDMP,RSRS,FPQ,VOLPU,VOLWU,TKCX,CNDT,VHCPX,&
          HFLWC1,PSILH,ILYR)
!
!     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
!
!     TKC=final estimate of canopy temperature TKCZ
!     TKA=current air temperature
!     DTKC=TKC-TKA for next hour
!
        TKC(NZ)=TKCZ(NZ)
        TCC(NZ)=TKC(NZ)-TC2K
        DTKC(NZ)=TKC(NZ)-TKA
!
!     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
!     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
!
        call HandlingDivergence(I,J,NN,NZ,PSIST1,DIFF,FDMP)

        call UpdateCanopyWater(NZ,PARHC,PSIST1,RSRT,RSSX,RSRS,TKCX,VHCPX,HFLWC1,UPRT,VFLXC,ILYR)
!
!     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
!
      ELSE
        call HandleBareSoil(NZ,PSIST1,FDMP)
      ENDIF

      call SetCanopyGrowthFuncs(NZ)


      call NutO2Uptake(NZ,FDMP,OSTRN,OSTRD,PATH,RRADL,&
        FPQ,FPP,FRTDPX,RTARR)

      TLEC=TLEC+EFLXC(NZ)*RA(NZ)
      TSHC=TSHC+SFLXC(NZ)*RA(NZ)
      IF(OSTRD.GT.ZEROP(NZ))THEN
        OSTR(NZ)=OSTRN/OSTRD
      ELSE
        OSTR(NZ)=0.0_r8
      ENDIF
    ENDIF
  ENDDO
  RETURN
  end associate
  END subroutine uptakes
!------------------------------------------------------------------------

  subroutine PrepUptake(PSIST1,WTRTG,VOLPU,VOLWU)
!
!     prepare for uptake calculation
  implicit none
  real(r8), intent(out) :: PSIST1(JZ1),WTRTG(JZ1),VOLPU(JZ1),VOLWU(JZ1)
  integer :: NZ, L, N
  real(r8) :: ARLSC

  associate(                          &
    ZERO   => plt_site%ZERO     , &
    NP     => plt_site%NP       , &
    NJ     => plt_site%NJ       , &
    VOLWM  => plt_site%VOLWM    , &
    ALT    => plt_site%ALT      , &
    NU     => plt_site%NU       , &
    NP0    => plt_site%NP0      , &
    PSIST  => plt_ew%PSIST      , &
    WTRTD  => plt_biom%WTRTD    , &
    VOLA   => plt_soilchem%VOLA , &
    VOLI   => plt_soilchem%VOLI , &
    THETY  => plt_soilchem%THETY, &
    BKDS   => plt_soilchem%BKDS , &
    VOLW   => plt_soilchem%VOLW , &
    VOLY   => plt_soilchem%VOLY , &
    ARSTP  => plt_morph%ARSTP   , &
    MY     => plt_morph%MY      , &
    ARLFP  => plt_morph%ARLFP   , &
    RAD1   => plt_rad%RAD1      , &
    THRM1  => plt_rad%THRM1       &
  )
!
!     RESET TOTAL UPTAKE ARRAYS
!
!     ARLFP,ARSTP=leaf,stalk areas
!

  ARLSC=0.0_r8
  DO 9984 NZ=1,NP0
!     TKC(NZ)=TKA+DTKC(NZ)
!     TCC(NZ)=TKC(NZ)-TC2K
    ARLSC=ARLSC+ARLFP(NZ)+ARSTP(NZ)
    RAD1(NZ)=0.0_r8
    plt_ew%EFLXC(NZ)=0.0_r8
    plt_ew%SFLXC(NZ)=0.0_r8
    plt_ew%HFLXC(NZ)=0.0_r8
    THRM1(NZ)=0.0_r8
    plt_ew%EP(NZ)=0.0_r8
    plt_ew%EVAPC(NZ)=0.0_r8
    plt_rbgc%UPOMC(NZ)=0.0_r8
    plt_rbgc%UPOMN(NZ)=0.0_r8
    plt_rbgc%UPOMP(NZ)=0.0_r8
    plt_rbgc%UPNH4(NZ)=0.0_r8
    plt_rbgc%UPNO3(NZ)=0.0_r8
    plt_rbgc%UPH2P(NZ)=0.0_r8
    plt_rbgc%UPH1P(NZ)=0.0_r8
    plt_rbgc%UPNF(NZ)=0.0_r8
!
!     RESET UPTAKE ARRAYS
!
    DO  L=NU,NJ
      DO  N=1,MY(NZ)
        plt_ew%UPWTR(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2P(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXP(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXS(N,L,NZ)=0.0_r8
        plt_rbgc%RUPCHS(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN2S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3B(N,L,NZ)=0.0_r8
        plt_rbgc%RUPHGS(N,L,NZ)=0.0_r8
        plt_rbgc%RCOFLA(N,L,NZ)=0.0_r8
        plt_rbgc%ROXFLA(N,L,NZ)=0.0_r8
        plt_rbgc%RCHFLA(N,L,NZ)=0.0_r8
        plt_rbgc%RN2FLA(N,L,NZ)=0.0_r8
        plt_rbgc%RNHFLA(N,L,NZ)=0.0_r8
        plt_rbgc%RHGFLA(N,L,NZ)=0.0_r8
        plt_rbgc%RCODFA(N,L,NZ)=0.0_r8
        plt_rbgc%ROXDFA(N,L,NZ)=0.0_r8
        plt_rbgc%RCHDFA(N,L,NZ)=0.0_r8
        plt_rbgc%RN2DFA(N,L,NZ)=0.0_r8
        plt_rbgc%RNHDFA(N,L,NZ)=0.0_r8
        plt_rbgc%RHGDFA(N,L,NZ)=0.0_r8
      enddo
    enddo
9984  CONTINUE
!
!     PSIST1=total soil water potential PSIST adjusted for surf elevn
!     ALT=surface elevation
!     VOLWU,VOLWM=water volume available for uptake,total water volume
!     THETY,VOLX=hygroscopic SWC,soil volume
!     VOLPU=air volume
!     WTRTG=total biome root mass
!
  DO 9000 L=NU,NJ
    PSIST1(L)=PSIST(L)-0.0098_r8*ALT
    IF(BKDS(L).GT.ZERO)THEN
      VOLWU(L)=VOLWM(NPH,L)-THETY(L)*VOLY(L)
      VOLPU(L)=AMAX1(0.0,VOLA(L)-VOLW(L)-VOLI(L))
    ELSE
      VOLWU(L)=VOLWM(NPH,L)
      VOLPU(L)=0.0_r8
    ENDIF
    WTRTG(L)=0.0_r8
    DO 9005 NZ=1,NP
      DO  N=1,MY(NZ)
!     IF(IFLGC(NZ).EQ.1.AND.PP(NZ).GT.0.0)THEN
      WTRTG(L)=WTRTG(L)+AMAX1(0.0,WTRTD(N,L,NZ))
!     ENDIF
      enddo
9005  CONTINUE
9000  CONTINUE
  end associate
  end subroutine PrepUptake
!------------------------------------------------------------------------
  subroutine UpdateCanopyProperty(NZ)
!
!     update canopy characterization
  implicit none
  integer, intent(in) :: NZ
  real(r8) :: ALFZ
  real(r8) :: TFRADP,RACZ(JP1)
  integer :: NB,K,L,N,NZZ

  associate(                          &
    ZERO   =>  plt_site%ZERO    , &
    NP     =>  plt_site%NP      , &
    IETYP  =>  plt_site%IETYP   , &
    RAB    =>  plt_ew%RAB       , &
    TKA    =>  plt_ew%TKA       , &
    DTKC   =>  plt_ew%DTKC      , &
    ZD     =>  plt_ew%ZD        , &
    ZR     =>  plt_ew%ZR        , &
    RAZ    =>  plt_ew%RAZ       , &
    TKCZ   =>  plt_ew%TKCZ      , &
    WSLF   =>  plt_biom%WSLF    , &
    ZEROP  =>  plt_biom%ZEROP   , &
    FRADP  =>  plt_rad%FRADP    , &
    SURFX  =>  plt_photo%SURFX  , &
    ARLF1  =>  plt_morph%ARLF1  , &
    KLEAFX =>  plt_morph%KLEAFX , &
    ZC     =>  plt_morph%ZC     , &
    CFX    =>  plt_morph%CFX    , &
    ARLFS  =>  plt_morph%ARLFS  , &
    NBR    =>  plt_morph%NBR    , &
    SURF   =>  plt_morph%SURF   , &
    ZT     =>  plt_morph%ZT       &
  )
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!
  DO 500 NB=1,NBR(NZ)
    DO 550 K=1,JNODS1
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     ARLF=leaf area
!     WSLF=leaf protein content
!     SURFX,SURF=unself-shaded,total leaf surface area
!     CFX=clumping factor from PFT file
!
      IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ) &
        .AND.WSLF(K,NB,NZ).GT.ZEROP(NZ))THEN
        KLEAFX(NB,NZ)=K
      ENDIF
      DO 600 L=JC1,1,-1
        DO 650 N=1,JLI1
          SURFX(N,L,K,NB,NZ)=SURF(N,L,K,NB,NZ)*CFX(NZ)
650     CONTINUE
600   CONTINUE
550 CONTINUE
500 CONTINUE
!
!     CANOPY HEIGHT FROM HEIGHT OF MAXIMUM LEAF LAYER
!
!     DIFFUSIVE RESISTANCE OF OTHER TALLER CANOPIES TO HEAT AND VAPOR
!     TRANSFER OF CURRENT CANOPY ADDED TO BOUNDARY LAYER RESISTANCE
!     OF TALLEST CANOPY CALCULATED IN 'HOUR1'
!
!     IETYP=Koppen climate zone
!     ZC,ZT,ZR=PFT canopy,biome,surface roughness height
!     FRADP=fraction of radiation received by each PFT canopy
!     ALFZ=shape parameter for windspeed attenuation in canopy
!     RAB=biome canopy isothermal boundary layer resistance
!     RACZ,RAZ=additional,total PFT canopy isothermal blr
!
  IF(ARLFS(NZ).GT.0.0_r8)THEN
    IF(IETYP.GE.0)THEN
      TFRADP=0.0_r8
      DO 700 NZZ=1,NP
        IF(ZC(NZZ).GT.ZC(NZ)+ZR)THEN
          TFRADP=TFRADP+FRADP(NZZ)
        ENDIF
700   CONTINUE
      ALFZ=2.0_r8*TFRADP
      IF(RAB.GT.ZERO.AND.ZT.GT.ZERO.AND.ALFZ.GT.ZERO)THEN
        RACZ(NZ)=AMIN1(RACX,AMAX1(0.0,ZT*EXP(ALFZ) &
          /(ALFZ/RAB)*(EXP(-ALFZ*ZC(NZ)/ZT) &
          -EXP(-ALFZ*(ZD+ZR)/ZT))))
      ELSE
        RACZ(NZ)=0.0_r8
      ENDIF
    ELSE
      RACZ(NZ)=0.0_r8
    ENDIF
  ELSE
    RACZ(NZ)=RACX
  ENDIF
  RAZ(NZ)=RAB+RACZ(NZ)
!
!     INITIALIZE CANOPY TEMPERATURE WITH CURRENT AIR TEMPERATURE AND
!     LAST HOUR'S CANOPY-AIR TEMPERATURE DIFFERENCE, AND CALL A
!     subroutine TO CALCULATE MINIMUM CANOPY STOMATAL RESISTANCE
!     FOR SUBSEQUENT USE IN ENERGY EXCHANGE CALCULATIONS
!
!     TKCZ=initial estimate of canopy temperature TKC
!     TKA=current air temperature
!     DTKC=TKC-TKA from previous hour
!
  TKCZ(NZ)=TKA+DTKC(NZ)
  end associate
  end subroutine UpdateCanopyProperty
!------------------------------------------------------------------------
  subroutine UpdateRootProperty(NZ,PATH,RRADL,WTRTG,FPQ,FPP,FRTDPX,RTARR)
!
!     update root characterization

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: WTRTG(JZ1)
  real(r8), intent(out) :: PATH(2,JZ1),RRADL(2,JZ1)
  real(r8), intent(out) :: FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(out) :: FRTDPX(JZ1,JP1),RTARR(2,JZ1)
  real(r8) :: RTDPZ,RTDPX
  integer :: N,L,NR

  associate(                           &
    WTRTD  =>  plt_biom%WTRTD    , &
    FMPR   =>  plt_site%FMPR     , &
    CDPTHZ =>  plt_site%CDPTHZ   , &
    DLYR3  =>  plt_site%DLYR3    , &
    ZEROS  =>  plt_site%ZEROS    , &
    ZERO   =>  plt_site%ZERO     , &
    PP     =>  plt_site%PP       , &
    NU     =>  plt_site%NU       , &
    NRT    =>  plt_morph%NRT     , &
    MY     =>  plt_morph%MY      , &
    RTDNP  =>  plt_morph%RTDNP   , &
    RRAD2X =>  plt_morph%RRAD2X  , &
    HTCTL  =>  plt_morph%HTCTL   , &
    RTLGP  =>  plt_morph%RTLGP   , &
    RTDP1  =>  plt_morph%RTDP1   , &
    PORT   =>  plt_morph%PORT    , &
    RTVLW  =>  plt_morph%RTVLW   , &
    RRAD2M =>  plt_morph%RRAD2M  , &
    SDPTH  =>  plt_morph%SDPTH   , &
    NI     =>  plt_morph%NI        &
  )
!     RTDPZ,RTDP1=primary root depth
!     FRTDPX=fraction of each soil layer with primary root
!     DLYR=layer thickness
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!     HTCTL=hypocotyledon height
!     FPQ=PFT fraction of biome root mass
!     RTDNP,RTLGP=root length density,root length per plant
!     RRADL=root radius
!     PORT=root porosity
!     FMPR=micropore fraction excluding macropore,rock
!     RTVLW=root aqueous volume
!     PP=plant population
!     PATH=path length of water and nutrient uptake
!     RTARR=root surface area/radius for uptake,diffusivity
!
  DO 2000 N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(N.EQ.1)THEN
        RTDPZ=0.0_r8
        DO 2005 NR=1,NRT(NZ)
          RTDPZ=AMAX1(RTDPZ,RTDP1(1,NR,NZ))
2005    CONTINUE
        IF(L.EQ.NU)THEN
          FRTDPX(L,NZ)=1.0_r8
        ELSE
          IF(DLYR3(L).GT.ZERO)THEN
            RTDPX=AMAX1(0.0_r8,RTDPZ-CDPTHZ(L-1))
            RTDPX=AMAX1(0.0_r8,AMIN1(DLYR3(L),RTDPX) &
              -AMAX1(0.0,SDPTH(NZ)-CDPTHZ(L-1)-HTCTL(NZ)))
            FRTDPX(L,NZ)=RTDPX/DLYR3(L)
          ELSE
            FRTDPX(L,NZ)=0.0_r8
          ENDIF
        ENDIF

      ENDIF
      IF(WTRTG(L).GT.ZEROS)THEN
        FPQ(N,L,NZ)=AMAX1(0.0,WTRTD(N,L,NZ))/WTRTG(L)
      ELSE
        FPQ(N,L,NZ)=1.0
      ENDIF
      FPP(N,L,NZ)=FMN*FPQ(N,L,NZ)
      IF(RTDNP(N,L,NZ).GT.ZERO.AND.FRTDPX(L,NZ).GT.ZERO)THEN
        RRADL(N,L)=AMAX1(RRAD2X(N,NZ),SQRT((RTVLW(N,L,NZ) &
          /(1.0-PORT(N,NZ)))/(PICON*PP(NZ)*RTLGP(N,L,NZ))))
        PATH(N,L)=AMAX1(1.001*RRADL(N,L) &
          ,1.0/(SQRT(PICON*(RTDNP(N,L,NZ)/FRTDPX(L,NZ))/FMPR(L))))
        RTARR(N,L)=6.283*RTLGP(N,L,NZ)/FRTDPX(L,NZ)
      ELSE
        RRADL(N,L)=RRAD2M(N,NZ)
        PATH(N,L)=1.001*RRADL(N,L)
        RTARR(N,L)=6.283*RTLGP(N,L,NZ)
      ENDIF
    enddo
2000  CONTINUE
  end associate
  end subroutine UpdateRootProperty
!------------------------------------------------------------------------

  subroutine HandlingDivergence(I,J,NN,NZ,PSIST1,DIFF,FDMP)

  implicit none
  integer  , intent(in) :: NN, I, J
  integer  , intent(in) :: NZ
  REAL(R8) ,INTENT(IN) :: PSIST1(JZ1),DIFF
  real(r8), intent(inout):: FDMP
  real(r8) :: APSILT,APSIRT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM,FDMR
  real(r8) :: OSWT,WFNC

  integer :: N,L
! begin_execution
  associate(                         &
   RAZ     => plt_ew%RAZ       , &
   DTKC    => plt_ew%DTKC      , &
   OSMO    => plt_ew%OSMO      , &
   PSILO   => plt_ew%PSILO     , &
   TKC     => plt_ew%TKC       , &
   TKA     => plt_ew%TKA       , &
   TKS     => plt_ew%TKS       , &
   PSIRO   => plt_ew%PSIRO     , &
   TCC     => plt_ew%TCC       , &
   UPWTR   => plt_ew%UPWTR     , &
   PSIRG   => plt_ew%PSIRG     , &
   PSILG   => plt_ew%PSILG     , &
   PSIRT   => plt_ew%PSIRT     , &
   VHCPC   => plt_ew%VHCPC     , &
   PSILT   => plt_ew%PSILT     , &
   NU      => plt_site%NU      , &
   AREA3   => plt_site%AREA3   , &
   IYRC    => plt_site%IYRC    , &
   NI      => plt_morph%NI     , &
   NG      => plt_morph%NG     , &
   MY      => plt_morph%MY     , &
   CCPOLR  => plt_biom%CCPOLR  , &
   CZPOLR  => plt_biom%CZPOLR  , &
   CPPOLR  => plt_biom%CPPOLR  , &
   CCPOLP  => plt_biom%CCPOLP  , &
   CZPOLP  => plt_biom%CZPOLP  , &
   CPPOLP  => plt_biom%CPPOLP  , &
   WTSHT   => plt_biom%WTSHT   , &
   RA      => plt_photo%RA     , &
   RSMN    => plt_photo%RSMN   , &
   RC      => plt_photo%RC     , &
   RSMH    => plt_photo%RSMH   , &
   RCS     => plt_photo%RCS    , &
   FRADP   => plt_rad%FRADP    , &
   THRM1   => plt_rad%THRM1      &
  )
  IF(NN.GE.MXN)THEN
    WRITE(*,9999)IYRC,I,J,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)
    IF(DIFF.GT.0.5)THEN
      plt_rad%RAD1(NZ)=0.0_r8
      plt_ew%EFLXC(NZ)=0.0_r8
      plt_ew%SFLXC(NZ)=0.0_r8
      plt_ew%HFLXC(NZ)=0.0_r8
      plt_ew%EVAPC(NZ)=0.0_r8
      plt_ew%EP(NZ)=0.0_r8
      TKC(NZ)=TKA+DTKC(NZ)
      TCC(NZ)=TKC(NZ)-TC2K
      FTHRM=EMMC*2.04E-10*FRADP(NZ)*AREA3(NU)
      THRM1(NZ)=FTHRM*TKC(NZ)**4
      PSILT(NZ)=PSIST1(NG(NZ))
      APSILT=ABS(PSILT(NZ))
      FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLP(NZ)+CZPOLP(NZ)+CPPOLP(NZ)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILO(NZ)=FDMP/0.16*OSMO(NZ)-RGAS*TKC(NZ)*FDMP*CCPOLT/OSWT
      PSILG(NZ)=AMAX1(0.0,PSILT(NZ)-PSILO(NZ))
      WFNC=EXP(RCS(NZ)*PSILG(NZ))
      RC(NZ)=RSMN(NZ)+(RSMH(NZ)-RSMN(NZ))*WFNC
      RA(NZ)=RAZ(NZ)
      VHCPC(NZ)=cpw*(WTSHT(NZ)*10.0E-06)
      DTKC(NZ)=0.0_r8
      DO 4290 N=1,MY(NZ)
        DO  L=NU,NI(NZ)
          PSIRT(N,L,NZ)=PSIST1(L)
          APSIRT=ABS(PSIRT(N,L,NZ))
          FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
          CCPOLT=CCPOLR(N,L,NZ)+CZPOLR(N,L,NZ)+CPPOLR(N,L,NZ)
          OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
          PSIRO(N,L,NZ)=FDMR/0.16*OSMO(NZ) &
            -RGAS*TKS(L)*FDMR*CCPOLT/OSWT
          PSIRG(N,L,NZ)=AMAX1(0.0,PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
          UPWTR(N,L,NZ)=0.0_r8
      enddo
4290  CONTINUE
    ENDIF
  ENDIF
  end associate
  end subroutine HandlingDivergence
!------------------------------------------------------------------------------
  function CanopyEnergyH2OIteration(I,J,NZ,FPC,WVPLT,&
    PSIST1,PARHC,DIFF,UPRT,VFLXC,FDMP,RSRS,FPQ,VOLPU,VOLWU,TKCX,CNDT,&
    VHCPX,HFLWC1,PSILH,ILYR) result(NN)
  implicit none
  integer  , intent(in) :: I, J
  integer  , intent(in) :: NZ
  real(r8) , intent(in) :: FPC,WVPLT,PSIST1(JZ1)
  real(r8) , intent(in) :: RSRS(2,JZ1),FPQ(2,JZ1,JP1),VOLPU(JZ1),VOLWU(JZ1)
  real(r8) , intent(in) :: TKCX,CNDT,VHCPX,HFLWC1,PSILH
  integer  , intent(in) :: ILYR(2,JZ1)
  real(r8) , intent(out):: PARHC,DIFF,UPRT,VFLXC,FDMP
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: DTHS1
  real(r8) :: DIFFZ,DIFFU,DPSI
  real(r8) :: EX,EPX
  real(r8) :: FTHRM,FDMR
  real(r8) :: FDTHS
  real(r8) :: HFLXS
  real(r8) :: OSWT
  real(r8) :: PAREX,PARHX,PSIL2,PAREC
  real(r8) :: PSILC
  real(r8) :: RSSUX,RSSU,RA1,RSSZ
  real(r8) :: TKC1,TKCY
  real(r8) :: UPRTX
  real(r8) :: VOLWPX,VPC,VOLWPZ
  real(r8) :: XC,WFNC
  real(r8) :: RI
  integer  :: IC,ICHK
!     return variables
  integer :: NN
!     local variables
  integer :: N,L
!     begin_execution

  associate(                          &
    NU      => plt_site%NU      , &
    AREA3   => plt_site%AREA3   , &
    PSILO   => plt_ew%PSILO     , &
    OSMO    => plt_ew%OSMO      , &
    RAZ     => plt_ew%RAZ       , &
    VPA     => plt_ew%VPA       , &
    TKA     => plt_ew%TKA       , &
    RIB     => plt_ew%RIB       , &
    EP      => plt_ew%EP        , &
    PSILT   => plt_ew%PSILT     , &
    VOLWP   => plt_ew%VOLWP     , &
    UPWTR   => plt_ew%UPWTR     , &
    TKCZ    => plt_ew%TKCZ      , &
    EVAPC   => plt_ew%EVAPC     , &
    VHCPC   => plt_ew%VHCPC     , &
    PSILG   => plt_ew%PSILG     , &
    VOLWC   => plt_ew%VOLWC     , &
    FLWC    => plt_ew%FLWC      , &
    EFLXC   => plt_ew%EFLXC     , &
    ZEROL   => plt_biom%ZEROL   , &
    ZEROP   => plt_biom%ZEROP   , &
    CCPOLP  => plt_biom%CCPOLP  , &
    CZPOLP  => plt_biom%CZPOLP  , &
    CPPOLP  => plt_biom%CPPOLP  , &
    NI      => plt_morph%NI     , &
    MY      => plt_morph%MY     , &
    RSMN    => plt_photo%RSMN   , &
    RCS     => plt_photo%RCS    , &
    RA      => plt_photo%RA     , &
    RSMH    => plt_photo%RSMH   , &
    RC      => plt_photo%RC     , &
    RAD1    => plt_rad%RAD1     , &
    RADC    => plt_rad%RADC     , &
    THS     => plt_rad%THS      , &
    FRADP   => plt_rad%FRADP    , &
    THRM1   => plt_rad%THRM1    , &
    THRMGX  => plt_rad%THRMGX     &
  )
  CCPOLT=CCPOLP(NZ)+CZPOLP(NZ)+CPPOLP(NZ)
  OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
  FTHRM=EMMC*2.04E-10*FRADP(NZ)*AREA3(NU)
  FDTHS=(THS+THRMGX)*FRADP(NZ)
!     RAZ=canopy isothermal boundary later resistance

  UPRT=0.0_r8
  PAREX=FPC*AREA3(NU)
  PARHX=FPC*AREA3(NU)*1.25E-03
  RA1=RAZ(NZ)

  IC=0
  XC=0.5
  ICHK=0
  PSIL2=0.0_r8
  EPX=0.0_r8
  UPRTX=0.0_r8
  VOLWPX=0.0_r8

  DO NN=1,MXN
!
!     NET RADIATION FROM ABSORBED SW AND NET LW
!
!     THRM1=LW emitted by canopy
!     DTHS1=net LW absorbed by canopy
!     RADC=total SW absorbed by canopy
!     RAD1=net SW+LW absorbed by canopy
!
    TKC1=TKCZ(NZ)
    THRM1(NZ)=FTHRM*TKC1**4
    DTHS1=FDTHS-THRM1(NZ)*2.0
    RAD1(NZ)=RADC(NZ)+DTHS1
!
!     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
!
!     RI=Ricardson's number
!     RA=canopy boundary layer resistance
!     PAREC,PARHC=canopy latent,sensible heat conductance
!
    RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB*(TKA-TKC1)))

    RA(NZ)=AMAX1(RACM,0.9_r8*RA1,AMIN1(1.1_r8*RA1,RAZ(NZ)/(1.0_r8-10.0_r8*RI)))
    RA1=RA(NZ)
    PAREC=PAREX/RA(NZ)
    PARHC=PARHX/RA(NZ)
!
!     CANOPY WATER AND OSMOTIC POTENTIALS
!
!     PSILT=canopy total water potential
!     FDMP=dry matter content
!     OSMO=osmotic potential at PSILT=0 from PFT file
!     PSILO,PSILG=canopy osmotic,turgor water potential
!
    APSILT=ABS(PSILT(NZ))
    FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
    PSILO(NZ)=FDMP/0.16*OSMO(NZ)-RGAS*TKC1*FDMP*CCPOLT/OSWT
    PSILG(NZ)=AMAX1(0.0,PSILT(NZ)-PSILO(NZ))
!
!     CANOPY STOMATAL RESISTANCE
!
!     RCS=shape parameter for RC vs PSILG from PFT file
!     RC=canopy stomatal resistance
!     RSMN=minimum RC at PSILT=0 from stomate.f
!     RSMX=cuticular resistance from PFT file
!
    WFNC=EXP(RCS(NZ)*PSILG(NZ))
    RC(NZ)=RSMN(NZ)+(RSMH(NZ)-RSMN(NZ))*WFNC
!
!     CANOPY VAPOR PRESSURE AND EVAPORATION OF INTERCEPTED WATER
!     OR TRANSPIRATION OF UPTAKEN WATER
!
!     VPC,VPA=vapor pressure inside canopy, in atmosphere
!     TKC1=canopy temperature
!     PSILT=canopy total water potential
!     EX=canopy-atmosphere water flux
!     RA,RZ=canopy boundary layer,surface resistance
!     EVAPC=water flux to,from canopy surfaces
!     VOLWC=water volume on canopy surfaces
!     EP=water flux to,from inside canopy
!     EFLXC=canopy latent heat flux
!     VFLXC=convective heat flux from EFLXC
!     VAP=latent heat of evaporation
!
      !VPC=2.173E-03/TKC1 &
      !*0.61*EXP(5360.0*(3.661E-03-1.0/TKC1)) &
    vpc=vapsat(tkc1)*EXP(18.0*PSILT(NZ)/(RGAS*TKC1))
    EX=PAREC*(VPA-VPC)
    IF(EX.GT.0.0)THEN
      EVAPC(NZ)=EX*RA(NZ)/(RA(NZ)+RZ)
      EX=0.0_r8
    ELSEIF(EX.LE.0.0.AND.VOLWC(NZ).GT.0.0)THEN
      EVAPC(NZ)=AMAX1(EX*RA(NZ)/(RA(NZ)+RZ),-VOLWC(NZ))
      EX=EX-EVAPC(NZ)
    ENDIF
    EP(NZ)=EX*RA(NZ)/(RA(NZ)+RC(NZ))
    EFLXC(NZ)=(EP(NZ)+EVAPC(NZ))*VAP
    VFLXC=EVAPC(NZ)*cpw*TKC1
!
!     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
!
!     HFLXS=initial estimate of sensible+storage heat flux
!     HFLWC1=convective heat flux from precip to canopy
!
    HFLXS=RAD1(NZ)+EFLXC(NZ)+VFLXC+HFLWC1
!
!     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
!
!     VHCPC=canopy heat capacity
!     TKCY=equilibrium canopy temperature for HFLXS
!
    VHCPC(NZ)=VHCPX+cpw*(EVAPC(NZ)+FLWC(NZ))
    TKCY=(TKCX*VHCPX+TKA*PARHC+HFLXS)/(VHCPC(NZ)+PARHC)
    TKCY=AMIN1(TKA+10.0,AMAX1(TKA-10.0,TKCY))
!
!     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
!
!     XC,IC=magnitude,direction of change in canopy temp for next cycle
!
    IF((IC.EQ.0.AND.TKCY.GT.TKC1).OR.(IC.EQ.1.AND.TKCY.LT.TKC1))THEN
      XC=0.5*XC
    ENDIF
    TKCZ(NZ)=TKC1+0.1*(TKCY-TKC1)
    IF(TKCY.GT.TKC1)THEN
      IC=1
    ELSE
      IC=0
    ENDIF

!
!     IF CONVERGENCE CRITERION IS MET OR ON EVERY TENTH ITERATION,
!     PROCEED TO WATER BALANCE
!
!     PSILC=canopy water potential adjusted for canopy height
!
    IF(ABS(TKCY-TKC1).LT.0.05_r8.OR.(NN/10)*10.EQ.NN)THEN
      UPRT=0.0_r8
      PSILC=PSILT(NZ)-PSILH
!
!     ROOT WATER UPTAKE FROM SOIL-CANOPY WATER POTENTIALS,
!     SOIL + ROOT HYDRAULIC RESISTANCES
!
!     ILYR=rooted layer flag
!     UPWTR=root water uptake from soil layer
!     VOLWU,VOLPU=water volume available for uptake,air volume
!     FPQ=PFT fraction of biome root mass
!     PSILC=canopy water potential adjusted for canopy height
!     PSIST1=total soil water potential PSIST adjusted for surf elevn
!     RSRS=total soil+root resistance
!     UPRT=total water uptake from soil profile
!
      DO N=1,MY(NZ)
        DO  L=NU,NI(NZ)
          IF(ILYR(N,L).EQ.1)THEN
            UPWTR(N,L,NZ)=AMAX1(AMIN1(0.0,-VOLWU(L)*FPQ(N,L,NZ)) &
              ,AMIN1((PSILC-PSIST1(L))/RSRS(N,L),VOLPU(L)*FPQ(N,L,NZ)))
            IF(UPWTR(N,L,NZ).GT.0.0)THEN
              UPWTR(N,L,NZ)=0.1*UPWTR(N,L,NZ)
            ENDIF
            UPRT=UPRT+UPWTR(N,L,NZ)
          ELSE
            UPWTR(N,L,NZ)=0.0_r8
          ENDIF
        enddo
      ENDDO
!
!     TEST TRANSPIRATION - ROOT WATER UPTAKE VS. CHANGE IN CANOPY
!     WATER STORAGE
!
!     VOLWPZ,VOLWP=canopy water content
!     DIFFZ,DIFFU=change in canopy water content,transpiration-uptake
!     DIFF=normalized difference between DIFFZ and DIFFU
!     5.0E-03=acceptance criterion for DIFF
!     RSSZ=change in canopy water potl vs change in canopy water cnt
!     RSSU=change in canopy water potl vs change in transpiration

      VOLWPZ=1.0E-06_r8*WVPLT/FDMP
      DIFFZ=VOLWPZ-VOLWP(NZ)
      DIFFU=EP(NZ)-UPRT
      IF(test_aneb(UPRT,0.0_r8))THEN
        DIFF=ABS((DIFFU-DIFFZ)/UPRT)
      ELSE
        DIFF=ABS((DIFFU-DIFFZ)/VOLWPZ)
      ENDIF
      IF(DIFF.LT.5.0E-03_r8)THEN
        IF(ICHK.EQ.1)EXIT
        ICHK=1
        CALL STOMATEs(I,J,NZ)
        CYCLE
      ENDIF
      IF(ABS(VOLWPZ-VOLWPX).GT.ZEROP(NZ))THEN
        RSSZ=ABS((PSILT(NZ)-PSIL2)/(VOLWPZ-VOLWPX))
      ELSEIF(CNDT.GT.ZEROP(NZ))THEN
        RSSZ=1.0_r8/CNDT
      ELSE
        RSSZ=ZEROL(NZ)
      ENDIF
      IF(ABS(EP(NZ)-EPX).GT.ZEROP(NZ))THEN
        RSSUX=ABS((PSILT(NZ)-PSIL2)/(EP(NZ)-EPX))
        IF(CNDT.GT.ZEROP(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(ABS(UPRT-UPRTX).GT.ZEROP(NZ))THEN
        RSSUX=ABS((PSILT(NZ)-PSIL2)/(UPRT-UPRTX))
        IF(CNDT.GT.ZEROP(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(CNDT.GT.ZEROP(NZ))THEN
        RSSU=1.0_r8/CNDT
      ELSE
        RSSU=ZEROL(NZ)
      ENDIF
!
!     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
!     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
!     WATER STORAGE
!
!     DPSI=change in PSILT for next convergence cycle
!     1.0E-03=acceptance criterion for DPSI
!
      DPSI=AMIN1(AMIN1(RSSZ,RSSU)*(DIFFU-DIFFZ),ABS(PSILT(NZ)))

!     IF CONVERGENCE CRITERION IS MET THEN FINISH,
!     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
!     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
!
      IF(.not.((NN.GE.30.AND.ABS(DPSI).LT.1.0E-03_r8).OR.NN.GE.MXN))then
        PSIL2=PSILT(NZ)
        EPX=EP(NZ)
        UPRTX=UPRT
        VOLWPX=VOLWPZ
        PSILT(NZ)=AMIN1(0.0_r8,PSILT(NZ)+0.5_r8*DPSI)
        XC=0.50_r8
        cycle
!
!     RESET MIN STOMATAL RESISTANCE IN STOMATE.F BEFORE FINAL ITERATION
!
      ELSE
        IF(ICHK.EQ.1)EXIT
        ICHK=1
        CALL STOMATEs(I,J,NZ)
      ENDIF
    ENDIF
  ENDDO
!4000  CONTINUE
!4500  CONTINUE
  end associate
  end function CanopyEnergyH2OIteration
!------------------------------------------------------------------------
  subroutine CalcResistance(NZ,PATH,RRADL,RTARR,&
    RSRT,RSRG,RSR1,RSR2,RSSX,RSRS,CNDT,PSILH,ILYR)

  implicit none
  integer, intent(in)   :: NZ
  real(r8), intent(in)  :: PATH(2,JZ1),RRADL(2,JZ1),RTARR(2,JZ1)
  real(r8), intent(out) :: RSRT(2,JZ1),RSRG(2,JZ1)
  real(r8), intent(out) :: RSR1(2,JZ1),RSR2(2,JZ1)
  real(r8), intent(out) :: RSSX(2,JZ1),RSRS(2,JZ1),CNDT
  real(r8), intent(out) :: PSILH
  integer, intent(out) :: ILYR(2,JZ1)
  real(r8) :: FRADW,FRAD1,FRAD2
  real(r8) :: RSSL,RTAR2
  integer :: N, L
  associate(                          &
    DPTHZ  => plt_site%DPTHZ    , &
    PP     => plt_site%PP       , &
    VOLWM  => plt_site%VOLWM    , &
    ZERO   => plt_site%ZERO     , &
    ZEROS2 => plt_site%ZEROS2   , &
    NU     => plt_site%NU       , &
    PSILT  => plt_ew%PSILT      , &
    ZEROP  => plt_biom%ZEROP    , &
    THETW  => plt_soilchem%THETW, &
    VOLA   => plt_soilchem%VOLA , &
    CNDU   => plt_soilchem%CNDU , &
    VOLX   => plt_soilchem%VOLX , &
    RTNL   => plt_morph%RTNL    , &
    RTN1   => plt_morph%RTN1    , &
    RRAD2M => plt_morph%RRAD2M  , &
    RSRA   => plt_morph%RSRA    , &
    RRAD2  => plt_morph%RRAD2   , &
    RTDNP  => plt_morph%RTDNP   , &
    HTSTZ  => plt_morph%HTSTZ   , &
    RTLGP  => plt_morph%RTLGP   , &
    RTLGA  => plt_morph%RTLGA   , &
    RRAD1  => plt_morph%RRAD1   , &
    MY     => plt_morph%MY      , &
    RSRR   => plt_morph%RSRR    , &
    ZC     => plt_morph%ZC      , &
    NG     => plt_morph%NG      , &
    NI     => plt_morph%NI        &
  )

  !     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
  !
  !     HTSTZ=canopy height for water uptake
  !     PSILH=gravimetric water potential at HTSTZ
  !     FRADW=conducting elements of stalk relative to those of primary root
  !     PSILT=canopy total water potential
  !     EMODW=wood modulus of elasticity (MPa)
!
  CNDT=0.0_r8
  HTSTZ(NZ)=0.80*ZC(NZ)
  PSILH=-0.0098*HTSTZ(NZ)
  FRADW=1.0E+04*(AMAX1(0.5,1.0+PSILT(NZ)/EMODW))**4
!
  !     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
  !
  !      VOLX,VOLWM,THETW=soil,water volume,content
  !     RTDNP,RTLGP=root length density,root length per plant
  !     CNDU=soil hydraulic conductivity for root uptake
  !     RTN1,RTNL=number of root,myco primary,secondary axes
  !     ILYR:1=rooted,0=not rooted
  !     N:1=root,2=mycorrhizae
!
  DO 3880 N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(VOLX(L).GT.ZEROS2 &
        .AND.VOLWM(NPH,L).GT.ZEROS2 &
        .AND.RTDNP(N,L,NZ).GT.ZERO &
        .AND.CNDU(L).GT.ZERO &
        .AND.RTN1(1,L,NZ).GT.ZEROP(NZ) &
        .AND.RTNL(N,L,NZ).GT.ZEROP(NZ) &
        .AND.THETW(L).GT.ZERO)THEN
        ILYR(N,L)=1
        !
        !     SOIL HYDRAULIC RESISTANCE FROM RADIAL UPTAKE GEOMETRY
        !     AND SOIL HYDRAULIC CONDUCTIVITY
        !
        !     RSSX=soil hydraulic resistance
        !     PP=plant population
        !     PATH=path length of water and nutrient uptake
        !     RRADL,RTARR=root radius,surface/radius area
        !
        RSSL=(LOG(PATH(N,L)/RRADL(N,L))/RTARR(N,L))/PP(NZ)
        RSSX(N,L)=RSSL/CNDU(L)
        !
        !     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
        !     ENTERED IN 'READQ'
        !
        !     RRAD2=secondary root radius
        !     RTLGP=root length per plant
        !     RSRG=radial resistance
        !     RSRR=radial resistivity from PFT file
        !     VOLA,VOLWM=soil micropore,water volume
        !
        RTAR2=6.283*RRAD2(N,L,NZ)*RTLGP(N,L,NZ)*PP(NZ)
        RSRG(N,L)=RSRR(N,NZ)/RTAR2*VOLA(L)/VOLWM(NPH,L)
!
        !     ROOT AXIAL RESISTANCE FROM RADII AND LENGTHS OF PRIMARY AND
        !     SECONDARY ROOTS AND FROM AXIAL RESISTIVITY ENTERED IN 'READQ'
        !
        !     FRAD1,FRAD2=primary,secondary root radius relative to maximum
        !     secondary radius from PFT file RRAD2M at which RSRA is defined
        !     RRAD1,RRAD2=primary,secondary root radius
        !     RSRA=axial resistivity from PFT file
        !     DPTHZ=depth of primary root from surface
        !     RSR1,RSR2=axial resistance of primary,secondary roots
        !     RTLGA=average secondary root length
        !     RTN1,RTNL=number of primary,secondary axes
!
        FRAD1=(RRAD1(N,L,NZ)/RRAD2M(N,NZ))**4
        RSR1(N,L)=RSRA(N,NZ)*DPTHZ(L)/(FRAD1*RTN1(1,L,NZ)) &
          +RSRA(1,NZ)*HTSTZ(NZ)/(FRADW*RTN1(1,L,NZ))
        FRAD2=(RRAD2(N,L,NZ)/RRAD2M(N,NZ))**4
        RSR2(N,L)=RSRA(N,NZ)*RTLGA(N,L,NZ)/(FRAD2*RTNL(N,L,NZ))
      ELSE
        ILYR(N,L)=0
      ENDIF
    enddo
3880  CONTINUE

  DO 3890 N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(ILYR(N,L).EQ.1)THEN
        !
        !     TOTAL ROOT RESISTANCE = SOIL + RADIAL + AXIAL
        !
        !     RSRT=root radial+axial resistance
        !     RSRS=total soil+root resistance
        !     CNDT=total soil+root conductance for all layers
        !
        RSRT(N,L)=RSRG(N,L)+RSR1(N,L)+RSR2(N,L)
        RSRS(N,L)=RSSX(N,L)+RSRT(N,L)
        CNDT=CNDT+1.0/RSRS(N,L)

      ENDIF
    enddo
3890  CONTINUE
  end associate
  end subroutine CalcResistance
!------------------------------------------------------------------------

  subroutine HandleBareSoil(NZ,PSIST1,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: PSIST1(JZ1)
  real(r8), intent(out):: FDMP
  integer :: N,L
  real(r8) :: APSILT,APSIRT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM,FDMR
  real(r8) :: OSWT,WFNC

! begin_execution
  associate(                         &
    TCC    =>  plt_ew%TCC       , &
    OSMO   =>  plt_ew%OSMO      , &
    TKW    =>  plt_ew%TKW       , &
    TKC    =>  plt_ew%TKC       , &
    TKS    =>  plt_ew%TKS       , &
    TKA    =>  plt_ew%TKA       , &
    DPTHS  =>  plt_ew%DPTHS     , &
    DTKC   =>  plt_ew%DTKC      , &
    RAZ    =>  plt_ew%RAZ       , &
    VHCPC  =>  plt_ew%VHCPC     , &
    PSILO  =>  plt_ew%PSILO     , &
    PSIRG  =>  plt_ew%PSIRG     , &
    UPWTR  =>  plt_ew%UPWTR     , &
    PSILG  =>  plt_ew%PSILG     , &
    PSILT  =>  plt_ew%PSILT     , &
    PSIRO  =>  plt_ew%PSIRO     , &
    PSIRT  =>  plt_ew%PSIRT     , &
    EP     =>  plt_ew%EP        , &
    HFLXC  =>  plt_ew%HFLXC     , &
    EFLXC  =>  plt_ew%EFLXC     , &
    SFLXC  =>  plt_ew%SFLXC     , &
    EVAPC  =>  plt_ew%EVAPC     , &
    NU     =>  plt_site%NU      , &
    ZERO   =>  plt_site%ZERO    , &
    AREA3  =>  plt_site%AREA3   , &
    CCPOLP =>  plt_biom%CCPOLP  , &
    CZPOLP =>  plt_biom%CZPOLP  , &
    CPPOLP =>  plt_biom%CPPOLP  , &
    CCPOLR =>  plt_biom%CCPOLR  , &
    CZPOLR =>  plt_biom%CZPOLR  , &
    CPPOLR =>  plt_biom%CPPOLR  , &
    WTSHT  =>  plt_biom%WTSHT   , &
    NI     =>  plt_morph%NI     , &
    ZC     =>  plt_morph%ZC     , &
    NG     =>  plt_morph%NG     , &
    MY     =>  plt_morph%MY     , &
    RCS    =>  plt_photo%RCS    , &
    RA     =>  plt_photo%RA     , &
    RC     =>  plt_photo%RC     , &
    RSMN   =>  plt_photo%RSMN   , &
    RSMH   =>  plt_photo%RSMH   , &
    FRADP  =>  plt_rad%FRADP    , &
    THRM1  =>  plt_rad%THRM1    , &
    RAD1   =>  plt_rad%RAD1       &
  )
  RAD1(NZ)=0.0_r8
  EFLXC(NZ)=0.0_r8
  SFLXC(NZ)=0.0_r8
  HFLXC(NZ)=0.0_r8
  EVAPC(NZ)=0.0_r8
  EP(NZ)=0.0_r8
  IF(ZC(NZ).GE.DPTHS-ZERO)THEN
    TKC(NZ)=TKA
  ELSE
    TKC(NZ)=TKW
  ENDIF
  TCC(NZ)=TKC(NZ)-TC2K
  FTHRM=EMMC*2.04E-10_r8*FRADP(NZ)*AREA3(NU)
  THRM1(NZ)=FTHRM*TKC(NZ)**4
  PSILT(NZ)=PSIST1(NG(NZ))
  APSILT=ABS(PSILT(NZ))
  FDMP=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
  CCPOLT=CCPOLP(NZ)+CZPOLP(NZ)+CPPOLP(NZ)
  OSWT=36.0_r8+840.0_r8*AMAX1(0.0_r8,CCPOLT)
  PSILO(NZ)=FDMP/0.16_r8*OSMO(NZ)-RGAS*TKC(NZ)*FDMP*CCPOLT/OSWT
  PSILG(NZ)=AMAX1(0.0_r8,PSILT(NZ)-PSILO(NZ))
  WFNC=EXP(RCS(NZ)*PSILG(NZ))
  RC(NZ)=RSMN(NZ)+(RSMH(NZ)-RSMN(NZ))*WFNC
  RA(NZ)=RAZ(NZ)
  VHCPC(NZ)=cpw*(WTSHT(NZ)*10.0E-06_r8)
  DTKC(NZ)=0.0_r8
  DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      PSIRT(N,L,NZ)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ))
      FDMR=0.16_r8+0.10_r8*APSIRT/(0.05*APSIRT+2.0_r8)
      CCPOLT=CCPOLR(N,L,NZ)+CZPOLR(N,L,NZ)+CPPOLR(N,L,NZ)
      OSWT=36.0_r8+840.0_r8*AMAX1(0.0_r8,CCPOLT)
      PSIRO(N,L,NZ)=FDMR/0.16_r8*OSMO(NZ)-RGAS*TKS(L)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ)=AMAX1(0.0_r8,PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      UPWTR(N,L,NZ)=0.0_r8
    enddo
  ENDDO
  end associate
  end subroutine HandleBareSoil
!------------------------------------------------------------------------

  subroutine UpdateCanopyWater(NZ,PARHC,PSIST1,RSRT,RSSX,&
    RSRS,TKCX,VHCPX,HFLWC1,UPRT,VFLXC,ILYR)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: PARHC,PSIST1(JZ1),RSRT(2,JZ1)
  real(r8), intent(in) :: RSSX(2,JZ1),RSRS(2,JZ1)
  real(r8), intent(in) :: TKCX,VHCPX,HFLWC1,UPRT,VFLXC
  integer , intent(in) :: ILYR(2,JZ1)
  real(r8) :: APSIRT
  real(r8) :: CCPOLT
  real(r8) :: FDMR
  real(r8) :: OSWT
  integer :: N,L
  associate(                           &
    NU       => plt_site%NU      , &
    OSMO     => plt_ew%OSMO      , &
    TKA      => plt_ew%TKA       , &
    TKC      => plt_ew%TKC       , &
    PSILZ    => plt_ew%PSILZ     , &
    TKS      => plt_ew%TKS       , &
    TKCZ     => plt_ew%TKCZ      , &
    VOLWC    => plt_ew%VOLWC     , &
    VOLWP    => plt_ew%VOLWP     , &
    VHCPC    => plt_ew%VHCPC     , &
    HFLXC    => plt_ew%HFLXC     , &
    FLWC     => plt_ew%FLWC      , &
    EP       => plt_ew%EP        , &
    EVAPC    => plt_ew%EVAPC     , &
    PSIRT    => plt_ew%PSIRT     , &
    SFLXC    => plt_ew%SFLXC     , &
    PSIRO    => plt_ew%PSIRO     , &
    PSIRG    => plt_ew%PSIRG     , &
    PSILT    => plt_ew%PSILT     , &
    CCPOLR   => plt_biom%CCPOLR  , &
    CZPOLR   => plt_biom%CZPOLR  , &
    CPPOLR   => plt_biom%CPPOLR  , &
    MY       => plt_morph%MY     , &
    NI       => plt_morph%NI       &
  )
  !
  !     CANOPY SURFACE WATER STORAGE, SENSIBLE AND STORAGE HEAT FLUXES
  !     (NOT EXPLICITLY CALCULATED IN CONVERGENCE SOLUTION)
  !
  !     VOLWP,VOLWC=water volume in canopy,on canopy surfaces
  !     SFLXC,HFLXC=canopy sensible,storage heat fluxes
  !     VHCPX,VHCPC=previous,current canopy heat capacity
  !     PARHC=canopy sensible heat conductance
  !     VFLXC=convective heat flux from latent heat flux
  !     HFLWC1=convective heat flux from precip to canopy
!
  VOLWP(NZ)=VOLWP(NZ)+EP(NZ)-UPRT
  VOLWC(NZ)=VOLWC(NZ)+FLWC(NZ)+EVAPC(NZ)
  SFLXC(NZ)=PARHC*(TKA-TKCZ(NZ))
  HFLXC(NZ)=TKCX*VHCPX-TKCZ(NZ)*VHCPC(NZ)+VFLXC+HFLWC1
  !
  !     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
  !
  !     PSIRT,PSILT=root,canopy total water potential
  !     PSIST1=total soil water potential PSIST adjusted for surf elevn
  !     RSSX,RSRS,RSRT=soil,soil+root,root radial+axial resistance
  !     PSIRO,PSIRG=root osmotic,turgor water potential
  !     FDMR=dry matter content
  !     OSMO=osmotic potential at PSIRT=0 from PFT file
!
  !
  DO 4505 N=1,MY(NZ)
    DO 4510 L=NU,NI(NZ)
      IF(ILYR(N,L).EQ.1)THEN
        PSIRT(N,L,NZ)=AMIN1(0.0,(PSIST1(L)*RSRT(N,L) &
          +PSILT(NZ)*RSSX(N,L))/RSRS(N,L))
        APSIRT=ABS(PSIRT(N,L,NZ))
        FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
        CCPOLT=CCPOLR(N,L,NZ)+CZPOLR(N,L,NZ)+CPPOLR(N,L,NZ)
        OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
        PSIRO(N,L,NZ)=FDMR/0.16*OSMO(NZ) &
          -RGAS*TKS(L)*FDMR*CCPOLT/OSWT
        PSIRG(N,L,NZ)=AMAX1(0.0,PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      ELSE
        PSIRT(N,L,NZ)=PSIST1(L)
        APSIRT=ABS(PSIRT(N,L,NZ))
        FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
        CCPOLT=CCPOLR(N,L,NZ)+CZPOLR(N,L,NZ)+CPPOLR(N,L,NZ)
        OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
        PSIRO(N,L,NZ)=FDMR/0.16*OSMO(NZ) &
          -RGAS*TKS(L)*FDMR*CCPOLT/OSWT
        PSIRG(N,L,NZ)=AMAX1(0.0,PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      ENDIF
4510  CONTINUE
4505  CONTINUE
  end associate
  end subroutine UpdateCanopyWater
!------------------------------------------------------------------------

  subroutine SetCanopyGrowthFuncs(NZ)

  implicit none
  integer, intent(in) :: NZ
  real(r8) :: ACTV,RTK,STK,TKGO,TKSO
  integer :: L
  associate(                          &
    TCC    =>  plt_ew%TCC       , &
    TKC    =>  plt_ew%TKC       , &
    TKS    =>  plt_ew%TKS       , &
    PSILZ  =>  plt_ew%PSILZ     , &
    PSILT  =>  plt_ew%PSILT     , &
    NU     =>  plt_site%NU      , &
    CHILL  =>  plt_photo%CHILL  , &
    OFFST  =>  plt_pheno%OFFST  , &
    CTC    =>  plt_pheno%CTC    , &
    TFN4   =>  plt_pheno%TFN4   , &
    TCG    =>  plt_pheno%TCG    , &
    TKG    =>  plt_pheno%TKG    , &
    IDAY   =>  plt_pheno%IDAY   , &
    TFN3   =>  plt_pheno%TFN3   , &
    NI     =>  plt_morph%NI     , &
    NB1    =>  plt_morph%NB1      &
  )
  !
  !     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
  !     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
  !
  IF(IDAY(1,NB1(NZ),NZ).EQ.0)THEN
    TKG(NZ)=TKS(NU)
    !     ELSEIF((IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)
    !    2.AND.IDAY(2,NB1(NZ),NZ).EQ.0)THEN
    !     TKG(NZ)=TKS(NU)
  ELSE
    TKG(NZ)=TKC(NZ)
  ENDIF
  TCG(NZ)=TKG(NZ)-TC2K
  !
  !     ARRHENIUS FUNCTION FOR CANOPY AND ROOT GROWTH WITH OFFSET
  !     FOR ZONE OF THERMAL ADAPTATION ENTERED IN 'READQ'
  !
  !     TKG,TKGO=canopy temperature,canopy temp used in Arrhenius eqn
  !     TKS,TKSO=soil temperature,soil temp used in Arrhenius eqn
  !     OFFST=shift in Arrhenius curve for thermal adaptation
  !     TFN3,TFN4=temperature function for canopy,root growth (25 oC =1)
  !     RGAS,710.0=gas constant,enthalpy
  !     62500,197500,222500=energy of activn,high,low temp inactivn(KJ mol-1)
  !     PSILZ=minimum daily canopy water potential
  !
  TKGO=TKG(NZ)+OFFST(NZ)
  RTK=RGAS*TKGO
  STK=710.0*TKGO
  ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
  TFN3(NZ)=EXP(25.229-62500/RTK)/ACTV
  DO 100 L=NU,NI(NZ)
    TKSO=TKS(L)+OFFST(NZ)
    RTK=RGAS*TKSO
    STK=710.0*TKSO
    ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
    TFN4(L,NZ)=EXP(25.229-62500/RTK)/ACTV
100   CONTINUE
  PSILZ(NZ)=AMIN1(PSILZ(NZ),PSILT(NZ))
  !
  !     DIURNAL CHILLING
  !
  !     CTC=chilling temperature from PFT file
  !     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
  !
  IF(TCC(NZ).LT.CTC(NZ))THEN
    CHILL(NZ)=AMIN1(24.0,CHILL(NZ)+1.0)
  ELSE
    CHILL(NZ)=AMAX1(0.0,CHILL(NZ)-1.0)
  ENDIF
  end associate
  end subroutine SetCanopyGrowthFuncs

end module UptakesMod
