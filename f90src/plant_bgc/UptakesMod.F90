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
    RAs1     => plt_photo%RAs1     , &
    VOLWPs1  => plt_ew%VOLWPs1     , &
    TCCs1    => plt_ew%TCCs1       , &
    TKCZs1   => plt_ew%TKCZs1      , &
    FLWCs1   => plt_ew%FLWCs1      , &
    EFLXCs1  => plt_ew%EFLXCs1     , &
    TKCs1    => plt_ew%TKCs1       , &
    PSILTs1  => plt_ew%PSILTs1     , &
    SFLXCs1  => plt_ew%SFLXCs1     , &
    TSHCs1   => plt_ew%TSHCs1      , &
    DTKCs1   => plt_ew%DTKCs1      , &
    EVAPCs1  => plt_ew%EVAPCs1     , &
    TKAs1    => plt_ew%TKAs1       , &
    TLECs1   => plt_ew%TLECs1      , &
    EPs1     => plt_ew%EPs1        , &
    VOLWCs1  => plt_ew%VOLWCs1     , &
    CDPTHZs1 => plt_site%CDPTHZs1  , &
    PPTs1    => plt_site%PPTs1     , &
    NPs1     => plt_site%NPs1      , &
    PPs1     => plt_site%PPs1      , &
    NUs1     => plt_site%NUs1      , &
    AREA3s1  => plt_site%AREA3s1   , &
    ZEROSs1  => plt_site%ZEROSs1   , &
    WTRTDs1  => plt_biom%WTRTDs1   , &
    ZEROLs1  => plt_biom%ZEROLs1   , &
    ZEROPs1  => plt_biom%ZEROPs1   , &
    WTLSs1   => plt_biom%WTLSs1    , &
    WVSTKs1  => plt_biom%WVSTKs1   , &
    IDAYs1   => plt_pheno%IDAYs1   , &
    OSTRs1   => plt_pheno%OSTRs1   , &
    IFLGCs1  => plt_pheno%IFLGCs1  , &
    ARLFCs1  => plt_morph%ARLFCs1  , &
    RTDP1s1  => plt_morph%RTDP1s1  , &
    ARLSSs1  => plt_morph%ARLSSs1  , &
    ARLFSs1  => plt_morph%ARLFSs1  , &
    SDPTHs1  => plt_morph%SDPTHs1  , &
    NB1s1    => plt_morph%NB1s1    , &
    FRADPs1  => plt_rad%FRADPs1      &
  )

  call PrepUptake(PSIST1,WTRTG,VOLPU,VOLWU)
!
!     IF PLANT SPECIES EXISTS
!
  DO NZ=1,NPs1
    OSTRN=0.0_r8
    OSTRD=0.0_r8
    IF(IFLGCs1(NZ).EQ.1.AND.PPs1(NZ).GT.0.0)THEN

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
      IF((IDAYs1(1,NB1s1(NZ),NZ).NE.0).AND.(ARLFSs1(NZ).GT.ZEROLs1(NZ) &
        .AND.FRADPs1(NZ).GT.0.0_r8).AND.(RTDP1s1(1,1,NZ).GT.SDPTHs1(NZ)+CDPTHZs1(0)))THEN
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
        PSILTs1(NZ)=AMIN1(-1.0E-06,0.667*PSILTs1(NZ))
        EPs1(NZ)=0.0_r8
        EVAPCs1(NZ)=0.0_r8
        HFLWC1=FLWCs1(NZ)*cpw*TKAs1

        IF(ARLSSs1.GT.ZEROSs1)THEN
          FPC=ARLFSs1(NZ)/ARLSSs1*AMIN1(1.0,0.5*ARLFCs1/AREA3s1(NUs1))
        ELSEIF(PPTs1.GT.ZEROSs1)THEN
          FPC=PPs1(NZ)/PPTs1
        ELSE
          FPC=1.0/NPs1
        ENDIF

        TKCX=TKCs1(NZ)
        WVPLT=AMAX1(0.0,WTLSs1(NZ)+WVSTKs1(NZ))
        VHCPX=cpw*(WVPLT*VSTK+VOLWCs1(NZ)+VOLWPs1(NZ))
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
        TKCs1(NZ)=TKCZs1(NZ)
        TCCs1(NZ)=TKCs1(NZ)-TC2K
        DTKCs1(NZ)=TKCs1(NZ)-TKAs1
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

      TLECs1=TLECs1+EFLXCs1(NZ)*RAs1(NZ)
      TSHCs1=TSHCs1+SFLXCs1(NZ)*RAs1(NZ)
      IF(OSTRD.GT.ZEROPs1(NZ))THEN
        OSTRs1(NZ)=OSTRN/OSTRD
      ELSE
        OSTRs1(NZ)=0.0_r8
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
    ZEROs1   => plt_site%ZEROs1     , &
    NPs1     => plt_site%NPs1       , &
    NJs1     => plt_site%NJs1       , &
    VOLWMs1  => plt_site%VOLWMs1    , &
    ALTs1    => plt_site%ALTs1      , &
    NUs1     => plt_site%NUs1       , &
    NP0s1    => plt_site%NP0s1      , &
    PSISTs1  => plt_ew%PSISTs1      , &
    WTRTDs1  => plt_biom%WTRTDs1    , &
    VOLAs1   => plt_soilchem%VOLAs1 , &
    VOLIs1   => plt_soilchem%VOLIs1 , &
    THETYs1  => plt_soilchem%THETYs1, &
    BKDSs1   => plt_soilchem%BKDSs1 , &
    VOLWs1   => plt_soilchem%VOLWs1 , &
    VOLYs1   => plt_soilchem%VOLYs1 , &
    ARSTPs1  => plt_morph%ARSTPs1   , &
    MYs1     => plt_morph%MYs1      , &
    ARLFPs1  => plt_morph%ARLFPs1   , &
    RAD1s1   => plt_rad%RAD1s1      , &
    THRM1s1  => plt_rad%THRM1s1       &
  )
!
!     RESET TOTAL UPTAKE ARRAYS
!
!     ARLFP,ARSTP=leaf,stalk areas
!

  ARLSC=0.0_r8
  DO 9984 NZ=1,NP0s1
!     TKCs1(NZ)=TKAs1+DTKCs1(NZ)
!     TCCs1(NZ)=TKCs1(NZ)-TC2K
    ARLSC=ARLSC+ARLFPs1(NZ)+ARSTPs1(NZ)
    RAD1s1(NZ)=0.0_r8
    plt_ew%EFLXCs1(NZ)=0.0_r8
    plt_ew%SFLXCs1(NZ)=0.0_r8
    plt_ew%HFLXCs1(NZ)=0.0_r8
    THRM1s1(NZ)=0.0_r8
    plt_ew%EPs1(NZ)=0.0_r8
    plt_ew%EVAPCs1(NZ)=0.0_r8
    plt_rbgc%UPOMCs1(NZ)=0.0_r8
    plt_rbgc%UPOMNs1(NZ)=0.0_r8
    plt_rbgc%UPOMPs1(NZ)=0.0_r8
    plt_rbgc%UPNH4s1(NZ)=0.0_r8
    plt_rbgc%UPNO3s1(NZ)=0.0_r8
    plt_rbgc%UPH2Ps1(NZ)=0.0_r8
    plt_rbgc%UPH1Ps1(NZ)=0.0_r8
    plt_rbgc%UPNFs1(NZ)=0.0_r8
!
!     RESET UPTAKE ARRAYS
!
    DO  L=NUs1,NJs1
      DO  N=1,MYs1(NZ)
        plt_ew%UPWTRs1(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2Ps1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXPs1(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2Ss1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXSs1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPCHSs1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN2Ss1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3Ss1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3Bs1(N,L,NZ)=0.0_r8
        plt_rbgc%RUPHGSs1(N,L,NZ)=0.0_r8
        plt_rbgc%RCOFLAs1(N,L,NZ)=0.0_r8
        plt_rbgc%ROXFLAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RCHFLAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RN2FLAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RNHFLAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RHGFLAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RCODFAs1(N,L,NZ)=0.0_r8
        plt_rbgc%ROXDFAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RCHDFAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RN2DFAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RNHDFAs1(N,L,NZ)=0.0_r8
        plt_rbgc%RHGDFAs1(N,L,NZ)=0.0_r8
      enddo
    enddo
9984  CONTINUE
!
!     PSIST1=total soil water potential PSIST adjusted for surf elevn
!     ALT=surface elevation
!     VOLWU,VOLWMs1=water volume available for uptake,total water volume
!     THETY,VOLX=hygroscopic SWC,soil volume
!     VOLPU=air volume
!     WTRTG=total biome root mass
!
  DO 9000 L=NUs1,NJs1
    PSIST1(L)=PSISTs1(L)-0.0098_r8*ALTs1
    IF(BKDSs1(L).GT.ZEROs1)THEN
      VOLWU(L)=VOLWMs1(NPH,L)-THETYs1(L)*VOLYs1(L)
      VOLPU(L)=AMAX1(0.0,VOLAs1(L)-VOLWs1(L)-VOLIs1(L))
    ELSE
      VOLWU(L)=VOLWMs1(NPH,L)
      VOLPU(L)=0.0_r8
    ENDIF
    WTRTG(L)=0.0_r8
    DO 9005 NZ=1,NPs1
      DO  N=1,MYs1(NZ)
!     IF(IFLGCs1(NZ).EQ.1.AND.PPs1(NZ).GT.0.0)THEN
      WTRTG(L)=WTRTG(L)+AMAX1(0.0,WTRTDs1(N,L,NZ))
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
    ZEROs1   =>  plt_site%ZEROs1    , &
    NPs1     =>  plt_site%NPs1      , &
    IETYPs1  =>  plt_site%IETYPs1   , &
    RABs1    =>  plt_ew%RABs1       , &
    TKAs1    =>  plt_ew%TKAs1       , &
    DTKCs1   =>  plt_ew%DTKCs1      , &
    ZDs1     =>  plt_ew%ZDs1        , &
    ZRs1     =>  plt_ew%ZRs1        , &
    RAZs1    =>  plt_ew%RAZs1       , &
    TKCZs1   =>  plt_ew%TKCZs1      , &
    WSLFs1   =>  plt_biom%WSLFs1    , &
    ZEROPs1  =>  plt_biom%ZEROPs1   , &
    FRADPs1  =>  plt_rad%FRADPs1    , &
    SURFXs1  =>  plt_photo%SURFXs1  , &
    ARLF1s1  =>  plt_morph%ARLF1s1  , &
    KLEAFXs1 =>  plt_morph%KLEAFXs1 , &
    ZCs1     =>  plt_morph%ZCs1     , &
    CFXs1    =>  plt_morph%CFXs1    , &
    ARLFSs1  =>  plt_morph%ARLFSs1  , &
    NBRs1    =>  plt_morph%NBRs1    , &
    SURFs1   =>  plt_morph%SURFs1   , &
    ZTs1     =>  plt_morph%ZTs1       &
  )
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!
  DO 500 NB=1,NBRs1(NZ)
    DO 550 K=1,JNODS1
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     ARLF=leaf area
!     WSLF=leaf protein content
!     SURFX,SURF=unself-shaded,total leaf surface area
!     CFX=clumping factor from PFT file
!
      IF(ARLF1s1(K,NB,NZ).GT.ZEROPs1(NZ) &
        .AND.WSLFs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
        KLEAFXs1(NB,NZ)=K
      ENDIF
      DO 600 L=JC1,1,-1
        DO 650 N=1,JLI1
          SURFXs1(N,L,K,NB,NZ)=SURFs1(N,L,K,NB,NZ)*CFXs1(NZ)
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
  IF(ARLFSs1(NZ).GT.0.0_r8)THEN
    IF(IETYPs1.GE.0)THEN
      TFRADP=0.0_r8
      DO 700 NZZ=1,NPs1
        IF(ZCs1(NZZ).GT.ZCs1(NZ)+ZRs1)THEN
          TFRADP=TFRADP+FRADPs1(NZZ)
        ENDIF
700   CONTINUE
      ALFZ=2.0_r8*TFRADP
      IF(RABs1.GT.ZEROs1.AND.ZTs1.GT.ZEROs1.AND.ALFZ.GT.ZEROs1)THEN
        RACZ(NZ)=AMIN1(RACX,AMAX1(0.0,ZTs1*EXP(ALFZ) &
          /(ALFZ/RABs1)*(EXP(-ALFZ*ZCs1(NZ)/ZTs1) &
          -EXP(-ALFZ*(ZDs1+ZRs1)/ZTs1))))
      ELSE
        RACZ(NZ)=0.0_r8
      ENDIF
    ELSE
      RACZ(NZ)=0.0_r8
    ENDIF
  ELSE
    RACZ(NZ)=RACX
  ENDIF
  RAZs1(NZ)=RABs1+RACZ(NZ)
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
  TKCZs1(NZ)=TKAs1+DTKCs1(NZ)
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
    WTRTDs1  =>  plt_biom%WTRTDs1    , &
    FMPRs1   =>  plt_site%FMPRs1     , &
    CDPTHZs1 =>  plt_site%CDPTHZs1   , &
    DLYR3s1  =>  plt_site%DLYR3s1    , &
    ZEROSs1  =>  plt_site%ZEROSs1    , &
    ZEROs1   =>  plt_site%ZEROs1     , &
    PPs1     =>  plt_site%PPs1       , &
    NUs1     =>  plt_site%NUs1       , &
    NRTs1    =>  plt_morph%NRTs1     , &
    MYs1     =>  plt_morph%MYs1      , &
    RTDNPs1  =>  plt_morph%RTDNPs1   , &
    RRAD2Xs1 =>  plt_morph%RRAD2Xs1  , &
    HTCTLs1  =>  plt_morph%HTCTLs1   , &
    RTLGPs1  =>  plt_morph%RTLGPs1   , &
    RTDP1s1  =>  plt_morph%RTDP1s1   , &
    PORTs1   =>  plt_morph%PORTs1    , &
    RTVLWs1  =>  plt_morph%RTVLWs1   , &
    RRAD2Ms1 =>  plt_morph%RRAD2Ms1  , &
    SDPTHs1  =>  plt_morph%SDPTHs1   , &
    NIs1     =>  plt_morph%NIs1        &
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
  DO 2000 N=1,MYs1(NZ)
    DO  L=NUs1,NIs1(NZ)
      IF(N.EQ.1)THEN
        RTDPZ=0.0_r8
        DO 2005 NR=1,NRTs1(NZ)
          RTDPZ=AMAX1(RTDPZ,RTDP1s1(1,NR,NZ))
2005    CONTINUE
        IF(L.EQ.NUs1)THEN
          FRTDPX(L,NZ)=1.0_r8
        ELSE
          IF(DLYR3s1(L).GT.ZEROs1)THEN
            RTDPX=AMAX1(0.0_r8,RTDPZ-CDPTHZs1(L-1))
            RTDPX=AMAX1(0.0_r8,AMIN1(DLYR3s1(L),RTDPX) &
              -AMAX1(0.0,SDPTHs1(NZ)-CDPTHZs1(L-1)-HTCTLs1(NZ)))
            FRTDPX(L,NZ)=RTDPX/DLYR3s1(L)
          ELSE
            FRTDPX(L,NZ)=0.0_r8
          ENDIF
        ENDIF

      ENDIF
      IF(WTRTG(L).GT.ZEROSs1)THEN
        FPQ(N,L,NZ)=AMAX1(0.0,WTRTDs1(N,L,NZ))/WTRTG(L)
      ELSE
        FPQ(N,L,NZ)=1.0
      ENDIF
      FPP(N,L,NZ)=FMN*FPQ(N,L,NZ)
      IF(RTDNPs1(N,L,NZ).GT.ZEROs1.AND.FRTDPX(L,NZ).GT.ZEROs1)THEN
        RRADL(N,L)=AMAX1(RRAD2Xs1(N,NZ),SQRT((RTVLWs1(N,L,NZ) &
          /(1.0-PORTs1(N,NZ)))/(PICON*PPs1(NZ)*RTLGPs1(N,L,NZ))))
        PATH(N,L)=AMAX1(1.001*RRADL(N,L) &
          ,1.0/(SQRT(PICON*(RTDNPs1(N,L,NZ)/FRTDPX(L,NZ))/FMPRs1(L))))
        RTARR(N,L)=6.283*RTLGPs1(N,L,NZ)/FRTDPX(L,NZ)
      ELSE
        RRADL(N,L)=RRAD2Ms1(N,NZ)
        PATH(N,L)=1.001*RRADL(N,L)
        RTARR(N,L)=6.283*RTLGPs1(N,L,NZ)
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
   RAZs1     => plt_ew%RAZs1       , &
   DTKCs1    => plt_ew%DTKCs1      , &
   OSMOs1    => plt_ew%OSMOs1      , &
   PSILOs1   => plt_ew%PSILOs1     , &
   TKCs1     => plt_ew%TKCs1       , &
   TKAs1     => plt_ew%TKAs1       , &
   TKSs1     => plt_ew%TKSs1       , &
   PSIROs1   => plt_ew%PSIROs1     , &
   TCCs1     => plt_ew%TCCs1       , &
   UPWTRs1   => plt_ew%UPWTRs1     , &
   PSIRGs1   => plt_ew%PSIRGs1     , &
   PSILGs1   => plt_ew%PSILGs1     , &
   PSIRTs1   => plt_ew%PSIRTs1     , &
   VHCPCs1   => plt_ew%VHCPCs1     , &
   PSILTs1   => plt_ew%PSILTs1     , &
   NUs1      => plt_site%NUs1      , &
   AREA3s1   => plt_site%AREA3s1   , &
   IYRCs1    => plt_site%IYRCs1    , &
   NIs1      => plt_morph%NIs1     , &
   NGs1      => plt_morph%NGs1     , &
   MYs1      => plt_morph%MYs1     , &
   CCPOLRs1  => plt_biom%CCPOLRs1  , &
   CZPOLRs1  => plt_biom%CZPOLRs1  , &
   CPPOLRs1  => plt_biom%CPPOLRs1  , &
   CCPOLPs1  => plt_biom%CCPOLPs1  , &
   CZPOLPs1  => plt_biom%CZPOLPs1  , &
   CPPOLPs1  => plt_biom%CPPOLPs1  , &
   WTSHTs1   => plt_biom%WTSHTs1   , &
   RAs1      => plt_photo%RAs1     , &
   RSMNs1    => plt_photo%RSMNs1   , &
   RCs1      => plt_photo%RCs1     , &
   RSMHs1    => plt_photo%RSMHs1   , &
   RCSs1     => plt_photo%RCSs1    , &
   FRADPs1   => plt_rad%FRADPs1    , &
   THRM1s1   => plt_rad%THRM1s1      &
  )
  IF(NN.GE.MXN)THEN
    WRITE(*,9999)IYRCs1,I,J,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)
    IF(DIFF.GT.0.5)THEN
      plt_rad%RAD1s1(NZ)=0.0_r8
      plt_ew%EFLXCs1(NZ)=0.0_r8
      plt_ew%SFLXCs1(NZ)=0.0_r8
      plt_ew%HFLXCs1(NZ)=0.0_r8
      plt_ew%EVAPCs1(NZ)=0.0_r8
      plt_ew%EPs1(NZ)=0.0_r8
      TKCs1(NZ)=TKAs1+DTKCs1(NZ)
      TCCs1(NZ)=TKCs1(NZ)-TC2K
      FTHRM=EMMC*2.04E-10*FRADPs1(NZ)*AREA3s1(NUs1)
      THRM1s1(NZ)=FTHRM*TKCs1(NZ)**4
      PSILTs1(NZ)=PSIST1(NGs1(NZ))
      APSILT=ABS(PSILTs1(NZ))
      FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLPs1(NZ)+CZPOLPs1(NZ)+CPPOLPs1(NZ)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILOs1(NZ)=FDMP/0.16*OSMOs1(NZ)-RGAS*TKCs1(NZ)*FDMP*CCPOLT/OSWT
      PSILGs1(NZ)=AMAX1(0.0,PSILTs1(NZ)-PSILOs1(NZ))
      WFNC=EXP(RCSs1(NZ)*PSILGs1(NZ))
      RCs1(NZ)=RSMNs1(NZ)+(RSMHs1(NZ)-RSMNs1(NZ))*WFNC
      RAs1(NZ)=RAZs1(NZ)
      VHCPCs1(NZ)=cpw*(WTSHTs1(NZ)*10.0E-06)
      DTKCs1(NZ)=0.0_r8
      DO 4290 N=1,MYs1(NZ)
        DO  L=NUs1,NIs1(NZ)
          PSIRTs1(N,L,NZ)=PSIST1(L)
          APSIRT=ABS(PSIRTs1(N,L,NZ))
          FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
          CCPOLT=CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)
          OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
          PSIROs1(N,L,NZ)=FDMR/0.16*OSMOs1(NZ) &
            -RGAS*TKSs1(L)*FDMR*CCPOLT/OSWT
          PSIRGs1(N,L,NZ)=AMAX1(0.0,PSIRTs1(N,L,NZ)-PSIROs1(N,L,NZ))
          UPWTRs1(N,L,NZ)=0.0_r8
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
    NUs1      => plt_site%NUs1      , &
    AREA3s1   => plt_site%AREA3s1   , &
    PSILOs1   => plt_ew%PSILOs1     , &
    OSMOs1    => plt_ew%OSMOs1      , &
    RAZs1     => plt_ew%RAZs1       , &
    VPAs1     => plt_ew%VPAs1       , &
    TKAs1     => plt_ew%TKAs1       , &
    RIBs1     => plt_ew%RIBs1       , &
    EPs1      => plt_ew%EPs1        , &
    PSILTs1   => plt_ew%PSILTs1     , &
    VOLWPs1   => plt_ew%VOLWPs1     , &
    UPWTRs1   => plt_ew%UPWTRs1     , &
    TKCZs1    => plt_ew%TKCZs1      , &
    EVAPCs1   => plt_ew%EVAPCs1     , &
    VHCPCs1   => plt_ew%VHCPCs1     , &
    PSILGs1   => plt_ew%PSILGs1     , &
    VOLWCs1   => plt_ew%VOLWCs1     , &
    FLWCs1    => plt_ew%FLWCs1      , &
    EFLXCs1   => plt_ew%EFLXCs1     , &
    ZEROLs1   => plt_biom%ZEROLs1   , &
    ZEROPs1   => plt_biom%ZEROPs1   , &
    CCPOLPs1  => plt_biom%CCPOLPs1  , &
    CZPOLPs1  => plt_biom%CZPOLPs1  , &
    CPPOLPs1  => plt_biom%CPPOLPs1  , &
    NIs1      => plt_morph%NIs1     , &
    MYs1      => plt_morph%MYs1     , &
    RSMNs1    => plt_photo%RSMNs1   , &
    RCSs1     => plt_photo%RCSs1    , &
    RAs1      => plt_photo%RAs1     , &
    RSMHs1    => plt_photo%RSMHs1   , &
    RCs1      => plt_photo%RCs1     , &
    RAD1s1    => plt_rad%RAD1s1     , &
    RADCs1    => plt_rad%RADCs1     , &
    THSs1     => plt_rad%THSs1      , &
    FRADPs1   => plt_rad%FRADPs1    , &
    THRM1s1   => plt_rad%THRM1s1    , &
    THRMGXs1  => plt_rad%THRMGXs1     &
  )
  CCPOLT=CCPOLPs1(NZ)+CZPOLPs1(NZ)+CPPOLPs1(NZ)
  OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
  FTHRM=EMMC*2.04E-10*FRADPs1(NZ)*AREA3s1(NUs1)
  FDTHS=(THSs1+THRMGXs1)*FRADPs1(NZ)
!     RAZ=canopy isothermal boundary later resistance

  UPRT=0.0_r8
  PAREX=FPC*AREA3s1(NUs1)
  PARHX=FPC*AREA3s1(NUs1)*1.25E-03
  RA1=RAZs1(NZ)

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
    TKC1=TKCZs1(NZ)
    THRM1s1(NZ)=FTHRM*TKC1**4
    DTHS1=FDTHS-THRM1s1(NZ)*2.0
    RAD1s1(NZ)=RADCs1(NZ)+DTHS1
!
!     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
!
!     RI=Ricardson's number
!     RA=canopy boundary layer resistance
!     PAREC,PARHC=canopy latent,sensible heat conductance
!
    RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIBs1*(TKAs1-TKC1)))

    RAs1(NZ)=AMAX1(RACM,0.9_r8*RA1,AMIN1(1.1_r8*RA1,RAZs1(NZ)/(1.0_r8-10.0_r8*RI)))
    RA1=RAs1(NZ)
    PAREC=PAREX/RAs1(NZ)
    PARHC=PARHX/RAs1(NZ)
!
!     CANOPY WATER AND OSMOTIC POTENTIALS
!
!     PSILT=canopy total water potential
!     FDMP=dry matter content
!     OSMO=osmotic potential at PSILT=0 from PFT file
!     PSILO,PSILG=canopy osmotic,turgor water potential
!
    APSILT=ABS(PSILTs1(NZ))
    FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
    PSILOs1(NZ)=FDMP/0.16*OSMOs1(NZ)-RGAS*TKC1*FDMP*CCPOLT/OSWT
    PSILGs1(NZ)=AMAX1(0.0,PSILTs1(NZ)-PSILOs1(NZ))
!
!     CANOPY STOMATAL RESISTANCE
!
!     RCS=shape parameter for RC vs PSILG from PFT file
!     RC=canopy stomatal resistance
!     RSMN=minimum RC at PSILT=0 from stomate.f
!     RSMX=cuticular resistance from PFT file
!
    WFNC=EXP(RCSs1(NZ)*PSILGs1(NZ))
    RCs1(NZ)=RSMNs1(NZ)+(RSMHs1(NZ)-RSMNs1(NZ))*WFNC
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
    vpc=vapsat(tkc1)*EXP(18.0*PSILTs1(NZ)/(RGAS*TKC1))
    EX=PAREC*(VPAs1-VPC)
    IF(EX.GT.0.0)THEN
      EVAPCs1(NZ)=EX*RAs1(NZ)/(RAs1(NZ)+RZ)
      EX=0.0_r8
    ELSEIF(EX.LE.0.0.AND.VOLWCs1(NZ).GT.0.0)THEN
      EVAPCs1(NZ)=AMAX1(EX*RAs1(NZ)/(RAs1(NZ)+RZ),-VOLWCs1(NZ))
      EX=EX-EVAPCs1(NZ)
    ENDIF
    EPs1(NZ)=EX*RAs1(NZ)/(RAs1(NZ)+RCs1(NZ))
    EFLXCs1(NZ)=(EPs1(NZ)+EVAPCs1(NZ))*VAP
    VFLXC=EVAPCs1(NZ)*cpw*TKC1
!
!     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
!
!     HFLXS=initial estimate of sensible+storage heat flux
!     HFLWC1=convective heat flux from precip to canopy
!
    HFLXS=RAD1s1(NZ)+EFLXCs1(NZ)+VFLXC+HFLWC1
!
!     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
!
!     VHCPC=canopy heat capacity
!     TKCY=equilibrium canopy temperature for HFLXS
!
    VHCPCs1(NZ)=VHCPX+cpw*(EVAPCs1(NZ)+FLWCs1(NZ))
    TKCY=(TKCX*VHCPX+TKAs1*PARHC+HFLXS)/(VHCPCs1(NZ)+PARHC)
    TKCY=AMIN1(TKAs1+10.0,AMAX1(TKAs1-10.0,TKCY))
!
!     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
!
!     XC,IC=magnitude,direction of change in canopy temp for next cycle
!
    IF((IC.EQ.0.AND.TKCY.GT.TKC1).OR.(IC.EQ.1.AND.TKCY.LT.TKC1))THEN
      XC=0.5*XC
    ENDIF
    TKCZs1(NZ)=TKC1+0.1*(TKCY-TKC1)
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
      PSILC=PSILTs1(NZ)-PSILH
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
      DO N=1,MYs1(NZ)
        DO  L=NUs1,NIs1(NZ)
          IF(ILYR(N,L).EQ.1)THEN
            UPWTRs1(N,L,NZ)=AMAX1(AMIN1(0.0,-VOLWU(L)*FPQ(N,L,NZ)) &
              ,AMIN1((PSILC-PSIST1(L))/RSRS(N,L),VOLPU(L)*FPQ(N,L,NZ)))
            IF(UPWTRs1(N,L,NZ).GT.0.0)THEN
              UPWTRs1(N,L,NZ)=0.1*UPWTRs1(N,L,NZ)
            ENDIF
            UPRT=UPRT+UPWTRs1(N,L,NZ)
          ELSE
            UPWTRs1(N,L,NZ)=0.0_r8
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
      DIFFZ=VOLWPZ-VOLWPs1(NZ)
      DIFFU=EPs1(NZ)-UPRT
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
      IF(ABS(VOLWPZ-VOLWPX).GT.ZEROPs1(NZ))THEN
        RSSZ=ABS((PSILTs1(NZ)-PSIL2)/(VOLWPZ-VOLWPX))
      ELSEIF(CNDT.GT.ZEROPs1(NZ))THEN
        RSSZ=1.0_r8/CNDT
      ELSE
        RSSZ=ZEROLs1(NZ)
      ENDIF
      IF(ABS(EPs1(NZ)-EPX).GT.ZEROPs1(NZ))THEN
        RSSUX=ABS((PSILTs1(NZ)-PSIL2)/(EPs1(NZ)-EPX))
        IF(CNDT.GT.ZEROPs1(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(ABS(UPRT-UPRTX).GT.ZEROPs1(NZ))THEN
        RSSUX=ABS((PSILTs1(NZ)-PSIL2)/(UPRT-UPRTX))
        IF(CNDT.GT.ZEROPs1(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(CNDT.GT.ZEROPs1(NZ))THEN
        RSSU=1.0_r8/CNDT
      ELSE
        RSSU=ZEROLs1(NZ)
      ENDIF
!
!     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
!     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
!     WATER STORAGE
!
!     DPSI=change in PSILT for next convergence cycle
!     1.0E-03=acceptance criterion for DPSI
!
      DPSI=AMIN1(AMIN1(RSSZ,RSSU)*(DIFFU-DIFFZ),ABS(PSILTs1(NZ)))

!     IF CONVERGENCE CRITERION IS MET THEN FINISH,
!     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
!     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
!
      IF(.not.((NN.GE.30.AND.ABS(DPSI).LT.1.0E-03_r8).OR.NN.GE.MXN))then
        PSIL2=PSILTs1(NZ)
        EPX=EPs1(NZ)
        UPRTX=UPRT
        VOLWPX=VOLWPZ
        PSILTs1(NZ)=AMIN1(0.0_r8,PSILTs1(NZ)+0.5_r8*DPSI)
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
    DPTHZs1  => plt_site%DPTHZs1    , &
    PPs1     => plt_site%PPs1       , &
    VOLWMs1  => plt_site%VOLWMs1    , &
    ZEROs1   => plt_site%ZEROs1     , &
    ZEROS2s1 => plt_site%ZEROS2s1   , &
    NUs1     => plt_site%NUs1       , &
    PSILTs1  => plt_ew%PSILTs1      , &
    ZEROPs1  => plt_biom%ZEROPs1    , &
    THETWs1  => plt_soilchem%THETWs1, &
    VOLAs1   => plt_soilchem%VOLAs1 , &
    CNDUs1   => plt_soilchem%CNDUs1 , &
    VOLXs1   => plt_soilchem%VOLXs1 , &
    RTNLs1   => plt_morph%RTNLs1    , &
    RTN1s1   => plt_morph%RTN1s1    , &
    RRAD2Ms1 => plt_morph%RRAD2Ms1  , &
    RSRAs1   => plt_morph%RSRAs1    , &
    RRAD2s1  => plt_morph%RRAD2s1   , &
    RTDNPs1  => plt_morph%RTDNPs1   , &
    HTSTZs1  => plt_morph%HTSTZs1   , &
    RTLGPs1  => plt_morph%RTLGPs1   , &
    RTLGAs1  => plt_morph%RTLGAs1   , &
    RRAD1s1  => plt_morph%RRAD1s1   , &
    MYs1     => plt_morph%MYs1      , &
    RSRRs1   => plt_morph%RSRRs1    , &
    ZCs1     => plt_morph%ZCs1      , &
    NGs1     => plt_morph%NGs1      , &
    NIs1     => plt_morph%NIs1        &
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
  HTSTZs1(NZ)=0.80*ZCs1(NZ)
  PSILH=-0.0098*HTSTZs1(NZ)
  FRADW=1.0E+04*(AMAX1(0.5,1.0+PSILTs1(NZ)/EMODW))**4
!
  !     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
  !
  !      VOLX,VOLWMs1,THETW=soil,water volume,content
  !     RTDNP,RTLGP=root length density,root length per plant
  !     CNDU=soil hydraulic conductivity for root uptake
  !     RTN1,RTNL=number of root,myco primary,secondary axes
  !     ILYR:1=rooted,0=not rooted
  !     N:1=root,2=mycorrhizae
!
  DO 3880 N=1,MYs1(NZ)
    DO  L=NUs1,NIs1(NZ)
      IF(VOLXs1(L).GT.ZEROS2s1 &
        .AND.VOLWMs1(NPH,L).GT.ZEROS2s1 &
        .AND.RTDNPs1(N,L,NZ).GT.ZEROs1 &
        .AND.CNDUs1(L).GT.ZEROs1 &
        .AND.RTN1s1(1,L,NZ).GT.ZEROPs1(NZ) &
        .AND.RTNLs1(N,L,NZ).GT.ZEROPs1(NZ) &
        .AND.THETWs1(L).GT.ZEROs1)THEN
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
        RSSL=(LOG(PATH(N,L)/RRADL(N,L))/RTARR(N,L))/PPs1(NZ)
        RSSX(N,L)=RSSL/CNDUs1(L)
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
        RTAR2=6.283*RRAD2s1(N,L,NZ)*RTLGPs1(N,L,NZ)*PPs1(NZ)
        RSRG(N,L)=RSRRs1(N,NZ)/RTAR2*VOLAs1(L)/VOLWMs1(NPH,L)
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
        FRAD1=(RRAD1s1(N,L,NZ)/RRAD2Ms1(N,NZ))**4
        RSR1(N,L)=RSRAs1(N,NZ)*DPTHZs1(L)/(FRAD1*RTN1s1(1,L,NZ)) &
          +RSRAs1(1,NZ)*HTSTZs1(NZ)/(FRADW*RTN1s1(1,L,NZ))
        FRAD2=(RRAD2s1(N,L,NZ)/RRAD2Ms1(N,NZ))**4
        RSR2(N,L)=RSRAs1(N,NZ)*RTLGAs1(N,L,NZ)/(FRAD2*RTNLs1(N,L,NZ))
      ELSE
        ILYR(N,L)=0
      ENDIF
    enddo
3880  CONTINUE

  DO 3890 N=1,MYs1(NZ)
    DO  L=NUs1,NIs1(NZ)
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
    TCCs1    =>  plt_ew%TCCs1       , &
    OSMOs1   =>  plt_ew%OSMOs1      , &
    TKWs1    =>  plt_ew%TKWs1       , &
    TKCs1    =>  plt_ew%TKCs1       , &
    TKSs1    =>  plt_ew%TKSs1       , &
    TKAs1    =>  plt_ew%TKAs1       , &
    DPTHSs1  =>  plt_ew%DPTHSs1     , &
    DTKCs1   =>  plt_ew%DTKCs1      , &
    RAZs1    =>  plt_ew%RAZs1       , &
    VHCPCs1  =>  plt_ew%VHCPCs1     , &
    PSILOs1  =>  plt_ew%PSILOs1     , &
    PSIRGs1  =>  plt_ew%PSIRGs1     , &
    UPWTRs1  =>  plt_ew%UPWTRs1     , &
    PSILGs1  =>  plt_ew%PSILGs1     , &
    PSILTs1  =>  plt_ew%PSILTs1     , &
    PSIROs1  =>  plt_ew%PSIROs1     , &
    PSIRTs1  =>  plt_ew%PSIRTs1     , &
    EPs1     =>  plt_ew%EPs1        , &
    HFLXCs1  =>  plt_ew%HFLXCs1     , &
    EFLXCs1  =>  plt_ew%EFLXCs1     , &
    SFLXCs1  =>  plt_ew%SFLXCs1     , &
    EVAPCs1  =>  plt_ew%EVAPCs1     , &
    NUs1     =>  plt_site%NUs1      , &
    ZEROs1   =>  plt_site%ZEROs1    , &
    AREA3s1  =>  plt_site%AREA3s1   , &
    CCPOLPs1 =>  plt_biom%CCPOLPs1  , &
    CZPOLPs1 =>  plt_biom%CZPOLPs1  , &
    CPPOLPs1 =>  plt_biom%CPPOLPs1  , &
    CCPOLRs1 =>  plt_biom%CCPOLRs1  , &
    CZPOLRs1 =>  plt_biom%CZPOLRs1  , &
    CPPOLRs1 =>  plt_biom%CPPOLRs1  , &
    WTSHTs1  =>  plt_biom%WTSHTs1   , &
    NIs1     =>  plt_morph%NIs1     , &
    ZCs1     =>  plt_morph%ZCs1     , &
    NGs1     =>  plt_morph%NGs1     , &
    MYs1     =>  plt_morph%MYs1     , &
    RCSs1    =>  plt_photo%RCSs1    , &
    RAs1     =>  plt_photo%RAs1     , &
    RCs1     =>  plt_photo%RCs1     , &
    RSMNs1   =>  plt_photo%RSMNs1   , &
    RSMHs1   =>  plt_photo%RSMHs1   , &
    FRADPs1  =>  plt_rad%FRADPs1    , &
    THRM1s1  =>  plt_rad%THRM1s1    , &
    RAD1s1   =>  plt_rad%RAD1s1       &
  )
  RAD1s1(NZ)=0.0_r8
  EFLXCs1(NZ)=0.0_r8
  SFLXCs1(NZ)=0.0_r8
  HFLXCs1(NZ)=0.0_r8
  EVAPCs1(NZ)=0.0_r8
  EPs1(NZ)=0.0_r8
  IF(ZCs1(NZ).GE.DPTHSs1-ZEROs1)THEN
    TKCs1(NZ)=TKAs1
  ELSE
    TKCs1(NZ)=TKWs1
  ENDIF
  TCCs1(NZ)=TKCs1(NZ)-TC2K
  FTHRM=EMMC*2.04E-10_r8*FRADPs1(NZ)*AREA3s1(NUs1)
  THRM1s1(NZ)=FTHRM*TKCs1(NZ)**4
  PSILTs1(NZ)=PSIST1(NGs1(NZ))
  APSILT=ABS(PSILTs1(NZ))
  FDMP=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
  CCPOLT=CCPOLPs1(NZ)+CZPOLPs1(NZ)+CPPOLPs1(NZ)
  OSWT=36.0_r8+840.0_r8*AMAX1(0.0_r8,CCPOLT)
  PSILOs1(NZ)=FDMP/0.16_r8*OSMOs1(NZ)-RGAS*TKCs1(NZ)*FDMP*CCPOLT/OSWT
  PSILGs1(NZ)=AMAX1(0.0_r8,PSILTs1(NZ)-PSILOs1(NZ))
  WFNC=EXP(RCSs1(NZ)*PSILGs1(NZ))
  RCs1(NZ)=RSMNs1(NZ)+(RSMHs1(NZ)-RSMNs1(NZ))*WFNC
  RAs1(NZ)=RAZs1(NZ)
  VHCPCs1(NZ)=cpw*(WTSHTs1(NZ)*10.0E-06_r8)
  DTKCs1(NZ)=0.0_r8
  DO N=1,MYs1(NZ)
    DO  L=NUs1,NIs1(NZ)
      PSIRTs1(N,L,NZ)=PSIST1(L)
      APSIRT=ABS(PSIRTs1(N,L,NZ))
      FDMR=0.16_r8+0.10_r8*APSIRT/(0.05*APSIRT+2.0_r8)
      CCPOLT=CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)
      OSWT=36.0_r8+840.0_r8*AMAX1(0.0_r8,CCPOLT)
      PSIROs1(N,L,NZ)=FDMR/0.16_r8*OSMOs1(NZ)-RGAS*TKSs1(L)*FDMR*CCPOLT/OSWT
      PSIRGs1(N,L,NZ)=AMAX1(0.0_r8,PSIRTs1(N,L,NZ)-PSIROs1(N,L,NZ))
      UPWTRs1(N,L,NZ)=0.0_r8
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
    NUs1       => plt_site%NUs1      , &
    OSMOs1     => plt_ew%OSMOs1      , &
    TKAs1      => plt_ew%TKAs1       , &
    TKCs1      => plt_ew%TKCs1       , &
    PSILZs1    => plt_ew%PSILZs1     , &
    TKSs1      => plt_ew%TKSs1       , &
    TKCZs1     => plt_ew%TKCZs1      , &
    VOLWCs1    => plt_ew%VOLWCs1     , &
    VOLWPs1    => plt_ew%VOLWPs1     , &
    VHCPCs1    => plt_ew%VHCPCs1     , &
    HFLXCs1    => plt_ew%HFLXCs1     , &
    FLWCs1     => plt_ew%FLWCs1      , &
    EPs1       => plt_ew%EPs1        , &
    EVAPCs1    => plt_ew%EVAPCs1     , &
    PSIRTs1    => plt_ew%PSIRTs1     , &
    SFLXCs1    => plt_ew%SFLXCs1     , &
    PSIROs1    => plt_ew%PSIROs1     , &
    PSIRGs1    => plt_ew%PSIRGs1     , &
    PSILTs1    => plt_ew%PSILTs1     , &
    CCPOLRs1   => plt_biom%CCPOLRs1  , &
    CZPOLRs1   => plt_biom%CZPOLRs1  , &
    CPPOLRs1   => plt_biom%CPPOLRs1  , &
    MYs1       => plt_morph%MYs1     , &
    NIs1       => plt_morph%NIs1       &
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
  VOLWPs1(NZ)=VOLWPs1(NZ)+EPs1(NZ)-UPRT
  VOLWCs1(NZ)=VOLWCs1(NZ)+FLWCs1(NZ)+EVAPCs1(NZ)
  SFLXCs1(NZ)=PARHC*(TKAs1-TKCZs1(NZ))
  HFLXCs1(NZ)=TKCX*VHCPX-TKCZs1(NZ)*VHCPCs1(NZ)+VFLXC+HFLWC1
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
  DO 4505 N=1,MYs1(NZ)
    DO 4510 L=NUs1,NIs1(NZ)
      IF(ILYR(N,L).EQ.1)THEN
        PSIRTs1(N,L,NZ)=AMIN1(0.0,(PSIST1(L)*RSRT(N,L) &
          +PSILTs1(NZ)*RSSX(N,L))/RSRS(N,L))
        APSIRT=ABS(PSIRTs1(N,L,NZ))
        FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
        CCPOLT=CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)
        OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
        PSIROs1(N,L,NZ)=FDMR/0.16*OSMOs1(NZ) &
          -RGAS*TKSs1(L)*FDMR*CCPOLT/OSWT
        PSIRGs1(N,L,NZ)=AMAX1(0.0,PSIRTs1(N,L,NZ)-PSIROs1(N,L,NZ))
      ELSE
        PSIRTs1(N,L,NZ)=PSIST1(L)
        APSIRT=ABS(PSIRTs1(N,L,NZ))
        FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
        CCPOLT=CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)
        OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
        PSIROs1(N,L,NZ)=FDMR/0.16*OSMOs1(NZ) &
          -RGAS*TKSs1(L)*FDMR*CCPOLT/OSWT
        PSIRGs1(N,L,NZ)=AMAX1(0.0,PSIRTs1(N,L,NZ)-PSIROs1(N,L,NZ))
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
    TCCs1    =>  plt_ew%TCCs1       , &
    TKCs1    =>  plt_ew%TKCs1       , &
    TKSs1    =>  plt_ew%TKSs1       , &
    PSILZs1  =>  plt_ew%PSILZs1     , &
    PSILTs1  =>  plt_ew%PSILTs1     , &
    NUs1     =>  plt_site%NUs1      , &
    CHILLs1  =>  plt_photo%CHILLs1  , &
    OFFSTs1  =>  plt_pheno%OFFSTs1  , &
    CTCs1    =>  plt_pheno%CTCs1    , &
    TFN4s1   =>  plt_pheno%TFN4s1   , &
    TCGs1    =>  plt_pheno%TCGs1    , &
    TKGs1    =>  plt_pheno%TKGs1    , &
    IDAYs1   =>  plt_pheno%IDAYs1   , &
    TFN3s1   =>  plt_pheno%TFN3s1   , &
    NIs1     =>  plt_morph%NIs1     , &
    NB1s1    =>  plt_morph%NB1s1      &
  )
  !
  !     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
  !     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
  !
  IF(IDAYs1(1,NB1s1(NZ),NZ).EQ.0)THEN
    TKGs1(NZ)=TKSs1(NUs1)
    !     ELSEIF((IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)
    !    2.AND.IDAYs1(2,NB1s1(NZ),NZ).EQ.0)THEN
    !     TKGs1(NZ)=TKSs1(NUs1)
  ELSE
    TKGs1(NZ)=TKCs1(NZ)
  ENDIF
  TCGs1(NZ)=TKGs1(NZ)-TC2K
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
  TKGO=TKGs1(NZ)+OFFSTs1(NZ)
  RTK=RGAS*TKGO
  STK=710.0*TKGO
  ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
  TFN3s1(NZ)=EXP(25.229-62500/RTK)/ACTV
  DO 100 L=NUs1,NIs1(NZ)
    TKSO=TKSs1(L)+OFFSTs1(NZ)
    RTK=RGAS*TKSO
    STK=710.0*TKSO
    ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
    TFN4s1(L,NZ)=EXP(25.229-62500/RTK)/ACTV
100   CONTINUE
  PSILZs1(NZ)=AMIN1(PSILZs1(NZ),PSILTs1(NZ))
  !
  !     DIURNAL CHILLING
  !
  !     CTC=chilling temperature from PFT file
  !     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
  !
  IF(TCCs1(NZ).LT.CTCs1(NZ))THEN
    CHILLs1(NZ)=AMIN1(24.0,CHILLs1(NZ)+1.0)
  ELSE
    CHILLs1(NZ)=AMAX1(0.0,CHILLs1(NZ)-1.0)
  ENDIF
  end associate
  end subroutine SetCanopyGrowthFuncs

end module UptakesMod
