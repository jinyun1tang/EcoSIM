module UptakesMod
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : GravAcceleration=>DAT_CONST_G
  use StomatesMod   , only : stomates
  use minimathmod   , only : safe_adb,vapsat,isclose,AZMAX1,AZMIN1
  use UnitMod       , only : units  
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

  real(r8), parameter :: mGravAcceleration=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000.  
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
  real(r8) :: PopPlantO2Uptake,PopPlantO2Demand,PARHC
  real(r8) :: WVPLT
  real(r8) :: adjTotalSoilH2OPSIMPa(JZ1)
  real(r8) :: PATH(2,JZ1)
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
    TairK    => plt_ew%TairK       , &
    TLEC   => plt_ew%TLEC      , &
    EP     => plt_ew%EP        , &
    VOLWC  => plt_ew%VOLWC     , &
    CDPTHZ => plt_site%CDPTHZ  , &
    PPT    => plt_site%PPT     , &
    NP     => plt_site%NP      , &
    pftPlantPopulation     => plt_site%pftPlantPopulation      , &
    NU     => plt_site%NU      , &
    AREA3  => plt_site%AREA3   , &
    ZEROS  => plt_site%ZEROS   , &
    WTRTD  => plt_biom%WTRTD   , &
    ZEROL  => plt_biom%ZEROL   , &
    ZEROP  => plt_biom%ZEROP   , &
    WTLS   => plt_biom%WTLS    , &
    WVSTK  => plt_biom%WVSTK   , &
    IDAY   => plt_pheno%IDAY   , &
    PlantO2Stress   => plt_pheno%PlantO2Stress   , &
    IFLGC  => plt_pheno%IFLGC  , &
    ARLFC  => plt_morph%ARLFC  , &
    RTDP1  => plt_morph%RTDP1  , &
    ARLSS  => plt_morph%ARLSS  , &
    ARLFS  => plt_morph%ARLFS  , &
    SDPTH  => plt_morph%SDPTH  , &
    NB1    => plt_morph%NB1    , &
    FRADP  => plt_rad%FRADP      &
  )

  call PrepH2ONutrientUptake(adjTotalSoilH2OPSIMPa,WTRTG,VOLPU,VOLWU)
!

!     IF PLANT SPECIES EXISTS
!
  DO NZ=1,NP
    PopPlantO2Uptake=0.0_r8
    PopPlantO2Demand=0.0_r8

    IF(IFLGC(NZ).EQ.PlantIsActive.AND.pftPlantPopulation(NZ).GT.0.0_r8)THEN

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
        !leaf area > 0, absorped par, and rooting depth > seeding depth
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
        PSILT(NZ)=AMIN1(-ppmc,0.667_r8*PSILT(NZ))
        EP(NZ)=0.0_r8
        EVAPC(NZ)=0.0_r8
        HFLWC1=FLWC(NZ)*cpw*TairK

        IF(ARLSS.GT.ZEROS)THEN
          FPC=ARLFS(NZ)/ARLSS*AMIN1(1.0_r8,0.5_r8*ARLFC/AREA3(NU))
        ELSEIF(PPT.GT.ZEROS)THEN
          FPC=pftPlantPopulation(NZ)/PPT
        ELSE
          FPC=1.0_r8/NP
        ENDIF

        TKCX=TKC(NZ)
        WVPLT=AZMAX1(WTLS(NZ)+WVSTK(NZ))
        VHCPX=cpw*(WVPLT*VSTK+VOLWC(NZ)+VOLWP(NZ))
!
!     CONVERGENCE SOLUTION
!
        NN=CanopyEnergyH2OIteration(I,J,NZ,FPC,WVPLT,&
          adjTotalSoilH2OPSIMPa,PARHC,DIFF,UPRT,VFLXC,FDMP,RSRS,FPQ,VOLPU,&
          VOLWU,TKCX,CNDT,VHCPX,HFLWC1,PSILH,ILYR)
!
!     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
!
!     TKC=final estimate of canopy temperature TKCZ
!     TairK=current air temperature
!     DTKC=TKC-TairK for next hour
!
        TKC(NZ)=TKCZ(NZ)
        TCC(NZ)=units%Kelvin2Celcius(TKC(NZ))
        DTKC(NZ)=TKC(NZ)-TairK
!
!     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
!     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
!
        call HandlingDivergence(I,J,NN,NZ,adjTotalSoilH2OPSIMPa,DIFF,FDMP)

        call UpdateCanopyWater(NZ,PARHC,adjTotalSoilH2OPSIMPa,RSRT,RSSX,RSRS,TKCX,VHCPX,HFLWC1,UPRT,VFLXC,ILYR)
!
!     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
!
      ELSE
        call HandleBareSoil(NZ,adjTotalSoilH2OPSIMPa,FDMP)
      ENDIF

      call SetCanopyGrowthFuncs(NZ)

      call PopPlantNutientO2Uptake(NZ,FDMP,PopPlantO2Uptake,PopPlantO2Demand,PATH,RRADL,FPQ,FPP,FRTDPX,RTARR)

      TLEC=TLEC+EFLXC(NZ)*RA(NZ)
      TSHC=TSHC+SFLXC(NZ)*RA(NZ)
      IF(PopPlantO2Demand.GT.ZEROP(NZ))THEN
        PlantO2Stress(NZ)=PopPlantO2Uptake/PopPlantO2Demand
      ELSE
        PlantO2Stress(NZ)=0.0_r8
      ENDIF
    ENDIF
  ENDDO
  RETURN
  end associate
  END subroutine uptakes
!------------------------------------------------------------------------

  subroutine PrepH2ONutrientUptake(adjTotalSoilH2OPSIMPa,WTRTG,VOLPU,VOLWU)
!
!     prepare for uptake calculation
  implicit none
  real(r8), intent(out) :: adjTotalSoilH2OPSIMPa(JZ1),WTRTG(JZ1),VOLPU(JZ1),VOLWU(JZ1)
  integer :: NZ, L, N
  real(r8) :: ARLSC

  associate(                          &
    ZERO   => plt_site%ZERO     , &
    NP     => plt_site%NP       , &
    NJ     => plt_site%NJ       , &
    VWatMicPM  => plt_site%VWatMicPM    , &
    ALT    => plt_site%ALT      , &
    NU     => plt_site%NU       , &
    NP0    => plt_site%NP0      , &
    TotalSoilH2OPSIMPa  => plt_ew%TotalSoilH2OPSIMPa      , &
    WTRTD  => plt_biom%WTRTD    , &
    VMicP   => plt_soilchem%VMicP , &
    ViceMicP   => plt_soilchem%ViceMicP , &
    THETY  => plt_soilchem%THETY, &
    BKDS   => plt_soilchem%BKDS , &
    VWatMicP   => plt_soilchem%VWatMicP , &
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
  D9984: DO NZ=1,NP0
!     TKC(NZ)=TairK+DTKC(NZ)
!     TCC(NZ)=TKC(NZ)-TC2K
    ARLSC=ARLSC+ARLFP(NZ)+ARSTP(NZ)
    RAD1(NZ)=0.0_r8
    plt_ew%EFLXC(NZ)=0.0_r8
    plt_ew%SFLXC(NZ)=0.0_r8
    plt_ew%HFLXC(NZ)=0.0_r8
    THRM1(NZ)=0.0_r8
    plt_ew%EP(NZ)=0.0_r8
    plt_ew%EVAPC(NZ)=0.0_r8
    plt_rbgc%UPOME(1:npelms,NZ)=0.0_r8
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
        plt_ew%PopPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2P(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXP(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXS(N,L,NZ)=0.0_r8
        plt_rbgc%RUPCHS(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN2S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3B(N,L,NZ)=0.0_r8
        plt_rbgc%RUPHGS(N,L,NZ)=0.0_r8
        plt_rbgc%trcg_RFLA(idg_beg:idg_end-1,N,L,NZ)=0.0_r8
        plt_rbgc%trcg_RDFA(idg_beg:idg_end-1,N,L,NZ)=0.0_r8
      enddo
    enddo
  ENDDO D9984
!
!     adjTotalSoilH2OPSIMPa=total soil water potential PSIST adjusted for surf elevn
!     ALT=surface elevation
!     VOLWU,VWatMicPM=water volume available for uptake,total water volume
!     THETY,VSoilPoreMicP=hygroscopic SWC,soil volume
!     VOLPU=air volume
!     WTRTG=total biome root mass
!
  D9000: DO L=NU,NJ
    adjTotalSoilH2OPSIMPa(L)=TotalSoilH2OPSIMPa(L)-mGravAcceleration*ALT
    IF(BKDS(L).GT.ZERO)THEN
      VOLWU(L)=VWatMicPM(NPH,L)-THETY(L)*VOLY(L)
      VOLPU(L)=AZMAX1(VMicP(L)-VWatMicP(L)-ViceMicP(L))
    ELSE
      VOLWU(L)=VWatMicPM(NPH,L)
      VOLPU(L)=0.0_r8
    ENDIF
    WTRTG(L)=0.0_r8
    D9005: DO NZ=1,NP
      DO  N=1,MY(NZ)
!     IF(IFLGC(NZ).EQ.PlantIsActive.AND.pftPlantPopulation(NZ).GT.0.0)THEN
      WTRTG(L)=WTRTG(L)+AZMAX1(WTRTD(N,L,NZ))
!     ENDIF
      enddo
    ENDDO D9005
  ENDDO D9000  
  end associate
  end subroutine PrepH2ONutrientUptake
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
    TairK    =>  plt_ew%TairK       , &
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
    CanopyHeight     =>  plt_morph%CanopyHeight     , &
    CFX    =>  plt_morph%CFX    , &
    ARLFS  =>  plt_morph%ARLFS  , &
    NBR    =>  plt_morph%NBR    , &
    SURF   =>  plt_morph%SURF   , &
    GridMaxCanopyHeight     =>  plt_morph%GridMaxCanopyHeight       &
  )
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!
  D500: DO NB=1,NBR(NZ)
    D550: DO K=1,JNODS1
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     ARLF=leaf area
!     WSLF=leaf protein content
!     SURFX,SURF=unself-shaded,total leaf surface area
!     CFX=clumping factor from PFT file
!
      IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ).AND.WSLF(K,NB,NZ).GT.ZEROP(NZ))THEN
        KLEAFX(NB,NZ)=K
      ENDIF
      D600: DO L=JC1,1,-1
        D650: DO N=1,JLI1
          SURFX(N,L,K,NB,NZ)=SURF(N,L,K,NB,NZ)*CFX(NZ)
        ENDDO D650
      ENDDO D600
    ENDDO D550 
  ENDDO D500 
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
      D700: DO NZZ=1,NP
        IF(CanopyHeight(NZZ).GT.CanopyHeight(NZ)+ZR)THEN
          TFRADP=TFRADP+FRADP(NZZ)
        ENDIF
      ENDDO D700
      ALFZ=2.0_r8*TFRADP
      IF(RAB.GT.ZERO.AND.GridMaxCanopyHeight.GT.ZERO.AND.ALFZ.GT.ZERO)THEN
        RACZ(NZ)=AMIN1(RACX,AZMAX1(GridMaxCanopyHeight*EXP(ALFZ) &
          /(ALFZ/RAB)*(EXP(-ALFZ*CanopyHeight(NZ)/GridMaxCanopyHeight) &
          -EXP(-ALFZ*(ZD+ZR)/GridMaxCanopyHeight))))
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
!     TairK=current air temperature
!     DTKC=TKC-TairK from previous hour
!
  TKCZ(NZ)=TairK+DTKC(NZ)
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
    pftPlantPopulation     =>  plt_site%pftPlantPopulation       , &
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
  D2000: DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(N.EQ.1)THEN
        RTDPZ=0.0_r8
        D2005: DO NR=1,NRT(NZ)
          RTDPZ=AMAX1(RTDPZ,RTDP1(1,NR,NZ))
        ENDDO D2005
        IF(L.EQ.NU)THEN
          FRTDPX(L,NZ)=1.0_r8
        ELSE
          IF(DLYR3(L).GT.ZERO)THEN
            RTDPX=AZMAX1(RTDPZ-CDPTHZ(L-1))
            RTDPX=AZMAX1(AMIN1(DLYR3(L),RTDPX)-AZMAX1(SDPTH(NZ)-CDPTHZ(L-1)-HTCTL(NZ)))
            FRTDPX(L,NZ)=RTDPX/DLYR3(L)
          ELSE
            FRTDPX(L,NZ)=0.0_r8
          ENDIF
        ENDIF

      ENDIF
      IF(WTRTG(L).GT.ZEROS)THEN
        FPQ(N,L,NZ)=AZMAX1(WTRTD(N,L,NZ))/WTRTG(L)
      ELSE
        FPQ(N,L,NZ)=1.0_r8
      ENDIF
      FPP(N,L,NZ)=FMN*FPQ(N,L,NZ)
      IF(RTDNP(N,L,NZ).GT.ZERO.AND.FRTDPX(L,NZ).GT.ZERO)THEN
        RRADL(N,L)=AMAX1(RRAD2X(N,NZ),SQRT((RTVLW(N,L,NZ) &
          /(1.0_r8-PORT(N,NZ)))/(PICON*pftPlantPopulation(NZ)*RTLGP(N,L,NZ))))
        PATH(N,L)=AMAX1(1.001_r8*RRADL(N,L),1.0_r8/(SQRT(PICON*(RTDNP(N,L,NZ)/FRTDPX(L,NZ))/FMPR(L))))
        RTARR(N,L)=PICON2s*RTLGP(N,L,NZ)/FRTDPX(L,NZ)
      ELSE
        RRADL(N,L)=RRAD2M(N,NZ)
        PATH(N,L)=1.001_r8*RRADL(N,L)
        RTARR(N,L)=PICON2s*RTLGP(N,L,NZ)
      ENDIF
    enddo
  ENDDO D2000
  end associate
  end subroutine UpdateRootProperty
!------------------------------------------------------------------------

  subroutine HandlingDivergence(I,J,NN,NZ,adjTotalSoilH2OPSIMPa,DIFF,FDMP)

  implicit none
  integer  , intent(in) :: NN, I, J
  integer  , intent(in) :: NZ
  REAL(R8) ,INTENT(IN) :: adjTotalSoilH2OPSIMPa(JZ1),DIFF
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
   TairK     => plt_ew%TairK       , &
   TKS     => plt_ew%TKS       , &
   PSIRO   => plt_ew%PSIRO     , &
   TCC     => plt_ew%TCC       , &
   PopPlantRootH2OUptake_vr   => plt_ew%PopPlantRootH2OUptake_vr     , &
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
   CEPOLR  => plt_biom%CEPOLR  , &
   CEPOLP  => plt_biom%CEPOLP  , &
   WTSHTE  => plt_biom%WTSHTE  , &
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
      TKC(NZ)=TairK+DTKC(NZ)
      TCC(NZ)=units%Kelvin2Celcius(TKC(NZ))
      FTHRM=EMMC*2.04E-10_r8*FRADP(NZ)*AREA3(NU)
      THRM1(NZ)=FTHRM*TKC(NZ)**4._r8
      PSILT(NZ)=adjTotalSoilH2OPSIMPa(NG(NZ))
      APSILT=ABS(PSILT(NZ))
      FDMP=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
      CCPOLT=CEPOLP(ielmc,NZ)+CEPOLP(ielmn,NZ)+CEPOLP(ielmp,NZ)
      OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
      PSILO(NZ)=FDMP/0.16*OSMO(NZ)-RGAS*TKC(NZ)*FDMP*CCPOLT/OSWT
      PSILG(NZ)=AZMAX1(PSILT(NZ)-PSILO(NZ))
      WFNC=EXP(RCS(NZ)*PSILG(NZ))
      RC(NZ)=RSMN(NZ)+(RSMH(NZ)-RSMN(NZ))*WFNC
      RA(NZ)=RAZ(NZ)
      VHCPC(NZ)=cpw*(WTSHTE(ielmc,NZ)*10.0E-06_r8)
      DTKC(NZ)=0.0_r8
      D4290: DO N=1,MY(NZ)
        DO  L=NU,NI(NZ)
          PSIRT(N,L,NZ)=adjTotalSoilH2OPSIMPa(L)
          APSIRT=ABS(PSIRT(N,L,NZ))
          FDMR=0.16_r8+0.10_r8*APSIRT/(0.05_r8*APSIRT+2.0_r8)
          CCPOLT=sum(CEPOLR(1:npelms,N,L,NZ))
          OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
          PSIRO(N,L,NZ)=FDMR/0.16_r8*OSMO(NZ)-RGAS*TKS(L)*FDMR*CCPOLT/OSWT
          PSIRG(N,L,NZ)=AZMAX1(PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
          PopPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
      enddo
      ENDDO D4290
    ENDIF
  ENDIF
  end associate
  end subroutine HandlingDivergence
!------------------------------------------------------------------------------
  function CanopyEnergyH2OIteration(I,J,NZ,FPC,WVPLT,&
    adjTotalSoilH2OPSIMPa,PARHC,DIFF,UPRT,VFLXC,FDMP,RSRS,FPQ,VOLPU,VOLWU,TKCX,CNDT,&
    VHCPX,HFLWC1,PSILH,ILYR) result(NN)
  implicit none
  integer  , intent(in) :: I, J
  integer  , intent(in) :: NZ
  real(r8) , intent(in) :: FPC,WVPLT,adjTotalSoilH2OPSIMPa(JZ1)
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
    TairK     => plt_ew%TairK       , &
    RIB     => plt_ew%RIB       , &
    EP      => plt_ew%EP        , &
    PSILT   => plt_ew%PSILT     , &
    VOLWP   => plt_ew%VOLWP     , &
    PopPlantRootH2OUptake_vr   => plt_ew%PopPlantRootH2OUptake_vr     , &
    TKCZ    => plt_ew%TKCZ      , &
    EVAPC   => plt_ew%EVAPC     , &
    VHCPC   => plt_ew%VHCPC     , &
    PSILG   => plt_ew%PSILG     , &
    VOLWC   => plt_ew%VOLWC     , &
    FLWC    => plt_ew%FLWC      , &
    EFLXC   => plt_ew%EFLXC     , &
    ZEROL   => plt_biom%ZEROL   , &
    ZEROP   => plt_biom%ZEROP   , &
    CEPOLP  => plt_biom%CEPOLP  , &
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
  CCPOLT=CEPOLP(ielmc,NZ)+CEPOLP(ielmn,NZ)+CEPOLP(ielmp,NZ)
  OSWT=36.0+840.0*AZMAX1(CCPOLT)
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
    THRM1(NZ)=FTHRM*TKC1**4._r8   !long wave radiation
    DTHS1=FDTHS-THRM1(NZ)*2.0_r8  !upper and down
    RAD1(NZ)=RADC(NZ)+DTHS1
!
!     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
!
!     RI=Ricardson's number
!     RA=canopy boundary layer resistance
!     PAREC,PARHC=canopy latent,sensible heat conductance
!
    RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB*(TairK-TKC1)))

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
    FDMP=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
    PSILO(NZ)=FDMP/0.16_r8*OSMO(NZ)-RGAS*TKC1*FDMP*CCPOLT/OSWT
    PSILG(NZ)=AZMAX1(PSILT(NZ)-PSILO(NZ))  !turgor pressure
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
!     PAREC=aerodynamic conductance

    VPC=vapsat(tkc1)*EXP(18.0_r8*PSILT(NZ)/(RGAS*TKC1))
    EX=PAREC*(VPA-VPC)   !evaporation demand
    IF(EX.GT.0.0_r8)THEN
      !condensation  > 0._r8, to canopy
      EVAPC(NZ)=EX*RA(NZ)/(RA(NZ)+RZ)
      EX=0.0_r8
    ELSEIF(EX.LE.0.0_r8.AND.VOLWC(NZ).GT.0.0_r8)THEN
      !evaporation, and there is water stored in canopy
      !<0._r8, off canopy
      EVAPC(NZ)=AMAX1(EX*RA(NZ)/(RA(NZ)+RZ),-VOLWC(NZ))
      EX=EX-EVAPC(NZ)
    ENDIF

    EP(NZ)=EX*RA(NZ)/(RA(NZ)+RC(NZ))
    EFLXC(NZ)=(EP(NZ)+EVAPC(NZ))*VAP   !latent heat flux
    VFLXC=EVAPC(NZ)*cpw*TKC1           !enthalpy of evaporated water
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
!   VHCPX= canopy heat capacity
    VHCPC(NZ)=VHCPX+cpw*(EVAPC(NZ)+FLWC(NZ))
    TKCY=(TKCX*VHCPX+TairK*PARHC+HFLXS)/(VHCPC(NZ)+PARHC)
    TKCY=AMIN1(TairK+10.0_r8,AMAX1(TairK-10.0_r8,TKCY))
!
!     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
!
!     XC,IC=magnitude,direction of change in canopy temp for next cycle
!
    IF((IC.EQ.0.AND.TKCY.GT.TKC1).OR.(IC.EQ.1.AND.TKCY.LT.TKC1))THEN
      XC=0.5_r8*XC
    ENDIF

    TKCZ(NZ)=TKC1+0.1_r8*(TKCY-TKC1)
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
!     PopPlantRootH2OUptake_vr=root water uptake from soil layer
!     VOLWU,VOLPU=water volume available for uptake,air volume
!     FPQ=PFT fraction of biome root mass
!     PSILC=canopy water potential adjusted for canopy height
!     adjTotalSoilH2OPSIMPa=total soil water potential PSIST adjusted for surf elevn
!     RSRS=total soil+root resistance
!     UPRT=total water uptake from soil profile
!
      DO N=1,MY(NZ)
        DO  L=NU,NI(NZ)
          IF(ILYR(N,L).EQ.1)THEN
            PopPlantRootH2OUptake_vr(N,L,NZ)=AMAX1(AZMIN1(-VOLWU(L)*FPQ(N,L,NZ)) &
              ,AMIN1((PSILC-adjTotalSoilH2OPSIMPa(L))/RSRS(N,L),VOLPU(L)*FPQ(N,L,NZ)))
            IF(PopPlantRootH2OUptake_vr(N,L,NZ).GT.0.0_r8)THEN
              PopPlantRootH2OUptake_vr(N,L,NZ)=0.1*PopPlantRootH2OUptake_vr(N,L,NZ)
            ENDIF
            UPRT=UPRT+PopPlantRootH2OUptake_vr(N,L,NZ)
          ELSE
            PopPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
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

      VOLWPZ=ppmc*WVPLT/FDMP
      DIFFZ=VOLWPZ-VOLWP(NZ)
      DIFFU=EP(NZ)-UPRT
      IF(.not.isclose(UPRT,0.0_r8))THEN
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
        PSILT(NZ)=AZMIN1(PSILT(NZ)+0.5_r8*DPSI)
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
    pftPlantPopulation     => plt_site%pftPlantPopulation       , &
    VWatMicPM  => plt_site%VWatMicPM    , &
    ZERO   => plt_site%ZERO     , &
    ZEROS2 => plt_site%ZEROS2   , &
    NU     => plt_site%NU       , &
    PSILT  => plt_ew%PSILT      , &
    ZEROP  => plt_biom%ZEROP    , &
    THETW  => plt_soilchem%THETW, &
    VMicP   => plt_soilchem%VMicP , &
    CNDU   => plt_soilchem%CNDU , &
    VSoilPoreMicP   => plt_soilchem%VSoilPoreMicP , &
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
    CanopyHeight     => plt_morph%CanopyHeight      , &
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
  HTSTZ(NZ)=0.80*CanopyHeight(NZ)
  PSILH=-mGravAcceleration*HTSTZ(NZ)
  FRADW=1.0E+04_r8*(AMAX1(0.5_r8,1.0_r8+PSILT(NZ)/EMODW))**4._r8
!
  !     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
  !
  !      VSoilPoreMicP,VWatMicPM,THETW=soil,water volume,content
  !     RTDNP,RTLGP=root length density,root length per plant
  !     CNDU=soil hydraulic conductivity for root uptake
  !     RTN1,RTNL=number of root,myco primary,secondary axes
  !     ILYR:1=rooted,0=not rooted
  !     N:1=root,2=mycorrhizae
!
  DO 3880 N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(VSoilPoreMicP(L).GT.ZEROS2 &
        .AND.VWatMicPM(NPH,L).GT.ZEROS2 &
        .AND.RTDNP(N,L,NZ).GT.ZERO &
        .AND.CNDU(L).GT.ZERO &
        .AND.RTN1(ipltroot,L,NZ).GT.ZEROP(NZ) &
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
        RSSL=(LOG(PATH(N,L)/RRADL(N,L))/RTARR(N,L))/pftPlantPopulation(NZ)
        RSSX(N,L)=RSSL/CNDU(L)
        !
        !     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
        !     ENTERED IN 'READQ'
        !
        !     RRAD2=secondary root radius
        !     RTLGP=root length per plant
        !     RSRG=radial resistance
        !     RSRR=radial resistivity from PFT file
        !     VMicP,VWatMicPM=soil micropore,water volume
        !
        RTAR2=PICON2s*RRAD2(N,L,NZ)*RTLGP(N,L,NZ)*pftPlantPopulation(NZ)
        RSRG(N,L)=RSRR(N,NZ)/RTAR2*VMicP(L)/VWatMicPM(NPH,L)
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
        FRAD1=(RRAD1(N,L,NZ)/RRAD2M(N,NZ))**4._r8
        RSR1(N,L)=RSRA(N,NZ)*DPTHZ(L)/(FRAD1*RTN1(ipltroot,L,NZ)) &
          +RSRA(ipltroot,NZ)*HTSTZ(NZ)/(FRADW*RTN1(ipltroot,L,NZ))
        FRAD2=(RRAD2(N,L,NZ)/RRAD2M(N,NZ))**4._r8
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

  subroutine HandleBareSoil(NZ,adjTotalSoilH2OPSIMPa,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: adjTotalSoilH2OPSIMPa(JZ1)
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
    TairK    =>  plt_ew%TairK       , &
    SnowDepth  =>  plt_ew%SnowDepth     , &
    DTKC   =>  plt_ew%DTKC      , &
    RAZ    =>  plt_ew%RAZ       , &
    VHCPC  =>  plt_ew%VHCPC     , &
    PSILO  =>  plt_ew%PSILO     , &
    PSIRG  =>  plt_ew%PSIRG     , &
    PopPlantRootH2OUptake_vr  =>  plt_ew%PopPlantRootH2OUptake_vr     , &
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
    CEPOLP =>  plt_biom%CEPOLP  , &
    CEPOLR =>  plt_biom%CEPOLR  , &
    WTSHTE =>  plt_biom%WTSHTE  , &
    NI     =>  plt_morph%NI     , &
    CanopyHeight     =>  plt_morph%CanopyHeight     , &
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
  IF(CanopyHeight(NZ).GE.SnowDepth-ZERO)THEN
    TKC(NZ)=TairK
  ELSE
    TKC(NZ)=TKW
  ENDIF
  TCC(NZ)=units%Kelvin2Celcius(TKC(NZ))
  FTHRM=EMMC*2.04E-10_r8*FRADP(NZ)*AREA3(NU)
  THRM1(NZ)=FTHRM*TKC(NZ)**4._r8
  PSILT(NZ)=adjTotalSoilH2OPSIMPa(NG(NZ))
  APSILT=ABS(PSILT(NZ))
  FDMP=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
  CCPOLT=CEPOLP(ielmc,NZ)+CEPOLP(ielmn,NZ)+CEPOLP(ielmp,NZ)
  OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
  PSILO(NZ)=FDMP/0.16_r8*OSMO(NZ)-RGAS*TKC(NZ)*FDMP*CCPOLT/OSWT
  PSILG(NZ)=AZMAX1(PSILT(NZ)-PSILO(NZ))
  WFNC=EXP(RCS(NZ)*PSILG(NZ))
  RC(NZ)=RSMN(NZ)+(RSMH(NZ)-RSMN(NZ))*WFNC
  RA(NZ)=RAZ(NZ)
  VHCPC(NZ)=cpw*(WTSHTE(ielmc,NZ)*10.0E-06_r8)
  DTKC(NZ)=0.0_r8
  DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      PSIRT(N,L,NZ)=adjTotalSoilH2OPSIMPa(L)
      APSIRT=ABS(PSIRT(N,L,NZ))
      FDMR=0.16_r8+0.10_r8*APSIRT/(0.05_r8*APSIRT+2.0_r8)
      CCPOLT=sum(CEPOLR(1:npelms,N,L,NZ))
      OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
      PSIRO(N,L,NZ)=FDMR/0.16_r8*OSMO(NZ)-RGAS*TKS(L)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ)=AZMAX1(PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      PopPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
    enddo
  ENDDO
  end associate
  end subroutine HandleBareSoil
!------------------------------------------------------------------------

  subroutine UpdateCanopyWater(NZ,PARHC,adjTotalSoilH2OPSIMPa,RSRT,RSSX,&
    RSRS,TKCX,VHCPX,HFLWC1,UPRT,VFLXC,ILYR)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: PARHC,adjTotalSoilH2OPSIMPa(JZ1),RSRT(2,JZ1)
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
    TairK      => plt_ew%TairK       , &
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
    CEPOLR   => plt_biom%CEPOLR  , &
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
  SFLXC(NZ)=PARHC*(TairK-TKCZ(NZ))
  HFLXC(NZ)=TKCX*VHCPX-TKCZ(NZ)*VHCPC(NZ)+VFLXC+HFLWC1
  !
  !     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
  !
  !     PSIRT,PSILT=root,canopy total water potential
  !     adjTotalSoilH2OPSIMPa=total soil water potential PSIST adjusted for surf elevn
  !     RSSX,RSRS,RSRT=soil,soil+root,root radial+axial resistance
  !     PSIRO,PSIRG=root osmotic,turgor water potential
  !     FDMR=dry matter content
  !     OSMO=osmotic potential at PSIRT=0 from PFT file
!
  !
  D4505: DO N=1,MY(NZ)
    D4510: DO L=NU,NI(NZ)
      IF(ILYR(N,L).EQ.1)THEN
        PSIRT(N,L,NZ)=AZMIN1((adjTotalSoilH2OPSIMPa(L)*RSRT(N,L) &
          +PSILT(NZ)*RSSX(N,L))/RSRS(N,L))
        APSIRT=ABS(PSIRT(N,L,NZ))
        FDMR=0.16_r8+0.10_r8*APSIRT/(0.05_r8*APSIRT+2.0_r8)
        CCPOLT=sum(CEPOLR(1:npelms,N,L,NZ))
        OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
        PSIRO(N,L,NZ)=FDMR/0.16_r8*OSMO(NZ)-RGAS*TKS(L)*FDMR*CCPOLT/OSWT
        PSIRG(N,L,NZ)=AZMAX1(PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      ELSE
        PSIRT(N,L,NZ)=adjTotalSoilH2OPSIMPa(L)
        APSIRT=ABS(PSIRT(N,L,NZ))
        FDMR=0.16_r8+0.10_r8*APSIRT/(0.05_r8*APSIRT+2.0_r8)
        CCPOLT=sum(CEPOLR(1:npelms,N,L,NZ))
        OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
        PSIRO(N,L,NZ)=FDMR/0.16_r8*OSMO(NZ)-RGAS*TKS(L)*FDMR*CCPOLT/OSWT
        PSIRG(N,L,NZ)=AZMAX1(PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      ENDIF
    ENDDO D4510
  ENDDO D4505
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
  TCG(NZ)=units%Kelvin2Celcius(TKG(NZ))
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
    CHILL(NZ)=AMIN1(24.0_r8,CHILL(NZ)+1.0_r8)
  ELSE
    CHILL(NZ)=AZMAX1(CHILL(NZ)-1.0_r8)
  ENDIF
  end associate
  end subroutine SetCanopyGrowthFuncs

end module UptakesMod
