module UptakesMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,test_aneb
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
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
    TCCs1    => plt_ew%TCCs1       , &
    TKCs1    => plt_ew%TKCs1       , &
    TSHCs1   => plt_ew%TSHCs1      , &
    DTKCs1   => plt_ew%DTKCs1      , &
    TKAs1    => plt_ew%TKAs1       , &
    TLECs1   => plt_ew%TLECs1      , &
    WTRTDs1  => plt_biom%WTRTDs1   , &
    WTLSs1   => plt_biom%WTLSs1    , &
    WVSTKs1  => plt_biom%WVSTKs1   , &
    IDAYs1   => plt_pheno%IDAYs1   , &
    IFLGCs1  => plt_pheno%IFLGCs1  , &
    ARLFCs1  => plt_morph%ARLFCs1  , &
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

      call CanopyNH3Flux(NZ,FDMP)
!
!     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
!
      call RootMycoO2NutrientUptake(NZ,OSTRN,OSTRD,PATH,RRADL,&
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
    WTRTDs1  => plt_biom%WTRTDs1    , &
    VOLAs1   => plt_soilchem%VOLAs1 , &
    PSISTs1  => plt_soilchem%PSISTs1, &
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
    EFLXCs1(NZ)=0.0_r8
    SFLXCs1(NZ)=0.0_r8
    HFLXCs1(NZ)=0.0_r8
    THRM1s1(NZ)=0.0_r8
    EPs1(NZ)=0.0_r8
    EVAPCs1(NZ)=0.0_r8
    UPOMCs1(NZ)=0.0_r8
    UPOMNs1(NZ)=0.0_r8
    UPOMPs1(NZ)=0.0_r8
    UPNH4s1(NZ)=0.0_r8
    UPNO3s1(NZ)=0.0_r8
    UPH2Ps1(NZ)=0.0_r8
    UPH1Ps1(NZ)=0.0_r8
    UPNFs1(NZ)=0.0_r8
!
!     RESET UPTAKE ARRAYS
!
    DO  L=NUs1,NJs1
      DO  N=1,MYs1(NZ)
        UPWTRs1(N,L,NZ)=0.0_r8
        RCO2Ps1(N,L,NZ)=0.0_r8
        RUPOXPs1(N,L,NZ)=0.0_r8
        RCO2Ss1(N,L,NZ)=0.0_r8
        RUPOXSs1(N,L,NZ)=0.0_r8
        RUPCHSs1(N,L,NZ)=0.0_r8
        RUPN2Ss1(N,L,NZ)=0.0_r8
        RUPN3Ss1(N,L,NZ)=0.0_r8
        RUPN3Bs1(N,L,NZ)=0.0_r8
        RUPHGSs1(N,L,NZ)=0.0_r8
        RCOFLAs1(N,L,NZ)=0.0_r8
        ROXFLAs1(N,L,NZ)=0.0_r8
        RCHFLAs1(N,L,NZ)=0.0_r8
        RN2FLAs1(N,L,NZ)=0.0_r8
        RNHFLAs1(N,L,NZ)=0.0_r8
        RHGFLAs1(N,L,NZ)=0.0_r8
        RCODFAs1(N,L,NZ)=0.0_r8
        ROXDFAs1(N,L,NZ)=0.0_r8
        RCHDFAs1(N,L,NZ)=0.0_r8
        RN2DFAs1(N,L,NZ)=0.0_r8
        RNHDFAs1(N,L,NZ)=0.0_r8
        RHGDFAs1(N,L,NZ)=0.0_r8
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
    RABs1    =>  plt_ew%RABs1       , &
    TKAs1    =>  plt_ew%TKAs1       , &
    DTKCs1   =>  plt_ew%DTKCs1      , &
    ZDs1     =>  plt_ew%ZDs1        , &
    ZRs1     =>  plt_ew%ZRs1        , &
    RAZs1    =>  plt_ew%RAZs1       , &
    WSLFs1   =>  plt_biom%WSLFs1    , &
    FRADPs1  =>  plt_rad%FRADPs1    , &
    ARLF1s1  =>  plt_morph%ARLF1s1  , &
    KLEAFXs1 =>  plt_morph%KLEAFXs1 , &
    ZCs1     => plt_morph%ZCs1      , &
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

  associate(                          &
    WTRTDs1  => plt_biom%WTRTDs1    , &
    NRTs1    =>  plt_morph%NRTs1    , &
    MYs1     => plt_morph%MYs1      , &
    SDPTHs1  =>  plt_morph%SDPTHs1  , &
    NIs1     =>  plt_morph%NIs1       &
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
          IF(DLYR3S1(L).GT.ZEROs1)THEN
            RTDPX=AMAX1(0.0_r8,RTDPZ-CDPTHZs1(L-1))
            RTDPX=AMAX1(0.0_r8,AMIN1(DLYR3S1(L),RTDPX) &
              -AMAX1(0.0,SDPTHs1(NZ)-CDPTHZs1(L-1)-HTCTLs1(NZ)))
            FRTDPX(L,NZ)=RTDPX/DLYR3S1(L)
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
   TCCs1     => plt_ew%TCCs1       , &
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
   FRADPs1   => plt_rad%FRADPs1    , &
   RSMHs1    => plt_photo%RSMHs1   , &
   RCSs1     => plt_photo%RCSs1    , &
   RAD1s1    => plt_rad%RAD1s1     , &
   THRM1s1   => plt_rad%THRM1s1      &
  )
  IF(NN.GE.MXN)THEN
    WRITE(*,9999)IYRCs1,I,J,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)
    IF(DIFF.GT.0.5)THEN
      RAD1s1(NZ)=0.0_r8
      EFLXCs1(NZ)=0.0_r8
      SFLXCs1(NZ)=0.0_r8
      HFLXCs1(NZ)=0.0_r8
      EVAPCs1(NZ)=0.0_r8
      EPs1(NZ)=0.0_r8
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
    PSILOs1   => plt_ew%PSILOs1     , &
    OSMOs1    => plt_ew%OSMOs1      , &
    RAZs1     => plt_ew%RAZs1       , &
    VPAs1     => plt_ew%VPAs1       , &
    TKAs1     => plt_ew%TKAs1       , &
    RIBs1     => plt_ew%RIBs1       , &
    CCPOLPs1  => plt_biom%CCPOLPs1  , &
    CZPOLPs1  => plt_biom%CZPOLPs1  , &
    CPPOLPs1  => plt_biom%CPPOLPs1  , &
    NIs1       => plt_morph%NIs1    , &
    MYs1       => plt_morph%MYs1    , &
    RCSs1      => plt_photo%RCSs1   , &
    RSMHs1     => plt_photo%RSMHs1  , &
    RAD1s1     => plt_rad%RAD1s1    , &
    RADCs1     => plt_rad%RADCs1    , &
    THSs1      => plt_rad%THSs1     , &
    FRADPs1    => plt_rad%FRADPs1   , &
    THRM1s1    => plt_rad%THRM1s1   , &
    THRMGXs1   => plt_rad%THRMGXs1    &
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
    THETWs1  => plt_soilchem%THETWs1, &
    VOLAs1   => plt_soilchem%VOLAs1 , &
    CNDUs1   => plt_soilchem%CNDUs1 , &
    VOLXs1   => plt_soilchem%VOLXs1 , &
    HTSTZs1  => plt_morph%HTSTZs1   , &
    MYs1     => plt_morph%MYs1      , &
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
    TCCs1    => plt_ew%TCCs1       , &
    OSMOs1   => plt_ew%OSMOs1      , &
    TKWs1    => plt_ew%TKWs1       , &
    TKCs1    => plt_ew%TKCs1       , &
    TKSs1    => plt_ew%TKSs1       , &
    TKAs1    => plt_ew%TKAs1       , &
    DPTHSs1  => plt_ew%DPTHSs1     , &
    DTKCs1   => plt_ew%DTKCs1      , &
    RAZs1    => plt_ew%RAZs1       , &
    PSILOs1  => plt_ew%PSILOs1     , &
    CCPOLPs1 =>  plt_biom%CCPOLPs1 , &
    CZPOLPs1 =>  plt_biom%CZPOLPs1 , &
    CPPOLPs1 =>  plt_biom%CPPOLPs1 , &
    CCPOLRs1 =>  plt_biom%CCPOLRs1 , &
    CZPOLRs1 =>  plt_biom%CZPOLRs1 , &
    CPPOLRs1 =>  plt_biom%CPPOLRs1 , &
    WTSHTs1  =>  plt_biom%WTSHTs1  , &
    NIs1     =>  plt_morph%NIs1    , &
    ZCs1     =>  plt_morph%ZCs1    , &
    NGs1     =>  plt_morph%NGs1    , &
    MYs1     =>  plt_morph%MYs1    , &
    RCSs1    =>  plt_photo%RCSs1   , &
    RSMHs1   =>  plt_photo%RSMHs1  , &
    FRADPs1  =>  plt_rad%FRADPs1   , &
    THRM1s1  =>  plt_rad%THRM1s1   , &
    RAD1s1   =>  plt_rad%RAD1s1      &
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
    OSMOs1     => plt_ew%OSMOs1      , &
    TKAs1      => plt_ew%TKAs1       , &
    TKCs1      => plt_ew%TKCs1       , &
    PSILZs1    => plt_ew%PSILZs1     , &
    TKSs1      => plt_ew%TKSs1       , &
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
    TCCs1    => plt_ew%TCCs1        , &
    TKCs1    => plt_ew%TKCs1        , &
    TKSs1    => plt_ew%TKSs1        , &
    PSILZs1  => plt_ew%PSILZs1      , &
    CHILLs1  =>  plt_photo%CHILLs1  , &
    IDAYs1   =>  plt_pheno%IDAYs1   , &
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
!------------------------------------------------------------------------

  subroutine CanopyNH3Flux(NZ,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP
  real(r8) :: FNH3P
  real(r8) :: CNH3P
  real(r8) :: SNH3P
  real(r8) :: ZPOOLB
  integer :: NB

  associate(                            &
    TCCs1     => plt_ew%TCCs1         , &
    WTLSBs1   =>  plt_biom%WTLSBs1    , &
    ZPOOLs1   =>  plt_biom%ZPOOLs1    , &
    CCPOLBs1  =>  plt_biom%CCPOLBs1   , &
    CZPOLBs1  =>  plt_biom%CZPOLBs1   , &
    CPPOLBs1  =>  plt_biom%CPPOLBs1   , &
    ARLFBs1   =>  plt_morph%ARLFBs1   , &
    NBRs1     =>  plt_morph%NBRs1     , &
    FRADPs1   =>  plt_rad%FRADPs1     , &
    ARLFPs1   =>  plt_morph%ARLFPs1     &
  )
  !
  !     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
  !     CONCENTRATION DIFFERENCES 'CNH3E' (ATMOSPHERE FROM 'READS') AND
  !     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
  !
  !     SNH3P,SNH3X=NH3 solubility at TCC, 25 oC
  !     TCC=canopy temperature (oC)
  !     FDMP,FNH3P=canopy dry matter content,NH3 concentration
  !     ARLFB,ARLFP=branch,canopy leaf area
  !     CNH3P,CNH3E=gaseous NH3 concentration in branch,atmosphere
  !     CZPOLB,ZPOOLB=nonstructural N concentration,content in branch
  !     RNH3B=NH3 flux between atmosphere and branch
  !     RA,RC=canopy boundary layer,stomatal resistance
  !     FRADP=fraction of radiation received by each PFT canopy
  !
  SNH3P=SNH3X*EXP(0.513-0.0171*TCCs1(NZ))
  FNH3P=1.0E-04*FDMP
  DO 105 NB=1,NBRs1(NZ)
    IF(WTLSBs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.ARLFBs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.ARLFPs1(NZ).GT.ZEROPs1(NZ))THEN
      CNH3P=AMAX1(0.0,FNH3P*CZPOLBs1(NB,NZ)/SNH3P)
      ZPOOLB=AMAX1(0.0,ZPOOLs1(NB,NZ))
      RNH3Bs1(NB,NZ)=AMIN1(0.1*ZPOOLB,AMAX1((CNH3Es1-CNH3P)/(RAs1(NZ)+RCs1(NZ)) &
      *FRADPs1(NZ)*AREA3s1(NUs1)*ARLFBs1(NB,NZ)/ARLFPs1(NZ),-0.1*ZPOOLB))
      ELSE
      RNH3Bs1(NB,NZ)=0.0_r8
      ENDIF

105   CONTINUE
  end associate
  end subroutine CanopyNH3Flux
!------------------------------------------------------------------------

  subroutine ZeroUptake(N,L,NZ)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ

  integer :: K
  !     begin_execution

  RCOFLAs1(N,L,NZ)=0.0_r8
  ROXFLAs1(N,L,NZ)=0.0_r8
  RCHFLAs1(N,L,NZ)=0.0_r8
  RN2FLAs1(N,L,NZ)=0.0_r8
  RNHFLAs1(N,L,NZ)=0.0_r8
  RCODFAs1(N,L,NZ)=0.0_r8
  ROXDFAs1(N,L,NZ)=0.0_r8
  RCHDFAs1(N,L,NZ)=0.0_r8
  RN2DFAs1(N,L,NZ)=0.0_r8
  RNHDFAs1(N,L,NZ)=0.0_r8
  RCO2Ss1(N,L,NZ)=0.0_r8
  RUPOXSs1(N,L,NZ)=0.0_r8
  RUPCHSs1(N,L,NZ)=0.0_r8
  RUPN2Ss1(N,L,NZ)=0.0_r8
  RUPN3Ss1(N,L,NZ)=0.0_r8
  RCO2Ps1(N,L,NZ)=0.0_r8
  RUPOXPs1(N,L,NZ)=0.0_r8
  DO 395 K=0,jcplx11
    RDFOMCs1(N,K,L,NZ)=0.0_r8
    RDFOMNs1(N,K,L,NZ)=0.0_r8
    RDFOMPs1(N,K,L,NZ)=0.0_r8
395   CONTINUE
  WFRs1(N,L,NZ)=1.0
  RUNNHPs1(N,L,NZ)=0.0_r8
  RUPNH4s1(N,L,NZ)=0.0_r8
  RUONH4s1(N,L,NZ)=0.0_r8
  RUCNH4s1(N,L,NZ)=0.0_r8
  RUNNBPs1(N,L,NZ)=0.0_r8
  RUPNHBs1(N,L,NZ)=0.0_r8
  RUONHBs1(N,L,NZ)=0.0_r8
  RUCNHBs1(N,L,NZ)=0.0_r8
  RUNNOPs1(N,L,NZ)=0.0_r8
  RUPNO3s1(N,L,NZ)=0.0_r8
  RUONO3s1(N,L,NZ)=0.0_r8
  RUCNO3s1(N,L,NZ)=0.0_r8
  RUNNXPs1(N,L,NZ)=0.0_r8
  RUPNOBs1(N,L,NZ)=0.0_r8
  RUONOBs1(N,L,NZ)=0.0_r8
  RUCNOBs1(N,L,NZ)=0.0_r8
  RUPP2Ps1(N,L,NZ)=0.0_r8
  RUPH2Ps1(N,L,NZ)=0.0_r8
  RUOH2Ps1(N,L,NZ)=0.0_r8
  RUCH2Ps1(N,L,NZ)=0.0_r8
  RUPP2Bs1(N,L,NZ)=0.0_r8
  RUPH2Bs1(N,L,NZ)=0.0_r8
  RUOH2Bs1(N,L,NZ)=0.0_r8
  RUCH2Bs1(N,L,NZ)=0.0_r8
  RUPP1Ps1(N,L,NZ)=0.0_r8
  RUPH1Ps1(N,L,NZ)=0.0_r8
  RUOH1Ps1(N,L,NZ)=0.0_r8
  RUCH1Ps1(N,L,NZ)=0.0_r8
  RUPP1Bs1(N,L,NZ)=0.0_r8
  RUPH1Bs1(N,L,NZ)=0.0_r8
  RUOH1Bs1(N,L,NZ)=0.0_r8
  RUCH1Bs1(N,L,NZ)=0.0_r8
  IF(N.EQ.1)RUPNFs1(L,NZ)=0.0_r8
  end subroutine ZeroUptake
!------------------------------------------------------------------------
  subroutine NoActiveNutrientUptake(N,L,NZ)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ
!     begin_execution

  RUNNHPs1(N,L,NZ)=0.0_r8
  RUPNH4s1(N,L,NZ)=0.0_r8
  RUONH4s1(N,L,NZ)=0.0_r8
  RUCNH4s1(N,L,NZ)=0.0_r8
  RUNNBPs1(N,L,NZ)=0.0_r8
  RUPNHBs1(N,L,NZ)=0.0_r8
  RUONHBs1(N,L,NZ)=0.0_r8
  RUCNHBs1(N,L,NZ)=0.0_r8
  RUNNOPs1(N,L,NZ)=0.0_r8
  RUPNO3s1(N,L,NZ)=0.0_r8
  RUONO3s1(N,L,NZ)=0.0_r8
  RUCNO3s1(N,L,NZ)=0.0_r8
  RUNNXPs1(N,L,NZ)=0.0_r8
  RUPNOBs1(N,L,NZ)=0.0_r8
  RUONOBs1(N,L,NZ)=0.0_r8
  RUCNOBs1(N,L,NZ)=0.0_r8
  RUPP2Ps1(N,L,NZ)=0.0_r8
  RUPH2Ps1(N,L,NZ)=0.0_r8
  RUOH2Ps1(N,L,NZ)=0.0_r8
  RUCH2Ps1(N,L,NZ)=0.0_r8
  RUPP2Bs1(N,L,NZ)=0.0_r8
  RUPH2Bs1(N,L,NZ)=0.0_r8
  RUOH2Bs1(N,L,NZ)=0.0_r8
  RUCH2Bs1(N,L,NZ)=0.0_r8
  RUPP1Ps1(N,L,NZ)=0.0_r8
  RUPH1Ps1(N,L,NZ)=0.0_r8
  RUOH1Ps1(N,L,NZ)=0.0_r8
  RUCH1Ps1(N,L,NZ)=0.0_r8
  RUPP1Bs1(N,L,NZ)=0.0_r8
  RUPH1Bs1(N,L,NZ)=0.0_r8
  RUOH1Bs1(N,L,NZ)=0.0_r8
  RUCH1Bs1(N,L,NZ)=0.0_r8
  end subroutine NoActiveNutrientUptake
!------------------------------------------------------------------------
  subroutine UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FP14X,FP1BX
  real(r8), intent(in) :: FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFH1P
  real(r8) :: DIFH1B
  real(r8) :: H1POM,H1POX,H1PXM,H1PXB
  real(r8) :: RMFH1P,RTKH1P,RTKHP1,RTKH1B,RTKHB1
  real(r8) :: RMFH2B
  real(r8) :: UPMX,UPMXP
  real(r8) :: X,Y
  !     begin_execution
  associate(                             &
    VOLWs1    => plt_soilchem%VOLWs1   , &
    VLPO4s1   => plt_soilchem%VLPO4s1  , &
    CH1P4s1   => plt_soilchem%CH1P4s1  , &
    H1PO4s1   => plt_soilchem%H1PO4s1  , &
    VLPOBs1   => plt_soilchem%VLPOBs1  , &
    CH1P4Bs1  => plt_soilchem%CH1P4Bs1 , &
    H1POBs1   => plt_soilchem%H1POBs1    &
  )
  !
  !     HPO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4=HPO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
  !     DIFH1P=soil-root HPO4 diffusion per plant in non-band
!
    IF(VLPO4s1(L).GT.ZEROs1.AND.CH1P4s1(L) &
      .GT.UPMNPOs1(N,NZ))THEN
      RMFH1P=UPWTRP*CH1P4s1(L)*VLPO4s1(L)
      DIFH1P=DIFFL*VLPO4s1(L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max HPO4 uptake in non-band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1*UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLPO4s1(L)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
    !     DIFH1P=soil-root HPO4 diffusion per plant in non-band
    !     CH1P4=HPO4 concentration in non-band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1P,RTKHP1=HPO4 uptake per plant in non-band lmtd,unlmtd by O2
    !     H1POM,H1POX=minimum,maximum HPO4 available for uptake in non-band
    !     FP14X=fraction of total HPO4 uptake in non-band by root,myco populn
    !     RUPP1P,RUPH1P=HPO4 uptake in non-band unlimited,limited by HPO4
    !     RUOH1P=HPO4 uptake in non-band unlimited by O2
    !     RUCH1P=HPO4 uptake in non-band unlimited by nonstructural C
!
    X=(DIFH1P+RMFH1P)*CH1P4s1(L)
    Y=DIFH1P*UPMNPOs1(N,NZ)
    B=-UPMX-DIFH1P*UPKMPOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKH1P=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH1P*UPKMPOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKHP1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H1POM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPO4s1(L)
    H1POX=AMAX1(0.0,FP14X*(H1PO4s1(L)-H1POM))
    RUPP1Ps1(N,L,NZ)=AMAX1(0.0,RTKH1P*PPs1(NZ))
    RUPH1Ps1(N,L,NZ)=AMIN1(H1POX,RUPP1Ps1(N,L,NZ))
    RUOH1Ps1(N,L,NZ)=AMIN1(H1POX,AMAX1(0.0,RTKHP1*PPs1(NZ)))
    RUCH1Ps1(N,L,NZ)=RUPH1Ps1(N,L,NZ)/FCUP
    !     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.3)THEN
    !     WRITE(*,2226)'UPPO4',I,J,NZ,L,N,RUPH2Ps1(N,L,NZ),FPO4X
    !    2,H2PO4s1(L),RUPP2Ps1(N,L,NZ),UPMX,DIFPO,UPKMPOs1(N,NZ)
    !    3,UPMNPOs1(N,NZ),RMFH2P,CH2P4s1(L),UPMXP,WFRs1(N,L,NZ)
    !    4,FCUP,FZUP,FPUP,UPMXPOs1(N,NZ),RTARPs1(N,L,NZ),FWSRT
    !    5,TFN4s1(L,NZ),DIFFL,FH2P,CPO4Ss1(L),CPOOLRs1(N,L,NZ)
    !    6,PPOOLRs1(N,L,NZ),RTKH2P,PPs1(NZ)
    !    2,RTLGPs1(N,L,NZ)
!2226  FORMAT(A8,5I4,40E12.4)
    !     ENDIF
  ELSE
    RUPP1Ps1(N,L,NZ)=0.0_r8
    RUPH1Ps1(N,L,NZ)=0.0_r8
    RUOH1Ps1(N,L,NZ)=0.0_r8
    RUCH1Ps1(N,L,NZ)=0.0_r8
  ENDIF
  !
  !     HPO4 UPTAKE IN BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4B=HPO4 concentration in band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH1B=soil-root convective HPO4 flux per plant in band
  !     DIFH1B=soil-root HPO4 diffusion per plant in band
  !
  IF(VLPOBs1(L).GT.ZEROs1.AND.CH1P4Bs1(L) &
    .GT.UPMNPOs1(N,NZ))THEN
    RMFH2B=UPWTRP*CH1P4Bs1(L)*VLPOBs1(L)
    DIFH1B=DIFFL*VLPOBs1(L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1*UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLPOBs1(L)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1B=soil-root convective HPO4 flux per plant in band
    !     DIFH1B=soil-root HPO4 diffusion per plant in band
    !     CH1P4B=HPO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1B,RTKHB1=HPO4 uptake per plant in band lmtd,unlmtd by O2
    !     H1PXM,H1PXB=minimum,maximum HPO4 available for uptake in band
    !     FP1BX=fraction of total HPO4 uptake in band by root,myco populn
    !     RUPP1B,RUPH1B=HPO4 uptake in band unlimited,limited by H2PO4
    !     RUOH1B=HPO4 uptake in band unlimited by O2
    !     RUCH1B=HPO4 uptake in band unlimited by nonstructural C
!
    X=(DIFH1B+RMFH2B)*CH1P4Bs1(L)
    Y=DIFH1B*UPMNPOs1(N,NZ)
    B=-UPMX-DIFH1B*UPKMPOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKH1B=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH1B*UPKMPOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKHB1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H1PXM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPOBs1(L)
    H1PXB=AMAX1(0.0,FP1BX*(H1POBs1(L)-H1PXM))
    RUPP1Bs1(N,L,NZ)=AMAX1(0.0,RTKH1B*PPs1(NZ))
    RUPH1Bs1(N,L,NZ)=AMIN1(H1PXB,RUPP1Bs1(N,L,NZ))
    RUOH1Bs1(N,L,NZ)=AMIN1(H1PXB,AMAX1(0.0,RTKHB1*PPs1(NZ)))
    RUCH1Bs1(N,L,NZ)=RUPH1Bs1(N,L,NZ)/FCUP
  ELSE
    RUPP1Bs1(N,L,NZ)=0.0_r8
    RUPH1Bs1(N,L,NZ)=0.0_r8
    RUOH1Bs1(N,L,NZ)=0.0_r8
    RUCH1Bs1(N,L,NZ)=0.0_r8
  ENDIF
  end associate
  end subroutine UptakeHPO4
!------------------------------------------------------------------------

  subroutine UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer,  intent(in) :: N,L
  integer,  intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FPO4X,FPOBX
  real(r8), intent(in) :: FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFH2P,DIFH2B
  real(r8) :: H2POM,H2POX,H2PXM,H2PXB
  real(r8) :: RTKHPB,RMFH2B,RMFH2P,RTKH2P,RTKHPP,RTKH2B
  real(r8) :: UPMX,UPMXP,X,Y
  !
  associate(                             &
    H2POBs1   => plt_soilchem%H2POBs1  , &
    VOLWs1    => plt_soilchem%VOLWs1   , &
    H2PO4s1   => plt_soilchem%H2PO4s1  , &
    CH2P4s1   => plt_soilchem%CH2P4s1  , &
    VLPO4s1   => plt_soilchem%VLPO4s1  , &
    CH2P4Bs1  => plt_soilchem%CH2P4Bs1 , &
    VLPOBs1   => plt_soilchem%VLPOBs1    &
  )
  !     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH2P4=H2PO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
  !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
!
    IF(VLPO4s1(L).GT.ZEROs1.AND.CH2P4s1(L) &
      .GT.UPMNPOs1(N,NZ))THEN
      RMFH2P=UPWTRP*CH2P4s1(L)*VLPO4s1(L)
      DIFH2P=DIFFL*VLPO4s1(L)
      !
      !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
      !     AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS CALCULATED ABOVE
      !
      !     UPMXP,UPMX=max H2PO4 uptake in non-band unlimited,limited by O2
      !     RTARP=root surface area per plant from grosub.f
      !     FWSRT=protein concentration relative to 5%
      !     TFN4=temperature function for root growth
      !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
      !     WFR=constraint by O2 consumption on all biological processes
!
      UPMXP=UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
        *FWSRT*TFN4s1(L,NZ)*VLPO4s1(L)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFRs1(N,L,NZ)
      !
      !     SOLUTION FOR MASS FLOW + DIFFUSION OF H2PO4 IN AQUEOUS PHASE OF
      !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY
      !     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
      !
      !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
      !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
      !     CH2P4=H2PO4 concentration in non-band
      !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
      !     RTKH2P,RTKHPP=H2PO4 uptake per plant in non-band lmtd,unlmtd by O2
      !     H2POM,H2POX=minimum,maximum H2PO4 available for uptake in non-band
      !     FPO4X=fraction of total H2PO4 uptake in non-band by root,myco populn
      !     RUPP2P,RUPH2P=H2PO4 uptake in non-band unlimited,limited by H2PO4
      !     RUOH2P=H2PO4 uptake in non-band unlimited by O2
      !     RUCH2P=H2PO4 uptake in non-band unlimited by nonstructural C
!
      X=(DIFH2P+RMFH2P)*CH2P4s1(L)
      Y=DIFH2P*UPMNPOs1(N,NZ)
      B=-UPMX-DIFH2P*UPKMPOs1(N,NZ)-X+Y
      C=(X-Y)*UPMX
      RTKH2P=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH2P*UPKMPOs1(N,NZ)-X+Y
      CP=(X-Y)*UPMXP
      RTKHPP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H2POM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPO4s1(L)
      H2POX=AMAX1(0.0,FPO4X*(H2PO4s1(L)-H2POM))
      RUPP2Ps1(N,L,NZ)=AMAX1(0.0,RTKH2P*PPs1(NZ))
      RUPH2Ps1(N,L,NZ)=AMIN1(H2POX,RUPP2Ps1(N,L,NZ))
      RUOH2Ps1(N,L,NZ)=AMIN1(H2POX,AMAX1(0.0_r8,RTKHPP*PPs1(NZ)))
      RUCH2Ps1(N,L,NZ)=RUPH2Ps1(N,L,NZ)/FCUP
    ELSE
      RUPP2Ps1(N,L,NZ)=0.0_r8
      RUPH2Ps1(N,L,NZ)=0.0_r8
      RUOH2Ps1(N,L,NZ)=0.0_r8
      RUCH2Ps1(N,L,NZ)=0.0_r8
    ENDIF
    !
    !     H2PO4 UPTAKE IN BAND SOIL ZONE
    !
    !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
    !     CH2P4B=H2PO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     UPWTRP=root water uptake per plant
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !

  IF(VLPOBs1(L).GT.ZEROs1.AND.CH2P4Bs1(L) &
    .GT.UPMNPOs1(N,NZ))THEN
    RMFH2B=UPWTRP*CH2P4Bs1(L)*VLPOBs1(L)
    DIFH2B=DIFFL*VLPOBs1(L)
    !
    !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLPOBs1(L)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF PO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !     CH2P4B=H2PO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH2B,RTKHPB=H2PO4 uptake per plant in band lmtd,unlmtd by O2
    !     H2PXM,H2PXB=minimum,maximum H2PO4 available for uptake in band
    !     FPOBX=fraction of total H2PO4 uptake in band by root,myco populn
    !     RUPP2B,RUPH2B=H2PO4 uptake in band unlimited,limited by H2PO4
    !     RUOH2B=H2PO4 uptake in band unlimited by O2
    !     RUCH2B=H2PO4 uptake in band unlimited by nonstructural C
!
    X=(DIFH2B+RMFH2B)*CH2P4Bs1(L)
    Y=DIFH2B*UPMNPOs1(N,NZ)
    B=-UPMX-DIFH2B*UPKMPOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKH2B=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH2B*UPKMPOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKHPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H2PXM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPOBs1(L)
    H2PXB=AMAX1(0.0,FPOBX*(H2POBs1(L)-H2PXM))
    RUPP2Bs1(N,L,NZ)=AMAX1(0.0,RTKH2B*PPs1(NZ))
    RUPH2Bs1(N,L,NZ)=AMIN1(H2PXB,RUPP2Bs1(N,L,NZ))
    RUOH2Bs1(N,L,NZ)=AMIN1(H2PXB,AMAX1(0.0,RTKHPB*PPs1(NZ)))
    RUCH2Bs1(N,L,NZ)=RUPH2Bs1(N,L,NZ)/FCUP
  ELSE
    RUPP2Bs1(N,L,NZ)=0.0_r8
    RUPH2Bs1(N,L,NZ)=0.0_r8
    RUOH2Bs1(N,L,NZ)=0.0_r8
    RUCH2Bs1(N,L,NZ)=0.0_r8
  ENDIF
  end associate
  end subroutine UptakeH2PO4
!------------------------------------------------------------------------

  subroutine UptakeMineralPhosporhus(N,L,NZ,PATH,RRADL,FPQ,FPP,&
    RTARR,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(in):: RTARR(2,JZ1),FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: TFPO4X,TFP14X,TFPOBX,TFP1BX
  real(r8) :: DIFFL
  real(r8) :: FP14X,FP1BX
  real(r8) :: FPO4X,FPOBX
  real(r8) :: POSGX
  real(r8) :: PATHL
  associate(                               &
    POSGLs1   =>   plt_soilchem%POSGLs1    &
  )
  TFPO4X=0.0_r8
  TFPOBX=0.0_r8
  TFP14X=0.0_r8
  TFP1BX=0.0_r8

!     begin_execution
  IF(RPO4Ys1(L).GT.ZEROSs1)THEN
    FPO4X=AMAX1(FPP(N,L,NZ),RUPP2Ps1(N,L,NZ)/RPO4Ys1(L))
  ELSE
    FPO4X=FPQ(N,L,NZ)
  ENDIF
  IF(RPOBYs1(L).GT.ZEROSs1)THEN
    FPOBX=AMAX1(FPP(N,L,NZ),RUPP2Bs1(N,L,NZ)/RPOBYs1(L))
  ELSE
    FPOBX=FPQ(N,L,NZ)
  ENDIF
  IF(RP14Ys1(L).GT.ZEROSs1)THEN
    FP14X=AMAX1(FPP(N,L,NZ),RUPP1Ps1(N,L,NZ)/RP14Ys1(L))
  ELSE
    FP14X=FPQ(N,L,NZ)
  ENDIF
  IF(RP1BYs1(L).GT.ZEROSs1)THEN
    FP1BX=AMAX1(FPP(N,L,NZ),RUPP1Bs1(N,L,NZ)/RP1BYs1(L))
  ELSE
    FP1BX=FPQ(N,L,NZ)
  ENDIF

  TFPO4X=TFPO4X+FPO4X
  TFP14X=TFP14X+FP14X

  TFPOBX=TFPOBX+FPOBX
  TFP1BX=TFP1BX+FP1BX

  IF(FPUP.GT.ZERO2s1)THEN
    !
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF H2PO4,HPO4
    !     FROM SOIL TO ROOT
    !
    !     POSGL=PO4 diffusivity
    !     TORT=soil tortuosity
    !     PATH=path length of water and nutrient uptake
    !     RRADL=root radius
    !     DIFFL=PO4 diffusion per plant
    !
    POSGX=POSGLs1(L)*TORTs1(NPH,L)
    PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*POSGX))
    DIFFL=POSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))

    call UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,UPWTRP)

    call UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,UPWTRP)

  ELSE
    RUPP2Ps1(N,L,NZ)=0.0_r8
    RUPH2Ps1(N,L,NZ)=0.0_r8
    RUOH2Ps1(N,L,NZ)=0.0_r8
    RUCH2Ps1(N,L,NZ)=0.0_r8
    RUPP2Bs1(N,L,NZ)=0.0_r8
    RUPH2Bs1(N,L,NZ)=0.0_r8
    RUOH2Bs1(N,L,NZ)=0.0_r8
    RUCH2Bs1(N,L,NZ)=0.0_r8
    RUPP1Ps1(N,L,NZ)=0.0_r8
    RUPH1Ps1(N,L,NZ)=0.0_r8
    RUOH1Ps1(N,L,NZ)=0.0_r8
    RUCH1Ps1(N,L,NZ)=0.0_r8
    RUPP1Bs1(N,L,NZ)=0.0_r8
    RUPH1Bs1(N,L,NZ)=0.0_r8
    RUOH1Bs1(N,L,NZ)=0.0_r8
    RUCH1Bs1(N,L,NZ)=0.0_r8
  ENDIF
  end associate
  end subroutine UptakeMineralPhosporhus
!------------------------------------------------------------------------

  subroutine UptakeNO3(N,L,NZ,FNO3X,FNOBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: FNO3X,FNOBX,PATH(2,JZ1),RRADL(2,JZ1),RTARR(2,JZ1)
  real(r8), intent(in):: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFFL
  real(r8) :: DIFNO3,DIFNOB
  real(r8) :: PATHL
  real(r8) :: RMFNO3,RTKNO3,RTKNOP,RMFNOB,RTKNOB,RTKNPB
  real(r8) :: UPMX,UPMXP,X,Y
  real(r8) :: ZOSGX,ZNO3M,ZNO3X,ZNOBM,ZNOBX
! begin_execution
  associate(                              &
    VLNOBs1   =>  plt_soilchem%VLNOBs1  , &
    CNO3Bs1   =>  plt_soilchem%CNO3Bs1  , &
    ZNO3Ss1   =>  plt_soilchem%ZNO3Ss1  , &
    ZOSGLs1   =>  plt_soilchem%ZOSGLs1  , &
    CNO3Ss1   =>  plt_soilchem%CNO3Ss1  , &
    VOLWs1    =>  plt_soilchem%VOLWs1   , &
    VLNO3s1   =>  plt_soilchem%VLNO3s1  , &
    ZNO3Bs1   =>  plt_soilchem%ZNO3Bs1    &
  )
!
! PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NO3
! FROM SOIL TO ROOT
!
! ZOSGL=NO3 diffusivity
! TORT=soil tortuosity
! RRADL=root radius
! PATH=path length of water and nutrient uptake
! DIFFL=NO3 diffusion per plant
!
  ZOSGX=ZOSGLs1(L)*TORTs1(NPH,L)
  PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*ZOSGX))
  DIFFL=ZOSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))
  !
  ! NO3 UPTAKE IN NON-BAND SOIL ZONE
  !
  !  VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3S=NO3 concentration in non-band
  !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFNO3=soil-root convective NO3 flux per plant in non-band
  !     DIFNO3=soil-root NO3 diffusion per plant in non-band
  !
  IF(VLNO3s1(L).GT.ZEROs1.AND.CNO3Ss1(L) &
    .GT.UPMNZOs1(N,NZ))THEN
    RMFNO3=UPWTRP*CNO3Ss1(L)*VLNO3s1(L)
    DIFNO3=DIFFL*VLNO3s1(L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max NO3 uptake in non-band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXZOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNO3s1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNO3=soil-root convective N03 flux per plant in non-band
    !     DIFNO3=soil-root N03 diffusion per plant in non-band
    !     CNO3S=NO3 concentration in non-band
    !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
    !     RTKNO3,RTKNOP=NO3 uptake per plant in non-band lmtd,unlmtd by O2
    !     ZNO3M,ZNO3X=minimum,maximum NO3 available for uptake in non-band
    !     FNO3X=fraction of total NH4 uptake in non-band by root,myco populn
    !     RUNNOP,RUPNO3=NO3 uptake in non-band unlimited,limited by NO3
    !     RUONO3=NO3 uptake in non-band unlimited by O2
    !     RUCNO3=NO3 uptake in non-band unlimited by nonstructural C
!
    X=(DIFNO3+RMFNO3)*CNO3Ss1(L)
    Y=DIFNO3*UPMNZOs1(N,NZ)
    B=-UPMX-DIFNO3*UPKMZOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNO3=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNO3*UPKMZOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNOP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNO3M=UPMNZOs1(N,NZ)*VOLWs1(L)*VLNO3s1(L)
    ZNO3X=AMAX1(0.0,FNO3X*(ZNO3Ss1(L)-ZNO3M))
    RUNNOPs1(N,L,NZ)=AMAX1(0.0,RTKNO3*PPs1(NZ))
    RUPNO3s1(N,L,NZ)=AMIN1(ZNO3X,RUNNOPs1(N,L,NZ))
    RUONO3s1(N,L,NZ)=AMIN1(ZNO3X,AMAX1(0.0 &
      ,RTKNOP*PPs1(NZ)))
    RUCNO3s1(N,L,NZ)=RUPNO3s1(N,L,NZ)/FCUP
  ELSE
    RUNNOPs1(N,L,NZ)=0.0_r8
    RUPNO3s1(N,L,NZ)=0.0_r8
    RUONO3s1(N,L,NZ)=0.0_r8
    RUCNO3s1(N,L,NZ)=0.0_r8
  ENDIF
  !
  !     NO3 UPTAKE IN BAND SOIL ZONE
  !
  !     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3B=NO3 concentration in band
  !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFNOB=soil-root convective NO3 flux per plant in band
  !     DIFNOB=soil-root NO3 diffusion per plant in band
  !

  IF(VLNOBs1(L).GT.ZEROs1.AND.CNO3Bs1(L) &
    .GT.UPMNZOs1(N,NZ))THEN
    RMFNOB=UPWTRP*CNO3Bs1(L)*VLNOBs1(L)
    DIFNOB=DIFFL*VLNOBs1(L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum NO3 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXZOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNOBs1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNOB=soil-root convective NO3 flux per plant in band
    !     DIFNOB=soil-root NO3 diffusion per plant in band
    !     CNO3B=NH4 concentration in band
    !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
    !     RTKNOB,RTKNPB=NO3 uptake per plant in band lmtd,unlmtd by O2
    !     ZNOBM,ZNOBX=minimum,maximum NO3 available for uptake in band
    !     FNOBX=fraction of total NO3 uptake in band by root,myco populn
    !     RUNNXP,RUPNOB=NO3 uptake in band unlimited,limited by NH4
    !     RUONOB=NO3 uptake in band unlimited by O2
    !     RUCNOB=NO3 uptake in band unlimited by nonstructural C
    !
    X=(DIFNOB+RMFNOB)*CNO3Bs1(L)
    Y=DIFNOB*UPMNZOs1(N,NZ)
    B=-UPMX-DIFNOB*UPKMZOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNOB=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNOB*UPKMZOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNOBM=UPMNZOs1(N,NZ)*VOLWs1(L)*VLNOBs1(L)
    ZNOBX=AMAX1(0.0,FNOBX*(ZNO3Bs1(L)-ZNOBM))
    RUNNXPs1(N,L,NZ)=AMAX1(0.0,RTKNOB*PPs1(NZ))
    RUPNOBs1(N,L,NZ)=AMIN1(ZNOBX,RUNNXPs1(N,L,NZ))
    RUONOBs1(N,L,NZ)=AMIN1(ZNOBX,AMAX1(0.0,RTKNPB*PPs1(NZ)))
    RUCNOBs1(N,L,NZ)=RUPNOBs1(N,L,NZ)/FCUP
  ELSE
    RUNNXPs1(N,L,NZ)=0.0_r8
    RUPNOBs1(N,L,NZ)=0.0_r8
    RUONOBs1(N,L,NZ)=0.0_r8
    RUCNOBs1(N,L,NZ)=0.0_r8
  ENDIF
  end associate
  end subroutine UptakeNO3
!------------------------------------------------------------------------

  subroutine UptakeNH4(N,L,NZ,FNH4X,FNHBX,PATH,RRADL,RTARR,&
    FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: FNH4X,FNHBX,PATH(2,JZ1),RRADL(2,JZ1)
  real(r8), intent(in) :: RTARR(2,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFFL
  real(r8) :: DIFNH4,DIFNHB
  real(r8) :: PATHL
  real(r8) :: RMFNH4,RTKNH4,RTKNHP,RMFNHB,RTKNHB,RTKNBP
  real(r8) :: UPMX,UPMXP,X,Y
  real(r8) :: ZNHBX,ZNSGX,ZNH4M,ZNH4X,ZNHBM
! begin_execution
  associate(                              &
    ZNSGLs1   =>  plt_soilchem%ZNSGLs1  , &
    VOLWs1    =>  plt_soilchem%VOLWs1   , &
    VLNHBs1   =>  plt_soilchem%VLNHBs1  , &
    ZNH4Bs1   =>  plt_soilchem%ZNH4Bs1  , &
    CNH4Bs1   =>  plt_soilchem%CNH4Bs1  , &
    VLNH4s1   =>  plt_soilchem%VLNH4s1  , &
    ZNH4Ss1   =>  plt_soilchem%ZNH4Ss1  , &
    CNH4Ss1   =>  plt_soilchem%CNH4Ss1    &
  )
! ZNSGL=NH4 diffusivity
! TORT=soil tortuosity
! PATH=path length of water and nutrient uptake
! RRADL=root radius
! DIFFL=NH4 diffusion per plant
!
  ZNSGX=ZNSGLs1(L)*TORTs1(NPH,L)
  PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*ZNSGX))
  DIFFL=ZNSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))
!
! NH4 UPTAKE IN NON-BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4S=NH4 concentration in non-band
! UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
! UPWTRP=root water uptake per plant
! RMFNH4=soil-root convective NH4 flux per plant in non-band
! DIFNH4=soil-root NH4 diffusion per plant in non-band
!
  IF(VLNH4s1(L).GT.ZEROs1.AND.CNH4Ss1(L).GT.UPMNZHs1(N,NZ))THEN
    RMFNH4=UPWTRP*CNH4Ss1(L)*VLNH4s1(L)
    DIFNH4=DIFFL*VLNH4s1(L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=max NH4 uptake in non-band unlimited,limited by O2
!   RTARP=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   TFN4=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   WFR=constraint by O2 consumption on all biological processes
!
    UPMXP=UPMXZHs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNH4s1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNH4=soil-root convective NH4 flux per plant in non-band
!   DIFNH4=soil-root NH4 diffusion per plant in non-band
!   CNH4S=NH4 concentration in non-band
!   UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
!   RTKNH4,RTKNHP=NH4 uptake per plant in non-band lmtd,unlmtd by O2
!   ZNH4M,ZNH4X=minimum,maximum NH4 available for uptake in non-band
!   FNH4X=fraction of total NH4 uptake in non-band by root,myco populn
!   RUNNHP,RUPNH4=NH4 uptake in non-band unlimited,limited by NH4
!   RUONH4=NH4 uptake in non-band unlimited by O2
!   RUCNH4=NH4 uptake in non-band unlimited by nonstructural C
!
    X=(DIFNH4+RMFNH4)*CNH4Ss1(L)
    Y=DIFNH4*UPMNZHs1(N,NZ)
    B=-UPMX-DIFNH4*UPKMZHs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNH4=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNH4*UPKMZHs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNHP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNH4M=UPMNZHs1(N,NZ)*VOLWs1(L)*VLNH4s1(L)
    ZNH4X=AMAX1(0.0,FNH4X*(ZNH4Ss1(L)-ZNH4M))
    RUNNHPs1(N,L,NZ)=AMAX1(0.0,RTKNH4*PPs1(NZ))
    RUPNH4s1(N,L,NZ)=AMIN1(ZNH4X,RUNNHPs1(N,L,NZ))
    RUONH4s1(N,L,NZ)=AMIN1(ZNH4X,AMAX1(0.0 &
      ,RTKNHP*PPs1(NZ)))
    RUCNH4s1(N,L,NZ)=RUPNH4s1(N,L,NZ)/FCUP
  ELSE
    RUNNHPs1(N,L,NZ)=0.0_r8
    RUPNH4s1(N,L,NZ)=0.0_r8
    RUONH4s1(N,L,NZ)=0.0_r8
    RUCNH4s1(N,L,NZ)=0.0_r8
  ENDIF
!
! NH4 UPTAKE IN BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4B=NH4 concentration in band
! UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
! UPWTRP=root water uptake per plant
! RMFNHB=soil-root convective NH4 flux per plant in band
! DIFNHB=soil-root NH4 diffusion per plant in band
!

  IF(VLNHBs1(L).GT.ZEROs1.AND.CNH4Bs1(L) &
    .GT.UPMNZHs1(N,NZ))THEN
    RMFNHB=UPWTRP*CNH4Bs1(L)*VLNHBs1(L)
    DIFNHB=DIFFL*VLNHBs1(L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=maximum NH4 uptake in band unlimited,limited by O2
!   RTARP=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   TFN4=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   WFR=constraint by O2 consumption on all biological processes
!
    UPMXP=UPMXZHs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNHBs1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNHB=soil-root convective NH4 flux per plant in band
!   DIFNHB=soil-root NH4 diffusion per plant in band
!   CNH4B=NH4 concentration in band
!   UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
!   RTKNHB,RTKNBP=NH4 uptake per plant in band lmtd,unlmtd by O2
!   ZNHBM,ZNHBX=minimum,maximum NH4 available for uptake in band
!   FNHBX=fraction of total NH4 uptake in band by root,myco populn
!   RUNNBP,RUPNHB=NH4 uptake in band unlimited,limited by NH4
!   RUONHB=NH4 uptake in band unlimited by O2
!   RUCNHB=NH4 uptake in band unlimited by nonstructural C
!
    X=(DIFNHB+RMFNHB)*CNH4Bs1(L)
    Y=DIFNHB*UPMNZHs1(N,NZ)
    B=-UPMX-DIFNHB*UPKMZHs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNHB=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNHB*UPKMZHs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNBP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNHBM=UPMNZHs1(N,NZ)*VOLWs1(L)*VLNHBs1(L)
    ZNHBX=AMAX1(0.0,FNHBX*(ZNH4Bs1(L)-ZNHBM))
    RUNNBPs1(N,L,NZ)=AMAX1(0.0,RTKNHB*PPs1(NZ))
    RUPNHBs1(N,L,NZ)=AMIN1(ZNHBX,RUNNBPs1(N,L,NZ))
    RUONHBs1(N,L,NZ)=AMIN1(ZNHBX,AMAX1(0.0 &
      ,RTKNBP*PPs1(NZ)))
    RUCNHBs1(N,L,NZ)=RUPNHBs1(N,L,NZ)/FCUP
  ELSE
    RUNNBPs1(N,L,NZ)=0.0_r8
    RUPNHBs1(N,L,NZ)=0.0_r8
    RUONHBs1(N,L,NZ)=0.0_r8
    RUCNHBs1(N,L,NZ)=0.0_r8
  ENDIF
  end associate
  end subroutine UptakeNH4
!------------------------------------------------------------------------

  subroutine UptakeMineralNitrogen(N,L,NZ,PATH,RRADL,FPQ,&
    FPP,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1)
  real(r8), intent(in) :: FPP(2,JZ1,JP1),RTARR(2,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: FNO3X,FNOBX,FNH4X,FNHBX
  real(r8) :: TFNH4X,TFNO3X,TFNHBX,TFNOBX

!     begin_execution

  TFNH4X=0.0_r8
  TFNHBX=0.0_r8
  TFNO3X=0.0_r8
  TFNOBX=0.0_r8


  IF(RNH4Ys1(L).GT.ZEROSs1)THEN
    FNH4X=AMAX1(FPP(N,L,NZ),RUNNHPs1(N,L,NZ)/RNH4Ys1(L))
  ELSE
    FNH4X=FPQ(N,L,NZ)
  ENDIF
  IF(RNHBYs1(L).GT.ZEROSs1)THEN
    FNHBX=AMAX1(FPP(N,L,NZ),RUNNBPs1(N,L,NZ)/RNHBYs1(L))
  ELSE
    FNHBX=FPQ(N,L,NZ)
  ENDIF

  IF(RNO3Ys1(L).GT.ZEROSs1)THEN
    FNO3X=AMAX1(FPP(N,L,NZ),RUNNOPs1(N,L,NZ)/RNO3Ys1(L))
  ELSE
    FNO3X=FPQ(N,L,NZ)
  ENDIF
  IF(RN3BYs1(L).GT.ZEROSs1)THEN
    FNOBX=AMAX1(FPP(N,L,NZ),RUNNXPs1(N,L,NZ)/RN3BYs1(L))
  ELSE
    FNOBX=FPQ(N,L,NZ)
  ENDIF

  TFNH4X=TFNH4X+FNH4X
  TFNO3X=TFNO3X+FNO3X
  TFNHBX=TFNHBX+FNHBX
  TFNOBX=TFNOBX+FNOBX

  IF(FZUP.GT.ZERO2s1)THEN
!
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
    !     FROM SOIL TO ROOT
    !
    call UptakeNH4(N,L,NZ,FNH4X,FNHBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

    call UptakeNO3(N,L,NZ,FNO3X,FNOBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  ELSE
    RUNNHPs1(N,L,NZ)=0.0_r8
    RUPNH4s1(N,L,NZ)=0.0_r8
    RUONH4s1(N,L,NZ)=0.0_r8
    RUCNH4s1(N,L,NZ)=0.0_r8
    RUNNBPs1(N,L,NZ)=0.0_r8
    RUPNHBs1(N,L,NZ)=0.0_r8
    RUONHBs1(N,L,NZ)=0.0_r8
    RUCNHBs1(N,L,NZ)=0.0_r8
    RUNNOPs1(N,L,NZ)=0.0_r8
    RUPNO3s1(N,L,NZ)=0.0_r8
    RUONO3s1(N,L,NZ)=0.0_r8
    RUCNO3s1(N,L,NZ)=0.0_r8
    RUNNXPs1(N,L,NZ)=0.0_r8
    RUPNOBs1(N,L,NZ)=0.0_r8
    RUONOBs1(N,L,NZ)=0.0_r8
    RUCNOBs1(N,L,NZ)=0.0_r8
  ENDIF
  end subroutine UptakeMineralNitrogen
!------------------------------------------------------------------------

  subroutine RootExudates(N,L,NZ)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ

  real(r8) :: CPOOLX,CPOOLT
  real(r8) :: PPOOLX,ZPOOLX
  real(r8) :: VOLWK,VOLWT
  real(r8) :: XFRC,XFRN,XFRP
  integer :: K
  !     begin_execution
  associate(                           &
    CPOOLRs1=>  plt_biom%CPOOLRs1    , &
    ZPOOLRs1=>  plt_biom%ZPOOLRs1    , &
    PPOOLRs1=>  plt_biom%PPOOLRs1    , &
    OQNs1   =>  plt_soilchem%OQNs1   , &
    OQPs1   =>  plt_soilchem%OQPs1   , &
    OQCs1   =>  plt_soilchem%OQCs1     &
  )
  !
  !     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
  !     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
  !
  !     VOLWMM=soil micropore water volume
  !     FOSRH=fraction of total SOC in each substrate K from nitro.f
  !     RTVLW=root aqueous volume
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,myco
  !     XFRC,XFRN,XFRP=nonstructural C,N,P exchg at root-soil DOC equilibrium
  !     OQC=soil DOC
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation
  !     TLEC,TSHC=total fluxes x blr for calculating canopy air temperature,
  !     vapor pressure in watsub.f
  !      EFLXC,SFLXC=canopylatent,sensible heat fluxes
  !      RA=canopy boundary layer resistance
  !     OSTR=O2 stress indicator
  !
  DO 195 K=0,jcplx11
    VOLWK=VOLWMs1(NPH,L)*FOSRHs1(K,L)
    IF(VOLWK.GT.ZEROS2s1 &
      .AND.RTVLWs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
      VOLWT=VOLWK+RTVLWs1(N,L,NZ)
      CPOOLX=AMIN1(1.25E+03*RTVLWs1(N,L,NZ),CPOOLRs1(N,L,NZ))
      XFRC=(OQCs1(K,L)*RTVLWs1(N,L,NZ)-CPOOLX*VOLWK)/VOLWT
      RDFOMCs1(N,K,L,NZ)=FEXUC*XFRC
      IF(OQCs1(K,L).GT.ZEROSs1 &
        .AND.CPOOLRs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
        CPOOLT=OQCs1(K,L)+CPOOLRs1(N,L,NZ)
        ZPOOLX=0.1*ZPOOLRs1(N,L,NZ)
        PPOOLX=0.1*PPOOLRs1(N,L,NZ)
        XFRN=(OQNs1(K,L)*CPOOLRs1(N,L,NZ)-ZPOOLX*OQCs1(K,L))/CPOOLT
        XFRP=(OQPs1(K,L)*CPOOLRs1(N,L,NZ)-PPOOLX*OQCs1(K,L))/CPOOLT
        RDFOMNs1(N,K,L,NZ)=FEXUN*XFRN
        RDFOMPs1(N,K,L,NZ)=FEXUP*XFRP
      ELSE
        RDFOMNs1(N,K,L,NZ)=0.0_r8
        RDFOMPs1(N,K,L,NZ)=0.0_r8
      ENDIF
    ELSE
      RDFOMCs1(N,K,L,NZ)=0.0_r8
      RDFOMNs1(N,K,L,NZ)=0.0_r8
      RDFOMPs1(N,K,L,NZ)=0.0_r8
    ENDIF

195   CONTINUE
  end associate
  end subroutine RootExudates
!------------------------------------------------------------------------

  subroutine SumupNutrientUptake(N,L,NZ)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ

  integer :: K
  !     begin_execution
  !
  !     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL
  !
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
  !     XOQCS,XOQNZ,XOQPS=accumulated change in DOC,DON,DOP from nitro.f
  !     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
  !     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
  !     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
  !
  DO 295 K=0,jcplx11
    UPOMCs1(NZ)=UPOMCs1(NZ)+RDFOMCs1(N,K,L,NZ)
    UPOMNs1(NZ)=UPOMNs1(NZ)+RDFOMNs1(N,K,L,NZ)
    UPOMPs1(NZ)=UPOMPs1(NZ)+RDFOMPs1(N,K,L,NZ)
    XOQCSs1(K,L)=XOQCSs1(K,L)-RDFOMCs1(N,K,L,NZ)
    XOQNSs1(K,L)=XOQNSs1(K,L)-RDFOMNs1(N,K,L,NZ)
    XOQPSs1(K,L)=XOQPSs1(K,L)-RDFOMPs1(N,K,L,NZ)
295   CONTINUE
  UPNH4s1(NZ)=UPNH4s1(NZ)+RUPNH4s1(N,L,NZ)+RUPNHBs1(N,L,NZ)
  UPNO3s1(NZ)=UPNO3s1(NZ)+RUPNO3s1(N,L,NZ)+RUPNOBs1(N,L,NZ)
  UPH2Ps1(NZ)=UPH2Ps1(NZ)+RUPH2Ps1(N,L,NZ)+RUPH2Bs1(N,L,NZ)
  UPH1Ps1(NZ)=UPH1Ps1(NZ)+RUPH1Ps1(N,L,NZ)+RUPH1Bs1(N,L,NZ)
  end subroutine SumupNutrientUptake
!------------------------------------------------------------------------

  subroutine GetUptakeCapcity(N,L,NZ,FPQ,FPP,FCUP,FZUP,FPUP,&
    FWSRT,UPWTRP,UPWTRH,FOXYX)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  REAL(R8), INTENT(IN):: FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(out):: FCUP,FZUP,FPUP,FWSRT,UPWTRP,UPWTRH,FOXYX

  associate(                          &
    CWSRTLs1  => plt_biom%CWSRTLs1  , &
    WSRTLs1   => plt_biom%WSRTLs1   , &
    CPPOLRs1  => plt_biom%CPPOLRs1  , &
    CZPOLRs1  => plt_biom%CZPOLRs1  , &
    CPOOLRs1  => plt_biom%CPOOLRs1  , &
    CCPOLRs1  => plt_biom%CCPOLRs1  , &
    WTRTLs1   => plt_biom%WTRTLs1     &
  )
  !
  !     UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL
  !     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
  !     PARAMETERS ARE DEFINED
  !
  !     CWSRTL,CWSRT=current,maximum protein concentration
  !     WSRTL,WTRTL=protein content,mass
  !     FWSRT=protein concentration relative to 5%
  !
  IF(WTRTLs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
    CWSRTLs1(N,L,NZ)=AMIN1(CWSRTs1(NZ),WSRTLs1(N,L,NZ)/WTRTLs1(N,L,NZ))
    FWSRT=CWSRTLs1(N,L,NZ)/0.05
  ELSE
    CWSRTLs1(N,L,NZ)=CWSRTs1(NZ)
    FWSRT=1.0
  ENDIF
  !
  !     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
  !
  !     RCO2N=total respiration from CPOOLR
  !     FCUP=limitation to active uptake respiration from CPOOLR
  !     CPOOLR=nonstructural C content
  !
  IF(RCO2Ns1(N,L,NZ).GT.ZEROPs1(NZ))THEN
    FCUP=AMAX1(0.0,AMIN1(1.0,0.25*safe_adb(CPOOLRs1(N,L,NZ) &
      ,RCO2Ns1(N,L,NZ))))
  ELSE
    FCUP=0.0_r8
  ENDIF
  !
  !     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
  !
  !     FZUP,FPUP=limitn to active uptake respiration from CZPOLR,CPPOLR
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake
  !     UPWTRH=water uptake at time step for gas flux calculations
  !
  IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
    FZUP=AMIN1(safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)/ZCKI) &
      ,safe_adb(CPPOLRs1(N,L,NZ),CPPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)/ZPKI))
    FPUP=AMIN1(safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)/PCKI) &
      ,safe_adb(CZPOLRs1(N,L,NZ),CZPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)/PZKI))
  ELSE
    FZUP=0.0_r8
    FPUP=0.0_r8
  ENDIF
  !NN=0
  UPWTRP=AMAX1(0.0,-UPWTRs1(N,L,NZ)/PPs1(NZ))
  UPWTRH=UPWTRP*XNPG
  !
  !     FACTORS CONSTRAINING O2 AND NUTRIENT UPTAKE AMONG
  !     COMPETING ROOT,MYCORRHIZAL AND MICROBIAL POPULATIONS
  !     IN BAND AND NON-BAND SOIL ZONES FROM DEMAND CALCULATED
  !     IN PREVIOUS HOUR
  !
  !     ROXYY=O2 demand by all microbial,root,myco populations
  !     ROXYP=O2 demand by each root,myco population
  !     FOXYX=fraction of ROXYY by each root,myco population
  !     RNH4Y=NH4 demand in non-band by all microbial,root,myco populations
  !     RUNNHP=NH4 demand in non-band by each root,myco population
  !     FNH4X=fraction of RNH4Y by each root,myco populn
  !     RNHBY=NH4 demand in band by all microbial,root,myco populations
  !     RUNNBP=NH4 demand in band by each root,myco population
  !     FNHBX=fraction of RNHBY by each root,myco populn
  !     RNO3Y=NO3 demand in non-band by all microbial,root,myco populations
  !     RUNNOP=NO3 demand in non-band by each root,myco population
  !     FNO3X=fraction of RNO3Y by each root,myco populn
  !     RN3BY=NO3 demand in band by all microbial,root,myco populations
  !     RUNNXB=NO3 demand in band by each root,myco population
  !     FNOBX=fraction of RN3BY by each root,myco populn
  !     RPO4Y=H2PO4 demand in non-band by all microbial,root,myco populations
  !     RUPP2P=H2PO4 demand in non-band by each root,myco population
  !     FPO4X=fraction of RPO4Y by each root,myco populn
  !     RPOBY=H2PO4 demand in band by all microbial,root,myco populations
  !     RUPP2B=H2PO4 demand in band by each root,myco population
  !     FPOBX=fraction of RPOBY by each root,myco populn
  !     RP14Y=HPO4 demand in non-band by all microbial,root,myco populations
  !     RUPP1P=HPO4 demand in non-band by each root,myco population
  !     FP14X=fraction of RP14Y by each root,myco populn
  !     RP1BY=HPO4 demand in band by all microbial,root,myco populations
  !     RUPP1B=HPO4 demand in band by each root,myco population
  !     FP1BX=fraction of RP1BY by each root,myco populn
  !     FPP=minimum uptake fraction
  !     FPQ=PFT fraction of biome root mass
  !
  IF(ROXYYs1(L).GT.ZEROSs1)THEN
    FOXYX=AMAX1(FPP(N,L,NZ),ROXYPs1(N,L,NZ)/ROXYYs1(L))
  ELSE
    FOXYX=FPQ(N,L,NZ)
  ENDIF
  end associate
  end subroutine GetUptakeCapcity
!------------------------------------------------------------------------

  subroutine RootSoilGasExchange(N,L,NZ,RRADL,FPQ,FRTDPX,RTARR,&
    UPWTRH,FOXYX,RUPOXT)

  implicit none
  integer , intent(in) :: N,L,NZ
  real(r8), intent(in) :: RRADL(2,JZ1),FPQ(2,JZ1,JP1),FRTDPX(JZ1,JP1)
  real(r8), intent(in) :: RTARR(2,JZ1),UPWTRH,FOXYX
  real(r8), intent(out):: RUPOXT
  integer :: M,MX
  real(r8) :: B,C
  real(r8) :: CO2A1,CO2P1,CO2G1,CO2S1,CH4A1
  real(r8) :: CH4P1,CH4S1,CCH4S1,CCH4P1
  real(r8) :: CN2OS1,CN2OP1,CNH3S1,CNH3B1,CNH3P1,CH2GS1,CH2GP1
  real(r8) :: CGSGL1,CHSGL1,CLSGL1
  real(r8) :: CQSGL1,CH4G1,CCO2S1,COXYS1,CCO2A1,COXYA1,CCH4A1
  real(r8) :: CZ2OA1,CNH3A1,CH2GA1,CCO2P1,COXYP1,COXYR,CO2PX,CH4PX
  real(r8) :: DIFOP,DIFCL,DIFZL,DIFNL,DIFNB,DIFHL,DIFOX,DFGSP
  real(r8) :: DFCOA,DFOXA,DFCHA,DFN2A,DFNHA,DFHGA,DFGP,DIFOL
  real(r8) :: H2GA1,H2GP1,H2GS1,HGSGL1,HLSGL1,H2GG1,H2GPX
  real(r8) :: OXYA1,OXYP1,OXYG1,OXYS1,OGSGL1
  real(r8) :: OLSGL1,OLSGLP,OXYPX,RTVLWA
  real(r8) :: RTVLWB,ROXYFX,RCO2FX,ROXYLX,ROXDFQ,RCHDFQ,RN2DFQ
  real(r8) :: RNHDFQ,RHGDFQ,ROXDF1,RCHDF1,RN2DF1,RNHDF1,RHGDF1
  real(r8) :: RTCR1,RTCR2,RTCRA,RTARRX,RCO2PX,RRADS,RMFCOS,RMFOXS
  real(r8) :: RMFCHS,RMFN2S,RMFN3S,RMFN3B,RMFHGS,RUPOXR,RDFOXS
  real(r8) :: RDFOXP,RUPOSX,RUPOPX,RDFCOS,RDXCOS,RCO2SX,RDFCHS
  real(r8) :: RDXCHS,RUPCSX,RDFN2S,RDXN2S,RUPZSX,RDFN3S,RDXNHS
  real(r8) :: RUPNSX,RDFN3B,RDXNHB,RUPNBX,RDFHGS,RDXHGS,RUPHGX
  real(r8) :: RCODFQ,RUPOST,RNBDFQ,RUPNTX,RCODF1,RCOFL1,ROXFL1
  real(r8) :: RCHFL1,RN2FL1,RNHFL1,RHGFL1,THETW1,THETM
  real(r8) :: UPMXP
  real(r8) :: VOLWCA,VOLWOA
  real(r8) :: VOLWC4,VOLWZA,VOLWNA,VOLWH2,VOLWMO,VOLWMM,VOLPMM
  real(r8) :: VOLWSP,VOLWMA,VOLWMB,VOLWSA,VOLWSB,VOLWCO,VOLWOX
  real(r8) :: VOLWCH,VOLWN2,VOLWNH,VOLWNB,VOLWHG,VOLPNH,VOLPNB
  real(r8) :: X
  real(r8) :: Z2OA1,Z2OP1,Z2OS1,ZH3A1,ZH3P1,ZH3S1
  real(r8) :: ZH3B1,Z2SGL1,ZHSGL1,ZVSGL1,ZNSGL1,Z2OG1,ZH3G1
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB,Z2OPX,ZH3PX

!     begin_execution
  associate(                           &
    WTRTSs1  =>  plt_biom%WTRTSs1    , &
    TFNDs1   =>  plt_soilchem%TFNDs1 , &
    CO2Gs1   =>  plt_soilchem%CO2Gs1 , &
    CO2Ss1   =>  plt_soilchem%CO2Ss1 , &
    ZNH3Ss1  =>  plt_soilchem%ZNH3Ss1, &
    SN2OLs1  =>  plt_soilchem%SN2OLs1, &
    H2GSs1   =>  plt_soilchem%H2GSs1 , &
    CZ2OSs1  =>  plt_soilchem%CZ2OSs1, &
    CNH3Gs1  =>  plt_soilchem%CNH3Gs1, &
    CH2GGs1  =>  plt_soilchem%CH2GGs1, &
    CH2GSs1  =>  plt_soilchem%CH2GSs1, &
    CZ2OGs1  =>  plt_soilchem%CZ2OGs1, &
    VLNH4s1  =>  plt_soilchem%VLNH4s1, &
    HGSGLs1  =>  plt_soilchem%HGSGLs1, &
    SCH4Ls1  =>  plt_soilchem%SCH4Ls1, &
    CCH4Gs1  =>  plt_soilchem%CCH4Gs1, &
    CLSGLs1  =>  plt_soilchem%CLSGLs1, &
    THETYs1  =>  plt_soilchem%THETYs1, &
    OLSGLs1  =>  plt_soilchem%OLSGLs1, &
    SOXYLs1  =>  plt_soilchem%SOXYLs1, &
    SNH3Ls1  =>  plt_soilchem%SNH3Ls1, &
    SH2GLs1  =>  plt_soilchem%SH2GLs1, &
    VOLYs1   =>  plt_soilchem%VOLYs1 , &
    ZNSGLs1  =>  plt_soilchem%ZNSGLs1, &
    Z2SGLs1  =>  plt_soilchem%Z2SGLs1, &
    ZHSGLs1  =>  plt_soilchem%ZHSGLs1, &
    ZVSGLs1  =>  plt_soilchem%ZVSGLs1, &
    SCO2Ls1  =>  plt_soilchem%SCO2Ls1, &
    HLSGLs1  =>  plt_soilchem%HLSGLs1, &
    CQSGLs1  =>  plt_soilchem%CQSGLs1, &
    VLNHBs1  =>  plt_soilchem%VLNHBs1, &
    ZNH3Bs1  =>  plt_soilchem%ZNH3Bs1, &
    CGSGLs1  =>  plt_soilchem%CGSGLs1, &
    CNH3Bs1  =>  plt_soilchem%CNH3Bs1, &
    CHSGLs1  =>  plt_soilchem%CHSGLs1, &
    OGSGLs1  =>  plt_soilchem%OGSGLs1, &
    Z2OSs1   =>  plt_soilchem%Z2OSs1 , &
    CNH3Ss1  =>  plt_soilchem%CNH3Ss1, &
    THETPMs1 =>  plt_soilchem%THETPMs1,&
    OXYGs1   =>  plt_soilchem%OXYGs1 , &
    DFGSs1   =>  plt_soilchem%DFGSs1 , &
    CCH4Ss1  =>  plt_soilchem%CCH4Ss1, &
    OXYSs1   =>  plt_soilchem%OXYSs1 , &
    CH4Ss1   =>  plt_soilchem%CH4Ss1 , &
    IDAYs1   =>  plt_pheno%IDAYs1    , &
    NGs1     =>   plt_morph%NGs1     , &
    NB1s1    =>  plt_morph%NB1s1       &
  )
  IF(RCO2Ms1(N,L,NZ).GT.ZEROPs1(NZ) &
    .AND.RTVLWs1(N,L,NZ).GT.ZEROPs1(NZ) &
    .AND.FOXYX.GT.ZEROQs1(NZ))THEN
!
!     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
!     (CO2, O2, CH4, N2, N2O, NH3, H2)
!
!     CO2A1,CO2P1,CO2G1,CO2S1=gaseous,aqueous CO2 in root,soil
!     OXYA1,OXYP1,OXYG1,OXYS1=gaseous,aqueous O2 in root,soil
!     CH4A1,CH4P1,CH4G1,CH4S1=gaseous,aqueous CH4 in root,soil
!     Z2OA1,Z2OP1,Z2OG1,Z2OS1=gaseous,aqueous N2O in root,soil
!     ZH3A1,ZH3P1,ZH3G1,ZH3S1=gaseous,aqueous NH3 in root,soil
!     H2GA1,H2GP1,H2GG1,H2GS1=gaseous,aqueous H2 in root,soil
!     CCH4S1,CCH4P1=aqueous CH4 concentration in soil,root
!     CN2OS1,CN2OP1=aqueous N2O concentration in soil,root
!     CNH3S1,CNH3B1,CNH3P1=aqueous NH3 concn in soil non-band,band,root
!     CH2GS1,CH2GP1=aqueous H2 concentration in soil,root
!     RTVLWA,RTVLWB=root aqueous volume in non-band,band
!     XNPG=time step of flux calculation
!     UPMXP=O2 demand per plant at time step of flux calculation
!     ROXYFX=net O2 gas flux at time step of flux calculation
!     RCO2FX=net CO2 gas flux at time step of flux calculation
!     ROXYLX=net O2 aqueous flux at time step of flux calculation
!
    CO2A1=AMAX1(ZEROPs1(NZ),CO2As1(N,L,NZ))
    CO2P1=AMAX1(ZEROPs1(NZ),CO2Ps1(N,L,NZ))
    CO2G1=AMAX1(ZEROPs1(NZ),CO2Gs1(L)*FPQ(N,L,NZ))
    CO2S1=AMAX1(ZEROPs1(NZ),CO2Ss1(L)*FPQ(N,L,NZ))
    OXYA1=AMAX1(ZEROPs1(NZ),OXYAs1(N,L,NZ))
    OXYP1=AMAX1(ZEROPs1(NZ),OXYPs1(N,L,NZ))
    OXYG1=AMAX1(ZEROPs1(NZ),OXYGs1(L)*FOXYX)
    OXYS1=OXYSs1(L)*FOXYX
    CH4A1=CH4As1(N,L,NZ)
    CH4P1=CH4Ps1(N,L,NZ)
    CH4S1=CH4Ss1(L)*FPQ(N,L,NZ)
    CCH4S1=CCH4Ss1(L)
    CCH4P1=AMAX1(0.0,CH4P1/RTVLWs1(N,L,NZ))
    Z2OA1=Z2OAs1(N,L,NZ)
    Z2OP1=Z2OPs1(N,L,NZ)
    Z2OS1=Z2OSs1(L)*FPQ(N,L,NZ)
    CN2OS1=CZ2OSs1(L)
    CN2OP1=AMAX1(0.0,Z2OP1/RTVLWs1(N,L,NZ))
    ZH3A1=ZH3As1(N,L,NZ)
    ZH3P1=ZH3Ps1(N,L,NZ)
    ZH3S1=ZNH3Ss1(L)*FPQ(N,L,NZ)
    ZH3B1=ZNH3Bs1(L)*FPQ(N,L,NZ)
    CNH3S1=CNH3Ss1(L)
    CNH3B1=CNH3Bs1(L)
    CNH3P1=AMAX1(0.0,ZH3P1/RTVLWs1(N,L,NZ))
    H2GA1=H2GAs1(N,L,NZ)
    H2GP1=H2GPs1(N,L,NZ)
    H2GS1=H2GSs1(L)*FPQ(N,L,NZ)
    CH2GS1=CH2GSs1(L)
    CH2GP1=AMAX1(0.0,H2GP1/RTVLWs1(N,L,NZ))
    RTVLWA=RTVLWs1(N,L,NZ)*VLNH4s1(L)
    RTVLWB=RTVLWs1(N,L,NZ)*VLNHBs1(L)
    UPMXP=ROXYPs1(N,L,NZ)*XNPG/PPs1(NZ)
    ROXYFX=ROXYFs1(L)*FOXYX*XNPG
    RCO2FX=RCO2Fs1(L)*FOXYX*XNPG
    ROXYLX=ROXYLs1(L)*FOXYX*XNPG
!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
!     *SGL1=diffusivity
!     PORTX=tortuosity effect of root porosity on diffusivity
!     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
!     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
!
    CGSGL1=CGSGLs1(L)*XNPG*PORTXs1(N,NZ)
    OGSGL1=OGSGLs1(L)*XNPG*PORTXs1(N,NZ)
    CHSGL1=CHSGLs1(L)*XNPG*PORTXs1(N,NZ)
    Z2SGL1=Z2SGLs1(L)*XNPG*PORTXs1(N,NZ)
    ZHSGL1=ZHSGLs1(L)*XNPG*PORTXs1(N,NZ)
    HGSGL1=HGSGLs1(L)*XNPG*PORTXs1(N,NZ)
    CLSGL1=CLSGLs1(L)*XNPG*FOXYX
    OLSGL1=OLSGLs1(L)*XNPG*FOXYX
    CQSGL1=CQSGLs1(L)*XNPG*FOXYX
    ZVSGL1=ZVSGLs1(L)*XNPG*FOXYX
    ZNSGL1=ZNSGLs1(L)*XNPG*FOXYX
    HLSGL1=HLSGLs1(L)*XNPG*FOXYX
    OLSGLP=OLSGLs1(L)*XNPG
    ROXDFQ=0.0_r8
    RCHDFQ=0.0_r8
    RN2DFQ=0.0_r8
    RNHDFQ=0.0_r8
    RHGDFQ=0.0_r8
    ROXDF1=0.0_r8
    RCHDF1=0.0_r8
    RN2DF1=0.0_r8
    RNHDF1=0.0_r8
    RHGDF1=0.0_r8
!
!     ROOT CONDUCTANCE TO GAS TRANSFER
!
!     WTRTS=total root,myco mass
!     FRTDPX=fraction of each soil layer with primary root
!     RTCR1,RTCR2,RTCRA=cross-sectional area/length of
!     primary,secondary,total root,myco system
!     RTN1,RTNL=number of root,myco primary,secondary axes
!     RRAD1,RRAD2=primary,secondary root radius
!     DPTHZ=depth of primary root from surface
!     RTLGA=average secondary root length
!
    IF(WTRTSs1(NZ).GT.ZEROPs1(NZ).AND.FRTDPX(L,NZ).GT.ZEROs1)THEN
      RTCR1=AMAX1(PPs1(NZ),RTN1s1(N,L,NZ)) &
        *PICON*RRAD1s1(N,L,NZ)**2/DPTHZs1(L)
      RTCR2=(RTNLs1(N,L,NZ)*PICON*RRAD2s1(N,L,NZ)**2 &
        /RTLGAs1(N,L,NZ))/FRTDPX(L,NZ)
      IF(RTCR2.GT.RTCR1)THEN
        RTCRA=RTCR1*RTCR2/(RTCR1+RTCR2)
      ELSE
        RTCRA=RTCR1
      ENDIF
    ELSE
      RTCRA=0.0_r8
    ENDIF
!
!     VARIABLES USED TO CALCULATE ROOT GAS TRANSFER
!     BETWEEN AQUEOUS AND GASEOUS PHASES
!
!     RTLGP=root,myco length per plant
!     IDAYs1(1,=emergence date
!     RTARR=root surface area/radius for uptake
!     RRADP=path length for radial diffusion within root
!     DIFOP=aqueous diffusivity of O2 within root
!     RTVLW=root,myco aqueous volume
!     S*L=solubility of gas in water from hour1.f:
!     CO2=CO2,OXY=O2,CH4=CH4,N2O=N2O,NH3=NH3,H2G=H2
!     DF*A=root-atmosphere gas conductance
!     DFGP=rate const for equilibrn of gas concn in gaseous-aqueous phases
!     RCO2PX=root CO2 gas flux at time step for gas flux calculations
!     RCO2A=root CO2 flux from grosub.f
!
    IF(N.EQ.1.AND.IDAYs1(1,NB1s1(NZ),NZ).GT.0 &
      .AND.RTLGPs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
      RTARRX=RTARR(N,L)/RRADPs1(N,NZ)
      DIFOP=OLSGLP*RTARRX
      VOLWCA=RTVLWs1(N,L,NZ)*SCO2Ls1(L)
      VOLWOA=RTVLWs1(N,L,NZ)*SOXYLs1(L)
      VOLWC4=RTVLWs1(N,L,NZ)*SCH4Ls1(L)
      VOLWZA=RTVLWs1(N,L,NZ)*SN2OLs1(L)
      VOLWNA=RTVLWs1(N,L,NZ)*SNH3Ls1(L)
      VOLWH2=RTVLWs1(N,L,NZ)*SH2GLs1(L)
      DFCOA=CGSGL1*RTCRA
      DFOXA=OGSGL1*RTCRA
      DFCHA=CHSGL1*RTCRA
      DFN2A=Z2SGL1*RTCRA
      DFNHA=ZHSGL1*RTCRA
      DFHGA=HGSGL1*RTCRA
    ELSE
      RTARRX=0.0_r8
      DIFOP=0.0_r8
      VOLWCA=0.0_r8
      VOLWOA=0.0_r8
      VOLWC4=0.0_r8
      VOLWZA=0.0_r8
      VOLWNA=0.0_r8
      VOLWH2=0.0_r8
      DFCOA=0.0_r8
      DFOXA=0.0_r8
      DFCHA=0.0_r8
      DFN2A=0.0_r8
      DFNHA=0.0_r8
      DFHGA=0.0_r8
    ENDIF
    DFGP=AMIN1(1.0,XNPD*SQRT(PORTs1(N,NZ))*TFNDs1(L))
    RCO2PX=-RCO2As1(N,L,NZ)*XNPG
!
!     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
!     AT SMALLER TIME STEP NPH
!
    DO 99 M=1,NPH
!
!     AQUEOUS GAS DIFFUSIVITY THROUGH SOIL WATER TO ROOT
!
!     gas code:CO2=CO2,OXY=O2,CH4=CH4,Z2O=N2O,NH3=NH3 non-band,
!     NHB=NH3 band,H2G=H2
!     VOLWMM,VOLPMM=soil micropore water,air volume
!     FOXYX=root fraction of total O2 demand from previous hour
!     FPQ=PFT fraction of biome root mass
!     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
!     VOLX=soil volume excluding rock,macropores
!     THETW1=soil water concentration
!     TORT=soil tortuosity
!     FILM=soil water film thickness
!     RRADL=root radius
!     RRADS=path length for radial diffusion from soil to root
!     DIF*=aqueous diffusivity from soil to root:OL=O2,CL=CH4
!     ZL=N2O,NL=NH3 non-band,NB=NH4 band,HL=H2
!     C*G=soil gaseous concentration
!     VOLW*,VOLP*=VOLWMM,VOLPMM*gas solubility
!
      VOLWMO=VOLWMs1(M,L)*FOXYX
      VOLWMM=VOLWMs1(M,L)*FPQ(N,L,NZ)
      VOLPMM=VOLPMs1(M,L)*FPQ(N,L,NZ)
      VOLWSP=RTVLWs1(N,L,NZ)+VOLWMM
      VOLWMA=VOLWMM*VLNH4s1(L)
      VOLWMB=VOLWMM*VLNHBs1(L)
      VOLWSA=RTVLWA+VOLWMA
      VOLWSB=RTVLWB+VOLWMB
      THETW1=AMAX1(0.0,VOLWMs1(M,L)/VOLYs1(L))
      IF(THETW1.GT.THETYs1(L).AND.FPQ(N,L,NZ).GT.ZEROQs1(NZ))THEN
        THETM=TORTs1(M,L)*THETW1
        RRADS=LOG((FILMs1(M,L)+RRADL(N,L))/RRADL(N,L))
        RTARRX=RTARR(N,L)/RRADS
        DIFOL=THETM*OLSGL1*RTARRX
        DIFCL=THETM*CQSGL1*RTARRX
        DIFZL=THETM*ZVSGL1*RTARRX
        DIFNL=THETM*ZNSGL1*RTARRX*VLNH4s1(L)
        DIFNB=THETM*ZNSGL1*RTARRX*VLNHBs1(L)
        DIFHL=THETM*HLSGL1*RTARRX
        CH4G1=CCH4Gs1(L)*VOLPMM
        Z2OG1=CZ2OGs1(L)*VOLPMM
        ZH3G1=CNH3Gs1(L)*VOLPMM
        H2GG1=CH2GGs1(L)*VOLPMM
        VOLWCO=VOLWMM*SCO2Ls1(L)
        VOLWOX=VOLWMM*SOXYLs1(L)
        VOLWCH=VOLWMM*SCH4Ls1(L)
        VOLWN2=VOLWMM*SN2OLs1(L)
        VOLWNH=VOLWMM*SNH3Ls1(L)*VLNH4s1(L)
        VOLWNB=VOLWMM*SNH3Ls1(L)*VLNHBs1(L)
        VOLWHG=VOLWMM*SH2GLs1(L)
        VOLPNH=VOLPMM*VLNH4s1(L)
        VOLPNB=VOLPMM*VLNHBs1(L)
!
!     MASS FLOW OF GAS FROM SOIL TO ROOT AT SHORTER TIME STEP NPT
!
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*A1=root gaseous concentration
!     C*P1=root aqueous concentration
!     ROXYLX=soil net O2 aqueous flux
!     VOLWMM=micropore water volume
!     RTVLW,RTVLP=root aqueous,gaseous volume
!     RMF*=soil convective solute flux:COS=CO2,OXS=O2,CHS=CH4,
!     N2S=N2O,NHS=NH3 non-band,NHB=NH3 band,HGS=H2
!     UPWTRH=water uptake
!
        DO 90 MX=1,NPT
          OXYS1=OXYS1+ROXYLX
          CCO2S1=AMAX1(0.0,CO2S1/VOLWMM)
          COXYS1=AMIN1(COXYEs1*SOXYLs1(L),AMAX1(0.0,OXYS1/VOLWMO))
          CCH4S1=AMAX1(0.0,CH4S1/VOLWMM)
          CN2OS1=AMAX1(0.0,Z2OS1/VOLWMM)
          CNH3S1=AMAX1(0.0,ZH3S1/VOLWMM)
          CNH3B1=AMAX1(0.0,ZH3B1/VOLWMM)
          CH2GS1=AMAX1(0.0,H2GS1/VOLWMM)
          IF(RTVLPs1(N,L,NZ).GT.ZEROs1)THEN
            CCO2A1=AMAX1(0.0,CO2A1/RTVLPs1(N,L,NZ))
            COXYA1=AMAX1(0.0,OXYA1/RTVLPs1(N,L,NZ))
            CCH4A1=AMAX1(0.0,CH4A1/RTVLPs1(N,L,NZ))
            CZ2OA1=AMAX1(0.0,Z2OA1/RTVLPs1(N,L,NZ))
            CNH3A1=AMAX1(0.0,ZH3A1/RTVLPs1(N,L,NZ))
            CH2GA1=AMAX1(0.0,H2GA1/RTVLPs1(N,L,NZ))
          ELSE
            CCO2A1=0.0_r8
            COXYA1=0.0_r8
            CCH4A1=0.0_r8
            CZ2OA1=0.0_r8
            CNH3A1=0.0_r8
            CH2GA1=0.0_r8
          ENDIF
          CCO2P1=AMAX1(0.0,CO2P1/RTVLWs1(N,L,NZ))
          COXYP1=AMIN1(COXYEs1*SOXYLs1(L),AMAX1(0.0,OXYP1/RTVLWs1(N,L,NZ)))
          CCH4P1=AMAX1(0.0,CH4P1/RTVLWs1(N,L,NZ))
          CN2OP1=AMAX1(0.0,Z2OP1/RTVLWs1(N,L,NZ))
          CNH3P1=AMAX1(0.0,ZH3P1/RTVLWs1(N,L,NZ))
          CH2GP1=AMAX1(0.0,H2GP1/RTVLWs1(N,L,NZ))
          DIFOX=DIFOL+DIFOP
          RMFCOS=UPWTRH*CCO2S1
          RMFOXS=UPWTRH*COXYS1
          RMFCHS=UPWTRH*CCH4S1
          RMFN2S=UPWTRH*CN2OS1
          RMFN3S=UPWTRH*CNH3S1*VLNH4s1(L)
          RMFN3B=UPWTRH*CNH3B1*VLNHBs1(L)
          RMFHGS=UPWTRH*CH2GS1
!
!     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
!     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
!
!     DIFOL=O2 aqueous diffusivity from soil to root
!     UPWTRH=water uptake
!     DIFOP=aqueous diffusivity of O2 within root
!     COXYS1,COXYP1=soil,root aqueous O2 concentration
!     UPMXP=O2 demand per plant
!     RUPOXR=root O2 uptake per plant
!     COXYR=aqueous O2 concentration at root surface
!     RDFOXS,RDFOXP=aqueous O2 diffusion per plant:soil-root,within root
!
          X=(DIFOL+UPWTRH)*COXYS1+DIFOP*COXYP1
          IF(X.GT.ZEROs1.AND.OXYS1.GT.ZEROPs1(NZ))THEN
            B=-UPMXP-DIFOX*OXKM-X
            C=X*UPMXP
            RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
            COXYR=(X-RUPOXR)/DIFOX
            RDFOXS=RMFOXS+DIFOL*(COXYS1-COXYR)
            RDFOXP=DIFOP*(COXYP1-COXYR)
          ELSE
            X=DIFOP*COXYP1
            IF(X.GT.ZEROs1.AND.OXYP1.GT.ZEROPs1(NZ))THEN
              B=-UPMXP-DIFOP*OXKM-X
              C=X*UPMXP
              RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
              COXYR=(X-RUPOXR)/DIFOP
              RDFOXS=0.0_r8
              RDFOXP=DIFOP*(COXYP1-COXYR)
            ELSE
              RUPOXR=0.0_r8
              COXYR=0.0_r8
              RDFOXS=0.0_r8
              RDFOXP=0.0_r8
            ENDIF
          ENDIF
!
!     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
!     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
!     WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!     RUPOSX,RUPOPX=aqueous O2 uptake from soil,root
!     PP=PFT population
!     RDFCOS,RCO2SX=aqueous CO2 soil-root diffusion,root uptake
!     RDFCHS,RUPCSX=aqueous CH4 soil-root diffusion,root uptake
!     RDFN2S,RUPZSX=aqueous N2O soil-root diffusion,root uptake
!     RDFN3S,RUPNSX=aqueous NH3 soil-root diffusion,root uptake:non-band
!     RDFN3B,RUPNBX=aqueous NH3 soil-root diffusion,root uptake:band
!     RDFHGS,RUPHGX=aqueous H2 soil-root diffusion,root uptake
!     RMF*=soil convective solute flux
!     DIF*=aqueous diffusivity from soil to root
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*P1=root aqueous concentration
!

          RUPOSX=RDFOXS*PPs1(NZ)
          RUPOPX=RDFOXP*PPs1(NZ)
          RDFCOS=RMFCOS+DIFCL*(CCO2S1-CCO2P1)
          RDXCOS=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),CO2S1) &
            -VOLWMM*AMAX1(ZEROPs1(NZ),CO2P1))/VOLWSP
          IF(RDFCOS.GT.0.0)THEN
            RCO2SX=AMIN1(AMAX1(0.0,RDXCOS),RDFCOS*PPs1(NZ))
          ELSE
            RCO2SX=AMAX1(AMIN1(0.0,RDXCOS),RDFCOS*PPs1(NZ))
          ENDIF
          IF(N.EQ.1)THEN
            RDFCHS=RMFCHS+DIFCL*(CCH4S1-CCH4P1)
            RDXCHS=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),CH4S1) &
              -VOLWMM*AMAX1(ZEROPs1(NZ),CH4P1))/VOLWSP
            IF(RDFCHS.GT.0.0)THEN
              RUPCSX=AMIN1(AMAX1(0.0,RDXCHS),RDFCHS*PPs1(NZ))
            ELSE
              RUPCSX=AMAX1(AMIN1(0.0,RDXCHS),RDFCHS*PPs1(NZ))
            ENDIF
            RDFN2S=RMFN2S+DIFZL*(CN2OS1-CN2OP1)
            RDXN2S=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),Z2OS1) &
              -VOLWMM*AMAX1(ZEROPs1(NZ),Z2OP1))/VOLWSP
            IF(RDFN2S.GT.0.0)THEN
              RUPZSX=AMIN1(AMAX1(0.0,RDXN2S),RDFN2S*PPs1(NZ))
            ELSE
              RUPZSX=AMAX1(AMIN1(0.0,RDXN2S),RDFN2S*PPs1(NZ))
            ENDIF
            RDFN3S=RMFN3S+DIFNL*(CNH3S1-CNH3P1)
            IF(VOLWSA.GT.ZEROPs1(NZ))THEN
              ZH3PA=ZH3P1*VLNH4s1(L)
              RDXNHS=(RTVLWA*AMAX1(ZEROPs1(NZ),ZH3S1) &
                -VOLWMA*AMAX1(ZEROPs1(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXNHS=0.0_r8
            ENDIF
            IF(RDFN3S.GT.0.0)THEN
              RUPNSX=AMIN1(AMAX1(0.0,RDXNHS),RDFN3S*PPs1(NZ))
            ELSE
              RUPNSX=AMAX1(AMIN1(0.0,RDXNHS),RDFN3S*PPs1(NZ))
            ENDIF
            RDFN3B=RMFN3B+DIFNB*(CNH3B1-CNH3P1)
            IF(VOLWSB.GT.ZEROPs1(NZ))THEN
              ZH3PB=ZH3P1*VLNHBs1(L)
              RDXNHB=(RTVLWB*AMAX1(ZEROPs1(NZ),ZH3B1) &
                -VOLWMB*AMAX1(ZEROPs1(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXNHB=0.0_r8
            ENDIF
            IF(RDFN3B.GT.0.0)THEN
              RUPNBX=AMIN1(AMAX1(0.0,RDXNHB),RDFN3B*PPs1(NZ))
            ELSE
              RUPNBX=AMAX1(AMIN1(0.0,RDXNHB),RDFN3B*PPs1(NZ))
            ENDIF
            RDFHGS=RMFHGS+DIFHL*(CH2GS1-CH2GP1)
            RDXHGS=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),H2GS1) &
              -VOLWMM*AMAX1(ZEROPs1(NZ),H2GP1))/VOLWSP
            IF(RDFHGS.GT.0.0)THEN
              RUPHGX=AMIN1(AMAX1(0.0,RDXHGS),RDFHGS*PPs1(NZ))
            ELSE
              RUPHGX=AMAX1(AMIN1(0.0,RDXHGS),RDFHGS*PPs1(NZ))
            ENDIF
          ELSE
            RUPCSX=0.0_r8
            RUPZSX=0.0_r8
            RUPNSX=0.0_r8
            RUPNBX=0.0_r8
            RUPHGX=0.0_r8
          ENDIF
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN SOIL
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENTS
!     FROM 'WATSUB'
!
!     THETPM,THETX=air-filled porosity,minimum THETPM
!     R*DFQ=soil gas exchange between gaseous-aqueous phases
!     DFGS=rate constant for soil gas exchange from watsub.f
!     CO2G1,CO2S1=gaseous,aqueous CO2 in soil
!     OXYG1,OXYS1=gaseous,aqueous O2 in soil
!     CH4G1,CH4S1=gaseous,aqueous CH4 in soil
!     Z2OG1,Z2OS1=gaseous,aqueous N2O in soil
!     ZH3G1,ZH3S1,ZH3B1=gaseous,aqueous NH3 in soil non-band,band
!     H2GG1,H2GS1=gaseous,aqueous H2 in soil
!     RUPOSX=root aqueous O2 uptake
!     ROXYLX=root net O2 aqueous flux
!     RCO2SX=root aqueous CO2 uptake
!     RUPCSX=root aqueous CH4 uptake
!     RUPZSX=root aqueous N2O uptake
!     RUPNSX=root aqueous NH3 uptake non-band
!     RUPNBX=root aqueous NH3 uptake band
!     RUPHGX=root aqueous H2 uptake
!     VOLWMM,VOLPMM=soil micropore water,air volume
!     VOLW*=VOLWMM*gas solubility
!
          IF(THETPMs1(M,L).GT.THETX)THEN
            DFGSP=FPQ(N,L,NZ)*DFGSs1(M,L)
            RCODFQ=DFGSP*(AMAX1(ZEROPs1(NZ),CO2G1)*VOLWCO &
              -(AMAX1(ZEROSs1,CO2S1)-RCO2SX)*VOLPMM)/(VOLWCO+VOLPMM)
            RUPOST=RUPOSX-ROXYLX
            ROXDFQ=DFGSP*(AMAX1(ZEROPs1(NZ),OXYG1)*VOLWOX &
              -(AMAX1(ZEROSs1,OXYS1)-RUPOST)*VOLPMM)/(VOLWOX+VOLPMM)
            IF(N.EQ.1)THEN
              RCHDFQ=DFGSP*(AMAX1(ZEROPs1(NZ),CH4G1)*VOLWCH &
                -(AMAX1(ZEROSs1,CH4S1)-RUPCSX)*VOLPMM)/(VOLWCH+VOLPMM)
              RN2DFQ=DFGSP*(AMAX1(ZEROPs1(NZ),Z2OG1)*VOLWN2 &
                -(AMAX1(ZEROSs1,Z2OS1)-RUPZSX)*VOLPMM)/(VOLWN2+VOLPMM)
              IF(VOLWNH+VOLPNH.GT.ZEROPs1(NZ))THEN
                ZH3GA=ZH3G1*VLNH4s1(L)
                RNHDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROPs1(NZ),ZH3GA)*VOLWNH &
                  -(AMAX1(ZEROSs1,ZH3S1)-RUPNSX)*VOLPNH)/(VOLWNH+VOLPNH)))
              ELSE
                RNHDFQ=0.0_r8
              ENDIF
              IF(VOLWNB+VOLPNB.GT.ZEROPs1(NZ))THEN
                ZH3GB=ZH3G1*VLNHBs1(L)
                RNBDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROPs1(NZ),ZH3GB)*VOLWNB &
                  -(AMAX1(ZEROSs1,ZH3B1)-RUPNBX)*VOLPNB)/(VOLWNB+VOLPNB)))
              ELSE
                RNBDFQ=0.0_r8
              ENDIF
              RHGDFQ=DFGSP*(AMAX1(ZEROPs1(NZ),H2GG1)*VOLWHG &
                -(AMAX1(ZEROSs1,H2GS1)-RUPHGX)*VOLPMM)/(VOLWHG+VOLPMM)
            ELSE
              RCHDFQ=0.0_r8
              RN2DFQ=0.0_r8
              RNHDFQ=0.0_r8
              RNBDFQ=0.0_r8
              RHGDFQ=0.0_r8
            ENDIF
          ELSE
            RCODFQ=0.0_r8
            ROXDFQ=0.0_r8
            RCHDFQ=0.0_r8
            RN2DFQ=0.0_r8
            RNHDFQ=0.0_r8
            RNBDFQ=0.0_r8
            RHGDFQ=0.0_r8
          ENDIF
!
!     UPDATE GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
!     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
!
          OXYG1=OXYG1-ROXDFQ+ROXYFX
          OXYS1=OXYS1+ROXDFQ-RUPOSX
          CO2G1=CO2G1-RCODFQ+RCO2FX
          CO2S1=CO2S1+RCODFQ-RCO2SX
          CH4S1=CH4S1+RCHDFQ-RUPCSX
          Z2OS1=Z2OS1+RN2DFQ-RUPZSX
          ZH3S1=ZH3S1+RNHDFQ-RUPNSX
          ZH3B1=ZH3B1+RNBDFQ-RUPNBX
          H2GS1=H2GS1+RHGDFQ-RUPHGX
!
!     GAS TRANSFER THROUGH ROOTS
!
          IF(N.EQ.1.AND.RTVLPs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
            RUPNTX=RUPNSX+RUPNBX
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
!
!     R*DF1=root gas exchange between gaseous-aqueous phases
!     R*FL1=root gas exchange with atmosphere
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     CO2A1,CO2P1=gaseous,aqueous CO2 in root
!     OXYA1,OXYP1=gaseous,aqueous O2 in root
!     CH4A1,CH4P1=gaseous,aqueous CH4 in root
!     Z2OA1,Z2OP1=gaseous,aqueous N2O in root
!     ZH3A1,ZH3P1=gaseous,aqueous NH3 in root
!     H2GA1,H2GP1=gaseous,aqueous H2 in root
!     RTVLW,RTVLP=root aqueous,gaseous volume
!     VOLW*=RTVLW*gas solubility
!     C*E,C*A1=atmosphere,root gas concentration
!     DF*A=root-atmosphere gas conductance
!
            CO2PX=CO2P1+RCO2PX
            RCODF1=AMAX1(-CO2PX,DFGP*(AMAX1(ZEROPs1(NZ),CO2A1)*VOLWCA &
              -CO2PX*RTVLPs1(N,L,NZ))/(VOLWCA+RTVLPs1(N,L,NZ)))
            OXYPX=OXYP1-RUPOPX
            ROXDF1=AMAX1(-OXYPX,DFGP*(AMAX1(ZEROPs1(NZ),OXYA1)*VOLWOA &
              -OXYPX*RTVLPs1(N,L,NZ))/(VOLWOA+RTVLPs1(N,L,NZ)))
            CH4PX=CH4P1+RUPCSX
            RCHDF1=AMAX1(-CH4PX,DFGP*(AMAX1(ZEROPs1(NZ),CH4A1)*VOLWC4 &
              -CH4PX*RTVLPs1(N,L,NZ))/(VOLWC4+RTVLPs1(N,L,NZ)))
            Z2OPX=Z2OP1+RUPZSX
            RN2DF1=AMAX1(-Z2OPX,DFGP*(AMAX1(ZEROPs1(NZ),Z2OA1)*VOLWZA &
              -Z2OPX*RTVLPs1(N,L,NZ))/(VOLWZA+RTVLPs1(N,L,NZ)))
            ZH3PX=ZH3P1+RUPNTX
            RNHDF1=AMAX1(-ZH3PX,DFGP*(AMAX1(ZEROPs1(NZ),ZH3A1)*VOLWNA &
              -ZH3PX*RTVLPs1(N,L,NZ))/(VOLWNA+RTVLPs1(N,L,NZ)))
            H2GPX=H2GP1+RUPHGX
            RHGDF1=AMAX1(-H2GPX,DFGP*(AMAX1(ZEROPs1(NZ),H2GA1)*VOLWH2 &
              -H2GPX*RTVLPs1(N,L,NZ))/(VOLWH2+RTVLPs1(N,L,NZ)))
            RCOFL1=AMIN1(DFCOA,RTVLPs1(N,L,NZ))*(CCO2Es1-CCO2A1)
            ROXFL1=AMIN1(DFOXA,RTVLPs1(N,L,NZ))*(COXYEs1-COXYA1)
            RCHFL1=AMIN1(DFCHA,RTVLPs1(N,L,NZ))*(CCH4Es1-CCH4A1)
            RN2FL1=AMIN1(DFN2A,RTVLPs1(N,L,NZ))*(CZ2OEs1-CZ2OA1)
            RNHFL1=AMIN1(DFNHA,RTVLPs1(N,L,NZ))*(CNH3Es1-CNH3A1)
            RHGFL1=AMIN1(DFHGA,RTVLPs1(N,L,NZ))*(CH2GEs1-CH2GA1)
          ELSE
            RCODF1=0.0_r8
            ROXDF1=0.0_r8
            RCHDF1=0.0_r8
            RN2DF1=0.0_r8
            RNHDF1=0.0_r8
            RHGDF1=0.0_r8
            RCOFL1=0.0_r8
            ROXFL1=0.0_r8
            RCHFL1=0.0_r8
            RN2FL1=0.0_r8
            RNHFL1=0.0_r8
            RHGFL1=0.0_r8
          ENDIF
!
!     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
!     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
!
          CO2A1=CO2A1-RCODF1+RCOFL1
          OXYA1=OXYA1-ROXDF1+ROXFL1
          CH4A1=CH4A1-RCHDF1+RCHFL1
          Z2OA1=Z2OA1-RN2DF1+RN2FL1
          ZH3A1=ZH3A1-RNHDF1+RNHFL1
          H2GA1=H2GA1-RHGDF1+RHGFL1
          CO2P1=CO2P1+RCODF1+RCO2SX+RCO2PX
          OXYP1=OXYP1+ROXDF1-RUPOPX
          CH4P1=CH4P1+RCHDF1+RUPCSX
          Z2OP1=Z2OP1+RN2DF1+RUPZSX
          ZH3P1=ZH3P1+RNHDF1+RUPNSX+RUPNBX
          H2GP1=H2GP1+RHGDF1+RUPHGX
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2S=soil-root CO2 exchange
!     RUPOXS=soil-root O2 exchange
!     RUPCHS=soil-root CH4 exchange
!     RUPN2S=soil-root N2O exchange
!     RUPN3S=soil-root NH3 exchange non-band
!     RUPN3B=soil-root NH3 exchange band
!     RUPHGS=soil-root H2 exchange
!
          RCO2Ss1(N,L,NZ)=RCO2Ss1(N,L,NZ)+RCO2SX
          RUPOXSs1(N,L,NZ)=RUPOXSs1(N,L,NZ)+RUPOSX
          RUPCHSs1(N,L,NZ)=RUPCHSs1(N,L,NZ)+RUPCSX
          RUPN2Ss1(N,L,NZ)=RUPN2Ss1(N,L,NZ)+RUPZSX
          RUPN3Ss1(N,L,NZ)=RUPN3Ss1(N,L,NZ)+RUPNSX
          RUPN3Bs1(N,L,NZ)=RUPN3Bs1(N,L,NZ)+RUPNBX
          RUPHGSs1(N,L,NZ)=RUPHGSs1(N,L,NZ)+RUPHGX
!
!     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE
!
!     R*DFA=root aqueous-gaseous CO2 exchange
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!
          RCODFAs1(N,L,NZ)=RCODFAs1(N,L,NZ)+RCODF1
          ROXDFAs1(N,L,NZ)=ROXDFAs1(N,L,NZ)+ROXDF1
          RCHDFAs1(N,L,NZ)=RCHDFAs1(N,L,NZ)+RCHDF1
          RN2DFAs1(N,L,NZ)=RN2DFAs1(N,L,NZ)+RN2DF1
          RNHDFAs1(N,L,NZ)=RNHDFAs1(N,L,NZ)+RNHDF1
          RHGDFAs1(N,L,NZ)=RHGDFAs1(N,L,NZ)+RHGDF1
          RCOFLAs1(N,L,NZ)=RCOFLAs1(N,L,NZ)+RCOFL1
          ROXFLAs1(N,L,NZ)=ROXFLAs1(N,L,NZ)+ROXFL1
          RCHFLAs1(N,L,NZ)=RCHFLAs1(N,L,NZ)+RCHFL1
          RN2FLAs1(N,L,NZ)=RN2FLAs1(N,L,NZ)+RN2FL1
          RNHFLAs1(N,L,NZ)=RNHFLAs1(N,L,NZ)+RNHFL1
          RHGFLAs1(N,L,NZ)=RHGFLAs1(N,L,NZ)+RHGFL1
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2P=root CO2 emission into root
!     RUPOXP=root O2 uptake from root
!     ROXSK=total O2 uptake from soil by all microbial,root popns
!
          RCO2Ps1(N,L,NZ)=RCO2Ps1(N,L,NZ)+RCO2PX+RCO2SX
          RUPOXPs1(N,L,NZ)=RUPOXPs1(N,L,NZ)+RUPOPX
          ROXSKs1(M,L)=ROXSKs1(M,L)+RUPOSX

90      CONTINUE
      ENDIF
99  CONTINUE
!
!     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
!     OF ROOT O2 UPTAKE 'RUPOXT' TO ROOT O2 DEMAND 'ROXYP'
!
!     RUPOXT=O2 uptake from soil+root by each root,myco population
!     ROXYP=O2 demand by each root,myco population
!     WFR=constraint by O2 consumption on all root processes
!     imposed by O2 uptake
!
    RUPOXT=RUPOXPs1(N,L,NZ)+RUPOXSs1(N,L,NZ)
    WFRs1(N,L,NZ)=AMIN1(1.0,AMAX1(0.0 &
      ,RUPOXT/ROXYPs1(N,L,NZ)))
  ELSE
    RUPOXT=0.0_r8
    IF(L.GT.NGs1(NZ))THEN
      WFRs1(N,L,NZ)=WFRs1(N,L-1,NZ)
    ELSE
      WFRs1(N,L,NZ)=1.0
    ENDIF
  ENDIF
  end associate
  end subroutine RootSoilGasExchange
!------------------------------------------------------------------------------------------

  subroutine RootMycoO2NutrientUptake(NZ,OSTRN,OSTRD,PATH,RRADL,&
    FPQ,FPP,FRTDPX,RTARR)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(in) :: FRTDPX(JZ1,JP1),RTARR(2,JZ1)
  real(r8), intent(inout) :: OSTRN,OSTRD
  real(r8) :: TFOXYX
  real(r8) :: FCUP,FZUP,FPUP,FWSRT,UPWTRP,UPWTRH,FOXYX,RUPOXT
  integer :: N,L
!     begin_execution
  associate(                             &
    THETWs1 =>  plt_soilchem%THETWs1   , &
    VOLXs1  =>  plt_soilchem%VOLXs1    , &
    MYs1    =>  plt_morph%MYs1         , &
    NIs1    =>  plt_morph%NIs1           &
  )
  DO 955 N=1,MYs1(NZ)
    DO 950 L=NUs1,NIs1(NZ)
      IF(VOLXs1(L).GT.ZEROS2s1 &
        .AND.RTDNPs1(N,L,NZ).GT.ZEROs1 &
        .AND.RTVLWs1(N,L,NZ).GT.ZEROPs1(NZ) &
        .AND.THETWs1(L).GT.ZEROs1)THEN
        TFOXYX=0.0_r8
        call GetUptakeCapcity(N,L,NZ,FPQ,FPP,FCUP,FZUP,FPUP,&
          FWSRT,UPWTRP,UPWTRH,FOXYX)

        TFOXYX=TFOXYX+FOXYX
!
!     ROOT O2 DEMAND CALCULATED FROM O2 NON-LIMITED RESPIRATION RATE
!
!     ROXYP=O2 demand
!     RCO2M=respiration unlimited by O2
!     RTVLW=root or myco aqueous volume
!     FOXYX=fraction of total O2 demand from previous hour
!
        ROXYPs1(N,L,NZ)=2.667*RCO2Ms1(N,L,NZ)

        call RootSoilGasExchange(N,L,NZ,RRADL,FPQ,FRTDPX,RTARR,UPWTRH,&
          FOXYX,RUPOXT)

        OSTRD=OSTRD+ROXYPs1(N,L,NZ)
        OSTRN=OSTRN+RUPOXT

        call RootExudates(N,L,NZ)
!
!     NUTRIENT UPTAKE
!
!     WFR=constraint by O2 consumption on all biological processes
!     FCUP=limitation to active uptake respiration from CPOOLR
!     FWSRT=protein concentration relative to 5%
!     RTLGP=root,myco length per plant
!
        IF(WFRs1(N,L,NZ).GT.ZEROs1.AND.FCUP.GT.ZEROs1.AND.FWSRT.GT.ZEROs1 &
          .AND.RTLGPs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
!
!     FZUP=limitn to active uptake respiration from CZPOLR
!
          call UptakeMineralNitrogen(N,L,NZ,PATH,RRADL,FPQ,FPP,&
            RTARR,FCUP,FZUP,FWSRT,UPWTRP)
!
!     FPUP=limitn to active uptake respiration from CPPOLR
!
          call UptakeMineralPhosporhus(N,L,NZ,PATH,RRADL,FPQ,FPP,&
            RTARR,FCUP,FPUP,FWSRT,UPWTRP)

        ELSE
          call NoActiveNutrientUptake(N,L,NZ)

        ENDIF
      ELSE
        call ZeroUptake(N,L,NZ)
      ENDIF

      call SumupNutrientUptake(N,L,NZ)

950 CONTINUE
955 CONTINUE
  end associate
  end subroutine RootMycoO2NutrientUptake

end module UptakesMod
