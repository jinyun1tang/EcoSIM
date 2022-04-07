module UptakeMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StomateMod   , only : stomate
  use minimathmod  , only : safe_adb,vapsat,test_aneb
  use EcosimConst
  use SOMDataType
  use ChemTranspDataType
  use UptakePars
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use PlantDataRateType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use LandSurfDataType
  use PlantTraitDataType
  use SnowDataType
  use RootDataType
  use CanopyDataType
  use EcosysBGCFluxType
  use SoilPropertyDataType
  use CanopyRadDataType
  use GridDataType
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__

  public :: uptake
  public :: InitUptake
  contains


  subroutine InitUptake

  implicit none

  call InitUptakePars
!------------------------------------------------------------------------------------------

  end subroutine InitUptake

  subroutine uptake(I,J,NHW,NHE,NVN,NVS)
!
!     THIS subroutine CALCULATES EXCHANGES OF ENERGY, C, N AND P
!     BETWEEN THE CANOPY AND THE ATMOSPHERE AND BETWEEN ROOTS AND THE SOIL
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NN,N,NX,NY,NZ,K,L
  real(r8) :: FPC
  real(r8) :: OSTRN,OSTRD,PARHC
  real(r8) :: WVPLT
  real(r8) :: PSIST1(JZ),PATH(2,JZ)
  real(r8) :: RRADL(2,JZ),RSRT(2,JZ)
  real(r8) :: RSRG(2,JZ),RSR1(2,JZ)
  real(r8) :: RSR2(2,JZ),RSSX(2,JZ)
  real(r8) :: RSRS(2,JZ),WTRTG(JZ)
  real(r8) :: FPQ(2,JZ,05),FPP(2,JZ,05)
  real(r8) :: FRTDPX(JZ,05),RTARR(2,JZ)
  real(r8) :: VOLPU(JZ),VOLWU(JZ)
  real(r8) :: TKCX,CNDT,HFLWC1,VHCPX
  real(r8) :: DIFF,PSILH,UPRT,VFLXC
  real(r8) :: FDMP
  integer :: ILYR(2,JZ)
!     begin_execution

  DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS

      call PrepUptake(NY,NX,PSIST1,WTRTG,VOLPU,VOLWU)
!
!     IF PLANT SPECIES EXISTS
!
      DO 9985 NZ=1,NP(NY,NX)
        OSTRN=0.0_r8
        OSTRD=0.0_r8
        IF(IFLGC(NZ,NY,NX).EQ.1.AND.PP(NZ,NY,NX).GT.0.0)THEN

          call UpdateCanopyProperty(NZ,NY,NX)

!     STOMATE=solve for minimum canopy stomatal resistance
          CALL STOMATE(I,J,NZ,NY,NX)
!
!     CALCULATE VARIABLES USED IN ROOT UPTAKE OF WATER AND NUTRIENTS
          call UpdateRootProperty(NZ,NY,NX,PATH,RRADL,WTRTG,FPQ,FPP, &
            FRTDPX,RTARR)
!
!     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
!     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER CONTENT
!
!     (AG: - originally this line had a N0B1 here )
          IF((IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0) &
            .AND.(ARLFS(NZ,NY,NX).GT.ZEROL(NZ,NY,NX) &
            .AND.FRADP(NZ,NY,NX).GT.0.0) &
            .AND.(RTDP1(1,1,NZ,NY,NX).GT.SDPTH(NZ,NY,NX)+CDPTHZ(0,NY,NX)))THEN
!
            call CalcResistance(NZ,NY,NX,PATH,RRADL,RTARR,RSRT,&
              RSRG,RSR1,RSR2,RSSX,RSRS,CNDT,PSILH,ILYR)
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
            PSILT(NZ,NY,NX)=AMIN1(-1.0E-06,0.667*PSILT(NZ,NY,NX))
            EP(NZ,NY,NX)=0.0_r8
            EVAPC(NZ,NY,NX)=0.0_r8
            HFLWC1=FLWC(NZ,NY,NX)*cpw*TKA(NY,NX)

            IF(ARLSS(NY,NX).GT.ZEROS(NY,NX))THEN
              FPC=ARLFS(NZ,NY,NX)/ARLSS(NY,NX)*AMIN1(1.0,0.5*ARLFC(NY,NX)/AREA(3,NU(NY,NX),NY,NX))
            ELSEIF(PPT(NY,NX).GT.ZEROS(NY,NX))THEN
              FPC=PP(NZ,NY,NX)/PPT(NY,NX)
            ELSE
              FPC=1.0/NP(NY,NX)
            ENDIF

            TKCX=TKC(NZ,NY,NX)
            WVPLT=AMAX1(0.0,WTLS(NZ,NY,NX)+WVSTK(NZ,NY,NX))
            VHCPX=cpw*(WVPLT*VSTK+VOLWC(NZ,NY,NX)+VOLWP(NZ,NY,NX))
!
!     CONVERGENCE SOLUTION
!
            NN=CanopyEnergyH2OIteration(I,J,NZ,NY,NX,FPC,WVPLT,&
              PSIST1,PARHC,DIFF,UPRT,VFLXC,FDMP,RSRS,FPQ,VOLPU,VOLWU,TKCX,CNDT,VHCPX,&
              HFLWC1,PSILH,ILYR)
!
!     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
!
!     TKC=final estimate of canopy temperature TKCZ
!     TKA=current air temperature
!     DTKC=TKC-TKA for next hour
!
            TKC(NZ,NY,NX)=TKCZ(NZ,NY,NX)
            TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-TC2K
            DTKC(NZ,NY,NX)=TKC(NZ,NY,NX)-TKA(NY,NX)
!
!     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
!     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
!
            call HandlingDivergence(I,J,NN,NZ,NY,NX,PSIST1,DIFF,FDMP)

            call UpdateCanopyWater(NZ,NY,NX,PARHC,PSIST1,RSRT,RSSX,&
              RSRS,TKCX,VHCPX,HFLWC1,UPRT,VFLXC,ILYR)
!
!     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
!
          ELSE
            call HandleBareSoil(NZ,NY,NX,PSIST1,FDMP)
          ENDIF

          call SetCanopyGrowthFuncs(NZ,NY,NX)

          call CanopyNH3Flux(NZ,NY,NX,FDMP)
!
!     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
!
          call RootMycoO2NutrientUptake(NZ,NY,NX,OSTRN,OSTRD,PATH,RRADL,&
            FPQ,FPP,FRTDPX,RTARR)

          TLEC(NY,NX)=TLEC(NY,NX)+EFLXC(NZ,NY,NX)*RA(NZ,NY,NX)
          TSHC(NY,NX)=TSHC(NY,NX)+SFLXC(NZ,NY,NX)*RA(NZ,NY,NX)
          IF(OSTRD.GT.ZEROP(NZ,NY,NX))THEN
            OSTR(NZ,NY,NX)=OSTRN/OSTRD
          ELSE
            OSTR(NZ,NY,NX)=0.0_r8
          ENDIF
        ENDIF
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
  RETURN
  END subroutine uptake
!------------------------------------------------------------------------

  subroutine PrepUptake(NY,NX,PSIST1,WTRTG,VOLPU,VOLWU)
!
!     prepare for uptake calculation
  implicit none
  integer, intent(in) :: NY, NX
  real(r8), intent(out) :: PSIST1(JZ),WTRTG(JZ),VOLPU(JZ),VOLWU(JZ)
  integer :: NZ, L, N
  real(r8) :: ARLSC

!
!     RESET TOTAL UPTAKE ARRAYS
!
!     ARLFP,ARSTP=leaf,stalk areas
!

  ARLSC=0.0_r8
  DO 9984 NZ=1,NP0(NY,NX)
!     TKC(NZ,NY,NX)=TKA(NY,NX)+DTKC(NZ,NY,NX)
!     TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-TC2K
    ARLSC=ARLSC+ARLFP(NZ,NY,NX)+ARSTP(NZ,NY,NX)
    RAD1(NZ,NY,NX)=0.0_r8
    EFLXC(NZ,NY,NX)=0.0_r8
    SFLXC(NZ,NY,NX)=0.0_r8
    HFLXC(NZ,NY,NX)=0.0_r8
    THRM1(NZ,NY,NX)=0.0_r8
    EP(NZ,NY,NX)=0.0_r8
    EVAPC(NZ,NY,NX)=0.0_r8
    UPOMC(NZ,NY,NX)=0.0_r8
    UPOMN(NZ,NY,NX)=0.0_r8
    UPOMP(NZ,NY,NX)=0.0_r8
    UPNH4(NZ,NY,NX)=0.0_r8
    UPNO3(NZ,NY,NX)=0.0_r8
    UPH2P(NZ,NY,NX)=0.0_r8
    UPH1P(NZ,NY,NX)=0.0_r8
    UPNF(NZ,NY,NX)=0.0_r8
!
!     RESET UPTAKE ARRAYS
!
    DO  L=NU(NY,NX),NJ(NY,NX)
      DO  N=1,MY(NZ,NY,NX)
        UPWTR(N,L,NZ,NY,NX)=0.0_r8
        RCO2P(N,L,NZ,NY,NX)=0.0_r8
        RUPOXP(N,L,NZ,NY,NX)=0.0_r8
        RCO2S(N,L,NZ,NY,NX)=0.0_r8
        RUPOXS(N,L,NZ,NY,NX)=0.0_r8
        RUPCHS(N,L,NZ,NY,NX)=0.0_r8
        RUPN2S(N,L,NZ,NY,NX)=0.0_r8
        RUPN3S(N,L,NZ,NY,NX)=0.0_r8
        RUPN3B(N,L,NZ,NY,NX)=0.0_r8
        RUPHGS(N,L,NZ,NY,NX)=0.0_r8
        RCOFLA(N,L,NZ,NY,NX)=0.0_r8
        ROXFLA(N,L,NZ,NY,NX)=0.0_r8
        RCHFLA(N,L,NZ,NY,NX)=0.0_r8
        RN2FLA(N,L,NZ,NY,NX)=0.0_r8
        RNHFLA(N,L,NZ,NY,NX)=0.0_r8
        RHGFLA(N,L,NZ,NY,NX)=0.0_r8
        RCODFA(N,L,NZ,NY,NX)=0.0_r8
        ROXDFA(N,L,NZ,NY,NX)=0.0_r8
        RCHDFA(N,L,NZ,NY,NX)=0.0_r8
        RN2DFA(N,L,NZ,NY,NX)=0.0_r8
        RNHDFA(N,L,NZ,NY,NX)=0.0_r8
        RHGDFA(N,L,NZ,NY,NX)=0.0_r8
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
  DO 9000 L=NU(NY,NX),NJ(NY,NX)
    PSIST1(L)=PSIST(L,NY,NX)-0.0098*ALT(NY,NX)
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLWU(L)=VOLWM(NPH,L,NY,NX)-THETY(L,NY,NX)*VOLY(L,NY,NX)
      VOLPU(L)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX)-VOLI(L,NY,NX))
    ELSE
      VOLWU(L)=VOLWM(NPH,L,NY,NX)
      VOLPU(L)=0.0_r8
    ENDIF
    WTRTG(L)=0.0_r8
    DO 9005 NZ=1,NP(NY,NX)
      DO  N=1,MY(NZ,NY,NX)
!     IF(IFLGC(NZ,NY,NX).EQ.1.AND.PP(NZ,NY,NX).GT.0.0)THEN
      WTRTG(L)=WTRTG(L)+AMAX1(0.0,WTRTD(N,L,NZ,NY,NX))
!     ENDIF
      enddo
9005  CONTINUE
9000  CONTINUE
  end subroutine PrepUptake
!------------------------------------------------------------------------
  subroutine UpdateCanopyProperty(NZ,NY,NX)
!
!     update canopy characterization
  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8) :: ALFZ
  real(r8) :: TFRADP,RACZ(JP,JY,JX)
  integer :: NB,K,L,N,NZZ
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!
  DO 500 NB=1,NBR(NZ,NY,NX)
    DO 550 K=1,JNODS
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     ARLF=leaf area
!     WSLF=leaf protein content
!     SURFX,SURF=unself-shaded,total leaf surface area
!     CFX=clumping factor from PFT file
!
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WSLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        KLEAFX(NB,NZ,NY,NX)=K
      ENDIF
      DO 600 L=JC,1,-1
        DO 650 N=1,JLI
          SURFX(N,L,K,NB,NZ,NY,NX)=SURF(N,L,K,NB,NZ,NY,NX)*CFX(NZ,NY,NX)
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
  IF(ARLFS(NZ,NY,NX).GT.0.0)THEN
    IF(IETYP(NY,NX).GE.0)THEN
      TFRADP=0.0_r8
      DO 700 NZZ=1,NP(NY,NX)
        IF(ZC(NZZ,NY,NX).GT.ZC(NZ,NY,NX)+ZR(NY,NX))THEN
          TFRADP=TFRADP+FRADP(NZZ,NY,NX)
        ENDIF
700   CONTINUE
      ALFZ=2.0*TFRADP
      IF(RAB(NY,NX).GT.ZERO.AND.ZT(NY,NX).GT.ZERO.AND.ALFZ.GT.ZERO)THEN
        RACZ(NZ,NY,NX)=AMIN1(RACX,AMAX1(0.0,ZT(NY,NX)*EXP(ALFZ) &
          /(ALFZ/RAB(NY,NX))*(EXP(-ALFZ*ZC(NZ,NY,NX)/ZT(NY,NX)) &
          -EXP(-ALFZ*(ZD(NY,NX)+ZR(NY,NX))/ZT(NY,NX)))))
      ELSE
        RACZ(NZ,NY,NX)=0.0_r8
      ENDIF
    ELSE
      RACZ(NZ,NY,NX)=0.0_r8
    ENDIF
  ELSE
    RACZ(NZ,NY,NX)=RACX
  ENDIF
  RAZ(NZ,NY,NX)=RAB(NY,NX)+RACZ(NZ,NY,NX)
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
  TKCZ(NZ,NY,NX)=TKA(NY,NX)+DTKC(NZ,NY,NX)
  end subroutine UpdateCanopyProperty
!------------------------------------------------------------------------
  subroutine UpdateRootProperty(NZ,NY,NX,PATH,RRADL,WTRTG,FPQ,FPP,FRTDPX,RTARR)
!
!     update root characterization

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: WTRTG(JZ)
  real(r8), intent(out) :: PATH(2,JZ),RRADL(2,JZ)
  real(r8), intent(out) :: FPQ(2,JZ,05),FPP(2,JZ,05)
  real(r8), intent(out) :: FRTDPX(JZ,05),RTARR(2,JZ)
  real(r8) :: RTDPZ,RTDPX
  integer :: N,L,NR
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
  DO 2000 N=1,MY(NZ,NY,NX)
    DO  L=NU(NY,NX),NI(NZ,NY,NX)
      IF(N.EQ.1)THEN
        RTDPZ=0.0_r8
        DO 2005 NR=1,NRT(NZ,NY,NX)
          RTDPZ=AMAX1(RTDPZ,RTDP1(1,NR,NZ,NY,NX))
2005    CONTINUE
        IF(L.EQ.NU(NY,NX))THEN
          FRTDPX(L,NZ)=1.0
        ELSE
          IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
            RTDPX=AMAX1(0.0,RTDPZ-CDPTHZ(L-1,NY,NX))
            RTDPX=AMAX1(0.0,AMIN1(DLYR(3,L,NY,NX),RTDPX) &
              -AMAX1(0.0,SDPTH(NZ,NY,NX)-CDPTHZ(L-1,NY,NX)-HTCTL(NZ,NY,NX)))
            FRTDPX(L,NZ)=RTDPX/DLYR(3,L,NY,NX)
          ELSE
            FRTDPX(L,NZ)=0.0_r8
          ENDIF
        ENDIF

      ENDIF
      IF(WTRTG(L).GT.ZEROS(NY,NX))THEN
        FPQ(N,L,NZ)=AMAX1(0.0,WTRTD(N,L,NZ,NY,NX))/WTRTG(L)
      ELSE
        FPQ(N,L,NZ)=1.0
      ENDIF
      FPP(N,L,NZ)=FMN*FPQ(N,L,NZ)
      IF(RTDNP(N,L,NZ,NY,NX).GT.ZERO.AND.FRTDPX(L,NZ).GT.ZERO)THEN
        RRADL(N,L)=AMAX1(RRAD2X(N,NZ,NY,NX),SQRT((RTVLW(N,L,NZ,NY,NX) &
          /(1.0-PORT(N,NZ,NY,NX)))/(PICON*PP(NZ,NY,NX)*RTLGP(N,L,NZ,NY,NX))))
        PATH(N,L)=AMAX1(1.001*RRADL(N,L) &
          ,1.0/(SQRT(PICON*(RTDNP(N,L,NZ,NY,NX)/FRTDPX(L,NZ))/FMPR(L,NY,NX))))
        RTARR(N,L)=6.283*RTLGP(N,L,NZ,NY,NX)/FRTDPX(L,NZ)
      ELSE
        RRADL(N,L)=RRAD2M(N,NZ,NY,NX)
        PATH(N,L)=1.001*RRADL(N,L)
        RTARR(N,L)=6.283*RTLGP(N,L,NZ,NY,NX)
      ENDIF
    enddo
2000  CONTINUE
  end subroutine UpdateRootProperty
!------------------------------------------------------------------------

  subroutine HandlingDivergence(I,J,NN,NZ,NY,NX,PSIST1,DIFF,FDMP)

  implicit none
  integer  , intent(in) :: NN, I, J
  integer  , intent(in) :: NZ,NY,NX
  REAL(R8) ,INTENT(IN) :: PSIST1(JZ),DIFF
  real(r8), intent(inout):: FDMP
  real(r8) :: APSILT,APSIRT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM,FDMR
  real(r8) :: OSWT,WFNC

! begin_execution
  integer :: N,L

  IF(NN.GE.MXN)THEN
    WRITE(*,9999)IYRC,I,J,NX,NY,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)
    IF(DIFF.GT.0.5)THEN
      RAD1(NZ,NY,NX)=0.0_r8
      EFLXC(NZ,NY,NX)=0.0_r8
      SFLXC(NZ,NY,NX)=0.0_r8
      HFLXC(NZ,NY,NX)=0.0_r8
      EVAPC(NZ,NY,NX)=0.0_r8
      EP(NZ,NY,NX)=0.0_r8
      TKC(NZ,NY,NX)=TKA(NY,NX)+DTKC(NZ,NY,NX)
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-TC2K
      FTHRM=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      THRM1(NZ,NY,NX)=FTHRM*TKC(NZ,NY,NX)**4
      PSILT(NZ,NY,NX)=PSIST1(NG(NZ,NY,NX))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILO(NZ,NY,NX)=FDMP/0.16*OSMO(NZ,NY,NX) &
        -8.3143*TKC(NZ,NY,NX)*FDMP*CCPOLT/OSWT
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)-RSMN(NZ,NY,NX))*WFNC
      RA(NZ,NY,NX)=RAZ(NZ,NY,NX)
      VHCPC(NZ,NY,NX)=cpw*(WTSHT(NZ,NY,NX)*10.0E-06)
      DTKC(NZ,NY,NX)=0.0_r8
      DO 4290 N=1,MY(NZ,NY,NX)
        DO  L=NU(NY,NX),NI(NZ,NY,NX)
          PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
          APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
          FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
          CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)+CPPOLR(N,L,NZ,NY,NX)
          OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
          PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX) &
            -8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
          PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)-PSIRO(N,L,NZ,NY,NX))
          UPWTR(N,L,NZ,NY,NX)=0.0_r8
      enddo
4290  CONTINUE
    ENDIF
  ENDIF
  end subroutine HandlingDivergence
!------------------------------------------------------------------------------
  function CanopyEnergyH2OIteration(I,J,NZ,NY,NX,FPC,WVPLT,&
    PSIST1,PARHC,DIFF,UPRT,VFLXC,FDMP,RSRS,FPQ,VOLPU,VOLWU,TKCX,CNDT,&
    VHCPX,HFLWC1,PSILH,ILYR) result(NN)
  implicit none
  integer  , intent(in) :: I, J
  integer  , intent(in) :: NZ,NY,NX
  real(r8) , intent(in) :: FPC,WVPLT,PSIST1(JZ)
  real(r8) , intent(in) :: RSRS(2,JZ),FPQ(2,JZ,05),VOLPU(JZ),VOLWU(JZ)
  real(r8) , intent(in) :: TKCX,CNDT,VHCPX,HFLWC1,PSILH
  integer  , intent(in) :: ILYR(2,JZ)
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

  CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
  OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
  FTHRM=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
  FDTHS=(THS(NY,NX)+THRMGX(NY,NX))*FRADP(NZ,NY,NX)
!     RAZ=canopy isothermal boundary later resistance

  UPRT=0.0_r8
  PAREX=FPC*AREA(3,NU(NY,NX),NY,NX)
  PARHX=FPC*AREA(3,NU(NY,NX),NY,NX)*1.25E-03
  RA1=RAZ(NZ,NY,NX)

  IC=0
  XC=0.5
  ICHK=0
  PSIL2=0.0_r8
  EPX=0.0_r8
  UPRTX=0.0_r8
  VOLWPX=0.0_r8

  DO 4000 NN=1,MXN
!
!     NET RADIATION FROM ABSORBED SW AND NET LW
!
!     THRM1=LW emitted by canopy
!     DTHS1=net LW absorbed by canopy
!     RADC=total SW absorbed by canopy
!     RAD1=net SW+LW absorbed by canopy
!
    TKC1=TKCZ(NZ,NY,NX)
    THRM1(NZ,NY,NX)=FTHRM*TKC1**4
    DTHS1=FDTHS-THRM1(NZ,NY,NX)*2.0
    RAD1(NZ,NY,NX)=RADC(NZ,NY,NX)+DTHS1
!
!     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
!
!     RI=Ricardson's number
!     RA=canopy boundary layer resistance
!     PAREC,PARHC=canopy latent,sensible heat conductance
!
    RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB(NY,NX)*(TKA(NY,NX)-TKC1)))

    RA(NZ,NY,NX)=AMAX1(RACM,0.9_r8*RA1,AMIN1(1.1_r8*RA1,RAZ(NZ,NY,NX)/(1.0_r8-10.0_r8*RI)))
    RA1=RA(NZ,NY,NX)
    PAREC=PAREX/RA(NZ,NY,NX)
    PARHC=PARHX/RA(NZ,NY,NX)
!
!     CANOPY WATER AND OSMOTIC POTENTIALS
!
!     PSILT=canopy total water potential
!     FDMP=dry matter content
!     OSMO=osmotic potential at PSILT=0 from PFT file
!     PSILO,PSILG=canopy osmotic,turgor water potential
!
    APSILT=ABS(PSILT(NZ,NY,NX))
    FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
    PSILO(NZ,NY,NX)=FDMP/0.16*OSMO(NZ,NY,NX)-8.3143*TKC1*FDMP*CCPOLT/OSWT
    PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
!
!     CANOPY STOMATAL RESISTANCE
!
!     RCS=shape parameter for RC vs PSILG from PFT file
!     RC=canopy stomatal resistance
!     RSMN=minimum RC at PSILT=0 from stomate.f
!     RSMX=cuticular resistance from PFT file
!
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)-RSMN(NZ,NY,NX))*WFNC
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
      vpc=vapsat(tkc1)*EXP(18.0*PSILT(NZ,NY,NX)/(8.3143*TKC1))
      EX=PAREC*(VPA(NY,NX)-VPC)
      IF(EX.GT.0.0)THEN
      EVAPC(NZ,NY,NX)=EX*RA(NZ,NY,NX)/(RA(NZ,NY,NX)+RZ)
      EX=0.0_r8
      ELSEIF(EX.LE.0.0.AND.VOLWC(NZ,NY,NX).GT.0.0)THEN
      EVAPC(NZ,NY,NX)=AMAX1(EX*RA(NZ,NY,NX)/(RA(NZ,NY,NX)+RZ) &
      ,-VOLWC(NZ,NY,NX))
      EX=EX-EVAPC(NZ,NY,NX)
      ENDIF
      EP(NZ,NY,NX)=EX*RA(NZ,NY,NX)/(RA(NZ,NY,NX)+RC(NZ,NY,NX))
      EFLXC(NZ,NY,NX)=(EP(NZ,NY,NX)+EVAPC(NZ,NY,NX))*VAP
      VFLXC=EVAPC(NZ,NY,NX)*cpw*TKC1
!
!     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
!
!     HFLXS=initial estimate of sensible+storage heat flux
!     HFLWC1=convective heat flux from precip to canopy
!
      HFLXS=RAD1(NZ,NY,NX)+EFLXC(NZ,NY,NX)+VFLXC+HFLWC1
!
!     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
!
!     VHCPC=canopy heat capacity
!     TKCY=equilibrium canopy temperature for HFLXS
!
      VHCPC(NZ,NY,NX)=VHCPX+cpw*(EVAPC(NZ,NY,NX)+FLWC(NZ,NY,NX))
      TKCY=(TKCX*VHCPX+TKA(NY,NX)*PARHC+HFLXS)/(VHCPC(NZ,NY,NX)+PARHC)
      TKCY=AMIN1(TKA(NY,NX)+10.0,AMAX1(TKA(NY,NX)-10.0,TKCY))
!
!     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
!
!     XC,IC=magnitude,direction of change in canopy temp for next cycle
!
      IF((IC.EQ.0.AND.TKCY.GT.TKC1).OR.(IC.EQ.1.AND.TKCY.LT.TKC1))THEN
      XC=0.5*XC
      ENDIF
      TKCZ(NZ,NY,NX)=TKC1+0.1*(TKCY-TKC1)
      IF(TKCY.GT.TKC1)THEN
      IC=1
      ELSE
      IC=0
      ENDIF
!     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
!     WRITE(*,4444)'TKZ',I,J,NX,NY,NZ,NN,XC,TKC1,TKCY,TKCZ(NZ,NY,NX)
!    2,TKA(NY,NX),TKCX,VHCPX,PARHC,HFLXS,VHCPC(NZ,NY,NX),WVPLT,EX
!    2,FLWC(NZ,NY,NX),VOLWC(NZ,NY,NX),VOLWP(NZ,NY,NX),EVAPC(NZ,NY,NX)
!    2,RAD1(NZ,NY,NX),EFLXC(NZ,NY,NX),RA(NZ,NY,NX),RC(NZ,NY,NX)
!    2,EP(NZ,NY,NX),HFLXS,VFLXC,HFLWC1,RADC(NZ,NY,NX),FRADP(NZ,NY,NX)
!    3,THS(NY,NX),THRMGX(NY,NX)
!    2,RSMN(NZ,NY,NX),CCPOLT,OSWT,CCPOLP(NZ,NY,NX),CPOOLP(NZ,NY,NX)
!    4,DCO2(NZ,NY,NX),AREA(3,NU(NY,NX),NY,NX),WTLS(NZ,NY,NX)
!    2,PSILT(NZ,NY,NX),PSILG(NZ,NY,NX),RACZ(NZ,NY,NX),RAZ(NZ,NY,NX),RI
!    3,RIB(NY,NX),RA1,ARLFV(1,NZ,NY,NX),ARSTV(1,NZ,NY,NX)
!    4,WVPLT,VSTK
!4444  FORMAT(A8,6I4,60E16.8)
!     ENDIF
!
!     IF CONVERGENCE CRITERION IS MET OR ON EVERY TENTH ITERATION,
!     PROCEED TO WATER BALANCE
!
!     PSILC=canopy water potential adjusted for canopy height
!
      IF(ABS(TKCY-TKC1).LT.0.05.OR.(NN/10)*10.EQ.NN)THEN
      UPRT=0.0_r8
      PSILC=PSILT(NZ,NY,NX)-PSILH
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
      DO 4200 N=1,MY(NZ,NY,NX)
      DO  L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
      UPWTR(N,L,NZ,NY,NX)=AMAX1(AMIN1(0.0,-VOLWU(L)*FPQ(N,L,NZ)) &
      ,AMIN1((PSILC-PSIST1(L))/RSRS(N,L),VOLPU(L)*FPQ(N,L,NZ)))
      IF(UPWTR(N,L,NZ,NY,NX).GT.0.0)THEN
      UPWTR(N,L,NZ,NY,NX)=0.1*UPWTR(N,L,NZ,NY,NX)
      ENDIF
      UPRT=UPRT+UPWTR(N,L,NZ,NY,NX)
      ELSE
      UPWTR(N,L,NZ,NY,NX)=0.0_r8
      ENDIF
!     IF(NZ.EQ.2)THEN
!     WRITE(*,6565)'UPRT',I,J,NX,NY,NZ,L,N,NN,ILYR(N,L)
!    2,UPRT,UPWTR(N,L,NZ,NY,NX)
!    2,PSILC,PSIST1(L),PSISM(L,NY,NX),RSRS(N,L),RSSX(N,L)
!    2,RSRT(N,L),RSRG(N,L),RSR1(N,L),RSR2(N,L),PSILH,RTAR2
!    3,RSRR(N,NZ,NY,NX),VOLA(L,NY,NX),VOLWM(NPH,L,NY,NX)
!    4,VOLWU(L),THETY(L,NY,NX),VOLY(L,NY,NX),FPQ(N,L,NZ)
!6565  FORMAT(A8,9I4,30E12.4)
!     ENDIF
      enddo
4200  CONTINUE
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

      VOLWPZ=1.0E-06*WVPLT/FDMP
      DIFFZ=VOLWPZ-VOLWP(NZ,NY,NX)
      DIFFU=EP(NZ,NY,NX)-UPRT
      IF(test_aneb(UPRT,0.0_r8))THEN
      DIFF=ABS((DIFFU-DIFFZ)/UPRT)
      ELSE
      DIFF=ABS((DIFFU-DIFFZ)/VOLWPZ)
      ENDIF
      IF(DIFF.LT.5.0E-03)GO TO 4250
      IF(ABS(VOLWPZ-VOLWPX).GT.ZEROP(NZ,NY,NX))THEN
      RSSZ=ABS((PSILT(NZ,NY,NX)-PSIL2)/(VOLWPZ-VOLWPX))
      ELSEIF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSZ=1.0/CNDT
      ELSE
      RSSZ=ZEROL(NZ,NY,NX)
      ENDIF
      IF(ABS(EP(NZ,NY,NX)-EPX).GT.ZEROP(NZ,NY,NX))THEN
      RSSUX=ABS((PSILT(NZ,NY,NX)-PSIL2)/(EP(NZ,NY,NX)-EPX))
      IF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSU=AMIN1(1.0/CNDT,RSSUX)
      ELSE
      RSSU=RSSUX
      ENDIF
      ELSEIF(ABS(UPRT-UPRTX).GT.ZEROP(NZ,NY,NX))THEN
      RSSUX=ABS((PSILT(NZ,NY,NX)-PSIL2)/(UPRT-UPRTX))
      IF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSU=AMIN1(1.0/CNDT,RSSUX)
      ELSE
      RSSU=RSSUX
      ENDIF
      ELSEIF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSU=1.0/CNDT
      ELSE
      RSSU=ZEROL(NZ,NY,NX)
      ENDIF
!
!     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
!     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
!     WATER STORAGE
!
!     DPSI=change in PSILT for next convergence cycle
!     1.0E-03=acceptance criterion for DPSI
!
      DPSI=AMIN1(AMIN1(RSSZ,RSSU)*(DIFFU-DIFFZ),ABS(PSILT(NZ,NY,NX)))
!     IF(NX.EQ.1.AND.NY.EQ.5.AND.NZ.EQ.3)THEN
!     WRITE(*,2222)'PSI',I,J,NX,NY,NZ,NN,PSILT(NZ,NY,NX),PSIL2,DPSI
!    2,RSSUX,RSSU,RSSZ,1.0/CNDT,UPRT,UPRTX,EP(NZ,NY,NX),EPX,EX
!    3,EVAPC(NZ,NY,NX),RC(NZ,NY,NX),RA(NZ,NY,NX),FRADP(NZ,NY,NX)
!    3,PAREC,VPA(NY,NX),VPC,TKA(NY,NX),TKC1,VOLWP(NZ,NY,NX)
!    4,VOLWPZ,VOLWPX,WVPLT,DIFF,WFNC,PSILG(NZ,NY,NX)
!    5,FDMP,CCPOLT,OSWT,RAZ(NZ,NY,NX),RI,RIB(NY,NX)
!    5,((UPWTR(N,L,NZ,NY,NX),L=1,8),N=1,1)
!    6,((RSRS(N,L),L=1,8),N=1,1)
!    7,(PSIST1(L),L=1,8)
!    4,DIFFZ,DIFFU,DIFF
!2222  FORMAT(A8,6I4,80E12.4)
!     ENDIF
!
!     IF CONVERGENCE CRITERION IS MET THEN FINISH,
!     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
!     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
!
      IF((NN.GE.30.AND.ABS(DPSI).LT.1.0E-03).OR.NN.GE.MXN)GO TO 4250
      PSIL2=PSILT(NZ,NY,NX)
      EPX=EP(NZ,NY,NX)
      UPRTX=UPRT
      VOLWPX=VOLWPZ
      PSILT(NZ,NY,NX)=AMIN1(0.0,PSILT(NZ,NY,NX)+0.5*DPSI)
      XC=0.50
      GO TO 4000
!
!     RESET MIN STOMATAL RESISTANCE IN STOMATE.F BEFORE FINAL ITERATION
!
4250  IF(ICHK.EQ.1)THEN
      GO TO 4500
      ELSE
      ICHK=1
      CALL STOMATE(I,J,NZ,NY,NX)
      ENDIF
      ENDIF
4000  CONTINUE
4500  CONTINUE
      return
      end function CanopyEnergyH2OIteration
!------------------------------------------------------------------------
  subroutine CalcResistance(NZ,NY,NX,PATH,RRADL,RTARR,&
    RSRT,RSRG,RSR1,RSR2,RSSX,RSRS,CNDT,PSILH,ILYR)

  implicit none
  integer, intent(in)   :: NZ, NY, NX
  real(r8), intent(in)  :: PATH(2,JZ),RRADL(2,JZ),RTARR(2,JZ)
  real(r8), intent(out) :: RSRT(2,JZ),RSRG(2,JZ)
  real(r8), intent(out) :: RSR1(2,JZ),RSR2(2,JZ)
  real(r8), intent(out) :: RSSX(2,JZ),RSRS(2,JZ),CNDT
  real(r8), intent(out) :: PSILH
  integer, intent(out) :: ILYR(2,JZ)
  real(r8) :: FRADW,FRAD1,FRAD2
  real(r8) :: RSSL,RTAR2
  integer :: N, L

  !     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
  !
  !     HTSTZ=canopy height for water uptake
  !     PSILH=gravimetric water potential at HTSTZ
  !     FRADW=conducting elements of stalk relative to those of primary root
  !     PSILT=canopy total water potential
  !     EMODW=wood modulus of elasticity (MPa)
!
  CNDT=0.0_r8
  HTSTZ(NZ,NY,NX)=0.80*ZC(NZ,NY,NX)
  PSILH=-0.0098*HTSTZ(NZ,NY,NX)
  FRADW=1.0E+04*(AMAX1(0.5,1.0+PSILT(NZ,NY,NX)/EMODW))**4
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
  DO 3880 N=1,MY(NZ,NY,NX)
    DO  L=NU(NY,NX),NI(NZ,NY,NX)
      !     IF(NZ.EQ.2)THEN
      !     WRITE(*,2124)'ILYR',I,J,NX,NY,NZ,L,N,RTDNP(N,L,NZ,NY,NX)
      !    2,CNDU(L,NY,NX),RTN1(1,L,NZ,NY,NX),RTNL(N,L,NZ,NY,NX)
      !    3,THETW(L,NY,NX),ZEROP(NZ,NY,NX)
!2124  FORMAT(A8,7I4,20E12.4)
!      ENDIF
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX) &
        .AND.VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX) &
        .AND.RTDNP(N,L,NZ,NY,NX).GT.ZERO &
        .AND.CNDU(L,NY,NX).GT.ZERO &
        .AND.RTN1(1,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.RTNL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.THETW(L,NY,NX).GT.ZERO)THEN
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
        RSSL=(LOG(PATH(N,L)/RRADL(N,L))/RTARR(N,L))/PP(NZ,NY,NX)
        RSSX(N,L)=RSSL/CNDU(L,NY,NX)
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
        RTAR2=6.283*RRAD2(N,L,NZ,NY,NX)*RTLGP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
        RSRG(N,L)=RSRR(N,NZ,NY,NX)/RTAR2*VOLA(L,NY,NX)/VOLWM(NPH,L,NY,NX)
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
        FRAD1=(RRAD1(N,L,NZ,NY,NX)/RRAD2M(N,NZ,NY,NX))**4
        RSR1(N,L)=RSRA(N,NZ,NY,NX)*DPTHZ(L,NY,NX)/(FRAD1*RTN1(1,L,NZ,NY,NX)) &
          +RSRA(1,NZ,NY,NX)*HTSTZ(NZ,NY,NX)/(FRADW*RTN1(1,L,NZ,NY,NX))
        FRAD2=(RRAD2(N,L,NZ,NY,NX)/RRAD2M(N,NZ,NY,NX))**4
        RSR2(N,L)=RSRA(N,NZ,NY,NX)*RTLGA(N,L,NZ,NY,NX)/(FRAD2*RTNL(N,L,NZ,NY,NX))
      ELSE
        ILYR(N,L)=0
      ENDIF
    enddo
3880  CONTINUE

  DO 3890 N=1,MY(NZ,NY,NX)
    DO  L=NU(NY,NX),NI(NZ,NY,NX)
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
      !     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
      !     WRITE(*,8855)'RSRT',I,J,NX,NY,NZ,L,N,RSRT(N,L),RSRG(N,L)
      !    2,RSR1(N,L),RSR2(N,L),RSSX(N,L),RSRS(N,L),DPTHZ(L,NY,NX)
      !    3,HTSTZ(NZ,NY,NX),RSRA(1,NZ,NY,NX)*HTSTZ(NZ,NY,NX)
      !    4/(FRADW*AMAX1(PP(NZ,NY,NX),RTN1(1,L,NZ,NY,NX)))
      !    4,RTNL(N,L,NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)
      !    7,RTLGA(N,L,NZ,NY,NX),FRAD1,PP(NZ,NY,NX)
      !    8,RTN1(1,L,NZ,NY,NX),FRADW,RTNL(N,L,NZ,NY,NX),CNDT
    !8855  FORMAT(A8,7I4,30E14.6)
      !     ENDIF
      ENDIF
    enddo
3890  CONTINUE
  end subroutine CalcResistance
!------------------------------------------------------------------------

  subroutine HandleBareSoil(NZ,NY,NX,PSIST1,FDMP)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: PSIST1(JZ)
  real(r8), intent(out):: FDMP
  integer :: N,L
  real(r8) :: APSILT,APSIRT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM,FDMR
  real(r8) :: OSWT,WFNC

! begin_execution
  RAD1(NZ,NY,NX)=0.0_r8
  EFLXC(NZ,NY,NX)=0.0_r8
  SFLXC(NZ,NY,NX)=0.0_r8
  HFLXC(NZ,NY,NX)=0.0_r8
  EVAPC(NZ,NY,NX)=0.0_r8
  EP(NZ,NY,NX)=0.0_r8
  IF(ZC(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO)THEN
    TKC(NZ,NY,NX)=TKA(NY,NX)
  ELSE
    TKC(NZ,NY,NX)=TKW(1,NY,NX)
  ENDIF
  TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-TC2K
  FTHRM=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
  THRM1(NZ,NY,NX)=FTHRM*TKC(NZ,NY,NX)**4
  PSILT(NZ,NY,NX)=PSIST1(NG(NZ,NY,NX))
  APSILT=ABS(PSILT(NZ,NY,NX))
  FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
  CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
  OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
  PSILO(NZ,NY,NX)=FDMP/0.16*OSMO(NZ,NY,NX) &
    -8.3143*TKC(NZ,NY,NX)*FDMP*CCPOLT/OSWT
  PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
  WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
  RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)-RSMN(NZ,NY,NX))*WFNC
  RA(NZ,NY,NX)=RAZ(NZ,NY,NX)
  VHCPC(NZ,NY,NX)=cpw*(WTSHT(NZ,NY,NX)*10.0E-06)
  DTKC(NZ,NY,NX)=0.0_r8
  DO 4300 N=1,MY(NZ,NY,NX)
    DO  L=NU(NY,NX),NI(NZ,NY,NX)
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX) &
        -8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)-PSIRO(N,L,NZ,NY,NX))
      UPWTR(N,L,NZ,NY,NX)=0.0_r8
    enddo
4300  CONTINUE
  end subroutine HandleBareSoil
!------------------------------------------------------------------------

  subroutine UpdateCanopyWater(NZ,NY,NX,PARHC,PSIST1,RSRT,RSSX,&
    RSRS,TKCX,VHCPX,HFLWC1,UPRT,VFLXC,ILYR)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: PARHC,PSIST1(JZ),RSRT(2,JZ)
  real(r8), intent(in) :: RSSX(2,JZ),RSRS(2,JZ)
  real(r8), intent(in) :: TKCX,VHCPX,HFLWC1,UPRT,VFLXC
  integer , intent(in) :: ILYR(2,JZ)
  real(r8) :: APSIRT
  real(r8) :: CCPOLT
  real(r8) :: FDMR
  real(r8) :: OSWT
  integer :: N,L
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
  VOLWP(NZ,NY,NX)=VOLWP(NZ,NY,NX)+EP(NZ,NY,NX)-UPRT
  VOLWC(NZ,NY,NX)=VOLWC(NZ,NY,NX)+FLWC(NZ,NY,NX)+EVAPC(NZ,NY,NX)
  SFLXC(NZ,NY,NX)=PARHC*(TKA(NY,NX)-TKCZ(NZ,NY,NX))
  HFLXC(NZ,NY,NX)=TKCX*VHCPX-TKCZ(NZ,NY,NX)*VHCPC(NZ,NY,NX)+VFLXC+HFLWC1
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
  DO 4505 N=1,MY(NZ,NY,NX)
    DO 4510 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
        PSIRT(N,L,NZ,NY,NX)=AMIN1(0.0,(PSIST1(L)*RSRT(N,L) &
          +PSILT(NZ,NY,NX)*RSSX(N,L))/RSRS(N,L))
        APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
        FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
        CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)+CPPOLR(N,L,NZ,NY,NX)
        OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
        PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX) &
          -8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
        PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)-PSIRO(N,L,NZ,NY,NX))
      ELSE
        PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
        APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
        FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
        CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)+CPPOLR(N,L,NZ,NY,NX)
        OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
        PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX) &
          -8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
        PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)-PSIRO(N,L,NZ,NY,NX))
      ENDIF
!     IF(I.EQ.284)THEN
!     WRITE(*,1256)'PSIRT',I,J,NX,NY,NZ,NN,PSIRT(N,L,NZ,NY,NX)
!    2,PSIST1(L),RSRT(N,L),PSILT(NZ,NY,NX),RSSX(N,L),RSRS(N,L)
!    3,RSRG(N,L),RSR1(N,L),RSR2(N,L),RTAR2,VOLWM(NPH,L,NY,NX)
!1256  FORMAT(A8,6I4,20E12.4)
!     ENDIF
4510  CONTINUE
4505  CONTINUE
  end subroutine UpdateCanopyWater
!------------------------------------------------------------------------

  subroutine SetCanopyGrowthFuncs(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8) :: ACTV,RTK,STK,TKGO,TKSO
  integer :: L
  !
  !     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
  !     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
  !
  IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
    TKG(NZ,NY,NX)=TKS(NU(NY,NX),NY,NX)
    !     ELSEIF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)
    !    2.AND.IDAY(2,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
    !     TKG(NZ,NY,NX)=TKS(NU(NY,NX),NY,NX)
  ELSE
    TKG(NZ,NY,NX)=TKC(NZ,NY,NX)
  ENDIF
  TCG(NZ,NY,NX)=TKG(NZ,NY,NX)-TC2K
  !
  !     ARRHENIUS FUNCTION FOR CANOPY AND ROOT GROWTH WITH OFFSET
  !     FOR ZONE OF THERMAL ADAPTATION ENTERED IN 'READQ'
  !
  !     TKG,TKGO=canopy temperature,canopy temp used in Arrhenius eqn
  !     TKS,TKSO=soil temperature,soil temp used in Arrhenius eqn
  !     OFFST=shift in Arrhenius curve for thermal adaptation
  !     TFN3,TFN4=temperature function for canopy,root growth (25 oC =1)
  !     8.3143,710.0=gas constant,enthalpy
  !     62500,197500,222500=energy of activn,high,low temp inactivn(KJ mol-1)
  !     PSILZ=minimum daily canopy water potential
  !
  TKGO=TKG(NZ,NY,NX)+OFFST(NZ,NY,NX)
  RTK=8.3143*TKGO
  STK=710.0*TKGO
  ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
  TFN3(NZ,NY,NX)=EXP(25.229-62500/RTK)/ACTV
  DO 100 L=NU(NY,NX),NI(NZ,NY,NX)
    TKSO=TKS(L,NY,NX)+OFFST(NZ,NY,NX)
    RTK=8.3143*TKSO
    STK=710.0*TKSO
    ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
    TFN4(L,NZ,NY,NX)=EXP(25.229-62500/RTK)/ACTV
100   CONTINUE
  PSILZ(NZ,NY,NX)=AMIN1(PSILZ(NZ,NY,NX),PSILT(NZ,NY,NX))
  !
  !     DIURNAL CHILLING
  !
  !     CTC=chilling temperature from PFT file
  !     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
  !
  IF(TCC(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
    CHILL(NZ,NY,NX)=AMIN1(24.0,CHILL(NZ,NY,NX)+1.0)
  ELSE
    CHILL(NZ,NY,NX)=AMAX1(0.0,CHILL(NZ,NY,NX)-1.0)
  ENDIF
  end subroutine SetCanopyGrowthFuncs
!------------------------------------------------------------------------

  subroutine CanopyNH3Flux(NZ,NY,NX,FDMP)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in):: FDMP
  real(r8) :: FNH3P
  real(r8) :: CNH3P
  real(r8) :: SNH3P
  real(r8) :: ZPOOLB
  integer :: NB
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
  SNH3P=SNH3X*EXP(0.513-0.0171*TCC(NZ,NY,NX))
  FNH3P=1.0E-04*FDMP
  DO 105 NB=1,NBR(NZ,NY,NX)
    IF(WTLSB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.ARLFB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
      .AND.ARLFP(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CNH3P=AMAX1(0.0,FNH3P*CZPOLB(NB,NZ,NY,NX)/SNH3P)
      ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      RNH3B(NB,NZ,NY,NX)=AMIN1(0.1*ZPOOLB,AMAX1((CNH3E(NY,NX)-CNH3P)/(RA(NZ,NY,NX)+RC(NZ,NY,NX)) &
      *FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX) &
      *ARLFB(NB,NZ,NY,NX)/ARLFP(NZ,NY,NX),-0.1*ZPOOLB))
      ELSE
      RNH3B(NB,NZ,NY,NX)=0.0_r8
      ENDIF

105   CONTINUE
  end subroutine CanopyNH3Flux
!------------------------------------------------------------------------

  subroutine ZeroUptake(N,L,NZ,NY,NX)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ,NY,NX

  integer :: K
  !     begin_execution

  RCOFLA(N,L,NZ,NY,NX)=0.0_r8
  ROXFLA(N,L,NZ,NY,NX)=0.0_r8
  RCHFLA(N,L,NZ,NY,NX)=0.0_r8
  RN2FLA(N,L,NZ,NY,NX)=0.0_r8
  RNHFLA(N,L,NZ,NY,NX)=0.0_r8
  RCODFA(N,L,NZ,NY,NX)=0.0_r8
  ROXDFA(N,L,NZ,NY,NX)=0.0_r8
  RCHDFA(N,L,NZ,NY,NX)=0.0_r8
  RN2DFA(N,L,NZ,NY,NX)=0.0_r8
  RNHDFA(N,L,NZ,NY,NX)=0.0_r8
  RCO2S(N,L,NZ,NY,NX)=0.0_r8
  RUPOXS(N,L,NZ,NY,NX)=0.0_r8
  RUPCHS(N,L,NZ,NY,NX)=0.0_r8
  RUPN2S(N,L,NZ,NY,NX)=0.0_r8
  RUPN3S(N,L,NZ,NY,NX)=0.0_r8
  RCO2P(N,L,NZ,NY,NX)=0.0_r8
  RUPOXP(N,L,NZ,NY,NX)=0.0_r8
  DO 395 K=0,4
    RDFOMC(N,K,L,NZ,NY,NX)=0.0_r8
    RDFOMN(N,K,L,NZ,NY,NX)=0.0_r8
    RDFOMP(N,K,L,NZ,NY,NX)=0.0_r8
395   CONTINUE
  WFR(N,L,NZ,NY,NX)=1.0
  RUNNHP(N,L,NZ,NY,NX)=0.0_r8
  RUPNH4(N,L,NZ,NY,NX)=0.0_r8
  RUONH4(N,L,NZ,NY,NX)=0.0_r8
  RUCNH4(N,L,NZ,NY,NX)=0.0_r8
  RUNNBP(N,L,NZ,NY,NX)=0.0_r8
  RUPNHB(N,L,NZ,NY,NX)=0.0_r8
  RUONHB(N,L,NZ,NY,NX)=0.0_r8
  RUCNHB(N,L,NZ,NY,NX)=0.0_r8
  RUNNOP(N,L,NZ,NY,NX)=0.0_r8
  RUPNO3(N,L,NZ,NY,NX)=0.0_r8
  RUONO3(N,L,NZ,NY,NX)=0.0_r8
  RUCNO3(N,L,NZ,NY,NX)=0.0_r8
  RUNNXP(N,L,NZ,NY,NX)=0.0_r8
  RUPNOB(N,L,NZ,NY,NX)=0.0_r8
  RUONOB(N,L,NZ,NY,NX)=0.0_r8
  RUCNOB(N,L,NZ,NY,NX)=0.0_r8
  RUPP2P(N,L,NZ,NY,NX)=0.0_r8
  RUPH2P(N,L,NZ,NY,NX)=0.0_r8
  RUOH2P(N,L,NZ,NY,NX)=0.0_r8
  RUCH2P(N,L,NZ,NY,NX)=0.0_r8
  RUPP2B(N,L,NZ,NY,NX)=0.0_r8
  RUPH2B(N,L,NZ,NY,NX)=0.0_r8
  RUOH2B(N,L,NZ,NY,NX)=0.0_r8
  RUCH2B(N,L,NZ,NY,NX)=0.0_r8
  RUPP1P(N,L,NZ,NY,NX)=0.0_r8
  RUPH1P(N,L,NZ,NY,NX)=0.0_r8
  RUOH1P(N,L,NZ,NY,NX)=0.0_r8
  RUCH1P(N,L,NZ,NY,NX)=0.0_r8
  RUPP1B(N,L,NZ,NY,NX)=0.0_r8
  RUPH1B(N,L,NZ,NY,NX)=0.0_r8
  RUOH1B(N,L,NZ,NY,NX)=0.0_r8
  RUCH1B(N,L,NZ,NY,NX)=0.0_r8
  IF(N.EQ.1)RUPNF(L,NZ,NY,NX)=0.0_r8
  end subroutine ZeroUptake
!------------------------------------------------------------------------
  subroutine NoActiveNutrientUptake(N,L,NZ,NY,NX)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ,NY,NX
!     begin_execution

  RUNNHP(N,L,NZ,NY,NX)=0.0_r8
  RUPNH4(N,L,NZ,NY,NX)=0.0_r8
  RUONH4(N,L,NZ,NY,NX)=0.0_r8
  RUCNH4(N,L,NZ,NY,NX)=0.0_r8
  RUNNBP(N,L,NZ,NY,NX)=0.0_r8
  RUPNHB(N,L,NZ,NY,NX)=0.0_r8
  RUONHB(N,L,NZ,NY,NX)=0.0_r8
  RUCNHB(N,L,NZ,NY,NX)=0.0_r8
  RUNNOP(N,L,NZ,NY,NX)=0.0_r8
  RUPNO3(N,L,NZ,NY,NX)=0.0_r8
  RUONO3(N,L,NZ,NY,NX)=0.0_r8
  RUCNO3(N,L,NZ,NY,NX)=0.0_r8
  RUNNXP(N,L,NZ,NY,NX)=0.0_r8
  RUPNOB(N,L,NZ,NY,NX)=0.0_r8
  RUONOB(N,L,NZ,NY,NX)=0.0_r8
  RUCNOB(N,L,NZ,NY,NX)=0.0_r8
  RUPP2P(N,L,NZ,NY,NX)=0.0_r8
  RUPH2P(N,L,NZ,NY,NX)=0.0_r8
  RUOH2P(N,L,NZ,NY,NX)=0.0_r8
  RUCH2P(N,L,NZ,NY,NX)=0.0_r8
  RUPP2B(N,L,NZ,NY,NX)=0.0_r8
  RUPH2B(N,L,NZ,NY,NX)=0.0_r8
  RUOH2B(N,L,NZ,NY,NX)=0.0_r8
  RUCH2B(N,L,NZ,NY,NX)=0.0_r8
  RUPP1P(N,L,NZ,NY,NX)=0.0_r8
  RUPH1P(N,L,NZ,NY,NX)=0.0_r8
  RUOH1P(N,L,NZ,NY,NX)=0.0_r8
  RUCH1P(N,L,NZ,NY,NX)=0.0_r8
  RUPP1B(N,L,NZ,NY,NX)=0.0_r8
  RUPH1B(N,L,NZ,NY,NX)=0.0_r8
  RUOH1B(N,L,NZ,NY,NX)=0.0_r8
  RUCH1B(N,L,NZ,NY,NX)=0.0_r8
  end subroutine NoActiveNutrientUptake
!------------------------------------------------------------------------
  subroutine UptakeHPO4(N,L,NZ,NY,NX,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ,NY,NX
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
    IF(VLPO4(L,NY,NX).GT.ZERO.AND.CH1P4(L,NY,NX) &
      .GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH1P=UPWTRP*CH1P4(L,NY,NX)*VLPO4(L,NY,NX)
      DIFH1P=DIFFL*VLPO4(L,NY,NX)
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
    UPMXP=0.1*UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLPO4(L,NY,NX)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFH1P+RMFH1P)*CH1P4(L,NY,NX)
    Y=DIFH1P*UPMNPO(N,NZ,NY,NX)
    B=-UPMX-DIFH1P*UPKMPO(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKH1P=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH1P*UPKMPO(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKHP1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H1POM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    H1POX=AMAX1(0.0,FP14X*(H1PO4(L,NY,NX)-H1POM))
    RUPP1P(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH1P*PP(NZ,NY,NX))
    RUPH1P(N,L,NZ,NY,NX)=AMIN1(H1POX,RUPP1P(N,L,NZ,NY,NX))
    RUOH1P(N,L,NZ,NY,NX)=AMIN1(H1POX,AMAX1(0.0,RTKHP1*PP(NZ,NY,NX)))
    RUCH1P(N,L,NZ,NY,NX)=RUPH1P(N,L,NZ,NY,NX)/FCUP
    !     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.3)THEN
    !     WRITE(*,2226)'UPPO4',I,J,NZ,L,N,RUPH2P(N,L,NZ,NY,NX),FPO4X
    !    2,H2PO4(L,NY,NX),RUPP2P(N,L,NZ,NY,NX),UPMX,DIFPO,UPKMPO(N,NZ,NY,NX)
    !    3,UPMNPO(N,NZ,NY,NX),RMFH2P,CH2P4(L,NY,NX),UPMXP,WFR(N,L,NZ,NY,NX)
    !    4,FCUP,FZUP,FPUP,UPMXPO(N,NZ,NY,NX),RTARP(N,L,NZ,NY,NX),FWSRT
    !    5,TFN4(L,NZ,NY,NX),DIFFL,FH2P,CPO4S(L,NY,NX),CPOOLR(N,L,NZ,NY,NX)
    !    6,PPOOLR(N,L,NZ,NY,NX),RTKH2P,PP(NZ,NY,NX)
    !    2,RTLGP(N,L,NZ,NY,NX)
!2226  FORMAT(A8,5I4,40E12.4)
    !     ENDIF
  ELSE
    RUPP1P(N,L,NZ,NY,NX)=0.0_r8
    RUPH1P(N,L,NZ,NY,NX)=0.0_r8
    RUOH1P(N,L,NZ,NY,NX)=0.0_r8
    RUCH1P(N,L,NZ,NY,NX)=0.0_r8
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
  IF(VLPOB(L,NY,NX).GT.ZERO.AND.CH1P4B(L,NY,NX) &
    .GT.UPMNPO(N,NZ,NY,NX))THEN
    RMFH2B=UPWTRP*CH1P4B(L,NY,NX)*VLPOB(L,NY,NX)
    DIFH1B=DIFFL*VLPOB(L,NY,NX)
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
    UPMXP=0.1*UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLPOB(L,NY,NX)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFH1B+RMFH2B)*CH1P4B(L,NY,NX)
    Y=DIFH1B*UPMNPO(N,NZ,NY,NX)
    B=-UPMX-DIFH1B*UPKMPO(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKH1B=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH1B*UPKMPO(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKHB1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H1PXM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    H1PXB=AMAX1(0.0,FP1BX*(H1POB(L,NY,NX)-H1PXM))
    RUPP1B(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH1B*PP(NZ,NY,NX))
    RUPH1B(N,L,NZ,NY,NX)=AMIN1(H1PXB,RUPP1B(N,L,NZ,NY,NX))
    RUOH1B(N,L,NZ,NY,NX)=AMIN1(H1PXB,AMAX1(0.0,RTKHB1*PP(NZ,NY,NX)))
    RUCH1B(N,L,NZ,NY,NX)=RUPH1B(N,L,NZ,NY,NX)/FCUP
  ELSE
    RUPP1B(N,L,NZ,NY,NX)=0.0_r8
    RUPH1B(N,L,NZ,NY,NX)=0.0_r8
    RUOH1B(N,L,NZ,NY,NX)=0.0_r8
    RUCH1B(N,L,NZ,NY,NX)=0.0_r8
  ENDIF
  end subroutine UptakeHPO4
!------------------------------------------------------------------------

  subroutine UptakeH2PO4(N,L,NZ,NY,NX,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer,  intent(in) :: N,L
  integer,  intent(in) :: NZ,NY,NX
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
  !     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH2P4=H2PO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
  !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
!
    IF(VLPO4(L,NY,NX).GT.ZERO.AND.CH2P4(L,NY,NX) &
      .GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2P=UPWTRP*CH2P4(L,NY,NX)*VLPO4(L,NY,NX)
      DIFH2P=DIFFL*VLPO4(L,NY,NX)
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
      UPMXP=UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
        *FWSRT*TFN4(L,NZ,NY,NX)*VLPO4(L,NY,NX)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
      X=(DIFH2P+RMFH2P)*CH2P4(L,NY,NX)
      Y=DIFH2P*UPMNPO(N,NZ,NY,NX)
      B=-UPMX-DIFH2P*UPKMPO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKH2P=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH2P*UPKMPO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKHPP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H2POM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
      H2POX=AMAX1(0.0,FPO4X*(H2PO4(L,NY,NX)-H2POM))
      RUPP2P(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH2P*PP(NZ,NY,NX))
      RUPH2P(N,L,NZ,NY,NX)=AMIN1(H2POX,RUPP2P(N,L,NZ,NY,NX))
      RUOH2P(N,L,NZ,NY,NX)=AMIN1(H2POX,AMAX1(0.0_r8,RTKHPP*PP(NZ,NY,NX)))
      RUCH2P(N,L,NZ,NY,NX)=RUPH2P(N,L,NZ,NY,NX)/FCUP
      !     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.3)THEN
      !     WRITE(*,2223)'UPPO4',I,J,NZ,L,N,RUPH2P(N,L,NZ,NY,NX),FPO4X
      !    2,H2PO4(L,NY,NX),RUPP2P(N,L,NZ,NY,NX),UPMX,DIFPO,UPKMPO(N,NZ,NY,NX)
      !    3,UPMNPO(N,NZ,NY,NX),RMFH2P,CH2P4(L,NY,NX),UPMXP,WFR(N,L,NZ,NY,NX)
      !    4,FCUP,FZUP,FPUP,UPMXPO(N,NZ,NY,NX),RTARP(N,L,NZ,NY,NX),FWSRT
      !    5,TFN4(L,NZ,NY,NX),DIFFL,CPO4S(L,NY,NX),CPOOLR(N,L,NZ,NY,NX)
      !    6,PPOOLR(N,L,NZ,NY,NX),RTKH2P,PP(NZ,NY,NX)
      !    2,RTLGP(N,L,NZ,NY,NX)
!2223  FORMAT(A8,5I4,40E12.4)
!     ENDIF
    ELSE
      RUPP2P(N,L,NZ,NY,NX)=0.0_r8
      RUPH2P(N,L,NZ,NY,NX)=0.0_r8
      RUOH2P(N,L,NZ,NY,NX)=0.0_r8
      RUCH2P(N,L,NZ,NY,NX)=0.0_r8
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

  IF(VLPOB(L,NY,NX).GT.ZERO.AND.CH2P4B(L,NY,NX) &
    .GT.UPMNPO(N,NZ,NY,NX))THEN
    RMFH2B=UPWTRP*CH2P4B(L,NY,NX)*VLPOB(L,NY,NX)
    DIFH2B=DIFFL*VLPOB(L,NY,NX)
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
    UPMXP=UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLPOB(L,NY,NX)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFH2B+RMFH2B)*CH2P4B(L,NY,NX)
    Y=DIFH2B*UPMNPO(N,NZ,NY,NX)
    B=-UPMX-DIFH2B*UPKMPO(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKH2B=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH2B*UPKMPO(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKHPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H2PXM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    H2PXB=AMAX1(0.0,FPOBX*(H2POB(L,NY,NX)-H2PXM))
    RUPP2B(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH2B*PP(NZ,NY,NX))
    RUPH2B(N,L,NZ,NY,NX)=AMIN1(H2PXB,RUPP2B(N,L,NZ,NY,NX))
    RUOH2B(N,L,NZ,NY,NX)=AMIN1(H2PXB,AMAX1(0.0,RTKHPB*PP(NZ,NY,NX)))
    RUCH2B(N,L,NZ,NY,NX)=RUPH2B(N,L,NZ,NY,NX)/FCUP
  ELSE
    RUPP2B(N,L,NZ,NY,NX)=0.0_r8
    RUPH2B(N,L,NZ,NY,NX)=0.0_r8
    RUOH2B(N,L,NZ,NY,NX)=0.0_r8
    RUCH2B(N,L,NZ,NY,NX)=0.0_r8
  ENDIF
  end subroutine UptakeH2PO4
!------------------------------------------------------------------------

  subroutine UptakeMineralPhosporhus(N,L,NZ,NY,NX,PATH,RRADL,FPQ,FPP,&
    RTARR,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in):: PATH(2,JZ),RRADL(2,JZ),FPQ(2,JZ,05),FPP(2,JZ,05)
  real(r8), intent(in):: RTARR(2,JZ),FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: TFPO4X,TFP14X,TFPOBX,TFP1BX
  real(r8) :: DIFFL
  real(r8) :: FP14X,FP1BX
  real(r8) :: FPO4X,FPOBX
  real(r8) :: POSGX
  real(r8) :: PATHL

  TFPO4X=0.0_r8
  TFPOBX=0.0_r8
  TFP14X=0.0_r8
  TFP1BX=0.0_r8

!     begin_execution
  IF(RPO4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FPO4X=AMAX1(FPP(N,L,NZ),RUPP2P(N,L,NZ,NY,NX)/RPO4Y(L,NY,NX))
  ELSE
    FPO4X=FPQ(N,L,NZ)
  ENDIF
  IF(RPOBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FPOBX=AMAX1(FPP(N,L,NZ),RUPP2B(N,L,NZ,NY,NX)/RPOBY(L,NY,NX))
  ELSE
    FPOBX=FPQ(N,L,NZ)
  ENDIF
  IF(RP14Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FP14X=AMAX1(FPP(N,L,NZ),RUPP1P(N,L,NZ,NY,NX)/RP14Y(L,NY,NX))
  ELSE
    FP14X=FPQ(N,L,NZ)
  ENDIF
  IF(RP1BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FP1BX=AMAX1(FPP(N,L,NZ),RUPP1B(N,L,NZ,NY,NX)/RP1BY(L,NY,NX))
  ELSE
    FP1BX=FPQ(N,L,NZ)
  ENDIF

  TFPO4X=TFPO4X+FPO4X
  TFP14X=TFP14X+FP14X

  TFPOBX=TFPOBX+FPOBX
  TFP1BX=TFP1BX+FP1BX

  IF(FPUP.GT.ZERO2)THEN
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
    POSGX=POSGL(L,NY,NX)*TORT(NPH,L,NY,NX)
    PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*POSGX))
    DIFFL=POSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))

    call UptakeH2PO4(N,L,NZ,NY,NX,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,UPWTRP)

    call UptakeHPO4(N,L,NZ,NY,NX,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,UPWTRP)

  ELSE
    RUPP2P(N,L,NZ,NY,NX)=0.0_r8
    RUPH2P(N,L,NZ,NY,NX)=0.0_r8
    RUOH2P(N,L,NZ,NY,NX)=0.0_r8
    RUCH2P(N,L,NZ,NY,NX)=0.0_r8
    RUPP2B(N,L,NZ,NY,NX)=0.0_r8
    RUPH2B(N,L,NZ,NY,NX)=0.0_r8
    RUOH2B(N,L,NZ,NY,NX)=0.0_r8
    RUCH2B(N,L,NZ,NY,NX)=0.0_r8
    RUPP1P(N,L,NZ,NY,NX)=0.0_r8
    RUPH1P(N,L,NZ,NY,NX)=0.0_r8
    RUOH1P(N,L,NZ,NY,NX)=0.0_r8
    RUCH1P(N,L,NZ,NY,NX)=0.0_r8
    RUPP1B(N,L,NZ,NY,NX)=0.0_r8
    RUPH1B(N,L,NZ,NY,NX)=0.0_r8
    RUOH1B(N,L,NZ,NY,NX)=0.0_r8
    RUCH1B(N,L,NZ,NY,NX)=0.0_r8
  ENDIF
  end subroutine UptakeMineralPhosporhus
!------------------------------------------------------------------------

  subroutine UptakeNO3(N,L,NZ,NY,NX,FNO3X,FNOBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in):: FNO3X,FNOBX,PATH(2,JZ),RRADL(2,JZ),RTARR(2,JZ)
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
  ZOSGX=ZOSGL(L,NY,NX)*TORT(NPH,L,NY,NX)
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
  IF(VLNO3(L,NY,NX).GT.ZERO.AND.CNO3S(L,NY,NX) &
    .GT.UPMNZO(N,NZ,NY,NX))THEN
    RMFNO3=UPWTRP*CNO3S(L,NY,NX)*VLNO3(L,NY,NX)
    DIFNO3=DIFFL*VLNO3(L,NY,NX)
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
    UPMXP=UPMXZO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLNO3(L,NY,NX)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFNO3+RMFNO3)*CNO3S(L,NY,NX)
    Y=DIFNO3*UPMNZO(N,NZ,NY,NX)
    B=-UPMX-DIFNO3*UPKMZO(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKNO3=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNO3*UPKMZO(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKNOP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNO3M=UPMNZO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNO3(L,NY,NX)
    ZNO3X=AMAX1(0.0,FNO3X*(ZNO3S(L,NY,NX)-ZNO3M))
    RUNNOP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNO3*PP(NZ,NY,NX))
    RUPNO3(N,L,NZ,NY,NX)=AMIN1(ZNO3X,RUNNOP(N,L,NZ,NY,NX))
    RUONO3(N,L,NZ,NY,NX)=AMIN1(ZNO3X,AMAX1(0.0 &
      ,RTKNOP*PP(NZ,NY,NX)))
    RUCNO3(N,L,NZ,NY,NX)=RUPNO3(N,L,NZ,NY,NX)/FCUP
    !     IF(NX.EQ.4.AND.NY.EQ.2)THEN
    !     WRITE(*,1111)'UPNO3',I,J,NZ,L,N,RUPNO3(N,L,NZ,NY,NX),FNO3X
    !    2,ZNO3S(L,NY,NX),ZNO3M,RTDNP(N,L,NZ,NY,NX),RTKNO3,RMFNO3,X,Y,B,C
    !    2,UPMX,CNO3S(L,NY,NX),DIFNO,RUONO3(N,L,NZ,NY,NX)
    !    3,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX)
    !    4,THETW(L,NY,NX),TKS(L,NY,NX),RSCS(L,NY,NX),UPMXP,FWSRT
    !    5,FZUP,FCUP,COXYS(L,NY,NX),COXYG(L,NY,NX),WFR(N,L,NZ,NY,NX)
    !    6,CCPOLP(NZ,NY,NX),CZPOLP(NZ,NY,NX),CPPOLP(NZ,NY,NX)
    !    7,FDBK(1,NZ,NY,NX),PSIST1(L),PSIRT(N,L,NZ,NY,NX)
    !    2,ZPOOLR(N,L,NZ,NY,NX),WTRTL(N,L,NZ,NY,NX)
    !    3,RUNNOP(N,L,NZ,NY,NX),RNO3Y(L,NY,NX)
!1111  FORMAT(A8,5I4,40E12.4)
!     ENDIF
  ELSE
    RUNNOP(N,L,NZ,NY,NX)=0.0_r8
    RUPNO3(N,L,NZ,NY,NX)=0.0_r8
    RUONO3(N,L,NZ,NY,NX)=0.0_r8
    RUCNO3(N,L,NZ,NY,NX)=0.0_r8
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

  IF(VLNOB(L,NY,NX).GT.ZERO.AND.CNO3B(L,NY,NX) &
    .GT.UPMNZO(N,NZ,NY,NX))THEN
    RMFNOB=UPWTRP*CNO3B(L,NY,NX)*VLNOB(L,NY,NX)
    DIFNOB=DIFFL*VLNOB(L,NY,NX)
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
    UPMXP=UPMXZO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLNOB(L,NY,NX)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFNOB+RMFNOB)*CNO3B(L,NY,NX)
    Y=DIFNOB*UPMNZO(N,NZ,NY,NX)
    B=-UPMX-DIFNOB*UPKMZO(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKNOB=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNOB*UPKMZO(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKNPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNOBM=UPMNZO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNOB(L,NY,NX)
    ZNOBX=AMAX1(0.0,FNOBX*(ZNO3B(L,NY,NX)-ZNOBM))
    RUNNXP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNOB*PP(NZ,NY,NX))
    RUPNOB(N,L,NZ,NY,NX)=AMIN1(ZNOBX,RUNNXP(N,L,NZ,NY,NX))
    RUONOB(N,L,NZ,NY,NX)=AMIN1(ZNOBX &
      ,AMAX1(0.0,RTKNPB*PP(NZ,NY,NX)))
    RUCNOB(N,L,NZ,NY,NX)=RUPNOB(N,L,NZ,NY,NX)/FCUP
  ELSE
    RUNNXP(N,L,NZ,NY,NX)=0.0_r8
    RUPNOB(N,L,NZ,NY,NX)=0.0_r8
    RUONOB(N,L,NZ,NY,NX)=0.0_r8
    RUCNOB(N,L,NZ,NY,NX)=0.0_r8
  ENDIF
  end subroutine UptakeNO3
!------------------------------------------------------------------------

  subroutine UptakeNH4(N,L,NZ,NY,NX,FNH4X,FNHBX,PATH,RRADL,RTARR,&
    FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: FNH4X,FNHBX,PATH(2,JZ),RRADL(2,JZ)
  real(r8), intent(in) :: RTARR(2,JZ)
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
! ZNSGL=NH4 diffusivity
! TORT=soil tortuosity
! PATH=path length of water and nutrient uptake
! RRADL=root radius
! DIFFL=NH4 diffusion per plant
!
  ZNSGX=ZNSGL(L,NY,NX)*TORT(NPH,L,NY,NX)
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
  IF(VLNH4(L,NY,NX).GT.ZERO.AND.CNH4S(L,NY,NX) &
    .GT.UPMNZH(N,NZ,NY,NX))THEN
    RMFNH4=UPWTRP*CNH4S(L,NY,NX)*VLNH4(L,NY,NX)
    DIFNH4=DIFFL*VLNH4(L,NY,NX)
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
    UPMXP=UPMXZH(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLNH4(L,NY,NX)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFNH4+RMFNH4)*CNH4S(L,NY,NX)
    Y=DIFNH4*UPMNZH(N,NZ,NY,NX)
    B=-UPMX-DIFNH4*UPKMZH(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKNH4=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNH4*UPKMZH(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKNHP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNH4M=UPMNZH(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNH4(L,NY,NX)
    ZNH4X=AMAX1(0.0,FNH4X*(ZNH4S(L,NY,NX)-ZNH4M))
    RUNNHP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNH4*PP(NZ,NY,NX))
    RUPNH4(N,L,NZ,NY,NX)=AMIN1(ZNH4X,RUNNHP(N,L,NZ,NY,NX))
    RUONH4(N,L,NZ,NY,NX)=AMIN1(ZNH4X,AMAX1(0.0 &
      ,RTKNHP*PP(NZ,NY,NX)))
    RUCNH4(N,L,NZ,NY,NX)=RUPNH4(N,L,NZ,NY,NX)/FCUP
!   IF(NX.EQ.1.OR.NZ.EQ.4)THEN
!     WRITE(*,1110)'UPNH4',I,J,NZ,L,N,RUNNHP(N,L,NZ,NY,NX)
!    2,RUPNH4(N,L,NZ,NY,NX),RTKNH4,RMFNH4,X,Y,B,C,UPMX,UPMXP
!    2,WFR(N,L,NZ,NY,NX),CNH4S(L,NY,NX),DIFNH,RTDNP(N,L,NZ,NY,NX)
!    2,WTRTD(N,L,NZ,NY,NX),CNH4S(L,NY,NX),RDFOMN(N,L,NZ,NY,NX)
!    3,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX)
!    4,THETW(L,NY,NX),TKS(L,NY,NX),RSCS(L,NY,NX),UPMXP,FWSRT
!    5,FZUP,FCUP,COXYS(L,NY,NX),COXYG(L,NY,NX)
!    6,CCPOLP(NZ,NY,NX)
!    7,CZPOLP(NZ,NY,NX),CPPOLP(NZ,NY,NX),FDBK(1,NZ,NY,NX),PSIST1(L)
!    2,PSIRT(N,L,NZ,NY,NX),ZPOOLR(N,L,NZ,NY,NX),WTRTL(N,L,NZ,NY,NX)
!    3,RTARP(N,L,NZ,NY,NX),RRADL(N,L),PATH(N,L)
!    4,DIFFL,ZNSGX,RTARR(N,L),PATHL,RRADL(N,L),VLNH4(L,NY,NX)
!1110  FORMAT(A8,5I4,100E24.16)
!     ENDIF
  ELSE
    RUNNHP(N,L,NZ,NY,NX)=0.0_r8
    RUPNH4(N,L,NZ,NY,NX)=0.0_r8
    RUONH4(N,L,NZ,NY,NX)=0.0_r8
    RUCNH4(N,L,NZ,NY,NX)=0.0_r8
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

  IF(VLNHB(L,NY,NX).GT.ZERO.AND.CNH4B(L,NY,NX) &
    .GT.UPMNZH(N,NZ,NY,NX))THEN
    RMFNHB=UPWTRP*CNH4B(L,NY,NX)*VLNHB(L,NY,NX)
    DIFNHB=DIFFL*VLNHB(L,NY,NX)
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
    UPMXP=UPMXZH(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) &
      *FWSRT*TFN4(L,NZ,NY,NX)*VLNHB(L,NY,NX)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
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
    X=(DIFNHB+RMFNHB)*CNH4B(L,NY,NX)
    Y=DIFNHB*UPMNZH(N,NZ,NY,NX)
    B=-UPMX-DIFNHB*UPKMZH(N,NZ,NY,NX)-X+Y
    C=(X-Y)*UPMX
    RTKNHB=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNHB*UPKMZH(N,NZ,NY,NX)-X+Y
    CP=(X-Y)*UPMXP
    RTKNBP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNHBM=UPMNZH(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNHB(L,NY,NX)
    ZNHBX=AMAX1(0.0,FNHBX*(ZNH4B(L,NY,NX)-ZNHBM))
    RUNNBP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNHB*PP(NZ,NY,NX))
    RUPNHB(N,L,NZ,NY,NX)=AMIN1(ZNHBX,RUNNBP(N,L,NZ,NY,NX))
    RUONHB(N,L,NZ,NY,NX)=AMIN1(ZNHBX,AMAX1(0.0 &
      ,RTKNBP*PP(NZ,NY,NX)))
    RUCNHB(N,L,NZ,NY,NX)=RUPNHB(N,L,NZ,NY,NX)/FCUP
  ELSE
    RUNNBP(N,L,NZ,NY,NX)=0.0_r8
    RUPNHB(N,L,NZ,NY,NX)=0.0_r8
    RUONHB(N,L,NZ,NY,NX)=0.0_r8
    RUCNHB(N,L,NZ,NY,NX)=0.0_r8
  ENDIF
  end subroutine UptakeNH4
!------------------------------------------------------------------------

  subroutine UptakeMineralNitrogen(N,L,NZ,NY,NX,PATH,RRADL,FPQ,&
    FPP,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: PATH(2,JZ),RRADL(2,JZ),FPQ(2,JZ,05)
  real(r8), intent(in) :: FPP(2,JZ,05),RTARR(2,JZ)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: FNO3X,FNOBX,FNH4X,FNHBX
  real(r8) :: TFNH4X,TFNO3X,TFNHBX,TFNOBX

!     begin_execution

  TFNH4X=0.0_r8
  TFNHBX=0.0_r8
  TFNO3X=0.0_r8
  TFNOBX=0.0_r8


  IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNH4X=AMAX1(FPP(N,L,NZ),RUNNHP(N,L,NZ,NY,NX)/RNH4Y(L,NY,NX))
  ELSE
    FNH4X=FPQ(N,L,NZ)
  ENDIF
  IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNHBX=AMAX1(FPP(N,L,NZ),RUNNBP(N,L,NZ,NY,NX)/RNHBY(L,NY,NX))
  ELSE
    FNHBX=FPQ(N,L,NZ)
  ENDIF

  IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO3X=AMAX1(FPP(N,L,NZ),RUNNOP(N,L,NZ,NY,NX)/RNO3Y(L,NY,NX))
  ELSE
    FNO3X=FPQ(N,L,NZ)
  ENDIF
  IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNOBX=AMAX1(FPP(N,L,NZ),RUNNXP(N,L,NZ,NY,NX)/RN3BY(L,NY,NX))
  ELSE
    FNOBX=FPQ(N,L,NZ)
  ENDIF

  TFNH4X=TFNH4X+FNH4X
  TFNO3X=TFNO3X+FNO3X
  TFNHBX=TFNHBX+FNHBX
  TFNOBX=TFNOBX+FNOBX

  IF(FZUP.GT.ZERO2)THEN
!
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
    !     FROM SOIL TO ROOT
    !
    call UptakeNH4(N,L,NZ,NY,NX,FNH4X,FNHBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

    call UptakeNO3(N,L,NZ,NY,NX,FNO3X,FNOBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  ELSE
    RUNNHP(N,L,NZ,NY,NX)=0.0_r8
    RUPNH4(N,L,NZ,NY,NX)=0.0_r8
    RUONH4(N,L,NZ,NY,NX)=0.0_r8
    RUCNH4(N,L,NZ,NY,NX)=0.0_r8
    RUNNBP(N,L,NZ,NY,NX)=0.0_r8
    RUPNHB(N,L,NZ,NY,NX)=0.0_r8
    RUONHB(N,L,NZ,NY,NX)=0.0_r8
    RUCNHB(N,L,NZ,NY,NX)=0.0_r8
    RUNNOP(N,L,NZ,NY,NX)=0.0_r8
    RUPNO3(N,L,NZ,NY,NX)=0.0_r8
    RUONO3(N,L,NZ,NY,NX)=0.0_r8
    RUCNO3(N,L,NZ,NY,NX)=0.0_r8
    RUNNXP(N,L,NZ,NY,NX)=0.0_r8
    RUPNOB(N,L,NZ,NY,NX)=0.0_r8
    RUONOB(N,L,NZ,NY,NX)=0.0_r8
    RUCNOB(N,L,NZ,NY,NX)=0.0_r8
  ENDIF
  end subroutine UptakeMineralNitrogen
!------------------------------------------------------------------------

  subroutine RootExudates(N,L,NZ,NY,NX)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ,NY,NX

  real(r8) :: CPOOLX,CPOOLT
  real(r8) :: PPOOLX,ZPOOLX
  real(r8) :: VOLWK,VOLWT
  real(r8) :: XFRC,XFRN,XFRP
  integer :: K
  !     begin_execution
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
  DO 195 K=0,4
    VOLWK=VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)
    IF(VOLWK.GT.ZEROS2(NY,NX) &
      .AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      VOLWT=VOLWK+RTVLW(N,L,NZ,NY,NX)
      CPOOLX=AMIN1(1.25E+03*RTVLW(N,L,NZ,NY,NX),CPOOLR(N,L,NZ,NY,NX))
      XFRC=(OQC(K,L,NY,NX)*RTVLW(N,L,NZ,NY,NX)-CPOOLX*VOLWK)/VOLWT
      RDFOMC(N,K,L,NZ,NY,NX)=FEXUC*XFRC
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX) &
        .AND.CPOOLR(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        CPOOLT=OQC(K,L,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
        ZPOOLX=0.1*ZPOOLR(N,L,NZ,NY,NX)
        PPOOLX=0.1*PPOOLR(N,L,NZ,NY,NX)
        XFRN=(OQN(K,L,NY,NX)*CPOOLR(N,L,NZ,NY,NX)-ZPOOLX*OQC(K,L,NY,NX))/CPOOLT
        XFRP=(OQP(K,L,NY,NX)*CPOOLR(N,L,NZ,NY,NX)-PPOOLX*OQC(K,L,NY,NX))/CPOOLT
        RDFOMN(N,K,L,NZ,NY,NX)=FEXUN*XFRN
        RDFOMP(N,K,L,NZ,NY,NX)=FEXUP*XFRP
      ELSE
        RDFOMN(N,K,L,NZ,NY,NX)=0.0_r8
        RDFOMP(N,K,L,NZ,NY,NX)=0.0_r8
      ENDIF
    ELSE
      RDFOMC(N,K,L,NZ,NY,NX)=0.0_r8
      RDFOMN(N,K,L,NZ,NY,NX)=0.0_r8
      RDFOMP(N,K,L,NZ,NY,NX)=0.0_r8
    ENDIF
    !     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.1)THEN
    !     WRITE(*,2224)'RDFOMC',I,J,NX,NY,L,NZ,K,N,RDFOMC(N,K,L,NZ,NY,NX)
    !    2,RDFOMN(N,K,L,NZ,NY,NX),RDFOMP(N,K,L,NZ,NY,NX)
    !    3,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
    !    2,CPOOLR(N,L,NZ,NY,NX),ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
    !    3,VOLWM(NPH,L,NY,NX),RTVLW(N,L,NZ,NY,NX),RTAR1X(N,NZ,NY,NX)
    !    4,RTAR2X(N,NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
    !    4,WTRTD(N,L,NZ,NY,NX)
    !    5,VOLWK,VOLWM(NPH,L,NY,NX),FOSRH(K,L,NY,NX)
    !    5,OQC(K,L,NY,NX)/VOLWK
    !    5,OQN(K,L,NY,NX)/OQC(K,L,NY,NX)
    !    5,OQP(K,L,NY,NX)/OQC(K,L,NY,NX)
    !    6,CPOOLR(N,L,NZ,NY,NX)/RTVLW(N,L,NZ,NY,NX)
    !    6,ZPOOLX/CPOOLR(N,L,NZ,NY,NX)
    !    6,PPOOLX/CPOOLR(N,L,NZ,NY,NX)
!2224  FORMAT(A8,8I4,30E12.4)
!     ENDIF
195   CONTINUE
  end subroutine RootExudates
!------------------------------------------------------------------------

  subroutine SumupNutrientUptake(N,L,NZ,NY,NX)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ,NY,NX

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
  DO 295 K=0,jcplx1
    UPOMC(NZ,NY,NX)=UPOMC(NZ,NY,NX)+RDFOMC(N,K,L,NZ,NY,NX)
    UPOMN(NZ,NY,NX)=UPOMN(NZ,NY,NX)+RDFOMN(N,K,L,NZ,NY,NX)
    UPOMP(NZ,NY,NX)=UPOMP(NZ,NY,NX)+RDFOMP(N,K,L,NZ,NY,NX)
    XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-RDFOMC(N,K,L,NZ,NY,NX)
    XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-RDFOMN(N,K,L,NZ,NY,NX)
    XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-RDFOMP(N,K,L,NZ,NY,NX)
295   CONTINUE
  UPNH4(NZ,NY,NX)=UPNH4(NZ,NY,NX)+RUPNH4(N,L,NZ,NY,NX)+RUPNHB(N,L,NZ,NY,NX)
  UPNO3(NZ,NY,NX)=UPNO3(NZ,NY,NX)+RUPNO3(N,L,NZ,NY,NX)+RUPNOB(N,L,NZ,NY,NX)
  UPH2P(NZ,NY,NX)=UPH2P(NZ,NY,NX)+RUPH2P(N,L,NZ,NY,NX)+RUPH2B(N,L,NZ,NY,NX)
  UPH1P(NZ,NY,NX)=UPH1P(NZ,NY,NX)+RUPH1P(N,L,NZ,NY,NX)+RUPH1B(N,L,NZ,NY,NX)
!     IF(J.EQ.12)THEN
!     WRITE(*,8765)'PLANT',I,J,NX,NY,L,NZ,N,TFOXYX,TFNH4X
!    2,TFNO3X,TFPO4X,TFNHBX,TFNOBX,TFPOBX
!8765  FORMAT(A8,7I4,7F15.6)
!     ENDIF
  end subroutine SumupNutrientUptake
!------------------------------------------------------------------------

  subroutine GetUptakeCapcity(N,L,NZ,NY,NX,FPQ,FPP,FCUP,FZUP,FPUP,&
    FWSRT,UPWTRP,UPWTRH,FOXYX)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ,NY,NX
  REAL(R8), INTENT(IN):: FPQ(2,JZ,05),FPP(2,JZ,05)
  real(r8), intent(out):: FCUP,FZUP,FPUP,FWSRT,UPWTRP,UPWTRH,FOXYX
  !
  !     UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL
  !     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
  !     PARAMETERS ARE DEFINED
  !
  !     CWSRTL,CWSRT=current,maximum protein concentration
  !     WSRTL,WTRTL=protein content,mass
  !     FWSRT=protein concentration relative to 5%
  !
  IF(WTRTL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
    CWSRTL(N,L,NZ,NY,NX)=AMIN1(CWSRT(NZ,NY,NX) &
      ,WSRTL(N,L,NZ,NY,NX)/WTRTL(N,L,NZ,NY,NX))
    FWSRT=CWSRTL(N,L,NZ,NY,NX)/0.05
  ELSE
    CWSRTL(N,L,NZ,NY,NX)=CWSRT(NZ,NY,NX)
    FWSRT=1.0
  ENDIF
  !
  !     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
  !
  !     RCO2N=total respiration from CPOOLR
  !     FCUP=limitation to active uptake respiration from CPOOLR
  !     CPOOLR=nonstructural C content
  !
  IF(RCO2N(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
    FCUP=AMAX1(0.0,AMIN1(1.0,0.25*safe_adb(CPOOLR(N,L,NZ,NY,NX) &
      ,RCO2N(N,L,NZ,NY,NX))))
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
  IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
    FZUP=AMIN1(safe_adb(CCPOLR(N,L,NZ,NY,NX),CCPOLR(N,L,NZ,NY,NX) &
      +CZPOLR(N,L,NZ,NY,NX)/ZCKI) &
      ,safe_adb(CPPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX) &
      +CZPOLR(N,L,NZ,NY,NX)/ZPKI))
    FPUP=AMIN1(safe_adb(CCPOLR(N,L,NZ,NY,NX),CCPOLR(N,L,NZ,NY,NX) &
      +CPPOLR(N,L,NZ,NY,NX)/PCKI) &
      ,safe_adb(CZPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX) &
      +CPPOLR(N,L,NZ,NY,NX)/PZKI))
  ELSE
    FZUP=0.0_r8
    FPUP=0.0_r8
  ENDIF
  !NN=0
  UPWTRP=AMAX1(0.0,-UPWTR(N,L,NZ,NY,NX)/PP(NZ,NY,NX))
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
  IF(ROXYY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FOXYX=AMAX1(FPP(N,L,NZ),ROXYP(N,L,NZ,NY,NX)/ROXYY(L,NY,NX))
  ELSE
    FOXYX=FPQ(N,L,NZ)
  ENDIF
  end subroutine GetUptakeCapcity
!------------------------------------------------------------------------

  subroutine RootSoilGasExchange(N,L,NZ,NY,NX,RRADL,FPQ,FRTDPX,RTARR,&
    UPWTRH,FOXYX,RUPOXT)

  implicit none
  integer , intent(in) :: N,L,NZ,NY,NX
  real(r8), intent(in) :: RRADL(2,JZ),FPQ(2,JZ,05),FRTDPX(JZ,05)
  real(r8), intent(in) :: RTARR(2,JZ),UPWTRH,FOXYX
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

  IF(RCO2M(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
    .AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
    .AND.FOXYX.GT.ZEROQ(NZ,NY,NX))THEN
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
    CO2A1=AMAX1(ZEROP(NZ,NY,NX),CO2A(N,L,NZ,NY,NX))
    CO2P1=AMAX1(ZEROP(NZ,NY,NX),CO2P(N,L,NZ,NY,NX))
    CO2G1=AMAX1(ZEROP(NZ,NY,NX),CO2G(L,NY,NX)*FPQ(N,L,NZ))
    CO2S1=AMAX1(ZEROP(NZ,NY,NX),CO2S(L,NY,NX)*FPQ(N,L,NZ))
    OXYA1=AMAX1(ZEROP(NZ,NY,NX),OXYA(N,L,NZ,NY,NX))
    OXYP1=AMAX1(ZEROP(NZ,NY,NX),OXYP(N,L,NZ,NY,NX))
    OXYG1=AMAX1(ZEROP(NZ,NY,NX),OXYG(L,NY,NX)*FOXYX)
    OXYS1=OXYS(L,NY,NX)*FOXYX
    CH4A1=CH4A(N,L,NZ,NY,NX)
    CH4P1=CH4P(N,L,NZ,NY,NX)
    CH4S1=CH4S(L,NY,NX)*FPQ(N,L,NZ)
    CCH4S1=CCH4S(L,NY,NX)
    CCH4P1=AMAX1(0.0,CH4P1/RTVLW(N,L,NZ,NY,NX))
    Z2OA1=Z2OA(N,L,NZ,NY,NX)
    Z2OP1=Z2OP(N,L,NZ,NY,NX)
    Z2OS1=Z2OS(L,NY,NX)*FPQ(N,L,NZ)
    CN2OS1=CZ2OS(L,NY,NX)
    CN2OP1=AMAX1(0.0,Z2OP1/RTVLW(N,L,NZ,NY,NX))
    ZH3A1=ZH3A(N,L,NZ,NY,NX)
    ZH3P1=ZH3P(N,L,NZ,NY,NX)
    ZH3S1=ZNH3S(L,NY,NX)*FPQ(N,L,NZ)
    ZH3B1=ZNH3B(L,NY,NX)*FPQ(N,L,NZ)
    CNH3S1=CNH3S(L,NY,NX)
    CNH3B1=CNH3B(L,NY,NX)
    CNH3P1=AMAX1(0.0,ZH3P1/RTVLW(N,L,NZ,NY,NX))
    H2GA1=H2GA(N,L,NZ,NY,NX)
    H2GP1=H2GP(N,L,NZ,NY,NX)
    H2GS1=H2GS(L,NY,NX)*FPQ(N,L,NZ)
    CH2GS1=CH2GS(L,NY,NX)
    CH2GP1=AMAX1(0.0,H2GP1/RTVLW(N,L,NZ,NY,NX))
    RTVLWA=RTVLW(N,L,NZ,NY,NX)*VLNH4(L,NY,NX)
    RTVLWB=RTVLW(N,L,NZ,NY,NX)*VLNHB(L,NY,NX)
    UPMXP=ROXYP(N,L,NZ,NY,NX)*XNPG/PP(NZ,NY,NX)
    ROXYFX=ROXYF(L,NY,NX)*FOXYX*XNPG
    RCO2FX=RCO2F(L,NY,NX)*FOXYX*XNPG
    ROXYLX=ROXYL(L,NY,NX)*FOXYX*XNPG
!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
!     *SGL1=diffusivity
!     PORTX=tortuosity effect of root porosity on diffusivity
!     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
!     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
!
    CGSGL1=CGSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
    OGSGL1=OGSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
    CHSGL1=CHSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
    Z2SGL1=Z2SGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
    ZHSGL1=ZHSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
    HGSGL1=HGSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
    CLSGL1=CLSGL(L,NY,NX)*XNPG*FOXYX
    OLSGL1=OLSGL(L,NY,NX)*XNPG*FOXYX
    CQSGL1=CQSGL(L,NY,NX)*XNPG*FOXYX
    ZVSGL1=ZVSGL(L,NY,NX)*XNPG*FOXYX
    ZNSGL1=ZNSGL(L,NY,NX)*XNPG*FOXYX
    HLSGL1=HLSGL(L,NY,NX)*XNPG*FOXYX
    OLSGLP=OLSGL(L,NY,NX)*XNPG
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
    IF(WTRTS(NZ,NY,NX).GT.ZEROP(NZ,NY,NX).AND.FRTDPX(L,NZ).GT.ZERO)THEN
      RTCR1=AMAX1(PP(NZ,NY,NX),RTN1(N,L,NZ,NY,NX)) &
        *PICON*RRAD1(N,L,NZ,NY,NX)**2/DPTHZ(L,NY,NX)
      RTCR2=(RTNL(N,L,NZ,NY,NX)*PICON*RRAD2(N,L,NZ,NY,NX)**2 &
        /RTLGA(N,L,NZ,NY,NX))/FRTDPX(L,NZ)
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
!     IDAY(1,=emergence date
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
    IF(N.EQ.1.AND.IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).GT.0 &
      .AND.RTLGP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RTARRX=RTARR(N,L)/RRADP(N,NZ,NY,NX)
      DIFOP=OLSGLP*RTARRX
      VOLWCA=RTVLW(N,L,NZ,NY,NX)*SCO2L(L,NY,NX)
      VOLWOA=RTVLW(N,L,NZ,NY,NX)*SOXYL(L,NY,NX)
      VOLWC4=RTVLW(N,L,NZ,NY,NX)*SCH4L(L,NY,NX)
      VOLWZA=RTVLW(N,L,NZ,NY,NX)*SN2OL(L,NY,NX)
      VOLWNA=RTVLW(N,L,NZ,NY,NX)*SNH3L(L,NY,NX)
      VOLWH2=RTVLW(N,L,NZ,NY,NX)*SH2GL(L,NY,NX)
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
    DFGP=AMIN1(1.0,XNPD*SQRT(PORT(N,NZ,NY,NX))*TFND(L,NY,NX))
    RCO2PX=-RCO2A(N,L,NZ,NY,NX)*XNPG
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
      VOLWMO=VOLWM(M,L,NY,NX)*FOXYX
      VOLWMM=VOLWM(M,L,NY,NX)*FPQ(N,L,NZ)
      VOLPMM=VOLPM(M,L,NY,NX)*FPQ(N,L,NZ)
      VOLWSP=RTVLW(N,L,NZ,NY,NX)+VOLWMM
      VOLWMA=VOLWMM*VLNH4(L,NY,NX)
      VOLWMB=VOLWMM*VLNHB(L,NY,NX)
      VOLWSA=RTVLWA+VOLWMA
      VOLWSB=RTVLWB+VOLWMB
      THETW1=AMAX1(0.0,VOLWM(M,L,NY,NX)/VOLY(L,NY,NX))
      IF(THETW1.GT.THETY(L,NY,NX).AND.FPQ(N,L,NZ).GT.ZEROQ(NZ,NY,NX))THEN
        THETM=TORT(M,L,NY,NX)*THETW1
        RRADS=LOG((FILM(M,L,NY,NX)+RRADL(N,L))/RRADL(N,L))
        RTARRX=RTARR(N,L)/RRADS
        DIFOL=THETM*OLSGL1*RTARRX
        DIFCL=THETM*CQSGL1*RTARRX
        DIFZL=THETM*ZVSGL1*RTARRX
        DIFNL=THETM*ZNSGL1*RTARRX*VLNH4(L,NY,NX)
        DIFNB=THETM*ZNSGL1*RTARRX*VLNHB(L,NY,NX)
        DIFHL=THETM*HLSGL1*RTARRX
        CH4G1=CCH4G(L,NY,NX)*VOLPMM
        Z2OG1=CZ2OG(L,NY,NX)*VOLPMM
        ZH3G1=CNH3G(L,NY,NX)*VOLPMM
        H2GG1=CH2GG(L,NY,NX)*VOLPMM
        VOLWCO=VOLWMM*SCO2L(L,NY,NX)
        VOLWOX=VOLWMM*SOXYL(L,NY,NX)
        VOLWCH=VOLWMM*SCH4L(L,NY,NX)
        VOLWN2=VOLWMM*SN2OL(L,NY,NX)
        VOLWNH=VOLWMM*SNH3L(L,NY,NX)*VLNH4(L,NY,NX)
        VOLWNB=VOLWMM*SNH3L(L,NY,NX)*VLNHB(L,NY,NX)
        VOLWHG=VOLWMM*SH2GL(L,NY,NX)
        VOLPNH=VOLPMM*VLNH4(L,NY,NX)
        VOLPNB=VOLPMM*VLNHB(L,NY,NX)
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
          COXYS1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX),AMAX1(0.0,OXYS1/VOLWMO))
          CCH4S1=AMAX1(0.0,CH4S1/VOLWMM)
          CN2OS1=AMAX1(0.0,Z2OS1/VOLWMM)
          CNH3S1=AMAX1(0.0,ZH3S1/VOLWMM)
          CNH3B1=AMAX1(0.0,ZH3B1/VOLWMM)
          CH2GS1=AMAX1(0.0,H2GS1/VOLWMM)
          IF(RTVLP(N,L,NZ,NY,NX).GT.ZERO)THEN
            CCO2A1=AMAX1(0.0,CO2A1/RTVLP(N,L,NZ,NY,NX))
            COXYA1=AMAX1(0.0,OXYA1/RTVLP(N,L,NZ,NY,NX))
            CCH4A1=AMAX1(0.0,CH4A1/RTVLP(N,L,NZ,NY,NX))
            CZ2OA1=AMAX1(0.0,Z2OA1/RTVLP(N,L,NZ,NY,NX))
            CNH3A1=AMAX1(0.0,ZH3A1/RTVLP(N,L,NZ,NY,NX))
            CH2GA1=AMAX1(0.0,H2GA1/RTVLP(N,L,NZ,NY,NX))
          ELSE
            CCO2A1=0.0_r8
            COXYA1=0.0_r8
            CCH4A1=0.0_r8
            CZ2OA1=0.0_r8
            CNH3A1=0.0_r8
            CH2GA1=0.0_r8
          ENDIF
          CCO2P1=AMAX1(0.0,CO2P1/RTVLW(N,L,NZ,NY,NX))
          COXYP1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX),AMAX1(0.0,OXYP1/RTVLW(N,L,NZ,NY,NX)))
          CCH4P1=AMAX1(0.0,CH4P1/RTVLW(N,L,NZ,NY,NX))
          CN2OP1=AMAX1(0.0,Z2OP1/RTVLW(N,L,NZ,NY,NX))
          CNH3P1=AMAX1(0.0,ZH3P1/RTVLW(N,L,NZ,NY,NX))
          CH2GP1=AMAX1(0.0,H2GP1/RTVLW(N,L,NZ,NY,NX))
          DIFOX=DIFOL+DIFOP
          RMFCOS=UPWTRH*CCO2S1
          RMFOXS=UPWTRH*COXYS1
          RMFCHS=UPWTRH*CCH4S1
          RMFN2S=UPWTRH*CN2OS1
          RMFN3S=UPWTRH*CNH3S1*VLNH4(L,NY,NX)
          RMFN3B=UPWTRH*CNH3B1*VLNHB(L,NY,NX)
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
          IF(X.GT.ZERO.AND.OXYS1.GT.ZEROP(NZ,NY,NX))THEN
            B=-UPMXP-DIFOX*OXKM-X
            C=X*UPMXP
            RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
            COXYR=(X-RUPOXR)/DIFOX
            RDFOXS=RMFOXS+DIFOL*(COXYS1-COXYR)
            RDFOXP=DIFOP*(COXYP1-COXYR)
          ELSE
            X=DIFOP*COXYP1
            IF(X.GT.ZERO.AND.OXYP1.GT.ZEROP(NZ,NY,NX))THEN
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

          RUPOSX=RDFOXS*PP(NZ,NY,NX)
          RUPOPX=RDFOXP*PP(NZ,NY,NX)
          RDFCOS=RMFCOS+DIFCL*(CCO2S1-CCO2P1)
          RDXCOS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),CO2S1) &
            -VOLWMM*AMAX1(ZEROP(NZ,NY,NX),CO2P1))/VOLWSP
          IF(RDFCOS.GT.0.0)THEN
            RCO2SX=AMIN1(AMAX1(0.0,RDXCOS),RDFCOS*PP(NZ,NY,NX))
          ELSE
            RCO2SX=AMAX1(AMIN1(0.0,RDXCOS),RDFCOS*PP(NZ,NY,NX))
          ENDIF
          IF(N.EQ.1)THEN
            RDFCHS=RMFCHS+DIFCL*(CCH4S1-CCH4P1)
            RDXCHS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),CH4S1) &
              -VOLWMM*AMAX1(ZEROP(NZ,NY,NX),CH4P1))/VOLWSP
            IF(RDFCHS.GT.0.0)THEN
              RUPCSX=AMIN1(AMAX1(0.0,RDXCHS),RDFCHS*PP(NZ,NY,NX))
            ELSE
              RUPCSX=AMAX1(AMIN1(0.0,RDXCHS),RDFCHS*PP(NZ,NY,NX))
            ENDIF
            RDFN2S=RMFN2S+DIFZL*(CN2OS1-CN2OP1)
            RDXN2S=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),Z2OS1) &
              -VOLWMM*AMAX1(ZEROP(NZ,NY,NX),Z2OP1))/VOLWSP
            IF(RDFN2S.GT.0.0)THEN
              RUPZSX=AMIN1(AMAX1(0.0,RDXN2S),RDFN2S*PP(NZ,NY,NX))
            ELSE
              RUPZSX=AMAX1(AMIN1(0.0,RDXN2S),RDFN2S*PP(NZ,NY,NX))
            ENDIF
            RDFN3S=RMFN3S+DIFNL*(CNH3S1-CNH3P1)
            IF(VOLWSA.GT.ZEROP(NZ,NY,NX))THEN
              ZH3PA=ZH3P1*VLNH4(L,NY,NX)
              RDXNHS=(RTVLWA*AMAX1(ZEROP(NZ,NY,NX),ZH3S1) &
                -VOLWMA*AMAX1(ZEROP(NZ,NY,NX),ZH3PA))/VOLWSA
            ELSE
              RDXNHS=0.0_r8
            ENDIF
            IF(RDFN3S.GT.0.0)THEN
              RUPNSX=AMIN1(AMAX1(0.0,RDXNHS),RDFN3S*PP(NZ,NY,NX))
            ELSE
              RUPNSX=AMAX1(AMIN1(0.0,RDXNHS),RDFN3S*PP(NZ,NY,NX))
            ENDIF
            RDFN3B=RMFN3B+DIFNB*(CNH3B1-CNH3P1)
            IF(VOLWSB.GT.ZEROP(NZ,NY,NX))THEN
              ZH3PB=ZH3P1*VLNHB(L,NY,NX)
              RDXNHB=(RTVLWB*AMAX1(ZEROP(NZ,NY,NX),ZH3B1) &
                -VOLWMB*AMAX1(ZEROP(NZ,NY,NX),ZH3PB))/VOLWSB
            ELSE
              RDXNHB=0.0_r8
            ENDIF
            IF(RDFN3B.GT.0.0)THEN
              RUPNBX=AMIN1(AMAX1(0.0,RDXNHB),RDFN3B*PP(NZ,NY,NX))
            ELSE
              RUPNBX=AMAX1(AMIN1(0.0,RDXNHB),RDFN3B*PP(NZ,NY,NX))
            ENDIF
            RDFHGS=RMFHGS+DIFHL*(CH2GS1-CH2GP1)
            RDXHGS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),H2GS1) &
              -VOLWMM*AMAX1(ZEROP(NZ,NY,NX),H2GP1))/VOLWSP
            IF(RDFHGS.GT.0.0)THEN
              RUPHGX=AMIN1(AMAX1(0.0,RDXHGS),RDFHGS*PP(NZ,NY,NX))
            ELSE
              RUPHGX=AMAX1(AMIN1(0.0,RDXHGS),RDFHGS*PP(NZ,NY,NX))
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
          IF(THETPM(M,L,NY,NX).GT.THETX)THEN
            DFGSP=FPQ(N,L,NZ)*DFGS(M,L,NY,NX)
            RCODFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),CO2G1)*VOLWCO &
              -(AMAX1(ZEROS(NY,NX),CO2S1)-RCO2SX)*VOLPMM)/(VOLWCO+VOLPMM)
            RUPOST=RUPOSX-ROXYLX
            ROXDFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),OXYG1)*VOLWOX &
              -(AMAX1(ZEROS(NY,NX),OXYS1)-RUPOST)*VOLPMM)/(VOLWOX+VOLPMM)
            IF(N.EQ.1)THEN
              RCHDFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),CH4G1)*VOLWCH &
                -(AMAX1(ZEROS(NY,NX),CH4S1)-RUPCSX)*VOLPMM)/(VOLWCH+VOLPMM)
              RN2DFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),Z2OG1)*VOLWN2 &
                -(AMAX1(ZEROS(NY,NX),Z2OS1)-RUPZSX)*VOLPMM)/(VOLWN2+VOLPMM)
              IF(VOLWNH+VOLPNH.GT.ZEROP(NZ,NY,NZ))THEN
                ZH3GA=ZH3G1*VLNH4(L,NY,NX)
                RNHDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROP(NZ,NY,NX),ZH3GA)*VOLWNH &
                  -(AMAX1(ZEROS(NY,NX),ZH3S1)-RUPNSX)*VOLPNH)/(VOLWNH+VOLPNH)))
              ELSE
                RNHDFQ=0.0_r8
              ENDIF
              IF(VOLWNB+VOLPNB.GT.ZEROP(NZ,NY,NZ))THEN
                ZH3GB=ZH3G1*VLNHB(L,NY,NX)
                RNBDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROP(NZ,NY,NX),ZH3GB)*VOLWNB &
                  -(AMAX1(ZEROS(NY,NX),ZH3B1)-RUPNBX)*VOLPNB)/(VOLWNB+VOLPNB)))
              ELSE
                RNBDFQ=0.0_r8
              ENDIF
              RHGDFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),H2GG1)*VOLWHG &
                -(AMAX1(ZEROS(NY,NX),H2GS1)-RUPHGX)*VOLPMM)/(VOLWHG+VOLPMM)
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
          IF(N.EQ.1.AND.RTVLP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
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
            RCODF1=AMAX1(-CO2PX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),CO2A1)*VOLWCA &
              -CO2PX*RTVLP(N,L,NZ,NY,NX))/(VOLWCA+RTVLP(N,L,NZ,NY,NX)))
            OXYPX=OXYP1-RUPOPX
            ROXDF1=AMAX1(-OXYPX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),OXYA1)*VOLWOA &
              -OXYPX*RTVLP(N,L,NZ,NY,NX))/(VOLWOA+RTVLP(N,L,NZ,NY,NX)))
            CH4PX=CH4P1+RUPCSX
            RCHDF1=AMAX1(-CH4PX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),CH4A1)*VOLWC4 &
              -CH4PX*RTVLP(N,L,NZ,NY,NX))/(VOLWC4+RTVLP(N,L,NZ,NY,NX)))
            Z2OPX=Z2OP1+RUPZSX
            RN2DF1=AMAX1(-Z2OPX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),Z2OA1)*VOLWZA &
              -Z2OPX*RTVLP(N,L,NZ,NY,NX))/(VOLWZA+RTVLP(N,L,NZ,NY,NX)))
            ZH3PX=ZH3P1+RUPNTX
            RNHDF1=AMAX1(-ZH3PX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),ZH3A1)*VOLWNA &
              -ZH3PX*RTVLP(N,L,NZ,NY,NX))/(VOLWNA+RTVLP(N,L,NZ,NY,NX)))
            H2GPX=H2GP1+RUPHGX
            RHGDF1=AMAX1(-H2GPX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),H2GA1)*VOLWH2 &
              -H2GPX*RTVLP(N,L,NZ,NY,NX))/(VOLWH2+RTVLP(N,L,NZ,NY,NX)))
            RCOFL1=AMIN1(DFCOA,RTVLP(N,L,NZ,NY,NX))*(CCO2E(NY,NX)-CCO2A1)
            ROXFL1=AMIN1(DFOXA,RTVLP(N,L,NZ,NY,NX))*(COXYE(NY,NX)-COXYA1)
            RCHFL1=AMIN1(DFCHA,RTVLP(N,L,NZ,NY,NX))*(CCH4E(NY,NX)-CCH4A1)
            RN2FL1=AMIN1(DFN2A,RTVLP(N,L,NZ,NY,NX))*(CZ2OE(NY,NX)-CZ2OA1)
            RNHFL1=AMIN1(DFNHA,RTVLP(N,L,NZ,NY,NX))*(CNH3E(NY,NX)-CNH3A1)
            RHGFL1=AMIN1(DFHGA,RTVLP(N,L,NZ,NY,NX))*(CH2GE(NY,NX)-CH2GA1)
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
          RCO2S(N,L,NZ,NY,NX)=RCO2S(N,L,NZ,NY,NX)+RCO2SX
          RUPOXS(N,L,NZ,NY,NX)=RUPOXS(N,L,NZ,NY,NX)+RUPOSX
          RUPCHS(N,L,NZ,NY,NX)=RUPCHS(N,L,NZ,NY,NX)+RUPCSX
          RUPN2S(N,L,NZ,NY,NX)=RUPN2S(N,L,NZ,NY,NX)+RUPZSX
          RUPN3S(N,L,NZ,NY,NX)=RUPN3S(N,L,NZ,NY,NX)+RUPNSX
          RUPN3B(N,L,NZ,NY,NX)=RUPN3B(N,L,NZ,NY,NX)+RUPNBX
          RUPHGS(N,L,NZ,NY,NX)=RUPHGS(N,L,NZ,NY,NX)+RUPHGX
!
!     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE
!
!     R*DFA=root aqueous-gaseous CO2 exchange
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!
          RCODFA(N,L,NZ,NY,NX)=RCODFA(N,L,NZ,NY,NX)+RCODF1
          ROXDFA(N,L,NZ,NY,NX)=ROXDFA(N,L,NZ,NY,NX)+ROXDF1
          RCHDFA(N,L,NZ,NY,NX)=RCHDFA(N,L,NZ,NY,NX)+RCHDF1
          RN2DFA(N,L,NZ,NY,NX)=RN2DFA(N,L,NZ,NY,NX)+RN2DF1
          RNHDFA(N,L,NZ,NY,NX)=RNHDFA(N,L,NZ,NY,NX)+RNHDF1
          RHGDFA(N,L,NZ,NY,NX)=RHGDFA(N,L,NZ,NY,NX)+RHGDF1
          RCOFLA(N,L,NZ,NY,NX)=RCOFLA(N,L,NZ,NY,NX)+RCOFL1
          ROXFLA(N,L,NZ,NY,NX)=ROXFLA(N,L,NZ,NY,NX)+ROXFL1
          RCHFLA(N,L,NZ,NY,NX)=RCHFLA(N,L,NZ,NY,NX)+RCHFL1
          RN2FLA(N,L,NZ,NY,NX)=RN2FLA(N,L,NZ,NY,NX)+RN2FL1
          RNHFLA(N,L,NZ,NY,NX)=RNHFLA(N,L,NZ,NY,NX)+RNHFL1
          RHGFLA(N,L,NZ,NY,NX)=RHGFLA(N,L,NZ,NY,NX)+RHGFL1
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2P=root CO2 emission into root
!     RUPOXP=root O2 uptake from root
!     ROXSK=total O2 uptake from soil by all microbial,root popns
!
          RCO2P(N,L,NZ,NY,NX)=RCO2P(N,L,NZ,NY,NX)+RCO2PX+RCO2SX
          RUPOXP(N,L,NZ,NY,NX)=RUPOXP(N,L,NZ,NY,NX)+RUPOPX
          ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RUPOSX

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
    RUPOXT=RUPOXP(N,L,NZ,NY,NX)+RUPOXS(N,L,NZ,NY,NX)
    WFR(N,L,NZ,NY,NX)=AMIN1(1.0,AMAX1(0.0 &
      ,RUPOXT/ROXYP(N,L,NZ,NY,NX)))
  ELSE
    RUPOXT=0.0_r8
    IF(L.GT.NG(NZ,NY,NX))THEN
      WFR(N,L,NZ,NY,NX)=WFR(N,L-1,NZ,NY,NX)
    ELSE
      WFR(N,L,NZ,NY,NX)=1.0
    ENDIF
  ENDIF
  end subroutine RootSoilGasExchange
!------------------------------------------------------------------------------------------

  subroutine RootMycoO2NutrientUptake(NZ,NY,NX,OSTRN,OSTRD,PATH,RRADL,&
    FPQ,FPP,FRTDPX,RTARR)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: PATH(2,JZ),RRADL(2,JZ),FPQ(2,JZ,05),FPP(2,JZ,05)
  real(r8), intent(in) :: FRTDPX(JZ,05),RTARR(2,JZ)
  real(r8), intent(inout) :: OSTRN,OSTRD
  real(r8) :: TFOXYX
  real(r8) :: FCUP,FZUP,FPUP,FWSRT,UPWTRP,UPWTRH,FOXYX,RUPOXT
  integer :: N,L
!     begin_execution

  DO 955 N=1,MY(NZ,NY,NX)
    DO 950 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX) &
        .AND.RTDNP(N,L,NZ,NY,NX).GT.ZERO &
        .AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.THETW(L,NY,NX).GT.ZERO)THEN
        TFOXYX=0.0_r8
        call GetUptakeCapcity(N,L,NZ,NY,NX,FPQ,FPP,FCUP,FZUP,FPUP,&
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
        ROXYP(N,L,NZ,NY,NX)=2.667*RCO2M(N,L,NZ,NY,NX)

        call RootSoilGasExchange(N,L,NZ,NY,NX,RRADL,FPQ,FRTDPX,RTARR,UPWTRH,&
          FOXYX,RUPOXT)

        OSTRD=OSTRD+ROXYP(N,L,NZ,NY,NX)
        OSTRN=OSTRN+RUPOXT

        call RootExudates(N,L,NZ,NY,NX)
!
!     NUTRIENT UPTAKE
!
!     WFR=constraint by O2 consumption on all biological processes
!     FCUP=limitation to active uptake respiration from CPOOLR
!     FWSRT=protein concentration relative to 5%
!     RTLGP=root,myco length per plant
!
        IF(WFR(N,L,NZ,NY,NX).GT.ZERO.AND.FCUP.GT.ZERO.AND.FWSRT.GT.ZERO &
          .AND.RTLGP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
!
!     FZUP=limitn to active uptake respiration from CZPOLR
!
          call UptakeMineralNitrogen(N,L,NZ,NY,NX,PATH,RRADL,FPQ,FPP,&
            RTARR,FCUP,FZUP,FWSRT,UPWTRP)
!
!     FPUP=limitn to active uptake respiration from CPPOLR
!
          call UptakeMineralPhosporhus(N,L,NZ,NY,NX,PATH,RRADL,FPQ,FPP,&
            RTARR,FCUP,FPUP,FWSRT,UPWTRP)

        ELSE
          call NoActiveNutrientUptake(N,L,NZ,NY,NX)

        ENDIF
      ELSE
        call ZeroUptake(N,L,NZ,NY,NX)
      ENDIF

      call SumupNutrientUptake(N,L,NZ,NY,NX)

950 CONTINUE
955 CONTINUE
  end subroutine RootMycoO2NutrientUptake

end module UptakeMod
