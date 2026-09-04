module MicrobMathFuncMod
  use EcosimConst
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use EcoSIMSolverPar,      only: dts_gas, NPH, NPT
  use MicFLuxTypeMod,       only: micfluxtype
  use MicForcTypeMod,       only: micforctype
  use MicStateTraitTypeMod, only: micsttype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use ElmIDMod,             only: ielmc, ielmn, ibiom_kinetic, ibiom_struct
  use MicrobeDiagTypes,     only: Cumlate_Flux_Diag_type, Microbe_Diag_type, &
                                  Microbe_Flux_type, Microbe_State_type
  use minimathmod,          only: AZMAX1, safe_adb, fixEXConsumpFlux, &
                                  real_truncate, SubstrateDribbling
  use NitroPars,            only: BIOS, FMN, ORAD, RMOM
  implicit none

  character(len=*), parameter, private :: mod_filename=&
  __FILE__

  real(r8), parameter, private :: ZERO=1.0E-15_r8

  contains
!------------------------------------------------------------------------

  subroutine MicrobPhysTempFun(TKSO, TSensGrowth, TSensMaintR)
  !
  !the physiological temperature dependence of microbes
  implicit none
  real(r8), intent(in) :: TKSO
  real(r8), intent(out):: TSensGrowth !temperature sensitivity for growth respiration
  real(r8), intent(out):: TSensMaintR !temperature sensitivity for maintenance respiration
  real(r8) :: RTK,STK
  real(r8) :: ACTV,ACTVM

  RTK=RGASC*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TSensGrowth=EXP(25.229_r8-62500._r8/RTK)/ACTV
  ACTVM=1+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TSensMaintR=EXP(25.214_r8-62500._r8/RTK)/ACTVM

  end subroutine MicrobPhysTempFun
!------------------------------------------------------------------------


  pure function TranspBasedsubstrateUptake(S_conc,diffusc, KM, V_max, zeros)result(uptake)
  !
  !transport based substrate uptake
  !assuming balance between MM-based uptake and supply by diffusion
  !q=V_max*s/(KM+s)=D(S_ext-S)
  !solve for q.
  implicit none
  real(r8), intent(in) :: S_conc   !substrate concentration
  real(r8), intent(in) :: diffusc !diffusion coefficient
  real(r8), intent(in) :: KM      !half saturation parameter
  real(r8), intent(in) :: V_max   !maximum uptake rate
  real(r8), optional, intent(in) :: zeros !threshold for active uptake
  real(r8) :: uptake
  real(r8) :: X,B,C
  real(r8) :: zero1

  if(present(zeros))then
    zero1=zeros
  else
    zero1=zero
  endif

  !obtain uptake flux
  X=diffusc*S_conc
  IF(X.GT.ZERO1)THEN
    B=-V_max-diffusc*KM-X
    C=X*V_max
    uptake=(-B-SQRT(B*B-4.0_r8*C))/2.0_r8
  ELSE
    uptake=0.0_r8
  ENDIF
  end function TranspBasedsubstrateUptake

!------------------------------------------------------------------------------------------

  subroutine AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: OXKX
  real(r8), intent(in) :: FOXYX  !preallocated O2 for uptake
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx
  integer  :: M,MX
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,O2AquaDiffusvity1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMAX
  real(r8) :: ROXYFX       !gas dissolution flux
  real(r8) :: ROXYLX       !aqueous flux transport flux
  real(r8) :: RRADO,RMPOX
  real(r8) :: ROXDFQ  !gas dissolution
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X
  real(r8) :: VOLWPM,VOLOXM
  real(r8) :: dsignO2,dribbling_flx
  ! begin_execution
  associate(                                                 &
    OxyLimterHeter         => nmics%OxyLimterHeter,          &
    OMActHeter             => nmics%OMActHeter,              &
    RO2UptkHeter           => nmicf%RO2UptkHeter,            &
    RespGrossHeter         => nmicf%RespGrossHeter,          &
    RO2DmndHeter           => nmicf%RO2DmndHeter,            &
    O2_irrig_conc          => micfor%O2_irrig_conc,          &
    O2_rain_conc           => micfor%O2_rain_conc,           &
    COXYE                  => micfor%COXYE,                  &
    RO2GasXchangePrev      => micfor%RO2GasXchangePrev,      &
    RO2AquaXchangePrev     => micfor%RO2AquaXchangePrev,     &
    Irrig2LitRSurf_col     => micfor%Irrig2LitRSurf_col,     &
    Rain2LitRSurf          => micfor%Rain2LitRSurf,          &
    litrm                  => micfor%litrm,                  &
    O2AquaDiffusvity       => micfor%O2AquaDiffusvity,       &
    VLSoilPoreMicP         => micfor%VLSoilPoreMicP,         &
    VLSoilMicP             => micfor%VLSoilMicP,             &
    ZERO                   => micfor%ZERO,                   &
    ZEROS                  => micfor%ZEROS,                  &
    VLsoiAirPM             => micfor%VLsoiAirPM,             &
    VLWatMicP              => micfor%VLWatMicP,              &
    VLWatMicPM             => micfor%VLWatMicPM,             &
    THETPM                 => micfor%THETPM,                 &
    DiffusivitySolutEff    => micfor%DiffusivitySolutEff,    &
    FILM                   => micfor%FILM,                   &
    TortMicPM              => micfor%TortMicPM,              &
    OXYG                   => micstt%OXYG,                   &
    OXYS                   => micstt%OXYS,                   &
    COXYS                  => micstt%COXYS,                  &
    O2GSolubility          => micstt%O2GSolubility,          &
    COXYG                  => micstt%COXYG,                  &
    RNO2DmndReduxSoilHeter => micflx%RNO2DmndReduxSoilHeter, &
    RNO2DmndReduxBandHeter => micflx%RNO2DmndReduxBandHeter, &
    REcoUptkSoilO2M        => micflx%REcoUptkSoilO2M         &
  )

  IF(RO2DmndHeter(NGL,K).GT.ZEROS .AND. FOXYX.GT.ZERO)THEN
    IF(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMAX            = RO2DmndHeter(NGL,K)*dts_gas
      ROXYFX            = -RO2GasXchangePrev*dts_gas*FOXYX
      O2AquaDiffusvity1 = O2AquaDiffusvity*dts_gas

      IF(.not.litrm)THEN
        OXYG1  = OXYG*FOXYX
        ROXYLX = -RO2AquaXchangePrev*dts_gas*FOXYX
      ELSE
        !litter layer
        OXYG1  = COXYG*VLsoiAirPM(1)*FOXYX
        ROXYLX = -(RO2AquaXchangePrev+Rain2LitRSurf*O2_rain_conc+Irrig2LitRSurf_col*O2_irrig_conc)*dts_gas*FOXYX
      ENDIF
      OXYS1=OXYS*FOXYX

      if(OXYS1 <= 0._r8 .and. ROXYLX > 0._r8)ROXYLX=0._r8
      IF(OXYG1 <=0._r8 .and. ROXYFX>0._R8)ROXYFX=0._r8
!
      !write(*,*)'O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP'
!     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
!
      dribbling_flx=0._r8
      D420: DO M=1,NPH
        !
        !     ACTUAL REDUCTION OF AQUEOUS BY AEROBES CALCULATED
        !     FROM MASS FLOW PLUS DIFFUSION = ACTIVE UPTAKE
        !     COUPLED WITH DISSOLUTION OF GASEOUS O2 DURING REDUCTION
        !     OF AQUEOUS O2 FROM DISSOLUTION RATE CONSTANT 'DiffusivitySolutEff'
        !     CALCULATED IN 'WATSUB'
        !
        !     VLWatMicPM,VLsoiAirPM,VLSoilPoreMicP=water, air and total volumes
        !     ORAD=microbial radius,FILM=water film thickness
        !     DIFOX=aqueous O2 diffusion, TortMicPM=tortuosity
        !     BIOS=microbial number, OMA=active biomass
        !     O2GSolubility=O2 solubility, OXKX=Km for O2 uptake
        !     OXYS,COXYS=aqueous O2 amount, concentration
        !     OXYG,COXYG=gaseous O2 amount, concentration
        !     RMPOX,REcoUptkSoilO2M=O2 uptake
        !
        THETW1 = AZMAX1(safe_adb(VLWatMicPM(M),VLSoilMicP))
        RRADO  = ORAD*(FILM(M)+ORAD)/FILM(M)
        DIFOX  = TortMicPM(M)*O2AquaDiffusvity1*12.57_r8*BIOS*OMActHeter(NGL,K)*RRADO
        VOLWOX = VLWatMicPM(M)*O2GSolubility
        VOLPOX = VLsoiAirPM(M)
        VOLWPM = VOLWOX+VOLPOX
        VOLOXM = VLWatMicPM(M)*FOXYX
        D425: DO MX=1,NPT
          call fixEXConsumpFlux(OXYG1,ROXYFX)
          call fixEXConsumpFlux(OXYS1,ROXYLX)

          COXYS1 = AMIN1(COXYE*O2GSolubility,safe_adb(OXYS1,VOLOXM))

          !obtain O2 uptake flux
          if(OXYS1<=ZEROS)then
            RMPOX=0.0_r8
          else
            RMPOX=TranspBasedsubstrateUptake(COXYS1,DIFOX, OXKX, RUPMAX, ZEROS)
          endif

          !apply the uptake
          call SubstrateDribbling(RMPOX,dribbling_flx,OXYS1)

          !apply dissolution-volatilization
          IF(THETPM(M).GT.AirFillPore_Min.AND.VOLPOX.GT.ZEROS)THEN
            ROXDFQ=DiffusivitySolutEff(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
            ROXDFQ=AMAX1(AMIN1(ROXDFQ,OXYG1),-OXYS1)
          ELSE
            ROXDFQ=0.0_r8
          ENDIF
          OXYG1 = OXYG1-ROXDFQ
          OXYS1 = OXYS1+ROXDFQ
          !accumulate upatke
          RO2UptkHeter(NGL,K) = RO2UptkHeter(NGL,K)+RMPOX
          REcoUptkSoilO2M(M)     = REcoUptkSoilO2M(M)+RMPOX
        ENDDO D425

      ENDDO D420
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (OxyLimterHeter)
      !
      !     OxyLimterHeter=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter,RNO2DmndReduxBandHeter=NH3,NO2 oxidation in non-band, band
      !
      OxyLimterHeter(NGL,K)=AMIN1(1.0_r8,AZMAX1(RO2UptkHeter(NGL,K)/RO2DmndHeter(NGL,K)))

    ELSE
      RO2UptkHeter(NGL,K)   = RO2DmndHeter(NGL,K)
      OxyLimterHeter(NGL,K) = 1.0_r8
    ENDIF
  ELSE
    RO2UptkHeter(NGL,K)   = 0.0_r8
    OxyLimterHeter(NGL,K) = 1.0_r8
  ENDIF
  !
  !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2ProdHeter,RAcetateProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
  !     RO2Uptk4RespHeter=O2-limited O2 uptake
  !     RSMetaOxidSoilAutor,RSMetaOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !

  end associate
  end subroutine AerobicHeterO2Uptake

!------------------------------------------------------------------------------------------
  subroutine StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC,FOQA,micfor,naqfdiag,nmicdiag,nmics)
  implicit none
  integer, intent(in) :: N   !functional group id
  integer, intent(in) :: NGL !functional guild id
  integer, intent(in) :: K   !complex id
  real(r8),intent(in) :: TotActMicrobiom         !total active microbial biomass
  type(micforctype), intent(in) :: micfor
  real(r8),intent(out) :: FOQC                   !fraction of DOC acetate demand over all microbial demand, soil/band
  real(r8),intent(out) :: FOQA                   !fraction of acetate demand over all microbial demand, soil/band
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type),intent(inout) :: nmicdiag
  type(Microbe_State_type), intent(inout):: nmics
  real(r8) :: WatStressMicb
  associate(                                                 &
    PSISoilMatricP         => micfor%PSISoilMatricP,         &
    ZEROS                  => micfor%ZEROS,                  &
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,        &
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,    &
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,  &
    mid_Aerob_Fungi        => micpar%mid_Aerob_Fungi,        &
    mid_Facult_DenitBacter => micpar%mid_Facult_DenitBacter, &
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,      &
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,     &
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,      &
    FracOMActHeter         => nmics%FracOMActHeter,          &
    OMActHeter             => nmics%OMActHeter,              &
    TempMaintRHeter        => nmics%TempMaintRHeter,         &
    TOMEK                  => nmicdiag%TOMEK          ,      &
    FracNO2ReduxHeter      => nmics%FracNO2ReduxHeter,       &
    WSensGroHeter          => nmics%WSensGroHeter         ,  &
    TSensGroHeter          => nmics%TSensGroHeter         ,  &
    TSensMaintR            => nmicdiag%TSensMaintR,          &
    TotBiomNO2Consumers    => nmicdiag%TotBiomNO2Consumers,  &
    TSensGrowth            => nmicdiag%TSensGrowth           &
  )

  ! WatStressMicb=water potential (PSISoilMatricP_vr) effect on microbial respiration
  ! OXKX=Km for O2 uptake
  ! OXKM=Km for heterotrophic O2 uptake set in starts.f
  !  GrowthEnvScalHeter=combined temp and water stress effect on growth respiration
  !  TempMaintRHeter=temperature effect on maintenance respiration
  !
  !different guilds can have different temperature and moisture sensitivity
  IF(N.EQ.mid_Aerob_Fungi)THEN
    WatStressMicb=EXP(0.1_r8*AMAX1(PSISoilMatricP,-500._r8))
  ELSE
    WatStressMicb=EXP(0.2_r8*AMAX1(PSISoilMatricP,-500._r8))
  ENDIF

  WSensGroHeter(NGL,K)=real_truncate(WatStressMicb,1.e-3_r8)
  TSensGroHeter(NGL,K)=TSensGrowth

  GrowthEnvScalHeter(NGL,K) = WSensGroHeter(NGL,K)*TSensGroHeter(NGL,K)
  TempMaintRHeter(NGL,K)    = TSensMaintR

  ! FracOMActHeter,FracNO2ReduxHeter=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActHeter(NGL,K)=OMActHeter(NGL,K)/TotActMicrobiom
  ELSE
    FracOMActHeter(NGL,K)=1.0_r8
  ENDIF

  IF(TotBiomNO2Consumers.GT.ZEROS .and. N.EQ.mid_Facult_DenitBacter)THEN
    FracNO2ReduxHeter(NGL,K)=OMActHeter(NGL,K)/TotBiomNO2Consumers
  ELSE
    FracNO2ReduxHeter(NGL,K)=1.0_r8
  ENDIF

  IF(TOMEK(ielmc,K).GT.ZEROS)THEN
    FracHeterBiomOfActK(NGL,K)=OMActHeter(NGL,K)/TOMEK(ielmc,K)
  ELSE
    FracHeterBiomOfActK(NGL,K)=1.0_r8
  ENDIF
  !
  IF(RDOMEcoDmndPrev(K).GT.ZEROS)THEN
    FOQC=AMAX1(FMN,RDOCUptkHeterPrev(NGL,K)/RDOMEcoDmndPrev(K))
  ELSE
    FOQC=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF

  naqfdiag%TFOQC=naqfdiag%TFOQC+FOQC

  IF(RAcetateEcoDmndPrev(K).GT.ZEROS)THEN
    FOQA=AMAX1(FMN,RAcetateUptkHeterPrev(NGL,K)/RAcetateEcoDmndPrev(K))
  ELSE
    FOQA=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF
  naqfdiag%TFOQA  = naqfdiag%TFOQA+FOQA

  end associate
  end subroutine StageFuncGuild

!------------------------------------------------------------------------------------------

  subroutine StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  real(r8) :: WatStressMicb

  associate(                                             &
    PSISoilMatricP      => micfor%PSISoilMatricP,        &
    OMActAutor          => nmics%OMActAutor,             &
    FracOMActAutor      => nmics%FracOMActAutor,         &
    FracNO2XupAutor     => nmics%FracNO2XupAutor,        &
    FracAutorBiomOfActK => nmics%FracAutorBiomOfActK,    &
    ZEROS               => micfor%ZEROS,                 &
    TSensMaintRAutor    => nmics%TSensMaintRAutor,       &
    GrowthEnvScalAutor  => nmics%GrowthEnvScalAutor,     &
    TotBiomNO2Consumers => nmicdiag%TotBiomNO2Consumers, &
    TotActMicrobiom     => nmicdiag%TotActMicrobiom,     &
    WSensGroAutor       => nmics%WSensGroAutor,          &
    TSensGroAutor       => nmics%TSensGroAutor,          &
    TSensGrowth         => nmicdiag%TSensGrowth,         &
    TSensMaintR         => nmicdiag%TSensMaintR          &
  )
  !replace with trait specific parameterization
  WatStressMicb           = EXP(0.2_r8*PSISoilMatricP)
  WSensGroAutor(NGL)      = WatStressMicb
  TSensGroAutor(NGL)      = TSensGrowth

  GrowthEnvScalAutor(NGL) = AZMAX1(TSensGroAutor(NGL)*WSensGroAutor(NGL))
  TSensMaintRAutor(NGL)   = AZMAX1(TSensMaintR)

! FracOMActHeter,FOMN=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActAutor(NGL)=OMActAutor(NGL)/TotActMicrobiom
  ELSE
    FracOMActAutor(NGL)=1.0_r8
  ENDIF

  IF(TotBiomNO2Consumers.GT.ZEROS.and.N.eq.micpar%mid_AutoAmmoniaOxidBacter)THEN
    FracNO2XupAutor(NGL)=OMActAutor(NGL)/TotBiomNO2Consumers
  ELSE
    FracNO2XupAutor(NGL)=1.0_r8
  ENDIF

  IF(TOMEAutoKC.GT.ZEROS)THEN
    FracAutorBiomOfActK(NGL)=OMActAutor(NGL)/TOMEAutoKC
  ELSE
    FracAutorBiomOfActK(NGL)=1.0_r8
  ENDIF
  end associate
  end subroutine StageAutotroph
!------------------------------------------------------------------------------------------
  subroutine CalcRespMaintAutor(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: I,J,NGL
  real(r8), intent(in) :: RMOMK(2)  !effect of low microbial C concentration on maintenance respiration
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  character(len=*),parameter :: subname='CalcRespMaintAutor'
  REAL(R8) :: FPH,RMOMX
  integer :: MID1

  associate(                                             &
    OMActAutor           => nmics%OMActAutor,            &
    OMN2Autor            => nmics%OMN2Autor,             &
    TSensMaintRAutor     => nmics%TSensMaintRAutor,      &
    RMaintDmndAutor      => nmicf%RMaintDmndAutor,       &
    RMaintRespAutor      => micflx%RMaintRespAutor,      &
    pH                   => micfor%pH,                   &
    mBiomeAutor          => micstt%mBiomeAutor,          &
    JGniA                => micpar%JGniA,                &
    JGnfA                => micpar%JGnfA                 &
  )
  call PrintInfo('beg '//subname)

  MID1                   = micpar%get_micb_id(ibiom_kinetic,NGL)
  FPH                    = 1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX                  = RMOM*TSensMaintRAutor(NGL)*FPH
  RMaintDmndAutor(ibiom_kinetic,NGL) = mBiomeAutor(ielmn,MID1)*RMOMX*RMOMK(ibiom_kinetic)
  RMaintDmndAutor(ibiom_struct,NGL)  = OMN2Autor(NGL)*RMOMX*RMOMK(ibiom_struct)
  !
  !     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
  !
  !     RMaintRespAutor=total maintenance respiration, as measured by C needs to be respired
  !     RGrowthRespAutor=growth respiration
  !     RMaintDefcitcitAutor=senescence respiration
  !
  RMaintRespAutor(NGL)      = RMaintDmndAutor(ibiom_kinetic,NGL)+RMaintDmndAutor(ibiom_struct,NGL)
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcRespMaintAutor
!------------------------------------------------------------------------------------------

  subroutine CalcRespMaintHeter(NGL,K,RMOMK,micfor,micstt,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,K
  real(r8), intent(in) :: RMOMK(2)  !effect of low microbial C concentration on maintenance respiration
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  character(len=*),parameter :: subname='CalcRespMaintHeter'
  REAL(R8) :: FPH,RMOMX
  integer :: MID1

  associate(                                             &
    OMN2                 => nmics%OMN2,                  &
    TempMaintRHeter      => nmics%TempMaintRHeter,       &
    RMaintDmndHeter      => nmicf%RMaintDmndHeter,       &
    RMaintRespHeter      => nmicf%RMaintRespHeter,       &
    pH                   => micfor%pH,                   &
    mBiomeHeter          => micstt%mBiomeHeter           &
  )
  call PrintInfo('beg '//subname)

  !     RMOMK=effect of low microbial C concentration on mntc respn

  FPH                                  = 1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX                                = RMOM*TempMaintRHeter(NGL,K)*FPH
  MID1                                 = micpar%get_micb_id(ibiom_kinetic,NGL)
  RMaintDmndHeter(ibiom_kinetic,NGL,K) = mBiomeHeter(ielmn,MID1,K)*RMOMX*RMOMK(ibiom_kinetic)
  RMaintDmndHeter(ibiom_struct,NGL,K)  = OMN2(NGL,K)*RMOMX*RMOMK(ibiom_struct)
  !
  !     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
  !
  !     RMaintRespHeter=total maintenance respiration
  !     RGrowthRespHeter=growth respiration
  !     RMaintDefcitcitHeter=senescence respiration
  !
  RMaintRespHeter(NGL,K)      = RMaintDmndHeter(ibiom_kinetic,NGL,K)+RMaintDmndHeter(ibiom_struct,NGL,K)
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcRespMaintHeter

end module MicrobMathFuncMod
