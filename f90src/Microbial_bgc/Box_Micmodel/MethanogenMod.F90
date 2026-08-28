module MethanogenMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use minimathmod,          only: AZMAX1, real_truncate
  use TracerIDMod
  use EcosimConst,          only: RGASC
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: CalcRespMaint, StageAutotroph, StageFuncGuild

  implicit none

  private
  public :: AcetoMethanogenCatabolism, H2MethanogensCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AcetoMethanogenCatabolism(N,K,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  implicit none
  integer, intent(in) :: N,K

  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(OMCplx_State_type),intent(in):: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  integer :: NGL
  reaL(r8) :: RGOMP         !substrate-limited potential respiration
  real(r8) :: GOMX,GOMM
  real(r8) :: RGOGY,RGOGZ
  real(r8) :: RGroMax   !kinetically unlimited acetate uptake
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQA
  real(r8)  :: WatStressMicb

! begin_execution
  associate(                                            &
    FBiomStoiScalarHeter => nmics%FBiomStoiScalarHeter, &
    OMActHeter           => nmics%OMActHeter,           &
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,   &
    FSBSTHeter           => nmicdiag%FSBSTHeter,        &
    RO2DmndHeter         => nmicf%RO2DmndHeter,         &
    RO2Dmnd4RespHeter    => nmicf%RO2Dmnd4RespHeter,    &
    ROQC4HeterMicrobAct  => nmicf%ROQC4HeterMicrobAct,  &
    RAcetateProdHeter    => nmicf%RAcetateProdHeter,    &
    RO2Uptk4RespHeter    => nmicf%RO2Uptk4RespHeter,    &
    RCO2ProdHeter        => nmicf%RCO2ProdHeter,        &
    RH2ProdHeter         => nmicf%RH2ProdHeter,         &
    RCH4ProdHeter        => nmicf%RCH4ProdHeter,        &
    FOQC                 => nmicf%FOQC,                 &
    ECHZHeter            => nmicf%ECHZHeter,            &
    FGOCP                => nmicf%FGOCP,                &
    RespGrossHeter       => nmicf%RespGrossHeter,       &
    FGOAP                => nmicf%FGOAP,                &
    TotActMicrobiom      => nmicdiag%TotActMicrobiom,   &
    CDOM                 => ncplxs%CDOM,                &
    DOM                  => micstt%DOM,                 &
    RO2DmndHetert        => micflx%RO2DmndHetert,       &
    RDOCUptkHeter        => micflx%RDOCUptkHeter,       &
    RAcetateUptkHeter    => micflx%RAcetateUptkHeter,   &
    PSISoilMatricP       => micfor%PSISoilMatricP,      &
    ZERO                 => micfor%ZERO,                &
    TKS                  => micfor%TKS                  &
  )
  !
  !     GOMX=acetate effect on energy yield
  !     ECHZHeter=growth respiration efficiency of aceto. methanogenesis
  !
  !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle
    WatStressMicb=real_truncate(EXP(0.2_r8*AMAX1(PSISoilMatricP,-500._r8)),1.e-3_r8)
    !prepare parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA,micfor,naqfdiag,nmicdiag,nmics)

    GOMX = RGASC*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI))
    GOMM = GOMX/24.0_r8
    ECHZHeter(NGL,K) = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GC4X+GOMM))/EOMH)))
    !
    !     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
    !     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
    !     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
    !     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
    !     MICROBIAL COMPETITION FACTOR
    !
    !     COQA=DOA concentration
    !     OQKAM=Km for acetate uptake,FBiomStoiScalarHeter=N,P limitation
    !     VMXCH4gAcet=specific respiration rate
    !     WatStressMicb=water stress effect, OMA=active biomass
    !     TSensGrowth=temp stress effect, FOQA= acetate limitation
    !     RGroMax=substrate-limited respiration of acetate
    !     RGroMax=competition-limited respiration of acetate
    !     OQA=acetate, FOQA=fraction of biological demand for acetate
    !     RGOMP=O2-unlimited respiration of acetate
    !     ROXY*=O2 demand, RDOCUptkHeter,ROQCA=DOC, acetate demand
    !     ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
    !
    FSBSTHeter(NGL,K)          = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKAM)
    RGOGY                      = FBiomStoiScalarHeter(NGL,K)*VMXCH4gAcet*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)
    RGOGZ                      = RGOGY*FSBSTHeter(NGL,K)
    RGroMax                    = AZMAX1(DOM(idom_acetate,K)*FOQA*ECHZHeter(NGL,K))
    RGOMP                      = AMIN1(RGroMax,RGOGZ)
    FGOCP(NGL,K)               = 0.0_r8
    FGOAP(NGL,K)               = 1.0_r8
    RO2Dmnd4RespHeter(NGL,K)   = 0.0_r8
    RO2DmndHeter(NGL,K)        = 0.0_r8
    RO2DmndHetert(NGL,K)       = 0.0_r8
    RDOCUptkHeter(NGL,K)       = 0.0_r8
    RAcetateUptkHeter(NGL,K)   = RGOGZ
    ROQC4HeterMicrobAct(NGL,K) = 0.0_r8

    !given CH3COOH -> CH4+CO2, 0.5 is into CH4.
    naqfdiag%tCH4ProdAceto=naqfdiag%tCH4ProdAceto+0.5_r8*RGOMP

    ! CH3COOH -> CO2 + CH4
    RespGrossHeter(NGL,K)   = RGOMP
    RCO2ProdHeter(NGL,K)    = 0.50_r8*RespGrossHeter(NGL,K)
    RAcetateProdHeter(NGL,K)  = 0.0_r8
    RCH4ProdHeter(NGL,K)    = AZMAX1(RespGrossHeter(NGL,K)-RespGrossHeter(NGL,K))
    RO2Uptk4RespHeter(NGL,K)= RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)     = 0.0_r8
  ENDDO
  end associate
  end subroutine AcetoMethanogenCatabolism

!------------------------------------------------------------------------------------------

  subroutine H2MethanogensCatabolism(I,J,N,RMOMK,TOMEAutoKC,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !Hydrogenotrophic CH4 production
  !CO2 + 4H2 -> CH4 + 2H2O
  !H2 is produced from fermentation
  !CO2+0.667H2-> CH4+ 3H2O
  !use CO2 for both energy generation and C biomass
  implicit none
  integer, intent(in) :: I,J, N
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(in) :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='H2MethanogensCatabolism'
  real(r8) :: GH2X,GH2H
  real(r8) :: H2GSX
  real(r8) :: VMAX
  REAL(R8) :: XCO2
  real(r8) :: RGOMP,RVOXP,GH2C
  real(r8) :: ECH2         !Efficiency of converting CO2 into biomass (CH2O) by hydrogenotrophic methanogen
  real(r8), parameter :: GCHA=38.9/12._r8 !Gibbs free energy for anabolic reaction, CO2(aq)+2H2(aq)->CH2O + H2O,[kJ (gC)-1]

  integer  :: NGL

  associate(                                                &
    GrowthEnvScalAutor     => nmics%GrowthEnvScalAutor,     &
    FBiomNutStoiScalAutor  => nmics%FBiomNutStoiScalAutor,  &
    FSBSTAutor             => nmicdiag%FSBSTAutor,          &
    OMActAutor             => nmics%OMActAutor,             &
    RO2Dmnd4GrossRespAutor => nmicf%RO2Dmnd4GrossRespAutor, &
    RO2Uptk4RespAutor      => nmicf%RO2Uptk4RespAutor,      &
    RespGrossAutor         => nmicf%RespGrossAutor,         &
    RO2UptkAutor           => nmicf%RO2UptkAutor,           &
    RCO2ProdAutor          => nmicf%RCO2ProdAutor,          &
    RCH4ProdAutor          => nmicf%RCH4ProdAutor,          &
    ECHZAutor              => nmicf%ECHZAutor,              &
    JGniA                  => micpar%JGniA,                 &
    JGnfA                  => micpar%JGnfA,                 &
    TKS                    => micfor%TKS,                   &
    CH2GS                  => micstt%CH2GS,                 &
    CCO2S                  => micstt%CCO2S,                 &
    H2GS                   => micstt%H2GS,                  &
    RH2UptkAutor           => nmicdiag%RH2UptkAutor,        &
    RO2MetaDmndAutor       => micflx%RO2MetaDmndAutor       &
  )
  !     begin_execution
  !
  !     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
  !     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
  !
  !     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
  !     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
  !     VMAX=substrate-unlimited H2 oxidation rate
  !     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (tCResp4H2Prod)
  !     CH2GS=H2 concentration, H2KM=Km for H2 uptake
  !     RGOMP=H2 oxidation, ROXY*=O2 demand
  !
  !     and energy yield of hydrogenotrophic
  !     methanogenesis GH2X at ambient H2 concentration CH2GS
  !     CCO2S=aqueous CO2 concentration
  !
  call PrintInfo('beg '//subname)
  XCO2         = CCO2S/(CCO2S+CCKM)
  RH2UptkAutor = 0.0_r8
  DO NGL=micpar%JGniA(N),micpar%JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle
    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaint(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    !Use catabolic reaction: CO2(aq)+4H2(aq) -> CH4 + 2H2O, 8/12=0.667, 1.5=12/8,
    !to drive anabolic reaction: CO2(aq)+2H2(aq) -> CH2O + H2O,

    GH2X = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**2)
    GH2C = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**4)/12._r8
    GH2H = GH2X/12.0_r8
    ECH2 = AMIN1(1._r8,AZMAX1((GCOX+GH2C)/(GCHA+GH2H)))
    !biomass yield as measured based on C, using respiration CH2O +2H2 -> CH4 + H2O for energy
    ECHZAutor(NGL)  = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCOX+GH2H))/EOMH)))
    VMAX            = OMActAutor(NGL)*VMXCH4gH2*GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*XCO2
    !0.111 is the stoichiometry from fermentation, C6H12O6 + 2H2O-> 2(C2H4O2)+ 4H2 + 2CO2, 8/72=0.111
    H2GSX           = AZMAX1(H2GS+0.111_r8*naqfdiag%tCResp4H2Prod)
    FSBSTAutor(NGL) = CH2GS/(CH2GS+H2KM)

    !first CO2 is partitioned into CH4 (catabolic) and CH2O (anabolic for respiration)
    !CO2+2H2 -> CH2O+H2O, 3CO2+H2->3CH2O+H2O potential C for respiration
    !CH2O +2H2-> CH4 + H2O
    !RGOMP is based on H2-driven methanogen respiration, which is used to support growth + (growth/maint resp)
    !assuming all electrons/reducing power are produced during catabolic reaction, so no more H2 is needed for biomass growth computed with RGOMP
    !H2 uptake rate
    RVOXP = AMIN1(1.5_r8*H2GSX/(1._r8+0.5_r8*ECH2*ECHZAutor(NGL)),VMAX*FSBSTAutor(NGL))
    RGOMP = RVOXP*ECH2*ECHZAutor(NGL)

    RO2Dmnd4GrossRespAutor(NGL) = 0.0_r8
    RO2MetaDmndAutor(NGL)       = 0.0_r8

    !obtains CO2 uptake for energy generation
    RespGrossAutor(NGL)    = RGOMP
    RCH4ProdAutor(NGL)     = RVOXP+RGOMP
    naqfdiag%tCH4ProdH2    = naqfdiag%tCH4ProdH2+RCH4ProdAutor(NGL)
    RO2Uptk4RespAutor(NGL) = 0._r8
    RH2UptkAutor           = RH2UptkAutor+0.667_r8*RCH4ProdAutor(NGL)
    RO2UptkAutor(NGL)      = 0.0_r8
  ENDDO
!
  call PrintInfo('end '//subname)
  end associate
  end subroutine H2MethanogensCatabolism
end module MethanogenMod
