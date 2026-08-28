module FermenterMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use minimathmod
  use TracerIDMod
  use EcosimConst
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: StageFuncGuild

  implicit none

  private
  public :: AcetogFermentCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AcetogFermentCatabolism(N,K,micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx,nmicdiag)
  !
  !Description:
  !Fermentation and acetogenic N2 fixers
  !(CH2O)6 +2H2O-> 2CO2 + 2(CH2O)2 + 4H2, mole based
  !(CH2O)6 -> 2CO2 + 2/3 (CH2O)2 + 8/(72)H2, mass based
  !fermenters only take up DOC/glucose
  !it can be fermenters or anaerobic N2 fixers
  implicit none
  integer, intent(in) :: N,K

  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  real(r8) :: WatStressMicb
  real(r8) :: RGOMP   !potential respiration for metabolism
  integer  :: NGL
  real(r8) :: GH2X,GH2F
  real(r8) :: GOAX,GOAF
  real(r8) :: GHAX
  REAL(R8) :: oxyi
  real(r8) :: RGOFX,RGOFY,RGOFZ

  real(r8), parameter :: GlucoseC=72._r8  !glucose has 72 gC/mol

  ! begin_execution
  associate(                                            &
    FBiomStoiScalarHeter => nmics%FBiomStoiScalarHeter, &
    OMActHeter           => nmics%OMActHeter,           &
    FSBSTHeter           => nmicdiag%FSBSTHeter,        &
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,   &
    RO2Dmnd4RespHeter    => nmicf%RO2Dmnd4RespHeter,    &
    RO2DmndHeter         => nmicf%RO2DmndHeter,         &
    ECHZHeter            => nmicf%ECHZHeter,            &
    FOQC                 => nmicf%FOQC,                 &
    FOQA                 => nmicf%FOQA,                 &
    RCH4ProdHeter        => nmicf%RCH4ProdHeter,        &
    RO2Uptk4RespHeter    => nmicf%RO2Uptk4RespHeter,    &
    RCO2ProdHeter        => nmicf%RCO2ProdHeter,        &
    RespGrossHeter       => nmicf%RespGrossHeter,       &
    RH2ProdHeter         => nmicf%RH2ProdHeter,         &
    ROQC4HeterMicrobAct  => nmicf%ROQC4HeterMicrobAct,  &
    RAcetateProdHeter    => nmicf%RAcetateProdHeter,    &
    TotActMicrobiom      => nmicdiag%TotActMicrobiom,   &
    FGOCP                => nmicf%FGOCP,                &
    FGOAP                => nmicf%FGOAP,                &
    TKS                  => micfor%TKS,                 &
    PSISoilMatricP       => micfor%PSISoilMatricP,      &
    ZERO                 => micfor%ZERO,                &
    DOM                  => micstt%DOM,                 &
    CH2GS                => micstt%CH2GS,               &
    COXYS                => micstt%COXYS,               &
    RO2DmndHetert        => micflx%RO2DmndHetert,       &
    RDOCUptkHeter        => micflx%RDOCUptkHeter,       &
    RAcetateUptkHeter    => micflx%RAcetateUptkHeter,   &
    mid_fermentor        => micpar%mid_fermentor,       &
    CDOM                 => ncplxs%CDOM                 &
  )

  !
  !     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND
  !     ACETATE CONCENTRATION
  !
  !     GH2F=energy yield of acetotrophic methanogenesis per g C
  !     GHAX=H2 effect on energy yield of fermentation
  !     GOAX=acetate effect on energy yield of fermentation
  !     ECHZHeter=growth respiration efficiency of fermentation

  !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle
    !prepare trait parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA(NGL,K),micfor,naqfdiag,nmicdiag,nmics)

    GH2X = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**4)
    GH2F = GH2X/GlucoseC    !
    GOAX = RGASC*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI)**2)
    GOAF = GOAX/GlucoseC
    GHAX = GH2F+GOAF
    IF(N.EQ.mid_fermentor)THEN
      ECHZHeter(NGL,K)=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCHX-GHAX))/EOMF)))
    ELSE
        !dizotrophs, i.e. N2 fixers
      ECHZHeter(NGL,K)=AMAX1(ENFX,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCHX-GHAX))/EOMN)))
    ENDIF
    !
    !     RESPIRATION RATES BY HETEROTROPHIC ANAEROBES 'RGOMP' FROM
    !     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
    !     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
    !     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
    !     MICROBIAL COMPETITION FACTOR
    !
    !     OXYI=O2 inhibition of fermentation
    !     FBiomStoiScalarHeter=N,P limitation on respiration
    !     VMXF=maximum respiration rate by fermenters
    !     WatStressMicb=water stress effect on respiration
    !     OMA=active fermenter biomass
    !     TSensGrowth=temp stress effect, FOQC=OQC limitation
    !     RFOMP=O2-unlimited respiration of DOC
    !     ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
    !
    OXYI  = 1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*AMAX1(-COXYS+2.5_r8,-50._r8)))
    FSBSTHeter(NGL,K) = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)*OXYI
    RGOFY             = AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))*VMXF*GrowthEnvScalHeter(NGL,K)
    RGOFZ             = RGOFY*FSBSTHeter(NGL,K)
    RGOFX             = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*ECHZHeter(NGL,K))

    !potential respiration to expense
    RGOMP                      = AMIN1(RGOFX,RGOFZ)
    FGOCP(NGL,K)               = 1.0_r8
    FGOAP(NGL,K)               = 0.0_r8
    RO2Dmnd4RespHeter(NGL,K)   = 0.0_r8   !demand no oxygen
    RO2DmndHeter(NGL,K)        = 0.0_r8
    RO2DmndHetert(NGL,K)       = 0.0_r8
    RDOCUptkHeter(NGL,K)       = RGOFZ    !potential DOC (unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = 0.0_r8
    ROQC4HeterMicrobAct(NGL,K) = RGOFY*OXYI    !DOC/temperature-unlimited fermentation rate for microbial activity calculation
    naqfdiag%tCResp4H2Prod     = naqfdiag%tCResp4H2Prod+RGOMP

    !fermentation  (CH2O)6 -> 2CO2 + 2(CH2O)2
    RespGrossHeter(NGL,K)    = RGOMP
    RAcetateProdHeter(NGL,K)   = 0.667_r8*RespGrossHeter(NGL,K)
    RCO2ProdHeter(NGL,K)     = AZMAX1(RespGrossHeter(NGL,K)-RAcetateProdHeter(NGL,K))
    RCH4ProdHeter(NGL,K)     = 0.0_r8
    RO2Uptk4RespHeter(NGL,K) = RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)      = 0.111_r8*RespGrossHeter(NGL,K)

  ENDDO
!
  end associate
  end subroutine AcetogFermentCatabolism
end module FermenterMod
