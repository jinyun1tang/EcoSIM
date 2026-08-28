module FungiMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use minimathmod
  use ElmIDMod
  use TracerIDMod
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: AerobicHeterO2Uptake, StageFuncGuild
  implicit none

  private
  public :: AerobicFungiCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains
!------------------------------------------------------------------------------------------

  subroutine AerobicFungiCatabolism(I,J,N,K,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N,K
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_State_type),intent(inout):: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  real(r8) :: WatStressMicb  !moisture sensivity of microbial activity
  integer  :: NGL
  real(r8) :: EO2Q  !respiraiton efficiency
  real(r8) :: OXKX
  real(r8) :: RGOMP  !total DOC/acetate C uptake for potential respiraiton
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX,FOXYX

!     begin_execution
  associate(                                                               &
    OMActHeter             => nmics%OMActHeter,                            &
    FBiomStoiScalarHeter   => nmics%FBiomStoiScalarHeter,                  &
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,                    &
    WSensGroHeter          => nmics%WSensGroHeter,                         &
    TSensGroHeter          => nmics%TSensGroHeter,                         &
    TempMaintRHeter        => nmics%TempMaintRHeter,                       &
    FracOMActHeter         => nmics%FracOMActHeter,                        &
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,                   &
    OxyLimterHeter         => nmics%OxyLimterHeter,                        &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,                     &
    RO2DmndHeter           => nmicf%RO2DmndHeter,                          &
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,                   &
    ECHZHeter              => nmicf%ECHZHeter,                             &
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,                         &
    FGOCP                  => nmicf%FGOCP,                                 &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,                     &
    RespGrossHeter         => nmicf%RespGrossHeter,                        &
    FGOAP                  => nmicf%FGOAP,                                 &
    RH2ProdHeter           => nmicf%RH2ProdHeter,                          &
    FOQC                   => nmicf%FOQC,                                  &
    FOQA                   => nmicf%FOQA,                                  &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,                         &
    RGOCP                  => nmicf%RGOCP,                                 &
    RGOAP                  => nmicf%RGOAP,                                 &
    RAcetateProdHeter      => nmicf%RAcetateProdHeter,                     &
    RO2EcoDmndPrev         => micfor%RO2EcoDmndPrev,                       &
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,                      &
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,                  &
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,                    &
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,                &
    ZEROS                  => micfor%ZEROS,                                &
    PSISoilMatricP         => micfor%PSISoilMatricP,                       &
    DOM                    => micstt%DOM,                                  &
    FSBSTHeter             => nmicdiag%FSBSTHeter,                         &
    TOMEK                  => nmicdiag%TOMEK,                              &
    TSensGrowth            => nmicdiag%TSensGrowth,                        &
    TSensMaintR            => nmicdiag%TSensMaintR,                        &
    RO2DmndHetertPrev      => micflx%RO2DmndHetertPrev,                    &
    tRespGrossHeterUlm     => naqfdiag%tRespGrossHeterUlm,                 &
    RO2DmndHetert          => micflx%RO2DmndHetert,                        &
    RDOCUptkHeter          => micflx%RDOCUptkHeter,                        &
    TotActMicrobiom        => nmicdiag%TotActMicrobiom,                    &
    RAcetateUptkHeter      => micflx%RAcetateUptkHeter,                    &
    tRGOXP                 => micflx%tRGOXP,                               &
    tRGOZP                 => micflx%tRGOZP,                               &
    FOCA                   => ncplxs%FOCA,                                 &
    FOAA                   => ncplxs%FOAA,                                 &
    CDOM                   => ncplxs%CDOM                                  &
  )
      !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle

    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA(NGL,K),micfor,naqfdiag,nmicdiag,nmics)

    OXKX  = OXKM
    IF(RO2EcoDmndPrev.GT.ZEROS)THEN
      FOXYX=AMAX1(FMN,RO2DmndHetertPrev(NGL,K)/RO2EcoDmndPrev)
    ELSE
      FOXYX=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    naqfdiag%TFOXYX = naqfdiag%TFOXYX+FOXYX

    EO2Q=EO2G
    FSBSTC = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
    FSBSTA = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
    FSBSTHeter(NGL,K)  = FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
    RGOCY  = AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))*VMXO*GrowthEnvScalHeter(NGL,K)
    RGOCZ  = RGOCY*FSBSTC*FOCA(K)
    RGOAZ  = RGOCY*FSBSTA*FOAA(K)

    !obtain kinetically unlimited DOM/acetate uptake
    RGOCX = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*EO2Q)     !potential DOC respiraiton
    RGOAX = AZMAX1(DOM(idom_acetate,K)*FOQA(NGL,K)*EO2A)        !potential acetate respiraiton
    !obtain the final uptake
    RGOCP(NGL,K) = AMIN1(RGOCX,RGOCZ)       !potential DOC(-limited) respiraiton
    RGOAP(NGL,K) = AMIN1(RGOAX,RGOAZ)
    RGOMP        = RGOCP(NGL,K)+RGOAP(NGL,K)
    tRGOXP       = tRGOXP+RGOCX+RGOAX
    tRGOZP       = tRGOZP+RGOCZ+RGOAZ

    IF(RGOMP.GT.ZEROS)THEN
      FGOCP(NGL,K) = RGOCP(NGL,K)/RGOMP
      FGOAP(NGL,K) = RGOAP(NGL,K)/RGOMP
    ELSE
      FGOCP(NGL,K) = 1.0_r8
      FGOAP(NGL,K) = 0.0_r8
    ENDIF
    !
    ! ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
    ! BY HETEROTROPHIC AEROBES

    ! ECHZHeter=growth respiration yield, averaged over acetate and DOC/glucose
    ! RO2Dmnd4RespHeter,RO2DmndHeter,RO2DmndHetert=O2 demand from DOC,DOA oxidation
    ! RDOCUptkHeter,RAcetateUptkHeter=DOC,DOA demand from DOC,DOA oxidation
    ! ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
    ! CH2O+O2 -> CO2 + H2O, (32/12.=2.667)
    ECHZHeter(NGL,K)         = EO2Q*FGOCP(NGL,K)+EO2A*FGOAP(NGL,K)
    RO2Dmnd4RespHeter(NGL,K) = 2.667_r8*RGOMP
    RO2DmndHeter(NGL,K)      = RO2Dmnd4RespHeter(NGL,K)
    !make a copy for flux limiter
    RO2DmndHetert(NGL,K)       = RO2DmndHeter(NGL,K)
    RDOCUptkHeter(NGL,K)       = RGOCZ    !potential DOC (unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = RGOAZ

    call AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)

    ROQC4HeterMicrobAct(NGL,K) = RGOCY*OxyLimterHeter(NGL,K)
    RespGrossHeter(NGL,K)      = RGOMP*OxyLimterHeter(NGL,K)
    RCO2ProdHeter(NGL,K)       = RespGrossHeter(NGL,K)
    RAcetateProdHeter(NGL,K)     = 0.0_r8
    RCH4ProdHeter(NGL,K)       = 0.0_r8
    RO2Uptk4RespHeter(NGL,K)   = RO2Dmnd4RespHeter(NGL,K)*OxyLimterHeter(NGL,K)
    RH2ProdHeter(NGL,K)        = 0.0_r8
    tRespGrossHeterUlm         = tRespGrossHeterUlm+RGOMP
  ENDDO
  end associate
  end subroutine AerobicFungiCatabolism

end module FungiMod
