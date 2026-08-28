module AerobicBacteriaMod
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
  public :: AerobicHeteroBactCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AerobicHeteroBactCatabolism(I,J,N,K,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  !
  !Description
  !catabolism of aerobic heterotrophs
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
  character(len=*), parameter :: subname='AerobicHeteroBactCatabolism'
  real(r8) :: WatStressMicb  !moisture sensivity of microbial activity
  integer :: NGL
  real(r8) :: RGOMP  !total DOC/acetate C uptake for potential respiraiton
  real(r8) :: EO2Q   !respiraiton efficiency, i.e. fraction of 1 gC (uptake) used for respiraiton
  real(r8) :: OXKX,FOXYX
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX

!     begin_execution
  associate(                                                               &
    OMActHeter             => nmics%OMActHeter,                            &
    FBiomStoiScalarHeter   => nmics%FBiomStoiScalarHeter,                  &
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,                    &
    WSensGroHeter          => nmics%WSensGroHeter,                         &
    TSensGroHeter          => nmics%TSensGroHeter,                         &
    TempMaintRHeter        => nmics%TempMaintRHeter,                       &
    FracOMActHeter         => nmics%FracOMActHeter,                        &
    FracNO2ReduxHeter      => nmics%FracNO2ReduxHeter,                     &
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,                   &
    FSBSTHeter             => nmicdiag%FSBSTHeter,                         &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,                     &
    OxyLimterHeter         => nmics%OxyLimterHeter,                        &
    RO2DmndHeter           => nmicf%RO2DmndHeter,                          &
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,                   &
    ECHZHeter              => nmicf%ECHZHeter,                             &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,                         &
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,                         &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,                     &
    RH2ProdHeter           => nmicf%RH2ProdHeter,                          &
    RGOCP                  => nmicf%RGOCP,                                 &
    RGOAP                  => nmicf%RGOAP,                                 &
    FOQC                   => nmicf%FOQC,                                  &
    FOQA                   => nmicf%FOQA,                                  &
    FGOCP                  => nmicf%FGOCP,                                 &
    RespGrossHeter         => nmicf%RespGrossHeter,                        &
    RAcetateProdHeter      => nmicf%RAcetateProdHeter,                     &
    FGOAP                  => nmicf%FGOAP,                                 &
    ZEROS                  => micfor%ZEROS,                                &
    RO2EcoDmndPrev         => micfor%RO2EcoDmndPrev,                       &
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,                      &
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,                  &
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,                    &
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,                &
    PSISoilMatricP         => micfor%PSISoilMatricP,                       &
    DOM                    => micstt%DOM,                                  &
    mid_HeterAerobBacter   => micpar%mid_HeterAerobBacter,                 &
    mid_Facult_DenitBacter => micpar%mid_Facult_DenitBacter,               &
    mid_HeterAerobN2Fixer  => micpar%mid_HeterAerobN2Fixer,                &
    TOMEK                  => nmicdiag%TOMEK,                              &
    TSensGrowth            => nmicdiag%TSensGrowth,                        &
    TSensMaintR            => nmicdiag%TSensMaintR,                        &
    TotBiomNO2Consumers    => nmicdiag%TotBiomNO2Consumers,                &
    RO2DmndHetert          => micflx%RO2DmndHetert,                        &
    RDOCUptkHeter          => micflx%RDOCUptkHeter,                        &
    RAcetateUptkHeter      => micflx%RAcetateUptkHeter,                    &
    RO2DmndHetertPrev      => micflx%RO2DmndHetertPrev,                    &
    tRespGrossHeterUlm     => naqfdiag%tRespGrossHeterUlm,                 &
    TotActMicrobiom        => nmicdiag%TotActMicrobiom,                    &
    tRGOXP                 => micflx%tRGOXP,                               &
    tRGOZP                 => micflx%tRGOZP,                               &
    FOCA                   => ncplxs%FOCA,                                 &
    FOAA                   => ncplxs%FOAA,                                 &
    CDOM                   => ncplxs%CDOM                                  &
  )
  call PrintInfo('beg '//subname)
  !     ENERGY YIELDS OF O2 REDOX REACTIONS
  !     E* = growth respiration efficiency calculated in PARAMETERS
  !
  !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE..0_r8)cycle

    !prepare trait parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA(NGL,K),micfor,naqfdiag,nmicdiag,nmics)

    OXKX  = OXKM
    IF(RO2EcoDmndPrev.GT.ZEROS)THEN
      FOXYX=AMAX1(FMN,RO2DmndHetertPrev(NGL,K)/RO2EcoDmndPrev)
    ELSE
      FOXYX=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    naqfdiag%TFOXYX = naqfdiag%TFOXYX+FOXYX

    ! N=OBLIGATE AEROBIC bacteria,
    IF(N.EQ.mid_HeterAerobBacter)THEN
      EO2Q=EO2X
    ! FACULTATIVE ANAEROBES,
    ELSEIF(N.EQ.mid_Facult_DenitBacter)THEN
      EO2Q=EO2D
    !aerobic N2 FIXERS
    ELSEIF(N.EQ.mid_HeterAerobN2Fixer)THEN
      EO2Q=ENFX
    ENDIF
    !
    ! O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
    ! 'RGO*Z'FROM SPECIFIC RESPIRATION RATE, ACTIVE BIOMASS, DOC OR
    ! ACETATE CONCENTRATION,MICROBIAL C:N:P FACTOR, AND TEMPERATURE
    ! FOLLOWED BY POTENTIAL RESPIRATION RATES 'RGO*P' WITH UNLIMITED
    ! SUBSTRATE USED FOR MICROBIAL COMPETITION FACTOR

    ! COQC,COQA=DOC,DOA concentration, FOCA,FOAA=DOC,DOA vs DOC+DOA
    ! FBiomStoiScalarHeter=N,P limitation,VMXO=specific respiration rate
    ! WatStressMicb=water stress effect, OMA=active biomass
    ! TSensGrowth=temp stress effect,FOQC,FOQA=OQC,OQA limitation
    ! RGOMP=O2-unlimited respiration of DOC+DOA
    ! RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
    !
    FSBSTC            = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
    FSBSTA            = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
    FSBSTHeter(NGL,K) = FOCA(K)*FSBSTC+FOAA(K)*FSBSTA

    RGOCY  = AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))*VMXO*GrowthEnvScalHeter(NGL,K)
    RGOCZ  = RGOCY*FSBSTC*FOCA(K) !MM uptake of DOC
    RGOAZ  = RGOCY*FSBSTA*FOAA(K) !MM uptake of acetate

    !obtain kinetically unlimited DOM/acetate uptake
    RGOCX = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*EO2Q)             !DOC respiration
    RGOAX = AZMAX1(DOM(idom_acetate,K)*FOQA(NGL,K)*EO2A)         !acetate respiraiton

    !obtain the final uptake
    RGOCP(NGL,K) = AMIN1(RGOCX,RGOCZ)      !DOC respiraiton
    RGOAP(NGL,K) = AMIN1(RGOAX,RGOAZ)      !acetate respiraiton
    RGOMP        = RGOCP(NGL,K)+RGOAP(NGL,K)      !total C respiration before O2 limitation

    tRGOXP = tRGOXP+RGOCX+RGOAX
    tRGOZP = tRGOZP+RGOCZ+RGOAZ            !potential C oxidation without C and O2 limitation
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
    RO2Dmnd4RespHeter(NGL,K) = 2.667_r8*RGOMP                 !O2 demand
    RO2DmndHeter(NGL,K)      = RO2Dmnd4RespHeter(NGL,K)

    !make a copy for flux limiter
    RO2DmndHetert(NGL,K)       = RO2DmndHeter(NGL,K)
    RDOCUptkHeter(NGL,K)       = RGOCZ               !potential DOC (O2-unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = RGOAZ               !potential acetate (O2-unlimited) uptake flux

    call AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)

    ROQC4HeterMicrobAct(NGL,K) = RGOCY*OxyLimterHeter(NGL,K)  !C demand for oxidation
    RespGrossHeter(NGL,K)      = RGOMP*OxyLimterHeter(NGL,K)  !actual respiration O2-limited
    RCO2ProdHeter(NGL,K)       = RespGrossHeter(NGL,K)
    RAcetateProdHeter(NGL,K)     = 0.0_r8
    RCH4ProdHeter(NGL,K)       = 0.0_r8
    RO2Uptk4RespHeter(NGL,K)   = RO2Dmnd4RespHeter(NGL,K)*OxyLimterHeter(NGL,K)
    RH2ProdHeter(NGL,K)        = 0.0_r8
    tRespGrossHeterUlm         = tRespGrossHeterUlm+RGOMP
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AerobicHeteroBactCatabolism

end module AerobicBacteriaMod
