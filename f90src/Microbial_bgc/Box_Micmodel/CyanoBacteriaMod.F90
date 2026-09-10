module CyanoBacteriaMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use abortutils,           only: endrun,   destroy
  use EcoSIMCtrlMod,        only: etimer
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: DebugPrint,PrintInfo
  use minimathmod
  use ElmIDMod
  use TracerIDMod
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: AerobicHeterO2Uptake, CalcRespMaintHeter, &
                                  StageFuncGuild
  implicit none

  private
  public :: CyanoBacteriaCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains
!------------------------------------------------------------------------------------------
  subroutine CyanoBacteriaCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  !

  implicit none
  integer, intent(in) :: I,J,N,K
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='CyanoBacteriaCatabolism'

  real(r8), parameter :: K_low_photo= 6.e-3_r8   ![h-1]
  integer :: NGL
  real(r8) :: EO2Q
  real(r8) :: fBIOM
  real(r8) :: fPAR,xCO2,r_photo,f_low_photo
  real(r8) :: RPhotoResp
  real(r8) :: RGOMP,RHeterRespUlm
  
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX
  real(r8) :: OXKX,FOXYX
  real(r8) :: WatStressMicb
  real(r8) :: f_cyno_heter
  associate(                                                               &
    OMActHeter             => nmics%OMActHeter,                            &
    FBiomStoiScalarHeter   => nmics%FBiomStoiScalarHeter,                  &
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,                    &
    WSensGroHeter          => nmics%WSensGroHeter,                         &
    TSensGroHeter          => nmics%TSensGroHeter,                         &
    TempMaintRHeter        => nmics%TempMaintRHeter,                       &
    FracOMActHeter         => nmics%FracOMActHeter,                        &
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,                   &
    FSBSTHeter             => nmicdiag%FSBSTHeter,                         &
    TotActMicrobiom        => nmicdiag%TotActMicrobiom,                    &
    TOMEK                  => nmicdiag%TOMEK,                              &
    TSensGrowth            => nmicdiag%TSensGrowth,                        &
    TSensMaintR            => nmicdiag%TSensMaintR,                        &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,                     &
    OxyLimterHeter         => nmics%OxyLimterHeter,                        &
    RO2DmndHeter           => nmicf%RO2DmndHeter,                          &
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,                   &
    ECHZHeter              => nmicf%ECHZHeter,                             &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,                         &
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,                         &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,                     &
    RCO2FixCyano           => nmicf%RCO2FixCyano   ,                       &
    RH2ProdHeter           => nmicf%RH2ProdHeter,                          &
    RGOCP                  => nmicf%RGOCP,                                 &
    RGOAP                  => nmicf%RGOAP,                                 &
    FOQC                   => nmicf%FOQC,                                  &
    FOQA                   => nmicf%FOQA,                                  &
    FGOCP                  => nmicf%FGOCP,                                 &
    FGOAP                  => nmicf%FGOAP,                                 &
    fPhotoR                => nmicf%fPhotoR,                               &
    RespGrossHeter         => nmicf%RespGrossHeter,                        &
    RAcetateProdHeter      => nmicf%RAcetateProdHeter,                     &
    RO2DmndHetert          => micflx%RO2DmndHetert,                        &
    RDOCUptkHeter          => micflx%RDOCUptkHeter,                        &
    RAcetateUptkHeter      => micflx%RAcetateUptkHeter,                    &
    RMaintRespHeter        => nmicf%RMaintRespHeter,                       &    
    tRGOXP                 => micflx%tRGOXP,                               &
    tRGOZP                 => micflx%tRGOZP,                               &
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,                      &
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,                  &
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,                    &
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,                &
    RO2EcoDmndPrev         => micfor%RO2EcoDmndPrev,                       &
    RO2DmndHetertPrev      => micflx%RO2DmndHetertPrev,                    &
    PSISoilMatricP         => micfor%PSISoilMatricP,                       &
    PAR_rad                => micfor%PAR_rad,                              &
    ZEROS                  => micfor%ZEROS,                                &
    CCO2S                  => micstt%CCO2S,                                &
    DOM                    => micstt%DOM,                                  &
    FOCA                   => ncplxs%FOCA,                                 &
    FOAA                   => ncplxs%FOAA,                                 &
    CDOM                   => ncplxs%CDOM,                                 &
    tRespGrossHeterUlm     => naqfdiag%tRespGrossHeterUlm                  &
  )
  call PrintInfo('beg '//subname)

  !reduction of heterotrophic rate compared to actual heterotrophs
  f_cyno_heter=0.1_r8
  if(PAR_rad.LE.ZEROS)then
    f_cyno_heter=0.01_r8
  endif

  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)

    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle

    !prepare trait parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA(NGL,K),micfor,naqfdiag,nmicdiag,nmics)

    call CalcRespMaintHeter(NGL,K,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    OXKX  = OXKM
    IF(RO2EcoDmndPrev.GT.ZEROS)THEN
      FOXYX=AMAX1(FMN,RO2DmndHetertPrev(NGL,K)/RO2EcoDmndPrev)
    ELSE
      FOXYX=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    naqfdiag%TFOXYX = naqfdiag%TFOXYX+FOXYX

    EO2Q=ENFX
    fBIOM=AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))

    !photosynthesis
    !6CO2 + 6H2O + light -> (CH2O)6 + 6O2

    if(PAR_rad.GT.ZEROS)then

      fPAR=fPAR_func(PAR_rad)
      xCO2 = CCO2S/(CCO2S+CCKM)
      r_photo=fPAR*xCO2*GrowthEnvScalHeter(NGL,K)
      RCO2FixCyano(NGL,K) = r_photo * fBIOM

      !photosynthate respiraiton
      RPhotoResp=RCO2FixCyano(NGL,K)

      !(CH2O)6  + 6O2 -> 6CO2 +6 H2O
      f_low_photo = K_low_photo / (r_photo + K_low_photo)
    else
      f_low_photo         = 1._r8
      RCO2FixCyano(NGL,K) = 0._r8
      RPhotoResp          = 0._r8
    endif

    !grow on DOC and acetate
    FSBSTC            = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
    FSBSTA            = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
    FSBSTHeter(NGL,K) = FOCA(K)*FSBSTC+FOAA(K)*FSBSTA

    RGOCY  = fBIOM*VMXO*GrowthEnvScalHeter(NGL,K)*f_low_photo*f_cyno_heter
    RGOCZ  = RGOCY*FSBSTC*FOCA(K)  !MM uptake of DOC
    RGOAZ  = RGOCY*FSBSTA*FOAA(K)  !MM uptake of acetate

    !obtain kinetically unlimited DOM/acetate uptake
    RGOCX = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*EO2Q)      !DOC respiration
    RGOAX = AZMAX1(DOM(idom_acetate,K)*FOQA(NGL,K)*EO2A)  !acetate respiraiton

    !obtain the final uptake
    RGOCP(NGL,K) = AMIN1(RGOCX,RGOCZ)             !DOC respiraiton
    RGOAP(NGL,K) = AMIN1(RGOAX,RGOAZ)             !acetate respiraiton
    RHeterRespUlm = RGOCP(NGL,K)+RGOAP(NGL,K)
    RGOMP        = RPhotoResp+RHeterRespUlm      !total C respiration before O2 limitation

    tRGOXP = tRGOXP+RGOCX+RGOAX
    tRGOZP = tRGOZP+RGOCZ+RGOAZ            !potential C oxidation without C and O2 limitation
    IF(RGOMP.GT.ZEROS)THEN
      FGOCP(NGL,K) = RGOCP(NGL,K)/RGOMP
      FGOAP(NGL,K) = RGOAP(NGL,K)/RGOMP      
    ELSE
      FGOCP(NGL,K) = 1.0_r8
      FGOAP(NGL,K) = 0.0_r8      
    ENDIF

    RO2Dmnd4RespHeter(NGL,K) = 2.667_r8*(RHeterRespUlm-RCO2FixCyano(NGL,K))         !external O2 demand
    RO2DmndHeter(NGL,K)      = RO2Dmnd4RespHeter(NGL,K)
    ECHZHeter(NGL,K)         = EO2Q*FGOCP(NGL,K)+EO2A*FGOAP(NGL,K)

    !make a copy for flux limiter
    RO2DmndHetert(NGL,K)       = RO2DmndHeter(NGL,K)
    RDOCUptkHeter(NGL,K)       = RGOCZ               !potential DOC (O2-unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = RGOAZ               !potential acetate (O2-unlimited) uptake flux

    call AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)

    ROQC4HeterMicrobAct(NGL,K) = RGOCY*OxyLimterHeter(NGL,K) !C demand for oxidation
    RespGrossHeter(NGL,K)      = RPhotoResp+RHeterRespUlm*OxyLimterHeter(NGL,K)
    
    IF(RespGrossHeter(NGL,K).GT.0._R8)THEN
      fPhotoR(NGL,K)   = RPhotoResp/RespGrossHeter(NGL,K) 
    ELSE
      fPhotoR(NGL,K)   = 0._R8
    ENDIF  

    RCO2ProdHeter(NGL,K)       = RespGrossHeter(NGL,K)-RCO2FixCyano(NGL,K)
    RAcetateProdHeter(NGL,K)   = 0.0_r8
    RCH4ProdHeter(NGL,K)       = 0.0_r8
    RO2Uptk4RespHeter(NGL,K)   = RO2Dmnd4RespHeter(NGL,K)*OxyLimterHeter(NGL,K)
    RH2ProdHeter(NGL,K)        = 0.0_r8
    tRespGrossHeterUlm         = tRespGrossHeterUlm+RGOMP

  ENDDO

  call PrintInfo('end '//subname)
  end associate

  end subroutine CyanoBacteriaCatabolism
!------------------------------------------------------------------------------------------
  function fPAR_func(PAR_RAD)result(ans)
  implicit none
  real(r8), intent(in) :: PAR_RAD                ![umol photon m-2 s-1]
  real(r8), parameter :: alpha_cyno = 1.1e-4_r8  ! [h-1] / [umol photon m-2 s-1]
  real(r8), parameter :: beta_cyno  = 3.3e-4_r8  ! inverse [PAR]
  real(r8), parameter :: gamma_cyno = 1.5e-3_r8  ! inverse [PAR]

  real(r8) :: ans   !photosynthesis rate [h-1]


  ans = alpha_cyno * PAR_rad * (1._r8 - beta_cyno*PAR_rad) &
        / (1._r8 + gamma_cyno*PAR_rad)

  end function fPAR_func

  end module CyanoBacteriaMod
