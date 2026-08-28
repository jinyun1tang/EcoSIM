module MicAutoCPLXMod
! USES:
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use minimathmod,          only: safe_adb, AZMAX1, fixEXConsumpFlux, SubstrateDribbling
  use MicForcTypeMod,       only: micforctype
  use MicFluxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use MicrobeDiagTypes
  use ElmIDMod
  use TracerIDMod
  use EcoSIMSolverPar
  use EcoSimConst
  use NitroPars
  use MicrobMathFuncMod
  use MethanogenMod,        only: H2MethanogensCatabolism
  use MethanotrophMod,      only: AeroMethanotrophCatabolism, AMONC10Catabolism, &
                                  AMOANME2dCatabolism
  use NitrifierMod,         only: AmmoniaOxiDenitCatabolism, &
                                  AmmoniaOxidizerCatabolism, NitriteOxidizerCatabolism
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ActiveAutotrophs
  public :: AutotrophAnabolicUpdate
  contains

!------------------------------------------------------------------------------------------
  subroutine ActiveAutotrophs(I,J,N,SPOMK, RMOMK, &
    micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs,nmicdiag)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: N
  real(r8), intent(in) :: SPOMK(2)
  real(r8), intent(in)  :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='ActiveAutotrophs'
  integer  :: M
  real(r8) :: COMC
  real(r8) :: FOQC,FOQA
  real(r8) :: RGOCP
  real(r8) :: RGOMP
  real(r8) :: RVOXP
  real(r8) :: RVOXPA    !oxidation in soil, for NH3, NO2 and CH4
  real(r8) :: RVOXPB    !oxidation in band, specifically for NH3 and NO2
  real(r8) :: RGrowthRespAutor
  real(r8) :: RMaintDefcitcitAutor
  real(r8) :: RMaintRespAutor

! begin_execution
  associate(                                                         &
    TOMEAutoK                  => ncplxs%TOMEAutoK,                  &
    ZNH4T                      => nmicdiag%ZNH4T,                    &
    ZNO3T                      => nmicdiag%ZNO3T,                    &
    ZNO2T                      => nmicdiag%ZNO2T,                    &
    H2P4T                      => nmicdiag%H2P4T,                    &
    H1P4T                      => nmicdiag%H1P4T,                    &
    VOLWZ                      => nmicdiag%VOLWZ,                    &
    mid_AutoAmmoniaOxidBacter  => micpar%mid_AutoAmmoniaOxidBacter,  &
    mid_AutoNitriteOxidBacter  => micpar%mid_AutoNitriteOxidBacter,  &
    mid_AutoH2GenoCH4GenArchea => micpar%mid_AutoH2GenoCH4GenArchea, &
    mid_AutoAeroCH4OxiBacter   => micpar%mid_AutoAeroCH4OxiBacter,   &
    mid_AutoAMOANME2D          => micpar%mid_AutoAMOANME2D       ,   &
    mid_AutoAMONC10            => micpar%mid_AutoAMONC10           , &
    ZEROS                      => micfor%ZEROS,                      &
    SoilMicPMassLayer          => micfor%SoilMicPMassLayer,          &
    litrm                      => micfor%litrm,                      &
    VLSoilPoreMicP             => micfor%VLSoilPoreMicP              &
  )
  !
  !
  !     RESPIRATION RATES BY AUTOTROPHS 'RGOMP' FROM SPECIFIC
  !     OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
  !     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
  !     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
  !     MICROBIAL COMPETITION FACTOR. N=(1) NH4 OXIDIZERS (2) NO2
  !     OXIDIZERS,(3) CH4 OXIDIZERS, (5) H2TROPHIC METHANOGENS
  !
  !
  call PrintInfo('beg '//subname)
  if (N.eq.mid_AutoAmmoniaOxidBacter)then
    ! NH3 OXIDIZERS
    call AmmoniaOxidizerCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)

  elseif (N.eq.mid_AutoNitriteOxidBacter)then
    ! NO2 OXIDIZERS
    call NitriteOxidizerCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)

  elseif (N.eq.mid_AutoH2GenoCH4GenArchea)then
    ! H2TROPHIC METHANOGENS
    call H2MethanogensCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)

  elseif (N.eq.mid_AutoAeroCH4OxiBacter)then
    ! METHANOTROPHS
    call AeroMethanotrophCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  elseif (N.eq.mid_AutoAMONC10)then
    call AMONC10Catabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)    
  elseif (N.eq.mid_AutoAMOANME2D)then
    call AMOANME2dCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
    
  ENDIF

  IF(micpar%is_aerobic_autor(N))then
    call AerobicAutorO2Uptake(I,J,N,micfor,micstt,nmicf,nmics,micflx,naqfdiag)
  ENDif
  !
  !
  ! RO2UptkHeter, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
  ! RUPMX=O2-unlimited rate of O2 uptake
  ! FOXYX=fraction of O2 uptake by N,K relative to total
  ! dts_gas=1/(NPH*NPT)
  ! ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
  ! O2AquaDiffusvity=aqueous O2 diffusivity
  ! OXYG,OXYS=gaseous, aqueous O2 amounts
  ! Rain2LitRSurf,Irrig2LitRSurf_col=surface water flux from precipitation, irrigation
  ! O2_rain_conc,O2_irrig_conc=O2 concentration in Rain2LitRSurf,Irrig2LitRSurf_col
  !
  !
  !     AUTOTROPHIC DENITRIFICATION by ammonia-oxidizer
  !
  IF(N.EQ.mid_AutoAmmoniaOxidBacter .AND. (.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS))THEN  
    call AmmoniaOxiDenitCatabolism(I,J,N,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  ENDIF
  !
  !     BIOMASS DECOMPOSITION AND MINERALIZATION
  !
  ! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
  ! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
  !

  call BiomNutMinerMobilAutor(I,J,N,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,nmics,naqfdiag)
  !
  call GatherAutotrophRespiration(I,J,N,micfor,micflx,nmicf,nmics)
  !
  call GatherAutotrophAnabolicFlux(I,J,N,micflx,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  call PrintInfo('end '//subname)
  end associate
  end subroutine ActiveAutotrophs


!------------------------------------------------------------------------------------------

  subroutine SubstrateCompetAuto(NGL,N,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    micfor,naqfdiag,nmicf,nmics,micflx)
  !
  !Description:
  !Substrate competition for autotrophs
  !  
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(out):: FNH4X
  real(r8),intent(out) :: FNB3X,FNB4X,FNO3X
  real(r8),intent(out) :: FPO4X,FPOBX,FP14X,FP1BX
  type(micforctype), intent(in) :: micfor
  type(Cumlate_Flux_Diag_type),INTENT(INOUT)::  naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                                                   &
    FracOMActAutor          => nmics%FracOMActAutor,           &
    FracAutorBiomOfActK     => nmics%FracAutorBiomOfActK,      &
    AttenfNH4Autor          => micflx%AttenfNH4Autor,          &
    AttenfNO3Autor          => micflx%AttenfNO3Autor,          &
    AttenfH1PO4Autor        => micflx%AttenfH1PO4Autor,        &
    AttenfH2PO4Autor        => micflx%AttenfH2PO4Autor,        &
    RNH4EcoDmndSoilPrev     => micfor%RNH4EcoDmndSoilPrev,     &
    RNH4EcoDmndBandPrev     => micfor%RNH4EcoDmndBandPrev,     &
    RNO3EcoDmndSoilPrev     => micfor%RNO3EcoDmndSoilPrev,     &
    RNO3EcoDmndBandPrev     => micfor%RNO3EcoDmndBandPrev,     &
    RH2PO4EcoDmndSoilPrev   => micfor%RH2PO4EcoDmndSoilPrev,   &
    RH2PO4EcoDmndBandPrev   => micfor%RH2PO4EcoDmndBandPrev,   &
    RH1PO4EcoDmndSoilPrev   => micfor%RH1PO4EcoDmndSoilPrev,   &
    RH1PO4EcoDmndBandPrev   => micfor%RH1PO4EcoDmndBandPrev,   &
    RNH4EcoDmndLitrPrev     => micfor%RNH4EcoDmndLitrPrev,     &
    RNO3EcoDmndLitrPrev     => micfor%RNO3EcoDmndLitrPrev,     &
    RH1PO4EcoDmndLitrPrev   => micfor%RH1PO4EcoDmndLitrPrev,   &
    RH2PO4EcoDmndLitrPrev   => micfor%RH2PO4EcoDmndLitrPrev,   &
    RDOMEcoDmndPrev         => micfor%RDOMEcoDmndPrev,         &
    RAcetateEcoDmndPrev     => micfor%RAcetateEcoDmndPrev,     &
    SoilMicPMassLayer0      => micfor%SoilMicPMassLayer0,      &
    Lsurf                   => micfor%Lsurf,                   &
    litrm                   => micfor%litrm,                   &
    ZEROS                   => micfor%ZEROS,                   &
    VLNO3                   => micfor%VLNO3,                   &
    VLNOB                   => micfor%VLNOB,                   &
    VLPO4                   => micfor%VLPO4,                   &
    VLPOB                   => micfor%VLPOB,                   &
    VLNHB                   => micfor%VLNHB,                   &
    VLNH4                   => micfor%VLNH4,                   &
    RH1PO4UptkLitrAutorPrev => micflx%RH1PO4UptkLitrAutorPrev, &
    RH2PO4UptkLitrAutorPrev => micflx%RH2PO4UptkLitrAutorPrev, &
    RNO3UptkLitrAutorPrev   => micflx%RNO3UptkLitrAutorPrev,   &
    RNH4UptkLitrAutorPrev   => micflx%RNH4UptkLitrAutorPrev,   &
    RH1PO4UptkBandAutorPrev => micflx%RH1PO4UptkBandAutorPrev, &
    RH1PO4UptkSoilAutorPrev => micflx%RH1PO4UptkSoilAutorPrev, &
    RH2PO4UptkBandAutorPrev => micflx%RH2PO4UptkBandAutorPrev, &
    RH2PO4UptkSoilAutorPrev => micflx%RH2PO4UptkSoilAutorPrev, &
    RNO3UptkBandAutorPrev   => micflx%RNO3UptkBandAutorPrev,   &
    RNO3UptkSoilAutorPrev   => micflx%RNO3UptkSoilAutorPrev,   &
    RNH4UptkBandAutorPrev   => micflx%RNH4UptkBandAutorPrev,   &
    RNH4UptkSoilAutorPrev   => micflx%RNH4UptkSoilAutorPrev    &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4UptkSoilAutorPrev(NGL)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4UptkBandAutorPrev(NGL)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNHB)
  ENDIF
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3UptkSoilAutorPrev(NGL)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3UptkBandAutorPrev(NGL)/RNO3EcoDmndBandPrev)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNOB)
  ENDIF
  IF(RH2PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4UptkSoilAutorPrev(NGL)/RH2PO4EcoDmndSoilPrev)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4UptkBandAutorPrev(NGL)/RH2PO4EcoDmndBandPrev)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4UptkSoilAutorPrev(NGL)/RH1PO4EcoDmndSoilPrev)
  ELSE
    FP14X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4UptkBandAutorPrev(NGL)/RH1PO4EcoDmndBandPrev)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF

!  naqfdiag%TFNH4X = naqfdiag%TFNH4X+FNH4X
!  naqfdiag%TFNO3X = naqfdiag%TFNO3X+FNO3X
!  naqfdiag%TFPO4X = naqfdiag%TFPO4X+FPO4X
!  naqfdiag%TFP14X = naqfdiag%TFP14X+FP14X
!  naqfdiag%TFNH4B = naqfdiag%TFNH4B+FNB4X
!  naqfdiag%TFNO3B = naqfdiag%TFNO3B+FNB3X
!  naqfdiag%TFPO4B = naqfdiag%TFPO4B+FPOBX
!  naqfdiag%TFP14B = naqfdiag%TFP14B+FP1BX
  !
  ! FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
  ! MICROBIAL POPULATIONS IN SURFACE RESIDUE
  ! F*=fraction of substrate uptake relative to total uptake from
  ! previous hour in surface litter, labels as for soil layers above
  !
  !litter layer
  IF(litrm)THEN
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNH4Autor(NGL)=AMAX1(FMN,RNH4UptkLitrAutorPrev(NGL)/RNH4EcoDmndLitrPrev)
    ELSE
      AttenfNH4Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNO3Autor(NGL)=AMAX1(FMN,RNO3UptkLitrAutorPrev(NGL)/RNO3EcoDmndLitrPrev)
    ELSE
      AttenfNO3Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH2PO4Autor(NGL)=AMAX1(FMN,RH2PO4UptkLitrAutorPrev(NGL)/RH2PO4EcoDmndLitrPrev)
    ELSE
      AttenfH2PO4Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH1PO4Autor(NGL)=AMAX1(FMN,RH1PO4UptkLitrAutorPrev(NGL)/RH1PO4EcoDmndLitrPrev)
    ELSE
      AttenfH1PO4Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
  ENDIF
  !top soil layer
  !diagnostics off
!  IF(Lsurf.AND.SoilMicPMassLayer0.GT.ZEROS)THEN
!    naqfdiag%TFNH4X=naqfdiag%TFNH4X+micfor%AttenfNH4AutorR(NGL)
!    naqfdiag%TFNO3X=naqfdiag%TFNO3X+micfor%AttenfNO3AutorR(NGL)
!    naqfdiag%TFPO4X=naqfdiag%TFPO4X+micfor%AttenfH2PO4AutorR(NGL)
!    naqfdiag%TFP14X=naqfdiag%TFP14X+micfor%AttenfH1PO4AutorR(NGL)
!  ENDIF
  end associate
  end subroutine SubstrateCompetAuto
!------------------------------------------------------------------------------------------

  subroutine GatherAutotrophAnabolicFlux(I,J,N,micflx,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: N
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx  
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs

  character(len=*), parameter :: subname='GatherAutotrophAnabolicFlux'
  integer :: M,K,MID3,MID,MID1,NE,idom,NGL
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: CXC,RCCE(NumPlantChemElms)
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FRM
!     begin_execution
  associate(                                                                    &
    rCNBiomeActAutor                 => nmics%rCNBiomeActAutor,                 &
    GrowthEnvScalAutor               => nmics%GrowthEnvScalAutor,               &
    OMActAutor                       => nmics%OMActAutor,                       &
    DOMuptk4GrothAutor               => nmicf%DOMuptk4GrothAutor,               &
    NonstX2stBiomAutor               => nmicf%NonstX2stBiomAutor,               &
    RespGrossAutor                   => nmicf%RespGrossAutor,                   &
    RNOxReduxRespAutorLim            => nmicf%RNOxReduxRespAutorLim,            & !respiration energy due to NO2 
    RMaintDmndAutor                  => nmicf%RMaintDmndAutor,                  &
    RkillLitfalOMAutor               => nmicf%RkillLitfalOMAutor,               &
    RkillLitrfal2HumOMAutor          => nmicf%RkillLitrfal2HumOMAutor,          &
    RkillLitrfal2ResduOMAutor        => nmicf%RkillLitrfal2ResduOMAutor,        &
    RMaintDefcitLitrfalOMAutor       => nmicf%RMaintDefcitLitrfalOMAutor,       &
    RMaintDefLitrfal2HumOMAutor      => nmicf%RMaintDefLitrfal2HumOMAutor,      &
    RMaintDefLitrfal2ResduOMAutor    => nmicf%RMaintDefLitrfal2ResduOMAutor,    &
    RKillOMAutor                     => nmicf%RKillOMAutor,                     &
    RkillRecycOMAutor                => nmicf%RkillRecycOMAutor,                &
    RMaintDefcitKillOMAutor          => nmicf%RMaintDefcitKillOMAutor,          &
    RMaintDefcitRecycOMAutor         => nmicf%RMaintDefcitRecycOMAutor,         &
    Resp4NFixAutor                   => nmicf%Resp4NFixAutor,                   &
    ECHZAutor                        => nmicf%ECHZAutor,                        & !respiraiton ratio
    rNCOMCAutor                      => micpar%rNCOMCAutor,                     &
    rPCOMCAutor                      => micpar%rPCOMCAutor,                     &
    JGniA                            => micpar%JGniA,                           &
    JGnfA                            => micpar%JGnfA,                           &
    FL                               => micpar%FL,                              &
    ZEROS                            => micfor%ZEROS,                           &
    ZERO                             => micfor%ZERO,                            &
    RGrowthRespAutor                 => micflx%RGrowthRespAutor,                &
    RMaintDefcitcitAutor             => micflx%RMaintDefcitcitAutor,            &
    RMaintRespAutor                  => micflx%RMaintRespAutor,                 &        
    mBiomeAutor                      => micstt%mBiomeAutor,                     &
    EHUM                             => micstt%EHUM                             &
  )
  call PrintInfo('beg '//subname)
  !     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
  !     FROM O2, NOX AND C REDUCTION
  !
  !     CGOMX=DOC+acetate uptake from aerobic growth respiration
  !     CGOMD=DOC+acetate uptake from denitrifier growth respiration
  !     RMaintRespAutor=maintenance respiration
  !     RespGrossHeter=total respiration
  !     RNOxReduxRespDenitLim=respiration for denitrifcation
  !     Resp4NFixHeter=respiration for N2 fixation
  !     ECHZ,ENOX=growth respiration efficiencies for CO2, O2, and NOx reduction
  !     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake heterotrophs
  !     CGOMC=total CO2,CH4 uptake (autotrophs)
  !     CGOMN,CGOMP=DON, DOP uptake
  !     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
  !     OQN,OPQ=DON,DOP
  !     FOMK=faction of OMA in total OMA
  !     FCN,FCP=limitation from N,P
  !
  DO NGL=JGniA(N),JGnfA(N)  
    IF(OMActAutor(NGL).LE.0.0_r8)cycle        

    !potential growth respiraiton-respiraiton for N2-fixation 
    DOMuptk4GrothAutor(idom_beg:idom_end,NGL)=0._r8
    CGOMX = AMIN1(RMaintRespAutor(NGL),RespGrossAutor(NGL))+Resp4NFixAutor(NGL)+(RGrowthRespAutor(NGL)-Resp4NFixAutor(NGL))/ECHZAutor(NGL)
    CGOMD = RNOxReduxRespAutorLim(NGL)/ENOX         !CO2 synthesis due to NO2(-) reduction by NH3

    !total C uptake, which could be CO2, or CH4, depending on the type of organism
    !for aerobic methanotrophs, the following equals to CH4 uptake for maintenance+growth respiraiton+biomass
    DOMuptk4GrothAutor(ielmc,NGL)=CGOMX+CGOMD

    !
    !     TRANSFER UPTAKEN C,N,P FROM STORAGE TO ACTIVE BIOMASS
    !
    !     OMC,OMN,OMP=nonstructural C,N,P
    !     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
    !     rNCOMC,rPCOMC=maximum microbial N:C, P:C ratios
    !     RCCC,RCCN,RCCP=C,N,P recycling fractions
    !     RCCZ,RCCY=min, max C recycling fractions
    !     RCCX,RCCQ=max N,P recycling fractions
    !
    MID1=micpar%get_micb_id(ibiom_kinetic,NGL);MID3=micpar%get_micb_id(ibiom_reserve,NGL)
    IF(mBiomeAutor(ielmc,MID3).GT.ZEROS.AND.mBiomeAutor(ielmc,MID1).GT.ZEROS)THEN
      CCC=AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeAutor(ielmn,MID3)/(mBiomeAutor(ielmn,MID3)+mBiomeAutor(ielmc,MID3)*rNCOMCAutor(ibiom_reserve,NGL)) &
        ,mBiomeAutor(ielmp,MID3)/(mBiomeAutor(ielmp,MID3)+mBiomeAutor(ielmc,MID3)*rPCOMCAutor(ibiom_reserve,NGL))))
      CXC = mBiomeAutor(ielmc,MID3)/mBiomeAutor(ielmc,MID1)
      C3C = 1.0_r8/(1.0_r8+CXC/CKC)
      CNC = AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeAutor(ielmc,MID3)/(mBiomeAutor(ielmc,MID3)+mBiomeAutor(ielmn,MID3)/rNCOMCAutor(ibiom_reserve,NGL))))
      CPC=AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeAutor(ielmc,MID3)/(mBiomeAutor(ielmc,MID3)+mBiomeAutor(ielmp,MID3)/rPCOMCAutor(ibiom_reserve,NGL))))
      RCCC = RCCZ+AMAX1(CCC,C3C)*RCCY
      RCCN = CNC*RCCX
      RCCP = CPC*RCCQ
    ELSE
      RCCC = RCCZ
      RCCN = 0.0_r8
      RCCP = 0.0_r8
    ENDIF
    RCCE(ielmc)=RCCC
    RCCE(ielmn)=(RCCN+(1.0_r8-RCCN)*RCCC)
    RCCE(ielmp)=(RCCP+(1.0_r8-RCCP)*RCCC)    
    !
    !     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
    !
    !     CGOMZ=transfer from nonstructural to structural microbial C
    !     TFNG=temperature+water stress function
    !     OMGR=rate constant for transferring nonstructural to structural C
    !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
    !     FL=partitioning between labile and resistant microbial components
    !     OMC,OMN,OMP=nonstructural microbial C,N,P
    !
    MID3  = micpar%get_micb_id(ibiom_reserve,NGL)
    CGOMZ = GrowthEnvScalAutor(NGL)*OMGR*AZMAX1(mBiomeAutor(ielmc,MID3))

    DO M = 1, 2
      NonstX2stBiomAutor(ielmc,M,NGL)=FL(M)*CGOMZ
      IF(mBiomeAutor(ielmc,MID3).GT.ZEROS)THEN
        NonstX2stBiomAutor(ielmn,M,NGL)=AMIN1(FL(M)*AZMAX1(mBiomeAutor(ielmn,MID3)) &
          ,NonstX2stBiomAutor(ielmc,M,NGL)*mBiomeAutor(ielmn,MID3)/mBiomeAutor(ielmc,MID3))
        NonstX2stBiomAutor(ielmp,M,NGL)=AMIN1(FL(M)*AZMAX1(mBiomeAutor(ielmp,MID3)) &
          ,NonstX2stBiomAutor(ielmc,M,NGL)*mBiomeAutor(ielmp,MID3)/mBiomeAutor(ielmc,MID3))
      ELSE
        NonstX2stBiomAutor(ielmn,M,NGL)=0.0_r8
        NonstX2stBiomAutor(ielmp,M,NGL)=0.0_r8
      ENDIF

    !
    !     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
    !     RATE, TEMPERATURE
    !
    !     SPOMX=rate constant for microbial decomposition
    !     SPOMC=basal decomposition rate
    !     SPOMK=effect of low microbial C concentration on microbial decay
    !     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
    !     RDOMC,RDOMN,RDOMP=microbial C,N,P LitrFall
    !     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
    !
      MID   = micpar%get_micb_id(M,NGL)
      SPOMX = SQRT(GrowthEnvScalAutor(NGL))*SPOMC(M)*SPOMK(M)
      
      DO NE=1,NumPlantChemElms
        RKillOMAutor(NE,M,NGL)=AZMAX1(mBiomeAutor(NE,MID)*SPOMX)
            
        RkillRecycOMAutor(NE,M,NGL)=RKillOMAutor(NE,M,NGL)*RCCE(NE)

        RkillLitfalOMAutor(NE,M,NGL)=AZMAX1(RKillOMAutor(NE,M,NGL)-RkillRecycOMAutor(NE,M,NGL))
    !
    !     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
    !     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
    !
    !     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
    !     EHUM=humus transfer fraction from hour1.f
    !     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P LitrFall to residue
    !
        RkillLitrfal2HumOMAutor(NE,M,NGL)=RkillLitfalOMAutor(NE,M,NGL)*EHUM
    !
    !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
    !
        RkillLitrfal2ResduOMAutor(NE,M,NGL)=RkillLitfalOMAutor(NE,M,NGL)-RkillLitrfal2HumOMAutor(NE,M,NGL)
      ENDDO
      
    ENDDO
    
    !
    !     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
    !     EXCEEDS UPTAKE
    !
    !     OMC,OMN,OMP=microbial C,N,P
    !     RMaintRespAutor=total maintenance respiration
    !     RMaintDefcitcitAutor=senescence respiration
    !     RCCC=C recycling fraction
    !     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
    !     RMaintDmndHeter=maintenance respiration
    !     CNOMA,CPOMA=N:C,P:C ratios of active biomass
    !     RDMMC,RDMMN,RDMMP=microbial C,N,P LitrFall from senescence
    !     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
    !
    IF(RMaintDefcitcitAutor(NGL).GT.ZEROS.AND.RMaintRespAutor(NGL).GT.ZEROS.AND.RCCC.GT.ZERO)THEN
      FRM=RMaintDefcitcitAutor(NGL)/RMaintRespAutor(NGL)
      DO  M=1,2
        RMaintDefcitKillOMAutor(ielmc,M,NGL)=AMIN1(mBiomeAutor(ielmc,MID),AZMAX1(FRM*RMaintDmndAutor(M,NGL)/RCCC))
        RMaintDefcitKillOMAutor(ielmn,M,NGL)=AMIN1(mBiomeAutor(ielmn,MID),AZMAX1(RMaintDefcitKillOMAutor(ielmc,M,NGL)*rCNBiomeActAutor(ielmn,NGL)))
        RMaintDefcitKillOMAutor(ielmp,M,NGL)=AMIN1(mBiomeAutor(ielmp,MID),AZMAX1(RMaintDefcitKillOMAutor(ielmc,M,NGL)*rCNBiomeActAutor(ielmp,NGL)))
        DO NE=1,NumPlantChemElms
          RMaintDefcitRecycOMAutor(NE,M,NGL)   = RMaintDefcitKillOMAutor(NE,M,NGL)*RCCE(NE)
          RMaintDefcitLitrfalOMAutor(NE,M,NGL) = AZMAX1(RMaintDefcitKillOMAutor(NE,M,NGL)-RMaintDefcitRecycOMAutor(NE,M,NGL))
          !
          !     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
          !     PRODUCTS
          !
          !     RHMMC,RHMMN,RHMMC=transfer of senesence LitrFall C,N,P to humus
          !     EHUM=humus transfer fraction
          !     RCMMC,RCMMN,RCMMC=transfer of senesence LitrFall C,N,P to residue
          !
          RMaintDefLitrfal2HumOMAutor(NE,M,NGL)   = RMaintDefcitLitrfalOMAutor(NE,M,NGL)*EHUM
          RMaintDefLitrfal2ResduOMAutor(NE,M,NGL) = RMaintDefcitLitrfalOMAutor(NE,M,NGL)-RMaintDefLitrfal2HumOMAutor(NE,M,NGL)
        ENDDO
      ENDDO
    ELSE
      DO  M=1,2
        DO NE=1,NumPlantChemElms
          RMaintDefcitKillOMAutor(NE,M,NGL)          = 0.0_r8
          RMaintDefcitLitrfalOMAutor(NE,M,NGL)       = 0.0_r8
          RMaintDefcitRecycOMAutor(NE,M,NGL)         = 0.0_r8
          RMaintDefLitrfal2HumOMAutor(NE,M,NGL)   = 0.0_r8
          RMaintDefLitrfal2ResduOMAutor(NE,M,NGL) = 0.0_r8
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherAutotrophAnabolicFlux
!------------------------------------------------------------------------------------------

  subroutine AerobicAutorO2Uptake(I,J,N,micfor,micstt,nmicf,nmics,micflx,naqfdiag)
  implicit none
  integer, intent(in) :: I,J,N     !functional group id
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag  
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FOXYX,OXKX
  integer  :: M,MX,NGL
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,O2AquaDiffusvity1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMX,RO2DmndX
  real(r8) :: ROXYLX
  real(r8) :: RRADO,RMPOX,ROXDFQ
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: dribbling_flx
  real(r8) :: X,VOLOXM
  real(r8) :: VOLWPM

  ! begin_execution
  associate(                                                       &
    fLimO2Autor               => nmics%fLimO2Autor,                &
    OMActAutor                => nmics%OMActAutor,                 &
    FracOMActAutor            => nmics%FracOMActAutor,             &
    RO2UptkAutor              => nmicf%RO2UptkAutor,               &
    RespGrossAutor            => nmicf%RespGrossAutor,             &
    RO2Dmnd4GrossRespAutor    => nmicf%RO2Dmnd4GrossRespAutor,     &
    RO2Uptk4RespAutor         => nmicf%RO2Uptk4RespAutor,          &
    RCO2ProdAutor             => nmicf%RCO2ProdAutor,              &
    RCH4ProdAutor             => nmicf%RCH4ProdAutor,              &
    RSMetaOxidSoilAutor       => nmicf%RSMetaOxidSoilAutor,        &
    RSMetaOxidBandAutor       => nmicf%RSMetaOxidBandAutor,        &
    mid_AutoAmmoniaOxidBacter => micpar%mid_AutoAmmoniaOxidBacter, &
    mid_AutoNitriteOxidBacter => micpar%mid_AutoNitriteOxidBacter, &
    RO2GasXchangePrev         => micfor%RO2GasXchangePrev,         &
    RO2MetaDmndAutorPrev      => micflx%RO2MetaDmndAutorPrev,      &
    RO2MetaDmndAutor          => micflx%RO2MetaDmndAutor,          &
    COXYE                     => micfor%COXYE,                     &
    RO2EcoDmndPrev            => micfor%RO2EcoDmndPrev,            &
    O2_rain_conc              => micfor%O2_rain_conc,              &
    O2_irrig_conc             => micfor%O2_irrig_conc,             &
    Irrig2LitRSurf_col        => micfor%Irrig2LitRSurf_col,        &
    Rain2LitRSurf             => micfor%Rain2LitRSurf,             &
    litrm                     => micfor%litrm,                     &
    O2AquaDiffusvity          => micfor%O2AquaDiffusvity,          &
    RO2AquaXchangePrev        => micfor%RO2AquaXchangePrev,        &
    VLSoilPoreMicP            => micfor%VLSoilPoreMicP,            &
    VLSoilMicP                => micfor%VLSoilMicP,                &
    VLsoiAirPM                => micfor%VLsoiAirPM,                &
    VLWatMicPM                => micfor%VLWatMicPM,                &
    FILM                      => micfor%FILM,                      &
    THETPM                    => micfor%THETPM,                    &
    TortMicPM                 => micfor%TortMicPM,                 &
    ZERO                      => micfor%ZERO,                      &
    ZEROS                     => micfor%ZEROS,                     &
    DiffusivitySolutEff       => micfor%DiffusivitySolutEff,       &
    O2GSolubility             => micstt%O2GSolubility,             &
    JGniA                     => micpar%JGniA,                     &
    JGnfA                     => micpar%JGnfA,                     &
    OXYG                      => micstt%OXYG,                      &
    OXYS                      => micstt%OXYS,                      &
    COXYG                     => micstt%COXYG,                     &
    REcoUptkSoilO2M           => micflx%REcoUptkSoilO2M,           &
    RNH3OxidAutor             => micflx%RNH3OxidAutor,             &
    RNH3OxidAutorBand         => micflx%RNH3OxidAutorBand,         &
    RNO2XupAutor              => micflx%RNO2XupAutor,              &
    RNO2XupAutorBand          => micflx%RNO2XupAutorBand           &
  )

  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle    
    IF(RO2EcoDmndPrev.GT.ZEROS)THEN
      FOXYX=AMAX1(FMN,RO2MetaDmndAutorPrev(NGL)/RO2EcoDmndPrev)
    ELSE
      FOXYX=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF
    naqfdiag%TFOXYX   = naqfdiag%TFOXYX+FOXYX
    OXKX              = OXKA
    RO2UptkAutor(NGL) = 0._r8

    IF(RO2MetaDmndAutor(NGL).GT.ZEROS .AND. FOXYX.GT.ZERO)THEN
      IF(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS)THEN
        !
        !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
        !     POPULATION
        !
        RUPMX             = RO2MetaDmndAutor(NGL)*dts_gas    !
        RO2DmndX          = -RO2GasXchangePrev*dts_gas*FOXYX    !O2 demand
        O2AquaDiffusvity1 = O2AquaDiffusvity*dts_gas
        IF(.not.litrm)THEN
          OXYG1  = OXYG*FOXYX
          ROXYLX = -RO2AquaXchangePrev*dts_gas*FOXYX
        ELSE
          OXYG1  = COXYG*VLsoiAirPM(1)*FOXYX
          ROXYLX = -(RO2AquaXchangePrev+Rain2LitRSurf*O2_rain_conc &
            +Irrig2LitRSurf_col*O2_irrig_conc)*dts_gas*FOXYX
        ENDIF
        if(OXYG1<=0._r8 .and. ROXYLX>0._r8)ROXYLX=0._r8
        OXYS1=OXYS*FOXYX
        !
            !write(*,*)'O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP'
        !     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
        !
        dribbling_flx=0._r8
        DO  M=1,NPH
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
          DIFOX  = TortMicPM(M)*O2AquaDiffusvity1*12.57_r8*BIOS*OMActAutor(NGL)*RRADO
          VOLWOX = VLWatMicPM(M)*O2GSolubility
          VOLPOX = VLsoiAirPM(M)
          VOLWPM = VOLWOX+VOLPOX
          VOLOXM = VLWatMicPM(M)*FOXYX
          !oxygen uptake in the layer
          DO  MX=1,NPT
            call fixEXConsumpFlux(OXYG1,RO2DmndX)
            call fixEXConsumpFlux(OXYS1,ROXYLX)
            COXYS1 = AMIN1(COXYE*O2GSolubility,safe_adb(OXYS1,VOLOXM))

            !solve for uptake flux
            IF(OXYS1<=ZEROS)THEN
              RMPOX=0.0_r8
            else
              RMPOX=TranspBasedsubstrateUptake(COXYS1,DIFOX, OXKX, RUPMX, ZEROS)
            ENDIF

            !apply the uptake flux
            !apply the uptake
            call SubstrateDribbling(RMPOX,dribbling_flx,OXYS1)

            !apply volatilization-dissolution
            IF(THETPM(M).GT.AirFillPore_Min.AND.VOLPOX.GT.ZEROS)THEN
              ROXDFQ=DiffusivitySolutEff(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
              ROXDFQ=AMAX1(AMIN1(OXYG1,ROXDFQ),-OXYS1)
            ELSE
              ROXDFQ=0.0_r8
            ENDIF
            OXYG1 = OXYG1-ROXDFQ
            OXYS1 = OXYS1+ROXDFQ
            !accumulate uptake flux
            RO2UptkAutor(NGL)  = RO2UptkAutor(NGL)+RMPOX
            REcoUptkSoilO2M(M) = REcoUptkSoilO2M(M)+RMPOX
          ENDDO
        ENDDO
        !
        !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (OxyLimterHeter)
        !
        !     OxyLimterHeter=ratio of O2-limited to O2-unlimited uptake
        !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter_vr,RNO2DmndReduxBandHeter_vr=NH3,NO2 oxidation in non-band, band
        !
        fLimO2Autor(NGL)=AMIN1(1.0,AZMAX1(RO2UptkAutor(NGL)/RO2MetaDmndAutor(NGL)))
        IF(N.EQ.mid_AutoAmmoniaOxidBacter)THEN
          RNH3OxidAutor(NGL)     = RNH3OxidAutor(NGL)*fLimO2Autor(NGL)
          RNH3OxidAutorBand(NGL) = RNH3OxidAutorBand(NGL)*fLimO2Autor(NGL)
        ELSEIF(N.EQ.mid_AutoNitriteOxidBacter)THEN
          RNO2XupAutor(NGL)     = RNO2XupAutor(NGL)*fLimO2Autor(NGL)
          RNO2XupAutorBand(NGL) = RNO2XupAutorBand(NGL)*fLimO2Autor(NGL)
        ENDIF
      ELSE
        RO2UptkAutor(NGL) = RO2MetaDmndAutor(NGL)
        fLimO2Autor(NGL)  = 1.0_r8
      ENDIF
    ELSE
      fLimO2Autor(NGL)  = 1.0_r8
    ENDIF
    !
    !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
    !     RCO2X,RAcetateProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
    !     RO2Uptk4RespHeter=O2-limited O2 uptake
    !     RSMetaOxidSoilAutor,RSMetaOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation

    !NH3 oxidizer assimilate CO2, CH4 oxidizer produces CO2, nitrite oxidizer assimilates CO2
    RespGrossAutor(NGL)      = RespGrossAutor(NGL)*fLimO2Autor(NGL)
    RCO2ProdAutor(NGL)       = RespGrossAutor(NGL)
    RCH4ProdAutor(NGL)       = 0.0_r8
    RO2Uptk4RespAutor(NGL)   = RO2Dmnd4GrossRespAutor(NGL)*fLimO2Autor(NGL)
    RSMetaOxidSoilAutor(NGL) = RSMetaOxidSoilAutor(NGL)*fLimO2Autor(NGL)
    RSMetaOxidBandAutor(NGL) = RSMetaOxidBandAutor(NGL)*fLimO2Autor(NGL)
  ENDDO
  end associate
  end subroutine AerobicAutorO2Uptake



!------------------------------------------------------------------------------------------
  subroutine BiomNutMinerMobilAutor(I,J,N,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,&
    nmics,naqfdiag)
  !
  !Description:
  !Do biomass nutrient mobilization and immobilization.   
  implicit none
  integer, intent(in) :: I,J,N
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag  
  real(r8)  :: FNH4X
  real(r8)  :: FNB3X,FNB4X,FNO3X
  real(r8)  :: FPO4X,FPOBX,FP14X,FP1BX  
  real(r8) :: CNH4X,CNH4Y,CNO3X,CNO3Y
  real(r8) :: CH2PX,CH2PY
  real(r8) :: CH1PX,CH1PY
  real(r8) :: FNH4S,FNHBS
  real(r8) :: FNO3S,FNO3B
  real(r8) :: FH1PS,FH1PB
  real(r8) :: FH2PS,FH2PB
  real(r8) :: H2POM,H2PBM
  real(r8) :: H1POM,H1PBM
  real(r8) :: H1P4M,H2P4M
  real(r8) :: RNetNH4MinPotent
  real(r8) :: RINHX,RNetNO3Dmnd,RINOX,RNetH2PO4MinPotent,RIPOX,RNetH1PO4Dmnd
  real(r8) :: RIP1X,RNetNH4MinPotentLitr,RNetNO3DmndLitr,RNetH2PO4MinPotentLitr,RNetH1PO4DmndLitr
  real(r8) :: ZNH4M,ZNHBM
  real(r8) :: ZNO3M
  real(r8) :: ZNOBM
  integer :: MID3,NGL

!     begin_execution
  associate(                                             &
   GrowthEnvScalAutor    => nmics%GrowthEnvScalAutor,    &
   OMActAutor            => nmics%OMActAutor,            &
   AttenfNH4Autor        => micflx%AttenfNH4Autor,       &
   AttenfNO3Autor        => micflx%AttenfNO3Autor,       &
   AttenfH2PO4Autor      => micflx%AttenfH2PO4Autor,     &
   AttenfH1PO4Autor      => micflx%AttenfH1PO4Autor,     &
   RNH4TransfSoilAutor   => nmicf%RNH4TransfSoilAutor,   &
   RNO3TransfSoilAutor   => nmicf%RNO3TransfSoilAutor,   &
   RH2PO4TransfSoilAutor => nmicf%RH2PO4TransfSoilAutor, &
   RNH4TransfBandAutor   => nmicf%RNH4TransfBandAutor,   &
   RNO3TransfBandAutor   => nmicf%RNO3TransfBandAutor,   &
   RH2PO4TransfBandAutor => nmicf%RH2PO4TransfBandAutor, &
   RNH4TransfLitrAutor   => nmicf%RNH4TransfLitrAutor,   &
   RNO3TransfLitrAutor   => nmicf%RNO3TransfLitrAutor,   &
   RH2PO4TransfLitrAutor => nmicf%RH2PO4TransfLitrAutor, &
   RH1PO4TransfSoilAutor => nmicf%RH1PO4TransfSoilAutor, &
   RH1PO4TransfBandAutor => nmicf%RH1PO4TransfBandAutor, &
   RH1PO4TransfLitrAutor => nmicf%RH1PO4TransfLitrAutor, &
   rNCOMCAutor           => micpar%rNCOMCAutor,          &
   rPCOMCAutor           => micpar%rPCOMCAutor,          &
   VLNH4                 => micfor%VLNH4,                &
   VLNHB                 => micfor%VLNHB,                &
   VLWatMicP             => micfor%VLWatMicP,            &
   VLNO3                 => micfor%VLNO3,                &
   VLNOB                 => micfor%VLNOB,                &
   VLPO4                 => micfor%VLPO4,                &
   VLPOB                 => micfor%VLPOB,                &
   litrm                 => micfor%litrm,                &
   mBiomeAutor           => micstt%mBiomeAutor,          &
   ZNH4S                 => micstt%ZNH4S,                &
   ZNH4B                 => micstt%ZNH4B,                &
   ZNO3S                 => micstt%ZNO3S,                &
   ZNO3B                 => micstt%ZNO3B,                &
   CNO3S                 => micstt%CNO3S,                &
   CNO3B                 => micstt%CNO3B,                &
   CNH4S                 => micstt%CNH4S,                &
   CNH4B                 => micstt%CNH4B,                &
   CH2P4                 => micstt%CH2P4,                &
   CH2P4B                => micstt%CH2P4B,               &
   H2PO4                 => micstt%H2PO4,                &
   H2POB                 => micstt%H2POB,                &
   CH1P4                 => micstt%CH1P4,                &
   CH1P4B                => micstt%CH1P4B,               &
   H1PO4                 => micstt%H1PO4,                &
   H1POB                 => micstt%H1POB,                &
   JGniA                 => micpar%JGniA,                &
   JGnfA                 => micpar%JGnfA,                &
   RNH4UptkSoilAutor     => micflx%RNH4UptkSoilAutor,    &
   RNH4UptkBandAutor     => micflx%RNH4UptkBandAutor,    &
   RNO3UptkSoilAutor     => micflx%RNO3UptkSoilAutor,    &
   RNO3UptkBandAutor     => micflx%RNO3UptkBandAutor,    &
   NetNH4Mineralize      => micflx%NetNH4Mineralize,     &
   RH2PO4UptkSoilAutor   => micflx%RH2PO4UptkSoilAutor,  &
   RH2PO4UptkBandAutor   => micflx%RH2PO4UptkBandAutor,  &
   NetPO4Mineralize      => micflx%NetPO4Mineralize,     &
   RH1PO4UptkSoilAutor   => micflx%RH1PO4UptkSoilAutor,  &
   RH1PO4UptkBandAutor   => micflx%RH1PO4UptkBandAutor,  &
   RH1PO4UptkLitrAutor   => micflx%RH1PO4UptkLitrAutor,  &
   RNH4UptkLitrAutor     => micflx%RNH4UptkLitrAutor,    &
   RNO3UptkLitrAutor     => micflx%RNO3UptkLitrAutor,    &
   RH2PO4UptkLitrAutor   => micflx%RH2PO4UptkLitrAutor   &
  )
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
!     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetNH4MinPotent=NH4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMN=microbial nonstructural C,N
!     rNCOMC=maximum microbial N:C ratio
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     RINHX=microbially limited NH4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     VLWatMicP=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RNH4imobilSoilHeter,RNH4imobilBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle      
    call SubstrateCompetAuto(NGL,N,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
        micfor,naqfdiag,nmicf,nmics,micflx)

    FNH4S=VLNH4
    FNHBS=VLNHB
    MID3=micpar%get_micb_id(ibiom_reserve,NGL)
    RNetNH4MinPotent=mBiomeAutor(ielmc,MID3)*rNCOMCAutor(ibiom_reserve,NGL)-mBiomeAutor(ielmn,MID3)
    IF(RNetNH4MinPotent.GT.0.0_r8)THEN
      CNH4X                    = AZMAX1(CNH4S-Z4MN)
      CNH4Y                    = AZMAX1(CNH4B-Z4MN)
      RINHX                    = AMIN1(RNetNH4MinPotent,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*Z4MX)
      RNH4UptkSoilAutor(NGL)   = FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
      RNH4UptkBandAutor(NGL)   = FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
      ZNH4M                    = Z4MN*VLWatMicP*FNH4S
      ZNHBM                    = Z4MN*VLWatMicP*FNHBS
      RNH4TransfSoilAutor(NGL) = AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4UptkSoilAutor(NGL))
      RNH4TransfBandAutor(NGL) = AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4UptkBandAutor(NGL))
    ELSE
      RNH4UptkSoilAutor(NGL)   = 0.0_r8
      RNH4UptkBandAutor(NGL)   = 0.0_r8
      RNH4TransfSoilAutor(NGL) = RNetNH4MinPotent*FNH4S
      RNH4TransfBandAutor(NGL) = RNetNH4MinPotent*FNHBS
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+(RNH4TransfSoilAutor(NGL)+RNH4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
!     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetNO3Dmnd=NO3 immobilization (+ve) demand
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     min NO3 concentration and Km for NO3 uptake
!     RINOX=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     VLWatMicP=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RNO3imobilSoilHeter,RNO3imobilBandHeter=substrate-limited NO3 immobiln in non-band, band
!     NetNH4Mineralize=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
    FNO3S=VLNO3
    FNO3B=VLNOB
    RNetNO3Dmnd=AZMAX1(RNetNH4MinPotent-RNH4TransfSoilAutor(NGL)-RNH4TransfBandAutor(NGL))
    IF(RNetNO3Dmnd.GT.0.0_r8)THEN
      CNO3X                    = AZMAX1(CNO3S-ZOMN)
      CNO3Y                    = AZMAX1(CNO3B-ZOMN)
      RINOX                    = AMIN1(RNetNO3Dmnd,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*ZOMX)
      RNO3UptkSoilAutor(NGL)   = FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
      RNO3UptkBandAutor(NGL)   = FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
      ZNO3M                    = ZOMN*VLWatMicP*FNO3S
      ZNOBM                    = ZOMN*VLWatMicP*FNO3B
      RNO3TransfSoilAutor(NGL) = AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)),RNO3UptkSoilAutor(NGL))
      RNO3TransfBandAutor(NGL) = AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)),RNO3UptkBandAutor(NGL))
    ELSE
      RNO3UptkSoilAutor(NGL)   = 0.0_r8
      RNO3UptkBandAutor(NGL)   = 0.0_r8
      RNO3TransfSoilAutor(NGL) = RNetNO3Dmnd*FNO3S
      RNO3TransfBandAutor(NGL) = RNetNO3Dmnd*FNO3B
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+(RNO3TransfSoilAutor(NGL)+RNO3TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetH2PO4MinPotent=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMP=microbial nonstructural C,P
!     rPCOMC=maximum microbial P:C ratio
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     min H2PO4 concentration and Km for H2PO4 uptake
!     RIPOX=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RH2PO4DmndSoilHeter_vr,RH2PO4DmndBandHeter_vr=substrate-unlimited H2PO4 mineraln-immobiln
!     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
!     VOLW=water content
!     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
!     RH2PO4imobilSoilHeter,RH2PO4imobilBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    FH2PS=VLPO4
    FH2PB=VLPOB
    MID3=micpar%get_micb_id(ibiom_reserve,NGL)
    RNetH2PO4MinPotent=(mBiomeAutor(ielmc,MID3)*rPCOMCAutor(ibiom_reserve,NGL)-mBiomeAutor(ielmp,MID3))
    IF(RNetH2PO4MinPotent.GT.0.0)THEN
      CH2PX=AZMAX1(CH2P4-HPMN)
      CH2PY=AZMAX1(CH2P4B-HPMN)
      RIPOX=AMIN1(RNetH2PO4MinPotent,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX)
      RH2PO4UptkSoilAutor(NGL)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
      RH2PO4UptkBandAutor(NGL)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
      H2POM=HPMN*VLWatMicP*FH2PS
      H2PBM=HPMN*VLWatMicP*FH2PB
      RH2PO4TransfSoilAutor(NGL)=AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RH2PO4UptkSoilAutor(NGL))
      RH2PO4TransfBandAutor(NGL)=AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RH2PO4UptkBandAutor(NGL))
    ELSE
      RH2PO4UptkSoilAutor(NGL)=0.0_r8
      RH2PO4UptkBandAutor(NGL)=0.0_r8
      RH2PO4TransfSoilAutor(NGL)=RNetH2PO4MinPotent*FH2PS
      RH2PO4TransfBandAutor(NGL)=RNetH2PO4MinPotent*FH2PB
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+(RH2PO4TransfSoilAutor(NGL)+RH2PO4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetH1PO4Dmnd=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     min HPO4 concentration and Km for HPO4 uptake
!     RIP1X=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RH1PO4DmndSoilHeter_vr,RH1PO4DmndBandHeter_vr=substrate-unlimited HPO4 mineraln-immobiln
!     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
!     VOLW=water content
!     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
!     RH1PO4imobilSoilHeter,RH1PO4imobilBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RNetH1PO4Dmnd=0.1_r8*AZMAX1(RNetH2PO4MinPotent-RH2PO4TransfSoilAutor(NGL)-RH2PO4TransfBandAutor(NGL))
    IF(RNetH1PO4Dmnd.GT.0.0_r8)THEN
      CH1PX=AZMAX1(CH1P4-HPMN)
      CH1PY=AZMAX1(CH1P4B-HPMN)
      RIP1X=AMIN1(RNetH1PO4Dmnd,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX)
      RH1PO4UptkSoilAutor(NGL)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
      RH1PO4UptkBandAutor(NGL)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
      H1POM=HPMN*VLWatMicP*FH1PS
      H1PBM=HPMN*VLWatMicP*FH1PB
      RH1PO4TransfSoilAutor(NGL)=AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RH1PO4UptkSoilAutor(NGL))
      RH1PO4TransfBandAutor(NGL)=AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RH1PO4UptkBandAutor(NGL))
    ELSE
      RH1PO4UptkSoilAutor(NGL)=0.0_r8
      RH1PO4UptkBandAutor(NGL)=0.0_r8
      RH1PO4TransfSoilAutor(NGL)=RNetH1PO4Dmnd*FH1PS
      RH1PO4TransfBandAutor(NGL)=RNetH1PO4Dmnd*FH1PB
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+(RH1PO4TransfSoilAutor(NGL)+RH1PO4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RNetNH4MinPotentLitr=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RNH4DmndLitrHeter=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M=NH4 not available for uptake
!     AttenfNH4Heter=fractions of biological NH4 demand
!     RNH4imobilLitrHeter=substrate-limited NH4 mineraln-immobiln
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
!
    IF(litrm)THEN
      RNetNH4MinPotentLitr=RNetNH4MinPotent-RNH4TransfSoilAutor(NGL)-RNO3TransfSoilAutor(NGL)
      IF(RNetNH4MinPotentLitr.GT.0.0_r8)THEN
        CNH4X=AZMAX1(CNH4S-Z4MN)
        CNH4Y=AZMAX1(CNH4B-Z4MN)
        RNH4UptkLitrAutor(NGL)=AMIN1(RNetNH4MinPotentLitr,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*Z4MX) &
            *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
        ZNH4M=Z4MN*VLWatMicP
        RNH4TransfLitrAutor(NGL)=AMIN1(AttenfNH4Autor(NGL)*AZMAX1((ZNH4T-ZNH4M)),RNH4UptkLitrAutor(NGL))
      ELSE
        RNH4UptkLitrAutor(NGL)=0.0_r8
        RNH4TransfLitrAutor(NGL)=RNetNH4MinPotentLitr
      ENDIF
      NetNH4Mineralize=NetNH4Mineralize+RNH4TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RNetNO3DmndLitr=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     minimum NO3 concentration and Km for NO3 uptake
!     RNO3DmndLitrHeter_col=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RNO3imobilLitrHeter=substrate-unlimited NO3 immobiln
!     VLWatMicP=water content
!     ZNO3M=NO3 not available for uptake
!     AttenfNO3Heter=fraction of biological NO3 demand
!     RNO3imobilLitrHeter=substrate-limited NO3 immobiln
!     NetNH4Mineralize=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
      RNetNO3DmndLitr=AZMAX1(RNetNH4MinPotentLitr-RNH4TransfLitrAutor(NGL))
      IF(RNetNO3DmndLitr.GT.0.0_r8)THEN
        CNO3X=AZMAX1(CNO3S-ZOMN)
        CNO3Y=AZMAX1(CNO3B-ZOMN)
        RNO3UptkLitrAutor(NGL)=AMAX1(RNetNO3DmndLitr,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*ZOMX) &
            *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
        ZNO3M=ZOMN*VLWatMicP
        RNO3TransfLitrAutor(NGL)=AMIN1(AttenfNO3Autor(NGL)*AZMAX1((ZNO3T-ZNO3M)),RNO3UptkLitrAutor(NGL))
      ELSE
        RNO3UptkLitrAutor(NGL)=0.0_r8
        RNO3TransfLitrAutor(NGL)=RNetNO3DmndLitr
      ENDIF
      NetNH4Mineralize=NetNH4Mineralize+RNO3TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RNetH2PO4MinPotentLitr=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     minimum H2PO4 concentration and Km for H2PO4 uptake
!     RH2PO4DmndLitrHeter=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RH2PO4DmndLitrHeter=substrate-unlimited H2PO4 mineraln-immobiln
!     VLWatMicP=water content
!     H2P4M=H2PO4 not available for uptake
!     AttenfH2PO4Heter=fractions of biological H2PO4 demand
!     RH2PO4imobilLitrHeter=substrate-limited H2PO4 mineraln-immobiln
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
      RNetH2PO4MinPotentLitr=RNetH2PO4MinPotent-RH2PO4TransfSoilAutor(NGL)
      IF(RNetH2PO4MinPotentLitr.GT.0.0_r8)THEN
        CH2PX=AZMAX1(CH2P4-HPMN)
        CH2PY=AZMAX1(CH2P4B-HPMN)
        RH2PO4UptkLitrAutor(NGL)=AMIN1(RNetH2PO4MinPotentLitr,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX) &
            *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
        H2P4M=HPMN*VLWatMicP
        RH2PO4TransfLitrAutor(NGL)=AMIN1(AttenfH2PO4Autor(NGL)*AZMAX1((H2P4T-H2P4M)),RH2PO4UptkLitrAutor(NGL))
      ELSE
        RH2PO4UptkLitrAutor(NGL)=0.0_r8
        RH2PO4TransfLitrAutor(NGL)=RNetH2PO4MinPotentLitr
      ENDIF
      NetPO4Mineralize=NetPO4Mineralize+RH2PO4TransfLitrAutor(NGL)
      !
      !     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
      !     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
      !     ZONES OF SOIL SURFACE
      !
      !     RNetH1PO4DmndLitr=HPO4 mineralization (-ve) or immobilization (+ve) demand
      !     NU=surface layer number
      !     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
      !     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
      !     minimum HPO4 concentration and Km for HPO4 uptake
      !     RH1PO4DmndLitrHeter_col=microbially limited HPO4 demand
      !     BIOA=microbial surface area, OMA=active biomass
      !     TFNG=temp+water stress
      !     FH1PS,FH1PB=fractions of HPO4 in non-band, band
      !     RH1PO4DmndLitrHeter_col=substrate-unlimited HPO4 mineraln-immobiln
      !     VLWatMicP=water content
      !     H1P4M=HPO4 not available for uptake
      !     AttenfH1PO4Heter=fraction of biological HPO4 demand
      !     RH1PO4imobilLitrHeter=substrate-limited HPO4 minereraln-immobiln
      !     NetPO4Mineralize=total HPO4 net mineraln (-ve) or immobiln (+ve)
      !
      FH1PS = VLPO4
      FH1PB = VLPOB
      RNetH1PO4DmndLitr=0.1_r8*AZMAX1(RNetH2PO4MinPotentLitr-RH2PO4TransfLitrAutor(NGL))
      IF(RNetH1PO4DmndLitr.GT.0.0_r8)THEN
        CH1PX=AZMAX1(CH1P4-HPMN)
        CH1PY=AZMAX1(CH1P4B-HPMN)
        RH1PO4UptkLitrAutor(NGL)=AMIN1(RNetH1PO4DmndLitr,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX) &
            *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
        H1P4M=HPMN*VLWatMicP
        RH1PO4TransfLitrAutor(NGL)=AMIN1(AttenfH1PO4Autor(NGL)*AZMAX1((H1P4T-H1P4M)),RH1PO4UptkLitrAutor(NGL))
      ELSE
        RH1PO4UptkLitrAutor(NGL)=0.0_r8
        RH1PO4TransfLitrAutor(NGL)=RNetH1PO4DmndLitr
      ENDIF
      NetPO4Mineralize=NetPO4Mineralize+RH1PO4TransfLitrAutor(NGL)
    ENDIF
  ENDDO
  end associate
  end subroutine BiomNutMinerMobilAutor
!------------------------------------------------------------------------------------------

  subroutine GatherAutotrophRespiration(I,J,N,micfor,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: I,J,N
  type(micforctype), intent(in) :: micfor

  type(micfluxtype), intent(inout) :: micflx  
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  character(len=*), parameter :: subname='GatherAutotrophRespiration'
  real(r8) :: RGN2P
  integer  :: NGL
!     begin_execution
  associate(                                             &
    OMActAutor           => nmics%OMActAutor,            &  
    RespGrossAutor       => nmicf%RespGrossAutor,        &
    Resp4NFixAutor       => nmicf%Resp4NFixAutor,        &
    RN2FixAutor          => nmicf%RN2FixAutor,           &
    RGrowthRespAutor     => micflx%RGrowthRespAutor,     &
    RMaintDefcitcitAutor => micflx%RMaintDefcitcitAutor, &
    JGniA                => micpar%JGniA,                &
    JGnfA                => micpar%JGnfA,                &
    RMaintRespAutor      => micflx%RMaintRespAutor       &
  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TempMaintRHeter=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!
  call PrintInfo('beg '//subname)
  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle      
    RGrowthRespAutor(NGL)     = AZMAX1(RespGrossAutor(NGL)-RMaintRespAutor(NGL))
    RMaintDefcitcitAutor(NGL) = AZMAX1(RMaintRespAutor(NGL)-RespGrossAutor(NGL))
    !
    !     N2 FIXATION: N=(6) AEROBIC, (7) ANAEROBIC
    !     FROM GROWTH RESPIRATION, FIXATION ENERGY REQUIREMENT,
    !     MICROBIAL N REQUIREMENT IN LABILE (1) AND
    !     RESISTANT (2) FRACTIONS
    !
    !     RGN2P=respiration to meet N2 fixation demand
    !     OMC,OMN=microbial nonstructural C,N
    !     rNCOMC=maximum microbial N:C ratio
    !     EN2F=N2 fixation yield per unit nonstructural C
    !     RGrowthRespAutor=growth respiration
    !     Resp4NFixHeter=respiration for N2 fixation
    !     CZ2GS=aqueous N2 concentration
    !     ZFKM=Km for N2 uptake
    !     OMGR*OMC(3,NGL,N,K)=nonstructural C limitation to Resp4NFixHeter
    !     RN2FixHeter=N2 fixation rate
    !
    RN2FixAutor(NGL)    = 0.0_r8
    Resp4NFixAutor(NGL) = 0.0_r8
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherAutotrophRespiration
!------------------------------------------------------------------------------------------

  subroutine AutotrophAnabolicUpdate(micfor,micstt,nmicf,nmicdiag)

  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag    
  character(len=*), parameter :: subname='AutotrophAnabolicUpdate'

  real(r8) :: CGROMC   !C for microbial biomass growth
  integer :: N,M,NGL,MID,MID3,NE
  associate(                                                                &
    DOMuptk4GrothAutor             => nmicf%DOMuptk4GrothAutor,             &
    NonstX2stBiomAutor             => nmicf%NonstX2stBiomAutor,             &
    Resp4NFixAutor                 => nmicf%Resp4NFixAutor,                 &
    RespGrossAutor                 => nmicf%RespGrossAutor,                 &
    RNOxReduxRespAutorLim          => nmicf%RNOxReduxRespAutorLim,          &
    RNO3TransfSoilAutor            => nmicf%RNO3TransfSoilAutor,            &
    RCO2ProdAutor                  => nmicf%RCO2ProdAutor,                  &
    ROMProdCO2Autor                => nmicf%ROMProdCO2Autor,                &
    RGrowthCAutor                  => nmicf%RGrowthCAutor,                  & !growth C biomass
    RCO2XumpAutor                  => nmicf%RCO2XumpAutor,                  & !CO2 uptake for biomass and respiration
    RH2PO4TransfSoilAutor          => nmicf%RH2PO4TransfSoilAutor,          &
    RNH4TransfBandAutor            => nmicf%RNH4TransfBandAutor,            &
    RNO3TransfBandAutor            => nmicf%RNO3TransfBandAutor,            &
    RH2PO4TransfBandAutor          => nmicf%RH2PO4TransfBandAutor,          &
    RkillLitrfal2HumOMAutor        => nmicf%RkillLitrfal2HumOMAutor,        &
    RMaintDefLitrfal2HumOMAutor    => nmicf%RMaintDefLitrfal2HumOMAutor,    &
    RN2FixAutor                    => nmicf%RN2FixAutor,                    &
    RKillOMAutor                   => nmicf%RKillOMAutor,                   &
    RkillRecycOMAutor              => nmicf%RkillRecycOMAutor,              &
    RMaintDefcitKillOMAutor        => nmicf%RMaintDefcitKillOMAutor,        &
    RMaintDefcitRecycOMAutor       => nmicf%RMaintDefcitRecycOMAutor,       &
    RNH4TransfLitrAutor            => nmicf%RNH4TransfLitrAutor,            &
    RNO3TransfLitrAutor            => nmicf%RNO3TransfLitrAutor,            &
    RH2PO4TransfLitrAutor          => nmicf%RH2PO4TransfLitrAutor,          &
    RH1PO4TransfSoilAutor          => nmicf%RH1PO4TransfSoilAutor,          &
    RH1PO4TransfBandAutor          => nmicf%RH1PO4TransfBandAutor,          &
    RH1PO4TransfLitrAutor          => nmicf%RH1PO4TransfLitrAutor,          &
    RNH4TransfSoilAutor            => nmicf%RNH4TransfSoilAutor,            &
    litrm                          => micfor%litrm,                         &
    ElmAllocmatMicrblitr2POM       => micfor%ElmAllocmatMicrblitr2POM,      &
    ElmAllocmatMicrblitr2POMU      => micfor%ElmAllocmatMicrblitr2POMU,     &
    SolidOM                        => micstt%SolidOM,                       &
    mBiomeAutor                    => micstt%mBiomeAutor,                   &
    SOMHumProtein                  => micstt%SOMHumProtein,                 &
    SOMHumCarbohyd                 => micstt%SOMHumCarbohyd,                &
    mid_AutoAmmoniaOxidBacter      => micpar%mid_AutoAmmoniaOxidBacter,     &
    mid_AutoAMONC10                => micpar%mid_AutoAMONC10          ,     &
    mid_AutoAMOANME2D              => micpar%mid_AutoAMOANME2D            , &
    mid_AutoNitriteOxidBacter      => micpar%mid_AutoNitriteOxidBacter,     &
    mid_AutoH2GenoCH4GenArchea     => micpar%mid_AutoH2GenoCH4GenArchea ,   &
    JGniA                          => micpar%JGniA,                         &
    JGnfA                          => micpar%JGnfA,                         &
    NumMicbAFunGrupsPerCmplx       => micpar%NumMicbAFunGrupsPerCmplx,      &
    icarbhyro                      => micpar%icarbhyro,                     &
    iprotein                       => micpar%iprotein,                      &
    k_POM                          => micpar%k_POM,                         &
    is_activeMicrbFungrpAutor      => micpar%is_activeMicrbFungrpAutor      &
  )
  call PrintInfo('beg '//subname)
  DO  N=1,NumMicbAFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        DO  M=1,2
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor(NE,MID)=mBiomeAutor(NE,MID)+NonstX2stBiomAutor(NE,M,NGL)-RKillOMAutor(NE,M,NGL)-RMaintDefcitKillOMAutor(NE,M,NGL)
          ENDDO

!     HUMIFICATION PRODUCTS
!
!     ElmAllocmatMicrblitr2POM=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence LitrFall C,N,P to humus
!
          IF(.not.litrm)THEN
            DO NE=1,NumPlantChemElms
              SolidOM(NE,iprotein,k_POM)=SolidOM(NE,iprotein,k_POM)+ElmAllocmatMicrblitr2POM(1) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
              SolidOM(NE,icarbhyro,k_POM)=SolidOM(NE,icarbhyro,k_POM)+ElmAllocmatMicrblitr2POM(2)&
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
            ENDDO
          ELSE
            DO NE=1,NumPlantChemElms
              SOMHumProtein(NE)=SOMHumProtein(NE)+ElmAllocmatMicrblitr2POMU(1) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
              SOMHumCarbohyd(NE)=SOMHumCarbohyd(NE)+ElmAllocmatMicrblitr2POMU(2) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
            ENDDO
          ENDIF
        ENDDO
        !
        !     INPUTS TO NONSTRUCTURAL POOLS
        !
        !     CGOMC=total DOC+acetate uptake
        !     RespGrossHeter=total respiration
        !     RNOxReduxRespDenitLim=respiration for denitrifcation
        !     Resp4NFixHeter=respiration for N2 fixation
        !     RCO2X=total CO2 emission
        !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
        !     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
        !     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
        !     CGOMN,CGOMP=DON, DOP uptake
        !
        CGROMC             = DOMuptk4GrothAutor(ielmc,NGL)-RespGrossAutor(NGL)-RNOxReduxRespAutorLim(NGL)-Resp4NFixAutor(NGL)
        RGrowthCAutor(NGL) = CGROMC
        RCO2ProdAutor(NGL) = RCO2ProdAutor(NGL)+Resp4NFixAutor(NGL)
        MID3               = micpar%get_micb_id(ibiom_reserve,NGL)

        if(N.eq.mid_AutoAMONC10)then
          !environmental CO2 is assimilated for C biomass
          RCO2XumpAutor(NGL)= CGROMC
        elseif(N.EQ.mid_AutoAMOANME2D)then
          !recyle some CO2
          RCO2ProdAutor(NGL)=RCO2ProdAutor(NGL)-CGROMC
        elseif(N.eq.mid_AutoH2GenoCH4GenArchea .or. N.eq.mid_AutoNitriteOxidBacter .or. N.eq.mid_AutoAmmoniaOxidBacter)then
          !for H2-methanogen, CO2 is used for CH4 and biomass
          ! for NH3/NO2 oxidizer, some CO2 is first converted into CH2O and respired, and some is fixed right away
          RCO2XumpAutor(NGL)= DOMuptk4GrothAutor(ielmc,NGL)          

          if(N.eq.mid_AutoH2GenoCH4GenArchea)then
            !CO2 + 2H2 -> CH2O + 2H2O, 4/12=0.333
            nmicdiag%RH2UptkAutor=nmicdiag%RH2UptkAutor+0.333_r8*CGROMC
          endif  
        endif

        DO M=1,2
          DO NE=1,NumPlantChemElms
            mBiomeAutor(NE,MID3)=mBiomeAutor(NE,MID3)-NonstX2stBiomAutor(NE,M,NGL)+RkillRecycOMAutor(NE,M,NGL)
          ENDDO

          !C is respired as CO2 while N and P are recycled.
          mBiomeAutor(ielmn,MID3) = mBiomeAutor(ielmn,MID3)+RMaintDefcitRecycOMAutor(ielmn,M,NGL)
          mBiomeAutor(ielmp,MID3) = mBiomeAutor(ielmp,MID3)+RMaintDefcitRecycOMAutor(ielmp,M,NGL)
          RCO2ProdAutor(NGL)      = RCO2ProdAutor(NGL)+RMaintDefcitRecycOMAutor(ielmc,M,NGL)
          ROMProdCO2Autor(NGL)    =  ROMProdCO2Autor(NGL)+RMaintDefcitRecycOMAutor(ielmc,M,NGL)
        ENDDO
        
        mBiomeAutor(ielmc,MID3)=mBiomeAutor(ielmc,MID3)+CGROMC
        mBiomeAutor(ielmn,MID3)=mBiomeAutor(ielmn,MID3)+DOMuptk4GrothAutor(ielmn,NGL) &
          +RNH4TransfSoilAutor(NGL)+RNH4TransfBandAutor(NGL)+RNO3TransfSoilAutor(NGL) &
          +RNO3TransfBandAutor(NGL)+RN2FixAutor(NGL)
        
        mBiomeAutor(ielmp,MID3)=mBiomeAutor(ielmp,MID3)+DOMuptk4GrothAutor(ielmp,NGL) &
          +RH2PO4TransfSoilAutor(NGL)+RH2PO4TransfBandAutor(NGL)+RH1PO4TransfSoilAutor(NGL)&
          +RH1PO4TransfBandAutor(NGL)
        IF(litrm)THEN
          mBiomeAutor(ielmn,MID3)=mBiomeAutor(ielmn,MID3)+RNH4TransfLitrAutor(NGL)+RNO3TransfLitrAutor(NGL)
          mBiomeAutor(ielmp,MID3)=mBiomeAutor(ielmp,MID3)+RH2PO4TransfLitrAutor(NGL)+RH1PO4TransfLitrAutor(NGL)
        ENDIF
      enddo
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AutotrophAnabolicUpdate

end module MicAutoCPLXMod
