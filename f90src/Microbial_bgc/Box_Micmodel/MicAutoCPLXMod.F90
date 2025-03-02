module MicAutoCPLXMod
! USES:
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use minimathmod,          only: safe_adb, AZMAX1, fixEXConsumpFlux
  use MicForcTypeMod,       only: micforctype
  use MicFluxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use EcoSiMParDataMod,     only: micpar
  use MicrobeDiagTypes
  use ElmIDMod
  use TracerIDMod
  use EcoSIMSolverPar
  use EcoSimConst
  use NitroPars
  use MicrobMathFuncMod
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ActiveAutotrophs
  public :: AutotrophAnabolicUpdate
  contains

!------------------------------------------------------------------------------------------
  subroutine ActiveAutotrophs(I,J,NGL,N,VOLWZ,XCO2,TSensGrowth,WatStressMicb,&
    SPOMK, RMOMK, OXKX,TotActMicrobiom,TotBiomNO2Consumers,RH2UptkAutor,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: VOLWZ
  real(r8), intent(in):: OXKX,WatStressMicb,TSensGrowth,XCO2
  real(r8), intent(in) :: TotActMicrobiom,TotBiomNO2Consumers
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  real(r8), intent(in) :: SPOMK(2)
  real(r8), intent(in)  :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  real(r8), intent(out):: RH2UptkAutor
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  integer  :: M
  real(r8) :: COMC
  real(r8) :: ECHZ
  real(r8) :: FOXYX
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX
  real(r8) :: FOQC,FOQA
  real(r8) :: FGOCP,FGOAP
  real(r8) :: RGOCP
  real(r8) :: RGOMP
  real(r8) :: RVOXP
  real(r8) :: RVOXPA    !oxidation in soil, for NH3, NO2 and CH4
  real(r8) :: RVOXPB    !oxidation in band, specifically for NH3 and NO2
  real(r8) :: RGrowthRespAutor
  real(r8) :: RMaintDefcitcitAutor
  real(r8) :: RMaintRespAutor

! begin_execution
  associate(                                                             &
    FracOMActAutor               => nmics%FracOMActAutor,                &
    FracNO2ReduxAutor            => nmics%FracNO2ReduxAutor,             &
    FracAutorBiomOfActK          => nmics%FracAutorBiomOfActK,           &
    OMActAutor                   => nmics%OMActAutor,                    &
    RO2UptkAutor                 => nmicf%RO2UptkAutor,                  &
    RespGrossAutor               => nmicf%RespGrossAutor,                &
    RO2Dmnd4RespAutor            => nmicf%RO2Dmnd4RespAutor,             &
    RO2DmndAutor                 => nmicf%RO2DmndAutor,                  &
    RO2Uptk4RespAutor            => nmicf%RO2Uptk4RespAutor,             &
    RNO3UptkAutor                => nmicf%RNO3UptkAutor,                 &
    RNO2ReduxAutorSoil           => nmicf%RNO2ReduxAutorSoil,            &
    RNO2ReduxAutorBand           => nmicf%RNO2ReduxAutorBand,            &
    RNOxReduxRespAutorLim        => nmicf%RNOxReduxRespAutorLim,         &
    RCO2ProdAutor                => nmicf%RCO2ProdAutor,                 &
    RCH4ProdAutor                => nmicf%RCH4ProdAutor,                 &
    TOMEAutoK                    => ncplxs%TOMEAutoK,                    &
    mid_AmmoniaOxidBacter        => micpar%mid_AmmoniaOxidBacter,        &
    mid_NitriteOxidBacter        => micpar%mid_NitriteOxidBacter,        &
    mid_H2GenoMethanogArchea     => micpar%mid_H2GenoMethanogArchea,     &
    mid_AerobicMethanotrofBacter => micpar%mid_AerobicMethanotrofBacter, &
    ZEROS                        => micfor%ZEROS,                        &
    SoilMicPMassLayer            => micfor%SoilMicPMassLayer,            &
    litrm                        => micfor%litrm,                        &
    VLSoilPoreMicP               => micfor%VLSoilPoreMicP,               &
    RO2DmndAutort                => micflx%RO2DmndAutort                 &
  )
! FracOMActHeter,FOMN=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActAutor(NGL)=OMActAutor(NGL)/TotActMicrobiom
  ELSE
    FracOMActAutor(NGL)=1.0_r8
  ENDIF

  IF(TotBiomNO2Consumers.GT.ZEROS.and.N.eq.mid_AmmoniaOxidBacter)THEN
    FracNO2ReduxAutor(NGL)=OMActAutor(NGL)/TotBiomNO2Consumers
  ELSE
    FracNO2ReduxAutor(NGL)=1.0_r8
  ENDIF

  IF(TOMEAutoK(ielmc).GT.ZEROS)THEN
    FracAutorBiomOfActK(NGL)=OMActAutor(NGL)/TOMEAutoK(ielmc)
  ELSE
    FracAutorBiomOfActK(NGL)=1.0_r8
  ENDIF
!
! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
            !write(*,*)'SubstrateCompetitionFactors'
  call SubstrateCompetitionFactorsff(NGL,N,FOXYX,FNH4X,&
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    micfor,naqfdiag,nmicf,nmics,micflx)
!
! HETEROTROPHIC BIOMASS RESPIRATION
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
  if (N.eq.mid_AmmoniaOxidBacter)then
    ! NH3 OXIDIZERS
    call NH3OxidizerCatabolism(NGL,N,XCO2,VOLWZ,TSensGrowth,ECHZ,RGOMP,RVOXP,&
      RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.mid_NitriteOxidBacter)then
    ! NO2 OXIDIZERS
    call NO2OxidizerCatabolism(NGL,N,XCO2,ECHZ,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.mid_H2GenoMethanogArchea)then
    ! H2TROPHIC METHANOGENS
    call H2MethanogensCatabolism(NGL,N,ECHZ,RGOMP,XCO2,micfor,micstt,&
      naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.mid_AerobicMethanotrofBacter)then
    ! METHANOTROPHS
    call MethanotrophCatabolism(I,J,NGL,N,ECHZ,RGOMP,&
      RVOXP,RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  else
    RGOMP                  = 0.0_r8
    RO2Dmnd4RespAutor(NGL) = 0.0_r8
    RO2DmndAutor(NGL)      = 0.0_r8
    RO2DmndAutort(NGL)     = 0.0_r8
  ENDif
!
!  write(*,*)'O2 UPTAKE BY AEROBES'
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
  RO2UptkAutor(NGL)=0.0_r8

  if (N.eq.mid_AmmoniaOxidBacter .or. N.eq.mid_NitriteOxidBacter .or. N.eq.mid_AerobicMethanotrofBacter)then
    call AerobicAutorO2Uptake(NGL,N,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,micfor,micstt,nmicf,nmics,micflx)
  elseif (N.eq.mid_H2GenoMethanogArchea)then
    RespGrossAutor(NGL)    = RGOMP
    RCO2ProdAutor(NGL)     = 0.0_r8
    RCH4ProdAutor(NGL)     = RespGrossAutor(NGL)
    RO2Uptk4RespAutor(NGL) = RO2Dmnd4RespAutor(NGL)
    RH2UptkAutor           = 0.667_r8*RespGrossAutor(NGL)
  ENDIF
!
!     AUTOTROPHIC DENITRIFICATION
!
  IF(N.EQ.mid_AmmoniaOxidBacter .AND. RO2Dmnd4RespAutor(NGL).GT.0.0_r8 .AND. (.not.litrm.OR.VLSoilPoreMicP.GT.ZEROS))THEN
    call AutotrophDenitrificCatabolism(NGL,N,XCO2,VOLWZ,micfor,micstt,&
      naqfdiag,nmicf,nmics,micflx)
  ELSE
    RNO3UptkAutor(NGL)         = 0.0_r8
    RNO2ReduxAutorSoil(NGL)    = 0.0_r8
    RNO2ReduxAutorBand(NGL)    = 0.0_r8
    RNOxReduxRespAutorLim(NGL) = 0.0_r8
  ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
  call BiomassMineralizationff(NGL,N,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,nmics)
!
  call GatherAutotrophRespiration(NGL,N,RMOMK,micfor,micstt,RGrowthRespAutor,&
    RMaintDefcitcitAutor,RMaintRespAutor,nmicf,nmics)
!
  call GatherAutotrophAnabolicFlux(I,J,NGL,N,ECHZ,FGOCP,FGOAP,RGrowthRespAutor,&
    RMaintDefcitcitAutor,RMaintRespAutor,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)

  end associate
  end subroutine ActiveAutotrophs


!------------------------------------------------------------------------------------------

  subroutine SubstrateCompetitionFactorsff(NGL,N,FOXYX,&
    FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    micfor,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(out):: FOXYX,FNH4X
  real(r8),intent(out) :: FNB3X,FNB4X,FNO3X
  real(r8),intent(out) :: FPO4X,FPOBX,FP14X,FP1BX
  type(micforctype), intent(in) :: micfor
  type(Cumlate_Flux_Diag_type),INTENT(INOUT)::  naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                                               &
    FracOMActAutor        => nmics%FracOMActAutor,         &
    FracAutorBiomOfActK   => nmics%FracAutorBiomOfActK,    &
    AttenfNH4Autor        => nmicf%AttenfNH4Autor,         &
    AttenfNO3Autor        => nmicf%AttenfNO3Autor,         &
    AttenfH1PO4Autor      => nmicf%AttenfH1PO4Autor,       &
    AttenfH2PO4Autor      => nmicf%AttenfH2PO4Autor,       &
    RO2EcoDmndPrev        => micfor%RO2EcoDmndPrev,        &
    RNH4EcoDmndSoilPrev   => micfor%RNH4EcoDmndSoilPrev,   &
    RNH4EcoDmndBandPrev   => micfor%RNH4EcoDmndBandPrev,   &
    RNO3EcoDmndSoilPrev   => micfor%RNO3EcoDmndSoilPrev,   &
    RNO3EcoDmndBandPrev   => micfor%RNO3EcoDmndBandPrev,   &
    RH2PO4EcoDmndSoilPrev => micfor%RH2PO4EcoDmndSoilPrev, &
    RH2PO4EcoDmndBandPrev => micfor%RH2PO4EcoDmndBandPrev, &
    RH1PO4EcoDmndSoilPrev => micfor%RH1PO4EcoDmndSoilPrev, &
    RH1PO4EcoDmndBandPrev => micfor%RH1PO4EcoDmndBandPrev, &
    RNH4EcoDmndLitrPrev   => micfor%RNH4EcoDmndLitrPrev,   &
    RNO3EcoDmndLitrPrev   => micfor%RNO3EcoDmndLitrPrev,   &
    RH1PO4EcoDmndLitrPrev => micfor%RH1PO4EcoDmndLitrPrev, &
    RH2PO4EcoDmndLitrPrev => micfor%RH2PO4EcoDmndLitrPrev, &
    RDOMEcoDmndPrev       => micfor%RDOMEcoDmndPrev,       &
    RAcetateEcoDmndPrev   => micfor%RAcetateEcoDmndPrev,   &
    SoilMicPMassLayer0    => micfor%SoilMicPMassLayer0,    &
    Lsurf                 => micfor%Lsurf,                 &
    litrm                 => micfor%litrm,                 &
    ZEROS                 => micfor%ZEROS,                 &
    VLNO3                 => micfor%VLNO3,                 &
    VLNOB                 => micfor%VLNOB,                 &
    VLPO4                 => micfor%VLPO4,                 &
    VLPOB                 => micfor%VLPOB,                 &
    VLNHB                 => micfor%VLNHB,                 &
    VLNH4                 => micfor%VLNH4,                 &
    RH1PO4UptkLitrAutor   => micflx%RH1PO4UptkLitrAutor,   &
    RH2PO4UptkLitrAutor   => micflx%RH2PO4UptkLitrAutor,   &
    RNO3UptkLitrAutor     => micflx%RNO3UptkLitrAutor,     &
    RNH4UptkLitrAutor     => micflx%RNH4UptkLitrAutor,     &
    RH1PO4UptkBandAutor   => micflx%RH1PO4UptkBandAutor,   &
    RH1PO4UptkSoilAutor   => micflx%RH1PO4UptkSoilAutor,   &
    RH2PO4UptkBandAutor   => micflx%RH2PO4UptkBandAutor,   &
    RH2PO4UptkSoilAutor   => micflx%RH2PO4UptkSoilAutor,   &
    RNO3UptkBandAutor     => micflx%RNO3UptkBandAutor,     &
    RNO3UptkSoilAutor     => micflx%RNO3UptkSoilAutor,     &
    RNH4UptkBandAutor     => micflx%RNH4UptkBandAutor,     &
    RNH4UptkSoilAutor     => micflx%RNH4UptkSoilAutor,     &
    RO2DmndAutort         => micflx%RO2DmndAutort          &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(RO2EcoDmndPrev.GT.ZEROS)THEN
    FOXYX=AMAX1(FMN,RO2DmndAutort(NGL)/RO2EcoDmndPrev)
  ELSE
    FOXYX=AMAX1(FMN,FracOMActAutor(NGL))
  ENDIF
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4UptkSoilAutor(NGL)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4UptkBandAutor(NGL)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNHB)
  ENDIF
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3UptkSoilAutor(NGL)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3UptkBandAutor(NGL)/RNO3EcoDmndBandPrev)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNOB)
  ENDIF
  IF(RH2PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4UptkSoilAutor(NGL)/RH2PO4EcoDmndSoilPrev)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4UptkBandAutor(NGL)/RH2PO4EcoDmndBandPrev)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4UptkSoilAutor(NGL)/RH1PO4EcoDmndSoilPrev)
  ELSE
    FP14X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4UptkBandAutor(NGL)/RH1PO4EcoDmndBandPrev)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF
  naqfdiag%TFOXYX=naqfdiag%TFOXYX+FOXYX
  naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4X
  naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3X
  naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4X
  naqfdiag%TFP14X=naqfdiag%TFP14X+FP14X
  naqfdiag%TFNH4B=naqfdiag%TFNH4B+FNB4X
  naqfdiag%TFNO3B=naqfdiag%TFNO3B+FNB3X
  naqfdiag%TFPO4B=naqfdiag%TFPO4B+FPOBX
  naqfdiag%TFP14B=naqfdiag%TFP14B+FP1BX
!
! FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
! MICROBIAL POPULATIONS IN SURFACE RESIDUE
! F*=fraction of substrate uptake relative to total uptake from
! previous hour in surface litter, labels as for soil layers above
!
  IF(litrm)THEN
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNH4Autor(NGL)=AMAX1(FMN,RNH4UptkLitrAutor(NGL)/RNH4EcoDmndLitrPrev)
    ELSE
      AttenfNH4Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNO3Autor(NGL)=AMAX1(FMN,RNO3UptkLitrAutor(NGL)/RNO3EcoDmndLitrPrev)
    ELSE
      AttenfNO3Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH2PO4Autor(NGL)=AMAX1(FMN,RH2PO4UptkLitrAutor(NGL)/RH2PO4EcoDmndLitrPrev)
    ELSE
      AttenfH2PO4Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH1PO4Autor(NGL)=AMAX1(FMN,RH1PO4UptkLitrAutor(NGL)/RH1PO4EcoDmndLitrPrev)
    ELSE
      AttenfH1PO4Autor(NGL)=AMAX1(FMN,FracAutorBiomOfActK(NGL))
    ENDIF
  ENDIF
  IF(Lsurf.AND.SoilMicPMassLayer0.GT.ZEROS)THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+AttenfNH4Autor(NGL)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+AttenfNO3Autor(NGL)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+AttenfH2PO4Autor(NGL)
    naqfdiag%TFP14X=naqfdiag%TFP14X+AttenfH1PO4Autor(NGL)
  ENDIF
  end associate
  end subroutine SubstrateCompetitionFactorsff
!------------------------------------------------------------------------------------------

  subroutine GatherAutotrophAnabolicFlux(I,J,NGL,N,ECHZ,FGOCP,FGOAP,&
    RGrowthRespAutor,RMaintDefcitcitAutor,RMaintRespAutor,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: ECHZ
  real(r8), intent(in) :: FGOCP,FGOAP
  real(r8), intent(in) :: RGrowthRespAutor,RMaintDefcitcitAutor,RMaintRespAutor
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  integer :: M,K,MID3,MID,MID1,NE,idom
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: CXC
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FRM
!     begin_execution
  associate(                                                                    &
    rCNBiomeActAutor                 => nmics%rCNBiomeActAutor,                 &
    GrowthEnvScalAutor               => nmics%GrowthEnvScalAutor,               &
    DOMuptk4GrothAutor               => nmicf%DOMuptk4GrothAutor,               &
    NonstX2stBiomAutor               => nmicf%NonstX2stBiomAutor,               &
    RespGrossAutor                   => nmicf%RespGrossAutor,                   &
    RNOxReduxRespAutorLim            => nmicf%RNOxReduxRespAutorLim,            &
    RMaintDmndAutor                  => nmicf%RMaintDmndAutor,                  &
    RkillLitfalOMAutor               => nmicf%RkillLitfalOMAutor,               &
    RkillLitrfal2HumOMAutor          => nmicf%RkillLitrfal2HumOMAutor,          &
    RkillLitrfal2ResduOMAutor        => nmicf%RkillLitrfal2ResduOMAutor,        &
    RMaintDefcitLitrfalOMAutor       => nmicf%RMaintDefcitLitrfalOMAutor,       &
    RMaintDefcitLitrfal2HumOMAutor   => nmicf%RMaintDefcitLitrfal2HumOMAutor,   &
    RMaintDefcitLitrfal2ResduOMAutor => nmicf%RMaintDefcitLitrfal2ResduOMAutor, &
    RKillOMAutor                     => nmicf%RKillOMAutor,                     &
    RkillRecycOMAutor                => nmicf%RkillRecycOMAutor,                &
    RMaintDefcitKillOMAutor          => nmicf%RMaintDefcitKillOMAutor,          &
    RMaintDefcitRecycOMAutor         => nmicf%RMaintDefcitRecycOMAutor,         &
    Resp4NFixAutor                   => nmicf%Resp4NFixAutor,                   &
    TDOMUptkHeter                    => ncplxf%TDOMUptkHeter,                   &
    rNCOMCAutor                      => micpar%rNCOMCAutor,                     &
    rPCOMCAutor                      => micpar%rPCOMCAutor,                     &
    FL                               => micpar%FL,                              &
    ZEROS                            => micfor%ZEROS,                           &
    ZERO                             => micfor%ZERO,                            &
    mBiomeAutor                      => micstt%mBiomeAutor,                     &
    EHUM                             => micstt%EHUM                             &
  )

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

  DOMuptk4GrothAutor(idom_beg:idom_end,NGL)=0._r8
  CGOMX = AMIN1(RMaintRespAutor,RespGrossAutor(NGL))+Resp4NFixAutor(NGL)+(RGrowthRespAutor-Resp4NFixAutor(NGL))/ECHZ  
  CGOMD = RNOxReduxRespAutorLim(NGL)/ENOX  !for NO2(-) reduction by NH3
  !total C uptake, which could be CO2, or CH4, depending on the type of organism
  !for aerobic methanotrophs, the following equals to CH4 uptake for biomass
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
  MID1=micpar%get_micb_id(1,NGL);MID3=micpar%get_micb_id(3,NGL)
  IF(mBiomeAutor(ielmc,MID3).GT.ZEROS.AND.mBiomeAutor(ielmc,MID1).GT.ZEROS)THEN
    CCC=AZMAX1(AMIN1(1.0_r8 &
      ,mBiomeAutor(ielmn,MID3)/(mBiomeAutor(ielmn,MID3)+mBiomeAutor(ielmc,MID3)*rNCOMCAutor(3,NGL)) &
      ,mBiomeAutor(ielmp,MID3)/(mBiomeAutor(ielmp,MID3)+mBiomeAutor(ielmc,MID3)*rPCOMCAutor(3,NGL))))
    CXC = mBiomeAutor(ielmc,MID3)/mBiomeAutor(ielmc,MID1)
    C3C = 1.0_r8/(1.0_r8+CXC/CKC)
    CNC = AZMAX1(AMIN1(1.0_r8 &
      ,mBiomeAutor(ielmc,MID3)/(mBiomeAutor(ielmc,MID3)+mBiomeAutor(ielmn,MID3)/rNCOMCAutor(3,NGL))))
    CPC=AZMAX1(AMIN1(1.0_r8 &
      ,mBiomeAutor(ielmc,MID3)/(mBiomeAutor(ielmc,MID3)+mBiomeAutor(ielmp,MID3)/rPCOMCAutor(3,NGL))))
    RCCC = RCCZ+AMAX1(CCC,C3C)*RCCY
    RCCN = CNC*RCCX
    RCCP = CPC*RCCQ
  ELSE
    RCCC = RCCZ
    RCCN = 0.0_r8
    RCCP = 0.0_r8
  ENDIF
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
  MID3  = micpar%get_micb_id(3,NGL)
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
    ENDDO

    RkillLitfalOMAutor(ielmc,M,NGL)=RKillOMAutor(ielmc,M,NGL)*(1.0_r8-RCCC)
    RkillLitfalOMAutor(ielmn,M,NGL)=RKillOMAutor(ielmn,M,NGL)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RkillLitfalOMAutor(ielmp,M,NGL)=RKillOMAutor(ielmp,M,NGL)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
    DO NE=1,NumPlantChemElms
      RkillRecycOMAutor(NE,M,NGL)=RKillOMAutor(NE,M,NGL)-RkillLitfalOMAutor(NE,M,NGL)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P LitrFall to residue
!
      RkillLitrfal2HumOMAutor(NE,M,NGL)=AZMAX1(RkillLitfalOMAutor(NE,M,NGL)*EHUM)
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
  IF(RMaintDefcitcitAutor.GT.ZEROS.AND.RMaintRespAutor.GT.ZEROS.AND.RCCC.GT.ZERO)THEN
    FRM=RMaintDefcitcitAutor/RMaintRespAutor
    DO  M=1,2
      RMaintDefcitKillOMAutor(ielmc,M,NGL)=AMIN1(mBiomeAutor(ielmc,MID),AZMAX1(FRM*RMaintDmndAutor(M,NGL)/RCCC))
      RMaintDefcitKillOMAutor(ielmn,M,NGL)=AMIN1(mBiomeAutor(ielmn,MID),&
        AZMAX1(RMaintDefcitKillOMAutor(ielmc,M,NGL)*rCNBiomeActAutor(ielmn,NGL)))
      RMaintDefcitKillOMAutor(ielmp,M,NGL)=AMIN1(mBiomeAutor(ielmp,MID),&
        AZMAX1(RMaintDefcitKillOMAutor(ielmc,M,NGL)*rCNBiomeActAutor(ielmp,NGL)))

      RMaintDefcitLitrfalOMAutor(ielmc,M,NGL)=RMaintDefcitKillOMAutor(ielmc,M,NGL)*(1.0_r8-RCCC)
      RMaintDefcitLitrfalOMAutor(ielmn,M,NGL)=RMaintDefcitKillOMAutor(ielmn,M,NGL)*(1.0_r8-RCCN)*(1.0_r8-RCCC)
      RMaintDefcitLitrfalOMAutor(ielmp,M,NGL)=RMaintDefcitKillOMAutor(ielmp,M,NGL)*(1.0_r8-RCCP)*(1.0_r8-RCCC)
      DO NE=1,NumPlantChemElms
        RMaintDefcitRecycOMAutor(NE,M,NGL)=RMaintDefcitKillOMAutor(NE,M,NGL)-RMaintDefcitLitrfalOMAutor(NE,M,NGL)
      ENDDO
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RHMMC,RHMMN,RHMMC=transfer of senesence LitrFall C,N,P to humus
!     EHUM=humus transfer fraction
!     RCMMC,RCMMN,RCMMC=transfer of senesence LitrFall C,N,P to residue
!
      DO NE=1,NumPlantChemElms
        RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL)   = AZMAX1(RMaintDefcitLitrfalOMAutor(NE,M,NGL)*EHUM)
        RMaintDefcitLitrfal2ResduOMAutor(NE,M,NGL) = RMaintDefcitLitrfalOMAutor(NE,M,NGL)-RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL)
      ENDDO
    ENDDO
  ELSE
    DO  M=1,2
      DO NE=1,NumPlantChemElms
        RMaintDefcitKillOMAutor(NE,M,NGL)          = 0.0_r8
        RMaintDefcitLitrfalOMAutor(NE,M,NGL)       = 0.0_r8
        RMaintDefcitRecycOMAutor(NE,M,NGL)         = 0.0_r8
        RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL)   = 0.0_r8
        RMaintDefcitLitrfal2ResduOMAutor(NE,M,NGL) = 0.0_r8
      ENDDO
    ENDDO
  ENDIF
  end associate
  end subroutine GatherAutotrophAnabolicFlux
!------------------------------------------------------------------------------------------

  subroutine AerobicAutorO2Uptake(NGL,N,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
    micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL   !guild id
  integer, intent(in) :: N     !functional group id
  real(r8), intent(in) :: OXKX,FOXYX
  real(r8), intent(in) :: RGOMP,RVOXP
  real(r8), intent(in) :: RVOXPA
  real(r8), intent(in) :: RVOXPB
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx

  integer  :: M,MX
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,O2AquaDiffusvity1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMX,RO2DmndX
  real(r8) :: ROXYLX
  real(r8) :: RRADO,RMPOX,ROXDFQ
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X,VOLOXM
  real(r8) :: VOLWPM

  ! begin_execution
  associate(                                               &
    fLimO2Autor           => nmics%fLimO2Autor,            &
    OMActAutor            => nmics%OMActAutor,             &
    RO2UptkAutor          => nmicf%RO2UptkAutor,           &
    RespGrossAutor        => nmicf%RespGrossAutor,         &
    RO2Dmnd4RespAutor     => nmicf%RO2Dmnd4RespAutor,      &
    RO2DmndAutor          => nmicf%RO2DmndAutor,           &
    RO2Uptk4RespAutor     => nmicf%RO2Uptk4RespAutor,      &
    RCO2ProdAutor         => nmicf%RCO2ProdAutor,          &
    RCH4ProdAutor         => nmicf%RCH4ProdAutor,          &
    RSOxidSoilAutor       => nmicf%RSOxidSoilAutor,        &
    RSOxidBandAutor       => nmicf%RSOxidBandAutor,        &
    mid_AmmoniaOxidBacter => micpar%mid_AmmoniaOxidBacter, &
    mid_NitriteOxidBacter => micpar%mid_NitriteOxidBacter, &
    RO2GasXchangePrev     => micfor%RO2GasXchangePrev,     &
    COXYE                 => micfor%COXYE,                 &
    O2_rain_conc          => micfor%O2_rain_conc,          &
    O2_irrig_conc         => micfor%O2_irrig_conc,         &
    Irrig2LitRSurf_col    => micfor%Irrig2LitRSurf_col,    &
    Rain2LitRSurf         => micfor%Rain2LitRSurf,         &
    litrm                 => micfor%litrm,                 &
    O2AquaDiffusvity      => micfor%O2AquaDiffusvity,      &
    RO2AquaXchangePrev    => micfor%RO2AquaXchangePrev,    &
    VLSoilPoreMicP        => micfor%VLSoilPoreMicP,        &
    VLSoilMicP            => micfor%VLSoilMicP,            &
    VLsoiAirPM            => micfor%VLsoiAirPM,            &
    VLWatMicPM            => micfor%VLWatMicPM,            &
    FILM                  => micfor%FILM,                  &
    THETPM                => micfor%THETPM,                &
    TortMicPM             => micfor%TortMicPM,             &
    ZERO                  => micfor%ZERO,                  &
    ZEROS                 => micfor%ZEROS,                 &
    DiffusivitySolutEff   => micfor%DiffusivitySolutEff,   &
    O2GSolubility         => micstt%O2GSolubility,         &
    OXYG                  => micstt%OXYG,                  &
    OXYS                  => micstt%OXYS,                  &
    COXYG                 => micstt%COXYG,                 &
    RO2UptkSoilM          => micflx%RO2UptkSoilM,          &
    RNH3OxidAutor         => micflx%RNH3OxidAutor,         &
    RNH3OxidAutorBand     => micflx%RNH3OxidAutorBand,     &
    RNO2OxidAutor         => micflx%RNO2OxidAutor,         &
    RNO2OxidAutorBand     => micflx%RNO2OxidAutorBand      &
  )

  IF(RO2DmndAutor(NGL).GT.ZEROS .AND. FOXYX.GT.ZERO)THEN
    IF(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX             = RO2DmndAutor(NGL)*dts_gas    !
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
        !     RMPOX,RO2UptkSoilM=O2 uptake
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
          OXYS1=OXYS1-RMPOX
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
          RO2UptkAutor(NGL) = RO2UptkAutor(NGL)+RMPOX
          RO2UptkSoilM(M)   = RO2UptkSoilM(M)+RMPOX
        ENDDO

      ENDDO
      !write(*,*)'420'
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (OxyLimterHeter)
      !
      !     OxyLimterHeter=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter_vr,RNO2DmndReduxBandHeter_vr=NH3,NO2 oxidation in non-band, band
      !
      fLimO2Autor(NGL)=AMIN1(1.0,AZMAX1(RO2UptkAutor(NGL)/RO2DmndAutor(NGL)))
      IF(N.EQ.mid_AmmoniaOxidBacter)THEN
        RNH3OxidAutor(NGL)     = RNH3OxidAutor(NGL)*fLimO2Autor(NGL)
        RNH3OxidAutorBand(NGL) = RNH3OxidAutorBand(NGL)*fLimO2Autor(NGL)
      ELSEIF(N.EQ.mid_NitriteOxidBacter)THEN
        RNO2OxidAutor(NGL)     = RNO2OxidAutor(NGL)*fLimO2Autor(NGL)
        RNO2OxidAutorBand(NGL) = RNO2OxidAutorBand(NGL)*fLimO2Autor(NGL)
      ENDIF
    ELSE
      RO2UptkAutor(NGL) = RO2DmndAutor(NGL)
      fLimO2Autor(NGL)  = 1.0_r8
    ENDIF
  ELSE
    RO2UptkAutor(NGL) = 0.0_r8
    fLimO2Autor(NGL)  = 1.0_r8
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2X,RAcettProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
  !     RO2Uptk4RespHeter=O2-limited O2 uptake
  !     RSOxidSoilAutor,RSOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !
  RespGrossAutor(NGL)    = RGOMP*fLimO2Autor(NGL)
  RCO2ProdAutor(NGL)     = RespGrossAutor(NGL)
  RCH4ProdAutor(NGL)     = 0.0_r8
  RO2Uptk4RespAutor(NGL) = RO2Dmnd4RespAutor(NGL)*fLimO2Autor(NGL)
  RSOxidSoilAutor(NGL)   = RVOXPA*fLimO2Autor(NGL)
  RSOxidBandAutor(NGL)   = RVOXPB*fLimO2Autor(NGL)

  end associate
  end subroutine AerobicAutorO2Uptake

!------------------------------------------------------------------------------------------

  subroutine AutotrophDenitrificCatabolism(NGL,N,XCO2,VOLWZ,micfor,micstt,&
    naqfdiag,nmicf,nmics, micflx)
  !
  !2NO2(-) + NH3 -> 1.5N2O + 2OH(-) + 0.5H2O, molar based
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: VOLWZ,XCO2
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FNO2S,FNO2B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: FVMXDX
  real(r8) :: ROXYD,RDNOT
  real(r8) :: VMXDXS
  real(r8) :: VMXDXB
  real(r8) :: VMXDXT
  real(r8) :: VMXD4
  real(r8) :: VMXD4S
  real(r8) :: VMXD4B
  real(r8) :: ZNO2SX,ZNO2BX

!     begin_execution
  associate(                                              &
    FracNO2ReduxAutor     => nmics%FracNO2ReduxAutor,     &
    RO2Dmnd4RespAutor     => nmicf%RO2Dmnd4RespAutor,     &
    RO2Uptk4RespAutor     => nmicf%RO2Uptk4RespAutor,     &
    RNO3UptkAutor         => nmicf%RNO3UptkAutor,         &
    RNO2ReduxAutorSoil    => nmicf%RNO2ReduxAutorSoil,    &
    RNO2ReduxAutorBand    => nmicf%RNO2ReduxAutorBand,    &
    RNOxReduxRespAutorLim => nmicf%RNOxReduxRespAutorLim, &
    RSOxidSoilAutor       => nmicf%RSOxidSoilAutor,       &
    RSOxidBandAutor       => nmicf%RSOxidBandAutor,       &
    RTotNH3OxidSoilAutor  => nmicf%RTotNH3OxidSoilAutor,  &
    RTotNH3OxidBandAutor  => nmicf%RTotNH3OxidBandAutor,  &
    RNO2EcoUptkSoilPrev   => micfor%RNO2EcoUptkSoilPrev,  &
    VLNO3                 => micfor%VLNO3,                &
    VLNOB                 => micfor%VLNOB,                &
    RNO2EcoUptkBandPrev   => micfor%RNO2EcoUptkBandPrev,  &
    ZEROS                 => micfor%ZEROS,                &
    ZEROS2                => micfor%ZEROS2,               &
    CNO2B                 => micstt%CNO2B,                &
    CNO2S                 => micstt%CNO2S,                &
    ZNO2B                 => micstt%ZNO2B,                &
    ZNO2S                 => micstt%ZNO2S,                &
    RNO2OxidAutor         => micflx%RNO2OxidAutor,        &
    RNO2OxidAutorBand     => micflx%RNO2OxidAutorBand     &
  )
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
  IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RNO2OxidAutor(NGL)/RNO2EcoUptkSoilPrev)
  ELSE
    FNO2=AMAX1(FMN,FracNO2ReduxAutor(NGL)*VLNO3)
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2OxidAutorBand(NGL)/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=AMAX1(FMN,FracNO2ReduxAutor(NGL)*VLNOB)
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2
!
!     ROXYD=O2 demand RO2Dmnd4RespHeter not met by O2 uptake RO2Uptk4RespHeter
!     VMXD4=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM=Km for NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD4S,VMXD4B=substrate-unlimited NO2 reduction in non-band,band
!     RNO2ReduxHeterSoil,RNO2ReduxHeterBand=substrate-limited NO2 reduction in non-band,band
!     RNOxReduxRespDenitUlm,RNOxReduxRespDenitLim=total substrate-unltd,-ltd respn from NO2 reduction
!     ECNO=efficiency CO2 conversion to biomass
!     ECHZ=growth respiration efficiency
!     RSOxidSoilAutor,RSOxidBandAutor=total O2-limited (1)NH4,(2)NO2,(3)CH4 oxidation
!
  FNO2S=VLNO3
  FNO2B=VLNOB
  ROXYD=AZMAX1(RO2Dmnd4RespAutor(NGL)-RO2Uptk4RespAutor(NGL))
  VMXD4=0.875_r8*ROXYD*XCO2
  VMXDXS=FNO2S*VMXD4*CNO2S/(CNO2S+Z2KM)
  VMXDXB=FNO2B*VMXD4*CNO2B/(CNO2B+Z2KM)
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD4S=VMXDXS*FVMXDX
  VMXD4B=VMXDXB*FVMXDX
  ZNO2SX=ZNO2S+RTotNH3OxidSoilAutor
  ZNO2BX=ZNO2B+RTotNH3OxidBandAutor
  RNO2ReduxAutorSoil(NGL)=AZMAX1(AMIN1(VMXD4S,ZNO2SX))
  RNO2ReduxAutorBand(NGL)=AZMAX1(AMIN1(VMXD4B,ZNO2BX))
  RDNOT=RNO2ReduxAutorSoil(NGL)+RNO2ReduxAutorBand(NGL)
  RNOxReduxRespAutorLim(NGL)=RDNOT*ECNO*ENOX
  RNO3UptkAutor(NGL)=0.0_r8
  RNO2OxidAutor(NGL)=VMXD4S
  RNO2OxidAutorBand(NGL)=VMXD4B
  !NH4 oxidation by NO2(-)
  RSOxidSoilAutor(NGL)=RSOxidSoilAutor(NGL)+0.333_r8*RNO2ReduxAutorSoil(NGL)
  RSOxidBandAutor(NGL)=RSOxidBandAutor(NGL)+0.333_r8*RNO2ReduxAutorBand(NGL)

!     TRN2ON=TRN2ON+RNO2ReduxAutorSoil(NGL)+RNO2ReduxAutorBand(NGL)
  end associate
  end subroutine AutotrophDenitrificCatabolism
!------------------------------------------------------------------------------------------
  subroutine NH3OxidizerCatabolism(NGL,N,XCO2,VOLWZ,TSensGrowth,ECHZ,RGOMP,RVOXP,RVOXPA,&
    RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N
  REAL(R8), INTENT(IN) :: VOLWZ,TSensGrowth
  REAL(R8), INTENT(IN) :: XCO2
  real(r8),intent(out) :: ECHZ
  real(r8),intent(out) :: RGOMP,RVOXP
  real(r8),intent(out) :: RVOXPA,RVOXPB
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type),INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  real(r8) :: FNH4S,FNHBS
  real(r8) :: FNH4,FNB4
  real(r8) :: FCN4S,FCN4B
  real(r8) :: RNNH4,RNNHB
  real(r8) :: VMXX,VMX4S
  real(r8) :: VMX4B
  real(r8) :: ZNFN4S,ZNFN4B
  real(r8) :: FSBST
  real(r8) :: VMXA

!     begin_execution
  associate(                                            &
    GrowthEnvScalAutor   => nmics%GrowthEnvScalAutor,   &
    FBiomStoiScalarAutor => nmics%FBiomStoiScalarAutor, &
    FracOMActAutor       => nmics%FracOMActAutor,       &
    OMActAutor           => nmics%OMActAutor,           &
    RO2Dmnd4RespAutor    => nmicf%RO2Dmnd4RespAutor,    &
    RO2DmndAutor         => nmicf%RO2DmndAutor,         &
    VLNH4                => micfor%VLNH4,               &
    VLNHB                => micfor%VLNHB,               &
    ZEROS                => micfor%ZEROS,               &
    ZEROS2               => micfor%ZEROS2,              &
    RNH4EcoDmndSoilPrev  => micfor%RNH4EcoDmndSoilPrev, &
    RNH4EcoDmndBandPrev  => micfor%RNH4EcoDmndBandPrev, &
    ZNFN0                => micstt%ZNFN0,               &
    ZNFNI                => micstt%ZNFNI,               &
    CNH4S                => micstt%CNH4S,               &
    CNH4B                => micstt%CNH4B,               &
    ZNH4S                => micstt%ZNH4S,               &
    ZNH4B                => micstt%ZNH4B,               &
    RNH3OxidAutor        => micflx%RNH3OxidAutor,       &
    RNH3OxidAutorBand    => micflx%RNH3OxidAutorBand,   &
    RO2DmndAutort        => micflx%RO2DmndAutort        &
  )
!
!     FACTOR TO REGULATE COMPETITION FOR NH4 AMONG DIFFERENT
!     MICROBIAL AND ROOT POPULATIONS FNH4
!
!     FNH4,FNB4=frac of total biol demand for NH4 in non-band, band
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4=AMAX1(FMN,RNH3OxidAutor(NGL)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4=AMAX1(FMN,VLNH4*FracOMActAutor(NGL))
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4=AMAX1(FMN,RNH3OxidAutorBand(NGL)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4=AMAX1(FMN,VLNHB*FracOMActAutor(NGL))
  ENDIF
  naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4
  naqfdiag%TFNH4B=naqfdiag%TFNH4B+FNB4
!
!     NITRIFICATION INHIBITION
!
!     ZNFN0=inhibition when fertilizer added
!     ZNFNI=reduction in inhibition since fertilizer added
!     CNH4S,CNH4B=NH4 concentrations in non-band, band
!     TSensGrowth=temperature effect
!     RNFNI=rate constant for inhibition decline
!     ZHKI=inhibition from high CNH4
!     ZNFN4S,ZNFN4B=inhibition in non-band, band
!
  IF(ZNFN0.GT.ZEROS)THEN
    ZNFNI=ZNFNI*(1.0_r8-RNFNI*TSensGrowth)
    ZNFN4S=ZNFN0-ZNFNI/(1.0_r8+CNH4S/ZHKI)
    ZNFN4B=ZNFN0-ZNFNI/(1.0_r8+CNH4B/ZHKI)
  ELSE
    ZNFN4S=1.0_r8
    ZNFN4B=1.0_r8
  ENDIF
!
!     NH3 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     NH3 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXX=potential NH3 oxidation, VMXH=specific oxidation
!     TFNG=temperature+water limitation, FBiomStoiScalarAutorr=N,P limitation
!     XCO2=aqueous CO2 limitation, OMA=active biomass
!     VMXA= non-substrate limited NH3 oxidation
!     VHKI=nonlinear increase in VMXA with VMXH
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     CNH4S,CNH4B=NH4 concentration in non-band, band
!     ZHKM=Km for NH4 uptake
!     FNH4,FNB4=fractions of total NH4 demand in non-band, band
!     ZNH4S,ZNH4B=NH4 amount in non-band, band
!     RNNH4,RNNHB=NH3 oxidation in non-band, band
!     RGOMP=O2-unlimited respiration
!     ECNH=efficiency CO2 conversion to biomass
!     RVMX4,RVMXB=nitrifier demand for NH4 in non-band, band
!
  ECHZ=EO2X
  VMXX=VMXH*GrowthEnvScalAutor(NGL)*FBiomStoiScalarAutor(NGL)*XCO2*OMActAutor(NGL)
  IF(VOLWZ.GT.ZEROS2)THEN
    VMXA=VMXX/(1.0_r8+VMXX/(VHKI*VOLWZ))
  ELSE
    VMXA=0.0_r8
  ENDIF
  FCN4S                  = FNH4S*CNH4S/(CNH4S+ZHKM)
  FCN4B                  = FNHBS*CNH4B/(CNH4B+ZHKM)
  FSBST                  = FCN4S+FCN4B
  VMX4S                  = VMXA*FCN4S
  VMX4B                  = VMXA*FCN4B
  RNNH4                  = AZMAX1(AMIN1(VMX4S,FNH4*ZNH4S))*ZNFN4S
  RNNHB                  = AZMAX1(AMIN1(VMX4B,FNB4*ZNH4B))*ZNFN4B
  RVOXP                  = RNNH4+RNNHB
  RVOXPA                 = RNNH4
  RVOXPB                 = RNNHB
  RGOMP                  = AZMAX1(RVOXP*ECNH*ECHZ)
  RNH3OxidAutor(NGL)     = VMX4S
  RNH3OxidAutorBand(NGL) = VMX4B
!
!     O2 DEMAND FROM NH3 OXIDATION
!
!     RO2Dmnd4RespHeter=O2 demand from respiration by nitrifiers
!     ROXYP,RO2Dmnd4RespHeter=O2 demand from respiration + NH3 oxidation
! C+O2 -> CO2,  respiration, 2.667=32./12.
! NH3+1.5O2-> NO2(-)+H2O+H(+), 1.5*32/14.=3.249

  RO2Dmnd4RespAutor(NGL)=2.667_r8*RGOMP
  RO2DmndAutor(NGL)=RO2Dmnd4RespAutor(NGL)+3.429_r8*RVOXP
  RO2DmndAutort(NGL)=RO2DmndAutor(NGL)
!
  end associate
  end subroutine NH3OxidizerCatabolism
!------------------------------------------------------------------------------------------
  subroutine NO2OxidizerCatabolism(NGL,N,XCO2,ECHZ,RGOMP,RVOXP,&
    RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  !
  !nitrite oxidation
  !NO2(-) + 0.5O2 -> NO3(-)   
  implicit none
  integer, intent(in) :: NGL,N
  REAL(r8), intent(in) :: XCO2
  real(r8), intent(out) :: ECHZ,RGOMP,RVOXP
  real(r8), intent(out) :: RVOXPA,RVOXPB
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  real(r8) :: FNH4S,FNHBS
  REAL(R8) :: fno2,FNB2
  real(r8) :: FCN2S,FCN2B
  real(r8) :: RNNO2,RNNOB
  real(r8) :: VMX2S,VMX2B
  real(r8) :: FSBST
  real(r8) :: VMXA

!     begin_execution
  associate(                                            &
    GrowthEnvScalAutor   => nmics%GrowthEnvScalAutor,   &
    FBiomStoiScalarAutor => nmics%FBiomStoiScalarAutor, &
    FracNO2ReduxAutor    => nmics%FracNO2ReduxAutor,    &
    OMActAutor           => nmics%OMActAutor,           &
    RO2Dmnd4RespAutor    => nmicf%RO2Dmnd4RespAutor,    &
    RO2DmndAutor         => nmicf%RO2DmndAutor,         &
    VLNH4                => micfor%VLNH4,               &
    VLNHB                => micfor%VLNHB,               &
    VLNO3                => micfor%VLNO3,               &
    VLNOB                => micfor%VLNOB,               &
    ZEROS                => micfor%ZEROS,               &
    RNO2EcoUptkSoilPrev  => micfor%RNO2EcoUptkSoilPrev, &
    RNO2EcoUptkBandPrev  => micfor%RNO2EcoUptkBandPrev, &
    CNO2S                => micstt%CNO2S,               &
    CNO2B                => micstt%CNO2B,               &
    ZNO2S                => micstt%ZNO2S,               &
    ZNO2B                => micstt%ZNO2B,               &
    RNO2OxidAutor        => micflx%RNO2OxidAutor,       &
    RNO2OxidAutorBand    => micflx%RNO2OxidAutorBand,   &
    RO2DmndAutort        => micflx%RO2DmndAutort        &
  )
!     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
!     MICROBIAL POPULATIONS
!
!     FNO2=fraction of total biological demand for NO2 in non-band, band
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RNO2OxidAutor(NGL)/RNO2EcoUptkSoilPrev)
  ELSE
    FNO2=AMAX1(FMN,FracNO2ReduxAutor(NGL)*VLNO3)
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2OxidAutorBand(NGL)/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=AMAX1(FMN,FracNO2ReduxAutor(NGL)*VLNOB)
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     NO2 CONCENTRATIONS
!
!     ECHZ=growth respiration efficiency
!     VMXA= non-substrate limited NH3 oxidation
!     VMXN=specific oxidation
!     TFNG=temperature+water limitation, FBiomStoiScalarAutorr=N,P limitation
!     XCO2=aqueous CO2 limitation, OMA=active biomass
!     OMA=active biomass
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     CNO2S,CNO2B=NO2 concentration in non-band, band
!     ZNKM=Km for NO2 uptake
!     FNO2,FNB2=fractions of total NO2 demand in non-band, band
!     ZNO2S,ZNO2B=NO2 amount in non-band, band
!     RNNO2,RNNOB=NO2 oxidation in non-band, band
!     RGOMP=O2-unlimited respiration
!     ECNO=efficiency CO2 conversion to biomass
!     RNO2DmndReduxSoilHeter_vr,RNO2DmndReduxBandHeter_vr=nitrifier demand for NO2 in non-band, band
!
  ECHZ=EO2X
  VMXA=GrowthEnvScalAutor(NGL)*FBiomStoiScalarAutor(NGL)*XCO2*OMActAutor(NGL)*VMXN
  FCN2S                  = FNH4S*CNO2S/(CNO2S+ZNKM)
  FCN2B                  = FNHBS*CNO2B/(CNO2B+ZNKM)
  FSBST                  = FCN2S+FCN2B
  VMX2S                  = VMXA*FCN2S
  VMX2B                  = VMXA*FCN2B
  RNNO2                  = AZMAX1(AMIN1(VMX2S,FNO2*ZNO2S))
  RNNOB                  = AZMAX1(AMIN1(VMX2B,FNB2*ZNO2B))
  RVOXP                  = RNNO2+RNNOB
  RVOXPA                 = RNNO2
  RVOXPB                 = RNNOB
  RGOMP                  = AZMAX1(RVOXP*ECNO*ECHZ)
  RNO2OxidAutor(NGL)     = VMX2S
  RNO2OxidAutorBand(NGL) = VMX2B
!
!     O2 DEMAND FROM NO2 OXIDATION
!
!     RO2Dmnd4RespHeter=O2 demand from respiration by nitrifiers
!     ROXYP,RO2Dmnd4RespHeter=O2 demand from respiration + NO2 oxidation
!
  RO2Dmnd4RespAutor(NGL)=2.667_r8*RGOMP
  RO2DmndAutor(NGL)=RO2Dmnd4RespAutor(NGL)+1.143_r8*RVOXP
  RO2DmndAutort(NGL)=RO2DmndAutor(NGL)

  end associate
  end subroutine NO2OxidizerCatabolism
!------------------------------------------------------------------------------------------


  subroutine H2MethanogensCatabolism(NGL,N,ECHZ,RGOMP,XCO2,micfor,micstt,&
    naqfdiag,nmicf,nmics,micflx)
  !
  !Hydrogenotrophic CH4 production
  !CO2 + 4H2 -> CH4 + 2H2O
  !H2 is produced from fermentation   

  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(out) :: ECHZ
  real(r8), intent(out) :: RGOMP
  REAL(R8), INTENT(IN) :: XCO2
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  real(r8) :: GH2X,GH2H
  real(r8) :: H2GSX
  real(r8) :: FSBST
  real(r8) :: VMXA

  associate(                                            &
    GrowthEnvScalAutor   => nmics%GrowthEnvScalAutor,   &
    FBiomStoiScalarAutor => nmics%FBiomStoiScalarAutor, &
    OMActAutor           => nmics%OMActAutor,           &
    RO2Dmnd4RespAutor    => nmicf%RO2Dmnd4RespAutor,    &
    RO2DmndAutor         => nmicf%RO2DmndAutor,         &
    TKS                  => micfor%TKS,                 &
    CH2GS                => micstt%CH2GS,               &
    H2GS                 => micstt%H2GS,                &
    RO2DmndAutort        => micflx%RO2DmndAutort        &
  )
!     begin_execution
!
!     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
!
!     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
!     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
!     VMXA=substrate-unlimited H2 oxidation rate
!     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (tCResp4H2Prod)
!     CH2GS=H2 concentration, H2KM=Km for H2 uptake
!     RGOMP=H2 oxidation, ROXY*=O2 demand
!
!     and energy yield of hydrogenotrophic
!     methanogenesis GH2X at ambient H2 concentration CH2GS

  GH2X=RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**4)
  GH2H=GH2X/12.08_r8
  ECHZ=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCOX+GH2H))/EOMH)))
  VMXA=GrowthEnvScalAutor(NGL)*FBiomStoiScalarAutor(NGL)*XCO2*OMActAutor(NGL)*VMXC
  H2GSX=H2GS+0.111_r8*naqfdiag%tCResp4H2Prod
  FSBST=CH2GS/(CH2GS+H2KM)
  !why 1.5? 
  RGOMP=AZMAX1(AMIN1(1.5_r8*H2GSX,VMXA*FSBST))
  RO2Dmnd4RespAutor(NGL)=0.0_r8
  RO2DmndAutor(NGL)=0.0_r8
  RO2DmndAutort(NGL)=0.0_r8
  naqfdiag%tCH4ProdH2=naqfdiag%tCH4ProdH2+RGOMP
!
  end associate
  end subroutine H2MethanogensCatabolism
!------------------------------------------------------------------------------------------

  subroutine MethanotrophCatabolism(I,J,NGL,N,ECHZ,RGOMP,&
    RCH4Oxid,RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NGL,N
  real(r8), intent(out) :: ECHZ
  real(r8), intent(out) :: RGOMP     !methane oxidized into CH2O
  real(r8), intent(out) :: RCH4Oxid
  real(r8), intent(out) :: RVOXPA    !methane oxidation
  real(r8), intent(out) :: RVOXPB
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(Cumlate_Flux_Diag_type), intent(in) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  integer  :: M,MM
  real(r8) :: CH4G1,CH4S1,CCH4S1
  real(r8) :: RCH4L1,RCH4F1,RCH4S1
  real(r8) :: RGOMP1,RCHDF
  real(r8) :: VMXA1,VOLWCH
  real(r8) :: FSBST
  real(r8) :: RVOXP1
  real(r8) :: VMXA
  real(r8) :: pscal
  REAL(R8) :: VOLWPM

  associate(                                            &
    GrowthEnvScalAutor   => nmics%GrowthEnvScalAutor,   &
    FBiomStoiScalarAutor => nmics%FBiomStoiScalarAutor, &
    OMActAutor           => nmics%OMActAutor,           &
    RO2Dmnd4RespAutor    => nmicf%RO2Dmnd4RespAutor,    &
    RO2DmndAutor         => nmicf%RO2DmndAutor,         &
    CCH4E                => micfor%CCH4E,               &
    VLsoiAirPM           => micfor%VLsoiAirPM,          &
    VLWatMicPM           => micfor%VLWatMicPM,          &
    ZEROS2               => micfor%ZEROS2,              &
    ZEROS                => micfor%ZEROS,               &
    THETPM               => micfor%THETPM,              &
    DiffusivitySolutEff  => micfor%DiffusivitySolutEff, &
    litrm                => micfor%litrm,               &
    CCH4G                => micstt%CCH4G,               &
    CH4S                 => micstt%CH4S,                &
    CH4AquaSolubility    => micstt%CH4AquaSolubility,   &
    RCH4PhysexchPrev     => micfor%RCH4PhysexchPrev,    &
    RCH4GasXchangePrev   => micfor%RCH4GasXchangePrev,  &
    RO2DmndAutort        => micflx%RO2DmndAutort        &
  )
!     begin_execution
!
!     CH4 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     CH4 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXA=potential oxidation
!     RCH4PhysexchPrev=total aqueous CH4 exchange from previous hour
!     RCH4GasXchangePrev=total gaseous CH4 exchange from previous hour
!     tCH4ProdAceto+tCH4ProdH2=total CH4 generated from methanogenesis
!     dts_gas=1.0/(NPH*NPT)
!     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
!     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
!     VLsoiAirPM,VLWatMicPM=air,water-filled porosity
!     CH4AquaSolubility=CH4 aqueous solubility
!     CCK4=Km for CH4 uptake
!     ECHO=efficiency CO2 conversion to biomass
!     RGOMP1=substrate-limited CH4 oxidation
!     RCHDF=gaseous-aqueous CH4 exchange
!     DiffusivitySolutEff=rate constant for gaseous-aqueous exchange
!
  ECHZ   = EH4X
  VMXA   = GrowthEnvScalAutor(NGL)*FBiomStoiScalarAutor(NGL)*OMActAutor(NGL)*VMX4
  RCH4L1 = RCH4PhysexchPrev*dts_gas
  RCH4F1 = RCH4GasXchangePrev*dts_gas
  RCH4S1 = (naqfdiag%tCH4ProdAceto+naqfdiag%tCH4ProdH2)*dts_gas

  IF(litrm)THEN
    !surface residue layer
    CH4G1 = CCH4E*VLsoiAirPM(1)
    VMXA1 = AZMAX1(AMIN1(VMXA*dts_gas,CH4G1))  !apparent vmax for uptake
  ELSE
    CH4G1=CCH4G*VLsoiAirPM(1)
    VMXA1    = VMXA*dts_gas  !apparent vmax for uptake    
  ENDIF
  CH4S1    = CH4S
  RCH4Oxid = 0.0_r8
  RGOMP    = 0.0_r8
!
!     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
!     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
! for aerobic methanotrophs, CH4 is oxiized to CO2 for energy, and also 
! to intracellular C for respiration. (Grant 1999), the C yield is approximated
! as the energy required for turning CH4 into organic C, and the energy released
! from turnining CH4 into CO2. 
! the catabolic reaction is (molar-basis)
!  CH4 + 2O2 -> CO2 + 2H2O
!  mass basis becomes (2*32/12=5.33)
!  CH4 + 5.33 O2 -> CO2 + 2H2O  
!  
!  C+ O2 -> CO2,  32/12=2.667
  D320: DO M=1,NPH
    IF(VLWatMicPM(M).GT.ZEROS2)THEN
      VOLWCH = VLWatMicPM(M)*CH4AquaSolubility
      VOLWPM = VOLWCH+VLsoiAirPM(M)

      !CH4 uptake by aerobic oxidation
      !RCH4F1: net gaseous CH4 flux from transport into the layer
      !RCH4L1: net aqueous CH4 flux from transport into the layer
      !RCH4S1: aquoues CH4 production from methanogenesis in the layer
      D325: DO MM=1,NPT
        CH4G1  = CH4G1+RCH4F1
        CH4S1  = CH4S1+RCH4L1+RCH4S1
        CCH4S1 = AZMAX1(safe_adb(CH4S1,VLWatMicPM(M)))
        FSBST  = CCH4S1/(CCH4S1+CCK4)
        !RVOXP1 energy from oxidizing CH4 into CO2
        RVOXP1=AMIN1(AZMAX1(CH4S1)/(1.0_r8+ECHO*ECHZ),VMXA1*FSBST)
        !the respiration yield of CH2O, CH4+O2 -> CH2O + H2O (molar basis)
        RGOMP1=RVOXP1*ECHO*ECHZ
        !
        if(abs(RVOXP1+RGOMP1)>0._r8)then
          pscal=CH4S1/(RVOXP1+RGOMP1)
          if(pscal < 1._r8)then
            RVOXP1=RVOXP1*pscal
            RGOMP1=RGOMP1*pscal
          endif
          CH4S1=CH4S1-RVOXP1-RGOMP1
        endif

        !dissolution-vaporization
        IF(THETPM(M).GT.AirFillPore_Min)THEN
          RCHDF=DiffusivitySolutEff(M)*(AMAX1(ZEROS,CH4G1)*VOLWCH-CH4S1*VLsoiAirPM(M))/VOLWPM
          RCHDF=AMAX1(AMIN1(CH4G1,RCHDF),-CH4S1)
        ELSE
          RCHDF=0.0_r8
        ENDIF
        CH4G1    = CH4G1-RCHDF
        CH4S1    = CH4S1+RCHDF
        RCH4Oxid = RCH4Oxid+RVOXP1
        RGOMP    = RGOMP+RGOMP1
      ENDDO D325
    ENDIF
  ENDDO D320
  RVOXPA = AZMAX1(RCH4Oxid)
  RVOXPB = 0.0_r8
!
!     O2 DEMAND FROM CH4 OXIDATION
!
!     RO2Dmnd4RespHeter=O2 demand from respiration
!     ROXYP=O2 demand from respiration + CH4 oxidation
! 
!  IF(litrm)THEN
!  write(114,*)I+J/24.,RGOMP,OMActAutor(NGL)
!  if(OMActAutor(NGL)>1.e2)stop
!  endif
  RO2Dmnd4RespAutor(NGL) = 2.667_r8*RGOMP
  RO2DmndAutor(NGL)      = RO2Dmnd4RespAutor(NGL)+5.333_r8*RCH4Oxid
  RO2DmndAutort(NGL)     = RO2DmndAutor(NGL)
  end associate
  end subroutine MethanotrophCatabolism
!------------------------------------------------------------------------------------------
  subroutine BiomassMineralizationff(NGL,N,FNH4X,&
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: FNH4X
  real(r8), intent(in) :: FNB3X,FNB4X,FNO3X
  real(r8), intent(in) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
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
  integer :: MID3

!     begin_execution
  associate(                                             &
   GrowthEnvScalAutor    => nmics%GrowthEnvScalAutor,    &
   OMActAutor            => nmics%OMActAutor,            &
   AttenfNH4Autor        => nmicf%AttenfNH4Autor,        &
   AttenfNO3Autor        => nmicf%AttenfNO3Autor,        &
   AttenfH2PO4Autor      => nmicf%AttenfH2PO4Autor,      &
   AttenfH1PO4Autor      => nmicf%AttenfH1PO4Autor,      &
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
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  MID3=micpar%get_micb_id(3,NGL)
  RNetNH4MinPotent=mBiomeAutor(ielmc,MID3)*rNCOMCAutor(3,NGL)-mBiomeAutor(ielmn,MID3)
  IF(RNetNH4MinPotent.GT.0.0_r8)THEN
    CNH4X=AZMAX1(CNH4S-Z4MN)
    CNH4Y=AZMAX1(CNH4B-Z4MN)
    RINHX=AMIN1(RNetNH4MinPotent,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*Z4MX)
    RNH4UptkSoilAutor(NGL)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RNH4UptkBandAutor(NGL)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M=Z4MN*VLWatMicP*FNH4S
    ZNHBM=Z4MN*VLWatMicP*FNHBS
    RNH4TransfSoilAutor(NGL)=AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4UptkSoilAutor(NGL))
    RNH4TransfBandAutor(NGL)=AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4UptkBandAutor(NGL))
  ELSE
    RNH4UptkSoilAutor(NGL)=0.0_r8
    RNH4UptkBandAutor(NGL)=0.0_r8
    RNH4TransfSoilAutor(NGL)=RNetNH4MinPotent*FNH4S
    RNH4TransfBandAutor(NGL)=RNetNH4MinPotent*FNHBS
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
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     NetNH4Mineralize=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
  FNO3S=VLNO3
  FNO3B=VLNOB
  RNetNO3Dmnd=AZMAX1(RNetNH4MinPotent-RNH4TransfSoilAutor(NGL)-RNH4TransfBandAutor(NGL))
  IF(RNetNO3Dmnd.GT.0.0)THEN
    CNO3X=AZMAX1(CNO3S-ZOMN)
    CNO3Y=AZMAX1(CNO3B-ZOMN)
    RINOX=AMIN1(RNetNO3Dmnd,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*ZOMX)
    RNO3UptkSoilAutor(NGL)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RNO3UptkBandAutor(NGL)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M=ZOMN*VLWatMicP*FNO3S
    ZNOBM=ZOMN*VLWatMicP*FNO3B
    RNO3TransfSoilAutor(NGL)=AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)) &
      ,RNO3UptkSoilAutor(NGL))
    RNO3TransfBandAutor(NGL)=AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)) &
      ,RNO3UptkBandAutor(NGL))
  ELSE
    RNO3UptkSoilAutor(NGL)=0.0_r8
    RNO3UptkBandAutor(NGL)=0.0_r8
    RNO3TransfSoilAutor(NGL)=RNetNO3Dmnd*FNO3S
    RNO3TransfBandAutor(NGL)=RNetNO3Dmnd*FNO3B
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
!     RH2PO4TransfSoilHeter,RH2PO4TransfBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
  FH2PS=VLPO4
  FH2PB=VLPOB
  MID3=micpar%get_micb_id(3,NGL)
  RNetH2PO4MinPotent=(mBiomeAutor(ielmc,MID3)*rPCOMCAutor(3,NGL)-mBiomeAutor(ielmp,MID3))
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
!     RH1PO4TransfSoilHeter,RH1PO4TransfBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
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
!     RNH4TransfLitrHeter=substrate-limited NH4 mineraln-immobiln
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
!     RNO3TransfLitrHeter=substrate-unlimited NO3 immobiln
!     VLWatMicP=water content
!     ZNO3M=NO3 not available for uptake
!     AttenfNO3Heter=fraction of biological NO3 demand
!     RNO3TransfLitrHeter=substrate-limited NO3 immobiln
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
!     RH2PO4TransfLitrHeter=substrate-limited H2PO4 mineraln-immobiln
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
!     RH1PO4TransfLitrHeter=substrate-limited HPO4 minereraln-immobiln
!     NetPO4Mineralize=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RNetH1PO4DmndLitr=0.1*AZMAX1(RNetH2PO4MinPotentLitr-RH2PO4TransfLitrAutor(NGL))
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
  end associate
  end subroutine BiomassMineralizationff
!------------------------------------------------------------------------------------------

  subroutine GatherAutotrophRespiration(NGL,N,RMOMK,micfor,micstt,RGrowthRespAutor,&
    RMaintDefcitcitAutor,RMaintRespAutor,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  real(r8), intent(out) :: RGrowthRespAutor,RMaintDefcitcitAutor
  real(r8), intent(out) :: RMaintRespAutor
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
  integer :: MID1
!     begin_execution
  associate(                                    &
    TSensMaintRAutor => nmics%TSensMaintRAutor, &
    OMN2Autor        => nmics%OMN2Autor,        &
    RMaintDmndAutor  => nmicf%RMaintDmndAutor,  &
    RespGrossAutor   => nmicf%RespGrossAutor,   &
    Resp4NFixAutor   => nmicf%Resp4NFixAutor,   &
    RN2FixAutor      => nmicf%RN2FixAutor,      &
    pH               => micfor%pH,              &
    mBiomeAutor      => micstt%mBiomeAutor      &

  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TempMaintRHeter=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  MID1=micpar%get_micb_id(1,NGL)
  FPH=1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX=RMOM*TSensMaintRAutor(NGL)*FPH
  RMaintDmndAutor(1,NGL)=mBiomeAutor(ielmn,MID1)*RMOMX*RMOMK(1)
  RMaintDmndAutor(2,NGL)=OMN2Autor(NGL)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMaintRespAutor=total maintenance respiration
!     RGrowthRespAutor=growth respiration
!     RMaintDefcitcitAutor=senescence respiration
!
  RMaintRespAutor=RMaintDmndAutor(1,NGL)+RMaintDmndAutor(2,NGL)
  RGrowthRespAutor=AZMAX1(RespGrossAutor(NGL)-RMaintRespAutor)
  RMaintDefcitcitAutor=AZMAX1(RMaintRespAutor-RespGrossAutor(NGL))

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

  RN2FixAutor(NGL)=0.0_r8
  Resp4NFixAutor(NGL)=0.0_r8
  end associate
  end subroutine GatherAutotrophRespiration

!------------------------------------------------------------------------------------------

  subroutine AutotrophAnabolicUpdate(micfor,micstt,nmicf)

  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  real(r8) :: CGROMC
  integer :: N,M,NGL,MID,MID3,NE
  associate(                                                                &
    DOMuptk4GrothAutor             => nmicf%DOMuptk4GrothAutor,             &
    NonstX2stBiomAutor             => nmicf%NonstX2stBiomAutor,             &
    Resp4NFixAutor                 => nmicf%Resp4NFixAutor,                 &
    RespGrossAutor                 => nmicf%RespGrossAutor,                 &
    RNOxReduxRespAutorLim          => nmicf%RNOxReduxRespAutorLim,          &
    RNO3TransfSoilAutor            => nmicf%RNO3TransfSoilAutor,            &
    RCO2ProdAutor                  => nmicf%RCO2ProdAutor,                  &
    RH2PO4TransfSoilAutor          => nmicf%RH2PO4TransfSoilAutor,          &
    RNH4TransfBandAutor            => nmicf%RNH4TransfBandAutor,            &
    RNO3TransfBandAutor            => nmicf%RNO3TransfBandAutor,            &
    RH2PO4TransfBandAutor          => nmicf%RH2PO4TransfBandAutor,          &
    RkillLitrfal2HumOMAutor        => nmicf%RkillLitrfal2HumOMAutor,        &
    RMaintDefcitLitrfal2HumOMAutor => nmicf%RMaintDefcitLitrfal2HumOMAutor, &
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
    JGniA                          => micpar%JGniA,                         &
    JGnfA                          => micpar%JGnfA,                         &
    NumMicbFunGrupsPerCmplx        => micpar%NumMicbFunGrupsPerCmplx,       &
    icarbhyro                      => micpar%icarbhyro,                     &
    iprotein                       => micpar%iprotein,                      &
    k_POM                          => micpar%k_POM,                         &
    is_activeMicrbFungrpAutor      => micpar%is_activeMicrbFungrpAutor      &
  )
  DO  N=1,NumMicbFunGrupsPerCmplx
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
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL))
              SolidOM(NE,icarbhyro,k_POM)=SolidOM(NE,icarbhyro,k_POM)+ElmAllocmatMicrblitr2POM(2)&
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL))
            ENDDO
          ELSE
            DO NE=1,NumPlantChemElms
              SOMHumProtein(NE)=SOMHumProtein(NE)+ElmAllocmatMicrblitr2POMU(1) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL))
              SOMHumCarbohyd(NE)=SOMHumCarbohyd(NE)+ElmAllocmatMicrblitr2POMU(2) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefcitLitrfal2HumOMAutor(NE,M,NGL))
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
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     RH2PO4TransfSoilHeter,RH2PO4TransfBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     RH1PO4TransfSoilHeter,RH1PO4TransfBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
!     RNH4TransfLitrHeter,RNO3TransfLitrHeter =substrate-limited NH4,NO3 mineraln-immobiln
!     RH2PO4TransfLitrHeter,RH1PO4TransfLitrHeter=substrate-limited H2PO4,HPO4 mineraln-immobiln
!
        CGROMC             = DOMuptk4GrothAutor(ielmc,NGL)-RespGrossAutor(NGL)-RNOxReduxRespAutorLim(NGL)-Resp4NFixAutor(NGL)
        RCO2ProdAutor(NGL) = RCO2ProdAutor(NGL)+Resp4NFixAutor(NGL)
        MID3               = micpar%get_micb_id(3,NGL)
        DO M=1,2
          DO NE=1,NumPlantChemElms
            mBiomeAutor(NE,MID3)=mBiomeAutor(NE,MID3)-NonstX2stBiomAutor(NE,M,NGL)+RkillRecycOMAutor(NE,M,NGL)
          ENDDO
          !C is respired as CO2 while N and P are recycled.
          mBiomeAutor(ielmn,MID3) = mBiomeAutor(ielmn,MID3)+RMaintDefcitRecycOMAutor(ielmn,M,NGL)
          mBiomeAutor(ielmp,MID3) = mBiomeAutor(ielmp,MID3)+RMaintDefcitRecycOMAutor(ielmp,M,NGL)
          RCO2ProdAutor(NGL)      = RCO2ProdAutor(NGL)+RMaintDefcitRecycOMAutor(ielmc,M,NGL)
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
  end associate
  end subroutine AutotrophAnabolicUpdate


end module MicAutoCPLXMod
