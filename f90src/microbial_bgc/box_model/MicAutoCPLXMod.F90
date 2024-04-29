module MicAutoCPLXMod
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : safe_adb,AZMAX1
  use MicForcTypeMod, only : micforctype
  use MicFluxTypeMod, only: micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use NitroDiagTypes
  use ElmIDMod
  use TracerIDMod
  use EcoSiMParDataMod, only : micpar
  use EcoSIMSolverPar
  use EcoSimConst
  use NitroPars
  use MicrobMathFuncMod
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ActiveMicrobesff
  public :: MicrobialAnabolicUpdateff
  contains

!------------------------------------------------------------------------------------------
  subroutine ActiveMicrobesff(NGL,N,VOLWZ,XCO2,TSensGrowth,WatStressMicb,TOMCNK,OXKX,TOMA,TOMN,RH2GZ,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: VOLWZ
  real(r8), intent(in):: OXKX,tomcnk(2),WatStressMicb,TSensGrowth,XCO2
  real(r8), intent(in) :: TOMA,TOMN
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  real(r8), intent(out):: RH2GZ
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  integer  :: M,K
  real(r8) :: COMC
  real(r8) :: ECHZ
  real(r8) :: ORGCL
  real(r8) :: FOXYX
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX
  real(r8) :: FOQC,FOQA
  real(r8) :: FGOCP,FGOAP
  real(r8) :: RGOCP
  real(r8) :: RGOMP
  real(r8) :: RVOXP
  real(r8) :: RVOXPA
  real(r8) :: RVOXPB
  real(r8) :: RGOMT
  real(r8) :: RXOMT
  real(r8) :: RMOMT
  real(r8) :: SPOMK(2)
  real(r8) :: RMOMK(2)
! begin_execution
  associate(                      &
    FracOMActAutor => nmics%FracOMActAutor, &
    FOMNff => nmics%FOMNff, &
    FOMKff => nmics%FOMKff, &
    OMActAutor => nmics%OMActAutor,         &
    RO2UptkAutor  => nmicf%RO2UptkAutor, &
    RGOMOff  => nmicf%RGOMOff, &
    RO2Dmnd4RespAutor => nmicf%RO2Dmnd4RespAutor , &
    RO2DmndAutor=> nmicf%RO2DmndAutor  , &
    RO2Uptk4RespAutor => nmicf%RO2Uptk4RespAutor , &
    RNO3UptkAutor => nmicf%RNO3UptkAutor , &
    RNO2ReduxAutorSoil => nmicf%RNO2ReduxAutorSoil , &
    RNO2ReduxAutorBand => nmicf%RNO2ReduxAutorBand , &
    RGOMDff  => nmicf%RGOMDff, &
    RH2GXff  => nmicf%RH2GXff, &
    RGOMYff  => nmicf%RGOMYff, &
    RCO2ProdAutor  => nmicf%RCO2ProdAutor, &
    RCH4ProdAutor  => nmicf%RCH4ProdAutor, &
    TOMK  => ncplxs%TOMK     , &
    jcplx  => micpar%jcplx   , &
    AmmoniaOxidBacter => micpar%AmmoniaOxidBacter, &
    NitriteOxidBacter  => micpar%NitriteOxidBacter, &
    H2GenoMethanogArchea  => micpar%H2GenoMethanogArchea, &
    AerobicMethanotrofBacter  => micpar%AerobicMethanotrofBacter, &
    ZEROS  => micfor%ZEROS  , &
    SoilMicPMassLayer  => micfor%SoilMicPMassLayer    , &
    litrm => micfor%litrm  , &
    VLSoilPoreMicP_vr  => micfor%VLSoilPoreMicP_vr, &
    ORGC   => micfor%ORGC  ,     &
    RO2DmndAutort => micflx%RO2DmndAutort  &
  )
! FracOMActHeter,FOMN=fraction of total active biomass C,N in each N and K

  IF(TOMA.GT.ZEROS)THEN
    FracOMActAutor(NGL)=OMActAutor(NGL)/TOMA
  ELSE
    FracOMActAutor(NGL)=1.0_r8
  ENDIF
  IF(TOMN.GT.ZEROS)THEN
    FOMNff(NGL)=OMActAutor(NGL)/TOMN
  ELSE
    FOMNff(NGL)=1.0_r8
  ENDIF

  K=micpar%jcplx
  IF(TOMK(K).GT.ZEROS)THEN
    FOMKff(NGL)=OMActAutor(NGL)/TOMK(K)
  ELSE
    FOMKff(NGL)=1.0_r8
  ENDIF
!
  !     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
  !
  !     COMC=microbial C concentration relative to substrate
  !     SPOMK=effect of microbial C concentration on microbial decay
  !     RMOMK=effect of microbial C concentration on maintenance respn
  !
  ORGCL=AMIN1(1.0E+05_r8*SoilMicPMassLayer,ORGC)
  IF(ORGCL.GT.ZEROS)THEN
    DO M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
    ENDDO
  ELSE
    DO  M=1,2
      SPOMK(M)=1.0
      RMOMK(M)=1.0
    ENDDO
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
  if (N.eq.AmmoniaOxidBacter)then
!   NH3 OXIDIZERS
    call NH3OxidizerCatabolism(NGL,N,XCO2,VOLWZ,TSensGrowth,ECHZ,RGOMP,RVOXP,&
      RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.NitriteOxidBacter)then
!     NO2 OXIDIZERS
    call NO2OxidizerCatabolism(NGL,N,XCO2,ECHZ,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.H2GenoMethanogArchea)then
!     H2TROPHIC METHANOGENS
    call H2MethanogensCatabolism(NGL,N,ECHZ,RGOMP,XCO2,micfor,micstt,&
      naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.AerobicMethanotrofBacter)then
!     METHANOTROPHS
    call MethanotrophCatabolism(NGL,N,ECHZ,RGOMP,&
      RVOXP,RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  else
    RGOMP=0.0_r8
    RO2Dmnd4RespAutor(NGL)=0.0_r8
    RO2DmndAutor(NGL)=0.0_r8
    RO2DmndAutort(NGL)=0.0_r8
  ENDif
!
!  write(*,*)'O2 UPTAKE BY AEROBES'
!
! RO2UptkHeter, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
! RUPMX=O2-unlimited rate of O2 uptake
! FOXYX=fraction of O2 uptake by N,K relative to total
! dts_gas=1/(NPH*NPT)
! ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
! OLSGL=aqueous O2 diffusivity
! OXYG,OXYS=gaseous, aqueous O2 amounts
! Rain2LitRSurf_col,Irrig2LitRSurf=surface water flux from precipitation, irrigation
! O2_rain_conc,O2_irrig_conc=O2 concentration in Rain2LitRSurf_col,Irrig2LitRSurf
!
  RO2UptkAutor(NGL)=0.0_r8

  if (N.eq.AmmoniaOxidBacter .or. N.eq.NitriteOxidBacter .or. N.eq.AerobicMethanotrofBacter)then
!   write(*,*)'AerobLeafO2Solubility_pftUptake'
    call AerobicAutorO2Uptake(NGL,N,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      micfor,micstt,nmicf,nmics,micflx)
  elseif (N.eq.H2GenoMethanogArchea)then
    RGOMOff(NGL)=RGOMP
    RCO2ProdAutor(NGL)=0.0_r8
    RCH4ProdAutor(NGL)=RGOMOff(NGL)
    RO2Uptk4RespAutor(NGL)=RO2Dmnd4RespAutor(NGL)
    RH2GXff(NGL)=0.0_r8
    RH2GZ=0.667_r8*RGOMOff(NGL)
  ENDIF
!
!     AUTOTROPHIC DENITRIFICATION
!
  IF(N.EQ.AmmoniaOxidBacter.AND.RO2Dmnd4RespAutor(NGL).GT.0.0_r8.AND.(.not.litrm.OR.VLSoilPoreMicP_vr.GT.ZEROS))THEN
    call AutotrophDenitrificCatabolism(NGL,N,XCO2,VOLWZ,micfor,micstt,&
      naqfdiag,nmicf,nmics,micflx)
  ELSE
    RNO3UptkAutor(NGL)=0.0_r8
    RNO2ReduxAutorSoil(NGL)=0.0_r8
    RNO2ReduxAutorBand(NGL)=0.0_r8
    RGOMYff(NGL)=0.0_r8
    RGOMDff(NGL)=0.0_r8
  ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
  call BiomassMineralizationff(NGL,N,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,nmics)
!
  call GatherMicrobialRespirationff(NGL,N,RMOMK,micfor,micstt,RGOMT,&
    RXOMT,RMOMT,nmicf,nmics)
!
  call GetMicrobialAnabolismFluxff(NGL,N,ECHZ,FGOCP,&
    FGOAP,RGOMT,RXOMT,RMOMT,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  end associate
  end subroutine ActiveMicrobesff


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
  type(NitroAQMFluxDiagType),INTENT(INOUT)::  naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                  &
    FracOMActAutor    => nmics%FracOMActAutor, &
    FOMKff    => nmics%FOMKff, &
    FNH4XRff  =>nmicf%FNH4XRff,   &
    FNO3XRff  =>nmicf%FNO3XRff,   &
    FP14XRff  => nmicf%FP14XRff,  &
    FPO4XRff  => nmicf%FPO4XRff,  &
    ROXYY     => micfor%ROXYY, &
    RNH4Y     => micfor%RNH4Y, &
    RNHBY     => micfor%RNHBY, &
    RNO3Y     => micfor%RNO3Y, &
    RN3BY     => micfor%RN3BY, &
    RPO4Y     => micfor%RPO4Y, &
    RPOBY     => micfor%RPOBY, &
    RP14Y     => micfor%RP14Y, &
    RP1BY     => micfor%RP1BY, &
    RNH4YU    => micfor%RNH4YU,   &
    RNO3YU    => micfor%RNO3YU,   &
    RP14YU    => micfor%RP14YU,   &
    RPO4YU    => micfor%RPO4YU,   &
    ROQCY     => micfor%ROQCY, &
    ROQAY     => micfor%ROQAY, &
    SoilMicPMassLayer0     => micfor%SoilMicPMassLayer0,   &
    Lsurf     => micfor%Lsurf,   &
    litrm     => micfor%litrm,   &
    ZEROS     => micfor%ZEROS,   &
    VLNO3     => micfor%VLNO3,   &
    VLNOB     => micfor%VLNOB,   &
    VLPO4     => micfor%VLPO4,   &
    VLPOB     => micfor%VLPOB,   &
    VLNHB     => micfor%VLNHB,   &
    VLNH4     => micfor%VLNH4,   &
    RH1PO4UptkLitrAutor  => micflx%RH1PO4UptkLitrAutor, &
    RH2PO4UptkLitrAutor  => micflx%RH2PO4UptkLitrAutor, &
    RNO3UptkLitrAutor  => micflx%RNO3UptkLitrAutor, &
    RNH4UptkLitrAutor  => micflx%RNH4UptkLitrAutor, &
    RH1PO4UptkBandAutor   => micflx%RH1PO4UptkBandAutor, &
    RH1PO4UptkSoilAutor   => micflx%RH1PO4UptkSoilAutor, &
    RH2PO4UptkBandAutor   => micflx%RH2PO4UptkBandAutor, &
    RH2PO4UptkSoilAutor   => micflx%RH2PO4UptkSoilAutor, &
    RNO3UptkBandAutor   => micflx%RNO3UptkBandAutor, &
    RNO3UptkSoilAutor   => micflx%RNO3UptkSoilAutor, &
    RNH4UptkBandAutor   => micflx%RNH4UptkBandAutor, &
    RNH4UptkSoilAutor   => micflx%RNH4UptkSoilAutor, &
    RO2DmndAutort   => micflx%RO2DmndAutort  &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(ROXYY.GT.ZEROS)THEN
    FOXYX=AMAX1(FMN,RO2DmndAutort(NGL)/ROXYY)
  ELSE
    FOXYX=AMAX1(FMN,FracOMActAutor(NGL))
  ENDIF
  IF(RNH4Y.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4UptkSoilAutor(NGL)/RNH4Y)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNH4)
  ENDIF
  IF(RNHBY.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4UptkBandAutor(NGL)/RNHBY)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNHB)
  ENDIF
  IF(RNO3Y.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3UptkSoilAutor(NGL)/RNO3Y)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNO3)
  ENDIF
  IF(RN3BY.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3UptkBandAutor(NGL)/RN3BY)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNOB)
  ENDIF
  IF(RPO4Y.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4UptkSoilAutor(NGL)/RPO4Y)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RPOBY.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4UptkBandAutor(NGL)/RPOBY)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF
  IF(RP14Y.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4UptkSoilAutor(NGL)/RP14Y)
  ELSE
    FP14X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RP1BY.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4UptkBandAutor(NGL)/RP1BY)
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
    IF(RNH4YU.GT.ZEROS)THEN
      FNH4XRff(NGL)=AMAX1(FMN,RNH4UptkLitrAutor(NGL)/RNH4YU)
    ELSE
      FNH4XRff(NGL)=AMAX1(FMN,FOMKff(NGL))
    ENDIF
    IF(RNO3YU.GT.ZEROS)THEN
      FNO3XRff(NGL)=AMAX1(FMN,RNO3UptkLitrAutor(NGL)/RNO3YU)
    ELSE
      FNO3XRff(NGL)=AMAX1(FMN,FOMKff(NGL))
    ENDIF
    IF(RPO4YU.GT.ZEROS)THEN
      FPO4XRff(NGL)=AMAX1(FMN,RH2PO4UptkLitrAutor(NGL)/RPO4YU)
    ELSE
      FPO4XRff(NGL)=AMAX1(FMN,FOMKff(NGL))
    ENDIF
    IF(RP14YU.GT.ZEROS)THEN
      FP14XRff(NGL)=AMAX1(FMN,RH1PO4UptkLitrAutor(NGL)/RP14YU)
    ELSE
      FP14XRff(NGL)=AMAX1(FMN,FOMKff(NGL))
    ENDIF
  ENDIF
  IF(Lsurf.AND.SoilMicPMassLayer0.GT.ZEROS)THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4XRff(NGL)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3XRff(NGL)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4XRff(NGL)
    naqfdiag%TFP14X=naqfdiag%TFP14X+FP14XRff(NGL)
  ENDIF
  end associate
  end subroutine SubstrateCompetitionFactorsff
!------------------------------------------------------------------------------------------

  subroutine GetMicrobialAnabolismFluxff(NGL,N,ECHZ,FGOCP,FGOAP,&
    RGOMT,RXOMT,RMOMT,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: ECHZ
  real(r8), intent(in) :: FGOCP,FGOAP,RGOMT,RXOMT,RMOMT
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
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
  associate(                      &
    CNOMActAutor  => nmics%CNOMActAutor, &
    CPOMActAutor  => nmics%CPOMActAutor, &
    TFNGff   => nmics%TFNGff , &
    CGOMEautor  => nmicf%CGOMEautor, &
    CGOSEautor  => nmicf%CGOSEautor, &
    RGOMOff  => nmicf%RGOMOff, &
    RGOMDff  => nmicf%RGOMDff, &
    RMaintCompAutor  => nmicf%RMaintCompAutor, &
    RDOMEautor => nmicf%RDOMEautor, &
    RHOMEautor  => nmicf%RHOMEautor, &
    RCOMEautor  => nmicf%RCOMEautor, &
    RDMMEautor  => nmicf%RDMMEautor, &
    RHMMEautor  => nmicf%RHMMEautor, &
    RCMMEautor  => nmicf%RCMMEautor, &
    RXOMEautor  => nmicf%RXOMEautor, &
    R3OMEautor  => nmicf%R3OMEautor, &
    RXMMEautor  => nmicf%RXMMEautor, &
    R3MMEautor  => nmicf%R3MMEautor, &
    RGN2Fff  => nmicf%RGN2Fff, &
    TCGOMEheter   => ncplxf%TCGOMEheter, &
    CNQ      => ncplxs%CNQ, &
    CPQ      => ncplxs%CPQ, &
    rNCOMCff  => micpar%rNCOMCff,   &
    rPCOMCff  => micpar%rPCOMCff,   &
    FL       => micpar%FL       , &
    ZEROS    => micfor%ZEROS    , &
    ZERO     => micfor%ZERO     , &
    OMEauto    => micstt%OMEauto    , &
    EHUM     => micstt%EHUM       &
  )

!     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
!     FROM O2, NOX AND C REDUCTION
!
!     CGOMX=DOC+acetate uptake from aerobic growth respiration
!     CGOMD=DOC+acetate uptake from denitrifier growth respiration
!     RMOMT=maintenance respiration
!     RGOMO=total respiration
!     RGOMD=respiration for denitrifcation
!     RGN2F=respiration for N2 fixation
!     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
!     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake(heterotrophs
!     CGOMC=total CO2,CH4 uptake (autotrophs)
!     CGOMN,CGOMP=DON, DOP uptake
!     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
!     OQN,OPQ=DON,DOP
!     FOMK=faction of OMA in total OMA
!     CNQ,CPQ=DON/DOC, DOP/DOC
!     FCN,FCP=limitation from N,P
!

  CGOMEautor(idom_beg:idom_end,NGL)=0._r8
  CGOMX=AMIN1(RMOMT,RGOMOff(NGL))+RGN2Fff(NGL)+(RGOMT-RGN2Fff(NGL))/ECHZ
  CGOMD=RGOMDff(NGL)/ENOX
  CGOMEautor(ielmc,NGL)=CGOMX+CGOMD
  K=micpar%jcplx
  DO idom=idom_beg,idom_end
    TCGOMEheter(idom,K)=TCGOMEheter(idom,K)+CGOMEautor(idom,NGL)
  ENDDO
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
  IF(OMEauto(ielmc,MID3).GT.ZEROS.AND.OMEauto(ielmc,MID1).GT.ZEROS)THEN
    CCC=AZMAX1(AMIN1(1.0_r8 &
      ,OMEauto(ielmn,MID3)/(OMEauto(ielmn,MID3)+OMEauto(ielmc,MID3)*rNCOMCff(3,NGL)) &
      ,OMEauto(ielmp,MID3)/(OMEauto(ielmp,MID3)+OMEauto(ielmc,MID3)*rPCOMCff(3,NGL))))
    CXC=OMEauto(ielmc,MID3)/OMEauto(ielmc,MID1)
    C3C=1.0_r8/(1.0_r8+CXC/CKC)
    CNC=AZMAX1(AMIN1(1.0_r8 &
      ,OMEauto(ielmc,MID3)/(OMEauto(ielmc,MID3)+OMEauto(ielmn,MID3)/rNCOMCff(3,NGL))))
    CPC=AZMAX1(AMIN1(1.0_r8 &
      ,OMEauto(ielmc,MID3)/(OMEauto(ielmc,MID3)+OMEauto(ielmp,MID3)/rPCOMCff(3,NGL))))
    RCCC=RCCZ+AMAX1(CCC,C3C)*RCCY
    RCCN=CNC*RCCX
    RCCP=CPC*RCCQ
  ELSE
    RCCC=RCCZ
    RCCN=0.0_r8
    RCCP=0.0_r8
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
  MID3=micpar%get_micb_id(3,NGL)
  CGOMZ=TFNGff(NGL)*OMGR*AZMAX1(OMEauto(ielmc,MID3))
  DO  M=1,2
    CGOSEautor(ielmc,M,NGL)=FL(M)*CGOMZ
    IF(OMEauto(ielmc,MID3).GT.ZEROS)THEN
      CGOSEautor(ielmn,M,NGL)=AMIN1(FL(M)*AZMAX1(OMEauto(ielmn,MID3)) &
        ,CGOSEautor(ielmc,M,NGL)*OMEauto(ielmn,MID3)/OMEauto(ielmc,MID3))
      CGOSEautor(ielmp,M,NGL)=AMIN1(FL(M)*AZMAX1(OMEauto(ielmp,MID3)) &
        ,CGOSEautor(ielmc,M,NGL)*OMEauto(ielmp,MID3)/OMEauto(ielmc,MID3))
    ELSE
      CGOSEautor(ielmn,M,NGL)=0.0_r8
      CGOSEautor(ielmp,M,NGL)=0.0_r8
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
    MID=micpar%get_micb_id(M,NGL)
    SPOMX=SQRT(TFNGff(NGL))*SPOMC(M)*SPOMK(M)
    RXOMEautor(ielmc,M,NGL)=AZMAX1(OMEauto(ielmc,MID)*SPOMX)
    RXOMEautor(ielmn,M,NGL)=AZMAX1(OMEauto(ielmn,MID)*SPOMX)
    RXOMEautor(ielmp,M,NGL)=AZMAX1(OMEauto(ielmp,MID)*SPOMX)

    RDOMEautor(ielmc,M,NGL)=RXOMEautor(ielmc,M,NGL)*(1.0_r8-RCCC)
    RDOMEautor(ielmn,M,NGL)=RXOMEautor(ielmn,M,NGL)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDOMEautor(ielmp,M,NGL)=RXOMEautor(ielmp,M,NGL)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
    DO NE=1,NumPlantChemElms
      R3OMEautor(NE,M,NGL)=RXOMEautor(NE,M,NGL)-RDOMEautor(NE,M,NGL)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P LitrFall to residue
!
      RHOMEautor(NE,M,NGL)=AZMAX1(RDOMEautor(NE,M,NGL)*EHUM)
  !
  !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
  !
      RCOMEautor(NE,M,NGL)=RDOMEautor(NE,M,NGL)-RHOMEautor(NE,M,NGL)
    ENDDO
  ENDDO
!
!     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
!     EXCEEDS UPTAKE
!
!     OMC,OMN,OMP=microbial C,N,P
!     RMOMT=total maintenance respiration
!     RXOMT=senescence respiration
!     RCCC=C recycling fraction
!     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
!     RMOMC=maintenance respiration
!     CNOMA,CPOMA=N:C,P:C ratios of active biomass
!     RDMMC,RDMMN,RDMMP=microbial C,N,P LitrFall from senescence
!     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!
  IF(RXOMT.GT.ZEROS.AND.RMOMT.GT.ZEROS.AND.RCCC.GT.ZERO)THEN
    FRM=RXOMT/RMOMT
    DO  M=1,2
      RXMMEautor(ielmc,M,NGL)=AMIN1(OMEauto(ielmc,MID),AZMAX1(FRM*RMaintCompAutor(M,NGL)/RCCC))
      RXMMEautor(ielmn,M,NGL)=AMIN1(OMEauto(ielmn,MID),AZMAX1(RXMMEautor(ielmc,M,NGL)*CNOMActAutor(NGL)))
      RXMMEautor(ielmp,M,NGL)=AMIN1(OMEauto(ielmp,MID),AZMAX1(RXMMEautor(ielmc,M,NGL)*CPOMActAutor(NGL)))

      RDMMEautor(ielmc,M,NGL)=RXMMEautor(ielmc,M,NGL)*(1.0_r8-RCCC)
      RDMMEautor(ielmn,M,NGL)=RXMMEautor(ielmn,M,NGL)*(1.0_r8-RCCN)*(1.0_r8-RCCC)
      RDMMEautor(ielmp,M,NGL)=RXMMEautor(ielmp,M,NGL)*(1.0_r8-RCCP)*(1.0_r8-RCCC)
      DO NE=1,NumPlantChemElms
        R3MMEautor(NE,M,NGL)=RXMMEautor(NE,M,NGL)-RDMMEautor(NE,M,NGL)
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
        RHMMEautor(NE,M,NGL)=AZMAX1(RDMMEautor(NE,M,NGL)*EHUM)
        RCMMEautor(NE,M,NGL)=RDMMEautor(NE,M,NGL)-RHMMEautor(NE,M,NGL)
      ENDDO
    ENDDO
  ELSE
    DO  M=1,2
      DO NE=1,NumPlantChemElms
        RXMMEautor(NE,M,NGL)=0.0_r8
        RDMMEautor(NE,M,NGL)=0.0_r8
        R3MMEautor(NE,M,NGL)=0.0_r8
        RHMMEautor(NE,M,NGL)=0.0_r8
        RCMMEautor(NE,M,NGL)=0.0_r8
      ENDDO
    ENDDO
  ENDIF
  end associate
  end subroutine GetMicrobialAnabolismFluxff
!------------------------------------------------------------------------------------------

  subroutine AerobicAutorO2Uptake(NGL,N,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
    micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL   !guild id
  integer, intent(in) :: N     !functional group id
  real(r8), intent(in) :: OXKX,FOXYX,RGOMP,RVOXP
  real(r8), intent(in) :: RVOXPA
  real(r8), intent(in) :: RVOXPB
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx

  integer  :: M,MX
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,OLSGL1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMX,ROXYFX
  real(r8) :: ROXYLX
  real(r8) :: RRADO,RMPOX,ROXDFQ
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X
  real(r8) :: VOLWPM

  ! begin_execution
  associate(                  &
    fLimO2Autor    => nmics%fLimO2Autor,      &
    OMActAutor    => nmics%OMActAutor,      &
    RO2UptkAutor  => nmicf%RO2UptkAutor, &
    RGOMOff  => nmicf%RGOMOff, &
    RO2Dmnd4RespAutor => nmicf%RO2Dmnd4RespAutor,     &
    RO2DmndAutor=> nmicf%RO2DmndAutor ,     &
    RO2Uptk4RespAutor => nmicf%RO2Uptk4RespAutor,     &
    RH2GXff  => nmicf%RH2GXff, &
    RCO2ProdAutor  => nmicf%RCO2ProdAutor, &
    RCH4ProdAutor  => nmicf%RCH4ProdAutor, &
    RSOxidSoilAutor  => nmicf%RSOxidSoilAutor, &
    RSOxidBandAutor  => nmicf%RSOxidBandAutor, &
    AmmoniaOxidBacter => micpar%AmmoniaOxidBacter, &
    NitriteOxidBacter => micpar%NitriteOxidBacter, &
    ROXYF  => micfor%ROXYF,  &
    COXYE  => micfor%COXYE  , &
    O2_rain_conc   => micfor%O2_rain_conc  , &
    O2_irrig_conc   => micfor%O2_irrig_conc  , &
    Irrig2LitRSurf  => micfor%Irrig2LitRSurf  , &
    Rain2LitRSurf_col  => micfor%Rain2LitRSurf_col  , &
    litrm  => micfor%litrm  , &
    OLSGL => micfor%OLSGL , &
    ROXYL => micfor%ROXYL  , &
    VLSoilPoreMicP_vr  => micfor%VLSoilPoreMicP_vr   , &
    VLSoilMicP  => micfor%VLSoilMicP  , &
    VLsoiAirPM  => micfor%VLsoiAirPM , &
    VLWatMicPM  => micfor%VLWatMicPM , &
    FILM  => micfor%FILM , &
    THETPM => micfor%THETPM ,&
    TortMicPM => micfor%TortMicPM , &
    ZERO => micfor%ZERO , &
    ZEROS => micfor%ZEROS , &
    DiffusivitySolutEff => micfor%DiffusivitySolutEff , &
    SOXYL => micstt%SOXYL  , &
    OXYG => micstt%OXYG , &
    OXYS => micstt%OXYS , &
    COXYG => micstt%COXYG,   &
    ROXSK=> micflx%ROXSK, &
    RNH3OxidAutor => micflx%RNH3OxidAutor,  &
    RNH3OxidAutorBand => micflx%RNH3OxidAutorBand, &
    RNO2OxidAutor => micflx%RNO2OxidAutor, &
    RNO2OxidAutorBand => micflx%RNO2OxidAutorBand  &
  )

  IF(RO2DmndAutor(NGL).GT.ZEROS.AND.FOXYX.GT.ZERO)THEN
    IF(.not.litrm.OR.VLSoilPoreMicP_vr.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX=RO2DmndAutor(NGL)*dts_gas    !
      ROXYFX=ROXYF*dts_gas*FOXYX    !O2 demand
      OLSGL1=OLSGL*dts_gas
      IF(.not.litrm)THEN
        OXYG1=OXYG*FOXYX
        ROXYLX=ROXYL*dts_gas*FOXYX
      ELSE
        OXYG1=COXYG*VLsoiAirPM(1)*FOXYX
        ROXYLX=(ROXYL+Rain2LitRSurf_col*O2_rain_conc &
          +Irrig2LitRSurf*O2_irrig_conc)*dts_gas*FOXYX
      ENDIF
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
        !     VLWatMicPM,VLsoiAirPM,VLSoilPoreMicP_vr=water, air and total volumes
        !     ORAD=microbial radius,FILM=water film thickness
        !     DIFOX=aqueous O2 diffusion, TortMicPM=tortuosity
        !     BIOS=microbial number, OMA=active biomass
        !     SOXYL=O2 solubility, OXKX=Km for O2 uptake
        !     OXYS,COXYS=aqueous O2 amount, concentration
        !     OXYG,COXYG=gaseous O2 amount, concentration
        !     RMPOX,ROXSK=O2 uptake
        !
        THETW1=AZMAX1(safe_adb(VLWatMicPM(M),VLSoilMicP))
        RRADO=ORAD*(FILM(M)+ORAD)/FILM(M)
        DIFOX=TortMicPM(M)*OLSGL1*12.57_r8*BIOS*OMActAutor(NGL)*RRADO
        VOLWOX=VLWatMicPM(M)*SOXYL
        VOLPOX=VLsoiAirPM(M)
        VOLWPM=VOLWOX+VOLPOX

!oxygen uptake in the layer
        DO  MX=1,NPT
          OXYG1=OXYG1+ROXYFX
          OXYS1=OXYS1+ROXYLX
          COXYS1=AMIN1(COXYE*SOXYL,AZMAX1(safe_adb(OXYS1,(VLWatMicPM(M)*FOXYX))))

          !solve for uptake flux
          IF(OXYS1<=ZEROS)THEN
            RMPOX=0.0_r8
          else
            RMPOX=TranspBasedsubstrateUptake(COXYS1,DIFOX, OXKX, RUPMX, ZEROS)
          ENDIF

          !apply the uptake flux
          OXYS1=OXYS1-RMPOX
          !apply volatilization-dissolution
          IF(THETPM(M).GT.THETX.AND.VOLPOX.GT.ZEROS)THEN
            ROXDFQ=DiffusivitySolutEff(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
          ELSE
            ROXDFQ=0.0_r8
          ENDIF
          OXYG1=OXYG1-ROXDFQ
          OXYS1=OXYS1+ROXDFQ
          !accumulate uptake flux
          RO2UptkAutor(NGL)=RO2UptkAutor(NGL)+RMPOX
          ROXSK(M)=ROXSK(M)+RMPOX
        ENDDO

      ENDDO
      !write(*,*)'420'
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
      !
      !     WFN=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RVMX2,RVMB2=NH3,NO2 oxidation in non-band, band
      !
      fLimO2Autor(NGL)=AMIN1(1.0,AZMAX1(RO2UptkAutor(NGL)/RO2DmndAutor(NGL)))
      IF(N.EQ.AmmoniaOxidBacter)THEN
        RNH3OxidAutor(NGL)=RNH3OxidAutor(NGL)*fLimO2Autor(NGL)
        RNH3OxidAutorBand(NGL)=RNH3OxidAutorBand(NGL)*fLimO2Autor(NGL)
      ELSEIF(N.EQ.NitriteOxidBacter)THEN
        RNO2OxidAutor(NGL)=RNO2OxidAutor(NGL)*fLimO2Autor(NGL)
        RNO2OxidAutorBand(NGL)=RNO2OxidAutorBand(NGL)*fLimO2Autor(NGL)
      ENDIF
    ELSE
      RO2UptkAutor(NGL)=RO2DmndAutor(NGL)
      fLimO2Autor(NGL)=1.0_r8
    ENDIF
  ELSE
    RO2UptkAutor(NGL)=0.0_r8
    fLimO2Autor(NGL)=1.0_r8
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RGOMO,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2X,RAcettProdHeter,RCH4ProdHeter,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
  !     RO2Uptk4RespHeter=O2-limited O2 uptake
  !     RSOxidSoilAutor,RSOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !
  RGOMOff(NGL)=RGOMP*fLimO2Autor(NGL)
  RCO2ProdAutor(NGL)=RGOMOff(NGL)
  RCH4ProdAutor(NGL)=0.0_r8
  RO2Uptk4RespAutor(NGL)=RO2Dmnd4RespAutor(NGL)*fLimO2Autor(NGL)
  RH2GXff(NGL)=0.0_r8
  RSOxidSoilAutor(NGL)=RVOXPA*fLimO2Autor(NGL)
  RSOxidBandAutor(NGL)=RVOXPB*fLimO2Autor(NGL)
  end associate
  end subroutine AerobicAutorO2Uptake

!------------------------------------------------------------------------------------------

  subroutine AutotrophDenitrificCatabolism(NGL,N,XCO2,VOLWZ,micfor,micstt,&
    naqfdiag,nmicf,nmics, micflx)

  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: VOLWZ,XCO2
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
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
  associate(                              &
    FOMNff               => nmics%FOMNff, &
    RO2Dmnd4RespAutor  => nmicf%RO2Dmnd4RespAutor, &
    RO2Uptk4RespAutor => nmicf%RO2Uptk4RespAutor , &
    RNO3UptkAutor => nmicf%RNO3UptkAutor , &
    RNO2ReduxAutorSoil => nmicf%RNO2ReduxAutorSoil , &
    RNO2ReduxAutorBand => nmicf%RNO2ReduxAutorBand , &
    RGOMDff  => nmicf%RGOMDff, &
    RGOMYff  => nmicf%RGOMYff, &
    RSOxidSoilAutor    => nmicf%RSOxidSoilAutor  , &
    RSOxidBandAutor    => nmicf%RSOxidBandAutor  , &
    RTotNH3OxidSoilAutor => nmicf%RTotNH3OxidSoilAutor,   &
    RTotNH3OxidBandAutor => nmicf%RTotNH3OxidBandAutor,   &
    RNO2EcoUptkSoilPrev  => micfor%RNO2EcoUptkSoilPrev   , &
    VLNO3  => micfor%VLNO3   , &
    VLNOB  => micfor%VLNOB   , &
    RNO2EcoUptkBandPrev  =>  micfor%RNO2EcoUptkBandPrev  , &
    ZEROS  => micfor%ZEROS   , &
    ZEROS2  => micfor%ZEROS2   , &
    CNO2B  => micstt%CNO2B   , &
    CNO2S  => micstt%CNO2S   , &
    ZNO2B  => micstt%ZNO2B   , &
    ZNO2S  => micstt%ZNO2S   , &
    RNO2OxidAutor => micflx%RNO2OxidAutor, &
    RNO2OxidAutorBand  => micflx%RNO2OxidAutorBand    &
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
    FNO2=AMAX1(FMN,FOMNff(NGL)*VLNO3)
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2OxidAutorBand(NGL)/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=AMAX1(FMN,FOMNff(NGL)*VLNOB)
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
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NO2 reduction
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
  RGOMYff(NGL)=0.0_r8
  RGOMDff(NGL)=RDNOT*ECNO*ENOX
  RNO3UptkAutor(NGL)=0.0_r8
  RNO2OxidAutor(NGL)=VMXD4S
  RNO2OxidAutorBand(NGL)=VMXD4B
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
  type(NitroAQMFluxDiagType),INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
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
  associate(             &
    TFNGff => nmics%TFNGff,  &
    FCNPff => nmics%FCNPff,  &
    FracOMActAutor => nmics%FracOMActAutor,  &
    OMActAutor  => nmics%OMActAutor ,  &
    RO2Dmnd4RespAutor=> nmicf%RO2Dmnd4RespAutor, &
    RO2DmndAutor=> nmicf%RO2DmndAutor, &
    VLNH4  => micfor%VLNH4 , &
    VLNHB  => micfor%VLNHB , &
    ZEROS  => micfor%ZEROS , &
    ZEROS2 => micfor%ZEROS2, &
    RNH4Y  => micfor%RNH4Y  , &
    RNHBY  => micfor%RNHBY  , &
    ZNFN0  => micstt%ZNFN0 , &
    ZNFNI  => micstt%ZNFNI , &
    CNH4S  => micstt%CNH4S , &
    CNH4B  => micstt%CNH4B , &
    ZNH4S  => micstt%ZNH4S , &
    ZNH4B  => micstt%ZNH4B , &
    RNH3OxidAutor=> micflx%RNH3OxidAutor, &
    RNH3OxidAutorBand=> micflx%RNH3OxidAutorBand,  &
    RO2DmndAutort => micflx%RO2DmndAutort   &
  )
!
!     FACTOR TO REGULATE COMPETITION FOR NH4 AMONG DIFFERENT
!     MICROBIAL AND ROOT POPULATIONS FNH4
!
!     FNH4,FNB4=frac of total biol demand for NH4 in non-band, band
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  IF(RNH4Y.GT.ZEROS)THEN
    FNH4=AMAX1(FMN,RNH3OxidAutor(NGL)/RNH4Y)
  ELSE
    FNH4=AMAX1(FMN,VLNH4*FracOMActAutor(NGL))
  ENDIF
  IF(RNHBY.GT.ZEROS)THEN
    FNB4=AMAX1(FMN,RNH3OxidAutorBand(NGL)/RNHBY)
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
!     TFNG=temperature+water limitation, FCNP=N,P limitation
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
  VMXX=VMXH*TFNGff(NGL)*FCNPff(NGL)*XCO2*OMActAutor(NGL)
  IF(VOLWZ.GT.ZEROS2)THEN
    VMXA=VMXX/(1.0_r8+VMXX/(VHKI*VOLWZ))
  ELSE
    VMXA=0.0_r8
  ENDIF
  FCN4S=FNH4S*CNH4S/(CNH4S+ZHKM)
  FCN4B=FNHBS*CNH4B/(CNH4B+ZHKM)
  FSBST=FCN4S+FCN4B
  VMX4S=VMXA*FCN4S
  VMX4B=VMXA*FCN4B
  RNNH4=AZMAX1(AMIN1(VMX4S,FNH4*ZNH4S))*ZNFN4S
  RNNHB=AZMAX1(AMIN1(VMX4B,FNB4*ZNH4B))*ZNFN4B
  RVOXP=RNNH4+RNNHB
  RVOXPA=RNNH4
  RVOXPB=RNNHB
  RGOMP=AZMAX1(RVOXP*ECNH*ECHZ)
  RNH3OxidAutor(NGL)=VMX4S
  RNH3OxidAutorBand(NGL)=VMX4B
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
  implicit none
  integer, intent(in) :: NGL,N
  REAL(r8), intent(in) :: XCO2
  real(r8), intent(out) :: ECHZ,RGOMP,RVOXP
  real(r8), intent(out) :: RVOXPA,RVOXPB
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: FNH4S,FNHBS
  REAL(R8) :: fno2,FNB2
  real(r8) :: FCN2S,FCN2B
  real(r8) :: RNNO2,RNNOB
  real(r8) :: VMX2S,VMX2B
  real(r8) :: FSBST
  real(r8) :: VMXA

!     begin_execution
  associate(               &
    TFNGff => nmics%TFNGff, &
    FCNPff => nmics%FCNPff, &
    FOMNff => nmics%FOMNff, &
    OMActAutor  => nmics%OMActAutor , &
    RO2Dmnd4RespAutor => nmicf%RO2Dmnd4RespAutor,  &
    RO2DmndAutor => nmicf%RO2DmndAutor,  &
    VLNH4   => micfor%VLNH4   , &
    VLNHB  => micfor%VLNHB   , &
    VLNO3  => micfor%VLNO3   , &
    VLNOB  => micfor%VLNOB   , &
    ZEROS  => micfor%ZEROS   , &
    RNO2EcoUptkSoilPrev  => micfor%RNO2EcoUptkSoilPrev   , &
    RNO2EcoUptkBandPrev  =>  micfor%RNO2EcoUptkBandPrev  , &
    CNO2S  => micstt%CNO2S   , &
    CNO2B  => micstt%CNO2B   , &
    ZNO2S  => micstt%ZNO2S   , &
    ZNO2B  => micstt%ZNO2B   , &
    RNO2OxidAutor=> micflx%RNO2OxidAutor  , &
    RNO2OxidAutorBand => micflx%RNO2OxidAutorBand , &
    RO2DmndAutort => micflx%RO2DmndAutort   &
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
    FNO2=AMAX1(FMN,FOMNff(NGL)*VLNO3)
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2OxidAutorBand(NGL)/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=AMAX1(FMN,FOMNff(NGL)*VLNOB)
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
!     TFNG=temperature+water limitation, FCNP=N,P limitation
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
!     RVMX2,RVMB2=nitrifier demand for NO2 in non-band, band
!
  ECHZ=EO2X
  VMXA=TFNGff(NGL)*FCNPff(NGL)*XCO2*OMActAutor(NGL)*VMXN
  FCN2S=FNH4S*CNO2S/(CNO2S+ZNKM)
  FCN2B=FNHBS*CNO2B/(CNO2B+ZNKM)
  FSBST=FCN2S+FCN2B
  VMX2S=VMXA*FCN2S
  VMX2B=VMXA*FCN2B
  RNNO2=AZMAX1(AMIN1(VMX2S,FNO2*ZNO2S))
  RNNOB=AZMAX1(AMIN1(VMX2B,FNB2*ZNO2B))
  RVOXP=RNNO2+RNNOB
  RVOXPA=RNNO2
  RVOXPB=RNNOB
  RGOMP=AZMAX1(RVOXP*ECNO*ECHZ)
  RNO2OxidAutor(NGL)=VMX2S
  RNO2OxidAutorBand(NGL)=VMX2B
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
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(out) :: ECHZ,RGOMP
  REAL(R8), INTENT(IN) :: XCO2
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: GH2X,GH2H
  real(r8) :: H2GSX
  real(r8) :: FSBST
  real(r8) :: VMXA

  associate(                  &
    TFNGff => nmics%TFNGff, &
    FCNPff => nmics%FCNPff, &
    OMActAutor  => nmics%OMActAutor , &
    RO2Dmnd4RespAutor => nmicf%RO2Dmnd4RespAutor,     &
    RO2DmndAutor => nmicf%RO2DmndAutor  , &
    TKS     => micfor%TKS       , &
    CH2GS   => micstt%CH2GS    , &
    H2GS    => micstt%H2GS     , &
    RO2DmndAutort  => micflx%RO2DmndAutort     &
  )
!     begin_execution
!
!     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
!
!     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
!     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
!     VMXA=substrate-unlimited H2 oxidation rate
!     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (TRH2G)
!     CH2GS=H2 concentration, H2KM=Km for H2 uptake
!     RGOMP=H2 oxidation, ROXY*=O2 demand
!
!     and energy yield of hydrogenotrophic
!     methanogenesis GH2X at ambient H2 concentration CH2GS

  GH2X=RGAS*1.E-3_r8*TKS*LOG((AMAX1(1.0E-03_r8,CH2GS)/H2KI)**4)
  GH2H=GH2X/12.08_r8
  ECHZ=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCOX+GH2H))/EOMH)))
  VMXA=TFNGff(NGL)*FCNPff(NGL)*XCO2*OMActAutor(NGL)*VMXC
  H2GSX=H2GS+0.111_r8*naqfdiag%TRH2G
  FSBST=CH2GS/(CH2GS+H2KM)
  RGOMP=AZMAX1(AMIN1(1.5*H2GSX,VMXA*FSBST))
  RO2Dmnd4RespAutor(NGL)=0.0_r8
  RO2DmndAutor(NGL)=0.0_r8
  RO2DmndAutort(NGL)=0.0_r8
  naqfdiag%TCH4A=naqfdiag%TCH4A+RGOMP
!
  end associate
  end subroutine H2MethanogensCatabolism
!------------------------------------------------------------------------------------------

  subroutine MethanotrophCatabolism(NGL,N,ECHZ,RGOMP,&
    RVOXP,RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(out) :: ECHZ,RGOMP,RVOXP
  real(r8), intent(out) :: RVOXPA,RVOXPB
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(NitroAQMFluxDiagType), intent(in) :: naqfdiag
  type(NitroMicFLuxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  integer  :: M,MM
  real(r8) :: CH4G1,CH4S1,CCH4S1
  real(r8) :: RCH4L1,RCH4F1,RCH4S1
  real(r8) :: RGOMP1,RCHDF
  real(r8) :: VMXA1,VOLWCH
  real(r8) :: FSBST
  real(r8) :: RVOXP1
  real(r8) :: VMXA
  REAL(R8) :: VOLWPM

  associate(             &
    TFNGff => nmics%TFNGff,  &
    FCNPff => nmics%FCNPff,  &
    OMActAutor  => nmics%OMActAutor ,  &
    RO2Dmnd4RespAutor=> nmicf%RO2Dmnd4RespAutor, &
    RO2DmndAutor=> nmicf%RO2DmndAutor, &
    CCH4E  => micfor%CCH4E  , &
    VLsoiAirPM  => micfor%VLsoiAirPM  , &
    VLWatMicPM  => micfor%VLWatMicPM  , &
    ZEROS2  => micfor%ZEROS2 , &
    ZEROS  => micfor%ZEROS , &
    THETPM  => micfor%THETPM, &
    DiffusivitySolutEff   => micfor%DiffusivitySolutEff   , &
    litrm  => micfor%litrm  , &
    CCH4G  => micstt%CCH4G  , &
    CH4S   => micstt%CH4S   , &
    SCH4L  => micstt%SCH4L  , &
    RCH4L  => micfor%RCH4L , &
    RCH4F  => micfor%RCH4F , &
    RO2DmndAutort  => micflx%RO2DmndAutort &
  )
!     begin_execution
!
!     CH4 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     CH4 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXA=potential oxidation
!     TFNG=temperature+water effect,FCNP=N,P limitation
!     OMA=active biomass,VMX4=specific respiration rate
!     RCH4L=total aqueous CH4 exchange from previous hour
!     RCH4F=total gaseous CH4 exchange from previous hour
!     TCH4H+TCH4A=total CH4 generated from methanogenesis
!     dts_gas=1.0/(NPH*NPT)
!     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
!     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
!     VLsoiAirPM,VLWatMicPM=air,water-filled porosity
!     SCH4L=CH4 aqueous solubility
!     CCK4=Km for CH4 uptake
!     ECHO=efficiency CO2 conversion to biomass
!     RGOMP1=substrate-limited CH4 oxidation
!     RCHDF=gaseous-aqueous CH4 exchange
!     DiffusivitySolutEff=rate constant for gaseous-aqueous exchange
!
  ECHZ=EH4X
  VMXA=TFNGff(NGL)*FCNPff(NGL)*OMActAutor(NGL)*VMX4
  RCH4L1=RCH4L*dts_gas
  RCH4F1=RCH4F*dts_gas
  RCH4S1=(naqfdiag%TCH4H+naqfdiag%TCH4A)*dts_gas

  IF(litrm)THEN
    !surface residue layer
    CH4G1=CCH4E*VLsoiAirPM(1)
  ELSE
    CH4G1=CCH4G*VLsoiAirPM(1)
  ENDIF
  CH4S1=CH4S
  !apparent vmax for uptake
  VMXA1=VMXA*dts_gas
  RVOXP=0.0_r8
  RGOMP=0.0_r8
!
!     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
!     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
! for aerobic methanotrophs, CH4 is oxiized to CO2 for energy, and also 
! to intracellular C for respiration. (Grant 1999), the C yield is approximated
! as the energy required for turning CH4 into organic C, and the energy released
! from turnining CH4 into CO2. 
! the catabolic reaction is
!  CH4 + 2O2 -> CO2 + 2H2O
!    
!  C+ O2 -> CO2,  32/12=2.667
  D320: DO M=1,NPH
    IF(VLWatMicPM(M).GT.ZEROS2)THEN
      VOLWCH=VLWatMicPM(M)*SCH4L
      VOLWPM=VOLWCH+VLsoiAirPM(M)

      !CH4 uptake by aerobic oxidation
      !RCH4F1: net gaseous CH4 flux from transport into the layer
      !RCH4L1: net aqueous CH4 flux from transport into the layer
      !RCH4S1: aquoues CH4 production from methanogenesis in the layer
      D325: DO MM=1,NPT
        CH4G1=CH4G1+RCH4F1
        CH4S1=CH4S1+RCH4L1+RCH4S1
        CCH4S1=AZMAX1(safe_adb(CH4S1,VLWatMicPM(M)))
        
        FSBST=CCH4S1/(CCH4S1+CCK4)
        !RVOXP1 energy from oxidizing CH4 into CO2
        RVOXP1=AMIN1(AZMAX1(CH4S1)/(1.0+ECHO*ECHZ),VMXA1*FSBST)
        !the respiration yield
        RGOMP1=RVOXP1*ECHO*ECHZ
        !
        CH4S1=CH4S1-RVOXP1-RGOMP1

        !dissolution-vaporization
        IF(THETPM(M).GT.THETX)THEN
          RCHDF=DiffusivitySolutEff(M)*(AMAX1(ZEROS,CH4G1)*VOLWCH-CH4S1*VLsoiAirPM(M))/VOLWPM
        ELSE
          RCHDF=0.0_r8
        ENDIF
        CH4G1=CH4G1-RCHDF
        CH4S1=CH4S1+RCHDF
        RVOXP=RVOXP+RVOXP1
        RGOMP=RGOMP+RGOMP1

      ENDDO D325
    ENDIF
  ENDDO D320
  RVOXPA=RVOXP
  RVOXPB=0.0_r8
!
!     O2 DEMAND FROM CH4 OXIDATION
!
!     RO2Dmnd4RespHeter=O2 demand from respiration
!     ROXYP=O2 demand from respiration + CH4 oxidation
!
  RO2Dmnd4RespAutor(NGL)=2.667_r8*RGOMP
  RO2DmndAutor(NGL)=RO2Dmnd4RespAutor(NGL)+5.333_r8*RVOXP
  RO2DmndAutort(NGL)=RO2DmndAutor(NGL)
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
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType),intent(inout) :: nmics
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
  real(r8) :: RINHP
  real(r8) :: RINHX,RINOP,RINOX,RIPOP,RIPOX,RIP1P
  real(r8) :: RIP1X,RINHPR,RINOPR,RIPOPR,RIP1PR
  real(r8) :: ZNH4M,ZNHBM
  real(r8) :: ZNO3M
  real(r8) :: ZNOBM
  integer :: MID3

!     begin_execution
  associate(                 &
   TFNGff    => nmics%TFNGff  ,  &
   OMActAutor     => nmics%OMActAutor  ,  &
   FNH4XRff  => nmicf%FNH4XRff,  &
   FNO3XRff  => nmicf%FNO3XRff,  &
   FPO4XRff  => nmicf%FPO4XRff,  &
   FP14XRff  => nmicf%FP14XRff,  &
   RNH4TransfSoilAutor => nmicf%RNH4TransfSoilAutor   ,  &
   RNO3TransfSoilAutor  => nmicf%RNO3TransfSoilAutor  ,  &
   RH2PO4TransfSoilAutor  => nmicf%RH2PO4TransfSoilAutor  ,  &
   RNH4TransfBandAutor => nmicf%RNH4TransfBandAutor   ,  &
   RNO3TransfBandAutor => nmicf%RNO3TransfBandAutor   ,  &
   RH2PO4TransfBandAutor => nmicf%RH2PO4TransfBandAutor   ,  &
   RNH4TransfLitrAutor  => nmicf%RNH4TransfLitrAutor,  &
   RNO3TransfLitrAutor  => nmicf%RNO3TransfLitrAutor,  &
   RH2PO4TransfLitrAutor  => nmicf%RH2PO4TransfLitrAutor,  &
   RH1PO4TransfSoilAutor  => nmicf%RH1PO4TransfSoilAutor, &
   RH1PO4TransfBandAutor  => nmicf%RH1PO4TransfBandAutor, &
   RH1PO4TransfLitrAutor  => nmicf%RH1PO4TransfLitrAutor,  &
   rNCOMCff  => micpar%rNCOMCff   , &
   rPCOMCff  => micpar%rPCOMCff   , &
   VLNH4    => micfor%VLNH4     , &
   VLNHB    => micfor%VLNHB     , &
   VLWatMicP     => micfor%VLWatMicP      , &
   VLNO3    => micfor%VLNO3     , &
   VLNOB    => micfor%VLNOB     , &
   VLPO4    => micfor%VLPO4     , &
   VLPOB    => micfor%VLPOB     , &
   litrm    => micfor%litrm     , &
   OMEauto    => micstt%OMEauto     , &
   ZNH4S    => micstt%ZNH4S     , &
   ZNH4B    => micstt%ZNH4B     , &
   ZNO3S    => micstt%ZNO3S     , &
   ZNO3B    => micstt%ZNO3B     , &
   CNO3S    => micstt%CNO3S     , &
   CNO3B    => micstt%CNO3B     , &
   CNH4S    => micstt%CNH4S     , &
   CNH4B    => micstt%CNH4B     , &
   CH2P4    => micstt%CH2P4     , &
   CH2P4B   => micstt%CH2P4B    , &
   H2PO4    => micstt%H2PO4     , &
   H2POB    => micstt%H2POB     , &
   CH1P4    => micstt%CH1P4     , &
   CH1P4B   => micstt%CH1P4B    , &
   H1PO4    => micstt%H1PO4     , &
   H1POB    => micstt%H1POB     , &
   RNH4UptkSoilAutor  => micflx%RNH4UptkSoilAutor   , &
   RNH4UptkBandAutor  => micflx%RNH4UptkBandAutor   , &
   RNO3UptkSoilAutor  => micflx%RNO3UptkSoilAutor   , &
   RNO3UptkBandAutor  => micflx%RNO3UptkBandAutor   , &
   NetNH4Mineralize_col   => micflx%NetNH4Mineralize_col    , &
   RH2PO4UptkSoilAutor  => micflx%RH2PO4UptkSoilAutor   , &
   RH2PO4UptkBandAutor  => micflx%RH2PO4UptkBandAutor   , &
   NetPO4Mineralize_col   => micflx%NetPO4Mineralize_col    , &
   RIPO1    => micflx%RIPO1     , &
   RIPB1    => micflx%RIPB1     , &
   RINHOR   => micflx%RINHOR    , &
   RINOOR   => micflx%RINOOR    , &
   RIPOOR   => micflx%RIPOOR    , &
   RH1PO4UptkSoilAutor  => micflx%RH1PO4UptkSoilAutor   , &
   RH1PO4UptkBandAutor  => micflx%RH1PO4UptkBandAutor   , &
   RH1PO4UptkLitrAutor   => micflx%RH1PO4UptkLitrAutor,  &
   RNH4UptkLitrAutor  => micflx%RNH4UptkLitrAutor , &
   RNO3UptkLitrAutor  => micflx%RNO3UptkLitrAutor , &
   RH2PO4UptkLitrAutor  => micflx%RH2PO4UptkLitrAutor   &
  )
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
!     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMN=microbial nonstructural C,N
!     rNCOMC=maximum microbial N:C ratio
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     RINHX=microbially limited NH4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHO,RINHB=substrate-unlimited NH4 mineraln-immobiln
!     VLWatMicP=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     NetNH4Mineralize_col=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  MID3=micpar%get_micb_id(3,NGL)
  RINHP=OMEauto(ielmc,MID3)*rNCOMCff(3,NGL)-OMEauto(ielmn,MID3)
  IF(RINHP.GT.0.0_r8)THEN
    CNH4X=AZMAX1(CNH4S-Z4MN)
    CNH4Y=AZMAX1(CNH4B-Z4MN)
    RINHX=AMIN1(RINHP,BIOA*OMActAutor(NGL)*TFNGff(NGL)*Z4MX)
    RNH4UptkSoilAutor(NGL)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RNH4UptkBandAutor(NGL)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M=Z4MN*VLWatMicP*FNH4S
    ZNHBM=Z4MN*VLWatMicP*FNHBS
    RNH4TransfSoilAutor(NGL)=AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4UptkSoilAutor(NGL))
    RNH4TransfBandAutor(NGL)=AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4UptkBandAutor(NGL))
  ELSE
    RNH4UptkSoilAutor(NGL)=0.0_r8
    RNH4UptkBandAutor(NGL)=0.0_r8
    RNH4TransfSoilAutor(NGL)=RINHP*FNH4S
    RNH4TransfBandAutor(NGL)=RINHP*FNHBS
  ENDIF
  NetNH4Mineralize_col=NetNH4Mineralize_col+(RNH4TransfSoilAutor(NGL)+RNH4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
!     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINOP=NO3 immobilization (+ve) demand
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     min NO3 concentration and Km for NO3 uptake
!     RINOX=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RINOO,RINOB=substrate-unlimited NO3 immobiln
!     VLWatMicP=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     NetNH4Mineralize_col=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
  FNO3S=VLNO3
  FNO3B=VLNOB
  RINOP=AZMAX1(RINHP-RNH4TransfSoilAutor(NGL)-RNH4TransfBandAutor(NGL))
  IF(RINOP.GT.0.0)THEN
    CNO3X=AZMAX1(CNO3S-ZOMN)
    CNO3Y=AZMAX1(CNO3B-ZOMN)
    RINOX=AMIN1(RINOP,BIOA*OMActAutor(NGL)*TFNGff(NGL)*ZOMX)
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
    RNO3TransfSoilAutor(NGL)=RINOP*FNO3S
    RNO3TransfBandAutor(NGL)=RINOP*FNO3B
  ENDIF
  NetNH4Mineralize_col=NetNH4Mineralize_col+(RNO3TransfSoilAutor(NGL)+RNO3TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMP=microbial nonstructural C,P
!     rPCOMC=maximum microbial P:C ratio
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     min H2PO4 concentration and Km for H2PO4 uptake
!     RIPOX=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RIPOO,RIPBO=substrate-unlimited H2PO4 mineraln-immobiln
!     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
!     VOLW=water content
!     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
!     RH2PO4TransfSoilHeter,RH2PO4TransfBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize_col=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
  FH2PS=VLPO4
  FH2PB=VLPOB
  MID3=micpar%get_micb_id(3,NGL)
  RIPOP=(OMEauto(ielmc,MID3)*rPCOMCff(3,NGL)-OMEauto(ielmp,MID3))
  IF(RIPOP.GT.0.0)THEN
    CH2PX=AZMAX1(CH2P4-HPMN)
    CH2PY=AZMAX1(CH2P4B-HPMN)
    RIPOX=AMIN1(RIPOP,BIOA*OMActAutor(NGL)*TFNGff(NGL)*HPMX)
    RH2PO4UptkSoilAutor(NGL)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RH2PO4UptkBandAutor(NGL)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM=HPMN*VLWatMicP*FH2PS
    H2PBM=HPMN*VLWatMicP*FH2PB
    RH2PO4TransfSoilAutor(NGL)=AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RH2PO4UptkSoilAutor(NGL))
    RH2PO4TransfBandAutor(NGL)=AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RH2PO4UptkBandAutor(NGL))
  ELSE
    RH2PO4UptkSoilAutor(NGL)=0.0_r8
    RH2PO4UptkBandAutor(NGL)=0.0_r8
    RH2PO4TransfSoilAutor(NGL)=RIPOP*FH2PS
    RH2PO4TransfBandAutor(NGL)=RIPOP*FH2PB
  ENDIF
  NetPO4Mineralize_col=NetPO4Mineralize_col+(RH2PO4TransfSoilAutor(NGL)+RH2PO4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIP1P=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     min HPO4 concentration and Km for HPO4 uptake
!     RIP1X=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RIPO1,RIPB1=substrate-unlimited HPO4 mineraln-immobiln
!     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
!     VOLW=water content
!     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
!     RH1PO4TransfSoilHeter,RH1PO4TransfBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize_col=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
  FH1PS=VLPO4
  FH1PB=VLPOB
  RIP1P=0.1_r8*AZMAX1(RIPOP-RH2PO4TransfSoilAutor(NGL)-RH2PO4TransfBandAutor(NGL))
  IF(RIP1P.GT.0.0_r8)THEN
    CH1PX=AZMAX1(CH1P4-HPMN)
    CH1PY=AZMAX1(CH1P4B-HPMN)
    RIP1X=AMIN1(RIP1P,BIOA*OMActAutor(NGL)*TFNGff(NGL)*HPMX)
    RH1PO4UptkSoilAutor(NGL)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RH1PO4UptkBandAutor(NGL)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM=HPMN*VLWatMicP*FH1PS
    H1PBM=HPMN*VLWatMicP*FH1PB
    RH1PO4TransfSoilAutor(NGL)=AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RH1PO4UptkSoilAutor(NGL))
    RH1PO4TransfBandAutor(NGL)=AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RH1PO4UptkBandAutor(NGL))
  ELSE
    RH1PO4UptkSoilAutor(NGL)=0.0_r8
    RH1PO4UptkBandAutor(NGL)=0.0_r8
    RH1PO4TransfSoilAutor(NGL)=RIP1P*FH1PS
    RH1PO4TransfBandAutor(NGL)=RIP1P*FH1PB
  ENDIF
  NetPO4Mineralize_col=NetPO4Mineralize_col+(RH1PO4TransfSoilAutor(NGL)+RH1PO4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RINHPR=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHOR=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M=NH4 not available for uptake
!     FNH4XR=fractions of biological NH4 demand
!     RNH4TransfLitrHeter=substrate-limited NH4 mineraln-immobiln
!     NetNH4Mineralize_col=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  IF(litrm)THEN
    RINHPR=RINHP-RNH4TransfSoilAutor(NGL)-RNO3TransfSoilAutor(NGL)
    IF(RINHPR.GT.0.0_r8)THEN
      CNH4X=AZMAX1(CNH4S-Z4MN)
      CNH4Y=AZMAX1(CNH4B-Z4MN)
      RNH4UptkLitrAutor(NGL)=AMIN1(RINHPR,BIOA*OMActAutor(NGL)*TFNGff(NGL)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VLWatMicP
      RNH4TransfLitrAutor(NGL)=AMIN1(FNH4XRff(NGL)*AZMAX1((ZNH4T-ZNH4M)),RNH4UptkLitrAutor(NGL))
    ELSE
      RNH4UptkLitrAutor(NGL)=0.0_r8
      RNH4TransfLitrAutor(NGL)=RINHPR
    ENDIF
    NetNH4Mineralize_col=NetNH4Mineralize_col+RNH4TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RINOPR=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     minimum NO3 concentration and Km for NO3 uptake
!     RINOOR=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RNO3TransfLitrHeter=substrate-unlimited NO3 immobiln
!     VLWatMicP=water content
!     ZNO3M=NO3 not available for uptake
!     FNO3XR=fraction of biological NO3 demand
!     RNO3TransfLitrHeter=substrate-limited NO3 immobiln
!     NetNH4Mineralize_col=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
    RINOPR=AZMAX1(RINHPR-RNH4TransfLitrAutor(NGL))
    IF(RINOPR.GT.0.0_r8)THEN
      CNO3X=AZMAX1(CNO3S-ZOMN)
      CNO3Y=AZMAX1(CNO3B-ZOMN)
      RNO3UptkLitrAutor(NGL)=AMAX1(RINOPR,BIOA*OMActAutor(NGL)*TFNGff(NGL)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VLWatMicP
      RNO3TransfLitrAutor(NGL)=AMIN1(FNO3XRff(NGL)*AZMAX1((ZNO3T-ZNO3M)),RNO3UptkLitrAutor(NGL))
    ELSE
      RNO3UptkLitrAutor(NGL)=0.0_r8
      RNO3TransfLitrAutor(NGL)=RINOPR
    ENDIF
    NetNH4Mineralize_col=NetNH4Mineralize_col+RNO3TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RIPOPR=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     minimum H2PO4 concentration and Km for H2PO4 uptake
!     RIPOOR=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RIPOOR=substrate-unlimited H2PO4 mineraln-immobiln
!     VLWatMicP=water content
!     H2P4M=H2PO4 not available for uptake
!     FPO4XR=fractions of biological H2PO4 demand
!     RH2PO4TransfLitrHeter=substrate-limited H2PO4 mineraln-immobiln
!     NetPO4Mineralize_col=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    RIPOPR=RIPOP-RH2PO4TransfSoilAutor(NGL)
    IF(RIPOPR.GT.0.0_r8)THEN
      CH2PX=AZMAX1(CH2P4-HPMN)
      CH2PY=AZMAX1(CH2P4B-HPMN)
      RH2PO4UptkLitrAutor(NGL)=AMIN1(RIPOPR,BIOA*OMActAutor(NGL)*TFNGff(NGL)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M=HPMN*VLWatMicP
      RH2PO4TransfLitrAutor(NGL)=AMIN1(FPO4XRff(NGL)*AZMAX1((H2P4T-H2P4M)),RH2PO4UptkLitrAutor(NGL))
    ELSE
      RH2PO4UptkLitrAutor(NGL)=0.0_r8
      RH2PO4TransfLitrAutor(NGL)=RIPOPR
    ENDIF
    NetPO4Mineralize_col=NetPO4Mineralize_col+RH2PO4TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RIP1PR=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     minimum HPO4 concentration and Km for HPO4 uptake
!     RIPO1R=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RIPO1R=substrate-unlimited HPO4 mineraln-immobiln
!     VLWatMicP=water content
!     H1P4M=HPO4 not available for uptake
!     FP14XR=fraction of biological HPO4 demand
!     RH1PO4TransfLitrHeter=substrate-limited HPO4 minereraln-immobiln
!     NetPO4Mineralize_col=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RIP1PR=0.1*AZMAX1(RIPOPR-RH2PO4TransfLitrAutor(NGL))
    IF(RIP1PR.GT.0.0_r8)THEN
      CH1PX=AZMAX1(CH1P4-HPMN)
      CH1PY=AZMAX1(CH1P4B-HPMN)
      RH1PO4UptkLitrAutor(NGL)=AMIN1(RIP1PR,BIOA*OMActAutor(NGL)*TFNGff(NGL)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M=HPMN*VLWatMicP
      RH1PO4TransfLitrAutor(NGL)=AMIN1(FP14XRff(NGL)*AZMAX1((H1P4T-H1P4M)),RH1PO4UptkLitrAutor(NGL))
    ELSE
      RH1PO4UptkLitrAutor(NGL)=0.0_r8
      RH1PO4TransfLitrAutor(NGL)=RIP1PR
    ENDIF
    NetPO4Mineralize_col=NetPO4Mineralize_col+RH1PO4TransfLitrAutor(NGL)
  ENDIF
  end associate
  end subroutine BiomassMineralizationff
!------------------------------------------------------------------------------------------

  subroutine GatherMicrobialRespirationff(NGL,N,RMOMK,micfor,micstt,RGOMT,RXOMT,RMOMT,&
    nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  real(r8), intent(out) :: RGOMT,RXOMT
  real(r8), intent(out) :: RMOMT
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
  integer :: MID1
!     begin_execution
  associate(                      &
    TFNRff  => nmics%TFNRff,      &
    OMN2ff  => nmics%OMN2ff,      &
    RMaintCompAutor => nmicf%RMaintCompAutor , &
    RGOMOff  => nmicf%RGOMOff, &
    RGN2Fff => nmicf%RGN2Fff,     &
    RN2FixAutor  => nmicf%RN2FixAutor, &
    pH      => micfor%pH    , &
    OMEauto   => micstt%OMEauto       &

  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TFNR=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  MID1=micpar%get_micb_id(1,NGL)
  FPH=1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX=RMOM*TFNRff(NGL)*FPH
  RMaintCompAutor(1,NGL)=OMEauto(ielmn,MID1)*RMOMX*RMOMK(1)
  RMaintCompAutor(2,NGL)=OMN2ff(NGL)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMOMT=total maintenance respiration
!     RGOMT=growth respiration
!     RXOMT=senescence respiration
!
  RMOMT=RMaintCompAutor(1,NGL)+RMaintCompAutor(2,NGL)
  RGOMT=AZMAX1(RGOMOff(NGL)-RMOMT)
  RXOMT=AZMAX1(RMOMT-RGOMOff(NGL))

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
!     RGOMT=growth respiration
!     RGN2F=respiration for N2 fixation
!     CZ2GS=aqueous N2 concentration
!     ZFKM=Km for N2 uptake
!     OMGR*OMC(3,NGL,N,K)=nonstructural C limitation to RGN2F
!     RN2FixHeter=N2 fixation rate
!

  RN2FixAutor(NGL)=0.0_r8
  RGN2Fff(NGL)=0.0_r8
  end associate
  end subroutine GatherMicrobialRespirationff

!------------------------------------------------------------------------------------------

  subroutine MicrobialAnabolicUpdateff(micfor,micstt,nmicf)

  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: CGROMC
  integer :: N,M,NGL,MID,MID3,NE
  associate(                      &
    CGOMEautor    => nmicf%CGOMEautor,  &
    CGOSEautor  => nmicf%CGOSEautor, &
    RGN2Fff  => nmicf%RGN2Fff, &
    RGOMOff  => nmicf%RGOMOff, &
    RGOMDff  => nmicf%RGOMDff, &
    RNO3TransfSoilAutor  => nmicf%RNO3TransfSoilAutor, &
    RCO2ProdAutor  => nmicf%RCO2ProdAutor, &
    RH2PO4TransfSoilAutor  => nmicf%RH2PO4TransfSoilAutor, &
    RNH4TransfBandAutor => nmicf%RNH4TransfBandAutor , &
    RNO3TransfBandAutor => nmicf%RNO3TransfBandAutor , &
    RH2PO4TransfBandAutor => nmicf%RH2PO4TransfBandAutor , &
    RHOMEautor => nmicf%RHOMEautor , &
    RHMMEautor  => nmicf%RHMMEautor, &
    RN2FixAutor  => nmicf%RN2FixAutor, &
    RXOMEautor  => nmicf%RXOMEautor, &
    R3OMEautor  => nmicf%R3OMEautor, &
    RXMMEautor  => nmicf%RXMMEautor, &
    R3MMEautor  => nmicf%R3MMEautor, &
    RNH4TransfLitrAutor  => nmicf%RNH4TransfLitrAutor,  &
    RNO3TransfLitrAutor  => nmicf%RNO3TransfLitrAutor,  &
    RH2PO4TransfLitrAutor  => nmicf%RH2PO4TransfLitrAutor,  &
    RH1PO4TransfSoilAutor   => nmicf%RH1PO4TransfSoilAutor ,  &
    RH1PO4TransfBandAutor   => nmicf%RH1PO4TransfBandAutor,  &
    RH1PO4TransfLitrAutor  => nmicf%RH1PO4TransfLitrAutor,  &
    RNH4TransfSoilAutor   => nmicf%RNH4TransfSoilAutor,  &
    litrm     => micfor%litrm  , &
    CFOMC     => micfor%CFOMC  , &
    CFOMCU    => micfor%CFOMCU , &
    SolidOM       => micstt%SolidOM    , &
    OMEauto     => micstt%OMEauto  , &
    OSC14U    => micstt%OSC14U , &
    OSN14U    => micstt%OSN14U , &
    OSP14U    => micstt%OSP14U , &
    OSC24U    => micstt%OSC24U , &
    OSN24U    => micstt%OSN24U , &
    OSP24U    => micstt%OSP24U , &
    JGniA     => micpar%JGniA  , &
    JGnfA     => micpar%JGnfA  , &
    NumMicbFunGroups      => micpar%NumMicbFunGroups   , &
    icarbhyro => micpar%icarbhyro, &
    iprotein  => micpar%iprotein , &
    k_POM     => micpar%k_POM  , &
    is_activef_micb => micpar%is_activef_micb  &
  )
  DO  N=1,NumMicbFunGroups
    IF(is_activef_micb(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        DO  M=1,2
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            OMEauto(NE,MID)=OMEauto(NE,MID)+CGOSEautor(NE,M,NGL)-RXOMEautor(NE,M,NGL)-RXMMEautor(NE,M,NGL)
          ENDDO

!     HUMIFICATION PRODUCTS
!
!     CFOMC=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence LitrFall C,N,P to humus
!
          IF(.not.litrm)THEN
            DO NE=1,NumPlantChemElms
              SolidOM(NE,iprotein,k_POM)=SolidOM(NE,iprotein,k_POM)+CFOMC(1)*(RHOMEautor(NE,M,NGL)+RHMMEautor(NE,M,NGL))
              SolidOM(NE,icarbhyro,k_POM)=SolidOM(NE,icarbhyro,k_POM)+CFOMC(2)*(RHOMEautor(NE,M,NGL)+RHMMEautor(NE,M,NGL))
            ENDDO
          ELSE
            OSC14U=OSC14U+CFOMCU(1)*(RHOMEautor(ielmc,M,NGL)+RHMMEautor(ielmc,M,NGL))
            OSN14U=OSN14U+CFOMCU(1)*(RHOMEautor(ielmn,M,NGL)+RHMMEautor(ielmn,M,NGL))
            OSP14U=OSP14U+CFOMCU(1)*(RHOMEautor(ielmp,M,NGL)+RHMMEautor(ielmp,M,NGL))
            OSC24U=OSC24U+CFOMCU(2)*(RHOMEautor(ielmc,M,NGL)+RHMMEautor(ielmc,M,NGL))
            OSN24U=OSN24U+CFOMCU(2)*(RHOMEautor(ielmn,M,NGL)+RHMMEautor(ielmn,M,NGL))
            OSP24U=OSP24U+CFOMCU(2)*(RHOMEautor(ielmp,M,NGL)+RHMMEautor(ielmp,M,NGL))
          ENDIF
        ENDDO
!
!     INPUTS TO NONSTRUCTURAL POOLS
!
!     CGOMC=total DOC+acetate uptake
!     RGOMO=total respiration
!     RGOMD=respiration for denitrifcation
!     RGN2F=respiration for N2 fixation
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
        CGROMC=CGOMEautor(ielmc,NGL)-RGOMOff(NGL)-RGOMDff(NGL)-RGN2Fff(NGL)
        RCO2ProdAutor(NGL)=RCO2ProdAutor(NGL)+RGN2Fff(NGL)
        MID3=micpar%get_micb_id(3,NGL)
        DO M=1,2
          OMEauto(ielmc,MID3)=OMEauto(ielmc,MID3)-CGOSEautor(ielmc,M,NGL)+R3OMEautor(ielmc,M,NGL)
          OMEauto(ielmn,MID3)=OMEauto(ielmn,MID3)-CGOSEautor(ielmn,M,NGL)+R3OMEautor(ielmn,M,NGL)+R3MMEautor(ielmn,M,NGL)
          OMEauto(ielmp,MID3)=OMEauto(ielmp,MID3)-CGOSEautor(ielmp,M,NGL)+R3OMEautor(ielmp,M,NGL)+R3MMEautor(ielmp,M,NGL)
          RCO2ProdAutor(NGL)=RCO2ProdAutor(NGL)+R3MMEautor(ielmc,M,NGL)
        ENDDO
        OMEauto(ielmc,MID3)=OMEauto(ielmc,MID3)+CGROMC
        OMEauto(ielmn,MID3)=OMEauto(ielmn,MID3)+CGOMEautor(ielmn,NGL) &
          +RNH4TransfSoilAutor(NGL)+RNH4TransfBandAutor(NGL)+RNO3TransfSoilAutor(NGL)+RNO3TransfBandAutor(NGL)+RN2FixAutor(NGL)
        OMEauto(ielmp,MID3)=OMEauto(ielmp,MID3)+CGOMEautor(ielmp,NGL) &
          +RH2PO4TransfSoilAutor(NGL)+RH2PO4TransfBandAutor(NGL)+RH1PO4TransfSoilAutor(NGL)+RH1PO4TransfBandAutor(NGL)
        IF(litrm)THEN
          OMEauto(ielmn,MID3)=OMEauto(ielmn,MID3)+RNH4TransfLitrAutor(NGL)+RNO3TransfLitrAutor(NGL)
          OMEauto(ielmp,MID3)=OMEauto(ielmp,MID3)+RH2PO4TransfLitrAutor(NGL)+RH1PO4TransfLitrAutor(NGL)
        ENDIF
      enddo
    ENDIF
  ENDDO
  end associate
  end subroutine MicrobialAnabolicUpdateff


end module MicAutoCPLXMod
