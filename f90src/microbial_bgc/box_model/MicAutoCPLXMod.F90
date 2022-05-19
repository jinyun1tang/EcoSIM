module MicAutoCPLXMod
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : safe_adb
  use MicForcTypeMod, only : micforctype
  use MicFluxTypeMod, only: micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use NitroDiagTypes
  use MicBGCPars, only : micpar
  use EcoSIMSolverPar
  use EcoSimConst
  use NitroPars
  implicit none

  private
  character(len=*), parameter :: mod_filename = __FILE__

  public :: ActiveMicrobesff
  public :: MicrobialAnabolicUpdateff
  contains

!------------------------------------------------------------------------------------------
  subroutine ActiveMicrobesff(NGL,N,VOLWZ,XCO2,TFNX,WFNG,TOMCNK,OXKX,TOMA,TOMN,RH2GZ,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: VOLWZ
  real(r8), intent(in):: OXKX,tomcnk(2),WFNG,TFNX,XCO2
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
    FOMAff => nmics%FOMAff,       &
    FOMNff => nmics%FOMNff,       &
    FOMKff => nmics%FOMKff,       &
    OMAff => nmics%OMAff,         &
    RUPOXff  => nmicf%RUPOXff,    &
    RGOMOff  => nmicf%RGOMOff,    &
    ROXYMff => nmicf%ROXYMff ,    &
    ROXYPff=> nmicf%ROXYPff  ,    &
    ROXYOff => nmicf%ROXYOff ,    &
    RDNO3ff => nmicf%RDNO3ff ,    &
    RDNOBff => nmicf%RDNOBff ,    &
    RDNO2ff => nmicf%RDNO2ff ,    &
    RDN2Bff => nmicf%RDN2Bff ,    &
    RDN2Off => nmicf%RDN2Off ,    &
    RGOMDff  => nmicf%RGOMDff,    &
    RH2GXff  => nmicf%RH2GXff,    &
    RGOMYff  => nmicf%RGOMYff,    &
    RCO2Xff  => nmicf%RCO2Xff,    &
    RCH3Xff  => nmicf%RCH3Xff,    &
    RCH4Xff  => nmicf%RCH4Xff,    &
    TOMK  => ncplxs%TOMK     ,    &
    jcplx  => micpar%jcplx   ,    &
    nf_amonia_oxi => micpar%nf_amonia_oxi, &
    nf_nitrite_oxi  => micpar%nf_nitrite_oxi, &
    nf_hydro_methang  => micpar%nf_hydro_methang, &
    nf_aero_methanot  => micpar%nf_aero_methanot, &
    ZEROS  => micfor%ZEROS  ,    &
    BKVL  => micfor%BKVL    ,    &
    litrm => micfor%litrm  , &
    VOLX  => micfor%VOLX, &
    ORGC   => micfor%ORGC  ,     &
    ROXYSff => micflx%ROXYSff  &
  )
! FOMA,FOMN=fraction of total active biomass C,N in each N and K

  IF(TOMA.GT.ZEROS)THEN
    FOMAff(NGL,N)=OMAff(NGL,N)/TOMA
  ELSE
    FOMAff(NGL,N)=1.0
  ENDIF
  IF(TOMN.GT.ZEROS)THEN
    FOMNff(NGL,N)=OMAff(NGL,N)/TOMN
  ELSE
    FOMNff(NGL,N)=1.0_r8
  ENDIF

  K=jcplx
  IF(TOMK(K).GT.ZEROS)THEN
    FOMKff(NGL,N)=OMAff(NGL,N)/TOMK(K)
  ELSE
    FOMKff(NGL,N)=1.0_r8
  ENDIF
!
  !     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
  !
  !     COMC=microbial C concentration relative to substrate
  !     SPOMK=effect of microbial C concentration on microbial decay
  !     RMOMK=effect of microbial C concentration on maintenance respn
  !
  ORGCL=AMIN1(1.0E+05_r8*BKVL,ORGC)
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
  if (N.eq.nf_amonia_oxi)then
!   NH3 OXIDIZERS
    call NH3OxidizerCatabolism(NGL,N,XCO2,VOLWZ,TFNX,ECHZ,RGOMP,RVOXP,&
      RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.nf_nitrite_oxi)then
!     NO2 OXIDIZERS
    call NO2OxidizerCatabolism(NGL,N,XCO2,ECHZ,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      micfor,micstt,naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.nf_hydro_methang)then
!     H2TROPHIC METHANOGENS
    call H2MethanogensCatabolism(NGL,N,ECHZ,RGOMP,XCO2,micfor,micstt,&
      naqfdiag,nmicf,nmics,micflx)

  elseif (N.eq.nf_aero_methanot)then
!     METHANOTROPHS
    call MethanotrophCatabolism(NGL,N,ECHZ,RGOMP,&
      RVOXP,RVOXPA,RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  else
    RGOMP=0.0_r8
    ROXYMff(NGL,N)=0.0_r8
    ROXYPff(NGL,N)=0.0_r8
    ROXYSff(NGL,N)=0.0_r8
  ENDif
!
!  write(*,*)'O2 UPTAKE BY AEROBES'
!
! RUPOX, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
! RUPMX=O2-unlimited rate of O2 uptake
! FOXYX=fraction of O2 uptake by N,K relative to total
! XNPG=1/(NPH*NPT)
! ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
! OLSGL=aqueous O2 diffusivity
! OXYG,OXYS=gaseous, aqueous O2 amounts
! FLQRQ,FLQRI=surface water flux from precipitation, irrigation
! COXR,COXQ=O2 concentration in FLQRQ,FLQRI
!
  RUPOXff(NGL,N)=0.0_r8

  if (N.eq.nf_amonia_oxi .or. N.eq.nf_nitrite_oxi .or. N.eq.nf_aero_methanot)then
!   write(*,*)'AerobsO2Uptake'
    call AerobsO2Uptakeff(NGL,N,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      micfor,micstt,nmicf,nmics,micflx)
  elseif (N.eq.nf_hydro_methang)then
    RGOMOff(NGL,N)=RGOMP
    RCO2Xff(NGL,N)=0.0_r8
    RCH3Xff(NGL,N)=0.0_r8
    RCH4Xff(NGL,N)=RGOMOff(NGL,N)
    ROXYOff(NGL,N)=ROXYMff(NGL,N)
    RH2GXff(NGL,N)=0.0_r8
    RH2GZ=0.667_r8*RGOMOff(NGL,N)
  ENDIF
!
!     AUTOTROPHIC DENITRIFICATION
!
  IF(N.EQ.nf_amonia_oxi.AND.ROXYMff(NGL,N).GT.0.0_r8 &
    .AND.(.not.litrm.OR.VOLX.GT.ZEROS))THEN
    call AutotrophDenitrificCatabolism(NGL,N,XCO2,VOLWZ,micfor,micstt,&
      naqfdiag,nmicf,nmics,micflx)
  ELSE
    RDNO3ff(NGL,N)=0.0_r8
    RDNOBff(NGL,N)=0.0_r8
    RDNO2ff(NGL,N)=0.0_r8
    RDN2Bff(NGL,N)=0.0_r8
    RDN2Off(NGL,N)=0.0_r8
    RGOMYff(NGL,N)=0.0_r8
    RGOMDff(NGL,N)=0.0_r8
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
    FOMAff    => nmics%FOMAff,    &
    FOMKff    => nmics%FOMKff,    &
    FNH4XRff  =>nmicf%FNH4XRff,   &
    FNO3XRff  =>nmicf%FNO3XRff,   &
    FP14XRff  => nmicf%FP14XRff,  &
    FPO4XRff  => nmicf%FPO4XRff,  &
    ROXYY     => micfor%ROXYY,    &
    RNH4Y     => micfor%RNH4Y,    &
    RNHBY     => micfor%RNHBY,    &
    RNO3Y     => micfor%RNO3Y,    &
    RN3BY     => micfor%RN3BY,    &
    RPO4Y     => micfor%RPO4Y,    &
    RPOBY     => micfor%RPOBY,    &
    RP14Y     => micfor%RP14Y,    &
    RP1BY     => micfor%RP1BY,    &
    ROQCY     => micfor%ROQCY,    &
    ROQAY     => micfor%ROQAY,    &
    BKVL0     => micfor%BKVL0,   &
    Lsurf     => micfor%Lsurf,   &
    litrm     => micfor%litrm,   &
    ZEROS     => micfor%ZEROS,   &
    VLNO3     => micfor%VLNO3,   &
    VLNOB     => micfor%VLNOB,   &
    VLPO4     => micfor%VLPO4,   &
    VLPOB     => micfor%VLPOB,   &
    VLNHB     => micfor%VLNHB,   &
    VLNH4     => micfor%VLNH4,   &
    RIPO1Rff  => micflx%RIPO1Rff, &
    RIPOORff  => micflx%RIPOORff, &
    RINOORff  => micflx%RINOORff, &
    RINHORff  => micflx%RINHORff, &
    RIPB1ff   => micflx%RIPB1ff, &
    RIPO1ff   => micflx%RIPO1ff, &
    RIPBOff   => micflx%RIPBOff, &
    RIPOOff   => micflx%RIPOOff, &
    RINOBff   => micflx%RINOBff, &
    RINOOff   => micflx%RINOOff, &
    RINHBff   => micflx%RINHBff, &
    RINHOff   => micflx%RINHOff, &
    ROXYSff   => micflx%ROXYSff  &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(ROXYY.GT.ZEROS)THEN
    FOXYX=AMAX1(FMN,ROXYSff(NGL,N)/ROXYY)
  ELSE
    FOXYX=AMAX1(FMN,FOMAff(NGL,N))
  ENDIF
  IF(RNH4Y.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RINHOff(NGL,N)/RNH4Y)
  ELSE
    FNH4X=AMAX1(FMN,FOMAff(NGL,N)*VLNH4)
  ENDIF
  IF(RNHBY.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RINHBff(NGL,N)/RNHBY)
  ELSE
    FNB4X=AMAX1(FMN,FOMAff(NGL,N)*VLNHB)
  ENDIF
  IF(RNO3Y.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RINOOff(NGL,N)/RNO3Y)
  ELSE
    FNO3X=AMAX1(FMN,FOMAff(NGL,N)*VLNO3)
  ENDIF
  IF(RN3BY.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RINOBff(NGL,N)/RN3BY)
  ELSE
    FNB3X=AMAX1(FMN,FOMAff(NGL,N)*VLNOB)
  ENDIF
  IF(RPO4Y.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RIPOOff(NGL,N)/RPO4Y)
  ELSE
    FPO4X=AMAX1(FMN,FOMAff(NGL,N)*VLPO4)
  ENDIF
  IF(RPOBY.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RIPBOff(NGL,N)/RPOBY)
  ELSE
    FPOBX=AMAX1(FMN,FOMAff(NGL,N)*VLPOB)
  ENDIF
  IF(RP14Y.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RIPO1ff(NGL,N)/RP14Y)
  ELSE
    FP14X=AMAX1(FMN,FOMAff(NGL,N)*VLPO4)
  ENDIF
  IF(RP1BY.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RIPB1ff(NGL,N)/RP1BY)
  ELSE
    FP1BX=AMAX1(FMN,FOMAff(NGL,N)*VLPOB)
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
    IF(RNH4Y.GT.ZEROS)THEN
      FNH4XRff(NGL,N)=AMAX1(FMN,RINHORff(NGL,N)/RNH4Y)
    ELSE
      FNH4XRff(NGL,N)=AMAX1(FMN,FOMKff(NGL,N))
    ENDIF
    IF(RNO3Y.GT.ZEROS)THEN
      FNO3XRff(NGL,N)=AMAX1(FMN,RINOORff(NGL,N)/RNO3Y)
    ELSE
      FNO3XRff(NGL,N)=AMAX1(FMN,FOMKff(NGL,N))
    ENDIF
    IF(RPO4Y.GT.ZEROS)THEN
      FPO4XRff(NGL,N)=AMAX1(FMN,RIPOORff(NGL,N)/RPO4Y)
    ELSE
      FPO4XRff(NGL,N)=AMAX1(FMN,FOMKff(NGL,N))
    ENDIF
    IF(RP14Y.GT.ZEROS)THEN
      FP14XRff(NGL,N)=AMAX1(FMN,RIPO1Rff(NGL,N)/RP14Y)
    ELSE
      FP14XRff(NGL,N)=AMAX1(FMN,FOMKff(NGL,N))
    ENDIF
  ENDIF
  IF(Lsurf.AND.BKVL0.GT.ZEROS)THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4XRff(NGL,N)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3XRff(NGL,N)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4XRff(NGL,N)
    naqfdiag%TFP14X=naqfdiag%TFP14X+FP14XRff(NGL,N)
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
  integer :: M,K
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: CXC
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FRM
!     begin_execution
  associate(                  &
    CNOMAff  => nmics%CNOMAff,    &
    CPOMAff  => nmics%CPOMAff,    &
    TFNGff   => nmics%TFNGff ,    &
    CGOMCff  => nmicf%CGOMCff,    &
    CGOMNff  => nmicf%CGOMNff,    &
    CGOMPff  => nmicf%CGOMPff,    &
    CGOQCff  => nmicf%CGOQCff,    &
    CGOACff  => nmicf%CGOACff,    &
    CGOMSff  => nmicf%CGOMSff,    &
    CGONSff  => nmicf%CGONSff,    &
    CGOPSff  => nmicf%CGOPSff,    &
    RGOMOff  => nmicf%RGOMOff,    &
    RGOMDff  => nmicf%RGOMDff,    &
    RMOMCff => nmicf%RMOMCff ,    &
    RDOMCff  => nmicf%RDOMCff,    &
    RDOMNff => nmicf%RDOMNff ,    &
    RDOMPff => nmicf%RDOMPff ,    &
    RHOMCff => nmicf%RHOMCff ,    &
    RHOMNff => nmicf%RHOMNff ,    &
    RHOMPff => nmicf%RHOMPff ,    &
    RCOMCff  => nmicf%RCOMCff,    &
    RCOMNff  => nmicf%RCOMNff,    &
    RCOMPff  => nmicf%RCOMPff,    &
    RDMMCff  => nmicf%RDMMCff,    &
    RHMMCff  => nmicf%RHMMCff,    &
    RCMMCff  => nmicf%RCMMCff,    &
    RDMMNff  => nmicf%RDMMNff,    &
    RHMMNff  => nmicf%RHMMNff,    &
    RCMMNff  => nmicf%RCMMNff,    &
    RDMMPff  => nmicf%RDMMPff,    &
    RHMMPff  => nmicf%RHMMPff,    &
    RCMMPff  => nmicf%RCMMPff,    &
    RXOMCff  => nmicf%RXOMCff,    &
    RXOMNff  => nmicf%RXOMNff,    &
    RXOMPff  => nmicf%RXOMPff,    &
    R3OMCff  => nmicf%R3OMCff,    &
    R3OMNff  => nmicf%R3OMNff,    &
    R3OMPff  => nmicf%R3OMPff,    &
    RXMMCff  => nmicf%RXMMCff,    &
    RXMMNff  => nmicf%RXMMNff,    &
    RXMMPff  => nmicf%RXMMPff,    &
    R3MMCff  => nmicf%R3MMCff,    &
    R3MMNff  => nmicf%R3MMNff,    &
    R3MMPff  => nmicf%R3MMPff,    &
    RGN2Fff  => nmicf%RGN2Fff,    &
    TCGOQC  => ncplxf%TCGOQC ,    &
    TCGOAC  => ncplxf%TCGOAC ,    &
    TCGOMN  => ncplxf%TCGOMN ,    &
    TCGOMP  => ncplxf%TCGOMP ,    &
    CNQ     => ncplxs%CNQ,        &
    CPQ     => ncplxs%CPQ,        &
    CNOMCff  => micpar%CNOMCff,   &
    CPOMCff  => micpar%CPOMCff,   &
    FL       => micpar%FL       , &
    ZEROS     => micfor%ZEROS   , &
    ZERO    => micfor%ZERO      , &
    OMCff    => micstt%OMCff    , &
    OMNff    => micstt%OMNff    , &
    OMPff    => micstt%OMPff    , &
    EHUM    => micstt%EHUM        &
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

  CGOMX=AMIN1(RMOMT,RGOMOff(NGL,N))+RGN2Fff(NGL,N) &
    +(RGOMT-RGN2Fff(NGL,N))/ECHZ
  CGOMD=RGOMDff(NGL,N)/ENOX
  CGOMCff(NGL,N)=CGOMX+CGOMD
  CGOQCff(NGL,N)=CGOMX+CGOMD
  CGOACff(NGL,N)=0.0_r8
  CGOMNff(NGL,N)=0.0_r8
  CGOMPff(NGL,N)=0.0_r8
  K=micpar%jcplx
  TCGOQC(K)=TCGOQC(K)+CGOQCff(NGL,N)
  TCGOAC(K)=TCGOAC(K)+CGOACff(NGL,N)
  TCGOMN(K)=TCGOMN(K)+CGOMNff(NGL,N)
  TCGOMP(K)=TCGOMP(K)+CGOMPff(NGL,N)
!
!     TRANSFER UPTAKEN C,N,P FROM STORAGE TO ACTIVE BIOMASS
!
!     OMC,OMN,OMP=nonstructural C,N,P
!     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
!     CNOMC,CPOMC=maximum microbial N:C, P:C ratios
!     RCCC,RCCN,RCCP=C,N,P recycling fractions
!     RCCZ,RCCY=min, max C recycling fractions
!     RCCX,RCCQ=max N,P recycling fractions
!
  IF(OMCff(3,NGL,N).GT.ZEROS &
    .AND.OMCff(1,NGL,N).GT.ZEROS)THEN
    CCC=AMAX1(0.0_r8,AMIN1(1.0_r8 &
      ,OMNff(3,NGL,N)/(OMNff(3,NGL,N) &
      +OMCff(3,NGL,N)*CNOMCff(3,NGL,N)) &
      ,OMPff(3,NGL,N)/(OMPff(3,NGL,N) &
      +OMCff(3,NGL,N)*CPOMCff(3,NGL,N))))
    CXC=OMCff(3,NGL,N)/OMCff(1,NGL,N)
    C3C=1.0_r8/(1.0_r8+CXC/CKC)
    CNC=AMAX1(0.0_r8,AMIN1(1.0_r8 &
      ,OMCff(3,NGL,N)/(OMCff(3,NGL,N) &
      +OMNff(3,NGL,N)/CNOMCff(3,NGL,N))))
    CPC=AMAX1(0.0_r8,AMIN1(1.0_r8 &
      ,OMCff(3,NGL,N)/(OMCff(3,NGL,N) &
      +OMPff(3,NGL,N)/CPOMCff(3,NGL,N))))
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
  CGOMZ=TFNGff(NGL,N)*OMGR*AMAX1(0.0_r8,OMCff(3,NGL,N))
  DO  M=1,2
    CGOMSff(M,NGL,N)=FL(M)*CGOMZ
    IF(OMCff(3,NGL,N).GT.ZEROS)THEN
      CGONSff(M,NGL,N)=AMIN1(FL(M)*AMAX1(0.0_r8,OMNff(3,NGL,N)) &
        ,CGOMSff(M,NGL,N)*OMNff(3,NGL,N)/OMCff(3,NGL,N))
      CGOPSff(M,NGL,N)=AMIN1(FL(M)*AMAX1(0.0_r8,OMPff(3,NGL,N)) &
        ,CGOMSff(M,NGL,N)*OMPff(3,NGL,N)/OMCff(3,NGL,N))
    ELSE
      CGONSff(M,NGL,N)=0.0_r8
      CGOPSff(M,NGL,N)=0.0_r8
    ENDIF
!
!     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
!     RATE, TEMPERATURE
!
!     SPOMX=rate constant for microbial decomposition
!     SPOMC=basal decomposition rate
!     SPOMK=effect of low microbial C concentration on microbial decay
!     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
!     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall
!     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
!
    SPOMX=SQRT(TFNGff(NGL,N))*SPOMC(M)*SPOMK(M)
    RXOMCff(M,NGL,N)=AMAX1(0.0_r8,OMCff(M,NGL,N)*SPOMX)
    RXOMNff(M,NGL,N)=AMAX1(0.0_r8,OMNff(M,NGL,N)*SPOMX)
    RXOMPff(M,NGL,N)=AMAX1(0.0_r8,OMPff(M,NGL,N)*SPOMX)
    RDOMCff(M,NGL,N)=RXOMCff(M,NGL,N)*(1.0_r8-RCCC)
    RDOMNff(M,NGL,N)=RXOMNff(M,NGL,N)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDOMPff(M,NGL,N)=RXOMPff(M,NGL,N)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
    R3OMCff(M,NGL,N)=RXOMCff(M,NGL,N)-RDOMCff(M,NGL,N)
    R3OMNff(M,NGL,N)=RXOMNff(M,NGL,N)-RDOMNff(M,NGL,N)
    R3OMPff(M,NGL,N)=RXOMPff(M,NGL,N)-RDOMPff(M,NGL,N)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!
    RHOMCff(M,NGL,N)=AMAX1(0.0_r8,RDOMCff(M,NGL,N)*EHUM)
    RHOMNff(M,NGL,N)=AMAX1(0.0_r8,RDOMNff(M,NGL,N)*EHUM)
    RHOMPff(M,NGL,N)=AMAX1(0.0_r8,RDOMPff(M,NGL,N)*EHUM)
!
!     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
!
    RCOMCff(M,NGL,N)=RDOMCff(M,NGL,N)-RHOMCff(M,NGL,N)
    RCOMNff(M,NGL,N)=RDOMNff(M,NGL,N)-RHOMNff(M,NGL,N)
    RCOMPff(M,NGL,N)=RDOMPff(M,NGL,N)-RHOMPff(M,NGL,N)
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
!     RDMMC,RDMMN,RDMMP=microbial C,N,P litterfall from senescence
!     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!
  IF(RXOMT.GT.ZEROS.AND.RMOMT.GT.ZEROS &
    .AND.RCCC.GT.ZERO)THEN
    FRM=RXOMT/RMOMT
    DO  M=1,2
      RXMMCff(M,NGL,N)=AMIN1(OMCff(M,NGL,N) &
        ,AMAX1(0.0_r8,FRM*RMOMCff(M,NGL,N)/RCCC))
      RXMMNff(M,NGL,N)=AMIN1(OMNff(M,NGL,N) &
        ,AMAX1(0.0_r8,RXMMCff(M,NGL,N)*CNOMAff(NGL,N)))
      RXMMPff(M,NGL,N)=AMIN1(OMPff(M,NGL,N) &
        ,AMAX1(0.0_r8,RXMMCff(M,NGL,N)*CPOMAff(NGL,N)))
      RDMMCff(M,NGL,N)=RXMMCff(M,NGL,N)*(1.0_r8-RCCC)
      RDMMNff(M,NGL,N)=RXMMNff(M,NGL,N)*(1.0_r8-RCCN)*(1.0_r8-RCCC)
      RDMMPff(M,NGL,N)=RXMMPff(M,NGL,N)*(1.0_r8-RCCP)*(1.0_r8-RCCC)
      R3MMCff(M,NGL,N)=RXMMCff(M,NGL,N)-RDMMCff(M,NGL,N)
      R3MMNff(M,NGL,N)=RXMMNff(M,NGL,N)-RDMMNff(M,NGL,N)
      R3MMPff(M,NGL,N)=RXMMPff(M,NGL,N)-RDMMPff(M,NGL,N)
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!     EHUM=humus transfer fraction
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
      RHMMCff(M,NGL,N)=AMAX1(0.0_r8,RDMMCff(M,NGL,N)*EHUM)
      RHMMNff(M,NGL,N)=AMAX1(0.0_r8,RDMMNff(M,NGL,N)*EHUM)
      RHMMPff(M,NGL,N)=AMAX1(0.0_r8,RDMMPff(M,NGL,N)*EHUM)
      RCMMCff(M,NGL,N)=RDMMCff(M,NGL,N)-RHMMCff(M,NGL,N)
      RCMMNff(M,NGL,N)=RDMMNff(M,NGL,N)-RHMMNff(M,NGL,N)
      RCMMPff(M,NGL,N)=RDMMPff(M,NGL,N)-RHMMPff(M,NGL,N)

    ENDDO
  ELSE
    DO  M=1,2
      RXMMCff(M,NGL,N)=0.0_r8
      RXMMNff(M,NGL,N)=0.0_r8
      RXMMPff(M,NGL,N)=0.0_r8
      RDMMCff(M,NGL,N)=0.0_r8
      RDMMNff(M,NGL,N)=0.0_r8
      RDMMPff(M,NGL,N)=0.0_r8
      R3MMCff(M,NGL,N)=0.0_r8
      R3MMNff(M,NGL,N)=0.0_r8
      R3MMPff(M,NGL,N)=0.0_r8
      RHMMCff(M,NGL,N)=0.0_r8
      RHMMNff(M,NGL,N)=0.0_r8
      RHMMPff(M,NGL,N)=0.0_r8
      RCMMCff(M,NGL,N)=0.0_r8
      RCMMNff(M,NGL,N)=0.0_r8
      RCMMPff(M,NGL,N)=0.0_r8
    ENDDO
  ENDIF
  end associate
  end subroutine GetMicrobialAnabolismFluxff
!------------------------------------------------------------------------------------------

  subroutine AerobsO2Uptakeff(NGL,N,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
    micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N
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
    WFNff    => nmics%WFNff,      &
    OMAff    => nmics%OMAff,      &
    RUPOXff  => nmicf%RUPOXff,    &
    RGOMOff  => nmicf%RGOMOff,    &
    ROXYMff => nmicf%ROXYMff,     &
    ROXYPff=> nmicf%ROXYPff ,     &
    ROXYOff => nmicf%ROXYOff,     &
    RH2GXff  => nmicf%RH2GXff,    &
    ROQCD  => nmicf%ROQCD,    &
    RCO2Xff  => nmicf%RCO2Xff,    &
    RCH3Xff  => nmicf%RCH3Xff,    &
    RCH4Xff  => nmicf%RCH4Xff,    &
    RVOXA  => nmicf%RVOXA,    &
    RVOXB  => nmicf%RVOXB,    &
    nf_amonia_oxi => micpar%nf_amonia_oxi, &
    nf_nitrite_oxi => micpar%nf_nitrite_oxi, &
    ROXYF  => micfor%ROXYF,  &
    COXYE  => micfor%COXYE  , &
    COXR   => micfor%COXR  , &
    COXQ   => micfor%COXQ  , &
    FLQRI  => micfor%FLQRI  , &
    FLQRQ  => micfor%FLQRQ  , &
    litrm  => micfor%litrm  , &
    OLSGL => micfor%OLSGL , &
    ROXYL => micfor%ROXYL  , &
    VOLX  => micfor%VOLX   , &
    VOLY  => micfor%VOLY  , &
    VOLPM  => micfor%VOLPM , &
    VOLWM  => micfor%VOLWM , &
    FILM  => micfor%FILM , &
    THETPM => micfor%THETPM ,&
    TORT => micfor%TORT , &
    ZERO => micfor%ZERO , &
    ZEROS => micfor%ZEROS , &
    DFGS => micfor%DFGS , &
    SOXYL => micstt%SOXYL  , &
    OXYG => micstt%OXYG , &
    OXYS => micstt%OXYS , &
    COXYG => micstt%COXYG,   &
    ROXSK=> micflx%ROXSK, &
    RVMX4ff => micflx%RVMX4ff,  &
    RVMB4ff => micflx%RVMB4ff, &
    RVMX2ff => micflx%RVMX2ff, &
    RVMB2ff => micflx%RVMB2ff  &
  )

  IF(ROXYPff(NGL,N).GT.ZEROS.AND.FOXYX.GT.ZERO)THEN
    IF(.not.litrm.OR.VOLX.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX=ROXYPff(NGL,N)*XNPG
      ROXYFX=ROXYF*XNPG*FOXYX
      OLSGL1=OLSGL*XNPG
      IF(.not.litrm)THEN
        OXYG1=OXYG*FOXYX
        ROXYLX=ROXYL*XNPG*FOXYX
      ELSE
        OXYG1=COXYG*VOLPM(1)*FOXYX
        ROXYLX=(ROXYL+FLQRQ*COXR &
          +FLQRI*COXQ)*XNPG*FOXYX
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
        !     OF AQUEOUS O2 FROM DISSOLUTION RATE CONSTANT 'DFGS'
        !     CALCULATED IN 'WATSUB'
        !
        !     VOLWM,VOLPM,VOLX=water, air and total volumes
        !     ORAD=microbial radius,FILM=water film thickness
        !     DIFOX=aqueous O2 diffusion, TORT=tortuosity
        !     BIOS=microbial number, OMA=active biomass
        !     SOXYL=O2 solubility, OXKX=Km for O2 uptake
        !     OXYS,COXYS=aqueous O2 amount, concentration
        !     OXYG,COXYG=gaseous O2 amount, concentration
        !     RMPOX,ROXSK=O2 uptake
        !
        THETW1=AMAX1(0.0_r8,safe_adb(VOLWM(M),VOLY))
        RRADO=ORAD*(FILM(M)+ORAD)/FILM(M)
        DIFOX=TORT(M)*OLSGL1*12.57_r8*BIOS*OMAff(NGL,N)*RRADO
        VOLWOX=VOLWM(M)*SOXYL
        VOLPOX=VOLPM(M)
        VOLWPM=VOLWOX+VOLPOX
        DO  MX=1,NPT
          OXYG1=OXYG1+ROXYFX
          OXYS1=OXYS1+ROXYLX
          COXYS1=AMIN1(COXYE*SOXYL &
            ,AMAX1(0.0,safe_adb(OXYS1,(VOLWM(M)*FOXYX))))
          X=DIFOX*COXYS1
          IF(X.GT.ZEROS.AND.OXYS1.GT.ZEROS)THEN
            B=-RUPMX-DIFOX*OXKX-X
            C=X*RUPMX
            RMPOX=(-B-SQRT(B*B-4.0*C))/2.0_r8
          ELSE
            RMPOX=0.0_r8
          ENDIF
          OXYS1=OXYS1-RMPOX
          IF(THETPM(M).GT.THETX.AND.VOLPOX.GT.ZEROS)THEN
            ROXDFQ=DFGS(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
          ELSE
            ROXDFQ=0.0_r8
          ENDIF
          OXYG1=OXYG1-ROXDFQ
          OXYS1=OXYS1+ROXDFQ
          RUPOXff(NGL,N)=RUPOXff(NGL,N)+RMPOX
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
      WFNff(NGL,N)=AMIN1(1.0,AMAX1(0.0,RUPOXff(NGL,N)/ROXYPff(NGL,N)))
      IF(N.EQ.nf_amonia_oxi)THEN
        RVMX4ff(NGL,N)=RVMX4ff(NGL,N)*WFNff(NGL,N)
        RVMB4ff(NGL,N)=RVMB4ff(NGL,N)*WFNff(NGL,N)
      ELSEIF(N.EQ.nf_nitrite_oxi)THEN
        RVMX2ff(NGL,N)=RVMX2ff(NGL,N)*WFNff(NGL,N)
        RVMB2ff(NGL,N)=RVMB2ff(NGL,N)*WFNff(NGL,N)
      ENDIF
    ELSE
      RUPOXff(NGL,N)=ROXYPff(NGL,N)
      WFNff(NGL,N)=1.0_r8
    ENDIF
  ELSE
    RUPOXff(NGL,N)=0.0_r8
    WFNff(NGL,N)=1.0
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RGOMO,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2X,RCH3X,RCH4X,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
  !     ROXYO=O2-limited O2 uptake
  !     RVOXA,RVOXB=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !
  RGOMOff(NGL,N)=RGOMP*WFNff(NGL,N)
  RCO2Xff(NGL,N)=RGOMOff(NGL,N)
  RCH3Xff(NGL,N)=0.0_r8
  RCH4Xff(NGL,N)=0.0_r8
  ROXYOff(NGL,N)=ROXYMff(NGL,N)*WFNff(NGL,N)
  RH2GXff(NGL,N)=0.0_r8
  RVOXA(NGL,N)=RVOXPA*WFNff(NGL,N)
  RVOXB(NGL,N)=RVOXPB*WFNff(NGL,N)
  end associate
  end subroutine AerobsO2Uptakeff

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
  associate(                  &
    FOMNff => nmics%FOMNff,       &
    ROXYMff  => nmicf%ROXYMff,    &
    ROXYOff => nmicf%ROXYOff ,    &
    RDNO3ff => nmicf%RDNO3ff ,    &
    RDNOBff => nmicf%RDNOBff ,    &
    RDNO2ff => nmicf%RDNO2ff ,    &
    RDN2Bff => nmicf%RDN2Bff ,    &
    RDN2Off => nmicf%RDN2Off ,    &
    RGOMDff  => nmicf%RGOMDff,    &
    RGOMYff  => nmicf%RGOMYff,    &
    RVOXA    => nmicf%RVOXA  ,    &
    RVOXB    => nmicf%RVOXB  ,    &
    RNO2Y  => micfor%RNO2Y   ,    &
    VLNO3  => micfor%VLNO3   ,    &
    VLNOB  => micfor%VLNOB   ,    &
    RN2BY  =>  micfor%RN2BY  ,    &
    ZEROS  => micfor%ZEROS   ,    &
    ZEROS2  => micfor%ZEROS2   ,    &
    CNO2B  => micstt%CNO2B   ,    &
    CNO2S  => micstt%CNO2S   ,    &
    ZNO2B  => micstt%ZNO2B   ,    &
    ZNO2S  => micstt%ZNO2S   ,    &
    RVMX2ff => micflx%RVMX2ff,    &
    RVMB2ff  => micflx%RVMB2ff    &
  )
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
  IF(RNO2Y.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RVMX2ff(NGL,N)/RNO2Y)
  ELSE
    FNO2=AMAX1(FMN,FOMNff(NGL,N)*VLNO3)
  ENDIF
  IF(RN2BY.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RVMB2ff(NGL,N)/RN2BY)
  ELSE
    FNB2=AMAX1(FMN,FOMNff(NGL,N)*VLNOB)
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2
!
!     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
!     VMXD4=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM=Km for NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD4S,VMXD4B=substrate-unlimited NO2 reduction in non-band,band
!     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NO2 reduction
!     ECNO=efficiency CO2 conversion to biomass
!     ECHZ=growth respiration efficiency
!     RVOXA,RVOXB=total O2-limited (1)NH4,(2)NO2,(3)CH4 oxidation
!
  FNO2S=VLNO3
  FNO2B=VLNOB
  ROXYD=AMAX1(0.0,ROXYMff(NGL,N)-ROXYOff(NGL,N))
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
  ZNO2SX=ZNO2S+RVOXA(NGL,1)
  ZNO2BX=ZNO2B+RVOXB(NGL,1)
  RDNO2ff(NGL,N)=AMAX1(0.0_r8,AMIN1(VMXD4S,ZNO2SX))
  RDN2Bff(NGL,N)=AMAX1(0.0_r8,AMIN1(VMXD4B,ZNO2BX))
  RDNOT=RDNO2ff(NGL,N)+RDN2Bff(NGL,N)
  RGOMYff(NGL,N)=0.0_r8
  RGOMDff(NGL,N)=RDNOT*ECNO*ENOX
  RDNO3ff(NGL,N)=0.0_r8
  RDNOBff(NGL,N)=0.0_r8
  RDN2Off(NGL,N)=0.0_r8
  RVMX2ff(NGL,N)=VMXD4S
  RVMB2ff(NGL,N)=VMXD4B
  RVOXA(NGL,N)=RVOXA(NGL,N)+0.333*RDNO2ff(NGL,N)
  RVOXB(NGL,N)=RVOXB(NGL,N)+0.333*RDN2Bff(NGL,N)
!     TRN2ON=TRN2ON+RDNO2ff(NGL,N)+RDN2Bff(NGL,N)
  end associate
  end subroutine AutotrophDenitrificCatabolism
!------------------------------------------------------------------------------------------
  subroutine NH3OxidizerCatabolism(NGL,N,XCO2,VOLWZ,TFNX,ECHZ,RGOMP,RVOXP,RVOXPA,&
    RVOXPB,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N
  REAL(R8), INTENT(IN) :: VOLWZ,TFNX
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
    FOMAff => nmics%FOMAff,  &
    OMAff  => nmics%OMAff ,  &
    ROXYMff=> nmicf%ROXYMff, &
    ROXYPff=> nmicf%ROXYPff, &
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
    RVMX4ff=> micflx%RVMX4ff, &
    RVMB4ff=> micflx%RVMB4ff,  &
    ROXYSff => micflx%ROXYSff   &
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
    FNH4=AMAX1(FMN,RVMX4ff(NGL,N)/RNH4Y)
  ELSE
    FNH4=AMAX1(FMN,VLNH4*FOMAff(NGL,N))
  ENDIF
  IF(RNHBY.GT.ZEROS)THEN
    FNB4=AMAX1(FMN,RVMB4ff(NGL,N)/RNHBY)
  ELSE
    FNB4=AMAX1(FMN,VLNHB*FOMAff(NGL,N))
  ENDIF
  naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4
  naqfdiag%TFNH4B=naqfdiag%TFNH4B+FNB4
!
!     NITRIFICATION INHIBITION
!
!     ZNFN0=inhibition when fertilizer added
!     ZNFNI=reduction in inhibition since fertilizer added
!     CNH4S,CNH4B=NH4 concentrations in non-band, band
!     TFNX=temperature effect
!     RNFNI=rate constant for inhibition decline
!     ZHKI=inhibition from high CNH4
!     ZNFN4S,ZNFN4B=inhibition in non-band, band
!
  IF(ZNFN0.GT.ZEROS)THEN
    ZNFNI=ZNFNI*(1.0_r8-RNFNI*TFNX)
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
  VMXX=VMXH*TFNGff(NGL,N)*FCNPff(NGL,N)*XCO2*OMAff(NGL,N)
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
  RNNH4=AMAX1(0.0,AMIN1(VMX4S,FNH4*ZNH4S))*ZNFN4S
  RNNHB=AMAX1(0.0,AMIN1(VMX4B,FNB4*ZNH4B))*ZNFN4B
  RVOXP=RNNH4+RNNHB
  RVOXPA=RNNH4
  RVOXPB=RNNHB
  RGOMP=AMAX1(0.0,RVOXP*ECNH*ECHZ)
  RVMX4ff(NGL,N)=VMX4S
  RVMB4ff(NGL,N)=VMX4B
!
!     O2 DEMAND FROM NH3 OXIDATION
!
!     ROXYM=O2 demand from respiration by nitrifiers
!     ROXYP,ROXYM=O2 demand from respiration + NH3 oxidation
!
  ROXYMff(NGL,N)=2.667_r8*RGOMP
  ROXYPff(NGL,N)=ROXYMff(NGL,N)+3.429_r8*RVOXP
  ROXYSff(NGL,N)=ROXYPff(NGL,N)
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
    TFNGff => nmics%TFNGff,    &
    FCNPff => nmics%FCNPff,    &
    FOMNff => nmics%FOMNff,    &
    OMAff  => nmics%OMAff ,    &
    ROXYMff => nmicf%ROXYMff,  &
    ROXYPff => nmicf%ROXYPff,  &
    VLNH4   => micfor%VLNH4   , &
    VLNHB  => micfor%VLNHB   , &
    VLNO3  => micfor%VLNO3   , &
    VLNOB  => micfor%VLNOB   , &
    ZEROS  => micfor%ZEROS   , &
    RNO2Y  => micfor%RNO2Y   , &
    RN2BY  =>  micfor%RN2BY  , &
    CNO2S  => micstt%CNO2S   , &
    CNO2B  => micstt%CNO2B   , &
    ZNO2S  => micstt%ZNO2S   , &
    ZNO2B  => micstt%ZNO2B   , &
    RVMX2ff=> micflx%RVMX2ff  , &
    RVMB2ff => micflx%RVMB2ff , &
    ROXYSff => micflx%ROXYSff   &
  )
!     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
!     MICROBIAL POPULATIONS
!
!     FNO2=fraction of total biological demand for NO2 in non-band, band
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  IF(RNO2Y.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RVMX2ff(NGL,N)/RNO2Y)
  ELSE
    FNO2=AMAX1(FMN,FOMNff(NGL,N)*VLNO3)
  ENDIF
  IF(RN2BY.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RVMB2ff(NGL,N)/RN2BY)
  ELSE
    FNB2=AMAX1(FMN,FOMNff(NGL,N)*VLNOB)
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
  VMXA=TFNGff(NGL,N)*FCNPff(NGL,N)*XCO2*OMAff(NGL,N)*VMXN
  FCN2S=FNH4S*CNO2S/(CNO2S+ZNKM)
  FCN2B=FNHBS*CNO2B/(CNO2B+ZNKM)
  FSBST=FCN2S+FCN2B
  VMX2S=VMXA*FCN2S
  VMX2B=VMXA*FCN2B
  RNNO2=AMAX1(0.0,AMIN1(VMX2S,FNO2*ZNO2S))
  RNNOB=AMAX1(0.0,AMIN1(VMX2B,FNB2*ZNO2B))
  RVOXP=RNNO2+RNNOB
  RVOXPA=RNNO2
  RVOXPB=RNNOB
  RGOMP=AMAX1(0.0,RVOXP*ECNO*ECHZ)
  RVMX2ff(NGL,N)=VMX2S
  RVMB2ff(NGL,N)=VMX2B
!
!     O2 DEMAND FROM NO2 OXIDATION
!
!     ROXYM=O2 demand from respiration by nitrifiers
!     ROXYP,ROXYM=O2 demand from respiration + NO2 oxidation
!
  ROXYMff(NGL,N)=2.667_r8*RGOMP
  ROXYPff(NGL,N)=ROXYMff(NGL,N)+1.143_r8*RVOXP
  ROXYSff(NGL,N)=ROXYPff(NGL,N)

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
    TFNGff => nmics%TFNGff,       &
    FCNPff => nmics%FCNPff,       &
    OMAff  => nmics%OMAff ,       &
    ROXYMff => nmicf%ROXYMff,     &
    ROXYPff => nmicf%ROXYPff  ,    &
    TKS     => micfor%TKS       , &
    CH2GS   => micstt%CH2GS    , &
    H2GS    => micstt%H2GS     , &
    ROXYSff  => micflx%ROXYSff     &
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

  GH2X=8.3143E-03_r8*TKS*LOG((AMAX1(1.0E-03_r8,CH2GS)/H2KI)**4)
  GH2H=GH2X/12.08_r8
  ECHZ=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AMAX1(0.0_r8,(GCOX+GH2H))/EOMH)))
  VMXA=TFNGff(NGL,N)*FCNPff(NGL,N)*XCO2*OMAff(NGL,N)*VMXC
  H2GSX=H2GS+0.111_r8*naqfdiag%TRH2G
  FSBST=CH2GS/(CH2GS+H2KM)
  RGOMP=AMAX1(0.0,AMIN1(1.5*H2GSX,VMXA*FSBST))
  ROXYMff(NGL,N)=0.0_r8
  ROXYPff(NGL,N)=0.0_r8
  ROXYSff(NGL,N)=0.0_r8
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
    OMAff  => nmics%OMAff ,  &
    ROXYMff=> nmicf%ROXYMff, &
    ROXYPff=> nmicf%ROXYPff, &
    CCH4E  => micfor%CCH4E  , &
    VOLPM  => micfor%VOLPM  , &
    VOLWM  => micfor%VOLWM  , &
    ZEROS2  => micfor%ZEROS2 , &
    ZEROS  => micfor%ZEROS , &
    THETPM  => micfor%THETPM, &
    DFGS   => micfor%DFGS   , &
    litrm  => micfor%litrm  , &
    CCH4G  => micstt%CCH4G  , &
    CH4S   => micstt%CH4S   , &
    SCH4L  => micstt%SCH4L  , &
    RCH4L  => micfor%RCH4L , &
    RCH4F  => micfor%RCH4F , &
    ROXYSff  => micflx%ROXYSff &
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
!     XNPG=1.0/(NPH*NPT)
!     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
!     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
!     VOLPM,VOLWM=air,water-filled porosity
!     SCH4L=CH4 aqueous solubility
!     CCK4=Km for CH4 uptake
!     ECHO=efficiency CO2 conversion to biomass
!     RGOMP1=substrate-limited CH4 oxidation
!     RCHDF=gaseous-aqueous CH4 exchange
!     DFGS=rate constant for gaseous-aqueous exchange
!
  ECHZ=EH4X
  VMXA=TFNGff(NGL,N)*FCNPff(NGL,N)*OMAff(NGL,N)*VMX4
  RCH4L1=RCH4L*XNPG
  RCH4F1=RCH4F*XNPG
  RCH4S1=(naqfdiag%TCH4H+naqfdiag%TCH4A)*XNPG
  IF(litrm)THEN
    CH4G1=CCH4E*VOLPM(1)
  ELSE
    CH4G1=CCH4G*VOLPM(1)
  ENDIF
  CH4S1=CH4S
  VMXA1=VMXA*XNPG
  RVOXP=0.0_r8
  RGOMP=0.0_r8
!
!     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
!     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
!
  DO 320 M=1,NPH
    IF(VOLWM(M).GT.ZEROS2)THEN
      VOLWCH=VOLWM(M)*SCH4L
      VOLWPM=VOLWCH+VOLPM(M)
      DO 325 MM=1,NPT
        CH4G1=CH4G1+RCH4F1
        CH4S1=CH4S1+RCH4L1+RCH4S1
        CCH4S1=AMAX1(0.0,safe_adb(CH4S1,VOLWM(M)))
        FSBST=CCH4S1/(CCH4S1+CCK4)
        RVOXP1=AMIN1(AMAX1(0.0,CH4S1)/(1.0+ECHO*ECHZ),VMXA1*FSBST)
        RGOMP1=RVOXP1*ECHO*ECHZ
        CH4S1=CH4S1-RVOXP1-RGOMP1
        IF(THETPM(M).GT.THETX)THEN
          RCHDF=DFGS(M)*(AMAX1(ZEROS,CH4G1)*VOLWCH &
            -CH4S1*VOLPM(M))/VOLWPM
        ELSE
          RCHDF=0.0_r8
        ENDIF
        CH4G1=CH4G1-RCHDF
        CH4S1=CH4S1+RCHDF
        RVOXP=RVOXP+RVOXP1
        RGOMP=RGOMP+RGOMP1

325   CONTINUE
    ENDIF
320   CONTINUE
  RVOXPA=RVOXP
  RVOXPB=0.0_r8
!
!     O2 DEMAND FROM CH4 OXIDATION
!
!     ROXYM=O2 demand from respiration
!     ROXYP=O2 demand from respiration + CH4 oxidation
!
  ROXYMff(NGL,N)=2.667_r8*RGOMP
  ROXYPff(NGL,N)=ROXYMff(NGL,N)+4.00_r8*RVOXP
  ROXYSff(NGL,N)=ROXYPff(NGL,N)
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

!     begin_execution
  associate(                 &
   TFNGff    => nmics%TFNGff  ,  &
   OMAff     => nmics%OMAff  ,  &
   FNH4XRff  => nmicf%FNH4XRff,  &
   FNO3XRff  => nmicf%FNO3XRff,  &
   FPO4XRff  => nmicf%FPO4XRff,  &
   FP14XRff  => nmicf%FP14XRff,  &
   RINH4ff => nmicf%RINH4ff   ,  &
   RINO3ff  => nmicf%RINO3ff  ,  &
   RIPO4ff  => nmicf%RIPO4ff  ,  &
   RINB4ff => nmicf%RINB4ff   ,  &
   RINB3ff => nmicf%RINB3ff   ,  &
   RIPOBff => nmicf%RIPOBff   ,  &
   RINH4Rff  => nmicf%RINH4Rff,  &
   RINO3Rff  => nmicf%RINO3Rff,  &
   RIPO4Rff  => nmicf%RIPO4Rff,  &
   RIP14ff  => nmicf%RIP14ff,    &
   RIP1Bff  => nmicf%RIP1Bff,    &
   RIP14Rff  => nmicf%RIP14Rff,  &
   CNOMCff  => micpar%CNOMCff   , &
   CPOMCff  => micpar%CPOMCff   , &
   VLNH4    => micfor%VLNH4     , &
   VLNHB    => micfor%VLNHB     , &
   VOLW     => micfor%VOLW      , &
   VLNO3    => micfor%VLNO3     , &
   VLNOB    => micfor%VLNOB     , &
   VLPO4    => micfor%VLPO4     , &
   VLPOB    => micfor%VLPOB     , &
   litrm    => micfor%litrm     , &
   OMCff    => micstt%OMCff     , &
   OMNff    => micstt%OMNff     , &
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
   OMPff    => micstt%OMPff     , &
   RINHOff  => micflx%RINHOff   , &
   RINHBff  => micflx%RINHBff   , &
   RINOOff  => micflx%RINOOff   , &
   RINOBff  => micflx%RINOBff   , &
   TRINH4   => micflx%TRINH4    , &
   RIPOOff  => micflx%RIPOOff   , &
   RIPBOff  => micflx%RIPBOff   , &
   TRIPO4   => micflx%TRIPO4    , &
   RIPO1    => micflx%RIPO1     , &
   RIPB1    => micflx%RIPB1     , &
   RINHOR   => micflx%RINHOR    , &
   RINOOR   => micflx%RINOOR    , &
   RIPOOR   => micflx%RIPOOR    , &
   RIPO1ff  => micflx%RIPO1ff   , &
   RIPB1ff  => micflx%RIPB1ff   , &
   RIPO1Rff   => micflx%RIPO1Rff,  &
   RINHORff  => micflx%RINHORff , &
   RINOORff  => micflx%RINOORff , &
   RIPOORff  => micflx%RIPOORff   &
  )
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
!     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMN=microbial nonstructural C,N
!     CNOMC=maximum microbial N:C ratio
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     RINHX=microbially limited NH4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHO,RINHB=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
!     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  FNH4S=VLNH4
  FNHBS=VLNHB
  RINHP=(OMCff(3,NGL,N)*CNOMCff(3,NGL,N)-OMNff(3,NGL,N))
  IF(RINHP.GT.0.0_r8)THEN
    CNH4X=AMAX1(0.0_r8,CNH4S-Z4MN)
    CNH4Y=AMAX1(0.0_r8,CNH4B-Z4MN)
    RINHX=AMIN1(RINHP,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*Z4MX)
    RINHOff(NGL,N)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RINHBff(NGL,N)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M=Z4MN*VOLW*FNH4S
    ZNHBM=Z4MN*VOLW*FNHBS
    RINH4ff(NGL,N)=AMIN1(FNH4X*AMAX1(0.0_r8,(ZNH4S-ZNH4M)) &
      ,RINHOff(NGL,N))
    RINB4ff(NGL,N)=AMIN1(FNB4X*AMAX1(0.0_r8,(ZNH4B-ZNHBM)) &
      ,RINHBff(NGL,N))
  ELSE
    RINHOff(NGL,N)=0.0_r8
    RINHBff(NGL,N)=0.0_r8
    RINH4ff(NGL,N)=RINHP*FNH4S
    RINB4ff(NGL,N)=RINHP*FNHBS
  ENDIF
  TRINH4=TRINH4+(RINH4ff(NGL,N)+RINB4ff(NGL,N))
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
!     VOLW=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
!     TRINH4=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
  FNO3S=VLNO3
  FNO3B=VLNOB
  RINOP=AMAX1(0.0_r8,RINHP-RINH4ff(NGL,N)-RINB4ff(NGL,N))
  IF(RINOP.GT.0.0)THEN
    CNO3X=AMAX1(0.0,CNO3S-ZOMN)
    CNO3Y=AMAX1(0.0,CNO3B-ZOMN)
    RINOX=AMIN1(RINOP,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*ZOMX)
    RINOOff(NGL,N)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RINOBff(NGL,N)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M=ZOMN*VOLW*FNO3S
    ZNOBM=ZOMN*VOLW*FNO3B
    RINO3ff(NGL,N)=AMIN1(FNO3X*AMAX1(0.0_r8,(ZNO3S-ZNO3M)) &
      ,RINOOff(NGL,N))
    RINB3ff(NGL,N)=AMIN1(FNB3X*AMAX1(0.0_r8,(ZNO3B-ZNOBM)) &
      ,RINOBff(NGL,N))
  ELSE
    RINOOff(NGL,N)=0.0_r8
    RINOBff(NGL,N)=0.0_r8
    RINO3ff(NGL,N)=RINOP*FNO3S
    RINB3ff(NGL,N)=RINOP*FNO3B
  ENDIF
  TRINH4=TRINH4+(RINO3ff(NGL,N)+RINB3ff(NGL,N))
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMP=microbial nonstructural C,P
!     CPOMC=maximum microbial P:C ratio
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
!     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
  FH2PS=VLPO4
  FH2PB=VLPOB
  RIPOP=(OMCff(3,NGL,N)*CPOMCff(3,NGL,N)-OMPff(3,NGL,N))
  IF(RIPOP.GT.0.0)THEN
    CH2PX=AMAX1(0.0,CH2P4-HPMN)
    CH2PY=AMAX1(0.0,CH2P4B-HPMN)
    RIPOX=AMIN1(RIPOP,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*HPMX)
    RIPOOff(NGL,N)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RIPBOff(NGL,N)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM=HPMN*VOLW*FH2PS
    H2PBM=HPMN*VOLW*FH2PB
    RIPO4ff(NGL,N)=AMIN1(FPO4X*AMAX1(0.0,(H2PO4-H2POM)) &
      ,RIPOOff(NGL,N))
    RIPOBff(NGL,N)=AMIN1(FPOBX*AMAX1(0.0,(H2POB-H2PBM)) &
      ,RIPBOff(NGL,N))
  ELSE
    RIPOOff(NGL,N)=0.0_r8
    RIPBOff(NGL,N)=0.0_r8
    RIPO4ff(NGL,N)=RIPOP*FH2PS
    RIPOBff(NGL,N)=RIPOP*FH2PB
  ENDIF
  TRIPO4=TRIPO4+(RIPO4ff(NGL,N)+RIPOBff(NGL,N))
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
!     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
!     TRIPO4=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
  FH1PS=VLPO4
  FH1PB=VLPOB
  RIP1P=0.1_r8*AMAX1(0.0_r8,RIPOP-RIPO4ff(NGL,N)-RIPOBff(NGL,N))
  IF(RIP1P.GT.0.0_r8)THEN
    CH1PX=AMAX1(0.0_r8,CH1P4-HPMN)
    CH1PY=AMAX1(0.0_r8,CH1P4B-HPMN)
    RIP1X=AMIN1(RIP1P,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*HPMX)
    RIPO1ff(NGL,N)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RIPB1ff(NGL,N)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM=HPMN*VOLW*FH1PS
    H1PBM=HPMN*VOLW*FH1PB
    RIP14ff(NGL,N)=AMIN1(FP14X*AMAX1(0.0_r8,(H1PO4-H1POM)) &
      ,RIPO1ff(NGL,N))
    RIP1Bff(NGL,N)=AMIN1(FP1BX*AMAX1(0.0_r8,(H1POB-H1PBM)) &
      ,RIPB1ff(NGL,N))
  ELSE
    RIPO1ff(NGL,N)=0.0_r8
    RIPB1ff(NGL,N)=0.0_r8
    RIP14ff(NGL,N)=RIP1P*FH1PS
    RIP1Bff(NGL,N)=RIP1P*FH1PB
  ENDIF
  TRIPO4=TRIPO4+(RIP14ff(NGL,N)+RIP1Bff(NGL,N))
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
!     RINH4R=substrate-limited NH4 mineraln-immobiln
!     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  IF(litrm)THEN
    RINHPR=RINHP-RINH4ff(NGL,N)-RINO3ff(NGL,N)
    IF(RINHPR.GT.0.0_r8)THEN
      CNH4X=AMAX1(0.0_r8,CNH4S-Z4MN)
      CNH4Y=AMAX1(0.0_r8,CNH4B-Z4MN)
      RINHORff(NGL,N)=AMIN1(RINHPR,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VOLW
      RINH4Rff(NGL,N)=AMIN1(FNH4XRff(NGL,N)*AMAX1(0.0_r8 &
        ,(ZNH4T-ZNH4M)),RINHORff(NGL,N))
    ELSE
      RINHORff(NGL,N)=0.0_r8
      RINH4Rff(NGL,N)=RINHPR
    ENDIF
    TRINH4=TRINH4+RINH4Rff(NGL,N)
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
!     RINO3R=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M=NO3 not available for uptake
!     FNO3XR=fraction of biological NO3 demand
!     RINO3R=substrate-limited NO3 immobiln
!     TRINH4=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
    RINOPR=AMAX1(0.0_r8,RINHPR-RINH4Rff(NGL,N))
    IF(RINOPR.GT.0.0_r8)THEN
      CNO3X=AMAX1(0.0_r8,CNO3S-ZOMN)
      CNO3Y=AMAX1(0.0_r8,CNO3B-ZOMN)
      RINOORff(NGL,N)=AMAX1(RINOPR,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VOLW
      RINO3Rff(NGL,N)=AMIN1(FNO3XRff(NGL,N)*AMAX1(0.0_r8 &
        ,(ZNO3T-ZNO3M)),RINOORff(NGL,N))
    ELSE
      RINOORff(NGL,N)=0.0_r8
      RINO3Rff(NGL,N)=RINOPR
    ENDIF
    TRINH4=TRINH4+RINO3Rff(NGL,N)
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
!     VOLW=water content
!     H2P4M=H2PO4 not available for uptake
!     FPO4XR=fractions of biological H2PO4 demand
!     RIPO4R=substrate-limited H2PO4 mineraln-immobiln
!     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    RIPOPR=RIPOP-RIPO4ff(NGL,N)
    IF(RIPOPR.GT.0.0_r8)THEN
      CH2PX=AMAX1(0.0_r8,CH2P4-HPMN)
      CH2PY=AMAX1(0.0_r8,CH2P4B-HPMN)
      RIPOORff(NGL,N)=AMIN1(RIPOPR,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M=HPMN*VOLW
      RIPO4Rff(NGL,N)=AMIN1(FPO4XRff(NGL,N)*AMAX1(0.0_r8 &
        ,(H2P4T-H2P4M)),RIPOORff(NGL,N))
    ELSE
      RIPOORff(NGL,N)=0.0_r8
      RIPO4Rff(NGL,N)=RIPOPR
    ENDIF
    TRIPO4=TRIPO4+RIPO4Rff(NGL,N)
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
!     VOLW=water content
!     H1P4M=HPO4 not available for uptake
!     FP14XR=fraction of biological HPO4 demand
!     RIP14R=substrate-limited HPO4 minereraln-immobiln
!     TRIPO4=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RIP1PR=0.1*AMAX1(0.0_r8,RIPOPR-RIPO4Rff(NGL,N))
    IF(RIP1PR.GT.0.0_r8)THEN
      CH1PX=AMAX1(0.0_r8,CH1P4-HPMN)
      CH1PY=AMAX1(0.0_r8,CH1P4B-HPMN)
      RIPO1Rff(NGL,N)=AMIN1(RIP1PR,BIOA*OMAff(NGL,N)*TFNGff(NGL,N)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M=HPMN*VOLW
      RIP14Rff(NGL,N)=AMIN1(FP14XRff(NGL,N)*AMAX1(0.0_r8 &
        ,(H1P4T-H1P4M)),RIPO1Rff(NGL,N))
    ELSE
      RIPO1Rff(NGL,N)=0.0_r8
      RIP14Rff(NGL,N)=RIP1PR
    ENDIF
    TRIPO4=TRIPO4+RIP14Rff(NGL,N)
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
!     begin_execution
  associate(                      &
    TFNRff  => nmics%TFNRff,      &
    OMN2ff  => nmics%OMN2ff,      &
    RMOMCff => nmicf%RMOMCff ,    &
    RGOMOff  => nmicf%RGOMOff,    &
    RGN2Fff => nmicf%RGN2Fff,     &
    RN2FXff  => nmicf%RN2FXff,    &
    pH      => micfor%pH    ,    &
    OMNff   => micstt%OMNFF       &

  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TFNR=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  FPH=1.0_r8+AMAX1(0.0_r8,0.25_r8*(6.5_r8-PH))
  RMOMX=RMOM*TFNRff(NGL,N)*FPH
  RMOMCff(1,NGL,N)=OMNff(1,NGL,N)*RMOMX*RMOMK(1)
  RMOMCff(2,NGL,N)=OMN2ff(NGL,N)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMOMT=total maintenance respiration
!     RGOMT=growth respiration
!     RXOMT=senescence respiration
!
  RMOMT=RMOMCff(1,NGL,N)+RMOMCff(2,NGL,N)
  RGOMT=AMAX1(0.0_r8,RGOMOff(NGL,N)-RMOMT)
  RXOMT=AMAX1(0.0_r8,RMOMT-RGOMOff(NGL,N))

!
!     N2 FIXATION: N=(6) AEROBIC, (7) ANAEROBIC
!     FROM GROWTH RESPIRATION, FIXATION ENERGY REQUIREMENT,
!     MICROBIAL N REQUIREMENT IN LABILE (1) AND
!     RESISTANT (2) FRACTIONS
!
!     RGN2P=respiration to meet N2 fixation demand
!     OMC,OMN=microbial nonstructural C,N
!     CNOMC=maximum microbial N:C ratio
!     EN2F=N2 fixation yield per unit nonstructural C
!     RGOMT=growth respiration
!     RGN2F=respiration for N2 fixation
!     CZ2GS=aqueous N2 concentration
!     ZFKM=Km for N2 uptake
!     OMGR*OMC(3,NGL,N,K)=nonstructural C limitation to RGN2F
!     RN2FX=N2 fixation rate
!

  RN2FXff(NGL,N)=0.0_r8
  RGN2Fff(NGL,N)=0.0_r8
  end associate
  end subroutine GatherMicrobialRespirationff

!------------------------------------------------------------------------------------------

  subroutine MicrobialAnabolicUpdateff(micfor,micstt,nmicf)

  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: CGROMC
  integer :: N,M,NGL
  associate(                      &
    CGOMCff    => nmicf%CGOMCff,  &
    CGOMNff  => nmicf%CGOMNff,    &
    CGOMPff  => nmicf%CGOMPff,    &
    CGOMSff  => nmicf%CGOMSff,    &
    CGONSff  => nmicf%CGONSff,    &
    CGOPSff  => nmicf%CGOPSff,    &
    RGN2Fff  => nmicf%RGN2Fff,    &
    RGOMOff  => nmicf%RGOMOff,    &
    RGOMDff  => nmicf%RGOMDff,    &
    RINO3ff  => nmicf%RINO3ff,    &
    RCO2Xff  => nmicf%RCO2Xff,    &
    RIPO4ff  => nmicf%RIPO4ff,    &
    RINB4ff => nmicf%RINB4ff ,    &
    RINB3ff => nmicf%RINB3ff ,    &
    RIPOBff => nmicf%RIPOBff ,    &
    RHOMCff => nmicf%RHOMCff ,    &
    RHOMNff => nmicf%RHOMNff ,    &
    RHOMPff => nmicf%RHOMPff ,    &
    RHMMCff  => nmicf%RHMMCff,    &
    RHMMNff  => nmicf%RHMMNff,    &
    RHMMPff  => nmicf%RHMMPff,    &
    RN2FXff  => nmicf%RN2FXff,    &
    RXOMCff  => nmicf%RXOMCff,    &
    RXOMNff  => nmicf%RXOMNff,    &
    RXOMPff  => nmicf%RXOMPff,    &
    R3OMCff  => nmicf%R3OMCff,    &
    R3OMNff  => nmicf%R3OMNff,    &
    R3OMPff  => nmicf%R3OMPff,    &
    RXMMCff  => nmicf%RXMMCff,    &
    RXMMNff  => nmicf%RXMMNff,    &
    RXMMPff  => nmicf%RXMMPff,    &
    R3MMCff  => nmicf%R3MMCff,    &
    R3MMNff  => nmicf%R3MMNff,    &
    R3MMPff  => nmicf%R3MMPff,    &
    RINH4Rff  => nmicf%RINH4Rff,  &
    RINO3Rff  => nmicf%RINO3Rff,  &
    RIPO4Rff  => nmicf%RIPO4Rff,  &
    RIP14ff   => nmicf%RIP14ff ,  &
    RIP1Bff   => nmicf%RIP1Bff,  &
    RIP14Rff  => nmicf%RIP14Rff,  &
    RINH4ff   => nmicf%RINH4ff,  &
    litrm => micfor%litrm,&
    CFOMC => micfor%CFOMC, &
    CFOMCU => micfor%CFOMCU, &
    OSC  => micstt%OSC, &
    OSN  => micstt%OSN, &
    OSP  => micstt%OSP, &
    OMCff => micstt%OMCff, &
    OMNff => micstt%OMNff, &
    OMPff => micstt%OMPff, &
    OSC14U => micstt%OSC14U, &
    OSN14U=> micstt%OSN14U, &
    OSP14U=> micstt%OSP14U, &
    OSC24U=> micstt%OSC24U, &
    OSN24U=> micstt%OSN24U, &
    OSP24U=> micstt%OSP24U, &
    NFGs=> micpar%NFGs, &
    JG=> micpar%jguilds, &
    is_activef_micb => micpar%is_activef_micb  &
  )
  DO  N=1,NFGs
    IF(is_activef_micb(N))THEN
      DO NGL=1,JG
        DO  M=1,2
          OMCff(M,NGL,N)=OMCff(M,NGL,N)+CGOMSff(M,NGL,N)-RXOMCff(M,NGL,N)-RXMMCff(M,NGL,N)
          OMNff(M,NGL,N)=OMNff(M,NGL,N)+CGONSff(M,NGL,N)-RXOMNff(M,NGL,N)-RXMMNff(M,NGL,N)
          OMPff(M,NGL,N)=OMPff(M,NGL,N)+CGOPSff(M,NGL,N)-RXOMPff(M,NGL,N)-RXMMPff(M,NGL,N)

!     HUMIFICATION PRODUCTS
!
!     CFOMC=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!
          IF(.not.litrm)THEN
            OSC(1,4)=OSC(1,4)+CFOMC(1)*(RHOMCff(M,NGL,N)+RHMMCff(M,NGL,N))
            OSN(1,4)=OSN(1,4)+CFOMC(1)*(RHOMNff(M,NGL,N)+RHMMNff(M,NGL,N))
            OSP(1,4)=OSP(1,4)+CFOMC(1)*(RHOMPff(M,NGL,N)+RHMMPff(M,NGL,N))
            OSC(2,4)=OSC(2,4)+CFOMC(2)*(RHOMCff(M,NGL,N)+RHMMCff(M,NGL,N))
            OSN(2,4)=OSN(2,4)+CFOMC(2)*(RHOMNff(M,NGL,N)+RHMMNff(M,NGL,N))
            OSP(2,4)=OSP(2,4)+CFOMC(2)*(RHOMPff(M,NGL,N)+RHMMPff(M,NGL,N))
          ELSE
            OSC14U=OSC14U+CFOMCU(1)*(RHOMCff(M,NGL,N)+RHMMCff(M,NGL,N))
            OSN14U=OSN14U+CFOMCU(1)*(RHOMNff(M,NGL,N)+RHMMNff(M,NGL,N))
            OSP14U=OSP14U+CFOMCU(1)*(RHOMPff(M,NGL,N)+RHMMPff(M,NGL,N))
            OSC24U=OSC24U+CFOMCU(2)*(RHOMCff(M,NGL,N)+RHMMCff(M,NGL,N))
            OSN24U=OSN24U+CFOMCU(2)*(RHOMNff(M,NGL,N)+RHMMNff(M,NGL,N))
            OSP24U=OSP24U+CFOMCU(2)*(RHOMPff(M,NGL,N)+RHMMPff(M,NGL,N))
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
!     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
!     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
!     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
!     RINH4R,RINO3R =substrate-limited NH4,NO3 mineraln-immobiln
!     RIPO4R,RIP14R=substrate-limited H2PO4,HPO4 mineraln-immobiln
!
        CGROMC=CGOMCff(NGL,N)-RGOMOff(NGL,N)-RGOMDff(NGL,N)-RGN2Fff(NGL,N)
        RCO2Xff(NGL,N)=RCO2Xff(NGL,N)+RGN2Fff(NGL,N)
        DO M=1,2
          OMCff(3,NGL,N)=OMCff(3,NGL,N)-CGOMSff(M,NGL,N)+R3OMCff(M,NGL,N)
          OMNff(3,NGL,N)=OMNff(3,NGL,N)-CGONSff(M,NGL,N)+R3OMNff(M,NGL,N)+R3MMNff(M,NGL,N)
          OMPff(3,NGL,N)=OMPff(3,NGL,N)-CGOPSff(M,NGL,N)+R3OMPff(M,NGL,N)+R3MMPff(M,NGL,N)
          RCO2Xff(NGL,N)=RCO2Xff(NGL,N)+R3MMCff(M,NGL,N)
        ENDDO
        OMCff(3,NGL,N)=OMCff(3,NGL,N)+CGROMC
        OMNff(3,NGL,N)=OMNff(3,NGL,N)+CGOMNff(NGL,N) &
          +RINH4ff(NGL,N)+RINB4ff(NGL,N)+RINO3ff(NGL,N)+RINB3ff(NGL,N)+RN2FXff(NGL,N)
        OMPff(3,NGL,N)=OMPff(3,NGL,N)+CGOMPff(NGL,N) &
          +RIPO4ff(NGL,N)+RIPOBff(NGL,N)+RIP14ff(NGL,N)+RIP1Bff(NGL,N)
        IF(litrm)THEN
          OMNff(3,NGL,N)=OMNff(3,NGL,N)+RINH4Rff(NGL,N)+RINO3Rff(NGL,N)
          OMPff(3,NGL,N)=OMPff(3,NGL,N)+RIPO4Rff(NGL,N)+RIP14Rff(NGL,N)
        ENDIF
      enddo
    ENDIF
  ENDDO
  end associate
  end subroutine MicrobialAnabolicUpdateff


end module MicAutoCPLXMod
