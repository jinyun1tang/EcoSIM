module MicBGCMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils  , only : endrun,destroy
  use MicAutoCPLXMod
  use minimathmod, only : safe_adb
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use NitroDiagTypes
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod, only : micforctype
  use MicBGCPars, only : micpar

  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = __FILE__

  integer :: jcplx,jcplx1,NFGs,JG,jsken,ndbiomcp,nlbiomcp
!
  public :: initNitro1Layer, SoilBGCOneLayer

  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro1Layer
!
! DESCRIPTION:
! initialize single layer microibal bgc model
  implicit none

  jcplx =micpar%jcplx
  jcplx1=micpar%jcplx1
  NFGs  =micpar%NFGs
  JG    =micpar%jguilds
  jsken =micpar%jsken
  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp
  call initNitroPars

  end subroutine initNitro1Layer

!------------------------------------------------------------------------------------------

  subroutine SoilBGCOneLayer(micfor,micstt,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
!local variables
  integer :: LL,K,KL,NGL
  integer :: M,N
  type(NitroMicDiagType) :: nmicdiag
  type(NitroAQMFluxDiagType) :: naqfdiag
  type(NitroMicStateType) :: nmics
  type(NitroMicFluxType) :: nmicf
  type(NitroOMcplxFluxType) :: ncplxf
  type(NitroOMcplxStateType) :: ncplxs

! begin_execution
  call nmicf%Init(JG,jcplx,NFGs)
  call nmics%Init(JG,jcplx,NFGs)
  call ncplxf%Init()
  call ncplxs%Init()
  call naqfdiag%ZeroOut()

  micflx%TRINH4=0._r8;micflx%TRIPO4=0._r8
! write(*,*)'StageBGCEnvironCondition'
  call StageBGCEnvironCondition(micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)
!
! write(*,*)'MicrobialCatabolism'
  call MicrobialCatabolism(micfor,micstt,micflx,nmicdiag,&
    naqfdiag,nmicf,nmics,ncplxf,ncplxs)
!
!write(*,*)'ChemoDenitrification'
  call ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)
!
!     DECOMPOSITION
!
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
  DO  K=0,KL
    DO  N=1,7
      DO NGL=1,JG
        ncplxf%ROQCK(K)=ncplxf%ROQCK(K)+nmicf%ROQCD(NGL,N,K)
      enddo
    ENDDO
  ENDDO
        !
        !write(*,*)'PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C'
  call OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
        !
        !     TRANSFER ALL PRIMING AMONG ALL K
        !
        !     TOQCK=total respiration of DOC+DOA in soil layer
        !     ROQCK=total respiration of DOC+DOA in substrate complex
        !     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
        !     OMC,OMN,OMP=microbial C,N,P
        !
  DO K=0,KL
!
          !write(*,*)'DECOMPOSITION OF ORGANIC SUBSTRATES'
!
    call SolidOMDecomposition(K,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs)
!
          !write(*,*)'DOC ADSORPTION - DESORPTION'
!
    call DOMSorption(K,micfor,micstt,nmicf,ncplxf,ncplxs)

  ENDDO
        !write(*,*)'RedistDecompositionProduct'
  call RedistDecompositionProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)
!
        !write(*,*)'MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION'

  call MicrobialAnabolicUpdate(micfor,micstt,nmicf)
!
        !write(*,*)'MICROBIAL COLONIZATION OF NEW LITTER'
!
  call MicrobialLitterColonization(KL,micfor,micstt,ncplxf,ncplxs)
!
!     AGGREGATE ALL TRANSFORMATIONS CALCULATED ABOVE FOR EACH N,K
!
  call AggregateTransformations(micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
!
  call nmics%destroy()
  call nmicf%destroy()
  call ncplxf%destroy()
  call ncplxs%destroy()
  end subroutine SoilBGCOneLayer
!------------------------------------------------------------------------------------------

  subroutine StageBGCEnvironCondition(micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)
  implicit none
  type(micforctype), intent(in) :: micfor
  integer, intent(out) :: KL
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroOMcplxStateType),intent(inout):: ncplxs
  integer  :: K
  integer  :: M,N,NGL
  real(r8) :: ACTV,ACTVM
  real(r8) :: ORGCL
  real(r8) :: RTK,STK,TKSO
  real(r8) :: TOSC,TOSA,TOHC
  real(r8) :: TSRH
!     begin_execution
  associate(                   &
    CNOMA  => nmics%CNOMA,     &
    CPOMA  => nmics%CPOMA,     &
    OMA    => nmics%OMA  ,     &
    OMC2   => nmics%OMC2 ,     &
    OMN2   => nmics%OMN2 ,     &
    FOM2   => nmics%FOM2 ,     &
    FCN    => nmics%FCN  ,     &
    FCP    => nmics%FCP  ,     &
    FCNP   => nmics%FCNP ,     &
    OSRH  => ncplxs%OSRH,      &
    TOMK  => ncplxs%TOMK,      &
    TONK  => ncplxs%TONK,      &
    TOPK  => ncplxs%TOPK,      &
    FOCA  => ncplxs%FOCA,      &
    FOAA  => ncplxs%FOAA,      &
    CNQ  => ncplxs%CNQ,        &
    CPQ  => ncplxs%CPQ,        &
    COQC    => ncplxs%COQC, &
    COQA    => ncplxs%COQA,  &
    ORCT  => ncplxs%ORCT,      &
    OSCT  => ncplxs%OSCT,      &
    OSAT  => ncplxs%OSAT,      &
    TONX  => ncplxs%TONX,      &
    TOPX  => ncplxs%TOPX,      &
    TORC  =>  nmicdiag%TORC,   &
    TOMA  =>  nmicdiag%TOMA,   &
    TOMN  =>  nmicdiag%TOMN,   &
    XCO2  =>  nmicdiag%XCO2,   &
    TFNX  =>  nmicdiag%TFNX,   &
    TFNY  =>  nmicdiag%TFNY,   &
    THETR => nmicdiag%THETR,   &
    THETZ => nmicdiag%THETZ,   &
    VOLWZ => nmicdiag%VOLWZ,   &
    ZNH4T  =>  nmicdiag%ZNH4T, &
    ZNO3T  =>  nmicdiag%ZNO3T, &
    ZNO2T  =>  nmicdiag%ZNO2T, &
    H2P4T  =>  nmicdiag%H2P4T, &
    H1P4T  =>  nmicdiag%H1P4T, &
    CNOMAff  => nmics%CNOMAff,     &
    CPOMAff  => nmics%CPOMAff,     &
    OMAff    => nmics%OMAff  ,     &
    OMC2ff   => nmics%OMC2ff ,     &
    OMN2ff   => nmics%OMN2ff ,     &
    FOM2ff   => nmics%FOM2ff ,     &
    FCNff    => nmics%FCNff  ,     &
    FCPff    => nmics%FCPff  ,     &
    FCNPff   => nmics%FCNPff ,     &
    CNOMCff  => micpar%CNOMCff,     &
    CPOMCff  => micpar%CPOMCff,    &
    CNOMC  => micpar%CNOMC,     &
    CPOMC  => micpar%CPOMC,    &
    FL       => micpar%FL   ,      &
    k_humus=>micpar%k_humus, &
    k_POM=>micpar%k_POM                     , &
    is_activef_micb=> micpar%is_activef_micb, &
    n_anero_faculb  => micpar%n_anero_faculb, &
    nf_amonia_oxi => micpar%nf_amonia_oxi, &
    litrm    => micfor%litrm  , &
    VOLWRX  => micfor%VOLWRX , &
    VOLW   => micfor%VOLW, &
    VOLW0   => micfor%VOLW0, &
    THETY  => micfor%THETY, &
    VOLR   => micfor%VOLR , &
    VOLY   => micfor%VOLY , &
    POROS  => micfor%POROS, &
    ZEROS => micfor%ZEROS, &
    FC     => micfor%FC   , &
    THETW   => micfor%THETW , &
    TKS    => micfor%TKS, &
    OFFSET  => micfor%OFFSET, &
    VOLWM  => micfor%VOLWM  , &
    ZEROS2  => micfor%ZEROS2 , &
    ZERO   => micfor%ZERO    , &
    CCO2S   => micstt%CCO2S, &
    OSC     => micstt%OSC  ,&
    OSA     => micstt%OSA  ,&
    ORC     => micstt%ORC   , &
    OHC     => micstt%OHC   , &
    OHA     => micstt%OHA   , &
    OMC     => micstt%OMC   , &
    OMN     => micstt%OMN   , &
    OMP     => micstt%OMP   , &
    OQC     => micstt%OQC   , &
    OQN     => micstt%OQN   , &
    OQP     => micstt%OQP   , &
    OQA     => micstt%OQA   , &
    H1PO4 => micstt%H1PO4, &
    H1POB => micstt%H1POB, &
    H2PO4 => micstt%H2PO4, &
    H2POB => micstt%H2POB, &
    ZNH4B => micstt%ZNH4B, &
    ZNH4S => micstt%ZNH4S, &
    ZNO2B => micstt%ZNO2B, &
    ZNO2S => micstt%ZNO2S, &
    ZNO3B => micstt%ZNO3B, &
    ZNO3S => micstt%ZNO3S, &
    OMCff => micstt%OMCff, &
    OMNff => micstt%OMNff, &
    OMPff => micstt%OMPff, &
    FOSRH   => micstt%FOSRH   &
  )

! get KL, the number of mic-om complexes

!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH OFFSET FOR THERMAL ADAPTATION
  IF(litrm)THEN
    ! surface litter layer
    KL=2
    IF(VOLWRX.GT.ZEROS2)THEN
      THETR=VOLW0/VOLR
      THETZ=AMAX1(0.0_r8,THETR-THETY)
      VOLWZ=THETZ*VOLR
    ELSE
      VOLWZ=0.0_r8
    ENDIF
  ELSE
!     non-surface layer
    KL=4
    THETZ=AMAX1(0.0_r8,(AMIN1(AMAX1(0.5_r8*POROS,FC),THETW)-THETY))
    VOLWZ=THETZ*VOLY
  ENDIF


!     TKS=soil temperature
!     OFFSET=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TFNX,TFNY=temperature function for growth,maintenance respiration
!
  TKSO=TKS+OFFSET
  RTK=RGAS*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TFNX=EXP(25.229_r8-62500._r8/RTK)/ACTV
  ACTVM=1+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TFNY=EXP(25.214_r8-62500._r8/RTK)/ACTVM
!
!     OXYI=inhibition of fermenters by O2
!     ORGCL=SOC used to calculate microbial concentration
!
!  OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS+2.5)))
!  ORGCL=AMIN1(1.0E+05*BKVL,ORGC)
!
!     TOTAL MINERAL NH4, NO3 AND PO4
!
!     allocate NH4, NO3, HPO4, H2PO4 to non-band and band fractions
!
  ZNH4T=AMAX1(0.0_r8,ZNH4S)+AMAX1(0.0_r8,ZNH4B)
  ZNO3T=AMAX1(0.0_r8,ZNO3S)+AMAX1(0.0_r8,ZNO3B)
  H1P4T=AMAX1(0.0_r8,H1PO4)+AMAX1(0.0_r8,H1POB)
  H2P4T=AMAX1(0.0_r8,H2PO4)+AMAX1(0.0_r8,H2POB)
  ZNO2T=AMAX1(0.0_r8,ZNO2S)+AMAX1(0.0_r8,ZNO2B)
!
!     CCO2S=aqueous CO2 concentration
!
  XCO2=CCO2S/(CCO2S+CCKM)
!
!     TOTAL SUBSTRATE
!
!     TOSC=total SOC, TOSA=total colonized SOC
!     TORC=total microbial residue, TOHC=total adsorbed C
!     in each K:
!     OSCT=total SOC n each K, OSAT=total colonized SOC
!     ORCT=total microbial residue, OHCT=total adsorbed C
!
  TOSC=0.0_r8
  TOSA=0.0_r8
  TORC=0.0_r8
  TOHC=0.0_r8
!
!     TOTAL SOLID SUBSTRATE
!
  DO  K=0,KL
    OSCT(K)=0.0_r8
    OSAT(K)=0.0_r8

    DO M=1,jsken
      OSCT(K)=OSCT(K)+OSC(M,K)
      OSAT(K)=OSAT(K)+OSA(M,K)
    enddo
    TOSC=TOSC+OSCT(K)
    TOSA=TOSA+OSAT(K)
  enddo
!
!     TOTAL BIORESIDUE
!
  DO  K=0,KL
    ORCT(K)=0.0_r8
    DO  M=1,ndbiomcp
      ORCT(K)=ORCT(K)+ORC(M,K)
    ENDDO
    TORC=TORC+ORCT(K)
!
!     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
!
!     OSRH=total SOC
!
    TOHC=TOHC+OHC(K)+OHA(K)
  enddo

  DO 860 K=0,KL
    OSRH(K)=OSAT(K)+ORCT(K)+OHC(K)+OHA(K)
860 CONTINUE
  TSRH=TOSA+TORC+TOHC
!
!     C:N AND C:P RATIOS OF TOTAL BIOMASS
!     CNOMA,CPOMA=N,P contents of active biomass OMA
!     FCN,FCP=effects of N,P limitations on biomass activity
!
  TOMA=0.0_r8
  TOMN=0.0_r8
  DO 890 K=0,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      IF(K.EQ.jcplx)THEN
! the abstract complex
        DO N=1,NFGs
          IF(is_activef_micb(N))THEN
            DO NGL=1,JG
              IF(OMCff(1,NGL,N).GT.ZEROS)THEN
                CNOMAff(NGL,N)=AMAX1(0.0_r8,OMNff(1,NGL,N)/OMCff(1,NGL,N))
                CPOMAff(NGL,N)=AMAX1(0.0_r8,OMPff(1,NGL,N)/OMCff(1,NGL,N))
              ELSE
                CNOMAff(NGL,N)=CNOMCff(1,NGL,N)
                CPOMAff(NGL,N)=CPOMCff(1,NGL,N)
              ENDIF
              OMAff(NGL,N)=AMAX1(0.0_r8,OMCff(1,NGL,N)/FL(1))
              FCNff(NGL,N)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CNOMAff(NGL,N)/CNOMCff(1,NGL,N))))
              FCPff(NGL,N)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CPOMAff(NGL,N)/CPOMCff(1,NGL,N))))
              FCNPff(NGL,N)=AMIN1(FCNff(NGL,N),FCPff(NGL,N))
            ENDDO
!
!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
            DO NGL=1,JG
              TOMA=TOMA+OMAff(NGL,N)
            ENDDO
            IF(N.EQ.nf_amonia_oxi)THEN
              DO NGL=1,JG
                TOMN=TOMN+OMAff(NGL,N)
              ENDDO
            ENDIF
            DO NGL=1,JG
              OMC2ff(NGL,N)=AMAX1(0.0_r8,AMIN1(OMAff(NGL,N)*FL(2),OMCff(2,NGL,N)))
              IF(OMCff(2,NGL,N).GT.ZEROS)THEN
                FOM2ff(NGL,N)=AMAX1(0.0_r8,OMC2ff(NGL,N)/OMCff(2,NGL,N))
                OMN2ff(NGL,N)=AMAX1(0.0_r8,FOM2ff(NGL,N)*OMNff(2,NGL,N))
              ELSE
                FOM2ff(NGL,N)=0.0_r8
                OMN2ff(NGL,N)=0.0_r8
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ELSE
! the omb complexes
        DO 895 N=1,NFGs
          DO NGL=1,JG
            IF(OMC(1,NGL,N,K).GT.ZEROS)THEN
              CNOMA(NGL,N,K)=AMAX1(0.0_r8,OMN(1,NGL,N,K)/OMC(1,NGL,N,K))
              CPOMA(NGL,N,K)=AMAX1(0.0_r8,OMP(1,NGL,N,K)/OMC(1,NGL,N,K))
            ELSE
              CNOMA(NGL,N,K)=CNOMC(1,NGL,N,K)
              CPOMA(NGL,N,K)=CPOMC(1,NGL,N,K)
            ENDIF
            OMA(NGL,N,K)=AMAX1(0.0_r8,OMC(1,NGL,N,K)/FL(1))
            FCN(NGL,N,K)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CNOMA(NGL,N,K)/CNOMC(1,NGL,N,K))))
            FCP(NGL,N,K)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CPOMA(NGL,N,K)/CPOMC(1,NGL,N,K))))
            FCNP(NGL,N,K)=AMIN1(FCN(NGL,N,K),FCP(NGL,N,K))
          ENDDO

!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
          DO NGL=1,JG
            TOMA=TOMA+OMA(NGL,N,K)
          ENDDO
          IF(N.EQ.n_anero_faculb)THEN
            DO NGL=1,JG
              TOMN=TOMN+OMA(NGL,N,K)
            ENDDO
          ENDIF
          DO NGL=1,JG
            OMC2(NGL,N,K)=AMAX1(0.0_r8,AMIN1(OMA(NGL,N,K)*FL(2),OMC(2,NGL,N,K)))
            IF(OMC(2,NGL,N,K).GT.ZEROS)THEN
              FOM2(NGL,N,K)=AMAX1(0.0_r8,OMC2(NGL,N,K)/OMC(2,NGL,N,K))
              OMN2(NGL,N,K)=AMAX1(0.0_r8,FOM2(NGL,N,K)*OMN(2,NGL,N,K))
            ELSE
              FOM2(NGL,N,K)=0.0_r8
              OMN2(NGL,N,K)=0.0_r8
            ENDIF
          ENDDO
!        ENDIF
895     CONTINUE
      ENDIF
    ENDIF
890 CONTINUE

  DO 690 K=0,KL
    TOMK(K)=0.0_r8
    TONK(K)=0.0_r8
    TOPK(K)=0.0_r8
    TONX(K)=0.0_r8
    TOPX(K)=0.0_r8
    DO 685 N=1,NFGs
      DO NGL=1,JG
        TOMK(K)=TOMK(K)+OMA(NGL,N,K)
        TONK(K)=TONK(K)+OMA(NGL,N,K)*CNOMA(NGL,N,K)
        TOPK(K)=TOPK(K)+OMA(NGL,N,K)*CPOMA(NGL,N,K)
        TONX(K)=TONX(K)+OMA(NGL,N,K)*CNOMC(1,NGL,N,K)   !maximum total N in active micb
        TOPX(K)=TOPX(K)+OMA(NGL,N,K)*CPOMC(1,NGL,N,K)   !maximum total P in active micb
      ENDDO
685 CONTINUE
690 CONTINUE
!
!     FOSRH=fraction of total SOC in each substrate complex K
!
  DO 790 K=0,KL
    IF(TSRH.GT.ZEROS)THEN
      FOSRH(K)=OSRH(K)/TSRH
    ELSE
      FOSRH(K)=1.0_r8
    ENDIF
    !
    !     DOC CONCENTRATIONS
    !
    !     COQC,COQA=aqueous DOC,acetate concentrations
    !     VOLWM=soil water content, FOSRH=fraction of total SOC
    !     occupied by each substrate complex K
    !
    IF(VOLWM(NPH).GT.ZEROS2)THEN
      IF(FOSRH(K).GT.ZERO)THEN
        COQC(K)=AMAX1(0.0_r8,OQC(K)/(VOLWM(NPH)*FOSRH(K)))
        COQA(K)=AMAX1(0.0_r8,OQA(K)/(VOLWM(NPH)*FOSRH(K)))
      ELSE
        COQC(K)=AMAX1(0.0_r8,OQC(K)/VOLWM(NPH))
        COQA(K)=AMAX1(0.0_r8,OQA(K)/VOLWM(NPH))
      ENDIF
    ELSE
      COQC(K)=0.0_r8
      COQA(K)=0.0_r8
    ENDIF
!
!     CNQ,CPQ=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
!
    IF(OQC(K).GT.ZEROS)THEN
      CNQ(K)=AMAX1(0.0_r8,OQN(K)/OQC(K))
      CPQ(K)=AMAX1(0.0_r8,OQP(K)/OQC(K))
    ELSE
      CNQ(K)=0.0_r8
      CPQ(K)=0.0_r8
    ENDIF
    IF(OQC(K).GT.ZEROS.AND.OQA(K).GT.ZEROS)THEN
      FOCA(K)=OQC(K)/(OQC(K)+OQA(K))
      FOAA(K)=1.0_r8-FOCA(K)
    ELSEIF(OQC(K).GT.ZEROS)THEN
      FOCA(K)=1.0_r8
      FOAA(K)=0.0_r8
    ELSE
      FOCA(K)=0.0_r8
      FOAA(K)=1.0_r8
    ENDIF
790 CONTINUE
!
  end associate
  end subroutine StageBGCEnvironCondition
!------------------------------------------------------------------------------------------

  subroutine MicrobialCatabolism(micfor,micstt,micflx,nmicdiag,naqfdiag, &
    nmicf, nmics,ncplxf,ncplxs)
  !
  !  Description:
  !
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout):: ncplxs
  integer :: K,M,N,NGL
  real(r8) :: TOMCNK(2)
  REAL(R8) :: OXKX
! begin_execution
  associate(                   &
    TFNR => nmics%TFNR,        &
    TFNG => nmics%TFNG,        &
    OMA  => nmics%OMA,         &
    TCGOQC  => ncplxf%TCGOQC,  &
    TCGOAC  => ncplxf%TCGOAC,  &
    TCGOMN  => ncplxf%TCGOMN,  &
    TCGOMP  => ncplxf%TCGOMP,  &
    XCO2   =>  nmicdiag%XCO2,   &
    TFNX   =>  nmicdiag%TFNX,   &
    TOMA   =>  nmicdiag%TOMA,   &
    TOMN   =>  nmicdiag%TOMN,   &
    RH2GZ  =>  nmicdiag%RH2GZ, &
    WFNG   =>  nmicdiag%WFNG,   &
    TFNY   =>  nmicdiag%TFNY,   &
    ZNH4T  =>  nmicdiag%ZNH4T, &
    ZNO3T  =>  nmicdiag%ZNO3T, &
    ZNO2T  =>  nmicdiag%ZNO2T, &
    H2P4T  =>  nmicdiag%H2P4T, &
    H1P4T  =>  nmicdiag%H1P4T, &
    VOLWZ  =>  nmicdiag%VOLWZ, &
    TFNGff => nmics%TFNGff   , &
    TFNRff => nmics%TFNRff  ,     &
    OMAff  => nmics%OMAff   ,      &
    k_humus => micpar%k_humus, &
    k_POM => micpar%k_POM, &
    n_aero_fungi => micpar%n_aero_fungi, &
    is_activef_micb=> micpar%is_activef_micb ,&
    PSISM  => micfor%PSISM, &
    litrm => micfor%litrm, &
    H1PO4 => micstt%H1PO4, &
    H1POB => micstt%H1POB, &
    H2PO4 => micstt%H2PO4, &
    H2POB => micstt%H2POB, &
    OMCff  => micstt%OMCff, &
    OMC    => micstt%OMC    &
  )
  RH2GZ=0.0_r8
  DO 760 K=0,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      TCGOQC(K)=0.0_r8
      TCGOAC(K)=0.0_r8
      TCGOMN(K)=0.0_r8
      TCGOMP(K)=0.0_r8
      IF(K.EQ.jcplx)THEN
        DO  N=1,NFGs
          IF(is_activef_micb(N))THEN
            TOMCNK(:)=0.0_r8
            DO NGL=1,JG
              DO M=1,2
                TOMCNK(M)=TOMCNK(M)+OMCff(M,NGL,N)
              ENDDO
            ENDDO
            DO NGL=1,JG
              WFNG=EXP(0.2_r8*PSISM)
              OXKX=OXKA
              TFNGff(NGL,N)=TFNX*WFNG
              TFNRff(NGL,N)=TFNY
              IF(OMAff(NGL,N).GT.0.0_r8)THEN
                call ActiveMicrobesff(NGL,N,VOLWZ,XCO2,TFNX,WFNG,TOMCNK,&
                  OXKX,TOMA,TOMN,RH2GZ,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,&
                  micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ELSE
        DO  N=1,NFGs
          TOMCNK(:)=0.0_r8
          DO NGL=1,JG
            DO M=1,2
              TOMCNK(M)=TOMCNK(M)+OMC(M,NGL,N,K)
            ENDDO
          ENDDO

          DO NGL=1,JG
!           WFNG=water potential (PSISM) effect on microbial respiration
!           OXKX=Km for O2 uptake
!           OXKM=Km for heterotrophic O2 uptake set in starts.f
!           TFNG=combined temp and water stress effect on growth respiration
!           TFNR=temperature effect on maintenance respiration
            IF(N.EQ.n_aero_fungi)THEN
              WFNG=EXP(0.1_r8*PSISM)
            ELSE
              WFNG=EXP(0.2_r8*PSISM)
            ENDIF
            OXKX=OXKM
            TFNG(NGL,N,K)=TFNX*WFNG
            TFNR(NGL,N,K)=TFNY
            IF(OMA(NGL,N,K).GT.0.0_r8)THEN
              call ActiveMicrobes(NGL,N,K,VOLWZ,XCO2,TFNX,WFNG,TOMCNK,&
                OXKX,TOMA,TOMN,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,&
                micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
760 CONTINUE
  end associate
  end subroutine MicrobialCatabolism
!------------------------------------------------------------------------------------------
  subroutine ActiveMicrobes(NGL,N,K,VOLWZ,XCO2,TFNX,WFNG,TOMCNK,OXKX,TOMA,TOMN,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: VOLWZ
  real(r8), intent(in):: OXKX,tomcnk(2),WFNG,TFNX,XCO2
  real(r8), intent(in) :: TOMA,TOMN
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  integer  :: M
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
  associate(                  &
    FOMA => nmics%FOMA,       &
    FOMN => nmics%FOMN,       &
    FOMK => nmics%FOMK,       &
    OMA => nmics%OMA,         &
    RUPOX  => nmicf%RUPOX,    &
    RGOMO  => nmicf%RGOMO,    &
    ROXYM => nmicf%ROXYM ,    &
    ROXYP=> nmicf%ROXYP  ,    &
    ROXYO => nmicf%ROXYO ,    &
    RDNO3 => nmicf%RDNO3 ,    &
    RDNOB => nmicf%RDNOB ,    &
    RDNO2 => nmicf%RDNO2 ,    &
    RDN2B => nmicf%RDN2B ,    &
    RDN2O => nmicf%RDN2O ,    &
    RGOMD  => nmicf%RGOMD,    &
    RH2GX  => nmicf%RH2GX,    &
    RGOMY  => nmicf%RGOMY,    &
    RCO2X  => nmicf%RCO2X,    &
    RCH3X  => nmicf%RCH3X,    &
    RCH4X  => nmicf%RCH4X,    &
    TOMK  => ncplxs%TOMK ,    &
    BKVL  => micfor%BKVL, &
    litrm => micfor%litrm, &
    ORGC=> micfor%ORGC, &
    ZEROS => micfor%ZEROS, &
    VOLX => micfor%VOLX, &
    RNO2Y  => micfor%RNO2Y    &

  )
! FOMA,FOMN=fraction of total active biomass C,N in each N and K

  IF(TOMA.GT.ZEROS)THEN
    FOMA(NGL,N,K)=OMA(NGL,N,K)/TOMA
  ELSE
    FOMA(NGL,N,K)=1.0_r8
  ENDIF
  IF(TOMN.GT.ZEROS)THEN
    FOMN(NGL,N,K)=OMA(NGL,N,K)/TOMN
  ELSE
    FOMN(NGL,N,K)=1.0_r8
  ENDIF
  IF(TOMK(K).GT.ZEROS)THEN
    FOMK(NGL,N,K)=OMA(NGL,N,K)/TOMK(K)
  ELSE
    FOMK(NGL,N,K)=1.0_r8
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
    DO 765 M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
765 CONTINUE
  ELSE
    DO 770 M=1,2
      SPOMK(M)=1.0
      RMOMK(M)=1.0
770 CONTINUE
  ENDIF
!
! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
! write(*,*)'SubstrateCompetitionFactors'
  call SubstrateCompetitionFactors(NGL,N,K,FOXYX,FNH4X,&
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA,&
    micfor,naqfdiag,nmicf,nmics,micflx)


  RGOMP=0.0_r8

!
! HETEROTROPHIC BIOMASS RESPIRATION
          !
!  IF(K.LE.4)THEN
!
!   RESPIRATION BY HETEROTROPHIC AEROBES:
!   N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!   (6)N2 FIXERS
!
  IF(N.LE.3.OR.N.EQ.6)THEN
!     write(*,*)'AerobicHeterotrophCatabolism'
    call AerobicHeterotrophCatabolism(NGL,N,K,TFNX,WFNG,FOQC,FOQA,&
      ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
!     RESPIRATION BY HETEROTROPHIC ANAEROBES:
!     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
!
!     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND
!     ACETATE CONCENTRATION
!
!     GH2F=energy yield of acetotrophic methanogenesis per g C
!     GHAX=H2 effect on energy yield of fermentation
!     GOAX=acetate effect on energy yield of fermentation
!     ECHZ=growth respiration efficiency of fermentation
!
  ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
!     write(*,*)'AnaerobCatabolism'
    call AnaerobCatabolism(NGL,N,K,TFNX,WFNG,FOQC,ECHZ,FGOCP,FGOAP,RGOMP,&
      micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx)
!     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
!
!     GOMX=acetate effect on energy yield
!     ECHZ=growth respiration efficiency of aceto. methanogenesis
!
  ELSEIF(N.EQ.5)THEN
!     write(*,*)'AcetoMethanogenCatabolism'
    call AcetoMethanogenCatabolism(NGL,N,K,TFNX,WFNG,FOQA,ECHZ,&
      FGOCP,FGOAP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  ENDIF
!
!     RESPIRATION RATES BY AUTOTROPHS 'RGOMP' FROM SPECIFIC
!     OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR. N=(1) NH4 OXIDIZERS (2) NO2
!     OXIDIZERS,(3) CH4 OXIDIZERS, (5) H2TROPHIC METHANOGENS
!
!  ENDIF
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
  RUPOX(NGL,N,K)=0.0_r8
  IF(N.LE.3.OR.N.EQ.6)THEN
!  N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!    (6)N2 FIXERS
!   write(*,*)'AerobsO2Uptake'
    call AerobsO2Uptake(NGL,N,K,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      micfor,micstt,nmicf,nmics,micflx)
  ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
    RGOMO(NGL,N,K)=RGOMP
    RCO2X(NGL,N,K)=0.333_r8*RGOMO(NGL,N,K)
    RCH3X(NGL,N,K)=0.667_r8*RGOMO(NGL,N,K)
    RCH4X(NGL,N,K)=0.0_r8
    ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
    RH2GX(NGL,N,K)=0.111*RGOMO(NGL,N,K)
  ELSEIF(N.EQ.5)THEN
    RGOMO(NGL,N,K)=RGOMP
    RCO2X(NGL,N,K)=0.50_r8*RGOMO(NGL,N,K)
    RCH3X(NGL,N,K)=0.0_r8
    RCH4X(NGL,N,K)=0.50_r8*RGOMO(NGL,N,K)
    ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
    RH2GX(NGL,N,K)=0.0_r8
  ENDIF
!
!  write(*,*)'HETEROTROPHIC DENITRIFICATION'
!
  IF(N.EQ.2.AND.ROXYM(NGL,N,K).GT.0.0_r8 .AND.(.not.litrm.OR.VOLX.GT.ZEROS))THEN
    call HeteroDenitrificCatabolism(NGL,N,K,FOQC,RGOCP,&
      VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  ELSE
    RDNO3(NGL,N,K)=0.0_r8
    RDNOB(NGL,N,K)=0.0_r8
    RDNO2(NGL,N,K)=0.0_r8
    RDN2B(NGL,N,K)=0.0_r8
    RDN2O(NGL,N,K)=0.0_r8
    RGOMY(NGL,N,K)=0.0_r8
    RGOMD(NGL,N,K)=0.0_r8
  ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
  call BiomassMineralization(NGL,N,K,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,&
    nmicf,nmics,micflx)
!
  call GatherMicrobialRespiration(NGL,N,K,RMOMK,RGOMT,RXOMT,RMOMT,&
    micfor,micstt,nmicf,nmics)
!
  call GetMicrobialAnabolismFlux(NGL,N,K,ECHZ,FGOCP,&
    FGOAP,RGOMT,RXOMT,RMOMT,spomk,rmomk,micfor,micstt,nmicf,&
    nmics,ncplxf,ncplxs)
  end associate
  end subroutine ActiveMicrobes

!------------------------------------------------------------------------------------------

  subroutine ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroAQMFluxDiagType), INTENT(INOUT):: naqfdiag
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: chy1,CHNO2,CHNOB
  real(r8) :: FNO3S,FNO3B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: VMXC4S,VMXC4B
!     begin_execution
  associate(                   &
    TFNX  =>  nmicdiag%TFNX,   &
    RCNO2  =>  nmicdiag%RCNO2, &
    RCNOB  =>  nmicdiag%RCNOB, &
    RCN2O  =>  nmicdiag%RCN2O, &
    RCN2B  =>  nmicdiag%RCN2B, &
    RCNO3  =>  nmicdiag%RCNO3, &
    RCN3B  =>  nmicdiag%RCN3B, &
    RCOQN  =>  nmicdiag%RCOQN,  &
    RNO2Y  =>  micfor%RNO2Y  , &
    RN2BY  =>  micfor%RN2BY ,  &
    ph => micfor%pH, &
    VOLWM => micfor%VOLWM, &
    ZEROS => micfor%ZEROS, &
    ZERO  => micfor%ZERO, &
    VLNOB => micfor%VLNOB, &
    VLNO3 => micfor%VLNO3, &
    CNO2S => micstt%CNO2S, &
    CNO2B => micstt%CNO2B, &
    ZNO2B => micstt%ZNO2B, &
    ZNO2S => micstt%ZNO2S, &
    RVMXC => micflx%RVMXC, &
    RVMBC => micflx%RVMBC &
  )
!
!     FNO2,FNB2=fraction of total NO2 demand in non-band,band
!     VMXC4S,VMXC4B=substrate-unlimited NO2 reduction in non-band,band
!     CHNO2,CHNOB=nitrous acid concentration in non-band,band
!     VOLWM=soil water content
!     FNO3S,FNO3B=fractions of NO2 in non-band,band
!     TFNX=temperature stress function
!     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
!     RCN2O,RCN2B=N2O production from nitrous acid reduction in non-band,band
!     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
!     RCOQN=DON production from nitrous acid reduction
!     RVMXC,RVMBC=demand for NO2 reduction in non-band,band
!     nitrous acid concn CHNO2
  CHY1=AMAX1(ZERO,10.0**(-(PH-3.0)))
  CHNO2=CNO2S*CHY1/0.5
  CHNOB=CNO2B*CHY1/0.5

  IF(RNO2Y.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RVMXC/RNO2Y)
  ELSE
    FNO2=FMN*VLNO3
  ENDIF
  IF(RN2BY.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RVMBC/RN2BY)
  ELSE
    FNB2=FMN*VLNOB
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
  FNO3S=VLNO3
  FNO3B=VLNOB
  VMXC4S=7.5E-02_r8*CHNO2*VOLWM(NPH)*FNO3S*TFNX
  VMXC4B=7.5E-02_r8*CHNOB*VOLWM(NPH)*FNO3B*TFNX
  RCNO2=AMAX1(0.0_r8,AMIN1(ZNO2S*FNO2,VMXC4S))
  RCNOB=AMAX1(0.0_r8,AMIN1(ZNO2B*FNB2,VMXC4B))
  RCN2O=0.10_r8*RCNO2
  RCN2B=0.10_r8*RCNOB
  RCNO3=0.80_r8*RCNO2
  RCN3B=0.80_r8*RCNOB
  RCOQN=0.10*(RCNO2+RCNOB)
  RVMXC=VMXC4S
  RVMBC=VMXC4B

  end associate
  end subroutine ChemoDenitrification
!------------------------------------------------------------------------------------------

  subroutine OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFLuxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroOMcplxFluxType), intent(inout):: ncplxf
  type(NitroOMcplxStateType),intent(inout):: ncplxs
  integer  :: K,M,N,KK,NGL
  real(r8) :: OSRT
  real(r8) :: XFRK,XFRC
  real(r8) :: XFRN,XFRP,XFRA
  real(r8) :: XFMC,XFMN,XFMP
!     begin_execution
  associate(                 &
    TFNG  => nmics%TFNG,     &
    XOMCZ  => nmicf%XOMCZ,   &
    XOMNZ  => nmicf%XOMNZ,   &
    XOMPZ  => nmicf%XOMPZ,   &
    ROQCK  => ncplxf%ROQCK,  &
    XOQCK  => ncplxf%XOQCK,  &
    XOQCZ  => ncplxf%XOQCZ,  &
    XOQNZ  => ncplxf%XOQNZ,  &
    XOQPZ  => ncplxf%XOQPZ,  &
    XOQAZ  => ncplxf%XOQAZ,  &
    OSRH  => ncplxs%OSRH  ,  &
    TOQCK => micstt%TOQCK, &
    OQC => micstt%OQC, &
    OQN => micstt%OQN, &
    OQP => micstt%OQP, &
    OQA => micstt%OQA, &
    OMC => micstt%OMC, &
    OMN => micstt%OMN, &
    OMP => micstt%OMP, &
    ZEROS => micfor%ZEROS , &
    TFND => micfor%TFND  &
  )
!
!     OSRH=total SOC in each K
!     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
!     between each K and KK, FPRIM=priming transfer rate constant
!     TFND=temperature effect on priming transfers
!     ROQCK,OQC,OQN,OQP=respiration,DOC,DON,DOP
!     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
!
  DO 795 K=0,KL
    IF(K.LE.KL-1)THEN
      DO 800 KK=K+1,KL
        OSRT=OSRH(K)+OSRH(KK)
        IF(OSRH(K).GT.ZEROS.AND.OSRH(KK).GT.ZEROS)THEN
          XFRK=FPRIM*TFND*(ROQCK(K)*OSRH(KK)-ROQCK(KK)*OSRH(K))/OSRT
          XFRC=FPRIM*TFND*(OQC(K)*OSRH(KK)-OQC(KK)*OSRH(K))/OSRT
          XFRN=FPRIM*TFND*(OQN(K)*OSRH(KK)-OQN(KK)*OSRH(K))/OSRT
          XFRP=FPRIM*TFND*(OQP(K)*OSRH(KK)-OQP(KK)*OSRH(K))/OSRT
          XFRA=FPRIM*TFND*(OQA(K)*OSRH(KK)-OQA(KK)*OSRH(K))/OSRT
          IF(ROQCK(K)+XOQCK(K)-XFRK.GT.0.0 &
            .AND.ROQCK(KK)+XOQCK(KK)+XFRK.GT.0.0)THEN
            XOQCK(K)=XOQCK(K)-XFRK
            XOQCK(KK)=XOQCK(KK)+XFRK
          ENDIF
          IF(OQC(K)+XOQCZ(K)-XFRC.GT.0.0 &
            .AND.OQC(KK)+XOQCZ(KK)+XFRC.GT.0.0)THEN
            XOQCZ(K)=XOQCZ(K)-XFRC
            XOQCZ(KK)=XOQCZ(KK)+XFRC

          ENDIF
          IF(OQN(K)+XOQNZ(K)-XFRN.GT.0.0 &
            .AND.OQN(KK)+XOQNZ(KK)+XFRN.GT.0.0)THEN
            XOQNZ(K)=XOQNZ(K)-XFRN
            XOQNZ(KK)=XOQNZ(KK)+XFRN

          ENDIF
          IF(OQP(K)+XOQPZ(K)-XFRP.GT.0.0 &
            .AND.OQP(KK)+XOQPZ(KK)+XFRP.GT.0.0)THEN
            XOQPZ(K)=XOQPZ(K)-XFRP
            XOQPZ(KK)=XOQPZ(KK)+XFRP
          ENDIF
          IF(OQA(K)+XOQAZ(K)-XFRA.GT.0.0 &
            .AND.OQA(KK)+XOQAZ(KK)+XFRA.GT.0.0)THEN
            XOQAZ(K)=XOQAZ(K)-XFRA
            XOQAZ(KK)=XOQAZ(KK)+XFRA
          ENDIF
!
!     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
!
!     XFMC,XFMN,XFMP=transfer of microbial C,N,P
!     between each K and KK, FPRIMM=priming transfer rate constant
!     TFNG=temperature+water effect
!     OMC,OMN,OMP=microbial C,N,P
!     OSRH=total SOC in each K
!     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
!
          DO 850 N=1,7
            DO  M=1,3
              DO NGL=1,JG
                XFMC=FPRIMM*TFNG(NGL,N,K)*(OMC(M,NGL,N,K)*OSRH(KK) &
                  -OMC(M,NGL,N,KK)*OSRH(K))/OSRT
                XFMN=FPRIMM*TFNG(NGL,N,K)*(OMN(M,NGL,N,K)*OSRH(KK) &
                  -OMN(M,NGL,N,KK)*OSRH(K))/OSRT
                XFMP=FPRIMM*TFNG(NGL,N,K)*(OMP(M,NGL,N,K)*OSRH(KK) &
                  -OMP(M,NGL,N,KK)*OSRH(K))/OSRT
                IF(OMC(M,NGL,N,K)+XOMCZ(M,NGL,N,K)-XFMC.GT.0.0 &
                  .AND.OMC(M,NGL,N,KK)+XOMCZ(M,NGL,N,KK)+XFMC.GT.0.0)THEN
                  XOMCZ(M,NGL,N,K)=XOMCZ(M,NGL,N,K)-XFMC
                  XOMCZ(M,NGL,N,KK)=XOMCZ(M,NGL,N,KK)+XFMC
                ENDIF
                IF(OMN(M,NGL,N,K)+XOMNZ(M,NGL,N,K)-XFMN.GT.0.0 &
                  .AND.OMN(M,NGL,N,KK)+XOMNZ(M,NGL,N,KK)+XFMN.GT.0.0)THEN
                  XOMNZ(M,NGL,N,K)=XOMNZ(M,NGL,N,K)-XFMN
                  XOMNZ(M,NGL,N,KK)=XOMNZ(M,NGL,N,KK)+XFMN
                ENDIF
                IF(OMP(M,NGL,N,K)+XOMPZ(M,NGL,N,K)-XFMP.GT.0.0 &
                  .AND.OMP(M,NGL,N,KK)+XOMPZ(M,NGL,N,KK)+XFMP.GT.0.0)THEN
                  XOMPZ(M,NGL,N,K)=XOMPZ(M,NGL,N,K)-XFMP
                  XOMPZ(M,NGL,N,KK)=XOMPZ(M,NGL,N,KK)+XFMP
                ENDIF
              enddo
            enddo
850       CONTINUE
        ENDIF
800   CONTINUE
    ENDIF
795   CONTINUE
!
!     TRANSFER ALL PRIMING AMONG ALL K
!
!     TOQCK=total respiration of DOC+DOA in soil layer
!     ROQCK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
  TOQCK=0.0_r8
  DO 840 K=0,KL
    ROQCK(K)=ROQCK(K)+XOQCK(K)
    TOQCK=TOQCK+ROQCK(K)
    OQC(K)=OQC(K)+XOQCZ(K)
    OQN(K)=OQN(K)+XOQNZ(K)
    OQP(K)=OQP(K)+XOQPZ(K)
    OQA(K)=OQA(K)+XOQAZ(K)
    DO  N=1,7
      DO  M=1,3
        do NGL=1,JG
          OMC(M,NGL,N,K)=OMC(M,NGL,N,K)+XOMCZ(M,NGL,N,K)
          OMN(M,NGL,N,K)=OMN(M,NGL,N,K)+XOMNZ(M,NGL,N,K)
          OMP(M,NGL,N,K)=OMP(M,NGL,N,K)+XOMPZ(M,NGL,N,K)

        enddo
      enddo
    enddo
840   CONTINUE
  end associate
  end subroutine OMTransferForPriming
!------------------------------------------------------------------------------------------

  subroutine DOMSorption(K,micfor,micstt,nmicf,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: K
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  real(r8) :: AECX
  real(r8) :: OQCX
  real(r8) :: OQNX,OQPX,OQAX
  real(r8) :: OHCX,OHNX,OHPX,OHAX
  real(r8) :: VOLXX
  real(r8) :: VOLXW,VOLCX,VOLCW,VOLAX,VOLAW
!     begin_execution
  associate(                   &
    CGOMN  => nmicf%CGOMN,     &
    CGOMP  => nmicf%CGOMP,     &
    CGOQC  => nmicf%CGOQC,     &
    CGOAC  => nmicf%CGOAC,     &
    CSORP  => ncplxf%CSORP,    &
    ZSORP  => ncplxf%ZSORP,    &
    PSORP  => ncplxf%PSORP,    &
    CSORPA  => ncplxf%CSORPA,  &
    TCGOQC  => ncplxf%TCGOQC,  &
    TCGOAC  => ncplxf%TCGOAC,  &
    TCGOMN  => ncplxf%TCGOMN,  &
    TCGOMP  => ncplxf%TCGOMP,  &
    OSRH  => ncplxs%OSRH,      &
    FOCA  => ncplxs%FOCA,      &
    FOAA  => ncplxs%FOAA,      &
    BKVL => micfor%BKVL, &
    ZERO => micfor%ZERO, &
    ZEROS2 => micfor%ZEROS2, &
    ZEROS => micfor%ZEROS, &
    litrm => micfor%litrm, &
    VOLWM => micfor%VOLWM, &
    FOSRH => micstt%FOSRH, &
    OQC => micstt%OQC, &
    OQN => micstt%OQN, &
    OQP => micstt%OQP, &
    OQA => micstt%OQA, &
    OHC => micstt%OHC, &
    OHN => micstt%OHN, &
    OHP => micstt%OHP, &
    OHA => micstt%OHA, &
    AEC => micfor%AEC  &
  )
!     VOLWM=soil water content, FOSRH=fraction of total SOC
!     AEC,AECX=anion exchange capacity
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     TCGOQC,TCGOMN,TCGOMP,TCGOAC=total uptake of DOC,DON,DOP,acetate
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     TSORP,HSORP=sorption rate constant and coefficient for OHC
!     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
  IF(VOLWM(NPH).GT.ZEROS2.AND.FOSRH(K).GT.ZERO)THEN
    IF(litrm)THEN
      AECX=0.5E+03_r8
    ELSE
      AECX=AEC
    ENDIF
    OQCX=AMAX1(ZEROS,OQC(K)-TCGOQC(K))
    OQNX=AMAX1(ZEROS,OQN(K)-TCGOMN(K))
    OQPX=AMAX1(ZEROS,OQP(K)-TCGOMP(K))
    OQAX=AMAX1(ZEROS,OQA(K)-TCGOAC(K))
    OHCX=AMAX1(ZEROS,OHC(K))
    OHNX=AMAX1(ZEROS,OHN(K))
    OHPX=AMAX1(ZEROS,OHP(K))
    OHAX=AMAX1(ZEROS,OHA(K))
    VOLXX=BKVL*AECX*HSORP*FOSRH(K)
    VOLXW=VOLWM(NPH)*FOSRH(K)
    IF(FOCA(K).GT.ZERO)THEN
      VOLCX=FOCA(K)*VOLXX
      VOLCW=FOCA(K)*VOLXW
      CSORP(K)=TSORP*(OQCX*VOLCX-OHCX*VOLCW)/(VOLCX+VOLCW)
    ELSE
      CSORP(K)=TSORP*(OQCX*VOLXX-OHCX*VOLXW)/(VOLXX+VOLXW)
    ENDIF
    IF(FOAA(K).GT.ZERO)THEN
      VOLAX=FOAA(K)*VOLXX
      VOLAW=FOAA(K)*VOLXW
      CSORPA(K)=TSORP*(OQAX*VOLAX-OHAX*VOLAW)/(VOLAX+VOLAW)
    ELSE
      CSORPA(K)=TSORP*(OQAX*VOLXX-OHAX*VOLXW)/(VOLXX+VOLXW)
    ENDIF
    ZSORP(K)=TSORP*(OQNX*VOLXX-OHNX*VOLXW)/(VOLXX+VOLXW)
    PSORP(K)=TSORP*(OQPX*VOLXX-OHPX*VOLXW)/(VOLXX+VOLXW)
  ELSE
    CSORP(K)=0.0_r8
    CSORPA(K)=0.0_r8
    ZSORP(K)=0.0_r8
    PSORP(K)=0.0_r8
  ENDIF
  end associate
  end subroutine DOMSorption
!------------------------------------------------------------------------------------------

  subroutine SolidOMDecomposition(K,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs)

  implicit none
  integer, intent(in) :: K
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  integer  :: M
  real(r8) :: CNOMX,CPOMX
  real(r8) :: COQCK,COSC
  real(r8) :: CPR,CNR
  real(r8) :: DCKD
  real(r8) :: DFNS
  real(r8) :: OQCI
  real(r8) :: RHOSCM
  real(r8) :: FCNK(0:jcplx1),FCPK(0:jcplx1)
  real(r8) :: CNS(4,0:jcplx1),CPS(4,0:jcplx1)

!     begin_execution
  associate(                 &
    RDOSC  => ncplxf%RDOSC,  &
    RDOSN  => ncplxf%RDOSN,  &
    RDOSP  => ncplxf%RDOSP,  &
    RHOSC  => ncplxf%RHOSC,  &
    RHOSN  => ncplxf%RHOSN,  &
    RHOSP  => ncplxf%RHOSP,  &
    RCOSC  => ncplxf%RCOSC,  &
    RCOSN  => ncplxf%RCOSN,  &
    RCOSP  => ncplxf%RCOSP,  &
    RDORC  => ncplxf%RDORC,  &
    RDORN  => ncplxf%RDORN,  &
    RDORP  => ncplxf%RDORP,  &
    RDOHC  => ncplxf%RDOHC,  &
    RDOHN  => ncplxf%RDOHN,  &
    RDOHP  => ncplxf%RDOHP,  &
    RDOHA  => ncplxf%RDOHA,  &
    ROQCK  => ncplxf%ROQCK,  &
    OSRH  => ncplxs%OSRH,    &
    TOMK  => ncplxs%TOMK,    &
    TONK  => ncplxs%TONK,    &
    TOPK  => ncplxs%TOPK,    &
    CNH  => ncplxs%CNH,      &
    CPH  => ncplxs%CPH,      &
    TONX  => ncplxs%TONX,    &
    TOPX  => ncplxs%TOPX,    &
    VOLWZ =>  nmicdiag%VOLWZ,&
    TFNX  =>  nmicdiag%TFNX ,&
    SPOSC => micpar%SPOSC   ,&
    CNRH  => micpar%CNRH    ,&
    CPRH  => micpar%CPRH    ,&
    COQC  => ncplxs%COQC  ,&
    EPOC => micstt%EPOC, &
    CNOSC => micstt%CNOSC, &
    CPOSC => micstt%CPOSC, &
    OHC => micstt%OHC, &
    OHN => micstt%OHN, &
    OHP => micstt%OHP, &
    OSC => micstt%OSC, &
    OSN => micstt%OSN, &
    OSP => micstt%OSP, &
    OSA => micstt%OSA, &
    ORC => micstt%ORC, &
    ORN => micstt%ORN, &
    ORP => micstt%ORP, &
    OHA => micstt%OHA, &
    VOLY => micfor%VOLY, &
    ZEROS => micfor%ZEROS, &
    ZEROS2 => micfor%ZEROS2, &
    litrm => micfor%litrm, &
    BKVL  => micfor%BKVL  &
  )
!     FCPK=N,P limitation to microbial activity in each K
!     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
!     COQCK=aqueous concentration of microbial activity
!     DCKD=Km for decomposition of SOC at current COQCK
!     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
!     DCKI=inhibition of decomposition by microbial concentration
!     OSRH=total SOC
!     COSC=concentration of total SOC
!     BKVL,VOLX=mass, volume of soil layer
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     OQKI=DOC product inhibition constant for decomposition
!
  IF(TOMK(K).GT.ZEROS)THEN
    CNOMX=TONK(K)/TONX(K)
    CPOMX=TOPK(K)/TOPX(K)
    FCNK(K)=AMIN1(1.0,AMAX1(0.50,CNOMX))
    FCPK(K)=AMIN1(1.0,AMAX1(0.50,CPOMX))
  ELSE
    FCNK(K)=1.0
    FCPK(K)=1.0
  ENDIF
!
!     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
!     CONSTANT FOR DECOMPOSITION
!
  IF(VOLWZ.GT.ZEROS2)THEN
    COQCK=AMIN1(0.1E+06,ROQCK(K)/VOLWZ)
  ELSE
    COQCK=0.1E+06
  ENDIF
  IF(litrm)THEN
    DCKD=DCKM0*(1.0+COQCK/DCKI)
  ELSE
    DCKD=DCKML*(1.0+COQCK/DCKI)
  ENDIF
  IF(OSRH(K).GT.ZEROS)THEN
    IF(BKVL.GT.ZEROS)THEN
      COSC=OSRH(K)/BKVL
    ELSE
      COSC=OSRH(K)/VOLY
    ENDIF
    DFNS=COSC/(COSC+DCKD)
    OQCI=1.0/(1.0+COQC(K)/OQKI)
!
!     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES 'RDOS*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     CNS,CPS=N:C,P:C ratios of SOC
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OSA,OSN,OSP=active biomass C,N,P
!     SPOSC=specific decomposition rate constant
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
    DO 785 M=1,jsken
      IF(OSC(M,K).GT.ZEROS)THEN
        CNS(M,K)=AMAX1(0.0_r8,OSN(M,K)/OSC(M,K))
        CPS(M,K)=AMAX1(0.0_r8,OSP(M,K)/OSC(M,K))
        RDOSC(M,K)=AMAX1(0.0_r8,AMIN1(0.5*OSA(M,K) &
          ,SPOSC(M,K)*ROQCK(K)*DFNS*OQCI*TFNX*OSA(M,K)/OSRH(K)))
        RDOSN(M,K)=AMAX1(0.0_r8,AMIN1(OSN(M,K),CNS(M,K)*RDOSC(M,K)))/FCNK(K)
        RDOSP(M,K)=AMAX1(0.0_r8,AMIN1(OSP(M,K),CPS(M,K)*RDOSC(M,K)))/FCPK(K)

      ELSE
        CNS(M,K)=CNOSC(M,K)
        CPS(M,K)=CPOSC(M,K)
        RDOSC(M,K)=0.0_r8
        RDOSN(M,K)=0.0_r8
        RDOSP(M,K)=0.0_r8
      ENDIF
785 CONTINUE
!
!     HUMIFICATION OF DECOMPOSED RESIDUE LIGNIN WITH PROTEIN,
!     CH2O AND CELLULOSE 'RHOS*' WITH REMAINDER 'RCOS*' TO DOC,DON,DOP
!
!     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP
!     CNRH,CPRH=N:C,P:C in POC
!     EPOC=fraction of RDOSC allocated to POC from hour1.f
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
    IF(K.LE.2)THEN
      RHOSC(4,K)=AMAX1(0.0_r8,AMIN1(RDOSN(4,K)/CNRH(3) &
        ,RDOSP(4,K)/CPRH(3),EPOC*RDOSC(4,K)))
      RHOSCM=0.10*RHOSC(4,K)
      RHOSC(1,K)=AMAX1(0.0_r8,AMIN1(RDOSC(1,K),RDOSN(1,K)/CNRH(3) &
        ,RDOSP(1,K)/CPRH(3),RHOSCM))
      RHOSC(2,K)=AMAX1(0.0_r8,AMIN1(RDOSC(2,K),RDOSN(2,K)/CNRH(3) &
        ,RDOSP(2,K)/CPRH(3),RHOSCM))
      RHOSC(3,K)=AMAX1(0.0_r8,AMIN1(RDOSC(3,K),RDOSN(3,K)/CNRH(3) &
        ,RDOSP(3,K)/CPRH(3),RHOSCM-RHOSC(2,K)))
      DO 805 M=1,jsken
        RHOSN(M,K)=AMIN1(RDOSN(M,K),RHOSC(M,K)*CNRH(3))
        RHOSP(M,K)=AMIN1(RDOSP(M,K),RHOSC(M,K)*CPRH(3))
        RCOSC(M,K)=RDOSC(M,K)-RHOSC(M,K)
        RCOSN(M,K)=RDOSN(M,K)-RHOSN(M,K)
        RCOSP(M,K)=RDOSP(M,K)-RHOSP(M,K)
805   CONTINUE
    ELSE
      DO 810 M=1,jsken
        RHOSC(M,K)=0.0_r8
        RHOSN(M,K)=0.0_r8
        RHOSP(M,K)=0.0_r8
        RCOSC(M,K)=RDOSC(M,K)
        RCOSN(M,K)=RDOSN(M,K)
        RCOSP(M,K)=RDOSP(M,K)
810   CONTINUE
    ENDIF
  ELSE
    DO 780 M=1,jsken
      RDOSC(M,K)=0.0_r8
      RDOSN(M,K)=0.0_r8
      RDOSP(M,K)=0.0_r8
      RHOSC(M,K)=0.0_r8
      RHOSN(M,K)=0.0_r8
      RHOSP(M,K)=0.0_r8
      RCOSC(M,K)=0.0_r8
      RCOSN(M,K)=0.0_r8
      RCOSP(M,K)=0.0_r8
780 CONTINUE
  ENDIF
!
!     C, N, P DECOMPOSITION RATE OF BIORESIDUE 'RDOR*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     CNR,CPR=N:C,P:C ratios of microbial residue
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     SPORC=specific decomposition rate constant for microbial residue
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
  IF(OSRH(K).GT.ZEROS)THEN
    DO 775 M=1,2
      IF(ORC(M,K).GT.ZEROS)THEN
        CNR=AMAX1(0.0_r8,ORN(M,K)/ORC(M,K))
        CPR=AMAX1(0.0_r8,ORP(M,K)/ORC(M,K))
        RDORC(M,K)=AMAX1(0.0_r8,AMIN1(ORC(M,K) &
          ,SPORC(M)*ROQCK(K)*DFNS*OQCI*TFNX*ORC(M,K)/OSRH(K)))
    !    3*AMIN1(FCNK(K),FCPK(K))
        RDORN(M,K)=AMAX1(0.0_r8,AMIN1(ORN(M,K),CNR*RDORC(M,K)))/FCNK(K)
        RDORP(M,K)=AMAX1(0.0_r8,AMIN1(ORP(M,K),CPR*RDORC(M,K)))/FCPK(K)
      ELSE
        RDORC(M,K)=0.0_r8
        RDORN(M,K)=0.0_r8
        RDORP(M,K)=0.0_r8
      ENDIF
775 CONTINUE
  ELSE
    DO 776 M=1,2
      RDORC(M,K)=0.0_r8
      RDORN(M,K)=0.0_r8
      RDORP(M,K)=0.0_r8
776 CONTINUE
  ENDIF
!
!     C, N, P DECOMPOSITION RATE OF SORBED SUBSTRATES 'RDOH*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     CNH,CPH=N:C,P:C ratios of adsorbed C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     SPOHC=specific decomposition rate constant for adsorbed C
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
  IF(OSRH(K).GT.ZEROS)THEN
    IF(OHC(K).GT.ZEROS)THEN
      CNH(K)=AMAX1(0.0_r8,OHN(K)/OHC(K))
      CPH(K)=AMAX1(0.0_r8,OHP(K)/OHC(K))
      RDOHC(K)=AMAX1(0.0_r8,AMIN1(OHC(K) &
        ,SPOHC*ROQCK(K)*DFNS*OQCI*TFNX*OHC(K)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      RDOHN(K)=AMAX1(0.0_r8,AMIN1(OHN(K),CNH(K)*RDOHC(K)))/FCNK(K)
      RDOHP(K)=AMAX1(0.0_r8,AMIN1(OHP(K),CPH(K)*RDOHC(K)))/FCPK(K)
      RDOHA(K)=AMAX1(0.0_r8,AMIN1(OHA(K) &
        ,SPOHA*ROQCK(K)*DFNS*TFNX*OHA(K)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
    ELSE
      CNH(K)=0.0_r8
      CPH(K)=0.0_r8
      RDOHC(K)=0.0_r8
      RDOHN(K)=0.0_r8
      RDOHP(K)=0.0_r8
      RDOHA(K)=0.0_r8
    ENDIF
  ELSE
    CNH(K)=0.0_r8
    CPH(K)=0.0_r8
    RDOHC(K)=0.0_r8
    RDOHN(K)=0.0_r8
    RDOHP(K)=0.0_r8
    RDOHA(K)=0.0_r8
  ENDIF
  end associate
  end subroutine SolidOMDecomposition
!------------------------------------------------------------------------------------------

  subroutine RedistDecompositionProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicDiagType), intent(in) :: nmicdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout):: ncplxs
  integer :: K,M,N,NGL
  real(r8) :: FORC(0:jcplx)
!     begin_execution
  associate(                   &
    CGOMN  => nmicf%CGOMN,     &
    CGOMP  => nmicf%CGOMP,     &
    CGOQC  => nmicf%CGOQC,     &
    CGOAC  => nmicf%CGOAC,     &
    RCOMC  => nmicf%RCOMC,     &
    RCOMN  => nmicf%RCOMN,     &
    RCOMP  => nmicf%RCOMP,     &
    RCMMC  => nmicf%RCMMC,     &
    RCMMN  => nmicf%RCMMN,     &
    RCMMP  => nmicf%RCMMP,     &
    RCCMC  => nmicf%RCCMC,     &
    RCCMN  => nmicf%RCCMN,     &
    RCCMP  => nmicf%RCCMP,     &
    RCH3X  => nmicf%RCH3X,     &
    RDOSC  => ncplxf%RDOSC,  &
    RDOSN  => ncplxf%RDOSN,  &
    RDOSP  => ncplxf%RDOSP,  &
    RHOSC  => ncplxf%RHOSC,  &
    RHOSN  => ncplxf%RHOSN,  &
    RHOSP  => ncplxf%RHOSP,  &
    RCOSC  => ncplxf%RCOSC,  &
    RCOSN  => ncplxf%RCOSN,  &
    RCOSP  => ncplxf%RCOSP,  &
    RDORC  => ncplxf%RDORC,  &
    RDORN  => ncplxf%RDORN,  &
    RDORP  => ncplxf%RDORP,  &
    RDOHC  => ncplxf%RDOHC,  &
    RDOHN  => ncplxf%RDOHN,  &
    RDOHP  => ncplxf%RDOHP,  &
    RDOHA  => ncplxf%RDOHA,  &
    CSORP  => ncplxf%CSORP,  &
    ZSORP  => ncplxf%ZSORP,  &
    PSORP  => ncplxf%PSORP,  &
    CSORPA  => ncplxf%CSORPA,&
    ORCT  => ncplxs%ORCT    ,&
    RCOQN =>  nmicdiag%RCOQN,&
    TORC  =>  nmicdiag%TORC , &
    RCMMCff  => nmicf%RCMMCff,     &
    RCMMNff  => nmicf%RCMMNff,     &
    RCMMPff  => nmicf%RCMMPff,     &
    RCOMCff  => nmicf%RCOMCff,     &
    RCOMNff  => nmicf%RCOMNff,     &
    RCOMPff  => nmicf%RCOMPff,     &
    OSC      => micstt%OSC   , &
    OSN      => micstt%OSN   , &
    OSP      => micstt%OSP   , &
!    OSA      => micstt%OSA   , &
    OSC13U      => micstt%OSC13U   , &
    OSN13U      => micstt%OSN13U   , &
    OSP13U      => micstt%OSP13U   , &
   OQC  => micstt%OQC, &
   OQN  => micstt%OQN, &
   OQP  => micstt%OQP, &
   OQA  => micstt%OQA, &
   ORC  => micstt%ORC, &
   ORN  => micstt%ORN, &
   ORP  => micstt%ORP, &
   OHC  => micstt%OHC, &
   OHN  => micstt%OHN, &
   OHP  => micstt%OHP, &
   OHA  => micstt%OHA, &
    ZEROS => micfor%ZEROS, &
    Litrm => micfor%litrm  &
  )
!
!     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
!     HETEROTROPHIC SUBSTRATE-MICROBE COMPLEXES
!
!     FORC=fraction of total microbial residue
!     ORCT=microbial residue
!     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
  DO 1690 K=0,KL
    IF(TORC.GT.ZEROS)THEN
      FORC(K)=ORCT(K)/TORC
    ELSE
      IF(K.EQ.3)THEN
        FORC(K)=1.0
      ELSE
        FORC(K)=0.0_r8
      ENDIF
    ENDIF
    DO 1685 N=1,7
      DO 1680 M=1,2
        DO NGL=1,JG
          RCCMC(M,NGL,N,K)=(RCOMCff(M,NGL,N)+RCMMCff(M,NGL,N))*FORC(K)
          RCCMN(M,NGL,N,K)=(RCOMNff(M,NGL,N)+RCMMNff(M,NGL,N))*FORC(K)
          RCCMP(M,NGL,N,K)=(RCOMPff(M,NGL,N)+RCMMPff(M,NGL,N))*FORC(K)
        ENDDO
1680  CONTINUE
1685  CONTINUE
1690  CONTINUE
!
!   REDISTRIBUTE C,N AND P TRANSFORMATIONS AMONG STATE
!   VARIABLES IN SUBSTRATE-MICROBE COMPLEXES
!

  DO 590 K=0,KL
    DO 580 M=1,jsken
!
!     SUBSTRATE DECOMPOSITION PRODUCTS
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OQC,OQN,OQP,OQA=DOC,DON,DOP
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
      OSC(M,K)=OSC(M,K)-RDOSC(M,K)
!     OSA(M,K)=OSA(M,K)-RDOSC(M,K)
      OSN(M,K)=OSN(M,K)-RDOSN(M,K)
      OSP(M,K)=OSP(M,K)-RDOSP(M,K)
      OQC(K)=OQC(K)+RCOSC(M,K)
      OQN(K)=OQN(K)+RCOSN(M,K)
      OQP(K)=OQP(K)+RCOSP(M,K)
!
!     LIGNIFICATION PRODUCTS
!
!       RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!
      IF(.not.litrm)THEN
        OSC(1,3)=OSC(1,3)+RHOSC(M,K)
 !      OSA(1,3)=OSA(1,3)+RHOSC(M,K)
        OSN(1,3)=OSN(1,3)+RHOSN(M,K)
        OSP(1,3)=OSP(1,3)+RHOSP(M,K)
      ELSE
        OSC13U=OSC13U+RHOSC(M,K)
 !      OSA13U=OSA13U+RHOSC(M,K)
        OSN13U=OSN13U+RHOSN(M,K)
        OSP13U=OSP13U+RHOSP(M,K)
      ENDIF

580 CONTINUE
!
!     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     RCOQN=DON production from nitrous acid reduction
!
    DO 575 M=1,2
      ORC(M,K)=ORC(M,K)-RDORC(M,K)
      ORN(M,K)=ORN(M,K)-RDORN(M,K)
      ORP(M,K)=ORP(M,K)-RDORP(M,K)
      OQC(K)=OQC(K)+RDORC(M,K)
      OQN(K)=OQN(K)+RDORN(M,K)
      OQP(K)=OQP(K)+RDORP(M,K)
575 CONTINUE
    OQC(K)=OQC(K)+RDOHC(K)
    OQN(K)=OQN(K)+RDOHN(K)+RCOQN*FORC(K)
    OQP(K)=OQP(K)+RDOHP(K)
    OQA(K)=OQA(K)+RDOHA(K)
    OHC(K)=OHC(K)-RDOHC(K)
    OHN(K)=OHN(K)-RDOHN(K)
    OHP(K)=OHP(K)-RDOHP(K)
    OHA(K)=OHA(K)-RDOHA(K)
!
!     MICROBIAL UPTAKE OF DISSOLVED C, N, P
!
!     CGOQC,CGOAC,CGOMN,CGOMP=DOC,acetate,DON,DOP uptake
!     RCH3X=acetate production from fermentation
!
    DO 570 N=1,7
      DO NGL=1,JG
        OQC(K)=OQC(K)-CGOQC(NGL,N,K)
        OQN(K)=OQN(K)-CGOMN(NGL,N,K)
        OQP(K)=OQP(K)-CGOMP(NGL,N,K)
        OQA(K)=OQA(K)-CGOAC(NGL,N,K)+RCH3X(NGL,N,K)
!
!     MICROBIAL DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
        DO 565 M=1,2
          ORC(M,K)=ORC(M,K)+RCOMC(M,NGL,N,K)+RCCMC(M,NGL,N,K) &
            +RCMMC(M,NGL,N,K)
          ORN(M,K)=ORN(M,K)+RCOMN(M,NGL,N,K)+RCCMN(M,NGL,N,K) &
            +RCMMN(M,NGL,N,K)
          ORP(M,K)=ORP(M,K)+RCOMP(M,NGL,N,K)+RCCMP(M,NGL,N,K) &
            +RCMMP(M,NGL,N,K)
565     CONTINUE
      enddo
570 CONTINUE
!
!     SORPTION PRODUCTS
!
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
    OQC(K)=OQC(K)-CSORP(K)
    OQN(K)=OQN(K)-ZSORP(K)
    OQP(K)=OQP(K)-PSORP(K)
    OQA(K)=OQA(K)-CSORPA(K)
    OHC(K)=OHC(K)+CSORP(K)
    OHN(K)=OHN(K)+ZSORP(K)
    OHP(K)=OHP(K)+PSORP(K)
    OHA(K)=OHA(K)+CSORPA(K)

590 CONTINUE
  end associate
  end subroutine RedistDecompositionProduct
!------------------------------------------------------------------------------------------

  subroutine MicrobialAnabolicUpdate(micfor,micstt,nmicf)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFluxType), intent(inout) :: nmicf
  integer  :: K,M,N,NGL
  real(r8) ::CGROMC
!     begin_execution
  associate(                  &
    CGOMC    => nmicf%CGOMC,  &
    CGOMN  => nmicf%CGOMN,    &
    CGOMP  => nmicf%CGOMP,    &
    CGOMS  => nmicf%CGOMS,    &
    CGONS  => nmicf%CGONS,    &
    CGOPS  => nmicf%CGOPS,    &
    RGN2F  => nmicf%RGN2F,    &
    RGOMO  => nmicf%RGOMO,    &
    RGOMD  => nmicf%RGOMD,    &
    RINO3  => nmicf%RINO3,    &
    RCO2X  => nmicf%RCO2X,    &
    RIPO4  => nmicf%RIPO4,    &
    RINB4 => nmicf%RINB4 ,    &
    RINB3 => nmicf%RINB3 ,    &
    RIPOB => nmicf%RIPOB ,    &
    RHOMC => nmicf%RHOMC ,    &
    RHOMN => nmicf%RHOMN ,    &
    RHOMP => nmicf%RHOMP ,    &
    RHMMC  => nmicf%RHMMC,    &
    RHMMN  => nmicf%RHMMN,    &
    RHMMP  => nmicf%RHMMP,    &
    RN2FX  => nmicf%RN2FX,    &
    RXOMC  => nmicf%RXOMC,    &
    RXOMN  => nmicf%RXOMN,    &
    RXOMP  => nmicf%RXOMP,    &
    R3OMC  => nmicf%R3OMC,    &
    R3OMN  => nmicf%R3OMN,    &
    R3OMP  => nmicf%R3OMP,    &
    RXMMC  => nmicf%RXMMC,    &
    RXMMN  => nmicf%RXMMN,    &
    RXMMP  => nmicf%RXMMP,    &
    R3MMC  => nmicf%R3MMC,    &
    R3MMN  => nmicf%R3MMN,    &
    R3MMP  => nmicf%R3MMP,    &
    RINH4R  => nmicf%RINH4R,  &
    RINO3R  => nmicf%RINO3R,  &
    RIPO4R  => nmicf%RIPO4R,  &
    RIP14   => nmicf%RIP14 ,  &
    RIP1B   => nmicf%RIP1B ,  &
    RIP14R  => nmicf%RIP14R,  &
    RINH4   => nmicf%RINH4,   &
    k_POM => micpar%k_POM, &
    k_humus => micpar%k_humus, &
    OMC => micstt%OMC, &
    OMN => micstt%OMN, &
    OMP => micstt%OMP, &
    OSC => micstt%OSC, &
    OSN => micstt%OSN, &
    OSP => micstt%OSP, &
    OSC14U => micstt%OSC14U, &
    OSN14U=> micstt%OSN14U, &
    OSP14U=> micstt%OSP14U, &
    OSC24U=> micstt%OSC24U, &
    OSN24U=> micstt%OSN24U, &
    OSP24U=> micstt%OSP24U, &
    CFOMC => micfor%CFOMC, &
    CFOMCU => micfor%CFOMCU, &
    Litrm => micfor%litrm  &
  )
!
!     OMC,OMN,OMP=microbial C,N,P
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
!     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
!

  call MicrobialAnabolicUpdateff(micfor,micstt,nmicf)

  DO 550 K=0,jcplx1
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      DO  N=1,NFGs
        DO NGL=1,JG
          DO 540 M=1,2
            OMC(M,NGL,N,K)=OMC(M,NGL,N,K)+CGOMS(M,NGL,N,K) &
              -RXOMC(M,NGL,N,K)-RXMMC(M,NGL,N,K)
            OMN(M,NGL,N,K)=OMN(M,NGL,N,K)+CGONS(M,NGL,N,K) &
              -RXOMN(M,NGL,N,K)-RXMMN(M,NGL,N,K)
            OMP(M,NGL,N,K)=OMP(M,NGL,N,K)+CGOPS(M,NGL,N,K) &
              -RXOMP(M,NGL,N,K)-RXMMP(M,NGL,N,K)
!
!     HUMIFICATION PRODUCTS
!
!     CFOMC=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!
            IF(.not.litrm)THEN
              OSC(1,4)=OSC(1,4)+CFOMC(1)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
              OSN(1,4)=OSN(1,4)+CFOMC(1)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
              OSP(1,4)=OSP(1,4)+CFOMC(1)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
              OSC(2,4)=OSC(2,4)+CFOMC(2)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
              OSN(2,4)=OSN(2,4)+CFOMC(2)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
              OSP(2,4)=OSP(2,4)+CFOMC(2)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
            ELSE
              OSC14U=OSC14U+CFOMCU(1)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
              OSN14U=OSN14U+CFOMCU(1)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
              OSP14U=OSP14U+CFOMCU(1)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
              OSC24U=OSC24U+CFOMC(2)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
              OSN24U=OSN24U+CFOMCU(2)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
              OSP24U=OSP24U+CFOMCU(2)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
            ENDIF
540       CONTINUE

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
          CGROMC=CGOMC(NGL,N,K)-RGOMO(NGL,N,K)-RGOMD(NGL,N,K)-RGN2F(NGL,N,K)
          RCO2X(NGL,N,K)=RCO2X(NGL,N,K)+RGN2F(NGL,N,K)
          DO 555 M=1,2
            OMC(3,NGL,N,K)=OMC(3,NGL,N,K)-CGOMS(M,NGL,N,K)+R3OMC(M,NGL,N,K)
            OMN(3,NGL,N,K)=OMN(3,NGL,N,K)-CGONS(M,NGL,N,K)+R3OMN(M,NGL,N,K)+R3MMN(M,NGL,N,K)
            OMP(3,NGL,N,K)=OMP(3,NGL,N,K)-CGOPS(M,NGL,N,K)+R3OMP(M,NGL,N,K)+R3MMP(M,NGL,N,K)
            RCO2X(NGL,N,K)=RCO2X(NGL,N,K)+R3MMC(M,NGL,N,K)
555       CONTINUE
          OMC(3,NGL,N,K)=OMC(3,NGL,N,K)+CGROMC
          OMN(3,NGL,N,K)=OMN(3,NGL,N,K)+CGOMN(NGL,N,K) &
            +RINH4(NGL,N,K)+RINB4(NGL,N,K)+RINO3(NGL,N,K)+RINB3(NGL,N,K)+RN2FX(NGL,N,K)
          OMP(3,NGL,N,K)=OMP(3,NGL,N,K)+CGOMP(NGL,N,K) &
            +RIPO4(NGL,N,K)+RIPOB(NGL,N,K)+RIP14(NGL,N,K)+RIP1B(NGL,N,K)
          IF(litrm)THEN
            OMN(3,NGL,N,K)=OMN(3,NGL,N,K)+RINH4R(NGL,N,K)+RINO3R(NGL,N,K)
            OMP(3,NGL,N,K)=OMP(3,NGL,N,K)+RIPO4R(NGL,N,K)+RIP14R(NGL,N,K)
          ENDIF
        enddo
      ENDDO
    ENDIF
550 CONTINUE
  end associate
  end subroutine MicrobialAnabolicUpdate
!------------------------------------------------------------------------------------------

  subroutine MicrobialLitterColonization(KL,micfor,micstt,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  integer  :: K,M
  real(r8) :: DOSAK
!     begin_execution
  associate(                 &
    ROQCK  => ncplxf%ROQCK,  &
    OSCT  => ncplxs%OSCT,    &
    OSAT  => ncplxs%OSAT,    &
    ZEROS => micfor%ZEROS, &
    OSA =>micstt%OSA, &
    OSC =>micstt%OSC, &
    DOSA  => micpar%DOSA     &
  )
!     OSCT,OSAT,OSCX=total,colonized,uncolonized SOC
!     OSA,OSC=colonized,total litter
!     DOSA=rate constant for litter colonization
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
  DO 475 K=0,KL
    OSCT(K)=0.0_r8
    OSAT(K)=0.0_r8
    DO  M=1,jsken
      OSCT(K)=OSCT(K)+OSC(M,K)
      OSAT(K)=OSAT(K)+OSA(M,K)
    enddo
475 CONTINUE
  DO 480 K=0,KL
    IF(OSCT(K).GT.ZEROS)THEN
      DOSAK=DOSA(K)*AMAX1(0.0_r8,ROQCK(K))
      DO 485 M=1,jsken
        OSA(M,K)=AMIN1(OSC(M,K) &
          ,OSA(M,K)+DOSAK*OSC(M,K)/OSCT(K))
485   CONTINUE
    ELSE
      DO 490 M=1,jsken
        OSA(M,K)=AMIN1(OSC(M,K),OSA(M,K))
490   CONTINUE
    ENDIF

480 CONTINUE
  end associate
  end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

  subroutine AggregateTransformations(micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicDiagType),intent(in) :: nmicdiag
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(micfluxtype), intent(inout) :: micflx
  integer  :: K,M,N,NGL
!     begin_execution
  associate(                   &
    CGOMC  => nmicf%CGOMC,     &
    CGOMN  => nmicf%CGOMN,     &
    CGOMP  => nmicf%CGOMP,     &
    CGOQC  => nmicf%CGOQC,     &
    CGOAC  => nmicf%CGOAC,     &
    RUPOX  => nmicf%RUPOX,     &
    RDNO3 => nmicf%RDNO3 ,     &
    RDNOB => nmicf%RDNOB ,     &
    RDNO2 => nmicf%RDNO2 ,     &
    RDN2B => nmicf%RDN2B ,     &
    RDN2O => nmicf%RDN2O ,     &
    RGOMD  => nmicf%RGOMD,     &
    RINH4 => nmicf%RINH4 ,     &
    RINO3  => nmicf%RINO3,     &
    RIPO4  => nmicf%RIPO4,     &
    RINB4 => nmicf%RINB4 ,     &
    RINB3 => nmicf%RINB3 ,     &
    RIPOB => nmicf%RIPOB ,     &
    RH2GX  => nmicf%RH2GX,     &
    RINH4R  => nmicf%RINH4R,   &
    RINO3R  => nmicf%RINO3R,   &
    RIPO4R  => nmicf%RIPO4R,   &
    RCO2X  => nmicf%RCO2X  ,   &
    RCH3X  => nmicf%RCH3X ,    &
    RCH4X  => nmicf%RCH4X ,    &
    RVOXA  => nmicf%RVOXA,     &
    RVOXB  => nmicf%RVOXB,     &
    RIP14  => nmicf%RIP14,     &
    RIP1B  => nmicf%RIP1B,     &
    RIP14R  => nmicf%RIP14R,   &
    RN2FX  => nmicf%RN2FX  ,   &
    RCOSC  => ncplxf%RCOSC,    &
    RCOSN  => ncplxf%RCOSN,    &
    RCOSP  => ncplxf%RCOSP,    &
    RDORC  => ncplxf%RDORC,    &
    RDORN  => ncplxf%RDORN,    &
    RDORP  => ncplxf%RDORP,    &
    RDOHC  => ncplxf%RDOHC,    &
    RDOHN  => ncplxf%RDOHN,    &
    RDOHP  => ncplxf%RDOHP,    &
    RDOHA  => ncplxf%RDOHA,    &
    CSORP  => ncplxf%CSORP,    &
    ZSORP  => ncplxf%ZSORP,    &
    PSORP  => ncplxf%PSORP,    &
    CSORPA  => ncplxf%CSORPA,  &
    TFNX  =>  nmicdiag%TFNX,   &
    RH2GZ  =>  nmicdiag%RH2GZ, &
    RCNO2  =>  nmicdiag%RCNO2, &
    RCNOB  =>  nmicdiag%RCNOB, &
    RCN2O  =>  nmicdiag%RCN2O, &
    RCN2B  =>  nmicdiag%RCN2B, &
    RCNO3  =>  nmicdiag%RCNO3, &
    RCN3B  =>  nmicdiag%RCN3B, &
    VOLWZ  =>  nmicdiag%VOLWZ, &
    CGOMCff  => nmicf%CGOMCff,     &
    RH2GXff  => nmicf%RH2GXff,     &
    RDN2Off => nmicf%RDN2Off,    &
    RDN2Bff => nmicf%RDN2Bff ,   &
    RDNO3ff => nmicf%RDNO3ff ,     &
    RCH4Xff  => nmicf%RCH4Xff ,    &
    RDNO2ff => nmicf%RDNO2ff ,     &
    RUPOXff  => nmicf%RUPOXff,     &
    RGOMDff  => nmicf%RGOMDff,     &
    RCO2Xff  => nmicf%RCO2Xff  ,   &
    RIP14Rff  => nmicf%RIP14Rff,   &
    RIPO4Rff  => nmicf%RIPO4Rff,   &
    RINO3Rff  => nmicf%RINO3Rff,   &
    RINH4Rff  => nmicf%RINH4Rff,   &
    RN2FXff  => nmicf%RN2FXff  ,   &
    RIP1Bff  => nmicf%RIP1Bff,     &
    RIPOBff => nmicf%RIPOBff ,     &
    RINB3ff => nmicf%RINB3ff ,     &
    RINB4ff => nmicf%RINB4ff ,     &
    RIPO4ff  => nmicf%RIPO4ff,     &
    RINH4ff => nmicf%RINH4ff ,     &
    RINO3ff  => nmicf%RINO3ff,     &
    RIP14ff  => nmicf%RIP14ff,     &
    RDNOBff => nmicf%RDNOBff ,   &
    TFNQ => micstt%TFNQ, &
    VOLQ => micstt%VOLQ, &
    litrm => micfor%litrm, &
    Lsurf => micfor%Lsurf, &
    k_POM =>micpar%k_POM, &
    k_humus => micpar%k_humus, &
    is_activef_micb => micpar%is_activef_micb, &
    RCH4O => micflx%RCH4O, &
    RCO2O => micflx%RCO2O, &
    RH2GO => micflx%RH2GO, &
    RN2G  => micflx%RN2G, &
    RN2O => micflx%RN2O, &
    RUPOXO => micflx%RUPOXO, &
    XH1BS => micflx%XH1BS, &
    XH1PS => micflx%XH1PS, &
    XH2BS => micflx%XH2BS, &
    XH2PS => micflx%XH2PS, &
    XN2GS => micflx%XN2GS, &
    XNH4B => micflx%XNH4B, &
    XNH4S => micflx%XNH4S, &
    XNO2B => micflx%XNO2B, &
    XNO2S => micflx%XNO2S, &
    XNO3B => micflx%XNO3B, &
    XNO3S => micflx%XNO3S, &
    XOQCS =>micflx%XOQCS     , &
    XOQNS =>micflx%XOQNS     , &
    XOQPS =>micflx%XOQPS     , &
    XOQAS =>micflx%XOQAS       &
  )
  DO 650 K=0,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      IF(K.EQ.jcplx)THEN
        DO  N=1,NFGs
          IF(is_activef_micb(N))THEN
          DO NGL=1,JG
            naqfdiag%TRINH=naqfdiag%TRINH+RINH4ff(NGL,N)
            naqfdiag%TRINO=naqfdiag%TRINO+RINO3ff(NGL,N)
            naqfdiag%TRIPO=naqfdiag%TRIPO+RIPO4ff(NGL,N)
            naqfdiag%TRIP1=naqfdiag%TRIP1+RIP14ff(NGL,N)
            naqfdiag%TRINB=naqfdiag%TRINB+RINB4ff(NGL,N)
            naqfdiag%TRIOB=naqfdiag%TRIOB+RINB3ff(NGL,N)
            naqfdiag%TRIPB=naqfdiag%TRIPB+RIPOBff(NGL,N)
            naqfdiag%TRIB1=naqfdiag%TRIB1+RIP1Bff(NGL,N)
            naqfdiag%TRN2F=naqfdiag%TRN2F+RN2FXff(NGL,N)
            IF(Lsurf)THEN
              naqfdiag%TRINH=naqfdiag%TRINH+RINH4Rff(NGL,N)
              naqfdiag%TRINO=naqfdiag%TRINO+RINO3Rff(NGL,N)
              naqfdiag%TRIPO=naqfdiag%TRIPO+RIPO4Rff(NGL,N)
              naqfdiag%TRIP1=naqfdiag%TRIP1+RIP14Rff(NGL,N)
            ENDIF
            naqfdiag%TRGOM=naqfdiag%TRGOM+RCO2Xff(NGL,N)
            naqfdiag%TRGOC=naqfdiag%TRGOC+RCH4Xff(NGL,N)
            naqfdiag%TRGOD=naqfdiag%TRGOD+RGOMDff(NGL,N)
            naqfdiag%TUPOX=naqfdiag%TUPOX+RUPOXff(NGL,N)
            naqfdiag%TRDN3=naqfdiag%TRDN3+RDNO3ff(NGL,N)
            naqfdiag%TRDNB=naqfdiag%TRDNB+RDNOBff(NGL,N)
            naqfdiag%TRDN2=naqfdiag%TRDN2+RDNO2ff(NGL,N)
            naqfdiag%TRD2B=naqfdiag%TRD2B+RDN2Bff(NGL,N)
            naqfdiag%TRDNO=naqfdiag%TRDNO+RDN2Off(NGL,N)
            naqfdiag%TRGOH=naqfdiag%TRGOH+RH2GXff(NGL,N)
          ENDDO
!          print*,'K',K,naqfdiag%TRGOM
          ENDIF
        ENDDO
      ELSE
        DO N=1,NFGs

          DO NGL=1,JG
            naqfdiag%TRINH=naqfdiag%TRINH+RINH4(NGL,N,K)
            naqfdiag%TRINO=naqfdiag%TRINO+RINO3(NGL,N,K)
            naqfdiag%TRIPO=naqfdiag%TRIPO+RIPO4(NGL,N,K)
            naqfdiag%TRIP1=naqfdiag%TRIP1+RIP14(NGL,N,K)
            naqfdiag%TRINB=naqfdiag%TRINB+RINB4(NGL,N,K)
            naqfdiag%TRIOB=naqfdiag%TRIOB+RINB3(NGL,N,K)
            naqfdiag%TRIPB=naqfdiag%TRIPB+RIPOB(NGL,N,K)
            naqfdiag%TRIB1=naqfdiag%TRIB1+RIP1B(NGL,N,K)
            naqfdiag%TRN2F=naqfdiag%TRN2F+RN2FX(NGL,N,K)
            IF(Lsurf)THEN
              naqfdiag%TRINH=naqfdiag%TRINH+RINH4R(NGL,N,K)
              naqfdiag%TRINO=naqfdiag%TRINO+RINO3R(NGL,N,K)
              naqfdiag%TRIPO=naqfdiag%TRIPO+RIPO4R(NGL,N,K)
              naqfdiag%TRIP1=naqfdiag%TRIP1+RIP14R(NGL,N,K)
            ENDIF
            naqfdiag%TRGOM=naqfdiag%TRGOM+RCO2X(NGL,N,K)
            naqfdiag%TRGOC=naqfdiag%TRGOC+RCH4X(NGL,N,K)
            naqfdiag%TRGOD=naqfdiag%TRGOD+RGOMD(NGL,N,K)
            naqfdiag%TUPOX=naqfdiag%TUPOX+RUPOX(NGL,N,K)
            naqfdiag%TRDN3=naqfdiag%TRDN3+RDNO3(NGL,N,K)
            naqfdiag%TRDNB=naqfdiag%TRDNB+RDNOB(NGL,N,K)
            naqfdiag%TRDN2=naqfdiag%TRDN2+RDNO2(NGL,N,K)
            naqfdiag%TRD2B=naqfdiag%TRD2B+RDN2B(NGL,N,K)
            naqfdiag%TRDNO=naqfdiag%TRDNO+RDN2O(NGL,N,K)
            naqfdiag%TRGOH=naqfdiag%TRGOH+RH2GX(NGL,N,K)
          ENDDO
!          print*,'K',K,naqfdiag%TRGOM
        ENDDO
      ENDIF
    ENDIF
650 CONTINUE
  DO 645 N=1,7
    IF(N.LT.3.OR.N.EQ.5)THEN
!      IF(N.NE.3)THEN
        DO NGL=1,JG
          naqfdiag%TRGOA=naqfdiag%TRGOA+CGOMCff(NGL,N)
        ENDDO
!      ENDIF
    ENDIF
645 CONTINUE
!
!     ALLOCATE AGGREGATED TRANSFORMATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2O=net CO2 uptake
!     TRGOA=total CO2 uptake by autotrophs
!     TRGOM total CO2 emission by heterotrophs reducing O2
!     TRGOD=total CO2 emission by denitrifiers reducing NOx
!     RVOXA(3)=CH4 oxidation
!     RCH4O=net CH4 uptake
!     CGOMC=total CH4 uptake by autotrophs
!     TRGOC=total CH4 emission
!     RH2GO=net H2 uptake
!     RH2GZ,TRGOH=total H2 uptake, emission
!     RUPOXO,TUPOX=total O2 uptake
!     RN2G=total N2 production
!     TRDNO=total N2O reduction
!     RN2O=total N2O uptake
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RCN2O,RCN2B=nitrous acid reduction in non-band,band
!
!  print*,'rco2o',naqfdiag%TRGOA,naqfdiag%TRGOM,naqfdiag%TRGOD
  RCO2O=naqfdiag%TRGOA-naqfdiag%TRGOM-naqfdiag%TRGOD
  RCH4O=-naqfdiag%TRGOC
  DO NGL=1,JG
    RCO2O=RCO2O-RVOXA(NGL,3)
    RCH4O=RCH4O+RVOXA(NGL,3)+CGOMCff(NGL,3)
  ENDDO
  RH2GO =RH2GZ-naqfdiag%TRGOH
  RUPOXO=naqfdiag%TUPOX
  RN2G  =-naqfdiag%TRDNO
  RN2O  =-naqfdiag%TRDN2-naqfdiag%TRD2B-RCN2O-RCN2B+naqfdiag%TRDNO
!
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate
!
  DO 655 K=0,jcplx1
    DO 660 M=1,jsken
      XOQCS(K)=XOQCS(K)+RCOSC(M,K)
      XOQNS(K)=XOQNS(K)+RCOSN(M,K)
      XOQPS(K)=XOQPS(K)+RCOSP(M,K)
660 CONTINUE
    DO 665 M=1,2
      XOQCS(K)=XOQCS(K)+RDORC(M,K)
      XOQNS(K)=XOQNS(K)+RDORN(M,K)
      XOQPS(K)=XOQPS(K)+RDORP(M,K)
665 CONTINUE
    XOQCS(K)=XOQCS(K)+RDOHC(K)
    XOQNS(K)=XOQNS(K)+RDOHN(K)
    XOQPS(K)=XOQPS(K)+RDOHP(K)
    XOQAS(K)=XOQAS(K)+RDOHA(K)
    DO 670 N=1,NFGs
      DO NGL=1,JG
        XOQCS(K)=XOQCS(K)-CGOQC(NGL,N,K)
        XOQNS(K)=XOQNS(K)-CGOMN(NGL,N,K)
        XOQPS(K)=XOQPS(K)-CGOMP(NGL,N,K)
        XOQAS(K)=XOQAS(K)-CGOAC(NGL,N,K)+RCH3X(NGL,N,K)
      ENDDO
670 CONTINUE
    XOQCS(K)=XOQCS(K)-CSORP(K)
    XOQNS(K)=XOQNS(K)-ZSORP(K)
    XOQPS(K)=XOQPS(K)-PSORP(K)
    XOQAS(K)=XOQAS(K)-CSORPA(K)
655 CONTINUE
!
!     XNH4S,XNH4B=net change in NH4 in band,non-band
!     TRINH,TRINB=total NH4 mineraln-immobn in non-band,band
!     RVOXA(1),RVOXB(1)=total NH4 oxidation in non-band,band
!     XNO3S,XNO3B=net change in NO3 in band,non-band
!     TRINO,TRIOB=total NO3 immobn in non-band,band
!     RVOXA(2),RVOXB(2)=total NO2 oxidation in non-band,band
!     TRDN3,TRDNB=total NO3 reduction in non-band,band
!     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
!     XNO2S,XNO2B=net change in NO3 in band,non-band
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
!     XH2PS,XH2BS=net change in H2PO4 in band,non-band
!     TRIPO,TRIPB=total H2PO4 mineraln-immobn in non-band,band
!     XH1PS,XH1BS=net change in HPO4 in band,non-band
!     TRIP1,TRIB1=total HPO4 mineraln-immobn in non-band,band
!     XN2GS=total N2 fixation
!     XZHYS=total H+ production
!     TRN2F=total N2 fixation
!
  XNH4S=-naqfdiag%TRINH
  XNO3S=-naqfdiag%TRINO-naqfdiag%TRDN3+RCNO3
  XNO2S=+naqfdiag%TRDN3-naqfdiag%TRDN2-RCNO2
  XH2PS=-naqfdiag%TRIPO
  XH1PS=-naqfdiag%TRIP1
  XNH4B=-naqfdiag%TRINB
  XNO3B=-naqfdiag%TRIOB-naqfdiag%TRDNB+RCN3B
  XNO2B=naqfdiag%TRDNB-naqfdiag%TRD2B-RCNOB
  DO NGL=1,JG
    XNH4S=XNH4S-RVOXA(NGL,1)
    XNO3S=XNO3S+RVOXA(NGL,2)
    XNO2S=XNO2S+RVOXA(NGL,1)-RVOXA(NGL,2)
    XNH4B=XNH4B-RVOXB(NGL,1)
    XNO3B=XNO3B+RVOXB(NGL,2)
    XNO2B=XNO2B+RVOXB(NGL,1)-RVOXB(NGL,2)
  ENDDO

  XH2BS=-naqfdiag%TRIPB
  XH1BS=-naqfdiag%TRIB1
  XN2GS=naqfdiag%TRN2F
  TFNQ=TFNX
  VOLQ=VOLWZ
  end associate
  end subroutine AggregateTransformations
!------------------------------------------------------------------------------------------

  subroutine SubstrateCompetitionFactors(NGL,N,K,FOXYX,&
    FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA,&
    micfor,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(out):: FOXYX,FNH4X
  real(r8),intent(out) :: FNB3X,FNB4X,FNO3X
  real(r8),intent(out) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8),intent(out) :: FOQC,FOQA
  type(micforctype), intent(in) :: micfor
  type(NitroAQMFluxDiagType),INTENT(INOUT)::  naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                  &
    FOMA    => nmics%FOMA,    &
    FOMK    => nmics%FOMK,    &
    FNH4XR  =>nmicf%FNH4XR,   &
    FNO3XR  =>nmicf%FNO3XR,   &
    FP14XR  => nmicf%FP14XR,  &
    FPO4XR  => nmicf%FPO4XR,  &
    ROXYS => micflx%ROXYS, &
    RINHO => micflx%RINHO, &
    RINHB => micflx%RINHB, &
    RINOO => micflx%RINOO, &
    RINOB => micflx%RINOB, &
    RIPOO => micflx%RIPOO, &
    RIPBO => micflx%RIPBO, &
    RIPO1  => micflx%RIPO1, &
    RIPB1 => micflx%RIPB1, &
    ROQCS  => micflx%ROQCS, &
    ROQAS => micflx%ROQAS, &
    RINHOR => micflx%RINHOR, &
    RINOOR => micflx%RINOOR, &
    RIPOOR => micflx%RIPOOR, &
    RIPO1R => micflx%RIPO1R, &
    litrm  => micfor%litrm, &
    VLNH4 => micfor%VLNH4, &
    VLNHB => micfor%VLNHB, &
    VLNOB  => micfor%VLNOB, &
    VLNO3  => micfor%VLNO3, &
    VLPOB  => micfor%VLPOB, &
    VLPO4  => micfor%VLPO4, &
    ZEROS => micfor%ZEROS, &
    ROXYY   => micfor%ROXYY,  &
    RNH4Y   => micfor%RNH4Y,  &
    RNHBY   => micfor%RNHBY,  &
    RNO3Y   => micfor%RNO3Y,  &
    RNH4YU  => micfor%RNH4YU, &
    RNO3YU  => micfor%RNO3YU, &
    RP14YU  => micfor%RP14YU, &
    RPO4YU  => micfor%RPO4YU, &
    RN3BY     => micfor%RN3BY,    &
    RPO4Y     => micfor%RPO4Y,    &
    RPOBY     => micfor%RPOBY,    &
    RP14Y     => micfor%RP14Y,    &
    RP1BY     => micfor%RP1BY,    &
    ROQCY     => micfor%ROQCY,    &
    ROQAY     => micfor%ROQAY,    &
    Lsurf     => micfor%Lsurf,    &
    BKVL0     => micfor%BKVL0     &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(ROXYY.GT.ZEROS)THEN
    FOXYX=AMAX1(FMN,ROXYS(NGL,N,K)/ROXYY)
  ELSE
    FOXYX=AMAX1(FMN,FOMA(NGL,N,K))
  ENDIF
  IF(RNH4Y.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RINHO(NGL,N,K)/RNH4Y)
  ELSE
    FNH4X=AMAX1(FMN,FOMA(NGL,N,K)*VLNH4)
  ENDIF
  IF(RNHBY.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RINHB(NGL,N,K)/RNHBY)
  ELSE
    FNB4X=AMAX1(FMN,FOMA(NGL,N,K)*VLNHB)
  ENDIF
  IF(RNO3Y.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RINOO(NGL,N,K)/RNO3Y)
  ELSE
    FNO3X=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3)
  ENDIF
  IF(RN3BY.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RINOB(NGL,N,K)/RN3BY)
  ELSE
    FNB3X=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB)
  ENDIF
  IF(RPO4Y.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RIPOO(NGL,N,K)/RPO4Y)
  ELSE
    FPO4X=AMAX1(FMN,FOMA(NGL,N,K)*VLPO4)
  ENDIF
  IF(RPOBY.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RIPBO(NGL,N,K)/RPOBY)
  ELSE
    FPOBX=AMAX1(FMN,FOMA(NGL,N,K)*VLPOB)
  ENDIF
  IF(RP14Y.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RIPO1(NGL,N,K)/RP14Y)
  ELSE
    FP14X=AMAX1(FMN,FOMA(NGL,N,K)*VLPO4)
  ENDIF
  IF(RP1BY.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RIPB1(NGL,N,K)/RP1BY)
  ELSE
    FP1BX=AMAX1(FMN,FOMA(NGL,N,K)*VLPOB)
  ENDIF

  IF(ROQCY(K).GT.ZEROS)THEN
    FOQC=AMAX1(FMN,ROQCS(NGL,N,K)/ROQCY(K))
  ELSE
    FOQC=AMAX1(FMN,FOMK(NGL,N,K))
  ENDIF
  naqfdiag%TFOQC=naqfdiag%TFOQC+FOQC
  IF(ROQAY(K).GT.ZEROS)THEN
    FOQA=AMAX1(FMN,ROQAS(NGL,N,K)/ROQAY(K))
  ELSE
    FOQA=AMAX1(FMN,FOMK(NGL,N,K))
  ENDIF
  naqfdiag%TFOQA=naqfdiag%TFOQA+FOQA
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
      FNH4XR(NGL,N,K)=AMAX1(FMN,RINHOR(NGL,N,K)/RNH4YU)
    ELSE
      FNH4XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    IF(RNO3YU.GT.ZEROS)THEN
      FNO3XR(NGL,N,K)=AMAX1(FMN,RINOOR(NGL,N,K)/RNO3YU)
    ELSE
      FNO3XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    IF(RPO4YU.GT.ZEROS)THEN
      FPO4XR(NGL,N,K)=AMAX1(FMN,RIPOOR(NGL,N,K)/RPO4YU)
    ELSE
      FPO4XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    IF(RP14YU.GT.ZEROS)THEN
      FP14XR(NGL,N,K)=AMAX1(FMN,RIPO1R(NGL,N,K)/RP14YU)
    ELSE
      FP14XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
  ENDIF
  IF(Lsurf.AND.K.NE.micpar%k_POM.AND.K.NE.micpar%k_humus &
    .AND.BKVL0.GT.ZEROS)THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4XR(NGL,N,K)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3XR(NGL,N,K)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4XR(NGL,N,K)
    naqfdiag%TFP14X=naqfdiag%TFP14X+FP14XR(NGL,N,K)
  ENDIF
  end associate
  end subroutine SubstrateCompetitionFactors
!------------------------------------------------------------------------------------------

  subroutine AcetoMethanogenCatabolism(NGL,N,K,TFNX,WFNG,FOQA,ECHZ,&
    FGOCP,FGOAP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: FOQA,TFNX,WFNG
  real(r8), intent(out) :: ECHZ
  REAL(R8), intent(out) :: FGOCP,FGOAP,RGOMP
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  type(NitroOMcplxStateType),intent(in):: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FSBST
  real(r8) :: GOMX,GOMM
  real(r8) :: RGOGY,RGOGZ,RGOGX

! begin_execution
  associate(                  &
    FCNP  => nmics%FCNP,      &
    OMA   => nmics%OMA ,      &
    ROXYP => nmicf%ROXYP,     &
    ROXYM => nmicf%ROXYM,     &
    ROQCD => nmicf%ROQCD,     &
    COQA  => ncplxs%COQA,   &
    OQA => micstt%OQA, &
    ROXYS => micflx%ROXYS, &
    ROQCS => micflx%ROQCS, &
    ROQAS => micflx%ROQAS, &
    ZERO  => micfor%ZERO, &
    TKS  => micfor%TKS &
  )
  GOMX=8.3143E-03*TKS*LOG((AMAX1(ZERO,COQA(K))/OAKI))
  GOMM=GOMX/24.0
  ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0/(1.0+AMAX1(0.0_r8,(GC4X+GOMM))/EOMH)))
!
!     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     COQA=DOA concentration
!     OQKAM=Km for acetate uptake,FCNP=N,P limitation
!     VMXM=specific respiration rate
!     WFNG=water stress effect, OMA=active biomass
!     TFNX=temp stress effect, FOQA= acetate limitation
!     RGOGX=substrate-limited respiration of acetate
!     RGOGX=competition-limited respiration of acetate
!     OQA=acetate, FOQA=fraction of biological demand for acetate
!     RGOMP=O2-unlimited respiration of acetate
!     ROXY*=O2 demand, ROQCS,ROQCA=DOC, acetate demand
!     ROQCD=microbial respiration used to represent microbial activity
!
  FSBST=COQA(K)/(COQA(K)+OQKAM)
  RGOGY=AMAX1(0.0_r8,FCNP(NGL,N,K)*VMXM*WFNG*OMA(NGL,N,K))
  RGOGZ=RGOGY*FSBST*TFNX
  RGOGX=AMAX1(0.0_r8,OQA(K)*FOQA*ECHZ)
  RGOMP=AMIN1(RGOGX,RGOGZ)
  FGOCP=0.0_r8
  FGOAP=1.0_r8
  ROXYM(NGL,N,K)=0.0_r8
  ROXYP(NGL,N,K)=0.0_r8
  ROXYS(NGL,N,K)=0.0_r8
  ROQCS(NGL,N,K)=0.0_r8
  ROQAS(NGL,N,K)=RGOGZ
  ROQCD(NGL,N,K)=0.0_r8
  naqfdiag%TCH4H=naqfdiag%TCH4H+0.5*RGOMP
  end associate
  end subroutine AcetoMethanogenCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobicHeterotrophCatabolism(NGL,N,K,TFNX,WFNG,FOQC,FOQA,&
    ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC,FOQA,WFNG,TFNX
  real(r8), intent(out) :: ECHZ,RGOMP
  REAL(R8), INTENT(OUT) :: FGOCP,FGOAP
  real(r8), intent(out) :: RGOCP
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroOMcplxStateType),intent(inout):: ncplxs
  type(micfluxtype), intent(inout) :: micflx

  real(r8) :: EO2Q
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX
  real(r8) :: RGOAP
  real(r8) :: ROXYSX,ROQCSX,ROQASX
  real(r8) :: FSBST
!     begin_execution
  associate(                  &
    OMA   => nmics%OMA  ,     &
    FCNP  => nmics%FCNP ,     &
    ROXYM => nmicf%ROXYM,     &
    ROXYP=> nmicf%ROXYP ,     &
    ROQCD => nmicf%ROQCD,     &
    ZEROS => micfor%ZEROS, &
    OQC  => micstt%OQC, &
    OQA => micstt%OQA, &
    n_aero_hetrophb => micpar%n_aero_hetrophb, &
    n_anero_faculb => micpar%n_anero_faculb, &
    n_aero_fungi => micpar%n_aero_fungi, &
    n_aero_n2fixer => micpar%n_aero_n2fixer, &
    ROXYS => micflx%ROXYS, &
    ROQCS => micflx%ROQCS, &
    ROQAS => micflx%ROQAS, &
    FOCA  => ncplxs%FOCA,     &
    FOAA  => ncplxs%FOAA,     &
    COQC  => ncplxs%COQC,   &
    COQA  => ncplxs%COQA    &
  )
!     ENERGY YIELDS OF O2 REDOX REACTIONS
!     E* = growth respiration efficiency calculated in PARAMETERS
!
  IF(N.EQ.n_aero_hetrophb)THEN
    EO2Q=EO2X
  ELSEIF(N.EQ.n_anero_faculb)THEN
    EO2Q=EO2D
  ELSEIF(N.EQ.n_aero_fungi)THEN
    EO2Q=EO2G
  ELSEIF(N.EQ.n_aero_n2fixer)THEN
    EO2Q=ENFX
  ENDIF
!
! O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
! 'RGO*Z'FROM SPECIFIC RESPIRATION RATE, ACTIVE BIOMASS, DOC OR
! ACETATE CONCENTRATION,MICROBIAL C:N:P FACTOR, AND TEMPERATURE
! FOLLOWED BY POTENTIAL RESPIRATION RATES 'RGO*P' WITH UNLIMITED
! SUBSTRATE USED FOR MICROBIAL COMPETITION FACTOR

! COQC,COQA=DOC,DOA concentration, FOCA,FOAA=DOC,DOA vs DOC+DOA
! FCNP=N,P limitation,VMXO=specific respiration rate
! WFNG=water stress effect, OMA=active biomass
! TFNX=temp stress effect,FOQC,FOQA=OQC,OQA limitation
! RGOMP=O2-unlimited respiration of DOC+DOA
! RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
!
  FSBSTC=COQC(K)/(COQC(K)+OQKM)
  FSBSTA=COQA(K)/(COQA(K)+OQKA)
  FSBST=FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
  RGOCY=AMAX1(0.0_r8,FCNP(NGL,N,K)*VMXO*WFNG*OMA(NGL,N,K))
  RGOCZ=RGOCY*FSBSTC*FOCA(K)*TFNX
  RGOAZ=RGOCY*FSBSTA*FOAA(K)*TFNX
  RGOCX=AMAX1(0.0_r8,OQC(K)*FOQC*EO2Q)
  RGOAX=AMAX1(0.0_r8,OQA(K)*FOQA*EO2A)
  RGOCP=AMIN1(RGOCX,RGOCZ)
  RGOAP=AMIN1(RGOAX,RGOAZ)
  RGOMP=RGOCP+RGOAP
  IF(RGOMP.GT.ZEROS)THEN
    FGOCP=RGOCP/RGOMP
    FGOAP=RGOAP/RGOMP
  ELSE
    FGOCP=1.0
    FGOAP=0.0_r8
  ENDIF
!
! ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
! BY HETEROTROPHIC AEROBES

! ECHZ=growth respiration yield
! ROXYM,ROXYP,ROXYS=O2 demand from DOC,DOA oxidation
! ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
! ROQCD=microbial respiration used to represent microbial activity
!
  ECHZ=EO2Q*FGOCP+EO2A*FGOAP
  ROXYM(NGL,N,K)=2.667*RGOMP
  ROXYP(NGL,N,K)=ROXYM(NGL,N,K)
  ROXYSX=ROXYS(NGL,N,K)
  ROQCSX=ROQCS(NGL,N,K)
  ROQASX=ROQAS(NGL,N,K)
  ROXYS(NGL,N,K)=ROXYP(NGL,N,K)
  ROQCS(NGL,N,K)=RGOCZ
  ROQAS(NGL,N,K)=RGOAZ
  ROQCD(NGL,N,K)=RGOCY
!
  end associate
  end subroutine AerobicHeterotrophCatabolism
!------------------------------------------------------------------------------------------

  subroutine AnaerobCatabolism(NGL,N,K,TFNX,WFNG,FOQC,ECHZ,FGOCP,FGOAP,RGOMP,&
    micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC,WFNG,TFNX
  real(r8), intent(out) :: ECHZ,RGOMP
  real(r8), intent(out) :: FGOCP,FGOAP
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: GH2X,GH2F
  real(r8) :: GOAX,GOAF
  real(r8) :: GHAX
  REAL(R8) :: oxyi
  real(r8) :: RGOFX,RGOFY,RGOFZ
  real(r8) :: FSBST
!     begin_execution
  associate(                    &
    FCNP   => nmics%FCNP,       &
    OMA    => nmics%OMA ,       &
    ROXYM  => nmicf%ROXYM,      &
    ROXYP  => nmicf%ROXYP  ,    &
    ROQCD  => nmicf%ROQCD  ,    &
    TKS => micfor%TKS, &
    ZERO => micfor%ZERO, &
    OQC  => micstt%OQC, &
    CH2GS => micstt%CH2GS , &
    COXYS => micstt%COXYS, &
    ROXYS  => micflx%ROXYS, &
    ROQCS => micflx%ROQCS , &
    ROQAS => micflx%ROQAS, &
    n_anaero_ferm => micpar%n_anaero_ferm, &
    COQA    => ncplxs%COQA,  &
    COQC   => ncplxs%COQC     &
  )
  GH2X=8.3143E-03*TKS*LOG((AMAX1(1.0E-03,CH2GS)/H2KI)**4)
  GH2F=GH2X/72.0
  GOAX=8.3143E-03*TKS*LOG((AMAX1(ZERO,COQA(K))/OAKI)**2)
  GOAF=GOAX/72.0
  GHAX=GH2F+GOAF
  IF(N.EQ.n_anaero_ferm)THEN
    ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0/(1.0+AMAX1(0.0_r8,(GCHX-GHAX))/EOMF)))
  ELSE
    ECHZ=AMAX1(ENFX,AMIN1(1.0,1.0/(1.0+AMAX1(0.0_r8,(GCHX-GHAX))/EOMN)))
  ENDIF
!
!     RESPIRATION RATES BY HETEROTROPHIC ANAEROBES 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     OXYI=O2 inhibition of fermentation
!     FCNP=N,P limitation on respiration
!     VMXF=maximum respiration rate by fermenters
!     WFNG=water stress effect on respiration
!     OMA=active fermenter biomass
!     TFNX=temp stress effect, FOQC=OQC limitation
!     RFOMP=O2-unlimited respiration of DOC
!     ROQCD=microbial respiration used to represent microbial activity
!
  OXYI=1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*(-COXYS+2.5_r8)))
  FSBST=COQC(K)/(COQC(K)+OQKM)*OXYI
  RGOFY=AMAX1(0.0_r8,FCNP(NGL,N,K)*VMXF*WFNG*OMA(NGL,N,K))
  RGOFZ=RGOFY*FSBST*TFNX
  RGOFX=AMAX1(0.0_r8,OQC(K)*FOQC*ECHZ)
  RGOMP=AMIN1(RGOFX,RGOFZ)
  FGOCP=1.0_r8
  FGOAP=0.0_r8
  ROXYM(NGL,N,K)=0.0_r8
  ROXYP(NGL,N,K)=0.0_r8
  ROXYS(NGL,N,K)=0.0_r8
  ROQCS(NGL,N,K)=RGOFZ
  ROQAS(NGL,N,K)=0.0_r8
  ROQCD(NGL,N,K)=RGOFY
  naqfdiag%TRH2G=naqfdiag%TRH2G+RGOMP
!
  end associate
  end subroutine AnaerobCatabolism
!------------------------------------------------------------------------------------------

  subroutine HeteroDenitrificCatabolism(NGL,N,K,FOQC,RGOCP,&
    VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC,RGOCP
  real(r8), intent(in) :: VOLWZ
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FNO3S,FNO3B
  real(r8) :: FNO2S,FNO2B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: FNO3,FNB3
  real(r8) :: FVMXDX
  REAL(R8) :: FN2O
  real(r8) :: OQCZ3
  real(r8) :: OQCD3
  real(r8) :: OQCD3S
  real(r8) :: OQCD3B
  real(r8) :: OQCZ2
  real(r8) :: OQCD2
  real(r8) :: OQCD2S
  real(r8) :: OQCD2B
  real(r8) :: OQCZ1
  real(r8) :: OQCD1
  real(r8) :: ROXYD
  real(r8) :: RDNO3X,RDNOBX,RDNOX
  real(r8) :: RDNOT
  real(r8) :: RGOM3X,RGOMD3
  real(r8) :: RDNO2X,RDN2X,RDN2T,RGOM2X,RGOMD2,RDN2OX,RGOM1X
  real(r8) :: RGOMD1
  real(r8) :: VMXD3
  real(r8) :: VMXDXS
  real(r8) :: VMXDXB
  real(r8) :: VMXDXT
  real(r8) :: VMXD3S,VMXD3B,VMXD2,VMXD2S,VMXD2B,VMXD1
  real(r8) :: VMXD1S
  real(r8) :: ZNO3SX,ZNO3BX
  real(r8) :: ZNO2SX,ZNO2BX
  real(r8) :: Z2OSX

! begin_execution
  associate(                  &
    WFN    => nmics%WFN,      &
    FOMA   => nmics%FOMA,     &
    ROXYM => nmicf%ROXYM,     &
    ROXYO => nmicf%ROXYO,     &
    RDNO3 => nmicf%RDNO3,     &
    RDNOB => nmicf%RDNOB,     &
    RDNO2 => nmicf%RDNO2,     &
    RDN2B => nmicf%RDN2B,     &
    RDN2O => nmicf%RDN2O,     &
    RGOMD  => nmicf%RGOMD,    &
    RGOMY  => nmicf%RGOMY,    &
    OSRH  => ncplxs%OSRH ,    &
    RNO2Y  => micfor%RNO2Y,   &
    RN2OY  => micfor%RN2OY,   &
    RNO3Y  => micfor%RNO3Y , &
    VLNO3  => micfor%VLNO3, &
    ZERO   => micfor%ZERO, &
    ZEROS  => micfor%ZEROS, &
    ZEROS2 => micfor%ZEROS2, &
    RN2BY => micfor%RN2BY, &
    RN3BY => micfor%RN3BY, &
    VLNOB => micfor%VLNOB, &
    CNO3B => micstt%CNO3B, &
    CNO3S => micstt%CNO3S, &
    CZ2OS => micstt%CZ2OS, &
    Z2OS => micstt%Z2OS, &
    ZNO2B => micstt%ZNO2B, &
    ZNO2S => micstt%ZNO2S, &
    ZNO3B => micstt%ZNO3B, &
    ZNO3S => micstt%ZNO3S, &
    CNO2B => micstt%CNO2B , &
    CNO2S => micstt%CNO2S, &
    FOSRH => micstt%FOSRH, &
    CH2GS => micstt%CH2GS, &
    OQC => micstt%OQC, &
    RVMX3 => micflx%RVMX3,  &
    RVMX2 => micflx%RVMX2, &
    RVMX1 => micflx%RVMX1, &
    RVMB2 => micflx%RVMB2, &
    RVMB3 => micflx%RVMB3   &
  )
!
! FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
! AND ROOT POPULATIONS
!
! FNO3,FNB3=fraction of total biological demand for NO3
!

  FNO3S=VLNO3
  FNO3B=VLNOB
  IF(RNO3Y.GT.ZEROS)THEN
    FNO3=AMAX1(FMN,RVMX3(NGL,N,K)/RNO3Y)
  ELSE
    FNO3=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3)
  ENDIF
  IF(RN3BY.GT.ZEROS)THEN
    FNB3=AMAX1(FMN,RVMB3(NGL,N,K)/RN3BY)
  ELSE
    FNB3=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB)
  ENDIF
  naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3
  naqfdiag%TFNO3B=naqfdiag%TFNO3B+FNB3
!
!     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
!
!     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
!     VMXD3=demand for NO3-N reduction
!     VMXDXS,VMXDXB=maximum NO3 reduction in non-band, band
!     FNO3S,FNO3B=fractions of total NO3 in non-band, band
!     CNO3S,CNO3B=NO3 concentrations in non-band, band
!     Z3KM,Z2KM=Km for NO3, NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD3S,VMXD3B=substrate-unlimited NO3 reduction in non-band,band
!     OQCD3S,OQCD3B=DOC limitation to NO3 reduction in non-band, band
!     RDNO3,RDNOB=substrate-limited NO3 reduction in non-band,band
!     RGOM3X,RGOMD3=substrate-unltd,-ltd respn from NO3 reduction
!     RVMX3,RVMB3=demand for NO3 reduction in non-band,band
!
  ROXYD=AMAX1(0.0_r8,ROXYM(NGL,N,K)-ROXYO(NGL,N,K))
  VMXD3=0.875*ROXYD
  IF(CNO3S.GT.ZERO)THEN
    VMXDXS=FNO3S*VMXD3*CNO3S/(CNO3S+Z3KM) &
      /(1.0+(CNO2S*Z3KM)/(CNO3S*Z2KM))
  ELSE
    VMXDXS=0.0_r8
  ENDIF
  IF(CNO3B.GT.ZERO)THEN
    VMXDXB=FNO3B*VMXD3*CNO3B/(CNO3B+Z3KM) &
      /(1.0+(CNO2B*Z3KM)/(CNO3B*Z2KM))
  ELSE
    VMXDXB=0.0_r8
  ENDIF
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2.AND.FOSRH(K).GT.ZERO)THEN
    FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD3S=VMXDXS*FVMXDX
  VMXD3B=VMXDXB*FVMXDX
  OQCZ3=AMAX1(0.0_r8,OQC(K)*FOQC-RGOCP*WFN(NGL,N,K))
  OQCD3=OQCZ3/ECN3
  OQCD3S=OQCD3*FNO3S
  OQCD3B=OQCD3*FNO3B
  ZNO3SX=ZNO3S*FNO3
  ZNO3BX=ZNO3B*FNB3
  RDNO3X=AMAX1(0.0_r8,AMIN1(ZNO3SX,VMXD3S))
  RDNOBX=AMAX1(0.0_r8,AMIN1(ZNO3BX,VMXD3B))
  RDNO3(NGL,N,K)=AMAX1(0.0_r8,AMIN1(VMXD3S,OQCD3S,ZNO3SX))
  RDNOB(NGL,N,K)=AMAX1(0.0_r8,AMIN1(VMXD3B,OQCD3B,ZNO3BX))
  RDNOX=RDNO3X+RDNOBX
  RDNOT=RDNO3(NGL,N,K)+RDNOB(NGL,N,K)
  RGOM3X=ECN3*RDNOX
  RGOMD3=ECN3*RDNOT
  RVMX3(NGL,N,K)=VMXD3S
  RVMB3(NGL,N,K)=VMXD3B
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
  FNO2S=VLNO3
  FNO2B=VLNOB
  IF(RNO2Y.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RVMX2(NGL,N,K)/RNO2Y)
  ELSE
    FNO2=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3)
  ENDIF
  IF(RN2BY.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RVMB2(NGL,N,K)/RN2BY)
  ELSE
    FNB2=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB)
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 AND NO3 IN BAND AND NON-BAND SOIL ZONES
!
!     VMXD2=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM,Z1KM=Km for NO2, N2O uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD2S,VMXD2B=substrate-unlimited NO2 reduction in non-band,band
!     OQCD2S,OQCD2B=DOC limitation to NO2 reduction in non-band, band
!     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
!     RGOM2X,RGOMD2=substrate-unltd,-ltd respn from NO2 reduction
!
  VMXD2=VMXD3-RDNOT
  IF(CNO2S.GT.ZERO)THEN
    VMXDXS=FNO2S*VMXD2*CNO2S/(CNO2S+Z2KM) &
      /(1.0+(CZ2OS*Z2KM)/(CNO2S*Z1KM))
  ELSE
    VMXDXS=0.0_r8
  ENDIF
  IF(CNO2B.GT.ZERO)THEN
    VMXDXB=FNO2B*VMXD2*CNO2B/(CNO2B+Z2KM) &
      /(1.0+(CZ2OS*Z2KM)/(CNO2B*Z1KM))
  ELSE
    VMXDXB=0.0_r8
  ENDIF
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2.AND.FOSRH(K).GT.ZERO)THEN
    FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD2S=VMXDXS*FVMXDX
  VMXD2B=VMXDXB*FVMXDX
  OQCZ2=AMAX1(0.0_r8,OQCZ3-RGOMD3)
  OQCD2=OQCZ2/ECN2
  OQCD2S=OQCD2*FNO3S
  OQCD2B=OQCD2*FNO3B
  ZNO2SX=(ZNO2S+RDNO3(NGL,N,K))*FNO2
  ZNO2BX=(ZNO2B+RDNOB(NGL,N,K))*FNB2
  RDNO2X=AMAX1(0.0_r8,AMIN1(ZNO2SX,VMXD2S))
  RDNOBX=AMAX1(0.0_r8,AMIN1(ZNO2BX,VMXD2B))
  RDNO2(NGL,N,K)=AMAX1(0.0_r8,AMIN1(VMXD2S,OQCD2S,ZNO2SX))
  RDN2B(NGL,N,K)=AMAX1(0.0_r8,AMIN1(VMXD2B,OQCD2B,ZNO2BX))
  RDN2X=RDNO2X+RDNOBX
  RDN2T=RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
  RGOM2X=ECN2*RDN2X
  RGOMD2=ECN2*RDN2T
  RVMX2(NGL,N,K)=VMXD2S
  RVMB2(NGL,N,K)=VMXD2B
!
!     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
!     AND ROOT POPULATIONS
!
!     FN2O=fraction of total biological demand for N2O
!
  IF(RN2OY.GT.ZEROS)THEN
    FN2O=AMAX1(FMN,RVMX1(NGL,N,K)/RN2OY)
  ELSE
    FN2O=AMAX1(FMN,FOMA(NGL,N,K))
  ENDIF
  naqfdiag%TFN2OX=naqfdiag%TFN2OX+FN2O
!
!     N2O REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS N2O
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2, NO3 AND NO2 IN BAND AND NON-BAND SOIL ZONES
!
!     VMXD1=demand for N2O-N reduction
!     VMXDXS=maximum N2O reduction
!     CZ2OS=N2O concentrations
!     Z1KM=Km for N2O uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD1S=substrate-unlimited N2O reduction
!     OQCD1=DOC limitation to N2O reduction
!     RDN2O=substrate-limited N2O reduction
!     RGOM1X,RGOMD1=substrate-unltd,-ltd  respn from N2O reduction
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NOx reduction
!     RVMX1=demand for N2O reduction
!
  VMXD1=(VMXD2-RDN2T)*2.0
  VMXDXS=VMXD1*CZ2OS/(CZ2OS+Z1KM)
  IF(VOLWZ.GT.ZEROS2.AND.FOSRH(K).GT.ZERO)THEN
    FVMXDX=1.0/(1.0+VMXDXS/(VMKI*VOLWZ*FOSRH(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD1S=VMXDXS*FVMXDX
  OQCZ1=AMAX1(0.0_r8,OQCZ2-RGOMD2)
  OQCD1=OQCZ1/ECN1
  Z2OSX=(Z2OS+RDN2T)*FN2O
  RDN2OX=AMAX1(0.0_r8,AMIN1(Z2OSX,VMXD1S))
  RDN2O(NGL,N,K)=AMAX1(0.0_r8,AMIN1(VMXD1S,OQCD1,Z2OSX))
  RGOM1X=ECN1*RDN2OX
  RGOMD1=ECN1*RDN2O(NGL,N,K)
  RGOMY(NGL,N,K)=RGOM3X+RGOM2X+RGOM1X
  RGOMD(NGL,N,K)=RGOMD3+RGOMD2+RGOMD1
  RVMX1(NGL,N,K)=VMXD1S
  end associate
  end subroutine HeteroDenitrificCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobsO2Uptake(NGL,N,K,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
    micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: OXKX,FOXYX,RGOMP,RVOXP
  real(r8), intent(in) :: RVOXPA
  real(r8), intent(in) :: RVOXPB
  type(micforctype), intent(in) :: micfor
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
    WFN    => nmics%WFN,      &
    OMA    => nmics%OMA,      &
    RUPOX  => nmicf%RUPOX,    &
    RGOMO  => nmicf%RGOMO,    &
    ROXYM => nmicf%ROXYM,     &
    ROXYP=> nmicf%ROXYP ,     &
    ROXYO => nmicf%ROXYO,     &
    RH2GX  => nmicf%RH2GX,    &
    ROQCD  => nmicf%ROQCD,    &
    RCO2X  => nmicf%RCO2X,    &
    RCH3X  => nmicf%RCH3X,    &
    RCH4X  => nmicf%RCH4X,    &
    RVOXA  => nmicf%RVOXA,    &
    RVOXB  => nmicf%RVOXB,    &
    COXQ   => micfor%COXQ, &
    COXR   => micfor%COXR, &
    COXYE => micfor%COXYE, &
    ROXYF  => micfor%ROXYF,   &
    ROXYL => micfor%ROXYL, &
    FLQRI  => micfor%FLQRI, &
    FLQRQ => micfor%FLQRQ , &
    litrm => micfor%litrm , &
    OLSGL  => micfor%OLSGL, &
    VOLX => micfor%VOLX, &
    VOLY  => micfor%VOLY, &
    ZERO => micfor%ZERO, &
    ZEROS  => micfor%ZEROS, &
    VOLPM => micfor%VOLPM, &
    VOLW  => micfor%VOLW, &
    VOLWM => micfor%VOLWM, &
    THETPM => micfor%THETPM, &
    DFGS => micfor%DFGS, &
    FILM => micfor%FILM, &
    TORT  => micfor%TORT, &
    OXYG => micstt%OXYG, &
    OXYS => micstt%OXYS, &
    COXYS => micstt%COXYS, &
    SOXYL => micstt%SOXYL, &
    COXYG => micstt%COXYG, &
    RVMX4 => micflx%RVMX4, &
    RVMB4 => micflx%RVMB4, &
    RVMX2 => micflx%RVMX2, &
    RVMB2 => micflx%RVMB2, &
    ROXSK => micflx%ROXSK &
  )
  IF(ROXYP(NGL,N,K).GT.ZEROS.AND.FOXYX.GT.ZERO)THEN
    IF(.not.litrm.OR.VOLX.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX=ROXYP(NGL,N,K)*XNPG
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
      DO 420 M=1,NPH
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
        !write(*,*)'VOLY=',VOLY
        THETW1=AMAX1(0.0_r8,safe_adb(VOLWM(M),VOLY))
        RRADO=ORAD*(FILM(M)+ORAD)/FILM(M)
        DIFOX=TORT(M)*OLSGL1*12.57_r8*BIOS*OMA(NGL,N,K)*RRADO
        VOLWOX=VOLWM(M)*SOXYL
        VOLPOX=VOLPM(M)
        VOLWPM=VOLWOX+VOLPOX
        DO 425 MX=1,NPT
          OXYG1=OXYG1+ROXYFX
          OXYS1=OXYS1+ROXYLX
          COXYS1=AMIN1(COXYE*SOXYL,AMAX1(0.0_r8,safe_adb(OXYS1,(VOLWM(M)*FOXYX))))
          X=DIFOX*COXYS1
          IF(X.GT.ZEROS.AND.OXYS1.GT.ZEROS)THEN
            B=-RUPMX-DIFOX*OXKX-X
            C=X*RUPMX
            RMPOX=(-B-SQRT(B*B-4.0_r8*C))/2.0_r8
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
          RUPOX(NGL,N,K)=RUPOX(NGL,N,K)+RMPOX
          ROXSK(M)=ROXSK(M)+RMPOX

425     CONTINUE
420   CONTINUE
      !write(*,*)'420'
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
      !
      !     WFN=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RVMX2,RVMB2=NH3,NO2 oxidation in non-band, band
      !
      WFN(NGL,N,K)=AMIN1(1.0,AMAX1(0.0_r8,RUPOX(NGL,N,K)/ROXYP(NGL,N,K)))
!     IF(K.LE.4)THEN
!       ROQCS(NGL,N,K)=ROQCS(NGL,N,K)*WFN(NGL,N,K)
!       ROQAS(NGL,N,K)=ROQAS(NGL,N,K)*WFN(NGL,N,K)
!       ROQCD(NGL,N,K)=ROQCD(NGL,N,K)*WFN(NGL,N,K)
!     ENDIF
      IF(K.EQ.5)THEN
        IF(N.EQ.1)THEN
          RVMX4(NGL,N,K)=RVMX4(NGL,N,K)*WFN(NGL,N,K)
          RVMB4(NGL,N,K)=RVMB4(NGL,N,K)*WFN(NGL,N,K)
        ELSEIF(N.EQ.2)THEN
          RVMX2(NGL,N,K)=RVMX2(NGL,N,K)*WFN(NGL,N,K)
          RVMB2(NGL,N,K)=RVMB2(NGL,N,K)*WFN(NGL,N,K)
        ENDIF
      ENDIF
    ELSE
      RUPOX(NGL,N,K)=ROXYP(NGL,N,K)
      WFN(NGL,N,K)=1.0
    ENDIF
  ELSE
    RUPOX(NGL,N,K)=0.0_r8
    WFN(NGL,N,K)=1.0
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RGOMO,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2X,RCH3X,RCH4X,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
  !     ROXYO=O2-limited O2 uptake
  !     RVOXA,RVOXB=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !
  RGOMO(NGL,N,K)=RGOMP*WFN(NGL,N,K)
  RCO2X(NGL,N,K)=RGOMO(NGL,N,K)
  RCH3X(NGL,N,K)=0.0_r8
  RCH4X(NGL,N,K)=0.0_r8
  ROXYO(NGL,N,K)=ROXYM(NGL,N,K)*WFN(NGL,N,K)
  RH2GX(NGL,N,K)=0.0_r8
  IF(K.EQ.5)THEN
    RVOXA(NGL,N)=RVOXPA*WFN(NGL,N,K)
    RVOXB(NGL,N)=RVOXPB*WFN(NGL,N,K)
  ENDIF
  !write(*,*)'finish AerobsO2Uptake'
  end associate
  end subroutine AerobsO2Uptake

!------------------------------------------------------------------------------------------

  subroutine BiomassMineralization(NGL,N,K,FNH4X,&
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: FNH4X
  real(r8), intent(in) :: FNB3X,FNB4X,FNO3X
  real(r8), intent(in) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
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
   TFNG    => nmics%TFNG  ,  &
   OMA     => nmics%OMA   ,  &
   FNH4XR  => nmicf%FNH4XR,  &
   FNO3XR  => nmicf%FNO3XR,  &
   FPO4XR  => nmicf%FPO4XR,  &
   FP14XR  => nmicf%FP14XR,  &
   RINH4 => nmicf%RINH4   ,  &
   RINO3  => nmicf%RINO3  ,  &
   RIPO4  => nmicf%RIPO4  ,  &
   RINB4 => nmicf%RINB4   ,  &
   RINB3 => nmicf%RINB3   ,  &
   RIPOB => nmicf%RIPOB   ,  &
   RINH4R  => nmicf%RINH4R,  &
   RINO3R  => nmicf%RINO3R,  &
   RIPO4R  => nmicf%RIPO4R,  &
   RIP14  => nmicf%RIP14,    &
   RIP1B  => nmicf%RIP1B,    &
   RIP14R  => nmicf%RIP14R,  &
   CNOMC  => micpar%CNOMC   , &
   CPOMC  => micpar%CPOMC ,   &
   litrm => micfor%litrm , &
   VLNH4  => micfor%VLNH4 , &
   VLNHB  => micfor%VLNHB , &
   VLNO3  => micfor%VLNO3, &
   VLNOB  => micfor%VLNOB, &
   VLPOB => micfor%VLPOB , &
   VOLW  => micfor%VOLW , &
   VOLWU => micfor%VOLWU, &
   VLPO4  => micfor%VLPO4 , &
   ZNH4B  => micstt%ZNH4B, &
   ZNH4S  => micstt%ZNH4S, &
   ZNO3B  => micstt%ZNO3B, &
   ZNO3S  => micstt%ZNO3S, &
   ZNH4TU => micstt%ZNH4TU, &
   ZNO3TU => micstt%ZNO3TU, &
   H1P4TU => micstt%H1P4TU, &
   H2P4TU  => micstt%H2P4TU, &
   CNH4BU => micstt%CNH4BU, &
   CNH4SU => micstt%CNH4SU, &
   CH2P4U => micstt%CH2P4U, &
   CH2P4BU => micstt%CH2P4BU, &
   CH1P4U => micstt%CH1P4U, &
   CH1P4BU => micstt%CH1P4BU, &
   CH2P4  => micstt%CH2P4, &
   CH2P4B  => micstt%CH2P4B, &
   CNH4B  => micstt%CNH4B , &
   CNH4S  => micstt%CNH4S, &
   CH1P4  => micstt%CH1P4 , &
   CH1P4B  => micstt%CH1P4B , &
   H1PO4  => micstt%H1PO4 , &
   H1POB => micstt%H1POB, &
   H2PO4  => micstt%H2PO4, &
   H2POB => micstt%H2POB, &
   CNO3B  => micstt%CNO3B, &
   CNO3S  => micstt%CNO3S, &
   CNO3SU => micstt%CNO3SU, &
   CNO3BU => micstt%CNO3BU, &
   OMC  => micstt%OMC, &
   OMN  => micstt%OMN, &
   OMP  => micstt%OMP, &
   RINOB => micflx%RINOB, &
   RINOO  => micflx%RINOO, &
   RINHO  => micflx%RINHO, &
   RINHB  => micflx%RINHB , &
   RIPOO  => micflx%RIPOO, &
   RIPBO   => micflx%RIPBO, &
   RIPO1  => micflx%RIPO1, &
   RIPB1 => micflx%RIPB1, &
   RINHOR => micflx%RINHOR, &
   RINOOR => micflx%RINOOR, &
   RIPOOR  => micflx%RIPOOR, &
   RIPO1R => micflx%RIPO1R, &
   TRINH4  => micflx%TRINH4, &
   TRIPO4 => micflx%TRIPO4 &
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
  RINHP=(OMC(3,NGL,N,K)*CNOMC(3,NGL,N,K)-OMN(3,NGL,N,K))
  IF(RINHP.GT.0.0)THEN
    CNH4X=AMAX1(0.0_r8,CNH4S-Z4MN)
    CNH4Y=AMAX1(0.0_r8,CNH4B-Z4MN)
    RINHX=AMIN1(RINHP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*Z4MX)
    RINHO(NGL,N,K)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RINHB(NGL,N,K)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M=Z4MN*VOLWU*FNH4S
    ZNHBM=Z4MN*VOLWU*FNHBS
    RINH4(NGL,N,K)=AMIN1(FNH4X*AMAX1(0.0_r8,(ZNH4S-ZNH4M)),RINHO(NGL,N,K))
    RINB4(NGL,N,K)=AMIN1(FNB4X*AMAX1(0.0_r8,(ZNH4B-ZNHBM)),RINHB(NGL,N,K))
  ELSE
    RINHO(NGL,N,K)=0.0_r8
    RINHB(NGL,N,K)=0.0_r8
    RINH4(NGL,N,K)=RINHP*FNH4S
    RINB4(NGL,N,K)=RINHP*FNHBS
  ENDIF
  TRINH4=TRINH4+(RINH4(NGL,N,K)+RINB4(NGL,N,K))
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
  RINOP=AMAX1(0.0_r8,RINHP-RINH4(NGL,N,K)-RINB4(NGL,N,K))
  IF(RINOP.GT.0.0)THEN
    CNO3X=AMAX1(0.0_r8,CNO3S-ZOMN)
    CNO3Y=AMAX1(0.0_r8,CNO3B-ZOMN)
    RINOX=AMIN1(RINOP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*ZOMX)
    RINOO(NGL,N,K)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RINOB(NGL,N,K)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M=ZOMN*VOLWU*FNO3S
    ZNOBM=ZOMN*VOLWU*FNO3B
    RINO3(NGL,N,K)=AMIN1(FNO3X*AMAX1(0.0_r8,(ZNO3S-ZNO3M)) &
      ,RINOO(NGL,N,K))
    RINB3(NGL,N,K)=AMIN1(FNB3X*AMAX1(0.0_r8,(ZNO3B-ZNOBM)) &
      ,RINOB(NGL,N,K))
  ELSE
    RINOO(NGL,N,K)=0.0_r8
    RINOB(NGL,N,K)=0.0_r8
    RINO3(NGL,N,K)=RINOP*FNO3S
    RINB3(NGL,N,K)=RINOP*FNO3B
  ENDIF
  TRINH4=TRINH4+(RINO3(NGL,N,K)+RINB3(NGL,N,K))
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
  RIPOP=(OMC(3,NGL,N,K)*CPOMC(3,NGL,N,K)-OMP(3,NGL,N,K))
  IF(RIPOP.GT.0.0)THEN
    CH2PX=AMAX1(0.0_r8,CH2P4-HPMN)
    CH2PY=AMAX1(0.0_r8,CH2P4B-HPMN)
    RIPOX=AMIN1(RIPOP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX)
    RIPOO(NGL,N,K)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RIPBO(NGL,N,K)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM=HPMN*VOLW*FH2PS
    H2PBM=HPMN*VOLW*FH2PB
    RIPO4(NGL,N,K)=AMIN1(FPO4X*AMAX1(0.0_r8,(H2PO4-H2POM)) &
      ,RIPOO(NGL,N,K))
    RIPOB(NGL,N,K)=AMIN1(FPOBX*AMAX1(0.0_r8,(H2POB-H2PBM)) &
      ,RIPBO(NGL,N,K))
  ELSE
    RIPOO(NGL,N,K)=0.0_r8
    RIPBO(NGL,N,K)=0.0_r8
    RIPO4(NGL,N,K)=RIPOP*FH2PS
    RIPOB(NGL,N,K)=RIPOP*FH2PB
  ENDIF
  TRIPO4=TRIPO4+(RIPO4(NGL,N,K)+RIPOB(NGL,N,K))
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
  RIP1P=0.1*AMAX1(0.0_r8,RIPOP-RIPO4(NGL,N,K)-RIPOB(NGL,N,K))
  IF(RIP1P.GT.0.0)THEN
    CH1PX=AMAX1(0.0_r8,CH1P4-HPMN)
    CH1PY=AMAX1(0.0_r8,CH1P4B-HPMN)
    RIP1X=AMIN1(RIP1P,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX)
    RIPO1(NGL,N,K)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RIPB1(NGL,N,K)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM=HPMN*VOLW*FH1PS
    H1PBM=HPMN*VOLW*FH1PB
    RIP14(NGL,N,K)=AMIN1(FP14X*AMAX1(0.0_r8,(H1PO4-H1POM)) &
      ,RIPO1(NGL,N,K))
    RIP1B(NGL,N,K)=AMIN1(FP1BX*AMAX1(0.0_r8,(H1POB-H1PBM)) &
      ,RIPB1(NGL,N,K))
  ELSE
    RIPO1(NGL,N,K)=0.0_r8
    RIPB1(NGL,N,K)=0.0_r8
    RIP14(NGL,N,K)=RIP1P*FH1PS
    RIP1B(NGL,N,K)=RIP1P*FH1PB
  ENDIF
  TRIPO4=TRIPO4+(RIP14(NGL,N,K)+RIP1B(NGL,N,K))
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
    RINHPR=RINHP-RINH4(NGL,N,K)-RINO3(NGL,N,K)
    IF(RINHPR.GT.0.0)THEN
      CNH4X=AMAX1(0.0_r8,CNH4SU-Z4MN)
      CNH4Y=AMAX1(0.0_r8,CNH4BU-Z4MN)
      RINHOR(NGL,N,K)=AMIN1(RINHPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VOLW
      RINH4R(NGL,N,K)=AMIN1(FNH4XR(NGL,N,K)*AMAX1(0.0 &
        ,(ZNH4TU-ZNH4M)),RINHOR(NGL,N,K))
    ELSE
      RINHOR(NGL,N,K)=0.0_r8
      RINH4R(NGL,N,K)=RINHPR
    ENDIF
    TRINH4=TRINH4+RINH4R(NGL,N,K)
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
    RINOPR=AMAX1(0.0_r8,RINHPR-RINH4R(NGL,N,K))
    IF(RINOPR.GT.0.0)THEN
      CNO3X=AMAX1(0.0_r8,CNO3SU-ZOMN)
      CNO3Y=AMAX1(0.0_r8,CNO3BU-ZOMN)
      RINOOR(NGL,N,K)=AMAX1(RINOPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VOLW
      RINO3R(NGL,N,K)=AMIN1(FNO3XR(NGL,N,K)*AMAX1(0.0 &
        ,(ZNO3TU-ZNO3M)),RINOOR(NGL,N,K))
    ELSE
      RINOOR(NGL,N,K)=0.0_r8
      RINO3R(NGL,N,K)=RINOPR
    ENDIF
    TRINH4=TRINH4+RINO3R(NGL,N,K)
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
    RIPOPR=RIPOP-RIPO4(NGL,N,K)
    IF(RIPOPR.GT.0.0)THEN
      CH2PX=AMAX1(0.0_r8,CH2P4U-HPMN)
      CH2PY=AMAX1(0.0_r8,CH2P4BU-HPMN)
      RIPOOR(NGL,N,K)=AMIN1(RIPOPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M=HPMN*VOLWU
      RIPO4R(NGL,N,K)=AMIN1(FPO4XR(NGL,N,K)*AMAX1(0.0 &
        ,(H2P4TU-H2P4M)),RIPOOR(NGL,N,K))
    ELSE
      RIPOOR(NGL,N,K)=0.0_r8
      RIPO4R(NGL,N,K)=RIPOPR
    ENDIF
    TRIPO4=TRIPO4+RIPO4R(NGL,N,K)
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
    RIP1PR=0.1*AMAX1(0.0_r8,RIPOPR-RIPO4R(NGL,N,K))
    IF(RIP1PR.GT.0.0)THEN
      CH1PX=AMAX1(0.0_r8,CH1P4U-HPMN)
      CH1PY=AMAX1(0.0_r8,CH1P4BU-HPMN)
      RIPO1R(NGL,N,K)=AMIN1(RIP1PR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M=HPMN*VOLWU
      RIP14R(NGL,N,K)=AMIN1(FP14XR(NGL,N,K)*AMAX1(0.0 &
        ,(H1P4TU-H1P4M)),RIPO1R(NGL,N,K))
    ELSE
      RIPO1R(NGL,N,K)=0.0_r8
      RIP14R(NGL,N,K)=RIP1PR
    ENDIF
    TRIPO4=TRIPO4+RIP14R(NGL,N,K)
  ENDIF
  end associate
  end subroutine BiomassMineralization
!------------------------------------------------------------------------------------------

  subroutine GatherMicrobialRespiration(NGL,N,K,RMOMK,RGOMT,RXOMT,RMOMT,&
    micfor,micstt,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(out) :: RGOMT,RXOMT
  real(r8), intent(out) :: RMOMT
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
!     begin_execution
  associate(                  &
    TFNR  => nmics%TFNR,      &
    OMN2  => nmics%OMN2,      &
    RGN2F => nmicf%RGN2F,     &
    RGOMO  => nmicf%RGOMO,    &
    RMOMC => nmicf%RMOMC ,    &
    RINH4 => nmicf%RINH4 ,    &
    RINO3  => nmicf%RINO3,    &
    RN2FX  => nmicf%RN2FX,    &
    CNOMC  => micpar%CNOMC   , &
    CPOMC  => micpar%CPOMC  ,  &
    pH  => micfor%pH, &
    ZEROS => micfor%ZEROS, &
    n_aero_n2fixer => micpar%n_aero_n2fixer, &
    n_anero_n2fixer => micpar%n_anero_n2fixer ,&
    CZ2GS => micstt%CZ2GS, &
    OMC => micstt%OMC , &
    OMN => micstt%OMN  &
  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TFNR=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  FPH=1.0+AMAX1(0.0_r8,0.25*(6.5-PH))

  RMOMX=RMOM*TFNR(NGL,N,K)*FPH
  RMOMC(1,NGL,N,K)=OMN(1,NGL,N,K)*RMOMX*RMOMK(1)
  RMOMC(2,NGL,N,K)=OMN2(NGL,N,K)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMOMT=total maintenance respiration
!     RGOMT=growth respiration
!     RXOMT=senescence respiration
!
  RMOMT=RMOMC(1,NGL,N,K)+RMOMC(2,NGL,N,K)
  RGOMT=AMAX1(0.0_r8,RGOMO(NGL,N,K)-RMOMT)
  RXOMT=AMAX1(0.0_r8,RMOMT-RGOMO(NGL,N,K))
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
  IF(N.EQ.n_aero_n2fixer.OR.N.EQ.n_anero_n2fixer)THEN
    RGN2P=AMAX1(0.0_r8,OMC(3,NGL,N,K)*CNOMC(3,NGL,N,K) &
      -OMN(3,NGL,N,K))/EN2F(N)
    IF(RGOMT.GT.ZEROS)THEN
      RGN2F(NGL,N,K)=AMIN1(RGOMT*RGN2P/(RGOMT+RGN2P) &
        *CZ2GS/(CZ2GS+ZFKM),OMGR*OMC(3,NGL,N,K))
    ELSE
      RGN2F(NGL,N,K)=0.0_r8
    ENDIF
    RN2FX(NGL,N,K)=RGN2F(NGL,N,K)*EN2F(N)
  ENDIF
  end associate
  end subroutine GatherMicrobialRespiration
!------------------------------------------------------------------------------------------

  subroutine GetMicrobialAnabolismFlux(NGL,N,K,ECHZ,FGOCP,FGOAP,&
    RGOMT,RXOMT,RMOMT,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: ECHZ
  real(r8), intent(in) :: FGOCP,FGOAP,RGOMT,RXOMT,RMOMT
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  integer :: M
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
    CNOMA  => nmics%CNOMA,    &
    CPOMA  => nmics%CPOMA,    &
    TFNG   => nmics%TFNG ,    &
    FCN    => nmics%FCN  ,    &
    FCP    => nmics%FCP  ,    &
    FOMK   => nmics%FOMK ,    &
    CGOMC  => nmicf%CGOMC,    &
    CGOMN  => nmicf%CGOMN,    &
    CGOMP  => nmicf%CGOMP,    &
    CGOQC  => nmicf%CGOQC,    &
    CGOAC  => nmicf%CGOAC,    &
    CGOMS  => nmicf%CGOMS,    &
    CGONS  => nmicf%CGONS,    &
    CGOPS  => nmicf%CGOPS,    &
    RGOMO  => nmicf%RGOMO,    &
    RGOMD  => nmicf%RGOMD,    &
    RMOMC => nmicf%RMOMC ,    &
    RDOMC  => nmicf%RDOMC,    &
    RDOMN => nmicf%RDOMN ,    &
    RDOMP => nmicf%RDOMP ,    &
    RHOMC => nmicf%RHOMC ,    &
    RHOMN => nmicf%RHOMN ,    &
    RHOMP => nmicf%RHOMP ,    &
    RCOMC  => nmicf%RCOMC,    &
    RCOMN  => nmicf%RCOMN,    &
    RCOMP  => nmicf%RCOMP,    &
    RDMMC  => nmicf%RDMMC,    &
    RHMMC  => nmicf%RHMMC,    &
    RCMMC  => nmicf%RCMMC,    &
    RDMMN  => nmicf%RDMMN,    &
    RHMMN  => nmicf%RHMMN,    &
    RCMMN  => nmicf%RCMMN,    &
    RDMMP  => nmicf%RDMMP,    &
    RHMMP  => nmicf%RHMMP,    &
    RCMMP  => nmicf%RCMMP,    &
    RXOMC  => nmicf%RXOMC,    &
    RXOMN  => nmicf%RXOMN,    &
    RXOMP  => nmicf%RXOMP,    &
    R3OMC  => nmicf%R3OMC,    &
    R3OMN  => nmicf%R3OMN,    &
    R3OMP  => nmicf%R3OMP,    &
    RXMMC  => nmicf%RXMMC,    &
    RXMMN  => nmicf%RXMMN,    &
    RXMMP  => nmicf%RXMMP,    &
    R3MMC  => nmicf%R3MMC,    &
    R3MMN  => nmicf%R3MMN,    &
    R3MMP  => nmicf%R3MMP,    &
    RCO2X  => nmicf%RCO2X,    &
    RGN2F  => nmicf%RGN2F,    &
    TCGOQC  => ncplxf%TCGOQC, &
    TCGOAC  => ncplxf%TCGOAC, &
    TCGOMN  => ncplxf%TCGOMN, &
    TCGOMP  => ncplxf%TCGOMP, &
    CNQ  => ncplxs%CNQ,       &
    CPQ  => ncplxs%CPQ,       &
    CNOMC  => micpar%CNOMC,     &
    CPOMC  => micpar%CPOMC,    &
    FL       => micpar%FL    ,     &
    EHUM     => micstt%EHUM,  &
    OMC  => micstt%OMC , &
    OMN  => micstt%OMN , &
    OMP  => micstt%OMP, &
    OQN   => micstt%OQN , &
    OQP  => micstt%OQP , &
    ZEROS  => micfor%ZEROS , &
    ZERO   => micfor%ZERO  &
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

  CGOMX=AMIN1(RMOMT,RGOMO(NGL,N,K))+RGN2F(NGL,N,K) &
    +(RGOMT-RGN2F(NGL,N,K))/ECHZ
  CGOMD=RGOMD(NGL,N,K)/ENOX
  CGOMC(NGL,N,K)=CGOMX+CGOMD

  CGOQC(NGL,N,K)=CGOMX*FGOCP+CGOMD
  CGOAC(NGL,N,K)=CGOMX*FGOAP
  CGOXC=CGOQC(NGL,N,K)+CGOAC(NGL,N,K)
  CGOMN(NGL,N,K)=AMAX1(0.0_r8,AMIN1(OQN(K)*FOMK(NGL,N,K) &
    ,CGOXC*CNQ(K)/FCN(NGL,N,K)))
  CGOMP(NGL,N,K)=AMAX1(0.0_r8,AMIN1(OQP(K)*FOMK(NGL,N,K) &
    ,CGOXC*CPQ(K)/FCP(NGL,N,K)))
  TCGOQC(K)=TCGOQC(K)+CGOQC(NGL,N,K)
  TCGOAC(K)=TCGOAC(K)+CGOAC(NGL,N,K)
  TCGOMN(K)=TCGOMN(K)+CGOMN(NGL,N,K)
  TCGOMP(K)=TCGOMP(K)+CGOMP(NGL,N,K)
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
  IF(OMC(3,NGL,N,K).GT.ZEROS &
    .AND.OMC(1,NGL,N,K).GT.ZEROS)THEN
    CCC=AMAX1(0.0_r8,AMIN1(1.0 &
      ,OMN(3,NGL,N,K)/(OMN(3,NGL,N,K) &
      +OMC(3,NGL,N,K)*CNOMC(3,NGL,N,K)) &
      ,OMP(3,NGL,N,K)/(OMP(3,NGL,N,K) &
      +OMC(3,NGL,N,K)*CPOMC(3,NGL,N,K))))
    CXC=OMC(3,NGL,N,K)/OMC(1,NGL,N,K)
    C3C=1.0/(1.0+CXC/CKC)
    CNC=AMAX1(0.0_r8,AMIN1(1.0 &
      ,OMC(3,NGL,N,K)/(OMC(3,NGL,N,K) &
      +OMN(3,NGL,N,K)/CNOMC(3,NGL,N,K))))
    CPC=AMAX1(0.0_r8,AMIN1(1.0 &
      ,OMC(3,NGL,N,K)/(OMC(3,NGL,N,K) &
      +OMP(3,NGL,N,K)/CPOMC(3,NGL,N,K))))
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
! M=1:labile, 2, resistant
  CGOMZ=TFNG(NGL,N,K)*OMGR*AMAX1(0.0_r8,OMC(3,NGL,N,K))
  DO 745 M=1,2
    CGOMS(M,NGL,N,K)=FL(M)*CGOMZ
    IF(OMC(3,NGL,N,K).GT.ZEROS)THEN
      CGONS(M,NGL,N,K)=AMIN1(FL(M)*AMAX1(0.0_r8,OMN(3,NGL,N,K)) &
        ,CGOMS(M,NGL,N,K)*OMN(3,NGL,N,K)/OMC(3,NGL,N,K))
      CGOPS(M,NGL,N,K)=AMIN1(FL(M)*AMAX1(0.0_r8,OMP(3,NGL,N,K)) &
        ,CGOMS(M,NGL,N,K)*OMP(3,NGL,N,K)/OMC(3,NGL,N,K))
    ELSE
      CGONS(M,NGL,N,K)=0.0_r8
      CGOPS(M,NGL,N,K)=0.0_r8
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
    SPOMX=SQRT(TFNG(NGL,N,K))*SPOMC(M)*SPOMK(M)
    RXOMC(M,NGL,N,K)=AMAX1(0.0_r8,OMC(M,NGL,N,K)*SPOMX)
    RXOMN(M,NGL,N,K)=AMAX1(0.0_r8,OMN(M,NGL,N,K)*SPOMX)
    RXOMP(M,NGL,N,K)=AMAX1(0.0_r8,OMP(M,NGL,N,K)*SPOMX)
    RDOMC(M,NGL,N,K)=RXOMC(M,NGL,N,K)*(1.0_r8-RCCC)
    RDOMN(M,NGL,N,K)=RXOMN(M,NGL,N,K)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDOMP(M,NGL,N,K)=RXOMP(M,NGL,N,K)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
    R3OMC(M,NGL,N,K)=RXOMC(M,NGL,N,K)-RDOMC(M,NGL,N,K)
    R3OMN(M,NGL,N,K)=RXOMN(M,NGL,N,K)-RDOMN(M,NGL,N,K)
    R3OMP(M,NGL,N,K)=RXOMP(M,NGL,N,K)-RDOMP(M,NGL,N,K)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!
    RHOMC(M,NGL,N,K)=AMAX1(0.0_r8,RDOMC(M,NGL,N,K)*EHUM)
    RHOMN(M,NGL,N,K)=AMAX1(0.0_r8,RDOMN(M,NGL,N,K)*EHUM)
    RHOMP(M,NGL,N,K)=AMAX1(0.0_r8,RDOMP(M,NGL,N,K)*EHUM)
!
!     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
!
    RCOMC(M,NGL,N,K)=RDOMC(M,NGL,N,K)-RHOMC(M,NGL,N,K)
    RCOMN(M,NGL,N,K)=RDOMN(M,NGL,N,K)-RHOMN(M,NGL,N,K)
    RCOMP(M,NGL,N,K)=RDOMP(M,NGL,N,K)-RHOMP(M,NGL,N,K)
745 CONTINUE
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
    DO 730 M=1,2
      RXMMC(M,NGL,N,K)=AMIN1(OMC(M,NGL,N,K) &
        ,AMAX1(0.0_r8,FRM*RMOMC(M,NGL,N,K)/RCCC))
      RXMMN(M,NGL,N,K)=AMIN1(OMN(M,NGL,N,K) &
        ,AMAX1(0.0_r8,RXMMC(M,NGL,N,K)*CNOMA(NGL,N,K)))
      RXMMP(M,NGL,N,K)=AMIN1(OMP(M,NGL,N,K) &
        ,AMAX1(0.0_r8,RXMMC(M,NGL,N,K)*CPOMA(NGL,N,K)))
      RDMMC(M,NGL,N,K)=RXMMC(M,NGL,N,K)*(1.0-RCCC)
      RDMMN(M,NGL,N,K)=RXMMN(M,NGL,N,K)*(1.0-RCCN)*(1.0-RCCC)
      RDMMP(M,NGL,N,K)=RXMMP(M,NGL,N,K)*(1.0-RCCP)*(1.0-RCCC)
      R3MMC(M,NGL,N,K)=RXMMC(M,NGL,N,K)-RDMMC(M,NGL,N,K)
      R3MMN(M,NGL,N,K)=RXMMN(M,NGL,N,K)-RDMMN(M,NGL,N,K)
      R3MMP(M,NGL,N,K)=RXMMP(M,NGL,N,K)-RDMMP(M,NGL,N,K)
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!     EHUM=humus transfer fraction
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
      RHMMC(M,NGL,N,K)=AMAX1(0.0_r8,RDMMC(M,NGL,N,K)*EHUM)
      RHMMN(M,NGL,N,K)=AMAX1(0.0_r8,RDMMN(M,NGL,N,K)*EHUM)
      RHMMP(M,NGL,N,K)=AMAX1(0.0_r8,RDMMP(M,NGL,N,K)*EHUM)
      RCMMC(M,NGL,N,K)=RDMMC(M,NGL,N,K)-RHMMC(M,NGL,N,K)
      RCMMN(M,NGL,N,K)=RDMMN(M,NGL,N,K)-RHMMN(M,NGL,N,K)
      RCMMP(M,NGL,N,K)=RDMMP(M,NGL,N,K)-RHMMP(M,NGL,N,K)

730 CONTINUE
  ELSE
    DO 720 M=1,2
      RXMMC(M,NGL,N,K)=0.0_r8
      RXMMN(M,NGL,N,K)=0.0_r8
      RXMMP(M,NGL,N,K)=0.0_r8
      RDMMC(M,NGL,N,K)=0.0_r8
      RDMMN(M,NGL,N,K)=0.0_r8
      RDMMP(M,NGL,N,K)=0.0_r8
      R3MMC(M,NGL,N,K)=0.0_r8
      R3MMN(M,NGL,N,K)=0.0_r8
      R3MMP(M,NGL,N,K)=0.0_r8
      RHMMC(M,NGL,N,K)=0.0_r8
      RHMMN(M,NGL,N,K)=0.0_r8
      RHMMP(M,NGL,N,K)=0.0_r8
      RCMMC(M,NGL,N,K)=0.0_r8
      RCMMN(M,NGL,N,K)=0.0_r8
      RCMMP(M,NGL,N,K)=0.0_r8
720 CONTINUE
  ENDIF
  end associate
  end subroutine GetMicrobialAnabolismFlux


end module MicBGCMod
