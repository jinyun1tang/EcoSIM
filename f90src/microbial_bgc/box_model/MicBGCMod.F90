module MicBGCMod
!!
! DESCRIPTION:
! codes to do soil biological transfOMBioResduations
!
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils  , only : endrun,destroy
  use TracerIDMod
  use MicAutoCPLXMod
  use minimathmod, only : safe_adb,AZMAX1
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use NitroDiagTypes
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod, only : micforctype
  use EcoSiMParDataMod, only : micpar
  use MicrobMathFuncMod
  use MicrobMathFuncMod
  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: jcplx,NumMicbFunGroups,jsken,ndbiomcp,nlbiomcp
  integer, pointer :: JGniA(:)
  integer, pointer :: JGnfA(:)
  integer, pointer :: JGnio(:)
  integer, pointer :: JGnfo(:)
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
  NumMicbFunGroups  =micpar%NumMicbFunGroups
  jsken =micpar%jsken
  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp

  JGniA => micpar%JGniA
  JGnfA => micpar%JGnfA
  JGnio => micpar%JGnio
  JGnfo => micpar%JGnfo

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
  call nmicf%Init(jcplx,NumMicbFunGroups)
  call nmics%Init(jcplx,NumMicbFunGroups)
  call ncplxf%Init()
  call ncplxs%Init()
  call naqfdiag%ZeroOut()

  micflx%NetNH4Mineralize_col=0._r8;micflx%NetPO4Mineralize_col=0._r8
! write(*,*)'StageBGCEnvironCondition'
  call StageBGCEnvironCondition(micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)
!
! write(*,*)'MicrobialCatabolism'
  call MicrobialCatabolism(micfor,micstt,micflx,nmicdiag, &
    naqfdiag,nmicf,nmics,ncplxf,ncplxs)
!
!write(*,*)'ChemoDenitrification'
  call ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)
!
!     DECOMPOSITION
!
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
  DO  K=1,KL
    DO  N=1,NumMicbFunGroups
      DO NGL=JGnio(N),JGnfo(N)
        ncplxf%ROQCK(K)=ncplxf%ROQCK(K)+nmicf%ROQCD(NGL,K)
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
  DO K=1,KL
!
    !write(*,*)'DECOMPOSITION OF ORGANIC SUBSTRATES'
!
    call SolidOMDecomposition(K,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs)
!
          !write(*,*)'DOC ADSORPTION - DESORPTION'
!
    call DOMSorption(K,micfor,micstt,nmicf,ncplxf,ncplxs)

  ENDDO
        !write(*,*)'RedistDecompProduct'
  call RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)
!
        !write(*,*)'MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION'

  call MicrobialAnabolicUpdate(micfor,micstt,nmicf)
!
        !write(*,*)'MICROBIAL COLONIZATION OF NEW LITTER'
!
  call MicrobialLitterColonization(KL,micfor,micstt,ncplxf,ncplxs)
!
!     AGGREGATE ALL TRANSFOMBioResduATIONS CALCULATED ABOVE FOR EACH N,K
!
  call AggregateTransfOMBioResduations(micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
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
  integer  :: M,N,NGL,MID1,MID2
  real(r8) :: ORGCL
  real(r8) :: TKSO
  real(r8) :: TOSC,TOSA,TOHC
  real(r8) :: TSRH
!     begin_execution
  associate(                                               &
    CNOMActHeter         => nmics%CNOMActHeter           , &
    CPOMActHeter         => nmics%CPOMActHeter           , &
    OMActHeter           => nmics%OMActHeter             , &
    OMC2                 => nmics%OMC2                   , &
    OMN2                 => nmics%OMN2                   , &
    FOM2                 => nmics%FOM2                   , &
    FCN                  => nmics%FCN                    , &
    FCP                  => nmics%FCP                    , &
    FCNP                 => nmics%FCNP                   , &
    BulkSOM              => ncplxs%BulkSOM               , &
    TOMK                 => ncplxs%TOMK                  , &
    TONK                 => ncplxs%TONK                  , &
    TOPK                 => ncplxs%TOPK                  , &
    FOCA                 => ncplxs%FOCA                  , &
    FOAA                 => ncplxs%FOAA                  , &
    CNQ                  => ncplxs%CNQ                   , &
    CPQ                  => ncplxs%CPQ                   , &
    CDOM                 => ncplxs%CDOM                  , &
    ORCT                 => ncplxs%ORCT                  , &
    OSCT                 => ncplxs%OSCT                  , &
    OSAT                 => ncplxs%OSAT                  , &
    TONX                 => ncplxs%TONX                  , &
    TOPX                 => ncplxs%TOPX                  , &
    TORC                 => nmicdiag%TORC                , &
    TotActMicrobiom      => nmicdiag%TotActMicrobiom     , &
    TotBiomNO2Consumers  => nmicdiag%TotBiomNO2Consumers , &
    XCO2                 => nmicdiag%XCO2                , &
    TSensGrowth          => nmicdiag%TSensGrowth         , &
    TSensMaintR          => nmicdiag%TSensMaintR         , &
    ThetaLitr            => nmicdiag%ThetaLitr           , &
    ThetaZ               => nmicdiag%ThetaZ              , &
    VOLWZ                => nmicdiag%VOLWZ               , &
    ZNH4T                => nmicdiag%ZNH4T               , &
    ZNO3T                => nmicdiag%ZNO3T               , &
    ZNO2T                => nmicdiag%ZNO2T               , &
    H2P4T                => nmicdiag%H2P4T               , &
    H1P4T                => nmicdiag%H1P4T               , &
    CNOMActAutor         => nmics%CNOMActAutor           , &
    CPOMActAutor         => nmics%CPOMActAutor           , &
    OMActAutor           => nmics%OMActAutor             , &
    OMC2ff               => nmics%OMC2ff                 , &
    OMN2ff               => nmics%OMN2ff                 , &
    FOM2ff               => nmics%FOM2ff                 , &
    FCNff                => nmics%FCNff                  , &
    FCPff                => nmics%FCPff                  , &
    FCNPff               => nmics%FCNPff                 , &
    rNCOMCAutor             => micpar%rNCOMCAutor        , &
    rPCOMCAutor             => micpar%rPCOMCAutor        , &
    rNCOMC               => micpar%rNCOMC                , &
    rPCOMC               => micpar%rPCOMC                , &
    FL                   => micpar%FL                    , &
    k_humus              => micpar%k_humus               , &
    k_POM                => micpar%k_POM                 , &
    is_activef_micb      => micpar%is_activef_micb       , &
    n_O2facult_bacter    => micpar%n_O2facult_bacter     , &
    AmmoniaOxidBacter    => micpar%AmmoniaOxidBacter     , &
    litrm                => micfor%litrm                 , &
    VWatLitRHoldCapcity  => micfor%VWatLitRHoldCapcity   , &
    VLWatMicP            => micfor%VLWatMicP             , &
    VOLW0                => micfor%VOLW0                 , &
    THETY                => micfor%THETY                 , &
    VLitR                => micfor%VLitR                 , &
    VLSoilMicP           => micfor%VLSoilMicP            , &
    POROS                => micfor%POROS                 , &
    ZEROS                => micfor%ZEROS                 , &
    FieldCapacity        => micfor%FieldCapacity         , &
    THETW                => micfor%THETW                 , &
    TKS                  => micfor%TKS                   , &
    OFFSET               => micfor%OFFSET                , &
    VLWatMicPM           => micfor%VLWatMicPM            , &
    ZEROS2               => micfor%ZEROS2                , &
    ZERO                 => micfor%ZERO                  , &
    CCO2S                => micstt%CCO2S                 , &
    SolidOM              => micstt%SolidOM               , &
    OSA                  => micstt%OSA                   , &
    OMBioResdu           => micstt%OMBioResdu            , &
    SorbedOM             => micstt%SorbedOM              , &
    OMEheter             => micstt%OMEheter              , &
    DOM                  => micstt%DOM                   , &
    H1PO4                => micstt%H1PO4                 , &
    H1POB                => micstt%H1POB                 , &
    H2PO4                => micstt%H2PO4                 , &
    H2POB                => micstt%H2POB                 , &
    ZNH4B                => micstt%ZNH4B                 , &
    ZNH4S                => micstt%ZNH4S                 , &
    ZNO2B                => micstt%ZNO2B                 , &
    ZNO2S                => micstt%ZNO2S                 , &
    ZNO3B                => micstt%ZNO3B                 , &
    ZNO3S                => micstt%ZNO3S                 , &
    OMEauto              => micstt%OMEauto               , &
    FracBulkSOM          => micstt%FracBulkSOM             &
  )

! get KL, the number of mic-om complexes

!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH OFFSET FOR THERMAL ADAPTATION
  IF(litrm)THEN
    ! surface litter layer
    KL=micpar%NumOfLitrCmplxs
    IF(VWatLitRHoldCapcity.GT.ZEROS2)THEN
      ThetaLitr=VOLW0/VLitR
      ThetaZ=AZMAX1(ThetaLitr-THETY)
      VOLWZ=ThetaZ*VLitR
    ELSE
      VOLWZ=0.0_r8
    ENDIF
  ELSE
!     non-surface layer
    KL=micpar%jcplx
    ThetaZ=AZMAX1((AMIN1(AMAX1(0.5_r8*POROS,FieldCapacity),THETW)-THETY))
    VOLWZ=ThetaZ*VLSoilMicP
  ENDIF


!     TKS=soil temperature
!     OFFSET=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TSensGrowth,TSensMaintR=temperature function for growth,maintenance respiration
!
  TKSO=TKS+OFFSET

  call MicrobPhysTempFun(TKSO, TSensGrowth, TSensMaintR)

!
!     OXYI=inhibition of fermenters by O2
!     ORGCL=SOC used to calculate microbial concentration
!
!  OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS+2.5)))
!  ORGCL=AMIN1(1.0E+05*SoilMicPMassLayer,ORGC)
!
!     TOTAL MINERAL NH4, NO3 AND PO4
!
!     allocate NH4, NO3, HPO4, H2PO4 to non-band and band fractions
!
  ZNH4T=AZMAX1(ZNH4S)+AZMAX1(ZNH4B)
  ZNO3T=AZMAX1(ZNO3S)+AZMAX1(ZNO3B)
  H1P4T=AZMAX1(H1PO4)+AZMAX1(H1POB)
  H2P4T=AZMAX1(H2PO4)+AZMAX1(H2POB)
  ZNO2T=AZMAX1(ZNO2S)+AZMAX1(ZNO2B)
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
  DO  K=1,KL
    OSCT(K)=0.0_r8
    OSAT(K)=0.0_r8

    DO M=1,jsken
      OSCT(K)=OSCT(K)+SolidOM(ielmc,M,K)
      OSAT(K)=OSAT(K)+OSA(M,K)
    enddo
    TOSC=TOSC+OSCT(K)
    TOSA=TOSA+OSAT(K)
  enddo
!
!     TOTAL BIORESIDUE
!
  DO  K=1,KL
    ORCT(K)=0.0_r8
    DO  M=1,ndbiomcp
      ORCT(K)=ORCT(K)+OMBioResdu(ielmc,M,K)
    ENDDO
    TORC=TORC+ORCT(K)
!
!     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
!
!     BulkSOM=total SOC
!
    TOHC=TOHC+SorbedOM(ielmc,K)+SorbedOM(idom_acetate,K)
  enddo

  D860: DO K=1,KL
    BulkSOM(K)=OSAT(K)+ORCT(K)+SorbedOM(ielmc,K)+SorbedOM(idom_acetate,K)
  ENDDO D860
  TSRH=TOSA+TORC+TOHC
!
!     C:N AND C:P RATIOS OF TOTAL BIOMASS
!     CNOMA,CPOMA=N,P contents of active biomass OMA
!     FCN,FCP=effects of N,P limitations on biomass activity
!
  TotActMicrobiom=0.0_r8
  TotBiomNO2Consumers=0.0_r8
  D890: DO K=1,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
! the omb complexes
      D895: DO N=1,NumMicbFunGroups
        DO NGL=JGnio(n),JGnfo(n)
          MID1=micpar%get_micb_id(1,NGL)
          IF(OMEheter(ielmc,MID1,K).GT.ZEROS)THEN
            CNOMActHeter(NGL,K)=AZMAX1(OMEheter(ielmn,MID1,K)/OMEheter(ielmc,MID1,K))
            CPOMActHeter(NGL,K)=AZMAX1(OMEheter(ielmp,MID1,K)/OMEheter(ielmc,MID1,K))
          ELSE
            CNOMActHeter(NGL,K)=rNCOMC(1,NGL,K)
            CPOMActHeter(NGL,K)=rPCOMC(1,NGL,K)
          ENDIF
          OMActHeter(NGL,K)=AZMAX1(OMEheter(ielmc,MID1,K)/FL(1))
          FCN(NGL,K)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CNOMActHeter(NGL,K)/rNCOMC(1,NGL,K))))
          FCP(NGL,K)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CPOMActHeter(NGL,K)/rPCOMC(1,NGL,K))))
          FCNP(NGL,K)=AMIN1(FCN(NGL,K),FCP(NGL,K))

!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
          TotActMicrobiom=TotActMicrobiom+OMActHeter(NGL,K)
          IF(N.EQ.n_O2facult_bacter)THEN
            TotBiomNO2Consumers=TotBiomNO2Consumers+OMActHeter(NGL,K)
          ENDIF
          MID2=micpar%get_micb_id(2,NGL)
          OMC2(NGL,K)=AZMAX1(AMIN1(OMActHeter(NGL,K)*FL(2),OMEheter(ielmc,MID2,K)))
          IF(OMEheter(ielmc,MID2,K).GT.ZEROS)THEN
            FOM2(NGL,K)=AZMAX1(OMC2(NGL,K)/OMEheter(ielmc,MID2,K))
            OMN2(NGL,K)=AZMAX1(FOM2(NGL,K)*OMEheter(ielmn,MID2,K))
          ELSE
            FOM2(NGL,K)=0.0_r8
            OMN2(NGL,K)=0.0_r8
          ENDIF
        ENDDO
      ENDDO D895
    ENDIF
  ENDDO D890

! the abstract complex
  DO N=1,NumMicbFunGroups
    IF(is_activef_micb(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        MID1=micpar%get_micb_id(1,NGL)
        IF(OMEauto(ielmc,MID1).GT.ZEROS)THEN
          CNOMActAutor(NGL)=AZMAX1(OMEauto(ielmn,MID1)/OMEauto(ielmc,MID1))
          CPOMActAutor(NGL)=AZMAX1(OMEauto(ielmp,MID1)/OMEauto(ielmc,MID1))
        ELSE
          CNOMActAutor(NGL)=rNCOMCAutor(1,NGL)
          CPOMActAutor(NGL)=rPCOMCAutor(1,NGL)
        ENDIF
        OMActAutor(NGL)=AZMAX1(OMEauto(ielmc,MID1)/FL(1))
        FCNff(NGL)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CNOMActAutor(NGL)/rNCOMCAutor(1,NGL))))
        FCPff(NGL)=AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(CPOMActAutor(NGL)/rPCOMCAutor(1,NGL))))
        FCNPff(NGL)=AMIN1(FCNff(NGL),FCPff(NGL))
!
!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
        TotActMicrobiom=TotActMicrobiom+OMActAutor(NGL)

        IF(N.EQ.AmmoniaOxidBacter)THEN
          TotBiomNO2Consumers=TotBiomNO2Consumers+OMActAutor(NGL)
        ENDIF
        MID2=micpar%get_micb_id(2,NGL)
        OMC2ff(NGL)=AZMAX1(AMIN1(OMActAutor(NGL)*FL(2),OMEauto(ielmc,MID2)))
        IF(OMEauto(ielmc,MID2).GT.ZEROS)THEN
          FOM2ff(NGL)=AZMAX1(OMC2ff(NGL)/OMEauto(ielmc,MID2))
          OMN2ff(NGL)=AZMAX1(FOM2ff(NGL)*OMEauto(ielmn,MID2))
        ELSE
          FOM2ff(NGL)=0.0_r8
          OMN2ff(NGL)=0.0_r8
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  D690: DO K=1,KL
    TOMK(K)=0.0_r8
    TONK(K)=0.0_r8
    TOPK(K)=0.0_r8
    TONX(K)=0.0_r8
    TOPX(K)=0.0_r8
    D685: DO N=1,NumMicbFunGroups
      DO NGL=JGnio(N),JGnfo(N)
        if(OMActHeter(NGL,K)>ZEROS)THEN
          TOMK(K)=TOMK(K)+OMActHeter(NGL,K)
          TONK(K)=TONK(K)+OMActHeter(NGL,K)*CNOMActHeter(NGL,K)
          TOPK(K)=TOPK(K)+OMActHeter(NGL,K)*CPOMActHeter(NGL,K)
          TONX(K)=TONX(K)+OMActHeter(NGL,K)*rNCOMC(1,NGL,K)   !maximum total N in active micb
          TOPX(K)=TOPX(K)+OMActHeter(NGL,K)*rPCOMC(1,NGL,K)   !maximum total P in active micb
        ENDIF
      ENDDO
    ENDDO D685
  ENDDO D690

  K=jcplx+1
  TOMK(K)=0.0_r8
  TONK(K)=0.0_r8
  TOPK(K)=0.0_r8
  TONX(K)=0.0_r8
  TOPX(K)=0.0_r8  
  DO N=1,NumMicbFunGroups
    DO NGL=JGniA(N),JGnfA(N)
      if(OMActAutor(NGL)>ZEROS)then
        TOMK(K)=TOMK(K)+OMActAutor(NGL)      
        TONK(K)=TONK(K)+OMActAutor(NGL)*CNOMActAutor(NGL)
        TOPK(K)=TOPK(K)+OMActAutor(NGL)*CPOMActAutor(NGL)
        TONX(K)=TONX(K)+OMActAutor(NGL)*rNCOMCAutor(1,NGL)   !maximum total N in active micb
        TOPX(K)=TOPX(K)+OMActAutor(NGL)*rPCOMCAutor(1,NGL)   !maximum total P in active micb
      endif
    ENDDO
  ENDDO

!
!     FracBulkSOM=fraction of total SOC in each substrate complex K
!
  D790: DO K=1,KL
    IF(TSRH.GT.ZEROS)THEN
      FracBulkSOM(K)=BulkSOM(K)/TSRH
    ELSE
      FracBulkSOM(K)=1.0_r8
    ENDIF
    !
    !     DOC CONCENTRATIONS
    !
    !     COQC,COQA=aqueous DOC,acetate concentrations
    !     VLWatMicPM=soil water content, FracBulkSOM=fraction of total SOC
    !     occupied by each substrate complex K
    !
    IF(VLWatMicPM(NPH).GT.ZEROS2)THEN
      IF(FracBulkSOM(K).GT.ZERO)THEN
        CDOM(idom_doc,K)=AZMAX1(DOM(idom_doc,K)/(VLWatMicPM(NPH)*FracBulkSOM(K)))
        CDOM(idom_acetate,K)=AZMAX1(DOM(idom_acetate,K)/(VLWatMicPM(NPH)*FracBulkSOM(K)))
      ELSE
        CDOM(idom_doc,K)=AZMAX1(DOM(idom_doc,K)/VLWatMicPM(NPH))
        CDOM(idom_acetate,K)=AZMAX1(DOM(idom_acetate,K)/VLWatMicPM(NPH))
      ENDIF
    ELSE
      CDOM(idom_doc,K)=0.0_r8
      CDOM(idom_acetate,K)=0.0_r8
    ENDIF
!
!     CNQ,CPQ=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
!
    IF(DOM(idom_doc,K).GT.ZEROS)THEN
      CNQ(K)=AZMAX1(DOM(idom_don,K)/DOM(idom_doc,K))
      CPQ(K)=AZMAX1(DOM(idom_dop,K)/DOM(idom_doc,K))
    ELSE
      CNQ(K)=0.0_r8
      CPQ(K)=0.0_r8
    ENDIF
    IF(DOM(idom_doc,K).GT.ZEROS.AND.DOM(idom_acetate,K).GT.ZEROS)THEN
      FOCA(K)=DOM(idom_doc,K)/(DOM(idom_doc,K)+DOM(idom_acetate,K))
      FOAA(K)=1.0_r8-FOCA(K)
    ELSEIF(DOM(idom_doc,K).GT.ZEROS)THEN
      FOCA(K)=1.0_r8
      FOAA(K)=0.0_r8
    ELSE
      FOCA(K)=0.0_r8
      FOAA(K)=1.0_r8
    ENDIF
  ENDDO D790
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
  integer :: K,M,N,NGL,MID
  real(r8) :: TOMCNK(2)
  REAL(R8) :: OXKX
! begin_execution
  associate(                                             &
    TFNR                => nmics%TFNR,                   &
    TFNG                => nmics%TFNG,                   &
    OMActHeter          => nmics%OMActHeter,             &
    TCGOMEheter         => ncplxf%TCGOMEheter,           &
    XCO2                => nmicdiag%XCO2,                &
    TSensGrowth         => nmicdiag%TSensGrowth,         &
    TotActMicrobiom     => nmicdiag%TotActMicrobiom,     &
    TotBiomNO2Consumers => nmicdiag%TotBiomNO2Consumers, &
    RH2GZ               => nmicdiag%RH2GZ,               &
    WatStressMicb       => nmicdiag%WatStressMicb,       &
    TSensMaintR         => nmicdiag%TSensMaintR,         &
    ZNH4T               => nmicdiag%ZNH4T,               &
    ZNO3T               => nmicdiag%ZNO3T,               &
    ZNO2T               => nmicdiag%ZNO2T,               &
    H2P4T               => nmicdiag%H2P4T,               &
    H1P4T               => nmicdiag%H1P4T,               &
    VOLWZ               => nmicdiag%VOLWZ,               &
    TFNGff              => nmics%TFNGff,                 &
    TFNRff              => nmics%TFNRff,                 &
    OMActAutor          => nmics%OMActAutor,             &
    k_humus             => micpar%k_humus,               &
    k_POM               => micpar%k_POM,                 &
    n_aero_fungi        => micpar%n_aero_fungi,          &
    is_activef_micb     => micpar%is_activef_micb,       &
    PSISoilMatricP      => micfor%PSISoilMatricP,        &
    litrm               => micfor%litrm,                 &
    H1PO4               => micstt%H1PO4,                 &
    H1POB               => micstt%H1POB,                 &
    H2PO4               => micstt%H2PO4,                 &
    H2POB               => micstt%H2POB,                 &
    OMEauto             => micstt%OMEauto,               &
    OMEheter            => micstt%OMEheter               &
  )
  RH2GZ=0.0_r8

  TCGOMEheter(:,:)=0.0_r8

  D760: DO K=1,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      DO  N=1,NumMicbFunGroups
        TOMCNK(:)=0.0_r8
        DO NGL=JGnio(N),JGnfo(N)
          DO M=1,2
            MID=micpar%get_micb_id(M,NGL)          
            TOMCNK(M)=TOMCNK(M)+OMEheter(ielmc,MID,K)
          ENDDO

!           WatStressMicb=water potential (PSISoilMatricP) effect on microbial respiration
!           OXKX=Km for O2 uptake
!           OXKM=Km for heterotrophic O2 uptake set in starts.f
!           TFNG=combined temp and water stress effect on growth respiration
!           TFNR=temperature effect on maintenance respiration
          IF(N.EQ.n_aero_fungi)THEN
            WatStressMicb=EXP(0.1_r8*PSISoilMatricP)
          ELSE
            WatStressMicb=EXP(0.2_r8*PSISoilMatricP)
          ENDIF
          OXKX=OXKM
          TFNG(NGL,K)=TSensGrowth*WatStressMicb
          TFNR(NGL,K)=TSensMaintR
          IF(OMActHeter(NGL,K).GT.0.0_r8)THEN
            call ActiveMicrobes(NGL,N,K,VOLWZ,XCO2,TSensGrowth,WatStressMicb,TOMCNK, &
              OXKX,TotActMicrobiom,TotBiomNO2Consumers,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T, &
              micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO D760

  N=micpar%AmmoniaOxidBacter
  nmicf%RTotNH3OxidSoilAutor=SUM(nmicf%RSOxidSoilAutor(JGniA(N):JGnfA(N)))
  nmicf%RTotNH3OxidBandAutor=SUM(nmicf%RSOxidBandAutor(JGniA(N):JGnfA(N)))
  DO  N=1,NumMicbFunGroups
    IF(is_activef_micb(N))THEN
      TOMCNK(:)=0.0_r8
      DO NGL=JGniA(N),JGnfA(N)
        DO M=1,2
          MID=micpar%get_micb_id(M,NGL)
          TOMCNK(M)=TOMCNK(M)+OMEauto(ielmc,MID)
        ENDDO

        WatStressMicb=EXP(0.2_r8*PSISoilMatricP)
        OXKX=OXKA
        TFNGff(NGL)=TSensGrowth*WatStressMicb
        TFNRff(NGL)=TSensMaintR
        IF(OMActAutor(NGL).GT.0.0_r8)THEN
          call ActiveMicrobAutotrophs(NGL,N,VOLWZ,XCO2,TSensGrowth,WatStressMicb,TOMCNK, &
            OXKX,TotActMicrobiom,TotBiomNO2Consumers,RH2GZ,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T, &
            micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  end associate
  end subroutine MicrobialCatabolism
!------------------------------------------------------------------------------------------
  subroutine ActiveMicrobes(NGL,N,K,VOLWZ,XCO2,TSensGrowth,WatStressMicb,TOMCNK,OXKX,TotActMicrobiom,TotBiomNO2Consumers, &
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,K,N
  real(r8), intent(in) :: VOLWZ
  real(r8), intent(in):: OXKX,tomcnk(2),WatStressMicb,TSensGrowth,XCO2
  real(r8), intent(in) :: TotActMicrobiom,TotBiomNO2Consumers
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
  real(r8) :: RGrowthRespHeter
  real(r8) :: RMaintDefcitAutor
  real(r8) :: RMaintRespHeter
  real(r8) :: SPOMK(2)
  real(r8) :: RMOMK(2)
! begin_execution
  associate(                                                    &
    FracOMActHeter         => nmics%FracOMActHeter            , &
    FracNO2ReduxHeter      => nmics%FracNO2ReduxHeter         , &
    n_O2facult_bacter      => micpar%n_O2facult_bacter        , &    
    FracHeterBiomOfActK                   => nmics%FracHeterBiomOfActK                      , &
    OMActHeter             => nmics%OMActHeter                , &
    RO2UptkHeter           => nmicf%RO2UptkHeter              , &
    RespGrossHeter         => nmicf%RespGrossHeter            , &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter         , &
    RO2DmndHeter           => nmicf%RO2DmndHeter              , &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter         , &
    RNO3ReduxHeterSoil     => nmicf%RNO3ReduxHeterSoil        , &
    RNO3ReduxHeterBand     => nmicf%RNO3ReduxHeterBand        , &
    RNO2ReduxHeterSoil     => nmicf%RNO2ReduxHeterSoil        , &
    RNO2ReduxHeterBand     => nmicf%RNO2ReduxHeterBand        , &
    RN2OReduxHeter         => nmicf%RN2OReduxHeter            , &
    RGOMD                  => nmicf%RGOMD                     , &
    RH2ProdHeter           => nmicf%RH2ProdHeter              , &
    RGOMY                  => nmicf%RGOMY                     , &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter             , &
    RAcettProdHeter        => nmicf%RAcettProdHeter           , &
    RCH4ProdHeter          => nmicf%RCH4ProdHeter             , &
    TOMK                   => ncplxs%TOMK                     , &
    SoilMicPMassLayer      => micfor%SoilMicPMassLayer        , &
    litrm                  => micfor%litrm                    , &
    ORGC                   => micfor%ORGC                     , &
    ZEROS                  => micfor%ZEROS                    , &
    VLSoilPoreMicP_vr      => micfor%VLSoilPoreMicP_vr        , &
    RNO2EcoUptkSoilPrev    => micfor%RNO2EcoUptkSoilPrev        &

  )
! FracOMActHeter,FracNO2ReduxHeter=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActHeter(NGL,K)=OMActHeter(NGL,K)/TotActMicrobiom
  ELSE
    FracOMActHeter(NGL,K)=1.0_r8
  ENDIF

  IF(TotBiomNO2Consumers.GT.ZEROS .and. N.EQ.n_O2facult_bacter)THEN
    FracNO2ReduxHeter(NGL,K)=OMActHeter(NGL,K)/TotBiomNO2Consumers
  ELSE
    FracNO2ReduxHeter(NGL,K)=1.0_r8
  ENDIF
  IF(TOMK(K).GT.ZEROS)THEN
    FracHeterBiomOfActK(NGL,K)=OMActHeter(NGL,K)/TOMK(K)
  ELSE
    FracHeterBiomOfActK(NGL,K)=1.0_r8
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
    D765: DO M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
    ENDDO D765
  ELSE
    D770: DO M=1,2
      SPOMK(M)=1.0_r8
      RMOMK(M)=1.0_r8
    ENDDO D770
  ENDIF
!
! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
! write(*,*)'SubstrateCompetitionFactors'
  call SubstrateCompetitionFactors(NGL,N,K,FOXYX,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA, &
    micfor,naqfdiag,nmicf,nmics,micflx)


  RGOMP=0.0_r8

!
! HETEROTROPHIC BIOMASS RESPIRATION
          !
!  IF(K.LE.jcplx1)THEN
!
!   RESPIRATION BY HETEROTROPHIC AEROBES:
!   N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!   (6)N2 FIXERS
!
  IF(micpar%is_aerobic_hetr(N))THEN
!     write(*,*)'AerobicHeterotrophCatabolism'
    call AerobicHeterotrophCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQC,FOQA, &
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
  ELSEIF(micpar%is_anerobic_hetr(N))THEN
!     write(*,*)'AnaerobCatabolism'
    call AnaerobCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQC,ECHZ,FGOCP,FGOAP,RGOMP, &
      micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx)
!     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
!
!     GOMX=acetate effect on energy yield
!     ECHZ=growth respiration efficiency of aceto. methanogenesis
!
  ELSEIF(N.EQ.micpar%AcetotroMethanogenArchea)THEN
!     write(*,*)'AcetoMethanogenCatabolism'
    call AcetoMethanogenCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQA,ECHZ, &
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
! RO2UptkHeter, RO2DmndHeter=O2-limited, O2-unlimited rates of O2 uptake
! RUPMX=O2-unlimited rate of O2 uptake
! FOXYX=fraction of O2 uptake by N,K relative to total
! dts_gas=1/(NPH*NPT)
! ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
! OLSGL=aqueous O2 diffusivity
! OXYG,OXYS=gaseous, aqueous O2 amounts
! Rain2LitRSurf_col,Irrig2LitRSurf=surface water flux from precipitation, irrigation
! O2_rain_conc,O2_irrig_conc=O2 concentration in Rain2LitRSurf_col,Irrig2LitRSurf
!
  RO2UptkHeter(NGL,K)=0.0_r8
  !aerobic heterotrophs
  IF(micpar%is_aerobic_hetr(N))THEN
!  N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!    (6)N2 FIXERS
!   write(*,*)'AerobicHeterO2Uptake'
    call AerobicHeterO2Uptake(NGL,N,K,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB, &
      micfor,micstt,nmicf,nmics,micflx)
  !anaerboic fertmenting heterotrophs
  ELSEIF(micpar%is_anerobic_hetr(N))THEN
    RespGrossHeter(NGL,K)=RGOMP
    RCO2ProdHeter(NGL,K)=0.333_r8*RespGrossHeter(NGL,K)
    RAcettProdHeter(NGL,K)=0.667_r8*RespGrossHeter(NGL,K)
    RCH4ProdHeter(NGL,K)=0.0_r8
    RO2Uptk4RespHeter(NGL,K)=RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)=0.111_r8*RespGrossHeter(NGL,K)
  ELSEIF(N.EQ.micpar%AcetotroMethanogenArchea)THEN
    RespGrossHeter(NGL,K)=RGOMP
    RCO2ProdHeter(NGL,K)=0.50_r8*RespGrossHeter(NGL,K)
    RAcettProdHeter(NGL,K)=0.0_r8
    RCH4ProdHeter(NGL,K)=0.50_r8*RespGrossHeter(NGL,K)
    RO2Uptk4RespHeter(NGL,K)=RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)=0.0_r8
  ENDIF
!
!  write(*,*)'HETEROTROPHIC DENITRIFICATION'
!
  IF(N.EQ.micpar%n_O2facult_bacter.AND.RO2Dmnd4RespHeter(NGL,K).GT.0.0_r8 &
    .AND.(.not.litrm.OR.VLSoilPoreMicP_vr.GT.ZEROS))THEN

    call HeteroDenitrificCatabolism(NGL,N,K,FOQC,RGOCP, &
      VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  ELSE
    RNO3ReduxHeterSoil(NGL,K)=0.0_r8
    RNO3ReduxHeterBand(NGL,K)=0.0_r8
    RNO2ReduxHeterSoil(NGL,K)=0.0_r8
    RNO2ReduxHeterBand(NGL,K)=0.0_r8
    RN2OReduxHeter(NGL,K)=0.0_r8
    RGOMY(NGL,K)=0.0_r8
    RGOMD(NGL,K)=0.0_r8
  ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
  call BiomassMineralization(NGL,N,K,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt, &
    nmicf,nmics,micflx)
!
  call GatherMicrobialRespiration(NGL,N,K,RMOMK,RGrowthRespHeter,RMaintDefcitAutor,RMaintRespHeter, &
    micfor,micstt,nmicf,nmics)
!
  call GetMicrobialAnabolismFlux(NGL,N,K,ECHZ,FGOCP, &
    FGOAP,RGrowthRespHeter,RMaintDefcitAutor,RMaintRespHeter,spomk,rmomk,micfor,micstt,nmicf, &
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
  real(r8) :: H_1p_conc,CHNO2,CHNOB
  real(r8) :: FNO3S,FNO3B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: VMXC4S,VMXC4B
!     begin_execution
  associate(                                               &
    TSensGrowth         =>  nmicdiag%TSensGrowth         , &
    RNO2ReduxSoilChemo  =>  nmicdiag%RNO2ReduxSoilChemo  , &
    RNO2ReduxBandChemo  =>  nmicdiag%RNO2ReduxBandChemo  , &
    RN2OProdSoilChemo   =>  nmicdiag%RN2OProdSoilChemo   , &
    RN2OProdBandChemo   =>  nmicdiag%RN2OProdBandChemo   , &
    RNO3ProdSoilChemo   =>  nmicdiag%RNO3ProdSoilChemo   , &
    RNO3ProdBandChemo   =>  nmicdiag%RNO3ProdBandChemo   , &
    RNO2ReduxChemo      =>  nmicdiag%RNO2ReduxChemo      , &
    RNO2EcoUptkSoilPrev =>  micfor%RNO2EcoUptkSoilPrev   , &
    RNO2EcoUptkBandPrev =>  micfor%RNO2EcoUptkBandPrev   , &
    ph                  => micfor%pH                     , &
    VLWatMicPM          => micfor%VLWatMicPM             , &
    ZEROS               => micfor%ZEROS                  , &
    ZERO                => micfor%ZERO                   , &
    VLNOB               => micfor%VLNOB                  , &
    VLNO3               => micfor%VLNO3                  , &
    CNO2S               => micstt%CNO2S                  , &
    CNO2B               => micstt%CNO2B                  , &
    ZNO2B               => micstt%ZNO2B                  , &
    ZNO2S               => micstt%ZNO2S                  , &
    RNO2DmndSoilChemo   => micflx%RNO2DmndSoilChemo      , &
    RNO2DmndBandChemo   => micflx%RNO2DmndBandChemo        &
  )
!
!     FNO2,FNB2=fraction of total NO2 demand in non-band,band
!     VMXC4S,VMXC4B=substrate-unlimited NO2 reduction in non-band,band
!     CHNO2,CHNOB=nitrous acid concentration in non-band,band
!     VLWatMicPM=soil water content
!     FNO3S,FNO3B=fractions of NO2 in non-band,band
!     TSensGrowth=temperature stress function
!     RNO2ReduxSoilChemo,RNO2ReduxBandChemo=substrate-limited nitrous acid reduction in non-band,band
!     RN2OProdSoilChemo,RN2OProdBandChemo=N2O production from nitrous acid reduction in non-band,band
!     RNO3ProdSoilChemo,RNO3ProdBandChemo=NO3 production from nitrous acid reduction in non-band,band
!     RNO2ReduxChemo=DON production from nitrous acid reduction
!     RNO2DmndSoilChemo,RNO2DmndBandChemo=demand for NO2 reduction in non-band,band
!     nitrous acid concn CHNO2
  H_1p_conc=AMAX1(ZERO,10.0**(-(PH-3.0_r8)))
  CHNO2=CNO2S*H_1p_conc/0.5_r8
  CHNOB=CNO2B*H_1p_conc/0.5_r8

  IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RNO2DmndSoilChemo/RNO2EcoUptkSoilPrev)
  ELSE
    FNO2=FMN*VLNO3
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2DmndBandChemo/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=FMN*VLNOB
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
  FNO3S=VLNO3
  FNO3B=VLNOB
  VMXC4S=7.5E-02_r8*CHNO2*VLWatMicPM(NPH)*FNO3S*TSensGrowth
  VMXC4B=7.5E-02_r8*CHNOB*VLWatMicPM(NPH)*FNO3B*TSensGrowth
  RNO2ReduxSoilChemo=AZMAX1(AMIN1(ZNO2S*FNO2,VMXC4S))
  RNO2ReduxBandChemo=AZMAX1(AMIN1(ZNO2B*FNB2,VMXC4B))
  RN2OProdSoilChemo=0.10_r8*RNO2ReduxSoilChemo
  RN2OProdBandChemo=0.10_r8*RNO2ReduxBandChemo
  RNO3ProdSoilChemo=0.80_r8*RNO2ReduxSoilChemo
  RNO3ProdBandChemo=0.80_r8*RNO2ReduxBandChemo
  RNO2ReduxChemo=0.10*(RNO2ReduxSoilChemo+RNO2ReduxBandChemo)
  RNO2DmndSoilChemo=VMXC4S
  RNO2DmndBandChemo=VMXC4B

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
  integer  :: K,M,N,KK,NGL,MID,idom,NE
  real(r8) :: OSRT
  real(r8) :: XFRK,XFROM(idom_beg:idom_end)
  real(r8) :: XFME
!     begin_execution
  associate(                         &
    TFNG        => nmics%TFNG      , &
    XOMZ        => nmicf%XOMZ      , &
    ROQCK       => ncplxf%ROQCK    , &
    XOQCK       => ncplxf%XOQCK    , &
    XOQMZ       => ncplxf%XOQMZ    , &
    BulkSOM     => ncplxs%BulkSOM  , &
    TOQCK       => micstt%TOQCK    , &
    DOM         => micstt%DOM      , &
    OMEheter    => micstt%OMEheter , &
    ZEROS       => micfor%ZEROS    , &
    TFND        => micfor%TFND       &
  )
!
!     BulkSOM=total SOC in each K
!     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
!     between each K and KK, FPRIM=priming transfer rate constant
!     TFND=temperature effect on priming transfers
!     ROQCK,OQC,OQN,OQP=respiration,DOC,DON,DOP
!     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
!
  D795: DO K=1,KL
    IF(K.LE.KL-1)THEN
      D800: DO KK=K+1,KL
        OSRT=BulkSOM(K)+BulkSOM(KK)
        IF(BulkSOM(K).GT.ZEROS.AND.BulkSOM(KK).GT.ZEROS)THEN
          XFRK=FPRIM*TFND*(ROQCK(K)*BulkSOM(KK)-ROQCK(KK)*BulkSOM(K))/OSRT
          DO idom=idom_beg,idom_end
            XFROM(idom)=FPRIM*TFND*(DOM(idom,K)*BulkSOM(KK)-DOM(idom,KK)*BulkSOM(K))/OSRT
          ENDDO
          IF(ROQCK(K)+XOQCK(K)-XFRK.GT.0.0_r8 .AND. ROQCK(KK)+XOQCK(KK)+XFRK.GT.0.0_r8)THEN
            XOQCK(K)=XOQCK(K)-XFRK
            XOQCK(KK)=XOQCK(KK)+XFRK
          ENDIF
          DO iDOM=idom_beg,idom_end
            IF(DOM(idom,K)+XOQMZ(idom,K)-XFROM(idom).GT.0.0_r8 &
              .AND.DOM(idom,KK)+XOQMZ(idom,KK)+XFROM(idom).GT.0.0_r8)THEN
              XOQMZ(idom,K)=XOQMZ(idom,K)-XFROM(idom)
              XOQMZ(idom,KK)=XOQMZ(idom,KK)+XFROM(idom)
            ENDIF
          ENDDO
!
!     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
!
!     XFMC,XFMN,XFMP=transfer of microbial C,N,P
!     between each K and KK, FPRIMM=priming transfer rate constant
!     TFNG=temperature+water effect
!     OMC,OMN,OMP=microbial C,N,P
!     BulkSOM=total SOC in each K
!     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
!
          D850: DO N=1,NumMicbFunGroups
            DO  M=1,nlbiomcp
              DO NGL=JGnio(N),JGnfo(N)
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  XFME=FPRIMM*TFNG(NGL,K)*(OMEheter(NE,MID,K)*BulkSOM(KK) &
                    -OMEheter(NE,MID,KK)*BulkSOM(K))/OSRT
                  IF(OMEheter(NE,MID,K)+XOMZ(NE,M,NGL,K)-XFME.GT.0.0_r8 &
                    .AND.OMEheter(NE,MID,KK)+XOMZ(NE,M,NGL,KK)+XFME.GT.0.0_r8)THEN
                    XOMZ(NE,M,NGL,K)=XOMZ(NE,M,NGL,K)-XFME
                    XOMZ(NE,M,NGL,KK)=XOMZ(NE,M,NGL,KK)+XFME
                  ENDIF
                ENDDO
              enddo
            enddo
          ENDDO D850
        ENDIF
      ENDDO D800
    ENDIF
  ENDDO D795
!
!     TRANSFER ALL PRIMING AMONG ALL K
!
!     TOQCK=total respiration of DOC+DOA in soil layer
!     ROQCK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
  TOQCK=0.0_r8
  D840: DO K=1,KL
    ROQCK(K)=ROQCK(K)+XOQCK(K)
    TOQCK=TOQCK+ROQCK(K)
    DO idom=idom_beg,idom_end
      DOM(idom,K)=DOM(idom,K)+XOQMZ(idom,K)
    ENDDO
    DO  N=1,NumMicbFunGroups
      DO  M=1,nlbiomcp
        do NGL=JGnio(N),JGnfo(N)
          MID=micpar%get_micb_id(M,NGL)        
          DO NE=1,NumPlantChemElms
            OMEheter(NE,MID,K)=OMEheter(NE,MID,K)+XOMZ(NE,M,NGL,K)
          ENDDO
        enddo
      enddo
    enddo
  ENDDO D840
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
  real(r8) :: OQEX(idom_beg:idom_end)
  real(r8) :: OHEX(idom_beg:idom_end)
  integer :: idom
  real(r8) :: VLSoilPoreMicP_vrX
  real(r8) :: VLSoilPoreMicP_vrW,VOLCX,VOLCW,VOLAX,VOLAW
!     begin_execution
  associate(                                         &
    CGOMEheter         => nmicf%CGOMEheter         , &
    CGOQC              => nmicf%CGOQC              , &
    CGOAC              => nmicf%CGOAC              , &
    DOMSorp            => ncplxf%DOMSorp           , &
    TCGOMEheter        => ncplxf%TCGOMEheter       , &
    BulkSOM            => ncplxs%BulkSOM           , &
    FOCA               => ncplxs%FOCA              , &
    FOAA               => ncplxs%FOAA              , &
    SoilMicPMassLayer  => micfor%SoilMicPMassLayer , &
    ZERO               => micfor%ZERO              , &
    ZEROS2             => micfor%ZEROS2            , &
    ZEROS              => micfor%ZEROS             , &
    litrm              => micfor%litrm             , &
    VLWatMicPM         => micfor%VLWatMicPM        , &
    FracBulkSOM        => micstt%FracBulkSOM       , &
    DOM                => micstt%DOM               , &
    SorbedOM           => micstt%SorbedOM          , &
    AEC                => micfor%AEC                 &
  )
!     VLWatMicPM=soil water content, FracBulkSOM=fraction of total SOC
!     AEC,AECX=anion exchange capacity
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     TCGOQC,TCGOMEheter,TCGOMEheter,TCGOAC=total uptake of DOC,DON,DOP,acetate
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     TSORP,HSORP=sorption rate constant and coefficient for OHC
!     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
  IF(VLWatMicPM(NPH).GT.ZEROS2.AND.FracBulkSOM(K).GT.ZERO)THEN
    IF(litrm)THEN
      AECX=0.5E+03_r8
    ELSE
      AECX=AEC
    ENDIF
    DO idom=idom_beg,idom_end
      OQEX(idom)=AMAX1(ZEROS,DOM(idom,K)-TCGOMEheter(idom,K))
      OHEX(idom)=AMAX1(ZEROS,SorbedOM(idom,K))
    ENDDO

    VLSoilPoreMicP_vrX=SoilMicPMassLayer*AECX*HSORP*FracBulkSOM(K)
    VLSoilPoreMicP_vrW=VLWatMicPM(NPH)*FracBulkSOM(K)
    IF(FOCA(K).GT.ZERO)THEN
      VOLCX=FOCA(K)*VLSoilPoreMicP_vrX
      VOLCW=FOCA(K)*VLSoilPoreMicP_vrW
      DOMSorp(idom_doc,K)=TSORP*(OQEX(idom_doc)*VOLCX-OHEX(idom_doc)*VOLCW)/(VOLCX+VOLCW)
    ELSE
      DOMSorp(idom_doc,K)=TSORP*(OQEX(idom_doc)*VLSoilPoreMicP_vrX  &
        -OHEX(idom_doc)*VLSoilPoreMicP_vrW)/(VLSoilPoreMicP_vrX+VLSoilPoreMicP_vrW)
    ENDIF

    IF(FOAA(K).GT.ZERO)THEN
      VOLAX=FOAA(K)*VLSoilPoreMicP_vrX
      VOLAW=FOAA(K)*VLSoilPoreMicP_vrW
      DOMSorp(idom_acetate,K)=TSORP*(OQEX(idom_acetate)*VOLAX-OHEX(idom_acetate)*VOLAW)/(VOLAX+VOLAW)
    ELSE
      DOMSorp(idom_acetate,K)=TSORP*(OQEX(idom_acetate)*VLSoilPoreMicP_vrX &
        -OHEX(idom_acetate)*VLSoilPoreMicP_vrW)/(VLSoilPoreMicP_vrX+VLSoilPoreMicP_vrW)
    ENDIF
    DOMSorp(idom_don,K)=TSORP*(OQEX(idom_don)*VLSoilPoreMicP_vrX &
      -OHEX(idom_don)*VLSoilPoreMicP_vrW)/(VLSoilPoreMicP_vrX+VLSoilPoreMicP_vrW)
    DOMSorp(idom_dop,K)=TSORP*(OQEX(idom_dop)*VLSoilPoreMicP_vrX &
      -OHEX(idom_dop)*VLSoilPoreMicP_vrW)/(VLSoilPoreMicP_vrX+VLSoilPoreMicP_vrW)
  ELSE
    DO idom=idom_beg,idom_end
      DOMSorp(idom,K)=0.0_r8
    enddo
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
  integer  :: M,NE,idom
  real(r8) :: CNOMX,CPOMX
  real(r8) :: COQCK,COSC
  real(r8) :: CPR,CNR
  real(r8) :: DCKD
  real(r8) :: DFNS
  real(r8) :: OQCI
  real(r8) :: RHOSCM
  real(r8) :: FCNK(1:jcplx),FCPK(1:jcplx)
  real(r8) :: CNS(4,1:jcplx),CPS(4,1:jcplx)

!     begin_execution
  associate(                                          &
    RHydlysSolidOM     => ncplxf%RHydlysSolidOM     , &
    RHumifySolidOM     => ncplxf%RHumifySolidOM     , &
    RDcmpProdDOM       => ncplxf%RDcmpProdDOM       , &
    RHydlysBioResduOM  => ncplxf%RHydlysBioResduOM  , &
    RHydlysSorptOM     => ncplxf%RHydlysSorptOM     , &
    ROQCK              => ncplxf%ROQCK              , &
    BulkSOM            => ncplxs%BulkSOM            , &
    TOMK               => ncplxs%TOMK               , &
    TONK               => ncplxs%TONK               , &
    TOPK               => ncplxs%TOPK               , &
    CNH                => ncplxs%CNH                , &
    CPH                => ncplxs%CPH                , &
    TONX               => ncplxs%TONX               , &
    TOPX               => ncplxs%TOPX               , &
    VOLWZ              =>  nmicdiag%VOLWZ           , &
    TSensGrowth        =>  nmicdiag%TSensGrowth     , &
    iprotein           => micpar%iprotein           , &
    icarbhyro          => micpar%icarbhyro          , &
    icellulos          => micpar%icellulos          , &
    ilignin            => micpar%ilignin            , &
    SPOSC              => micpar%SPOSC              , &
    k_POM              => micpar%k_POM              , &
    CNRH               => micpar%CNRH               , &
    CPRH               => micpar%CPRH               , &
    CDOM               => ncplxs%CDOM               , &
    EPOC               => micstt%EPOC               , &
    CNOSC              => micstt%CNOSC              , &
    CPOSC              => micstt%CPOSC              , &
    SorbedOM           => micstt%SorbedOM           , &
    SolidOM            => micstt%SolidOM            , &
    OSA                => micstt%OSA                , &
    OMBioResdu         => micstt%OMBioResdu         , &
    VLSoilMicP         => micfor%VLSoilMicP         , &
    ZEROS              => micfor%ZEROS              , &
    ZEROS2             => micfor%ZEROS2             , &
    litrm              => micfor%litrm              , &
    SoilMicPMassLayer  => micfor%SoilMicPMassLayer    &
  )
!     FCPK=N,P limitation to microbial activity in each K
!     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
!     COQCK=aqueous concentration of microbial activity
!     DCKD=Km for decomposition of SOC at current COQCK
!     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
!     DCKI=inhibition of decomposition by microbial concentration
!     BulkSOM=total SOC
!     COSC=concentration of total SOC
!     SoilMicPMassLayer,VLSoilPoreMicP_vr=mass, volume of soil layer
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     OQKI=DOC product inhibition constant for decomposition
!
  IF(TOMK(K).GT.ZEROS)THEN
    CNOMX=TONK(K)/TONX(K)
    CPOMX=TOPK(K)/TOPX(K)
    FCNK(K)=AMIN1(1.0_r8,AMAX1(0.50_r8,CNOMX))
    FCPK(K)=AMIN1(1.0_r8,AMAX1(0.50_r8,CPOMX))
  ELSE
    FCNK(K)=1.0_r8
    FCPK(K)=1.0_r8
  ENDIF
!
!     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
!     CONSTANT FOR DECOMPOSITION
!
  IF(VOLWZ.GT.ZEROS2)THEN
    COQCK=AMIN1(0.1E+06_r8,ROQCK(K)/VOLWZ)
  ELSE
    COQCK=0.1E+06_r8
  ENDIF
  IF(litrm)THEN
    DCKD=DCKM0*(1.0_r8+COQCK/DCKI)
  ELSE
    DCKD=DCKML*(1.0_r8+COQCK/DCKI)
  ENDIF
  IF(BulkSOM(K).GT.ZEROS)THEN
    IF(SoilMicPMassLayer.GT.ZEROS)THEN
      COSC=BulkSOM(K)/SoilMicPMassLayer
    ELSE
      COSC=BulkSOM(K)/VLSoilMicP
    ENDIF
    DFNS=COSC/(COSC+DCKD)
    OQCI=1.0_r8/(1.0_r8+CDOM(idom_doc,K)/OQKI)
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
!     TSensGrowth=temperature stress effect
!     BulkSOM=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
    D785: DO M=1,jsken
      IF(SolidOM(ielmc,M,K).GT.ZEROS)THEN
        CNS(M,K)=AZMAX1(SolidOM(ielmn,M,K)/SolidOM(ielmc,M,K))
        CPS(M,K)=AZMAX1(SolidOM(ielmp,M,K)/SolidOM(ielmc,M,K))
        RHydlysSolidOM(ielmc,M,K)=AZMAX1(AMIN1(0.5_r8*OSA(M,K) &
          ,SPOSC(M,K)*ROQCK(K)*DFNS*OQCI*TSensGrowth*OSA(M,K)/BulkSOM(K)))
        RHydlysSolidOM(ielmn,M,K)=AZMAX1(AMIN1(SolidOM(ielmn,M,K),CNS(M,K)*RHydlysSolidOM(ielmc,M,K)))/FCNK(K)
        RHydlysSolidOM(ielmp,M,K)=AZMAX1(AMIN1(SolidOM(ielmp,M,K),CPS(M,K)*RHydlysSolidOM(ielmc,M,K)))/FCPK(K)

      ELSE
        CNS(M,K)=CNOSC(M,K)
        CPS(M,K)=CPOSC(M,K)
        DO NE=1,NumPlantChemElms
          RHydlysSolidOM(NE,M,K)=0.0_r8
        ENDDO
      ENDIF
    ENDDO D785
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
    IF(K.LE.micpar%NumOfLitrCmplxs)THEN
      RHumifySolidOM(ielmc,ilignin,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmn,ilignin,K)/CNRH(k_POM) &
        ,RHydlysSolidOM(ielmp,ilignin,K)/CPRH(k_POM),EPOC*RHydlysSolidOM(ielmc,ilignin,K)))
      RHOSCM=0.10_r8*RHumifySolidOM(ielmc,ilignin,K)
      RHumifySolidOM(ielmc,iprotein,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,iprotein,K) &
        ,RHydlysSolidOM(ielmn,iprotein,K)/CNRH(k_POM) &
        ,RHydlysSolidOM(ielmp,iprotein,K)/CPRH(k_POM),RHOSCM))
      RHumifySolidOM(ielmc,icarbhyro,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,icarbhyro,K) &
        ,RHydlysSolidOM(ielmn,icarbhyro,K)/CNRH(k_POM) &
        ,RHydlysSolidOM(ielmp,icarbhyro,K)/CPRH(k_POM),RHOSCM))
      RHumifySolidOM(ielmc,icellulos,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,icellulos,K) &
        ,RHydlysSolidOM(ielmn,icellulos,K)/CNRH(k_POM) &
        ,RHydlysSolidOM(ielmp,icellulos,K)/CPRH(k_POM),RHOSCM-RHumifySolidOM(ielmc,icarbhyro,K)))

      D805: DO M=1,jsken
        RHumifySolidOM(ielmn,M,K)=AMIN1(RHydlysSolidOM(ielmn,M,K),RHumifySolidOM(ielmc,M,K)*CNRH(k_POM))
        RHumifySolidOM(ielmp,M,K)=AMIN1(RHydlysSolidOM(ielmp,M,K),RHumifySolidOM(ielmc,M,K)*CPRH(k_POM))
        DO NE=1,NumPlantChemElms
          RDcmpProdDOM(NE,M,K)=RHydlysSolidOM(NE,M,K)-RHumifySolidOM(NE,M,K)
        ENDDO
      ENDDO D805
    ELSE
      D810: DO M=1,jsken
        DO NE=1,NumPlantChemElms      
          RHumifySolidOM(NE,M,K)=0.0_r8
          RDcmpProdDOM(NE,M,K)=RHydlysSolidOM(NE,M,K)
        ENDDO
      ENDDO D810
    ENDIF
  ELSE
    D780: DO M=1,jsken
      DO NE=1,NumPlantChemElms    
        RHydlysSolidOM(NE,M,K)=0.0_r8
        RHumifySolidOM(NE,M,K)=0.0_r8
        RDcmpProdDOM(NE,M,K)=0.0_r8
      ENDDO  
    ENDDO D780
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
!     TSensGrowth=temperature stress effect
!     BulkSOM=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
  IF(BulkSOM(K).GT.ZEROS)THEN
    D775: DO M=1,ndbiomcp
      IF(OMBioResdu(ielmc,M,K).GT.ZEROS)THEN
        CNR=AZMAX1(OMBioResdu(ielmn,M,K)/OMBioResdu(ielmc,M,K))
        CPR=AZMAX1(OMBioResdu(ielmp,M,K)/OMBioResdu(ielmc,M,K))
        RHydlysBioResduOM(ielmc,M,K)=AZMAX1(AMIN1(OMBioResdu(ielmc,M,K) &
          ,SPORC(M)*ROQCK(K)*DFNS*OQCI*TSensGrowth*OMBioResdu(ielmc,M,K)/BulkSOM(K)))
    !    3*AMIN1(FCNK(K),FCPK(K))
        RHydlysBioResduOM(ielmn,M,K)=AZMAX1(AMIN1(OMBioResdu(ielmn,M,K),CNR*RHydlysBioResduOM(ielmc,M,K)))/FCNK(K)
        RHydlysBioResduOM(ielmp,M,K)=AZMAX1(AMIN1(OMBioResdu(ielmp,M,K),CPR*RHydlysBioResduOM(ielmc,M,K)))/FCPK(K)
      ELSE
        DO NE=1,NumPlantChemElms      
          RHydlysBioResduOM(NE,M,K)=0.0_r8
        ENDDO  
      ENDIF
    ENDDO D775
  ELSE
    D776: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms    
        RHydlysBioResduOM(NE,M,K)=0.0_r8
      ENDDO  
    ENDDO D776
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
!     TSensGrowth=temperature stress effect
!     BulkSOM=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
  IF(BulkSOM(K).GT.ZEROS)THEN
    IF(SorbedOM(ielmc,K).GT.ZEROS)THEN
      CNH(K)=AZMAX1(SorbedOM(ielmn,K)/SorbedOM(ielmc,K))
      CPH(K)=AZMAX1(SorbedOM(ielmp,K)/SorbedOM(ielmc,K))
      RHydlysSorptOM(ielmc,K)=AZMAX1(AMIN1(SorbedOM(ielmc,K) &
        ,SPOHC*ROQCK(K)*DFNS*OQCI*TSensGrowth*SorbedOM(ielmc,K)/BulkSOM(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      RHydlysSorptOM(ielmn,K)=AZMAX1(AMIN1(SorbedOM(ielmn,K),CNH(K)*RHydlysSorptOM(ielmc,K)))/FCNK(K)
      RHydlysSorptOM(ielmp,K)=AZMAX1(AMIN1(SorbedOM(ielmp,K),CPH(K)*RHydlysSorptOM(ielmc,K)))/FCPK(K)
      RHydlysSorptOM(idom_acetate,K)=AZMAX1(AMIN1(SorbedOM(idom_acetate,K) &
        ,SPOHA*ROQCK(K)*DFNS*TSensGrowth*SorbedOM(idom_acetate,K)/BulkSOM(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
    ELSE
      CNH(K)=0.0_r8
      CPH(K)=0.0_r8
      DO idom=idom_beg,idom_end
        RHydlysSorptOM(idom,K)=0.0_r8
      ENDDO
    ENDIF
  ELSE
    CNH(K)=0.0_r8
    CPH(K)=0.0_r8
    DO idom=idom_beg,idom_end          
      RHydlysSorptOM(idom,K)=0.0_r8
    ENDDO
  ENDIF
  end associate
  end subroutine SolidOMDecomposition
!------------------------------------------------------------------------------------------

  subroutine RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicDiagType), intent(in) :: nmicdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout):: ncplxs
  integer :: K,M,N,NGL,NE,idom
  real(r8) :: FORC(0:jcplx)
!     begin_execution
  associate(                                             &
    k_POM                  => micpar%k_POM             , &
    CGOMEheter             => nmicf%CGOMEheter         , &
    CGOQC                  => nmicf%CGOQC              , &
    CGOAC                  => nmicf%CGOAC              , &
    RCOMEheter             => nmicf%RCOMEheter         , &
    RCMMEheter             => nmicf%RCMMEheter         , &
    RCCMEheter             => nmicf%RCCMEheter         , &
    RAcettProdHeter        => nmicf%RAcettProdHeter    , &
    RHydlysSolidOM         => ncplxf%RHydlysSolidOM    , &
    RHumifySolidOM         => ncplxf%RHumifySolidOM    , &
    RDcmpProdDOM           => ncplxf%RDcmpProdDOM      , &
    RHydlysBioResduOM      => ncplxf%RHydlysBioResduOM , &
    RHydlysSorptOM         => ncplxf%RHydlysSorptOM    , &
    DOMSorp                => ncplxf%DOMSorp           , &
    ORCT                   => ncplxs%ORCT              , &
    RNO2ReduxChemo         =>  nmicdiag%RNO2ReduxChemo , &
    TORC                   =>  nmicdiag%TORC           , &
    RCMMEautor             => nmicf%RCMMEautor         , &
    RCOMEautor             => nmicf%RCOMEautor         , &
    SolidOM                => micstt%SolidOM           , &
    iprotein               => micpar%iprotein          , &
    OSA                    => micstt%OSA               , &
    OSC13U                 => micstt%OSC13U            , &
    OSN13U                 => micstt%OSN13U            , &
    OSP13U                 => micstt%OSP13U            , &
    DOM                    => micstt%DOM               , &
    OMBioResdu             => micstt%OMBioResdu        , &
    SorbedOM               => micstt%SorbedOM          , &
    ZEROS                  => micfor%ZEROS             , &
    Litrm                  => micfor%litrm               &
  )
!
!     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
!     HETEROTROPHIC SUBSTRATE-MICROBE complexES
!
!     FORC=fraction of total microbial residue
!     ORCT=microbial residue
!     RCCMEheter,RCCMN,RCCMP=transfer of auto LitrFall C,N,P to each hetero K
!     RCOMEheter,RCOMEheter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!     RCMMEheter,RCMMN,RCMMEheter=transfer of senesence LitrFall C,N,P to residue
!
  D1690: DO K=1,KL
    IF(TORC.GT.ZEROS)THEN
      FORC(K)=ORCT(K)/TORC
    ELSE
      IF(K.EQ.k_POM)THEN
        FORC(K)=1.0_r8
      ELSE
        FORC(K)=0.0_r8
      ENDIF
    ENDIF
    D1685: DO N=1,NumMicbFunGroups
      D1680: DO M=1,ndbiomcp
        DO NGL=JGniA(N),JGnfA(N)
        DO NE=1,NumPlantChemElms
          RCCMEheter(NE,M,NGL,K)=(RCOMEautor(NE,M,NGL)+RCMMEautor(NE,M,NGL))*FORC(K)
          ENDDO
        ENDDO
      ENDDO D1680
    ENDDO D1685
  ENDDO D1690
!
!   REDISTRIBUTE C,N AND P TRANSFOMBioResduATIONS AMONG STATE
!   VARIABLES IN SUBSTRATE-MICROBE complexES
!

  D590: DO K=1,KL
    D580: DO M=1,jsken
!
!     SUBSTRATE DECOMPOSITION PRODUCTS
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OQC,OQN,OQP,OQA=DOC,DON,DOP
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
      DO NE=1,NumPlantChemElms
        SolidOM(NE,M,K)=SolidOM(NE,M,K)-RHydlysSolidOM(NE,M,K)
      ENDDO

!     OSA(M,K)=OSA(M,K)-RHydlysSolidOM(ielmc,M,K)
      DO NE=1,NumPlantChemElms
      DOM(NE,K)=DOM(NE,K)+RDcmpProdDOM(NE,M,K)
      ENDDO
!
!     LIGNIFICATION PRODUCTS
!
!       RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!
      IF(.not.litrm)THEN
! add to POM carbonhydrate
 !      OSA(1,k_POM)=OSA(1,k_POM)+RHumifySolidOM(ielmc,M,K)
        DO NE=1,NumPlantChemElms
          SolidOM(NE,iprotein,k_POM)=SolidOM(NE,iprotein,k_POM)+RHumifySolidOM(NE,M,K)
        ENDDO
      ELSE
 !      OSA13U=OSA13U+RHumifySolidOM(ielmc,M,K)
        OSC13U=OSC13U+RHumifySolidOM(ielmc,M,K)
        OSN13U=OSN13U+RHumifySolidOM(ielmn,M,K)
        OSP13U=OSP13U+RHumifySolidOM(ielmp,M,K)
      ENDIF

    ENDDO D580
!
!     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     RNO2ReduxChemo=DON production from nitrous acid reduction
!
    D575: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms    
        OMBioResdu(NE,M,K)=OMBioResdu(NE,M,K)-RHydlysBioResduOM(NE,M,K)
        DOM(NE,K)=DOM(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D575
    DO idom=idom_beg,idom_end
      DOM(idom,K)=DOM(idom,K)+RHydlysSorptOM(idom,K)
      SorbedOM(idom,K)=SorbedOM(idom,K)-RHydlysSorptOM(idom,K)
    ENDDO
!
!     MICROBIAL UPTAKE OF DISSOLVED C, N, P
!
!     CGOQC,CGOAC,CGOMEheter,CGOMEheter=DOC,acetate,DON,DOP uptake
!     RAcettProdHeter=acetate production from fermentation
!
    D570: DO N=1,NumMicbFunGroups
      DO NGL=JGnio(N),JGnfo(N)
        DOM(idom_doc,K)=DOM(idom_doc,K)-CGOQC(NGL,K)
        DOM(idom_don,K)=DOM(idom_don,K)-CGOMEheter(ielmp,NGL,K)
        DOM(idom_dop,K)=DOM(idom_dop,K)-CGOMEheter(ielmp,NGL,K)
        DOM(idom_acetate,K)=DOM(idom_acetate,K)-CGOAC(NGL,K)+RAcettProdHeter(NGL,K)
!
!     MICROBIAL DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RCOMEheter,RCOMEheter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!     RCCMEheter,RCCMN,RCCMP=transfer of auto LitrFall C,N,P to each hetero K
!     RCMMEheter,RCMMN,RCMMEheter=transfer of senesence LitrFall C,N,P to residue
!
        D565: DO M=1,ndbiomcp
          DO NE=1,NumPlantChemElms
            OMBioResdu(NE,M,K)=OMBioResdu(NE,M,K)+RCOMEheter(NE,M,NGL,K)+RCCMEheter(NE,M,NGL,K)+RCMMEheter(NE,M,NGL,K)
          ENDDO
        ENDDO D565
      enddo
    ENDDO D570
!
!     SORPTION PRODUCTS
!
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
    DO idom=idom_beg,idom_end
      DOM(idom,K)=DOM(idom,K)-DOMSorp(idom,K)    
      SorbedOM(idom,K)=SorbedOM(idom,K)+DOMSorp(idom,K)
    ENDDO

  ENDDO D590
  end associate
  end subroutine RedistDecompProduct
!------------------------------------------------------------------------------------------

  subroutine MicrobialAnabolicUpdate(micfor,micstt,nmicf)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicFluxType), intent(inout) :: nmicf
  integer  :: K,M,N,NGL,MID3,MID,NE
  real(r8) ::CGROMC
!     begin_execution
  associate(                                                  &
    CGOMEheter             => nmicf%CGOMEheter              , &
    CGOMES                 => nmicf%CGOMES                  , &
    Resp4NFixHeter         => nmicf%Resp4NFixHeter          , &
    RespGrossHeter         => nmicf%RespGrossHeter          , &
    RGOMD                  => nmicf%RGOMD                   , &
    RNO3TransfSoilHeter    => nmicf%RNO3TransfSoilHeter     , &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter           , &
    RH2PO4TransfSoilHeter  => nmicf%RH2PO4TransfSoilHeter   , &
    RNH4TransfBandHeter    => nmicf%RNH4TransfBandHeter     , &
    RNO3TransfBandHeter    => nmicf%RNO3TransfBandHeter     , &
    RH2PO4TransfBandHeter  => nmicf%RH2PO4TransfBandHeter   , &
    RHOMEheter             => nmicf%RHOMEheter              , &
    RHMMEheter             => nmicf%RHMMEheter              , &
    RN2FixHeter            => nmicf%RN2FixHeter             , &
    RXOMEheter             => nmicf%RXOMEheter              , &
    R3OMEheter             => nmicf%R3OMEheter              , &
    RXMMEheter             => nmicf%RXMMEheter              , &
    R3MMEheter             => nmicf%R3MMEheter              , &
    RNH4TransfLitrHeter    => nmicf%RNH4TransfLitrHeter     , &
    RNO3TransfLitrHeter    => nmicf%RNO3TransfLitrHeter     , &
    RH2PO4TransfLitrHeter  => nmicf%RH2PO4TransfLitrHeter   , &
    RH1PO4TransfSoilHeter  => nmicf%RH1PO4TransfSoilHeter   , &
    RH1PO4TransfBandHeter  => nmicf%RH1PO4TransfBandHeter   , &
    RH1PO4TransfLitrHeter  => nmicf%RH1PO4TransfLitrHeter   , &
    RNH4TransfSoilHeter    => nmicf%RNH4TransfSoilHeter     , &
    k_POM                  => micpar%k_POM                  , &
    k_humus                => micpar%k_humus                , &
    OMEheter               => micstt%OMEheter               , &
    SolidOM                => micstt%SolidOM                , &
    OSC14U                 => micstt%OSC14U                 , &
    OSN14U                 => micstt%OSN14U                 , &
    OSP14U                 => micstt%OSP14U                 , &
    OSC24U                 => micstt%OSC24U                 , &
    OSN24U                 => micstt%OSN24U                 , &
    OSP24U                 => micstt%OSP24U                 , &
    CFOMC                  => micfor%CFOMC                  , &
    CFOMCU                 => micfor%CFOMCU                 , &
    icarbhyro              => micpar%icarbhyro              , &
    iprotein               => micpar%iprotein               , &
    Litrm                  => micfor%litrm                    &
  )
!
!     OMC,OMN,OMP=microbial C,N,P
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     RXOMC,RXOMN,RXOMEheter=microbial C,N,P decomposition
!     RXMMEheter,RXMMN,RXMMP=microbial C,N,P loss from senescence
!

  call MicrobAutotrophAnabolicUpdate(micfor,micstt,nmicf)

  D550: DO K=1,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      DO  N=1,NumMicbFunGroups
        DO NGL=JGnio(N),JGnfo(N)
          D540: DO M=1,2
            MID=micpar%get_micb_id(M,NGL)          
            OMEheter(ielmc,MID,K)=OMEheter(ielmc,MID,K)+CGOMES(ielmc,M,NGL,K) &
              -RXOMEheter(ielmc,M,NGL,K)-RXMMEheter(ielmc,M,NGL,K)
            OMEheter(ielmn,MID,K)=OMEheter(ielmn,MID,K)+CGOMES(ielmn,M,NGL,K) &
              -RXOMEheter(ielmn,M,NGL,K)-RXMMEheter(ielmn,M,NGL,K)
            OMEheter(ielmp,MID,K)=OMEheter(ielmp,MID,K)+CGOMES(ielmp,M,NGL,K) &
              -RXOMEheter(ielmp,M,NGL,K)-RXMMEheter(ielmp,M,NGL,K)
!
!     HUMIFICATION PRODUCTS
!
!     CFOMC=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMEheter=transfer of microbial C,N,P LitrFall to humus
!     RHMMEheter,RHMMN,RHMMEheter=transfer of senesence LitrFall C,N,P to humus
!
            IF(.not.litrm)THEN
!add as protein
              DO NE=1,NumPlantChemElms
                SolidOM(NE,iprotein,k_humus)=SolidOM(NE,iprotein,k_humus) &
                  +CFOMC(1)*(RHOMEheter(NE,M,NGL,K)+RHMMEheter(NE,M,NGL,K))
  !add as carbon hydro
                SolidOM(NE,icarbhyro,k_humus)=SolidOM(NE,icarbhyro,k_humus) &
                  +CFOMC(2)*(RHOMEheter(NE,M,NGL,K)+RHMMEheter(NE,M,NGL,K))
              ENDDO
            ELSE
              OSC14U=OSC14U+CFOMCU(1)*(RHOMEheter(ielmc,M,NGL,K)+RHMMEheter(ielmc,M,NGL,K))
              OSN14U=OSN14U+CFOMCU(1)*(RHOMEheter(ielmn,M,NGL,K)+RHMMEheter(ielmn,M,NGL,K))
              OSP14U=OSP14U+CFOMCU(1)*(RHOMEheter(ielmp,M,NGL,K)+RHMMEheter(ielmp,M,NGL,K))
              OSC24U=OSC24U+CFOMC(2)*(RHOMEheter(ielmc,M,NGL,K)+RHMMEheter(ielmc,M,NGL,K))
              OSN24U=OSN24U+CFOMCU(2)*(RHOMEheter(ielmn,M,NGL,K)+RHMMEheter(ielmn,M,NGL,K))
              OSP24U=OSP24U+CFOMCU(2)*(RHOMEheter(ielmp,M,NGL,K)+RHMMEheter(ielmp,M,NGL,K))
            ENDIF
          ENDDO D540

!
!     INPUTS TO NONSTRUCTURAL POOLS
!
!     CGOMEheter=total DOC+acetate uptake
!     RespGrossHeter=total respiration
!     RGOMD=respiration for denitrifcation
!     Resp4NFixHeter=respiration for N2 fixation
!     RCO2ProdHeter=total CO2 emission
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     R3OMC,R3OMN,R3OMEheter=microbial C,N,P recycling
!     R3MMEheter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!     CGOMEheter,CGOMEheter=DON, DOP uptake
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     RH2PO4TransfSoilHeter,RH2PO4TransfBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     RH1PO4TransfSoilHeter,RH1PO4TransfBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
!     RNH4TransfLitrHeter,RNO3TransfLitrHeter =substrate-limited NH4,NO3 mineraln-immobiln
!     RH2PO4TransfLitrHeter,RH1PO4TransfLitrHeter=substrate-limited H2PO4,HPO4 mineraln-immobiln
!
          CGROMC=CGOMEheter(ielmc,NGL,K)-RespGrossHeter(NGL,K)-RGOMD(NGL,K)-Resp4NFixHeter(NGL,K)
          RCO2ProdHeter(NGL,K)=RCO2ProdHeter(NGL,K)+Resp4NFixHeter(NGL,K)
          MID3=micpar%get_micb_id(3,NGL)
          D555: DO M=1,2
            DO NE=1,NumPlantChemElms
              OMEheter(NE,MID3,K)=OMEheter(NE,MID3,K)-CGOMES(NE,M,NGL,K)+R3OMEheter(NE,M,NGL,K)
            ENDDO
            OMEheter(ielmn,MID3,K)=OMEheter(ielmn,MID3,K)+R3MMEheter(ielmn,M,NGL,K)
            OMEheter(ielmp,MID3,K)=OMEheter(ielmp,MID3,K)+R3MMEheter(ielmp,M,NGL,K)

            RCO2ProdHeter(NGL,K)=RCO2ProdHeter(NGL,K)+R3MMEheter(ielmc,M,NGL,K)
          ENDDO D555
          OMEheter(ielmc,MID3,K)=OMEheter(ielmc,MID3,K)+CGROMC
          OMEheter(ielmn,MID3,K)=OMEheter(ielmn,MID3,K)+CGOMEheter(ielmp,NGL,K) &
            +RNH4TransfSoilHeter(NGL,K)+RNH4TransfBandHeter(NGL,K)+RNO3TransfSoilHeter(NGL,K) &
            +RNO3TransfBandHeter(NGL,K)+RN2FixHeter(NGL,K)
          OMEheter(ielmp,MID3,K)=OMEheter(ielmp,MID3,K)+CGOMEheter(ielmp,NGL,K) &
            +RH2PO4TransfSoilHeter(NGL,K)+RH2PO4TransfBandHeter(NGL,K)+RH1PO4TransfSoilHeter(NGL,K) &
            +RH1PO4TransfBandHeter(NGL,K)
          IF(litrm)THEN
            OMEheter(ielmn,MID3,K)=OMEheter(ielmn,MID3,K)+RNH4TransfLitrHeter(NGL,K)+RNO3TransfLitrHeter(NGL,K)
            OMEheter(ielmp,MID3,K)=OMEheter(ielmp,MID3,K)+RH2PO4TransfLitrHeter(NGL,K)+RH1PO4TransfLitrHeter(NGL,K)
          ENDIF
        enddo
      ENDDO
    ENDIF
  ENDDO D550
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
  associate(                   &
    ROQCK   => ncplxf%ROQCK  , &
    OSCT    => ncplxs%OSCT   , &
    OSAT    => ncplxs%OSAT   , &
    ZEROS   => micfor%ZEROS  , &
    OSA     => micstt%OSA    , &
    SolidOM => micstt%SolidOM, &
    DOSA    => micpar%DOSA     &
  )
!     OSCT,OSAT,OSCX=total,colonized,uncolonized SOC
!     OSA,OSC=colonized,total litter
!     DOSA=rate constant for litter colonization
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
  D475: DO K=1,KL
    OSCT(K)=0.0_r8
    OSAT(K)=0.0_r8
    DO  M=1,jsken
      OSCT(K)=OSCT(K)+SolidOM(ielmc,M,K)
      OSAT(K)=OSAT(K)+OSA(M,K)
    enddo
  ENDDO D475

  D480: DO K=1,KL
    IF(OSCT(K).GT.ZEROS)THEN
      DOSAK=DOSA(K)*AZMAX1(ROQCK(K))
      D485: DO M=1,jsken
        OSA(M,K)=AMIN1(SolidOM(ielmc,M,K),OSA(M,K)+DOSAK*SolidOM(ielmc,M,K)/OSCT(K))
      ENDDO D485
    ELSE
      D490: DO M=1,jsken
        OSA(M,K)=AMIN1(SolidOM(ielmc,M,K),OSA(M,K))
      ENDDO D490
    ENDIF
  ENDDO D480
  end associate
  end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

  subroutine AggregateTransfOMBioResduations(micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicDiagType),intent(in) :: nmicdiag
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(micfluxtype), intent(inout) :: micflx
  integer  :: K,M,N,NGL,NE
!     begin_execution
  associate(                                                     &
    CGOMEheter               => nmicf%CGOMEheter               , &
    CGOQC                    => nmicf%CGOQC                    , &
    CGOAC                    => nmicf%CGOAC                    , &
    RO2UptkHeter             => nmicf%RO2UptkHeter             , &
    RNO3ReduxHeterSoil       => nmicf%RNO3ReduxHeterSoil       , &
    RNO3ReduxHeterBand       => nmicf%RNO3ReduxHeterBand       , &
    RNO2ReduxHeterSoil       => nmicf%RNO2ReduxHeterSoil       , &
    RNO2ReduxHeterBand       => nmicf%RNO2ReduxHeterBand       , &
    RN2OReduxHeter           => nmicf%RN2OReduxHeter           , &
    RGOMD                    => nmicf%RGOMD                    , &
    RNH4TransfSoilHeter      => nmicf%RNH4TransfSoilHeter      , &
    RNO3TransfSoilHeter      => nmicf%RNO3TransfSoilHeter      , &
    RH2PO4TransfSoilHeter    => nmicf%RH2PO4TransfSoilHeter    , &
    RNH4TransfBandHeter      => nmicf%RNH4TransfBandHeter      , &
    RNO3TransfBandHeter      => nmicf%RNO3TransfBandHeter      , &
    RH2PO4TransfBandHeter    => nmicf%RH2PO4TransfBandHeter    , &
    RH2ProdHeter             => nmicf%RH2ProdHeter             , &
    RNH4TransfLitrHeter      => nmicf%RNH4TransfLitrHeter      , &
    RNO3TransfLitrHeter      => nmicf%RNO3TransfLitrHeter      , &
    RH2PO4TransfLitrHeter    => nmicf%RH2PO4TransfLitrHeter    , &
    RCO2ProdHeter            => nmicf%RCO2ProdHeter            , &
    RAcettProdHeter          => nmicf%RAcettProdHeter          , &
    RCH4ProdHeter            => nmicf%RCH4ProdHeter            , &
    RSOxidSoilAutor          => nmicf%RSOxidSoilAutor          , &
    RSOxidBandAutor          => nmicf%RSOxidBandAutor          , &
    RH1PO4TransfSoilHeter    => nmicf%RH1PO4TransfSoilHeter    , &
    RH1PO4TransfBandHeter    => nmicf%RH1PO4TransfBandHeter    , &
    RH1PO4TransfLitrHeter    => nmicf%RH1PO4TransfLitrHeter    , &
    RN2FixHeter              => nmicf%RN2FixHeter              , &
    RDcmpProdDOM             => ncplxf%RDcmpProdDOM            , &
    RHydlysBioResduOM        => ncplxf%RHydlysBioResduOM       , &
    RHydlysSorptOM           => ncplxf%RHydlysSorptOM          , &
    DOMSorp                  => ncplxf%DOMSorp                 , &
    TSensGrowth              => nmicdiag%TSensGrowth           , &
    RH2GZ                    => nmicdiag%RH2GZ                 , &
    RNO2ReduxSoilChemo       => nmicdiag%RNO2ReduxSoilChemo    , &
    RNO2ReduxBandChemo       => nmicdiag%RNO2ReduxBandChemo    , &
    RN2OProdSoilChemo        => nmicdiag%RN2OProdSoilChemo     , &
    RN2OProdBandChemo        => nmicdiag%RN2OProdBandChemo     , &
    RNO3ProdSoilChemo        => nmicdiag%RNO3ProdSoilChemo     , &
    RNO3ProdBandChemo        => nmicdiag%RNO3ProdBandChemo     , &
    VOLWZ                    => nmicdiag%VOLWZ                 , &
    CGOMEautor               => nmicf%CGOMEautor               , &
    RNO2ReduxAutorBand       => nmicf%RNO2ReduxAutorBand       , &
    RNO3UptkAutor            => nmicf%RNO3UptkAutor            , &
    RCH4ProdAutor            => nmicf%RCH4ProdAutor            , &
    RNO2ReduxAutorSoil       => nmicf%RNO2ReduxAutorSoil       , &
    RO2UptkAutor             => nmicf%RO2UptkAutor             , &
    RGOMDff                  => nmicf%RGOMDff                  , &
    RCO2ProdAutor            => nmicf%RCO2ProdAutor            , &
    RH1PO4TransfLitrAutor    => nmicf%RH1PO4TransfLitrAutor    , &
    RH2PO4TransfLitrAutor    => nmicf%RH2PO4TransfLitrAutor    , &
    RNO3TransfLitrAutor      => nmicf%RNO3TransfLitrAutor      , &
    RNH4TransfLitrAutor      => nmicf%RNH4TransfLitrAutor      , &
    RN2FixAutor              => nmicf%RN2FixAutor              , &
    RH1PO4TransfBandAutor    => nmicf%RH1PO4TransfBandAutor    , &
    RH2PO4TransfBandAutor    => nmicf%RH2PO4TransfBandAutor    , &
    RNO3TransfBandAutor      => nmicf%RNO3TransfBandAutor      , &
    RNH4TransfBandAutor      => nmicf%RNH4TransfBandAutor      , &
    RH2PO4TransfSoilAutor    => nmicf%RH2PO4TransfSoilAutor    , &
    RNH4TransfSoilAutor      => nmicf%RNH4TransfSoilAutor      , &
    RNO3TransfSoilAutor      => nmicf%RNO3TransfSoilAutor      , &
    RH1PO4TransfSoilAutor    => nmicf%RH1PO4TransfSoilAutor    , &
    TFNQ                     => micstt%TFNQ                    , &
    VOLQ                     => micstt%VOLQ                    , &
    litrm                    => micfor%litrm                   , &
    Lsurf                    => micfor%Lsurf                   , &
    k_POM                    => micpar%k_POM                   , &
    k_humus                  => micpar%k_humus                 , &
    AmmoniaOxidBacter        => micpar%AmmoniaOxidBacter       , &
    AerobicMethanotrofBacter => micpar%AerobicMethanotrofBacter, &
    NitriteOxidBacter        => micpar%NitriteOxidBacter       , &
    is_activef_micb          => micpar%is_activef_micb         , &
    RCH4O                    => micflx%RCH4O                   , &
    RCO2O                    => micflx%RCO2O                   , &
    RH2GO                    => micflx%RH2GO                   , &
    RN2G                     => micflx%RN2G                    , &
    RN2O                     => micflx%RN2O                    , &
    RO2UptkMicb              => micflx%RO2UptkMicb             , &
    XH1BS                    => micflx%XH1BS                   , &
    RH1PO4MicbTransf_vr      => micflx%RH1PO4MicbTransf_vr     , &
    XH2BS                    => micflx%XH2BS                   , &
    RH2PO4MicbTransf_vr      => micflx%RH2PO4MicbTransf_vr     , &
    XN2GS                    => micflx%XN2GS                   , &
    XNH4B                    => micflx%XNH4B                   , &
    RNH4MicbTransf_vr        => micflx%RNH4MicbTransf_vr       , &
    XNO2B                    => micflx%XNO2B                   , &
    RNO2MicbTransf_vr        => micflx%RNO2MicbTransf_vr       , &
    XNO3B                    => micflx%XNO3B                   , &
    RNO3MicbTransf_vr        => micflx%RNO3MicbTransf_vr       , &
    RDOM_micb_flx            => micflx%RDOM_micb_flx             &
  )
  D650: DO K=1,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      DO N=1,NumMicbFunGroups
        DO NGL=JGnio(N),JGnfo(N)
          naqfdiag%TRINH=naqfdiag%TRINH+RNH4TransfSoilHeter(NGL,K)
          naqfdiag%TRINO=naqfdiag%TRINO+RNO3TransfSoilHeter(NGL,K)
          naqfdiag%TRIPO=naqfdiag%TRIPO+RH2PO4TransfSoilHeter(NGL,K)
          naqfdiag%TRIP1=naqfdiag%TRIP1+RH1PO4TransfSoilHeter(NGL,K)
          naqfdiag%TRINB=naqfdiag%TRINB+RNH4TransfBandHeter(NGL,K)
          naqfdiag%TRIOB=naqfdiag%TRIOB+RNO3TransfBandHeter(NGL,K)
          naqfdiag%TRIPB=naqfdiag%TRIPB+RH2PO4TransfBandHeter(NGL,K)
          naqfdiag%TRIB1=naqfdiag%TRIB1+RH1PO4TransfBandHeter(NGL,K)
          naqfdiag%TRN2F=naqfdiag%TRN2F+RN2FixHeter(NGL,K)
          IF(Lsurf)THEN
            naqfdiag%TRINH=naqfdiag%TRINH+RNH4TransfLitrHeter(NGL,K)
            naqfdiag%TRINO=naqfdiag%TRINO+RNO3TransfLitrHeter(NGL,K)
            naqfdiag%TRIPO=naqfdiag%TRIPO+RH2PO4TransfLitrHeter(NGL,K)
            naqfdiag%TRIP1=naqfdiag%TRIP1+RH1PO4TransfLitrHeter(NGL,K)
          ENDIF
          naqfdiag%TRGOM=naqfdiag%TRGOM+RCO2ProdHeter(NGL,K)
          naqfdiag%TRGOC=naqfdiag%TRGOC+RCH4ProdHeter(NGL,K)
          naqfdiag%TRGOD=naqfdiag%TRGOD+RGOMD(NGL,K)
          naqfdiag%TUPOX=naqfdiag%TUPOX+RO2UptkHeter(NGL,K)
          naqfdiag%TRDN3=naqfdiag%TRDN3+RNO3ReduxHeterSoil(NGL,K)
          naqfdiag%TRDNB=naqfdiag%TRDNB+RNO3ReduxHeterBand(NGL,K)
          naqfdiag%TRDN2=naqfdiag%TRDN2+RNO2ReduxHeterSoil(NGL,K)
          naqfdiag%TRD2B=naqfdiag%TRD2B+RNO2ReduxHeterBand(NGL,K)
          naqfdiag%TRDNO=naqfdiag%TRDNO+RN2OReduxHeter(NGL,K)
          naqfdiag%TRGOH=naqfdiag%TRGOH+RH2ProdHeter(NGL,K)
        ENDDO
      ENDDO
    ENDIF
  ENDDO D650

  DO  N=1,NumMicbFunGroups
    IF(is_activef_micb(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%TRINH=naqfdiag%TRINH+RNH4TransfSoilAutor(NGL)
        naqfdiag%TRINO=naqfdiag%TRINO+RNO3TransfSoilAutor(NGL)
        naqfdiag%TRIPO=naqfdiag%TRIPO+RH2PO4TransfSoilAutor(NGL)
        naqfdiag%TRIP1=naqfdiag%TRIP1+RH1PO4TransfSoilAutor(NGL)
        naqfdiag%TRINB=naqfdiag%TRINB+RNH4TransfBandAutor(NGL)
        naqfdiag%TRIOB=naqfdiag%TRIOB+RNO3TransfBandAutor(NGL)
        naqfdiag%TRIPB=naqfdiag%TRIPB+RH2PO4TransfBandAutor(NGL)
        naqfdiag%TRIB1=naqfdiag%TRIB1+RH1PO4TransfBandAutor(NGL)
        naqfdiag%TRN2F=naqfdiag%TRN2F+RN2FixAutor(NGL)
        IF(Lsurf)THEN
          naqfdiag%TRINH=naqfdiag%TRINH+RNH4TransfLitrAutor(NGL)
          naqfdiag%TRINO=naqfdiag%TRINO+RNO3TransfLitrAutor(NGL)
          naqfdiag%TRIPO=naqfdiag%TRIPO+RH2PO4TransfLitrAutor(NGL)
          naqfdiag%TRIP1=naqfdiag%TRIP1+RH1PO4TransfLitrAutor(NGL)
        ENDIF
        naqfdiag%TRGOM=naqfdiag%TRGOM+RCO2ProdAutor(NGL)
        naqfdiag%TRGOC=naqfdiag%TRGOC+RCH4ProdAutor(NGL)
        naqfdiag%TRGOD=naqfdiag%TRGOD+RGOMDff(NGL)
        naqfdiag%TUPOX=naqfdiag%TUPOX+RO2UptkAutor(NGL)
        naqfdiag%TRDN3=naqfdiag%TRDN3+RNO3UptkAutor(NGL)
        naqfdiag%TRDN2=naqfdiag%TRDN2+RNO2ReduxAutorSoil(NGL)
        naqfdiag%TRD2B=naqfdiag%TRD2B+RNO2ReduxAutorBand(NGL)
      ENDDO
    ENDIF
  ENDDO

!     TRGOA=total CO2 uptake by autotrophs, ammonia oxidizer
!  nitrite oxidizer, and hydrogenotophic methanogens
  D645: DO N=1,NumMicbFunGroups
    IF(micpar%is_CO2_autotroph(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%TRGOA=naqfdiag%TRGOA+CGOMEautor(ielmc,NGL)
      ENDDO
    ENDIF
  ENDDO D645
!
!     ALLOCATE AGGREGATED TRANSFOMBioResduATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2O=net CO2 uptake
!     TRGOM total CO2 emission by heterotrophs reducing O2
!     TRGOD=total CO2 emission by denitrifiers reducing NOx
!     RSOxidSoilAutor(3)=CH4 oxidation
!     RCH4O=net CH4 uptake
!     CGOMEheter=total CH4 uptake by autotrophs
!     TRGOC=total CH4 emission
!     RH2GO=net H2 uptake
!     RH2GZ,TRGOH=total H2 uptake, emission
!     RO2UptkMicb,TUPOX=total O2 uptake
!     RN2G=total N2 production
!     TRDNO=total N2O reduction
!     RN2O=total N2O uptake
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RN2OProdSoilChemo,RN2OProdBandChemo=nitrous acid reduction in non-band,band
!
!  print*,'rco2o',naqfdiag%TRGOA,naqfdiag%TRGOM,naqfdiag%TRGOD
  RCO2O=naqfdiag%TRGOA-naqfdiag%TRGOM-naqfdiag%TRGOD
  RCH4O=-naqfdiag%TRGOC

  DO NGL=JGniA(AerobicMethanotrofBacter),JGnfA(AerobicMethanotrofBacter)
    RCO2O=RCO2O-RSOxidSoilAutor(NGL)
    RCH4O=RCH4O+RSOxidSoilAutor(NGL)+CGOMEautor(ielmc,NGL)
  ENDDO
  RH2GO =RH2GZ-naqfdiag%TRGOH
  RO2UptkMicb=naqfdiag%TUPOX
  RN2G  =-naqfdiag%TRDNO
  RN2O  =-naqfdiag%TRDN2-naqfdiag%TRD2B-RN2OProdSoilChemo-RN2OProdBandChemo+naqfdiag%TRDNO
!
  D655: DO K=1,jcplx
    D660: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        RDOM_micb_flx(NE,K)=RDOM_micb_flx(NE,K)+RDcmpProdDOM(NE,M,K)
      ENDDO
    ENDDO D660

    D665: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        RDOM_micb_flx(NE,K)=RDOM_micb_flx(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D665
    DO NE=1,NumPlantChemElms
      RDOM_micb_flx(NE,K)=RDOM_micb_flx(NE,K)+RHydlysSorptOM(NE,K)
    ENDDO
    RDOM_micb_flx(idom_acetate,K)=RDOM_micb_flx(idom_acetate,K)+RHydlysSorptOM(idom_acetate,K)
    D670: DO N=1,NumMicbFunGroups
      DO NGL=JGnio(N),JGnfo(N)
        RDOM_micb_flx(idom_doc,K)=RDOM_micb_flx(idom_doc,K)-CGOQC(NGL,K)
        RDOM_micb_flx(idom_don,K)=RDOM_micb_flx(idom_don,K)-CGOMEheter(ielmp,NGL,K)
        RDOM_micb_flx(idom_dop,K)=RDOM_micb_flx(idom_dop,K)-CGOMEheter(ielmp,NGL,K)
        RDOM_micb_flx(idom_acetate,K)=RDOM_micb_flx(idom_acetate,K)-CGOAC(NGL,K)+RAcettProdHeter(NGL,K)
      ENDDO
    ENDDO D670
    DO NE=1,NumPlantChemElms
      RDOM_micb_flx(NE,K)=RDOM_micb_flx(NE,K)-DOMSorp(NE,K)
    ENDDO
    RDOM_micb_flx(idom_acetate,K)=RDOM_micb_flx(idom_acetate,K)-DOMSorp(idom_acetate,K)
  ENDDO D655
!
!     RNH4MicbTransf_vr,XNH4B=net change in NH4 in band,non-band
!     TRINH,TRINB=total NH4 mineraln-immobn in non-band,band
!     RSOxidSoilAutor(1),RSOxidBandAutor(1)=total NH4 oxidation in non-band,band
!     RNO3MicbTransf_vr,XNO3B=net change in NO3 in band,non-band
!     TRINO,TRIOB=total NO3 immobn in non-band,band
!     RSOxidSoilAutor(2),RSOxidBandAutor(2)=total NO2 oxidation in non-band,band
!     TRDN3,TRDNB=total NO3 reduction in non-band,band
!     RNO3ProdSoilChemo,RNO3ProdBandChemo=NO3 production from nitrous acid reduction in non-band,band
!     RNO2MicbTransf_vr,XNO2B=net change in NO3 in band,non-band
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RNO2ReduxSoilChemo,RNO2ReduxBandChemo=substrate-limited nitrous acid reduction in non-band,band
!     RH2PO4MicbTransf_vr,XH2BS=net change in H2PO4 in band,non-band
!     TRIPO,TRIPB=total H2PO4 mineraln-immobn in non-band,band
!     RH1PO4MicbTransf_vr,XH1BS=net change in HPO4 in band,non-band
!     TRIP1,TRIB1=total HPO4 mineraln-immobn in non-band,band
!     XN2GS=total N2 fixation
!     XZHYS=total H+ production
!     TRN2F=total N2 fixation
!
  RNH4MicbTransf_vr=-naqfdiag%TRINH
  RNO3MicbTransf_vr=-naqfdiag%TRINO-naqfdiag%TRDN3+RNO3ProdSoilChemo
  RNO2MicbTransf_vr=+naqfdiag%TRDN3-naqfdiag%TRDN2-RNO2ReduxSoilChemo
  RH2PO4MicbTransf_vr=-naqfdiag%TRIPO
  RH1PO4MicbTransf_vr=-naqfdiag%TRIP1
  XNH4B=-naqfdiag%TRINB
  XNO3B=-naqfdiag%TRIOB-naqfdiag%TRDNB+RNO3ProdBandChemo
  XNO2B=naqfdiag%TRDNB-naqfdiag%TRD2B-RNO2ReduxBandChemo
  !AmmoniaOxidBacter=1, NitriteOxidBacter=2, AerobicMethanotrofBacter=3
  DO NGL=JGniA(AmmoniaOxidBacter),JGnfA(AmmoniaOxidBacter)
    RNH4MicbTransf_vr=RNH4MicbTransf_vr-RSOxidSoilAutor(NGL)
    RNO2MicbTransf_vr=RNO2MicbTransf_vr+RSOxidSoilAutor(NGL)
    XNH4B=XNH4B-RSOxidBandAutor(NGL)
  ENDDO
  DO NGL=JGniA(NitriteOxidBacter),JGnfA(NitriteOxidBacter)
    RNO3MicbTransf_vr=RNO3MicbTransf_vr+RSOxidSoilAutor(NGL)
    RNO2MicbTransf_vr=RNO2MicbTransf_vr-RSOxidSoilAutor(NGL)
    XNO3B=XNO3B+RSOxidBandAutor(NGL)
    XNO2B=XNO2B-RSOxidBandAutor(NGL)
  ENDDO

  XH2BS=-naqfdiag%TRIPB
  XH1BS=-naqfdiag%TRIB1
  XN2GS=naqfdiag%TRN2F
  TFNQ=TSensGrowth
  VOLQ=VOLWZ
  end associate
  end subroutine AggregateTransfOMBioResduations
!------------------------------------------------------------------------------------------

  subroutine SubstrateCompetitionFactors(NGL,N,K,FOXYX, &
    FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA, &
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
  associate(                                            &
    FracOMActHeter     => nmics%FracOMActHeter        , &
    FracHeterBiomOfActK               => nmics%FracHeterBiomOfActK                  , &
    FNH4XR             => nmicf%FNH4XR                , &
    FNO3XR             => nmicf%FNO3XR                , &
    FP14XR             => nmicf%FP14XR                , &
    FPO4XR             => nmicf%FPO4XR                , &
    RO2DmndHetert      => micflx%RO2DmndHetert        , &
    RINHO              => micflx%RINHO                , &
    RINHB              => micflx%RINHB                , &
    RINOO              => micflx%RINOO                , &
    RINOB              => micflx%RINOB                , &
    RIPOO              => micflx%RIPOO                , &
    RIPBO              => micflx%RIPBO                , &
    RIPO1              => micflx%RIPO1                , &
    RIPB1              => micflx%RIPB1                , &
    ROQCS              => micflx%ROQCS                , &
    ROQAS              => micflx%ROQAS                , &
    RINHOR             => micflx%RINHOR               , &
    RINOOR             => micflx%RINOOR               , &
    RIPOOR             => micflx%RIPOOR               , &
    RIPO1R             => micflx%RIPO1R               , &
    litrm              => micfor%litrm                , &
    VLNH4              => micfor%VLNH4                , &
    VLNHB              => micfor%VLNHB                , &
    VLNOB              => micfor%VLNOB                , &
    VLNO3              => micfor%VLNO3                , &
    VLPOB              => micfor%VLPOB                , &
    VLPO4              => micfor%VLPO4                , &
    ZEROS              => micfor%ZEROS                , &
    ROXYY              => micfor%ROXYY                , &
    RNH4EcoDmndSoilPrev              => micfor%RNH4EcoDmndSoilPrev                , &
    RNH4EcoDmndBandPrev              => micfor%RNH4EcoDmndBandPrev                , &
    RNO3Y              => micfor%RNO3Y                , &
    RNH4EcoDmndLitrPrev             => micfor%RNH4EcoDmndLitrPrev               , &
    RNO3EcoDmndLitrPrev             => micfor%RNO3EcoDmndLitrPrev               , &
    RH1PO4EcoDmndLitrPrev             => micfor%RH1PO4EcoDmndLitrPrev               , &
    RH2PO4EcoDmndLitrPrev             => micfor%RH2PO4EcoDmndLitrPrev               , &
    RN3BY              => micfor%RN3BY                , &
    RPO4Y              => micfor%RPO4Y                , &
    RPOBY              => micfor%RPOBY                , &
    RP14Y              => micfor%RP14Y                , &
    RP1BY              => micfor%RP1BY                , &
    ROQCY              => micfor%ROQCY                , &
    ROQAY              => micfor%ROQAY                , &
    Lsurf              => micfor%Lsurf                , &
    SoilMicPMassLayer0 => micfor%SoilMicPMassLayer0     &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(ROXYY.GT.ZEROS)THEN
    FOXYX=AMAX1(FMN,RO2DmndHetert(NGL,K)/ROXYY)
  ELSE
    FOXYX=AMAX1(FMN,FracOMActHeter(NGL,K))
  ENDIF
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RINHO(NGL,K)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RINHB(NGL,K)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNHB)
  ENDIF
  IF(RNO3Y.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RINOO(NGL,K)/RNO3Y)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RN3BY.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RINOB(NGL,K)/RN3BY)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
  ENDIF
  IF(RPO4Y.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RIPOO(NGL,K)/RPO4Y)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RPOBY.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RIPBO(NGL,K)/RPOBY)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
  ENDIF
  IF(RP14Y.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RIPO1(NGL,K)/RP14Y)
  ELSE
    FP14X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RP1BY.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RIPB1(NGL,K)/RP1BY)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
  ENDIF

  IF(ROQCY(K).GT.ZEROS)THEN
    FOQC=AMAX1(FMN,ROQCS(NGL,K)/ROQCY(K))
  ELSE
    FOQC=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF
  naqfdiag%TFOQC=naqfdiag%TFOQC+FOQC
  IF(ROQAY(K).GT.ZEROS)THEN
    FOQA=AMAX1(FMN,ROQAS(NGL,K)/ROQAY(K))
  ELSE
    FOQA=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
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
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      FNH4XR(NGL,K)=AMAX1(FMN,RINHOR(NGL,K)/RNH4EcoDmndLitrPrev)
    ELSE
      FNH4XR(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      FNO3XR(NGL,K)=AMAX1(FMN,RINOOR(NGL,K)/RNO3EcoDmndLitrPrev)
    ELSE
      FNO3XR(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      FPO4XR(NGL,K)=AMAX1(FMN,RIPOOR(NGL,K)/RH2PO4EcoDmndLitrPrev)
    ELSE
      FPO4XR(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      FP14XR(NGL,K)=AMAX1(FMN,RIPO1R(NGL,K)/RH1PO4EcoDmndLitrPrev)
    ELSE
      FP14XR(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
  ENDIF
  IF(Lsurf.AND.K.NE.micpar%k_POM.AND.K.NE.micpar%k_humus &
    .AND.SoilMicPMassLayer0.GT.ZEROS)THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4XR(NGL,K)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3XR(NGL,K)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4XR(NGL,K)
    naqfdiag%TFP14X=naqfdiag%TFP14X+FP14XR(NGL,K)
  ENDIF
  end associate
  end subroutine SubstrateCompetitionFactors
!------------------------------------------------------------------------------------------

  subroutine AcetoMethanogenCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQA,ECHZ, &
    FGOCP,FGOAP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: FOQA,TSensGrowth,WatStressMicb
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
  associate(                                      &
    FCNP              => nmics%FCNP             , &
    OMActHeter        => nmics%OMActHeter       , &
    RO2DmndHeter      => nmicf%RO2DmndHeter     , &
    RO2Dmnd4RespHeter => nmicf%RO2Dmnd4RespHeter, &
    ROQCD             => nmicf%ROQCD            , &
    CDOM              => ncplxs%CDOM            , &
    DOM               => micstt%DOM             , &
    RO2DmndHetert     => micflx%RO2DmndHetert   , &
    ROQCS             => micflx%ROQCS           , &
    ROQAS             => micflx%ROQAS           , &
    ZERO              => micfor%ZERO            , &
    TKS               => micfor%TKS               &
  )
  GOMX=RGAS*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI))
  GOMM=GOMX/24.0_r8
  ECHZ=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GC4X+GOMM))/EOMH)))
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
!     WatStressMicb=water stress effect, OMA=active biomass
!     TSensGrowth=temp stress effect, FOQA= acetate limitation
!     RGOGX=substrate-limited respiration of acetate
!     RGOGX=competition-limited respiration of acetate
!     OQA=acetate, FOQA=fraction of biological demand for acetate
!     RGOMP=O2-unlimited respiration of acetate
!     ROXY*=O2 demand, ROQCS,ROQCA=DOC, acetate demand
!     ROQCD=microbial respiration used to represent microbial activity
!
  FSBST=CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKAM)
  RGOGY=AZMAX1(FCNP(NGL,K)*VMXM*WatStressMicb*OMActHeter(NGL,K))
  RGOGZ=RGOGY*FSBST*TSensGrowth
  RGOGX=AZMAX1(DOM(idom_acetate,K)*FOQA*ECHZ)
  RGOMP=AMIN1(RGOGX,RGOGZ)
  FGOCP=0.0_r8
  FGOAP=1.0_r8
  RO2Dmnd4RespHeter(NGL,K)=0.0_r8
  RO2DmndHeter(NGL,K)=0.0_r8
  RO2DmndHetert(NGL,K)=0.0_r8
  ROQCS(NGL,K)=0.0_r8
  ROQAS(NGL,K)=RGOGZ
  ROQCD(NGL,K)=0.0_r8
  naqfdiag%TCH4H=naqfdiag%TCH4H+0.5_r8*RGOMP
  end associate
  end subroutine AcetoMethanogenCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobicHeterotrophCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQC,FOQA, &
    ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC,FOQA,WatStressMicb,TSensGrowth
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
  real(r8) :: RO2DmndHetertX,ROQCSX,ROQASX
  real(r8) :: FSBST
!     begin_execution
  associate(                                       &
    OMActHeter        => nmics%OMActHeter,         &
    FCNP              => nmics%FCNP,               &
    RO2Dmnd4RespHeter => nmicf%RO2Dmnd4RespHeter,  &
    RO2DmndHeter      => nmicf%RO2DmndHeter,       &
    ROQCD             => nmicf%ROQCD,              &
    ZEROS             => micfor%ZEROS,             &
    DOM               => micstt%DOM,               &
    n_aero_hetrophb   => micpar%n_aero_hetrophb,   &
    n_O2facult_bacter => micpar%n_O2facult_bacter, &
    n_aero_fungi      => micpar%n_aero_fungi,      &
    n_aero_n2fixer    => micpar%n_aero_n2fixer,    &
    RO2DmndHetert     => micflx%RO2DmndHetert,     &
    ROQCS             => micflx%ROQCS,             &
    ROQAS             => micflx%ROQAS,             &
    FOCA              => ncplxs%FOCA,              &
    FOAA              => ncplxs%FOAA,              &
    CDOM              => ncplxs%CDOM               &
  )
!     ENERGY YIELDS OF O2 REDOX REACTIONS
!     E* = growth respiration efficiency calculated in PARAMETERS
!
  IF(N.EQ.n_aero_hetrophb)THEN
    EO2Q=EO2X
  ELSEIF(N.EQ.n_O2facult_bacter)THEN
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
! WatStressMicb=water stress effect, OMA=active biomass
! TSensGrowth=temp stress effect,FOQC,FOQA=OQC,OQA limitation
! RGOMP=O2-unlimited respiration of DOC+DOA
! RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
!
  FSBSTC=CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
  FSBSTA=CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
  FSBST=FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
  RGOCY=AZMAX1(FCNP(NGL,K)*VMXO*WatStressMicb*OMActHeter(NGL,K))
  RGOCZ=RGOCY*FSBSTC*FOCA(K)*TSensGrowth
  RGOAZ=RGOCY*FSBSTA*FOAA(K)*TSensGrowth
  RGOCX=AZMAX1(DOM(idom_doc,K)*FOQC*EO2Q)
  RGOAX=AZMAX1(DOM(idom_acetate,K)*FOQA*EO2A)
  RGOCP=AMIN1(RGOCX,RGOCZ)
  RGOAP=AMIN1(RGOAX,RGOAZ)
  RGOMP=RGOCP+RGOAP
  IF(RGOMP.GT.ZEROS)THEN
    FGOCP=RGOCP/RGOMP
    FGOAP=RGOAP/RGOMP
  ELSE
    FGOCP=1.0_r8
    FGOAP=0.0_r8
  ENDIF
!
! ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
! BY HETEROTROPHIC AEROBES

! ECHZ=growth respiration yield
! RO2Dmnd4RespHeter,RO2DmndHeter,RO2DmndHetert=O2 demand from DOC,DOA oxidation
! ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
! ROQCD=microbial respiration used to represent microbial activity
!
  ECHZ=EO2Q*FGOCP+EO2A*FGOAP
  RO2Dmnd4RespHeter(NGL,K)=2.667_r8*RGOMP
  RO2DmndHeter(NGL,K)=RO2Dmnd4RespHeter(NGL,K)
  RO2DmndHetertX=RO2DmndHetert(NGL,K)
  ROQCSX=ROQCS(NGL,K)
  ROQASX=ROQAS(NGL,K)
  RO2DmndHetert(NGL,K)=RO2DmndHeter(NGL,K)
  ROQCS(NGL,K)=RGOCZ
  ROQAS(NGL,K)=RGOAZ
  ROQCD(NGL,K)=RGOCY
!
  end associate
  end subroutine AerobicHeterotrophCatabolism
!------------------------------------------------------------------------------------------

  subroutine AnaerobCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQC,ECHZ,FGOCP,FGOAP,RGOMP, &
    micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC,WatStressMicb,TSensGrowth
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
  associate(                                      &
    FCNP              => nmics%FCNP             , &
    OMActHeter        => nmics%OMActHeter       , &
    RO2Dmnd4RespHeter => nmicf%RO2Dmnd4RespHeter, &
    RO2DmndHeter      => nmicf%RO2DmndHeter     , &
    ROQCD             => nmicf%ROQCD            , &
    TKS               => micfor%TKS             , &
    ZERO              => micfor%ZERO            , &
    DOM               => micstt%DOM             , &
    CH2GS             => micstt%CH2GS           , &
    COXYS             => micstt%COXYS           , &
    RO2DmndHetert     => micflx%RO2DmndHetert   , &
    ROQCS             => micflx%ROQCS           , &
    ROQAS             => micflx%ROQAS           , &
    n_anaero_ferm     => micpar%n_anaero_ferm   , &
    CDOM              => ncplxs%CDOM              &
  )
  GH2X=RGAS*1.E-3_r8*TKS*LOG((AMAX1(1.0E-03,CH2GS)/H2KI)**4)
  GH2F=GH2X/72.0
  GOAX=RGAS*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI)**2)
  GOAF=GOAX/72.0
  GHAX=GH2F+GOAF
  IF(N.EQ.n_anaero_ferm)THEN
    ECHZ=AMAX1(EO2X,AMIN1(1.0_r8,1.0/(1.0+AZMAX1((GCHX-GHAX))/EOMF)))
  ELSE
    ECHZ=AMAX1(ENFX,AMIN1(1.0_r8,1.0/(1.0+AZMAX1((GCHX-GHAX))/EOMN)))
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
!     WatStressMicb=water stress effect on respiration
!     OMA=active fermenter biomass
!     TSensGrowth=temp stress effect, FOQC=OQC limitation
!     RFOMP=O2-unlimited respiration of DOC
!     ROQCD=microbial respiration used to represent microbial activity
!
  OXYI=1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*(-COXYS+2.5_r8)))
  FSBST=CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)*OXYI
  RGOFY=AZMAX1(FCNP(NGL,K)*VMXF*WatStressMicb*OMActHeter(NGL,K))
  RGOFZ=RGOFY*FSBST*TSensGrowth
  RGOFX=AZMAX1(DOM(idom_doc,K)*FOQC*ECHZ)
  RGOMP=AMIN1(RGOFX,RGOFZ)
  FGOCP=1.0_r8
  FGOAP=0.0_r8
  RO2Dmnd4RespHeter(NGL,K)=0.0_r8
  RO2DmndHeter(NGL,K)=0.0_r8
  RO2DmndHetert(NGL,K)=0.0_r8
  ROQCS(NGL,K)=RGOFZ
  ROQAS(NGL,K)=0.0_r8
  ROQCD(NGL,K)=RGOFY
  naqfdiag%TRH2G=naqfdiag%TRH2G+RGOMP
!
  end associate
  end subroutine AnaerobCatabolism
!------------------------------------------------------------------------------------------

  subroutine HeteroDenitrificCatabolism(NGL,N,K,FOQC,RGOCP, &
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
  real(r8) :: RNO3UptkSoil,RNO3UptkBand,RDNOX
  real(r8) :: RDNOT
  real(r8) :: RGOM3X,RGOMD3
  real(r8) :: RNO2UptkSoil,RNO2UptkBand,RDN2X,RDN2T,RGOM2X,RGOMD2,RDN2OX,RGOM1X
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
  associate(                                                 &
    WFN                    => nmics%WFN,                     &
    FracOMActHeter         => nmics%FracOMActHeter,          &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,       &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,       &
    RNO3ReduxHeterSoil     => nmicf%RNO3ReduxHeterSoil,      &
    RNO3ReduxHeterBand     => nmicf%RNO3ReduxHeterBand,      &
    RNO2ReduxHeterSoil     => nmicf%RNO2ReduxHeterSoil,      &
    RNO2ReduxHeterBand     => nmicf%RNO2ReduxHeterBand,      &
    RN2OReduxHeter         => nmicf%RN2OReduxHeter,          &
    RGOMD                  => nmicf%RGOMD,                   &
    RGOMY                  => nmicf%RGOMY,                   &
    BulkSOM                => ncplxs%BulkSOM,                &
    RNO2EcoUptkSoilPrev    => micfor%RNO2EcoUptkSoilPrev,    &
    RN2OEcoUptkSoilPrev    => micfor%RN2OEcoUptkSoilPrev,    &
    RNO3Y                  => micfor%RNO3Y,                  &
    VLNO3                  => micfor%VLNO3,                  &
    ZERO                   => micfor%ZERO,                   &
    ZEROS                  => micfor%ZEROS,                  &
    ZEROS2                 => micfor%ZEROS2,                 &
    RNO2EcoUptkBandPrev    => micfor%RNO2EcoUptkBandPrev,    &
    RN3BY                  => micfor%RN3BY,                  &
    VLNOB                  => micfor%VLNOB,                  &
    CNO3B                  => micstt%CNO3B,                  &
    CNO3S                  => micstt%CNO3S,                  &
    CZ2OS                  => micstt%CZ2OS,                  &
    Z2OS                   => micstt%Z2OS,                   &
    ZNO2B                  => micstt%ZNO2B,                  &
    ZNO2S                  => micstt%ZNO2S,                  &
    ZNO3B                  => micstt%ZNO3B,                  &
    ZNO3S                  => micstt%ZNO3S,                  &
    CNO2B                  => micstt%CNO2B,                  &
    CNO2S                  => micstt%CNO2S,                  &
    FracBulkSOM            => micstt%FracBulkSOM,            &
    CH2GS                  => micstt%CH2GS,                  &
    DOM                    => micstt%DOM,                    &
    RVMX3                  => micflx%RVMX3,                  &
    RNO2DmndReduxSoilHeter => micflx%RNO2DmndReduxSoilHeter, &
    RN2ODmndReduxHeter     => micflx%RN2ODmndReduxHeter,     &
    RNO2DmndReduxBandHeter => micflx%RNO2DmndReduxBandHeter, &
    RVMB3                  => micflx%RVMB3                   &
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
    FNO3=AMAX1(FMN,RVMX3(NGL,K)/RNO3Y)
  ELSE
    FNO3=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RN3BY.GT.ZEROS)THEN
    FNB3=AMAX1(FMN,RVMB3(NGL,K)/RN3BY)
  ELSE
    FNB3=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
  ENDIF
  naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3
  naqfdiag%TFNO3B=naqfdiag%TFNO3B+FNB3
!
!     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
!
!     ROXYD=O2 demand RO2Dmnd4RespHeter not met by O2 uptake RO2Uptk4RespHeter
!     VMXD3=demand for NO3-N reduction
!     VMXDXS,VMXDXB=maximum NO3 reduction in non-band, band
!     FNO3S,FNO3B=fractions of total NO3 in non-band, band
!     CNO3S,CNO3B=NO3 concentrations in non-band, band
!     Z3KM,Z2KM=Km for NO3, NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD3S,VMXD3B=substrate-unlimited NO3 reduction in non-band,band
!     OQCD3S,OQCD3B=DOC limitation to NO3 reduction in non-band, band
!     RNO3ReduxHeterSoil,RNO3ReduxHeterBand=substrate-limited NO3 reduction in non-band,band
!     RGOM3X,RGOMD3=substrate-unltd,-ltd respn from NO3 reduction
!     RVMX3,RVMB3=demand for NO3 reduction in non-band,band
!
  ROXYD=AZMAX1(RO2Dmnd4RespHeter(NGL,K)-RO2Uptk4RespHeter(NGL,K))
  VMXD3=0.875_r8*ROXYD
  IF(CNO3S.GT.ZERO)THEN
    VMXDXS=FNO3S*VMXD3*CNO3S/(CNO3S+Z3KM)/(1.0_r8+(CNO2S*Z3KM)/(CNO3S*Z2KM))
  ELSE
    VMXDXS=0.0_r8
  ENDIF
  IF(CNO3B.GT.ZERO)THEN
    VMXDXB=FNO3B*VMXD3*CNO3B/(CNO3B+Z3KM)/(1.0_r8+(CNO2B*Z3KM)/(CNO3B*Z2KM))
  ELSE
    VMXDXB=0.0_r8
  ENDIF
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOM(K).GT.ZERO)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ*FracBulkSOM(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD3S=VMXDXS*FVMXDX
  VMXD3B=VMXDXB*FVMXDX
  OQCZ3=AZMAX1(DOM(idom_doc,K)*FOQC-RGOCP*WFN(NGL,K))
  OQCD3=OQCZ3/eQNO3toOxy
  OQCD3S=OQCD3*FNO3S
  OQCD3B=OQCD3*FNO3B
  ZNO3SX=ZNO3S*FNO3
  ZNO3BX=ZNO3B*FNB3
  RNO3UptkSoil=AZMAX1(AMIN1(ZNO3SX,VMXD3S))
  RNO3UptkBand=AZMAX1(AMIN1(ZNO3BX,VMXD3B))
  RNO3ReduxHeterSoil(NGL,K)=AZMAX1(AMIN1(VMXD3S,OQCD3S,ZNO3SX))
  RNO3ReduxHeterBand(NGL,K)=AZMAX1(AMIN1(VMXD3B,OQCD3B,ZNO3BX))
  RDNOX=RNO3UptkSoil+RNO3UptkBand
  RDNOT=RNO3ReduxHeterSoil(NGL,K)+RNO3ReduxHeterBand(NGL,K)
  RGOM3X=eQNO3toOxy*RDNOX
  RGOMD3=eQNO3toOxy*RDNOT
  RVMX3(NGL,K)=VMXD3S
  RVMB3(NGL,K)=VMXD3B
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
  FNO2S=VLNO3
  FNO2B=VLNOB
  IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RNO2DmndReduxSoilHeter(NGL,K)/RNO2EcoUptkSoilPrev)
  ELSE
    FNO2=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2DmndReduxBandHeter(NGL,K)/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
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
!     RNO2ReduxHeterSoil,RNO2ReduxHeterBand=substrate-limited NO2 reduction in non-band,band
!     RGOM2X,RGOMD2=substrate-unltd,-ltd respn from NO2 reduction
!
  VMXD2=VMXD3-RDNOT
  IF(CNO2S.GT.ZERO)THEN
    VMXDXS=FNO2S*VMXD2*CNO2S/(CNO2S+Z2KM)/(1.0_r8+(CZ2OS*Z2KM)/(CNO2S*Z1KM))
  ELSE
    VMXDXS=0.0_r8
  ENDIF
  IF(CNO2B.GT.ZERO)THEN
    VMXDXB=FNO2B*VMXD2*CNO2B/(CNO2B+Z2KM)/(1.0_r8+(CZ2OS*Z2KM)/(CNO2B*Z1KM))
  ELSE
    VMXDXB=0.0_r8
  ENDIF
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOM(K).GT.ZERO)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ*FracBulkSOM(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD2S=VMXDXS*FVMXDX
  VMXD2B=VMXDXB*FVMXDX
  OQCZ2=AZMAX1(OQCZ3-RGOMD3)
  OQCD2=OQCZ2/eQNO2toOxy
  OQCD2S=OQCD2*FNO3S
  OQCD2B=OQCD2*FNO3B
  ZNO2SX=(ZNO2S+RNO3ReduxHeterSoil(NGL,K))*FNO2
  ZNO2BX=(ZNO2B+RNO3ReduxHeterBand(NGL,K))*FNB2
  RNO2UptkSoil=AZMAX1(AMIN1(ZNO2SX,VMXD2S))
  RNO2UptkBand=AZMAX1(AMIN1(ZNO2BX,VMXD2B))
  RNO2ReduxHeterSoil(NGL,K)=AZMAX1(AMIN1(VMXD2S,OQCD2S,ZNO2SX))
  RNO2ReduxHeterBand(NGL,K)=AZMAX1(AMIN1(VMXD2B,OQCD2B,ZNO2BX))
  RDN2X=RNO2UptkSoil+RNO2UptkBand
  RDN2T=RNO2ReduxHeterSoil(NGL,K)+RNO2ReduxHeterBand(NGL,K)
  RGOM2X=eQNO2toOxy*RDN2X
  RGOMD2=eQNO2toOxy*RDN2T
  RNO2DmndReduxSoilHeter(NGL,K)=VMXD2S
  RNO2DmndReduxBandHeter(NGL,K)=VMXD2B
!
!     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
!     AND ROOT POPULATIONS
!
!     FN2O=fraction of total biological demand for N2O
!
  IF(RN2OEcoUptkSoilPrev.GT.ZEROS)THEN
    FN2O=AMAX1(FMN,RN2ODmndReduxHeter(NGL,K)/RN2OEcoUptkSoilPrev)
  ELSE
    FN2O=AMAX1(FMN,FracOMActHeter(NGL,K))
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
!     RN2ODmndReduxHeter=demand for N2O reduction
!
  VMXD1=(VMXD2-RDN2T)*2.0_r8
  VMXDXS=VMXD1*CZ2OS/(CZ2OS+Z1KM)
  IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOM(K).GT.ZERO)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXS/(VMKI*VOLWZ*FracBulkSOM(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD1S=VMXDXS*FVMXDX
  OQCZ1=AZMAX1(OQCZ2-RGOMD2)
  OQCD1=OQCZ1/eQN2OtoOxy
  Z2OSX=(Z2OS+RDN2T)*FN2O
  RDN2OX=AZMAX1(AMIN1(Z2OSX,VMXD1S))
  RN2OReduxHeter(NGL,K)=AZMAX1(AMIN1(VMXD1S,OQCD1,Z2OSX))
  RGOM1X=eQN2OtoOxy*RDN2OX
  RGOMD1=eQN2OtoOxy*RN2OReduxHeter(NGL,K)
  RGOMY(NGL,K)=RGOM3X+RGOM2X+RGOM1X
  RGOMD(NGL,K)=RGOMD3+RGOMD2+RGOMD1
  RN2ODmndReduxHeter(NGL,K)=VMXD1S
  end associate
  end subroutine HeteroDenitrificCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobicHeterO2Uptake(NGL,N,K,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB, &
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
  associate(                                           &
    WFN                 => nmics%WFN                 , &
    OMActHeter          => nmics%OMActHeter          , &
    RO2UptkHeter        => nmicf%RO2UptkHeter        , &
    RespGrossHeter      => nmicf%RespGrossHeter      , &
    RO2Dmnd4RespHeter   => nmicf%RO2Dmnd4RespHeter   , &
    RO2DmndHeter        => nmicf%RO2DmndHeter        , &
    RO2Uptk4RespHeter   => nmicf%RO2Uptk4RespHeter   , &
    RH2ProdHeter        => nmicf%RH2ProdHeter        , &
    ROQCD               => nmicf%ROQCD               , &
    RCO2ProdHeter       => nmicf%RCO2ProdHeter       , &
    RAcettProdHeter     => nmicf%RAcettProdHeter     , &
    RCH4ProdHeter       => nmicf%RCH4ProdHeter       , &
    RSOxidSoilAutor     => nmicf%RSOxidSoilAutor     , &
    RSOxidBandAutor     => nmicf%RSOxidBandAutor     , &
    O2_irrig_conc       => micfor%O2_irrig_conc      , &
    O2_rain_conc        => micfor%O2_rain_conc       , &
    COXYE               => micfor%COXYE              , &
    ROXYF               => micfor%ROXYF              , &
    ROXYL               => micfor%ROXYL              , &
    Irrig2LitRSurf      => micfor%Irrig2LitRSurf     , &
    Rain2LitRSurf_col   => micfor%Rain2LitRSurf_col  , &
    litrm               => micfor%litrm              , &
    OLSGL               => micfor%OLSGL              , &
    VLSoilPoreMicP_vr   => micfor%VLSoilPoreMicP_vr  , &
    VLSoilMicP          => micfor%VLSoilMicP         , &
    ZERO                => micfor%ZERO               , &
    ZEROS               => micfor%ZEROS              , &
    VLsoiAirPM          => micfor%VLsoiAirPM         , &
    VLWatMicP           => micfor%VLWatMicP          , &
    VLWatMicPM          => micfor%VLWatMicPM         , &
    THETPM              => micfor%THETPM             , &
    DiffusivitySolutEff => micfor%DiffusivitySolutEff, &
    FILM                => micfor%FILM               , &
    TortMicPM           => micfor%TortMicPM          , &
    OXYG                => micstt%OXYG               , &
    OXYS                => micstt%OXYS               , &
    COXYS               => micstt%COXYS              , &
    SOXYL               => micstt%SOXYL              , &
    COXYG               => micstt%COXYG              , &
    RVMX4               => micflx%RVMX4              , &
    RVMB4               => micflx%RVMB4              , &
    RNO2DmndReduxSoilHeter  => micflx%RNO2DmndReduxSoilHeter , &
    RNO2DmndReduxBandHeter               => micflx%RNO2DmndReduxBandHeter              , &
    ROXSK               => micflx%ROXSK                &
  )
  IF(RO2DmndHeter(NGL,K).GT.ZEROS.AND.FOXYX.GT.ZERO)THEN
    IF(.not.litrm.OR.VLSoilPoreMicP_vr.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX=RO2DmndHeter(NGL,K)*dts_gas
      ROXYFX=ROXYF*dts_gas*FOXYX
      OLSGL1=OLSGL*dts_gas
      IF(.not.litrm)THEN
        OXYG1=OXYG*FOXYX
        ROXYLX=ROXYL*dts_gas*FOXYX
      ELSE
        OXYG1=COXYG*VLsoiAirPM(1)*FOXYX
        ROXYLX=(ROXYL+Rain2LitRSurf_col*O2_rain_conc+Irrig2LitRSurf*O2_irrig_conc)*dts_gas*FOXYX
      ENDIF
      OXYS1=OXYS*FOXYX
!
      !write(*,*)'O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP'
!     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
!
      D420: DO M=1,NPH
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
        !write(*,*)'VLSoilMicP=',VLSoilMicP
        THETW1=AZMAX1(safe_adb(VLWatMicPM(M),VLSoilMicP))
        RRADO=ORAD*(FILM(M)+ORAD)/FILM(M)
        DIFOX=TortMicPM(M)*OLSGL1*12.57_r8*BIOS*OMActHeter(NGL,K)*RRADO
        VOLWOX=VLWatMicPM(M)*SOXYL
        VOLPOX=VLsoiAirPM(M)
        VOLWPM=VOLWOX+VOLPOX

        D425: DO MX=1,NPT
          OXYG1=OXYG1+ROXYFX
          OXYS1=OXYS1+ROXYLX
          COXYS1=AMIN1(COXYE*SOXYL,AZMAX1(safe_adb(OXYS1,(VLWatMicPM(M)*FOXYX))))

          !obtain uptake flux
          if(OXYS1<=ZEROS)then
            RMPOX=0.0_r8            
          else
            RMPOX=TranspBasedsubstrateUptake(COXYS1,DIFOX, OXKX, RUPMX, ZEROS)
          endif  
  
          !apply the uptake
          OXYS1=OXYS1-RMPOX
          !apply dissolution-volatilization
          IF(THETPM(M).GT.THETX.AND.VOLPOX.GT.ZEROS)THEN
            ROXDFQ=DiffusivitySolutEff(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
          ELSE
            ROXDFQ=0.0_r8
          ENDIF
          OXYG1=OXYG1-ROXDFQ
          OXYS1=OXYS1+ROXDFQ
          !accumulate upatke
          RO2UptkHeter(NGL,K)=RO2UptkHeter(NGL,K)+RMPOX
          ROXSK(M)=ROXSK(M)+RMPOX
        ENDDO D425
        
      ENDDO D420
      !write(*,*)'420'
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
      !
      !     WFN=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter,RNO2DmndReduxBandHeter=NH3,NO2 oxidation in non-band, band
      !
      WFN(NGL,K)=AMIN1(1.0_r8,AZMAX1(RO2UptkHeter(NGL,K)/RO2DmndHeter(NGL,K)))
!     IF(K.LE.4)THEN
!       ROQCS(NGL,K)=ROQCS(NGL,K)*WFN(NGL,K)
!       ROQAS(NGL,K)=ROQAS(NGL,K)*WFN(NGL,K)
!       ROQCD(NGL,K)=ROQCD(NGL,K)*WFN(NGL,K)
!     ENDIF
    ELSE
      RO2UptkHeter(NGL,K)=RO2DmndHeter(NGL,K)
      WFN(NGL,K)=1.0_r8
    ENDIF
  ELSE
    RO2UptkHeter(NGL,K)=0.0_r8
    WFN(NGL,K)=1.0_r8
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2ProdHeter,RAcettProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
  !     RO2Uptk4RespHeter=O2-limited O2 uptake
  !     RSOxidSoilAutor,RSOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !
  RespGrossHeter(NGL,K)=RGOMP*WFN(NGL,K)
  RCO2ProdHeter(NGL,K)=RespGrossHeter(NGL,K)
  RAcettProdHeter(NGL,K)=0.0_r8
  RCH4ProdHeter(NGL,K)=0.0_r8
  RO2Uptk4RespHeter(NGL,K)=RO2Dmnd4RespHeter(NGL,K)*WFN(NGL,K)
  RH2ProdHeter(NGL,K)=0.0_r8
  !write(*,*)'finish AerobicHeterO2Uptake'
  end associate
  end subroutine AerobicHeterO2Uptake

!------------------------------------------------------------------------------------------

  subroutine BiomassMineralization(NGL,N,K,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
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
  integer :: MID3

!     begin_execution
  associate(                                             &
   TFNG                  => nmics%TFNG                 , &
   OMActHeter            => nmics%OMActHeter           , &
   FNH4XR                => nmicf%FNH4XR               , &
   FNO3XR                => nmicf%FNO3XR               , &
   FPO4XR                => nmicf%FPO4XR               , &
   FP14XR                => nmicf%FP14XR               , &
   RNH4TransfSoilHeter   => nmicf%RNH4TransfSoilHeter  , &
   RNO3TransfSoilHeter   => nmicf%RNO3TransfSoilHeter  , &
   RH2PO4TransfSoilHeter => nmicf%RH2PO4TransfSoilHeter, &
   RNH4TransfBandHeter   => nmicf%RNH4TransfBandHeter  , &
   RNO3TransfBandHeter   => nmicf%RNO3TransfBandHeter  , &
   RH2PO4TransfBandHeter => nmicf%RH2PO4TransfBandHeter, &
   RNH4TransfLitrHeter   => nmicf%RNH4TransfLitrHeter  , &
   RNO3TransfLitrHeter   => nmicf%RNO3TransfLitrHeter  , &
   RH2PO4TransfLitrHeter => nmicf%RH2PO4TransfLitrHeter, &
   RH1PO4TransfSoilHeter => nmicf%RH1PO4TransfSoilHeter, &
   RH1PO4TransfBandHeter => nmicf%RH1PO4TransfBandHeter, &
   RH1PO4TransfLitrHeter => nmicf%RH1PO4TransfLitrHeter, &
   rNCOMC                => micpar%rNCOMC              , &
   rPCOMC                => micpar%rPCOMC              , &
   litrm                 => micfor%litrm               , &
   VLNH4                 => micfor%VLNH4               , &
   VLNHB                 => micfor%VLNHB               , &
   VLNO3                 => micfor%VLNO3               , &
   VLNOB                 => micfor%VLNOB               , &
   VLPOB                 => micfor%VLPOB               , &
   VLWatMicP             => micfor%VLWatMicP           , &
   VOLWU                 => micfor%VOLWU               , &
   VLPO4                 => micfor%VLPO4               , &
   ZNH4B                 => micstt%ZNH4B               , &
   ZNH4S                 => micstt%ZNH4S               , &
   ZNO3B                 => micstt%ZNO3B               , &
   ZNO3S                 => micstt%ZNO3S               , &
   ZNH4TU                => micstt%ZNH4TU              , &
   ZNO3TU                => micstt%ZNO3TU              , &
   H1P4TU                => micstt%H1P4TU              , &
   H2P4TU                => micstt%H2P4TU              , &
   CNH4BU                => micstt%CNH4BU              , &
   CNH4SU                => micstt%CNH4SU              , &
   CH2P4U                => micstt%CH2P4U              , &
   CH2P4BU               => micstt%CH2P4BU             , &
   CH1P4U                => micstt%CH1P4U              , &
   CH1P4BU               => micstt%CH1P4BU             , &
   CH2P4                 => micstt%CH2P4               , &
   CH2P4B                => micstt%CH2P4B              , &
   CNH4B                 => micstt%CNH4B               , &
   CNH4S                 => micstt%CNH4S               , &
   CH1P4                 => micstt%CH1P4               , &
   CH1P4B                => micstt%CH1P4B              , &
   H1PO4                 => micstt%H1PO4               , &
   H1POB                 => micstt%H1POB               , &
   H2PO4                 => micstt%H2PO4               , &
   H2POB                 => micstt%H2POB               , &
   CNO3B                 => micstt%CNO3B               , &
   CNO3S                 => micstt%CNO3S               , &
   CNO3SU                => micstt%CNO3SU              , &
   CNO3BU                => micstt%CNO3BU              , &
   OMEheter              => micstt%OMEheter            , &
   RINOB                 => micflx%RINOB               , &
   RINOO                 => micflx%RINOO               , &
   RINHO                 => micflx%RINHO               , &
   RINHB                 => micflx%RINHB               , &
   RIPOO                 => micflx%RIPOO               , &
   RIPBO                 => micflx%RIPBO               , &
   RIPO1                 => micflx%RIPO1               , &
   RIPB1                 => micflx%RIPB1               , &
   RINHOR                => micflx%RINHOR              , &
   RINOOR                => micflx%RINOOR              , &
   RIPOOR                => micflx%RIPOOR              , &
   RIPO1R                => micflx%RIPO1R              , &
   NetNH4Mineralize_col  => micflx%NetNH4Mineralize_col, &
   NetPO4Mineralize_col  => micflx%NetPO4Mineralize_col  &
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
!     VOLW=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     NetNH4Mineralize_col=total NH4 net mineraln (-ve) or immobiln (+ve)
! update may be needed, May 17th, 2023, jyt.
  FNH4S=VLNH4
  FNHBS=VLNHB
  MID3=micpar%get_micb_id(3,NGL)
  RINHP=(OMEheter(ielmc,MID3,K)*rNCOMC(3,NGL,K)-OMEheter(ielmn,MID3,K))
  IF(RINHP.GT.0.0_r8)THEN
    CNH4X=AZMAX1(CNH4S-Z4MN)
    CNH4Y=AZMAX1(CNH4B-Z4MN)
    RINHX=AMIN1(RINHP,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*Z4MX)
    RINHO(NGL,K)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RINHB(NGL,K)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M=Z4MN*VOLWU*FNH4S
    ZNHBM=Z4MN*VOLWU*FNHBS
    RNH4TransfSoilHeter(NGL,K)=AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RINHO(NGL,K))
    RNH4TransfBandHeter(NGL,K)=AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RINHB(NGL,K))
  ELSE
    RINHO(NGL,K)=0.0_r8
    RINHB(NGL,K)=0.0_r8
    RNH4TransfSoilHeter(NGL,K)=RINHP*FNH4S
    RNH4TransfBandHeter(NGL,K)=RINHP*FNHBS
  ENDIF
  NetNH4Mineralize_col=NetNH4Mineralize_col+(RNH4TransfSoilHeter(NGL,K)+RNH4TransfBandHeter(NGL,K))
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
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     NetNH4Mineralize_col=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
  FNO3S=VLNO3
  FNO3B=VLNOB
  RINOP=AZMAX1(RINHP-RNH4TransfSoilHeter(NGL,K)-RNH4TransfBandHeter(NGL,K))
  IF(RINOP.GT.0.0_r8)THEN
    CNO3X=AZMAX1(CNO3S-ZOMN)
    CNO3Y=AZMAX1(CNO3B-ZOMN)
    RINOX=AMIN1(RINOP,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*ZOMX)
    RINOO(NGL,K)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RINOB(NGL,K)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M=ZOMN*VOLWU*FNO3S
    ZNOBM=ZOMN*VOLWU*FNO3B
    RNO3TransfSoilHeter(NGL,K)=AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)),RINOO(NGL,K))
    RNO3TransfBandHeter(NGL,K)=AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)),RINOB(NGL,K))
  ELSE
    RINOO(NGL,K)=0.0_r8
    RINOB(NGL,K)=0.0_r8
    RNO3TransfSoilHeter(NGL,K)=RINOP*FNO3S
    RNO3TransfBandHeter(NGL,K)=RINOP*FNO3B
  ENDIF
  NetNH4Mineralize_col=NetNH4Mineralize_col+(RNO3TransfSoilHeter(NGL,K)+RNO3TransfBandHeter(NGL,K))
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
  RIPOP=(OMEheter(ielmc,MID3,K)*rPCOMC(3,NGL,K)-OMEheter(ielmp,MID3,K))
  IF(RIPOP.GT.0.0_r8)THEN
    CH2PX=AZMAX1(CH2P4-HPMN)
    CH2PY=AZMAX1(CH2P4B-HPMN)
    RIPOX=AMIN1(RIPOP,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*HPMX)
    RIPOO(NGL,K)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RIPBO(NGL,K)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM=HPMN*VLWatMicP*FH2PS
    H2PBM=HPMN*VLWatMicP*FH2PB
    RH2PO4TransfSoilHeter(NGL,K)=AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RIPOO(NGL,K))
    RH2PO4TransfBandHeter(NGL,K)=AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RIPBO(NGL,K))
  ELSE
    RIPOO(NGL,K)=0.0_r8
    RIPBO(NGL,K)=0.0_r8
    RH2PO4TransfSoilHeter(NGL,K)=RIPOP*FH2PS
    RH2PO4TransfBandHeter(NGL,K)=RIPOP*FH2PB
  ENDIF
  NetPO4Mineralize_col=NetPO4Mineralize_col+(RH2PO4TransfSoilHeter(NGL,K)+RH2PO4TransfBandHeter(NGL,K))
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
  RIP1P=0.1_r8*AZMAX1(RIPOP-RH2PO4TransfSoilHeter(NGL,K)-RH2PO4TransfBandHeter(NGL,K))
  IF(RIP1P.GT.0.0_r8)THEN
    CH1PX=AZMAX1(CH1P4-HPMN)
    CH1PY=AZMAX1(CH1P4B-HPMN)
    RIP1X=AMIN1(RIP1P,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*HPMX)
    RIPO1(NGL,K)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RIPB1(NGL,K)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM=HPMN*VLWatMicP*FH1PS
    H1PBM=HPMN*VLWatMicP*FH1PB
    RH1PO4TransfSoilHeter(NGL,K)=AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RIPO1(NGL,K))
    RH1PO4TransfBandHeter(NGL,K)=AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RIPB1(NGL,K))
  ELSE
    RIPO1(NGL,K)=0.0_r8
    RIPB1(NGL,K)=0.0_r8
    RH1PO4TransfSoilHeter(NGL,K)=RIP1P*FH1PS
    RH1PO4TransfBandHeter(NGL,K)=RIP1P*FH1PB
  ENDIF
  NetPO4Mineralize_col=NetPO4Mineralize_col+(RH1PO4TransfSoilHeter(NGL,K)+RH1PO4TransfBandHeter(NGL,K))
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
    RINHPR=RINHP-RNH4TransfSoilHeter(NGL,K)-RNO3TransfSoilHeter(NGL,K)
    IF(RINHPR.GT.0.0_r8)THEN
      CNH4X=AZMAX1(CNH4SU-Z4MN)
      CNH4Y=AZMAX1(CNH4BU-Z4MN)
      RINHOR(NGL,K)=AMIN1(RINHPR,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VLWatMicP
      RNH4TransfLitrHeter(NGL,K)=AMIN1(FNH4XR(NGL,K)*AZMAX1((ZNH4TU-ZNH4M)),RINHOR(NGL,K))
    ELSE
      RINHOR(NGL,K)=0.0_r8
      RNH4TransfLitrHeter(NGL,K)=RINHPR
    ENDIF
    NetNH4Mineralize_col=NetNH4Mineralize_col+RNH4TransfLitrHeter(NGL,K)
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
!     VOLW=water content
!     ZNO3M=NO3 not available for uptake
!     FNO3XR=fraction of biological NO3 demand
!     RNO3TransfLitrHeter=substrate-limited NO3 immobiln
!     NetNH4Mineralize_col=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
    RINOPR=AZMAX1(RINHPR-RNH4TransfLitrHeter(NGL,K))
    IF(RINOPR.GT.0.0_r8)THEN
      CNO3X=AZMAX1(CNO3SU-ZOMN)
      CNO3Y=AZMAX1(CNO3BU-ZOMN)
      RINOOR(NGL,K)=AMAX1(RINOPR,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VLWatMicP
      RNO3TransfLitrHeter(NGL,K)=AMIN1(FNO3XR(NGL,K)*AZMAX1((ZNO3TU-ZNO3M)),RINOOR(NGL,K))
    ELSE
      RINOOR(NGL,K)=0.0_r8
      RNO3TransfLitrHeter(NGL,K)=RINOPR
    ENDIF
    NetNH4Mineralize_col=NetNH4Mineralize_col+RNO3TransfLitrHeter(NGL,K)
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
!     RH2PO4TransfLitrHeter=substrate-limited H2PO4 mineraln-immobiln
!     NetPO4Mineralize_col=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    RIPOPR=RIPOP-RH2PO4TransfSoilHeter(NGL,K)
    IF(RIPOPR.GT.0.0_r8)THEN
      CH2PX=AZMAX1(CH2P4U-HPMN)
      CH2PY=AZMAX1(CH2P4BU-HPMN)
      RIPOOR(NGL,K)=AMIN1(RIPOPR,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M=HPMN*VOLWU
      RH2PO4TransfLitrHeter(NGL,K)=AMIN1(FPO4XR(NGL,K)*AZMAX1((H2P4TU-H2P4M)),RIPOOR(NGL,K))
    ELSE
      RIPOOR(NGL,K)=0.0_r8
      RH2PO4TransfLitrHeter(NGL,K)=RIPOPR
    ENDIF
    NetPO4Mineralize_col=NetPO4Mineralize_col+RH2PO4TransfLitrHeter(NGL,K)
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
!     RH1PO4TransfLitrHeter=substrate-limited HPO4 minereraln-immobiln
!     NetPO4Mineralize_col=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RIP1PR=0.1_r8*AZMAX1(RIPOPR-RH2PO4TransfLitrHeter(NGL,K))
    IF(RIP1PR.GT.0.0_r8)THEN
      CH1PX=AZMAX1(CH1P4U-HPMN)
      CH1PY=AZMAX1(CH1P4BU-HPMN)
      RIPO1R(NGL,K)=AMIN1(RIP1PR,BIOA*OMActHeter(NGL,K)*TFNG(NGL,K)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M=HPMN*VOLWU
      RH1PO4TransfLitrHeter(NGL,K)=AMIN1(FP14XR(NGL,K)*AZMAX1((H1P4TU-H1P4M)),RIPO1R(NGL,K))
    ELSE
      RIPO1R(NGL,K)=0.0_r8
      RH1PO4TransfLitrHeter(NGL,K)=RIP1PR
    ENDIF
    NetPO4Mineralize_col=NetPO4Mineralize_col+RH1PO4TransfLitrHeter(NGL,K)
  ENDIF
  end associate
  end subroutine BiomassMineralization
!------------------------------------------------------------------------------------------

  subroutine GatherMicrobialRespiration(NGL,N,K,RMOMK,RGrowthRespHeter,RMaintDefcitAutor,RMaintRespHeter, &
    micfor,micstt,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(out) :: RGrowthRespHeter,RMaintDefcitAutor
  real(r8), intent(out) :: RMaintRespHeter
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  integer :: MID3,MID1
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
!     begin_execution
  associate(                                          &
    TFNR                => nmics%TFNR               , &
    OMN2                => nmics%OMN2               , &
    Resp4NFixHeter      => nmicf%Resp4NFixHeter     , &
    RespGrossHeter      => nmicf%RespGrossHeter     , &
    RMOMC               => nmicf%RMOMC              , &
    RNH4TransfSoilHeter => nmicf%RNH4TransfSoilHeter, &
    RNO3TransfSoilHeter => nmicf%RNO3TransfSoilHeter, &
    RN2FixHeter         => nmicf%RN2FixHeter        , &
    rNCOMC              => micpar%rNCOMC            , &
    rPCOMC              => micpar%rPCOMC            , &
    pH                  => micfor%pH                , &
    ZEROS               => micfor%ZEROS             , &
    n_aero_n2fixer      => micpar%n_aero_n2fixer    , &
    n_anero_n2fixer     => micpar%n_anero_n2fixer   , &
    CZ2GS               => micstt%CZ2GS             , &
    OMEheter            => micstt%OMEheter            &
  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TFNR=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  FPH=1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))

  RMOMX=RMOM*TFNR(NGL,K)*FPH
  MID1=micpar%get_micb_id(1,NGL)
  RMOMC(1,NGL,K)=OMEheter(ielmn,MID1,K)*RMOMX*RMOMK(1)
  RMOMC(2,NGL,K)=OMN2(NGL,K)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMaintRespHeter=total maintenance respiration
!     RGrowthRespHeter=growth respiration
!     RMaintDefcitAutor=senescence respiration
!
  RMaintRespHeter=RMOMC(1,NGL,K)+RMOMC(2,NGL,K)
  RGrowthRespHeter=AZMAX1(RespGrossHeter(NGL,K)-RMaintRespHeter)
  RMaintDefcitAutor=AZMAX1(RMaintRespHeter-RespGrossHeter(NGL,K))
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
!     RGrowthRespHeter=growth respiration
!     Resp4NFixHeter=respiration for N2 fixation
!     CZ2GS=aqueous N2 concentration
!     ZFKM=Km for N2 uptake
!     OMGR*OMC(3,NGL,N,K)=nonstructural C limitation to Resp4NFixHeter
!     RN2FixHeter=N2 fixation rate
!
  IF(N.EQ.n_aero_n2fixer.OR.N.EQ.n_anero_n2fixer)THEN
    MID3=micpar%get_micb_id(3,NGL)
    RGN2P=AZMAX1(OMEheter(ielmc,MID3,K)*rNCOMC(3,NGL,K)-OMEheter(ielmn,MID3,K))/EN2F(N)
    IF(RGrowthRespHeter.GT.ZEROS)THEN
      Resp4NFixHeter(NGL,K)=AMIN1(RGrowthRespHeter*RGN2P/(RGrowthRespHeter+RGN2P) &
        *CZ2GS/(CZ2GS+ZFKM),OMGR*OMEheter(ielmc,MID3,K))
    ELSE
      Resp4NFixHeter(NGL,K)=0.0_r8
    ENDIF
    RN2FixHeter(NGL,K)=Resp4NFixHeter(NGL,K)*EN2F(N)
  ENDIF
  end associate
  end subroutine GatherMicrobialRespiration
!------------------------------------------------------------------------------------------

  subroutine GetMicrobialAnabolismFlux(NGL,N,K,ECHZ,FGOCP,FGOAP, &
    RGrowthRespHeter,RMaintDefcitAutor,RMaintRespHeter,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: ECHZ
  real(r8), intent(in) :: FGOCP,FGOAP,RGrowthRespHeter,RMaintDefcitAutor,RMaintRespHeter
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  integer :: M,MID1,MID3,MID,NE
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: CXC
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FRM
!     begin_execution
  associate(                                          &
    CNOMActHeter        => nmics%CNOMActHeter,        &
    CPOMActHeter        => nmics%CPOMActHeter,        &
    TFNG                => nmics%TFNG,                &
    FCN                 => nmics%FCN,                 &
    FCP                 => nmics%FCP,                 &
    FracHeterBiomOfActK => nmics%FracHeterBiomOfActK, &
    CGOMEheter          => nmicf%CGOMEheter,          &
    CGOQC               => nmicf%CGOQC,               &
    CGOAC               => nmicf%CGOAC,               &
    CGOMES              => nmicf%CGOMES,              &
    RespGrossHeter      => nmicf%RespGrossHeter,      &
    RGOMD               => nmicf%RGOMD,               &
    RMOMC               => nmicf%RMOMC,               &
    RDOMEheter          => nmicf%RDOMEheter,          &
    RHOMEheter          => nmicf%RHOMEheter,          &
    RCOMEheter          => nmicf%RCOMEheter,          &
    RHMMEheter          => nmicf%RHMMEheter,          &
    RCMMEheter          => nmicf%RCMMEheter,          &
    RDMMEheter          => nmicf%RDMMEheter,          &
    RXOMEheter          => nmicf%RXOMEheter,          &
    R3OMEheter          => nmicf%R3OMEheter,          &
    RXMMEheter          => nmicf%RXMMEheter,          &
    R3MMEheter          => nmicf%R3MMEheter,          &
    Resp4NFixHeter      => nmicf%Resp4NFixHeter,      &
    TCGOMEheter         => ncplxf%TCGOMEheter,        &
    CNQ                 => ncplxs%CNQ,                &
    CPQ                 => ncplxs%CPQ,                &
    rNCOMC              => micpar%rNCOMC,             &
    rPCOMC              => micpar%rPCOMC,             &
    FL                  => micpar%FL,                 &
    EHUM                => micstt%EHUM,               &
    OMEheter            => micstt%OMEheter,           &
    DOM                 => micstt%DOM,                &
    ZEROS               => micfor%ZEROS,              &
    ZERO                => micfor%ZERO                &
  )

!     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
!     FROM O2, NOX AND C REDUCTION
!
!     CGOMX=DOC+acetate uptake from aerobic growth respiration
!     CGOMD=DOC+acetate uptake from denitrifier growth respiration
!     RMaintRespHeter=maintenance respiration
!     RespGrossHeter=total respiration
!     RGOMD=respiration for denitrifcation
!     Resp4NFixHeter=respiration for N2 fixation
!     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
!     CGOMEheter,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake(heterotrophs
!     CGOMEheter=total CO2,CH4 uptake (autotrophs)
!     CGOMEheter,CGOMEheter=DON, DOP uptake
!     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
!     OQN,OPQ=DON,DOP
!     FracHeterBiomOfActK=faction of OMActHeterin total OMA
!     CNQ,CPQ=DON/DOC, DOP/DOC
!     FCN,FCP=limitation from N,P
!

  CGOMX=AMIN1(RMaintRespHeter,RespGrossHeter(NGL,K))+Resp4NFixHeter(NGL,K)+(RGrowthRespHeter-Resp4NFixHeter(NGL,K))/ECHZ
  CGOMD=RGOMD(NGL,K)/ENOX
  CGOMEheter(ielmc,NGL,K)=CGOMX+CGOMD

  CGOQC(NGL,K)=CGOMX*FGOCP+CGOMD
  CGOAC(NGL,K)=CGOMX*FGOAP
  CGOXC=CGOQC(NGL,K)+CGOAC(NGL,K)
  CGOMEheter(ielmp,NGL,K)=AZMAX1(AMIN1(DOM(idom_don,K)*FracHeterBiomOfActK(NGL,K),CGOXC*CNQ(K)/FCN(NGL,K)))
  CGOMEheter(ielmp,NGL,K)=AZMAX1(AMIN1(DOM(idom_dop,K)*FracHeterBiomOfActK(NGL,K),CGOXC*CPQ(K)/FCP(NGL,K)))
  TCGOMEheter(ielmc,K)=TCGOMEheter(ielmc,K)+CGOQC(NGL,K)
  TCGOMEheter(idom_acetate,K)=TCGOMEheter(idom_acetate,K)+CGOAC(NGL,K)
  TCGOMEheter(ielmp,K)=TCGOMEheter(ielmp,K)+CGOMEheter(ielmp,NGL,K)
  TCGOMEheter(ielmp,K)=TCGOMEheter(ielmp,K)+CGOMEheter(ielmp,NGL,K)
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
  IF(OMEheter(ielmc,MID3,K).GT.ZEROS .AND.OMEheter(ielmc,MID1,K).GT.ZEROS)THEN
    CCC=AZMAX1(AMIN1(1.0_r8 &
      ,OMEheter(ielmn,MID3,K)/(OMEheter(ielmn,MID3,K)+OMEheter(ielmc,MID3,K)*rNCOMC(3,NGL,K)) &
      ,OMEheter(ielmp,MID3,K)/(OMEheter(ielmp,MID3,K)+OMEheter(ielmc,MID3,K)*rPCOMC(3,NGL,K))))
    CXC=OMEheter(ielmc,MID3,K)/OMEheter(ielmc,MID1,K)
    C3C=1.0_r8/(1.0_r8+CXC/CKC)
    CNC=AZMAX1(AMIN1(1.0_r8,OMEheter(ielmc,MID3,K)/(OMEheter(ielmc,MID3,K)+OMEheter(ielmn,MID3,K)/rNCOMC(3,NGL,K))))
    CPC=AZMAX1(AMIN1(1.0_r8,OMEheter(ielmc,MID3,K)/(OMEheter(ielmc,MID3,K)+OMEheter(ielmp,MID3,K)/rPCOMC(3,NGL,K))))
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
  CGOMZ=TFNG(NGL,K)*OMGR*AZMAX1(OMEheter(ielmc,MID3,K))
  D745: DO M=1,2
    CGOMES(ielmc,M,NGL,K)=FL(M)*CGOMZ
    IF(OMEheter(ielmc,MID3,K).GT.ZEROS)THEN
      CGOMES(ielmn,M,NGL,K)=AMIN1(FL(M)*AZMAX1(OMEheter(ielmn,MID3,K)) &
        ,CGOMES(ielmc,M,NGL,K)*OMEheter(ielmn,MID3,K)/OMEheter(ielmc,MID3,K))
      CGOMES(ielmp,M,NGL,K)=AMIN1(FL(M)*AZMAX1(OMEheter(ielmp,MID3,K)) &
        ,CGOMES(ielmc,M,NGL,K)*OMEheter(ielmp,MID3,K)/OMEheter(ielmc,MID3,K))
    ELSE
      CGOMES(ielmn,M,NGL,K)=0.0_r8
      CGOMES(ielmp,M,NGL,K)=0.0_r8
    ENDIF
!
!     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
!     RATE, TEMPERATURE
!
!     SPOMX=rate constant for microbial decomposition
!     SPOMC=basal decomposition rate
!     SPOMK=effect of low microbial C concentration on microbial decay
!     RXOMC,RXOMN,RXOMEheter=microbial C,N,P decomposition
!     RDOMEheter,RDOMN,RDOMP=microbial C,N,P LitrFall
!     R3OMC,R3OMN,R3OMEheter=microbial C,N,P recycling
!
    MID=micpar%get_micb_id(M,NGL)
    SPOMX=SQRT(TFNG(NGL,K))*SPOMC(M)*SPOMK(M)
    RXOMEheter(ielmc,M,NGL,K)=AZMAX1(OMEheter(ielmc,MID,K)*SPOMX)
    RXOMEheter(ielmn,M,NGL,K)=AZMAX1(OMEheter(ielmn,MID,K)*SPOMX)
    RXOMEheter(ielmp,M,NGL,K)=AZMAX1(OMEheter(ielmp,MID,K)*SPOMX)
    RDOMEheter(ielmc,M,NGL,K)=RXOMEheter(ielmc,M,NGL,K)*(1.0_r8-RCCC)
    RDOMEheter(ielmn,M,NGL,K)=RXOMEheter(ielmn,M,NGL,K)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDOMEheter(ielmp,M,NGL,K)=RXOMEheter(ielmp,M,NGL,K)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
    DO NE=1,NumPlantChemElms    
      R3OMEheter(NE,M,NGL,K)=RXOMEheter(NE,M,NGL,K)-RDOMEheter(NE,M,NGL,K)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMEheter=transfer of microbial C,N,P LitrFall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMEheter,RCOMEheter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!

      RHOMEheter(NE,M,NGL,K)=AZMAX1(RDOMEheter(NE,M,NGL,K)*EHUM)
  !
  !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
  !
      RCOMEheter(NE,M,NGL,K)=RDOMEheter(NE,M,NGL,K)-RHOMEheter(NE,M,NGL,K)
    ENDDO
  ENDDO D745
!
!     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
!     EXCEEDS UPTAKE
!
!     OMC,OMN,OMP=microbial C,N,P
!     RMaintRespHeter=total maintenance respiration
!     RMaintDefcitAutor=senescence respiration
!     RCCC=C recycling fraction
!     RXMMEheter,RXMMN,RXMMP=microbial C,N,P loss from senescence
!     RMOMC=maintenance respiration
!     CNOMA,CPOMA=N:C,P:C ratios of active biomass
!     RDMMC,RDMMEheter,RDMMP=microbial C,N,P LitrFall from senescence
!     R3MMEheter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!
  IF(RMaintDefcitAutor.GT.ZEROS.AND.RMaintRespHeter.GT.ZEROS.AND.RCCC.GT.ZERO)THEN
    FRM=RMaintDefcitAutor/RMaintRespHeter
    D730: DO M=1,2
      MID=micpar%get_micb_id(M,NGL)    
      RXMMEheter(ielmc,M,NGL,K)=AMIN1(OMEheter(ielmc,MID,K),AZMAX1(FRM*RMOMC(M,NGL,K)/RCCC))
      RXMMEheter(ielmn,M,NGL,K)=AMIN1(OMEheter(ielmn,MID,K),AZMAX1(RXMMEheter(ielmc,M,NGL,K)*CNOMActHeter(NGL,K)))
      RXMMEheter(ielmp,M,NGL,K)=AMIN1(OMEheter(ielmp,MID,K),AZMAX1(RXMMEheter(ielmc,M,NGL,K)*CPOMActHeter(NGL,K)))
      RDMMEheter(ielmc,M,NGL,K)=RXMMEheter(ielmc,M,NGL,K)*(1.0_r8-RCCC)
      RDMMEheter(ielmn,M,NGL,K)=RXMMEheter(ielmn,M,NGL,K)*(1.0_r8-RCCN)*(1.0_r8-RCCC)
      RDMMEheter(ielmp,M,NGL,K)=RXMMEheter(ielmp,M,NGL,K)*(1.0_r8-RCCP)*(1.0_r8-RCCC)
      DO NE=1,NumPlantChemElms
        R3MMEheter(NE,M,NGL,K)=RXMMEheter(NE,M,NGL,K)-RDMMEheter(NE,M,NGL,K)
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RHMMEheter,RHMMN,RHMMEheter=transfer of senesence LitrFall C,N,P to humus
!     EHUM=humus transfer fraction
!     RCMMEheter,RCMMN,RCMMEheter=transfer of senesence LitrFall C,N,P to residue
!

        RHMMEheter(NE,M,NGL,K)=AZMAX1(RDMMEheter(NE,M,NGL,K)*EHUM)
        RCMMEheter(NE,M,NGL,K)=RDMMEheter(NE,M,NGL,K)-RHMMEheter(NE,M,NGL,K)
      ENDDO
    ENDDO D730
  ELSE
    D720: DO M=1,2
      DO NE=1,NumPlantChemElms
        RXMMEheter(NE,M,NGL,K)=0.0_r8
        RDMMEheter(NE,M,NGL,K)=0.0_r8
        R3MMEheter(NE,M,NGL,K)=0.0_r8
        RHMMEheter(NE,M,NGL,K)=0.0_r8
        RCMMEheter(NE,M,NGL,K)=0.0_r8
      ENDDO

    ENDDO D720
  ENDIF
  end associate
  end subroutine GetMicrobialAnabolismFlux

end module MicBGCMod
