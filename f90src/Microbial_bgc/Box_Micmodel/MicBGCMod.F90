module MicBGCMod
!!
! DESCRIPTION:
! codes to do soil biological transfOMBioResduations
!
! USES:
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use abortutils,           only: endrun,   destroy
  use EcoSIMCtrlMod,        only: etimer
  use minimathmod,          only: safe_adb, AZMAX1
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: DebugPrint
  use TracerIDMod
  use MicAutoCPLXMod
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod
  use MicrobMathFuncMod
  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: jcplx,NumMicbFunGrupsPerCmplx,jsken,ndbiomcp,nlbiomcp
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
  NumMicbFunGrupsPerCmplx  =micpar%NumMicbFunGrupsPerCmplx
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

  subroutine SoilBGCOneLayer(I,J,micfor,micstt,micflx,naqfdiag,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(out) :: naqfdiag
  type(Microbe_Diag_type),intent(out) :: nmicdiag

!local variables
  integer :: LL,K,KL,NGL
  integer :: M,N

  type(Microbe_State_type) :: nmics
  type(Microbe_Flux_type) :: nmicf
  type(OMCplx_Flux_type) :: ncplxf
  type(OMCplx_State_type) :: ncplxs
  real(r8) :: totOMbeg(1:NumPlantChemElms)
  real(r8) :: totOMend(1:NumPlantChemElms)

! begin_execution
  call nmicf%Init(jcplx,NumMicbFunGrupsPerCmplx)
  call nmics%Init(jcplx,NumMicbFunGrupsPerCmplx)
  call ncplxf%Init()
  call ncplxs%Init()
  call naqfdiag%ZeroOut()

  micflx%tRHydlySOM = 0._r8;   micflx%tRHydlyBioReSOM = 0._r8; micflx%tRHydlySoprtOM = 0._r8
  micflx%NetCAssimhr = 0._r8;micflx%NetNH4Mineralize = 0._r8;micflx%NetPO4Mineralize = 0._r8
  micflx%RPiDemand   = 0._r8; micflx%RNiDemand       = 0._r8;micflx%GrosAssimhr      = 0._r8
  micflx%CDOMuptk1   = 0._r8;micflx%CDOMuptk2        = 0._r8;micflx%tROMT            = 0._r8  
  micflx%tGROMO     = 0._r8; micflx%tRGOMP            = 0._r8;micflx%tRGOXP          = 0._r8
  micflx%tRGOZP     = 0._r8

! write(*,*)'StageBGCEnvironCondition'
  call StageBGCEnvironCondition(I,J,micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)

!  call SumOneLayer(KL,micstt,totOMbeg)
!
! write(*,*)'ActiveMicrobes'
  call ActiveMicrobes(I,J,micfor,micstt,micflx,nmicdiag, &
    naqfdiag,nmicf,nmics,ncplxf,ncplxs)
!
!write(*,*)'ChemoDenitrification'
  call ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)
!
!     DECOMPOSITION
!
!     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
!
  DO  K=1,KL
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        ncplxf%ROQC4HeterMicActCmpK(K)=ncplxf%ROQC4HeterMicActCmpK(K)+nmicf%ROQC4HeterMicrobAct(NGL,K)
      enddo
    ENDDO
  ENDDO
        !
        !write(*,*)'PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C'
  call OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
        !
        !     TRANSFER ALL PRIMING AMONG ALL K
        !
        !     TMicHeterActivity=total respiration of DOC+DOA in soil layer
        !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA in substrate complex
        !     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
        !     OMC,OMN,OMP=microbial C,N,P
        !
!
    !write(*,*)'DECOMPOSITION OF ORGANIC SUBSTRATES'
!
    call SolidOMDecomposition(I,J,KL,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs,micflx)
!
          !write(*,*)'DOC ADSORPTION - DESORPTION'
!
    call RDOMSorption(KL,micfor,micstt,nmicf,ncplxf,ncplxs)


        !write(*,*)'RedistDecompProduct'
  call RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)
!
        !write(*,*)'MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION'

  call AutotrophAnabolicUpdate(micfor,micstt,nmicf)

  call HeterotrophAnabolicUpdate(I,J,micfor,micstt,nmicf,micflx)
!
        !write(*,*)'MICROBIAL COLONIZATION OF NEW LITTER'
!
  call MicrobialLitterColonization(I,J,KL,micfor,micstt,ncplxf,ncplxs)
!
!     AGGREGATE ALL TRANSFOMBioResduATIONS CALCULATED ABOVE FOR EACH N,K
!
  call AggregateTransfOMBioResdue(micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)

!  call SumOneLayer(KL,micstt,totOMend)

  micstt%TotActMicrobiom=nmicdiag%TotActMicrobiom
  call nmics%destroy()
  call nmicf%destroy()
  call ncplxf%destroy()
  call ncplxs%destroy()

  end subroutine SoilBGCOneLayer
!------------------------------------------------------------------------------------------
  subroutine SumOneLayer(KL,micstt,toms)
  implicit none
  integer, intent(in) :: KL
  type(micsttype), intent(in) :: micstt  
  real(r8), intent(out) :: toms(1:NumPlantChemElms)
  integer :: K,M,NBM,NE,N

!     begin_execution
  associate(                                      &  
    SolidOM           => micstt%SolidOM,          &
    OMBioResdu        => micstt%OMBioResdu,       &
    SorbedOM          => micstt%SorbedOM,         &
    DOM               => micstt%DOM,              &
    mBiomeHeter       => micstt%mBiomeHeter,      &
    mBiomeAutor       => micstt%mBiomeAutor,      &
    NumLiveAutoBioms  => micpar%NumLiveAutoBioms, &
    NumLiveHeterBioms => micpar%NumLiveHeterBioms &
  )

  toms=0._r8
  DO K=1,KL
    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        TOMS(NE)=TOMS(NE)+SolidOM(NE,M,K)
      ENDDO
    ENDDO

    DO NBM=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        TOMS(NE)=TOMS(NE)+OMBioResdu(NE,NBM,K)
      ENDDO
    ENDDO

    DO NE=1,NumPlantChemElms
      TOMS(NE)=TOMS(NE)+SorbedOM(NE,K)
      TOMS(NE)=TOMS(NE)+DOM(NE,K)
    ENDDO

    DO N=1,NumLiveHeterBioms
      DO NE=1,NumPlantChemElms
        TOMS(NE)=TOMS(NE)+mBiomeHeter(NE,N,K)
      ENDDO
    ENDDO
  ENDDO

  DO N=1,NumLiveAutoBioms
    DO NE=1,NumPlantChemElms
      TOMS(NE)=TOMS(NE)+mBiomeAutor(NE,N)      
    ENDDO
  ENDDO
  
  end associate

  end subroutine SumOneLayer
!------------------------------------------------------------------------------------------

  subroutine StageBGCEnvironCondition(I,J,micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  integer, intent(out) :: KL
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_State_type),intent(inout):: ncplxs
  real(r8) :: TBulkOMC  
  integer  :: K
  integer  :: M,N,NGL,MID1,MID2
  real(r8) :: ORGCL
  real(r8) :: TKSO
  real(r8) :: TSolidOMC,TSorbedOMC

!     begin_execution
  associate(                                                       &
    rCNBiomeActHeter          => nmics%rCNBiomeActHeter,           &
    OMActHeter                => nmics%OMActHeter,                 &
    OMC2                      => nmics%OMC2,                       &
    OMN2                      => nmics%OMN2,                       &
    FOM2                      => nmics%FOM2,                       &
    FCN                       => nmics%FCN,                        &
    FCP                       => nmics%FCP,                        &
    FBiomStoiScalarHeter      => nmics%FBiomStoiScalarHeter,       &
    BulkSOMC                  => ncplxs%BulkSOMC,                  &
    TOMEAutoK                 => ncplxs%TOMEAutoK,                 &
    TOMK                      => ncplxs%TOMK,                      &
    TONK                      => ncplxs%TONK,                      &
    TOPK                      => ncplxs%TOPK,                      &
    FOCA                      => ncplxs%FOCA,                      &
    FOAA                      => ncplxs%FOAA,                      &
    rCNDOM                    => ncplxs%rCNDOM,                    &
    rCPDOM                    => ncplxs%rCPDOM,                    &
    CDOM                      => ncplxs%CDOM,                      &
    OMBioResduK               => ncplxs%OMBioResduK,               &
    SolidOMCK                  => ncplxs%SolidOMCK,                  &
    SolidOMActK               => ncplxs%SolidOMActK,               &
    tMaxNActMicrbK            => ncplxs%tMaxNActMicrbK,            &
    tMaxPActMicrbK            => ncplxs%tMaxPActMicrbK,            &
    TOMBioResdu               => nmicdiag%TOMBioResdu,             &
    TotActMicrobiom           => nmicdiag%TotActMicrobiom,         &
    TotBiomNO2Consumers       => nmicdiag%TotBiomNO2Consumers,     &
    XCO2                      => nmicdiag%XCO2,                    &
    TSensGrowth               => nmicdiag%TSensGrowth,             &
    TSensMaintR               => nmicdiag%TSensMaintR,             &
    ThetaLitr                 => nmicdiag%ThetaLitr,               &
    ThetaZ                    => nmicdiag%ThetaZ,                  &
    VOLWZ                     => nmicdiag%VOLWZ,                   &
    ZNH4T                     => nmicdiag%ZNH4T,                   &
    ZNO3T                     => nmicdiag%ZNO3T,                   &
    ZNO2T                     => nmicdiag%ZNO2T,                   &
    H2P4T                     => nmicdiag%H2P4T,                   &
    H1P4T                     => nmicdiag%H1P4T,                   &
    rCNBiomeActAutor          => nmics%rCNBiomeActAutor,           &
    OMActAutor                => nmics%OMActAutor,                 &
    OMC2Autor                 => nmics%OMC2Autor,                  &
    OMN2Autor                 => nmics%OMN2Autor,                  &
    FOM2Autor                 => nmics%FOM2Autor,                  &
    FCNAutor                  => nmics%FCNAutor,                   &
    FCPAutor                  => nmics%FCPAutor,                   &
    FBiomStoiScalarAutor      => nmics%FBiomStoiScalarAutor,       &
    rNCOMCAutor               => micpar%rNCOMCAutor,               &
    rPCOMCAutor               => micpar%rPCOMCAutor,               &
    rNCOMC                    => micpar%rNCOMC,                    &
    rPCOMC                    => micpar%rPCOMC,                    &
    FL                        => micpar%FL,                        &
    k_humus                   => micpar%k_humus,                   &
    k_POM                     => micpar%k_POM,                     &
    is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, &
    mid_Facult_DenitBacter    => micpar%mid_Facult_DenitBacter,    &
    mid_AmmoniaOxidBacter     => micpar%mid_AmmoniaOxidBacter,     &
    litrm                     => micfor%litrm,                     &
    VWatLitRHoldCapcity       => micfor%VWatLitRHoldCapcity,       &
    VLWatMicP                 => micfor%VLWatMicP,                 &
    VOLW0                     => micfor%VOLW0,                     &
    THETY                     => micfor%THETY,                     &
    VLitR                     => micfor%VLitR,                     &
    VLSoilMicP                => micfor%VLSoilMicP,                &
    POROS                     => micfor%POROS,                     &
    ZEROS                     => micfor%ZEROS,                     &
    FieldCapacity             => micfor%FieldCapacity,             &
    THETW                     => micfor%THETW,                     &
    TKS                       => micfor%TKS,                       &
    TempOffset                => micfor%TempOffset,                &
    VLWatMicPM                => micfor%VLWatMicPM,                &
    ZEROS2                    => micfor%ZEROS2,                    &
    ZERO                      => micfor%ZERO,                      &
    CCO2S                     => micstt%CCO2S,                     &
    tOMActC                   => micstt%tOMActC                  ,  &
    SolidOM                   => micstt%SolidOM,                   &
    SolidOMAct                => micstt%SolidOMAct,                &
    TSolidOMActC              => micstt%TSolidOMActC  ,            &
    TSolidOMC                 => micstt%TSolidOMC  ,               &
    OMBioResdu                => micstt%OMBioResdu,                &
    SorbedOM                  => micstt%SorbedOM,                  &
    mBiomeHeter               => micstt%mBiomeHeter,               &
    DOM                       => micstt%DOM,                       &
    H1PO4                     => micstt%H1PO4,                     &
    H1POB                     => micstt%H1POB,                     &
    H2PO4                     => micstt%H2PO4,                     &
    H2POB                     => micstt%H2POB,                     &
    ZNH4B                     => micstt%ZNH4B,                     &
    ZNH4S                     => micstt%ZNH4S,                     &
    ZNO2B                     => micstt%ZNO2B,                     &
    ZNO2S                     => micstt%ZNO2S,                     &
    ZNO3B                     => micstt%ZNO3B,                     &
    ZNO3S                     => micstt%ZNO3S,                     &
    mBiomeAutor               => micstt%mBiomeAutor,               &
    FracBulkSOMC              => micstt%FracBulkSOMC               &
  )

! get KL, the number of mic-om complexes

!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH TempOffset FOR THERMAL ADAPTATION
  IF(litrm)THEN
    ! surface litter layer
    KL=micpar%NumOfLitrCmplxs
    IF(VWatLitRHoldCapcity.GT.ZEROS2)THEN
      ThetaLitr = AMIN1(VOLW0/VLitR,1._r8)
      ThetaZ    = AZMAX1(ThetaLitr-THETY)
      VOLWZ     = ThetaZ*VLitR    !effective water volume
!      if(etimer%get_curr_yearAD()==1980)then
!      write(116,*)I+J/24.,ThetaZ,VOLW0,VLitR
!      endif
    ELSE
      VOLWZ=0.0_r8
    ENDIF
  ELSE
!     non-surface layer
    KL     = micpar%jcplx
    ThetaZ = AZMAX1((AMIN1(AMAX1(0.5_r8*POROS,FieldCapacity),THETW)-THETY))
    VOLWZ  = ThetaZ*VLSoilMicP    
  ENDIF


!     TKS=soil temperature
!     TempOffset=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TSensGrowth,TSensMaintR=temperature function for growth,maintenance respiration
! the offset could be micobial guild/group specific
  TKSO=TKS+TempOffset

  call MicrobPhysTempFun(TKSO, TSensGrowth, TSensMaintR)
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
!     TSolidOMC=total SOC, TSolidOMActC=total colonized SOC
!     TOMBioResdu=total microbial residue, TSorbedOMC=total adsorbed C
!     in each K:
!     SolidOMCK=total SOC n each K, SolidOMActK=total colonized SOC
!     OMBioResduK=total microbial residue, OHCT=total adsorbed C
!
  TSolidOMC    = 0.0_r8
  TSolidOMActC = 0.0_r8
  TOMBioResdu  = 0.0_r8
  TSorbedOMC   = 0.0_r8
!
!     TOTAL SOLID SUBSTRATE
!
  DO  K=1,KL
    SolidOMCK(K)   =0.0_r8
    SolidOMActK(K) =0.0_r8
    DO M=1,jsken
      SolidOMCK(K)   = SolidOMCK(K)+SolidOM(ielmc,M,K)
      SolidOMActK(K) = SolidOMActK(K)+SolidOMAct(M,K)
    enddo
    TSolidOMC    = TSolidOMC+SolidOMCK(K)
    TSolidOMActC = TSolidOMActC+SolidOMActK(K)
  enddo
!
!     TOTAL BIORESIDUE
!
  DO  K=1,KL
    OMBioResduK(K)=0.0_r8
    DO  M=1,ndbiomcp
      OMBioResduK(K)=OMBioResduK(K)+OMBioResdu(ielmc,M,K)
    ENDDO
    TOMBioResdu=TOMBioResdu+OMBioResduK(K)
!
!     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
!
!     BulkSOMC=total SOC
!
    TSorbedOMC=TSorbedOMC+SorbedOM(ielmc,K)+SorbedOM(idom_acetate,K)
  enddo

  D860: DO K=1,KL
    BulkSOMC(K)=SolidOMActK(K)+OMBioResduK(K)+SorbedOM(ielmc,K)+SorbedOM(idom_acetate,K)
  ENDDO D860
  
  TBulkOMC=TSolidOMActC+TOMBioResdu+TSorbedOMC
!
!     C:N AND C:P RATIOS OF TOTAL BIOMASS
!     CNOMA,CPOMA=N,P contents of active biomass OMA
!     FCN,FCP=effects of N,P limitations on biomass activity
!
  TotActMicrobiom     = 0.0_r8
  TotBiomNO2Consumers = 0.0_r8
  D890: DO K = 1, jcplx
    IF(.not.litrm .OR. (K.NE.k_POM .AND. K.NE.k_humus))THEN
      ! the omb complexes, three biomass components, labile, recalcitrant and reserve
      D895: DO N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(n),JGnfo(n)
          MID1=micpar%get_micb_id(1,NGL)
          IF(mBiomeHeter(ielmc,MID1,K).GT.ZEROS)THEN
            rCNBiomeActHeter(ielmn,NGL,K)=AZMAX1(mBiomeHeter(ielmn,MID1,K)/mBiomeHeter(ielmc,MID1,K))
            rCNBiomeActHeter(ielmp,NGL,K)=AZMAX1(mBiomeHeter(ielmp,MID1,K)/mBiomeHeter(ielmc,MID1,K))
          ELSE
            rCNBiomeActHeter(ielmn,NGL,K)=rNCOMC(1,NGL,K)
            rCNBiomeActHeter(ielmp,NGL,K)=rPCOMC(1,NGL,K)
          ENDIF
          OMActHeter(NGL,K)           = AZMAX1(mBiomeHeter(ielmc,MID1,K)/FL(1))
          FCN(NGL,K)                  = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActHeter(ielmn,NGL,K)/rNCOMC(1,NGL,K))))
          FCP(NGL,K)                  = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActHeter(ielmp,NGL,K)/rPCOMC(1,NGL,K))))
          FBiomStoiScalarHeter(NGL,K) = AMIN1(FCN(NGL,K),FCP(NGL,K))

!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
          TotActMicrobiom=TotActMicrobiom+OMActHeter(NGL,K)
          IF(N.EQ.mid_Facult_DenitBacter)THEN
            TotBiomNO2Consumers=TotBiomNO2Consumers+OMActHeter(NGL,K)
          ENDIF
          MID2=micpar%get_micb_id(2,NGL)
          OMC2(NGL,K)=AZMAX1(AMIN1(OMActHeter(NGL,K)*FL(2),mBiomeHeter(ielmc,MID2,K)))
          IF(mBiomeHeter(ielmc,MID2,K).GT.ZEROS)THEN
            FOM2(NGL,K)=AZMAX1(OMC2(NGL,K)/mBiomeHeter(ielmc,MID2,K))
            OMN2(NGL,K)=AZMAX1(FOM2(NGL,K)*mBiomeHeter(ielmn,MID2,K))
          ELSE
            FOM2(NGL,K)=0.0_r8
            OMN2(NGL,K)=0.0_r8
          ENDIF
        ENDDO
      ENDDO D895
    ENDIF
  ENDDO D890

! the abstract complex
  DO N=1,NumMicbFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        MID1=micpar%get_micb_id(1,NGL)
        IF(mBiomeAutor(ielmc,MID1).GT.ZEROS)THEN
          rCNBiomeActAutor(ielmn,NGL)=AZMAX1(mBiomeAutor(ielmn,MID1)/mBiomeAutor(ielmc,MID1))
          rCNBiomeActAutor(ielmp,NGL)=AZMAX1(mBiomeAutor(ielmp,MID1)/mBiomeAutor(ielmc,MID1))
        ELSE
          rCNBiomeActAutor(ielmn,NGL)=rNCOMCAutor(1,NGL)
          rCNBiomeActAutor(ielmp,NGL)=rPCOMCAutor(1,NGL)
        ENDIF
        OMActAutor(NGL)           = AZMAX1(mBiomeAutor(ielmc,MID1)/FL(1))
        FCNAutor(NGL)             = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActAutor(ielmn,NGL)/rNCOMCAutor(1,NGL))))
        FCPAutor(NGL)             = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActAutor(ielmp,NGL)/rPCOMCAutor(1,NGL))))
        FBiomStoiScalarAutor(NGL) = AMIN1(FCNAutor(NGL),FCPAutor(NGL))
!
!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
        TotActMicrobiom=TotActMicrobiom+OMActAutor(NGL)

        IF(N.EQ.mid_AmmoniaOxidBacter)THEN
          TotBiomNO2Consumers=TotBiomNO2Consumers+OMActAutor(NGL)
        ENDIF
        MID2=micpar%get_micb_id(2,NGL)
        OMC2Autor(NGL)=AZMAX1(AMIN1(OMActAutor(NGL)*FL(2),mBiomeAutor(ielmc,MID2)))
        IF(mBiomeAutor(ielmc,MID2).GT.ZEROS)THEN
          FOM2Autor(NGL)=AZMAX1(OMC2Autor(NGL)/mBiomeAutor(ielmc,MID2))
          OMN2Autor(NGL)=AZMAX1(FOM2Autor(NGL)*mBiomeAutor(ielmn,MID2))
        ELSE
          FOM2Autor(NGL)=0.0_r8
          OMN2Autor(NGL)=0.0_r8
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  tOMActC=0._r8
  D690: DO K=1,KL
    TOMK(K)           = 0.0_r8
    TONK(K)           = 0.0_r8
    TOPK(K)           = 0.0_r8
    tMaxNActMicrbK(K) = 0.0_r8
    tMaxPActMicrbK(K) = 0.0_r8
    D685: DO N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        if(OMActHeter(NGL,K)>ZEROS)THEN
          TOMK(K)           = TOMK(K)+OMActHeter(NGL,K)
          TONK(K)           = TONK(K)+OMActHeter(NGL,K)*rCNBiomeActHeter(ielmn,NGL,K)
          TOPK(K)           = TOPK(K)+OMActHeter(NGL,K)*rCNBiomeActHeter(ielmp,NGL,K)
          tMaxNActMicrbK(K) = tMaxNActMicrbK(K)+OMActHeter(NGL,K)*rNCOMC(1,NGL,K)   !maximum total N in active micb
          tMaxPActMicrbK(K) = tMaxPActMicrbK(K)+OMActHeter(NGL,K)*rPCOMC(1,NGL,K)   !maximum total P in active micb
        ENDIF
      ENDDO
    ENDDO D685
    tOMActC=tOMActC+TOMK(K)
  ENDDO D690

  K                 = jcplx+1
  TOMEAutoK(:)      = 0._r8
  tMaxNActMicrbK(K) = 0.0_r8
  tMaxPActMicrbK(K) = 0.0_r8
  DO N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      if(OMActAutor(NGL)>ZEROS)then
        TOMEAutoK(ielmc)  = TOMEAutoK(ielmc)+OMActAutor(NGL)
        TOMEAutoK(ielmn)  = TOMEAutoK(ielmn)+OMActAutor(NGL)*rCNBiomeActAutor(ielmn,NGL)
        TOMEAutoK(ielmp)  = TOMEAutoK(ielmp)+OMActAutor(NGL)*rCNBiomeActAutor(ielmp,NGL)
        tMaxNActMicrbK(K) = tMaxNActMicrbK(K)+OMActAutor(NGL)*rNCOMCAutor(1,NGL)   !maximum total N in active micb
        tMaxPActMicrbK(K) = tMaxPActMicrbK(K)+OMActAutor(NGL)*rPCOMCAutor(1,NGL)   !maximum total P in active micb
      endif
    ENDDO
  ENDDO

!
!     FracBulkSOMC=fraction of total SOC in each substrate complex K
!
  D790: DO K=1,KL
    IF(TBulkOMC.GT.ZEROS)THEN
      FracBulkSOMC(K)=BulkSOMC(K)/TBulkOMC
    ELSE
      FracBulkSOMC(K)=1.0_r8
    ENDIF
    !
    !     DOC CONCENTRATIONS
    !
    !     COQC,COQA=aqueous DOC,acetate concentrations
    !     VLWatMicPM=soil water content, FracBulkSOMC=fraction of total SOC
    !     occupied by each substrate complex K
    !
    IF(VLWatMicPM(NPH).GT.ZEROS2)THEN
      IF(FracBulkSOMC(K).GT.ZERO)THEN
        CDOM(idom_doc,K)     = AZMAX1(DOM(idom_doc,K)/(VLWatMicPM(NPH)*FracBulkSOMC(K)))
        CDOM(idom_acetate,K) = AZMAX1(DOM(idom_acetate,K)/(VLWatMicPM(NPH)*FracBulkSOMC(K)))
      ELSE
        CDOM(idom_doc,K)     = AZMAX1(DOM(idom_doc,K)/VLWatMicPM(NPH))
        CDOM(idom_acetate,K) = AZMAX1(DOM(idom_acetate,K)/VLWatMicPM(NPH))
      ENDIF
    ELSE
      CDOM(idom_doc,K)     = 0.0_r8
      CDOM(idom_acetate,K) = 0.0_r8
    ENDIF
!
!     rCNDOM,rCPDOM=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
!
    IF(DOM(idom_doc,K).GT.ZEROS)THEN
      rCNDOM(K)=AZMAX1(DOM(idom_don,K)/DOM(idom_doc,K))
      rCPDOM(K)=AZMAX1(DOM(idom_dop,K)/DOM(idom_doc,K))
    ELSE
      rCNDOM(K)=0.0_r8
      rCPDOM(K)=0.0_r8
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
  subroutine GetMicrobDensFactorHeter(N,K,micfor, micstt, ORGCL, SPOMK, RMOMK)

  implicit none  
  integer, intent(in) :: N,K
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt    
  real(r8), intent(in) :: ORGCL
  real(r8), intent(out) :: SPOMK(2)
  real(r8), intent(out)  :: RMOMK(2)    

  integer :: NGL,M,MID  
  real(r8) :: TOMCNK(2),COMC

  associate(                             &
    ZEROS       => micfor%ZEROS,         &
    mBiomeHeter => micstt%mBiomeHeter    &
  )

  TOMCNK(:)=0.0_r8
  DO NGL=JGnio(N),JGnfo(N)
    DO M=1,2
      MID=micpar%get_micb_id(M,NGL)          
      TOMCNK(M)=TOMCNK(M)+mBiomeHeter(ielmc,MID,K)
    ENDDO
  ENDDO  

!     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
!
!     COMC=microbial C concentration relative to substrate
!     SPOMK=effect of microbial C concentration on microbial decay
!     RMOMK=effect of microbial C concentration on maintenance respn
!

  IF(ORGCL.GT.ZEROS)THEN
    DO M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
    ENDDO
  ELSE
    DO M=1,2
      SPOMK(M)=1.0_r8
      RMOMK(M)=1.0_r8
    ENDDO
  ENDIF
  end associate
  end subroutine GetMicrobDensFactorHeter
!------------------------------------------------------------------------------------------
  subroutine GetMicrobDensFactorAutor(N, micfor, micstt, ORGCL, SPOMK, RMOMK)

  implicit none  
  integer, intent(in) :: N
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt  
  real(r8), intent(in) :: ORGCL
  real(r8), intent(out) :: SPOMK(2)
  real(r8), intent(out)  :: RMOMK(2)    

  integer :: NGL,M,MID  
  real(r8) :: TOMCNK(2),COMC

  associate(                             &
    ZEROS       => micfor%ZEROS,         &
    mBiomeAutor => micstt%mBiomeAutor    &
  )

  TOMCNK(:)=0.0_r8
  DO NGL=JGnio(N),JGnfo(N)
    DO M=1,2
      MID=micpar%get_micb_id(M,NGL)          
      TOMCNK(M)=TOMCNK(M)+mBiomeAutor(ielmc,MID)
    ENDDO
  ENDDO  

!     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
!
!     COMC=microbial C concentration relative to substrate
!     SPOMK=effect of microbial C concentration on microbial decay
!     RMOMK=effect of microbial C concentration on maintenance respn
!

  IF(ORGCL.GT.ZEROS)THEN
    DO M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
    ENDDO
  ELSE
    DO M=1,2
      SPOMK(M)=1.0_r8
      RMOMK(M)=1.0_r8
    ENDDO
  ENDIF
  end associate
  end subroutine GetMicrobDensFactorAutor
!------------------------------------------------------------------------------------------

  subroutine ActiveMicrobes(I,J,micfor,micstt,micflx,nmicdiag,naqfdiag, &
    nmicf, nmics,ncplxf,ncplxs)
  !
  !  Description:
  !
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout):: ncplxs
  integer :: K,M,N,NGL,MID
  real(r8) :: SPOMK(2)
  real(r8) :: RMOMK(2)    

  REAL(R8) :: OXKX
  real(r8) :: ORGCL
! begin_execution
  associate(                                                     &
  TempMaintRHeter           => nmics%TempMaintRHeter,            &
  GrowthEnvScalHeter        => nmics%GrowthEnvScalHeter,         &
  OMActHeter                => nmics%OMActHeter,                 &
  TDOMUptkHeter             => ncplxf%TDOMUptkHeter,             &
  XCO2                      => nmicdiag%XCO2,                    &
  TSensGrowth               => nmicdiag%TSensGrowth,             &
  TotActMicrobiom           => nmicdiag%TotActMicrobiom,         &
  TotBiomNO2Consumers       => nmicdiag%TotBiomNO2Consumers,     &
  RH2UptkAutor              => nmicdiag%RH2UptkAutor,            &
  WatStressMicb             => nmicdiag%WatStressMicb,           &
  TSensMaintR               => nmicdiag%TSensMaintR,             &
  ZNH4T                     => nmicdiag%ZNH4T,                   &
  ZNO3T                     => nmicdiag%ZNO3T,                   &
  ZNO2T                     => nmicdiag%ZNO2T,                   &
  H2P4T                     => nmicdiag%H2P4T,                   &
  H1P4T                     => nmicdiag%H1P4T,                   &
  VOLWZ                     => nmicdiag%VOLWZ,                   &
  GrowthEnvScalAutor        => nmics%GrowthEnvScalAutor,         &
  TSensMaintRAutor          => nmics%TSensMaintRAutor,           &
  OMActAutor                => nmics%OMActAutor,                 &
  k_humus                   => micpar%k_humus,                   &
  k_POM                     => micpar%k_POM,                     &
  mid_Aerob_Fungi           => micpar%mid_Aerob_Fungi,           &
  is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, &
  PSISoilMatricP            => micfor%PSISoilMatricP,            &
  litrm                     => micfor%litrm,                     &
  H1PO4                     => micstt%H1PO4,                     &
  H1POB                     => micstt%H1POB,                     &
  H2PO4                     => micstt%H2PO4,                     &
  H2POB                     => micstt%H2POB,                     &
  ORGC                      => micfor%ORGC,                      &
  SoilMicPMassLayer         => micfor%SoilMicPMassLayer,         &
  mBiomeAutor               => micstt%mBiomeAutor,               &
  mBiomeHeter               => micstt%mBiomeHeter                &
  )

  !Heterotrophic microbes
  TDOMUptkHeter(:,:) = 0.0_r8
  ORGCL              = AMIN1(1.0E+05_r8*SoilMicPMassLayer,ORGC)

  D760: DO K=1,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM.AND.K.NE.k_humus))THEN
      DO  N=1,NumMicbFunGrupsPerCmplx

        call GetMicrobDensFactorHeter(N,K,micfor, micstt, ORGCL,SPOMK,RMOMK)

        DO NGL=JGnio(N),JGnfo(N)        

!           WatStressMicb=water potential (PSISoilMatricP_vr) effect on microbial respiration
!           OXKX=Km for O2 uptake
!           OXKM=Km for heterotrophic O2 uptake set in starts.f
!           GrowthEnvScalHeter=combined temp and water stress effect on growth respiration
!           TempMaintRHeter=temperature effect on maintenance respiration
          IF(N.EQ.mid_Aerob_Fungi)THEN
            WatStressMicb=EXP(0.1_r8*AMAX1(PSISoilMatricP,-500._r8))
          ELSE
            WatStressMicb=EXP(0.2_r8*AMAX1(PSISoilMatricP,-500._r8))
          ENDIF
          OXKX                      = OXKM
          GrowthEnvScalHeter(NGL,K) = TSensGrowth*WatStressMicb
          TempMaintRHeter(NGL,K)    = TSensMaintR
          IF(OMActHeter(NGL,K).GT.0.0_r8)THEN
            call ActiveHeterotrophs(I,J,NGL,N,K,VOLWZ,XCO2,TSensGrowth,WatStressMicb,SPOMK, RMOMK, &
              OXKX,TotActMicrobiom,TotBiomNO2Consumers,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T, &
              micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO D760

! Autotrophic microbes
  RH2UptkAutor               = 0.0_r8
  N                          = micpar%mid_AmmoniaOxidBacter
  nmicf%RTotNH3OxidSoilAutor = SUM(nmicf%RSOxidSoilAutor(JGniA(N):JGnfA(N)))
  nmicf%RTotNH3OxidBandAutor = SUM(nmicf%RSOxidBandAutor(JGniA(N):JGnfA(N)))

  DO  N=1,NumMicbFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN

      call GetMicrobDensFactorAutor(N,micfor, micstt, ORGCL,SPOMK,RMOMK)

      DO NGL=JGniA(N),JGnfA(N)
        WatStressMicb           = EXP(0.2_r8*PSISoilMatricP)
        OXKX                    = OXKA
        GrowthEnvScalAutor(NGL) = TSensGrowth*WatStressMicb
        TSensMaintRAutor(NGL)   = TSensMaintR
        IF(OMActAutor(NGL).GT.0.0_r8)THEN
          call ActiveAutotrophs(I,J,NGL,N,VOLWZ,XCO2,TSensGrowth,WatStressMicb,SPOMK, RMOMK, &
            OXKX,TotActMicrobiom,TotBiomNO2Consumers,RH2UptkAutor,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T, &
            micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  end associate
  end subroutine ActiveMicrobes
!------------------------------------------------------------------------------------------
  subroutine ActiveHeterotrophs(I,J,NGL,N,K,VOLWZ,XCO2,TSensGrowth,WatStressMicb,SPOMK,RMOMK,&
    OXKX,TotActMicrobiom,TotBiomNO2Consumers, ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,&
    micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,K,N
  real(r8), intent(in) :: VOLWZ
  real(r8), intent(in):: OXKX
  real(r8), intent(in):: WatStressMicb
  real(r8), intent(in):: TSensGrowth
  real(r8), intent(in):: XCO2
  real(r8), intent(in) :: TotActMicrobiom,TotBiomNO2Consumers
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  real(r8), intent(in) :: SPOMK(2)
  real(r8), intent(in) :: RMOMK(2)  
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx
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
  real(r8) :: RVOXPA
  real(r8) :: RVOXPB
  real(r8) :: RGrowthRespHeter
  real(r8) :: RMaintDefcitcitHeter
  real(r8) :: RMaintRespHeter

! begin_execution
  associate(                                                 &
    tRGOMP                 => micflx%tRGOMP               ,  & !aerobic oxidation
    OxyLimterHeter         => nmics%OxyLimterHeter,          &  
    FracOMActHeter         => nmics%FracOMActHeter,          &
    FracNO2ReduxHeter      => nmics%FracNO2ReduxHeter,       &
    mid_Facult_DenitBacter => micpar%mid_Facult_DenitBacter, &
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,     &
    OMActHeter             => nmics%OMActHeter,              &
    RO2UptkHeter           => nmicf%RO2UptkHeter,            &
    RespGrossHeter         => nmicf%RespGrossHeter,          &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,       &
    RO2DmndHeter           => nmicf%RO2DmndHeter,            &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,       &
    RNO3ReduxHeterSoil     => nmicf%RNO3ReduxHeterSoil,      &
    RNO3ReduxHeterBand     => nmicf%RNO3ReduxHeterBand,      &
    RNO2ReduxHeterSoil     => nmicf%RNO2ReduxHeterSoil,      &
    RNO2ReduxHeterBand     => nmicf%RNO2ReduxHeterBand,      &
    RN2OReduxHeter         => nmicf%RN2OReduxHeter,          &
    RNOxReduxRespDenitLim  => nmicf%RNOxReduxRespDenitLim,   &
    RH2ProdHeter           => nmicf%RH2ProdHeter,            &
    RNOxReduxRespDenitUlm  => nmicf%RNOxReduxRespDenitUlm,   &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,           &
    RAcettProdHeter        => nmicf%RAcettProdHeter,         &
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,           &
    TOMK                   => ncplxs%TOMK,                   &
    litrm                  => micfor%litrm,                  &
    ZEROS                  => micfor%ZEROS,                  &
    VLSoilPoreMicP         => micfor%VLSoilPoreMicP,         &
    RNO2EcoUptkSoilPrev    => micfor%RNO2EcoUptkSoilPrev     &
  )
! FracOMActHeter,FracNO2ReduxHeter=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActHeter(NGL,K)=OMActHeter(NGL,K)/TotActMicrobiom
  ELSE
    FracOMActHeter(NGL,K)=1.0_r8
  ENDIF

  IF(TotBiomNO2Consumers.GT.ZEROS .and. N.EQ.mid_Facult_DenitBacter)THEN
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
!
! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
! write(*,*)'SubstrateAttenf4Compet'
  call SubstrateAttenf4Compet(NGL,N,K,FOXYX,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA, &
    micfor,naqfdiag,nmicf,nmics,micflx)

! RO2UptkHeter, RO2DmndHeter=O2-limited, O2-unlimited rates of O2 uptake
  RGOMP               = 0.0_r8
  RO2UptkHeter(NGL,K) = 0.0_r8

!
! HETEROTROPHIC BIOMASS RESPIRATION

  IF(micpar%is_aerobic_hetr(N))THEN
!   RESPIRATION BY HETEROTROPHIC AEROBES:
!   N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!   (6)N2 FIXERS

    call AerobicHeterotrophCatabolism(I,J,NGL,N,K,TSensGrowth,WatStressMicb,FOQC,FOQA, &
      ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)

!   write(*,*)'AerobicHeterO2Uptake'
    call AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB, &
      micfor,micstt,nmicf,nmics,micflx)

  ELSEIF(micpar%is_anaerobic_hetr(N))THEN
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

!     write(*,*)'AnaerobAcetogenCatabolism'
    call AnaerobAcetogenCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQC,ECHZ,FGOCP,&
      FGOAP,RGOMP, micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx)
  ELSEIF(N.EQ.micpar%mid_AcetoMethanogArchea)THEN
!     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
!
!     GOMX=acetate effect on energy yield
!     ECHZ=growth respiration efficiency of aceto. methanogenesis
!
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
  IF(micpar%is_aerobic_hetr(N))THEN    
    RespGrossHeter(NGL,K)  = RGOMP*OxyLimterHeter(NGL,K)
    RCO2ProdHeter(NGL,K)   = RespGrossHeter(NGL,K)
    RAcettProdHeter(NGL,K) = 0.0_r8
    RCH4ProdHeter(NGL,K)   = 0.0_r8
    RO2Uptk4RespHeter(NGL,K) = RO2Dmnd4RespHeter(NGL,K)*OxyLimterHeter(NGL,K)
    RH2ProdHeter(NGL,K)      = 0.0_r8
    tRGOMP                   = tRGOMP+RGOMP
    !anaerboic fertmenting heterotrophs
  ELSEIF(micpar%is_anaerobic_hetr(N))THEN
    !fermentation  (CH2O)6 -> 2CO2 + 2(CH2O)2  
    RespGrossHeter(NGL,K)    = RGOMP
    RCO2ProdHeter(NGL,K)     = 0.333_r8*RespGrossHeter(NGL,K)
    RAcettProdHeter(NGL,K)   = 0.667_r8*RespGrossHeter(NGL,K)
    RCH4ProdHeter(NGL,K)     = 0.0_r8
    RO2Uptk4RespHeter(NGL,K) = RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)      = 0.111_r8*RespGrossHeter(NGL,K)

  ELSEIF(N.EQ.micpar%mid_AcetoMethanogArchea)THEN
    ! CH3COOH -> CO2 + CH4
    RespGrossHeter(NGL,K)   = RGOMP
    RCO2ProdHeter(NGL,K)    = 0.50_r8*RespGrossHeter(NGL,K)
    RAcettProdHeter(NGL,K)  = 0.0_r8
    RCH4ProdHeter(NGL,K)    = 0.50_r8*RespGrossHeter(NGL,K)
    RO2Uptk4RespHeter(NGL,K)= RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)     = 0.0_r8

  ENDIF
!
!  write(*,*)'HETEROTROPHIC DENITRIFICATION'
!
  IF(N.EQ.micpar%mid_Facult_DenitBacter .AND. RO2Dmnd4RespHeter(NGL,K).GT.0.0_r8 &
    .AND.(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS))THEN
  !no litter layer denitrifcation
    call HeteroDenitrificCatabolism(NGL,N,K,FOQC,RGOCP, &
      VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  ELSE
    RNO3ReduxHeterSoil(NGL,K)    = 0.0_r8
    RNO3ReduxHeterBand(NGL,K)    = 0.0_r8
    RNO2ReduxHeterSoil(NGL,K)    = 0.0_r8
    RNO2ReduxHeterBand(NGL,K)    = 0.0_r8
    RN2OReduxHeter(NGL,K)        = 0.0_r8
    RNOxReduxRespDenitUlm(NGL,K) = 0.0_r8
    RNOxReduxRespDenitLim(NGL,K) = 0.0_r8
  ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
  call BiomassMineralization(NGL,N,K,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt, &
    nmicf,nmics,micflx)
!
  call GatherHetertrophRespiration(I,J,NGL,N,K,RMOMK,RGrowthRespHeter,RMaintDefcitcitHeter,RMaintRespHeter, &
    micfor,micstt,nmicf,nmics)
!
  call GatherHetertrophAnabolicFlux(I,J,NGL,N,K,ECHZ,FGOCP,FGOAP,RGrowthRespHeter,&
    RMaintDefcitcitHeter,RMaintRespHeter,spomk,rmomk,micfor,micstt,nmicf, &
    nmics,ncplxf,ncplxs,micflx)
  end associate
  end subroutine ActiveHeterotrophs

!------------------------------------------------------------------------------------------

  subroutine ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Cumlate_Flux_Diag_type), INTENT(INOUT):: naqfdiag
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: H_1p_conc,CHNO2,CHNOB
  real(r8) :: FNO3S,FNO3B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: VMXC4S,VMXC4B
!     begin_execution
  associate(                                               &
    TSensGrowth         =>  nmicdiag%TSensGrowth         , &
    RNO2ReduxSoilChemo  =>  naqfdiag%RNO2ReduxSoilChemo  , &
    RNO2ReduxBandChemo  =>  naqfdiag%RNO2ReduxBandChemo  , &
    RN2OProdSoilChemo   =>  naqfdiag%RN2OProdSoilChemo   , &
    RN2OProdBandChemo   =>  naqfdiag%RN2OProdBandChemo   , &
    RNO3ProdSoilChemo   =>  naqfdiag%RNO3ProdSoilChemo   , &
    RNO3ProdBandChemo   =>  naqfdiag%RNO3ProdBandChemo   , &
    RNO2ReduxChemo      =>  naqfdiag%RNO2ReduxChemo      , &
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
  H_1p_conc = AMAX1(ZERO,10.0_r8**(-(PH-3.0_r8)))
  CHNO2     = CNO2S*H_1p_conc/0.5_r8
  CHNOB     = CNO2B*H_1p_conc/0.5_r8

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
  naqfdiag%TFNO2X    = naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B    = naqfdiag%TFNO2B+FNB2
  FNO3S              = VLNO3
  FNO3B              = VLNOB
  VMXC4S             = 7.5E-02_r8*CHNO2*VLWatMicPM(NPH)*FNO3S*TSensGrowth
  VMXC4B             = 7.5E-02_r8*CHNOB*VLWatMicPM(NPH)*FNO3B*TSensGrowth
  RNO2ReduxSoilChemo = AZMAX1(AMIN1(ZNO2S*FNO2,VMXC4S))
  RNO2ReduxBandChemo = AZMAX1(AMIN1(ZNO2B*FNB2,VMXC4B))
  RN2OProdSoilChemo  = 0.10_r8*RNO2ReduxSoilChemo
  RN2OProdBandChemo  = 0.10_r8*RNO2ReduxBandChemo
  RNO3ProdSoilChemo  = 0.80_r8*RNO2ReduxSoilChemo
  RNO3ProdBandChemo  = 0.80_r8*RNO2ReduxBandChemo
  RNO2ReduxChemo     = 0.10_r8*(RNO2ReduxSoilChemo+RNO2ReduxBandChemo)
  RNO2DmndSoilChemo  = VMXC4S
  RNO2DmndBandChemo  = VMXC4B

  end associate
  end subroutine ChemoDenitrification
!------------------------------------------------------------------------------------------

  subroutine OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_Flux_type), intent(inout):: ncplxf
  type(OMCplx_State_type),intent(inout):: ncplxs
  integer  :: K,M,N,KK,NGL,MID,idom,NE
  real(r8) :: OSRT
  real(r8) :: XFRK,XFROM(idom_beg:idom_end)
  real(r8) :: XFME
!     begin_execution
  associate(                                             &
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,    &
    XferBiomeHeterK      => nmicf%XferBiomeHeterK,       &
    ROQC4HeterMicActCmpK => ncplxf%ROQC4HeterMicActCmpK, &
    XferRespHeterK       => ncplxf%XferRespHeterK,       &
    XferDOMK             => ncplxf%XferDOMK,             &
    BulkSOMC             => ncplxs%BulkSOMC,             &
    TMicHeterActivity    => micstt%TMicHeterActivity,    &
    DOM                  => micstt%DOM,                  &
    mBiomeHeter          => micstt%mBiomeHeter,          &
    ZEROS                => micfor%ZEROS,                &
    TScal4Difsvity       => micfor%TScal4Difsvity        &
  )
!
!     BulkSOMC=total SOC in each K
!     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
!     between each K and KK, FPRIM=priming transfer rate constant
!     TScal4Difsvity=temperature effect on priming transfers
!     ROQC4HeterMicActCmpK,OQC,OQN,OQP=respiration,DOC,DON,DOP
!     XferRespHeterK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
!
  D795: DO K=1,KL
    IF(K.LE.KL-1)THEN
      D800: DO KK=K+1,KL
        OSRT=BulkSOMC(K)+BulkSOMC(KK)
        IF(BulkSOMC(K).GT.ZEROS.AND.BulkSOMC(KK).GT.ZEROS)THEN
          XFRK=FPRIM*TScal4Difsvity*(ROQC4HeterMicActCmpK(K)*BulkSOMC(KK)-ROQC4HeterMicActCmpK(KK)*BulkSOMC(K))/OSRT
          DO idom=idom_beg,idom_end
            XFROM(idom)=FPRIM*TScal4Difsvity*(DOM(idom,K)*BulkSOMC(KK)-DOM(idom,KK)*BulkSOMC(K))/OSRT
          ENDDO
          IF(ROQC4HeterMicActCmpK(K)+XferRespHeterK(K)-XFRK.GT.0.0_r8 .AND. ROQC4HeterMicActCmpK(KK)+XferRespHeterK(KK)+XFRK.GT.0.0_r8)THEN
            XferRespHeterK(K)  = XferRespHeterK(K)-XFRK
            XferRespHeterK(KK) = XferRespHeterK(KK)+XFRK
          ENDIF
          DO iDOM=idom_beg,idom_end
            IF(DOM(idom,K)+XferDOMK(idom,K)-XFROM(idom).GT.0.0_r8 &
              .AND.DOM(idom,KK)+XferDOMK(idom,KK)+XFROM(idom).GT.0.0_r8)THEN
              XferDOMK(idom,K)  = XferDOMK(idom,K)-XFROM(idom)
              XferDOMK(idom,KK) = XferDOMK(idom,KK)+XFROM(idom)
            ENDIF
          ENDDO
!
!     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
!
!     XFMC,XFMN,XFMP=transfer of microbial C,N,P
!     between each K and KK, FPRIMM=priming transfer rate constant
!     GrowthEnvScalHeter=temperature+water effect
!     OMC,OMN,OMP=microbial C,N,P
!     BulkSOMC=total SOC in each K
!     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
!
          D850: DO N=1,NumMicbFunGrupsPerCmplx
            DO  M=1,nlbiomcp
              DO NGL=JGnio(N),JGnfo(N)
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  XFME=FPRIMM*GrowthEnvScalHeter(NGL,K)*(mBiomeHeter(NE,MID,K)*BulkSOMC(KK) &
                    -mBiomeHeter(NE,MID,KK)*BulkSOMC(K))/OSRT
                  IF(mBiomeHeter(NE,MID,K)+XferBiomeHeterK(NE,M,NGL,K)-XFME.GT.0.0_r8 &
                    .AND.mBiomeHeter(NE,MID,KK)+XferBiomeHeterK(NE,M,NGL,KK)+XFME.GT.0.0_r8)THEN
                    XferBiomeHeterK(NE,M,NGL,K)  = XferBiomeHeterK(NE,M,NGL,K)-XFME
                    XferBiomeHeterK(NE,M,NGL,KK) = XferBiomeHeterK(NE,M,NGL,KK)+XFME
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
!     TMicHeterActivity=total respiration of DOC+DOA in soil layer
!     ROQC4HeterMicActCmpK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
  TMicHeterActivity=0.0_r8
  D840: DO K=1,KL
    ROQC4HeterMicActCmpK(K) = ROQC4HeterMicActCmpK(K)+XferRespHeterK(K)
    TMicHeterActivity       = TMicHeterActivity+ROQC4HeterMicActCmpK(K)
    DO idom=idom_beg,idom_end
      DOM(idom,K)=DOM(idom,K)+XferDOMK(idom,K)
    ENDDO
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        do NGL=JGnio(N),JGnfo(N)
          MID=micpar%get_micb_id(M,NGL)        
          DO NE=1,NumPlantChemElms
            mBiomeHeter(NE,MID,K)=mBiomeHeter(NE,MID,K)+XferBiomeHeterK(NE,M,NGL,K)
          ENDDO
        enddo
      enddo
    enddo
  ENDDO D840
  end associate
  end subroutine OMTransferForPriming
!------------------------------------------------------------------------------------------

  subroutine RDOMSorption(KL,micfor,micstt,nmicf,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: KL                 !from 1 to 5
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  real(r8) :: AECX
  real(r8) :: OQEX(idom_beg:idom_end)
  real(r8) :: OHEX(idom_beg:idom_end)
  integer :: idom,K
  real(r8) :: VLSoilPoreMicPX
  real(r8) :: VLSoilPoreMicPW,VOLCX,VOLCW,VOLAX,VOLAW
!     begin_execution
  associate(                                            &
    RDOMSorp             => ncplxf%RDOMSorp,            &
    TDOMUptkHeter        => ncplxf%TDOMUptkHeter,       &
    BulkSOMC             => ncplxs%BulkSOMC,            &
    FOCA                 => ncplxs%FOCA,                &
    FOAA                 => ncplxs%FOAA,                &
    SoilMicPMassLayer    => micfor%SoilMicPMassLayer,   &
    ZERO                 => micfor%ZERO,                &
    ZEROS2               => micfor%ZEROS2,              &
    ZEROS                => micfor%ZEROS,               &
    litrm                => micfor%litrm,               &
    VLWatMicPM           => micfor%VLWatMicPM,          &
    FracBulkSOMC         => micstt%FracBulkSOMC,        &
    DOM                  => micstt%DOM,                 &
    SorbedOM             => micstt%SorbedOM,            &
    AEC                  => micfor%AEC                  &
  )
!     VLWatMicPM=soil water content, FracBulkSOMC=fraction of total SOC
!     AEC,AECX=anion exchange capacity
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     TCGOQC,TDOMUptkHeter,TDOMUptkHeter,TCGOAC=total uptake of DOC,DON,DOP,acetate
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     TSORP,HSORP=sorption rate constant and coefficient for OHC
!     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
  DO K=1,KL
    IF(VLWatMicPM(NPH).GT.ZEROS2.AND.FracBulkSOMC(K).GT.ZERO)THEN
      IF(litrm)THEN
        AECX=0.5E+03_r8
      ELSE
        AECX=AEC
      ENDIF
      DO idom=idom_beg,idom_end
        OQEX(idom)=AMAX1(ZEROS,DOM(idom,K)-TDOMUptkHeter(idom,K))  !free DOM
        OHEX(idom)=AMAX1(ZEROS,SorbedOM(idom,K))                   !adsorbed DOM
      ENDDO

      VLSoilPoreMicPX=SoilMicPMassLayer*AECX*HSORP*FracBulkSOMC(K)
      VLSoilPoreMicPW=VLWatMicPM(NPH)*FracBulkSOMC(K)
      IF(FOCA(K).GT.ZERO)THEN
        VOLCX                = FOCA(K)*VLSoilPoreMicPX
        VOLCW                = FOCA(K)*VLSoilPoreMicPW
        RDOMSorp(idom_doc,K) = TSORP*(OQEX(idom_doc)*VOLCX-OHEX(idom_doc)*VOLCW)/(VOLCX+VOLCW)
      ELSE
        RDOMSorp(idom_doc,K)=TSORP*(OQEX(idom_doc)*VLSoilPoreMicPX  &
          -OHEX(idom_doc)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      ENDIF

      IF(FOAA(K).GT.ZERO)THEN
        VOLAX                    = FOAA(K)*VLSoilPoreMicPX
        VOLAW                    = FOAA(K)*VLSoilPoreMicPW
        RDOMSorp(idom_acetate,K) = TSORP*(OQEX(idom_acetate)*VOLAX-OHEX(idom_acetate)*VOLAW)/(VOLAX+VOLAW)
      ELSE
        RDOMSorp(idom_acetate,K)=TSORP*(OQEX(idom_acetate)*VLSoilPoreMicPX &
          -OHEX(idom_acetate)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      ENDIF
      RDOMSorp(idom_don,K)=TSORP*(OQEX(idom_don)*VLSoilPoreMicPX &
        -OHEX(idom_don)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      RDOMSorp(idom_dop,K)=TSORP*(OQEX(idom_dop)*VLSoilPoreMicPX &
        -OHEX(idom_dop)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
    ELSE
      DO idom=idom_beg,idom_end
        RDOMSorp(idom,K)=0.0_r8
      enddo
    ENDIF
  ENDDO
  end associate
  end subroutine RDOMSorption
!------------------------------------------------------------------------------------------

  subroutine SolidOMDecomposition(I,J,KL,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs,micflx)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx  
  integer  :: M,NE,idom,K
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
  associate(                                             &
    tRHydlySOM           => micflx%tRHydlySOM,           &
    tRHydlyBioReSOM      => micflx%tRHydlyBioReSOM,      &
    tRHydlySoprtOM       => micflx%tRHydlySoprtOM,       &
    RHydlysSolidOM       => ncplxf%RHydlysSolidOM,       &
    RHumifySolidOM       => ncplxf%RHumifySolidOM,       &
    RDcmpProdDOM         => ncplxf%RDcmpProdDOM,         &
    RHydlysBioResduOM    => ncplxf%RHydlysBioResduOM,    &
    RHydlysSorptOM       => ncplxf%RHydlysSorptOM,       &
    ROQC4HeterMicActCmpK => ncplxf%ROQC4HeterMicActCmpK, &
    BulkSOMC             => ncplxs%BulkSOMC,             &
    TOMK                 => ncplxs%TOMK,                 &
    TONK                 => ncplxs%TONK,                 &
    TOPK                 => ncplxs%TOPK,                 &
    rCNSorbOM            => ncplxs%rCNSorbOM,            &
    rCPSorbOM            => ncplxs%rCPSorbOM,            &
    tMaxNActMicrbK       => ncplxs%tMaxNActMicrbK,       &
    tMaxPActMicrbK       => ncplxs%tMaxPActMicrbK,       &
    VOLWZ                => nmicdiag%VOLWZ,              &
    TSensGrowth          => nmicdiag%TSensGrowth,        &
    iprotein             => micpar%iprotein,             &
    icarbhyro            => micpar%icarbhyro,            &
    icellulos            => micpar%icellulos,            &
    ilignin              => micpar%ilignin,              &
    SPOSC                => micpar%SPOSC,                &
    k_POM                => micpar%k_POM,                &
    CNRH                 => micpar%CNRH,                 &
    CPRH                 => micpar%CPRH,                 &
    CDOM                 => ncplxs%CDOM,                 &
    EPOC                 => micstt%EPOC,                 &
    CNOSC                => micstt%CNOSC,                &
    CPOSC                => micstt%CPOSC,                &
    SorbedOM             => micstt%SorbedOM,             &
    SolidOM              => micstt%SolidOM,              &
    SolidOMAct           => micstt%SolidOMAct,           &
    OMBioResdu           => micstt%OMBioResdu,           &
    VLSoilMicP           => micfor%VLSoilMicP,           &
    ZEROS                => micfor%ZEROS,                &
    ZEROS2               => micfor%ZEROS2,               &
    litrm                => micfor%litrm,                &
    SoilMicPMassLayer    => micfor%SoilMicPMassLayer     &
  )
!     FCPK=N,P limitation to microbial activity in each K
!     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
!     COQCK=aqueous concentration of microbial activity
!     DCKD=Km for decomposition of SOC at current COQCK
!     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
!     DCKI=inhibition of decomposition by microbial concentration
!     BulkSOMC=total SOC
!     COSC=concentration of total SOC
!     SoilMicPMassLayer,VLSoilPoreMicP=mass, volume of soil layer
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     OQKI=DOC product inhibition constant for decomposition
!
  DO K=1,KL
    IF(TOMK(K).GT.ZEROS)THEN
      CNOMX   = TONK(K)/tMaxNActMicrbK(K)
      CPOMX   = TOPK(K)/tMaxPActMicrbK(K)
      FCNK(K) = AMIN1(1.0_r8,AMAX1(0.50_r8,CNOMX))
      FCPK(K) = AMIN1(1.0_r8,AMAX1(0.50_r8,CPOMX))
    ELSE
      FCNK(K) = 1.0_r8
      FCPK(K) = 1.0_r8
    ENDIF
  !
  !     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
  !     CONSTANT FOR DECOMPOSITION
  !
    IF(VOLWZ.GT.ZEROS2)THEN
      COQCK=AMIN1(0.1E+06_r8,ROQC4HeterMicActCmpK(K)/VOLWZ)
    ELSE
      COQCK=0.1E+06_r8
    ENDIF
    IF(litrm)THEN
      DCKD=DCKM0*(1.0_r8+COQCK/DCKI)
    ELSE
      DCKD=DCKML*(1.0_r8+COQCK/DCKI)
    ENDIF

    IF(BulkSOMC(K).GT.ZEROS)THEN
      IF(SoilMicPMassLayer.GT.ZEROS)THEN
        COSC=BulkSOMC(K)/SoilMicPMassLayer
      ELSE
        COSC=BulkSOMC(K)/VLSoilMicP
      ENDIF
      DFNS = COSC/(COSC+DCKD)
      OQCI = 1.0_r8/(1.0_r8+CDOM(idom_doc,K)/OQKI)
  !
  !     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES 'RDOS*' FROM
  !     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
  !     TEMPERATURE, SUBSTRATE C:N, C:P
  !
  !     CNS,CPS=N:C,P:C ratios of SOC
  !     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
  !     OSA,OSN,OSP=active biomass C,N,P
  !     SPOSC=specific decomposition rate constant
  !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
  !     DFNS=effect of microbial concentration on decomposition
  !     OQCI=DOC product inhibition for decomposition
  !     TSensGrowth=temperature stress effect
  !     BulkSOMC=total SOC
  !     FCNK,FCPK=N,P limitation to microbial activity in each K
  !
!      if(etimer%get_curr_yearAD()==1980 .and. litrm)then
!      write(115,*)I+J/24.,K,ROQC4HeterMicActCmpK(K),VOLWZ,COSC,COQCK,CDOM(idom_doc,K),ZEROS2
!      endif
      D785: DO M=1,jsken
        IF(SolidOM(ielmc,M,K).GT.ZEROS)THEN
          CNS(M,K)                  = AZMAX1(SolidOM(ielmn,M,K)/SolidOM(ielmc,M,K))
          CPS(M,K)                  = AZMAX1(SolidOM(ielmp,M,K)/SolidOM(ielmc,M,K))
          RHydlysSolidOM(ielmc,M,K) = SolidOMAct(M,K)*AZMAX1(AMIN1(0.5_r8 &
            ,SPOSC(M,K)*ROQC4HeterMicActCmpK(K)*DFNS*OQCI*TSensGrowth/BulkSOMC(K)))
          RHydlysSolidOM(ielmn,M,K)=AZMAX1(AMIN1(SolidOM(ielmn,M,K),CNS(M,K)*RHydlysSolidOM(ielmc,M,K)))/FCNK(K)
          RHydlysSolidOM(ielmp,M,K)=AZMAX1(AMIN1(SolidOM(ielmp,M,K),CPS(M,K)*RHydlysSolidOM(ielmc,M,K)))/FCPK(K)
!          IF(litrm)THEN
!            write(113,*)I+J/24.,ROQC4HeterMicActCmpK(K)/TOMK(K),ncplxs%CDOM(idom_doc,K),TOMK(K)
!          endif

          DO NE=1,NumPlantChemElms
            tRHydlySOM(NE)=tRHydlySOM(NE)+RHydlysSolidOM(NE,M,K)
          ENDDO

        ELSE
          CNS(M,K) = CNOSC(M,K)
          CPS(M,K) = CPOSC(M,K)
          DO NE    = 1, NumPlantChemElms
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
        RHOSCM                           = 0.10_r8*RHumifySolidOM(ielmc,ilignin,K)
        RHumifySolidOM(ielmc,iprotein,K) = AZMAX1(AMIN1(RHydlysSolidOM(ielmc,iprotein,K) &
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
            RHumifySolidOM(NE,M,K) = 0.0_r8
            RDcmpProdDOM(NE,M,K)   = RHydlysSolidOM(NE,M,K)
          ENDDO
        ENDDO D810
      ENDIF

    ELSE
      D780: DO M=1,jsken
        DO NE=1,NumPlantChemElms    
          RHydlysSolidOM(NE,M,K) = 0.0_r8
          RHumifySolidOM(NE,M,K) = 0.0_r8
          RDcmpProdDOM(NE,M,K)   = 0.0_r8
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
  !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
  !     DFNS=effect of microbial concentration on decomposition
  !     OQCI=DOC product inhibition for decomposition
  !     TSensGrowth=temperature stress effect
  !     BulkSOMC=total SOC
  !     FCNK,FCPK=N,P limitation to microbial activity in each K
  !
    IF(BulkSOMC(K).GT.ZEROS)THEN
      D775: DO M=1,ndbiomcp
        IF(OMBioResdu(ielmc,M,K).GT.ZEROS)THEN
          CNR                          = AZMAX1(OMBioResdu(ielmn,M,K)/OMBioResdu(ielmc,M,K))
          CPR                          = AZMAX1(OMBioResdu(ielmp,M,K)/OMBioResdu(ielmc,M,K))
          RHydlysBioResduOM(ielmc,M,K) = OMBioResdu(ielmc,M,K)*AZMAX1(AMIN1(1._r8 &
            ,SPORC(M)*ROQC4HeterMicActCmpK(K)*DFNS*OQCI*TSensGrowth/BulkSOMC(K)))
          RHydlysBioResduOM(ielmn,M,K) = AZMAX1(AMIN1(OMBioResdu(ielmn,M,K),CNR*RHydlysBioResduOM(ielmc,M,K)))/FCNK(K)
          RHydlysBioResduOM(ielmp,M,K) = AZMAX1(AMIN1(OMBioResdu(ielmp,M,K),CPR*RHydlysBioResduOM(ielmc,M,K)))/FCPK(K)

          DO NE=1,NumPlantChemElms
            tRHydlyBioReSOM(NE)=tRHydlyBioReSOM(NE)+RHydlysBioResduOM(NE,M,K)
          ENDDO
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
  !     rCNSorbOM,rCPSorbOM=N:C,P:C ratios of adsorbed C,N,P
  !     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
  !     SPOHC=specific decomposition rate constant for adsorbed C
  !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
  !     DFNS=effect of microbial concentration on decomposition
  !     OQCI=DOC product inhibition for decomposition
  !     TSensGrowth=temperature stress effect
  !     BulkSOMC=total SOC
  !     FCNK,FCPK=N,P limitation to microbial activity in each K
  !
    IF(BulkSOMC(K).GT.ZEROS)THEN
      IF(SorbedOM(ielmc,K).GT.ZEROS)THEN
        rCNSorbOM(K)                   = AZMAX1(SorbedOM(ielmn,K)/SorbedOM(ielmc,K))
        rCPSorbOM(K)                   = AZMAX1(SorbedOM(ielmp,K)/SorbedOM(ielmc,K))
        RHydlysSorptOM(ielmc,K)        = SorbedOM(ielmc,K)*AZMAX1(AMIN1(1._r8 &
          , SPOHC*ROQC4HeterMicActCmpK(K)*DFNS*OQCI*TSensGrowth/BulkSOMC(K)))
        RHydlysSorptOM(idom_acetate,K) = SorbedOM(idom_acetate,K)*AZMAX1(AMIN1(1._r8 &
          , SPOHA*ROQC4HeterMicActCmpK(K)*DFNS*TSensGrowth/BulkSOMC(K)))

        RHydlysSorptOM(ielmn,K)        = AZMAX1(AMIN1(SorbedOM(ielmn,K),rCNSorbOM(K)*RHydlysSorptOM(ielmc,K)))/FCNK(K)
        RHydlysSorptOM(ielmp,K)        = AZMAX1(AMIN1(SorbedOM(ielmp,K),rCPSorbOM(K)*RHydlysSorptOM(ielmc,K)))/FCPK(K)

        DO NE=1,NumPlantChemElms
          tRHydlySoprtOM(NE)=tRHydlySoprtOM(NE)+RHydlysSorptOM(NE,K)
        ENDDO

      ELSE
        rCNSorbOM(K)=0.0_r8
        rCPSorbOM(K)=0.0_r8
        DO idom=idom_beg,idom_end
          RHydlysSorptOM(idom,K)=0.0_r8
        ENDDO
      ENDIF
    ELSE
      rCNSorbOM(K)=0.0_r8
      rCPSorbOM(K)=0.0_r8
      DO idom=idom_beg,idom_end          
        RHydlysSorptOM(idom,K)=0.0_r8
      ENDDO
    ENDIF
  ENDDO
  end associate
  end subroutine SolidOMDecomposition
!------------------------------------------------------------------------------------------

  subroutine RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type), intent(in) :: nmicdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout):: ncplxs
  integer :: K,M,N,NGL,NE,idom
  real(r8) :: FORC(0:jcplx)
!     begin_execution
  associate(                                                                    &
    k_POM                            => micpar%k_POM,                           &
    DOMuptk4GrothHeter               => nmicf%DOMuptk4GrothHeter,               &
    RAnabolDOCUptkHeter              => nmicf%RAnabolDOCUptkHeter,              &
    RAnabolAcetUptkHeter             => nmicf%RAnabolAcetUptkHeter,             &
    RkillLitrfal2ResduOMHeter        => nmicf%RkillLitrfal2ResduOMHeter,        &
    RMaintDefcitLitrfal2ResduOMHeter => nmicf%RMaintDefcitLitrfal2ResduOMHeter, &
    RCCMEheter                       => nmicf%RCCMEheter,                       &
    RAcettProdHeter                  => nmicf%RAcettProdHeter,                  &
    RHydlysSolidOM                   => ncplxf%RHydlysSolidOM,                  &
    RHumifySolidOM                   => ncplxf%RHumifySolidOM,                  &
    RDcmpProdDOM                     => ncplxf%RDcmpProdDOM,                    &
    RHydlysBioResduOM                => ncplxf%RHydlysBioResduOM,               &
    RHydlysSorptOM                   => ncplxf%RHydlysSorptOM,                  &
    RDOMSorp                         => ncplxf%RDOMSorp,                        &
    OMBioResduK                      => ncplxs%OMBioResduK,                     &
    TOMBioResdu                      => nmicdiag%TOMBioResdu,                   &
    RMaintDefcitLitrfal2ResduOMAutor => nmicf%RMaintDefcitLitrfal2ResduOMAutor, &
    RkillLitrfal2ResduOMAutor        => nmicf%RkillLitrfal2ResduOMAutor,        &
    SolidOM                          => micstt%SolidOM,                         &
    iprotein                         => micpar%iprotein,                        &
    SolidOMAct                       => micstt%SolidOMAct,                      &
    SOMPomProtein                    => micstt%SOMPomProtein,                   &
    DOM                              => micstt%DOM,                             &
    OMBioResdu                       => micstt%OMBioResdu,                      &
    SorbedOM                         => micstt%SorbedOM,                        &
    ZEROS                            => micfor%ZEROS,                           &
    Litrm                            => micfor%litrm                            &
  )
!
!     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
!     HETEROTROPHIC SUBSTRATE-MICROBE complexES
!
!     FORC=fraction of total microbial residue
!     OMBioResduK=microbial residue
!     RCCMEheter,RCCMN,RCCMP=transfer of auto LitrFall C,N,P to each hetero K
!     RkillLitrfal2ResduOMHeter,RkillLitrfal2ResduOMHeter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!     RMaintDefcitLitrfal2ResduOMHeter,RCMMN,RMaintDefcitLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
!
  D1690: DO K=1,KL
    IF(TOMBioResdu.GT.ZEROS)THEN
      FORC(K)=OMBioResduK(K)/TOMBioResdu
    ELSE
      IF(K.EQ.k_POM)THEN
        FORC(K)=1.0_r8
      ELSE
        FORC(K)=0.0_r8
      ENDIF
    ENDIF
    D1685: DO N=1,NumMicbFunGrupsPerCmplx
      D1680: DO M=1,ndbiomcp
        DO NGL=JGniA(N),JGnfA(N)
        DO NE=1,NumPlantChemElms
          RCCMEheter(NE,M,NGL,K)=(RkillLitrfal2ResduOMAutor(NE,M,NGL)+RMaintDefcitLitrfal2ResduOMAutor(NE,M,NGL))*FORC(K)
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

!     SolidOMAct(M,K)=SolidOMAct(M,K)-RHydlysSolidOM(ielmc,M,K)
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
 !      SolidOMAct(1,k_POM)=SolidOMAct(1,k_POM)+RHumifySolidOM(ielmc,M,K)
        DO NE=1,NumPlantChemElms
          SolidOM(NE,iprotein,k_POM)=SolidOM(NE,iprotein,k_POM)+RHumifySolidOM(NE,M,K)
        ENDDO
      ELSE
        DO NE=1,NumPlantChemElms
          SOMPomProtein(NE)=SOMPomProtein(NE)+RHumifySolidOM(NE,M,K)
        ENDDO
      ENDIF

    ENDDO D580
!
!     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!
    D575: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms    
        OMBioResdu(NE,M,K)=OMBioResdu(NE,M,K)-RHydlysBioResduOM(NE,M,K)
        DOM(NE,K)=DOM(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D575

    DO idom=idom_beg,idom_end
      DOM(idom,K)      = DOM(idom,K)+RHydlysSorptOM(idom,K)
      SorbedOM(idom,K) = SorbedOM(idom,K)-RHydlysSorptOM(idom,K)
    ENDDO
!
!     MICROBIAL UPTAKE OF DISSOLVED C, N, P
!
!     RAnabolDOCUptkHeter,RAnabolAcetUptkHeter,DOMuptk4GrothHeter,DOMuptk4GrothHeter=DOC,acetate,DON,DOP uptake
!     RAcettProdHeter=acetate production from fermentation
!
    D570: DO N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        DOM(idom_doc,K)     = DOM(idom_doc,K)-RAnabolDOCUptkHeter(NGL,K)
        DOM(idom_don,K)     = DOM(idom_don,K)-DOMuptk4GrothHeter(ielmn,NGL,K)
        DOM(idom_dop,K)     = DOM(idom_dop,K)-DOMuptk4GrothHeter(ielmp,NGL,K)
        DOM(idom_acetate,K) = DOM(idom_acetate,K)-RAnabolAcetUptkHeter(NGL,K)+RAcettProdHeter(NGL,K)
!
!     MICROBIAL DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RkillLitrfal2ResduOMHeter,RkillLitrfal2ResduOMHeter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!     RCCMEheter,RCCMN,RCCMP=transfer of auto LitrFall C,N,P to each hetero K
!     RMaintDefcitLitrfal2ResduOMHeter,RCMMN,RMaintDefcitLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
!
        D565: DO M=1,ndbiomcp
          DO NE=1,NumPlantChemElms
            OMBioResdu(NE,M,K)=OMBioResdu(NE,M,K)+RkillLitrfal2ResduOMHeter(NE,M,NGL,K) &
              +RCCMEheter(NE,M,NGL,K)+RMaintDefcitLitrfal2ResduOMHeter(NE,M,NGL,K)
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
      DOM(idom,K)=DOM(idom,K)-RDOMSorp(idom,K)    
      SorbedOM(idom,K)=SorbedOM(idom,K)+RDOMSorp(idom,K)
    ENDDO

  ENDDO D590
  end associate
  end subroutine RedistDecompProduct
!------------------------------------------------------------------------------------------

  subroutine HeterotrophAnabolicUpdate(I,J,micfor,micstt,nmicf,micflx)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
  integer  :: K,M,N,NGL,MID3,MID,NE
  real(r8) ::CGROMC,dmassC
!     begin_execution
  associate(                                                                &
    DOMuptk4GrothHeter             => nmicf%DOMuptk4GrothHeter,             &
    NonstX2stBiomHeter             => nmicf%NonstX2stBiomHeter,             &
    Resp4NFixHeter                 => nmicf%Resp4NFixHeter,                 &
    RespGrossHeter                 => nmicf%RespGrossHeter,                 &
    RNOxReduxRespDenitLim          => nmicf%RNOxReduxRespDenitLim,          &
    RNO3TransfSoilHeter            => nmicf%RNO3TransfSoilHeter,            &
    RCO2ProdHeter                  => nmicf%RCO2ProdHeter,                  &
    RH2PO4TransfSoilHeter          => nmicf%RH2PO4TransfSoilHeter,          &
    RNH4TransfBandHeter            => nmicf%RNH4TransfBandHeter,            &
    RNO3TransfBandHeter            => nmicf%RNO3TransfBandHeter,            &
    RH2PO4TransfBandHeter          => nmicf%RH2PO4TransfBandHeter,          &
    RkillLitrfal2HumOMHeter        => nmicf%RkillLitrfal2HumOMHeter,        &
    RMaintDefcitLitrfal2HumOMHeter => nmicf%RMaintDefcitLitrfal2HumOMHeter, &
    RN2FixHeter                    => nmicf%RN2FixHeter,                    &
    RKillOMHeter                   => nmicf%RKillOMHeter,                   &
    RkillRecycOMHeter              => nmicf%RkillRecycOMHeter,              &
    RMaintDefcitKillOMHeter        => nmicf%RMaintDefcitKillOMHeter,        &
    RMaintDefcitRecycOMHeter       => nmicf%RMaintDefcitRecycOMHeter,       &
    RNH4TransfLitrHeter            => nmicf%RNH4TransfLitrHeter,            &
    RNO3TransfLitrHeter            => nmicf%RNO3TransfLitrHeter,            &
    RH2PO4TransfLitrHeter          => nmicf%RH2PO4TransfLitrHeter,          &
    RH1PO4TransfSoilHeter          => nmicf%RH1PO4TransfSoilHeter,          &
    RH1PO4TransfBandHeter          => nmicf%RH1PO4TransfBandHeter,          &
    RH1PO4TransfLitrHeter          => nmicf%RH1PO4TransfLitrHeter,          &
    RNH4TransfSoilHeter            => nmicf%RNH4TransfSoilHeter,            &
    k_POM                          => micpar%k_POM,                         &
    k_humus                        => micpar%k_humus,                       &
    mBiomeHeter                    => micstt%mBiomeHeter,                   &
    SolidOM                        => micstt%SolidOM,                       &
    SOMHumProtein                  => micstt%SOMHumProtein,                 &
    SOMHumCarbohyd                 => micstt%SOMHumCarbohyd,                &
    ElmAllocmatMicrblitr2POM       => micfor%ElmAllocmatMicrblitr2POM,      &
    ElmAllocmatMicrblitr2POMU      => micfor%ElmAllocmatMicrblitr2POMU,     &
    icarbhyro                      => micpar%icarbhyro,                     &
    iprotein                       => micpar%iprotein,                      &
    Litrm                          => micfor%litrm,                         &
    NetCAssimhr                    => micflx%NetCAssimhr,                   &
    GrosAssimhr                    => micflx%GrosAssimhr,                   &
    NetNH4Mineralize               => micflx%NetNH4Mineralize,              &
    NetPO4Mineralize               => micflx%NetPO4Mineralize               &
  )
!
!     OMC,OMN,OMP=microbial C,N,P
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     RXOMC,RXOMN,RKillOMHeter=microbial C,N,P decomposition
!     RMaintDefcitKillOMHeter,RXMMN,RXMMP=microbial C,N,P loss from senescence
!

  D550: DO K=1,jcplx
    IF(.not.litrm .OR. (K.NE.k_POM .AND. K.NE.k_humus))THEN
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          D540: DO M=1,2
            MID=micpar%get_micb_id(M,NGL)     
            DO NE=1,NumPlantChemElms     
              mBiomeHeter(NE,MID,K)=mBiomeHeter(NE,MID,K)+NonstX2stBiomHeter(NE,M,NGL,K) &
                -RKillOMHeter(NE,M,NGL,K)-RMaintDefcitKillOMHeter(NE,M,NGL,K)
            ENDDO  
!
!     HUMIFICATION PRODUCTS
!
!     ElmAllocmatMicrblitr2POM=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RkillLitrfal2HumOMHeter=transfer of microbial C,N,P LitrFall to humus
!     RMaintDefcitLitrfal2HumOMHeter,RHMMN,RMaintDefcitLitrfal2HumOMHeter=transfer of senesence LitrFall C,N,P to humus
!
            IF(.not.litrm)THEN
!add as protein
              DO NE=1,NumPlantChemElms
                SolidOM(NE,iprotein,k_humus)=SolidOM(NE,iprotein,k_humus) &
                  +ElmAllocmatMicrblitr2POM(iprotein)*(RkillLitrfal2HumOMHeter(NE,M,NGL,K)&
                  +RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K))
  !add as carbon hydro
                SolidOM(NE,icarbhyro,k_humus)=SolidOM(NE,icarbhyro,k_humus) &
                  +ElmAllocmatMicrblitr2POM(icarbhyro)*(RkillLitrfal2HumOMHeter(NE,M,NGL,K)&
                  +RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K))
              ENDDO
            ELSE

              DO NE=1,NumPlantChemElms
                SOMHumProtein(NE)=SOMHumProtein(NE)+ElmAllocmatMicrblitr2POMU(iprotein) &
                  *(RkillLitrfal2HumOMHeter(NE,M,NGL,K)+RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K))
                SOMHumCarbohyd(NE)=SOMHumCarbohyd(NE)+ElmAllocmatMicrblitr2POMU(icarbhyro) &
                  *(RkillLitrfal2HumOMHeter(NE,M,NGL,K)+RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K))
              ENDDO
            ENDIF
          ENDDO D540

!
!     INPUTS TO NONSTRUCTURAL POOLS
!
!     DOMuptk4GrothHeter=total DOC+acetate uptake
!     RespGrossHeter=total respiration
!     RNOxReduxRespDenitLim=respiration for denitrifcation
!     Resp4NFixHeter=respiration for N2 fixation
!     RCO2ProdHeter=total CO2 emission
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     R3OMC,R3OMN,RkillRecycOMHeter=microbial C,N,P recycling
!     RMaintDefcitRecycOMHeter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!     DOMuptk4GrothHeter,DOMuptk4GrothHeter=DON, DOP uptake
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     RH2PO4TransfSoilHeter,RH2PO4TransfBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     RH1PO4TransfSoilHeter,RH1PO4TransfBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
!     RNH4TransfLitrHeter,RNO3TransfLitrHeter =substrate-limited NH4,NO3 mineraln-immobiln
!     RH2PO4TransfLitrHeter,RH1PO4TransfLitrHeter=substrate-limited H2PO4,HPO4 mineraln-immobiln
!
          CGROMC               = DOMuptk4GrothHeter(ielmc,NGL,K)-RespGrossHeter(NGL,K)-RNOxReduxRespDenitLim(NGL,K)-Resp4NFixHeter(NGL,K)
          RCO2ProdHeter(NGL,K) = RCO2ProdHeter(NGL,K)+Resp4NFixHeter(NGL,K)
          NetCAssimhr          = NetCAssimhr+CGROMC
          GrosAssimhr = GrosAssimhr+DOMuptk4GrothHeter(ielmc,NGL,K)
          MID3        = micpar%get_micb_id(3,NGL)
          D555: DO M       = 1, 2
            DO NE=1,NumPlantChemElms
              mBiomeHeter(NE,MID3,K)=mBiomeHeter(NE,MID3,K)-NonstX2stBiomHeter(NE,M,NGL,K)+RkillRecycOMHeter(NE,M,NGL,K)
            ENDDO
            DO NE=2,NumPlantChemElms
              mBiomeHeter(NE,MID3,K)=mBiomeHeter(NE,MID3,K)+RMaintDefcitRecycOMHeter(NE,M,NGL,K)
            ENDDO
            RCO2ProdHeter(NGL,K)=RCO2ProdHeter(NGL,K)+RMaintDefcitRecycOMHeter(ielmc,M,NGL,K)
          ENDDO D555
          mBiomeHeter(ielmc,MID3,K)=mBiomeHeter(ielmc,MID3,K)+CGROMC
          mBiomeHeter(ielmn,MID3,K)=mBiomeHeter(ielmn,MID3,K)+DOMuptk4GrothHeter(ielmn,NGL,K) &
            +RNH4TransfSoilHeter(NGL,K)+RNH4TransfBandHeter(NGL,K)+RNO3TransfSoilHeter(NGL,K) &
            +RNO3TransfBandHeter(NGL,K)+RN2FixHeter(NGL,K)
          mBiomeHeter(ielmp,MID3,K)=mBiomeHeter(ielmp,MID3,K)+DOMuptk4GrothHeter(ielmp,NGL,K) &
            +RH2PO4TransfSoilHeter(NGL,K)+RH2PO4TransfBandHeter(NGL,K)+RH1PO4TransfSoilHeter(NGL,K) &
            +RH1PO4TransfBandHeter(NGL,K)

          !fix negative microbial N  by immobilization
          if(mBiomeHeter(ielmn,MID3,K)<0._r8)then
            RNH4TransfSoilHeter(NGL,K) = RNH4TransfSoilHeter(NGL,K)-mBiomeHeter(ielmn,MID3,K)
            NetNH4Mineralize           = NetNH4Mineralize-mBiomeHeter(ielmn,MID3,K)
            mBiomeHeter(ielmn,MID3,K)  = 0._r8
          endif            
          !fix negative P biomass by immobilization
          if(mBiomeHeter(ielmp,MID3,K)<0._r8)then
            RH2PO4TransfSoilHeter(NGL,K) = RH2PO4TransfSoilHeter(NGL,K)-mBiomeHeter(ielmp,MID3,K)
            NetPO4Mineralize             = NetPO4Mineralize-mBiomeHeter(ielmp,MID3,K)
            mBiomeHeter(ielmp,MID3,K)    = 0._r8
          endif  
          IF(litrm)THEN
            mBiomeHeter(ielmn,MID3,K)=mBiomeHeter(ielmn,MID3,K)+RNH4TransfLitrHeter(NGL,K)+RNO3TransfLitrHeter(NGL,K)
            mBiomeHeter(ielmp,MID3,K)=mBiomeHeter(ielmp,MID3,K)+RH2PO4TransfLitrHeter(NGL,K)+RH1PO4TransfLitrHeter(NGL,K)
          ENDIF
        enddo
      ENDDO
    ENDIF
  ENDDO D550
  end associate
  end subroutine HeterotrophAnabolicUpdate
!------------------------------------------------------------------------------------------

  subroutine MicrobialLitterColonization(I,J,KL,micfor,micstt,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  integer  :: K,M
  real(r8) :: DOSAK

!     begin_execution
  associate(                                             &
    ROQC4HeterMicActCmpK => ncplxf%ROQC4HeterMicActCmpK, &
    SolidOMCK            => ncplxs%SolidOMCK,            &
    SolidOMActK          => ncplxs%SolidOMActK,          &
    ZEROS                => micfor%ZEROS,                &
    SolidOMAct           => micstt%SolidOMAct,           &
    SolidOM              => micstt%SolidOM,              &
    DOSA                 => micpar%DOSA                  &
  )
!     SolidOMCK,SolidOMActK,OSCX=total,colonized,uncolonized SOC
!     OSA,OSC=colonized,total litter
!     DOSA=rate constant for litter colonization
!     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
!
  D475: DO K=1,KL
    SolidOMCK(K)   = 0.0_r8
    SolidOMActK(K) = 0.0_r8
    DO  M=1,jsken
      SolidOMCK(K)    = SolidOMCK(K)+SolidOM(ielmc,M,K)
      SolidOMActK(K) = SolidOMActK(K)+SolidOMAct(M,K)
    enddo
  ENDDO D475

  D480: DO K=1,KL
    IF(SolidOMCK(K).GT.ZEROS)THEN
      DOSAK=DOSA(K)*AZMAX1(ROQC4HeterMicActCmpK(K))

      D485: DO M=1,jsken
        SolidOMAct(M,K)=AMIN1(SolidOM(ielmc,M,K),SolidOMAct(M,K)+DOSAK*SolidOM(ielmc,M,K)/SolidOMCK(K))
      ENDDO D485
    ELSE
      D490: DO M=1,jsken
        SolidOMAct(M,K)=AMIN1(SolidOM(ielmc,M,K),SolidOMAct(M,K))
      ENDDO D490
    ENDIF
  ENDDO D480
  end associate
  end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

  subroutine AggregateTransfOMBioResdue(micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type),intent(in) :: nmicdiag
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(micfluxtype), intent(inout) :: micflx
  integer  :: K,M,N,NGL,NE
!     begin_execution
  associate(                                                             &
    DOMuptk4GrothHeter           => nmicf%DOMuptk4GrothHeter,            &
    RAnabolDOCUptkHeter          => nmicf%RAnabolDOCUptkHeter,           &
    RAnabolAcetUptkHeter         => nmicf%RAnabolAcetUptkHeter,          &
    RO2UptkHeter                 => nmicf%RO2UptkHeter,                  &
    RO2DmndHeter                 => nmicf%RO2DmndHeter,                  &
    RNO3ReduxHeterSoil           => nmicf%RNO3ReduxHeterSoil,            &
    RNO3ReduxHeterBand           => nmicf%RNO3ReduxHeterBand,            &
    RNO2ReduxHeterSoil           => nmicf%RNO2ReduxHeterSoil,            &
    RNO2ReduxHeterBand           => nmicf%RNO2ReduxHeterBand,            &
    RN2OReduxHeter               => nmicf%RN2OReduxHeter,                &
    RNOxReduxRespDenitLim        => nmicf%RNOxReduxRespDenitLim,         &
    RNH4TransfSoilHeter          => nmicf%RNH4TransfSoilHeter,           &
    RNO3TransfSoilHeter          => nmicf%RNO3TransfSoilHeter,           &
    RH2PO4TransfSoilHeter        => nmicf%RH2PO4TransfSoilHeter,         &
    RNH4TransfBandHeter          => nmicf%RNH4TransfBandHeter,           &
    RNO3TransfBandHeter          => nmicf%RNO3TransfBandHeter,           &
    RH2PO4TransfBandHeter        => nmicf%RH2PO4TransfBandHeter,         &
    RH2ProdHeter                 => nmicf%RH2ProdHeter,                  &
    RNH4TransfLitrHeter          => nmicf%RNH4TransfLitrHeter,           &
    RNO3TransfLitrHeter          => nmicf%RNO3TransfLitrHeter,           &
    RH2PO4TransfLitrHeter        => nmicf%RH2PO4TransfLitrHeter,         &
    RCO2ProdHeter                => nmicf%RCO2ProdHeter,                 &
    RAcettProdHeter              => nmicf%RAcettProdHeter,               &
    RCH4ProdHeter                => nmicf%RCH4ProdHeter,                 &
    RSOxidSoilAutor              => nmicf%RSOxidSoilAutor,               &
    RSOxidBandAutor              => nmicf%RSOxidBandAutor,               &
    RH1PO4TransfSoilHeter        => nmicf%RH1PO4TransfSoilHeter,         &
    RH1PO4TransfBandHeter        => nmicf%RH1PO4TransfBandHeter,         &
    RH1PO4TransfLitrHeter        => nmicf%RH1PO4TransfLitrHeter,         &
    RN2FixHeter                  => nmicf%RN2FixHeter,                   &
    RDcmpProdDOM                 => ncplxf%RDcmpProdDOM,                 &
    RHydlysBioResduOM            => ncplxf%RHydlysBioResduOM,            &
    RHydlysSorptOM               => ncplxf%RHydlysSorptOM,               &
    RDOMSorp                     => ncplxf%RDOMSorp,                     &
    TSensGrowth                  => nmicdiag%TSensGrowth,                &
    RH2UptkAutor                 => nmicdiag%RH2UptkAutor,               &
    RNO2ReduxSoilChemo           => naqfdiag%RNO2ReduxSoilChemo,         &
    RNO2ReduxBandChemo           => naqfdiag%RNO2ReduxBandChemo,         &
    RN2OProdSoilChemo            => naqfdiag%RN2OProdSoilChemo,          &
    RN2OProdBandChemo            => naqfdiag%RN2OProdBandChemo,          &
    RNO3ProdSoilChemo            => naqfdiag%RNO3ProdSoilChemo,          &
    RNO3ProdBandChemo            => naqfdiag%RNO3ProdBandChemo,          &
    VOLWZ                        => nmicdiag%VOLWZ,                      &
    DOMuptk4GrothAutor           => nmicf%DOMuptk4GrothAutor,            &
    RNO2ReduxAutorBand           => nmicf%RNO2ReduxAutorBand,            &
    RNO3UptkAutor                => nmicf%RNO3UptkAutor,                 &
    RCH4ProdAutor                => nmicf%RCH4ProdAutor,                 &
    RNO2ReduxAutorSoil           => nmicf%RNO2ReduxAutorSoil,            &
    RO2UptkAutor                 => nmicf%RO2UptkAutor,                  &
    RNOxReduxRespAutorLim        => nmicf%RNOxReduxRespAutorLim,         &
    RCO2ProdAutor                => nmicf%RCO2ProdAutor,                 &
    RH1PO4TransfLitrAutor        => nmicf%RH1PO4TransfLitrAutor,         &
    RH2PO4TransfLitrAutor        => nmicf%RH2PO4TransfLitrAutor,         &
    RNO3TransfLitrAutor          => nmicf%RNO3TransfLitrAutor,           &
    RNH4TransfLitrAutor          => nmicf%RNH4TransfLitrAutor,           &
    RN2FixAutor                  => nmicf%RN2FixAutor,                   &
    RH1PO4TransfBandAutor        => nmicf%RH1PO4TransfBandAutor,         &
    RH2PO4TransfBandAutor        => nmicf%RH2PO4TransfBandAutor,         &
    RNO3TransfBandAutor          => nmicf%RNO3TransfBandAutor,           &
    RNH4TransfBandAutor          => nmicf%RNH4TransfBandAutor,           &
    RH2PO4TransfSoilAutor        => nmicf%RH2PO4TransfSoilAutor,         &
    RNH4TransfSoilAutor          => nmicf%RNH4TransfSoilAutor,           &
    RNO3TransfSoilAutor          => nmicf%RNO3TransfSoilAutor,           &
    RH1PO4TransfSoilAutor        => nmicf%RH1PO4TransfSoilAutor,         &
    TSens4MicbGrwoth             => micstt%TSens4MicbGrwoth,             &
    VWatMicrobAct                => micstt%VWatMicrobAct,                &
    litrm                        => micfor%litrm,                        &
    Lsurf                        => micfor%Lsurf,                        &
    k_POM                        => micpar%k_POM,                        &
    k_humus                      => micpar%k_humus,                      &
    mid_AmmoniaOxidBacter        => micpar%mid_AmmoniaOxidBacter,        &
    mid_AerobicMethanotrofBacter => micpar%mid_AerobicMethanotrofBacter, &
    mid_NitriteOxidBacter        => micpar%mid_NitriteOxidBacter,        &
    is_activeMicrbFungrpAutor    => micpar%is_activeMicrbFungrpAutor,    &
    RCH4UptkAutor                => micflx%RCH4UptkAutor,                &
    RCO2NetUptkMicb              => micflx%RCO2NetUptkMicb,              &
    RH2NetUptkMicb               => micflx%RH2NetUptkMicb,               &
    RN2NetUptkMicb               => micflx%RN2NetUptkMicb,               &
    RN2ONetUptkMicb              => micflx%RN2ONetUptkMicb,              &
    RO2UptkMicb                  => micflx%RO2UptkMicb,                  &
    RH1PO4MicbTransfBand         => micflx%RH1PO4MicbTransfBand,         &
    RH1PO4MicbTransfSoil         => micflx%RH1PO4MicbTransfSoil,         &
    RH2PO4MicbTransfBand         => micflx%RH2PO4MicbTransfBand,         &
    RH2PO4MicbTransfSoil         => micflx%RH2PO4MicbTransfSoil,         &
    MicrbN2Fix                   => micflx%MicrbN2Fix,                   &
    RNH4MicbTransfBand           => micflx%RNH4MicbTransfBand,           &
    RNH4MicbTransfSoil           => micflx%RNH4MicbTransfSoil,           &
    RNO2MicbTransfBand           => micflx%RNO2MicbTransfBand,           &
    RNO2MicbTransfSoil           => micflx%RNO2MicbTransfSoil,           &
    RNO3MicbTransfBand           => micflx%RNO3MicbTransfBand,           &
    RNO3MicbTransfSoil           => micflx%RNO3MicbTransfSoil,           &
    REcoDOMProd                  => micflx%REcoDOMProd                   &
  )
  D650: DO K=1,jcplx
    IF(.not.litrm.OR.(K.NE.k_POM .AND. K.NE.k_humus))THEN
      DO N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          naqfdiag%tRNH4MicrbTransfSoil   = naqfdiag%tRNH4MicrbTransfSoil+RNH4TransfSoilHeter(NGL,K)
          naqfdiag%tRNO3MicrbTransfSoil   = naqfdiag%tRNO3MicrbTransfSoil+RNO3TransfSoilHeter(NGL,K)
          naqfdiag%tRH2PO4MicrbTransfSoil = naqfdiag%tRH2PO4MicrbTransfSoil+RH2PO4TransfSoilHeter(NGL,K)
          naqfdiag%tRH1PO4MicrbTransfSoil = naqfdiag%tRH1PO4MicrbTransfSoil+RH1PO4TransfSoilHeter(NGL,K)  !> 0 uptake 
          naqfdiag%tRNH4MicrbTransfBand   = naqfdiag%tRNH4MicrbTransfBand+RNH4TransfBandHeter(NGL,K)
          naqfdiag%tRNO3MicrbTransfBand   = naqfdiag%tRNO3MicrbTransfBand+RNO3TransfBandHeter(NGL,K)
          naqfdiag%tRH2PO4MicrbTransfBand = naqfdiag%tRH2PO4MicrbTransfBand+RH2PO4TransfBandHeter(NGL,K)
          naqfdiag%tRH1PO4MicrbTransfBand = naqfdiag%tRH1PO4MicrbTransfBand+RH1PO4TransfBandHeter(NGL,K)
          naqfdiag%TFixN2                 = naqfdiag%TFixN2+RN2FixHeter(NGL,K)
          IF(Lsurf)THEN
            naqfdiag%tRNH4MicrbTransfSoil   = naqfdiag%tRNH4MicrbTransfSoil+RNH4TransfLitrHeter(NGL,K)
            naqfdiag%tRNO3MicrbTransfSoil   = naqfdiag%tRNO3MicrbTransfSoil+RNO3TransfLitrHeter(NGL,K)
            naqfdiag%tRH2PO4MicrbTransfSoil = naqfdiag%tRH2PO4MicrbTransfSoil+RH2PO4TransfLitrHeter(NGL,K)
            naqfdiag%tRH1PO4MicrbTransfSoil = naqfdiag%tRH1PO4MicrbTransfSoil+RH1PO4TransfLitrHeter(NGL,K)
          ENDIF
          naqfdiag%tRCO2MicrbProd    = naqfdiag%tRCO2MicrbProd+RCO2ProdHeter(NGL,K)
          naqfdiag%tRCH4MicrbProd    = naqfdiag%tRCH4MicrbProd+RCH4ProdHeter(NGL,K)
          naqfdiag%tRNOxMicrbRedux   = naqfdiag%tRNOxMicrbRedux+RNOxReduxRespDenitLim(NGL,K)
          naqfdiag%tRO2MicrbUptk     = naqfdiag%tRO2MicrbUptk+RO2UptkHeter(NGL,K)
          naqfdiag%TReduxNO3Soil     = naqfdiag%TReduxNO3Soil+RNO3ReduxHeterSoil(NGL,K)
          naqfdiag%TReduxNO3Band     = naqfdiag%TReduxNO3Band+RNO3ReduxHeterBand(NGL,K)
          naqfdiag%TDeniReduxNO2Soil = naqfdiag%TDeniReduxNO2Soil+RNO2ReduxHeterSoil(NGL,K)
          naqfdiag%TDeniReduxNO2Band = naqfdiag%TDeniReduxNO2Band+RNO2ReduxHeterBand(NGL,K)
          naqfdiag%TReduxNO2Soil     = naqfdiag%TReduxNO2Soil+RNO2ReduxHeterSoil(NGL,K)
          naqfdiag%TReduxNO2Band     = naqfdiag%TReduxNO2Band+RNO2ReduxHeterBand(NGL,K)
          naqfdiag%TReduxN2O         = naqfdiag%TReduxN2O+RN2OReduxHeter(NGL,K)
          naqfdiag%TProdH2           = naqfdiag%TProdH2+RH2ProdHeter(NGL,K)
          naqfdiag%tRO2UptkHeterG    = naqfdiag%tRO2UptkHeterG+RO2UptkHeter(NGL,K)
          naqfdiag%tRO2DmndHeterG    = naqfdiag%tRO2DmndHeterG + RO2DmndHeter(NGL,K)
          nmicf%RO2UptkHeterG(NGL)   = nmicf%RO2UptkHeterG(NGL)+RO2UptkHeter(NGL,K)
          nmicf%RO2DmndHeterG(NGL)   = nmicf%RO2DmndHeterG(NGL)+RO2DmndHeter(NGL,K)
          micflx%TRDOE2DIE(ielmc)    = micflx%TRDOE2DIE(ielmc)+RCO2ProdHeter(NGL,K)+RCH4ProdHeter(NGL,K)
        ENDDO
      ENDDO
    ENDIF
  ENDDO D650

  DO  N=1,NumMicbFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%tRNH4MicrbTransfSoil   = naqfdiag%tRNH4MicrbTransfSoil+RNH4TransfSoilAutor(NGL)
        naqfdiag%tRNO3MicrbTransfSoil   = naqfdiag%tRNO3MicrbTransfSoil+RNO3TransfSoilAutor(NGL)
        naqfdiag%tRH2PO4MicrbTransfSoil = naqfdiag%tRH2PO4MicrbTransfSoil+RH2PO4TransfSoilAutor(NGL)
        naqfdiag%tRH1PO4MicrbTransfSoil = naqfdiag%tRH1PO4MicrbTransfSoil+RH1PO4TransfSoilAutor(NGL)
        naqfdiag%tRNH4MicrbTransfBand   = naqfdiag%tRNH4MicrbTransfBand+RNH4TransfBandAutor(NGL)
        naqfdiag%tRNO3MicrbTransfBand   = naqfdiag%tRNO3MicrbTransfBand+RNO3TransfBandAutor(NGL)
        naqfdiag%tRH2PO4MicrbTransfBand = naqfdiag%tRH2PO4MicrbTransfBand+RH2PO4TransfBandAutor(NGL)
        naqfdiag%tRH1PO4MicrbTransfBand = naqfdiag%tRH1PO4MicrbTransfBand+RH1PO4TransfBandAutor(NGL)
        naqfdiag%TFixN2                 = naqfdiag%TFixN2+RN2FixAutor(NGL)
        IF(Lsurf)THEN
          naqfdiag%tRNH4MicrbTransfSoil   = naqfdiag%tRNH4MicrbTransfSoil+RNH4TransfLitrAutor(NGL)
          naqfdiag%tRNO3MicrbTransfSoil   = naqfdiag%tRNO3MicrbTransfSoil+RNO3TransfLitrAutor(NGL)
          naqfdiag%tRH2PO4MicrbTransfSoil = naqfdiag%tRH2PO4MicrbTransfSoil+RH2PO4TransfLitrAutor(NGL)
          naqfdiag%tRH1PO4MicrbTransfSoil = naqfdiag%tRH1PO4MicrbTransfSoil+RH1PO4TransfLitrAutor(NGL)
        ENDIF
        naqfdiag%tRCO2MicrbProd   = naqfdiag%tRCO2MicrbProd+RCO2ProdAutor(NGL)
        naqfdiag%tRCH4MicrbProd   = naqfdiag%tRCH4MicrbProd+RCH4ProdAutor(NGL)
        naqfdiag%tRNOxMicrbRedux  = naqfdiag%tRNOxMicrbRedux+RNOxReduxRespAutorLim(NGL)
        naqfdiag%tRO2MicrbUptk    = naqfdiag%tRO2MicrbUptk+RO2UptkAutor(NGL)
        naqfdiag%TReduxNO3Soil    = naqfdiag%TReduxNO3Soil+RNO3UptkAutor(NGL)
        naqfdiag%TNitReduxNO2Soil = naqfdiag%TNitReduxNO2Soil+RNO2ReduxAutorSoil(NGL)
        naqfdiag%TNitReduxNO2Band = naqfdiag%TNitReduxNO2Band+RNO2ReduxAutorBand(NGL)
        naqfdiag%TReduxNO2Soil    = naqfdiag%TReduxNO2Soil+RNO2ReduxAutorSoil(NGL)
        naqfdiag%TReduxNO2Band    = naqfdiag%TReduxNO2Band+RNO2ReduxAutorBand(NGL)
        if(micpar%is_CO2_autotroph(N))then
          micflx%TRDOE2DIE(ielmc)=micflx%TRDOE2DIE(ielmc)+RCO2ProdAutor(NGL)
        endif
      ENDDO
    ENDIF
  ENDDO

!     tRCO2GrothAutor=total CO2 uptake by autotrophs, ammonia oxidizer
!  nitrite oxidizer, and hydrogenotophic methanogens
  D645: DO N=1,NumMicbFunGrupsPerCmplx
    IF(micpar%is_CO2_autotroph(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%tRCO2GrothAutor=naqfdiag%tRCO2GrothAutor+DOMuptk4GrothAutor(ielmc,NGL)
      ENDDO
    ENDIF
  ENDDO D645

!
!     ALLOCATE AGGREGATED TRANSFOMBioResduATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2NetUptkMicb=net CO2 uptake
!     tRCO2MicrbProd total CO2 emission by heterotrophs reducing O2
!     tRNOxMicrbRedux=total CO2 emission by denitrifiers reducing NOx
!     RSOxidSoilAutor(3)=CH4 oxidation
!     RCH4UptkAutor=net CH4 uptake, >0, means uptake
!     DOMuptk4GrothHeter=total CH4 uptake by autotrophs
!     tRCH4MicrbProd=total CH4 emission
!     RH2NetUptkMicb=net H2 uptake
!     RH2UptkAutor,TProdH2=total H2 uptake, emission
!     RO2UptkMicb,tRO2MicrbUptk=total O2 uptake
!     RN2NetUptkMicb=total N2 production
!     TReduxN2O=total N2O reduction
!     RN2ONetUptkMicb=total N2O uptake
!     NO2(-) -> N2O
!     TReduxNO2Soil,TReduxNO2Band=total NO2 reduction in non-band,band
!     RN2OProdSoilChemo,RN2OProdBandChemo=nitrous acid reduction in non-band,band
!
!  print*,'RCO2NetUptkMicb',naqfdiag%tRCO2GrothAutor,naqfdiag%tRCO2MicrbProd,naqfdiag%tRNOxMicrbRedux
  RCO2NetUptkMicb = naqfdiag%tRCO2GrothAutor-naqfdiag%tRCO2MicrbProd-naqfdiag%tRNOxMicrbRedux
  RCH4UptkAutor   = -naqfdiag%tRCH4MicrbProd

  DO NGL=JGniA(mid_AerobicMethanotrofBacter),JGnfA(mid_AerobicMethanotrofBacter)
    RCO2NetUptkMicb         = RCO2NetUptkMicb-RSOxidSoilAutor(NGL)
    RCH4UptkAutor           = RCH4UptkAutor+RSOxidSoilAutor(NGL)+DOMuptk4GrothAutor(ielmc,NGL)
    micflx%TRDOE2DIE(ielmc) = micflx%TRDOE2DIE(ielmc)-RSOxidSoilAutor(NGL)-DOMuptk4GrothAutor(ielmc,NGL)
  ENDDO
  !>0. microbial uptake
  RH2NetUptkMicb  = RH2UptkAutor-naqfdiag%TProdH2
  RO2UptkMicb     = naqfdiag%tRO2MicrbUptk
  RN2NetUptkMicb  = -naqfdiag%TReduxN2O
  RN2ONetUptkMicb = -naqfdiag%TReduxNO2Soil-naqfdiag%TReduxNO2Band-RN2OProdSoilChemo &
    -RN2OProdBandChemo+naqfdiag%TReduxN2O
!
  D655: DO K=1,jcplx
    D660: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RDcmpProdDOM(NE,M,K)
      ENDDO
    ENDDO D660

    D665: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D665
    DO NE=1,NumPlantChemElms
      REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RHydlysSorptOM(NE,K)
    ENDDO
    REcoDOMProd(idom_acetate,K)=REcoDOMProd(idom_acetate,K)+RHydlysSorptOM(idom_acetate,K)
    D670: DO N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        REcoDOMProd(idom_doc,K)     = REcoDOMProd(idom_doc,K)-RAnabolDOCUptkHeter(NGL,K)
        REcoDOMProd(idom_don,K)     = REcoDOMProd(idom_don,K)-DOMuptk4GrothHeter(ielmn,NGL,K)
        REcoDOMProd(idom_dop,K)     = REcoDOMProd(idom_dop,K)-DOMuptk4GrothHeter(ielmp,NGL,K)
        REcoDOMProd(idom_acetate,K) = REcoDOMProd(idom_acetate,K)-RAnabolAcetUptkHeter(NGL,K)+RAcettProdHeter(NGL,K)
      ENDDO
    ENDDO D670
    DO NE=idom_beg,idom_end
      REcoDOMProd(NE,K)=REcoDOMProd(NE,K)-RDOMSorp(NE,K)
    ENDDO
  ENDDO D655
!
!     RNH4MicbTransfSoil,RNH4MicbTransfBand=net change in NH4 in band,non-band
!     tRNH4MicrbTransfSoil,tRNH4MicrbTransfBand=total NH4 mineraln-immobn in non-band,band
!     RSOxidSoilAutor(1),RSOxidBandAutor(1)=total NH4 oxidation in non-band,band
!     RNO3MicbTransfSoil,RNO3MicbTransfBand=net change in NO3 in band,non-band
!     tRNO3MicrbTransfSoil,tRNO3MicrbTransfBand=total NO3 immobn in non-band,band
!     RSOxidSoilAutor(2),RSOxidBandAutor(2)=total NO2 oxidation in non-band,band
!     TReduxNO3Soil,TReduxNO3Band=total NO3 reduction in non-band,band
!     RNO3ProdSoilChemo,RNO3ProdBandChemo=NO3 production from nitrous acid reduction in non-band,band
!     RNO2MicbTransfSoil,RNO2MicbTransfBand=net change in NO3 in band,non-band
!     TReduxNO2Soil,TReduxNO2Band=total NO2 reduction in non-band,band
!     RNO2ReduxSoilChemo,RNO2ReduxBandChemo=substrate-limited nitrous acid reduction in non-band,band
!     RH2PO4MicbTransfSoil,RH2PO4MicbTransfBand=net change in H2PO4 in band,non-band
!     tRH2PO4MicrbTransfSoil,tRH2PO4MicrbTransfBand=total H2PO4 mineraln-immobn in non-band,band
!     RH1PO4MicbTransfSoil,RH1PO4MicbTransfBand=net change in HPO4 in band,non-band
!     tRH1PO4MicrbTransfSoil,tRH1PO4MicrbTransfBand=total HPO4 mineraln-immobn in non-band,band
!     MicrbN2Fix=total N2 fixation
!     XZHYS=total H+ production
!     TFixN2=total N2 fixation
!

  micflx%TRDOE2DIE(ielmn)=micflx%TRDOE2DIE(ielmn)+naqfdiag%tRNH4MicrbTransfSoil &
    +naqfdiag%tRNH4MicrbTransfBand+naqfdiag%tRNO3MicrbTransfSoil &
    +naqfdiag%tRNO3MicrbTransfBand+naqfdiag%TFixN2
  micflx%TRDOE2DIE(ielmp)=micflx%TRDOE2DIE(ielmp)+naqfdiag%tRH1PO4MicrbTransfSoil &
    +naqfdiag%tRH2PO4MicrbTransfSoil+naqfdiag%tRH1PO4MicrbTransfBand &
    +naqfdiag%tRH2PO4MicrbTransfBand


  RNH4MicbTransfSoil=-naqfdiag%tRNH4MicrbTransfSoil
  RNO3MicbTransfSoil=-naqfdiag%tRNO3MicrbTransfSoil-naqfdiag%TReduxNO3Soil+RNO3ProdSoilChemo
  RNO2MicbTransfSoil=+naqfdiag%TReduxNO3Soil-naqfdiag%TReduxNO2Soil-RNO2ReduxSoilChemo
  RH2PO4MicbTransfSoil=-naqfdiag%tRH2PO4MicrbTransfSoil
  RH1PO4MicbTransfSoil=-naqfdiag%tRH1PO4MicrbTransfSoil     !< 0 uptake
  RNH4MicbTransfBand=-naqfdiag%tRNH4MicrbTransfBand
  RNO3MicbTransfBand=-naqfdiag%tRNO3MicrbTransfBand-naqfdiag%TReduxNO3Band+RNO3ProdBandChemo
  RNO2MicbTransfBand=naqfdiag%TReduxNO3Band-naqfdiag%TReduxNO2Band-RNO2ReduxBandChemo

  !mid_AmmoniaOxidBacter=1, mid_NitriteOxidBacter=2, mid_AerobicMethanotrofBacter=3
  DO NGL=JGniA(mid_AmmoniaOxidBacter),JGnfA(mid_AmmoniaOxidBacter)
    RNH4MicbTransfSoil=RNH4MicbTransfSoil-RSOxidSoilAutor(NGL)
    RNO2MicbTransfSoil=RNO2MicbTransfSoil+RSOxidSoilAutor(NGL)
    RNH4MicbTransfBand=RNH4MicbTransfBand-RSOxidBandAutor(NGL)

  ENDDO
  DO NGL=JGniA(mid_NitriteOxidBacter),JGnfA(mid_NitriteOxidBacter)
    RNO3MicbTransfSoil=RNO3MicbTransfSoil+RSOxidSoilAutor(NGL)
    RNO2MicbTransfSoil=RNO2MicbTransfSoil-RSOxidSoilAutor(NGL)
    RNO3MicbTransfBand=RNO3MicbTransfBand+RSOxidBandAutor(NGL)
    RNO2MicbTransfBand=RNO2MicbTransfBand-RSOxidBandAutor(NGL)
  ENDDO

  RH2PO4MicbTransfBand = -naqfdiag%tRH2PO4MicrbTransfBand
  RH1PO4MicbTransfBand = -naqfdiag%tRH1PO4MicrbTransfBand
  MicrbN2Fix           = naqfdiag%TFixN2
  TSens4MicbGrwoth     = TSensGrowth
  VWatMicrobAct        = VOLWZ

  DO NGL=JGNiA(mid_AerobicMethanotrofBacter),jGnfA(mid_AerobicMethanotrofBacter)
    naqfdiag%tCH4OxiAero=naqfdiag%tCH4OxiAero+RSOxidSoilAutor(NGL)
  ENDDO

  naqfdiag%tRNH3Oxi=nmicf%RTotNH3OxidSoilAutor+nmicf%RTotNH3OxidBandAutor
  end associate
  end subroutine AggregateTransfOMBioResdue
!------------------------------------------------------------------------------------------

  subroutine SubstrateAttenf4Compet(NGL,N,K,FOXYX, &
    FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA, &
    micfor,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(out):: FOXYX   !fraction of demand over all demand
  real(r8), intent(out):: FNH4X
  real(r8),intent(out) :: FNB3X,FNB4X,FNO3X
  real(r8),intent(out) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8),intent(out) :: FOQC,FOQA
  type(micforctype), intent(in) :: micfor
  type(Cumlate_Flux_Diag_type),INTENT(INOUT)::  naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                                               &
    FracOMActHeter        => nmics%FracOMActHeter,         &
    FracHeterBiomOfActK   => nmics%FracHeterBiomOfActK,    &
    AttenfNH4Heter        => nmicf%AttenfNH4Heter,         &
    AttenfNO3Heter        => nmicf%AttenfNO3Heter,         &
    AttenfH1PO4Heter      => nmicf%AttenfH1PO4Heter,       &
    AttenfH2PO4Heter      => nmicf%AttenfH2PO4Heter,       &
    RO2DmndHetert         => micflx%RO2DmndHetert,         &
    RNH4DmndSoilHeter     => micflx%RNH4DmndSoilHeter,     &
    RNH4DmndBandHeter     => micflx%RNH4DmndBandHeter,     &
    RNO3DmndSoilHeter     => micflx%RNO3DmndSoilHeter,     &
    RNO3DmndBandHeter     => micflx%RNO3DmndBandHeter,     &
    RH2PO4DmndSoilHeter   => micflx%RH2PO4DmndSoilHeter,   &
    RH2PO4DmndBandHeter   => micflx%RH2PO4DmndBandHeter,   &
    RH1PO4DmndSoilHeter   => micflx%RH1PO4DmndSoilHeter,   &
    RH1PO4DmndBandHeter   => micflx%RH1PO4DmndBandHeter,   &
    RDOCUptkHeter         => micflx%RDOCUptkHeter,         &
    RAcetateUptkHeter     => micflx%RAcetateUptkHeter,     &
    RNH4DmndLitrHeter     => micflx%RNH4DmndLitrHeter,     &
    RNO3DmndLitrHeter     => micflx%RNO3DmndLitrHeter,     &
    RH2PO4DmndLitrHeter   => micflx%RH2PO4DmndLitrHeter,   &
    RH1PO4DmndLitrHeter   => micflx%RH1PO4DmndLitrHeter,   &
    litrm                 => micfor%litrm,                 &
    VLNH4                 => micfor%VLNH4,                 &
    VLNHB                 => micfor%VLNHB,                 &
    VLNOB                 => micfor%VLNOB,                 &
    VLNO3                 => micfor%VLNO3,                 &
    VLPOB                 => micfor%VLPOB,                 &
    VLPO4                 => micfor%VLPO4,                 &
    ZEROS                 => micfor%ZEROS,                 &
    RO2EcoDmndPrev        => micfor%RO2EcoDmndPrev,        &
    RNH4EcoDmndSoilPrev   => micfor%RNH4EcoDmndSoilPrev,   &
    RNH4EcoDmndBandPrev   => micfor%RNH4EcoDmndBandPrev,   &
    RNO3EcoDmndSoilPrev   => micfor%RNO3EcoDmndSoilPrev,   &
    RNH4EcoDmndLitrPrev   => micfor%RNH4EcoDmndLitrPrev,   &
    RNO3EcoDmndLitrPrev   => micfor%RNO3EcoDmndLitrPrev,   &
    RH1PO4EcoDmndLitrPrev => micfor%RH1PO4EcoDmndLitrPrev, &
    RH2PO4EcoDmndLitrPrev => micfor%RH2PO4EcoDmndLitrPrev, &
    RNO3EcoDmndBandPrev   => micfor%RNO3EcoDmndBandPrev,   &
    RH2PO4EcoDmndSoilPrev => micfor%RH2PO4EcoDmndSoilPrev, &
    RH2PO4EcoDmndBandPrev => micfor%RH2PO4EcoDmndBandPrev, &
    RH1PO4EcoDmndSoilPrev => micfor%RH1PO4EcoDmndSoilPrev, &
    RH1PO4EcoDmndBandPrev => micfor%RH1PO4EcoDmndBandPrev, &
    RDOMEcoDmndPrev       => micfor%RDOMEcoDmndPrev,       &
    RAcetateEcoDmndPrev   => micfor%RAcetateEcoDmndPrev,   &
    Lsurf                 => micfor%Lsurf,                 &
    SoilMicPMassLayer0    => micfor%SoilMicPMassLayer0     &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(RO2EcoDmndPrev.GT.ZEROS)THEN
    FOXYX=AMAX1(FMN,RO2DmndHetert(NGL,K)/RO2EcoDmndPrev)
  ELSE
    FOXYX=AMAX1(FMN,FracOMActHeter(NGL,K))
  ENDIF
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4DmndSoilHeter(NGL,K)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4DmndBandHeter(NGL,K)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNHB)
  ENDIF
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3DmndSoilHeter(NGL,K)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3DmndBandHeter(NGL,K)/RNO3EcoDmndBandPrev)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
  ENDIF
  IF(RH2PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4DmndSoilHeter(NGL,K)/RH2PO4EcoDmndSoilPrev)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4DmndBandHeter(NGL,K)/RH2PO4EcoDmndBandPrev)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4DmndSoilHeter(NGL,K)/RH1PO4EcoDmndSoilPrev)
  ELSE
    FP14X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4DmndBandHeter(NGL,K)/RH1PO4EcoDmndBandPrev)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
  ENDIF

  IF(RDOMEcoDmndPrev(K).GT.ZEROS)THEN
    FOQC=AMAX1(FMN,RDOCUptkHeter(NGL,K)/RDOMEcoDmndPrev(K))
  ELSE
    FOQC=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF

  naqfdiag%TFOQC=naqfdiag%TFOQC+FOQC
  IF(RAcetateEcoDmndPrev(K).GT.ZEROS)THEN
    FOQA=AMAX1(FMN,RAcetateUptkHeter(NGL,K)/RAcetateEcoDmndPrev(K))
  ELSE
    FOQA=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF

  naqfdiag%TFOQA  = naqfdiag%TFOQA+FOQA
  naqfdiag%TFOXYX = naqfdiag%TFOXYX+FOXYX
  naqfdiag%TFNH4X = naqfdiag%TFNH4X+FNH4X
  naqfdiag%TFNO3X = naqfdiag%TFNO3X+FNO3X
  naqfdiag%TFPO4X = naqfdiag%TFPO4X+FPO4X
  naqfdiag%TFP14X = naqfdiag%TFP14X+FP14X
  naqfdiag%TFNH4B = naqfdiag%TFNH4B+FNB4X
  naqfdiag%TFNO3B = naqfdiag%TFNO3B+FNB3X
  naqfdiag%TFPO4B = naqfdiag%TFPO4B+FPOBX
  naqfdiag%TFP14B = naqfdiag%TFP14B+FP1BX
!
! FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
! MICROBIAL POPULATIONS IN SURFACE RESIDUE
! F*=fraction of substrate uptake relative to total uptake from
! previous hour in surface litter, labels as for soil layers above
!
  IF(litrm)THEN
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNH4Heter(NGL,K)=AMAX1(FMN,RNH4DmndLitrHeter(NGL,K)/RNH4EcoDmndLitrPrev)
    ELSE
      AttenfNH4Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNO3Heter(NGL,K)=AMAX1(FMN,RNO3DmndLitrHeter(NGL,K)/RNO3EcoDmndLitrPrev)
    ELSE
      AttenfNO3Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH2PO4Heter(NGL,K)=AMAX1(FMN,RH2PO4DmndLitrHeter(NGL,K)/RH2PO4EcoDmndLitrPrev)
    ELSE
      AttenfH2PO4Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH1PO4Heter(NGL,K)=AMAX1(FMN,RH1PO4DmndLitrHeter(NGL,K)/RH1PO4EcoDmndLitrPrev)
    ELSE
      AttenfH1PO4Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
  ENDIF

  IF(Lsurf.AND.K.NE.micpar%k_POM.AND.K.NE.micpar%k_humus .AND. SoilMicPMassLayer0.GT.ZEROS)THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+AttenfNH4Heter(NGL,K)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+AttenfNO3Heter(NGL,K)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+AttenfH2PO4Heter(NGL,K)
    naqfdiag%TFP14X=naqfdiag%TFP14X+AttenfH1PO4Heter(NGL,K)
  ENDIF
  end associate
  end subroutine SubstrateAttenf4Compet
!------------------------------------------------------------------------------------------

  subroutine AcetoMethanogenCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQA,ECHZ, &
    FGOCP,FGOAP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: FOQA
  real(r8), intent(in) :: TSensGrowth
  real(r8), intent(in) :: WatStressMicb
  real(r8), intent(out) :: ECHZ
  REAL(R8), intent(out) :: FGOCP,FGOAP
  reaL(r8), intent(out) :: RGOMP         !substrate-limited potential respiration 
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(OMCplx_State_type),intent(in):: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FSBST
  real(r8) :: GOMX,GOMM
  real(r8) :: RGOGY,RGOGZ
  real(r8) :: RGroMax   !kinetically unlimited acetate uptake

! begin_execution
  associate(                                            &
    FBiomStoiScalarHeter => nmics%FBiomStoiScalarHeter, &
    OMActHeter           => nmics%OMActHeter,           &
    RO2DmndHeter         => nmicf%RO2DmndHeter,         &
    RO2Dmnd4RespHeter    => nmicf%RO2Dmnd4RespHeter,    &
    ROQC4HeterMicrobAct  => nmicf%ROQC4HeterMicrobAct,  &
    CDOM                 => ncplxs%CDOM,                &
    DOM                  => micstt%DOM,                 &
    RO2DmndHetert        => micflx%RO2DmndHetert,       &
    RDOCUptkHeter        => micflx%RDOCUptkHeter,       &
    RAcetateUptkHeter    => micflx%RAcetateUptkHeter,   &
    ZERO                 => micfor%ZERO,                &
    TKS                  => micfor%TKS                  &
  )
  GOMX = RGASC*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI))
  GOMM = GOMX/24.0_r8
  ECHZ = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GC4X+GOMM))/EOMH)))
!
!     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     COQA=DOA concentration
!     OQKAM=Km for acetate uptake,FBiomStoiScalarHeter=N,P limitation
!     VMXM=specific respiration rate
!     WatStressMicb=water stress effect, OMA=active biomass
!     TSensGrowth=temp stress effect, FOQA= acetate limitation
!     RGroMax=substrate-limited respiration of acetate
!     RGroMax=competition-limited respiration of acetate
!     OQA=acetate, FOQA=fraction of biological demand for acetate
!     RGOMP=O2-unlimited respiration of acetate
!     ROXY*=O2 demand, RDOCUptkHeter,ROQCA=DOC, acetate demand
!     ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
!
  FSBST                      = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKAM)
  RGOGY                      = AZMAX1(FBiomStoiScalarHeter(NGL,K)*VMXM*WatStressMicb*OMActHeter(NGL,K))
  RGOGZ                      = RGOGY*FSBST*TSensGrowth
  RGroMax                    = AZMAX1(DOM(idom_acetate,K)*FOQA*ECHZ)
  RGOMP                      = AMIN1(RGroMax,RGOGZ)
  FGOCP                      = 0.0_r8
  FGOAP                      = 1.0_r8
  RO2Dmnd4RespHeter(NGL,K)   = 0.0_r8
  RO2DmndHeter(NGL,K)        = 0.0_r8
  RO2DmndHetert(NGL,K)       = 0.0_r8
  RDOCUptkHeter(NGL,K)       = 0.0_r8
  RAcetateUptkHeter(NGL,K)   = RGOGZ
  ROQC4HeterMicrobAct(NGL,K) = 0.0_r8
  !given CH3COOH -> CH4+CO2, 0.5 is into CH4.
  naqfdiag%tCH4ProdAceto=naqfdiag%tCH4ProdAceto+0.5_r8*RGOMP
  end associate
  end subroutine AcetoMethanogenCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobicHeterotrophCatabolism(I,J,NGL,N,K,TSensGrowth,WatStressMicb,FOQC,FOQA, &
    ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC,FOQA
  real(r8), intent(in) :: WatStressMicb  !moisture sensivity of microbial activity
  real(r8), intent(in) :: TSensGrowth    !temperature sensitivity of microbial activity
  real(r8), intent(out) :: ECHZ
  real(r8), intent(out) :: RGOMP  !total DOC/acetate C uptake for potential respiraiton
  REAL(R8), INTENT(OUT) :: FGOCP,FGOAP
  real(r8), intent(out) :: RGOCP  !oxygen-unlimited DOC-based respiraiton
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_State_type),intent(inout):: ncplxs
  type(micfluxtype), intent(inout) :: micflx

  real(r8) :: EO2Q
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX
  real(r8) :: RGOAP
  real(r8) :: FSBST
!     begin_execution
  associate(                                                 &
    OMActHeter             => nmics%OMActHeter,              &
    FBiomStoiScalarHeter   => nmics%FBiomStoiScalarHeter,    &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,       &
    RO2DmndHeter           => nmicf%RO2DmndHeter,            &
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,     &
    ZEROS                  => micfor%ZEROS,                  &
    DOM                    => micstt%DOM,                    &
    mid_Aerob_HeteroBacter => micpar%mid_Aerob_HeteroBacter, &
    mid_Facult_DenitBacter => micpar%mid_Facult_DenitBacter, &
    mid_Aerob_Fungi        => micpar%mid_Aerob_Fungi,        &
    mid_aerob_N2Fixer      => micpar%mid_aerob_N2Fixer,      &
    RO2DmndHetert          => micflx%RO2DmndHetert,          &
    RDOCUptkHeter          => micflx%RDOCUptkHeter,          &
    RAcetateUptkHeter      => micflx%RAcetateUptkHeter,      &
    tRGOXP                 => micflx%tRGOXP,                 &
    tRGOZP                 => micflx%tRGOZP,                 &
    FOCA                   => ncplxs%FOCA,                   &
    FOAA                   => ncplxs%FOAA,                   &
    CDOM                   => ncplxs%CDOM                    &
  )
!     ENERGY YIELDS OF O2 REDOX REACTIONS
!     E* = growth respiration efficiency calculated in PARAMETERS
!
  IF(N.EQ.mid_Aerob_HeteroBacter)THEN
    EO2Q=EO2X
  ELSEIF(N.EQ.mid_Facult_DenitBacter)THEN
    EO2Q=EO2D
  ELSEIF(N.EQ.mid_Aerob_Fungi)THEN
    EO2Q=EO2G
  ELSEIF(N.EQ.mid_aerob_N2Fixer)THEN
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
  FSBSTC = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
  FSBSTA = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
  FSBST  = FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
  RGOCY  = AZMAX1(FBiomStoiScalarHeter(NGL,K)*VMXO*WatStressMicb*OMActHeter(NGL,K))*TSensGrowth
  RGOCZ  = RGOCY*FSBSTC*FOCA(K)
  RGOAZ  = RGOCY*FSBSTA*FOAA(K)

  !obtain kinetically unlimited DOM/acetate uptake 
  RGOCX = AZMAX1(DOM(idom_doc,K)*FOQC*EO2Q)
  RGOAX = AZMAX1(DOM(idom_acetate,K)*FOQA*EO2A)
  !obtain the final uptake
  RGOCP  = AMIN1(RGOCX,RGOCZ)
  RGOAP  = AMIN1(RGOAX,RGOAZ)
  RGOMP  = RGOCP+RGOAP
  tRGOXP = tRGOXP+RGOCX+RGOAX
  tRGOZP = tRGOZP+RGOCZ+RGOAZ
  IF(RGOMP.GT.ZEROS)THEN
    FGOCP = RGOCP/RGOMP
    FGOAP = RGOAP/RGOMP
  ELSE
    FGOCP = 1.0_r8
    FGOAP = 0.0_r8
  ENDIF
!
! ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
! BY HETEROTROPHIC AEROBES

! ECHZ=growth respiration yield, averaged over acetate and DOC/glucose
! RO2Dmnd4RespHeter,RO2DmndHeter,RO2DmndHetert=O2 demand from DOC,DOA oxidation
! RDOCUptkHeter,RAcetateUptkHeter=DOC,DOA demand from DOC,DOA oxidation
! ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
! CH2O+O2 -> CO2 + H2O, (32/12.=2.667)
  ECHZ                     = EO2Q*FGOCP+EO2A*FGOAP
  RO2Dmnd4RespHeter(NGL,K) = 2.667_r8*RGOMP
  RO2DmndHeter(NGL,K)      = RO2Dmnd4RespHeter(NGL,K)
  !make a copy for flux limiter 
  RO2DmndHetert(NGL,K)       = RO2DmndHeter(NGL,K)
  RDOCUptkHeter(NGL,K)       = RGOCZ
  RAcetateUptkHeter(NGL,K)   = RGOAZ
  ROQC4HeterMicrobAct(NGL,K) = RGOCY
!
  end associate
  end subroutine AerobicHeterotrophCatabolism
!------------------------------------------------------------------------------------------

  subroutine AnaerobAcetogenCatabolism(NGL,N,K,TSensGrowth,WatStressMicb,FOQC,ECHZ,FGOCP,FGOAP,RGOMP, &
    micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx)
  !
  !fermentation
  !(CH2O)6 +2H2O-> 2CO2 + 2(CH2O)2 + 4H2, mole based
  !(CH2O)6 -> 2CO2 + 2/3 (CH2O)2 + 8/(72)H2, mass based
  !fermenters only take up DOC/glucose
  !it can be fermenters or anaerobic N2 fixers
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC
  real(r8), intent(in) :: WatStressMicb
  real(r8), intent(in) :: TSensGrowth
  real(r8), intent(out) :: ECHZ
  real(r8), intent(out) :: RGOMP   !potential respiration for metabolism
  real(r8), intent(out) :: FGOCP   !fraction of growth contributed by DOC
  real(r8), intent(out) :: FGOAP   !fraction of growth contributed by acetate
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: GH2X,GH2F
  real(r8) :: GOAX,GOAF
  real(r8) :: GHAX
  REAL(R8) :: oxyi
  real(r8) :: RGOFX,RGOFY,RGOFZ
  real(r8) :: FSBST
  real(r8), parameter :: GlucoseC=72._r8
!     begin_execution
  associate(                                            &
    FBiomStoiScalarHeter => nmics%FBiomStoiScalarHeter, &
    OMActHeter           => nmics%OMActHeter,           &
    RO2Dmnd4RespHeter    => nmicf%RO2Dmnd4RespHeter,    &
    RO2DmndHeter         => nmicf%RO2DmndHeter,         &
    ROQC4HeterMicrobAct  => nmicf%ROQC4HeterMicrobAct,  &
    TKS                  => micfor%TKS,                 &
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
  GH2X = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**4)
  GH2F = GH2X/GlucoseC    !
  GOAX = RGASC*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI)**2)
  GOAF = GOAX/GlucoseC
  GHAX = GH2F+GOAF
  IF(N.EQ.mid_fermentor)THEN
    ECHZ=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCHX-GHAX))/EOMF)))
  ELSE
    !denitrifier
    ECHZ=AMAX1(ENFX,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCHX-GHAX))/EOMN)))
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
  FSBST = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)*OXYI
  RGOFY = AZMAX1(FBiomStoiScalarHeter(NGL,K)*VMXF*WatStressMicb*OMActHeter(NGL,K))
  RGOFZ = RGOFY*FSBST*TSensGrowth
  RGOFX = AZMAX1(DOM(idom_doc,K)*FOQC*ECHZ)
  !potential respiration to expense
  RGOMP                      = AMIN1(RGOFX,RGOFZ)
  FGOCP                      = 1.0_r8
  FGOAP                      = 0.0_r8
  RO2Dmnd4RespHeter(NGL,K)   = 0.0_r8
  RO2DmndHeter(NGL,K)        = 0.0_r8
  RO2DmndHetert(NGL,K)       = 0.0_r8
  RDOCUptkHeter(NGL,K)       = RGOFZ
  RAcetateUptkHeter(NGL,K)   = 0.0_r8
  ROQC4HeterMicrobAct(NGL,K) = RGOFY    !DOC-unlimited fermentation rate
  naqfdiag%tCResp4H2Prod     = naqfdiag%tCResp4H2Prod+RGOMP
!
  end associate
  end subroutine AnaerobAcetogenCatabolism
!------------------------------------------------------------------------------------------

  subroutine HeteroDenitrificCatabolism(NGL,N,K,FOQC,RGOCP, &
    VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)

  !Description
  !FACULTATIVE denitrifcation
  !(CH2O)6  + 6O2 -> 6CO2 +6 H2O
  !(CH2O)6 + 12NO3(-) -> 6CO2 + 12NO2(-) + 6H2O, 12*14/(6*32) = 7/8=0.875
  !(CH2O)6 + 12NO2(-) -> 6CO2 + 6N2O + 12OH(-), 12*14/(6*32)=7/8=0.875
  !(CH2O)6 + 12N2O    -> 6CO2 + 12N2 + 6H2O,  24*14/(6*32)=7/4 = 1.75
  !the reduction of NO2 into NO is not considered
  !Ref: The microbial nitrogen-cycling network, Kuypers et al., 2018
  implicit none
  integer, intent(in) :: NGL,N,K
  REAL(R8), INTENT(IN) :: FOQC
  real(r8), intent(in) :: RGOCP
  real(r8), intent(in) :: VOLWZ            !volume of water to support biogeochemistry
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_State_type),intent(inout) :: ncplxs
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
  real(r8) :: RGOM3X,RNOxReduxRespDenitLim3
  real(r8) :: RNO2UptkSoil,RNO2UptkBand,RDN2X,RDN2T,RGOM2X,RNOxReduxRespDenitLim2,RDN2OX,RGOM1X
  real(r8) :: RNOxReduxRespDenitLim1
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
    OxyLimterHeter         => nmics%OxyLimterHeter,          &
    FracOMActHeter         => nmics%FracOMActHeter,          &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,       &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,       &
    RNO3ReduxHeterSoil     => nmicf%RNO3ReduxHeterSoil,      &
    RNO3ReduxHeterBand     => nmicf%RNO3ReduxHeterBand,      &
    RNO2ReduxHeterSoil     => nmicf%RNO2ReduxHeterSoil,      &
    RNO2ReduxHeterBand     => nmicf%RNO2ReduxHeterBand,      &
    RN2OReduxHeter         => nmicf%RN2OReduxHeter,          &
    RNOxReduxRespDenitLim  => nmicf%RNOxReduxRespDenitLim,   &
    RNOxReduxRespDenitUlm  => nmicf%RNOxReduxRespDenitUlm,   &
    BulkSOMC               => ncplxs%BulkSOMC,               &
    RNO2EcoUptkSoilPrev    => micfor%RNO2EcoUptkSoilPrev,    &
    RN2OEcoUptkSoilPrev    => micfor%RN2OEcoUptkSoilPrev,    &
    RNO3EcoDmndSoilPrev    => micfor%RNO3EcoDmndSoilPrev,    &
    VLNO3                  => micfor%VLNO3,                  &
    ZERO                   => micfor%ZERO,                   &
    ZEROS                  => micfor%ZEROS,                  &
    ZEROS2                 => micfor%ZEROS2,                 &
    RNO2EcoUptkBandPrev    => micfor%RNO2EcoUptkBandPrev,    &
    RNO3EcoDmndBandPrev    => micfor%RNO3EcoDmndBandPrev,    &
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
    FracBulkSOMC           => micstt%FracBulkSOMC,           &
    CH2GS                  => micstt%CH2GS,                  &
    DOM                    => micstt%DOM,                    &
    RNO3ReduxDmndSoilHeter => micflx%RNO3ReduxDmndSoilHeter, &
    RNO2DmndReduxSoilHeter => micflx%RNO2DmndReduxSoilHeter, &
    RN2ODmndReduxHeter     => micflx%RN2ODmndReduxHeter,     &
    RNO2DmndReduxBandHeter => micflx%RNO2DmndReduxBandHeter, &
    RNO3ReduxDmndBandHeter => micflx%RNO3ReduxDmndBandHeter  &
  )
!
! FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
! AND ROOT POPULATIONS
!
! FNO3,FNB3=fraction of total biological demand for NO3
!

  FNO3S=VLNO3
  FNO3B=VLNOB
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3=AMAX1(FMN,RNO3ReduxDmndSoilHeter(NGL,K)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3=AMAX1(FMN,RNO3ReduxDmndBandHeter(NGL,K)/RNO3EcoDmndBandPrev)
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
!     RGOM3X,RNOxReduxRespDenitLim3=substrate-unltd,-ltd respn from NO3 reduction
!     RNO3ReduxDmndSoilHeter,RNO3ReduxDmndBandHeter=demand for NO3 reduction in non-band,band
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
  IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOMC(K).GT.ZERO)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ*FracBulkSOMC(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD3S                        = VMXDXS*FVMXDX
  VMXD3B                        = VMXDXB*FVMXDX
  OQCZ3                         = AZMAX1(DOM(idom_doc,K)*FOQC-RGOCP*OxyLimterHeter(NGL,K))
  OQCD3                         = OQCZ3/eQNO3toOxy
  OQCD3S                        = OQCD3*FNO3S
  OQCD3B                        = OQCD3*FNO3B
  ZNO3SX                        = ZNO3S*FNO3
  ZNO3BX                        = ZNO3B*FNB3
  RNO3UptkSoil                  = AZMAX1(AMIN1(ZNO3SX,VMXD3S))
  RNO3UptkBand                  = AZMAX1(AMIN1(ZNO3BX,VMXD3B))
  RNO3ReduxHeterSoil(NGL,K)     = AZMAX1(AMIN1(VMXD3S,OQCD3S,ZNO3SX))
  RNO3ReduxHeterBand(NGL,K)     = AZMAX1(AMIN1(VMXD3B,OQCD3B,ZNO3BX))
  RDNOX                         = RNO3UptkSoil+RNO3UptkBand
  RDNOT                         = RNO3ReduxHeterSoil(NGL,K)+RNO3ReduxHeterBand(NGL,K)
  RGOM3X                        = eQNO3toOxy*RDNOX
  RNOxReduxRespDenitLim3        = eQNO3toOxy*RDNOT
  RNO3ReduxDmndSoilHeter(NGL,K) = VMXD3S
  RNO3ReduxDmndBandHeter(NGL,K) = VMXD3B
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
!     RGOM2X,RNOxReduxRespDenitLim2=substrate-unltd,-ltd respn from NO2 reduction
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
  IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOMC(K).GT.ZERO)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ*FracBulkSOMC(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD2S=VMXDXS*FVMXDX
  VMXD2B=VMXDXB*FVMXDX
  OQCZ2=AZMAX1(OQCZ3-RNOxReduxRespDenitLim3)
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
  RNOxReduxRespDenitLim2=eQNO2toOxy*RDN2T
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
!     RGOM1X,RNOxReduxRespDenitLim1=substrate-unltd,-ltd  respn from N2O reduction
!     RNOxReduxRespDenitUlm,RNOxReduxRespDenitLim=total substrate-unltd,-ltd respn from NOx reduction
!     RN2ODmndReduxHeter=demand for N2O reduction
!
  VMXD1=(VMXD2-RDN2T)*2.0_r8
  VMXDXS=VMXD1*CZ2OS/(CZ2OS+Z1KM)
  IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOMC(K).GT.ZERO)THEN
    FVMXDX=1.0_r8/(1.0_r8+VMXDXS/(VMKI*VOLWZ*FracBulkSOMC(K)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD1S=VMXDXS*FVMXDX
  OQCZ1=AZMAX1(OQCZ2-RNOxReduxRespDenitLim2)
  OQCD1=OQCZ1/eQN2OtoOxy
  Z2OSX=(Z2OS+RDN2T)*FN2O
  RDN2OX=AZMAX1(AMIN1(Z2OSX,VMXD1S))
  RN2OReduxHeter(NGL,K)=AZMAX1(AMIN1(VMXD1S,OQCD1,Z2OSX))
  RGOM1X=eQN2OtoOxy*RDN2OX
  RNOxReduxRespDenitLim1=eQN2OtoOxy*RN2OReduxHeter(NGL,K)
  RNOxReduxRespDenitUlm(NGL,K)=RGOM3X+RGOM2X+RGOM1X
  RNOxReduxRespDenitLim(NGL,K)=RNOxReduxRespDenitLim3+RNOxReduxRespDenitLim2+RNOxReduxRespDenitLim1
  RN2ODmndReduxHeter(NGL,K)=VMXD1S
  end associate
  end subroutine HeteroDenitrificCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB, &
    micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: OXKX
  real(r8), intent(in) :: FOXYX,RGOMP,RVOXP
  real(r8), intent(in) :: RVOXPA
  real(r8), intent(in) :: RVOXPB
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx
  integer  :: M,MX
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,O2AquaDiffusvity1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMX
  real(r8) :: ROXYFX !dissolution flux
  real(r8) :: ROXYLX
  real(r8) :: RRADO,RMPOX
  real(r8) :: ROXDFQ  !gas dissolution 
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X
  real(r8) :: VOLWPM
  real(r8) :: dsignO2
  ! begin_execution
  associate(                                                 &
    OxyLimterHeter         => nmics%OxyLimterHeter,          &
    OMActHeter             => nmics%OMActHeter,              &
    RO2UptkHeter           => nmicf%RO2UptkHeter,            &
    RespGrossHeter         => nmicf%RespGrossHeter,          &
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,       &
    RO2DmndHeter           => nmicf%RO2DmndHeter,            &
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,       &
    RH2ProdHeter           => nmicf%RH2ProdHeter,            &
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,     &
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,           &
    RAcettProdHeter        => nmicf%RAcettProdHeter,         &
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,           &
    RSOxidSoilAutor        => nmicf%RSOxidSoilAutor,         &
    RSOxidBandAutor        => nmicf%RSOxidBandAutor,         &
    O2_irrig_conc          => micfor%O2_irrig_conc,          &
    O2_rain_conc           => micfor%O2_rain_conc,           &
    COXYE                  => micfor%COXYE,                  &
    RO2GasXchangePrev      => micfor%RO2GasXchangePrev,      &
    RO2AquaXchangePrev     => micfor%RO2AquaXchangePrev,     &
    Irrig2LitRSurf_col         => micfor%Irrig2LitRSurf_col,         &
    Rain2LitRSurf          => micfor%Rain2LitRSurf,          &
    litrm                  => micfor%litrm,                  &
    O2AquaDiffusvity       => micfor%O2AquaDiffusvity,       &
    VLSoilPoreMicP         => micfor%VLSoilPoreMicP,         &
    VLSoilMicP             => micfor%VLSoilMicP,             &
    ZERO                   => micfor%ZERO,                   &
    ZEROS                  => micfor%ZEROS,                  &
    VLsoiAirPM             => micfor%VLsoiAirPM,             &
    VLWatMicP              => micfor%VLWatMicP,              &
    VLWatMicPM             => micfor%VLWatMicPM,             &
    THETPM                 => micfor%THETPM,                 &
    DiffusivitySolutEff    => micfor%DiffusivitySolutEff,    &
    FILM                   => micfor%FILM,                   &
    TortMicPM              => micfor%TortMicPM,              &
    OXYG                   => micstt%OXYG,                   &
    OXYS                   => micstt%OXYS,                   &
    COXYS                  => micstt%COXYS,                  &
    O2GSolubility          => micstt%O2GSolubility,          &
    COXYG                  => micstt%COXYG,                  &
    RNO2DmndReduxSoilHeter => micflx%RNO2DmndReduxSoilHeter, &
    RNO2DmndReduxBandHeter => micflx%RNO2DmndReduxBandHeter, &
    RO2UptkSoilM           => micflx%RO2UptkSoilM            &
  )

  IF(RO2DmndHeter(NGL,K).GT.ZEROS .AND. FOXYX.GT.ZERO)THEN
    IF(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX             = RO2DmndHeter(NGL,K)*dts_gas
      ROXYFX            = RO2GasXchangePrev*dts_gas*FOXYX
      O2AquaDiffusvity1 = O2AquaDiffusvity*dts_gas
      IF(.not.litrm)THEN
        OXYG1  = OXYG*FOXYX
        ROXYLX = RO2AquaXchangePrev*dts_gas*FOXYX
      ELSE
        OXYG1  = COXYG*VLsoiAirPM(1)*FOXYX
        ROXYLX = (RO2AquaXchangePrev+Rain2LitRSurf*O2_rain_conc+Irrig2LitRSurf_col*O2_irrig_conc)*dts_gas*FOXYX
      ENDIF
      if(OXYS <= 0._r8 .and. ROXYLX < 0._r8)then 
        ROXYLX=0._r8
      endif        
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
        DIFOX  = TortMicPM(M)*O2AquaDiffusvity1*12.57_r8*BIOS*OMActHeter(NGL,K)*RRADO
        VOLWOX = VLWatMicPM(M)*O2GSolubility
        VOLPOX = VLsoiAirPM(M)
        VOLWPM = VOLWOX+VOLPOX

        D425: DO MX=1,NPT
          OXYG1  = OXYG1+ROXYFX
          OXYS1  = OXYS1+ROXYLX
          COXYS1 = AMIN1(COXYE*O2GSolubility,AZMAX1(safe_adb(OXYS1,(VLWatMicPM(M)*FOXYX))))

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
          OXYG1 = OXYG1-ROXDFQ
          OXYS1 = OXYS1+ROXDFQ
          !accumulate upatke
          RO2UptkHeter(NGL,K) = RO2UptkHeter(NGL,K)+RMPOX
          RO2UptkSoilM(M)     = RO2UptkSoilM(M)+RMPOX
        ENDDO D425
        
      ENDDO D420
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (OxyLimterHeter)
      !
      !     OxyLimterHeter=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter,RNO2DmndReduxBandHeter=NH3,NO2 oxidation in non-band, band
      !
      OxyLimterHeter(NGL,K)=AMIN1(1.0_r8,AZMAX1(RO2UptkHeter(NGL,K)/RO2DmndHeter(NGL,K)))

    ELSE
      RO2UptkHeter(NGL,K)   = RO2DmndHeter(NGL,K)
      OxyLimterHeter(NGL,K) = 1.0_r8
    ENDIF
  ELSE
    RO2UptkHeter(NGL,K)   = 0.0_r8
    OxyLimterHeter(NGL,K) = 1.0_r8
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2ProdHeter,RAcettProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
  !     RO2Uptk4RespHeter=O2-limited O2 uptake
  !     RSOxidSoilAutor,RSOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !

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
  real(r8) :: RINHP
  real(r8) :: RINHX,RINOP,RINOX,RIPOP,RIPOX,RIP1P
  real(r8) :: RIP1X,RINHPR,RINOPR,RIPOPR,RIP1PR
  real(r8) :: ZNH4M,ZNHBM
  real(r8) :: ZNO3M
  real(r8) :: ZNOBM
  integer :: MID3

!     begin_execution
  associate(                                                       &
   GrowthEnvScalHeter              => nmics%GrowthEnvScalHeter,    &
   OMActHeter                      => nmics%OMActHeter,            &
   AttenfNH4Heter                  => nmicf%AttenfNH4Heter,        &
   AttenfNO3Heter                  => nmicf%AttenfNO3Heter,        &
   AttenfH2PO4Heter                => nmicf%AttenfH2PO4Heter,      &
   AttenfH1PO4Heter                => nmicf%AttenfH1PO4Heter,      &
   RNH4TransfSoilHeter             => nmicf%RNH4TransfSoilHeter,   &
   RNO3TransfSoilHeter             => nmicf%RNO3TransfSoilHeter,   &
   RH2PO4TransfSoilHeter           => nmicf%RH2PO4TransfSoilHeter, &
   RNH4TransfBandHeter             => nmicf%RNH4TransfBandHeter,   &
   RNO3TransfBandHeter             => nmicf%RNO3TransfBandHeter,   &
   RH2PO4TransfBandHeter           => nmicf%RH2PO4TransfBandHeter, &
   RNH4TransfLitrHeter             => nmicf%RNH4TransfLitrHeter,   &
   RNO3TransfLitrHeter             => nmicf%RNO3TransfLitrHeter,   &
   RH2PO4TransfLitrHeter           => nmicf%RH2PO4TransfLitrHeter, &
   RH1PO4TransfSoilHeter           => nmicf%RH1PO4TransfSoilHeter, &
   RH1PO4TransfBandHeter           => nmicf%RH1PO4TransfBandHeter, &
   RH1PO4TransfLitrHeter           => nmicf%RH1PO4TransfLitrHeter, &
   rNCOMC                          => micpar%rNCOMC,               &
   rPCOMC                          => micpar%rPCOMC,               &
   litrm                           => micfor%litrm,                &
   VLNH4                           => micfor%VLNH4,                &
   VLNHB                           => micfor%VLNHB,                &
   VLNO3                           => micfor%VLNO3,                &
   VLNOB                           => micfor%VLNOB,                &
   VLPOB                           => micfor%VLPOB,                &
   VLWatMicP                       => micfor%VLWatMicP,            &
   VOLWU                           => micfor%VOLWU,                &
   VLPO4                           => micfor%VLPO4,                &
   ZNH4B                           => micstt%ZNH4B,                &
   ZNH4S                           => micstt%ZNH4S,                &
   ZNO3B                           => micstt%ZNO3B,                &
   ZNO3S                           => micstt%ZNO3S,                &
   ZNH4TU                          => micstt%ZNH4TU,               &
   ZNO3TU                          => micstt%ZNO3TU,               &
   H1P4TU                          => micstt%H1P4TU,               &
   H2P4TU                          => micstt%H2P4TU,               &
   CNH4BU                          => micstt%CNH4BU,               &
   CNH4SU                          => micstt%CNH4SU,               &
   CH2P4U                          => micstt%CH2P4U,               &
   CH2P4BU                         => micstt%CH2P4BU,              &
   CH1P4U                          => micstt%CH1P4U,               &
   CH1P4BU                         => micstt%CH1P4BU,              &
   CH2P4                           => micstt%CH2P4,                &
   CH2P4B                          => micstt%CH2P4B,               &
   CNH4B                           => micstt%CNH4B,                &
   CNH4S                           => micstt%CNH4S,                &
   CH1P4                           => micstt%CH1P4,                &
   CH1P4B                          => micstt%CH1P4B,               &
   H1PO4                           => micstt%H1PO4,                &
   H1POB                           => micstt%H1POB,                &
   H2PO4                           => micstt%H2PO4,                &
   H2POB                           => micstt%H2POB,                &
   CNO3B                           => micstt%CNO3B,                &
   CNO3S                           => micstt%CNO3S,                &
   CNO3SU                          => micstt%CNO3SU,               &
   CNO3BU                          => micstt%CNO3BU,               &
   mBiomeHeter                     => micstt%mBiomeHeter,          &
   RNO3DmndBandHeter               => micflx%RNO3DmndBandHeter,    &
   RNO3DmndSoilHeter               => micflx%RNO3DmndSoilHeter,    &
   RNH4DmndSoilHeter               => micflx%RNH4DmndSoilHeter,    &
   RNH4DmndBandHeter               => micflx%RNH4DmndBandHeter,    &
   RH2PO4DmndSoilHeter             => micflx%RH2PO4DmndSoilHeter,  &
   RH2PO4DmndBandHeter             => micflx%RH2PO4DmndBandHeter,  &
   RH1PO4DmndSoilHeter             => micflx%RH1PO4DmndSoilHeter,  &
   RH1PO4DmndBandHeter             => micflx%RH1PO4DmndBandHeter,  &
   RNH4DmndLitrHeter               => micflx%RNH4DmndLitrHeter,    &
   RNO3DmndLitrHeter               => micflx%RNO3DmndLitrHeter,    &
   RH2PO4DmndLitrHeter             => micflx%RH2PO4DmndLitrHeter,  &
   RH1PO4DmndLitrHeter             => micflx%RH1PO4DmndLitrHeter,  &
   RNiDemand                       => micflx%RNiDemand          ,  &
   RPiDemand                       => micflx%RPiDemand          ,  &
   NetNH4Mineralize                => micflx%NetNH4Mineralize,     &
   NetPO4Mineralize                => micflx%NetPO4Mineralize      &
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
!     GrowthEnvScalHeter=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RNH4DmndSoilHeter,RNH4DmndBandHeter=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RNH4TransfSoilHeter,RNH4TransfBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
! update may be needed, May 17th, 2023, jyt.
  FNH4S=VLNH4
  FNHBS=VLNHB
  MID3=micpar%get_micb_id(3,NGL)
  RINHP=(mBiomeHeter(ielmc,MID3,K)*rNCOMC(3,NGL,K)-mBiomeHeter(ielmn,MID3,K))
  RNiDemand=RNiDemand+RINHP
  !immobilization
  IF(RINHP.GT.0.0_r8)THEN
    CNH4X                      = AZMAX1(CNH4S-Z4MN)
    CNH4Y                      = AZMAX1(CNH4B-Z4MN)
    RINHX                      = AMIN1(RINHP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*Z4MX)
    RNH4DmndSoilHeter(NGL,K)   = FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RNH4DmndBandHeter(NGL,K)   = FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M                      = Z4MN*VOLWU*FNH4S
    ZNHBM                      = Z4MN*VOLWU*FNHBS
    RNH4TransfSoilHeter(NGL,K) = AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4DmndSoilHeter(NGL,K))
    RNH4TransfBandHeter(NGL,K) = AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4DmndBandHeter(NGL,K))
    !mineralization  (<0._r8)
  ELSE
    RNH4DmndSoilHeter(NGL,K)   = 0.0_r8
    RNH4DmndBandHeter(NGL,K)   = 0.0_r8
    RNH4TransfSoilHeter(NGL,K) = RINHP*FNH4S
    RNH4TransfBandHeter(NGL,K) = RINHP*FNHBS
  ENDIF
  NetNH4Mineralize=NetNH4Mineralize+(RNH4TransfSoilHeter(NGL,K)+RNH4TransfBandHeter(NGL,K))
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
!     GrowthEnvScalHeter=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RNO3DmndSoilHeter,RNO3DmndBandHeter=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RNO3TransfSoilHeter,RNO3TransfBandHeter=substrate-limited NO3 immobiln in non-band, band
!     NetNH4Mineralize=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
  FNO3S=VLNO3
  FNO3B=VLNOB
  RINOP=AZMAX1(RINHP-RNH4TransfSoilHeter(NGL,K)-RNH4TransfBandHeter(NGL,K))
  !immobilization
  IF(RINOP.GT.0.0_r8)THEN
    CNO3X                      = AZMAX1(CNO3S-ZOMN)
    CNO3Y                      = AZMAX1(CNO3B-ZOMN)
    RINOX                      = AMIN1(RINOP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*ZOMX)
    RNO3DmndSoilHeter(NGL,K)   = FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RNO3DmndBandHeter(NGL,K)   = FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M                      = ZOMN*VOLWU*FNO3S
    ZNOBM                      = ZOMN*VOLWU*FNO3B
    RNO3TransfSoilHeter(NGL,K) = AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)),RNO3DmndSoilHeter(NGL,K))
    RNO3TransfBandHeter(NGL,K) = AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)),RNO3DmndBandHeter(NGL,K))
    !mineralization, 
  ELSE
    RNO3DmndSoilHeter(NGL,K)   = 0.0_r8
    RNO3DmndBandHeter(NGL,K)   = 0.0_r8
    RNO3TransfSoilHeter(NGL,K) = 0.0_r8
    RNO3TransfBandHeter(NGL,K) = 0.0_r8
  ENDIF
  NetNH4Mineralize=NetNH4Mineralize+(RNO3TransfSoilHeter(NGL,K)+RNO3TransfBandHeter(NGL,K))
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
!     GrowthEnvScalHeter=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RH2PO4DmndSoilHeter,RH2PO4DmndBandHeter=substrate-unlimited H2PO4 mineraln-immobiln
!     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
!     VOLW=water content
!     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
!     RH2PO4TransfSoilHeter,RH2PO4TransfBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
  FH2PS     = VLPO4
  FH2PB     = VLPOB
  MID3      = micpar%get_micb_id(3,NGL)
  RIPOP     = (mBiomeHeter(ielmc,MID3,K)*rPCOMC(3,NGL,K)-mBiomeHeter(ielmp,MID3,K))
  RPiDemand = RPiDemand+RIPOP
  !immobilization
  IF(RIPOP.GT.0.0_r8)THEN
    CH2PX                        = AZMAX1(CH2P4-HPMN)
    CH2PY                        = AZMAX1(CH2P4B-HPMN)
    RIPOX                        = AMIN1(RIPOP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX)
    RH2PO4DmndSoilHeter(NGL,K)   = FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RH2PO4DmndBandHeter(NGL,K)   = FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM                        = HPMN*VLWatMicP*FH2PS
    H2PBM                        = HPMN*VLWatMicP*FH2PB
    RH2PO4TransfSoilHeter(NGL,K) = AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RH2PO4DmndSoilHeter(NGL,K))
    RH2PO4TransfBandHeter(NGL,K) = AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RH2PO4DmndBandHeter(NGL,K))
  !mineralization  
  ELSE
    RH2PO4DmndSoilHeter(NGL,K)   = 0.0_r8
    RH2PO4DmndBandHeter(NGL,K)   = 0.0_r8
    RH2PO4TransfSoilHeter(NGL,K) = RIPOP*FH2PS
    RH2PO4TransfBandHeter(NGL,K) = RIPOP*FH2PB
  ENDIF
  NetPO4Mineralize=NetPO4Mineralize+(RH2PO4TransfSoilHeter(NGL,K)+RH2PO4TransfBandHeter(NGL,K))
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
!     GrowthEnvScalHeter=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RH1PO4DmndSoilHeter,RH1PO4DmndBandHeter=substrate-unlimited HPO4 mineraln-immobiln
!     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
!     VOLW=water content
!     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
!     RH1PO4TransfSoilHeter,RH1PO4TransfBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band uptake (> 0)
!     NetPO4Mineralize=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
  FH1PS=VLPO4
  FH1PB=VLPOB
! why 0.1 here?  
  RIP1P=0.1_r8*AZMAX1(RIPOP-RH2PO4TransfSoilHeter(NGL,K)-RH2PO4TransfBandHeter(NGL,K))
  !immobilization
  IF(RIP1P.GT.0.0_r8)THEN
    CH1PX                        = AZMAX1(CH1P4-HPMN)
    CH1PY                        = AZMAX1(CH1P4B-HPMN)
    RIP1X                        = AMIN1(RIP1P,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX)
    RH1PO4DmndSoilHeter(NGL,K)   = FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RH1PO4DmndBandHeter(NGL,K)   = FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM                        = HPMN*VLWatMicP*FH1PS
    H1PBM                        = HPMN*VLWatMicP*FH1PB
    RH1PO4TransfSoilHeter(NGL,K) = AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RH1PO4DmndSoilHeter(NGL,K))
    RH1PO4TransfBandHeter(NGL,K) = AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RH1PO4DmndBandHeter(NGL,K))
  !mineralization  
  ELSE
    RH1PO4DmndSoilHeter(NGL,K)   = 0.0_r8
    RH1PO4DmndBandHeter(NGL,K)   = 0.0_r8
    RH1PO4TransfSoilHeter(NGL,K) = 0.0_r8
    RH1PO4TransfBandHeter(NGL,K) = 0.0_r8
  ENDIF
  NetPO4Mineralize=NetPO4Mineralize+(RH1PO4TransfSoilHeter(NGL,K)+RH1PO4TransfBandHeter(NGL,K))
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
!     GrowthEnvScalHeter=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RNH4DmndLitrHeter=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M=NH4 not available for uptake
!     AttenfNH4Heter=fractions of biological NH4 demand
!     RNH4TransfLitrHeter=substrate-limited NH4 mineraln-immobiln
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  IF(litrm)THEN
    RINHPR=RINHP-RNH4TransfSoilHeter(NGL,K)-RNO3TransfSoilHeter(NGL,K)
    !immobilization by tap into the top soil layer
    IF(RINHPR.GT.0.0_r8)THEN
      CNH4X                    = AZMAX1(CNH4SU-Z4MN)
      CNH4Y                    = AZMAX1(CNH4BU-Z4MN)
      RNH4DmndLitrHeter(NGL,K) = AMIN1(RINHPR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M                      = Z4MN*VLWatMicP
      RNH4TransfLitrHeter(NGL,K) = AMIN1(AttenfNH4Heter(NGL,K)*AZMAX1((ZNH4TU-ZNH4M)),RNH4DmndLitrHeter(NGL,K))
    ELSE
      RNH4DmndLitrHeter(NGL,K)   = 0.0_r8
      RNH4TransfLitrHeter(NGL,K) = RINHPR
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+RNH4TransfLitrHeter(NGL,K)
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
!     RNO3DmndLitrHeter=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     GrowthEnvScalHeter=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RNO3TransfLitrHeter=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M=NO3 not available for uptake
!     AttenfNO3Heter=fraction of biological NO3 demand
!     RNO3TransfLitrHeter=substrate-limited NO3 immobiln
!     NetNH4Mineralize=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
    RINOPR=AZMAX1(RINHPR-RNH4TransfLitrHeter(NGL,K))
    !immobilization by tapping into the top soil layer
    IF(RINOPR.GT.0.0_r8)THEN
      CNO3X                    = AZMAX1(CNO3SU-ZOMN)
      CNO3Y                    = AZMAX1(CNO3BU-ZOMN)
      RNO3DmndLitrHeter(NGL,K) = AMAX1(RINOPR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M                      = ZOMN*VLWatMicP
      RNO3TransfLitrHeter(NGL,K) = AMIN1(AttenfNO3Heter(NGL,K)*AZMAX1((ZNO3TU-ZNO3M)),RNO3DmndLitrHeter(NGL,K))
    ELSE
      RNO3DmndLitrHeter(NGL,K)   = 0._r8
      RNO3TransfLitrHeter(NGL,K) = 0._r8
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+RNO3TransfLitrHeter(NGL,K)
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
!     RH2PO4DmndLitrHeter=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     GrowthEnvScalHeter=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RH2PO4DmndLitrHeter=substrate-unlimited H2PO4 mineraln-immobiln
!     VOLW=water content
!     H2P4M=H2PO4 not available for uptake
!     AttenfH2PO4Heter=fractions of biological H2PO4 demand
!     RH2PO4TransfLitrHeter=substrate-limited H2PO4 mineraln-immobiln
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    RIPOPR=RIPOP-RH2PO4TransfSoilHeter(NGL,K)
    !immobilization by tapping into top soil layer
    IF(RIPOPR.GT.0.0_r8)THEN
      CH2PX                      = AZMAX1(CH2P4U-HPMN)
      CH2PY                      = AZMAX1(CH2P4BU-HPMN)
      RH2PO4DmndLitrHeter(NGL,K) = AMIN1(RIPOPR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M                        = HPMN*VOLWU
      RH2PO4TransfLitrHeter(NGL,K) = AMIN1(AttenfH2PO4Heter(NGL,K)*AZMAX1((H2P4TU-H2P4M)),RH2PO4DmndLitrHeter(NGL,K))
    ELSE
      RH2PO4DmndLitrHeter(NGL,K)   = 0.0_r8
      RH2PO4TransfLitrHeter(NGL,K) = RIPOPR
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+RH2PO4TransfLitrHeter(NGL,K)
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
!     RH1PO4DmndLitrHeter=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     GrowthEnvScalHeter=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RH1PO4DmndLitrHeter=substrate-unlimited HPO4 mineraln-immobiln
!     VOLW=water content
!     H1P4M=HPO4 not available for uptake
!     AttenfH1PO4Heter=fraction of biological HPO4 demand
!     RH1PO4TransfLitrHeter=substrate-limited HPO4 minereraln-immobiln
!     NetPO4Mineralize=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RIP1PR=0.1_r8*AZMAX1(RIPOPR-RH2PO4TransfLitrHeter(NGL,K))
    !immobilization
    IF(RIP1PR.GT.0.0_r8)THEN
      CH1PX                      = AZMAX1(CH1P4U-HPMN)
      CH1PY                      = AZMAX1(CH1P4BU-HPMN)
      RH1PO4DmndLitrHeter(NGL,K) = AMIN1(RIP1PR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M                        = HPMN*VOLWU
      RH1PO4TransfLitrHeter(NGL,K) = AMIN1(AttenfH1PO4Heter(NGL,K)*AZMAX1((H1P4TU-H1P4M)),RH1PO4DmndLitrHeter(NGL,K))
    ELSE
      RH1PO4DmndLitrHeter(NGL,K)   = 0.0_r8
      RH1PO4TransfLitrHeter(NGL,K) = 0._r8
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+RH1PO4TransfLitrHeter(NGL,K)
  ENDIF
  end associate
  end subroutine BiomassMineralization
!------------------------------------------------------------------------------------------

  subroutine GatherHetertrophRespiration(I,J,NGL,N,K,RMOMK,RGrowthRespHeter,&
    RMaintDefcitcitHeter,RMaintRespHeter,micfor,micstt,nmicf,nmics)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(out) :: RGrowthRespHeter       !growth respiraiton
  real(r8), intent(out) :: RMaintDefcitcitHeter   !deficit for maintenance respiraiton
  real(r8), intent(out) :: RMaintRespHeter
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  integer :: MID3,MID1
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
!     begin_execution
  associate(                                           &
    TempMaintRHeter     => nmics%TempMaintRHeter,      &
    OMN2                => nmics%OMN2,                 &
    Resp4NFixHeter      => nmicf%Resp4NFixHeter,       &
    RespGrossHeter      => nmicf%RespGrossHeter,       &
    RMaintDmndHeter     => nmicf%RMaintDmndHeter,      &
    RNH4TransfSoilHeter => nmicf%RNH4TransfSoilHeter,  &
    RNO3TransfSoilHeter => nmicf%RNO3TransfSoilHeter,  &
    RN2FixHeter         => nmicf%RN2FixHeter,          &
    rNCOMC              => micpar%rNCOMC,              &
    rPCOMC              => micpar%rPCOMC,              &
    pH                  => micfor%pH,                  &
    ZEROS               => micfor%ZEROS,               &
    mid_aerob_N2Fixer   => micpar%mid_aerob_N2Fixer,   &
    mid_Anaerob_N2Fixer => micpar%mid_Anaerob_N2Fixer, &
    CZ2GS               => micstt%CZ2GS,               &
    mBiomeHeter         => micstt%mBiomeHeter          &
  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TempMaintRHeter=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  FPH                      = 1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX                    = RMOM*TempMaintRHeter(NGL,K)*FPH
  MID1                     = micpar%get_micb_id(1,NGL)
  RMaintDmndHeter(1,NGL,K) = mBiomeHeter(ielmn,MID1,K)*RMOMX*RMOMK(1)
  RMaintDmndHeter(2,NGL,K) = OMN2(NGL,K)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMaintRespHeter=total maintenance respiration
!     RGrowthRespHeter=growth respiration
!     RMaintDefcitcitHeter=senescence respiration
!
  RMaintRespHeter      = RMaintDmndHeter(1,NGL,K)+RMaintDmndHeter(2,NGL,K)
  RGrowthRespHeter     = AZMAX1(RespGrossHeter(NGL,K)-RMaintRespHeter)
  RMaintDefcitcitHeter = AZMAX1(RMaintRespHeter-RespGrossHeter(NGL,K))
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
  IF(N.EQ.mid_aerob_N2Fixer .OR. N.EQ.mid_Anaerob_N2Fixer)THEN
    MID3  = micpar%get_micb_id(3,NGL)
    RGN2P = AZMAX1(mBiomeHeter(ielmc,MID3,K)*rNCOMC(3,NGL,K)-mBiomeHeter(ielmn,MID3,K))/EN2F(N)
    IF(RGrowthRespHeter.GT.ZEROS)THEN
      Resp4NFixHeter(NGL,K)=AMIN1(RGrowthRespHeter*RGN2P/(RGrowthRespHeter+RGN2P) &
        *CZ2GS/(CZ2GS+ZFKM),OMGR*mBiomeHeter(ielmc,MID3,K))
      RN2FixHeter(NGL,K)=Resp4NFixHeter(NGL,K)*EN2F(N)  
    ELSE
      Resp4NFixHeter(NGL,K) = 0.0_r8
      RN2FixHeter(NGL,K)    = 0._r8
    ENDIF    
  ENDIF
  end associate
  end subroutine GatherHetertrophRespiration
!------------------------------------------------------------------------------------------

  subroutine GatherHetertrophAnabolicFlux(I,J,NGL,N,K,ECHZ,FGOCP,FGOAP, &
    RGrowthRespHeter,RMaintDefcitcitHeter,RMaintRespHeter,spomk,rmomk,micfor,&
    micstt,nmicf,nmics,ncplxf,ncplxs,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: ECHZ
  real(r8), intent(in) :: FGOCP,FGOAP
  real(r8), intent(in) :: RGrowthRespHeter      !growth respiraiton
  real(r8), intent(in) :: RMaintDefcitcitHeter  !maintenance deficit
  real(r8), intent(in) :: RMaintRespHeter       !respiraiton for maintenance
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx  
  integer :: M,MID1,MID3,MID,NE
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: CXC
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FracMaintDeficit
!     begin_execution
  associate(                                                                    &
    rCNBiomeActHeter                 => nmics%rCNBiomeActHeter,                 &
    GrowthEnvScalHeter               => nmics%GrowthEnvScalHeter,               &
    FCN                              => nmics%FCN,                              &
    FCP                              => nmics%FCP,                              &
    FracHeterBiomOfActK              => nmics%FracHeterBiomOfActK,              &
    DOMuptk4GrothHeter               => nmicf%DOMuptk4GrothHeter,               &
    RAnabolDOCUptkHeter              => nmicf%RAnabolDOCUptkHeter,              &
    RAnabolAcetUptkHeter             => nmicf%RAnabolAcetUptkHeter,             &
    NonstX2stBiomHeter               => nmicf%NonstX2stBiomHeter,               &
    RespGrossHeter                   => nmicf%RespGrossHeter,                   &
    RNOxReduxRespDenitLim            => nmicf%RNOxReduxRespDenitLim,            &
    RMaintDmndHeter                  => nmicf%RMaintDmndHeter,                  &
    RkillLitfalOMHeter               => nmicf%RkillLitfalOMHeter,               &
    RkillLitrfal2HumOMHeter          => nmicf%RkillLitrfal2HumOMHeter,          &
    RkillLitrfal2ResduOMHeter        => nmicf%RkillLitrfal2ResduOMHeter,        &
    RMaintDefcitLitrfal2HumOMHeter   => nmicf%RMaintDefcitLitrfal2HumOMHeter,   &
    RMaintDefcitLitrfal2ResduOMHeter => nmicf%RMaintDefcitLitrfal2ResduOMHeter, &
    RMaintDefcitLitrfalOMHeter       => nmicf%RMaintDefcitLitrfalOMHeter,       &
    RKillOMHeter                     => nmicf%RKillOMHeter,                     &
    RkillRecycOMHeter                => nmicf%RkillRecycOMHeter,                &
    RMaintDefcitKillOMHeter          => nmicf%RMaintDefcitKillOMHeter,          &
    RMaintDefcitRecycOMHeter         => nmicf%RMaintDefcitRecycOMHeter,         &
    Resp4NFixHeter                   => nmicf%Resp4NFixHeter,                   &
    TDOMUptkHeter                    => ncplxf%TDOMUptkHeter,                   &
    rCNDOM                           => ncplxs%rCNDOM,                          &
    rCPDOM                           => ncplxs%rCPDOM,                          &
    rNCOMC                           => micpar%rNCOMC,                          &
    rPCOMC                           => micpar%rPCOMC,                          &
    FL                               => micpar%FL,                              &
    EHUM                             => micstt%EHUM,                            &
    mBiomeHeter                      => micstt%mBiomeHeter,                     &
    DOM                              => micstt%DOM,                             &
    CDOMuptk1                        => micflx%CDOMuptk1,                       &
    CDOMuptk2                        => micflx%CDOMuptk2,                       &
    tROMT                            => micflx%tROMT,                           &
    tGROMO                           => micflx%tGROMO,                          &
    ZEROS                            => micfor%ZEROS,                           &
    ZERO                             => micfor%ZERO                             &
  )

!     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
!     FROM O2, NOX AND C REDUCTION
!
!     CGOMX=DOC+acetate uptake from aerobic growth respiration
!     CGOMD=DOC+acetate uptake from denitrifier growth respiration
!     RMaintRespHeter=maintenance respiration
!     RespGrossHeter=total respiration
!     RNOxReduxRespDenitLim=respiration for denitrifcation
!     Resp4NFixHeter=respiration for N2 fixation
!     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
!     DOMuptk4GrothHeter,RAnabolDOCUptkHeter,RAnabolAcetUptkHeter=total DOC+acetate, DOC, acetate uptake(heterotrophs
!     DOMuptk4GrothHeter=total CO2,CH4 uptake (autotrophs)
!     DOMuptk4GrothHeter,DOMuptk4GrothHeter=DON, DOP uptake
!     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
!     OQN,OPQ=DON,DOP
!     FracHeterBiomOfActK=faction of OMActHeterin total OMA
!     rCNDOM,rCPDOM=DON/DOC, DOP/DOC
!     FCN,FCP=limitation from N,P
!
!  write(*,*)'bf cgom',N,K,ECHZ,ENOX
  CGOMX=AMIN1(RMaintRespHeter,RespGrossHeter(NGL,K))+Resp4NFixHeter(NGL,K)+(RGrowthRespHeter-Resp4NFixHeter(NGL,K))/ECHZ
  CGOMD=RNOxReduxRespDenitLim(NGL,K)/ENOX
  CDOMuptk1 = CDOMuptk1+CGOMX
  CDOMuptk2 = CDOMuptk2+CGOMD
  tROMT     = tROMT+RMaintRespHeter
  tGROMO    = tGROMO+RespGrossHeter(NGL,K)
  DOMuptk4GrothHeter(ielmc,NGL,K)=CGOMX+CGOMD
!  write(*,*)'CGOMX',N,K,CGOMX
  RAnabolDOCUptkHeter(NGL,K)=CGOMX*FGOCP+CGOMD
  RAnabolAcetUptkHeter(NGL,K)=CGOMX*FGOAP
  CGOXC=RAnabolDOCUptkHeter(NGL,K)+RAnabolAcetUptkHeter(NGL,K)
  !obtain organic nutrient uptake
  DOMuptk4GrothHeter(ielmn,NGL,K)=AZMAX1(AMIN1(DOM(idom_don,K)*FracHeterBiomOfActK(NGL,K),CGOXC*rCNDOM(K)/FCN(NGL,K)))
  DOMuptk4GrothHeter(ielmp,NGL,K)=AZMAX1(AMIN1(DOM(idom_dop,K)*FracHeterBiomOfActK(NGL,K),CGOXC*rCPDOM(K)/FCP(NGL,K)))

  TDOMUptkHeter(idom_doc,K)=TDOMUptkHeter(idom_doc,K)+RAnabolDOCUptkHeter(NGL,K)
  TDOMUptkHeter(idom_acetate,K)=TDOMUptkHeter(idom_acetate,K)+RAnabolAcetUptkHeter(NGL,K)
  TDOMUptkHeter(idom_don,K)=TDOMUptkHeter(idom_don,K)+DOMuptk4GrothHeter(ielmn,NGL,K)
  TDOMUptkHeter(idom_dop,K)=TDOMUptkHeter(idom_dop,K)+DOMuptk4GrothHeter(ielmp,NGL,K)

!
!     TRANSFER UPTAKEN C,N,P FROM STORAGE/nonstructural TO ACTIVE BIOMASS
!
!     OMC,OMN,OMP=nonstructural C,N,P
!     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
!     rNCOMC,rPCOMC=maximum microbial N:C, P:C ratios
!     RCCC,RCCN,RCCP=C,N,P recycling fractions
!     RCCZ,RCCY=min, max C recycling fractions
!     RCCX,RCCQ=max N,P recycling fractions
!
  MID1=micpar%get_micb_id(1,NGL);MID3=micpar%get_micb_id(3,NGL)

  IF(mBiomeHeter(ielmc,MID3,K).GT.ZEROS .AND.mBiomeHeter(ielmc,MID1,K).GT.ZEROS)THEN
    CCC=AZMAX1(AMIN1(1.0_r8 &
      ,mBiomeHeter(ielmn,MID3,K)/(mBiomeHeter(ielmn,MID3,K)+mBiomeHeter(ielmc,MID3,K)*rNCOMC(3,NGL,K)) &
      ,mBiomeHeter(ielmp,MID3,K)/(mBiomeHeter(ielmp,MID3,K)+mBiomeHeter(ielmc,MID3,K)*rPCOMC(3,NGL,K))))
    CXC  = mBiomeHeter(ielmc,MID3,K)/mBiomeHeter(ielmc,MID1,K)
    C3C  = 1.0_r8/(1.0_r8+CXC/CKC)
    CNC  = AZMAX1(AMIN1(1.0_r8,mBiomeHeter(ielmc,MID3,K)/(mBiomeHeter(ielmc,MID3,K)+mBiomeHeter(ielmn,MID3,K)/rNCOMC(3,NGL,K))))
    CPC  = AZMAX1(AMIN1(1.0_r8,mBiomeHeter(ielmc,MID3,K)/(mBiomeHeter(ielmc,MID3,K)+mBiomeHeter(ielmp,MID3,K)/rPCOMC(3,NGL,K))))
    RCCC = RCCZ+AMAX1(CCC,C3C)*RCCY
    RCCN = CNC*RCCX
    RCCP = CPC*RCCQ
  ELSE
    RCCC=RCCZ
    RCCN=0.0_r8
    RCCP=0.0_r8
  ENDIF
!
!     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
!
!     CGOMZ=transfer from nonstructural to structural microbial C
!     GrowthEnvScalHeter=temperature+water stress function
!     OMGR=rate constant for transferring nonstructural to structural C
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     FL=partitioning between labile and resistant microbial components
!     OMC,OMN,OMP=nonstructural microbial C,N,P
! M=1:labile, 2, resistant
  CGOMZ=GrowthEnvScalHeter(NGL,K)*OMGR*AZMAX1(mBiomeHeter(ielmc,MID3,K))
  D745: DO M=1,2
    NonstX2stBiomHeter(ielmc,M,NGL,K)=FL(M)*CGOMZ
    IF(mBiomeHeter(ielmc,MID3,K).GT.ZEROS)THEN
      Do NE=2,NumPlantChemElms
        NonstX2stBiomHeter(NE,M,NGL,K)=AMIN1(FL(M)*AZMAX1(mBiomeHeter(NE,MID3,K)) &
          ,NonstX2stBiomHeter(ielmc,M,NGL,K)*mBiomeHeter(NE,MID3,K)/mBiomeHeter(ielmc,MID3,K))
      ENDDO  
    ELSE
      NonstX2stBiomHeter(2:NumPlantChemElms,M,NGL,K)=0.0_r8
    ENDIF
!
!     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
!     RATE, TEMPERATURE
!
!     SPOMX=rate constant for microbial decomposition
!     SPOMC=basal decomposition rate
!     SPOMK=effect of low microbial C concentration on microbial decay
!     RXOMC,RXOMN,RKillOMHeter=microbial C,N,P decomposition
!     RkillLitfalOMHeter,RDOMN,RDOMP=microbial C,N,P LitrFall
!     R3OMC,R3OMN,RkillRecycOMHeter=microbial C,N,P recycling
!
    MID   = micpar%get_micb_id(M,NGL)
    SPOMX = SQRT(GrowthEnvScalHeter(NGL,K))*SPOMC(M)*SPOMK(M)
    DO NE=1,NumPlantChemElms
      RKillOMHeter(NE,M,NGL,K)=AZMAX1(mBiomeHeter(NE,MID,K)*SPOMX)
    ENDDO

    RkillLitfalOMHeter(ielmc,M,NGL,K)=RKillOMHeter(ielmc,M,NGL,K)*(1.0_r8-RCCC)
    RkillLitfalOMHeter(ielmn,M,NGL,K)=RKillOMHeter(ielmn,M,NGL,K)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RkillLitfalOMHeter(ielmp,M,NGL,K)=RKillOMHeter(ielmp,M,NGL,K)*(1.0_r8-RCCC)*(1.0_r8-RCCP)

    DO NE=1,NumPlantChemElms    
      RkillRecycOMHeter(NE,M,NGL,K)=RKillOMHeter(NE,M,NGL,K)-RkillLitfalOMHeter(NE,M,NGL,K)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RkillLitrfal2HumOMHeter=transfer of microbial C,N,P LitrFall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RkillLitrfal2ResduOMHeter,RkillLitrfal2ResduOMHeter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!

      RkillLitrfal2HumOMHeter(NE,M,NGL,K)=AZMAX1(RkillLitfalOMHeter(NE,M,NGL,K)*EHUM)
  !
  !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
  !
      RkillLitrfal2ResduOMHeter(NE,M,NGL,K)=RkillLitfalOMHeter(NE,M,NGL,K)-RkillLitrfal2HumOMHeter(NE,M,NGL,K)
    ENDDO
  ENDDO D745
!
!     MICROBIAL DECOMPOSITION/renomalization WHEN MAINTENANCE RESPIRATION
!     EXCEEDS UPTAKE
!
!     OMC,OMN,OMP=microbial C,N,P
!     RMaintRespHeter=total maintenance respiration
!     RMaintDefcitcitHeter=senescence respiration
!     RCCC=C recycling fraction
!     RMaintDefcitKillOMHeter,RXMMN,RXMMP=microbial C,N,P loss from senescence
!     RMaintDmndHeter=maintenance respiration
!     CNOMA,CPOMA=N:C,P:C ratios of active biomass
!     RDMMC,RMaintDefcitLitrfalOMHeter,RDMMP=microbial C,N,P LitrFall from senescence
!     RMaintDefcitRecycOMHeter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!

  IF(RMaintDefcitcitHeter.GT.ZEROS.AND.RMaintRespHeter.GT.ZEROS.AND.RCCC.GT.ZERO)THEN
    FracMaintDeficit=RMaintDefcitcitHeter/RMaintRespHeter
    D730: DO M=1,2
      MID=micpar%get_micb_id(M,NGL)    
      RMaintDefcitKillOMHeter(ielmc,M,NGL,K)=AMIN1(mBiomeHeter(ielmc,MID,K),AZMAX1(FracMaintDeficit*RMaintDmndHeter(M,NGL,K)/RCCC))
      RMaintDefcitKillOMHeter(ielmn,M,NGL,K)=AMIN1(mBiomeHeter(ielmn,MID,K),&
        AZMAX1(RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*rCNBiomeActHeter(ielmn,NGL,K)))
      RMaintDefcitKillOMHeter(ielmp,M,NGL,K)=AMIN1(mBiomeHeter(ielmp,MID,K),&
        AZMAX1(RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*rCNBiomeActHeter(ielmp,NGL,K)))
        
      RMaintDefcitLitrfalOMHeter(ielmc,M,NGL,K)=RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*(1.0_r8-RCCC)
      RMaintDefcitLitrfalOMHeter(ielmn,M,NGL,K)=RMaintDefcitKillOMHeter(ielmn,M,NGL,K)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
      RMaintDefcitLitrfalOMHeter(ielmp,M,NGL,K)=RMaintDefcitKillOMHeter(ielmp,M,NGL,K)*(1.0_r8-RCCC)*(1.0_r8-RCCP)

      DO NE=1,NumPlantChemElms
        RMaintDefcitRecycOMHeter(NE,M,NGL,K)=RMaintDefcitKillOMHeter(NE,M,NGL,K)-RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RMaintDefcitLitrfal2HumOMHeter,RHMMN,RMaintDefcitLitrfal2HumOMHeter=transfer of senesence LitrFall C,N,P to humus
!     EHUM=humus transfer fraction
!     RMaintDefcitLitrfal2ResduOMHeter,RCMMN,RMaintDefcitLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
!

        RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K)   = AZMAX1(RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)*EHUM)
        RMaintDefcitLitrfal2ResduOMHeter(NE,M,NGL,K) = RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)-RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K)
      ENDDO
    ENDDO D730
  ELSE
    D720: DO M=1,2
      DO NE=1,NumPlantChemElms
        RMaintDefcitKillOMHeter(NE,M,NGL,K)          = 0.0_r8
        RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)       = 0.0_r8
        RMaintDefcitRecycOMHeter(NE,M,NGL,K)         = 0.0_r8
        RMaintDefcitLitrfal2HumOMHeter(NE,M,NGL,K)   = 0.0_r8
        RMaintDefcitLitrfal2ResduOMHeter(NE,M,NGL,K) = 0.0_r8
      ENDDO

    ENDDO D720
  ENDIF
  end associate
  end subroutine GatherHetertrophAnabolicFlux

end module MicBGCMod
