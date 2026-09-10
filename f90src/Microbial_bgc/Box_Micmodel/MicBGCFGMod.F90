module MicBGCMod
!!
! DESCRIPTION:
! codes to do soil biological transfOMBioResduations
!
! USES:
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
  use MicAutoCPLXMod
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod
  use AerobicBacteriaMod,   only: AerobicHeteroBactCatabolism
  use CyanoBacteriaMod,     only: CyanoBacteriaCatabolism
  use DenitrifierMod,       only: HeteroDenitrificCatabolism
  use FermenterMod,         only: AcetogFermentCatabolism
  use FungiMod,             only: AerobicFungiCatabolism
  use MethanogenMod,        only: AcetoMethanogenCatabolism
  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: jcplx,NumMicbHFunGrupsPerCmplx,NumMicbAFunGrupsPerCmplx,jsken,ndbiomcp,nlbiomcp
  integer, pointer :: JGniA(:)
  integer, pointer :: JGnfA(:)
  integer, pointer :: JGniH(:)
  integer, pointer :: JGnfH(:)
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
  NumMicbHFunGrupsPerCmplx  =micpar%NumMicbHFunGrupsPerCmplx
  NumMicbAFunGrupsPerCmplx  =micpar%NumMicbAFunGrupsPerCmplx  
  jsken =micpar%jsken
  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp

  JGniA => micpar%JGniA
  JGnfA => micpar%JGnfA
  JGniH => micpar%JGniH
  JGnfH => micpar%JGnfH

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
  type(Microbe_Diag_type),intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='SoilBGCOneLayer'
!local variables
  integer :: LL,K,KL,NGL
  integer :: M,N

  type(Microbe_State_type) :: nmics
  type(Microbe_Flux_type) :: nmicf
  type(OMCplx_Flux_type) :: ncplxf
  type(OMCplx_State_type) :: ncplxs
  real(r8) :: totOMbeg(1:NumPlantChemElms)
  real(r8) :: totOMend(1:NumPlantChemElms)
  real(r8) :: domsb(1:NumPlantChemElms),somb(1:NumPlantChemElms)
  real(r8) :: sorbomb(1:NumPlantChemElms),biomheterb(1:NumPlantChemElms)
  real(r8) :: biomautob(1:NumPlantChemElms),biomresb(1:NumPlantChemElms)
  real(r8) :: domse(1:NumPlantChemElms),some(1:NumPlantChemElms)
  real(r8) :: sorbome(1:NumPlantChemElms),biomhetere(1:NumPlantChemElms)
  real(r8) :: biomautoe(1:NumPlantChemElms),biomrese(1:NumPlantChemElms)

! begin_execution
  call PrintInfo('beg '//subname)
  call nmicf%Init(jcplx)
  call nmics%Init(jcplx)
  call ncplxf%Init()
  call ncplxs%Init()
  call naqfdiag%ZeroOut()
  call micflx%ZeroOut()
  call nmicdiag%ZeroOut()

  call StageBGCEnvironCondition(I,J,micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)

!  call SumOneLayer('b',KL,micstt,micfor%litrM,totOMbeg,domsb,somb,sorbomb,biomheterb,biomautob,biomresb)
!
  ! write(*,*)'ActiveMicrobeCatabolism'
  call ActiveMicrobeCatabolism(I,J,KL,micfor,micstt,micflx,nmicdiag, &
    naqfdiag,nmicf,nmics,ncplxf,ncplxs)
  !
  !write(*,*)'ChemoDenitrification'
  call ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)

        !
        !write(*,*)'PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C'
  call OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs,nmicdiag)
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
  call SolidOMDecomposition(I,J,KL,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs,micflx)

  call RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)

  call RDOMSorption(KL,micfor,micstt,nmicf,ncplxf,ncplxs)

  call AutotrophAnabolicUpdate(micfor,micstt,nmicf,nmicdiag)

  call HeterotrophAnabolicUpdate(I,J,micfor,micstt,nmicf,micflx)
  
  call MicrobialLitterColonization(I,J,KL,micfor,micstt,ncplxf,ncplxs,nmicdiag)
  
  !     AGGREGATE ALL TRANSFOMBioResduATIONS CALCULATED ABOVE FOR EACH N,K
  !
  call AggregateTransfOMBioResdue(KL,micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)

!  call SumOneLayer('e',KL,micstt,micfor%litrM,totOMend,domse,some,sorbome,biomhetere,biomautoe,biomrese)

  micstt%TotActMicrobiom=nmicdiag%TotActMicrobiom
  !obtain total gross respiration by heterotrophs
  naqfdiag%tRespGrossHeter =sum(nmicf%RespGrossHeter)
  call nmics%destroy()
  call nmicf%destroy()
  call ncplxf%destroy()
  call ncplxs%destroy()
  call PrintInfo('end '//subname)
  end subroutine SoilBGCOneLayer
!------------------------------------------------------------------------------------------
  subroutine SumOneLayer(tag,KL,micstt,litrM,toms,doms,som,sorbom,biomheter,biomauto,biomres)
  !
  !Description:
  !summarize all organic matter states in the specified layer
  implicit none
  character(len=*), intent(in) :: tag
  integer, intent(in) :: KL
  type(micsttype), intent(in) :: micstt  
  logical, intent(in) :: litrM
  real(r8), intent(out) :: toms(1:NumPlantChemElms)
  real(r8), intent(out) :: doms(1:NumPlantChemElms)
  real(r8), intent(out) :: sorbom(1:NumPlantChemElms)
  real(r8), intent(out) :: som(1:NumPlantChemElms)
  real(r8), intent(out) :: biomheter(1:NumPlantChemElms)
  real(r8), intent(out) :: biomauto(1:NumPlantChemElms)
  real(r8), intent(out) :: biomres(1:NumPlantChemElms)
  integer :: K,M,NBM,NE,N,NGL,NLB,MID

!     begin_execution
  associate(                                      &
    SOMPomProtein     => micstt%SOMPomProtein,    &
    SOMHumProtein     => micstt%SOMHumProtein,    &
    SOMHumCarbohyd    => micstt%SOMHumCarbohyd,   &
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
  som=0._r8;biomheter=0._r8;biomauto=0._r8;biomres=0._r8;doms=0._r8;sorbom=0._r8
  if(litrM)then
    DO NE=1,NumPlantChemElms
      toms(NE)=toms(NE)+SOMPomProtein(NE)+SOMHumProtein(NE)+SOMHumCarbohyd(NE)
      som(NE)=som(NE)+SOMPomProtein(NE)+SOMHumProtein(NE)+SOMHumCarbohyd(NE)
    ENDDO
  endif
  
  DO K=1,KL
   !add solid organic matter
    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        TOMS(NE)=TOMS(NE)+SolidOM(NE,M,K)
        som(NE)=som(NE)+SolidOM(NE,M,K)
      ENDDO
    ENDDO

    !add dead microbial biomass
    DO NBM=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        TOMS(NE)    = TOMS(NE)+OMBioResdu(NE,NBM,K)
        biomres(NE) = biomres(NE)+OMBioResdu(NE,NBM,K)
      ENDDO
    ENDDO

    !add sorbed om and DOM
    DO NE=1,NumPlantChemElms
      TOMS(NE)   = TOMS(NE)+SorbedOM(NE,K)
      TOMS(NE)   = TOMS(NE)+DOM(NE,K)
      sorbom(NE) = sorbom(NE)+SorbedOM(NE,K)
      doms(NE)   = doms(NE)+DOM(NE,K)
    ENDDO

    !add live heterotrophic biomass     
    DO N=1,NumMicbHFunGrupsPerCmplx
      if(.not.micpar%is_activeMicrbFungrpHeter(N))cycle
      DO NGL=JGniH(n),JGnfH(n)
        DO nlb=1,micpar%nlbiomcp
          MID=micpar%get_micb_id(nlb,NGL)                
          DO NE=1,NumPlantChemElms
            TOMS(NE)=TOMS(NE)+mBiomeHeter(NE,mid,K)
            biomheter(NE)=biomheter(NE)+mBiomeHeter(NE,mid,K)
          ENDDO
        ENDDO  
      ENDDO
    ENDDO 
    sorbom(ielmc) = sorbom(ielmc)+SorbedOM(idom_acetate,K)
    doms(ielmc)   = doms(ielmc)+DOM(idom_acetate,K)
    TOMS(ielmc)   = TOMS(ielmc)+DOM(idom_acetate,K)+SorbedOM(idom_acetate,K)
  ENDDO

  !add live autotrophic biomass
  DO N=1,NumMicbAFunGrupsPerCmplx
    if(.not.micpar%is_activeMicrbFungrpAutor(N))cycle
    DO NGL=JGniA(N),JGnfA(N)
      DO nlb=1,micpar%nlbiomcp
        MID=micpar%get_micb_id(nlb,NGL)    
        DO NE=1,NumPlantChemElms
          TOMS(NE)=TOMS(NE)+mBiomeAutor(NE,mid)      
          biomauto(NE)=biomauto(NE)+mBiomeAutor(NE,mid)      
        ENDDO 
      ENDDO    
    ENDDO
  ENDDO

  end associate

  end subroutine SumOneLayer
!------------------------------------------------------------------------------------------

  subroutine StageBGCEnvironCondition(I,J,micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  integer, intent(out) :: KL                     !total number of complexes to do calculation
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_State_type),intent(inout):: ncplxs
  character(len=*), parameter :: subname='StageBGCEnvironCondition'
  real(r8) :: TBulkOMC  
  integer  :: K,NE
  integer  :: M,N,NGL,MID1,MID2
  real(r8) :: ORGCL
  real(r8) :: TKSO
  real(r8) :: TSolidOMC,TSorbedOMC
  real(r8) :: theta0=0.5_r8            !tuning parameter for litter layer moist volume 
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
    TOMEK                     => nmicdiag%TOMEK,                   &
    FOCA                      => ncplxs%FOCA,                      &
    FOAA                      => ncplxs%FOAA,                      &
    rCNDOM                    => ncplxs%rCNDOM,                    &
    rCPDOM                    => ncplxs%rCPDOM,                    &
    CDOM                      => ncplxs%CDOM,                      &
    OMBioResduK               => ncplxs%OMBioResduK,               &
    SolidOMCK                 => ncplxs%SolidOMCK,                 &
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
    FBiomNutStoiScalAutor     => nmics%FBiomNutStoiScalAutor,      &
    rNCOMCAutor               => micpar%rNCOMCAutor,               &
    rPCOMCAutor               => micpar%rPCOMCAutor,               &
    rNCOMC                    => micpar%rNCOMC,                    &
    rPCOMC                    => micpar%rPCOMC,                    &
    FL                        => micpar%FL,                        &
    k_humus                   => micpar%k_humus,                   &
    k_POM                     => micpar%k_POM,                     &
    is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, &
    mid_Facult_DenitBacter    => micpar%mid_Facult_DenitBacter,    &
    mid_AutoAmmoniaOxidBacter => micpar%mid_AutoAmmoniaOxidBacter, &
    litrm                     => micfor%litrm,                     &
    VLSoilPoreMicP            => micfor%VLSoilPoreMicP,            &
    VWatLitRHoldCapcity       => micfor%VWatLitRHoldCapcity,       &
    VLWatMicP                 => micfor%VLWatMicP,                 &
    VOLW0                     => micfor%VOLW0,                     & !input: volumetric liquid water in litter layer, [m3 water m-3 pore]
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
    tOMActC                   => micstt%tOMActC,                   &
    SolidOM                   => micstt%SolidOM,                   &
    SolidOMAct                => micstt%SolidOMAct,                &
    TSolidOMActC              => micstt%TSolidOMActC,              &
    TSolidOMC                 => micstt%TSolidOMC,                 &
    OMBioResdu                => micstt%OMBioResdu,                &
    SorbedOM                  => micstt%SorbedOM,                  &
    mBiomeHeter               => micstt%mBiomeHeter,               &
    DOM                       => micstt%DOM,                       &
    H1PO4                     => micstt%H1PO4,                     &
    H1POB                     => micstt%H1POB,                     &
    H2PO4                     => micstt%H2PO4,                     &
    H2POB                     => micstt%H2POB,                     &
    ZNH4B                     => micstt%ZNH4B,                     & !NH4 concentration in banded soil
    ZNH4S                     => micstt%ZNH4S,                     & !NH4 concentration in non-banded soil
    ZNO2B                     => micstt%ZNO2B,                     &
    ZNO2S                     => micstt%ZNO2S,                     &
    ZNO3B                     => micstt%ZNO3B,                     &
    ZNO3S                     => micstt%ZNO3S,                     &
    mBiomeAutor               => micstt%mBiomeAutor,               &
    FracBulkSOMC              => micstt%FracBulkSOMC               &
  )
  call PrintInfo('beg '//subname)
! get KL, the number of mic-om complexes

!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH TempOffset FOR THERMAL ADAPTATION
  IF(litrm)THEN
    ! surface litter layer
    KL=micpar%NumOfLitrCmplxs
    IF(VWatLitRHoldCapcity.GT.ZEROS2)THEN
      ThetaLitr = AMIN1(VOLW0/VWatLitRHoldCapcity,1._r8)
      ThetaZ    = AZMAX1(AMIN1(ThetaLitr,FieldCapacity)-THETY)
!      VOLWZ     = ThetaZ**2/(theta0+ThetaZ)*VWatLitRHoldCapcity    !effective water volume, the definition seems not very accurate
!      VOLWZ     = ThetaZ**0.4_r8*VWatLitRHoldCapcity 
      VOLWZ     = ThetaZ*VWatLitRHoldCapcity 
    ELSE
      VOLWZ=0.0_r8
    ENDIF
  ELSE
!     non-surface layer
    KL     = micpar%jcplx
    ThetaZ = AZMAX1((AMIN1(AMAX1(0.5_r8*POROS,FieldCapacity),THETW)-THETY))
!    VOLWZ  = ThetaZ**2/(theta0+ThetaZ)*VLSoilPoreMicP  !effective water volume
!    VOLWZ  = ThetaZ**0.4_r8*VLSoilPoreMicP  !effective water volume
    VOLWZ  = ThetaZ*VLSoilPoreMicP  !effective water volume
  ENDIF
  VOLWZ=real_truncate(VOLWZ,1.e-5_r8)

!     TKS=soil temperature
!     TempOffset=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TSensGrowth,TSensMaintR=temperature function for growth,maintenance respiration
! the offset could be micobial guild/group specific
  TKSO=real_truncate(TKS+TempOffset,1.e-3_r8)

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
      D895: DO N=1,NumMicbHFunGrupsPerCmplx
        DO NGL=JGniH(n),JGnfH(n)
          MID1=micpar%get_micb_id(ibiom_kinetic,NGL)
          
          IF(mBiomeHeter(ielmc,MID1,K).GT.ZEROS)THEN
            rCNBiomeActHeter(ielmn,NGL,K)=AZMAX1(mBiomeHeter(ielmn,MID1,K)/mBiomeHeter(ielmc,MID1,K))
            rCNBiomeActHeter(ielmp,NGL,K)=AZMAX1(mBiomeHeter(ielmp,MID1,K)/mBiomeHeter(ielmc,MID1,K))
          ELSE
            rCNBiomeActHeter(ielmn,NGL,K)=rNCOMC(ibiom_kinetic,NGL,K)
            rCNBiomeActHeter(ielmp,NGL,K)=rPCOMC(ibiom_kinetic,NGL,K)
          ENDIF
          OMActHeter(NGL,K)           = AZMAX1(mBiomeHeter(ielmc,MID1,K)/FL(ibiom_kinetic))
          FCN(NGL,K)                  = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActHeter(ielmn,NGL,K)/rNCOMC(ibiom_kinetic,NGL,K))))
          FCP(NGL,K)                  = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActHeter(ielmp,NGL,K)/rPCOMC(ibiom_kinetic,NGL,K))))
          FBiomStoiScalarHeter(NGL,K) = AMIN1(FCN(NGL,K),FCP(NGL,K))

!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
          TotActMicrobiom=TotActMicrobiom+OMActHeter(NGL,K)
          IF(N.EQ.mid_Facult_DenitBacter)THEN
            TotBiomNO2Consumers=TotBiomNO2Consumers+OMActHeter(NGL,K)
          ENDIF
          MID2=micpar%get_micb_id(ibiom_struct,NGL)
          OMC2(NGL,K)=AZMAX1(AMIN1(OMActHeter(NGL,K)*FL(ibiom_struct),mBiomeHeter(ielmc,MID2,K)))
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

! the autotrohpic complex
  DO N=1,NumMicbAFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        MID1=micpar%get_micb_id(ibiom_kinetic,NGL)
        IF(mBiomeAutor(ielmc,MID1).GT.ZEROS)THEN
          rCNBiomeActAutor(ielmn,NGL)=AZMAX1(mBiomeAutor(ielmn,MID1)/mBiomeAutor(ielmc,MID1))
          rCNBiomeActAutor(ielmp,NGL)=AZMAX1(mBiomeAutor(ielmp,MID1)/mBiomeAutor(ielmc,MID1))
        ELSE
          rCNBiomeActAutor(ielmn,NGL)=rNCOMCAutor(ibiom_kinetic,NGL)
          rCNBiomeActAutor(ielmp,NGL)=rPCOMCAutor(ibiom_kinetic,NGL)
        ENDIF
        OMActAutor(NGL)           = AZMAX1(mBiomeAutor(ielmc,MID1)/FL(ibiom_kinetic))
        FCNAutor(NGL)             = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActAutor(ielmn,NGL)/rNCOMCAutor(ibiom_kinetic,NGL))))
        FCPAutor(NGL)             = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActAutor(ielmp,NGL)/rPCOMCAutor(ibiom_kinetic,NGL))))
        FBiomNutStoiScalAutor(NGL) = AMIN1(FCNAutor(NGL),FCPAutor(NGL))
      !
      !       TOTAL BIOMASS
      !       OMC2=active biomass in recalcitrant fraction
      !
        TotActMicrobiom=TotActMicrobiom+OMActAutor(NGL)

        IF(N.EQ.mid_AutoAmmoniaOxidBacter)THEN
          TotBiomNO2Consumers=TotBiomNO2Consumers+OMActAutor(NGL)
        ENDIF

        MID2           = micpar%get_micb_id(ibiom_struct,NGL)
        OMC2Autor(NGL) = AZMAX1(AMIN1(OMActAutor(NGL)*FL(ibiom_struct),mBiomeAutor(ielmc,MID2)))
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
    TOMEK(:,K)        = 0.0_r8
    tMaxNActMicrbK(K) = 0.0_r8
    tMaxPActMicrbK(K) = 0.0_r8
    D685: DO N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        if(OMActHeter(NGL,K)>ZEROS)THEN
          TOMEK(ielmc,K)    = TOMEK(ielmc,K)+OMActHeter(NGL,K)
          DO NE=2,NumPlantChemElms
            TOMEK(NE,K)    = TOMEK(NE,K)+OMActHeter(NGL,K)*rCNBiomeActHeter(NE,NGL,K)
          ENDDO
          tMaxNActMicrbK(K) = tMaxNActMicrbK(K)+OMActHeter(NGL,K)*rNCOMC(ibiom_kinetic,NGL,K)   !maximum total N in active micb
          tMaxPActMicrbK(K) = tMaxPActMicrbK(K)+OMActHeter(NGL,K)*rPCOMC(ibiom_kinetic,NGL,K)   !maximum total P in active micb
        ENDIF
      ENDDO
    ENDDO D685
    tOMActC=tOMActC+TOMEK(ielmc,K)
  ENDDO D690
  
  TOMEAutoK(:)      = 0._r8
  DO N=1,NumMicbAFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      if(OMActAutor(NGL)>ZEROS)then
        TOMEAutoK(ielmc)  = TOMEAutoK(ielmc)+OMActAutor(NGL)
        TOMEAutoK(ielmn)  = TOMEAutoK(ielmn)+OMActAutor(NGL)*rCNBiomeActAutor(ielmn,NGL)
        TOMEAutoK(ielmp)  = TOMEAutoK(ielmp)+OMActAutor(NGL)*rCNBiomeActAutor(ielmp,NGL)
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
  call PrintInfo('end '//subname)
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
  real(r8) :: TOMCNK(2)   !total microbial C in labile (M=1) and recalcitrant (M=2) biomass components 
  real(r8) :: COMC

  associate(                             &
    ZEROS       => micfor%ZEROS,         &
    mBiomeHeter => micstt%mBiomeHeter    &
  )

  TOMCNK(:)=0.0_r8
  DO NGL=JGniH(N),JGnfH(N)
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
      COMC     = TOMCNK(M)/ORGCL
      SPOMK(M) = COMC/(COMC+COMKI)
      RMOMK(M) = COMC/(COMC+COMKM)
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
  DO NGL=JGniA(N),JGnfA(N)
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

  subroutine ActiveMicrobeCatabolism(I,J,KL,micfor,micstt,micflx,nmicdiag,naqfdiag, &
    nmicf, nmics,ncplxf,ncplxs)
  !
  !  Description:
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL                !number of complexes involved in calculation
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
  character(len=*), parameter :: subname='ActiveMicrobeCatabolism'
  REAL(R8) :: OXKX
  real(r8) :: ORGCL
! begin_execution
  associate(                                                     &
  TempMaintRHeter           => nmics%TempMaintRHeter,            &
  GrowthEnvScalHeter        => nmics%GrowthEnvScalHeter,         &
  OMActHeter                => nmics%OMActHeter,                 &
  TDOMUptkHeter             => ncplxf%TDOMUptkHeter,             &
  k_humus                   => micpar%k_humus,                   &
  k_POM                     => micpar%k_POM,                     &
  is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, &
  is_activeMicrbFungrpHeter => micpar%is_activeMicrbFungrpHeter, &
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
  call PrintInfo('beg '//subname)
  !Heterotrophic microbes
  TDOMUptkHeter(:,:) = 0.0_r8
  ORGCL              = AMIN1(1.0E+05_r8*SoilMicPMassLayer,ORGC)

  !heterotrophs
  D760: DO K=1,KL

    DO  N=1,NumMicbHFunGrupsPerCmplx
      if(.not.is_activeMicrbFungrpHeter(N))cycle
      call GetMicrobDensFactorHeter(N,K,micfor, micstt, ORGCL,SPOMK,RMOMK)
      
      call ActiveHeterotrophsK(I,J,N,K,SPOMK,RMOMK,&
        micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx,nmicdiag)
    ENDDO
  ENDDO D760

! Autotrophic microbes
  
  DO  N=1,NumMicbAFunGrupsPerCmplx
    IF(.not.is_activeMicrbFungrpAutor(N))cycle

    call GetMicrobDensFactorAutor(N,micfor, micstt, ORGCL,SPOMK,RMOMK)

    call ActiveAutotrophs(I,J,N,SPOMK, RMOMK, &
      micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs,nmicdiag)

  ENDDO

  !summarize microbial activity as a proxy for hydrolysis
  DO  K=1,KL
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        nmicdiag%ROQC4HeterMicActCmpK(K)=nmicdiag%ROQC4HeterMicActCmpK(K)+nmicf%ROQC4HeterMicrobAct(NGL,K)
      enddo
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine ActiveMicrobeCatabolism
!------------------------------------------------------------------------------------------
  subroutine ActiveHeterotrophsK(I,J,N,K,SPOMK,RMOMK,&
     micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: K,N
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
  type(Microbe_Diag_type), intent(inout) :: nmicdiag  
  character(len=*), parameter :: subname='ActiveHeterotrophsK'
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX
  real(r8) :: RMaintDefcitcitHeter
  integer  :: NGL

! begin_execution
  associate(                                             &
    ZNH4T               => nmicdiag%ZNH4T,               &
    ZNO3T               => nmicdiag%ZNO3T,               &
    ZNO2T               => nmicdiag%ZNO2T,               &
    H2P4T               => nmicdiag%H2P4T,               &
    H1P4T               => nmicdiag%H1P4T,               &
    VOLWZ               => nmicdiag%VOLWZ,               &
    TotBiomNO2Consumers => nmicdiag%TotBiomNO2Consumers, &
    OMActHeter          => nmics%OMActHeter,             &
    litrm               => micfor%litrm,                 &
    ZEROS               => micfor%ZEROS,                 &
    VLSoilPoreMicP      => micfor%VLSoilPoreMicP         &
  )
  call PrintInfo('beg '//subname)
  !
  ! HETEROTROPHIC BIOMASS RESPIRATION

  IF(micpar%is_aerobic_hetr(N))THEN
    !   RESPIRATION BY HETEROTROPHIC AEROBES:
    !   N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI, (6)N2 FIXERS

    IF(N.EQ.micpar%mid_Aerob_Fungi)THEN
      !(3)FUNGI
      call AerobicFungiCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
    ELSEIF(N.EQ.micpar%mid_HeterMixtCynoBacter)THEN
      call CyanoBacteriaCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
    else  
      !heterotrophic aerboic bacteria, including facultative denitrifiers
      call AerobicHeteroBactCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
    endif

    IF(N.EQ.micpar%mid_Facult_DenitBacter .AND. (.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS))THEN
      !non-litter layer denitrifcation
      call HeteroDenitrificCatabolism(N,K,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
    ENDIF

  ELSEIF(micpar%is_anaerobic_hetr(N))THEN
    !     RESPIRATION BY HETEROTROPHIC ANAEROBES:
    !     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
    call AcetogFermentCatabolism(N,K,RMOMK,micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx,nmicdiag)

  ELSEIF(N.EQ.micpar%mid_HeterAcetoCH4GenArchea)THEN
    !     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
    call AcetoMethanogenCatabolism(N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  ENDIF

  !
  DO NGL=JGniH(N),JGnfH(N)        
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle
    !     BIOMASS DECOMPOSITION AND MINERALIZATION
    !
    call SubstrateAttenf4Compet(NGL,N,K,FNH4X, FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
      micfor,naqfdiag,nmicf,nmics,micflx)

    call BiomassMineralization(NGL,N,K,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
      ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt, nmicf,nmics,micflx)

    call GatherHeterotrophRespiration(I,J,NGL,N,K,RMaintDefcitcitHeter, &
      micfor,micstt,nmicf,nmics,micflx)

    call GatherHetertrophAnabolicFlux(I,J,NGL,N,K,&
      RMaintDefcitcitHeter,SPOMK,micfor,micstt,nmicf, nmics,ncplxf,ncplxs,micflx)
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine ActiveHeterotrophsK

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
  character(len=*), parameter :: subname='ChemoDenitrification'
!     begin_execution
  associate(                                               &
    TSensGrowth           => nmicdiag%TSensGrowth,         &
    RNO2ReduxSoilChemo    => naqfdiag%RNO2ReduxSoilChemo,  &
    RNO2ReduxBandChemo    => naqfdiag%RNO2ReduxBandChemo,  &
    RN2OProdSoilChemo     => naqfdiag%RN2OProdSoilChemo,   &
    RN2OProdBandChemo     => naqfdiag%RN2OProdBandChemo,   &
    RNO3ProdSoilChemo     => naqfdiag%RNO3ProdSoilChemo,   &
    RNO3ProdBandChemo     => naqfdiag%RNO3ProdBandChemo,   &
    RNO2ReduxChemo        => naqfdiag%RNO2ReduxChemo,      &
    RNO2EcoUptkSoilPrev   => micfor%RNO2EcoUptkSoilPrev,   &
    RNO2EcoUptkBandPrev   => micfor%RNO2EcoUptkBandPrev,   &
    ph                    => micfor%pH,                    &
    VLWatMicPM            => micfor%VLWatMicPM,            &
    ZEROS                 => micfor%ZEROS,                 &
    ZERO                  => micfor%ZERO,                  &
    VLNOB                 => micfor%VLNOB,                 &
    VLNO3                 => micfor%VLNO3,                 &
    CNO2S                 => micstt%CNO2S,                 &
    CNO2B                 => micstt%CNO2B,                 &
    ZNO2B                 => micstt%ZNO2B,                 &
    ZNO2S                 => micstt%ZNO2S,                 &
    RNO2DmndSoilChemoPrev => micflx%RNO2DmndSoilChemoPrev, &
    RNO2DmndBandChemoPrev => micflx%RNO2DmndBandChemoPrev, &
    RNO2DmndSoilChemo     => micflx%RNO2DmndSoilChemo,     &
    RNO2DmndBandChemo     => micflx%RNO2DmndBandChemo      &
  )
  call PrintInfo('beg '//subname)
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
    FNO2=AMAX1(FMN,RNO2DmndSoilChemoPrev/RNO2EcoUptkSoilPrev)
  ELSE
    FNO2=FMN*VLNO3
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2DmndBandChemoPrev/RNO2EcoUptkBandPrev)
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
  call PrintInfo('end '//subname)
  end associate
  end subroutine ChemoDenitrification
!------------------------------------------------------------------------------------------

  subroutine OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs,nmicdiag)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_Flux_type), intent(inout):: ncplxf
  type(OMCplx_State_type),intent(inout):: ncplxs
  type(Microbe_Diag_type), intent(inout) :: nmicdiag  
  integer  :: K,M,N,KK,NGL,MID,idom,NE
  real(r8) :: OSRT
  real(r8) :: XFRK,XFROM(idom_beg:idom_end)
  real(r8) :: XFME
  character(len=*), parameter :: subname='OMTransferForPriming'

!     begin_execution
  associate(                                               &
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,      &
    XferBiomeHeterK      => nmicf%XferBiomeHeterK,         &
    ROQC4HeterMicActCmpK => nmicdiag%ROQC4HeterMicActCmpK, &
    XferRespHeterK       => ncplxf%XferRespHeterK,         &
    XferDOMK             => ncplxf%XferDOMK,               &
    BulkSOMC             => ncplxs%BulkSOMC,               &
    TMicHeterActivity    => micstt%TMicHeterActivity,      &
    DOM                  => micstt%DOM,                    &
    mBiomeHeter          => micstt%mBiomeHeter,            &
    ZEROS                => micfor%ZEROS,                  &
    TScal4Difsvity       => micfor%TScal4Difsvity          &
  )
  call PrintInfo('beg '//subname)
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
            XFROM(idom)=AMAX1(AMIN1(FPRIM*TScal4Difsvity*(DOM(idom,K)*BulkSOMC(KK)-DOM(idom,KK)*BulkSOMC(K))/OSRT,DOM(idom,K)),-DOM(idom,KK))
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
          D850: DO N=1,NumMicbHFunGrupsPerCmplx
            DO  M=1,nlbiomcp
              DO NGL=JGniH(N),JGnfH(N)
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
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        do NGL=JGniH(N),JGnfH(N)
          MID=micpar%get_micb_id(M,NGL)        
          DO NE=1,NumPlantChemElms
            mBiomeHeter(NE,MID,K)=AZERO(mBiomeHeter(NE,MID,K)+XferBiomeHeterK(NE,M,NGL,K))
          ENDDO
        enddo
      enddo
    enddo
  ENDDO D840
  call PrintInfo('end '//subname)
  end associate
  end subroutine OMTransferForPriming
!------------------------------------------------------------------------------------------

  subroutine RDOMSorption(KL,micfor,micstt,nmicf,ncplxf,ncplxs)
  !
  ! Description:
  !DOC ADSORPTION - DESORPTION  
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
    RDOMSorp             => ncplxf%RDOMSorp,            &  !output: >0, adsorption of free DOM/aceate
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
  !TSORP=0._r8
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
        OQEX(idom) = AMAX1(ZEROS,DOM(idom,K)-TDOMUptkHeter(idom,K))  !free DOM
        OHEX(idom) = AMAX1(ZEROS,SorbedOM(idom,K))                   !adsorbed DOM
      ENDDO

      VLSoilPoreMicPX = SoilMicPMassLayer*AECX*HSORP*FracBulkSOMC(K)
      VLSoilPoreMicPW = VLWatMicPM(NPH)*FracBulkSOMC(K)
      IF(FOCA(K).GT.ZERO)THEN
        VOLCX                = FOCA(K)*VLSoilPoreMicPX
        VOLCW                = FOCA(K)*VLSoilPoreMicPW
        RDOMSorp(idom_doc,K) = TSORP*(OQEX(idom_doc)*VOLCX-OHEX(idom_doc)*VOLCW)/(VOLCX+VOLCW)
      ELSE
        RDOMSorp(idom_doc,K)=TSORP*(OQEX(idom_doc)*VLSoilPoreMicPX-OHEX(idom_doc)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      ENDIF

      IF(FOAA(K).GT.ZERO)THEN
        VOLAX                    = FOAA(K)*VLSoilPoreMicPX
        VOLAW                    = FOAA(K)*VLSoilPoreMicPW
        RDOMSorp(idom_acetate,K) = TSORP*(OQEX(idom_acetate)*VOLAX-OHEX(idom_acetate)*VOLAW)/(VOLAX+VOLAW) 
      ELSE
        RDOMSorp(idom_acetate,K)=TSORP*(OQEX(idom_acetate)*VLSoilPoreMicPX-OHEX(idom_acetate)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      ENDIF
      RDOMSorp(idom_don,K)=TSORP*(OQEX(idom_don)*VLSoilPoreMicPX-OHEX(idom_don)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      RDOMSorp(idom_dop,K)=TSORP*(OQEX(idom_dop)*VLSoilPoreMicPX-OHEX(idom_dop)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
    ELSE
      DO idom=idom_beg,idom_end
        RDOMSorp(idom,K)=0.0_r8
      enddo
    ENDIF
    DO idom=idom_beg,idom_end
      DOM(idom,K)      = DOM(idom,K)-RDOMSorp(idom,K)
      SorbedOM(idom,K) = SorbedOM(idom,K)+RDOMSorp(idom,K)
    ENDDO
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
  real(r8) :: RHOSCM,dHyd
  real(r8) :: FCNK(1:jcplx),FCPK(1:jcplx)

  character(len=*), parameter :: subname='SolidOMDecomposition'

!     begin_execution
  associate(                                               &
    tRHydlySOM           => micflx%tRHydlySOM,             &
    RHydlySOCK           => micflx%RHydlySOCK,             &
    tRHydlyBioReSOM      => micflx%tRHydlyBioReSOM,        &
    tRHydlySoprtOM       => micflx%tRHydlySoprtOM,         &
    RHydlysSolidOM       => ncplxf%RHydlysSolidOM,         &
    RHumifySolidOM       => ncplxf%RHumifySolidOM,         &
    RDecmpProdDOM        => ncplxf%RDecmpProdDOM,          &
    RHydlysBioResduOM    => ncplxf%RHydlysBioResduOM,      &
    RHydlysSorptOM       => ncplxf%RHydlysSorptOM,         &
    ROQC4HeterMicActCmpK => nmicdiag%ROQC4HeterMicActCmpK, &
    BulkSOMC             => ncplxs%BulkSOMC,               &
    TOMEK                => nmicdiag%TOMEK,                &
    tMaxNActMicrbK       => ncplxs%tMaxNActMicrbK,         &
    tMaxPActMicrbK       => ncplxs%tMaxPActMicrbK,         &
    VOLWZ                => nmicdiag%VOLWZ,                &
    TSensGrowth          => nmicdiag%TSensGrowth,          &
    RHydrolysisScalCmpK  => nmicdiag%RHydrolysisScalCmpK,  &
    iprotein             => micpar%iprotein,               &
    icarbhyro            => micpar%icarbhyro,              &
    icellulos            => micpar%icellulos,              &
    ilignin              => micpar%ilignin,                &
    SPOSC                => micpar%SPOSC,                  &
    k_POM                => micpar%k_POM,                  &
    CNRH                 => micpar%CNRH,                   &
    CPRH                 => micpar%CPRH,                   &
    CDOM                 => ncplxs%CDOM,                   &
    EPOC                 => micstt%EPOC,                   &
    CNOSC                => micstt%CNOSC,                  &
    CPOSC                => micstt%CPOSC,                  &
    SorbedOM             => micstt%SorbedOM,               &
    SolidOM              => micstt%SolidOM,                &
    SolidOMAct           => micstt%SolidOMAct,             &
    OMBioResdu           => micstt%OMBioResdu,             &
    VLSoilMicP           => micfor%VLSoilMicP,             &
    ZEROS                => micfor%ZEROS,                  &
    ZEROS2               => micfor%ZEROS2,                 &
    litrm                => micfor%litrm,                  &
    SoilMicPMassLayer    => micfor%SoilMicPMassLayer       &
  )
  call PrintInfo('beg '//subname)
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

    IF(TOMEK(ielmc,K).GT.ZEROS)THEN
      CNOMX   = TOMEK(ielmn,K)/tMaxNActMicrbK(K)
      CPOMX   = TOMEK(ielmp,K)/tMaxPActMicrbK(K)
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
      RHydrolysisScalCmpK(K)=AZMAX1(ROQC4HeterMicActCmpK(K))*DFNS*OQCI/BulkSOMC(K)
      
      D785: DO M=1,jsken
        IF(SolidOM(ielmc,M,K).GT.ZEROS)THEN                  
          RHydlysSolidOM(ielmc,M,K) = SolidOMAct(M,K)*AMIN1(0.5_r8,SPOSC(M,K)*RHydrolysisScalCmpK(K))            
          dHyd=AMIN1(1._r8,RHydlysSolidOM(ielmc,M,K)/SolidOM(ielmc,M,K))      
          RHydlysSolidOM(ielmn,M,K) = AZERO(SolidOM(ielmn,M,K))*dHyd/FCNK(K)          
          RHydlysSolidOM(ielmp,M,K) = AZERO(SolidOM(ielmp,M,K))*dHyd/FCPK(K)
                
          DO NE=1,NumPlantChemElms
            RHydlysSolidOM(NE,M,K) = AMIN1(RHydlysSolidOM(NE,M,K),SolidOM(NE,M,K))
            tRHydlySOM(NE)         = tRHydlySOM(NE)+RHydlysSolidOM(NE,M,K)
          ENDDO
          RHydlySOCK(K)=RHydlySOCK(K)+RHydlysSolidOM(ielmc,M,K)
        ELSE

          DO NE    = 1, NumPlantChemElms
            RHydlysSolidOM(NE,M,K)=0.0_r8
          ENDDO
        ENDIF
      ENDDO D785
      call PrintInfo('end 785')      
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
        !litter complexes
        RHumifySolidOM(ielmc,ilignin,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmn,ilignin,K)/CNRH(k_POM) &
          ,RHydlysSolidOM(ielmp,ilignin,K)/CPRH(k_POM),EPOC*RHydlysSolidOM(ielmc,ilignin,K)))
        RHOSCM                           = 0.10_r8*RHumifySolidOM(ielmc,ilignin,K)
        RHumifySolidOM(ielmc,iprotein,K) = AZMAX1(AMIN1(RHydlysSolidOM(ielmc,iprotein,K) &
          ,RHydlysSolidOM(ielmn,iprotein,K)/CNRH(k_POM),RHydlysSolidOM(ielmp,iprotein,K)/CPRH(k_POM),RHOSCM))
        RHumifySolidOM(ielmc,icarbhyro,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,icarbhyro,K) &
          ,RHydlysSolidOM(ielmn,icarbhyro,K)/CNRH(k_POM),RHydlysSolidOM(ielmp,icarbhyro,K)/CPRH(k_POM),RHOSCM))
        RHumifySolidOM(ielmc,icellulos,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,icellulos,K) &
          ,RHydlysSolidOM(ielmn,icellulos,K)/CNRH(k_POM),RHydlysSolidOM(ielmp,icellulos,K)/CPRH(k_POM),RHOSCM-RHumifySolidOM(ielmc,icarbhyro,K)))

        D805: DO M=1,jsken
          RHumifySolidOM(ielmn,M,K) = AMIN1(RHydlysSolidOM(ielmn,M,K),RHumifySolidOM(ielmc,M,K)*CNRH(k_POM))
          RHumifySolidOM(ielmp,M,K) = AMIN1(RHydlysSolidOM(ielmp,M,K),RHumifySolidOM(ielmc,M,K)*CPRH(k_POM))
          DO NE=1,NumPlantChemElms
            RDecmpProdDOM(NE,M,K)=RHydlysSolidOM(NE,M,K)-RHumifySolidOM(NE,M,K)
          ENDDO
        ENDDO D805
      ELSE
        !non-litter complexes
        D810: DO M=1,jsken
          DO NE=1,NumPlantChemElms      
            RHumifySolidOM(NE,M,K) = 0.0_r8
            RDecmpProdDOM(NE,M,K)  = RHydlysSolidOM(NE,M,K)
          ENDDO
        ENDDO D810
      ENDIF
      call PrintInfo('end 805')      
    ELSE
      D780: DO M=1,jsken
        DO NE=1,NumPlantChemElms    
          RHydlysSolidOM(NE,M,K) = 0.0_r8
          RHumifySolidOM(NE,M,K) = 0.0_r8
          RDecmpProdDOM(NE,M,K)  = 0.0_r8
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
      RHydrolysisScalCmpK(K)=AZMAX1(ROQC4HeterMicActCmpK(K))*DFNS*OQCI/BulkSOMC(K)
      D775: DO M=1,ndbiomcp
        IF(OMBioResdu(ielmc,M,K).GT.ZEROS)THEN
          dHyd=AMIN1(1._r8,SPORC(M)*RHydrolysisScalCmpK(K))
          RHydlysBioResduOM(ielmc,M,K) = AZERO(OMBioResdu(ielmc,M,K))*dHyd
          RHydlysBioResduOM(ielmn,M,K) = AZERO(OMBioResdu(ielmn,M,K))*dHyd/FCNK(K)
          RHydlysBioResduOM(ielmp,M,K) = AZERO(OMBioResdu(ielmp,M,K))*dHyd/FCPK(K)

          DO NE=1,NumPlantChemElms
            RHydlysBioResduOM(NE,M,K) = AMIN1(RHydlysBioResduOM(NE,M,K),OMBioResdu(NE,M,K))
            tRHydlyBioReSOM(NE)       = tRHydlyBioReSOM(NE)+RHydlysBioResduOM(NE,M,K)
          ENDDO
        ELSE
          DO NE=1,NumPlantChemElms      
            RHydlysBioResduOM(NE,M,K)=0.0_r8
          ENDDO  
        ENDIF
      ENDDO D775
      call PrintInfo('end 775')      

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
      RHydrolysisScalCmpK(K)=AZMAX1(ROQC4HeterMicActCmpK(K))*DFNS*OQCI/BulkSOMC(K)    
      IF(SorbedOM(ielmc,K).GT.ZEROS)THEN
        dHyd=AMIN1(1._r8, SPOHC*RHydrolysisScalCmpK(K))
        RHydlysSorptOM(ielmc,K)        = SorbedOM(ielmc,K)*dHyd
        RHydlysSorptOM(ielmn,K)        = AZERO(SorbedOM(ielmn,K))*dHyd/FCNK(K)
        RHydlysSorptOM(ielmp,K)        = AZERO(SorbedOM(ielmp,K))*dHyd/FCPK(K)
        RHydlysSorptOM(idom_acetate,K) = AZERO(SorbedOM(idom_acetate,K))*AZMAX1(AMIN1(1._r8,SPOHA*ROQC4HeterMicActCmpK(K)*DFNS/BulkSOMC(K)))

        DO NE=1,NumPlantChemElms
          RHydlysSorptOM(NE,K) = AMIN1(SorbedOM(NE,K),RHydlysSorptOM(NE,K))
          tRHydlySoprtOM(NE)   = tRHydlySoprtOM(NE)+RHydlysSorptOM(NE,K)
        ENDDO
      ELSE
        DO idom=idom_beg,idom_end
          RHydlysSorptOM(idom,K)=0.0_r8
        ENDDO
      ENDIF
    ELSE

      DO idom=idom_beg,idom_end          
        RHydlysSorptOM(idom,K)=0.0_r8
      ENDDO
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
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
  integer :: K,M,N,NGL,NE,idom,N1,N2
  real(r8) :: FORC(0:jcplx)
  real(r8) :: dflux,scal
!     begin_execution
  associate(                                                                    &
    k_POM                            => micpar%k_POM,                           &
    TDOMUptkHeter                    => ncplxf%TDOMUptkHeter,                   &    
    DOMuptk4GrothHeter               => nmicf%DOMuptk4GrothHeter,               &
    RMetabDOCUptkHeter               => nmicf%RMetabDOCUptkHeter,               &
    RMetabAcetUptkHeter              => nmicf%RMetabAcetUptkHeter,              &
    RkillLitrfal2ResduOMHeter        => nmicf%RkillLitrfal2ResduOMHeter,        &
    RMaintDefLitrfal2ResduOMHeter    => nmicf%RMaintDefLitrfal2ResduOMHeter,    &
    RCCMEAutor                       => nmicf%RCCMEAutor,                       &
    RAcetateProdHeter                => nmicf%RAcetateProdHeter,                &
    RCCMEHeter                       => nmicf%RCCMEHeter,                       &
    RHydlysSolidOM                   => ncplxf%RHydlysSolidOM,                  &
    RHumifySolidOM                   => ncplxf%RHumifySolidOM,                  &
    RDecmpProdDOM                    => ncplxf%RDecmpProdDOM,                   &
    RHydlysBioResduOM                => ncplxf%RHydlysBioResduOM,               &
    RHydlysSorptOM                   => ncplxf%RHydlysSorptOM,                  &
    RDOMSorp                         => ncplxf%RDOMSorp,                        &
    OMBioResduK                      => ncplxs%OMBioResduK,                     &
    TOMBioResdu                      => nmicdiag%TOMBioResdu,                   &
    RMaintDefLitrfal2ResduOMAutor    => nmicf%RMaintDefLitrfal2ResduOMAutor,    &
    RkillLitrfal2ResduOMAutor        => nmicf%RkillLitrfal2ResduOMAutor,        &
    SolidOM                          => micstt%SolidOM,                         &
    iprotein                         => micpar%iprotein,                        &
    SolidOMAct                       => micstt%SolidOMAct,                      &
    SOMPomProtein                    => micstt%SOMPomProtein,                   &
    DOM                              => micstt%DOM,                             &
    DOM_MicP_drib                    => micstt%DOM_MicP_drib ,                  &
    OMBioResdu                       => micstt%OMBioResdu,                      &
    SorbedOM                         => micstt%SorbedOM,                        &
    ZEROS                            => micfor%ZEROS,                           &
    NumHetetr1MicCmplx               => micpar%NumHetetr1MicCmplx ,             &
    is_activeMicrbFungrpHeter        => micpar%is_activeMicrbFungrpHeter,       &
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
!     RMaintDefLitrfal2ResduOMHeter,RCMMN,RMaintDefLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
!
  D1690: DO K=1,KL
    IF(TOMBioResdu.GT.ZEROS)THEN
      !partition dead autotrophs proportionally
      FORC(K)=OMBioResduK(K)/TOMBioResdu
    ELSE
      IF(K.EQ.k_POM)THEN
        FORC(K)=1.0_r8
      ELSE
        FORC(K)=0.0_r8
      ENDIF
    ENDIF
    D1685: DO N=1,NumMicbAFunGrupsPerCmplx
      if(.not.micpar%is_activeMicrbFungrpHeter(N))cycle
      D1680: DO M=1,ndbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          DO NE=1,NumPlantChemElms
            RCCMEAutor(NE,M,K)=RCCMEAutor(NE,M,K)+(RkillLitrfal2ResduOMAutor(NE,M,NGL)+RMaintDefLitrfal2ResduOMAutor(NE,M,NGL))*FORC(K)
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
        DOM(NE,K)=DOM(NE,K)+RDecmpProdDOM(NE,M,K)
      ENDDO
      !
      !     LIGNIFICATION PRODUCTS
      !
      !       RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
      !
      IF(.not.litrm)THEN
       ! add to POM carbonhydrate from humification
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
        OMBioResdu(NE,M,K) = OMBioResdu(NE,M,K)-RHydlysBioResduOM(NE,M,K)
        DOM(NE,K)          = DOM(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D575

    N1=1; N2=NumHetetr1MicCmplx
    DO  N=1,NumMicbHFunGrupsPerCmplx
      if(micpar%is_anaerobic_hetr(N))then
        DO NGL=JGniH(N),JGnfH(N)
          DOM(idom_acetate,K) = DOM(idom_acetate,K)+RAcetateProdHeter(NGL,K)
        ENDDO
      endif
    ENDDO

!    write(151,*)K,'b',DOM(idom_doc,K)+SorbedOM(idom_doc,K)+DOM(idom_acetate,K)+SorbedOM(idom_acetate,K),&
!      DOM(idom_doc,K),SorbedOM(idom_doc,K),DOM(idom_acetate,K),SorbedOM(idom_acetate,K)    
    !update adsorbed OM    
    DO idom=idom_beg,idom_end
      dflux            = RHydlysSorptOM(idom,K)
      DOM(idom,K)      = DOM(idom,K)+dflux
      SorbedOM(idom,K) = SorbedOM(idom,K)-dflux
    ENDDO
!    write(151,*)K,'e',DOM(idom_doc,K)+SorbedOM(idom_doc,K)+DOM(idom_acetate,K)+SorbedOM(idom_acetate,K), &
!      DOM(idom_doc,K),SorbedOM(idom_doc,K),DOM(idom_acetate,K),SorbedOM(idom_acetate,K)    
    !
    !     MICROBIAL UPTAKE OF DISSOLVED C, N, P
    !
    !     RMetabDOCUptkHeter,RMetabAcetUptkHeter,DOMuptk4GrothHeter,DOMuptk4GrothHeter=DOC,acetate,DON,DOP uptake
    !     RAcetateProdHeter=acetate production from fermentation
    
    call SubstrateDribbling(n1,n2,RMetabDOCUptkHeter(N1:N2,K),DOM_MicP_drib(idom_doc,K),DOM(idom_doc,K))

    call SubstrateDribbling(n1,n2,DOMuptk4GrothHeter(ielmn,N1:N2,K),DOM_MicP_drib(idom_don,K),DOM(idom_don,K))

    call SubstrateDribbling(n1,n2,DOMuptk4GrothHeter(ielmp,N1:N2,K),DOM_MicP_drib(idom_dop,K),DOM(idom_dop,K))

    call SubstrateDribbling(n1,n2,RMetabAcetUptkHeter(N1:N2,K),DOM_MicP_drib(idom_acetate,K),DOM(idom_acetate,K))
!   
!     MICROBIAL DECOMPOSITION PRODUCTS
!
    D570: DO N=1,NumMicbHFunGrupsPerCmplx
      if(.not.is_activeMicrbFungrpHeter(N))cycle
      DO NGL=JGniH(N),JGnfH(N)
        TDOMUptkHeter(idom_doc,K)     = TDOMUptkHeter(idom_doc,K)+RMetabDOCUptkHeter(NGL,K)
        TDOMUptkHeter(idom_don,K)     = TDOMUptkHeter(idom_don,K)+DOMuptk4GrothHeter(ielmn,NGL,K)
        TDOMUptkHeter(idom_dop,K)     = TDOMUptkHeter(idom_dop,K)+DOMuptk4GrothHeter(ielmp,NGL,K)
        TDOMUptkHeter(idom_acetate,K) = TDOMUptkHeter(idom_acetate,K)+RMetabAcetUptkHeter(NGL,K)        
      enddo
      D565: DO M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          OMBioResdu(NE,M,K)=OMBioResdu(NE,M,K)+RCCMEHeter(NE,M,K)+RCCMEAutor(NE,M,K)
        ENDDO
      ENDDO D565

    ENDDO D570

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
  character(len=*), parameter :: subname='HeterotrophAnabolicUpdate'
  integer  :: K,M,N,NGL,MID3,MID,NE
  real(r8) ::CGROMC,dmassC
!     begin_execution
  associate(                                                                &
    DOMuptk4GrothHeter             => nmicf%DOMuptk4GrothHeter,             &
    NonstX2stBiomHeter             => nmicf%NonstX2stBiomHeter,             &
    Resp4NFixHeter                 => nmicf%Resp4NFixHeter,                 &
    RespGrossHeter                 => nmicf%RespGrossHeter,                 &
    RNOxDOCReduxRespDenitLim       => nmicf%RNOxDOCReduxRespDenitLim,       &
    RNO3imobilSoilHeter            => nmicf%RNO3imobilSoilHeter,            &
    RCO2ProdHeter                  => nmicf%RCO2ProdHeter,                  &
    RH2PO4imobilSoilHeter          => nmicf%RH2PO4imobilSoilHeter,          &
    RNH4imobilBandHeter            => nmicf%RNH4imobilBandHeter,            &
    RNO3imobilBandHeter            => nmicf%RNO3imobilBandHeter,            &
    RH2PO4imobilBandHeter          => nmicf%RH2PO4imobilBandHeter,          &
    RkillLitrfal2HumOMHeter        => nmicf%RkillLitrfal2HumOMHeter,        &
    RMaintDefLitrfal2HumOMHeter    => nmicf%RMaintDefLitrfal2HumOMHeter,    &
    RN2FixHeter                    => nmicf%RN2FixHeter,                    &
    RKillOMHeter                   => nmicf%RKillOMHeter,                   &
    RkillRecycOMHeter              => nmicf%RkillRecycOMHeter,              &
    RMaintDefcitKillOMHeter        => nmicf%RMaintDefcitKillOMHeter,        &
    RMaintDefcitRecycOMHeter       => nmicf%RMaintDefcitRecycOMHeter,       &
    RNH4imobilLitrHeter            => nmicf%RNH4imobilLitrHeter,            &
    RNO3imobilLitrHeter            => nmicf%RNO3imobilLitrHeter,            &
    RH2PO4imobilLitrHeter          => nmicf%RH2PO4imobilLitrHeter,          &
    RH1PO4imobilSoilHeter          => nmicf%RH1PO4imobilSoilHeter,          &
    RH1PO4imobilBandHeter          => nmicf%RH1PO4imobilBandHeter,          &
    RH1PO4imobilLitrHeter          => nmicf%RH1PO4imobilLitrHeter,          &
    RNH4imobilSoilHeter            => nmicf%RNH4imobilSoilHeter,            &
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
    NetCAssimhr                    => micflx%NetCAssimhr,                   & !net C assimilation
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
      DO  N=1,NumMicbHFunGrupsPerCmplx
        if(.not.micpar%is_activeMicrbFungrpHeter(N))cycle
        DO NGL=JGniH(N),JGnfH(N)
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
!     RMaintDefLitrfal2HumOMHeter,RHMMN,RMaintDefLitrfal2HumOMHeter=transfer of senesence LitrFall C,N,P to humus
!
            IF(.not.litrm)THEN              
              DO NE=1,NumPlantChemElms
                !add as protein
                SolidOM(NE,iprotein,k_humus)=SolidOM(NE,iprotein,k_humus) &
                  +ElmAllocmatMicrblitr2POM(iprotein)*(RkillLitrfal2HumOMHeter(NE,M,NGL,K)&
                  +RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))
                !add as carbon hydro
                SolidOM(NE,icarbhyro,k_humus)=SolidOM(NE,icarbhyro,k_humus) &
                  +ElmAllocmatMicrblitr2POM(icarbhyro)*(RkillLitrfal2HumOMHeter(NE,M,NGL,K)&
                  +RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))
              ENDDO
            ELSE
              DO NE=1,NumPlantChemElms
                SOMHumProtein(NE)=SOMHumProtein(NE)+ElmAllocmatMicrblitr2POMU(iprotein) &
                  *(RkillLitrfal2HumOMHeter(NE,M,NGL,K)+RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))
                SOMHumCarbohyd(NE)=SOMHumCarbohyd(NE)+ElmAllocmatMicrblitr2POMU(icarbhyro) &
                  *(RkillLitrfal2HumOMHeter(NE,M,NGL,K)+RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))                  
              ENDDO
            ENDIF
          ENDDO D540

          !
          !     INPUTS TO NONSTRUCTURAL POOLS
          !
          !     DOMuptk4GrothHeter=total DOC+acetate uptake
          !     RespGrossHeter=total respiration, including produciton of CO2 and acetate
          !     RNOxDOCReduxRespDenitLim=DOC respiration for denitrifcation
          !     Resp4NFixHeter=respiration for N2 fixation
          !     RCO2ProdHeter=total CO2 emission
          !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
          !     R3OMC,R3OMN,RkillRecycOMHeter=microbial C,N,P recycling
          !     RMaintDefcitRecycOMHeter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
          !     DOMuptk4GrothHeter,DOMuptk4GrothHeter=DON, DOP uptake
          !     RNH4imobilSoilHeter,RNH4imobilBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
          !     RNO3imobilSoilHeter,RNO3imobilBandHeter=substrate-limited NO3 immobiln in non-band, band
          !     RH2PO4imobilSoilHeter,RH2PO4imobilBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
          !     RH1PO4imobilSoilHeter,RH1PO4imobilBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
          !     RNH4imobilLitrHeter,RNO3imobilLitrHeter =substrate-limited NH4,NO3 mineraln-immobiln
          !     RH2PO4imobilLitrHeter,RH1PO4imobilLitrHeter=substrate-limited H2PO4,HPO4 mineraln-immobiln
          !
          CGROMC               = DOMuptk4GrothHeter(ielmc,NGL,K)-RespGrossHeter(NGL,K)-RNOxDOCReduxRespDenitLim(NGL,K)-Resp4NFixHeter(NGL,K)
          RCO2ProdHeter(NGL,K) = RCO2ProdHeter(NGL,K)+Resp4NFixHeter(NGL,K)
          NetCAssimhr = NetCAssimhr+CGROMC
          GrosAssimhr = GrosAssimhr+DOMuptk4GrothHeter(ielmc,NGL,K)
          MID3        = micpar%get_micb_id(ibiom_reserve,NGL)

          D555: DO M = 1, 2
            DO NE=1,NumPlantChemElms
              mBiomeHeter(NE,MID3,K)=mBiomeHeter(NE,MID3,K)-NonstX2stBiomHeter(NE,M,NGL,K)+RkillRecycOMHeter(NE,M,NGL,K)
            ENDDO
            DO NE=2,NumPlantChemElms
              mBiomeHeter(NE,MID3,K)=mBiomeHeter(NE,MID3,K)+RMaintDefcitRecycOMHeter(NE,M,NGL,K)
            ENDDO
            !respire
            RCO2ProdHeter(NGL,K)=RCO2ProdHeter(NGL,K)+RMaintDefcitRecycOMHeter(ielmc,M,NGL,K)
          ENDDO D555

          mBiomeHeter(ielmc,MID3,K) = AZERO(mBiomeHeter(ielmc,MID3,K)+CGROMC)
          mBiomeHeter(ielmn,MID3,K) = mBiomeHeter(ielmn,MID3,K)+DOMuptk4GrothHeter(ielmn,NGL,K) &
            +RNH4imobilSoilHeter(NGL,K)+RNH4imobilBandHeter(NGL,K)+RNO3imobilSoilHeter(NGL,K) &
            +RNO3imobilBandHeter(NGL,K)+RN2FixHeter(NGL,K)
          mBiomeHeter(ielmp,MID3,K)=mBiomeHeter(ielmp,MID3,K)+DOMuptk4GrothHeter(ielmp,NGL,K) &
            +RH2PO4imobilSoilHeter(NGL,K)+RH2PO4imobilBandHeter(NGL,K)+RH1PO4imobilSoilHeter(NGL,K) &
            +RH1PO4imobilBandHeter(NGL,K)
         
          !fix negative microbial N  by immobilization
          if(mBiomeHeter(ielmn,MID3,K)<0._r8)then
            RNH4imobilSoilHeter(NGL,K) = RNH4imobilSoilHeter(NGL,K)-mBiomeHeter(ielmn,MID3,K)
            NetNH4Mineralize           = NetNH4Mineralize-mBiomeHeter(ielmn,MID3,K)
            mBiomeHeter(ielmn,MID3,K)  = 0._r8
          endif        

          !fix negative P biomass by immobilization
          if(mBiomeHeter(ielmp,MID3,K)<0._r8)then
            RH2PO4imobilSoilHeter(NGL,K) = RH2PO4imobilSoilHeter(NGL,K)-mBiomeHeter(ielmp,MID3,K)
            NetPO4Mineralize             = NetPO4Mineralize-mBiomeHeter(ielmp,MID3,K)
            mBiomeHeter(ielmp,MID3,K)    = 0._r8
          endif  
          IF(litrm)THEN
            mBiomeHeter(ielmn,MID3,K)=mBiomeHeter(ielmn,MID3,K)+RNH4imobilLitrHeter(NGL,K)+RNO3imobilLitrHeter(NGL,K)
            mBiomeHeter(ielmp,MID3,K)=mBiomeHeter(ielmp,MID3,K)+RH2PO4imobilLitrHeter(NGL,K)+RH1PO4imobilLitrHeter(NGL,K)
          ENDIF
        enddo
      ENDDO
    ENDIF
  ENDDO D550
  end associate
  end subroutine HeterotrophAnabolicUpdate
!------------------------------------------------------------------------------------------

  subroutine MicrobialLitterColonization(I,J,KL,micfor,micstt,ncplxf,ncplxs,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  type(Microbe_Diag_type),intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='MicrobialLitterColonization'

  integer  :: K,M
  real(r8) :: DOSAK

!     begin_execution
  associate(                                               &
    ROQC4HeterMicActCmpK => nmicdiag%ROQC4HeterMicActCmpK, &
    SolidOMCK            => ncplxs%SolidOMCK,              &
    SolidOMActK          => ncplxs%SolidOMActK,            &
    ZEROS                => micfor%ZEROS,                  &
    SolidOMAct           => micstt%SolidOMAct,             &
    SolidOM              => micstt%SolidOM,                &
    DOSA                 => micpar%DOSA                    &
  )

  call PrintInfo('beg '//subname)
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
  call PrintInfo('end '//subname)
  end associate
  end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

  subroutine AggregateTransfOMBioResdue(KL,micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
  implicit none
  integer,intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type),intent(in) :: nmicdiag
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(micfluxtype), intent(inout) :: micflx
  integer  :: K,M,N,NGL,NE,idom
!     begin_execution
  associate(                                                       &
    DOMuptk4GrothHeter        => nmicf%DOMuptk4GrothHeter,         &
    RMetabDOCUptkHeter        => nmicf%RMetabDOCUptkHeter,         &
    RMetabAcetUptkHeter       => nmicf%RMetabAcetUptkHeter,        &
    RO2UptkHeter              => nmicf%RO2UptkHeter,               &
    RO2DmndHeter              => nmicf%RO2DmndHeter,               &
    RNO3ReduxHeterSoil        => nmicf%RNO3ReduxHeterSoil,         &
    RNO3ReduxHeterBand        => nmicf%RNO3ReduxHeterBand,         &
    RNO2ReduxHeterSoil        => nmicf%RNO2ReduxHeterSoil,         &
    RNO2ReduxHeterBand        => nmicf%RNO2ReduxHeterBand,         &
    RN2OReduxHeter            => nmicf%RN2OReduxHeter,             &
    RNOxDOCReduxRespDenitLim  => nmicf%RNOxDOCReduxRespDenitLim,   &
    RNOxAcetReduxRespDenitLim => nmicf%RNOxAcetReduxRespDenitLim,  &
    RNH4imobilSoilHeter       => nmicf%RNH4imobilSoilHeter,        &
    RNO3imobilSoilHeter       => nmicf%RNO3imobilSoilHeter,        &
    RH2PO4imobilSoilHeter     => nmicf%RH2PO4imobilSoilHeter,      &
    RNH4imobilBandHeter       => nmicf%RNH4imobilBandHeter,        &
    RNO3imobilBandHeter       => nmicf%RNO3imobilBandHeter,        &
    RH2PO4imobilBandHeter     => nmicf%RH2PO4imobilBandHeter,      &
    RH2ProdHeter              => nmicf%RH2ProdHeter,               &
    RNH4imobilLitrHeter       => nmicf%RNH4imobilLitrHeter,        &
    RNO3imobilLitrHeter       => nmicf%RNO3imobilLitrHeter,        &
    RH2PO4imobilLitrHeter     => nmicf%RH2PO4imobilLitrHeter,      &
    RCO2ProdHeter             => nmicf%RCO2ProdHeter,              &
    RAcetateProdHeter         => nmicf%RAcetateProdHeter,          &
    RCH4ProdHeter             => nmicf%RCH4ProdHeter,              &
    TRDOM2DIE                 => micflx%TRDOM2DIE,                 & !total conversion flux between organic matter and inorganic matter
    RSMetaOxidSoilAutor       => nmicf%RSMetaOxidSoilAutor,        & !catbolic oxidation
    RSMetaOxidBandAutor       => nmicf%RSMetaOxidBandAutor,        &
    RH1PO4imobilSoilHeter     => nmicf%RH1PO4imobilSoilHeter,      &
    RH1PO4imobilBandHeter     => nmicf%RH1PO4imobilBandHeter,      &
    RH1PO4imobilLitrHeter     => nmicf%RH1PO4imobilLitrHeter,      &
    RCO2XumpAutor             => nmicf%RCO2XumpAutor,              & !CO2 uptake for biomass and respiration    
    RN2FixHeter               => nmicf%RN2FixHeter,                &
    RDecmpProdDOM             => ncplxf%RDecmpProdDOM,             &
    RHydlysBioResduOM         => ncplxf%RHydlysBioResduOM,         &
    RHydlysSorptOM            => ncplxf%RHydlysSorptOM,            &
    RDOMSorp                  => ncplxf%RDOMSorp,                  &
    TSensGrowth               => nmicdiag%TSensGrowth,             &
    RH2UptkAutor              => nmicdiag%RH2UptkAutor,            &
    RNO2ReduxSoilChemo        => naqfdiag%RNO2ReduxSoilChemo,      &
    RNO2ReduxBandChemo        => naqfdiag%RNO2ReduxBandChemo,      &
    RN2OProdSoilChemo         => naqfdiag%RN2OProdSoilChemo,       &
    RN2OProdBandChemo         => naqfdiag%RN2OProdBandChemo,       &
    RNO3ProdSoilChemo         => naqfdiag%RNO3ProdSoilChemo,       &
    RNO3ProdBandChemo         => naqfdiag%RNO3ProdBandChemo,       &
    VOLWZ                     => nmicdiag%VOLWZ,                   &
    DOMuptk4GrothAutor        => nmicf%DOMuptk4GrothAutor,         &
    RNOxReduxAutorBand        => nmicf%RNOxReduxAutorBand,         &
    RNO3UptkAutor             => nmicf%RNO3UptkAutor,              &
    RCH4ProdAutor             => nmicf%RCH4ProdAutor,              &
    RGrowthCAutor             => nmicf%RGrowthCAutor,              & !growth C biomass    
    RNOxReduxAutorSoil        => nmicf%RNOxReduxAutorSoil,         &
    RO2UptkAutor              => nmicf%RO2UptkAutor,               &
    RNOxReduxRespAutorLim     => nmicf%RNOxReduxRespAutorLim,      &
    RCO2ProdAutor             => nmicf%RCO2ProdAutor,              &
    RH1PO4TransfLitrAutor     => nmicf%RH1PO4TransfLitrAutor,      &
    RH2PO4TransfLitrAutor     => nmicf%RH2PO4TransfLitrAutor,      &
    RNO3TransfLitrAutor       => nmicf%RNO3TransfLitrAutor,        &
    RNH4TransfLitrAutor       => nmicf%RNH4TransfLitrAutor,        &
    RN2FixAutor               => nmicf%RN2FixAutor,                &
    RH1PO4TransfBandAutor     => nmicf%RH1PO4TransfBandAutor,      &
    RH2PO4TransfBandAutor     => nmicf%RH2PO4TransfBandAutor,      &
    RNO3TransfBandAutor       => nmicf%RNO3TransfBandAutor,        &
    RNH4TransfBandAutor       => nmicf%RNH4TransfBandAutor,        &
    RH2PO4TransfSoilAutor     => nmicf%RH2PO4TransfSoilAutor,      &
    RNH4TransfSoilAutor       => nmicf%RNH4TransfSoilAutor,        &
    RNO3TransfSoilAutor       => nmicf%RNO3TransfSoilAutor,        &
    RH1PO4TransfSoilAutor     => nmicf%RH1PO4TransfSoilAutor,      &
    TSens4MicbGrwoth          => micstt%TSens4MicbGrwoth,          &
    VWatMicrobAct             => micstt%VWatMicrobAct,             &
    litrm                     => micfor%litrm,                     &
    Lsurf                     => micfor%Lsurf,                     &
    k_POM                     => micpar%k_POM,                     &
    k_humus                   => micpar%k_humus,                   &
    mid_AutoAmmoniaOxidBacter => micpar%mid_AutoAmmoniaOxidBacter, &
    mid_AutoAeroCH4OxiBacter  => micpar%mid_AutoAeroCH4OxiBacter,  &
    mid_AutoNitriteOxidBacter => micpar%mid_AutoNitriteOxidBacter, &
    is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, &
    mid_AutoAMONC10           => micpar%mid_AutoAMONC10,           &
    mid_AutoAMOANME2D         => micpar%mid_AutoAMOANME2D    ,     &
    RCH4MetaDmndAutor         => micflx%RCH4MetaDmndAutor,         &
    RCH4UptkAutor             => micflx%RCH4UptkAutor,             &
    RCO2NetUptkMicb           => micflx%RCO2NetUptkMicb,           &
    RH2NetUptkMicb            => micflx%RH2NetUptkMicb,            &
    RN2NetUptkMicb            => micflx%RN2NetUptkMicb,            &
    RN2ONetUptkMicb           => micflx%RN2ONetUptkMicb,           &
    RO2UptkMicb               => micflx%RO2UptkMicb,               &
    RH1PO4MicbReliz2Band      => micflx%RH1PO4MicbReliz2Band,      &
    RH1PO4MicbReliz2Soil      => micflx%RH1PO4MicbReliz2Soil,      &
    RH2PO4MicbReliz2Band      => micflx%RH2PO4MicbReliz2Band,      &
    RH2PO4MicbReliz2Soil      => micflx%RH2PO4MicbReliz2Soil,      &
    MicrbN2Fix                => micflx%MicrbN2Fix,                &
    RNH4MicbReliz2Band        => micflx%RNH4MicbReliz2Band,        &
    RNH4MicbReliz2Soil        => micflx%RNH4MicbReliz2Soil,        &
    RNO2MicbReliz2Band        => micflx%RNO2MicbReliz2Band,        &
    RNO2MicbReliz2Soil        => micflx%RNO2MicbReliz2Soil,        &
    RNO3MicbReliz2Band        => micflx%RNO3MicbReliz2Band,        &
    RNO3MicbReliz2Soil        => micflx%RNO3MicbReliz2Soil,        &
    REcoDOMProd               => micflx%REcoDOMProd                &
  )
  D650: DO K=1,KL
    IF(.not.litrm .OR. (K.NE.k_POM .AND. K.NE.k_humus))THEN
      DO N=1,NumMicbHFunGrupsPerCmplx
        DO NGL=JGniH(N),JGnfH(N)
          naqfdiag%tRNH4MicrbImobilSoil   = naqfdiag%tRNH4MicrbImobilSoil+RNH4imobilSoilHeter(NGL,K)
          naqfdiag%tRNO3MicrbImobilSoil   = naqfdiag%tRNO3MicrbImobilSoil+RNO3imobilSoilHeter(NGL,K)
          naqfdiag%tRH2PO4MicrbImobilSoil = naqfdiag%tRH2PO4MicrbImobilSoil+RH2PO4imobilSoilHeter(NGL,K)
          naqfdiag%tRH1PO4MicrbImobilSoil = naqfdiag%tRH1PO4MicrbImobilSoil+RH1PO4imobilSoilHeter(NGL,K)  !> 0 uptake 
          naqfdiag%tRNH4MicrbImobilBand   = naqfdiag%tRNH4MicrbImobilBand+RNH4imobilBandHeter(NGL,K)
          naqfdiag%tRNO3MicrbImobilBand   = naqfdiag%tRNO3MicrbImobilBand+RNO3imobilBandHeter(NGL,K)
          naqfdiag%tRH2PO4MicrbImobilBand = naqfdiag%tRH2PO4MicrbImobilBand+RH2PO4imobilBandHeter(NGL,K)
          naqfdiag%tRH1PO4MicrbImobilBand = naqfdiag%tRH1PO4MicrbImobilBand+RH1PO4imobilBandHeter(NGL,K)
          naqfdiag%TFixN2                 = naqfdiag%TFixN2+RN2FixHeter(NGL,K)
          IF(litrm)THEN
            micflx%tRNH4MicrbImobilSoil   = micflx%tRNH4MicrbImobilSoil+RNH4imobilLitrHeter(NGL,K)
            micflx%tRNO3MicrbImobilSoil   = micflx%tRNO3MicrbImobilSoil+RNO3imobilLitrHeter(NGL,K)
            micflx%tRH2PO4MicrbImobilSoil = micflx%tRH2PO4MicrbImobilSoil+RH2PO4imobilLitrHeter(NGL,K)
            micflx%tRH1PO4MicrbImobilSoil = micflx%tRH1PO4MicrbImobilSoil+RH1PO4imobilLitrHeter(NGL,K)
          ENDIF
      !
          naqfdiag%tRCO2MicrbProd     = naqfdiag%tRCO2MicrbProd+RCO2ProdHeter(NGL,K)
          naqfdiag%tRCH4MicrbProd     = naqfdiag%tRCH4MicrbProd+RCH4ProdHeter(NGL,K)
          naqfdiag%tRNOxMicrbRedux    = naqfdiag%tRNOxMicrbRedux+RNOxDOCReduxRespDenitLim(NGL,K)+RNOxAcetReduxRespDenitLim(NGL,K)
          naqfdiag%tRO2MicrbUptk      = naqfdiag%tRO2MicrbUptk+RO2UptkHeter(NGL,K)
          naqfdiag%TReduxNO3Soil      = naqfdiag%TReduxNO3Soil+RNO3ReduxHeterSoil(NGL,K)       !NO3->NO2
          naqfdiag%TReduxNO3Band      = naqfdiag%TReduxNO3Band+RNO3ReduxHeterBand(NGL,K)       !NO3->NO2
          naqfdiag%TDeniReduxNO2Soil  = naqfdiag%TDeniReduxNO2Soil+RNO2ReduxHeterSoil(NGL,K)   !NO2->N2O
          naqfdiag%TDeniReduxNO2Band  = naqfdiag%TDeniReduxNO2Band+RNO2ReduxHeterBand(NGL,K)   !NO2->N2O
          naqfdiag%TReduxNO2toN2OSoil = naqfdiag%TReduxNO2toN2OSoil+RNO2ReduxHeterSoil(NGL,K)
          naqfdiag%TReduxNO2toN2OBand = naqfdiag%TReduxNO2toN2OBand+RNO2ReduxHeterBand(NGL,K)
          naqfdiag%TReduxN2OtoN2      = naqfdiag%TReduxN2OtoN2+RN2OReduxHeter(NGL,K)               !N2O -> N2
          naqfdiag%TProdH2            = naqfdiag%TProdH2+RH2ProdHeter(NGL,K)
          naqfdiag%tRO2UptkHeterG     = naqfdiag%tRO2UptkHeterG+RO2UptkHeter(NGL,K)
          naqfdiag%tRO2DmndHeterG     = naqfdiag%tRO2DmndHeterG + RO2DmndHeter(NGL,K)
          nmicf%RO2UptkHeterG(NGL)    = nmicf%RO2UptkHeterG(NGL)+RO2UptkHeter(NGL,K)
          nmicf%RO2DmndHeterG(NGL)    = nmicf%RO2DmndHeterG(NGL)+RO2DmndHeter(NGL,K)
          TRDOM2DIE(ielmc)            = TRDOM2DIE(ielmc)+RCO2ProdHeter(NGL,K)+RCH4ProdHeter(NGL,K)
        ENDDO
      ENDDO
    ENDIF
  ENDDO D650

  DO  N=1,NumMicbAFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%tRNH4MicrbImobilSoil   = naqfdiag%tRNH4MicrbImobilSoil+RNH4TransfSoilAutor(NGL)
        naqfdiag%tRNO3MicrbImobilSoil   = naqfdiag%tRNO3MicrbImobilSoil+RNO3TransfSoilAutor(NGL)
        naqfdiag%tRH2PO4MicrbImobilSoil = naqfdiag%tRH2PO4MicrbImobilSoil+RH2PO4TransfSoilAutor(NGL)
        naqfdiag%tRH1PO4MicrbImobilSoil = naqfdiag%tRH1PO4MicrbImobilSoil+RH1PO4TransfSoilAutor(NGL)
        naqfdiag%tRNH4MicrbImobilBand   = naqfdiag%tRNH4MicrbImobilBand+RNH4TransfBandAutor(NGL)
        naqfdiag%tRNO3MicrbImobilBand   = naqfdiag%tRNO3MicrbImobilBand+RNO3TransfBandAutor(NGL)
        naqfdiag%tRH2PO4MicrbImobilBand = naqfdiag%tRH2PO4MicrbImobilBand+RH2PO4TransfBandAutor(NGL)
        naqfdiag%tRH1PO4MicrbImobilBand = naqfdiag%tRH1PO4MicrbImobilBand+RH1PO4TransfBandAutor(NGL)
        naqfdiag%TFixN2                 = naqfdiag%TFixN2+RN2FixAutor(NGL)

        IF(litrm)then
          micflx%tRNH4MicrbImobilSoil   = micflx%tRNH4MicrbImobilSoil+RNH4TransfLitrAutor(NGL)
          micflx%tRNO3MicrbImobilSoil   = micflx%tRNO3MicrbImobilSoil+RNO3TransfLitrAutor(NGL)
          micflx%tRH2PO4MicrbImobilSoil = micflx%tRH2PO4MicrbImobilSoil+RH2PO4TransfLitrAutor(NGL)
          micflx%tRH1PO4MicrbImobilSoil = micflx%tRH1PO4MicrbImobilSoil+RH1PO4TransfLitrAutor(NGL)
        ENDIF

        naqfdiag%tRCO2MicrbProd       = naqfdiag%tRCO2MicrbProd+RCO2ProdAutor(NGL)
        naqfdiag%tRCH4MicrbProd       = naqfdiag%tRCH4MicrbProd+RCH4ProdAutor(NGL)
        naqfdiag%tRNOxMicrbRedux      = naqfdiag%tRNOxMicrbRedux+RNOxReduxRespAutorLim(NGL)
        naqfdiag%tRO2MicrbUptk        = naqfdiag%tRO2MicrbUptk+RO2UptkAutor(NGL)
        naqfdiag%TReduxNO3Soil        = naqfdiag%TReduxNO3Soil+RNO3UptkAutor(NGL)
        naqfdiag%TNitNO2Redux2N2OSoil = naqfdiag%TNitNO2Redux2N2OSoil+RNOxReduxAutorSoil(NGL)*1.5_r8 !NO2(-) -> N2O from aerobic denitrification
        naqfdiag%TNitNO2Redux2N2OBand = naqfdiag%TNitNO2Redux2N2OBand+RNOxReduxAutorBand(NGL)*1.5_r8 !NO2(-) -> N2O from aerobic denitrification
        naqfdiag%TReduxNO2toN2OSoil        = naqfdiag%TReduxNO2toN2OSoil+RNOxReduxAutorSoil(NGL)*1.5_r8
        naqfdiag%TReduxNO2toN2OBand   = naqfdiag%TReduxNO2toN2OBand+RNOxReduxAutorBand(NGL)*1.5_r8

      ENDDO
    ENDIF
  ENDDO

  IF(Lsurf)THEN
    naqfdiag%tRNH4MicrbImobilSoil   = naqfdiag%tRNH4MicrbImobilSoil+micfor%tRNH4MicrbImobilSoil
    naqfdiag%tRNO3MicrbImobilSoil   = naqfdiag%tRNO3MicrbImobilSoil+micfor%tRNO3MicrbImobilSoil
    naqfdiag%tRH2PO4MicrbImobilSoil = naqfdiag%tRH2PO4MicrbImobilSoil+micfor%tRH2PO4MicrbImobilSoil
    naqfdiag%tRH1PO4MicrbImobilSoil = naqfdiag%tRH1PO4MicrbImobilSoil+micfor%tRH1PO4MicrbImobilSoil
  ENDIF

  ! tRCO2GrothAutor=total CO2 uptake by autotrophs, ammonia oxidizer
  ! nitrite oxidizer, and hydrogenotrophic methanogens,
  ! all of which involves CO2 for both energy and C biomass.
  D645: DO N=1,NumMicbAFunGrupsPerCmplx
    IF(micpar%is_CO2_autotroph(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%tRCO2GrothAutor=naqfdiag%tRCO2GrothAutor+RCO2XumpAutor(NGL)
      ENDDO
    ENDIF
  ENDDO D645

!
!     ALLOCATE AGGREGATED TRANSFOMBioResduATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2NetUptkMicb=net CO2 uptake, < 0 producig CO2
!     tRCO2MicrbProd total CO2 emission by heterotrophs reducing O2
!     tRNOxMicrbRedux=total CO2 emission by denitrifiers reducing NOx
!     RSMetaOxidSoilAutor(3)=CH4 oxidation
!     RCH4UptkAutor=net CH4 uptake, >0, means uptake
!     DOMuptk4GrothHeter=total CH4 uptake by autotrophs
!     tRCH4MicrbProd=total CH4 emission
!     RH2NetUptkMicb=net H2 uptake
!     RH2UptkAutor,TProdH2=total H2 uptake, emission
!     RO2UptkMicb,tRO2MicrbUptk=total O2 uptake
!     RN2NetUptkMicb=total N2 production
!     TReduxN2OtoN2=total N2O reduction
!     RN2ONetUptkMicb=total N2O uptake
!     NO2(-) -> N2O
!     TReduxNO2toN2OSoil,TReduxNO2toN2OBand=total NO2 reduction in non-band,band
!     RN2OProdSoilChemo,RN2OProdBandChemo=nitrous acid reduction in non-band,band
!
  RCO2NetUptkMicb = naqfdiag%tRCO2GrothAutor-naqfdiag%tRCO2MicrbProd-naqfdiag%tRNOxMicrbRedux
  RCH4UptkAutor   = -naqfdiag%tRCH4MicrbProd

  DO  N=1,NumMicbAFunGrupsPerCmplx
    if(N.eq.mid_AutoAeroCH4OxiBacter)then
      DO NGL=JGniA(N),JGnfA(N)
        !use CH4 for both energy and biomass
        !CO2 is released from CH4 oxidation
        RCO2NetUptkMicb = RCO2NetUptkMicb-RSMetaOxidSoilAutor(NGL)
        !total CH4 uptake includes catabolic oxidation and biomass synthesis
        RCH4UptkAutor          = RCH4UptkAutor+RSMetaOxidSoilAutor(NGL)+DOMuptk4GrothAutor(ielmc,NGL)
        RCH4MetaDmndAutor(NGL) = RCH4MetaDmndAutor(NGL)+DOMuptk4GrothAutor(ielmc,NGL)
        TRDOM2DIE(ielmc)       = TRDOM2DIE(ielmc)-RGrowthCAutor(NGL)
      ENDDO
    elseif(N.eq.mid_AutoAMONC10)then
      !use CO2 from the environment for carbon biomass
      !CO2 is released from CH4 oxidation and taken up for biomass
      RCO2NetUptkMicb  = RCO2NetUptkMicb+RGrowthCAutor(NGL)-RSMetaOxidSoilAutor(NGL)
      RCH4UptkAutor    = RCH4UptkAutor+RSMetaOxidSoilAutor(NGL)
      TRDOM2DIE(ielmc) = TRDOM2DIE(ielmc)-RGrowthCAutor(NGL)
    elseif (N.eq.mid_AutoAMOANME2D)then  
      !CO2 is released from CH4 oxidation
      !CH4 is taken up for catabolic and anabolic reactions
      RCO2NetUptkMicb         = RCO2NetUptkMicb-RSMetaOxidSoilAutor(NGL)
      RCH4UptkAutor           = RCH4UptkAutor+RSMetaOxidSoilAutor(NGL)
      RCH4MetaDmndAutor(NGL)  = RCH4MetaDmndAutor(NGL)+RGrowthCAutor(NGL)
      TRDOM2DIE(ielmc)        = TRDOM2DIE(ielmc)-RGrowthCAutor(NGL)
    ENDIF
  ENDDO

  !>0. microbial uptake
  RH2NetUptkMicb  = RH2UptkAutor-naqfdiag%TProdH2
  RO2UptkMicb     = naqfdiag%tRO2MicrbUptk
  RN2NetUptkMicb  = -naqfdiag%TReduxN2OtoN2
  RN2ONetUptkMicb = -naqfdiag%TReduxNO2toN2OSoil-naqfdiag%TReduxNO2toN2OBand-RN2OProdSoilChemo &
    -RN2OProdBandChemo+naqfdiag%TReduxN2OtoN2
!
  D655: DO K=1,jcplx
    D660: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RDecmpProdDOM(NE,M,K)
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
    D670: DO N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        REcoDOMProd(idom_doc,K)     = REcoDOMProd(idom_doc,K)-RMetabDOCUptkHeter(NGL,K)
        REcoDOMProd(idom_don,K)     = REcoDOMProd(idom_don,K)-DOMuptk4GrothHeter(ielmn,NGL,K)
        REcoDOMProd(idom_dop,K)     = REcoDOMProd(idom_dop,K)-DOMuptk4GrothHeter(ielmp,NGL,K)
        REcoDOMProd(idom_acetate,K) = REcoDOMProd(idom_acetate,K)-RMetabAcetUptkHeter(NGL,K)+RAcetateProdHeter(NGL,K)
      ENDDO
    ENDDO D670
    DO idom=idom_beg,idom_end
      REcoDOMProd(idom,K)=REcoDOMProd(idom,K)-RDOMSorp(idom,K)
    ENDDO
  ENDDO D655
!
!     RNH4MicbReliz2Soil,RNH4MicbReliz2Band=net change in NH4 in band,non-band
!     tRNH4MicrbImobilSoil,tRNH4MicrbImobilBand=total NH4 mineraln-immobn in non-band,band
!     RSMetaOxidSoilAutor(1),RSMetaOxidBandAutor(1)=total NH4 oxidation in non-band,band
!     RNO3MicbReliz2Soil,RNO3MicbReliz2Band=net change in NO3 in band,non-band
!     tRNO3MicrbImobilSoil,tRNO3MicrbImobilBand=total NO3 immobn in non-band,band
!     RSMetaOxidSoilAutor(2),RSMetaOxidBandAutor(2)=total NO2 oxidation in non-band,band
!     TReduxNO3Soil,TReduxNO3Band=total NO3 reduction in non-band,band
!     RNO3ProdSoilChemo,RNO3ProdBandChemo=NO3 production from nitrous acid reduction in non-band,band
!     RNO2MicbReliz2Soil,RNO2MicbReliz2Band=net change in NO3 in band,non-band
!     TReduxNO2toN2OSoil,TReduxNO2toN2OBand=total NO2 reduction in non-band,band
!     RNO2ReduxSoilChemo,RNO2ReduxBandChemo=substrate-limited nitrous acid reduction in non-band,band
!     RH2PO4MicbReliz2Soil,RH2PO4MicbReliz2Band=net change in H2PO4 in band,non-band
!     tRH2PO4MicrbImobilSoil,tRH2PO4MicrbImobilBand=total H2PO4 mineraln-immobn in non-band,band
!     RH1PO4MicbReliz2Soil,RH1PO4MicbReliz2Band=net change in HPO4 in band,non-band
!     tRH1PO4MicrbImobilSoil,tRH1PO4MicrbImobilBand=total HPO4 mineraln-immobn in non-band,band
!     MicrbN2Fix=total N2 fixation
!     XZHYS=total H+ production
!     TFixN2=total N2 fixation
!

  micflx%TRDOM2DIE(ielmn)=micflx%TRDOM2DIE(ielmn)+naqfdiag%tRNH4MicrbImobilSoil &
    +naqfdiag%tRNH4MicrbImobilBand+naqfdiag%tRNO3MicrbImobilSoil &
    +naqfdiag%tRNO3MicrbImobilBand+naqfdiag%TFixN2
  micflx%TRDOM2DIE(ielmp)=micflx%TRDOM2DIE(ielmp)+naqfdiag%tRH1PO4MicrbImobilSoil &
    +naqfdiag%tRH2PO4MicrbImobilSoil+naqfdiag%tRH1PO4MicrbImobilBand &
    +naqfdiag%tRH2PO4MicrbImobilBand


  RNH4MicbReliz2Soil=-naqfdiag%tRNH4MicrbImobilSoil
  RNO3MicbReliz2Soil=-naqfdiag%tRNO3MicrbImobilSoil-naqfdiag%TReduxNO3Soil+RNO3ProdSoilChemo
  RNO2MicbReliz2Soil=+naqfdiag%TReduxNO3Soil-naqfdiag%TReduxNO2toN2OSoil-RNO2ReduxSoilChemo
  RH2PO4MicbReliz2Soil=-naqfdiag%tRH2PO4MicrbImobilSoil
  RH1PO4MicbReliz2Soil=-naqfdiag%tRH1PO4MicrbImobilSoil     !< 0 uptake
  RNH4MicbReliz2Band=-naqfdiag%tRNH4MicrbImobilBand
  RNO3MicbReliz2Band=-naqfdiag%tRNO3MicrbImobilBand-naqfdiag%TReduxNO3Band+RNO3ProdBandChemo
  RNO2MicbReliz2Band=naqfdiag%TReduxNO3Band-naqfdiag%TReduxNO2toN2OBand-RNO2ReduxBandChemo

  !mid_AutoAmmoniaOxidBacter=1, mid_AutoNitriteOxidBacter=2, mid_AutoAeroCH4OxiBacter=3
  DO NGL=JGniA(mid_AutoAmmoniaOxidBacter),JGnfA(mid_AutoAmmoniaOxidBacter)
    RNH4MicbReliz2Soil=RNH4MicbReliz2Soil-RSMetaOxidSoilAutor(NGL)   !some NH3 -> NO2, some NH3-> N2O
    RNH4MicbReliz2Band=RNH4MicbReliz2Band-RSMetaOxidBandAutor(NGL)

    !on mole-basis 2NO2(-) + NH3 -> 1.5N2O + 2OH(-) + 0.5H2O, RNOxReduxAutorSoil(NGL)/2._r8 goes to N2O
    RNO2MicbReliz2Soil=RNO2MicbReliz2Soil+RSMetaOxidSoilAutor(NGL)-RNOxReduxAutorSoil(NGL)/2._r8
    RNO2MicbReliz2Band=RNO2MicbReliz2Band+RSMetaOxidBandAutor(NGL)-RNOxReduxAutorBand(NGL)/2._r8
  ENDDO

  DO NGL=JGniA(mid_AutoNitriteOxidBacter),JGnfA(mid_AutoNitriteOxidBacter)
    naqfdiag%tNO2OxiAuto=naqfdiag%tNO2OxiAuto+RSMetaOxidSoilAutor(NGL)+RSMetaOxidBandAutor(NGL)
    RNO3MicbReliz2Soil=RNO3MicbReliz2Soil+RSMetaOxidSoilAutor(NGL)
    RNO2MicbReliz2Soil=RNO2MicbReliz2Soil-RSMetaOxidSoilAutor(NGL)
    RNO3MicbReliz2Band=RNO3MicbReliz2Band+RSMetaOxidBandAutor(NGL)
    RNO2MicbReliz2Band=RNO2MicbReliz2Band-RSMetaOxidBandAutor(NGL)
  ENDDO

  RH2PO4MicbReliz2Band = -naqfdiag%tRH2PO4MicrbImobilBand
  RH1PO4MicbReliz2Band = -naqfdiag%tRH1PO4MicrbImobilBand
  MicrbN2Fix           = naqfdiag%TFixN2
  !how to aggregate the mean T sensitivity?
  TSens4MicbGrwoth     = TSensGrowth
  VWatMicrobAct        = VOLWZ

  DO NGL=JGNiA(mid_AutoAeroCH4OxiBacter),jGnfA(mid_AutoAeroCH4OxiBacter)
    naqfdiag%tCH4OxiAero=naqfdiag%tCH4OxiAero+RSMetaOxidSoilAutor(NGL)
  ENDDO

  DO NGL=JGNiA(mid_AutoAMONC10),jGnfA((mid_AutoAMONC10))
    naqfdiag%tCH4OxiANMO=naqfdiag%tCH4OxiANMO+RSMetaOxidSoilAutor(NGL)
  ENDDO

  DO NGL=JGNiA(mid_AutoAMOANME2D),jGnfA(mid_AutoAMOANME2D)
    naqfdiag%tCH4OxiANMO=naqfdiag%tCH4OxiANMO+RSMetaOxidSoilAutor(NGL)    
  ENDDO  
  naqfdiag%tRNH3Oxi=nmicf%RTotNH3OxidSoilAutor+nmicf%RTotNH3OxidBandAutor

  end associate
  end subroutine AggregateTransfOMBioResdue
!------------------------------------------------------------------------------------------

  subroutine SubstrateAttenf4Compet(NGL,N,K,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    micfor,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(out):: FNH4X,FNB4X            !fraction of NH4 demand over all plant+microbial demand, soil/band 
  real(r8),intent(out) :: FNO3X,FNB3X            !fraction of NO3 demand over all plant+microbial demand, soil/band
  real(r8),intent(out) :: FPO4X,FPOBX            !fraction of H2PO4 demand over all plant+microbial demand, soil/band
  real(r8), intent(out):: FP14X,FP1BX            !fraction of H1PO4 demand over all plant+microbial demand, soil/band
 
  type(micforctype), intent(in) :: micfor
  type(Cumlate_Flux_Diag_type),INTENT(INOUT)::  naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                                                   &
    FracOMActHeter          => nmics%FracOMActHeter,           &
    FracHeterBiomOfActK     => nmics%FracHeterBiomOfActK,      &
    AttenfNH4Heter          => micflx%AttenfNH4Heter,          &
    AttenfNO3Heter          => micflx%AttenfNO3Heter,          &
    AttenfH1PO4Heter        => micflx%AttenfH1PO4Heter,        &
    AttenfH2PO4Heter        => micflx%AttenfH2PO4Heter,        &
    RO2DmndHetertPrev       => micflx%RO2DmndHetertPrev,       &
    RNH4DmndSoilHeterPrev   => micflx%RNH4DmndSoilHeterPrev,   &
    RNH4DmndBandHeterPrev   => micflx%RNH4DmndBandHeterPrev,   &
    RNO3DmndSoilHeterPrev   => micflx%RNO3DmndSoilHeterPrev,   &
    RNO3DmndBandHeterPrev   => micflx%RNO3DmndBandHeterPrev,   &
    RH2PO4DmndSoilHeterPrev => micflx%RH2PO4DmndSoilHeterPrev, &
    RH2PO4DmndBandHeterPrev => micflx%RH2PO4DmndBandHeterPrev, &
    RH1PO4DmndSoilHeterPrev => micflx%RH1PO4DmndSoilHeterPrev, &
    RH1PO4DmndBandHeterPrev => micflx%RH1PO4DmndBandHeterPrev, &
    RNH4DmndLitrHeterPrev   => micflx%RNH4DmndLitrHeterPrev,   &
    RNO3DmndLitrHeterPrev   => micflx%RNO3DmndLitrHeterPrev,   &
    RH2PO4DmndLitrHeterPrev => micflx%RH2PO4DmndLitrHeterPrev, &
    RH1PO4DmndLitrHeterPrev => micflx%RH1PO4DmndLitrHeterPrev, &
    litrm                   => micfor%litrm,                   &
    VLNH4                   => micfor%VLNH4,                   &
    VLNHB                   => micfor%VLNHB,                   &
    VLNOB                   => micfor%VLNOB,                   &
    VLNO3                   => micfor%VLNO3,                   &
    VLPOB                   => micfor%VLPOB,                   &
    VLPO4                   => micfor%VLPO4,                   &
    ZEROS                   => micfor%ZEROS,                   &
    RNH4EcoDmndSoilPrev     => micfor%RNH4EcoDmndSoilPrev,     &
    RNH4EcoDmndBandPrev     => micfor%RNH4EcoDmndBandPrev,     &
    RNO3EcoDmndSoilPrev     => micfor%RNO3EcoDmndSoilPrev,     &
    RNH4EcoDmndLitrPrev     => micfor%RNH4EcoDmndLitrPrev,     &
    RNO3EcoDmndLitrPrev     => micfor%RNO3EcoDmndLitrPrev,     &
    RH1PO4EcoDmndLitrPrev   => micfor%RH1PO4EcoDmndLitrPrev,   &
    RH2PO4EcoDmndLitrPrev   => micfor%RH2PO4EcoDmndLitrPrev,   &
    RNO3EcoDmndBandPrev     => micfor%RNO3EcoDmndBandPrev,     &
    RH2PO4EcoDmndSoilPrev   => micfor%RH2PO4EcoDmndSoilPrev,   &
    RH2PO4EcoDmndBandPrev   => micfor%RH2PO4EcoDmndBandPrev,   &
    RH1PO4EcoDmndSoilPrev   => micfor%RH1PO4EcoDmndSoilPrev,   &
    RH1PO4EcoDmndBandPrev   => micfor%RH1PO4EcoDmndBandPrev,   &
    RDOMEcoDmndPrev         => micfor%RDOMEcoDmndPrev,         &
    RAcetateEcoDmndPrev     => micfor%RAcetateEcoDmndPrev,     &
    Lsurf                   => micfor%Lsurf,                   &
    SoilMicPMassLayer0      => micfor%SoilMicPMassLayer0       &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4DmndSoilHeterPrev(NGL,K)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4DmndBandHeterPrev(NGL,K)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNHB)
  ENDIF
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3DmndSoilHeterPrev(NGL,K)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3DmndBandHeterPrev(NGL,K)/RNO3EcoDmndBandPrev)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
  ENDIF
  IF(RH2PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4DmndSoilHeterPrev(NGL,K)/RH2PO4EcoDmndSoilPrev)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4DmndBandHeterPrev(NGL,K)/RH2PO4EcoDmndBandPrev)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4DmndSoilHeterPrev(NGL,K)/RH1PO4EcoDmndSoilPrev)
  ELSE
    FP14X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4DmndBandHeterPrev(NGL,K)/RH1PO4EcoDmndBandPrev)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
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
  IF(litrm)THEN
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNH4Heter(NGL,K)=AMAX1(FMN,RNH4DmndLitrHeterPrev(NGL,K)/RNH4EcoDmndLitrPrev)
    ELSE
      AttenfNH4Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNO3Heter(NGL,K)=AMAX1(FMN,RNO3DmndLitrHeterPrev(NGL,K)/RNO3EcoDmndLitrPrev)
    ELSE
      AttenfNO3Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH2PO4Heter(NGL,K)=AMAX1(FMN,RH2PO4DmndLitrHeterPrev(NGL,K)/RH2PO4EcoDmndLitrPrev)
    ELSE
      AttenfH2PO4Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH1PO4Heter(NGL,K)=AMAX1(FMN,RH1PO4DmndLitrHeterPrev(NGL,K)/RH1PO4EcoDmndLitrPrev)
    ELSE
      AttenfH1PO4Heter(NGL,K)=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
    ENDIF
  ENDIF

!diagnostics off
!  IF(Lsurf .AND. K.NE.micpar%k_POM .AND. K.NE.micpar%k_humus .AND. SoilMicPMassLayer0.GT.ZEROS)THEN
!    naqfdiag%TFNH4X=naqfdiag%TFNH4X+micfor%AttenfNH4HeterR(NGL,K)
!    naqfdiag%TFNO3X=naqfdiag%TFNO3X+micfor%AttenfNO3HeterR(NGL,K)
!    naqfdiag%TFPO4X=naqfdiag%TFPO4X+micfor%AttenfH2PO4HeterR(NGL,K)
!    naqfdiag%TFP14X=naqfdiag%TFP14X+micfor%AttenfH1PO4HeterR(NGL,K)
!  ENDIF
  end associate
  end subroutine SubstrateAttenf4Compet

!------------------------------------------------------------------------------------------

  subroutine BiomassMineralization(NGL,N,K,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
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
  character(len=*), parameter :: subname='BiomassMineralization'

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
   AttenfNH4Heter                  => micflx%AttenfNH4Heter,       &
   AttenfNO3Heter                  => micflx%AttenfNO3Heter,       &
   AttenfH2PO4Heter                => micflx%AttenfH2PO4Heter,     &
   AttenfH1PO4Heter                => micflx%AttenfH1PO4Heter,     &
   RNH4imobilSoilHeter             => nmicf%RNH4imobilSoilHeter,   &
   RNO3imobilSoilHeter             => nmicf%RNO3imobilSoilHeter,   &
   RH2PO4imobilSoilHeter           => nmicf%RH2PO4imobilSoilHeter, &
   RNH4imobilBandHeter             => nmicf%RNH4imobilBandHeter,   &
   RNO3imobilBandHeter             => nmicf%RNO3imobilBandHeter,   &
   RH2PO4imobilBandHeter           => nmicf%RH2PO4imobilBandHeter, &
   RNH4imobilLitrHeter             => nmicf%RNH4imobilLitrHeter,   &
   RNO3imobilLitrHeter             => nmicf%RNO3imobilLitrHeter,   &
   RH2PO4imobilLitrHeter           => nmicf%RH2PO4imobilLitrHeter, &
   RH1PO4imobilSoilHeter           => nmicf%RH1PO4imobilSoilHeter, &
   RH1PO4imobilBandHeter           => nmicf%RH1PO4imobilBandHeter, &
   RH1PO4imobilLitrHeter           => nmicf%RH1PO4imobilLitrHeter, &
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
  call PrintInfo('beg '//subname)
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
  !     RNH4imobilSoilHeter,RNH4imobilBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
  !     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
  ! update may be needed, May 17th, 2023, jyt.
  FNH4S     = VLNH4
  FNHBS     = VLNHB
  MID3      = micpar%get_micb_id(ibiom_reserve,NGL)
  RINHP     = (mBiomeHeter(ielmc,MID3,K)*rNCOMC(ibiom_reserve,NGL,K)-mBiomeHeter(ielmn,MID3,K))
  RNiDemand = RNiDemand+RINHP
  
  IF(RINHP.GT.0.0_r8)THEN
    !immobilization
    CNH4X                      = AZMAX1(CNH4S-Z4MN)
    CNH4Y                      = AZMAX1(CNH4B-Z4MN)
    RINHX                      = AMIN1(RINHP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*Z4MX)
    RNH4DmndSoilHeter(NGL,K)   = FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RNH4DmndBandHeter(NGL,K)   = FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M                      = Z4MN*VOLWU*FNH4S
    ZNHBM                      = Z4MN*VOLWU*FNHBS
    RNH4imobilSoilHeter(NGL,K) = AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4DmndSoilHeter(NGL,K))
    RNH4imobilBandHeter(NGL,K) = AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4DmndBandHeter(NGL,K))
    !mineralization  (<0._r8)
  ELSE
    !mineralization
    RNH4DmndSoilHeter(NGL,K)   = 0.0_r8
    RNH4DmndBandHeter(NGL,K)   = 0.0_r8
    RNH4imobilSoilHeter(NGL,K) = RINHP*FNH4S
    RNH4imobilBandHeter(NGL,K) = RINHP*FNHBS
  ENDIF

  NetNH4Mineralize=NetNH4Mineralize+(RNH4imobilSoilHeter(NGL,K)+RNH4imobilBandHeter(NGL,K))
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
  !     RNO3imobilSoilHeter,RNO3imobilBandHeter=substrate-limited NO3 immobiln in non-band, band
  !     NetNH4Mineralize=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
  !
  FNO3S = VLNO3
  FNO3B = VLNOB
  RINOP = AZMAX1(RINHP-RNH4imobilSoilHeter(NGL,K)-RNH4imobilBandHeter(NGL,K))
  !immobilization
  IF(RINOP.GT.0.0_r8)THEN
    CNO3X                      = AZMAX1(CNO3S-ZOMN)
    CNO3Y                      = AZMAX1(CNO3B-ZOMN)
    RINOX                      = AMIN1(RINOP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*ZOMX)
    RNO3DmndSoilHeter(NGL,K)   = FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RNO3DmndBandHeter(NGL,K)   = FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M                      = ZOMN*VOLWU*FNO3S
    ZNOBM                      = ZOMN*VOLWU*FNO3B
    RNO3imobilSoilHeter(NGL,K) = AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)),RNO3DmndSoilHeter(NGL,K))
    RNO3imobilBandHeter(NGL,K) = AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)),RNO3DmndBandHeter(NGL,K))
    !mineralization, 
  ELSE
    RNO3DmndSoilHeter(NGL,K)   = 0.0_r8
    RNO3DmndBandHeter(NGL,K)   = 0.0_r8
    RNO3imobilSoilHeter(NGL,K) = 0.0_r8
    RNO3imobilBandHeter(NGL,K) = 0.0_r8
  ENDIF
  NetNH4Mineralize=NetNH4Mineralize+(RNO3imobilSoilHeter(NGL,K)+RNO3imobilBandHeter(NGL,K))
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
  !     RH2PO4imobilSoilHeter,RH2PO4imobilBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
  !     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
  !
  FH2PS     = VLPO4
  FH2PB     = VLPOB
  MID3      = micpar%get_micb_id(ibiom_reserve,NGL)
  RIPOP     = (mBiomeHeter(ielmc,MID3,K)*rPCOMC(ibiom_reserve,NGL,K)-mBiomeHeter(ielmp,MID3,K))
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
    RH2PO4imobilSoilHeter(NGL,K) = AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RH2PO4DmndSoilHeter(NGL,K))
    RH2PO4imobilBandHeter(NGL,K) = AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RH2PO4DmndBandHeter(NGL,K))
  !mineralization  
  ELSE
    RH2PO4DmndSoilHeter(NGL,K)   = 0.0_r8
    RH2PO4DmndBandHeter(NGL,K)   = 0.0_r8
    RH2PO4imobilSoilHeter(NGL,K) = RIPOP*FH2PS
    RH2PO4imobilBandHeter(NGL,K) = RIPOP*FH2PB
  ENDIF
  NetPO4Mineralize=NetPO4Mineralize+(RH2PO4imobilSoilHeter(NGL,K)+RH2PO4imobilBandHeter(NGL,K))
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
  !     RH1PO4imobilSoilHeter,RH1PO4imobilBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band uptake (> 0)
  !     NetPO4Mineralize=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
  !
  FH1PS=VLPO4
  FH1PB=VLPOB
  ! why 0.1 here?  
  RIP1P=0.1_r8*AZMAX1(RIPOP-RH2PO4imobilSoilHeter(NGL,K)-RH2PO4imobilBandHeter(NGL,K))
  !immobilization
  IF(RIP1P.GT.0.0_r8)THEN
    CH1PX                        = AZMAX1(CH1P4-HPMN)
    CH1PY                        = AZMAX1(CH1P4B-HPMN)
    RIP1X                        = AMIN1(RIP1P,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX)
    RH1PO4DmndSoilHeter(NGL,K)   = FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RH1PO4DmndBandHeter(NGL,K)   = FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM                        = HPMN*VLWatMicP*FH1PS
    H1PBM                        = HPMN*VLWatMicP*FH1PB
    RH1PO4imobilSoilHeter(NGL,K) = AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RH1PO4DmndSoilHeter(NGL,K))
    RH1PO4imobilBandHeter(NGL,K) = AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RH1PO4DmndBandHeter(NGL,K))
  !mineralization  
  ELSE
    RH1PO4DmndSoilHeter(NGL,K)   = 0.0_r8
    RH1PO4DmndBandHeter(NGL,K)   = 0.0_r8
    RH1PO4imobilSoilHeter(NGL,K) = 0.0_r8
    RH1PO4imobilBandHeter(NGL,K) = 0.0_r8
  ENDIF
  NetPO4Mineralize=NetPO4Mineralize+(RH1PO4imobilSoilHeter(NGL,K)+RH1PO4imobilBandHeter(NGL,K))
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
  !     RNH4imobilLitrHeter=substrate-limited NH4 mineraln-immobiln
  !     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
  !when there is not sufficient uptake from the litter layer, take it from soil
  IF(litrm)THEN
    RINHPR=RINHP-RNH4imobilSoilHeter(NGL,K)-RNO3imobilSoilHeter(NGL,K)
    !immobilization by tap into the top soil layer
    IF(RINHPR.GT.0.0_r8)THEN
      CNH4X                    = AZMAX1(CNH4SU-Z4MN)
      CNH4Y                    = AZMAX1(CNH4BU-Z4MN)
      RNH4DmndLitrHeter(NGL,K) = AMIN1(RINHPR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M                      = Z4MN*VLWatMicP
      RNH4imobilLitrHeter(NGL,K) = AMIN1(AttenfNH4Heter(NGL,K)*AZMAX1((ZNH4TU-ZNH4M)),RNH4DmndLitrHeter(NGL,K))
    ELSE
      RNH4DmndLitrHeter(NGL,K)   = 0.0_r8
      RNH4imobilLitrHeter(NGL,K) = RINHPR
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+RNH4imobilLitrHeter(NGL,K)
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
    !     RNO3imobilLitrHeter=substrate-unlimited NO3 immobiln
    !     VOLW=water content
    !     ZNO3M=NO3 not available for uptake
    !     AttenfNO3Heter=fraction of biological NO3 demand
    !     RNO3imobilLitrHeter=substrate-limited NO3 immobiln
    !     NetNH4Mineralize=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
    !
    RINOPR=AZMAX1(RINHPR-RNH4imobilLitrHeter(NGL,K))
    !immobilization by tapping into the top soil layer
    IF(RINOPR.GT.0.0_r8)THEN
      CNO3X                    = AZMAX1(CNO3SU-ZOMN)
      CNO3Y                    = AZMAX1(CNO3BU-ZOMN)
      RNO3DmndLitrHeter(NGL,K) = AMAX1(RINOPR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M                      = ZOMN*VLWatMicP
      RNO3imobilLitrHeter(NGL,K) = AMIN1(AttenfNO3Heter(NGL,K)*AZMAX1((ZNO3TU-ZNO3M)),RNO3DmndLitrHeter(NGL,K))
    ELSE
      RNO3DmndLitrHeter(NGL,K)   = 0._r8
      RNO3imobilLitrHeter(NGL,K) = 0._r8
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+RNO3imobilLitrHeter(NGL,K)
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
    !     RH2PO4imobilLitrHeter=substrate-limited H2PO4 mineraln-immobiln
    !     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
    !
    RIPOPR=RIPOP-RH2PO4imobilSoilHeter(NGL,K)
    !immobilization by tapping into top soil layer
    IF(RIPOPR.GT.0.0_r8)THEN
      CH2PX                      = AZMAX1(CH2P4U-HPMN)
      CH2PY                      = AZMAX1(CH2P4BU-HPMN)
      RH2PO4DmndLitrHeter(NGL,K) = AMIN1(RIPOPR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M                        = HPMN*VOLWU
      RH2PO4imobilLitrHeter(NGL,K) = AMIN1(AttenfH2PO4Heter(NGL,K)*AZMAX1((H2P4TU-H2P4M)),RH2PO4DmndLitrHeter(NGL,K))
    ELSE
      RH2PO4DmndLitrHeter(NGL,K)   = 0.0_r8
      RH2PO4imobilLitrHeter(NGL,K) = RIPOPR
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+RH2PO4imobilLitrHeter(NGL,K)
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
    !     RH1PO4imobilLitrHeter=substrate-limited HPO4 minereraln-immobiln
    !     NetPO4Mineralize=total HPO4 net mineraln (-ve) or immobiln (+ve)
    !
    FH1PS=VLPO4
    FH1PB=VLPOB
    RIP1PR=0.1_r8*AZMAX1(RIPOPR-RH2PO4imobilLitrHeter(NGL,K))
    !immobilization
    IF(RIP1PR.GT.0.0_r8)THEN
      CH1PX                      = AZMAX1(CH1P4U-HPMN)
      CH1PY                      = AZMAX1(CH1P4BU-HPMN)
      RH1PO4DmndLitrHeter(NGL,K) = AMIN1(RIP1PR,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M                        = HPMN*VOLWU
      RH1PO4imobilLitrHeter(NGL,K) = AMIN1(AttenfH1PO4Heter(NGL,K)*AZMAX1((H1P4TU-H1P4M)),RH1PO4DmndLitrHeter(NGL,K))
    ELSE
      RH1PO4DmndLitrHeter(NGL,K)   = 0.0_r8
      RH1PO4imobilLitrHeter(NGL,K) = 0._r8
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+RH1PO4imobilLitrHeter(NGL,K)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BiomassMineralization

!------------------------------------------------------------------------------------------

  subroutine GatherHeterotrophRespiration(I,J,NGL,N,K,RMaintDefcitcitHeter,micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(out) :: RMaintDefcitcitHeter   !deficit for maintenance respiraiton, [gC d-2 h-1]
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx    
  character(len=*), parameter :: subname='GatherHeterotrophRespiration'
  integer :: MID3,MID1
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
!     begin_execution
  associate(                                                   &
    TempMaintRHeter         => nmics%TempMaintRHeter,          &
    OMN2                    => nmics%OMN2,                     &
    Resp4NFixHeter          => nmicf%Resp4NFixHeter,           &
    RespGrossHeter          => nmicf%RespGrossHeter,           &
    RMaintRespHeter         => nmicf%RMaintRespHeter,          &
    RMaintDmndHeter         => nmicf%RMaintDmndHeter,          &
    RNH4imobilSoilHeter     => nmicf%RNH4imobilSoilHeter,      &
    RNO3imobilSoilHeter     => nmicf%RNO3imobilSoilHeter,      &
    RGrowthRespHeter        => micflx%RGrowthRespHeter    ,    &         !growth respiraiton, [gC d-2 h-1]
    RN2FixHeter             => nmicf%RN2FixHeter,              &
    rNCOMC                  => micpar%rNCOMC,                  &
    rPCOMC                  => micpar%rPCOMC,                  &
    pH                      => micfor%pH,                      &
    ZEROS                   => micfor%ZEROS,                   &
    mid_HeterMixtCynoBacter => micpar%mid_HeterMixtCynoBacter, &
    mid_HeterAerobN2Fixer   => micpar%mid_HeterAerobN2Fixer,   &
    mid_HeterAnaerobN2Fixer => micpar%mid_HeterAnaerobN2Fixer, &
    CZ2GS                   => micstt%CZ2GS,                   &
    mBiomeHeter             => micstt%mBiomeHeter              &
  )
  call PrintInfo('beg '//subname)
  !     pH EFFECT ON MAINTENANCE RESPIRATION
  !
  !     FPH=pH effect on maintenance respiration
  !     RMOM=specific maintenance respiration rate
  !     TempMaintRHeter=temperature effect on maintenance respiration
  !     OMN=microbial N biomass
  !
  !
  RGrowthRespHeter(NGL,K) = AZMAX1(RespGrossHeter(NGL,K)-RMaintRespHeter(NGL,K))
  RMaintDefcitcitHeter    = AZMAX1(RMaintRespHeter(NGL,K)-RespGrossHeter(NGL,K))
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
  !     RN2FixHeter=N2 fixation rate, [gN d-2 h-1]
  !
  IF(N.EQ.mid_HeterAerobN2Fixer .OR. N.EQ.mid_HeterAnaerobN2Fixer .or. N.eq.mid_HeterMixtCynoBacter)THEN
    MID3  = micpar%get_micb_id(ibiom_reserve,NGL)
    RGN2P = AZMAX1(mBiomeHeter(ielmc,MID3,K)*rNCOMC(ibiom_reserve,NGL,K)-mBiomeHeter(ielmn,MID3,K))/EN2F(N)
    IF(RGrowthRespHeter(NGL,K).GT.ZEROS)THEN
      Resp4NFixHeter(NGL,K)=AMIN1(RGrowthRespHeter(NGL,K)*RGN2P/(RGrowthRespHeter(NGL,K)+RGN2P) &
        *CZ2GS/(CZ2GS+ZFKM),OMGR*mBiomeHeter(ielmc,MID3,K))
      RN2FixHeter(NGL,K)=Resp4NFixHeter(NGL,K)*EN2F(N)  
    ELSE
      Resp4NFixHeter(NGL,K) = 0.0_r8
      RN2FixHeter(NGL,K)    = 0._r8
    ENDIF    
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherHeterotrophRespiration
!------------------------------------------------------------------------------------------

  subroutine GatherHetertrophAnabolicFlux(I,J,NGL,N,K,RMaintDefcitcitHeter,&
    SPOMK,micfor,micstt,nmicf,nmics,ncplxf,ncplxs,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: RMaintDefcitcitHeter  !maintenance deficit, [gC d-2 h-1]
  real(r8), intent(in) :: SPOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx  
  character(len=*), parameter :: subname='GatherHetertrophAnabolicFlux'

  integer :: M,MID1,MID3,MID,NE
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CCC,CXC,RCCE(NumPlantChemElms)  

  real(r8) :: RCCC,RCCN,RCCP  
  real(r8) :: CGOMX,CGOMD,AGOMD
  real(r8) :: CGOXC              !DOC+acetate uptake flux, [gC d-2 h-1]
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
    RMetabDOCUptkHeter               => nmicf%RMetabDOCUptkHeter,               &
    RMaintRespHeter                  => nmicf%RMaintRespHeter,                  &    
    FGOCP                            => nmicf%FGOCP,                            & !fraction of C uptake as DOC
    FGOAP                            => nmicf%FGOAP,                            & !fraction of C uptake as acetate
    RMetabAcetUptkHeter              => nmicf%RMetabAcetUptkHeter,              &
    ECHZHeter                        => nmicf%ECHZHeter,                        &
    NonstX2stBiomHeter               => nmicf%NonstX2stBiomHeter,               &
    RespGrossHeter                   => nmicf%RespGrossHeter,                   &
    RCO2FixCyano                     => nmicf%RCO2FixCyano,                     &
    RGrowthRespHeter                 => micflx%RGrowthRespHeter        ,        &
    RNOxDOCReduxRespDenitLim         => nmicf%RNOxDOCReduxRespDenitLim,         & !DOC contributed respiraiton by denitrifcation
    RNOxAcetReduxRespDenitLim        => nmicf%RNOxAcetReduxRespDenitLim,        & !acetate contributed respiraiton by denitrifcation   
    RMaintDmndHeter                  => nmicf%RMaintDmndHeter,                  &
    fPhotoR                          => nmicf%fPhotoR,                          &    
    RkillLitfalOMHeter               => nmicf%RkillLitfalOMHeter,               &
    RkillLitrfal2HumOMHeter          => nmicf%RkillLitrfal2HumOMHeter,          &
    RkillLitrfal2ResduOMHeter        => nmicf%RkillLitrfal2ResduOMHeter,        &
    RCCMEHeter                       => nmicf%RCCMEHeter,                       &
    RMaintDefLitrfal2HumOMHeter      => nmicf%RMaintDefLitrfal2HumOMHeter,      &
    RMaintDefLitrfal2ResduOMHeter    => nmicf%RMaintDefLitrfal2ResduOMHeter,    &
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
  call PrintInfo('beg '//subname)
  !     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
  !     FROM O2, NOX AND C REDUCTION
  !
  !     CGOMX=DOC+acetate uptake from aerobic growth respiration
  !     CGOMD=DOC  uptake from denitrifier growth respiration
  !     RMaintRespHeter=maintenance respiration
  !     RespGrossHeter=total respiration
  !     RNOxDOCReduxRespDenitLim=respiration generated by denitrifcation
  !     Resp4NFixHeter=respiration for N2 fixation
  !     ECHZHeter,ENOX=growth respiration efficiencies for O2, NOx reduction
  !     DOMuptk4GrothHeter,RMetabDOCUptkHeter,RMetabAcetUptkHeter=total DOC+acetate, DOC, acetate uptake(heterotrophs
  !     DOMuptk4GrothHeter=total CO2,CH4 uptake (autotrophs)
  !     DOMuptk4GrothHeter,DOMuptk4GrothHeter=DON, DOP uptake
  !     FracHeterBiomOfActK=faction of OMActHeterin total OMA
  !     rCNDOM,rCPDOM=DON/DOC, DOP/DOC
  !     FCN,FCP=limitation from N,P
  ! gross respiration equals to maintenance+respiraiton for N-fixation + growth respiraiton
  if(N.eq.micpar%mid_HeterMixtCynoBacter)THEN
    IF(ECHZHeter(NGL,K).GT.0._r8)then
    CGOMX     = AMIN1(RMaintRespHeter(NGL,K),RespGrossHeter(NGL,K))+Resp4NFixHeter(NGL,K)+(RGrowthRespHeter(NGL,K)-Resp4NFixHeter(NGL,K))*&
      (fPhotoR(NGL,K)+(1._r8-fPhotoR(NGL,K))/ECHZHeter(NGL,K))      
    ELSE
      CGOMX     = AMIN1(RMaintRespHeter(NGL,K),RespGrossHeter(NGL,K))+Resp4NFixHeter(NGL,K)+(RGrowthRespHeter(NGL,K)-Resp4NFixHeter(NGL,K))*fPhotoR(NGL,K)          
    ENDIF
    CGOMX = CGOMX-RCO2FixCyano(NGL,K)      
  ELSE
    CGOMX = AMIN1(RMaintRespHeter(NGL,K),RespGrossHeter(NGL,K))+Resp4NFixHeter(NGL,K)+(RGrowthRespHeter(NGL,K)-Resp4NFixHeter(NGL,K))/ECHZHeter(NGL,K)    
  endif

  CGOMD     = RNOxDOCReduxRespDenitLim(NGL,K)/ENOX
  AGOMD     = RNOxAcetReduxRespDenitLim(NGL,K)/ENOX
  CDOMuptk1 = CDOMuptk1+CGOMX !DOC used for growth
  CDOMuptk2 = CDOMuptk2+CGOMD !DOC used for denitrifcation
  tROMT     = tROMT+RMaintRespHeter(NGL,K)
  tGROMO    = tGROMO+RespGrossHeter(NGL,K)

  DOMuptk4GrothHeter(ielmc,NGL,K) = CGOMX+CGOMD
  if(N.eq.micpar%mid_HeterMixtCynoBacter .and. RGrowthRespHeter(NGL,K).GT.0._r8)then
    DOMuptk4GrothHeter(ielmc,NGL,K) = DOMuptk4GrothHeter(ielmc,NGL,K)+RCO2FixCyano(NGL,K)
  endif
  RMetabDOCUptkHeter(NGL,K)       = CGOMX*FGOCP(NGL,K)+CGOMD         !include DOC for respiraiton+denitrifcation
  RMetabAcetUptkHeter(NGL,K)      = CGOMX*FGOAP(NGL,K)+AGOMD         !acetate uptake for metabolism
  CGOXC                           = RMetabDOCUptkHeter(NGL,K)+RMetabAcetUptkHeter(NGL,K)

  !obtain organic nutrient uptake
  DOMuptk4GrothHeter(ielmn,NGL,K)=AZMAX1(AMIN1(DOM(idom_don,K)*FracHeterBiomOfActK(NGL,K),CGOXC*rCNDOM(K)/FCN(NGL,K)))
  DOMuptk4GrothHeter(ielmp,NGL,K)=AZMAX1(AMIN1(DOM(idom_dop,K)*FracHeterBiomOfActK(NGL,K),CGOXC*rCPDOM(K)/FCP(NGL,K)))

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
  !
  !     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
  !
  !     CGOMZ=transfer from nonstructural to structural microbial C
  !     GrowthEnvScalHeter=temperature+water stress function
  !     OMGR=rate constant for transferring nonstructural to structural C
  !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
  !     FL=partitioning between labile and resistant microbial components
  !     OMC,OMN,OMP=nonstructural microbial C,N,P
  !
  MID3  = micpar%get_micb_id(ibiom_reserve,NGL)
  CGOMZ = GrowthEnvScalHeter(NGL,K)*OMGR*AZMAX1(mBiomeHeter(ielmc,MID3,K))
  
  ! M=1:labile, 2, recalcitrant
  D745: DO M=1,2
    !reserve flux to kinetic and structural biomass
    NonstX2stBiomHeter(ielmc,M,NGL,K)=FL(M)*CGOMZ
    IF(mBiomeHeter(ielmc,MID3,K).GT.ZEROS)THEN
      Do NE=2,NumPlantChemElms
        NonstX2stBiomHeter(NE,M,NGL,K)=AZMAX1(mBiomeHeter(NE,MID3,K))*AMIN1(FL(M), &
          NonstX2stBiomHeter(ielmc,M,NGL,K)/mBiomeHeter(ielmc,MID3,K))
      ENDDO  
    ELSE
      NonstX2stBiomHeter(2:NumPlantChemElms,M,NGL,K)=0.0_r8
    ENDIF
    !
    !     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
    !     RATE, TEMPERATURE
    !
    !     SPOMX=rate constant for microbial decomposition
    !     SPOMC=basal decomposition/mortality rate
    !     SPOMK=effect of low microbial C concentration on microbial decay
    !     RXOMC,RXOMN,RKillOMHeter=microbial C,N,P decomposition
    !     RkillLitfalOMHeter,RDOMN,RDOMP=microbial C,N,P LitrFall
    !     R3OMC,R3OMN,RkillRecycOMHeter=microbial C,N,P recycling
    !

    MID1=micpar%get_micb_id(ibiom_kinetic,NGL)

    IF(mBiomeHeter(ielmc,MID3,K).GT.ZEROS .AND. mBiomeHeter(ielmc,MID1,K).GT.ZEROS)THEN
      CCC=AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeHeter(ielmn,MID3,K)/(mBiomeHeter(ielmn,MID3,K)+mBiomeHeter(ielmc,MID3,K)*rNCOMC(ibiom_reserve,NGL,K)) &
        ,mBiomeHeter(ielmp,MID3,K)/(mBiomeHeter(ielmp,MID3,K)+mBiomeHeter(ielmc,MID3,K)*rPCOMC(ibiom_reserve,NGL,K))))
      CXC  = mBiomeHeter(ielmc,MID3,K)/mBiomeHeter(ielmc,MID1,K)
      C3C  = 1.0_r8/(1.0_r8+CXC/CKC)

      CNC  = AZMAX1(AMIN1(1.0_r8,mBiomeHeter(ielmc,MID3,K)/(mBiomeHeter(ielmc,MID3,K)+mBiomeHeter(ielmn,MID3,K)/rNCOMC(ibiom_reserve,NGL,K))))
      CPC  = AZMAX1(AMIN1(1.0_r8,mBiomeHeter(ielmc,MID3,K)/(mBiomeHeter(ielmc,MID3,K)+mBiomeHeter(ielmp,MID3,K)/rPCOMC(ibiom_reserve,NGL,K))))
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

    MID   = micpar%get_micb_id(M,NGL)
    SPOMX = SQRT(GrowthEnvScalHeter(NGL,K))*SPOMC(M)*SPOMK(M)
    DO NE=1,NumPlantChemElms
      RKillOMHeter(NE,M,NGL,K)=AZMAX1(mBiomeHeter(NE,MID,K)*SPOMX)

      RkillRecycOMHeter(NE,M,NGL,K)= RKillOMHeter(NE,M,NGL,K)*RCCE(NE)

      RkillLitfalOMHeter(NE,M,NGL,K)=AZMAX1(RKillOMHeter(NE,M,NGL,K)-RkillRecycOMHeter(NE,M,NGL,K))
      !
      !     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
      !     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
      !
      !     RHOMC,RHOMN,RkillLitrfal2HumOMHeter=transfer of microbial C,N,P LitrFall to humus
      !     EHUM=humus transfer fraction from hour1.f
      !     RkillLitrfal2ResduOMHeter,RkillLitrfal2ResduOMHeter,RCOMP=transfer of microbial C,N,P LitrFall to residue
      !
      RkillLitrfal2HumOMHeter(NE,M,NGL,K)=RkillLitfalOMHeter(NE,M,NGL,K)*EHUM
      !
      !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
      !
      RkillLitrfal2ResduOMHeter(NE,M,NGL,K)=RkillLitfalOMHeter(NE,M,NGL,K)-RkillLitrfal2HumOMHeter(NE,M,NGL,K)

    ENDDO
  ENDDO D745

  !=========================================================================================================
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

  IF(RMaintDefcitcitHeter.GT.ZEROS .AND. RMaintRespHeter(NGL,K).GT.ZEROS .AND. RCCC.GT.ZERO)THEN
    FracMaintDeficit=RMaintDefcitcitHeter/RMaintRespHeter(NGL,K)
    D730: DO M=1,2
      MID                                    = micpar%get_micb_id(M,NGL)
      RMaintDefcitKillOMHeter(ielmc,M,NGL,K) = AMIN1(mBiomeHeter(ielmc,MID,K),AZMAX1(FracMaintDeficit*RMaintDmndHeter(M,NGL,K)/RCCC))
      RMaintDefcitKillOMHeter(ielmn,M,NGL,K) = AMIN1(mBiomeHeter(ielmn,MID,K),AZMAX1(RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*rCNBiomeActHeter(ielmn,NGL,K)))
      RMaintDefcitKillOMHeter(ielmp,M,NGL,K) = AMIN1(mBiomeHeter(ielmp,MID,K),AZMAX1(RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*rCNBiomeActHeter(ielmp,NGL,K)))
        
      DO NE=1,NumPlantChemElms
        RMaintDefcitRecycOMHeter(NE,M,NGL,K)   = RMaintDefcitKillOMHeter(NE,M,NGL,K)*RCCE(NE)
        RMaintDefcitLitrfalOMHeter(NE,M,NGL,K) = AZMAX1(RMaintDefcitKillOMHeter(NE,M,NGL,K)-RMaintDefcitRecycOMHeter(NE,M,NGL,K))
        !
        ! HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
        ! PRODUCTS
        !
        ! RMaintDefLitrfal2HumOMHeter,RHMMN,RMaintDefLitrfal2HumOMHeter=transfer of senesence LitrFall C,N,P to humus
        ! EHUM=humus transfer fraction
        ! RMaintDefLitrfal2ResduOMHeter,RCMMN,RMaintDefLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
        !

        RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K)   = RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)*EHUM
        RMaintDefLitrfal2ResduOMHeter(NE,M,NGL,K) = RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)-RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K)
        RCCMEheter(NE,M,K)                        = RCCMEheter(NE,M,K)+RkillLitrfal2ResduOMHeter(NE,M,NGL,K)+RMaintDefLitrfal2ResduOMHeter(NE,M,NGL,K)

      ENDDO
    ENDDO D730
  ELSE
    D720: DO M=1,2
      DO NE=1,NumPlantChemElms
        RMaintDefcitKillOMHeter(NE,M,NGL,K)       = 0.0_r8
        RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)    = 0.0_r8
        RMaintDefcitRecycOMHeter(NE,M,NGL,K)      = 0.0_r8
        RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K)   = 0.0_r8
        RMaintDefLitrfal2ResduOMHeter(NE,M,NGL,K) = 0.0_r8
      ENDDO

    ENDDO D720
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherHetertrophAnabolicFlux

end module MicBGCMod
