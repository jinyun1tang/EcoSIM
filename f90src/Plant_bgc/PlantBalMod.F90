module PlantBalMod
!
!Description
!code to do mass balance calculation for plant bgc

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only : pltpar
  use PlantMathFuncMod
  use PlantAPIData  
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: SumPlantBiom
  public :: SumPlantBiomStates
  public :: ZeroGrosub  
  public :: SumRootBiome
  public :: SumPlantBranchBiome
  public :: EnterPlantBalance
  public :: ExitPlantBalance
  public :: SumPlantRootGas
  logical,save  :: lfile(2)=.true.
  contains

!------------------------------------------------------------------------------------------

  subroutine SumPlantBiom(I,J,NZ,header)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), intent(in) :: header
  integer :: L,K,N,NE,NB,M
  real(r8) :: root1st,root2nd

!     begin_execution
  associate(                                                        &
    NU                       => plt_site%NU,                        &
    MaxNodesPerBranch1       => pltpar%MaxNodesPerBranch1,          &
    iPlantPhotosynthesisType => plt_photo%iPlantPhotosynthesisType, &
    NH3Dep2Can_pft           => plt_bgcr%NH3Dep2Can_pft,            &
    NH3Dep2Can_brch          => plt_rbgc%NH3Dep2Can_brch,           &    !
    GrossResp_pft            => plt_bgcr%GrossResp_pft,             &
    CanopyGrosRCO2_pft       => plt_bgcr%CanopyGrosRCO2_pft,        &
    RootGrosRCO2_pft         => plt_bgcr%RootGrosRCO2_pft,          &
    GrossCO2Fix_pft          => plt_bgcr%GrossCO2Fix_pft,           &
    NumOfBranches_pft        => plt_morph%NumOfBranches_pft,        &
    NodulInfectElms_pft      => plt_bgcr%NodulInfectElms_pft,       &
    LitrfalStrutElms_pft     => plt_bgcr%LitrfalStrutElms_pft,      &
    LitrfalStrutElms_pvr     => plt_bgcr%LitrfalStrutElms_pvr,      &
    RootMycoExudElms_pft     => plt_rbgc%RootMycoExudElms_pft,      &
    MY                       => plt_morph%MY,                       &
    MaxSoiL4Root_pft         => plt_morph%MaxSoiL4Root_pft,         &
    NumRootAxes_pft          => plt_morph%NumRootAxes_pft,          &
    MaxNumRootLays           => plt_site%MaxNumRootLays,            &
    RootElms_pft             => plt_biom%RootElms_pft,              &
    ShootC4NonstC_brch       => plt_biom%ShootC4NonstC_brch,        &
    RootCO2Autor_pvr         => plt_rbgc%RootCO2Autor_pvr,          &
    ShootElms_pft            => plt_biom%ShootElms_pft              &
  )

  CALL SumPlantBiomStates(I,J,NZ,header)
  
  if(RootElms_pft(ielmc,NZ)>1.e16 .or. ShootElms_pft(ielmc,NZ)>1.e16)then
    write(*,*)'rootshoot',RootElms_pft(ielmc,NZ),ShootElms_pft(ielmc,NZ)
    stop
  endif
  !sum fluxes
!     NH3Dep2Can_brch,NH3Dep2Can_pft=PFT NH3 flux between atmosphere and branch,canopy
  NH3Dep2Can_pft(NZ)=0._r8
  DO NB=1,NumOfBranches_pft(NZ)
    NH3Dep2Can_pft(NZ)=NH3Dep2Can_pft(NZ)+NH3Dep2Can_brch(NB,NZ)
  ENDDO    

  LitrfalStrutElms_pft(1:NumPlantChemElms,NZ)=0._r8
  DO L=0,MaxNumRootLays
    DO K=1,pltpar%NumOfPlantLitrCmplxs      
      DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_pft(NE,NZ)=LitrfalStrutElms_pft(NE,NZ)+LitrfalStrutElms_pvr(NE,M,K,L,NZ)
        ENDDO
      ENDDO
    ENDDO      
  ENDDO

  RootGrosRCO2_pft(NZ)=0._r8
  DO N=1,MY(NZ)
    DO L=NU,MaxSoiL4Root_pft(NZ)
      RootGrosRCO2_pft(NZ)=RootGrosRCO2_pft(NZ)+RootCO2Autor_pvr(N,L,NZ)
    ENDDO     
  ENDDO  

  GrossResp_pft(NZ)=RootGrosRCO2_pft(NZ)+CanopyGrosRCO2_pft(NZ)

  end associate
  end subroutine SumPlantBiom

!------------------------------------------------------------------------------------------

  subroutine SumPlantRootGas(I,J)

  implicit none
  integer, intent(in) :: I,J
  integer :: NZ ,L,N,idg,K
  real(r8) :: trcg(idg_beg:idg_NH3)
  associate(                                                       &
    NU                       => plt_site%NU,                       &  
    trcg_root_vr             => plt_rbgc%trcg_root_vr,             &
    trcg_rootml_pvr          => plt_rbgc%trcg_rootml_pvr,          &
    trcs_rootml_pvr          => plt_rbgc%trcs_rootml_pvr,          &
    MY                       => plt_morph%MY,                      &
    MaxSoiL4Root_pft         => plt_morph%MaxSoiL4Root_pft         &
  )
  trcg_root_vr(idg_beg:idg_NH3,:)         = 0._r8

  trcg(:)=0._r8
  DO NZ=1,plt_site%NP  
    IF(.not.plt_pheno%IsPlantActive_pft(NZ).EQ.iActive)cycle    
    DO L=NU,MaxSoiL4Root_pft(NZ)
      DO N=1,MY(NZ)  
        DO idg=idg_beg,idg_NH3
          trcg_root_vr(idg,L)=trcg_root_vr(idg,L)+trcs_rootml_pvr(idg,N,L,NZ)+trcg_rootml_pvr(idg,N,L,NZ)        
        ENDDO        
      ENDDO  
      DO idg=idg_beg,idg_NH3
        trcg(idg)=trcg(idg)+trcg_root_vr(idg,L)
      ENDDO
    ENDDO
!    if(I==140 .and. J<=2)write(116,*)I*1000+J,MaxSoiL4Root_pft(NZ),trcg(idg_CH4)
  ENDDO
!  if(I==140 .and. J>=20)write(116,*)trcg(idg_N2),'N2'
  end associate
  end subroutine SumPlantRootGas
!------------------------------------------------------------------------------------------

  subroutine SumPlantBiomStates(I,J,NZ,header)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  character(len=*),intent(in) :: header
  integer :: NB,NE,K,jk
  real(r8) :: root1st,root2nd
  real(r8) :: balc

!     begin_execution
  associate(                                                         &
    NU                        => plt_site%NU,                        &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,       &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,        &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,         &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,          &
    MaxNumRootLays            => plt_site%MaxNumRootLays,            &
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,     &
    RootElmsBeg_pft           => plt_biom%RootElmsBeg_pft,           &
    ShootElmsBeg_pft          => plt_biom%ShootElmsBeg_pft,          &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr,   &
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr,   &
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch, &
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch, &
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft,    &
    StandDeadKCompElms_pft    => plt_biom%StandDeadKCompElms_pft,    &
    ShootStrutElms_brch       => plt_biom%ShootStrutElms_brch,       &
    RootElms_pft              => plt_biom%RootElms_pft,              &
    ShootElms_pft             => plt_biom%ShootElms_pft              &
  )
  
  !shoots
  DO NB=1,NumOfBranches_pft(NZ)
    CALL SumPlantBranchBiome(NB,NZ)
  ENDDO

  DO NE=1,NumPlantChemElms
    ShootElms_pft(NE,NZ)=sum(ShootStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))    
    if(ShootElms_pft(NE,NZ)>1.e20)then
      write(*,*)'sumplantstates',NE,NumOfBranches_pft(NZ),ShootStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)
      stop
    endif
  ENDDO

  !add nodule, currently, a plant either has canopy or root N-fixing symbiosis, not both
  IF(is_plant_N2fix(iPlantNfixType_pft(NZ)))THEN
    IF(is_canopy_N2fix(iPlantNfixType_pft(NZ)))THEN
      DO NE=1,NumPlantChemElms
        ShootElms_pft(NE,NZ)=ShootElms_pft(NE,NZ)  &
          +sum(CanopyNodulStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)) &
          +sum(CanopyNodulNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
        if(ShootElms_pft(NE,NZ)>1.e20)then
          write(*,*)'sum2plantstates'
          stop
        endif
      ENDDO
    ENDIF
  ENDIF
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!add standing dead to shoot
  DO NE=1,NumPlantChemElms
    StandDeadStrutElms_pft(NE,NZ)=sum(StandDeadKCompElms_pft(NE,1:jsken,NZ))    
    if(StandDeadStrutElms_pft(NE,NZ)>1.e20)then
      write(*,*)'sum3plantstates',StandDeadStrutElms_pft(NE,NZ)
      DO jk=1,jsken
        write(*,*)NE,StandDeadKCompElms_pft(NE,jk,NZ)
      enddo
      stop
    endif
  ENDDO

  call SumRootBiome(NZ,RootElms_pft(:,NZ))

  return
  if(I>=154)then
  if(NZ==1)then
    if(lfile(NZ))then
      write(243,'(A14,X,A13,6(X,A16))')'header','doy','rootC','shootC','rootN','shootN','rootP','shootP'
    endif
    write(243,'(A14,X,F13.6,6(X,F16.6))')trim(header),I+J/24.,&
      (RootElms_pft(NE,NZ),ShootElms_pft(NE,NZ),NE=1,NumPlantChemElms)
    write(243,'(A14,X,F13.6,6(X,F16.6))')'x',I+J/24.,&
      (RootElmsbeg_pft(NE,NZ),ShootElmsbeg_pft(NE,NZ),NE=1,NumPlantChemElms)      
  else
    if(lfile(NZ))then
      write(244,'(A14,X,A13,6(X,A14))')'header','doy','rootC','shootC','rootN','shootN','rootP','shootP'
    endif  
    write(244,'(A14,X,F13.6,6(X,F14.6))')trim(header),I+J/24.,&
      (RootElms_pft(NE,NZ),ShootElms_pft(NE,NZ),NE=1,NumPlantChemElms)
    write(244,'(A14,X,F13.6,6(X,F14.6))')'x',I+J/24.,&
      (RootElmsbeg_pft(NE,NZ),ShootElmsbeg_pft(NE,NZ),NE=1,NumPlantChemElms)            
  endif
  lfile(NZ)=.false.
  endif
  end associate
  END subroutine SumPlantBiomStates

!------------------------------------------------------------------------------------------

  subroutine SumPlantBranchBiome(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  integer :: NE,K
  associate(                                                            &
    MaxNodesPerBranch1         => pltpar%MaxNodesPerBranch1,            &
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft,          &
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch,          &
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch,        &
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch,         &
    CanopyStalkC_pft           => plt_biom%CanopyStalkC_pft,            &
    StalkLiveBiomassC_brch     => plt_biom%StalkLiveBiomassC_brch,      &
    CanopyLeafShethC_pft       => plt_biom%CanopyLeafShethC_pft,        &
    LeafPetolBiomassC_brch     => plt_biom%LeafPetolBiomassC_brch,      &
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch,          &
    HuskStrutElms_brch         => plt_biom%HuskStrutElms_brch,          &
    EarStrutElms_brch          => plt_biom%EarStrutElms_brch,           &
    GrainStrutElms_brch        => plt_biom%GrainStrutElms_brch,         &
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch,        &
    ShootC4NonstC_brch         => plt_biom%ShootC4NonstC_brch,          &
    CPOOL3_node                => plt_photo%CPOOL3_node,                &
    CPOOL4_node                => plt_photo%CPOOL4_node,                &
    iPlantPhotosynthesisType   => plt_photo%iPlantPhotosynthesisType,   &
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node, &
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node,  &
    ShootStrutElms_brch        => plt_biom%ShootStrutElms_brch          &
  )  
  DO NE=1,NumPlantChemElms
    ShootStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ) &
      +PetoleStrutElms_brch(NE,NB,NZ)+StalkStrutElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ) &
      +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ) &
      +CanopyNonstElms_brch(NE,NB,NZ)
    
  ENDDO
!  CanopyStalkC_pft(NZ)     = sum(StalkLiveBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafShethC_pft(NZ) = sum(LeafPetolBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))

  !add C4 specific reserve carbon
  ShootC4NonstC_brch(NB,NZ)=0._r8
  IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN      
    D3251: DO K=1,MaxNodesPerBranch1
      ShootC4NonstC_brch(NB,NZ)=ShootC4NonstC_brch(NB,NZ)+CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D3251
    ShootStrutElms_brch(ielmc,NB,NZ)=ShootStrutElms_brch(ielmc,NB,NZ)+ShootC4NonstC_brch(NB,NZ)
  ENDIF

  end associate
  end subroutine SumPlantBranchBiome

!------------------------------------------------------------------------------------------

  subroutine ZeroGrosub()
  
  implicit none
  integer :: NZ,K,L,M,NE

  associate(                                                 &
    NP0                   => plt_site%NP0,                   &
    MY                    => plt_morph%MY,                   &
    NU                    => plt_site%NU,                    &
    MaxSoiL4Root_pft      => plt_morph%MaxSoiL4Root_pft,     &
    MaxNumRootLays        => plt_site%MaxNumRootLays,        &
    CO2NetFix_pft         => plt_bgcr%CO2NetFix_pft,         &
    CanopyRespC_CumYr_pft => plt_bgcr%CanopyRespC_CumYr_pft, &
    LitrfalStrutElms_pft  => plt_bgcr%LitrfalStrutElms_pft,  &
    NodulInfectElms_pft   => plt_bgcr%NodulInfectElms_pft,   &
    RootMycoExudElms_pft  => plt_rbgc%RootMycoExudElms_pft,  &
    NH3Dep2Can_pft        => plt_bgcr%NH3Dep2Can_pft,        &
    GrossResp_pft         => plt_bgcr%GrossResp_pft,         &
    GrossCO2Fix_pft       => plt_bgcr%GrossCO2Fix_pft,       &
    CO2FixCL_pft          => plt_rbgc%CO2FixCL_pft,          &
    CO2FixLL_pft          => plt_rbgc%CO2FixLL_pft,          &
    CanopyGrosRCO2_pft    => plt_bgcr%CanopyGrosRCO2_pft,    &
    RootCO2Autor_pvr      => plt_rbgc%RootCO2Autor_pvr,      &
    LitrfalStrutElms_pvr  => plt_bgcr%LitrfalStrutElms_pvr   &
  )  
  
  plt_rbgc%trcs_plant_uptake_vr=0._r8

  D9980: DO NZ=1,NP0
    D1: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pvr(NE,M,K,L,NZ)=0._r8
          ENDDO
        ENDDO
      ENDDO
    ENDDO D1
    RootCO2Autor_pvr(1:MY(NZ),NU:MaxSoiL4Root_pft(NZ),NZ) = 0._r8
    NH3Dep2Can_pft(NZ)                                    = 0._r8
    GrossResp_pft(NZ)                                     = 0._r8
    GrossCO2Fix_pft(NZ)                                   = 0._r8
    CO2NetFix_pft(NZ)                                     = 0._r8
    CanopyGrosRCO2_pft(NZ)                                = 0._r8
    RootMycoExudElms_pft(1:NumPlantChemElms,NZ)           = 0._r8
    NodulInfectElms_pft(1:NumPlantChemElms,NZ)            = 0._r8
    CO2FixCL_pft(NZ)                                      = 0._r8
    CO2FixLL_pft(NZ)                                      = 0._r8

  ENDDO D9980
  end associate
  end subroutine ZeroGrosub

!------------------------------------------------------------------------------------------

  subroutine SumRootBiome(NZ,mass_roots,massr1st,massr2nd,massnonst,massnodul)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(out) :: mass_roots(NumPlantChemElms)
  real(r8), optional, intent(out) :: massr1st(NumPlantChemElms)
  real(r8), optional, intent(out) :: massr2nd(NumPlantChemElms)
  real(r8), optional, intent(out) :: massnonst(NumPlantChemElms)
  real(r8), optional, intent(out) :: massnodul(NumPlantChemElms)
  integer :: NE,N,L
  real(r8) :: massr1st1(NumPlantChemElms)
  real(r8) :: massr2nd1(NumPlantChemElms)
  real(r8) :: massnodul1(NumPlantChemElms)
  real(r8) :: massnonst1(NumPlantChemElms)
  associate(                                                         &
    NU                        => plt_site%NU,                        &
    MY                        => plt_morph%MY,                       &
    MaxNumRootLays            => plt_site%MaxNumRootLays,            &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,         &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,       &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,          &
    RootNodulStrutElms_rpvr    => plt_biom%RootNodulStrutElms_rpvr,    &
    RootNodulNonstElms_rpvr    => plt_biom%RootNodulNonstElms_rpvr,    &
    RootMycoNonstElms_pft     => plt_biom%RootMycoNonstElms_pft,     &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    RootStrutElms_pft         => plt_biom%RootStrutElms_pft,         &
    RootMassElm_pvr           => plt_biom%RootMassElm_pvr ,          &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     &
  )
  massr1st1=0._r8;massr2nd1=0._r8;massnodul1=0._r8
  mass_roots=0._r8

  DO NE=1,NumPlantChemElms
    DO L=NU,MaxNumRootLays
      RootMassElm_pvr(NE,L,NZ)= sum(RootMyco1stStrutElms_rpvr(NE,1:MY(NZ),L,1:NumRootAxes_pft(NZ),NZ)) + &
        sum(RootMyco2ndStrutElms_rpvr(NE,1:MY(NZ),L,1:NumRootAxes_pft(NZ),NZ)) + &
        sum(RootMycoNonstElms_rpvr(NE,1:MY(NZ),L,NZ))
    ENDDO
    massr1st1(NE)=sum(RootMyco1stStrutElms_rpvr(NE,1:MY(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    massr2nd1(NE)=sum(RootMyco2ndStrutElms_rpvr(NE,1:MY(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    RootStrutElms_pft(NE,NZ)=massr1st1(NE)+massr2nd1(NE)
    DO N=1,MY(NZ)
      RootMycoNonstElms_pft(NE,N,NZ)=sum(RootMycoNonstElms_rpvr(NE,N,NU:MaxNumRootLays,NZ))
    enddo
    massnonst1(NE)=sum(RootMycoNonstElms_pft(NE,1:MY(NZ),NZ))
    mass_roots(NE)=massr1st1(NE)+massr2nd1(NE)+massnonst1(NE)
    !add reserve to struct
    if(is_plant_N2fix(iPlantNfixType_pft(NZ)) .and. is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
      massnodul1(NE)=sum(RootNodulStrutElms_rpvr(NE,NU:MaxSoiL4Root_pft(NZ),NZ))+&
          sum(RootNodulNonstElms_rpvr(NE,NU:MaxSoiL4Root_pft(NZ),NZ))
      mass_roots(NE)=mass_roots(NE)+massnodul1(NE)    
    endif      
  ENDDO
  if(present(massr1st))massr1st=massr1st1
  if(present(massr2nd))massr2nd=massr2nd1
  if(present(massnonst))massnonst=massnonst1
  if(present(massnodul))massnodul=massnodul1
  end associate
  end subroutine SumRootBiome

!------------------------------------------------------------------------------------------
  subroutine EnterPlantBalance(I,J,NP)
  !
  implicit none
  INTEGER, INTENT(IN) :: I,J
  integer, intent(in) :: NP
  INTEGER :: NZ

  DO NZ=1,NP
    plt_biom%RootElmsBeg_pft(:,NZ)           = plt_biom%RootElms_pft(:,NZ)
    plt_biom%StandDeadStrutElmsBeg_pft(:,NZ) = plt_biom%StandDeadStrutElms_pft(:,NZ)
    plt_biom%ShootElmsBeg_pft(:,NZ)          = plt_biom%ShootElms_pft(:,NZ)
    plt_biom%SeasonalNonstElmsbeg_pft(:,NZ)  = plt_biom%SeasonalNonstElms_pft(:,NZ)
  ENDDO
  end subroutine EnterPlantBalance
!------------------------------------------------------------------------------------------

  subroutine ExitPlantBalance(I,J,NP)
  !
  !Description:
  ! mass balance check for plants
  ! dplantE/dt=input-output
  !
  ! plantE = root+shoot+standing dead
  ! C input: photosynthesis, nodule infection, DOM uptake/exudation
  ! C output: respiration, harvest
  
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NP
  integer :: NZ
  real(r8) :: balE(1:NP)

  associate(                                                         &
    GrossCO2Fix_pft           => plt_bgcr%GrossCO2Fix_pft,           &  !>0 onto plant
    GrossResp_pft             => plt_bgcr%GrossResp_pft,             &  !<0 onto plant
    NodulInfectElms_pft       => plt_bgcr%NodulInfectElms_pft,       &  !>0 onto plant
    LitrfalStrutElms_pft      => plt_bgcr%LitrfalStrutElms_pft,      &  !>0 off plant
    RootMycoExudElms_pft      => plt_rbgc%RootMycoExudElms_pft,      &  !>0 onto plant
    EcoHavstElmnt_CumYr_pft   => plt_distb%EcoHavstElmnt_CumYr_pft,  &  !>0 off plant
    RootElmsBeg_pft           => plt_biom%RootElmsBeg_pft,           &
    RootElms_pft              => plt_biom%RootElms_pft,              &
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,     &
    SeasonalNonstElmsbeg_pft  => plt_biom%SeasonalNonstElmsbeg_pft,  &
    StandDeadStrutElmsBeg_pft => plt_biom%StandDeadStrutElmsBeg_pft, &
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft,    &
    ShootElmsBeg_pft          => plt_biom%ShootElmsBeg_pft,          &
    ShootElms_pft             => plt_biom%ShootElms_pft              &
  )

  DO NZ=1,NP
    balE(NZ)=RootElms_pft(ielmc,NZ)-RootElmsBeg_pft(ielmc,NZ)  &
      +ShootElms_pft(ielmc,NZ)-ShootElmsBeg_pft(ielmc,NZ)   &
      +StandDeadStrutElms_pft(ielmc,NZ)-StandDeadStrutElmsBeg_pft(ielmc,NZ)  &
      +SeasonalNonstElms_pft(ielmc,NZ)-SeasonalNonstElmsbeg_pft(ielmc,NZ) &
      -GrossCO2Fix_pft(NZ)-GrossResp_pft(NZ)   &
      -NodulInfectElms_pft(ielmc,NZ)  &      
      -RootMycoExudElms_pft(ielmc,NZ) &
      +LitrfalStrutElms_pft(ielmc,NZ) 
!    if(plt_rad%RadPARbyCanopy_pft(NZ)>0._r8)then
!      write(123,*)I+J/24.,GrossCO2Fix_pft(NZ)/plt_rad%RadPARbyCanopy_pft(NZ),plt_rad%RadPARbyCanopy_pft(NZ)
!    endif
  ENDDO
  end associate

  end subroutine ExitPlantBalance

end module PlantBalMod