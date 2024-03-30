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
  logical,save  :: lfile(2)=.true.
  contains

!------------------------------------------------------------------------------------------

  subroutine SumPlantBiom(I,J,NZ,header)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), intent(in) :: header
  integer :: L,K,N,NE,NB,M
  real(r8) :: root1st,root2nd
  real(r8) :: balc

!     begin_execution
  associate(                                                                &
    NU                             =>  plt_site%NU                        , &  
    MaxNodesPerBranch1             =>  pltpar%MaxNodesPerBranch1          , &
    iPlantPhotosynthesisType       =>  plt_photo%iPlantPhotosynthesisType , &    
    NH3Dep2Can_pft                 =>  plt_bgcr%NH3Dep2Can_pft            , &             
    NH3Dep2Can_brch                =>  plt_rbgc%NH3Dep2Can_brch           , &    !
    GrossResp_pft                  =>  plt_bgcr%GrossResp_pft             , &
    CanopyGrosRCO2_pft             =>  plt_bgcr%CanopyGrosRCO2_pft     , &           
    RootGrosRCO2_pft               =>  plt_bgcr%RootGrosRCO2_pft       , &
    GrossCO2Fix_pft                =>  plt_bgcr%GrossCO2Fix_pft        , &
    NumOfBranches_pft              =>  plt_morph%NumOfBranches_pft     , &       
    NodulInfectElms_pft            =>  plt_bgcr%NodulInfectElms_pft    , &    
    LitrfalStrutElms_pft           =>  plt_bgcr%LitrfalStrutElms_pft   , &
    LitrfalStrutElms_pvr           =>  plt_bgcr%LitrfalStrutElms_pvr   , &    
    RootMycoExudElms_pft           =>  plt_rbgc%RootMycoExudElms_pft   , &    
    MY                             =>  plt_morph%MY                    , &        
    MaxSoiL4Root                   =>  plt_morph%MaxSoiL4Root          , &    
    NumRootAxes_pft                =>  plt_morph%NumRootAxes_pft       , &
    MaxNumRootLays                 =>  plt_site%MaxNumRootLays         , &  
    SeasonalNonstElms_pft          =>  plt_biom%SeasonalNonstElms_pft  , &    
    RootElmsBeg_pft                =>  plt_biom%RootElmsBeg_pft        , &
    ShootElmsBeg_pft               =>  plt_biom%ShootElmsBeg_pft       , &
    RootElms_pft                   =>  plt_biom%RootElms_pft           , &    
    ShootElms_brch                 =>  plt_biom%ShootElms_brch         , &
    RootStrutElms_pft              =>  plt_biom%RootStrutElms_pft      , &
    ShootC4NonstC_brch             =>  plt_biom%ShootC4NonstC_brch     , &    
    RCO2A_pvr                      =>  plt_rbgc%RCO2A_pvr              , & 
    ShootElms_pft                  =>  plt_biom%ShootElms_pft            &
  )

  CALL SumPlantBiomStates(I,J,NZ,header)
  
  if(RootElms_pft(ielmc,NZ)>1.e16 .or. ShootElms_pft(ielmc,NZ)>1.e16)stop

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
    DO L=NU,MaxSoiL4Root(NZ)
      RootGrosRCO2_pft(NZ)=RootGrosRCO2_pft(NZ)+RCO2A_pvr(N,L,NZ)
    ENDDO     
  ENDDO  

  GrossResp_pft(NZ)=RootGrosRCO2_pft(NZ)+CanopyGrosRCO2_pft(NZ)

  balc=RootElms_pft(ielmc,NZ)+ShootElms_pft(ielmc,NZ)-RootElmsBeg_pft(ielmc,NZ)-ShootElmsBeg_pft(ielmc,NZ)&
    -GrossCO2Fix_pft(NZ)-GrossResp_pft(NZ)+LitrfalStrutElms_pft(ielmc,NZ)-NodulInfectElms_pft(ielmc,NZ) &
    -RootMycoExudElms_pft(ielmc,NZ)
  if(I>=154)then
  if(NZ==1)THEN
    WRITE(111,*)I+J/24.,trim(header)//' BALC=',BALC,SeasonalNonstElms_pft(ielmc,NZ)&
      ,RootElmsbeg_pft(ielmc,NZ),ShootElmsbeg_pft(ielmc,NZ)
    WRITE(111,*)'cfix=',GrossCO2Fix_pft(NZ),'Rauto=',GrossResp_pft(NZ),'exud=',RootMycoExudElms_pft(ielmc,NZ) 
    WRITE(111,*)'litrf=',LitrfalStrutElms_pft(ielmc,NZ),'infec=',NodulInfectElms_pft(ielmc,NZ)
    write(111,*)RootElms_pft(ielmc,NZ),ShootElms_pft(ielmc,NZ),RootElmsBeg_pft(ielmc,NZ),ShootElmsBeg_pft(ielmc,NZ),'|'
  else
    WRITE(112,*)I+J/24.,trim(header)//' BALC=',BALC,SeasonalNonstElms_pft(ielmc,NZ)&
      ,RootElmsbeg_pft(ielmc,NZ),ShootElmsbeg_pft(ielmc,NZ)
    WRITE(112,*)'cfix=',GrossCO2Fix_pft(NZ),'Rauto=',GrossResp_pft(NZ),'exud=',RootMycoExudElms_pft(ielmc,NZ) 
    WRITE(112,*)'litrf=',LitrfalStrutElms_pft(ielmc,NZ),'infec=',NodulInfectElms_pft(ielmc,NZ)
    write(112,*)RootElms_pft(ielmc,NZ),ShootElms_pft(ielmc,NZ),RootElmsBeg_pft(ielmc,NZ),ShootElmsBeg_pft(ielmc,NZ),'|'
  ENDIF  
  endif
  end associate
  end subroutine SumPlantBiom

!------------------------------------------------------------------------------------------


  subroutine SumPlantBiomStates(I,J,NZ,header)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  character(len=*),intent(in) :: header
  integer :: NB,NE,K
  real(r8) :: root1st,root2nd
  real(r8) :: balc

!     begin_execution
  associate(                                                     &
    NU                             =>  plt_site%NU      , &  
    iPlantNfixType                 =>  plt_morph%iPlantNfixType  , &
    MaxNodesPerBranch1             =>  pltpar%MaxNodesPerBranch1    , &
    CMassHCO3BundleSheath_node     =>  plt_photo%CMassHCO3BundleSheath_node      , &
    CMassCO2BundleSheath_node      =>  plt_photo%CMassCO2BundleSheath_node      , &
    CPOOL3_node                    =>  plt_photo%CPOOL3_node    , &
    CPOOL4_node                    =>  plt_photo%CPOOL4_node    , &
    iPlantPhotosynthesisType       =>  plt_photo%iPlantPhotosynthesisType   , &    
    NodulInfectElms_pft            =>  plt_bgcr%NodulInfectElms_pft, &    
    NumOfBranches_pft              =>  plt_morph%NumOfBranches_pft    , &   
    MaxSoiL4Root                   =>  plt_morph%MaxSoiL4Root     , &    
    NumRootAxes_pft                =>  plt_morph%NumRootAxes_pft   , &
    MaxNumRootLays                 =>  plt_site%MaxNumRootLays      , &  
    SeasonalNonstElms_pft          =>  plt_biom%SeasonalNonstElms_pft     , &    
    RootElmsBeg_pft                =>  plt_biom%RootElmsBeg_pft     , &
    ShootElmsBeg_pft               =>  plt_biom%ShootElmsBeg_pft    , &
    RootMyco1stStrutElms_rpvr      =>  plt_biom%RootMyco1stStrutElms_rpvr  , &
    RootMyco2ndStrutElms_rpvr      =>  plt_biom%RootMyco2ndStrutElms_rpvr  , &    
    RootMycoNonstElms_rpvr         =>  plt_biom%RootMycoNonstElms_rpvr , &
    ShootStrutElms_brch            =>  plt_biom%ShootStrutElms_brch    , &
    RootNodulStrutElms_pvr         =>  plt_biom%RootNodulStrutElms_pvr , &
    EarStrutElms_brch              =>  plt_biom%EarStrutElms_brch    , &
    StalkStrutElms_brch            =>  plt_biom%StalkStrutElms_brch    , &
    PetoleStrutElms_brch           =>  plt_biom%PetoleStrutElms_brch   , &
    HuskStrutElms_brch             =>  plt_biom%HuskStrutElms_brch    , &
    StalkRsrvElms_brch             =>  plt_biom%StalkRsrvElms_brch    , &
    CanopyNonstElms_brch           =>  plt_biom%CanopyNonstElms_brch     , &
    GrainStrutElms_brch            =>  plt_biom%GrainStrutElms_brch  , &    
    RootNodulNonstElms_pvr         =>  plt_biom%RootNodulNonstElms_pvr , &    
    LeafStrutElms_brch             =>  plt_biom%LeafStrutElms_brch , &    
    CanopyNodulNonstElms_brch      =>  plt_biom%CanopyNodulNonstElms_brch, &
    CanopyNodulStrutElms_brch      =>  plt_biom%CanopyNodulStrutElms_brch, &
    StandDeadStrutElms_pft         =>  plt_biom%StandDeadStrutElms_pft   , & 
    StandDeadKCompElms_pft         =>  plt_biom%StandDeadKCompElms_pft   , &       
    RootElms_pft                   =>  plt_biom%RootElms_pft    , &    
    ShootElms_brch                 =>  plt_biom%ShootElms_brch  , &
    RootStrutElms_pft              =>  plt_biom%RootStrutElms_pft, &
    ShootC4NonstC_brch             =>  plt_biom%ShootC4NonstC_brch, &            
    ShootElms_pft                  =>  plt_biom%ShootElms_pft     &
  )
  
  !shoots

  DO NB=1,NumOfBranches_pft(NZ)
    DO NE=1,NumPlantChemElms
      ShootStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ) &
        +PetoleStrutElms_brch(NE,NB,NZ)+StalkStrutElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ) &
        +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ) 
        
      ShootElms_brch(NE,NB,NZ)=ShootStrutElms_brch(NE,NB,NZ)+CanopyNonstElms_brch(NE,NB,NZ)

    ENDDO
  ENDDO

  !add C4 specific reserve carbon
  IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN  
    D3201: DO NB=1,NumOfBranches_pft(NZ)
      ShootC4NonstC_brch(NB,NZ)=0._r8
      D3251: DO K=1,MaxNodesPerBranch1
        ShootC4NonstC_brch(NB,NZ)=ShootC4NonstC_brch(NB,NZ)+CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
          +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ)
      ENDDO D3251
      ShootElms_brch(ielmc,NB,NZ)=ShootElms_brch(ielmc,NB,NZ)+ShootC4NonstC_brch(NB,NZ)
    ENDDO D3201
  ENDIF

  DO NE=1,NumPlantChemElms
    ShootElms_pft(NE,NZ)=sum(ShootElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))+SeasonalNonstElms_pft(NE,NZ)
  ENDDO


  !add nodule, currently, a plant either has canopy or root N-fixing symbiosis, not both
  IF(is_plant_N2fix(iPlantNfixType(NZ)))THEN
    IF(is_canopy_N2fix(iPlantNfixType(NZ)))THEN
      DO NE=1,NumPlantChemElms
        ShootElms_pft(NE,NZ)=ShootElms_pft(NE,NZ)  &
          +sum(CanopyNodulStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)) &
          +sum(CanopyNodulNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
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
    ShootElms_pft(NE,NZ)=ShootElms_pft(NE,NZ)+StandDeadStrutElms_pft(NE,NZ)
  ENDDO

  call SumRootBiome(NZ,RootElms_pft(:,NZ))

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


  subroutine ZeroGrosub()
  
  implicit none
  integer :: NZ,K,L,M,NE

  associate(                                                       &
    NP0                     =>  plt_site%NP0                     , &  
    MY                      =>  plt_morph%MY                     , &     
    NU                      =>  plt_site%NU                      , &      
    MaxSoiL4Root            =>  plt_morph%MaxSoiL4Root           , &        
    MaxNumRootLays          =>  plt_site%MaxNumRootLays          , &
    CO2NetFix_pft           =>  plt_bgcr%CO2NetFix_pft           , &
    LitrfalStrutElms_pft    =>  plt_bgcr%LitrfalStrutElms_pft    , &
    NodulInfectElms_pft     =>  plt_bgcr%NodulInfectElms_pft     , &
    RootMycoExudElms_pft    =>  plt_rbgc%RootMycoExudElms_pft    , &
    NH3Dep2Can_pft          =>  plt_bgcr%NH3Dep2Can_pft          , &             
    GrossResp_pft           =>  plt_bgcr%GrossResp_pft           , &
    GrossCO2Fix_pft         =>  plt_bgcr%GrossCO2Fix_pft         , &    
    CanopyGrosRCO2_pft      =>  plt_bgcr%CanopyGrosRCO2_pft      , &               
    RCO2A_pvr               =>  plt_rbgc%RCO2A_pvr               , &     
    LitrfalStrutElms_pvr    =>  plt_bgcr%LitrfalStrutElms_pvr      &
  )  
  
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
    RCO2A_pvr(1:MY(NZ),NU:MaxSoiL4Root(NZ),NZ)=0._r8
    NH3Dep2Can_pft(NZ)=0._r8
    GrossResp_pft(NZ)=0._r8
    GrossCO2Fix_pft(NZ)=0._r8    
    CO2NetFix_pft(NZ)=0._r8
    CanopyGrosRCO2_pft(NZ)=0._r8
    RootMycoExudElms_pft(1:NumPlantChemElms,NZ)=0._r8
    NodulInfectElms_pft(1:NumPlantChemElms,NZ)=0._r8
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
  integer :: NE
  real(r8) :: massr1st1(NumPlantChemElms)
  real(r8) :: massr2nd1(NumPlantChemElms)
  real(r8) :: massnonst1(NumPlantChemElms)
  real(r8) :: massnodul1(NumPlantChemElms)

  associate(                                                                   &
    NU                            =>  plt_site%NU                            , &  
    MY                            =>  plt_morph%MY                           , &                  
    MaxNumRootLays                =>  plt_site%MaxNumRootLays                , &  
    MaxSoiL4Root                  =>  plt_morph%MaxSoiL4Root                 , &        
    iPlantNfixType                =>  plt_morph%iPlantNfixType               , &    
    NumRootAxes_pft               =>  plt_morph%NumRootAxes_pft              , &  
    RootNodulStrutElms_pvr        =>  plt_biom%RootNodulStrutElms_pvr        , &   
    RootNodulNonstElms_pvr        =>  plt_biom%RootNodulNonstElms_pvr        , &
    RootMyco1stStrutElms_rpvr     =>  plt_biom%RootMyco1stStrutElms_rpvr     , &
    RootMyco2ndStrutElms_rpvr     =>  plt_biom%RootMyco2ndStrutElms_rpvr     , &    
    RootMycoNonstElms_rpvr        =>  plt_biom%RootMycoNonstElms_rpvr          &
  )
  massr1st1=0._r8;massr2nd1=0._r8;massnonst1=0._r8;massnodul1=0._r8
  mass_roots=0._r8
  DO NE=1,NumPlantChemElms
    massr1st1(NE)=sum(RootMyco1stStrutElms_rpvr(NE,1:MY(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    massr2nd1(NE)=sum(RootMyco2ndStrutElms_rpvr(NE,1:MY(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    massnonst1(NE)=sum(RootMycoNonstElms_rpvr(NE,1:MY(NZ),NU:MaxNumRootLays,NZ))
    mass_roots(NE)=massr1st1(NE)+massr2nd1(NE)+massnonst1(NE)
    !add reserve to struct
    if(is_plant_N2fix(iPlantNfixType(NZ)) .and. is_root_N2fix(iPlantNfixType(NZ)))THEN
      massnodul1(NE)=sum(RootNodulStrutElms_pvr(NE,NU:MaxSoiL4Root(NZ),NZ))+&
          sum(RootNodulNonstElms_pvr(NE,NU:MaxSoiL4Root(NZ),NZ))
      mass_roots(NE)=mass_roots(NE)+massnodul1(NE)    
    endif      
  ENDDO
  if(present(massr1st))massr1st=massr1st1
  if(present(massr2nd))massr2nd=massr2nd1
  if(present(massnonst))massnonst=massnonst1
  if(present(massnodul))massnodul=massnodul1
  end associate
  end subroutine SumRootBiome
end module PlantBalMod