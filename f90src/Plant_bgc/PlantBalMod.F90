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
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine SumPlantBiom(I,J,NZ,header)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), intent(in) :: header
  integer :: L,K,N,NE,NB,M
  real(r8) :: root1st,root2nd

!     begin_execution
  associate(                                                         &
    NU                       => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NH3Dep2Can_brch          => plt_rbgc%NH3Dep2Can_brch            ,& !input  :gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
    CanopyGrosRCO2_pft       => plt_bgcr%CanopyGrosRCO2_pft         ,& !input  :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    NumOfBranches_pft        => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    LitrfalStrutElms_pvr     => plt_bgcr%LitrfalStrutElms_pvr       ,& !input  :plant LitrFall element, [g d-2 h-1]
    Myco_pft                 => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    MaxSoiL4Root_pft         => plt_morph%MaxSoiL4Root_pft          ,& !input  :maximum soil layer number for all root axes,[-]
    MaxNumRootLays           => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    RootElms_pft             => plt_biom%RootElms_pft               ,& !input  :plant root element mass, [g d-2]
    RootCO2Autor_pvr         => plt_rbgc%RootCO2Autor_pvr           ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    ShootElms_pft            => plt_biom%ShootElms_pft              ,& !input  :current time whole plant shoot element mass, [g d-2]
    NH3Dep2Can_pft           => plt_bgcr%NH3Dep2Can_pft             ,& !inoput :canopy NH3 flux, [g d-2 h-1]
    RotoAutoCO2_pft          => plt_bgcr%RotoAutoCO2_pft            ,& !inoput :root autotrophic respiraiton, [gC d-2]
    LitrfalStrutElms_pft     => plt_bgcr%LitrfalStrutElms_pft       ,& !inoput :plant element LitrFall, [g d-2 h-1]
    GrossResp_pft            => plt_bgcr%GrossResp_pft               & !output :total plant respiration, [gC d-2 ]
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

  RotoAutoCO2_pft(NZ)=0._r8
  DO N=1,Myco_pft(NZ)
    DO L=NU,MaxSoiL4Root_pft(NZ)
      RotoAutoCO2_pft(NZ)=RotoAutoCO2_pft(NZ)+RootCO2Autor_pvr(N,L,NZ)
    ENDDO     
  ENDDO  

  GrossResp_pft(NZ)=RotoAutoCO2_pft(NZ)+CanopyGrosRCO2_pft(NZ)

  end associate
  end subroutine SumPlantBiom

!----------------------------------------------------------------------------------------------------
  subroutine SumPlantRootGas(I,J,NZ1)

  implicit none
  integer, intent(in) :: I,J
  INTEGER, optional, INTENT(IN) :: NZ1
  integer :: NZ ,L,N,idg,K
  real(r8) :: trcg(idg_beg:idg_NH3)
  associate(                                         &
    trcg_rootml_pvr  => plt_rbgc%trcg_rootml_pvr    ,& !input  :root gas content, [g d-2]
    trcs_rootml_pvr  => plt_rbgc%trcs_rootml_pvr    ,& !input  :root aqueous content, [g d-2]
    Myco_pft         => plt_morph%Myco_pft          ,& !input  :mycorrhizal type (no or yes),[-]
    MaxSoiL4Root_pft => plt_morph%MaxSoiL4Root_pft  ,& !input  :maximum soil layer number for all root axes,[-]
    trcg_root_vr     => plt_rbgc%trcg_root_vr        & !inoput :total root internal gas flux, [g d-2 h-1]
  )
  trcg_root_vr(idg_beg:idg_NH3,:)         = 0._r8

  trcg(:)=0._r8
  IF(present(NZ1))THEN
    NZ=NZ1
    DO L=1,MaxSoiL4Root_pft(NZ)
      DO N=1,Myco_pft(NZ)  
        DO idg=idg_beg,idg_NH3
          trcg_root_vr(idg,L)=trcg_root_vr(idg,L)+trcs_rootml_pvr(idg,N,L,NZ)+trcg_rootml_pvr(idg,N,L,NZ)        
        ENDDO        
      ENDDO  
      DO idg=idg_beg,idg_NH3
        trcg(idg)=trcg(idg)+trcg_root_vr(idg,L)
      ENDDO
    ENDDO
  ELSE
    DO NZ=1,plt_site%NP  

      DO L=1,MaxSoiL4Root_pft(NZ)
        DO N=1,Myco_pft(NZ)  
          DO idg=idg_beg,idg_NH3
            trcg_root_vr(idg,L)=trcg_root_vr(idg,L)+trcs_rootml_pvr(idg,N,L,NZ)+trcg_rootml_pvr(idg,N,L,NZ)        
          ENDDO        
        ENDDO  
        DO idg=idg_beg,idg_NH3
          trcg(idg)=trcg(idg)+trcg_root_vr(idg,L)
        ENDDO
      ENDDO
    ENDDO
  ENDIF  

  end associate
  end subroutine SumPlantRootGas

!----------------------------------------------------------------------------------------------------
  subroutine SumPlantBiomStates(I,J,NZ,header)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  character(len=*),intent(in) :: header
  integer :: NB,NE,K,jk
  real(r8) :: root1st,root2nd
  real(r8) :: balc

!     begin_execution
  associate(                                                          &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    RootElmsBeg_pft           => plt_biom%RootElmsBeg_pft            ,& !input  :plant root element at previous time step, [g d-2]    
    ShootElmsBeg_pft          => plt_biom%ShootElmsBeg_pft           ,& !input  :previous whole plant shoot element mass, [g d-2]    
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch  ,& !input  :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch  ,& !input  :branch nodule structural element, [g d-2]
    StandDeadKCompElms_pft    => plt_biom%StandDeadKCompElms_pft     ,& !input  :standing dead element fraction, [g d-2]
    ShootStrutElms_brch       => plt_biom%ShootStrutElms_brch        ,& !input  :branch shoot structural element mass, [g d-2]
    RootElms_pft              => plt_biom%RootElms_pft               ,& !input  :plant root element mass, [g d-2]
    ShootElms_pft             => plt_biom%ShootElms_pft              ,& !inoput :current time whole plant shoot element mass, [g d-2]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft      & !output :standing dead element, [g d-2]
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

!----------------------------------------------------------------------------------------------------
  subroutine SumPlantBranchBiome(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  integer :: NE,K
  associate(                                                             &
    MaxNodesPerBranch1         => pltpar%MaxNodesPerBranch1             ,& !input  :maximum number of canopy nodes, 25
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft           ,& !input  :number of branches,[-]
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch           ,& !input  :branch leaf structural element mass, [g d-2]
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch         ,& !input  :branch sheath structural element, [g d-2]
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch          ,& !input  :branch stalk structural element mass, [g d-2]
    LeafPetolBiomassC_brch     => plt_biom%LeafPetolBiomassC_brch       ,& !input  :plant branch leaf + sheath C, [g d-2]
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch           ,& !input  :branch reserve element mass, [g d-2]
    HuskStrutElms_brch         => plt_biom%HuskStrutElms_brch           ,& !input  :branch husk structural element mass, [g d-2]
    EarStrutElms_brch          => plt_biom%EarStrutElms_brch            ,& !input  :branch ear structural chemical element mass, [g d-2]
    GrainStrutElms_brch        => plt_biom%GrainStrutElms_brch          ,& !input  :branch grain structural element mass, [g d-2]
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch         ,& !input  :branch nonstructural element, [g d-2]
    CPOOL3_node                => plt_photo%CPOOL3_node                 ,& !input  :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                => plt_photo%CPOOL4_node                 ,& !input  :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    iPlantPhotosynthesisType   => plt_photo%iPlantPhotosynthesisType    ,& !input  :plant photosynthetic type (C3 or C4),[-]
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node  ,& !input  :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node   ,& !input  :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    C4PhotoShootNonstC_brch    => plt_biom%C4PhotoShootNonstC_brch      ,& !inoput :branch shoot nonstrucal elelment, [g d-2]
    ShootStrutElms_brch        => plt_biom%ShootStrutElms_brch          ,& !inoput :branch shoot structural element mass, [g d-2]
    CanopyLeafShethC_pft       => plt_biom%CanopyLeafShethC_pft          & !output :canopy leaf + sheath C, [g d-2]
  )
  DO NE=1,NumPlantChemElms
    ShootStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ) &
      +PetoleStrutElms_brch(NE,NB,NZ)+StalkStrutElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ) &
      +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ) &
      +CanopyNonstElms_brch(NE,NB,NZ)
    
  ENDDO
!  CanopySapwoodC_pft(NZ)     = sum(SapwoodBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafShethC_pft(NZ) = sum(LeafPetolBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))

  !add C4 specific reserve carbon
  C4PhotoShootNonstC_brch(NB,NZ)=0._r8
  IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN      
    D3251: DO K=1,MaxNodesPerBranch1
      C4PhotoShootNonstC_brch(NB,NZ)=C4PhotoShootNonstC_brch(NB,NZ)+CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D3251
    ShootStrutElms_brch(ielmc,NB,NZ)=ShootStrutElms_brch(ielmc,NB,NZ)+C4PhotoShootNonstC_brch(NB,NZ)
  ENDIF

  end associate
  end subroutine SumPlantBranchBiome

!----------------------------------------------------------------------------------------------------
  subroutine ZeroGrosub()
  
  implicit none
  integer :: NZ,K,L,M,NE

  associate(                                                  &
    NP0                   => plt_site%NP0                    ,& !input  :intitial number of plant species,[-]
    Myco_pft              => plt_morph%Myco_pft              ,& !input  :mycorrhizal type (no or yes),[-]
    NU                    => plt_site%NU                     ,& !input  :current soil surface layer number, [-]
    MaxSoiL4Root_pft      => plt_morph%MaxSoiL4Root_pft      ,& !input  :maximum soil layer number for all root axes,[-]
    MaxNumRootLays        => plt_site%MaxNumRootLays         ,& !input  :maximum root layer number,[-]
    RootCO2Autor_pvr      => plt_rbgc%RootCO2Autor_pvr       ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    CO2NetFix_pft         => plt_bgcr%CO2NetFix_pft          ,& !output :canopy net CO2 exchange, [gC d-2 h-1]
    RCanMaintDef_CO2_pft  => plt_bgcr%RCanMaintDef_CO2_pft   ,& !output :canopy maintenance respiraiton deficit as CO2, [gC d-2]    
    RootMaintDef_CO2_pvr  => plt_bgcr%RootMaintDef_CO2_pvr   ,& !output :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]        
    NodulInfectElms_pft   => plt_bgcr%NodulInfectElms_pft    ,& !output :nodule infection chemical element mass, [g d-2]
    RootMycoExudElms_pft  => plt_rbgc%RootMycoExudElms_pft   ,& !output :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    NH3Dep2Can_pft        => plt_bgcr%NH3Dep2Can_pft         ,& !output :canopy NH3 flux, [g d-2 h-1]
    GrossResp_pft         => plt_bgcr%GrossResp_pft          ,& !output :total plant respiration, [gC d-2 ]
    GrossCO2Fix_pft       => plt_bgcr%GrossCO2Fix_pft        ,& !output :total gross CO2 fixation, [gC d-2 ]
    CO2FixCL_pft          => plt_rbgc%CO2FixCL_pft           ,& !output :Rubisco-limited CO2 fixation, [gC d-2 h-1]
    CO2FixLL_pft          => plt_rbgc%CO2FixLL_pft           ,& !output :Light-limited CO2 fixation, [gC d-h2 h-1]
    CanopyGrosRCO2_pft    => plt_bgcr%CanopyGrosRCO2_pft     ,& !output :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    LitrfalStrutElms_pvr  => plt_bgcr%LitrfalStrutElms_pvr   ,& !output :plant LitrFall element, [g d-2 h-1]
    PARSunlit_pft         => plt_photo%PARSunlit_pft         ,& !output :PAR absorbed by sunlit leaf, [umol m-2 s-1]
    PARSunsha_pft         => plt_photo%PARSunsha_pft         ,& !output :PAR absorbed by sun-shaded leaf, [umol m-2 s-1]
    CH2OSunlit_pft        => plt_photo%CH2OSunlit_pft        ,& !output :carbon fixation by sun-lit leaf, [gC d-2 h-1]
    CH2OSunsha_pft        => plt_photo%CH2OSunsha_pft         & !output :carbon fixation by sun-shaded leaf, [gC d-2 h-1]    
  )
  
  plt_rbgc%trcs_Soil2plant_uptake_vr=0._r8

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
    PARSunsha_pft(NZ)                                           = 0._r8
    PARSunlit_pft(NZ)                                           = 0._r8
    CH2OSunlit_pft(NZ)                                          = 0._r8
    CH2OSunsha_pft(NZ)                                          = 0._r8
    RootMaintDef_CO2_pvr(:,:,NZ)                                = 0._r8
    RootCO2Autor_pvr(1:Myco_pft(NZ),NU:MaxSoiL4Root_pft(NZ),NZ) = 0._r8
    NH3Dep2Can_pft(NZ)                                          = 0._r8
    GrossResp_pft(NZ)                                           = 0._r8
    GrossCO2Fix_pft(NZ)                                         = 0._r8
    CO2NetFix_pft(NZ)                                           = 0._r8
    RCanMaintDef_CO2_pft(NZ)                                    = 0._r8
    CanopyGrosRCO2_pft(NZ)                                      = 0._r8
    RootMycoExudElms_pft(1:NumPlantChemElms,NZ)                 = 0._r8
    NodulInfectElms_pft(1:NumPlantChemElms,NZ)                  = 0._r8
    CO2FixCL_pft(NZ)                                            = 0._r8
    CO2FixLL_pft(NZ)                                            = 0._r8

  ENDDO D9980
  end associate
  end subroutine ZeroGrosub

!----------------------------------------------------------------------------------------------------
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
  associate(                                                          &
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft          ,& !input  :maximum soil layer number for all root axes,[-]
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft           ,& !input  :root primary axis number,[-]
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr    ,& !input  :root layer nodule element, [g d-2]
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr    ,& !input  :root layer nonstructural element, [g d-2]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !input  :root layer element primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !input  :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !input  :root layer nonstructural element, [g d-2]
    RootMycoNonstElms_pft     => plt_biom%RootMycoNonstElms_pft      ,& !output :nonstructural root-myco chemical element, [g d-2]
    RootStrutElms_pft         => plt_biom%RootStrutElms_pft          ,& !output :plant root structural element mass, [g d-2]
    RootMassElm_pvr           => plt_biom%RootMassElm_pvr             & !output :root biomass in chemical elements, [g d-2]
  )
  massr1st1=0._r8;massr2nd1=0._r8;massnodul1=0._r8
  mass_roots=0._r8

  DO NE=1,NumPlantChemElms
    DO L=NU,MaxNumRootLays
      RootMassElm_pvr(NE,L,NZ)= sum(RootMyco1stStrutElms_rpvr(NE,1:Myco_pft(NZ),L,1:NumRootAxes_pft(NZ),NZ)) + &
        sum(RootMyco2ndStrutElms_rpvr(NE,1:Myco_pft(NZ),L,1:NumRootAxes_pft(NZ),NZ)) + &
        sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),L,NZ))
    ENDDO
    massr1st1(NE)=sum(RootMyco1stStrutElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    massr2nd1(NE)=sum(RootMyco2ndStrutElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    RootStrutElms_pft(NE,NZ)=massr1st1(NE)+massr2nd1(NE)
    DO N=1,Myco_pft(NZ)
      RootMycoNonstElms_pft(NE,N,NZ)=sum(RootMycoNonstElms_rpvr(NE,N,NU:MaxNumRootLays,NZ))
    enddo
    massnonst1(NE)=sum(RootMycoNonstElms_pft(NE,1:Myco_pft(NZ),NZ))
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

!----------------------------------------------------------------------------------------------------
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


!----------------------------------------------------------------------------------------------------
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

  associate(                                                          &
    GrossCO2Fix_pft           => plt_bgcr%GrossCO2Fix_pft            ,& !input  :total gross CO2 fixation, [gC d-2 ]
    GrossResp_pft             => plt_bgcr%GrossResp_pft              ,& !input  :total plant respiration, [gC d-2 ]
    NodulInfectElms_pft       => plt_bgcr%NodulInfectElms_pft        ,& !input  :nodule infection chemical element mass, [g d-2]
    LitrfalStrutElms_pft      => plt_bgcr%LitrfalStrutElms_pft       ,& !input  :plant element LitrFall, [g d-2 h-1]
    RootMycoExudElms_pft      => plt_rbgc%RootMycoExudElms_pft       ,& !input  :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    RootElmsBeg_pft           => plt_biom%RootElmsBeg_pft            ,& !input  :plant root element at previous time step, [g d-2]
    RootElms_pft              => plt_biom%RootElms_pft               ,& !input  :plant root element mass, [g d-2]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft      ,& !input  :plant stored nonstructural element at current step, [g d-2]
    SeasonalNonstElmsbeg_pft  => plt_biom%SeasonalNonstElmsbeg_pft   ,& !input  :plant stored nonstructural element at prev step, [g d-2]
    StandDeadStrutElmsBeg_pft => plt_biom%StandDeadStrutElmsBeg_pft  ,& !input  :standing dead element at previous time step, [g d-2]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft     ,& !input  :standing dead element, [g d-2]
    ShootElmsBeg_pft          => plt_biom%ShootElmsBeg_pft           ,& !input  :previous whole plant shoot element mass, [g d-2]
    ShootElms_pft             => plt_biom%ShootElms_pft               & !input  :current time whole plant shoot element mass, [g d-2]
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
  ![tail]
end module PlantBalMod