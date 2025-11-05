module PlantPhenolMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : AZMAX1
  use InitPlantMod, only : StartPlants
  use EcoSIMCtrlMod, only : etimer
  use DebugToolMod
  use EcosimConst
  use PlantAPIData
  use PlantMathFuncMod

  implicit none

  private


  character(len=*), parameter :: mod_filename = &
  __FILE__

!
! PSIMin4LeafExpansion=
! PSIMin4LeafOut=
! PSIMin4LeafOff=
! GrowStageNorm4VegetaPheno,GrowStageNorm4ReprodPheno=,reproductive phenology
! BranchNumMax=
! MaxHour4LeafOutOff=maximum hours for leafout,leafoff
!
  real(r8), PARAMETER :: PSIMin4LeafExpansion=0.1_r8                             !minimum canopy turgor potential for leaf expansion, [MPa]
  real(r8), parameter :: PSIMin4LeafOut=-0.2_r8                                  !minimum canopy water potential for leafout of drought-deciduous PFT, [MPa]
  real(r8) ,PARAMETER :: GrowStageNorm4VegetaPheno=2.00_r8                       !normalized growth stage durations for vegetative phenology,[-]
  real(r8), PARAMETER :: GrowStageNorm4ReprodPheno=0.667_r8                      !normalized growth stage durations for reproductive phenology,[-]
  real(r8), PARAMETER :: MaxHour4LeafOutOff=3600.0_r8                            !maximum branch number for PFT defined by iPlantTurnoverPattern_pftin PFT file [h]
  real(r8), parameter :: PSIMin4LeafOff(0:3)=real((/-200.0,-2.0,-2.0,-2.0/),r8)  !minimum leaf water potential of leave off, [h]
  integer , parameter :: BranchNumMax(0:3)=(/5,1,1,1/)                           !maximum branch number

  public :: PhenologyUpdate
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine PhenologyUpdate(I,J)
!
!     THIS subroutine CALCULATES PLANT PHENOLOGY
!
  use PlantBalMod, only : SumPlantBiomStates
  implicit none

  integer, intent(in) :: I,J
  INTEGER :: NB, NZ
  integer :: NE
  character(len=*), parameter :: subname='PhenologyUpdate'
! begin_execution
  associate(                                                                   &
    doInitPlant_pft               => plt_pheno%doInitPlant_pft                ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    IsPlantActive_pft             => plt_pheno%IsPlantActive_pft              ,& !input  :flag for living pft, [-]
    iPlantCalendar_brch           => plt_pheno%iPlantCalendar_brch            ,& !input  :plant growth stage, [-]
    DATAP                         => plt_site%DATAP                           ,& !input  :parameter file name,[-]
    NP                            => plt_site%NP                              ,& !input  :current number of plant species,[-]
    PlantPopulation_pft           => plt_site%PlantPopulation_pft             ,& !input  :plant population, [d-2]
    MainBranchNum_pft             => plt_morph%MainBranchNum_pft              ,& !input  :id number of main branch,[-]
    PlantPopu_col                 => plt_site%PlantPopu_col                    & !inoput :total plant population, [plants d-2]
  )
  call PrintInfo('beg '//subname)
  D9985: DO NZ=1,NP

    IF(DATAP(NZ).NE.'NO')THEN
!
!     PlantPopu_col=total biome population
!
      PlantPopu_col=PlantPopu_col+PlantPopulation_pft(NZ)
!
!         SET CROP FLAG ACCORDINGTopRootLayer_pftTO PLANTING, HARVEST DATES, DEATH,
!         1 = ALIVE, 0 = NOT ALIVE
!         DATAP=PFT file name
!
      call set_plant_flags(I,J,NZ)

!
!         INITIALIZE VARIABLES IN ACTIVE PFT
!
      IF(IsPlantActive_pft(NZ).EQ.iActive)THEN

        call FindMainBranchNumber(NZ)

        call StagePlantPhenology(I,J,NZ)

        call TestPlantEmergence(I,J,NZ)

        call root_shoot_branching(I,J,NZ)

!
!           THE REST OF THE subroutine MODELS THE PHENOLOGY OF EACH BRANCH
!
!           doInitLeafOut_brch,doPlantLeafOut_brch=flags for initializing leafout,leafoff
!           Hours4Leafout_brch=leafout hours
!
        IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 .OR. doInitPlant_pft(NZ).EQ.itrue)THEN          
          call Emerged_plant_Phenology(I,J,NZ)
        ENDIF        
      ENDIF
      
    ENDIF
  ENDDO D9985

  call PrintInfo('end '//subname)
  RETURN
  end associate
  END subroutine PhenologyUpdate

!----------------------------------------------------------------------------------------------------
  subroutine Emerged_plant_Phenology(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='Emerged_plant_Phenology'
  integer :: NB
  integer :: LeafNumberGrowing
  associate(                                                               &
    iPlantBranchState_brch      => plt_pheno%iPlantBranchState_brch       ,& !input  :flag to detect branch death, [-]
    LeafNumberAtFloralInit_brch => plt_pheno%LeafNumberAtFloralInit_brch  ,& !input  :leaf number at floral initiation, [-]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                   ,& !input  :canopy total water potential, [Mpa]
    NumOfLeaves_brch            => plt_morph%NumOfLeaves_brch             ,& !input  :leaf number, [-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft        ,& !input  :plant growth type (vascular, non-vascular),[-]
    HoursTooLowPsiCan_pft       => plt_pheno%HoursTooLowPsiCan_pft        ,& !inoput :canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
    KHiestGroLeafNode_brch      => plt_pheno%KHiestGroLeafNode_brch       ,& !inoput :leaf growth stage counter, [-]
    doRemobilization_brch       => plt_pheno%doRemobilization_brch        ,& !output :branch phenology flag, [-]
    KLeafNumber_brch            => plt_morph%KLeafNumber_brch              & !output :leaf number, [-]
  )
  call PrintInfo('beg '//subname)
  D2010: DO NB=1,NumOfBranches_pft(NZ)

    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      call live_branch_phenology(I,J,NB,nz)
    ENDIF
!
!           KHiestGroLeafNode_brch=integer of most recent leaf number currently growing
!
    LeafNumberGrowing=KHiestGroLeafNode_brch(NB,NZ)
    IF(LeafNumberAtFloralInit_brch(NB,NZ).LE.ppmc)THEN
      KHiestGroLeafNode_brch(NB,NZ)=INT(NumOfLeaves_brch(NB,NZ))+1
    ELSE
      KHiestGroLeafNode_brch(NB,NZ)=INT(AMIN1(NumOfLeaves_brch(NB,NZ),LeafNumberAtFloralInit_brch(NB,NZ)))+1
    ENDIF
    KLeafNumber_brch(NB,NZ)=MIN(MaxNodesPerBranch1-1,KHiestGroLeafNode_brch(NB,NZ))
    IF(KHiestGroLeafNode_brch(NB,NZ).GT.LeafNumberGrowing)THEN
      doRemobilization_brch(NB,NZ)=itrue
    ELSE
      doRemobilization_brch(NB,NZ)=ifalse
    ENDIF
!
    call branch_specific_phenology(I,J,NB,NZ)

  ENDDO D2010
!
!             WATER STRESS INDICATOR
!
!             PSICanopy_pft=canopy total water potential
!             PSIMin4LeafOff=minimum canopy water potential for leafoff
!             HoursTooLowPsiCan_pft=number of hours PSICanopy_pft(< PSIMin4LeafOff (for output only)
!
  IF(PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iPlantRootProfile_pft(NZ)))THEN
    HoursTooLowPsiCan_pft(NZ)=HoursTooLowPsiCan_pft(NZ)+1.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine Emerged_plant_Phenology          

!----------------------------------------------------------------------------------------------------
  subroutine set_plant_flags(I,J,NZ)
!
! Description
!
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='set_plant_flags'
  INTEGER :: L
  logical :: HarvestChk,PlantingChk,postharvst_check

! begin_execution
  associate(                                                   &
    iYearPlanting_pft     => plt_distb%iYearPlanting_pft      ,& !input  :year of planting,[-]
    iYearPlantHarvest_pft => plt_distb%iYearPlantHarvest_pft  ,& !input  :year of harvest,[-]
    iDayPlanting_pft      => plt_distb%iDayPlanting_pft       ,& !input  :day of planting,[-]
    iDayPlantHarvest_pft  => plt_distb%iDayPlantHarvest_pft   ,& !input  :day of harvest,[-]
    DATAP                 => plt_site%DATAP                   ,& !input  :parameter file name,[-]
    iYearCurrent          => plt_site%iYearCurrent            ,& !input  :current year,[-]
    SeedPlantedElm_pft      => plt_biom%SeedPlantedElm_pft        ,& !input  :plant stored nonstructural C at planting, [gC d-2]
    NumActivePlants       => plt_site%NumActivePlants         ,& !inoput :number of active PFT in the grid, [-]
    Eco_NBP_CumYr_col     => plt_bgcr%Eco_NBP_CumYr_col       ,& !inoput :total NBP, [g d-2]
    iPlantState_pft       => plt_pheno%iPlantState_pft        ,& !inoput :flag for species death, [-]
    IsPlantActive_pft     => plt_pheno%IsPlantActive_pft       & !output :flag for living pft, [-]
  )
  call PrintInfo('beg '//subname)
  !first hour of day
  IF(J.EQ.1)THEN
    HarvestChk=iDayPlanting_pft(NZ).LE.iDayPlantHarvest_pft(NZ) .OR. iYearPlanting_pft(NZ).LT.iYearPlantHarvest_pft(NZ)

    !Before harvest  
    IF(HarvestChk)THEN
      !planting is feasible
      IF(I.GE.iDayPlanting_pft(NZ) .OR. iYearCurrent.GT.iYearPlanting_pft(NZ))THEN
        !planted 
        postharvst_check=I.GT.iDayPlantHarvest_pft(NZ) .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ) .AND. iPlantState_pft(NZ).EQ.iDead

        IF(postharvst_check)THEN
          !post harvest
          IsPlantActive_pft(NZ)=iDormant
        ELSE
          IF(I.EQ.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ))THEN
            !planting day of year
            IsPlantActive_pft(NZ) = iDormant
            iPlantState_pft(NZ)   = iLive
            CALL StartPlants(NZ,NZ)
            Eco_NBP_CumYr_col=Eco_NBP_CumYr_col+SeedPlantedElm_pft(ielmc,NZ)
          ENDIF

          !the living plant has actual properties set
          IF(DATAP(NZ).NE.'NO' .AND. iPlantState_pft(NZ).EQ.iLive)then
            IsPlantActive_pft(NZ)=iActive
          endif
        ENDIF
      ELSE
        IsPlantActive_pft(NZ)=iDormant
      ENDIF
      !After harvest
    ELSE
      HarvestChk  = I.LT.iDayPlanting_pft(NZ) .AND. I.GT.iDayPlantHarvest_pft(NZ) .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ)
      PlantingChk = I.LT.iDayPlanting_pft(NZ) .AND. iYearPlanting_pft(NZ).GT.iYearPlantHarvest_pft(NZ)

      !not planted

      IF((HarvestChk .AND. iPlantState_pft(NZ).EQ.iDead) .OR. PlantingChk)THEN
        IsPlantActive_pft(NZ)=iDormant
      ELSE
        !planting
        IF(I.EQ.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ))THEN
          IsPlantActive_pft(NZ) = iDormant
          iPlantState_pft(NZ)   = iLive
          CALL StartPlants(NZ,NZ)
          Eco_NBP_CumYr_col = Eco_NBP_CumYr_col+SeedPlantedElm_pft(ielmc,NZ)
        ENDIF
        IF(DATAP(NZ).NE.'NO' .AND. iPlantState_pft(NZ).EQ.iLive)then
          IsPlantActive_pft(NZ)=iActive
        endif
      ENDIF
    ENDIF
    NumActivePlants=NumActivePlants+IsPlantActive_pft(NZ)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine set_plant_flags

!----------------------------------------------------------------------------------------------------
  subroutine root_shoot_branching(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: NB
  character(len=*), parameter :: subname='root_shoot_branching'
  integer :: BranchNumber_new
! begin_execution
  associate(                                                                 &
    CanopyNonstElmConc_pft       => plt_biom%CanopyNonstElmConc_pft         ,& !input  :canopy nonstructural element concentration, [g d-2]
    SeasonalNonstElms_pft        => plt_biom%SeasonalNonstElms_pft          ,& !input  :plant stored nonstructural element at current step, [g d-2]
    iPlantCalendar_brch          => plt_pheno%iPlantCalendar_brch           ,& !input  :plant growth stage, [-]
    doInitPlant_pft              => plt_pheno%doInitPlant_pft               ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    iPlantPhenolPattern_pft      => plt_pheno%iPlantPhenolPattern_pft       ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantTurnoverPattern_pft    => plt_pheno%iPlantTurnoverPattern_pft     ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    MinNonstC2InitRoot_pft       => plt_pheno%MinNonstC2InitRoot_pft        ,& !input  :threshold root nonstructural C content for initiating new root axis, [gC gC-1]
    MatureGroup_pft              => plt_pheno%MatureGroup_pft               ,& !input  :acclimated plant maturity group, [-]
    NonstCMinConc2InitBranch_pft => plt_pheno%NonstCMinConc2InitBranch_pft  ,& !input  :branch nonstructural C content required for new branch, [gC gC-1]
    PlantPopulation_pft          => plt_site%PlantPopulation_pft            ,& !input  :plant population, [d-2]
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                   ,& !input  :root turgor water potential, [Mpa]
    FracGroth2Node_pft           => plt_allom%FracGroth2Node_pft            ,& !input  :parameter for allocation of growth to nodes, [-]
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft             ,& !input  :id number of main branch,[-]
    NumCogrowthNode_pft          => plt_morph%NumCogrowthNode_pft           ,& !input  :number of concurrently growing nodes,[-]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    ShootNodeNumAtPlanting_pft   => plt_morph%ShootNodeNumAtPlanting_pft    ,& !input  :number of nodes in seed, [-]
    ShootNodeNum_brch            => plt_morph%ShootNodeNum_brch             ,& !input  :shoot node number, [-]
    iPlantBranchState_brch       => plt_pheno%iPlantBranchState_brch        ,& !inoput :flag to detect branch death, [-]
    NumPrimeRootAxes_pft         => plt_morph%NumPrimeRootAxes_pft          ,& !inoput :root primary axis number,[-]
    NumOfBranches_pft            => plt_morph%NumOfBranches_pft             ,& !inoput :number of branches,[-]
    BranchNumber_pft             => plt_morph%BranchNumber_pft              ,& !inoput :main branch numeric id,[-]
    MatureGroup_brch             => plt_pheno%MatureGroup_brch              ,& !output :plant maturity group, [-]
    iPlantRootState_pft          => plt_pheno%iPlantRootState_pft           ,& !output :flag to detect root system death,[-]
    iPlantShootState_pft         => plt_pheno%iPlantShootState_pft          ,& !output :flag to detect canopy death,[-]
    Hours4Leafout_brch           => plt_pheno%Hours4Leafout_brch            ,& !output :heat requirement for spring leafout/dehardening, [h]
    BranchNumerID_brch            => plt_morph%BranchNumerID_brch              & !output :branch meric id, [-]
  )
  call PrintInfo('beg '//subname)
!
! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
! CONCENTRATION PERMIT
!
! doInitPlant_pft=PFT initialization flag:0=no,1=yes
! PSIRootTurg_vr=root turgor potential
! iPlantPhenolPattern_pft=growth habit from PFT file
! iPlantCalendar_brch(ipltcal_InitFloral,=floral initiation date
! NumOfBranches_pft=primary root axis number
! WTRVC=nonstructural C storage
! PB=nonstructural C concentration needed for branching
! iPlantBranchState_brch=branch life flag:0=living,1=dead
! PSTG=node number
! FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
! NumCogrowthNode_pft=number of concurrently growing nodes
! ShootNodeNumAtPlanting_pft,GROUP=node number at planting,floral initiation
! IBTYP: setup for phenologically-driven above-ground turnover

  IF(doInitPlant_pft(NZ).EQ.ifalse)THEN
    !plant initialized
    IF(J.EQ.1 .AND. PlantPopulation_pft(NZ).GT.0.0_r8)THEN
      !first hour of the day, population > 0
      IF(PSIRootTurg_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ).GT.PSIMin4LeafExpansion)THEN
        IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .OR. iPlantCalendar_brch(ipltcal_InitFloral,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
          !perennial plant or flower not initiated for annual plant 

          IF((NumOfBranches_pft(NZ).EQ.0 .AND. SeasonalNonstElms_pft(ielmc,NZ).GT.0.0_r8) &
            .OR. (CanopyNonstElmConc_pft(ielmc,NZ).GT.NonstCMinConc2InitBranch_pft(NZ) &
            .AND. NonstCMinConc2InitBranch_pft(NZ).GT.0.0_r8))THEN

            D120: DO NB=1,MaxNumBranches
              IF(iPlantBranchState_brch(NB,NZ).EQ.iDead)THEN
                BranchNumber_new=BranchNumber_pft(NZ)+NumCogrowthNode_pft(NZ)/FracGroth2Node_pft(NZ)+ShootNodeNumAtPlanting_pft(NZ)
                IF(NB.EQ.MainBranchNum_pft(NZ) .OR. ShootNodeNum_brch(MainBranchNum_pft(NZ),NZ).GT.BranchNumber_new)THEN
                  !initiate a new branch
                  BranchNumber_pft(NZ)          = BranchNumber_pft(NZ)+1
                  NumOfBranches_pft(NZ)         = MIN(BranchNumMax(iPlantTurnoverPattern_pft(NZ)),MAX(NB,NumOfBranches_pft(NZ)))
                  BranchNumerID_brch(NB,NZ)     = BranchNumber_pft(NZ)-1
                  iPlantShootState_pft(NZ)      = iLive
                  iPlantBranchState_brch(NB,NZ) = iLive
                  Hours4Leafout_brch(NB,NZ)     = 0.0_r8
                  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
                    MatureGroup_brch(NB,NZ)=AZMAX1(MatureGroup_pft(NZ)-BranchNumerID_brch(NB,NZ))
                  ELSE
                    MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
                  ENDIF
                  exit
                ENDIF
              ENDIF
            ENDDO D120            
          ENDIF
        ENDIF
      ENDIF
!
!     ADD AXIS TO ROOT IF PLANT GROWTH STAGE, ROOT NON-STRUCTURAL C
!     CONCENTRATION PERMIT
!
!     PR=nonstructural C concentration needed for root branching
!     FracGroth2Node_pft: parameter for allocation of growth to nodes
!     XLI: number of nodes in seed
!     PSTG: node number
!     MainBranchNum_pft: number of main branch
!     CanopyNonstElmConc_pft: canopy nonstructural element concentration
!     PSIRootTurg_vr: root turgor pressure
!     SeasonalNonstElms_pft: non-structural carbon
!     root axis initialization
      !   if(I>176)print*,'rootaxis'
      IF(PSIRootTurg_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ).GT.PSIMin4LeafExpansion)THEN
        IF(NumPrimeRootAxes_pft(NZ).EQ.0 .OR. ShootNodeNum_brch(MainBranchNum_pft(NZ),NZ) &
          .GT.NumPrimeRootAxes_pft(NZ)/FracGroth2Node_pft(NZ)+ShootNodeNumAtPlanting_pft(NZ))THEN
          IF((NumPrimeRootAxes_pft(NZ).EQ.0 .AND. SeasonalNonstElms_pft(ielmc,NZ).GT.0.0_r8) &
            .OR. (CanopyNonstElmConc_pft(ielmc,NZ).GT.MinNonstC2InitRoot_pft(NZ) & 
            .AND. MinNonstC2InitRoot_pft(NZ).GT.0.0_r8))THEN
            NumPrimeRootAxes_pft(NZ)=MIN(NumCanopyLayers1,NumPrimeRootAxes_pft(NZ)+1)
            iPlantRootState_pft(NZ)=iLive
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine root_shoot_branching

!----------------------------------------------------------------------------------------------------
  subroutine FindMainBranchNumber(NZ)
  implicit none
  integer, intent(in) :: NZ
  integer :: NB,BranchNumberX_pft
  character(len=*), parameter :: subname='FindMainBranchNumber'
  associate(                                                     &
    iPlantBranchState_brch => plt_pheno%iPlantBranchState_brch  ,& !input  :flag to detect branch death, [-]
    NumOfBranches_pft      => plt_morph%NumOfBranches_pft       ,& !input  :number of branches,[-]
    BranchNumerID_brch      => plt_morph%BranchNumerID_brch       ,& !input  :branch meric id, [-]
    MainBranchNum_pft      => plt_morph%MainBranchNum_pft        & !output :id number of main branch,[-]
  )
  call PrintInfo('beg '//subname)
  MainBranchNum_pft(NZ) = 1
  BranchNumberX_pft     = 1.0E+06_r8

  DD140: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      !find main branch number, which is the most recent live branch      
      IF(BranchNumerID_brch(NB,NZ).LT.BranchNumberX_pft)THEN
        MainBranchNum_pft(NZ) = NB
        BranchNumberX_pft     = BranchNumerID_brch(NB,NZ)
      ENDIF      
    ENDIF
  ENDDO DD140
  call PrintInfo('end '//subname)
  end associate
  end subroutine FindMainBranchNumber

!----------------------------------------------------------------------------------------------------
  subroutine StagePlantPhenology(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='StagePlantPhenology'
  integer :: NB,N,L,NE

  associate(                                                            &
    CanopyLeafSheathC_pft       => plt_biom%CanopyLeafSheathC_pft      ,& !input  :canopy leaf + sheath C, [g d-2]
    CanopyLeafSheathC_brch      => plt_biom%CanopyLeafSheathC_brch     ,& !input  :plant branch leaf + sheath C, [g d-2]
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch  ,& !input  :branch nodule nonstructural element, [g d-2]
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr     ,& !input  :root layer nonstructural element, [g d-2]
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch       ,& !input  :branch nonstructural element, [g d-2]
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft           ,& !input  :threshold zero for leaf calculation, [-]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    RootMycoActiveBiomC_pvr     => plt_biom%RootMycoActiveBiomC_pvr    ,& !input  :root layer structural C, [gC d-2]
    NU                          => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    iPlantBranchState_brch      => plt_pheno%iPlantBranchState_brch    ,& !input  :flag to detect branch death, [-]
    Myco_pft                    => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NMaxRootBotLayer_pft        => plt_morph%NMaxRootBotLayer_pft      ,& !input  :maximum soil layer number for all root axes, [-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    CanopyNonstElms_pft         => plt_biom%CanopyNonstElms_pft        ,& !inoput :canopy nonstructural element concentration, [g d-2]
    CanopyNodulNonstElms_pft    => plt_biom%CanopyNodulNonstElms_pft   ,& !inoput :canopy nodule nonstructural element, [g d-2]
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft        ,& !inoput :soil layer at planting depth, [-]
    CanopyNoduleNonstCConc_pft  => plt_biom%CanopyNoduleNonstCConc_pft ,& !output :canopy nodule nonstructural C concentration, [gC (gC)-1]
    LeafPetoNonstElmConc_brch   => plt_biom%LeafPetoNonstElmConc_brch  ,& !output :branch nonstructural C concentration, [g d-2]
    CanopyNonstElmConc_pft      => plt_biom%CanopyNonstElmConc_pft     ,& !output :canopy nonstructural element concentration, [g d-2]
    RootNonstructElmConc_rpvr   => plt_biom%RootNonstructElmConc_rpvr  ,& !output :root layer nonstructural C concentration, [g g-1]
    MaxSoiL4Root_pft            => plt_morph%MaxSoiL4Root_pft           & !output :maximum soil layer number for all root axes,[-]
  )
  call PrintInfo('beg '//subname)
  plt_bgcr%RootGasLossDisturb_pft(idg_beg:idg_NH3,NZ)=0.0_r8
  CanopyNonstElms_pft(1:NumPlantChemElms,NZ)=0.0_r8
  MaxSoiL4Root_pft(NZ)   = NMaxRootBotLayer_pft(NZ)
  NGTopRootLayer_pft(NZ) = MIN(MaxSoiL4Root_pft(NZ),MAX(NGTopRootLayer_pft(NZ),NU))
  
!
! TOTAL PLANT NON-STRUCTURAL C, N, P
!
! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
! MainBranchNum_pft=main branch number
!
  D140: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      DO NE=1,NumPlantChemElms      
        CanopyNonstElms_pft(NE,NZ)=CanopyNonstElms_pft(NE,NZ)+CanopyNonstElms_brch(NE,NB,NZ)
        CanopyNodulNonstElms_pft(NE,NZ)=CanopyNodulNonstElms_pft(NE,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)    
      ENDDO     
    ENDIF
  ENDDO D140
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN ROOT
!
! WTRTL=root mas(g)
! CPOOLR,ZPOOLR,PPOOLR=non-structl C,N,P in root(1),myco(2)(g)
! CCPOLR,CZPOLR,CPPOLR=non-structl C,N,P concn in root(1),myco(2)(g g-1)
!
  D180: DO N=1,Myco_pft(NZ)
    D160: DO L=NU,MaxSoiL4Root_pft(NZ)
      IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        DO NE=1,NumPlantChemElms
          RootNonstructElmConc_rpvr(NE,N,L,NZ)=AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ) &
            /RootMycoActiveBiomC_pvr(N,L,NZ))
        ENDDO
      ELSE
        DO NE=1,NumPlantChemElms
          RootNonstructElmConc_rpvr(NE,N,L,NZ)=1.0_r8
        ENDDO
      ENDIF
    ENDDO D160
  ENDDO D180
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN SHOOT
!
! CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy(g g-1)
! CanopyNoduleNonstCConc_pft=nonstructural C concentration in canopy nodules
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!
  IF(CanopyLeafSheathC_pft(NZ).GT.ZERO4LeafVar_pft(NZ))THEN
    DO NE=1,NumPlantChemElms
      CanopyNonstElmConc_pft(NE,NZ)=AZMAX1(AMIN1(1.0_r8,CanopyNonstElms_pft(NE,NZ)/CanopyLeafSheathC_pft(NZ)))
    ENDDO
    CanopyNoduleNonstCConc_pft(NZ)=AZMAX1(AMIN1(1.0_r8,CanopyNodulNonstElms_pft(ielmc,NZ)/CanopyLeafSheathC_pft(NZ)))
  ELSE
    CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ) = 1.0_r8
    CanopyNoduleNonstCConc_pft(NZ)                = 1.0_r8
  ENDIF
  
  D190: DO NB=1,NumOfBranches_pft(NZ)
    IF(CanopyLeafSheathC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      DO NE=1,NumPlantChemElms      
        LeafPetoNonstElmConc_brch(NE,NB,NZ)=AZMAX1(CanopyNonstElms_brch(NE,NB,NZ)/CanopyLeafSheathC_brch(NB,NZ))
      ENDDO  
    ELSE
      LeafPetoNonstElmConc_brch(1:NumPlantChemElms,NB,NZ)=1.0_r8
    ENDIF    
  ENDDO D190
  call PrintInfo('end '//subname)
  end associate
  end subroutine StagePlantPhenology

!----------------------------------------------------------------------------------------------------
  subroutine TestPlantEmergence(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8):: ShootArea
  logical :: CanopyChk,RootChk

  associate(                                               &
    MainBranchNum_pft   => plt_morph%MainBranchNum_pft    ,& !input  :number of main branch,[-]
    CanopyLeafArea_pft  => plt_morph%CanopyLeafArea_pft   ,& !input  :plant canopy leaf area, [m2 d-2]
    ShootElms_pft  => plt_biom%ShootElms_pft    ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    HypoctoHeight_pft   => plt_morph%HypoctoHeight_pft    ,& !input  :cotyledon height, [m]
    Root1stDepz_pft     => plt_morph%Root1stDepz_pft      ,& !input  :root layer depth, [m]
    ZERO4LeafVar_pft    => plt_biom%ZERO4LeafVar_pft      ,& !input  :threshold zero for leaf calculation, [-]
    WatHeldOnCanopy_pft => plt_ew%WatHeldOnCanopy_pft     ,& !input  :canopy surface water content, [m3 d-2]
    CanopyStemArea_pft  => plt_morph%CanopyStemArea_pft   ,& !input  :plant stem area, [m2 d-2]
    SeedDepth_pft       => plt_morph%SeedDepth_pft        ,& !input  :seeding depth, [m]
    iPlantCalendar_brch => plt_pheno%iPlantCalendar_brch  ,& !input  :plant growth stage, [-]
    VHeatCapCanopy_pft  => plt_ew%VHeatCapCanopy_pft       & !output :canopy heat capacity, [MJ d-2 K-1]
  )
!
! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
!
! iPlantCalendar_brch(ipltcal_Emerge,=emergence date
! CanopyLeafArea_pft,CanopyStemArea_pft=leaf,stalk areas
! HypoctoHeight_pft=hypocotyledon height
! SeedDepth_pft=seeding depth
! Root1stDepz_pft=primary root depth,avoid plant the seed at the grid interface
! VHeatCapCanopy_pft,WTSHT,WatHeldOnCanopy_pft=canopy heat capacity,mass,water content
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
    ShootArea = CanopyLeafArea_pft(NZ)+CanopyStemArea_pft(NZ)
    CanopyChk = (HypoctoHeight_pft(NZ).GT.SeedDepth_pft(NZ)).AND.(ShootArea.GT.ZERO4LeafVar_pft(NZ))
    RootChk   = Root1stDepz_pft(ipltroot,1,NZ).GT.(SeedDepth_pft(NZ)+ppmc)
!    write(666,*)I+J/24.,NZ,(HypoctoHeight_pft(NZ).GT.SeedDepth_pft(NZ)),(ShootArea.GT.ZERO4LeafVar_pft(NZ)),&
!      RootChk,CanopyLeafArea_pft(NZ),CanopyStemArea_pft(NZ),HypoctoHeight_pft(NZ),SeedDepth_pft(NZ),Root1stDepz_pft(ipltroot,1,NZ)
    IF(CanopyChk .AND. RootChk)THEN
      iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ)=I
      VHeatCapCanopy_pft(NZ)=cpw*(ShootElms_pft(ielmc,NZ)*10.0E-06_r8+WatHeldOnCanopy_pft(NZ))
    ENDIF
  ENDIF
  end associate
  end subroutine TestPlantEmergence

!----------------------------------------------------------------------------------------------------
  subroutine branch_specific_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J,NB,NZ

  associate(                                                                   &
    PSICanopy_pft                 => plt_ew%PSICanopy_pft                     ,& !input  :canopy total water potential, [Mpa]
    DayLenthCurrent               => plt_site%DayLenthCurrent                 ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                  => plt_site%DayLenthPrev                    ,& !input  :daylength of previous day, [h]
    iPlantCalendar_brch           => plt_pheno%iPlantCalendar_brch            ,& !input  :plant growth stage, [-]
    iPlantRootProfile_pft         => plt_pheno%iPlantRootProfile_pft          ,& !input  :plant growth type (vascular, non-vascular),[-]
    HourReq4LeafOut_brch          => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    iPlantPhenolType_pft          => plt_pheno%iPlantPhenolType_pft           ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    doInitPlant_pft               => plt_pheno%doInitPlant_pft                ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    iPlantBranchState_brch        => plt_pheno%iPlantBranchState_brch         ,& !input  :flag to detect branch death, [-]
    HourReq4LeafOff_brch          => plt_pheno%HourReq4LeafOff_brch           ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    doPlantLeaveOff_brch          => plt_pheno%doPlantLeaveOff_brch           ,& !input  :branch phenology flag, [-]
    doPlantLeafOut_brch           => plt_pheno%doPlantLeafOut_brch            ,& !inoput :branch phenology flag, [-]
    Hours4ShortenPhotoPeriod_brch => plt_pheno%Hours4ShortenPhotoPeriod_brch  ,& !inoput :initial cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch            => plt_pheno%Hours4Leafout_brch             ,& !inoput :heat requirement for spring leafout/dehardening, [h]
    Hours4LenthenPhotoPeriod_brch => plt_pheno%Hours4LenthenPhotoPeriod_brch  ,& !inoput :initial heat requirement for spring leafout/dehardening, [h]
    Hours4LeafOff_brch            => plt_pheno%Hours4LeafOff_brch              & !inoput :cold requirement for autumn leafoff/hardening, [h]
  )

!               PHENOLOGY
!
!           DayLenthPrev,DLYN=daylength of previous,current day
!           Hours4LenthenPhotoPeriod_brch,Hours4ShortenPhotoPeriod_brch
!            =hourly counter for lengthening,shortening photoperiods
!

  IF(iPlantBranchState_brch(NB,NZ).EQ.iLive .OR. doInitPlant_pft(NZ).EQ.itrue)THEN
    IF(DayLenthCurrent.GE.DayLenthPrev)THEN
      Hours4LenthenPhotoPeriod_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)+1.0_r8
      Hours4ShortenPhotoPeriod_brch(NB,NZ)=0.0_r8
    ELSE
      Hours4LenthenPhotoPeriod_brch(NB,NZ)=0.0_r8
      Hours4ShortenPhotoPeriod_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)+1.0_r8
    ENDIF
!
  !
    IF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen)THEN

      call CropEvergreenPhenology(I,J,NB,NZ)
  !
    ELSEIF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid)THEN

      call ColdDeciduousBranchPhenology(I,J,NB,NZ)

  !
  !     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS
  !     ABOVE SPECIFIED WATER POTENTIAL DURINGTopRootLayer_pftDORMANCY
  !
    ELSEIF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu &
      .OR. iPlantPhenolType_pft(NZ).EQ.4 &
      .OR. iPlantPhenolType_pft(NZ).EQ.5)THEN
      !type 4 and 5 are place holders for subtropical/tropical evergreen
      !

      IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable)THEN
        IF(PSICanopy_pft(NZ).GE.PSIMin4LeafOut)THEN
          Hours4Leafout_brch(NB,NZ)=Hours4Leafout_brch(NB,NZ)+1.0_r8
        ENDIF

        IF(Hours4Leafout_brch(NB,NZ).LT.HourReq4LeafOut_brch(NB,NZ))THEN
          IF(PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iPlantRootProfile_pft(NZ)))THEN
            Hours4Leafout_brch(NB,NZ)=AZMAX1(Hours4Leafout_brch(NB,NZ)-12.0_r8)
          ENDIF
        ENDIF

        IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
          Hours4LeafOff_brch(NB,NZ)=0.0_r8
          IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
        ENDIF
      ENDIF

      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
  !
  !     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS
  !     BELOW SPECIFIED WATER POTENTIAL DURINGTopRootLayer_pftGROWINGTopRootLayer_pftSEASON
  !
  !     Hours4Leafout_brch=leafout hours,hours required for leafout
  !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
  !     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
  !     PSICanopy_pft=canopy total water potential
  !     PSIMin4LeafOut,PSIMin4LeafOff=minimum canopy water potential for leafout,leafoff
  !     ALAT=latitude
  !     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
  !     Hours4LenthenPhotoPeriod_brch,Hours4ShortenPhotoPeriod_brch=hourly counter for lengthening,shortening photoperiods
  !     MaxHour4LeafOutOff=maximum hours for leafout,leafoff
  !
      IF(doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. doPlantLeaveOff_brch(NB,NZ).EQ.iEnable)THEN
        IF(PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iPlantRootProfile_pft(NZ)))THEN
          Hours4LeafOff_brch(NB,NZ)=Hours4LeafOff_brch(NB,NZ)+1.0_r8
        ENDIF

        IF(iPlantPhenolType_pft(NZ).EQ.4)THEN
          IF(Hours4ShortenPhotoPeriod_brch(NB,NZ).GT.MaxHour4LeafOutOff)THEN
            Hours4LeafOff_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)
          ENDIF
        ELSEIF(iPlantPhenolType_pft(NZ).EQ.5)THEN
          IF(Hours4LenthenPhotoPeriod_brch(NB,NZ).GT.MaxHour4LeafOutOff)THEN
            Hours4LeafOff_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)
          ENDIF
        ENDIF

        IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ) .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable)THEN
          Hours4Leafout_brch(NB,NZ)  = 0.0_r8
          doPlantLeafOut_brch(NB,NZ) = iEnable
        ENDIF
      ENDIF
  
    ELSEIF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid)THEN
      call ColdDroughtDeciduPhenology(I,J,NB,NZ)

    ENDIF
  ENDIF
  end associate
  end subroutine branch_specific_phenology

!----------------------------------------------------------------------------------------------------
  subroutine ColdDroughtDeciduPhenology(I,J,NB,NZ)

!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS ABOVE SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     LENGTHENINGTopRootLayer_pftPHOTOPERIODS
  implicit none
  integer, intent(in) :: I   !day
  integer, intent(in) :: J   !hour
  integer, intent(in) :: NB  !branch id
  integer, intent(in) :: NZ  !pft id

  associate(                                                   &
    DayLenthCurrent       => plt_site%DayLenthCurrent         ,& !input  :current daylength of the grid, [h]
    DayLenthMax_col       => plt_site%DayLenthMax_col         ,& !input  :maximum daylength of the grid, [h]
    PSICanopyTurg_pft     => plt_ew%PSICanopyTurg_pft         ,& !input  :plant canopy turgor water potential, [MPa]
    DayLenthPrev          => plt_site%DayLenthPrev            ,& !input  :daylength of previous day, [h]
    HourReq4LeafOut_brch  => plt_pheno%HourReq4LeafOut_brch   ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    iPlantCalendar_brch   => plt_pheno%iPlantCalendar_brch    ,& !input  :plant growth stage, [-]
    doPlantLeaveOff_brch  => plt_pheno%doPlantLeaveOff_brch   ,& !input  :branch phenology flag, [-]
    TCGroth_pft           => plt_pheno%TCGroth_pft            ,& !input  :canopy growth temperature, [oC]
    TC4LeafOff_pft        => plt_pheno%TC4LeafOff_pft         ,& !input  :threshold temperature for autumn leafoff/hardening, [oC]
    TCChill4Seed_pft      => plt_pheno%TCChill4Seed_pft       ,& !input  :temperature below which seed set is adversely affected, [oC]
    TC4LeafOut_pft        => plt_pheno%TC4LeafOut_pft         ,& !input  :threshold temperature for spring leafout/dehardening, [oC]
    PSICanopy_pft         => plt_ew%PSICanopy_pft             ,& !input  :canopy total water potential, [Mpa]
    HourReq4LeafOff_brch  => plt_pheno%HourReq4LeafOff_brch   ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantRootProfile_pft => plt_pheno%iPlantRootProfile_pft  ,& !input  :plant growth type (vascular, non-vascular),[-]
    doPlantLeafOut_brch   => plt_pheno%doPlantLeafOut_brch    ,& !inoput :branch phenology flag, [-]
    Hours4LeafOff_brch    => plt_pheno%Hours4LeafOff_brch     ,& !inoput :cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch    => plt_pheno%Hours4Leafout_brch      & !inoput :heat requirement for spring leafout/dehardening, [h]
  )

  IF((DayLenthCurrent.GE.DayLenthPrev .OR. DayLenthCurrent.GE.DayLenthMax_col-2.0_r8) &
    .AND. doPlantLeafOut_brch(NB,NZ).EQ.iEnable)THEN
    IF(TCGroth_pft(NZ).GE.TC4LeafOut_pft(NZ) .AND. PSICanopyTurg_pft(NZ).GT.PSIMin4LeafExpansion)THEN
      Hours4Leafout_brch(NB,NZ)=Hours4Leafout_brch(NB,NZ)+1.0_r8
    ENDIF
    IF(Hours4Leafout_brch(NB,NZ).LT.HourReq4LeafOut_brch(NB,NZ))THEN
      IF(TCGroth_pft(NZ).LT.TCChill4Seed_pft(NZ) .OR. PSICanopyTurg_pft(NZ).LT.PSIMin4LeafExpansion)THEN
        Hours4Leafout_brch(NB,NZ)=AZMAX1(Hours4Leafout_brch(NB,NZ)-1.5_r8)
      ENDIF
    ENDIF
    IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
      Hours4LeafOff_brch(NB,NZ)=0.0_r8
      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
    ENDIF
  ENDIF
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS BELOW SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     SHORTENINGTopRootLayer_pftPHOTOPERIODS
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCGroth_pft,TCZ,TCChill4Seed_pft=canopy temp,leafout threshold temp,chilling temp
!     PSICanopy_pft=canopy total water potential
!     PSIMin4LeafOut,PSIMin4LeafOff=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
  IF((DayLenthCurrent.LT.DayLenthPrev .OR. DayLenthCurrent.LT.24.0_r8-DayLenthMax_col+2.0_r8) &
    .AND. doPlantLeaveOff_brch(NB,NZ).EQ.iEnable)THEN
    IF(TCGroth_pft(NZ).LE.TC4LeafOff_pft(NZ) .OR. &
      PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iPlantRootProfile_pft(NZ)))THEN
      Hours4LeafOff_brch(NB,NZ)=Hours4LeafOff_brch(NB,NZ)+1.0_r8
    ENDIF
    IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ) &
      .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable)THEN
      Hours4Leafout_brch(NB,NZ)  = 0.0_r8
      doPlantLeafOut_brch(NB,NZ) = iEnable
    ENDIF
  ENDIF
  end associate
  end subroutine ColdDroughtDeciduPhenology

!----------------------------------------------------------------------------------------------------
  subroutine CropEvergreenPhenology(I,J,NB,NZ)

  ! CALCULATE EVERGREEN PHENOLOGY DURING TopRootLayer_pft LENGTHENING TopRootLayer_pft PHOTOPERIODS
  implicit none
  integer, intent(in) :: I   !day
  integer, intent(in) :: J   !hour
  integer, intent(in) :: NB  !branch id
  integer, intent(in) :: NZ  !pft id

  associate(                                                                   &
    DayLenthCurrent               => plt_site%DayLenthCurrent                 ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                  => plt_site%DayLenthPrev                    ,& !input  :daylength of previous day, [h]
    HourReq4LeafOff_brch          => plt_pheno%HourReq4LeafOff_brch           ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    Hours4LenthenPhotoPeriod_brch => plt_pheno%Hours4LenthenPhotoPeriod_brch  ,& !input  :initial heat requirement for spring leafout/dehardening, [h]
    HourReq4LeafOut_brch          => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    Hours4ShortenPhotoPeriod_brch => plt_pheno%Hours4ShortenPhotoPeriod_brch  ,& !input  :initial cold requirement for autumn leafoff/hardening, [h]
    ALAT                          => plt_site%ALAT                            ,& !input  :latitude, [degrees north]
    Hours4Leafout_brch            => plt_pheno%Hours4Leafout_brch             ,& !output :heat requirement for spring leafout/dehardening, [h]
    Hours4LeafOff_brch            => plt_pheno%Hours4LeafOff_brch             ,& !output :cold requirement for autumn leafoff/hardening, [h]
    doPlantLeafOut_brch           => plt_pheno%doPlantLeafOut_brch            ,& !output :branch phenology flag, [-]
    doPlantLeaveOff_brch          => plt_pheno%doPlantLeaveOff_brch            & !output :branch phenology flag, [-]
  )
  IF(DayLenthCurrent.GE.DayLenthPrev)THEN
    Hours4Leafout_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)
    IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ) &
      .OR.(ALAT.GT.0.0_r8 .AND. I.EQ.173) &     !northern hemisphere
      .OR.(ALAT.LT.0.0_r8 .AND. I.EQ.355))THEN  !southern hemisphere
      Hours4LeafOff_brch(NB,NZ)   = 0.0_r8
      doPlantLeaveOff_brch(NB,NZ) = iEnable
    ENDIF
  ENDIF
  !
  !   CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayer_pftSHORTENINGTopRootLayer_pftPHOTOPERIODS
  !
  !   Hours4Leafout_brch,Hours4LeafOff_brch=leafout,leafoff hours
  !   Hours4ShortenPhotoPeriod_brch=hourly counter for shortening photoperiods
  !   doPlantLeafOut_brch=flag for enabling leafout:0=enable,1=disable
  !   ALAT=latitude
  !
  IF(DayLenthCurrent.LT.DayLenthPrev)THEN
    Hours4LeafOff_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)
    IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ) &
      .OR.(ALAT.GT.0.0_r8 .AND. I.EQ.355) &     !northern hemisphere
      .OR.(ALAT.LT.0.0_r8 .AND. I.EQ.173))THEN  !southern hemisphere
      Hours4Leafout_brch(NB,NZ)  = 0.0_r8
      doPlantLeafOut_brch(NB,NZ) = iEnable
    ENDIF
  ENDIF
  end associate  
  end subroutine CropEvergreenPhenology

!----------------------------------------------------------------------------------------------------
  subroutine ColdDeciduousBranchPhenology(I,J,NB,NZ)
  !   CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS ABOVE
  !   SPECIFIED TEMPERATURE DURINGTopRootLayer_pftLENGTHENINGTopRootLayer_pftPHOTOPERIODS

  implicit none
  integer, intent(in) :: I   !day
  integer, intent(in) :: J   !hour
  integer, intent(in) :: NB  !branch id
  integer, intent(in) :: NZ  !pft id
  associate(                                                 &
    DayLenthPrev         => plt_site%DayLenthPrev           ,& !input  :daylength of previous day, [h]
    DayLenthCurrent      => plt_site%DayLenthCurrent        ,& !input  :current daylength of the grid, [h]
    HourReq4LeafOff_brch => plt_pheno%HourReq4LeafOff_brch  ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch  => plt_pheno%iPlantCalendar_brch   ,& !input  :plant growth stage, [-]
    TCGroth_pft          => plt_pheno%TCGroth_pft           ,& !input  :canopy growth temperature, [oC]
    TC4LeafOff_pft       => plt_pheno%TC4LeafOff_pft        ,& !input  :threshold temperature for autumn leafoff/hardening, [oC]
    TC4LeafOut_pft       => plt_pheno%TC4LeafOut_pft        ,& !input  :threshold temperature for spring leafout/dehardening, [oC]
    TCChill4Seed_pft     => plt_pheno%TCChill4Seed_pft      ,& !input  :temperature below which seed set is adversely affected, [oC]
    HourReq4LeafOut_brch => plt_pheno%HourReq4LeafOut_brch  ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    ALAT                 => plt_site%ALAT                   ,& !input  :latitude, [degrees north]
    doPlantLeafOut_brch  => plt_pheno%doPlantLeafOut_brch   ,& !inoput :branch phenology flag, [-]
    Hours4LeafOff_brch   => plt_pheno%Hours4LeafOff_brch    ,& !inoput :cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch   => plt_pheno%Hours4Leafout_brch    ,& !inoput :heat requirement for spring leafout/dehardening, [h]
    doPlantLeaveOff_brch => plt_pheno%doPlantLeaveOff_brch   & !output :branch phenology flag, [-]
  )

  IF((DayLenthCurrent.GE.DayLenthPrev .OR. (DayLenthCurrent.LT.DayLenthPrev .AND. &
    Hours4LeafOff_brch(NB,NZ).LT.HourReq4LeafOff_brch(NB,NZ))) &
    .AND.doPlantLeafOut_brch(NB,NZ).EQ.iEnable)THEN
    !favorable temperature for leafout
    IF(TCGroth_pft(NZ).GE.TC4LeafOut_pft(NZ))THEN
      Hours4Leafout_brch(NB,NZ)=Hours4Leafout_brch(NB,NZ)+1.0_r8
    ENDIF
    !not sufficient leafout hour
    IF(Hours4Leafout_brch(NB,NZ).LT.HourReq4LeafOut_brch(NB,NZ))THEN
      !unfavorable cold temperature
      IF(TCGroth_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
        Hours4Leafout_brch(NB,NZ)=AZMAX1(Hours4Leafout_brch(NB,NZ)-1.0_r8)
      ENDIF
    ENDIF
    IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ) &
      .OR.(ALAT.GT.0.0_r8 .AND.I.EQ.173) &     !northern hemisphere
      .OR.(ALAT.LT.0.0_r8 .AND.I.EQ.355))THEN  !southern hemisphere
      Hours4LeafOff_brch(NB,NZ)=0.0_r8
    ENDIF
  ENDIF

  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0 .OR. &
    (DayLenthCurrent.LT.DayLenthPrev .AND. DayLenthCurrent.LT.12.0_r8))THEN
    !has flower or (day length is decreasing, and current daytime is short)
    doPlantLeaveOff_brch(NB,NZ)=iEnable
  ENDIF
!
!     CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS BELOW
!     SPECIFIED TEMPERATURE DURINGTopRootLayer_pftSHORTENINGTopRootLayer_pftPHOTOPERIODS
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCGroth_pft,TCZ,TCChill4Seed_pft=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
  IF(DayLenthCurrent.LT.DayLenthPrev .AND. doPlantLeaveOff_brch(NB,NZ).EQ.iEnable &
    .AND.iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
    IF(TCGroth_pft(NZ).LE.TC4LeafOff_pft(NZ))THEN
      Hours4LeafOff_brch(NB,NZ)=Hours4LeafOff_brch(NB,NZ)+1.0_r8
    ENDIF
    IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ).AND.&
      doPlantLeafOut_brch(NB,NZ).EQ.iDisable)THEN
      Hours4Leafout_brch(NB,NZ)  = 0.0_r8
      doPlantLeafOut_brch(NB,NZ) = iEnable
    ENDIF
  ENDIF
  end associate
  end subroutine ColdDeciduousBranchPhenology

!----------------------------------------------------------------------------------------------------
  subroutine live_branch_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB  !plant branch id
  integer, intent(in) :: NZ  !plant species id
  character(len=*), parameter :: subname='live_branch_phenology'
  real(r8) :: TFNP,WFNG
  real(r8) :: ACTV,OFNG
  real(r8) :: RTK
  real(r8) :: STK,TKCO
  real(r8) :: NodeInitRate,LeafAppearRate,HourlyNodeNumNormByMatgrp_brch
  integer :: kk
  logical :: NodeNumChk,PlantDayChk,LeafOutChk,LeafOffChk,CalChk
  logical :: DayLenChk,CanopyHeightChk,PhenoChk1,PhenoChk2,PhotoPrdChk
! begin_execution
  associate(                                                                           &
    PSICanopy_pft                     => plt_ew%PSICanopy_pft                         ,& !input  :canopy total water potential, [Mpa]
    DayLenthCurrent                   => plt_site%DayLenthCurrent                     ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                      => plt_site%DayLenthPrev                        ,& !input  :daylength of previous day, [h]
    TKGroth_pft                       => plt_pheno%TKGroth_pft                        ,& !input  :canopy growth temperature, [K]
    RefLeafAppearRate_pft             => plt_pheno%RefLeafAppearRate_pft              ,& !input  :rate of leaf initiation, [h-1 at 25 oC]
    TempOffset_pft                    => plt_pheno%TempOffset_pft                     ,& !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft            ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft               ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                    ,& !input  :acclimated plant maturity group, [-]
    PlantO2Stress_pft                 => plt_pheno%PlantO2Stress_pft                  ,& !input  :plant O2 stress indicator, [-]
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch               ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    RefNodeInitRate_pft               => plt_pheno%RefNodeInitRate_pft                ,& !input  :rate of node initiation, [h-1 at 25 oC]
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch            ,& !input  :shoot node number at floral initiation, [-]
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch          ,& !input  :shoot node number at anthesis, [-]
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                  ,& !input  :id number of main branch,[-]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                ,& !inoput :plant growth stage, [-]
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch                   ,& !inoput :leaf number, [-]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !inoput :normalized node number during reproductive growth stages, [-]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !inoput :normalized node number during vegetative growth stages, [-]
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch                 ,& !inoput :heat requirement for spring leafout/dehardening, [h]
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                  ,& !inoput :shoot node number, [-]
    ReprodNodeNumNormByMatrgrp_brch   => plt_pheno%ReprodNodeNumNormByMatrgrp_brch    ,& !output :normalized node number during reproductive growth stages, [-]
    doPlantLeafOut_brch               => plt_pheno%doPlantLeafOut_brch                ,& !output :branch phenology flag, [-]
    NodeNumNormByMatgrp_brch          => plt_pheno%NodeNumNormByMatgrp_brch           ,& !output :normalized node number during vegetative growth stages, [-]
    dReproNodeNumNormByMatG_brch      => plt_pheno%dReproNodeNumNormByMatG_brch       ,& !output :gain in normalized node number during reproductive growth stages, [h-1]
    doSenescence_brch                 => plt_pheno%doSenescence_brch                  ,& !output :branch phenology flag, [-]
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch                  & !output :branch phenology flag, [-]
  )
  call PrintInfo('beg '//subname)
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).EQ.0)THEN
    !plant emergence
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ) = I
    doInitLeafOut_brch(NB,NZ)                 = iDisable
    doPlantLeafOut_brch(NB,NZ)                = iEnable
    Hours4Leafout_brch(NB,NZ)                 = 0.5_r8*Hours4Leafout_brch(MainBranchNum_pft(NZ),NZ)
  ENDIF
!
! CALCULATE NODE INITIATION AND LEAF APPEARANCE RATES
! FROM TEMPERATURE FUNCTION CALCULATED IN 'UPTAKE' AND
! RATES AT 25C ENTERED IN 'READQ' EXCEPT WHEN DORMANT
!
! iPlantPhenolType_pft=phenology type from PFT file, 0=evergreen,1=cold deciduous, 2=drought deciduous,3=1+2
! Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
! TKGroth_pft,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
! TempOffset_pft=shift in Arrhenius curve for thermal adaptation
! TFNP=temperature function for phenology (25 oC =1 )
! 8.313,710.0=gas constant,enthalpy
! 60000,197500,218500=energy of activn,high,low temp inactivn(KJ mol-1)
! NodeInitRate,LeafAppearRate=rates of node initiation,leaf appearance
! XRNI,XRLA=rate of node initiation,leaf appearance at 25 oC (h-1)
!
  call DebugPrint('TKGroth_pft(NZ)',TKGroth_pft(NZ))

  IF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen .OR. Hours4LeafOff_brch(NB,NZ).LT.HourReq4LeafOff_brch(NB,NZ))THEN
    TKCO           = TKGroth_pft(NZ)+TempOffset_pft(NZ)
    TFNP           = calc_leave_grow_tempf(TKCO)
    call DebugPrint('TFNP',TFNP)

    NodeInitRate   = AZMAX1(RefNodeInitRate_pft(NZ)*TFNP)
    LeafAppearRate = AZMAX1(RefLeafAppearRate_pft(NZ)*TFNP)
!
!   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
!
!   PSICanopyTurg_pft=leaf turgor potential
!   WFNG=water stress effect on phenology
!   only annual plants depends on moisture
    call DebugPrint('LeafAppearRate',LeafAppearRate)
    call DebugPrint('PSICanopy_pft(NZ)',PSICanopy_pft(NZ))

    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
      IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
        WFNG           = EXP(0.025_r8*AMAX1(PSICanopy_pft(NZ),-1000._r8))
        NodeInitRate   = NodeInitRate*WFNG
        LeafAppearRate = LeafAppearRate*WFNG
      ENDIF

      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
        OFNG           = SQRT(PlantO2Stress_pft(NZ))
        NodeInitRate   = NodeInitRate*OFNG
        LeafAppearRate = LeafAppearRate*OFNG
      ENDIF
!    ELSE 
      !perennial plants tend to maintain stronger antioxidant defenses and show delayed senescence compared to annuals    
!      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
!        OFNG           = PlantO2Stress_pft(NZ)**0.5_r8
!        NodeInitRate   = NodeInitRate*OFNG
!        LeafAppearRate = LeafAppearRate*OFNG
!      ENDIF      
    ENDIF
    call DebugPrint('NodeInitRate',NodeInitRate)
!
!   ACCUMULATE NODE INITIATION AND LEAF APPEARANCE RATES
!   INTO TOTAL NUMBER OF NODES AND LEAVES
!
!   PSTG,NumOfLeaves_brch=number of nodes initiated,leaves appeared

    ShootNodeNum_brch(NB,NZ) = ShootNodeNum_brch(NB,NZ)+NodeInitRate
    NumOfLeaves_brch(NB,NZ)  = NumOfLeaves_brch(NB,NZ)+LeafAppearRate
!
!   USE TOTAL NUMBER OF NODES TO CALCULATE PROGRESSION THROUGH
!   VEGETATIVE AND REPRODUCTIVE GROWTH STAGES. THIS PROGRESSION
!   IS USED TO SET START AND END DATES FOR GROWTH STAGES BELOW
!
!   NodeNumNormByMatgrp_brch=vegetative node number normalized for maturity group
!   MatureGroup_pft=node number required for floral initiation

    IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
      NodeNumNormByMatgrp_brch(NB,NZ)      = (ShootNodeNum_brch(NB,NZ)-NodeNum2InitFloral_brch(NB,NZ))/MatureGroup_pft(NZ)
      HourlyNodeNumNormByMatgrp_brch       = NodeInitRate/(MatureGroup_pft(NZ)*GrowStageNorm4VegetaPheno)
      TotalNodeNumNormByMatgrp_brch(NB,NZ) = TotalNodeNumNormByMatgrp_brch(NB,NZ)+HourlyNodeNumNormByMatgrp_brch
    ENDIF

    IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0)THEN
      ReprodNodeNumNormByMatrgrp_brch(NB,NZ)   = (ShootNodeNum_brch(NB,NZ)-NodeNumberAtAnthesis_brch(NB,NZ))/MatureGroup_pft(NZ)
      dReproNodeNumNormByMatG_brch(NB,NZ)      = NodeInitRate/(MatureGroup_pft(NZ)*GrowStageNorm4ReprodPheno)
      TotReproNodeNumNormByMatrgrp_brch(NB,NZ) = TotReproNodeNumNormByMatrgrp_brch(NB,NZ)+dReproNodeNumNormByMatG_brch(NB,NZ)
    ENDIF
    call DebugPrint('TotReproNodeNumNormByMatrgrp_brch(NB,NZ)',TotReproNodeNumNormByMatrgrp_brch(NB,NZ))
    doSenescence_brch(NB,NZ)=itrue
  ELSE
    doSenescence_brch(NB,NZ)=ifalse
  ENDIF
  call PrintInfo('before testing phenolgoy')
!
! REPRODUCTIVE GROWTH STAGES ADVANCE WHEN THRESHOLD NUMBER
! OF NODES HAVE BEEN INITIATED. FIRST DETERMINE PHOTOPERIOD
! AND TEMPERATURE EFFECTS ON FINAL VEG NODE NUMBER FROM
! NUMBER OF INITIATED NODES
!
! PSTG=number of nodes initiated
! NodeNum2InitFloral_brch=node number at floral initiation
! GROUP=node number required for floral initiation
! Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
! DayLenthPrev,DLYN=daylength of previous,current day
! iPlantPhenolPattern_pft=growth habit from PFT file
! iPlantPhenolType_pft=phenology type from PFT file
! ZC,SnowDepth=canopy height,snowpack depth
!
  DayLenChk=DayLenthCurrent.LT.DayLenthPrev

  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN    !2

    call InitiateBranchFlora(I,J,NB,NZ,DayLenChk)
  !
  ELSEIF(iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).EQ.0)THEN  !3

    call BranchStemJointing(I,J,NB,NZ,DayLenChk)
    !
    !   ipltcal_Elongation,=mid stem elongation
    !
  ELSEIF(iPlantCalendar_brch(ipltcal_Elongation,NB,NZ).EQ.0)THEN !4

    call BranchStemElongation(I,J,NB,NZ,DayLenChk)
    !
    !   ipltcal_Heading,=end of stem elongation and setting max seed number
    !
  ELSEIF(iPlantCalendar_brch(ipltcal_Heading,NB,NZ).EQ.0)THEN !5

    call BranchHeading(I,J,NB,NZ,DayLenChk) 
    !!
  ELSEIF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN !6
  
    call BranchAnthesis(I,J,NB,NZ,DayLenChk) 
!
  ELSEIF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN !7

    call InitBranchGrainFill(I,J,NB,NZ,DayLenChk) 
!
!   END SEED NUMBER SET PERIOD
!
!   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
!
  ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN !8
    call PrintInfo('set seed number')
    IF(ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.1.00_r8*GrowStageNorm4ReprodPheno)THEN
      iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ)=I
    ENDIF
!
!   END SEED SIZE SET PERIOD
!
!   iPlantCalendar_brch(ipltcal_SetSeedMass,=end of setting max seed size
!
  ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN  !9
    call PrintInfo('set seed number')
    IF(ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.1.50_r8*GrowStageNorm4ReprodPheno)THEN
      iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ)=I
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine live_branch_phenology

!----------------------------------------------------------------------------------------------------
  subroutine InitBranchGrainFill(I,J,NB,NZ,DayLenChk) 
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLenChk
  character(len=*), parameter :: subname='InitBranchGrainFill'
  logical :: NodeNumChk,LeafOffChk,PhenoChk1,PhenoChk2

  associate(                                                                       &
    ReprodNodeNumNormByMatrgrp_brch => plt_pheno%ReprodNodeNumNormByMatrgrp_brch  ,& !input  :normalized node number during reproductive growth stages, [-]
    Hours4LeafOff_brch              => plt_pheno%Hours4LeafOff_brch               ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft         => plt_pheno%iPlantPhenolPattern_pft          ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantPhenolType_pft            => plt_pheno%iPlantPhenolType_pft             ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    HourReq4LeafOff_brch            => plt_pheno%HourReq4LeafOff_brch             ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    doPlantLeafOut_brch             => plt_pheno%doPlantLeafOut_brch              ,& !input  :branch phenology flag, [-]
    iPlantPhotoperiodType_pft       => plt_pheno%iPlantPhotoperiodType_pft        ,& !input  :photoperiod type (neutral, long day, short day),[-]
    iPlantCalendar_brch             => plt_pheno%iPlantCalendar_brch               & !output :plant growth stage, [-]
  )
  call PrintInfo('beg '//subname)
  NodeNumChk = ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.0.50_r8*GrowStageNorm4ReprodPheno
  LeafOffChk = Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)
  PhenoChk1  = (iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid) &
    .AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short .AND. DayLenChk
  PhenoChk2 = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual

  IF(NodeNumChk .OR.(PhenoChk1.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable.AND.LeafOffChk) &
                .OR. PhenoChk2.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable.AND.LeafOffChk)THEN
      iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ)=I
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitBranchGrainFill

!----------------------------------------------------------------------------------------------------
  subroutine BranchAnthesis(I,J,NB,NZ,DayLenChk) 
!   ANTHESIS OCCURS WHEN THE NUMBER OF LEAVES THAT HAVE APPEARED
!   EQUALS THE NUMBER OF NODES INITIATED WHEN THE FINAL VEGETATIVE
!   NODE NUMBER WAS SET ABOVE

  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLenChk
  character(len=*), parameter :: subname='BranchAnthesis'
  logical :: NodeNumChk,LeafOffChk,PhenoChk1,PhenoChk2,CalChk

  associate(                                                           &
    NumOfLeaves_brch          => plt_morph%NumOfLeaves_brch           ,& !input  :leaf number, [-]
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    NodeNum2InitFloral_brch   => plt_morph%NodeNum2InitFloral_brch    ,& !input  :shoot node number at floral initiation, [-]
    doPlantLeafOut_brch       => plt_pheno%doPlantLeafOut_brch        ,& !input  :branch phenology flag, [-]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft          ,& !input  :id number of main branch,[-]
    ShootNodeNum_brch         => plt_morph%ShootNodeNum_brch          ,& !input  :shoot node number, [-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch        ,& !inoput :plant growth stage, [-]
    NodeNumberAtAnthesis_brch => plt_morph%NodeNumberAtAnthesis_brch   & !output :shoot node number at anthesis, [-]
  )
  call PrintInfo('beg '//subname)
  NodeNumChk = NumOfLeaves_brch(NB,NZ).GT.NodeNum2InitFloral_brch(NB,NZ)
  LeafOffChk = Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)
  PhenoChk1  = (iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid) &
    .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short .AND. DayLenChk
  PhenoChk2 = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual
  CalChk    = iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. iPlantCalendar_brch(ipltcal_Heading,NB,NZ).NE.0


  IF(NodeNumChk .OR. CalChk &
    .OR.(PhenoChk1 .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. LeafOffChk) &
    .OR. PhenoChk2 .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. LeafOffChk)THEN

    IF(NB.EQ.MainBranchNum_pft(NZ) .OR. iPlantCalendar_brch(ipltcal_Anthesis,MainBranchNum_pft(NZ),NZ).NE.0)THEN
      iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ) = I
      NodeNumberAtAnthesis_brch(NB,NZ)            = ShootNodeNum_brch(NB,NZ)
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine BranchAnthesis

!----------------------------------------------------------------------------------------------------
  subroutine BranchHeading(I,J,NB,NZ,DayLenChk)
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLenChk
  character(len=*), parameter :: subname='BranchHeading'
  logical :: NodeNumChk,LeafOffChk,PhenoChk1,PhenoChk2

  associate(                                                           &
    NodeNumNormByMatgrp_brch  => plt_pheno%NodeNumNormByMatgrp_brch   ,& !input  :normalized node number during vegetative growth stages, [-]
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    doPlantLeafOut_brch       => plt_pheno%doPlantLeafOut_brch        ,& !input  :branch phenology flag, [-]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch         & !output :plant growth stage, [-]
  )

  call PrintInfo('beg '//subname)
  NodeNumChk=NodeNumNormByMatgrp_brch(NB,NZ).GT.1.00_r8*GrowStageNorm4VegetaPheno 
  LeafOffChk=Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)
  PhenoChk1=(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid.OR.iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid) &
    .AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhotoperiodType_pft(NZ).NE.iphotop_short.AND.DayLenChk
  PhenoChk2=iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu.AND.iPlantPhenolPattern_pft(NZ).EQ.iplt_annual  
  
  IF(NodeNumChk .OR.(PhenoChk1.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND.LeafOffChk) &
                .OR. PhenoChk2 .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable.AND.LeafOffChk)THEN
    iPlantCalendar_brch(ipltcal_Heading,NB,NZ)=I
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BranchHeading

!----------------------------------------------------------------------------------------------------
  subroutine BranchStemElongation(I,J,NB,NZ,DayLenChk)
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLenChk
  character(len=*), parameter :: subname='BranchStemElongation'
  logical :: NodeNumChk,LeafOffChk,PhenoChk1,PhenoChk2

  associate(                                                               &
    NodeNumNormByMatgrp_brch    => plt_pheno%NodeNumNormByMatgrp_brch     ,& !input  :normalized node number during vegetative growth stages, [-]
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    Hours4LeafOff_brch          => plt_pheno%Hours4LeafOff_brch           ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft      ,& !input  :plant growth habit: annual or perennial,[-]
    doPlantLeafOut_brch         => plt_pheno%doPlantLeafOut_brch          ,& !input  :branch phenology flag, [-]
    ShootNodeNum_brch           => plt_morph%ShootNodeNum_brch            ,& !input  :shoot node number, [-]
    iPlantDevelopPattern_pft    => plt_pheno%iPlantDevelopPattern_pft     ,& !input  :plant growth habit (determinate or indeterminate),[-]
    iPlantPhotoperiodType_pft   => plt_pheno%iPlantPhotoperiodType_pft    ,& !input  :photoperiod type (neutral, long day, short day),[-]
    HourReq4LeafOff_brch        => plt_pheno%HourReq4LeafOff_brch         ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch         => plt_pheno%iPlantCalendar_brch          ,& !output :plant growth stage, [-]
    LeafNumberAtFloralInit_brch => plt_pheno%LeafNumberAtFloralInit_brch   & !output :leaf number at floral initiation, [-]
  )
  call PrintInfo('beg '//subname)
  NodeNumChk = NodeNumNormByMatgrp_brch(NB,NZ).GT.0.50_r8*GrowStageNorm4VegetaPheno
  LeafOffChk = Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)
  PhenoChk1=(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid) &
    .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short .AND. DayLenChk
  PhenoChk2=iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual
  
  IF(NodeNumChk .OR.(PhenoChk1 .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. LeafOffChk) &
                .OR. PhenoChk2 .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. LeafOffChk)THEN
    iPlantCalendar_brch(ipltcal_Elongation,NB,NZ)=I

    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantDevelopPattern_pft(NZ).NE.ideterminate)THEN
      LeafNumberAtFloralInit_brch(NB,NZ)=ShootNodeNum_brch(NB,NZ)
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine BranchStemElongation    

!----------------------------------------------------------------------------------------------------
  subroutine BranchStemJointing(I,J,NB,NZ,DayLenChk)
!   STEM ELONGATION
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLenChk
  character(len=*), parameter :: subname='BranchStemJointing'
  logical :: NodeNumChk,PhenoChk1,PhenoChk2,LeafOffChk

  associate(                                                           &
    NodeNumNormByMatgrp_brch  => plt_pheno%NodeNumNormByMatgrp_brch   ,& !input  :normalized node number during vegetative growth stages, [-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    doPlantLeafOut_brch       => plt_pheno%doPlantLeafOut_brch        ,& !input  :branch phenology flag, [-]
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch         & !output :plant growth stage, [-]
  )
  call PrintInfo('beg '//subname)
  NodeNumChk = NodeNumNormByMatgrp_brch(NB,NZ).GT.0.25_r8*GrowStageNorm4VegetaPheno
  PhenoChk1  = (iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid) &
    .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short .AND. DayLenChk
  PhenoChk2  = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual
  LeafOffChk = Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)
  
  IF(NodeNumChk .OR. (PhenoChk1 .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. LeafOffChk) &
                .OR. PhenoChk2 .AND. doPlantLeafOut_brch(NB,NZ).EQ.iDisable .AND. LeafOffChk)THEN
    
    iPlantCalendar_brch(ipltcal_Jointing,NB,NZ)=I
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BranchStemJointing

!----------------------------------------------------------------------------------------------------
  subroutine InitiateBranchFlora(I,J,NB,NZ,DayLenChk)
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLenChk
  character(len=*), parameter :: subname='InitiateBranchFlora'
  logical :: NodeNumChk,LeafOutChk,PlantDayChk,CanopyHeightChk,PhenoChk1,PhenoChk2,PhotoPrdChk
  real(r8) :: PPD

  associate(                                                               &
    ShootNodeNum_brch           => plt_morph%ShootNodeNum_brch            ,& !input  :shoot node number, [-]
    ZERO                        => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft      ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantDevelopPattern_pft    => plt_pheno%iPlantDevelopPattern_pft     ,& !input  :plant growth habit (determinate or indeterminate),[-]
    iYearCurrent                => plt_site%iYearCurrent                  ,& !input  :current year,[-]
    SnowDepth                   => plt_ew%SnowDepth                       ,& !input  :snowpack depth, [m]
    CanopyHeight_pft            => plt_morph%CanopyHeight_pft             ,& !input  :canopy height, [m]
    iDayPlanting_pft            => plt_distb%iDayPlanting_pft             ,& !input  :day of planting,[-]
    PhotoPeriodSens_pft         => plt_pheno%PhotoPeriodSens_pft          ,& !input  :difference between current and critical daylengths used to calculate phenological progress, [h]
    HourReq4LeafOut_brch        => plt_pheno%HourReq4LeafOut_brch         ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    DayLenthCurrent             => plt_site%DayLenthCurrent               ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                => plt_site%DayLenthPrev                  ,& !input  :daylength of previous day, [h]
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iYearPlanting_pft           => plt_distb%iYearPlanting_pft            ,& !input  :year of planting,[-]
    Hours4Leafout_brch          => plt_pheno%Hours4Leafout_brch           ,& !input  :heat requirement for spring leafout/dehardening, [h]
    CriticPhotoPeriod_pft       => plt_pheno%CriticPhotoPeriod_pft        ,& !input  :critical daylength for phenological progress, [h]
    iPlantPhotoperiodType_pft   => plt_pheno%iPlantPhotoperiodType_pft    ,& !input  :photoperiod type (neutral, long day, short day),[-]
    MatureGroup_brch            => plt_pheno%MatureGroup_brch             ,& !input  :plant maturity group, [-]
    NodeNum2InitFloral_brch     => plt_morph%NodeNum2InitFloral_brch      ,& !inoput :shoot node number at floral initiation, [-]
    LeafNumberAtFloralInit_brch => plt_pheno%LeafNumberAtFloralInit_brch  ,& !output :leaf number at floral initiation, [-]
    iPlantCalendar_brch         => plt_pheno%iPlantCalendar_brch           & !output :plant growth stage, [-]
  )

  call PrintInfo('beg '//subname)
  NodeNumChk      = ShootNodeNum_brch(NB,NZ).GT.MatureGroup_brch(NB,NZ)+NodeNum2InitFloral_brch(NB,NZ)
  LeafOutChk      = Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ)
  PlantDayChk     = I.GE.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ) .AND. DayLenthCurrent.GT.DayLenthPrev
  CanopyHeightChk = CanopyHeight_pft(NZ).GE.SnowDepth-ZERO
  PhenoChk1       = iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. &
    (iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid)

  !Annual crop plants (i.e. those seeded by human) are set as evergreen, if it is self-seeding, then 
  !iPlantPhenolType_pft should be set according to koppen climate zone.  
  PhenoChk2=iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen

  IF((NodeNumChk .AND. (LeafOutChk.OR.PlantDayChk)) .OR. &
    ((PhenoChk1.OR.PhenoChk2) .AND. CanopyHeightChk .AND. DayLenChk))THEN
!
!     FINAL VEGETATIVE NODE NUMBER DEPENDS ON PHOTOPERIOD FROM 'DAY'
!     AND ON MATURITY GROUP, CRITICAL PHOTOPERIOD AND PHOTOPERIOD
!     SENSITIVITY ENTERED IN 'READQ'
!
!     iPlantPhotoperiodType_pft=photoperiod type from PFT file
!     PPD=photoperiod sensitivity
!     CriticPhotoPeriod_pft=critical photoperiod from PFT file
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!     VSTGX=node number on date of floral initiation
!
    IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_neutral)THEN
      PPD=0.0_r8
    ELSE
      PPD=AZMAX1(CriticPhotoPeriod_pft(NZ)-DayLenthCurrent)
      IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_short .AND. DayLenthCurrent.GE.DayLenthPrev)PPD=0.0_r8
    ENDIF
    
    PhotoPrdChk=iPlantPhotoperiodType_pft(NZ).EQ.iphotop_neutral &
      .OR.(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_short.AND. PPD.GT.PhotoPeriodSens_pft(NZ)) &
      .OR.(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_long .AND. PPD.LT.PhotoPeriodSens_pft(NZ))

    IF( PhotoPrdChk .OR. ((PhenoChk1.OR.PhenoChk2) .AND. CanopyHeightChk.AND.DayLenChk))THEN
      iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ) = I
      NodeNum2InitFloral_brch(NB,NZ)                = ShootNodeNum_brch(NB,NZ)
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
        LeafNumberAtFloralInit_brch(NB,NZ)=ShootNodeNum_brch(NB,NZ)
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitiateBranchFlora
  ![tail]
end module PlantPhenolMod
