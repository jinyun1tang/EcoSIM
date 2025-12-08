module LitterFallMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use DebugToolMod,  only: PrintInfo
  use EcosimConst
  use PlantBGCPars
  use PlantAPIData
  use PlantMathFuncMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: ResetDeadPlant
  public :: ReSeedPlants
  contains
  ![header]  
!----------------------------------------------------------------------------------------------------
  subroutine ResetDeadPlant(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NumDeadBranches
  logical :: emerge_check,harvst_check,reset_check

!     begin_execution
  associate(                                                       &
    C4PhotoShootNonstC_brch => plt_biom%C4PhotoShootNonstC_brch   ,& !input  :branch shoot nonstrucal elelment, [g d-2]
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft        ,& !input  :number of main branch,[-]
    SolarNoonHour_col       => plt_site%SolarNoonHour_col         ,& !input  :time of solar noon, [h]
    iDayPlantHarvest_pft    => plt_distb%iDayPlantHarvest_pft     ,& !input  :day of harvest,[-]
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch      ,& !input  :plant growth stage, [-]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    iYearCurrent            => plt_site%iYearCurrent              ,& !input  :current year,[-]
    iYearPlantHarvest_pft   => plt_distb%iYearPlantHarvest_pft     & !input  :year of harvest,[-]
  )
!
!     SolarNoonHour_col=hour of solar noon
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iDayPlantHarvest_pft,iYearPlantHarvest_pft=day,year of harvesting
!     iYearCurrent=current year
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     GROUP=node number required for floral initiation
!     NodeNum2InitFloral_brch=node number at floral initiation
!     NodeNumberAtAnthesis_brch=node number at flowering
!     NumOfLeaves_brch=number of leaves appeared
!     KHiestGroLeafNode_brch=integer of most recent leaf number currently growing
!     VSTGX=leaf number on date of floral initiation
!     TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
!     TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
!     HourFailGrainFill_brch=number of hours with no grain fill
!     doInitLeafOut_brch=flag for initializing leafout
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     Hours2LeafOut_brch=hourly leafout counter
!     RubiscoActivity_brch,C4PhotosynDowreg_brch=N,P feedback inhibition on C3 CO2 fixation
!     doInitLeafOut_brch,doPlantLeafOut_brch=flag for initializing,enabling leafout
!     doPlantLeaveOff_brch=flag for enabling leafoff:0=enable,1=disable
!     Hours4LiterfalAftMature_brch=current hours after physl maturity until start of LitrFall
!
  plt_pheno%doReSeed_pft(NZ)=.FALSE.
  harvst_check = (I.GE.iDayPlantHarvest_pft(NZ) .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ))
  emerge_check = J.EQ.INT(SolarNoonHour_col) .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).GT.0
  reset_check  = emerge_check .AND. (iPlantPhenolPattern_pft(NZ).NE.iplt_annual .OR. harvst_check)

  IF(reset_check)THEN
    
    NumDeadBranches=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!

    call LiterFallDeadBranches(I,J,NZ,NumDeadBranches,C4PhotoShootNonstC_brch)

    call SetDeadPlant(I,J,NZ,NumDeadBranches)
!
!     DEAD ROOTS
!
!     LitrFall FROM DEAD ROOTS
!
    call LiterFallDeadRoots(I,J,NZ)

    call LiterFallRootShootStorage(I,J,NZ,C4PhotoShootNonstC_brch)

  ENDIF

  end associate
  end subroutine ResetDeadPlant

!----------------------------------------------------------------------------------------------------
  subroutine ReSeedPlants(I,J,NZ)
  !
  !Description:
  !Re-Seed plants if there is non-zero storage (sort of seed bank)
  !
  use InitPlantMod, only: InitPlantPhenoMorphoBio
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  logical :: reseed_check
  character(len=*), parameter :: subname='ReSeedPlants'
  associate(                                                  &
    SeasonalNonstElms_pft => plt_biom%SeasonalNonstElms_pft  ,& !input  :plant stored nonstructural element at current step, [g d-2]
    ZERO4Groth_pft        => plt_biom%ZERO4Groth_pft         ,& !input  :threshold zero for plang growth calculation, [-]
    doReSeed_pft          => plt_pheno%doReSeed_pft          ,& !input  :flag to do annual plant reseeding, [-]
    IsPlantActive_pft     => plt_pheno%IsPlantActive_pft      & !output :flag for living pft, [-]

  )
  call PrintInfo('beg '//subname)

  reseed_check=doReSeed_pft(NZ) .and. SeasonalNonstElms_pft(ielmc,NZ) .GT. ZERO4Groth_pft(NZ)
  if(doReSeed_pft(NZ) .and. SeasonalNonstElms_pft(ielmc,NZ) .GT. ZERO4Groth_pft(NZ))then

    IsPlantActive_pft(NZ) = iDormant

    call InitPlantPhenoMorphoBio(NZ)

  endif  
  call PrintInfo('end '//subname)

  end associate
  end subroutine ReSeedPlants

!----------------------------------------------------------------------------------------------------
  subroutine SetDeadPlant(I,J,NZ,NumDeadBranches)
  implicit none
  INTEGER, intent(in) :: i,J
  integer, intent(in) :: NZ
  integer, intent(in) :: NumDeadBranches

  associate(                                                       &
    PlantPopulation_pft     => plt_site%PlantPopulation_pft       ,& !input  :plant population, [d-2]
    RootElms_pft            => plt_biom%RootElms_pft              ,& !input  :plant root element mass, [g d-2]
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft     ,& !input  :plant stored nonstructural element at current step, [g d-2]
    doInitPlant_pft         => plt_pheno%doInitPlant_pft          ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    jHarvstType_pft         => plt_distb%jHarvstType_pft          ,& !input  :flag for stand replacing disturbance,[-]
    CanopyBiomWater_pft     => plt_ew%CanopyBiomWater_pft         ,& !inoput :canopy water content, [m3 d-2]
    H2OLoss_CumYr_col       => plt_ew%H2OLoss_CumYr_col           ,& !inoput :total subsurface water flux, [m3 d-2]
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft        ,& !inoput :number of branches,[-]
    QH2OLoss_lnds           => plt_site%QH2OLoss_lnds             ,& !inoput :total subsurface water loss flux over the landscape, [m3 d-2]
    iPlantRootState_pft     => plt_pheno%iPlantRootState_pft      ,& !output :flag to detect root system death,[-]
    BranchNumber_pft        => plt_morph%BranchNumber_pft         ,& !output :main branch numeric id,[-]
    HoursTooLowPsiCan_pft   => plt_pheno%HoursTooLowPsiCan_pft    ,& !output :canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
    HypoctoHeight_pft       => plt_morph%HypoctoHeight_pft        ,& !output :cotyledon height, [m]
    doReSeed_pft            => plt_pheno%doReSeed_pft             ,& !output :flag to do annual plant reseeding, [-]
    iPlantShootState_pft    => plt_pheno%iPlantShootState_pft      & !output :flag to detect canopy death,[-]
  )

  IF(NumDeadBranches.EQ.NumOfBranches_pft(NZ))THEN
    iPlantShootState_pft(NZ)  = iDead
    BranchNumber_pft(NZ)      = 0
    HoursTooLowPsiCan_pft(NZ) = 0._r8

    IF(doInitPlant_pft(NZ).EQ.itrue)THEN
      NumOfBranches_pft(NZ)=1
    ELSE
      NumOfBranches_pft(NZ)=0
    ENDIF
    HypoctoHeight_pft(NZ)   = 0._r8
    QH2OLoss_lnds           = QH2OLoss_lnds+CanopyBiomWater_pft(NZ)
    H2OLoss_CumYr_col       = H2OLoss_CumYr_col+CanopyBiomWater_pft(NZ)
    CanopyBiomWater_pft(NZ) = 0._r8
    !
    !     RESET LIVING FLAGS
    !
    !     WTRVC,WTRT=PFT storage,root C
    !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
    !     jHarvstType_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
    !     PP=PFT population
    !     iPlantShootState_pft,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
    !
    IF(SeasonalNonstElms_pft(ielmc,NZ).LT.1.0E-04_r8*RootElms_pft(ielmc,NZ) .AND. &
      iPlantPhenolPattern_pft(NZ).NE.iplt_annual)then
      iPlantRootState_pft(NZ)=iDead
    endif
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)then
      iPlantRootState_pft(NZ)=iDead
    endif
    IF(jHarvstType_pft(NZ).NE.jharvtyp_noaction)then
      iPlantRootState_pft(NZ)=iDead
    endif
    IF(PlantPopulation_pft(NZ).LE.0.0_r8)then
      iPlantRootState_pft(NZ)=iDead
    endif
    IF(iPlantRootState_pft(NZ).EQ.iDead)then
      iPlantShootState_pft(NZ)=iDead
    endif
    doReSeed_pft(NZ)=.true.    
  ENDIF
  end associate
  end subroutine SetDeadPlant

!----------------------------------------------------------------------------------------------------
  subroutine LiterFallRootShootStorage(I,J,NZ,C4PhotoShootNonstC_brch)

  use EcoSIMCtrlDataType, only : iYearCurrent
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: C4PhotoShootNonstC_brch(NumCanopyLayers1,JP1)
  integer :: L,M,NR,NB,N,NE
  REAL(R8) :: dRootMyco,dDeadE
!     begin_execution
  associate(                                                               &
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch     ,& !input  :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch   => plt_biom%CanopyNodulStrutElms_brch     ,& !input  :branch nodule structural element, [g d-2]
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch          ,& !input  :branch nonstructural element, [g d-2]
    DazCurrYear                 => plt_site%DazCurrYear                   ,& !input  :number of days in current year,[-]
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch             ,& !input  :branch ear structural chemical element mass, [g d-2]
    PlantElmAllocMat4Litr       => plt_soilchem%PlantElmAllocMat4Litr     ,& !input  :litter kinetic fraction, [-]
    FracRootElmAlloc2Litr       => plt_allom%FracRootElmAlloc2Litr        ,& !input  :C woody fraction in root,[-]
    FracWoodStalkElmAlloc2Litr  => plt_allom%FracWoodStalkElmAlloc2Litr   ,& !input  :woody element allocation,[-]
    FracShootLeafAlloc2Litr     => plt_allom%FracShootLeafAlloc2Litr      ,& !input  :woody element allocation, [-]
    FracShootPetolAlloc2Litr    => plt_allom%FracShootPetolAlloc2Litr     ,& !input  :leaf element allocation,[-]
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch           ,& !input  :branch grain structural element mass, [g d-2]
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch            ,& !input  :branch husk structural element mass, [g d-2]
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch            ,& !input  :branch leaf structural element mass, [g d-2]
    MaxNumRootLays              => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]
    Myco_pft                    => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft           ,& !input  :soil layer at planting depth, [-]
    NU                          => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    NumPrimeRootAxes_pft        => plt_morph%NumPrimeRootAxes_pft         ,& !input  :root primary axis number,[-]
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch          ,& !input  :branch sheath structural element, [g d-2]
    RootMyco1stStrutElms_rpvr   => plt_biom%RootMyco1stStrutElms_rpvr     ,& !input  :root layer element primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr   => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !input  :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr        ,& !input  :root layer nonstructural element, [g d-2]
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch            ,& !input  :branch reserve element mass, [g d-2]
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch           ,& !input  :branch stalk structural element mass, [g d-2]
    doInitPlant_pft             => plt_pheno%doInitPlant_pft              ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft      ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft        ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantRootState_pft         => plt_pheno%iPlantRootState_pft          ,& !input  :flag to detect root system death,[-]
    iPlantShootState_pft        => plt_pheno%iPlantShootState_pft         ,& !input  :flag to detect canopy death,[-]
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft    ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    icwood                      => pltpar%icwood                          ,& !input  :group id of coarse woody litter
    ifoliar                     => pltpar%ifoliar                         ,& !input  :group id of plant foliar litter
    inonfoliar                  => pltpar%inonfoliar                      ,& !input  :group id of plant non-foliar litter group
    inonstruct                  => pltpar%inonstruct                      ,& !input  :group id of plant nonstructural litter
    iroot                       => pltpar%iroot                           ,& !input  :group id of plant root litter
    istalk                      => pltpar%istalk                          ,& !input  :group id of plant stalk litter group
    jHarvstType_pft             => plt_distb%jHarvstType_pft              ,& !input  :flag for stand replacing disturbance,[-]
    k_fine_comp                 => pltpar%k_fine_comp                     ,& !input  :fine litter complex id
    k_woody_comp                => pltpar%k_woody_comp                    ,& !input  :woody litter complex id
    LitrfallElms_pvr        => plt_bgcr%LitrfallElms_pvr          ,& !inoput :plant LitrFall element, [g d-2 h-1]
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft         ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    StandDeadKCompElms_pft      => plt_biom%StandDeadKCompElms_pft        ,& !inoput :standing dead element fraction, [g d-2]
    iDayPlanting_pft            => plt_distb%iDayPlanting_pft             ,& !output :day of planting,[-]
    iYearPlanting_pft           => plt_distb%iYearPlanting_pft             & !output :year of planting,[-]
  )
!     LitrFall AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM SHOOT AT DEATH
!
!     iPlantShootState_pft,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!     doInitPlant_pft=PFT initialization flag:0=no,1=yes
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     C4PhotoShootNonstC_brch=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
! whole plant is dead
  IF(iPlantShootState_pft(NZ).EQ.iDead .AND. iPlantRootState_pft(NZ).EQ.iDead)THEN
    !both plant shoots and roots are dead
    IF(doInitPlant_pft(NZ).EQ.ifalse)THEN      
      !surface litterfall
      D8825: DO NB=1,NumOfBranches_pft(NZ)
        D6425: DO M=1,jsken
          LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ) &
            +PlantElmAllocMat4Litr(ielmc,inonstruct,M,NZ)*AZMAX1(C4PhotoShootNonstC_brch(NB,NZ))

          DO NE=1,NumPlantChemElms
            LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafAlloc2Litr(NE,k_woody_comp) &
              +PetoleStrutElms_brch(NE,NB,NZ)*FracShootPetolAlloc2Litr(NE,k_woody_comp))

            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*(CanopyNonstElms_brch(NE,NB,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ) &
              +StalkRsrvElms_brch(NE,NB,NZ)) &
              +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafAlloc2Litr(NE,k_fine_comp) &
              +CanopyNodulStrutElms_brch(NE,NB,NZ)) &
              +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootPetolAlloc2Litr(NE,k_fine_comp) &
              +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))

            IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
              SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ) &
                *AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
            ELSE
              LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
                +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
            ENDIF

            IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. .not.is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
              !all above ground, or not a tree
              LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)&
                +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
            ELSE
              !planting mortality adds to standing dead
              dDeadE                          = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
              StandDeadKCompElms_pft(NE,M,NZ) = StandDeadKCompElms_pft(NE,M,NZ)+dDeadE
            ENDIF
          ENDDO
        ENDDO D6425  
      ENDDO D8825      
      !
      !     LitrFall AND STATE VARIABLES FOR SEASONAL STORAGE
      !     RESERVES FROM ROOT AND STORGE AT DEATH
      !        
      D6415: DO L=NU,MaxNumRootLays
        DO N=1,Myco_pft(NZ)
          DO M=1,jsken
            DO NE=1,NumPlantChemElms
              LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
                +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ))
            ENDDO    

            DO NR=1,NumPrimeRootAxes_pft(NZ)
              DO NE=1,NumPlantChemElms
                dRootMyco=AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))
                LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)&
                  +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dRootMyco*FracRootElmAlloc2Litr(NE,k_woody_comp)

                LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
                  +PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dRootMyco*FracRootElmAlloc2Litr(NE,k_fine_comp)
              ENDDO
            ENDDO  
          ENDDO  
        ENDDO

        DO M=1,jsken
          DO NR=1,NumPrimeRootAxes_pft(NZ)
            DO NE=1,NumPlantChemElms
              dRootMyco=AZMAX1(RootMyco1stStrutElms_rpvr(NE,L,NR,NZ))     
              LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)&
                +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dRootMyco*FracRootElmAlloc2Litr(NE,k_woody_comp)

              LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
                +PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dRootMyco*FracRootElmAlloc2Litr(NE,k_fine_comp)
            ENDDO    
          ENDDO      
        ENDDO
      ENDDO D6415

      DO M=1,jsken
        DO NE=1,NumPlantChemElms  
          LitrfallElms_pvr(NE,M,k_woody_comp,NGTopRootLayer_pft(NZ),NZ)=&
             LitrfallElms_pvr(NE,M,k_woody_comp,NGTopRootLayer_pft(NZ),NZ) &
            +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(SeasonalNonstElms_pft(NE,NZ))*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)

          LitrfallElms_pvr(NE,M,k_fine_comp,NGTopRootLayer_pft(NZ),NZ)= &
             LitrfallElms_pvr(NE,M,k_fine_comp,NGTopRootLayer_pft(NZ),NZ) &
            +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(SeasonalNonstElms_pft(NE,NZ))*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
        ENDDO
      ENDDO
!
      call ResetBranchRootStates(NZ,C4PhotoShootNonstC_brch)
    ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     jHarvstType_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     DazCurrYear=number of days in current year
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. jHarvstType_pft(NZ).EQ.jharvtyp_noaction)THEN
      IF(I.LT.DazCurrYear)THEN
        iDayPlanting_pft(NZ)  = I+1
        iYearPlanting_pft(NZ) = iYearCurrent
      ELSE
        iDayPlanting_pft(NZ)  = 1
        iYearPlanting_pft(NZ) = iYearCurrent+1
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine LiterFallRootShootStorage

!----------------------------------------------------------------------------------------------------
  subroutine LiterFallDeadRoots(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,M,NR,N,NE,NTG
!     begin_execution
  associate(                                                          &
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr     ,& !input  :C woody fraction in root,[-]
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft        ,& !input  :soil layer at planting depth, [-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft      ,& !input  :maximum radius of primary roots, [m]
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft      ,& !input  :maximum radius of secondary roots, [m]
    SeedDepth_pft             => plt_morph%SeedDepth_pft             ,& !input  :seeding depth, [m]
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    iPlantRootState_pft       => plt_pheno%iPlantRootState_pft       ,& !input  :flag to detect root system death,[-]
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    inonstruct                => pltpar%inonstruct                   ,& !input  :group id of plant nonstructural litter
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    NActiveRootSegs_raxes     => plt_morph%NActiveRootSegs_raxes     ,& !number of active root segments    
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    NumPrimeRootAxes_pft      => plt_morph%NumPrimeRootAxes_pft      ,& !inoput :root primary axis number,[-]
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft     ,& !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !inoput :root layer element primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr    ,& !inoput :root layer nonstructural element, [g d-2]
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr    ,& !inoput :root layer nodule element, [g d-2]
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr            ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr            ,& !inoput :root aqueous content, [g d-2]
    RootSegBaseDepth_raxes    => plt_morph%RootSegBaseDepth_raxes    ,& !base depth of different root axes, [m]    
    NMaxRootBotLayer_pft      => plt_morph%NMaxRootBotLayer_pft      ,& !output :maximum soil layer number for all root axes, [-]
    NRoot1stTipLay_raxes      => plt_morph%NRoot1stTipLay_raxes      ,& !output :maximum soil layer number for root axes, [-]
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr         ,& !output :root layer C, [gC d-2]
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft           ,& !output :root layer depth, [m]
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr           ,& !output :root layer length primary axes, [m d-2]
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr         ,& !output :root layer diameter primary axes, [m]
    Root1stXNumL_rpvr         => plt_morph%Root1stXNumL_rpvr         ,& !output :root layer number primary axes, [d-2]
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr           ,& !output :root layer length secondary axes, [m d-2]
    Root2ndEffLen4uptk_rpvr   => plt_morph%Root2ndEffLen4uptk_rpvr   ,& !output :root layer average length, [m]
    Root2ndRadius_rpvr        => plt_morph%Root2ndRadius_rpvr        ,& !output :root layer diameter secondary axes, [m]
    Root2ndXNumL_rpvr         => plt_morph%Root2ndXNumL_rpvr         ,& !output :root layer number axes, [d-2]
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr          ,& !output :root layer number secondary axes, [d-2]
    RootAreaPerPlant_pvr      => plt_morph%RootAreaPerPlant_pvr      ,& !output :root layer area per plant, [m p-1]
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr   ,& !output :root layer length density, [m m-3]
    RootTotLenPerPlant_pvr    => plt_morph%RootTotLenPerPlant_pvr    ,& !output :root layer length per plant, [m p-1]
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs        ,& !output :root C primary axes, [g d-2]
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr    ,& !output :root layer structural C, [gC d-2]
    RootPoreVol_rpvr          => plt_morph%RootPoreVol_rpvr          ,& !output :root layer volume air, [m2 d-2]
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr           ,& !output :root layer protein C, [gC d-2]
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr               & !output :root layer volume water, [m2 d-2]
  )
!     IDTHR=PFT root living flag: 0=alive,1=dead
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!

  IF(iPlantRootState_pft(NZ).EQ.iDead)THEN
    !add root to LitrFall
    DO L=NU,MaxNumRootLays        
      DO N=1,Myco_pft(NZ)
        DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
              +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)* RootMycoNonstElms_rpvr(NE,N,L,NZ)
          ENDDO 
        ENDDO
      ENDDO  
    ENDDO  


    DO  NR=1,NumPrimeRootAxes_pft(NZ)
      DO L=NU,MaxNumRootLays        

        DO N=1,Myco_pft(NZ)
          DO M=1,jsken
            DO NE=1,NumPlantChemElms
              LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)&
                +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)&
                *FracRootElmAlloc2Litr(NE,k_woody_comp)

              LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
                +PlantElmAllocMat4Litr(NE,iroot,M,NZ)*RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)&
                *FracRootElmAlloc2Litr(NE,k_fine_comp)
            enddo
          ENDDO
        ENDDO 

        DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)&
              +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)&
              *FracRootElmAlloc2Litr(NE,k_woody_comp)

            LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
              +PlantElmAllocMat4Litr(NE,iroot,M,NZ)*RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)&
              *FracRootElmAlloc2Litr(NE,k_fine_comp)
          enddo
        ENDDO

      ENDDO
    ENDDO    


    DO L=NU,MaxNumRootLays             
      DO N=1,Myco_pft(NZ)
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        DO NTG=idg_beg,idg_NH3
          RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-trcg_rootml_pvr(NTG,N,L,NZ)&
            -trcs_rootml_pvr(NTG,N,L,NZ)
        ENDDO
        trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
        trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
      ENDDO
    ENDDO    
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!
    D8870: DO NR=1,NumPrimeRootAxes_pft(NZ)
      DO L=NU,MaxNumRootLays       
        RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = 0._r8    
        Root1stLen_rpvr(L,NR,NZ)                              = 0._r8                
        DO N=1,Myco_pft(NZ)        
          RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
          Root2ndLen_rpvr(N,L,NR,NZ)                              = 0._r8
          Root2ndXNum_rpvr(N,L,NR,NZ)                             = 0._r8
        ENDDO
      ENDDO          
      RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ)=0._r8            
    ENDDO D8870

    DO L=NU,MaxNumRootLays             
      DO N=1,Myco_pft(NZ)           
        RootMycoNonstElms_rpvr(:,N,L,NZ) = 0._r8
        RootMycoActiveBiomC_pvr(N,L,NZ)  = 0._r8
        PopuRootMycoC_pvr(N,L,NZ)        = 0._r8
        RootProteinC_pvr(N,L,NZ)         = 0._r8
        Root1stXNumL_rpvr(N,L,NZ)         = 0._r8
        Root2ndXNumL_rpvr(N,L,NZ)          = 0._r8
        RootTotLenPerPlant_pvr(N,L,NZ)      = 0._r8
        RootLenDensPerPlant_pvr(N,L,NZ)  = 0._r8
        RootPoreVol_rpvr(N,L,NZ)          = 0._r8
        RootVH2O_pvr(N,L,NZ)             = 0._r8
        Root1stRadius_pvr(N,L,NZ)        = Root1stMaxRadius_pft(N,NZ)
        Root2ndRadius_rpvr(N,L,NZ)        = Root2ndMaxRadius_pft(N,NZ)
        RootAreaPerPlant_pvr(N,L,NZ)     = 0._r8
        Root2ndEffLen4uptk_rpvr(N,L,NZ)        = Root2ndTipLen4uptk
      ENDDO
    ENDDO    
!
!     LitrFall AND STATE VARIABLES FROM DEAD NODULES
!
!     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
    DO N=1,Myco_pft(NZ)       
      IF(is_plant_N2fix(iPlantNfixType_pft(NZ)).AND.N.EQ.ipltroot)THEN
        DO L=NU,MaxNumRootLays
          D6420: DO M=1,jsken
            DO NE=1,NumPlantChemElms            
              LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+&
                PlantElmAllocMat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_rpvr(NE,L,NZ)+&
                PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_rpvr(NE,L,NZ)
            ENDDO
          ENDDO D6420
          RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ)=0._r8
          RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ)=0._r8  
        ENDDO          
      ENDIF
    ENDDO 
!
!     RESET DEPTH VARIABLES OF DEAD ROOTS
!
!   
    D8795: DO NR=1,NumPrimeRootAxes_pft(NZ)
     NRoot1stTipLay_raxes(NR,NZ)   = NGTopRootLayer_pft(NZ)
     Root1stDepz_pft(NR,NZ)        = SeedDepth_pft(NZ)
     RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)
     NActiveRootSegs_raxes(NR,NZ)  = 0._r8
    ENDDO D8795
    NMaxRootBotLayer_pft(NZ) = NGTopRootLayer_pft(NZ)
    NumPrimeRootAxes_pft(NZ)     = 0
  ENDIF
  end associate
  end subroutine LiterFallDeadRoots

!----------------------------------------------------------------------------------------------------
  subroutine LiterFallDeadBranches(I,J,NZ,NumDeadBranches,C4PhotoShootNonstC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: NumDeadBranches
  real(r8), intent(inout) :: C4PhotoShootNonstC_brch(NumCanopyLayers1,JP1)
  real(r8) :: dDeadE
  integer :: M,NE,NB
!     begin_execution
  associate(                                                                           &
    CanopyNodulNonstElms_brch         => plt_biom%CanopyNodulNonstElms_brch           ,& !input  :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch         => plt_biom%CanopyNodulStrutElms_brch           ,& !input  :branch nodule structural element, [g d-2]
    CanopyNonstElms_brch              => plt_biom%CanopyNonstElms_brch                ,& !input  :branch nonstructural element, [g d-2]
    EarStrutElms_brch                 => plt_biom%EarStrutElms_brch                   ,& !input  :branch ear structural chemical element mass, [g d-2]
    PlantElmAllocMat4Litr             => plt_soilchem%PlantElmAllocMat4Litr           ,& !input  :litter kinetic fraction, [-]
    FracShootLeafAlloc2Litr           => plt_allom%FracShootLeafAlloc2Litr            ,& !input  :woody element allocation, [-]
    FracShootPetolAlloc2Litr          => plt_allom%FracShootPetolAlloc2Litr           ,& !input  :leaf element allocation,[-]
    GrainStrutElms_brch               => plt_biom%GrainStrutElms_brch                 ,& !input  :branch grain structural element mass, [g d-2]
    HuskStrutElms_brch                => plt_biom%HuskStrutElms_brch                  ,& !input  :branch husk structural element mass, [g d-2]
    LeafStrutElms_brch                => plt_biom%LeafStrutElms_brch                  ,& !input  :branch leaf structural element mass, [g d-2]
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                    ,& !input  :acclimated plant maturity group, [-]
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft                  ,& !input  :number of branches,[-]
    PetoleStrutElms_brch              => plt_biom%PetoleStrutElms_brch                ,& !input  :branch sheath structural element, [g d-2]
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft         ,& !input  :number of nodes in seed, [-]
    StalkRsrvElms_brch                => plt_biom%StalkRsrvElms_brch                  ,& !input  :branch reserve element mass, [g d-2]
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch                 ,& !input  :branch stalk structural element mass, [g d-2]
    iHarvstType_pft                   => plt_distb%iHarvstType_pft                    ,& !input  :type of harvest,[-]
    iPlantBranchState_brch            => plt_pheno%iPlantBranchState_brch             ,& !input  :flag to detect branch death, [-]
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft            ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft               ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft              ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft          ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    icwood                            => pltpar%icwood                                ,& !input  :group id of coarse woody litter
    ifoliar                           => pltpar%ifoliar                               ,& !input  :group id of plant foliar litter
    inonfoliar                        => pltpar%inonfoliar                            ,& !input  :group id of plant non-foliar litter group
    inonstruct                        => pltpar%inonstruct                            ,& !input  :group id of plant nonstructural litter
    istalk                            => pltpar%istalk                                ,& !input  :group id of plant stalk litter group
    k_fine_comp                       => pltpar%k_fine_comp                           ,& !input  :fine litter complex id
    k_woody_comp                      => pltpar%k_woody_comp                          ,& !input  :woody litter complex id
    LitrfallElms_pvr              => plt_bgcr%LitrfallElms_pvr                ,& !inoput :plant LitrFall element, [g d-2 h-1]
    SeasonalNonstElms_pft             => plt_biom%SeasonalNonstElms_pft               ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    StandDeadKCompElms_pft            => plt_biom%StandDeadKCompElms_pft              ,& !inoput :standing dead element fraction, [g d-2]
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                  ,& !output :shoot node number, [-]
    BranchNumerID_brch                => plt_morph%BranchNumerID_brch                 ,& !output :branch meric id, [-]
    C4PhotosynDowreg_brch             => plt_photo%C4PhotosynDowreg_brch              ,& !output :down-regulation of C4 photosynthesis, [-]
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch             ,& !output :flag to detect physiological maturity from grain fill, [-]
    Hours2LeafOut_brch                => plt_pheno%Hours2LeafOut_brch                 ,& !output :counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 ,& !output :cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch                 ,& !output :heat requirement for spring leafout/dehardening, [h]
    Hours4LenthenPhotoPeriod_brch     => plt_pheno%Hours4LenthenPhotoPeriod_brch      ,& !output :initial heat requirement for spring leafout/dehardening, [h]
    Hours4LiterfalAftMature_brch      => plt_pheno%Hours4LiterfalAftMature_brch       ,& !output :branch phenology flag, [h]
    Hours4ShortenPhotoPeriod_brch     => plt_pheno%Hours4ShortenPhotoPeriod_brch      ,& !output :initial cold requirement for autumn leafoff/hardening, [h]
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch             ,& !output :leaf growth stage counter, [-]
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch                   ,& !output :leaf number, [-]
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch        ,& !output :leaf number at floral initiation, [-]
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch                   ,& !output :plant maturity group, [-]
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch            ,& !output :shoot node number at floral initiation, [-]
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch          ,& !output :shoot node number at anthesis, [-]
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch                   ,& !output :leaf number, [-]
    Prep4Literfall_brch               => plt_pheno%Prep4Literfall_brch                ,& !output :branch phenology flag, [-]
    RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch               ,& !output :branch down-regulation of CO2 fixation, [-]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !output :normalized node number during reproductive growth stages, [-]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !output :normalized node number during vegetative growth stages, [-]
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch                 ,& !output :branch phenology flag, [-]
    doPlantLeafOut_brch               => plt_pheno%doPlantLeafOut_brch                ,& !output :branch phenology flag, [-]
    doPlantLeaveOff_brch              => plt_pheno%doPlantLeaveOff_brch               ,& !output :branch phenology flag, [-]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                 & !output :plant growth stage, [-]
  )
  D8845: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iDead)THEN
      MatureGroup_brch(NB,NZ)                  = MatureGroup_pft(NZ)
      ShootNodeNum_brch(NB,NZ)                 = ShootNodeNumAtPlanting_pft(NZ)
      NodeNum2InitFloral_brch(NB,NZ)           = ShootNodeNum_brch(NB,NZ)
      NodeNumberAtAnthesis_brch(NB,NZ)         = 0.0_r8
      NumOfLeaves_brch(NB,NZ)                  = 0.0_r8
      LeafNumberAtFloralInit_brch(NB,NZ)       = 0.0_r8
      KLeafNumber_brch(NB,NZ)                  = 1
      KHiestGroLeafNode_brch(NB,NZ)            = 1
      TotalNodeNumNormByMatgrp_brch(NB,NZ)     = 0.0_r8
      TotReproNodeNumNormByMatrgrp_brch(NB,NZ) = 0.0_r8
      Hours4Leafout_brch(NB,NZ)                = 0.0_r8
      Hours4LeafOff_brch(NB,NZ)                = 0.0_r8
      Hours4LenthenPhotoPeriod_brch(NB,NZ)     = 0.0_r8
      Hours4ShortenPhotoPeriod_brch(NB,NZ)     = 0.0_r8
      Hours2LeafOut_brch(NB,NZ)                = 0.0_r8
      HourFailGrainFill_brch(NB,NZ)            = 0.0_r8
      RubiscoActivity_brch(NB,NZ)              = 1.0_r8
      C4PhotosynDowreg_brch(NB,NZ)             = 1.0_r8
      doInitLeafOut_brch(NB,NZ)                = iEnable
      doPlantLeafOut_brch(NB,NZ)               = iDisable
      doPlantLeaveOff_brch(NB,NZ)              = iEnable
      Prep4Literfall_brch(NB,NZ)               = ifalse
      Hours4LiterfalAftMature_brch(NB,NZ)      = 0
      BranchNumerID_brch(NB,NZ)                = 0
      D8850: DO M=1,pltpar%NumGrowthStages
        iPlantCalendar_brch(M,NB,NZ) = 0
      ENDDO D8850
      !
      !     LitrFall FROM DEAD BRANCHES
      !
      !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
      !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
      !     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
      !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
      !
      D6405: DO M=1,jsken
        DO NE=1,NumPlantChemElms        
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
            +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*CanopyNodulNonstElms_brch(NE,NB,NZ) &
            +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafAlloc2Litr(NE,k_fine_comp)+CanopyNodulStrutElms_brch(NE,NB,NZ)) &
            +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootPetolAlloc2Litr(NE,k_fine_comp)+HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))

          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) &
            +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafAlloc2Litr(NE,k_woody_comp)+PetoleStrutElms_brch(NE,NB,NZ)*FracShootPetolAlloc2Litr(NE,k_woody_comp))

          IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
            SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ) &
              +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
          ELSE
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
          ENDIF
          !terminate all-aboveground or plant is not tree
          IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. .not.is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)
          ELSE
            !standing dead
            dDeadE                          = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)
            StandDeadKCompElms_pft(NE,M,NZ) = StandDeadKCompElms_pft(NE,M,NZ)+dDeadE              
          ENDIF
        ENDDO
      ENDDO D6405

      !
      !     RECOVER NON-STRUCTURAL C,N,P FROM BRANCH TO
      !     SEASONAL STORAGE RESERVES
      !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
      !                       ,3=pruning,4=grazing,5=fire,6=herbivory
      !
      SeasonalNonstElms_pft(ielmc,NZ)=SeasonalNonstElms_pft(ielmc,NZ)+C4PhotoShootNonstC_brch(NB,NZ)

      DO NE=1,NumPlantChemElms
        SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+CanopyNonstElms_brch(NE,NB,NZ)
      ENDDO  

      IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
        D6406: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*StalkRsrvElms_brch(NE,NB,NZ)
          ENDDO    
        ENDDO D6406
      ELSE
        DO NE=1,NumPlantChemElms
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+StalkRsrvElms_brch(NE,NB,NZ)
        ENDDO  
      ENDIF

      call ResetDeadRootStates(NB,NZ,C4PhotoShootNonstC_brch)

      NumDeadBranches=NumDeadBranches+1
    ENDIF
  ENDDO D8845
  end associate
  end subroutine LiterFallDeadBranches

!----------------------------------------------------------------------------------------------------
  subroutine ResetBranchRootStates(NZ,C4PhotoShootNonstC_brch)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: C4PhotoShootNonstC_brch(NumCanopyLayers1,JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                                                          &
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    NumPrimeRootAxes_pft      => plt_morph%NumPrimeRootAxes_pft      ,& !input  :root primary axis number,[-]
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch  ,& !output :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch  ,& !output :branch nodule structural element, [g d-2]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch       ,& !output :branch nonstructural element, [g d-2]
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch          ,& !output :branch ear structural chemical element mass, [g d-2]
    GrainStrutElms_brch       => plt_biom%GrainStrutElms_brch        ,& !output :branch grain structural element mass, [g d-2]
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch         ,& !output :branch husk structural element mass, [g d-2]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch     ,& !output :plant branch leaf + sheath C, [g d-2]
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch         ,& !output :branch leaf structural element mass, [g d-2]
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch       ,& !output :branch sheath structural element, [g d-2]
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr           ,& !output :root layer length primary axes, [m d-2]
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr           ,& !output :root layer length secondary axes, [m d-2]
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr          ,& !output :root layer number secondary axes, [d-2]
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs        ,& !output :root C primary axes, [g d-2]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !output :root layer element primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !output :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !output :root layer nonstructural element, [g d-2]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft      ,& !output :plant stored nonstructural element at current step, [g d-2]
    SenecStalkStrutElms_brch  => plt_biom%SenecStalkStrutElms_brch   ,& !output :branch stalk structural element, [g d-2]
    ShootElms_brch            => plt_biom%ShootElms_brch             ,& !output :branch shoot structural element mass, [g d-2]
    SapwoodBiomassC_brch      => plt_biom%SapwoodBiomassC_brch       ,& !output :branch live stalk C, [gC d-2]
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch         ,& !output :branch reserve element mass, [g d-2]
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch        ,& !output :branch stalk structural element mass, [g d-2]
    iPlantState_pft           => plt_pheno%iPlantState_pft            & !output :flag for species death, [-]
  )
!     RESET BRANCH STATE VARIABLES
!
  
  DO NB=1,NumOfBranches_pft(NZ)
    DO NE=1,NumPlantChemElms
      CanopyNonstElms_brch(NE,NB,NZ)      = 0._r8
      CanopyNodulNonstElms_brch(NE,NB,NZ) = 0._r8
      ShootElms_brch(NE,NB,NZ)            = 0._r8
      LeafStrutElms_brch(NE,NB,NZ)        = 0._r8
      PetoleStrutElms_brch(NE,NB,NZ)      = 0._r8
      StalkStrutElms_brch(NE,NB,NZ)       = 0._r8
      StalkRsrvElms_brch(NE,NB,NZ)        = 0._r8
    ENDDO
  ENDDO

  D8835: DO NB=1,NumOfBranches_pft(NZ)
    C4PhotoShootNonstC_brch(NB,NZ)                      = 0._r8
    SapwoodBiomassC_brch(NB,NZ)                         = 0._r8
    CanopyNodulStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
    HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
    EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)         = 0._r8
    GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
    CanopyLeafSheathC_brch(NB,NZ)                       = 0._r8
    SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)  = 0._r8
  ENDDO D8835
!
!     RESET ROOT STATE VARIABLES
!
  DO  NR=1,NumPrimeRootAxes_pft(NZ)    
    RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ)=0._r8    
  ENDDO

  D6416: DO L=NU,MaxNumRootLays
    DO  NR=1,NumPrimeRootAxes_pft(NZ)
      RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = 0._r8
      Root1stLen_rpvr(L,NR,NZ)                              = 0._r8      
    ENDDO  
    DO  N=1,Myco_pft(NZ)
       RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ)=0._r8
      DO  NR=1,NumPrimeRootAxes_pft(NZ)
        RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
        Root2ndLen_rpvr(N,L,NR,NZ)                              = 0._r8
        Root2ndXNum_rpvr(N,L,NR,NZ)                             = 0._r8
      enddo
    enddo
  ENDDO D6416
  SeasonalNonstElms_pft(1:NumPlantChemElms,NZ) = 0._r8
  iPlantState_pft(NZ)                          = iDead
  end associate
  end subroutine ResetBranchRootStates

!----------------------------------------------------------------------------------------------------
  subroutine ResetDeadRootStates(NB,NZ,C4PhotoShootNonstC_brch)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: C4PhotoShootNonstC_brch(NumCanopyLayers1,JP1)
  integer :: L,K,N
!     begin_execution
  associate(                                                             &
    CanopyLeafAreaZ_pft        => plt_morph%CanopyLeafAreaZ_pft         ,& !inoput :canopy layer leaf area, [m2 d-2]
    CanopyLeafArea_lnode       => plt_morph%CanopyLeafArea_lnode        ,& !inoput :layer/node/branch leaf area, [m2 d-2]
    CanopyLeafCLyr_pft         => plt_biom%CanopyLeafCLyr_pft           ,& !inoput :canopy layer leaf C, [g d-2]
    LeafLayerElms_node   => plt_biom%LeafLayerElms_node     ,& !inoput :layer leaf element, [g d-2]
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node   ,& !output :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node  ,& !output :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CPOOL3_node                => plt_photo%CPOOL3_node                 ,& !output :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                => plt_photo%CPOOL4_node                 ,& !output :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    CanopyNodulNonstElms_brch  => plt_biom%CanopyNodulNonstElms_brch    ,& !output :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch  => plt_biom%CanopyNodulStrutElms_brch    ,& !output :branch nodule structural element, [g d-2]
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch         ,& !output :branch nonstructural element, [g d-2]
    CanopyStalkArea_lbrch      => plt_morph%CanopyStalkArea_lbrch       ,& !output :plant canopy layer branch stem area, [m2 d-2]
    EarStrutElms_brch          => plt_biom%EarStrutElms_brch            ,& !output :branch ear structural chemical element mass, [g d-2]
    GrainSeedBiomCMean_brch    => plt_allom%GrainSeedBiomCMean_brch     ,& !output :maximum grain C during grain fill, [g d-2]
    GrainStrutElms_brch        => plt_biom%GrainStrutElms_brch          ,& !output :branch grain structural element mass, [g d-2]
    HuskStrutElms_brch         => plt_biom%HuskStrutElms_brch           ,& !output :branch husk structural element mass, [g d-2]
    StalkNodeVertLength_brch   => plt_morph%StalkNodeVertLength_brch    ,& !output :internode height, [m]
    StructInternodeElms_brch   => plt_biom%StructInternodeElms_brch     ,& !output :internode C, [g d-2]
    LeafAreaLive_brch          => plt_morph%LeafAreaLive_brch           ,& !output :branch leaf area, [m2 d-2]
    LeafAreaZsec_brch          => plt_morph%LeafAreaZsec_brch           ,& !output :leaf surface area, [m2 d-2]
    LeafElmntNode_brch         => plt_biom%LeafElmntNode_brch           ,& !output :leaf element, [g d-2]
    LeafArea_node              => plt_morph%LeafArea_node               ,& !output :leaf area, [m2 d-2]
    CanopyLeafSheathC_brch     => plt_biom%CanopyLeafSheathC_brch       ,& !output :plant branch leaf + sheath C, [g d-2]
    LeafProteinC_node          => plt_biom%LeafProteinC_node            ,& !output :layer leaf protein C, [g d-2]
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch           ,& !output :branch leaf structural element mass, [g d-2]
    StalkNodeHeight_brch       => plt_morph%StalkNodeHeight_brch        ,& !output :internode height, [m]
    PetioleElmntNode_brch      => plt_biom%PetioleElmntNode_brch        ,& !output :sheath chemical element, [g d-2]
    PetoleLength_node        => plt_morph%PetoleLength_node         ,& !output :sheath height, [m]
    PetoleProteinC_node    => plt_biom%PetoleProteinC_node      ,& !output :layer sheath protein C, [g d-2]
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch         ,& !output :branch sheath structural element, [g d-2]
    PotentialSeedSites_brch    => plt_morph%PotentialSeedSites_brch     ,& !output :branch potential grain number, [d-2]
    SeedSitesSet_brch          => plt_morph%SeedSitesSet_brch           ,& !output :branch grain number, [d-2]
    SenecStalkStrutElms_brch   => plt_biom%SenecStalkStrutElms_brch     ,& !output :branch stalk structural element, [g d-2]
    ShootElms_brch             => plt_biom%ShootElms_brch               ,& !output :branch shoot structural element mass, [g d-2]
    SapwoodBiomassC_brch       => plt_biom%SapwoodBiomassC_brch         ,& !output :branch live stalk C, [gC d-2]
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch           ,& !output :branch reserve element mass, [g d-2]
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch          ,& !output :branch stalk structural element mass, [g d-2]
    StemAreaZsec_brch          => plt_morph%StemAreaZsec_brch            & !output :stem surface area, [m2 d-2]
  )
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     SeedSitesSet_brch=seed set number
!     PotentialSeedSites_brch=potential number of seed set sites
!     GrainSeedBiomCMean_brch=individual seed size
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     LeafProteinC_node=leaf protein mass
!     LeafAreaLive_brch=branch leaf area
!     WGLF,WGLFN,WGLFP,LeafProteinC_node=node leaf C,N,P,protein mass
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinC_node=node petiole C,N,P,protein mass
!     StructInternodeElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
  C4PhotoShootNonstC_brch(NB,NZ)                      = 0._r8
  CanopyNonstElms_brch(1:NumPlantChemElms,NB,NZ)      = 0._r8
  CanopyNodulNonstElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
  ShootElms_brch(1:NumPlantChemElms,NB,NZ)            = 0._r8
  LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
  CanopyNodulStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
  PetoleStrutElms_brch(1:NumPlantChemElms,NB,NZ)      = 0._r8
  StalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
  StalkRsrvElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
  HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
  EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)         = 0._r8
  GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
  SapwoodBiomassC_brch(NB,NZ)                         = 0._r8
  CanopyLeafSheathC_brch(NB,NZ)                       = 0._r8
  PotentialSeedSites_brch(NB,NZ)                      = 0._r8
  SeedSitesSet_brch(NB,NZ)                            = 0._r8
  GrainSeedBiomCMean_brch(NB,NZ)                      = 0._r8
  LeafAreaLive_brch(NB,NZ)                            = 0._r8
  SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)  = 0._r8

  D8855: DO K=0,MaxNodesPerBranch1
    IF(K.NE.0)THEN
      CPOOL3_node(K,NB,NZ)                = 0._r8
      CPOOL4_node(K,NB,NZ)                = 0._r8
      CMassCO2BundleSheath_node(K,NB,NZ)  = 0._r8
      CMassHCO3BundleSheath_node(K,NB,NZ) = 0._r8
    ENDIF
    LeafArea_node(K,NB,NZ)                          = 0._r8
    StalkNodeHeight_brch(K,NB,NZ)                    = 0._r8
    StalkNodeVertLength_brch(K,NB,NZ)                  = 0._r8
    PetoleLength_node(K,NB,NZ)                        = 0._r8
    LeafProteinC_node(K,NB,NZ)                      = 0._r8
    PetoleProteinC_node(K,NB,NZ)                    = 0._r8
    LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)      = 0._r8
    PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)   = 0._r8
    StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
    D8865: DO L=1,NumCanopyLayers1
      CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)-CanopyLeafArea_lnode(L,K,NB,NZ)
      CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)-LeafLayerElms_node(ielmc,L,K,NB,NZ)
      CanopyLeafArea_lnode(L,K,NB,NZ)                         = 0._r8
      LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ) = 0._r8
      IF(K.NE.0)THEN
        D8860: DO N=1,NumLeafZenithSectors1
          LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
        ENDDO D8860
      ENDIF
    ENDDO D8865
  ENDDO D8855
  D8875: DO L=1,NumCanopyLayers1
    CanopyStalkArea_lbrch(L,NB,NZ)=0._r8
    DO  N=1,NumLeafZenithSectors1
      StemAreaZsec_brch(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D8875
  end associate
  end subroutine ResetDeadRootStates
  ![tail]
end module LitterFallMod
