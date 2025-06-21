module LitterFallMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  use PlantMathFuncMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: ResetDeadPlant
  public :: ReSeedPlants
  contains
!------------------------------------------------------------------------------------------

  subroutine ResetDeadPlant(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NumDeadBranches

!     begin_execution
  associate(                                                       &
    C4PhotoShootNonstC_brch => plt_biom%C4PhotoShootNonstC_brch  , & !input  
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft       , & !input  
    SolarNoonHour_col       => plt_site%SolarNoonHour_col        , & !input  
    iDayPlantHarvest_pft    => plt_distb%iDayPlantHarvest_pft    , & !input  
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch     , & !input  
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft , & !input  
    iYearCurrent            => plt_site%iYearCurrent             , & !input  
    iYearPlantHarvest_pft   => plt_distb%iYearPlantHarvest_pft     & !input  
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
  IF(J.EQ.INT(SolarNoonHour_col) .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND. (iPlantPhenolPattern_pft(NZ).NE.iplt_annual .OR. (I.GE.iDayPlantHarvest_pft(NZ) &
    .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ))))THEN
    NumDeadBranches=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!

    call LiterfallFromDeadBranches(I,J,NZ,NumDeadBranches,C4PhotoShootNonstC_brch)

    call SetDeadPlant(i,j,NZ,NumDeadBranches)
!
!     DEAD ROOTS
!
!     LitrFall FROM DEAD ROOTS
!
    call LiterfallFromDeadRoots(I,J,NZ)

    call LiterfallFromRootShootStorage(I,J,NZ,C4PhotoShootNonstC_brch)

  ENDIF
  end associate
  end subroutine ResetDeadPlant
!------------------------------------------------------------------------------------------
  subroutine ReSeedPlants(I,J,NZ)
  !
  !Description:
  !Re-Seed plants if there is non-zero storage (sort of seed bank)
  !
  use InitPlantMod, only: InitPlantPhenoMorphoBio
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ

  associate(                                                  &
    SeasonalNonstElms_pft => plt_biom%SeasonalNonstElms_pft , & !input  
    ZERO4Groth_pft        => plt_biom%ZERO4Groth_pft        , & !input  
    doReSeed_pft          => plt_pheno%doReSeed_pft         , & !input  
    IsPlantActive_pft     => plt_pheno%IsPlantActive_pft    , & !output 
    SeedCPlanted_pft      => plt_biom%SeedCPlanted_pft        & !output 
  )
  if(doReSeed_pft(NZ) .and. SeasonalNonstElms_pft(ielmc,NZ) > ZERO4Groth_pft(NZ) )then

    IsPlantActive_pft(NZ) = iDormant
    SeedCPlanted_pft(NZ)  = SeasonalNonstElms_pft(ielmc,NZ)

    call InitPlantPhenoMorphoBio(NZ)

  endif  
  end associate
  end subroutine ReSeedPlants
!------------------------------------------------------------------------------------------
  subroutine SetDeadPlant(I,J,NZ,NumDeadBranches)
  implicit none
  INTEGER, intent(in) :: i,J
  integer, intent(in) :: NZ
  integer, intent(in) :: NumDeadBranches

  associate(                                                       &
    PlantPopulation_pft     => plt_site%PlantPopulation_pft      , & !input  
    RootElms_pft            => plt_biom%RootElms_pft             , & !input  
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft    , & !input  
    doInitPlant_pft         => plt_pheno%doInitPlant_pft         , & !input  
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft , & !input  
    jHarvst_pft             => plt_distb%jHarvst_pft             , & !input  
    CanopyBiomWater_pft     => plt_ew%CanopyBiomWater_pft        , & !inoput 
    H2OLoss_CumYr_col       => plt_ew%H2OLoss_CumYr_col          , & !inoput 
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft       , & !inoput 
    QH2OLoss_lnds           => plt_site%QH2OLoss_lnds            , & !inoput 
    iPlantRootState_pft     => plt_pheno%iPlantRootState_pft     , & !inoput 
    BranchNumber_pft        => plt_morph%BranchNumber_pft        , & !output 
    HoursTooLowPsiCan_pft   => plt_pheno%HoursTooLowPsiCan_pft   , & !output 
    HypoctoHeight_pft       => plt_morph%HypoctoHeight_pft       , & !output 
    doReSeed_pft            => plt_pheno%doReSeed_pft            , & !output 
    iPlantShootState_pft    => plt_pheno%iPlantShootState_pft      & !output 
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
    HypoctoHeight_pft(NZ) = 0._r8
    QH2OLoss_lnds         = QH2OLoss_lnds+CanopyBiomWater_pft(NZ)
    H2OLoss_CumYr_col     = H2OLoss_CumYr_col+CanopyBiomWater_pft(NZ)
    CanopyBiomWater_pft(NZ)   = 0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
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
    IF(jHarvst_pft(NZ).NE.jharvtyp_noaction)then
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
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromRootShootStorage(I,J,NZ,C4PhotoShootNonstC_brch)

  use EcoSIMCtrlDataType, only : iYearCurrent
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: C4PhotoShootNonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,M,NR,NB,N,NE
  REAL(R8) :: dRootMyco
!     begin_execution
  associate(                                                               &
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch    , & !input  
    CanopyNodulStrutElms_brch   => plt_biom%CanopyNodulStrutElms_brch    , & !input  
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch         , & !input  
    DazCurrYear                 => plt_site%DazCurrYear                  , & !input  
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch            , & !input  
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr         , & !input  
    FracRootElmAlloc2Litr       => plt_allom%FracRootElmAlloc2Litr       , & !input  
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr  , & !input  
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr  , & !input  
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr , & !input  
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch          , & !input  
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch           , & !input  
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch           , & !input  
    Myco_pft                      => plt_morph%Myco_pft                      , & !input  
    MaxNumRootLays              => plt_site%MaxNumRootLays               , & !input  
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft          , & !input  
    NU                          => plt_site%NU                           , & !input  
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft           , & !input  
    NumRootAxes_pft             => plt_morph%NumRootAxes_pft             , & !input  
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch         , & !input  
    RootMyco1stStrutElms_rpvr   => plt_biom%RootMyco1stStrutElms_rpvr    , & !input  
    RootMyco2ndStrutElms_rpvr   => plt_biom%RootMyco2ndStrutElms_rpvr    , & !input  
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr       , & !input  
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch           , & !input  
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch          , & !input  
    doInitPlant_pft             => plt_pheno%doInitPlant_pft             , & !input  
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft     , & !input  
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft        , & !input  
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft       , & !input  
    iPlantRootState_pft         => plt_pheno%iPlantRootState_pft         , & !input  
    iPlantShootState_pft        => plt_pheno%iPlantShootState_pft        , & !input  
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft   , & !input  
    icwood                      => pltpar%icwood                         , & !input  
    ifoliar                     => pltpar%ifoliar                        , & !input  
    inonfoliar                  => pltpar%inonfoliar                     , & !input  
    inonstruct                  => pltpar%inonstruct                     , & !input  
    iroot                       => pltpar%iroot                          , & !input  
    istalk                      => pltpar%istalk                         , & !input  
    jHarvst_pft                 => plt_distb%jHarvst_pft                 , & !input  
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr         , & !inoput 
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft        , & !inoput 
    StandDeadKCompElms_pft      => plt_biom%StandDeadKCompElms_pft       , & !inoput 
    k_fine_litr                 => pltpar%k_fine_litr                    , & !inoput 
    k_woody_litr                => pltpar%k_woody_litr                   , & !inoput 
    iDayPlanting_pft            => plt_distb%iDayPlanting_pft            , & !output 
    iYearPlanting_pft           => plt_distb%iYearPlanting_pft             & !output 
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
!
  IF(iPlantShootState_pft(NZ).EQ.iDead .AND. iPlantRootState_pft(NZ).EQ.iDead)THEN
    !both plant shoots and roots are dead
    IF(doInitPlant_pft(NZ).EQ.ifalse)THEN      
      D8825: DO NB=1,NumOfBranches_pft(NZ)
        D6425: DO M=1,jsken
          LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmc,inonstruct,M,NZ)*AZMAX1(C4PhotoShootNonstC_brch(NB,NZ))

          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr) &
              +PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr))

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*(CanopyNonstElms_brch(NE,NB,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ) &
              +StalkRsrvElms_brch(NE,NB,NZ)) &
              +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr) &
              +CanopyNodulStrutElms_brch(NE,NB,NZ)) &
              +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr) &
              +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))

            IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
              SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+ElmAllocmat4Litr(NE,inonfoliar,M,NZ) &
                *AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
            ELSE
              LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
                +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
            ENDIF

            IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. .not.is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
              !all above ground, or not a tree
              LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)&
                +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
            ELSE
              StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ)&
                +ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
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
              LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
                +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ))
            ENDDO    

            DO NR=1,NumRootAxes_pft(NZ)
              DO NE=1,NumPlantChemElms
                dRootMyco=AZMAX1(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))
                LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)&
                  +ElmAllocmat4Litr(NE,icwood,M,NZ)*dRootMyco*FracRootElmAlloc2Litr(NE,k_woody_litr)

                LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
                  +ElmAllocmat4Litr(NE,iroot,M,NZ)*dRootMyco*FracRootElmAlloc2Litr(NE,k_fine_litr)
              ENDDO
            ENDDO  
          ENDDO  
        ENDDO
      ENDDO D6415

      DO M=1,jsken
        DO NE=1,NumPlantChemElms  
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=&
             LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*AZMAX1(SeasonalNonstElms_pft(NE,NZ))*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

          LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)= &
             LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*AZMAX1(SeasonalNonstElms_pft(NE,NZ))*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
        ENDDO
      ENDDO
!
      call ResetBranchRootStates(NZ,C4PhotoShootNonstC_brch)
    ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     DazCurrYear=number of days in current year
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. jHarvst_pft(NZ).EQ.jharvtyp_noaction)THEN
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
  end subroutine LiterfallFromRootShootStorage
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadRoots(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,M,NR,N,NE,NTG
!     begin_execution
  associate(                                                         &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr      , & !input  
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr    , & !input  
    Myco_pft                    => plt_morph%Myco_pft                   , & !input  
    MaxNumRootLays            => plt_site%MaxNumRootLays            , & !input  
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft       , & !input  
    NU                        => plt_site%NU                        , & !input  
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft     , & !input  
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft     , & !input  
    SeedDepth_pft             => plt_morph%SeedDepth_pft            , & !input  
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft       , & !input  
    iPlantRootState_pft       => plt_pheno%iPlantRootState_pft      , & !input  
    icwood                    => pltpar%icwood                      , & !input  
    inonstruct                => pltpar%inonstruct                  , & !input  
    iroot                     => pltpar%iroot                       , & !input  
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr      , & !inoput 
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft          , & !inoput 
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft    , & !inoput 
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr , & !inoput 
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr , & !inoput 
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr    , & !inoput 
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr   , & !inoput 
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr   , & !inoput 
    k_fine_litr               => pltpar%k_fine_litr                 , & !inoput 
    k_woody_litr              => pltpar%k_woody_litr                , & !inoput 
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr           , & !inoput 
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr           , & !inoput 
    NIXBotRootLayer_pft       => plt_morph%NIXBotRootLayer_pft      , & !output 
    NIXBotRootLayer_rpft      => plt_morph%NIXBotRootLayer_rpft     , & !output 
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr        , & !output 
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft          , & !output 
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr          , & !output 
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr        , & !output 
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr         , & !output 
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr          , & !output 
    Root2ndMeanLens_pvr       => plt_morph%Root2ndMeanLens_pvr      , & !output 
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr        , & !output 
    Root2ndXNum_pvr           => plt_morph%Root2ndXNum_pvr          , & !output 
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr         , & !output 
    RootAreaPerPlant_pvr      => plt_morph%RootAreaPerPlant_pvr     , & !output 
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr  , & !output 
    RootLenPerPlant_pvr       => plt_morph%RootLenPerPlant_pvr      , & !output 
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs       , & !output 
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr   , & !output 
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr          , & !output 
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr          , & !output 
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr               & !output 
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
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
              +ElmAllocmat4Litr(NE,inonstruct,M,NZ)* RootMycoNonstElms_rpvr(NE,N,L,NZ)
          ENDDO 
        ENDDO
      ENDDO  
    ENDDO  


    DO  NR=1,NumRootAxes_pft(NZ)
      DO L=NU,MaxNumRootLays        
        DO N=1,Myco_pft(NZ)
          DO M=1,jsken
            DO NE=1,NumPlantChemElms
              LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)&
                +ElmAllocmat4Litr(NE,icwood,M,NZ) &
                *(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))&
                *FracRootElmAlloc2Litr(NE,k_woody_litr)

              LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
                +ElmAllocmat4Litr(NE,iroot,M,NZ) &
                *(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))&
                *FracRootElmAlloc2Litr(NE,k_fine_litr)
            enddo
          ENDDO
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
    D8870: DO NR=1,NumRootAxes_pft(NZ)
      DO L=NU,MaxNumRootLays             
        DO N=1,Myco_pft(NZ)        
          RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
          RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
          Root1stLen_rpvr(N,L,NR,NZ)                              = 0._r8
          Root2ndLen_rpvr(N,L,NR,NZ)                              = 0._r8
          Root2ndXNum_rpvr(N,L,NR,NZ)                             = 0._r8
        ENDDO
      ENDDO    
      DO N=1,Myco_pft(NZ)        
        RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)=0._r8      
      ENDDO    
    ENDDO D8870

    DO L=NU,MaxNumRootLays             
      DO N=1,Myco_pft(NZ)           
        RootMycoNonstElms_rpvr(:,N,L,NZ) = 0._r8
        RootMycoActiveBiomC_pvr(N,L,NZ)  = 0._r8
        PopuRootMycoC_pvr(N,L,NZ)        = 0._r8
        RootProteinC_pvr(N,L,NZ)         = 0._r8
        Root1stXNumL_pvr(N,L,NZ)         = 0._r8
        Root2ndXNum_pvr(N,L,NZ)          = 0._r8
        RootLenPerPlant_pvr(N,L,NZ)      = 0._r8
        RootLenDensPerPlant_pvr(N,L,NZ)  = 0._r8
        RootPoreVol_pvr(N,L,NZ)          = 0._r8
        RootVH2O_pvr(N,L,NZ)             = 0._r8
        Root1stRadius_pvr(N,L,NZ)        = Root1stMaxRadius_pft(N,NZ)
        Root2ndRadius_pvr(N,L,NZ)        = Root2ndMaxRadius_pft(N,NZ)
        RootAreaPerPlant_pvr(N,L,NZ)     = 0._r8
        Root2ndMeanLens_pvr(N,L,NZ)        = Root2ndMeanLensMin
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
              LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+&
                ElmAllocmat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_rpvr(NE,L,NZ)+&
                ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_rpvr(NE,L,NZ)
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
    D8795: DO NR=1,NumRootAxes_pft(NZ)
      NIXBotRootLayer_rpft(NR,NZ)=NGTopRootLayer_pft(NZ)
      D8790: DO N=1,Myco_pft(NZ)
        Root1stDepz_pft(N,NR,NZ)=SeedDepth_pft(NZ)
      ENDDO D8790
    ENDDO D8795
    NIXBotRootLayer_pft(NZ) = NGTopRootLayer_pft(NZ)
    NumRootAxes_pft(NZ)     = 0
  ENDIF
  end associate
  end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadBranches(I,J,NZ,NumDeadBranches,C4PhotoShootNonstC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: NumDeadBranches
  real(r8), intent(inout) :: C4PhotoShootNonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: M,NE,NB
!     begin_execution
  associate(                                                                           &
    CanopyNodulNonstElms_brch         => plt_biom%CanopyNodulNonstElms_brch          , & !input  
    CanopyNodulStrutElms_brch         => plt_biom%CanopyNodulStrutElms_brch          , & !input  
    CanopyNonstElms_brch              => plt_biom%CanopyNonstElms_brch               , & !input  
    EarStrutElms_brch                 => plt_biom%EarStrutElms_brch                  , & !input  
    ElmAllocmat4Litr                  => plt_soilchem%ElmAllocmat4Litr               , & !input  
    FracShootLeafElmAlloc2Litr        => plt_allom%FracShootLeafElmAlloc2Litr        , & !input  
    FracShootStalkElmAlloc2Litr       => plt_allom%FracShootStalkElmAlloc2Litr       , & !input  
    GrainStrutElms_brch               => plt_biom%GrainStrutElms_brch                , & !input  
    HuskStrutElms_brch                => plt_biom%HuskStrutElms_brch                 , & !input  
    LeafStrutElms_brch                => plt_biom%LeafStrutElms_brch                 , & !input  
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                   , & !input  
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft                 , & !input  
    PetoleStrutElms_brch              => plt_biom%PetoleStrutElms_brch               , & !input  
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft        , & !input  
    StalkRsrvElms_brch                => plt_biom%StalkRsrvElms_brch                 , & !input  
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch                , & !input  
    iHarvstType_pft                   => plt_distb%iHarvstType_pft                   , & !input  
    iPlantBranchState_brch            => plt_pheno%iPlantBranchState_brch            , & !input  
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft           , & !input  
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft              , & !input  
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft             , & !input  
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft         , & !input  
    icwood                            => pltpar%icwood                               , & !input  
    ifoliar                           => pltpar%ifoliar                              , & !input  
    inonfoliar                        => pltpar%inonfoliar                           , & !input  
    inonstruct                        => pltpar%inonstruct                           , & !input  
    istalk                            => pltpar%istalk                               , & !input  
    LitrfalStrutElms_pvr              => plt_bgcr%LitrfalStrutElms_pvr               , & !inoput 
    SeasonalNonstElms_pft             => plt_biom%SeasonalNonstElms_pft              , & !inoput 
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                 , & !inoput 
    StandDeadKCompElms_pft            => plt_biom%StandDeadKCompElms_pft             , & !inoput 
    k_fine_litr                       => pltpar%k_fine_litr                          , & !inoput 
    k_woody_litr                      => pltpar%k_woody_litr                         , & !inoput 
    BranchNumber_brch                 => plt_morph%BranchNumber_brch                 , & !output 
    C4PhotosynDowreg_brch             => plt_photo%C4PhotosynDowreg_brch             , & !output 
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch            , & !output 
    Hours2LeafOut_brch                => plt_pheno%Hours2LeafOut_brch                , & !output 
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                , & !output 
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch                , & !output 
    Hours4LenthenPhotoPeriod_brch     => plt_pheno%Hours4LenthenPhotoPeriod_brch     , & !output 
    Hours4LiterfalAftMature_brch      => plt_pheno%Hours4LiterfalAftMature_brch      , & !output 
    Hours4ShortenPhotoPeriod_brch     => plt_pheno%Hours4ShortenPhotoPeriod_brch     , & !output 
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch            , & !output 
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch                  , & !output 
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch       , & !output 
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch                  , & !output 
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch           , & !output 
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch         , & !output 
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch                  , & !output 
    Prep4Literfall_brch               => plt_pheno%Prep4Literfall_brch               , & !output 
    RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch              , & !output 
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch , & !output 
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch     , & !output 
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch                , & !output 
    doPlantLeafOut_brch               => plt_pheno%doPlantLeafOut_brch               , & !output 
    doPlantLeaveOff_brch              => plt_pheno%doPlantLeaveOff_brch              , & !output 
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                 & !output 
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
      BranchNumber_brch(NB,NZ)                 = 0
      D8850: DO M=1,pltpar%NumGrowthStages
        iPlantCalendar_brch(M,NB,NZ)=0
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
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*CanopyNodulNonstElms_brch(NE,NB,NZ) &
            +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr) &
            +CanopyNodulStrutElms_brch(NE,NB,NZ)) &
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr) &
            +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))

          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr) &
            +PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr))

          IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
            SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ) &
              +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
          ELSE
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
          ENDIF
          IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. iPlantRootProfile_pft(NZ).LE.1)THEN
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)
          ELSE
            StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)
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
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*StalkRsrvElms_brch(NE,NB,NZ)
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
  end subroutine LiterfallFromDeadBranches

!------------------------------------------------------------------------------------------

  subroutine ResetBranchRootStates(NZ,C4PhotoShootNonstC_brch)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: C4PhotoShootNonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                                                          &
    Myco_pft                    => plt_morph%Myco_pft                   , & !input  
    MaxNumRootLays            => plt_site%MaxNumRootLays            , & !input  
    NU                        => plt_site%NU                        , & !input  
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft        , & !input  
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft          , & !input  
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch , & !output 
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch , & !output 
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch      , & !output 
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch         , & !output 
    GrainStrutElms_brch       => plt_biom%GrainStrutElms_brch       , & !output 
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch        , & !output 
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch    , & !output 
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch        , & !output 
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch      , & !output 
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr          , & !output 
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr          , & !output 
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr         , & !output 
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs       , & !output 
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr , & !output 
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr , & !output 
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr    , & !output 
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft     , & !output 
    SenecStalkStrutElms_brch  => plt_biom%SenecStalkStrutElms_brch  , & !output 
    ShootStrutElms_brch       => plt_biom%ShootStrutElms_brch       , & !output 
    StalkLiveBiomassC_brch    => plt_biom%StalkLiveBiomassC_brch    , & !output 
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch        , & !output 
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch       , & !output 
    iPlantState_pft           => plt_pheno%iPlantState_pft            & !output 
  )
!     RESET BRANCH STATE VARIABLES
!
  
  DO NB=1,NumOfBranches_pft(NZ)
    DO NE=1,NumPlantChemElms
      CanopyNonstElms_brch(NE,NB,NZ)      = 0._r8
      CanopyNodulNonstElms_brch(NE,NB,NZ) = 0._r8
      ShootStrutElms_brch(NE,NB,NZ)       = 0._r8
      LeafStrutElms_brch(NE,NB,NZ)        = 0._r8
      PetoleStrutElms_brch(NE,NB,NZ)      = 0._r8
      StalkStrutElms_brch(NE,NB,NZ)       = 0._r8
      StalkRsrvElms_brch(NE,NB,NZ)        = 0._r8
    ENDDO
  ENDDO

  D8835: DO NB=1,NumOfBranches_pft(NZ)
    C4PhotoShootNonstC_brch(NB,NZ)                           = 0._r8
    StalkLiveBiomassC_brch(NB,NZ)                           = 0._r8
    CanopyNodulStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
    HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
    EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)         = 0._r8
    GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
    LeafPetolBiomassC_brch(NB,NZ)                       = 0._r8
    SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)  = 0._r8
  ENDDO D8835
!
!     RESET ROOT STATE VARIABLES
!
  DO  NR=1,NumRootAxes_pft(NZ)
    DO  N=1,Myco_pft(NZ)
      RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)=0._r8
    ENDDO
  ENDDO

  D6416: DO L=NU,MaxNumRootLays
    DO  N=1,Myco_pft(NZ)
       RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ)=0._r8
      DO  NR=1,NumRootAxes_pft(NZ)
        RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
        RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
        Root1stLen_rpvr(N,L,NR,NZ)                              = 0._r8
        Root2ndLen_rpvr(N,L,NR,NZ)                               = 0._r8
        Root2ndXNum_rpvr(N,L,NR,NZ)                             = 0._r8
      enddo
    enddo
  ENDDO D6416
  SeasonalNonstElms_pft(1:NumPlantChemElms,NZ) = 0._r8
  iPlantState_pft(NZ)                          = iDead
  end associate
  end subroutine ResetBranchRootStates

!------------------------------------------------------------------------------------------

  subroutine ResetDeadRootStates(NB,NZ,C4PhotoShootNonstC_brch)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: C4PhotoShootNonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,K,N
!     begin_execution
  associate(                                                             &
    CanopyLeafAreaZ_pft        => plt_morph%CanopyLeafAreaZ_pft        , & !inoput 
    CanopyLeafArea_lnode       => plt_morph%CanopyLeafArea_lnode       , & !inoput 
    CanopyLeafCLyr_pft         => plt_biom%CanopyLeafCLyr_pft          , & !inoput 
    LeafElmsByLayerNode_brch   => plt_biom%LeafElmsByLayerNode_brch    , & !inoput 
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node  , & !output 
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node , & !output 
    CPOOL3_node                => plt_photo%CPOOL3_node                , & !output 
    CPOOL4_node                => plt_photo%CPOOL4_node                , & !output 
    CanopyNodulNonstElms_brch  => plt_biom%CanopyNodulNonstElms_brch   , & !output 
    CanopyNodulStrutElms_brch  => plt_biom%CanopyNodulStrutElms_brch   , & !output 
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch        , & !output 
    CanopyStalkArea_lbrch      => plt_morph%CanopyStalkArea_lbrch      , & !output 
    EarStrutElms_brch          => plt_biom%EarStrutElms_brch           , & !output 
    GrainSeedBiomCMean_brch    => plt_allom%GrainSeedBiomCMean_brch    , & !output 
    GrainStrutElms_brch        => plt_biom%GrainStrutElms_brch         , & !output 
    HuskStrutElms_brch         => plt_biom%HuskStrutElms_brch          , & !output 
    InternodeHeightDead_brch   => plt_morph%InternodeHeightDead_brch   , & !output 
    InternodeStrutElms_brch    => plt_biom%InternodeStrutElms_brch     , & !output 
    LeafAreaLive_brch          => plt_morph%LeafAreaLive_brch          , & !output 
    LeafAreaZsec_brch          => plt_morph%LeafAreaZsec_brch          , & !output 
    LeafElmntNode_brch         => plt_biom%LeafElmntNode_brch          , & !output 
    LeafNodeArea_brch          => plt_morph%LeafNodeArea_brch          , & !output 
    LeafPetolBiomassC_brch     => plt_biom%LeafPetolBiomassC_brch      , & !output 
    LeafProteinCNode_brch      => plt_biom%LeafProteinCNode_brch       , & !output 
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch          , & !output 
    LiveInterNodeHight_brch    => plt_morph%LiveInterNodeHight_brch    , & !output 
    PetioleElmntNode_brch      => plt_biom%PetioleElmntNode_brch       , & !output 
    PetoleLensNode_brch        => plt_morph%PetoleLensNode_brch        , & !output 
    PetoleProteinCNode_brch    => plt_biom%PetoleProteinCNode_brch     , & !output 
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch        , & !output 
    PotentialSeedSites_brch    => plt_morph%PotentialSeedSites_brch    , & !output 
    SeedNumSet_brch            => plt_morph%SeedNumSet_brch            , & !output 
    SenecStalkStrutElms_brch   => plt_biom%SenecStalkStrutElms_brch    , & !output 
    ShootStrutElms_brch        => plt_biom%ShootStrutElms_brch         , & !output 
    StalkLiveBiomassC_brch     => plt_biom%StalkLiveBiomassC_brch      , & !output 
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch          , & !output 
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch         , & !output 
    StemAreaZsec_brch          => plt_morph%StemAreaZsec_brch            & !output 
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
!     SeedNumSet_brch=seed set number
!     PotentialSeedSites_brch=potential number of seed set sites
!     GrainSeedBiomCMean_brch=individual seed size
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     LeafProteinCNode_brch=leaf protein mass
!     LeafAreaLive_brch=branch leaf area
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinCNode_brch=node petiole C,N,P,protein mass
!     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
  C4PhotoShootNonstC_brch(NB,NZ)                           = 0._r8
  CanopyNonstElms_brch(1:NumPlantChemElms,NB,NZ)      = 0._r8
  CanopyNodulNonstElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
  ShootStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
  LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
  CanopyNodulStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
  PetoleStrutElms_brch(1:NumPlantChemElms,NB,NZ)      = 0._r8
  StalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
  StalkRsrvElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
  HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)        = 0._r8
  EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)         = 0._r8
  GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ)       = 0._r8
  StalkLiveBiomassC_brch(NB,NZ)                           = 0._r8
  LeafPetolBiomassC_brch(NB,NZ)                       = 0._r8
  PotentialSeedSites_brch(NB,NZ)                      = 0._r8
  SeedNumSet_brch(NB,NZ)                              = 0._r8
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
    LeafNodeArea_brch(K,NB,NZ)                          = 0._r8
    LiveInterNodeHight_brch(K,NB,NZ)                    = 0._r8
    InternodeHeightDead_brch(K,NB,NZ)                  = 0._r8
    PetoleLensNode_brch(K,NB,NZ)                        = 0._r8
    LeafProteinCNode_brch(K,NB,NZ)                      = 0._r8
    PetoleProteinCNode_brch(K,NB,NZ)                    = 0._r8
    LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)      = 0._r8
    PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)   = 0._r8
    InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
    D8865: DO L=1,NumOfCanopyLayers1
      CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)-CanopyLeafArea_lnode(L,K,NB,NZ)
      CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)-LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)
      CanopyLeafArea_lnode(L,K,NB,NZ)                         = 0._r8
      LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ) = 0._r8
      IF(K.NE.0)THEN
        D8860: DO N=1,NumOfLeafZenithSectors1
          LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
        ENDDO D8860
      ENDIF
    ENDDO D8865
  ENDDO D8855
  D8875: DO L=1,NumOfCanopyLayers1
    CanopyStalkArea_lbrch(L,NB,NZ)=0._r8
    DO  N=1,NumOfLeafZenithSectors1
      StemAreaZsec_brch(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D8875
  end associate
  end subroutine ResetDeadRootStates


end module LitterFallMod
