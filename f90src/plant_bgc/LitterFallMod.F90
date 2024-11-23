module LitrFallMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  use PlantMathFuncMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: ResetDeadBranch
  contains
!------------------------------------------------------------------------------------------

  subroutine ResetDeadBranch(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NumDeadBranches

!     begin_execution
  associate(                                                      &
    iYearPlantHarvest_pft   => plt_distb%iYearPlantHarvest_pft,   &
    iDayPlantHarvest_pft    => plt_distb%iDayPlantHarvest_pft,    &
    ShootC4NonstC_brch      => plt_biom%ShootC4NonstC_brch,       &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft, &
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch,     &
    iYearCurrent            => plt_site%iYearCurrent,             &
    SolarNoonHour_col       => plt_site%SolarNoonHour_col,        &
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft        &
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
  IF(J.EQ.INT(SolarNoonHour_col) .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND. (iPlantPhenolPattern_pft(NZ).NE.iplt_annual .OR. (I.GE.iDayPlantHarvest_pft(NZ) &
    .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ))))THEN
    NumDeadBranches=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,NumDeadBranches,ShootC4NonstC_brch)

    call SetDeadPlant(i,j,NZ,NumDeadBranches)
!
!     DEAD ROOTS
!
!     LitrFall FROM DEAD ROOTS
!
    call LiterfallFromDeadRoots(I,J,NZ)
!
    call LiterfallFromRootShootStorage(I,J,NZ,ShootC4NonstC_brch)
  ENDIF
  end associate
  end subroutine ResetDeadBranch
!------------------------------------------------------------------------------------------
  subroutine SetDeadPlant(I,J,NZ,NumDeadBranches)
  implicit none
  INTEGER, intent(in) :: i,J
  integer, intent(in) :: NZ
  integer, intent(in) :: NumDeadBranches

  associate(                                                      &
    doInitPlant_pft         => plt_pheno%doInitPlant_pft,         &
    iPlantShootState_pft    => plt_pheno%iPlantShootState_pft,    &
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft,       &
    HoursTooLowPsiCan_pft   => plt_pheno%HoursTooLowPsiCan_pft,   &
    HypoctoHeight_pft       => plt_morph%HypoctoHeight_pft,       &
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft,    &
    CanopyWater_pft         => plt_ew%CanopyWater_pft,            &
    RootElms_pft            => plt_biom%RootElms_pft,             &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft, &
    PlantPopulation_pft     => plt_site%PlantPopulation_pft,      &    
    jHarvst_pft             => plt_distb%jHarvst_pft,             &    
    iPlantRootState_pft     => plt_pheno%iPlantRootState_pft,     &
    QH2OLoss_lnds           => plt_site%QH2OLoss_lnds,            &
    H2OLoss_CumYr_col       => plt_ew%H2OLoss_CumYr_col,          &
    BranchNumber_pft        => plt_morph%BranchNumber_pft         &
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
    QH2OLoss_lnds         = QH2OLoss_lnds+CanopyWater_pft(NZ)
    H2OLoss_CumYr_col     = H2OLoss_CumYr_col+CanopyWater_pft(NZ)
    CanopyWater_pft(NZ)   = 0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     iPlantShootState_pft,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
    IF(SeasonalNonstElms_pft(ielmc,NZ).LT.1.0E-04_r8*RootElms_pft(ielmc,NZ).AND.&
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
  ENDIF
  end associate
  end subroutine SetDeadPlant    
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromRootShootStorage(I,J,NZ,ShootC4NonstC_brch)

  use EcoSIMCtrlDataType, only : iYearCurrent
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: ShootC4NonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,M,NR,NB,N,NE
  REAL(R8) :: dRootMyco
!     begin_execution
  associate(                                                              &
    jHarvst_pft                 => plt_distb%jHarvst_pft,                 &
    iYearPlanting_pft           => plt_distb%iYearPlanting_pft,           &
    iDayPlanting_pft            => plt_distb%iDayPlanting_pft,            &
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch,            &
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch,         &
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch,          &
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch,           &
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch,         &
    CanopyNodulStrutElms_brch   => plt_biom%CanopyNodulStrutElms_brch,    &
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch,           &
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch,    &
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft,        &
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch,          &
    RootMyco1stStrutElms_rpvr   => plt_biom%RootMyco1stStrutElms_rpvr,    &
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch,           &
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr,       &
    StandDeadKCompElms_pft      => plt_biom%StandDeadKCompElms_pft,       &
    RootMyco2ndStrutElms_rpvr   => plt_biom%RootMyco2ndStrutElms_rpvr,    &
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr,  &
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr,  &
    FracRootElmAlloc2Litr       => plt_allom%FracRootElmAlloc2Litr,       &
    doInitPlant_pft             => plt_pheno%doInitPlant_pft,             &
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft,     &
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft,        &
    iPlantRootState_pft         => plt_pheno%iPlantRootState_pft,         &
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft,       &
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft,   &
    iPlantShootState_pft        => plt_pheno%iPlantShootState_pft,        &
    icwood                      => pltpar%icwood,                         &
    ifoliar                     => pltpar%ifoliar,                        &
    k_fine_litr                 => pltpar%k_fine_litr,                    &
    k_woody_litr                => pltpar%k_woody_litr,                   &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,         &
    inonfoliar                  => pltpar%inonfoliar,                     &
    istalk                      => pltpar%istalk,                         &
    iroot                       => pltpar%iroot,                          &
    inonstruct                  => pltpar%inonstruct,                     &
    DazCurrYear                 => plt_site%DazCurrYear,                  &
    MaxNumRootLays              => plt_site%MaxNumRootLays,               &
    NU                          => plt_site%NU,                           &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,         &
    MY                          => plt_morph%MY,                          &
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft,          &
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft,           &
    NumRootAxes_pft             => plt_morph%NumRootAxes_pft              &
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
!     ShootC4NonstC_brch=total C4 nonstructural C in branch
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
            +ElmAllocmat4Litr(ielmc,inonstruct,M,NZ)*AZMAX1(ShootC4NonstC_brch(NB,NZ))

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
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        
      D6415: DO L=NU,MaxNumRootLays
        DO N=1,MY(NZ)
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
      call ResetBranchRootStates(NZ,ShootC4NonstC_brch)
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
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,       &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,        &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,   &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,          &
    RootNodulNonstElms_rpvr    => plt_biom%RootNodulNonstElms_rpvr,    &
    RootNodulStrutElms_rpvr    => plt_biom%RootNodulStrutElms_rpvr,    &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,    &
    iPlantRootState_pft       => plt_pheno%iPlantRootState_pft,      &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,           &
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,           &
    MaxNumRootLays            => plt_site%MaxNumRootLays,            &
    NU                        => plt_site%NU,                        &
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,    &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    icwood                    => pltpar%icwood,                      &
    iroot                     => pltpar%iroot,                       &
    inonstruct                => pltpar%inonstruct,                  &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    k_woody_litr              => pltpar%k_woody_litr,                &
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr,  &
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr,             &
    RootAreaPerPlant_pvr      => plt_morph%RootAreaPerPlant_pvr,     &
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr,          &
    Root2ndXNum_pvr           => plt_morph%Root2ndXNum_pvr,          &
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr,         &
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft,          &
    Root2ndAveLen_pvr         => plt_morph%Root2ndAveLen_pvr,        &
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr,        &
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr,        &
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft,     &
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft,     &
    RootLenPerPlant_pvr       => plt_morph%RootLenPerPlant_pvr,      &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    Root2ndLen_rpvr            => plt_morph%Root2ndLen_rpvr,           &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,       &
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr,         &
    MY                        => plt_morph%MY,                       &
    NIXBotRootLayer_pft       => plt_morph%NIXBotRootLayer_pft,      &
    NIXBotRootLayer_rpft      => plt_morph%NIXBotRootLayer_rpft,     &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,            &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,       &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft           &
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
      DO N=1,MY(NZ)
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
        DO N=1,MY(NZ)
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
      DO N=1,MY(NZ)
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        DO NTG=idg_beg,idg_end-1
          RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-trcg_rootml_pvr(NTG,N,L,NZ)&
            -trcs_rootml_pvr(NTG,N,L,NZ)
        ENDDO
        trcg_rootml_pvr(idg_beg:idg_end-1,N,L,NZ)=0._r8
        trcs_rootml_pvr(idg_beg:idg_end-1,N,L,NZ)=0._r8
      ENDDO
    ENDDO    
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stLen_rpvr,Root2ndLen_rpvr=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootMycoActiveBiomC_pvr,WTRTD=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,Root2ndXNum_pvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_pvr,RootPoreVol_pvr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_pvr=root surface area per plant
!     Root2ndAveLen_pvr=average secondary root length
!
    D8870: DO NR=1,NumRootAxes_pft(NZ)
      DO L=NU,MaxNumRootLays             
        DO N=1,MY(NZ)        
          RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
          RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
          Root1stLen_rpvr(N,L,NR,NZ)                              = 0._r8
          Root2ndLen_rpvr(N,L,NR,NZ)                               = 0._r8
          Root2ndXNum_rpvr(N,L,NR,NZ)                             = 0._r8
        ENDDO
      ENDDO    
      DO N=1,MY(NZ)        
        RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)=0._r8      
      ENDDO    
    ENDDO D8870

    DO L=NU,MaxNumRootLays             
      DO N=1,MY(NZ)           
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
        Root2ndAveLen_pvr(N,L,NZ)        = Root2ndAveLenMin
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
    DO N=1,MY(NZ)       
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
!     NIXBotRootLayer_rpft=deepest root layer
!     RTDP1=primary root depth from soil surface
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!   
    D8795: DO NR=1,NumRootAxes_pft(NZ)
      NIXBotRootLayer_rpft(NR,NZ)=NGTopRootLayer_pft(NZ)
      D8790: DO N=1,MY(NZ)
        Root1stDepz_pft(N,NR,NZ)=SeedDepth_pft(NZ)
      ENDDO D8790
    ENDDO D8795
    NIXBotRootLayer_pft(NZ) = NGTopRootLayer_pft(NZ)
    NumRootAxes_pft(NZ)     = 0
  ENDIF
  end associate
  end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadBranches(I,J,NZ,NumDeadBranches,ShootC4NonstC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: NumDeadBranches
  real(r8), intent(inout) :: ShootC4NonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: M,NE,NB
!     begin_execution
  associate(                                                                          &
    iHarvstType_pft                   => plt_distb%iHarvstType_pft,                   &
    FracShootLeafElmAlloc2Litr        => plt_allom%FracShootLeafElmAlloc2Litr,        &
    FracShootStalkElmAlloc2Litr       => plt_allom%FracShootStalkElmAlloc2Litr,       &
    CanopyNodulNonstElms_brch         => plt_biom%CanopyNodulNonstElms_brch,          &
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch,                &
    HuskStrutElms_brch                => plt_biom%HuskStrutElms_brch,                 &
    PetoleStrutElms_brch              => plt_biom%PetoleStrutElms_brch,               &
    CanopyNodulStrutElms_brch         => plt_biom%CanopyNodulStrutElms_brch,          &
    LeafStrutElms_brch                => plt_biom%LeafStrutElms_brch,                 &
    ZERO4Groth_pft                    => plt_biom%ZERO4Groth_pft,                     &
    EarStrutElms_brch                 => plt_biom%EarStrutElms_brch,                  &
    CanopyNonstElms_brch              => plt_biom%CanopyNonstElms_brch,               &
    GrainStrutElms_brch               => plt_biom%GrainStrutElms_brch,                &
    StalkRsrvElms_brch                => plt_biom%StalkRsrvElms_brch,                 &
    StandDeadKCompElms_pft            => plt_biom%StandDeadKCompElms_pft,             &
    SeasonalNonstElms_pft             => plt_biom%SeasonalNonstElms_pft,              &
    ElmAllocmat4Litr                  => plt_soilchem%ElmAllocmat4Litr,               &
    iPlantBranchState_brch            => plt_pheno%iPlantBranchState_brch,            &
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch,                  &
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch,       &
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch,            &
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch,     &
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch, &
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch,                &
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch,                &
    Hours4LenthenPhotoPeriod_brch     => plt_pheno%Hours4LenthenPhotoPeriod_brch,     &
    Hours4ShortenPhotoPeriod_brch     => plt_pheno%Hours4ShortenPhotoPeriod_brch,     &
    Hours2LeafOut_brch                => plt_pheno%Hours2LeafOut_brch,                &
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch,            &
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch,                &
    doPlantLeafOut_brch               => plt_pheno%doPlantLeafOut_brch,               &
    Prep4Literfall_brch               => plt_pheno%Prep4Literfall_brch,               &
    Hours4LiterfalAftMature_brch      => plt_pheno%Hours4LiterfalAftMature_brch,      &
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft,                   &
    doPlantLeaveOff_brch              => plt_pheno%doPlantLeaveOff_brch,              &
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch,               &
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft,         &
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft,             &
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft,              &
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft,           &
    LitrfalStrutElms_pvr              => plt_bgcr%LitrfalStrutElms_pvr,               &
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft,                 &
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch,           &
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch,                 &
    BranchNumber_brch                 => plt_morph%BranchNumber_brch,                 &
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch,         &
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft,        &
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch,                  &
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch,                  &
    istalk                            => pltpar%istalk,                               &
    k_fine_litr                       => pltpar%k_fine_litr,                          &
    k_woody_litr                      => pltpar%k_woody_litr,                         &
    inonstruct                        => pltpar%inonstruct,                           &
    icwood                            => pltpar%icwood,                               &
    inonfoliar                        => pltpar%inonfoliar,                           &
    ifoliar                           => pltpar%ifoliar,                              &
    RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch,              &
    C4PhotosynDowreg_brch             => plt_photo%C4PhotosynDowreg_brch              &
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
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
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
          IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantRootProfile_pft(NZ).LE.1)THEN
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
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     ShootC4NonstC_brch=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
      SeasonalNonstElms_pft(ielmc,NZ)=SeasonalNonstElms_pft(ielmc,NZ)+ShootC4NonstC_brch(NB,NZ)

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

      call ResetDeadRootStates(NB,NZ,ShootC4NonstC_brch)

      NumDeadBranches=NumDeadBranches+1
    ENDIF
  ENDDO D8845
  end associate
  end subroutine LiterfallFromDeadBranches

!------------------------------------------------------------------------------------------

  subroutine ResetBranchRootStates(NZ,ShootC4NonstC_brch)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: ShootC4NonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                                                         &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch, &
    ShootStrutElms_brch       => plt_biom%ShootStrutElms_brch,       &
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,        &
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch,        &
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,      &
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch,        &
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch,         &
    GrainStrutElms_brch       => plt_biom%GrainStrutElms_brch,       &
    StalkBiomassC_brch        => plt_biom%StalkBiomassC_brch,        &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    SenecStalkStrutElms_brch  => plt_biom%SenecStalkStrutElms_brch,  &
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,    &
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,     &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch, &
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch,       &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,       &
    MaxNumRootLays            => plt_site%MaxNumRootLays,            &
    NU                        => plt_site%NU,                        &
    iPlantState_pft           => plt_pheno%iPlantState_pft,          &
    Root2ndLen_rpvr            => plt_morph%Root2ndLen_rpvr,           &
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr,         &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    MY                        => plt_morph%MY,                       &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,        &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft           &
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
    ShootC4NonstC_brch(NB,NZ)                           = 0._r8
    StalkBiomassC_brch(NB,NZ)                           = 0._r8
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
    DO  N=1,MY(NZ)
      RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)=0._r8
    ENDDO
  ENDDO

  D6416: DO L=NU,MaxNumRootLays
    DO  N=1,MY(NZ)
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

  subroutine ResetDeadRootStates(NB,NZ,ShootC4NonstC_brch)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: ShootC4NonstC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,K,N
!     begin_execution
  associate(                                                            &
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch,        &
    CanopyNodulNonstElms_brch  => plt_biom%CanopyNodulNonstElms_brch,   &
    ShootStrutElms_brch        => plt_biom%ShootStrutElms_brch,         &
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch,          &
    CanopyNodulStrutElms_brch  => plt_biom%CanopyNodulStrutElms_brch,   &
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch,         &
    LeafElmntNode_brch         => plt_biom%LeafElmntNode_brch,          &
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch,          &
    SenecStalkStrutElms_brch   => plt_biom%SenecStalkStrutElms_brch,    &
    StalkBiomassC_brch         => plt_biom%StalkBiomassC_brch,          &
    GrainStrutElms_brch        => plt_biom%GrainStrutElms_brch,         &
    LeafPetolBiomassC_brch     => plt_biom%LeafPetolBiomassC_brch,      &
    HuskStrutElms_brch         => plt_biom%HuskStrutElms_brch,          &
    EarStrutElms_brch          => plt_biom%EarStrutElms_brch,           &
    PetoleProteinCNode_brch    => plt_biom%PetoleProteinCNode_brch,     &
    LeafProteinCNode_brch      => plt_biom%LeafProteinCNode_brch,       &
    InternodeStrutElms_brch    => plt_biom%InternodeStrutElms_brch,     &
    PetioleElmntNode_brch      => plt_biom%PetioleElmntNode_brch,       &
    LeafElmsByLayerNode_brch   => plt_biom%LeafElmsByLayerNode_brch,    &
    CanopyLeafCLyr_pft         => plt_biom%CanopyLeafCLyr_pft,          &
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch,        &
    CPOOL3_node                => plt_photo%CPOOL3_node,                &
    CPOOL4_node                => plt_photo%CPOOL4_node,                &
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node,  &
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node, &
    GrainSeedBiomCMean_brch    => plt_allom%GrainSeedBiomCMean_brch,    &
    PotentialSeedSites_brch    => plt_morph%PotentialSeedSites_brch,    &
    SeedNumSet_brch            => plt_morph%SeedNumSet_brch,            &
    LeafAreaLive_brch          => plt_morph%LeafAreaLive_brch,          &
    CanopyStalkArea_lbrch      => plt_morph%CanopyStalkArea_lbrch,      &
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft,          &
    LeafAreaNode_brch          => plt_morph%LeafAreaNode_brch,          &
    PetoleLensNode_brch        => plt_morph%PetoleLensNode_brch,        &
    InternodeHeightDying_brch  => plt_morph%InternodeHeightDying_brch,  &
    LeafAreaZsec_brch          => plt_morph%LeafAreaZsec_brch,          &
    LiveInterNodeHight_brch    => plt_morph%LiveInterNodeHight_brch,    &
    CanopyLeafAreaZ_pft        => plt_morph%CanopyLeafAreaZ_pft,        &
    StemAreaZsec_brch          => plt_morph%StemAreaZsec_brch,          &
    CanopyLeafArea_lpft        => plt_morph%CanopyLeafArea_lpft         &
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
  ShootC4NonstC_brch(NB,NZ)                           = 0._r8
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
  StalkBiomassC_brch(NB,NZ)                           = 0._r8
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
    LeafAreaNode_brch(K,NB,NZ)                          = 0._r8
    LiveInterNodeHight_brch(K,NB,NZ)                    = 0._r8
    InternodeHeightDying_brch(K,NB,NZ)                  = 0._r8
    PetoleLensNode_brch(K,NB,NZ)                        = 0._r8
    LeafProteinCNode_brch(K,NB,NZ)                      = 0._r8
    PetoleProteinCNode_brch(K,NB,NZ)                    = 0._r8
    LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)      = 0._r8
    PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)   = 0._r8
    InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
    D8865: DO L=1,NumOfCanopyLayers1
      CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)-CanopyLeafArea_lpft(L,K,NB,NZ)
      CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)-LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)
      CanopyLeafArea_lpft(L,K,NB,NZ)                         = 0._r8
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


end module LitrFallMod
