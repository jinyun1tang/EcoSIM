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

  subroutine ResetDeadBranch(I,J,NZ,ShootNonstructC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8),intent(inout) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: IDTHY

!     begin_execution
  associate(                                 &
    iYearPlantHarvest_pft       =>   plt_distb%iYearPlantHarvest_pft     , &
    jHarvst_pft                 =>   plt_distb%jHarvst_pft    , &
    iDayPlantHarvest_pft        =>   plt_distb%iDayPlantHarvest_pft    , &
    UVOLO                       =>   plt_ew%UVOLO       , &
    CanopyWater_pft             =>   plt_ew%CanopyWater_pft       , &
    RootElmnts_pft              =>   plt_biom%RootElmnts_pft     , &
    NonstructalElmnts_pft       =>   plt_biom%NonstructalElmnts_pft     , &
    iPlantPhenologyPattern_pft  =>   plt_pheno%iPlantPhenologyPattern_pft   , &
    iPlantRootState_pft         =>   plt_pheno%iPlantRootState_pft   , &
    iPlantShootState_pft        =>   plt_pheno%iPlantShootState_pft    , &
    doInitPlant_pft             =>   plt_pheno%doInitPlant_pft    , &
    HoursCanopyPSITooLow        =>   plt_pheno%HoursCanopyPSITooLow     , &
    iPlantCalendar_brch         =>   plt_pheno%iPlantCalendar_brch   , &
    PlantPopulation_pft         =>   plt_site%PlantPopulation_pft        , &
    iYearCurrent                =>   plt_site%iYearCurrent      , &
    VOLWOU                      =>   plt_site%VOLWOU    , &
    SolarNoonHour_col           =>   plt_site%SolarNoonHour_col    , &
    HypoctoHeight_pft           =>   plt_morph%HypoctoHeight_pft    , &
    NumOfBranches_pft           =>   plt_morph%NumOfBranches_pft      , &
    NumOfMainBranch_pft         =>   plt_morph%NumOfMainBranch_pft      , &
    BranchNumber_pft            =>   plt_morph%BranchNumber_pft        &
  )
!
!     SolarNoonHour_col=hour of solar noon
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iDayPlantHarvest_pft,iYearPlantHarvest_pft=day,year of harvesting
!     iYearCurrent=current year
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     GROUP=node number required for floral initiation
!     NodeNumberToInitFloral_brch=node number at floral initiation
!     NodeNumberAtAnthesis_brch=node number at flowering
!     NumOfLeaves_brch=number of leaves appeared
!     KLeafNodeNumber=integer of most recent leaf number currently growing
!     VSTGX=leaf number on date of floral initiation
!     TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
!     TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
!     HourFailGrainFill_brch=number of hours with no grain fill
!     doInitLeafOut_brch=flag for initializing leafout
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     HourCounter4LeafOut_brch=hourly leafout counter
!     RubiscoActivity_brch,C4PhotosynDowreg_brch=N,P feedback inhibition on C3 CO2 fixation
!     doInitLeafOut_brch,doPlantLeafOut_brch=flag for initializing,enabling leafout
!     doPlantLeaveOff_brch=flag for enabling leafoff:0=enable,1=disable
!     Hours4LiterfalAftMature_brch=current hours after physl maturity until start of LitrFall
!
  IF(J.EQ.INT(SolarNoonHour_col).AND.iPlantCalendar_brch(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).NE.0 &
    .AND.(iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.OR.(I.GE.iDayPlantHarvest_pft(NZ) &
    .AND.iYearCurrent.GE.iYearPlantHarvest_pft(NZ))))THEN
    IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,IDTHY,ShootNonstructC_brch)

    IF(IDTHY.EQ.NumOfBranches_pft(NZ))THEN
      iPlantShootState_pft(NZ)=iDead
      BranchNumber_pft(NZ)=0
      HoursCanopyPSITooLow(NZ)=0._r8
      IF(doInitPlant_pft(NZ).EQ.itrue)THEN
        NumOfBranches_pft(NZ)=1
      ELSE
        NumOfBranches_pft(NZ)=0
      ENDIF
      HypoctoHeight_pft(NZ)=0._r8
      VOLWOU=VOLWOU+CanopyWater_pft(NZ)
      UVOLO=UVOLO+CanopyWater_pft(NZ)
      CanopyWater_pft(NZ)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     iPlantShootState_pft,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(NonstructalElmnts_pft(ielmc,NZ).LT.1.0E-04_r8*RootElmnts_pft(ielmc,NZ).AND.&
        iPlantPhenologyPattern_pft(NZ).NE.iplt_annual)then
        iPlantRootState_pft(NZ)=iDead
      endif
      IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual)then
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
!
!     DEAD ROOTS
!
!
!     LitrFall FROM DEAD ROOTS
!
    call LiterfallFromDeadRoots(I,J,NZ)
!
    call LiterfallFromRootShootStorage(I,J,NZ,ShootNonstructC_brch)
  ENDIF
  end associate
  end subroutine ResetDeadBranch
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromRootShootStorage(I,J,NZ,ShootNonstructC_brch)

  use EcoSIMCtrlDataType, only : iYearCurrent
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,M,NR,NB,N,NE
!     begin_execution
  associate(                            &
    jHarvst_pft                      =>   plt_distb%jHarvst_pft   , &
    iYearPlanting_pft                =>   plt_distb%iYearPlanting_pft    , &
    iDayPlanting_pft                 =>   plt_distb%iDayPlanting_pft   , &
    EarChemElms_brch               =>   plt_biom%EarChemElms_brch  , &
    NonstructElmnt_brch              =>   plt_biom%NonstructElmnt_brch   , &
    StalkChemElms_brch             =>   plt_biom%StalkChemElms_brch  , &
    HuskChemElms_brch              =>   plt_biom%HuskChemElms_brch  , &
    PetoleChemElm_brch             =>   plt_biom%PetoleChemElm_brch , &
    CanopyNoduleChemElm_brch       =>   plt_biom%CanopyNoduleChemElm_brch   , &
    ReserveElmnts_brch               =>   plt_biom%ReserveElmnts_brch  , &
    NoduleNonstructElmnt_brch        =>   plt_biom%NoduleNonstructElmnt_brch   , &
    NonstructalElmnts_pft            =>   plt_biom%NonstructalElmnts_pft    , &
    GrainChemElms_brch             =>   plt_biom%GrainChemElms_brch   , &
    Root1stStructChemElm_pvr       =>   plt_biom%Root1stStructChemElm_pvr   , &
    LeafChemElms_brch              =>   plt_biom%LeafChemElms_brch  , &
     RootMycoNonstructElmnt_vr       =>   plt_biom%RootMycoNonstructElmnt_vr   , &
    StandingDeadKCompChemElms_pft  =>   plt_biom%StandingDeadKCompChemElms_pft   , &
    Root2ndStructChemElm_pvr       =>   plt_biom%Root2ndStructChemElm_pvr   , &
    FWODLE                           =>   plt_allom%FWODLE  , &
    FWODBE                           =>   plt_allom%FWODBE  , &
    FWOODE                           =>   plt_allom%FWOODE  , &
    FWODRE                           =>   plt_allom%FWODRE  , &
    doInitPlant_pft                  =>   plt_pheno%doInitPlant_pft   , &
    iPlantPhenologyPattern_pft       =>   plt_pheno%iPlantPhenologyPattern_pft  , &
    iPlantPhenolType_pft          =>   plt_pheno%iPlantPhenolType_pft  , &
    iPlantRootState_pft              =>   plt_pheno%iPlantRootState_pft  , &
    iPlantMorphologyType_pft         =>   plt_pheno%iPlantMorphologyType_pft  , &
    iPlantTurnoverPattern_pft        =>   plt_pheno%iPlantTurnoverPattern_pft  , &
    iPlantShootState_pft             =>   plt_pheno%iPlantShootState_pft   , &
    icwood                           =>   pltpar%icwood     , &
    ifoliar                          =>   pltpar%ifoliar    , &
    k_fine_litr                      =>   pltpar%k_fine_litr, &
    k_woody_litr                     =>  pltpar%k_woody_litr, &
    LitrFallChemElm_pvr        =>   plt_bgcr%LitrFallChemElm_pvr     , &
    inonfoliar                       =>   pltpar%inonfoliar   , &
    istalk                           =>   pltpar%istalk     , &
    iroot                            =>   pltpar%iroot      , &
    inonstruct                       =>   pltpar%inonstruct   , &
    LYRC                             =>   plt_site%LYRC     , &
    MaxNumRootLays                   =>   plt_site%MaxNumRootLays       , &
    NU                               =>   plt_site%NU       , &
    CFOPE                            =>   plt_soilchem%CFOPE, &
    MY                               =>   plt_morph%MY      , &
    NGTopRootLayer_pft               =>   plt_morph%NGTopRootLayer_pft      , &
    NumOfBranches_pft                =>   plt_morph%NumOfBranches_pft     , &
    NumRootAxes_pft                  =>   plt_morph%NumRootAxes_pft      &
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
!     ShootNonstructC_brch=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
  IF(iPlantShootState_pft(NZ).EQ.iDead.AND.iPlantRootState_pft(NZ).EQ.iDead)THEN
    !both plant shoots and roots are dead
    IF(doInitPlant_pft(NZ).EQ.ifalse)THEN
      D6425: DO M=1,jsken
        D8825: DO NB=1,NumOfBranches_pft(NZ)

          LitrFallChemElm_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(ielmc,M,k_fine_litr,0,NZ) &
            +CFOPE(ielmc,inonstruct,M,NZ)*ShootNonstructC_brch(NB,NZ)

          DO NE=1,NumPlantChemElms
            LitrFallChemElm_pvr(NE,M,k_woody_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_woody_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(LeafChemElms_brch(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
              +PetoleChemElm_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

            LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,inonstruct,M,NZ)*(NonstructElmnt_brch(NE,NB,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ) &
              +ReserveElmnts_brch(NE,NB,NZ)) &
              +CFOPE(NE,ifoliar,M,NZ)*(LeafChemElms_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
              +CanopyNoduleChemElm_brch(NE,NB,NZ)) &
              +CFOPE(NE,inonfoliar,M,NZ)*(PetoleChemElm_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
              +HuskChemElms_brch(NE,NB,NZ)+EarChemElms_brch(NE,NB,NZ))

            IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
              NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)+CFOPE(NE,inonfoliar,M,NZ) &
                *GrainChemElms_brch(NE,NB,NZ)
            ELSE
              LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
                +CFOPE(NE,inonfoliar,M,NZ)*GrainChemElms_brch(NE,NB,NZ)
            ENDIF
            IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantMorphologyType_pft(NZ).LE.1)THEN
              !all above ground
              LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)&
                +CFOPE(NE,istalk,M,NZ)*StalkChemElms_brch(NE,NB,NZ)
            ELSE
              StandingDeadKCompChemElms_pft(NE,M,NZ)=StandingDeadKCompChemElms_pft(NE,M,NZ)&
                +CFOPE(NE,icwood,M,NZ)*StalkChemElms_brch(NE,NB,NZ)
            ENDIF
          ENDDO
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
        DO NE=1,NumPlantChemElms
          D6415: DO L=NU,MaxNumRootLays
            DO N=1,MY(NZ)
              LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ) &
                +CFOPE(NE,inonstruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
              DO NR=1,NumRootAxes_pft(NZ)
                LitrFallChemElm_pvr(NE,M,k_woody_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_woody_litr,L,NZ)&
                  +CFOPE(NE,icwood,M,NZ) &
                  *(Root1stStructChemElm_pvr(NE,N,L,NR,NZ)+Root2ndStructChemElm_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

                LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                  *(Root1stStructChemElm_pvr(NE,N,L,NR,NZ)+Root2ndStructChemElm_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
              ENDDO
            ENDDO
          ENDDO D6415
          LitrFallChemElm_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=LitrFallChemElm_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
            +CFOPE(NE,inonstruct,M,NZ)*NonstructalElmnts_pft(NE,NZ)*FWOODE(NE,k_woody_litr)

          LitrFallChemElm_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
            +CFOPE(NE,inonstruct,M,NZ)*NonstructalElmnts_pft(NE,NZ)*FWOODE(NE,k_fine_litr)
        ENDDO
      ENDDO D6425
!
      call ResetBranchRootStates(NZ,ShootNonstructC_brch)
    ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     LYRC=number of days in current year
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!
    IF(iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.jHarvst_pft(NZ).EQ.jharvtyp_noaction)THEN
      IF(I.LT.LYRC)THEN
        iDayPlanting_pft(NZ)=I+1
        iYearPlanting_pft(NZ)=iYearCurrent
      ELSE
        iDayPlanting_pft(NZ)=1
        iYearPlanting_pft(NZ)=iYearCurrent+1
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
  associate(                                &
    Root2ndStructChemElm_pvr    =>   plt_biom%Root2ndStructChemElm_pvr    , &
    Root1stChemElm              =>   plt_biom%Root1stChemElm    , &
    Root1stStructChemElm_pvr    =>   plt_biom%Root1stStructChemElm_pvr    , &
     RootMycoNonstructElmnt_vr    =>   plt_biom%RootMycoNonstructElmnt_vr    , &
     PopuPlantRootC_vr            =>   plt_biom% PopuPlantRootC_vr    , &
    RootStructBiomC_vr            =>   plt_biom%RootStructBiomC_vr    , &
    RootProteinC_pvr              =>   plt_biom%RootProteinC_pvr     , &
    RootNoduleNonstructElmnt_vr   =>   plt_biom%RootNoduleNonstructElmnt_vr   , &
    RootNodueChemElm_pvr        =>   plt_biom%RootNodueChemElm_pvr    , &
    FWODRE                        =>   plt_allom%FWODRE   , &
    iPlantRootState_pft           =>   plt_pheno%iPlantRootState_pft   , &
    CFOPE                         =>   plt_soilchem%CFOPE , &
    trcg_rootml_vr                =>   plt_rbgc%trcg_rootml_vr  , &
    trcs_rootml_vr                => plt_rbgc%trcs_rootml_vr, &
    MaxNumRootLays                =>   plt_site%MaxNumRootLays        , &
    NU                            =>   plt_site%NU        , &
    RootGasLossDisturb_pft        =>   plt_bgcr%RootGasLossDisturb_pft     , &
    LitrFallChemElm_pvr     =>   plt_bgcr%LitrFallChemElm_pvr      , &
    icwood                        =>   pltpar%icwood      , &
    iroot                         =>   pltpar%iroot       , &
    inonstruct                    =>   pltpar%inonstruct   , &
    k_fine_litr                   =>  pltpar%k_fine_litr , &
    k_woody_litr                  => pltpar%k_woody_litr, &
    RootLenDensPerPlant_pvr       =>   plt_morph%RootLenDensPerPlant_pvr    , &
    RootVH2O_vr                   =>   plt_morph%RootVH2O_vr   , &
    RootAreaPerPlant_vr           =>   plt_morph%RootAreaPerPlant_vr    , &
    RootVolume_vr                 =>   plt_morph%RootVolume_vr    , &
    SecndRootXNum_pvr             =>   plt_morph%SecndRootXNum_pvr     , &
    PrimRootXNumL_pvr             =>   plt_morph%PrimRootXNumL_pvr    , &
    PrimRootDepth                 =>   plt_morph%PrimRootDepth   , &
    AveSecndRootLen               =>   plt_morph%AveSecndRootLen    , &
    PrimRootRadius_pvr            =>   plt_morph%PrimRootRadius_pvr    , &
    SecndRootRadius_pvr           =>   plt_morph%SecndRootRadius_pvr    , &
    Max1stRootRadius              =>   plt_morph%Max1stRootRadius   , &
    Max2ndRootRadius              =>   plt_morph%Max2ndRootRadius   , &
    RootLenPerPlant_pvr           =>   plt_morph%RootLenPerPlant_pvr    , &
    PrimRootLen                   =>   plt_morph%PrimRootLen    , &
    SecndRootLen                  =>   plt_morph%SecndRootLen    , &
    iPlantNfixType                =>   plt_morph%iPlantNfixType   , &
    SecndRootXNum_rpvr            =>   plt_morph%SecndRootXNum_rpvr    , &
    MY                            =>   plt_morph%MY       , &
    NIXBotRootLayer_pft           =>   plt_morph%NIXBotRootLayer_pft     , &
    NIXBotRootLayer_rpft          =>   plt_morph%NIXBotRootLayer_rpft     , &
    SeedDepth_pft                 =>   plt_morph%SeedDepth_pft    , &
    NGTopRootLayer_pft            =>   plt_morph%NGTopRootLayer_pft       , &
    NumRootAxes_pft               =>   plt_morph%NumRootAxes_pft       &
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
    D8900: DO N=1,MY(NZ)
      D8895: DO L=NU,MaxNumRootLays
        DO NE=1,NumPlantChemElms
          D6410: DO M=1,jsken
            LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ) &
              +CFOPE(NE,inonstruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
            DO  NR=1,NumRootAxes_pft(NZ)
              LitrFallChemElm_pvr(NE,M,k_woody_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
                *(Root1stStructChemElm_pvr(NE,N,L,NR,NZ)+Root2ndStructChemElm_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

              LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                *(Root1stStructChemElm_pvr(NE,N,L,NR,NZ)+Root2ndStructChemElm_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
            enddo
          ENDDO D6410
        ENDDO
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        DO NTG=idg_beg,idg_end-1
          RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-trcg_rootml_vr(NTG,N,L,NZ)&
            -trcs_rootml_vr(NTG,N,L,NZ)
        ENDDO
        trcg_rootml_vr(idg_beg:idg_end-1,N,L,NZ)=0._r8
        trcs_rootml_vr(idg_beg:idg_end-1,N,L,NZ)=0._r8
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     PrimRootLen,SecndRootLen=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootStructBiomC_vr,WTRTD=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,SecndRootXNum_pvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_vr,RootVolume_vr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_vr=root surface area per plant
!     AveSecndRootLen=average secondary root length
!
        D8870: DO NR=1,NumRootAxes_pft(NZ)
          Root1stStructChemElm_pvr(1:NumPlantChemElms,N,L,NR,NZ)=0._r8
          Root2ndStructChemElm_pvr(1:NumPlantChemElms,N,L,NR,NZ)=0._r8
          Root1stChemElm(1:NumPlantChemElms,N,NR,NZ)=0._r8
          PrimRootLen(N,L,NR,NZ)=0._r8
          SecndRootLen(N,L,NR,NZ)=0._r8
          SecndRootXNum_rpvr(N,L,NR,NZ)=0._r8
        ENDDO D8870
         RootMycoNonstructElmnt_vr(:,N,L,NZ)=0._r8
        RootStructBiomC_vr(N,L,NZ)=0._r8
         PopuPlantRootC_vr(N,L,NZ)=0._r8
        RootProteinC_pvr(N,L,NZ)=0._r8
        PrimRootXNumL_pvr(N,L,NZ)=0._r8
        SecndRootXNum_pvr(N,L,NZ)=0._r8
        RootLenPerPlant_pvr(N,L,NZ)=0._r8
        RootLenDensPerPlant_pvr(N,L,NZ)=0._r8
        RootVolume_vr(N,L,NZ)=0._r8
        RootVH2O_vr(N,L,NZ)=0._r8
        PrimRootRadius_pvr(N,L,NZ)=Max1stRootRadius(N,NZ)
        SecndRootRadius_pvr(N,L,NZ)=Max2ndRootRadius(N,NZ)
        RootAreaPerPlant_vr(N,L,NZ)=0._r8
        AveSecndRootLen(N,L,NZ)=MinAve2ndRootLen
!
!     LitrFall AND STATE VARIABLES FROM DEAD NODULES
!
!     iPlantNfixType=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(is_plant_N2fix(iPlantNfixType(NZ)).AND.N.EQ.ipltroot)THEN
          D6420: DO M=1,jsken
            DO NE=1,NumPlantChemElms            
              LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,L,NZ)+&
                CFOPE(NE,iroot,M,NZ)*RootNodueChemElm_pvr(NE,L,NZ)+&
                CFOPE(NE,inonstruct,M,NZ)*RootNoduleNonstructElmnt_vr(NE,L,NZ)
            ENDDO
          ENDDO D6420
          RootNodueChemElm_pvr(1:NumPlantChemElms,L,NZ)=0._r8
          RootNoduleNonstructElmnt_vr(1:NumPlantChemElms,L,NZ)=0._r8          
        ENDIF
      ENDDO D8895
    ENDDO D8900
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
        PrimRootDepth(N,NR,NZ)=SeedDepth_pft(NZ)
        Root1stChemElm(1:NumPlantChemElms,N,NR,NZ)=0._r8
      ENDDO D8790
    ENDDO D8795
    NIXBotRootLayer_pft(NZ)=NGTopRootLayer_pft(NZ)
    NumRootAxes_pft(NZ)=0
  ENDIF
  end associate
  end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadBranches(I,J,NZ,IDTHY,ShootNonstructC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: IDTHY
  real(r8), intent(inout) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: M,NE,NB
!     begin_execution
  associate(                                &
    iHarvstType_pft              =>  plt_distb%iHarvstType_pft     , &
    FWODBE                       =>  plt_allom%FWODBE    , &
    FWODLE                       =>  plt_allom%FWODLE    , &
    NoduleNonstructElmnt_brch    =>  plt_biom%NoduleNonstructElmnt_brch     , &
    StalkChemElms_brch         =>  plt_biom%StalkChemElms_brch    , &
    HuskChemElms_brch          =>  plt_biom%HuskChemElms_brch    , &
    PetoleChemElm_brch         =>  plt_biom%PetoleChemElm_brch   , &
    CanopyNoduleChemElm_brch   =>  plt_biom%CanopyNoduleChemElm_brch     , &
    LeafChemElms_brch   =>  plt_biom%LeafChemElms_brch    , &
    ZEROP     =>  plt_biom%ZEROP      , &
    EarChemElms_brch   =>  plt_biom%EarChemElms_brch    , &
    NonstructElmnt_brch    =>  plt_biom%NonstructElmnt_brch     , &
    GrainChemElms_brch    =>  plt_biom%GrainChemElms_brch     , &
    ReserveElmnts_brch   =>  plt_biom%ReserveElmnts_brch    , &
    StandingDeadKCompChemElms_pft    =>  plt_biom%StandingDeadKCompChemElms_pft     , &
    NonstructalElmnts_pft     =>  plt_biom%NonstructalElmnts_pft      , &
    CFOPE     =>  plt_soilchem%CFOPE  , &
    iPlantBranchState_brch     =>  plt_pheno%iPlantBranchState_brch     , &
    MatureGroup_brch    =>  plt_pheno%MatureGroup_brch    , &
    LeafNumberAtFloralInit_brch    =>  plt_pheno%LeafNumberAtFloralInit_brch    , &
    KLeafNodeNumber    =>  plt_pheno%KLeafNodeNumber    , &
    TotalNodeNumNormByMatgrp_brch    =>  plt_pheno%TotalNodeNumNormByMatgrp_brch    , &
    TotReproNodeNumNormByMatrgrp_brch    =>  plt_pheno%TotReproNodeNumNormByMatrgrp_brch    , &
    Hours4Leafout_brch      =>  plt_pheno%Hours4Leafout_brch      , &
    Hours4LeafOff_brch      =>  plt_pheno%Hours4LeafOff_brch      , &
    Hours4LenthenPhotoPeriod_brch      =>  plt_pheno%Hours4LenthenPhotoPeriod_brch      , &
    Hours4ShortenPhotoPeriod_brch      =>  plt_pheno%Hours4ShortenPhotoPeriod_brch      , &
    HourCounter4LeafOut_brch      =>  plt_pheno%HourCounter4LeafOut_brch      , &
    HourFailGrainFill_brch      =>  plt_pheno%HourFailGrainFill_brch      , &
    doInitLeafOut_brch     =>  plt_pheno%doInitLeafOut_brch     , &
    doPlantLeafOut_brch     =>  plt_pheno%doPlantLeafOut_brch     , &
    Prep4Literfall_brch     =>  plt_pheno%Prep4Literfall_brch     , &
    Hours4LiterfalAftMature_brch     =>  plt_pheno%Hours4LiterfalAftMature_brch     , &
    MatureGroup_pft   =>  plt_pheno%MatureGroup_pft   , &
    doPlantLeaveOff_brch     =>  plt_pheno%doPlantLeaveOff_brch     , &
    iPlantCalendar_brch    =>  plt_pheno%iPlantCalendar_brch    , &
    iPlantTurnoverPattern_pft    =>  plt_pheno%iPlantTurnoverPattern_pft    , &
    iPlantMorphologyType_pft    =>  plt_pheno%iPlantMorphologyType_pft    , &
    iPlantPhenolType_pft    =>  plt_pheno%iPlantPhenolType_pft    , &
    iPlantPhenologyPattern_pft    =>  plt_pheno%iPlantPhenologyPattern_pft    , &
    LitrFallChemElm_pvr      =>  plt_bgcr%LitrFallChemElm_pvr       , &
    NumOfBranches_pft       =>  plt_morph%NumOfBranches_pft       , &
    NodeNumberToInitFloral_brch     =>  plt_morph%NodeNumberToInitFloral_brch     , &
    ShootNodeNumber_brch     =>  plt_morph%ShootNodeNumber_brch     , &
    BranchNumber_brch      =>  plt_morph%BranchNumber_brch      , &
    NodeNumberAtAnthesis_brch     =>  plt_morph%NodeNumberAtAnthesis_brch     , &
    XTLI      =>  plt_morph%XTLI      , &
    NumOfLeaves_brch     =>  plt_morph%NumOfLeaves_brch     , &
    KLeafNumber_brch    =>  plt_morph%KLeafNumber_brch    , &
    istalk    =>  pltpar%istalk       , &
    k_fine_litr=>   pltpar%k_fine_litr, &
    k_woody_litr=> pltpar%k_woody_litr, &
    inonstruct   =>   pltpar%inonstruct   , &
    icwood     =>   pltpar%icwood     , &
    inonfoliar   =>  pltpar%inonfoliar    , &
    ifoliar    => pltpar%ifoliar      , &
    RubiscoActivity_brch      =>  plt_photo%RubiscoActivity_brch      , &
    C4PhotosynDowreg_brch     =>  plt_photo%C4PhotosynDowreg_brch       &
  )
  D8845: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iDead)THEN
      MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
      ShootNodeNumber_brch(NB,NZ)=XTLI(NZ)
      NodeNumberToInitFloral_brch(NB,NZ)=ShootNodeNumber_brch(NB,NZ)
      NodeNumberAtAnthesis_brch(NB,NZ)=0._r8
      NumOfLeaves_brch(NB,NZ)=0._r8
      LeafNumberAtFloralInit_brch(NB,NZ)=0._r8
      KLeafNumber_brch(NB,NZ)=1
      KLeafNodeNumber(NB,NZ)=1
      TotalNodeNumNormByMatgrp_brch(NB,NZ)=0._r8
      TotReproNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
      Hours4Leafout_brch(NB,NZ)=0._r8
      Hours4LeafOff_brch(NB,NZ)=0._r8
      Hours4LenthenPhotoPeriod_brch(NB,NZ)=0._r8
      Hours4ShortenPhotoPeriod_brch(NB,NZ)=0._r8
      HourCounter4LeafOut_brch(NB,NZ)=0._r8
      HourFailGrainFill_brch(NB,NZ)=0._r8
      RubiscoActivity_brch(NB,NZ)=1.0_r8
      C4PhotosynDowreg_brch(NB,NZ)=1.0_r8
      doInitLeafOut_brch(NB,NZ)=0
      doPlantLeafOut_brch(NB,NZ)=iDisable
      doPlantLeaveOff_brch(NB,NZ)=iEnable
      Prep4Literfall_brch(NB,NZ)=ifalse
      Hours4LiterfalAftMature_brch(NB,NZ)=0
      BranchNumber_brch(NB,NZ)=0
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
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      D6405: DO M=1,jsken
        DO NE=1,NumPlantChemElms        
          LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,inonstruct,M,NZ)*NoduleNonstructElmnt_brch(NE,NB,NZ) &
            +CFOPE(NE,ifoliar,M,NZ)*(LeafChemElms_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
            +CanopyNoduleChemElm_brch(NE,NB,NZ)) &
            +CFOPE(NE,inonfoliar,M,NZ)*(PetoleChemElm_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
            +HuskChemElms_brch(NE,NB,NZ)+EarChemElms_brch(NE,NB,NZ))

          LitrFallChemElm_pvr(NE,M,k_woody_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_woody_litr,0,NZ) &
            +CFOPE(NE,icwood,M,NZ)*(LeafChemElms_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr) &
            +PetoleChemElm_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

          IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
            NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)+CFOPE(NE,inonfoliar,M,NZ)*GrainChemElms_brch(NE,NB,NZ)
          ELSE
            LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,inonfoliar,M,NZ)*GrainChemElms_brch(NE,NB,NZ)
          ENDIF
          IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantMorphologyType_pft(NZ).LE.1)THEN
            LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,istalk,M,NZ)*StalkChemElms_brch(NE,NB,NZ)
          ELSE
            StandingDeadKCompChemElms_pft(NE,M,NZ)=StandingDeadKCompChemElms_pft(NE,M,NZ)+CFOPE(NE,icwood,M,NZ)*StalkChemElms_brch(NE,NB,NZ)
          ENDIF
        ENDDO
      ENDDO D6405

!
!     RECOVER NON-STRUCTURAL C,N,P FROM BRANCH TO
!     SEASONAL STORAGE RESERVES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     ShootNonstructC_brch=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
      NonstructalElmnts_pft(ielmc,NZ)=NonstructalElmnts_pft(ielmc,NZ)+ShootNonstructC_brch(NB,NZ)
      DO NE=1,NumPlantChemElms
        NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)+NonstructElmnt_brch(NE,NB,NZ)
        IF(iHarvstType_pft(NZ).NE.4.AND.iHarvstType_pft(NZ).NE.6)THEN
          D6406: DO M=1,jsken
            LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitrFallChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,inonstruct,M,NZ)*ReserveElmnts_brch(NE,NB,NZ)
          ENDDO D6406
        ELSE
          NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)+ReserveElmnts_brch(NE,NB,NZ)
        ENDIF
      ENDDO
!
      call ResetDeadRootStates(NB,NZ,ShootNonstructC_brch)

      IDTHY=IDTHY+1
    ENDIF
  ENDDO D8845
  end associate
  end subroutine LiterfallFromDeadBranches

!------------------------------------------------------------------------------------------

  subroutine ResetBranchRootStates(NZ,ShootNonstructC_brch)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                           &
    NonstructElmnt_brch => plt_biom%NonstructElmnt_brch        , &
    NoduleNonstructElmnt_brch => plt_biom%NoduleNonstructElmnt_brch        , &
    ShootChemElm_brch=> plt_biom%ShootChemElm_brch       , &
    LeafChemElms_brch=> plt_biom%LeafChemElms_brch       , &
    ReserveElmnts_brch=> plt_biom%ReserveElmnts_brch       , &
    PetoleChemElm_brch=> plt_biom%PetoleChemElm_brch      , &
    HuskChemElms_brch=> plt_biom%HuskChemElms_brch       , &
    EarChemElms_brch=> plt_biom%EarChemElms_brch       , &
    GrainChemElms_brch => plt_biom%GrainChemElms_brch        , &
    StalkBiomassC_brch => plt_biom%StalkBiomassC_brch        , &
     RootMycoNonstructElmnt_vr => plt_biom%RootMycoNonstructElmnt_vr        , &
    BranchStalkChemElms_pft=> plt_biom%BranchStalkChemElms_pft       , &
    LeafPetolBiomassC_brch  => plt_biom%LeafPetolBiomassC_brch         , &
    NonstructalElmnts_pft  => plt_biom%NonstructalElmnts_pft         , &
    Root1stStructChemElm_pvr => plt_biom%Root1stStructChemElm_pvr        , &
    Root2ndStructChemElm_pvr => plt_biom%Root2ndStructChemElm_pvr        , &
    CanopyNoduleChemElm_brch => plt_biom%CanopyNoduleChemElm_brch        , &
    StalkChemElms_brch=> plt_biom%StalkChemElms_brch       , &
    Root1stChemElm => plt_biom%Root1stChemElm        , &
    MaxNumRootLays     => plt_site%MaxNumRootLays            , &
    NU     => plt_site%NU            , &
    iPlantState_pft  => plt_pheno%iPlantState_pft        , &
    SecndRootLen  => plt_morph%SecndRootLen        , &
    SecndRootXNum_rpvr  => plt_morph%SecndRootXNum_rpvr        , &
    PrimRootLen  => plt_morph%PrimRootLen        , &
    MY     => plt_morph%MY           , &
    NumOfBranches_pft    => plt_morph%NumOfBranches_pft          , &
    NumRootAxes_pft   => plt_morph%NumRootAxes_pft           &
  )
!     RESET BRANCH STATE VARIABLES
!
  DO NE=1,NumPlantChemElms
    DO NB=1,NumOfBranches_pft(NZ)
      NonstructElmnt_brch(NE,NB,NZ)=0._r8
      NoduleNonstructElmnt_brch(NE,NB,NZ)=0._r8
      ShootChemElm_brch(NE,NB,NZ)=0._r8
      LeafChemElms_brch(NE,NB,NZ)=0._r8
      PetoleChemElm_brch(NE,NB,NZ)=0._r8
      StalkChemElms_brch(NE,NB,NZ)=0._r8
      ReserveElmnts_brch(NE,NB,NZ)=0._r8
    ENDDO
  ENDDO
  D8835: DO NB=1,NumOfBranches_pft(NZ)
    ShootNonstructC_brch(NB,NZ)=0._r8
    StalkBiomassC_brch(NB,NZ)=0._r8
    CanopyNoduleChemElm_brch(1:NumPlantChemElms,NB,NZ)=0._r8
    HuskChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
    EarChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
    GrainChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
    LeafPetolBiomassC_brch(NB,NZ)=0._r8
    BranchStalkChemElms_pft(1:NumPlantChemElms,NB,NZ)=0._r8
  ENDDO D8835
!
!     RESET ROOT STATE VARIABLES
!
  D6416: DO L=NU,MaxNumRootLays
    DO  N=1,MY(NZ)
       RootMycoNonstructElmnt_vr(1:NumPlantChemElms,N,L,NZ)=0._r8
      DO  NR=1,NumRootAxes_pft(NZ)
        Root1stStructChemElm_pvr(1:NumPlantChemElms,N,L,NR,NZ)=0._r8
        Root2ndStructChemElm_pvr(1:NumPlantChemElms,N,L,NR,NZ)=0._r8
        Root1stChemElm(1:NumPlantChemElms,N,NR,NZ)=0._r8
        PrimRootLen(N,L,NR,NZ)=0._r8
        SecndRootLen(N,L,NR,NZ)=0._r8
        SecndRootXNum_rpvr(N,L,NR,NZ)=0._r8
      enddo
    enddo
  ENDDO D6416
  NonstructalElmnts_pft(1:NumPlantChemElms,NZ)=0._r8
  iPlantState_pft(NZ)=1
  end associate
  end subroutine ResetBranchRootStates

!------------------------------------------------------------------------------------------

  subroutine ResetDeadRootStates(NB,NZ,ShootNonstructC_brch)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,K,N
!     begin_execution
  associate(                          &
    NonstructElmnt_brch   => plt_biom%NonstructElmnt_brch     , &
    NoduleNonstructElmnt_brch   => plt_biom%NoduleNonstructElmnt_brch     , &
    ShootChemElm_brch  => plt_biom%ShootChemElm_brch    , &
    LeafChemElms_brch  => plt_biom%LeafChemElms_brch    , &
    CanopyNoduleChemElm_brch   => plt_biom%CanopyNoduleChemElm_brch     , &
    StalkChemElms_brch  => plt_biom%StalkChemElms_brch    , &
    LeafElmntNode_brch    => plt_biom%LeafElmntNode_brch      , &
    ReserveElmnts_brch  => plt_biom%ReserveElmnts_brch    , &
    BranchStalkChemElms_pft  => plt_biom%BranchStalkChemElms_pft    , &
    StalkBiomassC_brch   => plt_biom%StalkBiomassC_brch     , &
    GrainChemElms_brch   => plt_biom%GrainChemElms_brch     , &
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch      , &
    HuskChemElms_brch  => plt_biom%HuskChemElms_brch    , &
    EarChemElms_brch  => plt_biom%EarChemElms_brch    , &
    PetioleProteinCNode_brch   => plt_biom%PetioleProteinCNode_brch     , &
    LeafProteinCNode_brch     => plt_biom%LeafProteinCNode_brch       , &
    InternodeChemElm_brch   => plt_biom%InternodeChemElm_brch     , &
    PetioleElmntNode_brch   => plt_biom%PetioleElmntNode_brch     , &
    LeafChemElmByLayer_pft   => plt_biom%LeafChemElmByLayer_pft     , &
    CanopyLeafCpft_lyr    => plt_biom%CanopyLeafCpft_lyr      , &
    PetoleChemElm_brch => plt_biom%PetoleChemElm_brch   , &
    CPOOL3   => plt_photo%CPOOL3    , &
    CPOOL4   => plt_photo%CPOOL4    , &
    CMassCO2BundleSheath_node     => plt_photo%CMassCO2BundleSheath_node      , &
    CMassHCO3BundleSheath_node     => plt_photo%CMassHCO3BundleSheath_node      , &
    GrainSeedBiomCMean_brch    => plt_allom%GrainSeedBiomCMean_brch     , &
    PotentialSeedSites_brch    => plt_morph%PotentialSeedSites_brch     , &
    SeedNumberSet_brch    => plt_morph%SeedNumberSet_brch     , &
    LeafAreaLive_brch    => plt_morph%LeafAreaLive_brch     , &
    CanopyBranchStemApft_lyr    => plt_morph%CanopyBranchStemApft_lyr     , &
    NumOfBranches_pft      => plt_morph%NumOfBranches_pft       , &
    LeafAreaNode_brch    => plt_morph%LeafAreaNode_brch     , &
    PetioleLengthNode_brch    => plt_morph%PetioleLengthNode_brch     , &
    InternodeHeightDying_brch   => plt_morph%InternodeHeightDying_brch    , &
    LeafAreaZsec_brch     => plt_morph%LeafAreaZsec_brch      , &
    InternodeHeightLive_brch   => plt_morph%InternodeHeightLive_brch    , &
    CanopyLeafApft_lyr    => plt_morph%CanopyLeafApft_lyr     , &
    StemAreaZsec_brch    => plt_morph%StemAreaZsec_brch     , &
    CanopyLeafAreaByLayer_pft    => plt_morph%CanopyLeafAreaByLayer_pft       &
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
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     SeedNumberSet_brch=seed set number
!     PotentialSeedSites_brch=potential number of seed set sites
!     GrainSeedBiomCMean_brch=individual seed size
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     LeafProteinCNode_brch=leaf protein mass
!     LeafAreaLive_brch=branch leaf area
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
!     InternodeChemElm_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
  ShootNonstructC_brch(NB,NZ)=0._r8
  NonstructElmnt_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  NoduleNonstructElmnt_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  ShootChemElm_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  LeafChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  CanopyNoduleChemElm_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  PetoleChemElm_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  StalkChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  ReserveElmnts_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  HuskChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  EarChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  GrainChemElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
  StalkBiomassC_brch(NB,NZ)=0._r8
  LeafPetolBiomassC_brch(NB,NZ)=0._r8
  PotentialSeedSites_brch(NB,NZ)=0._r8
  SeedNumberSet_brch(NB,NZ)=0._r8
  GrainSeedBiomCMean_brch(NB,NZ)=0._r8
  LeafAreaLive_brch(NB,NZ)=0._r8
  BranchStalkChemElms_pft(1:NumPlantChemElms,NB,NZ)=0._r8

  D8855: DO K=0,MaxNodesPerBranch1
    IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ)=0._r8
      CPOOL4(K,NB,NZ)=0._r8
      CMassCO2BundleSheath_node(K,NB,NZ)=0._r8
      CMassHCO3BundleSheath_node(K,NB,NZ)=0._r8
    ENDIF
    LeafAreaNode_brch(K,NB,NZ)=0._r8
    InternodeHeightLive_brch(K,NB,NZ)=0._r8
    InternodeHeightDying_brch(K,NB,NZ)=0._r8
    PetioleLengthNode_brch(K,NB,NZ)=0._r8
    LeafProteinCNode_brch(K,NB,NZ)=0._r8
    PetioleProteinCNode_brch(K,NB,NZ)=0._r8
    LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
    PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
    InternodeChemElm_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
    D8865: DO L=1,NumOfCanopyLayers1
      CanopyLeafApft_lyr(L,NZ)=CanopyLeafApft_lyr(L,NZ)-CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
      CanopyLeafCpft_lyr(L,NZ)=CanopyLeafCpft_lyr(L,NZ)-LeafChemElmByLayer_pft(ielmc,L,K,NB,NZ)
      CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=0._r8
      LeafChemElmByLayer_pft(1:NumPlantChemElms,L,K,NB,NZ)=0._r8
      IF(K.NE.0)THEN
        D8860: DO N=1,NumOfLeafZenithSectors1
          LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
        ENDDO D8860
      ENDIF
    ENDDO D8865
  ENDDO D8855
  D8875: DO L=1,NumOfCanopyLayers1
    CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
    DO  N=1,NumOfLeafZenithSectors1
      StemAreaZsec_brch(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D8875
  end associate
  end subroutine ResetDeadRootStates


end module LitrFallMod
