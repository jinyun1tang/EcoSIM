module PlantDisturbsMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use minimathmod,        only: isclose, AZMAX1
  use EcoSIMCtrlDataType, only: DazCurrYear
  use PlantBalMod,        only: SumPlantBranchBiome
  use ElmIDMod
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  use PlantMathFuncMod
  use PlantDisturbByFireMod
  use PlantDisturbByTillageMod
  use PlantDisturbByGrazingMod
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section

! disturbance variables
  real(r8) :: NonstructElmntRemoval(NumPlantChemElms)
  real(r8) :: LeafElmntRemoval(NumPlantChemElms)
  real(r8) :: FineNonleafElmntRemoval(NumPlantChemElms)
  real(r8) :: WoodyElmntRemoval(NumPlantChemElms)
  real(r8) :: StandeadElmntRemoval(NumPlantChemElms)
  real(r8) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8) :: StandeadElmnt2Litr(NumPlantChemElms)
  real(r8) :: LeafElmntHarv2Litr(NumPlantChemElms)
  real(r8) :: PetioleElmntHarv2Litr(NumPlantChemElms)
  real(r8) :: WoodyElmntHarv2Litr(NumPlantChemElms)
  real(r8) :: StandeadElmntHarv2Litr(NumPlantChemElms)
  real(r8) :: GrainHarvst(NumPlantChemElms)

  public :: RemoveBiomassByDisturbance
  public :: RemoveBiomByMgmt
  public :: InitPlantDisturbance
  contains

  subroutine InitPlantDisturbance
  implicit none
  call InitPlantFireMod
  end subroutine InitPlantDisturbance

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByMgmt(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
!     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
!
  NonstructElmntRemoval(1:NumPlantChemElms)   = 0._r8
  LeafElmntRemoval(1:NumPlantChemElms)        = 0._r8
  FineNonleafElmntRemoval(1:NumPlantChemElms) = 0._r8
  WoodyElmntRemoval(1:NumPlantChemElms)       = 0._r8
  LeafElmnt2Litr(1:NumPlantChemElms)          = 0._r8
  FineNonleafElmnt2Litr(1:NumPlantChemElms)   = 0._r8
  WoodyElmnt2Litr(1:NumPlantChemElms)         = 0._r8
  StandeadElmnt2Litr(1:NumPlantChemElms)      = 0._r8
  LeafElmntHarv2Litr(1:NumPlantChemElms)      = 0._r8
  PetioleElmntHarv2Litr(1:NumPlantChemElms)   = 0._r8
  WoodyElmntHarv2Litr(1:NumPlantChemElms)     = 0._r8
  GrainHarvst(1:NumPlantChemElms)             = 0._r8

  call RemoveBiomByHarvest(I,J,NZ)
!
!     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
!
  call RemoveBiomByTillage(I,J,NZ)

  end subroutine RemoveBiomByMgmt
!------------------------------------------------------------------------------------------

  subroutine RemoveStandingDead(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ

  real(r8) :: FHVSE
  real(r8) :: FHVSH
  integer :: M,NE
  associate(                                                   &
    NU                     => plt_site%NU,                     &
    ZERO                   => plt_site%ZERO,                   &
    AREA3                  => plt_site%AREA3,                  &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft,         &
    SolarNoonHour_col      => plt_site%SolarNoonHour_col,      &
    FracBiomHarvsted       => plt_distb%FracBiomHarvsted,      &
    THIN_pft               => plt_distb%THIN_pft,              &
    iHarvstType_pft        => plt_distb%iHarvstType_pft,       &
    StandDeadKCompElms_pft => plt_biom%StandDeadKCompElms_pft  &
  )
  StandeadElmntRemoval(1:NumPlantChemElms)=0._r8
  StandeadElmntHarv2Litr(1:NumPlantChemElms)=0._r8

  IF(J.EQ.INT(SolarNoonHour_col) .AND. iHarvstType_pft(NZ).NE.iharvtyp_grazing &
    .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(isclose(THIN_pft(NZ),0._r8))THEN
      FHVSE=AZMAX1(1._r8-FracBiomHarvsted(1,4,NZ))
      FHVSH=FHVSE
    ELSE
      FHVSE=AZMAX1(1._r8-THIN_pft(NZ))
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
        FHVSH=AZMAX1(1._r8-FracBiomHarvsted(1,4,NZ)*THIN_pft(NZ))
      ELSE
        FHVSH=FHVSE
      ENDIF
    ENDIF
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    call RemoveStandDeadByGrazing(I,J,NZ,FHVSE,FHVSH)
  ELSE
    FHVSE=1.0_r8
    FHVSH=1.0_r8
  ENDIF

  D6475: DO M=1,jsken
    DO NE=1,NumPlantChemElms
      StandeadElmntRemoval(NE)        = StandeadElmntRemoval(NE)+(1._r8-FHVSH)*StandDeadKCompElms_pft(NE,M,NZ)
      StandeadElmntHarv2Litr(NE)      = StandeadElmntHarv2Litr(NE)+(FHVSH-FHVSE)*StandDeadKCompElms_pft(NE,M,NZ)
      StandDeadKCompElms_pft(NE,M,NZ) = FHVSE*StandDeadKCompElms_pft(NE,M,NZ)
    ENDDO
  ENDDO D6475
  end associate
  end subroutine RemoveStandingDead  

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomassByDisturbance(I,J,NZ)
  implicit none
  integer , intent(in) :: I,J,NZ

!     begin_execution
  associate(                                             &
    NU                  => plt_site%NU,                  &
    ZERO                => plt_site%ZERO,                &
    AREA3               => plt_site%AREA3,               &
    iHarvstType_pft     => plt_distb%iHarvstType_pft,    &
    PlantPopulation_pft => plt_site%PlantPopulation_pft, &
    ZERO4Uptk_pft       => plt_rbgc%ZERO4Uptk_pft,       &
    ZERO4Groth_pft      => plt_biom%ZERO4Groth_pft,      &
    ZERO4LeafVar_pft    => plt_biom%ZERO4LeafVar_pft,    &
    ShootC4NonstC_brch  => plt_biom%ShootC4NonstC_brch   &
  )

!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     THIN_pft=thinning:fraction of population removed
!     FrcLeafMassNotHarvst(ielmc)=fraction of standing dead mass not harvested
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!
  IF(iHarvstType_pft(NZ).GE.iharvtyp_none)THEN
  !
    CALL RemoveStandingDead(I,J,NZ)
!
    call PlantDisturbance(I,J,NZ)

    ZERO4Groth_pft(NZ)=ZERO*PlantPopulation_pft(NZ)
    ZERO4Uptk_pft(NZ)=ZERO*PlantPopulation_pft(NZ)/AREA3(NU)
    ZERO4LeafVar_pft(NZ)=ZERO*PlantPopulation_pft(NZ)*1.0E+06_r8
  ENDIF
   
  end associate
  end subroutine RemoveBiomassByDisturbance

!------------------------------------------------------------------------------------------
  subroutine PlantDisturbance(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  real(r8) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8) :: HarvestElmnt2Litr(NumPlantChemElms)

  NonstructElmntOffEcosystem(1:NumPlantChemElms) = 0._r8
  LeafElmntOffEcosystem(1:NumPlantChemElms)      = 0._r8
  FineNonleafElmOffEcosystem(1:NumPlantChemElms) = 0._r8
  WoodyElmntOffEcosystem(1:NumPlantChemElms)     = 0._r8
  StandeadElmntOffEcosystem(1:NumPlantChemElms)  = 0._r8
  NonstructElmnt2Litr(1:NumPlantChemElms)        = 0._r8

  call ApplyDisturbanceBiomRemoval(I,J,NZ,NonstructElmnt2Litr,NonstructElmntOffEcosystem,&
    LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call AbvgBiomRemovalByDisturb(I,J,NZ,NonstructElmnt2Litr,HarvestElmnt2Litr,TotalElmnt2Litr)
!
!     ABOVE-GROUND LitrFall FROM HARVESTING
!
  call LiterfallByDisturbance(I,J,NZ,HarvestElmnt2Litr,TotalElmnt2Litr,NonstructElmnt2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
  end subroutine PlantDisturbance
!------------------------------------------------------------------------------------------

  subroutine LiterfallByDisturbance(I,J,NZ,HarvestElmnt2Litr,TotalElmnt2Litr,NonstructElmnt2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)

  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: HarvestElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  integer :: M,NE
  real(r8) :: dWoody
!     begin_execution
  associate(                                                                   &
    ilignin                        => pltpar%ilignin,                          &
    inonstruct                     => pltpar%inonstruct,                       &
    k_fine_litr                    => pltpar%k_fine_litr,                      &
    k_woody_litr                   => pltpar%k_woody_litr,                     &
    ifoliar                        => pltpar%ifoliar,                          &
    inonfoliar                     => pltpar%inonfoliar,                       &
    istalk                         => pltpar%istalk,                           &
    iroot                          => pltpar%iroot,                            &
    icwood                         => pltpar%icwood,                           &
    StandDeadKCompElms_pft         => plt_biom%StandDeadKCompElms_pft,         &
    iHarvstType_pft                => plt_distb%iHarvstType_pft,               &
    FracRootStalkElmAlloc2Litr     => plt_allom%FracRootStalkElmAlloc2Litr,    &
    ElmAllocmat4Litr               => plt_soilchem%ElmAllocmat4Litr,           &
    LitrfalStrutElms_pvr           => plt_bgcr%LitrfalStrutElms_pvr,           &
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft,     &
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft, &
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft,     &
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft          &
  )
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    !not by fire
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      D6375: DO M=1,jsken
        DO NE=1,NumPlantChemElms        
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*AZMAX1(NonstructElmnt2Litr(NE)) &
            +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafElmnt2Litr(NE)+LeafElmntHarv2Litr(NE)) &
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(NE)+PetioleElmntHarv2Litr(NE))

          IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmnt2Litr(NE)&
              +StandeadElmntHarv2Litr(NE))
          ELSE
            StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE))

            dWoody=AZMAX1(WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE))
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*dWoody*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*dWoody*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
          ENDIF
        ENDDO
      ENDDO D6375

!
!     ABOVE-GROUND LitrFall FROM FIRE
!
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
    ELSE
      call AbvGrndLiterFallByFire(I,J,NZ,NonstructElmnt2Litr,StandeadElmntOffEcosystem, &
        FineNonleafElmOffEcosystem,LeafElmnt2Litr,LeafElmntOffEcosystem,NonstructElmntOffEcosystem,&
        WoodyElmntOffEcosystem,WoodyElmnt2Litr,StandeadElmnt2Litr,PetioleElmntHarv2Litr,&
        FineNonleafElmnt2Litr,LeafElmntHarv2Litr,StandeadElmntHarv2Litr,WoodyElmntHarv2Litr)

    ENDIF
    !by grazing or herbivory
  ELSE
!
!     ABOVE-GROUND LitrFall FROM GRAZING
!
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P LitrFall
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P LitrFall
!
    DO NE=1,NumPlantChemElms
      LitrfalStrutElms_CumYr_pft(NE,NZ)=LitrfalStrutElms_CumYr_pft(NE,NZ)+TotalElmnt2Litr(NE)+HarvestElmnt2Litr(NE)
      SurfLitrfalStrutElms_CumYr_pft(NE,NZ)=SurfLitrfalStrutElms_CumYr_pft(NE,NZ)+TotalElmnt2Litr(NE)+HarvestElmnt2Litr(NE)
    ENDDO
  ENDIF
  end associate
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------
  subroutine AbvgBiomRemovalByDisturb(I,J,NZ,NonstructElmnt2Litr,HarvestElmnt2Litr,TotalElmnt2Litr)

  implicit none
  integer , intent(in)  :: I,J,NZ
  real(r8), intent(in)  :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: HarvestElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8) :: TotalElmntRemoval(NumPlantChemElms)
  integer :: NE
!     begin_execution
  associate(                                                      &
    iHarvstType_pft         => plt_distb%iHarvstType_pft,         &
    jHarvst_pft             => plt_distb%jHarvst_pft,             &
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft,    &
    EcoHavstElmnt_CumYr_pft => plt_distb%EcoHavstElmnt_CumYr_pft, &
    EcoHavstElmnt_CumYr_col => plt_distb%EcoHavstElmnt_CumYr_col, &
    CO2NetFix_pft           => plt_bgcr%CO2NetFix_pft,            &
    iYearCurrent            => plt_site%iYearCurrent,             &
    Eco_NBP_CumYr_col       => plt_bgcr%Eco_NBP_CumYr_col         &
  )
!
!     TotalElmntRemoval=total C,N,P removed
!     TotalElmnt2Litr=total C,N,P to litter
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  DO NE=1,NumPlantChemElms
    TotalElmntRemoval(NE) = NonstructElmntRemoval(NE)+LeafElmntRemoval(NE)+FineNonleafElmntRemoval(NE)+&
      WoodyElmntRemoval(NE)+StandeadElmntRemoval(NE)
    TotalElmnt2Litr(NE) = NonstructElmnt2Litr(NE)+LeafElmnt2Litr(NE)+FineNonleafElmnt2Litr(NE)+&
      WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE)
    HarvestElmnt2Litr(NE) = LeafElmntHarv2Litr(NE)+PetioleElmntHarv2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE)
  ENDDO

  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      !is not terminate and reseed
      IF(jHarvst_pft(NZ).NE.jharvtyp_tmareseed)THEN
        DO NE=1,NumPlantChemElms
          EcoHavstElmnt_CumYr_pft(NE,NZ) = EcoHavstElmnt_CumYr_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
          EcoHavstElmnt_CumYr_col(NE)    = EcoHavstElmnt_CumYr_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
        ENDDO
        Eco_NBP_CumYr_col=Eco_NBP_CumYr_col+TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc)
        !terminate and reseed
      ELSE
        DO NE=1,NumPlantChemElms
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
        ENDDO
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
    ELSE
      call AbvgBiomRemovalByFire(I,J,NZ,TotalElmnt2Litr,TotalElmntRemoval)
    ENDIF

  ELSE
!
!     C,N,P REMOVED FROM GRAZING
!  
     CALL AbvgBiomRemovalByGrazing(I,J,NZ,TotalElmnt2Litr,TotalElmntRemoval)

  ENDIF
  end associate
  end subroutine AbvgBiomRemovalByDisturb

!------------------------------------------------------------------------------------------
  subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,NonstructElmnt2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmntOffEcosystem(NumPlantChemElms)

  real(r8) :: EHVST21,EHVST22,EHVST23,EHVST24

  integer  :: NE
!     begin_execution
  associate(                                        &
    NU               => plt_site%NU,                &
    AREA3            => plt_site%AREA3,             &
    FracBiomHarvsted => plt_distb%FracBiomHarvsted, &
    FERT             => plt_distb%FERT,             &
    IYTYP            => plt_distb%IYTYP,            &
    iHarvstType_pft  => plt_distb%iHarvstType_pft   &
  )
!     IF NO PLANT C,N,P REMOVED AT HARVEST (ALL RESIDUE RETURNED)
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  EHVST21=1._r8-FracBiomHarvsted(2,iplthvst_leaf,NZ)
  EHVST22=1._r8-FracBiomHarvsted(2,iplthvst_finenonleaf,NZ)
  EHVST23=1._r8-FracBiomHarvsted(2,iplthvst_woody,NZ)
  EHVST24=1._r8-FracBiomHarvsted(2,iplthvst_stdead,NZ)

  IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)   = NonstructElmntRemoval(NE)*EHVST21     !non-structural
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*EHVST21          !leaf
      FineNonleafElmnt2Litr(NE) = FineNonleafElmntRemoval(NE)*EHVST22   !fine, non-woody
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*EHVST23         !woody
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*EHVST24      !standing dead
    ENDDO
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grain)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)   = NonstructElmntRemoval(NE)
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)
      FineNonleafElmnt2Litr(NE) = FineNonleafElmntRemoval(NE)-GrainHarvst(NE)*FracBiomHarvsted(2,iplthvst_finenonleaf,NZ)
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)
    ENDDO
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_allabv)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)   = NonstructElmntRemoval(NE)*EHVST21
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*EHVST21
      FineNonleafElmnt2Litr(NE) = FineNonleafElmntRemoval(NE)*EHVST22
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*EHVST23
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*EHVST24
    ENDDO
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)   = NonstructElmntRemoval(NE)*EHVST21
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*EHVST21
      FineNonleafElmnt2Litr(NE) = FineNonleafElmntRemoval(NE)*EHVST22
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*EHVST23
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*EHVST24
    ENDDO
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN

    call ApplyBiomRemovalByGrazing(I,J,NZ,EHVST21,EHVST22,EHVST23,EHVST24,&
    NonstructElmntRemoval,LeafElmntRemoval,FineNonleafElmntRemoval,WoodyElmntRemoval,StandeadElmntRemoval,&
    NonstructElmnt2Litr,LeafElmnt2Litr,FineNonleafElmnt2Litr,WoodyElmnt2Litr,StandeadElmnt2Litr)
!
!     REMOVALS BY FIRE
!
!     EFIRE=combustion  of N,P relative to C
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FracBiomHarvsted(2,1,FracBiomHarvsted(2,2,FracBiomHarvsted(2,3,FracBiomHarvsted(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_fire)THEN

    call ApplyBiomRemovalByFire(I,J,NZ,EHVST21,EHVST22, EHVST23, EHVST24,&
      StandeadElmntRemoval,NonstructElmntRemoval,LeafElmntRemoval,WoodyElmntRemoval,&
      FineNonleafElmntRemoval,NonstructElmnt2Litr,NonstructElmntOffEcosystem,&
      LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,WoodyElmntOffEcosystem,&
      StandeadElmntOffEcosystem,LeafElmnt2Litr,FineNonleafElmnt2Litr,&
      WoodyElmnt2Litr,StandeadElmnt2Litr)

  ENDIF
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!--------------------------------------------------------------------------------
  subroutine RemoveBiomByHarvest(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,K,M,NR,N,NB,NBX,NE
  real(r8) :: FracLeftThin
  real(r8) :: XHVST1
  REAL(R8) :: LeafC_lbrch(NumOfCanopyLayers1,JP1,JP1)
  real(r8) :: ARLFY,ARLFR
  real(r8) :: APSILT
  real(r8) :: FHGT
  real(r8) :: FHVSH
  real(r8) :: HvstedLeafC,HvstedShethC,HvstedEarC,HvstedGrainC,WHVSCP
  real(r8) :: HvstedStalkC,HvstedRsrvC
  real(r8) :: WHVSHH
  real(r8) :: WHVSNP
  real(r8) :: dLitR

!     begin_execution
  associate(                                                            &
    FracCanopyHeightCut_pft    => plt_distb%FracCanopyHeightCut_pft,    &
    THIN_pft                   => plt_distb%THIN_pft,                   &
    iHarvstType_pft            => plt_distb%iHarvstType_pft,            &
    jHarvst_pft                => plt_distb%jHarvst_pft,                &
    PlantPopulation_pft        => plt_site%PlantPopulation_pft,         &
    PPI_pft                    => plt_site%PPI_pft,                     &
    PPX_pft                    => plt_site%PPX_pft,                     &
    NU                         => plt_site%NU,                          &
    MaxNumRootLays             => plt_site%MaxNumRootLays,              &
    SolarNoonHour_col          => plt_site%SolarNoonHour_col,           &
    ZEROS                      => plt_site%ZEROS,                       &
    AREA3                      => plt_site%AREA3,                       &
    NumRootAxes_pft            => plt_morph%NumRootAxes_pft,            &
    ShootC4NonstC_brch         => plt_biom%ShootC4NonstC_brch,          &
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft,       &
    RootMyco1stElm_raxs        => plt_biom%RootMyco1stElm_raxs,         &
    CanopyStalkC_pft           => plt_biom%CanopyStalkC_pft,            &
    CanopyLeafShethC_pft       => plt_biom%CanopyLeafShethC_pft,        &
    StalkBiomassC_brch         => plt_biom%StalkBiomassC_brch,          &
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch,         &
    ShootStrutElms_brch        => plt_biom%ShootStrutElms_brch,         &
    LeafPetolBiomassC_brch     => plt_biom%LeafPetolBiomassC_brch,      &
    LeafElmsByLayerNode_brch   => plt_biom%LeafElmsByLayerNode_brch,    &
    StalkStrutElms_pft         => plt_biom%StalkStrutElms_pft,          &
    FracRootStalkElmAlloc2Litr => plt_allom%FracRootStalkElmAlloc2Litr, &
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft,    &
    iPlantPhenolType_pft       => plt_pheno%iPlantPhenolType_pft,       &
    ElmAllocmat4Litr           => plt_soilchem%ElmAllocmat4Litr,        &
    inonstruct                 => pltpar%inonstruct,                    &
    ifoliar                    => pltpar%ifoliar,                       &
    istalk                     => pltpar%istalk,                        &
    inonfoliar                 => pltpar%inonfoliar,                    &
    k_fine_litr                => pltpar%k_fine_litr,                   &
    k_woody_litr               => pltpar%k_woody_litr,                  &
    LitrfalStrutElms_pvr       => plt_bgcr%LitrfalStrutElms_pvr,        &
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft,         &
    MY                         => plt_morph%MY,                         &
    CanopyLeafAareZ_col        => plt_morph%CanopyLeafAareZ_col,        &
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col,          &
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft,          &
    CanopyStemArea_pft         => plt_morph%CanopyStemArea_pft,         &
    LiveInterNodeHight_brch    => plt_morph%LiveInterNodeHight_brch,    &
    CanopyStalkArea_lbrch      => plt_morph%CanopyStalkArea_lbrch,      &
    ClumpFactor_pft            => plt_morph%ClumpFactor_pft,            &
    CanopyLeafArea_col         => plt_morph%CanopyLeafArea_col          &
  )
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!

  IF((iHarvstType_pft(NZ).GE.iharvtyp_none .AND. J.EQ.INT(SolarNoonHour_col) &
    .AND. iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
    .OR. (iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo))THEN
!
!     ACCUMULATE ALL HARVESTED MATERIAL ABOVE CUTTING HEIGHT
!     ACCOUNTING FOR HARVEST EFFICIENCY ENTERED IN 'READQ'
!
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,and reseed
!     PPX,PP=PFT population per m2,grid cell
!     THIN_pft=thinning:fraction of population removed
!     CF=clumping factor
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     CanopyLeafArea_col,CanopyLeafAareZ_col=leaf area of combined canopy, canopy layer
!     ARLFR,ARLFY=leaf area harvested,remaining
!     ZL=height to bottom of each canopy layer
!
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(jHarvst_pft(NZ).NE.jharvtyp_tmareseed)THEN                
        PPX_pft(NZ)             = PPX_pft(NZ)*(1._r8-THIN_pft(NZ))
        PlantPopulation_pft(NZ) = PlantPopulation_pft(NZ)*(1._r8-THIN_pft(NZ))
        !terminate and reseed        
      ELSE
!     PPI_pft(NZ)=AMAX1(1.0_r8,0.5_r8*(PPI_pft(NZ)+CanopySeedNum_pft(NZ)/AREA3(NU)))
        PPX_pft(NZ)             = PPI_pft(NZ)
        PlantPopulation_pft(NZ) = PPX_pft(NZ)*AREA3(NU)
      ENDIF
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
        ClumpFactor_pft(NZ)=ClumpFactor_pft(NZ)*FracCanopyHeightCut_pft(NZ)
      ENDIF
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv .AND. FracCanopyHeightCut_pft(NZ).LT.0.0_r8)THEN
        ARLFY=(1._r8-ABS(FracCanopyHeightCut_pft(NZ)))*CanopyLeafArea_col
        ARLFR=0._r8
        D9875: DO L=1,NumOfCanopyLayers1
          IF(CanopyHeightZ_col(L).GT.CanopyHeightZ_col(L-1) &
            .AND. CanopyLeafAareZ_col(L).GT.ZEROS .AND. ARLFR.LT.ARLFY)THEN
            IF(ARLFR+CanopyLeafAareZ_col(L).GT.ARLFY)THEN
              FracCanopyHeightCut_pft(NZ)=CanopyHeightZ_col(L-1)+((ARLFY-ARLFR)/CanopyLeafAareZ_col(L))&
                *(CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1))
            ENDIF
          ELSE
            FracCanopyHeightCut_pft(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+CanopyLeafAareZ_col(L)
        ENDDO D9875
      ENDIF
      HvstedLeafC  = 0._r8
      HvstedShethC = 0._r8
      HvstedEarC   = 0._r8
      HvstedGrainC = 0._r8
      WHVSCP       = 0._r8
      HvstedStalkC = 0._r8
      HvstedRsrvC  = 0._r8
      LeafC_lbrch  = 0._r8          !it is a filler
    ELSE
!
!     GRAZING REMOVAL
      call GrazingPlant(I,J,NZ,HvstedLeafC,HvstedShethC,HvstedEarC,HvstedGrainC,&
        WHVSCP,HvstedStalkC,HvstedRsrvC,WHVSHH,WHVSNP,LeafC_lbrch)

    ENDIF
!
!     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
    call HarvestCanopy(I,J,NZ,HvstedLeafC,LeafC_lbrch)

    CALL CutPlant(I,J,NZ,WHVSHH,WHVSCP,WHVSNP,HvstedShethC,HvstedGrainC,HvstedEarC,&
      HvstedRsrvC,HvstedStalkC)

    CanopyLeafShethC_pft(NZ)     = 0._r8
    StalkStrutElms_pft(ielmc,NZ) = 0._r8
    CanopyStalkC_pft(NZ)         = 0._r8
    CanopyStemArea_pft(NZ)       = 0._r8
    D9840: DO NB=1,NumOfBranches_pft(NZ)
      CanopyLeafShethC_pft(NZ)     = CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)
      StalkStrutElms_pft(ielmc,NZ) = StalkStrutElms_pft(ielmc,NZ)+StalkStrutElms_brch(ielmc,NB,NZ)
      CanopyStalkC_pft(NZ)         = CanopyStalkC_pft(NZ)+StalkBiomassC_brch(NB,NZ)
      D9830: DO L=1,NumOfCanopyLayers1
        CanopyStemArea_pft(NZ)=CanopyStemArea_pft(NZ)+CanopyStalkArea_lbrch(L,NB,NZ)
      ENDDO D9830
    ENDDO D9840
!
!     ROOT LitrFall FROM HARVESTING OR FIRE
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     THETW=soil water concentration
!     CORGC=SOC concentration
!     iSoilDisturbType_col=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     EFIRE=combustion  of N,P relative to C
!     FrcLeafMassNotHarvst(ielmc),FrcLeafMassNotHarvst(ielmn),FrcLeafMassNotHarvst(ielmp)=fraction of root layer C,N,P not removed by disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     CO2ByFire_CumYr_pft,CH4ByFire_CumYr_pft,O2ByFire_CumYr_pft,NH3byFire_CumYr_pft,N2ObyFire_CumYr_pft,PO4byFire_CumYr_pft=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_CumYr_col=total net biome productivity
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      FracLeftThin=1.0_r8-THIN_pft(NZ)

      DO NR=1,NumRootAxes_pft(NZ)
        DO N=1,MY(NZ)
          DO NE=1,NumPlantChemElms
            RootMyco1stElm_raxs(NE,N,NR,NZ)=RootMyco1stElm_raxs(NE,N,NR,NZ)*FracLeftThin
          ENDDO
        ENDDO    
      ENDDO

      D3985: DO N=1,MY(NZ)
        D3980: DO L=NU,MaxNumRootLays
          CALL RootMaterialRemovalL(I,J,N,L,NZ,FracLeftThin,XHVST1)
          call HarvstUpdateRootStateL(I,J,N,L,NZ,FracLeftThin,XHVST1)
        ENDDO D3980
      ENDDO D3985
!
!     STORAGE LitrFall AND STATE VARIABLES DURING HARVESTING
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P

      IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
        XHVST1=1._r8-FracLeftThin
        DO NE=1,NumPlantChemElms
          D3400: DO M=1,jsken
            dLitR=XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*AZMAX1(SeasonalNonstElms_pft(NE,NZ))
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=&
                LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
              +dLitR*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)=&
                LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
              +dLitR*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)

          ENDDO D3400
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)*FracLeftThin
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

!--------------------------------------------------------------------------------

  subroutine ResetCutBranch(I,J,NZ,NB)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ,NB
  logical :: nonevgreenchk
  integer :: NBX,M

  associate(                                                                         &
   iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft,              &
   HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch,              &
   iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch,               &
   MatureGroup_pft                   => plt_pheno%MatureGroup_pft,                   &
   MatureGroup_brch                  => plt_pheno%MatureGroup_brch,                  &
   NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch,         &
   ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch,                 &
   TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch,     &
   NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch,           &
   FracHour4LeafoffRemob             => plt_allom%FracHour4LeafoffRemob,             &
   TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch, &
   doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch,                &
   NumOfBranches_pft                 => plt_morph%NumOfBranches_pft,                 &
   HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch,            &
   MainBranchNum_pft                 => plt_morph%MainBranchNum_pft,                 &
   LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch,       &
   Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 &
  )  
!     ZC=canopy height
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     GROUP=node number required for floral initiation
!     NodeNum2InitFloral_brch=node number at floral initiation
!     NodeNumberAtAnthesis_brch=node number at flowering
!     VSTGX=leaf number on date of floral initiation
!     TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
!     TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
!     HourFailGrainFill_brch=number of hours with no grain fill
!     doInitLeafOut_brch=flag for initializing leafout
!

  nonevgreenchk=iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen &
    .AND. Hours4LeafOff_brch(NB,NZ).LE.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)

  IF(nonevgreenchk .OR. (iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen .AND. iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0))THEN

    MatureGroup_brch(NB,NZ)                      = MatureGroup_pft(NZ)
    NodeNum2InitFloral_brch(NB,NZ)               = ShootNodeNum_brch(NB,NZ)
    NodeNumberAtAnthesis_brch(NB,NZ)             = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)           = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)         = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)     = 0._r8
    HourFailGrainFill_brch(NB,NZ)                = 0._r8
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)    = I
    iPlantCalendar_brch(2:NumGrowthStages,NB,NZ) = 0
    doInitLeafOut_brch(NB,NZ)                    = iEnable

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      D3010: DO NBX=1,NumOfBranches_pft(NZ)
        IF(NBX.NE.MainBranchNum_pft(NZ))THEN
          MatureGroup_brch(NBX,NZ)                   = MatureGroup_pft(NZ)
          NodeNum2InitFloral_brch(NBX,NZ)            = ShootNodeNum_brch(NBX,NZ)
          NodeNumberAtAnthesis_brch(NBX,NZ)          = 0._r8
          LeafNumberAtFloralInit_brch(NBX,NZ)        = 0._r8
          TotalNodeNumNormByMatgrp_brch(NBX,NZ)      = 0._r8
          TotReproNodeNumNormByMatrgrp_brch(NBX,NZ)  = 0._r8
          HourFailGrainFill_brch(NBX,NZ)             = 0._r8
          iPlantCalendar_brch(ipltcal_Emerge,NBX,NZ) = I
          D3015: DO M=2,NumGrowthStages
            iPlantCalendar_brch(M,NBX,NZ)=0
          ENDDO D3015
          doInitLeafOut_brch(NBX,NZ)=0
        ENDIF
      ENDDO D3010
    ENDIF
  ENDIF
  end associate
  end subroutine ResetCutBranch

!--------------------------------------------------------------------------------

  subroutine BranchCutPlantStalk(I,J,NB,NZ,RMedInternodeLen,HvstedStalkC,HvstedRsrvC)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: RMedInternodeLen
  real(r8), intent(in) :: HvstedStalkC,HvstedRsrvC
  integer :: NE,K  
  real(r8) :: FrcLeafMassLeft,FHGT,FHGTK,FHVSETS,FHVSH

  associate(                                                          &
    FracCanopyHeightCut_pft   => plt_distb%FracCanopyHeightCut_pft,   &
    THIN_pft                  => plt_distb%THIN_pft,                  &
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch,         &
    InternodeHeightDying_brch => plt_morph%InternodeHeightDying_brch, &
    FracBiomHarvsted          => plt_distb%FracBiomHarvsted,          &
    StalkBiomassC_brch        => plt_biom%StalkBiomassC_brch,         &
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch,        &
    LiveInterNodeHight_brch   => plt_morph%LiveInterNodeHight_brch,   &
    StalkRsrvElms_pft         => plt_biom%StalkRsrvElms_pft,          &
    SenecStalkStrutElms_brch  => plt_biom%SenecStalkStrutElms_brch,   &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,             &
    InternodeStrutElms_brch   => plt_biom%InternodeStrutElms_brch,    &
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft,           &
    ZERO                      => plt_site%ZERO,                       &
    StalkStrutElms_pft        => plt_biom%StalkStrutElms_pft,         &
    iHarvstType_pft           => plt_distb%iHarvstType_pft            &
  )
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     RMedInternodeLen=internode length removed
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     FHGT=fraction of canopy layer height not harvested
!     FrcLeafMassLeft=fraction of canopy layer mass not harvested
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WTSTK=stalk C mass
!
!
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(RMedInternodeLen.GT.ZERO)THEN
      IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
        FHGT=AZMAX1(AMIN1(1.0_r8,FracCanopyHeightCut_pft(NZ)/RMedInternodeLen))
      ELSE
        FHGT=0._r8
      ENDIF
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FrcLeafMassLeft=AZMAX1(1._r8-(1._r8-FHGT)*FracBiomHarvsted(1,iplthvst_woody,NZ))
        FHVSH=FrcLeafMassLeft
      ELSE
        FrcLeafMassLeft=AZMAX1(1._r8-THIN_pft(NZ))
        IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
          FHVSH=1.0_r8-(1._r8-FHGT)*FracBiomHarvsted(1,iplthvst_woody,NZ)*THIN_pft(NZ)
        ELSE
          FHVSH=FrcLeafMassLeft
        ENDIF
      ENDIF
    ELSE
      FrcLeafMassLeft=1.0_r8
      FHVSH=1.0_r8
    ENDIF
  ELSE
    !grazing or herbivory
    IF(StalkStrutElms_pft(ielmc,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
      FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,1._r8-HvstedStalkC/StalkStrutElms_pft(ielmc,NZ)))
      FHVSH=FrcLeafMassLeft
    ELSE
      FrcLeafMassLeft=1.0_r8
      FHVSH=1.0_r8
    ENDIF
  ENDIF
    !
!     HARVESTED STALK C,N,P
!
!     WoodyElmntRemoval=harvested stalk C,N,P
!     WoodyElmntHarv2Litr=harvested stalk C,N,P to litter
!     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in harvested stalk
!
  DO NE=1,NumPlantChemElms
    WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE)+(1._r8-FHVSH)*StalkStrutElms_brch(NE,NB,NZ)
    WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE)+(FHVSH-FrcLeafMassLeft)*StalkStrutElms_brch(NE,NB,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
    StalkStrutElms_brch(NE,NB,NZ)=FrcLeafMassLeft*StalkStrutElms_brch(NE,NB,NZ)
    SenecStalkStrutElms_brch(NE,NB,NZ)=FrcLeafMassLeft*SenecStalkStrutElms_brch(NE,NB,NZ)
  ENDDO

  StalkBiomassC_brch(NB,NZ)=FrcLeafMassLeft*StalkBiomassC_brch(NB,NZ)
!
!     CUT STALK NODES
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     InternodeHeightDying_brch,LiveInterNodeHight_brch=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTSTK=stalk C mass
!     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
  D9820: DO K=MaxNodesPerBranch1,0,-1
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(InternodeHeightDying_brch(K,NB,NZ).GT.ZERO)THEN
        IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
          FHGTK=AZMAX1(AMIN1(1.0_r8,(LiveInterNodeHight_brch(K,NB,NZ)-FracCanopyHeightCut_pft(NZ))/&
            InternodeHeightDying_brch(K,NB,NZ)))
        ELSE
          FHGTK=0._r8
        ENDIF
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FHVSETS=AZMAX1(1._r8-FHGTK*FracBiomHarvsted(1,iplthvst_woody,NZ))
        ELSE
          FHVSETS=AZMAX1(1._r8-THIN_pft(NZ))
        ENDIF
      ELSE
        FHVSETS=1.0_r8
      ENDIF
    ELSE
      IF(StalkStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FHVSETS=AZMAX1(AMIN1(1.0_r8,1._r8-HvstedStalkC/StalkStrutElms_pft(ielmc,NZ)))
      ELSE
        FHVSETS=1.0_r8
      ENDIF
    ENDIF
    DO NE=1,NumPlantChemElms
      InternodeStrutElms_brch(NE,K,NB,NZ)=FHVSETS*InternodeStrutElms_brch(NE,K,NB,NZ)
    ENDDO
    IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.isclose(THIN_pft(NZ),0._r8))THEN
      InternodeHeightDying_brch(K,NB,NZ)=FHVSETS*InternodeHeightDying_brch(K,NB,NZ)
      LiveInterNodeHight_brch(K,NB,NZ)=AMIN1(LiveInterNodeHight_brch(K,NB,NZ),FracCanopyHeightCut_pft(NZ))
    ENDIF

  ENDDO D9820
!
!     CUT STALK RESERVES
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSTKB=C mass remaining in harvested stalk
!     WTRSV=stalk reserve C mass
!     HvstedRsrvC=remaining stalk reserve C mass
!     FrcLeafMassLeft=fraction of reserve mass not harvested
! double check below
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(StalkStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=FrcLeafMassLeft
      FHVSH=FHVSH
    ELSE
      FrcLeafMassLeft=0._r8
      FHVSH=0._r8
    ENDIF
    !grazing or herbivory
  ELSE
    !
    IF(StalkRsrvElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,1._r8-HvstedRsrvC/StalkRsrvElms_pft(ielmc,NZ)))
      FHVSH=FrcLeafMassLeft
    ELSE
      FrcLeafMassLeft=0._r8
      FHVSH=0._r8
    ENDIF
  ENDIF
!
!     HARVESTED STALK RESERVE C,N,P
!
!     WoodyElmntRemoval=harvested stalk C,N,P
!     WoodyElmntHarv2Litr=harvested stalk C,N,P to litter
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!
  DO NE=1,NumPlantChemElms
    WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE)+(1._r8-FHVSH)*StalkRsrvElms_brch(NE,NB,NZ)
    WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(ielmc)+(FHVSH-FrcLeafMassLeft)*StalkRsrvElms_brch(NE,NB,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
    StalkRsrvElms_brch(NE,NB,NZ)=FrcLeafMassLeft*StalkRsrvElms_brch(NE,NB,NZ)
  ENDDO
  end associate
  end subroutine BranchCutPlantStalk

!--------------------------------------------------------------------------------
  subroutine BranchCutReprodOrgans(I,J,NB,NZ,RMedInternodeLen,HvstedShethC,HvstedGrainC,HvstedEarC)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: RMedInternodeLen
  real(r8), intent(in) :: HvstedShethC
  real(r8), intent(in) :: HvstedGrainC
  real(r8), intent(in) :: HvstedEarC

  integer :: NE
  real(r8) :: FracGrainNotHvsted,FracShethGrainNotHvsted
  real(r8) :: FracHuskNotHvsted,FracEarNotHvsted,FracShethHuskNotHvsted,FracShethNotHvsted
  associate(                                                      &
    FracCanopyHeightCut_pft => plt_distb%FracCanopyHeightCut_pft, &
    GrainStrutElms_brch     => plt_biom%GrainStrutElms_brch,      &
    EarStrutElms_brch       => plt_biom%EarStrutElms_brch,        &
    PotentialSeedSites_brch => plt_morph%PotentialSeedSites_brch, &
    FracBiomHarvsted        => plt_distb%FracBiomHarvsted,        &
    SeedNumSet_brch         => plt_morph%SeedNumSet_brch,         &
    GrainStrutElms_pft      => plt_biom%GrainStrutElms_pft,       &
    GrainSeedBiomCMean_brch => plt_allom%GrainSeedBiomCMean_brch, &
    THIN_pft                => plt_distb%THIN_pft,                &
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft,           &
    HuskStrutElms_pft       => plt_biom%HuskStrutElms_pft,        &
    HuskStrutElms_brch      => plt_biom%HuskStrutElms_brch,       &
    EarStrutElms_pft        => plt_biom%EarStrutElms_pft,         &
    iHarvstType_pft         => plt_distb%iHarvstType_pft          &
  )
! 
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FracGrainNotHvsted,FracHuskNotHvsted,FracEarNotHvsted=fraction of grain,husk,ear mass not harvested
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(FracCanopyHeightCut_pft(NZ).LT.RMedInternodeLen .OR. iHarvstType_pft(NZ).EQ.iharvtyp_grain &
      .OR. iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FracGrainNotHvsted=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)
        FracShethGrainNotHvsted=FracGrainNotHvsted
      ELSE
        FracGrainNotHvsted=1.0_r8-THIN_pft(NZ)
        FracShethGrainNotHvsted=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
      ENDIF
    ELSE
      FracGrainNotHvsted=1.0_r8-THIN_pft(NZ)
      FracShethGrainNotHvsted=FracGrainNotHvsted
    ENDIF
    FracHuskNotHvsted      = FracGrainNotHvsted
    FracEarNotHvsted       = FracGrainNotHvsted
    FracShethHuskNotHvsted = FracShethGrainNotHvsted
    FracShethNotHvsted     = FracShethGrainNotHvsted

    !grazing or herbivory
  ELSE

    IF(HuskStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracHuskNotHvsted      = AZMAX1(AMIN1(1.0_r8,1._r8-HvstedShethC/HuskStrutElms_pft(ielmc,NZ)))
      FracShethHuskNotHvsted = FracHuskNotHvsted
    ELSE
      FracHuskNotHvsted      = 1.0_r8
      FracShethHuskNotHvsted = 1.0_r8
    ENDIF

    IF(EarStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracEarNotHvsted   = AZMAX1(AMIN1(1.0_r8,1._r8-HvstedEarC/EarStrutElms_pft(ielmc,NZ)))
      FracShethNotHvsted = FracEarNotHvsted
    ELSE
      FracEarNotHvsted   = 1.0_r8
      FracShethNotHvsted = 1.0_r8
    ENDIF

    IF(GrainStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracGrainNotHvsted      = AZMAX1(AMIN1(1.0_r8,1._r8-HvstedGrainC/GrainStrutElms_pft(ielmc,NZ)))
      FracShethGrainNotHvsted = FracGrainNotHvsted
    ELSE
      FracGrainNotHvsted      = 1.0_r8
      FracShethGrainNotHvsted = 1.0_r8
    ENDIF
  ENDIF
!
!     HARVESTED REPRODUCTIVE C,N,P
!
!     FineNonleafElmntRemoval=reproductive C,N,P removed
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     GrainHarvst(ielmc),GrainHarvst(ielmn),GrainHarvst(ielmp)=grain harvested
!
  DO NE=1,NumPlantChemElms
    FineNonleafElmntRemoval(NE)=FineNonleafElmntRemoval(NE)+(1._r8-FracShethHuskNotHvsted)*HuskStrutElms_brch(NE,NB,NZ)&
      +(1._r8-FracShethNotHvsted)*EarStrutElms_brch(NE,NB,NZ)+(1._r8-FracShethGrainNotHvsted)*GrainStrutElms_brch(NE,NB,NZ)
    PetioleElmntHarv2Litr(NE)=PetioleElmntHarv2Litr(NE)+(FracShethHuskNotHvsted-FracHuskNotHvsted)*HuskStrutElms_brch(NE,NB,NZ) &
      +(FracShethNotHvsted-FracEarNotHvsted)*EarStrutElms_brch(NE,NB,NZ) &
      +(FracShethGrainNotHvsted-FracGrainNotHvsted)*GrainStrutElms_brch(NE,NB,NZ)
    GrainHarvst(NE)=GrainHarvst(NE)+(1._r8-FracGrainNotHvsted)*GrainStrutElms_brch(NE,NB,NZ)

!
!     REMAINING REPRODUCTIVE C,N,P
!
    HuskStrutElms_brch(NE,NB,NZ)  = FracHuskNotHvsted*HuskStrutElms_brch(NE,NB,NZ)
    EarStrutElms_brch(NE,NB,NZ)   = FracEarNotHvsted*EarStrutElms_brch(NE,NB,NZ)
    GrainStrutElms_brch(NE,NB,NZ) = FracGrainNotHvsted*GrainStrutElms_brch(NE,NB,NZ)
  ENDDO
  PotentialSeedSites_brch(NB,NZ) = FracGrainNotHvsted*PotentialSeedSites_brch(NB,NZ)
  SeedNumSet_brch(NB,NZ)         = FracGrainNotHvsted*SeedNumSet_brch(NB,NZ)
  GrainSeedBiomCMean_brch(NB,NZ) = FracGrainNotHvsted*GrainSeedBiomCMean_brch(NB,NZ)
  end associate
  END subroutine BranchCutReprodOrgans

!--------------------------------------------------------------------------------

  subroutine CutBranchNonstructural(I,J,NB,NZ,LeafCafCut_brch,WGSHGX,LeafCbfCut_brch,WGSHGY,WHVSCP,WHVSNP)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: LeafCafCut_brch,WGSHGX,LeafCbfCut_brch,WGSHGY
  real(r8), intent(in) :: WHVSNP,WHVSCP
  real(r8) :: CanopyNonstElmCopy(NumPlantChemElms)  
  real(r8) :: CanopyNodulNonstElmCopy(NumPlantChemElms)
  real(r8) :: CanopyNonstElmAfhvst(NumPlantChemElms)
  real(r8) :: CanopyNodulNonstElmAfhvst(NumPlantChemElms)  
  real(r8) :: CanopyNodulStrutElmAfhvst(NumPlantChemElms)
  real(r8) :: FHVST4,dFHVST4  
  integer  :: K,NE
  real(r8) :: FrcLeafMassLeft
  
  associate(                                                            &
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch,        &
    CanopyNodulStrutElms_brch  => plt_biom%CanopyNodulStrutElms_brch,   &
    CPOOL3_node                => plt_photo%CPOOL3_node,                &
    CPOOL4_node                => plt_photo%CPOOL4_node,                &
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node, &
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node,  &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft,              &
    CanopyNodulNonstElms_brch  => plt_biom%CanopyNodulNonstElms_brch,   &
    iPlantPhotosynthesisType   => plt_photo%iPlantPhotosynthesisType,   &
    iHarvstType_pft            => plt_distb%iHarvstType_pft             &
  )
!
!     CUT NON-STRUCTURAL C,N,P IN HARVESTED BRANCHES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     FrcLeafMassLeft=fraction of leaf+petiole node mass not harvested
!     CanopyNonstElmAfhvst(),=branch non-structural C,N,P mass after harvest
!     CanopyNodulNonstElmAfhvst=nonstructural C,N,P in bacteria after harvest
!     CanopyNodulStrutElmAfhvst=bacterial C,N,P mass after harvest
!     WTLS,LeafPetolBiomassC_brch=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
  DO NE=1,NumPlantChemElms
    CanopyNonstElmCopy(NE)      = AZMAX1(CanopyNonstElms_brch(NE,NB,NZ))
    CanopyNodulNonstElmCopy(NE) = AZMAX1(CanopyNodulNonstElms_brch(NE,NB,NZ))
  ENDDO

  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(LeafCbfCut_brch+WGSHGY.GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,(LeafCafCut_brch+WGSHGX)/(LeafCbfCut_brch+WGSHGY)))
      DO NE=1,NumPlantChemElms
        CanopyNonstElmAfhvst(NE)      = CanopyNonstElmCopy(NE)*FrcLeafMassLeft
        CanopyNodulNonstElmAfhvst(NE) = CanopyNodulNonstElmCopy(NE)*FrcLeafMassLeft
        CanopyNodulStrutElmAfhvst(NE) = CanopyNodulStrutElms_brch(NE,NB,NZ)*FrcLeafMassLeft
      ENDDO
    ELSE
      CanopyNonstElmAfhvst(:)      = 0._r8
      CanopyNodulNonstElmAfhvst(:) = 0._r8
      CanopyNodulStrutElmAfhvst(:) = 0._r8
    ENDIF
    !grazing
  ELSE
    call CutBranchNonstalByGrazing(I,J,NB,NZ,WHVSCP,CanopyNonstElmCopy,&
      WHVSNP,CanopyNodulNonstElmCopy,CanopyNonstElmAfhvst,CanopyNodulNonstElmAfhvst,&
      CanopyNodulStrutElmAfhvst)

  ENDIF
!
!     HARVESTED NON-STRUCTURAL C, N, P
!
!     NonstructElmntRemoval(ielmc),NonstructElmntRemoval(ielmn),NonstructElmntRemoval(ielmp)=nonstructural C,N,P removed
!
  DO NE=1,NumPlantChemElms
    NonstructElmntRemoval(NE)=NonstructElmntRemoval(NE)+CanopyNonstElmCopy(NE)-CanopyNonstElmAfhvst(NE)+&
      CanopyNodulNonstElmCopy(NE)-CanopyNodulNonstElmAfhvst(NE)
    NonstructElmntRemoval(NE)=NonstructElmntRemoval(NE)+CanopyNodulStrutElms_brch(NE,NB,NZ)-CanopyNodulStrutElmAfhvst(NE)
  ENDDO
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  CanopyNonstElms_brch(:,NB,NZ)=CanopyNonstElmAfhvst(:)
  CanopyNodulNonstElms_brch(:,NB,NZ)=CanopyNodulNonstElmAfhvst(:)
  CanopyNodulStrutElms_brch(:,NB,NZ)=CanopyNodulStrutElmAfhvst(:)
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CanopyNonstElmAfhvst(ielmc)=branch non-structural C mass after harvest
!     NonstructElmntRemoval()=nonstructural C,N,P removed
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!
  IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo .AND. CanopyNonstElmCopy(ielmc).GT.ZERO4Groth_pft(NZ))THEN
    FHVST4=CanopyNonstElmAfhvst(ielmc)/CanopyNonstElmCopy(ielmc)
    dFHVST4=1._r8-FHVST4
    D9810: DO K=1,MaxNodesPerBranch1
      NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc) &
        +dFHVST4*(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ))
      CPOOL3_node(K,NB,NZ)                = FHVST4*CPOOL3_node(K,NB,NZ)
      CPOOL4_node(K,NB,NZ)                = FHVST4*CPOOL4_node(K,NB,NZ)
      CMassCO2BundleSheath_node(K,NB,NZ)  = FHVST4*CMassCO2BundleSheath_node(K,NB,NZ)
      CMassHCO3BundleSheath_node(K,NB,NZ) = FHVST4*CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D9810
  ENDIF
  end associate
  end subroutine CutBranchNonstructural
!--------------------------------------------------------------------------------
  subroutine CutBranchSheathPetole(I,J,NB,NZ,WHVSHH,FracIntnodeNotHvsted,&
    FracNodeNotHvsted,RMedInternodeLen,WGSHGX,WGSHGY)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: WHVSHH
  real(r8), intent(inout) :: FracIntnodeNotHvsted(0:MaxNodesPerBranch1)
  real(r8), intent(inout) :: FracNodeNotHvsted(0:MaxNodesPerBranch1)  
  real(r8), intent(out) :: RMedInternodeLen
  real(r8), intent(out) :: WGSHGX,WGSHGY
  real(r8) :: PetoleCRmved_brnch,FHGT
  integer :: K,NE

  associate(                                                            &
    PetoleStrutElms_pft        => plt_biom%PetoleStrutElms_pft,         &
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch,        &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft,              &
    PetioleElmntNode_brch      => plt_biom%PetioleElmntNode_brch,       &
    FracCanopyHeightCut_pft    => plt_distb%FracCanopyHeightCut_pft,    &
    LiveInterNodeHight_brch    => plt_morph%LiveInterNodeHight_brch,    &
    PetoleLensNode_brch        => plt_morph%PetoleLensNode_brch,        &
    PetoleProteinCNode_brch    => plt_biom%PetoleProteinCNode_brch,     &
    k_fine_litr                => pltpar%k_fine_litr,                   &
    k_woody_litr               => pltpar%k_woody_litr,                  &
    FracShootLeafElmAlloc2Litr => plt_allom%FracShootLeafElmAlloc2Litr, &
    iHarvstType_pft            => plt_distb%iHarvstType_pft             &
  )
  WGSHGY=0._r8
  WGSHGX=0._r8
!
!     CUT SHEATHS OR PETIOLES AND STALKS HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSHE,WTSHEB=PFT,branch petiole C mass
!     PetoleCRmved_brnch,WHVSHH=branch, PFT petiole C mass removed
!     LiveInterNodeHight_brch=internode length
!     RMedInternodeLen=internode length removed
!
  RMedInternodeLen=0._r8
  IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
    .AND. PetoleStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    PetoleCRmved_brnch=WHVSHH*PetoleStrutElms_brch(ielmc,NB,NZ)/PetoleStrutElms_pft(ielmc,NZ)
  ELSE
    PetoleCRmved_brnch=0._r8
  ENDIF

  D9805: DO K=MaxNodesPerBranch1,0,-1
    IF(LiveInterNodeHight_brch(K,NB,NZ).GT.0.0) RMedInternodeLen=AMAX1(RMedInternodeLen,LiveInterNodeHight_brch(K,NB,NZ))
!
!     HARVESTED SHEATH OR PETIOLE C,N,P
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     PetoleCRmved_brnch=branch petiole C mass removed
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinCNode_brch=node petiole C,N,P,protein mass
!     FracIntnodeNotHvsted=fraction of internode layer mass not harvested
!     FineNonleafElmntRemoval=harvested petiole C,N,P
!     PetioleElmntHarv2Litr=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     PetoleLensNode_brch,LiveInterNodeHight_brch=petiole,internode length
!
    IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
      .OR. PetoleCRmved_brnch.GT.0.0_r8)THEN
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
        IF(PetioleElmntNode_brch(ielmc,K,NB,NZ).GT.PetoleCRmved_brnch)THEN
          FracIntnodeNotHvsted(K)=AZMAX1(AMIN1(1.0_r8,(PetioleElmntNode_brch(ielmc,K,NB,NZ)-PetoleCRmved_brnch) &
            /PetioleElmntNode_brch(ielmc,K,NB,NZ)))
          FracNodeNotHvsted(K)=FracIntnodeNotHvsted(K)
        ELSE
          FracIntnodeNotHvsted(K) = 0._r8
          FracNodeNotHvsted(K)    = 0._r8
        ENDIF
      ENDIF

      PetoleCRmved_brnch=PetoleCRmved_brnch-(1._r8-FracIntnodeNotHvsted(K))*PetioleElmntNode_brch(ielmc,K,NB,NZ)
      DO NE=1,NumPlantChemElms
        FineNonleafElmntRemoval(NE)=FineNonleafElmntRemoval(NE) &
          +(1._r8-FracNodeNotHvsted(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        PetioleElmntHarv2Litr(NE)=PetioleElmntHarv2Litr(NE) &
          +(FracNodeNotHvsted(K)-FracIntnodeNotHvsted(K))*PetioleElmntNode_brch(NE,K,NB,NZ)&
          *FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE) &
          +(1._r8-FracNodeNotHvsted(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
        WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE) &
          +(FracNodeNotHvsted(K)-FracIntnodeNotHvsted(K))*PetioleElmntNode_brch(NE,K,NB,NZ)&
          *FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
      ENDDO
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     PetioleElmntNode_brch=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     PetoleLensNode_brch=node petiole height
!     PetoleProteinCNode_brch=petiole protein mass
!
      WGSHGY=WGSHGY+PetioleElmntNode_brch(ielmc,K,NB,NZ)
      PetoleProteinCNode_brch(K,NB,NZ)=FracIntnodeNotHvsted(K)*PetoleProteinCNode_brch(K,NB,NZ)

      DO NE=1,NumPlantChemElms
        PetoleStrutElms_brch(NE,NB,NZ)=PetoleStrutElms_brch(NE,NB,NZ) &
          -(1._r8-FracIntnodeNotHvsted(K))*PetioleElmntNode_brch(NE,K,NB,NZ)
        PetioleElmntNode_brch(NE,K,NB,NZ)=FracIntnodeNotHvsted(K)*PetioleElmntNode_brch(NE,K,NB,NZ)
      ENDDO
!            PetoleProteinCNode_brch(K,NB,NZ)=FracIntnodeNotHvsted(K)*PetoleProteinCNode_brch(K,NB,NZ)
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.PetoleLensNode_brch(K,NB,NZ).GT.0.0_r8)THEN
        FHGT=AZMAX1(AMIN1(1.0_r8,(LiveInterNodeHight_brch(K,NB,NZ) &
          +PetoleLensNode_brch(K,NB,NZ)-FracCanopyHeightCut_pft(NZ))/PetoleLensNode_brch(K,NB,NZ)))
        PetoleLensNode_brch(K,NB,NZ)=(1._r8-FHGT)*PetoleLensNode_brch(K,NB,NZ)
      ELSE
        PetoleLensNode_brch(K,NB,NZ)=FracIntnodeNotHvsted(K)*PetoleLensNode_brch(K,NB,NZ)
      ENDIF
      WGSHGX=WGSHGX+PetioleElmntNode_brch(ielmc,K,NB,NZ)

!     IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
!     IF(LiveInterNodeHight_brch(K,NB,NZ).GT.FracCanopyHeightCut_pft(NZ)
!    2.OR.iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
!     IF(isclose(FracIntnodeNotHvsted(K),0._r8).AND.K.GT.0)THEN
!     IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
!     NumOfLeaves_brch(NB,NZ)=AZMAX1(NumOfLeaves_brch(NB,NZ)-1.0)
!     ELSE
!     NumOfLeaves_brch(NB,NZ)=AZMAX1(NumOfLeaves_brch(NB,NZ)-0.04)
!     ENDIF
!     ENDIF
!     ENDIF
!     ENDIF
    ENDIF
  ENDDO D9805
  end associate
  end subroutine CutBranchSheathPetole
!--------------------------------------------------------------------------------

  subroutine HarvestCanopy(I,J,NZ,HvstedLeafC,LeafC_lbrch)
  !
  !Canopy harvest
  !from layer to branch and to leaf node

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in):: HvstedLeafC
  REAL(R8), intent(in) :: LeafC_lbrch(NumOfCanopyLayers1,JP1,JP1)
  integer :: L,NB,K,NE
  real(r8) :: FHGT,FHVSH,FracHarvested
  real(r8) :: FrcLeafMassLeft
  real(r8) :: HvestedLeafCLayer_brch
  real(r8) :: FHVSHT

  associate(                                                              &
    CanopyHeightZ_col           => plt_morph%CanopyHeightZ_col,           &
    FracCanopyHeightCut_pft     => plt_distb%FracCanopyHeightCut_pft,     &
    FracBiomHarvsted            => plt_distb%FracBiomHarvsted,            &
    LeafElmsByLayerNode_brch    => plt_biom%LeafElmsByLayerNode_brch,     &
    CanopyStalkArea_lbrch       => plt_morph%CanopyStalkArea_lbrch,       &
    k_fine_litr                 => pltpar%k_fine_litr,                    &
    k_woody_litr                => pltpar%k_woody_litr,                   &
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft,           &
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft,             &
    CanopyLeafAreaZ_pft         => plt_morph%CanopyLeafAreaZ_pft,         &
    THIN_pft                    => plt_distb%THIN_pft,                    &
    LeafStrutElms_pft           => plt_biom%LeafStrutElms_pft,            &
    CanopyLeafCLyr_pft          => plt_biom%CanopyLeafCLyr_pft,           &
    CanopyLeafArea_lpft         => plt_morph%CanopyLeafArea_lpft,         &
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &
    CanopyStemAreaZ_pft         => plt_morph%CanopyStemAreaZ_pft,         &
    iHarvstType_pft             => plt_distb%iHarvstType_pft              &
  )
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     ZL=height to bottom of each canopy layer
!     FHGT=fraction of canopy layer height not harvested
!     FrcLeafMassLeft=fraction of canopy layer mass not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  D9865: DO L=NumOfCanopyLayers1,1,-1
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
        IF(CanopyHeightZ_col(L).GT.CanopyHeightZ_col(L-1))THEN
          FHGT=AZMAX1(AMIN1(1.0_r8,1._r8-((CanopyHeightZ_col(L))-FracCanopyHeightCut_pft(NZ))/ &
            (CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1))))
        ELSE
          FHGT=1.0_r8
        ENDIF
      ELSE
        FHGT=0._r8
      ENDIF
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FrcLeafMassLeft=AZMAX1(1._r8-(1._r8-FHGT)*FracBiomHarvsted(1,iplthvst_leaf,NZ))
        FHVSH=FrcLeafMassLeft
      ELSE
        FrcLeafMassLeft=AZMAX1(1._r8-THIN_pft(NZ))
        IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
          FHVSH=1.0_r8-(1._r8-FHGT)*FracBiomHarvsted(1,iplthvst_leaf,NZ)*THIN_pft(NZ)
        ELSE
          FHVSH=FrcLeafMassLeft
        ENDIF
      ENDIF
    ELSE
      FrcLeafMassLeft = 0._r8
      FHVSH           = 0._r8
    ENDIF
!
!     CUT LEAVES AT HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTLF=PFT leaf C mass
!     LeafC_lbrch=branch leaf C mass in canopy layer
!     WHVBSL,HvstedLeafC=layer,total leaf C mass removed
!     WGLFL=leaf node C in canopy layer
!     FrcLeafMassLeft=fraction of leaf node mass not harvested
! 
    D9855: DO NB=1,NumOfBranches_pft(NZ)
      IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
        .AND. LeafStrutElms_pft(ielmc,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        HvestedLeafCLayer_brch=HvstedLeafC*AZMAX1(LeafC_lbrch(L,NB,NZ))/LeafStrutElms_pft(ielmc,NZ)
      ELSE
        HvestedLeafCLayer_brch=0._r8
      ENDIF

      !distribute the harvest to leaf nodes on branch NB
      D9845: DO K=MaxNodesPerBranch1,0,-1
        IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
          .OR. HvestedLeafCLayer_brch.GT.0.0_r8)THEN
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
            IF(LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ).GT.HvestedLeafCLayer_brch)THEN
              FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,(LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)-HvestedLeafCLayer_brch) &
                /LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)))
              FHVSH=FrcLeafMassLeft
            ELSE
              FrcLeafMassLeft=1.0_r8
              FHVSH=1.0_r8
            ENDIF
          ENDIF
      !
!     HARVESTED LEAF AREA, C, N, P
!
!     FrcLeafMassLeft=fraction of leaf node mass not harvested
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanopyLeafArea_lpft,CanopyStalkArea_lbrch=leaf,stalk node area in canopy layer
!     LeafElmntRemoval(ielmc),LeafElmntRemoval(ielmn),LeafElmntRemoval(ielmp)=harvested leaf C,N,P
!     LeafElmntHarv2Litr(ielmc),LeafElmntHarv2Litr(ielmn),LeafElmntHarv2Litr(ielmp)=harvested leaf C,N,P to litter
!     WoodyElmntRemoval(ielmc),WoodyElmntRemoval(ielmn),WoodyElmntRemoval(ielmp)=harvested woody C,N,P
!     WoodyElmntHarv2Litr(ielmc),WoodyElmntHarv2Litr(ielmn),WoodyElmntHarv2Litr(ielmp)=harvested woody C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!

          HvestedLeafCLayer_brch=HvestedLeafCLayer_brch-(1._r8-FrcLeafMassLeft)*LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)
          FracHarvested=1._r8-FHVSH
          FHVSHT=FHVSH-FrcLeafMassLeft
          DO NE=1,NumPlantChemElms
            LeafElmntRemoval(NE)=LeafElmntRemoval(NE) &
              +FracHarvested*LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)
            LeafElmntHarv2Litr(NE)=LeafElmntHarv2Litr(NE) &
              +FHVSHT*LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)

            WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE) &
              +FracHarvested*LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
            WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE) &
              +FHVSHT*LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
            LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)=FrcLeafMassLeft*LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)
          ENDDO
!
!     REMAINING LEAF C,N,P AND AREA
!
          CanopyLeafArea_lpft(L,K,NB,NZ)=FrcLeafMassLeft*CanopyLeafArea_lpft(L,K,NB,NZ)
          IF(K.EQ.1)THEN
            CanopyStalkArea_lbrch(L,NB,NZ)=FrcLeafMassLeft*CanopyStalkArea_lbrch(L,NB,NZ)
          ENDIF
        ENDIF
      ENDDO D9845
    ENDDO D9855
    CanopyLeafAreaZ_pft(L,NZ) = 0._r8
    CanopyLeafCLyr_pft(L,NZ)  = 0._r8
    CanopyStemAreaZ_pft(L,NZ) = CanopyStemAreaZ_pft(L,NZ)*FrcLeafMassLeft
  ENDDO D9865
  end associate
  end subroutine HarvestCanopy
!--------------------------------------------------------------------------------

  subroutine StageBranch4Cut(I,J,NB,NZ,LeafCafCut_brch,LeafCbfCut_brch,FracIntnodeNotHvsted,FracNodeNotHvsted)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(out) :: LeafCafCut_brch,LeafCbfCut_brch
  real(r8), intent(out) :: FracIntnodeNotHvsted(0:MaxNodesPerBranch1)  
  real(r8), intent(out) :: FracNodeNotHvsted(0:MaxNodesPerBranch1)
  integer :: K,NE,L
  real(r8) :: ARLFG
  real(r8) :: LeafElmNodeK_brch(NumPlantChemElms)  

  associate(                                                       &
    CanopyLeafArea_lpft      => plt_morph%CanopyLeafArea_lpft,     &
    LeafElmntNode_brch       => plt_biom%LeafElmntNode_brch,       &
    ZERO4Groth_pft           => plt_biom%ZERO4Groth_pft,           &
    THIN_pft                 => plt_distb%THIN_pft,                &
    FracBiomHarvsted         => plt_distb%FracBiomHarvsted,        &
    LeafElmsByLayerNode_brch => plt_biom%LeafElmsByLayerNode_brch, &
    CanopyLeafCLyr_pft       => plt_biom%CanopyLeafCLyr_pft,       &
    LeafAreaLive_brch        => plt_morph%LeafAreaLive_brch,       &
    LeafAreaNode_brch        => plt_morph%LeafAreaNode_brch,       &
    LeafProteinCNode_brch    => plt_biom%LeafProteinCNode_brch,    &
    LeafStrutElms_brch       => plt_biom%LeafStrutElms_brch,       &
    CanopyLeafAreaZ_pft      => plt_morph%CanopyLeafAreaZ_pft,     &
    iHarvstType_pft          => plt_distb%iHarvstType_pft          &
  )
  LeafCafCut_brch=0._r8
  LeafCbfCut_brch=0._r8

  D9825: DO K=0,MaxNodesPerBranch1
    ARLFG=0._r8
    LeafElmNodeK_brch(1:NumPlantChemElms)=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanopyLeafArea_lpft,CanopyLeafAreaZ_pft=leaf node,total area in canopy layer
!
    D9815: DO L=1,NumOfCanopyLayers1
      ARLFG=ARLFG+CanopyLeafArea_lpft(L,K,NB,NZ)
      DO NE=1,NumPlantChemElms
        LeafElmNodeK_brch(NE)=LeafElmNodeK_brch(NE)+LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)
      ENDDO
      CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)+CanopyLeafArea_lpft(L,K,NB,NZ)
      CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)+LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)
    ENDDO D9815
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FracIntnodeNotHvsted=fraction of internode layer mass not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ).AND.FracBiomHarvsted(1,iplthvst_leaf,NZ).GT.0.0)THEN
        FracIntnodeNotHvsted(K)=AZMAX1(AMIN1(1.0_r8,(1._r8-(1._r8-AZMAX1(LeafElmNodeK_brch(ielmc)) &
          /LeafElmntNode_brch(ielmc,K,NB,NZ))*FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)&
          /FracBiomHarvsted(1,iplthvst_leaf,NZ))))
        FracNodeNotHvsted(K)=FracIntnodeNotHvsted(K)
      ELSE
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FracIntnodeNotHvsted(K) = 1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)
          FracNodeNotHvsted(K)    = FracIntnodeNotHvsted(K)
        ELSE
          FracIntnodeNotHvsted(K)=1.0_r8-THIN_pft(NZ)
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
            FracNodeNotHvsted(K)=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
          ELSE
            FracNodeNotHvsted(K)=FracIntnodeNotHvsted(K)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      FracIntnodeNotHvsted(K) = 0._r8
      FracNodeNotHvsted(K)    = 0._r8
    ENDIF
!
!     ACCUMULATE REMAINING BRANCH LEAF AREA, C, N, P
!
!     WGLF=leaf node C mass
!     LeafAreaLive_brch,LeafAreaNode_brch=branch,node leaf area
!     LeafProteinCNode_brch=leaf protein mass
!
    LeafCbfCut_brch=LeafCbfCut_brch+LeafElmntNode_brch(ielmc,K,NB,NZ)
    DO NE=1,NumPlantChemElms
      LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ)+LeafElmNodeK_brch(NE)
    ENDDO
    LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-LeafAreaNode_brch(K,NB,NZ)+ARLFG
    IF(LeafAreaNode_brch(K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*ARLFG/LeafAreaNode_brch(K,NB,NZ)
    ELSE
      LeafProteinCNode_brch(K,NB,NZ)=0._r8
    ENDIF
    LeafAreaNode_brch(K,NB,NZ)=ARLFG

    DO NE=1,NumPlantChemElms
      LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmNodeK_brch(NE)
    ENDDO
    LeafCafCut_brch=LeafCafCut_brch+LeafElmntNode_brch(ielmc,K,NB,NZ)
  ENDDO D9825
  end associate
  end subroutine StageBranch4Cut    

!--------------------------------------------------------------------------------
  subroutine CutPlant(I,J,NZ,WHVSHH,WHVSCP,WHVSNP,HvstedShethC,HvstedGrainC,HvstedEarC,&
    HvstedRsrvC,HvstedStalkC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: WHVSHH,WHVSCP,WHVSNP
  real(r8), intent(in) :: HvstedShethC,HvstedGrainC,HvstedEarC,HvstedRsrvC,HvstedStalkC

  integer :: NB
  real(r8) :: FracNodeNotHvsted(0:MaxNodesPerBranch1),FracIntnodeNotHvsted(0:MaxNodesPerBranch1)
  real(r8) :: LeafCafCut_brch,WGSHGX,LeafCbfCut_brch,WGSHGY
  real(r8) :: VOLWPX,WVPLT  
  real(r8) :: RMedInternodeLen
  real(r8) :: FDM  
  associate(                                                          &
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,     &
    CanopyWater_pft           => plt_ew%CanopyWater_pft,              &
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,        &
    jHarvst_pft               => plt_distb%jHarvst_pft,               &
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft, &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,     &
    iPlantBranchState_brch    => plt_pheno%iPlantBranchState_brch,    &
    QH2OLoss_lnds             => plt_site%QH2OLoss_lnds,              &
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,       &
    FracCanopyHeightCut_pft   => plt_distb%FracCanopyHeightCut_pft,   &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,                &
    CanopyStalkC_pft          => plt_biom%CanopyStalkC_pft,           &
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,       &
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,         &
    H2OLoss_CumYr_col         => plt_ew%H2OLoss_CumYr_col,            &
    CanopyHeight_pft          => plt_morph%CanopyHeight_pft,          &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,         &
    iHarvstType_pft           => plt_distb%iHarvstType_pft            &
  )

  D9835: DO NB=1,NumOfBranches_pft(NZ)

    CALL StageBranch4Cut(I,J,NB,NZ,LeafCafCut_brch,LeafCbfCut_brch,FracIntnodeNotHvsted,FracNodeNotHvsted)

    call CutBranchSheathPetole(I,J,NB,NZ,WHVSHH,FracIntnodeNotHvsted,FracNodeNotHvsted,RMedInternodeLen,WGSHGX,WGSHGY)

    call CutBranchNonstructural(I,J,NB,NZ,LeafCafCut_brch,WGSHGX,LeafCbfCut_brch,WGSHGY,WHVSCP,WHVSNP)

  !
  !     CUT STALKS
    call BranchCutPlantStalk(I,J,NB,NZ,RMedInternodeLen,HvstedStalkC,HvstedRsrvC)
!
!     CUT REPRODUCTIVE ORGANS FracHuskNotHvsted
    call BranchCutReprodOrgans(I,J,NB,NZ,RMedInternodeLen,HvstedShethC,HvstedGrainC,HvstedEarC)
!
!     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
!
!     ShootC4NonstC_brch=total C4 nonstructural C in branch
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     StalkBiomassC_brch=stalk sapwood mass
!     PSICanopy_pft=canopy water potential
!     CanopyWater_pft=water volume in canopy
!     QH2OLoss_lnds,H2OLoss_CumYr_col=accumulated water loss for water balance calculation
!
    LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))

    call SumPlantBranchBiome(NB,NZ)

    VOLWPX              = CanopyWater_pft(NZ)
    WVPLT               = AZMAX1(CanopyLeafShethC_pft(NZ)+CanopyStalkC_pft(NZ))
    FDM                 = get_FDM(PSICanopy_pft(NZ))    !drymatter/water=fdm
    CanopyWater_pft(NZ) = ppmc*WVPLT/FDM
    QH2OLoss_lnds       = QH2OLoss_lnds+VOLWPX-CanopyWater_pft(NZ)
    H2OLoss_CumYr_col   = H2OLoss_CumYr_col+VOLWPX-CanopyWater_pft(NZ)
!
!     RESET PHENOLOGY, GROWTH STAGE IF STALKS ARE CUT
!
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire

    IF((iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))) &
      .AND. (iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
      .AND. CanopyHeight_pft(NZ).GT.FracCanopyHeightCut_pft(NZ))THEN
      call ResetCutBranch(I,J,NZ,NB)
    ENDIF
!
!     DEATH OF BRANCH IF KILLING HARVEST ENTERED IN 'READQ'
!
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,and reseed
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     PP=PFT population
!     WTLS=total PFT leaf+petiole C mass
!     WTSTK=total PFT stalk C mass
!     WVSTK=total PFT sapwood C mass
!     CanopyStalkArea_lbrch=total PFT stalk surface area
!
    IF(jHarvst_pft(NZ).NE.jharvtyp_noaction)then
      iPlantBranchState_brch(NB,NZ)=iDead      
    endif

    IF(PlantPopulation_pft(NZ).LE.0.0_r8)then
      iPlantBranchState_brch(NB,NZ)=iDead
    endif
  ENDDO D9835
  end associate
  end subroutine CutPlant

!--------------------------------------------------------------------------------
  subroutine RootMaterialRemovalL(I,J,N,L,NZ,FracLeftThin,XHVST1)

  implicit none
  integer , intent(in) :: I,J,N,L,NZ
  real(r8), intent(out) :: FracLeftThin
  real(r8), intent(out):: XHVST1
  real(r8) :: FrcMassNotHarvst(NumPlantChemElms)
  real(r8) :: FFIRE(NumPlantChemElms)
  integer  :: NR,NTG,M,NE

  associate(                                                         &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    Eco_NBP_CumYr_col         => plt_bgcr%Eco_NBP_CumYr_col,         &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,    &
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,           &
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,           &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    THIN_pft                  => plt_distb%THIN_pft,                 &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,          &
    THETW_vr                  => plt_soilchem%THETW_vr,              &
    DCORP                     => plt_distb%DCORP,                    &
    FracBiomHarvsted          => plt_distb%FracBiomHarvsted,         &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,    &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    k_woody_litr              => pltpar%k_woody_litr,                &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    iroot                     => pltpar%iroot,                       &
    inonstruct                => pltpar%inonstruct,                  &
    icwood                    => pltpar%icwood,                      &
    iHarvstType_pft           => plt_distb%iHarvstType_pft           &
  )

  IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
    FracLeftThin=1.0_r8-THIN_pft(NZ)
    FFIRE(1:NumPlantChemElms)=0._r8
  ELSE
    !fire
    call StageRootRemovalByFire(I,J,NZ,L,FFIRE,DCORP,FracLeftThin)
  ENDIF

  XHVST1=1._r8-FracLeftThin
  D3385: DO M=1,jsken
    DO NE=1,NumPlantChemElms
      FrcMassNotHarvst(NE)=XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ))
      LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
        +(1._r8-FFIRE(NE))*FrcMassNotHarvst(NE)
    ENDDO
    !nonstructural root biomass
    call RemoveRootByFire(I,J,NZ,FrcMassNotHarvst,FFIRE)

    DO NR=1,NumRootAxes_pft(NZ)
      DO NE=1,NumPlantChemElms
        FrcMassNotHarvst(NE)=XHVST1*ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
          +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_woody_litr)
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+&
          (1._r8-FFIRE(NE))*FrcMassNotHarvst(NE)       
      ENDDO
      !woody roots
      call RemoveRootByFire(I,J,NZ,FrcMassNotHarvst,FFIRE)

      DO NE=1,NumPlantChemElms
        FrcMassNotHarvst(NE)=XHVST1*ElmAllocmat4Litr(NE,iroot,M,NZ)*AZMAX1(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
          +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_fine_litr)
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
          +(1._r8-FFIRE(NE))*FrcMassNotHarvst(NE)
      ENDDO
      !fine roots
      CALL RemoveRootByFire(I,J,NZ,FrcMassNotHarvst,FFIRE)

    enddo
  ENDDO D3385
!
!     RELEASE ROOT GAS CONTENTS DURING HARVESTING
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
  DO NTG=idg_beg,idg_end-1
    RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-XHVST1 &
      *(trcg_rootml_pvr(idg_CO2,N,L,NZ)+trcs_rootml_pvr(idg_CO2,N,L,NZ))
    trcg_rootml_pvr(NTG,N,L,NZ)=FracLeftThin*trcg_rootml_pvr(NTG,N,L,NZ)
    trcs_rootml_pvr(NTG,N,L,NZ)=FracLeftThin*trcs_rootml_pvr(NTG,N,L,NZ)
  ENDDO
  end associate          
  end subroutine RootMaterialRemovalL

!--------------------------------------------------------------------------------
  subroutine HarvstUpdateRootStateL(I,J,N,L,NZ,FracLeftThin,XHVST1)            
  implicit none
  integer,  intent(in) :: I,J,N,L,NZ
  real(r8), intent(in) :: FracLeftThin,XHVST1
  integer :: NE,NR,M
  associate(                                                         &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,          &
    inonstruct                => pltpar%inonstruct,                  &
    iroot                     => pltpar%iroot,                       &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,       &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    k_woody_litr              => pltpar%k_woody_litr,                &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    Root2ndLen_rpvr            => plt_morph%Root2ndLen_rpvr,           &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr,         &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,          &
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr,         &
    Root2ndXNum_pvr           => plt_morph%Root2ndXNum_pvr,          &
    RootLenPerPlant_pvr       => plt_morph%RootLenPerPlant_pvr,      &
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr,  &
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr,          &
    RootAreaPerPlant_pvr      => plt_morph%RootAreaPerPlant_pvr,     &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,        &
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr,             &
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr,          &
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr,        &
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr,        &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    RootNodulStrutElms_rpvr    => plt_biom%RootNodulStrutElms_rpvr,    &
    RootNodulNonstElms_rpvr    => plt_biom%RootNodulNonstElms_rpvr,    &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr    &
  )
!
!     REDUCE ROOT STATE VARIABLES DURING HARVESTING
!
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stLen_rpvr,Root2ndLen_rpvr=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootMycoActiveBiomC_pvr, PopuRootMycoC_pvr=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,Root2ndXNum_pvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_pvr,RootPoreVol_pvr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_pvr=root surface area per plant
!     RootRespPotent_pvr,RootCO2EmisPot_pvr,RootCO2Autor_pvr unlimited by O2,nonstructural C
!

    D3960: DO NR=1,NumRootAxes_pft(NZ)
      DO NE=1,NumPlantChemElms
        RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)*FracLeftThin
        RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*FracLeftThin
      ENDDO
      Root1stLen_rpvr(N,L,NR,NZ)=Root1stLen_rpvr(N,L,NR,NZ)*FracLeftThin
      Root2ndLen_rpvr(N,L,NR,NZ)=Root2ndLen_rpvr(N,L,NR,NZ)*FracLeftThin
      Root2ndXNum_rpvr(N,L,NR,NZ)=Root2ndXNum_rpvr(N,L,NR,NZ)*FracLeftThin
    ENDDO D3960

    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)*FracLeftThin
    ENDDO
    RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ)*FracLeftThin
    PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ)*FracLeftThin
    RootProteinC_pvr(N,L,NZ)        = RootProteinC_pvr(N,L,NZ)*FracLeftThin
    Root1stXNumL_pvr(N,L,NZ)        = Root1stXNumL_pvr(N,L,NZ)*FracLeftThin
    Root2ndXNum_pvr(N,L,NZ)         = Root2ndXNum_pvr(N,L,NZ)*FracLeftThin
    RootLenPerPlant_pvr(N,L,NZ)     = RootLenPerPlant_pvr(N,L,NZ)*FracLeftThin
    RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ)*FracLeftThin
    RootPoreVol_pvr(N,L,NZ)         = RootPoreVol_pvr(N,L,NZ)*FracLeftThin
    RootVH2O_pvr(N,L,NZ)            = RootVH2O_pvr(N,L,NZ)*FracLeftThin
    RootAreaPerPlant_pvr(N,L,NZ)    = RootAreaPerPlant_pvr(N,L,NZ)*FracLeftThin
    RootRespPotent_pvr(N,L,NZ)      = RootRespPotent_pvr(N,L,NZ)*FracLeftThin
    RootCO2EmisPot_pvr(N,L,NZ)      = RootCO2EmisPot_pvr(N,L,NZ)*FracLeftThin
    RootCO2Autor_pvr(N,L,NZ)        = RootCO2Autor_pvr(N,L,NZ)*FracLeftThin
!
!     NODULE LitrFall AND STATE VARIABLES DURING HARVESTING
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
    IF(is_plant_N2fix(iPlantNfixType_pft(NZ)).AND.N.EQ.ipltroot)THEN
      DO NE=1,NumPlantChemElms
        D3395: DO M=1,jsken
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+ &
            XHVST1*AZMAX1(ElmAllocmat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_rpvr(NE,L,NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_rpvr(NE,L,NZ))
        ENDDO D3395
        RootNodulStrutElms_rpvr(NE,L,NZ)=RootNodulStrutElms_rpvr(NE,L,NZ)*FracLeftThin
        RootNodulNonstElms_rpvr(NE,L,NZ)=RootNodulNonstElms_rpvr(NE,L,NZ)*FracLeftThin
      ENDDO
    ENDIF
  end associate          
  end subroutine HarvstUpdateRootStateL

end module PlantDisturbsMod
