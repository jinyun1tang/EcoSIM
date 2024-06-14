module PlantDisturbsMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  use PlantMathFuncMod
  use PlantBalMod, only : SumPlantBranchBiome
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
  real(r8) :: WTHTGE(NumPlantChemElms)
  real(r8) :: EFIRE(2,5:5)

  public :: RemoveBiomassByDisturbance
  public :: RemoveBiomByMgmt
  public :: InitPlantDisturbance
  contains

  subroutine InitPlantDisturbance
  implicit none
  EFIRE=reshape(real((/0.917,0.167/),r8),shape(EFIRE))

  end subroutine InitPlantDisturbance

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByMgmt(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

!     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
!

  call RemoveBiomByHarvest(I,J,NZ)
!
!     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
!
  call RemoveBiomByTillage(I,J,NZ)
  end subroutine RemoveBiomByMgmt
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomassByDisturbance(I,J,NZ)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8) :: FHVSE
  real(r8) :: FHVSH
  real(r8) :: WHVSTD
  integer :: M,NE

!     begin_execution
  associate(                                                    &
    FracBiomHarvsted       => plt_distb%FracBiomHarvsted,       &
    CutHeightORFrac_pft => plt_distb%CutHeightORFrac_pft, &
    iHarvstType_pft        => plt_distb%iHarvstType_pft,        &
    THIN_pft               => plt_distb%THIN_pft,               &
    PlantPopulation_pft    => plt_site%PlantPopulation_pft,     &
    NU                     => plt_site%NU,                      &
    ZERO                   => plt_site%ZERO,                    &
    SolarNoonHour_col      => plt_site%SolarNoonHour_col,       &
    AREA3                  => plt_site%AREA3,                   &
    ZERO4Uptk_pft          => plt_rbgc%ZERO4Uptk_pft,           &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft,          &
    ZERO4LeafVar_pft       => plt_biom%ZERO4LeafVar_pft,        &
    ShootC4NonstC_brch     => plt_biom%ShootC4NonstC_brch,      &
    StandDeadKCompElms_pft => plt_biom%StandDeadKCompElms_pft,  &
    StandDeadStrutElms_pft => plt_biom%StandDeadStrutElms_pft   &
  )
!  write(102,*)'iHarvstType_pft',I,iHarvstType_pft(NZ),plt_distb%jHarvst_pft(NZ),NZ
  NonstructElmntRemoval(1:NumPlantChemElms)=0._r8
  LeafElmntRemoval(1:NumPlantChemElms)=0._r8
  FineNonleafElmntRemoval(1:NumPlantChemElms)=0._r8
  WoodyElmntRemoval(1:NumPlantChemElms)=0._r8
  StandeadElmntRemoval(1:NumPlantChemElms)=0._r8
  LeafElmnt2Litr(1:NumPlantChemElms)=0._r8
  FineNonleafElmnt2Litr(1:NumPlantChemElms)=0._r8
  WoodyElmnt2Litr(1:NumPlantChemElms)=0._r8
  StandeadElmnt2Litr(1:NumPlantChemElms)=0._r8
  LeafElmntHarv2Litr(1:NumPlantChemElms)=0._r8
  PetioleElmntHarv2Litr(1:NumPlantChemElms)=0._r8
  WoodyElmntHarv2Litr(1:NumPlantChemElms)=0._r8
  StandeadElmntHarv2Litr(1:NumPlantChemElms)=0._r8
  WTHTGE(1:NumPlantChemElms)=0._r8
!     ENDIF

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
    IF(J.EQ.INT(SolarNoonHour_col).AND.iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.&
      iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
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
    ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
      IF(StandDeadStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
        WHVSTD=CutHeightORFrac_pft(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8*AREA3(NU)*FracBiomHarvsted(1,4,NZ)
        FHVSE=AZMAX1(1._r8-WHVSTD/StandDeadStrutElms_pft(ielmc,NZ))
        FHVSH=FHVSE
      ELSE
        FHVSE=1.0_r8
        FHVSH=1.0_r8
      ENDIF
    ELSE
      FHVSE=1.0_r8
      FHVSH=1.0_r8
    ENDIF
    D6475: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        StandeadElmntRemoval(NE)=StandeadElmntRemoval(NE)+(1._r8-FHVSH)*StandDeadKCompElms_pft(NE,M,NZ)
        StandeadElmntHarv2Litr(NE)=StandeadElmntHarv2Litr(NE)+(FHVSH-FHVSE)*StandDeadKCompElms_pft(NE,M,NZ)
        StandDeadKCompElms_pft(NE,M,NZ)=FHVSE*StandDeadKCompElms_pft(NE,M,NZ)
      ENDDO
    ENDDO D6475
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
  real(r8) :: FineNonleafElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8) :: HarvestElmnt2Litr(NumPlantChemElms)

  NonstructElmntOffEcosystem(1:NumPlantChemElms)=0._r8
  LeafElmntOffEcosystem(1:NumPlantChemElms)=0._r8
  FineNonleafElmntOffEcosystem(1:NumPlantChemElms)=0._r8
  WoodyElmntOffEcosystem(1:NumPlantChemElms)=0._r8
  StandeadElmntOffEcosystem(1:NumPlantChemElms)=0._r8
  NonstructElmnt2Litr(1:NumPlantChemElms)=0._r8

  call ApplyDisturbanceBiomRemoval(I,J,NZ,NonstructElmnt2Litr,NonstructElmntOffEcosystem,&
    LeafElmntOffEcosystem,FineNonleafElmntOffEcosystem,WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call TotalBiomRemovalByDisturbance(I,J,NZ,NonstructElmnt2Litr,HarvestElmnt2Litr,TotalElmnt2Litr)
!
!     ABOVE-GROUND LitrFall FROM HARVESTING
!
  call LiterfallByDisturbance(I,J,NZ,HarvestElmnt2Litr,TotalElmnt2Litr,NonstructElmnt2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmntOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
  end subroutine PlantDisturbance
!------------------------------------------------------------------------------------------

  subroutine LiterfallByDisturbance(I,J,NZ,HarvestElmnt2Litr,TotalElmnt2Litr,NonstructElmnt2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmntOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)

  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: HarvestElmnt2Litr(NumPlantChemElms),TotalElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  integer :: M,NE
!     begin_execution
  associate(                                                             &
    ilignin                     => pltpar%ilignin,                       &
    inonstruct                  => pltpar%inonstruct,                    &
    k_fine_litr                 => pltpar%k_fine_litr,                   &
    k_woody_litr                => pltpar%k_woody_litr,                  &
    ifoliar                     => pltpar%ifoliar,                       &
    inonfoliar                  => pltpar%inonfoliar,                    &
    istalk                      => pltpar%istalk,                        &
    iroot                       => pltpar%iroot,                         &
    icwood                      => pltpar%icwood,                        &
    StandDeadKCompElms_pft      => plt_biom%StandDeadKCompElms_pft,      &
    iHarvstType_pft             => plt_distb%iHarvstType_pft,            &
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr, &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,        &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,        &
    LitrfalStrutElmsCum_pft     => plt_bgcr%LitrfalStrutElmsCum_pft,     &
    SurfLitrfalStrutElmsCum_pft => plt_bgcr%SurfLitrfalStrutElmsCum_pft, &
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft,  &
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft       &
  )
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      D6375: DO M=1,jsken
        DO NE=1,NumPlantChemElms        
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*(NonstructElmnt2Litr(NE)) &
            +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafElmnt2Litr(NE)+LeafElmntHarv2Litr(NE)) &
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(NE)+PetioleElmntHarv2Litr(NE))

          IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,istalk,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmnt2Litr(NE)&
              +StandeadElmntHarv2Litr(NE))
          ELSE
            StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*(WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE))

            LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE) &
              +StandeadElmnt2Litr(NE))*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE) &
              +StandeadElmnt2Litr(NE))*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
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
      D6485: DO M=1,jsken
        LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(ielmc,inonstruct,M,NZ)*NonstructElmnt2Litr(ielmc) &
          +ElmAllocmat4Litr(ielmc,ifoliar,M,NZ)   *LeafElmnt2Litr(ielmc)+LeafElmntHarv2Litr(ielmc) &
          +ElmAllocmat4Litr(ielmc,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmc)+PetioleElmntHarv2Litr(ielmc))

        LitrfalStrutElms_pvr(ielmn,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,M,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(ielmn,inonstruct,M,NZ)*NonstructElmntOffEcosystem(ielmn) &
          +ElmAllocmat4Litr(ielmn,ifoliar,M,NZ)   *LeafElmntOffEcosystem(ielmn) &
          +ElmAllocmat4Litr(ielmn,inonfoliar,M,NZ)*FineNonleafElmntOffEcosystem(ielmn)

        LitrfalStrutElms_pvr(ielmp,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,M,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(ielmp,inonstruct,M,NZ)*NonstructElmntOffEcosystem(ielmp) &
          +ElmAllocmat4Litr(ielmp,ifoliar,M,NZ)   *LeafElmntOffEcosystem(ielmp) &
          +ElmAllocmat4Litr(ielmp,inonfoliar,M,NZ)*FineNonleafElmntOffEcosystem(ielmp)

        LitrfalStrutElms_pvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,ilignin,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(ielmn,inonstruct,M,NZ)*(NonstructElmnt2Litr(ielmn)-NonstructElmntOffEcosystem(ielmn)) &
          +ElmAllocmat4Litr(ielmn,ifoliar,M,NZ)   *(LeafElmnt2Litr(ielmn)+LeafElmntHarv2Litr(ielmn)-LeafElmntOffEcosystem(ielmn)) &
          +ElmAllocmat4Litr(ielmn,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmn)+PetioleElmntHarv2Litr(ielmn) &
          -FineNonleafElmntOffEcosystem(ielmn))

        LitrfalStrutElms_pvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,ilignin,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(ielmp,inonstruct,M,NZ)*(NonstructElmnt2Litr(ielmp)-NonstructElmntOffEcosystem(ielmp)) &
          +ElmAllocmat4Litr(ielmp,ifoliar,M,NZ)*(LeafElmnt2Litr(ielmp)+LeafElmntHarv2Litr(ielmp)-LeafElmntOffEcosystem(ielmp)) &
          +ElmAllocmat4Litr(ielmp,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmp)+PetioleElmntHarv2Litr(ielmp)&
          -FineNonleafElmntOffEcosystem(ielmp))

        IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
          LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)+&
            ElmAllocmat4Litr(ielmc,istalk,M,NZ)*(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc) &
            +StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))

          LitrfalStrutElms_pvr(ielmn,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,M,k_fine_litr,0,NZ)+&
            ElmAllocmat4Litr(ielmn,istalk,M,NZ)*(WoodyElmntOffEcosystem(ielmn)+StandeadElmntOffEcosystem(ielmn))
          LitrfalStrutElms_pvr(ielmp,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,M,k_fine_litr,0,NZ)+&
            ElmAllocmat4Litr(ielmp,istalk,M,NZ)*(WoodyElmntOffEcosystem(ielmp)+StandeadElmntOffEcosystem(ielmp))

          LitrfalStrutElms_pvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,ilignin,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmn,istalk,M,NZ)*(WoodyElmnt2Litr(ielmn)+WoodyElmntHarv2Litr(ielmn) &
            -WoodyElmntOffEcosystem(ielmn)+StandeadElmnt2Litr(ielmn)+StandeadElmntHarv2Litr(ielmn)&
            -StandeadElmntOffEcosystem(ielmn))

          LitrfalStrutElms_pvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,ilignin,k_fine_litr,0,NZ)+&
            ElmAllocmat4Litr(ielmp,istalk,M,NZ)*(WoodyElmnt2Litr(ielmp)+WoodyElmntHarv2Litr(ielmp)- &
            WoodyElmntOffEcosystem(ielmp)+StandeadElmnt2Litr(ielmp)+StandeadElmntHarv2Litr(ielmp)-&
            StandeadElmntOffEcosystem(ielmp))
        ELSE
          StandDeadKCompElms_pft(ielmc,M,NZ)=StandDeadKCompElms_pft(ielmc,M,NZ)+ElmAllocmat4Litr(ielmc,icwood,M,NZ) &
            *(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc))
          StandDeadKCompElms_pft(ielmn,M,NZ)=StandDeadKCompElms_pft(ielmn,M,NZ)+ElmAllocmat4Litr(ielmn,icwood,M,NZ) &
            *WoodyElmntOffEcosystem(ielmn)
          StandDeadKCompElms_pft(ielmp,M,NZ)=StandDeadKCompElms_pft(ielmp,M,NZ)+ElmAllocmat4Litr(ielmp,icwood,M,NZ) &
            *WoodyElmntOffEcosystem(ielmp)
            
          LitrfalStrutElms_pvr(ielmc,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_woody_litr,0,NZ) &
            *ElmAllocmat4Litr(ielmc,istalk,M,NZ)&
            *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FracRootStalkElmAlloc2Litr(ielmc,k_woody_litr)
          LitrfalStrutElms_pvr(ielmn,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmn,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmn)*FracRootStalkElmAlloc2Litr(ielmn,k_woody_litr)
          LitrfalStrutElms_pvr(ielmp,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmp,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmp)*FracRootStalkElmAlloc2Litr(ielmp,k_woody_litr)

          LitrfalStrutElms_pvr(ielmn,ilignin,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,ilignin,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmn,icwood,M,NZ) &
            *(WoodyElmnt2Litr(ielmn)+WoodyElmntHarv2Litr(ielmn)-WoodyElmntOffEcosystem(ielmn) &
            +StandeadElmnt2Litr(ielmn)+StandeadElmntHarv2Litr(ielmn)-StandeadElmntOffEcosystem(ielmn)) &
            *FracRootStalkElmAlloc2Litr(ielmn,k_woody_litr)
          LitrfalStrutElms_pvr(ielmp,ilignin,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,ilignin,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmp,icwood,M,NZ) &
            *(WoodyElmnt2Litr(ielmp)+WoodyElmntHarv2Litr(ielmp)-WoodyElmntOffEcosystem(ielmp) &
            +StandeadElmnt2Litr(ielmp)+StandeadElmntHarv2Litr(ielmp)-StandeadElmntOffEcosystem(ielmp)) &
            *FracRootStalkElmAlloc2Litr(ielmp,k_woody_litr)

          LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmc,istalk,M,NZ) &
            *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FracRootStalkElmAlloc2Litr(ielmc,k_fine_litr)
          LitrfalStrutElms_pvr(ielmn,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmn,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmn)*FracRootStalkElmAlloc2Litr(ielmn,k_fine_litr)
          LitrfalStrutElms_pvr(ielmp,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmp,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmp)*FracRootStalkElmAlloc2Litr(ielmp,k_fine_litr)
            
          LitrfalStrutElms_pvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmn,ilignin,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmn,icwood,M,NZ)*(WoodyElmnt2Litr(ielmn)+WoodyElmntHarv2Litr(ielmn)-WoodyElmntOffEcosystem(ielmn) &
            +StandeadElmnt2Litr(ielmn)+StandeadElmntHarv2Litr(ielmn)-StandeadElmntOffEcosystem(ielmn)) &
            *FracRootStalkElmAlloc2Litr(ielmn,k_fine_litr)
          LitrfalStrutElms_pvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmp,ilignin,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(ielmp,icwood,M,NZ)*(WoodyElmnt2Litr(ielmp)+WoodyElmntHarv2Litr(ielmp)-WoodyElmntOffEcosystem(ielmp) &
            +StandeadElmnt2Litr(ielmp)+StandeadElmntHarv2Litr(ielmp)-StandeadElmntOffEcosystem(ielmp))&
            *FracRootStalkElmAlloc2Litr(ielmp,k_fine_litr)
        ENDIF
      ENDDO D6485
    ENDIF
  ELSE
!
!     ABOVE-GROUND LitrFall FROM GRAZING
!
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P LitrFall
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P LitrFall
!
    DO NE=1,NumPlantChemElms
      LitrfalStrutElmsCum_pft(NE,NZ)=LitrfalStrutElmsCum_pft(NE,NZ)+TotalElmnt2Litr(NE)+HarvestElmnt2Litr(NE)
      SurfLitrfalStrutElmsCum_pft(NE,NZ)=SurfLitrfalStrutElmsCum_pft(NE,NZ)+TotalElmnt2Litr(NE)+HarvestElmnt2Litr(NE)
    ENDDO
  ENDIF
  end associate
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

  subroutine TotalBiomRemovalByDisturbance(I,J,NZ,NonstructElmnt2Litr,HarvestElmnt2Litr,TotalElmnt2Litr)
  implicit none
  integer , intent(in)  :: I,J,NZ
  real(r8), intent(in)  :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: HarvestElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8) :: TotalElmntRemoval(NumPlantChemElms)
  integer :: NE
!     begin_execution
  associate(                                              &
    iHarvstType_pft       => plt_distb%iHarvstType_pft,   &
    jHarvst_pft           => plt_distb%jHarvst_pft,       &
    EcoHavstElmnt_pft     => plt_distb%EcoHavstElmnt_pft, &
    EcoHavstElmnt_col     => plt_distb%EcoHavstElmnt_col, &
    NH3byFire_pft         => plt_distb%NH3byFire_pft,     &
    PO4byFire_pft         => plt_distb%PO4byFire_pft,     &
    CH4ByFire_pft         => plt_distb%CH4ByFire_pft,     &
    O2ByFire_pft          => plt_distb%O2ByFire_pft,      &
    N2ObyFire_pft         => plt_distb%N2ObyFire_pft,     &
    CO2ByFire_pft         => plt_distb%CO2ByFire_pft,     &
    CO2NetFix_pft         => plt_bgcr%CO2NetFix_pft,      &
    Eco_NBP_col           => plt_bgcr%Eco_NBP_col,        &
    Eco_AutoR_col         => plt_bgcr%Eco_AutoR_col,      &
    ECO_ER_col            => plt_bgcr%ECO_ER_col,         &
    CanopyRespC_pft       => plt_bgcr%CanopyRespC_pft,    &
    GrossResp_pft         => plt_bgcr%GrossResp_pft,      &
    SeasonalNonstElms_pft => plt_biom%SeasonalNonstElms_pft     &
  )
!
!     TotalElmntRemoval(ielmc),TotalElmntRemoval(ielmn),TotalElmntRemoval(ielmp)=total C,N,P removed
!     TotalElmnt2Litr(ielmc),TotalElmnt2Litr(ielmn),TotalElmnt2Litr(ielmp)=total C,N,P to litter
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  DO NE=1,NumPlantChemElms
    TotalElmntRemoval(NE)=NonstructElmntRemoval(NE)+LeafElmntRemoval(NE)+FineNonleafElmntRemoval(NE)+&
      WoodyElmntRemoval(NE)+StandeadElmntRemoval(NE)
    TotalElmnt2Litr(NE)=NonstructElmnt2Litr(NE)+LeafElmnt2Litr(NE)+FineNonleafElmnt2Litr(NE)+&
      WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE)
    HarvestElmnt2Litr(NE)=LeafElmntHarv2Litr(NE)+PetioleElmntHarv2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE)
  ENDDO

  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      IF(jHarvst_pft(NZ).NE.jharvtyp_tmareseed)THEN
        DO NE=1,NumPlantChemElms
          EcoHavstElmnt_pft(NE,NZ)=EcoHavstElmnt_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
          EcoHavstElmnt_col(NE)=EcoHavstElmnt_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
        ENDDO
        Eco_NBP_col=Eco_NBP_col+TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc)
      ELSE
        DO NE=1,NumPlantChemElms
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
        ENDDO
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     CO2ByFire_pft,CH4ByFire_pft,O2ByFire_pft,NH3byFire_pft,N2ObyFire_pft,PO4byFire_pft=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_col=total net biome productivity
!
    ELSE
      CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
      CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FrcAsCH4byFire*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
      O2ByFire_pft(NZ)=O2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))*2.667_r8
      NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-TotalElmntRemoval(ielmn)+TotalElmnt2Litr(ielmn)
      N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0_r8
      PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-TotalElmntRemoval(ielmp)+TotalElmnt2Litr(ielmp)
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
      Eco_NBP_col=Eco_NBP_col-FrcAsCH4byFire*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
    ENDIF

  ELSE
!
!     C,N,P REMOVED FROM GRAZING
!
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     GY=growth yield of grazers
!     TotalElmntRemoval(ielmc),TotalElmntRemoval(ielmn),TotalElmntRemoval(ielmp)=total C,N,P removed
!     TotalElmnt2Litr(ielmc),TotalElmnt2Litr(ielmn),TotalElmnt2Litr(ielmp)=total C,N,P to litter
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!
    EcoHavstElmnt_pft(ielmc,NZ)=EcoHavstElmnt_pft(ielmc,NZ)+GY*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
    EcoHavstElmnt_col(ielmc)=EcoHavstElmnt_col(ielmc)+GY*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
    DO NE=2,NumPlantChemElms
      EcoHavstElmnt_pft(NE,NZ)=EcoHavstElmnt_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
      EcoHavstElmnt_col(NE)=EcoHavstElmnt_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
    ENDDO
    GrossResp_pft(NZ)=GrossResp_pft(NZ)-GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
    CanopyRespC_pft(NZ)=CanopyRespC_pft(NZ)-GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
!     Eco_NBP_col=Eco_NBP_col+GY*(TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc))
!     CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)+GZ*(TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc))
    ECO_ER_col=ECO_ER_col-GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
    Eco_AutoR_col=Eco_AutoR_col-GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  ENDIF
  end associate
  end subroutine TotalBiomRemovalByDisturbance
!------------------------------------------------------------------------------------------

  subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,NonstructElmnt2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmntOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmntOffEcosystem(NumPlantChemElms)

  real(r8) :: EHVST21,EHVST22,EHVST23,EHVST24
  real(r8) :: EHVST21h,EHVST22h,EHVST23h,EHVST24h
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
      NonstructElmnt2Litr(NE)=NonstructElmntRemoval(NE)*EHVST21    !non-structural
      LeafElmnt2Litr(NE)=LeafElmntRemoval(NE)*EHVST21   !leaf
      FineNonleafElmnt2Litr(NE)=FineNonleafElmntRemoval(NE)*EHVST22   !fine, non-woody
      WoodyElmnt2Litr(NE)=WoodyElmntRemoval(NE)*EHVST23   !woody
      StandeadElmnt2Litr(NE)=StandeadElmntRemoval(NE)*EHVST24   !standing dead
    ENDDO
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grain)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)=NonstructElmntRemoval(NE)
      LeafElmnt2Litr(NE)=LeafElmntRemoval(NE)
      FineNonleafElmnt2Litr(NE)=FineNonleafElmntRemoval(NE)-WTHTGE(NE)*FracBiomHarvsted(2,iplthvst_finenonleaf,NZ)
      WoodyElmnt2Litr(NE)=WoodyElmntRemoval(NE)
      StandeadElmnt2Litr(NE)=StandeadElmntRemoval(NE)
    ENDDO
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_allabv)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)=NonstructElmntRemoval(NE)*EHVST21
      LeafElmnt2Litr(NE)=LeafElmntRemoval(NE)*EHVST21
      FineNonleafElmnt2Litr(NE)=FineNonleafElmntRemoval(NE)*EHVST22
      WoodyElmnt2Litr(NE)=WoodyElmntRemoval(NE)*EHVST23
      StandeadElmnt2Litr(NE)=StandeadElmntRemoval(NE)*EHVST24
    ENDDO
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
    DO NE=1,NumPlantChemElms
      NonstructElmnt2Litr(NE)=NonstructElmntRemoval(NE)*EHVST21
      LeafElmnt2Litr(NE)=LeafElmntRemoval(NE)*EHVST21
      FineNonleafElmnt2Litr(NE)=FineNonleafElmntRemoval(NE)*EHVST22
      WoodyElmnt2Litr(NE)=WoodyElmntRemoval(NE)*EHVST23
      StandeadElmnt2Litr(NE)=StandeadElmntRemoval(NE)*EHVST24
    ENDDO
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    EHVST21h=1._r8-FracBiomHarvsted(2,iplthvst_leaf,NZ)*0.5_r8
    EHVST22h=1._r8-FracBiomHarvsted(2,iplthvst_finenonleaf,NZ)*0.5_r8
    EHVST23h=1._r8-FracBiomHarvsted(2,iplthvst_woody,NZ)*0.5_r8
    EHVST24h=1._r8-FracBiomHarvsted(2,iplthvst_stdead,NZ)*0.5_r8

    NonstructElmnt2Litr(ielmc)=NonstructElmntRemoval(ielmc)*EHVST21
    LeafElmnt2Litr(ielmc)=LeafElmntRemoval(ielmc)*EHVST21
    FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
    WoodyElmnt2Litr(ielmc)=WoodyElmntRemoval(ielmc)*EHVST23
    StandeadElmnt2Litr(ielmc)=StandeadElmntRemoval(ielmc)*EHVST24

    DO NE=2,NumPlantChemElms
      NonstructElmnt2Litr(NE)=NonstructElmntRemoval(NE)*EHVST21h
      LeafElmnt2Litr(NE)=LeafElmntRemoval(NE)*EHVST21h
      FineNonleafElmnt2Litr(NE)=FineNonleafElmntRemoval(NE)*EHVST22h
      WoodyElmnt2Litr(NE)=WoodyElmntRemoval(NE)*EHVST23h
      StandeadElmnt2Litr(NE)=StandeadElmntRemoval(NE)*EHVST24h
    ENDDO
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
    FERT(17)=FERT(17)+(NonstructElmnt2Litr(ielmc)+LeafElmnt2Litr(ielmc)+&
      FineNonleafElmnt2Litr(ielmc)+WoodyElmnt2Litr(ielmc)+StandeadElmnt2Litr(ielmc))/AREA3(NU)
    FERT(18)=FERT(18)+(NonstructElmnt2Litr(ielmn)+LeafElmnt2Litr(ielmn)+&
      FineNonleafElmnt2Litr(ielmn)+WoodyElmnt2Litr(ielmn)+StandeadElmnt2Litr(ielmn))/AREA3(NU)*0.5_r8
    FERT(3)=FERT(3)+(NonstructElmnt2Litr(ielmn)+LeafElmnt2Litr(ielmn)+&
      FineNonleafElmnt2Litr(ielmn)+WoodyElmnt2Litr(ielmn)+StandeadElmnt2Litr(ielmn))/AREA3(NU)*0.5_r8
    FERT(19)=FERT(19)+(NonstructElmnt2Litr(ielmp)+LeafElmnt2Litr(ielmp)+&
      FineNonleafElmnt2Litr(ielmp)+WoodyElmnt2Litr(ielmp)+StandeadElmnt2Litr(ielmp))/AREA3(NU)
    IYTYP=3
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

    NonstructElmnt2Litr(ielmc)=NonstructElmntRemoval(ielmc)*EHVST21
    NonstructElmnt2Litr(ielmn)=NonstructElmntRemoval(ielmn)*&
      (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
    NonstructElmnt2Litr(ielmp)=NonstructElmntRemoval(ielmp)*&
      (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
    NonstructElmntOffEcosystem(ielmn)=NonstructElmntRemoval(ielmn)*EHVST21
    NonstructElmntOffEcosystem(ielmp)=NonstructElmntRemoval(ielmp)*EHVST21

    LeafElmnt2Litr(ielmc)=LeafElmntRemoval(ielmc)*EHVST21
    LeafElmnt2Litr(ielmn)=LeafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
    LeafElmnt2Litr(ielmp)=LeafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
    LeafElmntOffEcosystem(ielmn)=LeafElmntRemoval(ielmn)*EHVST21
    LeafElmntOffEcosystem(ielmp)=LeafElmntRemoval(ielmp)*EHVST21

    FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
    FineNonleafElmnt2Litr(ielmn)=FineNonleafElmntRemoval(ielmn)*&
      (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_finenonleaf,NZ))
    FineNonleafElmnt2Litr(ielmp)=FineNonleafElmntRemoval(ielmp)*&
      (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_finenonleaf,NZ))
    FineNonleafElmntOffEcosystem(ielmn)=FineNonleafElmntRemoval(ielmn)*EHVST22
    FineNonleafElmntOffEcosystem(ielmp)=FineNonleafElmntRemoval(ielmp)*EHVST22

    WoodyElmnt2Litr(ielmc)=WoodyElmntRemoval(ielmc)*EHVST23
    WoodyElmnt2Litr(ielmn)=WoodyElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_woody,NZ))
    WoodyElmnt2Litr(ielmp)=WoodyElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_woody,NZ))
    WoodyElmntOffEcosystem(ielmn)=WoodyElmntRemoval(ielmn)*EHVST23
    WoodyElmntOffEcosystem(ielmp)=WoodyElmntRemoval(ielmp)*EHVST23

    StandeadElmnt2Litr(ielmc)=StandeadElmntRemoval(ielmc)*EHVST24
    StandeadElmnt2Litr(ielmn)=StandeadElmntRemoval(ielmn)*&
      (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_stdead,NZ))
    StandeadElmnt2Litr(ielmp)=StandeadElmntRemoval(ielmp)*&
      (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_stdead,NZ))
    StandeadElmntOffEcosystem(ielmn)=StandeadElmntRemoval(ielmn)*EHVST24
    StandeadElmntOffEcosystem(ielmp)=StandeadElmntRemoval(ielmp)*EHVST24
  ENDIF
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByTillage(I,J,NZ)
  !REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
  implicit none
  integer , intent(in) :: I,J,NZ
  integer :: L,K,M,N,NR,NE,NB,NTG
  real(r8) :: XHVST,XHVST1
  REAL(R8) :: APSILT
  real(r8) :: FDM,VOLWPX
  real(r8) :: WVPLT
!     begin_execution
  associate(                                                              &
    jHarvst_pft                 => plt_distb%jHarvst_pft,                 &
    iDayPlantHarvest_pft        => plt_distb%iDayPlantHarvest_pft,        &
    iDayPlanting_pft            => plt_distb%iDayPlanting_pft,            &
    iSoilDisturbType_col        => plt_distb%iSoilDisturbType_col,        &
    iYearPlanting_pft           => plt_distb%iYearPlanting_pft,           &
    iYearPlantHarvest_pft       => plt_distb%iYearPlantHarvest_pft,       &
    XCORP                       => plt_distb%XCORP,                       &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,         &
    trcg_rootml_pvr             => plt_rbgc%trcg_rootml_pvr,              &
    trcs_rootml_pvr             => plt_rbgc%trcs_rootml_pvr,              &
    UVOLO                       => plt_ew%UVOLO,                          &
    CanopyWater_pft             => plt_ew%CanopyWater_pft,                &
    VHeatCapCanP_pft            => plt_ew%VHeatCapCanP_pft,               &
    PSICanopy_pft               => plt_ew%PSICanopy_pft,                  &
    PPX_pft                     => plt_site%PPX_pft,                      &
    ShootC4NonstC_brch          => plt_biom%ShootC4NonstC_brch,           &
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr,       &
    RootProteinC_pvr            => plt_biom%RootProteinC_pvr,             &
    PopuRootMycoC_pvr           => plt_biom% PopuRootMycoC_pvr,           &
    RootMycoActiveBiomC_pvr     => plt_biom%RootMycoActiveBiomC_pvr,      &
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch,           &
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch,          &
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch,            &
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch,    &
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch,         &
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch,           &
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch,           &
    CanopyNodulStrutElms_brch   => plt_biom%CanopyNodulStrutElms_brch,    &
    ShootStrutElms_brch         => plt_biom%ShootStrutElms_brch,          &
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch,          &
    LeafChemElmByLayerNode_brch => plt_biom%LeafChemElmByLayerNode_brch,  &
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch,         &
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch,           &
    LeafPetolBiomassC_brch      => plt_biom%LeafPetolBiomassC_brch,       &
    SenecStalkStrutElms_brch    => plt_biom%SenecStalkStrutElms_brch,     &
    StalkBiomassC_brch          => plt_biom%StalkBiomassC_brch,           &
    PetoleProteinCNode_brch     => plt_biom%PetoleProteinCNode_brch,      &
    PetioleElmntNode_brch       => plt_biom%PetioleElmntNode_brch,        &
    RootMyco1stStrutElms_rpvr   => plt_biom%RootMyco1stStrutElms_rpvr,    &
    LeafProteinCNode_brch       => plt_biom%LeafProteinCNode_brch,        &
    InternodeStrutElms_brch     => plt_biom%InternodeStrutElms_brch,      &
    CanopyStalkC_pft            => plt_biom%CanopyStalkC_pft,             &
    Root1stElm_raxs             => plt_biom%Root1stElm_raxs,              &
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft,        &
    CanopyLeafShethC_pft        => plt_biom%CanopyLeafShethC_pft,         &
    RootMyco2ndStrutElms_rpvr   => plt_biom%RootMyco2ndStrutElms_rpvr,    &
    RootNodulNonstElms_pvr      => plt_biom%RootNodulNonstElms_pvr,       &
    RootNodulStrutElms_pvr      => plt_biom%RootNodulStrutElms_pvr,       &
    GrainSeedBiomCMean_brch     => plt_allom%GrainSeedBiomCMean_brch,     &
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr,  &
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr,  &
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &
    FracRootElmAlloc2Litr       => plt_allom%FracRootElmAlloc2Litr,       &
    iPlantBranchState_brch      => plt_pheno%iPlantBranchState_brch,      &
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft,     &
    iPlantState_pft             => plt_pheno%iPlantState_pft,             &
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft,       &
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft,   &
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft,        &
    iPlantRootState_pft         => plt_pheno%iPlantRootState_pft,         &
    iPlantShootState_pft        => plt_pheno%iPlantShootState_pft,        &
    CMassHCO3BundleSheath_node  => plt_photo%CMassHCO3BundleSheath_node,  &
    CMassCO2BundleSheath_node   => plt_photo%CMassCO2BundleSheath_node,   &
    CPOOL3_node                 => plt_photo%CPOOL3_node,                 &
    CPOOL4_node                 => plt_photo%CPOOL4_node,                 &
    MaxNumRootLays              => plt_site%MaxNumRootLays,               &
    PlantPopulation_pft         => plt_site%PlantPopulation_pft,          &
    iYearCurrent                => plt_site%iYearCurrent,                 &
    SolarNoonHour_col           => plt_site%SolarNoonHour_col,            &
    VOLWOU                      => plt_site%VOLWOU,                       &
    NU                          => plt_site%NU,                           &
    k_fine_litr                 => pltpar%k_fine_litr,                    &
    k_woody_litr                => pltpar%k_woody_litr,                   &
    inonstruct                  => pltpar%inonstruct,                     &
    ifoliar                     => pltpar%ifoliar,                        &
    istalk                      => pltpar%istalk,                         &
    iroot                       => pltpar%iroot,                          &
    inonfoliar                  => pltpar%inonfoliar,                     &
    icwood                      => pltpar%icwood,                         &
    RootGasLossDisturb_pft      => plt_bgcr%RootGasLossDisturb_pft,       &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,         &
    RootCO2Autor_pvr            => plt_rbgc%RootCO2Autor_pvr,             &
    RootRespPotent_pvr          => plt_rbgc%RootRespPotent_pvr,           &
    RootCO2EmisPot_pvr          => plt_rbgc%RootCO2EmisPot_pvr,           &
    FracPARRadbyCanopy_pft      => plt_rad%FracPARRadbyCanopy_pft,        &
    Root1stLen_rpvr             => plt_morph%Root1stLen_rpvr,             &
    RootVH2O_pvr                => plt_morph%RootVH2O_pvr,                &
    RootAreaPerPlant_pvr        => plt_morph%RootAreaPerPlant_pvr,        &
    RootPoreVol_pvr             => plt_morph%RootPoreVol_pvr,             &
    RootLenDensPerPlant_pvr     => plt_morph%RootLenDensPerPlant_pvr,     &
    Root1stXNumL_pvr            => plt_morph%Root1stXNumL_pvr,            &
    iPlantNfixType_pft          => plt_morph%iPlantNfixType_pft,          &
    RootLenPerPlant_pvr         => plt_morph%RootLenPerPlant_pvr,         &
    Root2ndXNum_rpvr            => plt_morph%Root2ndXNum_rpvr,            &
    Root2ndLen_pvr              => plt_morph%Root2ndLen_pvr,              &
    Root2ndXNum_pvr             => plt_morph%Root2ndXNum_pvr,             &
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft,          &
    MY                          => plt_morph%MY,                          &
    NumRootAxes_pft             => plt_morph%NumRootAxes_pft,             &
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft,           &
    LeafAreaNode_brch           => plt_morph%LeafAreaNode_brch,           &
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch,           &
    PotentialSeedSites_brch     => plt_morph%PotentialSeedSites_brch,     &
    SeedNumSet_brch             => plt_morph%SeedNumSet_brch,             &
    CanopyLeafArea_lpft         => plt_morph%CanopyLeafArea_lpft          &
  )
!     SolarNoonHour_col=hour of solar noon
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!     iYearCurrent=current year
!     iSoilDisturbType_col=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!     PPX,PP=PFT population per m2,grid cell
!     FracPARRadbyCanopy_pft=fraction of radiation received by each PFT canopy
!     VHeatCapCanP_pft=canopy heat capacity
!
  IF(J.EQ.INT(SolarNoonHour_col) .AND. (iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))) &
    .AND. (I.NE.iDayPlanting_pft(NZ) .OR. iYearCurrent.NE.iYearPlanting_pft(NZ)))THEN

    IF(iSoilDisturbType_col.LE.10 .OR. NZ.NE.1)THEN
      IF(I.GT.iDayPlanting_pft(NZ) .OR. iYearCurrent.GT.iYearPlanting_pft(NZ))THEN
        XHVST=XCORP
        PPX_pft(NZ)=PPX_pft(NZ)*XHVST
        PlantPopulation_pft(NZ)=PlantPopulation_pft(NZ)*XHVST
        FracPARRadbyCanopy_pft(NZ)=FracPARRadbyCanopy_pft(NZ)*XHVST
        VHeatCapCanP_pft(NZ)=VHeatCapCanP_pft(NZ)*XHVST
        CanopyLeafShethC_pft(NZ)=0._r8
        CanopyStalkC_pft(NZ)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
        D8975: DO NB=1,NumOfBranches_pft(NZ)
          IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
            IF(PlantPopulation_pft(NZ).LE.0.0)then
              iPlantBranchState_brch(NB,NZ)=iDead
            endif
!
!     LitrFall FROM BRANCHES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
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
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
            XHVST1=1._r8-XHVST
            D6380: DO M=1,jsken
              LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
                +XHVST1*(ElmAllocmat4Litr(ielmc,inonstruct,M,NZ)*(CanopyNonstElms_brch(ielmc,NB,NZ) &
                +CanopyNodulNonstElms_brch(ielmc,NB,NZ) &
                +ShootC4NonstC_brch(NB,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)) &
                +ElmAllocmat4Litr(ielmc,ifoliar,M,NZ)*(LeafStrutElms_brch(ielmc,NB,NZ)*FracShootStalkElmAlloc2Litr(ielmc,k_fine_litr) &
                +CanopyNodulStrutElms_brch(ielmc,NB,NZ)) &
                +ElmAllocmat4Litr(ielmc,inonfoliar,M,NZ)*(PetoleStrutElms_brch(ielmc,NB,NZ)*FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr) &
                +HuskStrutElms_brch(ielmc,NB,NZ)+EarStrutElms_brch(ielmc,NB,NZ)))

              DO NE=2,NumPlantChemElms
                LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *(ElmAllocmat4Litr(NE,inonstruct,M,NZ)*(CanopyNonstElms_brch(NE,NB,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)&
                  +StalkRsrvElms_brch(NE,NB,NZ)) &
                  +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr) &
                  +CanopyNodulStrutElms_brch(NE,NB,NZ)) &
                  +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr) &
                  +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)))
              ENDDO
            ENDDO D6380

            DO M=1,jsken
              DO NE=1,NumPlantChemElms
                LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *ElmAllocmat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr) &
                  +PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr))

                IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
                  SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+XHVST1 &
                    *ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
                ELSE
                  LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                    *ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
                ENDIF
                LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *ElmAllocmat4Litr(NE,icwood,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

                LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *ElmAllocmat4Litr(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
              ENDDO
            ENDDO
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!

            ShootC4NonstC_brch(NB,NZ)=ShootC4NonstC_brch(NB,NZ)*XHVST
            StalkBiomassC_brch(NB,NZ)=StalkBiomassC_brch(NB,NZ)*XHVST
            DO NE=1,NumPlantChemElms
              CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)*XHVST
              CanopyNodulNonstElms_brch(NE,NB,NZ)=CanopyNodulNonstElms_brch(NE,NB,NZ)*XHVST
              ShootStrutElms_brch(NE,NB,NZ)=ShootStrutElms_brch(NE,NB,NZ)*XHVST
              StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)*XHVST
              HuskStrutElms_brch(NE,NB,NZ)=HuskStrutElms_brch(NE,NB,NZ)*XHVST
              EarStrutElms_brch(NE,NB,NZ)=EarStrutElms_brch(NE,NB,NZ)*XHVST
              GrainStrutElms_brch(NE,NB,NZ)=GrainStrutElms_brch(NE,NB,NZ)*XHVST
              LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)*XHVST
              CanopyNodulStrutElms_brch(NE,NB,NZ)=CanopyNodulStrutElms_brch(NE,NB,NZ)*XHVST
              PetoleStrutElms_brch(NE,NB,NZ)=PetoleStrutElms_brch(NE,NB,NZ)*XHVST
              StalkStrutElms_brch(NE,NB,NZ)=StalkStrutElms_brch(NE,NB,NZ)*XHVST
              SenecStalkStrutElms_brch(NE,NB,NZ)=SenecStalkStrutElms_brch(NE,NB,NZ)*XHVST
            ENDDO

            PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)*XHVST
            SeedNumSet_brch(NB,NZ)=SeedNumSet_brch(NB,NZ)*XHVST
            GrainSeedBiomCMean_brch(NB,NZ)=GrainSeedBiomCMean_brch(NB,NZ)*XHVST
            LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)*XHVST
            LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))
            CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)

            CanopyStalkC_pft(NZ)=CanopyStalkC_pft(NZ)+StalkBiomassC_brch(NB,NZ)
            D8970: DO K=0,MaxNodesPerBranch1
              IF(K.NE.0)THEN
                CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)*XHVST
                CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)*XHVST
                CMassCO2BundleSheath_node(K,NB,NZ)=CMassCO2BundleSheath_node(K,NB,NZ)*XHVST
                CMassHCO3BundleSheath_node(K,NB,NZ)=CMassHCO3BundleSheath_node(K,NB,NZ)*XHVST
              ENDIF
              LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)*XHVST

              LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*XHVST
!     PetoleLensNode_brch(K,NB,NZ)=PetoleLensNode_brch(K,NB,NZ)*XHVST

              PetoleProteinCNode_brch(K,NB,NZ)=PetoleProteinCNode_brch(K,NB,NZ)*XHVST
!     LiveInterNodeHight_brch(K,NB,NZ)=LiveInterNodeHight_brch(K,NB,NZ)*XHVST
!     InternodeHeightDying_brch(K,NB,NZ)=InternodeHeightDying_brch(K,NB,NZ)*XHVST
              DO NE=1,NumPlantChemElms
                InternodeStrutElms_brch(NE,K,NB,NZ)=InternodeStrutElms_brch(NE,K,NB,NZ)*XHVST
                LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)*XHVST
                PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)*XHVST
                DO L=1,NumOfCanopyLayers1
                  LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)=LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)*XHVST
                ENDDO
              ENDDO
              D8965: DO L=1,NumOfCanopyLayers1
                CanopyLeafArea_lpft(L,K,NB,NZ)=CanopyLeafArea_lpft(L,K,NB,NZ)*XHVST
              ENDDO D8965
            ENDDO D8970
          ENDIF
        ENDDO D8975
!
!     PSICanopy_pft=canopy water potential
!     CanopyWater_pft=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        VOLWPX=CanopyWater_pft(NZ)
        WVPLT=AZMAX1(CanopyLeafShethC_pft(NZ)+CanopyStalkC_pft(NZ))

        FDM=get_FDM(PSICanopy_pft(NZ))
!        APSILT=ABS(PSICanopy_pft(NZ))
!        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)

        CanopyWater_pft(NZ)=ppmc*WVPLT/FDM
        VOLWOU=VOLWOU+VOLWPX-CanopyWater_pft(NZ)
        UVOLO=UVOLO+VOLWPX-CanopyWater_pft(NZ)
!
!     TERMINATE ROOTS IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     PP=PFT population
!     IDTHR,iPlantShootState_pft=PFT root,shoot living flag: 0=alive,1=dead
!     IDTH=PFT living flag: 0=alive,1=dead
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,and reseed
!     iDayPlantHarvest_pft,iYearPlantHarvest_pft=day,year of harvesting
!     iYearCurrent=current year
!
        IF(PlantPopulation_pft(NZ).LE.0.0_r8)THEN
          iPlantRootState_pft(NZ)=iDead
          iPlantShootState_pft(NZ)=iDead
          iPlantState_pft(NZ)=iDead
          jHarvst_pft(NZ)=jharvtyp_terminate
          iDayPlantHarvest_pft(NZ)=I
          iYearPlantHarvest_pft(NZ)=iYearCurrent
        ENDIF
!
!     LitrFall FROM ROOTS DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!
        D8980: DO L=NU,MaxNumRootLays
          D8985: DO N=1,MY(NZ)

            D6385: DO M=1,jsken
                DO NE=1,NumPlantChemElms
                  LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *ElmAllocmat4Litr(NE,inonstruct,M,NZ)* RootMycoNonstElms_rpvr(NE,N,L,NZ)
                ENDDO

              DO NR=1,NumRootAxes_pft(NZ)
                DO NE=1,NumPlantChemElms
                  LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)+XHVST1 &
                    *ElmAllocmat4Litr(NE,icwood,M,NZ)*(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
                    +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_woody_litr)

                  LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *ElmAllocmat4Litr(NE,iroot,M,NZ)*(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
                    +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_fine_litr)
                ENDDO
              ENDDO
            ENDDO D6385
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!

            DO NTG=idg_beg,idg_end-1
              RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-XHVST1 &
                *(trcg_rootml_pvr(NTG,N,L,NZ)+trcs_rootml_pvr(NTG,N,L,NZ))
              trcg_rootml_pvr(NTG,N,L,NZ)=XHVST*trcg_rootml_pvr(NTG,N,L,NZ)
              trcs_rootml_pvr(NTG,N,L,NZ)=XHVST*trcs_rootml_pvr(NTG,N,L,NZ)
            ENDDO
!
!     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stLen_rpvr,Root2ndLen_pvr=primary,secondary root length
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
            D8960: DO NR=1,NumRootAxes_pft(NZ)
              DO NE=1,NumPlantChemElms
                RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)*XHVST
                RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*XHVST
                Root1stElm_raxs(NE,N,NR,NZ)=Root1stElm_raxs(NE,N,NR,NZ)*XHVST
              ENDDO
              Root1stLen_rpvr(N,L,NR,NZ)=Root1stLen_rpvr(N,L,NR,NZ)*XHVST
              Root2ndLen_pvr(N,L,NR,NZ)=Root2ndLen_pvr(N,L,NR,NZ)*XHVST
              Root2ndXNum_rpvr(N,L,NR,NZ)=Root2ndXNum_rpvr(N,L,NR,NZ)*XHVST
            ENDDO D8960
            DO NE=1,NumPlantChemElms
               RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)*XHVST
            ENDDO
            RootMycoActiveBiomC_pvr(N,L,NZ)=RootMycoActiveBiomC_pvr(N,L,NZ)*XHVST
             PopuRootMycoC_pvr(N,L,NZ)= PopuRootMycoC_pvr(N,L,NZ)*XHVST
            RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)*XHVST
            Root1stXNumL_pvr(N,L,NZ)=Root1stXNumL_pvr(N,L,NZ)*XHVST
            Root2ndXNum_pvr(N,L,NZ)=Root2ndXNum_pvr(N,L,NZ)*XHVST
            RootLenPerPlant_pvr(N,L,NZ)=RootLenPerPlant_pvr(N,L,NZ)*XHVST
            RootLenDensPerPlant_pvr(N,L,NZ)=RootLenDensPerPlant_pvr(N,L,NZ)*XHVST
            RootPoreVol_pvr(N,L,NZ)=RootPoreVol_pvr(N,L,NZ)*XHVST
            RootVH2O_pvr(N,L,NZ)=RootVH2O_pvr(N,L,NZ)*XHVST
            RootAreaPerPlant_pvr(N,L,NZ)=RootAreaPerPlant_pvr(N,L,NZ)*XHVST
            RootRespPotent_pvr(N,L,NZ)=RootRespPotent_pvr(N,L,NZ)*XHVST
            RootCO2EmisPot_pvr(N,L,NZ)=RootCO2EmisPot_pvr(N,L,NZ)*XHVST
            RootCO2Autor_pvr(N,L,NZ)=RootCO2Autor_pvr(N,L,NZ)*XHVST
!
!     LitrFall AND STATE VARIABLES FOR NODULES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(is_plant_N2fix(iPlantNfixType_pft(NZ)).AND.N.EQ.ipltroot)THEN
              DO NE=1,NumPlantChemElms
                D6395: DO M=1,jsken
                  LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+&
                    XHVST1*(ElmAllocmat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_pvr(NE,L,NZ) &
                    +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_pvr(NE,L,NZ))
                ENDDO D6395
                RootNodulStrutElms_pvr(NE,L,NZ)=RootNodulStrutElms_pvr(NE,L,NZ)*XHVST
                RootNodulNonstElms_pvr(NE,L,NZ)=RootNodulNonstElms_pvr(NE,L,NZ)*XHVST
              ENDDO
            ENDIF
          ENDDO D8985
        ENDDO D8980
!
!     LitrFall AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
!     DURING TILLAGE
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        DO NE=1,NumPlantChemElms
          D6400: DO M=1,jsken
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=&
              LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
              +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ))*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)=&
              LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
              +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ))*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
          ENDDO D6400
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)*XHVST
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------

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
  real(r8) :: WHVSNX
  integer :: NTG
!     begin_execution
  associate(                                                                          &
    CutHeightORFrac_pft               => plt_distb%CutHeightORFrac_pft,               &
    THIN_pft                          => plt_distb%THIN_pft,                          &
    iHarvstType_pft                   => plt_distb%iHarvstType_pft,                   &
    jHarvst_pft                       => plt_distb%jHarvst_pft,                       &
    PlantPopulation_pft               => plt_site%PlantPopulation_pft,                &
    PPI_pft                           => plt_site%PPI_pft,                            &
    PPX_pft                           => plt_site%PPX_pft,                            &
    NU                                => plt_site%NU,                                 &
    MaxNumRootLays                    => plt_site%MaxNumRootLays,                     &
    SolarNoonHour_col                 => plt_site%SolarNoonHour_col,                  &
    ZEROS                             => plt_site%ZEROS,                              &
    AREA3                             => plt_site%AREA3,                              &
    ShootC4NonstC_brch                => plt_biom%ShootC4NonstC_brch,                 &
    SeasonalNonstElms_pft             => plt_biom%SeasonalNonstElms_pft,              &
    CanopyStalkC_pft                  => plt_biom%CanopyStalkC_pft,                   &
    CanopyLeafShethC_pft              => plt_biom%CanopyLeafShethC_pft,               &
    StalkBiomassC_brch                => plt_biom%StalkBiomassC_brch,                 &
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch,                &
    ShootStrutElms_brch               => plt_biom%ShootStrutElms_brch,                &
    LeafPetolBiomassC_brch            => plt_biom%LeafPetolBiomassC_brch,             &
    LeafChemElmByLayerNode_brch       => plt_biom%LeafChemElmByLayerNode_brch,        &
    StalkStrutElms_pft                => plt_biom%StalkStrutElms_pft,                 &
    FracRootStalkElmAlloc2Litr        => plt_allom%FracRootStalkElmAlloc2Litr,        &
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft,           &
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft,              &
    ElmAllocmat4Litr                  => plt_soilchem%ElmAllocmat4Litr,               &
    inonstruct                        => pltpar%inonstruct,                           &
    ifoliar                           => pltpar%ifoliar,                              &
    istalk                            => pltpar%istalk,                               &
    inonfoliar                        => pltpar%inonfoliar,                           &
    k_fine_litr                       => pltpar%k_fine_litr,                          &
    k_woody_litr                      => pltpar%k_woody_litr,                         &
    LitrfalStrutElms_pvr              => plt_bgcr%LitrfalStrutElms_pvr,               &
    NGTopRootLayer_pft                => plt_morph%NGTopRootLayer_pft,                &
    MY                                => plt_morph%MY,                                &
    CanopyLeafAareZ_col               => plt_morph%CanopyLeafAareZ_col,               &
    CanopyHeightZ_col                 => plt_morph%CanopyHeightZ_col,                 &
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft,                 &
    CanopyStemArea_pft                => plt_morph%CanopyStemArea_pft,                &
    LiveInterNodeHight_brch           => plt_morph%LiveInterNodeHight_brch,           &
    CanopyStalkArea_lbrch             => plt_morph%CanopyStalkArea_lbrch,             &
    ClumpFactor_pft                   => plt_morph%ClumpFactor_pft,                   &
    CanopyLeafArea_col                => plt_morph%CanopyLeafArea_col                 &
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
        !terminate and reseed
        PPX_pft(NZ)=PPX_pft(NZ)*(1._r8-THIN_pft(NZ))
        PlantPopulation_pft(NZ)=PlantPopulation_pft(NZ)*(1._r8-THIN_pft(NZ))
      ELSE
!     PPI_pft(NZ)=AMAX1(1.0_r8,0.5_r8*(PPI_pft(NZ)+CanopySeedNum_pft(NZ)/AREA3(NU)))
        PPX_pft(NZ)=PPI_pft(NZ)
        PlantPopulation_pft(NZ)=PPX_pft(NZ)*AREA3(NU)
      ENDIF
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
        ClumpFactor_pft(NZ)=ClumpFactor_pft(NZ)*CutHeightORFrac_pft(NZ)
      ENDIF
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.CutHeightORFrac_pft(NZ).LT.0.0)THEN
        ARLFY=(1._r8-ABS(CutHeightORFrac_pft(NZ)))*CanopyLeafArea_col
        ARLFR=0._r8
        D9875: DO L=1,NumOfCanopyLayers1
          IF(CanopyHeightZ_col(L).GT.CanopyHeightZ_col(L-1) &
            .AND. CanopyLeafAareZ_col(L).GT.ZEROS .AND. ARLFR.LT.ARLFY)THEN
            IF(ARLFR+CanopyLeafAareZ_col(L).GT.ARLFY)THEN
              CutHeightORFrac_pft(NZ)=CanopyHeightZ_col(L-1)+((ARLFY-ARLFR)/CanopyLeafAareZ_col(L))&
                *(CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1))
            ENDIF
          ELSE
            CutHeightORFrac_pft(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+CanopyLeafAareZ_col(L)
        ENDDO D9875
      ENDIF
      HvstedLeafC=0._r8
      HvstedShethC=0._r8
      HvstedEarC=0._r8
      HvstedGrainC=0._r8
      WHVSCP=0._r8
      HvstedStalkC=0._r8
      HvstedRsrvC=0._r8
      LeafC_lbrch=0._r8          !it is a filler 
    ELSE
!
!     GRAZING REMOVAL
      call GrazingPlant(I,J,NZ,HvstedLeafC,HvstedShethC,HvstedEarC,HvstedGrainC,&
        WHVSCP,HvstedStalkC,HvstedRsrvC,WHVSHH,WHVSNP)

!     ALL HARVEST REMOVALS
!
!     LeafC_lbrch=branch leaf C mass in canopy layer
!
      D9860: DO NB=1,NumOfBranches_pft(NZ)
        DO  L=1,NumOfCanopyLayers1
          DO  K=0,MaxNodesPerBranch1
            LeafC_lbrch(L,NB,NZ)=0._r8
          enddo
        enddo
      ENDDO D9860

      D9870: DO NB=1,NumOfBranches_pft(NZ)
        DO  L=1,NumOfCanopyLayers1
          DO  K=0,MaxNodesPerBranch1
            LeafC_lbrch(L,NB,NZ)=LeafC_lbrch(L,NB,NZ)+LeafChemElmByLayerNode_brch(ielmc,L,K,NB,NZ)
          enddo
        enddo
      ENDDO D9870
    ENDIF
!
!     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
    call HarvestCanopy(I,J,NZ,HvstedLeafC,LeafC_lbrch)

    CALL CutPlant(I,J,NZ,WHVSHH,WHVSCP,WHVSNP,HvstedShethC,HvstedGrainC,HvstedEarC,HvstedRsrvC,HvstedStalkC)

    CanopyLeafShethC_pft(NZ)=0._r8
    StalkStrutElms_pft(ielmc,NZ)=0._r8
    CanopyStalkC_pft(NZ)=0._r8
    CanopyStemArea_pft(NZ)=0._r8
    D9840: DO NB=1,NumOfBranches_pft(NZ)
      CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)
      StalkStrutElms_pft(ielmc,NZ)=StalkStrutElms_pft(ielmc,NZ)+StalkStrutElms_brch(ielmc,NB,NZ)
      CanopyStalkC_pft(NZ)=CanopyStalkC_pft(NZ)+StalkBiomassC_brch(NB,NZ)
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
!     CO2ByFire_pft,CH4ByFire_pft,O2ByFire_pft,NH3byFire_pft,N2ObyFire_pft,PO4byFire_pft=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_col=total net biome productivity
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      FracLeftThin=1.0_r8-THIN_pft(NZ)

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
        DO NE=1,NumPlantChemElms
          D3400: DO M=1,jsken
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=&
                LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
              +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ))&
              *FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)=&
                LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
              +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ))&
              *FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
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
    MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
    NodeNum2InitFloral_brch(NB,NZ)=ShootNodeNum_brch(NB,NZ)
    NodeNumberAtAnthesis_brch(NB,NZ)=0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)=0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)=0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
    HourFailGrainFill_brch(NB,NZ)=0._r8
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)=I
    D3005: DO M=2,NumGrowthStages
      iPlantCalendar_brch(M,NB,NZ)=0
    ENDDO D3005
    doInitLeafOut_brch(NB,NZ)=iEnable

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      D3010: DO NBX=1,NumOfBranches_pft(NZ)
        IF(NBX.NE.MainBranchNum_pft(NZ))THEN
          MatureGroup_brch(NBX,NZ)=MatureGroup_pft(NZ)
          NodeNum2InitFloral_brch(NBX,NZ)=ShootNodeNum_brch(NBX,NZ)
          NodeNumberAtAnthesis_brch(NBX,NZ)=0._r8
          LeafNumberAtFloralInit_brch(NBX,NZ)=0._r8
          TotalNodeNumNormByMatgrp_brch(NBX,NZ)=0._r8
          TotReproNodeNumNormByMatrgrp_brch(NBX,NZ)=0._r8
          HourFailGrainFill_brch(NBX,NZ)=0._r8
          iPlantCalendar_brch(ipltcal_Emerge,NBX,NZ)=I
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
    CutHeightORFrac_pft       => plt_distb%CutHeightORFrac_pft,       &
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
        FHGT=AZMAX1(AMIN1(1.0_r8,CutHeightORFrac_pft(NZ)/RMedInternodeLen))
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
!     WoodyElmntRemoval(ielmc),WoodyElmntRemoval(ielmn),WoodyElmntRemoval(ielmp)=harvested stalk C,N,P
!     WoodyElmntHarv2Litr(ielmc),WoodyElmntHarv2Litr(ielmn),WoodyElmntHarv2Litr(ielmp)=harvested stalk C,N,P to litter
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
          FHGTK=AZMAX1(AMIN1(1.0_r8,(LiveInterNodeHight_brch(K,NB,NZ)-CutHeightORFrac_pft(NZ))/&
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
      LiveInterNodeHight_brch(K,NB,NZ)=AMIN1(LiveInterNodeHight_brch(K,NB,NZ),CutHeightORFrac_pft(NZ))
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
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(StalkStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=FrcLeafMassLeft
      FHVSH=FHVSH
    ELSE
      FrcLeafMassLeft=0._r8
      FHVSH=0._r8
    ENDIF
  ELSE
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
!     WoodyElmntRemoval(ielmc),WoodyElmntRemoval(ielmn),WoodyElmntRemoval(ielmp)=harvested stalk C,N,P
!     WoodyElmntHarv2Litr(ielmc),WoodyElmntHarv2Litr(ielmn),WoodyElmntHarv2Litr(ielmp)=harvested stalk C,N,P to litter
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
  real(r8) :: FHVSETG,FHVSHG,FHVSETH,FHVSETE,FHVSHH,FHVSHE
  associate(                                                      &
    CutHeightORFrac_pft     => plt_distb%CutHeightORFrac_pft,     &
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
!     FHVSETG,FHVSETH,FHVSETE=fraction of grain,husk,ear mass not harvested
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(CutHeightORFrac_pft(NZ).LT.RMedInternodeLen .OR. iHarvstType_pft(NZ).EQ.iharvtyp_grain &
      .OR. iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FHVSETG=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)
        FHVSHG=FHVSETG
      ELSE
        FHVSETG=1.0_r8-THIN_pft(NZ)
        FHVSHG=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
      ENDIF
    ELSE
      FHVSETG=1.0_r8-THIN_pft(NZ)
      FHVSHG=FHVSETG
    ENDIF
    FHVSETH=FHVSETG
    FHVSETE=FHVSETG
    FHVSHH=FHVSHG
    FHVSHE=FHVSHG
  ELSE
    IF(HuskStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FHVSETH=AZMAX1(AMIN1(1.0_r8,1._r8-HvstedShethC/HuskStrutElms_pft(ielmc,NZ)))
      FHVSHH=FHVSETH
    ELSE
      FHVSETH=1.0_r8
      FHVSHH=1.0_r8
    ENDIF
    IF(EarStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FHVSETE=AZMAX1(AMIN1(1.0_r8,1._r8-HvstedEarC/EarStrutElms_pft(ielmc,NZ)))
      FHVSHE=FHVSETE
    ELSE
      FHVSETE=1.0_r8
      FHVSHE=1.0_r8
    ENDIF
    IF(GrainStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FHVSETG=AZMAX1(AMIN1(1.0_r8,1._r8-HvstedGrainC/GrainStrutElms_pft(ielmc,NZ)))
      FHVSHG=FHVSETG
    ELSE
      FHVSETG=1.0_r8
      FHVSHG=1.0_r8
    ENDIF
  ENDIF
!
!     HARVESTED REPRODUCTIVE C,N,P
!
!     FineNonleafElmntRemoval(ielmc),FineNonleafElmntRemoval(ielmn),FineNonleafElmntRemoval(ielmp)=reproductive C,N,P removed
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTHTGE(ielmc),WTHTGE(ielmn),WTHTGE(ielmp)=grain harvested
!
  DO NE=1,NumPlantChemElms
    FineNonleafElmntRemoval(NE)=FineNonleafElmntRemoval(NE)+(1._r8-FHVSHH)*HuskStrutElms_brch(NE,NB,NZ)&
      +(1._r8-FHVSHE)*EarStrutElms_brch(NE,NB,NZ)+(1._r8-FHVSHG)*GrainStrutElms_brch(NE,NB,NZ)
    PetioleElmntHarv2Litr(NE)=PetioleElmntHarv2Litr(NE)+(FHVSHH-FHVSETH)*HuskStrutElms_brch(NE,NB,NZ) &
      +(FHVSHE-FHVSETE)*EarStrutElms_brch(NE,NB,NZ)+(FHVSHG-FHVSETG)*GrainStrutElms_brch(NE,NB,NZ)
    WTHTGE(NE)=WTHTGE(NE)+(1._r8-FHVSETG)*GrainStrutElms_brch(NE,NB,NZ)

!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
    HuskStrutElms_brch(NE,NB,NZ)=FHVSETH*HuskStrutElms_brch(NE,NB,NZ)
    EarStrutElms_brch(NE,NB,NZ)=FHVSETE*EarStrutElms_brch(NE,NB,NZ)
    GrainStrutElms_brch(NE,NB,NZ)=FHVSETG*GrainStrutElms_brch(NE,NB,NZ)
  ENDDO
  PotentialSeedSites_brch(NB,NZ)=FHVSETG*PotentialSeedSites_brch(NB,NZ)
  SeedNumSet_brch(NB,NZ)=FHVSETG*SeedNumSet_brch(NB,NZ)
  GrainSeedBiomCMean_brch(NB,NZ)=FHVSETG*GrainSeedBiomCMean_brch(NB,NZ)
  end associate
  END subroutine BranchCutReprodOrgans

!--------------------------------------------------------------------------------

  subroutine BranchCutNonstructural(I,J,NB,NZ,WGLFGX,WGSHGX,WGLFGY,WGSHGY,WHVSCP,WHVSNP)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: WGLFGX,WGSHGX,WGLFGY,WGSHGY
  real(r8), intent(in) :: WHVSNP,WHVSCP
  real(r8) :: WHVSNX
  real(r8) :: CPOOLX,ZPOOLX,PPOOLX  
  real(r8) :: CPOLNX,ZPOLNX,PPOLNX
  real(r8) :: CPOOLG,ZPOOLG
  real(r8) :: PPOOLG
  real(r8) :: CPOLNG  
  real(r8) :: ZPOLNG  
  real(r8) :: PPOLNG,WTNDG,WTNDNG  
  real(r8) :: WTNDPG,WHVSCX
  real(r8) :: WTLSBX  
  real(r8) :: FHVST4,dFHVST4  
  integer  :: K,NE
  real(r8) :: FrcLeafMassLeft
  
  associate(                                                            &
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch,        &
    CanopyNodulStrutElms_brch  => plt_biom%CanopyNodulStrutElms_brch,   &
    LeafPetolBiomassC_brch     => plt_biom%LeafPetolBiomassC_brch,      &
    CanopyLeafShethC_pft       => plt_biom%CanopyLeafShethC_pft,        &
    CPOOL3_node                => plt_photo%CPOOL3_node,                &
    CPOOL4_node                => plt_photo%CPOOL4_node,                &
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node, &
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node,  &
    ZERO4LeafVar_pft           => plt_biom%ZERO4LeafVar_pft,            &
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
!     CPOOLG,ZPOOLG,PPOOLG=branch non-structural C,N,P mass after harvest
!     CPOLNG,ZPOLNG,PPOLNG=nonstructural C,N,P in bacteria after harvest
!     WTNDG,WTNDNG,WTNDPG=bacterial C,N,P mass after harvest
!     WTLS,LeafPetolBiomassC_brch=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
  CPOOLX=AZMAX1(CanopyNonstElms_brch(ielmc,NB,NZ))
  ZPOOLX=AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
  PPOOLX=AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
  CPOLNX=AZMAX1(CanopyNodulNonstElms_brch(ielmc,NB,NZ))
  ZPOLNX=AZMAX1(CanopyNodulNonstElms_brch(ielmn,NB,NZ))
  PPOLNX=AZMAX1(CanopyNodulNonstElms_brch(ielmp,NB,NZ))
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(WGLFGY+WGSHGY.GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
      CPOOLG=CPOOLX*FrcLeafMassLeft
      ZPOOLG=ZPOOLX*FrcLeafMassLeft
      PPOOLG=PPOOLX*FrcLeafMassLeft
      CPOLNG=CPOLNX*FrcLeafMassLeft
      ZPOLNG=ZPOLNX*FrcLeafMassLeft
      PPOLNG=PPOLNX*FrcLeafMassLeft
      
      WTNDG=CanopyNodulStrutElms_brch(ielmc,NB,NZ)*FrcLeafMassLeft
      WTNDNG=CanopyNodulStrutElms_brch(ielmn,NB,NZ)*FrcLeafMassLeft
      WTNDPG=CanopyNodulStrutElms_brch(ielmp,NB,NZ)*FrcLeafMassLeft
    ELSE
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
    ENDIF
  ELSE
    IF(CanopyLeafShethC_pft(NZ).GT.ZERO4LeafVar_pft(NZ))THEN
      WTLSBX=AZMAX1(LeafPetolBiomassC_brch(NB,NZ))
      IF(CanopyNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        WHVSCX=AZMAX1(WHVSCP)*WTLSBX/CanopyLeafShethC_pft(NZ)
        CPOOLG=AZMAX1(CPOOLX-WHVSCX)
        ZPOOLG=AZMAX1(ZPOOLX-WHVSCX*ZPOOLX/CanopyNonstElms_brch(ielmc,NB,NZ))
        PPOOLG=AZMAX1(PPOOLX-WHVSCX*PPOOLX/CanopyNonstElms_brch(ielmc,NB,NZ))
      ELSE
        CPOOLG=0._r8
        ZPOOLG=0._r8
        PPOOLG=0._r8
      ENDIF
      IF(CanopyNodulNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        WHVSNX=AZMAX1(WHVSNP)*WTLSBX/CanopyLeafShethC_pft(NZ)
        CPOLNG=AZMAX1(CPOLNX-WHVSNX)
        ZPOLNG=AZMAX1(ZPOLNX-WHVSNX*ZPOLNX/CanopyNodulNonstElms_brch(ielmc,NB,NZ))
        PPOLNG=AZMAX1(PPOLNX-WHVSNX*PPOLNX/CanopyNodulNonstElms_brch(ielmc,NB,NZ))

        WTNDG=CanopyNodulStrutElms_brch(ielmc,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
        WTNDNG=CanopyNodulStrutElms_brch(ielmn,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
        WTNDPG=CanopyNodulStrutElms_brch(ielmp,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
      ELSE
        CPOLNG=0._r8
        ZPOLNG=0._r8
        PPOLNG=0._r8

        WTNDG=0._r8
        WTNDNG=0._r8
        WTNDPG=0._r8
      ENDIF
    ELSE
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
    ENDIF
  ENDIF
!
!     HARVESTED NON-STRUCTURAL C, N, P
!
!     NonstructElmntRemoval(ielmc),NonstructElmntRemoval(ielmn),NonstructElmntRemoval(ielmp)=nonstructural C,N,P removed
!
  NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+CPOOLX-CPOOLG+CPOLNX-CPOLNG
  NonstructElmntRemoval(ielmn)=NonstructElmntRemoval(ielmn)+ZPOOLX-ZPOOLG+ZPOLNX-ZPOLNG
  NonstructElmntRemoval(ielmp)=NonstructElmntRemoval(ielmp)+PPOOLX-PPOOLG+PPOLNX-PPOLNG
  NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+CanopyNodulStrutElms_brch(ielmc,NB,NZ)-WTNDG
  NonstructElmntRemoval(ielmn)=NonstructElmntRemoval(ielmn)+CanopyNodulStrutElms_brch(ielmn,NB,NZ)-WTNDNG
  NonstructElmntRemoval(ielmp)=NonstructElmntRemoval(ielmp)+CanopyNodulStrutElms_brch(ielmp,NB,NZ)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  CanopyNonstElms_brch(ielmc,NB,NZ)=CPOOLG
  CanopyNonstElms_brch(ielmn,NB,NZ)=ZPOOLG
  CanopyNonstElms_brch(ielmp,NB,NZ)=PPOOLG
  CanopyNodulNonstElms_brch(ielmc,NB,NZ)=CPOLNG
  CanopyNodulNonstElms_brch(ielmn,NB,NZ)=ZPOLNG
  CanopyNodulNonstElms_brch(ielmp,NB,NZ)=PPOLNG
  CanopyNodulStrutElms_brch(ielmc,NB,NZ)=WTNDG
  CanopyNodulStrutElms_brch(ielmn,NB,NZ)=WTNDNG
  CanopyNodulStrutElms_brch(ielmp,NB,NZ)=WTNDPG
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CPOOLG=branch non-structural C mass after harvest
!     NonstructElmntRemoval(ielmc),NonstructElmntRemoval(ielmn),NonstructElmntRemoval(ielmp)=nonstructural C,N,P removed
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!
  IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo.AND.CPOOLX.GT.ZERO4Groth_pft(NZ))THEN
    FHVST4=CPOOLG/CPOOLX
    dFHVST4=1._r8-FHVST4
    D9810: DO K=1,MaxNodesPerBranch1
      NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc) &
        +dFHVST4*(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ))
      CPOOL3_node(K,NB,NZ)=FHVST4*CPOOL3_node(K,NB,NZ)
      CPOOL4_node(K,NB,NZ)=FHVST4*CPOOL4_node(K,NB,NZ)
      CMassCO2BundleSheath_node(K,NB,NZ)=FHVST4*CMassCO2BundleSheath_node(K,NB,NZ)
      CMassHCO3BundleSheath_node(K,NB,NZ)=FHVST4*CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D9810
  ENDIF
  end associate
  end subroutine BranchCutNonstructural
!--------------------------------------------------------------------------------
  subroutine BranchCutSheathPetole(I,J,NB,NZ,WHVSHH,FHVSETK,FHVSHK,RMedInternodeLen,WGSHGX,WGSHGY)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: WHVSHH
  real(r8), intent(inout) :: FHVSETK(0:MaxNodesPerBranch1)
  real(r8), intent(inout) :: FHVSHK(0:MaxNodesPerBranch1)  
  real(r8), intent(out) :: RMedInternodeLen
  real(r8), intent(out) :: WGSHGX,WGSHGY
  real(r8) :: WHVSBS,FHGT
  integer :: K,NE

  associate(                                                            &
    PetoleStrutElms_pft       => plt_biom%PetoleStrutElms_pft,        &
    PetoleStrutElms_brch       => plt_biom%PetoleStrutElms_brch,        &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft,              &
    PetioleElmntNode_brch      => plt_biom%PetioleElmntNode_brch,       &
    CutHeightORFrac_pft        => plt_distb%CutHeightORFrac_pft,        &
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
!     WHVSBS,WHVSHH=branch, PFT petiole C mass removed
!     LiveInterNodeHight_brch=internode length
!     RMedInternodeLen=internode length removed
!
  RMedInternodeLen=0._r8
  IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
    .AND. PetoleStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    WHVSBS=WHVSHH*PetoleStrutElms_brch(ielmc,NB,NZ)/PetoleStrutElms_pft(ielmc,NZ)
  ELSE
    WHVSBS=0._r8
  ENDIF

  D9805: DO K=MaxNodesPerBranch1,0,-1
    IF(LiveInterNodeHight_brch(K,NB,NZ).GT.0.0) RMedInternodeLen=AMAX1(RMedInternodeLen,LiveInterNodeHight_brch(K,NB,NZ))
!
!     HARVESTED SHEATH OR PETIOLE C,N,P
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WHVSBS=branch petiole C mass removed
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinCNode_brch=node petiole C,N,P,protein mass
!     FHVSETK=fraction of internode layer mass not harvested
!     FineNonleafElmntRemoval(ielmc),FineNonleafElmntRemoval(ielmn),FineNonleafElmntRemoval(ielmp)=harvested petiole C,N,P
!     PetioleElmntHarv2Litr(ielmc),PetioleElmntHarv2Litr(ielmn),PetioleElmntHarv2Litr(ielmp)=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     PetoleLensNode_brch,LiveInterNodeHight_brch=petiole,internode length
!
    IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) .OR. WHVSBS.GT.0.0_r8)THEN

      IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
        IF(PetioleElmntNode_brch(ielmc,K,NB,NZ).GT.WHVSBS)THEN
          FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(PetioleElmntNode_brch(ielmc,K,NB,NZ)-WHVSBS)/PetioleElmntNode_brch(ielmc,K,NB,NZ)))
          FHVSHK(K)=FHVSETK(K)
        ELSE
          FHVSETK(K)=0._r8
          FHVSHK(K)=0._r8
        ENDIF
      ENDIF
      WHVSBS=WHVSBS-(1._r8-FHVSETK(K))*PetioleElmntNode_brch(ielmc,K,NB,NZ)
      DO NE=1,NumPlantChemElms
        FineNonleafElmntRemoval(NE)=FineNonleafElmntRemoval(NE) &
          +(1._r8-FHVSHK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        PetioleElmntHarv2Litr(NE)=PetioleElmntHarv2Litr(NE) &
          +(FHVSHK(K)-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE) &
          +(1._r8-FHVSHK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
        WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE) &
          +(FHVSHK(K)-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
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
      PetoleProteinCNode_brch(K,NB,NZ)=FHVSETK(K)*PetoleProteinCNode_brch(K,NB,NZ)

      DO NE=1,NumPlantChemElms
        PetoleStrutElms_brch(NE,NB,NZ)=PetoleStrutElms_brch(NE,NB,NZ) &
          -(1._r8-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)
        PetioleElmntNode_brch(NE,K,NB,NZ)=FHVSETK(K)*PetioleElmntNode_brch(NE,K,NB,NZ)
      ENDDO
!            PetoleProteinCNode_brch(K,NB,NZ)=FHVSETK(K)*PetoleProteinCNode_brch(K,NB,NZ)
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.PetoleLensNode_brch(K,NB,NZ).GT.0.0_r8)THEN
        FHGT=AZMAX1(AMIN1(1.0_r8,(LiveInterNodeHight_brch(K,NB,NZ) &
          +PetoleLensNode_brch(K,NB,NZ)-CutHeightORFrac_pft(NZ))/PetoleLensNode_brch(K,NB,NZ)))
        PetoleLensNode_brch(K,NB,NZ)=(1._r8-FHGT)*PetoleLensNode_brch(K,NB,NZ)
      ELSE
        PetoleLensNode_brch(K,NB,NZ)=FHVSETK(K)*PetoleLensNode_brch(K,NB,NZ)
      ENDIF
      WGSHGX=WGSHGX+PetioleElmntNode_brch(ielmc,K,NB,NZ)

!     IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
!     IF(LiveInterNodeHight_brch(K,NB,NZ).GT.CutHeightORFrac_pft(NZ)
!    2.OR.iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
!     IF(isclose(FHVSETK(K),0._r8).AND.K.GT.0)THEN
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
  end subroutine BranchCutSheathPetole
!--------------------------------------------------------------------------------

  subroutine HarvestCanopy(I,J,NZ,HvstedLeafC,LeafC_lbrch)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in):: HvstedLeafC
  REAL(R8), intent(in) :: LeafC_lbrch(NumOfCanopyLayers1,JP1,JP1)
  integer :: L,NB,K,NE
  real(r8) :: FHGT,FHVSH,FHVSH1
  real(r8) :: FrcLeafMassLeft
  real(r8) :: WHVSBL
  real(r8) :: FHVSHT

  associate(                                                              &
    CanopyHeightZ_col           => plt_morph%CanopyHeightZ_col,           &
    CutHeightORFrac_pft         => plt_distb%CutHeightORFrac_pft,         &
    FracBiomHarvsted            => plt_distb%FracBiomHarvsted,            &
    LeafChemElmByLayerNode_brch => plt_biom%LeafChemElmByLayerNode_brch,  &
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
          FHGT=AZMAX1(AMIN1(1.0_r8,1._r8-((CanopyHeightZ_col(L))-CutHeightORFrac_pft(NZ))/ &
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
      FrcLeafMassLeft=0._r8
      FHVSH=0._r8
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
        WHVSBL=HvstedLeafC*AZMAX1(LeafC_lbrch(L,NB,NZ))/LeafStrutElms_pft(ielmc,NZ)
      ELSE
        WHVSBL=0._r8
      ENDIF
      D9845: DO K=MaxNodesPerBranch1,0,-1
        IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
          .OR. WHVSBL.GT.0.0_r8)THEN
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
            IF(LeafChemElmByLayerNode_brch(ielmc,L,K,NB,NZ).GT.WHVSBL)THEN
              FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,(LeafChemElmByLayerNode_brch(ielmc,L,K,NB,NZ)-WHVSBL) &
                /LeafChemElmByLayerNode_brch(ielmc,L,K,NB,NZ)))
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

          WHVSBL=WHVSBL-(1._r8-FrcLeafMassLeft)*LeafChemElmByLayerNode_brch(ielmc,L,K,NB,NZ)
          FHVSH1=1._r8-FHVSH
          FHVSHT=FHVSH-FrcLeafMassLeft
          DO NE=1,NumPlantChemElms
            LeafElmntRemoval(NE)=LeafElmntRemoval(NE) &
              +FHVSH1*LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)
            LeafElmntHarv2Litr(NE)=LeafElmntHarv2Litr(NE) &
              +FHVSHT*LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)
            WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE) &
              +FHVSH1*LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
            WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE) &
              +FHVSHT*LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
            LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)=FrcLeafMassLeft*LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)
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
    CanopyLeafAreaZ_pft(L,NZ)=0._r8
    CanopyLeafCLyr_pft(L,NZ)=0._r8
    CanopyStemAreaZ_pft(L,NZ)=CanopyStemAreaZ_pft(L,NZ)*FrcLeafMassLeft
  ENDDO D9865
  end associate
  end subroutine HarvestCanopy
!--------------------------------------------------------------------------------

  subroutine GrazingPlant(I,J,NZ,HvstedLeafC,HvstedShethC,HvstedEarC,HvstedGrainC,&
    WHVSCP,HvstedStalkC,HvstedRsrvC,WHVSHH,WHVSNP)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out):: HvstedLeafC
  real(r8), intent(out):: HvstedShethC
  real(r8), intent(out):: HvstedEarC
  real(r8), intent(out):: HvstedGrainC
  real(r8), intent(out):: WHVSCP
  real(r8), intent(out):: HvstedStalkC
  real(r8), intent(out):: HvstedRsrvC
  real(r8), intent(out) :: WHVSHH
  real(r8), intent(out) :: WHVSNP
  real(r8) :: totShootC,TotPhytomassRemoval
  real(r8) :: WHVSLX,WHVSLY,WHVSCL,WHVSNL,WHVXXX,WHVSSX
  real(r8) :: WHVSTY  
  real(r8) :: WTSTKT  
  real(r8) :: WHVSHY,WHVSCS,WHVSNS,WHVEAX,WHVEAY,WHVGRX,WHVGRY  
  real(r8) :: WHVSHX,WHVHSX,WHVHSY,WHVRVX,WHVRVY,WHVSKX,WHVSTX
  real(r8) :: CCPOLX  
  real(r8) :: CCPLNX  

  associate(                                                       &
    ZERO4Groth_pft           => plt_biom%ZERO4Groth_pft,           &
    CutHeightORFrac_pft      => plt_distb%CutHeightORFrac_pft,     &
    iHarvstType_pft          => plt_distb%iHarvstType_pft,         &
    LeafStrutElms_pft        => plt_biom%LeafStrutElms_pft,        &
    HuskStrutElms_pft        => plt_biom%HuskStrutElms_pft,        &
    ShootStrutElms_pft       => plt_biom%ShootStrutElms_pft,       &
    AREA3                    => plt_site%AREA3,                    &
    NU                       => plt_site%NU,                       &
    FracBiomHarvsted         => plt_distb%FracBiomHarvsted,        &
    THIN_pft                 => plt_distb%THIN_pft,                &
    NoduleNonstructCconc_pft => plt_biom%NoduleNonstructCconc_pft, &
    GrainStrutElms_pft       => plt_biom%GrainStrutElms_pft,       &
    StalkStrutElms_pft       => plt_biom%StalkStrutElms_pft,       &
    StalkRsrvElms_pft        => plt_biom%StalkRsrvElms_pft,        &
    EarStrutElms_pft         => plt_biom%EarStrutElms_pft,         &
    PetoleStrutElms_pft      => plt_biom%PetoleStrutElms_pft,      &
    CanopyNonstElmConc_pft   => plt_biom%CanopyNonstElmConc_pft,   &
    fTCanopyGroth_pft        => plt_pheno%fTCanopyGroth_pft,       &
    AvgCanopyBiomC2Graze_pft => plt_biom%AvgCanopyBiomC2Graze_pft  &
  )
!
!     AvgCanopyBiomC2Graze_pft=average biomass in landscape grazing section
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     TotPhytomassRemoval=total phytomass grazed, removed
!     fTCanopyGroth_pft=temperature function for canopy growth
!     CCPOLP=nonstructural C concentration in canopy
!     NoduleNonstructCconc_pft=nonstructural C concentration in canopy nodules
!
  IF(AvgCanopyBiomC2Graze_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
    TotPhytomassRemoval=CutHeightORFrac_pft(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8 &
      *AREA3(NU)*ShootStrutElms_pft(ielmc,NZ)/AvgCanopyBiomC2Graze_pft(NZ)
  ELSE
    TotPhytomassRemoval=0._r8
  ENDIF
  IF(iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    TotPhytomassRemoval=TotPhytomassRemoval*fTCanopyGroth_pft(NZ)
  ENDIF
  CCPOLX=CanopyNonstElmConc_pft(ielmc,NZ)/(1.0_r8+CanopyNonstElmConc_pft(ielmc,NZ))
  CCPLNX=NoduleNonstructCconc_pft(NZ)/(1.0_r8+NoduleNonstructCconc_pft(NZ))
!
!     LEAF,BACTERIA GRAZED,REMOVED
!
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FracBiomHarvsted(2,1,FracBiomHarvsted(2,2,FracBiomHarvsted(2,3,FracBiomHarvsted(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosyst
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WTLF=PFT leaf C mass
!     WHVXXX=grazing requirement unmet by leaf
!
  WHVSLX=TotPhytomassRemoval*FracBiomHarvsted(1,iplthvst_leaf,NZ)
  WHVSLY=AMIN1(LeafStrutElms_pft(ielmc,NZ),WHVSLX)
  HvstedLeafC=WHVSLY*(1._r8-CCPOLX)
  WHVSCL=WHVSLY*CCPOLX
  WHVSNL=WHVSLY*CCPLNX
  WHVXXX=AZMAX1(WHVSLX-WHVSLY)
  WHVSSX=TotPhytomassRemoval*FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     WHVXXX=grazing requirement unmet by non-foliar removal
!
  totShootC=PetoleStrutElms_pft(ielmc,NZ)+HuskStrutElms_pft(ielmc,NZ)+EarStrutElms_pft(ielmc,NZ)+&
    GrainStrutElms_pft(ielmc,NZ)
  IF(totShootC.GT.ZERO4Groth_pft(NZ))THEN
    WHVSHX=WHVSSX*PetoleStrutElms_pft(ielmc,NZ)/totShootC+WHVXXX
    WHVSHY=AMIN1(PetoleStrutElms_pft(ielmc,NZ),WHVSHX)
    WHVSHH=WHVSHY*(1._r8-CCPOLX)
    WHVSCS=WHVSHY*CCPOLX
    WHVSNS=WHVSHY*CCPLNX
    WHVXXX=AZMAX1(WHVSHX-WHVSHY)
    WHVHSX=WHVSSX*HuskStrutElms_pft(ielmc,NZ)/totShootC+WHVXXX
    WHVHSY=AMIN1(HuskStrutElms_pft(ielmc,NZ),WHVHSX)
    HvstedShethC=WHVHSY
    WHVXXX=AZMAX1(WHVHSX-WHVHSY)
    WHVEAX=WHVSSX*EarStrutElms_pft(ielmc,NZ)/totShootC+WHVXXX
    WHVEAY=AMIN1(EarStrutElms_pft(ielmc,NZ),WHVEAX)
    HvstedEarC=WHVEAY
    WHVXXX=AZMAX1(WHVEAX-WHVEAY)
    WHVGRX=WHVSSX*GrainStrutElms_pft(ielmc,NZ)/totShootC+WHVXXX
    WHVGRY=AMIN1(GrainStrutElms_pft(ielmc,NZ),WHVGRX)
    HvstedGrainC=WHVGRY
    WHVXXX=AZMAX1(WHVGRX-WHVGRY)
  ELSE
    WHVSHH=0._r8
    WHVSCS=0._r8
    WHVSNS=0._r8
    HvstedShethC=0._r8
    HvstedEarC=0._r8
    HvstedGrainC=0._r8
    WHVXXX=WHVXXX+WHVSSX
  ENDIF
  WHVSCP=WHVSCL+WHVSCS
  WHVSNP=WHVSNL+WHVSNS
  WHVSKX=TotPhytomassRemoval*FracBiomHarvsted(1,iplthvst_woody,NZ)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     WHVXXX=grazing requirement unmet by stalk,reserve
!
  WTSTKT=StalkStrutElms_pft(ielmc,NZ)+StalkRsrvElms_pft(ielmc,NZ)
  IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
    WHVSTX=WHVSKX*StalkStrutElms_pft(ielmc,NZ)/WTSTKT+WHVXXX
    WHVSTY=AMIN1(StalkStrutElms_pft(ielmc,NZ),WHVSTX)
    HvstedStalkC=WHVSTY
    WHVXXX=AZMAX1(WHVSTX-WHVSTY)
    WHVRVX=WHVSKX*StalkRsrvElms_pft(ielmc,NZ)/WTSTKT+WHVXXX
    WHVRVY=AMIN1(StalkRsrvElms_pft(ielmc,NZ),WHVRVX)
    HvstedRsrvC=WHVRVY
    WHVXXX=AZMAX1(WHVRVX-WHVRVY)
  ELSE
    HvstedStalkC=0._r8
    HvstedRsrvC=0._r8
    WHVXXX=AZMAX1(WHVSKX)
!
!     ALLOCATE UNMET DEMAND FOR GRAZING TO LEAF,PETIOLE,HUSK
!     EAR,GRAIN
!
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WHVSH*,WHVHS,WHVEA,WHVGR,WHVSC=
!            petiole,husk,ear,grain,nonstructural C removed
!
    IF(WHVXXX.GT.0.0_r8)THEN
      WHVSLY=AMIN1(LeafStrutElms_pft(ielmc,NZ)-HvstedLeafC-WHVSCL,WHVXXX)
      HvstedLeafC=HvstedLeafC+WHVSLY*(1._r8-CCPOLX)
      WHVSCL=WHVSCL+WHVSLY*CCPOLX
      WHVSNL=WHVSNL+WHVSLY*CCPLNX
      WHVXXX=AZMAX1(WHVXXX-WHVSLY)
      IF(totShootC.GT.ZERO4Groth_pft(NZ))THEN
        WHVSHX=WHVXXX*PetoleStrutElms_pft(ielmc,NZ)/totShootC
        WHVSHY=AMIN1(PetoleStrutElms_pft(ielmc,NZ),WHVSHX)
        WHVSHH=WHVSHH+WHVSHY*(1._r8-CCPOLX)
        WHVSCS=WHVSCS+WHVSHY*CCPOLX
        WHVSNS=WHVSNS+WHVSHY*CCPLNX
        WHVXXX=AZMAX1(WHVXXX-WHVSHY)
        WHVHSX=WHVXXX*HuskStrutElms_pft(ielmc,NZ)/totShootC
        WHVHSY=AMIN1(HuskStrutElms_pft(ielmc,NZ),WHVHSX)
        HvstedShethC=HvstedShethC+WHVHSY
        WHVXXX=AZMAX1(WHVXXX-WHVHSY)
        WHVEAX=WHVXXX*EarStrutElms_pft(ielmc,NZ)/totShootC
        WHVEAY=AMIN1(EarStrutElms_pft(ielmc,NZ),WHVEAX)
        HvstedEarC=HvstedEarC+WHVEAY
        WHVXXX=AZMAX1(WHVEAX-WHVEAY)
        WHVGRX=WHVXXX*GrainStrutElms_pft(ielmc,NZ)/totShootC
        WHVGRY=AMIN1(GrainStrutElms_pft(ielmc,NZ),WHVGRX)
        HvstedGrainC=HvstedGrainC+WHVGRY
        WHVXXX=AZMAX1(WHVGRX-WHVGRY)
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine GrazingPlant
!--------------------------------------------------------------------------------

  subroutine PrepBranch4Cut(I,J,NB,NZ,WGLFGX,WGLFGY,FHVSETK,FHVSHK)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(out) :: WGLFGX,WGLFGY
  real(r8), intent(out) :: FHVSETK(0:MaxNodesPerBranch1)  
  real(r8), intent(out) :: FHVSHK(0:MaxNodesPerBranch1)
  integer :: K,NE,L
  real(r8) :: ARLFG
  real(r8) :: WGLFGE(NumPlantChemElms)  

  associate(                                                             &
    CanopyLeafArea_lpft         => plt_morph%CanopyLeafArea_lpft,        &
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch,          &
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,              &
    THIN_pft                    => plt_distb%THIN_pft,                   &
    FracBiomHarvsted            => plt_distb%FracBiomHarvsted,           &
    LeafChemElmByLayerNode_brch => plt_biom%LeafChemElmByLayerNode_brch, &
    CanopyLeafCLyr_pft          => plt_biom%CanopyLeafCLyr_pft,          &
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch,          &
    LeafAreaNode_brch           => plt_morph%LeafAreaNode_brch,          &
    LeafProteinCNode_brch       => plt_biom%LeafProteinCNode_brch,       &
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch,          &
    CanopyLeafAreaZ_pft         => plt_morph%CanopyLeafAreaZ_pft,        &
    iHarvstType_pft             => plt_distb%iHarvstType_pft             &
  )
  WGLFGX=0._r8
  WGLFGY=0._r8

  D9825: DO K=0,MaxNodesPerBranch1
    ARLFG=0._r8
    WGLFGE(1:NumPlantChemElms)=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanopyLeafArea_lpft,CanopyLeafAreaZ_pft=leaf node,total area in canopy layer
!
    D9815: DO L=1,NumOfCanopyLayers1
      ARLFG=ARLFG+CanopyLeafArea_lpft(L,K,NB,NZ)
      DO NE=1,NumPlantChemElms
        WGLFGE(NE)=WGLFGE(NE)+LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)
      ENDDO
      CanopyLeafAreaZ_pft(L,NZ)=CanopyLeafAreaZ_pft(L,NZ)+CanopyLeafArea_lpft(L,K,NB,NZ)
      CanopyLeafCLyr_pft(L,NZ)=CanopyLeafCLyr_pft(L,NZ)+LeafChemElmByLayerNode_brch(ielmc,L,K,NB,NZ)
    ENDDO D9815
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FHVSETK=fraction of internode layer mass not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ).AND.FracBiomHarvsted(1,iplthvst_leaf,NZ).GT.0.0)THEN
        FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(1._r8-(1._r8-AZMAX1(WGLFGE(ielmc)) &
          /LeafElmntNode_brch(ielmc,K,NB,NZ))*FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)/FracBiomHarvsted(1,iplthvst_leaf,NZ))))
        FHVSHK(K)=FHVSETK(K)
      ELSE
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FHVSETK(K)=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)
          FHVSHK(K)=FHVSETK(K)
        ELSE
          FHVSETK(K)=1.0_r8-THIN_pft(NZ)
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
            FHVSHK(K)=1.0_r8-FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
          ELSE
            FHVSHK(K)=FHVSETK(K)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      FHVSETK(K)=0._r8
      FHVSHK(K)=0._r8
    ENDIF
!
!     ACCUMULATE REMAINING BRANCH LEAF AREA, C, N, P
!
!     WGLF=leaf node C mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     LeafAreaLive_brch,LeafAreaNode_brch=branch,node leaf area
!     LeafProteinCNode_brch=leaf protein mass
!
    WGLFGY=WGLFGY+LeafElmntNode_brch(ielmc,K,NB,NZ)
    DO NE=1,NumPlantChemElms
      LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ)+WGLFGE(NE)
    ENDDO
    LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-LeafAreaNode_brch(K,NB,NZ)+ARLFG
    IF(LeafAreaNode_brch(K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*ARLFG/LeafAreaNode_brch(K,NB,NZ)
    ELSE
      LeafProteinCNode_brch(K,NB,NZ)=0._r8
    ENDIF
    LeafAreaNode_brch(K,NB,NZ)=ARLFG
    DO NE=1,NumPlantChemElms
      LeafElmntNode_brch(NE,K,NB,NZ)=WGLFGE(NE)
    ENDDO
    WGLFGX=WGLFGX+LeafElmntNode_brch(ielmc,K,NB,NZ)
  ENDDO D9825
  end associate
  end subroutine PrepBranch4Cut    

!--------------------------------------------------------------------------------
  subroutine CutPlant(I,J,NZ,WHVSHH,WHVSCP,WHVSNP,HvstedShethC,HvstedGrainC,HvstedEarC,HvstedRsrvC,HvstedStalkC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: WHVSHH,WHVSCP,WHVSNP
  real(r8), intent(in) :: HvstedShethC,HvstedGrainC,HvstedEarC,HvstedRsrvC,HvstedStalkC
  integer :: NB
  real(r8) :: FHVSHK(0:MaxNodesPerBranch1),FHVSETK(0:MaxNodesPerBranch1)
  real(r8) :: WGLFGX,WGSHGX,WGLFGY,WGSHGY
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
    VOLWOU                    => plt_site%VOLWOU,                     &
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,       &
    CutHeightORFrac_pft       => plt_distb%CutHeightORFrac_pft,       &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,                &
    CanopyStalkC_pft          => plt_biom%CanopyStalkC_pft,           &
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,       &
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,         &
    UVOLO                     => plt_ew%UVOLO,                        &
    CanopyHeight_pft          => plt_morph%CanopyHeight_pft,          &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,         &
    iHarvstType_pft           => plt_distb%iHarvstType_pft            &
  )
  D9835: DO NB=1,NumOfBranches_pft(NZ)

    CALL PrepBranch4Cut(I,J,NB,NZ,WGLFGX,WGLFGY,FHVSETK,FHVSHK)

    call BranchCutSheathPetole(I,J,NB,NZ,WHVSHH,FHVSETK,FHVSHK,RMedInternodeLen,WGSHGX,WGSHGY)

    call BranchCutNonstructural(I,J,NB,NZ,WGLFGX,WGSHGX,WGLFGY,WGSHGY,WHVSCP,WHVSNP)

  !
  !     CUT STALKS
    call BranchCutPlantStalk(I,J,NB,NZ,RMedInternodeLen,HvstedStalkC,HvstedRsrvC)
!
!     CUT REPRODUCTIVE ORGANS FHVSETH
    call BranchCutReprodOrgans(I,J,NB,NZ,RMedInternodeLen,HvstedShethC,HvstedGrainC,HvstedEarC)

!
!     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
!
!     ShootC4NonstC_brch=total C4 nonstructural C in branch
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     WTLSB=leaf+petiole mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     StalkBiomassC_brch=stalk sapwood mass
!     PSICanopy_pft=canopy water potential
!     CanopyWater_pft=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
    LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))

    call SumPlantBranchBiome(NB,NZ)

    VOLWPX=CanopyWater_pft(NZ)
    WVPLT=AZMAX1(CanopyLeafShethC_pft(NZ)+CanopyStalkC_pft(NZ))

    FDM=get_FDM(PSICanopy_pft(NZ))
!        APSILT=ABS(PSICanopy_pft(NZ))
!        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
    CanopyWater_pft(NZ)=ppmc*WVPLT/FDM

    VOLWOU=VOLWOU+VOLWPX-CanopyWater_pft(NZ)
    UVOLO=UVOLO+VOLWPX-CanopyWater_pft(NZ)
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
      .AND. CanopyHeight_pft(NZ).GT.CutHeightORFrac_pft(NZ))THEN
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
  real(r8), intent(inout) :: FracLeftThin
  real(r8), intent(out):: XHVST1
  real(r8) :: FrcLeafMassNotHarvst(NumPlantChemElms)
  real(r8) :: FFIRE(NumPlantChemElms)
  integer  :: NR,NTG,M,NE

  associate(                                                         &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    Eco_NBP_col               => plt_bgcr%Eco_NBP_col,               &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,    &
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,           &
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,           &
    PO4byFire_pft             => plt_distb%PO4byFire_pft,            &
    N2ObyFire_pft             => plt_distb%N2ObyFire_pft,            &
    NH3byFire_pft             => plt_distb%NH3byFire_pft,            &
    O2ByFire_pft              => plt_distb%O2ByFire_pft,             &
    CH4ByFire_pft             => plt_distb%CH4ByFire_pft,            &
    CO2ByFire_pft             => plt_distb%CO2ByFire_pft,            &
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft,             &
    CSoilOrgM_vr              => plt_soilchem%CSoilOrgM_vr,          &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    THIN_pft                  => plt_distb%THIN_pft,                 &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,          &
    iSoilDisturbType_col      => plt_distb%iSoilDisturbType_col,     &
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
    IF(THETW_vr(L).GT.VolMaxSoilMoist4Fire .OR. CSoilOrgM_vr(ielmc,L).LE.FORGC .OR. iSoilDisturbType_col.NE.22)THEN
      FracLeftThin=1.0_r8
      FFIRE(1:NumPlantChemElms)=0._r8
    ELSE
      FracLeftThin=1.0_r8-DCORP*FracBiomHarvsted(1,iplthvst_woody,NZ) &
        *AMIN1(1.0_r8,(CSoilOrgM_vr(ielmc,L)-FORGC)/(orgcden-FORGC))
      FFIRE(ielmc)=FracBiomHarvsted(2,iplthvst_woody,NZ)
      FFIRE(ielmn)=FFIRE(ielmc)*EFIRE(1,iHarvstType_pft(NZ))
      FFIRE(ielmp)=FFIRE(ielmc)*EFIRE(2,iHarvstType_pft(NZ))
    ENDIF
  ENDIF

  XHVST1=1._r8-FracLeftThin
  D3385: DO M=1,jsken
    DO NE=1,NumPlantChemElms
      FrcLeafMassNotHarvst(NE)=XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootMycoNonstElms_rpvr(NE,N,L,NZ)
      LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
        +(1._r8-FFIRE(NE))*FrcLeafMassNotHarvst(NE)
    ENDDO

    CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
    CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FrcAsCH4byFire*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
    O2ByFire_pft(NZ)=O2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)*2.667_r8
    NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-FFIRE(ielmn)*FrcLeafMassNotHarvst(ielmn)
    N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0_r8
    PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-FFIRE(ielmp)*FrcLeafMassNotHarvst(ielmp)
    CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
    Eco_NBP_col=Eco_NBP_col-FrcAsCH4byFire*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)

    DO NR=1,NumRootAxes_pft(NZ)
      DO NE=1,NumPlantChemElms
        FrcLeafMassNotHarvst(NE)=XHVST1*ElmAllocmat4Litr(NE,icwood,M,NZ)*(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
          +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_woody_litr)
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+&
          (1._r8-FFIRE(NE))*FrcLeafMassNotHarvst(NE)
      ENDDO
      CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
      CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FrcAsCH4byFire*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
      O2ByFire_pft(NZ)=O2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)*2.667_r8
      NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-FFIRE(ielmn)*FrcLeafMassNotHarvst(ielmn)
      N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0
      PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-FFIRE(ielmp)*FrcLeafMassNotHarvst(ielmp)
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
      Eco_NBP_col=Eco_NBP_col-FrcAsCH4byFire*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)

      DO NE=1,NumPlantChemElms
        FrcLeafMassNotHarvst(NE)=XHVST1*ElmAllocmat4Litr(NE,iroot,M,NZ)*(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
          +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_fine_litr)
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
          +(1._r8-FFIRE(NE))*FrcLeafMassNotHarvst(NE)
      ENDDO
      CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
      CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FrcAsCH4byFire*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
      O2ByFire_pft(NZ)=O2ByFire_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)*2.667_r8
      NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-FFIRE(ielmn)*FrcLeafMassNotHarvst(ielmn)
      N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0_r8
      PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-FFIRE(ielmp)*FrcLeafMassNotHarvst(ielmp)
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
      Eco_NBP_col=Eco_NBP_col-FrcAsCH4byFire*FFIRE(ielmc)*FrcLeafMassNotHarvst(ielmc)
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
    Root1stElm_raxs           => plt_biom%Root1stElm_raxs,           &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    Root2ndLen_pvr            => plt_morph%Root2ndLen_pvr,           &
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
    RootNodulStrutElms_pvr    => plt_biom%RootNodulStrutElms_pvr,    &
    RootNodulNonstElms_pvr    => plt_biom%RootNodulNonstElms_pvr,    &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr    &
  )
!
!     REDUCE ROOT STATE VARIABLES DURING HARVESTING
!
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stLen_rpvr,Root2ndLen_pvr=primary,secondary root length
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
        Root1stElm_raxs(NE,N,NR,NZ)=Root1stElm_raxs(NE,N,NR,NZ)*FracLeftThin
      ENDDO
      Root1stLen_rpvr(N,L,NR,NZ)=Root1stLen_rpvr(N,L,NR,NZ)*FracLeftThin
      Root2ndLen_pvr(N,L,NR,NZ)=Root2ndLen_pvr(N,L,NR,NZ)*FracLeftThin
      Root2ndXNum_rpvr(N,L,NR,NZ)=Root2ndXNum_rpvr(N,L,NR,NZ)*FracLeftThin
    ENDDO D3960

    DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)*FracLeftThin
    ENDDO
    RootMycoActiveBiomC_pvr(N,L,NZ)=RootMycoActiveBiomC_pvr(N,L,NZ)*FracLeftThin
      PopuRootMycoC_pvr(N,L,NZ)= PopuRootMycoC_pvr(N,L,NZ)*FracLeftThin
    RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)*FracLeftThin
    Root1stXNumL_pvr(N,L,NZ)=Root1stXNumL_pvr(N,L,NZ)*FracLeftThin
    Root2ndXNum_pvr(N,L,NZ)=Root2ndXNum_pvr(N,L,NZ)*FracLeftThin
    RootLenPerPlant_pvr(N,L,NZ)=RootLenPerPlant_pvr(N,L,NZ)*FracLeftThin
    RootLenDensPerPlant_pvr(N,L,NZ)=RootLenDensPerPlant_pvr(N,L,NZ)*FracLeftThin
    RootPoreVol_pvr(N,L,NZ)=RootPoreVol_pvr(N,L,NZ)*FracLeftThin
    RootVH2O_pvr(N,L,NZ)=RootVH2O_pvr(N,L,NZ)*FracLeftThin
    RootAreaPerPlant_pvr(N,L,NZ)=RootAreaPerPlant_pvr(N,L,NZ)*FracLeftThin
    RootRespPotent_pvr(N,L,NZ)=RootRespPotent_pvr(N,L,NZ)*FracLeftThin
    RootCO2EmisPot_pvr(N,L,NZ)=RootCO2EmisPot_pvr(N,L,NZ)*FracLeftThin
    RootCO2Autor_pvr(N,L,NZ)=RootCO2Autor_pvr(N,L,NZ)*FracLeftThin
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
            XHVST1*(ElmAllocmat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_pvr(NE,L,NZ) &
            +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_pvr(NE,L,NZ))
        ENDDO D3395
        RootNodulStrutElms_pvr(NE,L,NZ)=RootNodulStrutElms_pvr(NE,L,NZ)*FracLeftThin
        RootNodulNonstElms_pvr(NE,L,NZ)=RootNodulNonstElms_pvr(NE,L,NZ)*FracLeftThin
      ENDDO
    ENDIF
  end associate          
  end subroutine HarvstUpdateRootStateL

end module PlantDisturbsMod
