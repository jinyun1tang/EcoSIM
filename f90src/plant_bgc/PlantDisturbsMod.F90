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
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section

! disturbance variables
  real(r8) :: NonstructElmntRemoval(NumPlantChemElmnts)
  real(r8) :: LeafElmntRemoval(NumPlantChemElmnts)
  real(r8) :: FineNonleafElmntRemoval(NumPlantChemElmnts)
  real(r8) :: WoodyElmntRemoval(NumPlantChemElmnts)
  real(r8) :: StandeadElmntRemoval(NumPlantChemElmnts)
  real(r8) :: LeafElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: FineNonleafElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: WoodyElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: StandeadElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: LeafElmntHarv2Litr(NumPlantChemElmnts)
  real(r8) :: PetioleElmntHarv2Litr(NumPlantChemElmnts)
  real(r8) :: WoodyElmntHarv2Litr(NumPlantChemElmnts)
  real(r8) :: StandeadElmntHarv2Litr(NumPlantChemElmnts)
  real(r8) :: WTHTGE(NumPlantChemElmnts)
  real(r8) :: EFIRE(2,5:5)

  public :: RemoveBiomassByDisturbance
  public :: RemoveBiomByManagement
  public :: InitPlantDisturbance
  contains

  subroutine InitPlantDisturbance
  implicit none
  EFIRE=reshape(real((/0.917,0.167/),r8),shape(EFIRE))

  end subroutine InitPlantDisturbance

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByManagement(I,J,NZ,ShootNonstructC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(inout) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
!     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
!

  call RemoveBiomByHarvest(I,J,NZ,ShootNonstructC_brch)
!
!     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
!
  call RemoveBiomByTillage(I,J,NZ,ShootNonstructC_brch)
  end subroutine RemoveBiomByManagement
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomassByDisturbance(I,J,NZ,ShootNonstructC_brch)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), INTENT(INOUT) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  real(r8) :: FHVSE(ielmc)
  real(r8) :: FHVSH
  real(r8) :: WHVSTD
  integer :: M

!     begin_execution
  associate(                            &
    EHVST    =>  plt_distb%EHVST  , &
    HVST     =>  plt_distb%HVST   , &
    iHarvstType_pft    =>  plt_distb%iHarvstType_pft  , &
    THIN_pft     =>  plt_distb%THIN_pft   , &
    PlantPopulation_pft       =>  plt_site%PlantPopulation_pft      , &
    NU       =>  plt_site%NU      , &
    ZERO     =>  plt_site%ZERO    , &
    SolarNoonHour_col   =>  plt_site%SolarNoonHour_col  , &
    AREA3    =>  plt_site%AREA3   , &
    ZEROQ    =>  plt_rbgc%ZEROQ   , &
    ZEROP    =>  plt_biom%ZEROP   , &
    ZEROL    =>  plt_biom%ZEROL   , &
    StandingDeadKCompChemElmnts_pft   =>  plt_biom%StandingDeadKCompChemElmnts_pft  , &
    StandingDeadChemElmnts_pft   =>  plt_biom%StandingDeadChemElmnts_pft    &
  )
!  write(102,*)'iHarvstType_pft',I,iHarvstType_pft(NZ),plt_distb%jHarvst_pft(NZ),NZ
  NonstructElmntRemoval(ielmc)=0._r8;NonstructElmntRemoval(ielmn)=0._r8;NonstructElmntRemoval(ielmp)=0._r8
  LeafElmntRemoval(ielmc)=0._r8;LeafElmntRemoval(ielmn)=0._r8;LeafElmntRemoval(ielmp)=0._r8
  FineNonleafElmntRemoval(ielmc)=0._r8;FineNonleafElmntRemoval(ielmn)=0._r8;FineNonleafElmntRemoval(ielmp)=0._r8
  WoodyElmntRemoval(ielmc)=0._r8;WoodyElmntRemoval(ielmn)=0._r8;WoodyElmntRemoval(ielmp)=0._r8
  StandeadElmntRemoval(ielmc)=0._r8;StandeadElmntRemoval(ielmn)=0._r8;StandeadElmntRemoval(ielmp)=0._r8
  LeafElmnt2Litr(ielmc)=0._r8;LeafElmnt2Litr(ielmn)=0._r8;LeafElmnt2Litr(ielmp)=0._r8
  FineNonleafElmnt2Litr(ielmc)=0._r8;FineNonleafElmnt2Litr(ielmn)=0._r8;FineNonleafElmnt2Litr(ielmp)=0._r8
  WoodyElmnt2Litr(ielmc)=0._r8;WoodyElmnt2Litr(ielmn)=0._r8;WoodyElmnt2Litr(ielmp)=0._r8
  StandeadElmnt2Litr(ielmc)=0._r8;StandeadElmnt2Litr(ielmn)=0._r8;StandeadElmnt2Litr(ielmp)=0._r8
  LeafElmntHarv2Litr(ielmc)=0._r8;LeafElmntHarv2Litr(ielmn)=0._r8;LeafElmntHarv2Litr(ielmp)=0._r8
  PetioleElmntHarv2Litr(ielmc)=0._r8;PetioleElmntHarv2Litr(ielmn)=0._r8;PetioleElmntHarv2Litr(ielmp)=0._r8
  WoodyElmntHarv2Litr(ielmc)=0._r8;WoodyElmntHarv2Litr(ielmn)=0._r8;WoodyElmntHarv2Litr(ielmp)=0._r8
  StandeadElmntHarv2Litr(ielmc)=0._r8;StandeadElmntHarv2Litr(ielmn)=0._r8;StandeadElmntHarv2Litr(ielmp)=0._r8
  WTHTGE(ielmc)=0._r8;WTHTGE(ielmn)=0._r8;WTHTGE(ielmp)=0._r8
!     ENDIF

!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     THIN_pft=thinning:fraction of population removed
!     FHVSE(ielmc)=fraction of standing dead mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!
  IF(iHarvstType_pft(NZ).GE.iharvtyp_none)THEN
    IF(J.EQ.INT(SolarNoonHour_col).AND.iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.&
      iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FHVSE(ielmc)=AZMAX1(1._r8-EHVST(1,4,NZ))
        FHVSH=FHVSE(ielmc)
      ELSE
        FHVSE(ielmc)=AZMAX1(1._r8-THIN_pft(NZ))
        IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
          FHVSH=AZMAX1(1._r8-EHVST(1,4,NZ)*THIN_pft(NZ))
        ELSE
          FHVSH=FHVSE(ielmc)
        ENDIF
      ENDIF
    ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
      IF(StandingDeadChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
        WHVSTD=HVST(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8*AREA3(NU)*EHVST(1,4,NZ)
        FHVSE(ielmc)=AZMAX1(1._r8-WHVSTD/StandingDeadChemElmnts_pft(ielmc,NZ))
        FHVSH=FHVSE(ielmc)
      ELSE
        FHVSE(ielmc)=1.0_r8
        FHVSH=1.0_r8
      ENDIF
    ELSE
      FHVSE(ielmc)=1.0_r8
      FHVSH=1.0_r8
    ENDIF
    D6475: DO M=1,jsken
      StandeadElmntRemoval(ielmc)=StandeadElmntRemoval(ielmc)+(1._r8-FHVSH)*StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)
      StandeadElmntRemoval(ielmn)=StandeadElmntRemoval(ielmn)+(1._r8-FHVSH)*StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)
      StandeadElmntRemoval(ielmp)=StandeadElmntRemoval(ielmp)+(1._r8-FHVSH)*StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)
      StandeadElmntHarv2Litr(ielmc)=StandeadElmntHarv2Litr(ielmc)+(FHVSH-FHVSE(ielmc))*StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)
      StandeadElmntHarv2Litr(ielmn)=StandeadElmntHarv2Litr(ielmn)+(FHVSH-FHVSE(ielmc))*StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)
      StandeadElmntHarv2Litr(ielmp)=StandeadElmntHarv2Litr(ielmp)+(FHVSH-FHVSE(ielmc))*StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)=FHVSE(ielmc)*StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)=FHVSE(ielmc)*StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)=FHVSE(ielmc)*StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)
    ENDDO D6475
!
    call PlantDisturbance(I,J,NZ)

    ZEROP(NZ)=ZERO*PlantPopulation_pft(NZ)
    ZEROQ(NZ)=ZERO*PlantPopulation_pft(NZ)/AREA3(NU)
    ZEROL(NZ)=ZERO*PlantPopulation_pft(NZ)*1.0E+06_r8
  ENDIF
  end associate
  end subroutine RemoveBiomassByDisturbance


!------------------------------------------------------------------------------------------
  subroutine PlantDisturbance(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  real(r8) :: NonstructElmntOffEcosystem(NumPlantChemElmnts)
  real(r8) :: LeafElmntOffEcosystem(NumPlantChemElmnts)
  real(r8) :: FineNonleafElmntOffEcosystem(NumPlantChemElmnts)
  real(r8) :: WoodyElmntOffEcosystem(NumPlantChemElmnts)
  real(r8) :: StandeadElmntOffEcosystem(NumPlantChemElmnts)
  real(r8) :: NonstructElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: TotalElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: HarvestElmnt2Litr(NumPlantChemElmnts)

  NonstructElmntOffEcosystem(1:NumPlantChemElmnts)=0._r8
  LeafElmntOffEcosystem(1:NumPlantChemElmnts)=0._r8
  FineNonleafElmntOffEcosystem(1:NumPlantChemElmnts)=0._r8
  WoodyElmntOffEcosystem(1:NumPlantChemElmnts)=0._r8
  StandeadElmntOffEcosystem(1:NumPlantChemElmnts)=0._r8
  NonstructElmnt2Litr(1:NumPlantChemElmnts)=0._r8

  call ApplyDisturbanceBiomRemoval(I,J,NZ,NonstructElmnt2Litr,NonstructElmntOffEcosystem,&
    LeafElmntOffEcosystem,FineNonleafElmntOffEcosystem,WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call TotalBiomRemovalByDisturbance(I,J,NZ,NonstructElmnt2Litr,HarvestElmnt2Litr,TotalElmnt2Litr)
!
!     ABOVE-GROUND LITTERFALL FROM HARVESTING
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
  real(r8), intent(in) :: HarvestElmnt2Litr(NumPlantChemElmnts),TotalElmnt2Litr(NumPlantChemElmnts)
  real(r8), intent(in) :: NonstructElmnt2Litr(NumPlantChemElmnts)
  real(r8), intent(in) :: NonstructElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(in) :: LeafElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(in) :: FineNonleafElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(in) :: WoodyElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(in) :: StandeadElmntOffEcosystem(NumPlantChemElmnts)
  integer :: M,NE
!     begin_execution
  associate(                        &
    ilignin  =>  pltpar%ilignin   , &
    inonstruct =>  pltpar%inonstruct  , &
    k_fine_litr=> pltpar%k_fine_litr,&
    k_woody_litr=> pltpar%k_woody_litr,&
    ifoliar  =>  pltpar%ifoliar   , &
    inonfoliar =>  pltpar%inonfoliar  , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    icwood   =>  pltpar%icwood    , &
    StandingDeadKCompChemElmnts_pft   =>  plt_biom%StandingDeadKCompChemElmnts_pft  , &
    iHarvstType_pft    =>  plt_distb%iHarvstType_pft  , &
    FWOODE   =>  plt_allom%FWOODE , &
    CFOPE    =>  plt_soilchem%CFOPE, &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr    , &
    LitrfallChemElmnts_pft    =>  plt_bgcr%LitrfallChemElmnts_pft   , &
    SurfLitrfallChemElmnts_pft    =>  plt_bgcr%SurfLitrfallChemElmnts_pft   , &
    iPlantTurnoverPattern_pft   =>  plt_pheno%iPlantTurnoverPattern_pft , &
    iPlantMorphologyType_pft   =>  plt_pheno%iPlantMorphologyType_pft   &
  )
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      D6375: DO M=1,jsken
        DO NE=1,NumPlantChemElmnts        
          LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,inonstruct,M,NZ)*(NonstructElmnt2Litr(NE)) &
            +CFOPE(NE,ifoliar,M,NZ)*(LeafElmnt2Litr(NE)+LeafElmntHarv2Litr(NE)) &
            +CFOPE(NE,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(NE)+PetioleElmntHarv2Litr(NE))

          IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType_pft(NZ))))THEN
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,istalk,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmnt2Litr(NE)&
              +StandeadElmntHarv2Litr(NE))
          ELSE
            StandingDeadKCompChemElmnts_pft(NE,M,NZ)=StandingDeadKCompChemElmnts_pft(NE,M,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE))

            LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE))*FWOODE(NE,k_woody_litr)

            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE))*FWOODE(NE,k_fine_litr)
          ENDIF
        ENDDO
      ENDDO D6375

!
!     ABOVE-GROUND LITTERFALL FROM FIRE
!
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
    ELSE
      D6485: DO M=1,jsken
        LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmc,inonstruct,M,NZ)*(NonstructElmnt2Litr(ielmc)) &
          +CFOPE(ielmc,ifoliar,M,NZ)*(LeafElmnt2Litr(ielmc)+LeafElmntHarv2Litr(ielmc)) &
          +CFOPE(ielmc,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmc)+PetioleElmntHarv2Litr(ielmc))
        LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmn,inonstruct,M,NZ)*NonstructElmntOffEcosystem(ielmn) &
          +CFOPE(ielmn,ifoliar,M,NZ)*LeafElmntOffEcosystem(ielmn) &
          +CFOPE(ielmn,inonfoliar,M,NZ)*FineNonleafElmntOffEcosystem(ielmn)
        LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmp,inonstruct,M,NZ)*NonstructElmntOffEcosystem(ielmp) &
          +CFOPE(ielmp,ifoliar,M,NZ)*LeafElmntOffEcosystem(ielmp) &
          +CFOPE(ielmp,inonfoliar,M,NZ)*FineNonleafElmntOffEcosystem(ielmp)

        LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ) &
          +CFOPE(ielmn,inonstruct,M,NZ)*(NonstructElmnt2Litr(ielmn)-NonstructElmntOffEcosystem(ielmn)) &
          +CFOPE(ielmn,ifoliar,M,NZ)*(LeafElmnt2Litr(ielmn)+LeafElmntHarv2Litr(ielmn)-LeafElmntOffEcosystem(ielmn)) &
          +CFOPE(ielmn,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmn)+PetioleElmntHarv2Litr(ielmn)&
          -FineNonleafElmntOffEcosystem(ielmn))
        LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ) &
          +CFOPE(ielmp,inonstruct,M,NZ)*(NonstructElmnt2Litr(ielmp)-NonstructElmntOffEcosystem(ielmp)) &
          +CFOPE(ielmp,ifoliar,M,NZ)*(LeafElmnt2Litr(ielmp)+LeafElmntHarv2Litr(ielmp)-LeafElmntOffEcosystem(ielmp)) &
          +CFOPE(ielmp,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmp)+PetioleElmntHarv2Litr(ielmp)&
          -FineNonleafElmntOffEcosystem(ielmp))

        IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType_pft(NZ))))THEN
          LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)+&
            CFOPE(ielmc,istalk,M,NZ)*(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc)+StandeadElmnt2Litr(ielmc)+&
            StandeadElmntHarv2Litr(ielmc))
          LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)+&
            CFOPE(ielmn,istalk,M,NZ)*(WoodyElmntOffEcosystem(ielmn)+StandeadElmntOffEcosystem(ielmn))
          LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)+&
            CFOPE(ielmp,istalk,M,NZ)*(WoodyElmntOffEcosystem(ielmp)+StandeadElmntOffEcosystem(ielmp))
          LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ) &
            +CFOPE(ielmn,istalk,M,NZ)*(WoodyElmnt2Litr(ielmn)+WoodyElmntHarv2Litr(ielmn) &
            -WoodyElmntOffEcosystem(ielmn)+StandeadElmnt2Litr(ielmn)+StandeadElmntHarv2Litr(ielmn)&
            -StandeadElmntOffEcosystem(ielmn))
          LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)+&
            CFOPE(ielmp,istalk,M,NZ)*(WoodyElmnt2Litr(ielmp)+WoodyElmntHarv2Litr(ielmp)- &
            WoodyElmntOffEcosystem(ielmp)+StandeadElmnt2Litr(ielmp)+StandeadElmntHarv2Litr(ielmp)-&
            StandeadElmntOffEcosystem(ielmp))
        ELSE
          StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)=StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)+CFOPE(ielmc,icwood,M,NZ) &
            *(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc))
          StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)=StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)+CFOPE(ielmn,icwood,M,NZ) &
            *WoodyElmntOffEcosystem(ielmn)
          StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)=StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)+CFOPE(ielmp,icwood,M,NZ) &
            *WoodyElmntOffEcosystem(ielmp)
          LitterFallChemElmnt_pftvr(ielmc,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_woody_litr,0,NZ) &
            *CFOPE(ielmc,istalk,M,NZ)&
            *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FWOODE(ielmc,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmn,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_woody_litr,0,NZ) &
            +CFOPE(ielmn,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmn)*FWOODE(ielmn,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmp,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_woody_litr,0,NZ) &
            +CFOPE(ielmp,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmp)*FWOODE(ielmp,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmn,ilignin,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_woody_litr,0,NZ) &
            +CFOPE(ielmn,icwood,M,NZ) &
            *(WoodyElmnt2Litr(ielmn)+WoodyElmntHarv2Litr(ielmn)-WoodyElmntOffEcosystem(ielmn) &
            +StandeadElmnt2Litr(ielmn)+StandeadElmntHarv2Litr(ielmn)-StandeadElmntOffEcosystem(ielmn))*FWOODE(ielmn,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmp,ilignin,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_woody_litr,0,NZ) &
            +CFOPE(ielmp,icwood,M,NZ) &
            *(WoodyElmnt2Litr(ielmp)+WoodyElmntHarv2Litr(ielmp)-WoodyElmntOffEcosystem(ielmp) &
            +StandeadElmnt2Litr(ielmp)+StandeadElmntHarv2Litr(ielmp)-StandeadElmntOffEcosystem(ielmp)) &
            *FWOODE(ielmp,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ) &
            +CFOPE(ielmc,istalk,M,NZ) &
            *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FWOODE(ielmc,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ) &
            +CFOPE(ielmn,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmn)*FWOODE(ielmn,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ) &
            +CFOPE(ielmp,istalk,M,NZ)*StandeadElmntOffEcosystem(ielmp)*FWOODE(ielmp,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ) &
            +CFOPE(ielmn,icwood,M,NZ)*(WoodyElmnt2Litr(ielmn)+WoodyElmntHarv2Litr(ielmn)-WoodyElmntOffEcosystem(ielmn) &
            +StandeadElmnt2Litr(ielmn)+StandeadElmntHarv2Litr(ielmn)-StandeadElmntOffEcosystem(ielmn))*FWOODE(ielmn,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ) &
            +CFOPE(ielmp,icwood,M,NZ)*(WoodyElmnt2Litr(ielmp)+WoodyElmntHarv2Litr(ielmp)-WoodyElmntOffEcosystem(ielmp) &
            +StandeadElmnt2Litr(ielmp)+StandeadElmntHarv2Litr(ielmp)-StandeadElmntOffEcosystem(ielmp))*FWOODE(ielmp,k_fine_litr)
        ENDIF
      ENDDO D6485
    ENDIF
  ELSE
!
!     ABOVE-GROUND LITTERFALL FROM GRAZING
!
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!
    DO NE=1,NumPlantChemElmnts
      LitrfallChemElmnts_pft(NE,NZ)=LitrfallChemElmnts_pft(NE,NZ)+TotalElmnt2Litr(NE)+HarvestElmnt2Litr(NE)
      SurfLitrfallChemElmnts_pft(NE,NZ)=SurfLitrfallChemElmnts_pft(NE,NZ)+TotalElmnt2Litr(NE)+HarvestElmnt2Litr(NE)
    ENDDO
  ENDIF
  end associate
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

  subroutine TotalBiomRemovalByDisturbance(I,J,NZ,NonstructElmnt2Litr,HarvestElmnt2Litr,TotalElmnt2Litr)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: NonstructElmnt2Litr(NumPlantChemElmnts)
  real(r8), intent(out):: HarvestElmnt2Litr(NumPlantChemElmnts),TotalElmnt2Litr(NumPlantChemElmnts)
  real(r8) :: TotalElmntRemoval(NumPlantChemElmnts)
  integer :: NE
!     begin_execution
  associate(                            &
    iHarvstType_pft    =>  plt_distb%iHarvstType_pft  , &
    jHarvst_pft    =>  plt_distb%jHarvst_pft  , &
    EcoHavstElmnt_pft    =>  plt_distb%EcoHavstElmnt_pft  , &
    EcoHavstElmnt_col   =>  plt_distb%EcoHavstElmnt_col , &
    NH3byFire_pft    =>  plt_distb%NH3byFire_pft  , &
    PO4byFire_pft    =>  plt_distb%PO4byFire_pft  , &
    CH4ByFire_pft    =>  plt_distb%CH4ByFire_pft  , &
    VOXYF    =>  plt_distb%VOXYF  , &
    N2ObyFire_pft    =>  plt_distb%N2ObyFire_pft  , &
    CO2ByFire_pft    =>  plt_distb%CO2ByFire_pft  , &
    CO2NetFix_pft     =>  plt_bgcr%CO2NetFix_pft    , &
    Eco_NBP_col     =>  plt_bgcr%Eco_NBP_col    , &
    Eco_AutoR_col     =>  plt_bgcr%Eco_AutoR_col    , &
    ECO_ER_col     =>  plt_bgcr%ECO_ER_col    , &
    CanopyPlusNoduRespC_pft    =>  plt_bgcr%CanopyPlusNoduRespC_pft   , &
    GrossResp_pft    =>  plt_bgcr%GrossResp_pft   , &
    NonstructalElmnts_pft    =>  plt_biom%NonstructalElmnts_pft     &
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
  DO NE=1,NumPlantChemElmnts
    TotalElmntRemoval(NE)=NonstructElmntRemoval(NE)+LeafElmntRemoval(NE)+FineNonleafElmntRemoval(NE)+&
      WoodyElmntRemoval(NE)+StandeadElmntRemoval(NE)
    TotalElmnt2Litr(NE)=NonstructElmnt2Litr(NE)+LeafElmnt2Litr(NE)+FineNonleafElmnt2Litr(NE)+&
      WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE)
    HarvestElmnt2Litr(NE)=LeafElmntHarv2Litr(NE)+PetioleElmntHarv2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE)
  ENDDO
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      IF(jHarvst_pft(NZ).NE.jharvtyp_tmareseed)THEN
        DO NE=1,NumPlantChemElmnts
          EcoHavstElmnt_pft(NE,NZ)=EcoHavstElmnt_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
          EcoHavstElmnt_col(NE)=EcoHavstElmnt_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
        ENDDO
        Eco_NBP_col=Eco_NBP_col+TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc)
      ELSE
        DO NE=1,NumPlantChemElmnts
          NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
        ENDDO
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     CO2ByFire_pft,CH4ByFire_pft,VOXYF,NH3byFire_pft,N2ObyFire_pft,PO4byFire_pft=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_col=total net biome productivity
!
    ELSE
      CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
      CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
      VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))*2.667
      NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-TotalElmntRemoval(ielmn)+TotalElmnt2Litr(ielmn)
      N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0_r8
      PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-TotalElmntRemoval(ielmp)+TotalElmnt2Litr(ielmp)
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
      Eco_NBP_col=Eco_NBP_col-FCH4F*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
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
    DO NE=2,NumPlantChemElmnts
      EcoHavstElmnt_pft(NE,NZ)=EcoHavstElmnt_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
      EcoHavstElmnt_col(NE)=EcoHavstElmnt_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
    ENDDO
    GrossResp_pft(NZ)=GrossResp_pft(NZ)-GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
    CanopyPlusNoduRespC_pft(NZ)=CanopyPlusNoduRespC_pft(NZ)-GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
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
  real(r8), intent(out) :: NonstructElmnt2Litr(NumPlantChemElmnts)
  real(r8), intent(out) :: NonstructElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(out) :: LeafElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(out) :: FineNonleafElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(out) :: WoodyElmntOffEcosystem(NumPlantChemElmnts)
  real(r8), intent(out) :: StandeadElmntOffEcosystem(NumPlantChemElmnts)

  real(r8) :: EHVST21,EHVST22,EHVST23,EHVST24
  real(r8) :: EHVST21h,EHVST22h,EHVST23h,EHVST24h
  integer  :: NE
!     begin_execution
  associate(                            &
    NU       =>  plt_site%NU      , &
    AREA3    =>  plt_site%AREA3   , &
    EHVST    =>  plt_distb%EHVST  , &
    FERT     =>  plt_distb%FERT   , &
    IYTYP    =>  plt_distb%IYTYP  , &
    iHarvstType_pft    =>  plt_distb%iHarvstType_pft    &
  )
!     IF NO PLANT C,N,P REMOVED AT HARVEST (ALL RESIDUE RETURNED)
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  EHVST21=1._r8-EHVST(2,iplthvst_leaf,NZ)
  EHVST22=1._r8-EHVST(2,iplthvst_finenonleaf,NZ)
  EHVST23=1._r8-EHVST(2,iplthvst_woody,NZ)
  EHVST24=1._r8-EHVST(2,iplthvst_stdead,NZ)

  IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
    DO NE=1,NumPlantChemElmnts
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
    DO NE=1,NumPlantChemElmnts
      NonstructElmnt2Litr(NE)=NonstructElmntRemoval(NE)
      LeafElmnt2Litr(NE)=LeafElmntRemoval(NE)
      FineNonleafElmnt2Litr(NE)=FineNonleafElmntRemoval(NE)-WTHTGE(NE)*EHVST(2,iplthvst_finenonleaf,NZ)
      WoodyElmnt2Litr(NE)=WoodyElmntRemoval(NE)
      StandeadElmnt2Litr(NE)=StandeadElmntRemoval(NE)
    ENDDO
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_allabv)THEN
    DO NE=1,NumPlantChemElmnts
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
    DO NE=1,NumPlantChemElmnts
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
    EHVST21h=1._r8-EHVST(2,iplthvst_leaf,NZ)*0.5_r8
    EHVST22h=1._r8-EHVST(2,iplthvst_finenonleaf,NZ)*0.5_r8
    EHVST23h=1._r8-EHVST(2,iplthvst_woody,NZ)*0.5_r8
    EHVST24h=1._r8-EHVST(2,iplthvst_stdead,NZ)*0.5_r8

    NonstructElmnt2Litr(ielmc)=NonstructElmntRemoval(ielmc)*EHVST21
    LeafElmnt2Litr(ielmc)=LeafElmntRemoval(ielmc)*EHVST21
    FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
    WoodyElmnt2Litr(ielmc)=WoodyElmntRemoval(ielmc)*EHVST23
    StandeadElmnt2Litr(ielmc)=StandeadElmntRemoval(ielmc)*EHVST24

    DO NE=2,NumPlantChemElmnts
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
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_fire)THEN

    NonstructElmnt2Litr(ielmc)=NonstructElmntRemoval(ielmc)*EHVST21
    NonstructElmnt2Litr(ielmn)=NonstructElmntRemoval(ielmn)*&
      (1._r8-EFIRE(1,iHarvstType_pft(NZ))*EHVST(2,iplthvst_leaf,NZ))
    NonstructElmnt2Litr(ielmp)=NonstructElmntRemoval(ielmp)*&
      (1._r8-EFIRE(2,iHarvstType_pft(NZ))*EHVST(2,iplthvst_leaf,NZ))
    NonstructElmntOffEcosystem(ielmn)=NonstructElmntRemoval(ielmn)*EHVST21
    NonstructElmntOffEcosystem(ielmp)=NonstructElmntRemoval(ielmp)*EHVST21

    LeafElmnt2Litr(ielmc)=LeafElmntRemoval(ielmc)*EHVST21
    LeafElmnt2Litr(ielmn)=LeafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*EHVST(2,iplthvst_leaf,NZ))
    LeafElmnt2Litr(ielmp)=LeafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*EHVST(2,iplthvst_leaf,NZ))
    LeafElmntOffEcosystem(ielmn)=LeafElmntRemoval(ielmn)*EHVST21
    LeafElmntOffEcosystem(ielmp)=LeafElmntRemoval(ielmp)*EHVST21

    FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
    FineNonleafElmnt2Litr(ielmn)=FineNonleafElmntRemoval(ielmn)*&
      (1._r8-EFIRE(1,iHarvstType_pft(NZ))*EHVST(2,iplthvst_finenonleaf,NZ))
    FineNonleafElmnt2Litr(ielmp)=FineNonleafElmntRemoval(ielmp)*&
      (1._r8-EFIRE(2,iHarvstType_pft(NZ))*EHVST(2,iplthvst_finenonleaf,NZ))
    FineNonleafElmntOffEcosystem(ielmn)=FineNonleafElmntRemoval(ielmn)*EHVST22
    FineNonleafElmntOffEcosystem(ielmp)=FineNonleafElmntRemoval(ielmp)*EHVST22

    WoodyElmnt2Litr(ielmc)=WoodyElmntRemoval(ielmc)*EHVST23
    WoodyElmnt2Litr(ielmn)=WoodyElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*EHVST(2,iplthvst_woody,NZ))
    WoodyElmnt2Litr(ielmp)=WoodyElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*EHVST(2,iplthvst_woody,NZ))
    WoodyElmntOffEcosystem(ielmn)=WoodyElmntRemoval(ielmn)*EHVST23
    WoodyElmntOffEcosystem(ielmp)=WoodyElmntRemoval(ielmp)*EHVST23

    StandeadElmnt2Litr(ielmc)=StandeadElmntRemoval(ielmc)*EHVST24
    StandeadElmnt2Litr(ielmn)=StandeadElmntRemoval(ielmn)*&
      (1._r8-EFIRE(1,iHarvstType_pft(NZ))*EHVST(2,iplthvst_stdead,NZ))
    StandeadElmnt2Litr(ielmp)=StandeadElmntRemoval(ielmp)*&
      (1._r8-EFIRE(2,iHarvstType_pft(NZ))*EHVST(2,iplthvst_stdead,NZ))
    StandeadElmntOffEcosystem(ielmn)=StandeadElmntRemoval(ielmn)*EHVST24
    StandeadElmntOffEcosystem(ielmp)=StandeadElmntRemoval(ielmp)*EHVST24
  ENDIF
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByTillage(I,J,NZ,ShootNonstructC_brch)
  use SurfLitterDataType, only : XCORP
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(inout) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,K,M,N,NR,NE,NB,NTG
  real(r8) :: XHVST,XHVST1
  REAL(R8) :: APSILT
  real(r8) :: FDM,VOLWPX
  real(r8) :: WVPLT
!     begin_execution
  associate(                               &
    jHarvst_pft    =>  plt_distb%jHarvst_pft     , &
    iDayPlantHarvest_pft    =>  plt_distb%iDayPlantHarvest_pft     , &
    iDayPlanting_pft    =>  plt_distb%iDayPlanting_pft     , &
    ITILL    =>  plt_distb%ITILL     , &
    iYearPlanting_pft     =>  plt_distb%iYearPlanting_pft      , &
    iYearPlantHarvest_pft     =>  plt_distb%iYearPlantHarvest_pft      , &
    XCORP    =>  plt_distb%XCORP     , &
    CFOPE    =>  plt_soilchem%CFOPE  , &
    trcg_rootml_vr     =>  plt_rbgc%trcg_rootml_vr       , &
    trcs_rootml_vr => plt_rbgc%trcs_rootml_vr , &
    UVOLO    =>  plt_ew%UVOLO        , &
    CanopyWater_pft    =>  plt_ew%CanopyWater_pft        , &
    VHeatCapCanP    =>  plt_ew%VHeatCapCanP        , &
    PSICanopy_pft   =>  plt_ew%PSICanopy_pft       , &
    PPX      =>  plt_site%PPX        , &
     RootMycoNonstructElmnt_vr   =>  plt_biom%RootMycoNonstructElmnt_vr     , &
    RootProteinC_pvr    =>  plt_biom%RootProteinC_pvr      , &
     PopuPlantRootC_vr    =>  plt_biom% PopuPlantRootC_vr      , &
    RootStructBiomC_vr   =>  plt_biom%RootStructBiomC_vr     , &
    LeafChemElmnts_brch  =>  plt_biom%LeafChemElmnts_brch    , &
    GrainChemElmnts_brch   =>  plt_biom%GrainChemElmnts_brch     , &
    EarChemElmnts_brch  =>  plt_biom%EarChemElmnts_brch    , &
    NoduleNonstructElmnt_brch   =>  plt_biom%NoduleNonstructElmnt_brch     , &
    NonstructElmnt_brch   =>  plt_biom%NonstructElmnt_brch     , &
    HuskChemElmnts_brch  =>  plt_biom%HuskChemElmnts_brch    , &
    ReserveElmnts_brch  =>  plt_biom%ReserveElmnts_brch    , &
    CanopyNoduleChemElmnt_brch   =>  plt_biom%CanopyNoduleChemElmnt_brch     , &
    ShootChemElmnt_brch  =>  plt_biom%ShootChemElmnt_brch    , &
    StalkChemElmnts_brch  =>  plt_biom%StalkChemElmnts_brch    , &
    LeafChemElmntByLayer_pft   =>  plt_biom%LeafChemElmntByLayer_pft     , &
    PetoleChemElmnt_brch =>  plt_biom%PetoleChemElmnt_brch   , &
    LeafElmntNode_brch    =>  plt_biom%LeafElmntNode_brch      , &
    LeafPetolBiomassC_brch    =>  plt_biom%LeafPetolBiomassC_brch      , &
    BranchStalkChemElmnts_pft  =>  plt_biom%BranchStalkChemElmnts_pft    , &
    StalkBiomassC_brch   =>  plt_biom%StalkBiomassC_brch     , &
    PetioleProteinCNode_brch   =>  plt_biom%PetioleProteinCNode_brch     , &
    PetioleElmntNode_brch   =>  plt_biom%PetioleElmntNode_brch     , &
    Root1stStructChemElmnt_pvr   =>  plt_biom%Root1stStructChemElmnt_pvr     , &
    LeafProteinCNode_brch     =>  plt_biom%LeafProteinCNode_brch       , &
    InternodeChemElmnt_brch   =>  plt_biom%InternodeChemElmnt_brch     , &
    CanopyStalkC_pft    =>  plt_biom%CanopyStalkC_pft      , &
    Root1stChemElmnt   =>  plt_biom%Root1stChemElmnt     , &
    NonstructalElmnts_pft    =>  plt_biom%NonstructalElmnts_pft      , &
    CanopyLeafShethC_pft     =>  plt_biom%CanopyLeafShethC_pft       , &
    Root2ndStructChemElmnt_pvr   =>  plt_biom%Root2ndStructChemElmnt_pvr     , &
    RootNoduleNonstructElmnt_vr  =>  plt_biom%RootNoduleNonstructElmnt_vr    , &
    RootNodueChemElmnt_pvr   =>  plt_biom%RootNodueChemElmnt_pvr     , &
    GrainSeedBiomCMean_brch    =>  plt_allom%GrainSeedBiomCMean_brch     , &
    FWOODE   =>  plt_allom%FWOODE    , &
    FWODBE   =>  plt_allom%FWODBE    , &
    FWODLE   =>  plt_allom%FWODLE    , &
    FWODRE   =>  plt_allom%FWODRE    , &
    iPlantBranchState_brch    =>  plt_pheno%iPlantBranchState_brch     , &
    iPlantPhenologyPattern_pft   =>  plt_pheno%iPlantPhenologyPattern_pft    , &
    iPlantState_pft    =>  plt_pheno%iPlantState_pft     , &
    iPlantMorphologyType_pft   =>  plt_pheno%iPlantMorphologyType_pft    , &
    iPlantTurnoverPattern_pft   =>  plt_pheno%iPlantTurnoverPattern_pft    , &
    iPlantPhenologyType_pft   =>  plt_pheno%iPlantPhenologyType_pft    , &
    iPlantRootState_pft   =>  plt_pheno%iPlantRootState_pft    , &
    iPlantShootState_pft    =>  plt_pheno%iPlantShootState_pft     , &
    CMassHCO3BundleSheath_node     =>  plt_photo%CMassHCO3BundleSheath_node      , &
    CMassCO2BundleSheath_node     =>  plt_photo%CMassCO2BundleSheath_node      , &
    CPOOL3   =>  plt_photo%CPOOL3    , &
    CPOOL4   =>  plt_photo%CPOOL4    , &
    MaxNumRootLays       =>  plt_site%MaxNumRootLays         , &
    PlantPopulation_pft       =>  plt_site%PlantPopulation_pft         , &
    iYearCurrent     =>  plt_site%iYearCurrent       , &
    SolarNoonHour_col   =>  plt_site%SolarNoonHour_col     , &
    VOLWOU   =>  plt_site%VOLWOU     , &
    NU       =>  plt_site%NU         , &
    k_fine_litr=> pltpar%k_fine_litr ,&
    k_woody_litr=> pltpar%k_woody_litr,&
    inonstruct =>  pltpar%inonstruct     , &
    ifoliar  =>  pltpar%ifoliar      , &
    istalk   =>  pltpar%istalk       , &
    iroot    =>  pltpar%iroot        , &
    inonfoliar =>  pltpar%inonfoliar     , &
    icwood   =>  pltpar%icwood       , &
    RootGasLossDisturb_pft    =>   plt_bgcr%RootGasLossDisturb_pft    , &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr       , &
    RCO2A    =>  plt_rbgc%RCO2A      , &
    RootRespPotential_vr    =>  plt_rbgc%RootRespPotential_vr      , &
    RCO2N    =>  plt_rbgc%RCO2N      , &
    FracRadPARbyCanopy_pft    =>  plt_rad%FracRadPARbyCanopy_pft       , &
    PrimRootLen    =>  plt_morph%PrimRootLen     , &
    RootVH2O_vr   =>  plt_morph%RootVH2O_vr    , &
    RootAreaPerPlant_vr    =>  plt_morph%RootAreaPerPlant_vr     , &
    RootVolume_vr    =>  plt_morph%RootVolume_vr     , &
    RootLenDensPerPlant_pvr    =>  plt_morph%RootLenDensPerPlant_pvr     , &
    PrimRootXNumL_pvr    =>  plt_morph%PrimRootXNumL_pvr     , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType    , &
    RootLenPerPlant_pvr    =>  plt_morph%RootLenPerPlant_pvr     , &
    SecndRootXNum_rpvr    =>  plt_morph%SecndRootXNum_rpvr     , &
    SecndRootLen    =>  plt_morph%SecndRootLen     , &
    SecndRootXNum_pvr     =>  plt_morph%SecndRootXNum_pvr      , &
    NGTopRootLayer_pft      =>  plt_morph%NGTopRootLayer_pft       , &
    MY       =>  plt_morph%MY        , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft      , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft       , &
    LeafAreaNode_brch    =>  plt_morph%LeafAreaNode_brch     , &
    LeafAreaLive_brch    =>  plt_morph%LeafAreaLive_brch     , &
    PotentialSeedSites_brch    =>  plt_morph%PotentialSeedSites_brch     , &
    SeedNumberSet_brch    =>  plt_morph%SeedNumberSet_brch     , &
    CanopyLeafAreaByLayer_pft    =>  plt_morph%CanopyLeafAreaByLayer_pft       &
  )
!     SolarNoonHour_col=hour of solar noon
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!     iYearCurrent=current year
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!     PPX,PP=PFT population per m2,grid cell
!     FracRadPARbyCanopy_pft=fraction of radiation received by each PFT canopy
!     VHeatCapCanP=canopy heat capacity
!
  IF(J.EQ.INT(SolarNoonHour_col).AND.(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantMorphologyType_pft(NZ)))) &
    .AND.(I.NE.iDayPlanting_pft(NZ) &
    .OR.iYearCurrent.NE.iYearPlanting_pft(NZ)))THEN
    IF(ITILL.LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.iDayPlanting_pft(NZ).OR.iYearCurrent.GT.iYearPlanting_pft(NZ))THEN
        XHVST=XCORP
        PPX(NZ)=PPX(NZ)*XHVST
        PlantPopulation_pft(NZ)=PlantPopulation_pft(NZ)*XHVST
        FracRadPARbyCanopy_pft(NZ)=FracRadPARbyCanopy_pft(NZ)*XHVST
        VHeatCapCanP(NZ)=VHeatCapCanP(NZ)*XHVST
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
!     LITTERFALL FROM BRANCHES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
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
!     iPlantPhenologyType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
            XHVST1=1._r8-XHVST
            D6380: DO M=1,jsken
              LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ) &
                +XHVST1*(CFOPE(ielmc,inonstruct,M,NZ)*(NonstructElmnt_brch(ielmc,NB,NZ) &
                +NoduleNonstructElmnt_brch(ielmc,NB,NZ) &
                +ShootNonstructC_brch(NB,NZ)+ReserveElmnts_brch(ielmc,NB,NZ)) &
                +CFOPE(ielmc,ifoliar,M,NZ)*(LeafChemElmnts_brch(ielmc,NB,NZ)*FWODLE(ielmc,k_fine_litr) &
                +CanopyNoduleChemElmnt_brch(ielmc,NB,NZ)) &
                +CFOPE(ielmc,inonfoliar,M,NZ)*(PetoleChemElmnt_brch(ielmc,NB,NZ)*FWODBE(ielmc,k_fine_litr) &
                +HuskChemElmnts_brch(ielmc,NB,NZ)+EarChemElmnts_brch(ielmc,NB,NZ)))

              DO NE=2,NumPlantChemElmnts
                LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *(CFOPE(NE,inonstruct,M,NZ)*(NonstructElmnt_brch(NE,NB,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ)&
                  +ReserveElmnts_brch(NE,NB,NZ)) &
                  +CFOPE(NE,ifoliar,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
                  +CanopyNoduleChemElmnt_brch(NE,NB,NZ)) &
                  +CFOPE(NE,inonfoliar,M,NZ)*(PetoleChemElmnt_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
                  +HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ)))
              ENDDO
            ENDDO D6380

            DO M=1,jsken
              DO NE=1,NumPlantChemElmnts
                LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,icwood,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
                  +PetoleChemElmnt_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

                IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenologyType_pft(NZ).NE.0)THEN
                  NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)+XHVST1*CFOPE(NE,inonfoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
                ELSE
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                    *CFOPE(NE,inonfoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
                ENDIF
                LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,icwood,M,NZ)*StalkChemElmnts_brch(NE,NB,NZ)*FWOODE(NE,k_woody_litr)

                LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,istalk,M,NZ)*StalkChemElmnts_brch(NE,NB,NZ)*FWOODE(NE,k_fine_litr)
              ENDDO
            ENDDO
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!

            ShootNonstructC_brch(NB,NZ)=ShootNonstructC_brch(NB,NZ)*XHVST
            StalkBiomassC_brch(NB,NZ)=StalkBiomassC_brch(NB,NZ)*XHVST
            DO NE=1,NumPlantChemElmnts
              NonstructElmnt_brch(NE,NB,NZ)=NonstructElmnt_brch(NE,NB,NZ)*XHVST
              NoduleNonstructElmnt_brch(NE,NB,NZ)=NoduleNonstructElmnt_brch(NE,NB,NZ)*XHVST
              ShootChemElmnt_brch(NE,NB,NZ)=ShootChemElmnt_brch(NE,NB,NZ)*XHVST
              ReserveElmnts_brch(NE,NB,NZ)=ReserveElmnts_brch(NE,NB,NZ)*XHVST
              HuskChemElmnts_brch(NE,NB,NZ)=HuskChemElmnts_brch(NE,NB,NZ)*XHVST
              EarChemElmnts_brch(NE,NB,NZ)=EarChemElmnts_brch(NE,NB,NZ)*XHVST
              GrainChemElmnts_brch(NE,NB,NZ)=GrainChemElmnts_brch(NE,NB,NZ)*XHVST
              LeafChemElmnts_brch(NE,NB,NZ)=LeafChemElmnts_brch(NE,NB,NZ)*XHVST
              CanopyNoduleChemElmnt_brch(NE,NB,NZ)=CanopyNoduleChemElmnt_brch(NE,NB,NZ)*XHVST
              PetoleChemElmnt_brch(NE,NB,NZ)=PetoleChemElmnt_brch(NE,NB,NZ)*XHVST
              StalkChemElmnts_brch(NE,NB,NZ)=StalkChemElmnts_brch(NE,NB,NZ)*XHVST
              BranchStalkChemElmnts_pft(NE,NB,NZ)=BranchStalkChemElmnts_pft(NE,NB,NZ)*XHVST
            ENDDO

            PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)*XHVST
            SeedNumberSet_brch(NB,NZ)=SeedNumberSet_brch(NB,NZ)*XHVST
            GrainSeedBiomCMean_brch(NB,NZ)=GrainSeedBiomCMean_brch(NB,NZ)*XHVST
            LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)*XHVST
            LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafChemElmnts_brch(ielmc,NB,NZ)+PetoleChemElmnt_brch(ielmc,NB,NZ))
            CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)

            CanopyStalkC_pft(NZ)=CanopyStalkC_pft(NZ)+StalkBiomassC_brch(NB,NZ)
            D8970: DO K=0,MaxNodesPerBranch1
              IF(K.NE.0)THEN
                CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)*XHVST
                CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)*XHVST
                CMassCO2BundleSheath_node(K,NB,NZ)=CMassCO2BundleSheath_node(K,NB,NZ)*XHVST
                CMassHCO3BundleSheath_node(K,NB,NZ)=CMassHCO3BundleSheath_node(K,NB,NZ)*XHVST
              ENDIF
              LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)*XHVST

              LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*XHVST
!     PetioleLengthNode_brch(K,NB,NZ)=PetioleLengthNode_brch(K,NB,NZ)*XHVST

              PetioleProteinCNode_brch(K,NB,NZ)=PetioleProteinCNode_brch(K,NB,NZ)*XHVST
!     InternodeHeightLive_brch(K,NB,NZ)=InternodeHeightLive_brch(K,NB,NZ)*XHVST
!     InternodeHeightDying_brch(K,NB,NZ)=InternodeHeightDying_brch(K,NB,NZ)*XHVST
              DO NE=1,NumPlantChemElmnts
                InternodeChemElmnt_brch(NE,K,NB,NZ)=InternodeChemElmnt_brch(NE,K,NB,NZ)*XHVST
                LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)*XHVST
                PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)*XHVST
                DO L=1,NumOfCanopyLayers1
                  LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)=LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)*XHVST
                ENDDO
              ENDDO
              D8965: DO L=1,NumOfCanopyLayers1
                CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=CanopyLeafAreaByLayer_pft(L,K,NB,NZ)*XHVST
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
!     LITTERFALL FROM ROOTS DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!
        D8980: DO L=NU,MaxNumRootLays
          D8985: DO N=1,MY(NZ)

            D6385: DO M=1,jsken
                DO NE=1,NumPlantChemElmnts
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,inonstruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
                ENDDO

              DO NR=1,NumRootAxes_pft(NZ)
                DO NE=1,NumPlantChemElmnts
                  LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,icwood,M,NZ)*(Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ) &
                    +Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,iroot,M,NZ)*(Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ) &
                    +Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
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
                *(trcg_rootml_vr(NTG,N,L,NZ)+trcs_rootml_vr(NTG,N,L,NZ))
              trcg_rootml_vr(NTG,N,L,NZ)=XHVST*trcg_rootml_vr(NTG,N,L,NZ)
              trcs_rootml_vr(NTG,N,L,NZ)=XHVST*trcs_rootml_vr(NTG,N,L,NZ)
            ENDDO
!
!     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     PrimRootLen,SecndRootLen=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootStructBiomC_vr, PopuPlantRootC_vr=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,SecndRootXNum_pvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_vr,RootVolume_vr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_vr=root surface area per plant
!     RootRespPotential_vr,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D8960: DO NR=1,NumRootAxes_pft(NZ)
              DO NE=1,NumPlantChemElmnts
                Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ)=Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ)*XHVST
                Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)=Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)*XHVST
                Root1stChemElmnt(NE,N,NR,NZ)=Root1stChemElmnt(NE,N,NR,NZ)*XHVST
              ENDDO
              PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)*XHVST
              SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)*XHVST
              SecndRootXNum_rpvr(N,L,NR,NZ)=SecndRootXNum_rpvr(N,L,NR,NZ)*XHVST
            ENDDO D8960
            DO NE=1,NumPlantChemElmnts
               RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)*XHVST
            ENDDO
            RootStructBiomC_vr(N,L,NZ)=RootStructBiomC_vr(N,L,NZ)*XHVST
             PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)*XHVST
            RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)*XHVST
            PrimRootXNumL_pvr(N,L,NZ)=PrimRootXNumL_pvr(N,L,NZ)*XHVST
            SecndRootXNum_pvr(N,L,NZ)=SecndRootXNum_pvr(N,L,NZ)*XHVST
            RootLenPerPlant_pvr(N,L,NZ)=RootLenPerPlant_pvr(N,L,NZ)*XHVST
            RootLenDensPerPlant_pvr(N,L,NZ)=RootLenDensPerPlant_pvr(N,L,NZ)*XHVST
            RootVolume_vr(N,L,NZ)=RootVolume_vr(N,L,NZ)*XHVST
            RootVH2O_vr(N,L,NZ)=RootVH2O_vr(N,L,NZ)*XHVST
            RootAreaPerPlant_vr(N,L,NZ)=RootAreaPerPlant_vr(N,L,NZ)*XHVST
            RootRespPotential_vr(N,L,NZ)=RootRespPotential_vr(N,L,NZ)*XHVST
            RCO2N(N,L,NZ)=RCO2N(N,L,NZ)*XHVST
            RCO2A(N,L,NZ)=RCO2A(N,L,NZ)*XHVST
!
!     LITTERFALL AND STATE VARIABLES FOR NODULES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(is_plant_N2fix(iPlantNfixType(NZ)).AND.N.EQ.ipltroot)THEN
              DO NE=1,NumPlantChemElmnts
                D6395: DO M=1,jsken
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+&
                    XHVST1*(CFOPE(NE,iroot,M,NZ)*RootNodueChemElmnt_pvr(NE,L,NZ) &
                    +CFOPE(NE,inonstruct,M,NZ)*RootNoduleNonstructElmnt_vr(NE,L,NZ))
                ENDDO D6395
                RootNodueChemElmnt_pvr(NE,L,NZ)=RootNodueChemElmnt_pvr(NE,L,NZ)*XHVST
                RootNoduleNonstructElmnt_vr(NE,L,NZ)=RootNoduleNonstructElmnt_vr(NE,L,NZ)*XHVST
              ENDDO
            ENDIF
          ENDDO D8985
        ENDDO D8980
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
!     DURING TILLAGE
!
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        DO NE=1,NumPlantChemElmnts
          D6400: DO M=1,jsken
            LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=&
              LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
              +(XHVST1*CFOPE(NE,inonstruct,M,NZ)*NonstructalElmnts_pft(NE,NZ))*FWOODE(NE,k_woody_litr)

            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)=&
              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
              +(XHVST1*CFOPE(NE,inonstruct,M,NZ)*NonstructalElmnts_pft(NE,NZ))*FWOODE(NE,k_fine_litr)
          ENDDO D6400
          NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)*XHVST
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByHarvest(I,J,NZ,ShootNonstructC_brch)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(inout) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,K,M,NR,N,NB,NBX,NE
  real(r8):: ZPOOLG,ZPOLNG,ZPOOLX
  real(r8) :: ZPOLNX,XHVST(NumPlantChemElmnts)
  real(r8) :: XHVST1(NumPlantChemElmnts)
  REAL(R8) :: WGLFBL(NumOfCanopyLayers1,JP1,JP1)
  real(r8) :: FHVSHK(0:MaxNodesPerBranch1),FHVSETK(0:MaxNodesPerBranch1)
  real(r8) :: ARLFY,ARLFR,ARLFG
  real(r8) :: APSILT
  real(r8) :: CPOOLX
  real(r8) :: CCPOLX
  real(r8) :: CPOOLG
  real(r8) :: CPOLNG
  real(r8) :: CPOLNX
  real(r8) :: CCPLNX
  real(r8) :: FHGT
  real(r8) :: FHVSH
  real(r8) :: FHVST4
  real(r8) :: FHGTK
  real(r8) :: FHVSETS
  real(r8) :: FHVSETG
  real(r8) :: FHVSHG
  real(r8) :: FHVSETH
  real(r8) :: FHVSETE
  real(r8) :: FHVSHH
  real(r8) :: FHVSHE
  real(r8) :: FDM
  real(r8) :: FFIRE(NumPlantChemElmnts)
  real(r8) :: FHVSE(NumPlantChemElmnts)
  real(r8) :: HTSTKX
  real(r8) :: PPOOLG
  real(r8) :: PPOLNG,PPOOLX,PPOLNX
  real(r8) :: VOLWPX
  real(r8) :: WHVSBL
  real(r8) :: WTSTKT
  real(r8) :: WTLSBX
  real(r8) :: TotPhytomassRemoval,WHVSLF,WHVHSH,WHVEAH,WHVGRH,WHVSCP
  real(r8) :: WHVSTH,WHVRVH,WHVSLX,WHVSLY,WHVSCL,WHVSNL,WHVXXX
  real(r8) :: WHVSSX,WTSHTT,WHVSHX,WHVSHY,WHVSHH,WHVSCS,WHVSNS
  real(r8) :: WHVHSX,WHVHSY,WHVEAX,WHVEAY,WHVGRX,WHVGRY,WHVSNP
  real(r8) :: WHVSKX,WHVSTX,WHVSTY,WHVRVX,WHVRVY,WTNDG,WTNDNG
  real(r8) :: WTNDPG,WGLFGX,WGSHGX,WGLFGY,WGSHGY
  real(r8) :: WHVSBS,WHVSCX,WHVSNX,WVPLT
  real(r8) :: FHVSH1,FHVSHT
  real(r8) :: WGLFGE(NumPlantChemElmnts)
  integer :: NTG
!     begin_execution
  associate(                             &
    HVST     =>  plt_distb%HVST    , &
    EHVST    =>  plt_distb%EHVST   , &
    DCORP    =>  plt_distb%DCORP   , &
    THIN_pft     =>  plt_distb%THIN_pft    , &
    ITILL    =>  plt_distb%ITILL   , &
    iHarvstType_pft    =>  plt_distb%iHarvstType_pft   , &
    jHarvst_pft    =>  plt_distb%jHarvst_pft   , &
    PO4byFire_pft    =>  plt_distb%PO4byFire_pft   , &
    N2ObyFire_pft    =>  plt_distb%N2ObyFire_pft   , &
    NH3byFire_pft    =>  plt_distb%NH3byFire_pft   , &
    VOXYF    =>  plt_distb%VOXYF   , &
    CH4ByFire_pft    =>  plt_distb%CH4ByFire_pft   , &
    CO2ByFire_pft    =>  plt_distb%CO2ByFire_pft   , &
    UVOLO    =>  plt_ew%UVOLO      , &
    CanopyWater_pft    =>  plt_ew%CanopyWater_pft      , &
    PSICanopy_pft   =>  plt_ew%PSICanopy_pft     , &
    CMassHCO3BundleSheath_node     =>  plt_photo%CMassHCO3BundleSheath_node    , &
    CMassCO2BundleSheath_node     =>  plt_photo%CMassCO2BundleSheath_node    , &
    CPOOL3   =>  plt_photo%CPOOL3  , &
    CPOOL4   =>  plt_photo%CPOOL4  , &
    PlantPopulation_pft       =>  plt_site%PlantPopulation_pft       , &
    PPI      =>  plt_site%PPI      , &
    PPX      =>  plt_site%PPX      , &
    NU       =>  plt_site%NU       , &
    MaxNumRootLays       => plt_site%MaxNumRootLays        , &
    SolarNoonHour_col   => plt_site%SolarNoonHour_col    , &
    ZEROS    => plt_site%ZEROS     , &
    ZERO     => plt_site%ZERO      , &
    AREA3    => plt_site%AREA3     , &
    VOLWOU   => plt_site%VOLWOU    , &
    RootMycoNonstructElmnt_vr   => plt_biom%RootMycoNonstructElmnt_vr    , &
    RootProteinC_pvr    => plt_biom%RootProteinC_pvr     , &
    PopuPlantRootC_vr    => plt_biom% PopuPlantRootC_vr     , &
    RootStructBiomC_vr   => plt_biom%RootStructBiomC_vr    , &
    NonstructalElmnts_pft    => plt_biom%NonstructalElmnts_pft     , &
    CanopyStalkC_pft    => plt_biom%CanopyStalkC_pft     , &
    CanopyLeafShethC_pft     => plt_biom%CanopyLeafShethC_pft      , &
    StalkBiomassC_brch   => plt_biom%StalkBiomassC_brch    , &
    CanopyNoduleChemElmnt_brch   => plt_biom%CanopyNoduleChemElmnt_brch    , &
    ReserveChemElmnts_pft   => plt_biom%ReserveChemElmnts_pft    , &
    GrainChemElmnts_brch   => plt_biom%GrainChemElmnts_brch    , &
    StalkChemElmnts_brch  => plt_biom%StalkChemElmnts_brch   , &
    ShootChemElmnt_brch  => plt_biom%ShootChemElmnt_brch   , &
    HuskChemElmnts_brch  => plt_biom%HuskChemElmnts_brch   , &
    EarChemElmnts_brch  => plt_biom%EarChemElmnts_brch   , &
    InternodeChemElmnt_brch   => plt_biom%InternodeChemElmnt_brch    , &
    LeafPetolBiomassC_brch  => plt_biom%LeafPetolBiomassC_brch     , &
    ReserveElmnts_brch  => plt_biom%ReserveElmnts_brch   , &
    BranchStalkChemElmnts_pft  => plt_biom%BranchStalkChemElmnts_pft   , &
    NoduleNonstructElmnt_brch   => plt_biom%NoduleNonstructElmnt_brch     , &
    LeafChemElmntByLayer_pft   => plt_biom%LeafChemElmntByLayer_pft    , &
    PetoleChemElmnt_brch => plt_biom%PetoleChemElmnt_brch  , &
    NonstructElmnt_brch   => plt_biom%NonstructElmnt_brch    , &
    PetioleElmntNode_brch   => plt_biom%PetioleElmntNode_brch    , &
    LeafProteinCNode_brch     => plt_biom%LeafProteinCNode_brch      , &
    LeafElmntNode_brch    => plt_biom%LeafElmntNode_brch     , &
    LeafChemElmnts_brch  => plt_biom%LeafChemElmnts_brch   , &
    PetioleProteinCNode_brch   => plt_biom%PetioleProteinCNode_brch    , &
    StalkChemElmnts_pft   => plt_biom%StalkChemElmnts_pft    , &
    CanopyNonstructElementConc_pft   => plt_biom%CanopyNonstructElementConc_pft    , &
    NoduleNonstructCconc_pft   => plt_biom%NoduleNonstructCconc_pft    , &
    LeafChemElmnts_pft    => plt_biom%LeafChemElmnts_pft     , &
    GrainChemElmnts_pft    => plt_biom%GrainChemElmnts_pft     , &
    ShootChemElmnts_pft   => plt_biom%ShootChemElmnts_pft     , &
    HuskChemElmnts_pft   => plt_biom%HuskChemElmnts_pft    , &
    EarChemElmnts_pft   => plt_biom%EarChemElmnts_pft    , &
    PetioleChemElmnts_pft   => plt_biom%PetioleChemElmnts_pft    , &
    AvgCanopyBiomC2Graze_pft   => plt_biom%AvgCanopyBiomC2Graze_pft    , &
    Root1stStructChemElmnt_pvr   => plt_biom%Root1stStructChemElmnt_pvr    , &
    Root1stChemElmnt   => plt_biom%Root1stChemElmnt    , &
    Root2ndStructChemElmnt_pvr   => plt_biom%Root2ndStructChemElmnt_pvr    , &
    RootNodueChemElmnt_pvr   => plt_biom%RootNodueChemElmnt_pvr    , &
    RootNoduleNonstructElmnt_vr  => plt_biom%RootNoduleNonstructElmnt_vr   , &
    ZEROP    => plt_biom%ZEROP     , &
    ZEROL    => plt_biom%ZEROL     , &
    CanopyLeafCpft_lyr    => plt_biom%CanopyLeafCpft_lyr     , &
    FracHour4LeafoffRemob     => plt_allom%FracHour4LeafoffRemob     , &
    FWODRE   => plt_allom%FWODRE   , &
    FWOODE   => plt_allom%FWOODE   , &
    FWODBE   => plt_allom%FWODBE   , &
    FWODLE   => plt_allom%FWODLE   , &
    GrainSeedBiomCMean_brch    => plt_allom%GrainSeedBiomCMean_brch    , &
    iPlantBranchState_brch    =>  plt_pheno%iPlantBranchState_brch   , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP    , &
    iPlantCalendar_brch   =>  plt_pheno%iPlantCalendar_brch  , &
    MatureGroup_brch   =>  plt_pheno%MatureGroup_brch  , &
    LeafNumberAtFloralInit_brch   =>  plt_pheno%LeafNumberAtFloralInit_brch  , &
    iPlantMorphologyType_pft   =>  plt_pheno%iPlantMorphologyType_pft  , &
    iPlantTurnoverPattern_pft   =>  plt_pheno%iPlantTurnoverPattern_pft  , &
    doInitLeafOut_brch    =>  plt_pheno%doInitLeafOut_brch   , &
    iPlantPhenologyPattern_pft   =>  plt_pheno%iPlantPhenologyPattern_pft  , &
    HourThreshold4LeafOff_brch    =>  plt_pheno%HourThreshold4LeafOff_brch   , &
    Hours4LeafOff_brch     =>  plt_pheno%Hours4LeafOff_brch    , &
    iPlantPhenologyType_pft   =>  plt_pheno%iPlantPhenologyType_pft  , &
    TotalNodeNumNormByMatgrp_brch   =>  plt_pheno%TotalNodeNumNormByMatgrp_brch  , &
    TotReproNodeNumNormByMatrgrp_brch   =>  plt_pheno%TotReproNodeNumNormByMatrgrp_brch  , &
    HourFailGrainFill_brch     =>  plt_pheno%HourFailGrainFill_brch    , &
    MatureGroup_pft  =>  plt_pheno%MatureGroup_pft , &
    CORGC    =>  plt_soilchem%CORGC, &
    THETW    =>  plt_soilchem%THETW, &
    CFOPE    =>  plt_soilchem%CFOPE, &
    inonstruct =>  pltpar%inonstruct  , &
    ifoliar  =>  pltpar%ifoliar   , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    inonfoliar =>  pltpar%inonfoliar  , &
    k_fine_litr=> pltpar%k_fine_litr ,&
    k_woody_litr=> pltpar%k_woody_litr,&
    icwood   =>  pltpar%icwood    , &
    trcg_rootml_vr     =>  plt_rbgc%trcg_rootml_vr , &
    trcs_rootml_vr     =>  plt_rbgc%trcs_rootml_vr , &
    RootGasLossDisturb_pft    =>   plt_bgcr%RootGasLossDisturb_pft    , &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr     , &
    Eco_NBP_col     =>  plt_bgcr%Eco_NBP_col     , &
    CO2NetFix_pft     =>  plt_bgcr%CO2NetFix_pft     , &
    RCO2A    =>  plt_rbgc%RCO2A    , &
    RootRespPotential_vr    =>  plt_rbgc%RootRespPotential_vr    , &
    RCO2N    =>  plt_rbgc%RCO2N    , &
    SecndRootXNum_pvr     =>  plt_morph%SecndRootXNum_pvr    , &
    SecndRootXNum_rpvr    =>  plt_morph%SecndRootXNum_rpvr   , &
    RootLenPerPlant_pvr    =>  plt_morph%RootLenPerPlant_pvr   , &
    RootAreaPerPlant_vr    =>  plt_morph%RootAreaPerPlant_vr   , &
    RootVolume_vr    =>  plt_morph%RootVolume_vr   , &
    NGTopRootLayer_pft      =>  plt_morph%NGTopRootLayer_pft     , &
    MY       =>  plt_morph%MY      , &
    CanopyHeight_pft      =>  plt_morph%CanopyHeight_pft     , &
    RootVH2O_vr   =>  plt_morph%RootVH2O_vr  , &
    RootLenDensPerPlant_pvr    =>  plt_morph%RootLenDensPerPlant_pvr   , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType  , &
    CanopyLAgrid_lyr    =>  plt_morph%CanopyLAgrid_lyr   , &
    CanopyHeightz_col       =>  plt_morph%CanopyHeightz_col      , &
    LeafAreaLive_brch    =>  plt_morph%LeafAreaLive_brch   , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft     , &
    CanopyStemA_pft    =>  plt_morph%CanopyStemA_pft   , &
    InternodeHeightDying_brch   =>  plt_morph%InternodeHeightDying_brch  , &
    PrimRootXNumL_pvr    =>  plt_morph%PrimRootXNumL_pvr   , &
    SecndRootLen    =>  plt_morph%SecndRootLen   , &
    PrimRootLen    =>  plt_morph%PrimRootLen   , &
    InternodeHeightLive_brch   =>  plt_morph%InternodeHeightLive_brch  , &
    PotentialSeedSites_brch    =>  plt_morph%PotentialSeedSites_brch   , &
    SeedNumberSet_brch    =>  plt_morph%SeedNumberSet_brch   , &
    PetioleLengthNode_brch    =>  plt_morph%PetioleLengthNode_brch   , &
    LeafAreaNode_brch    =>  plt_morph%LeafAreaNode_brch   , &
    CanopyLeafApft_lyr    =>  plt_morph%CanopyLeafApft_lyr   , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr   , &
    CanopyLeafAreaByLayer_pft    =>  plt_morph%CanopyLeafAreaByLayer_pft   , &
    CanopyBranchStemApft_lyr    =>  plt_morph%CanopyBranchStemApft_lyr   , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft    , &
    NodeNumberAtAnthesis_brch    =>  plt_morph%NodeNumberAtAnthesis_brch   , &
    NumOfMainBranch_pft      =>  plt_morph%NumOfMainBranch_pft     , &
    NodeNumberToInitFloral_brch    =>  plt_morph%NodeNumberToInitFloral_brch   , &
    ShootNodeNumber_brch    =>  plt_morph%ShootNodeNumber_brch   , &
    ClumpFactor      =>  plt_morph%ClumpFactor     , &
    CanopyLeafArea_grd    =>  plt_morph%CanopyLeafArea_grd   , &
    iPlantPhotosynthesisType   =>  plt_photo%iPlantPhotosynthesisType    &
  )
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
  IF((iHarvstType_pft(NZ).GE.iharvtyp_none.AND.J.EQ.INT(SolarNoonHour_col) &
    .AND.iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
    .OR.(iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo))THEN
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
!     CanopyLeafArea_grd,CanopyLAgrid_lyr=leaf area of combined canopy, canopy layer
!     ARLFR,ARLFY=leaf area harvested,remaining
!     ZL=height to bottom of each canopy layer
!
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(jHarvst_pft(NZ).NE.jharvtyp_tmareseed)THEN
        !terminate and reseed
        PPX(NZ)=PPX(NZ)*(1._r8-THIN_pft(NZ))
        PlantPopulation_pft(NZ)=PlantPopulation_pft(NZ)*(1._r8-THIN_pft(NZ))
      ELSE
!     PPI(NZ)=AMAX1(1.0_r8,0.5_r8*(PPI(NZ)+CanopySeedNumber_pft(NZ)/AREA3(NU)))
        PPX(NZ)=PPI(NZ)
        PlantPopulation_pft(NZ)=PPX(NZ)*AREA3(NU)
      ENDIF
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
        ClumpFactor(NZ)=ClumpFactor(NZ)*HVST(NZ)
      ENDIF
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.HVST(NZ).LT.0.0)THEN
        ARLFY=(1._r8-ABS(HVST(NZ)))*CanopyLeafArea_grd
        ARLFR=0._r8
        D9875: DO L=1,NumOfCanopyLayers1
          IF(CanopyHeightz_col(L).GT.CanopyHeightz_col(L-1).AND.&
            CanopyLAgrid_lyr(L).GT.ZEROS.AND.ARLFR.LT.ARLFY)THEN
            IF(ARLFR+CanopyLAgrid_lyr(L).GT.ARLFY)THEN
              HVST(NZ)=CanopyHeightz_col(L-1)+((ARLFY-ARLFR)/CanopyLAgrid_lyr(L))&
                *(CanopyHeightz_col(L)-CanopyHeightz_col(L-1))
            ENDIF
          ELSE
            HVST(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+CanopyLAgrid_lyr(L)
        ENDDO D9875
      ENDIF
      TotPhytomassRemoval=0._r8
      WHVSLF=0._r8
      WHVHSH=0._r8
      WHVEAH=0._r8
      WHVGRH=0._r8
      WHVSCP=0._r8
      WHVSTH=0._r8
      WHVRVH=0._r8
    ELSE
!
!     GRAZING REMOVAL
!
!     AvgCanopyBiomC2Graze_pft=average biomass in landscape grazing section
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     TotPhytomassRemoval=total phytomass grazed, removed
!     fTgrowCanP=temperature function for canopy growth
!     CCPOLP=nonstructural C concentration in canopy
!     NoduleNonstructCconc_pft=nonstructural C concentration in canopy nodules
!
      IF(AvgCanopyBiomC2Graze_pft(NZ).GT.ZEROP(NZ))THEN
        TotPhytomassRemoval=HVST(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8 &
          *AREA3(NU)*ShootChemElmnts_pft(ielmc,NZ)/AvgCanopyBiomC2Graze_pft(NZ)
      ELSE
        TotPhytomassRemoval=0._r8
      ENDIF
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
        TotPhytomassRemoval=TotPhytomassRemoval*fTgrowCanP(NZ)
      ENDIF
      CCPOLX=CanopyNonstructElementConc_pft(ielmc,NZ)/(1.0_r8+CanopyNonstructElementConc_pft(ielmc,NZ))
      CCPLNX=NoduleNonstructCconc_pft(NZ)/(1.0_r8+NoduleNonstructCconc_pft(NZ))
!
!     LEAF,BACTERIA GRAZED,REMOVED
!
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosyst
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WTLF=PFT leaf C mass
!     WHVXXX=grazing requirement unmet by leaf
!
      WHVSLX=TotPhytomassRemoval*EHVST(1,iplthvst_leaf,NZ)
      WHVSLY=AMIN1(LeafChemElmnts_pft(ielmc,NZ),WHVSLX)
      WHVSLF=WHVSLY*(1._r8-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AZMAX1(WHVSLX-WHVSLY)
      WHVSSX=TotPhytomassRemoval*EHVST(1,iplthvst_finenonleaf,NZ)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     WHVXXX=grazing requirement unmet by non-foliar removal
!
      WTSHTT=PetioleChemElmnts_pft(ielmc,NZ)+HuskChemElmnts_pft(ielmc,NZ)+EarChemElmnts_pft(ielmc,NZ)+&
        GrainChemElmnts_pft(ielmc,NZ)
      IF(WTSHTT.GT.ZEROP(NZ))THEN
        WHVSHX=WHVSSX*PetioleChemElmnts_pft(ielmc,NZ)/WTSHTT+WHVXXX
        WHVSHY=AMIN1(PetioleChemElmnts_pft(ielmc,NZ),WHVSHX)
        WHVSHH=WHVSHY*(1._r8-CCPOLX)
        WHVSCS=WHVSHY*CCPOLX
        WHVSNS=WHVSHY*CCPLNX
        WHVXXX=AZMAX1(WHVSHX-WHVSHY)
        WHVHSX=WHVSSX*HuskChemElmnts_pft(ielmc,NZ)/WTSHTT+WHVXXX
        WHVHSY=AMIN1(HuskChemElmnts_pft(ielmc,NZ),WHVHSX)
        WHVHSH=WHVHSY
        WHVXXX=AZMAX1(WHVHSX-WHVHSY)
        WHVEAX=WHVSSX*EarChemElmnts_pft(ielmc,NZ)/WTSHTT+WHVXXX
        WHVEAY=AMIN1(EarChemElmnts_pft(ielmc,NZ),WHVEAX)
        WHVEAH=WHVEAY
        WHVXXX=AZMAX1(WHVEAX-WHVEAY)
        WHVGRX=WHVSSX*GrainChemElmnts_pft(ielmc,NZ)/WTSHTT+WHVXXX
        WHVGRY=AMIN1(GrainChemElmnts_pft(ielmc,NZ),WHVGRX)
        WHVGRH=WHVGRY
        WHVXXX=AZMAX1(WHVGRX-WHVGRY)
      ELSE
        WHVSHH=0._r8
        WHVSCS=0._r8
        WHVSNS=0._r8
        WHVHSH=0._r8
        WHVEAH=0._r8
        WHVGRH=0._r8
        WHVXXX=WHVXXX+WHVSSX
      ENDIF
      WHVSCP=WHVSCL+WHVSCS
      WHVSNP=WHVSNL+WHVSNS
      WHVSKX=TotPhytomassRemoval*EHVST(1,iplthvst_woody,NZ)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     WHVXXX=grazing requirement unmet by stalk,reserve
!
      WTSTKT=StalkChemElmnts_pft(ielmc,NZ)+ReserveChemElmnts_pft(ielmc,NZ)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
        WHVSTX=WHVSKX*StalkChemElmnts_pft(ielmc,NZ)/WTSTKT+WHVXXX
        WHVSTY=AMIN1(StalkChemElmnts_pft(ielmc,NZ),WHVSTX)
        WHVSTH=WHVSTY
        WHVXXX=AZMAX1(WHVSTX-WHVSTY)
        WHVRVX=WHVSKX*ReserveChemElmnts_pft(ielmc,NZ)/WTSTKT+WHVXXX
        WHVRVY=AMIN1(ReserveChemElmnts_pft(ielmc,NZ),WHVRVX)
        WHVRVH=WHVRVY
        WHVXXX=AZMAX1(WHVRVX-WHVRVY)
      ELSE
        WHVSTH=0._r8
        WHVRVH=0._r8
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
          WHVSLY=AMIN1(LeafChemElmnts_pft(ielmc,NZ)-WHVSLF-WHVSCL,WHVXXX)
          WHVSLF=WHVSLF+WHVSLY*(1._r8-CCPOLX)
          WHVSCL=WHVSCL+WHVSLY*CCPOLX
          WHVSNL=WHVSNL+WHVSLY*CCPLNX
          WHVXXX=AZMAX1(WHVXXX-WHVSLY)
          IF(WTSHTT.GT.ZEROP(NZ))THEN
            WHVSHX=WHVXXX*PetioleChemElmnts_pft(ielmc,NZ)/WTSHTT
            WHVSHY=AMIN1(PetioleChemElmnts_pft(ielmc,NZ),WHVSHX)
            WHVSHH=WHVSHH+WHVSHY*(1._r8-CCPOLX)
            WHVSCS=WHVSCS+WHVSHY*CCPOLX
            WHVSNS=WHVSNS+WHVSHY*CCPLNX
            WHVXXX=AZMAX1(WHVXXX-WHVSHY)
            WHVHSX=WHVXXX*HuskChemElmnts_pft(ielmc,NZ)/WTSHTT
            WHVHSY=AMIN1(HuskChemElmnts_pft(ielmc,NZ),WHVHSX)
            WHVHSH=WHVHSH+WHVHSY
            WHVXXX=AZMAX1(WHVXXX-WHVHSY)
            WHVEAX=WHVXXX*EarChemElmnts_pft(ielmc,NZ)/WTSHTT
            WHVEAY=AMIN1(EarChemElmnts_pft(ielmc,NZ),WHVEAX)
            WHVEAH=WHVEAH+WHVEAY
            WHVXXX=AZMAX1(WHVEAX-WHVEAY)
            WHVGRX=WHVXXX*GrainChemElmnts_pft(ielmc,NZ)/WTSHTT
            WHVGRY=AMIN1(GrainChemElmnts_pft(ielmc,NZ),WHVGRX)
            WHVGRH=WHVGRH+WHVGRY
            WHVXXX=AZMAX1(WHVGRX-WHVGRY)
          ENDIF
        ENDIF
      ENDIF
!
!     ALL HARVEST REMOVALS
!
!     WGLFBL=branch leaf C mass in canopy layer
!
      D9860: DO NB=1,NumOfBranches_pft(NZ)
        DO  L=1,NumOfCanopyLayers1
          DO  K=0,MaxNodesPerBranch1
            WGLFBL(L,NB,NZ)=0._r8
          enddo
        enddo
      ENDDO D9860
      D9870: DO NB=1,NumOfBranches_pft(NZ)
        DO  L=1,NumOfCanopyLayers1
          DO  K=0,MaxNodesPerBranch1
            WGLFBL(L,NB,NZ)=WGLFBL(L,NB,NZ)+LeafChemElmntByLayer_pft(ielmc,L,K,NB,NZ)
          enddo
        enddo
      ENDDO D9870
    ENDIF
!
!     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     ZL=height to bottom of each canopy layer
!     FHGT=fraction of canopy layer height not harvested
!     FHVSE(ielmc)=fraction of canopy layer mass not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
    D9865: DO L=NumOfCanopyLayers1,1,-1
      IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
        IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
          IF(CanopyHeightz_col(L).GT.CanopyHeightz_col(L-1))THEN
            FHGT=AZMAX1(AMIN1(1.0_r8,1._r8-((CanopyHeightz_col(L))-HVST(NZ))/ &
              (CanopyHeightz_col(L)-CanopyHeightz_col(L-1))))
          ELSE
            FHGT=1.0_r8
          ENDIF
        ELSE
          FHGT=0._r8
        ENDIF
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FHVSE(ielmc)=AZMAX1(1._r8-(1._r8-FHGT)*EHVST(1,iplthvst_leaf,NZ))
          FHVSH=FHVSE(ielmc)
        ELSE
          FHVSE(ielmc)=AZMAX1(1._r8-THIN_pft(NZ))
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
            FHVSH=1.0_r8-(1._r8-FHGT)*EHVST(1,iplthvst_leaf,NZ)*THIN_pft(NZ)
          ELSE
            FHVSH=FHVSE(ielmc)
          ENDIF
        ENDIF
      ELSE
        FHVSE(ielmc)=0._r8
        FHVSH=0._r8
      ENDIF
!
!     CUT LEAVES AT HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTLF=PFT leaf C mass
!     WGLFBL=branch leaf C mass in canopy layer
!     WHVBSL,WHVSLF=layer,total leaf C mass removed
!     WGLFL=leaf node C in canopy layer
!     FHVSE(ielmc)=fraction of leaf node mass not harvested
!
      D9855: DO NB=1,NumOfBranches_pft(NZ)
        IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
          .AND.LeafChemElmnts_pft(ielmc,NZ).GT.ZEROL(NZ))THEN
          WHVSBL=WHVSLF*AZMAX1(WGLFBL(L,NB,NZ))/LeafChemElmnts_pft(ielmc,NZ)
        ELSE
          WHVSBL=0._r8
        ENDIF
        D9845: DO K=MaxNodesPerBranch1,0,-1
          IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo).OR.WHVSBL.GT.0.0_r8)THEN
            IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
              IF(LeafChemElmntByLayer_pft(ielmc,L,K,NB,NZ).GT.WHVSBL)THEN
                FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,(LeafChemElmntByLayer_pft(ielmc,L,K,NB,NZ)-WHVSBL)/LeafChemElmntByLayer_pft(ielmc,L,K,NB,NZ)))
                FHVSH=FHVSE(ielmc)
              ELSE
                FHVSE(ielmc)=1.0_r8
                FHVSH=1.0_r8
              ENDIF
            ENDIF
        !
!     HARVESTED LEAF AREA, C, N, P
!
!     FHVSE(ielmc)=fraction of leaf node mass not harvested
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanopyLeafAreaByLayer_pft,CanopyBranchStemApft_lyr=leaf,stalk node area in canopy layer
!     LeafElmntRemoval(ielmc),LeafElmntRemoval(ielmn),LeafElmntRemoval(ielmp)=harvested leaf C,N,P
!     LeafElmntHarv2Litr(ielmc),LeafElmntHarv2Litr(ielmn),LeafElmntHarv2Litr(ielmp)=harvested leaf C,N,P to litter
!     WoodyElmntRemoval(ielmc),WoodyElmntRemoval(ielmn),WoodyElmntRemoval(ielmp)=harvested woody C,N,P
!     WoodyElmntHarv2Litr(ielmc),WoodyElmntHarv2Litr(ielmn),WoodyElmntHarv2Litr(ielmp)=harvested woody C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!

            WHVSBL=WHVSBL-(1._r8-FHVSE(ielmc))*LeafChemElmntByLayer_pft(ielmc,L,K,NB,NZ)
            FHVSH1=1._r8-FHVSH
            FHVSHT=FHVSH-FHVSE(ielmc)
            DO NE=1,NumPlantChemElmnts
              LeafElmntRemoval(NE)=LeafElmntRemoval(NE)+FHVSH1*LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)*FWODLE(NE,k_fine_litr)
              LeafElmntHarv2Litr(NE)=LeafElmntHarv2Litr(NE)+FHVSHT*LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)*FWODLE(NE,k_fine_litr)
              WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE)+FHVSH1*LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)*FWODLE(NE,k_woody_litr)
              WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE)+FHVSHT*LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)*FWODLE(NE,k_woody_litr)
              LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)=FHVSE(ielmc)*LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)
            ENDDO
!
!     REMAINING LEAF C,N,P AND AREA
!
            CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=FHVSE(ielmc)*CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
            IF(K.EQ.1)THEN
              CanopyBranchStemApft_lyr(L,NB,NZ)=FHVSE(ielmc)*CanopyBranchStemApft_lyr(L,NB,NZ)
            ENDIF
          ENDIF

        ENDDO D9845
      ENDDO D9855
      CanopyLeafApft_lyr(L,NZ)=0._r8
      CanopyLeafCpft_lyr(L,NZ)=0._r8
      CanopyStemApft_lyr(L,NZ)=CanopyStemApft_lyr(L,NZ)*FHVSE(ielmc)
    ENDDO D9865

    D9835: DO NB=1,NumOfBranches_pft(NZ)
      CPOOLG=0._r8
      ZPOOLG=0._r8
      PPOOLG=0._r8
      CPOLNG=0._r8
      ZPOLNG=0._r8
      PPOLNG=0._r8
      WTNDG=0._r8
      WTNDNG=0._r8
      WTNDPG=0._r8
      WGLFGX=0._r8
      WGSHGX=0._r8
      WGLFGY=0._r8
      WGSHGY=0._r8
      D9825: DO K=0,MaxNodesPerBranch1
        ARLFG=0._r8
        WGLFGE(1:NumPlantChemElmnts)=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanopyLeafAreaByLayer_pft,CanopyLeafApft_lyr=leaf node,total area in canopy layer
!
        D9815: DO L=1,NumOfCanopyLayers1
          ARLFG=ARLFG+CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
          DO NE=1,NumPlantChemElmnts
            WGLFGE(NE)=WGLFGE(NE)+LeafChemElmntByLayer_pft(NE,L,K,NB,NZ)
          ENDDO
          CanopyLeafApft_lyr(L,NZ)=CanopyLeafApft_lyr(L,NZ)+CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
          CanopyLeafCpft_lyr(L,NZ)=CanopyLeafCpft_lyr(L,NZ)+LeafChemElmntByLayer_pft(ielmc,L,K,NB,NZ)
        ENDDO D9815
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FHVSETK=fraction of internode layer mass not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
        IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
          IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ).AND.EHVST(1,iplthvst_leaf,NZ).GT.0.0)THEN
            FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(1._r8-(1._r8-AZMAX1(WGLFGE(ielmc)) &
              /LeafElmntNode_brch(ielmc,K,NB,NZ))*EHVST(1,iplthvst_finenonleaf,NZ)/EHVST(1,iplthvst_leaf,NZ))))
            FHVSHK(K)=FHVSETK(K)
        ELSE
          IF(isclose(THIN_pft(NZ),0._r8))THEN
            FHVSETK(K)=1.0_r8-EHVST(1,iplthvst_finenonleaf,NZ)
            FHVSHK(K)=FHVSETK(K)
          ELSE
            FHVSETK(K)=1.0_r8-THIN_pft(NZ)
            IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
              FHVSHK(K)=1.0_r8-EHVST(1,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
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
      DO NE=1,NumPlantChemElmnts
        LeafChemElmnts_brch(NE,NB,NZ)=LeafChemElmnts_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ)+WGLFGE(NE)
      ENDDO
      LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-LeafAreaNode_brch(K,NB,NZ)+ARLFG
      IF(LeafAreaNode_brch(K,NB,NZ).GT.ZEROP(NZ))THEN
        LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*ARLFG/LeafAreaNode_brch(K,NB,NZ)
      ELSE
        LeafProteinCNode_brch(K,NB,NZ)=0._r8
      ENDIF
      LeafAreaNode_brch(K,NB,NZ)=ARLFG
      DO NE=1,NumPlantChemElmnts
        LeafElmntNode_brch(NE,K,NB,NZ)=WGLFGE(NE)
      ENDDO
      WGLFGX=WGLFGX+LeafElmntNode_brch(ielmc,K,NB,NZ)
    ENDDO D9825
!
!     CUT SHEATHS OR PETIOLES AND STALKS HARVESTED NODES AND LAYERS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSHE,WTSHEB=PFT,branch petiole C mass
!     WHVSBS,WHVSHH=branch, PFT petiole C mass removed
!     InternodeHeightLive_brch=internode length
!     HTSTKX=internode length removed
!
      HTSTKX=0._r8
      IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
        .AND.PetioleChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
        WHVSBS=WHVSHH*PetoleChemElmnt_brch(ielmc,NB,NZ)/PetioleChemElmnts_pft(ielmc,NZ)
      ELSE
        WHVSBS=0._r8
      ENDIF
      D9805: DO K=MaxNodesPerBranch1,0,-1
!112   FORMAT(A8,8I4,12E12.4)
        IF(InternodeHeightLive_brch(K,NB,NZ).GT.0.0) &
          HTSTKX=AMAX1(HTSTKX,InternodeHeightLive_brch(K,NB,NZ))
!
!     HARVESTED SHEATH OR PETIOLE C,N,P
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WHVSBS=branch petiole C mass removed
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
!     FHVSETK=fraction of internode layer mass not harvested
!     FineNonleafElmntRemoval(ielmc),FineNonleafElmntRemoval(ielmn),FineNonleafElmntRemoval(ielmp)=harvested petiole C,N,P
!     PetioleElmntHarv2Litr(ielmc),PetioleElmntHarv2Litr(ielmn),PetioleElmntHarv2Litr(ielmp)=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     PetioleLengthNode_brch,InternodeHeightLive_brch=petiole,internode length
!
          IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)&
            .OR.WHVSBS.GT.0.0_r8)THEN
            IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing.OR.iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
              IF(PetioleElmntNode_brch(ielmc,K,NB,NZ).GT.WHVSBS)THEN
                FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(PetioleElmntNode_brch(ielmc,K,NB,NZ)-WHVSBS)/&
                  PetioleElmntNode_brch(ielmc,K,NB,NZ)))
                FHVSHK(K)=FHVSETK(K)
              ELSE
                FHVSETK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1._r8-FHVSETK(K))*PetioleElmntNode_brch(ielmc,K,NB,NZ)
            DO NE=1,NumPlantChemElmnts
              FineNonleafElmntRemoval(NE)=FineNonleafElmntRemoval(NE)+(1._r8-FHVSHK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)
              PetioleElmntHarv2Litr(NE)=PetioleElmntHarv2Litr(NE)+(FHVSHK(K)-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)
              WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE)+(1._r8-FHVSHK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
              WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE)+(FHVSHK(K)-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
            ENDDO
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     PetioleElmntNode_brch=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     PetioleLengthNode_brch=node petiole height
!     PetioleProteinCNode_brch=petiole protein mass
!
            WGSHGY=WGSHGY+PetioleElmntNode_brch(ielmc,K,NB,NZ)
            PetioleProteinCNode_brch(K,NB,NZ)=FHVSETK(K)*PetioleProteinCNode_brch(K,NB,NZ)

            DO NE=1,NumPlantChemElmnts
              PetoleChemElmnt_brch(NE,NB,NZ)=PetoleChemElmnt_brch(NE,NB,NZ) &
                -(1._r8-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)
              PetioleElmntNode_brch(NE,K,NB,NZ)=FHVSETK(K)*PetioleElmntNode_brch(NE,K,NB,NZ)
            ENDDO
!            PetioleProteinCNode_brch(K,NB,NZ)=FHVSETK(K)*PetioleProteinCNode_brch(K,NB,NZ)
            IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.PetioleLengthNode_brch(K,NB,NZ).GT.0.0_r8)THEN
              FHGT=AZMAX1(AMIN1(1.0_r8,(InternodeHeightLive_brch(K,NB,NZ) &
                +PetioleLengthNode_brch(K,NB,NZ)-HVST(NZ))/PetioleLengthNode_brch(K,NB,NZ)))
              PetioleLengthNode_brch(K,NB,NZ)=(1._r8-FHGT)*PetioleLengthNode_brch(K,NB,NZ)
            ELSE
              PetioleLengthNode_brch(K,NB,NZ)=FHVSETK(K)*PetioleLengthNode_brch(K,NB,NZ)
            ENDIF
            WGSHGX=WGSHGX+PetioleElmntNode_brch(ielmc,K,NB,NZ)
!     IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
!     IF(InternodeHeightLive_brch(K,NB,NZ).GT.HVST(NZ)
!    2.OR.iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
!     IF(isclose(FHVSETK(K),0._r8).AND.K.GT.0)THEN
!     IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType_pft(NZ)))THEN
!     NumOfLeaves_brch(NB,NZ)=AZMAX1(NumOfLeaves_brch(NB,NZ)-1.0)
!     ELSE
!     NumOfLeaves_brch(NB,NZ)=AZMAX1(NumOfLeaves_brch(NB,NZ)-0.04)
!     ENDIF
!     ENDIF
!     ENDIF
!     ENDIF
          ENDIF
    ENDDO D9805
!
!     CUT NON-STRUCTURAL C,N,P IN HARVESTED BRANCHES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     FHVSE(ielmc)=fraction of leaf+petiole node mass not harvested
!     CPOOLG,ZPOOLG,PPOOLG=branch non-structural C,N,P mass after harvest
!     CPOLNG,ZPOLNG,PPOLNG=nonstructural C,N,P in bacteria after harvest
!     WTNDG,WTNDNG,WTNDPG=bacterial C,N,P mass after harvest
!     WTLS,LeafPetolBiomassC_brch=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
        CPOOLX=AZMAX1(NonstructElmnt_brch(ielmc,NB,NZ))
        ZPOOLX=AZMAX1(NonstructElmnt_brch(ielmn,NB,NZ))
        PPOOLX=AZMAX1(NonstructElmnt_brch(ielmp,NB,NZ))
        CPOLNX=AZMAX1(NoduleNonstructElmnt_brch(ielmc,NB,NZ))
        ZPOLNX=AZMAX1(NoduleNonstructElmnt_brch(ielmn,NB,NZ))
        PPOLNX=AZMAX1(NoduleNonstructElmnt_brch(ielmp,NB,NZ))
        IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
          IF(WGLFGY+WGSHGY.GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
            CPOOLG=CPOOLX*FHVSE(ielmc)
            ZPOOLG=ZPOOLX*FHVSE(ielmc)
            PPOOLG=PPOOLX*FHVSE(ielmc)
            CPOLNG=CPOLNX*FHVSE(ielmc)
            ZPOLNG=ZPOLNX*FHVSE(ielmc)
            PPOLNG=PPOLNX*FHVSE(ielmc)
            WTNDG=CanopyNoduleChemElmnt_brch(ielmc,NB,NZ)*FHVSE(ielmc)
            WTNDNG=CanopyNoduleChemElmnt_brch(ielmn,NB,NZ)*FHVSE(ielmc)
            WTNDPG=CanopyNoduleChemElmnt_brch(ielmp,NB,NZ)*FHVSE(ielmc)
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
          IF(CanopyLeafShethC_pft(NZ).GT.ZEROL(NZ))THEN
            WTLSBX=AZMAX1(LeafPetolBiomassC_brch(NB,NZ))
            IF(NonstructElmnt_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
              WHVSCX=AZMAX1(WHVSCP)*WTLSBX/CanopyLeafShethC_pft(NZ)
              CPOOLG=AZMAX1(CPOOLX-WHVSCX)
              ZPOOLG=AZMAX1(ZPOOLX-WHVSCX*ZPOOLX/NonstructElmnt_brch(ielmc,NB,NZ))
              PPOOLG=AZMAX1(PPOOLX-WHVSCX*PPOOLX/NonstructElmnt_brch(ielmc,NB,NZ))
            ELSE
              CPOOLG=0._r8
              ZPOOLG=0._r8
              PPOOLG=0._r8
            ENDIF
            IF(NoduleNonstructElmnt_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
              WHVSNX=AZMAX1(WHVSNP)*WTLSBX/CanopyLeafShethC_pft(NZ)
              CPOLNG=AZMAX1(CPOLNX-WHVSNX)
              ZPOLNG=AZMAX1(ZPOLNX-WHVSNX*ZPOLNX/NoduleNonstructElmnt_brch(ielmc,NB,NZ))
              PPOLNG=AZMAX1(PPOLNX-WHVSNX*PPOLNX/NoduleNonstructElmnt_brch(ielmc,NB,NZ))
              WTNDG=CanopyNoduleChemElmnt_brch(ielmc,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDNG=CanopyNoduleChemElmnt_brch(ielmn,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDPG=CanopyNoduleChemElmnt_brch(ielmp,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
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
        NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+CanopyNoduleChemElmnt_brch(ielmc,NB,NZ)-WTNDG
        NonstructElmntRemoval(ielmn)=NonstructElmntRemoval(ielmn)+CanopyNoduleChemElmnt_brch(ielmn,NB,NZ)-WTNDNG
        NonstructElmntRemoval(ielmp)=NonstructElmntRemoval(ielmp)+CanopyNoduleChemElmnt_brch(ielmp,NB,NZ)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
        NonstructElmnt_brch(ielmc,NB,NZ)=CPOOLG
        NonstructElmnt_brch(ielmn,NB,NZ)=ZPOOLG
        NonstructElmnt_brch(ielmp,NB,NZ)=PPOOLG
        NoduleNonstructElmnt_brch(ielmc,NB,NZ)=CPOLNG
        NoduleNonstructElmnt_brch(ielmn,NB,NZ)=ZPOLNG
        NoduleNonstructElmnt_brch(ielmp,NB,NZ)=PPOLNG
        CanopyNoduleChemElmnt_brch(ielmc,NB,NZ)=WTNDG
        CanopyNoduleChemElmnt_brch(ielmn,NB,NZ)=WTNDNG
        CanopyNoduleChemElmnt_brch(ielmp,NB,NZ)=WTNDPG
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CPOOLG=branch non-structural C mass after harvest
!     NonstructElmntRemoval(ielmc),NonstructElmntRemoval(ielmn),NonstructElmntRemoval(ielmp)=nonstructural C,N,P removed
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!
        IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo.AND.CPOOLX.GT.ZEROP(NZ))THEN
          FHVST4=CPOOLG/CPOOLX
          D9810: DO K=1,MaxNodesPerBranch1
            NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+(1._r8-FHVST4)*CPOOL3(K,NB,NZ)
            NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+(1._r8-FHVST4)*CPOOL4(K,NB,NZ)
            NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+(1._r8-FHVST4)*CMassCO2BundleSheath_node(K,NB,NZ)
            NonstructElmntRemoval(ielmc)=NonstructElmntRemoval(ielmc)+(1._r8-FHVST4)*CMassHCO3BundleSheath_node(K,NB,NZ)
            CPOOL3(K,NB,NZ)=FHVST4*CPOOL3(K,NB,NZ)
            CPOOL4(K,NB,NZ)=FHVST4*CPOOL4(K,NB,NZ)
            CMassCO2BundleSheath_node(K,NB,NZ)=FHVST4*CMassCO2BundleSheath_node(K,NB,NZ)
            CMassHCO3BundleSheath_node(K,NB,NZ)=FHVST4*CMassHCO3BundleSheath_node(K,NB,NZ)
          ENDDO D9810
        ENDIF
!
!     CUT STALKS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTSTKX=internode length removed
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     FHGT=fraction of canopy layer height not harvested
!     FHVSE(ielmc)=fraction of canopy layer mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WTSTK=stalk C mass
!
!
        IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
          IF(HTSTKX.GT.ZERO)THEN
            IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
              FHGT=AZMAX1(AMIN1(1.0_r8,HVST(NZ)/HTSTKX))
            ELSE
              FHGT=0._r8
            ENDIF
            IF(isclose(THIN_pft(NZ),0._r8))THEN
              FHVSE(ielmc)=AZMAX1(1._r8-(1._r8-FHGT)*EHVST(1,iplthvst_woody,NZ))
              FHVSH=FHVSE(ielmc)
            ELSE
              FHVSE(ielmc)=AZMAX1(1._r8-THIN_pft(NZ))
              IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
                FHVSH=1.0_r8-(1._r8-FHGT)*EHVST(1,iplthvst_woody,NZ)*THIN_pft(NZ)
              ELSE
                FHVSH=FHVSE(ielmc)
              ENDIF
            ENDIF
          ELSE
            FHVSE(ielmc)=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ELSE
          IF(StalkChemElmnts_pft(ielmc,NZ).GT.ZEROL(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/StalkChemElmnts_pft(ielmc,NZ)))
            FHVSH=FHVSE(ielmc)
          ELSE
            FHVSE(ielmc)=1.0_r8
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
        DO NE=1,NumPlantChemElmnts
          WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE)+(1._r8-FHVSH)*StalkChemElmnts_brch(NE,NB,NZ)
          WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE)+(FHVSH-FHVSE(ielmc))*StalkChemElmnts_brch(NE,NB,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
          StalkChemElmnts_brch(NE,NB,NZ)=FHVSE(ielmc)*StalkChemElmnts_brch(NE,NB,NZ)
          BranchStalkChemElmnts_pft(NE,NB,NZ)=FHVSE(ielmc)*BranchStalkChemElmnts_pft(NE,NB,NZ)
        ENDDO

        StalkBiomassC_brch(NB,NZ)=FHVSE(ielmc)*StalkBiomassC_brch(NB,NZ)
!
!     CUT STALK NODES
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     InternodeHeightDying_brch,InternodeHeightLive_brch=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTSTK=stalk C mass
!     InternodeChemElmnt_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
        D9820: DO K=MaxNodesPerBranch1,0,-1
          IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
            IF(InternodeHeightDying_brch(K,NB,NZ).GT.ZERO)THEN
              IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
                FHGTK=AZMAX1(AMIN1(1.0_r8,(InternodeHeightLive_brch(K,NB,NZ)-HVST(NZ))/&
                  InternodeHeightDying_brch(K,NB,NZ)))
              ELSE
                FHGTK=0._r8
              ENDIF
              IF(isclose(THIN_pft(NZ),0._r8))THEN
                FHVSETS=AZMAX1(1._r8-FHGTK*EHVST(1,iplthvst_woody,NZ))
              ELSE
                FHVSETS=AZMAX1(1._r8-THIN_pft(NZ))
              ENDIF
            ELSE
              FHVSETS=1.0_r8
            ENDIF
          ELSE
            IF(StalkChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
              FHVSETS=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/StalkChemElmnts_pft(ielmc,NZ)))
            ELSE
              FHVSETS=1.0_r8
            ENDIF
          ENDIF
          DO NE=1,NumPlantChemElmnts
            InternodeChemElmnt_brch(NE,K,NB,NZ)=FHVSETS*InternodeChemElmnt_brch(NE,K,NB,NZ)
          ENDDO
          IF(iHarvstType_pft(NZ).LE.iharvtyp_allabv.AND.isclose(THIN_pft(NZ),0._r8))THEN
            InternodeHeightDying_brch(K,NB,NZ)=FHVSETS*InternodeHeightDying_brch(K,NB,NZ)
            InternodeHeightLive_brch(K,NB,NZ)=AMIN1(InternodeHeightLive_brch(K,NB,NZ),HVST(NZ))
          ENDIF

        ENDDO D9820
!
!     CUT STALK RESERVES
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSTKB=C mass remaining in harvested stalk
!     WTRSV=stalk reserve C mass
!     WHVRVH=remaining stalk reserve C mass
!     FHVSE(ielmc)=fraction of reserve mass not harvested
!
        IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
          IF(StalkChemElmnts_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=FHVSE(ielmc)
            FHVSH=FHVSH
          ELSE
            FHVSE(ielmc)=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(ReserveChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,1._r8-WHVRVH/ReserveChemElmnts_pft(ielmc,NZ)))
            FHVSH=FHVSE(ielmc)
          ELSE
            FHVSE(ielmc)=0._r8
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
        DO NE=1,NumPlantChemElmnts
          WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE)+(1._r8-FHVSH)*ReserveElmnts_brch(NE,NB,NZ)
          WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(ielmc)+(FHVSH-FHVSE(ielmc))*ReserveElmnts_brch(NE,NB,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
          ReserveElmnts_brch(NE,NB,NZ)=FHVSE(ielmc)*ReserveElmnts_brch(NE,NB,NZ)
        ENDDO
!
!     CUT REPRODUCTIVE ORGANS
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FHVSETG,FHVSETH,FHVSETE=fraction of grain,husk,ear mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
        IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
          IF(HVST(NZ).LT.HTSTKX.OR.iHarvstType_pft(NZ).EQ.iharvtyp_grain.OR.iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
            IF(isclose(THIN_pft(NZ),0._r8))THEN
              FHVSETG=1.0_r8-EHVST(1,iplthvst_finenonleaf,NZ)
              FHVSHG=FHVSETG
            ELSE
              FHVSETG=1.0_r8-THIN_pft(NZ)
              FHVSHG=1.0_r8-EHVST(1,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
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
          IF(HuskChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETH=AZMAX1(AMIN1(1.0_r8,1._r8-WHVHSH/HuskChemElmnts_pft(ielmc,NZ)))
            FHVSHH=FHVSETH
          ELSE
            FHVSETH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(EarChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETE=AZMAX1(AMIN1(1.0_r8,1._r8-WHVEAH/EarChemElmnts_pft(ielmc,NZ)))
            FHVSHE=FHVSETE
          ELSE
            FHVSETE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(GrainChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETG=AZMAX1(AMIN1(1.0_r8,1._r8-WHVGRH/GrainChemElmnts_pft(ielmc,NZ)))
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
        DO NE=1,NumPlantChemElmnts
          FineNonleafElmntRemoval(NE)=FineNonleafElmntRemoval(NE)+(1._r8-FHVSHH)*HuskChemElmnts_brch(NE,NB,NZ)+(1._r8-FHVSHE) &
            *EarChemElmnts_brch(NE,NB,NZ)+(1._r8-FHVSHG)*GrainChemElmnts_brch(NE,NB,NZ)
          PetioleElmntHarv2Litr(NE)=PetioleElmntHarv2Litr(NE)+(FHVSHH-FHVSETH)*HuskChemElmnts_brch(NE,NB,NZ)+(FHVSHE-FHVSETE) &
            *EarChemElmnts_brch(NE,NB,NZ)+(FHVSHG-FHVSETG)*GrainChemElmnts_brch(NE,NB,NZ)
          WTHTGE(NE)=WTHTGE(NE)+(1._r8-FHVSETG)*GrainChemElmnts_brch(NE,NB,NZ)

!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
          HuskChemElmnts_brch(NE,NB,NZ)=FHVSETH*HuskChemElmnts_brch(NE,NB,NZ)
          EarChemElmnts_brch(NE,NB,NZ)=FHVSETE*EarChemElmnts_brch(NE,NB,NZ)
          GrainChemElmnts_brch(NE,NB,NZ)=FHVSETG*GrainChemElmnts_brch(NE,NB,NZ)

        ENDDO
        PotentialSeedSites_brch(NB,NZ)=FHVSETG*PotentialSeedSites_brch(NB,NZ)
        SeedNumberSet_brch(NB,NZ)=FHVSETG*SeedNumberSet_brch(NB,NZ)
        GrainSeedBiomCMean_brch(NB,NZ)=FHVSETG*GrainSeedBiomCMean_brch(NB,NZ)
!
!     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
!
!     ShootNonstructC_brch=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
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
        ShootNonstructC_brch(NB,NZ)=0._r8
        D1325: DO K=1,MaxNodesPerBranch1
          ShootNonstructC_brch(NB,NZ)=ShootNonstructC_brch(NB,NZ) &
            +CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
            +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ)
        ENDDO D1325

        LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafChemElmnts_brch(ielmc,NB,NZ) &
          +PetoleChemElmnt_brch(ielmc,NB,NZ))

        ShootChemElmnt_brch(ielmc,NB,NZ)=ShootChemElmnt_brch(ielmc,NB,NZ)+ShootNonstructC_brch(NB,NZ)
        DO NE=1,NumPlantChemElmnts
          ShootChemElmnt_brch(NE,NB,NZ)=AZMAX1(ShootChemElmnt_brch(ielmc,NB,NZ)+LeafChemElmnts_brch(NE,NB,NZ) &
            +PetoleChemElmnt_brch(NE,NB,NZ)+StalkChemElmnts_brch(NE,NB,NZ)+ReserveElmnts_brch(NE,NB,NZ) &
            +HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ)+GrainChemElmnts_brch(NE,NB,NZ) &
            +NonstructElmnt_brch(NE,NB,NZ))
        ENDDO
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
!     iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     ZC=canopy height
!     iPlantPhenologyType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     GROUP=node number required for floral initiation
!     NodeNumberToInitFloral_brch=node number at floral initiation
!     NodeNumberAtAnthesis_brch=node number at flowering
!     VSTGX=leaf number on date of floral initiation
!     TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
!     TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
!     HourFailGrainFill_brch=number of hours with no grain fill
!     doInitLeafOut_brch=flag for initializing leafout
!
        IF((iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType_pft(NZ)))) &
          .AND.(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
          .AND.CanopyHeight_pft(NZ).GT.HVST(NZ))THEN
          IF((iPlantPhenologyType_pft(NZ).NE.0.AND.Hours4LeafOff_brch(NB,NZ) &
            .LE.FracHour4LeafoffRemob(iPlantPhenologyType_pft(NZ))*HourThreshold4LeafOff_brch(NB,NZ)) &
            .OR.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_evgreen.AND.&
            iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0))THEN
            MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
            NodeNumberToInitFloral_brch(NB,NZ)=ShootNodeNumber_brch(NB,NZ)
            NodeNumberAtAnthesis_brch(NB,NZ)=0._r8
            LeafNumberAtFloralInit_brch(NB,NZ)=0._r8
            TotalNodeNumNormByMatgrp_brch(NB,NZ)=0._r8
            TotReproNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
            HourFailGrainFill_brch(NB,NZ)=0._r8
            iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)=I
            D3005: DO M=2,NumGrowthStages
              iPlantCalendar_brch(M,NB,NZ)=0
            ENDDO D3005
            doInitLeafOut_brch(NB,NZ)=0
            IF(NB.EQ.NumOfMainBranch_pft(NZ))THEN
              D3010: DO NBX=1,NumOfBranches_pft(NZ)
                IF(NBX.NE.NumOfMainBranch_pft(NZ))THEN
                  MatureGroup_brch(NBX,NZ)=MatureGroup_pft(NZ)
                  NodeNumberToInitFloral_brch(NBX,NZ)=ShootNodeNumber_brch(NBX,NZ)
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
!     CanopyBranchStemApft_lyr=total PFT stalk surface area
!
        IF(jHarvst_pft(NZ).NE.jharvtyp_noaction)then
          iPlantBranchState_brch(NB,NZ)=iDead
        endif
        IF(PlantPopulation_pft(NZ).LE.0.0)then
          iPlantBranchState_brch(NB,NZ)=iDead
        endif
      ENDDO D9835
      CanopyLeafShethC_pft(NZ)=0._r8
      StalkChemElmnts_pft(ielmc,NZ)=0._r8
      CanopyStalkC_pft(NZ)=0._r8
      CanopyStemA_pft(NZ)=0._r8
      D9840: DO NB=1,NumOfBranches_pft(NZ)
        CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)
        StalkChemElmnts_pft(ielmc,NZ)=StalkChemElmnts_pft(ielmc,NZ)+StalkChemElmnts_brch(ielmc,NB,NZ)
        CanopyStalkC_pft(NZ)=CanopyStalkC_pft(NZ)+StalkBiomassC_brch(NB,NZ)
        D9830: DO L=1,NumOfCanopyLayers1
          CanopyStemA_pft(NZ)=CanopyStemA_pft(NZ)+CanopyBranchStemApft_lyr(L,NB,NZ)
        ENDDO D9830
      ENDDO D9840
!
!     ROOT LITTERFALL FROM HARVESTING OR FIRE
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     THETW=soil water concentration
!     CORGC=SOC concentration
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     EFIRE=combustion  of N,P relative to C
!     FHVSE(ielmc),FHVSE(ielmn),FHVSE(ielmp)=fraction of root layer C,N,P not removed by disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CO2ByFire_pft,CH4ByFire_pft,VOXYF,NH3byFire_pft,N2ObyFire_pft,PO4byFire_pft=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_col=total net biome productivity
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
        XHVST(ielmc)=1.0_r8-THIN_pft(NZ)
        D3985: DO N=1,MY(NZ)
          D3980: DO L=NU,MaxNumRootLays
            IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
              XHVST(ielmc)=1.0_r8-THIN_pft(NZ)
              XHVST(ielmn)=XHVST(ielmc)
              XHVST(ielmp)=XHVST(ielmc)
              FFIRE(1:NumPlantChemElmnts)=0._r8
            ELSE
              IF(THETW(L).GT.FVLWB.OR.CORGC(L).LE.FORGC.OR.ITILL.NE.22)THEN
                XHVST(ielmc)=1.0_r8
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE(1:NumPlantChemElmnts)=0._r8
              ELSE
                XHVST(ielmc)=1.0_r8-DCORP*EHVST(1,iplthvst_woody,NZ) &
                  *AMIN1(1.0_r8,(CORGC(L)-FORGC)/(orgcden-FORGC))
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE(ielmc)=EHVST(2,iplthvst_woody,NZ)
                FFIRE(ielmn)=FFIRE(ielmc)*EFIRE(1,iHarvstType_pft(NZ))
                FFIRE(ielmp)=FFIRE(ielmc)*EFIRE(2,iHarvstType_pft(NZ))
              ENDIF
            ENDIF
            XHVST1=1._r8-XHVST
            D3385: DO M=1,jsken
              DO NE=1,NumPlantChemElmnts
                FHVSE(NE)=XHVST1(NE)*CFOPE(NE,inonstruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
                LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+(1._r8-FFIRE(NE))*FHVSE(NE)
              ENDDO
              CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
              CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
              VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667
              NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
              N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0_r8
              PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
              CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
              Eco_NBP_col=Eco_NBP_col-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
              DO NR=1,NumRootAxes_pft(NZ)
                DO NE=1,NumPlantChemElmnts
                  FHVSE(NE)=XHVST1(NE)*CFOPE(NE,icwood,M,NZ)*(Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ) &
                    +Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+&
                    (1._r8-FFIRE(NE))*FHVSE(NE)
                ENDDO
                CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667
                NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
                N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0
                PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
                CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                Eco_NBP_col=Eco_NBP_col-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)

                DO NE=1,NumPlantChemElmnts
                  FHVSE(NE)=XHVST1(NE)*CFOPE(NE,iroot,M,NZ)*(Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ) &
                    +Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ) &
                    +(1._r8-FFIRE(NE))*FHVSE(NE)
                ENDDO
                CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667_r8
                NH3byFire_pft(NZ)=NH3byFire_pft(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
                N2ObyFire_pft(NZ)=N2ObyFire_pft(NZ)-0.0_r8
                PO4byFire_pft(NZ)=PO4byFire_pft(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
                CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                Eco_NBP_col=Eco_NBP_col-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
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
              RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-XHVST1(ielmc) &
                *(trcg_rootml_vr(idg_CO2,N,L,NZ)+trcs_rootml_vr(idg_CO2,N,L,NZ))
              trcg_rootml_vr(NTG,N,L,NZ)=XHVST(ielmc)*trcg_rootml_vr(NTG,N,L,NZ)
              trcs_rootml_vr(NTG,N,L,NZ)=XHVST(ielmc)*trcs_rootml_vr(NTG,N,L,NZ)
            ENDDO

!
!     REDUCE ROOT STATE VARIABLES DURING HARVESTING
!
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     PrimRootLen,SecndRootLen=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootStructBiomC_vr, PopuPlantRootC_vr=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,SecndRootXNum_pvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_vr,RootVolume_vr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_vr=root surface area per plant
!     RootRespPotential_vr,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D3960: DO NR=1,NumRootAxes_pft(NZ)
              DO NE=1,NumPlantChemElmnts
                Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ)=Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ)*XHVST(NE)
                Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)=Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)*XHVST(NE)
                Root1stChemElmnt(NE,N,NR,NZ)=Root1stChemElmnt(NE,N,NR,NZ)*XHVST(NE)
              ENDDO
              PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)*XHVST(ielmc)
              SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)*XHVST(ielmc)
              SecndRootXNum_rpvr(N,L,NR,NZ)=SecndRootXNum_rpvr(N,L,NR,NZ)*XHVST(ielmc)
            ENDDO D3960
            DO NE=1,NumPlantChemElmnts
               RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)*XHVST(NE)
            ENDDO
            RootStructBiomC_vr(N,L,NZ)=RootStructBiomC_vr(N,L,NZ)*XHVST(ielmc)
             PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)*XHVST(ielmc)
            RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)*XHVST(ielmc)
            PrimRootXNumL_pvr(N,L,NZ)=PrimRootXNumL_pvr(N,L,NZ)*XHVST(ielmc)
            SecndRootXNum_pvr(N,L,NZ)=SecndRootXNum_pvr(N,L,NZ)*XHVST(ielmc)
            RootLenPerPlant_pvr(N,L,NZ)=RootLenPerPlant_pvr(N,L,NZ)*XHVST(ielmc)
            RootLenDensPerPlant_pvr(N,L,NZ)=RootLenDensPerPlant_pvr(N,L,NZ)*XHVST(ielmc)
            RootVolume_vr(N,L,NZ)=RootVolume_vr(N,L,NZ)*XHVST(ielmc)
            RootVH2O_vr(N,L,NZ)=RootVH2O_vr(N,L,NZ)*XHVST(ielmc)
            RootAreaPerPlant_vr(N,L,NZ)=RootAreaPerPlant_vr(N,L,NZ)*XHVST(ielmc)
            RootRespPotential_vr(N,L,NZ)=RootRespPotential_vr(N,L,NZ)*XHVST(ielmc)
            RCO2N(N,L,NZ)=RCO2N(N,L,NZ)*XHVST(ielmc)
            RCO2A(N,L,NZ)=RCO2A(N,L,NZ)*XHVST(ielmc)
!
!     NODULE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(is_plant_N2fix(iPlantNfixType(NZ)).AND.N.EQ.ipltroot)THEN
              DO NE=1,NumPlantChemElmnts
                D3395: DO M=1,jsken
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+ &
                    XHVST1(NE)*(CFOPE(NE,iroot,M,NZ)*RootNodueChemElmnt_pvr(NE,L,NZ) &
                    +CFOPE(NE,inonstruct,M,NZ)*RootNoduleNonstructElmnt_vr(NE,L,NZ))
                ENDDO D3395
                RootNodueChemElmnt_pvr(NE,L,NZ)=RootNodueChemElmnt_pvr(NE,L,NZ)*XHVST(NE)
                RootNoduleNonstructElmnt_vr(NE,L,NZ)=RootNoduleNonstructElmnt_vr(NE,L,NZ)*XHVST(NE)
              ENDDO
            ENDIF
          ENDDO D3980
        ENDDO D3985
!
!     STORAGE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        IF(iPlantPhenologyPattern_pft(NZ).NE.iplt_annual)THEN
          DO NE=1,NumPlantChemElmnts
            D3400: DO M=1,jsken
              LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
                +(XHVST1(NE)*CFOPE(NE,inonstruct,M,NZ)*NonstructalElmnts_pft(NE,NZ))*FWOODE(NE,k_woody_litr)

              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
                +(XHVST1(NE)*CFOPE(NE,inonstruct,M,NZ)*NonstructalElmnts_pft(NE,NZ))*FWOODE(NE,k_fine_litr)
            ENDDO D3400
            NonstructalElmnts_pft(NE,NZ)=NonstructalElmnts_pft(NE,NZ)*XHVST(NE)
          ENDDO
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

end module PlantDisturbsMod
