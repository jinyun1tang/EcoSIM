module PlantDisturbByGrazingMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use ElmIDMod
  use EcosimConst
  use PlantAPIData
  use GrosubPars
implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section
  public :: AbvgBiomRemovalByGrazing
  public :: RemoveStandDeadByGrazing
  public :: ApplyBiomRemovalByGrazing
  public :: CutBranchNonstalByGrazing
  public :: GrazingPlant

!     GY,GZ=partitioning of grazed material to removal,respiration
  real(r8),PARAMETER :: GY=1._r8
  real(r8),parameter :: GZ=1._r8-GY
contains
!------------------------------------------------------------------------------------------

  subroutine AbvgBiomRemovalByGrazing(I,J,NZ,TotalElmnt2Litr,TotalElmntRemoval)
  
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalElmntRemoval(NumPlantChemElms)  
  real(r8), intent(in) :: TotalElmnt2Litr(NumPlantChemElms)  
  integer :: NE
  real(r8) :: dRespC
  associate(                                                      &
    Eco_AutoR_CumYr_col     => plt_bgcr%Eco_AutoR_CumYr_col,      &
    CanopyRespC_CumYr_pft   => plt_bgcr%CanopyRespC_CumYr_pft,    &
    ECO_ER_col              => plt_bgcr%ECO_ER_col,               &
    GrossResp_pft           => plt_bgcr%GrossResp_pft,            &
    EcoHavstElmnt_CumYr_col => plt_distb%EcoHavstElmnt_CumYr_col, &
    EcoHavstElmnt_CumYr_pft => plt_distb%EcoHavstElmnt_CumYr_pft  &
  )
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     GY=growth yield of grazers
!     TotalElmntRemoval(ielmc),TotalElmntRemoval(ielmn),TotalElmntRemoval(ielmp)=total C,N,P removed
!     TotalElmnt2Litr(ielmc),TotalElmnt2Litr(ielmn),TotalElmnt2Litr(ielmp)=total C,N,P to litter
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_CumYr_col=total autotrophic respiration

  EcoHavstElmnt_CumYr_pft(ielmc,NZ)=EcoHavstElmnt_CumYr_pft(ielmc,NZ)+GY*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  EcoHavstElmnt_CumYr_col(ielmc)=EcoHavstElmnt_CumYr_col(ielmc)+GY*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))

  DO NE=2,NumPlantChemElms
    EcoHavstElmnt_CumYr_pft(NE,NZ)=EcoHavstElmnt_CumYr_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
    EcoHavstElmnt_CumYr_col(NE)=EcoHavstElmnt_CumYr_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
  ENDDO
  !GZ : fraction of removal into respiration by herbivory/grazing
  dRespC=GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  GrossResp_pft(NZ)=GrossResp_pft(NZ)-dRespC
  CanopyRespC_CumYr_pft(NZ)=CanopyRespC_CumYr_pft(NZ)-dRespC
!     Eco_NBP_CumYr_col=Eco_NBP_CumYr_col+GY*(TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc))
!     CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-dRespC
  ECO_ER_col=ECO_ER_col-dRespC
  Eco_AutoR_CumYr_col=Eco_AutoR_CumYr_col-dRespC
  end associate
  end subroutine AbvgBiomRemovalByGrazing

!------------------------------------------------------------------------------------------

  subroutine RemoveStandDeadByGrazing(I,J,NZ,FHVSE,FHVSH)
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: NZ
  real(r8), intent(out) :: FHVSE
  real(r8), intent(out) :: FHVSH
  real(r8) :: WHVSTD

  associate(   &
    FracCanopyHeightCut_pft => plt_distb%FracCanopyHeightCut_pft, &
    NU                      => plt_site%NU,                       &
    AREA3                   => plt_site%AREA3,                    &
    StandDeadStrutElms_pft  => plt_biom%StandDeadStrutElms_pft,   &
    FracBiomHarvsted        => plt_distb%FracBiomHarvsted,        &
    THIN_pft                => plt_distb%THIN_pft,                &
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            &
    )
  IF(StandDeadStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    WHVSTD=FracCanopyHeightCut_pft(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8*AREA3(NU)*FracBiomHarvsted(1,4,NZ)
    FHVSE=AZMAX1(1._r8-WHVSTD/StandDeadStrutElms_pft(ielmc,NZ))
    FHVSH=FHVSE
  ELSE
    FHVSE=1.0_r8
    FHVSH=1.0_r8
  ENDIF
  end associate  
  end subroutine RemoveStandDeadByGrazing
!------------------------------------------------------------------------------------------

  subroutine ApplyBiomRemovalByGrazing(I,J,NZ,EHVST21,EHVST22,EHVST23,EHVST24,&
    NonstructElmntRemoval,LeafElmntRemoval,FineNonleafElmntRemoval,WoodyElmntRemoval,StandeadElmntRemoval,&
    NonstructElmnt2Litr,LeafElmnt2Litr,FineNonleafElmnt2Litr,WoodyElmnt2Litr,StandeadElmnt2Litr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: EHVST21
  real(r8), intent(in) :: EHVST22
  real(r8), intent(in) :: EHVST23
  real(r8), intent(in) :: EHVST24
  real(r8), intent(in) :: NonstructElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntRemoval(NumPlantChemElms)
  real(r8), intent(out) :: NonstructElmnt2Litr(NumPlantChemElms)  
  real(r8), intent(out) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmnt2Litr(NumPlantChemElms)
  real(r8) :: EHVST21h,EHVST22h,EHVST23h,EHVST24h
  integer :: NE
  associate(                                        &
    NU               => plt_site%NU,                &
    AREA3            => plt_site%AREA3,             &
    FERT             => plt_distb%FERT,             &
    FracBiomHarvsted => plt_distb%FracBiomHarvsted, &
    IYTYP            => plt_distb%IYTYP             &
  )
  EHVST21h=1._r8-FracBiomHarvsted(2,iplthvst_leaf,NZ)*0.5_r8
  EHVST22h=1._r8-FracBiomHarvsted(2,iplthvst_finenonleaf,NZ)*0.5_r8
  EHVST23h=1._r8-FracBiomHarvsted(2,iplthvst_woody,NZ)*0.5_r8
  EHVST24h=1._r8-FracBiomHarvsted(2,iplthvst_stdead,NZ)*0.5_r8

  NonstructElmnt2Litr(ielmc)   = NonstructElmntRemoval(ielmc)*EHVST21
  LeafElmnt2Litr(ielmc)        = LeafElmntRemoval(ielmc)*EHVST21
  FineNonleafElmnt2Litr(ielmc) = FineNonleafElmntRemoval(ielmc)*EHVST22
  WoodyElmnt2Litr(ielmc)       = WoodyElmntRemoval(ielmc)*EHVST23
  StandeadElmnt2Litr(ielmc)    = StandeadElmntRemoval(ielmc)*EHVST24

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
  end associate
  end subroutine ApplyBiomRemovalByGrazing


!--------------------------------------------------------------------------------

  subroutine GrazingPlant(I,J,NZ,HvstedLeafC,HvstedShethC,HvstedEarC,HvstedGrainC,&
    WHVSCP,HvstedStalkC,HvstedRsrvC,WHVSHH,WHVSNP,LeafC_lbrch)
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
  REAL(R8), intent(out) :: LeafC_lbrch(NumOfCanopyLayers1,JP1,JP1)  
  real(r8) :: totShootC,TotPhytomassRemoval
  real(r8) :: WHVSLX,WHVSLY,WHVSCL,WHVSNL,CGrazedDeficit,WHVSSX
  real(r8) :: WHVSTY  
  real(r8) :: TotStalkC  
  real(r8) :: WHVSHY,WHVSCS,WHVSNS,WHVEAX,WHVEAY,WHVGRX,WHVGRY  
  real(r8) :: WHVSHX,WHVHSX,WHVHSY,WHVRVX,WHVRVY,WHVSKX,WHVSTX
  real(r8) :: CCPOLX  
  real(r8) :: CCPLNX  
  integer :: K,NB,L

  associate(                                                             &
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft,          &
    LeafElmsByLayerNode_brch => plt_biom%LeafElmsByLayerNode_brch, &
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,              &
    FracCanopyHeightCut_pft     => plt_distb%FracCanopyHeightCut_pft,    &
    iHarvstType_pft             => plt_distb%iHarvstType_pft,            &
    LeafStrutElms_pft           => plt_biom%LeafStrutElms_pft,           &
    HuskStrutElms_pft           => plt_biom%HuskStrutElms_pft,           &
    ShootStrutElms_pft          => plt_biom%ShootStrutElms_pft,          &
    AREA3                       => plt_site%AREA3,                       &
    NU                          => plt_site%NU,                          &
    FracBiomHarvsted            => plt_distb%FracBiomHarvsted,           &
    THIN_pft                    => plt_distb%THIN_pft,                   &
    NoduleNonstructCconc_pft    => plt_biom%NoduleNonstructCconc_pft,    &
    GrainStrutElms_pft          => plt_biom%GrainStrutElms_pft,          &
    StalkStrutElms_pft          => plt_biom%StalkStrutElms_pft,          &
    StalkRsrvElms_pft           => plt_biom%StalkRsrvElms_pft,           &
    EarStrutElms_pft            => plt_biom%EarStrutElms_pft,            &
    PetoleStrutElms_pft         => plt_biom%PetoleStrutElms_pft,         &
    CanopyNonstElmConc_pft      => plt_biom%CanopyNonstElmConc_pft,      &
    fTCanopyGroth_pft           => plt_pheno%fTCanopyGroth_pft,          &
    AvgCanopyBiomC2Graze_pft    => plt_biom%AvgCanopyBiomC2Graze_pft     &
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
! 24. has unit of hour
  IF(AvgCanopyBiomC2Graze_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
    TotPhytomassRemoval=FracCanopyHeightCut_pft(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8 &
      *AREA3(NU)*ShootStrutElms_pft(ielmc,NZ)/AvgCanopyBiomC2Graze_pft(NZ)
  ELSE
    TotPhytomassRemoval=0._r8
  ENDIF
  !herbivory, means fauna/insect cutting?
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
!     CGrazedDeficit=grazing requirement unmet by leaf
! grazing from leaf, 
  WHVSLX=TotPhytomassRemoval*FracBiomHarvsted(1,iplthvst_leaf,NZ)
  WHVSLY=AMIN1(LeafStrutElms_pft(ielmc,NZ),WHVSLX)
  HvstedLeafC=WHVSLY*(1._r8-CCPOLX)
  
  WHVSCL=WHVSLY*CCPOLX
  WHVSNL=WHVSLY*CCPLNX
  CGrazedDeficit=AZMAX1(WHVSLX-WHVSLY)
  WHVSSX=TotPhytomassRemoval*FracBiomHarvsted(1,iplthvst_finenonleaf,NZ)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     CGrazedDeficit=grazing requirement unmet by non-foliar removal
!
  totShootC=PetoleStrutElms_pft(ielmc,NZ)+HuskStrutElms_pft(ielmc,NZ) &
    +EarStrutElms_pft(ielmc,NZ)+GrainStrutElms_pft(ielmc,NZ)
  
  IF(totShootC.GT.ZERO4Groth_pft(NZ))THEN
    WHVSHX = WHVSSX*PetoleStrutElms_pft(ielmc,NZ)/totShootC+CGrazedDeficit
    WHVSHY = AMIN1(PetoleStrutElms_pft(ielmc,NZ),WHVSHX)
    WHVSHH = WHVSHY*(1._r8-CCPOLX)
    WHVSCS = WHVSHY*CCPOLX
    WHVSNS = WHVSHY*CCPLNX

    CGrazedDeficit = AZMAX1(WHVSHX-WHVSHY)
    WHVHSX         = WHVSSX*HuskStrutElms_pft(ielmc,NZ)/totShootC+CGrazedDeficit
    WHVHSY         = AMIN1(HuskStrutElms_pft(ielmc,NZ),WHVHSX)
    HvstedShethC   = WHVHSY

    CGrazedDeficit = AZMAX1(WHVHSX-WHVHSY)
    WHVEAX         = WHVSSX*EarStrutElms_pft(ielmc,NZ)/totShootC+CGrazedDeficit
    WHVEAY         = AMIN1(EarStrutElms_pft(ielmc,NZ),WHVEAX)
    HvstedEarC     = WHVEAY

    CGrazedDeficit = AZMAX1(WHVEAX-WHVEAY)
    WHVGRX         = WHVSSX*GrainStrutElms_pft(ielmc,NZ)/totShootC+CGrazedDeficit
    WHVGRY         = AMIN1(GrainStrutElms_pft(ielmc,NZ),WHVGRX)
    HvstedGrainC   = WHVGRY
    CGrazedDeficit = AZMAX1(WHVGRX-WHVGRY)
  ELSE
    WHVSHH         = 0._r8
    WHVSCS         = 0._r8
    WHVSNS         = 0._r8
    HvstedShethC   = 0._r8
    HvstedEarC     = 0._r8
    HvstedGrainC   = 0._r8
    CGrazedDeficit = CGrazedDeficit+WHVSSX
  ENDIF
  WHVSCP=WHVSCL+WHVSCS
  WHVSNP=WHVSNL+WHVSNS
  WHVSKX=TotPhytomassRemoval*FracBiomHarvsted(1,iplthvst_woody,NZ)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     CGrazedDeficit=grazing requirement unmet by stalk,reserve
!
  TotStalkC=StalkStrutElms_pft(ielmc,NZ)+StalkRsrvElms_pft(ielmc,NZ)
  IF(TotStalkC.GT.WHVSKX+CGrazedDeficit)THEN
    WHVSTX=WHVSKX*StalkStrutElms_pft(ielmc,NZ)/TotStalkC+CGrazedDeficit
    WHVSTY=AMIN1(StalkStrutElms_pft(ielmc,NZ),WHVSTX)
    HvstedStalkC=WHVSTY

    CGrazedDeficit=AZMAX1(WHVSTX-WHVSTY)
    WHVRVX=WHVSKX*StalkRsrvElms_pft(ielmc,NZ)/TotStalkC+CGrazedDeficit
    WHVRVY=AMIN1(StalkRsrvElms_pft(ielmc,NZ),WHVRVX)
    HvstedRsrvC=WHVRVY
    CGrazedDeficit=AZMAX1(WHVRVX-WHVRVY)
  ELSE
    HvstedStalkC=0._r8
    HvstedRsrvC=0._r8
    CGrazedDeficit=AZMAX1(WHVSKX)
!
!     ALLOCATE UNMET DEMAND FOR GRAZING TO LEAF,PETIOLE,HUSK
!     EAR,GRAIN
!
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WHVSH*,WHVHS,WHVEA,WHVGR,WHVSC=
!            petiole,husk,ear,grain,nonstructural C removed
!
    IF(CGrazedDeficit.GT.0.0_r8)THEN
      WHVSLY=AMIN1(LeafStrutElms_pft(ielmc,NZ)-HvstedLeafC-WHVSCL,CGrazedDeficit)
      HvstedLeafC=HvstedLeafC+WHVSLY*(1._r8-CCPOLX)

      WHVSCL=WHVSCL+WHVSLY*CCPOLX
      WHVSNL=WHVSNL+WHVSLY*CCPLNX

      CGrazedDeficit=AZMAX1(CGrazedDeficit-WHVSLY)

      IF(totShootC.GT.ZERO4Groth_pft(NZ))THEN
        WHVSHX = CGrazedDeficit*PetoleStrutElms_pft(ielmc,NZ)/totShootC
        WHVSHY = AMIN1(PetoleStrutElms_pft(ielmc,NZ),WHVSHX)
        WHVSHH = WHVSHH+WHVSHY*(1._r8-CCPOLX)
        WHVSCS = WHVSCS+WHVSHY*CCPOLX
        WHVSNS = WHVSNS+WHVSHY*CCPLNX

        CGrazedDeficit = AZMAX1(CGrazedDeficit-WHVSHY)
        WHVHSX         = CGrazedDeficit*HuskStrutElms_pft(ielmc,NZ)/totShootC
        WHVHSY         = AMIN1(HuskStrutElms_pft(ielmc,NZ),WHVHSX)
        HvstedShethC   = HvstedShethC+WHVHSY

        CGrazedDeficit = AZMAX1(CGrazedDeficit-WHVHSY)
        WHVEAX         = CGrazedDeficit*EarStrutElms_pft(ielmc,NZ)/totShootC
        WHVEAY         = AMIN1(EarStrutElms_pft(ielmc,NZ),WHVEAX)
        HvstedEarC     = HvstedEarC+WHVEAY

        CGrazedDeficit = AZMAX1(WHVEAX-WHVEAY)
        WHVGRX         = CGrazedDeficit*GrainStrutElms_pft(ielmc,NZ)/totShootC
        WHVGRY         = AMIN1(GrainStrutElms_pft(ielmc,NZ),WHVGRX)
        HvstedGrainC   = HvstedGrainC+WHVGRY
        
        CGrazedDeficit = AZMAX1(WHVGRX-WHVGRY)
      ENDIF
    ENDIF
  ENDIF

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
        LeafC_lbrch(L,NB,NZ)=LeafC_lbrch(L,NB,NZ)+LeafElmsByLayerNode_brch(ielmc,L,K,NB,NZ)
      enddo
    enddo
  ENDDO D9870
  end associate
  end subroutine GrazingPlant

!--------------------------------------------------------------------------------

  subroutine CutBranchNonstalByGrazing(I,J,NB,NZ,WHVSCP,CanopyNonstElmCopy,&
    WHVSNP,CanopyNodulNonstElmCopy,CanopyNonstElmAfhvst,CanopyNodulNonstElmAfhvst,&
    CanopyNodulStrutElmAfhvst)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: WHVSCP
  real(r8), intent(in) :: CanopyNonstElmCopy(NumPlantChemElms)
  real(r8), intent(in) :: WHVSNP
  real(r8), intent(in) :: CanopyNodulNonstElmCopy(NumPlantChemElms)
  REAL(R8),intent(out):: CanopyNonstElmAfhvst(NumPlantChemElms)
  REAL(R8),intent(out):: CanopyNodulNonstElmAfhvst(NumPlantChemElms)
  real(r8),intent(out) :: CanopyNodulStrutElmAfhvst(NumPlantChemElms)
  real(r8) :: WHVSNX
  real(r8) :: WHVSCX
  real(r8) :: WTLSBX  
  integer  :: NE
  associate(                                                         &
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,      &
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch, &
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,    &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch, &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft           &
  )
  IF(CanopyLeafShethC_pft(NZ).GT.ZERO4LeafVar_pft(NZ))THEN
    WTLSBX=AZMAX1(LeafPetolBiomassC_brch(NB,NZ))
    IF(CanopyNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      WHVSCX=AZMAX1(WHVSCP)*WTLSBX/CanopyLeafShethC_pft(NZ)
      CanopyNonstElmAfhvst(ielmc)=AZMAX1(CanopyNonstElmCopy(ielmc)-WHVSCX)
      DO NE=2,NumPlantChemElms
        CanopyNonstElmAfhvst(NE)=AZMAX1(CanopyNonstElmCopy(NE) &
          -WHVSCX*CanopyNonstElmCopy(NE)/CanopyNonstElms_brch(ielmc,NB,NZ))
      ENDDO
    ELSE
      CanopyNonstElmAfhvst(:)=0._r8
    ENDIF

    IF(CanopyNodulNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      WHVSNX=AZMAX1(WHVSNP)*WTLSBX/CanopyLeafShethC_pft(NZ)
      CanopyNodulNonstElmAfhvst(ielmc)=AZMAX1(CanopyNodulNonstElmCopy(ielmc)-WHVSNX)
      DO NE=2,NumPlantChemElms
        CanopyNodulNonstElmAfhvst(NE)=AZMAX1(CanopyNodulNonstElmCopy(NE) &
          -WHVSNX*CanopyNodulNonstElmCopy(NE)/CanopyNodulNonstElms_brch(ielmc,NB,NZ))
      ENDDO

      DO NE=1,NumPlantChemElms
        CanopyNodulStrutElmAfhvst(NE)=CanopyNodulStrutElms_brch(NE,NB,NZ)*(1._r8-WHVSNX/CanopyNodulNonstElmCopy(NE))
      ENDDO
    ELSE
      CanopyNodulNonstElmAfhvst(:)=0._r8
      CanopyNodulStrutElmAfhvst(:)=0._r8
    ENDIF
  ELSE
    CanopyNonstElmAfhvst(:)=0._r8
    CanopyNodulNonstElmAfhvst(:)=0._r8
    CanopyNodulStrutElmAfhvst(:)=0._r8
  ENDIF
  end associate
  end subroutine CutBranchNonstalByGrazing
end module PlantDisturbByGrazingMod