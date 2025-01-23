module PlantDisturbByFireMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ElmIDMod
  use PlantAPIData
  use GrosubPars
  use EcoSimConst
  use PlantMathFuncMod
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: StageRootRemovalByFire
  public :: RemoveRootByFire
  public :: AbvGrndLiterFallByFire
  public :: AbvgBiomRemovalByFire
  public :: ApplyBiomRemovalByFire
  public :: InitPlantFireMod

  real(r8) :: EFIRE(2,5:5)

contains
!--------------------------------------------------------------------------------

  subroutine InitPlantFireMod
  implicit none
  EFIRE=reshape(real((/0.917,0.167/),r8),shape(EFIRE))
  end subroutine InitPlantFireMod

!--------------------------------------------------------------------------------
  subroutine StageRootRemovalByFire(I,J,NZ,L,FFIRE,DCORP,FracLeftThin)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  integer, intent(in) :: L
  real(r8), intent(in):: DCORP
  real(r8),intent(out) :: FFIRE(NumPlantChemElms)
  real(r8), intent(out):: FracLeftThin

  associate(                                                &
    iSoilDisturbType_col => plt_distb%iSoilDisturbType_col, &
    CSoilOrgM_vr         => plt_soilchem%CSoilOrgM_vr,      &
    FracBiomHarvsted     => plt_distb%FracBiomHarvsted,     &
    iHarvstType_pft      => plt_distb%iHarvstType_pft,      &
    THETW_vr             => plt_soilchem%THETW_vr           &
  )
  IF(THETW_vr(L).GT.VolMaxSoilMoist4Fire .OR. CSoilOrgM_vr(ielmc,L).LE.FORGC .OR. iSoilDisturbType_col.NE.22)THEN
    FracLeftThin=1.0_r8
    FFIRE(1:NumPlantChemElms)=0._r8
  ELSE
    FracLeftThin=1.0_r8-DCORP*FracBiomHarvsted(1,iplthvst_woody,NZ) &
      *AMIN1(1.0_r8,(CSoilOrgM_vr(ielmc,L)-FORGC)/(orgcden-FORGC))
    FFIRE(ielmc) = FracBiomHarvsted(2,iplthvst_woody,NZ)
    FFIRE(ielmn) = FFIRE(ielmc)*EFIRE(1,iHarvstType_pft(NZ))
    FFIRE(ielmp) = FFIRE(ielmc)*EFIRE(2,iHarvstType_pft(NZ))
  ENDIF
  end associate
  end subroutine StageRootRemovalByFire


!--------------------------------------------------------------------------------

  subroutine RemoveRootByFire(I,J,NZ,FrcMassNotHarvst,FFIRE)
  !
  !upon combustion:
  !plant C is turned into CO2, or CH4
  !plant N is turned into NH3, or N2O (which is set to zero)
  !plant P is turned in PO4.
  implicit none

  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), INTENT(IN) :: FrcMassNotHarvst(NumPlantChemElms)
  real(r8), intent(in) :: FFIRE(NumPlantChemElms)
  associate(                                              &
    iHarvstType_pft     => plt_distb%iHarvstType_pft,     &
    CO2NetFix_pft       => plt_bgcr%CO2NetFix_pft,        &
    PO4byFire_CumYr_pft => plt_distb%PO4byFire_CumYr_pft, &
    N2ObyFire_CumYr_pft => plt_distb%N2ObyFire_CumYr_pft, &
    NH3byFire_CumYr_pft => plt_distb%NH3byFire_CumYr_pft, &
    O2ByFire_CumYr_pft  => plt_distb%O2ByFire_CumYr_pft,  &
    CH4ByFire_CumYr_pft => plt_distb%CH4ByFire_CumYr_pft, &
    CO2ByFire_CumYr_pft => plt_distb%CO2ByFire_CumYr_pft, &
    Eco_NBP_CumYr_col   => plt_bgcr%Eco_NBP_CumYr_col     &
  )

  CO2ByFire_CumYr_pft(NZ) = CO2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcMassNotHarvst(ielmc)
  CH4ByFire_CumYr_pft(NZ) = CH4ByFire_CumYr_pft(NZ)-FrcAsCH4byFire*FFIRE(ielmc)*FrcMassNotHarvst(ielmc)
  O2ByFire_CumYr_pft(NZ)  = O2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcMassNotHarvst(ielmc)*2.667_r8
  NH3byFire_CumYr_pft(NZ) = NH3byFire_CumYr_pft(NZ)-FFIRE(ielmn)*FrcMassNotHarvst(ielmn)
  N2ObyFire_CumYr_pft(NZ) = N2ObyFire_CumYr_pft(NZ)-0.0_r8
  PO4byFire_CumYr_pft(NZ) = PO4byFire_CumYr_pft(NZ)-FFIRE(ielmp)*FrcMassNotHarvst(ielmp)
  CO2NetFix_pft(NZ)       = CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*FFIRE(ielmc)*FrcMassNotHarvst(ielmc)
  Eco_NBP_CumYr_col       = Eco_NBP_CumYr_col-FrcAsCH4byFire*FFIRE(ielmc)*FrcMassNotHarvst(ielmc)
  end associate
  end subroutine RemoveRootByFire

!------------------------------------------------------------------------------------------
  subroutine AbvGrndLiterFallByFire(I,J,NZ,NonstructElmnt2Litr,StandeadElmntOffEcosystem, &
    FineNonleafElmOffEcosystem,LeafElmnt2Litr,LeafElmntOffEcosystem,NonstructElmntOffEcosystem,&
    WoodyElmntOffEcosystem,WoodyElmnt2Litr,StandeadElmnt2Litr,PetioleElmntHarv2Litr,&
    FineNonleafElmnt2Litr,LeafElmntHarv2Litr,StandeadElmntHarv2Litr,WoodyElmntHarv2Litr)
  implicit none 
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: NonstructElmnt2Litr(NumPlantChemElms)  
  real(r8), intent(in) :: StandeadElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: FineNonleafElmOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: LeafElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: WoodyElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: PetioleElmntHarv2Litr(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntHarv2Litr(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntHarv2Litr(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntHarv2Litr(NumPlantChemElms)
  integer :: M,NE

  associate(                                                             &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,        &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,        &
    StandDeadKCompElms_pft      => plt_biom%StandDeadKCompElms_pft,      &
    icwood                      => pltpar%icwood,                        &
    istalk                      => pltpar%istalk,                        &
    iroot                       => pltpar%iroot,                         &
    ilignin                     => pltpar%ilignin,                       &
    inonstruct                  => pltpar%inonstruct,                    &
    k_fine_litr                 => pltpar%k_fine_litr,                   &
    k_woody_litr                => pltpar%k_woody_litr,                  &
    ifoliar                     => pltpar%ifoliar,                       &
    inonfoliar                  => pltpar%inonfoliar,                    &
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr, &    
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft,  &    
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft       &    
  )

  D6485: DO M=1,jsken
    LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
      +ElmAllocmat4Litr(ielmc,inonstruct,M,NZ)*NonstructElmnt2Litr(ielmc) &
      +ElmAllocmat4Litr(ielmc,ifoliar,M,NZ)   *LeafElmnt2Litr(ielmc)+LeafElmntHarv2Litr(ielmc) &
      +ElmAllocmat4Litr(ielmc,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmc)+PetioleElmntHarv2Litr(ielmc))

    DO NE=2,NumPlantChemElms
      LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
        +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*NonstructElmntOffEcosystem(NE) &
        +ElmAllocmat4Litr(NE,ifoliar,M,NZ)   *LeafElmntOffEcosystem(NE) &
        +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*FineNonleafElmOffEcosystem(NE)

      LitrfalStrutElms_pvr(NE,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,ilignin,k_fine_litr,0,NZ) &
        +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*(NonstructElmnt2Litr(NE)-NonstructElmntOffEcosystem(NE)) &
        +ElmAllocmat4Litr(NE,ifoliar,M,NZ)   *(LeafElmnt2Litr(NE)+LeafElmntHarv2Litr(NE)-LeafElmntOffEcosystem(NE)) &
        +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(NE)+PetioleElmntHarv2Litr(NE) &
        -FineNonleafElmOffEcosystem(NE))
    ENDDO

    !all above ground is subject to fire
    IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
      LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)+&
        ElmAllocmat4Litr(ielmc,istalk,M,NZ)*(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc) &
        +StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))

      DO NE=2,NumPlantChemElms
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+&
          ElmAllocmat4Litr(NE,istalk,M,NZ)*(WoodyElmntOffEcosystem(NE)+StandeadElmntOffEcosystem(NE))

        LitrfalStrutElms_pvr(NE,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,ilignin,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(NE,istalk,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE) &
          -WoodyElmntOffEcosystem(NE)+StandeadElmnt2Litr(NE)+StandeadElmntHarv2Litr(NE)&
          -StandeadElmntOffEcosystem(NE))
      ENDDO
    !not grass, has woody component 
    ELSE
      StandDeadKCompElms_pft(ielmc,M,NZ)=StandDeadKCompElms_pft(ielmc,M,NZ)+ElmAllocmat4Litr(ielmc,icwood,M,NZ) &
        *(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc))
      DO NE=2,NumPlantChemElms  
        StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ)+ElmAllocmat4Litr(NE,icwood,M,NZ) &
          *WoodyElmntOffEcosystem(NE)
      ENDDO
        
      LitrfalStrutElms_pvr(ielmc,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_woody_litr,0,NZ) &
        *ElmAllocmat4Litr(ielmc,istalk,M,NZ)&
        *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FracRootStalkElmAlloc2Litr(ielmc,k_woody_litr)

      DO NE=2,NumPlantChemElms  
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
          +ElmAllocmat4Litr(NE,istalk,M,NZ)*StandeadElmntOffEcosystem(NE)*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

        LitrfalStrutElms_pvr(NE,ilignin,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,ilignin,k_woody_litr,0,NZ) &
          +ElmAllocmat4Litr(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)-WoodyElmntOffEcosystem(NE) &
          +StandeadElmnt2Litr(NE)+StandeadElmntHarv2Litr(NE)-StandeadElmntOffEcosystem(NE)) &
          *FracRootStalkElmAlloc2Litr(NE,k_woody_litr)
      ENDDO    

      LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
        +ElmAllocmat4Litr(ielmc,istalk,M,NZ) &
        *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FracRootStalkElmAlloc2Litr(ielmc,k_fine_litr)

      DO NE=2,NumPlantChemElms    
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(NE,istalk,M,NZ)*StandeadElmntOffEcosystem(NE)*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
        
        LitrfalStrutElms_pvr(NE,ilignin,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,ilignin,k_fine_litr,0,NZ) &
          +ElmAllocmat4Litr(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)-WoodyElmntOffEcosystem(NE) &
          +StandeadElmnt2Litr(NE)+StandeadElmntHarv2Litr(NE)-StandeadElmntOffEcosystem(NE)) &
          *FracRootStalkElmAlloc2Litr(NE,k_fine_litr)

      ENDDO  
    ENDIF
  ENDDO D6485
  end associate
  end subroutine AbvGrndLiterFallByFire

!------------------------------------------------------------------------------------------
  subroutine AbvgBiomRemovalByFire(I,J,NZ,TotalElmnt2Litr,TotalElmntRemoval)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: TotalElmntRemoval(NumPlantChemElms)

  associate(                                              &
    Eco_NBP_CumYr_col   => plt_bgcr%Eco_NBP_CumYr_col,    &
    CO2NetFix_pft       => plt_bgcr%CO2NetFix_pft,        &
    NH3byFire_CumYr_pft => plt_distb%NH3byFire_CumYr_pft, &
    PO4byFire_CumYr_pft => plt_distb%PO4byFire_CumYr_pft, &
    CH4ByFire_CumYr_pft => plt_distb%CH4ByFire_CumYr_pft, &
    O2ByFire_CumYr_pft  => plt_distb%O2ByFire_CumYr_pft,  &
    N2ObyFire_CumYr_pft => plt_distb%N2ObyFire_CumYr_pft, &
    CO2ByFire_CumYr_pft => plt_distb%CO2ByFire_CumYr_pft  &
  )
  CO2ByFire_CumYr_pft(NZ) = CO2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  CH4ByFire_CumYr_pft(NZ) = CH4ByFire_CumYr_pft(NZ)-FrcAsCH4byFire*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  O2ByFire_CumYr_pft(NZ)  = O2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))*2.667_r8
  NH3byFire_CumYr_pft(NZ) = NH3byFire_CumYr_pft(NZ)-TotalElmntRemoval(ielmn)+TotalElmnt2Litr(ielmn)
  N2ObyFire_CumYr_pft(NZ) = N2ObyFire_CumYr_pft(NZ)-0.0_r8
  PO4byFire_CumYr_pft(NZ) = PO4byFire_CumYr_pft(NZ)-TotalElmntRemoval(ielmp)+TotalElmnt2Litr(ielmp)
  CO2NetFix_pft(NZ)       = CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  Eco_NBP_CumYr_col       = Eco_NBP_CumYr_col-FrcAsCH4byFire*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  end associate
  end subroutine AbvgBiomRemovalByFire

!------------------------------------------------------------------------------------------
  subroutine ApplyBiomRemovalByFire(I,J,NZ,EHVST21,EHVST22, EHVST23, EHVST24,&
    StandeadElmntRemoval,NonstructElmntRemoval,LeafElmntRemoval,WoodyElmntRemoval,&
    FineNonleafElmntRemoval,NonstructElmnt2Litr,NonstructElmntOffEcosystem,&
    LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,WoodyElmntOffEcosystem,&
    StandeadElmntOffEcosystem,LeafElmnt2Litr,FineNonleafElmnt2Litr,&
    WoodyElmnt2Litr,StandeadElmnt2Litr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: EHVST21,EHVST22, EHVST23, EHVST24
  real(r8), intent(in) :: StandeadElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmntRemoval(NumPlantChemElms)  
  real(r8), intent(out) :: NonstructElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmnt2Litr(NumPlantChemElms)

  associate(                                       &
    iHarvstType_pft  => plt_distb%iHarvstType_pft, &
    FracBiomHarvsted => plt_distb%FracBiomHarvsted &
  )
  NonstructElmnt2Litr(ielmc) = NonstructElmntRemoval(ielmc)*EHVST21
  NonstructElmnt2Litr(ielmn) = NonstructElmntRemoval(ielmn)*&
    (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
  NonstructElmnt2Litr(ielmp)=NonstructElmntRemoval(ielmp)*&
    (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
  NonstructElmntOffEcosystem(ielmn) = NonstructElmntRemoval(ielmn)*EHVST21
  NonstructElmntOffEcosystem(ielmp) = NonstructElmntRemoval(ielmp)*EHVST21

  LeafElmnt2Litr(ielmc)        = LeafElmntRemoval(ielmc)*EHVST21
  LeafElmnt2Litr(ielmn)        = LeafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
  LeafElmnt2Litr(ielmp)        = LeafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_leaf,NZ))
  LeafElmntOffEcosystem(ielmn) = LeafElmntRemoval(ielmn)*EHVST21
  LeafElmntOffEcosystem(ielmp) = LeafElmntRemoval(ielmp)*EHVST21

  FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
  FineNonleafElmnt2Litr(ielmn)=FineNonleafElmntRemoval(ielmn)*&
    (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_finenonleaf,NZ))
  FineNonleafElmnt2Litr(ielmp)=FineNonleafElmntRemoval(ielmp)*&
    (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_finenonleaf,NZ))
  FineNonleafElmOffEcosystem(ielmn)=FineNonleafElmntRemoval(ielmn)*EHVST22
  FineNonleafElmOffEcosystem(ielmp)=FineNonleafElmntRemoval(ielmp)*EHVST22

  WoodyElmnt2Litr(ielmc)        = WoodyElmntRemoval(ielmc)*EHVST23
  WoodyElmnt2Litr(ielmn)        = WoodyElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_woody,NZ))
  WoodyElmnt2Litr(ielmp)        = WoodyElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_woody,NZ))
  WoodyElmntOffEcosystem(ielmn) = WoodyElmntRemoval(ielmn)*EHVST23
  WoodyElmntOffEcosystem(ielmp) = WoodyElmntRemoval(ielmp)*EHVST23

  StandeadElmnt2Litr(ielmc) = StandeadElmntRemoval(ielmc)*EHVST24
  StandeadElmnt2Litr(ielmn) = StandeadElmntRemoval(ielmn)*&
    (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_stdead,NZ))
  StandeadElmnt2Litr(ielmp)=StandeadElmntRemoval(ielmp)*&
    (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(2,iplthvst_stdead,NZ))
  StandeadElmntOffEcosystem(ielmn) = StandeadElmntRemoval(ielmn)*EHVST24
  StandeadElmntOffEcosystem(ielmp) = StandeadElmntRemoval(ielmp)*EHVST24
  end associate
  end subroutine ApplyBiomRemovalByFire

  end module PlantDisturbByFireMod