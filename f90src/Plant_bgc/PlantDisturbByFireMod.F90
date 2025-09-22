module PlantDisturbByFireMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ElmIDMod
  use PlantAPIData
  use PlantBGCPars
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
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine InitPlantFireMod
  implicit none
  EFIRE=reshape(real((/0.917,0.167/),r8),shape(EFIRE))
  end subroutine InitPlantFireMod


!----------------------------------------------------------------------------------------------------
  subroutine StageRootRemovalByFire(I,J,NZ,L,FFIRE,DCORP,FracLeftThin)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  integer, intent(in) :: L
  real(r8), intent(in):: DCORP
  real(r8),intent(out) :: FFIRE(NumPlantChemElms)
  real(r8), intent(out):: FracLeftThin

  associate(                                                 &
    iSoilDisturbType_col => plt_distb%iSoilDisturbType_col  ,& !input  :soil disturbance type, [-]
    CSoilOrgM_vr         => plt_soilchem%CSoilOrgM_vr       ,& !input  :soil organic C content, [gC kg soil-1]
    FracBiomHarvsted     => plt_distb%FracBiomHarvsted      ,& !input  :harvest efficiency, [-]
    iHarvstType_pft      => plt_distb%iHarvstType_pft       ,& !input  :type of harvest,[-]
    THETW_vr             => plt_soilchem%THETW_vr            & !input  :volumetric water content, [m3 m-3]
  )
  IF(THETW_vr(L).GT.VolMaxSoilMoist4Fire .OR. CSoilOrgM_vr(ielmc,L).LE.FORGC .OR. iSoilDisturbType_col.NE.22)THEN
    FracLeftThin              = 1.0_r8
    FFIRE(1:NumPlantChemElms) = 0._r8
  ELSE
    FracLeftThin=1.0_r8-DCORP*FracBiomHarvsted(iHarvst_pft,iplthvst_woody,NZ) &
      *AMIN1(1.0_r8,(CSoilOrgM_vr(ielmc,L)-FORGC)/(orgcden-FORGC))
    FFIRE(ielmc) = FracBiomHarvsted(iHarvst_col,iplthvst_woody,NZ)
    FFIRE(ielmn) = FFIRE(ielmc)*EFIRE(1,iHarvstType_pft(NZ))
    FFIRE(ielmp) = FFIRE(ielmc)*EFIRE(2,iHarvstType_pft(NZ))
  ENDIF
  end associate
  end subroutine StageRootRemovalByFire

!----------------------------------------------------------------------------------------------------
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
  real(r8) :: dFireE(NumPlantChemElms),dCO2,dCH4,dNH3
  integer :: NE

  associate(                                                     &
    CO2NetFix_pft        => plt_bgcr%CO2NetFix_pft,              & !inoput :canopy net CO2 exchange,      [gC d-2 h-1]
    PO4byFire_CumYr_pft  => plt_distb%PO4byFire_CumYr_pft,       & !inoput :plant PO4 emission from fire, [g d-2 ]
    N2ObyFire_CumYr_pft  => plt_distb%N2ObyFire_CumYr_pft,       & !inoput :plant N2O emission from fire, [g d-2 ]
    NH3byFire_CumYr_pft  => plt_distb%NH3byFire_CumYr_pft,       & !inoput :plant NH3 emission from fire, [g d-2 ]
    O2ByFire_CumYr_pft   => plt_distb%O2ByFire_CumYr_pft,        & !inoput :plant O2 uptake from fire,    [g d-2 ]
    CH4ByFire_CumYr_pft  => plt_distb%CH4ByFire_CumYr_pft,       & !inoput :plant CH4 emission from fire, [g d-2 ]
    CO2ByFire_CumYr_pft  => plt_distb%CO2ByFire_CumYr_pft,       & !inoput :plant CO2 emission from fire, [g d-2 ]
    Eco_NBP_CumYr_col    => plt_bgcr%Eco_NBP_CumYr_col           & !inoput :total NBP, >0 into ecosystem,  [g d-2]
  )
  DO NE=1,NumPlantChemElms
    dFireE(NE)    = FFIRE(NE)*FrcMassNotHarvst(NE)
  ENDDO

  dCO2                    = (1._r8-FrcAsCH4byFire)*dFireE(ielmc)
  dCH4                    = FrcAsCH4byFire*dFireE(ielmc)
  dNH3                    = dFireE(ielmn)
  CO2ByFire_CumYr_pft(NZ) = CO2ByFire_CumYr_pft(NZ)-dCO2
  CH4ByFire_CumYr_pft(NZ) = CH4ByFire_CumYr_pft(NZ)-dCH4
  O2ByFire_CumYr_pft(NZ)  = O2ByFire_CumYr_pft(NZ)-dCO2*2.667_r8
  NH3byFire_CumYr_pft(NZ) = NH3byFire_CumYr_pft(NZ)-dNH3
  N2ObyFire_CumYr_pft(NZ) = N2ObyFire_CumYr_pft(NZ)-0.0_r8
  PO4byFire_CumYr_pft(NZ) = PO4byFire_CumYr_pft(NZ)-dFireE(ielmp)
  CO2NetFix_pft(NZ)       = CO2NetFix_pft(NZ)-dCO2
  !should CO2 be removed here?
  Eco_NBP_CumYr_col       = Eco_NBP_CumYr_col-dCH4
  end associate
  end subroutine RemoveRootByFire

!----------------------------------------------------------------------------------------------------
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
  real(r8) :: dDeadE

  associate(                                                             &
    PlantElmAllocMat4Litr      => plt_soilchem%PlantElmAllocMat4Litr    ,& !input  :litter kinetic fraction, [-]
    icwood                     => pltpar%icwood                         ,& !input  :group id of coarse woody litter
    istalk                     => pltpar%istalk                         ,& !input  :group id of plant stalk litter group
    ilignin                    => pltpar%ilignin                        ,& !input  :kinetic id of litter component as lignin
    inonstruct                 => pltpar%inonstruct                     ,& !input  :group id of plant nonstructural litter
    k_fine_litr                => pltpar%k_fine_litr                    ,& !input  :fine litter complex id
    k_woody_litr               => pltpar%k_woody_litr                   ,& !input  :woody litter complex id
    ifoliar                    => pltpar%ifoliar                        ,& !input  :group id of plant foliar litter
    inonfoliar                 => pltpar%inonfoliar                     ,& !input  :group id of plant non-foliar litter group
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr  ,& !input  :woody element allocation,[-]
    iPlantTurnoverPattern_pft  => plt_pheno%iPlantTurnoverPattern_pft   ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iPlantRootProfile_pft      => plt_pheno%iPlantRootProfile_pft       ,& !input  :plant growth type (vascular, non-vascular),[-]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr             ,& !inoput :plant LitrFall element, [g d-2 h-1]
    StandDeadKCompElms_pft     => plt_biom%StandDeadKCompElms_pft        & !inoput :standing dead element fraction, [g d-2]
  )

  D6485: DO M=1,jsken
    LitrfallElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
      +PlantElmAllocMat4Litr(ielmc,inonstruct,M,NZ)*NonstructElmnt2Litr(ielmc) &
      +PlantElmAllocMat4Litr(ielmc,ifoliar,M,NZ)   *LeafElmnt2Litr(ielmc)+LeafElmntHarv2Litr(ielmc) &
      +PlantElmAllocMat4Litr(ielmc,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(ielmc)+PetioleElmntHarv2Litr(ielmc))

    DO NE=2,NumPlantChemElms
      LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ) &
        +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*NonstructElmntOffEcosystem(NE) &
        +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)   *LeafElmntOffEcosystem(NE) &
        +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*FineNonleafElmOffEcosystem(NE)

      LitrfallElms_pvr(NE,ilignin,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_fine_litr,0,NZ) &
        +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*(NonstructElmnt2Litr(NE)-NonstructElmntOffEcosystem(NE)) &
        +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)   *(LeafElmnt2Litr(NE)+LeafElmntHarv2Litr(NE)-LeafElmntOffEcosystem(NE)) &
        +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*(FineNonleafElmnt2Litr(NE)+PetioleElmntHarv2Litr(NE) &
        -FineNonleafElmOffEcosystem(NE))
    ENDDO

    !all above ground is subject to fire
    IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
      LitrfallElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_litr,0,NZ)+&
        PlantElmAllocMat4Litr(ielmc,istalk,M,NZ)*(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc) &
        +StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))

      DO NE=2,NumPlantChemElms
        LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ)+&
          PlantElmAllocMat4Litr(NE,istalk,M,NZ)*(WoodyElmntOffEcosystem(NE)+StandeadElmntOffEcosystem(NE))

        LitrfallElms_pvr(NE,ilignin,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_fine_litr,0,NZ) &
          +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE) &
          -WoodyElmntOffEcosystem(NE)+StandeadElmnt2Litr(NE)+StandeadElmntHarv2Litr(NE)&
          -StandeadElmntOffEcosystem(NE))
      ENDDO
      !not grass, has woody component 
    ELSE
      dDeadE=PlantElmAllocMat4Litr(ielmc,icwood,M,NZ)*(WoodyElmnt2Litr(ielmc)+WoodyElmntHarv2Litr(ielmc))
      StandDeadKCompElms_pft(ielmc,M,NZ)=StandDeadKCompElms_pft(ielmc,M,NZ)+dDeadE

      DO NE=2,NumPlantChemElms  
        dDeadE                          = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*WoodyElmntOffEcosystem(NE)
        StandDeadKCompElms_pft(NE,M,NZ) = StandDeadKCompElms_pft(NE,M,NZ)+dDeadE
      ENDDO
        
      LitrfallElms_pvr(ielmc,M,k_woody_litr,0,NZ)=LitrfallElms_pvr(ielmc,M,k_woody_litr,0,NZ) &
        *PlantElmAllocMat4Litr(ielmc,istalk,M,NZ)&
        *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FracWoodStalkElmAlloc2Litr(ielmc,k_woody_litr)

      DO NE=2,NumPlantChemElms  
        LitrfallElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_litr,0,NZ) &
          +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*StandeadElmntOffEcosystem(NE)*FracWoodStalkElmAlloc2Litr(NE,k_woody_litr)

        LitrfallElms_pvr(NE,ilignin,k_woody_litr,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_woody_litr,0,NZ) &
          +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)-WoodyElmntOffEcosystem(NE) &
          +StandeadElmnt2Litr(NE)+StandeadElmntHarv2Litr(NE)-StandeadElmntOffEcosystem(NE)) &
          *FracWoodStalkElmAlloc2Litr(NE,k_woody_litr)
      ENDDO    

      LitrfallElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
        +PlantElmAllocMat4Litr(ielmc,istalk,M,NZ) &
        *(StandeadElmnt2Litr(ielmc)+StandeadElmntHarv2Litr(ielmc))*FracWoodStalkElmAlloc2Litr(ielmc,k_fine_litr)

      DO NE=2,NumPlantChemElms    
        LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ) &
          +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*StandeadElmntOffEcosystem(NE)*FracWoodStalkElmAlloc2Litr(NE,k_fine_litr)
        
        LitrfallElms_pvr(NE,ilignin,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_fine_litr,0,NZ) &
          +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*(WoodyElmnt2Litr(NE)+WoodyElmntHarv2Litr(NE)-WoodyElmntOffEcosystem(NE) &
          +StandeadElmnt2Litr(NE)+StandeadElmntHarv2Litr(NE)-StandeadElmntOffEcosystem(NE)) &
          *FracWoodStalkElmAlloc2Litr(NE,k_fine_litr)

      ENDDO  
    ENDIF
  ENDDO D6485
  end associate
  end subroutine AbvGrndLiterFallByFire

!----------------------------------------------------------------------------------------------------
  subroutine AbvgBiomRemovalByFire(I,J,NZ,TotalElmnt2Litr,TotalElmntRemoval)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: TotalElmntRemoval(NumPlantChemElms)
  real(r8) :: dFireE(NumPlantChemElms)
  integer :: NE

  associate(                                                      &
    Eco_NBP_CumYr_col        => plt_bgcr%Eco_NBP_CumYr_col       ,& !inoput :total NBP, [g d-2]
    CO2NetFix_pft            => plt_bgcr%CO2NetFix_pft           ,& !inoput :canopy net CO2 exchange, [gC d-2 h-1]
    NH3byFire_CumYr_pft      => plt_distb%NH3byFire_CumYr_pft    ,& !inoput :plant NH3 emission from fire, [g d-2 ]
    PO4byFire_CumYr_pft      => plt_distb%PO4byFire_CumYr_pft    ,& !inoput :plant PO4 emission from fire, [g d-2 ]
    CH4ByFire_CumYr_pft      => plt_distb%CH4ByFire_CumYr_pft    ,& !inoput :plant CH4 emission from fire, [g d-2 ]
    O2ByFire_CumYr_pft       => plt_distb%O2ByFire_CumYr_pft     ,& !inoput :plant O2 uptake from fire, [g d-2 ]
    N2ObyFire_CumYr_pft      => plt_distb%N2ObyFire_CumYr_pft    ,& !inoput :plant N2O emission from fire, [g d-2 ]
    CO2ByFire_CumYr_pft      => plt_distb%CO2ByFire_CumYr_pft     & !inoput :plant CO2 emission from fire, [g d-2 ]
  )
  DO NE=1,NumPlantChemElms
    dFireE(NE)=TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
  ENDDO

  CO2ByFire_CumYr_pft(NZ) = CO2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*dFireE(ielmc)
  CH4ByFire_CumYr_pft(NZ) = CH4ByFire_CumYr_pft(NZ)-FrcAsCH4byFire*dFireE(ielmc)
  O2ByFire_CumYr_pft(NZ)  = O2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*dFireE(ielmc)*2.667_r8
  NH3byFire_CumYr_pft(NZ) = NH3byFire_CumYr_pft(NZ)-dFireE(ielmn)
  N2ObyFire_CumYr_pft(NZ) = N2ObyFire_CumYr_pft(NZ)-0.0_r8
  PO4byFire_CumYr_pft(NZ) = PO4byFire_CumYr_pft(NZ)-dFireE(ielmp)
  CO2NetFix_pft(NZ)       = CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*dFireE(ielmc)
  Eco_NBP_CumYr_col       = Eco_NBP_CumYr_col-FrcAsCH4byFire*dFireE(ielmc)
  end associate
  end subroutine AbvgBiomRemovalByFire

!----------------------------------------------------------------------------------------------------
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

  associate(                                         &
    iHarvstType_pft  => plt_distb%iHarvstType_pft   ,& !input  :type of harvest,[-]
    FracBiomHarvsted => plt_distb%FracBiomHarvsted   & !input  :harvest efficiency, [-]
  )

  NonstructElmnt2Litr(ielmc) = NonstructElmntRemoval(ielmc)*EHVST21
  NonstructElmnt2Litr(ielmn) = NonstructElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))
  NonstructElmnt2Litr(ielmp) = NonstructElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))

  NonstructElmntOffEcosystem(ielmn) = NonstructElmntRemoval(ielmn)*EHVST21
  NonstructElmntOffEcosystem(ielmp) = NonstructElmntRemoval(ielmp)*EHVST21

  LeafElmnt2Litr(ielmc)        = LeafElmntRemoval(ielmc)*EHVST21
  LeafElmnt2Litr(ielmn)        = LeafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))
  LeafElmnt2Litr(ielmp)        = LeafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))
  LeafElmntOffEcosystem(ielmn) = LeafElmntRemoval(ielmn)*EHVST21
  LeafElmntOffEcosystem(ielmp) = LeafElmntRemoval(ielmp)*EHVST21

  FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
  FineNonleafElmnt2Litr(ielmn)=FineNonleafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ))
  FineNonleafElmnt2Litr(ielmp)=FineNonleafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ))

  FineNonleafElmOffEcosystem(ielmn)=FineNonleafElmntRemoval(ielmn)*EHVST22
  FineNonleafElmOffEcosystem(ielmp)=FineNonleafElmntRemoval(ielmp)*EHVST22

  WoodyElmnt2Litr(ielmc)        = WoodyElmntRemoval(ielmc)*EHVST23
  WoodyElmnt2Litr(ielmn)        = WoodyElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_woody,NZ))
  WoodyElmnt2Litr(ielmp)        = WoodyElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_woody,NZ))
  WoodyElmntOffEcosystem(ielmn) = WoodyElmntRemoval(ielmn)*EHVST23
  WoodyElmntOffEcosystem(ielmp) = WoodyElmntRemoval(ielmp)*EHVST23

  StandeadElmnt2Litr(ielmc) = StandeadElmntRemoval(ielmc)*EHVST24
  StandeadElmnt2Litr(ielmn) = StandeadElmntRemoval(ielmn)*&
    (1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ))
  StandeadElmnt2Litr(ielmp)=StandeadElmntRemoval(ielmp)*&
    (1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ))
  StandeadElmntOffEcosystem(ielmn) = StandeadElmntRemoval(ielmn)*EHVST24
  StandeadElmntOffEcosystem(ielmp) = StandeadElmntRemoval(ielmp)*EHVST24
  end associate
  end subroutine ApplyBiomRemovalByFire
  ![tail]
  end module PlantDisturbByFireMod