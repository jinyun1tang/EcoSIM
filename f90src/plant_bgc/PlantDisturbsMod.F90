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
  real(r8) :: WTHTH0E(NumOfPlantChemElmnts)
  real(r8) :: WTHTH1E(NumOfPlantChemElmnts)
  real(r8) :: WTHTH2E(NumOfPlantChemElmnts)
  real(r8) :: WTHTH3E(NumOfPlantChemElmnts)
  real(r8) :: WTHTH4E(NumOfPlantChemElmnts)
  real(r8) :: WTHTR1E(NumOfPlantChemElmnts)
  real(r8) :: WTHTR2E(NumOfPlantChemElmnts)
  real(r8) :: WTHTR3E(NumOfPlantChemElmnts)
  real(r8) :: WTHTR4E(NumOfPlantChemElmnts)
  real(r8) :: WTHTX0E(NumOfPlantChemElmnts)
  real(r8) :: WTHTX1E(NumOfPlantChemElmnts)
  real(r8) :: WTHTX2E(NumOfPlantChemElmnts)
  real(r8) :: WTHTX3E(NumOfPlantChemElmnts)
  real(r8) :: WTHTX4E(NumOfPlantChemElmnts)
  real(r8) :: WTHTGE(NumOfPlantChemElmnts)
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

  subroutine RemoveBiomByManagement(I,J,NZ,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(inout) :: CPOOLK(NumOfCanopyLayers1,JP1)
!     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
!

  call RemoveBiomByHarvest(I,J,NZ,CPOOLK)
!
!     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
!
  call RemoveBiomByTillage(I,J,NZ,CPOOLK)
  end subroutine RemoveBiomByManagement
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomassByDisturbance(I,J,NZ,CPOOLK)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), INTENT(INOUT) :: CPOOLK(NumOfCanopyLayers1,JP1)
  real(r8) :: FHVSE(ielmc)
  real(r8) :: FHVSH
  real(r8) :: WHVSTD
  integer :: M

!     begin_execution
  associate(                            &
    EHVST    =>  plt_distb%EHVST  , &
    HVST     =>  plt_distb%HVST   , &
    IHVST    =>  plt_distb%IHVST  , &
    THIN_pft     =>  plt_distb%THIN_pft   , &
    pftPlantPopulation       =>  plt_site%pftPlantPopulation      , &
    NU       =>  plt_site%NU      , &
    ZERO     =>  plt_site%ZERO    , &
    ZNOON    =>  plt_site%ZNOON   , &
    AREA3    =>  plt_site%AREA3   , &
    ZEROQ    =>  plt_rbgc%ZEROQ   , &
    ZEROP    =>  plt_biom%ZEROP   , &
    ZEROL    =>  plt_biom%ZEROL   , &
    StandingDeadKCompChemElmnts_pft   =>  plt_biom%StandingDeadKCompChemElmnts_pft  , &
    StandingDeadChemElmnts_pft   =>  plt_biom%StandingDeadChemElmnts_pft    &
  )

  WTHTH0E(ielmc)=0._r8;WTHTH0E(ielmn)=0._r8;WTHTH0E(ielmp)=0._r8
  WTHTH1E(ielmc)=0._r8;WTHTH1E(ielmn)=0._r8;WTHTH1E(ielmp)=0._r8
  WTHTH2E(ielmc)=0._r8;WTHTH2E(ielmn)=0._r8;WTHTH2E(ielmp)=0._r8
  WTHTH3E(ielmc)=0._r8;WTHTH3E(ielmn)=0._r8;WTHTH3E(ielmp)=0._r8
  WTHTH4E(ielmc)=0._r8;WTHTH4E(ielmn)=0._r8;WTHTH4E(ielmp)=0._r8
  WTHTR1E(ielmc)=0._r8;WTHTR1E(ielmn)=0._r8;WTHTR1E(ielmp)=0._r8
  WTHTR2E(ielmc)=0._r8;WTHTR2E(ielmn)=0._r8;WTHTR2E(ielmp)=0._r8
  WTHTR3E(ielmc)=0._r8;WTHTR3E(ielmn)=0._r8;WTHTR3E(ielmp)=0._r8
  WTHTR4E(ielmc)=0._r8;WTHTR4E(ielmn)=0._r8;WTHTR4E(ielmp)=0._r8
  WTHTX0E(ielmc)=0._r8;WTHTX0E(ielmn)=0._r8;WTHTX0E(ielmp)=0._r8
  WTHTX1E(ielmc)=0._r8;WTHTX1E(ielmn)=0._r8;WTHTX1E(ielmp)=0._r8
  WTHTX2E(ielmc)=0._r8;WTHTX2E(ielmn)=0._r8;WTHTX2E(ielmp)=0._r8
  WTHTX3E(ielmc)=0._r8;WTHTX3E(ielmn)=0._r8;WTHTX3E(ielmp)=0._r8
  WTHTX4E(ielmc)=0._r8;WTHTX4E(ielmn)=0._r8;WTHTX4E(ielmp)=0._r8
  WTHTGE(ielmc)=0._r8;WTHTGE(ielmn)=0._r8;WTHTGE(ielmp)=0._r8
!     ENDIF

!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     THIN_pft=thinning:fraction of population removed
!     FHVSE(ielmc)=fraction of standing dead mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     WTHTH4E(ielmc),WTHTH4E(ielmn),WTHTH4E(ielmp)=harvested standing dead C,N,P
!     WTHTX4E(ielmc),WTHTX4E(ielmn),WTHTX4E(ielmp)=harvested standing dead C,N,P to litter
!
  IF(IHVST(NZ).GE.0)THEN
    IF(J.EQ.INT(ZNOON).AND.IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FHVSE(ielmc)=AZMAX1(1._r8-EHVST(1,4,NZ))
        FHVSH=FHVSE(ielmc)
      ELSE
        FHVSE(ielmc)=AZMAX1(1._r8-THIN_pft(NZ))
        IF(IHVST(NZ).EQ.0)THEN
          FHVSH=AZMAX1(1._r8-EHVST(1,4,NZ)*THIN_pft(NZ))
        ELSE
          FHVSH=FHVSE(ielmc)
        ENDIF
      ENDIF
    ELSEIF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
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
      WTHTH4E(ielmc)=WTHTH4E(ielmc)+(1._r8-FHVSH)*StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)
      WTHTH4E(ielmn)=WTHTH4E(ielmn)+(1._r8-FHVSH)*StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)
      WTHTH4E(ielmp)=WTHTH4E(ielmp)+(1._r8-FHVSH)*StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)
      WTHTX4E(ielmc)=WTHTX4E(ielmc)+(FHVSH-FHVSE(ielmc))*StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)
      WTHTX4E(ielmn)=WTHTX4E(ielmn)+(FHVSH-FHVSE(ielmc))*StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)
      WTHTX4E(ielmp)=WTHTX4E(ielmp)+(FHVSH-FHVSE(ielmc))*StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)=FHVSE(ielmc)*StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)=FHVSE(ielmc)*StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)=FHVSE(ielmc)*StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)
    ENDDO D6475
!
    call PlantDisturbance(I,J,NZ)

    ZEROP(NZ)=ZERO*pftPlantPopulation(NZ)
    ZEROQ(NZ)=ZERO*pftPlantPopulation(NZ)/AREA3(NU)
    ZEROL(NZ)=ZERO*pftPlantPopulation(NZ)*1.0E+06_r8
  ENDIF
  end associate
  end subroutine RemoveBiomassByDisturbance


!------------------------------------------------------------------------------------------
  subroutine PlantDisturbance(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  real(r8) :: WTHEL0(NumOfPlantChemElmnts)
  real(r8) :: WTHEL1(NumOfPlantChemElmnts)
  real(r8) :: WTHEL2(NumOfPlantChemElmnts)
  real(r8) :: WTHEL3(NumOfPlantChemElmnts)
  real(r8) :: WTHEL4(NumOfPlantChemElmnts)
  real(r8) :: WTHER0(NumOfPlantChemElmnts)
  real(r8) :: WTHERT(NumOfPlantChemElmnts)
  real(r8) :: WTHEXT(NumOfPlantChemElmnts)

  WTHEL0(1:NumOfPlantChemElmnts)=0._r8
  WTHEL1(1:NumOfPlantChemElmnts)=0._r8
  WTHEL2(1:NumOfPlantChemElmnts)=0._r8
  WTHEL3(1:NumOfPlantChemElmnts)=0._r8
  WTHEL4(1:NumOfPlantChemElmnts)=0._r8
  WTHER0(1:NumOfPlantChemElmnts)=0._r8

  call ApplyDisturbanceBiomRemoval(I,J,NZ,WTHER0,WTHEL0,WTHEL1,WTHEL2,WTHEL3,WTHEL4)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call TotalBiomRemovalByDisturbance(I,J,NZ,WTHER0,WTHEXT,WTHERT)
!
!     ABOVE-GROUND LITTERFALL FROM HARVESTING
!
  call LiterfallByDisturbance(I,J,NZ,WTHEXT,WTHERT,WTHER0,WTHEL0,WTHEL1,WTHEL2,WTHEL3,WTHEL4)
  end subroutine PlantDisturbance
!------------------------------------------------------------------------------------------

  subroutine LiterfallByDisturbance(I,J,NZ,WTHEXT,WTHERT,WTHER0,&
    WTHEL0,WTHEL1,WTHEL2,WTHEL3,WTHEL4)

  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: WTHEXT(NumOfPlantChemElmnts),WTHERT(NumOfPlantChemElmnts)
  real(r8), intent(in) :: WTHER0(NumOfPlantChemElmnts)
  real(r8), intent(in) :: WTHEL0(NumOfPlantChemElmnts)
  real(r8), intent(in) :: WTHEL1(NumOfPlantChemElmnts)
  real(r8), intent(in) :: WTHEL2(NumOfPlantChemElmnts)
  real(r8), intent(in) :: WTHEL3(NumOfPlantChemElmnts)
  real(r8), intent(in) :: WTHEL4(NumOfPlantChemElmnts)
  integer :: M,NE
!     begin_execution
  associate(                        &
    ilignin  =>  pltpar%ilignin   , &
    instruct =>  pltpar%instruct  , &
    k_fine_litr=> pltpar%k_fine_litr,&
    k_woody_litr=> pltpar%k_woody_litr,&
    ifoliar  =>  pltpar%ifoliar   , &
    infoliar =>  pltpar%infoliar  , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    icwood   =>  pltpar%icwood    , &
    StandingDeadKCompChemElmnts_pft   =>  plt_biom%StandingDeadKCompChemElmnts_pft  , &
    IHVST    =>  plt_distb%IHVST  , &
    FWOODE   =>  plt_allom%FWOODE , &
    CFOPE    =>  plt_soilchem%CFOPE, &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr    , &
    LitrfallChemElmnts_pft    =>  plt_bgcr%LitrfallChemElmnts_pft   , &
    SurfLitrfallChemElmnts_pft    =>  plt_bgcr%SurfLitrfallChemElmnts_pft   , &
    iPlantTurnoverPattern   =>  plt_pheno%iPlantTurnoverPattern , &
    iPlantMorphologyType   =>  plt_pheno%iPlantMorphologyType   &
  )
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTHER0(ielmc),WTHER0(ielmn),WTHER0(ielmp)=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     WTHTX1E(ielmc),WTHTX1E(ielmn),WTHTX1E(ielmp)=harvested leaf C,N,P to litter
!     WTHTX2E(ielmc),WTHTX2E(ielmn),WTHTX2E(ielmp)=harvested petiole C,N,P to litter
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested woody C,N,P to litter
!     WTHTX4E(ielmc),WTHTX4E(ielmn),WTHTX4E(ielmp)=harvested standing dead C,N,P to litter
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      D6375: DO M=1,jsken
        DO NE=1,NumOfPlantChemElmnts        
          LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,instruct,M,NZ)*(WTHER0(NE)+WTHTX0E(NE)) &
            +CFOPE(NE,ifoliar,M,NZ)*(WTHTR1E(NE)+WTHTX1E(NE)) &
            +CFOPE(NE,infoliar,M,NZ)*(WTHTR2E(NE)+WTHTX2E(NE))

          IF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ))))THEN
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,istalk,M,NZ)*(WTHTR3E(NE)+WTHTX3E(NE)+WTHTR4E(NE)+WTHTX4E(NE))
          ELSE
            StandingDeadKCompChemElmnts_pft(NE,M,NZ)=StandingDeadKCompChemElmnts_pft(NE,M,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WTHTX3E(NE)+WTHTX4E(NE))

            LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WTHTR3E(NE)+WTHTR4E(NE))*FWOODE(NE,k_woody_litr)

            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WTHTR3E(NE)+WTHTR4E(NE))*FWOODE(NE,k_fine_litr)
          ENDIF
        ENDDO
      ENDDO D6375

!
!     ABOVE-GROUND LITTERFALL FROM FIRE
!
!     WTHER0(ielmc),WTHER0(ielmn),WTHER0(ielmp)=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     WTHTX1E(ielmc),WTHTX1E(ielmn),WTHTX1E(ielmp)=harvested leaf C,N,P to litter
!     WTHTX2E(ielmc),WTHTX2E(ielmn),WTHTX2E(ielmp)=harvested petiole C,N,P to litter
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested woody C,N,P to litter
!     WTHTX4E(ielmc),WTHTX4E(ielmn),WTHTX4E(ielmp)=harvested standing dead C,N,P to litter
!     WTHEL0(ielmn),WTHEL0(ielmp)=nonstructural N,P to litter
!     WTHEL1(ielmn),WTHEL1(ielmp)=leaf N,P to litter
!     WTHEL2(ielmn),WTHEL2(ielmp)=fine,non-leaf N,P to litter
!     WTHEL3(ielmn),WTHEL3(ielmp)=woody N,P to litter
!     WTHEL4(ielmn),WTHEL4(ielmp)=standing dead N,P to litter
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    ELSE
      D6485: DO M=1,jsken
        LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmc,instruct,M,NZ)*(WTHER0(ielmc)+WTHTX0E(ielmc)) &
          +CFOPE(ielmc,ifoliar,M,NZ)*(WTHTR1E(ielmc)+WTHTX1E(ielmc)) &
          +CFOPE(ielmc,infoliar,M,NZ)*(WTHTR2E(ielmc)+WTHTX2E(ielmc))
        LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmn,instruct,M,NZ)*WTHEL0(ielmn) &
          +CFOPE(ielmn,ifoliar,M,NZ)*WTHEL1(ielmn) &
          +CFOPE(ielmn,infoliar,M,NZ)*WTHEL2(ielmn)
        LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmp,instruct,M,NZ)*WTHEL0(ielmp) &
          +CFOPE(ielmp,ifoliar,M,NZ)*WTHEL1(ielmp) &
          +CFOPE(ielmp,infoliar,M,NZ)*WTHEL2(ielmp)

        LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ) &
          +CFOPE(ielmn,instruct,M,NZ)*(WTHER0(ielmn)+WTHTX0E(ielmn)-WTHEL0(ielmn)) &
          +CFOPE(ielmn,ifoliar,M,NZ)*(WTHTR1E(ielmn)+WTHTX1E(ielmn)-WTHEL1(ielmn)) &
          +CFOPE(ielmn,infoliar,M,NZ)*(WTHTR2E(ielmn)+WTHTX2E(ielmn)-WTHEL2(ielmn))
        LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ) &
          +CFOPE(ielmp,instruct,M,NZ)*(WTHER0(ielmp)+WTHTX0E(ielmp)-WTHEL0(ielmp)) &
          +CFOPE(ielmp,ifoliar,M,NZ)*(WTHTR1E(ielmp)+WTHTX1E(ielmp)-WTHEL1(ielmp)) &
          +CFOPE(ielmp,infoliar,M,NZ)*(WTHTR2E(ielmp)+WTHTX2E(ielmp)-WTHEL2(ielmp))

        IF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ))))THEN
          LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)+CFOPE(ielmc,istalk,M,NZ) &
            * (WTHTR3E(ielmc)+WTHTX3E(ielmc)+WTHTR4E(ielmc)+WTHTX4E(ielmc))
          LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)+CFOPE(ielmn,istalk,M,NZ) &
            *(WTHEL3(ielmn)+WTHEL4(ielmn))
          LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)+CFOPE(ielmp,istalk,M,NZ) &
            *(WTHEL3(ielmp)+WTHEL4(ielmp))
          LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ) &
            +CFOPE(ielmn,istalk,M,NZ)*(WTHTR3E(ielmn)+WTHTX3E(ielmn) &
            -WTHEL3(ielmn)+WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHEL4(ielmn))
          LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)+&
            CFOPE(ielmp,istalk,M,NZ)*(WTHTR3E(ielmp)+WTHTX3E(ielmp)- &
            WTHEL3(ielmp)+WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHEL4(ielmp))
        ELSE
          StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)=StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)+CFOPE(ielmc,icwood,M,NZ) &
            *(WTHTR3E(ielmc)+WTHTX3E(ielmc))
          StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)=StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)+CFOPE(ielmn,icwood,M,NZ)*WTHEL3(ielmn)
          StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)=StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)+CFOPE(ielmp,icwood,M,NZ)*WTHEL3(ielmp)
          LitterFallChemElmnt_pftvr(ielmc,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_woody_litr,0,NZ)*CFOPE(ielmc,istalk,M,NZ)&
            *(WTHTR4E(ielmc)+WTHTX4E(ielmc))*FWOODE(ielmc,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmn,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_woody_litr,0,NZ)+CFOPE(ielmn,istalk,M,NZ)*WTHEL4(ielmn)*FWOODE(ielmn,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmp,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_woody_litr,0,NZ)+CFOPE(ielmp,istalk,M,NZ)*WTHEL4(ielmp)*FWOODE(ielmp,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmn,ilignin,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_woody_litr,0,NZ)+CFOPE(ielmn,icwood,M,NZ) &
            *(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHEL3(ielmn) &
            +WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHEL4(ielmn))*FWOODE(ielmn,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmp,ilignin,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_woody_litr,0,NZ)+CFOPE(ielmp,icwood,M,NZ) &
            *(WTHTR3E(ielmp)+WTHTX3E(ielmp)-WTHEL3(ielmp) &
            +WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHEL4(ielmp))*FWOODE(ielmp,k_woody_litr)
          LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)+CFOPE(ielmc,istalk,M,NZ) &
            *(WTHTR4E(ielmc)+WTHTX4E(ielmc))*FWOODE(ielmc,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,M,k_fine_litr,0,NZ)+CFOPE(ielmn,istalk,M,NZ)*WTHEL4(ielmn)*FWOODE(ielmn,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,M,k_fine_litr,0,NZ)+CFOPE(ielmp,istalk,M,NZ)*WTHEL4(ielmp)*FWOODE(ielmp,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmn,ilignin,k_fine_litr,0,NZ)+CFOPE(ielmn,icwood,M,NZ) &
            *(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHEL3(ielmn) &
            +WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHEL4(ielmn))*FWOODE(ielmn,k_fine_litr)
          LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmp,ilignin,k_fine_litr,0,NZ)+CFOPE(ielmp,icwood,M,NZ) &
            *(WTHTR3E(ielmp)+WTHTX3E(ielmp)-WTHEL3(ielmp) &
            +WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHEL4(ielmp))*FWOODE(ielmp,k_fine_litr)
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
    DO NE=1,NumOfPlantChemElmnts
      LitrfallChemElmnts_pft(NE,NZ)=LitrfallChemElmnts_pft(NE,NZ)+WTHERT(NE)+WTHEXT(NE)
      SurfLitrfallChemElmnts_pft(NE,NZ)=SurfLitrfallChemElmnts_pft(NE,NZ)+WTHERT(NE)+WTHEXT(NE)
    ENDDO
  ENDIF
  end associate
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

  subroutine TotalBiomRemovalByDisturbance(I,J,NZ,WTHER0,WTHEXT,WTHERT)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: WTHER0(NumOfPlantChemElmnts)
  real(r8), intent(out):: WTHEXT(NumOfPlantChemElmnts),WTHERT(NumOfPlantChemElmnts)
  real(r8) :: WTHEHT(NumOfPlantChemElmnts)
  integer :: NE
!     begin_execution
  associate(                            &
    IHVST    =>  plt_distb%IHVST  , &
    JHVST    =>  plt_distb%JHVST  , &
    HVSTE    =>  plt_distb%HVSTE  , &
    XHVSTE   =>  plt_distb%XHVSTE , &
    VNH3F    =>  plt_distb%VNH3F  , &
    VPO4F    =>  plt_distb%VPO4F  , &
    CH4ByFire_pft    =>  plt_distb%CH4ByFire_pft  , &
    VOXYF    =>  plt_distb%VOXYF  , &
    VN2OF    =>  plt_distb%VN2OF  , &
    CO2ByFire_pft    =>  plt_distb%CO2ByFire_pft  , &
    CO2NetFix_pft     =>  plt_bgcr%CO2NetFix_pft    , &
    Eco_NBP_col     =>  plt_bgcr%Eco_NBP_col    , &
    Eco_AutoR_col     =>  plt_bgcr%Eco_AutoR_col    , &
    ECO_ER_col     =>  plt_bgcr%ECO_ER_col    , &
    TCO2A    =>  plt_bgcr%TCO2A   , &
    GrossResp_pft    =>  plt_bgcr%GrossResp_pft   , &
    NonstructalChemElmnts_pft    =>  plt_biom%NonstructalChemElmnts_pft     &
  )
!
!     WTHEHT(ielmc),WTHEHT(ielmn),WTHEHT(ielmp)=total C,N,P removed
!     WTHERT(ielmc),WTHERT(ielmn),WTHERT(ielmp)=total C,N,P to litter
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  DO NE=1,NumOfPlantChemElmnts
    WTHEHT(NE)=WTHTH0E(NE)+WTHTH1E(NE)+WTHTH2E(NE)+WTHTH3E(NE)+WTHTH4E(NE)
    WTHERT(NE)=WTHER0(NE)+WTHTR1E(NE)+WTHTR2E(NE)+WTHTR3E(NE)+WTHTR4E(NE)
    WTHEXT(NE)=WTHTX0E(NE)+WTHTX1E(NE)+WTHTX2E(NE)+WTHTX3E(NE)+WTHTX4E(NE)
  ENDDO
  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      IF(JHVST(NZ).NE.2)THEN
        DO NE=1,NumOfPlantChemElmnts
          HVSTE(NE,NZ)=HVSTE(NE,NZ)+WTHEHT(NE)-WTHERT(NE)
          XHVSTE(NE)=XHVSTE(NE)+WTHEHT(NE)-WTHERT(NE)
        ENDDO
        Eco_NBP_col=Eco_NBP_col+WTHERT(ielmc)-WTHEHT(ielmc)
      ELSE
        DO NE=1,NumOfPlantChemElmnts
          NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+WTHEHT(NE)-WTHERT(NE)
        ENDDO
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     CO2ByFire_pft,CH4ByFire_pft,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_col=total net biome productivity
!
    ELSE
      CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*(WTHEHT(ielmc)-WTHERT(ielmc))
      CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*(WTHEHT(ielmc)-WTHERT(ielmc))
      VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*(WTHEHT(ielmc)-WTHERT(ielmc))*2.667
      VNH3F(NZ)=VNH3F(NZ)-WTHEHT(ielmn)+WTHERT(ielmn)
      VN2OF(NZ)=VN2OF(NZ)-0.0_r8
      VPO4F(NZ)=VPO4F(NZ)-WTHEHT(ielmp)+WTHERT(ielmp)
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*(WTHEHT(ielmc)-WTHERT(ielmc))
      Eco_NBP_col=Eco_NBP_col-FCH4F*(WTHEHT(ielmc)-WTHERT(ielmc))
    ENDIF
!
!     C,N,P REMOVED FROM GRAZING
!
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     GY=growth yield of grazers
!     WTHEHT(ielmc),WTHEHT(ielmn),WTHEHT(ielmp)=total C,N,P removed
!     WTHERT(ielmc),WTHERT(ielmn),WTHERT(ielmp)=total C,N,P to litter
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!
  ELSE
    HVSTE(ielmc,NZ)=HVSTE(ielmc,NZ)+GY*(WTHEHT(ielmc)-WTHERT(ielmc))
    XHVSTE(ielmc)=XHVSTE(ielmc)+GY*(WTHEHT(ielmc)-WTHERT(ielmc))
    DO NE=2,NumOfPlantChemElmnts
      HVSTE(NE,NZ)=HVSTE(NE,NZ)+WTHEHT(NE)-WTHERT(NE)
      XHVSTE(NE)=XHVSTE(NE)+WTHEHT(NE)-WTHERT(NE)
    ENDDO
    GrossResp_pft(NZ)=GrossResp_pft(NZ)-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
    TCO2A(NZ)=TCO2A(NZ)-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
!     Eco_NBP_col=Eco_NBP_col+GY*(WTHERT(ielmc)-WTHEHT(ielmc))
!     CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)+GZ*(WTHERT(ielmc)-WTHEHT(ielmc))
    ECO_ER_col=ECO_ER_col-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
    Eco_AutoR_col=Eco_AutoR_col-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
  ENDIF
  end associate
  end subroutine TotalBiomRemovalByDisturbance
!------------------------------------------------------------------------------------------

  subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,WTHER0,&
    WTHEL0,WTHEL1,WTHEL2,WTHEL3,WTHEL4)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: WTHER0(NumOfPlantChemElmnts)
  real(r8), intent(out) :: WTHEL0(NumOfPlantChemElmnts)
  real(r8), intent(out) :: WTHEL1(NumOfPlantChemElmnts)
  real(r8), intent(out) :: WTHEL2(NumOfPlantChemElmnts)
  real(r8), intent(out) :: WTHEL3(NumOfPlantChemElmnts)
  real(r8), intent(out) :: WTHEL4(NumOfPlantChemElmnts)

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
    IHVST    =>  plt_distb%IHVST    &
  )
!     IF NO PLANT C,N,P REMOVED AT HARVEST (ALL RESIDUE RETURNED)
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTHTH0E(ielmc),WTHTH0E(ielmn),WTHTH0E(ielmp)=nonstructural C,N,P removed
!     WTHTH1E(ielmc),WTHTH1E(ielmn),WTHTH1E(ielmp)=leaf C,N,P removed
!     WTHTH2E(ielmc),WTHTH2E(ielmn),WTHTH2E(ielmp)=fine,non-leaf C,N,P removed
!     WTHTH3E(ielmc),WTHTH3E(ielmn),WTHTH3E(ielmp)=woody C,N,P removed
!     WTHTH4E(ielmc),WTHTH4E(ielmn),WTHTH4E(ielmp)=standing dead C,N,P removed
!     WTHER0(ielmc),WTHER0(ielmn),WTHER0(ielmp)=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  EHVST21=1._r8-EHVST(2,ipld_leaf,NZ)
  EHVST22=1._r8-EHVST(2,ipld_nofoliar,NZ)
  EHVST23=1._r8-EHVST(2,ipld_woody,NZ)
  EHVST24=1._r8-EHVST(2,ipld_stdead,NZ)

  IF(IHVST(NZ).EQ.0)THEN
    DO NE=1,NumOfPlantChemElmnts
      WTHER0(NE)=WTHTH0E(NE)*EHVST21    !non-structural
      WTHTR1E(NE)=WTHTH1E(NE)*EHVST21   !leaf
      WTHTR2E(NE)=WTHTH2E(NE)*EHVST22   !fine, non-woody
      WTHTR3E(NE)=WTHTH3E(NE)*EHVST23   !woody
      WTHTR4E(NE)=WTHTH4E(NE)*EHVST24   !standing dead
    ENDDO
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ).EQ.1)THEN
    DO NE=1,NumOfPlantChemElmnts
      WTHER0(NE)=WTHTH0E(NE)
      WTHTR1E(NE)=WTHTH1E(NE)
      WTHTR2E(NE)=WTHTH2E(NE)-WTHTGE(NE)*EHVST(2,ipld_nofoliar,NZ)
      WTHTR3E(NE)=WTHTH3E(NE)
      WTHTR4E(NE)=WTHTH4E(NE)
    ENDDO
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ).EQ.2)THEN
    DO NE=1,NumOfPlantChemElmnts
      WTHER0(NE)=WTHTH0E(NE)*EHVST21
      WTHTR1E(NE)=WTHTH1E(NE)*EHVST21
      WTHTR2E(NE)=WTHTH2E(NE)*EHVST22
      WTHTR3E(NE)=WTHTH3E(NE)*EHVST23
      WTHTR4E(NE)=WTHTH4E(NE)*EHVST24
    ENDDO
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(IHVST(NZ).EQ.3)THEN
    DO NE=1,NumOfPlantChemElmnts
      WTHER0(NE)=WTHTH0E(NE)*EHVST21
      WTHTR1E(NE)=WTHTH1E(NE)*EHVST21
      WTHTR2E(NE)=WTHTH2E(NE)*EHVST22
      WTHTR3E(NE)=WTHTH3E(NE)*EHVST23
      WTHTR4E(NE)=WTHTH4E(NE)*EHVST24
    ENDDO
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
    EHVST21h=1._r8-EHVST(2,ipld_leaf,NZ)*0.5_r8
    EHVST22h=1._r8-EHVST(2,ipld_nofoliar,NZ)*0.5_r8
    EHVST23h=1._r8-EHVST(2,ipld_woody,NZ)*0.5_r8
    EHVST24h=1._r8-EHVST(2,ipld_stdead,NZ)*0.5_r8

    WTHER0(ielmc)=WTHTH0E(ielmc)*EHVST21
    WTHTR1E(ielmc)=WTHTH1E(ielmc)*EHVST21
    WTHTR2E(ielmc)=WTHTH2E(ielmc)*EHVST22
    WTHTR3E(ielmc)=WTHTH3E(ielmc)*EHVST23
    WTHTR4E(ielmc)=WTHTH4E(ielmc)*EHVST24

    DO NE=2,NumOfPlantChemElmnts
      WTHER0(NE)=WTHTH0E(NE)*EHVST21h
      WTHTR1E(NE)=WTHTH1E(NE)*EHVST21h
      WTHTR2E(NE)=WTHTH2E(NE)*EHVST22h
      WTHTR3E(NE)=WTHTH3E(NE)*EHVST23h
      WTHTR4E(NE)=WTHTH4E(NE)*EHVST24h
    ENDDO
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
    FERT(17)=FERT(17)+(WTHER0(ielmc)+WTHTR1E(ielmc)+WTHTR2E(ielmc)+WTHTR3E(ielmc)+WTHTR4E(ielmc))/AREA3(NU)
    FERT(18)=FERT(18)+(WTHER0(ielmn)+WTHTR1E(ielmn)+WTHTR2E(ielmn)+WTHTR3E(ielmn)+WTHTR4E(ielmn))/AREA3(NU)*0.5_r8
    FERT(3)=FERT(3)+(WTHER0(ielmn)+WTHTR1E(ielmn)+WTHTR2E(ielmn)+WTHTR3E(ielmn)+WTHTR4E(ielmn))/AREA3(NU)*0.5_r8
    FERT(19)=FERT(19)+(WTHER0(ielmp)+WTHTR1E(ielmp)+WTHTR2E(ielmp)+WTHTR3E(ielmp)+WTHTR4E(ielmp))/AREA3(NU)
    IYTYP=3
!
!     REMOVALS BY FIRE
!
!     EFIRE=combustion  of N,P relative to C
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     WTHTH0E(ielmc),WTHTH0E(ielmn),WTHTH0E(ielmp)=nonstructural C,N,P removed
!     WTHTH1E(ielmc),WTHTH1E(ielmn),WTHTH1E(ielmp)=leaf C,N,P removed
!     WTHTH2E(ielmc),WTHTH2E(ielmn),WTHTH2E(ielmp)=fine,non-leaf C,N,P removed
!     WTHTH3E(ielmc),WTHTH3E(ielmn),WTHTH3E(ielmp)=woody C,N,P removed
!     WTHTH4E(ielmc),WTHTH4E(ielmn),WTHTH4E(ielmp)=standing dead C,N,P removed
!     WTHER0(ielmc),WTHER0(ielmn),WTHER0(ielmp)=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     WTHTL0,WTHEL0(ielmn),WTHEL0(ielmp)=nonstructural C,N,P removed from ecosystem
!     WTHTL1,WTHEL1(ielmn),WTHEL1(ielmp)=leaf C,N,P removed from ecosystem
!     WTHTL2,WTHEL2(ielmn),WTHEL2(ielmp)=fine,non-leaf C,N,P removed from ecosystem
!     WTHTL3,WTHEL3(ielmn),WTHEL3(ielmp)=woody C,N,P removed from ecosystem
!     WTHTL4,WTHEL4(ielmn),WTHEL4(ielmp)=standing dead C,N,P removed from ecosystem
!
  ELSEIF(IHVST(NZ).EQ.5)THEN

    WTHER0(ielmc)=WTHTH0E(ielmc)*EHVST21
    WTHER0(ielmn)=WTHTH0E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,ipld_leaf,NZ))
    WTHER0(ielmp)=WTHTH0E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,ipld_leaf,NZ))
    WTHEL0(ielmn)=WTHTH0E(ielmn)*EHVST21
    WTHEL0(ielmp)=WTHTH0E(ielmp)*EHVST21

    WTHTR1E(ielmc)=WTHTH1E(ielmc)*EHVST21
    WTHTR1E(ielmn)=WTHTH1E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,ipld_leaf,NZ))
    WTHTR1E(ielmp)=WTHTH1E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,ipld_leaf,NZ))
    WTHEL1(ielmn)=WTHTH1E(ielmn)*EHVST21
    WTHEL1(ielmp)=WTHTH1E(ielmp)*EHVST21

    WTHTR2E(ielmc)=WTHTH2E(ielmc)*EHVST22
    WTHTR2E(ielmn)=WTHTH2E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,ipld_nofoliar,NZ))
    WTHTR2E(ielmp)=WTHTH2E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,ipld_nofoliar,NZ))
    WTHEL2(ielmn)=WTHTH2E(ielmn)*EHVST22
    WTHEL2(ielmp)=WTHTH2E(ielmp)*EHVST22

    WTHTR3E(ielmc)=WTHTH3E(ielmc)*EHVST23
    WTHTR3E(ielmn)=WTHTH3E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,ipld_woody,NZ))
    WTHTR3E(ielmp)=WTHTH3E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,ipld_woody,NZ))
    WTHEL3(ielmn)=WTHTH3E(ielmn)*EHVST23
    WTHEL3(ielmp)=WTHTH3E(ielmp)*EHVST23

    WTHTR4E(ielmc)=WTHTH4E(ielmc)*EHVST24
    WTHTR4E(ielmn)=WTHTH4E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,ipld_stdead,NZ))
    WTHTR4E(ielmp)=WTHTH4E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,ipld_stdead,NZ))
    WTHEL4(ielmn)=WTHTH4E(ielmn)*EHVST24
    WTHEL4(ielmp)=WTHTH4E(ielmp)*EHVST24
  ENDIF
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByTillage(I,J,NZ,CPOOLK)
  use SurfLitterDataType, only : XCORP
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(inout) :: CPOOLK(10,JP1)
  integer :: L,K,M,N,NR,NE,NB,NTG
  real(r8) :: XHVST,XHVST1
  REAL(R8) :: APSILT
  real(r8) :: FDM,VOLWPX
  real(r8) :: WVPLT
!     begin_execution
  associate(                               &
    JHVST    =>  plt_distb%JHVST     , &
    iDayPlantHarvest    =>  plt_distb%iDayPlantHarvest     , &
    iDayPlanting    =>  plt_distb%iDayPlanting     , &
    ITILL    =>  plt_distb%ITILL     , &
    iYearPlanting     =>  plt_distb%iYearPlanting      , &
    iYearPlantHarvest     =>  plt_distb%iYearPlantHarvest      , &
    XCORP    =>  plt_distb%XCORP     , &
    CFOPE    =>  plt_soilchem%CFOPE  , &
    trcg_rootml     =>  plt_rbgc%trcg_rootml       , &
    trcs_rootml => plt_rbgc%trcs_rootml , &
    UVOLO    =>  plt_ew%UVOLO        , &
    CanWatP    =>  plt_ew%CanWatP        , &
    VHeatCapCanP    =>  plt_ew%VHeatCapCanP        , &
    PSICanP    =>  plt_ew%PSICanP        , &
    PPX      =>  plt_site%PPX        , &
     RootMycoNonstructElmnt_vr   =>  plt_biom% RootMycoNonstructElmnt_vr     , &
    WSRTL    =>  plt_biom%WSRTL      , &
     PopuPlantRootC_vr    =>  plt_biom% PopuPlantRootC_vr      , &
    RootStructBiomC_vr   =>  plt_biom%RootStructBiomC_vr     , &
    LeafChemElmnts_brch  =>  plt_biom%LeafChemElmnts_brch    , &
    GrainChemElmnts_brch   =>  plt_biom%GrainChemElmnts_brch     , &
    EarChemElmnts_brch  =>  plt_biom%EarChemElmnts_brch    , &
    NoduleNonstructElmnt_brch   =>  plt_biom%NoduleNonstructElmnt_brch     , &
    NonstructElmnt_brch   =>  plt_biom%NonstructElmnt_brch     , &
    HuskChemElmnts_brch  =>  plt_biom%HuskChemElmnts_brch    , &
    ReserveChemElmnts_brch  =>  plt_biom%ReserveChemElmnts_brch    , &
    WTNDBE   =>  plt_biom%WTNDBE     , &
    WTSHTBE  =>  plt_biom%WTSHTBE    , &
    StalkChemElmnts_brch  =>  plt_biom%StalkChemElmnts_brch    , &
    WGLFLE   =>  plt_biom%WGLFLE     , &
    PetioleChemElmnts_brch =>  plt_biom%PetioleChemElmnts_brch   , &
    LeafChemElmntNode_brch    =>  plt_biom%LeafChemElmntNode_brch      , &
    LeafPetioleBiomassC_brch    =>  plt_biom%LeafPetioleBiomassC_brch      , &
    WTSTXBE  =>  plt_biom%WTSTXBE    , &
    StalkBiomassC_brch   =>  plt_biom%StalkBiomassC_brch     , &
    PetioleProteinCNode_brch   =>  plt_biom%PetioleProteinCNode_brch     , &
    PetioleElmntNode_brch   =>  plt_biom%PetioleElmntNode_brch     , &
    WTRT1E   =>  plt_biom%WTRT1E     , &
    LeafProteinCNode_brch     =>  plt_biom%LeafProteinCNode_brch       , &
    InternodeChemElmnt_brch   =>  plt_biom%InternodeChemElmnt_brch     , &
    CanPStalkC    =>  plt_biom%CanPStalkC      , &
    RTWT1E   =>  plt_biom%RTWT1E     , &
    NonstructalChemElmnts_pft    =>  plt_biom%NonstructalChemElmnts_pft      , &
    CanopyLeafShethC_pft     =>  plt_biom%CanopyLeafShethC_pft       , &
    WTRT2E   =>  plt_biom%WTRT2E     , &
    RootNoduleNonstructElmnt_vr  =>  plt_biom%RootNoduleNonstructElmnt_vr    , &
    WTNDLE   =>  plt_biom%WTNDLE     , &
    GRWTB    =>  plt_allom%GRWTB     , &
    FWOODE   =>  plt_allom%FWOODE    , &
    FWODBE   =>  plt_allom%FWODBE    , &
    FWODLE   =>  plt_allom%FWODLE    , &
    FWODRE   =>  plt_allom%FWODRE    , &
    iPlantBranchState    =>  plt_pheno%iPlantBranchState     , &
    iPlantPhenologyPattern   =>  plt_pheno%iPlantPhenologyPattern    , &
    iPlantState    =>  plt_pheno%iPlantState     , &
    iPlantMorphologyType   =>  plt_pheno%iPlantMorphologyType    , &
    iPlantTurnoverPattern   =>  plt_pheno%iPlantTurnoverPattern    , &
    iPlantPhenologyType   =>  plt_pheno%iPlantPhenologyType    , &
    iPlantRootState   =>  plt_pheno%iPlantRootState    , &
    iPlantShootState    =>  plt_pheno%iPlantShootState     , &
    HCOB     =>  plt_photo%HCOB      , &
    CO2B     =>  plt_photo%CO2B      , &
    CPOOL3   =>  plt_photo%CPOOL3    , &
    CPOOL4   =>  plt_photo%CPOOL4    , &
    NJ       =>  plt_site%NJ         , &
    pftPlantPopulation       =>  plt_site%pftPlantPopulation         , &
    iYearCurrent     =>  plt_site%iYearCurrent       , &
    ZNOON    =>  plt_site%ZNOON      , &
    VOLWOU   =>  plt_site%VOLWOU     , &
    NU       =>  plt_site%NU         , &
    k_fine_litr=> pltpar%k_fine_litr ,&
    k_woody_litr=> pltpar%k_woody_litr,&
    instruct =>  pltpar%instruct     , &
    ifoliar  =>  pltpar%ifoliar      , &
    istalk   =>  pltpar%istalk       , &
    iroot    =>  pltpar%iroot        , &
    infoliar =>  pltpar%infoliar     , &
    icwood   =>  pltpar%icwood       , &
    RootGasLoss_disturb    =>   plt_bgcr%RootGasLoss_disturb    , &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr       , &
    RCO2A    =>  plt_rbgc%RCO2A      , &
    RCO2M    =>  plt_rbgc%RCO2M      , &
    RCO2N    =>  plt_rbgc%RCO2N      , &
    FracPARByCanP    =>  plt_rad%FracPARByCanP       , &
    PrimRootLen    =>  plt_morph%PrimRootLen     , &
    RTVLW    =>  plt_morph%RTVLW     , &
    RTARP    =>  plt_morph%RTARP     , &
    RTVLP    =>  plt_morph%RTVLP     , &
    RootLenDensNLP    =>  plt_morph%RootLenDensNLP     , &
    PrimRootXNumL    =>  plt_morph%PrimRootXNumL     , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType    , &
    RootLenPerP    =>  plt_morph%RootLenPerP     , &
    RTN2     =>  plt_morph%RTN2      , &
    SecndRootLen    =>  plt_morph%SecndRootLen     , &
    SecndRootXNumL     =>  plt_morph%SecndRootXNumL      , &
    NGTopRootLayer      =>  plt_morph%NGTopRootLayer       , &
    MY       =>  plt_morph%MY        , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft      , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft       , &
    LeafAreaNode_brch    =>  plt_morph%LeafAreaNode_brch     , &
    LeafAreaLive_brch    =>  plt_morph%LeafAreaLive_brch     , &
    GRNXB    =>  plt_morph%GRNXB     , &
    GRNOB    =>  plt_morph%GRNOB     , &
    CanPLNBLA    =>  plt_morph%CanPLNBLA       &
  )
!     ZNOON=hour of solar noon
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iDayPlanting,iYearPlanting=day,year of planting
!     iYearCurrent=current year
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!     PPX,PP=PFT population per m2,grid cell
!     FracPARByCanP=fraction of radiation received by each PFT canopy
!     VHeatCapCanP=canopy heat capacity
!
  IF(J.EQ.INT(ZNOON).AND.(iPlantTurnoverPattern(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ)))) &
    .AND.(I.NE.iDayPlanting(NZ) &
    .OR.iYearCurrent.NE.iYearPlanting(NZ)))THEN
    IF(ITILL.LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.iDayPlanting(NZ).OR.iYearCurrent.GT.iYearPlanting(NZ))THEN
        XHVST=XCORP
        PPX(NZ)=PPX(NZ)*XHVST
        pftPlantPopulation(NZ)=pftPlantPopulation(NZ)*XHVST
        FracPARByCanP(NZ)=FracPARByCanP(NZ)*XHVST
        VHeatCapCanP(NZ)=VHeatCapCanP(NZ)*XHVST
        CanopyLeafShethC_pft(NZ)=0._r8
        CanPStalkC(NZ)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
        D8975: DO NB=1,NumOfBranches_pft(NZ)
          IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
            IF(pftPlantPopulation(NZ).LE.0.0)then
              iPlantBranchState(NB,NZ)=iDead
            endif
!
!     LITTERFALL FROM BRANCHES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenologyType=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
            XHVST1=1._r8-XHVST
            D6380: DO M=1,jsken
              LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)+XHVST1 &
                *(CFOPE(ielmc,instruct,M,NZ)*(NonstructElmnt_brch(ielmc,NB,NZ)+NoduleNonstructElmnt_brch(ielmc,NB,NZ) &
                +CPOOLK(NB,NZ)+ReserveChemElmnts_brch(ielmc,NB,NZ)) &
                +CFOPE(ielmc,ifoliar,M,NZ)*(LeafChemElmnts_brch(ielmc,NB,NZ)*FWODLE(ielmc,k_fine_litr) &
                +WTNDBE(ielmc,NB,NZ)) &
                +CFOPE(ielmc,infoliar,M,NZ)*(PetioleChemElmnts_brch(ielmc,NB,NZ)*FWODBE(ielmc,k_fine_litr) &
                +HuskChemElmnts_brch(ielmc,NB,NZ)+EarChemElmnts_brch(ielmc,NB,NZ)))

              DO NE=2,NumOfPlantChemElmnts
                LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *(CFOPE(NE,instruct,M,NZ)*(NonstructElmnt_brch(NE,NB,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ)&
                  +ReserveChemElmnts_brch(NE,NB,NZ)) &
                  + CFOPE(NE,ifoliar,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr)+WTNDBE(NE,NB,NZ)) &
                  + CFOPE(NE,infoliar,M,NZ)*(PetioleChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
                  +   HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ)))
              ENDDO
            ENDDO D6380

            DO M=1,jsken
              DO NE=1,NumOfPlantChemElmnts
                LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,icwood,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
                  +PetioleChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

                IF(iPlantPhenologyPattern(NZ).EQ.iplt_annual.AND.iPlantPhenologyType(NZ).NE.0)THEN
                  NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+XHVST1*CFOPE(NE,infoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
                ELSE
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                    *CFOPE(NE,infoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
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

            CPOOLK(NB,NZ)=CPOOLK(NB,NZ)*XHVST
            StalkBiomassC_brch(NB,NZ)=StalkBiomassC_brch(NB,NZ)*XHVST
            DO NE=1,NumOfPlantChemElmnts
              NonstructElmnt_brch(NE,NB,NZ)=NonstructElmnt_brch(NE,NB,NZ)*XHVST
              NoduleNonstructElmnt_brch(NE,NB,NZ)=NoduleNonstructElmnt_brch(NE,NB,NZ)*XHVST
              WTSHTBE(NE,NB,NZ)=WTSHTBE(NE,NB,NZ)*XHVST
              ReserveChemElmnts_brch(NE,NB,NZ)=ReserveChemElmnts_brch(NE,NB,NZ)*XHVST
              HuskChemElmnts_brch(NE,NB,NZ)=HuskChemElmnts_brch(NE,NB,NZ)*XHVST
              EarChemElmnts_brch(NE,NB,NZ)=EarChemElmnts_brch(NE,NB,NZ)*XHVST
              GrainChemElmnts_brch(NE,NB,NZ)=GrainChemElmnts_brch(NE,NB,NZ)*XHVST
              LeafChemElmnts_brch(NE,NB,NZ)=LeafChemElmnts_brch(NE,NB,NZ)*XHVST
              WTNDBE(NE,NB,NZ)=WTNDBE(NE,NB,NZ)*XHVST
              PetioleChemElmnts_brch(NE,NB,NZ)=PetioleChemElmnts_brch(NE,NB,NZ)*XHVST
              StalkChemElmnts_brch(NE,NB,NZ)=StalkChemElmnts_brch(NE,NB,NZ)*XHVST
              WTSTXBE(NE,NB,NZ)=WTSTXBE(NE,NB,NZ)*XHVST
            ENDDO

            GRNXB(NB,NZ)=GRNXB(NB,NZ)*XHVST
            GRNOB(NB,NZ)=GRNOB(NB,NZ)*XHVST
            GRWTB(NB,NZ)=GRWTB(NB,NZ)*XHVST
            LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)*XHVST
            LeafPetioleBiomassC_brch(NB,NZ)=AZMAX1(LeafChemElmnts_brch(ielmc,NB,NZ)+PetioleChemElmnts_brch(ielmc,NB,NZ))
            CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetioleBiomassC_brch(NB,NZ)

            CanPStalkC(NZ)=CanPStalkC(NZ)+StalkBiomassC_brch(NB,NZ)
            D8970: DO K=0,MaxNodesPerBranch1
              IF(K.NE.0)THEN
                CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)*XHVST
                CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)*XHVST
                CO2B(K,NB,NZ)=CO2B(K,NB,NZ)*XHVST
                HCOB(K,NB,NZ)=HCOB(K,NB,NZ)*XHVST
              ENDIF
              LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)*XHVST

              LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*XHVST
!     PetioleLengthNode_brch(K,NB,NZ)=PetioleLengthNode_brch(K,NB,NZ)*XHVST

              PetioleProteinCNode_brch(K,NB,NZ)=PetioleProteinCNode_brch(K,NB,NZ)*XHVST
!     InternodeHeightLive_brch(K,NB,NZ)=InternodeHeightLive_brch(K,NB,NZ)*XHVST
!     InternodeHeightDying_brch(K,NB,NZ)=InternodeHeightDying_brch(K,NB,NZ)*XHVST
              DO NE=1,NumOfPlantChemElmnts
                InternodeChemElmnt_brch(NE,K,NB,NZ)=InternodeChemElmnt_brch(NE,K,NB,NZ)*XHVST
                LeafChemElmntNode_brch(NE,K,NB,NZ)=LeafChemElmntNode_brch(NE,K,NB,NZ)*XHVST
                PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)*XHVST
                DO L=1,NumOfCanopyLayers1
                  WGLFLE(NE,L,K,NB,NZ)=WGLFLE(NE,L,K,NB,NZ)*XHVST
                ENDDO
              ENDDO
              D8965: DO L=1,NumOfCanopyLayers1
                CanPLNBLA(L,K,NB,NZ)=CanPLNBLA(L,K,NB,NZ)*XHVST
              ENDDO D8965
            ENDDO D8970
          ENDIF
        ENDDO D8975
!
!     PSICanP=canopy water potential
!     CanWatP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        VOLWPX=CanWatP(NZ)
        WVPLT=AZMAX1(CanopyLeafShethC_pft(NZ)+CanPStalkC(NZ))

        FDM=get_FDM(PSICanP(NZ))
!        APSILT=ABS(PSICanP(NZ))
!        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)

        CanWatP(NZ)=ppmc*WVPLT/FDM
        VOLWOU=VOLWOU+VOLWPX-CanWatP(NZ)
        UVOLO=UVOLO+VOLWPX-CanWatP(NZ)
!
!     TERMINATE ROOTS IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     PP=PFT population
!     IDTHR,iPlantShootState=PFT root,shoot living flag: 0=alive,1=dead
!     IDTH=PFT living flag: 0=alive,1=dead
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     iDayPlantHarvest,iYearPlantHarvest=day,year of harvesting
!     iYearCurrent=current year
!
        IF(pftPlantPopulation(NZ).LE.0.0_r8)THEN
          iPlantRootState(NZ)=iDead
          iPlantShootState(NZ)=iDead
          iPlantState(NZ)=iDead
          JHVST(NZ)=ihv_terminate
          iDayPlantHarvest(NZ)=I
          iYearPlantHarvest(NZ)=iYearCurrent
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
        D8980: DO L=NU,NJ
          D8985: DO N=1,MY(NZ)

            D6385: DO M=1,jsken
                DO NE=1,NumOfPlantChemElmnts
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,instruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
                ENDDO

              DO NR=1,NumRootAxes_pft(NZ)
                DO NE=1,NumOfPlantChemElmnts
                  LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,icwood,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,iroot,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
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
              RootGasLoss_disturb(NTG,NZ)=RootGasLoss_disturb(NTG,NZ)-XHVST1 &
                *(trcg_rootml(NTG,N,L,NZ)+trcs_rootml(NTG,N,L,NZ))
              trcg_rootml(NTG,N,L,NZ)=XHVST*trcg_rootml(NTG,N,L,NZ)
              trcs_rootml(NTG,N,L,NZ)=XHVST*trcs_rootml(NTG,N,L,NZ)
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
!     WSRTL=root protein C mass
!     RTN1,SecndRootXNumL=number of primary,secondary root axes
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D8960: DO NR=1,NumRootAxes_pft(NZ)
              DO NE=1,NumOfPlantChemElmnts
                WTRT1E(NE,N,L,NR,NZ)=WTRT1E(NE,N,L,NR,NZ)*XHVST
                WTRT2E(NE,N,L,NR,NZ)=WTRT2E(NE,N,L,NR,NZ)*XHVST
                RTWT1E(NE,N,NR,NZ)=RTWT1E(NE,N,NR,NZ)*XHVST
              ENDDO
              PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)*XHVST
              SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)*XHVST
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST
            ENDDO D8960
            DO NE=1,NumOfPlantChemElmnts
               RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)*XHVST
            ENDDO
            RootStructBiomC_vr(N,L,NZ)=RootStructBiomC_vr(N,L,NZ)*XHVST
             PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)*XHVST
            WSRTL(N,L,NZ)=WSRTL(N,L,NZ)*XHVST
            PrimRootXNumL(N,L,NZ)=PrimRootXNumL(N,L,NZ)*XHVST
            SecndRootXNumL(N,L,NZ)=SecndRootXNumL(N,L,NZ)*XHVST
            RootLenPerP(N,L,NZ)=RootLenPerP(N,L,NZ)*XHVST
            RootLenDensNLP(N,L,NZ)=RootLenDensNLP(N,L,NZ)*XHVST
            RTVLP(N,L,NZ)=RTVLP(N,L,NZ)*XHVST
            RTVLW(N,L,NZ)=RTVLW(N,L,NZ)*XHVST
            RTARP(N,L,NZ)=RTARP(N,L,NZ)*XHVST
            RCO2M(N,L,NZ)=RCO2M(N,L,NZ)*XHVST
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
            IF(iPlantNfixType(NZ).NE.0.AND.N.EQ.1)THEN
              DO NE=1,NumOfPlantChemElmnts
                D6395: DO M=1,jsken
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *(CFOPE(NE,iroot,M,NZ)*WTNDLE(NE,L,NZ) &
                    +CFOPE(NE,instruct,M,NZ)*RootNoduleNonstructElmnt_vr(NE,L,NZ))
                ENDDO D6395
                WTNDLE(NE,L,NZ)=WTNDLE(NE,L,NZ)*XHVST
                RootNoduleNonstructElmnt_vr(NE,L,NZ)=RootNoduleNonstructElmnt_vr(NE,L,NZ)*XHVST
              ENDDO
            ENDIF
          ENDDO D8985
        ENDDO D8980
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
!     DURING TILLAGE
!
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        DO NE=1,NumOfPlantChemElmnts
          D6400: DO M=1,jsken
            LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ) &
              +(XHVST1*CFOPE(NE,instruct,M,NZ)*NonstructalChemElmnts_pft(NE,NZ))*FWOODE(NE,k_woody_litr)

            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ) &
              +(XHVST1*CFOPE(NE,instruct,M,NZ)*NonstructalChemElmnts_pft(NE,NZ))*FWOODE(NE,k_fine_litr)
          ENDDO D6400
          NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)*XHVST
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByHarvest(I,J,NZ,CPOOLK)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(inout) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: L,K,M,NR,N,NB,NBX,NE
  real(r8):: ZPOOLG,ZPOLNG,ZPOOLX
  real(r8) :: ZPOLNX,XHVST(NumOfPlantChemElmnts)
  real(r8) :: XHVST1(NumOfPlantChemElmnts)
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
  real(r8) :: FFIRE(NumOfPlantChemElmnts)
  real(r8) :: FHVSE(NumOfPlantChemElmnts)
  real(r8) :: HTSTKX
  real(r8) :: PPOOLG
  real(r8) :: PPOLNG,PPOOLX,PPOLNX
  real(r8) :: VOLWPX
  real(r8) :: WHVSBL
  real(r8) :: WTSTKT
  real(r8) :: WTLSBX
  real(r8) :: WHVSTT,WHVSLF,WHVHSH,WHVEAH,WHVGRH,WHVSCP
  real(r8) :: WHVSTH,WHVRVH,WHVSLX,WHVSLY,WHVSCL,WHVSNL,WHVXXX
  real(r8) :: WHVSSX,WTSHTT,WHVSHX,WHVSHY,WHVSHH,WHVSCS,WHVSNS
  real(r8) :: WHVHSX,WHVHSY,WHVEAX,WHVEAY,WHVGRX,WHVGRY,WHVSNP
  real(r8) :: WHVSKX,WHVSTX,WHVSTY,WHVRVX,WHVRVY,WTNDG,WTNDNG
  real(r8) :: WTNDPG,WGLFGX,WGSHGX,WGLFGY,WGSHGY
  real(r8) :: WHVSBS,WHVSCX,WHVSNX,WVPLT
  real(r8) :: FHVSH1,FHVSHT
  real(r8) :: WGLFGE(NumOfPlantChemElmnts)
  integer :: NTG
!     begin_execution
  associate(                             &
    HVST     =>  plt_distb%HVST    , &
    EHVST    =>  plt_distb%EHVST   , &
    DCORP    =>  plt_distb%DCORP   , &
    THIN_pft     =>  plt_distb%THIN_pft    , &
    ITILL    =>  plt_distb%ITILL   , &
    IHVST    =>  plt_distb%IHVST   , &
    JHVST    =>  plt_distb%JHVST   , &
    VPO4F    =>  plt_distb%VPO4F   , &
    VN2OF    =>  plt_distb%VN2OF   , &
    VNH3F    =>  plt_distb%VNH3F   , &
    VOXYF    =>  plt_distb%VOXYF   , &
    CH4ByFire_pft    =>  plt_distb%CH4ByFire_pft   , &
    CO2ByFire_pft    =>  plt_distb%CO2ByFire_pft   , &
    UVOLO    =>  plt_ew%UVOLO      , &
    CanWatP    =>  plt_ew%CanWatP      , &
    PSICanP    =>  plt_ew%PSICanP      , &
    HCOB     =>  plt_photo%HCOB    , &
    CO2B     =>  plt_photo%CO2B    , &
    CPOOL3   =>  plt_photo%CPOOL3  , &
    CPOOL4   =>  plt_photo%CPOOL4  , &
    pftPlantPopulation       =>  plt_site%pftPlantPopulation       , &
    PPI      =>  plt_site%PPI      , &
    PPX      =>  plt_site%PPX      , &
    NU       =>  plt_site%NU       , &
    NJ       => plt_site%NJ        , &
    ZNOON    => plt_site%ZNOON     , &
    ZEROS    => plt_site%ZEROS     , &
    ZERO     => plt_site%ZERO      , &
    AREA3    => plt_site%AREA3     , &
    VOLWOU   => plt_site%VOLWOU    , &
     RootMycoNonstructElmnt_vr   => plt_biom% RootMycoNonstructElmnt_vr    , &
    WSRTL    => plt_biom%WSRTL     , &
     PopuPlantRootC_vr    => plt_biom% PopuPlantRootC_vr     , &
    RootStructBiomC_vr   => plt_biom%RootStructBiomC_vr    , &
    NonstructalChemElmnts_pft    => plt_biom%NonstructalChemElmnts_pft     , &
    CanPStalkC    => plt_biom%CanPStalkC     , &
    CanopyLeafShethC_pft     => plt_biom%CanopyLeafShethC_pft      , &
    StalkBiomassC_brch   => plt_biom%StalkBiomassC_brch    , &
    WTNDBE   => plt_biom%WTNDBE    , &
    CanopyReserveChemElmnts   => plt_biom%CanopyReserveChemElmnts    , &
    GrainChemElmnts_brch   => plt_biom%GrainChemElmnts_brch    , &
    StalkChemElmnts_brch  => plt_biom%StalkChemElmnts_brch   , &
    WTSHTBE  => plt_biom%WTSHTBE   , &
    HuskChemElmnts_brch  => plt_biom%HuskChemElmnts_brch   , &
    EarChemElmnts_brch  => plt_biom%EarChemElmnts_brch   , &
    InternodeChemElmnt_brch   => plt_biom%InternodeChemElmnt_brch    , &
    LeafPetioleBiomassC_brch  => plt_biom%LeafPetioleBiomassC_brch     , &
    ReserveChemElmnts_brch  => plt_biom%ReserveChemElmnts_brch   , &
    WTSTXBE  => plt_biom%WTSTXBE   , &
    NoduleNonstructElmnt_brch   => plt_biom%NoduleNonstructElmnt_brch     , &
    WGLFLE   => plt_biom%WGLFLE    , &
    PetioleChemElmnts_brch => plt_biom%PetioleChemElmnts_brch  , &
    NonstructElmnt_brch   => plt_biom%NonstructElmnt_brch    , &
    PetioleElmntNode_brch   => plt_biom%PetioleElmntNode_brch    , &
    LeafProteinCNode_brch     => plt_biom%LeafProteinCNode_brch      , &
    LeafChemElmntNode_brch    => plt_biom%LeafChemElmntNode_brch     , &
    LeafChemElmnts_brch  => plt_biom%LeafChemElmnts_brch   , &
    PetioleProteinCNode_brch   => plt_biom%PetioleProteinCNode_brch    , &
    StalkChemElmnts   => plt_biom%StalkChemElmnts    , &
    CanopyNonstructElementConc_pft   => plt_biom%CanopyNonstructElementConc_pft    , &
    NoduleNonstructCconc_pft   => plt_biom%NoduleNonstructCconc_pft    , &
    LeafChemElmnts    => plt_biom%LeafChemElmnts     , &
    GrainChemElmnts    => plt_biom%GrainChemElmnts     , &
    ShootChemElmnts_pft   => plt_biom%ShootChemElmnts_pft     , &
    HuskChemElmnts   => plt_biom%HuskChemElmnts    , &
    EarChemElmnts   => plt_biom%EarChemElmnts    , &
    SheathChemElmnts   => plt_biom%SheathChemElmnts    , &
    WTSHTA   => plt_biom%WTSHTA    , &
    WTRT1E   => plt_biom%WTRT1E    , &
    RTWT1E   => plt_biom%RTWT1E    , &
    WTRT2E   => plt_biom%WTRT2E    , &
    WTNDLE   => plt_biom%WTNDLE    , &
    RootNoduleNonstructElmnt_vr  => plt_biom%RootNoduleNonstructElmnt_vr   , &
    ZEROP    => plt_biom%ZEROP     , &
    ZEROL    => plt_biom%ZEROL     , &
    CanopyLeafCpft_lyr    => plt_biom%CanopyLeafCpft_lyr     , &
    FVRN     => plt_allom%FVRN     , &
    FWODRE   => plt_allom%FWODRE   , &
    FWOODE   => plt_allom%FWOODE   , &
    FWODBE   => plt_allom%FWODBE   , &
    FWODLE   => plt_allom%FWODLE   , &
    GRWTB    => plt_allom%GRWTB    , &
    iPlantBranchState    =>  plt_pheno%iPlantBranchState   , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP    , &
    iPlantCalendar    =>  plt_pheno%iPlantCalendar   , &
    MatureGroup_brch   =>  plt_pheno%MatureGroup_brch  , &
    VSTGX    =>  plt_pheno%VSTGX   , &
    iPlantMorphologyType   =>  plt_pheno%iPlantMorphologyType  , &
    iPlantTurnoverPattern   =>  plt_pheno%iPlantTurnoverPattern  , &
    doInitLeafOut    =>  plt_pheno%doInitLeafOut   , &
    iPlantPhenologyPattern   =>  plt_pheno%iPlantPhenologyPattern  , &
    HourThreshold4LeafOff    =>  plt_pheno%HourThreshold4LeafOff   , &
    Hours4LeafOff     =>  plt_pheno%Hours4LeafOff    , &
    iPlantPhenologyType   =>  plt_pheno%iPlantPhenologyType  , &
    TGSTGI   =>  plt_pheno%TGSTGI  , &
    TGSTGF   =>  plt_pheno%TGSTGF  , &
    FLG4     =>  plt_pheno%FLG4    , &
    MatureGroup_pft  =>  plt_pheno%MatureGroup_pft , &
    CORGC    =>  plt_soilchem%CORGC, &
    THETW    =>  plt_soilchem%THETW, &
    CFOPE    =>  plt_soilchem%CFOPE, &
    instruct =>  pltpar%instruct  , &
    ifoliar  =>  pltpar%ifoliar   , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    infoliar =>  pltpar%infoliar  , &
    k_fine_litr=> pltpar%k_fine_litr ,&
    k_woody_litr=> pltpar%k_woody_litr,&
    icwood   =>  pltpar%icwood    , &
    trcg_rootml     =>  plt_rbgc%trcg_rootml , &
    trcs_rootml     =>  plt_rbgc%trcs_rootml , &
    RootGasLoss_disturb    =>   plt_bgcr%RootGasLoss_disturb    , &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr     , &
    Eco_NBP_col     =>  plt_bgcr%Eco_NBP_col     , &
    CO2NetFix_pft     =>  plt_bgcr%CO2NetFix_pft     , &
    RCO2A    =>  plt_rbgc%RCO2A    , &
    RCO2M    =>  plt_rbgc%RCO2M    , &
    RCO2N    =>  plt_rbgc%RCO2N    , &
    SecndRootXNumL     =>  plt_morph%SecndRootXNumL    , &
    RTN2     =>  plt_morph%RTN2    , &
    RootLenPerP    =>  plt_morph%RootLenPerP   , &
    RTARP    =>  plt_morph%RTARP   , &
    RTVLP    =>  plt_morph%RTVLP   , &
    NGTopRootLayer      =>  plt_morph%NGTopRootLayer     , &
    MY       =>  plt_morph%MY      , &
    CanopyHeight       =>  plt_morph%CanopyHeight      , &
    RTVLW    =>  plt_morph%RTVLW   , &
    RootLenDensNLP    =>  plt_morph%RootLenDensNLP   , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType  , &
    CanopyLAgrid_lyr    =>  plt_morph%CanopyLAgrid_lyr   , &
    CanopyHeightz       =>  plt_morph%CanopyHeightz      , &
    LeafAreaLive_brch    =>  plt_morph%LeafAreaLive_brch   , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft     , &
    CanopyStemA_pft    =>  plt_morph%CanopyStemA_pft   , &
    InternodeHeightDying_brch   =>  plt_morph%InternodeHeightDying_brch  , &
    PrimRootXNumL    =>  plt_morph%PrimRootXNumL   , &
    SecndRootLen    =>  plt_morph%SecndRootLen   , &
    PrimRootLen    =>  plt_morph%PrimRootLen   , &
    InternodeHeightLive_brch   =>  plt_morph%InternodeHeightLive_brch  , &
    GRNXB    =>  plt_morph%GRNXB   , &
    GRNOB    =>  plt_morph%GRNOB   , &
    PetioleLengthNode_brch    =>  plt_morph%PetioleLengthNode_brch   , &
    LeafAreaNode_brch    =>  plt_morph%LeafAreaNode_brch   , &
    CanopyLeafApft_lyr    =>  plt_morph%CanopyLeafApft_lyr   , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr   , &
    CanPLNBLA    =>  plt_morph%CanPLNBLA   , &
    CanopyBranchStemApft_lyr    =>  plt_morph%CanopyBranchStemApft_lyr   , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft    , &
    NodeNumberAtAnthesis    =>  plt_morph%NodeNumberAtAnthesis   , &
    NumOfMainBranch_pft      =>  plt_morph%NumOfMainBranch_pft     , &
    NodeNumberToInitFloral    =>  plt_morph%NodeNumberToInitFloral   , &
    ShootNodeNumber    =>  plt_morph%ShootNodeNumber   , &
    ClumpFactor      =>  plt_morph%ClumpFactor     , &
    CanopyLA_grd    =>  plt_morph%CanopyLA_grd   , &
    iPlantPhotosynthesisType   =>  plt_photo%iPlantPhotosynthesisType    &
  )
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
  IF((IHVST(NZ).GE.0.AND.J.EQ.INT(ZNOON) &
    .AND.IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6) &
    .OR.(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6))THEN
!
!     ACCUMULATE ALL HARVESTED MATERIAL ABOVE CUTTING HEIGHT
!     ACCOUNTING FOR HARVEST EFFICIENCY ENTERED IN 'READQ'
!
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     PPX,PP=PFT population per m2,grid cell
!     THIN_pft=thinning:fraction of population removed
!     CF=clumping factor
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     CanopyLA_grd,CanopyLAgrid_lyr=leaf area of combined canopy, canopy layer
!     ARLFR,ARLFY=leaf area harvested,remaining
!     ZL=height to bottom of each canopy layer
!
    IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
      IF(JHVST(NZ).NE.ihv_tmareseed)THEN
        !terminate and reseed
        PPX(NZ)=PPX(NZ)*(1._r8-THIN_pft(NZ))
        pftPlantPopulation(NZ)=pftPlantPopulation(NZ)*(1._r8-THIN_pft(NZ))
      ELSE
!     PPI(NZ)=AMAX1(1.0_r8,0.5_r8*(PPI(NZ)+GRNO(NZ)/AREA3(NU)))
        PPX(NZ)=PPI(NZ)
        pftPlantPopulation(NZ)=PPX(NZ)*AREA3(NU)
      ENDIF
      IF(IHVST(NZ).EQ.3)THEN
        ClumpFactor(NZ)=ClumpFactor(NZ)*HVST(NZ)
      ENDIF
      IF(IHVST(NZ).LE.2.AND.HVST(NZ).LT.0.0)THEN
        ARLFY=(1._r8-ABS(HVST(NZ)))*CanopyLA_grd
        ARLFR=0._r8
        D9875: DO L=1,NumOfCanopyLayers1
          IF(CanopyHeightz(L).GT.CanopyHeightz(L-1).AND.CanopyLAgrid_lyr(L).GT.ZEROS.AND.ARLFR.LT.ARLFY)THEN
            IF(ARLFR+CanopyLAgrid_lyr(L).GT.ARLFY)THEN
              HVST(NZ)=CanopyHeightz(L-1)+((ARLFY-ARLFR)/CanopyLAgrid_lyr(L))*(CanopyHeightz(L)-CanopyHeightz(L-1))
            ENDIF
          ELSE
            HVST(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+CanopyLAgrid_lyr(L)
        ENDDO D9875
      ENDIF
      WHVSTT=0._r8
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
!     WTSHTA=average biomass in landscape grazing section
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WHVSTT=total phytomass grazed, removed
!     fTgrowCanP=temperature function for canopy growth
!     CCPOLP=nonstructural C concentration in canopy
!     NoduleNonstructCconc_pft=nonstructural C concentration in canopy nodules
!
      IF(WTSHTA(NZ).GT.ZEROP(NZ))THEN
        WHVSTT=HVST(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8 &
          *AREA3(NU)*ShootChemElmnts_pft(ielmc,NZ)/WTSHTA(NZ)
      ELSE
        WHVSTT=0._r8
      ENDIF
      IF(IHVST(NZ).EQ.6)THEN
        WHVSTT=WHVSTT*fTgrowCanP(NZ)
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
      WHVSLX=WHVSTT*EHVST(1,ipld_leaf,NZ)
      WHVSLY=AMIN1(LeafChemElmnts(ielmc,NZ),WHVSLX)
      WHVSLF=WHVSLY*(1._r8-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AZMAX1(WHVSLX-WHVSLY)
      WHVSSX=WHVSTT*EHVST(1,ipld_nofoliar,NZ)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     WHVXXX=grazing requirement unmet by non-foliar removal
!
      WTSHTT=SheathChemElmnts(ielmc,NZ)+HuskChemElmnts(ielmc,NZ)+EarChemElmnts(ielmc,NZ)+GrainChemElmnts(ielmc,NZ)
      IF(WTSHTT.GT.ZEROP(NZ))THEN
        WHVSHX=WHVSSX*SheathChemElmnts(ielmc,NZ)/WTSHTT+WHVXXX
        WHVSHY=AMIN1(SheathChemElmnts(ielmc,NZ),WHVSHX)
        WHVSHH=WHVSHY*(1._r8-CCPOLX)
        WHVSCS=WHVSHY*CCPOLX
        WHVSNS=WHVSHY*CCPLNX
        WHVXXX=AZMAX1(WHVSHX-WHVSHY)
        WHVHSX=WHVSSX*HuskChemElmnts(ielmc,NZ)/WTSHTT+WHVXXX
        WHVHSY=AMIN1(HuskChemElmnts(ielmc,NZ),WHVHSX)
        WHVHSH=WHVHSY
        WHVXXX=AZMAX1(WHVHSX-WHVHSY)
        WHVEAX=WHVSSX*EarChemElmnts(ielmc,NZ)/WTSHTT+WHVXXX
        WHVEAY=AMIN1(EarChemElmnts(ielmc,NZ),WHVEAX)
        WHVEAH=WHVEAY
        WHVXXX=AZMAX1(WHVEAX-WHVEAY)
        WHVGRX=WHVSSX*GrainChemElmnts(ielmc,NZ)/WTSHTT+WHVXXX
        WHVGRY=AMIN1(GrainChemElmnts(ielmc,NZ),WHVGRX)
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
      WHVSKX=WHVSTT*EHVST(1,ipld_woody,NZ)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     WHVXXX=grazing requirement unmet by stalk,reserve
!
      WTSTKT=StalkChemElmnts(ielmc,NZ)+CanopyReserveChemElmnts(ielmc,NZ)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
        WHVSTX=WHVSKX*StalkChemElmnts(ielmc,NZ)/WTSTKT+WHVXXX
        WHVSTY=AMIN1(StalkChemElmnts(ielmc,NZ),WHVSTX)
        WHVSTH=WHVSTY
        WHVXXX=AZMAX1(WHVSTX-WHVSTY)
        WHVRVX=WHVSKX*CanopyReserveChemElmnts(ielmc,NZ)/WTSTKT+WHVXXX
        WHVRVY=AMIN1(CanopyReserveChemElmnts(ielmc,NZ),WHVRVX)
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
          WHVSLY=AMIN1(LeafChemElmnts(ielmc,NZ)-WHVSLF-WHVSCL,WHVXXX)
          WHVSLF=WHVSLF+WHVSLY*(1._r8-CCPOLX)
          WHVSCL=WHVSCL+WHVSLY*CCPOLX
          WHVSNL=WHVSNL+WHVSLY*CCPLNX
          WHVXXX=AZMAX1(WHVXXX-WHVSLY)
          IF(WTSHTT.GT.ZEROP(NZ))THEN
            WHVSHX=WHVXXX*SheathChemElmnts(ielmc,NZ)/WTSHTT
            WHVSHY=AMIN1(SheathChemElmnts(ielmc,NZ),WHVSHX)
            WHVSHH=WHVSHH+WHVSHY*(1._r8-CCPOLX)
            WHVSCS=WHVSCS+WHVSHY*CCPOLX
            WHVSNS=WHVSNS+WHVSHY*CCPLNX
            WHVXXX=AZMAX1(WHVXXX-WHVSHY)
            WHVHSX=WHVXXX*HuskChemElmnts(ielmc,NZ)/WTSHTT
            WHVHSY=AMIN1(HuskChemElmnts(ielmc,NZ),WHVHSX)
            WHVHSH=WHVHSH+WHVHSY
            WHVXXX=AZMAX1(WHVXXX-WHVHSY)
            WHVEAX=WHVXXX*EarChemElmnts(ielmc,NZ)/WTSHTT
            WHVEAY=AMIN1(EarChemElmnts(ielmc,NZ),WHVEAX)
            WHVEAH=WHVEAH+WHVEAY
            WHVXXX=AZMAX1(WHVEAX-WHVEAY)
            WHVGRX=WHVXXX*GrainChemElmnts(ielmc,NZ)/WTSHTT
            WHVGRY=AMIN1(GrainChemElmnts(ielmc,NZ),WHVGRX)
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
            WGLFBL(L,NB,NZ)=WGLFBL(L,NB,NZ)+WGLFLE(ielmc,L,K,NB,NZ)
          enddo
        enddo
      ENDDO D9870
    ENDIF
!
!     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     ZL=height to bottom of each canopy layer
!     FHGT=fraction of canopy layer height not harvested
!     FHVSE(ielmc)=fraction of canopy layer mass not harvested
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
    D9865: DO L=NumOfCanopyLayers1,1,-1
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        IF(IHVST(NZ).NE.3)THEN
          IF(CanopyHeightz(L).GT.CanopyHeightz(L-1))THEN
            FHGT=AZMAX1(AMIN1(1.0_r8,1._r8-((CanopyHeightz(L))-HVST(NZ))/(CanopyHeightz(L)-CanopyHeightz(L-1))))
          ELSE
            FHGT=1.0_r8
          ENDIF
        ELSE
          FHGT=0._r8
        ENDIF
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FHVSE(ielmc)=AZMAX1(1._r8-(1._r8-FHGT)*EHVST(1,ipld_leaf,NZ))
          FHVSH=FHVSE(ielmc)
        ELSE
          FHVSE(ielmc)=AZMAX1(1._r8-THIN_pft(NZ))
          IF(IHVST(NZ).EQ.0)THEN
            FHVSH=1.0_r8-(1._r8-FHGT)*EHVST(1,ipld_leaf,NZ)*THIN_pft(NZ)
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
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTLF=PFT leaf C mass
!     WGLFBL=branch leaf C mass in canopy layer
!     WHVBSL,WHVSLF=layer,total leaf C mass removed
!     WGLFL=leaf node C in canopy layer
!     FHVSE(ielmc)=fraction of leaf node mass not harvested
!
      D9855: DO NB=1,NumOfBranches_pft(NZ)
        IF((IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6) &
          .AND.LeafChemElmnts(ielmc,NZ).GT.ZEROL(NZ))THEN
          WHVSBL=WHVSLF*AZMAX1(WGLFBL(L,NB,NZ))/LeafChemElmnts(ielmc,NZ)
        ELSE
          WHVSBL=0._r8
        ENDIF
        D9845: DO K=MaxNodesPerBranch1,0,-1
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBL.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(WGLFLE(ielmc,L,K,NB,NZ).GT.WHVSBL)THEN
                FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,(WGLFLE(ielmc,L,K,NB,NZ)-WHVSBL)/WGLFLE(ielmc,L,K,NB,NZ)))
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
!     CanPLNBLA,CanopyBranchStemApft_lyr=leaf,stalk node area in canopy layer
!     WTHTH1E(ielmc),WTHTH1E(ielmn),WTHTH1E(ielmp)=harvested leaf C,N,P
!     WTHTX1E(ielmc),WTHTX1E(ielmn),WTHTX1E(ielmp)=harvested leaf C,N,P to litter
!     WTHTH3E(ielmc),WTHTH3E(ielmn),WTHTH3E(ielmp)=harvested woody C,N,P
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested woody C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!

            WHVSBL=WHVSBL-(1._r8-FHVSE(ielmc))*WGLFLE(ielmc,L,K,NB,NZ)
            FHVSH1=1._r8-FHVSH
            FHVSHT=FHVSH-FHVSE(ielmc)
            DO NE=1,NumOfPlantChemElmnts
              WTHTH1E(NE)=WTHTH1E(NE)+FHVSH1*WGLFLE(NE,L,K,NB,NZ)*FWODLE(NE,k_fine_litr)
              WTHTX1E(NE)=WTHTX1E(NE)+FHVSHT*WGLFLE(NE,L,K,NB,NZ)*FWODLE(NE,k_fine_litr)
              WTHTH3E(NE)=WTHTH3E(NE)+FHVSH1*WGLFLE(NE,L,K,NB,NZ)*FWODLE(NE,k_woody_litr)
              WTHTX3E(NE)=WTHTX3E(NE)+FHVSHT*WGLFLE(NE,L,K,NB,NZ)*FWODLE(NE,k_woody_litr)
              WGLFLE(NE,L,K,NB,NZ)=FHVSE(ielmc)*WGLFLE(NE,L,K,NB,NZ)
            ENDDO
!
!     REMAINING LEAF C,N,P AND AREA
!
            CanPLNBLA(L,K,NB,NZ)=FHVSE(ielmc)*CanPLNBLA(L,K,NB,NZ)
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
        WGLFGE(1:NumOfPlantChemElmnts)=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanPLNBLA,CanopyLeafApft_lyr=leaf node,total area in canopy layer
!
        D9815: DO L=1,NumOfCanopyLayers1
          ARLFG=ARLFG+CanPLNBLA(L,K,NB,NZ)
          DO NE=1,NumOfPlantChemElmnts
            WGLFGE(NE)=WGLFGE(NE)+WGLFLE(NE,L,K,NB,NZ)
          ENDDO
          CanopyLeafApft_lyr(L,NZ)=CanopyLeafApft_lyr(L,NZ)+CanPLNBLA(L,K,NB,NZ)
          CanopyLeafCpft_lyr(L,NZ)=CanopyLeafCpft_lyr(L,NZ)+WGLFLE(ielmc,L,K,NB,NZ)
        ENDDO D9815
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FHVSETK=fraction of internode layer mass not harvested
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(LeafChemElmntNode_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ).AND.EHVST(1,ipld_leaf,NZ).GT.0.0)THEN
            FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(1._r8-(1._r8-AZMAX1(WGLFGE(ielmc)) &
              /LeafChemElmntNode_brch(ielmc,K,NB,NZ))*EHVST(1,ipld_nofoliar,NZ)/EHVST(1,ipld_leaf,NZ))))
            FHVSHK(K)=FHVSETK(K)
        ELSE
          IF(isclose(THIN_pft(NZ),0._r8))THEN
            FHVSETK(K)=1.0_r8-EHVST(1,ipld_nofoliar,NZ)
            FHVSHK(K)=FHVSETK(K)
          ELSE
            FHVSETK(K)=1.0_r8-THIN_pft(NZ)
            IF(IHVST(NZ).EQ.0)THEN
              FHVSHK(K)=1.0_r8-EHVST(1,ipld_nofoliar,NZ)*THIN_pft(NZ)
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
      WGLFGY=WGLFGY+LeafChemElmntNode_brch(ielmc,K,NB,NZ)
      DO NE=1,NumOfPlantChemElmnts
        LeafChemElmnts_brch(NE,NB,NZ)=LeafChemElmnts_brch(NE,NB,NZ)-LeafChemElmntNode_brch(NE,K,NB,NZ)+WGLFGE(NE)
      ENDDO
      LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-LeafAreaNode_brch(K,NB,NZ)+ARLFG
      IF(LeafAreaNode_brch(K,NB,NZ).GT.ZEROP(NZ))THEN
        LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*ARLFG/LeafAreaNode_brch(K,NB,NZ)
      ELSE
        LeafProteinCNode_brch(K,NB,NZ)=0._r8
      ENDIF
      LeafAreaNode_brch(K,NB,NZ)=ARLFG
      DO NE=1,NumOfPlantChemElmnts
        LeafChemElmntNode_brch(NE,K,NB,NZ)=WGLFGE(NE)
      ENDDO
      WGLFGX=WGLFGX+LeafChemElmntNode_brch(ielmc,K,NB,NZ)
    ENDDO D9825
!
!     CUT SHEATHS OR PETIOLES AND STALKS HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSHE,WTSHEB=PFT,branch petiole C mass
!     WHVSBS,WHVSHH=branch, PFT petiole C mass removed
!     InternodeHeightLive_brch=internode length
!     HTSTKX=internode length removed
!
      HTSTKX=0._r8
      IF((IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6) &
        .AND.SheathChemElmnts(ielmc,NZ).GT.ZEROP(NZ))THEN
        WHVSBS=WHVSHH*PetioleChemElmnts_brch(ielmc,NB,NZ)/SheathChemElmnts(ielmc,NZ)
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
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WHVSBS=branch petiole C mass removed
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
!     FHVSETK=fraction of internode layer mass not harvested
!     WTHTH2E(ielmc),WTHTH2E(ielmn),WTHTH2E(ielmp)=harvested petiole C,N,P
!     WTHTX2E(ielmc),WTHTX2E(ielmn),WTHTX2E(ielmp)=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     PetioleLengthNode_brch,InternodeHeightLive_brch=petiole,internode length
!
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBS.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(PetioleElmntNode_brch(ielmc,K,NB,NZ).GT.WHVSBS)THEN
                FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(PetioleElmntNode_brch(ielmc,K,NB,NZ)-WHVSBS)/PetioleElmntNode_brch(ielmc,K,NB,NZ)))
                FHVSHK(K)=FHVSETK(K)
              ELSE
                FHVSETK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1._r8-FHVSETK(K))*PetioleElmntNode_brch(ielmc,K,NB,NZ)
            DO NE=1,NumOfPlantChemElmnts
              WTHTH2E(NE)=WTHTH2E(NE)+(1._r8-FHVSHK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)
              WTHTX2E(NE)=WTHTX2E(NE)+(FHVSHK(K)-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)

              WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSHK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
              WTHTX3E(NE)=WTHTX3E(NE)+(FHVSHK(K)-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
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

            DO NE=1,NumOfPlantChemElmnts
              PetioleChemElmnts_brch(NE,NB,NZ)=PetioleChemElmnts_brch(NE,NB,NZ) &
                -(1._r8-FHVSETK(K))*PetioleElmntNode_brch(NE,K,NB,NZ)
              PetioleElmntNode_brch(NE,K,NB,NZ)=FHVSETK(K)*PetioleElmntNode_brch(NE,K,NB,NZ)
            ENDDO
!            PetioleProteinCNode_brch(K,NB,NZ)=FHVSETK(K)*PetioleProteinCNode_brch(K,NB,NZ)
            IF(IHVST(NZ).LE.2.AND.PetioleLengthNode_brch(K,NB,NZ).GT.0.0_r8)THEN
              FHGT=AZMAX1(AMIN1(1.0_r8,(InternodeHeightLive_brch(K,NB,NZ) &
                +PetioleLengthNode_brch(K,NB,NZ)-HVST(NZ))/PetioleLengthNode_brch(K,NB,NZ)))
              PetioleLengthNode_brch(K,NB,NZ)=(1._r8-FHGT)*PetioleLengthNode_brch(K,NB,NZ)
            ELSE
              PetioleLengthNode_brch(K,NB,NZ)=FHVSETK(K)*PetioleLengthNode_brch(K,NB,NZ)
            ENDIF
            WGSHGX=WGSHGX+PetioleElmntNode_brch(ielmc,K,NB,NZ)
!     IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
!     IF(InternodeHeightLive_brch(K,NB,NZ).GT.HVST(NZ)
!    2.OR.IHVST(NZ).EQ.3)THEN
!     IF(isclose(FHVSETK(K),0._r8).AND.K.GT.0)THEN
!     IF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ)))THEN
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
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     FHVSE(ielmc)=fraction of leaf+petiole node mass not harvested
!     CPOOLG,ZPOOLG,PPOOLG=branch non-structural C,N,P mass after harvest
!     CPOLNG,ZPOLNG,PPOLNG=nonstructural C,N,P in bacteria after harvest
!     WTNDG,WTNDNG,WTNDPG=bacterial C,N,P mass after harvest
!     WTLS,LeafPetioleBiomassC_brch=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
        CPOOLX=AZMAX1(NonstructElmnt_brch(ielmc,NB,NZ))
        ZPOOLX=AZMAX1(NonstructElmnt_brch(ielmn,NB,NZ))
        PPOOLX=AZMAX1(NonstructElmnt_brch(ielmp,NB,NZ))
        CPOLNX=AZMAX1(NoduleNonstructElmnt_brch(ielmc,NB,NZ))
        ZPOLNX=AZMAX1(NoduleNonstructElmnt_brch(ielmn,NB,NZ))
        PPOLNX=AZMAX1(NoduleNonstructElmnt_brch(ielmp,NB,NZ))
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(WGLFGY+WGSHGY.GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
            CPOOLG=CPOOLX*FHVSE(ielmc)
            ZPOOLG=ZPOOLX*FHVSE(ielmc)
            PPOOLG=PPOOLX*FHVSE(ielmc)
            CPOLNG=CPOLNX*FHVSE(ielmc)
            ZPOLNG=ZPOLNX*FHVSE(ielmc)
            PPOLNG=PPOLNX*FHVSE(ielmc)
            WTNDG=WTNDBE(ielmc,NB,NZ)*FHVSE(ielmc)
            WTNDNG=WTNDBE(ielmn,NB,NZ)*FHVSE(ielmc)
            WTNDPG=WTNDBE(ielmp,NB,NZ)*FHVSE(ielmc)
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
            WTLSBX=AZMAX1(LeafPetioleBiomassC_brch(NB,NZ))
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
              WTNDG=WTNDBE(ielmc,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDNG=WTNDBE(ielmn,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDPG=WTNDBE(ielmp,NB,NZ)*(1._r8-WHVSNX/CPOLNX)
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
!     WTHTH0E(ielmc),WTHTH0E(ielmn),WTHTH0E(ielmp)=nonstructural C,N,P removed
!
        WTHTH0E(ielmc)=WTHTH0E(ielmc)+CPOOLX-CPOOLG+CPOLNX-CPOLNG
        WTHTH0E(ielmn)=WTHTH0E(ielmn)+ZPOOLX-ZPOOLG+ZPOLNX-ZPOLNG
        WTHTH0E(ielmp)=WTHTH0E(ielmp)+PPOOLX-PPOOLG+PPOLNX-PPOLNG
        WTHTH0E(ielmc)=WTHTH0E(ielmc)+WTNDBE(ielmc,NB,NZ)-WTNDG
        WTHTH0E(ielmn)=WTHTH0E(ielmn)+WTNDBE(ielmn,NB,NZ)-WTNDNG
        WTHTH0E(ielmp)=WTHTH0E(ielmp)+WTNDBE(ielmp,NB,NZ)-WTNDPG
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
        WTNDBE(ielmc,NB,NZ)=WTNDG
        WTNDBE(ielmn,NB,NZ)=WTNDNG
        WTNDBE(ielmp,NB,NZ)=WTNDPG
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CPOOLG=branch non-structural C mass after harvest
!     WTHTH0E(ielmc),WTHTH0E(ielmn),WTHTH0E(ielmp)=nonstructural C,N,P removed
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
        IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo.AND.CPOOLX.GT.ZEROP(NZ))THEN
          FHVST4=CPOOLG/CPOOLX
          D9810: DO K=1,MaxNodesPerBranch1
            WTHTH0E(ielmc)=WTHTH0E(ielmc)+(1._r8-FHVST4)*CPOOL3(K,NB,NZ)
            WTHTH0E(ielmc)=WTHTH0E(ielmc)+(1._r8-FHVST4)*CPOOL4(K,NB,NZ)
            WTHTH0E(ielmc)=WTHTH0E(ielmc)+(1._r8-FHVST4)*CO2B(K,NB,NZ)
            WTHTH0E(ielmc)=WTHTH0E(ielmc)+(1._r8-FHVST4)*HCOB(K,NB,NZ)
            CPOOL3(K,NB,NZ)=FHVST4*CPOOL3(K,NB,NZ)
            CPOOL4(K,NB,NZ)=FHVST4*CPOOL4(K,NB,NZ)
            CO2B(K,NB,NZ)=FHVST4*CO2B(K,NB,NZ)
            HCOB(K,NB,NZ)=FHVST4*HCOB(K,NB,NZ)
          ENDDO D9810
        ENDIF
!
!     CUT STALKS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTSTKX=internode length removed
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     FHGT=fraction of canopy layer height not harvested
!     FHVSE(ielmc)=fraction of canopy layer mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WTSTK=stalk C mass
!
!
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(HTSTKX.GT.ZERO)THEN
            IF(IHVST(NZ).NE.3)THEN
              FHGT=AZMAX1(AMIN1(1.0_r8,HVST(NZ)/HTSTKX))
            ELSE
              FHGT=0._r8
            ENDIF
            IF(isclose(THIN_pft(NZ),0._r8))THEN
              FHVSE(ielmc)=AZMAX1(1._r8-(1._r8-FHGT)*EHVST(1,ipld_woody,NZ))
              FHVSH=FHVSE(ielmc)
            ELSE
              FHVSE(ielmc)=AZMAX1(1._r8-THIN_pft(NZ))
              IF(IHVST(NZ).EQ.0)THEN
                FHVSH=1.0_r8-(1._r8-FHGT)*EHVST(1,ipld_woody,NZ)*THIN_pft(NZ)
              ELSE
                FHVSH=FHVSE(ielmc)
              ENDIF
            ENDIF
          ELSE
            FHVSE(ielmc)=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ELSE
          IF(StalkChemElmnts(ielmc,NZ).GT.ZEROL(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/StalkChemElmnts(ielmc,NZ)))
            FHVSH=FHVSE(ielmc)
          ELSE
            FHVSE(ielmc)=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ENDIF
    !
!     HARVESTED STALK C,N,P
!
!     WTHTH3E(ielmc),WTHTH3E(ielmn),WTHTH3E(ielmp)=harvested stalk C,N,P
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested stalk C,N,P to litter
!     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in harvested stalk
!
        DO NE=1,NumOfPlantChemElmnts
          WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSH)*StalkChemElmnts_brch(NE,NB,NZ)
          WTHTX3E(NE)=WTHTX3E(NE)+(FHVSH-FHVSE(ielmc))*StalkChemElmnts_brch(NE,NB,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
          StalkChemElmnts_brch(NE,NB,NZ)=FHVSE(ielmc)*StalkChemElmnts_brch(NE,NB,NZ)
          WTSTXBE(NE,NB,NZ)=FHVSE(ielmc)*WTSTXBE(NE,NB,NZ)
        ENDDO

        StalkBiomassC_brch(NB,NZ)=FHVSE(ielmc)*StalkBiomassC_brch(NB,NZ)
!
!     CUT STALK NODES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     InternodeHeightDying_brch,InternodeHeightLive_brch=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTSTK=stalk C mass
!     InternodeChemElmnt_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
        D9820: DO K=MaxNodesPerBranch1,0,-1
          IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
            IF(InternodeHeightDying_brch(K,NB,NZ).GT.ZERO)THEN
              IF(IHVST(NZ).NE.3)THEN
                FHGTK=AZMAX1(AMIN1(1.0_r8,(InternodeHeightLive_brch(K,NB,NZ)-HVST(NZ))/InternodeHeightDying_brch(K,NB,NZ)))
              ELSE
                FHGTK=0._r8
              ENDIF
              IF(isclose(THIN_pft(NZ),0._r8))THEN
                FHVSETS=AZMAX1(1._r8-FHGTK*EHVST(1,ipld_woody,NZ))
              ELSE
                FHVSETS=AZMAX1(1._r8-THIN_pft(NZ))
              ENDIF
            ELSE
              FHVSETS=1.0_r8
            ENDIF
          ELSE
            IF(StalkChemElmnts(ielmc,NZ).GT.ZEROP(NZ))THEN
              FHVSETS=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/StalkChemElmnts(ielmc,NZ)))
            ELSE
              FHVSETS=1.0_r8
            ENDIF
          ENDIF
          DO NE=1,NumOfPlantChemElmnts
            InternodeChemElmnt_brch(NE,K,NB,NZ)=FHVSETS*InternodeChemElmnt_brch(NE,K,NB,NZ)
          ENDDO
          IF(IHVST(NZ).LE.2.AND.isclose(THIN_pft(NZ),0._r8))THEN
            InternodeHeightDying_brch(K,NB,NZ)=FHVSETS*InternodeHeightDying_brch(K,NB,NZ)
            InternodeHeightLive_brch(K,NB,NZ)=AMIN1(InternodeHeightLive_brch(K,NB,NZ),HVST(NZ))
          ENDIF

        ENDDO D9820
!
!     CUT STALK RESERVES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSTKB=C mass remaining in harvested stalk
!     WTRSV=stalk reserve C mass
!     WHVRVH=remaining stalk reserve C mass
!     FHVSE(ielmc)=fraction of reserve mass not harvested
!
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(StalkChemElmnts_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=FHVSE(ielmc)
            FHVSH=FHVSH
          ELSE
            FHVSE(ielmc)=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(CanopyReserveChemElmnts(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,1._r8-WHVRVH/CanopyReserveChemElmnts(ielmc,NZ)))
            FHVSH=FHVSE(ielmc)
          ELSE
            FHVSE(ielmc)=0._r8
            FHVSH=0._r8
          ENDIF
        ENDIF
!
!     HARVESTED STALK RESERVE C,N,P
!
!     WTHTH3E(ielmc),WTHTH3E(ielmn),WTHTH3E(ielmp)=harvested stalk C,N,P
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested stalk C,N,P to litter
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!
        DO NE=1,NumOfPlantChemElmnts
          WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSH)*ReserveChemElmnts_brch(NE,NB,NZ)
          WTHTX3E(NE)=WTHTX3E(ielmc)+(FHVSH-FHVSE(ielmc))*ReserveChemElmnts_brch(NE,NB,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
          ReserveChemElmnts_brch(NE,NB,NZ)=FHVSE(ielmc)*ReserveChemElmnts_brch(NE,NB,NZ)
        ENDDO
!
!     CUT REPRODUCTIVE ORGANS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FHVSETG,FHVSETH,FHVSETE=fraction of grain,husk,ear mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(HVST(NZ).LT.HTSTKX.OR.IHVST(NZ).EQ.1.OR.IHVST(NZ).EQ.3)THEN
            IF(isclose(THIN_pft(NZ),0._r8))THEN
              FHVSETG=1.0_r8-EHVST(1,ipld_nofoliar,NZ)
              FHVSHG=FHVSETG
            ELSE
              FHVSETG=1.0_r8-THIN_pft(NZ)
              FHVSHG=1.0_r8-EHVST(1,ipld_nofoliar,NZ)*THIN_pft(NZ)
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
          IF(HuskChemElmnts(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETH=AZMAX1(AMIN1(1.0_r8,1._r8-WHVHSH/HuskChemElmnts(ielmc,NZ)))
            FHVSHH=FHVSETH
          ELSE
            FHVSETH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(EarChemElmnts(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETE=AZMAX1(AMIN1(1.0_r8,1._r8-WHVEAH/EarChemElmnts(ielmc,NZ)))
            FHVSHE=FHVSETE
          ELSE
            FHVSETE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(GrainChemElmnts(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETG=AZMAX1(AMIN1(1.0_r8,1._r8-WHVGRH/GrainChemElmnts(ielmc,NZ)))
            FHVSHG=FHVSETG
          ELSE
            FHVSETG=1.0_r8
            FHVSHG=1.0_r8
          ENDIF
        ENDIF
!
!     HARVESTED REPRODUCTIVE C,N,P
!
!     WTHTH2E(ielmc),WTHTH2E(ielmn),WTHTH2E(ielmp)=reproductive C,N,P removed
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTHTGE(ielmc),WTHTGE(ielmn),WTHTGE(ielmp)=grain harvested
!
        DO NE=1,NumOfPlantChemElmnts
          WTHTH2E(NE)=WTHTH2E(NE)+(1._r8-FHVSHH)*HuskChemElmnts_brch(NE,NB,NZ)+(1._r8-FHVSHE) &
            *EarChemElmnts_brch(NE,NB,NZ)+(1._r8-FHVSHG)*GrainChemElmnts_brch(NE,NB,NZ)
          WTHTX2E(NE)=WTHTX2E(NE)+(FHVSHH-FHVSETH)*HuskChemElmnts_brch(NE,NB,NZ)+(FHVSHE-FHVSETE) &
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
        GRNXB(NB,NZ)=FHVSETG*GRNXB(NB,NZ)
        GRNOB(NB,NZ)=FHVSETG*GRNOB(NB,NZ)
        GRWTB(NB,NZ)=FHVSETG*GRWTB(NB,NZ)
!
!     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
!
!     CPOOLK=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WTLSB=leaf+petiole mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     StalkBiomassC_brch=stalk sapwood mass
!     PSICanP=canopy water potential
!     CanWatP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        CPOOLK(NB,NZ)=0._r8
        D1325: DO K=1,MaxNodesPerBranch1
          CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
            +CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
            +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
        ENDDO D1325

        LeafPetioleBiomassC_brch(NB,NZ)=AZMAX1(LeafChemElmnts_brch(ielmc,NB,NZ) &
          +PetioleChemElmnts_brch(ielmc,NB,NZ))

        WTSHTBE(ielmc,NB,NZ)=WTSHTBE(ielmc,NB,NZ)+CPOOLK(NB,NZ)
        DO NE=1,NumOfPlantChemElmnts
          WTSHTBE(NE,NB,NZ)=AZMAX1(WTSHTBE(ielmc,NB,NZ)+LeafChemElmnts_brch(NE,NB,NZ) &
            +PetioleChemElmnts_brch(NE,NB,NZ)+StalkChemElmnts_brch(NE,NB,NZ)+ReserveChemElmnts_brch(NE,NB,NZ) &
            +HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ)+GrainChemElmnts_brch(NE,NB,NZ) &
            +NonstructElmnt_brch(NE,NB,NZ))
        ENDDO
        VOLWPX=CanWatP(NZ)
        WVPLT=AZMAX1(CanopyLeafShethC_pft(NZ)+CanPStalkC(NZ))

        FDM=get_FDM(PSICanP(NZ))
!        APSILT=ABS(PSICanP(NZ))
!        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        CanWatP(NZ)=ppmc*WVPLT/FDM

        VOLWOU=VOLWOU+VOLWPX-CanWatP(NZ)
        UVOLO=UVOLO+VOLWPX-CanWatP(NZ)
!
!     RESET PHENOLOGY, GROWTH STAGE IF STALKS ARE CUT
!
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     ZC=canopy height
!     iPlantPhenologyType=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     iPlantCalendar(ipltcal_Emerge,=emergence date
!     GROUP=node number required for floral initiation
!     NodeNumberToInitFloral=node number at floral initiation
!     NodeNumberAtAnthesis=node number at flowering
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     doInitLeafOut=flag for initializing leafout
!
        IF((iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ)))) &
          .AND.(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6) &
          .AND.CanopyHeight(NZ).GT.HVST(NZ))THEN
          IF((iPlantPhenologyType(NZ).NE.0.AND.Hours4LeafOff(NB,NZ) &
            .LE.FVRN(iPlantPhenologyType(NZ))*HourThreshold4LeafOff(NB,NZ)) &
            .OR.(iPlantPhenologyType(NZ).EQ.0.AND.iPlantCalendar(ipltcal_Emerge,NB,NZ).NE.0))THEN
            MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
            NodeNumberToInitFloral(NB,NZ)=ShootNodeNumber(NB,NZ)
            NodeNumberAtAnthesis(NB,NZ)=0._r8
            VSTGX(NB,NZ)=0._r8
            TGSTGI(NB,NZ)=0._r8
            TGSTGF(NB,NZ)=0._r8
            FLG4(NB,NZ)=0._r8
            iPlantCalendar(ipltcal_Emerge,NB,NZ)=I
            D3005: DO M=2,10
              iPlantCalendar(M,NB,NZ)=0
            ENDDO D3005
            doInitLeafOut(NB,NZ)=0
            IF(NB.EQ.NumOfMainBranch_pft(NZ))THEN
              D3010: DO NBX=1,NumOfBranches_pft(NZ)
                IF(NBX.NE.NumOfMainBranch_pft(NZ))THEN
                  MatureGroup_brch(NBX,NZ)=MatureGroup_pft(NZ)
                  NodeNumberToInitFloral(NBX,NZ)=ShootNodeNumber(NBX,NZ)
                  NodeNumberAtAnthesis(NBX,NZ)=0._r8
                  VSTGX(NBX,NZ)=0._r8
                  TGSTGI(NBX,NZ)=0._r8
                  TGSTGF(NBX,NZ)=0._r8
                  FLG4(NBX,NZ)=0._r8
                  iPlantCalendar(ipltcal_Emerge,NBX,NZ)=I
                  D3015: DO M=2,10
                    iPlantCalendar(M,NBX,NZ)=0
                  ENDDO D3015
                  doInitLeafOut(NBX,NZ)=0
                ENDIF
              ENDDO D3010
            ENDIF
          ENDIF
        ENDIF
!
!     DEATH OF BRANCH IF KILLING HARVEST ENTERED IN 'READQ'
!
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!     PP=PFT population
!     WTLS=total PFT leaf+petiole C mass
!     WTSTK=total PFT stalk C mass
!     WVSTK=total PFT sapwood C mass
!     CanopyBranchStemApft_lyr=total PFT stalk surface area
!
        IF(JHVST(NZ).NE.ihv_noaction)then
          iPlantBranchState(NB,NZ)=iDead
        endif
        IF(pftPlantPopulation(NZ).LE.0.0)then
          iPlantBranchState(NB,NZ)=iDead
        endif
      ENDDO D9835
      CanopyLeafShethC_pft(NZ)=0._r8
      StalkChemElmnts(ielmc,NZ)=0._r8
      CanPStalkC(NZ)=0._r8
      CanopyStemA_pft(NZ)=0._r8
      D9840: DO NB=1,NumOfBranches_pft(NZ)
        CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetioleBiomassC_brch(NB,NZ)
        StalkChemElmnts(ielmc,NZ)=StalkChemElmnts(ielmc,NZ)+StalkChemElmnts_brch(ielmc,NB,NZ)
        CanPStalkC(NZ)=CanPStalkC(NZ)+StalkBiomassC_brch(NB,NZ)
        D9830: DO L=1,NumOfCanopyLayers1
          CanopyStemA_pft(NZ)=CanopyStemA_pft(NZ)+CanopyBranchStemApft_lyr(L,NB,NZ)
        ENDDO D9830
      ENDDO D9840
!
!     ROOT LITTERFALL FROM HARVESTING OR FIRE
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     THETW=soil water concentration
!     CORGC=SOC concentration
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     EFIRE=combustion  of N,P relative to C
!     FHVSE(ielmc),FHVSE(ielmn),FHVSE(ielmp)=fraction of root layer C,N,P not removed by disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CO2ByFire_pft,CH4ByFire_pft,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     Eco_NBP_col=total net biome productivity
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        XHVST(ielmc)=1.0_r8-THIN_pft(NZ)
        D3985: DO N=1,MY(NZ)
          D3980: DO L=NU,NJ
            IF(IHVST(NZ).NE.5)THEN
              XHVST(ielmc)=1.0_r8-THIN_pft(NZ)
              XHVST(ielmn)=XHVST(ielmc)
              XHVST(ielmp)=XHVST(ielmc)
              FFIRE(1:NumOfPlantChemElmnts)=0._r8
            ELSE
              IF(THETW(L).GT.FVLWB.OR.CORGC(L).LE.FORGC.OR.ITILL.NE.22)THEN
                XHVST(ielmc)=1.0_r8
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE(1:NumOfPlantChemElmnts)=0._r8
              ELSE
                XHVST(ielmc)=1.0_r8-DCORP*EHVST(1,ipld_woody,NZ) &
                  *AMIN1(1.0_r8,(CORGC(L)-FORGC)/(orgcden-FORGC))
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE(ielmc)=EHVST(2,ipld_woody,NZ)
                FFIRE(ielmn)=FFIRE(ielmc)*EFIRE(1,IHVST(NZ))
                FFIRE(ielmp)=FFIRE(ielmc)*EFIRE(2,IHVST(NZ))
              ENDIF
            ENDIF
            XHVST1=1._r8-XHVST
            D3385: DO M=1,jsken
              DO NE=1,NumOfPlantChemElmnts
                FHVSE(NE)=XHVST1(NE)*CFOPE(NE,instruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
                LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+(1._r8-FFIRE(NE))*FHVSE(NE)
              ENDDO
              CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
              CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
              VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667
              VNH3F(NZ)=VNH3F(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
              VN2OF(NZ)=VN2OF(NZ)-0.0_r8
              VPO4F(NZ)=VPO4F(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
              CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
              Eco_NBP_col=Eco_NBP_col-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
              DO NR=1,NumRootAxes_pft(NZ)
                DO NE=1,NumOfPlantChemElmnts
                  FHVSE(NE)=XHVST1(NE)*CFOPE(NE,icwood,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+(1._r8-FFIRE(NE))*FHVSE(NE)
                ENDDO
                CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667
                VNH3F(NZ)=VNH3F(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
                VN2OF(NZ)=VN2OF(NZ)-0.0
                VPO4F(NZ)=VPO4F(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
                CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                Eco_NBP_col=Eco_NBP_col-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)

                DO NE=1,NumOfPlantChemElmnts
                  FHVSE(NE)=XHVST1(NE)*CFOPE(NE,iroot,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ) &
                    +(1._r8-FFIRE(NE))*FHVSE(NE)
                ENDDO
                CO2ByFire_pft(NZ)=CO2ByFire_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                CH4ByFire_pft(NZ)=CH4ByFire_pft(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667_r8
                VNH3F(NZ)=VNH3F(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
                VN2OF(NZ)=VN2OF(NZ)-0.0_r8
                VPO4F(NZ)=VPO4F(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
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
              RootGasLoss_disturb(NTG,NZ)=RootGasLoss_disturb(NTG,NZ)-XHVST1(ielmc) &
                *(trcg_rootml(idg_CO2,N,L,NZ)+trcs_rootml(idg_CO2,N,L,NZ))
              trcg_rootml(NTG,N,L,NZ)=XHVST(ielmc)*trcg_rootml(NTG,N,L,NZ)
              trcs_rootml(NTG,N,L,NZ)=XHVST(ielmc)*trcs_rootml(NTG,N,L,NZ)
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
!     WSRTL=root protein C mass
!     RTN1,SecndRootXNumL=number of primary,secondary root axes
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D3960: DO NR=1,NumRootAxes_pft(NZ)
              DO NE=1,NumOfPlantChemElmnts
                WTRT1E(NE,N,L,NR,NZ)=WTRT1E(NE,N,L,NR,NZ)*XHVST(NE)
                WTRT2E(NE,N,L,NR,NZ)=WTRT2E(NE,N,L,NR,NZ)*XHVST(NE)
                RTWT1E(NE,N,NR,NZ)=RTWT1E(NE,N,NR,NZ)*XHVST(NE)
              ENDDO
              PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)*XHVST(ielmc)
              SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)*XHVST(ielmc)
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST(ielmc)
            ENDDO D3960
            DO NE=1,NumOfPlantChemElmnts
               RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)*XHVST(NE)
            ENDDO
            RootStructBiomC_vr(N,L,NZ)=RootStructBiomC_vr(N,L,NZ)*XHVST(ielmc)
             PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)*XHVST(ielmc)
            WSRTL(N,L,NZ)=WSRTL(N,L,NZ)*XHVST(ielmc)
            PrimRootXNumL(N,L,NZ)=PrimRootXNumL(N,L,NZ)*XHVST(ielmc)
            SecndRootXNumL(N,L,NZ)=SecndRootXNumL(N,L,NZ)*XHVST(ielmc)
            RootLenPerP(N,L,NZ)=RootLenPerP(N,L,NZ)*XHVST(ielmc)
            RootLenDensNLP(N,L,NZ)=RootLenDensNLP(N,L,NZ)*XHVST(ielmc)
            RTVLP(N,L,NZ)=RTVLP(N,L,NZ)*XHVST(ielmc)
            RTVLW(N,L,NZ)=RTVLW(N,L,NZ)*XHVST(ielmc)
            RTARP(N,L,NZ)=RTARP(N,L,NZ)*XHVST(ielmc)
            RCO2M(N,L,NZ)=RCO2M(N,L,NZ)*XHVST(ielmc)
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
            IF(iPlantNfixType(NZ).NE.0.AND.N.EQ.1)THEN
              DO NE=1,NumOfPlantChemElmnts
                D3395: DO M=1,jsken
                  LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+XHVST1(NE) &
                    *(CFOPE(NE,iroot,M,NZ)*WTNDLE(NE,L,NZ) &
                    +CFOPE(NE,instruct,M,NZ)*RootNoduleNonstructElmnt_vr(NE,L,NZ))
                ENDDO D3395
                WTNDLE(NE,L,NZ)=WTNDLE(NE,L,NZ)*XHVST(NE)
                RootNoduleNonstructElmnt_vr(NE,L,NZ)=RootNoduleNonstructElmnt_vr(NE,L,NZ)*XHVST(NE)
              ENDDO
            ENDIF
          ENDDO D3980
        ENDDO D3985
!
!     STORAGE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        IF(iPlantPhenologyPattern(NZ).NE.iplt_annual)THEN
          DO NE=1,NumOfPlantChemElmnts
            D3400: DO M=1,jsken
              LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ) &
                +(XHVST1(NE)*CFOPE(NE,instruct,M,NZ)*NonstructalChemElmnts_pft(NE,NZ))*FWOODE(NE,k_woody_litr)

              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ) &
                +(XHVST1(NE)*CFOPE(NE,instruct,M,NZ)*NonstructalChemElmnts_pft(NE,NZ))*FWOODE(NE,k_fine_litr)
            ENDDO D3400
            NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)*XHVST(NE)
          ENDDO
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

end module PlantDisturbsMod
