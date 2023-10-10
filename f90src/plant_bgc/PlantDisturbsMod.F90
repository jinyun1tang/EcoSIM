module PlantDisturbsMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  use PlantMathFuncMod, only : get_FDM
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section

! disturbance variables
  real(r8) :: WTHTH0E(NumOfPlantChemElements)
  real(r8) :: WTHTH1E(NumOfPlantChemElements)
  real(r8) :: WTHTH2E(NumOfPlantChemElements)
  real(r8) :: WTHTH3E(NumOfPlantChemElements)
  real(r8) :: WTHTH4E(NumOfPlantChemElements)
  real(r8) :: WTHTR1E(NumOfPlantChemElements)
  real(r8) :: WTHTR2E(NumOfPlantChemElements)
  real(r8) :: WTHTR3E(NumOfPlantChemElements)
  real(r8) :: WTHTR4E(NumOfPlantChemElements)
  real(r8) :: WTHTX0E(NumOfPlantChemElements)
  real(r8) :: WTHTX1E(NumOfPlantChemElements)
  real(r8) :: WTHTX2E(NumOfPlantChemElements)
  real(r8) :: WTHTX3E(NumOfPlantChemElements)
  real(r8) :: WTHTX4E(NumOfPlantChemElements)
  real(r8) :: WTHTGE(NumOfPlantChemElements)
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
    WTSTDE   =>  plt_biom%WTSTDE  , &
    WTSTGE   =>  plt_biom%WTSTGE    &
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
      IF(WTSTGE(ielmc,NZ).GT.ZEROP(NZ))THEN
        WHVSTD=HVST(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8*AREA3(NU)*EHVST(1,4,NZ)
        FHVSE(ielmc)=AZMAX1(1._r8-WHVSTD/WTSTGE(ielmc,NZ))
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
      WTHTH4E(ielmc)=WTHTH4E(ielmc)+(1._r8-FHVSH)*WTSTDE(ielmc,M,NZ)
      WTHTH4E(ielmn)=WTHTH4E(ielmn)+(1._r8-FHVSH)*WTSTDE(ielmn,M,NZ)
      WTHTH4E(ielmp)=WTHTH4E(ielmp)+(1._r8-FHVSH)*WTSTDE(ielmp,M,NZ)
      WTHTX4E(ielmc)=WTHTX4E(ielmc)+(FHVSH-FHVSE(ielmc))*WTSTDE(ielmc,M,NZ)
      WTHTX4E(ielmn)=WTHTX4E(ielmn)+(FHVSH-FHVSE(ielmc))*WTSTDE(ielmn,M,NZ)
      WTHTX4E(ielmp)=WTHTX4E(ielmp)+(FHVSH-FHVSE(ielmc))*WTSTDE(ielmp,M,NZ)
      WTSTDE(ielmc,M,NZ)=FHVSE(ielmc)*WTSTDE(ielmc,M,NZ)
      WTSTDE(ielmn,M,NZ)=FHVSE(ielmc)*WTSTDE(ielmn,M,NZ)
      WTSTDE(ielmp,M,NZ)=FHVSE(ielmc)*WTSTDE(ielmp,M,NZ)
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

  real(r8) :: WTHEL0(NumOfPlantChemElements)
  real(r8) :: WTHEL1(NumOfPlantChemElements)
  real(r8) :: WTHEL2(NumOfPlantChemElements)
  real(r8) :: WTHEL3(NumOfPlantChemElements)
  real(r8) :: WTHEL4(NumOfPlantChemElements)
  real(r8) :: WTHER0(NumOfPlantChemElements)
  real(r8) :: WTHERT(NumOfPlantChemElements)
  real(r8) :: WTHEXT(NumOfPlantChemElements)

  WTHEL0(1:NumOfPlantChemElements)=0._r8
  WTHEL1(1:NumOfPlantChemElements)=0._r8
  WTHEL2(1:NumOfPlantChemElements)=0._r8
  WTHEL3(1:NumOfPlantChemElements)=0._r8
  WTHEL4(1:NumOfPlantChemElements)=0._r8
  WTHER0(1:NumOfPlantChemElements)=0._r8

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
  real(r8), intent(in) :: WTHEXT(NumOfPlantChemElements),WTHERT(NumOfPlantChemElements)
  real(r8), intent(in) :: WTHER0(NumOfPlantChemElements)
  real(r8), intent(in) :: WTHEL0(NumOfPlantChemElements)
  real(r8), intent(in) :: WTHEL1(NumOfPlantChemElements)
  real(r8), intent(in) :: WTHEL2(NumOfPlantChemElements)
  real(r8), intent(in) :: WTHEL3(NumOfPlantChemElements)
  real(r8), intent(in) :: WTHEL4(NumOfPlantChemElements)
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
    WTSTDE   =>  plt_biom%WTSTDE  , &
    IHVST    =>  plt_distb%IHVST  , &
    FWOODE   =>  plt_allom%FWOODE , &
    CFOPE    =>  plt_soilchem%CFOPE, &
    ESNC     =>  plt_bgcr%ESNC    , &
    TESNC    =>  plt_bgcr%TESNC   , &
    TESN0    =>  plt_bgcr%TESN0   , &
    IBTYP    =>  plt_pheno%IBTYP  , &
    IGTYP    =>  plt_pheno%IGTYP    &
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
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      D6375: DO M=1,jsken
        DO NE=1,NumOfPlantChemElements        
          ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,instruct,M,NZ)*(WTHER0(NE)+WTHTX0E(NE)) &
            +CFOPE(NE,ifoliar,M,NZ)*(WTHTR1E(NE)+WTHTX1E(NE)) &
            +CFOPE(NE,infoliar,M,NZ)*(WTHTR2E(NE)+WTHTX2E(NE))

          IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,istalk,M,NZ)*(WTHTR3E(NE)+WTHTX3E(NE)+WTHTR4E(NE)+WTHTX4E(NE))
          ELSE
            WTSTDE(NE,M,NZ)=WTSTDE(NE,M,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WTHTX3E(NE)+WTHTX4E(NE))

            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WTHTR3E(NE)+WTHTR4E(NE))*FWOODE(NE,k_woody_litr)

            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
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
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    ELSE
      D6485: DO M=1,jsken
        ESNC(ielmc,M,k_fine_litr,0,NZ)=ESNC(ielmc,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmc,instruct,M,NZ)*(WTHER0(ielmc)+WTHTX0E(ielmc)) &
          +CFOPE(ielmc,ifoliar,M,NZ)*(WTHTR1E(ielmc)+WTHTX1E(ielmc)) &
          +CFOPE(ielmc,infoliar,M,NZ)*(WTHTR2E(ielmc)+WTHTX2E(ielmc))
        ESNC(ielmn,M,k_fine_litr,0,NZ)=ESNC(ielmn,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmn,instruct,M,NZ)*WTHEL0(ielmn) &
          +CFOPE(ielmn,ifoliar,M,NZ)*WTHEL1(ielmn) &
          +CFOPE(ielmn,infoliar,M,NZ)*WTHEL2(ielmn)
        ESNC(ielmp,M,k_fine_litr,0,NZ)=ESNC(ielmp,M,k_fine_litr,0,NZ) &
          +CFOPE(ielmp,instruct,M,NZ)*WTHEL0(ielmp) &
          +CFOPE(ielmp,ifoliar,M,NZ)*WTHEL1(ielmp) &
          +CFOPE(ielmp,infoliar,M,NZ)*WTHEL2(ielmp)

        ESNC(ielmn,ilignin,k_fine_litr,0,NZ)=ESNC(ielmn,ilignin,k_fine_litr,0,NZ) &
          +CFOPE(ielmn,instruct,M,NZ)*(WTHER0(ielmn)+WTHTX0E(ielmn)-WTHEL0(ielmn)) &
          +CFOPE(ielmn,ifoliar,M,NZ)*(WTHTR1E(ielmn)+WTHTX1E(ielmn)-WTHEL1(ielmn)) &
          +CFOPE(ielmn,infoliar,M,NZ)*(WTHTR2E(ielmn)+WTHTX2E(ielmn)-WTHEL2(ielmn))
        ESNC(ielmp,ilignin,k_fine_litr,0,NZ)=ESNC(ielmp,ilignin,k_fine_litr,0,NZ) &
          +CFOPE(ielmp,instruct,M,NZ)*(WTHER0(ielmp)+WTHTX0E(ielmp)-WTHEL0(ielmp)) &
          +CFOPE(ielmp,ifoliar,M,NZ)*(WTHTR1E(ielmp)+WTHTX1E(ielmp)-WTHEL1(ielmp)) &
          +CFOPE(ielmp,infoliar,M,NZ)*(WTHTR2E(ielmp)+WTHTX2E(ielmp)-WTHEL2(ielmp))

        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          ESNC(ielmc,M,k_fine_litr,0,NZ)=ESNC(ielmc,M,k_fine_litr,0,NZ)+CFOPE(ielmc,istalk,M,NZ) &
            * (WTHTR3E(ielmc)+WTHTX3E(ielmc)+WTHTR4E(ielmc)+WTHTX4E(ielmc))
          ESNC(ielmn,M,k_fine_litr,0,NZ)=ESNC(ielmn,M,k_fine_litr,0,NZ)+CFOPE(ielmn,istalk,M,NZ) &
            *(WTHEL3(ielmn)+WTHEL4(ielmn))
          ESNC(ielmp,M,k_fine_litr,0,NZ)=ESNC(ielmp,M,k_fine_litr,0,NZ)+CFOPE(ielmp,istalk,M,NZ) &
            *(WTHEL3(ielmp)+WTHEL4(ielmp))
          ESNC(ielmn,ilignin,k_fine_litr,0,NZ)=ESNC(ielmn,ilignin,k_fine_litr,0,NZ) &
            +CFOPE(ielmn,istalk,M,NZ)*(WTHTR3E(ielmn)+WTHTX3E(ielmn) &
            -WTHEL3(ielmn)+WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHEL4(ielmn))
          ESNC(ielmp,ilignin,k_fine_litr,0,NZ)=ESNC(ielmp,ilignin,k_fine_litr,0,NZ)+&
            CFOPE(ielmp,istalk,M,NZ)*(WTHTR3E(ielmp)+WTHTX3E(ielmp)- &
            WTHEL3(ielmp)+WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHEL4(ielmp))
        ELSE
          WTSTDE(ielmc,M,NZ)=WTSTDE(ielmc,M,NZ)+CFOPE(ielmc,icwood,M,NZ) &
            *(WTHTR3E(ielmc)+WTHTX3E(ielmc))
          WTSTDE(ielmn,M,NZ)=WTSTDE(ielmn,M,NZ)+CFOPE(ielmn,icwood,M,NZ)*WTHEL3(ielmn)
          WTSTDE(ielmp,M,NZ)=WTSTDE(ielmp,M,NZ)+CFOPE(ielmp,icwood,M,NZ)*WTHEL3(ielmp)
          ESNC(ielmc,M,k_woody_litr,0,NZ)=ESNC(ielmc,M,k_woody_litr,0,NZ)*CFOPE(ielmc,istalk,M,NZ)&
            *(WTHTR4E(ielmc)+WTHTX4E(ielmc))*FWOODE(ielmc,k_woody_litr)
          ESNC(ielmn,M,k_woody_litr,0,NZ)=ESNC(ielmn,M,k_woody_litr,0,NZ)+CFOPE(ielmn,istalk,M,NZ)*WTHEL4(ielmn)*FWOODE(ielmn,k_woody_litr)
          ESNC(ielmp,M,k_woody_litr,0,NZ)=ESNC(ielmp,M,k_woody_litr,0,NZ)+CFOPE(ielmp,istalk,M,NZ)*WTHEL4(ielmp)*FWOODE(ielmp,k_woody_litr)
          ESNC(ielmn,ilignin,k_woody_litr,0,NZ)=ESNC(ielmn,ilignin,k_woody_litr,0,NZ)+CFOPE(ielmn,icwood,M,NZ) &
            *(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHEL3(ielmn) &
            +WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHEL4(ielmn))*FWOODE(ielmn,k_woody_litr)
          ESNC(ielmp,ilignin,k_woody_litr,0,NZ)=ESNC(ielmp,ilignin,k_woody_litr,0,NZ)+CFOPE(ielmp,icwood,M,NZ) &
            *(WTHTR3E(ielmp)+WTHTX3E(ielmp)-WTHEL3(ielmp) &
            +WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHEL4(ielmp))*FWOODE(ielmp,k_woody_litr)
          ESNC(ielmc,M,k_fine_litr,0,NZ)=ESNC(ielmc,M,k_fine_litr,0,NZ)+CFOPE(ielmc,istalk,M,NZ) &
            *(WTHTR4E(ielmc)+WTHTX4E(ielmc))*FWOODE(ielmc,k_fine_litr)
          ESNC(ielmn,M,k_fine_litr,0,NZ)=ESNC(ielmn,M,k_fine_litr,0,NZ)+CFOPE(ielmn,istalk,M,NZ)*WTHEL4(ielmn)*FWOODE(ielmn,k_fine_litr)
          ESNC(ielmp,M,k_fine_litr,0,NZ)=ESNC(ielmp,M,k_fine_litr,0,NZ)+CFOPE(ielmp,istalk,M,NZ)*WTHEL4(ielmp)*FWOODE(ielmp,k_fine_litr)
          ESNC(ielmn,ilignin,k_fine_litr,0,NZ)=ESNC(ielmn,ilignin,k_fine_litr,0,NZ)+CFOPE(ielmn,icwood,M,NZ) &
            *(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHEL3(ielmn) &
            +WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHEL4(ielmn))*FWOODE(ielmn,k_fine_litr)
          ESNC(ielmp,ilignin,k_fine_litr,0,NZ)=ESNC(ielmp,ilignin,k_fine_litr,0,NZ)+CFOPE(ielmp,icwood,M,NZ) &
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
    DO NE=1,NumOfPlantChemElements
      TESNC(NE,NZ)=TESNC(NE,NZ)+WTHERT(NE)+WTHEXT(NE)
      TESN0(NE,NZ)=TESN0(NE,NZ)+WTHERT(NE)+WTHEXT(NE)
    ENDDO
  ENDIF
  end associate
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

  subroutine TotalBiomRemovalByDisturbance(I,J,NZ,WTHER0,WTHEXT,WTHERT)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: WTHER0(NumOfPlantChemElements)
  real(r8), intent(out):: WTHEXT(NumOfPlantChemElements),WTHERT(NumOfPlantChemElements)
  real(r8) :: WTHEHT(NumOfPlantChemElements)
  integer :: NE
!     begin_execution
  associate(                            &
    IHVST    =>  plt_distb%IHVST  , &
    JHVST    =>  plt_distb%JHVST  , &
    HVSTE    =>  plt_distb%HVSTE  , &
    XHVSTE   =>  plt_distb%XHVSTE , &
    VNH3F    =>  plt_distb%VNH3F  , &
    VPO4F    =>  plt_distb%VPO4F  , &
    VCH4F    =>  plt_distb%VCH4F  , &
    VOXYF    =>  plt_distb%VOXYF  , &
    VN2OF    =>  plt_distb%VN2OF  , &
    VCO2F    =>  plt_distb%VCO2F  , &
    CO2NetFix_pft     =>  plt_bgcr%CO2NetFix_pft    , &
    TNBP     =>  plt_bgcr%TNBP    , &
    TRAU     =>  plt_bgcr%TRAU    , &
    RECO     =>  plt_bgcr%RECO    , &
    TCO2A    =>  plt_bgcr%TCO2A   , &
    TCO2T    =>  plt_bgcr%TCO2T   , &
    WTRVE    =>  plt_biom%WTRVE     &
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
  DO NE=1,NumOfPlantChemElements
    WTHEHT(NE)=WTHTH0E(NE)+WTHTH1E(NE)+WTHTH2E(NE)+WTHTH3E(NE)+WTHTH4E(NE)
    WTHERT(NE)=WTHER0(NE)+WTHTR1E(NE)+WTHTR2E(NE)+WTHTR3E(NE)+WTHTR4E(NE)
    WTHEXT(NE)=WTHTX0E(NE)+WTHTX1E(NE)+WTHTX2E(NE)+WTHTX3E(NE)+WTHTX4E(NE)
  ENDDO
  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      IF(JHVST(NZ).NE.2)THEN
        DO NE=1,NumOfPlantChemElements
          HVSTE(NE,NZ)=HVSTE(NE,NZ)+WTHEHT(NE)-WTHERT(NE)
          XHVSTE(NE)=XHVSTE(NE)+WTHEHT(NE)-WTHERT(NE)
        ENDDO
        TNBP=TNBP+WTHERT(ielmc)-WTHEHT(ielmc)
      ELSE
        DO NE=1,NumOfPlantChemElements
          WTRVE(NE,NZ)=WTRVE(NE,NZ)+WTHEHT(NE)-WTHERT(NE)
        ENDDO
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     TNBP=total net biome productivity
!
    ELSE
      VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*(WTHEHT(ielmc)-WTHERT(ielmc))
      VCH4F(NZ)=VCH4F(NZ)-FCH4F*(WTHEHT(ielmc)-WTHERT(ielmc))
      VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*(WTHEHT(ielmc)-WTHERT(ielmc))*2.667
      VNH3F(NZ)=VNH3F(NZ)-WTHEHT(ielmn)+WTHERT(ielmn)
      VN2OF(NZ)=VN2OF(NZ)-0.0_r8
      VPO4F(NZ)=VPO4F(NZ)-WTHEHT(ielmp)+WTHERT(ielmp)
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*(WTHEHT(ielmc)-WTHERT(ielmc))
      TNBP=TNBP-FCH4F*(WTHEHT(ielmc)-WTHERT(ielmc))
    ENDIF
!
!     C,N,P REMOVED FROM GRAZING
!
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     GY=growth yield of grazers
!     WTHEHT(ielmc),WTHEHT(ielmn),WTHEHT(ielmp)=total C,N,P removed
!     WTHERT(ielmc),WTHERT(ielmn),WTHERT(ielmp)=total C,N,P to litter
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
  ELSE
    HVSTE(ielmc,NZ)=HVSTE(ielmc,NZ)+GY*(WTHEHT(ielmc)-WTHERT(ielmc))
    XHVSTE(ielmc)=XHVSTE(ielmc)+GY*(WTHEHT(ielmc)-WTHERT(ielmc))
    DO NE=2,NumOfPlantChemElements
      HVSTE(NE,NZ)=HVSTE(NE,NZ)+WTHEHT(NE)-WTHERT(NE)
      XHVSTE(NE)=XHVSTE(NE)+WTHEHT(NE)-WTHERT(NE)
    ENDDO
    TCO2T(NZ)=TCO2T(NZ)-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
    TCO2A(NZ)=TCO2A(NZ)-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
!     TNBP=TNBP+GY*(WTHERT(ielmc)-WTHEHT(ielmc))
!     CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)+GZ*(WTHERT(ielmc)-WTHEHT(ielmc))
    RECO=RECO-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
    TRAU=TRAU-GZ*(WTHEHT(ielmc)-WTHERT(ielmc))
  ENDIF
  end associate
  end subroutine TotalBiomRemovalByDisturbance
!------------------------------------------------------------------------------------------

  subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,WTHER0,&
    WTHEL0,WTHEL1,WTHEL2,WTHEL3,WTHEL4)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: WTHER0(NumOfPlantChemElements)
  real(r8), intent(out) :: WTHEL0(NumOfPlantChemElements)
  real(r8), intent(out) :: WTHEL1(NumOfPlantChemElements)
  real(r8), intent(out) :: WTHEL2(NumOfPlantChemElements)
  real(r8), intent(out) :: WTHEL3(NumOfPlantChemElements)
  real(r8), intent(out) :: WTHEL4(NumOfPlantChemElements)

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
    DO NE=1,NumOfPlantChemElements
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
    DO NE=1,NumOfPlantChemElements
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
    DO NE=1,NumOfPlantChemElements
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
    DO NE=1,NumOfPlantChemElements
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

    DO NE=2,NumOfPlantChemElements
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
    IDAYH    =>  plt_distb%IDAYH     , &
    IDAY0    =>  plt_distb%IDAY0     , &
    ITILL    =>  plt_distb%ITILL     , &
    IYR0     =>  plt_distb%IYR0      , &
    IYRH     =>  plt_distb%IYRH      , &
    XCORP    =>  plt_distb%XCORP     , &
    CFOPE    =>  plt_soilchem%CFOPE  , &
    trcg_rootml     =>  plt_rbgc%trcg_rootml       , &
    trcs_rootml => plt_rbgc%trcs_rootml , &
    UVOLO    =>  plt_ew%UVOLO        , &
    CanWatP    =>  plt_ew%CanWatP        , &
    VHeatCapCanP    =>  plt_ew%VHeatCapCanP        , &
    PSICanP    =>  plt_ew%PSICanP        , &
    PPX      =>  plt_site%PPX        , &
    EPOOLR   =>  plt_biom%EPOOLR     , &
    WSRTL    =>  plt_biom%WSRTL      , &
    PopPlantRootC_vr    =>  plt_biom%PopPlantRootC_vr      , &
    WTRTL    =>  plt_biom%WTRTL      , &
    WTLFBE   =>  plt_biom%WTLFBE     , &
    WTGRBE   =>  plt_biom%WTGRBE     , &
    WTEARBE  =>  plt_biom%WTEARBE    , &
    EPOLNB   =>  plt_biom%EPOLNB     , &
    EPOOL    =>  plt_biom%EPOOL      , &
    WTHSKBE  =>  plt_biom%WTHSKBE    , &
    WTRSVBE  =>  plt_biom%WTRSVBE    , &
    WTNDBE   =>  plt_biom%WTNDBE     , &
    WTSHTBE  =>  plt_biom%WTSHTBE    , &
    WTSTKBE  =>  plt_biom%WTSTKBE    , &
    WGLFLE   =>  plt_biom%WGLFLE     , &
    WTSHEBE  =>  plt_biom%WTSHEBE    , &
    WGLFE    =>  plt_biom%WGLFE      , &
    CanPBLeafShethC    =>  plt_biom%CanPBLeafShethC      , &
    WTSTXBE  =>  plt_biom%WTSTXBE    , &
    CanPBStalkC   =>  plt_biom%CanPBStalkC     , &
    WSSHE    =>  plt_biom%WSSHE      , &
    WGSHE    =>  plt_biom%WGSHE      , &
    WTRT1E   =>  plt_biom%WTRT1E     , &
    WSLF     =>  plt_biom%WSLF       , &
    WGNODE   =>  plt_biom%WGNODE     , &
    CanPStalkC    =>  plt_biom%CanPStalkC      , &
    RTWT1E   =>  plt_biom%RTWT1E     , &
    WTRVE    =>  plt_biom%WTRVE      , &
    CanopyLeafShethC_pft     =>  plt_biom%CanopyLeafShethC_pft       , &
    WTRT2E   =>  plt_biom%WTRT2E     , &
    EPOOLN   =>  plt_biom%EPOOLN     , &
    WTNDLE   =>  plt_biom%WTNDLE     , &
    GRWTB    =>  plt_allom%GRWTB     , &
    FWOODE   =>  plt_allom%FWOODE    , &
    FWODBE   =>  plt_allom%FWODBE    , &
    FWODLE   =>  plt_allom%FWODLE    , &
    FWODRE   =>  plt_allom%FWODRE    , &
    IDTHB    =>  plt_pheno%IDTHB     , &
    ISTYP    =>  plt_pheno%ISTYP     , &
    IDTH     =>  plt_pheno%IDTH      , &
    IGTYP    =>  plt_pheno%IGTYP     , &
    IBTYP    =>  plt_pheno%IBTYP     , &
    IWTYP    =>  plt_pheno%IWTYP     , &
    IDTHR    =>  plt_pheno%IDTHR     , &
    IDTHP    =>  plt_pheno%IDTHP     , &
    HCOB     =>  plt_photo%HCOB      , &
    CO2B     =>  plt_photo%CO2B      , &
    CPOOL3   =>  plt_photo%CPOOL3    , &
    CPOOL4   =>  plt_photo%CPOOL4    , &
    NJ       =>  plt_site%NJ         , &
    pftPlantPopulation       =>  plt_site%pftPlantPopulation         , &
    IYRC     =>  plt_site%IYRC       , &
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
    ESNC     =>  plt_bgcr%ESNC       , &
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
    INTYP    =>  plt_morph%INTYP     , &
    RootLenPerP    =>  plt_morph%RootLenPerP     , &
    RTN2     =>  plt_morph%RTN2      , &
    SecndRootLen    =>  plt_morph%SecndRootLen     , &
    SecndRootXNumL     =>  plt_morph%SecndRootXNumL      , &
    NGTopRootLayer      =>  plt_morph%NGTopRootLayer       , &
    MY       =>  plt_morph%MY        , &
    NRT      =>  plt_morph%NRT       , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft       , &
    ARLF1    =>  plt_morph%ARLF1     , &
    CanopyBranchLeafA_pft    =>  plt_morph%CanopyBranchLeafA_pft     , &
    GRNXB    =>  plt_morph%GRNXB     , &
    GRNOB    =>  plt_morph%GRNOB     , &
    CanPLNBLA    =>  plt_morph%CanPLNBLA       &
  )
!     ZNOON=hour of solar noon
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IDAY0,IYR0=day,year of planting
!     IYRC=current year
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!     PPX,PP=PFT population per m2,grid cell
!     FracPARByCanP=fraction of radiation received by each PFT canopy
!     VHeatCapCanP=canopy heat capacity
!
  IF(J.EQ.INT(ZNOON).AND.(IBTYP(NZ).EQ.0 &
    .OR.IGTYP(NZ).LE.1).AND.(I.NE.IDAY0(NZ) &
    .OR.IYRC.NE.IYR0(NZ)))THEN
    IF(ITILL.LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.IDAY0(NZ).OR.IYRC.GT.IYR0(NZ))THEN
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
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
        D8975: DO NB=1,NumOfBranches_pft(NZ)
          IF(IDTHB(NB,NZ).EQ.ibralive)THEN
            IF(pftPlantPopulation(NZ).LE.0.0)then
              IDTHB(NB,NZ)=ibrdead
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
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
            XHVST1=1._r8-XHVST
            D6380: DO M=1,jsken
              ESNC(ielmc,M,k_fine_litr,0,NZ)=ESNC(ielmc,M,k_fine_litr,0,NZ)+XHVST1 &
                *(CFOPE(ielmc,instruct,M,NZ)*(EPOOL(ielmc,NB,NZ)+EPOLNB(ielmc,NB,NZ) &
                +CPOOLK(NB,NZ)+WTRSVBE(ielmc,NB,NZ)) &
                +CFOPE(ielmc,ifoliar,M,NZ)*(WTLFBE(ielmc,NB,NZ)*FWODLE(ielmc,k_fine_litr) &
                +WTNDBE(ielmc,NB,NZ)) &
                +CFOPE(ielmc,infoliar,M,NZ)*(WTSHEBE(ielmc,NB,NZ)*FWODBE(ielmc,k_fine_litr) &
                +WTHSKBE(ielmc,NB,NZ)+WTEARBE(ielmc,NB,NZ)))

              DO NE=2,NumOfPlantChemElements
                ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *(CFOPE(NE,instruct,M,NZ)*(EPOOL(NE,NB,NZ)+EPOLNB(NE,NB,NZ)&
                  +WTRSVBE(NE,NB,NZ)) &
                  + CFOPE(NE,ifoliar,M,NZ)*(WTLFBE(NE,NB,NZ)*FWODLE(NE,k_fine_litr)+WTNDBE(NE,NB,NZ)) &
                  + CFOPE(NE,infoliar,M,NZ)*(WTSHEBE(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
                  +   WTHSKBE(NE,NB,NZ)+WTEARBE(NE,NB,NZ)))
              ENDDO
            ENDDO D6380

            DO M=1,jsken
              DO NE=1,NumOfPlantChemElements
                ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,icwood,M,NZ)*(WTLFBE(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
                  +WTSHEBE(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

                IF(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
                  WTRVE(NE,NZ)=WTRVE(NE,NZ)+XHVST1*CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
                ELSE
                  ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                    *CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
                ENDIF
                ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,icwood,M,NZ)*WTSTKBE(NE,NB,NZ)*FWOODE(NE,k_woody_litr)

                ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+XHVST1 &
                  *CFOPE(NE,istalk,M,NZ)*WTSTKBE(NE,NB,NZ)*FWOODE(NE,k_fine_litr)
              ENDDO
            ENDDO
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!

            CPOOLK(NB,NZ)=CPOOLK(NB,NZ)*XHVST
            CanPBStalkC(NB,NZ)=CanPBStalkC(NB,NZ)*XHVST
            DO NE=1,NumOfPlantChemElements
              EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)*XHVST
              EPOLNB(NE,NB,NZ)=EPOLNB(NE,NB,NZ)*XHVST
              WTSHTBE(NE,NB,NZ)=WTSHTBE(NE,NB,NZ)*XHVST
              WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ)*XHVST
              WTHSKBE(NE,NB,NZ)=WTHSKBE(NE,NB,NZ)*XHVST
              WTEARBE(NE,NB,NZ)=WTEARBE(NE,NB,NZ)*XHVST
              WTGRBE(NE,NB,NZ)=WTGRBE(NE,NB,NZ)*XHVST
              WTLFBE(NE,NB,NZ)=WTLFBE(NE,NB,NZ)*XHVST
              WTNDBE(NE,NB,NZ)=WTNDBE(NE,NB,NZ)*XHVST
              WTSHEBE(NE,NB,NZ)=WTSHEBE(NE,NB,NZ)*XHVST
              WTSTKBE(NE,NB,NZ)=WTSTKBE(NE,NB,NZ)*XHVST
              WTSTXBE(NE,NB,NZ)=WTSTXBE(NE,NB,NZ)*XHVST
            ENDDO

            GRNXB(NB,NZ)=GRNXB(NB,NZ)*XHVST
            GRNOB(NB,NZ)=GRNOB(NB,NZ)*XHVST
            GRWTB(NB,NZ)=GRWTB(NB,NZ)*XHVST
            CanopyBranchLeafA_pft(NB,NZ)=CanopyBranchLeafA_pft(NB,NZ)*XHVST
            CanPBLeafShethC(NB,NZ)=AZMAX1(WTLFBE(ielmc,NB,NZ)+WTSHEBE(ielmc,NB,NZ))
            CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+CanPBLeafShethC(NB,NZ)

            CanPStalkC(NZ)=CanPStalkC(NZ)+CanPBStalkC(NB,NZ)
            D8970: DO K=0,JNODS1
              IF(K.NE.0)THEN
                CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)*XHVST
                CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)*XHVST
                CO2B(K,NB,NZ)=CO2B(K,NB,NZ)*XHVST
                HCOB(K,NB,NZ)=HCOB(K,NB,NZ)*XHVST
              ENDIF
              ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)*XHVST

              WSLF(K,NB,NZ)=WSLF(K,NB,NZ)*XHVST
!     CanPSheathHeight(K,NB,NZ)=CanPSheathHeight(K,NB,NZ)*XHVST

              WSSHE(K,NB,NZ)=WSSHE(K,NB,NZ)*XHVST
!     HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)*XHVST
!     HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ)*XHVST
              DO NE=1,NumOfPlantChemElements
                WGNODE(NE,K,NB,NZ)=WGNODE(NE,K,NB,NZ)*XHVST
                WGLFE(NE,K,NB,NZ)=WGLFE(NE,K,NB,NZ)*XHVST
                WGSHE(NE,K,NB,NZ)=WGSHE(NE,K,NB,NZ)*XHVST
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
!     IDTHR,IDTHP=PFT root,shoot living flag: 0=alive,1=dead
!     IDTH=PFT living flag: 0=alive,1=dead
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     IDAYH,IYRH=day,year of harvesting
!     IYRC=current year
!
        IF(pftPlantPopulation(NZ).LE.0.0)THEN
          IDTHR(NZ)=ibrdead
          IDTHP(NZ)=ibrdead
          IDTH(NZ)=ibrdead
          JHVST(NZ)=ihv_terminate
          IDAYH(NZ)=I
          IYRH(NZ)=IYRC
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
                DO NE=1,NumOfPlantChemElements
                  ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,instruct,M,NZ)*EPOOLR(NE,N,L,NZ)
                ENDDO

              DO NR=1,NRT(NZ)
                DO NE=1,NumOfPlantChemElements
                  ESNC(NE,M,k_woody_litr,L,NZ)=ESNC(NE,M,k_woody_litr,L,NZ)+XHVST1 &
                    *CFOPE(NE,icwood,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

                  ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+XHVST1 &
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
!     WTRTL,PopPlantRootC_vr=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,SecndRootXNumL=number of primary,secondary root axes
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D8960: DO NR=1,NRT(NZ)
              DO NE=1,NumOfPlantChemElements
                WTRT1E(NE,N,L,NR,NZ)=WTRT1E(NE,N,L,NR,NZ)*XHVST
                WTRT2E(NE,N,L,NR,NZ)=WTRT2E(NE,N,L,NR,NZ)*XHVST
                RTWT1E(NE,N,NR,NZ)=RTWT1E(NE,N,NR,NZ)*XHVST
              ENDDO
              PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)*XHVST
              SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)*XHVST
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST
            ENDDO D8960
            DO NE=1,NumOfPlantChemElements
              EPOOLR(NE,N,L,NZ)=EPOOLR(NE,N,L,NZ)*XHVST
            ENDDO
            WTRTL(N,L,NZ)=WTRTL(N,L,NZ)*XHVST
            PopPlantRootC_vr(N,L,NZ)=PopPlantRootC_vr(N,L,NZ)*XHVST
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
            IF(INTYP(NZ).NE.0.AND.N.EQ.1)THEN
              DO NE=1,NumOfPlantChemElements
                D6395: DO M=1,jsken
                  ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+XHVST1 &
                    *(CFOPE(NE,iroot,M,NZ)*WTNDLE(NE,L,NZ) &
                    +CFOPE(NE,instruct,M,NZ)*EPOOLN(NE,L,NZ))
                ENDDO D6395
                WTNDLE(NE,L,NZ)=WTNDLE(NE,L,NZ)*XHVST
                EPOOLN(NE,L,NZ)=EPOOLN(NE,L,NZ)*XHVST
              ENDDO
            ENDIF
          ENDDO D8985
        ENDDO D8980
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
!     DURING TILLAGE
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        DO NE=1,NumOfPlantChemElements
          D6400: DO M=1,jsken
            ESNC(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ)=ESNC(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ) &
              +(XHVST1*CFOPE(NE,instruct,M,NZ)*WTRVE(NE,NZ))*FWOODE(NE,k_woody_litr)

            ESNC(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ)=ESNC(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ) &
              +(XHVST1*CFOPE(NE,instruct,M,NZ)*WTRVE(NE,NZ))*FWOODE(NE,k_fine_litr)
          ENDDO D6400
          WTRVE(NE,NZ)=WTRVE(NE,NZ)*XHVST
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
  real(r8) :: ZPOLNX,XHVST(NumOfPlantChemElements)
  real(r8) :: XHVST1(NumOfPlantChemElements)
  REAL(R8) :: WGLFBL(NumOfCanopyLayers1,JP1,JP1)
  real(r8) :: FHVSHK(0:JNODS1),FHVSETK(0:JNODS1)
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
  real(r8) :: FFIRE(NumOfPlantChemElements)
  real(r8) :: FHVSE(NumOfPlantChemElements)
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
  real(r8) :: WGLFGE(NumOfPlantChemElements)
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
    VCH4F    =>  plt_distb%VCH4F   , &
    VCO2F    =>  plt_distb%VCO2F   , &
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
    EPOOLR   => plt_biom%EPOOLR    , &
    WSRTL    => plt_biom%WSRTL     , &
    PopPlantRootC_vr    => plt_biom%PopPlantRootC_vr     , &
    WTRTL    => plt_biom%WTRTL     , &
    WTRVE    => plt_biom%WTRVE     , &
    CanPStalkC    => plt_biom%CanPStalkC     , &
    CanopyLeafShethC_pft     => plt_biom%CanopyLeafShethC_pft      , &
    CanPBStalkC   => plt_biom%CanPBStalkC    , &
    WTNDBE   => plt_biom%WTNDBE    , &
    WTRSVE   => plt_biom%WTRSVE    , &
    WTGRBE   => plt_biom%WTGRBE    , &
    WTSTKBE  => plt_biom%WTSTKBE   , &
    WTSHTBE  => plt_biom%WTSHTBE   , &
    WTHSKBE  => plt_biom%WTHSKBE   , &
    WTEARBE  => plt_biom%WTEARBE   , &
    WGNODE   => plt_biom%WGNODE    , &
    CanPBLeafShethC  => plt_biom%CanPBLeafShethC     , &
    WTRSVBE  => plt_biom%WTRSVBE   , &
    WTSTXBE  => plt_biom%WTSTXBE   , &
    EPOLNB   => plt_biom%EPOLNB     , &
    WGLFLE   => plt_biom%WGLFLE    , &
    WTSHEBE  => plt_biom%WTSHEBE   , &
    EPOOL    => plt_biom%EPOOL     , &
    WGSHE    => plt_biom%WGSHE     , &
    WSLF     => plt_biom%WSLF      , &
    WGLFE    => plt_biom%WGLFE     , &
    WTLFBE   => plt_biom%WTLFBE    , &
    WSSHE    => plt_biom%WSSHE     , &
    WTSTKE   => plt_biom%WTSTKE    , &
    CanopyNonstructElementConc_pft   => plt_biom%CanopyNonstructElementConc_pft    , &
    NoduleNonstructCconc_pft   => plt_biom%NoduleNonstructCconc_pft    , &
    WTLFE    => plt_biom%WTLFE     , &
    WTGRE    => plt_biom%WTGRE     , &
    CanPShootElmMass   => plt_biom%CanPShootElmMass     , &
    WTHSKE   => plt_biom%WTHSKE    , &
    WTEARE   => plt_biom%WTEARE    , &
    WTSHEE   => plt_biom%WTSHEE    , &
    WTSHTA   => plt_biom%WTSHTA    , &
    WTRT1E   => plt_biom%WTRT1E    , &
    RTWT1E   => plt_biom%RTWT1E    , &
    WTRT2E   => plt_biom%WTRT2E    , &
    WTNDLE   => plt_biom%WTNDLE    , &
    EPOOLN   => plt_biom%EPOOLN    , &
    ZEROP    => plt_biom%ZEROP     , &
    ZEROL    => plt_biom%ZEROL     , &
    CanopyLeafCpft_lyr    => plt_biom%CanopyLeafCpft_lyr     , &
    FVRN     => plt_allom%FVRN     , &
    FWODRE   => plt_allom%FWODRE   , &
    FWOODE   => plt_allom%FWOODE   , &
    FWODBE   => plt_allom%FWODBE   , &
    FWODLE   => plt_allom%FWODLE   , &
    GRWTB    => plt_allom%GRWTB    , &
    IDTHB    =>  plt_pheno%IDTHB   , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP    , &
    IDAY     =>  plt_pheno%IDAY    , &
    GROUP    =>  plt_pheno%GROUP   , &
    VSTGX    =>  plt_pheno%VSTGX   , &
    IGTYP    =>  plt_pheno%IGTYP   , &
    IBTYP    =>  plt_pheno%IBTYP   , &
    IFLGA    =>  plt_pheno%IFLGA   , &
    ISTYP    =>  plt_pheno%ISTYP   , &
    VRNX     =>  plt_pheno%VRNX    , &
    VRNF     =>  plt_pheno%VRNF    , &
    IWTYP    =>  plt_pheno%IWTYP   , &
    TGSTGI   =>  plt_pheno%TGSTGI  , &
    TGSTGF   =>  plt_pheno%TGSTGF  , &
    FLG4     =>  plt_pheno%FLG4    , &
    GROUPI   =>  plt_pheno%GROUPI  , &
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
    ESNC     =>  plt_bgcr%ESNC     , &
    TNBP     =>  plt_bgcr%TNBP     , &
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
    INTYP    =>  plt_morph%INTYP   , &
    CanopyLAgrid_lyr    =>  plt_morph%CanopyLAgrid_lyr   , &
    CanopyHeightz       =>  plt_morph%CanopyHeightz      , &
    CanopyBranchLeafA_pft    =>  plt_morph%CanopyBranchLeafA_pft   , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft     , &
    CanPSA    =>  plt_morph%CanPSA   , &
    HTNODX   =>  plt_morph%HTNODX  , &
    PrimRootXNumL    =>  plt_morph%PrimRootXNumL   , &
    SecndRootLen    =>  plt_morph%SecndRootLen   , &
    PrimRootLen    =>  plt_morph%PrimRootLen   , &
    HTNODE   =>  plt_morph%HTNODE  , &
    GRNXB    =>  plt_morph%GRNXB   , &
    GRNOB    =>  plt_morph%GRNOB   , &
    CanPSheathHeight    =>  plt_morph%CanPSheathHeight   , &
    ARLF1    =>  plt_morph%ARLF1   , &
    CanopyLeafApft_lyr    =>  plt_morph%CanopyLeafApft_lyr   , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr   , &
    CanPLNBLA    =>  plt_morph%CanPLNBLA   , &
    CanopyBranchStemApft_lyr    =>  plt_morph%CanopyBranchStemApft_lyr   , &
    NRT      =>  plt_morph%NRT     , &
    PSTGF    =>  plt_morph%PSTGF   , &
    NB1      =>  plt_morph%NB1     , &
    PSTGI    =>  plt_morph%PSTGI   , &
    PSTG     =>  plt_morph%PSTG    , &
    ClumpFactor      =>  plt_morph%ClumpFactor     , &
    CanopyLA_grd    =>  plt_morph%CanopyLA_grd   , &
    ICTYP    =>  plt_photo%ICTYP     &
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
          *AREA3(NU)*CanPShootElmMass(ielmc,NZ)/WTSHTA(NZ)
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
      WHVSLY=AMIN1(WTLFE(ielmc,NZ),WHVSLX)
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
      WTSHTT=WTSHEE(ielmc,NZ)+WTHSKE(ielmc,NZ)+WTEARE(ielmc,NZ)+WTGRE(ielmc,NZ)
      IF(WTSHTT.GT.ZEROP(NZ))THEN
        WHVSHX=WHVSSX*WTSHEE(ielmc,NZ)/WTSHTT+WHVXXX
        WHVSHY=AMIN1(WTSHEE(ielmc,NZ),WHVSHX)
        WHVSHH=WHVSHY*(1._r8-CCPOLX)
        WHVSCS=WHVSHY*CCPOLX
        WHVSNS=WHVSHY*CCPLNX
        WHVXXX=AZMAX1(WHVSHX-WHVSHY)
        WHVHSX=WHVSSX*WTHSKE(ielmc,NZ)/WTSHTT+WHVXXX
        WHVHSY=AMIN1(WTHSKE(ielmc,NZ),WHVHSX)
        WHVHSH=WHVHSY
        WHVXXX=AZMAX1(WHVHSX-WHVHSY)
        WHVEAX=WHVSSX*WTEARE(ielmc,NZ)/WTSHTT+WHVXXX
        WHVEAY=AMIN1(WTEARE(ielmc,NZ),WHVEAX)
        WHVEAH=WHVEAY
        WHVXXX=AZMAX1(WHVEAX-WHVEAY)
        WHVGRX=WHVSSX*WTGRE(ielmc,NZ)/WTSHTT+WHVXXX
        WHVGRY=AMIN1(WTGRE(ielmc,NZ),WHVGRX)
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
      WTSTKT=WTSTKE(ielmc,NZ)+WTRSVE(ielmc,NZ)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
        WHVSTX=WHVSKX*WTSTKE(ielmc,NZ)/WTSTKT+WHVXXX
        WHVSTY=AMIN1(WTSTKE(ielmc,NZ),WHVSTX)
        WHVSTH=WHVSTY
        WHVXXX=AZMAX1(WHVSTX-WHVSTY)
        WHVRVX=WHVSKX*WTRSVE(ielmc,NZ)/WTSTKT+WHVXXX
        WHVRVY=AMIN1(WTRSVE(ielmc,NZ),WHVRVX)
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
          WHVSLY=AMIN1(WTLFE(ielmc,NZ)-WHVSLF-WHVSCL,WHVXXX)
          WHVSLF=WHVSLF+WHVSLY*(1._r8-CCPOLX)
          WHVSCL=WHVSCL+WHVSLY*CCPOLX
          WHVSNL=WHVSNL+WHVSLY*CCPLNX
          WHVXXX=AZMAX1(WHVXXX-WHVSLY)
          IF(WTSHTT.GT.ZEROP(NZ))THEN
            WHVSHX=WHVXXX*WTSHEE(ielmc,NZ)/WTSHTT
            WHVSHY=AMIN1(WTSHEE(ielmc,NZ),WHVSHX)
            WHVSHH=WHVSHH+WHVSHY*(1._r8-CCPOLX)
            WHVSCS=WHVSCS+WHVSHY*CCPOLX
            WHVSNS=WHVSNS+WHVSHY*CCPLNX
            WHVXXX=AZMAX1(WHVXXX-WHVSHY)
            WHVHSX=WHVXXX*WTHSKE(ielmc,NZ)/WTSHTT
            WHVHSY=AMIN1(WTHSKE(ielmc,NZ),WHVHSX)
            WHVHSH=WHVHSH+WHVHSY
            WHVXXX=AZMAX1(WHVXXX-WHVHSY)
            WHVEAX=WHVXXX*WTEARE(ielmc,NZ)/WTSHTT
            WHVEAY=AMIN1(WTEARE(ielmc,NZ),WHVEAX)
            WHVEAH=WHVEAH+WHVEAY
            WHVXXX=AZMAX1(WHVEAX-WHVEAY)
            WHVGRX=WHVXXX*WTGRE(ielmc,NZ)/WTSHTT
            WHVGRY=AMIN1(WTGRE(ielmc,NZ),WHVGRX)
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
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=0._r8
          enddo
        enddo
      ENDDO D9860
      D9870: DO NB=1,NumOfBranches_pft(NZ)
        DO  L=1,NumOfCanopyLayers1
          DO  K=0,JNODS1
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
          .AND.WTLFE(ielmc,NZ).GT.ZEROL(NZ))THEN
          WHVSBL=WHVSLF*AZMAX1(WGLFBL(L,NB,NZ))/WTLFE(ielmc,NZ)
        ELSE
          WHVSBL=0._r8
        ENDIF
        D9845: DO K=JNODS1,0,-1
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
            DO NE=1,NumOfPlantChemElements
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
      D9825: DO K=0,JNODS1
        ARLFG=0._r8
        WGLFGE(1:NumOfPlantChemElements)=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     CanPLNBLA,CanopyLeafApft_lyr=leaf node,total area in canopy layer
!
        D9815: DO L=1,NumOfCanopyLayers1
          ARLFG=ARLFG+CanPLNBLA(L,K,NB,NZ)
          DO NE=1,NumOfPlantChemElements
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
          IF(WGLFE(ielmc,K,NB,NZ).GT.ZEROP(NZ).AND.EHVST(1,ipld_leaf,NZ).GT.0.0)THEN
            FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(1._r8-(1._r8-AZMAX1(WGLFGE(ielmc)) &
              /WGLFE(ielmc,K,NB,NZ))*EHVST(1,ipld_nofoliar,NZ)/EHVST(1,ipld_leaf,NZ))))
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
!     CanopyBranchLeafA_pft,ARLF=branch,node leaf area
!     WSLF=leaf protein mass
!
      WGLFGY=WGLFGY+WGLFE(ielmc,K,NB,NZ)
      DO NE=1,NumOfPlantChemElements
        WTLFBE(NE,NB,NZ)=WTLFBE(NE,NB,NZ)-WGLFE(NE,K,NB,NZ)+WGLFGE(NE)
      ENDDO
      CanopyBranchLeafA_pft(NB,NZ)=CanopyBranchLeafA_pft(NB,NZ)-ARLF1(K,NB,NZ)+ARLFG
      IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ))THEN
        WSLF(K,NB,NZ)=WSLF(K,NB,NZ)*ARLFG/ARLF1(K,NB,NZ)
      ELSE
        WSLF(K,NB,NZ)=0._r8
      ENDIF
      ARLF1(K,NB,NZ)=ARLFG
      DO NE=1,NumOfPlantChemElements
        WGLFE(NE,K,NB,NZ)=WGLFGE(NE)
      ENDDO
      WGLFGX=WGLFGX+WGLFE(ielmc,K,NB,NZ)
    ENDDO D9825
!
!     CUT SHEATHS OR PETIOLES AND STALKS HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WTSHE,WTSHEB=PFT,branch petiole C mass
!     WHVSBS,WHVSHH=branch, PFT petiole C mass removed
!     HTNODE=internode length
!     HTSTKX=internode length removed
!
      HTSTKX=0._r8
      IF((IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6) &
        .AND.WTSHEE(ielmc,NZ).GT.ZEROP(NZ))THEN
        WHVSBS=WHVSHH*WTSHEBE(ielmc,NB,NZ)/WTSHEE(ielmc,NZ)
      ELSE
        WHVSBS=0._r8
      ENDIF
      D9805: DO K=JNODS1,0,-1
!112   FORMAT(A8,8I4,12E12.4)
        IF(HTNODE(K,NB,NZ).GT.0.0) &
          HTSTKX=AMAX1(HTSTKX,HTNODE(K,NB,NZ))
!
!     HARVESTED SHEATH OR PETIOLE C,N,P
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WHVSBS=branch petiole C mass removed
!     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!     FHVSETK=fraction of internode layer mass not harvested
!     WTHTH2E(ielmc),WTHTH2E(ielmn),WTHTH2E(ielmp)=harvested petiole C,N,P
!     WTHTX2E(ielmc),WTHTX2E(ielmn),WTHTX2E(ielmp)=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     CanPSheathHeight,HTNODE=petiole,internode length
!
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBS.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(WGSHE(ielmc,K,NB,NZ).GT.WHVSBS)THEN
                FHVSETK(K)=AZMAX1(AMIN1(1.0_r8,(WGSHE(ielmc,K,NB,NZ)-WHVSBS)/WGSHE(ielmc,K,NB,NZ)))
                FHVSHK(K)=FHVSETK(K)
              ELSE
                FHVSETK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1._r8-FHVSETK(K))*WGSHE(ielmc,K,NB,NZ)
            DO NE=1,NumOfPlantChemElements
              WTHTH2E(NE)=WTHTH2E(NE)+(1._r8-FHVSHK(K))*WGSHE(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)
              WTHTX2E(NE)=WTHTX2E(NE)+(FHVSHK(K)-FHVSETK(K))*WGSHE(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)

              WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSHK(K))*WGSHE(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
              WTHTX3E(NE)=WTHTX3E(NE)+(FHVSHK(K)-FHVSETK(K))*WGSHE(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
            ENDDO
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     WGSHE=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     CanPSheathHeight=node petiole height
!     WSSHE=petiole protein mass
!
            WGSHGY=WGSHGY+WGSHE(ielmc,K,NB,NZ)
            WSSHE(K,NB,NZ)=FHVSETK(K)*WSSHE(K,NB,NZ)

            DO NE=1,NumOfPlantChemElements
              WTSHEBE(NE,NB,NZ)=WTSHEBE(NE,NB,NZ) &
                -(1._r8-FHVSETK(K))*WGSHE(NE,K,NB,NZ)
              WGSHE(NE,K,NB,NZ)=FHVSETK(K)*WGSHE(NE,K,NB,NZ)
            ENDDO
!            WSSHE(K,NB,NZ)=FHVSETK(K)*WSSHE(K,NB,NZ)
            IF(IHVST(NZ).LE.2.AND.CanPSheathHeight(K,NB,NZ).GT.0.0_r8)THEN
              FHGT=AZMAX1(AMIN1(1.0_r8,(HTNODE(K,NB,NZ) &
                +CanPSheathHeight(K,NB,NZ)-HVST(NZ))/CanPSheathHeight(K,NB,NZ)))
              CanPSheathHeight(K,NB,NZ)=(1._r8-FHGT)*CanPSheathHeight(K,NB,NZ)
            ELSE
              CanPSheathHeight(K,NB,NZ)=FHVSETK(K)*CanPSheathHeight(K,NB,NZ)
            ENDIF
            WGSHGX=WGSHGX+WGSHE(ielmc,K,NB,NZ)
!     IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
!     IF(HTNODE(K,NB,NZ).GT.HVST(NZ)
!    2.OR.IHVST(NZ).EQ.3)THEN
!     IF(isclose(FHVSETK(K),0._r8).AND.K.GT.0)THEN
!     IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
!     VSTG(NB,NZ)=AZMAX1(VSTG(NB,NZ)-1.0)
!     ELSE
!     VSTG(NB,NZ)=AZMAX1(VSTG(NB,NZ)-0.04)
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
!     WTLS,CanPBLeafShethC=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
        CPOOLX=AZMAX1(EPOOL(ielmc,NB,NZ))
        ZPOOLX=AZMAX1(EPOOL(ielmn,NB,NZ))
        PPOOLX=AZMAX1(EPOOL(ielmp,NB,NZ))
        CPOLNX=AZMAX1(EPOLNB(ielmc,NB,NZ))
        ZPOLNX=AZMAX1(EPOLNB(ielmn,NB,NZ))
        PPOLNX=AZMAX1(EPOLNB(ielmp,NB,NZ))
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
            WTLSBX=AZMAX1(CanPBLeafShethC(NB,NZ))
            IF(EPOOL(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
              WHVSCX=AZMAX1(WHVSCP)*WTLSBX/CanopyLeafShethC_pft(NZ)
              CPOOLG=AZMAX1(CPOOLX-WHVSCX)
              ZPOOLG=AZMAX1(ZPOOLX-WHVSCX*ZPOOLX/EPOOL(ielmc,NB,NZ))
              PPOOLG=AZMAX1(PPOOLX-WHVSCX*PPOOLX/EPOOL(ielmc,NB,NZ))
            ELSE
              CPOOLG=0._r8
              ZPOOLG=0._r8
              PPOOLG=0._r8
            ENDIF
            IF(EPOLNB(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
              WHVSNX=AZMAX1(WHVSNP)*WTLSBX/CanopyLeafShethC_pft(NZ)
              CPOLNG=AZMAX1(CPOLNX-WHVSNX)
              ZPOLNG=AZMAX1(ZPOLNX-WHVSNX*ZPOLNX/EPOLNB(ielmc,NB,NZ))
              PPOLNG=AZMAX1(PPOLNX-WHVSNX*PPOLNX/EPOLNB(ielmc,NB,NZ))
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
        EPOOL(ielmc,NB,NZ)=CPOOLG
        EPOOL(ielmn,NB,NZ)=ZPOOLG
        EPOOL(ielmp,NB,NZ)=PPOOLG
        EPOLNB(ielmc,NB,NZ)=CPOLNG
        EPOLNB(ielmn,NB,NZ)=ZPOLNG
        EPOLNB(ielmp,NB,NZ)=PPOLNG
        WTNDBE(ielmc,NB,NZ)=WTNDG
        WTNDBE(ielmn,NB,NZ)=WTNDNG
        WTNDBE(ielmp,NB,NZ)=WTNDPG
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CPOOLG=branch non-structural C mass after harvest
!     WTHTH0E(ielmc),WTHTH0E(ielmn),WTHTH0E(ielmp)=nonstructural C,N,P removed
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
        IF(ICTYP(NZ).EQ.ic4_photo.AND.CPOOLX.GT.ZEROP(NZ))THEN
          FHVST4=CPOOLG/CPOOLX
          D9810: DO K=1,JNODS1
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
          IF(WTSTKE(ielmc,NZ).GT.ZEROL(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/WTSTKE(ielmc,NZ)))
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
        DO NE=1,NumOfPlantChemElements
          WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSH)*WTSTKBE(NE,NB,NZ)
          WTHTX3E(NE)=WTHTX3E(NE)+(FHVSH-FHVSE(ielmc))*WTSTKBE(NE,NB,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
          WTSTKBE(NE,NB,NZ)=FHVSE(ielmc)*WTSTKBE(NE,NB,NZ)
          WTSTXBE(NE,NB,NZ)=FHVSE(ielmc)*WTSTXBE(NE,NB,NZ)
        ENDDO

        CanPBStalkC(NB,NZ)=FHVSE(ielmc)*CanPBStalkC(NB,NZ)
!
!     CUT STALK NODES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTNODX,HTNODE=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN_pft=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTSTK=stalk C mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
        D9820: DO K=JNODS1,0,-1
          IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
            IF(HTNODX(K,NB,NZ).GT.ZERO)THEN
              IF(IHVST(NZ).NE.3)THEN
                FHGTK=AZMAX1(AMIN1(1.0_r8,(HTNODE(K,NB,NZ)-HVST(NZ))/HTNODX(K,NB,NZ)))
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
            IF(WTSTKE(ielmc,NZ).GT.ZEROP(NZ))THEN
              FHVSETS=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/WTSTKE(ielmc,NZ)))
            ELSE
              FHVSETS=1.0_r8
            ENDIF
          ENDIF
          DO NE=1,NumOfPlantChemElements
            WGNODE(NE,K,NB,NZ)=FHVSETS*WGNODE(NE,K,NB,NZ)
          ENDDO
          IF(IHVST(NZ).LE.2.AND.isclose(THIN_pft(NZ),0._r8))THEN
            HTNODX(K,NB,NZ)=FHVSETS*HTNODX(K,NB,NZ)
            HTNODE(K,NB,NZ)=AMIN1(HTNODE(K,NB,NZ),HVST(NZ))
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
          IF(WTSTKBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=FHVSE(ielmc)
            FHVSH=FHVSH
          ELSE
            FHVSE(ielmc)=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(WTRSVE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSE(ielmc)=AZMAX1(AMIN1(1.0_r8,1._r8-WHVRVH/WTRSVE(ielmc,NZ)))
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
        DO NE=1,NumOfPlantChemElements
          WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSH)*WTRSVBE(NE,NB,NZ)
          WTHTX3E(NE)=WTHTX3E(ielmc)+(FHVSH-FHVSE(ielmc))*WTRSVBE(NE,NB,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
          WTRSVBE(NE,NB,NZ)=FHVSE(ielmc)*WTRSVBE(NE,NB,NZ)
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
          IF(WTHSKE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETH=AZMAX1(AMIN1(1.0_r8,1._r8-WHVHSH/WTHSKE(ielmc,NZ)))
            FHVSHH=FHVSETH
          ELSE
            FHVSETH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(WTEARE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETE=AZMAX1(AMIN1(1.0_r8,1._r8-WHVEAH/WTEARE(ielmc,NZ)))
            FHVSHE=FHVSETE
          ELSE
            FHVSETE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(WTGRE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSETG=AZMAX1(AMIN1(1.0_r8,1._r8-WHVGRH/WTGRE(ielmc,NZ)))
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
        DO NE=1,NumOfPlantChemElements
          WTHTH2E(NE)=WTHTH2E(NE)+(1._r8-FHVSHH)*WTHSKBE(NE,NB,NZ)+(1._r8-FHVSHE) &
            *WTEARBE(NE,NB,NZ)+(1._r8-FHVSHG)*WTGRBE(NE,NB,NZ)
          WTHTX2E(NE)=WTHTX2E(NE)+(FHVSHH-FHVSETH)*WTHSKBE(NE,NB,NZ)+(FHVSHE-FHVSETE) &
            *WTEARBE(NE,NB,NZ)+(FHVSHG-FHVSETG)*WTGRBE(NE,NB,NZ)
          WTHTGE(NE)=WTHTGE(NE)+(1._r8-FHVSETG)*WTGRBE(NE,NB,NZ)

!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
          WTHSKBE(NE,NB,NZ)=FHVSETH*WTHSKBE(NE,NB,NZ)
          WTEARBE(NE,NB,NZ)=FHVSETE*WTEARBE(NE,NB,NZ)
          WTGRBE(NE,NB,NZ)=FHVSETG*WTGRBE(NE,NB,NZ)

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
!     CanPBStalkC=stalk sapwood mass
!     PSICanP=canopy water potential
!     CanWatP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        CPOOLK(NB,NZ)=0._r8
        D1325: DO K=1,JNODS1
          CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
            +CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
            +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
        ENDDO D1325

        CanPBLeafShethC(NB,NZ)=AZMAX1(WTLFBE(ielmc,NB,NZ) &
          +WTSHEBE(ielmc,NB,NZ))

        WTSHTBE(ielmc,NB,NZ)=WTSHTBE(ielmc,NB,NZ)+CPOOLK(NB,NZ)
        DO NE=1,NumOfPlantChemElements
          WTSHTBE(NE,NB,NZ)=AZMAX1(WTSHTBE(ielmc,NB,NZ)+WTLFBE(NE,NB,NZ) &
            +WTSHEBE(NE,NB,NZ)+WTSTKBE(NE,NB,NZ)+WTRSVBE(NE,NB,NZ) &
            +WTHSKBE(NE,NB,NZ)+WTEARBE(NE,NB,NZ)+WTGRBE(NE,NB,NZ) &
            +EPOOL(NE,NB,NZ))
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
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     ZC=canopy height
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     IDAY(1,=emergence date
!     GROUP=node number required for floral initiation
!     PSTGI=node number at floral initiation
!     PSTGF=node number at flowering
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     IFLGA=flag for initializing leafout
!
        IF((IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1) &
          .AND.(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6) &
          .AND.CanopyHeight(NZ).GT.HVST(NZ))THEN
          IF((IWTYP(NZ).NE.0.AND.VRNF(NB,NZ) &
            .LE.FVRN(IWTYP(NZ))*VRNX(NB,NZ)) &
            .OR.(IWTYP(NZ).EQ.0.AND.IDAY(1,NB,NZ).NE.0))THEN
            GROUP(NB,NZ)=GROUPI(NZ)
            PSTGI(NB,NZ)=PSTG(NB,NZ)
            PSTGF(NB,NZ)=0._r8
            VSTGX(NB,NZ)=0._r8
            TGSTGI(NB,NZ)=0._r8
            TGSTGF(NB,NZ)=0._r8
            FLG4(NB,NZ)=0._r8
            IDAY(1,NB,NZ)=I
            D3005: DO M=2,10
              IDAY(M,NB,NZ)=0
            ENDDO D3005
            IFLGA(NB,NZ)=0
            IF(NB.EQ.NB1(NZ))THEN
              D3010: DO NBX=1,NumOfBranches_pft(NZ)
                IF(NBX.NE.NB1(NZ))THEN
                  GROUP(NBX,NZ)=GROUPI(NZ)
                  PSTGI(NBX,NZ)=PSTG(NBX,NZ)
                  PSTGF(NBX,NZ)=0._r8
                  VSTGX(NBX,NZ)=0._r8
                  TGSTGI(NBX,NZ)=0._r8
                  TGSTGF(NBX,NZ)=0._r8
                  FLG4(NBX,NZ)=0._r8
                  IDAY(1,NBX,NZ)=I
                  D3015: DO M=2,10
                    IDAY(M,NBX,NZ)=0
                  ENDDO D3015
                  IFLGA(NBX,NZ)=0
                ENDIF
              ENDDO D3010
            ENDIF
          ENDIF
        ENDIF
!
!     DEATH OF BRANCH IF KILLING HARVEST ENTERED IN 'READQ'
!
!     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!     WTLS=total PFT leaf+petiole C mass
!     WTSTK=total PFT stalk C mass
!     WVSTK=total PFT sapwood C mass
!     CanopyBranchStemApft_lyr=total PFT stalk surface area
!
        IF(JHVST(NZ).NE.ihv_noaction)then
          IDTHB(NB,NZ)=ibrdead
        endif
        IF(pftPlantPopulation(NZ).LE.0.0)then
          IDTHB(NB,NZ)=ibrdead
        endif
      ENDDO D9835
      CanopyLeafShethC_pft(NZ)=0._r8
      WTSTKE(ielmc,NZ)=0._r8
      CanPStalkC(NZ)=0._r8
      CanPSA(NZ)=0._r8
      D9840: DO NB=1,NumOfBranches_pft(NZ)
        CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+CanPBLeafShethC(NB,NZ)
        WTSTKE(ielmc,NZ)=WTSTKE(ielmc,NZ)+WTSTKBE(ielmc,NB,NZ)
        CanPStalkC(NZ)=CanPStalkC(NZ)+CanPBStalkC(NB,NZ)
        D9830: DO L=1,NumOfCanopyLayers1
          CanPSA(NZ)=CanPSA(NZ)+CanopyBranchStemApft_lyr(L,NB,NZ)
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
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CO2NetFix_pft=PFT net CO2 fixation
!     TNBP=total net biome productivity
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
              FFIRE(1:NumOfPlantChemElements)=0._r8
            ELSE
              IF(THETW(L).GT.FVLWB.OR.CORGC(L).LE.FORGC.OR.ITILL.NE.22)THEN
                XHVST(ielmc)=1.0_r8
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE(1:NumOfPlantChemElements)=0._r8
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
              DO NE=1,NumOfPlantChemElements
                FHVSE(NE)=XHVST1(NE)*CFOPE(NE,instruct,M,NZ)*EPOOLR(NE,N,L,NZ)
                ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+(1._r8-FFIRE(NE))*FHVSE(NE)
              ENDDO
              VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
              VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
              VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667
              VNH3F(NZ)=VNH3F(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
              VN2OF(NZ)=VN2OF(NZ)-0.0_r8
              VPO4F(NZ)=VPO4F(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
              CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
              TNBP=TNBP-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
              DO NR=1,NRT(NZ)
                DO NE=1,NumOfPlantChemElements
                  FHVSE(NE)=XHVST1(NE)*CFOPE(NE,icwood,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)
                  ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+(1._r8-FFIRE(NE))*FHVSE(NE)
                ENDDO
                VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667
                VNH3F(NZ)=VNH3F(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
                VN2OF(NZ)=VN2OF(NZ)-0.0
                VPO4F(NZ)=VPO4F(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
                CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                TNBP=TNBP-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)

                DO NE=1,NumOfPlantChemElements
                  FHVSE(NE)=XHVST1(NE)*CFOPE(NE,iroot,M,NZ)*(WTRT1E(NE,N,L,NR,NZ) &
                    +WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
                  ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ) &
                    +(1._r8-FFIRE(NE))*FHVSE(NE)
                ENDDO
                VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)*2.667_r8
                VNH3F(NZ)=VNH3F(NZ)-FFIRE(ielmn)*FHVSE(ielmn)
                VN2OF(NZ)=VN2OF(NZ)-0.0_r8
                VPO4F(NZ)=VPO4F(NZ)-FFIRE(ielmp)*FHVSE(ielmp)
                CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-(1._r8-FCH4F)*FFIRE(ielmc)*FHVSE(ielmc)
                TNBP=TNBP-FCH4F*FFIRE(ielmc)*FHVSE(ielmc)
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
!     WTRTL,PopPlantRootC_vr=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,SecndRootXNumL=number of primary,secondary root axes
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D3960: DO NR=1,NRT(NZ)
              DO NE=1,NumOfPlantChemElements
                WTRT1E(NE,N,L,NR,NZ)=WTRT1E(NE,N,L,NR,NZ)*XHVST(NE)
                WTRT2E(NE,N,L,NR,NZ)=WTRT2E(NE,N,L,NR,NZ)*XHVST(NE)
                RTWT1E(NE,N,NR,NZ)=RTWT1E(NE,N,NR,NZ)*XHVST(NE)
              ENDDO
              PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)*XHVST(ielmc)
              SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)*XHVST(ielmc)
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST(ielmc)
            ENDDO D3960
            DO NE=1,NumOfPlantChemElements
              EPOOLR(NE,N,L,NZ)=EPOOLR(NE,N,L,NZ)*XHVST(NE)
            ENDDO
            WTRTL(N,L,NZ)=WTRTL(N,L,NZ)*XHVST(ielmc)
            PopPlantRootC_vr(N,L,NZ)=PopPlantRootC_vr(N,L,NZ)*XHVST(ielmc)
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
            IF(INTYP(NZ).NE.0.AND.N.EQ.1)THEN
              DO NE=1,NumOfPlantChemElements
                D3395: DO M=1,jsken
                  ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+XHVST1(NE) &
                    *(CFOPE(NE,iroot,M,NZ)*WTNDLE(NE,L,NZ) &
                    +CFOPE(NE,instruct,M,NZ)*EPOOLN(NE,L,NZ))
                ENDDO D3395
                WTNDLE(NE,L,NZ)=WTNDLE(NE,L,NZ)*XHVST(NE)
                EPOOLN(NE,L,NZ)=EPOOLN(NE,L,NZ)*XHVST(NE)
              ENDDO
            ENDIF
          ENDDO D3980
        ENDDO D3985
!
!     STORAGE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        IF(ISTYP(NZ).NE.iplt_annual)THEN
          DO NE=1,NumOfPlantChemElements
            D3400: DO M=1,jsken
              ESNC(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ)=ESNC(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ) &
                +(XHVST1(NE)*CFOPE(NE,instruct,M,NZ)*WTRVE(NE,NZ))*FWOODE(NE,k_woody_litr)

              ESNC(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ)=ESNC(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ) &
                +(XHVST1(NE)*CFOPE(NE,instruct,M,NZ)*WTRVE(NE,NZ))*FWOODE(NE,k_fine_litr)
            ENDDO D3400
            WTRVE(NE,NZ)=WTRVE(NE,NZ)*XHVST(NE)
          ENDDO
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

end module PlantDisturbsMod
