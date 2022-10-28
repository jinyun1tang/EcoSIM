module PlantDisturbsMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb,AZMAX1
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = __FILE__
! end_include_section

! disturbance variables
  real(r8) :: WTHTH0E(npelms)
  real(r8) :: WTHTH1E(npelms)
  real(r8) :: WTHTH2E(npelms)
  real(r8) :: WTHTH3E(npelms)
  real(r8) :: WTHTH4E(npelms)
  real(r8) :: WTHTR1E(npelms)
  real(r8) :: WTHTR2E(npelms)
  real(r8) :: WTHTR3E(npelms)
  real(r8) :: WTHTR4E(npelms)
  real(r8) :: WTHTX0E(npelms)
  real(r8) :: WTHTX1E(npelms)
  real(r8) :: WTHTX2E(npelms)
  real(r8) :: WTHTX3E(npelms)
  real(r8) :: WTHTX4E(npelms)
  real(r8) :: WTHTGE(npelms)
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
  real(r8), intent(inout) :: CPOOLK(JC1,JP1)
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
  real(r8), INTENT(INOUT) :: CPOOLK(JC1,JP1)
  real(r8) :: FHVST
  real(r8) :: FHVSH
  real(r8) :: WHVSTD
  integer :: M

!     begin_execution
  associate(                            &
    EHVST    =>  plt_distb%EHVST  , &
    HVST     =>  plt_distb%HVST   , &
    IHVST    =>  plt_distb%IHVST  , &
    THIN     =>  plt_distb%THIN   , &
    PP       =>  plt_site%PP      , &
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
!     THIN=thinning:fraction of population removed
!     FHVST=fraction of standing dead mass not harvested
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
      IF(test_aeqb(THIN(NZ),0._r8))THEN
        FHVST=AZMAX1(1._r8-EHVST(1,4,NZ))
        FHVSH=FHVST
      ELSE
        FHVST=AZMAX1(1._r8-THIN(NZ))
        IF(IHVST(NZ).EQ.0)THEN
          FHVSH=AZMAX1(1._r8-EHVST(1,4,NZ)*THIN(NZ))
        ELSE
          FHVSH=FHVST
        ENDIF
      ENDIF
    ELSEIF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
      IF(WTSTGE(ielmc,NZ).GT.ZEROP(NZ))THEN
        WHVSTD=HVST(NZ)*THIN(NZ)*0.45_r8/24.0_r8*AREA3(NU)*EHVST(1,4,NZ)
        FHVST=AZMAX1(1._r8-WHVSTD/WTSTGE(ielmc,NZ))
        FHVSH=FHVST
      ELSE
        FHVST=1.0_r8
        FHVSH=1.0_r8
      ENDIF
    ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
    ENDIF
    D6475: DO M=1,jsken
      WTHTH4E(ielmc)=WTHTH4E(ielmc)+(1._r8-FHVSH)*WTSTDE(M,ielmc,NZ)
      WTHTH4E(ielmn)=WTHTH4E(ielmn)+(1._r8-FHVSH)*WTSTDE(M,ielmn,NZ)
      WTHTH4E(ielmp)=WTHTH4E(ielmp)+(1._r8-FHVSH)*WTSTDE(M,ielmp,NZ)
      WTHTX4E(ielmc)=WTHTX4E(ielmc)+(FHVSH-FHVST)*WTSTDE(M,ielmc,NZ)
      WTHTX4E(ielmn)=WTHTX4E(ielmn)+(FHVSH-FHVST)*WTSTDE(M,ielmn,NZ)
      WTHTX4E(ielmp)=WTHTX4E(ielmp)+(FHVSH-FHVST)*WTSTDE(M,ielmp,NZ)
      WTSTDE(M,ielmc,NZ)=FHVST*WTSTDE(M,ielmc,NZ)
      WTSTDE(M,ielmn,NZ)=FHVST*WTSTDE(M,ielmn,NZ)
      WTSTDE(M,ielmp,NZ)=FHVST*WTSTDE(M,ielmp,NZ)
    ENDDO D6475
!
    call PlantDisturbance(I,J,NZ)

    ZEROP(NZ)=ZERO*PP(NZ)
    ZEROQ(NZ)=ZERO*PP(NZ)/AREA3(NU)
    ZEROL(NZ)=ZERO*PP(NZ)*1.0E+06
  ENDIF
  end associate
  end subroutine RemoveBiomassByDisturbance


!------------------------------------------------------------------------------------------
  subroutine PlantDisturbance(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  real(r8) :: WTHNL0,WTHPL0
  real(r8) :: WTHNL1,WTHPL1
  real(r8) :: WTHNL2,WTHPL2
  real(r8) :: WTHNL3,WTHPL3
  real(r8) :: WTHNL4,WTHPL4
  real(r8) :: WTHTR0,WTHNR0,WTHPR0
  real(r8) :: WTHNRT,WTHPRT,WTHTRT
  real(r8) :: WTHTXT,WTHNXT,WTHPXT

  WTHNL0=0._r8;WTHPL0=0._r8
  WTHNL1=0._r8;WTHPL1=0._r8
  WTHNL2=0._r8;WTHPL2=0._r8
  WTHNL3=0._r8;WTHPL3=0._r8
  WTHNL4=0._r8;WTHPL4=0._r8
  WTHTR0=0._r8;WTHNR0=0._r8;WTHPR0=0._r8

  call ApplyDisturbanceBiomRemoval(I,J,NZ,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call TotalBiomRemovalByDisturbance(I,J,NZ,WTHTR0,WTHNR0,WTHPR0,WTHPXT,WTHTRT,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)
!
!     ABOVE-GROUND LITTERFALL FROM HARVESTING
!
  call LiterfallByDisturbance(I,J,NZ,WTHPXT,WTHTRT,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)
  end subroutine PlantDisturbance
!------------------------------------------------------------------------------------------

  subroutine LiterfallByDisturbance(I,J,NZ,WTHPXT,WTHTRT,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)

  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: WTHPXT,WTHTRT
  real(r8), intent(in) :: WTHTR0,WTHNR0,WTHPR0,WTHNL0,WTHPL0
  real(r8), intent(in) :: WTHNL1,WTHPL1
  real(r8), intent(in) :: WTHNL2,WTHPL2
  real(r8), intent(in) :: WTHNL3,WTHPL3
  real(r8), intent(in) :: WTHNL4,WTHPL4
  real(r8), intent(in) :: WTHNRT,WTHPRT,WTHTXT,WTHNXT
  integer :: m
!     begin_execution
  associate(                        &
    instruct =>  pltpar%instruct  , &
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
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
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
        ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
          +CFOPE(instruct,M,ielmc,NZ)*(WTHTR0+WTHTX0E(ielmc)) &
          +CFOPE(ifoliar,M,ielmc,NZ)*(WTHTR1E(ielmc)+WTHTX1E(ielmc)) &
          +CFOPE(infoliar,M,ielmc,NZ)*(WTHTR2E(ielmc)+WTHTX2E(ielmc))
        ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
          +CFOPE(instruct,M,ielmn,NZ)*(WTHNR0+WTHTX0E(ielmn)) &
          +CFOPE(ifoliar,M,ielmn,NZ)*(WTHTR1E(ielmn)+WTHTX1E(ielmn)) &
          +CFOPE(infoliar,M,ielmn,NZ)*(WTHTR2E(ielmn)+WTHTX2E(ielmn))
        ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
          +CFOPE(instruct,M,ielmp,NZ)*(WTHPR0+WTHTX0E(ielmp)) &
          +CFOPE(ifoliar,M,ielmp,NZ)*(WTHTR1E(ielmp)+WTHTX1E(ielmp)) &
          +CFOPE(infoliar,M,ielmp,NZ)*(WTHTR2E(ielmp)+WTHTX2E(ielmp))
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
            +CFOPE(istalk,M,ielmc,NZ)*(WTHTR3E(ielmc)+WTHTX3E(ielmc)+WTHTR4E(ielmc)+WTHTX4E(ielmc))
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
            +CFOPE(istalk,M,ielmn,NZ)*(WTHTR3E(ielmn)+WTHTX3E(ielmn)+WTHTR4E(ielmn)+WTHTX4E(ielmn))
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
            +CFOPE(istalk,M,ielmp,NZ)*(WTHTR3E(ielmp)+WTHTX3E(ielmp)+WTHTR4E(ielmp)+WTHTX4E(ielmp))
        ELSE
          WTSTDE(M,ielmc,NZ)=WTSTDE(M,ielmc,NZ) &
            +CFOPE(icwood,M,ielmc,NZ)*(WTHTX3E(ielmc)+WTHTX4E(ielmc))
          WTSTDE(M,ielmn,NZ)=WTSTDE(M,ielmn,NZ) &
            +CFOPE(icwood,M,ielmn,NZ)*(WTHTX3E(ielmn)+WTHTX4E(ielmn))
          WTSTDE(M,ielmp,NZ)=WTSTDE(M,ielmp,NZ) &
            +CFOPE(icwood,M,ielmp,NZ)*(WTHTX3E(ielmp)+WTHTX4E(ielmp))
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ) &
            +CFOPE(icwood,M,ielmc,NZ)*(WTHTR3E(ielmc)+WTHTR4E(ielmc))*FWOODE(ielmc,0)
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ) &
            +CFOPE(icwood,M,ielmn,NZ)*(WTHTR3E(ielmn)+WTHTR4E(ielmn))*FWOODE(ielmn,0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ) &
            +CFOPE(icwood,M,ielmp,NZ)*(WTHTR3E(ielmp)+WTHTR4E(ielmp))*FWOODE(ielmp,0)
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
            +CFOPE(icwood,M,ielmc,NZ)*(WTHTR3E(ielmc)+WTHTR4E(ielmc))*FWOODE(ielmc,1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
            +CFOPE(icwood,M,ielmn,NZ)*(WTHTR3E(ielmn)+WTHTR4E(ielmn))*FWOODE(ielmn,1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
            +CFOPE(icwood,M,ielmp,NZ)*(WTHTR3E(ielmp)+WTHTR4E(ielmp))*FWOODE(ielmp,0)
        ENDIF
      ENDDO D6375
!
!     ABOVE-GROUND LITTERFALL FROM FIRE
!
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     WTHTX1E(ielmc),WTHTX1E(ielmn),WTHTX1E(ielmp)=harvested leaf C,N,P to litter
!     WTHTX2E(ielmc),WTHTX2E(ielmn),WTHTX2E(ielmp)=harvested petiole C,N,P to litter
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested woody C,N,P to litter
!     WTHTX4E(ielmc),WTHTX4E(ielmn),WTHTX4E(ielmp)=harvested standing dead C,N,P to litter
!     WTHNL0,WTHPL0=nonstructural N,P to litter
!     WTHNL1,WTHPL1=leaf N,P to litter
!     WTHNL2,WTHPL2=fine,non-leaf N,P to litter
!     WTHNL3,WTHPL3=woody N,P to litter
!     WTHNL4,WTHPL4=standing dead N,P to litter
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    ELSE
      D6485: DO M=1,jsken
        ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
          +CFOPE(instruct,M,ielmc,NZ)*(WTHTR0+WTHTX0E(ielmc)) &
          +CFOPE(ifoliar,M,ielmc,NZ)*(WTHTR1E(ielmc)+WTHTX1E(ielmc)) &
          +CFOPE(infoliar,M,ielmc,NZ)*(WTHTR2E(ielmc)+WTHTX2E(ielmc))
        ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
          +CFOPE(instruct,M,ielmn,NZ)*WTHNL0 &
          +CFOPE(ifoliar,M,ielmn,NZ)*WTHNL1 &
          +CFOPE(infoliar,M,ielmn,NZ)*WTHNL2
        ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
          +CFOPE(instruct,M,ielmp,NZ)*WTHPL0 &
          +CFOPE(ifoliar,M,ielmp,NZ)*WTHPL1 &
          +CFOPE(infoliar,M,ielmp,NZ)*WTHPL2
        ESNC(4,ielmc,1,0,NZ)=ESNC(4,ielmc,1,0,NZ) &
          +CFOPE(instruct,M,ielmn,NZ)*(WTHNR0+WTHTX0E(ielmn)-WTHNL0) &
          +CFOPE(ifoliar,M,ielmn,NZ)*(WTHTR1E(ielmn)+WTHTX1E(ielmn)-WTHNL1) &
          +CFOPE(infoliar,M,ielmn,NZ)*(WTHTR2E(ielmn)+WTHTX2E(ielmn)-WTHNL2)
        ESNC(4,ielmp,1,0,NZ)=ESNC(4,ielmp,1,0,NZ) &
          +CFOPE(instruct,M,ielmp,NZ)*(WTHPR0+WTHTX0E(ielmp)-WTHPL0) &
          +CFOPE(ifoliar,M,ielmp,NZ)*(WTHTR1E(ielmp)+WTHTX1E(ielmp)-WTHPL1) &
          +CFOPE(infoliar,M,ielmp,NZ)*(WTHTR2E(ielmp)+WTHTX2E(ielmp)-WTHPL2)
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPE(istalk,M,ielmc,NZ)*(WTHTR3E(ielmc)+WTHTX3E(ielmc)+WTHTR4E(ielmc)+WTHTX4E(ielmc))
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPE(istalk,M,ielmn,NZ)*(WTHNL3+WTHNL4)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPE(istalk,M,ielmp,NZ)*(WTHPL3+WTHPL4)
          ESNC(4,ielmn,1,0,NZ)=ESNC(4,ielmn,1,0,NZ)+CFOPE(istalk,M,ielmn,NZ)*(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHNL3+WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHNL4)
          ESNC(4,ielmp,1,0,NZ)=ESNC(4,ielmp,1,0,NZ)+CFOPE(istalk,M,ielmp,NZ)*(WTHTR3E(ielmp)+WTHTX3E(ielmp)-WTHPL3+WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHPL4)
        ELSE
          WTSTDE(M,ielmc,NZ)=WTSTDE(M,ielmc,NZ)+CFOPE(icwood,M,ielmc,NZ)*(WTHTR3E(ielmc)+WTHTX3E(ielmc))
          WTSTDE(M,ielmn,NZ)=WTSTDE(M,ielmn,NZ)+CFOPE(icwood,M,ielmn,NZ)*WTHNL3
          WTSTDE(M,ielmp,NZ)=WTSTDE(M,ielmp,NZ)+CFOPE(icwood,M,ielmp,NZ)*WTHPL3
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)*CFOPE(istalk,M,ielmc,NZ)*(WTHTR4E(ielmc)+WTHTX4E(ielmc))*FWOODE(ielmc,0)
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+CFOPE(istalk,M,ielmn,NZ)*WTHNL4*FWOODE(ielmn,0)
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+CFOPE(istalk,M,ielmp,NZ)*WTHPL4*FWOODE(ielmp,0)
          ESNC(4,ielmn,0,0,NZ)=ESNC(4,ielmn,0,0,NZ)+CFOPE(icwood,M,ielmn,NZ)*(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHNL3 &
            +WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHNL4)*FWOODE(ielmn,0)
          ESNC(4,ielmp,0,0,NZ)=ESNC(4,ielmp,0,0,NZ)+CFOPE(icwood,M,ielmp,NZ)*(WTHTR3E(ielmp)+WTHTX3E(ielmp)-WTHPL3 &
            +WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHPL4)*FWOODE(ielmp,0)
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPE(istalk,M,ielmc,NZ)*(WTHTR4E(ielmc)+WTHTX4E(ielmc))*FWOODE(ielmc,1)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPE(istalk,M,ielmn,NZ)*WTHNL4*FWOODE(ielmn,1)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPE(istalk,M,ielmp,NZ)*WTHPL4*FWOODE(ielmp,1)
          ESNC(4,ielmn,1,0,NZ)=ESNC(4,ielmn,1,0,NZ)+CFOPE(icwood,M,ielmn,NZ)*(WTHTR3E(ielmn)+WTHTX3E(ielmn)-WTHNL3 &
            +WTHTR4E(ielmn)+WTHTX4E(ielmn)-WTHNL4)*FWOODE(ielmn,1)
          ESNC(4,ielmp,1,0,NZ)=ESNC(4,ielmp,1,0,NZ)+CFOPE(icwood,M,ielmp,NZ)*(WTHTR3E(ielmp)+WTHTX3E(ielmp)-WTHPL3 &
            +WTHTR4E(ielmp)+WTHTX4E(ielmp)-WTHPL4)*FWOODE(ielmp,1)
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
    TESNC(ielmc,NZ)=TESNC(ielmc,NZ)+WTHTRT+WTHTXT
    TESNC(ielmn,NZ)=TESNC(ielmn,NZ)+WTHNRT+WTHNXT
    TESNC(ielmp,NZ)=TESNC(ielmp,NZ)+WTHPRT+WTHPXT
    TESN0(ielmc,NZ)=TESN0(ielmc,NZ)+WTHTRT+WTHTXT
    TESN0(ielmn,NZ)=TESNC(ielmn,NZ)+WTHNRT+WTHNXT
    TESN0(ielmp,NZ)=TESNC(ielmp,NZ)+WTHPRT+WTHPXT
  ENDIF
  end associate
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

  subroutine TotalBiomRemovalByDisturbance(I,J,NZ,WTHTR0,WTHNR0,WTHPR0,WTHPXT,WTHTRT,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: WTHTR0,WTHNR0,WTHPR0
  real(r8), intent(out):: WTHPXT,WTHTRT,WTHNRT,WTHPRT,WTHTXT,WTHNXT
  real(r8) :: WTHTHT,WTHNHT,WTHPHT
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
    CNET     =>  plt_bgcr%CNET    , &
    TNBP     =>  plt_bgcr%TNBP    , &
    TRAU     =>  plt_bgcr%TRAU    , &
    RECO     =>  plt_bgcr%RECO    , &
    TCO2A    =>  plt_bgcr%TCO2A   , &
    TCO2T    =>  plt_bgcr%TCO2T   , &
    WTRVE    =>  plt_biom%WTRVE     &
  )
!
!     WTHTHT,WTHNHT,WTHPHT=total C,N,P removed
!     WTHTRT,WTHNRT,WTHPRT=total C,N,P to litter
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  WTHTHT=WTHTH0E(ielmc)+WTHTH1E(ielmc)+WTHTH2E(ielmc)+WTHTH3E(ielmc)+WTHTH4E(ielmc)
  WTHTRT=WTHTR0+WTHTR1E(ielmc)+WTHTR2E(ielmc)+WTHTR3E(ielmc)+WTHTR4E(ielmc)
  WTHNHT=WTHTH0E(ielmn)+WTHTH1E(ielmn)+WTHTH2E(ielmn)+WTHTH3E(ielmn)+WTHTH4E(ielmn)
  WTHNRT=WTHNR0+WTHTR1E(ielmn)+WTHTR2E(ielmn)+WTHTR3E(ielmn)+WTHTR4E(ielmn)
  WTHPHT=WTHTH0E(ielmp)+WTHTH1E(ielmp)+WTHTH2E(ielmp)+WTHTH3E(ielmp)+WTHTH4E(ielmp)
  WTHPRT=WTHPR0+WTHTR1E(ielmp)+WTHTR2E(ielmp)+WTHTR3E(ielmp)+WTHTR4E(ielmp)
  WTHTXT=WTHTX0E(ielmc)+WTHTX1E(ielmc)+WTHTX2E(ielmc)+WTHTX3E(ielmc)+WTHTX4E(ielmc)
  WTHNXT=WTHTX0E(ielmn)+WTHTX1E(ielmn)+WTHTX2E(ielmn)+WTHTX3E(ielmn)+WTHTX4E(ielmn)
  WTHPXT=WTHTX0E(ielmp)+WTHTX1E(ielmp)+WTHTX2E(ielmp)+WTHTX3E(ielmp)+WTHTX4E(ielmp)

  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      IF(JHVST(NZ).NE.2)THEN
        HVSTE(ielmc,NZ)=HVSTE(ielmc,NZ)+WTHTHT-WTHTRT
        HVSTE(ielmn,NZ)=HVSTE(ielmn,NZ)+WTHNHT-WTHNRT
        HVSTE(ielmp,NZ)=HVSTE(ielmp,NZ)+WTHPHT-WTHPRT
        TNBP=TNBP+WTHTRT-WTHTHT
        XHVSTE(ielmc)=XHVSTE(ielmc)+WTHTHT-WTHTRT
        XHVSTE(ielmn)=XHVSTE(ielmn)+WTHNHT-WTHNRT
        XHVSTE(ielmp)=XHVSTE(ielmp)+WTHPHT-WTHPRT
      ELSE
        WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+WTHTHT-WTHTRT
        WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+WTHNHT-WTHNRT
        WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+WTHPHT-WTHPRT
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CNET=PFT net CO2 fixation
!     TNBP=total net biome productivity
!
    ELSE
      VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*(WTHTHT-WTHTRT)
      VCH4F(NZ)=VCH4F(NZ)-FCH4F*(WTHTHT-WTHTRT)
      VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*(WTHTHT-WTHTRT)*2.667
      VNH3F(NZ)=VNH3F(NZ)-WTHNHT+WTHNRT
      VN2OF(NZ)=VN2OF(NZ)-0.0
      VPO4F(NZ)=VPO4F(NZ)-WTHPHT+WTHPRT
      CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*(WTHTHT-WTHTRT)
      TNBP=TNBP-FCH4F*(WTHTHT-WTHTRT)
    ENDIF
!
!     C,N,P REMOVED FROM GRAZING
!
!     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
!     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
!     GY=growth yield of grazers
!     WTHTHT,WTHNHT,WTHPHT=total C,N,P removed
!     WTHTRT,WTHNRT,WTHPRT=total C,N,P to litter
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
  ELSE
    HVSTE(ielmc,NZ)=HVSTE(ielmc,NZ)+GY*(WTHTHT-WTHTRT)
    HVSTE(ielmn,NZ)=HVSTE(ielmn,NZ)+WTHNHT-WTHNRT
    HVSTE(ielmp,NZ)=HVSTE(ielmp,NZ)+WTHPHT-WTHPRT
    TCO2T(NZ)=TCO2T(NZ)-GZ*(WTHTHT-WTHTRT)
    TCO2A(NZ)=TCO2A(NZ)-GZ*(WTHTHT-WTHTRT)
!     TNBP=TNBP+GY*(WTHTRT-WTHTHT)
!     CNET(NZ)=CNET(NZ)+GZ*(WTHTRT-WTHTHT)
    XHVSTE(ielmc)=XHVSTE(ielmc)+GY*(WTHTHT-WTHTRT)
    XHVSTE(ielmn)=XHVSTE(ielmn)+WTHNHT-WTHNRT
    XHVSTE(ielmp)=XHVSTE(ielmp)+WTHPHT-WTHPRT
    RECO=RECO-GZ*(WTHTHT-WTHTRT)
    TRAU=TRAU-GZ*(WTHTHT-WTHTRT)
  ENDIF
  end associate
  end subroutine TotalBiomRemovalByDisturbance
!------------------------------------------------------------------------------------------

  subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: WTHTR0,WTHNR0,WTHPR0
  real(r8), intent(out) :: WTHNL0,WTHPL0
  real(r8), intent(out) :: WTHNL1,WTHPL1
  real(r8), intent(out) :: WTHNL2,WTHPL2
  real(r8), intent(out) :: WTHNL3,WTHPL3
  real(r8), intent(out) :: WTHNL4,WTHPL4

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
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  IF(IHVST(NZ).EQ.0)THEN
    WTHTR0=WTHTH0E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHTH0E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHPR0=WTHTH0E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmc)=WTHTH1E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmn)=WTHTH1E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmp)=WTHTH1E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR2E(ielmc)=WTHTH2E(ielmc)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmn)=WTHTH2E(ielmn)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmp)=WTHTH2E(ielmp)*(1._r8-EHVST(2,2,NZ))
    WTHTR3E(ielmc)=WTHTH3E(ielmc)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmn)=WTHTH3E(ielmn)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmp)=WTHTH3E(ielmp)*(1._r8-EHVST(2,3,NZ))
    WTHTR4E(ielmc)=WTHTH4E(ielmc)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmn)=WTHTH4E(ielmn)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmp)=WTHTH4E(ielmp)*(1._r8-EHVST(2,4,NZ))
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ).EQ.1)THEN
    WTHTR0=WTHTH0E(ielmc)
    WTHNR0=WTHTH0E(ielmn)
    WTHPR0=WTHTH0E(ielmp)
    WTHTR1E(ielmc)=WTHTH1E(ielmc)
    WTHTR1E(ielmn)=WTHTH1E(ielmn)
    WTHTR1E(ielmp)=WTHTH1E(ielmp)
    WTHTR2E(ielmc)=WTHTH2E(ielmc)-WTHTGE(ielmc)*EHVST(2,2,NZ)
    WTHTR2E(ielmn)=WTHTH2E(ielmn)-WTHTGE(ielmn)*EHVST(2,2,NZ)
    WTHTR2E(ielmp)=WTHTH2E(ielmp)-WTHTGE(ielmp)*EHVST(2,2,NZ)
    WTHTR3E(ielmc)=WTHTH3E(ielmc)
    WTHTR3E(ielmn)=WTHTH3E(ielmn)
    WTHTR3E(ielmp)=WTHTH3E(ielmp)
    WTHTR4E(ielmc)=WTHTH4E(ielmc)
    WTHTR4E(ielmn)=WTHTH4E(ielmn)
    WTHTR4E(ielmp)=WTHTH4E(ielmp)
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ).EQ.2)THEN
    WTHTR0=WTHTH0E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHTH0E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHPR0=WTHTH0E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmc)=WTHTH1E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmn)=WTHTH1E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmp)=WTHTH1E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR2E(ielmc)=WTHTH2E(ielmc)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmn)=WTHTH2E(ielmn)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmp)=WTHTH2E(ielmp)*(1._r8-EHVST(2,2,NZ))
    WTHTR3E(ielmc)=WTHTH3E(ielmc)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmn)=WTHTH3E(ielmn)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmp)=WTHTH3E(ielmp)*(1._r8-EHVST(2,3,NZ))
    WTHTR4E(ielmc)=WTHTH4E(ielmc)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmn)=WTHTH4E(ielmn)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmp)=WTHTH4E(ielmp)*(1._r8-EHVST(2,4,NZ))
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(IHVST(NZ).EQ.3)THEN
    WTHTR0=WTHTH0E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHTH0E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHPR0=WTHTH0E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmc)=WTHTH1E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmn)=WTHTH1E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmp)=WTHTH1E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR2E(ielmc)=WTHTH2E(ielmc)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmn)=WTHTH2E(ielmn)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmp)=WTHTH2E(ielmp)*(1._r8-EHVST(2,2,NZ))
    WTHTR3E(ielmc)=WTHTH3E(ielmc)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmn)=WTHTH3E(ielmn)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmp)=WTHTH3E(ielmp)*(1._r8-EHVST(2,3,NZ))
    WTHTR4E(ielmc)=WTHTH4E(ielmc)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmn)=WTHTH4E(ielmn)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmp)=WTHTH4E(ielmp)*(1._r8-EHVST(2,4,NZ))
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
    WTHTR0=WTHTH0E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHTH0E(ielmn)*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHPR0=WTHTH0E(ielmp)*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHTR1E(ielmc)=WTHTH1E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmn)=WTHTH1E(ielmn)*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHTR1E(ielmp)=WTHTH1E(ielmp)*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHTR2E(ielmc)=WTHTH2E(ielmc)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmn)=WTHTH2E(ielmn)*(1._r8-EHVST(2,2,NZ)*0.5)
    WTHTR2E(ielmp)=WTHTH2E(ielmp)*(1._r8-EHVST(2,2,NZ)*0.5)
    WTHTR3E(ielmc)=WTHTH3E(ielmc)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmn)=WTHTH3E(ielmn)*(1._r8-EHVST(2,3,NZ)*0.5)
    WTHTR3E(ielmp)=WTHTH3E(ielmp)*(1._r8-EHVST(2,3,NZ)*0.5)
    WTHTR4E(ielmc)=WTHTH4E(ielmc)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmn)=WTHTH4E(ielmn)*(1._r8-EHVST(2,4,NZ)*0.5)
    WTHTR4E(ielmp)=WTHTH4E(ielmp)*(1._r8-EHVST(2,4,NZ)*0.5)
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
    FERT(17)=FERT(17)+(WTHTR0+WTHTR1E(ielmc)+WTHTR2E(ielmc)+WTHTR3E(ielmc)+WTHTR4E(ielmc))/AREA3(NU)
    FERT(18)=FERT(18)+(WTHNR0+WTHTR1E(ielmn)+WTHTR2E(ielmn)+WTHTR3E(ielmn)+WTHTR4E(ielmn))/AREA3(NU)*0.5_r8
    FERT(3)=FERT(3)+(WTHNR0+WTHTR1E(ielmn)+WTHTR2E(ielmn)+WTHTR3E(ielmn)+WTHTR4E(ielmn))/AREA3(NU)*0.5_r8
    FERT(19)=FERT(19)+(WTHPR0+WTHTR1E(ielmp)+WTHTR2E(ielmp)+WTHTR3E(ielmp)+WTHTR4E(ielmp))/AREA3(NU)
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
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1E(ielmc),WTHTR1E(ielmn),WTHTR1E(ielmp)=leaf C,N,P to litter
!     WTHTR2E(ielmc),WTHTR2E(ielmn),WTHTR2E(ielmp)=fine,non-leaf C,N,P to litter
!     WTHTR3E(ielmc),WTHTR3E(ielmn),WTHTR3E(ielmp)=woody C,N,P to litter
!     WTHTR4E(ielmc),WTHTR4E(ielmn),WTHTR4E(ielmp)=standing dead C,N,P to litter
!     WTHTL0,WTHNL0,WTHPL0=nonstructural C,N,P removed from ecosystem
!     WTHTL1,WTHNL1,WTHPL1=leaf C,N,P removed from ecosystem
!     WTHTL2,WTHNL2,WTHPL2=fine,non-leaf C,N,P removed from ecosystem
!     WTHTL3,WTHNL3,WTHPL3=woody C,N,P removed from ecosystem
!     WTHTL4,WTHNL4,WTHPL4=standing dead C,N,P removed from ecosystem
!
  ELSEIF(IHVST(NZ).EQ.5)THEN
    WTHTR0=WTHTH0E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHTH0E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,1,NZ))
    WTHPR0=WTHTH0E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,1,NZ))
    WTHNL0=WTHTH0E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHPL0=WTHTH0E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmc)=WTHTH1E(ielmc)*(1._r8-EHVST(2,1,NZ))
    WTHTR1E(ielmn)=WTHTH1E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,1,NZ))
    WTHTR1E(ielmp)=WTHTH1E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,1,NZ))
    WTHNL1=WTHTH1E(ielmn)*(1._r8-EHVST(2,1,NZ))
    WTHPL1=WTHTH1E(ielmp)*(1._r8-EHVST(2,1,NZ))
    WTHTR2E(ielmc)=WTHTH2E(ielmc)*(1._r8-EHVST(2,2,NZ))
    WTHTR2E(ielmn)=WTHTH2E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,2,NZ))
    WTHTR2E(ielmp)=WTHTH2E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,2,NZ))
    WTHNL2=WTHTH2E(ielmn)*(1._r8-EHVST(2,2,NZ))
    WTHPL2=WTHTH2E(ielmp)*(1._r8-EHVST(2,2,NZ))
    WTHTR3E(ielmc)=WTHTH3E(ielmc)*(1._r8-EHVST(2,3,NZ))
    WTHTR3E(ielmn)=WTHTH3E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,3,NZ))
    WTHTR3E(ielmp)=WTHTH3E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,3,NZ))
    WTHNL3=WTHTH3E(ielmn)*(1._r8-EHVST(2,3,NZ))
    WTHPL3=WTHTH3E(ielmp)*(1._r8-EHVST(2,3,NZ))
    WTHTR4E(ielmc)=WTHTH4E(ielmc)*(1._r8-EHVST(2,4,NZ))
    WTHTR4E(ielmn)=WTHTH4E(ielmn)*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,4,NZ))
    WTHTR4E(ielmp)=WTHTH4E(ielmp)*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,4,NZ))
    WTHNL4=WTHTH4E(ielmn)*(1._r8-EHVST(2,4,NZ))
    WTHPL4=WTHTH4E(ielmp)*(1._r8-EHVST(2,4,NZ))
  ENDIF
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByTillage(I,J,NZ,CPOOLK)
  use SurfLitterDataType, only : XCORP
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(inout) :: CPOOLK(10,JP1)
  integer :: L,K,M,N,NR,NB,NE
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
    CH4P     =>  plt_rbgc%CH4P       , &
    CH4A     =>  plt_rbgc%CH4A       , &
    H2GP     =>  plt_rbgc%H2GP       , &
    H2GA     =>  plt_rbgc%H2GA       , &
    OXYP     =>  plt_rbgc%OXYP       , &
    CO2P     =>  plt_rbgc%CO2P       , &
    OXYA     =>  plt_rbgc%OXYA       , &
    CO2A     =>  plt_rbgc%CO2A       , &
    UVOLO    =>  plt_ew%UVOLO        , &
    VOLWP    =>  plt_ew%VOLWP        , &
    VHCPC    =>  plt_ew%VHCPC        , &
    PSILT    =>  plt_ew%PSILT        , &
    PPX      =>  plt_site%PPX        , &
    EPOOLR   =>  plt_biom%EPOOLR     , &
    WSRTL    =>  plt_biom%WSRTL      , &
    WTRTD    =>  plt_biom%WTRTD      , &
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
    WTLSB    =>  plt_biom%WTLSB      , &
    WTSTXBE  =>  plt_biom%WTSTXBE    , &
    WVSTKB   =>  plt_biom%WVSTKB     , &
    WSSHE    =>  plt_biom%WSSHE      , &
    WGSHE    =>  plt_biom%WGSHE      , &
    WTRT1E   =>  plt_biom%WTRT1E     , &
    WSLF     =>  plt_biom%WSLF       , &
    WGNODE   =>  plt_biom%WGNODE     , &
    WVSTK    =>  plt_biom%WVSTK      , &
    RTWT1E   =>  plt_biom%RTWT1E     , &
    WTRVE    =>  plt_biom%WTRVE      , &
    WTLS     =>  plt_biom%WTLS       , &
    WTRT2E   =>  plt_biom%WTRT2E     , &
    EPOOLN   =>  plt_biom%EPOOLN     , &
    WTNDLE   =>  plt_biom%WTNDLE     , &
    GRWTB    =>  plt_allom%GRWTB     , &
    FWOODE   =>  plt_allom%FWOODE    , &
    FWODLP   =>  plt_allom%FWODLP    , &
    FWODSN   =>  plt_allom%FWODSN    , &
    FWODB    =>  plt_allom%FWODB     , &
    FWODRP   =>  plt_allom%FWODRP    , &
    FWODSP   =>  plt_allom%FWODSP    , &
    FWODLN   =>  plt_allom%FWODLN    , &
    FWODR    =>  plt_allom%FWODR     , &
    FWODRN   =>  plt_allom%FWODRN    , &
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
    PP       =>  plt_site%PP         , &
    IYRC     =>  plt_site%IYRC       , &
    ZNOON    =>  plt_site%ZNOON      , &
    VOLWOU   =>  plt_site%VOLWOU     , &
    NU       =>  plt_site%NU         , &
    instruct =>  pltpar%instruct  , &
    ifoliar  =>  pltpar%ifoliar   , &
    istalk   =>  pltpar%istalk    , &
    iroot    =>  pltpar%iroot     , &
    infoliar =>  pltpar%infoliar  , &
    icwood   =>  pltpar%icwood    , &
    ESNC     =>  plt_bgcr%ESNC       , &
    RH2GZ    =>  plt_bgcr%RH2GZ      , &
    RNH3Z    =>  plt_bgcr%RNH3Z      , &
    RN2OZ    =>  plt_bgcr%RN2OZ      , &
    RCO2Z    =>  plt_bgcr%RCO2Z      , &
    RCH4Z    =>  plt_bgcr%RCH4Z      , &
    ROXYZ    =>  plt_bgcr%ROXYZ      , &
    ZH3P     =>  plt_rbgc%ZH3P       , &
    Z2OP     =>  plt_rbgc%Z2OP       , &
    ZH3A     =>  plt_rbgc%ZH3A       , &
    RCO2A    =>  plt_rbgc%RCO2A      , &
    Z2OA     =>  plt_rbgc%Z2OA       , &
    RCO2M    =>  plt_rbgc%RCO2M      , &
    RCO2N    =>  plt_rbgc%RCO2N      , &
    FRADP    =>  plt_rad%FRADP       , &
    RTLG1    =>  plt_morph%RTLG1     , &
    RTVLW    =>  plt_morph%RTVLW     , &
    RTARP    =>  plt_morph%RTARP     , &
    RTVLP    =>  plt_morph%RTVLP     , &
    RTDNP    =>  plt_morph%RTDNP     , &
    RTN1     =>  plt_morph%RTN1      , &
    INTYP    =>  plt_morph%INTYP     , &
    RTLGP    =>  plt_morph%RTLGP     , &
    RTN2     =>  plt_morph%RTN2      , &
    RTLG2    =>  plt_morph%RTLG2     , &
    RTNL     =>  plt_morph%RTNL      , &
    NG       =>  plt_morph%NG        , &
    MY       =>  plt_morph%MY        , &
    NRT      =>  plt_morph%NRT       , &
    NBR      =>  plt_morph%NBR       , &
    ARLF1    =>  plt_morph%ARLF1     , &
    ARLFB    =>  plt_morph%ARLFB     , &
    GRNXB    =>  plt_morph%GRNXB     , &
    GRNOB    =>  plt_morph%GRNOB     , &
    ARLFL    =>  plt_morph%ARLFL       &
  )
!     ZNOON=hour of solar noon
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IDAY0,IYR0=day,year of planting
!     IYRC=current year
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!     PPX,PP=PFT population per m2,grid cell
!     FRADP=fraction of radiation received by each PFT canopy
!     VHCPC=canopy heat capacity
!
  IF(J.EQ.INT(ZNOON).AND.(IBTYP(NZ).EQ.0 &
    .OR.IGTYP(NZ).LE.1).AND.(I.NE.IDAY0(NZ) &
    .OR.IYRC.NE.IYR0(NZ)))THEN
    IF(ITILL.LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.IDAY0(NZ).OR.IYRC.GT.IYR0(NZ))THEN
        XHVST=XCORP
        PPX(NZ)=PPX(NZ)*XHVST
        PP(NZ)=PP(NZ)*XHVST
        FRADP(NZ)=FRADP(NZ)*XHVST
        VHCPC(NZ)=VHCPC(NZ)*XHVST
        WTLS(NZ)=0._r8
        WVSTK(NZ)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
        D8975: DO NB=1,NBR(NZ)
          IF(IDTHB(NB,NZ).EQ.0)THEN
            IF(PP(NZ).LE.0.0)IDTHB(NB,NZ)=1
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
              ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+XHVST1 &
                *(CFOPE(instruct,M,ielmc,NZ)*(EPOOL(NB,ielmc,NZ)+EPOLNB(NB,ielmc,NZ) &
                +CPOOLK(NB,NZ)+WTRSVBE(NB,ielmc,NZ)) &
                +CFOPE(ifoliar,M,ielmc,NZ)*(WTLFBE(NB,ielmc,NZ)*FWODB(1) &
                +WTNDBE(NB,ielmc,NZ)) &
                +CFOPE(infoliar,M,ielmc,NZ)*(WTSHEBE(NB,ielmc,NZ)*FWODB(1) &
                +WTHSKBE(NB,ielmc,NZ)+WTEARBE(NB,ielmc,NZ)))
              ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+XHVST1 &
                *CFOPE(icwood,M,ielmc,NZ)*(WTLFBE(NB,ielmc,NZ)*FWODB(0) &
                +WTSHEBE(NB,ielmc,NZ)*FWODB(0))

              ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+XHVST1 &
                *(CFOPE(instruct,M,ielmn,NZ)*(EPOOL(NB,ielmn,NZ)+EPOLNB(NB,ielmn,NZ) &
                +WTRSVBE(NB,ielmn,NZ)) &
                +CFOPE(ifoliar,M,ielmn,NZ)*(WTLFBE(NB,ielmn,NZ)*FWODLN(1) &
                +WTNDBE(NB,ielmn,NZ)) &
                +CFOPE(infoliar,M,ielmn,NZ)*(WTSHEBE(NB,ielmn,NZ)*FWODSN(1) &
                +WTHSKBE(NB,ielmn,NZ)+WTEARBE(NB,ielmn,NZ)))
              ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+XHVST1 &
                *CFOPE(icwood,M,ielmn,NZ)*(WTLFBE(NB,ielmn,NZ)*FWODLN(0) &
                +WTSHEBE(NB,ielmn,NZ)*FWODSN(0))

              ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+XHVST1 &
                *(CFOPE(instruct,M,ielmp,NZ)*(EPOOL(NB,ielmp,NZ)+EPOLNB(NB,ielmp,NZ) &
                +WTRSVBE(NB,ielmp,NZ)) &
                +CFOPE(ifoliar,M,ielmp,NZ)*(WTLFBE(NB,ielmp,NZ)*FWODLP(1) &
                +WTNDBE(NB,ielmp,NZ)) &
                +CFOPE(infoliar,M,ielmp,NZ)*(WTSHEBE(NB,ielmp,NZ)*FWODSP(1) &
                +WTHSKBE(NB,ielmp,NZ)+WTEARBE(NB,ielmp,NZ)))
              ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+XHVST1 &
                *CFOPE(icwood,M,ielmp,NZ)*(WTLFBE(NB,ielmp,NZ)*FWODLP(0) &
                +WTSHEBE(NB,ielmp,NZ)*FWODSP(0))
              IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
                WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+XHVST1 &
                  *CFOPE(infoliar,M,ielmc,NZ)*WTGRBE(NB,ielmc,NZ)
                WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+XHVST1 &
                  *CFOPE(infoliar,M,ielmn,NZ)*WTGRBE(NB,ielmn,NZ)
                WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+XHVST1 &
                  *CFOPE(infoliar,M,ielmp,NZ)*WTGRBE(NB,ielmp,NZ)
              ELSE
                ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+XHVST1 &
                  *CFOPE(infoliar,M,ielmc,NZ)*WTGRBE(NB,ielmc,NZ)
                ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+XHVST1 &
                  *CFOPE(infoliar,M,ielmn,NZ)*WTGRBE(NB,ielmn,NZ)
                ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+XHVST1 &
                  *CFOPE(infoliar,M,ielmp,NZ)*WTGRBE(NB,ielmp,NZ)
              ENDIF
              ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ)+XHVST1 &
                *CFOPE(icwood,M,ielmc,NZ)*WTSTKBE(NB,ielmc,NZ)*FWOODE(ielmc,0)
              ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ)+XHVST1 &
                *CFOPE(icwood,M,ielmn,NZ)*WTSTKBE(NB,ielmn,NZ)*FWOODE(ielmn,0)
              ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ)+XHVST1 &
                *CFOPE(icwood,M,ielmp,NZ)*WTSTKBE(NB,ielmp,NZ)*FWOODE(ielmp,0)
              ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+XHVST1 &
                *CFOPE(istalk,M,ielmc,NZ)*WTSTKBE(NB,ielmc,NZ)*FWOODE(ielmc,1)
              ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+XHVST1 &
                *CFOPE(istalk,M,ielmn,NZ)*WTSTKBE(NB,ielmn,NZ)*FWOODE(ielmn,1)
              ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+XHVST1 &
                *CFOPE(istalk,M,ielmp,NZ)*WTSTKBE(NB,ielmp,NZ)*FWOODE(ielmp,1)
            ENDDO D6380
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!

            CPOOLK(NB,NZ)=CPOOLK(NB,NZ)*XHVST
            EPOOL(NB,ielmc,NZ)=EPOOL(NB,ielmc,NZ)*XHVST
            EPOOL(NB,ielmn,NZ)=EPOOL(NB,ielmn,NZ)*XHVST
            EPOOL(NB,ielmp,NZ)=EPOOL(NB,ielmp,NZ)*XHVST
            EPOLNB(NB,ielmc,NZ)=EPOLNB(NB,ielmc,NZ)*XHVST
            EPOLNB(NB,ielmn,NZ)=EPOLNB(NB,ielmn,NZ)*XHVST
            EPOLNB(NB,ielmp,NZ)=EPOLNB(NB,ielmp,NZ)*XHVST
            WTSHTBE(NB,ielmc,NZ)=WTSHTBE(NB,ielmc,NZ)*XHVST
            WTSHTBE(NB,ielmn,NZ)=WTSHTBE(NB,ielmn,NZ)*XHVST
            WTSHTBE(NB,ielmp,NZ)=WTSHTBE(NB,ielmp,NZ)*XHVST

            WVSTKB(NB,NZ)=WVSTKB(NB,NZ)*XHVST
            DO NE=1,npelms
              WTRSVBE(NB,NE,NZ)=WTRSVBE(NB,NE,NZ)*XHVST
              WTHSKBE(NB,NE,NZ)=WTHSKBE(NB,NE,NZ)*XHVST
              WTEARBE(NB,NE,NZ)=WTEARBE(NB,NE,NZ)*XHVST
              WTGRBE(NB,NE,NZ)=WTGRBE(NB,NE,NZ)*XHVST
              WTLFBE(NB,NE,NZ)=WTLFBE(NB,NE,NZ)*XHVST
              WTNDBE(NB,NE,NZ)=WTNDBE(NB,NE,NZ)*XHVST
              WTSHEBE(NB,NE,NZ)=WTSHEBE(NB,NE,NZ)*XHVST
              WTSTKBE(NB,NE,NZ)=WTSTKBE(NB,NE,NZ)*XHVST
              WTSTXBE(NB,NE,NZ)=WTSTXBE(NB,NE,NZ)*XHVST
            ENDDO

            GRNXB(NB,NZ)=GRNXB(NB,NZ)*XHVST
            GRNOB(NB,NZ)=GRNOB(NB,NZ)*XHVST
            GRWTB(NB,NZ)=GRWTB(NB,NZ)*XHVST
            ARLFB(NB,NZ)=ARLFB(NB,NZ)*XHVST
            WTLSB(NB,NZ)=AZMAX1(WTLFBE(NB,ielmc,NZ)+WTSHEBE(NB,ielmc,NZ))
            WTLS(NZ)=WTLS(NZ)+WTLSB(NB,NZ)

            WVSTK(NZ)=WVSTK(NZ)+WVSTKB(NB,NZ)
            D8970: DO K=0,JNODS1
              IF(K.NE.0)THEN
                CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)*XHVST
                CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)*XHVST
                CO2B(K,NB,NZ)=CO2B(K,NB,NZ)*XHVST
                HCOB(K,NB,NZ)=HCOB(K,NB,NZ)*XHVST
              ENDIF
              ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)*XHVST

              WSLF(K,NB,NZ)=WSLF(K,NB,NZ)*XHVST
!     HTSHE(K,NB,NZ)=HTSHE(K,NB,NZ)*XHVST

              WSSHE(K,NB,NZ)=WSSHE(K,NB,NZ)*XHVST
!     HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)*XHVST
!     HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ)*XHVST
              DO NE=1,npelms
                WGNODE(K,NB,NE,NZ)=WGNODE(K,NB,NE,NZ)*XHVST
                WGLFE(K,NB,NE,NZ)=WGLFE(K,NB,NE,NZ)*XHVST
                WGSHE(K,NB,NE,NZ)=WGSHE(K,NB,NE,NZ)*XHVST
                DO L=1,JC1
                  WGLFLE(L,K,NB,NE,NZ)=WGLFLE(L,K,NB,NE,NZ)*XHVST
                ENDDO
              ENDDO
              D8965: DO L=1,JC1
                ARLFL(L,K,NB,NZ)=ARLFL(L,K,NB,NZ)*XHVST
              ENDDO D8965
            ENDDO D8970
          ENDIF
        ENDDO D8975
!
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        VOLWPX=VOLWP(NZ)
        WVPLT=AZMAX1(WTLS(NZ)+WVSTK(NZ))
        APSILT=ABS(PSILT(NZ))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWP(NZ)=ppmc*WVPLT/FDM
        VOLWOU=VOLWOU+VOLWPX-VOLWP(NZ)
        UVOLO=UVOLO+VOLWPX-VOLWP(NZ)
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
        IF(PP(NZ).LE.0.0)THEN
          IDTHR(NZ)=1
          IDTHP(NZ)=1
          IDTH(NZ)=1
          JHVST(NZ)=1
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
        D8985: DO N=1,MY(NZ)
          D8980: DO L=NU,NJ
            D6385: DO M=1,jsken
              ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+XHVST1 &
                *CFOPE(instruct,M,ielmc,NZ)*EPOOLR(ielmc,N,L,NZ)
              ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+XHVST1 &
                *CFOPE(instruct,M,ielmn,NZ)*EPOOLR(ielmn,N,L,NZ)
              ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+XHVST1 &
                *CFOPE(instruct,M,ielmp,NZ)*EPOOLR(ielmp,N,L,NZ)

              DO NR=1,NRT(NZ)
                ESNC(M,ielmc,0,L,NZ)=ESNC(M,ielmc,0,L,NZ)+XHVST1 &
                  *CFOPE(icwood,M,ielmc,NZ)*(WTRT1E(ielmc,N,L,NR,NZ) &
                  +WTRT2E(ielmc,N,L,NR,NZ))*FWODR(0)
                ESNC(M,ielmn,0,L,NZ)=ESNC(M,ielmn,0,L,NZ)+XHVST1 &
                  *CFOPE(icwood,M,ielmn,NZ)*(WTRT1E(ielmn,N,L,NR,NZ) &
                  +WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(0)
                ESNC(M,ielmp,0,L,NZ)=ESNC(M,ielmp,0,L,NZ)+XHVST1 &
                  *CFOPE(icwood,M,ielmp,NZ)*(WTRT1E(ielmp,N,L,NR,NZ) &
                  +WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(0)
                ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+XHVST1 &
                  *CFOPE(iroot,M,ielmc,NZ)*(WTRT1E(ielmc,N,L,NR,NZ) &
                  +WTRT2E(ielmc,N,L,NR,NZ))*FWODR(1)
                ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+XHVST1 &
                  *CFOPE(iroot,M,ielmn,NZ)*(WTRT1E(ielmn,N,L,NR,NZ) &
                  +WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(1)
                ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+XHVST1 &
                  *CFOPE(iroot,M,ielmp,NZ)*(WTRT1E(ielmp,N,L,NR,NZ) &
                  +WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(1)
              ENDDO
            ENDDO D6385
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Z(NZ)=RCO2Z(NZ)-XHVST1 &
              *(CO2A(N,L,NZ)+CO2P(N,L,NZ))
            ROXYZ(NZ)=ROXYZ(NZ)-XHVST1 &
              *(OXYA(N,L,NZ)+OXYP(N,L,NZ))
            RCH4Z(NZ)=RCH4Z(NZ)-XHVST1 &
              *(CH4A(N,L,NZ)+CH4P(N,L,NZ))
            RN2OZ(NZ)=RN2OZ(NZ)-XHVST1 &
              *(Z2OA(N,L,NZ)+Z2OP(N,L,NZ))
            RNH3Z(NZ)=RNH3Z(NZ)-XHVST1 &
              *(ZH3A(N,L,NZ)+ZH3P(N,L,NZ))
            RH2GZ(NZ)=RH2GZ(NZ)-XHVST1 &
              *(H2GA(N,L,NZ)+H2GP(N,L,NZ))
            CO2A(N,L,NZ)=XHVST*CO2A(N,L,NZ)
            OXYA(N,L,NZ)=XHVST*OXYA(N,L,NZ)
            CH4A(N,L,NZ)=XHVST*CH4A(N,L,NZ)
            Z2OA(N,L,NZ)=XHVST*Z2OA(N,L,NZ)
            ZH3A(N,L,NZ)=XHVST*ZH3A(N,L,NZ)
            H2GA(N,L,NZ)=XHVST*H2GA(N,L,NZ)
            CO2P(N,L,NZ)=XHVST*CO2P(N,L,NZ)
            OXYP(N,L,NZ)=XHVST*OXYP(N,L,NZ)
            CH4P(N,L,NZ)=XHVST*CH4P(N,L,NZ)
            Z2OP(N,L,NZ)=XHVST*Z2OP(N,L,NZ)
            ZH3P(N,L,NZ)=XHVST*ZH3P(N,L,NZ)
            H2GP(N,L,NZ)=XHVST*H2GP(N,L,NZ)
!
!     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D8960: DO NR=1,NRT(NZ)
              DO NE=1,npelms
                WTRT1E(NE,N,L,NR,NZ)=WTRT1E(NE,N,L,NR,NZ)*XHVST
                WTRT2E(NE,N,L,NR,NZ)=WTRT2E(NE,N,L,NR,NZ)*XHVST
                RTWT1E(N,NR,NE,NZ)=RTWT1E(N,NR,NE,NZ)*XHVST
              ENDDO
              RTLG1(N,L,NR,NZ)=RTLG1(N,L,NR,NZ)*XHVST
              RTLG2(N,L,NR,NZ)=RTLG2(N,L,NR,NZ)*XHVST
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST
            ENDDO D8960
            DO NE=1,npelms
              EPOOLR(NE,N,L,NZ)=EPOOLR(NE,N,L,NZ)*XHVST
            ENDDO
            WTRTL(N,L,NZ)=WTRTL(N,L,NZ)*XHVST
            WTRTD(N,L,NZ)=WTRTD(N,L,NZ)*XHVST
            WSRTL(N,L,NZ)=WSRTL(N,L,NZ)*XHVST
            RTN1(N,L,NZ)=RTN1(N,L,NZ)*XHVST
            RTNL(N,L,NZ)=RTNL(N,L,NZ)*XHVST
            RTLGP(N,L,NZ)=RTLGP(N,L,NZ)*XHVST
            RTDNP(N,L,NZ)=RTDNP(N,L,NZ)*XHVST
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
              DO NE=1,npelms
                D6395: DO M=1,jsken
                  ESNC(M,NE,1,L,NZ)=ESNC(M,NE,1,L,NZ)+XHVST1 &
                    *(CFOPE(iroot,M,NE,NZ)*WTNDLE(L,NE,NZ) &
                    +CFOPE(instruct,M,NE,NZ)*EPOOLN(L,NE,NZ))
                ENDDO D6395
                WTNDLE(L,NE,NZ)=WTNDLE(L,NE,NZ)*XHVST
                EPOOLN(L,NE,NZ)=EPOOLN(L,NE,NZ)*XHVST
              ENDDO
            ENDIF
          ENDDO D8980
        ENDDO D8985
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
        D6400: DO M=1,jsken
          ESNC(M,ielmc,0,NG(NZ),NZ)=ESNC(M,ielmc,0,NG(NZ),NZ) &
            +(XHVST1*CFOPE(instruct,M,ielmc,NZ)*WTRVE(ielmc,NZ))*FWOODE(ielmc,0)
          ESNC(M,ielmn,0,NG(NZ),NZ)=ESNC(M,ielmn,0,NG(NZ),NZ) &
            +(XHVST1*CFOPE(instruct,M,ielmn,NZ)*WTRVE(ielmn,NZ))*FWOODE(ielmn,0)
          ESNC(M,ielmp,0,NG(NZ),NZ)=ESNC(M,ielmp,0,NG(NZ),NZ) &
            +(XHVST1*CFOPE(instruct,M,ielmp,NZ)*WTRVE(ielmp,NZ))*FWOODE(ielmp,0)
          ESNC(M,ielmc,1,NG(NZ),NZ)=ESNC(M,ielmc,1,NG(NZ),NZ) &
            +(XHVST1*CFOPE(instruct,M,ielmc,NZ)*WTRVE(ielmc,NZ))*FWOODE(ielmc,1)
          ESNC(M,ielmn,1,NG(NZ),NZ)=ESNC(M,ielmn,1,NG(NZ),NZ) &
            +(XHVST1*CFOPE(instruct,M,ielmn,NZ)*WTRVE(ielmn,NZ))*FWOODE(ielmn,1)
          ESNC(M,ielmp,1,NG(NZ),NZ)=ESNC(M,ielmp,1,NG(NZ),NZ) &
            +(XHVST1*CFOPE(instruct,M,ielmp,NZ)*WTRVE(ielmp,NZ))*FWOODE(ielmp,1)
        ENDDO D6400
        WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)*XHVST
        WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)*XHVST
        WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)*XHVST
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByHarvest(I,J,NZ,CPOOLK)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(inout) :: CPOOLK(JC1,JP1)
  integer :: L,K,M,NR,N,NB,NBX,NE
  real(r8):: ZPOOLG,ZPOLNG,ZPOOLX
  real(r8) :: ZPOLNX,XHVST(npelms)
  real(r8) :: XHVST1(npelms)
  REAL(R8) :: WGLFBL(JC1,JP1,JP1)
  real(r8) :: FHVSHK(0:JNODS1),FHVSTK(0:JNODS1)
  real(r8) :: ARLFY,ARLFR,ARLFG
  real(r8) :: APSILT
  real(r8) :: CPOOLX
  real(r8) :: CCPOLX
  real(r8) :: CPOOLG
  real(r8) :: CPOLNG
  real(r8) :: CPOLNX
  real(r8) :: CCPLNX
  real(r8) :: FHGT
  real(r8) :: FHVST
  real(r8) :: FHVSH
  real(r8) :: FHVST4
  real(r8) :: FHGTK
  real(r8) :: FHVSTS
  real(r8) :: FHVSTG
  real(r8) :: FHVSHG
  real(r8) :: FHVSTH
  real(r8) :: FHVSTE
  real(r8) :: FHVSHH
  real(r8) :: FHVSHE
  real(r8) :: FDM
  real(r8) :: FFIRE,FFIRN,FFIRP
  real(r8) :: FHVSN,FHVSP
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
  real(r8) :: WTNDPG,WGLFGX,WGSHGX,WGLFGY,WGSHGY,WGLFG,WGLFNG
  real(r8) :: WGLFPG,WHVSBS,WHVSCX,WHVSNX,WVPLT
!     begin_execution
  associate(                             &
    HVST     =>  plt_distb%HVST    , &
    EHVST    =>  plt_distb%EHVST   , &
    DCORP    =>  plt_distb%DCORP   , &
    THIN     =>  plt_distb%THIN    , &
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
    VOLWP    =>  plt_ew%VOLWP      , &
    PSILT    =>  plt_ew%PSILT      , &
    HCOB     =>  plt_photo%HCOB    , &
    CO2B     =>  plt_photo%CO2B    , &
    CPOOL3   =>  plt_photo%CPOOL3  , &
    CPOOL4   =>  plt_photo%CPOOL4  , &
    PP       =>  plt_site%PP       , &
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
    WTRTD    => plt_biom%WTRTD     , &
    WTRTL    => plt_biom%WTRTL     , &
    WTRVE    => plt_biom%WTRVE     , &
    WVSTK    => plt_biom%WVSTK     , &
    WTLS     => plt_biom%WTLS      , &
    WVSTKB   => plt_biom%WVSTKB    , &
    WTNDBE   => plt_biom%WTNDBE    , &
    WTRSVE   => plt_biom%WTRSVE    , &
    WTGRBE   => plt_biom%WTGRBE    , &
    WTSTKBE  => plt_biom%WTSTKBE   , &
    WTSHTBE  => plt_biom%WTSHTBE   , &
    WTHSKBE  => plt_biom%WTHSKBE   , &
    WTEARBE  => plt_biom%WTEARBE   , &
    WGNODE   => plt_biom%WGNODE    , &
    WTLSB    => plt_biom%WTLSB     , &
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
    CEPOLP   => plt_biom%CEPOLP    , &
    CCPLNP   => plt_biom%CCPLNP    , &
    WTLFE    => plt_biom%WTLFE     , &
    WTGRE    => plt_biom%WTGRE     , &
    WTSHTE   => plt_biom%WTSHTE     , &
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
    WGLFV    => plt_biom%WGLFV     , &
    FVRN     => plt_allom%FVRN     , &
    FWODR    => plt_allom%FWODR    , &
    FWODRP   => plt_allom%FWODRP   , &
    FWODRN   => plt_allom%FWODRN   , &
    FWODSN   => plt_allom%FWODSN   , &
    FWODSP   => plt_allom%FWODSP   , &
    FWOODE   => plt_allom%FWOODE   , &
    FWODB    => plt_allom%FWODB    , &
    FWODLN   => plt_allom%FWODLN   , &
    FWODLP   => plt_allom%FWODLP   , &
    GRWTB    => plt_allom%GRWTB    , &
    IDTHB    =>  plt_pheno%IDTHB   , &
    TFN3     =>  plt_pheno%TFN3    , &
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
    icwood   =>  pltpar%icwood    , &
    H2GP     =>  plt_rbgc%H2GP     , &
    CO2P     =>  plt_rbgc%CO2P     , &
    CH4P     =>  plt_rbgc%CH4P     , &
    OXYP     =>  plt_rbgc%OXYP     , &
    H2GA     =>  plt_rbgc%H2GA     , &
    CH4A     =>  plt_rbgc%CH4A     , &
    OXYA     =>  plt_rbgc%OXYA     , &
    CO2A     =>  plt_rbgc%CO2A     , &
    ESNC     =>  plt_bgcr%ESNC     , &
    TNBP     =>  plt_bgcr%TNBP     , &
    RCH4Z    =>  plt_bgcr%RCH4Z    , &
    RCO2Z    =>  plt_bgcr%RCO2Z    , &
    ROXYZ    =>  plt_bgcr%ROXYZ    , &
    RN2OZ    =>  plt_bgcr%RN2OZ    , &
    RNH3Z    =>  plt_bgcr%RNH3Z    , &
    RH2GZ    =>  plt_bgcr%RH2GZ    , &
    CNET     =>  plt_bgcr%CNET     , &
    ZH3P     =>  plt_rbgc%ZH3P     , &
    Z2OP     =>  plt_rbgc%Z2OP     , &
    Z2OA     =>  plt_rbgc%Z2OA     , &
    RCO2A    =>  plt_rbgc%RCO2A    , &
    RCO2M    =>  plt_rbgc%RCO2M    , &
    RCO2N    =>  plt_rbgc%RCO2N    , &
    ZH3A     =>  plt_rbgc%ZH3A     , &
    RTNL     =>  plt_morph%RTNL    , &
    RTN2     =>  plt_morph%RTN2    , &
    RTLGP    =>  plt_morph%RTLGP   , &
    RTARP    =>  plt_morph%RTARP   , &
    RTVLP    =>  plt_morph%RTVLP   , &
    NG       =>  plt_morph%NG      , &
    MY       =>  plt_morph%MY      , &
    ZC       =>  plt_morph%ZC      , &
    RTVLW    =>  plt_morph%RTVLW   , &
    RTDNP    =>  plt_morph%RTDNP   , &
    INTYP    =>  plt_morph%INTYP   , &
    ARLFT    =>  plt_morph%ARLFT   , &
    ZL       =>  plt_morph%ZL      , &
    ARLFB    =>  plt_morph%ARLFB   , &
    NBR      =>  plt_morph%NBR     , &
    ARSTP    =>  plt_morph%ARSTP   , &
    HTNODX   =>  plt_morph%HTNODX  , &
    RTN1     =>  plt_morph%RTN1    , &
    RTLG2    =>  plt_morph%RTLG2   , &
    RTLG1    =>  plt_morph%RTLG1   , &
    HTNODE   =>  plt_morph%HTNODE  , &
    GRNXB    =>  plt_morph%GRNXB   , &
    GRNOB    =>  plt_morph%GRNOB   , &
    HTSHE    =>  plt_morph%HTSHE   , &
    ARLF1    =>  plt_morph%ARLF1   , &
    ARLFV    =>  plt_morph%ARLFV   , &
    ARSTV    =>  plt_morph%ARSTV   , &
    ARLFL    =>  plt_morph%ARLFL   , &
    ARSTK    =>  plt_morph%ARSTK   , &
    NRT      =>  plt_morph%NRT     , &
    PSTGF    =>  plt_morph%PSTGF   , &
    NB1      =>  plt_morph%NB1     , &
    PSTGI    =>  plt_morph%PSTGI   , &
    PSTG     =>  plt_morph%PSTG    , &
    CF       =>  plt_morph%CF      , &
    ARLFC    =>  plt_morph%ARLFC   , &
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
!     THIN=thinning:fraction of population removed
!     CF=clumping factor
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     ARLFC,ARLFT=leaf area of combined canopy, canopy layer
!     ARLFR,ARLFY=leaf area harvested,remaining
!     ZL=height to bottom of each canopy layer
!
    IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
      IF(JHVST(NZ).NE.2)THEN
        PPX(NZ)=PPX(NZ)*(1._r8-THIN(NZ))
        PP(NZ)=PP(NZ)*(1._r8-THIN(NZ))
      ELSE
!     PPI(NZ)=AMAX1(1.0_r8,0.5*(PPI(NZ)+GRNO(NZ)/AREA3(NU)))
        PPX(NZ)=PPI(NZ)
        PP(NZ)=PPX(NZ)*AREA3(NU)
      ENDIF
      IF(IHVST(NZ).EQ.3)THEN
        CF(NZ)=CF(NZ)*HVST(NZ)
      ENDIF
      IF(IHVST(NZ).LE.2.AND.HVST(NZ).LT.0.0)THEN
        ARLFY=(1._r8-ABS(HVST(NZ)))*ARLFC
        ARLFR=0._r8
        D9875: DO L=1,JC1
          IF(ZL(L).GT.ZL(L-1).AND.ARLFT(L).GT.ZEROS.AND.ARLFR.LT.ARLFY)THEN
            IF(ARLFR+ARLFT(L).GT.ARLFY)THEN
              HVST(NZ)=ZL(L-1)+((ARLFY-ARLFR)/ARLFT(L))*(ZL(L)-ZL(L-1))
            ENDIF
          ELSE
            HVST(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+ARLFT(L)
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
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WHVSTT=total phytomass grazed, removed
!     TFN3=temperature function for canopy growth
!     CCPOLP=nonstructural C concentration in canopy
!     CCPLNP=nonstructural C concentration in canopy nodules
!
      IF(WTSHTA(NZ).GT.ZEROP(NZ))THEN
        WHVSTT=HVST(NZ)*THIN(NZ)*0.45_r8/24.0_r8 &
          *AREA3(NU)*WTSHTE(ielmc,NZ)/WTSHTA(NZ)
      ELSE
        WHVSTT=0._r8
      ENDIF
      IF(IHVST(NZ).EQ.6)THEN
        WHVSTT=WHVSTT*TFN3(NZ)
      ENDIF
      CCPOLX=CEPOLP(ielmc,NZ)/(1.0_r8+CEPOLP(ielmc,NZ))
      CCPLNX=CCPLNP(NZ)/(1.0_r8+CCPLNP(NZ))
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
      WHVSLX=WHVSTT*EHVST(1,1,NZ)
      WHVSLY=AMIN1(WTLFE(ielmc,NZ),WHVSLX)
      WHVSLF=WHVSLY*(1._r8-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AZMAX1(WHVSLX-WHVSLY)
      WHVSSX=WHVSTT*EHVST(1,2,NZ)
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
      WHVSKX=WHVSTT*EHVST(1,3,NZ)
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
        IF(WHVXXX.GT.0.0)THEN
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
      D9860: DO NB=1,NBR(NZ)
        DO  L=1,JC1
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=0._r8
          enddo
        enddo
      ENDDO D9860
      D9870: DO NB=1,NBR(NZ)
        DO  L=1,JC1
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=WGLFBL(L,NB,NZ)+WGLFLE(L,K,NB,ielmc,NZ)
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
!     FHVST=fraction of canopy layer mass not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
    D9865: DO L=JC1,1,-1
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        IF(IHVST(NZ).NE.3)THEN
          IF(ZL(L).GT.ZL(L-1))THEN
            FHGT=AZMAX1(AMIN1(1.0_r8,1._r8-((ZL(L))-HVST(NZ))/(ZL(L)-ZL(L-1))))
          ELSE
            FHGT=1.0_r8
          ENDIF
        ELSE
          FHGT=0._r8
        ENDIF
        IF(test_aeqb(THIN(NZ),0._r8))THEN
          FHVST=AZMAX1(1._r8-(1._r8-FHGT)*EHVST(1,1,NZ))
          FHVSH=FHVST
        ELSE
          FHVST=AZMAX1(1._r8-THIN(NZ))
          IF(IHVST(NZ).EQ.0)THEN
            FHVSH=1.0_r8-(1._r8-FHGT)*EHVST(1,1,NZ)*THIN(NZ)
          ELSE
            FHVSH=FHVST
          ENDIF
        ENDIF
      ELSE
        FHVST=0._r8
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
!     FHVST=fraction of leaf node mass not harvested
!
      DO 9855 NB=1,NBR(NZ)
        IF((IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6) &
          .AND.WTLFE(ielmc,NZ).GT.ZEROL(NZ))THEN
          WHVSBL=WHVSLF*AZMAX1(WGLFBL(L,NB,NZ))/WTLFE(ielmc,NZ)
        ELSE
          WHVSBL=0._r8
        ENDIF
        DO 9845 K=JNODS1,0,-1
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBL.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(WGLFLE(L,K,NB,ielmc,NZ).GT.WHVSBL)THEN
                FHVST=AZMAX1(AMIN1(1.0_r8,(WGLFLE(L,K,NB,ielmc,NZ)-WHVSBL)/WGLFLE(L,K,NB,ielmc,NZ)))
                FHVSH=FHVST
              ELSE
                FHVST=1.0_r8
                FHVSH=1.0_r8
              ENDIF
            ENDIF
        !
!     HARVESTED LEAF AREA, C, N, P
!
!     FHVST=fraction of leaf node mass not harvested
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     ARLFL,ARSTK=leaf,stalk node area in canopy layer
!     WTHTH1E(ielmc),WTHTH1E(ielmn),WTHTH1E(ielmp)=harvested leaf C,N,P
!     WTHTX1E(ielmc),WTHTX1E(ielmn),WTHTX1E(ielmp)=harvested leaf C,N,P to litter
!     WTHTH3E(ielmc),WTHTH3E(ielmn),WTHTH3E(ielmp)=harvested woody C,N,P
!     WTHTX3E(ielmc),WTHTX3E(ielmn),WTHTX3E(ielmp)=harvested woody C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
            WHVSBL=WHVSBL-(1._r8-FHVST)*WGLFLE(L,K,NB,ielmc,NZ)
            WTHTH1E(ielmc)=WTHTH1E(ielmc)+(1._r8-FHVSH)*WGLFLE(L,K,NB,ielmc,NZ)*FWODB(1)
            WTHTH1E(ielmn)=WTHTH1E(ielmn)+(1._r8-FHVSH)*WGLFLE(L,K,NB,ielmn,NZ)*FWODLN(1)
            WTHTH1E(ielmp)=WTHTH1E(ielmp)+(1._r8-FHVSH)*WGLFLE(L,K,NB,ielmp,NZ)*FWODLP(1)
            WTHTX1E(ielmc)=WTHTX1E(ielmc)+(FHVSH-FHVST)*WGLFLE(L,K,NB,ielmc,NZ)*FWODB(1)
            WTHTX1E(ielmn)=WTHTX1E(ielmn)+(FHVSH-FHVST)*WGLFLE(L,K,NB,ielmn,NZ)*FWODLN(1)
            WTHTX1E(ielmp)=WTHTX1E(ielmp)+(FHVSH-FHVST)*WGLFLE(L,K,NB,ielmp,NZ)*FWODLP(1)
            WTHTH3E(ielmc)=WTHTH3E(ielmc)+(1._r8-FHVSH)*WGLFLE(L,K,NB,ielmc,NZ)*FWODB(0)
            WTHTH3E(ielmn)=WTHTH3E(ielmn)+(1._r8-FHVSH)*WGLFLE(L,K,NB,ielmn,NZ)*FWODLN(0)
            WTHTH3E(ielmp)=WTHTH3E(ielmp)+(1._r8-FHVSH)*WGLFLE(L,K,NB,ielmp,NZ)*FWODLP(0)
            WTHTX3E(ielmc)=WTHTX3E(ielmc)+(FHVSH-FHVST)*WGLFLE(L,K,NB,ielmc,NZ)*FWODB(0)
            WTHTX3E(ielmn)=WTHTX3E(ielmn)+(FHVSH-FHVST)*WGLFLE(L,K,NB,ielmn,NZ)*FWODLN(0)
            WTHTX3E(ielmp)=WTHTX3E(ielmp)+(FHVSH-FHVST)*WGLFLE(L,K,NB,ielmp,NZ)*FWODLP(0)
!
!     REMAINING LEAF C,N,P AND AREA
!
            WGLFLE(L,K,NB,ielmc,NZ)=FHVST*WGLFLE(L,K,NB,ielmc,NZ)
            WGLFLE(L,K,NB,ielmn,NZ)=FHVST*WGLFLE(L,K,NB,ielmn,NZ)
            WGLFLE(L,K,NB,ielmp,NZ)=FHVST*WGLFLE(L,K,NB,ielmp,NZ)
            ARLFL(L,K,NB,NZ)=FHVST*ARLFL(L,K,NB,NZ)
            IF(K.EQ.1)THEN
              ARSTK(L,NB,NZ)=FHVST*ARSTK(L,NB,NZ)
            ENDIF
          ENDIF

9845      CONTINUE
9855    CONTINUE
        ARLFV(L,NZ)=0._r8
        WGLFV(L,NZ)=0._r8
        ARSTV(L,NZ)=ARSTV(L,NZ)*FHVST
      ENDDO D9865
      D9835: DO NB=1,NBR(NZ)
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
          WGLFG=0._r8
          WGLFNG=0._r8
          WGLFPG=0._r8
!
!     ACCUMULATE REMAINING LEAF AREA, C, N, P
!
!     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
!     ARLFL,ARLFV=leaf node,total area in canopy layer
!
          D9815: DO L=1,JC1
            ARLFG=ARLFG+ARLFL(L,K,NB,NZ)
            WGLFG=WGLFG+WGLFLE(L,K,NB,ielmc,NZ)
            WGLFNG=WGLFNG+WGLFLE(L,K,NB,ielmn,NZ)
            WGLFPG=WGLFPG+WGLFLE(L,K,NB,ielmp,NZ)
            ARLFV(L,NZ)=ARLFV(L,NZ)+ARLFL(L,K,NB,NZ)
            WGLFV(L,NZ)=WGLFV(L,NZ)+WGLFLE(L,K,NB,ielmc,NZ)
          ENDDO D9815
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FHVSTK=fraction of internode layer mass not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
          IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
            IF(WGLFE(K,NB,ielmc,NZ).GT.ZEROP(NZ).AND.EHVST(1,1,NZ).GT.0.0)THEN
              FHVSTK(K)=AZMAX1(AMIN1(1.0_r8,(1._r8-(1._r8-AZMAX1(WGLFG) &
                /WGLFE(K,NB,ielmc,NZ))*EHVST(1,2,NZ)/EHVST(1,1,NZ))))
              FHVSHK(K)=FHVSTK(K)
          ELSE
            IF(test_aeqb(THIN(NZ),0._r8))THEN
              FHVSTK(K)=1.0_r8-EHVST(1,2,NZ)
              FHVSHK(K)=FHVSTK(K)
            ELSE
              FHVSTK(K)=1.0_r8-THIN(NZ)
              IF(IHVST(NZ).EQ.0)THEN
                FHVSHK(K)=1.0_r8-EHVST(1,2,NZ)*THIN(NZ)
              ELSE
                FHVSHK(K)=FHVSTK(K)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          FHVSTK(K)=0._r8
          FHVSHK(K)=0._r8
        ENDIF
!
!     ACCUMULATE REMAINING BRANCH LEAF AREA, C, N, P
!
!     WGLF=leaf node C mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     ARLFB,ARLF=branch,node leaf area
!     WSLF=leaf protein mass
!
        WGLFGY=WGLFGY+WGLFE(K,NB,ielmc,NZ)
        WTLFBE(NB,ielmc,NZ)=WTLFBE(NB,ielmc,NZ)-WGLFE(K,NB,ielmc,NZ)+WGLFG
        WTLFBE(NB,ielmn,NZ)=WTLFBE(NB,ielmn,NZ)-WGLFE(K,NB,ielmn,NZ)+WGLFNG
        WTLFBE(NB,ielmp,NZ)=WTLFBE(NB,ielmp,NZ)-WGLFE(K,NB,ielmp,NZ)+WGLFPG
        ARLFB(NB,NZ)=ARLFB(NB,NZ)-ARLF1(K,NB,NZ)+ARLFG
        IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ))THEN
          WSLF(K,NB,NZ)=WSLF(K,NB,NZ)*ARLFG/ARLF1(K,NB,NZ)
        ELSE
          WSLF(K,NB,NZ)=0._r8
        ENDIF
        ARLF1(K,NB,NZ)=ARLFG
        WGLFE(K,NB,ielmc,NZ)=WGLFG
        WGLFE(K,NB,ielmn,NZ)=WGLFNG
        WGLFE(K,NB,ielmp,NZ)=WGLFPG
        WGLFGX=WGLFGX+WGLFE(K,NB,ielmc,NZ)
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
        WHVSBS=WHVSHH*WTSHEBE(NB,ielmc,NZ)/WTSHEE(ielmc,NZ)
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
!     FHVSTK=fraction of internode layer mass not harvested
!     WTHTH2E(ielmc),WTHTH2E(ielmn),WTHTH2E(ielmp)=harvested petiole C,N,P
!     WTHTX2E(ielmc),WTHTX2E(ielmn),WTHTX2E(ielmp)=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     HTSHE,HTNODE=petiole,internode length
!
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBS.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(WGSHE(K,NB,ielmc,NZ).GT.WHVSBS)THEN
                FHVSTK(K)=AZMAX1(AMIN1(1.0_r8,(WGSHE(K,NB,ielmc,NZ)-WHVSBS)/WGSHE(K,NB,ielmc,NZ)))
                FHVSHK(K)=FHVSTK(K)
              ELSE
                FHVSTK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1._r8-FHVSTK(K))*WGSHE(K,NB,ielmc,NZ)
            WTHTH2E(ielmc)=WTHTH2E(ielmc)+(1._r8-FHVSHK(K))*WGSHE(K,NB,ielmc,NZ)*FWODB(1)
            WTHTH2E(ielmn)=WTHTH2E(ielmn)+(1._r8-FHVSHK(K))*WGSHE(K,NB,ielmn,NZ)*FWODSN(1)
            WTHTH2E(ielmp)=WTHTH2E(ielmp)+(1._r8-FHVSHK(K))*WGSHE(K,NB,ielmp,NZ)*FWODSP(1)
            WTHTX2E(ielmc)=WTHTX2E(ielmc)+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,ielmc,NZ)*FWODB(1)
            WTHTX2E(ielmn)=WTHTX2E(ielmn)+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,ielmn,NZ)*FWODSN(1)
            WTHTX2E(ielmp)=WTHTX2E(ielmp)+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,ielmp,NZ)*FWODSP(1)
            WTHTH3E(ielmc)=WTHTH3E(ielmc)+(1._r8-FHVSHK(K))*WGSHE(K,NB,ielmc,NZ)*FWODB(0)
            WTHTH3E(ielmn)=WTHTH3E(ielmn)+(1._r8-FHVSHK(K))*WGSHE(K,NB,ielmn,NZ)*FWODSN(0)
            WTHTH3E(ielmp)=WTHTH3E(ielmp)+(1._r8-FHVSHK(K))*WGSHE(K,NB,ielmp,NZ)*FWODSP(0)
            WTHTX3E(ielmc)=WTHTX3E(ielmc)+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,ielmc,NZ)*FWODB(0)
            WTHTX3E(ielmn)=WTHTX3E(ielmn)+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,ielmn,NZ)*FWODSN(0)
            WTHTX3E(ielmp)=WTHTX3E(ielmp)+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,ielmp,NZ)*FWODSP(0)
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     WGSHE=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     HTSHE=node petiole height
!     WSSHE=petiole protein mass
!
            WGSHGY=WGSHGY+WGSHE(K,NB,ielmc,NZ)
            WTSHEBE(NB,ielmc,NZ)=WTSHEBE(NB,ielmc,NZ) &
              -(1._r8-FHVSTK(K))*WGSHE(K,NB,ielmc,NZ)
            WTSHEBE(NB,ielmn,NZ)=WTSHEBE(NB,ielmn,NZ) &
              -(1._r8-FHVSTK(K))*WGSHE(K,NB,ielmn,NZ)
            WTSHEBE(NB,ielmp,NZ)=WTSHEBE(NB,ielmp,NZ) &
              -(1._r8-FHVSTK(K))*WGSHE(K,NB,ielmp,NZ)
            WGSHE(K,NB,ielmc,NZ)=FHVSTK(K)*WGSHE(K,NB,ielmc,NZ)
            WSSHE(K,NB,NZ)=FHVSTK(K)*WSSHE(K,NB,NZ)
            WGSHE(K,NB,ielmn,NZ)=FHVSTK(K)*WGSHE(K,NB,ielmn,NZ)
            WGSHE(K,NB,ielmp,NZ)=FHVSTK(K)*WGSHE(K,NB,ielmp,NZ)
            WSSHE(K,NB,NZ)=FHVSTK(K)*WSSHE(K,NB,NZ)
            IF(IHVST(NZ).LE.2.AND.HTSHE(K,NB,NZ).GT.0.0)THEN
              FHGT=AZMAX1(AMIN1(1.0_r8,(HTNODE(K,NB,NZ) &
                +HTSHE(K,NB,NZ)-HVST(NZ))/HTSHE(K,NB,NZ)))
              HTSHE(K,NB,NZ)=(1._r8-FHGT)*HTSHE(K,NB,NZ)
            ELSE
              HTSHE(K,NB,NZ)=FHVSTK(K)*HTSHE(K,NB,NZ)
            ENDIF
            WGSHGX=WGSHGX+WGSHE(K,NB,ielmc,NZ)
!     IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
!     IF(HTNODE(K,NB,NZ).GT.HVST(NZ)
!    2.OR.IHVST(NZ).EQ.3)THEN
!     IF(test_aeqb(FHVSTK(K),0._r8).AND.K.GT.0)THEN
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
!     FHVST=fraction of leaf+petiole node mass not harvested
!     CPOOLG,ZPOOLG,PPOOLG=branch non-structural C,N,P mass after harvest
!     CPOLNG,ZPOLNG,PPOLNG=nonstructural C,N,P in bacteria after harvest
!     WTNDG,WTNDNG,WTNDPG=bacterial C,N,P mass after harvest
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!     WHVSC*=nonstructural C removed
!
        CPOOLX=AZMAX1(EPOOL(NB,ielmc,NZ))
        ZPOOLX=AZMAX1(EPOOL(NB,ielmn,NZ))
        PPOOLX=AZMAX1(EPOOL(NB,ielmp,NZ))
        CPOLNX=AZMAX1(EPOLNB(NB,ielmc,NZ))
        ZPOLNX=AZMAX1(EPOLNB(NB,ielmn,NZ))
        PPOLNX=AZMAX1(EPOLNB(NB,ielmp,NZ))
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(WGLFGY+WGSHGY.GT.ZEROP(NZ))THEN
            FHVST=AZMAX1(AMIN1(1.0_r8,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
            CPOOLG=CPOOLX*FHVST
            ZPOOLG=ZPOOLX*FHVST
            PPOOLG=PPOOLX*FHVST
            CPOLNG=CPOLNX*FHVST
            ZPOLNG=ZPOLNX*FHVST
            PPOLNG=PPOLNX*FHVST
            WTNDG=WTNDBE(NB,ielmc,NZ)*FHVST
            WTNDNG=WTNDBE(NB,ielmn,NZ)*FHVST
            WTNDPG=WTNDBE(NB,ielmp,NZ)*FHVST
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
          IF(WTLS(NZ).GT.ZEROL(NZ))THEN
            WTLSBX=AZMAX1(WTLSB(NB,NZ))
            IF(EPOOL(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
              WHVSCX=AZMAX1(WHVSCP)*WTLSBX/WTLS(NZ)
              CPOOLG=AZMAX1(CPOOLX-WHVSCX)
              ZPOOLG=AZMAX1(ZPOOLX-WHVSCX*ZPOOLX/EPOOL(NB,ielmc,NZ))
              PPOOLG=AZMAX1(PPOOLX-WHVSCX*PPOOLX/EPOOL(NB,ielmc,NZ))
            ELSE
              CPOOLG=0._r8
              ZPOOLG=0._r8
              PPOOLG=0._r8
            ENDIF
            IF(EPOLNB(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
              WHVSNX=AZMAX1(WHVSNP)*WTLSBX/WTLS(NZ)
              CPOLNG=AZMAX1(CPOLNX-WHVSNX)
              ZPOLNG=AZMAX1(ZPOLNX-WHVSNX*ZPOLNX/EPOLNB(NB,ielmc,NZ))
              PPOLNG=AZMAX1(PPOLNX-WHVSNX*PPOLNX/EPOLNB(NB,ielmc,NZ))
              WTNDG=WTNDBE(NB,ielmc,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDNG=WTNDBE(NB,ielmn,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDPG=WTNDBE(NB,ielmp,NZ)*(1._r8-WHVSNX/CPOLNX)
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
        WTHTH0E(ielmc)=WTHTH0E(ielmc)+WTNDBE(NB,ielmc,NZ)-WTNDG
        WTHTH0E(ielmn)=WTHTH0E(ielmn)+WTNDBE(NB,ielmn,NZ)-WTNDNG
        WTHTH0E(ielmp)=WTHTH0E(ielmp)+WTNDBE(NB,ielmp,NZ)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
        EPOOL(NB,ielmc,NZ)=CPOOLG
        EPOOL(NB,ielmn,NZ)=ZPOOLG
        EPOOL(NB,ielmp,NZ)=PPOOLG
        EPOLNB(NB,ielmc,NZ)=CPOLNG
        EPOLNB(NB,ielmn,NZ)=ZPOLNG
        EPOLNB(NB,ielmp,NZ)=PPOLNG
        WTNDBE(NB,ielmc,NZ)=WTNDG
        WTNDBE(NB,ielmn,NZ)=WTNDNG
        WTNDBE(NB,ielmp,NZ)=WTNDPG
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
        IF(ICTYP(NZ).EQ.4.AND.CPOOLX.GT.ZEROP(NZ))THEN
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
!     FHVST=fraction of canopy layer mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     THIN=IHVST=0-3,5: fraction of population removed,
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
            IF(test_aeqb(THIN(NZ),0._r8))THEN
              FHVST=AZMAX1(1._r8-(1._r8-FHGT)*EHVST(1,3,NZ))
              FHVSH=FHVST
            ELSE
              FHVST=AZMAX1(1._r8-THIN(NZ))
              IF(IHVST(NZ).EQ.0)THEN
                FHVSH=1.0_r8-(1._r8-FHGT)*EHVST(1,3,NZ)*THIN(NZ)
              ELSE
                FHVSH=FHVST
              ENDIF
            ENDIF
          ELSE
            FHVST=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ELSE
          IF(WTSTKE(ielmc,NZ).GT.ZEROL(NZ))THEN
            FHVST=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/WTSTKE(ielmc,NZ)))
            FHVSH=FHVST
          ELSE
            FHVST=1.0_r8
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
        DO NE=1,npelms
          WTHTH3E(NE)=WTHTH3E(NE)+(1._r8-FHVSH)*WTSTKBE(NB,NE,NZ)
          WTHTX3E(NE)=WTHTX3E(NE)+(FHVSH-FHVST)*WTSTKBE(NB,NE,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
          WTSTKBE(NB,NE,NZ)=FHVST*WTSTKBE(NB,NE,NZ)
          WTSTXBE(NB,NE,NZ)=FHVST*WTSTXBE(NB,NE,NZ)
        ENDDO

        WVSTKB(NB,NZ)=FHVST*WVSTKB(NB,NZ)
!
!     CUT STALK NODES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTNODX,HTNODE=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
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
                FHGTK=AZMAX1(AMIN1(1.0_r8,(HTNODE(K,NB,NZ) &
                  -HVST(NZ))/HTNODX(K,NB,NZ)))
              ELSE
                FHGTK=0._r8
              ENDIF
              IF(test_aeqb(THIN(NZ),0._r8))THEN
                FHVSTS=AZMAX1(1._r8-FHGTK*EHVST(1,3,NZ))
              ELSE
                FHVSTS=AZMAX1(1._r8-THIN(NZ))
              ENDIF
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ELSE
            IF(WTSTKE(ielmc,NZ).GT.ZEROP(NZ))THEN
              FHVSTS=AZMAX1(AMIN1(1.0_r8,1._r8-WHVSTH/WTSTKE(ielmc,NZ)))
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ENDIF
          WGNODE(K,NB,ielmc,NZ)=FHVSTS*WGNODE(K,NB,ielmc,NZ)
          WGNODE(K,NB,ielmn,NZ)=FHVSTS*WGNODE(K,NB,ielmn,NZ)
          WGNODE(K,NB,ielmp,NZ)=FHVSTS*WGNODE(K,NB,ielmp,NZ)
          IF(IHVST(NZ).LE.2.AND.test_aeqb(THIN(NZ),0._r8))THEN
            HTNODX(K,NB,NZ)=FHVSTS*HTNODX(K,NB,NZ)
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
!     FHVST=fraction of reserve mass not harvested
!
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(WTSTKBE(NB,ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVST=FHVST
            FHVSH=FHVSH
          ELSE
            FHVST=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(WTRSVE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVST=AZMAX1(AMIN1(1.0_r8,1._r8-WHVRVH/WTRSVE(ielmc,NZ)))
            FHVSH=FHVST
          ELSE
            FHVST=0._r8
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
        WTHTH3E(ielmc)=WTHTH3E(ielmc)+(1._r8-FHVSH)*WTRSVBE(NB,ielmc,NZ)
        WTHTH3E(ielmn)=WTHTH3E(ielmn)+(1._r8-FHVSH)*WTRSVBE(NB,ielmn,NZ)
        WTHTH3E(ielmp)=WTHTH3E(ielmp)+(1._r8-FHVSH)*WTRSVBE(NB,ielmp,NZ)
        WTHTX3E(ielmc)=WTHTX3E(ielmc)+(FHVSH-FHVST)*WTRSVBE(NB,ielmc,NZ)
        WTHTX3E(ielmn)=WTHTX3E(ielmn)+(FHVSH-FHVST)*WTRSVBE(NB,ielmn,NZ)
        WTHTX3E(ielmp)=WTHTX3E(ielmp)+(FHVSH-FHVST)*WTRSVBE(NB,ielmp,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
        WTRSVBE(NB,ielmc,NZ)=FHVST*WTRSVBE(NB,ielmc,NZ)
        WTRSVBE(NB,ielmn,NZ)=FHVST*WTRSVBE(NB,ielmn,NZ)
        WTRSVBE(NB,ielmp,NZ)=FHVST*WTRSVBE(NB,ielmp,NZ)
!
!     CUT REPRODUCTIVE ORGANS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FHVSTG,FHVSTH,FHVSTE=fraction of grain,husk,ear mass not harvested
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(HVST(NZ).LT.HTSTKX.OR.IHVST(NZ).EQ.1.OR.IHVST(NZ).EQ.3)THEN
            IF(test_aeqb(THIN(NZ),0._r8))THEN
              FHVSTG=1.0_r8-EHVST(1,2,NZ)
              FHVSHG=FHVSTG
            ELSE
              FHVSTG=1.0_r8-THIN(NZ)
              FHVSHG=1.0_r8-EHVST(1,2,NZ)*THIN(NZ)
            ENDIF
          ELSE
            FHVSTG=1.0_r8-THIN(NZ)
            FHVSHG=FHVSTG
          ENDIF
          FHVSTH=FHVSTG
          FHVSTE=FHVSTG
          FHVSHH=FHVSHG
          FHVSHE=FHVSHG
        ELSE
          IF(WTHSKE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSTH=AZMAX1(AMIN1(1.0_r8,1._r8-WHVHSH/WTHSKE(ielmc,NZ)))
            FHVSHH=FHVSTH
          ELSE
            FHVSTH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(WTEARE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSTE=AZMAX1(AMIN1(1.0_r8,1._r8-WHVEAH/WTEARE(ielmc,NZ)))
            FHVSHE=FHVSTE
          ELSE
            FHVSTE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(WTGRE(ielmc,NZ).GT.ZEROP(NZ))THEN
            FHVSTG=AZMAX1(AMIN1(1.0_r8,1._r8-WHVGRH/WTGRE(ielmc,NZ)))
            FHVSHG=FHVSTG
          ELSE
            FHVSTG=1.0_r8
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
        WTHTH2E(ielmc)=WTHTH2E(ielmc)+(1._r8-FHVSHH)*WTHSKBE(NB,ielmc,NZ)+(1._r8-FHVSHE) &
          *WTEARBE(NB,ielmc,NZ)+(1._r8-FHVSHG)*WTGRBE(NB,ielmc,NZ)
        WTHTH2E(ielmn)=WTHTH2E(ielmn)+(1._r8-FHVSHH)*WTHSKBE(NB,ielmn,NZ)+(1._r8-FHVSHE) &
          *WTEARBE(NB,ielmn,NZ)+(1._r8-FHVSHG)*WTGRBE(NB,ielmn,NZ)
        WTHTH2E(ielmp)=WTHTH2E(ielmp)+(1._r8-FHVSHH)*WTHSKBE(NB,ielmp,NZ)+(1._r8-FHVSHE) &
          *WTEARBE(NB,ielmp,NZ)+(1._r8-FHVSHG)*WTGRBE(NB,ielmp,NZ)
        WTHTX2E(ielmc)=WTHTX2E(ielmc)+(FHVSHH-FHVSTH)*WTHSKBE(NB,ielmc,NZ)+(FHVSHE-FHVSTE) &
          *WTEARBE(NB,ielmc,NZ)+(FHVSHG-FHVSTG)*WTGRBE(NB,ielmc,NZ)
        WTHTX2E(ielmn)=WTHTX2E(ielmn)+(FHVSHH-FHVSTH)*WTHSKBE(NB,ielmn,NZ)+(FHVSHE-FHVSTE) &
          *WTEARBE(NB,ielmn,NZ)+(FHVSHG-FHVSTG)*WTGRBE(NB,ielmn,NZ)
        WTHTX2E(ielmp)=WTHTX2E(ielmp)+(FHVSHH-FHVSTH)*WTHSKBE(NB,ielmp,NZ)+(FHVSHE-FHVSTE) &
          *WTEARBE(NB,ielmp,NZ)+(FHVSHG-FHVSTG)*WTGRBE(NB,ielmp,NZ)
        WTHTGE(ielmc)=WTHTGE(ielmc)+(1._r8-FHVSTG)*WTGRBE(NB,ielmc,NZ)
        WTHTGE(ielmn)=WTHTGE(ielmn)+(1._r8-FHVSTG)*WTGRBE(NB,ielmn,NZ)
        WTHTGE(ielmp)=WTHTGE(ielmp)+(1._r8-FHVSTG)*WTGRBE(NB,ielmp,NZ)
!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
        WTHSKBE(NB,ielmc,NZ)=FHVSTH*WTHSKBE(NB,ielmc,NZ)
        WTEARBE(NB,ielmc,NZ)=FHVSTE*WTEARBE(NB,ielmc,NZ)
        WTGRBE(NB,ielmc,NZ)=FHVSTG*WTGRBE(NB,ielmc,NZ)
        WTHSKBE(NB,ielmn,NZ)=FHVSTH*WTHSKBE(NB,ielmn,NZ)
        WTEARBE(NB,ielmn,NZ)=FHVSTE*WTEARBE(NB,ielmn,NZ)
        WTGRBE(NB,ielmn,NZ)=FHVSTG*WTGRBE(NB,ielmn,NZ)
        WTHSKBE(NB,ielmp,NZ)=FHVSTH*WTHSKBE(NB,ielmp,NZ)
        WTEARBE(NB,ielmp,NZ)=FHVSTE*WTEARBE(NB,ielmp,NZ)
        WTGRBE(NB,ielmp,NZ)=FHVSTG*WTGRBE(NB,ielmp,NZ)
        GRNXB(NB,NZ)=FHVSTG*GRNXB(NB,NZ)
        GRNOB(NB,NZ)=FHVSTG*GRNOB(NB,NZ)
        GRWTB(NB,NZ)=FHVSTG*GRWTB(NB,NZ)
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
!     WVSTKB=stalk sapwood mass
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        CPOOLK(NB,NZ)=0._r8
        D1325: DO K=1,JNODS1
          CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
            +CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
            +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
        ENDDO D1325
        WTLSB(NB,NZ)=AZMAX1(WTLFBE(NB,ielmc,NZ) &
          +WTSHEBE(NB,ielmc,NZ))
        WTSHTBE(NB,ielmc,NZ)=AZMAX1(WTLFBE(NB,ielmc,NZ) &
          +WTSHEBE(NB,ielmc,NZ)+WTSTKBE(NB,ielmc,NZ)+WTRSVBE(NB,ielmc,NZ) &
          +WTHSKBE(NB,ielmc,NZ)+WTEARBE(NB,ielmc,NZ)+WTGRBE(NB,ielmc,NZ) &
          +EPOOL(NB,ielmc,NZ)+CPOOLK(NB,NZ))
        WTSHTBE(NB,ielmn,NZ)=AZMAX1(WTLFBE(NB,ielmn,NZ) &
          +WTSHEBE(NB,ielmn,NZ)+WTSTKBE(NB,ielmn,NZ)+WTRSVBE(NB,ielmn,NZ) &
          +WTHSKBE(NB,ielmn,NZ)+WTEARBE(NB,ielmn,NZ)+WTGRBE(NB,ielmn,NZ) &
          +EPOOL(NB,ielmn,NZ))
        WTSHTBE(NB,ielmp,NZ)=AZMAX1(WTLFBE(NB,ielmp,NZ) &
          +WTSHEBE(NB,ielmp,NZ)+WTSTKBE(NB,ielmp,NZ)+WTRSVBE(NB,ielmp,NZ) &
          +WTHSKBE(NB,ielmp,NZ)+WTEARBE(NB,ielmp,NZ)+WTGRBE(NB,ielmp,NZ) &
          +EPOOL(NB,ielmp,NZ))
        VOLWPX=VOLWP(NZ)
        WVPLT=AZMAX1(WTLS(NZ)+WVSTK(NZ))
        APSILT=ABS(PSILT(NZ))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWP(NZ)=ppmc*WVPLT/FDM
        VOLWOU=VOLWOU+VOLWPX-VOLWP(NZ)
        UVOLO=UVOLO+VOLWPX-VOLWP(NZ)
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
          .AND.ZC(NZ).GT.HVST(NZ))THEN
          IF((IWTYP(NZ).NE.0.AND.VRNF(NB,NZ) &
            .LE.FVRN(IWTYP(NZ))*VRNX(NB,NZ)) &
            .OR.(IWTYP(NZ).EQ.0 &
            .AND.IDAY(1,NB,NZ).NE.0))THEN
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
              D3010: DO NBX=1,NBR(NZ)
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
!     ARSTK=total PFT stalk surface area
!
      IF(JHVST(NZ).NE.0)IDTHB(NB,NZ)=1
      IF(PP(NZ).LE.0.0)IDTHB(NB,NZ)=1
      ENDDO D9835
      WTLS(NZ)=0._r8
      WTSTKE(ielmc,NZ)=0._r8
      WVSTK(NZ)=0._r8
      ARSTP(NZ)=0._r8
      D9840: DO NB=1,NBR(NZ)
        WTLS(NZ)=WTLS(NZ)+WTLSB(NB,NZ)
        WTSTKE(ielmc,NZ)=WTSTKE(ielmc,NZ)+WTSTKBE(NB,ielmc,NZ)
        WVSTK(NZ)=WVSTK(NZ)+WVSTKB(NB,NZ)
        D9830: DO L=1,JC1
          ARSTP(NZ)=ARSTP(NZ)+ARSTK(L,NB,NZ)
        ENDDO D9830
      ENDDO D9840
!
!     ROOT LITTERFALL FROM HARVESTING OR FIRE
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     THETW=soil water concentration
!     CORGC=SOC concentration
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     EFIRE=combustion  of N,P relative to C
!     FHVST,FHVSN,FHVSP=fraction of root layer C,N,P not removed by disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CNET=PFT net CO2 fixation
!     TNBP=total net biome productivity
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        XHVST(ielmc)=1.0_r8-THIN(NZ)
        D3985: DO N=1,MY(NZ)
          D3980: DO L=NU,NJ
            IF(IHVST(NZ).NE.5)THEN
              XHVST(ielmc)=1.0_r8-THIN(NZ)
              XHVST(ielmn)=XHVST(ielmc)
              XHVST(ielmp)=XHVST(ielmc)
              FFIRE=0._r8
              FFIRN=0._r8
              FFIRP=0._r8
            ELSE
              IF(THETW(L).GT.FVLWB.OR.CORGC(L).LE.FORGC.OR.ITILL.NE.22)THEN
                XHVST(ielmc)=1.0_r8
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE=0._r8
                FFIRN=0._r8
                FFIRP=0._r8
              ELSE
                XHVST(ielmc)=1.0_r8-DCORP*EHVST(1,3,NZ) &
                  *AMIN1(1.0_r8,(CORGC(L)-FORGC)/(0.55E+06-FORGC))
                XHVST(ielmn)=XHVST(ielmc)
                XHVST(ielmp)=XHVST(ielmc)
                FFIRE=EHVST(2,3,NZ)
                FFIRN=FFIRE*EFIRE(1,IHVST(NZ))
                FFIRP=FFIRE*EFIRE(2,IHVST(NZ))
              ENDIF
            ENDIF
            XHVST1=1._r8-XHVST
            D3385: DO M=1,jsken
              FHVST=XHVST1(ielmc)*CFOPE(instruct,M,ielmc,NZ)*EPOOLR(ielmc,N,L,NZ)
              FHVSN=XHVST1(ielmn)*CFOPE(instruct,M,ielmn,NZ)*EPOOLR(ielmn,N,L,NZ)
              FHVSP=XHVST1(ielmp)*CFOPE(instruct,M,ielmp,NZ)*EPOOLR(ielmp,N,L,NZ)
              ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+(1._r8-FFIRE)*FHVST
              ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+(1._r8-FFIRN)*FHVSN
              ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+(1._r8-FFIRP)*FHVSP
              VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
              VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE*FHVST
              VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
              VNH3F(NZ)=VNH3F(NZ)-FFIRN*FHVSN
              VN2OF(NZ)=VN2OF(NZ)-0.0_r8
              VPO4F(NZ)=VPO4F(NZ)-FFIRP*FHVSP
              CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
              TNBP=TNBP-FCH4F*FFIRE*FHVST
              DO NR=1,NRT(NZ)
                FHVST=XHVST1(ielmc)*CFOPE(icwood,M,ielmc,NZ)*(WTRT1E(ielmc,N,L,NR,NZ) &
                  +WTRT2E(ielmc,N,L,NR,NZ))*FWODR(0)
                FHVSN=XHVST1(ielmn)*CFOPE(icwood,M,ielmn,NZ)*(WTRT1E(ielmn,N,L,NR,NZ) &
                  +WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(0)
                FHVSP=XHVST1(ielmp)*CFOPE(icwood,M,ielmp,NZ)*(WTRT1E(ielmp,N,L,NR,NZ) &
                  +WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(0)
                ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+(1._r8-FFIRE)*FHVST
                ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+(1._r8-FFIRN)*FHVSN
                ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+(1._r8-FFIRP)*FHVSP
                VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE*FHVST
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
                VNH3F(NZ)=VNH3F(NZ)-FFIRN*FHVSN
                VN2OF(NZ)=VN2OF(NZ)-0.0
                VPO4F(NZ)=VPO4F(NZ)-FFIRP*FHVSP
                CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                TNBP=TNBP-FCH4F*FFIRE*FHVST
                FHVST=XHVST1(ielmc)*CFOPE(iroot,M,ielmc,NZ)*(WTRT1E(ielmc,N,L,NR,NZ) &
                  +WTRT2E(ielmc,N,L,NR,NZ))*FWODR(1)
                FHVSN=XHVST1(ielmn)*CFOPE(iroot,M,ielmn,NZ)*(WTRT1E(ielmn,N,L,NR,NZ) &
                  +WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(1)
                FHVSP=XHVST1(ielmp)*CFOPE(iroot,M,ielmp,NZ)*(WTRT1E(ielmp,N,L,NR,NZ) &
                  +WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(1)
                ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+(1._r8-FFIRE)*FHVST
                ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+(1._r8-FFIRN)*FHVSN
                ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+(1._r8-FFIRP)*FHVSP
                VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE*FHVST
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667_r8
                VNH3F(NZ)=VNH3F(NZ)-FFIRN*FHVSN
                VN2OF(NZ)=VN2OF(NZ)-0.0_r8
                VPO4F(NZ)=VPO4F(NZ)-FFIRP*FHVSP
                CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                TNBP=TNBP-FCH4F*FFIRE*FHVST
              enddo
            ENDDO D3385
!
!     RELEASE ROOT GAS CONTENTS DURING HARVESTING
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Z(NZ)=RCO2Z(NZ)-XHVST1(ielmc)*(CO2A(N,L,NZ)+CO2P(N,L,NZ))
            ROXYZ(NZ)=ROXYZ(NZ)-XHVST1(ielmc)*(OXYA(N,L,NZ)+OXYP(N,L,NZ))
            RCH4Z(NZ)=RCH4Z(NZ)-XHVST1(ielmc)*(CH4A(N,L,NZ)+CH4P(N,L,NZ))
            RN2OZ(NZ)=RN2OZ(NZ)-XHVST1(ielmc)*(Z2OA(N,L,NZ)+Z2OP(N,L,NZ))
            RNH3Z(NZ)=RNH3Z(NZ)-XHVST1(ielmc)*(ZH3A(N,L,NZ)+ZH3P(N,L,NZ))
            RH2GZ(NZ)=RH2GZ(NZ)-XHVST1(ielmc)*(H2GA(N,L,NZ)+H2GP(N,L,NZ))
            CO2A(N,L,NZ)=XHVST(ielmc)*CO2A(N,L,NZ)
            OXYA(N,L,NZ)=XHVST(ielmc)*OXYA(N,L,NZ)
            CH4A(N,L,NZ)=XHVST(ielmc)*CH4A(N,L,NZ)
            Z2OA(N,L,NZ)=XHVST(ielmc)*Z2OA(N,L,NZ)
            ZH3A(N,L,NZ)=XHVST(ielmc)*ZH3A(N,L,NZ)
            H2GA(N,L,NZ)=XHVST(ielmc)*H2GA(N,L,NZ)
            CO2P(N,L,NZ)=XHVST(ielmc)*CO2P(N,L,NZ)
            OXYP(N,L,NZ)=XHVST(ielmc)*OXYP(N,L,NZ)
            CH4P(N,L,NZ)=XHVST(ielmc)*CH4P(N,L,NZ)
            Z2OP(N,L,NZ)=XHVST(ielmc)*Z2OP(N,L,NZ)
            ZH3P(N,L,NZ)=XHVST(ielmc)*ZH3P(N,L,NZ)
            H2GP(N,L,NZ)=XHVST(ielmc)*H2GP(N,L,NZ)
!
!     REDUCE ROOT STATE VARIABLES DURING HARVESTING
!
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
            D3960: DO NR=1,NRT(NZ)
              WTRT1E(ielmc,N,L,NR,NZ)=WTRT1E(ielmc,N,L,NR,NZ)*XHVST(ielmc)
              WTRT2E(ielmc,N,L,NR,NZ)=WTRT2E(ielmc,N,L,NR,NZ)*XHVST(ielmc)
              WTRT1E(ielmn,N,L,NR,NZ)=WTRT1E(ielmn,N,L,NR,NZ)*XHVST(ielmn)
              WTRT2E(ielmn,N,L,NR,NZ)=WTRT2E(ielmn,N,L,NR,NZ)*XHVST(ielmn)
              WTRT1E(ielmp,N,L,NR,NZ)=WTRT1E(ielmp,N,L,NR,NZ)*XHVST(ielmp)
              WTRT2E(ielmp,N,L,NR,NZ)=WTRT2E(ielmp,N,L,NR,NZ)*XHVST(ielmp)
              RTWT1E(N,NR,ielmc,NZ)=RTWT1E(N,NR,ielmc,NZ)*XHVST(ielmc)
              RTWT1E(N,NR,ielmn,NZ)=RTWT1E(N,NR,ielmn,NZ)*XHVST(ielmc)
              RTWT1E(N,NR,ielmp,NZ)=RTWT1E(N,NR,ielmp,NZ)*XHVST(ielmc)
              RTLG1(N,L,NR,NZ)=RTLG1(N,L,NR,NZ)*XHVST(ielmc)
              RTLG2(N,L,NR,NZ)=RTLG2(N,L,NR,NZ)*XHVST(ielmc)
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST(ielmc)
            ENDDO D3960
            EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)*XHVST(ielmc)
            EPOOLR(ielmn,N,L,NZ)=EPOOLR(ielmn,N,L,NZ)*XHVST(ielmn)
            EPOOLR(ielmp,N,L,NZ)=EPOOLR(ielmp,N,L,NZ)*XHVST(ielmp)
            WTRTL(N,L,NZ)=WTRTL(N,L,NZ)*XHVST(ielmc)
            WTRTD(N,L,NZ)=WTRTD(N,L,NZ)*XHVST(ielmc)
            WSRTL(N,L,NZ)=WSRTL(N,L,NZ)*XHVST(ielmc)
            RTN1(N,L,NZ)=RTN1(N,L,NZ)*XHVST(ielmc)
            RTNL(N,L,NZ)=RTNL(N,L,NZ)*XHVST(ielmc)
            RTLGP(N,L,NZ)=RTLGP(N,L,NZ)*XHVST(ielmc)
            RTDNP(N,L,NZ)=RTDNP(N,L,NZ)*XHVST(ielmc)
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
              DO NE=1,npelms
                D3395: DO M=1,jsken
                  ESNC(M,NE,1,L,NZ)=ESNC(M,NE,1,L,NZ)+XHVST1(NE) &
                    *(CFOPE(iroot,M,NE,NZ)*WTNDLE(L,NE,NZ) &
                    +CFOPE(instruct,M,NE,NZ)*EPOOLN(L,NE,NZ))
                ENDDO D3395
                WTNDLE(L,NE,NZ)=WTNDLE(L,NE,NZ)*XHVST(NE)
                EPOOLN(L,NE,NZ)=EPOOLN(L,NE,NZ)*XHVST(NE)
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
        IF(ISTYP(NZ).NE.0)THEN
          D3400: DO M=1,jsken
            ESNC(M,ielmc,0,NG(NZ),NZ)=ESNC(M,ielmc,0,NG(NZ),NZ) &
              +(XHVST1(ielmc)*CFOPE(instruct,M,ielmc,NZ)*WTRVE(ielmc,NZ))*FWOODE(ielmc,0)
            ESNC(M,ielmn,0,NG(NZ),NZ)=ESNC(M,ielmn,0,NG(NZ),NZ) &
              +(XHVST1(ielmn)*CFOPE(instruct,M,ielmn,NZ)*WTRVE(ielmn,NZ))*FWOODE(ielmn,0)
            ESNC(M,ielmp,0,NG(NZ),NZ)=ESNC(M,ielmp,0,NG(NZ),NZ) &
              +(XHVST1(ielmp)*CFOPE(instruct,M,ielmp,NZ)*WTRVE(ielmp,NZ))*FWOODE(ielmp,0)
            ESNC(M,ielmc,1,NG(NZ),NZ)=ESNC(M,ielmc,1,NG(NZ),NZ) &
              +(XHVST1(ielmc)*CFOPE(instruct,M,ielmc,NZ)*WTRVE(ielmc,NZ))*FWOODE(ielmc,1)
            ESNC(M,ielmn,1,NG(NZ),NZ)=ESNC(M,ielmn,1,NG(NZ),NZ) &
              +(XHVST1(ielmn)*CFOPE(instruct,M,ielmn,NZ)*WTRVE(ielmn,NZ))*FWOODE(ielmn,1)
            ESNC(M,ielmp,1,NG(NZ),NZ)=ESNC(M,ielmp,1,NG(NZ),NZ) &
              +(XHVST1(ielmp)*CFOPE(instruct,M,ielmp,NZ)*WTRVE(ielmp,NZ))*FWOODE(ielmp,1)
          ENDDO D3400
          WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)*XHVST(ielmc)
          WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)*XHVST(ielmn)
          WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)*XHVST(ielmp)
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

end module PlantDisturbsMod
