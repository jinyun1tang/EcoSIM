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
  character(len=*),private, parameter :: mod_filename = __FILE__
! end_include_section

! disturbance variables
  real(r8) :: WTHTH0,WTHNH0,WTHPH0
  real(r8) :: WTHTH1,WTHNH1,WTHPH1
  real(r8) :: WTHTH2,WTHNH2,WTHPH2
  real(r8) :: WTHTH3,WTHNH3,WTHPH3
  real(r8) :: WTHTH4,WTHNH4,WTHPH4
  real(r8) :: WTHTR1,WTHNR1,WTHPR1
  real(r8) :: WTHTR2,WTHNR2,WTHPR2
  real(r8) :: WTHTR3,WTHNR3,WTHPR3
  real(r8) :: WTHTR4,WTHNR4,WTHPR4
  real(r8) :: WTHTX0,WTHNX0,WTHPX0
  real(r8) :: WTHTX1,WTHNX1,WTHPX1
  real(r8) :: WTHTX2,WTHNX2,WTHPX2
  real(r8) :: WTHTX3,WTHNX3,WTHPX3
  real(r8) :: WTHTX4,WTHNX4,WTHPX4
  real(r8) :: WTHTG,WTHNG,WTHPG
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
    WTSTDG   =>  plt_biom%WTSTDG  , &
    WTSTDN   =>  plt_biom%WTSTDN  , &
    WTSTDP   =>  plt_biom%WTSTDP  , &
    WTSTG    =>  plt_biom%WTSTG     &
  )

  WTHTH0=0._r8;WTHNH0=0._r8;WTHPH0=0._r8
  WTHTH1=0._r8;WTHNH1=0._r8;WTHPH1=0._r8
  WTHTH2=0._r8;WTHNH2=0._r8;WTHPH2=0._r8
  WTHTH3=0._r8;WTHNH3=0._r8;WTHPH3=0._r8
  WTHTH4=0._r8;WTHNH4=0._r8;WTHPH4=0._r8
  WTHTR1=0._r8;WTHNR1=0._r8;WTHPR1=0._r8
  WTHTR2=0._r8;WTHNR2=0._r8;WTHPR2=0._r8
  WTHTR3=0._r8;WTHNR3=0._r8;WTHPR3=0._r8
  WTHTR4=0._r8;WTHNR4=0._r8;WTHPR4=0._r8
  WTHTX0=0._r8;WTHNX0=0._r8;WTHPX0=0._r8
  WTHTX1=0._r8;WTHNX1=0._r8;WTHPX1=0._r8
  WTHTX2=0._r8;WTHNX2=0._r8;WTHPX2=0._r8
  WTHTX3=0._r8;WTHNX3=0._r8;WTHPX3=0._r8
  WTHTX4=0._r8;WTHNX4=0._r8;WTHPX4=0._r8
  WTHTG=0._r8;WTHNG=0._r8;WTHPG=0._r8
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
!     WTHTH4,WTHNH4,WTHPH4=harvested standing dead C,N,P
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
!
  IF(IHVST(NZ).GE.0)THEN
    IF(J.EQ.INT(ZNOON).AND.IHVST(NZ).NE.4 &
      .AND.IHVST(NZ).NE.6)THEN
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
      IF(WTSTG(NZ).GT.ZEROP(NZ))THEN
        WHVSTD=HVST(NZ)*THIN(NZ)*0.45/24.0 &
          *AREA3(NU)*EHVST(1,4,NZ)
        FHVST=AZMAX1(1._r8-WHVSTD/WTSTG(NZ))
        FHVSH=FHVST
      ELSE
        FHVST=1.0_r8
        FHVSH=1.0_r8
      ENDIF
    ELSE
      FHVST=1.0_r8
      FHVSH=1.0_r8
    ENDIF
    DO 6475 M=1,4
      WTHTH4=WTHTH4+(1._r8-FHVSH)*WTSTDG(M,NZ)
      WTHNH4=WTHNH4+(1._r8-FHVSH)*WTSTDN(M,NZ)
      WTHPH4=WTHPH4+(1._r8-FHVSH)*WTSTDP(M,NZ)
      WTHTX4=WTHTX4+(FHVSH-FHVST)*WTSTDG(M,NZ)
      WTHNX4=WTHNX4+(FHVSH-FHVST)*WTSTDN(M,NZ)
      WTHPX4=WTHPX4+(FHVSH-FHVST)*WTSTDP(M,NZ)
      WTSTDG(M,NZ)=FHVST*WTSTDG(M,NZ)
      WTSTDN(M,NZ)=FHVST*WTSTDN(M,NZ)
      WTSTDP(M,NZ)=FHVST*WTSTDP(M,NZ)
6475  CONTINUE
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
  associate(                            &
    WTSTDG   =>  plt_biom%WTSTDG  , &
    WTSTDN   =>  plt_biom%WTSTDN  , &
    WTSTDP   =>  plt_biom%WTSTDP  , &
    IHVST    =>  plt_distb%IHVST  , &
    FWOOD    =>  plt_allom%FWOOD  , &
    FWOODP   =>  plt_allom%FWOODP , &
    FWOODN   =>  plt_allom%FWOODN , &
    CFOPC    =>  plt_soilchem%CFOPC  , &
    CFOPN    =>  plt_soilchem%CFOPN  , &
    CFOPP    =>  plt_soilchem%CFOPP  , &
    CSNC     =>  plt_bgcr%CSNC    , &
    ZSNC     =>  plt_bgcr%ZSNC    , &
    PSNC     =>  plt_bgcr%PSNC    , &
    TPSNC    =>  plt_bgcr%TPSNC   , &
    TCSNC    =>  plt_bgcr%TCSNC   , &
    TPSN0    =>  plt_bgcr%TPSN0   , &
    TZSN0    =>  plt_bgcr%TZSN0   , &
    TCSN0    =>  plt_bgcr%TCSN0   , &
    TZSNC    =>  plt_bgcr%TZSNC   , &
    IBTYP    =>  plt_pheno%IBTYP  , &
    IGTYP    =>  plt_pheno%IGTYP    &
  )
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
!     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
!     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      DO 6375 M=1,jsken
        CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ) &
          +CFOPC(0,M,NZ)*(WTHTR0+WTHTX0) &
          +CFOPC(1,M,NZ)*(WTHTR1+WTHTX1) &
          +CFOPC(2,M,NZ)*(WTHTR2+WTHTX2)
        ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ) &
          +CFOPN(0,M,NZ)*(WTHNR0+WTHNX0) &
          +CFOPN(1,M,NZ)*(WTHNR1+WTHNX1) &
          +CFOPN(2,M,NZ)*(WTHNR2+WTHNX2)
        PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ) &
          +CFOPP(0,M,NZ)*(WTHPR0+WTHPX0) &
          +CFOPP(1,M,NZ)*(WTHPR1+WTHPX1) &
          +CFOPP(2,M,NZ)*(WTHPR2+WTHPX2)
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ) &
            +CFOPC(3,M,NZ)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
          ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ) &
            +CFOPN(3,M,NZ)*(WTHNR3+WTHNX3+WTHNR4+WTHNX4)
          PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ) &
            +CFOPP(3,M,NZ)*(WTHPR3+WTHPX3+WTHPR4+WTHPX4)
        ELSE
          WTSTDG(M,NZ)=WTSTDG(M,NZ) &
            +CFOPC(5,M,NZ)*(WTHTX3+WTHTX4)
          WTSTDN(M,NZ)=WTSTDN(M,NZ) &
            +CFOPN(5,M,NZ)*(WTHNX3+WTHNX4)
          WTSTDP(M,NZ)=WTSTDP(M,NZ) &
            +CFOPP(5,M,NZ)*(WTHPX3+WTHPX4)
          CSNC(M,0,0,NZ)=CSNC(M,0,0,NZ) &
            +CFOPC(5,M,NZ)*(WTHTR3+WTHTR4)*FWOOD(0)
          ZSNC(M,0,0,NZ)=ZSNC(M,0,0,NZ) &
            +CFOPN(5,M,NZ)*(WTHNR3+WTHNR4)*FWOODN(0)
          PSNC(M,0,0,NZ)=PSNC(M,0,0,NZ) &
            +CFOPP(5,M,NZ)*(WTHPR3+WTHPR4)*FWOODP(0)
          CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ) &
            +CFOPC(5,M,NZ)*(WTHTR3+WTHTR4)*FWOOD(1)
          ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ) &
            +CFOPN(5,M,NZ)*(WTHNR3+WTHNR4)*FWOODN(1)
          PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ) &
            +CFOPP(5,M,NZ)*(WTHPR3+WTHPR4)*FWOODP(0)
        ENDIF
6375  CONTINUE
!
!     ABOVE-GROUND LITTERFALL FROM FIRE
!
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
!     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
!     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
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
      DO 6485 M=1,jsken
        CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ) &
          +CFOPC(0,M,NZ)*(WTHTR0+WTHTX0) &
          +CFOPC(1,M,NZ)*(WTHTR1+WTHTX1) &
          +CFOPC(2,M,NZ)*(WTHTR2+WTHTX2)
        ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ) &
          +CFOPN(0,M,NZ)*WTHNL0 &
          +CFOPN(1,M,NZ)*WTHNL1 &
          +CFOPN(2,M,NZ)*WTHNL2
        PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ) &
          +CFOPP(0,M,NZ)*WTHPL0 &
          +CFOPP(1,M,NZ)*WTHPL1 &
          +CFOPP(2,M,NZ)*WTHPL2
        ZSNC(4,1,0,NZ)=ZSNC(4,1,0,NZ) &
          +CFOPN(0,M,NZ)*(WTHNR0+WTHNX0-WTHNL0) &
          +CFOPN(1,M,NZ)*(WTHNR1+WTHNX1-WTHNL1) &
          +CFOPN(2,M,NZ)*(WTHNR2+WTHNX2-WTHNL2)
        PSNC(4,1,0,NZ)=PSNC(4,1,0,NZ) &
          +CFOPP(0,M,NZ)*(WTHPR0+WTHPX0-WTHPL0) &
          +CFOPP(1,M,NZ)*(WTHPR1+WTHPX1-WTHPL1) &
          +CFOPP(2,M,NZ)*(WTHPR2+WTHPX2-WTHPL2)
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ)+CFOPC(3,M,NZ)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
          ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ)+CFOPN(3,M,NZ)*(WTHNL3+WTHNL4)
          PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ)+CFOPP(3,M,NZ)*(WTHPL3+WTHPL4)
          ZSNC(4,1,0,NZ)=ZSNC(4,1,0,NZ)+CFOPN(3,M,NZ)*(WTHNR3+WTHNX3-WTHNL3+WTHNR4+WTHNX4-WTHNL4)
          PSNC(4,1,0,NZ)=PSNC(4,1,0,NZ)+CFOPP(3,M,NZ)*(WTHPR3+WTHPX3-WTHPL3+WTHPR4+WTHPX4-WTHPL4)
        ELSE
          WTSTDG(M,NZ)=WTSTDG(M,NZ)+CFOPC(5,M,NZ)*(WTHTR3+WTHTX3)
          WTSTDN(M,NZ)=WTSTDN(M,NZ)+CFOPN(5,M,NZ)*WTHNL3
          WTSTDP(M,NZ)=WTSTDP(M,NZ)+CFOPP(5,M,NZ)*WTHPL3
          CSNC(M,0,0,NZ)=CSNC(M,0,0,NZ)*CFOPC(3,M,NZ)*(WTHTR4+WTHTX4)*FWOOD(0)
          ZSNC(M,0,0,NZ)=ZSNC(M,0,0,NZ)+CFOPN(3,M,NZ)*WTHNL4*FWOODN(0)
          PSNC(M,0,0,NZ)=PSNC(M,0,0,NZ)+CFOPP(3,M,NZ)*WTHPL4*FWOODP(0)
          ZSNC(4,0,0,NZ)=ZSNC(4,0,0,NZ)+CFOPN(5,M,NZ)*(WTHNR3+WTHNX3-WTHNL3 &
            +WTHNR4+WTHNX4-WTHNL4)*FWOODN(0)
          PSNC(4,0,0,NZ)=PSNC(4,0,0,NZ)+CFOPP(5,M,NZ)*(WTHPR3+WTHPX3-WTHPL3 &
            +WTHPR4+WTHPX4-WTHPL4)*FWOODP(0)
          CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ)+CFOPC(3,M,NZ)*(WTHTR4+WTHTX4)*FWOOD(1)
          ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ)+CFOPN(3,M,NZ)*WTHNL4*FWOODN(1)
          PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ)+CFOPP(3,M,NZ)*WTHPL4*FWOODP(1)
          ZSNC(4,1,0,NZ)=ZSNC(4,1,0,NZ)+CFOPN(5,M,NZ)*(WTHNR3+WTHNX3-WTHNL3 &
            +WTHNR4+WTHNX4-WTHNL4)*FWOODN(1)
          PSNC(4,1,0,NZ)=PSNC(4,1,0,NZ)+CFOPP(5,M,NZ)*(WTHPR3+WTHPX3-WTHPL3 &
            +WTHPR4+WTHPX4-WTHPL4)*FWOODP(1)
        ENDIF
6485  CONTINUE
    ENDIF
  ELSE
!
!     ABOVE-GROUND LITTERFALL FROM GRAZING
!
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!
    TCSNC(NZ)=TCSNC(NZ)+WTHTRT+WTHTXT
    TZSNC(NZ)=TZSNC(NZ)+WTHNRT+WTHNXT
    TPSNC(NZ)=TPSNC(NZ)+WTHPRT+WTHPXT
    TCSN0(NZ)=TCSN0(NZ)+WTHTRT+WTHTXT
    TZSN0(NZ)=TZSNC(NZ)+WTHNRT+WTHNXT
    TPSN0(NZ)=TPSNC(NZ)+WTHPRT+WTHPXT
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
    HVSTC    =>  plt_distb%HVSTC  , &
    HVSTN    =>  plt_distb%HVSTN  , &
    HVSTP    =>  plt_distb%HVSTP  , &
    XHVSTC   =>  plt_distb%XHVSTC , &
    XHVSTN   =>  plt_distb%XHVSTN , &
    XHVSTP   =>  plt_distb%XHVSTP , &
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
    WTRVN    =>  plt_biom%WTRVN   , &
    WTRVC    =>  plt_biom%WTRVC   , &
    WTRVP    =>  plt_biom%WTRVP     &
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
  WTHTHT=WTHTH0+WTHTH1+WTHTH2+WTHTH3+WTHTH4
  WTHTRT=WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4
  WTHNHT=WTHNH0+WTHNH1+WTHNH2+WTHNH3+WTHNH4
  WTHNRT=WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4
  WTHPHT=WTHPH0+WTHPH1+WTHPH2+WTHPH3+WTHPH4
  WTHPRT=WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4
  WTHTXT=WTHTX0+WTHTX1+WTHTX2+WTHTX3+WTHTX4
  WTHNXT=WTHNX0+WTHNX1+WTHNX2+WTHNX3+WTHNX4
  WTHPXT=WTHPX0+WTHPX1+WTHPX2+WTHPX3+WTHPX4

  IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
    IF(IHVST(NZ).NE.5)THEN
      IF(JHVST(NZ).NE.2)THEN
        HVSTC(NZ)=HVSTC(NZ)+WTHTHT-WTHTRT
        HVSTN(NZ)=HVSTN(NZ)+WTHNHT-WTHNRT
        HVSTP(NZ)=HVSTP(NZ)+WTHPHT-WTHPRT
        TNBP=TNBP+WTHTRT-WTHTHT
        XHVSTC=XHVSTC+WTHTHT-WTHTRT
        XHVSTN=XHVSTN+WTHNHT-WTHNRT
        XHVSTP=XHVSTP+WTHPHT-WTHPRT
      ELSE
        WTRVC(NZ)=WTRVC(NZ)+WTHTHT-WTHTRT
        WTRVN(NZ)=WTRVN(NZ)+WTHNHT-WTHNRT
        WTRVP(NZ)=WTRVP(NZ)+WTHPHT-WTHPRT
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
    HVSTC(NZ)=HVSTC(NZ)+GY*(WTHTHT-WTHTRT)
    HVSTN(NZ)=HVSTN(NZ)+WTHNHT-WTHNRT
    HVSTP(NZ)=HVSTP(NZ)+WTHPHT-WTHPRT
    TCO2T(NZ)=TCO2T(NZ)-GZ*(WTHTHT-WTHTRT)
    TCO2A(NZ)=TCO2A(NZ)-GZ*(WTHTHT-WTHTRT)
!     TNBP=TNBP+GY*(WTHTRT-WTHTHT)
!     CNET(NZ)=CNET(NZ)+GZ*(WTHTRT-WTHTHT)
    XHVSTC=XHVSTC+GY*(WTHTHT-WTHTRT)
    XHVSTN=XHVSTN+WTHNHT-WTHNRT
    XHVSTP=XHVSTP+WTHPHT-WTHPRT
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
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!     WTHTH1,WTHNH1,WTHPH1=leaf C,N,P removed
!     WTHTH2,WTHNH2,WTHPH2=fine,non-leaf C,N,P removed
!     WTHTH3,WTHNH3,WTHPH3=woody C,N,P removed
!     WTHTH4,WTHNH4,WTHPH4=standing dead C,N,P removed
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  IF(IHVST(NZ).EQ.0)THEN
    WTHTR0=WTHTH0*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVST(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EHVST(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVST(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVST(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EHVST(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVST(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVST(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EHVST(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVST(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVST(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EHVST(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVST(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVST(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EHVST(2,4,NZ))
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ).EQ.1)THEN
    WTHTR0=WTHTH0
    WTHNR0=WTHNH0
    WTHPR0=WTHPH0
    WTHTR1=WTHTH1
    WTHNR1=WTHNH1
    WTHPR1=WTHPH1
    WTHTR2=WTHTH2-WTHTG*EHVST(2,2,NZ)
    WTHNR2=WTHNH2-WTHNG*EHVST(2,2,NZ)
    WTHPR2=WTHPH2-WTHPG*EHVST(2,2,NZ)
    WTHTR3=WTHTH3
    WTHNR3=WTHNH3
    WTHPR3=WTHPH3
    WTHTR4=WTHTH4
    WTHNR4=WTHNH4
    WTHPR4=WTHPH4
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ).EQ.2)THEN
    WTHTR0=WTHTH0*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVST(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EHVST(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVST(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVST(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EHVST(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVST(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVST(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EHVST(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVST(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVST(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EHVST(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVST(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVST(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EHVST(2,4,NZ))
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(IHVST(NZ).EQ.3)THEN
    WTHTR0=WTHTH0*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVST(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EHVST(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVST(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVST(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EHVST(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVST(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVST(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EHVST(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVST(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVST(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EHVST(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVST(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVST(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EHVST(2,4,NZ))
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
    WTHTR0=WTHTH0*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHPR0=WTHPH0*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHTR1=WTHTH1*(1._r8-EHVST(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHPR1=WTHPH1*(1._r8-EHVST(2,1,NZ)*0.5)
    WTHTR2=WTHTH2*(1._r8-EHVST(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVST(2,2,NZ)*0.5)
    WTHPR2=WTHPH2*(1._r8-EHVST(2,2,NZ)*0.5)
    WTHTR3=WTHTH3*(1._r8-EHVST(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVST(2,3,NZ)*0.5)
    WTHPR3=WTHPH3*(1._r8-EHVST(2,3,NZ)*0.5)
    WTHTR4=WTHTH4*(1._r8-EHVST(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVST(2,4,NZ)*0.5)
    WTHPR4=WTHPH4*(1._r8-EHVST(2,4,NZ)*0.5)
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
    FERT(17)=FERT(17)+(WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4)/AREA3(NU)
    FERT(18)=FERT(18)+(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA3(NU)*0.5_r8
    FERT(3)=FERT(3)+(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA3(NU)*0.5_r8
    FERT(19)=FERT(19)+(WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4)/AREA3(NU)
    IYTYP=3
!
!     REMOVALS BY FIRE
!
!     EFIRE=combustion  of N,P relative to C
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!     WTHTH1,WTHNH1,WTHPH1=leaf C,N,P removed
!     WTHTH2,WTHNH2,WTHPH2=fine,non-leaf C,N,P removed
!     WTHTH3,WTHNH3,WTHPH3=woody C,N,P removed
!     WTHTH4,WTHNH4,WTHPH4=standing dead C,N,P removed
!     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
!     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
!     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
!     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
!     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
!     WTHTL0,WTHNL0,WTHPL0=nonstructural C,N,P removed from ecosystem
!     WTHTL1,WTHNL1,WTHPL1=leaf C,N,P removed from ecosystem
!     WTHTL2,WTHNL2,WTHPL2=fine,non-leaf C,N,P removed from ecosystem
!     WTHTL3,WTHNL3,WTHPL3=woody C,N,P removed from ecosystem
!     WTHTL4,WTHNL4,WTHPL4=standing dead C,N,P removed from ecosystem
!
  ELSEIF(IHVST(NZ).EQ.5)THEN
    WTHTR0=WTHTH0*(1._r8-EHVST(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,1,NZ))
    WTHNL0=WTHNH0*(1._r8-EHVST(2,1,NZ))
    WTHPL0=WTHPH0*(1._r8-EHVST(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVST(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,1,NZ))
    WTHNL1=WTHNH1*(1._r8-EHVST(2,1,NZ))
    WTHPL1=WTHPH1*(1._r8-EHVST(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVST(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,2,NZ))
    WTHNL2=WTHNH2*(1._r8-EHVST(2,2,NZ))
    WTHPL2=WTHPH2*(1._r8-EHVST(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVST(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,3,NZ))
    WTHNL3=WTHNH3*(1._r8-EHVST(2,3,NZ))
    WTHPL3=WTHPH3*(1._r8-EHVST(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVST(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EFIRE(1,IHVST(NZ))*EHVST(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EFIRE(2,IHVST(NZ))*EHVST(2,4,NZ))
    WTHNL4=WTHNH4*(1._r8-EHVST(2,4,NZ))
    WTHPL4=WTHPH4*(1._r8-EHVST(2,4,NZ))
  ENDIF
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByTillage(I,J,NZ,CPOOLK)
  use SurfLitterDataType, only : XCORP
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(inout) :: CPOOLK(10,JP1)
  integer :: L,K,M,N,NR,NB
  real(r8) :: XHVST
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
    CFOPC    =>  plt_soilchem%CFOPC  , &
    CFOPN    =>  plt_soilchem%CFOPN  , &
    CFOPP    =>  plt_soilchem%CFOPP  , &
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
    CPOOLR   =>  plt_biom%CPOOLR     , &
    ZPOOLR   =>  plt_biom%ZPOOLR     , &
    PPOOLR   =>  plt_biom%PPOOLR     , &
    WSRTL    =>  plt_biom%WSRTL      , &
    WTRTD    =>  plt_biom%WTRTD      , &
    WTRTL    =>  plt_biom%WTRTL      , &
    WGLFLP   =>  plt_biom%WGLFLP     , &
    WTLFBP   =>  plt_biom%WTLFBP     , &
    WTLFB    =>  plt_biom%WTLFB      , &
    WTGRBN   =>  plt_biom%WTGRBN     , &
    WTGRB    =>  plt_biom%WTGRB      , &
    WTEARB   =>  plt_biom%WTEARB     , &
    CPOLNB   =>  plt_biom%CPOLNB     , &
    PPOOL    =>  plt_biom%PPOOL      , &
    ZPOOL    =>  plt_biom%ZPOOL      , &
    CPOOL    =>  plt_biom%CPOOL      , &
    WTLFBN   =>  plt_biom%WTLFBN     , &
    WTSHTN   =>  plt_biom%WTSHTN     , &
    WTHSKB   =>  plt_biom%WTHSKB     , &
    WTRSVB   =>  plt_biom%WTRSVB     , &
    WTNDB    =>  plt_biom%WTNDB      , &
    WTSHTB   =>  plt_biom%WTSHTB     , &
    PPOLNB   =>  plt_biom%PPOLNB     , &
    ZPOLNB   =>  plt_biom%ZPOLNB     , &
    WTSTKB   =>  plt_biom%WTSTKB     , &
    WTEABN   =>  plt_biom%WTEABN     , &
    WTSHBN   =>  plt_biom%WTSHBN     , &
    WTSTBN   =>  plt_biom%WTSTBN     , &
    WTHSBN   =>  plt_biom%WTHSBN     , &
    WTNDBN   =>  plt_biom%WTNDBN     , &
    WTRSBN   =>  plt_biom%WTRSBN     , &
    WTSHTP   =>  plt_biom%WTSHTP     , &
    WGLFLN   =>  plt_biom%WGLFLN     , &
    WTSHBP   =>  plt_biom%WTSHBP     , &
    WTNDBP   =>  plt_biom%WTNDBP     , &
    WGLFL    =>  plt_biom%WGLFL      , &
    WGSHP    =>  plt_biom%WGSHP      , &
    WTSHEB   =>  plt_biom%WTSHEB     , &
    WTEABP   =>  plt_biom%WTEABP     , &
    WTHSBP   =>  plt_biom%WTHSBP     , &
    WTRSBP   =>  plt_biom%WTRSBP     , &
    WTSTBP   =>  plt_biom%WTSTBP     , &
    WGLF     =>  plt_biom%WGLF       , &
    WTSTXN   =>  plt_biom%WTSTXN     , &
    WTLSB    =>  plt_biom%WTLSB      , &
    WTGRBP   =>  plt_biom%WTGRBP     , &
    WTSTXB   =>  plt_biom%WTSTXB     , &
    WVSTKB   =>  plt_biom%WVSTKB     , &
    WTSTXP   =>  plt_biom%WTSTXP     , &
    WGLFP    =>  plt_biom%WGLFP      , &
    WSSHE    =>  plt_biom%WSSHE      , &
    WGSHE    =>  plt_biom%WGSHE      , &
    WTRT1    =>  plt_biom%WTRT1      , &
    WTRT1N   =>  plt_biom%WTRT1N     , &
    WTRT1P   =>  plt_biom%WTRT1P     , &
    WSLF     =>  plt_biom%WSLF       , &
    WGSHN    =>  plt_biom%WGSHN      , &
    WGNODE   =>  plt_biom%WGNODE     , &
    WGLFN    =>  plt_biom%WGLFN      , &
    WGNODN   =>  plt_biom%WGNODN     , &
    WGNODP   =>  plt_biom%WGNODP     , &
    WVSTK    =>  plt_biom%WVSTK      , &
    RTWT1    =>  plt_biom%RTWT1      , &
    RTWT1N   =>  plt_biom%RTWT1N     , &
    RTWT1P   =>  plt_biom%RTWT1P     , &
    WTRVC    =>  plt_biom%WTRVC      , &
    WTRVN    =>  plt_biom%WTRVN      , &
    WTRVP    =>  plt_biom%WTRVP      , &
    WTLS     =>  plt_biom%WTLS       , &
    WTRT2    =>  plt_biom%WTRT2      , &
    WTRT2N   =>  plt_biom%WTRT2N     , &
    WTRT2P   =>  plt_biom%WTRT2P     , &
    CPOOLN   =>  plt_biom%CPOOLN     , &
    ZPOOLN   =>  plt_biom%ZPOOLN     , &
    PPOOLN   =>  plt_biom%PPOOLN     , &
    WTNDL    =>  plt_biom%WTNDL      , &
    WTNDLN   =>  plt_biom%WTNDLN     , &
    WTNDLP   =>  plt_biom%WTNDLP     , &
    GRWTB    =>  plt_allom%GRWTB     , &
    FWOODP   =>  plt_allom%FWOODP    , &
    FWOODN   =>  plt_allom%FWOODN    , &
    FWOOD    =>  plt_allom%FWOOD     , &
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
    CSNC     =>  plt_bgcr%CSNC       , &
    ZSNC     =>  plt_bgcr%ZSNC       , &
    PSNC     =>  plt_bgcr%PSNC       , &
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
        DO 8975 NB=1,NBR(NZ)
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
            DO 6380 M=1,jsken
              CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ)+(1._r8-XHVST) &
                *(CFOPC(0,M,NZ)*(CPOOL(NB,NZ)+CPOLNB(NB,NZ) &
                +CPOOLK(NB,NZ)+WTRSVB(NB,NZ)) &
                +CFOPC(1,M,NZ)*(WTLFB(NB,NZ)*FWODB(1) &
                +WTNDB(NB,NZ)) &
                +CFOPC(2,M,NZ)*(WTSHEB(NB,NZ)*FWODB(1) &
                +WTHSKB(NB,NZ)+WTEARB(NB,NZ)))
              CSNC(M,0,0,NZ)=CSNC(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPC(5,M,NZ)*(WTLFB(NB,NZ)*FWODB(0) &
                +WTSHEB(NB,NZ)*FWODB(0))
              ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ)+(1._r8-XHVST) &
                *(CFOPN(0,M,NZ)*(ZPOOL(NB,NZ)+ZPOLNB(NB,NZ) &
                +WTRSBN(NB,NZ)) &
                +CFOPN(1,M,NZ)*(WTLFBN(NB,NZ)*FWODLN(1) &
                +WTNDBN(NB,NZ)) &
                +CFOPN(2,M,NZ)*(WTSHBN(NB,NZ)*FWODSN(1) &
                +WTHSBN(NB,NZ)+WTEABN(NB,NZ)))
              ZSNC(M,0,0,NZ)=ZSNC(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPN(5,M,NZ)*(WTLFBN(NB,NZ)*FWODLN(0) &
                +WTSHBN(NB,NZ)*FWODSN(0))
              PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ)+(1._r8-XHVST) &
                *(CFOPP(0,M,NZ)*(PPOOL(NB,NZ)+PPOLNB(NB,NZ) &
                +WTRSBP(NB,NZ)) &
                +CFOPP(1,M,NZ)*(WTLFBP(NB,NZ)*FWODLP(1) &
                +WTNDBP(NB,NZ)) &
                +CFOPP(2,M,NZ)*(WTSHBP(NB,NZ)*FWODSP(1) &
                +WTHSBP(NB,NZ)+WTEABP(NB,NZ)))
              PSNC(M,0,0,NZ)=PSNC(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPP(5,M,NZ)*(WTLFBP(NB,NZ)*FWODLP(0) &
                +WTSHBP(NB,NZ)*FWODSP(0))
              IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
                WTRVC(NZ)=WTRVC(NZ)+(1._r8-XHVST) &
                  *CFOPC(2,M,NZ)*WTGRB(NB,NZ)
                WTRVN(NZ)=WTRVN(NZ)+(1._r8-XHVST) &
                  *CFOPN(2,M,NZ)*WTGRBN(NB,NZ)
                WTRVP(NZ)=WTRVP(NZ)+(1._r8-XHVST) &
                  *CFOPP(2,M,NZ)*WTGRBP(NB,NZ)
              ELSE
                CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ)+(1._r8-XHVST) &
                  *CFOPC(2,M,NZ)*WTGRB(NB,NZ)
                ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ)+(1._r8-XHVST) &
                  *CFOPN(2,M,NZ)*WTGRBN(NB,NZ)
                PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ)+(1._r8-XHVST) &
                  *CFOPP(2,M,NZ)*WTGRBP(NB,NZ)
              ENDIF
              CSNC(M,0,0,NZ)=CSNC(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPC(5,M,NZ)*WTSTKB(NB,NZ)*FWOOD(0)
              ZSNC(M,0,0,NZ)=ZSNC(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPN(5,M,NZ)*WTSTBN(NB,NZ)*FWOODN(0)
              PSNC(M,0,0,NZ)=PSNC(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPP(5,M,NZ)*WTSTBP(NB,NZ)*FWOODP(0)
              CSNC(M,1,0,NZ)=CSNC(M,1,0,NZ)+(1._r8-XHVST) &
                *CFOPC(3,M,NZ)*WTSTKB(NB,NZ)*FWOOD(1)
              ZSNC(M,1,0,NZ)=ZSNC(M,1,0,NZ)+(1._r8-XHVST) &
                *CFOPN(3,M,NZ)*WTSTBN(NB,NZ)*FWOODN(1)
              PSNC(M,1,0,NZ)=PSNC(M,1,0,NZ)+(1._r8-XHVST) &
                *CFOPP(3,M,NZ)*WTSTBP(NB,NZ)*FWOODP(1)
6380        CONTINUE
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!
            CPOOL(NB,NZ)=CPOOL(NB,NZ)*XHVST
            CPOOLK(NB,NZ)=CPOOLK(NB,NZ)*XHVST
            ZPOOL(NB,NZ)=ZPOOL(NB,NZ)*XHVST
            PPOOL(NB,NZ)=PPOOL(NB,NZ)*XHVST
            CPOLNB(NB,NZ)=CPOLNB(NB,NZ)*XHVST
            ZPOLNB(NB,NZ)=ZPOLNB(NB,NZ)*XHVST
            PPOLNB(NB,NZ)=PPOLNB(NB,NZ)*XHVST
            WTSHTB(NB,NZ)=WTSHTB(NB,NZ)*XHVST
            WTLFB(NB,NZ)=WTLFB(NB,NZ)*XHVST
            WTNDB(NB,NZ)=WTNDB(NB,NZ)*XHVST
            WTSHEB(NB,NZ)=WTSHEB(NB,NZ)*XHVST
            WTSTKB(NB,NZ)=WTSTKB(NB,NZ)*XHVST
            WVSTKB(NB,NZ)=WVSTKB(NB,NZ)*XHVST
            WTRSVB(NB,NZ)=WTRSVB(NB,NZ)*XHVST
            WTHSKB(NB,NZ)=WTHSKB(NB,NZ)*XHVST
            WTEARB(NB,NZ)=WTEARB(NB,NZ)*XHVST
            WTGRB(NB,NZ)=WTGRB(NB,NZ)*XHVST
            WTSHTN(NB,NZ)=WTSHTN(NB,NZ)*XHVST
            WTLFBN(NB,NZ)=WTLFBN(NB,NZ)*XHVST
            WTNDBN(NB,NZ)=WTNDBN(NB,NZ)*XHVST
            WTSHBN(NB,NZ)=WTSHBN(NB,NZ)*XHVST
            WTSTBN(NB,NZ)=WTSTBN(NB,NZ)*XHVST
            WTRSBN(NB,NZ)=WTRSBN(NB,NZ)*XHVST
            WTHSBN(NB,NZ)=WTHSBN(NB,NZ)*XHVST
            WTEABN(NB,NZ)=WTEABN(NB,NZ)*XHVST
            WTGRBN(NB,NZ)=WTGRBN(NB,NZ)*XHVST
            WTSHTP(NB,NZ)=WTSHTP(NB,NZ)*XHVST
            WTLFBP(NB,NZ)=WTLFBP(NB,NZ)*XHVST
            WTNDBP(NB,NZ)=WTNDBP(NB,NZ)*XHVST
            WTSHBP(NB,NZ)=WTSHBP(NB,NZ)*XHVST
            WTSTBP(NB,NZ)=WTSTBP(NB,NZ)*XHVST
            WTRSBP(NB,NZ)=WTRSBP(NB,NZ)*XHVST
            WTHSBP(NB,NZ)=WTHSBP(NB,NZ)*XHVST
            WTEABP(NB,NZ)=WTEABP(NB,NZ)*XHVST
            WTGRBP(NB,NZ)=WTGRBP(NB,NZ)*XHVST
            GRNXB(NB,NZ)=GRNXB(NB,NZ)*XHVST
            GRNOB(NB,NZ)=GRNOB(NB,NZ)*XHVST
            GRWTB(NB,NZ)=GRWTB(NB,NZ)*XHVST
            ARLFB(NB,NZ)=ARLFB(NB,NZ)*XHVST
            WTLSB(NB,NZ)=AZMAX1(WTLFB(NB,NZ)+WTSHEB(NB,NZ))
            WTLS(NZ)=WTLS(NZ)+WTLSB(NB,NZ)
            WTSTXB(NB,NZ)=WTSTXB(NB,NZ)*XHVST
            WTSTXN(NB,NZ)=WTSTXN(NB,NZ)*XHVST
            WTSTXP(NB,NZ)=WTSTXP(NB,NZ)*XHVST
            WVSTK(NZ)=WVSTK(NZ)+WVSTKB(NB,NZ)
            DO 8970 K=0,JNODS1
              IF(K.NE.0)THEN
                CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ)*XHVST
                CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ)*XHVST
                CO2B(K,NB,NZ)=CO2B(K,NB,NZ)*XHVST
                HCOB(K,NB,NZ)=HCOB(K,NB,NZ)*XHVST
              ENDIF
              ARLF1(K,NB,NZ)=ARLF1(K,NB,NZ)*XHVST
              WGLF(K,NB,NZ)=WGLF(K,NB,NZ)*XHVST
              WSLF(K,NB,NZ)=WSLF(K,NB,NZ)*XHVST
!     HTSHE(K,NB,NZ)=HTSHE(K,NB,NZ)*XHVST
              WGSHE(K,NB,NZ)=WGSHE(K,NB,NZ)*XHVST
              WSSHE(K,NB,NZ)=WSSHE(K,NB,NZ)*XHVST
!     HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ)*XHVST
!     HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ)*XHVST
              WGNODE(K,NB,NZ)=WGNODE(K,NB,NZ)*XHVST
              WGLFN(K,NB,NZ)=WGLFN(K,NB,NZ)*XHVST
              WGSHN(K,NB,NZ)=WGSHN(K,NB,NZ)*XHVST
              WGNODN(K,NB,NZ)=WGNODN(K,NB,NZ)*XHVST
              WGLFP(K,NB,NZ)=WGLFP(K,NB,NZ)*XHVST
              WGSHP(K,NB,NZ)=WGSHP(K,NB,NZ)*XHVST
              WGNODP(K,NB,NZ)=WGNODP(K,NB,NZ)*XHVST
              DO 8965 L=1,JC1
                ARLFL(L,K,NB,NZ)=ARLFL(L,K,NB,NZ)*XHVST
                WGLFL(L,K,NB,NZ)=WGLFL(L,K,NB,NZ)*XHVST
                WGLFLN(L,K,NB,NZ)=WGLFLN(L,K,NB,NZ)*XHVST
                WGLFLP(L,K,NB,NZ)=WGLFLP(L,K,NB,NZ)*XHVST
8965          CONTINUE
8970        CONTINUE
          ENDIF
8975    CONTINUE
!
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        VOLWPX=VOLWP(NZ)
        WVPLT=AZMAX1(WTLS(NZ)+WVSTK(NZ))
        APSILT=ABS(PSILT(NZ))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWP(NZ)=1.0E-06_r8*WVPLT/FDM
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
        DO 8985 N=1,MY(NZ)
          DO 8980 L=NU,NJ
            DO 6385 M=1,jsken
              CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-XHVST) &
                *CFOPC(0,M,NZ)*CPOOLR(N,L,NZ)
              ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-XHVST) &
                *CFOPN(0,M,NZ)*ZPOOLR(N,L,NZ)
              PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-XHVST) &
                *CFOPP(0,M,NZ)*PPOOLR(N,L,NZ)
              DO NR=1,NRT(NZ)
                CSNC(M,0,L,NZ)=CSNC(M,0,L,NZ)+(1._r8-XHVST) &
                  *CFOPC(5,M,NZ)*(WTRT1(N,L,NR,NZ) &
                  +WTRT2(N,L,NR,NZ))*FWODR(0)
                ZSNC(M,0,L,NZ)=ZSNC(M,0,L,NZ)+(1._r8-XHVST) &
                  *CFOPN(5,M,NZ)*(WTRT1N(N,L,NR,NZ) &
                  +WTRT2N(N,L,NR,NZ))*FWODRN(0)
                PSNC(M,0,L,NZ)=PSNC(M,0,L,NZ)+(1._r8-XHVST) &
                  *CFOPP(5,M,NZ)*(WTRT1P(N,L,NR,NZ) &
                  +WTRT2P(N,L,NR,NZ))*FWODRP(0)
                CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *CFOPC(4,M,NZ)*(WTRT1(N,L,NR,NZ) &
                  +WTRT2(N,L,NR,NZ))*FWODR(1)
                ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *CFOPN(4,M,NZ)*(WTRT1N(N,L,NR,NZ) &
                  +WTRT2N(N,L,NR,NZ))*FWODRN(1)
                PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *CFOPP(4,M,NZ)*(WTRT1P(N,L,NR,NZ) &
                  +WTRT2P(N,L,NR,NZ))*FWODRP(1)
              ENDDO
6385        CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Z(NZ)=RCO2Z(NZ)-(1._r8-XHVST) &
              *(CO2A(N,L,NZ)+CO2P(N,L,NZ))
            ROXYZ(NZ)=ROXYZ(NZ)-(1._r8-XHVST) &
              *(OXYA(N,L,NZ)+OXYP(N,L,NZ))
            RCH4Z(NZ)=RCH4Z(NZ)-(1._r8-XHVST) &
              *(CH4A(N,L,NZ)+CH4P(N,L,NZ))
            RN2OZ(NZ)=RN2OZ(NZ)-(1._r8-XHVST) &
              *(Z2OA(N,L,NZ)+Z2OP(N,L,NZ))
            RNH3Z(NZ)=RNH3Z(NZ)-(1._r8-XHVST) &
              *(ZH3A(N,L,NZ)+ZH3P(N,L,NZ))
            RH2GZ(NZ)=RH2GZ(NZ)-(1._r8-XHVST) &
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
            DO 8960 NR=1,NRT(NZ)
              WTRT1(N,L,NR,NZ)=WTRT1(N,L,NR,NZ)*XHVST
              WTRT2(N,L,NR,NZ)=WTRT2(N,L,NR,NZ)*XHVST
              WTRT1N(N,L,NR,NZ)=WTRT1N(N,L,NR,NZ)*XHVST
              WTRT2N(N,L,NR,NZ)=WTRT2N(N,L,NR,NZ)*XHVST
              WTRT1P(N,L,NR,NZ)=WTRT1P(N,L,NR,NZ)*XHVST
              WTRT2P(N,L,NR,NZ)=WTRT2P(N,L,NR,NZ)*XHVST
              RTWT1(N,NR,NZ)=RTWT1(N,NR,NZ)*XHVST
              RTWT1N(N,NR,NZ)=RTWT1N(N,NR,NZ)*XHVST
              RTWT1P(N,NR,NZ)=RTWT1P(N,NR,NZ)*XHVST
              RTLG1(N,L,NR,NZ)=RTLG1(N,L,NR,NZ)*XHVST
              RTLG2(N,L,NR,NZ)=RTLG2(N,L,NR,NZ)*XHVST
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST
8960        CONTINUE
            CPOOLR(N,L,NZ)=CPOOLR(N,L,NZ)*XHVST
            ZPOOLR(N,L,NZ)=ZPOOLR(N,L,NZ)*XHVST
            PPOOLR(N,L,NZ)=PPOOLR(N,L,NZ)*XHVST
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
              DO 6395 M=1,jsken
                CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPC(4,M,NZ)*WTNDL(L,NZ) &
                  +CFOPC(0,M,NZ)*CPOOLN(L,NZ))
                ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPN(4,M,NZ)*WTNDLN(L,NZ) &
                  +CFOPN(0,M,NZ)*ZPOOLN(L,NZ))
                PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPP(4,M,NZ)*WTNDLP(L,NZ) &
                  +CFOPP(0,M,NZ)*PPOOLN(L,NZ))
6395          CONTINUE
              WTNDL(L,NZ)=WTNDL(L,NZ)*XHVST
              WTNDLN(L,NZ)=WTNDLN(L,NZ)*XHVST
              WTNDLP(L,NZ)=WTNDLP(L,NZ)*XHVST
              CPOOLN(L,NZ)=CPOOLN(L,NZ)*XHVST
              ZPOOLN(L,NZ)=ZPOOLN(L,NZ)*XHVST
              PPOOLN(L,NZ)=PPOOLN(L,NZ)*XHVST
            ENDIF
8980      CONTINUE
8985    CONTINUE
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
        DO 6400 M=1,jsken
          CSNC(M,0,NG(NZ),NZ)=CSNC(M,0,NG(NZ),NZ) &
            +((1._r8-XHVST)*CFOPC(0,M,NZ)*WTRVC(NZ))*FWOOD(0)
          ZSNC(M,0,NG(NZ),NZ)=ZSNC(M,0,NG(NZ),NZ) &
            +((1._r8-XHVST)*CFOPN(0,M,NZ)*WTRVN(NZ))*FWOODN(0)
          PSNC(M,0,NG(NZ),NZ)=PSNC(M,0,NG(NZ),NZ) &
            +((1._r8-XHVST)*CFOPP(0,M,NZ)*WTRVP(NZ))*FWOODP(0)
          CSNC(M,1,NG(NZ),NZ)=CSNC(M,1,NG(NZ),NZ) &
            +((1._r8-XHVST)*CFOPC(0,M,NZ)*WTRVC(NZ))*FWOOD(1)
          ZSNC(M,1,NG(NZ),NZ)=ZSNC(M,1,NG(NZ),NZ) &
            +((1._r8-XHVST)*CFOPN(0,M,NZ)*WTRVN(NZ))*FWOODN(1)
          PSNC(M,1,NG(NZ),NZ)=PSNC(M,1,NG(NZ),NZ) &
            +((1._r8-XHVST)*CFOPP(0,M,NZ)*WTRVP(NZ))*FWOODP(1)
6400    CONTINUE
        WTRVC(NZ)=WTRVC(NZ)*XHVST
        WTRVN(NZ)=WTRVN(NZ)*XHVST
        WTRVP(NZ)=WTRVP(NZ)*XHVST
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
  integer :: L,K,M,NR,N,NB,NBX
  real(r8):: ZPOOLG,ZPOLNG,ZPOOLX
  real(r8) :: ZPOLNX,XHVSN,XHVSP,XHVST
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
    CPOOLR   => plt_biom%CPOOLR    , &
    ZPOOLR   => plt_biom%ZPOOLR    , &
    PPOOLR   => plt_biom%PPOOLR    , &
    WSRTL    => plt_biom%WSRTL     , &
    WTRTD    => plt_biom%WTRTD     , &
    WTRTL    => plt_biom%WTRTL     , &
    WTRVC    => plt_biom%WTRVC     , &
    WVSTK    => plt_biom%WVSTK     , &
    WTLS     => plt_biom%WTLS      , &
    WTRVP    => plt_biom%WTRVP     , &
    WTRVN    => plt_biom%WTRVN     , &
    WVSTKB   => plt_biom%WVSTKB    , &
    WTNDB    => plt_biom%WTNDB     , &
    WTNDBN   => plt_biom%WTNDBN    , &
    WTRSV    => plt_biom%WTRSV     , &
    WTGRB    => plt_biom%WTGRB     , &
    WTGRBN   => plt_biom%WTGRBN    , &
    WTSTKB   => plt_biom%WTSTKB    , &
    WTSHTN   => plt_biom%WTSHTN    , &
    WTSHTP   => plt_biom%WTSHTP    , &
    WTHSBP   => plt_biom%WTHSBP    , &
    WTSHTB   => plt_biom%WTSHTB    , &
    WTSTBN   => plt_biom%WTSTBN    , &
    WTHSKB   => plt_biom%WTHSKB    , &
    WTHSBN   => plt_biom%WTHSBN    , &
    WTEABP   => plt_biom%WTEABP    , &
    WTEARB   => plt_biom%WTEARB    , &
    WTEABN   => plt_biom%WTEABN    , &
    WTGRBP   => plt_biom%WTGRBP    , &
    WTSTXN   => plt_biom%WTSTXN    , &
    WGNODE   => plt_biom%WGNODE    , &
    WTRSBN   => plt_biom%WTRSBN    , &
    WTSTXP   => plt_biom%WTSTXP    , &
    WTNDBP   => plt_biom%WTNDBP    , &
    WTRSBP   => plt_biom%WTRSBP    , &
    WTLSB    => plt_biom%WTLSB     , &
    WTSTBP   => plt_biom%WTSTBP    , &
    WGNODN   => plt_biom%WGNODN    , &
    WGNODP   => plt_biom%WGNODP    , &
    WTRSVB   => plt_biom%WTRSVB    , &
    WTSTXB   => plt_biom%WTSTXB    , &
    CPOLNB   => plt_biom%CPOLNB    , &
    ZPOLNB   => plt_biom%ZPOLNB    , &
    PPOLNB   => plt_biom%PPOLNB    , &
    WGLFL    => plt_biom%WGLFL     , &
    WGLFLN   => plt_biom%WGLFLN    , &
    WGLFLP   => plt_biom%WGLFLP    , &
    WGLFP    => plt_biom%WGLFP     , &
    WTSHEB   => plt_biom%WTSHEB    , &
    CPOOL    => plt_biom%CPOOL     , &
    ZPOOL    => plt_biom%ZPOOL     , &
    PPOOL    => plt_biom%PPOOL     , &
    WGSHE    => plt_biom%WGSHE     , &
    WTLFBN   => plt_biom%WTLFBN    , &
    WGLFN    => plt_biom%WGLFN     , &
    WTLFBP   => plt_biom%WTLFBP    , &
    WSLF     => plt_biom%WSLF      , &
    WGLF     => plt_biom%WGLF      , &
    WTLFB    => plt_biom%WTLFB     , &
    WSSHE    => plt_biom%WSSHE     , &
    WTSHBN   => plt_biom%WTSHBN    , &
    WTSHBP   => plt_biom%WTSHBP    , &
    WGSHN    => plt_biom%WGSHN     , &
    WGSHP    => plt_biom%WGSHP     , &
    WTSTK    => plt_biom%WTSTK     , &
    CCPOLP   => plt_biom%CCPOLP    , &
    CCPLNP   => plt_biom%CCPLNP    , &
    WTLF     => plt_biom%WTLF      , &
    WTGR     => plt_biom%WTGR      , &
    WTSHT    => plt_biom%WTSHT     , &
    WTHSK    => plt_biom%WTHSK     , &
    WTEAR    => plt_biom%WTEAR     , &
    WTSHE    => plt_biom%WTSHE     , &
    WTSHTA   => plt_biom%WTSHTA    , &
    WTRT1    => plt_biom%WTRT1     , &
    WTRT1N   => plt_biom%WTRT1N    , &
    WTRT1P   => plt_biom%WTRT1P    , &
    RTWT1    => plt_biom%RTWT1     , &
    RTWT1N   => plt_biom%RTWT1N    , &
    RTWT1P   => plt_biom%RTWT1P    , &
    WTRT2    => plt_biom%WTRT2     , &
    WTRT2N   => plt_biom%WTRT2N    , &
    WTRT2P   => plt_biom%WTRT2P    , &
    WTNDLN   => plt_biom%WTNDLN    , &
    WTNDLP   => plt_biom%WTNDLP    , &
    WTNDL    => plt_biom%WTNDL     , &
    CPOOLN   => plt_biom%CPOOLN    , &
    ZPOOLN   => plt_biom%ZPOOLN    , &
    PPOOLN   => plt_biom%PPOOLN    , &
    ZEROP    => plt_biom%ZEROP     , &
    ZEROL    => plt_biom%ZEROL     , &
    WGLFV    => plt_biom%WGLFV     , &
    FWOODN   => plt_allom%FWOODN   , &
    FVRN     => plt_allom%FVRN     , &
    FWODR    => plt_allom%FWODR    , &
    FWODRP   => plt_allom%FWODRP   , &
    FWODRN   => plt_allom%FWODRN   , &
    FWOODP   => plt_allom%FWOODP   , &
    FWODSN   => plt_allom%FWODSN   , &
    FWODSP   => plt_allom%FWODSP   , &
    FWOOD    => plt_allom%FWOOD    , &
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
    CFOPC    =>  plt_soilchem%CFOPC, &
    CFOPN    =>  plt_soilchem%CFOPN, &
    CFOPP    =>  plt_soilchem%CFOPP, &
    H2GP     =>  plt_rbgc%H2GP     , &
    CO2P     =>  plt_rbgc%CO2P     , &
    CH4P     =>  plt_rbgc%CH4P     , &
    OXYP     =>  plt_rbgc%OXYP     , &
    H2GA     =>  plt_rbgc%H2GA     , &
    CH4A     =>  plt_rbgc%CH4A     , &
    OXYA     =>  plt_rbgc%OXYA     , &
    CO2A     =>  plt_rbgc%CO2A     , &
    CSNC     =>  plt_bgcr%CSNC     , &
    CNET     =>  plt_bgcr%CNET     , &
    ZSNC     =>  plt_bgcr%ZSNC     , &
    PSNC     =>  plt_bgcr%PSNC     , &
    TNBP     =>  plt_bgcr%TNBP     , &
    RCH4Z    =>  plt_bgcr%RCH4Z    , &
    RCO2Z    =>  plt_bgcr%RCO2Z    , &
    ROXYZ    =>  plt_bgcr%ROXYZ    , &
    RN2OZ    =>  plt_bgcr%RN2OZ    , &
    RNH3Z    =>  plt_bgcr%RNH3Z    , &
    RH2GZ    =>  plt_bgcr%RH2GZ    , &
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
!     PPI(NZ)=AMAX1(1.0,0.5*(PPI(NZ)+GRNO(NZ)/AREA3(NU)))
        PPX(NZ)=PPI(NZ)
        PP(NZ)=PPX(NZ)*AREA3(NU)
      ENDIF
      IF(IHVST(NZ).EQ.3)THEN
        CF(NZ)=CF(NZ)*HVST(NZ)
      ENDIF
      IF(IHVST(NZ).LE.2.AND.HVST(NZ).LT.0.0)THEN
        ARLFY=(1._r8-ABS(HVST(NZ)))*ARLFC
        ARLFR=0._r8
        DO 9875 L=1,JC1
          IF(ZL(L).GT.ZL(L-1).AND.ARLFT(L).GT.ZEROS &
            .AND.ARLFR.LT.ARLFY)THEN
            IF(ARLFR+ARLFT(L).GT.ARLFY)THEN
              HVST(NZ)=ZL(L-1)+((ARLFY-ARLFR) &
                /ARLFT(L))*(ZL(L)-ZL(L-1))
            ENDIF
          ELSE
            HVST(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+ARLFT(L)
9875    CONTINUE
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
        WHVSTT=HVST(NZ)*THIN(NZ)*0.45/24.0 &
          *AREA3(NU)*WTSHT(NZ)/WTSHTA(NZ)
      ELSE
        WHVSTT=0._r8
      ENDIF
      IF(IHVST(NZ).EQ.6)THEN
        WHVSTT=WHVSTT*TFN3(NZ)
      ENDIF
      CCPOLX=CCPOLP(NZ)/(1.0+CCPOLP(NZ))
      CCPLNX=CCPLNP(NZ)/(1.0+CCPLNP(NZ))
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
      WHVSLY=AMIN1(WTLF(NZ),WHVSLX)
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
      WTSHTT=WTSHE(NZ)+WTHSK(NZ)+WTEAR(NZ)+WTGR(NZ)
      IF(WTSHTT.GT.ZEROP(NZ))THEN
        WHVSHX=WHVSSX*WTSHE(NZ)/WTSHTT+WHVXXX
        WHVSHY=AMIN1(WTSHE(NZ),WHVSHX)
        WHVSHH=WHVSHY*(1._r8-CCPOLX)
        WHVSCS=WHVSHY*CCPOLX
        WHVSNS=WHVSHY*CCPLNX
        WHVXXX=AZMAX1(WHVSHX-WHVSHY)
        WHVHSX=WHVSSX*WTHSK(NZ)/WTSHTT+WHVXXX
        WHVHSY=AMIN1(WTHSK(NZ),WHVHSX)
        WHVHSH=WHVHSY
        WHVXXX=AZMAX1(WHVHSX-WHVHSY)
        WHVEAX=WHVSSX*WTEAR(NZ)/WTSHTT+WHVXXX
        WHVEAY=AMIN1(WTEAR(NZ),WHVEAX)
        WHVEAH=WHVEAY
        WHVXXX=AZMAX1(WHVEAX-WHVEAY)
        WHVGRX=WHVSSX*WTGR(NZ)/WTSHTT+WHVXXX
        WHVGRY=AMIN1(WTGR(NZ),WHVGRX)
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
      WTSTKT=WTSTK(NZ)+WTRSV(NZ)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
        WHVSTX=WHVSKX*WTSTK(NZ)/WTSTKT+WHVXXX
        WHVSTY=AMIN1(WTSTK(NZ),WHVSTX)
        WHVSTH=WHVSTY
        WHVXXX=AZMAX1(WHVSTX-WHVSTY)
        WHVRVX=WHVSKX*WTRSV(NZ)/WTSTKT+WHVXXX
        WHVRVY=AMIN1(WTRSV(NZ),WHVRVX)
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
          WHVSLY=AMIN1(WTLF(NZ)-WHVSLF-WHVSCL,WHVXXX)
          WHVSLF=WHVSLF+WHVSLY*(1._r8-CCPOLX)
          WHVSCL=WHVSCL+WHVSLY*CCPOLX
          WHVSNL=WHVSNL+WHVSLY*CCPLNX
          WHVXXX=AZMAX1(WHVXXX-WHVSLY)
          IF(WTSHTT.GT.ZEROP(NZ))THEN
            WHVSHX=WHVXXX*WTSHE(NZ)/WTSHTT
            WHVSHY=AMIN1(WTSHE(NZ),WHVSHX)
            WHVSHH=WHVSHH+WHVSHY*(1._r8-CCPOLX)
            WHVSCS=WHVSCS+WHVSHY*CCPOLX
            WHVSNS=WHVSNS+WHVSHY*CCPLNX
            WHVXXX=AZMAX1(WHVXXX-WHVSHY)
            WHVHSX=WHVXXX*WTHSK(NZ)/WTSHTT
            WHVHSY=AMIN1(WTHSK(NZ),WHVHSX)
            WHVHSH=WHVHSH+WHVHSY
            WHVXXX=AZMAX1(WHVXXX-WHVHSY)
            WHVEAX=WHVXXX*WTEAR(NZ)/WTSHTT
            WHVEAY=AMIN1(WTEAR(NZ),WHVEAX)
            WHVEAH=WHVEAH+WHVEAY
            WHVXXX=AZMAX1(WHVEAX-WHVEAY)
            WHVGRX=WHVXXX*WTGR(NZ)/WTSHTT
            WHVGRY=AMIN1(WTGR(NZ),WHVGRX)
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
      DO 9860 NB=1,NBR(NZ)
        DO  L=1,JC1
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=0._r8
          enddo
        enddo
9860  CONTINUE
      DO 9870 NB=1,NBR(NZ)
        DO  L=1,JC1
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=WGLFBL(L,NB,NZ)+WGLFL(L,K,NB,NZ)
          enddo
        enddo
9870  CONTINUE
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
    DO 9865 L=JC1,1,-1
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        IF(IHVST(NZ).NE.3)THEN
          IF(ZL(L).GT.ZL(L-1))THEN
            FHGT=AZMAX1(AMIN1(1.0,1._r8-((ZL(L)) &
              -HVST(NZ))/(ZL(L)-ZL(L-1))))
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
          .AND.WTLF(NZ).GT.ZEROL(NZ))THEN
          WHVSBL=WHVSLF*AZMAX1(WGLFBL(L,NB,NZ))/WTLF(NZ)
        ELSE
          WHVSBL=0._r8
        ENDIF
        DO 9845 K=JNODS1,0,-1
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBL.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(WGLFL(L,K,NB,NZ).GT.WHVSBL)THEN
                FHVST=AZMAX1(AMIN1(1.0,(WGLFL(L,K,NB,NZ)-WHVSBL)/WGLFL(L,K,NB,NZ)))
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
!     WTHTH1,WTHNH1,WTHPH1=harvested leaf C,N,P
!     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
!     WTHTH3,WTHNH3,WTHPH3=harvested woody C,N,P
!     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
            WHVSBL=WHVSBL-(1._r8-FHVST)*WGLFL(L,K,NB,NZ)
            WTHTH1=WTHTH1+(1._r8-FHVSH)*WGLFL(L,K,NB,NZ)*FWODB(1)
            WTHNH1=WTHNH1+(1._r8-FHVSH)*WGLFLN(L,K,NB,NZ)*FWODLN(1)
            WTHPH1=WTHPH1+(1._r8-FHVSH)*WGLFLP(L,K,NB,NZ)*FWODLP(1)
            WTHTX1=WTHTX1+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ)*FWODB(1)
            WTHNX1=WTHNX1+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ)*FWODLN(1)
            WTHPX1=WTHPX1+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ)*FWODLP(1)
            WTHTH3=WTHTH3+(1._r8-FHVSH)*WGLFL(L,K,NB,NZ)*FWODB(0)
            WTHNH3=WTHNH3+(1._r8-FHVSH)*WGLFLN(L,K,NB,NZ)*FWODLN(0)
            WTHPH3=WTHPH3+(1._r8-FHVSH)*WGLFLP(L,K,NB,NZ)*FWODLP(0)
            WTHTX3=WTHTX3+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ)*FWODB(0)
            WTHNX3=WTHNX3+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ)*FWODLN(0)
            WTHPX3=WTHPX3+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ)*FWODLP(0)
!
!     REMAINING LEAF C,N,P AND AREA
!
            WGLFL(L,K,NB,NZ)=FHVST*WGLFL(L,K,NB,NZ)
            WGLFLN(L,K,NB,NZ)=FHVST*WGLFLN(L,K,NB,NZ)
            WGLFLP(L,K,NB,NZ)=FHVST*WGLFLP(L,K,NB,NZ)
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
9865  CONTINUE
      DO 9835 NB=1,NBR(NZ)
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
        DO 9825 K=0,JNODS1
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
          DO 9815 L=1,JC1
            ARLFG=ARLFG+ARLFL(L,K,NB,NZ)
            WGLFG=WGLFG+WGLFL(L,K,NB,NZ)
            WGLFNG=WGLFNG+WGLFLN(L,K,NB,NZ)
            WGLFPG=WGLFPG+WGLFLP(L,K,NB,NZ)
            ARLFV(L,NZ)=ARLFV(L,NZ)+ARLFL(L,K,NB,NZ)
            WGLFV(L,NZ)=WGLFV(L,NZ)+WGLFL(L,K,NB,NZ)
9815      CONTINUE
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
            IF(WGLF(K,NB,NZ).GT.ZEROP(NZ) &
              .AND.EHVST(1,1,NZ).GT.0.0)THEN
              FHVSTK(K)=AZMAX1(AMIN1(1.0,(1._r8-(1._r8-AZMAX1(WGLFG) &
                /WGLF(K,NB,NZ))*EHVST(1,2,NZ)/EHVST(1,1,NZ))))
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
        WGLFGY=WGLFGY+WGLF(K,NB,NZ)
        WTLFB(NB,NZ)=WTLFB(NB,NZ)-WGLF(K,NB,NZ)+WGLFG
        WTLFBN(NB,NZ)=WTLFBN(NB,NZ)-WGLFN(K,NB,NZ)+WGLFNG
        WTLFBP(NB,NZ)=WTLFBP(NB,NZ)-WGLFP(K,NB,NZ)+WGLFPG
        ARLFB(NB,NZ)=ARLFB(NB,NZ)-ARLF1(K,NB,NZ)+ARLFG
        IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ))THEN
          WSLF(K,NB,NZ)=WSLF(K,NB,NZ)*ARLFG/ARLF1(K,NB,NZ)
        ELSE
          WSLF(K,NB,NZ)=0._r8
        ENDIF
        ARLF1(K,NB,NZ)=ARLFG
        WGLF(K,NB,NZ)=WGLFG
        WGLFN(K,NB,NZ)=WGLFNG
        WGLFP(K,NB,NZ)=WGLFPG
        WGLFGX=WGLFGX+WGLF(K,NB,NZ)
9825  CONTINUE
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
        .AND.WTSHE(NZ).GT.ZEROP(NZ))THEN
        WHVSBS=WHVSHH*WTSHEB(NB,NZ)/WTSHE(NZ)
      ELSE
        WHVSBS=0._r8
      ENDIF
      DO 9805 K=JNODS1,0,-1
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
!     WTHTH2,WTHNH2,WTHPH2=harvested petiole C,N,P
!     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     HTSHE,HTNODE=petiole,internode length
!
          IF((IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6).OR.WHVSBS.GT.0.0)THEN
            IF(IHVST(NZ).EQ.4.OR.IHVST(NZ).EQ.6)THEN
              IF(WGSHE(K,NB,NZ).GT.WHVSBS)THEN
                FHVSTK(K)=AZMAX1(AMIN1(1.0,(WGSHE(K,NB,NZ)-WHVSBS)/WGSHE(K,NB,NZ)))
                FHVSHK(K)=FHVSTK(K)
              ELSE
                FHVSTK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1._r8-FHVSTK(K))*WGSHE(K,NB,NZ)
            WTHTH2=WTHTH2+(1._r8-FHVSHK(K))*WGSHE(K,NB,NZ)*FWODB(1)
            WTHNH2=WTHNH2+(1._r8-FHVSHK(K))*WGSHN(K,NB,NZ)*FWODSN(1)
            WTHPH2=WTHPH2+(1._r8-FHVSHK(K))*WGSHP(K,NB,NZ)*FWODSP(1)
            WTHTX2=WTHTX2+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ)*FWODB(1)
            WTHNX2=WTHNX2+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ)*FWODSN(1)
            WTHPX2=WTHPX2+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ)*FWODSP(1)
            WTHTH3=WTHTH3+(1._r8-FHVSHK(K))*WGSHE(K,NB,NZ)*FWODB(0)
            WTHNH3=WTHNH3+(1._r8-FHVSHK(K))*WGSHN(K,NB,NZ)*FWODSN(0)
            WTHPH3=WTHPH3+(1._r8-FHVSHK(K))*WGSHP(K,NB,NZ)*FWODSP(0)
            WTHTX3=WTHTX3+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ)*FWODB(0)
            WTHNX3=WTHNX3+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ)*FWODSN(0)
            WTHPX3=WTHPX3+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ)*FWODSP(0)
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     WGSHE=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     HTSHE=node petiole height
!     WSSHE=petiole protein mass
!
            WGSHGY=WGSHGY+WGSHE(K,NB,NZ)
            WTSHEB(NB,NZ)=WTSHEB(NB,NZ) &
              -(1._r8-FHVSTK(K))*WGSHE(K,NB,NZ)
            WTSHBN(NB,NZ)=WTSHBN(NB,NZ) &
              -(1._r8-FHVSTK(K))*WGSHN(K,NB,NZ)
            WTSHBP(NB,NZ)=WTSHBP(NB,NZ) &
              -(1._r8-FHVSTK(K))*WGSHP(K,NB,NZ)
            WGSHE(K,NB,NZ)=FHVSTK(K)*WGSHE(K,NB,NZ)
            WSSHE(K,NB,NZ)=FHVSTK(K)*WSSHE(K,NB,NZ)
            WGSHN(K,NB,NZ)=FHVSTK(K)*WGSHN(K,NB,NZ)
            WGSHP(K,NB,NZ)=FHVSTK(K)*WGSHP(K,NB,NZ)
            WSSHE(K,NB,NZ)=FHVSTK(K)*WSSHE(K,NB,NZ)
            IF(IHVST(NZ).LE.2 &
              .AND.HTSHE(K,NB,NZ).GT.0.0)THEN
              FHGT=AZMAX1(AMIN1(1.0,(HTNODE(K,NB,NZ) &
                +HTSHE(K,NB,NZ)-HVST(NZ))/HTSHE(K,NB,NZ)))
              HTSHE(K,NB,NZ)=(1._r8-FHGT)*HTSHE(K,NB,NZ)
            ELSE
              HTSHE(K,NB,NZ)=FHVSTK(K)*HTSHE(K,NB,NZ)
            ENDIF
            WGSHGX=WGSHGX+WGSHE(K,NB,NZ)
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
9805    CONTINUE
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
        CPOOLX=AZMAX1(CPOOL(NB,NZ))
        ZPOOLX=AZMAX1(ZPOOL(NB,NZ))
        PPOOLX=AZMAX1(PPOOL(NB,NZ))
        CPOLNX=AZMAX1(CPOLNB(NB,NZ))
        ZPOLNX=AZMAX1(ZPOLNB(NB,NZ))
        PPOLNX=AZMAX1(PPOLNB(NB,NZ))
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          IF(WGLFGY+WGSHGY.GT.ZEROP(NZ))THEN
            FHVST=AZMAX1(AMIN1(1.0,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
            CPOOLG=CPOOLX*FHVST
            ZPOOLG=ZPOOLX*FHVST
            PPOOLG=PPOOLX*FHVST
            CPOLNG=CPOLNX*FHVST
            ZPOLNG=ZPOLNX*FHVST
            PPOLNG=PPOLNX*FHVST
            WTNDG=WTNDB(NB,NZ)*FHVST
            WTNDNG=WTNDBN(NB,NZ)*FHVST
            WTNDPG=WTNDBP(NB,NZ)*FHVST
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
            IF(CPOOL(NB,NZ).GT.ZEROP(NZ))THEN
              WHVSCX=AZMAX1(WHVSCP)*WTLSBX/WTLS(NZ)
              CPOOLG=AZMAX1(CPOOLX-WHVSCX)
              ZPOOLG=AZMAX1(ZPOOLX-WHVSCX*ZPOOLX/CPOOL(NB,NZ))
              PPOOLG=AZMAX1(PPOOLX-WHVSCX*PPOOLX/CPOOL(NB,NZ))
            ELSE
              CPOOLG=0._r8
              ZPOOLG=0._r8
              PPOOLG=0._r8
            ENDIF
            IF(CPOLNB(NB,NZ).GT.ZEROP(NZ))THEN
              WHVSNX=AZMAX1(WHVSNP)*WTLSBX/WTLS(NZ)
              CPOLNG=AZMAX1(CPOLNX-WHVSNX)
              ZPOLNG=AZMAX1(ZPOLNX-WHVSNX*ZPOLNX/CPOLNB(NB,NZ))
              PPOLNG=AZMAX1(PPOLNX-WHVSNX*PPOLNX/CPOLNB(NB,NZ))
              WTNDG=WTNDB(NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDNG=WTNDBN(NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDPG=WTNDBP(NB,NZ)*(1._r8-WHVSNX/CPOLNX)
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
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!
        WTHTH0=WTHTH0+CPOOLX-CPOOLG+CPOLNX-CPOLNG
        WTHNH0=WTHNH0+ZPOOLX-ZPOOLG+ZPOLNX-ZPOLNG
        WTHPH0=WTHPH0+PPOOLX-PPOOLG+PPOLNX-PPOLNG
        WTHTH0=WTHTH0+WTNDB(NB,NZ)-WTNDG
        WTHNH0=WTHNH0+WTNDBN(NB,NZ)-WTNDNG
        WTHPH0=WTHPH0+WTNDBP(NB,NZ)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
        CPOOL(NB,NZ)=CPOOLG
        ZPOOL(NB,NZ)=ZPOOLG
        PPOOL(NB,NZ)=PPOOLG
        CPOLNB(NB,NZ)=CPOLNG
        ZPOLNB(NB,NZ)=ZPOLNG
        PPOLNB(NB,NZ)=PPOLNG
        WTNDB(NB,NZ)=WTNDG
        WTNDBN(NB,NZ)=WTNDNG
        WTNDBP(NB,NZ)=WTNDPG
!
!     REMOVE C4 NON-STRUCTURAL C
!
!     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!     FHVST4=fraction of nonstructural mass not harvested
!     CPOOLG=branch non-structural C mass after harvest
!     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!
        IF(ICTYP(NZ).EQ.4.AND.CPOOLX.GT.ZEROP(NZ))THEN
          FHVST4=CPOOLG/CPOOLX
          DO 9810 K=1,JNODS1
            WTHTH0=WTHTH0+(1._r8-FHVST4)*CPOOL3(K,NB,NZ)
            WTHTH0=WTHTH0+(1._r8-FHVST4)*CPOOL4(K,NB,NZ)
            WTHTH0=WTHTH0+(1._r8-FHVST4)*CO2B(K,NB,NZ)
            WTHTH0=WTHTH0+(1._r8-FHVST4)*HCOB(K,NB,NZ)
            CPOOL3(K,NB,NZ)=FHVST4*CPOOL3(K,NB,NZ)
            CPOOL4(K,NB,NZ)=FHVST4*CPOOL4(K,NB,NZ)
            CO2B(K,NB,NZ)=FHVST4*CO2B(K,NB,NZ)
            HCOB(K,NB,NZ)=FHVST4*HCOB(K,NB,NZ)
9810      CONTINUE
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
              FHGT=AZMAX1(AMIN1(1.0,HVST(NZ)/HTSTKX))
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
          IF(WTSTK(NZ).GT.ZEROL(NZ))THEN
            FHVST=AZMAX1(AMIN1(1.0,1._r8-WHVSTH/WTSTK(NZ)))
            FHVSH=FHVST
          ELSE
            FHVST=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ENDIF
    !
!     HARVESTED STALK C,N,P
!
!     WTHTH3,WTHNH3,WTHPH3=harvested stalk C,N,P
!     WTHTX3,WTHNX3,WTHPX3=harvested stalk C,N,P to litter
!     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in harvested stalk
!
        WTHTH3=WTHTH3+(1._r8-FHVSH)*WTSTKB(NB,NZ)
        WTHNH3=WTHNH3+(1._r8-FHVSH)*WTSTBN(NB,NZ)
        WTHPH3=WTHPH3+(1._r8-FHVSH)*WTSTBP(NB,NZ)
        WTHTX3=WTHTX3+(FHVSH-FHVST)*WTSTKB(NB,NZ)
        WTHNX3=WTHNX3+(FHVSH-FHVST)*WTSTBN(NB,NZ)
        WTHPX3=WTHPX3+(FHVSH-FHVST)*WTSTBP(NB,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
        WTSTKB(NB,NZ)=FHVST*WTSTKB(NB,NZ)
        WTSTBN(NB,NZ)=FHVST*WTSTBN(NB,NZ)
        WTSTBP(NB,NZ)=FHVST*WTSTBP(NB,NZ)
        WVSTKB(NB,NZ)=FHVST*WVSTKB(NB,NZ)
        WTSTXB(NB,NZ)=FHVST*WTSTXB(NB,NZ)
        WTSTXN(NB,NZ)=FHVST*WTSTXN(NB,NZ)
        WTSTXP(NB,NZ)=FHVST*WTSTXP(NB,NZ)
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
        DO 9820 K=JNODS1,0,-1
          IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
            IF(HTNODX(K,NB,NZ).GT.ZERO)THEN
              IF(IHVST(NZ).NE.3)THEN
                FHGTK=AZMAX1(AMIN1(1.0,(HTNODE(K,NB,NZ) &
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
            IF(WTSTK(NZ).GT.ZEROP(NZ))THEN
              FHVSTS=AZMAX1(AMIN1(1.0,1._r8-WHVSTH/WTSTK(NZ)))
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ENDIF
          WGNODE(K,NB,NZ)=FHVSTS*WGNODE(K,NB,NZ)
          WGNODN(K,NB,NZ)=FHVSTS*WGNODN(K,NB,NZ)
          WGNODP(K,NB,NZ)=FHVSTS*WGNODP(K,NB,NZ)
          IF(IHVST(NZ).LE.2.AND.test_aeqb(THIN(NZ),0._r8))THEN
            HTNODX(K,NB,NZ)=FHVSTS*HTNODX(K,NB,NZ)
            HTNODE(K,NB,NZ)=AMIN1(HTNODE(K,NB,NZ),HVST(NZ))
          ENDIF

9820    CONTINUE
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
          IF(WTSTKB(NB,NZ).GT.ZEROP(NZ))THEN
            FHVST=FHVST
            FHVSH=FHVSH
          ELSE
            FHVST=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(WTRSV(NZ).GT.ZEROP(NZ))THEN
            FHVST=AZMAX1(AMIN1(1.0,1._r8-WHVRVH/WTRSV(NZ)))
            FHVSH=FHVST
          ELSE
            FHVST=0._r8
            FHVSH=0._r8
          ENDIF
        ENDIF
!
!     HARVESTED STALK RESERVE C,N,P
!
!     WTHTH3,WTHNH3,WTHPH3=harvested stalk C,N,P
!     WTHTX3,WTHNX3,WTHPX3=harvested stalk C,N,P to litter
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!
        WTHTH3=WTHTH3+(1._r8-FHVSH)*WTRSVB(NB,NZ)
        WTHNH3=WTHNH3+(1._r8-FHVSH)*WTRSBN(NB,NZ)
        WTHPH3=WTHPH3+(1._r8-FHVSH)*WTRSBP(NB,NZ)
        WTHTX3=WTHTX3+(FHVSH-FHVST)*WTRSVB(NB,NZ)
        WTHNX3=WTHNX3+(FHVSH-FHVST)*WTRSBN(NB,NZ)
        WTHPX3=WTHPX3+(FHVSH-FHVST)*WTRSBP(NB,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
        WTRSVB(NB,NZ)=FHVST*WTRSVB(NB,NZ)
        WTRSBN(NB,NZ)=FHVST*WTRSBN(NB,NZ)
        WTRSBP(NB,NZ)=FHVST*WTRSBP(NB,NZ)
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
          IF(HVST(NZ).LT.HTSTKX &
            .OR.IHVST(NZ).EQ.1 &
            .OR.IHVST(NZ).EQ.3)THEN
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
          IF(WTHSK(NZ).GT.ZEROP(NZ))THEN
            FHVSTH=AZMAX1(AMIN1(1.0,1._r8-WHVHSH/WTHSK(NZ)))
            FHVSHH=FHVSTH
          ELSE
            FHVSTH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(WTEAR(NZ).GT.ZEROP(NZ))THEN
            FHVSTE=AZMAX1(AMIN1(1.0,1._r8-WHVEAH/WTEAR(NZ)))
            FHVSHE=FHVSTE
          ELSE
            FHVSTE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(WTGR(NZ).GT.ZEROP(NZ))THEN
            FHVSTG=AZMAX1(AMIN1(1.0,1._r8-WHVGRH/WTGR(NZ)))
            FHVSHG=FHVSTG
          ELSE
            FHVSTG=1.0_r8
            FHVSHG=1.0_r8
          ENDIF
        ENDIF
!
!     HARVESTED REPRODUCTIVE C,N,P
!
!     WTHTH2,WTHNH2,WTHPH2=reproductive C,N,P removed
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTHTG,WTHNG,WTHPG=grain harvested
!
        WTHTH2=WTHTH2+(1._r8-FHVSHH)*WTHSKB(NB,NZ)+(1._r8-FHVSHE) &
          *WTEARB(NB,NZ)+(1._r8-FHVSHG)*WTGRB(NB,NZ)
        WTHNH2=WTHNH2+(1._r8-FHVSHH)*WTHSBN(NB,NZ)+(1._r8-FHVSHE) &
          *WTEABN(NB,NZ)+(1._r8-FHVSHG)*WTGRBN(NB,NZ)
        WTHPH2=WTHPH2+(1._r8-FHVSHH)*WTHSBP(NB,NZ)+(1._r8-FHVSHE) &
          *WTEABP(NB,NZ)+(1._r8-FHVSHG)*WTGRBP(NB,NZ)
        WTHTX2=WTHTX2+(FHVSHH-FHVSTH)*WTHSKB(NB,NZ)+(FHVSHE-FHVSTE) &
          *WTEARB(NB,NZ)+(FHVSHG-FHVSTG)*WTGRB(NB,NZ)
        WTHNX2=WTHNX2+(FHVSHH-FHVSTH)*WTHSBN(NB,NZ)+(FHVSHE-FHVSTE) &
          *WTEABN(NB,NZ)+(FHVSHG-FHVSTG)*WTGRBN(NB,NZ)
        WTHPX2=WTHPX2+(FHVSHH-FHVSTH)*WTHSBP(NB,NZ)+(FHVSHE-FHVSTE) &
          *WTEABP(NB,NZ)+(FHVSHG-FHVSTG)*WTGRBP(NB,NZ)
        WTHTG=WTHTG+(1._r8-FHVSTG)*WTGRB(NB,NZ)
        WTHNG=WTHNG+(1._r8-FHVSTG)*WTGRBN(NB,NZ)
        WTHPG=WTHPG+(1._r8-FHVSTG)*WTGRBP(NB,NZ)
!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
        WTHSKB(NB,NZ)=FHVSTH*WTHSKB(NB,NZ)
        WTEARB(NB,NZ)=FHVSTE*WTEARB(NB,NZ)
        WTGRB(NB,NZ)=FHVSTG*WTGRB(NB,NZ)
        WTHSBN(NB,NZ)=FHVSTH*WTHSBN(NB,NZ)
        WTEABN(NB,NZ)=FHVSTE*WTEABN(NB,NZ)
        WTGRBN(NB,NZ)=FHVSTG*WTGRBN(NB,NZ)
        WTHSBP(NB,NZ)=FHVSTH*WTHSBP(NB,NZ)
        WTEABP(NB,NZ)=FHVSTE*WTEABP(NB,NZ)
        WTGRBP(NB,NZ)=FHVSTG*WTGRBP(NB,NZ)
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
        DO 1325 K=1,JNODS1
          CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
            +CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
            +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
1325    CONTINUE
        WTLSB(NB,NZ)=AZMAX1(WTLFB(NB,NZ) &
          +WTSHEB(NB,NZ))
        WTSHTB(NB,NZ)=AZMAX1(WTLFB(NB,NZ) &
          +WTSHEB(NB,NZ)+WTSTKB(NB,NZ)+WTRSVB(NB,NZ) &
          +WTHSKB(NB,NZ)+WTEARB(NB,NZ)+WTGRB(NB,NZ) &
          +CPOOL(NB,NZ)+CPOOLK(NB,NZ))
        WTSHTN(NB,NZ)=AZMAX1(WTLFBN(NB,NZ) &
          +WTSHBN(NB,NZ)+WTSTBN(NB,NZ)+WTRSBN(NB,NZ) &
          +WTHSBN(NB,NZ)+WTEABN(NB,NZ)+WTGRBN(NB,NZ) &
          +ZPOOL(NB,NZ))
        WTSHTP(NB,NZ)=AZMAX1(WTLFBP(NB,NZ) &
          +WTSHBP(NB,NZ)+WTSTBP(NB,NZ)+WTRSBP(NB,NZ) &
          +WTHSBP(NB,NZ)+WTEABP(NB,NZ)+WTGRBP(NB,NZ) &
          +PPOOL(NB,NZ))
        VOLWPX=VOLWP(NZ)
        WVPLT=AZMAX1(WTLS(NZ)+WVSTK(NZ))
        APSILT=ABS(PSILT(NZ))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWP(NZ)=1.0E-06_r8*WVPLT/FDM
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
            DO 3005 M=2,10
              IDAY(M,NB,NZ)=0
3005        CONTINUE
            IFLGA(NB,NZ)=0
            IF(NB.EQ.NB1(NZ))THEN
              DO 3010 NBX=1,NBR(NZ)
                IF(NBX.NE.NB1(NZ))THEN
                  GROUP(NBX,NZ)=GROUPI(NZ)
                  PSTGI(NBX,NZ)=PSTG(NBX,NZ)
                  PSTGF(NBX,NZ)=0._r8
                  VSTGX(NBX,NZ)=0._r8
                  TGSTGI(NBX,NZ)=0._r8
                  TGSTGF(NBX,NZ)=0._r8
                  FLG4(NBX,NZ)=0._r8
                  IDAY(1,NBX,NZ)=I
                  DO 3015 M=2,10
                    IDAY(M,NBX,NZ)=0
3015              CONTINUE
                  IFLGA(NBX,NZ)=0
                ENDIF
3010          CONTINUE
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
9835  CONTINUE
      WTLS(NZ)=0._r8
      WTSTK(NZ)=0._r8
      WVSTK(NZ)=0._r8
      ARSTP(NZ)=0._r8
      DO 9840 NB=1,NBR(NZ)
        WTLS(NZ)=WTLS(NZ)+WTLSB(NB,NZ)
        WTSTK(NZ)=WTSTK(NZ)+WTSTKB(NB,NZ)
        WVSTK(NZ)=WVSTK(NZ)+WVSTKB(NB,NZ)
        DO 9830 L=1,JC1
          ARSTP(NZ)=ARSTP(NZ)+ARSTK(L,NB,NZ)
9830    CONTINUE
9840  CONTINUE
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
        XHVST=1.0_r8-THIN(NZ)
        DO 3985 N=1,MY(NZ)
          DO 3980 L=NU,NJ
            IF(IHVST(NZ).NE.5)THEN
              XHVST=1.0_r8-THIN(NZ)
              XHVSN=XHVST
              XHVSP=XHVST
              FFIRE=0._r8
              FFIRN=0._r8
              FFIRP=0._r8
            ELSE
              IF(THETW(L).GT.FVLWB.OR.CORGC(L).LE.FORGC &
                .OR.ITILL.NE.22)THEN
                XHVST=1.0_r8
                XHVSN=XHVST
                XHVSP=XHVST
                FFIRE=0._r8
                FFIRN=0._r8
                FFIRP=0._r8
              ELSE
                XHVST=1.0_r8-DCORP*EHVST(1,3,NZ) &
                  *AMIN1(1.0,(CORGC(L)-FORGC)/(0.55E+06-FORGC))
                XHVSN=XHVST
                XHVSP=XHVST
                FFIRE=EHVST(2,3,NZ)
                FFIRN=FFIRE*EFIRE(1,IHVST(NZ))
                FFIRP=FFIRE*EFIRE(2,IHVST(NZ))
              ENDIF
            ENDIF
            DO 3385 M=1,jsken
              FHVST=(1._r8-XHVST)*CFOPC(0,M,NZ)*CPOOLR(N,L,NZ)
              FHVSN=(1._r8-XHVSN)*CFOPN(0,M,NZ)*ZPOOLR(N,L,NZ)
              FHVSP=(1._r8-XHVSP)*CFOPP(0,M,NZ)*PPOOLR(N,L,NZ)
              CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-FFIRE)*FHVST
              ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-FFIRN)*FHVSN
              PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-FFIRP)*FHVSP
              VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
              VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE*FHVST
              VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
              VNH3F(NZ)=VNH3F(NZ)-FFIRN*FHVSN
              VN2OF(NZ)=VN2OF(NZ)-0.0
              VPO4F(NZ)=VPO4F(NZ)-FFIRP*FHVSP
              CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
              TNBP=TNBP-FCH4F*FFIRE*FHVST
              DO NR=1,NRT(NZ)
                FHVST=(1._r8-XHVST)*CFOPC(5,M,NZ)*(WTRT1(N,L,NR,NZ) &
                  +WTRT2(N,L,NR,NZ))*FWODR(0)
                FHVSN=(1._r8-XHVSN)*CFOPN(5,M,NZ)*(WTRT1N(N,L,NR,NZ) &
                  +WTRT2N(N,L,NR,NZ))*FWODRN(0)
                FHVSP=(1._r8-XHVSP)*CFOPP(5,M,NZ)*(WTRT1P(N,L,NR,NZ) &
                  +WTRT2P(N,L,NR,NZ))*FWODRP(0)
                CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-FFIRE)*FHVST
                ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-FFIRN)*FHVSN
                PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-FFIRP)*FHVSP
                VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE*FHVST
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
                VNH3F(NZ)=VNH3F(NZ)-FFIRN*FHVSN
                VN2OF(NZ)=VN2OF(NZ)-0.0
                VPO4F(NZ)=VPO4F(NZ)-FFIRP*FHVSP
                CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                TNBP=TNBP-FCH4F*FFIRE*FHVST
                FHVST=(1._r8-XHVST)*CFOPC(4,M,NZ)*(WTRT1(N,L,NR,NZ) &
                  +WTRT2(N,L,NR,NZ))*FWODR(1)
                FHVSN=(1._r8-XHVSN)*CFOPN(4,M,NZ)*(WTRT1N(N,L,NR,NZ) &
                  +WTRT2N(N,L,NR,NZ))*FWODRN(1)
                FHVSP=(1._r8-XHVSP)*CFOPP(4,M,NZ)*(WTRT1P(N,L,NR,NZ) &
                  +WTRT2P(N,L,NR,NZ))*FWODRP(1)
                CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-FFIRE)*FHVST
                ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-FFIRN)*FHVSN
                PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-FFIRP)*FHVSP
                VCO2F(NZ)=VCO2F(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                VCH4F(NZ)=VCH4F(NZ)-FCH4F*FFIRE*FHVST
                VOXYF(NZ)=VOXYF(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
                VNH3F(NZ)=VNH3F(NZ)-FFIRN*FHVSN
                VN2OF(NZ)=VN2OF(NZ)-0.0
                VPO4F(NZ)=VPO4F(NZ)-FFIRP*FHVSP
                CNET(NZ)=CNET(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                TNBP=TNBP-FCH4F*FFIRE*FHVST
              enddo
3385        CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING HARVESTING
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Z(NZ)=RCO2Z(NZ)-(1._r8-XHVST) &
              *(CO2A(N,L,NZ)+CO2P(N,L,NZ))
            ROXYZ(NZ)=ROXYZ(NZ)-(1._r8-XHVST) &
              *(OXYA(N,L,NZ)+OXYP(N,L,NZ))
            RCH4Z(NZ)=RCH4Z(NZ)-(1._r8-XHVST) &
              *(CH4A(N,L,NZ)+CH4P(N,L,NZ))
            RN2OZ(NZ)=RN2OZ(NZ)-(1._r8-XHVST) &
              *(Z2OA(N,L,NZ)+Z2OP(N,L,NZ))
            RNH3Z(NZ)=RNH3Z(NZ)-(1._r8-XHVST) &
              *(ZH3A(N,L,NZ)+ZH3P(N,L,NZ))
            RH2GZ(NZ)=RH2GZ(NZ)-(1._r8-XHVST) &
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
            DO 3960 NR=1,NRT(NZ)
              WTRT1(N,L,NR,NZ)=WTRT1(N,L,NR,NZ)*XHVST
              WTRT2(N,L,NR,NZ)=WTRT2(N,L,NR,NZ)*XHVST
              WTRT1N(N,L,NR,NZ)=WTRT1N(N,L,NR,NZ)*XHVSN
              WTRT2N(N,L,NR,NZ)=WTRT2N(N,L,NR,NZ)*XHVSN
              WTRT1P(N,L,NR,NZ)=WTRT1P(N,L,NR,NZ)*XHVSP
              WTRT2P(N,L,NR,NZ)=WTRT2P(N,L,NR,NZ)*XHVSP
              RTWT1(N,NR,NZ)=RTWT1(N,NR,NZ)*XHVST
              RTWT1N(N,NR,NZ)=RTWT1N(N,NR,NZ)*XHVST
              RTWT1P(N,NR,NZ)=RTWT1P(N,NR,NZ)*XHVST
              RTLG1(N,L,NR,NZ)=RTLG1(N,L,NR,NZ)*XHVST
              RTLG2(N,L,NR,NZ)=RTLG2(N,L,NR,NZ)*XHVST
              RTN2(N,L,NR,NZ)=RTN2(N,L,NR,NZ)*XHVST
3960        CONTINUE
            CPOOLR(N,L,NZ)=CPOOLR(N,L,NZ)*XHVST
            ZPOOLR(N,L,NZ)=ZPOOLR(N,L,NZ)*XHVSN
            PPOOLR(N,L,NZ)=PPOOLR(N,L,NZ)*XHVSP
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
!     NODULE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(INTYP(NZ).NE.0.AND.N.EQ.1)THEN
              DO 3395 M=1,jsken
                CSNC(M,1,L,NZ)=CSNC(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPC(4,M,NZ)*WTNDL(L,NZ) &
                  +CFOPC(0,M,NZ)*CPOOLN(L,NZ))
                ZSNC(M,1,L,NZ)=ZSNC(M,1,L,NZ)+(1._r8-XHVSN) &
                  *(CFOPN(4,M,NZ)*WTNDLN(L,NZ) &
                  +CFOPN(0,M,NZ)*ZPOOLN(L,NZ))
                PSNC(M,1,L,NZ)=PSNC(M,1,L,NZ)+(1._r8-XHVSP) &
                  *(CFOPP(4,M,NZ)*WTNDLP(L,NZ) &
                  +CFOPP(0,M,NZ)*PPOOLN(L,NZ))
3395          CONTINUE
              WTNDL(L,NZ)=WTNDL(L,NZ)*XHVST
              WTNDLN(L,NZ)=WTNDLN(L,NZ)*XHVSN
              WTNDLP(L,NZ)=WTNDLP(L,NZ)*XHVSP
              CPOOLN(L,NZ)=CPOOLN(L,NZ)*XHVST
              ZPOOLN(L,NZ)=ZPOOLN(L,NZ)*XHVSN
              PPOOLN(L,NZ)=PPOOLN(L,NZ)*XHVSP
            ENDIF
3980      CONTINUE
3985    CONTINUE
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
          DO 3400 M=1,jsken
            CSNC(M,0,NG(NZ),NZ)=CSNC(M,0,NG(NZ),NZ) &
              +((1._r8-XHVST)*CFOPC(0,M,NZ)*WTRVC(NZ))*FWOOD(0)
            ZSNC(M,0,NG(NZ),NZ)=ZSNC(M,0,NG(NZ),NZ) &
              +((1._r8-XHVSN)*CFOPN(0,M,NZ)*WTRVN(NZ))*FWOODN(0)
            PSNC(M,0,NG(NZ),NZ)=PSNC(M,0,NG(NZ),NZ) &
              +((1._r8-XHVSP)*CFOPP(0,M,NZ)*WTRVP(NZ))*FWOODP(0)
            CSNC(M,1,NG(NZ),NZ)=CSNC(M,1,NG(NZ),NZ) &
              +((1._r8-XHVST)*CFOPC(0,M,NZ)*WTRVC(NZ))*FWOOD(1)
            ZSNC(M,1,NG(NZ),NZ)=ZSNC(M,1,NG(NZ),NZ) &
              +((1._r8-XHVSN)*CFOPN(0,M,NZ)*WTRVN(NZ))*FWOODN(1)
            PSNC(M,1,NG(NZ),NZ)=PSNC(M,1,NG(NZ),NZ) &
              +((1._r8-XHVSP)*CFOPP(0,M,NZ)*WTRVP(NZ))*FWOODP(1)
3400      CONTINUE
          WTRVC(NZ)=WTRVC(NZ)*XHVST
          WTRVN(NZ)=WTRVN(NZ)*XHVSN
          WTRVP(NZ)=WTRVP(NZ)*XHVSP
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

end module PlantDisturbsMod
