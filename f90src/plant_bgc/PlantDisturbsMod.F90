module PlantDisturbsMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  implicit none
  private

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
    EHVSTs1    =>  plt_distb%EHVSTs1  , &
    HVSTs1     =>  plt_distb%HVSTs1   , &
    IHVSTs1    =>  plt_distb%IHVSTs1  , &
    THINs1     =>  plt_distb%THINs1   , &
    PPs1       =>  plt_site%PPs1      , &
    NUs1       =>  plt_site%NUs1      , &
    ZEROs1     =>  plt_site%ZEROs1    , &
    ZNOONs1    =>  plt_site%ZNOONs1   , &
    AREA3s1    =>  plt_site%AREA3s1   , &
    ZEROQs1    =>  plt_rbgc%ZEROQs1   , &
    ZEROPs1    =>  plt_biom%ZEROPs1   , &
    ZEROLs1    =>  plt_biom%ZEROLs1   , &
    WTSTDGs1   =>  plt_biom%WTSTDGs1  , &
    WTSTDNs1   =>  plt_biom%WTSTDNs1  , &
    WTSTDPs1   =>  plt_biom%WTSTDPs1  , &
    WTSTGs1    =>  plt_biom%WTSTGs1     &
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
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
!          IHVST=3:reduction of clumping factor
!          IHVST=4 or 6:animal or insect biomass(g LM m-2),IHVST=5:fire
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     WTHTH4,WTHNH4,WTHPH4=harvested standing dead C,N,P
!     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
!
  IF(IHVSTs1(NZ).GE.0)THEN
    IF(J.EQ.INT(ZNOONs1).AND.IHVSTs1(NZ).NE.4 &
      .AND.IHVSTs1(NZ).NE.6)THEN
      IF(test_aeqb(THINs1(NZ),0._r8))THEN
        FHVST=AMAX1(0.0,1._r8-EHVSTs1(1,4,NZ))
        FHVSH=FHVST
      ELSE
        FHVST=AMAX1(0.0,1._r8-THINs1(NZ))
        IF(IHVSTs1(NZ).EQ.0)THEN
          FHVSH=AMAX1(0.0,1._r8-EHVSTs1(1,4,NZ)*THINs1(NZ))
        ELSE
          FHVSH=FHVST
        ENDIF
      ENDIF
    ELSEIF(IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6)THEN
      IF(WTSTGs1(NZ).GT.ZEROPs1(NZ))THEN
        WHVSTD=HVSTs1(NZ)*THINs1(NZ)*0.45/24.0 &
          *AREA3s1(NUs1)*EHVSTs1(1,4,NZ)
        FHVST=AMAX1(0.0,1._r8-WHVSTD/WTSTGs1(NZ))
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
      WTHTH4=WTHTH4+(1._r8-FHVSH)*WTSTDGs1(M,NZ)
      WTHNH4=WTHNH4+(1._r8-FHVSH)*WTSTDNs1(M,NZ)
      WTHPH4=WTHPH4+(1._r8-FHVSH)*WTSTDPs1(M,NZ)
      WTHTX4=WTHTX4+(FHVSH-FHVST)*WTSTDGs1(M,NZ)
      WTHNX4=WTHNX4+(FHVSH-FHVST)*WTSTDNs1(M,NZ)
      WTHPX4=WTHPX4+(FHVSH-FHVST)*WTSTDPs1(M,NZ)
      WTSTDGs1(M,NZ)=FHVST*WTSTDGs1(M,NZ)
      WTSTDNs1(M,NZ)=FHVST*WTSTDNs1(M,NZ)
      WTSTDPs1(M,NZ)=FHVST*WTSTDPs1(M,NZ)
6475  CONTINUE
!
    call PlantDisturbance(I,J,NZ)

    ZEROPs1(NZ)=ZEROs1*PPs1(NZ)
    ZEROQs1(NZ)=ZEROs1*PPs1(NZ)/AREA3s1(NUs1)
    ZEROLs1(NZ)=ZEROs1*PPs1(NZ)*1.0E+06
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
    WTSTDGs1   =>  plt_biom%WTSTDGs1  , &
    WTSTDNs1   =>  plt_biom%WTSTDNs1  , &
    WTSTDPs1   =>  plt_biom%WTSTDPs1  , &
    IHVSTs1    =>  plt_distb%IHVSTs1  , &
    FWOODs1    =>  plt_allom%FWOODs1  , &
    FWOODPs1   =>  plt_allom%FWOODPs1 , &
    FWOODNs1   =>  plt_allom%FWOODNs1 , &
    CFOPCs1    =>  plt_soilchem%CFOPCs1  , &
    CFOPNs1    =>  plt_soilchem%CFOPNs1  , &
    CFOPPs1    =>  plt_soilchem%CFOPPs1  , &
    CSNCs1     =>  plt_bgcr%CSNCs1    , &
    ZSNCs1     =>  plt_bgcr%ZSNCs1    , &
    PSNCs1     =>  plt_bgcr%PSNCs1    , &
    TPSNCs1    =>  plt_bgcr%TPSNCs1   , &
    TCSNCs1    =>  plt_bgcr%TCSNCs1   , &
    TPSN0s1    =>  plt_bgcr%TPSN0s1   , &
    TZSN0s1    =>  plt_bgcr%TZSN0s1   , &
    TCSN0s1    =>  plt_bgcr%TCSN0s1   , &
    TZSNCs1    =>  plt_bgcr%TZSNCs1   , &
    IBTYPs1    =>  plt_pheno%IBTYPs1  , &
    IGTYPs1    =>  plt_pheno%IGTYPs1    &
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
  IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
    IF(IHVSTs1(NZ).NE.5)THEN
      DO 6375 M=1,4
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
          +CFOPCs1(0,M,NZ)*(WTHTR0+WTHTX0) &
          +CFOPCs1(1,M,NZ)*(WTHTR1+WTHTX1) &
          +CFOPCs1(2,M,NZ)*(WTHTR2+WTHTX2)
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
          +CFOPNs1(0,M,NZ)*(WTHNR0+WTHNX0) &
          +CFOPNs1(1,M,NZ)*(WTHNR1+WTHNX1) &
          +CFOPNs1(2,M,NZ)*(WTHNR2+WTHNX2)
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
          +CFOPPs1(0,M,NZ)*(WTHPR0+WTHPX0) &
          +CFOPPs1(1,M,NZ)*(WTHPR1+WTHPX1) &
          +CFOPPs1(2,M,NZ)*(WTHPR2+WTHPX2)
        IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(3,M,NZ)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(3,M,NZ)*(WTHNR3+WTHNX3+WTHNR4+WTHNX4)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(3,M,NZ)*(WTHPR3+WTHPX3+WTHPR4+WTHPX4)
        ELSE
          WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ) &
            +CFOPCs1(5,M,NZ)*(WTHTX3+WTHTX4)
          WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ) &
            +CFOPNs1(5,M,NZ)*(WTHNX3+WTHNX4)
          WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ) &
            +CFOPPs1(5,M,NZ)*(WTHPX3+WTHPX4)
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ) &
            +CFOPCs1(5,M,NZ)*(WTHTR3+WTHTR4)*FWOODs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ) &
            +CFOPNs1(5,M,NZ)*(WTHNR3+WTHNR4)*FWOODNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ) &
            +CFOPPs1(5,M,NZ)*(WTHPR3+WTHPR4)*FWOODPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(5,M,NZ)*(WTHTR3+WTHTR4)*FWOODs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(5,M,NZ)*(WTHNR3+WTHNR4)*FWOODNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(5,M,NZ)*(WTHPR3+WTHPR4)*FWOODPs1(0)
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
      DO 6485 M=1,4
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
          +CFOPCs1(0,M,NZ)*(WTHTR0+WTHTX0) &
          +CFOPCs1(1,M,NZ)*(WTHTR1+WTHTX1) &
          +CFOPCs1(2,M,NZ)*(WTHTR2+WTHTX2)
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
          +CFOPNs1(0,M,NZ)*WTHNL0 &
          +CFOPNs1(1,M,NZ)*WTHNL1 &
          +CFOPNs1(2,M,NZ)*WTHNL2
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
          +CFOPPs1(0,M,NZ)*WTHPL0 &
          +CFOPPs1(1,M,NZ)*WTHPL1 &
          +CFOPPs1(2,M,NZ)*WTHPL2
        ZSNCs1(4,1,0,NZ)=ZSNCs1(4,1,0,NZ) &
          +CFOPNs1(0,M,NZ)*(WTHNR0+WTHNX0-WTHNL0) &
          +CFOPNs1(1,M,NZ)*(WTHNR1+WTHNX1-WTHNL1) &
          +CFOPNs1(2,M,NZ)*(WTHNR2+WTHNX2-WTHNL2)
        PSNCs1(4,1,0,NZ)=PSNCs1(4,1,0,NZ) &
          +CFOPPs1(0,M,NZ)*(WTHPR0+WTHPX0-WTHPL0) &
          +CFOPPs1(1,M,NZ)*(WTHPR1+WTHPX1-WTHPL1) &
          +CFOPPs1(2,M,NZ)*(WTHPR2+WTHPX2-WTHPL2)
        IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(3,M,NZ)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(3,M,NZ)*(WTHNL3+WTHNL4)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(3,M,NZ)*(WTHPL3+WTHPL4)
          ZSNCs1(4,1,0,NZ)=ZSNCs1(4,1,0,NZ)+CFOPNs1(3,M,NZ)*(WTHNR3+WTHNX3-WTHNL3+WTHNR4+WTHNX4-WTHNL4)
          PSNCs1(4,1,0,NZ)=PSNCs1(4,1,0,NZ)+CFOPPs1(3,M,NZ)*(WTHPR3+WTHPX3-WTHPL3+WTHPR4+WTHPX4-WTHPL4)
        ELSE
          WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ)+CFOPCs1(5,M,NZ)*(WTHTR3+WTHTX3)
          WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ)+CFOPNs1(5,M,NZ)*WTHNL3
          WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ)+CFOPPs1(5,M,NZ)*WTHPL3
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)*CFOPCs1(3,M,NZ)*(WTHTR4+WTHTX4)*FWOODs1(0)
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+CFOPNs1(3,M,NZ)*WTHNL4*FWOODNs1(0)
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+CFOPPs1(3,M,NZ)*WTHPL4*FWOODPs1(0)
          ZSNCs1(4,0,0,NZ)=ZSNCs1(4,0,0,NZ)+CFOPNs1(5,M,NZ)*(WTHNR3+WTHNX3-WTHNL3 &
            +WTHNR4+WTHNX4-WTHNL4)*FWOODNs1(0)
          PSNCs1(4,0,0,NZ)=PSNCs1(4,0,0,NZ)+CFOPPs1(5,M,NZ)*(WTHPR3+WTHPX3-WTHPL3 &
            +WTHPR4+WTHPX4-WTHPL4)*FWOODPs1(0)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(3,M,NZ)*(WTHTR4+WTHTX4)*FWOODs1(1)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(3,M,NZ)*WTHNL4*FWOODNs1(1)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(3,M,NZ)*WTHPL4*FWOODPs1(1)
          ZSNCs1(4,1,0,NZ)=ZSNCs1(4,1,0,NZ)+CFOPNs1(5,M,NZ)*(WTHNR3+WTHNX3-WTHNL3 &
            +WTHNR4+WTHNX4-WTHNL4)*FWOODNs1(1)
          PSNCs1(4,1,0,NZ)=PSNCs1(4,1,0,NZ)+CFOPPs1(5,M,NZ)*(WTHPR3+WTHPX3-WTHPL3 &
            +WTHPR4+WTHPX4-WTHPL4)*FWOODPs1(1)
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
    TCSNCs1(NZ)=TCSNCs1(NZ)+WTHTRT+WTHTXT
    TZSNCs1(NZ)=TZSNCs1(NZ)+WTHNRT+WTHNXT
    TPSNCs1(NZ)=TPSNCs1(NZ)+WTHPRT+WTHPXT
    TCSN0s1(NZ)=TCSN0s1(NZ)+WTHTRT+WTHTXT
    TZSN0s1(NZ)=TZSNCs1(NZ)+WTHNRT+WTHNXT
    TPSN0s1(NZ)=TPSNCs1(NZ)+WTHPRT+WTHPXT
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
    IHVSTs1    =>  plt_distb%IHVSTs1  , &
    JHVSTs1    =>  plt_distb%JHVSTs1  , &
    HVSTCs1    =>  plt_distb%HVSTCs1  , &
    HVSTNs1    =>  plt_distb%HVSTNs1  , &
    HVSTPs1    =>  plt_distb%HVSTPs1  , &
    XHVSTCs1   =>  plt_distb%XHVSTCs1 , &
    XHVSTNs1   =>  plt_distb%XHVSTNs1 , &
    XHVSTPs1   =>  plt_distb%XHVSTPs1 , &
    VNH3Fs1    =>  plt_distb%VNH3Fs1  , &
    VPO4Fs1    =>  plt_distb%VPO4Fs1  , &
    VCH4Fs1    =>  plt_distb%VCH4Fs1  , &
    VOXYFs1    =>  plt_distb%VOXYFs1  , &
    VN2OFs1    =>  plt_distb%VN2OFs1  , &
    VCO2Fs1    =>  plt_distb%VCO2Fs1  , &
    CNETs1     =>  plt_bgcr%CNETs1    , &
    TNBPs1     =>  plt_bgcr%TNBPs1    , &
    TRAUs1     =>  plt_bgcr%TRAUs1    , &
    RECOs1     =>  plt_bgcr%RECOs1    , &
    TCO2As1    =>  plt_bgcr%TCO2As1   , &
    TCO2Ts1    =>  plt_bgcr%TCO2Ts1   , &
    WTRVNs1    =>  plt_biom%WTRVNs1   , &
    WTRVCs1    =>  plt_biom%WTRVCs1   , &
    WTRVPs1    =>  plt_biom%WTRVPs1     &
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

  IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
    IF(IHVSTs1(NZ).NE.5)THEN
      IF(JHVSTs1(NZ).NE.2)THEN
        HVSTCs1(NZ)=HVSTCs1(NZ)+WTHTHT-WTHTRT
        HVSTNs1(NZ)=HVSTNs1(NZ)+WTHNHT-WTHNRT
        HVSTPs1(NZ)=HVSTPs1(NZ)+WTHPHT-WTHPRT
        TNBPs1=TNBPs1+WTHTRT-WTHTHT
        XHVSTCs1=XHVSTCs1+WTHTHT-WTHTRT
        XHVSTNs1=XHVSTNs1+WTHNHT-WTHNRT
        XHVSTPs1=XHVSTPs1+WTHPHT-WTHPRT
      ELSE
        WTRVCs1(NZ)=WTRVCs1(NZ)+WTHTHT-WTHTRT
        WTRVNs1(NZ)=WTRVNs1(NZ)+WTHNHT-WTHNRT
        WTRVPs1(NZ)=WTRVPs1(NZ)+WTHPHT-WTHPRT
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CNET=PFT net CO2 fixation
!     TNBP=total net biome productivity
!
    ELSE
      VCO2Fs1(NZ)=VCO2Fs1(NZ)-(1._r8-FCH4F)*(WTHTHT-WTHTRT)
      VCH4Fs1(NZ)=VCH4Fs1(NZ)-FCH4F*(WTHTHT-WTHTRT)
      VOXYFs1(NZ)=VOXYFs1(NZ)-(1._r8-FCH4F)*(WTHTHT-WTHTRT)*2.667
      VNH3Fs1(NZ)=VNH3Fs1(NZ)-WTHNHT+WTHNRT
      VN2OFs1(NZ)=VN2OFs1(NZ)-0.0
      VPO4Fs1(NZ)=VPO4Fs1(NZ)-WTHPHT+WTHPRT
      CNETs1(NZ)=CNETs1(NZ)-(1._r8-FCH4F)*(WTHTHT-WTHTRT)
      TNBPs1=TNBPs1-FCH4F*(WTHTHT-WTHTRT)
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
    HVSTCs1(NZ)=HVSTCs1(NZ)+GY*(WTHTHT-WTHTRT)
    HVSTNs1(NZ)=HVSTNs1(NZ)+WTHNHT-WTHNRT
    HVSTPs1(NZ)=HVSTPs1(NZ)+WTHPHT-WTHPRT
    TCO2Ts1(NZ)=TCO2Ts1(NZ)-GZ*(WTHTHT-WTHTRT)
    TCO2As1(NZ)=TCO2As1(NZ)-GZ*(WTHTHT-WTHTRT)
!     TNBPs1=TNBPs1+GY*(WTHTRT-WTHTHT)
!     CNETs1(NZ)=CNETs1(NZ)+GZ*(WTHTRT-WTHTHT)
    XHVSTCs1=XHVSTCs1+GY*(WTHTHT-WTHTRT)
    XHVSTNs1=XHVSTNs1+WTHNHT-WTHNRT
    XHVSTPs1=XHVSTPs1+WTHPHT-WTHPRT
    RECOs1=RECOs1-GZ*(WTHTHT-WTHTRT)
    TRAUs1=TRAUs1-GZ*(WTHTHT-WTHTRT)
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
    NUs1       =>  plt_site%NUs1      , &
    AREA3s1    =>  plt_site%AREA3s1   , &
    EHVSTs1    =>  plt_distb%EHVSTs1  , &
    FERTs1     =>  plt_distb%FERTs1   , &
    IYTYPs1    =>  plt_distb%IYTYPs1  , &
    IHVSTs1    =>  plt_distb%IHVSTs1    &
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
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
  IF(IHVSTs1(NZ).EQ.0)THEN
    WTHTR0=WTHTH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EHVSTs1(2,4,NZ))
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVSTs1(NZ).EQ.1)THEN
    WTHTR0=WTHTH0
    WTHNR0=WTHNH0
    WTHPR0=WTHPH0
    WTHTR1=WTHTH1
    WTHNR1=WTHNH1
    WTHPR1=WTHPH1
    WTHTR2=WTHTH2-WTHTG*EHVSTs1(2,2,NZ)
    WTHNR2=WTHNH2-WTHNG*EHVSTs1(2,2,NZ)
    WTHPR2=WTHPH2-WTHPG*EHVSTs1(2,2,NZ)
    WTHTR3=WTHTH3
    WTHNR3=WTHNH3
    WTHPR3=WTHPH3
    WTHTR4=WTHTH4
    WTHNR4=WTHNH4
    WTHPR4=WTHPH4
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVSTs1(NZ).EQ.2)THEN
    WTHTR0=WTHTH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EHVSTs1(2,4,NZ))
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(IHVSTs1(NZ).EQ.3)THEN
    WTHTR0=WTHTH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EHVSTs1(2,4,NZ))
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6)THEN
    WTHTR0=WTHTH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EHVSTs1(2,1,NZ)*0.5)
    WTHPR0=WTHPH0*(1._r8-EHVSTs1(2,1,NZ)*0.5)
    WTHTR1=WTHTH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EHVSTs1(2,1,NZ)*0.5)
    WTHPR1=WTHPH1*(1._r8-EHVSTs1(2,1,NZ)*0.5)
    WTHTR2=WTHTH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EHVSTs1(2,2,NZ)*0.5)
    WTHPR2=WTHPH2*(1._r8-EHVSTs1(2,2,NZ)*0.5)
    WTHTR3=WTHTH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EHVSTs1(2,3,NZ)*0.5)
    WTHPR3=WTHPH3*(1._r8-EHVSTs1(2,3,NZ)*0.5)
    WTHTR4=WTHTH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EHVSTs1(2,4,NZ)*0.5)
    WTHPR4=WTHPH4*(1._r8-EHVSTs1(2,4,NZ)*0.5)
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
    FERTs1(17)=FERTs1(17)+(WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4)/AREA3s1(NUs1)
    FERTs1(18)=FERTs1(18)+(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA3s1(NUs1)*0.5_r8
    FERTs1(3)=FERTs1(3)+(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA3s1(NUs1)*0.5_r8
    FERTs1(19)=FERTs1(19)+(WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4)/AREA3s1(NUs1)
    IYTYPs1=3
!
!     REMOVALS BY FIRE
!
!     EFIRE=combustion  of N,P relative to C
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVSTs1(2,1,EHVSTs1(2,2,EHVSTs1(2,3,EHVSTs1(2,4=fraction of
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
  ELSEIF(IHVSTs1(NZ).EQ.5)THEN
    WTHTR0=WTHTH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR0=WTHNH0*(1._r8-EFIRE(1,IHVSTs1(NZ))*EHVSTs1(2,1,NZ))
    WTHPR0=WTHPH0*(1._r8-EFIRE(2,IHVSTs1(NZ))*EHVSTs1(2,1,NZ))
    WTHNL0=WTHNH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHPL0=WTHPH0*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR1=WTHTH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHNR1=WTHNH1*(1._r8-EFIRE(1,IHVSTs1(NZ))*EHVSTs1(2,1,NZ))
    WTHPR1=WTHPH1*(1._r8-EFIRE(2,IHVSTs1(NZ))*EHVSTs1(2,1,NZ))
    WTHNL1=WTHNH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHPL1=WTHPH1*(1._r8-EHVSTs1(2,1,NZ))
    WTHTR2=WTHTH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHNR2=WTHNH2*(1._r8-EFIRE(1,IHVSTs1(NZ))*EHVSTs1(2,2,NZ))
    WTHPR2=WTHPH2*(1._r8-EFIRE(2,IHVSTs1(NZ))*EHVSTs1(2,2,NZ))
    WTHNL2=WTHNH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHPL2=WTHPH2*(1._r8-EHVSTs1(2,2,NZ))
    WTHTR3=WTHTH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHNR3=WTHNH3*(1._r8-EFIRE(1,IHVSTs1(NZ))*EHVSTs1(2,3,NZ))
    WTHPR3=WTHPH3*(1._r8-EFIRE(2,IHVSTs1(NZ))*EHVSTs1(2,3,NZ))
    WTHNL3=WTHNH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHPL3=WTHPH3*(1._r8-EHVSTs1(2,3,NZ))
    WTHTR4=WTHTH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHNR4=WTHNH4*(1._r8-EFIRE(1,IHVSTs1(NZ))*EHVSTs1(2,4,NZ))
    WTHPR4=WTHPH4*(1._r8-EFIRE(2,IHVSTs1(NZ))*EHVSTs1(2,4,NZ))
    WTHNL4=WTHNH4*(1._r8-EHVSTs1(2,4,NZ))
    WTHPL4=WTHPH4*(1._r8-EHVSTs1(2,4,NZ))
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
    JHVSTs1    =>  plt_distb%JHVSTs1     , &
    IDAYHs1    =>  plt_distb%IDAYHs1     , &
    IDAY0s1    =>  plt_distb%IDAY0s1     , &
    ITILLs1    =>  plt_distb%ITILLs1     , &
    IYR0s1     =>  plt_distb%IYR0s1      , &
    IYRHs1     =>  plt_distb%IYRHs1      , &
    XCORPs1    =>  plt_distb%XCORPs1     , &
    CFOPCs1    =>  plt_soilchem%CFOPCs1  , &
    CFOPNs1    =>  plt_soilchem%CFOPNs1  , &
    CFOPPs1    =>  plt_soilchem%CFOPPs1  , &
    CH4Ps1     =>  plt_rbgc%CH4Ps1       , &
    CH4As1     =>  plt_rbgc%CH4As1       , &
    H2GPs1     =>  plt_rbgc%H2GPs1       , &
    H2GAs1     =>  plt_rbgc%H2GAs1       , &
    OXYPs1     =>  plt_rbgc%OXYPs1       , &
    CO2Ps1     =>  plt_rbgc%CO2Ps1       , &
    OXYAs1     =>  plt_rbgc%OXYAs1       , &
    CO2As1     =>  plt_rbgc%CO2As1       , &
    UVOLOs1    =>  plt_ew%UVOLOs1        , &
    VOLWPs1    =>  plt_ew%VOLWPs1        , &
    VHCPCs1    =>  plt_ew%VHCPCs1        , &
    PSILTs1    =>  plt_ew%PSILTs1        , &
    PPXs1      =>  plt_site%PPXs1        , &
    CPOOLRs1   =>  plt_biom%CPOOLRs1     , &
    ZPOOLRs1   =>  plt_biom%ZPOOLRs1     , &
    PPOOLRs1   =>  plt_biom%PPOOLRs1     , &
    WSRTLs1    =>  plt_biom%WSRTLs1      , &
    WTRTDs1    =>  plt_biom%WTRTDs1      , &
    WTRTLs1    =>  plt_biom%WTRTLs1      , &
    WGLFLPs1   =>  plt_biom%WGLFLPs1     , &
    WTLFBPs1   =>  plt_biom%WTLFBPs1     , &
    WTLFBs1    =>  plt_biom%WTLFBs1      , &
    WTGRBNs1   =>  plt_biom%WTGRBNs1     , &
    WTGRBs1    =>  plt_biom%WTGRBs1      , &
    WTEARBs1   =>  plt_biom%WTEARBs1     , &
    CPOLNBs1   =>  plt_biom%CPOLNBs1     , &
    PPOOLs1    =>  plt_biom%PPOOLs1      , &
    ZPOOLs1    =>  plt_biom%ZPOOLs1      , &
    CPOOLs1    =>  plt_biom%CPOOLs1      , &
    WTLFBNs1   =>  plt_biom%WTLFBNs1     , &
    WTSHTNs1   =>  plt_biom%WTSHTNs1     , &
    WTHSKBs1   =>  plt_biom%WTHSKBs1     , &
    WTRSVBs1   =>  plt_biom%WTRSVBs1     , &
    WTNDBs1    =>  plt_biom%WTNDBs1      , &
    WTSHTBs1   =>  plt_biom%WTSHTBs1     , &
    PPOLNBs1   =>  plt_biom%PPOLNBs1     , &
    ZPOLNBs1   =>  plt_biom%ZPOLNBs1     , &
    WTSTKBs1   =>  plt_biom%WTSTKBs1     , &
    WTEABNs1   =>  plt_biom%WTEABNs1     , &
    WTSHBNs1   =>  plt_biom%WTSHBNs1     , &
    WTSTBNs1   =>  plt_biom%WTSTBNs1     , &
    WTHSBNs1   =>  plt_biom%WTHSBNs1     , &
    WTNDBNs1   =>  plt_biom%WTNDBNs1     , &
    WTRSBNs1   =>  plt_biom%WTRSBNs1     , &
    WTSHTPs1   =>  plt_biom%WTSHTPs1     , &
    WGLFLNs1   =>  plt_biom%WGLFLNs1     , &
    WTSHBPs1   =>  plt_biom%WTSHBPs1     , &
    WTNDBPs1   =>  plt_biom%WTNDBPs1     , &
    WGLFLs1    =>  plt_biom%WGLFLs1      , &
    WGSHPs1    =>  plt_biom%WGSHPs1      , &
    WTSHEBs1   =>  plt_biom%WTSHEBs1     , &
    WTEABPs1   =>  plt_biom%WTEABPs1     , &
    WTHSBPs1   =>  plt_biom%WTHSBPs1     , &
    WTRSBPs1   =>  plt_biom%WTRSBPs1     , &
    WTSTBPs1   =>  plt_biom%WTSTBPs1     , &
    WGLFs1     =>  plt_biom%WGLFs1       , &
    WTSTXNs1   =>  plt_biom%WTSTXNs1     , &
    WTLSBs1    =>  plt_biom%WTLSBs1      , &
    WTGRBPs1   =>  plt_biom%WTGRBPs1     , &
    WTSTXBs1   =>  plt_biom%WTSTXBs1     , &
    WVSTKBs1   =>  plt_biom%WVSTKBs1     , &
    WTSTXPs1   =>  plt_biom%WTSTXPs1     , &
    WGLFPs1    =>  plt_biom%WGLFPs1      , &
    WSSHEs1    =>  plt_biom%WSSHEs1      , &
    WGSHEs1    =>  plt_biom%WGSHEs1      , &
    WTRT1s1    =>  plt_biom%WTRT1s1      , &
    WTRT1Ns1   =>  plt_biom%WTRT1Ns1     , &
    WTRT1Ps1   =>  plt_biom%WTRT1Ps1     , &
    WSLFs1     =>  plt_biom%WSLFs1       , &
    WGSHNs1    =>  plt_biom%WGSHNs1      , &
    WGNODEs1   =>  plt_biom%WGNODEs1     , &
    WGLFNs1    =>  plt_biom%WGLFNs1      , &
    WGNODNs1   =>  plt_biom%WGNODNs1     , &
    WGNODPs1   =>  plt_biom%WGNODPs1     , &
    WVSTKs1    =>  plt_biom%WVSTKs1      , &
    RTWT1s1    =>  plt_biom%RTWT1s1      , &
    RTWT1Ns1   =>  plt_biom%RTWT1Ns1     , &
    RTWT1Ps1   =>  plt_biom%RTWT1Ps1     , &
    WTRVCs1    =>  plt_biom%WTRVCs1      , &
    WTRVNs1    =>  plt_biom%WTRVNs1      , &
    WTRVPs1    =>  plt_biom%WTRVPs1      , &
    WTLSs1     =>  plt_biom%WTLSs1       , &
    WTRT2s1    =>  plt_biom%WTRT2s1      , &
    WTRT2Ns1   =>  plt_biom%WTRT2Ns1     , &
    WTRT2Ps1   =>  plt_biom%WTRT2Ps1     , &
    CPOOLNs1   =>  plt_biom%CPOOLNs1     , &
    ZPOOLNs1   =>  plt_biom%ZPOOLNs1     , &
    PPOOLNs1   =>  plt_biom%PPOOLNs1     , &
    WTNDLs1    =>  plt_biom%WTNDLs1      , &
    WTNDLNs1   =>  plt_biom%WTNDLNs1     , &
    WTNDLPs1   =>  plt_biom%WTNDLPs1     , &
    GRWTBs1    =>  plt_allom%GRWTBs1     , &
    FWOODPs1   =>  plt_allom%FWOODPs1    , &
    FWOODNs1   =>  plt_allom%FWOODNs1    , &
    FWOODs1    =>  plt_allom%FWOODs1     , &
    FWODLPs1   =>  plt_allom%FWODLPs1    , &
    FWODSNs1   =>  plt_allom%FWODSNs1    , &
    FWODBs1    =>  plt_allom%FWODBs1     , &
    FWODRPs1   =>  plt_allom%FWODRPs1    , &
    FWODSPs1   =>  plt_allom%FWODSPs1    , &
    FWODLNs1   =>  plt_allom%FWODLNs1    , &
    FWODRs1    =>  plt_allom%FWODRs1     , &
    FWODRNs1   =>  plt_allom%FWODRNs1    , &
    IDTHBs1    =>  plt_pheno%IDTHBs1     , &
    ISTYPs1    =>  plt_pheno%ISTYPs1     , &
    IDTHs1     =>  plt_pheno%IDTHs1      , &
    IGTYPs1    =>  plt_pheno%IGTYPs1     , &
    IBTYPs1    =>  plt_pheno%IBTYPs1     , &
    IWTYPs1    =>  plt_pheno%IWTYPs1     , &
    IDTHRs1    =>  plt_pheno%IDTHRs1     , &
    IDTHPs1    =>  plt_pheno%IDTHPs1     , &
    HCOBs1     =>  plt_photo%HCOBs1      , &
    CO2Bs1     =>  plt_photo%CO2Bs1      , &
    CPOOL3s1   =>  plt_photo%CPOOL3s1    , &
    CPOOL4s1   =>  plt_photo%CPOOL4s1    , &
    NJs1       =>  plt_site%NJs1         , &
    PPs1       =>  plt_site%PPs1         , &
    IYRCs1     =>  plt_site%IYRCs1       , &
    ZNOONs1    =>  plt_site%ZNOONs1      , &
    VOLWOUs1   =>  plt_site%VOLWOUs1     , &
    NUs1       =>  plt_site%NUs1         , &
    CSNCs1     =>  plt_bgcr%CSNCs1       , &
    ZSNCs1     =>  plt_bgcr%ZSNCs1       , &
    PSNCs1     =>  plt_bgcr%PSNCs1       , &
    RH2GZs1    =>  plt_bgcr%RH2GZs1      , &
    RNH3Zs1    =>  plt_bgcr%RNH3Zs1      , &
    RN2OZs1    =>  plt_bgcr%RN2OZs1      , &
    RCO2Zs1    =>  plt_bgcr%RCO2Zs1      , &
    RCH4Zs1    =>  plt_bgcr%RCH4Zs1      , &
    ROXYZs1    =>  plt_bgcr%ROXYZs1      , &
    ZH3Ps1     =>  plt_rbgc%ZH3Ps1       , &
    Z2OPs1     =>  plt_rbgc%Z2OPs1       , &
    ZH3As1     =>  plt_rbgc%ZH3As1       , &
    RCO2As1    =>  plt_rbgc%RCO2As1      , &
    Z2OAs1     =>  plt_rbgc%Z2OAs1       , &
    RCO2Ms1    =>  plt_rbgc%RCO2Ms1      , &
    RCO2Ns1    =>  plt_rbgc%RCO2Ns1      , &
    FRADPs1    =>  plt_rad%FRADPs1       , &
    RTLG1s1    =>  plt_morph%RTLG1s1     , &
    RTVLWs1    =>  plt_morph%RTVLWs1     , &
    RTARPs1    =>  plt_morph%RTARPs1     , &
    RTVLPs1    =>  plt_morph%RTVLPs1     , &
    RTDNPs1    =>  plt_morph%RTDNPs1     , &
    RTN1s1     =>  plt_morph%RTN1s1      , &
    INTYPs1    =>  plt_morph%INTYPs1     , &
    RTLGPs1    =>  plt_morph%RTLGPs1     , &
    RTN2s1     =>  plt_morph%RTN2s1      , &
    RTLG2s1    =>  plt_morph%RTLG2s1     , &
    RTNLs1     =>  plt_morph%RTNLs1      , &
    NGs1       =>  plt_morph%NGs1        , &
    MYs1       =>  plt_morph%MYs1        , &
    NRTs1      =>  plt_morph%NRTs1       , &
    NBRs1      =>  plt_morph%NBRs1       , &
    ARLF1s1    =>  plt_morph%ARLF1s1     , &
    ARLFBs1    =>  plt_morph%ARLFBs1     , &
    GRNXBs1    =>  plt_morph%GRNXBs1     , &
    GRNOBs1    =>  plt_morph%GRNOBs1     , &
    ARLFLs1    =>  plt_morph%ARLFLs1       &
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
  IF(J.EQ.INT(ZNOONs1).AND.(IBTYPs1(NZ).EQ.0 &
    .OR.IGTYPs1(NZ).LE.1).AND.(I.NE.IDAY0s1(NZ) &
    .OR.IYRCs1.NE.IYR0s1(NZ)))THEN
    IF(ITILLs1.LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.IDAY0s1(NZ).OR.IYRCs1.GT.IYR0s1(NZ))THEN
        XHVST=XCORPs1
        PPXs1(NZ)=PPXs1(NZ)*XHVST
        PPs1(NZ)=PPs1(NZ)*XHVST
        FRADPs1(NZ)=FRADPs1(NZ)*XHVST
        VHCPCs1(NZ)=VHCPCs1(NZ)*XHVST
        WTLSs1(NZ)=0._r8
        WVSTKs1(NZ)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
        DO 8975 NB=1,NBRs1(NZ)
          IF(IDTHBs1(NB,NZ).EQ.0)THEN
            IF(PPs1(NZ).LE.0.0)IDTHBs1(NB,NZ)=1
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
            DO 6380 M=1,4
              CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                *(CFOPCs1(0,M,NZ)*(CPOOLs1(NB,NZ)+CPOLNBs1(NB,NZ) &
                +CPOOLK(NB,NZ)+WTRSVBs1(NB,NZ)) &
                +CFOPCs1(1,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(1) &
                +WTNDBs1(NB,NZ)) &
                +CFOPCs1(2,M,NZ)*(WTSHEBs1(NB,NZ)*FWODBs1(1) &
                +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ)))
              CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPCs1(5,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(0) &
                +WTSHEBs1(NB,NZ)*FWODBs1(0))
              ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                *(CFOPNs1(0,M,NZ)*(ZPOOLs1(NB,NZ)+ZPOLNBs1(NB,NZ) &
                +WTRSBNs1(NB,NZ)) &
                +CFOPNs1(1,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(1) &
                +WTNDBNs1(NB,NZ)) &
                +CFOPNs1(2,M,NZ)*(WTSHBNs1(NB,NZ)*FWODSNs1(1) &
                +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ)))
              ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPNs1(5,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(0) &
                +WTSHBNs1(NB,NZ)*FWODSNs1(0))
              PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                *(CFOPPs1(0,M,NZ)*(PPOOLs1(NB,NZ)+PPOLNBs1(NB,NZ) &
                +WTRSBPs1(NB,NZ)) &
                +CFOPPs1(1,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(1) &
                +WTNDBPs1(NB,NZ)) &
                +CFOPPs1(2,M,NZ)*(WTSHBPs1(NB,NZ)*FWODSPs1(1) &
                +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ)))
              PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPPs1(5,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(0) &
                +WTSHBPs1(NB,NZ)*FWODSPs1(0))
              IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
                WTRVCs1(NZ)=WTRVCs1(NZ)+(1._r8-XHVST) &
                  *CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
                WTRVNs1(NZ)=WTRVNs1(NZ)+(1._r8-XHVST) &
                  *CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
                WTRVPs1(NZ)=WTRVPs1(NZ)+(1._r8-XHVST) &
                  *CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
              ELSE
                CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                  *CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
                ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                  *CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
                PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                  *CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
              ENDIF
              CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPCs1(5,M,NZ)*WTSTKBs1(NB,NZ)*FWOODs1(0)
              ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPNs1(5,M,NZ)*WTSTBNs1(NB,NZ)*FWOODNs1(0)
              PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ)+(1._r8-XHVST) &
                *CFOPPs1(5,M,NZ)*WTSTBPs1(NB,NZ)*FWOODPs1(0)
              CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                *CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)*FWOODs1(1)
              ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                *CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)*FWOODNs1(1)
              PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+(1._r8-XHVST) &
                *CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)*FWOODPs1(1)
6380        CONTINUE
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!
            CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)*XHVST
            CPOOLK(NB,NZ)=CPOOLK(NB,NZ)*XHVST
            ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)*XHVST
            PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)*XHVST
            CPOLNBs1(NB,NZ)=CPOLNBs1(NB,NZ)*XHVST
            ZPOLNBs1(NB,NZ)=ZPOLNBs1(NB,NZ)*XHVST
            PPOLNBs1(NB,NZ)=PPOLNBs1(NB,NZ)*XHVST
            WTSHTBs1(NB,NZ)=WTSHTBs1(NB,NZ)*XHVST
            WTLFBs1(NB,NZ)=WTLFBs1(NB,NZ)*XHVST
            WTNDBs1(NB,NZ)=WTNDBs1(NB,NZ)*XHVST
            WTSHEBs1(NB,NZ)=WTSHEBs1(NB,NZ)*XHVST
            WTSTKBs1(NB,NZ)=WTSTKBs1(NB,NZ)*XHVST
            WVSTKBs1(NB,NZ)=WVSTKBs1(NB,NZ)*XHVST
            WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)*XHVST
            WTHSKBs1(NB,NZ)=WTHSKBs1(NB,NZ)*XHVST
            WTEARBs1(NB,NZ)=WTEARBs1(NB,NZ)*XHVST
            WTGRBs1(NB,NZ)=WTGRBs1(NB,NZ)*XHVST
            WTSHTNs1(NB,NZ)=WTSHTNs1(NB,NZ)*XHVST
            WTLFBNs1(NB,NZ)=WTLFBNs1(NB,NZ)*XHVST
            WTNDBNs1(NB,NZ)=WTNDBNs1(NB,NZ)*XHVST
            WTSHBNs1(NB,NZ)=WTSHBNs1(NB,NZ)*XHVST
            WTSTBNs1(NB,NZ)=WTSTBNs1(NB,NZ)*XHVST
            WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)*XHVST
            WTHSBNs1(NB,NZ)=WTHSBNs1(NB,NZ)*XHVST
            WTEABNs1(NB,NZ)=WTEABNs1(NB,NZ)*XHVST
            WTGRBNs1(NB,NZ)=WTGRBNs1(NB,NZ)*XHVST
            WTSHTPs1(NB,NZ)=WTSHTPs1(NB,NZ)*XHVST
            WTLFBPs1(NB,NZ)=WTLFBPs1(NB,NZ)*XHVST
            WTNDBPs1(NB,NZ)=WTNDBPs1(NB,NZ)*XHVST
            WTSHBPs1(NB,NZ)=WTSHBPs1(NB,NZ)*XHVST
            WTSTBPs1(NB,NZ)=WTSTBPs1(NB,NZ)*XHVST
            WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)*XHVST
            WTHSBPs1(NB,NZ)=WTHSBPs1(NB,NZ)*XHVST
            WTEABPs1(NB,NZ)=WTEABPs1(NB,NZ)*XHVST
            WTGRBPs1(NB,NZ)=WTGRBPs1(NB,NZ)*XHVST
            GRNXBs1(NB,NZ)=GRNXBs1(NB,NZ)*XHVST
            GRNOBs1(NB,NZ)=GRNOBs1(NB,NZ)*XHVST
            GRWTBs1(NB,NZ)=GRWTBs1(NB,NZ)*XHVST
            ARLFBs1(NB,NZ)=ARLFBs1(NB,NZ)*XHVST
            WTLSBs1(NB,NZ)=AMAX1(0.0,WTLFBs1(NB,NZ)+WTSHEBs1(NB,NZ))
            WTLSs1(NZ)=WTLSs1(NZ)+WTLSBs1(NB,NZ)
            WTSTXBs1(NB,NZ)=WTSTXBs1(NB,NZ)*XHVST
            WTSTXNs1(NB,NZ)=WTSTXNs1(NB,NZ)*XHVST
            WTSTXPs1(NB,NZ)=WTSTXPs1(NB,NZ)*XHVST
            WVSTKs1(NZ)=WVSTKs1(NZ)+WVSTKBs1(NB,NZ)
            DO 8970 K=0,JNODS1
              IF(K.NE.0)THEN
                CPOOL3s1(K,NB,NZ)=CPOOL3s1(K,NB,NZ)*XHVST
                CPOOL4s1(K,NB,NZ)=CPOOL4s1(K,NB,NZ)*XHVST
                CO2Bs1(K,NB,NZ)=CO2Bs1(K,NB,NZ)*XHVST
                HCOBs1(K,NB,NZ)=HCOBs1(K,NB,NZ)*XHVST
              ENDIF
              ARLF1s1(K,NB,NZ)=ARLF1s1(K,NB,NZ)*XHVST
              WGLFs1(K,NB,NZ)=WGLFs1(K,NB,NZ)*XHVST
              WSLFs1(K,NB,NZ)=WSLFs1(K,NB,NZ)*XHVST
!     HTSHEs1(K,NB,NZ)=HTSHEs1(K,NB,NZ)*XHVST
              WGSHEs1(K,NB,NZ)=WGSHEs1(K,NB,NZ)*XHVST
              WSSHEs1(K,NB,NZ)=WSSHEs1(K,NB,NZ)*XHVST
!     HTNODEs1(K,NB,NZ)=HTNODEs1(K,NB,NZ)*XHVST
!     HTNODXs1(K,NB,NZ)=HTNODXs1(K,NB,NZ)*XHVST
              WGNODEs1(K,NB,NZ)=WGNODEs1(K,NB,NZ)*XHVST
              WGLFNs1(K,NB,NZ)=WGLFNs1(K,NB,NZ)*XHVST
              WGSHNs1(K,NB,NZ)=WGSHNs1(K,NB,NZ)*XHVST
              WGNODNs1(K,NB,NZ)=WGNODNs1(K,NB,NZ)*XHVST
              WGLFPs1(K,NB,NZ)=WGLFPs1(K,NB,NZ)*XHVST
              WGSHPs1(K,NB,NZ)=WGSHPs1(K,NB,NZ)*XHVST
              WGNODPs1(K,NB,NZ)=WGNODPs1(K,NB,NZ)*XHVST
              DO 8965 L=1,JC1
                ARLFLs1(L,K,NB,NZ)=ARLFLs1(L,K,NB,NZ)*XHVST
                WGLFLs1(L,K,NB,NZ)=WGLFLs1(L,K,NB,NZ)*XHVST
                WGLFLNs1(L,K,NB,NZ)=WGLFLNs1(L,K,NB,NZ)*XHVST
                WGLFLPs1(L,K,NB,NZ)=WGLFLPs1(L,K,NB,NZ)*XHVST
8965          CONTINUE
8970        CONTINUE
          ENDIF
8975    CONTINUE
!
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOUs1,UVOLO=accumulated water loss for water balance calculation
!
        VOLWPX=VOLWPs1(NZ)
        WVPLT=AMAX1(0.0_r8,WTLSs1(NZ)+WVSTKs1(NZ))
        APSILT=ABS(PSILTs1(NZ))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWPs1(NZ)=1.0E-06_r8*WVPLT/FDM
        VOLWOUs1=VOLWOUs1+VOLWPX-VOLWPs1(NZ)
        UVOLOs1=UVOLOs1+VOLWPX-VOLWPs1(NZ)
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
        IF(PPs1(NZ).LE.0.0)THEN
          IDTHRs1(NZ)=1
          IDTHPs1(NZ)=1
          IDTHs1(NZ)=1
          JHVSTs1(NZ)=1
          IDAYHs1(NZ)=I
          IYRHs1(NZ)=IYRCs1
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
        DO 8985 N=1,MYs1(NZ)
          DO 8980 L=NUs1,NJs1
            DO 6385 M=1,4
              CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                *CFOPCs1(0,M,NZ)*CPOOLRs1(N,L,NZ)
              ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                *CFOPNs1(0,M,NZ)*ZPOOLRs1(N,L,NZ)
              PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                *CFOPPs1(0,M,NZ)*PPOOLRs1(N,L,NZ)
              DO NR=1,NRTs1(NZ)
                CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+(1._r8-XHVST) &
                  *CFOPCs1(5,M,NZ)*(WTRT1s1(N,L,NR,NZ) &
                  +WTRT2s1(N,L,NR,NZ))*FWODRs1(0)
                ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+(1._r8-XHVST) &
                  *CFOPNs1(5,M,NZ)*(WTRT1Ns1(N,L,NR,NZ) &
                  +WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(0)
                PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+(1._r8-XHVST) &
                  *CFOPPs1(5,M,NZ)*(WTRT1Ps1(N,L,NR,NZ) &
                  +WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(0)
                CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *CFOPCs1(4,M,NZ)*(WTRT1s1(N,L,NR,NZ) &
                  +WTRT2s1(N,L,NR,NZ))*FWODRs1(1)
                ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *CFOPNs1(4,M,NZ)*(WTRT1Ns1(N,L,NR,NZ) &
                  +WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(1)
                PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *CFOPPs1(4,M,NZ)*(WTRT1Ps1(N,L,NR,NZ) &
                  +WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(1)
              ENDDO
6385        CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Zs1(NZ)=RCO2Zs1(NZ)-(1._r8-XHVST) &
              *(CO2As1(N,L,NZ)+CO2Ps1(N,L,NZ))
            ROXYZs1(NZ)=ROXYZs1(NZ)-(1._r8-XHVST) &
              *(OXYAs1(N,L,NZ)+OXYPs1(N,L,NZ))
            RCH4Zs1(NZ)=RCH4Zs1(NZ)-(1._r8-XHVST) &
              *(CH4As1(N,L,NZ)+CH4Ps1(N,L,NZ))
            RN2OZs1(NZ)=RN2OZs1(NZ)-(1._r8-XHVST) &
              *(Z2OAs1(N,L,NZ)+Z2OPs1(N,L,NZ))
            RNH3Zs1(NZ)=RNH3Zs1(NZ)-(1._r8-XHVST) &
              *(ZH3As1(N,L,NZ)+ZH3Ps1(N,L,NZ))
            RH2GZs1(NZ)=RH2GZs1(NZ)-(1._r8-XHVST) &
              *(H2GAs1(N,L,NZ)+H2GPs1(N,L,NZ))
            CO2As1(N,L,NZ)=XHVST*CO2As1(N,L,NZ)
            OXYAs1(N,L,NZ)=XHVST*OXYAs1(N,L,NZ)
            CH4As1(N,L,NZ)=XHVST*CH4As1(N,L,NZ)
            Z2OAs1(N,L,NZ)=XHVST*Z2OAs1(N,L,NZ)
            ZH3As1(N,L,NZ)=XHVST*ZH3As1(N,L,NZ)
            H2GAs1(N,L,NZ)=XHVST*H2GAs1(N,L,NZ)
            CO2Ps1(N,L,NZ)=XHVST*CO2Ps1(N,L,NZ)
            OXYPs1(N,L,NZ)=XHVST*OXYPs1(N,L,NZ)
            CH4Ps1(N,L,NZ)=XHVST*CH4Ps1(N,L,NZ)
            Z2OPs1(N,L,NZ)=XHVST*Z2OPs1(N,L,NZ)
            ZH3Ps1(N,L,NZ)=XHVST*ZH3Ps1(N,L,NZ)
            H2GPs1(N,L,NZ)=XHVST*H2GPs1(N,L,NZ)
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
            DO 8960 NR=1,NRTs1(NZ)
              WTRT1s1(N,L,NR,NZ)=WTRT1s1(N,L,NR,NZ)*XHVST
              WTRT2s1(N,L,NR,NZ)=WTRT2s1(N,L,NR,NZ)*XHVST
              WTRT1Ns1(N,L,NR,NZ)=WTRT1Ns1(N,L,NR,NZ)*XHVST
              WTRT2Ns1(N,L,NR,NZ)=WTRT2Ns1(N,L,NR,NZ)*XHVST
              WTRT1Ps1(N,L,NR,NZ)=WTRT1Ps1(N,L,NR,NZ)*XHVST
              WTRT2Ps1(N,L,NR,NZ)=WTRT2Ps1(N,L,NR,NZ)*XHVST
              RTWT1s1(N,NR,NZ)=RTWT1s1(N,NR,NZ)*XHVST
              RTWT1Ns1(N,NR,NZ)=RTWT1Ns1(N,NR,NZ)*XHVST
              RTWT1Ps1(N,NR,NZ)=RTWT1Ps1(N,NR,NZ)*XHVST
              RTLG1s1(N,L,NR,NZ)=RTLG1s1(N,L,NR,NZ)*XHVST
              RTLG2s1(N,L,NR,NZ)=RTLG2s1(N,L,NR,NZ)*XHVST
              RTN2s1(N,L,NR,NZ)=RTN2s1(N,L,NR,NZ)*XHVST
8960        CONTINUE
            CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)*XHVST
            ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)*XHVST
            PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)*XHVST
            WTRTLs1(N,L,NZ)=WTRTLs1(N,L,NZ)*XHVST
            WTRTDs1(N,L,NZ)=WTRTDs1(N,L,NZ)*XHVST
            WSRTLs1(N,L,NZ)=WSRTLs1(N,L,NZ)*XHVST
            RTN1s1(N,L,NZ)=RTN1s1(N,L,NZ)*XHVST
            RTNLs1(N,L,NZ)=RTNLs1(N,L,NZ)*XHVST
            RTLGPs1(N,L,NZ)=RTLGPs1(N,L,NZ)*XHVST
            RTDNPs1(N,L,NZ)=RTDNPs1(N,L,NZ)*XHVST
            RTVLPs1(N,L,NZ)=RTVLPs1(N,L,NZ)*XHVST
            RTVLWs1(N,L,NZ)=RTVLWs1(N,L,NZ)*XHVST
            RTARPs1(N,L,NZ)=RTARPs1(N,L,NZ)*XHVST
            RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)*XHVST
            RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)*XHVST
            RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)*XHVST
!
!     LITTERFALL AND STATE VARIABLES FOR NODULES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(INTYPs1(NZ).NE.0.AND.N.EQ.1)THEN
              DO 6395 M=1,4
                CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPCs1(4,M,NZ)*WTNDLs1(L,NZ) &
                  +CFOPCs1(0,M,NZ)*CPOOLNs1(L,NZ))
                ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPNs1(4,M,NZ)*WTNDLNs1(L,NZ) &
                  +CFOPNs1(0,M,NZ)*ZPOOLNs1(L,NZ))
                PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPPs1(4,M,NZ)*WTNDLPs1(L,NZ) &
                  +CFOPPs1(0,M,NZ)*PPOOLNs1(L,NZ))
6395          CONTINUE
              WTNDLs1(L,NZ)=WTNDLs1(L,NZ)*XHVST
              WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)*XHVST
              WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)*XHVST
              CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)*XHVST
              ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)*XHVST
              PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)*XHVST
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
        DO 6400 M=1,4
          CSNCs1(M,0,NGs1(NZ),NZ)=CSNCs1(M,0,NGs1(NZ),NZ) &
            +((1._r8-XHVST)*CFOPCs1(0,M,NZ)*WTRVCs1(NZ))*FWOODs1(0)
          ZSNCs1(M,0,NGs1(NZ),NZ)=ZSNCs1(M,0,NGs1(NZ),NZ) &
            +((1._r8-XHVST)*CFOPNs1(0,M,NZ)*WTRVNs1(NZ))*FWOODNs1(0)
          PSNCs1(M,0,NGs1(NZ),NZ)=PSNCs1(M,0,NGs1(NZ),NZ) &
            +((1._r8-XHVST)*CFOPPs1(0,M,NZ)*WTRVPs1(NZ))*FWOODPs1(0)
          CSNCs1(M,1,NGs1(NZ),NZ)=CSNCs1(M,1,NGs1(NZ),NZ) &
            +((1._r8-XHVST)*CFOPCs1(0,M,NZ)*WTRVCs1(NZ))*FWOODs1(1)
          ZSNCs1(M,1,NGs1(NZ),NZ)=ZSNCs1(M,1,NGs1(NZ),NZ) &
            +((1._r8-XHVST)*CFOPNs1(0,M,NZ)*WTRVNs1(NZ))*FWOODNs1(1)
          PSNCs1(M,1,NGs1(NZ),NZ)=PSNCs1(M,1,NGs1(NZ),NZ) &
            +((1._r8-XHVST)*CFOPPs1(0,M,NZ)*WTRVPs1(NZ))*FWOODPs1(1)
6400    CONTINUE
        WTRVCs1(NZ)=WTRVCs1(NZ)*XHVST
        WTRVNs1(NZ)=WTRVNs1(NZ)*XHVST
        WTRVPs1(NZ)=WTRVPs1(NZ)*XHVST
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
    HVSTs1     =>  plt_distb%HVSTs1    , &
    EHVSTs1    =>  plt_distb%EHVSTs1   , &
    DCORPs1    =>  plt_distb%DCORPs1   , &
    THINs1     =>  plt_distb%THINs1    , &
    ITILLs1    =>  plt_distb%ITILLs1   , &
    IHVSTs1    =>  plt_distb%IHVSTs1   , &
    JHVSTs1    =>  plt_distb%JHVSTs1   , &
    VPO4Fs1    =>  plt_distb%VPO4Fs1   , &
    VN2OFs1    =>  plt_distb%VN2OFs1   , &
    VNH3Fs1    =>  plt_distb%VNH3Fs1   , &
    VOXYFs1    =>  plt_distb%VOXYFs1   , &
    VCH4Fs1    =>  plt_distb%VCH4Fs1   , &
    VCO2Fs1    =>  plt_distb%VCO2Fs1   , &
    UVOLOs1    =>  plt_ew%UVOLOs1      , &
    VOLWPs1    =>  plt_ew%VOLWPs1      , &
    PSILTs1    =>  plt_ew%PSILTs1      , &
    HCOBs1     =>  plt_photo%HCOBs1    , &
    CO2Bs1     =>  plt_photo%CO2Bs1    , &
    CPOOL3s1   =>  plt_photo%CPOOL3s1  , &
    CPOOL4s1   =>  plt_photo%CPOOL4s1  , &
    PPs1       =>  plt_site%PPs1       , &
    PPIs1      =>  plt_site%PPIs1      , &
    PPXs1      =>  plt_site%PPXs1      , &
    NUs1       =>  plt_site%NUs1       , &
    NJs1       => plt_site%NJs1        , &
    ZNOONs1    => plt_site%ZNOONs1     , &
    ZEROSs1    => plt_site%ZEROSs1     , &
    ZEROs1     => plt_site%ZEROs1      , &
    AREA3s1    => plt_site%AREA3s1     , &
    VOLWOUs1   => plt_site%VOLWOUs1    , &
    CPOOLRs1   => plt_biom%CPOOLRs1    , &
    ZPOOLRs1   => plt_biom%ZPOOLRs1    , &
    PPOOLRs1   => plt_biom%PPOOLRs1    , &
    WSRTLs1    => plt_biom%WSRTLs1     , &
    WTRTDs1    => plt_biom%WTRTDs1     , &
    WTRTLs1    => plt_biom%WTRTLs1     , &
    WTRVCs1    => plt_biom%WTRVCs1     , &
    WVSTKs1    => plt_biom%WVSTKs1     , &
    WTLSs1     => plt_biom%WTLSs1      , &
    WTRVPs1    => plt_biom%WTRVPs1     , &
    WTRVNs1    => plt_biom%WTRVNs1     , &
    WVSTKBs1   => plt_biom%WVSTKBs1    , &
    WTNDBs1    => plt_biom%WTNDBs1     , &
    WTNDBNs1   => plt_biom%WTNDBNs1    , &
    WTRSVs1    => plt_biom%WTRSVs1     , &
    WTGRBs1    => plt_biom%WTGRBs1     , &
    WTGRBNs1   => plt_biom%WTGRBNs1    , &
    WTSTKBs1   => plt_biom%WTSTKBs1    , &
    WTSHTNs1   => plt_biom%WTSHTNs1    , &
    WTSHTPs1   => plt_biom%WTSHTPs1    , &
    WTHSBPs1   => plt_biom%WTHSBPs1    , &
    WTSHTBs1   => plt_biom%WTSHTBs1    , &
    WTSTBNs1   => plt_biom%WTSTBNs1    , &
    WTHSKBs1   => plt_biom%WTHSKBs1    , &
    WTHSBNs1   => plt_biom%WTHSBNs1    , &
    WTEABPs1   => plt_biom%WTEABPs1    , &
    WTEARBs1   => plt_biom%WTEARBs1    , &
    WTEABNs1   => plt_biom%WTEABNs1    , &
    WTGRBPs1   => plt_biom%WTGRBPs1    , &
    WTSTXNs1   => plt_biom%WTSTXNs1    , &
    WGNODEs1   => plt_biom%WGNODEs1    , &
    WTRSBNs1   => plt_biom%WTRSBNs1    , &
    WTSTXPs1   => plt_biom%WTSTXPs1    , &
    WTNDBPs1   => plt_biom%WTNDBPs1    , &
    WTRSBPs1   => plt_biom%WTRSBPs1    , &
    WTLSBs1    => plt_biom%WTLSBs1     , &
    WTSTBPs1   => plt_biom%WTSTBPs1    , &
    WGNODNs1   => plt_biom%WGNODNs1    , &
    WGNODPs1   => plt_biom%WGNODPs1    , &
    WTRSVBs1   => plt_biom%WTRSVBs1    , &
    WTSTXBs1   => plt_biom%WTSTXBs1    , &
    CPOLNBs1   => plt_biom%CPOLNBs1    , &
    ZPOLNBs1   => plt_biom%ZPOLNBs1    , &
    PPOLNBs1   => plt_biom%PPOLNBs1    , &
    WGLFLs1    => plt_biom%WGLFLs1     , &
    WGLFLNs1   => plt_biom%WGLFLNs1    , &
    WGLFLPs1   => plt_biom%WGLFLPs1    , &
    WGLFPs1    => plt_biom%WGLFPs1     , &
    WTSHEBs1   => plt_biom%WTSHEBs1    , &
    CPOOLs1    => plt_biom%CPOOLs1     , &
    ZPOOLs1    => plt_biom%ZPOOLs1     , &
    PPOOLs1    => plt_biom%PPOOLs1     , &
    WGSHEs1    => plt_biom%WGSHEs1     , &
    WTLFBNs1   => plt_biom%WTLFBNs1    , &
    WGLFNs1    => plt_biom%WGLFNs1     , &
    WTLFBPs1   => plt_biom%WTLFBPs1    , &
    WSLFs1     => plt_biom%WSLFs1      , &
    WGLFs1     => plt_biom%WGLFs1      , &
    WTLFBs1    => plt_biom%WTLFBs1     , &
    WSSHEs1    => plt_biom%WSSHEs1     , &
    WTSHBNs1   => plt_biom%WTSHBNs1    , &
    WTSHBPs1   => plt_biom%WTSHBPs1    , &
    WGSHNs1    => plt_biom%WGSHNs1     , &
    WGSHPs1    => plt_biom%WGSHPs1     , &
    WTSTKs1    => plt_biom%WTSTKs1     , &
    CCPOLPs1   => plt_biom%CCPOLPs1    , &
    CCPLNPs1   => plt_biom%CCPLNPs1    , &
    WTLFs1     => plt_biom%WTLFs1      , &
    WTGRs1     => plt_biom%WTGRs1      , &
    WTSHTs1    => plt_biom%WTSHTs1     , &
    WTHSKs1    => plt_biom%WTHSKs1     , &
    WTEARs1    => plt_biom%WTEARs1     , &
    WTSHEs1    => plt_biom%WTSHEs1     , &
    WTSHTAs1   => plt_biom%WTSHTAs1    , &
    WTRT1s1    => plt_biom%WTRT1s1     , &
    WTRT1Ns1   => plt_biom%WTRT1Ns1    , &
    WTRT1Ps1   => plt_biom%WTRT1Ps1    , &
    RTWT1s1    => plt_biom%RTWT1s1     , &
    RTWT1Ns1   => plt_biom%RTWT1Ns1    , &
    RTWT1Ps1   => plt_biom%RTWT1Ps1    , &
    WTRT2s1    => plt_biom%WTRT2s1     , &
    WTRT2Ns1   => plt_biom%WTRT2Ns1    , &
    WTRT2Ps1   => plt_biom%WTRT2Ps1    , &
    WTNDLNs1   => plt_biom%WTNDLNs1    , &
    WTNDLPs1   => plt_biom%WTNDLPs1    , &
    WTNDLs1    => plt_biom%WTNDLs1     , &
    CPOOLNs1   => plt_biom%CPOOLNs1    , &
    ZPOOLNs1   => plt_biom%ZPOOLNs1    , &
    PPOOLNs1   => plt_biom%PPOOLNs1    , &
    ZEROPs1    => plt_biom%ZEROPs1     , &
    ZEROLs1    => plt_biom%ZEROLs1     , &
    WGLFVs1    => plt_biom%WGLFVs1     , &
    FWOODNs1   => plt_allom%FWOODNs1   , &
    FVRNs1     => plt_allom%FVRNs1     , &
    FWODRs1    => plt_allom%FWODRs1    , &
    FWODRPs1   => plt_allom%FWODRPs1   , &
    FWODRNs1   => plt_allom%FWODRNs1   , &
    FWOODPs1   => plt_allom%FWOODPs1   , &
    FWODSNs1   => plt_allom%FWODSNs1   , &
    FWODSPs1   => plt_allom%FWODSPs1   , &
    FWOODs1    => plt_allom%FWOODs1    , &
    FWODBs1    => plt_allom%FWODBs1    , &
    FWODLNs1   => plt_allom%FWODLNs1   , &
    FWODLPs1   => plt_allom%FWODLPs1   , &
    GRWTBs1    => plt_allom%GRWTBs1    , &
    IDTHBs1    =>  plt_pheno%IDTHBs1   , &
    TFN3s1     =>  plt_pheno%TFN3s1    , &
    IDAYs1     =>  plt_pheno%IDAYs1    , &
    GROUPs1    =>  plt_pheno%GROUPs1   , &
    VSTGXs1    =>  plt_pheno%VSTGXs1   , &
    IGTYPs1    =>  plt_pheno%IGTYPs1   , &
    IBTYPs1    =>  plt_pheno%IBTYPs1   , &
    IFLGAs1    =>  plt_pheno%IFLGAs1   , &
    ISTYPs1    =>  plt_pheno%ISTYPs1   , &
    VRNXs1     =>  plt_pheno%VRNXs1    , &
    VRNFs1     =>  plt_pheno%VRNFs1    , &
    IWTYPs1    =>  plt_pheno%IWTYPs1   , &
    TGSTGIs1   =>  plt_pheno%TGSTGIs1  , &
    TGSTGFs1   =>  plt_pheno%TGSTGFs1  , &
    FLG4s1     =>  plt_pheno%FLG4s1    , &
    GROUPIs1   =>  plt_pheno%GROUPIs1  , &
    CORGCs1    =>  plt_soilchem%CORGCs1, &
    THETWs1    =>  plt_soilchem%THETWs1, &
    CFOPCs1    =>  plt_soilchem%CFOPCs1, &
    CFOPNs1    =>  plt_soilchem%CFOPNs1, &
    CFOPPs1    =>  plt_soilchem%CFOPPs1, &
    H2GPs1     =>  plt_rbgc%H2GPs1     , &
    CO2Ps1     =>  plt_rbgc%CO2Ps1     , &
    CH4Ps1     =>  plt_rbgc%CH4Ps1     , &
    OXYPs1     =>  plt_rbgc%OXYPs1     , &
    H2GAs1     =>  plt_rbgc%H2GAs1     , &
    CH4As1     =>  plt_rbgc%CH4As1     , &
    OXYAs1     =>  plt_rbgc%OXYAs1     , &
    CO2As1     =>  plt_rbgc%CO2As1     , &
    CSNCs1     =>  plt_bgcr%CSNCs1     , &
    CNETs1     =>  plt_bgcr%CNETs1     , &
    ZSNCs1     =>  plt_bgcr%ZSNCs1     , &
    PSNCs1     =>  plt_bgcr%PSNCs1     , &
    TNBPs1     =>  plt_bgcr%TNBPs1     , &
    RCH4Zs1    =>  plt_bgcr%RCH4Zs1    , &
    RCO2Zs1    =>  plt_bgcr%RCO2Zs1    , &
    ROXYZs1    =>  plt_bgcr%ROXYZs1    , &
    RN2OZs1    =>  plt_bgcr%RN2OZs1    , &
    RNH3Zs1    =>  plt_bgcr%RNH3Zs1    , &
    RH2GZs1    =>  plt_bgcr%RH2GZs1    , &
    ZH3Ps1     =>  plt_rbgc%ZH3Ps1     , &
    Z2OPs1     =>  plt_rbgc%Z2OPs1     , &
    Z2OAs1     =>  plt_rbgc%Z2OAs1     , &
    RCO2As1    =>  plt_rbgc%RCO2As1    , &
    RCO2Ms1    =>  plt_rbgc%RCO2Ms1    , &
    RCO2Ns1    =>  plt_rbgc%RCO2Ns1    , &
    ZH3As1     =>  plt_rbgc%ZH3As1     , &
    RTNLs1     =>  plt_morph%RTNLs1    , &
    RTN2s1     =>  plt_morph%RTN2s1    , &
    RTLGPs1    =>  plt_morph%RTLGPs1   , &
    RTARPs1    =>  plt_morph%RTARPs1   , &
    RTVLPs1    =>  plt_morph%RTVLPs1   , &
    NGs1       =>  plt_morph%NGs1      , &
    MYs1       =>  plt_morph%MYs1      , &
    ZCs1       =>  plt_morph%ZCs1      , &
    RTVLWs1    =>  plt_morph%RTVLWs1   , &
    RTDNPs1    =>  plt_morph%RTDNPs1   , &
    INTYPs1    =>  plt_morph%INTYPs1   , &
    ARLFTs1    =>  plt_morph%ARLFTs1   , &
    ZLs1       =>  plt_morph%ZLs1      , &
    ARLFBs1    =>  plt_morph%ARLFBs1   , &
    NBRs1      =>  plt_morph%NBRs1     , &
    ARSTPs1    =>  plt_morph%ARSTPs1   , &
    HTNODXs1   =>  plt_morph%HTNODXs1  , &
    RTN1s1     =>  plt_morph%RTN1s1    , &
    RTLG2s1    =>  plt_morph%RTLG2s1   , &
    RTLG1s1    =>  plt_morph%RTLG1s1   , &
    HTNODEs1   =>  plt_morph%HTNODEs1  , &
    GRNXBs1    =>  plt_morph%GRNXBs1   , &
    GRNOBs1    =>  plt_morph%GRNOBs1   , &
    HTSHEs1    =>  plt_morph%HTSHEs1   , &
    ARLF1s1    =>  plt_morph%ARLF1s1   , &
    ARLFVs1    =>  plt_morph%ARLFVs1   , &
    ARSTVs1    =>  plt_morph%ARSTVs1   , &
    ARLFLs1    =>  plt_morph%ARLFLs1   , &
    ARSTKs1    =>  plt_morph%ARSTKs1   , &
    NRTs1      =>  plt_morph%NRTs1     , &
    PSTGFs1    =>  plt_morph%PSTGFs1   , &
    NB1s1      =>  plt_morph%NB1s1     , &
    PSTGIs1    =>  plt_morph%PSTGIs1   , &
    PSTGs1     =>  plt_morph%PSTGs1    , &
    CFs1       =>  plt_morph%CFs1      , &
    ARLFCs1    =>  plt_morph%ARLFCs1   , &
    ICTYPs1    =>  plt_photo%ICTYPs1     &
  )
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
  IF((IHVSTs1(NZ).GE.0.AND.J.EQ.INT(ZNOONs1) &
    .AND.IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6) &
    .OR.(IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6))THEN
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
    IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
      IF(JHVSTs1(NZ).NE.2)THEN
        PPXs1(NZ)=PPXs1(NZ)*(1._r8-THINs1(NZ))
        PPs1(NZ)=PPs1(NZ)*(1._r8-THINs1(NZ))
      ELSE
!     PPIs1(NZ)=AMAX1(1.0,0.5*(PPIs1(NZ)+GRNOs1(NZ)/AREA3s1(NUs1)))
        PPXs1(NZ)=PPIs1(NZ)
        PPs1(NZ)=PPXs1(NZ)*AREA3s1(NUs1)
      ENDIF
      IF(IHVSTs1(NZ).EQ.3)THEN
        CFs1(NZ)=CFs1(NZ)*HVSTs1(NZ)
      ENDIF
      IF(IHVSTs1(NZ).LE.2.AND.HVSTs1(NZ).LT.0.0)THEN
        ARLFY=(1._r8-ABS(HVSTs1(NZ)))*ARLFCs1
        ARLFR=0._r8
        DO 9875 L=1,JC1
          IF(ZLs1(L).GT.ZLs1(L-1).AND.ARLFTs1(L).GT.ZEROSs1 &
            .AND.ARLFR.LT.ARLFY)THEN
            IF(ARLFR+ARLFTs1(L).GT.ARLFY)THEN
              HVSTs1(NZ)=ZLs1(L-1)+((ARLFY-ARLFR) &
                /ARLFTs1(L))*(ZLs1(L)-ZLs1(L-1))
            ENDIF
          ELSE
            HVSTs1(NZ)=0._r8
          ENDIF
          ARLFR=ARLFR+ARLFTs1(L)
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
!     WTSHTAs1=average biomass in landscape grazing section
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
      IF(WTSHTAs1(NZ).GT.ZEROPs1(NZ))THEN
        WHVSTT=HVSTs1(NZ)*THINs1(NZ)*0.45/24.0 &
          *AREA3s1(NUs1)*WTSHTs1(NZ)/WTSHTAs1(NZ)
      ELSE
        WHVSTT=0._r8
      ENDIF
      IF(IHVSTs1(NZ).EQ.6)THEN
        WHVSTT=WHVSTT*TFN3s1(NZ)
      ENDIF
      CCPOLX=CCPOLPs1(NZ)/(1.0+CCPOLPs1(NZ))
      CCPLNX=CCPLNPs1(NZ)/(1.0+CCPLNPs1(NZ))
!
!     LEAF,BACTERIA GRAZED,REMOVED
!
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVSTs1(2,1,EHVSTs1(2,2,EHVSTs1(2,3,EHVSTs1(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosyst
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WTLF=PFT leaf C mass
!     WHVXXX=grazing requirement unmet by leaf
!
      WHVSLX=WHVSTT*EHVSTs1(1,1,NZ)
      WHVSLY=AMIN1(WTLFs1(NZ),WHVSLX)
      WHVSLF=WHVSLY*(1._r8-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVSLX-WHVSLY)
      WHVSSX=WHVSTT*EHVSTs1(1,2,NZ)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     WHVXXX=grazing requirement unmet by non-foliar removal
!
      WTSHTT=WTSHEs1(NZ)+WTHSKs1(NZ)+WTEARs1(NZ)+WTGRs1(NZ)
      IF(WTSHTT.GT.ZEROPs1(NZ))THEN
        WHVSHX=WHVSSX*WTSHEs1(NZ)/WTSHTT+WHVXXX
        WHVSHY=AMIN1(WTSHEs1(NZ),WHVSHX)
        WHVSHH=WHVSHY*(1._r8-CCPOLX)
        WHVSCS=WHVSHY*CCPOLX
        WHVSNS=WHVSHY*CCPLNX
        WHVXXX=AMAX1(0.0,WHVSHX-WHVSHY)
        WHVHSX=WHVSSX*WTHSKs1(NZ)/WTSHTT+WHVXXX
        WHVHSY=AMIN1(WTHSKs1(NZ),WHVHSX)
        WHVHSH=WHVHSY
        WHVXXX=AMAX1(0.0,WHVHSX-WHVHSY)
        WHVEAX=WHVSSX*WTEARs1(NZ)/WTSHTT+WHVXXX
        WHVEAY=AMIN1(WTEARs1(NZ),WHVEAX)
        WHVEAH=WHVEAY
        WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
        WHVGRX=WHVSSX*WTGRs1(NZ)/WTSHTT+WHVXXX
        WHVGRY=AMIN1(WTGRs1(NZ),WHVGRX)
        WHVGRH=WHVGRY
        WHVXXX=AMAX1(0.0,WHVGRX-WHVGRY)
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
      WHVSKX=WHVSTT*EHVSTs1(1,3,NZ)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     WHVXXX=grazing requirement unmet by stalk,reserve
!
      WTSTKT=WTSTKs1(NZ)+WTRSVs1(NZ)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
        WHVSTX=WHVSKX*WTSTKs1(NZ)/WTSTKT+WHVXXX
        WHVSTY=AMIN1(WTSTKs1(NZ),WHVSTX)
        WHVSTH=WHVSTY
        WHVXXX=AMAX1(0.0,WHVSTX-WHVSTY)
        WHVRVX=WHVSKX*WTRSVs1(NZ)/WTSTKT+WHVXXX
        WHVRVY=AMIN1(WTRSVs1(NZ),WHVRVX)
        WHVRVH=WHVRVY
        WHVXXX=AMAX1(0.0,WHVRVX-WHVRVY)
      ELSE
        WHVSTH=0._r8
        WHVRVH=0._r8
        WHVXXX=AMAX1(0.0,WHVSKX)
!
!     ALLOCATE UNMET DEMAND FOR GRAZING TO LEAF,PETIOLE,HUSK
!     EAR,GRAIN
!
!     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
!     WHVSH*,WHVHS,WHVEA,WHVGR,WHVSC=
!            petiole,husk,ear,grain,nonstructural C removed
!
        IF(WHVXXX.GT.0.0)THEN
          WHVSLY=AMIN1(WTLFs1(NZ)-WHVSLF-WHVSCL,WHVXXX)
          WHVSLF=WHVSLF+WHVSLY*(1._r8-CCPOLX)
          WHVSCL=WHVSCL+WHVSLY*CCPOLX
          WHVSNL=WHVSNL+WHVSLY*CCPLNX
          WHVXXX=AMAX1(0.0,WHVXXX-WHVSLY)
          IF(WTSHTT.GT.ZEROPs1(NZ))THEN
            WHVSHX=WHVXXX*WTSHEs1(NZ)/WTSHTT
            WHVSHY=AMIN1(WTSHEs1(NZ),WHVSHX)
            WHVSHH=WHVSHH+WHVSHY*(1._r8-CCPOLX)
            WHVSCS=WHVSCS+WHVSHY*CCPOLX
            WHVSNS=WHVSNS+WHVSHY*CCPLNX
            WHVXXX=AMAX1(0.0,WHVXXX-WHVSHY)
            WHVHSX=WHVXXX*WTHSKs1(NZ)/WTSHTT
            WHVHSY=AMIN1(WTHSKs1(NZ),WHVHSX)
            WHVHSH=WHVHSH+WHVHSY
            WHVXXX=AMAX1(0.0,WHVXXX-WHVHSY)
            WHVEAX=WHVXXX*WTEARs1(NZ)/WTSHTT
            WHVEAY=AMIN1(WTEARs1(NZ),WHVEAX)
            WHVEAH=WHVEAH+WHVEAY
            WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
            WHVGRX=WHVXXX*WTGRs1(NZ)/WTSHTT
            WHVGRY=AMIN1(WTGRs1(NZ),WHVGRX)
            WHVGRH=WHVGRH+WHVGRY
            WHVXXX=AMAX1(0.0,WHVGRX-WHVGRY)
          ENDIF
        ENDIF
      ENDIF
!
!     ALL HARVEST REMOVALS
!
!     WGLFBL=branch leaf C mass in canopy layer
!
      DO 9860 NB=1,NBRs1(NZ)
        DO  L=1,JC1
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=0._r8
          enddo
        enddo
9860  CONTINUE
      DO 9870 NB=1,NBRs1(NZ)
        DO  L=1,JC1
          DO  K=0,JNODS1
            WGLFBL(L,NB,NZ)=WGLFBL(L,NB,NZ)+WGLFLs1(L,K,NB,NZ)
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
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!
    DO 9865 L=JC1,1,-1
      IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
        IF(IHVSTs1(NZ).NE.3)THEN
          IF(ZLs1(L).GT.ZLs1(L-1))THEN
            FHGT=AMAX1(0.0,AMIN1(1.0,1._r8-((ZLs1(L)) &
              -HVSTs1(NZ))/(ZLs1(L)-ZLs1(L-1))))
          ELSE
            FHGT=1.0_r8
          ENDIF
        ELSE
          FHGT=0._r8
        ENDIF
        IF(test_aeqb(THINs1(NZ),0._r8))THEN
          FHVST=AMAX1(0.0,1._r8-(1._r8-FHGT)*EHVSTs1(1,1,NZ))
          FHVSH=FHVST
        ELSE
          FHVST=AMAX1(0.0,1._r8-THINs1(NZ))
          IF(IHVSTs1(NZ).EQ.0)THEN
            FHVSH=1.0_r8-(1._r8-FHGT)*EHVSTs1(1,1,NZ)*THINs1(NZ)
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
      DO 9855 NB=1,NBRs1(NZ)
        IF((IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6) &
          .AND.WTLFs1(NZ).GT.ZEROLs1(NZ))THEN
          WHVSBL=WHVSLF*AMAX1(0.0,WGLFBL(L,NB,NZ))/WTLFs1(NZ)
        ELSE
          WHVSBL=0._r8
        ENDIF
        DO 9845 K=JNODS1,0,-1
          IF((IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6).OR.WHVSBL.GT.0.0)THEN
            IF(IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6)THEN
              IF(WGLFLs1(L,K,NB,NZ).GT.WHVSBL)THEN
                FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFLs1(L,K,NB,NZ)-WHVSBL)/WGLFLs1(L,K,NB,NZ)))
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
            WHVSBL=WHVSBL-(1._r8-FHVST)*WGLFLs1(L,K,NB,NZ)
            WTHTH1=WTHTH1+(1._r8-FHVSH)*WGLFLs1(L,K,NB,NZ)*FWODBs1(1)
            WTHNH1=WTHNH1+(1._r8-FHVSH)*WGLFLNs1(L,K,NB,NZ)*FWODLNs1(1)
            WTHPH1=WTHPH1+(1._r8-FHVSH)*WGLFLPs1(L,K,NB,NZ)*FWODLPs1(1)
            WTHTX1=WTHTX1+(FHVSH-FHVST)*WGLFLs1(L,K,NB,NZ)*FWODBs1(1)
            WTHNX1=WTHNX1+(FHVSH-FHVST)*WGLFLNs1(L,K,NB,NZ)*FWODLNs1(1)
            WTHPX1=WTHPX1+(FHVSH-FHVST)*WGLFLPs1(L,K,NB,NZ)*FWODLPs1(1)
            WTHTH3=WTHTH3+(1._r8-FHVSH)*WGLFLs1(L,K,NB,NZ)*FWODBs1(0)
            WTHNH3=WTHNH3+(1._r8-FHVSH)*WGLFLNs1(L,K,NB,NZ)*FWODLNs1(0)
            WTHPH3=WTHPH3+(1._r8-FHVSH)*WGLFLPs1(L,K,NB,NZ)*FWODLPs1(0)
            WTHTX3=WTHTX3+(FHVSH-FHVST)*WGLFLs1(L,K,NB,NZ)*FWODBs1(0)
            WTHNX3=WTHNX3+(FHVSH-FHVST)*WGLFLNs1(L,K,NB,NZ)*FWODLNs1(0)
            WTHPX3=WTHPX3+(FHVSH-FHVST)*WGLFLPs1(L,K,NB,NZ)*FWODLPs1(0)
!
!     REMAINING LEAF C,N,P AND AREA
!
            WGLFLs1(L,K,NB,NZ)=FHVST*WGLFLs1(L,K,NB,NZ)
            WGLFLNs1(L,K,NB,NZ)=FHVST*WGLFLNs1(L,K,NB,NZ)
            WGLFLPs1(L,K,NB,NZ)=FHVST*WGLFLPs1(L,K,NB,NZ)
            ARLFLs1(L,K,NB,NZ)=FHVST*ARLFLs1(L,K,NB,NZ)
            IF(K.EQ.1)THEN
              ARSTKs1(L,NB,NZ)=FHVST*ARSTKs1(L,NB,NZ)
            ENDIF
          ENDIF

9845      CONTINUE
9855    CONTINUE
        ARLFVs1(L,NZ)=0._r8
        WGLFVs1(L,NZ)=0._r8
        ARSTVs1(L,NZ)=ARSTVs1(L,NZ)*FHVST
9865  CONTINUE
      DO 9835 NB=1,NBRs1(NZ)
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
            ARLFG=ARLFG+ARLFLs1(L,K,NB,NZ)
            WGLFG=WGLFG+WGLFLs1(L,K,NB,NZ)
            WGLFNG=WGLFNG+WGLFLNs1(L,K,NB,NZ)
            WGLFPG=WGLFPG+WGLFLPs1(L,K,NB,NZ)
            ARLFVs1(L,NZ)=ARLFVs1(L,NZ)+ARLFLs1(L,K,NB,NZ)
            WGLFVs1(L,NZ)=WGLFVs1(L,NZ)+WGLFLs1(L,K,NB,NZ)
9815      CONTINUE
!
!     CUT STALK AT HARVESTED NODES AND LAYERS
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     WGLF=leaf node C mass
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FHVSTK=fraction of internode layer mass not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!
          IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
            IF(WGLFs1(K,NB,NZ).GT.ZEROPs1(NZ) &
              .AND.EHVSTs1(1,1,NZ).GT.0.0)THEN
              FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(1._r8-(1._r8-AMAX1(0.0,WGLFG) &
                /WGLFs1(K,NB,NZ))*EHVSTs1(1,2,NZ)/EHVSTs1(1,1,NZ))))
              FHVSHK(K)=FHVSTK(K)
          ELSE
            IF(test_aeqb(THINs1(NZ),0._r8))THEN
              FHVSTK(K)=1.0_r8-EHVSTs1(1,2,NZ)
              FHVSHK(K)=FHVSTK(K)
            ELSE
              FHVSTK(K)=1.0_r8-THINs1(NZ)
              IF(IHVSTs1(NZ).EQ.0)THEN
                FHVSHK(K)=1.0_r8-EHVSTs1(1,2,NZ)*THINs1(NZ)
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
        WGLFGY=WGLFGY+WGLFs1(K,NB,NZ)
        WTLFBs1(NB,NZ)=WTLFBs1(NB,NZ)-WGLFs1(K,NB,NZ)+WGLFG
        WTLFBNs1(NB,NZ)=WTLFBNs1(NB,NZ)-WGLFNs1(K,NB,NZ)+WGLFNG
        WTLFBPs1(NB,NZ)=WTLFBPs1(NB,NZ)-WGLFPs1(K,NB,NZ)+WGLFPG
        ARLFBs1(NB,NZ)=ARLFBs1(NB,NZ)-ARLF1s1(K,NB,NZ)+ARLFG
        IF(ARLF1s1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
          WSLFs1(K,NB,NZ)=WSLFs1(K,NB,NZ)*ARLFG/ARLF1s1(K,NB,NZ)
        ELSE
          WSLFs1(K,NB,NZ)=0._r8
        ENDIF
        ARLF1s1(K,NB,NZ)=ARLFG
        WGLFs1(K,NB,NZ)=WGLFG
        WGLFNs1(K,NB,NZ)=WGLFNG
        WGLFPs1(K,NB,NZ)=WGLFPG
        WGLFGX=WGLFGX+WGLFs1(K,NB,NZ)
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
      IF((IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6) &
        .AND.WTSHEs1(NZ).GT.ZEROPs1(NZ))THEN
        WHVSBS=WHVSHH*WTSHEBs1(NB,NZ)/WTSHEs1(NZ)
      ELSE
        WHVSBS=0._r8
      ENDIF
      DO 9805 K=JNODS1,0,-1
!112   FORMAT(A8,8I4,12E12.4)
        IF(HTNODEs1(K,NB,NZ).GT.0.0) &
          HTSTKX=AMAX1(HTSTKX,HTNODEs1(K,NB,NZ))
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
          IF((IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6).OR.WHVSBS.GT.0.0)THEN
            IF(IHVSTs1(NZ).EQ.4.OR.IHVSTs1(NZ).EQ.6)THEN
              IF(WGSHEs1(K,NB,NZ).GT.WHVSBS)THEN
                FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(WGSHEs1(K,NB,NZ)-WHVSBS)/WGSHEs1(K,NB,NZ)))
                FHVSHK(K)=FHVSTK(K)
              ELSE
                FHVSTK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1._r8-FHVSTK(K))*WGSHEs1(K,NB,NZ)
            WTHTH2=WTHTH2+(1._r8-FHVSHK(K))*WGSHEs1(K,NB,NZ)*FWODBs1(1)
            WTHNH2=WTHNH2+(1._r8-FHVSHK(K))*WGSHNs1(K,NB,NZ)*FWODSNs1(1)
            WTHPH2=WTHPH2+(1._r8-FHVSHK(K))*WGSHPs1(K,NB,NZ)*FWODSPs1(1)
            WTHTX2=WTHTX2+(FHVSHK(K)-FHVSTK(K))*WGSHEs1(K,NB,NZ)*FWODBs1(1)
            WTHNX2=WTHNX2+(FHVSHK(K)-FHVSTK(K))*WGSHNs1(K,NB,NZ)*FWODSNs1(1)
            WTHPX2=WTHPX2+(FHVSHK(K)-FHVSTK(K))*WGSHPs1(K,NB,NZ)*FWODSPs1(1)
            WTHTH3=WTHTH3+(1._r8-FHVSHK(K))*WGSHEs1(K,NB,NZ)*FWODBs1(0)
            WTHNH3=WTHNH3+(1._r8-FHVSHK(K))*WGSHNs1(K,NB,NZ)*FWODSNs1(0)
            WTHPH3=WTHPH3+(1._r8-FHVSHK(K))*WGSHPs1(K,NB,NZ)*FWODSPs1(0)
            WTHTX3=WTHTX3+(FHVSHK(K)-FHVSTK(K))*WGSHEs1(K,NB,NZ)*FWODBs1(0)
            WTHNX3=WTHNX3+(FHVSHK(K)-FHVSTK(K))*WGSHNs1(K,NB,NZ)*FWODSNs1(0)
            WTHPX3=WTHPX3+(FHVSHK(K)-FHVSTK(K))*WGSHPs1(K,NB,NZ)*FWODSPs1(0)
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     WGSHE=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     HTSHE=node petiole height
!     WSSHE=petiole protein mass
!
            WGSHGY=WGSHGY+WGSHEs1(K,NB,NZ)
            WTSHEBs1(NB,NZ)=WTSHEBs1(NB,NZ) &
              -(1._r8-FHVSTK(K))*WGSHEs1(K,NB,NZ)
            WTSHBNs1(NB,NZ)=WTSHBNs1(NB,NZ) &
              -(1._r8-FHVSTK(K))*WGSHNs1(K,NB,NZ)
            WTSHBPs1(NB,NZ)=WTSHBPs1(NB,NZ) &
              -(1._r8-FHVSTK(K))*WGSHPs1(K,NB,NZ)
            WGSHEs1(K,NB,NZ)=FHVSTK(K)*WGSHEs1(K,NB,NZ)
            WSSHEs1(K,NB,NZ)=FHVSTK(K)*WSSHEs1(K,NB,NZ)
            WGSHNs1(K,NB,NZ)=FHVSTK(K)*WGSHNs1(K,NB,NZ)
            WGSHPs1(K,NB,NZ)=FHVSTK(K)*WGSHPs1(K,NB,NZ)
            WSSHEs1(K,NB,NZ)=FHVSTK(K)*WSSHEs1(K,NB,NZ)
            IF(IHVSTs1(NZ).LE.2 &
              .AND.HTSHEs1(K,NB,NZ).GT.0.0)THEN
              FHGT=AMAX1(0.0,AMIN1(1.0,(HTNODEs1(K,NB,NZ) &
                +HTSHEs1(K,NB,NZ)-HVSTs1(NZ))/HTSHEs1(K,NB,NZ)))
              HTSHEs1(K,NB,NZ)=(1._r8-FHGT)*HTSHEs1(K,NB,NZ)
            ELSE
              HTSHEs1(K,NB,NZ)=FHVSTK(K)*HTSHEs1(K,NB,NZ)
            ENDIF
            WGSHGX=WGSHGX+WGSHEs1(K,NB,NZ)
!     IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
!     IF(HTNODEs1(K,NB,NZ).GT.HVSTs1(NZ)
!    2.OR.IHVSTs1(NZ).EQ.3)THEN
!     IF(test_aeqb(FHVSTK(K),0._r8).AND.K.GT.0)THEN
!     IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
!     VSTGs1(NB,NZ)=AMAX1(0.0,VSTGs1(NB,NZ)-1.0)
!     ELSE
!     VSTGs1(NB,NZ)=AMAX1(0.0,VSTGs1(NB,NZ)-0.04)
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
        CPOOLX=AMAX1(0.0,CPOOLs1(NB,NZ))
        ZPOOLX=AMAX1(0.0,ZPOOLs1(NB,NZ))
        PPOOLX=AMAX1(0.0,PPOOLs1(NB,NZ))
        CPOLNX=AMAX1(0.0,CPOLNBs1(NB,NZ))
        ZPOLNX=AMAX1(0.0,ZPOLNBs1(NB,NZ))
        PPOLNX=AMAX1(0.0,PPOLNBs1(NB,NZ))
        IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
          IF(WGLFGY+WGSHGY.GT.ZEROPs1(NZ))THEN
            FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
            CPOOLG=CPOOLX*FHVST
            ZPOOLG=ZPOOLX*FHVST
            PPOOLG=PPOOLX*FHVST
            CPOLNG=CPOLNX*FHVST
            ZPOLNG=ZPOLNX*FHVST
            PPOLNG=PPOLNX*FHVST
            WTNDG=WTNDBs1(NB,NZ)*FHVST
            WTNDNG=WTNDBNs1(NB,NZ)*FHVST
            WTNDPG=WTNDBPs1(NB,NZ)*FHVST
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
          IF(WTLSs1(NZ).GT.ZEROLs1(NZ))THEN
            WTLSBX=AMAX1(0.0,WTLSBs1(NB,NZ))
            IF(CPOOLs1(NB,NZ).GT.ZEROPs1(NZ))THEN
              WHVSCX=AMAX1(0.0,WHVSCP)*WTLSBX/WTLSs1(NZ)
              CPOOLG=AMAX1(0.0,CPOOLX-WHVSCX)
              ZPOOLG=AMAX1(0.0,ZPOOLX-WHVSCX*ZPOOLX/CPOOLs1(NB,NZ))
              PPOOLG=AMAX1(0.0,PPOOLX-WHVSCX*PPOOLX/CPOOLs1(NB,NZ))
            ELSE
              CPOOLG=0._r8
              ZPOOLG=0._r8
              PPOOLG=0._r8
            ENDIF
            IF(CPOLNBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
              WHVSNX=AMAX1(0.0,WHVSNP)*WTLSBX/WTLSs1(NZ)
              CPOLNG=AMAX1(0.0,CPOLNX-WHVSNX)
              ZPOLNG=AMAX1(0.0,ZPOLNX-WHVSNX*ZPOLNX/CPOLNBs1(NB,NZ))
              PPOLNG=AMAX1(0.0,PPOLNX-WHVSNX*PPOLNX/CPOLNBs1(NB,NZ))
              WTNDG=WTNDBs1(NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDNG=WTNDBNs1(NB,NZ)*(1._r8-WHVSNX/CPOLNX)
              WTNDPG=WTNDBPs1(NB,NZ)*(1._r8-WHVSNX/CPOLNX)
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
        WTHTH0=WTHTH0+WTNDBs1(NB,NZ)-WTNDG
        WTHNH0=WTHNH0+WTNDBNs1(NB,NZ)-WTNDNG
        WTHPH0=WTHPH0+WTNDBPs1(NB,NZ)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
        CPOOLs1(NB,NZ)=CPOOLG
        ZPOOLs1(NB,NZ)=ZPOOLG
        PPOOLs1(NB,NZ)=PPOOLG
        CPOLNBs1(NB,NZ)=CPOLNG
        ZPOLNBs1(NB,NZ)=ZPOLNG
        PPOLNBs1(NB,NZ)=PPOLNG
        WTNDBs1(NB,NZ)=WTNDG
        WTNDBNs1(NB,NZ)=WTNDNG
        WTNDBPs1(NB,NZ)=WTNDPG
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
        IF(ICTYPs1(NZ).EQ.4.AND.CPOOLX.GT.ZEROPs1(NZ))THEN
          FHVST4=CPOOLG/CPOOLX
          DO 9810 K=1,JNODS1
            WTHTH0=WTHTH0+(1._r8-FHVST4)*CPOOL3s1(K,NB,NZ)
            WTHTH0=WTHTH0+(1._r8-FHVST4)*CPOOL4s1(K,NB,NZ)
            WTHTH0=WTHTH0+(1._r8-FHVST4)*CO2Bs1(K,NB,NZ)
            WTHTH0=WTHTH0+(1._r8-FHVST4)*HCOBs1(K,NB,NZ)
            CPOOL3s1(K,NB,NZ)=FHVST4*CPOOL3s1(K,NB,NZ)
            CPOOL4s1(K,NB,NZ)=FHVST4*CPOOL4s1(K,NB,NZ)
            CO2Bs1(K,NB,NZ)=FHVST4*CO2Bs1(K,NB,NZ)
            HCOBs1(K,NB,NZ)=FHVST4*HCOBs1(K,NB,NZ)
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
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     WTSTK=stalk C mass
!
!
        IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
          IF(HTSTKX.GT.ZEROs1)THEN
            IF(IHVSTs1(NZ).NE.3)THEN
              FHGT=AMAX1(0.0,AMIN1(1.0,HVSTs1(NZ)/HTSTKX))
            ELSE
              FHGT=0._r8
            ENDIF
            IF(test_aeqb(THINs1(NZ),0._r8))THEN
              FHVST=AMAX1(0.0,1._r8-(1._r8-FHGT)*EHVSTs1(1,3,NZ))
              FHVSH=FHVST
            ELSE
              FHVST=AMAX1(0.0,1._r8-THINs1(NZ))
              IF(IHVSTs1(NZ).EQ.0)THEN
                FHVSH=1.0_r8-(1._r8-FHGT)*EHVSTs1(1,3,NZ)*THINs1(NZ)
              ELSE
                FHVSH=FHVST
              ENDIF
            ENDIF
          ELSE
            FHVST=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ELSE
          IF(WTSTKs1(NZ).GT.ZEROLs1(NZ))THEN
            FHVST=AMAX1(0.0,AMIN1(1.0,1._r8-WHVSTH/WTSTKs1(NZ)))
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
        WTHTH3=WTHTH3+(1._r8-FHVSH)*WTSTKBs1(NB,NZ)
        WTHNH3=WTHNH3+(1._r8-FHVSH)*WTSTBNs1(NB,NZ)
        WTHPH3=WTHPH3+(1._r8-FHVSH)*WTSTBPs1(NB,NZ)
        WTHTX3=WTHTX3+(FHVSH-FHVST)*WTSTKBs1(NB,NZ)
        WTHNX3=WTHNX3+(FHVSH-FHVST)*WTSTBNs1(NB,NZ)
        WTHPX3=WTHPX3+(FHVSH-FHVST)*WTSTBPs1(NB,NZ)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
        WTSTKBs1(NB,NZ)=FHVST*WTSTKBs1(NB,NZ)
        WTSTBNs1(NB,NZ)=FHVST*WTSTBNs1(NB,NZ)
        WTSTBPs1(NB,NZ)=FHVST*WTSTBPs1(NB,NZ)
        WVSTKBs1(NB,NZ)=FHVST*WVSTKBs1(NB,NZ)
        WTSTXBs1(NB,NZ)=FHVST*WTSTXBs1(NB,NZ)
        WTSTXNs1(NB,NZ)=FHVST*WTSTXNs1(NB,NZ)
        WTSTXPs1(NB,NZ)=FHVST*WTSTXPs1(NB,NZ)
!
!     CUT STALK NODES
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     HTNODX,HTNODE=stalk height,stalk internode length
!     FHGTK=fraction of internode length not harvested
!     THIN=IHVST=0-3,5: fraction of population removed,
!          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTSTK=stalk C mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
        DO 9820 K=JNODS1,0,-1
          IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
            IF(HTNODXs1(K,NB,NZ).GT.ZEROs1)THEN
              IF(IHVSTs1(NZ).NE.3)THEN
                FHGTK=AMAX1(0.0,AMIN1(1.0,(HTNODEs1(K,NB,NZ) &
                  -HVSTs1(NZ))/HTNODXs1(K,NB,NZ)))
              ELSE
                FHGTK=0._r8
              ENDIF
              IF(test_aeqb(THINs1(NZ),0._r8))THEN
                FHVSTS=AMAX1(0.0,1._r8-FHGTK*EHVSTs1(1,3,NZ))
              ELSE
                FHVSTS=AMAX1(0.0,1._r8-THINs1(NZ))
              ENDIF
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ELSE
            IF(WTSTKs1(NZ).GT.ZEROPs1(NZ))THEN
              FHVSTS=AMAX1(0.0,AMIN1(1.0,1._r8-WHVSTH/WTSTKs1(NZ)))
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ENDIF
          WGNODEs1(K,NB,NZ)=FHVSTS*WGNODEs1(K,NB,NZ)
          WGNODNs1(K,NB,NZ)=FHVSTS*WGNODNs1(K,NB,NZ)
          WGNODPs1(K,NB,NZ)=FHVSTS*WGNODPs1(K,NB,NZ)
          IF(IHVSTs1(NZ).LE.2.AND.test_aeqb(THINs1(NZ),0._r8))THEN
            HTNODXs1(K,NB,NZ)=FHVSTS*HTNODXs1(K,NB,NZ)
            HTNODEs1(K,NB,NZ)=AMIN1(HTNODEs1(K,NB,NZ),HVSTs1(NZ))
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
        IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
          IF(WTSTKBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
            FHVST=FHVST
            FHVSH=FHVSH
          ELSE
            FHVST=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(WTRSVs1(NZ).GT.ZEROPs1(NZ))THEN
            FHVST=AMAX1(0.0,AMIN1(1.0,1._r8-WHVRVH/WTRSVs1(NZ)))
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
        WTHTH3=WTHTH3+(1._r8-FHVSH)*WTRSVBs1(NB,NZ)
        WTHNH3=WTHNH3+(1._r8-FHVSH)*WTRSBNs1(NB,NZ)
        WTHPH3=WTHPH3+(1._r8-FHVSH)*WTRSBPs1(NB,NZ)
        WTHTX3=WTHTX3+(FHVSH-FHVST)*WTRSVBs1(NB,NZ)
        WTHNX3=WTHNX3+(FHVSH-FHVST)*WTRSBNs1(NB,NZ)
        WTHPX3=WTHPX3+(FHVSH-FHVST)*WTRSBPs1(NB,NZ)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
        WTRSVBs1(NB,NZ)=FHVST*WTRSVBs1(NB,NZ)
        WTRSBNs1(NB,NZ)=FHVST*WTRSBNs1(NB,NZ)
        WTRSBPs1(NB,NZ)=FHVST*WTRSBPs1(NB,NZ)
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
!     EHVSTs1(1,1,EHVSTs1(1,2,EHVSTs1(1,3,EHVSTs1(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
!
        IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
          IF(HVSTs1(NZ).LT.HTSTKX &
            .OR.IHVSTs1(NZ).EQ.1 &
            .OR.IHVSTs1(NZ).EQ.3)THEN
            IF(test_aeqb(THINs1(NZ),0._r8))THEN
              FHVSTG=1.0_r8-EHVSTs1(1,2,NZ)
              FHVSHG=FHVSTG
            ELSE
              FHVSTG=1.0_r8-THINs1(NZ)
              FHVSHG=1.0_r8-EHVSTs1(1,2,NZ)*THINs1(NZ)
            ENDIF
          ELSE
            FHVSTG=1.0_r8-THINs1(NZ)
            FHVSHG=FHVSTG
          ENDIF
          FHVSTH=FHVSTG
          FHVSTE=FHVSTG
          FHVSHH=FHVSHG
          FHVSHE=FHVSHG
        ELSE
          IF(WTHSKs1(NZ).GT.ZEROPs1(NZ))THEN
            FHVSTH=AMAX1(0.0,AMIN1(1.0,1._r8-WHVHSH/WTHSKs1(NZ)))
            FHVSHH=FHVSTH
          ELSE
            FHVSTH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(WTEARs1(NZ).GT.ZEROPs1(NZ))THEN
            FHVSTE=AMAX1(0.0,AMIN1(1.0,1._r8-WHVEAH/WTEARs1(NZ)))
            FHVSHE=FHVSTE
          ELSE
            FHVSTE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(WTGRs1(NZ).GT.ZEROPs1(NZ))THEN
            FHVSTG=AMAX1(0.0,AMIN1(1.0,1._r8-WHVGRH/WTGRs1(NZ)))
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
        WTHTH2=WTHTH2+(1._r8-FHVSHH)*WTHSKBs1(NB,NZ)+(1._r8-FHVSHE) &
          *WTEARBs1(NB,NZ)+(1._r8-FHVSHG)*WTGRBs1(NB,NZ)
        WTHNH2=WTHNH2+(1._r8-FHVSHH)*WTHSBNs1(NB,NZ)+(1._r8-FHVSHE) &
          *WTEABNs1(NB,NZ)+(1._r8-FHVSHG)*WTGRBNs1(NB,NZ)
        WTHPH2=WTHPH2+(1._r8-FHVSHH)*WTHSBPs1(NB,NZ)+(1._r8-FHVSHE) &
          *WTEABPs1(NB,NZ)+(1._r8-FHVSHG)*WTGRBPs1(NB,NZ)
        WTHTX2=WTHTX2+(FHVSHH-FHVSTH)*WTHSKBs1(NB,NZ)+(FHVSHE-FHVSTE) &
          *WTEARBs1(NB,NZ)+(FHVSHG-FHVSTG)*WTGRBs1(NB,NZ)
        WTHNX2=WTHNX2+(FHVSHH-FHVSTH)*WTHSBNs1(NB,NZ)+(FHVSHE-FHVSTE) &
          *WTEABNs1(NB,NZ)+(FHVSHG-FHVSTG)*WTGRBNs1(NB,NZ)
        WTHPX2=WTHPX2+(FHVSHH-FHVSTH)*WTHSBPs1(NB,NZ)+(FHVSHE-FHVSTE) &
          *WTEABPs1(NB,NZ)+(FHVSHG-FHVSTG)*WTGRBPs1(NB,NZ)
        WTHTG=WTHTG+(1._r8-FHVSTG)*WTGRBs1(NB,NZ)
        WTHNG=WTHNG+(1._r8-FHVSTG)*WTGRBNs1(NB,NZ)
        WTHPG=WTHPG+(1._r8-FHVSTG)*WTGRBPs1(NB,NZ)
!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
        WTHSKBs1(NB,NZ)=FHVSTH*WTHSKBs1(NB,NZ)
        WTEARBs1(NB,NZ)=FHVSTE*WTEARBs1(NB,NZ)
        WTGRBs1(NB,NZ)=FHVSTG*WTGRBs1(NB,NZ)
        WTHSBNs1(NB,NZ)=FHVSTH*WTHSBNs1(NB,NZ)
        WTEABNs1(NB,NZ)=FHVSTE*WTEABNs1(NB,NZ)
        WTGRBNs1(NB,NZ)=FHVSTG*WTGRBNs1(NB,NZ)
        WTHSBPs1(NB,NZ)=FHVSTH*WTHSBPs1(NB,NZ)
        WTEABPs1(NB,NZ)=FHVSTE*WTEABPs1(NB,NZ)
        WTGRBPs1(NB,NZ)=FHVSTG*WTGRBPs1(NB,NZ)
        GRNXBs1(NB,NZ)=FHVSTG*GRNXBs1(NB,NZ)
        GRNOBs1(NB,NZ)=FHVSTG*GRNOBs1(NB,NZ)
        GRWTBs1(NB,NZ)=FHVSTG*GRWTBs1(NB,NZ)
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
!     VOLWOUs1,UVOLO=accumulated water loss for water balance calculation
!
        CPOOLK(NB,NZ)=0._r8
        DO 1325 K=1,JNODS1
          CPOOLK(NB,NZ)=CPOOLK(NB,NZ) &
            +CPOOL3s1(K,NB,NZ)+CPOOL4s1(K,NB,NZ) &
            +CO2Bs1(K,NB,NZ)+HCOBs1(K,NB,NZ)
1325    CONTINUE
        WTLSBs1(NB,NZ)=AMAX1(0.0,WTLFBs1(NB,NZ) &
          +WTSHEBs1(NB,NZ))
        WTSHTBs1(NB,NZ)=AMAX1(0.0,WTLFBs1(NB,NZ) &
          +WTSHEBs1(NB,NZ)+WTSTKBs1(NB,NZ)+WTRSVBs1(NB,NZ) &
          +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ)+WTGRBs1(NB,NZ) &
          +CPOOLs1(NB,NZ)+CPOOLK(NB,NZ))
        WTSHTNs1(NB,NZ)=AMAX1(0.0,WTLFBNs1(NB,NZ) &
          +WTSHBNs1(NB,NZ)+WTSTBNs1(NB,NZ)+WTRSBNs1(NB,NZ) &
          +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ)+WTGRBNs1(NB,NZ) &
          +ZPOOLs1(NB,NZ))
        WTSHTPs1(NB,NZ)=AMAX1(0.0,WTLFBPs1(NB,NZ) &
          +WTSHBPs1(NB,NZ)+WTSTBPs1(NB,NZ)+WTRSBPs1(NB,NZ) &
          +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ)+WTGRBPs1(NB,NZ) &
          +PPOOLs1(NB,NZ))
        VOLWPX=VOLWPs1(NZ)
        WVPLT=AMAX1(0.0_r8,WTLSs1(NZ)+WVSTKs1(NZ))
        APSILT=ABS(PSILTs1(NZ))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWPs1(NZ)=1.0E-06_r8*WVPLT/FDM
        VOLWOUs1=VOLWOUs1+VOLWPX-VOLWPs1(NZ)
        UVOLOs1=UVOLOs1+VOLWPX-VOLWPs1(NZ)
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
!     IDAYs1(1,=emergence date
!     GROUP=node number required for floral initiation
!     PSTGI=node number at floral initiation
!     PSTGF=node number at flowering
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     IFLGA=flag for initializing leafout
!
        IF((IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1) &
          .AND.(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6) &
          .AND.ZCs1(NZ).GT.HVSTs1(NZ))THEN
          IF((IWTYPs1(NZ).NE.0.AND.VRNFs1(NB,NZ) &
            .LE.FVRNs1(IWTYPs1(NZ))*VRNXs1(NB,NZ)) &
            .OR.(IWTYPs1(NZ).EQ.0 &
            .AND.IDAYs1(1,NB,NZ).NE.0))THEN
            GROUPs1(NB,NZ)=GROUPIs1(NZ)
            PSTGIs1(NB,NZ)=PSTGs1(NB,NZ)
            PSTGFs1(NB,NZ)=0._r8
            VSTGXs1(NB,NZ)=0._r8
            TGSTGIs1(NB,NZ)=0._r8
            TGSTGFs1(NB,NZ)=0._r8
            FLG4s1(NB,NZ)=0._r8
            IDAYs1(1,NB,NZ)=I
            DO 3005 M=2,10
              IDAYs1(M,NB,NZ)=0
3005        CONTINUE
            IFLGAs1(NB,NZ)=0
            IF(NB.EQ.NB1s1(NZ))THEN
              DO 3010 NBX=1,NBRs1(NZ)
                IF(NBX.NE.NB1s1(NZ))THEN
                  GROUPs1(NBX,NZ)=GROUPIs1(NZ)
                  PSTGIs1(NBX,NZ)=PSTGs1(NBX,NZ)
                  PSTGFs1(NBX,NZ)=0._r8
                  VSTGXs1(NBX,NZ)=0._r8
                  TGSTGIs1(NBX,NZ)=0._r8
                  TGSTGFs1(NBX,NZ)=0._r8
                  FLG4s1(NBX,NZ)=0._r8
                  IDAYs1(1,NBX,NZ)=I
                  DO 3015 M=2,10
                    IDAYs1(M,NBX,NZ)=0
3015              CONTINUE
                  IFLGAs1(NBX,NZ)=0
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
      IF(JHVSTs1(NZ).NE.0)IDTHBs1(NB,NZ)=1
      IF(PPs1(NZ).LE.0.0)IDTHBs1(NB,NZ)=1
9835  CONTINUE
      WTLSs1(NZ)=0._r8
      WTSTKs1(NZ)=0._r8
      WVSTKs1(NZ)=0._r8
      ARSTPs1(NZ)=0._r8
      DO 9840 NB=1,NBRs1(NZ)
        WTLSs1(NZ)=WTLSs1(NZ)+WTLSBs1(NB,NZ)
        WTSTKs1(NZ)=WTSTKs1(NZ)+WTSTKBs1(NB,NZ)
        WVSTKs1(NZ)=WVSTKs1(NZ)+WVSTKBs1(NB,NZ)
        DO 9830 L=1,JC1
          ARSTPs1(NZ)=ARSTPs1(NZ)+ARSTKs1(L,NB,NZ)
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
      IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
        XHVST=1.0_r8-THINs1(NZ)
        DO 3985 N=1,MYs1(NZ)
          DO 3980 L=NUs1,NJs1
            IF(IHVSTs1(NZ).NE.5)THEN
              XHVST=1.0_r8-THINs1(NZ)
              XHVSN=XHVST
              XHVSP=XHVST
              FFIRE=0._r8
              FFIRN=0._r8
              FFIRP=0._r8
            ELSE
              IF(THETWs1(L).GT.FVLWB.OR.CORGCs1(L).LE.FORGC &
                .OR.ITILLs1.NE.22)THEN
                XHVST=1.0_r8
                XHVSN=XHVST
                XHVSP=XHVST
                FFIRE=0._r8
                FFIRN=0._r8
                FFIRP=0._r8
              ELSE
                XHVST=1.0_r8-DCORPs1*EHVSTs1(1,3,NZ) &
                  *AMIN1(1.0,(CORGCs1(L)-FORGC)/(0.55E+06-FORGC))
                XHVSN=XHVST
                XHVSP=XHVST
                FFIRE=EHVSTs1(2,3,NZ)
                FFIRN=FFIRE*EFIRE(1,IHVSTs1(NZ))
                FFIRP=FFIRE*EFIRE(2,IHVSTs1(NZ))
              ENDIF
            ENDIF
            DO 3385 M=1,4
              FHVST=(1._r8-XHVST)*CFOPCs1(0,M,NZ)*CPOOLRs1(N,L,NZ)
              FHVSN=(1._r8-XHVSN)*CFOPNs1(0,M,NZ)*ZPOOLRs1(N,L,NZ)
              FHVSP=(1._r8-XHVSP)*CFOPPs1(0,M,NZ)*PPOOLRs1(N,L,NZ)
              CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-FFIRE)*FHVST
              ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-FFIRN)*FHVSN
              PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-FFIRP)*FHVSP
              VCO2Fs1(NZ)=VCO2Fs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
              VCH4Fs1(NZ)=VCH4Fs1(NZ)-FCH4F*FFIRE*FHVST
              VOXYFs1(NZ)=VOXYFs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
              VNH3Fs1(NZ)=VNH3Fs1(NZ)-FFIRN*FHVSN
              VN2OFs1(NZ)=VN2OFs1(NZ)-0.0
              VPO4Fs1(NZ)=VPO4Fs1(NZ)-FFIRP*FHVSP
              CNETs1(NZ)=CNETs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
              TNBPs1=TNBPs1-FCH4F*FFIRE*FHVST
              DO NR=1,NRTs1(NZ)
                FHVST=(1._r8-XHVST)*CFOPCs1(5,M,NZ)*(WTRT1s1(N,L,NR,NZ) &
                  +WTRT2s1(N,L,NR,NZ))*FWODRs1(0)
                FHVSN=(1._r8-XHVSN)*CFOPNs1(5,M,NZ)*(WTRT1Ns1(N,L,NR,NZ) &
                  +WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(0)
                FHVSP=(1._r8-XHVSP)*CFOPPs1(5,M,NZ)*(WTRT1Ps1(N,L,NR,NZ) &
                  +WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(0)
                CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-FFIRE)*FHVST
                ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-FFIRN)*FHVSN
                PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-FFIRP)*FHVSP
                VCO2Fs1(NZ)=VCO2Fs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                VCH4Fs1(NZ)=VCH4Fs1(NZ)-FCH4F*FFIRE*FHVST
                VOXYFs1(NZ)=VOXYFs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
                VNH3Fs1(NZ)=VNH3Fs1(NZ)-FFIRN*FHVSN
                VN2OFs1(NZ)=VN2OFs1(NZ)-0.0
                VPO4Fs1(NZ)=VPO4Fs1(NZ)-FFIRP*FHVSP
                CNETs1(NZ)=CNETs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                TNBPs1=TNBPs1-FCH4F*FFIRE*FHVST
                FHVST=(1._r8-XHVST)*CFOPCs1(4,M,NZ)*(WTRT1s1(N,L,NR,NZ) &
                  +WTRT2s1(N,L,NR,NZ))*FWODRs1(1)
                FHVSN=(1._r8-XHVSN)*CFOPNs1(4,M,NZ)*(WTRT1Ns1(N,L,NR,NZ) &
                  +WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(1)
                FHVSP=(1._r8-XHVSP)*CFOPPs1(4,M,NZ)*(WTRT1Ps1(N,L,NR,NZ) &
                  +WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(1)
                CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-FFIRE)*FHVST
                ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-FFIRN)*FHVSN
                PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-FFIRP)*FHVSP
                VCO2Fs1(NZ)=VCO2Fs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                VCH4Fs1(NZ)=VCH4Fs1(NZ)-FCH4F*FFIRE*FHVST
                VOXYFs1(NZ)=VOXYFs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST*2.667
                VNH3Fs1(NZ)=VNH3Fs1(NZ)-FFIRN*FHVSN
                VN2OFs1(NZ)=VN2OFs1(NZ)-0.0
                VPO4Fs1(NZ)=VPO4Fs1(NZ)-FFIRP*FHVSP
                CNETs1(NZ)=CNETs1(NZ)-(1._r8-FCH4F)*FFIRE*FHVST
                TNBPs1=TNBPs1-FCH4F*FFIRE*FHVST
              enddo
3385        CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING HARVESTING
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Zs1(NZ)=RCO2Zs1(NZ)-(1._r8-XHVST) &
              *(CO2As1(N,L,NZ)+CO2Ps1(N,L,NZ))
            ROXYZs1(NZ)=ROXYZs1(NZ)-(1._r8-XHVST) &
              *(OXYAs1(N,L,NZ)+OXYPs1(N,L,NZ))
            RCH4Zs1(NZ)=RCH4Zs1(NZ)-(1._r8-XHVST) &
              *(CH4As1(N,L,NZ)+CH4Ps1(N,L,NZ))
            RN2OZs1(NZ)=RN2OZs1(NZ)-(1._r8-XHVST) &
              *(Z2OAs1(N,L,NZ)+Z2OPs1(N,L,NZ))
            RNH3Zs1(NZ)=RNH3Zs1(NZ)-(1._r8-XHVST) &
              *(ZH3As1(N,L,NZ)+ZH3Ps1(N,L,NZ))
            RH2GZs1(NZ)=RH2GZs1(NZ)-(1._r8-XHVST) &
              *(H2GAs1(N,L,NZ)+H2GPs1(N,L,NZ))
            CO2As1(N,L,NZ)=XHVST*CO2As1(N,L,NZ)
            OXYAs1(N,L,NZ)=XHVST*OXYAs1(N,L,NZ)
            CH4As1(N,L,NZ)=XHVST*CH4As1(N,L,NZ)
            Z2OAs1(N,L,NZ)=XHVST*Z2OAs1(N,L,NZ)
            ZH3As1(N,L,NZ)=XHVST*ZH3As1(N,L,NZ)
            H2GAs1(N,L,NZ)=XHVST*H2GAs1(N,L,NZ)
            CO2Ps1(N,L,NZ)=XHVST*CO2Ps1(N,L,NZ)
            OXYPs1(N,L,NZ)=XHVST*OXYPs1(N,L,NZ)
            CH4Ps1(N,L,NZ)=XHVST*CH4Ps1(N,L,NZ)
            Z2OPs1(N,L,NZ)=XHVST*Z2OPs1(N,L,NZ)
            ZH3Ps1(N,L,NZ)=XHVST*ZH3Ps1(N,L,NZ)
            H2GPs1(N,L,NZ)=XHVST*H2GPs1(N,L,NZ)
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
            DO 3960 NR=1,NRTs1(NZ)
              WTRT1s1(N,L,NR,NZ)=WTRT1s1(N,L,NR,NZ)*XHVST
              WTRT2s1(N,L,NR,NZ)=WTRT2s1(N,L,NR,NZ)*XHVST
              WTRT1Ns1(N,L,NR,NZ)=WTRT1Ns1(N,L,NR,NZ)*XHVSN
              WTRT2Ns1(N,L,NR,NZ)=WTRT2Ns1(N,L,NR,NZ)*XHVSN
              WTRT1Ps1(N,L,NR,NZ)=WTRT1Ps1(N,L,NR,NZ)*XHVSP
              WTRT2Ps1(N,L,NR,NZ)=WTRT2Ps1(N,L,NR,NZ)*XHVSP
              RTWT1s1(N,NR,NZ)=RTWT1s1(N,NR,NZ)*XHVST
              RTWT1Ns1(N,NR,NZ)=RTWT1Ns1(N,NR,NZ)*XHVST
              RTWT1Ps1(N,NR,NZ)=RTWT1Ps1(N,NR,NZ)*XHVST
              RTLG1s1(N,L,NR,NZ)=RTLG1s1(N,L,NR,NZ)*XHVST
              RTLG2s1(N,L,NR,NZ)=RTLG2s1(N,L,NR,NZ)*XHVST
              RTN2s1(N,L,NR,NZ)=RTN2s1(N,L,NR,NZ)*XHVST
3960        CONTINUE
            CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)*XHVST
            ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)*XHVSN
            PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)*XHVSP
            WTRTLs1(N,L,NZ)=WTRTLs1(N,L,NZ)*XHVST
            WTRTDs1(N,L,NZ)=WTRTDs1(N,L,NZ)*XHVST
            WSRTLs1(N,L,NZ)=WSRTLs1(N,L,NZ)*XHVST
            RTN1s1(N,L,NZ)=RTN1s1(N,L,NZ)*XHVST
            RTNLs1(N,L,NZ)=RTNLs1(N,L,NZ)*XHVST
            RTLGPs1(N,L,NZ)=RTLGPs1(N,L,NZ)*XHVST
            RTDNPs1(N,L,NZ)=RTDNPs1(N,L,NZ)*XHVST
            RTVLPs1(N,L,NZ)=RTVLPs1(N,L,NZ)*XHVST
            RTVLWs1(N,L,NZ)=RTVLWs1(N,L,NZ)*XHVST
            RTARPs1(N,L,NZ)=RTARPs1(N,L,NZ)*XHVST
            RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)*XHVST
            RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)*XHVST
            RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)*XHVST
!
!     NODULE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(INTYPs1(NZ).NE.0.AND.N.EQ.1)THEN
              DO 3395 M=1,4
                CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+(1._r8-XHVST) &
                  *(CFOPCs1(4,M,NZ)*WTNDLs1(L,NZ) &
                  +CFOPCs1(0,M,NZ)*CPOOLNs1(L,NZ))
                ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+(1._r8-XHVSN) &
                  *(CFOPNs1(4,M,NZ)*WTNDLNs1(L,NZ) &
                  +CFOPNs1(0,M,NZ)*ZPOOLNs1(L,NZ))
                PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+(1._r8-XHVSP) &
                  *(CFOPPs1(4,M,NZ)*WTNDLPs1(L,NZ) &
                  +CFOPPs1(0,M,NZ)*PPOOLNs1(L,NZ))
3395          CONTINUE
              WTNDLs1(L,NZ)=WTNDLs1(L,NZ)*XHVST
              WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)*XHVSN
              WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)*XHVSP
              CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)*XHVST
              ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)*XHVSN
              PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)*XHVSP
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
        IF(ISTYPs1(NZ).NE.0)THEN
          DO 3400 M=1,4
            CSNCs1(M,0,NGs1(NZ),NZ)=CSNCs1(M,0,NGs1(NZ),NZ) &
              +((1._r8-XHVST)*CFOPCs1(0,M,NZ)*WTRVCs1(NZ))*FWOODs1(0)
            ZSNCs1(M,0,NGs1(NZ),NZ)=ZSNCs1(M,0,NGs1(NZ),NZ) &
              +((1._r8-XHVSN)*CFOPNs1(0,M,NZ)*WTRVNs1(NZ))*FWOODNs1(0)
            PSNCs1(M,0,NGs1(NZ),NZ)=PSNCs1(M,0,NGs1(NZ),NZ) &
              +((1._r8-XHVSP)*CFOPPs1(0,M,NZ)*WTRVPs1(NZ))*FWOODPs1(0)
            CSNCs1(M,1,NGs1(NZ),NZ)=CSNCs1(M,1,NGs1(NZ),NZ) &
              +((1._r8-XHVST)*CFOPCs1(0,M,NZ)*WTRVCs1(NZ))*FWOODs1(1)
            ZSNCs1(M,1,NGs1(NZ),NZ)=ZSNCs1(M,1,NGs1(NZ),NZ) &
              +((1._r8-XHVSN)*CFOPNs1(0,M,NZ)*WTRVNs1(NZ))*FWOODNs1(1)
            PSNCs1(M,1,NGs1(NZ),NZ)=PSNCs1(M,1,NGs1(NZ),NZ) &
              +((1._r8-XHVSP)*CFOPPs1(0,M,NZ)*WTRVPs1(NZ))*FWOODPs1(1)
3400      CONTINUE
          WTRVCs1(NZ)=WTRVCs1(NZ)*XHVST
          WTRVNs1(NZ)=WTRVNs1(NZ)*XHVSN
          WTRVPs1(NZ)=WTRVPs1(NZ)*XHVSP
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByHarvest

end module PlantDisturbsMod
