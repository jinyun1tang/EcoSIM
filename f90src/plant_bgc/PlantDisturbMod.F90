module PlantDisturbMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
  use SOMDataType
  use GrosubPars
  use PhenologyDataType
  use GridDataType
implicit none
  private

  include "blkc.h"
  include "blk1cp.h"
  include "blk1cr.h"
  include "blk1g.h"
  include "blk1n.h"
  include "blk1p.h"
  include "blk1s.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk3.h"
  include "blk5.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk9a.h"
  include "blk9b.h"
  include "blk9c.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk12a.h"
  include "blk12b.h"
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk14.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"
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

  subroutine RemoveBiomByManagement(I,J,NZ,NY,NX,WTSHTA,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX
  REAL(R8), intent(in) :: WTSHTA(JZ,JY,JX)
  real(r8), intent(inout) :: CPOOLK(10,JP,JY,JX)
!     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
!

  call RemoveBiomByHarvest(I,J,NZ,NY,NX,WTSHTA,CPOOLK)
!
!     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
!
  call RemoveBiomByTillage(I,J,NZ,NY,NX,CPOOLK)
  end subroutine RemoveBiomByManagement
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomassByDisturbance(I,J,NZ,NY,NX,CPOOLK)
  implicit none
  integer , intent(in) :: I,J,NZ,NY,NX
  real(r8), INTENT(INOUT) :: CPOOLK(10,JP,JY,JX)
  real(r8) :: FHVST
  real(r8) :: FHVSH
  real(r8) :: WHVSTD
  integer :: M

!     begin_execution
!     IF(J.EQ.INT(ZNOON(NY,NX)))THEN
!        XHVST=1.0_r8
!        WHVSBL=0._r8
  WTHTH0=0._r8
  WTHNH0=0._r8
  WTHPH0=0._r8
  WTHTH1=0._r8
  WTHNH1=0._r8
  WTHPH1=0._r8
  WTHTH2=0._r8
  WTHNH2=0._r8
  WTHPH2=0._r8
  WTHTH3=0._r8
  WTHNH3=0._r8
  WTHPH3=0._r8
  WTHTH4=0._r8
  WTHNH4=0._r8
  WTHPH4=0._r8
  WTHTR1=0._r8
  WTHNR1=0._r8
  WTHPR1=0._r8
  WTHTR2=0._r8
  WTHNR2=0._r8
  WTHPR2=0._r8
  WTHTR3=0._r8
  WTHNR3=0._r8
  WTHPR3=0._r8
  WTHTR4=0._r8
  WTHNR4=0._r8
  WTHPR4=0._r8
  WTHTX0=0._r8
  WTHNX0=0._r8
  WTHPX0=0._r8
  WTHTX1=0._r8
  WTHNX1=0._r8
  WTHPX1=0._r8
  WTHTX2=0._r8
  WTHNX2=0._r8
  WTHPX2=0._r8
  WTHTX3=0._r8
  WTHNX3=0._r8
  WTHPX3=0._r8
  WTHTX4=0._r8
  WTHNX4=0._r8
  WTHPX4=0._r8
  WTHTG=0._r8
  WTHNG=0._r8
  WTHPG=0._r8
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
  IF(IHVST(NZ,I,NY,NX).GE.0)THEN
    IF(J.EQ.INT(ZNOON(NY,NX)).AND.IHVST(NZ,I,NY,NX).NE.4 &
      .AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
        FHVST=AMAX1(0.0,1.0-EHVST(1,4,NZ,I,NY,NX))
        FHVSH=FHVST
      ELSE
        FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
        IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
          FHVSH=AMAX1(0.0,1.0-EHVST(1,4,NZ,I,NY,NX)*THIN(NZ,I,NY,NX))
        ELSE
          FHVSH=FHVST
        ENDIF
      ENDIF
    ELSEIF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WTSTG(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        WHVSTD=HVST(NZ,I,NY,NX)*THIN(NZ,I,NY,NX)*0.45/24.0 &
          *AREA(3,NU(NY,NX),NY,NX)*EHVST(1,4,NZ,I,NY,NX)
        FHVST=AMAX1(0.0,1.0-WHVSTD/WTSTG(NZ,NY,NX))
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
      WTHTH4=WTHTH4+(1.0-FHVSH)*WTSTDG(M,NZ,NY,NX)
      WTHNH4=WTHNH4+(1.0-FHVSH)*WTSTDN(M,NZ,NY,NX)
      WTHPH4=WTHPH4+(1.0-FHVSH)*WTSTDP(M,NZ,NY,NX)
      WTHTX4=WTHTX4+(FHVSH-FHVST)*WTSTDG(M,NZ,NY,NX)
      WTHNX4=WTHNX4+(FHVSH-FHVST)*WTSTDN(M,NZ,NY,NX)
      WTHPX4=WTHPX4+(FHVSH-FHVST)*WTSTDP(M,NZ,NY,NX)
      WTSTDG(M,NZ,NY,NX)=FHVST*WTSTDG(M,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=FHVST*WTSTDN(M,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=FHVST*WTSTDP(M,NZ,NY,NX)
6475  CONTINUE
!
    call PlantDisturbance(I,J,NZ,NY,NX)

    ZEROP(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)
    ZEROQ(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
    ZEROL(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)*1.0E+06
  ENDIF
  end subroutine RemoveBiomassByDisturbance


!------------------------------------------------------------------------------------------
  subroutine PlantDisturbance(I,J,NZ,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX

  real(r8) :: WTHNL0,WTHPL0
  real(r8) :: WTHNL1,WTHPL1
  real(r8) :: WTHNL2,WTHPL2
  real(r8) :: WTHNL3,WTHPL3
  real(r8) :: WTHNL4,WTHPL4
  real(r8) :: WTHNRT,WTHPRT,WTHTXT,WTHNXT
  real(r8) :: WTHTR0,WTHNR0,WTHPR0
  real(r8) :: WTHTRT
  real(r8) :: WTHPXT

  call ApplyDisturbanceBiomRemoval(I,J,NZ,NY,NX,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4)
!
!     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call TotalBiomRemovalByDisturbance(I,J,NZ,NY,NX,WTHTR0,WTHNR0,WTHPR0,WTHPXT,WTHTRT,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)
!
!     ABOVE-GROUND LITTERFALL FROM HARVESTING
!
  call LiterfallByDisturbance(I,J,NZ,NY,NX,WTHPXT,WTHTRT,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)
  end subroutine PlantDisturbance
!------------------------------------------------------------------------------------------

  subroutine LiterfallByDisturbance(I,J,NZ,NY,NX,WTHPXT,WTHTRT,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)

  implicit none
  integer , intent(in) :: I,J,NZ,NY,NX
  real(r8), intent(in) :: WTHPXT,WTHTRT
  real(r8), intent(in) :: WTHTR0,WTHNR0,WTHPR0,WTHNL0,WTHPL0
  real(r8), intent(in) :: WTHNL1,WTHPL1
  real(r8), intent(in) :: WTHNL2,WTHPL2
  real(r8), intent(in) :: WTHNL3,WTHPL3
  real(r8), intent(in) :: WTHNL4,WTHPL4
  real(r8), intent(in) :: WTHNRT,WTHPRT,WTHTXT,WTHNXT
  integer :: m
!     begin_execution
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
  IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
    IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      DO 6375 M=1,4
        CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
          +CFOPC(0,M,NZ,NY,NX)*(WTHTR0+WTHTX0) &
          +CFOPC(1,M,NZ,NY,NX)*(WTHTR1+WTHTX1) &
          +CFOPC(2,M,NZ,NY,NX)*(WTHTR2+WTHTX2)
        ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
          +CFOPN(0,M,NZ,NY,NX)*(WTHNR0+WTHNX0) &
          +CFOPN(1,M,NZ,NY,NX)*(WTHNR1+WTHNX1) &
          +CFOPN(2,M,NZ,NY,NX)*(WTHNR2+WTHNX2)
        PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
          +CFOPP(0,M,NZ,NY,NX)*(WTHPR0+WTHPX0) &
          +CFOPP(1,M,NZ,NY,NX)*(WTHPR1+WTHPX1) &
          +CFOPP(2,M,NZ,NY,NX)*(WTHPR2+WTHPX2)
        IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
          CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
            +CFOPC(3,M,NZ,NY,NX)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
          ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
            +CFOPN(3,M,NZ,NY,NX)*(WTHNR3+WTHNX3+WTHNR4+WTHNX4)
          PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
            +CFOPP(3,M,NZ,NY,NX)*(WTHPR3+WTHPX3+WTHPR4+WTHPX4)
        ELSE
          WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX) &
            +CFOPC(5,M,NZ,NY,NX)*(WTHTX3+WTHTX4)
          WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX) &
            +CFOPN(5,M,NZ,NY,NX)*(WTHNX3+WTHNX4)
          WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX) &
            +CFOPP(5,M,NZ,NY,NX)*(WTHPX3+WTHPX4)
          CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX) &
            +CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTR4)*FWOOD(0)
          ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX) &
            +CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNR4)*FWOODN(0)
          PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX) &
            +CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPR4)*FWOODP(0)
          CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
            +CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTR4)*FWOOD(1)
          ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
            +CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNR4)*FWOODN(1)
          PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
            +CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPR4)*FWOODP(0)
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
        CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
          +CFOPC(0,M,NZ,NY,NX)*(WTHTR0+WTHTX0) &
          +CFOPC(1,M,NZ,NY,NX)*(WTHTR1+WTHTX1) &
          +CFOPC(2,M,NZ,NY,NX)*(WTHTR2+WTHTX2)
        ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
          +CFOPN(0,M,NZ,NY,NX)*WTHNL0 &
          +CFOPN(1,M,NZ,NY,NX)*WTHNL1 &
          +CFOPN(2,M,NZ,NY,NX)*WTHNL2
        PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
          +CFOPP(0,M,NZ,NY,NX)*WTHPL0 &
          +CFOPP(1,M,NZ,NY,NX)*WTHPL1 &
          +CFOPP(2,M,NZ,NY,NX)*WTHPL2
        ZSNC(4,1,0,NZ,NY,NX)=ZSNC(4,1,0,NZ,NY,NX) &
          +CFOPN(0,M,NZ,NY,NX)*(WTHNR0+WTHNX0-WTHNL0) &
          +CFOPN(1,M,NZ,NY,NX)*(WTHNR1+WTHNX1-WTHNL1) &
          +CFOPN(2,M,NZ,NY,NX)*(WTHNR2+WTHNX2-WTHNL2)
        PSNC(4,1,0,NZ,NY,NX)=PSNC(4,1,0,NZ,NY,NX) &
          +CFOPP(0,M,NZ,NY,NX)*(WTHPR0+WTHPX0-WTHPL0) &
          +CFOPP(1,M,NZ,NY,NX)*(WTHPR1+WTHPX1-WTHPL1) &
          +CFOPP(2,M,NZ,NY,NX)*(WTHPR2+WTHPX2-WTHPL2)
        IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
          CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX) &
            +CFOPC(3,M,NZ,NY,NX)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
          ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX) &
            +CFOPN(3,M,NZ,NY,NX)*(WTHNL3+WTHNL4)
          PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX) &
            +CFOPP(3,M,NZ,NY,NX)*(WTHPL3+WTHPL4)
          ZSNC(4,1,0,NZ,NY,NX)=ZSNC(4,1,0,NZ,NY,NX) &
            +CFOPN(3,M,NZ,NY,NX)*(WTHNR3+WTHNX3-WTHNL3+WTHNR4+WTHNX4-WTHNL4)
          PSNC(4,1,0,NZ,NY,NX)=PSNC(4,1,0,NZ,NY,NX) &
            +CFOPP(3,M,NZ,NY,NX)*(WTHPR3+WTHPX3-WTHPL3+WTHPR4+WTHPX4-WTHPL4)
        ELSE
          WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTX3)
          WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)*WTHNL3
          WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)*WTHPL3
          CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)*CFOPC(3,M,NZ,NY,NX)*(WTHTR4+WTHTX4)*FWOOD(0)
          ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)*WTHNL4*FWOODN(0)
          PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)*WTHPL4*FWOODP(0)
          ZSNC(4,0,0,NZ,NY,NX)=ZSNC(4,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNX3-WTHNL3 &
            +WTHNR4+WTHNX4-WTHNL4)*FWOODN(0)
          PSNC(4,0,0,NZ,NY,NX)=PSNC(4,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPX3-WTHPL3 &
            +WTHPR4+WTHPX4-WTHPL4)*FWOODP(0)
          CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)*(WTHTR4+WTHTX4)*FWOOD(1)
          ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)*WTHNL4*FWOODN(1)
          PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)*WTHPL4*FWOODP(1)
          ZSNC(4,1,0,NZ,NY,NX)=ZSNC(4,1,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNX3-WTHNL3 &
            +WTHNR4+WTHNX4-WTHNL4)*FWOODN(1)
          PSNC(4,1,0,NZ,NY,NX)=PSNC(4,1,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPX3-WTHPL3 &
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
    TCSNC(NZ,NY,NX)=TCSNC(NZ,NY,NX)+WTHTRT+WTHTXT
    TZSNC(NZ,NY,NX)=TZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
    TPSNC(NZ,NY,NX)=TPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
    TCSN0(NZ,NY,NX)=TCSN0(NZ,NY,NX)+WTHTRT+WTHTXT
    TZSN0(NZ,NY,NX)=TZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
    TPSN0(NZ,NY,NX)=TPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
  ENDIF
  end subroutine LiterfallByDisturbance
!------------------------------------------------------------------------------------------

  subroutine TotalBiomRemovalByDisturbance(I,J,NZ,NY,NX,WTHTR0,WTHNR0,WTHPR0,WTHPXT,WTHTRT,&
    WTHNRT,WTHPRT,WTHTXT,WTHNXT)
  implicit none
  integer , intent(in) :: I,J,NZ,NY,NX
  real(r8), intent(in) :: WTHTR0,WTHNR0,WTHPR0
  real(r8), intent(out):: WTHPXT,WTHTRT,WTHNRT,WTHPRT,WTHTXT,WTHNXT
  real(r8) :: WTHTHT,WTHNHT,WTHPHT
!     begin_execution
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

  IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
    IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      IF(JHVST(NZ,I,NY,NX).NE.2)THEN
        HVSTC(NZ,NY,NX)=HVSTC(NZ,NY,NX)+WTHTHT-WTHTRT
        HVSTN(NZ,NY,NX)=HVSTN(NZ,NY,NX)+WTHNHT-WTHNRT
        HVSTP(NZ,NY,NX)=HVSTP(NZ,NY,NX)+WTHPHT-WTHPRT
        TNBP(NY,NX)=TNBP(NY,NX)+WTHTRT-WTHTHT
        XHVSTC(NY,NX)=XHVSTC(NY,NX)+WTHTHT-WTHTRT
        XHVSTN(NY,NX)=XHVSTN(NY,NX)+WTHNHT-WTHNRT
        XHVSTP(NY,NX)=XHVSTP(NY,NX)+WTHPHT-WTHPRT
      ELSE
        WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+WTHTHT-WTHTRT
        WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+WTHNHT-WTHNRT
        WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+WTHPHT-WTHPRT
      ENDIF
!
!     C,N,P LOST AS GAS IF FIRE
!
!     VCO2F,VCH4F,VOXYF,VNH3F,VN2OF,VPO4F=CO2,CH4,O2,NH3,N2O,PO4 emission from disturbance
!     CNET=PFT net CO2 fixation
!     TNBP=total net biome productivity
!
    ELSE
      VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*(WTHTHT-WTHTRT)
      VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*(WTHTHT-WTHTRT)
      VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*(WTHTHT-WTHTRT)*2.667
      VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-WTHNHT+WTHNRT
      VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
      VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-WTHPHT+WTHPRT
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*(WTHTHT-WTHTRT)
      TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*(WTHTHT-WTHTRT)
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
      HVSTC(NZ,NY,NX)=HVSTC(NZ,NY,NX)+GY*(WTHTHT-WTHTRT)
      HVSTN(NZ,NY,NX)=HVSTN(NZ,NY,NX)+WTHNHT-WTHNRT
      HVSTP(NZ,NY,NX)=HVSTP(NZ,NY,NX)+WTHPHT-WTHPRT
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-GZ*(WTHTHT-WTHTRT)
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-GZ*(WTHTHT-WTHTRT)
!     TNBP(NY,NX)=TNBP(NY,NX)+GY*(WTHTRT-WTHTHT)
!     CNET(NZ,NY,NX)=CNET(NZ,NY,NX)+GZ*(WTHTRT-WTHTHT)
      XHVSTC(NY,NX)=XHVSTC(NY,NX)+GY*(WTHTHT-WTHTRT)
      XHVSTN(NY,NX)=XHVSTN(NY,NX)+WTHNHT-WTHNRT
      XHVSTP(NY,NX)=XHVSTP(NY,NX)+WTHPHT-WTHPRT
      RECO(NY,NX)=RECO(NY,NX)-GZ*(WTHTHT-WTHTRT)
      TRAU(NY,NX)=TRAU(NY,NX)-GZ*(WTHTHT-WTHTRT)
  ENDIF
  end subroutine TotalBiomRemovalByDisturbance
!------------------------------------------------------------------------------------------

  subroutine ApplyDisturbanceBiomRemoval(I,J,NZ,NY,NX,WTHTR0,WTHNR0,WTHPR0,&
    WTHNL0,WTHPL0,WTHNL1,WTHPL1,WTHNL2,WTHPL2,WTHNL3,WTHPL3,WTHNL4,WTHPL4)
  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX
  real(r8), intent(out) :: WTHTR0,WTHNR0,WTHPR0
  real(r8), intent(out) :: WTHNL0,WTHPL0
  real(r8), intent(out) :: WTHNL1,WTHPL1
  real(r8), intent(out) :: WTHNL2,WTHPL2
  real(r8), intent(out) :: WTHNL3,WTHPL3
  real(r8), intent(out) :: WTHNL4,WTHPL4

!     begin_execution
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
  IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
    WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
!
!     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ,I,NY,NX).EQ.1)THEN
    WTHTR0=WTHTH0
    WTHNR0=WTHNH0
    WTHPR0=WTHPH0
    WTHTR1=WTHTH1
    WTHNR1=WTHNH1
    WTHPR1=WTHPH1
    WTHTR2=WTHTH2-WTHTG*EHVST(2,2,NZ,I,NY,NX)
    WTHNR2=WTHNH2-WTHNG*EHVST(2,2,NZ,I,NY,NX)
    WTHPR2=WTHPH2-WTHPG*EHVST(2,2,NZ,I,NY,NX)
    WTHTR3=WTHTH3
    WTHNR3=WTHNH3
    WTHPR3=WTHPH3
    WTHTR4=WTHTH4
    WTHNR4=WTHNH4
    WTHPR4=WTHPH4
!
!     IF ONLY WOOD C,N,P REMOVED AT HARVEST
!
  ELSEIF(IHVST(NZ,I,NY,NX).EQ.2)THEN
    WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
!
!     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
!
  ELSEIF(IHVST(NZ,I,NY,NX).EQ.3)THEN
    WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
!
!     IF PLANT C,N,P REMOVED BY GRAZING
!
  ELSEIF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
    WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
    WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
    WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
    WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX)*0.5)
    WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX)*0.5)
    WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX)*0.5)
    WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX)*0.5)
    WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX)*0.5)
    WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX)*0.5)
    WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX)*0.5)
!
!     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
!
!     FERT=fertilizer type from fertilizer input file
!     IYTYP=fertilizer release type from fertilizer input file
!
    FERT(17,I+1,NY,NX)=FERT(17,I+1,NY,NX) &
      +(WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4)/AREA(3,NU(NY,NX),NY,NX)
    FERT(18,I+1,NY,NX)=FERT(18,I+1,NY,NX) &
      +(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA(3,NU(NY,NX),NY,NX)*0.5
    FERT(3,I+1,NY,NX)=FERT(3,I+1,NY,NX) &
      +(WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4)/AREA(3,NU(NY,NX),NY,NX)*0.5
    FERT(19,I+1,NY,NX)=FERT(19,I+1,NY,NX) &
      +(WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4)/AREA(3,NU(NY,NX),NY,NX)
    IYTYP(2,I+1,NY,NX)=3
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
  ELSEIF(IHVST(NZ,I,NY,NX).EQ.5)THEN
    WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR0=WTHNH0*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX))*EHVST(2,1,NZ,I,NY,NX))
    WTHPR0=WTHPH0*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX))*EHVST(2,1,NZ,I,NY,NX))
    WTHNL0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPL0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHNR1=WTHNH1*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX))*EHVST(2,1,NZ,I,NY,NX))
    WTHPR1=WTHPH1*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX))*EHVST(2,1,NZ,I,NY,NX))
    WTHNL1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHPL1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
    WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHNR2=WTHNH2*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX))*EHVST(2,2,NZ,I,NY,NX))
    WTHPR2=WTHPH2*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX))*EHVST(2,2,NZ,I,NY,NX))
    WTHNL2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHPL2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
    WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHNR3=WTHNH3*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX))*EHVST(2,3,NZ,I,NY,NX))
    WTHPR3=WTHPH3*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX))*EHVST(2,3,NZ,I,NY,NX))
    WTHNL3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHPL3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
    WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHNR4=WTHNH4*(1.0-EFIRE(1,IHVST(NZ,I,NY,NX))*EHVST(2,4,NZ,I,NY,NX))
    WTHPR4=WTHPH4*(1.0-EFIRE(2,IHVST(NZ,I,NY,NX))*EHVST(2,4,NZ,I,NY,NX))
    WTHNL4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
    WTHPL4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
  ENDIF
  end subroutine ApplyDisturbanceBiomRemoval

!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByTillage(I,J,NZ,NY,NX,CPOOLK)

  implicit none
  integer , intent(in) :: I,J,NZ,NY,NX
  real(r8), intent(inout) :: CPOOLK(10,JP,JY,JX)
  integer :: L,K,M,N,NR,NB
  real(r8) :: XHVST
  REAL(R8) :: APSILT
  real(r8) :: FDM,VOLWPX
  real(r8) :: WVPLT
!     begin_execution
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
  IF(J.EQ.INT(ZNOON(NY,NX)).AND.(IBTYP(NZ,NY,NX).EQ.0 &
    .OR.IGTYP(NZ,NY,NX).LE.1).AND.(I.NE.IDAY0(NZ,NY,NX) &
    .OR.IYRC.NE.IYR0(NZ,NY,NX)))THEN
    IF(ITILL(I,NY,NX).LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.IDAY0(NZ,NY,NX).OR.IYRC.GT.IYR0(NZ,NY,NX))THEN
        XHVST=XCORP(NY,NX)
        PPX(NZ,NY,NX)=PPX(NZ,NY,NX)*XHVST
        PP(NZ,NY,NX)=PP(NZ,NY,NX)*XHVST
        FRADP(NZ,NY,NX)=FRADP(NZ,NY,NX)*XHVST
        VHCPC(NZ,NY,NX)=VHCPC(NZ,NY,NX)*XHVST
        WTLS(NZ,NY,NX)=0._r8
        WVSTK(NZ,NY,NX)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
        DO 8975 NB=1,NBR(NZ,NY,NX)
          IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
            IF(PP(NZ,NY,NX).LE.0.0)IDTHB(NB,NZ,NY,NX)=1
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
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                *(CFOPC(0,M,NZ,NY,NX)*(CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX) &
                +CPOOLK(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)) &
                +CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1) &
                +WTNDB(NB,NZ,NY,NX)) &
                +CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1) &
                +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)))
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0) &
                +WTSHEB(NB,NZ,NY,NX)*FWODB(0))
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                *(CFOPN(0,M,NZ,NY,NX)*(ZPOOL(NB,NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX) &
                +WTRSBN(NB,NZ,NY,NX)) &
                +CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1) &
                +WTNDBN(NB,NZ,NY,NX)) &
                +CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1) &
                +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)))
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0) &
                +WTSHBN(NB,NZ,NY,NX)*FWODSN(0))
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                *(CFOPP(0,M,NZ,NY,NX)*(PPOOL(NB,NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX) &
                +WTRSBP(NB,NZ,NY,NX)) &
                +CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1) &
                +WTNDBP(NB,NZ,NY,NX)) &
                +CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1) &
                +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)))
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0) &
                +WTSHBP(NB,NZ,NY,NX)*FWODSP(0))
              IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
                WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
                WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
                WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
              ELSE
                CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
                ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
                PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
              ENDIF
              CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPC(5,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)*FWOOD(0)
              ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPN(5,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)*FWOODN(0)
              PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPP(5,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)*FWOODP(0)
              CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)*FWOOD(1)
              ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)*FWOODN(1)
              PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)*FWOODP(1)
6380        CONTINUE
!
!     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
!
            CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)*XHVST
            CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX)*XHVST
            ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)*XHVST
            PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)*XHVST
            CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)*XHVST
            ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)*XHVST
            PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)*XHVST
            WTSHTB(NB,NZ,NY,NX)=WTSHTB(NB,NZ,NY,NX)*XHVST
            WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)*XHVST
            WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)*XHVST
            WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)*XHVST
            WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)*XHVST
            WVSTKB(NB,NZ,NY,NX)=WVSTKB(NB,NZ,NY,NX)*XHVST
            WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)*XHVST
            WTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)*XHVST
            WTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)*XHVST
            WTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)*XHVST
            WTSHTN(NB,NZ,NY,NX)=WTSHTN(NB,NZ,NY,NX)*XHVST
            WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)*XHVST
            WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)*XHVST
            WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)*XHVST
            WTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)*XHVST
            WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)*XHVST
            WTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)*XHVST
            WTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)*XHVST
            WTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)*XHVST
            WTSHTP(NB,NZ,NY,NX)=WTSHTP(NB,NZ,NY,NX)*XHVST
            WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)*XHVST
            WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)*XHVST
            WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)*XHVST
            WTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)*XHVST
            WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)*XHVST
            WTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)*XHVST
            WTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)*XHVST
            WTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)*XHVST
            GRNXB(NB,NZ,NY,NX)=GRNXB(NB,NZ,NY,NX)*XHVST
            GRNOB(NB,NZ,NY,NX)=GRNOB(NB,NZ,NY,NX)*XHVST
            GRWTB(NB,NZ,NY,NX)=GRWTB(NB,NZ,NY,NX)*XHVST
            ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)*XHVST
            WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)+WTSHEB(NB,NZ,NY,NX))
            WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
            WTSTXB(NB,NZ,NY,NX)=WTSTXB(NB,NZ,NY,NX)*XHVST
            WTSTXN(NB,NZ,NY,NX)=WTSTXN(NB,NZ,NY,NX)*XHVST
            WTSTXP(NB,NZ,NY,NX)=WTSTXP(NB,NZ,NY,NX)*XHVST
            WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
            DO 8970 K=0,25
              IF(K.NE.0)THEN
                CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)*XHVST
                CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)*XHVST
                CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)*XHVST
                HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)*XHVST
              ENDIF
              ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)*XHVST
              WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)*XHVST
              WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)*XHVST
!     HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)*XHVST
              WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)*XHVST
              WSSHE(K,NB,NZ,NY,NX)=WSSHE(K,NB,NZ,NY,NX)*XHVST
!     HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX)*XHVST
!     HTNODX(K,NB,NZ,NY,NX)=HTNODX(K,NB,NZ,NY,NX)*XHVST
              WGNODE(K,NB,NZ,NY,NX)=WGNODE(K,NB,NZ,NY,NX)*XHVST
              WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)*XHVST
              WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)*XHVST
              WGNODN(K,NB,NZ,NY,NX)=WGNODN(K,NB,NZ,NY,NX)*XHVST
              WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)*XHVST
              WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)*XHVST
              WGNODP(K,NB,NZ,NY,NX)=WGNODP(K,NB,NZ,NY,NX)*XHVST
              DO 8965 L=1,JC
                ARLFL(L,K,NB,NZ,NY,NX)=ARLFL(L,K,NB,NZ,NY,NX)*XHVST
                WGLFL(L,K,NB,NZ,NY,NX)=WGLFL(L,K,NB,NZ,NY,NX)*XHVST
                WGLFLN(L,K,NB,NZ,NY,NX)=WGLFLN(L,K,NB,NZ,NY,NX)*XHVST
                WGLFLP(L,K,NB,NZ,NY,NX)=WGLFLP(L,K,NB,NZ,NY,NX)*XHVST
8965          CONTINUE
8970        CONTINUE
          ENDIF
8975    CONTINUE
!
!     PSILT=canopy water potential
!     VOLWP=water volume in canopy
!     VOLWOU,UVOLO=accumulated water loss for water balance calculation
!
        VOLWPX=VOLWP(NZ,NY,NX)
        WVPLT=AMAX1(0.0_r8,WTLS(NZ,NY,NX)+WVSTK(NZ,NY,NX))
        APSILT=ABS(PSILT(NZ,NY,NX))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWP(NZ,NY,NX)=1.0E-06_r8*WVPLT/FDM
        VOLWOU=VOLWOU+VOLWPX-VOLWP(NZ,NY,NX)
        UVOLO(NY,NX)=UVOLO(NY,NX)+VOLWPX-VOLWP(NZ,NY,NX)
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
        IF(PP(NZ,NY,NX).LE.0.0)THEN
          IDTHR(NZ,NY,NX)=1
          IDTHP(NZ,NY,NX)=1
          IDTH(NZ,NY,NX)=1
          JHVST(NZ,I,NY,NX)=1
          IDAYH(NZ,NY,NX)=I
          IYRH(NZ,NY,NX)=IYRC
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
        DO 8985 N=1,MY(NZ,NY,NX)
          DO 8980 L=NU(NY,NX),NJ(NY,NX)
            DO 6385 M=1,4
              CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPC(0,M,NZ,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
              ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPN(0,M,NZ,NY,NX)*ZPOOLR(N,L,NZ,NY,NX)
              PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                *CFOPP(0,M,NZ,NY,NX)*PPOOLR(N,L,NZ,NY,NX)
              DO NR=1,NRT(NZ,NY,NX)
                CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPC(5,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
                  +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0)
                ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPN(5,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
                  +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0)
                PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPP(5,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
                  +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0)
                CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPC(4,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
                  +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1)
                ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPN(4,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
                  +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1)
                PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *CFOPP(4,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
                  +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1)
              ENDDO
6385        CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(1.0-XHVST) &
              *(CO2A(N,L,NZ,NY,NX)+CO2P(N,L,NZ,NY,NX))
            ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(1.0-XHVST) &
              *(OXYA(N,L,NZ,NY,NX)+OXYP(N,L,NZ,NY,NX))
            RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(1.0-XHVST) &
              *(CH4A(N,L,NZ,NY,NX)+CH4P(N,L,NZ,NY,NX))
            RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(1.0-XHVST) &
              *(Z2OA(N,L,NZ,NY,NX)+Z2OP(N,L,NZ,NY,NX))
            RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(1.0-XHVST) &
              *(ZH3A(N,L,NZ,NY,NX)+ZH3P(N,L,NZ,NY,NX))
            RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(1.0-XHVST) &
              *(H2GA(N,L,NZ,NY,NX)+H2GP(N,L,NZ,NY,NX))
            CO2A(N,L,NZ,NY,NX)=XHVST*CO2A(N,L,NZ,NY,NX)
            OXYA(N,L,NZ,NY,NX)=XHVST*OXYA(N,L,NZ,NY,NX)
            CH4A(N,L,NZ,NY,NX)=XHVST*CH4A(N,L,NZ,NY,NX)
            Z2OA(N,L,NZ,NY,NX)=XHVST*Z2OA(N,L,NZ,NY,NX)
            ZH3A(N,L,NZ,NY,NX)=XHVST*ZH3A(N,L,NZ,NY,NX)
            H2GA(N,L,NZ,NY,NX)=XHVST*H2GA(N,L,NZ,NY,NX)
            CO2P(N,L,NZ,NY,NX)=XHVST*CO2P(N,L,NZ,NY,NX)
            OXYP(N,L,NZ,NY,NX)=XHVST*OXYP(N,L,NZ,NY,NX)
            CH4P(N,L,NZ,NY,NX)=XHVST*CH4P(N,L,NZ,NY,NX)
            Z2OP(N,L,NZ,NY,NX)=XHVST*Z2OP(N,L,NZ,NY,NX)
            ZH3P(N,L,NZ,NY,NX)=XHVST*ZH3P(N,L,NZ,NY,NX)
            H2GP(N,L,NZ,NY,NX)=XHVST*H2GP(N,L,NZ,NY,NX)
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
            DO 8960 NR=1,NRT(NZ,NY,NX)
              WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*XHVST
              WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*XHVST
              WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*XHVST
              WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*XHVST
              WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*XHVST
              WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*XHVST
              RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*XHVST
              RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*XHVST
              RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*XHVST
              RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*XHVST
              RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*XHVST
              RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*XHVST
8960        CONTINUE
            CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*XHVST
            ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*XHVST
            PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*XHVST
            WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*XHVST
            WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*XHVST
            WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*XHVST
            RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*XHVST
            RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*XHVST
            RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*XHVST
            RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*XHVST
            RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*XHVST
            RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*XHVST
            RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*XHVST
            RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)*XHVST
            RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)*XHVST
            RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)*XHVST
!
!     LITTERFALL AND STATE VARIABLES FOR NODULES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
              DO 6395 M=1,4
                CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *(CFOPC(4,M,NZ,NY,NX)*WTNDL(L,NZ,NY,NX) &
                  +CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX))
                ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *(CFOPN(4,M,NZ,NY,NX)*WTNDLN(L,NZ,NY,NX) &
                  +CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX))
                PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *(CFOPP(4,M,NZ,NY,NX)*WTNDLP(L,NZ,NY,NX) &
                  +CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX))
6395          CONTINUE
              WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*XHVST
              WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*XHVST
              WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*XHVST
              CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*XHVST
              ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*XHVST
              PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*XHVST
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
          CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
            +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(0)
          ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
            +((1.0-XHVST)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(0)
          PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
            +((1.0-XHVST)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(0)
          CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
            +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(1)
          ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
            +((1.0-XHVST)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(1)
          PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
            +((1.0-XHVST)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(1)
6400    CONTINUE
        WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*XHVST
        WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*XHVST
        WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*XHVST
      ENDIF
    ENDIF
  ENDIF
  end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------

  subroutine RemoveBiomByHarvest(I,J,NZ,NY,NX,WTSHTA,CPOOLK)

  implicit none
  integer, intent(in) :: I,J,NZ,NY,NX
  REAL(R8), intent(in) :: WTSHTA(JZ,JY,JX)
  real(r8), intent(inout) :: CPOOLK(10,JP,JY,JX)
  integer :: L,K,M,NR,N,NB,NBX
  real(r8):: ZPOOLG,ZPOLNG,ZPOOLX
  real(r8) :: ZPOLNX,XHVSN,XHVSP,XHVST
  REAL(R8) :: WGLFBL(JZ,10,JP,JY,JX)
  real(r8) :: FHVSHK(0:25),FHVSTK(0:25)
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
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
  IF((IHVST(NZ,I,NY,NX).GE.0.AND.J.EQ.INT(ZNOON(NY,NX)) &
    .AND.IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6) &
    .OR.(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6))THEN
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
    IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(JHVST(NZ,I,NY,NX).NE.2)THEN
        PPX(NZ,NY,NX)=PPX(NZ,NY,NX)*(1.0-THIN(NZ,I,NY,NX))
        PP(NZ,NY,NX)=PP(NZ,NY,NX)*(1.0-THIN(NZ,I,NY,NX))
      ELSE
!     PPI(NZ,NY,NX)=AMAX1(1.0,0.5*(PPI(NZ,NY,NX)+GRNO(NZ,NY,NX)
!    2/AREA(3,NU(NY,NX),NY,NX)))
        PPX(NZ,NY,NX)=PPI(NZ,NY,NX)
        PP(NZ,NY,NX)=PPX(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      ENDIF
      IF(IHVST(NZ,I,NY,NX).EQ.3)THEN
        CF(NZ,NY,NX)=CF(NZ,NY,NX)*HVST(NZ,I,NY,NX)
      ENDIF
      IF(IHVST(NZ,I,NY,NX).LE.2.AND.HVST(NZ,I,NY,NX).LT.0.0)THEN
        ARLFY=(1.0-ABS(HVST(NZ,I,NY,NX)))*ARLFC(NY,NX)
        ARLFR=0._r8
        DO 9875 L=1,JC
          IF(ZL(L,NY,NX).GT.ZL(L-1,NY,NX).AND.ARLFT(L,NY,NX).GT.ZEROS(NY,NX) &
            .AND.ARLFR.LT.ARLFY)THEN
            IF(ARLFR+ARLFT(L,NY,NX).GT.ARLFY)THEN
              HVST(NZ,I,NY,NX)=ZL(L-1,NY,NX)+((ARLFY-ARLFR) &
                /ARLFT(L,NY,NX))*(ZL(L,NY,NX)-ZL(L-1,NY,NX))
            ENDIF
          ELSE
            HVST(NZ,I,NY,NX)=0._r8
          ENDIF
          ARLFR=ARLFR+ARLFT(L,NY,NX)
!     WRITE(*,6544)'HVST',I,J,L,NZ,IHVST(NZ,I,NY,NX),ARLFC(NY,NX)
!    2,ARLFT(L,NY,NX),ARLFY,ARLFR,ZL(L,NY,NX),ZL(L-1,NY,NX)
!    3,ARLFV(L,NZ,NY,NX),HVST(NZ,I,NY,NX)
!6544  FORMAT(A8,5I4,20E12.4)
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
      IF(WTSHTA(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        WHVSTT=HVST(NZ,I,NY,NX)*THIN(NZ,I,NY,NX)*0.45/24.0 &
          *AREA(3,NU(NY,NX),NY,NX)*WTSHT(NZ,NY,NX)/WTSHTA(NZ,NY,NX)
      ELSE
        WHVSTT=0._r8
      ENDIF
      IF(IHVST(NZ,I,NY,NX).EQ.6)THEN
        WHVSTT=WHVSTT*TFN3(NZ,NY,NX)
      ENDIF
      CCPOLX=CCPOLP(NZ,NY,NX)/(1.0+CCPOLP(NZ,NY,NX))
      CCPLNX=CCPLNP(NZ,NY,NX)/(1.0+CCPLNP(NZ,NY,NX))
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
      WHVSLX=WHVSTT*EHVST(1,1,NZ,I,NY,NX)
      WHVSLY=AMIN1(WTLF(NZ,NY,NX),WHVSLX)
      WHVSLF=WHVSLY*(1.0-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVSLX-WHVSLY)
      WHVSSX=WHVSTT*EHVST(1,2,NZ,I,NY,NX)
!
!     OTHER NON-FOLIAR GRAZED,REMOVED
!
!     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
!     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
!            petiole,husk,ear,grain,nonstructural C removed
!     WHVXXX=grazing requirement unmet by non-foliar removal
!
      WTSHTT=WTSHE(NZ,NY,NX)+WTHSK(NZ,NY,NX)+WTEAR(NZ,NY,NX)+WTGR(NZ,NY,NX)
      IF(WTSHTT.GT.ZEROP(NZ,NY,NX))THEN
        WHVSHX=WHVSSX*WTSHE(NZ,NY,NX)/WTSHTT+WHVXXX
        WHVSHY=AMIN1(WTSHE(NZ,NY,NX),WHVSHX)
        WHVSHH=WHVSHY*(1.0-CCPOLX)
        WHVSCS=WHVSHY*CCPOLX
        WHVSNS=WHVSHY*CCPLNX
        WHVXXX=AMAX1(0.0,WHVSHX-WHVSHY)
        WHVHSX=WHVSSX*WTHSK(NZ,NY,NX)/WTSHTT+WHVXXX
        WHVHSY=AMIN1(WTHSK(NZ,NY,NX),WHVHSX)
        WHVHSH=WHVHSY
        WHVXXX=AMAX1(0.0,WHVHSX-WHVHSY)
        WHVEAX=WHVSSX*WTEAR(NZ,NY,NX)/WTSHTT+WHVXXX
        WHVEAY=AMIN1(WTEAR(NZ,NY,NX),WHVEAX)
        WHVEAH=WHVEAY
        WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
        WHVGRX=WHVSSX*WTGR(NZ,NY,NX)/WTSHTT+WHVXXX
        WHVGRY=AMIN1(WTGR(NZ,NY,NX),WHVGRX)
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
      WHVSKX=WHVSTT*EHVST(1,3,NZ,I,NY,NX)
!
!     STALK GRAZED, REMOVED
!
!     WTSTK,WTRSV=stalk,reserve C mass
!     WHVST*,WHVRV*=stalk,reserve C removed
!     WHVXXX=grazing requirement unmet by stalk,reserve
!
      WTSTKT=WTSTK(NZ,NY,NX)+WTRSV(NZ,NY,NX)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
        WHVSTX=WHVSKX*WTSTK(NZ,NY,NX)/WTSTKT+WHVXXX
        WHVSTY=AMIN1(WTSTK(NZ,NY,NX),WHVSTX)
        WHVSTH=WHVSTY
        WHVXXX=AMAX1(0.0,WHVSTX-WHVSTY)
        WHVRVX=WHVSKX*WTRSV(NZ,NY,NX)/WTSTKT+WHVXXX
        WHVRVY=AMIN1(WTRSV(NZ,NY,NX),WHVRVX)
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
          WHVSLY=AMIN1(WTLF(NZ,NY,NX)-WHVSLF-WHVSCL,WHVXXX)
          WHVSLF=WHVSLF+WHVSLY*(1.0-CCPOLX)
          WHVSCL=WHVSCL+WHVSLY*CCPOLX
          WHVSNL=WHVSNL+WHVSLY*CCPLNX
          WHVXXX=AMAX1(0.0,WHVXXX-WHVSLY)
          IF(WTSHTT.GT.ZEROP(NZ,NY,NX))THEN
            WHVSHX=WHVXXX*WTSHE(NZ,NY,NX)/WTSHTT
            WHVSHY=AMIN1(WTSHE(NZ,NY,NX),WHVSHX)
            WHVSHH=WHVSHH+WHVSHY*(1.0-CCPOLX)
            WHVSCS=WHVSCS+WHVSHY*CCPOLX
            WHVSNS=WHVSNS+WHVSHY*CCPLNX
            WHVXXX=AMAX1(0.0,WHVXXX-WHVSHY)
            WHVHSX=WHVXXX*WTHSK(NZ,NY,NX)/WTSHTT
            WHVHSY=AMIN1(WTHSK(NZ,NY,NX),WHVHSX)
            WHVHSH=WHVHSH+WHVHSY
            WHVXXX=AMAX1(0.0,WHVXXX-WHVHSY)
            WHVEAX=WHVXXX*WTEAR(NZ,NY,NX)/WTSHTT
            WHVEAY=AMIN1(WTEAR(NZ,NY,NX),WHVEAX)
            WHVEAH=WHVEAH+WHVEAY
            WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
            WHVGRX=WHVXXX*WTGR(NZ,NY,NX)/WTSHTT
            WHVGRY=AMIN1(WTGR(NZ,NY,NX),WHVGRX)
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
      DO 9860 NB=1,NBR(NZ,NY,NX)
        DO  L=1,JC
          DO  K=0,25
            WGLFBL(L,NB,NZ,NY,NX)=0._r8
          enddo
        enddo
9860  CONTINUE
      DO 9870 NB=1,NBR(NZ,NY,NX)
        DO  L=1,JC
          DO  K=0,25
            WGLFBL(L,NB,NZ,NY,NX)=WGLFBL(L,NB,NZ,NY,NX)+WGLFL(L,K,NB,NZ,NY,NX)
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
    DO 9865 L=JC,1,-1
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
        IF(IHVST(NZ,I,NY,NX).NE.3)THEN
          IF(ZL(L,NY,NX).GT.ZL(L-1,NY,NX))THEN
            FHGT=AMAX1(0.0,AMIN1(1.0,1.0-((ZL(L,NY,NX)) &
              -HVST(NZ,I,NY,NX))/(ZL(L,NY,NX)-ZL(L-1,NY,NX))))
          ELSE
            FHGT=1.0_r8
          ENDIF
        ELSE
          FHGT=0._r8
        ENDIF
        IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
          FHVST=AMAX1(0.0,1.0-(1.0-FHGT)*EHVST(1,1,NZ,I,NY,NX))
          FHVSH=FHVST
        ELSE
          FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
          IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
            FHVSH=1.0_r8-(1.0-FHGT)*EHVST(1,1,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
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
      DO 9855 NB=1,NBR(NZ,NY,NX)
        IF((IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6) &
          .AND.WTLF(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
          WHVSBL=WHVSLF*AMAX1(0.0,WGLFBL(L,NB,NZ,NY,NX))/WTLF(NZ,NY,NX)
        ELSE
          WHVSBL=0._r8
        ENDIF
        DO 9845 K=25,0,-1
          IF((IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6).OR.WHVSBL.GT.0.0)THEN
            IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
              IF(WGLFL(L,K,NB,NZ,NY,NX).GT.WHVSBL)THEN
                FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFL(L,K,NB,NZ,NY,NX)-WHVSBL)/WGLFL(L,K,NB,NZ,NY,NX)))
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
            WHVSBL=WHVSBL-(1.0-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)
            WTHTH1=WTHTH1+(1.0-FHVSH)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(1)
            WTHNH1=WTHNH1+(1.0-FHVSH)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(1)
            WTHPH1=WTHPH1+(1.0-FHVSH)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(1)
            WTHTX1=WTHTX1+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(1)
            WTHNX1=WTHNX1+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(1)
            WTHPX1=WTHPX1+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(1)
            WTHTH3=WTHTH3+(1.0-FHVSH)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(0)
            WTHNH3=WTHNH3+(1.0-FHVSH)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(0)
            WTHPH3=WTHPH3+(1.0-FHVSH)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(0)
            WTHTX3=WTHTX3+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(0)
            WTHNX3=WTHNX3+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(0)
            WTHPX3=WTHPX3+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(0)
!
!     REMAINING LEAF C,N,P AND AREA
!
            WGLFL(L,K,NB,NZ,NY,NX)=FHVST*WGLFL(L,K,NB,NZ,NY,NX)
            WGLFLN(L,K,NB,NZ,NY,NX)=FHVST*WGLFLN(L,K,NB,NZ,NY,NX)
            WGLFLP(L,K,NB,NZ,NY,NX)=FHVST*WGLFLP(L,K,NB,NZ,NY,NX)
            ARLFL(L,K,NB,NZ,NY,NX)=FHVST*ARLFL(L,K,NB,NZ,NY,NX)
            IF(K.EQ.1)THEN
              ARSTK(L,NB,NZ,NY,NX)=FHVST*ARSTK(L,NB,NZ,NY,NX)
            ENDIF
          ENDIF
    !     IF(I.EQ.262.AND.K.EQ.5)THEN
!     WRITE(*,6543)'GRAZ',I,J,NZ,NB,K,L,IHVST(NZ,I,NY,NX)
!    2,ZL(L,NY,NX),ZL(L-1,NY,NX),HVST(NZ,I,NY,NX),FHVST,FHVSH
!    5,WGLFBL(L,NB,NZ,NY,NX),WTLF(NZ,NY,NX),CPOOLP(NZ,NY,NX)
!    6,ARLFL(L,K,NB,NZ,NY,NX),WGLF(K,NB,NZ,NY,NX),ARLF(K,NB,NZ,NY,NX)
!    7,HTNODE(K,NB,NZ,NY,NX)
!    7,WTSHTA(NZ,NY,NX),WHVSBL,WHVSTT,WHVSLF,WHVSHH
!    3,WHVHSH,WHVEAH,WHVGRH,WHVSCP,WHVSTH,WHVRVH,WHVXXX
!    4,WTSHTT,WHVSSX,CCPOLX
!6543  FORMAT(A8,7I4,30E12.4)
!     ENDIF
9845      CONTINUE
9855    CONTINUE
        ARLFV(L,NZ,NY,NX)=0._r8
        WGLFV(L,NZ,NY,NX)=0._r8
        ARSTV(L,NZ,NY,NX)=ARSTV(L,NZ,NY,NX)*FHVST
9865  CONTINUE
      DO 9835 NB=1,NBR(NZ,NY,NX)
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
        DO 9825 K=0,25
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
          DO 9815 L=1,JC
            ARLFG=ARLFG+ARLFL(L,K,NB,NZ,NY,NX)
            WGLFG=WGLFG+WGLFL(L,K,NB,NZ,NY,NX)
            WGLFNG=WGLFNG+WGLFLN(L,K,NB,NZ,NY,NX)
            WGLFPG=WGLFPG+WGLFLP(L,K,NB,NZ,NY,NX)
            ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
            WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)+WGLFL(L,K,NB,NZ,NY,NX)
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
          IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
            IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
              .AND.EHVST(1,1,NZ,I,NY,NX).GT.0.0)THEN
              FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(1.0-(1.0-AMAX1(0.0,WGLFG) &
                /WGLF(K,NB,NZ,NY,NX))*EHVST(1,2,NZ,I,NY,NX)/EHVST(1,1,NZ,I,NY,NX))))
              FHVSHK(K)=FHVSTK(K)
          ELSE
            IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
              FHVSTK(K)=1.0_r8-EHVST(1,2,NZ,I,NY,NX)
              FHVSHK(K)=FHVSTK(K)
            ELSE
              FHVSTK(K)=1.0_r8-THIN(NZ,I,NY,NX)
              IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
                FHVSHK(K)=1.0_r8-EHVST(1,2,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
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
        WGLFGY=WGLFGY+WGLF(K,NB,NZ,NY,NX)
        WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)-WGLF(K,NB,NZ,NY,NX)+WGLFG
        WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)-WGLFN(K,NB,NZ,NY,NX)+WGLFNG
        WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)-WGLFP(K,NB,NZ,NY,NX)+WGLFPG
        ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)-ARLF(K,NB,NZ,NY,NX)+ARLFG
        IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
          WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)*ARLFG/ARLF(K,NB,NZ,NY,NX)
        ELSE
          WSLF(K,NB,NZ,NY,NX)=0._r8
        ENDIF
        ARLF(K,NB,NZ,NY,NX)=ARLFG
        WGLF(K,NB,NZ,NY,NX)=WGLFG
        WGLFN(K,NB,NZ,NY,NX)=WGLFNG
        WGLFP(K,NB,NZ,NY,NX)=WGLFPG
        WGLFGX=WGLFGX+WGLF(K,NB,NZ,NY,NX)
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
      IF((IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6) &
        .AND.WTSHE(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        WHVSBS=WHVSHH*WTSHEB(NB,NZ,NY,NX)/WTSHE(NZ,NY,NX)
      ELSE
        WHVSBS=0._r8
      ENDIF
      DO 9805 K=25,0,-1
!112   FORMAT(A8,8I4,12E12.4)
        IF(HTNODE(K,NB,NZ,NY,NX).GT.0.0) &
          HTSTKX=AMAX1(HTSTKX,HTNODE(K,NB,NZ,NY,NX))
!     WRITE(*,112)'VSTG',I,J,NX,NY,NZ,NB,K,IDTHB(NB,NZ,NY,NX)
!    2,VSTG(NB,NZ,NY,NX),FHVSTK(K),HTSTKX,HTNODE(K,NB,NZ,NY,NX)
!    3,HVST(NZ,I,NY,NX)
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
          IF((IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6).OR.WHVSBS.GT.0.0)THEN
            IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
              IF(WGSHE(K,NB,NZ,NY,NX).GT.WHVSBS)THEN
                FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(WGSHE(K,NB,NZ,NY,NX)-WHVSBS)/WGSHE(K,NB,NZ,NY,NX)))
                FHVSHK(K)=FHVSTK(K)
              ELSE
                FHVSTK(K)=0._r8
                FHVSHK(K)=0._r8
              ENDIF
            ENDIF
            WHVSBS=WHVSBS-(1.0-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
            WTHTH2=WTHTH2+(1.0-FHVSHK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(1)
            WTHNH2=WTHNH2+(1.0-FHVSHK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(1)
            WTHPH2=WTHPH2+(1.0-FHVSHK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(1)
            WTHTX2=WTHTX2+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(1)
            WTHNX2=WTHNX2+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(1)
            WTHPX2=WTHPX2+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(1)
            WTHTH3=WTHTH3+(1.0-FHVSHK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(0)
            WTHNH3=WTHNH3+(1.0-FHVSHK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(0)
            WTHPH3=WTHPH3+(1.0-FHVSHK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(0)
            WTHTX3=WTHTX3+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(0)
            WTHNX3=WTHNX3+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(0)
            WTHPX3=WTHPX3+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(0)
!
!     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
!
!     WGSHE=petiole node C mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     HTSHE=node petiole height
!     WSSHE=petiole protein mass
!
            WGSHGY=WGSHGY+WGSHE(K,NB,NZ,NY,NX)
            WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX) &
              -(1.0-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
            WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX) &
              -(1.0-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)
            WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX) &
              -(1.0-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)
            WGSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHE(K,NB,NZ,NY,NX)
            WSSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WSSHE(K,NB,NZ,NY,NX)
            WGSHN(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHN(K,NB,NZ,NY,NX)
            WGSHP(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHP(K,NB,NZ,NY,NX)
            WSSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WSSHE(K,NB,NZ,NY,NX)
            IF(IHVST(NZ,I,NY,NX).LE.2 &
              .AND.HTSHE(K,NB,NZ,NY,NX).GT.0.0)THEN
              FHGT=AMAX1(0.0,AMIN1(1.0,(HTNODE(K,NB,NZ,NY,NX) &
                +HTSHE(K,NB,NZ,NY,NX)-HVST(NZ,I,NY,NX))/HTSHE(K,NB,NZ,NY,NX)))
              HTSHE(K,NB,NZ,NY,NX)=(1.0-FHGT)*HTSHE(K,NB,NZ,NY,NX)
            ELSE
              HTSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*HTSHE(K,NB,NZ,NY,NX)
            ENDIF
            WGSHGX=WGSHGX+WGSHE(K,NB,NZ,NY,NX)
!     IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
!     IF(HTNODE(K,NB,NZ,NY,NX).GT.HVST(NZ,I,NY,NX)
!    2.OR.IHVST(NZ,I,NY,NX).EQ.3)THEN
!     IF(test_aeqb(FHVSTK(K),0._r8).AND.K.GT.0)THEN
!     IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
!     VSTG(NB,NZ,NY,NX)=AMAX1(0.0,VSTG(NB,NZ,NY,NX)-1.0)
!     ELSE
!     VSTG(NB,NZ,NY,NX)=AMAX1(0.0,VSTG(NB,NZ,NY,NX)-0.04)
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
        CPOOLX=AMAX1(0.0,CPOOL(NB,NZ,NY,NX))
        ZPOOLX=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
        PPOOLX=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
        CPOLNX=AMAX1(0.0,CPOLNB(NB,NZ,NY,NX))
        ZPOLNX=AMAX1(0.0,ZPOLNB(NB,NZ,NY,NX))
        PPOLNX=AMAX1(0.0,PPOLNB(NB,NZ,NY,NX))
        IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
          IF(WGLFGY+WGSHGY.GT.ZEROP(NZ,NY,NX))THEN
            FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFGX+WGSHGX)/(WGLFGY+WGSHGY)))
            CPOOLG=CPOOLX*FHVST
            ZPOOLG=ZPOOLX*FHVST
            PPOOLG=PPOOLX*FHVST
            CPOLNG=CPOLNX*FHVST
            ZPOLNG=ZPOLNX*FHVST
            PPOLNG=PPOLNX*FHVST
            WTNDG=WTNDB(NB,NZ,NY,NX)*FHVST
            WTNDNG=WTNDBN(NB,NZ,NY,NX)*FHVST
            WTNDPG=WTNDBP(NB,NZ,NY,NX)*FHVST
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
          IF(WTLS(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
            WTLSBX=AMAX1(0.0,WTLSB(NB,NZ,NY,NX))
            IF(CPOOL(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
              WHVSCX=AMAX1(0.0,WHVSCP)*WTLSBX/WTLS(NZ,NY,NX)
              CPOOLG=AMAX1(0.0,CPOOLX-WHVSCX)
              ZPOOLG=AMAX1(0.0,ZPOOLX-WHVSCX*ZPOOLX/CPOOL(NB,NZ,NY,NX))
              PPOOLG=AMAX1(0.0,PPOOLX-WHVSCX*PPOOLX/CPOOL(NB,NZ,NY,NX))
            ELSE
              CPOOLG=0._r8
              ZPOOLG=0._r8
              PPOOLG=0._r8
            ENDIF
            IF(CPOLNB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
              WHVSNX=AMAX1(0.0,WHVSNP)*WTLSBX/WTLS(NZ,NY,NX)
              CPOLNG=AMAX1(0.0,CPOLNX-WHVSNX)
              ZPOLNG=AMAX1(0.0,ZPOLNX-WHVSNX*ZPOLNX/CPOLNB(NB,NZ,NY,NX))
              PPOLNG=AMAX1(0.0,PPOLNX-WHVSNX*PPOLNX/CPOLNB(NB,NZ,NY,NX))
              WTNDG=WTNDB(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
              WTNDNG=WTNDBN(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
              WTNDPG=WTNDBP(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
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
        WTHTH0=WTHTH0+WTNDB(NB,NZ,NY,NX)-WTNDG
        WTHNH0=WTHNH0+WTNDBN(NB,NZ,NY,NX)-WTNDNG
        WTHPH0=WTHPH0+WTNDBP(NB,NZ,NY,NX)-WTNDPG
!
!     REMAINING NON-STRUCTURAL C, N, P
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
        CPOOL(NB,NZ,NY,NX)=CPOOLG
        ZPOOL(NB,NZ,NY,NX)=ZPOOLG
        PPOOL(NB,NZ,NY,NX)=PPOOLG
        CPOLNB(NB,NZ,NY,NX)=CPOLNG
        ZPOLNB(NB,NZ,NY,NX)=ZPOLNG
        PPOLNB(NB,NZ,NY,NX)=PPOLNG
        WTNDB(NB,NZ,NY,NX)=WTNDG
        WTNDBN(NB,NZ,NY,NX)=WTNDNG
        WTNDBP(NB,NZ,NY,NX)=WTNDPG
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
        IF(ICTYP(NZ,NY,NX).EQ.4.AND.CPOOLX.GT.ZEROP(NZ,NY,NX))THEN
          FHVST4=CPOOLG/CPOOLX
          DO 9810 K=1,25
            WTHTH0=WTHTH0+(1.0-FHVST4)*CPOOL3(K,NB,NZ,NY,NX)
            WTHTH0=WTHTH0+(1.0-FHVST4)*CPOOL4(K,NB,NZ,NY,NX)
            WTHTH0=WTHTH0+(1.0-FHVST4)*CO2B(K,NB,NZ,NY,NX)
            WTHTH0=WTHTH0+(1.0-FHVST4)*HCOB(K,NB,NZ,NY,NX)
            CPOOL3(K,NB,NZ,NY,NX)=FHVST4*CPOOL3(K,NB,NZ,NY,NX)
            CPOOL4(K,NB,NZ,NY,NX)=FHVST4*CPOOL4(K,NB,NZ,NY,NX)
            CO2B(K,NB,NZ,NY,NX)=FHVST4*CO2B(K,NB,NZ,NY,NX)
            HCOB(K,NB,NZ,NY,NX)=FHVST4*HCOB(K,NB,NZ,NY,NX)
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
        IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
          IF(HTSTKX.GT.ZERO)THEN
            IF(IHVST(NZ,I,NY,NX).NE.3)THEN
              FHGT=AMAX1(0.0,AMIN1(1.0,HVST(NZ,I,NY,NX)/HTSTKX))
            ELSE
              FHGT=0._r8
            ENDIF
            IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
              FHVST=AMAX1(0.0,1.0-(1.0-FHGT)*EHVST(1,3,NZ,I,NY,NX))
              FHVSH=FHVST
            ELSE
              FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
              IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
                FHVSH=1.0_r8-(1.0-FHGT)*EHVST(1,3,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
              ELSE
                FHVSH=FHVST
              ENDIF
            ENDIF
          ELSE
            FHVST=1.0_r8
            FHVSH=1.0_r8
          ENDIF
        ELSE
          IF(WTSTK(NZ,NY,NX).GT.ZEROL(NZ,NY,NX))THEN
            FHVST=AMAX1(0.0,AMIN1(1.0,1.0-WHVSTH/WTSTK(NZ,NY,NX)))
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
        WTHTH3=WTHTH3+(1.0-FHVSH)*WTSTKB(NB,NZ,NY,NX)
        WTHNH3=WTHNH3+(1.0-FHVSH)*WTSTBN(NB,NZ,NY,NX)
        WTHPH3=WTHPH3+(1.0-FHVSH)*WTSTBP(NB,NZ,NY,NX)
        WTHTX3=WTHTX3+(FHVSH-FHVST)*WTSTKB(NB,NZ,NY,NX)
        WTHNX3=WTHNX3+(FHVSH-FHVST)*WTSTBN(NB,NZ,NY,NX)
        WTHPX3=WTHPX3+(FHVSH-FHVST)*WTSTBP(NB,NZ,NY,NX)
!
!     REMAINING STALK C,N,P
!
!     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
!
        WTSTKB(NB,NZ,NY,NX)=FHVST*WTSTKB(NB,NZ,NY,NX)
        WTSTBN(NB,NZ,NY,NX)=FHVST*WTSTBN(NB,NZ,NY,NX)
        WTSTBP(NB,NZ,NY,NX)=FHVST*WTSTBP(NB,NZ,NY,NX)
        WVSTKB(NB,NZ,NY,NX)=FHVST*WVSTKB(NB,NZ,NY,NX)
        WTSTXB(NB,NZ,NY,NX)=FHVST*WTSTXB(NB,NZ,NY,NX)
        WTSTXN(NB,NZ,NY,NX)=FHVST*WTSTXN(NB,NZ,NY,NX)
        WTSTXP(NB,NZ,NY,NX)=FHVST*WTSTXP(NB,NZ,NY,NX)
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
        DO 9820 K=25,0,-1
          IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
            IF(HTNODX(K,NB,NZ,NY,NX).GT.ZERO)THEN
              IF(IHVST(NZ,I,NY,NX).NE.3)THEN
                FHGTK=AMAX1(0.0,AMIN1(1.0,(HTNODE(K,NB,NZ,NY,NX) &
                  -HVST(NZ,I,NY,NX))/HTNODX(K,NB,NZ,NY,NX)))
              ELSE
                FHGTK=0._r8
              ENDIF
              IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
                FHVSTS=AMAX1(0.0,1.0-FHGTK*EHVST(1,3,NZ,I,NY,NX))
              ELSE
                FHVSTS=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
              ENDIF
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ELSE
            IF(WTSTK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
              FHVSTS=AMAX1(0.0,AMIN1(1.0,1.0-WHVSTH/WTSTK(NZ,NY,NX)))
            ELSE
              FHVSTS=1.0_r8
            ENDIF
          ENDIF
          WGNODE(K,NB,NZ,NY,NX)=FHVSTS*WGNODE(K,NB,NZ,NY,NX)
          WGNODN(K,NB,NZ,NY,NX)=FHVSTS*WGNODN(K,NB,NZ,NY,NX)
          WGNODP(K,NB,NZ,NY,NX)=FHVSTS*WGNODP(K,NB,NZ,NY,NX)
          IF(IHVST(NZ,I,NY,NX).LE.2.AND.test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
            HTNODX(K,NB,NZ,NY,NX)=FHVSTS*HTNODX(K,NB,NZ,NY,NX)
            HTNODE(K,NB,NZ,NY,NX)=AMIN1(HTNODE(K,NB,NZ,NY,NX),HVST(NZ,I,NY,NX))
          ENDIF
!     IF(NZ.EQ.2)THEN
!     WRITE(*,4811)'STK2',I,J,NX,NY,NZ,NB,K,IHVST(NZ,I,NY,NX)
!    2,HTNODX(K,NB,NZ,NY,NX),HTNODE(K,NB,NZ,NY,NX)
!    3,HVST(NZ,I,NY,NX),FHGTK,FHVSTS,ARLF(K,NB,NZ,NY,NX)
!    4,EHVST(1,3,NZ,I,NY,NX),THIN(NZ,I,NY,NX)
!4811  FORMAT(A8,8I4,12E12.4)
!     ENDIF
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
        IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
          IF(WTSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            FHVST=FHVST
            FHVSH=FHVSH
          ELSE
            FHVST=0._r8
            FHVSH=0._r8
          ENDIF
        ELSE
          IF(WTRSV(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            FHVST=AMAX1(0.0,AMIN1(1.0,1.0-WHVRVH/WTRSV(NZ,NY,NX)))
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
        WTHTH3=WTHTH3+(1.0-FHVSH)*WTRSVB(NB,NZ,NY,NX)
        WTHNH3=WTHNH3+(1.0-FHVSH)*WTRSBN(NB,NZ,NY,NX)
        WTHPH3=WTHPH3+(1.0-FHVSH)*WTRSBP(NB,NZ,NY,NX)
        WTHTX3=WTHTX3+(FHVSH-FHVST)*WTRSVB(NB,NZ,NY,NX)
        WTHNX3=WTHNX3+(FHVSH-FHVST)*WTRSBN(NB,NZ,NY,NX)
        WTHPX3=WTHPX3+(FHVSH-FHVST)*WTRSBP(NB,NZ,NY,NX)
!
!     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
!
        WTRSVB(NB,NZ,NY,NX)=FHVST*WTRSVB(NB,NZ,NY,NX)
        WTRSBN(NB,NZ,NY,NX)=FHVST*WTRSBN(NB,NZ,NY,NX)
        WTRSBP(NB,NZ,NY,NX)=FHVST*WTRSBP(NB,NZ,NY,NX)
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
        IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
          IF(HVST(NZ,I,NY,NX).LT.HTSTKX &
            .OR.IHVST(NZ,I,NY,NX).EQ.1 &
            .OR.IHVST(NZ,I,NY,NX).EQ.3)THEN
            IF(test_aeqb(THIN(NZ,I,NY,NX),0._r8))THEN
              FHVSTG=1.0_r8-EHVST(1,2,NZ,I,NY,NX)
              FHVSHG=FHVSTG
            ELSE
              FHVSTG=1.0_r8-THIN(NZ,I,NY,NX)
              FHVSHG=1.0_r8-EHVST(1,2,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
            ENDIF
          ELSE
            FHVSTG=1.0_r8-THIN(NZ,I,NY,NX)
            FHVSHG=FHVSTG
          ENDIF
          FHVSTH=FHVSTG
          FHVSTE=FHVSTG
          FHVSHH=FHVSHG
          FHVSHE=FHVSHG
        ELSE
          IF(WTHSK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            FHVSTH=AMAX1(0.0,AMIN1(1.0,1.0-WHVHSH/WTHSK(NZ,NY,NX)))
            FHVSHH=FHVSTH
          ELSE
            FHVSTH=1.0_r8
            FHVSHH=1.0_r8
          ENDIF
          IF(WTEAR(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            FHVSTE=AMAX1(0.0,AMIN1(1.0,1.0-WHVEAH/WTEAR(NZ,NY,NX)))
            FHVSHE=FHVSTE
          ELSE
            FHVSTE=1.0_r8
            FHVSHE=1.0_r8
          ENDIF
          IF(WTGR(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
            FHVSTG=AMAX1(0.0,AMIN1(1.0,1.0-WHVGRH/WTGR(NZ,NY,NX)))
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
        WTHTH2=WTHTH2+(1.0-FHVSHH)*WTHSKB(NB,NZ,NY,NX)+(1.0-FHVSHE) &
          *WTEARB(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRB(NB,NZ,NY,NX)
        WTHNH2=WTHNH2+(1.0-FHVSHH)*WTHSBN(NB,NZ,NY,NX)+(1.0-FHVSHE) &
          *WTEABN(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRBN(NB,NZ,NY,NX)
        WTHPH2=WTHPH2+(1.0-FHVSHH)*WTHSBP(NB,NZ,NY,NX)+(1.0-FHVSHE) &
          *WTEABP(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRBP(NB,NZ,NY,NX)
        WTHTX2=WTHTX2+(FHVSHH-FHVSTH)*WTHSKB(NB,NZ,NY,NX)+(FHVSHE-FHVSTE) &
          *WTEARB(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRB(NB,NZ,NY,NX)
        WTHNX2=WTHNX2+(FHVSHH-FHVSTH)*WTHSBN(NB,NZ,NY,NX)+(FHVSHE-FHVSTE) &
          *WTEABN(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRBN(NB,NZ,NY,NX)
        WTHPX2=WTHPX2+(FHVSHH-FHVSTH)*WTHSBP(NB,NZ,NY,NX)+(FHVSHE-FHVSTE) &
          *WTEABP(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRBP(NB,NZ,NY,NX)
        WTHTG=WTHTG+(1.0-FHVSTG)*WTGRB(NB,NZ,NY,NX)
        WTHNG=WTHNG+(1.0-FHVSTG)*WTGRBN(NB,NZ,NY,NX)
        WTHPG=WTHPG+(1.0-FHVSTG)*WTGRBP(NB,NZ,NY,NX)
!
!     REMAINING REPRODUCTIVE C,N,P
!
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!
        WTHSKB(NB,NZ,NY,NX)=FHVSTH*WTHSKB(NB,NZ,NY,NX)
        WTEARB(NB,NZ,NY,NX)=FHVSTE*WTEARB(NB,NZ,NY,NX)
        WTGRB(NB,NZ,NY,NX)=FHVSTG*WTGRB(NB,NZ,NY,NX)
        WTHSBN(NB,NZ,NY,NX)=FHVSTH*WTHSBN(NB,NZ,NY,NX)
        WTEABN(NB,NZ,NY,NX)=FHVSTE*WTEABN(NB,NZ,NY,NX)
        WTGRBN(NB,NZ,NY,NX)=FHVSTG*WTGRBN(NB,NZ,NY,NX)
        WTHSBP(NB,NZ,NY,NX)=FHVSTH*WTHSBP(NB,NZ,NY,NX)
        WTEABP(NB,NZ,NY,NX)=FHVSTE*WTEABP(NB,NZ,NY,NX)
        WTGRBP(NB,NZ,NY,NX)=FHVSTG*WTGRBP(NB,NZ,NY,NX)
        GRNXB(NB,NZ,NY,NX)=FHVSTG*GRNXB(NB,NZ,NY,NX)
        GRNOB(NB,NZ,NY,NX)=FHVSTG*GRNOB(NB,NZ,NY,NX)
        GRWTB(NB,NZ,NY,NX)=FHVSTG*GRWTB(NB,NZ,NY,NX)
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
        CPOOLK(NB,NZ,NY,NX)=0._r8
        DO 1325 K=1,25
          CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX) &
            +CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX) &
            +CO2B(K,NB,NZ,NY,NX)+HCOB(K,NB,NZ,NY,NX)
1325    CONTINUE
        WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX) &
          +WTSHEB(NB,NZ,NY,NX))
        WTSHTB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX) &
          +WTSHEB(NB,NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX) &
          +WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX) &
          +CPOOL(NB,NZ,NY,NX)+CPOOLK(NB,NZ,NY,NX))
        WTSHTN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX) &
          +WTSHBN(NB,NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX) &
          +WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX) &
          +ZPOOL(NB,NZ,NY,NX))
        WTSHTP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX) &
          +WTSHBP(NB,NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX) &
          +WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX) &
          +PPOOL(NB,NZ,NY,NX))
        VOLWPX=VOLWP(NZ,NY,NX)
        WVPLT=AMAX1(0.0_r8,WTLS(NZ,NY,NX)+WVSTK(NZ,NY,NX))
        APSILT=ABS(PSILT(NZ,NY,NX))
        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
        VOLWP(NZ,NY,NX)=1.0E-06_r8*WVPLT/FDM
        VOLWOU=VOLWOU+VOLWPX-VOLWP(NZ,NY,NX)
        UVOLO(NY,NX)=UVOLO(NY,NX)+VOLWPX-VOLWP(NZ,NY,NX)
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
        IF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1) &
          .AND.(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6) &
          .AND.ZC(NZ,NY,NX).GT.HVST(NZ,I,NY,NX))THEN
          IF((IWTYP(NZ,NY,NX).NE.0.AND.VRNF(NB,NZ,NY,NX) &
            .LE.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX)) &
            .OR.(IWTYP(NZ,NY,NX).EQ.0 &
            .AND.IDAY(1,NB,NZ,NY,NX).NE.0))THEN
            GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
            PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
            PSTGF(NB,NZ,NY,NX)=0._r8
            VSTGX(NB,NZ,NY,NX)=0._r8
            TGSTGI(NB,NZ,NY,NX)=0._r8
            TGSTGF(NB,NZ,NY,NX)=0._r8
            FLG4(NB,NZ,NY,NX)=0._r8
            IDAY(1,NB,NZ,NY,NX)=I
            DO 3005 M=2,10
              IDAY(M,NB,NZ,NY,NX)=0
3005        CONTINUE
            IFLGA(NB,NZ,NY,NX)=0
            IF(NB.EQ.NB1(NZ,NY,NX))THEN
              DO 3010 NBX=1,NBR(NZ,NY,NX)
                IF(NBX.NE.NB1(NZ,NY,NX))THEN
                  GROUP(NBX,NZ,NY,NX)=GROUPI(NZ,NY,NX)
                  PSTGI(NBX,NZ,NY,NX)=PSTG(NBX,NZ,NY,NX)
                  PSTGF(NBX,NZ,NY,NX)=0._r8
                  VSTGX(NBX,NZ,NY,NX)=0._r8
                  TGSTGI(NBX,NZ,NY,NX)=0._r8
                  TGSTGF(NBX,NZ,NY,NX)=0._r8
                  FLG4(NBX,NZ,NY,NX)=0._r8
                  IDAY(1,NBX,NZ,NY,NX)=I
                  DO 3015 M=2,10
                    IDAY(M,NBX,NZ,NY,NX)=0
3015              CONTINUE
                  IFLGA(NBX,NZ,NY,NX)=0
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
      IF(JHVST(NZ,I,NY,NX).NE.0)IDTHB(NB,NZ,NY,NX)=1
      IF(PP(NZ,NY,NX).LE.0.0)IDTHB(NB,NZ,NY,NX)=1
9835  CONTINUE
      WTLS(NZ,NY,NX)=0._r8
      WTSTK(NZ,NY,NX)=0._r8
      WVSTK(NZ,NY,NX)=0._r8
      ARSTP(NZ,NY,NX)=0._r8
      DO 9840 NB=1,NBR(NZ,NY,NX)
        WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
        WTSTK(NZ,NY,NX)=WTSTK(NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)
        WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
        DO 9830 L=1,JC
          ARSTP(NZ,NY,NX)=ARSTP(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
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
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
        XHVST=1.0_r8-THIN(NZ,I,NY,NX)
        DO 3985 N=1,MY(NZ,NY,NX)
          DO 3980 L=NU(NY,NX),NJ(NY,NX)
            IF(IHVST(NZ,I,NY,NX).NE.5)THEN
              XHVST=1.0_r8-THIN(NZ,I,NY,NX)
              XHVSN=XHVST
              XHVSP=XHVST
              FFIRE=0._r8
              FFIRN=0._r8
              FFIRP=0._r8
            ELSE
              IF(THETW(L,NY,NX).GT.FVLWB.OR.CORGC(L,NY,NX).LE.FORGC &
                .OR.ITILL(I,NY,NX).NE.22)THEN
                XHVST=1.0_r8
                XHVSN=XHVST
                XHVSP=XHVST
                FFIRE=0._r8
                FFIRN=0._r8
                FFIRP=0._r8
              ELSE
                XHVST=1.0_r8-DCORP(I,NY,NX)*EHVST(1,3,NZ,I,NY,NX) &
                  *AMIN1(1.0,(CORGC(L,NY,NX)-FORGC)/(0.55E+06-FORGC))
                XHVSN=XHVST
                XHVSP=XHVST
                FFIRE=EHVST(2,3,NZ,I,NY,NX)
                FFIRN=FFIRE*EFIRE(1,IHVST(NZ,I,NY,NX))
                FFIRP=FFIRE*EFIRE(2,IHVST(NZ,I,NY,NX))
              ENDIF
            ENDIF
            DO 3385 M=1,4
              FHVST=(1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
              FHVSN=(1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*ZPOOLR(N,L,NZ,NY,NX)
              FHVSP=(1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*PPOOLR(N,L,NZ,NY,NX)
              CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRE)*FHVST
              ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRN)*FHVSN
              PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRP)*FHVSP
              VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
              VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*FFIRE*FHVST
              VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST*2.667
              VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-FFIRN*FHVSN
              VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
              VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-FFIRP*FHVSP
              CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
              TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*FFIRE*FHVST
              DO NR=1,NRT(NZ,NY,NX)
                FHVST=(1.0-XHVST)*CFOPC(5,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
                  +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0)
                FHVSN=(1.0-XHVSN)*CFOPN(5,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
                  +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0)
                FHVSP=(1.0-XHVSP)*CFOPP(5,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
                  +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0)
                CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRE)*FHVST
                ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRN)*FHVSN
                PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRP)*FHVSP
                VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
                VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*FFIRE*FHVST
                VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST*2.667
                VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-FFIRN*FHVSN
                VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
                VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-FFIRP*FHVSP
                CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
                TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*FFIRE*FHVST
                FHVST=(1.0-XHVST)*CFOPC(4,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX) &
                  +WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1)
                FHVSN=(1.0-XHVSN)*CFOPN(4,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX) &
                  +WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1)
                FHVSP=(1.0-XHVSP)*CFOPP(4,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX) &
                  +WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1)
                CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRE)*FHVST
                ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRN)*FHVSN
                PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-FFIRP)*FHVSP
                VCO2F(NZ,NY,NX)=VCO2F(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
                VCH4F(NZ,NY,NX)=VCH4F(NZ,NY,NX)-FCH4F*FFIRE*FHVST
                VOXYF(NZ,NY,NX)=VOXYF(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST*2.667
                VNH3F(NZ,NY,NX)=VNH3F(NZ,NY,NX)-FFIRN*FHVSN
                VN2OF(NZ,NY,NX)=VN2OF(NZ,NY,NX)-0.0
                VPO4F(NZ,NY,NX)=VPO4F(NZ,NY,NX)-FFIRP*FHVSP
                CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-(1.0-FCH4F)*FFIRE*FHVST
                TNBP(NY,NX)=TNBP(NY,NX)-FCH4F*FFIRE*FHVST
              enddo
3385        CONTINUE
!
!     RELEASE ROOT GAS CONTENTS DURING HARVESTING
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!
            RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(1.0-XHVST) &
              *(CO2A(N,L,NZ,NY,NX)+CO2P(N,L,NZ,NY,NX))
            ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(1.0-XHVST) &
              *(OXYA(N,L,NZ,NY,NX)+OXYP(N,L,NZ,NY,NX))
            RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(1.0-XHVST) &
              *(CH4A(N,L,NZ,NY,NX)+CH4P(N,L,NZ,NY,NX))
            RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(1.0-XHVST) &
              *(Z2OA(N,L,NZ,NY,NX)+Z2OP(N,L,NZ,NY,NX))
            RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(1.0-XHVST) &
              *(ZH3A(N,L,NZ,NY,NX)+ZH3P(N,L,NZ,NY,NX))
            RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(1.0-XHVST) &
              *(H2GA(N,L,NZ,NY,NX)+H2GP(N,L,NZ,NY,NX))
            CO2A(N,L,NZ,NY,NX)=XHVST*CO2A(N,L,NZ,NY,NX)
            OXYA(N,L,NZ,NY,NX)=XHVST*OXYA(N,L,NZ,NY,NX)
            CH4A(N,L,NZ,NY,NX)=XHVST*CH4A(N,L,NZ,NY,NX)
            Z2OA(N,L,NZ,NY,NX)=XHVST*Z2OA(N,L,NZ,NY,NX)
            ZH3A(N,L,NZ,NY,NX)=XHVST*ZH3A(N,L,NZ,NY,NX)
            H2GA(N,L,NZ,NY,NX)=XHVST*H2GA(N,L,NZ,NY,NX)
            CO2P(N,L,NZ,NY,NX)=XHVST*CO2P(N,L,NZ,NY,NX)
            OXYP(N,L,NZ,NY,NX)=XHVST*OXYP(N,L,NZ,NY,NX)
            CH4P(N,L,NZ,NY,NX)=XHVST*CH4P(N,L,NZ,NY,NX)
            Z2OP(N,L,NZ,NY,NX)=XHVST*Z2OP(N,L,NZ,NY,NX)
            ZH3P(N,L,NZ,NY,NX)=XHVST*ZH3P(N,L,NZ,NY,NX)
            H2GP(N,L,NZ,NY,NX)=XHVST*H2GP(N,L,NZ,NY,NX)
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
            DO 3960 NR=1,NRT(NZ,NY,NX)
              WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*XHVST
              WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*XHVST
              WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*XHVSN
              WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*XHVSN
              WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*XHVSP
              WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*XHVSP
              RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*XHVST
              RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*XHVST
              RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*XHVST
              RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*XHVST
              RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*XHVST
              RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*XHVST
3960        CONTINUE
            CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*XHVST
            ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*XHVSN
            PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*XHVSP
            WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*XHVST
            WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*XHVST
            WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*XHVST
            RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*XHVST
            RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*XHVST
            RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*XHVST
            RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*XHVST
            RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*XHVST
            RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*XHVST
            RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*XHVST
            RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)*XHVST
            RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)*XHVST
            RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)*XHVST
!
!     NODULE LITTERFALL AND STATE VARIABLES DURING HARVESTING
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
            IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
              DO 3395 M=1,4
                CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST) &
                  *(CFOPC(4,M,NZ,NY,NX)*WTNDL(L,NZ,NY,NX) &
                  +CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX))
                ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVSN) &
                  *(CFOPN(4,M,NZ,NY,NX)*WTNDLN(L,NZ,NY,NX) &
                  +CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX))
                PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVSP) &
                  *(CFOPP(4,M,NZ,NY,NX)*WTNDLP(L,NZ,NY,NX) &
                  +CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX))
3395          CONTINUE
              WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*XHVST
              WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*XHVSN
              WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*XHVSP
              CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*XHVST
              ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*XHVSN
              PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*XHVSP
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
        IF(ISTYP(NZ,NY,NX).NE.0)THEN
          DO 3400 M=1,4
            CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
              +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(0)
            ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
              +((1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(0)
            PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX) &
              +((1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(0)
            CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
              +((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(1)
            ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
              +((1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(1)
            PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX) &
              +((1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(1)
3400      CONTINUE
          WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*XHVST
          WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*XHVSN
          WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*XHVSP
        ENDIF
    ENDIF
  ENDIF
  end subroutine RemoveBiomByHarvest
!------------------------------------------------------------------------------------------

end module PlantDisturbMod
