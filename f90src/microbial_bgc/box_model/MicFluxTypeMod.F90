module MicFluxTypeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL     
  use EcoSiMParDataMod, only : micpar
  use EcoSIMSolverPar, only : NPH
  use abortutils, only : destroy
  use TracerIDMod
implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  
  type, public :: micfluxtype

  real(r8) :: RCO2O
  real(r8) :: RCH4O
  real(r8) :: RH2GO
  real(r8) :: RO2UptkMicb
  real(r8) :: RN2G
  real(r8) :: RN2O
  real(r8) :: RNH4MicbTransf_vr
  real(r8) :: RNO3MicbTransf_vr
  real(r8) :: RNO2MicbTransf_vr
  real(r8) :: RH2PO4MicbTransf_vr
  real(r8) :: RH1PO4MicbTransf_vr
  real(r8) :: XNH4B
  real(r8) :: XNO3B
  real(r8) :: XNO2B
  real(r8) :: XH2BS
  real(r8) :: XH1BS
  real(r8) :: XN2GS
  real(r8) :: RVMXC
  real(r8) :: RVMBC
  real(r8) :: NetNH4Mineralize_col
  real(r8) :: NetPO4Mineralize_col
  real(r8), allocatable :: RDOM_micb_flx(:,:)
  real(r8), allocatable :: RO2DmndAutort(:)
  real(r8), allocatable :: RNH3OxidAutor(:)
  real(r8), allocatable :: RNH3OxidAutorBand(:)
  real(r8), allocatable :: RNO2OxidAutor(:)
  real(r8), allocatable :: RNO2OxidAutorBand(:)
  real(r8), allocatable :: RO2DmndHetert(:,:)
  real(r8), allocatable :: ROQCS(:,:)
  real(r8), allocatable :: ROQAS(:,:)
  real(r8), allocatable :: RVMX3(:,:)
  real(r8), allocatable :: RVMB3(:,:)
  real(r8), allocatable :: RVMX2(:,:)
  real(r8), allocatable :: RVMB2(:,:)
  real(r8), allocatable :: RVMX1(:,:)
  real(r8), allocatable :: ROXSK(:)
  real(r8), allocatable :: RVMX4(:,:)
  real(r8), allocatable :: RVMB4(:,:)
  real(r8), allocatable :: RINHO(:,:)
  real(r8), allocatable :: RINHB(:,:)
  real(r8), allocatable :: RINOO(:,:)
  real(r8), allocatable :: RINOB(:,:)
  real(r8), allocatable :: RIPOO(:,:)
  real(r8), allocatable :: RIPBO(:,:)
  real(r8), allocatable :: RIPO1(:,:)
  real(r8), allocatable :: RIPB1(:,:)
  real(r8), allocatable :: RINHOR(:,:)
  real(r8), allocatable :: RINOOR(:,:)
  real(r8), allocatable :: RIPOOR(:,:)
  real(r8), allocatable :: RIPO1R(:,:)
  real(r8), allocatable :: RNH4UptkSoilAutor(:)
  real(r8), allocatable :: RNH4UptkBandAutor(:)
  real(r8), allocatable :: RNO3UptkSoilAutor(:)
  real(r8), allocatable :: RNO3UptkBandAutor(:)
  real(r8), allocatable :: RH2PO4UptkSoilAutor(:)
  real(r8), allocatable :: RH2PO4UptkBandAutor(:)
  real(r8), allocatable :: RH1PO4UptkSoilAutor(:)
  real(r8), allocatable :: RH1PO4UptkBandAutor(:)
  real(r8), allocatable :: RNH4UptkLitrAutor(:)
  real(r8), allocatable :: RNO3UptkLitrAutor(:)
  real(r8), allocatable :: RH2PO4UptkLitrAutor(:)
  real(r8), allocatable :: RH1PO4UptkLitrAutor(:)
  contains
   procedure, public :: Init
   procedure, public :: Destroy=> Destruct
   procedure, public :: ZeroOut
  end type micfluxtype

  contains

  subroutine Init(this)

  implicit none
  class(micfluxtype) :: this
  integer :: jcplx,JG,NumMicbFunGroups
  integer :: NumMicrbHetetrophCmplx, NumMicrobAutrophCmplx

  jcplx=micpar%jcplx
  JG=micpar%jguilds
  NumMicbFunGroups=micpar%NumMicbFunGroups
  NumMicrbHetetrophCmplx=micpar%NumMicrbHetetrophCmplx
  NumMicrobAutrophCmplx=micpar%NumMicrobAutrophCmplx

  allocate(this%ROXSK(NPH));this%ROXSK = spval
  allocate(this%RDOM_micb_flx(idom_beg:idom_end,1:jcplx));this%RDOM_micb_flx=spval
  allocate(this%RO2DmndHetert(NumMicrbHetetrophCmplx,1:jcplx));this%RO2DmndHetert=spval
  allocate(this%ROQCS(NumMicrbHetetrophCmplx,1:jcplx));this%ROQCS=spval
  allocate(this%ROQAS(NumMicrbHetetrophCmplx,1:jcplx));this%ROQAS=spval
  allocate(this%RVMX3(NumMicrbHetetrophCmplx,1:jcplx));this%RVMX3=spval
  allocate(this%RVMB3(NumMicrbHetetrophCmplx,1:jcplx));this%RVMB3=spval
  allocate(this%RVMX2(NumMicrbHetetrophCmplx,1:jcplx));this%RVMX2=spval
  allocate(this%RVMB2(NumMicrbHetetrophCmplx,1:jcplx));this%RVMB2=spval
  allocate(this%RVMX1(NumMicrbHetetrophCmplx,1:jcplx));this%RVMX1=spval
  allocate(this%RVMX4(NumMicrbHetetrophCmplx,1:jcplx));this%RVMX4=spval
  allocate(this%RVMB4(NumMicrbHetetrophCmplx,1:jcplx));this%RVMB4=spval
  allocate(this%RINHO(NumMicrbHetetrophCmplx,1:jcplx));this%RINHO=spval
  allocate(this%RINHB(NumMicrbHetetrophCmplx,1:jcplx));this%RINHB=spval
  allocate(this%RINOO(NumMicrbHetetrophCmplx,1:jcplx));this%RINOO=spval
  allocate(this%RINOB(NumMicrbHetetrophCmplx,1:jcplx));this%RINOB=spval
  allocate(this%RIPOO(NumMicrbHetetrophCmplx,1:jcplx));this%RIPOO=spval
  allocate(this%RIPBO(NumMicrbHetetrophCmplx,1:jcplx));this%RIPBO=spval
  allocate(this%RIPO1(NumMicrbHetetrophCmplx,1:jcplx));this%RIPO1=spval
  allocate(this%RIPB1(NumMicrbHetetrophCmplx,1:jcplx));this%RIPB1=spval
  allocate(this%RINHOR(NumMicrbHetetrophCmplx,1:jcplx));this%RINHOR=spval
  allocate(this%RINOOR(NumMicrbHetetrophCmplx,1:jcplx));this%RINOOR=spval
  allocate(this%RIPOOR(NumMicrbHetetrophCmplx,1:jcplx));this%RIPOOR=spval
  allocate(this%RIPO1R(NumMicrbHetetrophCmplx,1:jcplx));this%RIPO1R=spval
  allocate(this%RO2DmndAutort(NumMicrobAutrophCmplx));this%RO2DmndAutort=spval
  allocate(this%RNH4UptkSoilAutor(NumMicrobAutrophCmplx));this%RNH4UptkSoilAutor=spval
  allocate(this%RNH4UptkBandAutor(NumMicrobAutrophCmplx));this%RNH4UptkBandAutor=spval
  allocate(this%RNO3UptkSoilAutor(NumMicrobAutrophCmplx));this%RNO3UptkSoilAutor=spval
  allocate(this%RNO3UptkBandAutor(NumMicrobAutrophCmplx));this%RNO3UptkBandAutor=spval
  allocate(this%RH2PO4UptkSoilAutor(NumMicrobAutrophCmplx));this%RH2PO4UptkSoilAutor=spval
  allocate(this%RH2PO4UptkBandAutor(NumMicrobAutrophCmplx));this%RH2PO4UptkBandAutor=spval
  allocate(this%RH1PO4UptkSoilAutor(NumMicrobAutrophCmplx));this%RH1PO4UptkSoilAutor=spval
  allocate(this%RH1PO4UptkBandAutor(NumMicrobAutrophCmplx));this%RH1PO4UptkBandAutor=spval
  allocate(this%RNH4UptkLitrAutor(NumMicrobAutrophCmplx));this%RNH4UptkLitrAutor=spval
  allocate(this%RNO3UptkLitrAutor(NumMicrobAutrophCmplx));this%RNO3UptkLitrAutor=spval
  allocate(this%RH2PO4UptkLitrAutor(NumMicrobAutrophCmplx));this%RH2PO4UptkLitrAutor=spval
  allocate(this%RH1PO4UptkLitrAutor(NumMicrobAutrophCmplx));this%RH1PO4UptkLitrAutor=spval
  allocate(this%RNH3OxidAutor(NumMicrobAutrophCmplx));this%RNH3OxidAutor=spval
  allocate(this%RNH3OxidAutorBand(NumMicrobAutrophCmplx));this%RNH3OxidAutorBand=spval
  allocate(this%RNO2OxidAutor(NumMicrobAutrophCmplx));this%RNO2OxidAutor=spval
  allocate(this%RNO2OxidAutorBand(NumMicrobAutrophCmplx));this%RNO2OxidAutorBand=spval

  end subroutine Init
!------------------------------------------------------------------------------------------
  subroutine ZeroOut(this)

  implicit none
  class(micfluxtype) :: this
  integer :: jcplx,JG,NumMicbFunGroups

  this%ROXSK = 0._r8
  this%RDOM_micb_flx=0._r8
  this%RO2DmndHetert=0._r8
  this%ROQCS=0._r8
  this%ROQAS=0._r8
  this%RVMX3=0._r8
  this%RVMB3=0._r8
  this%RVMX2=0._r8
  this%RVMB2=0._r8
  this%RVMX1=0._r8
  this%RVMX4=0._r8
  this%RVMB4=0._r8
  this%RINHO=0._r8
  this%RINHB=0._r8
  this%RINOO=0._r8
  this%RINOB=0._r8
  this%RIPOO=0._r8
  this%RIPBO=0._r8
  this%RIPO1=0._r8
  this%RIPB1=0._r8
  this%RINHOR=0._r8
  this%RINOOR=0._r8
  this%RIPOOR=0._r8
  this%RIPO1R=0._r8
  this%RO2DmndAutort=0._r8
  this%RNH4UptkSoilAutor=0._r8
  this%RNH4UptkBandAutor=0._r8
  this%RNO3UptkSoilAutor=0._r8
  this%RNO3UptkBandAutor=0._r8
  this%RH2PO4UptkSoilAutor=0._r8
  this%RH2PO4UptkBandAutor=0._r8
  this%RH1PO4UptkSoilAutor=0._r8
  this%RH1PO4UptkBandAutor=0._r8
  this%RNH4UptkLitrAutor=0._r8
  this%RNO3UptkLitrAutor=0._r8
  this%RH2PO4UptkLitrAutor=0._r8
  this%RH1PO4UptkLitrAutor=0._r8
  this%RNH3OxidAutor=0._r8
  this%RNH3OxidAutorBand=0._r8
  this%RNO2OxidAutor=0._r8
  this%RNO2OxidAutorBand=0._r8

  end subroutine ZeroOut
!------------------------------------------------------------------------------------------
  subroutine Destruct(this)
  use abortutils, only : destroy
  implicit none
  class(micfluxtype) :: this

  call destroy(this%RDOM_micb_flx)
  call destroy(this%RNH3OxidAutor)
  call destroy(this%RNH3OxidAutorBand)
  call destroy(this%RNO2OxidAutor)
  call destroy(this%RNO2OxidAutorBand)
  call destroy(this%RO2DmndHetert)
  call destroy(this%ROQCS)
  call destroy(this%ROQAS)
  call destroy(this%RVMX3)
  call destroy(this%RVMB3)
  call destroy(this%RVMX2)
  call destroy(this%RVMB2)
  call destroy(this%RVMX1)
  call destroy(this%ROXSK)
  call destroy(this%RVMX4)
  call destroy(this%RVMB4)
  call destroy(this%RINHO)
  call destroy(this%RINHB)
  call destroy(this%RINOO)
  call destroy(this%RINOB)
  call destroy(this%RIPOO)
  call destroy(this%RIPBO)
  call destroy(this%RIPO1)
  call destroy(this%RIPB1)
  call destroy(this%RINHOR)
  call destroy(this%RINOOR)
  call destroy(this%RIPOOR)
  call destroy(this%RIPO1R)
  call destroy(this%RNH4UptkSoilAutor)
  call destroy(this%RNH4UptkBandAutor)
  call destroy(this%RNO3UptkSoilAutor)
  call destroy(this%RNO3UptkBandAutor)
  call destroy(this%RH2PO4UptkSoilAutor)
  call destroy(this%RH2PO4UptkBandAutor)
  call destroy(this%RH1PO4UptkSoilAutor)
  call destroy(this%RH1PO4UptkBandAutor)
  call destroy(this%RNH4UptkLitrAutor)
  call destroy(this%RNO3UptkLitrAutor)
  call destroy(this%RH2PO4UptkLitrAutor)
  call destroy(this%RH1PO4UptkLitrAutor)
  call destroy(this%RO2DmndAutort)

  end subroutine Destruct
end module MicFluxTypeMod
