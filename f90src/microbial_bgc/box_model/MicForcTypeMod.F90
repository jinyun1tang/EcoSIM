module MicForcTypeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only : micpar
  use EcoSIMSolverPar, only : NPH
  use abortutils, only : destroy
  implicit none

  private
  character(len=*), parameter :: mod_filename = __FILE__

  type, public :: micforctype

  real(r8) :: CCH4E
  real(r8) :: COXYE
  real(r8) :: COXQ
  real(r8) :: COXR
  real(r8) :: FLQRI
  real(r8) :: FLQRQ
  real(r8) :: offset
  real(r8) :: VLitR
  real(r8) :: VWatLitrX
  real(r8) :: ZEROS2
  real(r8) :: ZEROS
  real(r8) :: VLWatMicP
  real(r8) :: VLsoiAirP
  real(r8) :: VOLW0
  real(r8) :: VLSoilMicP
  real(r8) :: THETY
  real(r8) :: POROS
  real(r8) :: TKS
  real(r8) :: FieldCapacity
  real(r8) :: THETW
  real(r8) :: pH
  real(r8) :: ZERO
  real(r8) :: SoilMicPMassLayer
  real(r8) :: VLSoilPoreMicP
  real(r8) :: TFND
  real(r8) :: VLNOB
  real(r8) :: VLNO3
  real(r8) :: VLNH4
  real(r8) :: VLNHB
  real(r8) :: VLPO4
  real(r8) :: VLPOB
  real(r8) :: AEC
  real(r8) :: PSISoilMatricP
  real(r8) :: OLSGL
  real(r8) :: ORGC
  real(r8) :: RNO2Y
  real(r8) :: RN2BY
  real(r8) :: RN2OY
  real(r8) :: ROXYY
  real(r8) :: ROXYF
  real(r8) :: RCH4L
  real(r8) :: RCH4F
  real(r8) :: RNH4Y
  real(r8) :: RNHBY
  real(r8) :: RNO3Y
  real(r8) :: RN3BY
  real(r8) :: RPO4Y
  real(r8) :: RPOBY
  real(r8) :: RP14Y
  real(r8) :: RP1BY
  real(r8) :: SoilMicPMassLayer0
  logical  :: LitrM
  logical  :: Lsurf
  real(r8) :: RNH4YU
  real(r8) :: RNO3YU
  real(r8) :: RPO4YU
  real(r8) :: RP14YU
  real(r8) :: ROXYL
  real(r8) :: VOLWU
  real(r8), allocatable :: CFOMC(:)
  real(r8), allocatable :: CFOMCU(:)
  real(r8), allocatable :: ROQCY(:)
  real(r8), allocatable :: ROQAY(:)
  real(r8), allocatable :: DFGS(:)  !rate constant for air-water gas exchange
  real(r8), allocatable :: FILM(:)
  real(r8), allocatable :: THETPM(:)
  real(r8), allocatable :: VLWatMicPM(:)
  real(r8), allocatable :: TORT(:)
  real(r8), allocatable :: VLsoiAirPM(:)
  contains
   procedure, public :: Init
   procedure, public :: destroy=>Destruct
  end type micforctype

  contains
!------------------------------------------------------------------------------------------

  subroutine Init(this)

  implicit none
  class(micforctype) :: this
  integer :: jcplx
  jcplx=micpar%jcplx
  allocate(this%CFOMC(1:micpar%ndbiomcp))
  allocate(this%CFOMCU(1:micpar%ndbiomcp))
  allocate(this%ROQCY(1:jcplx))
  allocate(this%ROQAY(1:jcplx))
  allocate(this%VLWatMicPM(NPH))
  allocate(this%THETPM(NPH))
  allocate(this%FILM(NPH))
  allocate(this%TORT(NPH))
  allocate(this%VLsoiAirPM(NPH))
  allocate(this%DFGS(NPH))
  end subroutine Init
!------------------------------------------------------------------------------------------


  subroutine Destruct(this)
  use abortutils, only : destroy
  implicit none
  class(micforctype) :: this

  call destroy(this%VLWatMicPM)
  call destroy(this%THETPM)
  call destroy(this%FILM)
  call destroy(this%TORT)
  call destroy(this%VLsoiAirPM)
  call destroy(this%DFGS)
  call destroy(this%ROQCY)
  call destroy(this%ROQAY)
  call destroy(this%CFOMC)
  call destroy(this%CFOMCU)

  end subroutine Destruct

end module MicForcTypeMod
