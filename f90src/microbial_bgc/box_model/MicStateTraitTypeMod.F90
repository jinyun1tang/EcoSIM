module MicStateTraitTypeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only : micpar
  use TracerIDMod
  use abortutils, only : destroy
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  type, public :: micsttype
  real(r8) :: EPOC
  real(r8) :: EHUM
  real(r8) :: TOQCK
  real(r8) :: ZNH4B
  real(r8) :: ZNH4S
  real(r8) :: ZNO3B
  real(r8) :: ZNO3S
  real(r8) :: H1POB
  real(r8) :: H1PO4
  real(r8) :: ZNO2B
  real(r8) :: ZNO2S
  real(r8) :: H2POB
  real(r8) :: H2PO4
  real(r8) :: CCO2S
  real(r8) :: CNO2S
  real(r8) :: CNO2B
  real(r8) :: CNO3S
  real(r8) :: CNO3B
  real(r8) :: CNH4S
  real(r8) :: CNH4B
  real(r8) :: CZ2OS
  real(r8) :: Z2OS
  real(r8) :: OXYG
  real(r8) :: OXYS
  real(r8) :: COXYS
  real(r8) :: SOXYL
  real(r8) :: COXYG
  real(r8) :: CH2P4
  real(r8) :: CH2P4B
  real(r8) :: CH1P4
  real(r8) :: CH1P4B
  real(r8) :: CZ2GS
  real(r8) :: CH2GS
  real(r8) :: H2GS
  real(r8) :: CCH4G
  real(r8) :: CH4S
  real(r8) :: SCH4L
  real(r8) :: ZNFN0
  real(r8) :: ZNFNI
  real(r8) ::ZNH4TU
  real(r8) ::ZNO3TU
  real(r8) ::H1P4TU
  real(r8) ::H2P4TU
  real(r8) ::CNH4BU
  real(r8) ::CNH4SU
  real(r8) ::CH2P4U
  real(r8) ::CH2P4BU
  real(r8) ::CH1P4U
  real(r8) ::CH1P4BU
  real(r8) ::CNO3SU
  real(r8) ::CNO3BU
  real(r8) :: TFNQ
  real(r8) :: VOLQ
  real(r8),allocatable :: FOSRH(:)
  real(r8),allocatable :: DOM(:,:)
  real(r8),allocatable :: OHM(:,:)
  real(r8),allocatable :: OSA(:,:)
  real(r8),allocatable :: OSM(:,:,:)
  real(r8),allocatable :: ORM(:,:,:)
  real(r8),allocatable :: CNOSC(:,:)
  real(r8),allocatable :: CPOSC(:,:)
  real(r8),allocatable :: OMEhetr(:,:,:)
  real(r8),allocatable :: OMEauto(:,:)
  real(r8) :: OSC13U,OSC14U,OSC24U
  real(r8) :: OSN13U,OSN14U,OSN24U
  real(r8) :: OSP13U,OSP14U,OSP24U
  contains
  procedure, public :: Init
  procedure, public :: Destroy=> Destruct
  end type micsttype
  contains

!------------------------------------------------------------------------------------------
  subroutine Init(this)
  implicit none
  class(micsttype) :: this
  integer :: jcplx,NumMicbFunGroups,jsken
  integer, pointer :: ndbiomcp, nlbiomcp
  integer, pointer :: NumMicrobAutotrophCmplx, NumMicrbHetetrophCmplx
  integer, pointer :: NumLiveHeterBioms
  integer, pointer :: NumLiveAutoBioms


  jcplx=micpar%jcplx
  NumLiveAutoBioms => micpar%NumLiveAutoBioms
  NumMicbFunGroups=micpar%NumMicbFunGroups
  jsken=micpar%jsken
  ndbiomcp =>micpar%ndbiomcp
  nlbiomcp =>micpar%nlbiomcp
  NumMicrobAutotrophCmplx=>micpar%NumMicrobAutotrophCmplx
  NumMicrbHetetrophCmplx=>micpar%NumMicrbHetetrophCmplx
  NumLiveHeterBioms => micpar%NumLiveHeterBioms

  allocate(this%FOSRH(1:jcplx))
  allocate(this%DOM(idom_beg:idom_end,1:jcplx))
  allocate(this%OHM(idom_beg:idom_end,1:jcplx))
  allocate(this%OSA(jsken,1:jcplx))
  allocate(this%OSM(NumPlantChemElmnts,jsken,1:jcplx))
  allocate(this%ORM(1:NumPlantChemElmnts,ndbiomcp,1:jcplx))
  allocate(this%CNOSC(jsken,1:jcplx))
  allocate(this%CPOSC(jsken,1:jcplx))
  allocate(this%OMEhetr(NumPlantChemElmnts,NumLiveHeterBioms,1:jcplx))
  allocate(this%OMEauto(NumPlantChemElmnts,NumLiveAutoBioms))

  end subroutine Init

!------------------------------------------------------------------------------------------


  subroutine Destruct(this)

  implicit none
  class(micsttype) :: this

  call destroy(this%FOSRH)
  call destroy(this%DOM)
  call destroy(this%OHM)
  call destroy(this%OSA)
  call destroy(this%OSM)
  call destroy(this%ORM)
  call destroy(this%CNOSC)
  call destroy(this%CPOSC)
  call destroy(this%OMEhetr)
  call destroy(this%OMEauto)
  end subroutine Destruct
end module MicStateTraitTypeMod
