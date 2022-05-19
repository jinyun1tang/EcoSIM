module MicBGCPars

! define microbial parameters
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
implicit none
  private

  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public :: MicParType
  integer :: jcplx   !# of microbe-substrate complexes
  integer :: jcplx1  !jcplx - 1
  integer :: jsken   !# of kinetic components of the substrates
  integer :: jguilds !# of guilds
  integer :: NFGs    !# of functional groups
  real(r8), allocatable :: OMCI(:,:)             !initializes microbial biomass
  real(r8), allocatable :: OHCK(:)    !fractions of SOC in adsorbed C
  real(r8), allocatable :: OMCK(:)    !fractions of SOC in biomass
  real(r8), allocatable :: ORCK(:)    !fractions of SOC in litter
  real(r8), allocatable :: OQCK(:)    !fractions of SOC in DOC
  real(r8), allocatable :: ORCI(:,:)   !allocation of residue to kinetic components
  real(r8), allocatable :: FL(:)       !allocation to microbial kinetic fractions
  real(r8),allocatable :: CNOMC(:,:,:,:)        !maximum/minimum microbial N:C
  real(r8),allocatable :: CPOMC(:,:,:,:)        !maximum/minimum microbial P:C
  real(r8),allocatable :: CNOMCff(:,:,:)        !maximum/minimum microbial N:C
  real(r8),allocatable :: CPOMCff(:,:,:)        !maximum/minimum microbial P:C
  real(r8), allocatable :: DOSA(:)
  real(r8), allocatable :: SPOSC(:,:)
  real(r8),allocatable :: CNOFC(:,:)                         !fractions to allocate N to kinetic components
  real(r8),allocatable :: CPOFC(:,:)                         !fractions to allocate P to kinetic components
  real(r8),allocatable :: CNRH(:)                            !default N:C ratios in SOC complexes
  real(r8),allocatable :: CPRH(:)                            !default P:C ratios in SOC complexes
  real(r8),allocatable :: OMCF(:)                            !hetero microbial biomass composition in SOC
  real(r8),allocatable :: OMCA(:)                            !autotrophic microbial biomass composition in SOC
  logical, allocatable :: is_activef_micb(:)
  integer :: k_woody_litr
  integer :: k_non_woody_litr
  integer :: k_manure
  integer :: k_POM
  integer :: k_humus
  integer :: nf_amonia_oxi
  integer :: nf_nitrite_oxi
  integer :: nf_aero_methanot
  integer :: nf_hydro_methang
  integer :: n_aero_hetrophb
  integer :: n_anero_faculb
  integer :: n_aero_fungi
  integer :: n_anaero_ferm
  integer :: n_aceto_methang
  integer :: n_aero_n2fixer
  integer :: n_anero_n2fixer
  contains
    procedure, public  :: Init
    procedure, public  :: SetPars
    procedure, private :: InitAllocate
    procedure, public  :: destroy      =>DestructMicBGCPar
  end type MicParType

  type(MicParType), public :: micpar
contains

  subroutine Init(this,nmicbguilds)
  implicit none
  class(MicParType) :: this
  integer, intent(in) :: nmicbguilds
  !organic matter is grouped into five complexes, including woody(0),
  ! non-woody(1), manure(2), litter, POC(3) and humus(4) (g Mg-1)

  this%jcplx=5 !# of microbe-substrate complexes
  this%jcplx1=this%jcplx-1
  this%jsken=4 !# of kinetic components of the substrates
  this%jguilds=nmicbguilds
  this%NFGs=7

  this%k_woody_litr=0
  this%k_non_woody_litr=1
  this%k_manure=2
  this%k_POM=3
  this%k_humus=4


  call this%Initallocate()

!set up functional group ids
!five om-complexes
  this%n_aero_hetrophb=1
  this%n_anero_faculb=2
  this%n_aero_fungi=3
  this%n_anaero_ferm=4
  this%n_aceto_methang=5
  this%n_aero_n2fixer=6
  this%n_anero_n2fixer=7

!the abstract complex
  this%nf_amonia_oxi=1
  this%nf_nitrite_oxi=2
  this%nf_aero_methanot=3
  this%nf_hydro_methang=5

  this%is_activef_micb(this%nf_amonia_oxi)=.True.
  this%is_activef_micb(this%nf_nitrite_oxi)=.True.
  this%is_activef_micb(this%nf_aero_methanot)=.True.
  this%is_activef_micb(this%nf_hydro_methang)=.True.
  end subroutine Init
!------------------------------------------------------------------------------------------

  subroutine SetPars(this)

  implicit none
  class(MicParType) :: this

  real(r8) :: COMCI(3,0:this%jcplx1)
  real(r8) :: OMCI1(3,0:this%jcplx1)  !allocation of biomass to kinetic components
  integer :: K,M,NGL,N
  associate(               &
    OHCK  => this%OHCK     , &
    OMCK  => this%OMCK     , &
    ORCK  => this%ORCK     , &
    OQCK  => this%OQCK     , &
    ORCI  => this%ORCI     , &
    OMCI  => this%OMCI     , &
    CNOFC  => this%CNOFC     , &
    CPOFC  => this%CPOFC     , &
    CNOMC  => this%CNOMC     , &
    CPOMC  => this%CPOMC     , &
    CNOMCff  => this%CNOMCff     , &
    CPOMCff  => this%CPOMCff     , &
    FL      => this%FL          , &
    DOSA    => this%DOSA       , &
    SPOSC   => this%SPOSC      , &
    CNRH    => this%CNRH       , &
    CPRH    => this%CPRH       , &
    OMCF    => this%OMCF       , &
    OMCA    => this%OMCA       , &
    JG    => this%jguilds      &
  )
  OHCK=real((/0.05,0.05,0.05,0.05,0.05/),r8)
  OMCK=real((/0.01,0.01,0.01,0.01,0.01/),r8)
  ORCK=real((/0.25,0.25,0.25,0.25,0.25/),r8)
  OQCK=real((/0.005,0.005,0.005,0.005,0.005/),r8)
  OMCI1=reshape(real((/0.010,0.050,0.005,0.050,0.050,0.005,0.050,0.050,0.005, &
     0.010,0.050,0.005,0.010,0.050,0.005/),r8),shape(OMCI1))
  ORCI=reshape(real((/0.01,0.05,0.01,0.05,0.01,0.05 &
     ,0.001,0.005,0.001,0.005/),r8),shape(ORCI))

  DOSA=(/0.25E-03_r8,0.25_r8,0.25_r8,0.25_r8,0.25_r8/)
  SPOSC=reshape((/7.5_r8,7.5_r8,1.5_r8,0.5_r8,7.5_r8,7.5_r8,1.5_r8,0.5_r8 &
    ,7.5_r8,7.5_r8,1.5_r8,0.5_r8,0.05_r8,0.00_r8,0.00_r8,0.00_r8 &
    ,0.05_r8,0.0167_r8,0.00_r8,0.00_r8/),shape(sposc))

  CNRH=(/3.33E-02_r8,3.33E-02_r8,3.33E-02_r8,5.00E-02_r8,12.50E-02_r8/)
  CPRH=(/3.33E-03_r8,3.33E-03_r8,3.33E-03_r8,5.00E-03_r8,12.50E-03_r8/)
  OMCF=(/0.20_r8,0.20_r8,0.30_r8,0.20_r8,0.050_r8,0.025_r8,0.025_r8/)
  OMCA=(/0.06_r8,0.02_r8,0.01_r8,0.0_r8,0.01_r8,0.0_r8,0.0_r8/)

  OMCI(1:3,:)=OMCI1/real(JG,r8)

  if(this%jguilds.GT.1)then
    COMCI=OMCI(1:3,:)
    DO K=0,4
      DO NGL=2,this%jguilds-1
        DO M=1,3
          OMCI(M+(NGL-1)*3,K)=OMCI(M,K)
          COMCI(M,K)=COMCI(M,K)+OMCI(M,K)
        enddo
      enddo
      DO M=1,3
        OMCI(M+(JG-1)*3,K)=OMCI1(M,K)-COMCI(M,K)
      ENDDO
    enddo
  endif


! CNOFC,CPOFC=fractions to allocate N,P to kinetic components
! CNOMC,CPOMC=maximum N:C and P:C ratios in microbial biomass

  CNOFC(1:4,0)=real((/0.0050,0.0050,0.0050,0.0200/),r8)
  CPOFC(1:4,0)=real((/0.0005,0.0005,0.0005,0.0020/),r8)
  CNOFC(1:4,1)=real((/0.0200,0.0200,0.0200,0.0200/),r8)
  CPOFC(1:4,1)=real((/0.0020,0.0020,0.0020,0.0020/),r8)
  CNOFC(1:4,2)=real((/0.0200,0.0200,0.0200,0.0200/),r8)
  CPOFC(1:4,2)=real((/0.0020,0.0020,0.0020,0.0020/),r8)
  FL(1:2)=real((/0.55,0.45/),r8)

  DO 95 K=0,4
    DO  N=1,7
      IF(N.EQ.3)THEN
        DO NGL=1,JG
          CNOMC(1,NGL,N,K)=0.15_r8           !maximum
          CNOMC(2,NGL,N,K)=0.09_r8           !minimum
          CPOMC(1,NGL,N,K)=0.015_r8
          CPOMC(2,NGL,N,K)=0.009_r8
        ENDDO
      ELSE
        do NGL=1,JG
          CNOMC(1,NGL,N,K)=0.225_r8
          CNOMC(2,NGL,N,K)=0.135_r8
          CPOMC(1,NGL,N,K)=0.0225_r8
          CPOMC(2,NGL,N,K)=0.0135_r8
        enddo
      ENDIF
      do NGL=1,JG
        CNOMC(3,NGL,N,K)=FL(1)*CNOMC(1,NGL,N,K)+FL(2)*CNOMC(2,NGL,N,K)
        CPOMC(3,NGL,N,K)=FL(1)*CPOMC(1,NGL,N,K)+FL(2)*CPOMC(2,NGL,N,K)
      enddo
     enddo
95  CONTINUE

    DO  N=1,7
        do NGL=1,JG
          CNOMCff(1,NGL,N)=0.225_r8
          CNOMCff(2,NGL,N)=0.135_r8
          CPOMCff(1,NGL,N)=0.0225_r8
          CPOMCff(2,NGL,N)=0.0135_r8
        enddo
      do NGL=1,JG
        CNOMCff(3,NGL,N)=FL(1)*CNOMCff(1,NGL,N)+FL(2)*CNOMCff(2,NGL,N)
        CPOMCff(3,NGL,N)=FL(1)*CPOMCff(1,NGL,N)+FL(2)*CPOMCff(2,NGL,N)
      enddo
     enddo

  end associate
  end subroutine SetPars
!------------------------------------------------------------------------------------------

  subroutine Initallocate(this)
  implicit none
  class(MicParType) :: this

  integer :: jguilds
  integer :: NFGs
  integer :: jcplx1,jcplx
  integer :: jsken

  jguilds = this%jguilds
  NFGs  =this%NFGs
  jcplx =this%jcplx
  jcplx1=this%jcplx1
  jsken =this%jsken
  allocate(this%DOSA(0:jcplx1))
  allocate(this%SPOSC(jsken,0:jcplx1))
  allocate(this%OHCK(0:jcplx1))
  allocate(this%OMCK(0:jcplx1))
  allocate(this%ORCK(0:jcplx1))
  allocate(this%OQCK(0:jcplx1))
  allocate(this%ORCI(2,0:jcplx1))
  allocate(this%OMCI(3*jguilds,0:jcplx1))
  allocate(this%CNOMC(3,jguilds,NFGs,0:jcplx1))
  allocate(this%CPOMC(3,jguilds,NFGs,0:jcplx1))
  allocate(this%CNOMCff(3,jguilds,NFGs))
  allocate(this%CPOMCff(3,jguilds,NFGs))
  allocate(this%CNOFC(jsken,0:2))
  allocate(this%CPOFC(jsken,0:2))
  allocate(this%CNRH(0:jcplx1))
  allocate(this%CPRH(0:jcplx1))
  allocate(this%OMCF(NFGs))
  allocate(this%OMCA(NFGs))
  allocate(this%FL(2))
  allocate(this%is_activef_micb(NFGs)); this%is_activef_micb=.false.
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructMicBGCPar(this)
  implicit none
  class(MicParType) :: this

  call destroy(this%OMCI)
  call destroy(this%OHCK)
  call destroy(this%OMCK)
  call destroy(this%ORCK)
  call destroy(this%OQCK)
  call destroy(this%ORCI)
  call destroy(this%OMCI)
  call destroy(this%CNOMC)
  call destroy(this%CPOMC)
  call destroy(this%CNOMCff)
  call destroy(this%CPOMCff)
  call destroy(this%CNOFC)
  call destroy(this%CPOFC)
  call destroy(this%FL)
  call destroy(this%DOSA)
  call destroy(this%SPOSC)
  call destroy(this%CNOFC)
  call destroy(this%CPOFC)
  call destroy(this%CNRH)
  call destroy(this%CPRH)
  call destroy(this%OMCF)
  call destroy(this%OMCA)
  call destroy(this%is_activef_micb)
  end subroutine DestructMicBGCPar

end module MicBGCPars
