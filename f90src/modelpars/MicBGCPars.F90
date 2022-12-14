module MicBGCPars

! define microbial parameters
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
  use EcoSIMConfig
implicit none
  private
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public :: MicParType
  integer :: jcplx   !# of microbe-substrate complexes
  integer :: jsken   !# of kinetic components of the substrates
  integer :: jguilds !# of guilds
  integer :: NFGs    !# of functional groups
  integer :: k_woody_litr
  integer :: k_fine_litr
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
  integer :: ndbiomcp   !number of necrobiomass components
  integer :: nlbiomcp   !number of living biomass components

  real(r8), pointer :: OMCI(:,:)             !initializes microbial biomass
  real(r8), pointer :: OHCK(:)    !fractions of SOC in adsorbed C
  real(r8), pointer :: OMCK(:)    !fractions of SOC in biomass
  real(r8), pointer :: ORCK(:)    !fractions of SOC in litter
  real(r8), pointer :: OQCK(:)    !fractions of SOC in DOC
  real(r8), pointer :: ORCI(:,:)   !allocation of residue to kinetic components
  real(r8), pointer :: FL(:)       !allocation to microbial kinetic fractions
  real(r8), pointer :: CNOMC(:,:,:)        !maximum/minimum microbial N:C
  real(r8), pointer :: CPOMC(:,:,:)        !maximum/minimum microbial P:C
  real(r8), pointer :: CNOMCff(:,:)        !maximum/minimum microbial N:C
  real(r8), pointer :: CPOMCff(:,:)        !maximum/minimum microbial P:C
  real(r8), pointer :: CNOMCa(:,:,:)         !average maximum/minimum microbial N:C
  real(r8), pointer :: CPOMCa(:,:,:)         !average maximum/minimum microbial P:C
  real(r8), pointer :: CNOMCffa(:,:)         !average maximum/minimum microbial N:C
  real(r8), pointer :: CPOMCffa(:,:)         !average maximum/minimum microbial P:C
  real(r8), pointer :: DOSA(:)
  real(r8), pointer :: SPOSC(:,:)
  real(r8), pointer :: CNOFC(:,:)                         !fractions to allocate N to kinetic components
  real(r8), pointer :: CPOFC(:,:)                         !fractions to allocate P to kinetic components
  real(r8), pointer :: CNRH(:)                            !default N:C ratios in SOC complexes
  real(r8), pointer :: CPRH(:)                            !default P:C ratios in SOC complexes
  real(r8), pointer :: OMCF(:)                            !hetero microbial biomass composition in SOC
  real(r8), pointer :: OMCA(:)                            !autotrophic microbial biomass composition in SOC
  logical,  pointer :: is_activef_micb(:)
  logical,  pointer :: is_aerobic_hetr(:)
  logical,  pointer :: is_anerobic_hetr(:)
  logical,  pointer :: is_litter(:)
  logical,  pointer :: is_finelitter(:)
  logical,  pointer :: is_CO2_autotroph(:)
  character(len=16) :: kiname(0:jskenc-1)
  character(len=16) :: cplxname(1:jcplxc)
  character(len=16) :: hmicname(NFGsc)
  character(len=16) :: amicname(NFGsc)
  character(len=16) :: micresb(0:ndbiomcpc-1)      !residual biomass name
  character(len=16) :: micbiom(1:nlbiomcpc)        !microbial biomass pool name
  integer, pointer :: JGnio(:)   !guid indices for organic-microbial complex
  integer, pointer :: JGnfo(:)   !guid indices for organic-microbial complex
  integer, pointer :: JGniA(:)   !guid indices for autotrophic-microbial complex
  integer, pointer :: JGnfA(:)   !guid indices for autotrophic-microbial complex
  integer :: NMICBSA             !total number of microbial guilds in the autotrophic complex
  integer :: NMICBSO             !total number of microbial guilds in one organic-microbial complex
  integer :: n_litrsfk
  integer :: n_pltlitrk
  integer :: iprotein
  integer :: icarbhyro
  integer :: icellulos
  integer :: ilignin
  contains
    procedure, public  :: Init
    procedure, public  :: SetPars
    procedure, private :: InitAllocate
    procedure, public  :: destroy      =>DestructMicBGCPar
  end type MicParType


contains

  subroutine Init(this,nmicbguilds)
  implicit none
  class(MicParType) :: this
  integer, intent(in) :: nmicbguilds

  !organic matter is grouped into five complexes, including woody(1),
  ! non-woody(2), manure(3), litter, POC(4) and humus(5) (g Mg-1)

  this%ndbiomcp=ndbiomcpc  !number of necrobiomass components
  this%nlbiomcp=nlbiomcpc  !number of living biomass components

  this%jcplx=jcplxc         !# of microbe-substrate complexes
  this%jsken=jskenc        !# of kinetic components of the substrates
  this%jguilds=nmicbguilds
  this%NFGs=NFGsc
  !woody, non_woody litter and manure are defined as litter
  allocate(this%is_litter(1:this%jcplx));this%is_litter(:)=.false.
  allocate(this%is_finelitter(1:this%jcplx));this%is_finelitter(:)=.false.

  this%k_woody_litr=1;     this%is_litter(this%k_woody_litr)=.true.
  this%k_fine_litr=this%k_woody_litr+1; this%is_litter(this%k_fine_litr)=.true.
  this%n_pltlitrk=this%k_fine_litr
  this%k_manure=this%k_fine_litr+1;   this%is_litter(this%k_manure)=.true.
  this%is_litter(this%k_fine_litr)=.true.
  this%is_litter(this%k_manure)=.true.

  this%iprotein =1
  this%icarbhyro=2
  this%icellulos=3
  this%ilignin  =4

  this%n_litrsfk=this%k_manure
  this%k_POM=this%k_manure+1
  this%k_humus=this%k_POM+1
  this%kiname(0)='protein'
  this%kiname(1)='carbhydro'
  this%kiname(2)='cellulose'
  this%kiname(3)='lignin'
  this%cplxname(1)='woodylitr'
  this%cplxname(2)='nwoodylit'
  this%cplxname(3)='manure'
  this%cplxname(4)='pom'
  this%cplxname(5)='humus'
  this%hmicname(1)='aerohetrob'
  this%hmicname(2)='anerofaclb'
  this%hmicname(3)='aerofungi'
  this%hmicname(4)='aneroferm'
  this%hmicname(5)='acetomethg'
  this%hmicname(6)='aeron2fix'
  this%hmicname(7)='aneron2fix'
  this%amicname(1)='amoniaoxib'
  this%amicname(2)='nititeoxib'
  this%amicname(3)='aeromethtp'
  this%amicname(5)='hydromethg'
  this%amicname(4)='null'
  this%amicname(6)='null'
  this%amicname(7)='null'
  this%micresb(0)='labile'
  this%micresb(1)='resist'
  this%micbiom(1)='labile'
  this%micbiom(2)='resist'
  this%micbiom(3)='active'
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

  this%is_aerobic_hetr(this%n_aero_hetrophb)=.true.
  this%is_aerobic_hetr(this%n_anero_faculb)=.true.
  this%is_aerobic_hetr(this%n_aero_fungi)=.true.
  this%is_aerobic_hetr(this%n_aero_n2fixer)=.true.

  this%is_anerobic_hetr(this%n_anaero_ferm)=.true.
  this%is_anerobic_hetr(this%n_anero_n2fixer)=.true.
!the abstract complex
  this%nf_amonia_oxi=1
  this%nf_nitrite_oxi=2
  this%nf_aero_methanot=3
  this%nf_hydro_methang=5

  this%is_activef_micb(this%nf_amonia_oxi)=.True.
  this%is_activef_micb(this%nf_nitrite_oxi)=.True.
  this%is_activef_micb(this%nf_aero_methanot)=.True.
  this%is_activef_micb(this%nf_hydro_methang)=.True.

  this%is_CO2_autotroph(this%nf_amonia_oxi)=.true.
  this%is_CO2_autotroph(this%nf_nitrite_oxi)=.true.
  this%is_CO2_autotroph(this%nf_hydro_methang)=.true.


  end subroutine Init
!------------------------------------------------------------------------------------------

  subroutine SetPars(this)

  implicit none
  class(MicParType) :: this

  real(r8) :: COMCI(nlbiomcpc,1:this%jcplx)
  real(r8) :: OMCI1(nlbiomcpc,1:this%jcplx)  !allocation of biomass to kinetic components
  integer :: K,M,NGL,N
  associate(                    &
    OHCK     => this%OHCK     , &
    OMCK     => this%OMCK     , &
    ORCK     => this%ORCK     , &
    OQCK     => this%OQCK     , &
    ORCI     => this%ORCI     , &
    OMCI     => this%OMCI     , &
    CNOFC    => this%CNOFC    , &
    CPOFC    => this%CPOFC    , &
    CNOMC    => this%CNOMC    , &
    CPOMC    => this%CPOMC    , &
    CNOMCff  => this%CNOMCff  , &
    CPOMCff  => this%CPOMCff  , &
    FL       => this%FL       , &
    DOSA     => this%DOSA     , &
    SPOSC    => this%SPOSC    , &
    CNRH     => this%CNRH     , &
    CPRH     => this%CPRH     , &
    OMCF     => this%OMCF     , &
    OMCA     => this%OMCA     , &
    JG       => this%jguilds    &
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

  OMCI(1:3,:)=OMCI1

!  if(this%jguilds.GT.1)then
!    COMCI=OMCI(1:nlbiomcpc,:)
!    DO K=1,jcplxc
!      DO NGL=2,this%jguilds-1
!        DO M=1,nlbiomcpc
!          OMCI(M+(NGL-1)*nlbiomcpc,K)=OMCI(M,K)
!          COMCI(M,K)=COMCI(M,K)+OMCI(M,K)
!        enddo
!      enddo
!      DO M=1,nlbiomcpc
!        OMCI(M+(JG-1)*3,K)=OMCI1(M,K)-COMCI(M,K)
!      ENDDO
!    enddo
!  endif


! CNOFC,CPOFC=fractions to allocate N,P to kinetic components
! CNOMC,CPOMC=maximum N:C and P:C ratios in microbial biomass

  CNOFC(1:jskenc,this%k_woody_litr)=real((/0.0050,0.0050,0.0050,0.0200/),r8)  !woody
  CPOFC(1:jskenc,this%k_woody_litr)=real((/0.0005,0.0005,0.0005,0.0020/),r8)  !woody
  CNOFC(1:jskenc,this%k_fine_litr)=real((/0.0200,0.0200,0.0200,0.0200/),r8)  !non-woody
  CPOFC(1:jskenc,this%k_fine_litr)=real((/0.0020,0.0020,0.0020,0.0020/),r8)  !non-woody
  CNOFC(1:jskenc,this%k_manure)=real((/0.0200,0.0200,0.0200,0.0200/),r8)   !manure
  CPOFC(1:jskenc,this%k_manure)=real((/0.0020,0.0020,0.0020,0.0020/),r8)   !manure
  FL(1:2)=real((/0.55,0.45/),r8)

  D95: DO K=1,this%jcplx
    DO  N=1,this%NFGs
      IF(N.EQ.this%n_aero_fungi)THEN

        DO NGL=this%JGnio(n),this%JGnfo(n)
          CNOMC(1,NGL,K)=0.15_r8           !maximum
          CNOMC(2,NGL,K)=0.09_r8           !minimum
          CPOMC(1,NGL,K)=0.015_r8
          CPOMC(2,NGL,K)=0.009_r8
        ENDDO
        this%CNOMCa(1,N,K)=0.15_r8           !maximum
        this%CNOMCa(2,N,K)=0.09_r8           !minimum
        this%CPOMCa(1,N,K)=0.015_r8
        this%CPOMCa(2,N,K)=0.009_r8

      ELSE
        do NGL=this%JGnio(n),this%JGnfo(n)
          CNOMC(1,NGL,K)=0.225_r8
          CNOMC(2,NGL,K)=0.135_r8
          CPOMC(1,NGL,K)=0.0225_r8
          CPOMC(2,NGL,K)=0.0135_r8
        enddo
        this%CNOMCa(1,N,K)=0.225_r8
        this%CNOMCa(2,N,K)=0.135_r8
        this%CPOMCa(1,N,K)=0.0225_r8
        this%CPOMCa(2,N,K)=0.0135_r8
      ENDIF
      do NGL=this%JGnio(n),this%JGnfo(n)
        CNOMC(3,NGL,K)=DOT_PRODUCT(FL,CNOMC(1:2,NGL,K))
        CPOMC(3,NGL,K)=DOT_PRODUCT(FL,CPOMC(1:2,NGL,K))
      enddo
      this%CNOMCa(3,N,K)=DOT_PRODUCT(FL,this%CNOMCa(1:2,N,K))
      this%CPOMCa(3,N,K)=DOT_PRODUCT(FL,this%CPOMCa(1:2,N,K))
    enddo
  ENDDO D95

  DO  N=1,this%NFGs
    do NGL=this%JGniA(n),this%JGnfA(n)
      CNOMCff(1,NGL)=0.225_r8
      CNOMCff(2,NGL)=0.135_r8
      CPOMCff(1,NGL)=0.0225_r8
      CPOMCff(2,NGL)=0.0135_r8
    enddo
    this%CNOMCffa(1,N)=0.225_r8
    this%CNOMCffa(2,N)=0.135_r8
    this%CPOMCffa(1,N)=0.0225_r8
    this%CPOMCffa(2,N)=0.0135_r8
    do NGL=this%JGniA(n),this%JGnfA(n)
      CNOMCff(3,NGL)=DOT_PRODUCT(FL,CNOMCff(1:2,NGL))
      CPOMCff(3,NGL)=DOT_PRODUCT(FL,CPOMCff(1:2,NGL))
    enddo
    this%CNOMCffa(3,N)=DOT_PRODUCT(FL,this%CNOMCffa(1:2,N))
    this%CPOMCffa(3,N)=DOT_PRODUCT(FL,this%CPOMCffa(1:2,N))
  enddo

  end associate
  end subroutine SetPars
!------------------------------------------------------------------------------------------

  subroutine Initallocate(this)
  implicit none
  class(MicParType) :: this

  integer :: jguilds
  integer :: NFGs
  integer :: jcplx
  integer :: jsken
  integer :: n, k

  jguilds = this%jguilds
  NFGs  =this%NFGs
  jcplx =this%jcplx
  jsken =this%jsken
  allocate(this%JGnio(NFGs))
  allocate(this%JGnfo(NFGs))

  allocate(this%JGniA(NFGsc))
  allocate(this%JGnfA(NFGsc))

  k=1
  this%NMICBSA=0
  this%NMICBSO=0
  do n=1,NFGsc
    this%JGnio(n)=k
    this%JGniA(n)=k
    k=k+jguilds
    this%JGnfo(n)=k-1
    this%JGnfA(n)=k-1
    this%NMICBSA=this%NMICBSA+this%JGnfA(n)-this%JGniA(n)+1
    this%NMICBSO=this%NMICBSO+this%JGnfo(n)-this%JGnio(n)+1
  enddo


  allocate(this%DOSA(1:jcplx))
  allocate(this%SPOSC(jsken,1:jcplx))
  allocate(this%OHCK(1:jcplx))
  allocate(this%OMCK(1:jcplx))
  allocate(this%ORCK(1:jcplx))
  allocate(this%OQCK(1:jcplx))
  allocate(this%ORCI(ndbiomcpc,1:jcplx))
  allocate(this%OMCI(nlbiomcpc,1:jcplx))
  allocate(this%CNOMC(nlbiomcpc,this%NMICBSO,1:jcplx))
  allocate(this%CPOMC(nlbiomcpc,this%NMICBSO,1:jcplx))
  allocate(this%CNOMCff(nlbiomcpc,this%NMICBSA))
  allocate(this%CPOMCff(nlbiomcpc,this%NMICBSA))
  allocate(this%CNOMCa(nlbiomcpc,NFGs,1:jcplx))
  allocate(this%CPOMCa(nlbiomcpc,NFGs,1:jcplx))
  allocate(this%CNOMCffa(nlbiomcpc,NFGs))
  allocate(this%CPOMCffa(nlbiomcpc,NFGs))

  allocate(this%CNOFC(jsken,1:this%n_litrsfk))
  allocate(this%CPOFC(jsken,1:this%n_litrsfk))
  allocate(this%CNRH(1:jcplx))
  allocate(this%CPRH(1:jcplx))
  allocate(this%OMCF(NFGs))
  allocate(this%OMCA(NFGs))
  allocate(this%FL(2))
  allocate(this%is_activef_micb(NFGs)); this%is_activef_micb=.false.
  allocate(this%is_CO2_autotroph(NFGs)); this%is_CO2_autotroph=.false.
  allocate(this%is_aerobic_hetr(NFGs)); this%is_aerobic_hetr=.false.
  allocate(this%is_anerobic_hetr(NFGs));this%is_anerobic_hetr=.false.
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
  call destroy(this%is_CO2_autotroph)
  call destroy(this%is_aerobic_hetr)
  call destroy(this%is_anerobic_hetr)
  end subroutine DestructMicBGCPar

end module MicBGCPars
