module MicBGCPars

! define microbial parameters
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use EcoSIMConfig
implicit none
  private
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  type, public :: MicParType
  integer :: jcplx   !# of microbe-substrate complexes
  integer :: jsken   !# of kinetic components of the substrates
  integer :: jguilds !# of guilds
  integer :: NumMicbFunGrupsPerCmplx    !# of functional groups
  integer :: k_woody_litr
  integer :: k_fine_litr
  integer :: k_manure
  integer :: k_POM
  integer :: k_humus
  integer :: mid_AmmoniaOxidBacter
  integer :: mid_NitriteOxidBacter
  integer :: mid_AerobicMethanotrofBacter
  integer :: mid_H2GenoMethanogArchea
  integer :: mid_Aerob_HeteroBacter
  integer :: mid_Facult_DenitBacter
  integer :: mid_Aerob_Fungi
  integer :: mid_fermentor
  integer :: mid_AcetoMethanogArchea
  integer :: mid_aerob_N2Fixer
  integer :: mid_Anaerob_N2Fixer
  integer :: ndbiomcp   !number of necrobiomass components
  integer :: nlbiomcp   !number of living biomass components

  real(r8), pointer :: OMCI(:,:)             !initializes microbial biomass
  real(r8), pointer :: OHCK(:)    !fractions of SOC in adsorbed C
  real(r8), pointer :: OMCK(:)    !fractions of SOC in biomass
  real(r8), pointer :: ORCK(:)    !fractions of SOC in litter
  real(r8), pointer :: OQCK(:)    !fractions of SOC in DOC
  real(r8), pointer :: ORCI(:,:)   !allocation of residue to kinetic components
  real(r8), pointer :: FL(:)       !allocation to microbial kinetic fractions
  real(r8), pointer :: rNCOMC(:,:,:)        !maximum/minimum mass based microbial N:C
  real(r8), pointer :: rPCOMC(:,:,:)        !maximum/minimum mass based microbial P:C
  real(r8), pointer :: rNCOMCAutor(:,:)        !maximum/minimum mass based microbial N:C
  real(r8), pointer :: rPCOMCAutor(:,:)        !maximum/minimum mass based microbial P:C
  real(r8), pointer :: rNCOMCa(:,:,:)         !average maximum/minimum mass based microbial N:C
  real(r8), pointer :: rPCOMCa(:,:,:)         !average maximum/minimum mass based microbial P:C
  real(r8), pointer :: rNCOMCAutora(:,:)         !average maximum/minimum mass based microbial N:C
  real(r8), pointer :: rPCOMCAutora(:,:)         !average maximum/minimum mass based microbial P:C
  real(r8), pointer :: DOSA(:)
  real(r8), pointer :: SPOSC(:,:)
  real(r8), pointer :: CNOFC(:,:)                         !fractions to allocate N to kinetic components
  real(r8), pointer :: CPOFC(:,:)                         !fractions to allocate P to kinetic components
  real(r8), pointer :: CNRH(:)                            !default N:C ratios in SOC complexes
  real(r8), pointer :: CPRH(:)                            !default P:C ratios in SOC complexes
  real(r8), pointer :: OMCF(:)                            !hetero microbial biomass composition in SOC
  real(r8), pointer :: OMCA(:)                            !autotrophic microbial biomass composition in SOC
  logical,  pointer :: is_activeMicrbFungrpAutor(:)
  logical,  pointer :: is_aerobic_hetr(:)
  logical,  pointer :: is_anaerobic_hetr(:)
  logical,  pointer :: is_litter(:)
  logical,  pointer :: is_finelitter(:)
  logical,  pointer :: is_CO2_autotroph(:)
  character(len=16) :: kiname(0:jskenc-1)
  character(len=16) :: cplxname(1:jcplxc)
  character(len=16) :: hmicname(NumMicbFunGrupsPerCmplx)
  character(len=16) :: amicname(NumMicbFunGrupsPerCmplx)
  character(len=16) :: micresb(0:NumDeadMicrbCompts-1)      !residual biomass name
  character(len=16) :: micbiom(1:NumLiveMicrbCompts)        !microbial biomass pool name
  integer, pointer :: JGnio(:)   !guid indices for organic-microbial complex
  integer, pointer :: JGnfo(:)   !guid indices for organic-microbial complex
  integer, pointer :: JGniA(:)   !guid indices for autotrophic-microbial complex
  integer, pointer :: JGnfA(:)   !guid indices for autotrophic-microbial complex
  integer :: NumMicrobAutrophCmplx        !total number of microbial guilds in the autotrophic complex
  integer :: NumHetetrMicCmplx         !total number of microbial guilds in one organic-microbial complex
  integer :: NumOfLitrCmplxs                !number of litter organo-microbial complexes, plant litter + manure
  integer :: NumOfPlantLitrCmplxs           !number of plant litter complexs, woody + fine litter
  integer :: NumLiveHeterBioms              !total number of biomass components in one heterotrophic OM complexes
  integer :: NumLiveAutoBioms               !total number of biomass components in the autotroph complex
  integer :: iprotein
  integer :: icarbhyro
  integer :: icellulos
  integer :: ilignin
  contains
    procedure, public  :: Init
    procedure, public  :: SetPars
    procedure, public  :: get_micb_id
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

  this%ndbiomcp = NumDeadMicrbCompts  !number of necrobiomass components
  this%nlbiomcp = NumLiveMicrbCompts  !number of living biomass components
  
  this%jcplx                   = jcplxc         !# of microbe-substrate complexes
  this%jsken                   = jskenc         !# of kinetic components of the substrates
  this%jguilds                 = nmicbguilds
  this%NumMicbFunGrupsPerCmplx = NumMicbFunGrupsPerCmplx

  !woody, non_woody litter and manure are defined as litter
  allocate(this%is_litter(1:this%jcplx));this%is_litter(:)=.false.
  allocate(this%is_finelitter(1:this%jcplx));this%is_finelitter(:)=.false.

  this%k_woody_litr                = 1;                   this%is_litter(this%k_woody_litr) = .true.
  this%k_fine_litr                 = this%k_woody_litr+1; this%is_litter(this%k_fine_litr) = .true.
  this%NumOfPlantLitrCmplxs        = this%k_fine_litr
  this%k_manure                    = this%k_fine_litr+1;   this%is_litter(this%k_manure)   = .true.
  this%is_litter(this%k_fine_litr) = .true.
  this%is_litter(this%k_manure)    = .true.

  this%iprotein  = 1
  this%icarbhyro = 2
  this%icellulos = 3
  this%ilignin   = 4

  this%NumOfLitrCmplxs = this%k_manure
  this%k_POM           = this%k_manure+1
  this%k_humus         = this%k_POM+1
  
  this%kiname(0)   = 'protein'
  this%kiname(1)   = 'carbhydro'
  this%kiname(2)   = 'cellulose'
  this%kiname(3)   = 'lignin'
  this%cplxname(1) = 'woodylitr'
  this%cplxname(2) = 'nwoodylit'
  this%cplxname(3) = 'manure'
  this%cplxname(4) = 'pom'
  this%cplxname(5) = 'humus'
  this%hmicname(1) = 'aerohetrob'
  this%hmicname(2) = 'anerofaclb'
  this%hmicname(3) = 'aerofungi'
  this%hmicname(4) = 'aneroferm'
  this%hmicname(5) = 'acetMicBiome_colhg'
  this%hmicname(6) = 'aeron2fix'
  this%hmicname(7) = 'aneron2fix'
  this%amicname(1) = 'amoniaoxib'
  this%amicname(2) = 'nititeoxib'
  this%amicname(3) = 'aeromethtp'
  this%amicname(5) = 'hydromethg'
  this%amicname(4) = 'null'
  this%amicname(6) = 'null'
  this%amicname(7) = 'null'
  this%micresb(0)  = 'labile'
  this%micresb(1)  = 'resist'
  this%micbiom(1)  = 'labile'
  this%micbiom(2)  = 'resist'
  this%micbiom(3)  = 'active'

  call this%Initallocate()

!set up functional group ids
! five om-complexes
  this%mid_Aerob_HeteroBacter  = 1
  this%mid_Facult_DenitBacter  = 2
  this%mid_Aerob_Fungi         = 3
  this%mid_fermentor           = 4
  this%mid_AcetoMethanogArchea = 5
  this%mid_aerob_N2Fixer       = 6
  this%mid_Anaerob_N2Fixer     = 7

  this%is_aerobic_hetr(this%mid_Aerob_HeteroBacter) = .true.
  this%is_aerobic_hetr(this%mid_Facult_DenitBacter) = .true.
  this%is_aerobic_hetr(this%mid_Aerob_Fungi)        = .true.
  this%is_aerobic_hetr(this%mid_aerob_N2Fixer)      = .true.

  this%is_anaerobic_hetr(this%mid_fermentor)       = .true.
  this%is_anaerobic_hetr(this%mid_Anaerob_N2Fixer) = .true.
!the autotrophic complex
  this%mid_AmmoniaOxidBacter        = 1
  this%mid_NitriteOxidBacter        = 2
  this%mid_AerobicMethanotrofBacter = 3
  this%mid_H2GenoMethanogArchea     = 5

  this%is_activeMicrbFungrpAutor(this%mid_AmmoniaOxidBacter)        = .True.
  this%is_activeMicrbFungrpAutor(this%mid_NitriteOxidBacter)        = .True.
  this%is_activeMicrbFungrpAutor(this%mid_AerobicMethanotrofBacter) = .True.
  this%is_activeMicrbFungrpAutor(this%mid_H2GenoMethanogArchea)     = .True.

  this%is_CO2_autotroph(this%mid_AmmoniaOxidBacter)    = .true.
  this%is_CO2_autotroph(this%mid_NitriteOxidBacter)    = .true.
  this%is_CO2_autotroph(this%mid_H2GenoMethanogArchea) = .true.

  end subroutine Init
!------------------------------------------------------------------------------------------

  subroutine SetPars(this)

  implicit none
  class(MicParType) :: this

  real(r8) :: COMCI(NumLiveMicrbCompts,1:this%jcplx)
  real(r8) :: OMCI1(NumLiveMicrbCompts,1:this%jcplx)  !allocation of biomass to kinetic components
  integer :: K,M,NGL,N
  associate(                         &
    OHCK        => this%OHCK,        &
    OMCK        => this%OMCK,        &
    ORCK        => this%ORCK,        &
    OQCK        => this%OQCK,        &
    ORCI        => this%ORCI,        &
    OMCI        => this%OMCI,        &
    CNOFC       => this%CNOFC,       &
    CPOFC       => this%CPOFC,       &
    rNCOMC      => this%rNCOMC,      &
    rPCOMC      => this%rPCOMC,      &
    rNCOMCAutor => this%rNCOMCAutor, &
    rPCOMCAutor => this%rPCOMCAutor, &
    FL          => this%FL,          &
    DOSA        => this%DOSA,        &
    SPOSC       => this%SPOSC,       &
    CNRH        => this%CNRH,        &
    CPRH        => this%CPRH,        &
    OMCF        => this%OMCF,        &
    OMCA        => this%OMCA,        &
    JG          => this%jguilds      &
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

  OMCI(1:NumLiveMicrbCompts,:)=OMCI1

!  if(this%jguilds.GT.1)then
!    COMCI=OMCI(1:NumLiveMicrbCompts,:)
!    DO K=1,jcplxc
!      DO NGL=2,this%jguilds-1
!        DO M=1,NumLiveMicrbCompts
!          OMCI(M+(NGL-1)*NumLiveMicrbCompts,K)=OMCI(M,K)
!          COMCI(M,K)=COMCI(M,K)+OMCI(M,K)
!        enddo
!      enddo
!      DO M=1,NumLiveMicrbCompts
!        OMCI(M+(JG-1)*3,K)=OMCI1(M,K)-COMCI(M,K)
!      ENDDO
!    enddo
!  endif


! CNOFC,CPOFC=fractions to allocate N,P to kinetic components
! rNCOMC,rPCOMC=maximum N:C and P:C ratios in microbial biomass

  CNOFC(1:jskenc,this%k_woody_litr) = real((/0.0050,0.0050,0.0050,0.0200/),r8)  !woody
  CPOFC(1:jskenc,this%k_woody_litr) = real((/0.0005,0.0005,0.0005,0.0020/),r8)  !woody
  CNOFC(1:jskenc,this%k_fine_litr)  = real((/0.0200,0.0200,0.0200,0.0200/),r8)  !non-woody
  CPOFC(1:jskenc,this%k_fine_litr)  = real((/0.0020,0.0020,0.0020,0.0020/),r8)  !non-woody
  CNOFC(1:jskenc,this%k_manure)     = real((/0.0200,0.0200,0.0200,0.0200/),r8)   !manure
  CPOFC(1:jskenc,this%k_manure)     = real((/0.0020,0.0020,0.0020,0.0020/),r8)   !manure
  FL(1:2)=real((/0.55,0.45/),r8)

  D95: DO K=1,this%jcplx
    DO  N=1,this%NumMicbFunGrupsPerCmplx
      IF(N.EQ.this%mid_Aerob_Fungi)THEN
        DO NGL=this%JGnio(n),this%JGnfo(n)
          rNCOMC(1,NGL,K)=0.15_r8           !maximum
          rNCOMC(2,NGL,K)=0.09_r8           !minimum
          rPCOMC(1,NGL,K)=0.015_r8
          rPCOMC(2,NGL,K)=0.009_r8
        ENDDO
        this%rNCOMCa(1,N,K)=0.15_r8           !maximum
        this%rNCOMCa(2,N,K)=0.09_r8           !minimum
        this%rPCOMCa(1,N,K)=0.015_r8
        this%rPCOMCa(2,N,K)=0.009_r8
      ELSE
        do NGL=this%JGnio(n),this%JGnfo(n)
          rNCOMC(1,NGL,K)=0.225_r8
          rNCOMC(2,NGL,K)=0.135_r8
          rPCOMC(1,NGL,K)=0.0225_r8
          rPCOMC(2,NGL,K)=0.0135_r8
        enddo
        this%rNCOMCa(1,N,K)=0.225_r8
        this%rNCOMCa(2,N,K)=0.135_r8
        this%rPCOMCa(1,N,K)=0.0225_r8
        this%rPCOMCa(2,N,K)=0.0135_r8
      ENDIF
      do NGL=this%JGnio(n),this%JGnfo(n)
        rNCOMC(3,NGL,K)=DOT_PRODUCT(FL,rNCOMC(1:2,NGL,K))
        rPCOMC(3,NGL,K)=DOT_PRODUCT(FL,rPCOMC(1:2,NGL,K))
      enddo
      this%rNCOMCa(3,N,K)=DOT_PRODUCT(FL,this%rNCOMCa(1:2,N,K))
      this%rPCOMCa(3,N,K)=DOT_PRODUCT(FL,this%rPCOMCa(1:2,N,K))
    enddo
  ENDDO D95
  DO  N=1,this%NumMicbFunGrupsPerCmplx
    do NGL=this%JGniA(n),this%JGnfA(n)
      rNCOMCAutor(1,NGL)=0.225_r8
      rNCOMCAutor(2,NGL)=0.135_r8
      rPCOMCAutor(1,NGL)=0.0225_r8
      rPCOMCAutor(2,NGL)=0.0135_r8
    enddo
    this%rNCOMCAutora(1,N)=0.225_r8
    this%rNCOMCAutora(2,N)=0.135_r8
    this%rPCOMCAutora(1,N)=0.0225_r8
    this%rPCOMCAutora(2,N)=0.0135_r8
    do NGL=this%JGniA(n),this%JGnfA(n)
      rNCOMCAutor(3,NGL)=DOT_PRODUCT(FL,rNCOMCAutor(1:2,NGL))
      rPCOMCAutor(3,NGL)=DOT_PRODUCT(FL,rPCOMCAutor(1:2,NGL))
    enddo
    this%rNCOMCAutora(3,N)=DOT_PRODUCT(FL,this%rNCOMCAutora(1:2,N))
    this%rPCOMCAutora(3,N)=DOT_PRODUCT(FL,this%rPCOMCAutora(1:2,N))
  enddo

  end associate
  end subroutine SetPars
!------------------------------------------------------------------------------------------

  subroutine Initallocate(this)
  implicit none
  class(MicParType) :: this

  integer :: jguilds
  integer :: NumMicbFunGrupsPerCmplx
  integer :: jcplx
  integer :: jsken
  integer :: n, k

  jguilds = this%jguilds
  NumMicbFunGrupsPerCmplx  =this%NumMicbFunGrupsPerCmplx
  jcplx =this%jcplx
  jsken =this%jsken
  allocate(this%JGnio(NumMicbFunGrupsPerCmplx))
  allocate(this%JGnfo(NumMicbFunGrupsPerCmplx))

  allocate(this%JGniA(NumMicbFunGrupsPerCmplx))
  allocate(this%JGnfA(NumMicbFunGrupsPerCmplx))

  k=1
  this%NumMicrobAutrophCmplx=0
  this%NumHetetrMicCmplx=0
  !replace the functional group specification with external input later
  do n=1,NumMicbFunGrupsPerCmplx
    this%JGnio(n)              = k
    this%JGniA(n)              = k
    k                          = k+jguilds
    this%JGnfo(n)              = k-1
    this%JGnfA(n)              = k-1
    this%NumMicrobAutrophCmplx = this%NumMicrobAutrophCmplx+this%JGnfA(n)-this%JGniA(n)+1
    this%NumHetetrMicCmplx     = this%NumHetetrMicCmplx+this%JGnfo(n)-this%JGnio(n)+1
  enddo
  this%NumLiveHeterBioms=this%nlbiomcp*this%NumHetetrMicCmplx
  this%NumLiveAutoBioms=this%nlbiomcp*this%NumMicrobAutrophCmplx

  allocate(this%DOSA(1:jcplx))
  allocate(this%SPOSC(jsken,1:jcplx))
  allocate(this%OHCK(1:jcplx))
  allocate(this%OMCK(1:jcplx))
  allocate(this%ORCK(1:jcplx))
  allocate(this%OQCK(1:jcplx))
  allocate(this%ORCI(NumDeadMicrbCompts,1:jcplx))
  allocate(this%OMCI(NumLiveMicrbCompts,1:jcplx))
  allocate(this%rNCOMC(NumLiveMicrbCompts,this%NumHetetrMicCmplx,1:jcplx))
  allocate(this%rPCOMC(NumLiveMicrbCompts,this%NumHetetrMicCmplx,1:jcplx))
  allocate(this%rNCOMCAutor(NumLiveMicrbCompts,this%NumMicrobAutrophCmplx))
  allocate(this%rPCOMCAutor(NumLiveMicrbCompts,this%NumMicrobAutrophCmplx))
  allocate(this%rNCOMCa(NumLiveMicrbCompts,NumMicbFunGrupsPerCmplx,1:jcplx))
  allocate(this%rPCOMCa(NumLiveMicrbCompts,NumMicbFunGrupsPerCmplx,1:jcplx))
  allocate(this%rNCOMCAutora(NumLiveMicrbCompts,NumMicbFunGrupsPerCmplx))
  allocate(this%rPCOMCAutora(NumLiveMicrbCompts,NumMicbFunGrupsPerCmplx))

  allocate(this%CNOFC(jsken,1:this%NumOfLitrCmplxs))
  allocate(this%CPOFC(jsken,1:this%NumOfLitrCmplxs))
  allocate(this%CNRH(1:jcplx))
  allocate(this%CPRH(1:jcplx))
  allocate(this%OMCF(NumMicbFunGrupsPerCmplx))
  allocate(this%OMCA(NumMicbFunGrupsPerCmplx))
  allocate(this%FL(2))
  allocate(this%is_activeMicrbFungrpAutor(NumMicbFunGrupsPerCmplx)); this%is_activeMicrbFungrpAutor=.false.
  allocate(this%is_CO2_autotroph(NumMicbFunGrupsPerCmplx)); this%is_CO2_autotroph=.false.
  allocate(this%is_aerobic_hetr(NumMicbFunGrupsPerCmplx)); this%is_aerobic_hetr=.false.
  allocate(this%is_anaerobic_hetr(NumMicbFunGrupsPerCmplx));this%is_anaerobic_hetr=.false.
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
  call destroy(this%rNCOMC)
  call destroy(this%rPCOMC)
  call destroy(this%rNCOMCAutor)
  call destroy(this%rPCOMCAutor)
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
  call destroy(this%is_activeMicrbFungrpAutor)
  call destroy(this%is_CO2_autotroph)
  call destroy(this%is_aerobic_hetr)
  call destroy(this%is_anaerobic_hetr)
  end subroutine DestructMicBGCPar

!------------------------------------------------------------------------------------------
  pure function get_micb_id(this,M,NGL)result(id)
  !return the id of biomass component M for
  !hetetroph guild NGL
  !1,2,3: labile, recalcitrant, reserve
  implicit none
  class(MicParType), intent(in) :: this
  integer, intent(in) :: M   !biomass component
  integer, intent(in) :: NGL !guild id
  integer :: id

  id=this%nlbiomcp*(NGL-1)+M

  end function get_micb_id

end module MicBGCPars
