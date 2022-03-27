module NitroDiagTypes
  ! !PUBLIC TYPES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  save
  private

! accumulative flux diagnostics
type, public :: NitroAQMFluxDiagType
! ratios
    real(r8) :: TFNH4B
    real(r8) :: TFNO3B
    real(r8) :: TFNO2B
    real(r8) :: TFP14B
    real(r8) :: TFPO4B
    real(r8) :: TCH4H
    real(r8) :: TCH4A
    real(r8) :: TFOQC
    real(r8) :: TFOQA
    real(r8) :: TFOXYX
    real(r8) :: TFNH4X
    real(r8) :: TFNO3X
    real(r8) :: TFNO2X
    real(r8) :: TFN2OX
    real(r8) :: TFP14X
    real(r8) :: TFPO4X
!fluxes
    real(r8) :: TRH2G
    real(r8) :: TRINH
    real(r8) :: TRIPO
    real(r8) :: TRINO
    real(r8) :: TRIP1
    real(r8) :: TRINB
    real(r8) :: TRIOB
    real(r8) :: TRIPB
    real(r8) :: TRIB1
    real(r8) :: TRGOM
    real(r8) :: TRGOC
    real(r8) :: TRGOD
    real(r8) :: TRGOA
    real(r8) :: TRGOH
    real(r8) :: TUPOX
    real(r8) :: TRDN3
    real(r8) :: TRDNB
    real(r8) :: TRDN2
    real(r8) :: TRD2B
    real(r8) :: TRDNO
    real(r8) :: TRN2F
  contains
    procedure, public :: ZeroOut => nit_aqmf_diag
  end type NitroAQMFluxDiagType

  type, public :: NitroMicStateType
  real(r8),allocatable :: CNOMA(:,:,:)
  real(r8),allocatable :: CPOMA(:,:,:)
  real(r8),allocatable :: OMA(:,:,:)

  real(r8),allocatable :: OMC2(:,:,:)
  real(r8),allocatable :: TFNG(:,:,:)
  real(r8),allocatable :: TFNR(:,:,:)
  real(r8),allocatable :: OMN2(:,:,:)
  real(r8),allocatable :: FOM2(:,:,:)
  real(r8),allocatable :: WFN(:,:,:)

  real(r8),allocatable :: FCN(:,:,:)
  real(r8),allocatable :: FCP(:,:,:)
  real(r8),allocatable :: FCNP(:,:,:)
  real(r8),allocatable :: FOMA(:,:,:)
  real(r8),allocatable :: FOMN(:,:,:)
  real(r8),allocatable :: FOMK(:,:,:)
  contains
   procedure, public :: Init => nit_mics_init
   procedure, public :: Destroy => nit_mics_destroy
  end type NitroMicStateType

  type, public :: NitroMicFluxType
! flux ratios
  real(r8),allocatable :: FNH4XR(:,:,:)
  real(r8),allocatable :: FNO3XR(:,:,:)
  real(r8),allocatable :: FPO4XR(:,:,:)
  real(r8),allocatable :: FP14XR(:,:,:)
! fluxes
  real(r8),allocatable :: CGOMC(:,:,:)
  real(r8),allocatable :: CGOMN(:,:,:)
  real(r8),allocatable :: CGOMP(:,:,:)
  real(r8),allocatable :: CGOQC(:,:,:)
  real(r8),allocatable :: CGOAC(:,:,:)
  real(r8),allocatable :: CGOMS(:,:,:,:)
  real(r8),allocatable :: CGONS(:,:,:,:)
  real(r8),allocatable :: CGOPS(:,:,:,:)
  real(r8),allocatable :: RUPOX(:,:,:)
  real(r8),allocatable :: RGN2F(:,:,:)
  real(r8),allocatable :: RGOMO(:,:,:)
  real(r8),allocatable :: ROXYM(:,:,:)
  real(r8),allocatable :: ROXYP(:,:,:)
  real(r8),allocatable :: ROXYO(:,:,:)
  real(r8),allocatable :: RDNO3(:,:,:)
  real(r8),allocatable :: RDNOB(:,:,:)
  real(r8),allocatable :: RDNO2(:,:,:)
  real(r8),allocatable :: RDN2B(:,:,:)
  real(r8),allocatable :: RDN2O(:,:,:)
  real(r8),allocatable :: RGOMD(:,:,:)
  real(r8),allocatable :: RMOMC(:,:,:,:)
  real(r8),allocatable :: RINH4(:,:,:)
  real(r8),allocatable :: RINO3(:,:,:)
  real(r8),allocatable :: RIPO4(:,:,:)
  real(r8),allocatable :: RINB4(:,:,:)
  real(r8),allocatable :: RINB3(:,:,:)
  real(r8),allocatable :: RIPOB(:,:,:)
  real(r8),allocatable :: RDOMC(:,:,:,:)
  real(r8),allocatable :: RDOMN(:,:,:,:)
  real(r8),allocatable :: RDOMP(:,:,:,:)
  real(r8),allocatable :: RHOMC(:,:,:,:)
  real(r8),allocatable :: RHOMN(:,:,:,:)
  real(r8),allocatable :: RHOMP(:,:,:,:)
  real(r8),allocatable :: RCOMC(:,:,:,:)
  real(r8),allocatable :: RCOMN(:,:,:,:)
  real(r8),allocatable :: RCOMP(:,:,:,:)
  real(r8),allocatable :: RH2GX(:,:,:)
  real(r8),allocatable :: RDMMC(:,:,:,:)
  real(r8),allocatable :: RHMMC(:,:,:,:)
  real(r8),allocatable :: RCMMC(:,:,:,:)
  real(r8),allocatable :: RDMMN(:,:,:,:)
  real(r8),allocatable :: RHMMN(:,:,:,:)
  real(r8),allocatable :: RCMMN(:,:,:,:)
  real(r8),allocatable :: RDMMP(:,:,:,:)
  real(r8),allocatable :: RHMMP(:,:,:,:)
  real(r8),allocatable :: RCMMP(:,:,:,:)
  real(r8),allocatable :: RCCMC(:,:,:,:)
  real(r8),allocatable :: RCCMN(:,:,:,:)
  real(r8),allocatable :: RCCMP(:,:,:,:)
  real(r8),allocatable :: RN2FX(:,:,:)
  real(r8),allocatable :: RXOMC(:,:,:,:)
  real(r8),allocatable :: RXOMN(:,:,:,:)
  real(r8),allocatable :: RXOMP(:,:,:,:)
  real(r8),allocatable :: R3OMC(:,:,:,:)
  real(r8),allocatable :: R3OMN(:,:,:,:)
  real(r8),allocatable :: R3OMP(:,:,:,:)
  real(r8),allocatable :: RXMMC(:,:,:,:)
  real(r8),allocatable :: RXMMN(:,:,:,:)
  real(r8),allocatable :: RXMMP(:,:,:,:)
  real(r8),allocatable :: R3MMC(:,:,:,:)
  real(r8),allocatable :: R3MMN(:,:,:,:)
  real(r8),allocatable :: R3MMP(:,:,:,:)
  real(r8),allocatable :: RINH4R(:,:,:)
  real(r8),allocatable :: RINO3R(:,:,:)
  real(r8),allocatable :: RIPO4R(:,:,:)
  real(r8),allocatable :: RGOMY(:,:,:)
  real(r8),allocatable :: ROQCD(:,:,:)
  real(r8),allocatable :: RCO2X(:,:,:)
  real(r8),allocatable :: RCH3X(:,:,:)
  real(r8),allocatable :: RCH4X(:,:,:)
  real(r8),allocatable :: RVOXA(:,:)
  real(r8),allocatable :: RVOXB(:,:)
  real(r8),allocatable :: XOMCZ(:,:,:,:)
  real(r8),allocatable :: XOMNZ(:,:,:,:)
  real(r8),allocatable :: XOMPZ(:,:,:,:)
  real(r8),allocatable :: RIP14(:,:,:)
  real(r8),allocatable :: RIP1B(:,:,:)
  real(r8),allocatable :: RIP14R(:,:,:)
  contains
    procedure, public :: Init => nit_micf_init
    procedure, public :: ZeroOut => nit_micf_zero
    procedure, public :: destroy => nit_micf_destroy
  end type NitroMicFluxType

  type, public :: NitroOMcplxFluxType
    real(r8),allocatable :: RDOSC(:,:)
    real(r8),allocatable :: RDOSN(:,:)
    real(r8),allocatable :: RDOSP(:,:)
    real(r8),allocatable :: RHOSC(:,:)
    real(r8),allocatable :: RHOSN(:,:)
    real(r8),allocatable :: RHOSP(:,:)
    real(r8),allocatable :: RCOSC(:,:)
    real(r8),allocatable :: RCOSN(:,:)
    real(r8),allocatable :: RCOSP(:,:)
    real(r8),allocatable :: RDORC(:,:)
    real(r8),allocatable :: RDORN(:,:)
    real(r8),allocatable :: RDORP(:,:)
    real(r8),allocatable :: RDOHC(:)
    real(r8),allocatable :: RDOHN(:)
    real(r8),allocatable :: RDOHP(:)
    real(r8),allocatable :: RDOHA(:)
    real(r8),allocatable :: CSORP(:)
    real(r8),allocatable :: ZSORP(:)
    real(r8),allocatable :: PSORP(:)
    real(r8),allocatable :: CSORPA(:)
    real(r8),allocatable :: TCGOQC(:)
    real(r8),allocatable :: TCGOAC(:)
    real(r8),allocatable :: TCGOMN(:)
    real(r8),allocatable :: TCGOMP(:)
    real(r8),allocatable :: ROQCK(:)
    real(r8),allocatable :: XOQCK(:)
    real(r8),allocatable :: XOQCZ(:)
    real(r8),allocatable :: XOQNZ(:)
    real(r8),allocatable :: XOQPZ(:)
    real(r8),allocatable :: XOQAZ(:)
  contains
    procedure, public :: Init => nit_omcplxf_init
    procedure, public :: ZeroOut => nit_omcplxf_zero
    procedure, public :: Destroy => nit_omcplxf_destroy
  end type NitroOMcplxFluxType

  type, public :: NitroOMcplxStateType
    real(r8),allocatable :: OSRH(:)
    real(r8),allocatable :: TOMK(:)
    real(r8),allocatable :: TONK(:)
    real(r8),allocatable :: TOPK(:)
    real(r8),allocatable :: FOCA(:)
    real(r8),allocatable :: FOAA(:)
    real(r8),allocatable :: CNQ(:)
    real(r8),allocatable :: CPQ(:)
    real(r8),allocatable :: CNH(:)
    real(r8),allocatable :: CPH(:)
    real(r8),allocatable :: ORCT(:)
    real(r8),allocatable :: OSCT(:)
    real(r8),allocatable :: OSAT(:)
    real(r8),allocatable :: TONX(:)
    real(r8),allocatable :: TOPX(:)
  contains
    procedure, public :: Init => nit_omcplxs_init
    procedure, public :: ZeroOut => nit_omcplxs_zero
    procedure, public :: Destroy => nit_omcplxs_destroy
  end type NitroOMcplxStateType

  type, public :: NitroMicDiagType
  real(r8) :: H1P4T
  real(r8) :: H2P4T
  real(r8) :: RH2GZ
  real(r8) :: RCNO2
  real(r8) :: RCNOB
  real(r8) :: RCN2O
  real(r8) :: RCN2B
  real(r8) :: RCNO3
  real(r8) :: RCN3B
  real(r8) :: RCOQN
  real(r8) :: THETR
  real(r8) :: THETZ
  real(r8) :: TORC
  real(r8) :: TOMA
  real(r8) :: TOMN
  real(r8) :: TFNX
  real(r8) :: TFNY
  real(r8) :: VOLWZ
  real(r8) :: WFNG
  real(r8) :: XCO2
  real(r8) :: ZNH4T
  real(r8) :: ZNO3T
  real(r8) :: ZNO2T

  end type NitroMicDiagType

  contains
!------------------------------------------------------------------------------------------

  subroutine nit_aqmf_diag(this)
  implicit none
  class(NitroAQMFluxDiagType) :: this

    this%TFNH4B = 0._r8
    this%TFNO3B = 0._r8
    this%TFNO2B = 0._r8
    this%TFP14B = 0._r8
    this%TFPO4B = 0._r8
    this%TCH4H = 0._r8
    this%TCH4A = 0._r8
    this%TFOQC = 0._r8
    this%TFOQA = 0._r8
    this%TFOXYX = 0._r8
    this%TFNH4X = 0._r8
    this%TFNO3X = 0._r8
    this%TFNO2X = 0._r8
    this%TFN2OX = 0._r8
    this%TFP14X = 0._r8
    this%TFPO4X = 0._r8

    this%TRH2G = 0._r8
    this%TRINH = 0._r8
    this%TRIPO = 0._r8
    this%TRINO=0.0_r8
    this%TRIPO=0.0_r8
    this%TRIP1=0.0_r8
    this%TRINB=0.0_r8
    this%TRIOB=0.0_r8
    this%TRIPB=0.0_r8
    this%TRIB1=0.0_r8
    this%TRGOM=0.0_r8
    this%TRGOC=0.0_r8
    this%TRGOD=0.0_r8
    this%TRGOA=0.0_r8
    this%TRGOH=0.0_r8
    this%TUPOX=0.0_r8
    this%TRDN3=0.0_r8
    this%TRDNB=0.0_r8
    this%TRDN2=0.0_r8
    this%TRD2B=0.0_r8
    this%TRDNO=0.0_r8
    this%TRN2F=0.0_r8
  end subroutine nit_aqmf_diag
!------------------------------------------------------------------------------------------

  subroutine nit_micf_init(this,JG,jcplx)
  implicit none
  class(NitroMicFluxType) :: this
  integer, intent(in) :: JG,jcplx
  integer :: jcplx1
  jcplx1=jcplx-1
  allocate(this%RUPOX(JG,7,0:jcplx))
  allocate(this%RGN2F(JG,7,0:jcplx))
  allocate(this%RGOMO(JG,7,0:jcplx))
  allocate(this%ROXYM(JG,7,0:jcplx))
  allocate(this%ROXYP(JG,7,0:jcplx))
  allocate(this%ROXYO(JG,7,0:jcplx))
  allocate(this%RDNO3(JG,7,0:jcplx))
  allocate(this%RDNOB(JG,7,0:jcplx))
  allocate(this%RDNO2(JG,7,0:jcplx))
  allocate(this%RDN2B(JG,7,0:jcplx))
  allocate(this%RDN2O(JG,7,0:jcplx))
  allocate(this%RGOMD(JG,7,0:jcplx))
  allocate(this%RMOMC(2,JG,7,0:jcplx))
  allocate(this%RINH4(JG,7,0:jcplx))
  allocate(this%RINO3(JG,7,0:jcplx))
  allocate(this%RIPO4(JG,7,0:jcplx))
  allocate(this%RINB4(JG,7,0:jcplx))
  allocate(this%RINB3(JG,7,0:jcplx))
  allocate(this%RIPOB(JG,7,0:jcplx))
  allocate(this%RDOMC(2,JG,7,0:jcplx))
  allocate(this%RDOMN(2,JG,7,0:jcplx))
  allocate(this%RDOMP(2,JG,7,0:jcplx))
  allocate(this%RHOMC(2,JG,7,0:jcplx))
  allocate(this%RHOMN(2,JG,7,0:jcplx))
  allocate(this%RHOMP(2,JG,7,0:jcplx))
  allocate(this%RCOMC(2,JG,7,0:jcplx))
  allocate(this%RCOMN(2,JG,7,0:jcplx))
  allocate(this%RCOMP(2,JG,7,0:jcplx))
  allocate(this%CGOMC(JG,7,0:jcplx))
  allocate(this%CGOMN(JG,7,0:jcplx))
  allocate(this%RH2GX(JG,7,0:jcplx))
  allocate(this%CGOMP(JG,7,0:jcplx))
  allocate(this%RDMMC(2,JG,7,0:jcplx))
  allocate(this%RHMMC(2,JG,7,0:jcplx))
  allocate(this%RCMMC(2,JG,7,0:jcplx))
  allocate(this%RDMMN(2,JG,7,0:jcplx))
  allocate(this%RHMMN(2,JG,7,0:jcplx))
  allocate(this%RCMMN(2,JG,7,0:jcplx))
  allocate(this%RDMMP(2,JG,7,0:jcplx))
  allocate(this%RHMMP(2,JG,7,0:jcplx))
  allocate(this%RCMMP(2,JG,7,0:jcplx))
  allocate(this%RCCMC(2,JG,7,0:jcplx1))
  allocate(this%RCCMN(2,JG,7,0:jcplx1))
  allocate(this%RCCMP(2,JG,7,0:jcplx1))
  allocate(this%RN2FX(JG,7,0:jcplx))

  allocate(this%RXOMC(2,JG,7,0:jcplx))
  allocate(this%RXOMN(2,JG,7,0:jcplx))
  allocate(this%RXOMP(2,JG,7,0:jcplx))
  allocate(this%R3OMC(2,JG,7,0:jcplx))
  allocate(this%R3OMN(2,JG,7,0:jcplx))
  allocate(this%R3OMP(2,JG,7,0:jcplx))
  allocate(this%RXMMC(2,JG,7,0:jcplx))
  allocate(this%RXMMN(2,JG,7,0:jcplx))
  allocate(this%RXMMP(2,JG,7,0:jcplx))
  allocate(this%R3MMC(2,JG,7,0:jcplx))
  allocate(this%R3MMN(2,JG,7,0:jcplx))
  allocate(this%R3MMP(2,JG,7,0:jcplx))
  allocate(this%CGOQC(JG,7,0:jcplx))
  allocate(this%CGOAC(JG,7,0:jcplx))
  allocate(this%RINH4R(JG,7,0:jcplx))
  allocate(this%RINO3R(JG,7,0:jcplx))
  allocate(this%RIPO4R(JG,7,0:jcplx))
  allocate(this%FNH4XR(JG,7,0:jcplx))
  allocate(this%FNO3XR(JG,7,0:jcplx))
  allocate(this%FPO4XR(JG,7,0:jcplx))
  allocate(this%RGOMY(JG,7,0:jcplx))
  allocate(this%ROQCD(JG,7,0:jcplx1))
  allocate(this%CGOMS(2,JG,7,0:jcplx))
  allocate(this%CGONS(2,JG,7,0:jcplx))
  allocate(this%CGOPS(2,JG,7,0:jcplx))
  allocate(this%FP14XR(JG,7,0:jcplx))
  allocate(this%RCO2X(JG,7,0:jcplx))
  allocate(this%RCH3X(JG,7,0:jcplx))
  allocate(this%RCH4X(JG,7,0:jcplx))
  allocate(this%RVOXA(JG,7))
  allocate(this%RVOXB(JG,7))
  allocate(this%XOMCZ(3,JG,7,0:jcplx1))
  allocate(this%XOMNZ(3,JG,7,0:jcplx1))
  allocate(this%XOMPZ(3,JG,7,0:jcplx1))
  allocate(this%RIP14(JG,7,0:jcplx))
  allocate(this%RIP1B(JG,7,0:jcplx))
  allocate(this%RIP14R(JG,7,0:jcplx))

  call this%ZeroOut()
  end subroutine nit_micf_init
!------------------------------------------------------------------------------------------

  subroutine nit_mics_init(this, JG, jcplx)

  implicit none
  class(NitroMicStateType) :: this
  integer, intent(in) :: JG, jcplx

  allocate(this%CNOMA(JG,7,0:jcplx))
  allocate(this%CPOMA(JG,7,0:jcplx))
  allocate(this%OMA(JG,7,0:jcplx))
  allocate(this%FOMA(JG,7,0:jcplx))
  allocate(this%FOMN(JG,7,0:jcplx))
  allocate(this%FOMK(JG,7,0:jcplx))
  allocate(this%OMC2(JG,7,0:jcplx))
  allocate(this%TFNG(JG,7,0:jcplx))
  allocate(this%TFNR(JG,7,0:jcplx))
  allocate(this%OMN2(JG,7,0:jcplx))
  allocate(this%FOM2(JG,7,0:jcplx))
  allocate(this%WFN(JG,7,0:jcplx))
  allocate(this%FCN(JG,7,0:jcplx))
  allocate(this%FCP(JG,7,0:jcplx))
  allocate(this%FCNP(JG,7,0:jcplx))
  end subroutine nit_mics_init
!------------------------------------------------------------------------------------------

  subroutine nit_micf_zero(this)
  implicit none
  class(NitroMicFLuxType) :: this

  this%RUPOX = 0._r8
  this%RGN2F = 0._r8
  this%RGOMO = 0._r8
  this%ROXYM = 0._r8
  this%ROXYP = 0._r8
  this%ROXYO = 0._r8
  this%RDNO3 = 0._r8
  this%RDNOB = 0._r8
  this%RDNO2 = 0._r8
  this%RDN2B = 0._r8
  this%RDN2O = 0._r8
  this%RGOMD = 0._r8
  this%RMOMC = 0._r8
  this%RINH4 = 0._r8
  this%RINO3 = 0._r8
  this%RIPO4 = 0._r8
  this%RINB4 = 0._r8
  this%RINB3 = 0._r8
  this%RIPOB = 0._r8
  this%RDOMC = 0._r8
  this%RDOMN = 0._r8
  this%RDOMP = 0._r8
  this%RHOMC = 0._r8
  this%RHOMN = 0._r8
  this%RHOMP = 0._r8
  this%RCOMC = 0._r8
  this%RCOMN = 0._r8
  this%RCOMP = 0._r8
  this%CGOMC = 0._r8
  this%CGOMN = 0._r8
  this%RH2GX = 0._r8
  this%CGOMP = 0._r8
  this%RDMMC = 0._r8
  this%RHMMC = 0._r8
  this%RCMMC = 0._r8
  this%RDMMN = 0._r8
  this%RHMMN = 0._r8
  this%RCMMN = 0._r8
  this%RDMMP = 0._r8
  this%RHMMP = 0._r8
  this%RCMMP = 0._r8
  this%RCCMC = 0._r8
  this%RCCMN = 0._r8
  this%RCCMP = 0._r8
  this%RN2FX = 0._r8

  this%RXOMC = 0._r8
  this%RXOMN = 0._r8
  this%RXOMP = 0._r8
  this%R3OMC = 0._r8
  this%R3OMN = 0._r8
  this%R3OMP = 0._r8
  this%RXMMC = 0._r8
  this%RXMMN = 0._r8
  this%RXMMP = 0._r8
  this%R3MMC = 0._r8
  this%R3MMN = 0._r8
  this%R3MMP = 0._r8
  this%CGOQC = 0._r8
  this%CGOAC = 0._r8
  this%RINH4R = 0._r8
  this%RINO3R = 0._r8
  this%RIPO4R = 0._r8
  this%FNH4XR = 0._r8
  this%FNO3XR = 0._r8
  this%FPO4XR = 0._r8
  this%RGOMY = 0._r8
  this%ROQCD = 0._r8
  this%CGOMS = 0._r8
  this%CGONS = 0._r8
  this%CGOPS = 0._r8
  this%FP14XR = 0._r8
  this%RCO2X = 0._r8
  this%RCH3X = 0._r8
  this%RCH4X = 0._r8
  this%RVOXA = 0._r8
  this%RVOXB = 0._r8
  this%XOMCZ = 0._r8
  this%XOMNZ = 0._r8
  this%XOMPZ = 0._r8
  this%RIP14 = 0._r8
  this%RIP1B = 0._r8
  this%RIP14R = 0._r8

  end subroutine nit_micf_zero


!------------------------------------------------------------------------------------------

  subroutine nit_micf_destroy(this)

  implicit none
  class(NitroMicFLuxType) :: this

  if(allocated(this%RUPOX))deallocate(this%RUPOX)
  if(allocated(this%RGN2F))deallocate(this%RGN2F)
  if(allocated(this%RGOMO))deallocate(this%RGOMO)
  if(allocated(this%ROXYM))deallocate(this%ROXYM)
  if(allocated(this%ROXYP))deallocate(this%ROXYP)
  if(allocated(this%ROXYO))deallocate(this%ROXYO)
  if(allocated(this%RDNO3))deallocate(this%RDNO3)
  if(allocated(this%RDNOB))deallocate(this%RDNOB)
  if(allocated(this%RDNO2))deallocate(this%RDNO2)
  if(allocated(this%RDN2B))deallocate(this%RDN2B)
  if(allocated(this%RDN2O))deallocate(this%RDN2O)
  if(allocated(this%RGOMD))deallocate(this%RGOMD)
  if(allocated(this%RMOMC))deallocate(this%RMOMC)
  if(allocated(this%RINH4))deallocate(this%RINH4)
  if(allocated(this%RINO3))deallocate(this%RINO3)
  if(allocated(this%RIPO4))deallocate(this%RIPO4)
  if(allocated(this%RINB4))deallocate(this%RINB4)
  if(allocated(this%RINB3))deallocate(this%RINB3)
  if(allocated(this%RIPOB))deallocate(this%RIPOB)
  if(allocated(this%RDOMC))deallocate(this%RDOMC)
  if(allocated(this%RDOMN))deallocate(this%RDOMN)
  if(allocated(this%RDOMP))deallocate(this%RDOMP)
  if(allocated(this%RHOMC))deallocate(this%RHOMC)
  if(allocated(this%RHOMN))deallocate(this%RHOMN)
  if(allocated(this%RHOMP))deallocate(this%RHOMP)
  if(allocated(this%RCOMC))deallocate(this%RCOMC)
  if(allocated(this%RCOMN))deallocate(this%RCOMN)
  if(allocated(this%RCOMP))deallocate(this%RCOMP)
  if(allocated(this%CGOMC))deallocate(this%CGOMC)
  if(allocated(this%CGOMN))deallocate(this%CGOMN)
  if(allocated(this%RH2GX))deallocate(this%RH2GX)
  if(allocated(this%CGOMP))deallocate(this%CGOMP)
  if(allocated(this%RDMMC))deallocate(this%RDMMC)
  if(allocated(this%RHMMC))deallocate(this%RHMMC)
  if(allocated(this%RCMMC))deallocate(this%RCMMC)
  if(allocated(this%RDMMN))deallocate(this%RDMMN)
  if(allocated(this%RHMMN))deallocate(this%RHMMN)
  if(allocated(this%RCMMN))deallocate(this%RCMMN)
  if(allocated(this%RDMMP))deallocate(this%RDMMP)
  if(allocated(this%RHMMP))deallocate(this%RHMMP)
  if(allocated(this%RCMMP))deallocate(this%RCMMP)
  if(allocated(this%RCCMC))deallocate(this%RCCMC)
  if(allocated(this%RCCMN))deallocate(this%RCCMN)
  if(allocated(this%RCCMP))deallocate(this%RCCMP)
  if(allocated(this%RN2FX))deallocate(this%RN2FX)
  if(allocated(this%RXOMC))deallocate(this%RXOMC)
  if(allocated(this%RXOMN))deallocate(this%RXOMN)
  if(allocated(this%RXOMP))deallocate(this%RXOMP)
  if(allocated(this%R3OMC))deallocate(this%R3OMC)
  if(allocated(this%R3OMN))deallocate(this%R3OMN)
  if(allocated(this%R3OMP))deallocate(this%R3OMP)
  if(allocated(this%RXMMC))deallocate(this%RXMMC)
  if(allocated(this%RXMMN))deallocate(this%RXMMN)
  if(allocated(this%RXMMP))deallocate(this%RXMMP)
  if(allocated(this%R3MMC))deallocate(this%R3MMC)
  if(allocated(this%R3MMN))deallocate(this%R3MMN)
  if(allocated(this%R3MMP))deallocate(this%R3MMP)
  if(allocated(this%CGOQC))deallocate(this%CGOQC)
  if(allocated(this%CGOAC))deallocate(this%CGOAC)
  if(allocated(this%RINH4R))deallocate(this%RINH4R)
  if(allocated(this%RINO3R))deallocate(this%RINO3R)
  if(allocated(this%RIPO4R))deallocate(this%RIPO4R)
  if(allocated(this%FNH4XR))deallocate(this%FNH4XR)
  if(allocated(this%FNO3XR))deallocate(this%FNO3XR)
  if(allocated(this%FPO4XR))deallocate(this%FPO4XR)
  if(allocated(this%RGOMY))deallocate(this%RGOMY)
  if(allocated(this%ROQCD))deallocate(this%ROQCD)
  if(allocated(this%CGOMS))deallocate(this%CGOMS)
  if(allocated(this%CGONS))deallocate(this%CGONS)
  if(allocated(this%CGOPS))deallocate(this%CGOPS)
  if(allocated(this%FP14XR))deallocate(this%FP14XR)
  if(allocated(this%RCO2X))deallocate(this%RCO2X)
  if(allocated(this%RCH3X))deallocate(this%RCH3X)
  if(allocated(this%RCH4X))deallocate(this%RCH4X)
  if(allocated(this%RVOXA))deallocate(this%RVOXA)
  if(allocated(this%RVOXB))deallocate(this%RVOXB)
  if(allocated(this%XOMCZ))deallocate(this%XOMCZ)
  if(allocated(this%XOMNZ))deallocate(this%XOMNZ)
  if(allocated(this%XOMPZ))deallocate(this%XOMPZ)
  if(allocated(this%RIP14))deallocate(this%RIP14)
  if(allocated(this%RIP1B))deallocate(this%RIP1B)
  if(allocated(this%RIP14R))deallocate(this%RIP14R)

  end subroutine nit_micf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_mics_destroy(this)
  implicit none
  class(NitroMicStateType) :: this

  if(allocated(this%CNOMA))deallocate(this%CNOMA)
  if(allocated(this%CPOMA))deallocate(this%CPOMA)
  if(allocated(this%OMA))deallocate(this%OMA)
  if(allocated(this%FOMA))deallocate(this%FOMA)
  if(allocated(this%FOMN))deallocate(this%FOMN)
  if(allocated(this%FOMK))deallocate(this%FOMK)
  if(allocated(this%OMC2))deallocate(this%OMC2)
  if(allocated(this%TFNG))deallocate(this%TFNG)
  if(allocated(this%TFNR))deallocate(this%TFNR)
  if(allocated(this%OMN2))deallocate(this%OMN2)
  if(allocated(this%FOM2))deallocate(this%FOM2)
  if(allocated(this%WFN))deallocate(this%WFN)
  if(allocated(this%FCN))deallocate(this%FCN)
  if(allocated(this%FCP))deallocate(this%FCP)
  if(allocated(this%FCNP))deallocate(this%FCNP)

  end subroutine nit_mics_destroy
!------------------------------------------------------------------------------------------
  subroutine nit_omcplxf_init(this)
  implicit none
  class(NitroOMcplxFluxType) :: this
  integer, parameter :: nkinets=4
  integer, parameter :: ncplx =4

  allocate(this%RDOSC(nkinets,0:ncplx))
  allocate(this%RDOSN(nkinets,0:ncplx))
  allocate(this%RDOSP(nkinets,0:ncplx))
  allocate(this%RHOSC(nkinets,0:ncplx))
  allocate(this%RHOSN(nkinets,0:ncplx))
  allocate(this%RHOSP(nkinets,0:ncplx))
  allocate(this%RCOSC(nkinets,0:ncplx))
  allocate(this%RCOSN(nkinets,0:ncplx))
  allocate(this%RCOSP(nkinets,0:ncplx))
  allocate(this%RDORC(2,0:ncplx))
  allocate(this%RDORN(2,0:ncplx))
  allocate(this%RDORP(2,0:ncplx))
  allocate(this%RDOHC(0:ncplx))
  allocate(this%RDOHN(0:ncplx))
  allocate(this%RDOHP(0:ncplx))
  allocate(this%RDOHA(0:ncplx))
  allocate(this%CSORP(0:ncplx))
  allocate(this%ZSORP(0:ncplx))
  allocate(this%PSORP(0:ncplx))
  allocate(this%CSORPA(0:ncplx))
  allocate(this%TCGOQC(0:ncplx+1))
  allocate(this%TCGOAC(0:ncplx+1))
  allocate(this%TCGOMN(0:ncplx+1))
  allocate(this%TCGOMP(0:ncplx+1))
  allocate(this%ROQCK(0:ncplx))
  allocate(this%XOQCK(0:ncplx))
  allocate(this%XOQCZ(0:ncplx))
  allocate(this%XOQNZ(0:ncplx))
  allocate(this%XOQPZ(0:ncplx))
  allocate(this%XOQAZ(0:ncplx))

  call this%ZeroOut()
  end subroutine nit_omcplxf_init
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_destroy(this)
  implicit none
  class(NitroOMcplxFluxType) :: this

  if(allocated(this%RDOSC))deallocate(this%RDOSC)
  if(allocated(this%RDOSN))deallocate(this%RDOSN)
  if(allocated(this%RDOSP))deallocate(this%RDOSP)
  if(allocated(this%RHOSC))deallocate(this%RHOSC)
  if(allocated(this%RHOSN))deallocate(this%RHOSN)
  if(allocated(this%RHOSP))deallocate(this%RHOSP)
  if(allocated(this%RCOSC))deallocate(this%RCOSC)
  if(allocated(this%RCOSN))deallocate(this%RCOSN)
  if(allocated(this%RCOSP))deallocate(this%RCOSP)
  if(allocated(this%RDORC))deallocate(this%RDORC)
  if(allocated(this%RDORN))deallocate(this%RDORN)
  if(allocated(this%RDORP))deallocate(this%RDORP)
  if(allocated(this%RDOHC))deallocate(this%RDOHC)
  if(allocated(this%RDOHN))deallocate(this%RDOHN)
  if(allocated(this%RDOHP))deallocate(this%RDOHP)
  if(allocated(this%RDOHA))deallocate(this%RDOHA)
  if(allocated(this%CSORP))deallocate(this%CSORP)
  if(allocated(this%ZSORP))deallocate(this%ZSORP)
  if(allocated(this%PSORP))deallocate(this%PSORP)
  if(allocated(this%CSORPA))deallocate(this%CSORPA)
  if(allocated(this%TCGOQC))deallocate(this%TCGOQC)
  if(allocated(this%TCGOAC))deallocate(this%TCGOAC)
  if(allocated(this%TCGOMN))deallocate(this%TCGOMN)
  if(allocated(this%TCGOMP))deallocate(this%TCGOMP)
  if(allocated(this%ROQCK))deallocate(this%ROQCK)
  if(allocated(this%XOQCK))deallocate(this%XOQCK)
  if(allocated(this%XOQCZ))deallocate(this%XOQCZ)
  if(allocated(this%XOQNZ))deallocate(this%XOQNZ)
  if(allocated(this%XOQPZ))deallocate(this%XOQPZ)
  if(allocated(this%XOQAZ))deallocate(this%XOQAZ)

  end subroutine nit_omcplxf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_zero(this)
  implicit none
  class(NitroOMcplxFluxType) :: this

  this%RDOSC=0._r8
  this%RDOSN=0._r8
  this%RDOSP=0._r8
  this%RHOSC=0._r8
  this%RHOSN=0._r8
  this%RHOSP=0._r8
  this%RCOSC=0._r8
  this%RCOSN=0._r8
  this%RCOSP=0._r8
  this%RDORC=0._r8
  this%RDORN=0._r8
  this%RDORP=0._r8
  this%RDOHC=0._r8
  this%RDOHN=0._r8
  this%RDOHP=0._r8
  this%RDOHA=0._r8
  this%CSORP=0._r8
  this%ZSORP=0._r8
  this%PSORP=0._r8
  this%CSORPA=0._r8
  this%TCGOQC=0._r8
  this%TCGOAC=0._r8
  this%TCGOMN=0._r8
  this%TCGOMP=0._r8
  this%ROQCK=0._r8
  this%XOQCK=0._r8
  this%XOQCZ=0._r8
  this%XOQNZ=0._r8
  this%XOQPZ=0._r8
  this%XOQAZ=0._r8
  end subroutine nit_omcplxf_zero
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxs_init(this)

  implicit none
  class(NitroOMcplxStateType) :: this

  integer, parameter :: ncplx=4

  allocate(this%OSRH(0:ncplx))
  allocate(this%TOMK(0:ncplx+1))
  allocate(this%TONK(0:ncplx+1))
  allocate(this%TOPK(0:ncplx+1))
  allocate(this%FOCA(0:ncplx))
  allocate(this%FOAA(0:ncplx))
  allocate(this%CNQ(0:ncplx))
  allocate(this%CPQ(0:ncplx))
  allocate(this%CNH(0:ncplx))
  allocate(this%CPH(0:ncplx))
  allocate(this%ORCT(0:ncplx))
  allocate(this%OSCT(0:ncplx))
  allocate(this%OSAT(0:ncplx))
  allocate(this%TONX(0:ncplx+1))
  allocate(this%TOPX(0:ncplx+1))

  call this%ZeroOut()
  end subroutine nit_omcplxs_init
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxs_zero(this)

  implicit none
  class(NitroOMcplxStateType) :: this

  this%OSRH=0._r8
  this%TOMK=0._r8
  this%TONK=0._r8
  this%TOPK=0._r8
  this%FOCA=0._r8
  this%FOAA=0._r8
  this%CNQ=0._r8
  this%CPQ=0._r8
  this%CNH=0._r8
  this%CPH=0._r8
  this%ORCT=0._r8
  this%OSCT=0._r8
  this%OSAT=0._r8
  this%TONX=0._r8
  this%TOPX=0._r8

  end subroutine nit_omcplxs_zero


!------------------------------------------------------------------------------------------
  subroutine nit_omcplxs_destroy(this)
  implicit none
  class(NitroOMcplxStateType) :: this

  if(allocated(this%OSRH))deallocate(this%OSRH)
  if(allocated(this%TOMK))deallocate(this%TOMK)
  if(allocated(this%TONK))deallocate(this%TONK)
  if(allocated(this%TOPK))deallocate(this%TOPK)
  if(allocated(this%FOCA))deallocate(this%FOCA)
  if(allocated(this%FOAA))deallocate(this%FOAA)
  if(allocated(this%CNQ))deallocate(this%CNQ)
  if(allocated(this%CPQ))deallocate(this%CPQ)
  if(allocated(this%CNH))deallocate(this%CNH)
  if(allocated(this%CPH))deallocate(this%CPH)
  if(allocated(this%ORCT))deallocate(this%ORCT)
  if(allocated(this%OSCT))deallocate(this%OSCT)
  if(allocated(this%OSAT))deallocate(this%OSAT)
  if(allocated(this%TONX))deallocate(this%TONX)
  if(allocated(this%TOPX))deallocate(this%TOPX)

  end subroutine nit_omcplxs_destroy
end module NitroDiagTypes
