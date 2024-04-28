module NitroDiagTypes
  ! !PUBLIC TYPES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL       
  USE abortutils, only : destroy
  use TracerIDMod
  use EcoSiMParDataMod, only : micpar
  implicit none
  save
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

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
  real(r8),allocatable :: CNOMActHeter(:,:)
  real(r8),allocatable :: CPOMActHeter(:,:)
  real(r8),allocatable :: OMActHeter(:,:)

  real(r8),allocatable :: OMC2(:,:)
  real(r8),allocatable :: TFNG(:,:)
  real(r8),allocatable :: TFNR(:,:)
  real(r8),allocatable :: OMN2(:,:)
  real(r8),allocatable :: FOM2(:,:)
  real(r8),allocatable :: WFN(:,:)

  real(r8),allocatable :: FCN(:,:)
  real(r8),allocatable :: FCP(:,:)
  real(r8),allocatable :: FCNP(:,:)
  real(r8),allocatable :: FracOMActHeter(:,:)
  real(r8),allocatable :: FOMN(:,:)
  real(r8),allocatable :: FOMK(:,:)

  real(r8),allocatable :: CNOMActAutor(:)
  real(r8),allocatable :: CPOMActAutor(:)
  real(r8),allocatable :: OMActAutor(:)
  real(r8),allocatable :: FracOMActAutor(:)
  real(r8),allocatable :: FOMNff(:)
  real(r8),allocatable :: FOMKff(:)
  real(r8),allocatable :: OMC2ff(:)
  real(r8),allocatable :: TFNGff(:)
  real(r8),allocatable :: TFNRff(:)
  real(r8),allocatable :: OMN2ff(:)
  real(r8),allocatable :: FOM2ff(:)
  real(r8),allocatable :: fLimO2Autor(:)
  real(r8),allocatable :: FCNff(:)
  real(r8),allocatable :: FCPff(:)
  real(r8),allocatable :: FCNPff(:)
  contains
   procedure, public :: Init => nit_mics_init
   procedure, public :: Destroy => nit_mics_destroy
  end type NitroMicStateType

  type, public :: NitroMicFluxType
! flux ratios
  real(r8),allocatable :: FNH4XR(:,:)
  real(r8),allocatable :: FNO3XR(:,:)
  real(r8),allocatable :: FPO4XR(:,:)
  real(r8),allocatable :: FP14XR(:,:)
! fluxes
  real(r8),allocatable :: CGOMEheter(:,:,:)
  real(r8),allocatable :: CGOQC(:,:)
  real(r8),allocatable :: CGOAC(:,:)
  real(r8),allocatable :: CGOMES(:,:,:,:)
  real(r8),allocatable :: RO2UptkHeter(:,:)
  real(r8),allocatable :: RGN2F(:,:)
  real(r8),allocatable :: RGOMO(:,:)
  real(r8),allocatable :: RO2Dmnd4RespHeter(:,:)
  real(r8),allocatable :: RO2DmndHeter(:,:)
  real(r8),allocatable :: RO2Uptk4RespHeter(:,:)
  real(r8),allocatable :: RNO3ReduxHeterSoil(:,:)
  real(r8),allocatable :: RNO3ReduxHeterBand(:,:)
  real(r8),allocatable :: RNO2ReduxHeterSoil(:,:)
  real(r8),allocatable :: RNO2ReduxHeterBand(:,:)
  real(r8),allocatable :: RN2OReduxHeter(:,:)
  real(r8),allocatable :: RGOMD(:,:)
  real(r8),allocatable :: RMOMC(:,:,:)
  real(r8),allocatable :: RNH4TransfSoilHeter(:,:)
  real(r8),allocatable :: RNO3TransfSoilHeter(:,:)
  real(r8),allocatable :: RH2PO4TransfSoilHeter(:,:)
  real(r8),allocatable :: RNH4TransfBandHeter(:,:)
  real(r8),allocatable :: RNO3TransfBandHeter(:,:)
  real(r8),allocatable :: RH2PO4TransfBandHeter(:,:)
  real(r8),allocatable :: RDOMEheter(:,:,:,:)
  real(r8),allocatable :: RHOMEheter(:,:,:,:)
  real(r8),allocatable :: RCOMEheter(:,:,:,:)
  real(r8),allocatable :: RH2GX(:,:)
  real(r8),allocatable :: RDMMEheter(:,:,:,:)
  real(r8),allocatable :: RHMMEheter(:,:,:,:)
  real(r8),allocatable :: RCMMEheter(:,:,:,:)
  real(r8),allocatable :: RCCMEheter(:,:,:,:)
  real(r8),allocatable :: RN2FixHeter(:,:)
  real(r8),allocatable :: RXOMEheter(:,:,:,:)
  real(r8),allocatable :: R3OMEheter(:,:,:,:)
  real(r8),allocatable :: RXMMEheter(:,:,:,:)
  real(r8),allocatable :: R3MMEheter(:,:,:,:)
  real(r8),allocatable :: RNH4TransfLitrHeter(:,:)
  real(r8),allocatable :: RNO3TransfLitrHeter(:,:)
  real(r8),allocatable :: RH2PO4TransfLitrHeter(:,:)
  real(r8),allocatable :: RGOMY(:,:)
  real(r8),allocatable :: ROQCD(:,:)
  real(r8),allocatable :: RCO2ProdHeter(:,:)
  real(r8),allocatable :: RAcettProdHeter(:,:)
  real(r8),allocatable :: RCH4ProdHeter(:,:)
  real(r8),allocatable :: RSOxidSoilAutor(:)
  real(r8),allocatable :: RSOxidBandAutor(:)
  real(r8) :: RTotNH3OxidSoilAutor
  real(r8) :: RTotNH3OxidBandAutor
  real(r8),allocatable :: XOMZ(:,:,:,:)
  real(r8),allocatable :: RH1PO4TransfSoilHeter(:,:)
  real(r8),allocatable :: RH1PO4TransfBandHeter(:,:)
  real(r8),allocatable :: RH1PO4TransfLitrHeter(:,:)

  real(r8),allocatable :: RO2UptkAutor(:)
  real(r8),allocatable :: RGN2Fff(:)
  real(r8),allocatable :: RGOMOff(:)
  real(r8),allocatable :: RO2Dmnd4RespAutor(:)
  real(r8),allocatable :: RO2DmndAutor(:)
  real(r8),allocatable :: RO2Uptk4RespAutor(:)
  real(r8),allocatable :: RNO3UptkAutor(:)
  real(r8),allocatable :: RNO2ReduxAutorSoil(:)
  real(r8),allocatable :: RNO2ReduxAutorBand(:)
  real(r8),allocatable :: RGOMDff(:)
  real(r8),allocatable :: RMaintCompAutor(:,:)
  real(r8),allocatable :: RNH4TransfSoilAutor(:)
  real(r8),allocatable :: RNO3TransfSoilAutor(:)
  real(r8),allocatable :: RH2PO4TransfSoilAutor(:)
  real(r8),allocatable :: RNH4TransfBandAutor(:)
  real(r8),allocatable :: RNO3TransfBandAutor(:)
  real(r8),allocatable :: RH2PO4TransfBandAutor(:)
  real(r8),allocatable :: RDOMEautor(:,:,:)
  real(r8),allocatable :: RHOMEautor(:,:,:)
  real(r8),allocatable :: RCOMEautor(:,:,:)
  real(r8),allocatable :: CGOMEautor(:,:)
  real(r8),allocatable :: RH2GXff(:)
  real(r8),allocatable :: RDMMEautor(:,:,:)
  real(r8),allocatable :: RHMMEautor(:,:,:)
  real(r8),allocatable :: RCMMEautor(:,:,:)
  real(r8),allocatable :: RN2FixAutor(:)
  real(r8),allocatable :: RXOMEautor(:,:,:)
  real(r8),allocatable :: R3OMEautor(:,:,:)
  real(r8),allocatable :: RXMMEautor(:,:,:)
  real(r8),allocatable :: R3MMEautor(:,:,:)
  real(r8),allocatable :: RNH4TransfLitrAutor(:)
  real(r8),allocatable :: RNO3TransfLitrAutor(:)
  real(r8),allocatable :: RH2PO4TransfLitrAutor(:)
  real(r8),allocatable :: FNH4XRff(:)
  real(r8),allocatable :: FNO3XRff(:)
  real(r8),allocatable :: FPO4XRff(:)
  real(r8),allocatable :: RGOMYff(:)
  real(r8),allocatable :: CGOSEautor(:,:,:)
  real(r8),allocatable :: FP14XRff(:)
  real(r8),allocatable :: RCO2ProdAutor(:)
  real(r8),allocatable :: RCH4ProdAutor(:)
  real(r8),allocatable :: RH1PO4TransfSoilAutor(:)
  real(r8),allocatable :: RH1PO4TransfBandAutor(:)
  real(r8),allocatable :: RH1PO4TransfLitrAutor(:)
  contains
    procedure, public :: Init => nit_micf_init
    procedure, public :: ZeroOut => nit_micf_zero
    procedure, public :: destroy => nit_micf_destroy
  end type NitroMicFluxType

  type, public :: NitroOMcplxFluxType
    real(r8),allocatable :: RDOSM(:,:,:)
    real(r8),allocatable :: RHOSM(:,:,:)
    real(r8),allocatable :: RCOSM(:,:,:)
    real(r8),allocatable :: RDORM(:,:,:)
    real(r8),allocatable :: RDOHM(:,:)
    real(r8),allocatable :: OMSORP(:,:)
    real(r8),allocatable :: TCGOMEheter(:,:)
    real(r8),allocatable :: ROQCK(:)
    real(r8),allocatable :: XOQCK(:)
    real(r8),allocatable :: XOQMZ(:,:)
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
    real(r8),allocatable :: CDOM(:,:)
  contains
    procedure, public :: Init => nit_omcplxs_init
    procedure, public :: ZeroOut => nit_omcplxs_zero
    procedure, public :: Destroy => nit_omcplxs_destroy
  end type NitroOMcplxStateType

  type, public :: NitroMicDiagType
  real(r8) :: H1P4T
  real(r8) :: H2P4T
  real(r8) :: RH2GZ
  real(r8) :: RNO2ReduxSoilChemo
  real(r8) :: RNO2ReduxBandChemo
  real(r8) :: RN2OProdSoilChemo
  real(r8) :: RN2OProdBandChemo
  real(r8) :: RNO3ProdSoilChemo
  real(r8) :: RNO3ProdBandChemo
  real(r8) :: RNO2ReduxChemo
  real(r8) :: ThetaLitr
  real(r8) :: ThetaZ
  real(r8) :: TORC
  real(r8) :: TOMA
  real(r8) :: TOMN
  real(r8) :: TSensGrowth
  real(r8) :: TSensMaintR
  real(r8) :: VOLWZ
  real(r8) :: WatStressMicb
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

  subroutine nit_micf_init(this,jcplx,NumMicbFunGroups)
  implicit none
  class(NitroMicFluxType) :: this
  integer, intent(in) :: jcplx,NumMicbFunGroups
  integer :: ndbiomcp
  integer :: NumMicrobAutrophCmplx
  integer :: NumMicrbHetetrophCmplx
  ndbiomcp=micpar%ndbiomcp
  NumMicrobAutrophCmplx=micpar%NumMicrobAutrophCmplx
  NumMicrbHetetrophCmplx=micpar%NumMicrbHetetrophCmplx

  allocate(this%RO2UptkHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RO2UptkHeter=spval
  allocate(this%RGN2F(NumMicrbHetetrophCmplx,1:jcplx));this%RGN2F=spval
  allocate(this%RGOMO(NumMicrbHetetrophCmplx,1:jcplx));this%RGOMO=spval
  allocate(this%RO2Dmnd4RespHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RO2Dmnd4RespHeter=spval
  allocate(this%RO2DmndHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RO2DmndHeter=spval
  allocate(this%RO2Uptk4RespHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RO2Uptk4RespHeter=spval
  allocate(this%RNO3ReduxHeterSoil(NumMicrbHetetrophCmplx,1:jcplx));this%RNO3ReduxHeterSoil=spval
  allocate(this%RNO3ReduxHeterBand(NumMicrbHetetrophCmplx,1:jcplx));this%RNO3ReduxHeterBand=spval
  allocate(this%RNO2ReduxHeterSoil(NumMicrbHetetrophCmplx,1:jcplx));this%RNO2ReduxHeterSoil=spval
  allocate(this%RNO2ReduxHeterBand(NumMicrbHetetrophCmplx,1:jcplx));this%RNO2ReduxHeterBand=spval
  allocate(this%RN2OReduxHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RN2OReduxHeter=spval
  allocate(this%RGOMD(NumMicrbHetetrophCmplx,1:jcplx));this%RGOMD=spval
  allocate(this%RMOMC(2,NumMicrbHetetrophCmplx,1:jcplx));this%RMOMC=spval
  allocate(this%RNH4TransfSoilHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RNH4TransfSoilHeter=spval
  allocate(this%RNO3TransfSoilHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RNO3TransfSoilHeter=spval
  allocate(this%RH2PO4TransfSoilHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RH2PO4TransfSoilHeter=spval
  allocate(this%RNH4TransfBandHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RNH4TransfBandHeter=spval
  allocate(this%RNO3TransfBandHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RNO3TransfBandHeter=spval
  allocate(this%RH2PO4TransfBandHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RH2PO4TransfBandHeter=spval
  allocate(this%RDOMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RDOMEheter=spval
  allocate(this%RHOMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RHOMEheter=spval
  allocate(this%RCOMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RCOMEheter=spval
  allocate(this%CGOMEheter(NumPlantChemElms,NumMicrbHetetrophCmplx,1:jcplx));this%CGOMEheter=spval
  allocate(this%RH2GX(NumMicrbHetetrophCmplx,1:jcplx));this%RH2GX=spval
  allocate(this%RDMMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RDMMEheter=spval
  allocate(this%RHMMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RHMMEheter=spval
  allocate(this%RCMMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RCMMEheter=spval
  allocate(this%RN2FixHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RN2FixHeter=spval
  allocate(this%RXOMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RXOMEheter=spval
  allocate(this%R3OMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%R3OMEheter=spval
  allocate(this%RXMMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%RXMMEheter=spval
  allocate(this%R3MMEheter(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%R3MMEheter=spval
  allocate(this%CGOQC(NumMicrbHetetrophCmplx,1:jcplx));this%CGOQC=spval
  allocate(this%CGOAC(NumMicrbHetetrophCmplx,1:jcplx));this%CGOAC=spval
  allocate(this%RNH4TransfLitrHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RNH4TransfLitrHeter=spval
  allocate(this%RNO3TransfLitrHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RNO3TransfLitrHeter=spval
  allocate(this%RH2PO4TransfLitrHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RH2PO4TransfLitrHeter=spval
  allocate(this%FNH4XR(NumMicrbHetetrophCmplx,1:jcplx));this%FNH4XR=spval
  allocate(this%FNO3XR(NumMicrbHetetrophCmplx,1:jcplx));this%FNO3XR=spval
  allocate(this%FPO4XR(NumMicrbHetetrophCmplx,1:jcplx));this%FPO4XR=spval
  allocate(this%RGOMY(NumMicrbHetetrophCmplx,1:jcplx));this%RGOMY=spval
  allocate(this%CGOMES(NumPlantChemElms,2,NumMicrbHetetrophCmplx,1:jcplx));this%CGOMES=spval
  allocate(this%FP14XR(NumMicrbHetetrophCmplx,1:jcplx));this%FP14XR=spval
  allocate(this%RCO2ProdHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RCO2ProdHeter=spval
  allocate(this%RAcettProdHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RAcettProdHeter=spval
  allocate(this%RCH4ProdHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RCH4ProdHeter=spval
  allocate(this%RH1PO4TransfSoilHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RH1PO4TransfSoilHeter=spval
  allocate(this%RH1PO4TransfBandHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RH1PO4TransfBandHeter=spval
  allocate(this%RH1PO4TransfLitrHeter(NumMicrbHetetrophCmplx,1:jcplx));this%RH1PO4TransfLitrHeter=spval

  allocate(this%RSOxidSoilAutor(NumMicrobAutrophCmplx));this%RSOxidSoilAutor=spval
  allocate(this%RSOxidBandAutor(NumMicrobAutrophCmplx));this%RSOxidBandAutor=spval
  allocate(this%XOMZ(1:NumPlantChemElms,3,NumMicrbHetetrophCmplx,1:jcplx));this%XOMZ=spval
  allocate(this%ROQCD(NumMicrbHetetrophCmplx,1:jcplx));this%ROQCD=spval
  allocate(this%RCCMEheter(NumPlantChemElms,ndbiomcp,NumMicrbHetetrophCmplx,1:jcplx));this%RCCMEheter=spval

  allocate(this%RO2UptkAutor(NumMicrobAutrophCmplx));this%RO2UptkAutor=spval
  allocate(this%RGN2Fff(NumMicrobAutrophCmplx));this%RGN2Fff=spval
  allocate(this%RGOMOff(NumMicrobAutrophCmplx));this%RGOMOff=spval
  allocate(this%RO2Dmnd4RespAutor(NumMicrobAutrophCmplx));this%RO2Dmnd4RespAutor=spval
  allocate(this%RO2DmndAutor(NumMicrobAutrophCmplx));this%RO2DmndAutor=spval
  allocate(this%RO2Uptk4RespAutor(NumMicrobAutrophCmplx));this%RO2Uptk4RespAutor=spval
  allocate(this%RNO3UptkAutor(NumMicrobAutrophCmplx));this%RNO3UptkAutor=spval
  allocate(this%RNO2ReduxAutorSoil(NumMicrobAutrophCmplx));this%RNO2ReduxAutorSoil=spval
  allocate(this%RNO2ReduxAutorBand(NumMicrobAutrophCmplx));this%RNO2ReduxAutorBand=spval
  allocate(this%RGOMDff(NumMicrobAutrophCmplx));this%RGOMDff=spval
  allocate(this%RMaintCompAutor(2,NumMicrobAutrophCmplx));this%RMaintCompAutor=spval
  allocate(this%RNH4TransfSoilAutor(NumMicrobAutrophCmplx));this%RNH4TransfSoilAutor=spval
  allocate(this%RNO3TransfSoilAutor(NumMicrobAutrophCmplx));this%RNO3TransfSoilAutor=spval
  allocate(this%RH2PO4TransfSoilAutor(NumMicrobAutrophCmplx));this%RH2PO4TransfSoilAutor=spval
  allocate(this%RNH4TransfBandAutor(NumMicrobAutrophCmplx));this%RNH4TransfBandAutor=spval
  allocate(this%RNO3TransfBandAutor(NumMicrobAutrophCmplx));this%RNO3TransfBandAutor=spval
  allocate(this%RH2PO4TransfBandAutor(NumMicrobAutrophCmplx));this%RH2PO4TransfBandAutor=spval
  allocate(this%RDOMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RDOMEautor=spval
  allocate(this%RHOMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RHOMEautor=spval
  allocate(this%RCOMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RCOMEautor=spval
  allocate(this%CGOMEautor(idom_beg:idom_end,NumMicrobAutrophCmplx));this%CGOMEautor=spval
  allocate(this%RH2GXff(NumMicrobAutrophCmplx));this%RH2GXff=spval
  allocate(this%RDMMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RDMMEautor=spval
  allocate(this%RHMMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RHMMEautor=spval
  allocate(this%RCMMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RCMMEautor=spval
  allocate(this%RN2FixAutor(NumMicrobAutrophCmplx));this%RN2FixAutor=spval
  allocate(this%RXOMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RXOMEautor=spval
  allocate(this%R3OMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%R3OMEautor=spval
  allocate(this%RXMMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RXMMEautor=spval
  allocate(this%R3MMEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%R3MMEautor=spval
  allocate(this%RNH4TransfLitrAutor(NumMicrobAutrophCmplx));this%RNH4TransfLitrAutor=spval
  allocate(this%RNO3TransfLitrAutor(NumMicrobAutrophCmplx));this%RNO3TransfLitrAutor=spval
  allocate(this%RH2PO4TransfLitrAutor(NumMicrobAutrophCmplx));this%RH2PO4TransfLitrAutor=spval
  allocate(this%FNH4XRff(NumMicrobAutrophCmplx));this%FNH4XRff=spval
  allocate(this%FNO3XRff(NumMicrobAutrophCmplx));this%FNO3XRff=spval
  allocate(this%FPO4XRff(NumMicrobAutrophCmplx));this%FPO4XRff=spval
  allocate(this%RGOMYff(NumMicrobAutrophCmplx));this%RGOMYff=spval
  allocate(this%CGOSEautor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%CGOSEautor=spval
  allocate(this%FP14XRff(NumMicrobAutrophCmplx));this%FP14XRff=spval
  allocate(this%RCO2ProdAutor(NumMicrobAutrophCmplx));this%RCO2ProdAutor=spval
  allocate(this%RCH4ProdAutor(NumMicrobAutrophCmplx));this%RCH4ProdAutor=spval
  allocate(this%RH1PO4TransfSoilAutor(NumMicrobAutrophCmplx));this%RH1PO4TransfSoilAutor=spval
  allocate(this%RH1PO4TransfBandAutor(NumMicrobAutrophCmplx));this%RH1PO4TransfBandAutor=spval
  allocate(this%RH1PO4TransfLitrAutor(NumMicrobAutrophCmplx));this%RH1PO4TransfLitrAutor=spval

  call this%ZeroOut()
  end subroutine nit_micf_init
!------------------------------------------------------------------------------------------

  subroutine nit_mics_init(this, jcplx,NumMicbFunGroups)

  implicit none
  class(NitroMicStateType) :: this
  integer, intent(in) :: jcplx,NumMicbFunGroups
  integer :: NumMicrobAutrophCmplx,NumMicrbHetetrophCmplx
  NumMicrobAutrophCmplx=micpar%NumMicrobAutrophCmplx
  NumMicrbHetetrophCmplx=micpar%NumMicrbHetetrophCmplx

  allocate(this%CNOMActHeter(NumMicrbHetetrophCmplx,1:jcplx));this%CNOMActHeter=spval
  allocate(this%CPOMActHeter(NumMicrbHetetrophCmplx,1:jcplx));this%CPOMActHeter=spval
  allocate(this%OMActHeter(NumMicrbHetetrophCmplx,1:jcplx));this%OMActHeter=spval
  allocate(this%FracOMActHeter(NumMicrbHetetrophCmplx,1:jcplx));this%FracOMActHeter=spval
  allocate(this%FOMN(NumMicrbHetetrophCmplx,1:jcplx));this%FOMN=spval
  allocate(this%FOMK(NumMicrbHetetrophCmplx,1:jcplx));this%FOMK=spval
  allocate(this%OMC2(NumMicrbHetetrophCmplx,1:jcplx));this%OMC2=spval
  allocate(this%TFNG(NumMicrbHetetrophCmplx,1:jcplx));this%TFNG=spval
  allocate(this%TFNR(NumMicrbHetetrophCmplx,1:jcplx));this%TFNR=spval
  allocate(this%OMN2(NumMicrbHetetrophCmplx,1:jcplx));this%OMN2=spval
  allocate(this%FOM2(NumMicrbHetetrophCmplx,1:jcplx));this%FOM2=spval
  allocate(this%WFN(NumMicrbHetetrophCmplx,1:jcplx));this%WFN=spval
  allocate(this%FCN(NumMicrbHetetrophCmplx,1:jcplx));this%FCN=spval
  allocate(this%FCP(NumMicrbHetetrophCmplx,1:jcplx));this%FCP=spval
  allocate(this%FCNP(NumMicrbHetetrophCmplx,1:jcplx));this%FCNP=spval

  allocate(this%CNOMActAutor(NumMicrobAutrophCmplx));this%CNOMActAutor=spval
  allocate(this%CPOMActAutor(NumMicrobAutrophCmplx));this%CPOMActAutor=spval
  allocate(this%OMActAutor(NumMicrobAutrophCmplx));this%OMActAutor=spval
  allocate(this%FracOMActAutor(NumMicrobAutrophCmplx));this%FracOMActAutor=spval
  allocate(this%FOMNff(NumMicrobAutrophCmplx));this%FOMNff=spval
  allocate(this%FOMKff(NumMicrobAutrophCmplx));this%FOMKff=spval
  allocate(this%OMC2ff(NumMicrobAutrophCmplx));this%OMC2ff=spval
  allocate(this%TFNGff(NumMicrobAutrophCmplx));this%TFNGff=spval
  allocate(this%TFNRff(NumMicrobAutrophCmplx));this%TFNRff=spval
  allocate(this%OMN2ff(NumMicrobAutrophCmplx));this%OMN2ff=spval
  allocate(this%FOM2ff(NumMicrobAutrophCmplx));this%FOM2ff=spval
  allocate(this%fLimO2Autor(NumMicrobAutrophCmplx));this%fLimO2Autor=spval
  allocate(this%FCNff(NumMicrobAutrophCmplx));this%FCNff=spval
  allocate(this%FCPff(NumMicrobAutrophCmplx));this%FCPff=spval
  allocate(this%FCNPff(NumMicrobAutrophCmplx));this%FCNPff=spval
  end subroutine nit_mics_init
!------------------------------------------------------------------------------------------

  subroutine nit_micf_zero(this)
  implicit none
  class(NitroMicFLuxType) :: this

  this%RO2UptkHeter = 0._r8
  this%RGN2F = 0._r8
  this%RGOMO = 0._r8
  this%RO2Dmnd4RespHeter = 0._r8
  this%RO2DmndHeter = 0._r8
  this%RO2Uptk4RespHeter = 0._r8
  this%RNO3ReduxHeterSoil = 0._r8
  this%RNO3ReduxHeterBand = 0._r8
  this%RNO2ReduxHeterSoil = 0._r8
  this%RNO2ReduxHeterBand = 0._r8
  this%RN2OReduxHeter = 0._r8
  this%RGOMD = 0._r8
  this%RMOMC = 0._r8
  this%RNH4TransfSoilHeter = 0._r8
  this%RNO3TransfSoilHeter = 0._r8
  this%RH2PO4TransfSoilHeter = 0._r8
  this%RNH4TransfBandHeter = 0._r8
  this%RNO3TransfBandHeter = 0._r8
  this%RH2PO4TransfBandHeter = 0._r8
  this%RDOMEheter = 0._r8
  this%RHOMEheter = 0._r8
  this%RCOMEheter = 0._r8
  this%CGOMEheter = 0._r8
  this%RH2GX = 0._r8
  this%RHMMEheter = 0._r8
  this%RCMMEheter = 0._r8
  this%RDMMEheter = 0._r8
  this%RCCMEheter = 0._r8
  this%RN2FixHeter = 0._r8

  this%RXOMEheter = 0._r8
  this%R3OMEheter = 0._r8
  this%RXMMEheter = 0._r8
  this%R3MMEheter = 0._r8
  this%CGOQC = 0._r8
  this%CGOAC = 0._r8
  this%RNH4TransfLitrHeter = 0._r8
  this%RNO3TransfLitrHeter = 0._r8
  this%RH2PO4TransfLitrHeter = 0._r8
  this%FNH4XR = 0._r8
  this%FNO3XR = 0._r8
  this%FPO4XR = 0._r8
  this%RGOMY = 0._r8
  this%ROQCD = 0._r8
  this%CGOMES = 0._r8
  this%FP14XR = 0._r8
  this%RCO2ProdHeter = 0._r8
  this%RAcettProdHeter = 0._r8
  this%RCH4ProdHeter = 0._r8
  this%RSOxidSoilAutor = 0._r8
  this%RSOxidBandAutor = 0._r8
  this%XOMZ = 0._r8
  this%RH1PO4TransfSoilHeter = 0._r8
  this%RH1PO4TransfBandHeter = 0._r8
  this%RH1PO4TransfLitrHeter = 0._r8

  this%RO2UptkAutor = 0._r8
  this%RGN2Fff = 0._r8
  this%RGOMOff = 0._r8
  this%RO2Dmnd4RespAutor = 0._r8
  this%RO2DmndAutor = 0._r8
  this%RO2Uptk4RespAutor = 0._r8
  this%RNO3UptkAutor = 0._r8
  this%RNO2ReduxAutorSoil = 0._r8
  this%RNO2ReduxAutorBand = 0._r8
  this%RGOMDff = 0._r8
  this%RMaintCompAutor = 0._r8
  this%RNH4TransfSoilAutor = 0._r8
  this%RNO3TransfSoilAutor = 0._r8
  this%RH2PO4TransfSoilAutor = 0._r8
  this%RNH4TransfBandAutor = 0._r8
  this%RNO3TransfBandAutor = 0._r8
  this%RH2PO4TransfBandAutor = 0._r8
  this%RDOMEautor = 0._r8
  this%RHOMEautor = 0._r8
  this%RCOMEautor = 0._r8
  this%CGOMEautor = 0._r8
  this%RH2GXff = 0._r8
  this%RDMMEautor = 0._r8
  this%RHMMEautor = 0._r8
  this%RCMMEautor = 0._r8
  this%RN2FixAutor = 0._r8

  this%RXOMEautor = 0._r8
  this%R3OMEautor = 0._r8
  this%RXMMEautor = 0._r8
  this%R3MMEautor = 0._r8
  this%RNH4TransfLitrAutor = 0._r8
  this%RNO3TransfLitrAutor = 0._r8
  this%RH2PO4TransfLitrAutor = 0._r8
  this%FNH4XRff = 0._r8
  this%FNO3XRff = 0._r8
  this%FPO4XRff = 0._r8
  this%RGOMYff = 0._r8
  this%CGOSEautor = 0._r8
  this%FP14XRff = 0._r8
  this%RCO2ProdAutor = 0._r8
  this%RCH4ProdAutor = 0._r8
  this%RH1PO4TransfSoilAutor = 0._r8
  this%RH1PO4TransfBandAutor = 0._r8
  this%RH1PO4TransfLitrAutor = 0._r8

  end subroutine nit_micf_zero


!------------------------------------------------------------------------------------------

  subroutine nit_micf_destroy(this)

  implicit none
  class(NitroMicFLuxType) :: this

  call destroy(this%RO2UptkHeter)
  call destroy(this%RGN2F)
  call destroy(this%RGOMO)
  call destroy(this%RO2Dmnd4RespHeter)
  call destroy(this%RO2DmndHeter)
  call destroy(this%RO2Uptk4RespHeter)
  call destroy(this%RNO3ReduxHeterSoil)
  call destroy(this%RNO3ReduxHeterBand)
  call destroy(this%RNO2ReduxHeterSoil)
  call destroy(this%RNO2ReduxHeterBand)
  call destroy(this%RN2OReduxHeter)
  call destroy(this%RGOMD)
  call destroy(this%RMOMC)
  call destroy(this%RNH4TransfSoilHeter)
  call destroy(this%RNO3TransfSoilHeter)
  call destroy(this%RH2PO4TransfSoilHeter)
  call destroy(this%RNH4TransfBandHeter)
  call destroy(this%RNO3TransfBandHeter)
  call destroy(this%RH2PO4TransfBandHeter)
  call destroy(this%RDOMEheter)
  call destroy(this%RHOMEheter)
  call destroy(this%RCOMEheter)
  call destroy(this%CGOMEheter)
  call destroy(this%RH2GX)
  call destroy(this%RHMMEheter)
  call destroy(this%RCMMEheter)
  call destroy(this%RDMMEheter)
  call destroy(this%RCCMEheter)
  call destroy(this%RN2FixHeter)
  call destroy(this%RXOMEheter)
  call destroy(this%R3OMEheter)
  call destroy(this%RXMMEheter)
  call destroy(this%R3MMEheter)
  call destroy(this%CGOQC)
  call destroy(this%CGOAC)
  call destroy(this%RNH4TransfLitrHeter)
  call destroy(this%RNO3TransfLitrHeter)
  call destroy(this%RH2PO4TransfLitrHeter)
  call destroy(this%FNH4XR)
  call destroy(this%FNO3XR)
  call destroy(this%FPO4XR)
  call destroy(this%RGOMY)
  call destroy(this%ROQCD)
  call destroy(this%CGOMES)
  call destroy(this%FP14XR)
  call destroy(this%RCO2ProdHeter)
  call destroy(this%RAcettProdHeter)
  call destroy(this%RCH4ProdHeter)
  call destroy(this%RSOxidSoilAutor)
  call destroy(this%RSOxidBandAutor)
  call destroy(this%XOMZ)
  call destroy(this%RH1PO4TransfSoilHeter)
  call destroy(this%RH1PO4TransfBandHeter)
  call destroy(this%RH1PO4TransfLitrHeter)

  call destroy(this%RO2UptkAutor)
  call destroy(this%RGN2Fff)
  call destroy(this%RGOMOff)
  call destroy(this%RO2Dmnd4RespAutor)
  call destroy(this%RO2DmndAutor)
  call destroy(this%RO2Uptk4RespAutor)
  call destroy(this%RNO3UptkAutor)
  call destroy(this%RNO2ReduxAutorSoil)
  call destroy(this%RNO2ReduxAutorBand)
  call destroy(this%RGOMDff)
  call destroy(this%RMaintCompAutor)
  call destroy(this%RNH4TransfSoilAutor)
  call destroy(this%RNO3TransfSoilAutor)
  call destroy(this%RH2PO4TransfSoilAutor)
  call destroy(this%RNH4TransfBandAutor)
  call destroy(this%RNO3TransfBandAutor)
  call destroy(this%RH2PO4TransfBandAutor)
  call destroy(this%RDOMEautor)
  call destroy(this%RHOMEautor)
  call destroy(this%RCOMEautor)
  call destroy(this%CGOMEautor)
  call destroy(this%RH2GXff)
  call destroy(this%RDMMEautor)
  call destroy(this%RHMMEautor)
  call destroy(this%RCMMEautor)
  call destroy(this%RN2FixAutor)
  call destroy(this%RXOMEautor)
  call destroy(this%R3OMEautor)
  call destroy(this%RXMMEautor)
  call destroy(this%R3MMEautor)
  call destroy(this%RNH4TransfLitrAutor)
  call destroy(this%RNO3TransfLitrAutor)
  call destroy(this%RH2PO4TransfLitrAutor)
  call destroy(this%FNH4XRff)
  call destroy(this%FNO3XRff)
  call destroy(this%FPO4XRff)
  call destroy(this%RGOMYff)
  call destroy(this%CGOSEautor)
  call destroy(this%FP14XRff)
  call destroy(this%RCO2ProdAutor)
  call destroy(this%RCH4ProdAutor)
  call destroy(this%RH1PO4TransfSoilAutor)
  call destroy(this%RH1PO4TransfBandAutor)
  call destroy(this%RH1PO4TransfLitrAutor)

  end subroutine nit_micf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_mics_destroy(this)
  implicit none
  class(NitroMicStateType) :: this

  call destroy(this%CNOMActHeter)
  call destroy(this%CPOMActHeter)
  call destroy(this%OMActHeter)
  call destroy(this%FracOMActHeter)
  call destroy(this%FOMN)
  call destroy(this%FOMK)
  call destroy(this%OMC2)
  call destroy(this%TFNG)
  call destroy(this%TFNR)
  call destroy(this%OMN2)
  call destroy(this%FOM2)
  call destroy(this%WFN)
  call destroy(this%FCN)
  call destroy(this%FCP)
  call destroy(this%FCNP)

  call destroy(this%CNOMActAutor)
  call destroy(this%CPOMActAutor)
  call destroy(this%OMActAutor)
  call destroy(this%FracOMActAutor)
  call destroy(this%FOMNff)
  call destroy(this%FOMKff)
  call destroy(this%OMC2ff)
  call destroy(this%TFNGff)
  call destroy(this%TFNRff)
  call destroy(this%OMN2ff)
  call destroy(this%FOM2ff)
  call destroy(this%fLimO2Autor)
  call destroy(this%FCNff)
  call destroy(this%FCPff)
  call destroy(this%FCNPff)
  end subroutine nit_mics_destroy
!------------------------------------------------------------------------------------------
  subroutine nit_omcplxf_init(this)
  implicit none
  class(NitroOMcplxFluxType) :: this
  integer :: nkinets
  integer :: ncplx
  integer :: ndbiomcp

  nkinets=micpar%jsken
  ncplx=micpar%jcplx
  ndbiomcp=micpar%ndbiomcp
  allocate(this%RDOSM(NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RHOSM(1:NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RCOSM(1:NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RDORM(1:NumPlantChemElms,ndbiomcp,1:ncplx))
  allocate(this%RDOHM(idom_beg:idom_end,1:ncplx))
  allocate(this%OMSORP(idom_beg:idom_end,1:ncplx))
  allocate(this%TCGOMEheter(idom_beg:idom_end,1:ncplx+1))
  allocate(this%ROQCK(1:ncplx))
  allocate(this%XOQCK(1:ncplx))
  allocate(this%XOQMZ(idom_beg:idom_end,1:ncplx))

  call this%ZeroOut()
  end subroutine nit_omcplxf_init
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_destroy(this)
  implicit none
  class(NitroOMcplxFluxType) :: this

  call destroy(this%RDOSM)
  call destroy(this%RHOSM)
  call destroy(this%RCOSM)
  call destroy(this%RDORM)
  call destroy(this%RDOHM)
  call destroy(this%OMSORP)
  call destroy(this%ROQCK)
  call destroy(this%XOQCK)
  call destroy(this%XOQMZ)

  end subroutine nit_omcplxf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_zero(this)
  implicit none
  class(NitroOMcplxFluxType) :: this

  this%RDOSM=0._r8
  this%RHOSM=0._r8
  this%RCOSM=0._r8
  this%RDORM=0._r8
  this%RDOHM=0._r8
  this%OMSORP=0._r8
  this%ROQCK=0._r8
  this%XOQCK=0._r8
  this%XOQMZ=0._r8
  end subroutine nit_omcplxf_zero
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxs_init(this)
  implicit none
  class(NitroOMcplxStateType) :: this
  integer :: ncplx

  ncplx=micpar%jcplx
  allocate(this%OSRH(1:ncplx));this%OSRH=spval
  allocate(this%TOMK(1:ncplx+1));this%TOMK=spval
  allocate(this%TONK(1:ncplx+1));this%TONK=spval
  allocate(this%TOPK(1:ncplx+1));this%TOPK=spval
  allocate(this%FOCA(1:ncplx));this%FOCA=spval
  allocate(this%FOAA(1:ncplx));this%FOAA=spval
  allocate(this%CNQ(1:ncplx));this%CNQ=spval
  allocate(this%CPQ(1:ncplx));this%CPQ=spval
  allocate(this%CNH(1:ncplx));this%CNH=spval
  allocate(this%CPH(1:ncplx));this%CPH=spval
  allocate(this%ORCT(1:ncplx));this%ORCT=spval
  allocate(this%OSCT(1:ncplx));this%OSCT=spval
  allocate(this%OSAT(1:ncplx));this%OSAT=spval
  allocate(this%TONX(1:ncplx+1));this%TONX=spval
  allocate(this%TOPX(1:ncplx+1));this%TOPX=spval
  allocate(this%CDOM(idom_beg:idom_end,1:ncplx));this%CDOM=spval
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

  call destroy(this%OSRH)
  call destroy(this%TOMK)
  call destroy(this%TONK)
  call destroy(this%TOPK)
  call destroy(this%FOCA)
  call destroy(this%FOAA)
  call destroy(this%CNQ)
  call destroy(this%CPQ)
  call destroy(this%CNH)
  call destroy(this%CPH)
  call destroy(this%ORCT)
  call destroy(this%OSCT)
  call destroy(this%OSAT)
  call destroy(this%TONX)
  call destroy(this%TOPX)
  call destroy(this%CDOM)

  end subroutine nit_omcplxs_destroy
end module NitroDiagTypes
