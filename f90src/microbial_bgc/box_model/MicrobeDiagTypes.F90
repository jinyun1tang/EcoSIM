module MicrobeDiagTypes
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
type, public :: Cumlate_Flux_Diag_type
! fraction diagnostics, used for substrate competition/uptake
    real(r8) :: TFNH4B
    real(r8) :: TFNO3B
    real(r8) :: TFNO2B
    real(r8) :: TFP14B
    real(r8) :: TFPO4B
    real(r8) :: tCH4ProdAceto    !acetoclastic CH4 production 
    real(r8) :: tCH4ProdH2       !Hydrogenotrophic methane production
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
    real(r8) :: tCResp4H2Prod
    real(r8) :: tRNH4MicrbTransfSoil
    real(r8) :: tRH2PO4MicrbTransfSoil
    real(r8) :: tRNO3MicrbTransfSoil
    real(r8) :: tRH1PO4MicrbTransfSoil
    real(r8) :: tRNH4MicrbTransfBand
    real(r8) :: tRNO3MicrbTransfBand
    real(r8) :: tRH2PO4MicrbTransfBand
    real(r8) :: tRH1PO4MicrbTransfBand
    real(r8) :: tRCO2MicrbProd
    real(r8) :: tRCH4MicrbProd
    real(r8) :: tRNOxMicrbRedux
    real(r8) :: tRCO2GrothAutor       !CO2 taken up for autotrophic biomass growth
    real(r8) :: TProdH2               !H2 production by fermenters
    real(r8) :: tRO2MicrbUptk         !O2 uptake by microbes
    real(r8) :: TReduxNO3Soil         !NO3 reduction in non-band soil by denitrifiers
    real(r8) :: TReduxNO3Band         !NO3 reduction in banded soil by denitrifiers
    real(r8) :: TReduxNO2Soil         !NO2 reduction in non-band soil, by ammonia oxidizers, and denitrifiers
    real(r8) :: TReduxNO2Band         !NO2 reduction in banded soil, by ammonia oxidizers, and denitrifiers
    real(r8) :: TDeniReduxNO2Band
    real(r8) :: TDeniReduxNO2Soil
    real(r8) :: TNitReduxNO2Band
    real(r8) :: TNitReduxNO2Soil

    real(r8) :: TReduxN2O             !N2O reduction by detnitrifiers
    real(r8) :: TFixN2                !N2 fixation by aerobic and anaerobic N2 fixers
    real(r8) :: tCH4OxiAero           !aerobic CH4 oxidation
    real(r8) :: tRNH3Oxi              !NH3 oxidation by nitrifiers, soil+band

    real(r8) :: RNO2ReduxSoilChemo
    real(r8) :: RNO2ReduxBandChemo
    real(r8) :: RN2OProdSoilChemo
    real(r8) :: RN2OProdBandChemo
    real(r8) :: RNO3ProdSoilChemo
    real(r8) :: RNO3ProdBandChemo
    real(r8) :: RNO2ReduxChemo

  contains
    procedure, public :: ZeroOut => nit_aqmf_diag
  end type Cumlate_Flux_Diag_type

  type, public :: Microbe_State_type
  real(r8),allocatable :: rCNBiomeActHeter(:,:,:)  !Nutrient to carbon ratio for active heterotrophs
  real(r8),allocatable :: OMActHeter(:,:)

  real(r8),allocatable :: OMC2(:,:)
  real(r8),allocatable :: GrowthEnvScalHeter(:,:)
  real(r8),allocatable :: TempMaintRHeter(:,:)
  real(r8),allocatable :: OMN2(:,:)
  real(r8),allocatable :: FOM2(:,:)
  real(r8),allocatable :: OxyLimterHeter(:,:)

  real(r8),allocatable :: FCN(:,:)
  real(r8),allocatable :: FCP(:,:)
  real(r8),allocatable :: FBiomStoiScalarHeter(:,:)
  real(r8),allocatable :: FracOMActHeter(:,:)
  real(r8),allocatable :: FracNO2ReduxHeter(:,:)
  real(r8),allocatable :: FracHeterBiomOfActK(:,:)

  real(r8),allocatable :: rCNBiomeActAutor(:,:)   !nutrient to carbon ratio for active autotrophs
  real(r8),allocatable :: OMActAutor(:)
  real(r8),allocatable :: FracOMActAutor(:)
  real(r8),allocatable :: FracNO2ReduxAutor(:)
  real(r8),allocatable :: FracAutorBiomOfActK(:)
  real(r8),allocatable :: OMC2Autor(:)
  real(r8),allocatable :: GrowthEnvScalAutor(:)
  real(r8),allocatable :: TSensMaintRAutor(:)
  real(r8),allocatable :: OMN2Autor(:)
  real(r8),allocatable :: FOM2Autor(:)
  real(r8),allocatable :: fLimO2Autor(:)
  real(r8),allocatable :: FCNAutor(:)
  real(r8),allocatable :: FCPAutor(:)
  real(r8),allocatable :: FBiomStoiScalarAutor(:)
  contains
   procedure, public :: Init => nit_mics_init
   procedure, public :: Destroy => nit_mics_destroy
  end type Microbe_State_type

  type, public :: Microbe_Flux_type
  !fluxes
  real(r8) :: RTotNH3OxidSoilAutor
  real(r8) :: RTotNH3OxidBandAutor  
! allocatable flux ratios
  real(r8),allocatable :: AttenfNH4Heter(:,:)
  real(r8),allocatable :: AttenfNO3Heter(:,:)
  real(r8),allocatable :: AttenfH2PO4Heter(:,:)
  real(r8),allocatable :: AttenfH1PO4Heter(:,:)
! allocatable fluxes
  real(r8),allocatable :: DOMuptk4GrothHeter(:,:,:)
  real(r8),allocatable :: RAnabolDOCUptkHeter(:,:)
  real(r8),allocatable :: RAnabolAcetUptkHeter(:,:)
  real(r8),allocatable :: NonstX2stBiomHeter(:,:,:,:)
  real(r8),allocatable :: RO2UptkHeter(:,:)
  real(r8),allocatable :: Resp4NFixHeter(:,:)
  real(r8),allocatable :: RespGrossHeter(:,:)
  real(r8),allocatable :: RO2Dmnd4RespHeter(:,:)
  real(r8),allocatable :: RO2DmndHeter(:,:)
  real(r8),allocatable :: RO2Uptk4RespHeter(:,:)
  real(r8),allocatable :: RNO3ReduxHeterSoil(:,:)
  real(r8),allocatable :: RNO3ReduxHeterBand(:,:)
  real(r8),allocatable :: RNO2ReduxHeterSoil(:,:)
  real(r8),allocatable :: RNO2ReduxHeterBand(:,:)
  real(r8),allocatable :: RN2OReduxHeter(:,:)
  real(r8),allocatable :: RNOxReduxRespDenitLim(:,:)
  real(r8),allocatable :: RMaintDmndHeter(:,:,:)
  real(r8),allocatable :: RNH4TransfSoilHeter(:,:)
  real(r8),allocatable :: RNO3TransfSoilHeter(:,:)
  real(r8),allocatable :: RH2PO4TransfSoilHeter(:,:)
  real(r8),allocatable :: RNH4TransfBandHeter(:,:)
  real(r8),allocatable :: RNO3TransfBandHeter(:,:)
  real(r8),allocatable :: RH2PO4TransfBandHeter(:,:)
  real(r8),allocatable :: RkillLitfalOMHeter(:,:,:,:)
  real(r8),allocatable :: RkillLitrfal2HumOMHeter(:,:,:,:)
  real(r8),allocatable :: RkillLitrfal2ResduOMHeter(:,:,:,:)
  real(r8),allocatable :: RH2ProdHeter(:,:)
  real(r8),allocatable :: RMaintDefcitLitrfalOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefcitLitrfal2HumOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefcitLitrfal2ResduOMHeter(:,:,:,:)
  real(r8),allocatable :: RCCMEheter(:,:,:,:)
  real(r8),allocatable :: RN2FixHeter(:,:)
  real(r8),allocatable :: RKillOMHeter(:,:,:,:)
  real(r8),allocatable :: RkillRecycOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefcitKillOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefcitRecycOMHeter(:,:,:,:)
  real(r8),allocatable :: RNH4TransfLitrHeter(:,:)
  real(r8),allocatable :: RNO3TransfLitrHeter(:,:)
  real(r8),allocatable :: RH2PO4TransfLitrHeter(:,:)
  real(r8),allocatable :: RNOxReduxRespDenitUlm(:,:)
  real(r8),allocatable :: ROQC4HeterMicrobAct(:,:)
  real(r8),allocatable :: RCO2ProdHeter(:,:)
  real(r8),allocatable :: RAcettProdHeter(:,:)
  real(r8),allocatable :: RCH4ProdHeter(:,:)
  real(r8),allocatable :: RSOxidSoilAutor(:)
  real(r8),allocatable :: RSOxidBandAutor(:)

  real(r8),allocatable :: XferBiomeHeterK(:,:,:,:)
  real(r8),allocatable :: RH1PO4TransfSoilHeter(:,:)
  real(r8),allocatable :: RH1PO4TransfBandHeter(:,:)
  real(r8),allocatable :: RH1PO4TransfLitrHeter(:,:)

  real(r8),allocatable :: RO2UptkAutor(:)
  real(r8),allocatable :: Resp4NFixAutor(:)
  real(r8),allocatable :: RespGrossAutor(:)
  real(r8),allocatable :: RO2Dmnd4RespAutor(:)
  real(r8),allocatable :: RO2DmndAutor(:)
  real(r8),allocatable :: RO2Uptk4RespAutor(:)
  real(r8),allocatable :: RNO3UptkAutor(:)
  real(r8),allocatable :: RNO2ReduxAutorSoil(:)
  real(r8),allocatable :: RNO2ReduxAutorBand(:)
  real(r8),allocatable :: RNOxReduxRespAutorLim(:)
  real(r8),allocatable :: RMaintDmndAutor(:,:)
  real(r8),allocatable :: RNH4TransfSoilAutor(:)
  real(r8),allocatable :: RNO3TransfSoilAutor(:)
  real(r8),allocatable :: RH2PO4TransfSoilAutor(:)
  real(r8),allocatable :: RNH4TransfBandAutor(:)
  real(r8),allocatable :: RNO3TransfBandAutor(:)
  real(r8),allocatable :: RH2PO4TransfBandAutor(:)
  real(r8),allocatable :: RkillLitfalOMAutor(:,:,:)
  real(r8),allocatable :: RkillLitrfal2HumOMAutor(:,:,:)
  real(r8),allocatable :: RkillLitrfal2ResduOMAutor(:,:,:)
  real(r8),allocatable :: DOMuptk4GrothAutor(:,:)
  real(r8),allocatable :: RMaintDefcitLitrfalOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefcitLitrfal2HumOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefcitLitrfal2ResduOMAutor(:,:,:)
  real(r8),allocatable :: RN2FixAutor(:)
  real(r8),allocatable :: RKillOMAutor(:,:,:)
  real(r8),allocatable :: RkillRecycOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefcitKillOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefcitRecycOMAutor(:,:,:)
  real(r8),allocatable :: RNH4TransfLitrAutor(:)
  real(r8),allocatable :: RNO3TransfLitrAutor(:)
  real(r8),allocatable :: RH2PO4TransfLitrAutor(:)
  real(r8),allocatable :: AttenfNH4Autor(:)
  real(r8),allocatable :: AttenfNO3Autor(:)
  real(r8),allocatable :: AttenfH2PO4Autor(:)
  real(r8),allocatable :: NonstX2stBiomAutor(:,:,:)
  real(r8),allocatable :: AttenfH1PO4Autor(:)
  real(r8),allocatable :: RCO2ProdAutor(:)
  real(r8),allocatable :: RCH4ProdAutor(:)
  real(r8),allocatable :: RH1PO4TransfSoilAutor(:)
  real(r8),allocatable :: RH1PO4TransfBandAutor(:)
  real(r8),allocatable :: RH1PO4TransfLitrAutor(:)

  contains
    procedure, public :: Init => nit_micf_init
    procedure, public :: ZeroOut => nit_micf_zero
    procedure, public :: destroy => nit_micf_destroy
  end type Microbe_Flux_type

  type, public :: OMCplx_Flux_type
    real(r8),allocatable :: RHydlysSolidOM(:,:,:)
    real(r8),allocatable :: RHumifySolidOM(:,:,:)
    real(r8),allocatable :: RDcmpProdDOM(:,:,:)
    real(r8),allocatable :: RHydlysBioResduOM(:,:,:)
    real(r8),allocatable :: RHydlysSorptOM(:,:)
    real(r8),allocatable :: RDOMSorp(:,:)                !rate of DOM adsorption
    real(r8),allocatable :: TDOMUptkHeter(:,:)
    real(r8),allocatable :: ROQC4HeterMicActCmpK(:)
    real(r8),allocatable :: XferRespHeterK(:)
    real(r8),allocatable :: XferDOMK(:,:)
  contains
    procedure, public :: Init => nit_omcplxf_init
    procedure, public :: ZeroOut => nit_omcplxf_zero
    procedure, public :: Destroy => nit_omcplxf_destroy
  end type OMCplx_Flux_type

  type, public :: OMCplx_State_type
    real(r8),allocatable :: TOMEAutoK(:)
    real(r8),allocatable :: BulkSOMC(:)
    real(r8),allocatable :: TOMK(:)
    real(r8),allocatable :: TONK(:)
    real(r8),allocatable :: TOPK(:)
    real(r8),allocatable :: FOCA(:)
    real(r8),allocatable :: FOAA(:)
    real(r8),allocatable :: rCNDOM(:)
    real(r8),allocatable :: rCPDOM(:)
    real(r8),allocatable :: rCNSorbOM(:)
    real(r8),allocatable :: rCPSorbOM(:)
    real(r8),allocatable :: OMBioResduK(:)
    real(r8),allocatable :: SolidOMK(:)
    real(r8),allocatable :: SolidOMActK(:)
    real(r8),allocatable :: tMaxNActMicrbK(:)
    real(r8),allocatable :: tMaxPActMicrbK(:)
    real(r8),allocatable :: CDOM(:,:)            !dom concentration
  contains
    procedure, public :: Init => nit_omcplxs_init
    procedure, public :: ZeroOut => nit_omcplxs_zero
    procedure, public :: Destroy => nit_omcplxs_destroy
  end type OMCplx_State_type

  type, public :: Microbe_Diag_type
  real(r8) :: H1P4T
  real(r8) :: H2P4T
  real(r8) :: RH2UptkAutor
  real(r8) :: ThetaLitr
  real(r8) :: ThetaZ
  real(r8) :: TOMBioResdu
  real(r8) :: TotActMicrobiom
  real(r8) :: TotBiomNO2Consumers
  real(r8) :: TSensGrowth
  real(r8) :: TSensMaintR
  real(r8) :: VOLWZ
  real(r8) :: WatStressMicb
  real(r8) :: XCO2
  real(r8) :: ZNH4T
  real(r8) :: ZNO3T
  real(r8) :: ZNO2T
  end type Microbe_Diag_type

  contains
!------------------------------------------------------------------------------------------

  subroutine nit_aqmf_diag(this)
  implicit none
  class(Cumlate_Flux_Diag_type) :: this

  this%TFNH4B = 0._r8
  this%TFNO3B = 0._r8
  this%TFNO2B = 0._r8
  this%TFP14B = 0._r8
  this%TFPO4B = 0._r8
  this%tCH4ProdAceto = 0._r8
  this%tCH4ProdH2 = 0._r8
  this%TFOQC = 0._r8
  this%TFOQA = 0._r8
  this%TFOXYX = 0._r8
  this%TFNH4X = 0._r8
  this%TFNO3X = 0._r8
  this%TFNO2X = 0._r8
  this%TFN2OX = 0._r8
  this%TFP14X = 0._r8
  this%TFPO4X = 0._r8

  this%tCResp4H2Prod = 0._r8
  this%tRNH4MicrbTransfSoil = 0._r8
  this%tRH2PO4MicrbTransfSoil = 0._r8
  this%tRNO3MicrbTransfSoil=0.0_r8
  this%tRH2PO4MicrbTransfSoil=0.0_r8
  this%tRH1PO4MicrbTransfSoil=0.0_r8
  this%tRNH4MicrbTransfBand=0.0_r8
  this%tRNO3MicrbTransfBand=0.0_r8
  this%tRH2PO4MicrbTransfBand=0.0_r8
  this%tRH1PO4MicrbTransfBand=0.0_r8
  this%tRCO2MicrbProd=0.0_r8
  this%tRCH4MicrbProd=0.0_r8
  this%tRNOxMicrbRedux=0.0_r8
  this%tRCO2GrothAutor=0.0_r8
  this%TProdH2=0.0_r8
  this%tRO2MicrbUptk=0.0_r8
  this%TReduxNO3Soil=0.0_r8
  this%TReduxNO3Band=0.0_r8
  this%TReduxNO2Soil=0.0_r8
  this%TReduxNO2Band=0.0_r8
  this%TReduxN2O=0.0_r8
  this%TFixN2=0.0_r8
  this%tCH4OxiAero=0._r8
  end subroutine nit_aqmf_diag
!------------------------------------------------------------------------------------------

  subroutine nit_micf_init(this,jcplx,NumMicbFunGrupsPerCmplx)
  implicit none
  class(Microbe_Flux_type) :: this
  integer, intent(in) :: jcplx,NumMicbFunGrupsPerCmplx
  integer :: ndbiomcp
  integer :: NumMicrobAutrophCmplx
  integer :: NumHetetrMicCmplx
  ndbiomcp=micpar%ndbiomcp
  NumMicrobAutrophCmplx=micpar%NumMicrobAutrophCmplx
  NumHetetrMicCmplx=micpar%NumHetetrMicCmplx

  allocate(this%RO2UptkHeter(NumHetetrMicCmplx,1:jcplx));this%RO2UptkHeter=spval
  allocate(this%Resp4NFixHeter(NumHetetrMicCmplx,1:jcplx));this%Resp4NFixHeter=spval
  allocate(this%RespGrossHeter(NumHetetrMicCmplx,1:jcplx));this%RespGrossHeter=spval
  allocate(this%RO2Dmnd4RespHeter(NumHetetrMicCmplx,1:jcplx));this%RO2Dmnd4RespHeter=spval
  allocate(this%RO2DmndHeter(NumHetetrMicCmplx,1:jcplx));this%RO2DmndHeter=spval
  allocate(this%RO2Uptk4RespHeter(NumHetetrMicCmplx,1:jcplx));this%RO2Uptk4RespHeter=spval
  allocate(this%RNO3ReduxHeterSoil(NumHetetrMicCmplx,1:jcplx));this%RNO3ReduxHeterSoil=spval
  allocate(this%RNO3ReduxHeterBand(NumHetetrMicCmplx,1:jcplx));this%RNO3ReduxHeterBand=spval
  allocate(this%RNO2ReduxHeterSoil(NumHetetrMicCmplx,1:jcplx));this%RNO2ReduxHeterSoil=spval
  allocate(this%RNO2ReduxHeterBand(NumHetetrMicCmplx,1:jcplx));this%RNO2ReduxHeterBand=spval
  allocate(this%RN2OReduxHeter(NumHetetrMicCmplx,1:jcplx));this%RN2OReduxHeter=spval
  allocate(this%RNOxReduxRespDenitLim(NumHetetrMicCmplx,1:jcplx));this%RNOxReduxRespDenitLim=spval
  allocate(this%RMaintDmndHeter(2,NumHetetrMicCmplx,1:jcplx));this%RMaintDmndHeter=spval
  allocate(this%RNH4TransfSoilHeter(NumHetetrMicCmplx,1:jcplx));this%RNH4TransfSoilHeter=spval
  allocate(this%RNO3TransfSoilHeter(NumHetetrMicCmplx,1:jcplx));this%RNO3TransfSoilHeter=spval
  allocate(this%RH2PO4TransfSoilHeter(NumHetetrMicCmplx,1:jcplx));this%RH2PO4TransfSoilHeter=spval
  allocate(this%RNH4TransfBandHeter(NumHetetrMicCmplx,1:jcplx));this%RNH4TransfBandHeter=spval
  allocate(this%RNO3TransfBandHeter(NumHetetrMicCmplx,1:jcplx));this%RNO3TransfBandHeter=spval
  allocate(this%RH2PO4TransfBandHeter(NumHetetrMicCmplx,1:jcplx));this%RH2PO4TransfBandHeter=spval
  allocate(this%RkillLitfalOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RkillLitfalOMHeter=spval
  allocate(this%RkillLitrfal2HumOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RkillLitrfal2HumOMHeter=spval
  allocate(this%RkillLitrfal2ResduOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RkillLitrfal2ResduOMHeter=spval
  allocate(this%DOMuptk4GrothHeter(NumPlantChemElms,NumHetetrMicCmplx,1:jcplx));this%DOMuptk4GrothHeter=spval
  allocate(this%RH2ProdHeter(NumHetetrMicCmplx,1:jcplx));this%RH2ProdHeter=spval
  allocate(this%RMaintDefcitLitrfalOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RMaintDefcitLitrfalOMHeter=spval
  allocate(this%RMaintDefcitLitrfal2HumOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RMaintDefcitLitrfal2HumOMHeter=spval
  allocate(this%RMaintDefcitLitrfal2ResduOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RMaintDefcitLitrfal2ResduOMHeter=spval
  allocate(this%RN2FixHeter(NumHetetrMicCmplx,1:jcplx));this%RN2FixHeter=spval
  allocate(this%RKillOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RKillOMHeter=spval
  allocate(this%RkillRecycOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RkillRecycOMHeter=spval
  allocate(this%RMaintDefcitKillOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RMaintDefcitKillOMHeter=spval
  allocate(this%RMaintDefcitRecycOMHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%RMaintDefcitRecycOMHeter=spval
  allocate(this%RAnabolDOCUptkHeter(NumHetetrMicCmplx,1:jcplx));this%RAnabolDOCUptkHeter=spval
  allocate(this%RAnabolAcetUptkHeter(NumHetetrMicCmplx,1:jcplx));this%RAnabolAcetUptkHeter=spval
  allocate(this%RNH4TransfLitrHeter(NumHetetrMicCmplx,1:jcplx));this%RNH4TransfLitrHeter=spval
  allocate(this%RNO3TransfLitrHeter(NumHetetrMicCmplx,1:jcplx));this%RNO3TransfLitrHeter=spval
  allocate(this%RH2PO4TransfLitrHeter(NumHetetrMicCmplx,1:jcplx));this%RH2PO4TransfLitrHeter=spval
  allocate(this%AttenfNH4Heter(NumHetetrMicCmplx,1:jcplx));this%AttenfNH4Heter=spval
  allocate(this%AttenfNO3Heter(NumHetetrMicCmplx,1:jcplx));this%AttenfNO3Heter=spval
  allocate(this%AttenfH2PO4Heter(NumHetetrMicCmplx,1:jcplx));this%AttenfH2PO4Heter=spval
  allocate(this%RNOxReduxRespDenitUlm(NumHetetrMicCmplx,1:jcplx));this%RNOxReduxRespDenitUlm=spval
  allocate(this%NonstX2stBiomHeter(NumPlantChemElms,2,NumHetetrMicCmplx,1:jcplx));this%NonstX2stBiomHeter=spval
  allocate(this%AttenfH1PO4Heter(NumHetetrMicCmplx,1:jcplx));this%AttenfH1PO4Heter=spval
  allocate(this%RCO2ProdHeter(NumHetetrMicCmplx,1:jcplx));this%RCO2ProdHeter=spval
  allocate(this%RAcettProdHeter(NumHetetrMicCmplx,1:jcplx));this%RAcettProdHeter=spval
  allocate(this%RCH4ProdHeter(NumHetetrMicCmplx,1:jcplx));this%RCH4ProdHeter=spval
  allocate(this%RH1PO4TransfSoilHeter(NumHetetrMicCmplx,1:jcplx));this%RH1PO4TransfSoilHeter=spval
  allocate(this%RH1PO4TransfBandHeter(NumHetetrMicCmplx,1:jcplx));this%RH1PO4TransfBandHeter=spval
  allocate(this%RH1PO4TransfLitrHeter(NumHetetrMicCmplx,1:jcplx));this%RH1PO4TransfLitrHeter=spval

  allocate(this%RSOxidSoilAutor(NumMicrobAutrophCmplx));this%RSOxidSoilAutor=spval
  allocate(this%RSOxidBandAutor(NumMicrobAutrophCmplx));this%RSOxidBandAutor=spval
  allocate(this%XferBiomeHeterK(1:NumPlantChemElms,3,NumHetetrMicCmplx,1:jcplx));this%XferBiomeHeterK=spval
  allocate(this%ROQC4HeterMicrobAct(NumHetetrMicCmplx,1:jcplx));this%ROQC4HeterMicrobAct=spval
  allocate(this%RCCMEheter(NumPlantChemElms,ndbiomcp,NumHetetrMicCmplx,1:jcplx));this%RCCMEheter=spval

  allocate(this%RO2UptkAutor(NumMicrobAutrophCmplx));this%RO2UptkAutor=spval
  allocate(this%Resp4NFixAutor(NumMicrobAutrophCmplx));this%Resp4NFixAutor=spval
  allocate(this%RespGrossAutor(NumMicrobAutrophCmplx));this%RespGrossAutor=spval
  allocate(this%RO2Dmnd4RespAutor(NumMicrobAutrophCmplx));this%RO2Dmnd4RespAutor=spval
  allocate(this%RO2DmndAutor(NumMicrobAutrophCmplx));this%RO2DmndAutor=spval
  allocate(this%RO2Uptk4RespAutor(NumMicrobAutrophCmplx));this%RO2Uptk4RespAutor=spval
  allocate(this%RNO3UptkAutor(NumMicrobAutrophCmplx));this%RNO3UptkAutor=spval
  allocate(this%RNO2ReduxAutorSoil(NumMicrobAutrophCmplx));this%RNO2ReduxAutorSoil=spval
  allocate(this%RNO2ReduxAutorBand(NumMicrobAutrophCmplx));this%RNO2ReduxAutorBand=spval
  allocate(this%RNOxReduxRespAutorLim(NumMicrobAutrophCmplx));this%RNOxReduxRespAutorLim=spval
  allocate(this%RMaintDmndAutor(2,NumMicrobAutrophCmplx));this%RMaintDmndAutor=spval
  allocate(this%RNH4TransfSoilAutor(NumMicrobAutrophCmplx));this%RNH4TransfSoilAutor=spval
  allocate(this%RNO3TransfSoilAutor(NumMicrobAutrophCmplx));this%RNO3TransfSoilAutor=spval
  allocate(this%RH2PO4TransfSoilAutor(NumMicrobAutrophCmplx));this%RH2PO4TransfSoilAutor=spval
  allocate(this%RNH4TransfBandAutor(NumMicrobAutrophCmplx));this%RNH4TransfBandAutor=spval
  allocate(this%RNO3TransfBandAutor(NumMicrobAutrophCmplx));this%RNO3TransfBandAutor=spval
  allocate(this%RH2PO4TransfBandAutor(NumMicrobAutrophCmplx));this%RH2PO4TransfBandAutor=spval
  allocate(this%RkillLitfalOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RkillLitfalOMAutor=spval
  allocate(this%RkillLitrfal2HumOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RkillLitrfal2HumOMAutor=spval
  allocate(this%RkillLitrfal2ResduOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RkillLitrfal2ResduOMAutor=spval
  allocate(this%DOMuptk4GrothAutor(idom_beg:idom_end,NumMicrobAutrophCmplx));this%DOMuptk4GrothAutor=spval
  allocate(this%RMaintDefcitLitrfalOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RMaintDefcitLitrfalOMAutor=spval
  allocate(this%RMaintDefcitLitrfal2HumOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RMaintDefcitLitrfal2HumOMAutor=spval
  allocate(this%RMaintDefcitLitrfal2ResduOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RMaintDefcitLitrfal2ResduOMAutor=spval
  allocate(this%RN2FixAutor(NumMicrobAutrophCmplx));this%RN2FixAutor=spval
  allocate(this%RKillOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RKillOMAutor=spval
  allocate(this%RkillRecycOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RkillRecycOMAutor=spval
  allocate(this%RMaintDefcitKillOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RMaintDefcitKillOMAutor=spval
  allocate(this%RMaintDefcitRecycOMAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%RMaintDefcitRecycOMAutor=spval
  allocate(this%RNH4TransfLitrAutor(NumMicrobAutrophCmplx));this%RNH4TransfLitrAutor=spval
  allocate(this%RNO3TransfLitrAutor(NumMicrobAutrophCmplx));this%RNO3TransfLitrAutor=spval
  allocate(this%RH2PO4TransfLitrAutor(NumMicrobAutrophCmplx));this%RH2PO4TransfLitrAutor=spval
  allocate(this%AttenfNH4Autor(NumMicrobAutrophCmplx));this%AttenfNH4Autor=spval
  allocate(this%AttenfNO3Autor(NumMicrobAutrophCmplx));this%AttenfNO3Autor=spval
  allocate(this%AttenfH2PO4Autor(NumMicrobAutrophCmplx));this%AttenfH2PO4Autor=spval
  allocate(this%NonstX2stBiomAutor(NumPlantChemElms,2,NumMicrobAutrophCmplx));this%NonstX2stBiomAutor=spval
  allocate(this%AttenfH1PO4Autor(NumMicrobAutrophCmplx));this%AttenfH1PO4Autor=spval
  allocate(this%RCO2ProdAutor(NumMicrobAutrophCmplx));this%RCO2ProdAutor=spval
  allocate(this%RCH4ProdAutor(NumMicrobAutrophCmplx));this%RCH4ProdAutor=spval
  allocate(this%RH1PO4TransfSoilAutor(NumMicrobAutrophCmplx));this%RH1PO4TransfSoilAutor=spval
  allocate(this%RH1PO4TransfBandAutor(NumMicrobAutrophCmplx));this%RH1PO4TransfBandAutor=spval
  allocate(this%RH1PO4TransfLitrAutor(NumMicrobAutrophCmplx));this%RH1PO4TransfLitrAutor=spval

  call this%ZeroOut()
  end subroutine nit_micf_init
!------------------------------------------------------------------------------------------

  subroutine nit_mics_init(this, jcplx,NumMicbFunGrupsPerCmplx)

  implicit none
  class(Microbe_State_type) :: this
  integer, intent(in) :: jcplx,NumMicbFunGrupsPerCmplx
  integer :: NumMicrobAutrophCmplx,NumHetetrMicCmplx
  NumMicrobAutrophCmplx=micpar%NumMicrobAutrophCmplx
  NumHetetrMicCmplx=micpar%NumHetetrMicCmplx

  allocate(this%rCNBiomeActHeter(NumPlantChemElms,NumHetetrMicCmplx,1:jcplx));this%rCNBiomeActHeter=1._r8
  allocate(this%OMActHeter(NumHetetrMicCmplx,1:jcplx));this%OMActHeter=spval
  allocate(this%FracOMActHeter(NumHetetrMicCmplx,1:jcplx));this%FracOMActHeter=spval
  allocate(this%FracNO2ReduxHeter(NumHetetrMicCmplx,1:jcplx));this%FracNO2ReduxHeter=spval
  allocate(this%FracHeterBiomOfActK(NumHetetrMicCmplx,1:jcplx));this%FracHeterBiomOfActK=spval
  allocate(this%OMC2(NumHetetrMicCmplx,1:jcplx));this%OMC2=spval
  allocate(this%GrowthEnvScalHeter(NumHetetrMicCmplx,1:jcplx));this%GrowthEnvScalHeter=spval
  allocate(this%TempMaintRHeter(NumHetetrMicCmplx,1:jcplx));this%TempMaintRHeter=spval
  allocate(this%OMN2(NumHetetrMicCmplx,1:jcplx));this%OMN2=spval
  allocate(this%FOM2(NumHetetrMicCmplx,1:jcplx));this%FOM2=spval
  allocate(this%OxyLimterHeter(NumHetetrMicCmplx,1:jcplx));this%OxyLimterHeter=spval
  allocate(this%FCN(NumHetetrMicCmplx,1:jcplx));this%FCN=spval
  allocate(this%FCP(NumHetetrMicCmplx,1:jcplx));this%FCP=spval
  allocate(this%FBiomStoiScalarHeter(NumHetetrMicCmplx,1:jcplx));this%FBiomStoiScalarHeter=spval

  allocate(this%rCNBiomeActAutor(NumPlantChemElms,NumMicrobAutrophCmplx));this%rCNBiomeActAutor=1.0_r8
  allocate(this%OMActAutor(NumMicrobAutrophCmplx));this%OMActAutor=spval
  allocate(this%FracOMActAutor(NumMicrobAutrophCmplx));this%FracOMActAutor=spval
  allocate(this%FracNO2ReduxAutor(NumMicrobAutrophCmplx));this%FracNO2ReduxAutor=spval
  allocate(this%FracAutorBiomOfActK(NumMicrobAutrophCmplx));this%FracAutorBiomOfActK=spval
  allocate(this%OMC2Autor(NumMicrobAutrophCmplx));this%OMC2Autor=spval
  allocate(this%GrowthEnvScalAutor(NumMicrobAutrophCmplx));this%GrowthEnvScalAutor=spval
  allocate(this%TSensMaintRAutor(NumMicrobAutrophCmplx));this%TSensMaintRAutor=spval
  allocate(this%OMN2Autor(NumMicrobAutrophCmplx));this%OMN2Autor=spval
  allocate(this%FOM2Autor(NumMicrobAutrophCmplx));this%FOM2Autor=spval
  allocate(this%fLimO2Autor(NumMicrobAutrophCmplx));this%fLimO2Autor=spval
  allocate(this%FCNAutor(NumMicrobAutrophCmplx));this%FCNAutor=spval
  allocate(this%FCPAutor(NumMicrobAutrophCmplx));this%FCPAutor=spval
  allocate(this%FBiomStoiScalarAutor(NumMicrobAutrophCmplx));this%FBiomStoiScalarAutor=spval
  end subroutine nit_mics_init
!------------------------------------------------------------------------------------------

  subroutine nit_micf_zero(this)
  implicit none
  class(Microbe_Flux_type) :: this

  this%RO2UptkHeter = 0._r8
  this%Resp4NFixHeter = 0._r8
  this%RespGrossHeter = 0._r8
  this%RO2Dmnd4RespHeter = 0._r8
  this%RO2DmndHeter = 0._r8
  this%RO2Uptk4RespHeter = 0._r8
  this%RNO3ReduxHeterSoil = 0._r8
  this%RNO3ReduxHeterBand = 0._r8
  this%RNO2ReduxHeterSoil = 0._r8
  this%RNO2ReduxHeterBand = 0._r8
  this%RN2OReduxHeter = 0._r8
  this%RNOxReduxRespDenitLim = 0._r8
  this%RMaintDmndHeter = 0._r8
  this%RNH4TransfSoilHeter = 0._r8
  this%RNO3TransfSoilHeter = 0._r8
  this%RH2PO4TransfSoilHeter = 0._r8
  this%RNH4TransfBandHeter = 0._r8
  this%RNO3TransfBandHeter = 0._r8
  this%RH2PO4TransfBandHeter = 0._r8
  this%RkillLitfalOMHeter = 0._r8
  this%RkillLitrfal2HumOMHeter = 0._r8
  this%RkillLitrfal2ResduOMHeter = 0._r8
  this%DOMuptk4GrothHeter = 0._r8
  this%RH2ProdHeter = 0._r8
  this%RMaintDefcitLitrfal2HumOMHeter = 0._r8
  this%RMaintDefcitLitrfal2ResduOMHeter = 0._r8
  this%RMaintDefcitLitrfalOMHeter = 0._r8
  this%RCCMEheter = 0._r8
  this%RN2FixHeter = 0._r8

  this%RKillOMHeter = 0._r8
  this%RkillRecycOMHeter = 0._r8
  this%RMaintDefcitKillOMHeter = 0._r8
  this%RMaintDefcitRecycOMHeter = 0._r8
  this%RAnabolDOCUptkHeter = 0._r8
  this%RAnabolAcetUptkHeter = 0._r8
  this%RNH4TransfLitrHeter = 0._r8
  this%RNO3TransfLitrHeter = 0._r8
  this%RH2PO4TransfLitrHeter = 0._r8
  this%AttenfNH4Heter = 0._r8
  this%AttenfNO3Heter = 0._r8
  this%AttenfH2PO4Heter = 0._r8
  this%RNOxReduxRespDenitUlm = 0._r8
  this%ROQC4HeterMicrobAct = 0._r8
  this%NonstX2stBiomHeter = 0._r8
  this%AttenfH1PO4Heter = 0._r8
  this%RCO2ProdHeter = 0._r8
  this%RAcettProdHeter = 0._r8
  this%RCH4ProdHeter = 0._r8
  this%RSOxidSoilAutor = 0._r8
  this%RSOxidBandAutor = 0._r8
  this%XferBiomeHeterK = 0._r8
  this%RH1PO4TransfSoilHeter = 0._r8
  this%RH1PO4TransfBandHeter = 0._r8
  this%RH1PO4TransfLitrHeter = 0._r8

  this%RO2UptkAutor = 0._r8
  this%Resp4NFixAutor = 0._r8
  this%RespGrossAutor = 0._r8
  this%RO2Dmnd4RespAutor = 0._r8
  this%RO2DmndAutor = 0._r8
  this%RO2Uptk4RespAutor = 0._r8
  this%RNO3UptkAutor = 0._r8
  this%RNO2ReduxAutorSoil = 0._r8
  this%RNO2ReduxAutorBand = 0._r8
  this%RNOxReduxRespAutorLim = 0._r8
  this%RMaintDmndAutor = 0._r8
  this%RNH4TransfSoilAutor = 0._r8
  this%RNO3TransfSoilAutor = 0._r8
  this%RH2PO4TransfSoilAutor = 0._r8
  this%RNH4TransfBandAutor = 0._r8
  this%RNO3TransfBandAutor = 0._r8
  this%RH2PO4TransfBandAutor = 0._r8
  this%RkillLitfalOMAutor = 0._r8
  this%RkillLitrfal2HumOMAutor = 0._r8
  this%RkillLitrfal2ResduOMAutor = 0._r8
  this%DOMuptk4GrothAutor = 0._r8
  this%RMaintDefcitLitrfalOMAutor = 0._r8
  this%RMaintDefcitLitrfal2HumOMAutor = 0._r8
  this%RMaintDefcitLitrfal2ResduOMAutor = 0._r8
  this%RN2FixAutor = 0._r8

  this%RKillOMAutor = 0._r8
  this%RkillRecycOMAutor = 0._r8
  this%RMaintDefcitKillOMAutor = 0._r8
  this%RMaintDefcitRecycOMAutor = 0._r8
  this%RNH4TransfLitrAutor = 0._r8
  this%RNO3TransfLitrAutor = 0._r8
  this%RH2PO4TransfLitrAutor = 0._r8
  this%AttenfNH4Autor = 0._r8
  this%AttenfNO3Autor = 0._r8
  this%AttenfH2PO4Autor = 0._r8
  this%NonstX2stBiomAutor = 0._r8
  this%AttenfH1PO4Autor = 0._r8
  this%RCO2ProdAutor = 0._r8
  this%RCH4ProdAutor = 0._r8
  this%RH1PO4TransfSoilAutor = 0._r8
  this%RH1PO4TransfBandAutor = 0._r8
  this%RH1PO4TransfLitrAutor = 0._r8

  end subroutine nit_micf_zero


!------------------------------------------------------------------------------------------

  subroutine nit_micf_destroy(this)

  implicit none
  class(Microbe_Flux_type) :: this

  call destroy(this%RO2UptkHeter)
  call destroy(this%Resp4NFixHeter)
  call destroy(this%RespGrossHeter)
  call destroy(this%RO2Dmnd4RespHeter)
  call destroy(this%RO2DmndHeter)
  call destroy(this%RO2Uptk4RespHeter)
  call destroy(this%RNO3ReduxHeterSoil)
  call destroy(this%RNO3ReduxHeterBand)
  call destroy(this%RNO2ReduxHeterSoil)
  call destroy(this%RNO2ReduxHeterBand)
  call destroy(this%RN2OReduxHeter)
  call destroy(this%RNOxReduxRespDenitLim)
  call destroy(this%RMaintDmndHeter)
  call destroy(this%RNH4TransfSoilHeter)
  call destroy(this%RNO3TransfSoilHeter)
  call destroy(this%RH2PO4TransfSoilHeter)
  call destroy(this%RNH4TransfBandHeter)
  call destroy(this%RNO3TransfBandHeter)
  call destroy(this%RH2PO4TransfBandHeter)
  call destroy(this%RkillLitfalOMHeter)
  call destroy(this%RkillLitrfal2HumOMHeter)
  call destroy(this%RkillLitrfal2ResduOMHeter)
  call destroy(this%DOMuptk4GrothHeter)
  call destroy(this%RH2ProdHeter)
  call destroy(this%RMaintDefcitLitrfal2HumOMHeter)
  call destroy(this%RMaintDefcitLitrfal2ResduOMHeter)
  call destroy(this%RMaintDefcitLitrfalOMHeter)
  call destroy(this%RCCMEheter)
  call destroy(this%RN2FixHeter)
  call destroy(this%RKillOMHeter)
  call destroy(this%RkillRecycOMHeter)
  call destroy(this%RMaintDefcitKillOMHeter)
  call destroy(this%RMaintDefcitRecycOMHeter)
  call destroy(this%RAnabolDOCUptkHeter)
  call destroy(this%RAnabolAcetUptkHeter)
  call destroy(this%RNH4TransfLitrHeter)
  call destroy(this%RNO3TransfLitrHeter)
  call destroy(this%RH2PO4TransfLitrHeter)
  call destroy(this%AttenfNH4Heter)
  call destroy(this%AttenfNO3Heter)
  call destroy(this%AttenfH2PO4Heter)
  call destroy(this%RNOxReduxRespDenitUlm)
  call destroy(this%ROQC4HeterMicrobAct)
  call destroy(this%NonstX2stBiomHeter)
  call destroy(this%AttenfH1PO4Heter)
  call destroy(this%RCO2ProdHeter)
  call destroy(this%RAcettProdHeter)
  call destroy(this%RCH4ProdHeter)
  call destroy(this%RSOxidSoilAutor)
  call destroy(this%RSOxidBandAutor)
  call destroy(this%XferBiomeHeterK)
  call destroy(this%RH1PO4TransfSoilHeter)
  call destroy(this%RH1PO4TransfBandHeter)
  call destroy(this%RH1PO4TransfLitrHeter)

  call destroy(this%RO2UptkAutor)
  call destroy(this%Resp4NFixAutor)
  call destroy(this%RespGrossAutor)
  call destroy(this%RO2Dmnd4RespAutor)
  call destroy(this%RO2DmndAutor)
  call destroy(this%RO2Uptk4RespAutor)
  call destroy(this%RNO3UptkAutor)
  call destroy(this%RNO2ReduxAutorSoil)
  call destroy(this%RNO2ReduxAutorBand)
  call destroy(this%RNOxReduxRespAutorLim)
  call destroy(this%RMaintDmndAutor)
  call destroy(this%RNH4TransfSoilAutor)
  call destroy(this%RNO3TransfSoilAutor)
  call destroy(this%RH2PO4TransfSoilAutor)
  call destroy(this%RNH4TransfBandAutor)
  call destroy(this%RNO3TransfBandAutor)
  call destroy(this%RH2PO4TransfBandAutor)
  call destroy(this%RkillLitfalOMAutor)
  call destroy(this%RkillLitrfal2HumOMAutor)
  call destroy(this%RkillLitrfal2ResduOMAutor)
  call destroy(this%DOMuptk4GrothAutor)
  call destroy(this%RMaintDefcitLitrfalOMAutor)
  call destroy(this%RMaintDefcitLitrfal2HumOMAutor)
  call destroy(this%RMaintDefcitLitrfal2ResduOMAutor)
  call destroy(this%RN2FixAutor)
  call destroy(this%RKillOMAutor)
  call destroy(this%RkillRecycOMAutor)
  call destroy(this%RMaintDefcitKillOMAutor)
  call destroy(this%RMaintDefcitRecycOMAutor)
  call destroy(this%RNH4TransfLitrAutor)
  call destroy(this%RNO3TransfLitrAutor)
  call destroy(this%RH2PO4TransfLitrAutor)
  call destroy(this%AttenfNH4Autor)
  call destroy(this%AttenfNO3Autor)
  call destroy(this%AttenfH2PO4Autor)
  call destroy(this%NonstX2stBiomAutor)
  call destroy(this%AttenfH1PO4Autor)
  call destroy(this%RCO2ProdAutor)
  call destroy(this%RCH4ProdAutor)
  call destroy(this%RH1PO4TransfSoilAutor)
  call destroy(this%RH1PO4TransfBandAutor)
  call destroy(this%RH1PO4TransfLitrAutor)

  end subroutine nit_micf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_mics_destroy(this)
  implicit none
  class(Microbe_State_type) :: this

  call destroy(this%rCNBiomeActHeter)
  call destroy(this%OMActHeter)
  call destroy(this%FracOMActHeter)
  call destroy(this%FracNO2ReduxHeter)
  call destroy(this%FracHeterBiomOfActK)
  call destroy(this%OMC2)
  call destroy(this%GrowthEnvScalHeter)
  call destroy(this%TempMaintRHeter)
  call destroy(this%OMN2)
  call destroy(this%FOM2)
  call destroy(this%OxyLimterHeter)
  call destroy(this%FCN)
  call destroy(this%FCP)
  call destroy(this%FBiomStoiScalarHeter)

  call destroy(this%rCNBiomeActAutor)
  call destroy(this%OMActAutor)
  call destroy(this%FracOMActAutor)
  call destroy(this%FracNO2ReduxAutor)
  call destroy(this%FracAutorBiomOfActK)
  call destroy(this%OMC2Autor)
  call destroy(this%GrowthEnvScalAutor)
  call destroy(this%TSensMaintRAutor)
  call destroy(this%OMN2Autor)
  call destroy(this%FOM2Autor)
  call destroy(this%fLimO2Autor)
  call destroy(this%FCNAutor)
  call destroy(this%FCPAutor)
  call destroy(this%FBiomStoiScalarAutor)
  end subroutine nit_mics_destroy
!------------------------------------------------------------------------------------------
  subroutine nit_omcplxf_init(this)
  implicit none
  class(OMCplx_Flux_type) :: this
  integer :: nkinets
  integer :: ncplx
  integer :: ndbiomcp

  nkinets=micpar%jsken
  ncplx=micpar%jcplx
  ndbiomcp=micpar%ndbiomcp
  allocate(this%RHydlysSolidOM(NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RHumifySolidOM(1:NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RDcmpProdDOM(1:NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RHydlysBioResduOM(1:NumPlantChemElms,ndbiomcp,1:ncplx))
  allocate(this%RHydlysSorptOM(idom_beg:idom_end,1:ncplx))
  allocate(this%RDOMSorp(idom_beg:idom_end,1:ncplx))
  allocate(this%TDOMUptkHeter(idom_beg:idom_end,1:ncplx+1))
  allocate(this%ROQC4HeterMicActCmpK(1:ncplx))
  allocate(this%XferRespHeterK(1:ncplx))
  allocate(this%XferDOMK(idom_beg:idom_end,1:ncplx))

  call this%ZeroOut()
  end subroutine nit_omcplxf_init
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_destroy(this)
  implicit none
  class(OMCplx_Flux_type) :: this

  call destroy(this%RHydlysSolidOM)
  call destroy(this%RHumifySolidOM)
  call destroy(this%RDcmpProdDOM)
  call destroy(this%RHydlysBioResduOM)
  call destroy(this%RHydlysSorptOM)
  call destroy(this%RDOMSorp)
  call destroy(this%ROQC4HeterMicActCmpK)
  call destroy(this%XferRespHeterK)
  call destroy(this%XferDOMK)

  end subroutine nit_omcplxf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_zero(this)
  implicit none
  class(OMCplx_Flux_type) :: this

  this%RHydlysSolidOM=0._r8
  this%RHumifySolidOM=0._r8
  this%RDcmpProdDOM=0._r8
  this%RHydlysBioResduOM=0._r8
  this%RHydlysSorptOM=0._r8
  this%RDOMSorp=0._r8
  this%ROQC4HeterMicActCmpK=0._r8
  this%XferRespHeterK=0._r8
  this%XferDOMK=0._r8
  end subroutine nit_omcplxf_zero
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxs_init(this)
  implicit none
  class(OMCplx_State_type) :: this
  integer :: ncplx

  ncplx=micpar%jcplx
  allocate(this%BulkSOMC(1:ncplx));this%BulkSOMC=spval
  allocate(this%TOMEAutoK(1:NumPlantChemElms)); this%TOMEAutoK=spval
  allocate(this%TOMK(1:ncplx+1));this%TOMK=spval
  allocate(this%TONK(1:ncplx+1));this%TONK=spval
  allocate(this%TOPK(1:ncplx+1));this%TOPK=spval
  allocate(this%FOCA(1:ncplx));this%FOCA=spval
  allocate(this%FOAA(1:ncplx));this%FOAA=spval
  allocate(this%rCNDOM(1:ncplx));this%rCNDOM=spval
  allocate(this%rCPDOM(1:ncplx));this%rCPDOM=spval
  allocate(this%rCNSorbOM(1:ncplx));this%rCNSorbOM=spval
  allocate(this%rCPSorbOM(1:ncplx));this%rCPSorbOM=spval
  allocate(this%OMBioResduK(1:ncplx));this%OMBioResduK=spval
  allocate(this%SolidOMK(1:ncplx));this%SolidOMK=spval
  allocate(this%SolidOMActK(1:ncplx));this%SolidOMActK=spval
  allocate(this%tMaxNActMicrbK(1:ncplx+1));this%tMaxNActMicrbK=spval
  allocate(this%tMaxPActMicrbK(1:ncplx+1));this%tMaxPActMicrbK=spval
  allocate(this%CDOM(idom_beg:idom_end,1:ncplx));this%CDOM=spval
  call this%ZeroOut()
  end subroutine nit_omcplxs_init
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxs_zero(this)

  implicit none
  class(OMCplx_State_type) :: this

  this%BulkSOMC=0._r8
  this%TOMK=0._r8
  this%TONK=0._r8
  this%TOPK=0._r8
  this%FOCA=0._r8
  this%FOAA=0._r8
  this%rCNDOM=0._r8
  this%rCPDOM=0._r8
  this%rCNSorbOM=0._r8
  this%rCPSorbOM=0._r8
  this%OMBioResduK=0._r8
  this%SolidOMK=0._r8
  this%SolidOMActK=0._r8
  this%tMaxNActMicrbK=0._r8
  this%tMaxPActMicrbK=0._r8

  end subroutine nit_omcplxs_zero


!------------------------------------------------------------------------------------------
  subroutine nit_omcplxs_destroy(this)
  implicit none
  class(OMCplx_State_type) :: this

  call destroy(this%BulkSOMC)
  call destroy(this%TOMK)
  call destroy(this%TONK)
  call destroy(this%TOPK)
  call destroy(this%FOCA)
  call destroy(this%FOAA)
  call destroy(this%rCNDOM)
  call destroy(this%rCPDOM)
  call destroy(this%rCNSorbOM)
  call destroy(this%rCPSorbOM)
  call destroy(this%OMBioResduK)
  call destroy(this%SolidOMK)
  call destroy(this%SolidOMActK)
  call destroy(this%tMaxNActMicrbK)
  call destroy(this%tMaxPActMicrbK)
  call destroy(this%CDOM)

  end subroutine nit_omcplxs_destroy
end module MicrobeDiagTypes
