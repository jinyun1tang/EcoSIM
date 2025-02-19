module TFlxTypeMod

  use GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc,NumMicbFunGrupsPerCmplx=>NumMicbFunGrupsPerCmplx
  use EcoSIMConfig, only : nlbiomcp=>NumLiveMicrbCompts,ndbiomcp=>NumDeadMicrbCompts
implicit none

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public

  real(r8),allocatable ::  trcs_TransptMicP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcs_TransptMacP_vr(:,:,:,:)                      !

  real(r8),allocatable ::  trcSalt_Flo2MicP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcSalt_Flo2MacP_vr(:,:,:,:)                      !

  real(r8),allocatable ::  TSandEros_col(:,:)                        !
  real(r8),allocatable ::  TSiltEros_col(:,:)                        !
  real(r8),allocatable ::  TCLAYEros_col(:,:)                        !
  real(r8),allocatable ::  TNH4Eros_col(:,:)                        !
  real(r8),allocatable ::  TNH3Eros_col(:,:)                        !
  real(r8),allocatable ::  TNUreaEros_col(:,:)                        !
  real(r8),allocatable ::  TNO3Eros_col(:,:)                        !
  real(r8),allocatable ::  TNH4ErosBand_col(:,:)                        !
  real(r8),allocatable ::  TNH3ErosBand_col(:,:)                        !
  real(r8),allocatable ::  TNUreaErosBand_col(:,:)                        !
  real(r8),allocatable ::  TNO3ErosBand_col(:,:)                        !
  real(r8),allocatable ::  trcx_TER_col(:,:,:)                         !
  real(r8),allocatable ::  trcp_TER_col(:,:,:)                        !
  real(r8),allocatable ::  tErosionSedmLoss_col(:,:)                        !
  real(r8),allocatable ::  TWatFlowCellMicP_vr(:,:,:)                        !
  real(r8),allocatable ::  TWatFlowCellMicPX_vr(:,:,:)                       !
  real(r8),allocatable ::  TWatFlowCellMacP_vr(:,:,:)                       !

  real(r8),allocatable ::  Gas_AdvDif_Flx_vr(:,:,:,:)                      !

  real(r8),allocatable ::  WatIceThawMicP_vr(:,:,:)                       !
  real(r8),allocatable ::  THeatSoiThaw_vr(:,:,:)                      !
  real(r8),allocatable ::  WatIceThawMacP_vr(:,:,:)                      !
  real(r8),allocatable ::  VLWatMicP2_vr(:,:,:)                       !local copy of soil water in micropore
  real(r8),allocatable ::  VLiceMicP2_vr(:,:,:)                       !local copy of soil ice in micropore
  real(r8),allocatable ::  VLWatMacP2_vr(:,:,:)                       !local copy of soil water in macropore
  real(r8),allocatable ::  VLiceMacP2_vr(:,:,:)                       !local copy of soil ice in macropore

  real(r8),allocatable :: TOMEERhetr_col(:,:,:,:,:)

  real(r8),allocatable :: TOMEERauto_col(:,:,:,:)
 
  real(r8),allocatable ::  DOM_Transp2Micp_vr(:,:,:,:,:)
  real(r8),allocatable ::  DOM_Transp2Macp_flx(:,:,:,:,:)
  real(r8),allocatable ::  TOMQRS_col(:,:,:,:)
  real(r8),allocatable ::  TORMER_col(:,:,:,:,:)
  real(r8),allocatable ::  TOHMER_col(:,:,:,:)
  real(r8),allocatable ::  TOSMER_col(:,:,:,:,:)
  real(r8),allocatable ::  TOSAER_col(:,:,:,:)

  real(r8) :: TDLYXF,TDayLenthPrevC,TDVOLI,TDORGC
  DATA TDORGC,TDayLenthPrevC/0.0,0.0/
  DATA TDVOLI,TDLYXF/0.0,0.0/


  contains

!------------------------------------------------------------------------------------------

  subroutine InitTflxType()

  implicit none

  allocate(trcs_TransptMacP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_TransptMacP_vr=0._r8
  allocate(Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,JZ,JY,JX));   Gas_AdvDif_Flx_vr=0._r8

  allocate(trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_Flo2MicP_vr=0._r8
  allocate(trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_Flo2MacP_vr=0._r8

  allocate(TSandEros_col(JY,JX));      TSandEros_col=0._r8
  allocate(TSiltEros_col(JY,JX));      TSiltEros_col=0._r8
  allocate(TCLAYEros_col(JY,JX));      TCLAYEros_col=0._r8
  allocate(TNH4Eros_col(JY,JX));      TNH4Eros_col=0._r8
  allocate(TNH3Eros_col(JY,JX));      TNH3Eros_col=0._r8
  allocate(TNUreaEros_col(JY,JX));      TNUreaEros_col=0._r8
  allocate(TNO3Eros_col(JY,JX));      TNO3Eros_col=0._r8
  allocate(TNH4ErosBand_col(JY,JX));      TNH4ErosBand_col=0._r8
  allocate(TNH3ErosBand_col(JY,JX));      TNH3ErosBand_col=0._r8
  allocate(TNUreaErosBand_col(JY,JX));      TNUreaErosBand_col=0._r8
  allocate(TNO3ErosBand_col(JY,JX));      TNO3ErosBand_col=0._r8

  allocate(trcx_TER_col(idx_beg:idx_end,JY,JX));    trcx_TER_col=0._r8
  allocate(trcp_TER_col(idsp_beg:idsp_end,JY,JX));      trcp_TER_col=0._r8
  allocate(trcs_TransptMicP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_TransptMicP_vr=0._r8

  allocate(tErosionSedmLoss_col(JY,JX));      tErosionSedmLoss_col=0._r8
  allocate(TWatFlowCellMicP_vr(JZ,JY,JX));     TWatFlowCellMicP_vr=0._r8
  allocate(TWatFlowCellMicPX_vr(JZ,JY,JX));    TWatFlowCellMicPX_vr=0._r8
  allocate(TWatFlowCellMacP_vr(JZ,JY,JX));    TWatFlowCellMacP_vr=0._r8
  allocate(WatIceThawMicP_vr(JZ,JY,JX));    WatIceThawMicP_vr=0._r8
  allocate(THeatSoiThaw_vr(JZ,JY,JX));   THeatSoiThaw_vr=0._r8
  allocate(WatIceThawMacP_vr(JZ,JY,JX));   WatIceThawMacP_vr=0._r8
  allocate(VLWatMicP2_vr(JZ,JY,JX));    VLWatMicP2_vr=0._r8
  allocate(VLiceMicP2_vr(JZ,JY,JX));    VLiceMicP2_vr=0._r8
  allocate(VLWatMacP2_vr(JZ,JY,JX));   VLWatMacP2_vr=0._r8
  allocate(VLiceMacP2_vr(JZ,JY,JX));   VLiceMacP2_vr=0._r8
  allocate(TOMEERhetr_col(NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,JY,JX)); TOMEERhetr_col=0._r8

  allocate(TOMEERauto_col(NumPlantChemElms,1:NumLiveAutoBioms,JY,JX));TOMEERauto_col=0._r8

  allocate(DOM_Transp2Micp_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Micp_vr=0._r8
  allocate(DOM_Transp2Macp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flx=0._r8
  allocate(TOMQRS_col(idom_beg:idom_end,1:jcplx,JY,JX));TOMQRS_col=0._r8
  allocate(TORMER_col(NumPlantChemElms,ndbiomcp,1:jcplx,JY,JX));TORMER_col=0._r8
  allocate(TOHMER_col(idom_beg:idom_end,1:jcplx,JY,JX));TOHMER_col=0._r8
  allocate(TOSMER_col(NumPlantChemElms,jsken,1:jcplx,JY,JX));TOSMER_col=0._r8
  allocate(TOSAER_col(jsken,1:jcplx,JY,JX));TOSAER_col=0._r8

  end subroutine InitTflxType


!------------------------------------------------------------------------------------------

  subroutine DestructTflxType

  implicit none

  call destroy(trcx_TER_col)
  call destroy(TOMEERhetr_col)
  call destroy(TOMEERauto_col)
  call destroy(TOMQRS_col)
  call destroy(TORMER_col)
  call destroy(TOHMER_col)
  call destroy(TOSMER_col)
  call destroy(TOSAER_col)
  call destroy(trcSalt_Flo2MicP_vr)
  call destroy(trcs_TransptMicP_vr)
  call destroy(trcSalt_Flo2MacP_vr)
!  call destroy(TFLWS)
!  call destroy(TFLWW)
!  call destroy(TFLWI)
!  call destroy(THFLWW)

  call destroy(trcs_TransptMacP_vr)

  call destroy(TSandEros_col)
  call destroy(TSiltEros_col)
  call destroy(TCLAYEros_col)
  call destroy(TNH4Eros_col)
  call destroy(TNH3Eros_col)
  call destroy(TNUreaEros_col)
  call destroy(TNO3Eros_col)
  call destroy(TNH4ErosBand_col)
  call destroy(TNH3ErosBand_col)
  call destroy(TNUreaErosBand_col)
  call destroy(TNO3ErosBand_col)
  call destroy(tErosionSedmLoss_col)
  call destroy(TWatFlowCellMicP_vr)
  call destroy(TWatFlowCellMicPX_vr)
  call destroy(TWatFlowCellMacP_vr)
  call destroy(WatIceThawMicP_vr)
  call destroy(THeatSoiThaw_vr)
  call destroy(WatIceThawMacP_vr)
  call destroy(VLWatMicP2_vr)
  call destroy(VLiceMicP2_vr)
  call destroy(VLWatMacP2_vr)
  call destroy(VLiceMacP2_vr)
  call destroy(DOM_Transp2Micp_vr)
  call destroy(DOM_Transp2Macp_flx)
  call destroy(trcp_TER_col)
  call destroy(Gas_AdvDif_Flx_vr)

  end subroutine DestructTflxType
end module TFlxTypeMod
