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

  real(r8),allocatable ::  trcs_Transp2MicP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcs_Transp2MacP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  TXGridSurfRunoff_2DH(:,:)                           !
  real(r8),allocatable ::  THeatXGridBySurfRunoff_2DH(:,:)                          !

  real(r8),allocatable ::  TCOQSS(:,:)                        !
  real(r8),allocatable ::  TCHQSS(:,:)                        !
  real(r8),allocatable ::  TOXQSS(:,:)                        !
  real(r8),allocatable ::  TNGQSS(:,:)                        !
  real(r8),allocatable ::  TN2QSS(:,:)                        !
  real(r8),allocatable ::  TN4QSS(:,:)                        !
  real(r8),allocatable ::  TN3QSS(:,:)                        !
  real(r8),allocatable ::  TNOQSS(:,:)                        !
  real(r8),allocatable ::  TPOQSS(:,:)                        !
  real(r8),allocatable ::  TP1QSS(:,:)                        !

  real(r8),allocatable ::  trcSalt_Flo2MicP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcSalt_Flo2MacP_vr(:,:,:,:)                      !

  real(r8),allocatable ::  TSANER(:,:)                        !
  real(r8),allocatable ::  TSILER(:,:)                        !
  real(r8),allocatable ::  TCLAER(:,:)                        !
  real(r8),allocatable ::  TNH4Eros_col(:,:)                        !
  real(r8),allocatable ::  TNH3Eros_col(:,:)                        !
  real(r8),allocatable ::  TNUreaEros_col(:,:)                        !
  real(r8),allocatable ::  TNO3Eros_col(:,:)                        !
  real(r8),allocatable ::  TNH4ErosBand_col(:,:)                        !
  real(r8),allocatable ::  TNH3ErosBand_col(:,:)                        !
  real(r8),allocatable ::  TNUreaErosBand_col(:,:)                        !
  real(r8),allocatable ::  TNO3ErosBand_col(:,:)                        !
  real(r8),allocatable ::  trcx_TER(:,:,:)                         !
  real(r8),allocatable ::  trcp_TER(:,:,:)                        !
  real(r8),allocatable ::  tErosionSedmLoss(:,:)                        !
  real(r8),allocatable ::  TWatFlowCellMicP_vr(:,:,:)                        !
  real(r8),allocatable ::  TWatFlowCellMicPX_vr(:,:,:)                       !
  real(r8),allocatable ::  TWaterFlowMacP_vr(:,:,:)                       !

  real(r8),allocatable ::  Gas_AdvDif_Flx_vr(:,:,:,:)                      !

  real(r8),allocatable ::  WatIceThawMicP_vr(:,:,:)                       !
  real(r8),allocatable ::  THeatSoiThaw_vr(:,:,:)                      !
  real(r8),allocatable ::  WatIceThawMacP_vr(:,:,:)                      !
  real(r8),allocatable ::  VLWatMicP1_vr(:,:,:)                       !
  real(r8),allocatable ::  VLiceMicP1_vr(:,:,:)                       !
  real(r8),allocatable ::  VLWatMacP1_vr(:,:,:)                      !
  real(r8),allocatable ::  VLiceMacP1_vr(:,:,:)                      !

  real(r8),allocatable :: TOMEERhetr(:,:,:,:,:)

  real(r8),allocatable :: TOMEERauto(:,:,:,:)
 
  real(r8),allocatable ::  DOM_Transp2Micp_vr(:,:,:,:,:)
  real(r8),allocatable ::  DOM_Transp2Macp_flx(:,:,:,:,:)
  real(r8),allocatable ::  TOMQRS(:,:,:,:)
  real(r8),allocatable ::  TORMER(:,:,:,:,:)
  real(r8),allocatable ::  TOHMER(:,:,:,:)
  real(r8),allocatable ::  TOSMER(:,:,:,:,:)
  real(r8),allocatable ::  TOSAER(:,:,:,:)

  real(r8) :: TDLYXF,TDayLenthPrevC,TDVOLI,TDORGC
  DATA TDORGC,TDayLenthPrevC/0.0,0.0/
  DATA TDVOLI,TDLYXF/0.0,0.0/


  contains

!------------------------------------------------------------------------------------------

  subroutine InitTflxType()

  implicit none

  allocate(trcs_Transp2MacP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_Transp2MacP_vr=0._r8
  allocate(Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,JZ,JY,JX));   Gas_AdvDif_Flx_vr=0._r8

  allocate(TXGridSurfRunoff_2DH(JY,JX));         TXGridSurfRunoff_2DH=0._r8
  allocate(THeatXGridBySurfRunoff_2DH(JY,JX));        THeatXGridBySurfRunoff_2DH=0._r8
 ! allocate(TXGridSurfRunoff_2DH(JY,JX));         TXGridSurfRunoff_2DH=0._r8
 ! allocate(THeatXGridBySurfRunoff_2DH(JY,JX));        THeatXGridBySurfRunoff_2DH=0._r8

!  allocate(TFLWS(JS,JY,JX));    TFLWS=0._r8
!  allocate(TFLWW(JS,JY,JX));    TFLWW=0._r8
!  allocate(TFLWI(JS,JY,JX));    TFLWI=0._r8
!  allocate(THFLWW(JS,JY,JX));   THFLWW=0._r8

  allocate(TCOQSS(JY,JX));      TCOQSS=0._r8
  allocate(TCHQSS(JY,JX));      TCHQSS=0._r8
  allocate(TOXQSS(JY,JX));      TOXQSS=0._r8
  allocate(TNGQSS(JY,JX));      TNGQSS=0._r8
  allocate(TN2QSS(JY,JX));      TN2QSS=0._r8
  allocate(TN4QSS(JY,JX));      TN4QSS=0._r8
  allocate(TN3QSS(JY,JX));      TN3QSS=0._r8
  allocate(TNOQSS(JY,JX));      TNOQSS=0._r8
  allocate(TPOQSS(JY,JX));      TPOQSS=0._r8
  allocate(TP1QSS(JY,JX));      TP1QSS=0._r8

  allocate(trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_Flo2MicP_vr=0._r8
  allocate(trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_Flo2MacP_vr=0._r8

  allocate(TSANER(JY,JX));      TSANER=0._r8
  allocate(TSILER(JY,JX));      TSILER=0._r8
  allocate(TCLAER(JY,JX));      TCLAER=0._r8
  allocate(TNH4Eros_col(JY,JX));      TNH4Eros_col=0._r8
  allocate(TNH3Eros_col(JY,JX));      TNH3Eros_col=0._r8
  allocate(TNUreaEros_col(JY,JX));      TNUreaEros_col=0._r8
  allocate(TNO3Eros_col(JY,JX));      TNO3Eros_col=0._r8
  allocate(TNH4ErosBand_col(JY,JX));      TNH4ErosBand_col=0._r8
  allocate(TNH3ErosBand_col(JY,JX));      TNH3ErosBand_col=0._r8
  allocate(TNUreaErosBand_col(JY,JX));      TNUreaErosBand_col=0._r8
  allocate(TNO3ErosBand_col(JY,JX));      TNO3ErosBand_col=0._r8

  allocate(trcx_TER(idx_beg:idx_end,JY,JX));    trcx_TER=0._r8
  allocate(trcp_TER(idsp_beg:idsp_end,JY,JX));      trcp_TER=0._r8
  allocate(trcs_Transp2MicP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_Transp2MicP_vr=0._r8

  allocate(tErosionSedmLoss(JY,JX));      tErosionSedmLoss=0._r8
  allocate(TWatFlowCellMicP_vr(JZ,JY,JX));     TWatFlowCellMicP_vr=0._r8
  allocate(TWatFlowCellMicPX_vr(JZ,JY,JX));    TWatFlowCellMicPX_vr=0._r8
  allocate(TWaterFlowMacP_vr(JZ,JY,JX));    TWaterFlowMacP_vr=0._r8
  allocate(WatIceThawMicP_vr(JZ,JY,JX));    WatIceThawMicP_vr=0._r8
  allocate(THeatSoiThaw_vr(JZ,JY,JX));   THeatSoiThaw_vr=0._r8
  allocate(WatIceThawMacP_vr(JZ,JY,JX));   WatIceThawMacP_vr=0._r8
  allocate(VLWatMicP1_vr(JZ,JY,JX));    VLWatMicP1_vr=0._r8
  allocate(VLiceMicP1_vr(JZ,JY,JX));    VLiceMicP1_vr=0._r8
  allocate(VLWatMacP1_vr(JZ,JY,JX));   VLWatMacP1_vr=0._r8
  allocate(VLiceMacP1_vr(JZ,JY,JX));   VLiceMacP1_vr=0._r8
  allocate(TOMEERhetr(NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,JY,JX)); TOMEERhetr=0._r8

  allocate(TOMEERauto(NumPlantChemElms,1:NumLiveAutoBioms,JY,JX));TOMEERauto=0._r8

  allocate(DOM_Transp2Micp_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Micp_vr=0._r8
  allocate(DOM_Transp2Macp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flx=0._r8
  allocate(TOMQRS(idom_beg:idom_end,1:jcplx,JY,JX));TOMQRS=0._r8
  allocate(TORMER(NumPlantChemElms,ndbiomcp,1:jcplx,JY,JX));TORMER=0._r8
  allocate(TOHMER(idom_beg:idom_end,1:jcplx,JY,JX));TOHMER=0._r8
  allocate(TOSMER(NumPlantChemElms,jsken,1:jcplx,JY,JX));TOSMER=0._r8
  allocate(TOSAER(jsken,1:jcplx,JY,JX));TOSAER=0._r8

  end subroutine InitTflxType


!------------------------------------------------------------------------------------------

  subroutine DestructTflxType

  implicit none

  call destroy(trcx_TER)

  call destroy(TOMEERhetr)
  call destroy(TOMEERauto)
  call destroy(TOMQRS)
  call destroy(TORMER)
  call destroy(TOHMER)
  call destroy(TOSMER)
  call destroy(TOSAER)
  call destroy(trcSalt_Flo2MicP_vr)
  call destroy(trcs_Transp2MicP_vr)
  call destroy(trcSalt_Flo2MacP_vr)
  call destroy(TXGridSurfRunoff_2DH)
  call destroy(THeatXGridBySurfRunoff_2DH)
!  call destroy(TFLWS)
!  call destroy(TFLWW)
!  call destroy(TFLWI)
!  call destroy(THFLWW)

  call destroy(trcs_Transp2MacP_vr)

  call destroy(TCOQSS)
  call destroy(TCHQSS)
  call destroy(TOXQSS)
  call destroy(TNGQSS)
  call destroy(TN2QSS)
  call destroy(TN4QSS)
  call destroy(TN3QSS)
  call destroy(TNOQSS)
  call destroy(TPOQSS)
  call destroy(TP1QSS)

  call destroy(TSANER)
  call destroy(TSILER)
  call destroy(TCLAER)
  call destroy(TNH4Eros_col)
  call destroy(TNH3Eros_col)
  call destroy(TNUreaEros_col)
  call destroy(TNO3Eros_col)
  call destroy(TNH4ErosBand_col)
  call destroy(TNH3ErosBand_col)
  call destroy(TNUreaErosBand_col)
  call destroy(TNO3ErosBand_col)
  call destroy(tErosionSedmLoss)
  call destroy(TWatFlowCellMicP_vr)
  call destroy(TWatFlowCellMicPX_vr)
  call destroy(TWaterFlowMacP_vr)
  call destroy(WatIceThawMicP_vr)
  call destroy(THeatSoiThaw_vr)
  call destroy(WatIceThawMacP_vr)
  call destroy(VLWatMicP1_vr)
  call destroy(VLiceMicP1_vr)
  call destroy(VLWatMacP1_vr)
  call destroy(VLiceMacP1_vr)
  call destroy(DOM_Transp2Micp_vr)
  call destroy(DOM_Transp2Macp_flx)
  call destroy(trcp_TER)
  call destroy(Gas_AdvDif_Flx_vr)

  end subroutine DestructTflxType
end module TFlxTypeMod
