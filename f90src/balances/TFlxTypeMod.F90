module TFlxTypeMod

  use GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc,NumMicbFunGroups=>NumMicbFunGroups
  use EcoSIMConfig, only : nlbiomcp=>NumLiveMicrbCompts,ndbiomcp=>NumDeadMicrbCompts
implicit none

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public

  real(r8),allocatable ::  trcs_Transp2MicP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcs_Transp2MacP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  TWat2GridBySurfRunoff(:,:)                           !
  real(r8),allocatable ::  THeat2GridBySurfRunoff(:,:)                          !

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
  real(r8),allocatable ::  TNH4ER(:,:)                        !
  real(r8),allocatable ::  TNH3ER(:,:)                        !
  real(r8),allocatable ::  TNHUER(:,:)                        !
  real(r8),allocatable ::  TNO3ER(:,:)                        !
  real(r8),allocatable ::  TNH4EB(:,:)                        !
  real(r8),allocatable ::  TNH3EB(:,:)                        !
  real(r8),allocatable ::  TNHUEB(:,:)                        !
  real(r8),allocatable ::  TNO3EB(:,:)                        !
  real(r8),allocatable ::  trcx_TER(:,:,:)                         !
  real(r8),allocatable ::  trcp_TER(:,:,:)                        !
  real(r8),allocatable ::  tErosionSedmLoss(:,:)                        !
  real(r8),allocatable ::  TWatFlowCellMicP(:,:,:)                        !
  real(r8),allocatable ::  TWatFlowCellMicPX(:,:,:)                       !
  real(r8),allocatable ::  THeatFlow2Soil(:,:,:)                       !
  real(r8),allocatable ::  TWaterFlowMacP(:,:,:)                       !

  real(r8),allocatable ::  Gas_AdvDif_Flx_vr(:,:,:,:)                      !

  real(r8),allocatable ::  WatIceThawMicP(:,:,:)                       !
  real(r8),allocatable ::  THeatSoiThaw(:,:,:)                      !
  real(r8),allocatable ::  WatIceThawMacP(:,:,:)                      !
  real(r8),allocatable ::  VLWatMicP1(:,:,:)                       !
  real(r8),allocatable ::  VLiceMicP1(:,:,:)                       !
  real(r8),allocatable ::  VLWatMacP1(:,:,:)                      !
  real(r8),allocatable ::  VLiceMacP1(:,:,:)                      !

  real(r8),allocatable :: TOMEERhetr(:,:,:,:,:)

  real(r8),allocatable :: TOMEERauto(:,:,:,:)
 
  real(r8),allocatable ::  DOM_Transp2Micp_flx(:,:,:,:,:)
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

  allocate(TWat2GridBySurfRunoff(JY,JX));         TWat2GridBySurfRunoff=0._r8
  allocate(THeat2GridBySurfRunoff(JY,JX));        THeat2GridBySurfRunoff=0._r8
  allocate(TWat2GridBySurfRunoff(JY,JX));         TWat2GridBySurfRunoff=0._r8
  allocate(THeat2GridBySurfRunoff(JY,JX));        THeat2GridBySurfRunoff=0._r8

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
  allocate(TNH4ER(JY,JX));      TNH4ER=0._r8
  allocate(TNH3ER(JY,JX));      TNH3ER=0._r8
  allocate(TNHUER(JY,JX));      TNHUER=0._r8
  allocate(TNO3ER(JY,JX));      TNO3ER=0._r8
  allocate(TNH4EB(JY,JX));      TNH4EB=0._r8
  allocate(TNH3EB(JY,JX));      TNH3EB=0._r8
  allocate(TNHUEB(JY,JX));      TNHUEB=0._r8
  allocate(TNO3EB(JY,JX));      TNO3EB=0._r8

  allocate(trcx_TER(idx_beg:idx_end,JY,JX));    trcx_TER=0._r8
  allocate(trcp_TER(idsp_beg:idsp_end,JY,JX));      trcp_TER=0._r8
  allocate(trcs_Transp2MicP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_Transp2MicP_vr=0._r8

  allocate(tErosionSedmLoss(JY,JX));      tErosionSedmLoss=0._r8
  allocate(TWatFlowCellMicP(JZ,JY,JX));     TWatFlowCellMicP=0._r8
  allocate(TWatFlowCellMicPX(JZ,JY,JX));    TWatFlowCellMicPX=0._r8
  allocate(THeatFlow2Soil(JZ,JY,JX));    THeatFlow2Soil=0._r8
  allocate(TWaterFlowMacP(JZ,JY,JX));    TWaterFlowMacP=0._r8
  allocate(WatIceThawMicP(JZ,JY,JX));    WatIceThawMicP=0._r8
  allocate(THeatSoiThaw(JZ,JY,JX));   THeatSoiThaw=0._r8
  allocate(WatIceThawMacP(JZ,JY,JX));   WatIceThawMacP=0._r8
  allocate(VLWatMicP1(JZ,JY,JX));    VLWatMicP1=0._r8
  allocate(VLiceMicP1(JZ,JY,JX));    VLiceMicP1=0._r8
  allocate(VLWatMacP1(JZ,JY,JX));   VLWatMacP1=0._r8
  allocate(VLiceMacP1(JZ,JY,JX));   VLiceMacP1=0._r8
  allocate(TOMEERhetr(NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,JY,JX)); TOMEERhetr=0._r8

  allocate(TOMEERauto(NumPlantChemElms,1:NumLiveAutoBioms,JY,JX));TOMEERauto=0._r8

  allocate(DOM_Transp2Micp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Micp_flx=0._r8
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
  call destroy(TWat2GridBySurfRunoff)
  call destroy(THeat2GridBySurfRunoff)
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
  call destroy(TNH4ER)
  call destroy(TNH3ER)
  call destroy(TNHUER)
  call destroy(TNO3ER)
  call destroy(TNH4EB)
  call destroy(TNH3EB)
  call destroy(TNHUEB)
  call destroy(TNO3EB)
  call destroy(tErosionSedmLoss)
  call destroy(TWatFlowCellMicP)
  call destroy(TWatFlowCellMicPX)
  call destroy(THeatFlow2Soil)
  call destroy(TWaterFlowMacP)
  call destroy(WatIceThawMicP)
  call destroy(THeatSoiThaw)
  call destroy(WatIceThawMacP)
  call destroy(VLWatMicP1)
  call destroy(VLiceMicP1)
  call destroy(VLWatMacP1)
  call destroy(VLiceMacP1)
  call destroy(DOM_Transp2Micp_flx)
  call destroy(DOM_Transp2Macp_flx)
  call destroy(trcp_TER)
  call destroy(Gas_AdvDif_Flx_vr)

  end subroutine DestructTflxType
end module TFlxTypeMod
