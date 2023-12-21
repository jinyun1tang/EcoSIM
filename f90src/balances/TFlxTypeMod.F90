module TFlxTypeMod

  use GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc,NFGs=>NFGsc
  use EcoSIMConfig, only : nlbiomcp=>NumOfLiveMicrobiomComponents,ndbiomcp=>NumOfDeadMicrobiomComponents
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
  real(r8),allocatable ::  TSEDER(:,:)                        !
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

  real(r8),allocatable :: TOMCER(:,:,:,:,:)
  real(r8),allocatable :: TOMNER(:,:,:,:,:)
  real(r8),allocatable :: TOMPER(:,:,:,:,:)


  real(r8),allocatable :: TOMCERff(:,:,:,:)
  real(r8),allocatable :: TOMNERff(:,:,:,:)
  real(r8),allocatable :: TOMPERff(:,:,:,:)

  real(r8),allocatable ::  DOM_Transp2Micp_flx(:,:,:,:,:)
  real(r8),allocatable ::  DOM_Transp2Macp_flx(:,:,:,:,:)
  real(r8),allocatable ::  TOCQRS(:,:,:)
  real(r8),allocatable ::  TONQRS(:,:,:)
  real(r8),allocatable ::  TOPQRS(:,:,:)
  real(r8),allocatable ::  TOAQRS(:,:,:)
  real(r8),allocatable ::  TORCER(:,:,:,:)
  real(r8),allocatable ::  TORNER(:,:,:,:)
  real(r8),allocatable ::  TORPER(:,:,:,:)
  real(r8),allocatable ::  TOHCER(:,:,:)
  real(r8),allocatable ::  TOHNER(:,:,:)
  real(r8),allocatable ::  TOHPER(:,:,:)
  real(r8),allocatable ::  TOHAER(:,:,:)
  real(r8),allocatable ::  TOSCER(:,:,:,:)
  real(r8),allocatable ::  TOSAER(:,:,:,:)
  real(r8),allocatable ::  TOSNER(:,:,:,:)
  real(r8),allocatable ::  TOSPER(:,:,:,:)

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

  allocate(TSEDER(JY,JX));      TSEDER=0._r8
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
  allocate(TOMCER(nlbiomcp,NumOfMicrobs1HetertrophCmplx,1:jcplx,JY,JX)); TOMCER=0._r8
  allocate(TOMNER(nlbiomcp,NumOfMicrobs1HetertrophCmplx,1:jcplx,JY,JX)); TOMNER=0._r8
  allocate(TOMPER(nlbiomcp,NumOfMicrobs1HetertrophCmplx,1:jcplx,JY,JX)); TOMPER=0._r8

  allocate(TOMCERff(nlbiomcp,NumOfMicrobsInAutotrophCmplx,JY,JX));TOMCERff=0._r8
  allocate(TOMNERff(nlbiomcp,NumOfMicrobsInAutotrophCmplx,JY,JX));TOMNERff=0._r8
  allocate(TOMPERff(nlbiomcp,NumOfMicrobsInAutotrophCmplx,JY,JX));TOMPERff=0._r8

  allocate(DOM_Transp2Micp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Micp_flx=0._r8
  allocate(DOM_Transp2Macp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flx=0._r8
  allocate(TOCQRS(1:jcplx,JY,JX));TOCQRS=0._r8
  allocate(TONQRS(1:jcplx,JY,JX));TONQRS=0._r8
  allocate(TOPQRS(1:jcplx,JY,JX));TOPQRS=0._r8
  allocate(TOAQRS(1:jcplx,JY,JX));TOAQRS=0._r8
  allocate(TORCER(ndbiomcp,1:jcplx,JY,JX));TORCER=0._r8
  allocate(TORNER(ndbiomcp,1:jcplx,JY,JX));TORNER=0._r8
  allocate(TORPER(ndbiomcp,1:jcplx,JY,JX));TORPER=0._r8
  allocate(TOHCER(1:jcplx,JY,JX));TOHCER=0._r8
  allocate(TOHNER(1:jcplx,JY,JX));TOHNER=0._r8
  allocate(TOHPER(1:jcplx,JY,JX));TOHPER=0._r8
  allocate(TOHAER(1:jcplx,JY,JX));TOHAER=0._r8
  allocate(TOSCER(jsken,1:jcplx,JY,JX));TOSCER=0._r8
  allocate(TOSAER(jsken,1:jcplx,JY,JX));TOSAER=0._r8
  allocate(TOSNER(jsken,1:jcplx,JY,JX));TOSNER=0._r8
  allocate(TOSPER(jsken,1:jcplx,JY,JX));TOSPER=0._r8

  end subroutine InitTflxType


!------------------------------------------------------------------------------------------

  subroutine DestructTflxType

  implicit none

  call destroy(trcx_TER)

  call destroy(TOMCER)
  call destroy(TOMNER)
  call destroy(TOMPER)
  call destroy(TOMCERff)
  call destroy(TOMNERff)
  call destroy(TOMPERff)
  call destroy(TOCQRS)
  call destroy(TONQRS)
  call destroy(TOPQRS)
  call destroy(TOAQRS)
  call destroy(TORCER)
  call destroy(TORNER)
  call destroy(TORPER)
  call destroy(TOHCER)
  call destroy(TOHNER)
  call destroy(TOHPER)
  call destroy(TOHAER)
  call destroy(TOSCER)
  call destroy(TOSAER)
  call destroy(TOSNER)
  call destroy(TOSPER)
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
  call destroy(TSEDER)
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
