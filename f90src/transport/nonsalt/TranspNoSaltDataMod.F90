module TranspNoSaltDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc
implicit none
  CHARACTER(LEN=*), Private,PARAMETER :: MOD_FILENAME=&
  __FILE__
  real(r8), parameter :: XFRS=0.05_r8
  real(r8), Parameter :: VFLWX=0.5

  real(r8), allocatable :: trcg_VFloSnow(:)
  real(r8), allocatable :: trcn_band_VFloSnow(:)
  real(r8), allocatable :: trcn_soil_VFloSnow(:)
  real(r8),allocatable :: DOM_Adv2MacP_flx(:,:)

  real(r8),allocatable :: CDOM_MacP1(:,:)
  real(r8),allocatable :: CDOM_MacP2(:,:)
  
  real(r8), allocatable :: DOM_Adv2MicP_flx(:,:)
  real(r8), allocatable :: Difus_Micp_flx_DOM(:,:)
  real(r8), allocatable :: Difus_Macp_flx_DOM(:,:)
  real(r8), allocatable :: CDOM_MicP1(:,:)
  real(r8), allocatable :: CDOM_MicP2(:,:)
  real(r8), allocatable :: DOM_MicP2(:,:,:,:,:)
  real(r8), allocatable :: DOM_3DMicp_Transp_flxM(:,:,:,:,:,:)
  real(r8), allocatable :: DOM_3DMacp_Transp_flxM(:,:,:,:,:,:)
  real(r8), allocatable :: DOM_Transp2Micp_flx(:,:,:,:,:)
  real(r8), allocatable :: RDOM_micb_cumflx(:,:,:,:,:)

  real(r8), allocatable :: dom_FloXSurRunoff(:,:,:,:)
  real(r8), allocatable :: RBGCSinkG_vr(:,:,:,:)   !BGC sink for gaseous tracers
  real(r8), allocatable :: RBGCSinkS(:,:,:,:)   !BGC sink for solute tracers
  real(r8), allocatable :: trcg_RBLS(:,:,:,:)                      !
  real(r8), allocatable :: trcn_RBLS(:,:,:,:)                      !
  real(r8), allocatable :: RDOMFL0(:,:,:,:)                      !
  real(r8), allocatable :: RDOMFL1(:,:,:,:)                      !
  real(r8), allocatable :: trcg_RFL0(:,:,:)                        !
  real(r8), allocatable :: trcn_RFL0(:,:,:)                        !
  real(r8), allocatable :: trcs_RFL1(:,:,:)                        !
  real(r8), allocatable :: DOM_XPoreTransp_flx(:,:,:,:,:)                    !

  real(r8), allocatable ::  POSGL2(:,:,:)                      !
  real(r8), allocatable ::  O2AquaDiffusvity2(:,:,:)                      !
  real(r8), allocatable ::  dom_TFloXSurRunoff(:,:,:,:)                       !

  real(r8), allocatable ::  trcg_SnowDrift(:,:,:)                        !
  real(r8), allocatable ::  trcn_SnowDrift(:,:,:)
  real(r8), allocatable ::  trcg_TFloXSurRunoff(:,:,:)
  real(r8), allocatable ::  trcn_TFloXSurRunoff(:,:,:)
  real(r8), allocatable ::  R3PorTSolFlx(:,:,:,:)          !total 3D micropore flux
  real(r8), allocatable ::  R3PorTSoHFlx(:,:,:,:)          !total 3D macropore flux

  real(r8), allocatable ::  GasDifc_vrc(:,:,:,:)
  real(r8), allocatable ::  SoluteDifusvty_vrc(:,:,:,:)
  real(r8), allocatable ::  DifuscG_vr(:,:,:,:,:)
  real(r8), allocatable ::  trcg_VLWatMicP(:,:,:,:)                      !

  real(r8), allocatable ::  DH2GG(:,:,:,:)                     !
  real(r8), allocatable ::  RHGFXS(:,:,:)                      !
  real(r8), allocatable ::  THGFLG(:,:,:)                      !
  real(r8), allocatable ::  CumReductVLsoiAirPM(:,:,:)                        !
  real(r8), allocatable ::  THETHL(:,:,:)                      !
  real(r8), allocatable ::  VLsoiAirPMA(:,:,:)                      !
  real(r8), allocatable ::  VLsoiAirPMB(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPMA(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPMB(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPXA(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPXB(:,:,:)                      !
  real(r8), allocatable ::  PARG_cef(:,:,:)                        !

  real(r8), allocatable ::  trcg_Ebu_vr(:,:,:,:)                      !
  real(r8), allocatable ::  trcg_FloXSurRunoff(:,:,:)                       !

  real(r8), allocatable ::  trcn_FloXSurRunoff(:,:,:)                       !
  real(r8), allocatable ::  DOM_MacP2(:,:,:,:,:)                     !
  real(r8), allocatable ::  DOM_Transp2Macp_flx(:,:,:,:,:)                    !
  real(r8), allocatable ::  DOMdiffusivity2_vr(:,:,:,:)                      !

  real(r8), allocatable :: trcg_solsml2(:,:,:,:)
  real(r8), allocatable :: trcn_solsml2(:,:,:,:)

  real(r8), allocatable ::  trcg_TBLS(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_TBLS(:,:,:,:)
  real(r8), allocatable ::  dom_2DFloXSurRunoffM(:,:,:,:,:,:)                   !
  real(r8), allocatable ::  trcg_2DFloXSurRunoffM(:,:,:,:,:)                    !

  real(r8), allocatable ::  RDFR_gas(:,:,:)                        !

  real(r8), allocatable :: trc_gasml_vr2(:,:,:,:)
  real(r8), allocatable :: trc_solml_vr2(:,:,:,:)
  real(r8), allocatable :: trc_soHml2_vr(:,:,:,:)

  real(r8), allocatable ::  trcn_2DSnowDrift(:,:,:,:)                      !
  real(r8), allocatable ::  trcg_2DSnowDrift(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_2DFloXSurRunoffM(:,:,:,:,:)                    !
                  !
  real(r8), allocatable ::  RGasSSVol(:,:,:)     !soil surface gas volatization
  real(r8), allocatable ::  RGasDSFlx(:,:,:,:)   !gas dissolution-volatilization
  real(r8), allocatable ::  R3GasADFlx(:,:,:,:,:) !3D gas flux advection + diffusion
  real(r8), allocatable ::  Gas_AdvDif_Flx_vr(:,:,:,:)  !total 3D gas flux advection + diffusion

  real(r8), allocatable :: R3PoreSoHFlx(:,:,:,:,:)        !3D macropore flux
  real(r8), allocatable :: R3PoreSolFlx_vr(:,:,:,:,:)        !3D micropore flux
  real(r8), allocatable :: RporeSoXFlx(:,:,:,:)        !Mac-mic pore exchange flux
!----------------------------------------------------------------------

contains
  subroutine InitTransfrData
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer :: NumOfLitrCmplxs
  NumOfLitrCmplxs=micpar%NumOfLitrCmplxs
  allocate(CDOM_MacP1(idom_beg:idom_end,1:jcplx));CDOM_MacP1=0._r8
  allocate(CDOM_MacP2(idom_beg:idom_end,1:jcplx));CDOM_MacP2=0._r8

  allocate(dom_FloXSurRunoff(idom_beg:idom_end,1:jcplx,JV,JH));    dom_FloXSurRunoff=0._r8
  allocate(DOM_MicP2(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));    DOM_MicP2=0._r8
  allocate(RDOM_micb_cumflx(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));  RDOM_micb_cumflx=0._r8
  allocate(DOM_3DMicp_Transp_flxM(idom_beg:idom_end,1:jcplx,3,0:JD,JV,JH));DOM_3DMicp_Transp_flxM=0._r8
  allocate(DOM_3DMacp_Transp_flxM(idom_beg:idom_end,1:jcplx,3,JD,JV,JH));  DOM_3DMacp_Transp_flxM=0._r8
  allocate(DOM_Transp2Micp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));    DOM_Transp2Micp_flx=0._r8

  allocate(trcg_RBLS(idg_beg:idg_NH3,JS,JY,JX));   trcg_RBLS=0._r8
  allocate(trcn_RBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX));   trcn_RBLS=0._r8
  allocate(RDOMFL0(idom_beg:idom_end,1:NumOfLitrCmplxs,JY,JX));  RDOMFL0=0._r8
  allocate(RDOMFL1(idom_beg:idom_end,1:NumOfLitrCmplxs,JY,JX));  RDOMFL1=0._r8
  allocate(trcg_RFL0(idg_beg:idg_NH3,JY,JX));      trcg_RFL0=0._r8
  allocate(trcn_RFL0(ids_nut_beg:ids_nuts_end,JY,JX));      trcn_RFL0=0._r8
  allocate(trcs_RFL1(ids_beg:ids_end,JY,JX));      trcs_RFL1=0._r8
  allocate(DOM_XPoreTransp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_XPoreTransp_flx=0._r8

  allocate(POSGL2(0:JZ,JY,JX)); POSGL2=0._r8
  allocate(O2AquaDiffusvity2(0:JZ,JY,JX)); O2AquaDiffusvity2=0._r8

  allocate(dom_TFloXSurRunoff(idom_beg:idom_end,1:jcplx,JY,JX)); dom_TFloXSurRunoff=0._r8
  allocate(Difus_Macp_flx_DOM(idom_beg:idom_end,1:jcplx));Difus_Macp_flx_DOM=0._r8
  allocate(trcg_SnowDrift(idg_beg:idg_NH3,JY,JX));      trcg_SnowDrift=0._r8
  allocate(trcn_SnowDrift(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_SnowDrift=0._r8
  allocate(trcg_TFloXSurRunoff(idg_beg:idg_NH3,JY,JX)); trcg_TFloXSurRunoff=0._r8
  allocate(trcn_TFloXSurRunoff(ids_nut_beg:ids_nuts_end,JY,JX));trcn_TFloXSurRunoff=0._r8

  allocate(CDOM_MicP1(idom_beg:idom_end,1:jcplx)); CDOM_MicP1=0._r8
  allocate(CDOM_MicP2(idom_beg:idom_end,1:jcplx)); CDOM_MicP2=0._r8

  allocate(GasDifc_vrc(idg_beg:idg_end,JZ,JY,JX))
  allocate(SoluteDifusvty_vrc(ids_beg:ids_end,0:JZ,JY,JX));SoluteDifusvty_vrc=0._r8
  allocate(DifuscG_vr(idg_beg:idg_end,3,JZ,JY,JX)); DifuscG_vr=0._r8

  allocate(trcg_VLWatMicP(idg_beg:idg_end,0:JZ,JY,JX)); trcg_VLWatMicP=0._r8

  allocate(DH2GG(3,JZ,JY,JX));  DH2GG=0._r8
  allocate(RHGFXS(JZ,JY,JX));   RHGFXS=0._r8
  allocate(THGFLG(JZ,JY,JX));   THGFLG=0._r8
  allocate(CumReductVLsoiAirPM(JZ,JY,JX));     CumReductVLsoiAirPM=0._r8
  allocate(THETHL(JZ,JY,JX));   THETHL=0._r8
  allocate(VLsoiAirPMA(JZ,JY,JX));   VLsoiAirPMA=0._r8
  allocate(VLsoiAirPMB(JZ,JY,JX));   VLsoiAirPMB=0._r8
  allocate(VLWatMicPMA(JZ,JY,JX));   VLWatMicPMA=0._r8
  allocate(VLWatMicPMB(JZ,JY,JX));   VLWatMicPMB=0._r8
  allocate(VLWatMicPXA(0:JZ,JY,JX)); VLWatMicPXA=0._r8
  allocate(VLWatMicPXB(JZ,JY,JX));   VLWatMicPXB=0._r8
  allocate(PARG_cef(idg_beg:idg_NH3,JY,JX));      PARG_cef=0._r8

  allocate(RBGCSinkG_vr(idg_beg:idg_end,0:JZ,JY,JX));RBGCSinkG_vr=0._r8
  allocate(RBGCSinkS(ids_nuts_beg:ids_nuts_end,0:JZ,JY,JX));RBGCSinkS=0._r8

  allocate(trcg_Ebu_vr(idg_beg:idg_end,JZ,JY,JX));   trcg_Ebu_vr=0._r8
  allocate(trcg_FloXSurRunoff(idg_beg:idg_end,JV,JH));     trcg_FloXSurRunoff=0._r8
  allocate(trcn_FloXSurRunoff(ids_nut_beg:ids_nuts_end,JV,JH));     trcn_FloXSurRunoff=0._r8
  allocate(DOM_MacP2(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_MacP2=0._r8
  allocate(DOM_Transp2Macp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flx=0._r8

  allocate(DOMdiffusivity2_vr(idom_beg:idom_end,0:JZ,JY,JX)); DOMdiffusivity2_vr=0._r8

  allocate(trcg_solsml2(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml2=0._r8
  allocate(trcn_solsml2(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml2=0._r8

  allocate(trcg_TBLS(idg_beg:idg_NH3,JS,JY,JX));   trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_TBLS=0._r8

  allocate(dom_2DFloXSurRunoffM(idom_beg:idom_end,1:jcplx,2,2,JV,JH));dom_2DFloXSurRunoffM=0._r8
  allocate(trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,2,2,JV,JH));  trcg_2DFloXSurRunoffM=0._r8

  allocate(RDFR_gas(idg_beg:idg_NH3,JY,JX));      RDFR_gas=0._r8

  allocate(trc_gasml_vr2(idg_beg:idg_end,0:JZ,JY,JX)); trc_gasml_vr2(:,:,:,:) = 0._r8
  allocate(trc_solml_vr2(ids_beg:ids_end,0:JZ,JY,JX)); trc_solml_vr2(:,:,:,:) = 0._r8
  allocate(trc_soHml2_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_soHml2_vr(:,:,:,:)       = 0._r8

  allocate(trcg_2DSnowDrift(idg_beg:idg_NH3,2,JV,JH));    trcg_2DSnowDrift                    = 0._r8
  allocate(trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_2DFloXSurRunoffM = 0._r8
  allocate(R3GasADFlx(idg_beg:idg_NH3,3,JD,JV,JH));R3GasADFlx                                 = 0._r8
  allocate(Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,JZ,JY,JH));Gas_AdvDif_Flx_vr                     = 0._r8

  allocate(DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx));DOM_Adv2MicP_flx = 0._r8

  allocate(RGasDSFlx(idg_beg:idg_end,0:JZ,JY,JX)); RGasDSFlx = 0._r8
  allocate(RGasSSVol(idg_beg:idg_end,JY,JX)); RGasSSVol      = 0._r8

  allocate(trcn_2DSnowDrift(ids_nut_beg:ids_nuts_end,2,JV,JH)); trcn_2DSnowDrift = 0._r8


  allocate(R3PoreSoHFlx(ids_beg:ids_end,3,JD,JV,JH));R3PoreSoHFlx            = 0._r8
  allocate(R3PoreSolFlx_vr(ids_beg:ids_end,3,0:JD,JV,JH));R3PoreSolFlx_vr    = 0._r8
  allocate(RporeSoXFlx(ids_beg:ids_end,JZ,JY,JX));RporeSoXFlx                = 0._r8
  allocate(R3PorTSolFlx(ids_beg:ids_end,JZ,JY,JX));R3PorTSolFlx              = 0._r8
  allocate(R3PorTSoHFlx(ids_beg:ids_end,JZ,JY,JX));R3PorTSoHFlx              = 0._r8
  allocate(trcg_VFloSnow(idg_beg:idg_end));trcg_VFloSnow                     = 0._r8
  allocate(trcn_band_VFloSnow(ids_nutb_beg:ids_nutb_end));trcn_band_VFloSnow = 0._r8
  allocate(trcn_soil_VFloSnow(ids_nut_beg:ids_nuts_end));trcn_soil_VFloSnow  = 0._r8
  allocate(DOM_Adv2MacP_flx(idom_beg:idom_end,1:jcplx));
  allocate(Difus_Micp_flx_DOM(idom_beg:idom_end,1:jcplx));

  end subroutine InitTransfrData

!----------------------------------------------------------------------
  subroutine DestructTransfrData
  use abortutils, only : destroy
  implicit none

  call destroy(DOM_MicP2)
  call destroy(CDOM_MicP1)
  call destroy(CDOM_MicP2)
  call destroy(RDOM_micb_cumflx)
  call destroy(DOM_3DMicp_Transp_flxM)
  call destroy(DOM_3DMacp_Transp_flxM)

  call destroy(dom_FloXSurRunoff)
  call destroy(dom_2DFloXSurRunoffM)
  call destroy(RDOMFL0)
  call destroy(RDOMFL1)

  call destroy(trcs_RFL1)

  call destroy(DOM_XPoreTransp_flx)

  call destroy(trcg_Ebu_vr)

  call destroy(POSGL2)
  call destroy(O2AquaDiffusvity2)

  call destroy(dom_TFloXSurRunoff)
  call destroy(trcg_TFloXSurRunoff)
  call destroy(trcn_TFloXSurRunoff)
  call destroy(trcg_SnowDrift)
  call destroy(trcn_SnowDrift)

  call destroy(trcg_VLWatMicP)

  call destroy(DH2GG)
  call destroy(RHGFXS)
  call destroy(THGFLG)
  call destroy(CumReductVLsoiAirPM)
  call destroy(THETHL)
  call destroy(VLsoiAirPMA)
  call destroy(VLsoiAirPMB)
  call destroy(VLWatMicPMA)
  call destroy(VLWatMicPMB)
  call destroy(VLWatMicPXA)
  call destroy(VLWatMicPXB)
  call destroy(PARG_cef)

  call destroy(DOMdiffusivity2_vr)

  call destroy(trcg_solsml2)
  call destroy(trcn_solsml2)

  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)
  call destroy(trcg_RFL0)
  call destroy(trcg_2DSnowDrift)
  call destroy(Difus_Micp_flx_DOM)
  call destroy(Difus_Macp_flx_DOM)


  call destroy(RDFR_gas)
  call destroy(DOM_MacP2)

  call destroy(trc_gasml_vr2)
  call destroy(trc_solml_vr2)
  call destroy(trc_soHml2_vr)
  call destroy(trcn_band_VFloSnow)
  call destroy(trcn_soil_VFloSnow)

  call destroy(RGasDSFlx)
  call destroy(trcn_2DSnowDrift)

  call destroy(R3PoreSoHFlx)
  call destroy(R3PoreSolFlx_vr)
  call destroy(RporeSoXFlx)
  call destroy(CDOM_MacP1)
  call destroy(CDOM_MacP2)
  call destroy(DOM_Adv2MicP_flx)
  call destroy(DOM_Adv2MacP_flx)
  call destroy(DOM_Transp2Macp_flx)
  call destroy(DOM_Transp2Micp_flx)
  end subroutine DestructTransfrData

end module TranspNoSaltDataMod
