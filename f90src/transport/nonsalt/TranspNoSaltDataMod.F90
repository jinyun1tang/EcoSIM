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
  real(r8), allocatable :: RBGCSinkG(:,:,:,:)   !BGC sink for gaseous tracers
  real(r8), allocatable :: RBGCSinkS(:,:,:,:)   !BGC sink for solute tracers
  real(r8), allocatable ::  trcg_RBLS(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_RBLS(:,:,:,:)                      !
  real(r8), allocatable ::  ROCFL0(:,:,:)                      !
  real(r8), allocatable ::  RONFL0(:,:,:)                      !
  real(r8), allocatable ::  ROPFL0(:,:,:)                      !
  real(r8), allocatable ::  ROAFL0(:,:,:)                      !
  real(r8), allocatable ::  ROCFL1(:,:,:)                      !
  real(r8), allocatable ::  RONFL1(:,:,:)                      !
  real(r8), allocatable ::  ROPFL1(:,:,:)                      !
  real(r8), allocatable ::  ROAFL1(:,:,:)                      !
  real(r8), allocatable ::  trcg_RFL0(:,:,:)                        !
  real(r8), allocatable ::  trcn_RFL0(:,:,:)                        !
  real(r8), allocatable ::  trcs_RFL1(:,:,:)                        !
  real(r8), allocatable ::  DOM_XPoreTransp_flx(:,:,:,:,:)                    !
  real(r8), allocatable ::  RCOFXS(:,:,:)                      !
  real(r8), allocatable ::  RCHFXS(:,:,:)                      !
  real(r8), allocatable ::  ROXFXS(:,:,:)                      !
  real(r8), allocatable ::  RNGFXS(:,:,:)                      !
  real(r8), allocatable ::  RN2FXS(:,:,:)                      !
  real(r8), allocatable ::  RN4FXW(:,:,:)                      !
  real(r8), allocatable ::  RN3FXW(:,:,:)                      !
  real(r8), allocatable ::  RNOFXW(:,:,:)                      !
  real(r8), allocatable ::  RH2PXS(:,:,:)                      !
  real(r8), allocatable ::  RN4FXB(:,:,:)                      !
  real(r8), allocatable ::  RN3FXB(:,:,:)                      !
  real(r8), allocatable ::  RNOFXB(:,:,:)                      !
  real(r8), allocatable ::  RH2BXB(:,:,:)                      !
  real(r8), allocatable ::  RNXFXS(:,:,:)                      !
  real(r8), allocatable ::  RNXFXB(:,:,:)                      !
  real(r8), allocatable ::  RH1PXS(:,:,:)                      !
  real(r8), allocatable ::  RH1BXB(:,:,:)                      !
  real(r8), allocatable ::  CLSGL2(:,:,:)                      !
  real(r8), allocatable ::  CQSGL2(:,:,:)                      !
  real(r8), allocatable ::  POSGL2(:,:,:)                      !
  real(r8), allocatable ::  OLSGL2(:,:,:)                      !
  real(r8), allocatable ::  ZNSGL2(:,:,:)                      !
  real(r8), allocatable ::  ZLSGL2(:,:,:)                      !
  real(r8), allocatable ::  ZVSGL2(:,:,:)                      !
  real(r8), allocatable ::  HLSGL2(:,:,:)                      !
  real(r8), allocatable ::  ZOSGL2(:,:,:)                      !
  real(r8), allocatable ::  TQRCOS(:,:)
  real(r8), allocatable ::  TQSCOS(:,:)
  real(r8), allocatable ::  RNBDFG(:,:,:)                      !NH3 volatilization/dissolution within soil layer band	g d-2 t-1
  real(r8), allocatable ::  dom_TFloXSurRunoff(:,:,:,:)                       !

  real(r8), allocatable ::  trcg_SnowDrift(:,:,:)                        !
  real(r8), allocatable ::  trcn_SnowDrift(:,:,:)
  real(r8), allocatable ::  trcg_TFloXSurRunoff(:,:,:)
  real(r8), allocatable ::  trcn_TFloXSurRunoff(:,:,:)
  real(r8), allocatable ::  R3PorTSolFlx(:,:,:,:) !total 3D micropore flux
  real(r8), allocatable ::  R3PorTSoHFlx(:,:,:,:) !total 3D macropore flux

  real(r8), allocatable ::  TN4FLW(:,:,:)                      !
  real(r8), allocatable ::  TN3FLW(:,:,:)                      !
  real(r8), allocatable ::  TNOFLW(:,:,:)                      !
  real(r8), allocatable ::  TH2PFS(:,:,:)                      !
  real(r8), allocatable ::  TQRH1P(:,:)                        !
  real(r8), allocatable ::  TH1PFS(:,:,:)                      !

  real(r8), allocatable ::  GasDifc_vrc(:,:,:,:)
  real(r8), allocatable ::  SolDifc_vrc(:,:,:,:)
  real(r8), allocatable ::  DifuscG(:,:,:,:,:)
  real(r8), allocatable ::  DCO2G(:,:,:,:)                     !
  real(r8), allocatable ::  DCH4G(:,:,:,:)                     !
  real(r8), allocatable ::  DOXYG(:,:,:,:)                     !
  real(r8), allocatable ::  DZ2GG(:,:,:,:)                     !
  real(r8), allocatable ::  DZ2OG(:,:,:,:)                     !
  real(r8), allocatable ::  DNH3G(:,:,:,:)                     !
  real(r8), allocatable ::  trcg_VLWatMicP(:,:,:,:)                      !
  real(r8), allocatable ::  HGSGL2(:,:,:)                      !
  real(r8), allocatable ::  DH2GG(:,:,:,:)                     !
  real(r8), allocatable ::  RHGFXS(:,:,:)                      !
  real(r8), allocatable ::  THGFLG(:,:,:)                      !
  real(r8), allocatable ::  CumReductVLsoiAirPM(:,:,:)                        !
  real(r8), allocatable ::  THETH2(:,:,:)                      !
  real(r8), allocatable ::  THETHL(:,:,:)                      !
  real(r8), allocatable ::  VLsoiAirPMA(:,:,:)                      !
  real(r8), allocatable ::  VLsoiAirPMB(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPMA(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPMB(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPXA(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPXB(:,:,:)                      !
  real(r8), allocatable ::  PARG_cef(:,:,:)                        !

  real(r8), allocatable ::  RCOSK2(:,:,:)                      !
  real(r8), allocatable ::  ROXSK2(:,:,:)                      !
  real(r8), allocatable ::  RCHSK2(:,:,:)                      !
  real(r8), allocatable ::  RNGSK2(:,:,:)                      !
  real(r8), allocatable ::  RN2SK2(:,:,:)                      !
  real(r8), allocatable ::  RN4SK2(:,:,:)                      !
  real(r8), allocatable ::  RN3SK2(:,:,:)                      !
  real(r8), allocatable ::  RNOSK2(:,:,:)                      !
  real(r8), allocatable ::  RHPSK2(:,:,:)                      !
  real(r8), allocatable ::  R4BSK2(:,:,:)                      !
  real(r8), allocatable ::  R3BSK2(:,:,:)                      !
  real(r8), allocatable ::  RNBSK2(:,:,:)                      !
  real(r8), allocatable ::  RHBSK2(:,:,:)                      !
  real(r8), allocatable ::  RNXSK2(:,:,:)                      !
  real(r8), allocatable ::  RNZSK2(:,:,:)                      !
  real(r8), allocatable ::  RHGSK2(:,:,:)                      !
  real(r8), allocatable ::  RNHSK2(:,:,:)                      !
  real(r8), allocatable ::  R1PSK2(:,:,:)                      !

  real(r8), allocatable ::  TN3FLG(:,:,:)                      !
  real(r8), allocatable ::  trcg_BBL(:,:,:,:)                      !
  real(r8), allocatable ::  trcg_FloXSurRunoff(:,:,:)                       !

  real(r8), allocatable ::  trcn_FloXSurRunoff(:,:,:)                       !
  real(r8), allocatable ::  DOM_MacP2(:,:,:,:,:)                     !
  real(r8), allocatable ::  DOM_Transp2Macp_flx(:,:,:,:,:)                    !
  real(r8), allocatable ::  TN4FHW(:,:,:)                      !
  real(r8), allocatable ::  TN3FHW(:,:,:)                      !
  real(r8), allocatable ::  TNOFHW(:,:,:)                      !
  real(r8), allocatable ::  TH2PHS(:,:,:)                      !
  real(r8), allocatable ::  ZNO2B2(:,:,:)                      !
  real(r8), allocatable ::  ZN2BH2(:,:,:)                      !
  real(r8), allocatable ::  H1P4H2(:,:,:)                      !
  real(r8), allocatable ::  H1PBH2(:,:,:)                      !
  real(r8), allocatable ::  TH1PHS(:,:,:)                      !
  REAL(R8), allocatable ::  H1PO42(:,:,:)                      !
  real(r8), allocatable ::  H1POB2(:,:,:)                      !
  real(r8), allocatable ::  OCSGL2(:,:,:)                      !
  real(r8), allocatable ::  ONSGL2(:,:,:)                      !
  real(r8), allocatable ::  OPSGL2(:,:,:)                      !
  real(r8), allocatable ::  OASGL2(:,:,:)                      !

  real(r8), allocatable :: trcg_solsml2(:,:,:,:)
  real(r8), allocatable :: trcn_solsml2(:,:,:,:)
  real(r8), allocatable ::  CO2W2(:,:,:)                       !
  real(r8), allocatable ::  CH4W2(:,:,:)                       !
  real(r8), allocatable ::  OXYW2(:,:,:)                       !
  real(r8), allocatable ::  ZNGW2(:,:,:)                       !
  real(r8), allocatable ::  ZN2W2(:,:,:)                       !
  real(r8), allocatable ::  ZN4W2(:,:,:)                       !
  real(r8), allocatable ::  ZN3W2(:,:,:)                       !
  real(r8), allocatable ::  ZNOW2(:,:,:)                       !
  real(r8), allocatable ::  ZHPW2(:,:,:)                       !
  real(r8), allocatable ::  Z1PW2(:,:,:)                       !
  real(r8), allocatable ::  trcg_TBLS(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_TBLS(:,:,:,:)
  real(r8), allocatable ::  dom_2DFloXSurRunoffM(:,:,:,:,:,:)                   !
  real(r8), allocatable ::  trcg_2DFloXSurRunoffM(:,:,:,:,:)                    !
  real(r8), allocatable ::  RQRHGS(:,:,:,:)                    !
  real(r8), allocatable ::  RQRNH4(:,:,:,:)                    !
  real(r8), allocatable ::  RCODFS(:,:)                        !
  real(r8), allocatable ::  RCHDFS(:,:)                        !
  real(r8), allocatable ::  ROXDFS(:,:)                        !
  real(r8), allocatable ::  RNGDFS(:,:)                        !
  real(r8), allocatable ::  RN2DFS(:,:)                        !
  real(r8), allocatable ::  RN3DFS(:,:)                        !
  real(r8), allocatable ::  RNBDFS(:,:)                        !
  real(r8), allocatable ::  RHGDFS(:,:)                        !
  real(r8), allocatable ::  RDFR_gas(:,:,:)                        !
  real(r8), allocatable ::  R1BSK2(:,:,:)                      !

  real(r8), allocatable ::  ZNH4H2(:,:,:)                      !
  real(r8), allocatable ::  ZN4BH2(:,:,:)                      !
  real(r8), allocatable ::  ZNH3H2(:,:,:)                      !
  real(r8), allocatable ::  ZN3BH2(:,:,:)                      !
  real(r8), allocatable ::  ZNO3H2(:,:,:)                      !
  real(r8), allocatable ::  ZNOBH2(:,:,:)                      !
  real(r8), allocatable ::  H2P4H2(:,:,:)                      !
  real(r8), allocatable ::  H2PBH2(:,:,:)                      !
  real(r8), allocatable ::  ZNO2H2(:,:,:)                      !

  real(r8), allocatable :: trc_gasml_vr2(:,:,:,:)
  real(r8), allocatable :: trc_solml_vr2(:,:,:,:)
  real(r8), allocatable :: trc_soHml2(:,:,:,:)

  real(r8), allocatable ::  trcn_2DSnowDrift(:,:,:,:)                      !
  real(r8), allocatable ::  trcg_2DSnowDrift(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_2DFloXSurRunoffM(:,:,:,:,:)                    !

  real(r8), allocatable ::  RNXFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RNXFHB(:,:,:,:)                    !
  real(r8), allocatable ::  RH1PFS(:,:,:,:)                    !
  real(r8), allocatable ::  RH1BFB(:,:,:,:)                    !
  real(r8), allocatable ::  RH1PHS(:,:,:,:)                    !
  real(r8), allocatable ::  RH1BHB(:,:,:,:)                    !
  real(r8), allocatable ::  RN3FHB(:,:,:,:)                    !
  real(r8), allocatable ::  RNOFHB(:,:,:,:)                    !
  real(r8), allocatable ::  RH2BHB(:,:,:,:)                    !
  real(r8), allocatable ::  RNOFHW(:,:,:,:)                    !
  real(r8), allocatable ::  RH2PHS(:,:,:,:)                    !
  real(r8), allocatable ::  RN4FHB(:,:,:,:)                    !
  real(r8), allocatable ::  RN2FHS(:,:,:,:)                    !
  real(r8), allocatable ::  RN4FHW(:,:,:,:)                    !
  real(r8), allocatable ::  RN3FHW(:,:,:,:)
                  !
  real(r8), allocatable ::  RGasSSVol(:,:,:)     !soil surface gas volatization
  real(r8), allocatable ::  RGasDSFlx(:,:,:,:)   !gas dissolution-volatilization
  real(r8), allocatable ::  R3GasADFlx(:,:,:,:,:) !3D gas flux advection + diffusion
  real(r8), allocatable ::  Gas_AdvDif_Flx_vr(:,:,:,:)  !total 3D gas flux advection + diffusion

  real(r8), allocatable ::  RHGFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RCHFHS(:,:,:,:)                    !
  real(r8), allocatable ::  ROXFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RNGFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RCOFLS(:,:,:,:)                    !
  real(r8), allocatable ::  RCHFLS(:,:,:,:)                    !
  real(r8), allocatable ::  ROXFLS(:,:,:,:)                    !
  real(r8), allocatable ::  RNGFLS(:,:,:,:)                    !
  real(r8), allocatable ::  RN2FLS(:,:,:,:)                    !
  real(r8), allocatable ::  RHGFLS(:,:,:,:)                    !
  real(r8), allocatable ::  RN4FLW(:,:,:,:)                    !
  real(r8), allocatable ::  RN3FLW(:,:,:,:)                    !
  real(r8), allocatable ::  RNOFLW(:,:,:,:)                    !
  real(r8), allocatable ::  RNXFLS(:,:,:,:)                    !
  real(r8), allocatable ::  RH2PFS(:,:,:,:)                    !
  real(r8), allocatable ::  RN4FLB(:,:,:,:)                    !
  real(r8), allocatable ::  RN3FLB(:,:,:,:)                    !
  real(r8), allocatable ::  RNOFLB(:,:,:,:)                    !
  real(r8), allocatable ::  RNXFLB(:,:,:,:)                    !
  real(r8), allocatable ::  RH2BFB(:,:,:,:)                    !
  real(r8), allocatable ::  RCOFHS(:,:,:,:)                    !

  real(r8), allocatable :: R3PoreSoHFlx(:,:,:,:,:)        !3D macropore flux
  real(r8), allocatable :: R3PoreSolFlx(:,:,:,:,:)        !3D micropore flux
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
  allocate(ROCFL0(1:NumOfLitrCmplxs,JY,JX));  ROCFL0=0._r8
  allocate(RONFL0(1:NumOfLitrCmplxs,JY,JX));  RONFL0=0._r8
  allocate(ROPFL0(1:NumOfLitrCmplxs,JY,JX));  ROPFL0=0._r8
  allocate(ROAFL0(1:NumOfLitrCmplxs,JY,JX));  ROAFL0=0._r8
  allocate(ROCFL1(1:NumOfLitrCmplxs,JY,JX));  ROCFL1=0._r8
  allocate(RONFL1(1:NumOfLitrCmplxs,JY,JX));  RONFL1=0._r8
  allocate(ROPFL1(1:NumOfLitrCmplxs,JY,JX));  ROPFL1=0._r8
  allocate(ROAFL1(1:NumOfLitrCmplxs,JY,JX));  ROAFL1=0._r8
  allocate(trcg_RFL0(idg_beg:idg_NH3,JY,JX));      trcg_RFL0=0._r8
  allocate(trcn_RFL0(ids_nut_beg:ids_nuts_end,JY,JX));      trcn_RFL0=0._r8
  allocate(trcs_RFL1(ids_beg:ids_end,JY,JX));      trcs_RFL1=0._r8
  allocate(DOM_XPoreTransp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_XPoreTransp_flx=0._r8
  allocate(RCOFXS(JZ,JY,JX));   RCOFXS=0._r8
  allocate(RCHFXS(JZ,JY,JX));   RCHFXS=0._r8
  allocate(ROXFXS(JZ,JY,JX));   ROXFXS=0._r8
  allocate(RNGFXS(JZ,JY,JX));   RNGFXS=0._r8
  allocate(RN2FXS(JZ,JY,JX));   RN2FXS=0._r8
  allocate(RN4FXW(JZ,JY,JX));   RN4FXW=0._r8
  allocate(RN3FXW(JZ,JY,JX));   RN3FXW=0._r8
  allocate(RNOFXW(JZ,JY,JX));   RNOFXW=0._r8
  allocate(RH2PXS(JZ,JY,JX));   RH2PXS=0._r8
  allocate(RN4FXB(JZ,JY,JX));   RN4FXB=0._r8
  allocate(RN3FXB(JZ,JY,JX));   RN3FXB=0._r8
  allocate(RNOFXB(JZ,JY,JX));   RNOFXB=0._r8
  allocate(RH2BXB(JZ,JY,JX));   RH2BXB=0._r8
  allocate(RNXFXS(JZ,JY,JX));   RNXFXS=0._r8
  allocate(RNXFXB(JZ,JY,JX));   RNXFXB=0._r8
  allocate(RH1PXS(JZ,JY,JX));   RH1PXS=0._r8
  allocate(RH1BXB(JZ,JY,JX));   RH1BXB=0._r8
  allocate(CLSGL2(0:JZ,JY,JX)); CLSGL2=0._r8
  allocate(CQSGL2(0:JZ,JY,JX)); CQSGL2=0._r8
  allocate(POSGL2(0:JZ,JY,JX)); POSGL2=0._r8
  allocate(OLSGL2(0:JZ,JY,JX)); OLSGL2=0._r8
  allocate(ZNSGL2(0:JZ,JY,JX)); ZNSGL2=0._r8
  allocate(ZLSGL2(0:JZ,JY,JX)); ZLSGL2=0._r8
  allocate(ZVSGL2(0:JZ,JY,JX)); ZVSGL2=0._r8
  allocate(HLSGL2(0:JZ,JY,JX)); HLSGL2=0._r8
  allocate(ZOSGL2(0:JZ,JY,JX)); ZOSGL2=0._r8
  allocate(RNBDFG(0:JZ,JY,JX)); RNBDFG=0._r8

  allocate(dom_TFloXSurRunoff(idom_beg:idom_end,1:jcplx,JY,JX)); dom_TFloXSurRunoff=0._r8
  allocate(Difus_Macp_flx_DOM(idom_beg:idom_end,1:jcplx));Difus_Macp_flx_DOM=0._r8
  allocate(TQSCOS(JY,JX));      TQSCOS=0._r8
  allocate(trcg_SnowDrift(idg_beg:idg_NH3,JY,JX));      trcg_SnowDrift=0._r8
  allocate(trcn_SnowDrift(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_SnowDrift=0._r8
  allocate(trcg_TFloXSurRunoff(idg_beg:idg_NH3,JY,JX)); trcg_TFloXSurRunoff=0._r8
  allocate(trcn_TFloXSurRunoff(ids_nut_beg:ids_nuts_end,JY,JX));trcn_TFloXSurRunoff=0._r8
  allocate(TN4FLW(JZ,JY,JX));   TN4FLW=0._r8
  allocate(TN3FLW(JZ,JY,JX));   TN3FLW=0._r8
  allocate(TNOFLW(JZ,JY,JX));   TNOFLW=0._r8
  allocate(TH2PFS(JZ,JY,JX));   TH2PFS=0._r8
  allocate(TQRH1P(JY,JX));      TQRH1P=0._r8
  allocate(TH1PFS(JZ,JY,JX));   TH1PFS=0._r8
  allocate(CDOM_MicP1(idom_beg:idom_end,1:jcplx)); CDOM_MicP1=0._r8
  allocate(CDOM_MicP2(idom_beg:idom_end,1:jcplx)); CDOM_MicP2=0._r8

  allocate(GasDifc_vrc(idg_beg:idg_end,JZ,JY,JX))
  allocate(SolDifc_vrc(ids_beg:ids_end,0:JZ,JY,JX));SolDifc_vrc=0._r8
  allocate(DifuscG(idg_beg:idg_end,3,JZ,JY,JX)); DifuscG=0._r8
  allocate(DCO2G(3,JZ,JY,JX));  DCO2G=0._r8
  allocate(DCH4G(3,JZ,JY,JX));  DCH4G=0._r8
  allocate(DOXYG(3,JZ,JY,JX));  DOXYG=0._r8
  allocate(DZ2GG(3,JZ,JY,JX));  DZ2GG=0._r8
  allocate(DZ2OG(3,JZ,JY,JX));  DZ2OG=0._r8
  allocate(DNH3G(3,JZ,JY,JX));  DNH3G=0._r8
  allocate(trcg_VLWatMicP(idg_beg:idg_end,0:JZ,JY,JX)); trcg_VLWatMicP=0._r8
  allocate(HGSGL2(JZ,JY,JX));   HGSGL2=0._r8
  allocate(DH2GG(3,JZ,JY,JX));  DH2GG=0._r8
  allocate(RHGFXS(JZ,JY,JX));   RHGFXS=0._r8
  allocate(THGFLG(JZ,JY,JX));   THGFLG=0._r8
  allocate(CumReductVLsoiAirPM(JZ,JY,JX));     CumReductVLsoiAirPM=0._r8
  allocate(THETH2(JZ,JY,JX));   THETH2=0._r8
  allocate(THETHL(JZ,JY,JX));   THETHL=0._r8
  allocate(VLsoiAirPMA(JZ,JY,JX));   VLsoiAirPMA=0._r8
  allocate(VLsoiAirPMB(JZ,JY,JX));   VLsoiAirPMB=0._r8
  allocate(VLWatMicPMA(JZ,JY,JX));   VLWatMicPMA=0._r8
  allocate(VLWatMicPMB(JZ,JY,JX));   VLWatMicPMB=0._r8
  allocate(VLWatMicPXA(0:JZ,JY,JX)); VLWatMicPXA=0._r8
  allocate(VLWatMicPXB(JZ,JY,JX));   VLWatMicPXB=0._r8
  allocate(PARG_cef(idg_beg:idg_NH3,JY,JX));      PARG_cef=0._r8

  allocate(RBGCSinkG(idg_beg:idg_end,0:JZ,JY,JX));RBGCSinkG=0._r8
  allocate(RBGCSinkS(ids_nuts_beg:ids_nuts_end,0:JZ,JY,JX));RBGCSinkS=0._r8

  allocate(RCOSK2(0:JZ,JY,JX)); RCOSK2=0._r8
  allocate(ROXSK2(0:JZ,JY,JX)); ROXSK2=0._r8
  allocate(RCHSK2(0:JZ,JY,JX)); RCHSK2=0._r8
  allocate(RNGSK2(0:JZ,JY,JX)); RNGSK2=0._r8
  allocate(RN2SK2(0:JZ,JY,JX)); RN2SK2=0._r8
  allocate(RN4SK2(0:JZ,JY,JX)); RN4SK2=0._r8
  allocate(RN3SK2(0:JZ,JY,JX)); RN3SK2=0._r8
  allocate(RNOSK2(0:JZ,JY,JX)); RNOSK2=0._r8
  allocate(RHPSK2(0:JZ,JY,JX)); RHPSK2=0._r8
  allocate(R4BSK2(JZ,JY,JX));   R4BSK2=0._r8
  allocate(R3BSK2(JZ,JY,JX));   R3BSK2=0._r8
  allocate(RNBSK2(JZ,JY,JX));   RNBSK2=0._r8
  allocate(RHBSK2(JZ,JY,JX));   RHBSK2=0._r8
  allocate(RNXSK2(0:JZ,JY,JX)); RNXSK2=0._r8
  allocate(RNZSK2(JZ,JY,JX));   RNZSK2=0._r8
  allocate(RHGSK2(0:JZ,JY,JX)); RHGSK2=0._r8
  allocate(RNHSK2(0:JZ,JY,JX)); RNHSK2=0._r8
  allocate(R1PSK2(0:JZ,JY,JX)); R1PSK2=0._r8

  allocate(TN3FLG(JZ,JY,JX));   TN3FLG=0._r8
  allocate(trcg_BBL(idg_beg:idg_end,JZ,JY,JX));   trcg_BBL=0._r8
  allocate(trcg_FloXSurRunoff(idg_beg:idg_end,JV,JH));     trcg_FloXSurRunoff=0._r8
  allocate(trcn_FloXSurRunoff(ids_nut_beg:ids_nuts_end,JV,JH));     trcn_FloXSurRunoff=0._r8
  allocate(DOM_MacP2(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_MacP2=0._r8
  allocate(DOM_Transp2Macp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flx=0._r8
  allocate(TN4FHW(JZ,JY,JX));   TN4FHW=0._r8
  allocate(TN3FHW(JZ,JY,JX));   TN3FHW=0._r8
  allocate(TNOFHW(JZ,JY,JX));   TNOFHW=0._r8
  allocate(TH2PHS(JZ,JY,JX));   TH2PHS=0._r8
  allocate(ZNO2B2(JZ,JY,JX));   ZNO2B2=0._r8
  allocate(ZN2BH2(JZ,JY,JX));   ZN2BH2=0._r8
  allocate(H1P4H2(JZ,JY,JX));   H1P4H2=0._r8
  allocate(H1PBH2(JZ,JY,JX));   H1PBH2=0._r8
  allocate(TH1PHS(JZ,JY,JX));   TH1PHS=0._r8
  allocate(H1PO42(0:JZ,JY,JX)); H1PO42=0._r8
  allocate(H1POB2(0:JZ,JY,JX)); H1POB2=0._r8
  allocate(OCSGL2(0:JZ,JY,JX)); OCSGL2=0._r8
  allocate(ONSGL2(0:JZ,JY,JX)); ONSGL2=0._r8
  allocate(OPSGL2(0:JZ,JY,JX)); OPSGL2=0._r8
  allocate(OASGL2(0:JZ,JY,JX)); OASGL2=0._r8

  allocate(trcg_solsml2(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml2=0._r8
  allocate(trcn_solsml2(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml2=0._r8

  allocate(CO2W2(JS,JY,JX));    CO2W2=0._r8
  allocate(CH4W2(JS,JY,JX));    CH4W2=0._r8
  allocate(OXYW2(JS,JY,JX));    OXYW2=0._r8
  allocate(ZNGW2(JS,JY,JX));    ZNGW2=0._r8
  allocate(ZN2W2(JS,JY,JX));    ZN2W2=0._r8
  allocate(ZN4W2(JS,JY,JX));    ZN4W2=0._r8
  allocate(ZN3W2(JS,JY,JX));    ZN3W2=0._r8
  allocate(ZNOW2(JS,JY,JX));    ZNOW2=0._r8
  allocate(ZHPW2(JS,JY,JX));    ZHPW2=0._r8
  allocate(Z1PW2(JS,JY,JX));    Z1PW2=0._r8

  allocate(trcg_TBLS(idg_beg:idg_NH3,JS,JY,JX));   trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_TBLS=0._r8

  allocate(dom_2DFloXSurRunoffM(idom_beg:idom_end,1:jcplx,2,2,JV,JH));dom_2DFloXSurRunoffM=0._r8
  allocate(trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,2,2,JV,JH));  trcg_2DFloXSurRunoffM=0._r8
  allocate(RQRHGS(2,2,JV,JH));  RQRHGS=0._r8
  allocate(RQRNH4(2,2,JV,JH));  RQRNH4=0._r8
  allocate(RCODFS(JY,JX));      RCODFS=0._r8
  allocate(RCHDFS(JY,JX));      RCHDFS=0._r8
  allocate(ROXDFS(JY,JX));      ROXDFS=0._r8
  allocate(RNGDFS(JY,JX));      RNGDFS=0._r8
  allocate(RN2DFS(JY,JX));      RN2DFS=0._r8
  allocate(RN3DFS(JY,JX));      RN3DFS=0._r8
  allocate(RNBDFS(JY,JX));      RNBDFS=0._r8
  allocate(RHGDFS(JY,JX));      RHGDFS=0._r8
  allocate(RDFR_gas(idg_beg:idg_NH3,JY,JX));      RDFR_gas=0._r8
  allocate(R1BSK2(JZ,JY,JX));   R1BSK2=0._r8

  allocate(ZNH4H2(JZ,JY,JX));   ZNH4H2=0._r8
  allocate(ZN4BH2(JZ,JY,JX));   ZN4BH2=0._r8
  allocate(ZNH3H2(JZ,JY,JX));   ZNH3H2=0._r8
  allocate(ZN3BH2(JZ,JY,JX));   ZN3BH2=0._r8
  allocate(ZNO3H2(JZ,JY,JX));   ZNO3H2=0._r8
  allocate(ZNOBH2(JZ,JY,JX));   ZNOBH2=0._r8
  allocate(H2P4H2(JZ,JY,JX));   H2P4H2=0._r8
  allocate(H2PBH2(JZ,JY,JX));   H2PBH2=0._r8
  allocate(ZNO2H2(JZ,JY,JX));   ZNO2H2=0._r8

  allocate(trc_gasml_vr2(idg_beg:idg_end,0:JZ,JY,JX)); trc_gasml_vr2(:,:,:,:)=0._r8
  allocate(trc_solml_vr2(ids_beg:ids_end,0:JZ,JY,JX)); trc_solml_vr2(:,:,:,:)=0._r8
  allocate(trc_soHml2(ids_beg:ids_end,0:JZ,JY,JX)); trc_soHml2(:,:,:,:)=0._r8

  allocate(trcg_2DSnowDrift(idg_beg:idg_NH3,2,JV,JH));    trcg_2DSnowDrift=0._r8
  allocate(trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_2DFloXSurRunoffM=0._r8
  allocate(R3GasADFlx(idg_beg:idg_NH3,3,JD,JV,JH));R3GasADFlx=0._r8
  allocate(Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,JZ,JY,JH));Gas_AdvDif_Flx_vr=0._r8

  allocate(RNXFHS(3,JD,JV,JH)); RNXFHS=0._r8
  allocate(RNXFHB(3,JD,JV,JH)); RNXFHB=0._r8
  allocate(RH1PFS(3,0:JD,JV,JH));RH1PFS=0._r8
  allocate(RH1BFB(3,0:JD,JV,JH));RH1BFB=0._r8
  allocate(RH1PHS(3,JD,JV,JH)); RH1PHS=0._r8
  allocate(RH1BHB(3,JD,JV,JH)); RH1BHB=0._r8
  allocate(RN3FHB(3,JD,JV,JH)); RN3FHB=0._r8
  allocate(RNOFHB(3,JD,JV,JH)); RNOFHB=0._r8
  allocate(RH2BHB(3,JD,JV,JH)); RH2BHB=0._r8
  allocate(RNOFHW(3,JD,JV,JH)); RNOFHW=0._r8
  allocate(RH2PHS(3,JD,JV,JH)); RH2PHS=0._r8
  allocate(RN4FHB(3,JD,JV,JH)); RN4FHB=0._r8
  allocate(RN2FHS(3,JD,JV,JH)); RN2FHS=0._r8
  allocate(RN4FHW(3,JD,JV,JH)); RN4FHW=0._r8
  allocate(RN3FHW(3,JD,JV,JH)); RN3FHW=0._r8
  allocate(DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx));DOM_Adv2MicP_flx=0._r8

  allocate(RGasDSFlx(idg_beg:idg_end,0:JZ,JY,JX)); RGasDSFlx=0._r8
  allocate(RGasSSVol(idg_beg:idg_end,JY,JX)); RGasSSVol=0._r8
  allocate(RHGFHS(3,JD,JV,JH)); RHGFHS=0._r8
  allocate(RCHFHS(3,JD,JV,JH)); RCHFHS=0._r8
  allocate(ROXFHS(3,JD,JV,JH)); ROXFHS=0._r8
  allocate(RNGFHS(3,JD,JV,JH)); RNGFHS=0._r8
  allocate(trcn_2DSnowDrift(ids_nut_beg:ids_nuts_end,2,JV,JH)); trcn_2DSnowDrift=0._r8
  allocate(RCOFLS(3,0:JD,JV,JH));RCOFLS=0._r8
  allocate(RCHFLS(3,0:JD,JV,JH));RCHFLS=0._r8
  allocate(ROXFLS(3,0:JD,JV,JH));ROXFLS=0._r8
  allocate(RNGFLS(3,0:JD,JV,JH));RNGFLS=0._r8
  allocate(RN2FLS(3,0:JD,JV,JH));RN2FLS=0._r8
  allocate(RHGFLS(3,0:JD,JV,JH));RHGFLS=0._r8
  allocate(RN4FLW(3,0:JD,JV,JH));RN4FLW=0._r8
  allocate(RN3FLW(3,0:JD,JV,JH));RN3FLW=0._r8
  allocate(RNOFLW(3,0:JD,JV,JH));RNOFLW=0._r8
  allocate(RNXFLS(3,0:JD,JV,JH));RNXFLS=0._r8
  allocate(RH2PFS(3,0:JD,JV,JH));RH2PFS=0._r8
  allocate(RN4FLB(3,0:JD,JV,JH));RN4FLB=0._r8
  allocate(RN3FLB(3,0:JD,JV,JH));RN3FLB=0._r8
  allocate(RNOFLB(3,0:JD,JV,JH));RNOFLB=0._r8
  allocate(RNXFLB(3,0:JD,JV,JH));RNXFLB=0._r8
  allocate(RH2BFB(3,0:JD,JV,JH));RH2BFB=0._r8
  allocate(RCOFHS(3,JD,JV,JH)); RCOFHS=0._r8

  allocate(R3PoreSoHFlx(ids_beg:ids_end,3,JD,JV,JH));R3PoreSoHFlx=0._r8
  allocate(R3PoreSolFlx(ids_beg:ids_end,3,0:JD,JV,JH));R3PoreSolFlx=0._r8
  allocate(RporeSoXFlx(ids_beg:ids_end,JZ,JY,JX));RporeSoXFlx=0._r8
  allocate(R3PorTSolFlx(ids_beg:ids_end,JZ,JY,JX));R3PorTSolFlx=0._r8
  allocate(R3PorTSoHFlx(ids_beg:ids_end,JZ,JY,JX));R3PorTSoHFlx=0._r8
  allocate(trcg_VFloSnow(idg_beg:idg_end));trcg_VFloSnow=0._r8
  allocate(trcn_band_VFloSnow(ids_nutb_beg:ids_nutb_end));trcn_band_VFloSnow=0._r8
  allocate(trcn_soil_VFloSnow(ids_nuts_beg:ids_nuts_end));trcn_soil_VFloSnow=0._r8
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
  call destroy(ROCFL0)
  call destroy(RONFL0)
  call destroy(ROPFL0)
  call destroy(ROAFL0)
  call destroy(ROCFL1)
  call destroy(RONFL1)
  call destroy(ROPFL1)
  call destroy(ROAFL1)

  call destroy(trcs_RFL1)

  call destroy(DOM_XPoreTransp_flx)
  call destroy(RCOFXS)
  call destroy(RCHFXS)
  call destroy(ROXFXS)
  call destroy(RNGFXS)
  call destroy(RN2FXS)
  call destroy(RN4FXW)
  call destroy(RN3FXW)
  call destroy(RNOFXW)
  call destroy(RH2PXS)
  call destroy(RN4FXB)
  call destroy(RN3FXB)
  call destroy(RNOFXB)
  call destroy(RH2BXB)
  call destroy(RNXFXS)
  call destroy(RNXFXB)
  call destroy(RH1PXS)
  call destroy(RH1BXB)

  call destroy(trcg_BBL)
  call destroy(CLSGL2)
  call destroy(CQSGL2)
  call destroy(POSGL2)
  call destroy(OLSGL2)
  call destroy(ZNSGL2)
  call destroy(ZLSGL2)
  call destroy(ZVSGL2)
  call destroy(HLSGL2)
  call destroy(ZOSGL2)

  call destroy(RNBDFG)
  call destroy(dom_TFloXSurRunoff)
  call destroy(trcg_TFloXSurRunoff)
  call destroy(trcn_TFloXSurRunoff)
  call destroy(trcg_SnowDrift)
  call destroy(trcn_SnowDrift)

  call destroy(TN4FLW)
  call destroy(TN3FLW)
  call destroy(TNOFLW)
  call destroy(TH2PFS)
  call destroy(TQRH1P)
  call destroy(TH1PFS)
  call destroy(DCO2G)
  call destroy(DCH4G)
  call destroy(DOXYG)
  call destroy(DZ2GG)
  call destroy(DZ2OG)
  call destroy(DNH3G)
  call destroy(trcg_VLWatMicP)
  call destroy(HGSGL2)
  call destroy(DH2GG)
  call destroy(RHGFXS)
  call destroy(THGFLG)
  call destroy(CumReductVLsoiAirPM)
  call destroy(THETH2)
  call destroy(THETHL)
  call destroy(VLsoiAirPMA)
  call destroy(VLsoiAirPMB)
  call destroy(VLWatMicPMA)
  call destroy(VLWatMicPMB)
  call destroy(VLWatMicPXA)
  call destroy(VLWatMicPXB)
  call destroy(PARG_cef)
  call destroy(RCOSK2)
  call destroy(ROXSK2)
  call destroy(RCHSK2)
  call destroy(RNGSK2)
  call destroy(RN2SK2)
  call destroy(RN4SK2)
  call destroy(RN3SK2)
  call destroy(RNOSK2)
  call destroy(RHPSK2)
  call destroy(R4BSK2)
  call destroy(R3BSK2)
  call destroy(RNBSK2)
  call destroy(RHBSK2)
  call destroy(RNXSK2)
  call destroy(RNZSK2)
  call destroy(RHGSK2)
  call destroy(RNHSK2)
  call destroy(R1PSK2)
  call destroy(TN3FLG)
  call destroy(TN4FHW)
  call destroy(TN3FHW)
  call destroy(TNOFHW)
  call destroy(TH2PHS)
  call destroy(ZNO2B2)
  call destroy(ZN2BH2)
  call destroy(H1P4H2)
  call destroy(H1PBH2)
  call destroy(TH1PHS)
  call destroy(H1PO42)
  call destroy(H1POB2)
  call destroy(OCSGL2)
  call destroy(ONSGL2)
  call destroy(OPSGL2)
  call destroy(OASGL2)

  call destroy(trcg_solsml2)
  call destroy(trcn_solsml2)

  call destroy(CO2W2)
  call destroy(CH4W2)
  call destroy(OXYW2)
  call destroy(ZNGW2)
  call destroy(ZN2W2)
  call destroy(ZN4W2)
  call destroy(ZN3W2)
  call destroy(ZNOW2)
  call destroy(ZHPW2)
  call destroy(Z1PW2)

  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)
  call destroy(trcg_RFL0)
  call destroy(trcg_2DSnowDrift)
  call destroy(Difus_Micp_flx_DOM)
  call destroy(Difus_Macp_flx_DOM)

  call destroy(RQRHGS)
  call destroy(RQRNH4)
  call destroy(RCODFS)
  call destroy(RCHDFS)
  call destroy(ROXDFS)
  call destroy(RNGDFS)
  call destroy(RN2DFS)
  call destroy(RN3DFS)
  call destroy(RNBDFS)
  call destroy(RHGDFS)
  call destroy(RDFR_gas)
  call destroy(R1BSK2)

  call destroy(ZNH4H2)
  call destroy(ZN4BH2)
  call destroy(ZNH3H2)
  call destroy(ZN3BH2)
  call destroy(ZNO3H2)
  call destroy(ZNOBH2)
  call destroy(H2P4H2)
  call destroy(H2PBH2)
  call destroy(ZNO2H2)
  call destroy(DOM_MacP2)

  call destroy(trc_gasml_vr2)
  call destroy(trc_solml_vr2)
  call destroy(trc_soHml2)
  call destroy(trcn_band_VFloSnow)
  call destroy(trcn_soil_VFloSnow)

  call destroy(RNXFHS)
  call destroy(RNXFHB)
  call destroy(RH1PFS)
  call destroy(RH1BFB)
  call destroy(RH1PHS)
  call destroy(RH1BHB)
  call destroy(RN3FHB)
  call destroy(RNOFHB)
  call destroy(RH2BHB)
  call destroy(RNOFHW)
  call destroy(RH2PHS)
  call destroy(RN4FHB)
  call destroy(RN2FHS)
  call destroy(RN4FHW)
  call destroy(RN3FHW)

  call destroy(RGasDSFlx)
  call destroy(RHGFHS)
  call destroy(RCHFHS)
  call destroy(ROXFHS)
  call destroy(RNGFHS)
  call destroy(trcn_2DSnowDrift)
  call destroy(RCOFLS)
  call destroy(RCHFLS)
  call destroy(ROXFLS)
  call destroy(RNGFLS)
  call destroy(RN2FLS)
  call destroy(RHGFLS)
  call destroy(RN4FLW)
  call destroy(RN3FLW)
  call destroy(RNOFLW)
  call destroy(RNXFLS)
  call destroy(RH2PFS)
  call destroy(RN4FLB)
  call destroy(RN3FLB)
  call destroy(RNOFLB)
  call destroy(RNXFLB)
  call destroy(RH2BFB)
  call destroy(RCOFHS)
  call destroy(R3PoreSoHFlx)
  call destroy(R3PoreSolFlx)
  call destroy(RporeSoXFlx)
  call destroy(CDOM_MacP1)
  call destroy(CDOM_MacP2)
  call destroy(DOM_Adv2MicP_flx)
  call destroy(DOM_Adv2MacP_flx)
  call destroy(DOM_Transp2Macp_flx)
  call destroy(DOM_Transp2Micp_flx)
  end subroutine DestructTransfrData

end module TranspNoSaltDataMod
