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

  real(r8), allocatable :: RFLOC(:)
  real(r8), allocatable :: RFLON(:)
  real(r8), allocatable :: RFLOP(:)
  real(r8), allocatable :: RFLOA(:)
  real(r8), allocatable :: DFVOC(:)
  real(r8), allocatable :: DFVON(:)
  real(r8), allocatable :: DFVOP(:)
  real(r8), allocatable :: DFVOA(:)
  real(r8), allocatable :: DFHOC(:)
  real(r8), allocatable :: DFHON(:)
  real(r8), allocatable :: DFHOP(:)
  real(r8), allocatable :: DFHOA(:)
  real(r8), allocatable :: COQC1(:)
  real(r8), allocatable :: COQC2(:)
  real(r8), allocatable :: COQN1(:)
  real(r8), allocatable :: COQN2(:)
  real(r8), allocatable :: COQP1(:)
  real(r8), allocatable :: COQP2(:)
  real(r8), allocatable :: COQA1(:)
  real(r8), allocatable :: COQA2(:)
  real(r8), allocatable :: OQC2(:,:,:,:)
  real(r8), allocatable :: OQN2(:,:,:,:)
  real(r8), allocatable :: OQP2(:,:,:,:)
  real(r8), allocatable :: OQA2(:,:,:,:)
  real(r8), allocatable :: ROCFLS(:,:,:,:,:)
  real(r8), allocatable :: RONFLS(:,:,:,:,:)
  real(r8), allocatable :: ROPFLS(:,:,:,:,:)
  real(r8), allocatable :: ROAFLS(:,:,:,:,:)
  real(r8), allocatable :: ROCFHS(:,:,:,:,:)
  real(r8), allocatable :: RONFHS(:,:,:,:,:)
  real(r8), allocatable :: ROPFHS(:,:,:,:,:)
  real(r8), allocatable :: ROAFHS(:,:,:,:,:)
  real(r8), allocatable :: TOCFLS(:,:,:,:)
  real(r8), allocatable :: TONFLS(:,:,:,:)
  real(r8), allocatable :: TOPFLS(:,:,:,:)
  real(r8), allocatable :: TOAFLS(:,:,:,:)
  real(r8), allocatable :: ROCSK2(:,:,:,:)
  real(r8), allocatable :: RONSK2(:,:,:,:)
  real(r8), allocatable :: ROPSK2(:,:,:,:)
  real(r8), allocatable :: ROASK2(:,:,:,:)

  real(r8), allocatable :: RQROC0(:,:,:)
  real(r8), allocatable :: RQRON0(:,:,:)
  real(r8), allocatable :: RQROA0(:,:,:)
  real(r8), allocatable :: RQROP0(:,:,:)

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
  real(r8), allocatable ::  ROCFXS(:,:,:,:)                    !
  real(r8), allocatable ::  RONFXS(:,:,:,:)                    !
  real(r8), allocatable ::  ROPFXS(:,:,:,:)                    !
  real(r8), allocatable ::  ROAFXS(:,:,:,:)                    !
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
  real(r8), allocatable ::  TQROC(:,:,:)                       !
  real(r8), allocatable ::  TQRON(:,:,:)                       !
  real(r8), allocatable ::  TQROP(:,:,:)                       !
  real(r8), allocatable ::  TQROA(:,:,:)                       !

  real(r8), allocatable ::  trcg_TQ(:,:,:)                        !
  real(r8), allocatable ::  trcn_TQ(:,:,:)
  real(r8), allocatable ::  trcg_TQR(:,:,:)
  real(r8), allocatable ::  trcn_TQR(:,:,:)
  real(r8), allocatable ::  R3PorTSolFlx(:,:,:,:) !total 3D micropore flux
  real(r8), allocatable ::  R3PorTSoHFlx(:,:,:,:) !total 3D macropore flux

  real(r8), allocatable ::  TN4FLW(:,:,:)                      !
  real(r8), allocatable ::  TN3FLW(:,:,:)                      !
  real(r8), allocatable ::  TNOFLW(:,:,:)                      !
  real(r8), allocatable ::  TH2PFS(:,:,:)                      !
  real(r8), allocatable ::  TQRH1P(:,:)                        !
  real(r8), allocatable ::  TH1PFS(:,:,:)                      !

  real(r8), allocatable ::  GasDifcc(:,:,:,:)
  real(r8), allocatable ::  SolDifcc(:,:,:,:)
  real(r8), allocatable ::  DifuscG(:,:,:,:,:)
  real(r8), allocatable ::  DCO2G(:,:,:,:)                     !
  real(r8), allocatable ::  DCH4G(:,:,:,:)                     !
  real(r8), allocatable ::  DOXYG(:,:,:,:)                     !
  real(r8), allocatable ::  DZ2GG(:,:,:,:)                     !
  real(r8), allocatable ::  DZ2OG(:,:,:,:)                     !
  real(r8), allocatable ::  DNH3G(:,:,:,:)                     !
  real(r8), allocatable ::  VLWatMicPCO(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPCH(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPOX(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPNG(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPN2(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPN3(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPNB(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPHG(:,:,:)                      !
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
  real(r8), allocatable ::  trcg_RQR0(:,:,:)                       !

  real(r8), allocatable ::  trcn_RQR0(:,:,:)                       !
  real(r8), allocatable ::  OQNH2(:,:,:,:)                     !
  real(r8), allocatable ::  OQPH2(:,:,:,:)                     !
  real(r8), allocatable ::  OQAH2(:,:,:,:)                     !
  real(r8), allocatable ::  TOCFHS(:,:,:,:)                    !
  real(r8), allocatable ::  TONFHS(:,:,:,:)                    !
  real(r8), allocatable ::  TOPFHS(:,:,:,:)                    !
  real(r8), allocatable ::  TOAFHS(:,:,:,:)                    !
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
  real(r8), allocatable ::  RQROP(:,:,:,:,:)                   !
  real(r8), allocatable ::  RQROA(:,:,:,:,:)                   !
  real(r8), allocatable ::  trcg_RQR(:,:,:,:,:)                    !
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
  real(r8), allocatable ::  RQROC(:,:,:,:,:)                   !
  real(r8), allocatable ::  RQRON(:,:,:,:,:)                   !


  real(r8), allocatable ::  ZNH4H2(:,:,:)                      !
  real(r8), allocatable ::  ZN4BH2(:,:,:)                      !
  real(r8), allocatable ::  ZNH3H2(:,:,:)                      !
  real(r8), allocatable ::  ZN3BH2(:,:,:)                      !
  real(r8), allocatable ::  ZNO3H2(:,:,:)                      !
  real(r8), allocatable ::  ZNOBH2(:,:,:)                      !
  real(r8), allocatable ::  H2P4H2(:,:,:)                      !
  real(r8), allocatable ::  H2PBH2(:,:,:)                      !
  real(r8), allocatable ::  ZNO2H2(:,:,:)                      !
  real(r8), allocatable ::  OQCH2(:,:,:,:)                     !

  real(r8), allocatable :: trc_gasml2(:,:,:,:)
  real(r8), allocatable :: trc_solml2(:,:,:,:)
  real(r8), allocatable :: trc_soHml2(:,:,:,:)

  real(r8), allocatable ::  RQSH2P(:,:,:)                      !
  real(r8), allocatable ::  RQSH1P(:,:,:)                      !
  real(r8), allocatable ::  trcg_RQS(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_RQR(:,:,:,:,:)                    !

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
  real(r8), allocatable ::  RTGasADFlx(:,:,:,:)  !total 3D gas flux advection + diffusion

  real(r8), allocatable ::  RHGFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RCHFHS(:,:,:,:)                    !
  real(r8), allocatable ::  ROXFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RNGFHS(:,:,:,:)                    !
  real(r8), allocatable ::  RQSNH4(:,:,:)                      !
  real(r8), allocatable ::  RQSNO3(:,:,:)                      !
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
  allocate(RQROC0(1:jcplx,JV,JH));       RQROC0=0._r8
  allocate(RQRON0(1:jcplx,JV,JH));       RQRON0=0._r8
  allocate(RQROA0(1:jcplx,JV,JH));       RQROA0=0._r8
  allocate(RQROP0(1:jcplx,JV,JH));       RQROP0=0._r8
  allocate(OQC2(1:jcplx,0:JZ,JY,JX));    OQC2=0._r8
  allocate(OQN2(1:jcplx,0:JZ,JY,JX));    OQN2=0._r8
  allocate(OQP2(1:jcplx,0:JZ,JY,JX));    OQP2=0._r8
  allocate(OQA2(1:jcplx,0:JZ,JY,JX));    OQA2=0._r8
  allocate(ROCSK2(1:jcplx,0:JZ,JY,JX));  ROCSK2=0._r8
  allocate(RONSK2(1:jcplx,0:JZ,JY,JX));  RONSK2=0._r8
  allocate(ROPSK2(1:jcplx,0:JZ,JY,JX));  ROPSK2=0._r8
  allocate(ROASK2(1:jcplx,0:JZ,JY,JX));  ROASK2=0._r8
  allocate(ROCFLS(1:jcplx,3,0:JD,JV,JH));ROCFLS=0._r8
  allocate(RONFLS(1:jcplx,3,0:JD,JV,JH));RONFLS=0._r8
  allocate(ROPFLS(1:jcplx,3,0:JD,JV,JH));ROPFLS=0._r8
  allocate(ROAFLS(1:jcplx,3,0:JD,JV,JH));ROAFLS=0._r8
  allocate(ROCFHS(1:jcplx,3,JD,JV,JH));  ROCFHS=0._r8
  allocate(RONFHS(1:jcplx,3,JD,JV,JH));  RONFHS=0._r8
  allocate(ROPFHS(1:jcplx,3,JD,JV,JH));  ROPFHS=0._r8
  allocate(ROAFHS(1:jcplx,3,JD,JV,JH));  ROAFHS=0._r8
  allocate(TOCFLS(1:jcplx,JZ,JY,JX));    TOCFLS=0._r8
  allocate(TONFLS(1:jcplx,JZ,JY,JX));    TONFLS=0._r8
  allocate(TOPFLS(1:jcplx,JZ,JY,JX));    TOPFLS=0._r8
  allocate(TOAFLS(1:jcplx,JZ,JY,JX));    TOAFLS=0._r8

  allocate(trcg_RBLS(idg_beg:idg_end-1,JS,JY,JX));   trcg_RBLS=0._r8
  allocate(trcn_RBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX));   trcn_RBLS=0._r8
  allocate(ROCFL0(1:NumOfLitrCmplxs,JY,JX));  ROCFL0=0._r8
  allocate(RONFL0(1:NumOfLitrCmplxs,JY,JX));  RONFL0=0._r8
  allocate(ROPFL0(1:NumOfLitrCmplxs,JY,JX));  ROPFL0=0._r8
  allocate(ROAFL0(1:NumOfLitrCmplxs,JY,JX));  ROAFL0=0._r8
  allocate(ROCFL1(1:NumOfLitrCmplxs,JY,JX));  ROCFL1=0._r8
  allocate(RONFL1(1:NumOfLitrCmplxs,JY,JX));  RONFL1=0._r8
  allocate(ROPFL1(1:NumOfLitrCmplxs,JY,JX));  ROPFL1=0._r8
  allocate(ROAFL1(1:NumOfLitrCmplxs,JY,JX));  ROAFL1=0._r8
  allocate(trcg_RFL0(idg_beg:idg_end-1,JY,JX));      trcg_RFL0=0._r8
  allocate(trcn_RFL0(ids_nut_beg:ids_nuts_end,JY,JX));      trcn_RFL0=0._r8
  allocate(trcs_RFL1(ids_beg:ids_end,JY,JX));      trcs_RFL1=0._r8
  allocate(ROCFXS(1:jcplx,JZ,JY,JX));ROCFXS=0._r8
  allocate(RONFXS(1:jcplx,JZ,JY,JX));RONFXS=0._r8
  allocate(ROPFXS(1:jcplx,JZ,JY,JX));ROPFXS=0._r8
  allocate(ROAFXS(1:jcplx,JZ,JY,JX));ROAFXS=0._r8
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

  allocate(TQROC(1:jcplx,JY,JX));   TQROC=0._r8
  allocate(TQRON(1:jcplx,JY,JX));   TQRON=0._r8
  allocate(TQROP(1:jcplx,JY,JX));   TQROP=0._r8
  allocate(TQROA(1:jcplx,JY,JX));   TQROA=0._r8

  allocate(TQSCOS(JY,JX));      TQSCOS=0._r8
  allocate(trcg_TQ(idg_beg:idg_end-1,JY,JX));      trcg_TQ=0._r8
  allocate(trcn_TQ(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_TQ=0._r8
  allocate(trcg_TQR(idg_beg:idg_end-1,JY,JX)); trcg_TQR=0._r8
  allocate(trcn_TQR(ids_nut_beg:ids_nuts_end,JY,JX));trcn_TQR=0._r8
  allocate(TN4FLW(JZ,JY,JX));   TN4FLW=0._r8
  allocate(TN3FLW(JZ,JY,JX));   TN3FLW=0._r8
  allocate(TNOFLW(JZ,JY,JX));   TNOFLW=0._r8
  allocate(TH2PFS(JZ,JY,JX));   TH2PFS=0._r8
  allocate(TQRH1P(JY,JX));      TQRH1P=0._r8
  allocate(TH1PFS(JZ,JY,JX));   TH1PFS=0._r8

  allocate(GasDifcc(idg_beg:idg_end,JZ,JY,JX))
  allocate(SolDifcc(ids_beg:ids_end,0:JZ,JY,JX));SolDifcc=0._r8
  allocate(DifuscG(idg_beg:idg_end,3,JZ,JY,JX)); DifuscG=0._r8
  allocate(DCO2G(3,JZ,JY,JX));  DCO2G=0._r8
  allocate(DCH4G(3,JZ,JY,JX));  DCH4G=0._r8
  allocate(DOXYG(3,JZ,JY,JX));  DOXYG=0._r8
  allocate(DZ2GG(3,JZ,JY,JX));  DZ2GG=0._r8
  allocate(DZ2OG(3,JZ,JY,JX));  DZ2OG=0._r8
  allocate(DNH3G(3,JZ,JY,JX));  DNH3G=0._r8
  allocate(VLWatMicPCO(0:JZ,JY,JX)); VLWatMicPCO=0._r8
  allocate(VLWatMicPCH(0:JZ,JY,JX)); VLWatMicPCH=0._r8
  allocate(VLWatMicPOX(0:JZ,JY,JX)); VLWatMicPOX=0._r8
  allocate(VLWatMicPNG(0:JZ,JY,JX)); VLWatMicPNG=0._r8
  allocate(VLWatMicPN2(0:JZ,JY,JX)); VLWatMicPN2=0._r8
  allocate(VLWatMicPN3(0:JZ,JY,JX)); VLWatMicPN3=0._r8
  allocate(VLWatMicPNB(0:JZ,JY,JX)); VLWatMicPNB=0._r8
  allocate(VLWatMicPHG(0:JZ,JY,JX)); VLWatMicPHG=0._r8
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
  allocate(PARG_cef(idg_beg:idg_end-1,JY,JX));      PARG_cef=0._r8

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
  allocate(trcg_RQR0(idg_beg:idg_end,JV,JH));     trcg_RQR0=0._r8
  allocate(trcn_RQR0(ids_nut_beg:ids_nuts_end,JV,JH));     trcn_RQR0=0._r8
  allocate(OQNH2(1:jcplx,JZ,JY,JX));OQNH2=0._r8
  allocate(OQPH2(1:jcplx,JZ,JY,JX));OQPH2=0._r8
  allocate(OQAH2(1:jcplx,JZ,JY,JX));OQAH2=0._r8
  allocate(TOCFHS(1:jcplx,JZ,JY,JX));TOCFHS=0._r8
  allocate(TONFHS(1:jcplx,JZ,JY,JX));TONFHS=0._r8
  allocate(TOPFHS(1:jcplx,JZ,JY,JX));TOPFHS=0._r8
  allocate(TOAFHS(1:jcplx,JZ,JY,JX));TOAFHS=0._r8
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

  allocate(trcg_solsml2(idg_beg:idg_end-1,JS,JY,JX));trcg_solsml2=0._r8
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

  allocate(trcg_TBLS(idg_beg:idg_end-1,JS,JY,JX));   trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_TBLS=0._r8

  allocate(RQROP(1:jcplx,2,2,JV,JH));RQROP=0._r8
  allocate(RQROA(1:jcplx,2,2,JV,JH));RQROA=0._r8
  allocate(trcg_RQR(idg_beg:idg_end-1,2,2,JV,JH));  trcg_RQR=0._r8
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
  allocate(RDFR_gas(idg_beg:idg_end-1,JY,JX));      RDFR_gas=0._r8
  allocate(R1BSK2(JZ,JY,JX));   R1BSK2=0._r8
  allocate(RQROC(1:jcplx,2,2,JV,JH));RQROC=0._r8
  allocate(RQRON(1:jcplx,2,2,JV,JH));RQRON=0._r8

  allocate(ZNH4H2(JZ,JY,JX));   ZNH4H2=0._r8
  allocate(ZN4BH2(JZ,JY,JX));   ZN4BH2=0._r8
  allocate(ZNH3H2(JZ,JY,JX));   ZNH3H2=0._r8
  allocate(ZN3BH2(JZ,JY,JX));   ZN3BH2=0._r8
  allocate(ZNO3H2(JZ,JY,JX));   ZNO3H2=0._r8
  allocate(ZNOBH2(JZ,JY,JX));   ZNOBH2=0._r8
  allocate(H2P4H2(JZ,JY,JX));   H2P4H2=0._r8
  allocate(H2PBH2(JZ,JY,JX));   H2PBH2=0._r8
  allocate(ZNO2H2(JZ,JY,JX));   ZNO2H2=0._r8
  allocate(OQCH2(1:jcplx,JZ,JY,JX));OQCH2=0._r8

  allocate(trc_gasml2(idg_beg:idg_end,0:JZ,JY,JX)); trc_gasml2(:,:,:,:)=0._r8
  allocate(trc_solml2(ids_beg:ids_end,0:JZ,JY,JX)); trc_solml2(:,:,:,:)=0._r8
  allocate(trc_soHml2(ids_beg:ids_end,0:JZ,JY,JX)); trc_soHml2(:,:,:,:)=0._r8

  allocate(RQSH2P(2,JV,JH));    RQSH2P=0._r8
  allocate(RQSH1P(2,JV,JH));    RQSH1P=0._r8
  allocate(trcg_RQS(idg_beg:idg_end-1,2,JV,JH));    trcg_RQS=0._r8
  allocate(trcn_RQR(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_RQR=0._r8
  allocate(R3GasADFlx(idg_beg:idg_end-1,3,JD,JV,JH));R3GasADFlx=0._r8
  allocate(RTGasADFlx(idg_beg:idg_end-1,JZ,JY,JH));RTGasADFlx=0._r8

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

  allocate(RGasDSFlx(idg_beg:idg_end,0:JZ,JY,JX)); RGasDSFlx=0._r8
  allocate(RGasSSVol(idg_beg:idg_end,JY,JX)); RGasSSVol=0._r8
  allocate(RHGFHS(3,JD,JV,JH)); RHGFHS=0._r8
  allocate(RCHFHS(3,JD,JV,JH)); RCHFHS=0._r8
  allocate(ROXFHS(3,JD,JV,JH)); ROXFHS=0._r8
  allocate(RNGFHS(3,JD,JV,JH)); RNGFHS=0._r8
  allocate(RQSNH4(2,JV,JH));    RQSNH4=0._r8
  allocate(RQSNO3(2,JV,JH));    RQSNO3=0._r8
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
  end subroutine InitTransfrData

!----------------------------------------------------------------------
  subroutine DestructTransfrData
  use abortutils, only : destroy
  implicit none

  call destroy(OQC2)
  call destroy(OQN2)
  call destroy(OQP2)
  call destroy(OQA2)
  call destroy(ROCSK2)
  call destroy(RONSK2)
  call destroy(ROPSK2)
  call destroy(ROASK2)
  call destroy(ROCFLS)
  call destroy(RONFLS)
  call destroy(ROPFLS)
  call destroy(ROAFLS)
  call destroy(ROCFHS)
  call destroy(RONFHS)
  call destroy(ROPFHS)
  call destroy(ROAFHS)
  call destroy(TOCFLS)
  call destroy(TONFLS)
  call destroy(TOPFLS)
  call destroy(TOAFLS)

  call destroy(RQROC0)
  call destroy(RQRON0)
  call destroy(RQROA0)
  call destroy(RQROP0)
  call destroy(ROCFL0)
  call destroy(RONFL0)
  call destroy(ROPFL0)
  call destroy(ROAFL0)
  call destroy(ROCFL1)
  call destroy(RONFL1)
  call destroy(ROPFL1)
  call destroy(ROAFL1)

  call destroy(trcs_RFL1)

  call destroy(ROCFXS)
  call destroy(RONFXS)
  call destroy(ROPFXS)
  call destroy(ROAFXS)
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
  call destroy(TQROC)
  call destroy(TQRON)
  call destroy(TQROP)
  call destroy(TQROA)
  call destroy(trcg_TQR)
  call destroy(trcn_TQR)
  call destroy(trcg_TQ)
  call destroy(trcn_TQ)

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
  call destroy(VLWatMicPCO)
  call destroy(VLWatMicPCH)
  call destroy(VLWatMicPOX)
  call destroy(VLWatMicPNG)
  call destroy(VLWatMicPN2)
  call destroy(VLWatMicPN3)
  call destroy(VLWatMicPNB)
  call destroy(VLWatMicPHG)
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
  call destroy(OQNH2)
  call destroy(OQPH2)
  call destroy(OQAH2)
  call destroy(TOCFHS)
  call destroy(TONFHS)
  call destroy(TOPFHS)
  call destroy(TOAFHS)
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
  call destroy(trcg_RQS)

  call destroy(RQROP)
  call destroy(RQROA)
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
  call destroy(RQROC)
  call destroy(RQRON)

  call destroy(ZNH4H2)
  call destroy(ZN4BH2)
  call destroy(ZNH3H2)
  call destroy(ZN3BH2)
  call destroy(ZNO3H2)
  call destroy(ZNOBH2)
  call destroy(H2P4H2)
  call destroy(H2PBH2)
  call destroy(ZNO2H2)
  call destroy(OQCH2)

  call destroy(trc_gasml2)
  call destroy(trc_solml2)
  call destroy(trc_soHml2)

  call destroy(RQSH2P)
  call destroy(RQSH1P)

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
  call destroy(RQSNH4)
  call destroy(RQSNO3)
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
  end subroutine DestructTransfrData

end module TranspNoSaltDataMod
