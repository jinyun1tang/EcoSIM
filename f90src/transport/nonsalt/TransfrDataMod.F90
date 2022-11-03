module TransfrDataMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc
implicit none
  CHARACTER(LEN=*), Private,PARAMETER :: MOD_FILENAME=__FILE__
  real(r8), parameter :: XFRS=0.05_r8
  real(r8), Parameter :: VFLWX=0.5

  real(r8),allocatable :: RFLOC(:)
  real(r8),allocatable :: RFLON(:)
  real(r8),allocatable :: RFLOP(:)
  real(r8),allocatable :: RFLOA(:)
  real(r8),allocatable :: DFVOC(:)
  real(r8),allocatable :: DFVON(:)
  real(r8),allocatable :: DFVOP(:)
  real(r8),allocatable :: DFVOA(:)
  real(r8),allocatable :: DFHOC(:)
  real(r8),allocatable :: DFHON(:)
  real(r8),allocatable :: DFHOP(:)
  real(r8),allocatable :: DFHOA(:)
  real(r8),allocatable :: COQC1(:)
  real(r8),allocatable :: COQC2(:)
  real(r8),allocatable :: COQN1(:)
  real(r8),allocatable :: COQN2(:)
  real(r8),allocatable :: COQP1(:)
  real(r8),allocatable :: COQP2(:)
  real(r8),allocatable :: COQA1(:)
  real(r8),allocatable :: COQA2(:)
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

  real(r8),allocatable ::  RCOBLS(:,:,:)                      !
  real(r8),allocatable ::  RCHBLS(:,:,:)                      !
  real(r8),allocatable ::  ROXBLS(:,:,:)                      !
  real(r8),allocatable ::  RNGBLS(:,:,:)                      !
  real(r8),allocatable ::  RN2BLS(:,:,:)                      !
  real(r8),allocatable ::  RN4BLW(:,:,:)                      !
  real(r8),allocatable ::  RN3BLW(:,:,:)                      !
  real(r8),allocatable ::  RNOBLW(:,:,:)                      !
  real(r8),allocatable ::  RH1PBS(:,:,:)                      !
  real(r8),allocatable ::  RH2PBS(:,:,:)                      !
  real(r8),allocatable ::  ROCFL0(:,:,:)                      !
  real(r8),allocatable ::  RONFL0(:,:,:)                      !
  real(r8),allocatable ::  ROPFL0(:,:,:)                      !
  real(r8),allocatable ::  ROAFL0(:,:,:)                      !
  real(r8),allocatable ::  ROCFL1(:,:,:)                      !
  real(r8),allocatable ::  RONFL1(:,:,:)                      !
  real(r8),allocatable ::  ROPFL1(:,:,:)                      !
  real(r8),allocatable ::  ROAFL1(:,:,:)                      !
  real(r8),allocatable ::  RCOFL0(:,:)                        !
  real(r8),allocatable ::  RCHFL0(:,:)                        !
  real(r8),allocatable ::  ROXFL0(:,:)                        !
  real(r8),allocatable ::  RNGFL0(:,:)                        !
  real(r8),allocatable ::  RN2FL0(:,:)                        !
  real(r8),allocatable ::  RHGFL0(:,:)                        !
  real(r8),allocatable ::  RN4FL0(:,:)                        !
  real(r8),allocatable ::  RN3FL0(:,:)                        !
  real(r8),allocatable ::  RNOFL0(:,:)                        !
  real(r8),allocatable ::  RNXFL0(:,:)                        !
  real(r8),allocatable ::  RH2PF0(:,:)                        !
  real(r8),allocatable ::  RCOFL1(:,:)                        !
  real(r8),allocatable ::  RCHFL1(:,:)                        !
  real(r8),allocatable ::  ROXFL1(:,:)                        !
  real(r8),allocatable ::  RNGFL1(:,:)                        !
  real(r8),allocatable ::  RN2FL1(:,:)                        !
  real(r8),allocatable ::  RHGFL1(:,:)                        !
  real(r8),allocatable ::  RN4FL1(:,:)                        !
  real(r8),allocatable ::  RN3FL1(:,:)                        !
  real(r8),allocatable ::  RNOFL1(:,:)                        !
  real(r8),allocatable ::  RNXFL1(:,:)                        !
  real(r8),allocatable ::  RH2PF1(:,:)                        !
  real(r8),allocatable ::  RN4FL2(:,:)                        !
  real(r8),allocatable ::  RN3FL2(:,:)                        !
  real(r8),allocatable ::  RNOFL2(:,:)                        !
  real(r8),allocatable ::  RNXFL2(:,:)                        !
  real(r8),allocatable ::  RH2BF2(:,:)                        !
  real(r8),allocatable ::  RH1PF0(:,:)                        !
  real(r8),allocatable ::  RH1PF1(:,:)                        !
  real(r8),allocatable ::  RH1BF2(:,:)                        !
  real(r8),allocatable ::  ROCFXS(:,:,:,:)                    !
  real(r8),allocatable ::  RONFXS(:,:,:,:)                    !
  real(r8),allocatable ::  ROPFXS(:,:,:,:)                    !
  real(r8),allocatable ::  ROAFXS(:,:,:,:)                    !
  real(r8),allocatable ::  RCOFXS(:,:,:)                      !
  real(r8),allocatable ::  RCHFXS(:,:,:)                      !
  real(r8),allocatable ::  ROXFXS(:,:,:)                      !
  real(r8),allocatable ::  RNGFXS(:,:,:)                      !
  real(r8),allocatable ::  RN2FXS(:,:,:)                      !
  real(r8),allocatable ::  RN4FXW(:,:,:)                      !
  real(r8),allocatable ::  RN3FXW(:,:,:)                      !
  real(r8),allocatable ::  RNOFXW(:,:,:)                      !
  real(r8),allocatable ::  RH2PXS(:,:,:)                      !
  real(r8),allocatable ::  RN4FXB(:,:,:)                      !
  real(r8),allocatable ::  RN3FXB(:,:,:)                      !
  real(r8),allocatable ::  RNOFXB(:,:,:)                      !
  real(r8),allocatable ::  RH2BXB(:,:,:)                      !
  real(r8),allocatable ::  RNXFXS(:,:,:)                      !
  real(r8),allocatable ::  RNXFXB(:,:,:)                      !
  real(r8),allocatable ::  RH1PXS(:,:,:)                      !
  real(r8),allocatable ::  RH1BXB(:,:,:)                      !
  real(r8),allocatable ::  CLSGL2(:,:,:)                      !
  real(r8),allocatable ::  CQSGL2(:,:,:)                      !
  real(r8),allocatable ::  POSGL2(:,:,:)                      !
  real(r8),allocatable ::  OLSGL2(:,:,:)                      !
  real(r8),allocatable ::  ZNSGL2(:,:,:)                      !
  real(r8),allocatable ::  ZLSGL2(:,:,:)                      !
  real(r8),allocatable ::  ZVSGL2(:,:,:)                      !
  real(r8),allocatable ::  HLSGL2(:,:,:)                      !
  real(r8),allocatable ::  ZOSGL2(:,:,:)                      !
  real(r8),allocatable ::  RCODFG(:,:,:)                      !CO2 volatilization/dissolution within soil layer	g d-2 t-1
  real(r8),allocatable ::  RCHDFG(:,:,:)                      !CH4 volatilization/dissolution within soil layer	g d-2 t-1
  real(r8),allocatable ::  ROXDFG(:,:,:)                      !O2 volatilization/dissolution within soil layer	g d-2 t-1
  real(r8),allocatable ::  RNGDFG(:,:,:)                      !N2 volatilization/dissolution within soil layer	g d-2 t-1
  real(r8),allocatable ::  RN2DFG(:,:,:)                      !N2O volatilization/dissolution within soil layer	g d-2 t-1
  real(r8),allocatable ::  RN3DFG(:,:,:)                      !NH3 volatilization/dissolution within soil layer non-band	g d-2 t-1
  real(r8),allocatable ::  RNBDFG(:,:,:)                      !NH3 volatilization/dissolution within soil layer band	g d-2 t-1
  real(r8),allocatable ::  TQROC(:,:,:)                       !
  real(r8),allocatable ::  TQRON(:,:,:)                       !
  real(r8),allocatable ::  TQROP(:,:,:)                       !
  real(r8),allocatable ::  TQROA(:,:,:)                       !
  real(r8),allocatable ::  TQRCHS(:,:)                        !
  real(r8),allocatable ::  TQROXS(:,:)                        !
  real(r8),allocatable ::  TQRNGS(:,:)                        !
  real(r8),allocatable ::  TQRN2S(:,:)                        !
  real(r8),allocatable ::  TQRNH4(:,:)                        !
  real(r8),allocatable ::  TQRNH3(:,:)                        !
  real(r8),allocatable ::  TQRNO3(:,:)                        !
  real(r8),allocatable ::  TQRH2P(:,:)                        !
  real(r8),allocatable ::  TQRNO2(:,:)                        !
  real(r8),allocatable ::  TQRHGS(:,:)                        !
  real(r8),allocatable ::  TQSCOS(:,:)                        !
  real(r8),allocatable ::  TQRCOS(:,:)                        !
  real(r8),allocatable ::  TQSCHS(:,:)                        !
  real(r8),allocatable ::  TQSOXS(:,:)                        !
  real(r8),allocatable ::  TQSNGS(:,:)                        !
  real(r8),allocatable ::  TQSN2S(:,:)                        !
  real(r8),allocatable ::  TQSNH4(:,:)                        !
  real(r8),allocatable ::  TQSNH3(:,:)                        !
  real(r8),allocatable ::  TQSNO3(:,:)                        !
  real(r8),allocatable ::  TQSH1P(:,:)                        !
  real(r8),allocatable ::  TQSH2P(:,:)                        !
  real(r8),allocatable ::  TCOFLS(:,:,:)                      !
  real(r8),allocatable ::  TCHFLS(:,:,:)                      !
  real(r8),allocatable ::  TOXFLS(:,:,:)                      !
  real(r8),allocatable ::  TNGFLS(:,:,:)                      !
  real(r8),allocatable ::  TN2FLS(:,:,:)                      !
  real(r8),allocatable ::  TN4FLW(:,:,:)                      !
  real(r8),allocatable ::  TN3FLW(:,:,:)                      !
  real(r8),allocatable ::  TNOFLW(:,:,:)                      !
  real(r8),allocatable ::  TH2PFS(:,:,:)                      !
  real(r8),allocatable ::  TN4FLB(:,:,:)                      !
  real(r8),allocatable ::  TN3FLB(:,:,:)                      !
  real(r8),allocatable ::  TNOFLB(:,:,:)                      !
  real(r8),allocatable ::  TH2BFB(:,:,:)                      !
  real(r8),allocatable ::  TNXFLS(:,:,:)                      !
  real(r8),allocatable ::  TCOFLG(:,:,:)                      !
  real(r8),allocatable ::  TCHFLG(:,:,:)                      !
  real(r8),allocatable ::  TOXFLG(:,:,:)                      !
  real(r8),allocatable ::  TNGFLG(:,:,:)                      !
  real(r8),allocatable ::  TN2FLG(:,:,:)                      !
  real(r8),allocatable ::  TQRH1P(:,:)                        !
  real(r8),allocatable ::  TH1PFS(:,:,:)                      !
  real(r8),allocatable ::  TH1BFB(:,:,:)                      !
  real(r8),allocatable ::  DCO2G(:,:,:,:)                     !
  real(r8),allocatable ::  DCH4G(:,:,:,:)                     !
  real(r8),allocatable ::  DOXYG(:,:,:,:)                     !
  real(r8),allocatable ::  DZ2GG(:,:,:,:)                     !
  real(r8),allocatable ::  DZ2OG(:,:,:,:)                     !
  real(r8),allocatable ::  DNH3G(:,:,:,:)                     !
  real(r8),allocatable ::  VOLWCO(:,:,:)                      !
  real(r8),allocatable ::  VOLWCH(:,:,:)                      !
  real(r8),allocatable ::  VOLWOX(:,:,:)                      !
  real(r8),allocatable ::  VOLWNG(:,:,:)                      !
  real(r8),allocatable ::  VOLWN2(:,:,:)                      !
  real(r8),allocatable ::  VOLWN3(:,:,:)                      !
  real(r8),allocatable ::  VOLWNB(:,:,:)                      !
  real(r8),allocatable ::  VOLWHG(:,:,:)                      !
  real(r8),allocatable ::  HGSGL2(:,:,:)                      !
  real(r8),allocatable ::  DH2GG(:,:,:,:)                     !
  real(r8),allocatable ::  RHGFXS(:,:,:)                      !
  real(r8),allocatable ::  THGFHS(:,:,:)                      !
  real(r8),allocatable ::  THGFLG(:,:,:)                      !
  real(r8),allocatable ::  FLVM(:,:,:)                        !
  real(r8),allocatable ::  THETH2(:,:,:)                      !
  real(r8),allocatable ::  THETHL(:,:,:)                      !
  real(r8),allocatable ::  VOLPMA(:,:,:)                      !
  real(r8),allocatable ::  VOLPMB(:,:,:)                      !
  real(r8),allocatable ::  VOLWMA(:,:,:)                      !
  real(r8),allocatable ::  VOLWMB(:,:,:)                      !
  real(r8),allocatable ::  VOLWXA(:,:,:)                      !
  real(r8),allocatable ::  VOLWXB(:,:,:)                      !
  real(r8),allocatable ::  PARGCO(:,:)                        !
  real(r8),allocatable ::  PARGCH(:,:)                        !
  real(r8),allocatable ::  PARGOX(:,:)                        !
  real(r8),allocatable ::  PARGNG(:,:)                        !
  real(r8),allocatable ::  PARGN2(:,:)                        !
  real(r8),allocatable ::  PARGN3(:,:)                        !
  real(r8),allocatable ::  PARGH2(:,:)                        !
  real(r8),allocatable ::  RCOSK2(:,:,:)                      !
  real(r8),allocatable ::  ROXSK2(:,:,:)                      !
  real(r8),allocatable ::  RCHSK2(:,:,:)                      !
  real(r8),allocatable ::  RNGSK2(:,:,:)                      !
  real(r8),allocatable ::  RN2SK2(:,:,:)                      !
  real(r8),allocatable ::  RN4SK2(:,:,:)                      !
  real(r8),allocatable ::  RN3SK2(:,:,:)                      !
  real(r8),allocatable ::  RNOSK2(:,:,:)                      !
  real(r8),allocatable ::  RHPSK2(:,:,:)                      !
  real(r8),allocatable ::  R4BSK2(:,:,:)                      !
  real(r8),allocatable ::  R3BSK2(:,:,:)                      !
  real(r8),allocatable ::  RNBSK2(:,:,:)                      !
  real(r8),allocatable ::  RHBSK2(:,:,:)                      !
  real(r8),allocatable ::  RNXSK2(:,:,:)                      !
  real(r8),allocatable ::  RNZSK2(:,:,:)                      !
  real(r8),allocatable ::  RHGSK2(:,:,:)                      !
  real(r8),allocatable ::  RNHSK2(:,:,:)                      !
  real(r8),allocatable ::  R1PSK2(:,:,:)                      !
  real(r8),allocatable ::  TN3FLG(:,:,:)                      !
  real(r8),allocatable ::  RCOBBL(:,:,:)                      !
  real(r8),allocatable ::  RCHBBL(:,:,:)                      !
  real(r8),allocatable ::  ROXBBL(:,:,:)                      !
  real(r8),allocatable ::  RNGBBL(:,:,:)                      !
  real(r8),allocatable ::  RN2BBL(:,:,:)                      !
  real(r8),allocatable ::  RN3BBL(:,:,:)                      !
  real(r8),allocatable ::  RNBBBL(:,:,:)                      !
  real(r8),allocatable ::  RHGBBL(:,:,:)                      !
  real(r8),allocatable ::  RQRCOS0(:,:)                       !
  real(r8),allocatable ::  RQRCHS0(:,:)                       !
  real(r8),allocatable ::  RQROXS0(:,:)                       !
  real(r8),allocatable ::  RQRNGS0(:,:)                       !
  real(r8),allocatable ::  RQRN2S0(:,:)                       !
  real(r8),allocatable ::  RQRHGS0(:,:)                       !
  real(r8),allocatable ::  RQRNH40(:,:)                       !
  real(r8),allocatable ::  RQRNH30(:,:)                       !
  real(r8),allocatable ::  RQRNO30(:,:)                       !
  real(r8),allocatable ::  RQRNO20(:,:)                       !
  real(r8),allocatable ::  RQRH2P0(:,:)                       !
  real(r8),allocatable ::  RQRH1P0(:,:)                       !
  real(r8),allocatable ::  OQNH2(:,:,:,:)                     !
  real(r8),allocatable ::  OQPH2(:,:,:,:)                     !
  real(r8),allocatable ::  OQAH2(:,:,:,:)                     !
  real(r8),allocatable ::  TOCFHS(:,:,:,:)                    !
  real(r8),allocatable ::  TONFHS(:,:,:,:)                    !
  real(r8),allocatable ::  TOPFHS(:,:,:,:)                    !
  real(r8),allocatable ::  TOAFHS(:,:,:,:)                    !
  real(r8),allocatable ::  TCOFHS(:,:,:)                      !
  real(r8),allocatable ::  TCHFHS(:,:,:)                      !
  real(r8),allocatable ::  TOXFHS(:,:,:)                      !
  real(r8),allocatable ::  TNGFHS(:,:,:)                      !
  real(r8),allocatable ::  TN2FHS(:,:,:)                      !
  real(r8),allocatable ::  TN4FHW(:,:,:)                      !
  real(r8),allocatable ::  TN3FHW(:,:,:)                      !
  real(r8),allocatable ::  TNOFHW(:,:,:)                      !
  real(r8),allocatable ::  TH2PHS(:,:,:)                      !
  real(r8),allocatable ::  TN4FHB(:,:,:)                      !
  real(r8),allocatable ::  TN3FHB(:,:,:)                      !
  real(r8),allocatable ::  TNOFHB(:,:,:)                      !
  real(r8),allocatable ::  TH2BHB(:,:,:)                      !
  real(r8),allocatable ::  TNXFHS(:,:,:)                      !
  real(r8),allocatable ::  ZNO2B2(:,:,:)                      !
  real(r8),allocatable ::  ZN2BH2(:,:,:)                      !
  real(r8),allocatable ::  TNXFLB(:,:,:)                      !
  real(r8),allocatable ::  TNXFHB(:,:,:)                      !
  real(r8),allocatable ::  H1P4H2(:,:,:)                      !
  real(r8),allocatable ::  H1PBH2(:,:,:)                      !
  real(r8),allocatable ::  TH1PHS(:,:,:)                      !
  real(r8),allocatable ::  TH1BHB(:,:,:)                      !
  REAL(R8),allocatable ::  H1PO42(:,:,:)                      !
  real(r8),allocatable ::  H1POB2(:,:,:)                      !
  real(r8),allocatable ::  OCSGL2(:,:,:)                      !
  real(r8),allocatable ::  ONSGL2(:,:,:)                      !
  real(r8),allocatable ::  OPSGL2(:,:,:)                      !
  real(r8),allocatable ::  OASGL2(:,:,:)                      !
  real(r8),allocatable ::  CO2W2(:,:,:)                       !
  real(r8),allocatable ::  CH4W2(:,:,:)                       !
  real(r8),allocatable ::  OXYW2(:,:,:)                       !
  real(r8),allocatable ::  ZNGW2(:,:,:)                       !
  real(r8),allocatable ::  ZN2W2(:,:,:)                       !
  real(r8),allocatable ::  ZN4W2(:,:,:)                       !
  real(r8),allocatable ::  ZN3W2(:,:,:)                       !
  real(r8),allocatable ::  ZNOW2(:,:,:)                       !
  real(r8),allocatable ::  ZHPW2(:,:,:)                       !
  real(r8),allocatable ::  Z1PW2(:,:,:)                       !
  real(r8),allocatable ::  TCOBLS(:,:,:)                      !
  real(r8),allocatable ::  TCHBLS(:,:,:)                      !
  real(r8),allocatable ::  TOXBLS(:,:,:)                      !
  real(r8),allocatable ::  TNGBLS(:,:,:)                      !
  real(r8),allocatable ::  TN2BLS(:,:,:)                      !
  real(r8),allocatable ::  TN4BLW(:,:,:)                      !
  real(r8),allocatable ::  TN3BLW(:,:,:)                      !
  real(r8),allocatable ::  TNOBLW(:,:,:)                      !
  real(r8),allocatable ::  TH1PBS(:,:,:)                      !
  real(r8),allocatable ::  TH2PBS(:,:,:)                      !
  real(r8),allocatable ::  RQROP(:,:,:,:,:)                   !
  real(r8),allocatable ::  RQROA(:,:,:,:,:)                   !
  real(r8),allocatable ::  RQRCOS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRCHS(:,:,:,:)                    !
  real(r8),allocatable ::  RQROXS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNGS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRN2S(:,:,:,:)                    !
  real(r8),allocatable ::  RQRHGS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNH4(:,:,:,:)                    !
  real(r8),allocatable ::  RCODFS(:,:)                        !
  real(r8),allocatable ::  RCHDFS(:,:)                        !
  real(r8),allocatable ::  ROXDFS(:,:)                        !
  real(r8),allocatable ::  RNGDFS(:,:)                        !
  real(r8),allocatable ::  RN2DFS(:,:)                        !
  real(r8),allocatable ::  RN3DFS(:,:)                        !
  real(r8),allocatable ::  RNBDFS(:,:)                        !
  real(r8),allocatable ::  RHGDFS(:,:)                        !
  real(r8),allocatable ::  RCODFR(:,:)                        !
  real(r8),allocatable ::  RCHDFR(:,:)                        !
  real(r8),allocatable ::  ROXDFR(:,:)                        !
  real(r8),allocatable ::  RNGDFR(:,:)                        !
  real(r8),allocatable ::  RN2DFR(:,:)                        !
  real(r8),allocatable ::  RN3DFR(:,:)                        !
  real(r8),allocatable ::  RHGDFR(:,:)                        !
  real(r8),allocatable ::  R1BSK2(:,:,:)                      !
  real(r8),allocatable ::  RQROC(:,:,:,:,:)                   !
  real(r8),allocatable ::  RQRON(:,:,:,:,:)                   !
  real(r8),allocatable ::  ZNH4S2(:,:,:)                      !
  real(r8),allocatable ::  ZNH4B2(:,:,:)                      !
  real(r8),allocatable ::  ZNH3S2(:,:,:)                      !
  real(r8),allocatable ::  ZNH3B2(:,:,:)                      !
  real(r8),allocatable ::  ZNO3S2(:,:,:)                      !
  real(r8),allocatable ::  ZNO3B2(:,:,:)                      !
  real(r8),allocatable ::  H2PO42(:,:,:)                      !
  real(r8),allocatable ::  H2POB2(:,:,:)                      !
  real(r8),allocatable ::  ZNO2S2(:,:,:)                      !
  real(r8),allocatable ::  CGSGL2(:,:,:)                      !
  real(r8),allocatable ::  CHSGL2(:,:,:)                      !
  real(r8),allocatable ::  OGSGL2(:,:,:)                      !
  real(r8),allocatable ::  ZGSGL2(:,:,:)                      !
  real(r8),allocatable ::  Z2SGL2(:,:,:)                      !
  real(r8),allocatable ::  ZHSGL2(:,:,:)                      !
  real(r8),allocatable ::  CO2SH2(:,:,:)                      !
  real(r8),allocatable ::  CH4SH2(:,:,:)                      !
  real(r8),allocatable ::  OXYSH2(:,:,:)                      !
  real(r8),allocatable ::  Z2GSH2(:,:,:)                      !
  real(r8),allocatable ::  Z2OSH2(:,:,:)                      !
  real(r8),allocatable ::  ZNH4H2(:,:,:)                      !
  real(r8),allocatable ::  ZN4BH2(:,:,:)                      !
  real(r8),allocatable ::  ZNH3H2(:,:,:)                      !
  real(r8),allocatable ::  ZN3BH2(:,:,:)                      !
  real(r8),allocatable ::  ZNO3H2(:,:,:)                      !
  real(r8),allocatable ::  ZNOBH2(:,:,:)                      !
  real(r8),allocatable ::  H2P4H2(:,:,:)                      !
  real(r8),allocatable ::  H2PBH2(:,:,:)                      !
  real(r8),allocatable ::  ZNO2H2(:,:,:)                      !
  real(r8),allocatable ::  OQCH2(:,:,:,:)                     !
  real(r8),allocatable ::  CH4G2(:,:,:)                       !
  real(r8),allocatable ::  CH4S2(:,:,:)                       !
  real(r8),allocatable ::  OXYG2(:,:,:)                       !
  real(r8),allocatable ::  OXYS2(:,:,:)                       !
  real(r8),allocatable ::  Z2GG2(:,:,:)                       !
  real(r8),allocatable ::  Z2GS2(:,:,:)                       !
  real(r8),allocatable ::  Z2OG2(:,:,:)                       !
  real(r8),allocatable ::  Z2OS2(:,:,:)                       !
  real(r8),allocatable ::  ZN3G2(:,:,:)                       !
  real(r8),allocatable ::  CO2G2(:,:,:)                       !
  real(r8),allocatable ::  CO2S2(:,:,:)                       !
  real(r8),allocatable ::  H2GG2(:,:,:)                       !
  real(r8),allocatable ::  H2GS2(:,:,:)                       !
  real(r8),allocatable ::  H2GSH2(:,:,:)                      !
  real(r8),allocatable ::  RQSH2P(:,:,:)                      !
  real(r8),allocatable ::  RQSH1P(:,:,:)                      !
  real(r8),allocatable ::  RQSCOS(:,:,:)                      !
  real(r8),allocatable ::  RQSCHS(:,:,:)                      !
  real(r8),allocatable ::  RQSOXS(:,:,:)                      !
  real(r8),allocatable ::  RQRH2P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRH1P(:,:,:,:)                    !
  real(r8),allocatable ::  RHGFLG(:,:,:,:)                    !
  real(r8),allocatable ::  THGFLS(:,:,:)                      !
  real(r8),allocatable ::  ROXFLG(:,:,:,:)                    !
  real(r8),allocatable ::  RN3FLG(:,:,:,:)                    !
  real(r8),allocatable ::  RCOFLG(:,:,:,:)                    !
  real(r8),allocatable ::  RNXFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RNXFHB(:,:,:,:)                    !
  real(r8),allocatable ::  RH1PFS(:,:,:,:)                    !
  real(r8),allocatable ::  RH1BFB(:,:,:,:)                    !
  real(r8),allocatable ::  RH1PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RH1BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RCHFLG(:,:,:,:)                    !
  real(r8),allocatable ::  RNGFLG(:,:,:,:)                    !
  real(r8),allocatable ::  RN2FLG(:,:,:,:)                    !
  real(r8),allocatable ::  RN3FHB(:,:,:,:)                    !
  real(r8),allocatable ::  RNOFHB(:,:,:,:)                    !
  real(r8),allocatable ::  RH2BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RNOFHW(:,:,:,:)                    !
  real(r8),allocatable ::  RH2PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RN4FHB(:,:,:,:)                    !
  real(r8),allocatable ::  RN2FHS(:,:,:,:)                    !
  real(r8),allocatable ::  RN4FHW(:,:,:,:)                    !
  real(r8),allocatable ::  RN3FHW(:,:,:,:)                    !
  real(r8),allocatable ::  RHGDFG(:,:,:)                      !
  real(r8),allocatable ::  RHGFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCHFHS(:,:,:,:)                    !
  real(r8),allocatable ::  ROXFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RNGFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNH3(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNO3(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNO2(:,:,:,:)                    !
  real(r8),allocatable ::  RQSNGS(:,:,:)                      !
  real(r8),allocatable ::  RQSN2S(:,:,:)                      !
  real(r8),allocatable ::  RQSNH4(:,:,:)                      !
  real(r8),allocatable ::  RQSNH3(:,:,:)                      !
  real(r8),allocatable ::  RQSNO3(:,:,:)                      !
  real(r8),allocatable ::  RCOFLS(:,:,:,:)                    !
  real(r8),allocatable ::  RCHFLS(:,:,:,:)                    !
  real(r8),allocatable ::  ROXFLS(:,:,:,:)                    !
  real(r8),allocatable ::  RNGFLS(:,:,:,:)                    !
  real(r8),allocatable ::  RN2FLS(:,:,:,:)                    !
  real(r8),allocatable ::  RHGFLS(:,:,:,:)                    !
  real(r8),allocatable ::  RN4FLW(:,:,:,:)                    !
  real(r8),allocatable ::  RN3FLW(:,:,:,:)                    !
  real(r8),allocatable ::  RNOFLW(:,:,:,:)                    !
  real(r8),allocatable ::  RNXFLS(:,:,:,:)                    !
  real(r8),allocatable ::  RH2PFS(:,:,:,:)                    !
  real(r8),allocatable ::  RN4FLB(:,:,:,:)                    !
  real(r8),allocatable ::  RN3FLB(:,:,:,:)                    !
  real(r8),allocatable ::  RNOFLB(:,:,:,:)                    !
  real(r8),allocatable ::  RNXFLB(:,:,:,:)                    !
  real(r8),allocatable ::  RH2BFB(:,:,:,:)                    !
  real(r8),allocatable ::  RCOFHS(:,:,:,:)                    !
!----------------------------------------------------------------------

contains
  subroutine InitTransfrData
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer :: n_litrsfk
  n_litrsfk=micpar%n_litrsfk
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

  allocate(RCOBLS(JS,JY,JX));   RCOBLS=0._r8
  allocate(RCHBLS(JS,JY,JX));   RCHBLS=0._r8
  allocate(ROXBLS(JS,JY,JX));   ROXBLS=0._r8
  allocate(RNGBLS(JS,JY,JX));   RNGBLS=0._r8
  allocate(RN2BLS(JS,JY,JX));   RN2BLS=0._r8
  allocate(RN4BLW(JS,JY,JX));   RN4BLW=0._r8
  allocate(RN3BLW(JS,JY,JX));   RN3BLW=0._r8
  allocate(RNOBLW(JS,JY,JX));   RNOBLW=0._r8
  allocate(RH1PBS(JS,JY,JX));   RH1PBS=0._r8
  allocate(RH2PBS(JS,JY,JX));   RH2PBS=0._r8
  allocate(ROCFL0(1:n_litrsfk,JY,JX));  ROCFL0=0._r8
  allocate(RONFL0(1:n_litrsfk,JY,JX));  RONFL0=0._r8
  allocate(ROPFL0(1:n_litrsfk,JY,JX));  ROPFL0=0._r8
  allocate(ROAFL0(1:n_litrsfk,JY,JX));  ROAFL0=0._r8
  allocate(ROCFL1(1:n_litrsfk,JY,JX));  ROCFL1=0._r8
  allocate(RONFL1(1:n_litrsfk,JY,JX));  RONFL1=0._r8
  allocate(ROPFL1(1:n_litrsfk,JY,JX));  ROPFL1=0._r8
  allocate(ROAFL1(1:n_litrsfk,JY,JX));  ROAFL1=0._r8
  allocate(RCOFL0(JY,JX));      RCOFL0=0._r8
  allocate(RCHFL0(JY,JX));      RCHFL0=0._r8
  allocate(ROXFL0(JY,JX));      ROXFL0=0._r8
  allocate(RNGFL0(JY,JX));      RNGFL0=0._r8
  allocate(RN2FL0(JY,JX));      RN2FL0=0._r8
  allocate(RHGFL0(JY,JX));      RHGFL0=0._r8
  allocate(RN4FL0(JY,JX));      RN4FL0=0._r8
  allocate(RN3FL0(JY,JX));      RN3FL0=0._r8
  allocate(RNOFL0(JY,JX));      RNOFL0=0._r8
  allocate(RNXFL0(JY,JX));      RNXFL0=0._r8
  allocate(RH2PF0(JY,JX));      RH2PF0=0._r8
  allocate(RCOFL1(JY,JX));      RCOFL1=0._r8
  allocate(RCHFL1(JY,JX));      RCHFL1=0._r8
  allocate(ROXFL1(JY,JX));      ROXFL1=0._r8
  allocate(RNGFL1(JY,JX));      RNGFL1=0._r8
  allocate(RN2FL1(JY,JX));      RN2FL1=0._r8
  allocate(RHGFL1(JY,JX));      RHGFL1=0._r8
  allocate(RN4FL1(JY,JX));      RN4FL1=0._r8
  allocate(RN3FL1(JY,JX));      RN3FL1=0._r8
  allocate(RNOFL1(JY,JX));      RNOFL1=0._r8
  allocate(RNXFL1(JY,JX));      RNXFL1=0._r8
  allocate(RH2PF1(JY,JX));      RH2PF1=0._r8
  allocate(RN4FL2(JY,JX));      RN4FL2=0._r8
  allocate(RN3FL2(JY,JX));      RN3FL2=0._r8
  allocate(RNOFL2(JY,JX));      RNOFL2=0._r8
  allocate(RNXFL2(JY,JX));      RNXFL2=0._r8
  allocate(RH2BF2(JY,JX));      RH2BF2=0._r8
  allocate(RH1PF0(JY,JX));      RH1PF0=0._r8
  allocate(RH1PF1(JY,JX));      RH1PF1=0._r8
  allocate(RH1BF2(JY,JX));      RH1BF2=0._r8
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
  allocate(RCODFG(0:JZ,JY,JX)); RCODFG=0._r8
  allocate(RCHDFG(0:JZ,JY,JX)); RCHDFG=0._r8
  allocate(ROXDFG(0:JZ,JY,JX)); ROXDFG=0._r8
  allocate(RNGDFG(0:JZ,JY,JX)); RNGDFG=0._r8
  allocate(RN2DFG(0:JZ,JY,JX)); RN2DFG=0._r8
  allocate(RN3DFG(0:JZ,JY,JX)); RN3DFG=0._r8
  allocate(RNBDFG(0:JZ,JY,JX)); RNBDFG=0._r8
  allocate(TQROC(1:jcplx,JY,JX));   TQROC=0._r8
  allocate(TQRON(1:jcplx,JY,JX));   TQRON=0._r8
  allocate(TQROP(1:jcplx,JY,JX));   TQROP=0._r8
  allocate(TQROA(1:jcplx,JY,JX));   TQROA=0._r8
  allocate(TQRCHS(JY,JX));      TQRCHS=0._r8
  allocate(TQROXS(JY,JX));      TQROXS=0._r8
  allocate(TQRNGS(JY,JX));      TQRNGS=0._r8
  allocate(TQRN2S(JY,JX));      TQRN2S=0._r8
  allocate(TQRNH4(JY,JX));      TQRNH4=0._r8
  allocate(TQRNH3(JY,JX));      TQRNH3=0._r8
  allocate(TQRNO3(JY,JX));      TQRNO3=0._r8
  allocate(TQRH2P(JY,JX));      TQRH2P=0._r8
  allocate(TQRNO2(JY,JX));      TQRNO2=0._r8
  allocate(TQRHGS(JY,JX));      TQRHGS=0._r8
  allocate(TQSCOS(JY,JX));      TQSCOS=0._r8
  allocate(TQRCOS(JY,JX));      TQRCOS=0._r8
  allocate(TQSCHS(JY,JX));      TQSCHS=0._r8
  allocate(TQSOXS(JY,JX));      TQSOXS=0._r8
  allocate(TQSNGS(JY,JX));      TQSNGS=0._r8
  allocate(TQSN2S(JY,JX));      TQSN2S=0._r8
  allocate(TQSNH4(JY,JX));      TQSNH4=0._r8
  allocate(TQSNH3(JY,JX));      TQSNH3=0._r8
  allocate(TQSNO3(JY,JX));      TQSNO3=0._r8
  allocate(TQSH1P(JY,JX));      TQSH1P=0._r8
  allocate(TQSH2P(JY,JX));      TQSH2P=0._r8
  allocate(TCOFLS(JZ,JY,JX));   TCOFLS=0._r8
  allocate(TCHFLS(JZ,JY,JX));   TCHFLS=0._r8
  allocate(TOXFLS(JZ,JY,JX));   TOXFLS=0._r8
  allocate(TNGFLS(JZ,JY,JX));   TNGFLS=0._r8
  allocate(TN2FLS(JZ,JY,JX));   TN2FLS=0._r8
  allocate(TN4FLW(JZ,JY,JX));   TN4FLW=0._r8
  allocate(TN3FLW(JZ,JY,JX));   TN3FLW=0._r8
  allocate(TNOFLW(JZ,JY,JX));   TNOFLW=0._r8
  allocate(TH2PFS(JZ,JY,JX));   TH2PFS=0._r8
  allocate(TN4FLB(JZ,JY,JX));   TN4FLB=0._r8
  allocate(TN3FLB(JZ,JY,JX));   TN3FLB=0._r8
  allocate(TNOFLB(JZ,JY,JX));   TNOFLB=0._r8
  allocate(TH2BFB(JZ,JY,JX));   TH2BFB=0._r8
  allocate(TNXFLS(JZ,JY,JX));   TNXFLS=0._r8
  allocate(TCOFLG(JZ,JY,JX));   TCOFLG=0._r8
  allocate(TCHFLG(JZ,JY,JX));   TCHFLG=0._r8
  allocate(TOXFLG(JZ,JY,JX));   TOXFLG=0._r8
  allocate(TNGFLG(JZ,JY,JX));   TNGFLG=0._r8
  allocate(TN2FLG(JZ,JY,JX));   TN2FLG=0._r8
  allocate(TQRH1P(JY,JX));      TQRH1P=0._r8
  allocate(TH1PFS(JZ,JY,JX));   TH1PFS=0._r8
  allocate(TH1BFB(JZ,JY,JX));   TH1BFB=0._r8
  allocate(DCO2G(3,JZ,JY,JX));  DCO2G=0._r8
  allocate(DCH4G(3,JZ,JY,JX));  DCH4G=0._r8
  allocate(DOXYG(3,JZ,JY,JX));  DOXYG=0._r8
  allocate(DZ2GG(3,JZ,JY,JX));  DZ2GG=0._r8
  allocate(DZ2OG(3,JZ,JY,JX));  DZ2OG=0._r8
  allocate(DNH3G(3,JZ,JY,JX));  DNH3G=0._r8
  allocate(VOLWCO(0:JZ,JY,JX)); VOLWCO=0._r8
  allocate(VOLWCH(0:JZ,JY,JX)); VOLWCH=0._r8
  allocate(VOLWOX(0:JZ,JY,JX)); VOLWOX=0._r8
  allocate(VOLWNG(0:JZ,JY,JX)); VOLWNG=0._r8
  allocate(VOLWN2(0:JZ,JY,JX)); VOLWN2=0._r8
  allocate(VOLWN3(0:JZ,JY,JX)); VOLWN3=0._r8
  allocate(VOLWNB(0:JZ,JY,JX)); VOLWNB=0._r8
  allocate(VOLWHG(0:JZ,JY,JX)); VOLWHG=0._r8
  allocate(HGSGL2(JZ,JY,JX));   HGSGL2=0._r8
  allocate(DH2GG(3,JZ,JY,JX));  DH2GG=0._r8
  allocate(RHGFXS(JZ,JY,JX));   RHGFXS=0._r8
  allocate(THGFHS(JZ,JY,JX));   THGFHS=0._r8
  allocate(THGFLG(JZ,JY,JX));   THGFLG=0._r8
  allocate(FLVM(JZ,JY,JX));     FLVM=0._r8
  allocate(THETH2(JZ,JY,JX));   THETH2=0._r8
  allocate(THETHL(JZ,JY,JX));   THETHL=0._r8
  allocate(VOLPMA(JZ,JY,JX));   VOLPMA=0._r8
  allocate(VOLPMB(JZ,JY,JX));   VOLPMB=0._r8
  allocate(VOLWMA(JZ,JY,JX));   VOLWMA=0._r8
  allocate(VOLWMB(JZ,JY,JX));   VOLWMB=0._r8
  allocate(VOLWXA(0:JZ,JY,JX)); VOLWXA=0._r8
  allocate(VOLWXB(JZ,JY,JX));   VOLWXB=0._r8
  allocate(PARGCO(JY,JX));      PARGCO=0._r8
  allocate(PARGCH(JY,JX));      PARGCH=0._r8
  allocate(PARGOX(JY,JX));      PARGOX=0._r8
  allocate(PARGNG(JY,JX));      PARGNG=0._r8
  allocate(PARGN2(JY,JX));      PARGN2=0._r8
  allocate(PARGN3(JY,JX));      PARGN3=0._r8
  allocate(PARGH2(JY,JX));      PARGH2=0._r8
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
  allocate(RCOBBL(JZ,JY,JX));   RCOBBL=0._r8
  allocate(RCHBBL(JZ,JY,JX));   RCHBBL=0._r8
  allocate(ROXBBL(JZ,JY,JX));   ROXBBL=0._r8
  allocate(RNGBBL(JZ,JY,JX));   RNGBBL=0._r8
  allocate(RN2BBL(JZ,JY,JX));   RN2BBL=0._r8
  allocate(RN3BBL(JZ,JY,JX));   RN3BBL=0._r8
  allocate(RNBBBL(JZ,JY,JX));   RNBBBL=0._r8
  allocate(RHGBBL(JZ,JY,JX));   RHGBBL=0._r8
  allocate(RQRCOS0(JV,JH));     RQRCOS0=0._r8
  allocate(RQRCHS0(JV,JH));     RQRCHS0=0._r8
  allocate(RQROXS0(JV,JH));     RQROXS0=0._r8
  allocate(RQRNGS0(JV,JH));     RQRNGS0=0._r8
  allocate(RQRN2S0(JV,JH));     RQRN2S0=0._r8
  allocate(RQRHGS0(JV,JH));     RQRHGS0=0._r8
  allocate(RQRNH40(JV,JH));     RQRNH40=0._r8
  allocate(RQRNH30(JV,JH));     RQRNH30=0._r8
  allocate(RQRNO30(JV,JH));     RQRNO30=0._r8
  allocate(RQRNO20(JV,JH));     RQRNO20=0._r8
  allocate(RQRH2P0(JV,JH));     RQRH2P0=0._r8
  allocate(RQRH1P0(JV,JH));     RQRH1P0=0._r8
  allocate(OQNH2(1:jcplx,JZ,JY,JX));OQNH2=0._r8
  allocate(OQPH2(1:jcplx,JZ,JY,JX));OQPH2=0._r8
  allocate(OQAH2(1:jcplx,JZ,JY,JX));OQAH2=0._r8
  allocate(TOCFHS(1:jcplx,JZ,JY,JX));TOCFHS=0._r8
  allocate(TONFHS(1:jcplx,JZ,JY,JX));TONFHS=0._r8
  allocate(TOPFHS(1:jcplx,JZ,JY,JX));TOPFHS=0._r8
  allocate(TOAFHS(1:jcplx,JZ,JY,JX));TOAFHS=0._r8
  allocate(TCOFHS(JZ,JY,JX));   TCOFHS=0._r8
  allocate(TCHFHS(JZ,JY,JX));   TCHFHS=0._r8
  allocate(TOXFHS(JZ,JY,JX));   TOXFHS=0._r8
  allocate(TNGFHS(JZ,JY,JX));   TNGFHS=0._r8
  allocate(TN2FHS(JZ,JY,JX));   TN2FHS=0._r8
  allocate(TN4FHW(JZ,JY,JX));   TN4FHW=0._r8
  allocate(TN3FHW(JZ,JY,JX));   TN3FHW=0._r8
  allocate(TNOFHW(JZ,JY,JX));   TNOFHW=0._r8
  allocate(TH2PHS(JZ,JY,JX));   TH2PHS=0._r8
  allocate(TN4FHB(JZ,JY,JX));   TN4FHB=0._r8
  allocate(TN3FHB(JZ,JY,JX));   TN3FHB=0._r8
  allocate(TNOFHB(JZ,JY,JX));   TNOFHB=0._r8
  allocate(TH2BHB(JZ,JY,JX));   TH2BHB=0._r8
  allocate(TNXFHS(JZ,JY,JX));   TNXFHS=0._r8
  allocate(ZNO2B2(JZ,JY,JX));   ZNO2B2=0._r8
  allocate(ZN2BH2(JZ,JY,JX));   ZN2BH2=0._r8
  allocate(TNXFLB(JZ,JY,JX));   TNXFLB=0._r8
  allocate(TNXFHB(JZ,JY,JX));   TNXFHB=0._r8
  allocate(H1P4H2(JZ,JY,JX));   H1P4H2=0._r8
  allocate(H1PBH2(JZ,JY,JX));   H1PBH2=0._r8
  allocate(TH1PHS(JZ,JY,JX));   TH1PHS=0._r8
  allocate(TH1BHB(JZ,JY,JX));   TH1BHB=0._r8
  allocate(H1PO42(0:JZ,JY,JX)); H1PO42=0._r8
  allocate(H1POB2(0:JZ,JY,JX)); H1POB2=0._r8
  allocate(OCSGL2(0:JZ,JY,JX)); OCSGL2=0._r8
  allocate(ONSGL2(0:JZ,JY,JX)); ONSGL2=0._r8
  allocate(OPSGL2(0:JZ,JY,JX)); OPSGL2=0._r8
  allocate(OASGL2(0:JZ,JY,JX)); OASGL2=0._r8
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
  allocate(TCOBLS(JS,JY,JX));   TCOBLS=0._r8
  allocate(TCHBLS(JS,JY,JX));   TCHBLS=0._r8
  allocate(TOXBLS(JS,JY,JX));   TOXBLS=0._r8
  allocate(TNGBLS(JS,JY,JX));   TNGBLS=0._r8
  allocate(TN2BLS(JS,JY,JX));   TN2BLS=0._r8
  allocate(TN4BLW(JS,JY,JX));   TN4BLW=0._r8
  allocate(TN3BLW(JS,JY,JX));   TN3BLW=0._r8
  allocate(TNOBLW(JS,JY,JX));   TNOBLW=0._r8
  allocate(TH1PBS(JS,JY,JX));   TH1PBS=0._r8
  allocate(TH2PBS(JS,JY,JX));   TH2PBS=0._r8
  allocate(RQROP(1:jcplx,2,2,JV,JH));RQROP=0._r8
  allocate(RQROA(1:jcplx,2,2,JV,JH));RQROA=0._r8
  allocate(RQRCOS(2,2,JV,JH));  RQRCOS=0._r8
  allocate(RQRCHS(2,2,JV,JH));  RQRCHS=0._r8
  allocate(RQROXS(2,2,JV,JH));  RQROXS=0._r8
  allocate(RQRNGS(2,2,JV,JH));  RQRNGS=0._r8
  allocate(RQRN2S(2,2,JV,JH));  RQRN2S=0._r8
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
  allocate(RCODFR(JY,JX));      RCODFR=0._r8
  allocate(RCHDFR(JY,JX));      RCHDFR=0._r8
  allocate(ROXDFR(JY,JX));      ROXDFR=0._r8
  allocate(RNGDFR(JY,JX));      RNGDFR=0._r8
  allocate(RN2DFR(JY,JX));      RN2DFR=0._r8
  allocate(RN3DFR(JY,JX));      RN3DFR=0._r8
  allocate(RHGDFR(JY,JX));      RHGDFR=0._r8
  allocate(R1BSK2(JZ,JY,JX));   R1BSK2=0._r8
  allocate(RQROC(1:jcplx,2,2,JV,JH));RQROC=0._r8
  allocate(RQRON(1:jcplx,2,2,JV,JH));RQRON=0._r8
  allocate(ZNH4S2(0:JZ,JY,JX)); ZNH4S2=0._r8
  allocate(ZNH4B2(0:JZ,JY,JX)); ZNH4B2=0._r8
  allocate(ZNH3S2(0:JZ,JY,JX)); ZNH3S2=0._r8
  allocate(ZNH3B2(0:JZ,JY,JX)); ZNH3B2=0._r8
  allocate(ZNO3S2(0:JZ,JY,JX)); ZNO3S2=0._r8
  allocate(ZNO3B2(0:JZ,JY,JX)); ZNO3B2=0._r8
  allocate(H2PO42(0:JZ,JY,JX)); H2PO42=0._r8
  allocate(H2POB2(0:JZ,JY,JX)); H2POB2=0._r8
  allocate(ZNO2S2(0:JZ,JY,JX)); ZNO2S2=0._r8
  allocate(CGSGL2(JZ,JY,JX));   CGSGL2=0._r8
  allocate(CHSGL2(JZ,JY,JX));   CHSGL2=0._r8
  allocate(OGSGL2(JZ,JY,JX));   OGSGL2=0._r8
  allocate(ZGSGL2(JZ,JY,JX));   ZGSGL2=0._r8
  allocate(Z2SGL2(JZ,JY,JX));   Z2SGL2=0._r8
  allocate(ZHSGL2(JZ,JY,JX));   ZHSGL2=0._r8
  allocate(CO2SH2(JZ,JY,JX));   CO2SH2=0._r8
  allocate(CH4SH2(JZ,JY,JX));   CH4SH2=0._r8
  allocate(OXYSH2(JZ,JY,JX));   OXYSH2=0._r8
  allocate(Z2GSH2(JZ,JY,JX));   Z2GSH2=0._r8
  allocate(Z2OSH2(JZ,JY,JX));   Z2OSH2=0._r8
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
  allocate(CH4G2(JZ,JY,JX));    CH4G2=0._r8
  allocate(CH4S2(0:JZ,JY,JX));  CH4S2=0._r8
  allocate(OXYG2(JZ,JY,JX));    OXYG2=0._r8
  allocate(OXYS2(0:JZ,JY,JX));  OXYS2=0._r8
  allocate(Z2GG2(JZ,JY,JX));    Z2GG2=0._r8
  allocate(Z2GS2(0:JZ,JY,JX));  Z2GS2=0._r8
  allocate(Z2OG2(JZ,JY,JX));    Z2OG2=0._r8
  allocate(Z2OS2(0:JZ,JY,JX));  Z2OS2=0._r8
  allocate(ZN3G2(0:JZ,JY,JX));  ZN3G2=0._r8
  allocate(CO2G2(JZ,JY,JX));    CO2G2=0._r8
  allocate(CO2S2(0:JZ,JY,JX));  CO2S2=0._r8
  allocate(H2GG2(JZ,JY,JX));    H2GG2=0._r8
  allocate(H2GS2(0:JZ,JY,JX));  H2GS2=0._r8
  allocate(H2GSH2(JZ,JY,JX));   H2GSH2=0._r8
  allocate(RQSH2P(2,JV,JH));    RQSH2P=0._r8
  allocate(RQSH1P(2,JV,JH));    RQSH1P=0._r8
  allocate(RQSCOS(2,JV,JH));    RQSCOS=0._r8
  allocate(RQSCHS(2,JV,JH));    RQSCHS=0._r8
  allocate(RQSOXS(2,JV,JH));    RQSOXS=0._r8
  allocate(RQRH2P(2,2,JV,JH));  RQRH2P=0._r8
  allocate(RQRH1P(2,2,JV,JH));  RQRH1P=0._r8
  allocate(RHGFLG(3,JD,JV,JH)); RHGFLG=0._r8
  allocate(THGFLS(JZ,JY,JX));   THGFLS=0._r8
  allocate(ROXFLG(3,JD,JV,JH)); ROXFLG=0._r8
  allocate(RN3FLG(3,JD,JV,JH)); RN3FLG=0._r8
  allocate(RCOFLG(3,JD,JV,JH)); RCOFLG=0._r8
  allocate(RNXFHS(3,JD,JV,JH)); RNXFHS=0._r8
  allocate(RNXFHB(3,JD,JV,JH)); RNXFHB=0._r8
  allocate(RH1PFS(3,0:JD,JV,JH));RH1PFS=0._r8
  allocate(RH1BFB(3,0:JD,JV,JH));RH1BFB=0._r8
  allocate(RH1PHS(3,JD,JV,JH)); RH1PHS=0._r8
  allocate(RH1BHB(3,JD,JV,JH)); RH1BHB=0._r8
  allocate(RCHFLG(3,JD,JV,JH)); RCHFLG=0._r8
  allocate(RNGFLG(3,JD,JV,JH)); RNGFLG=0._r8
  allocate(RN2FLG(3,JD,JV,JH)); RN2FLG=0._r8
  allocate(RN3FHB(3,JD,JV,JH)); RN3FHB=0._r8
  allocate(RNOFHB(3,JD,JV,JH)); RNOFHB=0._r8
  allocate(RH2BHB(3,JD,JV,JH)); RH2BHB=0._r8
  allocate(RNOFHW(3,JD,JV,JH)); RNOFHW=0._r8
  allocate(RH2PHS(3,JD,JV,JH)); RH2PHS=0._r8
  allocate(RN4FHB(3,JD,JV,JH)); RN4FHB=0._r8
  allocate(RN2FHS(3,JD,JV,JH)); RN2FHS=0._r8
  allocate(RN4FHW(3,JD,JV,JH)); RN4FHW=0._r8
  allocate(RN3FHW(3,JD,JV,JH)); RN3FHW=0._r8
  allocate(RHGDFG(0:JZ,JY,JX)); RHGDFG=0._r8
  allocate(RHGFHS(3,JD,JV,JH)); RHGFHS=0._r8
  allocate(RCHFHS(3,JD,JV,JH)); RCHFHS=0._r8
  allocate(ROXFHS(3,JD,JV,JH)); ROXFHS=0._r8
  allocate(RNGFHS(3,JD,JV,JH)); RNGFHS=0._r8
  allocate(RQRNH3(2,2,JV,JH));  RQRNH3=0._r8
  allocate(RQRNO3(2,2,JV,JH));  RQRNO3=0._r8
  allocate(RQRNO2(2,2,JV,JH));  RQRNO2=0._r8
  allocate(RQSNGS(2,JV,JH));    RQSNGS=0._r8
  allocate(RQSN2S(2,JV,JH));    RQSN2S=0._r8
  allocate(RQSNH4(2,JV,JH));    RQSNH4=0._r8
  allocate(RQSNH3(2,JV,JH));    RQSNH3=0._r8
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
  call destroy(RCOBLS)
  call destroy(RCHBLS)
  call destroy(ROXBLS)
  call destroy(RNGBLS)
  call destroy(RN2BLS)
  call destroy(RN4BLW)
  call destroy(RN3BLW)
  call destroy(RNOBLW)
  call destroy(RH1PBS)
  call destroy(RH2PBS)
  call destroy(ROCFL0)
  call destroy(RONFL0)
  call destroy(ROPFL0)
  call destroy(ROAFL0)
  call destroy(ROCFL1)
  call destroy(RONFL1)
  call destroy(ROPFL1)
  call destroy(ROAFL1)
  call destroy(RCOFL0)
  call destroy(RCHFL0)
  call destroy(ROXFL0)
  call destroy(RNGFL0)
  call destroy(RN2FL0)
  call destroy(RHGFL0)
  call destroy(RN4FL0)
  call destroy(RN3FL0)
  call destroy(RNOFL0)
  call destroy(RNXFL0)
  call destroy(RH2PF0)
  call destroy(RCOFL1)
  call destroy(RCHFL1)
  call destroy(ROXFL1)
  call destroy(RNGFL1)
  call destroy(RN2FL1)
  call destroy(RHGFL1)
  call destroy(RN4FL1)
  call destroy(RN3FL1)
  call destroy(RNOFL1)
  call destroy(RNXFL1)
  call destroy(RH2PF1)
  call destroy(RN4FL2)
  call destroy(RN3FL2)
  call destroy(RNOFL2)
  call destroy(RNXFL2)
  call destroy(RH2BF2)
  call destroy(RH1PF0)
  call destroy(RH1PF1)
  call destroy(RH1BF2)
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
  call destroy(CLSGL2)
  call destroy(CQSGL2)
  call destroy(POSGL2)
  call destroy(OLSGL2)
  call destroy(ZNSGL2)
  call destroy(ZLSGL2)
  call destroy(ZVSGL2)
  call destroy(HLSGL2)
  call destroy(ZOSGL2)
  call destroy(RCODFG)
  call destroy(RCHDFG)
  call destroy(ROXDFG)
  call destroy(RNGDFG)
  call destroy(RN2DFG)
  call destroy(RN3DFG)
  call destroy(RNBDFG)
  call destroy(TQROC)
  call destroy(TQRON)
  call destroy(TQROP)
  call destroy(TQROA)
  call destroy(TQRCHS)
  call destroy(TQROXS)
  call destroy(TQRNGS)
  call destroy(TQRN2S)
  call destroy(TQRNH4)
  call destroy(TQRNH3)
  call destroy(TQRNO3)
  call destroy(TQRH2P)
  call destroy(TQRNO2)
  call destroy(TQRHGS)
  call destroy(TQSCOS)
  call destroy(TQRCOS)
  call destroy(TQSCHS)
  call destroy(TQSOXS)
  call destroy(TQSNGS)
  call destroy(TQSN2S)
  call destroy(TQSNH4)
  call destroy(TQSNH3)
  call destroy(TQSNO3)
  call destroy(TQSH1P)
  call destroy(TQSH2P)
  call destroy(TCOFLS)
  call destroy(TCHFLS)
  call destroy(TOXFLS)
  call destroy(TNGFLS)
  call destroy(TN2FLS)
  call destroy(TN4FLW)
  call destroy(TN3FLW)
  call destroy(TNOFLW)
  call destroy(TH2PFS)
  call destroy(TN4FLB)
  call destroy(TN3FLB)
  call destroy(TNOFLB)
  call destroy(TH2BFB)
  call destroy(TNXFLS)
  call destroy(TCOFLG)
  call destroy(TCHFLG)
  call destroy(TOXFLG)
  call destroy(TNGFLG)
  call destroy(TN2FLG)
  call destroy(TQRH1P)
  call destroy(TH1PFS)
  call destroy(TH1BFB)
  call destroy(DCO2G)
  call destroy(DCH4G)
  call destroy(DOXYG)
  call destroy(DZ2GG)
  call destroy(DZ2OG)
  call destroy(DNH3G)
  call destroy(VOLWCO)
  call destroy(VOLWCH)
  call destroy(VOLWOX)
  call destroy(VOLWNG)
  call destroy(VOLWN2)
  call destroy(VOLWN3)
  call destroy(VOLWNB)
  call destroy(VOLWHG)
  call destroy(HGSGL2)
  call destroy(DH2GG)
  call destroy(RHGFXS)
  call destroy(THGFHS)
  call destroy(THGFLG)
  call destroy(FLVM)
  call destroy(THETH2)
  call destroy(THETHL)
  call destroy(VOLPMA)
  call destroy(VOLPMB)
  call destroy(VOLWMA)
  call destroy(VOLWMB)
  call destroy(VOLWXA)
  call destroy(VOLWXB)
  call destroy(PARGCO)
  call destroy(PARGCH)
  call destroy(PARGOX)
  call destroy(PARGNG)
  call destroy(PARGN2)
  call destroy(PARGN3)
  call destroy(PARGH2)
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
  call destroy(RCOBBL)
  call destroy(RCHBBL)
  call destroy(ROXBBL)
  call destroy(RNGBBL)
  call destroy(RN2BBL)
  call destroy(RN3BBL)
  call destroy(RNBBBL)
  call destroy(RHGBBL)
  call destroy(RQRCOS0)
  call destroy(RQRCHS0)
  call destroy(RQROXS0)
  call destroy(RQRNGS0)
  call destroy(RQRN2S0)
  call destroy(RQRHGS0)
  call destroy(RQRNH40)
  call destroy(RQRNH30)
  call destroy(RQRNO30)
  call destroy(RQRNO20)
  call destroy(RQRH2P0)
  call destroy(RQRH1P0)
  call destroy(OQNH2)
  call destroy(OQPH2)
  call destroy(OQAH2)
  call destroy(TOCFHS)
  call destroy(TONFHS)
  call destroy(TOPFHS)
  call destroy(TOAFHS)
  call destroy(TCOFHS)
  call destroy(TCHFHS)
  call destroy(TOXFHS)
  call destroy(TNGFHS)
  call destroy(TN2FHS)
  call destroy(TN4FHW)
  call destroy(TN3FHW)
  call destroy(TNOFHW)
  call destroy(TH2PHS)
  call destroy(TN4FHB)
  call destroy(TN3FHB)
  call destroy(TNOFHB)
  call destroy(TH2BHB)
  call destroy(TNXFHS)
  call destroy(ZNO2B2)
  call destroy(ZN2BH2)
  call destroy(TNXFLB)
  call destroy(TNXFHB)
  call destroy(H1P4H2)
  call destroy(H1PBH2)
  call destroy(TH1PHS)
  call destroy(TH1BHB)
  call destroy(H1PO42)
  call destroy(H1POB2)
  call destroy(OCSGL2)
  call destroy(ONSGL2)
  call destroy(OPSGL2)
  call destroy(OASGL2)
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
  call destroy(TCOBLS)
  call destroy(TCHBLS)
  call destroy(TOXBLS)
  call destroy(TNGBLS)
  call destroy(TN2BLS)
  call destroy(TN4BLW)
  call destroy(TN3BLW)
  call destroy(TNOBLW)
  call destroy(TH1PBS)
  call destroy(TH2PBS)
  call destroy(RQROP)
  call destroy(RQROA)
  call destroy(RQRCOS)
  call destroy(RQRCHS)
  call destroy(RQROXS)
  call destroy(RQRNGS)
  call destroy(RQRN2S)
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
  call destroy(RCODFR)
  call destroy(RCHDFR)
  call destroy(ROXDFR)
  call destroy(RNGDFR)
  call destroy(RN2DFR)
  call destroy(RN3DFR)
  call destroy(RHGDFR)
  call destroy(R1BSK2)
  call destroy(RQROC)
  call destroy(RQRON)
  call destroy(ZNH4S2)
  call destroy(ZNH4B2)
  call destroy(ZNH3S2)
  call destroy(ZNH3B2)
  call destroy(ZNO3S2)
  call destroy(ZNO3B2)
  call destroy(H2PO42)
  call destroy(H2POB2)
  call destroy(ZNO2S2)
  call destroy(CGSGL2)
  call destroy(CHSGL2)
  call destroy(OGSGL2)
  call destroy(ZGSGL2)
  call destroy(Z2SGL2)
  call destroy(ZHSGL2)
  call destroy(CO2SH2)
  call destroy(CH4SH2)
  call destroy(OXYSH2)
  call destroy(Z2GSH2)
  call destroy(Z2OSH2)
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
  call destroy(CH4G2)
  call destroy(CH4S2)
  call destroy(OXYG2)
  call destroy(OXYS2)
  call destroy(Z2GG2)
  call destroy(Z2GS2)
  call destroy(Z2OG2)
  call destroy(Z2OS2)
  call destroy(ZN3G2)
  call destroy(CO2G2)
  call destroy(CO2S2)
  call destroy(H2GG2)
  call destroy(H2GS2)
  call destroy(H2GSH2)
  call destroy(RQSH2P)
  call destroy(RQSH1P)
  call destroy(RQSCOS)
  call destroy(RQSCHS)
  call destroy(RQSOXS)
  call destroy(RQRH2P)
  call destroy(RQRH1P)
  call destroy(RHGFLG)
  call destroy(THGFLS)
  call destroy(ROXFLG)
  call destroy(RN3FLG)
  call destroy(RCOFLG)
  call destroy(RNXFHS)
  call destroy(RNXFHB)
  call destroy(RH1PFS)
  call destroy(RH1BFB)
  call destroy(RH1PHS)
  call destroy(RH1BHB)
  call destroy(RCHFLG)
  call destroy(RNGFLG)
  call destroy(RN2FLG)
  call destroy(RN3FHB)
  call destroy(RNOFHB)
  call destroy(RH2BHB)
  call destroy(RNOFHW)
  call destroy(RH2PHS)
  call destroy(RN4FHB)
  call destroy(RN2FHS)
  call destroy(RN4FHW)
  call destroy(RN3FHW)
  call destroy(RHGDFG)
  call destroy(RHGFHS)
  call destroy(RCHFHS)
  call destroy(ROXFHS)
  call destroy(RNGFHS)
  call destroy(RQRNH3)
  call destroy(RQRNO3)
  call destroy(RQRNO2)
  call destroy(RQSNGS)
  call destroy(RQSN2S)
  call destroy(RQSNH4)
  call destroy(RQSNH3)
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
  end subroutine DestructTransfrData


end module TransfrDataMod
