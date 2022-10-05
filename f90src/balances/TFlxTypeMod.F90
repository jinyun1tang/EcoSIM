module TFlxTypeMod

  use GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
  use EcoSIMConfig, only : jcplx1 => jcplx1c,jsken=>jskenc,NFGs=>NFGsc
  use EcoSIMConfig, only : nlbiomcp=>nlbiomcpc,ndbiomcp=>ndbiomcpc
implicit none

  character(len=*), private, parameter :: mod_filename = __FILE__

  public
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
  real(r8),allocatable ::  TALBLS(:,:,:)                      !
  real(r8),allocatable ::  TFEBLS(:,:,:)                      !
  real(r8),allocatable ::  THYBLS(:,:,:)                      !
  real(r8),allocatable ::  TCABLS(:,:,:)                      !
  real(r8),allocatable ::  TMGBLS(:,:,:)                      !
  real(r8),allocatable ::  TNABLS(:,:,:)                      !
  real(r8),allocatable ::  TKABLS(:,:,:)                      !
  real(r8),allocatable ::  TOHBLS(:,:,:)                      !
  real(r8),allocatable ::  TSOBLS(:,:,:)                      !
  real(r8),allocatable ::  TCLBLS(:,:,:)                      !
  real(r8),allocatable ::  TC3BLS(:,:,:)                      !
  real(r8),allocatable ::  THCBLS(:,:,:)                      !
  real(r8),allocatable ::  TAL1BS(:,:,:)                      !
  real(r8),allocatable ::  TAL2BS(:,:,:)                      !
  real(r8),allocatable ::  TAL3BS(:,:,:)                      !
  real(r8),allocatable ::  TAL4BS(:,:,:)                      !
  real(r8),allocatable ::  TALSBS(:,:,:)                      !
  real(r8),allocatable ::  TFE1BS(:,:,:)                      !
  real(r8),allocatable ::  TFE2BS(:,:,:)                      !
  real(r8),allocatable ::  TFE3BS(:,:,:)                      !
  real(r8),allocatable ::  TFE4BS(:,:,:)                      !
  real(r8),allocatable ::  TFESBS(:,:,:)                      !
  real(r8),allocatable ::  TCAOBS(:,:,:)                      !
  real(r8),allocatable ::  TCACBS(:,:,:)                      !
  real(r8),allocatable ::  TCAHBS(:,:,:)                      !
  real(r8),allocatable ::  TCASBS(:,:,:)                      !
  real(r8),allocatable ::  TMGOBS(:,:,:)                      !
  real(r8),allocatable ::  TMGCBS(:,:,:)                      !
  real(r8),allocatable ::  TMGHBS(:,:,:)                      !
  real(r8),allocatable ::  TMGSBS(:,:,:)                      !
  real(r8),allocatable ::  TNACBS(:,:,:)                      !
  real(r8),allocatable ::  TNASBS(:,:,:)                      !
  real(r8),allocatable ::  TKASBS(:,:,:)                      !
  real(r8),allocatable ::  TH0PBS(:,:,:)                      !
  real(r8),allocatable ::  TH3PBS(:,:,:)                      !
  real(r8),allocatable ::  TF1PBS(:,:,:)                      !
  real(r8),allocatable ::  TF2PBS(:,:,:)                      !
  real(r8),allocatable ::  TC0PBS(:,:,:)                      !
  real(r8),allocatable ::  TC1PBS(:,:,:)                      !
  real(r8),allocatable ::  TC2PBS(:,:,:)                      !
  real(r8),allocatable ::  TM1PBS(:,:,:)                      !
  real(r8),allocatable ::  THGFLG(:,:,:)                      !
  real(r8),allocatable ::  THGFLS(:,:,:)                      !
  real(r8),allocatable ::  THGFHS(:,:,:)                      !
  real(r8),allocatable ::  TQR(:,:)                           !
  real(r8),allocatable ::  THQR(:,:)                          !
  real(r8),allocatable ::  TQS(:,:)                           !
  real(r8),allocatable ::  TQW(:,:)                           !
  real(r8),allocatable ::  TQI(:,:)                           !
  real(r8),allocatable ::  THQS(:,:)                          !
  real(r8),allocatable ::  TFLWS(:,:,:)                       !
  real(r8),allocatable ::  TFLWW(:,:,:)                       !
  real(r8),allocatable ::  TFLWI(:,:,:)                       !
  real(r8),allocatable ::  THFLWW(:,:,:)                      !
  real(r8),allocatable ::  THGQRS(:,:)                        !
  real(r8),allocatable ::  TCOQRS(:,:)                        !
  real(r8),allocatable ::  TCHQRS(:,:)                        !
  real(r8),allocatable ::  TOXQRS(:,:)                        !
  real(r8),allocatable ::  TNGQRS(:,:)                        !
  real(r8),allocatable ::  TN2QRS(:,:)                        !
  real(r8),allocatable ::  TN4QRS(:,:)                        !
  real(r8),allocatable ::  TN3QRS(:,:)                        !
  real(r8),allocatable ::  TNOQRS(:,:)                        !
  real(r8),allocatable ::  TPOQRS(:,:)                        !
  real(r8),allocatable ::  TNXQRS(:,:)                        !
  real(r8),allocatable ::  TQRAL(:,:)                         !
  real(r8),allocatable ::  TQRFE(:,:)                         !
  real(r8),allocatable ::  TQRHY(:,:)                         !
  real(r8),allocatable ::  TQRCA(:,:)                         !
  real(r8),allocatable ::  TQRMG(:,:)                         !
  real(r8),allocatable ::  TQRNA(:,:)                         !
  real(r8),allocatable ::  TQRKA(:,:)                         !
  real(r8),allocatable ::  TQROH(:,:)                         !
  real(r8),allocatable ::  TQRSO(:,:)                         !
  real(r8),allocatable ::  TQRCL(:,:)                         !
  real(r8),allocatable ::  TQRC3(:,:)                         !
  real(r8),allocatable ::  TQRHC(:,:)                         !
  real(r8),allocatable ::  TQRAL1(:,:)                        !
  real(r8),allocatable ::  TQRAL2(:,:)                        !
  real(r8),allocatable ::  TQRAL3(:,:)                        !
  real(r8),allocatable ::  TQRAL4(:,:)                        !
  real(r8),allocatable ::  TQRALS(:,:)                        !
  real(r8),allocatable ::  TQRFE1(:,:)                        !
  real(r8),allocatable ::  TQRFE2(:,:)                        !
  real(r8),allocatable ::  TQRFE3(:,:)                        !
  real(r8),allocatable ::  TQRFE4(:,:)                        !
  real(r8),allocatable ::  TQRFES(:,:)                        !
  real(r8),allocatable ::  TQRCAO(:,:)                        !
  real(r8),allocatable ::  TQRCAC(:,:)                        !
  real(r8),allocatable ::  TQRCAH(:,:)                        !
  real(r8),allocatable ::  TQRCAS(:,:)                        !
  real(r8),allocatable ::  TQRMGO(:,:)                        !
  real(r8),allocatable ::  TQRMGC(:,:)                        !
  real(r8),allocatable ::  TQRMGH(:,:)                        !
  real(r8),allocatable ::  TQRMGS(:,:)                        !
  real(r8),allocatable ::  TQRNAC(:,:)                        !
  real(r8),allocatable ::  TQRNAS(:,:)                        !
  real(r8),allocatable ::  TQRKAS(:,:)                        !
  real(r8),allocatable ::  TQRH0P(:,:)                        !
  real(r8),allocatable ::  TQRH3P(:,:)                        !
  real(r8),allocatable ::  TQRF1P(:,:)                        !
  real(r8),allocatable ::  TP1QRS(:,:)                        !
  real(r8),allocatable ::  TQRF2P(:,:)                        !
  real(r8),allocatable ::  TQRC0P(:,:)                        !
  real(r8),allocatable ::  TQRC1P(:,:)                        !
  real(r8),allocatable ::  TQRC2P(:,:)                        !
  real(r8),allocatable ::  TQRM1P(:,:)                        !
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
  real(r8),allocatable ::  TQSAL(:,:)                         !
  real(r8),allocatable ::  TQSFE(:,:)                         !
  real(r8),allocatable ::  TQSHY(:,:)                         !
  real(r8),allocatable ::  TQSCA(:,:)                         !
  real(r8),allocatable ::  TQSMG(:,:)                         !
  real(r8),allocatable ::  TQSNA(:,:)                         !
  real(r8),allocatable ::  TQSKA(:,:)                         !
  real(r8),allocatable ::  TQSOH(:,:)                         !
  real(r8),allocatable ::  TQSSO(:,:)                         !
  real(r8),allocatable ::  TQSCL(:,:)                         !
  real(r8),allocatable ::  TQSC3(:,:)                         !
  real(r8),allocatable ::  TQSHC(:,:)                         !
  real(r8),allocatable ::  TQSAL1(:,:)                        !
  real(r8),allocatable ::  TQSAL2(:,:)                        !
  real(r8),allocatable ::  TQSAL3(:,:)                        !
  real(r8),allocatable ::  TQSAL4(:,:)                        !
  real(r8),allocatable ::  TQSALS(:,:)                        !
  real(r8),allocatable ::  TQSFE1(:,:)                        !
  real(r8),allocatable ::  TQSFE2(:,:)                        !
  real(r8),allocatable ::  TQSFE3(:,:)                        !
  real(r8),allocatable ::  TQSFE4(:,:)                        !
  real(r8),allocatable ::  TQSFES(:,:)                        !
  real(r8),allocatable ::  TQSCAO(:,:)                        !
  real(r8),allocatable ::  TQSCAC(:,:)                        !
  real(r8),allocatable ::  TQSCAH(:,:)                        !
  real(r8),allocatable ::  TQSCAS(:,:)                        !
  real(r8),allocatable ::  TQSMGO(:,:)                        !
  real(r8),allocatable ::  TQSMGC(:,:)                        !
  real(r8),allocatable ::  TQSMGH(:,:)                        !
  real(r8),allocatable ::  TQSMGS(:,:)                        !
  real(r8),allocatable ::  TQSNAC(:,:)                        !
  real(r8),allocatable ::  TQSNAS(:,:)                        !
  real(r8),allocatable ::  TQSKAS(:,:)                        !
  real(r8),allocatable ::  TQSH0P(:,:)                        !
  real(r8),allocatable ::  TQSH3P(:,:)                        !
  real(r8),allocatable ::  TQSF1P(:,:)                        !
  real(r8),allocatable ::  TQSF2P(:,:)                        !
  real(r8),allocatable ::  TQSC0P(:,:)                        !
  real(r8),allocatable ::  TQSC1P(:,:)                        !
  real(r8),allocatable ::  TQSC2P(:,:)                        !
  real(r8),allocatable ::  TQSM1P(:,:)                        !
  real(r8),allocatable ::  TALFLS(:,:,:)                      !
  real(r8),allocatable ::  TFEFLS(:,:,:)                      !
  real(r8),allocatable ::  TCAFLS(:,:,:)                      !
  real(r8),allocatable ::  THYFLS(:,:,:)                      !
  real(r8),allocatable ::  TMGFLS(:,:,:)                      !
  real(r8),allocatable ::  TNAFLS(:,:,:)                      !
  real(r8),allocatable ::  TKAFLS(:,:,:)                      !
  real(r8),allocatable ::  TOHFLS(:,:,:)                      !
  real(r8),allocatable ::  TSOFLS(:,:,:)                      !
  real(r8),allocatable ::  TCLFLS(:,:,:)                      !
  real(r8),allocatable ::  TC3FLS(:,:,:)                      !
  real(r8),allocatable ::  THCFLS(:,:,:)                      !
  real(r8),allocatable ::  TAL1FS(:,:,:)                      !
  real(r8),allocatable ::  TAL2FS(:,:,:)                      !
  real(r8),allocatable ::  TAL3FS(:,:,:)                      !
  real(r8),allocatable ::  TAL4FS(:,:,:)                      !
  real(r8),allocatable ::  TALSFS(:,:,:)                      !
  real(r8),allocatable ::  TFE1FS(:,:,:)                      !
  real(r8),allocatable ::  TFE2FS(:,:,:)                      !
  real(r8),allocatable ::  TFE3FS(:,:,:)                      !
  real(r8),allocatable ::  TFE4FS(:,:,:)                      !
  real(r8),allocatable ::  TFESFS(:,:,:)                      !
  real(r8),allocatable ::  TCAOFS(:,:,:)                      !
  real(r8),allocatable ::  TCACFS(:,:,:)                      !
  real(r8),allocatable ::  TCAHFS(:,:,:)                      !
  real(r8),allocatable ::  TCASFS(:,:,:)                      !
  real(r8),allocatable ::  TMGOFS(:,:,:)                      !
  real(r8),allocatable ::  TMGCFS(:,:,:)                      !
  real(r8),allocatable ::  TMGHFS(:,:,:)                      !
  real(r8),allocatable ::  TMGSFS(:,:,:)                      !
  real(r8),allocatable ::  TNACFS(:,:,:)                      !
  real(r8),allocatable ::  TNASFS(:,:,:)                      !
  real(r8),allocatable ::  TKASFS(:,:,:)                      !
  real(r8),allocatable ::  TH0PFS(:,:,:)                      !
  real(r8),allocatable ::  TH3PFS(:,:,:)                      !
  real(r8),allocatable ::  TF1PFS(:,:,:)                      !
  real(r8),allocatable ::  TF2PFS(:,:,:)                      !
  real(r8),allocatable ::  TC0PFS(:,:,:)                      !
  real(r8),allocatable ::  TC1PFS(:,:,:)                      !
  real(r8),allocatable ::  TC2PFS(:,:,:)                      !
  real(r8),allocatable ::  TM1PFS(:,:,:)                      !
  real(r8),allocatable ::  TH0BFB(:,:,:)                      !
  real(r8),allocatable ::  TH3BFB(:,:,:)                      !
  real(r8),allocatable ::  TF1BFB(:,:,:)                      !
  real(r8),allocatable ::  TF2BFB(:,:,:)                      !
  real(r8),allocatable ::  TC0BFB(:,:,:)                      !
  real(r8),allocatable ::  TC1BFB(:,:,:)                      !
  real(r8),allocatable ::  TC2BFB(:,:,:)                      !
  real(r8),allocatable ::  TM1BFB(:,:,:)                      !
  real(r8),allocatable ::  TALFHS(:,:,:)                      !
  real(r8),allocatable ::  TFEFHS(:,:,:)                      !
  real(r8),allocatable ::  THYFHS(:,:,:)                      !
  real(r8),allocatable ::  TCAFHS(:,:,:)                      !
  real(r8),allocatable ::  TMGFHS(:,:,:)                      !
  real(r8),allocatable ::  TNAFHS(:,:,:)                      !
  real(r8),allocatable ::  TKAFHS(:,:,:)                      !
  real(r8),allocatable ::  TOHFHS(:,:,:)                      !
  real(r8),allocatable ::  TSOFHS(:,:,:)                      !
  real(r8),allocatable ::  TCLFHS(:,:,:)                      !
  real(r8),allocatable ::  TC3FHS(:,:,:)                      !
  real(r8),allocatable ::  THCFHS(:,:,:)                      !
  real(r8),allocatable ::  TAL1HS(:,:,:)                      !
  real(r8),allocatable ::  TAL2HS(:,:,:)                      !
  real(r8),allocatable ::  TAL3HS(:,:,:)                      !
  real(r8),allocatable ::  TAL4HS(:,:,:)                      !
  real(r8),allocatable ::  TALSHS(:,:,:)                      !
  real(r8),allocatable ::  TFE1HS(:,:,:)                      !
  real(r8),allocatable ::  TFE2HS(:,:,:)                      !
  real(r8),allocatable ::  TFE3HS(:,:,:)                      !
  real(r8),allocatable ::  TFE4HS(:,:,:)                      !
  real(r8),allocatable ::  TFESHS(:,:,:)                      !
  real(r8),allocatable ::  TCAOHS(:,:,:)                      !
  real(r8),allocatable ::  TCACHS(:,:,:)                      !
  real(r8),allocatable ::  TCAHHS(:,:,:)                      !
  real(r8),allocatable ::  TCASHS(:,:,:)                      !
  real(r8),allocatable ::  TMGOHS(:,:,:)                      !
  real(r8),allocatable ::  TMGCHS(:,:,:)                      !
  real(r8),allocatable ::  TMGHHS(:,:,:)                      !
  real(r8),allocatable ::  TMGSHS(:,:,:)                      !
  real(r8),allocatable ::  TNACHS(:,:,:)                      !
  real(r8),allocatable ::  TNASHS(:,:,:)                      !
  real(r8),allocatable ::  TKASHS(:,:,:)                      !
  real(r8),allocatable ::  TH0PHS(:,:,:)                      !
  real(r8),allocatable ::  TH3PHS(:,:,:)                      !
  real(r8),allocatable ::  TF1PHS(:,:,:)                      !
  real(r8),allocatable ::  TF2PHS(:,:,:)                      !
  real(r8),allocatable ::  TC0PHS(:,:,:)                      !
  real(r8),allocatable ::  TC1PHS(:,:,:)                      !
  real(r8),allocatable ::  TC2PHS(:,:,:)                      !
  real(r8),allocatable ::  TM1PHS(:,:,:)                      !
  real(r8),allocatable ::  TH0BHB(:,:,:)                      !
  real(r8),allocatable ::  TH3BHB(:,:,:)                      !
  real(r8),allocatable ::  TF1BHB(:,:,:)                      !
  real(r8),allocatable ::  TF2BHB(:,:,:)                      !
  real(r8),allocatable ::  TC0BHB(:,:,:)                      !
  real(r8),allocatable ::  TC1BHB(:,:,:)                      !
  real(r8),allocatable ::  TC2BHB(:,:,:)                      !
  real(r8),allocatable ::  TM1BHB(:,:,:)                      !
  real(r8),allocatable ::  TSANER(:,:)                        !
  real(r8),allocatable ::  TSILER(:,:)                        !
  real(r8),allocatable ::  TCLAER(:,:)                        !
  real(r8),allocatable ::  TCECER(:,:)                        !
  real(r8),allocatable ::  TAECER(:,:)                        !
  real(r8),allocatable ::  TNH4ER(:,:)                        !
  real(r8),allocatable ::  TNH3ER(:,:)                        !
  real(r8),allocatable ::  TNHUER(:,:)                        !
  real(r8),allocatable ::  TNO3ER(:,:)                        !
  real(r8),allocatable ::  TNH4EB(:,:)                        !
  real(r8),allocatable ::  TNH3EB(:,:)                        !
  real(r8),allocatable ::  TNHUEB(:,:)                        !
  real(r8),allocatable ::  TNO3EB(:,:)                        !
  real(r8),allocatable ::  TN4ER(:,:)                         !
  real(r8),allocatable ::  TNBER(:,:)                         !
  real(r8),allocatable ::  THYER(:,:)                         !
  real(r8),allocatable ::  TALER(:,:)                         !
  real(r8),allocatable ::  TCAER(:,:)                         !
  real(r8),allocatable ::  TMGER(:,:)                         !
  real(r8),allocatable ::  TNAER(:,:)                         !
  real(r8),allocatable ::  TKAER(:,:)                         !
  real(r8),allocatable ::  THCER(:,:)                         !
  real(r8),allocatable ::  TAL2ER(:,:)                        !
  real(r8),allocatable ::  TOH0ER(:,:)                        !
  real(r8),allocatable ::  TOH1ER(:,:)                        !
  real(r8),allocatable ::  TOH2ER(:,:)                        !
  real(r8),allocatable ::  TH1PER(:,:)                        !
  real(r8),allocatable ::  TH2PER(:,:)                        !
  real(r8),allocatable ::  TOH0EB(:,:)                        !
  real(r8),allocatable ::  TOH1EB(:,:)                        !
  real(r8),allocatable ::  TOH2EB(:,:)                        !
  real(r8),allocatable ::  TH1PEB(:,:)                        !
  real(r8),allocatable ::  TH2PEB(:,:)                        !
  real(r8),allocatable ::  TALOER(:,:)                        !
  real(r8),allocatable ::  TFEOER(:,:)                        !
  real(r8),allocatable ::  TCACER(:,:)                        !
  real(r8),allocatable ::  TCASER(:,:)                        !
  real(r8),allocatable ::  TALPER(:,:)                        !
  real(r8),allocatable ::  TFEPER(:,:)                        !
  real(r8),allocatable ::  TCPDER(:,:)                        !
  real(r8),allocatable ::  TCPHER(:,:)                        !
  real(r8),allocatable ::  TCPMER(:,:)                        !
  real(r8),allocatable ::  TALPEB(:,:)                        !
  real(r8),allocatable ::  TFEPEB(:,:)                        !
  real(r8),allocatable ::  TCPDEB(:,:)                        !
  real(r8),allocatable ::  TCPHEB(:,:)                        !
  real(r8),allocatable ::  TCPMEB(:,:)                        !
  real(r8),allocatable ::  TFEER(:,:)                         !
  real(r8),allocatable ::  TFE2ER(:,:)                        !
  real(r8),allocatable ::  TSEDER(:,:)                        !
  real(r8),allocatable ::  TFLW(:,:,:)                        !
  real(r8),allocatable ::  TFLWX(:,:,:)                       !
  real(r8),allocatable ::  THFLW(:,:,:)                       !
  real(r8),allocatable ::  TFLWH(:,:,:)                       !
  real(r8),allocatable ::  TCOFLS(:,:,:)                      !
  real(r8),allocatable ::  TCHFLS(:,:,:)                      !
  real(r8),allocatable ::  TOXFLS(:,:,:)                      !
  real(r8),allocatable ::  TNXFLB(:,:,:)                      !
  real(r8),allocatable ::  TNGFLS(:,:,:)                      !
  real(r8),allocatable ::  TN2FLS(:,:,:)                      !
  real(r8),allocatable ::  TN4FLS(:,:,:)                      !
  real(r8),allocatable ::  TN4FLB(:,:,:)                      !
  real(r8),allocatable ::  TN3FLS(:,:,:)                      !
  real(r8),allocatable ::  TN3FLB(:,:,:)                      !
  real(r8),allocatable ::  TNOFLS(:,:,:)                      !
  real(r8),allocatable ::  TNOFLB(:,:,:)                      !
  real(r8),allocatable ::  TPOFLS(:,:,:)                      !
  real(r8),allocatable ::  TH2BFB(:,:,:)                      !
  real(r8),allocatable ::  TNXFLS(:,:,:)                      !
  real(r8),allocatable ::  TCOFHS(:,:,:)                      !
  real(r8),allocatable ::  TCHFHS(:,:,:)                      !
  real(r8),allocatable ::  TNXFHB(:,:,:)                      !
  real(r8),allocatable ::  TOXFHS(:,:,:)                      !
  real(r8),allocatable ::  TNGFHS(:,:,:)                      !
  real(r8),allocatable ::  TN2FHS(:,:,:)                      !
  real(r8),allocatable ::  TN4FHS(:,:,:)                      !
  real(r8),allocatable ::  TN4FHB(:,:,:)                      !
  real(r8),allocatable ::  TN3FHS(:,:,:)                      !
  real(r8),allocatable ::  TN3FHB(:,:,:)                      !
  real(r8),allocatable ::  TNOFHS(:,:,:)                      !
  real(r8),allocatable ::  TNOFHB(:,:,:)                      !
  real(r8),allocatable ::  TPOFHS(:,:,:)                      !
  real(r8),allocatable ::  TH2BHB(:,:,:)                      !
  real(r8),allocatable ::  TNXFHS(:,:,:)                      !
  real(r8),allocatable ::  TCOFLG(:,:,:)                      !
  real(r8),allocatable ::  TCHFLG(:,:,:)                      !
  real(r8),allocatable ::  TOXFLG(:,:,:)                      !
  real(r8),allocatable ::  TNGFLG(:,:,:)                      !
  real(r8),allocatable ::  TN2FLG(:,:,:)                      !
  real(r8),allocatable ::  TNHFLG(:,:,:)                      !
  real(r8),allocatable ::  TTHAW(:,:,:)                       !
  real(r8),allocatable ::  THTHAW(:,:,:)                      !
  real(r8),allocatable ::  TTHAWH(:,:,:)                      !
  real(r8),allocatable ::  TP1FLS(:,:,:)                      !
  real(r8),allocatable ::  TP1FHS(:,:,:)                      !
  real(r8),allocatable ::  TH1BFB(:,:,:)                      !
  real(r8),allocatable ::  TH1BHB(:,:,:)                      !
  real(r8),allocatable ::  VOLW1(:,:,:)                       !
  real(r8),allocatable ::  VOLI1(:,:,:)                       !
  real(r8),allocatable ::  VOLWH1(:,:,:)                      !
  real(r8),allocatable ::  VOLIH1(:,:,:)                      !

  real(r8),allocatable :: TOMCER(:,:,:,:,:,:)
  real(r8),allocatable :: TOMNER(:,:,:,:,:,:)
  real(r8),allocatable :: TOMPER(:,:,:,:,:,:)


  real(r8),allocatable :: TOMCERff(:,:,:,:,:)
  real(r8),allocatable :: TOMNERff(:,:,:,:,:)
  real(r8),allocatable :: TOMPERff(:,:,:,:,:)

  real(r8),allocatable ::  TOCFLS(:,:,:,:)
  real(r8),allocatable ::  TONFLS(:,:,:,:)
  real(r8),allocatable ::  TOPFLS(:,:,:,:)
  real(r8),allocatable ::  TOAFLS(:,:,:,:)
  real(r8),allocatable ::  TOCFHS(:,:,:,:)
  real(r8),allocatable ::  TONFHS(:,:,:,:)
  real(r8),allocatable ::  TOPFHS(:,:,:,:)
  real(r8),allocatable ::  TOAFHS(:,:,:,:)
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

  real(r8) :: TDLYXF,TDYLXC,TDVOLI,TDORGC
  DATA TDORGC,TDYLXC/0.0,0.0/
  DATA TDVOLI,TDLYXF/0.0,0.0/


  contains

!------------------------------------------------------------------------------------------

  subroutine InitTflxType()

  implicit none
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
  allocate(TALBLS(JS,JY,JX));   TALBLS=0._r8
  allocate(TFEBLS(JS,JY,JX));   TFEBLS=0._r8
  allocate(THYBLS(JS,JY,JX));   THYBLS=0._r8
  allocate(TCABLS(JS,JY,JX));   TCABLS=0._r8
  allocate(TMGBLS(JS,JY,JX));   TMGBLS=0._r8
  allocate(TNABLS(JS,JY,JX));   TNABLS=0._r8
  allocate(TKABLS(JS,JY,JX));   TKABLS=0._r8
  allocate(TOHBLS(JS,JY,JX));   TOHBLS=0._r8
  allocate(TSOBLS(JS,JY,JX));   TSOBLS=0._r8
  allocate(TCLBLS(JS,JY,JX));   TCLBLS=0._r8
  allocate(TC3BLS(JS,JY,JX));   TC3BLS=0._r8
  allocate(THCBLS(JS,JY,JX));   THCBLS=0._r8
  allocate(TAL1BS(JS,JY,JX));   TAL1BS=0._r8
  allocate(TAL2BS(JS,JY,JX));   TAL2BS=0._r8
  allocate(TAL3BS(JS,JY,JX));   TAL3BS=0._r8
  allocate(TAL4BS(JS,JY,JX));   TAL4BS=0._r8
  allocate(TALSBS(JS,JY,JX));   TALSBS=0._r8
  allocate(TFE1BS(JS,JY,JX));   TFE1BS=0._r8
  allocate(TFE2BS(JS,JY,JX));   TFE2BS=0._r8
  allocate(TFE3BS(JS,JY,JX));   TFE3BS=0._r8
  allocate(TFE4BS(JS,JY,JX));   TFE4BS=0._r8
  allocate(TFESBS(JS,JY,JX));   TFESBS=0._r8
  allocate(TCAOBS(JS,JY,JX));   TCAOBS=0._r8
  allocate(TCACBS(JS,JY,JX));   TCACBS=0._r8
  allocate(TCAHBS(JS,JY,JX));   TCAHBS=0._r8
  allocate(TCASBS(JS,JY,JX));   TCASBS=0._r8
  allocate(TMGOBS(JS,JY,JX));   TMGOBS=0._r8
  allocate(TMGCBS(JS,JY,JX));   TMGCBS=0._r8
  allocate(TMGHBS(JS,JY,JX));   TMGHBS=0._r8
  allocate(TMGSBS(JS,JY,JX));   TMGSBS=0._r8
  allocate(TNACBS(JS,JY,JX));   TNACBS=0._r8
  allocate(TNASBS(JS,JY,JX));   TNASBS=0._r8
  allocate(TKASBS(JS,JY,JX));   TKASBS=0._r8
  allocate(TH0PBS(JS,JY,JX));   TH0PBS=0._r8
  allocate(TH3PBS(JS,JY,JX));   TH3PBS=0._r8
  allocate(TF1PBS(JS,JY,JX));   TF1PBS=0._r8
  allocate(TF2PBS(JS,JY,JX));   TF2PBS=0._r8
  allocate(TC0PBS(JS,JY,JX));   TC0PBS=0._r8
  allocate(TC1PBS(JS,JY,JX));   TC1PBS=0._r8
  allocate(TC2PBS(JS,JY,JX));   TC2PBS=0._r8
  allocate(TM1PBS(JS,JY,JX));   TM1PBS=0._r8
  allocate(THGFLG(JZ,JY,JX));   THGFLG=0._r8
  allocate(THGFLS(JZ,JY,JX));   THGFLS=0._r8
  allocate(THGFHS(JZ,JY,JX));   THGFHS=0._r8
  allocate(TQR(JY,JX));         TQR=0._r8
  allocate(THQR(JY,JX));        THQR=0._r8
  allocate(TQS(JY,JX));         TQS=0._r8
  allocate(TQW(JY,JX));         TQW=0._r8
  allocate(TQI(JY,JX));         TQI=0._r8
  allocate(THQS(JY,JX));        THQS=0._r8
  allocate(TFLWS(JS,JY,JX));    TFLWS=0._r8
  allocate(TFLWW(JS,JY,JX));    TFLWW=0._r8
  allocate(TFLWI(JS,JY,JX));    TFLWI=0._r8
  allocate(THFLWW(JS,JY,JX));   THFLWW=0._r8
  allocate(THGQRS(JY,JX));      THGQRS=0._r8
  allocate(TCOQRS(JY,JX));      TCOQRS=0._r8
  allocate(TCHQRS(JY,JX));      TCHQRS=0._r8
  allocate(TOXQRS(JY,JX));      TOXQRS=0._r8
  allocate(TNGQRS(JY,JX));      TNGQRS=0._r8
  allocate(TN2QRS(JY,JX));      TN2QRS=0._r8
  allocate(TN4QRS(JY,JX));      TN4QRS=0._r8
  allocate(TN3QRS(JY,JX));      TN3QRS=0._r8
  allocate(TNOQRS(JY,JX));      TNOQRS=0._r8
  allocate(TPOQRS(JY,JX));      TPOQRS=0._r8
  allocate(TNXQRS(JY,JX));      TNXQRS=0._r8
  allocate(TQRAL(JY,JX));       TQRAL=0._r8
  allocate(TQRFE(JY,JX));       TQRFE=0._r8
  allocate(TQRHY(JY,JX));       TQRHY=0._r8
  allocate(TQRCA(JY,JX));       TQRCA=0._r8
  allocate(TQRMG(JY,JX));       TQRMG=0._r8
  allocate(TQRNA(JY,JX));       TQRNA=0._r8
  allocate(TQRKA(JY,JX));       TQRKA=0._r8
  allocate(TQROH(JY,JX));       TQROH=0._r8
  allocate(TQRSO(JY,JX));       TQRSO=0._r8
  allocate(TQRCL(JY,JX));       TQRCL=0._r8
  allocate(TQRC3(JY,JX));       TQRC3=0._r8
  allocate(TQRHC(JY,JX));       TQRHC=0._r8
  allocate(TQRAL1(JY,JX));      TQRAL1=0._r8
  allocate(TQRAL2(JY,JX));      TQRAL2=0._r8
  allocate(TQRAL3(JY,JX));      TQRAL3=0._r8
  allocate(TQRAL4(JY,JX));      TQRAL4=0._r8
  allocate(TQRALS(JY,JX));      TQRALS=0._r8
  allocate(TQRFE1(JY,JX));      TQRFE1=0._r8
  allocate(TQRFE2(JY,JX));      TQRFE2=0._r8
  allocate(TQRFE3(JY,JX));      TQRFE3=0._r8
  allocate(TQRFE4(JY,JX));      TQRFE4=0._r8
  allocate(TQRFES(JY,JX));      TQRFES=0._r8
  allocate(TQRCAO(JY,JX));      TQRCAO=0._r8
  allocate(TQRCAC(JY,JX));      TQRCAC=0._r8
  allocate(TQRCAH(JY,JX));      TQRCAH=0._r8
  allocate(TQRCAS(JY,JX));      TQRCAS=0._r8
  allocate(TQRMGO(JY,JX));      TQRMGO=0._r8
  allocate(TQRMGC(JY,JX));      TQRMGC=0._r8
  allocate(TQRMGH(JY,JX));      TQRMGH=0._r8
  allocate(TQRMGS(JY,JX));      TQRMGS=0._r8
  allocate(TQRNAC(JY,JX));      TQRNAC=0._r8
  allocate(TQRNAS(JY,JX));      TQRNAS=0._r8
  allocate(TQRKAS(JY,JX));      TQRKAS=0._r8
  allocate(TQRH0P(JY,JX));      TQRH0P=0._r8
  allocate(TQRH3P(JY,JX));      TQRH3P=0._r8
  allocate(TQRF1P(JY,JX));      TQRF1P=0._r8
  allocate(TP1QRS(JY,JX));      TP1QRS=0._r8
  allocate(TQRF2P(JY,JX));      TQRF2P=0._r8
  allocate(TQRC0P(JY,JX));      TQRC0P=0._r8
  allocate(TQRC1P(JY,JX));      TQRC1P=0._r8
  allocate(TQRC2P(JY,JX));      TQRC2P=0._r8
  allocate(TQRM1P(JY,JX));      TQRM1P=0._r8
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
  allocate(TQSAL(JY,JX));       TQSAL=0._r8
  allocate(TQSFE(JY,JX));       TQSFE=0._r8
  allocate(TQSHY(JY,JX));       TQSHY=0._r8
  allocate(TQSCA(JY,JX));       TQSCA=0._r8
  allocate(TQSMG(JY,JX));       TQSMG=0._r8
  allocate(TQSNA(JY,JX));       TQSNA=0._r8
  allocate(TQSKA(JY,JX));       TQSKA=0._r8
  allocate(TQSOH(JY,JX));       TQSOH=0._r8
  allocate(TQSSO(JY,JX));       TQSSO=0._r8
  allocate(TQSCL(JY,JX));       TQSCL=0._r8
  allocate(TQSC3(JY,JX));       TQSC3=0._r8
  allocate(TQSHC(JY,JX));       TQSHC=0._r8
  allocate(TQSAL1(JY,JX));      TQSAL1=0._r8
  allocate(TQSAL2(JY,JX));      TQSAL2=0._r8
  allocate(TQSAL3(JY,JX));      TQSAL3=0._r8
  allocate(TQSAL4(JY,JX));      TQSAL4=0._r8
  allocate(TQSALS(JY,JX));      TQSALS=0._r8
  allocate(TQSFE1(JY,JX));      TQSFE1=0._r8
  allocate(TQSFE2(JY,JX));      TQSFE2=0._r8
  allocate(TQSFE3(JY,JX));      TQSFE3=0._r8
  allocate(TQSFE4(JY,JX));      TQSFE4=0._r8
  allocate(TQSFES(JY,JX));      TQSFES=0._r8
  allocate(TQSCAO(JY,JX));      TQSCAO=0._r8
  allocate(TQSCAC(JY,JX));      TQSCAC=0._r8
  allocate(TQSCAH(JY,JX));      TQSCAH=0._r8
  allocate(TQSCAS(JY,JX));      TQSCAS=0._r8
  allocate(TQSMGO(JY,JX));      TQSMGO=0._r8
  allocate(TQSMGC(JY,JX));      TQSMGC=0._r8
  allocate(TQSMGH(JY,JX));      TQSMGH=0._r8
  allocate(TQSMGS(JY,JX));      TQSMGS=0._r8
  allocate(TQSNAC(JY,JX));      TQSNAC=0._r8
  allocate(TQSNAS(JY,JX));      TQSNAS=0._r8
  allocate(TQSKAS(JY,JX));      TQSKAS=0._r8
  allocate(TQSH0P(JY,JX));      TQSH0P=0._r8
  allocate(TQSH3P(JY,JX));      TQSH3P=0._r8
  allocate(TQSF1P(JY,JX));      TQSF1P=0._r8
  allocate(TQSF2P(JY,JX));      TQSF2P=0._r8
  allocate(TQSC0P(JY,JX));      TQSC0P=0._r8
  allocate(TQSC1P(JY,JX));      TQSC1P=0._r8
  allocate(TQSC2P(JY,JX));      TQSC2P=0._r8
  allocate(TQSM1P(JY,JX));      TQSM1P=0._r8
  allocate(TALFLS(JZ,JY,JX));   TALFLS=0._r8
  allocate(TFEFLS(JZ,JY,JX));   TFEFLS=0._r8
  allocate(TCAFLS(JZ,JY,JX));   TCAFLS=0._r8
  allocate(THYFLS(JZ,JY,JX));   THYFLS=0._r8
  allocate(TMGFLS(JZ,JY,JX));   TMGFLS=0._r8
  allocate(TNAFLS(JZ,JY,JX));   TNAFLS=0._r8
  allocate(TKAFLS(JZ,JY,JX));   TKAFLS=0._r8
  allocate(TOHFLS(JZ,JY,JX));   TOHFLS=0._r8
  allocate(TSOFLS(JZ,JY,JX));   TSOFLS=0._r8
  allocate(TCLFLS(JZ,JY,JX));   TCLFLS=0._r8
  allocate(TC3FLS(JZ,JY,JX));   TC3FLS=0._r8
  allocate(THCFLS(JZ,JY,JX));   THCFLS=0._r8
  allocate(TAL1FS(JZ,JY,JX));   TAL1FS=0._r8
  allocate(TAL2FS(JZ,JY,JX));   TAL2FS=0._r8
  allocate(TAL3FS(JZ,JY,JX));   TAL3FS=0._r8
  allocate(TAL4FS(JZ,JY,JX));   TAL4FS=0._r8
  allocate(TALSFS(JZ,JY,JX));   TALSFS=0._r8
  allocate(TFE1FS(JZ,JY,JX));   TFE1FS=0._r8
  allocate(TFE2FS(JZ,JY,JX));   TFE2FS=0._r8
  allocate(TFE3FS(JZ,JY,JX));   TFE3FS=0._r8
  allocate(TFE4FS(JZ,JY,JX));   TFE4FS=0._r8
  allocate(TFESFS(JZ,JY,JX));   TFESFS=0._r8
  allocate(TCAOFS(JZ,JY,JX));   TCAOFS=0._r8
  allocate(TCACFS(JZ,JY,JX));   TCACFS=0._r8
  allocate(TCAHFS(JZ,JY,JX));   TCAHFS=0._r8
  allocate(TCASFS(JZ,JY,JX));   TCASFS=0._r8
  allocate(TMGOFS(JZ,JY,JX));   TMGOFS=0._r8
  allocate(TMGCFS(JZ,JY,JX));   TMGCFS=0._r8
  allocate(TMGHFS(JZ,JY,JX));   TMGHFS=0._r8
  allocate(TMGSFS(JZ,JY,JX));   TMGSFS=0._r8
  allocate(TNACFS(JZ,JY,JX));   TNACFS=0._r8
  allocate(TNASFS(JZ,JY,JX));   TNASFS=0._r8
  allocate(TKASFS(JZ,JY,JX));   TKASFS=0._r8
  allocate(TH0PFS(JZ,JY,JX));   TH0PFS=0._r8
  allocate(TH3PFS(JZ,JY,JX));   TH3PFS=0._r8
  allocate(TF1PFS(JZ,JY,JX));   TF1PFS=0._r8
  allocate(TF2PFS(JZ,JY,JX));   TF2PFS=0._r8
  allocate(TC0PFS(JZ,JY,JX));   TC0PFS=0._r8
  allocate(TC1PFS(JZ,JY,JX));   TC1PFS=0._r8
  allocate(TC2PFS(JZ,JY,JX));   TC2PFS=0._r8
  allocate(TM1PFS(JZ,JY,JX));   TM1PFS=0._r8
  allocate(TH0BFB(JZ,JY,JX));   TH0BFB=0._r8
  allocate(TH3BFB(JZ,JY,JX));   TH3BFB=0._r8
  allocate(TF1BFB(JZ,JY,JX));   TF1BFB=0._r8
  allocate(TF2BFB(JZ,JY,JX));   TF2BFB=0._r8
  allocate(TC0BFB(JZ,JY,JX));   TC0BFB=0._r8
  allocate(TC1BFB(JZ,JY,JX));   TC1BFB=0._r8
  allocate(TC2BFB(JZ,JY,JX));   TC2BFB=0._r8
  allocate(TM1BFB(JZ,JY,JX));   TM1BFB=0._r8
  allocate(TALFHS(JZ,JY,JX));   TALFHS=0._r8
  allocate(TFEFHS(JZ,JY,JX));   TFEFHS=0._r8
  allocate(THYFHS(JZ,JY,JX));   THYFHS=0._r8
  allocate(TCAFHS(JZ,JY,JX));   TCAFHS=0._r8
  allocate(TMGFHS(JZ,JY,JX));   TMGFHS=0._r8
  allocate(TNAFHS(JZ,JY,JX));   TNAFHS=0._r8
  allocate(TKAFHS(JZ,JY,JX));   TKAFHS=0._r8
  allocate(TOHFHS(JZ,JY,JX));   TOHFHS=0._r8
  allocate(TSOFHS(JZ,JY,JX));   TSOFHS=0._r8
  allocate(TCLFHS(JZ,JY,JX));   TCLFHS=0._r8
  allocate(TC3FHS(JZ,JY,JX));   TC3FHS=0._r8
  allocate(THCFHS(JZ,JY,JX));   THCFHS=0._r8
  allocate(TAL1HS(JZ,JY,JX));   TAL1HS=0._r8
  allocate(TAL2HS(JZ,JY,JX));   TAL2HS=0._r8
  allocate(TAL3HS(JZ,JY,JX));   TAL3HS=0._r8
  allocate(TAL4HS(JZ,JY,JX));   TAL4HS=0._r8
  allocate(TALSHS(JZ,JY,JX));   TALSHS=0._r8
  allocate(TFE1HS(JZ,JY,JX));   TFE1HS=0._r8
  allocate(TFE2HS(JZ,JY,JX));   TFE2HS=0._r8
  allocate(TFE3HS(JZ,JY,JX));   TFE3HS=0._r8
  allocate(TFE4HS(JZ,JY,JX));   TFE4HS=0._r8
  allocate(TFESHS(JZ,JY,JX));   TFESHS=0._r8
  allocate(TCAOHS(JZ,JY,JX));   TCAOHS=0._r8
  allocate(TCACHS(JZ,JY,JX));   TCACHS=0._r8
  allocate(TCAHHS(JZ,JY,JX));   TCAHHS=0._r8
  allocate(TCASHS(JZ,JY,JX));   TCASHS=0._r8
  allocate(TMGOHS(JZ,JY,JX));   TMGOHS=0._r8
  allocate(TMGCHS(JZ,JY,JX));   TMGCHS=0._r8
  allocate(TMGHHS(JZ,JY,JX));   TMGHHS=0._r8
  allocate(TMGSHS(JZ,JY,JX));   TMGSHS=0._r8
  allocate(TNACHS(JZ,JY,JX));   TNACHS=0._r8
  allocate(TNASHS(JZ,JY,JX));   TNASHS=0._r8
  allocate(TKASHS(JZ,JY,JX));   TKASHS=0._r8
  allocate(TH0PHS(JZ,JY,JX));   TH0PHS=0._r8
  allocate(TH3PHS(JZ,JY,JX));   TH3PHS=0._r8
  allocate(TF1PHS(JZ,JY,JX));   TF1PHS=0._r8
  allocate(TF2PHS(JZ,JY,JX));   TF2PHS=0._r8
  allocate(TC0PHS(JZ,JY,JX));   TC0PHS=0._r8
  allocate(TC1PHS(JZ,JY,JX));   TC1PHS=0._r8
  allocate(TC2PHS(JZ,JY,JX));   TC2PHS=0._r8
  allocate(TM1PHS(JZ,JY,JX));   TM1PHS=0._r8
  allocate(TH0BHB(JZ,JY,JX));   TH0BHB=0._r8
  allocate(TH3BHB(JZ,JY,JX));   TH3BHB=0._r8
  allocate(TF1BHB(JZ,JY,JX));   TF1BHB=0._r8
  allocate(TF2BHB(JZ,JY,JX));   TF2BHB=0._r8
  allocate(TC0BHB(JZ,JY,JX));   TC0BHB=0._r8
  allocate(TC1BHB(JZ,JY,JX));   TC1BHB=0._r8
  allocate(TC2BHB(JZ,JY,JX));   TC2BHB=0._r8
  allocate(TM1BHB(JZ,JY,JX));   TM1BHB=0._r8
  allocate(TSANER(JY,JX));      TSANER=0._r8
  allocate(TSILER(JY,JX));      TSILER=0._r8
  allocate(TCLAER(JY,JX));      TCLAER=0._r8
  allocate(TCECER(JY,JX));      TCECER=0._r8
  allocate(TAECER(JY,JX));      TAECER=0._r8
  allocate(TNH4ER(JY,JX));      TNH4ER=0._r8
  allocate(TNH3ER(JY,JX));      TNH3ER=0._r8
  allocate(TNHUER(JY,JX));      TNHUER=0._r8
  allocate(TNO3ER(JY,JX));      TNO3ER=0._r8
  allocate(TNH4EB(JY,JX));      TNH4EB=0._r8
  allocate(TNH3EB(JY,JX));      TNH3EB=0._r8
  allocate(TNHUEB(JY,JX));      TNHUEB=0._r8
  allocate(TNO3EB(JY,JX));      TNO3EB=0._r8
  allocate(TN4ER(JY,JX));       TN4ER=0._r8
  allocate(TNBER(JY,JX));       TNBER=0._r8
  allocate(THYER(JY,JX));       THYER=0._r8
  allocate(TALER(JY,JX));       TALER=0._r8
  allocate(TCAER(JY,JX));       TCAER=0._r8
  allocate(TMGER(JY,JX));       TMGER=0._r8
  allocate(TNAER(JY,JX));       TNAER=0._r8
  allocate(TKAER(JY,JX));       TKAER=0._r8
  allocate(THCER(JY,JX));       THCER=0._r8
  allocate(TAL2ER(JY,JX));      TAL2ER=0._r8
  allocate(TOH0ER(JY,JX));      TOH0ER=0._r8
  allocate(TOH1ER(JY,JX));      TOH1ER=0._r8
  allocate(TOH2ER(JY,JX));      TOH2ER=0._r8
  allocate(TH1PER(JY,JX));      TH1PER=0._r8
  allocate(TH2PER(JY,JX));      TH2PER=0._r8
  allocate(TOH0EB(JY,JX));      TOH0EB=0._r8
  allocate(TOH1EB(JY,JX));      TOH1EB=0._r8
  allocate(TOH2EB(JY,JX));      TOH2EB=0._r8
  allocate(TH1PEB(JY,JX));      TH1PEB=0._r8
  allocate(TH2PEB(JY,JX));      TH2PEB=0._r8
  allocate(TALOER(JY,JX));      TALOER=0._r8
  allocate(TFEOER(JY,JX));      TFEOER=0._r8
  allocate(TCACER(JY,JX));      TCACER=0._r8
  allocate(TCASER(JY,JX));      TCASER=0._r8
  allocate(TALPER(JY,JX));      TALPER=0._r8
  allocate(TFEPER(JY,JX));      TFEPER=0._r8
  allocate(TCPDER(JY,JX));      TCPDER=0._r8
  allocate(TCPHER(JY,JX));      TCPHER=0._r8
  allocate(TCPMER(JY,JX));      TCPMER=0._r8
  allocate(TALPEB(JY,JX));      TALPEB=0._r8
  allocate(TFEPEB(JY,JX));      TFEPEB=0._r8
  allocate(TCPDEB(JY,JX));      TCPDEB=0._r8
  allocate(TCPHEB(JY,JX));      TCPHEB=0._r8
  allocate(TCPMEB(JY,JX));      TCPMEB=0._r8
  allocate(TFEER(JY,JX));       TFEER=0._r8
  allocate(TFE2ER(JY,JX));      TFE2ER=0._r8
  allocate(TSEDER(JY,JX));      TSEDER=0._r8
  allocate(TFLW(JZ,JY,JX));     TFLW=0._r8
  allocate(TFLWX(JZ,JY,JX));    TFLWX=0._r8
  allocate(THFLW(JZ,JY,JX));    THFLW=0._r8
  allocate(TFLWH(JZ,JY,JX));    TFLWH=0._r8
  allocate(TCOFLS(JZ,JY,JX));   TCOFLS=0._r8
  allocate(TCHFLS(JZ,JY,JX));   TCHFLS=0._r8
  allocate(TOXFLS(JZ,JY,JX));   TOXFLS=0._r8
  allocate(TNXFLB(JZ,JY,JX));   TNXFLB=0._r8
  allocate(TNGFLS(JZ,JY,JX));   TNGFLS=0._r8
  allocate(TN2FLS(JZ,JY,JX));   TN2FLS=0._r8
  allocate(TN4FLS(JZ,JY,JX));   TN4FLS=0._r8
  allocate(TN4FLB(JZ,JY,JX));   TN4FLB=0._r8
  allocate(TN3FLS(JZ,JY,JX));   TN3FLS=0._r8
  allocate(TN3FLB(JZ,JY,JX));   TN3FLB=0._r8
  allocate(TNOFLS(JZ,JY,JX));   TNOFLS=0._r8
  allocate(TNOFLB(JZ,JY,JX));   TNOFLB=0._r8
  allocate(TPOFLS(JZ,JY,JX));   TPOFLS=0._r8
  allocate(TH2BFB(JZ,JY,JX));   TH2BFB=0._r8
  allocate(TNXFLS(JZ,JY,JX));   TNXFLS=0._r8
  allocate(TCOFHS(JZ,JY,JX));   TCOFHS=0._r8
  allocate(TCHFHS(JZ,JY,JX));   TCHFHS=0._r8
  allocate(TNXFHB(JZ,JY,JX));   TNXFHB=0._r8
  allocate(TOXFHS(JZ,JY,JX));   TOXFHS=0._r8
  allocate(TNGFHS(JZ,JY,JX));   TNGFHS=0._r8
  allocate(TN2FHS(JZ,JY,JX));   TN2FHS=0._r8
  allocate(TN4FHS(JZ,JY,JX));   TN4FHS=0._r8
  allocate(TN4FHB(JZ,JY,JX));   TN4FHB=0._r8
  allocate(TN3FHS(JZ,JY,JX));   TN3FHS=0._r8
  allocate(TN3FHB(JZ,JY,JX));   TN3FHB=0._r8
  allocate(TNOFHS(JZ,JY,JX));   TNOFHS=0._r8
  allocate(TNOFHB(JZ,JY,JX));   TNOFHB=0._r8
  allocate(TPOFHS(JZ,JY,JX));   TPOFHS=0._r8
  allocate(TH2BHB(JZ,JY,JX));   TH2BHB=0._r8
  allocate(TNXFHS(JZ,JY,JX));   TNXFHS=0._r8
  allocate(TCOFLG(JZ,JY,JX));   TCOFLG=0._r8
  allocate(TCHFLG(JZ,JY,JX));   TCHFLG=0._r8
  allocate(TOXFLG(JZ,JY,JX));   TOXFLG=0._r8
  allocate(TNGFLG(JZ,JY,JX));   TNGFLG=0._r8
  allocate(TN2FLG(JZ,JY,JX));   TN2FLG=0._r8
  allocate(TNHFLG(JZ,JY,JX));   TNHFLG=0._r8
  allocate(TTHAW(JZ,JY,JX));    TTHAW=0._r8
  allocate(THTHAW(JZ,JY,JX));   THTHAW=0._r8
  allocate(TTHAWH(JZ,JY,JX));   TTHAWH=0._r8
  allocate(TP1FLS(JZ,JY,JX));   TP1FLS=0._r8
  allocate(TP1FHS(JZ,JY,JX));   TP1FHS=0._r8
  allocate(TH1BFB(JZ,JY,JX));   TH1BFB=0._r8
  allocate(TH1BHB(JZ,JY,JX));   TH1BHB=0._r8
  allocate(VOLW1(JZ,JY,JX));    VOLW1=0._r8
  allocate(VOLI1(JZ,JY,JX));    VOLI1=0._r8
  allocate(VOLWH1(JZ,JY,JX));   VOLWH1=0._r8
  allocate(VOLIH1(JZ,JY,JX));   VOLIH1=0._r8
  allocate(TOMCER(nlbiomcp,JG,NFGs,0:jcplx1,JY,JX)); TOMCER=0._r8
  allocate(TOMNER(nlbiomcp,JG,NFGs,0:jcplx1,JY,JX)); TOMNER=0._r8
  allocate(TOMPER(nlbiomcp,JG,NFGs,0:jcplx1,JY,JX)); TOMPER=0._r8

  allocate(TOMCERff(nlbiomcp,JG,NFGs,JY,JX));TOMCERff=0._r8
  allocate(TOMNERff(nlbiomcp,JG,NFGs,JY,JX));TOMNERff=0._r8
  allocate(TOMPERff(nlbiomcp,JG,NFGs,JY,JX));TOMPERff=0._r8

  allocate(TOCFLS(0:jcplx1,JZ,JY,JX));TOCFLS=0._r8
  allocate(TONFLS(0:jcplx1,JZ,JY,JX));TONFLS=0._r8
  allocate(TOPFLS(0:jcplx1,JZ,JY,JX));TOPFLS=0._r8
  allocate(TOAFLS(0:jcplx1,JZ,JY,JX));TOAFLS=0._r8
  allocate(TOCFHS(0:jcplx1,JZ,JY,JX));TOCFHS=0._r8
  allocate(TONFHS(0:jcplx1,JZ,JY,JX));TONFHS=0._r8
  allocate(TOPFHS(0:jcplx1,JZ,JY,JX));TOPFHS=0._r8
  allocate(TOAFHS(0:jcplx1,JZ,JY,JX));TOAFHS=0._r8
  allocate(TOCQRS(0:jcplx1,JY,JX));TOCQRS=0._r8
  allocate(TONQRS(0:jcplx1,JY,JX));TONQRS=0._r8
  allocate(TOPQRS(0:jcplx1,JY,JX));TOPQRS=0._r8
  allocate(TOAQRS(0:jcplx1,JY,JX));TOAQRS=0._r8
  allocate(TORCER(ndbiomcp,0:jcplx1,JY,JX));TORCER=0._r8
  allocate(TORNER(ndbiomcp,0:jcplx1,JY,JX));TORNER=0._r8
  allocate(TORPER(ndbiomcp,0:jcplx1,JY,JX));TORPER=0._r8
  allocate(TOHCER(0:jcplx1,JY,JX));TOHCER=0._r8
  allocate(TOHNER(0:jcplx1,JY,JX));TOHNER=0._r8
  allocate(TOHPER(0:jcplx1,JY,JX));TOHPER=0._r8
  allocate(TOHAER(0:jcplx1,JY,JX));TOHAER=0._r8
  allocate(TOSCER(jsken,0:jcplx1,JY,JX));TOSCER=0._r8
  allocate(TOSAER(jsken,0:jcplx1,JY,JX));TOSAER=0._r8
  allocate(TOSNER(jsken,0:jcplx1,JY,JX));TOSNER=0._r8
  allocate(TOSPER(jsken,0:jcplx1,JY,JX));TOSPER=0._r8
  end subroutine InitTflxType


!------------------------------------------------------------------------------------------

  subroutine DestructTflxType

  implicit none

  call destroy(TOMCER)
  call destroy(TOMNER)
  call destroy(TOMPER)
  call destroy(TOMCERff)
  call destroy(TOMNERff)
  call destroy(TOMPERff)
  call destroy(TOCFLS)
  call destroy(TONFLS)
  call destroy(TOPFLS)
  call destroy(TOAFLS)
  call destroy(TOCFHS)
  call destroy(TONFHS)
  call destroy(TOPFHS)
  call destroy(TOAFHS)
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
  call destroy(TALBLS)
  call destroy(TFEBLS)
  call destroy(THYBLS)
  call destroy(TCABLS)
  call destroy(TMGBLS)
  call destroy(TNABLS)
  call destroy(TKABLS)
  call destroy(TOHBLS)
  call destroy(TSOBLS)
  call destroy(TCLBLS)
  call destroy(TC3BLS)
  call destroy(THCBLS)
  call destroy(TAL1BS)
  call destroy(TAL2BS)
  call destroy(TAL3BS)
  call destroy(TAL4BS)
  call destroy(TALSBS)
  call destroy(TFE1BS)
  call destroy(TFE2BS)
  call destroy(TFE3BS)
  call destroy(TFE4BS)
  call destroy(TFESBS)
  call destroy(TCAOBS)
  call destroy(TCACBS)
  call destroy(TCAHBS)
  call destroy(TCASBS)
  call destroy(TMGOBS)
  call destroy(TMGCBS)
  call destroy(TMGHBS)
  call destroy(TMGSBS)
  call destroy(TNACBS)
  call destroy(TNASBS)
  call destroy(TKASBS)
  call destroy(TH0PBS)
  call destroy(TH3PBS)
  call destroy(TF1PBS)
  call destroy(TF2PBS)
  call destroy(TC0PBS)
  call destroy(TC1PBS)
  call destroy(TC2PBS)
  call destroy(TM1PBS)
  call destroy(THGFLG)
  call destroy(THGFLS)
  call destroy(THGFHS)
  call destroy(TQR)
  call destroy(THQR)
  call destroy(TQS)
  call destroy(TQW)
  call destroy(TQI)
  call destroy(THQS)
  call destroy(TFLWS)
  call destroy(TFLWW)
  call destroy(TFLWI)
  call destroy(THFLWW)
  call destroy(THGQRS)
  call destroy(TCOQRS)
  call destroy(TCHQRS)
  call destroy(TOXQRS)
  call destroy(TNGQRS)
  call destroy(TN2QRS)
  call destroy(TN4QRS)
  call destroy(TN3QRS)
  call destroy(TNOQRS)
  call destroy(TPOQRS)
  call destroy(TNXQRS)
  call destroy(TQRAL)
  call destroy(TQRFE)
  call destroy(TQRHY)
  call destroy(TQRCA)
  call destroy(TQRMG)
  call destroy(TQRNA)
  call destroy(TQRKA)
  call destroy(TQROH)
  call destroy(TQRSO)
  call destroy(TQRCL)
  call destroy(TQRC3)
  call destroy(TQRHC)
  call destroy(TQRAL1)
  call destroy(TQRAL2)
  call destroy(TQRAL3)
  call destroy(TQRAL4)
  call destroy(TQRALS)
  call destroy(TQRFE1)
  call destroy(TQRFE2)
  call destroy(TQRFE3)
  call destroy(TQRFE4)
  call destroy(TQRFES)
  call destroy(TQRCAO)
  call destroy(TQRCAC)
  call destroy(TQRCAH)
  call destroy(TQRCAS)
  call destroy(TQRMGO)
  call destroy(TQRMGC)
  call destroy(TQRMGH)
  call destroy(TQRMGS)
  call destroy(TQRNAC)
  call destroy(TQRNAS)
  call destroy(TQRKAS)
  call destroy(TQRH0P)
  call destroy(TQRH3P)
  call destroy(TQRF1P)
  call destroy(TP1QRS)
  call destroy(TQRF2P)
  call destroy(TQRC0P)
  call destroy(TQRC1P)
  call destroy(TQRC2P)
  call destroy(TQRM1P)
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
  call destroy(TQSAL)
  call destroy(TQSFE)
  call destroy(TQSHY)
  call destroy(TQSCA)
  call destroy(TQSMG)
  call destroy(TQSNA)
  call destroy(TQSKA)
  call destroy(TQSOH)
  call destroy(TQSSO)
  call destroy(TQSCL)
  call destroy(TQSC3)
  call destroy(TQSHC)
  call destroy(TQSAL1)
  call destroy(TQSAL2)
  call destroy(TQSAL3)
  call destroy(TQSAL4)
  call destroy(TQSALS)
  call destroy(TQSFE1)
  call destroy(TQSFE2)
  call destroy(TQSFE3)
  call destroy(TQSFE4)
  call destroy(TQSFES)
  call destroy(TQSCAO)
  call destroy(TQSCAC)
  call destroy(TQSCAH)
  call destroy(TQSCAS)
  call destroy(TQSMGO)
  call destroy(TQSMGC)
  call destroy(TQSMGH)
  call destroy(TQSMGS)
  call destroy(TQSNAC)
  call destroy(TQSNAS)
  call destroy(TQSKAS)
  call destroy(TQSH0P)
  call destroy(TQSH3P)
  call destroy(TQSF1P)
  call destroy(TQSF2P)
  call destroy(TQSC0P)
  call destroy(TQSC1P)
  call destroy(TQSC2P)
  call destroy(TQSM1P)
  call destroy(TALFLS)
  call destroy(TFEFLS)
  call destroy(TCAFLS)
  call destroy(THYFLS)
  call destroy(TMGFLS)
  call destroy(TNAFLS)
  call destroy(TKAFLS)
  call destroy(TOHFLS)
  call destroy(TSOFLS)
  call destroy(TCLFLS)
  call destroy(TC3FLS)
  call destroy(THCFLS)
  call destroy(TAL1FS)
  call destroy(TAL2FS)
  call destroy(TAL3FS)
  call destroy(TAL4FS)
  call destroy(TALSFS)
  call destroy(TFE1FS)
  call destroy(TFE2FS)
  call destroy(TFE3FS)
  call destroy(TFE4FS)
  call destroy(TFESFS)
  call destroy(TCAOFS)
  call destroy(TCACFS)
  call destroy(TCAHFS)
  call destroy(TCASFS)
  call destroy(TMGOFS)
  call destroy(TMGCFS)
  call destroy(TMGHFS)
  call destroy(TMGSFS)
  call destroy(TNACFS)
  call destroy(TNASFS)
  call destroy(TKASFS)
  call destroy(TH0PFS)
  call destroy(TH3PFS)
  call destroy(TF1PFS)
  call destroy(TF2PFS)
  call destroy(TC0PFS)
  call destroy(TC1PFS)
  call destroy(TC2PFS)
  call destroy(TM1PFS)
  call destroy(TH0BFB)
  call destroy(TH3BFB)
  call destroy(TF1BFB)
  call destroy(TF2BFB)
  call destroy(TC0BFB)
  call destroy(TC1BFB)
  call destroy(TC2BFB)
  call destroy(TM1BFB)
  call destroy(TALFHS)
  call destroy(TFEFHS)
  call destroy(THYFHS)
  call destroy(TCAFHS)
  call destroy(TMGFHS)
  call destroy(TNAFHS)
  call destroy(TKAFHS)
  call destroy(TOHFHS)
  call destroy(TSOFHS)
  call destroy(TCLFHS)
  call destroy(TC3FHS)
  call destroy(THCFHS)
  call destroy(TAL1HS)
  call destroy(TAL2HS)
  call destroy(TAL3HS)
  call destroy(TAL4HS)
  call destroy(TALSHS)
  call destroy(TFE1HS)
  call destroy(TFE2HS)
  call destroy(TFE3HS)
  call destroy(TFE4HS)
  call destroy(TFESHS)
  call destroy(TCAOHS)
  call destroy(TCACHS)
  call destroy(TCAHHS)
  call destroy(TCASHS)
  call destroy(TMGOHS)
  call destroy(TMGCHS)
  call destroy(TMGHHS)
  call destroy(TMGSHS)
  call destroy(TNACHS)
  call destroy(TNASHS)
  call destroy(TKASHS)
  call destroy(TH0PHS)
  call destroy(TH3PHS)
  call destroy(TF1PHS)
  call destroy(TF2PHS)
  call destroy(TC0PHS)
  call destroy(TC1PHS)
  call destroy(TC2PHS)
  call destroy(TM1PHS)
  call destroy(TH0BHB)
  call destroy(TH3BHB)
  call destroy(TF1BHB)
  call destroy(TF2BHB)
  call destroy(TC0BHB)
  call destroy(TC1BHB)
  call destroy(TC2BHB)
  call destroy(TM1BHB)
  call destroy(TSANER)
  call destroy(TSILER)
  call destroy(TCLAER)
  call destroy(TCECER)
  call destroy(TAECER)
  call destroy(TNH4ER)
  call destroy(TNH3ER)
  call destroy(TNHUER)
  call destroy(TNO3ER)
  call destroy(TNH4EB)
  call destroy(TNH3EB)
  call destroy(TNHUEB)
  call destroy(TNO3EB)
  call destroy(TN4ER)
  call destroy(TNBER)
  call destroy(THYER)
  call destroy(TALER)
  call destroy(TCAER)
  call destroy(TMGER)
  call destroy(TNAER)
  call destroy(TKAER)
  call destroy(THCER)
  call destroy(TAL2ER)
  call destroy(TOH0ER)
  call destroy(TOH1ER)
  call destroy(TOH2ER)
  call destroy(TH1PER)
  call destroy(TH2PER)
  call destroy(TOH0EB)
  call destroy(TOH1EB)
  call destroy(TOH2EB)
  call destroy(TH1PEB)
  call destroy(TH2PEB)
  call destroy(TALOER)
  call destroy(TFEOER)
  call destroy(TCACER)
  call destroy(TCASER)
  call destroy(TALPER)
  call destroy(TFEPER)
  call destroy(TCPDER)
  call destroy(TCPHER)
  call destroy(TCPMER)
  call destroy(TALPEB)
  call destroy(TFEPEB)
  call destroy(TCPDEB)
  call destroy(TCPHEB)
  call destroy(TCPMEB)
  call destroy(TFEER)
  call destroy(TFE2ER)
  call destroy(TSEDER)
  call destroy(TFLW)
  call destroy(TFLWX)
  call destroy(THFLW)
  call destroy(TFLWH)
  call destroy(TCOFLS)
  call destroy(TCHFLS)
  call destroy(TOXFLS)
  call destroy(TNXFLB)
  call destroy(TNGFLS)
  call destroy(TN2FLS)
  call destroy(TN4FLS)
  call destroy(TN4FLB)
  call destroy(TN3FLS)
  call destroy(TN3FLB)
  call destroy(TNOFLS)
  call destroy(TNOFLB)
  call destroy(TPOFLS)
  call destroy(TH2BFB)
  call destroy(TNXFLS)
  call destroy(TCOFHS)
  call destroy(TCHFHS)
  call destroy(TNXFHB)
  call destroy(TOXFHS)
  call destroy(TNGFHS)
  call destroy(TN2FHS)
  call destroy(TN4FHS)
  call destroy(TN4FHB)
  call destroy(TN3FHS)
  call destroy(TN3FHB)
  call destroy(TNOFHS)
  call destroy(TNOFHB)
  call destroy(TPOFHS)
  call destroy(TH2BHB)
  call destroy(TNXFHS)
  call destroy(TCOFLG)
  call destroy(TCHFLG)
  call destroy(TOXFLG)
  call destroy(TNGFLG)
  call destroy(TN2FLG)
  call destroy(TNHFLG)
  call destroy(TTHAW)
  call destroy(THTHAW)
  call destroy(TTHAWH)
  call destroy(TP1FLS)
  call destroy(TP1FHS)
  call destroy(TH1BFB)
  call destroy(TH1BHB)
  call destroy(VOLW1)
  call destroy(VOLI1)
  call destroy(VOLWH1)
  call destroy(VOLIH1)

  end subroutine DestructTflxType
end module TFlxTypeMod
