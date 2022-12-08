module TFlxTypeMod

  use GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc,NFGs=>NFGsc
  use EcoSIMConfig, only : nlbiomcp=>nlbiomcpc,ndbiomcp=>ndbiomcpc
implicit none

  character(len=*), private, parameter :: mod_filename = __FILE__

  public
  real(r8),allocatable ::  trcg_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcn_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcsa_TBLS(:,:,:,:)                      !

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
  real(r8),allocatable ::  TQRM1P(:,:)
                    !
  real(r8),allocatable ::  trcg_QSS(:,:,:)
  real(r8),allocatable ::  trcn_QSS(:,:,:)
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

  real(r8),allocatable ::  trcsa_TFLS(:,:,:,:)                      !
  real(r8),allocatable ::  trcsa_TFHS(:,:,:,:)                      !

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

  real(r8),allocatable :: TOMCER(:,:,:,:,:)
  real(r8),allocatable :: TOMNER(:,:,:,:,:)
  real(r8),allocatable :: TOMPER(:,:,:,:,:)

  real(r8),allocatable :: trcsa_TQS(:,:,:)
  real(r8),allocatable :: TOMCERff(:,:,:,:)
  real(r8),allocatable :: TOMNERff(:,:,:,:)
  real(r8),allocatable :: TOMPERff(:,:,:,:)

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

  allocate(trcg_TBLS(idg_beg:idg_end-1,JS,JY,JX)); trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_TBLS=0._r8
  allocate(trcsa_TBLS(idsa_beg:idsa_end,JS,JY,JX));   trcsa_TBLS=0._r8

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

  allocate(trcg_QSS(idg_beg:idg_end-1,JY,JX));trcg_QSS=0._r8
  allocate(trcn_QSS(ids_nut_beg:ids_nuts_end,JY,JX));trcn_QSS=0._r8

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

  allocate(trcsa_TFLS(idsa_beg:idsab_end,JZ,JY,JX)); trcsa_TFLS=0._r8
  allocate(trcsa_TFHS(idsa_beg:idsab_end,JZ,JY,JX)); trcsa_TFHS=0._r8

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
  allocate(TOMCER(nlbiomcp,NMICBSO,1:jcplx,JY,JX)); TOMCER=0._r8
  allocate(TOMNER(nlbiomcp,NMICBSO,1:jcplx,JY,JX)); TOMNER=0._r8
  allocate(TOMPER(nlbiomcp,NMICBSO,1:jcplx,JY,JX)); TOMPER=0._r8

  allocate(TOMCERff(nlbiomcp,NMICBSA,JY,JX));TOMCERff=0._r8
  allocate(TOMNERff(nlbiomcp,NMICBSA,JY,JX));TOMNERff=0._r8
  allocate(TOMPERff(nlbiomcp,NMICBSA,JY,JX));TOMPERff=0._r8

  allocate(TOCFLS(1:jcplx,JZ,JY,JX));TOCFLS=0._r8
  allocate(TONFLS(1:jcplx,JZ,JY,JX));TONFLS=0._r8
  allocate(TOPFLS(1:jcplx,JZ,JY,JX));TOPFLS=0._r8
  allocate(TOAFLS(1:jcplx,JZ,JY,JX));TOAFLS=0._r8
  allocate(TOCFHS(1:jcplx,JZ,JY,JX));TOCFHS=0._r8
  allocate(TONFHS(1:jcplx,JZ,JY,JX));TONFHS=0._r8
  allocate(TOPFHS(1:jcplx,JZ,JY,JX));TOPFHS=0._r8
  allocate(TOAFHS(1:jcplx,JZ,JY,JX));TOAFHS=0._r8
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
  allocate(trcsa_TQS(idsa_beg:idsa_end,JY,JX));trcsa_TQS=0._r8
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

  call destroy(trcsa_TQS)
  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)

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

  call destroy(trcg_QSS)
  call destroy(trcn_QSS)
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
