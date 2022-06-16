module TFlxTypeMod

  use GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
implicit none

  character(len=*), private, parameter :: mod_filename = __FILE__

  public
  real(r8) :: TCOBLS(JS,JY,JX),TCHBLS(JS,JY,JX),TOXBLS(JS,JY,JX) &
    ,TNGBLS(JS,JY,JX),TN2BLS(JS,JY,JX),TN4BLW(JS,JY,JX) &
    ,TN3BLW(JS,JY,JX),TNOBLW(JS,JY,JX),TH1PBS(JS,JY,JX) &
    ,TH2PBS(JS,JY,JX),TALBLS(JS,JY,JX),TFEBLS(JS,JY,JX) &
    ,THYBLS(JS,JY,JX),TCABLS(JS,JY,JX),TMGBLS(JS,JY,JX) &
    ,TNABLS(JS,JY,JX),TKABLS(JS,JY,JX),TOHBLS(JS,JY,JX) &
    ,TSOBLS(JS,JY,JX),TCLBLS(JS,JY,JX),TC3BLS(JS,JY,JX) &
    ,THCBLS(JS,JY,JX),TAL1BS(JS,JY,JX),TAL2BS(JS,JY,JX) &
    ,TAL3BS(JS,JY,JX),TAL4BS(JS,JY,JX),TALSBS(JS,JY,JX) &
    ,TFE1BS(JS,JY,JX),TFE2BS(JS,JY,JX),TFE3BS(JS,JY,JX) &
    ,TFE4BS(JS,JY,JX),TFESBS(JS,JY,JX),TCAOBS(JS,JY,JX) &
    ,TCACBS(JS,JY,JX),TCAHBS(JS,JY,JX),TCASBS(JS,JY,JX) &
    ,TMGOBS(JS,JY,JX),TMGCBS(JS,JY,JX),TMGHBS(JS,JY,JX) &
    ,TMGSBS(JS,JY,JX),TNACBS(JS,JY,JX),TNASBS(JS,JY,JX) &
    ,TKASBS(JS,JY,JX),TH0PBS(JS,JY,JX),TH3PBS(JS,JY,JX) &
    ,TF1PBS(JS,JY,JX),TF2PBS(JS,JY,JX),TC0PBS(JS,JY,JX) &
    ,TC1PBS(JS,JY,JX),TC2PBS(JS,JY,JX),TM1PBS(JS,JY,JX) &
    ,THGFLG(JZ,JY,JX),THGFLS(JZ,JY,JX),THGFHS(JZ,JY,JX)

!   real(r8) :: TQRH1P(JY,JX),TQSH1P(JY,JX),TH1PFS(JZ,JY,JX)

  real(r8) :: TQR(JY,JX),THQR(JY,JX),TQS(JY,JX),TQW(JY,JX) &
    ,TQI(JY,JX),THQS(JY,JX),TFLWS(JS,JY,JX),TFLWW(JS,JY,JX) &
    ,TFLWI(JS,JY,JX),THFLWW(JS,JY,JX),THGQRS(JY,JX) &
    ,TCOQRS(JY,JX),TCHQRS(JY,JX),TOXQRS(JY,JX) &
    ,TNGQRS(JY,JX),TN2QRS(JY,JX),TN4QRS(JY,JX),TN3QRS(JY,JX) &
    ,TNOQRS(JY,JX),TPOQRS(JY,JX),TNXQRS(JY,JX),TQRAL(JY,JX) &
    ,TQRFE(JY,JX),TQRHY(JY,JX),TQRCA(JY,JX),TQRMG(JY,JX) &
    ,TQRNA(JY,JX),TQRKA(JY,JX),TQROH(JY,JX),TQRSO(JY,JX) &
    ,TQRCL(JY,JX),TQRC3(JY,JX),TQRHC(JY,JX),TQRAL1(JY,JX) &
    ,TQRAL2(JY,JX),TQRAL3(JY,JX),TQRAL4(JY,JX),TQRALS(JY,JX) &
    ,TQRFE1(JY,JX),TQRFE2(JY,JX),TQRFE3(JY,JX),TQRFE4(JY,JX) &
    ,TQRFES(JY,JX),TQRCAO(JY,JX),TQRCAC(JY,JX),TQRCAH(JY,JX) &
    ,TQRCAS(JY,JX),TQRMGO(JY,JX),TQRMGC(JY,JX),TQRMGH(JY,JX) &
    ,TQRMGS(JY,JX),TQRNAC(JY,JX),TQRNAS(JY,JX),TQRKAS(JY,JX) &
    ,TQRH0P(JY,JX),TQRH3P(JY,JX),TQRF1P(JY,JX),TP1QRS(JY,JX) &
    ,TQRF2P(JY,JX),TQRC0P(JY,JX),TQRC1P(JY,JX),TQRC2P(JY,JX) &
    ,TQRM1P(JY,JX),TCOQSS(JY,JX),TCHQSS(JY,JX),TOXQSS(JY,JX) &
    ,TNGQSS(JY,JX),TN2QSS(JY,JX),TN4QSS(JY,JX),TN3QSS(JY,JX) &
    ,TNOQSS(JY,JX),TPOQSS(JY,JX),TP1QSS(JY,JX),TQSAL(JY,JX) &
    ,TQSFE(JY,JX),TQSHY(JY,JX),TQSCA(JY,JX),TQSMG(JY,JX) &
    ,TQSNA(JY,JX),TQSKA(JY,JX),TQSOH(JY,JX),TQSSO(JY,JX) &
    ,TQSCL(JY,JX),TQSC3(JY,JX),TQSHC(JY,JX),TQSAL1(JY,JX) &
    ,TQSAL2(JY,JX),TQSAL3(JY,JX),TQSAL4(JY,JX),TQSALS(JY,JX) &
    ,TQSFE1(JY,JX),TQSFE2(JY,JX),TQSFE3(JY,JX),TQSFE4(JY,JX) &
    ,TQSFES(JY,JX),TQSCAO(JY,JX),TQSCAC(JY,JX),TQSCAH(JY,JX) &
    ,TQSCAS(JY,JX),TQSMGO(JY,JX),TQSMGC(JY,JX),TQSMGH(JY,JX) &
    ,TQSMGS(JY,JX),TQSNAC(JY,JX),TQSNAS(JY,JX),TQSKAS(JY,JX) &
    ,TQSH0P(JY,JX),TQSH3P(JY,JX),TQSF1P(JY,JX) &
    ,TQSF2P(JY,JX),TQSC0P(JY,JX),TQSC1P(JY,JX),TQSC2P(JY,JX) &
    ,TQSM1P(JY,JX)

  real(r8) :: TALFLS(JZ,JY,JX),TFEFLS(JZ,JY,JX) &
    ,TCAFLS(JZ,JY,JX),THYFLS(JZ,JY,JX),TMGFLS(JZ,JY,JX) &
    ,TNAFLS(JZ,JY,JX),TKAFLS(JZ,JY,JX),TOHFLS(JZ,JY,JX) &
    ,TSOFLS(JZ,JY,JX),TCLFLS(JZ,JY,JX),TC3FLS(JZ,JY,JX) &
    ,THCFLS(JZ,JY,JX),TAL1FS(JZ,JY,JX),TAL2FS(JZ,JY,JX) &
    ,TAL3FS(JZ,JY,JX),TAL4FS(JZ,JY,JX),TALSFS(JZ,JY,JX) &
    ,TFE1FS(JZ,JY,JX),TFE2FS(JZ,JY,JX) &
    ,TFE3FS(JZ,JY,JX),TFE4FS(JZ,JY,JX),TFESFS(JZ,JY,JX) &
    ,TCAOFS(JZ,JY,JX),TCACFS(JZ,JY,JX),TCAHFS(JZ,JY,JX) &
    ,TCASFS(JZ,JY,JX),TMGOFS(JZ,JY,JX),TMGCFS(JZ,JY,JX) &
    ,TMGHFS(JZ,JY,JX),TMGSFS(JZ,JY,JX),TNACFS(JZ,JY,JX) &
    ,TNASFS(JZ,JY,JX),TKASFS(JZ,JY,JX),TH0PFS(JZ,JY,JX) &
    ,TH3PFS(JZ,JY,JX),TF1PFS(JZ,JY,JX) &
    ,TF2PFS(JZ,JY,JX),TC0PFS(JZ,JY,JX),TC1PFS(JZ,JY,JX) &
    ,TC2PFS(JZ,JY,JX),TM1PFS(JZ,JY,JX),TH0BFB(JZ,JY,JX) &
    ,TH3BFB(JZ,JY,JX),TF1BFB(JZ,JY,JX) &
    ,TF2BFB(JZ,JY,JX),TC0BFB(JZ,JY,JX),TC1BFB(JZ,JY,JX) &
    ,TC2BFB(JZ,JY,JX),TM1BFB(JZ,JY,JX)

  real(r8) :: TALFHS(JZ,JY,JX),TFEFHS(JZ,JY,JX) &
    ,THYFHS(JZ,JY,JX),TCAFHS(JZ,JY,JX),TMGFHS(JZ,JY,JX) &
    ,TNAFHS(JZ,JY,JX),TKAFHS(JZ,JY,JX),TOHFHS(JZ,JY,JX) &
    ,TSOFHS(JZ,JY,JX),TCLFHS(JZ,JY,JX),TC3FHS(JZ,JY,JX) &
    ,THCFHS(JZ,JY,JX),TAL1HS(JZ,JY,JX),TAL2HS(JZ,JY,JX) &
    ,TAL3HS(JZ,JY,JX),TAL4HS(JZ,JY,JX),TALSHS(JZ,JY,JX) &
    ,TFE1HS(JZ,JY,JX),TFE2HS(JZ,JY,JX) &
    ,TFE3HS(JZ,JY,JX),TFE4HS(JZ,JY,JX),TFESHS(JZ,JY,JX) &
    ,TCAOHS(JZ,JY,JX),TCACHS(JZ,JY,JX),TCAHHS(JZ,JY,JX) &
    ,TCASHS(JZ,JY,JX),TMGOHS(JZ,JY,JX),TMGCHS(JZ,JY,JX) &
    ,TMGHHS(JZ,JY,JX),TMGSHS(JZ,JY,JX),TNACHS(JZ,JY,JX) &
    ,TNASHS(JZ,JY,JX),TKASHS(JZ,JY,JX),TH0PHS(JZ,JY,JX) &
    ,TH3PHS(JZ,JY,JX),TF1PHS(JZ,JY,JX) &
    ,TF2PHS(JZ,JY,JX),TC0PHS(JZ,JY,JX),TC1PHS(JZ,JY,JX) &
    ,TC2PHS(JZ,JY,JX),TM1PHS(JZ,JY,JX),TH0BHB(JZ,JY,JX) &
    ,TH3BHB(JZ,JY,JX),TF1BHB(JZ,JY,JX) &
    ,TF2BHB(JZ,JY,JX),TC0BHB(JZ,JY,JX),TC1BHB(JZ,JY,JX) &
    ,TC2BHB(JZ,JY,JX),TM1BHB(JZ,JY,JX)

  real(r8) :: TSANER(JY,JX),TSILER(JY,JX),TCLAER(JY,JX) &
    ,TCECER(JY,JX),TAECER(JY,JX),TNH4ER(JY,JX),TNH3ER(JY,JX) &
    ,TNHUER(JY,JX),TNO3ER(JY,JX),TNH4EB(JY,JX),TNH3EB(JY,JX) &
    ,TNHUEB(JY,JX),TNO3EB(JY,JX),TN4ER(JY,JX),TNBER(JY,JX) &
    ,THYER(JY,JX),TALER(JY,JX),TCAER(JY,JX),TMGER(JY,JX) &
    ,TNAER(JY,JX),TKAER(JY,JX),THCER(JY,JX),TAL2ER(JY,JX) &
    ,TOH0ER(JY,JX),TOH1ER(JY,JX),TOH2ER(JY,JX),TH1PER(JY,JX) &
    ,TH2PER(JY,JX),TOH0EB(JY,JX),TOH1EB(JY,JX),TOH2EB(JY,JX) &
    ,TH1PEB(JY,JX),TH2PEB(JY,JX),TALOER(JY,JX),TFEOER(JY,JX) &
    ,TCACER(JY,JX),TCASER(JY,JX),TALPER(JY,JX),TFEPER(JY,JX) &
    ,TCPDER(JY,JX),TCPHER(JY,JX),TCPMER(JY,JX),TALPEB(JY,JX) &
    ,TFEPEB(JY,JX),TCPDEB(JY,JX),TCPHEB(JY,JX),TCPMEB(JY,JX) &
    ,TFEER(JY,JX),TFE2ER(JY,JX),TSEDER(JY,JX)

  real(r8) :: TFLW(JZ,JY,JX),TFLWX(JZ,JY,JX),THFLW(JZ,JY,JX)
  real(r8) :: TFLWH(JZ,JY,JX)

  real(r8) :: TCOFLS(JZ,JY,JX) &
    ,TCHFLS(JZ,JY,JX),TOXFLS(JZ,JY,JX),TNXFLB(JZ,JY,JX) &
    ,TNGFLS(JZ,JY,JX),TN2FLS(JZ,JY,JX),TN4FLS(JZ,JY,JX) &
    ,TN4FLB(JZ,JY,JX),TN3FLS(JZ,JY,JX),TN3FLB(JZ,JY,JX) &
    ,TNOFLS(JZ,JY,JX),TNOFLB(JZ,JY,JX),TPOFLS(JZ,JY,JX) &
    ,TH2BFB(JZ,JY,JX),TNXFLS(JZ,JY,JX) &
    ,TCOFHS(JZ,JY,JX),TCHFHS(JZ,JY,JX),TNXFHB(JZ,JY,JX) &
    ,TOXFHS(JZ,JY,JX),TNGFHS(JZ,JY,JX),TN2FHS(JZ,JY,JX) &
    ,TN4FHS(JZ,JY,JX),TN4FHB(JZ,JY,JX),TN3FHS(JZ,JY,JX) &
    ,TN3FHB(JZ,JY,JX),TNOFHS(JZ,JY,JX),TNOFHB(JZ,JY,JX) &
    ,TPOFHS(JZ,JY,JX),TH2BHB(JZ,JY,JX),TNXFHS(JZ,JY,JX) &
    ,TCOFLG(JZ,JY,JX),TCHFLG(JZ,JY,JX),TOXFLG(JZ,JY,JX) &
    ,TNGFLG(JZ,JY,JX),TN2FLG(JZ,JY,JX),TNHFLG(JZ,JY,JX) &
    ,TTHAW(JZ,JY,JX),THTHAW(JZ,JY,JX),TTHAWH(JZ,JY,JX) &
    ,TP1FLS(JZ,JY,JX),TP1FHS(JZ,JY,JX),TH1BFB(JZ,JY,JX) &
    ,TH1BHB(JZ,JY,JX)
  real(r8) :: VOLW1(JZ,JY,JX),VOLI1(JZ,JY,JX) &
    ,VOLWH1(JZ,JY,JX),VOLIH1(JZ,JY,JX)

  real(r8),allocatable :: TOMCER(:,:,:,:,:,:)
  real(r8),allocatable :: TOMNER(:,:,:,:,:,:)
  real(r8),allocatable :: TOMPER(:,:,:,:,:,:)
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

  allocate(TOMCER(3,JG,7,0:jcplx1,JY,JX))
  allocate(TOMNER(3,JG,7,0:jcplx1,JY,JX))
  allocate(TOMPER(3,JG,7,0:jcplx1,JY,JX))
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
  allocate(TORCER(2,0:jcplx1,JY,JX));TORCER=0._r8
  allocate(TORNER(2,0:jcplx1,JY,JX));TORNER=0._r8
  allocate(TORPER(2,0:jcplx1,JY,JX));TORPER=0._r8
  allocate(TOHCER(0:jcplx1,JY,JX));TOHCER=0._r8
  allocate(TOHNER(0:jcplx1,JY,JX));TOHNER=0._r8
  allocate(TOHPER(0:jcplx1,JY,JX));TOHPER=0._r8
  allocate(TOHAER(0:jcplx1,JY,JX));TOHAER=0._r8
  allocate(TOSCER(jsken,0:jcplx1,JY,JX));TOSCER=0._r8
  allocate(TOSAER(4,0:jcplx1,JY,JX));TOSAER=0._r8
  allocate(TOSNER(4,0:jcplx1,JY,JX));TOSNER=0._r8
  allocate(TOSPER(4,0:jcplx1,JY,JX));TOSPER=0._r8
  end subroutine InitTflxType


!------------------------------------------------------------------------------------------

  subroutine DestructTflxType

  implicit none

  call destroy(TOMCER)
  call destroy(TOMNER)
  call destroy(TOMPER)
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


  end subroutine DestructTflxType
end module TFlxTypeMod
