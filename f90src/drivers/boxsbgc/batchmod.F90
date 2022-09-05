module batchmod

  use abortutils     , only : endrun
  use MiscUtilMod    , only : addone
  use data_kind_mod  , only : r8 => SHR_KIND_R8
  use ModelStatusType, only : model_status_type
  use MicBGCPars     , only : micpar
  use MicForcTypeMod , only : micforctype
  use MicFLuxTypeMod , only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use fileUtil
  use MicIDMod
  use ChemIDMod
  use ChemIDMod  , only : getchemvarlist => getvarlist_nosalt
implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: BatchModelConfig
  public :: getvarllen, getvarlist,initmodel
  logical :: Litlayer


contains

  function getvarllen()result(nvars)
  implicit none
  integer :: nvars

  call Initboxbgc(nvars)

  end function getvarllen
! ----------------------------------------------------------------------
  subroutine initmodel(nvars, ystates0l, err_status)
!
! set initial conditions for the boxsbgc

  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystates0l(nvars)
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()


  end subroutine initmodel

! ----------------------------------------------------------------------

  subroutine BatchModelConfig(nvars,ystates0l,forc,micfor,micstt,micflx,err_status)
!
! DESCRIPTION:
! configure the batch mode of the soil bgc
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod      , only : micforctype
  use MicBGCPars          , only : micpar
  use ForcTypeMod         , only : forc_type
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(in) :: ystates0l(nvars)
  type(forc_type), intent(in) :: forc
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(model_status_type), intent(out) :: err_status

  real(r8), parameter :: ZERO=1.0E-15_r8
  real(r8), parameter :: ZEROS2=1.0E-08_r8
  integer :: kk

  call err_status%reset()

  associate(                      &
    nlbiomcp => micpar%nlbiomcp , &
    ndbiomcp => micpar%ndbiomcp , &
    NFGs    => micpar%NFGs      , &
    jcplx1  => micpar%jcplx1    , &
    JG      => micpar%jguilds     &
  )
  micfor%ZERO  =ZERO
  micfor%ZEROS2=ZEROS2
  micfor%ZEROS =ZERO

!  micfor%CCH4E =CCH4E
!  micfor%COXYE =COXYE
!  micfor%COXQ  =COXQ
!  micfor%COXR  =COXR
!  micfor%FLQRI =FLQRI
!  micfor%FLQRQ =FLQRQ
!  micfor%OFFSET=OFFSET
!  micfor%VOLR  =VOLR
!  micfor%VOLWRX=VOLWRX
!  micfor%VOLY  =VOLY
!  micfor%THETY =THETY
!  micfor%POROS =POROS
!  micfor%FC    =FC
!  micfor%TKS   =TKS
!  micfor%THETW =THETW
!  micfor%PH    =PH
!  micfor%BKVL  =BKVL
!  micfor%VOLX  =VOLX
!  micfor%TFND  =TFND
!  micfor%VLNOB =VLNOB
!  micfor%VLNO3 =VLNO3
!  micfor%VLNH4 =VLNH4
!  micfor%VLNHB =VLNHB
!  micfor%VLPO4 =VLPO4
!  micfor%VLPOB =VLPOB
!  micfor%PSISM =PSISM
!  micfor%OLSGL =OLSGL
!  micfor%ORGC  =ORGC
!  micfor%RNO2Y =ystates0l(fid_RNO2Y)
!  micfor%RN2OY =ystates0l(fid_RN2OY)
!  micfor%RN2BY =ystates0l(fid_RN2BY)
!  micfor%ROXYY =ystates0l(fid_ROXYY)
!  micfor%ROXYF =ystates0l(fid_ROXYF)
!  micfor%RNHBY =ystates0l(fid_RNHBY)
!  micfor%RN3BY =ystates0l(fid_RN3BY)
!  micfor%RPOBY =ystates0l(fid_RPOBY)
!  micfor%RP1BY =ystates0l(fid_RP1BY)
!  micfor%ROQCY(0:jcplx1)=ROQCY(0:jcplx1)
!  micfor%ROQAY(0:jcplx1)=ROQAY(0:jcplx1)
!  micfor%RCH4L =RCH4L
!  micfor%ROXYL = ROXYL
!  micfor%CFOMC =CFOMC(1:ndbiomcp)
!  micfor%litrm=(L==0)
!  micfor%Lsurf=.True.
!  if(micfor%litrm)then
!  !  micstt%ZNH4TU=AMAX1(0.0,ZNH4S)+AMAX1(0.0,ZNH4B)
!  !  micstt%ZNO3TU=AMAX1(0.0,ZNO3S)+AMAX1(0.0,ZNO3B)
!  !  micstt%H1P4TU=AMAX1(0.0,H1PO4)+AMAX1(0.0,H1POB)
!  !  micstt%H2P4TU=AMAX1(0.0,H2PO4)+AMAX1(0.0,H2POB)
!  !  micstt%CNH4BU=CNH4B
!  !  micstt%CNH4SU=CNH4S
!  !  micstt%CH2P4U=CH2P4
!  !  micstt%CH2P4BU=CH2P4B
!  !  micstt%CH1P4U=CH1P4
!  !  micstt%CH1P4BU=CH1P4B
!  !  micstt%CNO3SU=CNO3S
!  !  micstt%CNO3BU=CNO3B
!  !  micstt%OSC13U=OSC(1,3)
!  !  micstt%OSN13U=OSN(1,3)
!  !  micstt%OSP13U=OSP(1,3)
!  !  micstt%OSC14U=OSC(1,4)
!  !  micstt%OSN14U=OSN(1,4)
!  !  micstt%OSP14U=OSP(1,4)
!  !  micstt%OSC24U=OSC(2,4)
!  !  micstt%OSN24U=OSN(2,4)
!  !  micstt%OSP24U=OSP(2,4)
  !  micfor%RNH4YU =RNH4Y
  !  micfor%RNO3YU =RNO3Y
  !  micfor%RPO4YU =RPO4Y
  !  micfor%RP14YU =RP14Y
  !  micfor%VOLWU  =VOLW
  !  micfor%CFOMCU=CFOMC(1:ndbiomcp)
!  else
  !  micfor%AEC=AEC
!  !  micstt%OXYG=OXYG
!  endif
!  micstt%CNH4B =CNH4B
!  micstt%CNH4S =CNH4S
!  micstt%CNO3S =CNO3S
!  micstt%CNO3B =CNO3B
!  micstt%CH2P4 =CH2P4
!  micstt%CH2P4B=CH2P4B
!  micstt%CH1P4 =CH1P4
!  micstt%CH1P4B=CH1P4B
!  micfor%RNH4Y =RNH4Y
!  micfor%RNO3Y =RNO3Y
!  micfor%RPO4Y =RPO4Y
!  micfor%RP14Y =RP14Y
!  micfor%VOLW  =VOLW

!  if(micfor%Lsurf)then
  !  micfor%BKVL0=BKVL
!  endif
!  micfor%DFGS(1:NPH)=DFGS(1:NPH)
!  micfor%FILM(1:NPH)=FILM(1:NPH)
!  micfor%THETPM(1:NPH)=THETPM(1:NPH)
!  micfor%VOLWM(1:NPH)=VOLWM(1:NPH)
!  micfor%TORT(1:NPH)=TORT(1:NPH)
!  micfor%VOLPM(1:NPH)=VOLPM(1:NPH)

!  micstt%EPOC=EPOC
!  micstt%EHUM=EHUM
!  micstt%ZNH4B=ZNH4B
!  micstt%ZNH4S=ZNH4S
!  micstt%ZNO3B=ZNO3B
!  micstt%ZNO3S=ZNO3S
!  micstt%H1POB=H1POB
!  micstt%H1PO4=H1PO4
!  micstt%ZNO2B=ZNO2B
!  micstt%ZNO2S=ZNO2S
!  micstt%H2POB=H2POB
!  micstt%H2PO4=H2PO4
!  micstt%CCO2S=CCO2S
!  micstt%CNO2S=CNO2S
!  micstt%CNO2B=CNO2B
!  micstt%CZ2OS=CZ2OS
!  micstt%Z2OS=Z2OS
!  micstt%COXYS=COXYS
!  micstt%OXYS=OXYS
!  micstt%SOXYL=SOXYL
!  micstt%COXYG=COXYG
!  micstt%CZ2GS=CZ2GS
!  micstt%CH2GS=CH2GS
!  micstt%H2GS=H2GS
!  micstt%CCH4G=CCH4G
!  micstt%CH4S=CH4S
!  micstt%SCH4L=SCH4L
!  micstt%ZNFN0=ZNFN0
!  micstt%ZNFNI=ZNFNI

!  micstt%FOSRH(0:jcplx1)=FOSRH(0:jcplx1)
!  micstt%OQC(0:jcplx1)=OQC(0:jcplx1)
!  micstt%OQN(0:jcplx1)=OQN(0:jcplx1)
!  micstt%OQP(0:jcplx1)=OQP(0:jcplx1)
!  micstt%OQA(0:jcplx1)=OQA(0:jcplx1)
!  micstt%OHC(0:jcplx1)=OHC(0:jcplx1)
!  micstt%OHN(0:jcplx1)=OHN(0:jcplx1)
!  micstt%OHP(0:jcplx1)=OHP(0:jcplx1)
!  micstt%OHA(0:jcplx1)=OHA(0:jcplx1)

!  micstt%OSC(1:jsken,0:jcplx1)=OSC(1:jsken,0:jcplx1)
!  micstt%OSA(1:jsken,0:jcplx1)=OSA(1:jsken,0:jcplx1)
!  micstt%OSN(1:jsken,0:jcplx1)=OSN(1:jsken,0:jcplx1)
!  micstt%OSP(1:jsken,0:jcplx1)=OSP(1:jsken,0:jcplx1)
!  micstt%ORC(1:ndbiomcp,0:jcplx1)=ORC(1:ndbiomcp,0:jcplx1)
!  micstt%ORN(1:ndbiomcp,0:jcplx1)=ORN(1:ndbiomcp,0:jcplx1)
!  micstt%ORP(1:ndbiomcp,0:jcplx1)=ORP(1:ndbiomcp,0:jcplx1)
!  micstt%CNOSC(1:jsken,0:jcplx1)=CNOSC(1:jsken,0:jcplx1)
!  micstt%CPOSC(1:jsken,0:jcplx1)=CPOSC(1:jsken,0:jcplx1)
!  micstt%OMC(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)=OMC(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)
!  micstt%OMN(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)=OMN(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)
!  micstt%OMP(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)=OMP(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)
!  micstt%OMCff(1:nlbiomcp,1:JG,1:NFGs)=OMCff(1:nlbiomcp,1:JG,1:NFGs)
!  micstt%OMNff(1:nlbiomcp,1:JG,1:NFGs)=OMNff(1:nlbiomcp,1:JG,1:NFGs)
!  micstt%OMPff(1:nlbiomcp,1:JG,1:NFGs)=OMPff(1:nlbiomcp,1:JG,1:NFGs)

!  micflx%RVMXC=RVMXC
!  micflx%RVMBC=RVMBC
!  micflx%RINHOff(1:JG,1:NFGs)=RINHOff(1:JG,1:NFGs)
!  micflx%RINHBff(1:JG,1:NFGs)=RINHBff(1:JG,1:NFGs)
!  micflx%RINOOff(1:JG,1:NFGs)=RINOOff(1:JG,1:NFGs)
!  micflx%RINOBff(1:JG,1:NFGs)=RINOBff(1:JG,1:NFGs)
!  micflx%RIPOOff(1:JG,1:NFGs)=RIPOOff(1:JG,1:NFGs)
!  micflx%RIPBOff(1:JG,1:NFGs)=RIPBOff(1:JG,1:NFGs)
!  micflx%RIPO1ff(1:JG,1:NFGs)=RIPO1ff(1:JG,1:NFGs)
!  micflx%RIPB1ff(1:JG,1:NFGs)=RIPB1ff(1:JG,1:NFGs)
!  micflx%ROXYSff(1:JG,1:NFGs)=ROXYSff(1:JG,1:NFGs)

!  micflx%RINHO(1:JG,1:NFGs,0:JCPLX1)=RINHO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RINHB(1:JG,1:NFGs,0:JCPLX1)=RINHB(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RINOO(1:JG,1:NFGs,0:JCPLX1)=RINOO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RINOB(1:JG,1:NFGs,0:JCPLX1)=RINOB(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPOO(1:JG,1:NFGs,0:JCPLX1)=RIPOO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPBO(1:JG,1:NFGs,0:JCPLX1)=RIPBO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPO1(1:JG,1:NFGs,0:JCPLX1)=RIPO1(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPB1(1:JG,1:NFGs,0:JCPLX1)=RIPB1(1:JG,1:NFGs,0:JCPLX1)
!  micflx%ROXYS(1:JG,1:NFGs,0:JCPLX1)=ROXYS(1:JG,1:NFGs,0:JCPLX1)
  end associate
  end subroutine BatchModelConfig

! ----------------------------------------------------------------------
!  subroutine ReadForc(forc)
!
! DESCRIPTION
! read forcing data

!  use MicForcTypeMod, only : micforctype
!  implicit none
!  type(forc_type), intent(in) :: forc


!  end subroutine ReadForc
! ----------------------------------------------------------------------


  subroutine Initboxbgc(nvars)
!
! DESCRIPTION
! Initialize the id of relevant model variables &
! obtain total number of variables
  implicit none

  integer, intent(out) :: nvars
  integer :: itemp

  associate(                        &
    jcplx    => micpar%jcplx      , &
    jcplx1   => micpar%jcplx1     , &
    jsken    => micpar%jsken      , &
    NFGs     => micpar%NFGs       , &
    JG       => micpar%jguilds    , &
    ndbiomcp => micpar%ndbiomcp   , &
    nlbiomcp => micpar%nlbiomcp     &
  )
  itemp=0
  cid_ZMG    =addone(itemp)
  cid_ZNA    =addone(itemp)
  cid_ZKA    =addone(itemp)
  cid_CO2S   =addone(itemp)
  cid_CH1P1  =addone(itemp)
  cid_CH1PB  =addone(itemp)
  cid_CH2P1  =addone(itemp)
  cid_CH2PB  =addone(itemp)
  cid_CN31   =addone(itemp)
  cid_CN3B   =addone(itemp)
  cid_CN41   =addone(itemp)
  cid_CN4B   =addone(itemp)
  cid_XN41   =addone(itemp)
  cid_XN4B   =addone(itemp)
  cid_X1P1B  =addone(itemp)
  cid_X2P1B  =addone(itemp)
  cid_XH11B  =addone(itemp)
  cid_XH1P1  =addone(itemp)
  cid_XH21B  =addone(itemp)
  cid_XH2P1  =addone(itemp)
  cid_XOH11  =addone(itemp)
  cid_XOH21  =addone(itemp)
  cid_PALPO1 =addone(itemp)
  cid_PALPOB =addone(itemp)
  cid_PCAPD1 =addone(itemp)
  cid_PCAPDB =addone(itemp)
  cid_PCAPH1 =addone(itemp)
  cid_PCAPHB =addone(itemp)
  cid_PCAPM1 =addone(itemp)
  cid_PCAPMB =addone(itemp)
  cid_PFEPO1 =addone(itemp)
  cid_PFEPOB =addone(itemp)

  fid_TRN4S = addone(itemp)
  fid_TRN4B = addone(itemp)
  fid_TRN3S = addone(itemp)
  fid_TRN3B = addone(itemp)
  fid_TRH1P = addone(itemp)
  fid_TRH2P = addone(itemp)
  fid_TRH1B = addone(itemp)
  fid_TRH2B = addone(itemp)
  fid_TRXN4 = addone(itemp)
  fid_TRXNB = addone(itemp)
  fid_TRXH1 = addone(itemp)
  fid_TRXH2 = addone(itemp)
  fid_TRX1P = addone(itemp)
  fid_TRX2P = addone(itemp)
  fid_TRBH1 = addone(itemp)
  fid_TRBH2 = addone(itemp)
  fid_TRB1P = addone(itemp)
  fid_TRB2P = addone(itemp)
  fid_TRALPO= addone(itemp)
  fid_TRFEPO= addone(itemp)
  fid_TRCAPD= addone(itemp)
  fid_TRCAPH= addone(itemp)
  fid_TRCAPM= addone(itemp)
  fid_TRALPB= addone(itemp)
  fid_TRFEPB= addone(itemp)
  fid_TRCPDB= addone(itemp)
  fid_TRCPHB= addone(itemp)
  fid_TRCPMB= addone(itemp)
  fid_TRAL  = addone(itemp)

  cid_oqc_b=addone(itemp);cid_oqc_e=cid_oqc_b+jcplx1;itemp=cid_oqc_e
  cid_oqn_b=addone(itemp);cid_oqn_e=cid_oqn_b+jcplx1;itemp=cid_oqn_e
  cid_oqp_b=addone(itemp);cid_oqp_e=cid_oqp_b+jcplx1;itemp=cid_oqp_e
  cid_oqa_b=addone(itemp);cid_oqa_e=cid_oqa_b+jcplx1;itemp=cid_oqa_e
  cid_ohc_b=addone(itemp);cid_ohc_e=cid_ohc_b+jcplx1;itemp=cid_ohc_e
  cid_ohn_b=addone(itemp);cid_ohn_e=cid_ohn_b+jcplx1;itemp=cid_ohn_e
  cid_ohp_b=addone(itemp);cid_ohp_e=cid_ohp_b+jcplx1;itemp=cid_ohp_e
  cid_oha_b=addone(itemp);cid_oha_e=cid_oha_b+jcplx1;itemp=cid_oha_e
  cid_osc_b=addone(itemp);cid_osc_e=cid_osc_b+jsken*jcplx-1;itemp=cid_osc_e
  cid_osa_b=addone(itemp);cid_osa_e=cid_osa_b+jsken*jcplx-1;itemp=cid_osa_e
  cid_osn_b=addone(itemp);cid_osn_e=cid_osn_b+jsken*jcplx-1;itemp=cid_osn_e
  cid_osp_b=addone(itemp);cid_osp_e=cid_osp_b+jsken*jcplx-1;itemp=cid_osp_e
  cid_orc_b=addone(itemp);cid_orc_e=cid_orc_b+ndbiomcp*jcplx-1;itemp=cid_orc_e
  cid_orn_b=addone(itemp);cid_orn_e=cid_orn_b+ndbiomcp*jcplx-1;itemp=cid_orn_e
  cid_orp_b=addone(itemp);cid_orp_e=cid_orp_b+ndbiomcp*jcplx-1;itemp=cid_orp_e
  cid_omc_b=addone(itemp);cid_omc_e=cid_omc_b+nlbiomcp*JG*NFGs*jcplx-1;itemp=cid_omc_e
  cid_omn_b=addone(itemp);cid_omn_e=cid_omn_b+nlbiomcp*JG*NFGs*jcplx-1;itemp=cid_omn_e
  cid_omp_b=addone(itemp);cid_omp_e=cid_omp_b+nlbiomcp*JG*NFGs*jcplx-1;itemp=cid_omp_e
  cid_omcff_b=addone(itemp);cid_omcff_e=cid_omcff_b+nlbiomcp*JG*NFGs-1;itemp=cid_omcff_e
  cid_omnff_b=addone(itemp);cid_omnff_e=cid_omnff_b+nlbiomcp*JG*NFGs-1;itemp=cid_omnff_e
  cid_ompff_b=addone(itemp);cid_ompff_e=cid_ompff_b+nlbiomcp*JG*NFGs-1;itemp=cid_ompff_e

  fid_ROXYY=addone(itemp)
  fid_RNH4Y=addone(itemp)
  fid_RNO3Y=addone(itemp)
  fid_RNO2Y=addone(itemp)
  fid_RN2OY=addone(itemp)
  fid_RPO4Y=addone(itemp)
  fid_RP14Y=addone(itemp)
  fid_RNHBY=addone(itemp)
  fid_RN3BY=addone(itemp)
  fid_RN2BY=addone(itemp)
  fid_RPOBY=addone(itemp)
  fid_RP1BY=addone(itemp)
  fid_ROQCY_b=addone(itemp);fid_ROQCY_e=fid_ROQCY_b+jcplx1;itemp=fid_ROQCY_e
  fid_ROQAY_b=addone(itemp);fid_ROQAY_e=fid_ROQAY_b+jcplx1;itemp=fid_ROQAY_e

  cid_ZNH4B=addone(itemp)
  cid_ZNH4S=addone(itemp)
  cid_ZNO3B=addone(itemp)
  cid_ZNO3S=addone(itemp)
  cid_H1POB=addone(itemp)
  cid_H1PO4=addone(itemp)
  cid_ZNO2B=addone(itemp)
  cid_ZNO2S=addone(itemp)
  cid_H2POB=addone(itemp)
  cid_H2PO4=addone(itemp)
  cid_CCO2S=addone(itemp)
  cid_CNO2S=addone(itemp)
  cid_CNO2B=addone(itemp)
  cid_CZ2OS=addone(itemp)
  cid_Z2OS =addone(itemp)
  cid_COXYS=addone(itemp)
  cid_OXYS =addone(itemp)
  cid_COXYG=addone(itemp)
  cid_CZ2GS=addone(itemp)
  cid_CH2GS=addone(itemp)
  cid_H2GS =addone(itemp)
  cid_CCH4G=addone(itemp)
  cid_CH4S =addone(itemp)
  cid_ZNFN0=addone(itemp)
  cid_ZNFNI=addone(itemp)

  nvars=itemp
  end associate
  end subroutine Initboxbgc
! ----------------------------------------------------------------------
  subroutine UpdateStateVars(micfor, micstt,micflx,nvars,ystatesfl)
!
! DESCRIPTION
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micfluxtype), intent(in) :: micflx
  type(micsttype)  , intent(in) :: micstt
  integer , intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)

  real(r8) :: DC,DN,DP,OC,ON,OP
  real(r8) :: ORGC,ORGN,ORGR
  integer :: K,N,NGL,M
  associate(                         &
    jcplx1    => micpar%jcplx1     , &
    jcplx     => micpar%jcplx      , &
    JG        => micpar%jguilds    , &
    NFGs      => micpar%NFGs       , &
    jsken     => micpar%jsken      , &
    nlbiomcp  => micpar%nlbiomcp   , &
    ndbiomcp  => micpar%ndbiomcp   , &
    is_litter => micpar%is_litter  , &
    VOLW      => micfor%VOLW         &
  )
!atmospheric gaseous CO2,CH4,O2,NH3,N2,N2O,H2

  ystatesfl(cid_ZNH4B)=micstt%ZNH4B
  ystatesfl(cid_ZNH4S)=micstt%ZNH4S
  ystatesfl(cid_ZNO3B)=micstt%ZNO3B
  ystatesfl(cid_ZNO3S)=micstt%ZNO3S
  ystatesfl(cid_H1POB)=micstt%H1POB
  ystatesfl(cid_H1PO4)=micstt%H1PO4
  ystatesfl(cid_ZNO2B)=micstt%ZNO2B
  ystatesfl(cid_ZNO2S)=micstt%ZNO2S
  ystatesfl(cid_H2POB)=micstt%H2POB
  ystatesfl(cid_H2PO4)=micstt%H2PO4
  ystatesfl(cid_CCO2S)=micstt%CCO2S !-RCO2O/VOLW
  ystatesfl(cid_CNO2S)=micstt%ZNO2S/(VOLW*micfor%VLNO3)
  ystatesfl(cid_CNO2B)=micstt%ZNO2B/(VOLW*micfor%VLNOB)
  ystatesfl(cid_CZ2OS)=micstt%Z2OS/VOLW
  ystatesfl(cid_Z2OS) =micstt%Z2OS
  ystatesfl(cid_COXYS)=micstt%OXYS/VOLW
  ystatesfl(cid_OXYS) =micstt%OXYS  !-RUPOXO
  ystatesfl(cid_COXYG)=micstt%COXYG
  ystatesfl(cid_CZ2GS)=micstt%CZ2GS
  ystatesfl(cid_CH2GS)=micstt%CH2GS
  ystatesfl(cid_H2GS) =micstt%H2GS
  ystatesfl(cid_CCH4G)=micstt%CCH4G
  ystatesfl(cid_CH4S) =micstt%CH4S  !-RCH4O
  ystatesfl(cid_ZNFN0)=micstt%ZNFN0
  ystatesfl(cid_ZNFNI)=micstt%ZNFNI

! the following variables are updated in the microbial model
  ystatesfl(cid_oqc_b:cid_oqc_e)=micstt%OQC(0:jcplx1)
  ystatesfl(cid_oqn_b:cid_oqn_e)=micstt%OQN(0:jcplx1)
  ystatesfl(cid_oqp_b:cid_oqp_e)=micstt%OQP(0:jcplx1)
  ystatesfl(cid_oqa_b:cid_oqa_e)=micstt%OQA(0:jcplx1)
  ystatesfl(cid_ohc_b:cid_ohc_e)=micstt%OHC(0:jcplx1)
  ystatesfl(cid_ohn_b:cid_ohn_e)=micstt%OHN(0:jcplx1)
  ystatesfl(cid_ohp_b:cid_ohp_e)=micstt%OHP(0:jcplx1)
  ystatesfl(cid_oha_b:cid_oha_e)=micstt%OHA(0:jcplx1)
  ystatesfl(cid_osc_b:cid_osc_e)=reshape(micstt%OSC(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(cid_osa_b:cid_osa_e)=reshape(micstt%OSA(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(cid_osn_b:cid_osn_e)=reshape(micstt%OSN(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(cid_osp_b:cid_osp_e)=reshape(micstt%OSP(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(cid_orc_b:cid_orc_e)=reshape(micstt%ORC(1:ndbiomcp,0:jcplx1),&
    (/ndbiomcp*jcplx/))
  ystatesfl(cid_orn_b:cid_orn_e)=reshape(micstt%ORN(1:ndbiomcp,0:jcplx1),&
    (/ndbiomcp*jcplx/))
  ystatesfl(cid_orp_b:cid_orp_e)=reshape(micstt%ORP(1:ndbiomcp,0:jcplx1),&
    (/ndbiomcp*jcplx/))
  ystatesfl(cid_omc_b:cid_omc_e)=reshape(micstt%OMC(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1),&
    (/nlbiomcp*JG*NFGs*jcplx/))
  ystatesfl(cid_omn_b:cid_omn_e)=reshape(micstt%OMN(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1),&
    (/nlbiomcp*JG*NFGs*jcplx/))
  ystatesfl(cid_omp_b:cid_omp_e)=reshape(micstt%OMP(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1),&
    (/nlbiomcp*JG*NFGs*jcplx/))
  ystatesfl(cid_omcff_b:cid_omcff_e)=reshape(micstt%OMCff(1:nlbiomcp,1:JG,1:NFGs),&
    (/nlbiomcp*JG*NFGs/))
  ystatesfl(cid_omnff_b:cid_omnff_e)=reshape(micstt%OMNff(1:nlbiomcp,1:JG,1:NFGs),&
    (/nlbiomcp*JG*NFGs/))
  ystatesfl(cid_ompff_b:cid_ompff_e)=reshape(micstt%OMPff(1:nlbiomcp,1:JG,1:NFGs),&
    (/nlbiomcp*JG*NFGs/))

! summarize diagnostic fluxes
  DO K=0,jcplx1
    IF(.not.micfor%litrm.or.(micpar%is_litter(K)))THEN
      DO N=1,NFGs
        DO NGL=1,JG
          ystatesfl(fid_ROXYY)=ystatesfl(fid_ROXYY)!+ROXYS(NGL,N,K)
          ystatesfl(fid_RNH4Y)=ystatesfl(fid_RNH4Y)!+RVMX4(NGL,N,K)+RINHO(NGL,N,K)
          ystatesfl(fid_RNO3Y)=ystatesfl(fid_RNO3Y)!+RVMX3(NGL,N,K)+RINOO(NGL,N,K)
          ystatesfl(fid_RNO2Y)=ystatesfl(fid_RNO2Y)!+RVMX2(NGL,N,K)
          ystatesfl(fid_RN2OY)=ystatesfl(fid_RN2OY)!+RVMX1(NGL,N,K)
          ystatesfl(fid_RPO4Y)=ystatesfl(fid_RPO4Y)!+RIPOO(NGL,N,K)
          ystatesfl(fid_RP14Y)=ystatesfl(fid_RP14Y)!+RIPO1(NGL,N,K)
          ystatesfl(fid_RNHBY)=ystatesfl(fid_RNHBY)!+RVMB4(NGL,N,K)+RINHB(NGL,N,K)
          ystatesfl(fid_RN3BY)=ystatesfl(fid_RN3BY)!+RVMB3(NGL,N,K)+RINOB(NGL,N,K)
          ystatesfl(fid_RN2BY)=ystatesfl(fid_RN2BY)!+RVMB2(NGL,N,K)
          ystatesfl(fid_RPOBY)=ystatesfl(fid_RPOBY)!+RIPBO(NGL,N,K)
          ystatesfl(fid_RP1BY)=ystatesfl(fid_RP1BY)!+RIPB1(NGL,N,K)
          ystatesfl(fid_ROQCY_b+K)=ystatesfl(fid_ROQCY_b+K)!+ROQCS(NGL,N,K)
          ystatesfl(fid_ROQAY_b+K)=ystatesfl(fid_ROQAY_b+K)!+ROQAS(NGL,N,K)
        enddo
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NFGs
    DO NGL=1,JG
      ystatesfl(fid_ROXYY)=ystatesfl(fid_ROXYY) !+ROXYSff(NGL,N)
      ystatesfl(fid_RNH4Y)=ystatesfl(fid_RNH4Y) !+RVMX4ff(NGL,N)+RINHOff(NGL,N)
      ystatesfl(fid_RNO3Y)=ystatesfl(fid_RNO3Y) !+RVMX3ff(NGL,N)+RINOOff(NGL,N)
      ystatesfl(fid_RNO2Y)=ystatesfl(fid_RNO2Y) !+RVMX2ff(NGL,N)
      ystatesfl(fid_RN2OY)=ystatesfl(fid_RN2OY) !+RVMX1ff(NGL,N)
      ystatesfl(fid_RPO4Y)=ystatesfl(fid_RPO4Y) !+RIPOOff(NGL,N)
      ystatesfl(fid_RP14Y)=ystatesfl(fid_RP14Y) !+RIPO1ff(NGL,N)
      ystatesfl(fid_RNHBY)=ystatesfl(fid_RNHBY) !+RVMB4ff(NGL,N)+RINHBff(NGL,N)
      ystatesfl(fid_RN3BY)=ystatesfl(fid_RN3BY) !+RVMB3ff(NGL,N)+RINOBff(NGL,N)
      ystatesfl(fid_RN2BY)=ystatesfl(fid_RN2BY) !+RVMB2ff(NGL,N)
      ystatesfl(fid_RPOBY)=ystatesfl(fid_RPOBY) !+RIPBOff(NGL,N)
      ystatesfl(fid_RP1BY)=ystatesfl(fid_RP1BY) !+RIPB1ff(NGL,N)
    enddo
  ENDDO

  ystatesfl(fid_RNO2Y)=ystatesfl(fid_RNO2Y) !+RVMXC
  ystatesfl(fid_RN2BY)=ystatesfl(fid_RN2BY) !+RVMBC

  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8

  DO K=0,jcplx1
    IF(is_litter(K))THEN
      DO N=1,NFGs
        DO  M=1,nlbiomcp
          DO NGL=1,JG
!            DC=DC+OMC(M,NGL,N,K)
!            DN=DN+OMN(M,NGL,N,K)
!            DP=DP+OMP(M,NGL,N,K)
          ENDDO
        enddo
      ENDDO
    ELSE
      DO N=1,NFGs
        DO  M=1,nlbiomcp
          DO NGL=1,JG
!            OC=OC+OMC(M,NGL,N,K)
!            ON=ON+OMN(M,NGL,N,K)
!            OP=OP+OMP(M,NGL,N,K)
          enddo
        enddo
      ENDDO
    ENDIF
  ENDDO
! abstract complex
  DO  N=1,NFGs
    DO  M=1,nlbiomcp
      DO NGL=1,JG
!        OC=OC+OMCff(M,NGL,N)
!        ON=ON+OMNff(M,NGL,N)
!        OP=OP+OMPff(M,NGL,N)
      enddo
    enddo
  ENDDO
! microbial residue
  DO K=0,jcplx1
    IF(is_litter(K))THEN
      DO M=1,ndbiomcp
!        DC=DC+ORC(M,K)
!        DN=DN+ORN(M,K)
!        DP=DP+ORP(M,K)
      ENDDO
!      DC=DC+OQC(K)+OQCH(K)+OHC(K)+OQA(K)+OQAH(K)+OHA(K)
!      DN=DN+OQN(K)+OQNH(K)+OHN(K)
!      DP=DP+OQP(K)+OQPH(K)+OHP(K)
      DO M=1,jsken
!        DC=DC+OSC(M,K)
!        DN=DN+OSN(M,K)
!        DP=DP+OSP(M,K)
      ENDDO
    ELSE
      DO M=1,ndbiomcp
!        OC=OC+ORC(M,K)
!        ON=ON+ORN(M,K)
!        OP=OP+ORP(M,K)
      ENDDO
!      OC=OC+OQC(K)+OQCH(K)+OHC(K)+OQA(K)+OQAH(K)+OHA(K)
!      ON=ON+OQN(K)+OQNH(K)+OHN(K)
!      OP=OP+OQP(K)+OQPH(K)+OHP(K)
      DO M=1,jsken
!        OC=OC+OSC(M,K)
!        ON=ON+OSN(M,K)
!        OP=OP+OSP(M,K)
      ENDDO
    ENDIF
  ENDDO
! DC is for litter complex, and OC is for POM and humus complex
  ORGC=DC+OC
  ORGN=DN+ON
  ORGR=DC

  end associate
  end subroutine UpdateStateVars

! ----------------------------------------------------------------------

  subroutine getvarlist(nvars, varl, varlnml, unitl, vartypes)

  use histMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(nvars)            !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(nvars)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(nvars)           !variable unit
  integer                         ,intent(out) :: vartypes(nvars)        !variable type, flx or state

  integer :: iknen,icplx
  integer :: jj,ll,k,m,n,ngl

  associate(                        &
    jcplx1    => micpar%jcplx1    , &
    JG        => micpar%jguilds   , &
    jsken     => micpar%jsken     , &
    NFGs      => micpar%NFGs      , &
    nlbiomcp  => micpar%nlbiomcp  , &
    ndbiomcp  => micpar%ndbiomcp    &
  )

  call getchemvarlist(nvars, varl, varlnml, unitl, vartypes)

  varl(cid_ZMG) ='ZMG';varlnml(cid_ZMG)='soil aqueous Mg content micropore'
  unitl(cid_ZMG)='mol d-2';vartypes(cid_ZMG)=var_state_type

  varl(cid_ZNA) ='ZNA';varlnml(cid_ZNA)='soil aqueous Na content micropore'
  unitl(cid_ZNA) ='mol d-2';vartypes(cid_ZNA)=var_state_type

  varl(cid_ZKA) ='ZKA';varlnml(cid_ZKA)='soil aqueous K content micropore'
  unitl(cid_ZKA) ='mol d-2';vartypes(cid_ZKA)=var_state_type

  varl(cid_CO2S)='CO2S';varlnml(cid_CO2S)='aqueous CO2 concentration micropore';
  unitl(cid_CO2S)='gC d-2';vartypes(cid_CO2S)=var_state_type

  varl(cid_ZNH4B)='ZNH4B';varlnml(cid_ZNH4B)='band soil micropore NH4(+) mass'
  unitl(cid_ZNH4B)='gN d-2';vartypes(cid_ZNH4B)=var_state_type

  varl(cid_ZNH4S)='ZNH4S';varlnml(cid_ZNH4S)='non-band soil micropore NH4(+) mass'
  unitl(cid_ZNH4S)='gN d-2';vartypes(cid_ZNH4S)=var_state_type

  varl(cid_ZNO3B)='ZNO3B';varlnml(cid_ZNO3B)='band soil micropore NO3(-) mass'
  unitl(cid_ZNO3B)='gN d-2';vartypes(cid_ZNO3B)=var_state_type

  varl(cid_ZNO3S)='ZNO3S';varlnml(cid_ZNO3S)='non-band soil micropore NO3(-) mass'
  unitl(cid_ZNO3S)='gN d-2';vartypes(cid_ZNO3S)=var_state_type

  varl(cid_H1POB)='H1POB';varlnml(cid_H1POB)='band soil micropore aqueous HPO4 content'
  unitl(cid_H1POB)='gP m-2';vartypes(cid_H1POB)=var_state_type

  varl(cid_H1PO4)='H1PO4';varlnml(cid_H1PO4)='non-band soil micropore aqueous HPO4(--) content';
  unitl(cid_H1PO4)='gP m-2';vartypes(cid_H1PO4)=var_state_type

  varl(cid_ZNO2B)='ZNO2B';varlnml(cid_ZNO2B)='band soil micropore NO2(-) mass'
  unitl(cid_ZNO2B)='gN m-2';vartypes(cid_ZNO2B)=var_state_type

  varl(cid_ZNO2S)='ZNO2S';varlnml(cid_ZNO2S)='non-band soil micropore NO2(-) mass'
  unitl(cid_ZNO2S)='gN m-2';vartypes(cid_ZNO2S)=var_state_type

  varl(cid_H2POB)='H2POB';varlnml(cid_H2POB)='band soil micropore H2PO4 mass'
  unitl(cid_H2POB)='gP m-2';vartypes(cid_H2POB)=var_state_type

  varl(cid_H2PO4)='H2PO4';varlnml(cid_H2PO4)='non-band soil micropore H2PO4 mass'
  unitl(cid_H2PO4)='gP m-2';vartypes(cid_H2PO4)=var_state_type

  varl(cid_CCO2S)='CCO2S';varlnml(cid_CCO2S)='soil micropore aqueous CO2 concentration'
  unitl(cid_CCO2S)='gC m-3';vartypes(cid_CCO2S)=var_state_type

  varl(cid_CNO2S)='CNO2S';varlnml(cid_CNO2S)='non-band soil micropore NO2 concentration'
  unitl(cid_CNO2S)='gN m-3';vartypes(cid_CNO2S)=var_state_type

  varl(cid_CNO2B)='CNO2B';varlnml(cid_CNO2B)='band soil micropore NO2 concentration'
  unitl(cid_CNO2B)='gN m-3';vartypes(cid_CNO2B)=var_state_type

  varl(cid_CZ2OS)='CZ2OS';varlnml(cid_CZ2OS)='soil micropore aqueous N2O concentration'
  unitl(cid_CZ2OS)='gN m-3';vartypes(cid_CZ2OS)=var_state_type

  varl(cid_Z2OS) ='Z2OS';varlnml(cid_Z2OS)='soil micropore aqueous N2O mass'
  unitl(cid_Z2OS)='gN d-2';vartypes(cid_Z2OS)=var_state_type

  varl(cid_COXYS)='COXYS';varlnml(cid_COXYS)='soil micropore aqueous O2 concentration'
  unitl(cid_COXYS)='g m-3';vartypes(cid_COXYS)=var_state_type

  varl(cid_OXYS) ='OXYS';varlnml(cid_OXYS)='soil micropore aqueous O2 mass'
  unitl(cid_OXYS)='g d-2';vartypes(cid_OXYS)=var_state_type

  varl(cid_COXYG)='COXYG';varlnml(cid_COXYG)='soil micropore gaseous O2 concentration'
  unitl(cid_COXYG)='g m-3';vartypes(cid_COXYG)=var_state_type

  varl(cid_CZ2GS)='CZ2GS';varlnml(cid_CZ2GS)='soil micropore aqueous N2 concentration'
  unitl(cid_CZ2GS)='gN m-3';vartypes(cid_CZ2GS)=var_state_type

  varl(cid_CH2GS)='CH2GS';varlnml(cid_CH2GS)='soil micropore aqueous H2 concentration'
  unitl(cid_CH2GS)='g m-3';vartypes(cid_CH2GS)=var_state_type

  varl(cid_H2GS) ='H2GS';varlnml(cid_H2GS)='soil micropore aqueous H2 mass'
  unitl(cid_H2GS)='g d-2';vartypes(cid_H2GS)=var_state_type

  varl(cid_CCH4G)='CCH4G';varlnml(cid_CCH4G)='soil micropore gaseous CH4 concentration'
  unitl(cid_CCH4G)='gC m-3';vartypes(cid_CCH4G)=var_state_type

  varl(cid_CH4S) ='CH4S';varlnml(cid_CH4S)='soil micropore aqueous CH4 mass'
  unitl(cid_CH4S)='gC d-2';vartypes(cid_CH4S)=var_state_type

  varl(cid_ZNFN0)='ZNFN0';varlnml(cid_ZNFN0)='initial nitrification inhibition activity'
  unitl(cid_ZNFN0)='none';vartypes(cid_ZNFN0)=var_state_type

  varl(cid_ZNFNI)='ZNFNI';varlnml(cid_ZNFNI)='current nitrification inhibition activity'
  unitl(cid_ZNFNI)='none';vartypes(cid_ZNFNI)=var_state_type

  do jj=cid_oqc_b,cid_oqc_e
    write(varl(jj),'(A,I1)')'OQC',jj-cid_oqc_b
    varlnml(jj)='micropore dissolved organic C mass in complex '//trim(micpar%cplxname(jj-cid_oqc_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oqn_b,cid_oqn_e
    write(varl(jj),'(A,I1)')'OQN',jj-cid_oqn_b
    varlnml(jj)='micropore dissolved N mass in complex '//trim(micpar%cplxname(jj-cid_oqn_b))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oqp_b,cid_oqp_e
    write(varl(jj),'(A,I1)')'OQP',jj-cid_oqn_b
    varlnml(jj)='micropore dissolved N mass in complex '//trim(micpar%cplxname(jj-cid_oqp_b))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oqa_b,cid_oqa_e
    write(varl(jj),'(A,I1)')'OQA',jj-cid_oqa_b
    varlnml(jj)='micropore dissolved acetate mass in complex '//trim(micpar%cplxname(jj-cid_oqa_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_ohc_b,cid_ohc_e
    write(varl(jj),'(A,I1)')'OHC',jj-cid_ohc_b
    varlnml(jj)='adsorbed soil C mass in complex'//trim(micpar%cplxname(jj-cid_ohc_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_ohn_b,cid_ohn_e
    write(varl(jj),'(A,I1)')'OHN',jj-cid_ohn_b
    varlnml(jj)='adsorbed soil N mass in complex'//trim(micpar%cplxname(jj-cid_ohn_b))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_ohp_b,cid_ohp_e
    write(varl(jj),'(A,I1)')'OHP',jj-cid_ohp_b
    varlnml(jj)='adsorbed soil P mass in complex'//trim(micpar%cplxname(jj-cid_ohp_b))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oha_b,cid_oha_e
    write(varl(jj),'(A,I1)')'OHA',jj-cid_ohp_b
    varlnml(jj)='adsorbed soil acetate mass in complex'//trim(micpar%cplxname(jj-cid_oha_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_osc_b,cid_osc_e
    iknen=jj-cid_osc_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSC',iknen+1,icplx
    varlnml(jj)='humus soil C as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_osn_b,cid_osn_e
    iknen=jj-cid_osn_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSN',iknen+1,icplx
    varlnml(jj)='humus soil N as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_osp_b,cid_osp_e
    iknen=jj-cid_osp_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSP',iknen+1,icplx
    varlnml(jj)='humus soil P as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_osa_b,cid_osa_e
    iknen=jj-cid_osc_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSA',iknen+1,icplx
    varlnml(jj)='colonized humus soil C as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_orc_b,cid_orc_e
    iknen=jj-cid_orc_b
    icplx=floor((iknen-1.e-3_r8)/ndbiomcp)
    iknen=mod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORC',iknen+1,icplx
    varlnml(jj)='microbial residue C as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_orn_b,cid_orn_e
    iknen=jj-cid_orc_b
    icplx=floor((iknen-1.e-3_r8)/ndbiomcp)
    iknen=mod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORN',iknen+1,icplx
    varlnml(jj)='microbial residue N as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_orp_b,cid_orp_e
    iknen=jj-cid_orc_b
    icplx=floor((iknen-1.e-3_r8)/ndbiomcp)
    iknen=mod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORP',iknen+1,icplx
    varlnml(jj)='microbial residue P as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo

  jj=0
  DO k=0,jcplx1
  DO N=1,NFGs
  DO NGL=1,JG
  DO M=1,nlbiomcp
    ll=cid_omc_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMC'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass C in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gC d-2'
    vartypes(ll)=var_state_type

    ll=cid_omn_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMN'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass N in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gN d-2'
    vartypes(ll)=var_state_type

    ll=cid_omp_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMP'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass P in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gP d-2'
    vartypes(ll)=var_state_type
    jj=jj+1
  enddo
  ENDDO
  ENDDO
  ENDDO

  jj=0
  DO N=1,NFGs
  DO NGL=1,JG
  DO M=1,nlbiomcp
    ll=cid_omcff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMC'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass C in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gC d-2'
    vartypes(ll)=var_state_type

    ll=cid_omnff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMN'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass N in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gN d-2'
    vartypes(ll)=var_state_type

    ll=cid_ompff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMP'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass P in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gP d-2'
    vartypes(ll)=var_state_type

    jj=jj+1
  enddo
  ENDDO
  ENDDO

  varl(fid_ROXYY)='ROXYY';varlnml(fid_ROXYY)='total root + microbial O2 uptake potential'
  unitl(fid_ROXYY)='g d-2 h-1'; vartypes(fid_ROXYY)=var_flux_type

  varl(fid_RNH4Y)='RNH4Y';varlnml(fid_RNH4Y)='total root + microbial NH4 uptake potential non-band soil'
  unitl(fid_RNH4Y)='gN d-2 h-1'; vartypes(fid_RNH4Y)=var_flux_type

  varl(fid_RNO3Y)='RNO3Y';varlnml(fid_RNO3Y)='total root + microbial NO3 uptake potential non-band soil'
  unitl(fid_RNO3Y)='gN d-2 h-1'; vartypes(fid_RNO3Y)=var_flux_type

  varl(fid_RNO2Y)='RNO2Y';varlnml(fid_RNO2Y)='total root + microbial NO2 uptake potential non-band soil'
  unitl(fid_RNO2Y)='gN d-2 h-1'; vartypes(fid_RNO2Y)=var_flux_type

  varl(fid_RN2OY)='RN2OY';varlnml(fid_RN2OY)='total root + microbial N2O uptake potential';
  unitl(fid_RN2OY)='gN d-2 h-1'; vartypes(fid_RN2OY)=var_flux_type

  varl(fid_RPO4Y)='RPO4Y';varlnml(fid_RPO4Y)='total root + microbial PO4 uptake potential non-band soil'
  unitl(fid_RPO4Y)='gP d-2 h-1'; vartypes(fid_RPO4Y)=var_flux_type

  varl(fid_RP14Y)='RP14Y';varlnml(fid_RP14Y)='total root + microbial HPO4 uptake non-band soil'
  unitl(fid_RP14Y)='gP d-2 h-1'; vartypes(fid_RP14Y)=var_flux_type

  varl(fid_RNHBY)='RNHBY';varlnml(fid_RNHBY)='total root + microbial NH4 uptake potential band soil'
  unitl(fid_RNHBY)='gN d-2 h-1'; vartypes(fid_RNHBY)=var_flux_type

  varl(fid_RN3BY)='RN3BY';varlnml(fid_RN3BY)='total root + microbial NO3 uptake potential band soil'
  unitl(fid_RN3BY)='gN d-2 h-1'; vartypes(fid_RN3BY)=var_flux_type

  varl(fid_RN2BY)='RN2BY';varlnml(fid_RN2BY)='total root + microbial NO2 uptake potential band soil'
  unitl(fid_RN2BY)='gN d-2 h-1'; vartypes(fid_RN2BY)=var_flux_type

  varl(fid_RPOBY)='RPOBY';varlnml(fid_RPOBY)='total root + microbial PO4 uptake potential band soil'
  unitl(fid_RPOBY)='gP d-2 h-1'; vartypes(fid_RPOBY)=var_flux_type

  varl(fid_RP1BY)='RP1BY';varlnml(fid_RP1BY)='total root + microbial HPO4 uptake potential band soil';
  unitl(fid_RP1BY)='gP d-2 h-1'; vartypes(fid_RP1BY)=var_flux_type

  do jj =fid_ROQCY_b,fid_ROQCY_e
    write(varl(jj),'(A,I2.2)')'ROQCY',jj-fid_ROQCY_b
    varlnml(jj)='total root + microbial DOC uptake in complex ' &
      //micpar%cplxname(jj-fid_ROQCY_b)
    vartypes(jj)=var_flux_type
    unitl(jj)='gC d-2 h-1'
  enddo
  do jj =fid_ROQAY_b,fid_ROQAY_e
    write(varl(jj),'(A,I2.2)')'ROQAY',jj-fid_ROQAY_b
    varlnml(jj)='total root + microbial acetate uptake in complex ' &
      //micpar%cplxname(jj-fid_ROQCY_b)
    vartypes(jj)=var_flux_type
    unitl(jj)='gC d-2 h-1'
  enddo
  end associate
  end subroutine getvarlist
! ----------------------------------------------------------------------
  subroutine RunMicBGC(nvars, ystates0l, ystatesfl, micfor,micstt,micflx, err_status)
!
!
  use MicBGCMod, only : SoilBGCOneLayer
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(in) :: ystates0l(nvars)
  real(r8), intent(out) :: ystatesfl(nvars)
  type(micforctype), intent(in)    :: micfor
  type(micsttype)  , intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()

  call SoilBGCOneLayer(micfor,micstt,micflx)

  call RunModel_nosalt(micfor,nvars,ystates0l, ystatesfl, err_status)

  call UpdateStateVars(micfor,micstt,micflx,nvars,ystatesfl)
!
  end subroutine RunMicBGC


end module batchmod
