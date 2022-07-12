module batchmod

  use abortutils, only : endrun
implicit none
  private
  public :: BatchModelConfig
  public :: ReadNML
  logical :: Litlayer

  integer :: id_ZNH4B,id_ZNH4S,id_ZNO3B
  integer :: id_ZNO3S,id_H1POB,id_H1PO4
  integer :: id_ZNO2B,id_ZNO2S,id_H2POB
  integer :: id_H2PO4,id_CCO2S,id_CNO2S
  integer :: id_CNO2B,id_CZ2OS,id_Z2OS
  integer :: id_COXYS,id_OXYS,id_COXYG
  integer :: id_CZ2GS,id_CH2GS,id_H2GS
  integer :: id_CCH4G,id_CH4S ,id_ZNFN0
  integer :: id_ZNFNI
  integer :: id_oqc_b,id_oqc_e
  integer :: id_oqn_b,id_oqn_e
  integer :: id_oqp_b,id_oqp_e
  integer :: id_oqa_b,id_oqa_e
  integer :: id_ohc_b,id_ohc_e
  integer :: id_ohn_b,id_ohn_e
  integer :: id_ohp_b,id_ohp_e
  integer :: id_oha_b,id_oha_e
  integer :: id_osc_b,id_osc_e
  integer :: id_osa_b,id_osa_e
  integer :: id_osn_b,id_osn_e
  integer :: id_osp_b,id_osp_e
  integer :: id_orc_b,id_orc_e
  integer :: id_orn_b,id_orn_e
  integer :: id_orp_b,id_orp_e
  integer :: id_omc_b,id_omc_e
  integer :: id_omn_b,id_omn_e
  integer :: id_omp_b,id_omp_e
  integer :: id_omcff_b,id_omcff_e
  integer :: id_omnff_b,id_omnff_e
  integer :: id_ompff_b,id_ompff_e
  integer :: id_ROXYY,id_RNH4Y,id_RNO3Y
  integer :: id_RNO2Y,id_RN2OY,id_RPO4Y
  integer :: id_RP14Y,id_RNHBY,id_RN3BY
  integer :: id_RN2BY,id_RPOBY,id_RP1BY
  integer :: id_ROQCY_b,id_ROQCY_e
  integer :: id_ROQAY_b,id_ROQAY_e

contains

  subroutine ReadNML(nmlfile)

  implicit none
  character(len=*), intent(in) :: nmlfile

  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
  character(len=hist_freq_str_len) :: hist_freq
  namelist /batchmdel/LitLayer,hist_freq

  hist_freq='day'

  inquire (file=nmlfile, iostat=rc)
  if (rc /= 0) then
    write (stdout, '(3a)') 'Error: input file ', trim(nmlfile), &
      ' does not exist.'
    call endrun('stopped in readnml', 25)
  end if

  open (action='read', file=nmlfile, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (stdout, '(2a)') 'Error openning input file "', &
    trim(nmlfile)
    call endrun('stopped in readnml', 32)
  end if

  read(unit=fu, nml=batchmodel, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(stdout,'(a)')"ERROR reading ecosys namelist "
     call endrun('stopped in readnml', 38)
  end if

  close(fu)
  if (.true.) then
    write(stdout, *)
    write(stdout, *) '--------------------'
    write(stdout,ecosys)
    write(stdout, *)
    write(stdout, *) '--------------------'
  endif

  end subroutine ReadNML

! ----------------------------------------------------------------------

  subroutine BatchModelConfig(nvars,ystatesf0l,forc,micfor,micstt,micflx)
!
! DESCRIPTION:
! configure the batch mode of the soil bgc
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod, only : micforctype
  use MicBGCPars, only : micpar
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(in) :: ystatesf0l(nvars)
  type(forc_type), intent(in) :: forc
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx


  real(r8), parameter :: ZERO=1.0E-15_r8
  real(r8), parameter :: ZEROS2=1.0E-08_r8
  integer :: NFGs, jcplx1, JG
  integer :: kk

  NFGs=micpar%NFGs
  jcplx1=micpar%jcplx1
  JG=micpar%jguilds

  micfor%ZERO  =ZERO
  micfor%ZEROS2=ZERO2
  micfor%ZEROS =ZERO

  micfor%CCH4E =CCH4E
  micfor%COXYE =COXYE
  micfor%COXQ  =COXQ
  micfor%COXR  =COXR
  micfor%FLQRI =FLQRI
  micfor%FLQRQ =FLQRQ
  micfor%OFFSET=OFFSET
  micfor%VOLR  =VOLR
  micfor%VOLWRX=VOLWRX
  micfor%VOLY  =VOLY
  micfor%THETY =THETY
  micfor%POROS =POROS
  micfor%FC    =FC
  micfor%TKS   =TKS
  micfor%THETW =THETW
  micfor%PH    =PH
  micfor%BKVL  =BKVL
  micfor%VOLX  =VOLX
  micfor%TFND  =TFND
  micfor%VLNOB =VLNOB
  micfor%VLNO3 =VLNO3
  micfor%VLNH4 =VLNH4
  micfor%VLNHB =VLNHB
  micfor%VLPO4 =VLPO4
  micfor%VLPOB =VLPOB
  micfor%PSISM =PSISM
  micfor%OLSGL =OLSGL
  micfor%ORGC  =ORGC
  micfor%RNO2Y =ystatesf0l(id_RNO2Y)
  micfor%RN2OY =ystatesf0l(id_RN2OY)
  micfor%RN2BY =ystatesf0l(id_RN2BY)
  micfor%ROXYY =ystatesf0l(id_ROXYY)
  micfor%ROXYF =ystatesf0l(id_ROXYF)
  micfor%RNHBY =ystatesf0l(id_RNHBY)
  micfor%RN3BY =ystatesf0l(id_RN3BY)
  micfor%RPOBY =ystatesf0l(id_RPOBY)
  micfor%RP1BY =ystatesf0l(id_RP1BY)
  micfor%ROQCY(0:jcplx1)=ROQCY(0:jcplx1)
  micfor%ROQAY(0:jcplx1)=ROQAY(0:jcplx1)
  micfor%RCH4L =RCH4L
  micfor%ROXYL = ROXYL
  micfor%CFOMC =CFOMC(1:2)
  micfor%litrm=(L==0)
  micfor%Lsurf=.True.
  if(micfor%litrm)then
    micstt%ZNH4TU=AMAX1(0.0,ZNH4S)+AMAX1(0.0,ZNH4B)
    micstt%ZNO3TU=AMAX1(0.0,ZNO3S)+AMAX1(0.0,ZNO3B)
    micstt%H1P4TU=AMAX1(0.0,H1PO4)+AMAX1(0.0,H1POB)
    micstt%H2P4TU=AMAX1(0.0,H2PO4)+AMAX1(0.0,H2POB)
    micstt%CNH4BU=CNH4B
    micstt%CNH4SU=CNH4S
    micstt%CH2P4U=CH2P4
    micstt%CH2P4BU=CH2P4B
    micstt%CH1P4U=CH1P4
    micstt%CH1P4BU=CH1P4B
    micstt%CNO3SU=CNO3S
    micstt%CNO3BU=CNO3B
    micstt%OSC13U=OSC(1,3)
    micstt%OSN13U=OSN(1,3)
    micstt%OSP13U=OSP(1,3)
    micstt%OSC14U=OSC(1,4)
    micstt%OSN14U=OSN(1,4)
    micstt%OSP14U=OSP(1,4)
    micstt%OSC24U=OSC(2,4)
    micstt%OSN24U=OSN(2,4)
    micstt%OSP24U=OSP(2,4)
    micfor%RNH4YU =RNH4Y
    micfor%RNO3YU =RNO3Y
    micfor%RPO4YU =RPO4Y
    micfor%RP14YU =RP14Y
    micfor%VOLWU  =VOLW
    micfor%CFOMCU=CFOMC(1:2)
  else
    micfor%AEC=AEC
    micstt%OXYG=OXYG
  endif
  micstt%CNH4B =CNH4B
  micstt%CNH4S =CNH4S
  micstt%CNO3S =CNO3S
  micstt%CNO3B =CNO3B
  micstt%CH2P4 =CH2P4
  micstt%CH2P4B=CH2P4B
  micstt%CH1P4 =CH1P4
  micstt%CH1P4B=CH1P4B
  micfor%RNH4Y =RNH4Y
  micfor%RNO3Y =RNO3Y
  micfor%RPO4Y =RPO4Y
  micfor%RP14Y =RP14Y
  micfor%VOLW  =VOLW

  if(micfor%Lsurf)then
    micfor%BKVL0=BKVL
  endif
  micfor%DFGS(1:NPH)=DFGS(1:NPH)
  micfor%FILM(1:NPH)=FILM(1:NPH)
  micfor%THETPM(1:NPH)=THETPM(1:NPH)
  micfor%VOLWM(1:NPH)=VOLWM(1:NPH)
  micfor%TORT(1:NPH)=TORT(1:NPH)
  micfor%VOLPM(1:NPH)=VOLPM(1:NPH)

  micstt%EPOC=EPOC
  micstt%EHUM=EHUM
  micstt%ZNH4B=ZNH4B
  micstt%ZNH4S=ZNH4S
  micstt%ZNO3B=ZNO3B
  micstt%ZNO3S=ZNO3S
  micstt%H1POB=H1POB
  micstt%H1PO4=H1PO4
  micstt%ZNO2B=ZNO2B
  micstt%ZNO2S=ZNO2S
  micstt%H2POB=H2POB
  micstt%H2PO4=H2PO4
  micstt%CCO2S=CCO2S
  micstt%CNO2S=CNO2S
  micstt%CNO2B=CNO2B
  micstt%CZ2OS=CZ2OS
  micstt%Z2OS=Z2OS
  micstt%COXYS=COXYS
  micstt%OXYS=OXYS
  micstt%SOXYL=SOXYL
  micstt%COXYG=COXYG
  micstt%CZ2GS=CZ2GS
  micstt%CH2GS=CH2GS
  micstt%H2GS=H2GS
  micstt%CCH4G=CCH4G
  micstt%CH4S=CH4S
  micstt%SCH4L=SCH4L
  micstt%ZNFN0=ZNFN0
  micstt%ZNFNI=ZNFNI

  micstt%FOSRH(0:jcplx1)=FOSRH(0:jcplx1)
  micstt%OQC(0:jcplx1)=OQC(0:jcplx1)
  micstt%OQN(0:jcplx1)=OQN(0:jcplx1)
  micstt%OQP(0:jcplx1)=OQP(0:jcplx1)
  micstt%OQA(0:jcplx1)=OQA(0:jcplx1)
  micstt%OHC(0:jcplx1)=OHC(0:jcplx1)
  micstt%OHN(0:jcplx1)=OHN(0:jcplx1)
  micstt%OHP(0:jcplx1)=OHP(0:jcplx1)
  micstt%OHA(0:jcplx1)=OHA(0:jcplx1)

  micstt%OSC(1:jsken,0:jcplx1)=OSC(1:jsken,0:jcplx1)
  micstt%OSA(1:jsken,0:jcplx1)=OSA(1:jsken,0:jcplx1)
  micstt%OSN(1:jsken,0:jcplx1)=OSN(1:jsken,0:jcplx1)
  micstt%OSP(1:jsken,0:jcplx1)=OSP(1:jsken,0:jcplx1)
  micstt%ORC(1:2,0:jcplx1)=ORC(1:2,0:jcplx1)
  micstt%ORN(1:2,0:jcplx1)=ORN(1:2,0:jcplx1)
  micstt%ORP(1:2,0:jcplx1)=ORP(1:2,0:jcplx1)
  micstt%CNOSC(1:jsken,0:jcplx1)=CNOSC(1:jsken,0:jcplx1)
  micstt%CPOSC(1:jsken,0:jcplx1)=CPOSC(1:jsken,0:jcplx1)
  micstt%OMC(1:3,1:JG,1:NFGs,0:jcplx1)=OMC(1:3,1:JG,1:NFGs,0:jcplx1)
  micstt%OMN(1:3,1:JG,1:NFGs,0:jcplx1)=OMN(1:3,1:JG,1:NFGs,0:jcplx1)
  micstt%OMP(1:3,1:JG,1:NFGs,0:jcplx1)=OMP(1:3,1:JG,1:NFGs,0:jcplx1)
  micstt%OMCff(1:3,1:JG,1:NFGs)=OMCff(1:3,1:JG,1:NFGs)
  micstt%OMNff(1:3,1:JG,1:NFGs)=OMNff(1:3,1:JG,1:NFGs)
  micstt%OMPff(1:3,1:JG,1:NFGs)=OMPff(1:3,1:JG,1:NFGs)

  micflx%RVMXC=RVMXC
  micflx%RVMBC=RVMBC
  micflx%RINHOff(1:JG,1:NFGs)=RINHOff(1:JG,1:NFGs)
  micflx%RINHBff(1:JG,1:NFGs)=RINHBff(1:JG,1:NFGs)
  micflx%RINOOff(1:JG,1:NFGs)=RINOOff(1:JG,1:NFGs)
  micflx%RINOBff(1:JG,1:NFGs)=RINOBff(1:JG,1:NFGs)
  micflx%RIPOOff(1:JG,1:NFGs)=RIPOOff(1:JG,1:NFGs)
  micflx%RIPBOff(1:JG,1:NFGs)=RIPBOff(1:JG,1:NFGs)
  micflx%RIPO1ff(1:JG,1:NFGs)=RIPO1ff(1:JG,1:NFGs)
  micflx%RIPB1ff(1:JG,1:NFGs)=RIPB1ff(1:JG,1:NFGs)
  micflx%ROXYSff(1:JG,1:NFGs)=ROXYSff(1:JG,1:NFGs)

  micflx%RINHO(1:JG,1:NFGs,0:JCPLX1)=RINHO(1:JG,1:NFGs,0:JCPLX1)
  micflx%RINHB(1:JG,1:NFGs,0:JCPLX1)=RINHB(1:JG,1:NFGs,0:JCPLX1)
  micflx%RINOO(1:JG,1:NFGs,0:JCPLX1)=RINOO(1:JG,1:NFGs,0:JCPLX1)
  micflx%RINOB(1:JG,1:NFGs,0:JCPLX1)=RINOB(1:JG,1:NFGs,0:JCPLX1)
  micflx%RIPOO(1:JG,1:NFGs,0:JCPLX1)=RIPOO(1:JG,1:NFGs,0:JCPLX1)
  micflx%RIPBO(1:JG,1:NFGs,0:JCPLX1)=RIPBO(1:JG,1:NFGs,0:JCPLX1)
  micflx%RIPO1(1:JG,1:NFGs,0:JCPLX1)=RIPO1(1:JG,1:NFGs,0:JCPLX1)
  micflx%RIPB1(1:JG,1:NFGs,0:JCPLX1)=RIPB1(1:JG,1:NFGs,0:JCPLX1)
  micflx%ROXYS(1:JG,1:NFGs,0:JCPLX1)=ROXYS(1:JG,1:NFGs,0:JCPLX1)

  end subroutine BatchModelConfig

! ----------------------------------------------------------------------
  subroutine ReadForc(forc)
!
! DESCRIPTION
! read forcing data
  use MicForcTypeMod, only : micforctype
  implicit none
  type(forc_type), intent(in) :: forc



  end subroutine ReadForc
! ----------------------------------------------------------------------


  subroutine Initboxbgc(nvars)
  implicit none

  integer, intent(out) :: nvars
  integer, itemp
  itemp=0
  id_ZNH4B=addone(itemp)
  id_ZNH4S=addone(itemp)
  id_ZNO3B=addone(itemp)
  id_ZNO3S=addone(itemp)
  id_H1POB=addone(itemp)
  id_H1PO4=addone(itemp)
  id_ZNO2B=addone(itemp)
  id_ZNO2S=addone(itemp)
  id_H2POB=addone(itemp)
  id_H2PO4=addone(itemp)
  id_CCO2S=addone(itemp)
  id_CNO2S=addone(itemp)
  id_CNO2B=addone(itemp)
  id_CZ2OS=addone(itemp)
  id_Z2OS =addone(itemp)
  id_COXYS=addone(itemp)
  id_OXYS =addone(itemp)
  id_COXYG=addone(itemp)
  id_CZ2GS=addone(itemp)
  id_CH2GS=addone(itemp)
  id_H2GS =addone(itemp)
  id_CCH4G=addone(itemp)
  id_CH4S =addone(itemp)
  id_ZNFN0=addone(itemp)
  id_ZNFNI=addone(itemp)

  id_oqc_b=itemp+1
  id_oqc_e=id_oqc_b+jcplx1;itemp=id_oqc_e
  id_oqn_b=itemp+1
  id_oqn_e=id_oqn_b+jcplx1;itemp=id_oqn_e
  id_oqp_b=itemp+1
  id_oqp_e=id_oqp_b+jcplx1;itemp=id_oqp_e
  id_oqa_b=itemp+1
  id_oqa_e=id_oqa_b+jcplx1;itemp=id_oqa_e
  id_ohc_b=itemp+1
  id_ohc_e=id_ohc_b+jcplx1;itemp=id_ohc_e
  id_ohn_b=itemp+1
  id_ohn_e=id_ohn_b+jcplx1;itemp=id_ohn_e
  id_ohp_b=itemp+1
  id_ohp_e=id_ohp_b+jcplx1;itemp=id_ohp_e
  id_oha_b=itemp+1
  id_oha_e=id_oha_b+jcplx1;itemp=id_oha_e
  id_osc_b=itemp+1
  id_osc_e=id_osc_b+jsken*jcplx-1;itemp=id_osc_e
  id_osa_b=itemp+1
  id_osa_e=id_osa_b+jsken*jcplx-1;itemp=id_osa_e
  id_osn_b=itemp+1
  id_osn_e=id_osn_b+jsken*jcplx-1;itemp=id_osn_e
  id_osp_b=itemp+1
  id_osp_e=id_osp_b+jsken*jcplx-1;itemp=id_osp_e
  id_orc_b=itemp+1
  id_orc_e=id_orc_b+2*jcplx-1;itemp=id_orc_e
  id_orn_b=itemp+1
  id_orn_e=id_orn_b+2*jcplx-1;itemp=id_orn_e
  id_orp_b=itemp+1
  id_orp_e=id_orp_b+2*jcplx-1;itemp=id_orp_e
  id_omc_b=itemp+1
  id_omc_e=id_omc_b+3*JG*NFGs*jcplx-1;itemp=id_omc_e
  id_omn_b=itemp+1
  id_omn_e=id_omn_b+3*JG*NFGs*jcplx-1;itemp=id_omn_e
  id_omp_b=itemp+1
  id_omp_e=id_omp_b+3*JG*NFGs*jcplx-1;itemp=id_omp_e
  id_omcff_b=itemp+1
  id_omcff_e=id_omcff_b+3*JG*NFGs-1;itemp=id_omcff_e
  id_omnff_b=itemp+1
  id_omnff_e=id_omnff_b+3*JG*NFGs-1;itemp=id_omnff_e
  id_ompff_b=itemp+1
  id_ompff_e=id_ompff_b+3*JG*NFGs-1;itemp=id_ompff_e
  id_ROXYY=addone(itemp)
  id_RNH4Y=addone(itemp)
  id_RNO3Y=addone(itemp)
  id_RNO2Y=addone(itemp)
  id_RN2OY=addone(itemp)
  id_RPO4Y=addone(itemp)
  id_RP14Y=addone(itemp)
  id_RNHBY=addone(itemp)
  id_RN3BY=addone(itemp)
  id_RN2BY=addone(itemp)
  id_RPOBY=addone(itemp)
  id_RP1BY=addone(itemp)
  id_ROQCY_b=itemp+1
  id_ROQCY_e=id_ROQCY_b+jcplx1;itemp=id_ROQCY_e
  id_ROQAY_b=itemp+1
  id_ROQAY_e=id_ROQAY_b+jcplx1;itemp=id_ROQAY_e
  nvars=itemp
  end subroutine Initboxbgc
! ----------------------------------------------------------------------

  subroutine UpdateStateVars(micfor,micstt,micflx,ystatesfl)

  implicit none

  real(r8) :: DC,DN,DP,OC,ON,OP
  integer :: K,N,NGL,M
  associate(         &
    VOLW  => micfor%VOLW, &
  )
!atmospheric gaseous CO2,CH4,O2,NH3,N2,N2O,H2

  ystatesfl(id_ZNH4B)=micstt%ZNH4B+
  ystatesfl(id_ZNH4S)=micstt%ZNH4S+
  ystatesfl(id_ZNO3B)=micstt%ZNO3B+
  ystatesfl(id_ZNO3S)=micstt%ZNO3S+
  ystatesfl(id_H1POB)=micstt%H1POB+
  ystatesfl(id_H1PO4)=micstt%H1PO4+
  ystatesfl(id_ZNO2B)=micstt%ZNO2B+
  ystatesfl(id_ZNO2S)=micstt%ZNO2S+
  ystatesfl(id_H2POB)=micstt%H2POB+
  ystatesfl(id_H2PO4)=micstt%H2PO4+
  ystatesfl(id_CCO2S)=micstt%CCO2S-RCO2O/VOLW
  ystatesfl(id_CNO2S)=micstt%ZNO2S/(VOLW*micfor%VLNO3)
  ystatesfl(id_CNO2B)=micstt%ZNO2B/(VOLW*micfor%VLNOB)
  ystatesfl(id_CZ2OS)=micstt%Z2OS/VOLW
  ystatesfl(id_Z2OS) =micstt%Z2OS
  ystatesfl(id_COXYS)=micstt%OXYS/VOLW
  ystatesfl(id_OXYS) =micstt%OXYS-RUPOXO
  ystatesfl(id_COXYG)=micstt%COXYG
  ystatesfl(id_CZ2GS)=micstt%CZ2GS
  ystatesfl(id_CH2GS)=micstt%CH2GS
  ystatesfl(id_H2GS) =micstt%H2GS
  ystatesfl(id_CCH4G)=micstt%CCH4G
  ystatesfl(id_CH4S) =micstt%CH4S-RCH4O
  ystatesfl(id_ZNFN0)=micstt%ZNFN0
  ystatesfl(id_ZNFNI)=micstt%ZNFNI

! the following variables are updated in the microbial model
  ystatesfl(id_oqc_b:id_oqc_e)=micstt%OQC(0:jcplx1)
  ystatesfl(id_oqn_b:id_oqn_e)=micstt%OQN(0:jcplx1)
  ystatesfl(id_oqp_b:id_oqp_e)=micstt%OQP(0:jcplx1)
  ystatesfl(id_oqa_b:id_oqa_e)=micstt%OQA(0:jcplx1)
  ystatesfl(id_ohc_b:id_ohc_e)=micstt%OHC(0:jcplx1)
  ystatesfl(id_ohn_b:id_ohn_e)=micstt%OHN(0:jcplx1)
  ystatesfl(id_ohp_b:id_ohp_e)=micstt%OHP(0:jcplx1)
  ystatesfl(id_oha_b:id_oha_e)=micstt%OHA(0:jcplx1)
  ystatesfl(id_osc_b:id_osc_e)=micstt%OSC(1:jsken,0:jcplx1)
  ystatesfl(id_osa_b:id_osa_e)=micstt%OSA(1:jsken,0:jcplx1)
  ystatesfl(id_osn_b:id_osn_e)=micstt%OSN(1:jsken,0:jcplx1)
  ystatesfl(id_osp_b:id_osp_e)=micstt%OSP(1:jsken,0:jcplx1)
  ystatesfl(id_orc_b:id_orc_e)=reshape(micstt%ORC(1:2,0:jcplx1),&
    shape(size(micstt%ORC(1:2,0:jcplx1))))
  ystatesfl(id_orn_b:id_orn_e)=reshape(micstt%ORN(1:2,0:jcplx1),&
    shape(size(micstt%ORN(1:2,0:jcplx1))))
  ystatesfl(id_orp_b:id_orp_e)=reshape(micstt%ORP(1:2,0:jcplx1),&
    shape(size(micstt%ORP(1:2,0:jcplx1))))
  ystatesfl(id_omc_b:id_omc_e)=reshape(micstt%OMC(1:3,1:JG,1:NFGs,0:jcplx1),&
    shape(size(micstt%OMC(1:3,1:JG,1:NFGs,0:jcplx1))))
  ystatesfl(id_omn_b:id_omn_e)=reshape(micstt%OMN(1:3,1:JG,1:NFGs,0:jcplx1),&
    shape(size(micstt%OMN(1:3,1:JG,1:NFGs,0:jcplx1))))
  ystatesfl(id_omp_b:id_omp_e)=reshape(micstt%OMP(1:3,1:JG,1:NFGs,0:jcplx1),&
    shape(size(micstt%OMP(1:3,1:JG,1:NFGs,0:jcplx1))))
  ystatesfl(id_omcff_b:id_omcff_e)=reshape(micstt%OMCff(1:3,1:JG,1:NFGs),&
    shape(size(micstt%OMCff(1:3,1:JG,1:NFGs))))
  ystatesfl(id_omnff_b:id_omnff_e)=reshape(micstt%OMNff(1:3,1:JG,1:NFGs),&
    shape(size(micstt%OMNff(1:3,1:JG,1:NFGs))))
  ystatesfl(id_ompff_b:id_ompff_e)=reshape(micstt%OMPff(1:3,1:JG,1:NFGs),&
    shape(size(micstt%OMPff(1:3,1:JG,1:NFGs))))

! summarize diagnostic fluxes
  DO K=0,jcplx1
    IF(not.micfor%litrm.or.(K.NE.micpar%k_POM.AND.K.NE.micpar%k_humus))THEN
      DO N=1,NFGs
        DO NGL=1,JG
          ystatesfl(id_ROXYY)=ystatesfl(id_ROXYY)+ROXYS(NGL,N,K)
          ystatesfl(id_RNH4Y)=ystatesfl(id_RNH4Y)+RVMX4(NGL,N,K)+RINHO(NGL,N,K)
          ystatesfl(id_RNO3Y)=ystatesfl(id_RNO3Y)+RVMX3(NGL,N,K)+RINOO(NGL,N,K)
          ystatesfl(id_RNO2Y)=ystatesfl(id_RNO2Y)+RVMX2(NGL,N,K)
          ystatesfl(id_RN2OY)=ystatesfl(id_RN2OY)+RVMX1(NGL,N,K)
          ystatesfl(id_RPO4Y)=ystatesfl(id_RPO4Y)+RIPOO(NGL,N,K)
          ystatesfl(id_RP14Y)=ystatesfl(id_RP14Y)+RIPO1(NGL,N,K)
          ystatesfl(id_RNHBY)=ystatesfl(id_RNHBY)+RVMB4(NGL,N,K)+RINHB(NGL,N,K)
          ystatesfl(id_RN3BY)=ystatesfl(id_RN3BY)+RVMB3(NGL,N,K)+RINOB(NGL,N,K)
          ystatesfl(id_RN2BY)=ystatesfl(id_RN2BY)+RVMB2(NGL,N,K)
          ystatesfl(id_RPOBY)=ystatesfl(id_RPOBY)+RIPBO(NGL,N,K)
          ystatesfl(id_RP1BY)=ystatesfl(id_RP1BY)+RIPB1(NGL,N,K)
          ystatesfl(id_ROQCY_b+K)=ystatesfl(id_ROQCY_b+K)+ROQCS(NGL,N,K)
          ystatesfl(id_ROQAY_b+K)=ystatesfl(id_ROQAY_b+K)+ROQAS(NGL,N,K)
        enddo
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NFGs
    DO NGL=1,JG
      ystatesfl(id_ROXYY)=ystatesfl(id_ROXYY)+ROXYSff(NGL,N)
      ystatesfl(id_RNH4Y)=ystatesfl(id_RNH4Y)+RVMX4ff(NGL,N)+RINHOff(NGL,N)
      ystatesfl(id_RNO3Y)=ystatesfl(id_RNO3Y)+RVMX3ff(NGL,N)+RINOOff(NGL,N)
      ystatesfl(id_RNO2Y)=ystatesfl(id_RNO2Y)+RVMX2ff(NGL,N)
      ystatesfl(id_RN2OY)=ystatesfl(id_RN2OY)+RVMX1ff(NGL,N)
      ystatesfl(id_RPO4Y)=ystatesfl(id_RPO4Y)+RIPOOff(NGL,N)
      ystatesfl(id_RP14Y)=ystatesfl(id_RP14Y)+RIPO1ff(NGL,N)
      ystatesfl(id_RNHBY)=ystatesfl(id_RNHBY)+RVMB4ff(NGL,N)+RINHBff(NGL,N)
      ystatesfl(id_RN3BY)=ystatesfl(id_RN3BY)+RVMB3ff(NGL,N)+RINOBff(NGL,N)
      ystatesfl(id_RN2BY)=ystatesfl(id_RN2BY)+RVMB2ff(NGL,N)
      ystatesfl(id_RPOBY)=ystatesfl(id_RPOBY)+RIPBOff(NGL,N)
      ystatesfl(id_RP1BY)=ystatesfl(id_RP1BY)+RIPB1ff(NGL,N)
    enddo
  ENDDO

  ystatesfl(id_RNO2Y)=ystatesfl(id_RNO2Y)+RVMXC
  ystatesfl(id_RN2BY)=ystatesfl(id_RN2BY)+RVMBC

  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8

  DO K=0,jcplx1
    IF(K.LE.2)THEN
      DO N=1,NFGs
        DO  M=1,3
          DO NGL=1,JG
            DC=DC+OMC(M,NGL,N,K)
            DN=DN+OMN(M,NGL,N,K)
            DP=DP+OMP(M,NGL,N,K)
          ENDDO
        enddo
      ENDDO
    ELSE
      DO N=1,NFGs
        DO  M=1,3
          DO NGL=1,JG
            OC=OC+OMC(M,NGL,N,K)
            ON=ON+OMN(M,NGL,N,K)
            OP=OP+OMP(M,NGL,N,K)
          enddo
        enddo
      ENDDO
    ENDIF
  ENDDO
! abstract complex
  DO  N=1,NFGs
    DO  M=1,3
      DO NGL=1,JG
        OC=OC+OMCff(M,NGL,N)
        ON=ON+OMNff(M,NGL,N)
        OP=OP+OMPff(M,NGL,N)
      enddo
    enddo
  ENDDO
! microbial residue
  DO K=0,jcplx1
    IF(K.LE.2)THEN
      DO M=1,2
        DC=DC+ORC(M,K)
        DN=DN+ORN(M,K)
        DP=DP+ORP(M,K)
      ENDDO
      DC=DC+OQC(K)+OQCH(K)+OHC(K)+OQA(K)+OQAH(K)+OHA(K)
      DN=DN+OQN(K)+OQNH(K)+OHN(K)
      DP=DP+OQP(K)+OQPH(K)+OHP(K)
      DO M=1,jsken
        DC=DC+OSC(M,K)
        DN=DN+OSN(M,K)
        DP=DP+OSP(M,K)
      ENDDO
    ELSE
      DO M=1,2
        OC=OC+ORC(M,K)
        ON=ON+ORN(M,K)
        OP=OP+ORP(M,K)
      ENDDO
      OC=OC+OQC(K)+OQCH(K)+OHC(K)+OQA(K)+OQAH(K)+OHA(K)
      ON=ON+OQN(K)+OQNH(K)+OHN(K)
      OP=OP+OQP(K)+OQPH(K)+OHP(K)
      DO M=1,jsken
        OC=OC+OSC(M,K)
        ON=ON+OSN(M,K)
        OP=OP+OSP(M,K)
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

  function addone(itemp)result(ans)
!
!  DESCRIPTION
! increase itemp by one
  implicit none
  integer, intent(inout) :: itemp

  integer :: ans

  itemp=itemp+1
  ans=itemp
  end function addone
! ----------------------------------------------------------------------

  subroutine getvarlist(nvars, varl, unitl, vartypes)

  use histMod, only : hist_var_str_len,hist_unit_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)

  integer :: iknen,icplx
  integer :: jj,ll,k,m,n,ngl
  vartypes(:)=var_state_type
  varl(id_ZNH4B)='NH4p1_mass_band';unitl(id_ZNH4B)
  varl(id_ZNH4S)='NH4p1_mass';unitl(id_ZNH4S)
  varl(id_ZNO3B)='NO3e1_mass_band';unitl(id_ZNO3B)
  varl(id_ZNO3S)='NO3e1_mass';unitl(id_ZNO3S)
  varl(id_H1POB)='H1PO4e2_mass_band';unitl(id_H1POB)
  varl(id_H1PO4)='H1PO4e2_mass';unitl(id_H1PO4)
  varl(id_ZNO2B)='NO2e1_mass_band';unitl(id_ZNO2B)
  varl(id_ZNO2S)='NO2e1_mass';unitl(id_ZNO2S)
  varl(id_H2POB)='H2PO4e1_mass_band';unitl(id_H2POB)
  varl(id_H2PO4)='H2PO4e1_mass';unitl(id_H2PO4)
  varl(id_CCO2S)='CO2aq_conc';unitl(id_CCO2S)
  varl(id_CNO2S)='NO2e1_conc';unitl(id_CNO2S)
  varl(id_CNO2B)='NO2e1_conc_band';unitl(id_CNO2B)
  varl(id_CZ2OS)='N2Oaq_conc';unitl(id_CZ2OS)
  varl(id_Z2OS) ='N2Oaq_mass';unitl(id_Z2OS)
  varl(id_COXYS)='O2aq_conc';unitl(id_COXYS)
  varl(id_OXYS) ='O2aq_mass';unitl(id_OXYS)
  varl(id_COXYG)='O2gs_conc';unitl(id_COXYG)
  varl(id_CZ2GS)='N2gs_conc';unitl(id_CZ2GS)
  varl(id_CH2GS)='H2gs_conc';unitl(id_CH2GS)
  varl(id_H2GS) ='H2gs_mass';unitl(id_H2GS)
  varl(id_CCH4G)='CH4gs_conc';unitl(id_CCH4G)
  varl(id_CH4S) ='CH4aq_conc';unitl(id_CH4S)
  varl(id_ZNFN0)='NitIhb_0';unitl(id_ZNFN0)
  varl(id_ZNFNI)='NitIhb_t';unitl(id_ZNFNI)
  do jj=id_oqc_b,id_oqc_e
    varl(jj)='DOC_mass_cplx'//trim(micpar%cplxname(jj-id_oqc_b))
    unitl(jj)=
  enddo
  do jj=id_oqn_b,id_oqn_e
    varl(jj)='DON_mass_cplx'//trim(micpar%cplxname(jj-id_oqn_b))
    unitl(jj)=
  enddo
  do jj=id_oqp_b,id_oqp_e
    varl(jj)='DOP_mass_cplx'//trim(micpar%cplxname(jj-id_oqp_b))
    unitl(jj)=
  enddo
  do jj=id_oqa_b,id_oqa_e
    varl(jj)='Acetate_mass_cplx'//trim(micpar%cplxname(jj-id_oqa_b))
    unitl(jj)=
  enddo
  do jj=id_ohc_b,id_ohc_e
    varl(jj)='SORBC_mass_cplx'//trim(micpar%cplxname(jj-id_ohc_b))
    unitl(jj)=
  enddo
  do jj=id_ohn_b,id_ohn_e
    varl(jj)='SORBN_mass_cplx'//trim(micpar%cplxname(jj-id_ohn_b))
    unitl(jj)=
  enddo
  do jj=id_ohp_b,id_ohp_e
    varl(jj)='SORBP_mass_cplx'//trim(micpar%cplxname(jj-id_ohp_b))
    unitl(jj)=
  enddo
  do jj=id_oha_b,id_oha_e
    varl(jj)='SORBAcet_mass_cplx'//trim(micpar%cplxname(jj-id_oha_b))
    unitl(jj)=
  enddo
  do jj=id_osc_b,id_osc_e
    iknen=iknen
    icplx=floor((jj-1.e-8_r8)/jsken)
    varl(jj)='OSC_mass_'//trim(micpar%kiname(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo
  do jj=id_osn_b,id_osn_e
    iknen=iknen
    icplx=floor((jj-1.e-8_r8)/jsken)
    varl(jj)='OSN_mass_'//trim(micpar%kiname(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo

  do jj=id_osp_b,id_osp_e
    iknen=iknen
    icplx=floor((jj-1.e-8_r8)/jsken)
    varl(jj)='OSP_mass_'//trim(micpar%kiname(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo

  do jj=id_osa_b,id_osa_e
    iknen=iknen
    icplx=floor((jj-1.e-8_r8)/jsken)
    varl(jj)='OSA_mass_'//trim(micpar%kiname(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo

  do jj=id_orc_b,id_orc_e
    iknen=mod(jj,2)
    icplx=floor((jj-1.e-8_r8)/2)
    varl(jj)='mbresdC_mass_'//trim(micpar%micresb(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo

  do jj=id_orn_b,id_orn_e
    iknen=mod(jj,2)
    icplx=floor((jj-1.e-8_r8)/2)
    varl(jj)='mbresdN_mass_'//trim(micpar%micresb(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo

  do jj=id_orp_b,id_orp_e
    iknen=mod(jj,2)
    icplx=floor((jj-1.e-8_r8)/2)
    varl(jj)='mbresdP_mass_'//trim(micpar%micresb(iknen))//trim(micpar%cplxname(icplx))
    unitl(jj)=
  enddo


  jj=0
  DO k=0,jcplx1
  DO N=1,NFGs
  DO NGL=1,JG
  DO M=1,3
    ll=id_omc_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMC_mass_'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    unitl(ll)=
    ll=id_omn_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMN_mass_'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    unitl(ll)=
    ll=id_omp_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMP_mass_'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    unitl(ll)=
    jj=jj+1
  enddo
  ENDDO
  ENDDO
  ENDDO

  jj=0
  DO N=1,NFGs
  DO NGL=1,JG
  DO M=1,3
    ll=id_omcff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMC_mass_'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%amicname(N))
    unitl(ll)=
    ll=id_omnff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMN_mass_'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%amicname(N))
    unitl(ll)=
    ll=id_ompff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMP_mass_'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%amicname(N))
    unitl(ll)=
    jj=jj+1
  enddo
  ENDDO
  ENDDO

  varl(id_ROXYY)='oxygen_upatke'; unitl(id_ROXYY); vartypes(id_ROXYY)=var_flux_type
  varl(id_RNH4Y)='NH4p1_uptake'; unitl(id_RNH4Y); vartypes(id_RNH4Y)=var_flux_type
  varl(id_RNO3Y)='NO3e1_uptake'; unitl(id_RNO3Y); vartypes(id_RNO3Y)=var_flux_type
  varl(id_RNO2Y)='NO2e1_uptake'; unitl(id_RNO2Y); vartypes(id_RNO2Y)=var_flux_type
  varl(id_RN2OY)='N2O_uptake'; unitl(id_RN2OY); vartypes(id_RN2OY)=var_flux_type
  varl(id_RPO4Y)='PO4_uptake'; unitl(id_RPO4Y); vartypes(id_RPO4Y)=var_flux_type
  varl(id_RP14Y)='HPO4_demand'; unitl(id_RP14Y); vartypes(id_RP14Y)=var_flux_type
  varl(id_RNHBY)='NH4p1_uptake_band'; unitl(id_RNHBY); vartypes(id_RNHBY)=var_flux_type
  varl(id_RN3BY)='NO3e1_uptake_band'; unitl(id_RN3BY); vartypes(id_RN3BY)=var_flux_type
  varl(id_RN2BY)='NO2e1_uptake_band'; unitl(id_RN2BY); vartypes(id_RN2BY)=var_flux_type
  varl(id_RPOBY)='PO4_uptake_band'; unitl(id_RPOBY); vartypes(id_RPOBY)=var_flux_type
  varl(id_RP1BY)='HPO4_uptake_band'; unitl(id_RP1BY); vartypes(id_RP1BY)=var_flux_type
  do jj =id_ROQCY_b,id_ROQCY_e
    varl(jj)='DOC_uptake_'//micpar%cplxname(jj-id_ROQCY_b)
    vartypes(jj)=var_flux_type
    unitl(jj)
  enddo
  do jj =id_ROQAY_b,id_ROQAY_e
    varl(jj)='acetate_uptake_'//micpar%cplxname(jj-id_ROQCY_b)
    vartypes(jj)=var_flux_type
    unitl(jj)
  enddo
  end subroutine getvarlist
! ----------------------------------------------------------------------

end module batchmod
