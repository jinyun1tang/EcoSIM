module MicBGCAPI

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use NitrosMod      , only : VerticalLitterMixLvsLL
  use NitroDisturbMod, only : SOMRemovalByDisturbance
  use EcoSIMSolverPar
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod, only : micforctype
  use minimathmod, only : AZMAX1
  use SoilBGCDataType
  USE PlantDataRateType
  USE SoilWaterDataType
  USE SurfLitterDataType
  USE EcosimBGCFluxType
  USE EcoSIMCtrlDataType
  USE ClimForcDataType
  USE GridDataType
  USE SoilPropertyDataType
  use SoilPhysDataType
  use SOMDataType
  use ChemTranspDataType
  use MicrobialDataType
  use IrrigationDataType
  use SoilHeatDataType
  use MicBGCPars, only : micpar
  use MicBGCMod, only : SoilBGCOneLayer
implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = __FILE__

  type(micforctype) :: micfor
  type(micsttype) :: micstt
  type(micfluxtype) :: micflx
  integer :: curI,curJ


  public :: MicrobeModel
  public :: MicAPI_Init
  public :: MicAPI_cleanup
  contains

!------------------------------------------------------------------------------------------

  subroutine MicAPI_Init()

  implicit none

  call micfor%Init()
  call micstt%Init()
  call micflx%Init()

  end subroutine MicAPI_Init
!------------------------------------------------------------------------------------------
  subroutine MicAPI_cleanup()

  implicit none


  call micfor%destroy()
  call micflx%destroy()
  call micstt%destroy()
  end subroutine MicAPI_cleanup
!------------------------------------------------------------------------------------------
  subroutine MicrobeModel(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOIL BIOLOGICAL TRANSFORMATIONS
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: L,NX,NY

!   begin_execution
  curI=I; curJ=J
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!       VOLWZ=water volume used to calculate aqueous microbial
!       concentrations that drive microbial density effects on
!       decomposition
      D998: DO L=0,NL(NY,NX)
        IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
             call MicBGC1Layer(I,J,L,NY,NX)
          ELSE
            RCO2O(L,NY,NX)=0.0_r8
            RCH4O(L,NY,NX)=0.0_r8
            RH2GO(L,NY,NX)=0.0_r8
            RUPOXO(L,NY,NX)=0.0_r8
            RN2G(L,NY,NX)=0.0_r8
            RN2O(L,NY,NX)=0.0_r8
            XNH4S(L,NY,NX)=0.0_r8
            XNO3S(L,NY,NX)=0.0_r8
            XNO2S(L,NY,NX)=0.0_r8
            XH2PS(L,NY,NX)=0.0_r8
            XH1PS(L,NY,NX)=0.0_r8
            XNH4B(L,NY,NX)=0.0_r8
            XNO3B(L,NY,NX)=0.0_r8
            XNO2B(L,NY,NX)=0.0_r8
            XH2BS(L,NY,NX)=0.0_r8
            XH1BS(L,NY,NX)=0.0_r8
            XN2GS(L,NY,NX)=0.0_r8
          ENDIF
  !     MIX LITTER C BETWEEN ADJACENT SOIL LAYERS L AND LL
          call VerticalLitterMixLvsLL(I,J,L,NY,NX)

        ELSE
          RCO2O(L,NY,NX)=0.0_r8
          RCH4O(L,NY,NX)=0.0_r8
          RH2GO(L,NY,NX)=0.0_r8
          RUPOXO(L,NY,NX)=0.0_r8
          RN2G(L,NY,NX)=0.0_r8
          RN2O(L,NY,NX)=0.0_r8
          XNH4S(L,NY,NX)=0.0_r8
          XNO3S(L,NY,NX)=0.0_r8
          XNO2S(L,NY,NX)=0.0_r8
          XH2PS(L,NY,NX)=0.0_r8
          XH1PS(L,NY,NX)=0.0_r8
          XNH4B(L,NY,NX)=0.0_r8
          XNO3B(L,NY,NX)=0.0_r8
          XNO2B(L,NY,NX)=0.0_r8
          XH2BS(L,NY,NX)=0.0_r8
          XH1BS(L,NY,NX)=0.0_r8
          XN2GS(L,NY,NX)=0.0_r8
        ENDIF
      ENDDO D998
!
!       SOC LOSS IF FIRE OR REMOVAL EVENT IS ENTERED IN DISTURBANCE FILE
!
      call SOMRemovalByDisturbance(I,J,NY,NX)
    ENDDO D9990
  ENDDO D9995
  RETURN
  END subroutine MicrobeModel

!------------------------------------------------------------------------------------------

  subroutine MicBGC1Layer(I,J,L,NY,NX)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX

  call micflx%ZeroOut()

  call MicAPISend(L,NY,NX,micfor,micstt,micflx)

  call SoilBGCOneLayer(micfor,micstt,micflx)

  call MicAPIRecv(L,NY,NX,micfor%litrm,micstt,micflx)

  end subroutine Micbgc1Layer
!------------------------------------------------------------------------------------------

  subroutine MicAPISend(L,NY,NX,micfor,micstt,micflx)
  implicit none
  integer, intent(in) :: L,NY,NX
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx

  integer :: NFGs, jcplx1, JG, k_POM, k_humus
  integer :: kk, ndbiomcp, nlbiomcp, NMICBSA, NMICBSO
  NFGs=micpar%NFGs
  jcplx1=micpar%jcplx1
  JG=micpar%jguilds
  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp
  NMICBSA  = micpar%NMICBSA
  NMICBSO  = micpar%NMICBSO
  k_humus  = micpar%k_humus
  k_POM    = micpar%k_POM

  micfor%ZERO  =ZERO
  micfor%CCH4E =CCH4E(NY,NX)
  micfor%COXYE =COXYE(NY,NX)
  micfor%COXQ  =COXQ(NY,NX)
  micfor%COXR  =COXR(NY,NX)
  micfor%FLQRI =FLQRI(NY,NX)
  micfor%FLQRQ =FLQRQ(NY,NX)
  micfor%OFFSET=OFFSET(NY,NX)
  micfor%VOLR  =VOLR(NY,NX)
  micfor%VOLWRX=VOLWRX(NY,NX)
  micfor%ZEROS2=ZEROS2(NY,NX)
  micfor%ZEROS =ZEROS(NY,NX)
  micfor%VOLY  =VOLY(L,NY,NX)
  micfor%THETY =THETY(L,NY,NX)
  micfor%POROS =POROS(L,NY,NX)
  micfor%FC    =FC(L,NY,NX)
  micfor%TKS   =TKS(L,NY,NX)
  micfor%THETW =THETW(L,NY,NX)
  micfor%PH    =PH(L,NY,NX)
  micfor%BKVL  =BKVL(L,NY,NX)
  micfor%VOLX  =VOLX(L,NY,NX)
  micfor%TFND  =TFND(L,NY,NX)
  micfor%VLNOB =VLNOB(L,NY,NX)
  micfor%VLNO3 =VLNO3(L,NY,NX)
  micfor%VLNH4 =VLNH4(L,NY,NX)
  micfor%VLNHB =VLNHB(L,NY,NX)
  micfor%VLPO4 =VLPO4(L,NY,NX)
  micfor%VLPOB =VLPOB(L,NY,NX)
  micfor%PSISM =PSISM(L,NY,NX)
  micfor%OLSGL =OLSGL(L,NY,NX)
  micfor%ORGC  =ORGC(L,NY,NX)
  micfor%RNO2Y =RNO2Y(L,NY,NX)
  micfor%RN2OY =RN2OY(L,NY,NX)
  micfor%RN2BY =RN2BY(L,NY,NX)
  micfor%ROXYY =ROXYY(L,NY,NX)
  micfor%ROXYF =ROXYF(L,NY,NX)
  micfor%RCH4L =RCH4L(L,NY,NX)
  micfor%RNHBY =RNHBY(L,NY,NX)
  micfor%RN3BY =RN3BY(L,NY,NX)
  micfor%RPOBY =RPOBY(L,NY,NX)
  micfor%RP1BY =RP1BY(L,NY,NX)
  micfor%ROXYL =ROXYL(L,NY,NX)
  micfor%ROQCY(0:jcplx1)=ROQCY(0:jcplx1,L,NY,NX)
  micfor%ROQAY(0:jcplx1)=ROQAY(0:jcplx1,L,NY,NX)
  micfor%litrm=(L==0)
  micfor%Lsurf=(L==NU(NY,NX))
  if(micfor%litrm)then
    micstt%ZNH4TU=AZMAX1(ZNH4S(NU(NY,NX),NY,NX))+AZMAX1(ZNH4B(NU(NY,NX),NY,NX))
    micstt%ZNO3TU=AZMAX1(ZNO3S(NU(NY,NX),NY,NX))+AZMAX1(ZNO3B(NU(NY,NX),NY,NX))
    micstt%H1P4TU=AZMAX1(H1PO4(NU(NY,NX),NY,NX))+AZMAX1(H1POB(NU(NY,NX),NY,NX))
    micstt%H2P4TU=AZMAX1(H2PO4(NU(NY,NX),NY,NX))+AZMAX1(H2POB(NU(NY,NX),NY,NX))
    micstt%CNH4BU=CNH4B(NU(NY,NX),NY,NX)
    micstt%CNH4SU=CNH4S(NU(NY,NX),NY,NX)
    micstt%CH2P4U=CH2P4(NU(NY,NX),NY,NX)
    micstt%CH2P4BU=CH2P4B(NU(NY,NX),NY,NX)
    micstt%CH1P4U=CH1P4(NU(NY,NX),NY,NX)
    micstt%CH1P4BU=CH1P4B(NU(NY,NX),NY,NX)
    micstt%CNO3SU=CNO3S(NU(NY,NX),NY,NX)
    micstt%CNO3BU=CNO3B(NU(NY,NX),NY,NX)
    micstt%OSC13U=OSC(1,k_POM,NU(NY,NX),NY,NX)
    micstt%OSN13U=OSN(1,k_POM,NU(NY,NX),NY,NX)
    micstt%OSP13U=OSP(1,k_POM,NU(NY,NX),NY,NX)
    micstt%OSC14U=OSC(1,k_humus,NU(NY,NX),NY,NX)
    micstt%OSN14U=OSN(1,k_humus,NU(NY,NX),NY,NX)
    micstt%OSP14U=OSP(1,k_humus,NU(NY,NX),NY,NX)
    micstt%OSC24U=OSC(2,k_humus,NU(NY,NX),NY,NX)
    micstt%OSN24U=OSN(2,k_humus,NU(NY,NX),NY,NX)
    micstt%OSP24U=OSP(2,k_humus,NU(NY,NX),NY,NX)
    micfor%RNH4YU =RNH4Y(NU(NY,NX),NY,NX)
    micfor%RNO3YU =RNO3Y(NU(NY,NX),NY,NX)
    micfor%RPO4YU =RPO4Y(NU(NY,NX),NY,NX)
    micfor%RP14YU =RP14Y(NU(NY,NX),NY,NX)
    micfor%VOLWU =VOLW(NU(NY,NX),NY,NX)
    micfor%CFOMCU=CFOMC(1:2,NU(NY,NX),NY,NX)
  else
    micfor%CFOMC =CFOMC(1:2,L,NY,NX)
  endif
  micstt%CNH4B =CNH4B(L,NY,NX)
  micstt%CNH4S =CNH4S(L,NY,NX)
  micstt%CH2P4 =CH2P4(L,NY,NX)
  micstt%CH2P4B=CH2P4B(L,NY,NX)
  micstt%CH1P4=CH1P4(L,NY,NX)
  micstt%CH1P4B=CH1P4B(L,NY,NX)
  micstt%CNO3S=CNO3S(L,NY,NX)
  micstt%CNO3B=CNO3B(L,NY,NX)
  micfor%RNH4Y =RNH4Y(L,NY,NX)
  micfor%RNO3Y =RNO3Y(L,NY,NX)
  micfor%RPO4Y =RPO4Y(L,NY,NX)
  micfor%RP14Y =RP14Y(L,NY,NX)
  micfor%VOLW  =VOLW(L,NY,NX)

  if(micfor%Lsurf)then
    micfor%BKVL0=BKVL(0,NY,NX)
  endif
  micfor%DFGS(1:NPH)=DFGS(1:NPH,L,NY,NX)
  micfor%FILM(1:NPH)=FILM(1:NPH,L,NY,NX)
  micfor%THETPM(1:NPH)=THETPM(1:NPH,L,NY,NX)
  micfor%VOLWM(1:NPH)=VOLWM(1:NPH,L,NY,NX)
  micfor%TORT(1:NPH)=TORT(1:NPH,L,NY,NX)
  micfor%VOLPM(1:NPH)=VOLPM(1:NPH,L,NY,NX)
  micfor%VOLP=VOLP(L,NY,NX)
  micstt%EPOC=EPOC(L,NY,NX)
  micstt%EHUM=EHUM(L,NY,NX)
  micstt%ZNH4B=ZNH4B(L,NY,NX)
  micstt%ZNH4S=ZNH4S(L,NY,NX)
  micstt%ZNO3B=ZNO3B(L,NY,NX)
  micstt%ZNO3S=ZNO3S(L,NY,NX)
  micstt%H1POB=H1POB(L,NY,NX)
  micstt%H1PO4=H1PO4(L,NY,NX)
  micstt%ZNO2B=ZNO2B(L,NY,NX)
  micstt%ZNO2S=ZNO2S(L,NY,NX)
  micstt%H2POB=H2POB(L,NY,NX)
  micstt%H2PO4=H2PO4(L,NY,NX)
  micstt%CCO2S=CCO2S(L,NY,NX)
  micstt%CNO2S=CNO2S(L,NY,NX)
  micstt%CNO2B=CNO2B(L,NY,NX)
  micstt%CZ2OS=CZ2OS(L,NY,NX)
  micstt%Z2OS=Z2OS(L,NY,NX)
  micstt%COXYS=COXYS(L,NY,NX)
  micstt%OXYS=OXYS(L,NY,NX)
  micstt%SOXYL=SOXYL(L,NY,NX)
  micstt%COXYG=COXYG(L,NY,NX)
  micstt%CZ2GS=CZ2GS(L,NY,NX)
  micstt%CH2GS=CH2GS(L,NY,NX)
  micstt%H2GS=H2GS(L,NY,NX)
  micstt%CCH4G=CCH4G(L,NY,NX)
  micstt%CH4S=CH4S(L,NY,NX)
  micstt%SCH4L=SCH4L(L,NY,NX)
  micstt%ZNFN0=ZNFN0(L,NY,NX)
  micstt%ZNFNI=ZNFNI(L,NY,NX)
  micstt%FOSRH(0:jcplx1)=FOSRH(0:jcplx1,L,NY,NX)
  micstt%OQC(0:jcplx1)=OQC(0:jcplx1,L,NY,NX)
  micstt%OQN(0:jcplx1)=OQN(0:jcplx1,L,NY,NX)
  micstt%OQP(0:jcplx1)=OQP(0:jcplx1,L,NY,NX)
  micstt%OQA(0:jcplx1)=OQA(0:jcplx1,L,NY,NX)
  micstt%OHC(0:jcplx1)=OHC(0:jcplx1,L,NY,NX)
  micstt%OHN(0:jcplx1)=OHN(0:jcplx1,L,NY,NX)
  micstt%OHP(0:jcplx1)=OHP(0:jcplx1,L,NY,NX)
  micstt%OHA(0:jcplx1)=OHA(0:jcplx1,L,NY,NX)

  micstt%OSC(1:jsken,0:jcplx1)=OSC(1:jsken,0:jcplx1,L,NY,NX)
  micstt%OSA(1:jsken,0:jcplx1)=OSA(1:jsken,0:jcplx1,L,NY,NX)
  micstt%OSN(1:jsken,0:jcplx1)=OSN(1:jsken,0:jcplx1,L,NY,NX)
  micstt%OSP(1:jsken,0:jcplx1)=OSP(1:jsken,0:jcplx1,L,NY,NX)
  micstt%ORC(1:ndbiomcp,0:jcplx1)=ORC(1:ndbiomcp,0:jcplx1,L,NY,NX)
  micstt%ORN(1:ndbiomcp,0:jcplx1)=ORN(1:ndbiomcp,0:jcplx1,L,NY,NX)
  micstt%ORP(1:ndbiomcp,0:jcplx1)=ORP(1:ndbiomcp,0:jcplx1,L,NY,NX)
  micstt%CNOSC(1:jsken,0:jcplx1)=CNOSC(1:jsken,0:jcplx1,L,NY,NX)
  micstt%CPOSC(1:jsken,0:jcplx1)=CPOSC(1:jsken,0:jcplx1,L,NY,NX)
  micstt%OMC(1:nlbiomcp,1:NMICBSO,0:jcplx1)=OMC(1:nlbiomcp,1:NMICBSO,0:jcplx1,L,NY,NX)
  micstt%OMN(1:nlbiomcp,1:NMICBSO,0:jcplx1)=OMN(1:nlbiomcp,1:NMICBSO,0:jcplx1,L,NY,NX)
  micstt%OMP(1:nlbiomcp,1:NMICBSO,0:jcplx1)=OMP(1:nlbiomcp,1:NMICBSO,0:jcplx1,L,NY,NX)
  micstt%OMCff(1:nlbiomcp,1:NMICBSA)=OMCff(1:nlbiomcp,1:NMICBSA,L,NY,NX)
  micstt%OMNff(1:nlbiomcp,1:NMICBSA)=OMNff(1:nlbiomcp,1:NMICBSA,L,NY,NX)
  micstt%OMPff(1:nlbiomcp,1:NMICBSA)=OMPff(1:nlbiomcp,1:NMICBSA,L,NY,NX)
  if(.not.micfor%litrm)then
    micfor%AEC=AEC(L,NY,NX)
    micstt%OXYG=OXYG(L,NY,NX)
  endif
  micflx%RVMXC=RVMXC(L,NY,NX)
  micflx%RVMBC=RVMBC(L,NY,NX)
  micflx%RINHOff(1:NMICBSA)=RINHOff(1:NMICBSA,L,NY,NX)
  micflx%RINHBff(1:NMICBSA)=RINHBff(1:NMICBSA,L,NY,NX)
  micflx%RINOOff(1:NMICBSA)=RINOOff(1:NMICBSA,L,NY,NX)
  micflx%RINOBff(1:NMICBSA)=RINOBff(1:NMICBSA,L,NY,NX)
  micflx%RIPOOff(1:NMICBSA)=RIPOOff(1:NMICBSA,L,NY,NX)
  micflx%RIPBOff(1:NMICBSA)=RIPBOff(1:NMICBSA,L,NY,NX)
  micflx%RIPO1ff(1:NMICBSA)=RIPO1ff(1:NMICBSA,L,NY,NX)
  micflx%RIPB1ff(1:NMICBSA)=RIPB1ff(1:NMICBSA,L,NY,NX)
  micflx%ROXYSff(1:NMICBSA)=ROXYSff(1:NMICBSA,L,NY,NX)

  micflx%RINHO(1:NMICBSO,0:JCPLX1)=RINHO(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RINHB(1:NMICBSO,0:JCPLX1)=RINHB(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RINOO(1:NMICBSO,0:JCPLX1)=RINOO(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RINOB(1:NMICBSO,0:JCPLX1)=RINOB(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RIPOO(1:NMICBSO,0:JCPLX1)=RIPOO(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RIPBO(1:NMICBSO,0:JCPLX1)=RIPBO(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RIPO1(1:NMICBSO,0:JCPLX1)=RIPO1(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%RIPB1(1:NMICBSO,0:JCPLX1)=RIPB1(1:NMICBSO,0:JCPLX1,L,NY,NX)
  micflx%ROXYS(1:NMICBSO,0:JCPLX1)=ROXYS(1:NMICBSO,0:JCPLX1,L,NY,NX)

  end subroutine MicAPISend

!------------------------------------------------------------------------------------------


  subroutine MicAPIRecv(L,NY,NX,litrM,micstt,micflx)
  implicit none
  integer, intent(in) :: L,NY,NX
  logical, intent(in) :: litrM
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx
  integer :: NFGs, jcplx1,JG, NMICBSA
  NMICBSA = micpar%NMICBSA
  NFGs=micpar%NFGs
  jcplx1=micpar%jcplx1
  JG=micpar%jguilds
  RCO2O(L,NY,NX) =micflx%RCO2O
  RCH4O(L,NY,NX) =micflx%RCH4O
  RH2GO(L,NY,NX) =micflx%RH2GO
  RUPOXO(L,NY,NX)=micflx%RUPOXO
  RN2G(L,NY,NX)  =micflx%RN2G
  RN2O(L,NY,NX)  =micflx%RN2O
  XNH4S(L,NY,NX) =micflx%XNH4S
  XNO3S(L,NY,NX) =micflx%XNO3S
  XNO2S(L,NY,NX) =micflx%XNO2S
  XH2PS(L,NY,NX) =micflx%XH2PS
  XH1PS(L,NY,NX) =micflx%XH1PS
  XNH4B(L,NY,NX) =micflx%XNH4B
  XNO3B(L,NY,NX) =micflx%XNO3B
  XNO2B(L,NY,NX) =micflx%XNO2B
  XH2BS(L,NY,NX) =micflx%XH2BS
  XH1BS(L,NY,NX) =micflx%XH1BS
  XN2GS(L,NY,NX) =micflx%XN2GS
  RVMXC(L,NY,NX)=micflx%RVMXC
  RVMBC(L,NY,NX)=micflx%RVMBC
  TRINH4(NY,NX)=TRINH4(NY,NX)+micflx%TRINH4
  TRIPO4(NY,NX)=TRIPO4(NY,NX)+micflx%TRIPO4
  XOQCS(0:jcplx1,L,NY,NX)=micflx%XOQCS(0:jcplx1)
  XOQNS(0:jcplx1,L,NY,NX)=micflx%XOQNS(0:jcplx1)
  XOQPS(0:jcplx1,L,NY,NX)=micflx%XOQPS(0:jcplx1)
  XOQAS(0:jcplx1,L,NY,NX)=micflx%XOQAS(0:jcplx1)

  ROXYSff(1:NMICBSA,L,NY,NX)=micflx%ROXYSff(1:NMICBSA)
  RVMX4ff(1:NMICBSA,L,NY,NX)=micflx%RVMX4ff(1:NMICBSA)
  RVMB4ff(1:NMICBSA,L,NY,NX)=micflx%RVMB4ff(1:NMICBSA)
  RVMX2ff(1:NMICBSA,L,NY,NX)=micflx%RVMX2ff(1:NMICBSA)
  RVMB2ff(1:NMICBSA,L,NY,NX)=micflx%RVMB2ff(1:NMICBSA)
  ROXYS(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%ROXYS(1:NMICBSO,0:JCPLX1)
  ROQCS(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%ROQCS(1:NMICBSO,0:JCPLX1)
  ROQAS(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%ROQAS(1:NMICBSO,0:JCPLX1)
  RVMX3(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMX3(1:NMICBSO,0:JCPLX1)
  RVMB3(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMB3(1:NMICBSO,0:JCPLX1)
  RVMX2(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMX2(1:NMICBSO,0:JCPLX1)
  RVMB2(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMB2(1:NMICBSO,0:JCPLX1)
  RVMX1(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMX1(1:NMICBSO,0:JCPLX1)
  RVMX4(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMX4(1:NMICBSO,0:JCPLX1)
  RVMB4(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RVMB4(1:NMICBSO,0:JCPLX1)
  RINHO(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RINHO(1:NMICBSO,0:JCPLX1)
  RINHB(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RINHB(1:NMICBSO,0:JCPLX1)
  RINOO(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RINOO(1:NMICBSO,0:JCPLX1)
  RINOB(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RINOB(1:NMICBSO,0:JCPLX1)
  RIPOO(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RIPOO(1:NMICBSO,0:JCPLX1)
  RIPBO(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RIPBO(1:NMICBSO,0:JCPLX1)
  RIPO1(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RIPO1(1:NMICBSO,0:JCPLX1)
  RIPB1(1:NMICBSO,0:JCPLX1,L,NY,NX)=micflx%RIPB1(1:NMICBSO,0:JCPLX1)

  OSC(1:jsken,0:jcplx1,L,NY,NX)=micstt%OSC(1:jsken,0:jcplx1)
  OSA(1:jsken,0:jcplx1,L,NY,NX)=micstt%OSA(1:jsken,0:jcplx1)
  OSN(1:jsken,0:jcplx1,L,NY,NX)=micstt%OSN(1:jsken,0:jcplx1)
  OSP(1:jsken,0:jcplx1,L,NY,NX)=micstt%OSP(1:jsken,0:jcplx1)

  if(litrm)then
    RINHOR(1:NMICBSO,0:JCPLX1,NY,NX)=micflx%RINHOR(1:NMICBSO,0:JCPLX1)
    RINOOR(1:NMICBSO,0:JCPLX1,NY,NX)=micflx%RINOOR(1:NMICBSO,0:JCPLX1)
    RIPOOR(1:NMICBSO,0:JCPLX1,NY,NX)=micflx%RIPOOR(1:NMICBSO,0:JCPLX1)
    RIPO1R(1:NMICBSO,0:JCPLX1,NY,NX)=micflx%RIPO1R(1:NMICBSO,0:JCPLX1)
    RINHORff(1:NMICBSA,NY,NX)=micflx%RINHORff(1:NMICBSA)
    RINOORff(1:NMICBSA,NY,NX)=micflx%RINOORff(1:NMICBSA)
    RIPOORff(1:NMICBSA,NY,NX)=micflx%RIPOORff(1:NMICBSA)
    RIPO1Rff(1:NMICBSA,NY,NX)=micflx%RIPO1Rff(1:NMICBSA)
    OSC(1,3,NU(NY,NX),NY,NX)=micstt%OSC13U
    OSC(1,4,NU(NY,NX),NY,NX)=micstt%OSC14U
    OSC(2,4,NU(NY,NX),NY,NX)=micstt%OSC24U
    OSN(1,3,NU(NY,NX),NY,NX)=micstt%OSN13U
    OSN(1,4,NU(NY,NX),NY,NX)=micstt%OSN14U
    OSN(2,4,NU(NY,NX),NY,NX)=micstt%OSN24U
    OSP(1,3,NU(NY,NX),NY,NX)=micstt%OSP13U
    OSP(1,4,NU(NY,NX),NY,NX)=micstt%OSP14U
    OSP(2,4,NU(NY,NX),NY,NX)=micstt%OSP24U
  endif

  RINHOff(1:NMICBSA,L,NY,NX)=micflx%RINHOff(1:NMICBSA)
  RINHBff(1:NMICBSA,L,NY,NX)=micflx%RINHBff(1:NMICBSA)
  RINOOff(1:NMICBSA,L,NY,NX)=micflx%RINOOff(1:NMICBSA)
  RINOBff(1:NMICBSA,L,NY,NX)=micflx%RINOBff(1:NMICBSA)
  RIPOOff(1:NMICBSA,L,NY,NX)=micflx%RIPOOff(1:NMICBSA)
  RIPBOff(1:NMICBSA,L,NY,NX)=micflx%RIPBOff(1:NMICBSA)
  RIPO1ff(1:NMICBSA,L,NY,NX)=micflx%RIPO1ff(1:NMICBSA)
  RIPB1ff(1:NMICBSA,L,NY,NX)=micflx%RIPB1ff(1:NMICBSA)
  ROXSK(1:NPH,L,NY,NX)=micflx%ROXSK(1:NPH)

  TFNQ(L,NY,NX)=micstt%TFNQ
  VOLQ(L,NY,NX)=micstt%VOLQ
  TOQCK(L,NY,NX)=micstt%TOQCK
  ZNFNI(L,NY,NX)=micstt%ZNFNI
  FOSRH(0:jcplx1,L,NY,NX)=micstt%FOSRH(0:jcplx1)
  OQC(0:jcplx1,L,NY,NX)=micstt%OQC(0:jcplx1)
  OQN(0:jcplx1,L,NY,NX)=micstt%OQN(0:jcplx1)
  OQP(0:jcplx1,L,NY,NX)=micstt%OQP(0:jcplx1)
  OQA(0:jcplx1,L,NY,NX)=micstt%OQA(0:jcplx1)
  OHC(0:jcplx1,L,NY,NX)=micstt%OHC(0:jcplx1)
  OHN(0:jcplx1,L,NY,NX)=micstt%OHN(0:jcplx1)
  OHP(0:jcplx1,L,NY,NX)=micstt%OHP(0:jcplx1)
  OHA(0:jcplx1,L,NY,NX)=micstt%OHA(0:jcplx1)

  ORC(1:ndbiomcp,0:jcplx1,L,NY,NX)=micstt%ORC(1:ndbiomcp,0:jcplx1)
  ORN(1:ndbiomcp,0:jcplx1,L,NY,NX)=micstt%ORN(1:ndbiomcp,0:jcplx1)
  ORP(1:ndbiomcp,0:jcplx1,L,NY,NX)=micstt%ORP(1:ndbiomcp,0:jcplx1)
  OMC(1:nlbiomcp,1:NMICBSO,0:JCPLX1,L,NY,NX)=micstt%OMC(1:nlbiomcp,1:NMICBSO,0:JCPLX1)
  OMN(1:nlbiomcp,1:NMICBSO,0:JCPLX1,L,NY,NX)=micstt%OMN(1:nlbiomcp,1:NMICBSO,0:JCPLX1)
  OMP(1:nlbiomcp,1:NMICBSO,0:JCPLX1,L,NY,NX)=micstt%OMP(1:nlbiomcp,1:NMICBSO,0:JCPLX1)
  OMCff(1:nlbiomcp,1:NMICBSA,L,NY,NX)=micstt%OMCff(1:nlbiomcp,1:NMICBSA)
  OMNff(1:nlbiomcp,1:NMICBSA,L,NY,NX)=micstt%OMNff(1:nlbiomcp,1:NMICBSA)
  OMPff(1:nlbiomcp,1:NMICBSA,L,NY,NX)=micstt%OMPff(1:nlbiomcp,1:NMICBSA)

  end subroutine MicAPIRecv
end module MicBGCAPI
