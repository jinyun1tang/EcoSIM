module MicBGCAPI

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use NitrosMod      , only : VerticalLitterMixLvsLL
  use NitroDisturbMod, only : SOMRemovalByDisturbance
  use EcoSIMSolverPar
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod, only : micforctype
  use minimathmod, only : AZMAX1
  use TracerIDMod
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
  use EcoSiMParDataMod, only : micpar
  use MicBGCMod, only : SoilBGCOneLayer
implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

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
        IF(VLSoilPoreMicP(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
             call MicBGC1Layer(I,J,L,NY,NX)
          ELSE
            RCO2O(L,NY,NX)=0.0_r8
            RCH4O(L,NY,NX)=0.0_r8
            RH2GO(L,NY,NX)=0.0_r8
            RUPOXO(L,NY,NX)=0.0_r8
            RN2G(L,NY,NX)=0.0_r8
            RN2O(L,NY,NX)=0.0_r8
            RNH4MicbTransf_vr(L,NY,NX)=0.0_r8
            RNO3MicbTransf_vr(L,NY,NX)=0.0_r8
            RNO2MicbTransf_vr(L,NY,NX)=0.0_r8
            RH2PO4MicbTransf_vr(L,NY,NX)=0.0_r8
            RH1PO4MicbTransf_vr(L,NY,NX)=0.0_r8
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
          RNH4MicbTransf_vr(L,NY,NX)=0.0_r8
          RNO3MicbTransf_vr(L,NY,NX)=0.0_r8
          RNO2MicbTransf_vr(L,NY,NX)=0.0_r8
          RH2PO4MicbTransf_vr(L,NY,NX)=0.0_r8
          RH1PO4MicbTransf_vr(L,NY,NX)=0.0_r8
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

  integer :: NFGs, jcplx, k_POM, k_humus
  integer :: kk, ndbiomcp, nlbiomcp, NumOfMicrobsInAutotrophCmplx, NumOfMicrobs1HetertrophCmplx
  NFGs=micpar%NFGs
  jcplx=micpar%jcplx

  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp
  NumOfMicrobsInAutotrophCmplx  = micpar%NumOfMicrobsInAutotrophCmplx
  NumOfMicrobs1HetertrophCmplx  = micpar%NumOfMicrobs1HetertrophCmplx
  k_humus  = micpar%k_humus
  k_POM    = micpar%k_POM

  micfor%ZERO  =ZERO
  micfor%CCH4E =AtmGgms(idg_CH4,NY,NX)
  micfor%COXYE =AtmGgms(idg_O2,NY,NX)
  micfor%O2_irrig_conc  =O2_irrig_conc(NY,NX)
  micfor%O2_rain_conc  =O2_rain_conc(NY,NX)
  micfor%Irrig2LitRSurf =Irrig2LitRSurf(NY,NX)
  micfor%Rain2LitRSurf =Rain2LitRSurf(NY,NX)
  micfor%OFFSET=OFFSET(NY,NX)
  micfor%VLitR  =VLitR(NY,NX)
  micfor%VWatLitRHoldCapcity=VWatLitRHoldCapcity(NY,NX)
  micfor%ZEROS2=ZEROS2(NY,NX)
  micfor%ZEROS =ZEROS(NY,NX)
  micfor%VLSoilMicP  =VLSoilMicP(L,NY,NX)
  micfor%THETY =THETY(L,NY,NX)
  micfor%POROS =POROS(L,NY,NX)
  micfor%FieldCapacity    =FieldCapacity(L,NY,NX)
  micfor%TKS   =TKS(L,NY,NX)
  micfor%THETW =THETW(L,NY,NX)
  micfor%PH    =PH(L,NY,NX)
  micfor%SoilMicPMassLayer  =SoilMicPMassLayer(L,NY,NX)
  micfor%VLSoilPoreMicP  =VLSoilPoreMicP(L,NY,NX)
  micfor%TFND  =TFND(L,NY,NX)
  micfor%VLNOB =trcs_VLN(ids_NO3B,L,NY,NX)
  micfor%VLNO3 =trcs_VLN(ids_NO3,L,NY,NX)
  micfor%VLNH4 =trcs_VLN(ids_NH4,L,NY,NX)
  micfor%VLNHB =trcs_VLN(ids_NH4B,L,NY,NX)
  micfor%VLPO4 =trcs_VLN(ids_H1PO4,L,NY,NX)
  micfor%VLPOB =trcs_VLN(ids_H1PO4B,L,NY,NX)
  micfor%PSISoilMatricP =PSISoilMatricP(L,NY,NX)
  micfor%OLSGL =SolDifc(idg_O2,L,NY,NX)
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
  micfor%ROQCY(1:jcplx)=ROQCY(1:jcplx,L,NY,NX)
  micfor%ROQAY(1:jcplx)=ROQAY(1:jcplx,L,NY,NX)
  micfor%litrm=(L==0)
  micfor%Lsurf=(L==NU(NY,NX))
  if(micfor%litrm)then
    micstt%ZNH4TU=AZMAX1(trc_solml(ids_NH4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml(ids_NH4B,NU(NY,NX),NY,NX))
    micstt%ZNO3TU=AZMAX1(trc_solml(ids_NO3,NU(NY,NX),NY,NX))+AZMAX1(trc_solml(ids_NO3B,NU(NY,NX),NY,NX))
    micstt%H1P4TU=AZMAX1(trc_solml(ids_H1PO4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml(ids_H1PO4B,NU(NY,NX),NY,NX))
    micstt%H2P4TU=AZMAX1(trc_solml(ids_H2PO4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml(ids_H2PO4B,NU(NY,NX),NY,NX))
    micstt%CNH4BU=trc_solcl(ids_NH4B,NU(NY,NX),NY,NX)
    micstt%CNH4SU=trc_solcl(ids_NH4,NU(NY,NX),NY,NX)
    micstt%CH2P4U=trc_solcl(ids_H2PO4,NU(NY,NX),NY,NX)
    micstt%CH2P4BU=trc_solcl(ids_H2PO4B,NU(NY,NX),NY,NX)
    micstt%CH1P4U=trc_solcl(ids_H1PO4,NU(NY,NX),NY,NX)
    micstt%CH1P4BU=trc_solcl(ids_H1PO4B,NU(NY,NX),NY,NX)
    micstt%CNO3SU=trc_solcl(ids_NO3,NU(NY,NX),NY,NX)
    micstt%CNO3BU=trc_solcl(ids_NO3B,NU(NY,NX),NY,NX)
    micstt%OSC13U=OSC(micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
    micstt%OSN13U=OSN(micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
    micstt%OSP13U=OSP(micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
    micstt%OSC14U=OSC(micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSN14U=OSN(micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSP14U=OSP(micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSC24U=OSC(micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSN24U=OSN(micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSP24U=OSP(micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    micfor%RNH4YU =RNH4Y(NU(NY,NX),NY,NX)
    micfor%RNO3YU =RNO3Y(NU(NY,NX),NY,NX)
    micfor%RPO4YU =RPO4Y(NU(NY,NX),NY,NX)
    micfor%RP14YU =RP14Y(NU(NY,NX),NY,NX)
    micfor%VOLWU =VLWatMicP(NU(NY,NX),NY,NX)
    micfor%CFOMCU=CFOMC(1:2,NU(NY,NX),NY,NX)
  else
    micfor%CFOMC =CFOMC(1:2,L,NY,NX)
  endif
  micstt%CNH4B =trc_solcl(ids_NH4B,L,NY,NX)
  micstt%CNH4S =trc_solcl(ids_NH4,L,NY,NX)
  micstt%CH2P4 =trc_solcl(ids_H2PO4,L,NY,NX)
  micstt%CH2P4B=trc_solcl(ids_H2PO4B,L,NY,NX)
  micstt%CH1P4=trc_solcl(ids_H1PO4,L,NY,NX)
  micstt%CH1P4B=trc_solcl(ids_H1PO4B,L,NY,NX)
  micstt%CNO3S=trc_solcl(ids_NO3,L,NY,NX)
  micstt%CNO3B=trc_solcl(ids_NO3B,L,NY,NX)
  micfor%RNH4Y =RNH4Y(L,NY,NX)
  micfor%RNO3Y =RNO3Y(L,NY,NX)
  micfor%RPO4Y =RPO4Y(L,NY,NX)
  micfor%RP14Y =RP14Y(L,NY,NX)
  micfor%VLWatMicP  =VLWatMicP(L,NY,NX)

  if(micfor%Lsurf)then
    micfor%SoilMicPMassLayer0=SoilMicPMassLayer(0,NY,NX)
  endif
  micfor%DiffusivitySolutEff(1:NPH)=DiffusivitySolutEff(1:NPH,L,NY,NX)
  micfor%FILM(1:NPH)=FILM(1:NPH,L,NY,NX)
  micfor%THETPM(1:NPH)=THETPM(1:NPH,L,NY,NX)
  micfor%VLWatMicPM(1:NPH)=VLWatMicPM(1:NPH,L,NY,NX)
  micfor%TortMicPM(1:NPH)=TortMicPM(1:NPH,L,NY,NX)
  micfor%VLsoiAirPM(1:NPH)=VLsoiAirPM(1:NPH,L,NY,NX)
  micfor%VLsoiAirP=VLsoiAirP(L,NY,NX)
  micstt%EPOC=EPOC(L,NY,NX)
  micstt%EHUM=EHUM(L,NY,NX)
  micstt%ZNH4B=trc_solml(ids_NH4B,L,NY,NX)
  micstt%ZNH4S=trc_solml(ids_NH4,L,NY,NX)
  micstt%ZNO3B=trc_solml(ids_NO3B,L,NY,NX)
  micstt%ZNO3S=trc_solml(ids_NO3,L,NY,NX)
  micstt%H1POB=trc_solml(ids_H1PO4B,L,NY,NX)
  micstt%H1PO4=trc_solml(ids_H1PO4,L,NY,NX)
  micstt%ZNO2B=trc_solml(ids_NO2B,L,NY,NX)
  micstt%ZNO2S=trc_solml(ids_NO2,L,NY,NX)
  micstt%H2POB=trc_solml(ids_H2PO4B,L,NY,NX)
  micstt%H2PO4=trc_solml(ids_H2PO4,L,NY,NX)
  micstt%CCO2S=trc_solcl(idg_CO2,L,NY,NX)
  micstt%CNO2S=trc_solcl(ids_NO2,L,NY,NX)
  micstt%CNO2B=trc_solcl(ids_NO2B,L,NY,NX)
  micstt%CZ2OS=trc_solcl(idg_N2O,L,NY,NX)
  micstt%Z2OS=trc_solml(idg_N2O,L,NY,NX)
  micstt%COXYS=trc_solcl(idg_O2,L,NY,NX)
  micstt%OXYS=trc_solml(idg_O2,L,NY,NX)
  micstt%SOXYL=GasSolbility(idg_O2,L,NY,NX)
  micstt%COXYG=trc_gascl(idg_O2,L,NY,NX)
  micstt%CZ2GS=trc_solcl(idg_N2,L,NY,NX)
  micstt%CH2GS=trc_solcl(idg_H2,L,NY,NX)
  micstt%H2GS=trc_solml(idg_H2,L,NY,NX)
  micstt%CCH4G=trc_gascl(idg_CH4,L,NY,NX)
  micstt%CH4S=trc_solml(idg_CH4,L,NY,NX)
  micstt%SCH4L=GasSolbility(idg_CH4,L,NY,NX)
  micstt%ZNFN0=ZNFN0(L,NY,NX)
  micstt%ZNFNI=ZNFNI(L,NY,NX)
  micstt%FOSRH(1:jcplx)=FOSRH(1:jcplx,L,NY,NX)
  micstt%DOM(idom_doc,1:jcplx)=DOM(idom_doc,1:jcplx,L,NY,NX)
  micstt%DOM(idom_don,1:jcplx)=DOM(idom_don,1:jcplx,L,NY,NX)
  micstt%DOM(idom_dop,1:jcplx)=DOM(idom_dop,1:jcplx,L,NY,NX)
  micstt%DOM(idom_acetate,1:jcplx)=DOM(idom_acetate,1:jcplx,L,NY,NX)
  micstt%OHC(1:jcplx)=OHC(1:jcplx,L,NY,NX)
  micstt%OHN(1:jcplx)=OHN(1:jcplx,L,NY,NX)
  micstt%OHP(1:jcplx)=OHP(1:jcplx,L,NY,NX)
  micstt%OHA(1:jcplx)=OHA(1:jcplx,L,NY,NX)

  micstt%OSC(1:jsken,1:jcplx)=OSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%OSA(1:jsken,1:jcplx)=OSA(1:jsken,1:jcplx,L,NY,NX)
  micstt%OSN(1:jsken,1:jcplx)=OSN(1:jsken,1:jcplx,L,NY,NX)
  micstt%OSP(1:jsken,1:jcplx)=OSP(1:jsken,1:jcplx,L,NY,NX)
  micstt%ORC(1:ndbiomcp,1:jcplx)=ORC(1:ndbiomcp,1:jcplx,L,NY,NX)
  micstt%ORN(1:ndbiomcp,1:jcplx)=ORN(1:ndbiomcp,1:jcplx,L,NY,NX)
  micstt%ORP(1:ndbiomcp,1:jcplx)=ORP(1:ndbiomcp,1:jcplx,L,NY,NX)
  micstt%CNOSC(1:jsken,1:jcplx)=CNOSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%CPOSC(1:jsken,1:jcplx)=CPOSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%OMC(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=OMC(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micstt%OMN(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=OMN(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micstt%OMP(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=OMP(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micstt%OMCff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx)=OMCff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micstt%OMNff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx)=OMNff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micstt%OMPff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx)=OMPff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  if(.not.micfor%litrm)then
    micfor%AEC=AEC(L,NY,NX)
    micstt%OXYG=trc_gasml(idg_O2,L,NY,NX)
  endif
  micflx%RVMXC=RVMXC(L,NY,NX)
  micflx%RVMBC=RVMBC(L,NY,NX)
  micflx%RINHOff(1:NumOfMicrobsInAutotrophCmplx)=RINHOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RINHBff(1:NumOfMicrobsInAutotrophCmplx)=RINHBff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RINOOff(1:NumOfMicrobsInAutotrophCmplx)=RINOOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RINOBff(1:NumOfMicrobsInAutotrophCmplx)=RINOBff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RIPOOff(1:NumOfMicrobsInAutotrophCmplx)=RIPOOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RIPBOff(1:NumOfMicrobsInAutotrophCmplx)=RIPBOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RIPO1ff(1:NumOfMicrobsInAutotrophCmplx)=RIPO1ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%RIPB1ff(1:NumOfMicrobsInAutotrophCmplx)=RIPB1ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)
  micflx%ROXYSff(1:NumOfMicrobsInAutotrophCmplx)=ROXYSff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)

  micflx%RINHO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RINHO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RINHB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RINHB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RINOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RINOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RINOB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RINOB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RIPOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RIPOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RIPBO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RIPBO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RIPO1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RIPO1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%RIPB1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=RIPB1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)
  micflx%ROXYS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)=ROXYS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)

  end subroutine MicAPISend

!------------------------------------------------------------------------------------------


  subroutine MicAPIRecv(L,NY,NX,litrM,micstt,micflx)
  implicit none
  integer, intent(in) :: L,NY,NX
  logical, intent(in) :: litrM
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx
  integer :: NFGs, jcplx, NumOfMicrobsInAutotrophCmplx
  NumOfMicrobsInAutotrophCmplx = micpar%NumOfMicrobsInAutotrophCmplx
  NFGs=micpar%NFGs
  jcplx=micpar%jcplx

  RCO2O(L,NY,NX) =micflx%RCO2O
  RCH4O(L,NY,NX) =micflx%RCH4O
  RH2GO(L,NY,NX) =micflx%RH2GO
  RUPOXO(L,NY,NX)=micflx%RUPOXO
  RN2G(L,NY,NX)  =micflx%RN2G
  RN2O(L,NY,NX)  =micflx%RN2O
  RNH4MicbTransf_vr(L,NY,NX) =micflx%RNH4MicbTransf_vr
  RNO3MicbTransf_vr(L,NY,NX) =micflx%RNO3MicbTransf_vr
  RNO2MicbTransf_vr(L,NY,NX) =micflx%RNO2MicbTransf_vr
  RH2PO4MicbTransf_vr(L,NY,NX) =micflx%RH2PO4MicbTransf_vr
  RH1PO4MicbTransf_vr(L,NY,NX) =micflx%RH1PO4MicbTransf_vr
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
  RDOM_micb_flx(idom_doc,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_doc,1:jcplx)
  RDOM_micb_flx(idom_don,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_don,1:jcplx)
  RDOM_micb_flx(idom_dop,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_dop,1:jcplx)
  RDOM_micb_flx(idom_acetate,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_acetate,1:jcplx)

  ROXYSff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%ROXYSff(1:NumOfMicrobsInAutotrophCmplx)
  RVMX4ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RVMX4ff(1:NumOfMicrobsInAutotrophCmplx)
  RVMB4ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RVMB4ff(1:NumOfMicrobsInAutotrophCmplx)
  RVMX2ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RVMX2ff(1:NumOfMicrobsInAutotrophCmplx)
  RVMB2ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RVMB2ff(1:NumOfMicrobsInAutotrophCmplx)
  ROXYS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%ROXYS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  ROQCS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%ROQCS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  ROQAS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%ROQAS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMX3(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMX3(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMB3(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMB3(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMX2(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMX2(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMB2(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMB2(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMX1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMX1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMX4(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMX4(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RVMB4(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RVMB4(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RINHO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RINHO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RINHB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RINHB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RINOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RINOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RINOB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RINOB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RIPOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RIPOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RIPBO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RIPBO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RIPO1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RIPO1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  RIPB1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micflx%RIPB1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)

  OSC(1:jsken,1:jcplx,L,NY,NX)=micstt%OSC(1:jsken,1:jcplx)
  OSA(1:jsken,1:jcplx,L,NY,NX)=micstt%OSA(1:jsken,1:jcplx)
  OSN(1:jsken,1:jcplx,L,NY,NX)=micstt%OSN(1:jsken,1:jcplx)
  OSP(1:jsken,1:jcplx,L,NY,NX)=micstt%OSP(1:jsken,1:jcplx)

  if(litrm)then
    RINHOR(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,NY,NX)=micflx%RINHOR(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
    RINOOR(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,NY,NX)=micflx%RINOOR(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
    RIPOOR(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,NY,NX)=micflx%RIPOOR(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
    RIPO1R(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,NY,NX)=micflx%RIPO1R(1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
    RINHORff(1:NumOfMicrobsInAutotrophCmplx,NY,NX)=micflx%RINHORff(1:NumOfMicrobsInAutotrophCmplx)
    RINOORff(1:NumOfMicrobsInAutotrophCmplx,NY,NX)=micflx%RINOORff(1:NumOfMicrobsInAutotrophCmplx)
    RIPOORff(1:NumOfMicrobsInAutotrophCmplx,NY,NX)=micflx%RIPOORff(1:NumOfMicrobsInAutotrophCmplx)
    RIPO1Rff(1:NumOfMicrobsInAutotrophCmplx,NY,NX)=micflx%RIPO1Rff(1:NumOfMicrobsInAutotrophCmplx)
    OSC(micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)=micstt%OSC13U
    OSC(micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSC14U
    OSC(micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSC24U

    OSN(micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)=micstt%OSN13U
    OSN(micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSN14U
    OSN(micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSN24U

    OSP(micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)=micstt%OSP13U
    OSP(micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSP14U
    OSP(micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSP24U
  endif

  RINHOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RINHOff(1:NumOfMicrobsInAutotrophCmplx)
  RINHBff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RINHBff(1:NumOfMicrobsInAutotrophCmplx)
  RINOOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RINOOff(1:NumOfMicrobsInAutotrophCmplx)
  RINOBff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RINOBff(1:NumOfMicrobsInAutotrophCmplx)
  RIPOOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RIPOOff(1:NumOfMicrobsInAutotrophCmplx)
  RIPBOff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RIPBOff(1:NumOfMicrobsInAutotrophCmplx)
  RIPO1ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RIPO1ff(1:NumOfMicrobsInAutotrophCmplx)
  RIPB1ff(1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micflx%RIPB1ff(1:NumOfMicrobsInAutotrophCmplx)
  ROXSK(1:NPH,L,NY,NX)=micflx%ROXSK(1:NPH)

  TFNQ(L,NY,NX)=micstt%TFNQ
  VOLQ(L,NY,NX)=micstt%VOLQ
  TOQCK(L,NY,NX)=micstt%TOQCK
  ZNFNI(L,NY,NX)=micstt%ZNFNI
  FOSRH(1:jcplx,L,NY,NX)=micstt%FOSRH(1:jcplx)
  DOM(idom_doc,1:jcplx,L,NY,NX)=micstt%DOM(idom_doc,1:jcplx)
  DOM(idom_don,1:jcplx,L,NY,NX)=micstt%DOM(idom_don,1:jcplx)
  DOM(idom_dop,1:jcplx,L,NY,NX)=micstt%DOM(idom_dop,1:jcplx)
  DOM(idom_acetate,1:jcplx,L,NY,NX)=micstt%DOM(idom_acetate,1:jcplx)
  OHC(1:jcplx,L,NY,NX)=micstt%OHC(1:jcplx)
  OHN(1:jcplx,L,NY,NX)=micstt%OHN(1:jcplx)
  OHP(1:jcplx,L,NY,NX)=micstt%OHP(1:jcplx)
  OHA(1:jcplx,L,NY,NX)=micstt%OHA(1:jcplx)

  ORC(1:ndbiomcp,1:jcplx,L,NY,NX)=micstt%ORC(1:ndbiomcp,1:jcplx)
  ORN(1:ndbiomcp,1:jcplx,L,NY,NX)=micstt%ORN(1:ndbiomcp,1:jcplx)
  ORP(1:ndbiomcp,1:jcplx,L,NY,NX)=micstt%ORP(1:ndbiomcp,1:jcplx)
  OMC(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micstt%OMC(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  OMN(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micstt%OMN(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  OMP(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)=micstt%OMP(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,1:jcplx)
  OMCff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micstt%OMCff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx)
  OMNff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micstt%OMNff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx)
  OMPff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX)=micstt%OMPff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx)

  end subroutine MicAPIRecv
end module MicBGCAPI
