module TrnsfrMod
!!
! Description:
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
  use SOMDataType
  use ChemTranspDataType
  use GridConsts
  use EcoSIMSolverPar
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use EcosimConst
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use LandSurfDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use BoundaryTranspMod
  use InsideTranspMod
  use TransfrDataMod
  use IrrigationDataType
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8) :: XN

  real(r8) ::  CHY0(0:JZ,JY,JX),RHGFLZ(JZ,JY,JX)

  real(r8) :: RCOFLZ(JZ,JY,JX),RCHFLZ(JZ,JY,JX)
  real(r8) :: ROXFLZ(JZ,JY,JX),RNGFLZ(JZ,JY,JX)
  real(r8) :: RN2FLZ(JZ,JY,JX),RN4FLZ(JZ,JY,JX),RN3FLZ(JZ,JY,JX)
  real(r8) :: RNOFLZ(JZ,JY,JX),RH2PFZ(JZ,JY,JX),RN4FBZ(JZ,JY,JX)
  real(r8) :: RN3FBZ(JZ,JY,JX),RNOFBZ(JZ,JY,JX),RH2BBZ(JZ,JY,JX)
  real(r8) :: RH1PFZ(JZ,JY,JX),RH1BBZ(JZ,JY,JX)

  real(r8), PARAMETER :: DPN4=5.7E-07
  REAL(r8) :: CCO2SQ,CCH4SQ,COXYSQ,CZ2GSQ,CZ2OSQ,CNH3SQ &
  ,CNH3BQ,CH2GSQ

  public :: trnsfr
  public :: InitTrnsfr,DestructTrnsfr
  contains

  subroutine InitTrnsfr()
  implicit none

  allocate(RQROC0(0:jcplx1,JV,JH))
  allocate(RQRON0(0:jcplx1,JV,JH))
  allocate(RQROA0(0:jcplx1,JV,JH))
  allocate(RQROP0(0:jcplx1,JV,JH))
  allocate(OQC2(0:jcplx1,0:JZ,JY,JX))
  allocate(OQN2(0:jcplx1,0:JZ,JY,JX))
  allocate(OQP2(0:jcplx1,0:JZ,JY,JX))
  allocate(OQA2(0:jcplx1,0:JZ,JY,JX))
  allocate(ROCSK2(0:jcplx1,0:JZ,JY,JX))
  allocate(RONSK2(0:jcplx1,0:JZ,JY,JX))
  allocate(ROPSK2(0:jcplx1,0:JZ,JY,JX))
  allocate(ROASK2(0:jcplx1,0:JZ,JY,JX))
  allocate(ROCFLS(0:jcplx1,3,0:JD,JV,JH))
  allocate(RONFLS(0:jcplx1,3,0:JD,JV,JH))
  allocate(ROPFLS(0:jcplx1,3,0:JD,JV,JH))
  allocate(ROAFLS(0:jcplx1,3,0:JD,JV,JH))
  allocate(ROCFHS(0:jcplx1,3,JD,JV,JH))
  allocate(RONFHS(0:jcplx1,3,JD,JV,JH))
  allocate(ROPFHS(0:jcplx1,3,JD,JV,JH))
  allocate(ROAFHS(0:jcplx1,3,JD,JV,JH))
  allocate(TOCFLS(0:jcplx1,JZ,JY,JX))
  allocate(TONFLS(0:jcplx1,JZ,JY,JX))
  allocate(TOPFLS(0:jcplx1,JZ,JY,JX))
  allocate(TOAFLS(0:jcplx1,JZ,JY,JX))

  call InitInsTp
  end subroutine InitTrnsfr
!------------------------------------------------------------------------------------------
  subroutine DestructTrnsfr
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


  call DestructInsTp
  end subroutine DestructTrnsfr

!------------------------------------------------------------------------------------------

  SUBROUTINE trnsfr(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     NON-SALT SOLUTES AND GASES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: MX,MM,M

!     execution begins here
!
!     TIME STEPS FOR SOLUTE AND GAS FLUX CALCULATIONS
!

  call InitFluxandStateVariables(I,NHW,NHE,NVN,NVS)


!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=number of cycles per hour for gas flux calculations
! XNPT=1/number of cycles NPH-1 for gas flux calculations
  MX=0
  DO  MM=1,NPG
    M=MIN(NPH,INT((MM-1)*XNPT)+1)

    call ModelSoluteHydroFlux(M,MX,NHW, NHE, NVN, NVS)
!
!     BOUNDARY SOLUTE AND GAS FLUXES
!
    call BoundaryFlux(M,MX,NHW,NHE,NVN,NVS)
!
!     UPDATE STATE VARIABLES FROM TOTAL FLUXES CALCULATED ABOVE
!
    call UpdateStateVar(M,MM,MX,NPG,NHW,NHE,NVN,NVS)
    MX=M
  ENDDO
  RETURN
  END subroutine trnsfr
!------------------------------------------------------------------------------------------

  subroutine InitFluxandStateVariables(I,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: I, NHW,NHE,NVN,NVS

  integer :: NX,NY
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

!
!     GAS AND SOLUTE SINKS AND SOURCES IN SURFACE RESIDUE FROM MICROBIAL
!     TRANSFORMATIONS IN 'NITRO' + ROOT EXCHANGE IN 'EXTRACT'
!     + EQUILIBRIA REACTIONS IN 'SOLUTE' AT SUB-HOURLY TIME STEP
!
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
!     XNPG=1/number of cycles h-1 for gas flux calculations
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     RCO2O=net soil CO2 uptake from nitro.f
!     RCH4O=net soil CH4 uptake from nitro.f
!     RN2G=total soil N2 production from nitro.f
!     XN2GS=total N2 fixation from nitro.f
!     RN2O=net soil N2O uptake from nitro.f
!     RH2GO=net H2 uptake from nitro.f
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate from nitro.f
!     XNH4S=net change in NH4 from nitro.f
!     TRN4S,TRN3S=NH4,NH3 dissolution from solute.f
!     XNO3S=net change in NO3 from nitro.f
!     TRNO3=NO3 dissolution from solute.f
!     XNO2S=net change in NO2 from nitro.f
!     TRNO2=NO2 dissolution from solute.f
!     XH2PS=net change in H2PO4 from nitro.f
!     TRH2P=H2PO4 dissolution from solute.f
!     XH1PS=net change in HPO4 from nitro.f
!     TRH1P=HPO4 dissolution from solute.f
!
      call SurfaceSinksandSources(NY,NX)
!
!     INITIALIZE STATE VARIABLES FOR USE IN GAS, SOLUTE FLUX CALCULATIONS
!
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 content
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate from nitro.f
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4
!     CHY0=H concentration
!     PH=pH
!
      call StateVarforGasandSolute(NY,NX)
!
!     INITIALIZE SURFACE SOLUTE FLUXES FROM ATMOSPHERE
!
!     X*FLS,X*FHS=hourly solute flux in macropores,micropores
!
      call SurfaceSolutefromAtms(NY,NX)
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
!     PRECW,PRECR=snow,rain
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     X*BLS=hourly solute flux to snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
      call HourlySoluteFluxes(I,NY,NX)
!     GAS AND SOLUTE FLUXES AT SUB-HOURLY FLUX TIME STEP
!     ENTERED IN SITE FILE
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     R*BLS,R*FL0,R*FL1,R*FL2=solute flux to snowpack,surface litter,soil surface non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
      call SubHourlyFluxesFromSiteFile(NY,NX)
!
!     SOLUTE FLUXES FROM WATSUB.F, NITRO.F, UPTAKE.F, SOLUTE.F
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     CHY0=H concentration
!     PH=pH
!     FLWU,TUPWTR=total root water uptake from extract.f
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     RCO2O=net soil CO2 uptake from nitro.f
!     RCH4O=net soil CH4 uptake from nitro.f
!     RN2G=total soil N2 production from nitro.f
!     XN2GS=total N2 fixation from nitro.f
!     RN2O=net soil N2O uptake from nitro.f
!     RH2GO=net H2 uptake from nitro.f
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate from nitro.f
!     XNH4S=net change in NH4 from nitro.f
!     TRN4S,TRN3S=NH4,NH3 dissolution from solute.f
!     XNO3S=net change in NO3 from nitro.f
!     TRNO3=NO3 dissolution from solute.f
!     XNO2S=net change in NO2 from nitro.f
!     TRNO2=NO2 dissolution from solute.f
!     XH2PS=net change in H2PO4 from nitro.f
!     TRH2P=H2PO4 dissolution from solute.f
!     XH1PS=net change in HPO4 from nitro.f
!     TRH1P=HPO4 dissolution from solute.f
!     TUPNH4,TUPNHB=root NH4 uptake in non-band,band from extract.f
!     TUPNO3,TUPNOB=root NO3 uptake in non-band,band from extract.f
!     TUPH2P,TUPH2B=root H2PO4 uptake in non-band,band from extract.f
!     TUPH1P,TUPH1B=root HPO4 uptake in non-band,band from extract.f
!
      call ImportFluxFromOutsideModules(I,NY,NX)
    ENDDO
  ENDDO
  end subroutine InitFluxandStateVariables
!------------------------------------------------------------------------------------------

  subroutine UpdateStateVar(M,MM,MX,NPG,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) ::M, MM, MX,NPG, NHW, NHE, NVN, NVS
  integer :: NY,NX,K,L

  IF(MM.NE.NPG)THEN
    DO 9695 NX=NHW,NHE
      DO 9690 NY=NVN,NVS
        IF(M.NE.MX)THEN
!
!     STATE VARIABLES FOR SOLUTES IN SNOWPACK
!
!     *W2=solute content of snowpack
!     TQS*=net overland solute flux in snow
!     T*BLS=net solute flux in snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
!
          CO2W2(1,NY,NX)=CO2W2(1,NY,NX)+TQSCOS(NY,NX)
          CH4W2(1,NY,NX)=CH4W2(1,NY,NX)+TQSCHS(NY,NX)
          OXYW2(1,NY,NX)=OXYW2(1,NY,NX)+TQSOXS(NY,NX)
          ZNGW2(1,NY,NX)=ZNGW2(1,NY,NX)+TQSNGS(NY,NX)
          ZN2W2(1,NY,NX)=ZN2W2(1,NY,NX)+TQSN2S(NY,NX)
          ZN4W2(1,NY,NX)=ZN4W2(1,NY,NX)+TQSNH4(NY,NX)
          ZN3W2(1,NY,NX)=ZN3W2(1,NY,NX)+TQSNH3(NY,NX)
          ZNOW2(1,NY,NX)=ZNOW2(1,NY,NX)+TQSNO3(NY,NX)
          Z1PW2(1,NY,NX)=Z1PW2(1,NY,NX)+TQSH1P(NY,NX)
          ZHPW2(1,NY,NX)=ZHPW2(1,NY,NX)+TQSH2P(NY,NX)
          DO  L=1,JS
            CO2W2(L,NY,NX)=CO2W2(L,NY,NX)+TCOBLS(L,NY,NX)
            CH4W2(L,NY,NX)=CH4W2(L,NY,NX)+TCHBLS(L,NY,NX)
            OXYW2(L,NY,NX)=OXYW2(L,NY,NX)+TOXBLS(L,NY,NX)
            ZNGW2(L,NY,NX)=ZNGW2(L,NY,NX)+TNGBLS(L,NY,NX)
            ZN2W2(L,NY,NX)=ZN2W2(L,NY,NX)+TN2BLS(L,NY,NX)
            ZN4W2(L,NY,NX)=ZN4W2(L,NY,NX)+TN4BLW(L,NY,NX)
            ZN3W2(L,NY,NX)=ZN3W2(L,NY,NX)+TN3BLW(L,NY,NX)
            ZNOW2(L,NY,NX)=ZNOW2(L,NY,NX)+TNOBLW(L,NY,NX)
            Z1PW2(L,NY,NX)=Z1PW2(L,NY,NX)+TH1PBS(L,NY,NX)
            ZHPW2(L,NY,NX)=ZHPW2(L,NY,NX)+TH2PBS(L,NY,NX)
          ENDDO
!
!     STATE VARIABLES FOR SOLUTES IN SURFACE RESIDUE AND IN
!     MICROPORES AND MACROPORES IN SOIL SURFACE LAYER FROM OVERLAND
!     FLOW AND SURFACE VOLATILIZATION-DISSOLUTION
!
!     *S2=litter solute content
!     R*DFR=gas exchange between atmosphere and surface litter water
!     R*DFS=gas exchange between atmosphere and soil surface water
!     R*FLS=convective + diffusive solute flux into litter,soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into litter from non-band,band
!     TQR*=net overland solute flux
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
          DO  K=0,jcplx1
            OQC2(K,0,NY,NX)=OQC2(K,0,NY,NX)+ROCFLS(K,3,0,NY,NX)
            OQN2(K,0,NY,NX)=OQN2(K,0,NY,NX)+RONFLS(K,3,0,NY,NX)
            OQP2(K,0,NY,NX)=OQP2(K,0,NY,NX)+ROPFLS(K,3,0,NY,NX)
            OQA2(K,0,NY,NX)=OQA2(K,0,NY,NX)+ROAFLS(K,3,0,NY,NX)
          ENDDO
          CO2S2(0,NY,NX)=CO2S2(0,NY,NX)+RCODFR(NY,NX)+RCOFLS(3,0,NY,NX)
          CH4S2(0,NY,NX)=CH4S2(0,NY,NX)+RCHDFR(NY,NX)+RCHFLS(3,0,NY,NX)
          OXYS2(0,NY,NX)=OXYS2(0,NY,NX)+ROXDFR(NY,NX)+ROXFLS(3,0,NY,NX)
          Z2GS2(0,NY,NX)=Z2GS2(0,NY,NX)+RNGDFR(NY,NX)+RNGFLS(3,0,NY,NX)
          Z2OS2(0,NY,NX)=Z2OS2(0,NY,NX)+RN2DFR(NY,NX)+RN2FLS(3,0,NY,NX)
          H2GS2(0,NY,NX)=H2GS2(0,NY,NX)+RHGDFR(NY,NX)+RHGFLS(3,0,NY,NX)
          ZNH4S2(0,NY,NX)=ZNH4S2(0,NY,NX)+RN4FLW(3,0,NY,NX)
          ZNH3S2(0,NY,NX)=ZNH3S2(0,NY,NX)+RN3DFR(NY,NX)+RN3FLW(3,0,NY,NX)
          ZNO3S2(0,NY,NX)=ZNO3S2(0,NY,NX)+RNOFLW(3,0,NY,NX)
          ZNO2S2(0,NY,NX)=ZNO2S2(0,NY,NX)+RNXFLS(3,0,NY,NX)
          H1PO42(0,NY,NX)=H1PO42(0,NY,NX)+RH1PFS(3,0,NY,NX)
          H2PO42(0,NY,NX)=H2PO42(0,NY,NX)+RH2PFS(3,0,NY,NX)
          CO2S2(NU(NY,NX),NY,NX)=CO2S2(NU(NY,NX),NY,NX)+RCODFS(NY,NX)
          CH4S2(NU(NY,NX),NY,NX)=CH4S2(NU(NY,NX),NY,NX)+RCHDFS(NY,NX)
          OXYS2(NU(NY,NX),NY,NX)=OXYS2(NU(NY,NX),NY,NX)+ROXDFS(NY,NX)
          Z2GS2(NU(NY,NX),NY,NX)=Z2GS2(NU(NY,NX),NY,NX)+RNGDFS(NY,NX)
          Z2OS2(NU(NY,NX),NY,NX)=Z2OS2(NU(NY,NX),NY,NX)+RN2DFS(NY,NX)
          ZNH3S2(NU(NY,NX),NY,NX)=ZNH3S2(NU(NY,NX),NY,NX)+RN3DFS(NY,NX)
          ZNH3B2(NU(NY,NX),NY,NX)=ZNH3B2(NU(NY,NX),NY,NX)+RNBDFS(NY,NX)
          H2GS2(NU(NY,NX),NY,NX)=H2GS2(NU(NY,NX),NY,NX)+RHGDFS(NY,NX)
          DO 9680 K=0,jcplx1
            OQC2(K,0,NY,NX)=OQC2(K,0,NY,NX)+TQROC(K,NY,NX)
            OQN2(K,0,NY,NX)=OQN2(K,0,NY,NX)+TQRON(K,NY,NX)
            OQP2(K,0,NY,NX)=OQP2(K,0,NY,NX)+TQROP(K,NY,NX)
            OQA2(K,0,NY,NX)=OQA2(K,0,NY,NX)+TQROA(K,NY,NX)
9680  CONTINUE
          CO2S2(0,NY,NX)=CO2S2(0,NY,NX)+TQRCOS(NY,NX)
          CH4S2(0,NY,NX)=CH4S2(0,NY,NX)+TQRCHS(NY,NX)
          OXYS2(0,NY,NX)=OXYS2(0,NY,NX)+TQROXS(NY,NX)
          Z2GS2(0,NY,NX)=Z2GS2(0,NY,NX)+TQRNGS(NY,NX)
          Z2OS2(0,NY,NX)=Z2OS2(0,NY,NX)+TQRN2S(NY,NX)
          H2GS2(0,NY,NX)=H2GS2(0,NY,NX)+TQRHGS(NY,NX)
          ZNH4S2(0,NY,NX)=ZNH4S2(0,NY,NX)+TQRNH4(NY,NX)
          ZNH3S2(0,NY,NX)=ZNH3S2(0,NY,NX)+TQRNH3(NY,NX)
          ZNO3S2(0,NY,NX)=ZNO3S2(0,NY,NX)+TQRNO3(NY,NX)
          ZNO2S2(0,NY,NX)=ZNO2S2(0,NY,NX)+TQRNO2(NY,NX)
          H1PO42(0,NY,NX)=H1PO42(0,NY,NX)+TQRH1P(NY,NX)
          H2PO42(0,NY,NX)=H2PO42(0,NY,NX)+TQRH2P(NY,NX)
        ENDIF
!
        call UpdateSolutesInSoilLayers(M,MX,NY,NX)
!
!     GAS EXCHANGE IN SURFACE LITTER
!
!     *S2=litter aqueous gas content
!     R*DFG=litter-atmosphere gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
        CO2S2(0,NY,NX)=CO2S2(0,NY,NX)+RCODFG(0,NY,NX)
        CH4S2(0,NY,NX)=CH4S2(0,NY,NX)+RCHDFG(0,NY,NX)
        OXYS2(0,NY,NX)=OXYS2(0,NY,NX)+ROXDFG(0,NY,NX)
        Z2GS2(0,NY,NX)=Z2GS2(0,NY,NX)+RNGDFG(0,NY,NX)
        Z2OS2(0,NY,NX)=Z2OS2(0,NY,NX)+RN2DFG(0,NY,NX)
        ZNH3S2(0,NY,NX)=ZNH3S2(0,NY,NX)+RN3DFG(0,NY,NX)
        H2GS2(0,NY,NX)=H2GS2(0,NY,NX)+RHGDFG(0,NY,NX)
9690  CONTINUE
9695  CONTINUE
  ENDIF
  end subroutine UpdateStateVar

!------------------------------------------------------------------------------------------

  subroutine UpdateSolutesInSoilLayers(M,MX,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX,MX

  integer :: L,K
!     STATE VARIABLES FOR SOLUTES IN MICROPORES AND
!     MACROPORES IN SOIL LAYERS FROM SUBSURFACE FLOW, MICROBIAL
!     AND ROOT EXCHANGE IN 'NITRO' AND 'UPTAKE', AND EQUILIBRIUM
!     REACTIONS IN 'SOLUTE'
!
!     *S2,*B2=micropore solute content in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     T*FLS=net convective + diffusive solute flux through micropores
!     T*FHS=net convective + diffusive solute flux through macropores
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     R*BBL=bubble flux
!
  DO 9685 L=NU(NY,NX),NL(NY,NX)
    IF(M.NE.MX)THEN
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        CO2S2(L,NY,NX)=CO2S2(L,NY,NX)+TCOFLS(L,NY,NX)+RCOFXS(L,NY,NX) &
          +RCOFLZ(L,NY,NX)+RCOBBL(L,NY,NX)
        CH4S2(L,NY,NX)=CH4S2(L,NY,NX)+TCHFLS(L,NY,NX)+RCHFXS(L,NY,NX) &
          +RCHFLZ(L,NY,NX)+RCHBBL(L,NY,NX)
        OXYS2(L,NY,NX)=OXYS2(L,NY,NX)+TOXFLS(L,NY,NX)+ROXFXS(L,NY,NX) &
          +ROXFLZ(L,NY,NX)+ROXBBL(L,NY,NX)
        Z2GS2(L,NY,NX)=Z2GS2(L,NY,NX)+TNGFLS(L,NY,NX)+RNGFXS(L,NY,NX) &
          +RNGFLZ(L,NY,NX)+RNGBBL(L,NY,NX)
        Z2OS2(L,NY,NX)=Z2OS2(L,NY,NX)+TN2FLS(L,NY,NX)+RN2FXS(L,NY,NX) &
          +RN2FLZ(L,NY,NX)+RN2BBL(L,NY,NX)
        ZNH3S2(L,NY,NX)=ZNH3S2(L,NY,NX)+TN3FLW(L,NY,NX)+RN3FXW(L,NY,NX) &
          +RN3FLZ(L,NY,NX)+RN3BBL(L,NY,NX)
        ZNH3B2(L,NY,NX)=ZNH3B2(L,NY,NX)+TN3FLB(L,NY,NX)+RN3FXB(L,NY,NX) &
          +RN3FBZ(L,NY,NX)+RNBBBL(L,NY,NX)
        H2GS2(L,NY,NX)=H2GS2(L,NY,NX)+THGFLS(L,NY,NX)+RHGFXS(L,NY,NX) &
          +RHGFLZ(L,NY,NX)+RHGBBL(L,NY,NX)
        DO 9675 K=0,jcplx1
          OQC2(K,L,NY,NX)=OQC2(K,L,NY,NX)+TOCFLS(K,L,NY,NX)+ROCFXS(K,L,NY,NX)
          OQN2(K,L,NY,NX)=OQN2(K,L,NY,NX)+TONFLS(K,L,NY,NX)+RONFXS(K,L,NY,NX)
          OQP2(K,L,NY,NX)=OQP2(K,L,NY,NX)+TOPFLS(K,L,NY,NX)+ROPFXS(K,L,NY,NX)
          OQA2(K,L,NY,NX)=OQA2(K,L,NY,NX)+TOAFLS(K,L,NY,NX)+ROAFXS(K,L,NY,NX)
          OQCH2(K,L,NY,NX)=OQCH2(K,L,NY,NX)+TOCFHS(K,L,NY,NX)-ROCFXS(K,L,NY,NX)
          OQNH2(K,L,NY,NX)=OQNH2(K,L,NY,NX)+TONFHS(K,L,NY,NX)-RONFXS(K,L,NY,NX)
          OQPH2(K,L,NY,NX)=OQPH2(K,L,NY,NX)+TOPFHS(K,L,NY,NX)-ROPFXS(K,L,NY,NX)
          OQAH2(K,L,NY,NX)=OQAH2(K,L,NY,NX)+TOAFHS(K,L,NY,NX)-ROAFXS(K,L,NY,NX)
9675    CONTINUE
        ZNH4S2(L,NY,NX)=ZNH4S2(L,NY,NX)+TN4FLW(L,NY,NX)+RN4FXW(L,NY,NX)+RN4FLZ(L,NY,NX)
        ZNO3S2(L,NY,NX)=ZNO3S2(L,NY,NX)+TNOFLW(L,NY,NX)+RNOFXW(L,NY,NX)+RNOFLZ(L,NY,NX)
        ZNO2S2(L,NY,NX)=ZNO2S2(L,NY,NX)+TNXFLS(L,NY,NX)+RNXFXS(L,NY,NX)
        H1PO42(L,NY,NX)=H1PO42(L,NY,NX)+TH1PFS(L,NY,NX)+RH1PXS(L,NY,NX)+RH1PFZ(L,NY,NX)
        H2PO42(L,NY,NX)=H2PO42(L,NY,NX)+TH2PFS(L,NY,NX)+RH2PXS(L,NY,NX)+RH2PFZ(L,NY,NX)
        ZNH4B2(L,NY,NX)=ZNH4B2(L,NY,NX)+TN4FLB(L,NY,NX)+RN4FXB(L,NY,NX)+RN4FBZ(L,NY,NX)
        ZNO3B2(L,NY,NX)=ZNO3B2(L,NY,NX)+TNOFLB(L,NY,NX)+RNOFXB(L,NY,NX)+RNOFBZ(L,NY,NX)
        ZNO2B2(L,NY,NX)=ZNO2B2(L,NY,NX)+TNXFLB(L,NY,NX)+RNXFXB(L,NY,NX)
        H1POB2(L,NY,NX)=H1POB2(L,NY,NX)+TH1BFB(L,NY,NX)+RH1BXB(L,NY,NX)+RH1BBZ(L,NY,NX)
        H2POB2(L,NY,NX)=H2POB2(L,NY,NX)+TH2BFB(L,NY,NX)+RH2BXB(L,NY,NX)+RH2BBZ(L,NY,NX)
        CO2SH2(L,NY,NX)=CO2SH2(L,NY,NX)+TCOFHS(L,NY,NX)-RCOFXS(L,NY,NX)
        CH4SH2(L,NY,NX)=CH4SH2(L,NY,NX)+TCHFHS(L,NY,NX)-RCHFXS(L,NY,NX)
        OXYSH2(L,NY,NX)=OXYSH2(L,NY,NX)+TOXFHS(L,NY,NX)-ROXFXS(L,NY,NX)
        Z2GSH2(L,NY,NX)=Z2GSH2(L,NY,NX)+TNGFHS(L,NY,NX)-RNGFXS(L,NY,NX)
        Z2OSH2(L,NY,NX)=Z2OSH2(L,NY,NX)+TN2FHS(L,NY,NX)-RN2FXS(L,NY,NX)
        H2GSH2(L,NY,NX)=H2GSH2(L,NY,NX)+THGFHS(L,NY,NX)-RHGFXS(L,NY,NX)
        ZNH4H2(L,NY,NX)=ZNH4H2(L,NY,NX)+TN4FHW(L,NY,NX)-RN4FXW(L,NY,NX)
        ZNH3H2(L,NY,NX)=ZNH3H2(L,NY,NX)+TN3FHW(L,NY,NX)-RN3FXW(L,NY,NX)
        ZNO3H2(L,NY,NX)=ZNO3H2(L,NY,NX)+TNOFHW(L,NY,NX)-RNOFXW(L,NY,NX)
        ZNO2H2(L,NY,NX)=ZNO2H2(L,NY,NX)+TNXFHS(L,NY,NX)-RNXFXS(L,NY,NX)
        H1P4H2(L,NY,NX)=H1P4H2(L,NY,NX)+TH1PHS(L,NY,NX)-RH1PXS(L,NY,NX)
        H2P4H2(L,NY,NX)=H2P4H2(L,NY,NX)+TH2PHS(L,NY,NX)-RH2PXS(L,NY,NX)
        ZN4BH2(L,NY,NX)=ZN4BH2(L,NY,NX)+TN4FHB(L,NY,NX)-RN4FXB(L,NY,NX)
        ZN3BH2(L,NY,NX)=ZN3BH2(L,NY,NX)+TN3FHB(L,NY,NX)-RN3FXB(L,NY,NX)
        ZNOBH2(L,NY,NX)=ZNOBH2(L,NY,NX)+TNOFHB(L,NY,NX)-RNOFXB(L,NY,NX)
        ZN2BH2(L,NY,NX)=ZN2BH2(L,NY,NX)+TNXFHB(L,NY,NX)-RNXFXB(L,NY,NX)
        H1PBH2(L,NY,NX)=H1PBH2(L,NY,NX)+TH1BHB(L,NY,NX)-RH1BXB(L,NY,NX)
        H2PBH2(L,NY,NX)=H2PBH2(L,NY,NX)+TH2BHB(L,NY,NX)-RH2BXB(L,NY,NX)
      ENDIF
    ENDIF
!
!     STATE VARIABLES FOR GASES IN SOIL LAYERS FROM SUBSURFACE FLOW,
!     MICROBIAL AND ROOT EXCHANGE IN 'NITRO' AND 'UPTAKE'
!
!     *G2,*S2=soil gas, solute content
!     R*DFG=water-air gas flux
!     T*FLG=net convective+diffusive gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
    IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CO2S2(L,NY,NX)=CO2S2(L,NY,NX)+RCODFG(L,NY,NX)
      CH4S2(L,NY,NX)=CH4S2(L,NY,NX)+RCHDFG(L,NY,NX)
      OXYS2(L,NY,NX)=OXYS2(L,NY,NX)+ROXDFG(L,NY,NX)
      Z2GS2(L,NY,NX)=Z2GS2(L,NY,NX)+RNGDFG(L,NY,NX)
      Z2OS2(L,NY,NX)=Z2OS2(L,NY,NX)+RN2DFG(L,NY,NX)
      ZNH3S2(L,NY,NX)=ZNH3S2(L,NY,NX)+RN3DFG(L,NY,NX)
      ZNH3B2(L,NY,NX)=ZNH3B2(L,NY,NX)+RNBDFG(L,NY,NX)
      H2GS2(L,NY,NX)=H2GS2(L,NY,NX)+RHGDFG(L,NY,NX)
      CO2G2(L,NY,NX)=CO2G2(L,NY,NX)+TCOFLG(L,NY,NX)-RCODFG(L,NY,NX)
      CH4G2(L,NY,NX)=CH4G2(L,NY,NX)+TCHFLG(L,NY,NX)-RCHDFG(L,NY,NX)
      OXYG2(L,NY,NX)=OXYG2(L,NY,NX)+TOXFLG(L,NY,NX)-ROXDFG(L,NY,NX)
      Z2GG2(L,NY,NX)=Z2GG2(L,NY,NX)+TNGFLG(L,NY,NX)-RNGDFG(L,NY,NX)
      Z2OG2(L,NY,NX)=Z2OG2(L,NY,NX)+TN2FLG(L,NY,NX)-RN2DFG(L,NY,NX)
      ZN3G2(L,NY,NX)=ZN3G2(L,NY,NX)+TN3FLG(L,NY,NX)-RN3DFG(L,NY,NX)-RNBDFG(L,NY,NX)
      H2GG2(L,NY,NX)=H2GG2(L,NY,NX)+THGFLG(L,NY,NX)-RHGDFG(L,NY,NX)
    ENDIF

9685  CONTINUE
  end subroutine UpdateSolutesInSoilLayers


!------------------------------------------------------------------------------------------

  subroutine SurfaceSinksandSources(NY, NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K

  RCOSK2(0,NY,NX)=RCO2O(0,NY,NX)*XNPG
  RCHSK2(0,NY,NX)=RCH4O(0,NY,NX)*XNPG
  RNGSK2(0,NY,NX)=(RN2G(0,NY,NX)+XN2GS(0,NY,NX))*XNPG
  RN2SK2(0,NY,NX)=RN2O(0,NY,NX)*XNPG
  RNHSK2(0,NY,NX)=0.0_r8
  RHGSK2(0,NY,NX)=RH2GO(0,NY,NX)*XNPG
  DO  K=0,jcplx1
    ROCSK2(K,0,NY,NX)=-XOQCS(K,0,NY,NX)*XNPH
    RONSK2(K,0,NY,NX)=-XOQNS(K,0,NY,NX)*XNPH
    ROPSK2(K,0,NY,NX)=-XOQPS(K,0,NY,NX)*XNPH
    ROASK2(K,0,NY,NX)=-XOQAS(K,0,NY,NX)*XNPH
  ENDDO
  RN4SK2(0,NY,NX)=(-XNH4S(0,NY,NX)-TRN4S(0,NY,NX))*XNPH
  RN3SK2(0,NY,NX)=-TRN3S(0,NY,NX)*XNPH
  RNOSK2(0,NY,NX)=(-XNO3S(0,NY,NX)-TRNO3(0,NY,NX))*XNPH
  RNXSK2(0,NY,NX)=(-XNO2S(0,NY,NX)-TRNO2(0,NY,NX))*XNPH
  RHPSK2(0,NY,NX)=(-XH2PS(0,NY,NX)-TRH2P(0,NY,NX))*XNPH
  R1PSK2(0,NY,NX)=(-XH1PS(0,NY,NX)-TRH1P(0,NY,NX))*XNPH
  end subroutine SurfaceSinksandSources
!------------------------------------------------------------------------------------------

  subroutine StateVarforGasandSolute(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K

  CO2S2(0,NY,NX)=CO2S(0,NY,NX)
  CH4S2(0,NY,NX)=CH4S(0,NY,NX)
  OXYS2(0,NY,NX)=OXYS(0,NY,NX)
  Z2GS2(0,NY,NX)=Z2GS(0,NY,NX)
  Z2OS2(0,NY,NX)=Z2OS(0,NY,NX)
  H2GS2(0,NY,NX)=H2GS(0,NY,NX)
  DO 9979 K=0,jcplx1
    OQC2(K,0,NY,NX)=OQC(K,0,NY,NX)-XOQCS(K,0,NY,NX)
    OQN2(K,0,NY,NX)=OQN(K,0,NY,NX)-XOQNS(K,0,NY,NX)
    OQP2(K,0,NY,NX)=OQP(K,0,NY,NX)-XOQPS(K,0,NY,NX)
    OQA2(K,0,NY,NX)=OQA(K,0,NY,NX)-XOQAS(K,0,NY,NX)
9979  CONTINUE
  ZNH4S2(0,NY,NX)=ZNH4S(0,NY,NX)
  ZNH3S2(0,NY,NX)=ZNH3S(0,NY,NX)
  ZNO3S2(0,NY,NX)=ZNO3S(0,NY,NX)
  ZNO2S2(0,NY,NX)=ZNO2S(0,NY,NX)
  H1PO42(0,NY,NX)=H1PO4(0,NY,NX)
  H2PO42(0,NY,NX)=H2PO4(0,NY,NX)
  CHY0(0,NY,NX)=10.0**(-(PH(0,NY,NX)-3.0))
  end subroutine StateVarforGasandSolute
!------------------------------------------------------------------------------------------

  subroutine SurfaceSolutefromAtms(NY,NX)
  implicit none
  integer, intent(in) :: NY, NX
  integer :: K

  DO 8855 K=0,jcplx1
    IF(K.LE.2)THEN
      XOCFLS(K,3,0,NY,NX)=0.0
      XONFLS(K,3,0,NY,NX)=0.0
      XOPFLS(K,3,0,NY,NX)=0.0
      XOAFLS(K,3,0,NY,NX)=0.0
    ENDIF
    XOCFLS(K,3,NU(NY,NX),NY,NX)=0.0
    XONFLS(K,3,NU(NY,NX),NY,NX)=0.0
    XOPFLS(K,3,NU(NY,NX),NY,NX)=0.0
    XOAFLS(K,3,NU(NY,NX),NY,NX)=0.0
    XOCFHS(K,3,NU(NY,NX),NY,NX)=0.0
    XONFHS(K,3,NU(NY,NX),NY,NX)=0.0
    XOPFHS(K,3,NU(NY,NX),NY,NX)=0.0
    XOAFHS(K,3,NU(NY,NX),NY,NX)=0.0
8855  CONTINUE
  end subroutine SurfaceSolutefromAtms
!------------------------------------------------------------------------------------------

  subroutine HourlySoluteFluxes(I,NY,NX)
  implicit none

  integer, intent(in) :: I
  integer, intent(in) :: NY,NX

  IF(PRECW(NY,NX).GT.0.0.OR.(PRECR(NY,NX).GT.0.0 &
    .AND.VHCPWM(1,1,NY,NX).GT.VHCPWX(NY,NX)))THEN
    XCOBLS(1,NY,NX)=FLQGQ(NY,NX)*CCOR(NY,NX)+FLQGI(NY,NX)*CCOQ(NY,NX)
    XCHBLS(1,NY,NX)=FLQGQ(NY,NX)*CCHR(NY,NX)+FLQGI(NY,NX)*CCHQ(NY,NX)
    XOXBLS(1,NY,NX)=FLQGQ(NY,NX)*COXR(NY,NX)+FLQGI(NY,NX)*COXQ(NY,NX)
    XNGBLS(1,NY,NX)=FLQGQ(NY,NX)*CNNR(NY,NX)+FLQGI(NY,NX)*CNNQ(NY,NX)
    XN2BLS(1,NY,NX)=FLQGQ(NY,NX)*CN2R(NY,NX)+FLQGI(NY,NX)*CN2Q(NY,NX)
    XN4BLW(1,NY,NX)=(FLQGQ(NY,NX)*CN4R(NY,NX)+FLQGI(NY,NX)*CN4Q(I,NY,NX))*natomw
    XN3BLW(1,NY,NX)=(FLQGQ(NY,NX)*CN3R(NY,NX)+FLQGI(NY,NX)*CN3Q(I,NY,NX))*natomw
    XNOBLW(1,NY,NX)=(FLQGQ(NY,NX)*CNOR(NY,NX)+FLQGI(NY,NX)*CNOQ(I,NY,NX))*natomw
    XH1PBS(1,NY,NX)=(FLQGQ(NY,NX)*CH1PR(NY,NX)+FLQGI(NY,NX)*CH1PQ(I,NY,NX))*patomw
    XH2PBS(1,NY,NX)=(FLQGQ(NY,NX)*CPOR(NY,NX)+FLQGI(NY,NX)*CPOQ(I,NY,NX))*patomw
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IF RAINFALL AND IRRIGATION IS ZERO IF SNOWPACK IS PRESENT
!
!     X*FLS,X*FLB=hourly solute flux to micropores in non-band,band
!
    XCOFLS(3,0,NY,NX)=0.0_r8
    XCHFLS(3,0,NY,NX)=0.0_r8
    XOXFLS(3,0,NY,NX)=0.0_r8
    XNGFLS(3,0,NY,NX)=0.0_r8
    XN2FLS(3,0,NY,NX)=0.0_r8
    XHGFLS(3,0,NY,NX)=0.0_r8
    XN4FLW(3,0,NY,NX)=0.0_r8
    XN3FLW(3,0,NY,NX)=0.0_r8
    XNOFLW(3,0,NY,NX)=0.0_r8
    XNXFLS(3,0,NY,NX)=0.0_r8
    XH1PFS(3,0,NY,NX)=0.0_r8
    XH2PFS(3,0,NY,NX)=0.0_r8
    XCOFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XCHFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XOXFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XNGFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN2FLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XHGFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN4FLW(3,NU(NY,NX),NY,NX)=0.0_r8
    XN3FLW(3,NU(NY,NX),NY,NX)=0.0_r8
    XNOFLW(3,NU(NY,NX),NY,NX)=0.0_r8
    XNXFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XH1PFS(3,NU(NY,NX),NY,NX)=0.0_r8
    XH2PFS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN4FLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XN3FLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XNOFLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XNXFLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XH1BFB(3,NU(NY,NX),NY,NX)=0.0_r8
    XH2BFB(3,NU(NY,NX),NY,NX)=0.0_r8
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
  ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0) &
    .AND.VHCPWM(1,1,NY,NX).LE.VHCPWX(NY,NX))THEN
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IF SNOWFALL AND IRRIGATION IS ZERO AND SNOWPACK IS ABSENT
!
!     PRECQ,PRECI=snow+rain,irrigation
!     X*BLS=hourly solute flux to snowpack
!     X*FLS,X*FLB=hourly solute flux to surface litter,soil surface micropore non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to soil surface from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
    XCOBLS(1,NY,NX)=0.0_r8
    XCHBLS(1,NY,NX)=0.0_r8
    XOXBLS(1,NY,NX)=0.0_r8
    XNGBLS(1,NY,NX)=0.0_r8
    XN2BLS(1,NY,NX)=0.0_r8
    XN4BLW(1,NY,NX)=0.0_r8
    XN3BLW(1,NY,NX)=0.0_r8
    XNOBLW(1,NY,NX)=0.0_r8
    XH1PBS(1,NY,NX)=0.0_r8
    XH2PBS(1,NY,NX)=0.0_r8
    XCOFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CCOR(NY,NX)+FLQRI(NY,NX)*CCOQ(NY,NX)
    XCHFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CCHR(NY,NX)+FLQRI(NY,NX)*CCHQ(NY,NX)
    XOXFLS(3,0,NY,NX)=FLQRQ(NY,NX)*COXR(NY,NX)+FLQRI(NY,NX)*COXQ(NY,NX)
    XNGFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CNNR(NY,NX)+FLQRI(NY,NX)*CNNQ(NY,NX)
    XN2FLS(3,0,NY,NX)=FLQRQ(NY,NX)*CN2R(NY,NX)+FLQRI(NY,NX)*CN2Q(NY,NX)
    XHGFLS(3,0,NY,NX)=0.0_r8
    XN4FLW(3,0,NY,NX)=(FLQRQ(NY,NX)*CN4R(NY,NX)+FLQRI(NY,NX)*CN4Q(I,NY,NX))*natomw
    XN3FLW(3,0,NY,NX)=(FLQRQ(NY,NX)*CN3R(NY,NX)+FLQRI(NY,NX)*CN3Q(I,NY,NX))*natomw
    XNOFLW(3,0,NY,NX)=(FLQRQ(NY,NX)*CNOR(NY,NX)+FLQRI(NY,NX)*CNOQ(I,NY,NX))*natomw
    XNXFLS(3,0,NY,NX)=0.0_r8
    XH1PFS(3,0,NY,NX)=(FLQRQ(NY,NX)*CH1PR(NY,NX)+FLQRI(NY,NX)*CH1PQ(I,NY,NX))*patomw
    XH2PFS(3,0,NY,NX)=(FLQRQ(NY,NX)*CPOR(NY,NX)+FLQRI(NY,NX)*CPOQ(I,NY,NX))*patomw
    XCOFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCOR(NY,NX)+FLQGI(NY,NX)*CCOQ(NY,NX)
    XCHFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCHR(NY,NX)+FLQGI(NY,NX)*CCHQ(NY,NX)
    XOXFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*COXR(NY,NX)+FLQGI(NY,NX)*COXQ(NY,NX)
    XNGFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNNR(NY,NX)+FLQGI(NY,NX)*CNNQ(NY,NX)
    XN2FLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CN2R(NY,NX)+FLQGI(NY,NX)*CN2Q(NY,NX)
    XHGFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN4FLW(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN4R(NY,NX) &
      +FLQGI(NY,NX)*CN4Q(I,NY,NX))*natomw)*VLNH4(NU(NY,NX),NY,NX)
    XN3FLW(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN3R(NY,NX) &
      +FLQGI(NY,NX)*CN3Q(I,NY,NX))*natomw)*VLNH4(NU(NY,NX),NY,NX)
    XNOFLW(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CNOR(NY,NX) &
      +FLQGI(NY,NX)*CNOQ(I,NY,NX))*natomw)*VLNO3(NU(NY,NX),NY,NX)
    XNXFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XH1PFS(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CH1PR(NY,NX) &
      +FLQGI(NY,NX)*CH1PQ(I,NY,NX))*patomw)*VLPO4(NU(NY,NX),NY,NX)
    XH2PFS(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CPOR(NY,NX) &
      +FLQGI(NY,NX)*CPOQ(I,NY,NX))*patomw)*VLPO4(NU(NY,NX),NY,NX)
    XN4FLB(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN4R(NY,NX) &
      +FLQGI(NY,NX)*CN4Q(I,NY,NX))*natomw)*VLNHB(NU(NY,NX),NY,NX)
    XN3FLB(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN3R(NY,NX) &
      +FLQGI(NY,NX)*CN3Q(I,NY,NX))*natomw)*VLNHB(NU(NY,NX),NY,NX)
    XNOFLB(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CNOR(NY,NX) &
      +FLQGI(NY,NX)*CNOQ(I,NY,NX))*natomw)*VLNOB(NU(NY,NX),NY,NX)
    XNXFLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XH1BFB(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CH1PR(NY,NX) &
      +FLQGI(NY,NX)*CH1PQ(I,NY,NX))*patomw)*VLPOB(NU(NY,NX),NY,NX)
    XH2BFB(3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CPOR(NY,NX) &
      +FLQGI(NY,NX)*CPOQ(I,NY,NX))*patomw)*VLPOB(NU(NY,NX),NY,NX)
!
!     NO SOLUTE FLUXES FROM ATMOSPHERE
!
  ELSE
    XCOBLS(1,NY,NX)=0.0_r8
    XCHBLS(1,NY,NX)=0.0_r8
    XOXBLS(1,NY,NX)=0.0_r8
    XNGBLS(1,NY,NX)=0.0_r8
    XN2BLS(1,NY,NX)=0.0_r8
    XN4BLW(1,NY,NX)=0.0_r8
    XN3BLW(1,NY,NX)=0.0_r8
    XNOBLW(1,NY,NX)=0.0_r8
    XH1PBS(1,NY,NX)=0.0_r8
    XH2PBS(1,NY,NX)=0.0_r8
    XCOFLS(3,0,NY,NX)=0.0_r8
    XCHFLS(3,0,NY,NX)=0.0_r8
    XOXFLS(3,0,NY,NX)=0.0_r8
    XNGFLS(3,0,NY,NX)=0.0_r8
    XN2FLS(3,0,NY,NX)=0.0_r8
    XHGFLS(3,0,NY,NX)=0.0_r8
    XN4FLW(3,0,NY,NX)=0.0_r8
    XN3FLW(3,0,NY,NX)=0.0_r8
    XNOFLW(3,0,NY,NX)=0.0_r8
    XNXFLS(3,0,NY,NX)=0.0_r8
    XH1PFS(3,0,NY,NX)=0.0_r8
    XH2PFS(3,0,NY,NX)=0.0_r8
    XCOFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XCHFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XOXFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XNGFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN2FLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XHGFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN4FLW(3,NU(NY,NX),NY,NX)=0.0_r8
    XN3FLW(3,NU(NY,NX),NY,NX)=0.0_r8
    XNOFLW(3,NU(NY,NX),NY,NX)=0.0_r8
    XNXFLS(3,NU(NY,NX),NY,NX)=0.0_r8
    XH1PFS(3,NU(NY,NX),NY,NX)=0.0_r8
    XH2PFS(3,NU(NY,NX),NY,NX)=0.0_r8
    XN4FLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XN3FLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XNOFLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XNXFLB(3,NU(NY,NX),NY,NX)=0.0_r8
    XH1BFB(3,NU(NY,NX),NY,NX)=0.0_r8
    XH2BFB(3,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  XCOFHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XCHFHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XOXFHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XNGFHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XN2FHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XHGFHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XN4FHW(3,NU(NY,NX),NY,NX)=0.0_r8
  XN3FHW(3,NU(NY,NX),NY,NX)=0.0_r8
  XNOFHW(3,NU(NY,NX),NY,NX)=0.0_r8
  XH1PHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XH2PHS(3,NU(NY,NX),NY,NX)=0.0_r8
  XN4FHB(3,NU(NY,NX),NY,NX)=0.0_r8
  XN3FHB(3,NU(NY,NX),NY,NX)=0.0_r8
  XNOFHB(3,NU(NY,NX),NY,NX)=0.0_r8
  XNXFHB(3,NU(NY,NX),NY,NX)=0.0_r8
  XH1BHB(3,NU(NY,NX),NY,NX)=0.0_r8
  XH2BHB(3,NU(NY,NX),NY,NX)=0.0_r8
  XNXFHS(3,NU(NY,NX),NY,NX)=0.0_r8
  end subroutine HourlySoluteFluxes
!------------------------------------------------------------------------------------------

  subroutine SubHourlyFluxesFromSiteFile(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX
  INTEGER :: K,L

  DO  K=0,2
    ROCFL0(K,NY,NX)=XOCFLS(K,3,0,NY,NX)*XNPH
    RONFL0(K,NY,NX)=XONFLS(K,3,0,NY,NX)*XNPH
    ROPFL0(K,NY,NX)=XOPFLS(K,3,0,NY,NX)*XNPH
    ROAFL0(K,NY,NX)=XOAFLS(K,3,0,NY,NX)*XNPH
    ROCFL1(K,NY,NX)=XOCFLS(K,3,NU(NY,NX),NY,NX)*XNPH
    RONFL1(K,NY,NX)=XONFLS(K,3,NU(NY,NX),NY,NX)*XNPH
    ROPFL1(K,NY,NX)=XOPFLS(K,3,NU(NY,NX),NY,NX)*XNPH
    ROAFL1(K,NY,NX)=XOAFLS(K,3,NU(NY,NX),NY,NX)*XNPH
  enddo
  RCOBLS(1,NY,NX)=XCOBLS(1,NY,NX)*XNPH
  RCHBLS(1,NY,NX)=XCHBLS(1,NY,NX)*XNPH
  ROXBLS(1,NY,NX)=XOXBLS(1,NY,NX)*XNPH
  RNGBLS(1,NY,NX)=XNGBLS(1,NY,NX)*XNPH
  RN2BLS(1,NY,NX)=XN2BLS(1,NY,NX)*XNPH
  RN4BLW(1,NY,NX)=XN4BLW(1,NY,NX)*XNPH
  RN3BLW(1,NY,NX)=XN3BLW(1,NY,NX)*XNPH
  RNOBLW(1,NY,NX)=XNOBLW(1,NY,NX)*XNPH
  RH1PBS(1,NY,NX)=XH1PBS(1,NY,NX)*XNPH
  RH2PBS(1,NY,NX)=XH2PBS(1,NY,NX)*XNPH
  RCOFL0(NY,NX)=XCOFLS(3,0,NY,NX)*XNPH
  RCHFL0(NY,NX)=XCHFLS(3,0,NY,NX)*XNPH
  ROXFL0(NY,NX)=XOXFLS(3,0,NY,NX)*XNPH
  RNGFL0(NY,NX)=XNGFLS(3,0,NY,NX)*XNPH
  RN2FL0(NY,NX)=XN2FLS(3,0,NY,NX)*XNPH
  RHGFL0(NY,NX)=XHGFLS(3,0,NY,NX)*XNPH
  RN4FL0(NY,NX)=XN4FLW(3,0,NY,NX)*XNPH
  RN3FL0(NY,NX)=XN3FLW(3,0,NY,NX)*XNPH
  RNOFL0(NY,NX)=XNOFLW(3,0,NY,NX)*XNPH
  RNXFL0(NY,NX)=XNXFLS(3,0,NY,NX)*XNPH
  RH1PF0(NY,NX)=XH1PFS(3,0,NY,NX)*XNPH
  RH2PF0(NY,NX)=XH2PFS(3,0,NY,NX)*XNPH
  RCOFL1(NY,NX)=XCOFLS(3,NU(NY,NX),NY,NX)*XNPH
  RCHFL1(NY,NX)=XCHFLS(3,NU(NY,NX),NY,NX)*XNPH
  ROXFL1(NY,NX)=XOXFLS(3,NU(NY,NX),NY,NX)*XNPH
  RNGFL1(NY,NX)=XNGFLS(3,NU(NY,NX),NY,NX)*XNPH
  RN2FL1(NY,NX)=XN2FLS(3,NU(NY,NX),NY,NX)*XNPH
  RHGFL1(NY,NX)=XHGFLS(3,NU(NY,NX),NY,NX)*XNPH
  RN4FL1(NY,NX)=XN4FLW(3,NU(NY,NX),NY,NX)*XNPH
  RN3FL1(NY,NX)=XN3FLW(3,NU(NY,NX),NY,NX)*XNPH
  RNOFL1(NY,NX)=XNOFLW(3,NU(NY,NX),NY,NX)*XNPH
  RNXFL1(NY,NX)=XNXFLS(3,NU(NY,NX),NY,NX)*XNPH
  RH1PF1(NY,NX)=XH1PFS(3,NU(NY,NX),NY,NX)*XNPH
  RH2PF1(NY,NX)=XH2PFS(3,NU(NY,NX),NY,NX)*XNPH
  RN4FL2(NY,NX)=XN4FLB(3,NU(NY,NX),NY,NX)*XNPH
  RN3FL2(NY,NX)=XN3FLB(3,NU(NY,NX),NY,NX)*XNPH
  RNOFL2(NY,NX)=XNOFLB(3,NU(NY,NX),NY,NX)*XNPH
  RNXFL2(NY,NX)=XNXFLB(3,NU(NY,NX),NY,NX)*XNPH
  RH1BF2(NY,NX)=XH1BFB(3,NU(NY,NX),NY,NX)*XNPH
  RH2BF2(NY,NX)=XH2BFB(3,NU(NY,NX),NY,NX)*XNPH
!
!     GAS AND SOLUTE SINKS AND SOURCES IN SOIL LAYERS FROM MICROBIAL
!     TRANSFORMATIONS IN 'NITRO' + ROOT EXCHANGE IN 'EXTRACT'
!     + EQUILIBRIA REACTIONS IN 'SOLUTE' AT SUB-HOURLY TIME STEP
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     PARR=boundary layer conductance above litter from watsub.f
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
  CLSGL2(0,NY,NX)=CLSGL(0,NY,NX)*XNPH
  CQSGL2(0,NY,NX)=CQSGL(0,NY,NX)*XNPH
  OLSGL2(0,NY,NX)=OLSGL(0,NY,NX)*XNPH
  ZLSGL2(0,NY,NX)=ZLSGL(0,NY,NX)*XNPH
  ZNSGL2(0,NY,NX)=ZNSGL(0,NY,NX)*XNPH
  ZVSGL2(0,NY,NX)=ZVSGL(0,NY,NX)*XNPH
  HLSGL2(0,NY,NX)=HLSGL(0,NY,NX)*XNPH
  OCSGL2(0,NY,NX)=OCSGL(0,NY,NX)*XNPH
  ONSGL2(0,NY,NX)=ONSGL(0,NY,NX)*XNPH
  OPSGL2(0,NY,NX)=OPSGL(0,NY,NX)*XNPH
  OASGL2(0,NY,NX)=OASGL(0,NY,NX)*XNPH
  ZOSGL2(0,NY,NX)=ZOSGL(0,NY,NX)*XNPH
  POSGL2(0,NY,NX)=POSGL(0,NY,NX)*XNPH
!
!     INITIAL SOLUTES IN SNOWPACK
!
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  DO L=1,JS
    CO2W2(L,NY,NX)=CO2W(L,NY,NX)
    CH4W2(L,NY,NX)=CH4W(L,NY,NX)
    OXYW2(L,NY,NX)=OXYW(L,NY,NX)
    ZNGW2(L,NY,NX)=ZNGW(L,NY,NX)
    ZN2W2(L,NY,NX)=ZN2W(L,NY,NX)
    ZN4W2(L,NY,NX)=ZN4W(L,NY,NX)
    ZN3W2(L,NY,NX)=ZN3W(L,NY,NX)
    ZNOW2(L,NY,NX)=ZNOW(L,NY,NX)
    Z1PW2(L,NY,NX)=Z1PW(L,NY,NX)
    ZHPW2(L,NY,NX)=ZHPW(L,NY,NX)
  enddo
  end subroutine SubHourlyFluxesFromSiteFile
!------------------------------------------------------------------------------------------

  subroutine ImportFluxFromOutsideModules(I,NY,NX)
  implicit none

  integer, intent(in) ::  I,NY, NX
  real(r8) :: FLWU(JZ,JY,JX)
  integer :: L,K

  DO L=NU(NY,NX),NL(NY,NX)
    CHY0(L,NY,NX)=10.0_r8**(-(PH(L,NY,NX)-3.0_r8))
    FLWU(L,NY,NX)=TUPWTR(L,NY,NX)*XNPH
    RCOSK2(L,NY,NX)=(RCO2O(L,NY,NX)+TCO2S(L,NY,NX)-TRCO2(L,NY,NX))*XNPG
    RCHSK2(L,NY,NX)=(RCH4O(L,NY,NX)+TUPCHS(L,NY,NX))*XNPG
    RNGSK2(L,NY,NX)=(RN2G(L,NY,NX)+XN2GS(L,NY,NX)+TUPNF(L,NY,NX))*XNPG
    RN2SK2(L,NY,NX)=(RN2O(L,NY,NX)+TUPN2S(L,NY,NX))*XNPG
    RNHSK2(L,NY,NX)=-TRN3G(L,NY,NX)*XNPG
    RHGSK2(L,NY,NX)=(RH2GO(L,NY,NX)+TUPHGS(L,NY,NX))*XNPG
    DO  K=0,jcplx1
      ROCSK2(K,L,NY,NX)=-XOQCS(K,L,NY,NX)*XNPH
      RONSK2(K,L,NY,NX)=-XOQNS(K,L,NY,NX)*XNPH
      ROPSK2(K,L,NY,NX)=-XOQPS(K,L,NY,NX)*XNPH
      ROASK2(K,L,NY,NX)=-XOQAS(K,L,NY,NX)*XNPH
    ENDDO
    RN4SK2(L,NY,NX)=(-XNH4S(L,NY,NX)-TRN4S(L,NY,NX)+TUPNH4(L,NY,NX))*XNPH
    RN3SK2(L,NY,NX)=(-TRN3S(L,NY,NX)+TUPN3S(L,NY,NX))*XNPH
    RNOSK2(L,NY,NX)=(-XNO3S(L,NY,NX)-TRNO3(L,NY,NX)+TUPNO3(L,NY,NX))*XNPH
    RNXSK2(L,NY,NX)=(-XNO2S(L,NY,NX)-TRNO2(L,NY,NX))*XNPH
    RHPSK2(L,NY,NX)=(-XH2PS(L,NY,NX)-TRH2P(L,NY,NX)+TUPH2P(L,NY,NX))*XNPH
    R1PSK2(L,NY,NX)=(-XH1PS(L,NY,NX)-TRH1P(L,NY,NX)+TUPH1P(L,NY,NX))*XNPH
    R4BSK2(L,NY,NX)=(-XNH4B(L,NY,NX)-TRN4B(L,NY,NX)+TUPNHB(L,NY,NX))*XNPH
    R3BSK2(L,NY,NX)=(-TRN3B(L,NY,NX)+TUPN3B(L,NY,NX))*XNPH
    RNBSK2(L,NY,NX)=(-XNO3B(L,NY,NX)-TRNOB(L,NY,NX)+TUPNOB(L,NY,NX))*XNPH
    RNZSK2(L,NY,NX)=(-XNO2B(L,NY,NX)-TRN2B(L,NY,NX))*XNPH
    RHBSK2(L,NY,NX)=(-XH2BS(L,NY,NX)-TRH2B(L,NY,NX)+TUPH2B(L,NY,NX))*XNPH
    R1BSK2(L,NY,NX)=(-XH1BS(L,NY,NX)-TRH1B(L,NY,NX)+TUPH1B(L,NY,NX))*XNPH
!
!     SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     FLU=subsurface water flux from watsub.f
!     R*FLU,R*FBU=subsurface solute flux in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     C*Q=irrigation solute concentrations
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    RCOFLU(L,NY,NX)=FLU(L,NY,NX)*CCOQ(NY,NX)
    RCHFLU(L,NY,NX)=FLU(L,NY,NX)*CCHQ(NY,NX)
    ROXFLU(L,NY,NX)=FLU(L,NY,NX)*COXQ(NY,NX)
    RNGFLU(L,NY,NX)=FLU(L,NY,NX)*CNNQ(NY,NX)
    RN2FLU(L,NY,NX)=FLU(L,NY,NX)*CN2Q(NY,NX)
    RHGFLU(L,NY,NX)=0.0_r8
    RN4FLU(L,NY,NX)=FLU(L,NY,NX)*CN4Q(I,NY,NX)*VLNH4(L,NY,NX)*natomw
    RN3FLU(L,NY,NX)=FLU(L,NY,NX)*CN3Q(I,NY,NX)*VLNH4(L,NY,NX)*natomw
    RNOFLU(L,NY,NX)=FLU(L,NY,NX)*CNOQ(I,NY,NX)*VLNO3(L,NY,NX)*natomw
    RH1PFU(L,NY,NX)=FLU(L,NY,NX)*CH1PQ(I,NY,NX)*VLPO4(L,NY,NX)*patomw
    RH2PFU(L,NY,NX)=FLU(L,NY,NX)*CPOQ(I,NY,NX)*VLPO4(L,NY,NX)*patomw
    RN4FBU(L,NY,NX)=FLU(L,NY,NX)*CN4Q(I,NY,NX)*VLNHB(L,NY,NX)*natomw
    RN3FBU(L,NY,NX)=FLU(L,NY,NX)*CN3Q(I,NY,NX)*VLNHB(L,NY,NX)*natomw
    RNOFBU(L,NY,NX)=FLU(L,NY,NX)*CNOQ(I,NY,NX)*VLNOB(L,NY,NX)*natomw
    RH1BBU(L,NY,NX)=FLU(L,NY,NX)*CH1PQ(I,NY,NX)*VLPOB(L,NY,NX)*patomw
    RH2BBU(L,NY,NX)=FLU(L,NY,NX)*CPOQ(I,NY,NX)*VLPOB(L,NY,NX)*patomw
!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!
    RCOFLZ(L,NY,NX)=RCOFLU(L,NY,NX)*XNPH
    RCHFLZ(L,NY,NX)=RCHFLU(L,NY,NX)*XNPH
    ROXFLZ(L,NY,NX)=ROXFLU(L,NY,NX)*XNPH
    RNGFLZ(L,NY,NX)=RNGFLU(L,NY,NX)*XNPH
    RN2FLZ(L,NY,NX)=RN2FLU(L,NY,NX)*XNPH
    RHGFLZ(L,NY,NX)=RHGFLU(L,NY,NX)*XNPH
    RN4FLZ(L,NY,NX)=RN4FLU(L,NY,NX)*XNPH
    RN3FLZ(L,NY,NX)=RN3FLU(L,NY,NX)*XNPH
    RNOFLZ(L,NY,NX)=RNOFLU(L,NY,NX)*XNPH
    RH1PFZ(L,NY,NX)=RH1PFU(L,NY,NX)*XNPH
    RH2PFZ(L,NY,NX)=RH2PFU(L,NY,NX)*XNPH
    RN4FBZ(L,NY,NX)=RN4FBU(L,NY,NX)*XNPH
    RN3FBZ(L,NY,NX)=RN3FBU(L,NY,NX)*XNPH
    RNOFBZ(L,NY,NX)=RNOFBU(L,NY,NX)*XNPH
    RH1BBZ(L,NY,NX)=RH1BBU(L,NY,NX)*XNPH
    RH2BBZ(L,NY,NX)=RH2BBU(L,NY,NX)*XNPH
!
!     GAS AND SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     XNPG=1/number of cycles h-1 for gas flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
    OCSGL2(L,NY,NX)=OCSGL(L,NY,NX)*XNPH
    ONSGL2(L,NY,NX)=ONSGL(L,NY,NX)*XNPH
    OPSGL2(L,NY,NX)=OPSGL(L,NY,NX)*XNPH
    OASGL2(L,NY,NX)=OASGL(L,NY,NX)*XNPH
    CLSGL2(L,NY,NX)=CLSGL(L,NY,NX)*XNPH
    CQSGL2(L,NY,NX)=CQSGL(L,NY,NX)*XNPH
    OLSGL2(L,NY,NX)=OLSGL(L,NY,NX)*XNPH
    ZLSGL2(L,NY,NX)=ZLSGL(L,NY,NX)*XNPH
    ZVSGL2(L,NY,NX)=ZVSGL(L,NY,NX)*XNPH
    ZNSGL2(L,NY,NX)=ZNSGL(L,NY,NX)*XNPH
    HLSGL2(L,NY,NX)=HLSGL(L,NY,NX)*XNPH
    ZOSGL2(L,NY,NX)=ZOSGL(L,NY,NX)*XNPH
    POSGL2(L,NY,NX)=POSGL(L,NY,NX)*XNPH
    CGSGL2(L,NY,NX)=CGSGL(L,NY,NX)*XNPG
    CHSGL2(L,NY,NX)=CHSGL(L,NY,NX)*XNPG
    OGSGL2(L,NY,NX)=OGSGL(L,NY,NX)*XNPG
    ZGSGL2(L,NY,NX)=ZGSGL(L,NY,NX)*XNPG
    Z2SGL2(L,NY,NX)=Z2SGL(L,NY,NX)*XNPG
    ZHSGL2(L,NY,NX)=ZHSGL(L,NY,NX)*XNPG
    HGSGL2(L,NY,NX)=HGSGL(L,NY,NX)*XNPG
!
!     STATE VARIABLES FOR GASES AND SOLUTES USED IN 'TRNSFR'
!     TO STORE SUB-HOURLY CHANGES DURING FLUX CALCULATIONS
!     INCLUDING TRANSFORMATIONS FROM NITRO.F, UPTAKE.F AND SOLUTE.F
!
!     CO2G,CH4G,OXYG,ZN3G,Z2GG,Z2OG,H2GG=gaseous CO2,CH4,O2,NH3,N2,N2O,H2
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in micropores
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate from nitro.f
!     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!     CO2SH,CH4SH,OXYSH,Z2GSH,Z2OSH,H2GSH=aqueous CO2,CH4,O2,N2,N2O,H2 content in macropores
!     ZNH4SH,ZNH3SH,ZNO3SH,ZNO2SH,H1PO4H,H2PO4H=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band macropores
!     ZNH4BH,ZNH3BH,ZNO3BH,ZNO2BH,H1POBH,H2POBH=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band macropores
!
    CO2G2(L,NY,NX)=CO2G(L,NY,NX)
    CH4G2(L,NY,NX)=CH4G(L,NY,NX)
    OXYG2(L,NY,NX)=OXYG(L,NY,NX)
    ZN3G2(L,NY,NX)=ZNH3G(L,NY,NX)
    Z2GG2(L,NY,NX)=Z2GG(L,NY,NX)
    Z2OG2(L,NY,NX)=Z2OG(L,NY,NX)
    H2GG2(L,NY,NX)=H2GG(L,NY,NX)
    CO2S2(L,NY,NX)=CO2S(L,NY,NX)
    CH4S2(L,NY,NX)=CH4S(L,NY,NX)
    OXYS2(L,NY,NX)=OXYS(L,NY,NX)
    Z2GS2(L,NY,NX)=Z2GS(L,NY,NX)
    Z2OS2(L,NY,NX)=Z2OS(L,NY,NX)
    H2GS2(L,NY,NX)=H2GS(L,NY,NX)
    DO  K=0,jcplx1
      OQC2(K,L,NY,NX)=OQC(K,L,NY,NX)-XOQCS(K,L,NY,NX)
      OQN2(K,L,NY,NX)=OQN(K,L,NY,NX)-XOQNS(K,L,NY,NX)
      OQP2(K,L,NY,NX)=OQP(K,L,NY,NX)-XOQPS(K,L,NY,NX)
      OQA2(K,L,NY,NX)=OQA(K,L,NY,NX)-XOQAS(K,L,NY,NX)
      OQCH2(K,L,NY,NX)=OQCH(K,L,NY,NX)
      OQNH2(K,L,NY,NX)=OQNH(K,L,NY,NX)
      OQPH2(K,L,NY,NX)=OQPH(K,L,NY,NX)
      OQAH2(K,L,NY,NX)=OQAH(K,L,NY,NX)
    enddo
    ZNH4S2(L,NY,NX)=ZNH4S(L,NY,NX)
    ZNH3S2(L,NY,NX)=ZNH3S(L,NY,NX)
    ZNO3S2(L,NY,NX)=ZNO3S(L,NY,NX)
    ZNO2S2(L,NY,NX)=ZNO2S(L,NY,NX)
    H1PO42(L,NY,NX)=H1PO4(L,NY,NX)
    H2PO42(L,NY,NX)=H2PO4(L,NY,NX)
    ZNH4B2(L,NY,NX)=ZNH4B(L,NY,NX)
    ZNH3B2(L,NY,NX)=ZNH3B(L,NY,NX)
    ZNO3B2(L,NY,NX)=ZNO3B(L,NY,NX)
    ZNO2B2(L,NY,NX)=ZNO2B(L,NY,NX)
    H1POB2(L,NY,NX)=H1POB(L,NY,NX)
    H2POB2(L,NY,NX)=H2POB(L,NY,NX)
    CO2SH2(L,NY,NX)=CO2SH(L,NY,NX)
    CH4SH2(L,NY,NX)=CH4SH(L,NY,NX)
    OXYSH2(L,NY,NX)=OXYSH(L,NY,NX)
    Z2GSH2(L,NY,NX)=Z2GSH(L,NY,NX)
    Z2OSH2(L,NY,NX)=Z2OSH(L,NY,NX)
    H2GSH2(L,NY,NX)=H2GSH(L,NY,NX)
    ZNH4H2(L,NY,NX)=ZNH4SH(L,NY,NX)
    ZNH3H2(L,NY,NX)=ZNH3SH(L,NY,NX)
    ZNO3H2(L,NY,NX)=ZNO3SH(L,NY,NX)
    ZNO2H2(L,NY,NX)=ZNO2SH(L,NY,NX)
    H1P4H2(L,NY,NX)=H1PO4H(L,NY,NX)
    H2P4H2(L,NY,NX)=H2PO4H(L,NY,NX)
    ZN4BH2(L,NY,NX)=ZNH4BH(L,NY,NX)
    ZN3BH2(L,NY,NX)=ZNH3BH(L,NY,NX)
    ZNOBH2(L,NY,NX)=ZNO3BH(L,NY,NX)
    ZN2BH2(L,NY,NX)=ZNO2BH(L,NY,NX)
    H1PBH2(L,NY,NX)=H1POBH(L,NY,NX)
    H2PBH2(L,NY,NX)=H2POBH(L,NY,NX)
  enddo
  end subroutine ImportFluxFromOutsideModules
!

end module TrnsfrMod
