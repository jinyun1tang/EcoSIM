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

  real(r8), allocatable :: CHY0(:,:,:)
  real(r8), allocatable :: RHGFLZ(:,:,:)
  real(r8), allocatable :: RCOFLZ(:,:,:)
  real(r8), allocatable :: RCHFLZ(:,:,:)
  real(r8), allocatable :: ROXFLZ(:,:,:)
  real(r8), allocatable :: RNGFLZ(:,:,:)
  real(r8), allocatable :: RN2FLZ(:,:,:)
  real(r8), allocatable :: RN4FLZ(:,:,:)
  real(r8), allocatable :: RN3FLZ(:,:,:)
  real(r8), allocatable :: RNOFLZ(:,:,:)
  real(r8), allocatable :: RH2PFZ(:,:,:)
  real(r8), allocatable :: RN4FBZ(:,:,:)
  real(r8), allocatable :: RN3FBZ(:,:,:)
  real(r8), allocatable :: RNOFBZ(:,:,:)
  real(r8), allocatable :: RH2BBZ(:,:,:)
  real(r8), allocatable :: RH1PFZ(:,:,:)
  real(r8), allocatable :: RH1BBZ(:,:,:)

  real(r8), PARAMETER :: DPN4=5.7E-07
  REAL(r8) :: CCO2SQ,CCH4SQ,COXYSQ,CZ2GSQ,CZ2OSQ,CNH3SQ,CNH3BQ,CH2GSQ

  public :: trnsfr
  public :: InitTrnsfr,DestructTrnsfr
  contains

  subroutine InitTrnsfr()
  implicit none

  allocate(CHY0(0:JZ,JY,JX))
  allocate(RHGFLZ(JZ,JY,JX))
  allocate(RCOFLZ(JZ,JY,JX))
  allocate(RCHFLZ(JZ,JY,JX))
  allocate(ROXFLZ(JZ,JY,JX))
  allocate(RNGFLZ(JZ,JY,JX))
  allocate(RN2FLZ(JZ,JY,JX))
  allocate(RN4FLZ(JZ,JY,JX))
  allocate(RN3FLZ(JZ,JY,JX))
  allocate(RNOFLZ(JZ,JY,JX))
  allocate(RH2PFZ(JZ,JY,JX))
  allocate(RN4FBZ(JZ,JY,JX))
  allocate(RN3FBZ(JZ,JY,JX))
  allocate(RNOFBZ(JZ,JY,JX))
  allocate(RH2BBZ(JZ,JY,JX))
  allocate(RH1PFZ(JZ,JY,JX))
  allocate(RH1BBZ(JZ,JY,JX))

  call InitTransfrData

  call InitInsTp
  end subroutine InitTrnsfr
!------------------------------------------------------------------------------------------
  subroutine DestructTrnsfr
  implicit none


  call destroy(CHY0)
  call destroy(RHGFLZ)
  call destroy(RCOFLZ)
  call destroy(RCHFLZ)
  call destroy(ROXFLZ)
  call destroy(RNGFLZ)
  call destroy(RN2FLZ)
  call destroy(RN4FLZ)
  call destroy(RN3FLZ)
  call destroy(RNOFLZ)
  call destroy(RH2PFZ)
  call destroy(RN4FBZ)
  call destroy(RN3FBZ)
  call destroy(RNOFBZ)
  call destroy(RH2BBZ)
  call destroy(RH1PFZ)
  call destroy(RH1BBZ)

  call DestructTransfrData

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
  real(r8) :: FLQM(3,JD,JV,JH)
!     execution begins here
!
!     TIME STEPS FOR SOLUTE AND GAS FLUX CALCULATIONS
!

  FLQM(:,:,:,:)=0._r8
  call InitFluxandStateVariables(I,NHW,NHE,NVN,NVS)

!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=number of cycles per hour for gas flux calculations
! XNPT=1/number of cycles NPH-1 for gas flux calculations
  MX=0
  DO  MM=1,NPG
    M=MIN(NPH,INT((MM-1)*XNPT)+1)

    call ModelTracerHydroFlux(M,MX,NHW, NHE, NVN, NVS,FLQM)
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
  integer :: NY,NX,K,L,NGT

  IF(MM.NE.NPG)THEN
    D9695: DO NX=NHW,NHE
      D9690: DO NY=NVN,NVS
        IF(M.NE.MX)THEN
!
!     STATE VARIABLES FOR SOLUTES IN SNOWPACK
          call UpdateStateVarInSnowpack(NY,NX)
!
          call UpdateStateVarInResidue(NY,NX)

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
        DO NGT=idg_beg,idg_end
          trc_solml2(NGT,0,NY,NX)=trc_solml2(NGT,0,NY,NX)+RGasDSFlx(NGT,0,NY,NX)
        ENDDO
      ENDDO D9690
    ENDDO D9695
  ENDIF
  end subroutine UpdateStateVar
!------------------------------------------------------------------------------------------

  subroutine UpdateStateVarInSnowpack(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,NTG,NTN
!
!     *W2=solute content of snowpack
!     TQS*=net overland solute flux in snow
!     T*BLS=net solute flux in snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
!
  DO NTG=idg_beg,idg_end-1
    trcg_solsml2(NTG,1,NY,NX)=trcg_solsml2(NTG,1,NY,NX)+trcg_TQ(NTG,NY,NX)
  ENDDO

  DO NTN=ids_nut_beg,ids_nuts_end
    trcn_solsml2(NTN,1,NY,NX)=trcn_solsml2(NTN,1,NY,NX)+trcn_TQ(NTN,NY,NX)
  ENDDO

  DO  L=1,JS
    DO NTG=idg_beg,idg_end-1
      trcg_solsml2(NTG,L,NY,NX)=trcg_solsml2(NTG,L,NY,NX)+trcg_TBLS(NTG,L,NY,NX)
    ENDDO

    DO NTN=ids_nut_beg,ids_nuts_end
      trcn_solsml2(NTN,L,NY,NX)=trcn_solsml2(NTN,L,NY,NX)+trcn_TBLS(NTN,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine UpdateStateVarInSnowpack

!------------------------------------------------------------------------------------------

  subroutine UpdateStateVarInResidue(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,NTG,NTS
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
  DO  K=1,jcplx
    OQC2(K,0,NY,NX)=OQC2(K,0,NY,NX)+ROCFLS(K,3,0,NY,NX)
    OQN2(K,0,NY,NX)=OQN2(K,0,NY,NX)+RONFLS(K,3,0,NY,NX)
    OQP2(K,0,NY,NX)=OQP2(K,0,NY,NX)+ROPFLS(K,3,0,NY,NX)
    OQA2(K,0,NY,NX)=OQA2(K,0,NY,NX)+ROAFLS(K,3,0,NY,NX)
  ENDDO
  trc_solml2(idg_CO2,0,NY,NX)=trc_solml2(idg_CO2,0,NY,NX)+RCODFR(NY,NX)
  trc_solml2(idg_CH4,0,NY,NX)=trc_solml2(idg_CH4,0,NY,NX)+RCHDFR(NY,NX)
  trc_solml2(idg_O2,0,NY,NX)=trc_solml2(idg_O2,0,NY,NX)+ROXDFR(NY,NX)
  trc_solml2(idg_N2,0,NY,NX)=trc_solml2(idg_N2,0,NY,NX)+RNGDFR(NY,NX)
  trc_solml2(idg_N2O,0,NY,NX)=trc_solml2(idg_N2O,0,NY,NX)+RN2DFR(NY,NX)
  trc_solml2(idg_H2,0,NY,NX)=trc_solml2(idg_H2,0,NY,NX)+RHGDFR(NY,NX)
  trc_solml2(idg_NH3,0,NY,NX)=trc_solml2(idg_NH3,0,NY,NX)+RN3DFR(NY,NX)

! band does not exist in litter layer
  DO NTS=ids_beg,ids_end
    trc_solml2(NTS,0,NY,NX)=trc_solml2(NTS,0,NY,NX)+R3PoreSolFlx(NTS,3,0,NY,NX)
  ENDDO

  DO NTG=idg_beg,idg_end
    trc_solml2(NTG,NU(NY,NX),NY,NX)=trc_solml2(NTG,NU(NY,NX),NY,NX)+RGasSSVol(NTG,NY,NX)
  ENDDO
  D9680: DO K=1,jcplx
    OQC2(K,0,NY,NX)=OQC2(K,0,NY,NX)+TQROC(K,NY,NX)
    OQN2(K,0,NY,NX)=OQN2(K,0,NY,NX)+TQRON(K,NY,NX)
    OQP2(K,0,NY,NX)=OQP2(K,0,NY,NX)+TQROP(K,NY,NX)
    OQA2(K,0,NY,NX)=OQA2(K,0,NY,NX)+TQROA(K,NY,NX)
  ENDDO D9680


  trc_solml2(idg_CO2,0,NY,NX)=trc_solml2(idg_CO2,0,NY,NX)+trcg_TQR(idg_CO2,NY,NX)
  trc_solml2(idg_CH4,0,NY,NX)=trc_solml2(idg_CH4,0,NY,NX)+trcg_TQR(idg_CH4,NY,NX)
  trc_solml2(idg_O2,0,NY,NX)=trc_solml2(idg_O2,0,NY,NX)+trcg_TQR(idg_O2,NY,NX)
  trc_solml2(idg_N2,0,NY,NX)=trc_solml2(idg_N2,0,NY,NX)+trcg_TQR(idg_N2,NY,NX)
  trc_solml2(idg_N2O,0,NY,NX)=trc_solml2(idg_N2O,0,NY,NX)+trcg_TQR(idg_N2O,NY,NX)
  trc_solml2(idg_H2,0,NY,NX)=trc_solml2(idg_H2,0,NY,NX)+trcg_TQR(idg_H2,NY,NX)
  trc_solml2(idg_NH3,0,NY,NX)=trc_solml2(idg_NH3,0,NY,NX)+trcg_TQR(idg_NH3,NY,NX)

  trc_solml2(ids_NH4,0,NY,NX)=trc_solml2(ids_NH4,0,NY,NX)+trcn_TQR(ids_NH4,NY,NX)
  trc_solml2(ids_NO3,0,NY,NX)=trc_solml2(ids_NO3,0,NY,NX)+trcn_TQR(ids_NO3,NY,NX)
  trc_solml2(ids_NO2,0,NY,NX)=trc_solml2(ids_NO2,0,NY,NX)+trcn_TQR(ids_NO2,NY,NX)
  trc_solml2(ids_H1PO4,0,NY,NX)=trc_solml2(ids_H1PO4,0,NY,NX)+trcn_TQR(ids_H1PO4,NY,NX)
  trc_solml2(ids_H2PO4,0,NY,NX)=trc_solml2(ids_H2PO4,0,NY,NX)+trcn_TQR(ids_H2PO4,NY,NX)

  end subroutine UpdateStateVarInResidue

!------------------------------------------------------------------------------------------

  subroutine UpdateSolutesInSoilLayers(M,MX,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX,MX

  integer :: L,K,NTG,NTS

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
  D9685: DO L=NU(NY,NX),NL(NY,NX)
    IF(M.NE.MX)THEN
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN

        trc_solml2(idg_CO2,L,NY,NX)=trc_solml2(idg_CO2,L,NY,NX)+RCOFLZ(L,NY,NX)+RCOBBL(L,NY,NX)
        trc_solml2(idg_CH4,L,NY,NX)=trc_solml2(idg_CH4,L,NY,NX)+RCHFLZ(L,NY,NX)+RCHBBL(L,NY,NX)
        trc_solml2(idg_O2,L,NY,NX)=trc_solml2(idg_O2,L,NY,NX)+ROXFLZ(L,NY,NX)+ROXBBL(L,NY,NX)
        trc_solml2(idg_N2,L,NY,NX)=trc_solml2(idg_N2,L,NY,NX)+RNGFLZ(L,NY,NX)+RNGBBL(L,NY,NX)
        trc_solml2(idg_N2O,L,NY,NX)=trc_solml2(idg_N2O,L,NY,NX)+RN2FLZ(L,NY,NX)+RN2BBL(L,NY,NX)
        trc_solml2(idg_NH3,L,NY,NX)=trc_solml2(idg_NH3,L,NY,NX)+RN3FLZ(L,NY,NX)+RN3BBL(L,NY,NX)
        trc_solml2(idg_NH3B,L,NY,NX)=trc_solml2(idg_NH3B,L,NY,NX)+RN3FBZ(L,NY,NX)+RNBBBL(L,NY,NX)
        trc_solml2(idg_H2,L,NY,NX)=trc_solml2(idg_H2,L,NY,NX)+RHGFLZ(L,NY,NX)+RHGBBL(L,NY,NX)

        trc_solml2(ids_NH4,L,NY,NX)=trc_solml2(ids_NH4,L,NY,NX)+RN4FLZ(L,NY,NX)
        trc_solml2(ids_NO3,L,NY,NX)=trc_solml2(ids_NO3,L,NY,NX)+RNOFLZ(L,NY,NX)
        trc_solml2(ids_H1PO4,L,NY,NX)=trc_solml2(ids_H1PO4,L,NY,NX)+RH1PFZ(L,NY,NX)
        trc_solml2(ids_H2PO4,L,NY,NX)=trc_solml2(ids_H2PO4,L,NY,NX)+RH2PFZ(L,NY,NX)
        trc_solml2(ids_NH4B,L,NY,NX)=trc_solml2(ids_NH4B,L,NY,NX)+RN4FBZ(L,NY,NX)
        trc_solml2(ids_NO3B,L,NY,NX)=trc_solml2(ids_NO3B,L,NY,NX)+RNOFBZ(L,NY,NX)
        trc_solml2(ids_H1PO4B,L,NY,NX)=trc_solml2(ids_H1PO4B,L,NY,NX)+RH1BBZ(L,NY,NX)
        trc_solml2(ids_H2PO4B,L,NY,NX)=trc_solml2(ids_H2PO4B,L,NY,NX)+RH2BBZ(L,NY,NX)

        DO  K=1,jcplx
          OQC2(K,L,NY,NX)=OQC2(K,L,NY,NX)+TOCFLS(K,L,NY,NX)+ROCFXS(K,L,NY,NX)
          OQN2(K,L,NY,NX)=OQN2(K,L,NY,NX)+TONFLS(K,L,NY,NX)+RONFXS(K,L,NY,NX)
          OQP2(K,L,NY,NX)=OQP2(K,L,NY,NX)+TOPFLS(K,L,NY,NX)+ROPFXS(K,L,NY,NX)
          OQA2(K,L,NY,NX)=OQA2(K,L,NY,NX)+TOAFLS(K,L,NY,NX)+ROAFXS(K,L,NY,NX)
          OQCH2(K,L,NY,NX)=OQCH2(K,L,NY,NX)+TOCFHS(K,L,NY,NX)-ROCFXS(K,L,NY,NX)
          OQNH2(K,L,NY,NX)=OQNH2(K,L,NY,NX)+TONFHS(K,L,NY,NX)-RONFXS(K,L,NY,NX)
          OQPH2(K,L,NY,NX)=OQPH2(K,L,NY,NX)+TOPFHS(K,L,NY,NX)-ROPFXS(K,L,NY,NX)
          OQAH2(K,L,NY,NX)=OQAH2(K,L,NY,NX)+TOAFHS(K,L,NY,NX)-ROAFXS(K,L,NY,NX)
        ENDDO
        DO NTS=ids_beg,ids_end
          trc_solml2(NTS,L,NY,NX)=trc_solml2(NTS,L,NY,NX)+R3PorTSolFlx(NTS,L,NY,NX)+RporeSoXFlx(NTS,L,NY,NX)
          trc_soHml2(NTS,L,NY,NX)=trc_soHml2(NTS,L,NY,NX)+R3PorTSoHFlx(NTS,L,NY,NX)-RporeSoXFlx(NTS,L,NY,NX)
        ENDDO

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
      DO NTG=idg_beg,idg_end
        trc_solml2(NTG,L,NY,NX)=trc_solml2(NTG,L,NY,NX)+RGasDSFlx(NTG,L,NY,NX)
      ENDDO

      DO NTG=idg_beg,idg_end-1
        trc_gasml2(NTG,L,NY,NX)=trc_gasml2(NTG,L,NY,NX)+RTGasADFlx(NTG,L,NY,NX)-RGasDSFlx(NTG,L,NY,NX)
      ENDDO
      trc_gasml2(idg_NH3,L,NY,NX)=trc_gasml2(idg_NH3,L,NY,NX)-RGasDSFlx(idg_NH3B,L,NY,NX)
    ENDIF

  ENDDO D9685
  end subroutine UpdateSolutesInSoilLayers


!------------------------------------------------------------------------------------------

  subroutine SurfaceSinksandSources(NY, NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K

  RBGCSinkG(idg_CO2,0,NY,NX)=RCO2O(0,NY,NX)*XNPG
  RBGCSinkG(idg_CH4,0,NY,NX)=RCH4O(0,NY,NX)*XNPG
  RBGCSinkG(idg_N2,0,NY,NX)=(RN2G(0,NY,NX)+XN2GS(0,NY,NX))*XNPG
  RBGCSinkG(idg_N2O,0,NY,NX)=RN2O(0,NY,NX)*XNPG
  RBGCSinkG(idg_NH3,0,NY,NX)=0.0_r8
  RBGCSinkG(idg_H2,0,NY,NX)=RH2GO(0,NY,NX)*XNPG
  DO  K=1,jcplx
    ROCSK2(K,0,NY,NX)=-XOQCS(K,0,NY,NX)*XNPH
    RONSK2(K,0,NY,NX)=-XOQNS(K,0,NY,NX)*XNPH
    ROPSK2(K,0,NY,NX)=-XOQPS(K,0,NY,NX)*XNPH
    ROASK2(K,0,NY,NX)=-XOQAS(K,0,NY,NX)*XNPH
  ENDDO
  RBGCSinkS(ids_NH4,0,NY,NX)=(-XNH4S(0,NY,NX)-TRN4S(0,NY,NX))*XNPH
  RBGCSinkS(idg_NH3,0,NY,NX)=-TRN3S(0,NY,NX)*XNPH
  RBGCSinkS(ids_NO3,0,NY,NX)=(-XNO3S(0,NY,NX)-TRNO3(0,NY,NX))*XNPH
  RBGCSinkS(ids_NO2,0,NY,NX)=(-XNO2S(0,NY,NX)-TRNO2(0,NY,NX))*XNPH
  RBGCSinkS(ids_H2PO4,0,NY,NX)=(-XH2PS(0,NY,NX)-TRH2P(0,NY,NX))*XNPH
  RBGCSinkS(ids_H1PO4,0,NY,NX)=(-XH1PS(0,NY,NX)-TRH1P(0,NY,NX))*XNPH
  end subroutine SurfaceSinksandSources
!------------------------------------------------------------------------------------------

  subroutine StateVarforGasandSolute(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K,NTG,NTS

! exclude banded nutrient
  DO NTG=idg_beg,idg_end-1
    trc_solml2(NTG,0,NY,NX)=trc_solml(NTG,0,NY,NX)
  ENDDO

  D9979: DO K=1,jcplx
    OQC2(K,0,NY,NX)=OQC(K,0,NY,NX)-XOQCS(K,0,NY,NX)
    OQN2(K,0,NY,NX)=OQN(K,0,NY,NX)-XOQNS(K,0,NY,NX)
    OQP2(K,0,NY,NX)=OQP(K,0,NY,NX)-XOQPS(K,0,NY,NX)
    OQA2(K,0,NY,NX)=OQA(K,0,NY,NX)-XOQAS(K,0,NY,NX)
  ENDDO D9979

! exclude banded nutrient
  DO NTS=ids_nut_beg,ids_nuts_end
    trc_solml2(NTS,0,NY,NX)=trc_solml(NTS,0,NY,NX)
  ENDDO

  CHY0(0,NY,NX)=10.0_r8**(-(PH(0,NY,NX)-3.0_r8))
  end subroutine StateVarforGasandSolute
!------------------------------------------------------------------------------------------

  subroutine SurfaceSolutefromAtms(NY,NX)
  implicit none
  integer, intent(in) :: NY, NX
  integer :: K

  D8855: DO K=1,jcplx
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
  ENDDO D8855
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
    trcs_XFLS(idg_CO2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_CH4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_O2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_N2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_N2O,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_H2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_NH3,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NO3,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NO2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_H2PO4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_CO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_CH4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_O2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_N2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_N2O,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_H2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_NH3,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO3,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H2PO4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_NH3B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO3B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO2B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H2PO4B,3,NU(NY,NX),NY,NX)=0.0_r8
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
    trcs_XFLS(idg_CO2,3,0,NY,NX)=FLQRQ(NY,NX)*CCOR(NY,NX)+FLQRI(NY,NX)*CCOQ(NY,NX)
    trcs_XFLS(idg_CH4,3,0,NY,NX)=FLQRQ(NY,NX)*CCHR(NY,NX)+FLQRI(NY,NX)*CCHQ(NY,NX)
    trcs_XFLS(idg_O2,3,0,NY,NX)=FLQRQ(NY,NX)*COXR(NY,NX)+FLQRI(NY,NX)*COXQ(NY,NX)
    trcs_XFLS(idg_N2,3,0,NY,NX)=FLQRQ(NY,NX)*CNNR(NY,NX)+FLQRI(NY,NX)*CNNQ(NY,NX)
    trcs_XFLS(idg_N2O,3,0,NY,NX)=FLQRQ(NY,NX)*CN2R(NY,NX)+FLQRI(NY,NX)*CN2Q(NY,NX)
    trcs_XFLS(idg_H2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4,3,0,NY,NX)=(FLQRQ(NY,NX)*CN4R(NY,NX)+FLQRI(NY,NX)*CN4Q(I,NY,NX))*natomw
    trcs_XFLS(idg_NH3,3,0,NY,NX)=(FLQRQ(NY,NX)*CN3R(NY,NX)+FLQRI(NY,NX)*CN3Q(I,NY,NX))*natomw
    trcs_XFLS(ids_NO3,3,0,NY,NX)=(FLQRQ(NY,NX)*CNOR(NY,NX)+FLQRI(NY,NX)*CNOQ(I,NY,NX))*natomw
    trcs_XFLS(ids_NO2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4,3,0,NY,NX)=(FLQRQ(NY,NX)*CH1PR(NY,NX)+FLQRI(NY,NX)*CH1PQ(I,NY,NX))*patomw
    trcs_XFLS(ids_H2PO4,3,0,NY,NX)=(FLQRQ(NY,NX)*CPOR(NY,NX)+FLQRI(NY,NX)*CPOQ(I,NY,NX))*patomw
    trcs_XFLS(idg_CO2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCOR(NY,NX)+FLQGI(NY,NX)*CCOQ(NY,NX)
    trcs_XFLS(idg_CH4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCHR(NY,NX)+FLQGI(NY,NX)*CCHQ(NY,NX)
    trcs_XFLS(idg_O2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*COXR(NY,NX)+FLQGI(NY,NX)*COXQ(NY,NX)
    trcs_XFLS(idg_N2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNNR(NY,NX)+FLQGI(NY,NX)*CNNQ(NY,NX)
    trcs_XFLS(idg_N2O,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CN2R(NY,NX)+FLQGI(NY,NX)*CN2Q(NY,NX)
    trcs_XFLS(idg_H2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN4R(NY,NX) &
      +FLQGI(NY,NX)*CN4Q(I,NY,NX))*natomw)*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    trcs_XFLS(idg_NH3,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN3R(NY,NX) &
      +FLQGI(NY,NX)*CN3Q(I,NY,NX))*natomw)*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_NO3,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CNOR(NY,NX) &
      +FLQGI(NY,NX)*CNOQ(I,NY,NX))*natomw)*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_NO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CH1PR(NY,NX) &
      +FLQGI(NY,NX)*CH1PQ(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_H2PO4,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CPOR(NY,NX) &
      +FLQGI(NY,NX)*CPOQ(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_NH4B,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN4R(NY,NX) &
      +FLQGI(NY,NX)*CN4Q(I,NY,NX))*natomw)*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_XFLS(idg_NH3B,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CN3R(NY,NX) &
      +FLQGI(NY,NX)*CN3Q(I,NY,NX))*natomw)*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_NO3B,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CNOR(NY,NX) &
      +FLQGI(NY,NX)*CNOQ(I,NY,NX))*natomw)*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_NO2B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4B,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CH1PR(NY,NX) &
      +FLQGI(NY,NX)*CH1PQ(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
    trcs_XFLS(ids_H2PO4B,3,NU(NY,NX),NY,NX)=((FLQGQ(NY,NX)*CPOR(NY,NX) &
      +FLQGI(NY,NX)*CPOQ(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
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
    trcs_XFLS(idg_CO2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_CH4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_O2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_N2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_N2O,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_H2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_NH3,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NO3,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_NO2,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(ids_H2PO4,3,0,NY,NX)=0.0_r8
    trcs_XFLS(idg_CO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_CH4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_O2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_N2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_N2O,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_H2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_NH3,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO3,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H2PO4,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NH4B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(idg_NH3B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO3B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_NO2B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H1PO4B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_XFLS(ids_H2PO4B,3,NU(NY,NX),NY,NX)=0.0_r8
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
  use EcoSiMParDataMod, only : micpar
  implicit none

  integer, intent(in) :: NY, NX
  INTEGER :: K,L,NTS

  DO  K=1,micpar%n_litrsfk
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
  RCOFL0(NY,NX)=trcs_XFLS(idg_CO2,3,0,NY,NX)*XNPH
  RCHFL0(NY,NX)=trcs_XFLS(idg_CH4,3,0,NY,NX)*XNPH
  ROXFL0(NY,NX)=trcs_XFLS(idg_O2,3,0,NY,NX)*XNPH
  RNGFL0(NY,NX)=trcs_XFLS(idg_N2,3,0,NY,NX)*XNPH
  RN2FL0(NY,NX)=trcs_XFLS(idg_N2O,3,0,NY,NX)*XNPH
  RHGFL0(NY,NX)=trcs_XFLS(idg_H2,3,0,NY,NX)*XNPH
  RN4FL0(NY,NX)=trcs_XFLS(ids_NH4,3,0,NY,NX)*XNPH
  RN3FL0(NY,NX)=trcs_XFLS(idg_NH3,3,0,NY,NX)*XNPH
  RNOFL0(NY,NX)=trcs_XFLS(ids_NO3,3,0,NY,NX)*XNPH
  RNXFL0(NY,NX)=trcs_XFLS(ids_NO2,3,0,NY,NX)*XNPH
  RH1PF0(NY,NX)=trcs_XFLS(ids_H1PO4,3,0,NY,NX)*XNPH
  RH2PF0(NY,NX)=trcs_XFLS(ids_H2PO4,3,0,NY,NX)*XNPH
  RCOFL1(NY,NX)=trcs_XFLS(idg_CO2,3,NU(NY,NX),NY,NX)*XNPH
  RCHFL1(NY,NX)=trcs_XFLS(idg_CH4,3,NU(NY,NX),NY,NX)*XNPH
  ROXFL1(NY,NX)=trcs_XFLS(idg_O2,3,NU(NY,NX),NY,NX)*XNPH
  RNGFL1(NY,NX)=trcs_XFLS(idg_N2,3,NU(NY,NX),NY,NX)*XNPH
  RN2FL1(NY,NX)=trcs_XFLS(idg_N2O,3,NU(NY,NX),NY,NX)*XNPH
  RHGFL1(NY,NX)=trcs_XFLS(idg_H2,3,NU(NY,NX),NY,NX)*XNPH
  RN4FL1(NY,NX)=trcs_XFLS(ids_NH4,3,NU(NY,NX),NY,NX)*XNPH
  RN3FL1(NY,NX)=trcs_XFLS(idg_NH3,3,NU(NY,NX),NY,NX)*XNPH
  RNOFL1(NY,NX)=trcs_XFLS(ids_NO3,3,NU(NY,NX),NY,NX)*XNPH
  RNXFL1(NY,NX)=trcs_XFLS(ids_NO2,3,NU(NY,NX),NY,NX)*XNPH
  RH1PF1(NY,NX)=trcs_XFLS(ids_H1PO4,3,NU(NY,NX),NY,NX)*XNPH
  RH2PF1(NY,NX)=trcs_XFLS(ids_H2PO4,3,NU(NY,NX),NY,NX)*XNPH
  RN4FL2(NY,NX)=trcs_XFLS(ids_NH4B,3,NU(NY,NX),NY,NX)*XNPH
  RN3FL2(NY,NX)=trcs_XFLS(idg_NH3B,3,NU(NY,NX),NY,NX)*XNPH
  RNOFL2(NY,NX)=trcs_XFLS(ids_NO3B,3,NU(NY,NX),NY,NX)*XNPH
  RNXFL2(NY,NX)=trcs_XFLS(ids_NO2B,3,NU(NY,NX),NY,NX)*XNPH
  RH1BF2(NY,NX)=trcs_XFLS(ids_H1PO4B,3,NU(NY,NX),NY,NX)*XNPH
  RH2BF2(NY,NX)=trcs_XFLS(ids_H2PO4B,3,NU(NY,NX),NY,NX)*XNPH
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
  DO NTS=ids_beg,ids_end
    SolDifcc(NTS,0,NY,NX)=SolDifc(NTS,0,NY,NX)*XNPH
  ENDDO

  OCSGL2(0,NY,NX)=OCSGL(0,NY,NX)*XNPH
  ONSGL2(0,NY,NX)=ONSGL(0,NY,NX)*XNPH
  OPSGL2(0,NY,NX)=OPSGL(0,NY,NX)*XNPH
  OASGL2(0,NY,NX)=OASGL(0,NY,NX)*XNPH
!
!     INITIAL SOLUTES IN SNOWPACK
!
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  DO L=1,JS

    trcg_solsml2(idg_beg:idg_end-1,L,NY,NX)=trcg_solsml(idg_beg:idg_end-1,L,NY,NX)

    trcn_solsml2(ids_nut_beg:ids_nuts_end,L,NY,NX)=trcn_solsml(ids_nut_beg:ids_nuts_end,L,NY,NX)
  enddo
  end subroutine SubHourlyFluxesFromSiteFile
!------------------------------------------------------------------------------------------

  subroutine ImportFluxFromOutsideModules(I,NY,NX)
  implicit none

  integer, intent(in) ::  I,NY, NX

  real(r8) :: FLWU(JZ,JY,JX)
  integer :: L,K,NTG,NTS

  DO L=NU(NY,NX),NL(NY,NX)
    CHY0(L,NY,NX)=10.0_r8**(-(PH(L,NY,NX)-3.0_r8))
    FLWU(L,NY,NX)=TUPWTR(L,NY,NX)*XNPH
    RBGCSinkG(idg_CO2,L,NY,NX)=(RCO2O(L,NY,NX)+TCO2S(L,NY,NX)-TRCO2(L,NY,NX))*XNPG
    RBGCSinkG(idg_CH4,L,NY,NX)=(RCH4O(L,NY,NX)+TUPCHS(L,NY,NX))*XNPG
    RBGCSinkG(idg_N2,L,NY,NX)=(RN2G(L,NY,NX)+XN2GS(L,NY,NX)+TUPNF(L,NY,NX))*XNPG
    RBGCSinkG(idg_N2O,L,NY,NX)=(RN2O(L,NY,NX)+TUPN2S(L,NY,NX))*XNPG
    RBGCSinkG(idg_NH3,L,NY,NX)=-TRN3G(L,NY,NX)*XNPG
    RBGCSinkG(idg_H2,L,NY,NX)=(RH2GO(L,NY,NX)+TUPHGS(L,NY,NX))*XNPG
    DO  K=1,jcplx
      ROCSK2(K,L,NY,NX)=-XOQCS(K,L,NY,NX)*XNPH
      RONSK2(K,L,NY,NX)=-XOQNS(K,L,NY,NX)*XNPH
      ROPSK2(K,L,NY,NX)=-XOQPS(K,L,NY,NX)*XNPH
      ROASK2(K,L,NY,NX)=-XOQAS(K,L,NY,NX)*XNPH
    ENDDO
    RBGCSinkS(ids_NH4,L,NY,NX)=(-XNH4S(L,NY,NX)-TRN4S(L,NY,NX)+TUPNH4(L,NY,NX))*XNPH
    RBGCSinkS(idg_NH3,L,NY,NX)=(-TRN3S(L,NY,NX)+TUPN3S(L,NY,NX))*XNPH
    RBGCSinkS(ids_NO3,L,NY,NX)=(-XNO3S(L,NY,NX)-TRNO3(L,NY,NX)+TUPNO3(L,NY,NX))*XNPH
    RBGCSinkS(ids_NO2,L,NY,NX)=(-XNO2S(L,NY,NX)-TRNO2(L,NY,NX))*XNPH
    RBGCSinkS(ids_H2PO4,L,NY,NX)=(-XH2PS(L,NY,NX)-TRH2P(L,NY,NX)+TUPH2P(L,NY,NX))*XNPH
    RBGCSinkS(ids_H1PO4,L,NY,NX)=(-XH1PS(L,NY,NX)-TRH1P(L,NY,NX)+TUPH1P(L,NY,NX))*XNPH
    RBGCSinkS(ids_NH4B,L,NY,NX)=(-XNH4B(L,NY,NX)-TRN4B(L,NY,NX)+TUPNHB(L,NY,NX))*XNPH
    RBGCSinkS(idg_NH3B,L,NY,NX)=(-TRN3B(L,NY,NX)+TUPN3B(L,NY,NX))*XNPH
    RBGCSinkS(ids_NO3B,L,NY,NX)=(-XNO3B(L,NY,NX)-TRNOB(L,NY,NX)+TUPNOB(L,NY,NX))*XNPH
    RBGCSinkS(ids_NO2B,L,NY,NX)=(-XNO2B(L,NY,NX)-TRN2B(L,NY,NX))*XNPH
    RBGCSinkS(ids_H2PO4B,L,NY,NX)=(-XH2BS(L,NY,NX)-TRH2B(L,NY,NX)+TUPH2B(L,NY,NX))*XNPH
    RBGCSinkS(ids_H1PO4B,L,NY,NX)=(-XH1BS(L,NY,NX)-TRH1B(L,NY,NX)+TUPH1B(L,NY,NX))*XNPH
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
    RN4FLU(L,NY,NX)=FLU(L,NY,NX)*CN4Q(I,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    RN3FLU(L,NY,NX)=FLU(L,NY,NX)*CN3Q(I,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    RNOFLU(L,NY,NX)=FLU(L,NY,NX)*CNOQ(I,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)*natomw
    RH1PFU(L,NY,NX)=FLU(L,NY,NX)*CH1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    RH2PFU(L,NY,NX)=FLU(L,NY,NX)*CPOQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    RN4FBU(L,NY,NX)=FLU(L,NY,NX)*CN4Q(I,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    RN3FBU(L,NY,NX)=FLU(L,NY,NX)*CN3Q(I,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    RNOFBU(L,NY,NX)=FLU(L,NY,NX)*CNOQ(I,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)*natomw
    RH1BBU(L,NY,NX)=FLU(L,NY,NX)*CH1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw
    RH2BBU(L,NY,NX)=FLU(L,NY,NX)*CPOQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw
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

    DO NTS=ids_beg,ids_end
      SolDifcc(NTS,L,NY,NX)=SolDifc(NTS,L,NY,NX)*XNPH
    ENDDO

    DO NTG=idg_beg,idg_end
      GasDifcc(NTG,L,NY,NX)=GasDifc(NTG,L,NY,NX)*XNPG
    ENDDO
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
! exclude NH3B
    DO NTG=idg_beg,idg_end-1
      trc_gasml2(NTG,L,NY,NX)=trc_gasml(NTG,L,NY,NX)
    ENDDO

    DO  K=1,jcplx
      OQC2(K,L,NY,NX)=OQC(K,L,NY,NX)-XOQCS(K,L,NY,NX)
      OQN2(K,L,NY,NX)=OQN(K,L,NY,NX)-XOQNS(K,L,NY,NX)
      OQP2(K,L,NY,NX)=OQP(K,L,NY,NX)-XOQPS(K,L,NY,NX)
      OQA2(K,L,NY,NX)=OQA(K,L,NY,NX)-XOQAS(K,L,NY,NX)
      OQCH2(K,L,NY,NX)=OQCH(K,L,NY,NX)
      OQNH2(K,L,NY,NX)=OQNH(K,L,NY,NX)
      OQPH2(K,L,NY,NX)=OQPH(K,L,NY,NX)
      OQAH2(K,L,NY,NX)=OQAH(K,L,NY,NX)
    enddo
    DO NTS=ids_beg,ids_end
      trc_solml2(NTS,L,NY,NX)=trc_solml(NTS,L,NY,NX)
      trc_soHml2(NTS,L,NY,NX)=trc_soHml(NTS,L,NY,NX)
    ENDDO

  enddo
  end subroutine ImportFluxFromOutsideModules
!

end module TrnsfrMod
