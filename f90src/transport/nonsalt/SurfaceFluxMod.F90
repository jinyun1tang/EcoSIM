module SurfaceFluxMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridDataType
  use TransfrDataMod
  use SoilBGCDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SurfLitterDataType
  use AqueChemDatatype
  use SnowDataType
  use EcoSimConst
  use TracerIDMod
  use minimathmod, only : AZMAX1,AZMIN1
  use ChemTranspDataType
  use ClimForcDataType
  use EcoSIMSolverPar
  use LandSurfDataType
  use SurfSoilDataType
  use SoilPropertyDataType
implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__
  real(r8) :: RCODXR,RCHDXR,ROXDXR,RNGDXR,RN2DXR,RN3DXR,RHGDXR

  real(r8) :: RCOFLS1,RCHFLS1,RNGFLS1,RN2FLS1,ROXFLS1
  real(r8) :: RN4FLW1,RN3FLW1,RNOFLW1
  real(r8) :: RN4FLB1,RN3FLB1,RNOFLB1,RH1BFB1,RH2BFB1
  real(r8) :: RH1PFS1,RH2PFS1

  public :: SoluteFluxSurface
  public :: LitterGasVolatilDissol
  public :: SurfSoilFluxGasDifAdv
  public :: SoluteFluxSnowpackDisch
contains


!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSurface(M,NY,NX,NHE,NHW,NVS,NVN,&
    FLQM,trcg_RFLS0,trcn_RFLW0,RDXS_gas)
  implicit none
  integer, intent(in) :: NY,NX,M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8),intent(in) :: trcg_RFLS0(idg_beg:idg_end-1)
  real(r8),intent(in) :: trcn_RFLW0(ids_nut_beg:ids_nuts_end)
  real(r8),intent(inout) :: FLQM(3,JD,JV,JH)
  real(r8),intent(out) :: RDXS_gas(idg_beg:idg_end)
  real(r8) :: FLWRM1
  real(r8) :: trcs_cl1(ids_beg:ids_end)
  real(r8) :: trcs_cl2(ids_beg:ids_end)
  real(r8) :: RFLs_adv(ids_beg:ids_end)
  real(r8) :: SDifFlx(ids_beg:ids_end)
!     VLWatMicPM,VLWatMacPM,VLsoiAirPM,FLPM=micropore,macropore water volume, air volume and change in air volume
!     FLWM,WaterFlowMacPi=water flux into soil micropore,macropore from watsub.f
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     FLVM,FLM=air,water flux in gas flux calculations
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
  VLWatMicPMA(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
  VLWatMicPMB(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
  VLWatMicPXA(NU(NY,NX),NY,NX)=natomw*VLWatMicPMA(NU(NY,NX),NY,NX)
  VLWatMicPXB(NU(NY,NX),NY,NX)=natomw*VLWatMicPMB(NU(NY,NX),NY,NX)
  VLsoiAirPMA(NU(NY,NX),NY,NX)=VLsoiAirPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
  VLsoiAirPMB(NU(NY,NX),NY,NX)=VLsoiAirPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
  FLVM(NU(NY,NX),NY,NX)=FLPM(M,NU(NY,NX),NY,NX)*XNPT
  FLQM(3,NU(NY,NX),NY,NX)=(FLWM(M,3,NU(NY,NX),NY,NX)+WaterFlowMacPi(M,3,NU(NY,NX),NY,NX))*XNPT
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE RESIDUE
!
!     VOLT,DLYR,AREA=litter volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     TORT=tortuosity from hour1.f
!     *SGL*=solute diffusivity
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     DFGS*=effective solute diffusivity
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in litter
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in litter
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in litter
!     C*1=solute concentration in litter
!
  call LitterAtmosExchange(M,NY,NX,trcs_cl1)

  call SoilAtmosExchange(M,NY,NX,trcs_cl2,RDXS_gas)

!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
  call ConvectiveSurfaceSoluteFlux(M,NY,NX,FLWRM1,RFLs_adv)
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  call LitterSurfSoilExchange(M,NY,NX,FLWRM1,trcs_cl1,trcs_cl2,SDifFlx)

!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES
!
  call TotalPoreFluxAdjacentCell(NY,NX,SDifFlx,RFLs_adv,trcg_RFLS0,trcn_RFLW0)

!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     SURFACE LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
  call MacMicPoresTransfer(M,NY,NX)

!
!     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
!     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE LAYER
!
  call OverlandFlowSnowdriftTransport(M,NY,NX,NHE,NHW,NVS,NVN)
  end subroutine SoluteFluxSurface

!------------------------------------------------------------------------------------------

  subroutine ConvectiveSurfaceSoluteFlux(M,NY,NX,FLWRM1,RFLs_adv)
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8),intent(out) :: RFLs_adv(ids_beg:ids_end)
  real(r8),intent(out) :: FLWRM1
  REAL(R8) :: VFLW
  integer :: K,NTS

  FLWRM1=FLWRM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!     VLWatMicPM=litter water volume
!     RFL*=soil-litter convective solute flux
!     *S2=litter solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  IF(FLWRM1.GT.0.0_r8)THEN
    IF(VLWatMicPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VLWatMicPM(M,0,NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,0,NY,NX))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,0,NY,NX))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,0,NY,NX))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,0,NY,NX))
    ENDDO

    DO NTS=ids_beg,ids_end
      RFLs_adv(NTS)=VFLW*AZMAX1(trc_solml2(NTS,0,NY,NX))*trcs_VLN(NTS,NU(NY,NX),NY,NX)
    ENDDO
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     RFL*=soil-litter convective solute flux
!     *S2=soil solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  ELSE
    IF(VLWatMicPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,NU(NY,NX),NY,NX))
    ENDDO

    DO NTS=ids_beg,ids_end
      RFLs_adv(NTS)=VFLW*AZMAX1(trc_solml2(NTS,NU(NY,NX),NY,NX))
    ENDDO
  ENDIF
  end subroutine ConvectiveSurfaceSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine LitterSurfSoilExchange(M,NY,NX,FLWRM1,trcs_cl1,trcs_cl2,SDifFlx)
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8), intent(in) :: FLWRM1
  real(r8), intent(in) :: trcs_cl1(ids_beg:ids_end)
  real(r8), intent(in) :: trcs_cl2(ids_beg:ids_end)
  real(r8), intent(out):: SDifFlx(ids_beg:ids_end)

  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: DIFOC,DIFON,DIFOP,DIFOA

  real(r8) :: DIFOC0,DIFON0,DIFOP0,DIFOA0
  real(r8) :: DIFOC1,DIFON1,DIFOP1,DIFOA1
  real(r8) :: DISPN,DIF0,DIF1
  real(r8) :: SDifc(ids_beg:ids_end)
  integer  :: K,NTS,NTG,NTN
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!
  IF((VGeomLayer(0,NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPM(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND.(VLWatMicPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
!
!     DLYR0,DLYR1=litter, soil surface thickness
!     TORT=tortuosity from hour1.f
!     CVRD,BARE=litter cover fraction,1-CVRD
!     DISP=dispersivity parameter
!     FLWRM=litter-soil water flux from watsub.f
!     DIF*0,DIF*1=aqueous diffusivity-dispersivity in litter, soil surface
!     *SGL*=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!
    DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX))
    TORT0=TORT(M,0,NY,NX)/DLYR0*CVRD(NY,NX)
    DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
    TORT1=TORT(M,NU(NY,NX),NY,NX)/DLYR1
    DISPN=DISP(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))
    DIFOC0=(OCSGL2(0,NY,NX)*TORT0+DISPN)
    DIFON0=(ONSGL2(0,NY,NX)*TORT0+DISPN)
    DIFOP0=(OPSGL2(0,NY,NX)*TORT0+DISPN)
    DIFOA0=(OASGL2(0,NY,NX)*TORT0+DISPN)

    DIFOC1=(OCSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFON1=(ONSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFOP1=(OPSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFOA1=(OASGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)

    DIFOC=DIFOC0*DIFOC1/(DIFOC0+DIFOC1)*AREA(3,NU(NY,NX),NY,NX)
    DIFON=DIFON0*DIFON1/(DIFON0+DIFON1)*AREA(3,NU(NY,NX),NY,NX)
    DIFOP=DIFOP0*DIFOP1/(DIFOP0+DIFOP1)*AREA(3,NU(NY,NX),NY,NX)
    DIFOA=DIFOA0*DIFOA1/(DIFOA0+DIFOA1)*AREA(3,NU(NY,NX),NY,NX)

    DO NTS=ids_beg,ids_end
      DIF0=(SolDifcc(NTS,0,NY,NX)*TORT0+DISPN)
      DIF1=(SolDifcc(NTS,NU(NY,NX),NY,NX)*TORT1+DISPN)
      SDifc(NTS)=DIF0*DIF1/(DIF0+DIF1)*AREA(3,NU(NY,NX),NY,NX)
    ENDDO

!     DFV*S,DFV*B=diffusive solute flux between litter and soil surface in non-band,band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     C*1,C*2=solute concentration in litter,soil surface
!
    DO  K=1,jcplx
      DFVOC(K)=DIFOC*(COQC1(K)-COQC2(K))
      DFVON(K)=DIFON*(COQN1(K)-COQN2(K))
      DFVOP(K)=DIFOP*(COQP1(K)-COQP2(K))
      DFVOA(K)=DIFOA*(COQA1(K)-COQA2(K))
    ENDDO
!exclude NH3B and NH3
    DO NTG=idg_beg,idg_end-2
        SDifFlx(NTG)=SDifc(NTG)*(trcs_cl1(NTG)-trcs_cl2(NTG))
    ENDDO

    DO NTN=ids_nuts_beg,ids_nuts_end
      SDifFlx(NTN)=SDifc(NTN)*(trcs_cl1(NTN)-trcs_cl2(NTN))*trcs_VLN(NTN,NU(NY,NX),NY,NX)
    ENDDO
  ELSE
    DO  K=1,jcplx
      DFVOC(K)=0.0_r8
      DFVON(K)=0.0_r8
      DFVOP(K)=0.0_r8
      DFVOA(K)=0.0_r8
    ENDDO
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF
  end subroutine LitterSurfSoilExchange


!------------------------------------------------------------------------------------------

  subroutine TotalPoreFluxAdjacentCell(NY,NX,SDifFlx,RFLs_adv,trcg_RFLS0,trcn_RFLW0)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: NY, NX
  real(r8), intent(in) :: SDifFlx(ids_beg:ids_end)
  real(r8), intent(in) :: RFLs_adv(ids_beg:ids_end)
  real(r8), intent(in) :: trcg_RFLS0(idg_beg:idg_end-1)
  real(r8), intent(in) :: trcn_RFLW0(ids_nut_beg:ids_nuts_end)
  integer :: K,NTG

!     R*FLS=convective + diffusive solute flux between litter, soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into soil in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     R*FL0,R*FL1=convective flux into surface litter, soil surface
!     RFL*=convective flux between surface litter and soil surface
!     DFV*=diffusive solute flux between litter and soil surface
!
  DO K=1,micpar%n_litrsfk
    ROCFLS(K,3,0,NY,NX)=ROCFL0(K,NY,NX)-RFLOC(K)-DFVOC(K)
    RONFLS(K,3,0,NY,NX)=RONFL0(K,NY,NX)-RFLON(K)-DFVON(K)
    ROPFLS(K,3,0,NY,NX)=ROPFL0(K,NY,NX)-RFLOP(K)-DFVOP(K)
    ROAFLS(K,3,0,NY,NX)=ROAFL0(K,NY,NX)-RFLOA(K)-DFVOA(K)
    ROCFLS(K,3,NU(NY,NX),NY,NX)=ROCFL1(K,NY,NX)+RFLOC(K)+DFVOC(K)
    RONFLS(K,3,NU(NY,NX),NY,NX)=RONFL1(K,NY,NX)+RFLON(K)+DFVON(K)
    ROPFLS(K,3,NU(NY,NX),NY,NX)=ROPFL1(K,NY,NX)+RFLOP(K)+DFVOP(K)
    ROAFLS(K,3,NU(NY,NX),NY,NX)=ROAFL1(K,NY,NX)+RFLOA(K)+DFVOA(K)
  ENDDO

  DO NTG=idg_beg,idg_end-1
    R3PoreSolFlx(NTG,3,0,NY,NX)=trcg_RFL0(NTG,NY,NX)+trcg_RFLS0(NTG)-RFLs_adv(NTG)-SDifFlx(NTG)
  ENDDO

  R3PoreSolFlx(idg_NH3,3,0,NY,NX)=R3PoreSolFlx(idg_NH3,3,0,NY,NX)-RFLs_adv(idg_NH3B)-SDifFlx(idg_NH3B)

  R3PoreSolFlx(ids_NH4,3,0,NY,NX)=trcn_RFL0(ids_NH4,NY,NX)+trcn_RFLW0(ids_NH4)-RFLs_adv(ids_NH4)-SDifFlx(ids_NH4)-RFLs_adv(ids_NH4B)-SDifFlx(ids_NH4B)
  R3PoreSolFlx(ids_NO3,3,0,NY,NX)=trcn_RFL0(ids_NO3,NY,NX)+trcn_RFLW0(ids_NO3)-RFLs_adv(ids_NO3)-SDifFlx(ids_NO3)-RFLs_adv(ids_NO3B)-SDifFlx(ids_NO3B)
  R3PoreSolFlx(ids_NO2,3,0,NY,NX)=trcn_RFL0(ids_NO2,NY,NX)+trcn_RFLW0(ids_NO2)-RFLs_adv(ids_NO2)-SDifFlx(ids_NO2)-RFLs_adv(ids_NO2B)-SDifFlx(ids_NO2B)
  R3PoreSolFlx(ids_H1PO4,3,0,NY,NX)=trcn_RFL0(ids_H1PO4,NY,NX)+trcn_RFLW0(ids_H1PO4)-RFLs_adv(ids_H1PO4)-SDifFlx(ids_H1PO4)-RFLs_adv(ids_H1PO4B)-SDifFlx(ids_H1PO4B)
  R3PoreSolFlx(ids_H2PO4,3,0,NY,NX)=trcn_RFL0(ids_H2PO4,NY,NX)+trcn_RFLW0(ids_H2PO4)-RFLs_adv(ids_H2PO4)-SDifFlx(ids_H2PO4)-RFLs_adv(ids_H2PO4B)-SDifFlx(ids_H2PO4B)

  R3PoreSolFlx(idg_CO2,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_CO2,NY,NX)+RCOFLS1+RFLs_adv(idg_CO2)+SDifFlx(idg_CO2)
  R3PoreSolFlx(idg_CH4,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_CH4,NY,NX)+RCHFLS1+RFLs_adv(idg_CH4)+SDifFlx(idg_CH4)
  R3PoreSolFlx(idg_O2,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_O2,NY,NX)+ROXFLS1+RFLs_adv(idg_O2)+SDifFlx(idg_O2)
  R3PoreSolFlx(idg_N2,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_N2,NY,NX)+RNGFLS1+RFLs_adv(idg_N2)+SDifFlx(idg_N2)
  R3PoreSolFlx(idg_N2O,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_N2O,NY,NX)+RN2FLS1+RFLs_adv(idg_N2O)+SDifFlx(idg_N2O)
  R3PoreSolFlx(idg_H2,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_H2,NY,NX)+RFLs_adv(idg_H2)+SDifFlx(idg_H2)
  R3PoreSolFlx(idg_NH3,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_NH3,NY,NX)+RN3FLW1+RFLs_adv(idg_NH3)+SDifFlx(idg_NH3)

  R3PoreSolFlx(ids_NH4,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_NH4,NY,NX)+RN4FLW1+RFLs_adv(ids_NH4)+SDifFlx(ids_NH4)
  R3PoreSolFlx(ids_NO3,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_NO3,NY,NX)+RNOFLW1+RFLs_adv(ids_NO3)+SDifFlx(ids_NO3)
  R3PoreSolFlx(ids_NO2,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_NO2,NY,NX)+RFLs_adv(ids_NO2)+SDifFlx(ids_NO2)
  R3PoreSolFlx(ids_H1PO4,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_H1PO4,NY,NX)+RH1PFS1+RFLs_adv(ids_H1PO4)+SDifFlx(ids_H1PO4)
  R3PoreSolFlx(ids_H2PO4,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_H2PO4,NY,NX)+RH2PFS1+RFLs_adv(ids_H2PO4)+SDifFlx(ids_H2PO4)

  R3PoreSolFlx(ids_NH4B,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_NH4B,NY,NX)+RN4FLB1+RFLs_adv(ids_NH4B)+SDifFlx(ids_NH4B)
  R3PoreSolFlx(idg_NH3B,3,NU(NY,NX),NY,NX)=trcs_RFL1(idg_NH3B,NY,NX)+RN3FLB1+RFLs_adv(idg_NH3B)+SDifFlx(idg_NH3B)
  R3PoreSolFlx(ids_NO3B,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_NO3B,NY,NX)+RNOFLB1+RFLs_adv(ids_NO3B)+SDifFlx(ids_NO3B)
  R3PoreSolFlx(ids_NO2B,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_NO2B,NY,NX)+RFLs_adv(ids_NO2B)+SDifFlx(ids_NO2B)
  R3PoreSolFlx(ids_H1PO4B,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_H1PO4B,NY,NX)+RH1BFB1+RFLs_adv(ids_H1PO4B)+SDifFlx(ids_H1PO4B)
  R3PoreSolFlx(ids_H2PO4B,3,NU(NY,NX),NY,NX)=trcs_RFL1(ids_H2PO4B,NY,NX)+RH2BFB1+RFLs_adv(ids_H2PO4B)+SDifFlx(ids_H2PO4B)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  DO K=1,micpar%n_litrsfk
    XOCFLS(K,3,0,NY,NX)=XOCFLS(K,3,0,NY,NX)-RFLOC(K)-DFVOC(K)
    XONFLS(K,3,0,NY,NX)=XONFLS(K,3,0,NY,NX)-RFLON(K)-DFVON(K)
    XOPFLS(K,3,0,NY,NX)=XOPFLS(K,3,0,NY,NX)-RFLOP(K)-DFVOP(K)
    XOAFLS(K,3,0,NY,NX)=XOAFLS(K,3,0,NY,NX)-RFLOA(K)-DFVOA(K)
    XOCFLS(K,3,NU(NY,NX),NY,NX)=XOCFLS(K,3,NU(NY,NX),NY,NX)+RFLOC(K)+DFVOC(K)
    XONFLS(K,3,NU(NY,NX),NY,NX)=XONFLS(K,3,NU(NY,NX),NY,NX)+RFLON(K)+DFVON(K)
    XOPFLS(K,3,NU(NY,NX),NY,NX)=XOPFLS(K,3,NU(NY,NX),NY,NX)+RFLOP(K)+DFVOP(K)
    XOAFLS(K,3,NU(NY,NX),NY,NX)=XOAFLS(K,3,NU(NY,NX),NY,NX)+RFLOA(K)+DFVOA(K)
  ENDDO
  trcs_XFLS(idg_CO2,3,0,NY,NX)=trcs_XFLS(idg_CO2,3,0,NY,NX)+trcg_RFLS0(idg_CO2)-RFLs_adv(idg_CO2)-SDifFlx(idg_CO2)
  trcs_XFLS(idg_CH4,3,0,NY,NX)=trcs_XFLS(idg_CH4,3,0,NY,NX)+trcg_RFLS0(idg_CH4)-RFLs_adv(idg_CH4)-SDifFlx(idg_CH4)
  trcs_XFLS(idg_O2,3,0,NY,NX)=trcs_XFLS(idg_O2,3,0,NY,NX)+trcg_RFLS0(idg_O2)-RFLs_adv(idg_O2)-SDifFlx(idg_O2)
  trcs_XFLS(idg_N2,3,0,NY,NX)=trcs_XFLS(idg_N2,3,0,NY,NX)+trcg_RFLS0(idg_N2)-RFLs_adv(idg_N2)-SDifFlx(idg_N2)
  trcs_XFLS(idg_N2O,3,0,NY,NX)=trcs_XFLS(idg_N2O,3,0,NY,NX)+trcg_RFLS0(idg_N2O)-RFLs_adv(idg_N2O)-SDifFlx(idg_N2O)
  trcs_XFLS(idg_H2,3,0,NY,NX)=trcs_XFLS(idg_H2,3,0,NY,NX)-RFLs_adv(idg_H2)-SDifFlx(idg_H2)
  trcs_XFLS(ids_NH4,3,0,NY,NX)=trcs_XFLS(ids_NH4,3,0,NY,NX)+trcn_RFLW0(ids_NH4)-RFLs_adv(ids_NH4)-SDifFlx(ids_NH4)-RFLs_adv(ids_NH4B)-SDifFlx(ids_NH4B)
  trcs_XFLS(idg_NH3,3,0,NY,NX)=trcs_XFLS(idg_NH3,3,0,NY,NX)+trcg_RFLS0(idg_NH3)-RFLs_adv(idg_NH3)-SDifFlx(idg_NH3)-RFLs_adv(idg_NH3B)-SDifFlx(idg_NH3B)
  trcs_XFLS(ids_NO3,3,0,NY,NX)=trcs_XFLS(ids_NO3,3,0,NY,NX)+trcn_RFLW0(ids_NO3)-RFLs_adv(ids_NO3)-SDifFlx(ids_NO3)-RFLs_adv(ids_NO3B)-SDifFlx(ids_NO3B)
  trcs_XFLS(ids_NO2,3,0,NY,NX)=trcs_XFLS(ids_NO2,3,0,NY,NX)-RFLs_adv(ids_NO2)-SDifFlx(ids_NO2)-RFLs_adv(ids_NO2B)-SDifFlx(ids_NO2B)
  trcs_XFLS(ids_H1PO4,3,0,NY,NX)=trcs_XFLS(ids_H1PO4,3,0,NY,NX)+trcn_RFLW0(ids_H1PO4)-RFLs_adv(ids_H1PO4)-SDifFlx(ids_H1PO4)-RFLs_adv(ids_H1PO4B)-SDifFlx(ids_H1PO4B)
  trcs_XFLS(ids_H2PO4,3,0,NY,NX)=trcs_XFLS(ids_H2PO4,3,0,NY,NX)+trcn_RFLW0(ids_H2PO4)-RFLs_adv(ids_H2PO4)-SDifFlx(ids_H2PO4)-RFLs_adv(ids_H2PO4B)-SDifFlx(ids_H2PO4B)
  trcs_XFLS(idg_CO2,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_CO2,3,NU(NY,NX),NY,NX)+RCOFLS1+RFLs_adv(idg_CO2)+SDifFlx(idg_CO2)
  trcs_XFLS(idg_CH4,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_CH4,3,NU(NY,NX),NY,NX)+RCHFLS1+RFLs_adv(idg_CH4)+SDifFlx(idg_CH4)
  trcs_XFLS(idg_O2,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_O2,3,NU(NY,NX),NY,NX)+ROXFLS1+RFLs_adv(idg_O2)+SDifFlx(idg_O2)
  trcs_XFLS(idg_N2,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_N2,3,NU(NY,NX),NY,NX)+RNGFLS1+RFLs_adv(idg_N2)+SDifFlx(idg_N2)
  trcs_XFLS(idg_N2O,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_N2O,3,NU(NY,NX),NY,NX)+RN2FLS1+RFLs_adv(idg_N2O)+SDifFlx(idg_N2O)
  trcs_XFLS(idg_H2,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_H2,3,NU(NY,NX),NY,NX)+RFLs_adv(idg_H2)+SDifFlx(idg_H2)
  trcs_XFLS(ids_NH4,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_NH4,3,NU(NY,NX),NY,NX)+RN4FLW1+RFLs_adv(ids_NH4)+SDifFlx(ids_NH4)
  trcs_XFLS(idg_NH3,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_NH3,3,NU(NY,NX),NY,NX)+RN3FLW1+RFLs_adv(idg_NH3)+SDifFlx(idg_NH3)
  trcs_XFLS(ids_NO3,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_NO3,3,NU(NY,NX),NY,NX)+RNOFLW1+RFLs_adv(ids_NO3)+SDifFlx(ids_NO3)
  trcs_XFLS(ids_NO2,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_NO2,3,NU(NY,NX),NY,NX)+RFLs_adv(ids_NO2)+SDifFlx(ids_NO2)
  trcs_XFLS(ids_H1PO4,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_H1PO4,3,NU(NY,NX),NY,NX)+RH1PFS1+RFLs_adv(ids_H1PO4)+SDifFlx(ids_H1PO4)
  trcs_XFLS(ids_H2PO4,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_H2PO4,3,NU(NY,NX),NY,NX)+RH2PFS1+RFLs_adv(ids_H2PO4)+SDifFlx(ids_H2PO4)
  trcs_XFLS(ids_NH4B,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_NH4B,3,NU(NY,NX),NY,NX)+RN4FLB1+RFLs_adv(ids_NH4B)+SDifFlx(ids_NH4B)
  trcs_XFLS(idg_NH3B,3,NU(NY,NX),NY,NX)=trcs_XFLS(idg_NH3B,3,NU(NY,NX),NY,NX)+RN3FLB1+RFLs_adv(idg_NH3B)+SDifFlx(idg_NH3B)
  trcs_XFLS(ids_NO3B,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_NO3B,3,NU(NY,NX),NY,NX)+RNOFLB1+RFLs_adv(ids_NO3B)+SDifFlx(ids_NO3B)
  trcs_XFLS(ids_NO2B,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_NO2B,3,NU(NY,NX),NY,NX)+RFLs_adv(ids_NO2B)+SDifFlx(ids_NO2B)
  trcs_XFLS(ids_H1PO4B,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_H1PO4B,3,NU(NY,NX),NY,NX)+RH1BFB1+RFLs_adv(ids_H1PO4B)+SDifFlx(ids_H1PO4B)
  trcs_XFLS(ids_H2PO4B,3,NU(NY,NX),NY,NX)=trcs_XFLS(ids_H2PO4B,3,NU(NY,NX),NY,NX)+RH2BFB1+RFLs_adv(ids_H2PO4B)+SDifFlx(ids_H2PO4B)
  end subroutine TotalPoreFluxAdjacentCell
!------------------------------------------------------------------------------------------

  subroutine MacMicPoresTransfer(M,NY,NX)
  implicit none

  integer, intent(in) :: M, NY, NX
  integer :: K,NTS
  real(r8) :: VFLW
  real(r8) :: VLWatMacPS,VOLWT
  real(r8) :: SDifFlx(ids_beg:ids_end)
  real(r8) :: SAdvFlx(ids_beg:ids_end)
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VLWatMicPM,VLWatMacPM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FINHM(M,NU(NY,NX),NY,NX).GT.0.0_r8)THEN
    IF(VLWatMacPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,NU(NY,NX),NY,NX) &
      /VLWatMacPM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQCH2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AZMAX1(OQNH2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AZMAX1(OQPH2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AZMAX1(OQAH2(K,NU(NY,NX),NY,NX))
    ENDDO
    SAdvFlx(idg_CO2)=VFLW*AZMAX1(trc_soHml2(idg_CO2,NU(NY,NX),NY,NX))
    SAdvFlx(idg_CH4)=VFLW*AZMAX1(trc_soHml2(idg_CH4,NU(NY,NX),NY,NX))
    SAdvFlx(idg_O2)=VFLW*AZMAX1(trc_soHml2(idg_O2,NU(NY,NX),NY,NX))
    SAdvFlx(idg_N2)=VFLW*AZMAX1(trc_soHml2(idg_N2,NU(NY,NX),NY,NX))
    SAdvFlx(idg_N2O)=VFLW*AZMAX1(trc_soHml2(idg_N2O,NU(NY,NX),NY,NX))
    SAdvFlx(idg_H2)=VFLW*AZMAX1(trc_soHml2(idg_H2,NU(NY,NX),NY,NX))
    SAdvFlx(ids_NH4)=VFLW*AZMAX1(trc_soHml2(ids_NH4,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3)=VFLW*AZMAX1(trc_soHml2(idg_NH3,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3)=VFLW*AZMAX1(trc_soHml2(ids_NO3,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2)=VFLW*AZMAX1(trc_soHml2(ids_NO2,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4)=VFLW*AZMAX1(trc_soHml2(ids_H1PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4)=VFLW*AZMAX1(trc_soHml2(ids_H2PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NH4B)=VFLW*AZMAX1(trc_soHml2(ids_NH4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3B)=VFLW*AZMAX1(trc_soHml2(idg_NH3B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3B)=VFLW*AZMAX1(trc_soHml2(ids_NO3B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2B)=VFLW*AZMAX1(trc_soHml2(ids_NO2B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4B)=VFLW*AZMAX1(trc_soHml2(ids_H1PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4B)=VFLW*AZMAX1(trc_soHml2(ids_H2PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FINHM(M,NU(NY,NX),NY,NX).LT.0.0_r8)THEN
    IF(VLWatMicPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO  K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,NU(NY,NX),NY,NX))
    ENDDO
    SAdvFlx(idg_CO2)=VFLW*AZMAX1(trc_solml2(idg_CO2,NU(NY,NX),NY,NX))
    SAdvFlx(idg_CH4)=VFLW*AZMAX1(trc_solml2(idg_CH4,NU(NY,NX),NY,NX))
    SAdvFlx(idg_O2)=VFLW*AZMAX1(trc_solml2(idg_O2,NU(NY,NX),NY,NX))
    SAdvFlx(idg_N2)=VFLW*AZMAX1(trc_solml2(idg_N2,NU(NY,NX),NY,NX))
    SAdvFlx(idg_N2O)=VFLW*AZMAX1(trc_solml2(idg_N2O,NU(NY,NX),NY,NX))
    SAdvFlx(idg_H2)=VFLW*AZMAX1(trc_solml2(idg_H2,NU(NY,NX),NY,NX))

    SAdvFlx(ids_NH4)=VFLW*AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3)=VFLW*AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3)=VFLW*AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2)=VFLW*AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4)=VFLW*AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4)=VFLW*AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NH4B)=VFLW*AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3B)=VFLW*AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3B)=VFLW*AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2B)=VFLW*AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX))*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4B)=VFLW*AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4B)=VFLW*AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    DO  K=1,jcplx
      RFLOC(K)=0.0_r8
      RFLON(K)=0.0_r8
      RFLOP(K)=0.0_r8
      RFLOA(K)=0.0_r8
    ENDDO
    SAdvFlx(ids_beg:ids_end)=0._r8
  ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
!     VLWatMicPM,VLWatMacPM=micropore,macropore water volume
!     XFRS*VOLT=maximum macropore volume for solute transfer
!     DFV*=diffusive macropore-micropore solute transfer
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *H2,*2=macropore,micropore solute content
!
  IF(VLWatMacPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMacPS=AMIN1(XFRS*VGeomLayer(NU(NY,NX),NY,NX),VLWatMacPM(M,NU(NY,NX),NY,NX))
    VOLWT=VLWatMicPM(M,NU(NY,NX),NY,NX)+VLWatMacPS
    DO  K=1,jcplx
      DFVOC(K)=XNPH*(AZMAX1(OQCH2(K,NU(NY,NX),NY,NX))*VLWatMicPM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQC2(K,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
      DFVON(K)=XNPH*(AZMAX1(OQNH2(K,NU(NY,NX),NY,NX))*VLWatMicPM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQN2(K,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
      DFVOP(K)=XNPH*(AZMAX1(OQPH2(K,NU(NY,NX),NY,NX))*VLWatMicPM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQP2(K,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
      DFVOA(K)=XNPH*(AZMAX1(OQAH2(K,NU(NY,NX),NY,NX))*VLWatMicPM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQA2(K,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    ENDDO

    SDifFlx(idg_CO2)=XNPH*(AZMAX1(trc_soHml2(idg_CO2,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_CO2,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    SDifFlx(idg_CH4)=XNPH*(AZMAX1(trc_soHml2(idg_CH4,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_CH4,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    SDifFlx(idg_O2)=XNPH*(AZMAX1(trc_soHml2(idg_O2,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_O2,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    SDifFlx(idg_N2)=XNPH*(AZMAX1(trc_soHml2(idg_N2,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_N2,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    SDifFlx(idg_N2O)=XNPH*(AZMAX1(trc_soHml2(idg_N2O,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_N2O,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    SDifFlx(idg_H2)=XNPH*(AZMAX1(trc_soHml2(idg_H2,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_H2,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    SDifFlx(ids_NH4)=XNPH*(AZMAX1(trc_soHml2(ids_NH4,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3)=XNPH*(AZMAX1(trc_soHml2(idg_NH3,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    SDifFlx(ids_NO3)=XNPH*(AZMAX1(trc_soHml2(ids_NO3,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    SDifFlx(ids_NO2)=XNPH*(AZMAX1(trc_soHml2(ids_NO2,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    SDifFlx(ids_H1PO4)=XNPH*(AZMAX1(trc_soHml2(ids_H1PO4,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    SDifFlx(ids_H2PO4)=XNPH*(AZMAX1(trc_soHml2(ids_H2PO4,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    SDifFlx(ids_NH4B)=XNPH*(AZMAX1(trc_soHml2(ids_NH4B,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3B)=XNPH*(AZMAX1(trc_soHml2(idg_NH3B,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    SDifFlx(ids_NO3B)=XNPH*(AZMAX1(trc_soHml2(ids_NO3B,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    SDifFlx(ids_NO2B)=XNPH*(AZMAX1(trc_soHml2(ids_NO2B,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    SDifFlx(ids_H1PO4B)=XNPH*(AZMAX1(trc_soHml2(ids_H1PO4B,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
    SDifFlx(ids_H2PO4B)=XNPH*(AZMAX1(trc_soHml2(ids_H2PO4B,NU(NY,NX),NY,NX)) &
      *VLWatMicPM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  ELSE
    DO  K=1,jcplx
      DFVOC(K)=0.0_r8
      DFVON(K)=0.0_r8
      DFVOP(K)=0.0_r8
      DFVOA(K)=0.0_r8
    ENDDO
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  DO  K=1,jcplx
    ROCFXS(K,NU(NY,NX),NY,NX)=RFLOC(K)+DFVOC(K)
    RONFXS(K,NU(NY,NX),NY,NX)=RFLON(K)+DFVON(K)
    ROPFXS(K,NU(NY,NX),NY,NX)=RFLOP(K)+DFVOP(K)
    ROAFXS(K,NU(NY,NX),NY,NX)=RFLOA(K)+DFVOA(K)
  ENDDO

  DO NTS=ids_beg,ids_end
    RporeSoXFlx(NTS,NU(NY,NX),NY,NX)=SAdvFlx(NTS)+SDifFlx(NTS)
  ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FXS=hourly convective + diffusive solute flux between macropores and micropores
!     R*FXS=total convective + diffusive solute flux between macropores and micropores
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  DO  K=1,jcplx
    XOCFXS(K,NU(NY,NX),NY,NX)=XOCFXS(K,NU(NY,NX),NY,NX)+ROCFXS(K,NU(NY,NX),NY,NX)
    XONFXS(K,NU(NY,NX),NY,NX)=XONFXS(K,NU(NY,NX),NY,NX)+RONFXS(K,NU(NY,NX),NY,NX)
    XOPFXS(K,NU(NY,NX),NY,NX)=XOPFXS(K,NU(NY,NX),NY,NX)+ROPFXS(K,NU(NY,NX),NY,NX)
    XOAFXS(K,NU(NY,NX),NY,NX)=XOAFXS(K,NU(NY,NX),NY,NX)+ROAFXS(K,NU(NY,NX),NY,NX)
  ENDDO

  DO NTS=ids_beg,ids_end
    trcs_XFXS(NTS,NU(NY,NX),NY,NX)=trcs_XFXS(NTS,NU(NY,NX),NY,NX)&
      +RporeSoXFlx(NTS,NU(NY,NX),NY,NX)
  ENDDO

!
  end subroutine MacMicPoresTransfer
!------------------------------------------------------------------------------------------

  subroutine OverlandFlowSnowdriftTransport(M,NY,NX,NHE,NHW,NVS,NVN)
  implicit none

  integer, intent(in) :: M, NY, NX, NHE, NHW, NVS, NVN
  real(r8) :: FQRM,VFLW
  integer :: K,N,NN,N1,N2,N4,N5,N4B,N5B,NTG,NTN
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLWatMicPM=litter water volume from watsub.f
!     *S2=litter solute content
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     X*QRS=accumulated hourly solute in runoff
!
  N1=NX
  N2=NY
  IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VLWatMicPM(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,QRM(M,N2,N1)/VLWatMicPM(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      RQROC0(K,N2,N1)=VFLW*AZMAX1(OQC2(K,0,N2,N1))
      RQRON0(K,N2,N1)=VFLW*AZMAX1(OQN2(K,0,N2,N1))
      RQROP0(K,N2,N1)=VFLW*AZMAX1(OQP2(K,0,N2,N1))
      RQROA0(K,N2,N1)=VFLW*AZMAX1(OQA2(K,0,N2,N1))
    ENDDO
    DO NTG=idg_beg,idg_end-1
      trcg_RQR0(NTG,N2,N1)=VFLW*AZMAX1(trc_solml2(NTG,0,N2,N1))
    ENDDO

    DO NTN=ids_nut_beg,ids_nuts_end
      trcn_RQR0(NTN,N2,N1)=VFLW*AZMAX1(trc_solml2(NTN,0,N2,N1))
    ENDDO
  ELSE
    DO K=1,jcplx
      RQROC0(K,N2,N1)=0.0_r8
      RQRON0(K,N2,N1)=0.0_r8
      RQROP0(K,N2,N1)=0.0_r8
      RQROA0(K,N2,N1)=0.0_r8
    ENDDO
    trcg_RQR0(idg_beg:idg_end-1,N2,N1)=0.0_r8

    trcn_RQR0(ids_nut_beg:ids_nuts_end,N2,N1)=0.0_r8
  ENDIF
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  D4310: DO N=1,2
    D4305: DO NN=1,2
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE.AND.NN.EQ.1.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS.AND.NN.EQ.1.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N4B=NX
          N5B=NY-1
        ENDIF
      ENDIF
!
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLWatMicPM=litter water volume from watsub.f
!     *S2=litter solute content
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.1)THEN
          FQRM=QRMN(M,N,2,N5,N4)/QRM(M,N2,N1)
          DO  K=1,jcplx
            RQROC(K,N,2,N5,N4)=RQROC0(K,N2,N1)*FQRM
            RQRON(K,N,2,N5,N4)=RQRON0(K,N2,N1)*FQRM
            RQROP(K,N,2,N5,N4)=RQROP0(K,N2,N1)*FQRM
            RQROA(K,N,2,N5,N4)=RQROA0(K,N2,N1)*FQRM
          ENDDO

          DO NTG=idg_beg,idg_end-1
            trcg_RQR(NTG,N,2,N5,N4)=trcg_RQR0(NTG,N2,N1)*FQRM
          ENDDO

          DO NTN=ids_nut_beg,ids_nuts_end
            trcn_RQR(NTN,N,2,N5,N4)=trcn_RQR0(NTN,N2,N1)*FQRM
          ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
          DO  K=1,jcplx
            XOCQRS(K,N,2,N5,N4)=XOCQRS(K,N,2,N5,N4)+RQROC(K,N,2,N5,N4)
            XONQRS(K,N,2,N5,N4)=XONQRS(K,N,2,N5,N4)+RQRON(K,N,2,N5,N4)
            XOPQRS(K,N,2,N5,N4)=XOPQRS(K,N,2,N5,N4)+RQROP(K,N,2,N5,N4)
            XOAQRS(K,N,2,N5,N4)=XOAQRS(K,N,2,N5,N4)+RQROA(K,N,2,N5,N4)
          ENDDO
          DO NTG=idg_beg,idg_end-1
            trcg_XRS(NTG,N,2,N5,N4)=trcg_XRS(NTG,N,2,N5,N4)+trcg_RQR(NTG,N,2,N5,N4)
          ENDDO

          DO NTN=ids_nut_beg,ids_nuts_end
            trcn_XRS(NTN,N,2,N5,N4)=trcn_XRS(NTN,N,2,N5,N4)+trcn_RQR(NTN,N,2,N5,N4)
          ENDDO
        ELSE
          DO K=1,jcplx
            RQROC(K,N,2,N5,N4)=0.0_r8
            RQRON(K,N,2,N5,N4)=0.0_r8
            RQROP(K,N,2,N5,N4)=0.0_r8
            RQROA(K,N,2,N5,N4)=0.0_r8
          ENDDO
          trcg_RQR(idg_beg:idg_end-1,N,2,N5,N4)=0.0_r8

          trcn_RQR(ids_nut_beg:ids_nuts_end,N,2,N5,N4)=0.0_r8
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            FQRM=QRMN(M,N,1,N5B,N4B)/QRM(M,N2,N1)
            DO  K=1,jcplx
              RQROC(K,N,1,N5B,N4B)=RQROC0(K,N2,N1)*FQRM
              RQRON(K,N,1,N5B,N4B)=RQRON0(K,N2,N1)*FQRM
              RQROP(K,N,1,N5B,N4B)=RQROP0(K,N2,N1)*FQRM
              RQROA(K,N,1,N5B,N4B)=RQROA0(K,N2,N1)*FQRM
            ENDDO
            DO NTG=idg_beg,idg_end-1
              trcg_RQR(NTG,N,1,N5B,N4B)=trcg_RQR0(NTG,N2,N1)*FQRM
            ENDDO

            DO NTN=ids_nut_beg,ids_nuts_end
              trcn_RQR(NTN,N,1,N5B,N4B)=trcn_RQR0(NTN,N2,N1)*FQRM
            ENDDO

            DO K=1,jcplx
              XOCQRS(K,N,1,N5B,N4B)=XOCQRS(K,N,1,N5B,N4B)+RQROC(K,N,1,N5B,N4B)
              XONQRS(K,N,1,N5B,N4B)=XONQRS(K,N,1,N5B,N4B)+RQRON(K,N,1,N5B,N4B)
              XOPQRS(K,N,1,N5B,N4B)=XOPQRS(K,N,1,N5B,N4B)+RQROP(K,N,1,N5B,N4B)
              XOAQRS(K,N,1,N5B,N4B)=XOAQRS(K,N,1,N5B,N4B)+RQROA(K,N,1,N5B,N4B)
            ENDDO
            DO NTG=idg_beg,idg_end-1
              trcg_XRS(NTG,N,1,N5B,N4B)=trcg_XRS(NTG,N,1,N5B,N4B)+trcg_RQR(NTG,N,1,N5B,N4B)
            ENDDO

            DO NTN=ids_nut_beg,ids_nuts_end
              trcn_XRS(NTN,N,1,N5B,N4B)=trcn_XRS(NTN,N,1,N5B,N4B)+trcn_RQR(NTN,N,1,N5B,N4B)
            ENDDO
          ELSE
            DO  K=1,jcplx
              RQROC(K,N,1,N5B,N4B)=0.0_r8
              RQRON(K,N,1,N5B,N4B)=0.0_r8
              RQROP(K,N,1,N5B,N4B)=0.0_r8
              RQROA(K,N,1,N5B,N4B)=0.0_r8
            ENDDO
            trcg_RQR(idg_beg:idg_end,N,1,N5B,N4B)=0.0_r8
            trcn_RQR(ids_nut_beg:ids_nuts_end,N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        DO K=1,jcplx
          RQROC(K,N,2,N5,N4)=0.0_r8
          RQRON(K,N,2,N5,N4)=0.0_r8
          RQROP(K,N,2,N5,N4)=0.0_r8
          RQROA(K,N,2,N5,N4)=0.0_r8
        ENDDO
        trcg_RQR(idg_beg:idg_end-1,N,2,N5,N4)=0.0_r8
        trcn_RQR(ids_nut_beg:ids_nuts_end,N,2,N5,N4)=0.0_r8

        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          DO  K=1,jcplx
            RQROC(K,N,1,N5B,N4B)=0.0_r8
            RQRON(K,N,1,N5B,N4B)=0.0_r8
            RQROP(K,N,1,N5B,N4B)=0.0_r8
            RQROA(K,N,1,N5B,N4B)=0.0_r8
          ENDDO
          trcg_RQR(idg_beg:idg_end-1,N,1,N5B,N4B)=0.0_r8
          trcn_RQR(ids_nut_beg:ids_nuts_end,N,1,N5B,N4B)=0.0_r8
        ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     subroutine SnowdriftTransport(M)
!
!     QSM=snow transfer from watsub.f
!     RQS*=solute flux in snow transfer
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VOLS=volume of snowpack from watsub.f
!     *W2=solute content of snowpack
!
      IF(NN.EQ.1)THEN
!
!     IF NO SNOW DRIFT THEN NO TRANSPORT
!
        IF(ABS(QSM(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
          trcg_RQS(idg_beg:idg_end-1,N,N5,N4)=0.0_r8

          RQSNH4(N,N5,N4)=0.0_r8
          RQSNO3(N,N5,N4)=0.0_r8
          RQSH1P(N,N5,N4)=0.0_r8
          RQSH2P(N,N5,N4)=0.0_r8
    !
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).GT.ZEROS2(N2,N1))THEN
          IF(VOLS(N2,N1).GT.ZEROS2(N2,N1))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,QSM(M,N,N5,N4)/VOLS(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF

          DO NTG=idg_beg,idg_end-1
            trcg_RQS(NTG,N,N5,N4)=VFLW*AZMAX1(trcg_solsml2(NTG,1,N2,N1))
          ENDDO

          RQSNH4(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_NH4,1,N2,N1))
          RQSNO3(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_NO3,1,N2,N1))
          RQSH1P(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_H1PO4,1,N2,N1))
          RQSH2P(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_H2PO4,1,N2,N1))
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
          IF(VOLS(N5,N4).GT.ZEROS2(N5,N4))THEN
            VFLW=AZMIN1(AMAX1(-VFLWX,QSM(M,N,N5,N4)/VOLS(N5,N4)))
          ELSE
            VFLW=-VFLWX
          ENDIF

          DO NTG=idg_beg,idg_end-1
            trcg_RQS(NTG,N,N5,N4)=VFLW*AZMAX1(trcg_solsml2(NTG,1,N5,N4))
          ENDDO

          RQSNH4(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_NH4,1,N5,N4))
          RQSNO3(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_NO3,1,N5,N4))
          RQSH1P(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_H1PO4,1,N5,N4))
          RQSH2P(N,N5,N4)=VFLW*AZMAX1(trcn_solsml2(ids_H2PO4,1,N5,N4))
        ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QSS=hourly solute in snow flux
!     RQS*=solute in snow flux
!
        XCOQSS(N,N5,N4)=XCOQSS(N,N5,N4)+trcg_RQS(idg_CO2,N,N5,N4)
        XCHQSS(N,N5,N4)=XCHQSS(N,N5,N4)+trcg_RQS(idg_CH4,N,N5,N4)
        XOXQSS(N,N5,N4)=XOXQSS(N,N5,N4)+trcg_RQS(idg_O2,N,N5,N4)
        XNGQSS(N,N5,N4)=XNGQSS(N,N5,N4)+trcg_RQS(idg_N2,N,N5,N4)
        XN2QSS(N,N5,N4)=XN2QSS(N,N5,N4)+trcg_RQS(idg_N2O,N,N5,N4)
        XN3QSS(N,N5,N4)=XN3QSS(N,N5,N4)+trcg_RQS(idg_NH3,N,N5,N4)

        XN4QSS(N,N5,N4)=XN4QSS(N,N5,N4)+RQSNH4(N,N5,N4)
        XNOQSS(N,N5,N4)=XNOQSS(N,N5,N4)+RQSNO3(N,N5,N4)
        XP1QSS(N,N5,N4)=XP1QSS(N,N5,N4)+RQSH1P(N,N5,N4)
        XP4QSS(N,N5,N4)=XP4QSS(N,N5,N4)+RQSH2P(N,N5,N4)
      ENDIF
    ENDDO D4305
  ENDDO D4310
  end subroutine OverlandFlowSnowdriftTransport
!------------------------------------------------------------------------------------------

  subroutine LitterGasVolatilDissol(M,NY,NX)
  implicit none

  integer, intent(in) :: M,NY, NX
  real(r8) :: VOLGas(idg_beg:idg_end-1)
  real(r8) :: trc_gascl0(idg_beg:idg_end-1)
  integer  :: NTG
!
!     VOLT=litter volume from hour1.f
!     VLWatMicPM,VLsoiAirPM=micropore water volume, air volume from watsub.f
!     VLWatMicP*=equivalent aqueous volume for gas
!     S*L=solubility of gas in water from hour1.f
!     C*G=gas concentration
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     R*DFG=surface gas volatilization
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *S2=litter solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     R*DXR=gas exchange between atmosphere and surface litter water for gas flux calculations
!
  IF(VGeomLayer(0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VLsoiAirPM(M,0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VLWatMicPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN

    VLWatMicPCO(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_CO2,0,NY,NX)
    VLWatMicPCH(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_CH4,0,NY,NX)
    VLWatMicPOX(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_O2,0,NY,NX)
    VLWatMicPNG(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_N2,0,NY,NX)
    VLWatMicPN2(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_N2O,0,NY,NX)
    VLWatMicPN3(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_NH3,0,NY,NX)
    VLWatMicPHG(0,NY,NX)=VLWatMicPM(M,0,NY,NX)*GSolbility(idg_H2,0,NY,NX)

    VLWatMicPXA(0,NY,NX)=natomw*VLWatMicPM(M,0,NY,NX)
    DO NTG=idg_beg,idg_end-1
      trc_gascl0(NTG)=trc_gascl(NTG,0,NY,NX)*VLsoiAirPM(M,0,NY,NX)
    ENDDO

    VOLGas(idg_CO2)=VLWatMicPCO(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)
    VOLGas(idg_CH4)=VLWatMicPCH(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)
    VOLGas(idg_O2)=VLWatMicPOX(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)
    VOLGas(idg_N2)=VLWatMicPNG(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)
    VOLGas(idg_N2O)=VLWatMicPN2(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)
    VOLGas(idg_NH3)=VLWatMicPN3(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)
    VOLGas(idg_H2)=VLWatMicPHG(0,NY,NX)+VLsoiAirPM(M,0,NY,NX)

    RGasDSFlx(idg_CO2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_CO2))*VLWatMicPCO(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_CO2,0,NY,NX)+RCODXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_CO2)
    RGasDSFlx(idg_CH4,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_CH4))*VLWatMicPCH(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_CH4,0,NY,NX)+RCHDXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_CH4)
    RGasDSFlx(idg_O2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_O2))*VLWatMicPOX(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_O2,0,NY,NX)+ROXDXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_O2)
    RGasDSFlx(idg_N2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_N2))*VLWatMicPNG(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_N2,0,NY,NX)+RNGDXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_N2)
    RGasDSFlx(idg_N2O,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_N2O))*VLWatMicPN2(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_N2O,0,NY,NX)+RN2DXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_N2O)
    RGasDSFlx(idg_NH3,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_NH3))*VLWatMicPN3(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_NH3,0,NY,NX)+RN3DXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_NH3)
    RGasDSFlx(idg_H2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),trc_gascl0(idg_H2))*VLWatMicPHG(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_H2,0,NY,NX)+RHGDXR)*VLsoiAirPM(M,0,NY,NX))/VOLGas(idg_H2)

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly surface gas volatilization
!     R*DFG=surface gas volatilization
!   does not include band
    DO NTG=idg_beg,idg_end-1
      GasDisFlx(NTG,0,NY,NX)=GasDisFlx(NTG,0,NY,NX)+RGasDSFlx(NTG,0,NY,NX)
    ENDDO

  ELSE
    RGasDSFlx(idg_beg:idg_end,0,NY,NX)=0.0_r8
  ENDIF
  end subroutine LitterGasVolatilDissol

!------------------------------------------------------------------------------------------

  subroutine SurfSoilFluxGasDifAdv(M,NY,NX,FLQM,RDXS_gas)
!
! DESCRIPTION:
! surface soil gaseous diffusion, advection, dissolution & volatilization
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8),intent(in) :: FLQM(3,JD,JV,JH)
  real(r8),intent(in) :: RDXS_gas(idg_beg:idg_end)
  real(r8) :: VFLW
  integer :: NTG
!
!     THETPM=air-filled porosity from watsub.f
!     SoiBulkDensity=bulk density
!
  IF(THETPM(M,NU(NY,NX),NY,NX).GT.THETX.AND.SoiBulkDensity(NU(NY,NX),NY,NX).GT.ZERO)THEN
!
    call SurfSoilDifFlux(M,NY,NX)

!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SURFACE WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE SOIL SURFACE
!     OR THE ATMOSPHERE DEPENDING ON WATER FLUX DIRECTION
!
    call SurfSoillAdvFlux(M,NY,NX,FLQM)
!
!
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLG=hourly convective+diffusive gas flux
!
    DO NTG=idg_beg,idg_end-1
      R3GasADTFlx(NTG,3,NU(NY,NX),NY,NX)=R3GasADTFlx(NTG,3,NU(NY,NX),NY,NX) &
        +R3GasADFlx(NTG,3,NU(NY,NX),NY,NX)
    ENDDO

    IF(VLWatMicPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      call SurfSoilFluxDisolVapor(M,NY,NX,RDXS_gas)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      DO NTG=idg_beg,idg_end
        GasDisFlx(NTG,NU(NY,NX),NY,NX)=GasDisFlx(NTG,NU(NY,NX),NY,NX) &
          +RGasDSFlx(NTG,NU(NY,NX),NY,NX)
      ENDDO
    ELSE
      RGasDSFlx(idg_beg:idg_end,NU(NY,NX),NY,NX)=0.0_r8
    ENDIF
  ELSE
    R3GasADFlx(idg_beg:idg_end-1,3,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_beg:idg_end,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  end subroutine SurfSoilFluxGasDifAdv

!------------------------------------------------------------------------------------------

  subroutine LitterAtmosExchange(M,NY,NX,trcs_cl1)
  implicit none
  integer, intent(in) :: NY,NX,M
  real(r8),intent(out) :: trcs_cl1(ids_beg:ids_end)
  integer :: K,NTS,NTG

  real(r8) :: trc_gasq(idg_beg:idg_end-1)
  real(r8) :: DFGcc(idg_beg:idg_end-1)
  real(r8) :: DLYR0,TORT0

  IF(VGeomLayer(0,NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX)) !vertical layer thickness
    TORT0=TORT(M,0,NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR0)*CVRD(NY,NX)

    DO NTG=idg_beg,idg_end-1
      DFGcc(NTG)=SolDifcc(NTG,0,NY,NX)*TORT0
    ENDDO

    DO  K=1,jcplx
      COQC1(K)=AZMAX1(OQC2(K,0,NY,NX)/VLWatMicPM(M,0,NY,NX))
      COQN1(K)=AZMAX1(OQN2(K,0,NY,NX)/VLWatMicPM(M,0,NY,NX))
      COQP1(K)=AZMAX1(OQP2(K,0,NY,NX)/VLWatMicPM(M,0,NY,NX))
      COQA1(K)=AZMAX1(OQA2(K,0,NY,NX)/VLWatMicPM(M,0,NY,NX))
    ENDDO


    DO NTS=ids_beg,ids_end
      trcs_cl1(NTS)=AZMAX1(trc_solml2(NTS,0,NY,NX)/VLWatMicPM(M,0,NY,NX))
    ENDDO
!
!     EQUILIBRIUM CONCENTRATIONS AT RESIDUE SURFACE AT WHICH
!     AQUEOUS DIFFUSION THROUGH RESIDUE SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM AQUEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     PARR=boundary layer conductance above litter surface from watsub.f
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at litter surface
!     S*L=solubility of gas in water from hour1.f
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DFGS*=effective solute diffusivity
!
    DO NTG=idg_beg,idg_end-1
      trc_gasq(NTG)=(PARR(NY,NX)*AtmGgms(NTG,NY,NX)*GSolbility(NTG,0,NY,NX) &
        +DFGcc(NTG)*trcs_cl1(NTG))/(DFGcc(NTG)+PARR(NY,NX))

!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND RESIDUE SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!
!     R*DFR,X*DFR=current,hourly gas exchange between atmosphere and surface litter water
!     gas code: *CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at litter surface
!     VLWatMicPM=litter water volume
!     DFGS*=effective solute diffusivity
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     R*DXR=R*DFR for gas flux calculations
!
      RDFR_gas(NTG,NY,NX)=(trc_gasq(NTG)-trcs_cl1(NTG))*AMIN1(VLWatMicPM(M,0,NY,NX),DFGcc(NTG))
    ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
    DO NTG=idg_beg,idg_end-1
      trcg_XDFR(NTG,NY,NX)=trcg_XDFR(NTG,NY,NX)+RDFR_gas(NTG,NY,NX)
    ENDDO

  ELSE
    RDFR_gas(idg_beg:idg_end-1,NY,NX)=0.0_r8
  ENDIF

  RCODXR=RDFR_gas(idg_CO2,NY,NX)*XNPT
  RCHDXR=RDFR_gas(idg_CH4,NY,NX)*XNPT
  ROXDXR=RDFR_gas(idg_O2,NY,NX)*XNPT
  RNGDXR=RDFR_gas(idg_N2,NY,NX)*XNPT
  RN2DXR=RDFR_gas(idg_N2O,NY,NX)*XNPT
  RN3DXR=RDFR_gas(idg_NH3,NY,NX)*XNPT
  RHGDXR=RDFR_gas(idg_H2,NY,NX)*XNPT
  end subroutine LitterAtmosExchange
!------------------------------------------------------------------------------------------
  subroutine SoilAtmosExchange(M,NY,NX,trcs_cl2,RDXS_gas)

  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcs_cl2(ids_beg:ids_end)
  real(r8), intent(out) :: RDXS_gas(idg_beg:idg_end)
  real(r8) :: DLYR1,TORT1,VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: DFGS
  real(r8) :: trc_gsolc,trc_gsolc2
  integer  :: K,NTG
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE SOIL LAYER
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     TORT=tortuosity from hour1.f
!     *SGL*=solute diffusivity
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     DFGS*=effective solute diffusivity
!     C*2=solute concentration in soil surface
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in soil
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in litter
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!
  IF(VLWatMicPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMicPOA=VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    VLWatMicPOB=VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    VLWatMicPPA=VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    VLWatMicPPB=VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
    DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
    TORT1=TORT(M,NU(NY,NX),NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR1)

    D8910: DO  K=1,jcplx
      COQC2(K)=AZMAX1(OQC2(K,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX))
      COQN2(K)=AZMAX1(OQN2(K,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX))
      COQP2(K)=AZMAX1(OQP2(K,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX))
      COQA2(K)=AZMAX1(OQA2(K,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX))
    ENDDO D8910

!not include NH3 and NH3B
    DO NTG=idg_beg,idg_end-2
      trcs_cl2(NTG)=AZMAX1(trc_solml2(NTG,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX))
    ENDDO

    IF(VLWatMicPMA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl2(idg_NH3)=AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX)/VLWatMicPMA(NU(NY,NX),NY,NX))
      trcs_cl2(ids_NH4)=AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX)/VLWatMicPMA(NU(NY,NX),NY,NX))
    ELSE
      trcs_cl2(idg_NH3)=0.0_r8
      trcs_cl2(ids_NH4)=0.0_r8
    ENDIF

    IF(VLWatMicPOA.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_NO3)=AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX)/VLWatMicPOA)
      trcs_cl2(ids_NO2)=AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX)/VLWatMicPOA)
    ELSE
      trcs_cl2(ids_NO3)=0.0_r8
      trcs_cl2(ids_NO2)=0.0_r8
    ENDIF
    IF(VLWatMicPPA.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_H1PO4)=AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX)/VLWatMicPPA)
      trcs_cl2(ids_H2PO4)=AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX)/VLWatMicPPA)
    ELSE
      trcs_cl2(ids_H1PO4)=0.0_r8
      trcs_cl2(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VLWatMicPMB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl2(idg_NH3B)=AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX)/VLWatMicPMB(NU(NY,NX),NY,NX))
      trcs_cl2(ids_NH4B)=AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX)/VLWatMicPMB(NU(NY,NX),NY,NX))
    ELSE
      trcs_cl2(idg_NH3B)=trcs_cl2(idg_NH3)
      trcs_cl2(ids_NH4B)=trcs_cl2(ids_NH4)
    ENDIF

    IF(VLWatMicPOB.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_NO3B)=AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX)/VLWatMicPOB)
      trcs_cl2(ids_NO2)=AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX)/VLWatMicPOB)
    ELSE
      trcs_cl2(ids_NO3B)=trcs_cl2(ids_NO3)
      trcs_cl2(ids_NO2)=trcs_cl2(ids_NO2)
    ENDIF
    IF(VLWatMicPPB.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_H1PO4B)=AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX)/VLWatMicPPB)
      trcs_cl2(ids_H2PO4B)=AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX)/VLWatMicPPB)
    ELSE
      trcs_cl2(ids_H1PO4B)=trcs_cl2(ids_H1PO4)
      trcs_cl2(ids_H2PO4B)=trcs_cl2(ids_H2PO4)
    ENDIF
!
!     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
!     AQUEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM AQUEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     PARG=boundary layer conductance above soil surface from watsub.f
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at soil surface
!     S*L=solubility of gas in water from hour1.f
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DFGS*=effective solute diffusivity
!
!include NH3B
    DO NTG=idg_beg,idg_end
      DFGS=SolDifcc(NTG,NU(NY,NX),NY,NX)*TORT1
      trc_gsolc2=AZMAX1(trc_solml2(NTG,NU(NY,NX),NY,NX)/VLWatMicPM(M,NU(NY,NX),NY,NX))
      trc_gsolc=(PARG(M,NY,NX)*AtmGgms(NTG,NY,NX)*GSolbility(NTG,NU(NY,NX),NY,NX) &
        +DFGS*trc_gsolc2)/(DFGS+PARG(M,NY,NX))
!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND SOIL SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!     R*DFS, dissolution/volatilization
!     R*DFS,X*DFS=current,cumulative gas exchange between atmosphere and soil surface water
!     gas code:*CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at soil surface
!     VLWatMicPM=litter water volume
!     DFGS*=effective solute diffusivity
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     R*DXS=R*DFS for gas flux calculations
!  include NH3B

      RGasSSVol(NTG,NY,NX)=(trc_gsolc-trcs_cl2(NTG))&
        *AMIN1(VLWatMicPM(M,NU(NY,NX),NY,NX)*trcs_VLN(NTG,NU(NY,NX),NY,NX),DFGS)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!   seven gas species plus aqueous NH3 in band

      GasSfAtmFlx(NTG,NY,NX)=GasSfAtmFlx(NTG,NY,NX)+RGasSSVol(NTG,NY,NX)
    ENDDO

  ELSE
    RGasSSVol(idg_beg:idg_end,NY,NX)=0.0_r8
  ENDIF

  DO NTG=idg_beg,idg_end
    RDXS_gas(NTG)=RGasSSVol(NTG,NY,NX)*XNPT
  ENDDO
  end subroutine SoilAtmosExchange
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSnowpackDisch(M,NY,NX,trcg_RFLS0,trcn_RFLW0)

  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8), intent(out) :: trcg_RFLS0(idg_beg:idg_end-1)
  real(r8), intent(out) :: trcn_RFLW0(ids_nut_beg:ids_nuts_end)
  integer :: ICHKL,L,L2
  real(r8) :: VFLWS,VFLWW,VFLWR
  real(r8) :: VFLWNOB,VFLWNO3,VFLWNHB
  real(r8) :: VFLWNH4,VFLWPO4,VFLWPOB
  integer  :: NTG,NTS,NTN
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     FLQWM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  ICHKL=0
  DO  L=1,JS
    IF(VHCPWM(M,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWM(M,L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0,FLQWM(M,L2,NY,NX)/VOLWSL(L,NY,NX)))
        ELSE
          VFLWW=1.0
        ENDIF
        DO NTG=idg_beg,idg_end-1
          trcg_RBLS(NTG,L2,NY,NX)=trcg_solsml2(NTG,L,NY,NX)*VFLWW
          trcg_XBLS(NTG,L2,NY,NX)=trcg_XBLS(NTG,L2,NY,NX)+trcg_RBLS(NTG,L2,NY,NX)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_RBLS(NTN,L2,NY,NX)=trcn_solsml2(NTN,L,NY,NX)*VFLWW
          trcn_XBLS(NTN,L2,NY,NX)=trcn_XBLS(NTN,L2,NY,NX)+trcn_RBLS(NTN,L2,NY,NX)
        ENDDO
      ELSE
        IF(L.LT.JS)THEN
          trcg_RBLS(idg_beg:idg_end-1,L2,NY,NX)=0.0_r8
          trcn_RBLS(ids_nut_beg:ids_nuts_end,L2,NY,NX)=0.0_r8
        ENDIF
!
!     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
!
!     VOLWSL=snowpack water content
!     FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
!     CVRD,BARE=litter cover fraction,1-CVRD
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     R*S0=solute flux to surface litter
!     R*S1,R*B1=solute flux to soil surface non-band,band
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
        IF(ICHKL.EQ.0)THEN
          IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,FLQRM(M,NY,NX)/VOLWSL(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(FLQSM(M,NY,NX)+FLQHM(M,NY,NX))/VOLWSL(L,NY,NX)))
          ELSE
            VFLWR=CVRD(NY,NX)
            VFLWS=BARE(NY,NX)
          ENDIF
          VFLWNH4=VFLWS*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
          VFLWNHB=VFLWS*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
          VFLWNO3=VFLWS*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
          VFLWNOB=VFLWS*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
          VFLWPO4=VFLWS*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)

          DO NTG=idg_beg,idg_end-1
            trcg_RFLS0(NTG)=trcg_solsml2(NTG,L,NY,NX)*VFLWR
          ENDDO

          DO NTS=ids_nut_beg,ids_nuts_end
            trcn_RFLW0(NTS)=trcn_solsml2(NTS,L,NY,NX)*VFLWR
          ENDDO

          RCOFLS1=trcg_solsml2(idg_CO2,L,NY,NX)*VFLWS
          RCHFLS1=trcg_solsml2(idg_CH4,L,NY,NX)*VFLWS
          ROXFLS1=trcg_solsml2(idg_O2,L,NY,NX)*VFLWS
          RNGFLS1=trcg_solsml2(idg_N2,L,NY,NX)*VFLWS
          RN2FLS1=trcg_solsml2(idg_N2O,L,NY,NX)*VFLWS
          RN3FLW1=trcg_solsml2(idg_NH3,L,NY,NX)*VFLWNH4

          RN4FLW1=trcn_solsml2(ids_NH4,L,NY,NX)*VFLWNH4
          RNOFLW1=trcn_solsml2(ids_NO3,L,NY,NX)*VFLWNO3
          RH1PFS1=trcn_solsml2(ids_H1PO4,L,NY,NX)*VFLWPO4
          RH2PFS1=trcn_solsml2(ids_H2PO4,L,NY,NX)*VFLWPO4

          RN4FLB1=trcn_solsml2(ids_NH4,L,NY,NX)*VFLWNHB
          RN3FLB1=trcg_solsml2(idg_NH3,L,NY,NX)*VFLWNHB
          RNOFLB1=trcn_solsml2(ids_NO3,L,NY,NX)*VFLWNOB
          RH1BFB1=trcn_solsml2(ids_H1PO4,L,NY,NX)*VFLWPOB
          RH2BFB1=trcn_solsml2(ids_H2PO4,L,NY,NX)*VFLWPOB
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      trcg_RFLS0(idg_beg:idg_end-1)=0.0_r8
      trcn_RFLW0(ids_nut_beg:ids_nuts_end)=0.0_r8
      RCOFLS1=0.0_r8
      RCHFLS1=0.0_r8
      ROXFLS1=0.0_r8
      RNGFLS1=0.0_r8
      RN2FLS1=0.0_r8
      RN4FLW1=0.0_r8
      RN3FLW1=0.0_r8
      RNOFLW1=0.0_r8
      RH1PFS1=0.0_r8
      RH2PFS1=0.0_r8
      RN4FLB1=0.0_r8
      RN3FLB1=0.0_r8
      RNOFLB1=0.0_r8
      RH1BFB1=0.0_r8
      RH2BFB1=0.0_r8
    ENDIF
  ENDDO
  end subroutine SoluteFluxSnowpackDisch

! ----------------------------------------------------------------------
  subroutine SurfSoillAdvFlux(M,NY,NX,FLQM)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLQM(3,JD,JV,JH)
  real(r8) :: VFLW
  real(r8) :: RFLg_ADV(idg_beg:idg_end-1)
  integer :: NTG

!     FLQM=total water flux into soil micropore+macropore from watsub.f
!     VLsoiAirPM=air-filled porosity
!     RFL*G=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     *G2=gaseous content
!     C*E=atmospheric gas concentration from hour1.f
!

  IF(FLQM(3,NU(NY,NX),NY,NX).GT.0.0_r8)THEN
    IF(VLsoiAirPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQM(3,NU(NY,NX),NY,NX)/VLsoiAirPM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO NTG=idg_beg,idg_end-1
      RFLg_ADV(NTG)=VFLW*AZMAX1(trc_gasml2(NTG,NU(NY,NX),NY,NX))
    ENDDO
  ELSE
    DO NTG=idg_beg,idg_end-1
      RFLg_ADV(NTG)=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(NTG,NY,NX)
    ENDDO
  ENDIF
!     TOTAL SOIL GAS FLUX + CONVECTIVE FLUX
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux
  DO NTG=idg_beg,idg_end-1
    R3GasADFlx(NTG,3,NU(NY,NX),NY,NX)=R3GasADFlx(NTG,3,NU(NY,NX),NY,NX)+RFLg_ADV(NTG)
  ENDDO
  end subroutine SurfSoillAdvFlux

! ----------------------------------------------------------------------
  subroutine SurfSoilDifFlux(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: DFV_g
  real(r8) :: DFLG2
  real(r8) :: trcg_cl2
  real(r8) :: DGQ_cef
  integer  :: NTG
!     GASEOUS DIFFUSIVITIES
!
!     DFLG2=air-filled porosity effect on gaseous diffusivity
!     POROQ=Penman Water Linear Reduction tortuosity from starts.f
!     POROS=total porosity
!     DLYR=soil layer thickness
!     D*G=gaseous diffusivity in soil
!     *SGL2= gaseous diffusivity in air
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
  DFLG2=AZMAX1(THETPM(M,NU(NY,NX),NY,NX))*POROQ &
    *THETPM(M,NU(NY,NX),NY,NX)/POROS(NU(NY,NX),NY,NX) &
    *AREA(3,NU(NY,NX),NY,NX)/AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))

  DO NTG=idg_beg,idg_end-1
    DifuscG(NTG,3,NU(NY,NX),NY,NX)=DFLG2*GasDifcc(NTG,NU(NY,NX),NY,NX)
!
!     SURFACE GAS CONCENTRATIONS
!
!     C*G2=gaseous concentration
!     *G2=gaseous content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     VLsoiAirPM=air-filled porosity
!
    trcg_cl2=AZMAX1(trc_gasml2(NTG,NU(NY,NX),NY,NX)/VLsoiAirPM(M,NU(NY,NX),NY,NX))
!
!     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
!     GASEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM GASEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     D*GQ=gaseous diffusivity between soil surface and atmosphere
!     D*G=gaseous diffusivity in soil
!     PARG*=boundary layer conductance above soil surface
!     DFV*G=diffusive gas flux
!     C*E=atmospheric gas concentration from hour1.f
!     C*G2=gaseous concentration
!
    DGQ_cef=DifuscG(NTG,3,NU(NY,NX),NY,NX)*PARG_cef(NTG,NY,NX) &
      /(DifuscG(NTG,3,NU(NY,NX),NY,NX)+PARG_cef(NTG,NY,NX))

    DFV_g=DGQ_cef*(AtmGgms(NTG,NY,NX)-trcg_cl2)

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux

    R3GasADFlx(NTG,3,NU(NY,NX),NY,NX)=DFV_g
  ENDDO

  end subroutine SurfSoilDifFlux
! ----------------------------------------------------------------------

  subroutine SurfSoilFluxDisolVapor(M,NY,NX,RDXS_gas)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: RDXS_gas(idg_beg:idg_end)
  real(r8) :: CNH3B0,CNH4B0
  real(r8) :: VOLHGT(JY,JX),VOLNBT(JY,JX)
  real(r8) :: VOLN3T(JY,JX),VOLN2T(JY,JX)
  real(r8) :: VOLCOT(JY,JX),VOLCHT(JY,JX)
  real(r8) :: VOLOXT(JY,JX),VOLNGT(JY,JX)

!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     VLWatMicP*=equivalent aqueous volume for gas
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     S*L=solubility of gas in water from hour1.f
!     VLsoiAirPM=air-filled porosity
!     R*DFG=water-air gas flux
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     R*DXS=gas exchange between atmosphere and soil surface water
!

  VLWatMicPCO(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*GSolbility(idg_CO2,NU(NY,NX),NY,NX)
  VLWatMicPCH(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*GSolbility(idg_CH4,NU(NY,NX),NY,NX)
  VLWatMicPOX(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*GSolbility(idg_O2,NU(NY,NX),NY,NX)
  VLWatMicPNG(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*GSolbility(idg_N2,NU(NY,NX),NY,NX)
  VLWatMicPN2(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*GSolbility(idg_N2O,NU(NY,NX),NY,NX)
  VLWatMicPN3(NU(NY,NX),NY,NX)=VLWatMicPMA(NU(NY,NX),NY,NX)*GSolbility(idg_NH3,NU(NY,NX),NY,NX)
  VLWatMicPNB(NU(NY,NX),NY,NX)=VLWatMicPMB(NU(NY,NX),NY,NX)*GSolbility(idg_NH3,NU(NY,NX),NY,NX)
  VLWatMicPHG(NU(NY,NX),NY,NX)=VLWatMicPM(M,NU(NY,NX),NY,NX)*GSolbility(idg_H2,NU(NY,NX),NY,NX)

  VOLCOT(NY,NX)=VLWatMicPCO(NU(NY,NX),NY,NX)+VLsoiAirPM(M,NU(NY,NX),NY,NX)
  VOLCHT(NY,NX)=VLWatMicPCH(NU(NY,NX),NY,NX)+VLsoiAirPM(M,NU(NY,NX),NY,NX)
  VOLOXT(NY,NX)=VLWatMicPOX(NU(NY,NX),NY,NX)+VLsoiAirPM(M,NU(NY,NX),NY,NX)
  VOLNGT(NY,NX)=VLWatMicPNG(NU(NY,NX),NY,NX)+VLsoiAirPM(M,NU(NY,NX),NY,NX)
  VOLN2T(NY,NX)=VLWatMicPN2(NU(NY,NX),NY,NX)+VLsoiAirPM(M,NU(NY,NX),NY,NX)
  VOLN3T(NY,NX)=VLWatMicPN3(NU(NY,NX),NY,NX)+VLsoiAirPMA(NU(NY,NX),NY,NX)
  VOLNBT(NY,NX)=VLWatMicPNB(NU(NY,NX),NY,NX)+VLsoiAirPMB(NU(NY,NX),NY,NX)
  VOLHGT(NY,NX)=VLWatMicPHG(NU(NY,NX),NY,NX)+VLsoiAirPM(M,NU(NY,NX),NY,NX)

  RGasDSFlx(idg_CO2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_CO2,NU(NY,NX),NY,NX)) &
    *VLWatMicPCO(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_CO2,NU(NY,NX),NY,NX)+RDXS_gas(idg_CO2)) &
    *VLsoiAirPM(M,NU(NY,NX),NY,NX))/VOLCOT(NY,NX)
  RGasDSFlx(idg_CH4,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_CH4,NU(NY,NX),NY,NX)) &
    *VLWatMicPCH(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_CH4,NU(NY,NX),NY,NX)+RDXS_gas(idg_CH4)) &
    *VLsoiAirPM(M,NU(NY,NX),NY,NX))/VOLCHT(NY,NX)
  RGasDSFlx(idg_O2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_O2,NU(NY,NX),NY,NX)) &
    *VLWatMicPOX(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_O2,NU(NY,NX),NY,NX)+RDXS_gas(idg_O2)) &
    *VLsoiAirPM(M,NU(NY,NX),NY,NX))/VOLOXT(NY,NX)
  RGasDSFlx(idg_N2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_N2,NU(NY,NX),NY,NX)) &
    *VLWatMicPNG(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_N2,NU(NY,NX),NY,NX)+RDXS_gas(idg_N2)) &
  *VLsoiAirPM(M,NU(NY,NX),NY,NX))/VOLNGT(NY,NX)
  RGasDSFlx(idg_N2O,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_N2O,NU(NY,NX),NY,NX)) &
    *VLWatMicPN2(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_N2O,NU(NY,NX),NY,NX)+RDXS_gas(idg_N2O)) &
    *VLsoiAirPM(M,NU(NY,NX),NY,NX))/VOLN2T(NY,NX)
  RGasDSFlx(idg_H2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_H2,NU(NY,NX),NY,NX)) &
    *VLWatMicPHG(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_H2,NU(NY,NX),NY,NX)+RDXS_gas(idg_H2)) &
    *VLsoiAirPM(M,NU(NY,NX),NY,NX))/VOLHGT(NY,NX)

  IF(VOLN3T(NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPXA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_NH3,NU(NY,NX),NY,NX)) &
      *VLWatMicPN3(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,trc_solml2(idg_NH3,NU(NY,NX),NY,NX)+RDXS_gas(idg_NH3)) &
      *VLsoiAirPMA(NU(NY,NX),NY,NX))/VOLN3T(NY,NX)
  ELSE
    RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  IF(VOLNBT(NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPXB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    RGasDSFlx(idg_NH3B,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_NH3,NU(NY,NX),NY,NX)) &
      *VLWatMicPNB(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,trc_solml2(idg_NH3B,NU(NY,NX),NY,NX)+RDXS_gas(idg_NH3B)) &
      *VLsoiAirPMB(NU(NY,NX),NY,NX))/VOLNBT(NY,NX)
    CNH3B0=AZMAX1((trc_solml2(idg_NH3B,NU(NY,NX),NY,NX) &
      +RGasDSFlx(idg_NH3B,NU(NY,NX),NY,NX))/VLWatMicPXB(NU(NY,NX),NY,NX))
    CNH4B0=AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))/VLWatMicPXB(NU(NY,NX),NY,NX)
  ELSE
    RGasDSFlx(idg_NH3B,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  end subroutine SurfSoilFluxDisolVapor
! ----------------------------------------------------------------------
end module SurfaceFluxMod
