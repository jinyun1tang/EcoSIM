module SurfaceFluxMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
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
  real(r8) :: COXYS2,CZ2GS2,CZ2OS2,CH2GS2,CNH3S2,CNH4S2,CNO3S2
  real(r8) :: CNO2S2,CP14S2,CPO4S2,CCH4S2,CCO2S2
  real(r8) :: CNH3B2,CNH4B2,CNO2B2,CP14B2,CPO4B2,CNO3B2
  real(r8) :: RCODXS,RCHDXS,ROXDXS,RNGDXS,RN2DXS,RN3DXS,RNBDXS,RHGDXS
  real(r8) :: CCH4S1,COXYS1,CZ2GS1,CZ2OS1,CH2GS1,CNH4S1
  real(r8) :: CNH3S1,CNO3S1,CNO2S1,CP14S1,CPO4S1,CCO2S1
  real(r8) :: RCODXR,RCHDXR,ROXDXR,RNGDXR,RN2DXR,RN3DXR,RHGDXR
  real(r8) :: RFLCOS,RFLCHS,RFLNGS,RFLN2S,RFLOXS
  real(r8) :: RFLHGS,RFLNH4,RFLNH3,RFLNO3,RFLNO2,RFLP14,RFLPO4
  real(r8) :: RFLN4B,RFLN3B,RFLNOB,RFLN2B,RFLP1B,RFLPOB
  real(r8) :: RCOFLS1,RCHFLS1,RNGFLS1,RN2FLS1,ROXFLS1
  real(r8) :: RN4FLW1,RN3FLW1,RNOFLW1
  real(r8) :: RCOFLS0,RCHFLS0,RN2FLS0,ROXFLS0,RNGFLS0
  real(r8) :: RN4FLB1,RN3FLB1,RNOFLB1,RH1BFB1,RH2BFB1
  real(r8) :: RN4FLW0,RN3FLW0,RNOFLW0,RH1PFS0,RH2PFS0
  real(r8) :: RH1PFS1,RH2PFS1

  public :: SoluteFluxSurface
  public :: LitterGasVolatilDissol
  public :: SurfSoilFluxGasDifAdv
  public :: SoluteFluxSnowpackDisch
contains


!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSurface(M,NY,NX,NHE,NHW,NVS,NVN,FLQM)
  implicit none
  integer, intent(in) :: NY,NX,M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8),intent(inout) :: FLQM(3,JD,JV,JH)
  real(r8) :: FLWRM1
  real(r8) :: SDifFlx(ids_beg:ids_end)
!     VOLWM,VOLWHM,VOLPM,FLPM=micropore,macropore water volume, air volume and change in air volume
!     FLWM,FLWHM=water flux into soil micropore,macropore from watsub.f
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     FLVM,FLM=air,water flux in gas flux calculations
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
  VOLWMA(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*VLNH4(NU(NY,NX),NY,NX)
  VOLWMB(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*VLNHB(NU(NY,NX),NY,NX)
  VOLWXA(NU(NY,NX),NY,NX)=natomw*VOLWMA(NU(NY,NX),NY,NX)
  VOLWXB(NU(NY,NX),NY,NX)=natomw*VOLWMB(NU(NY,NX),NY,NX)
  VOLPMA(NU(NY,NX),NY,NX)=VOLPM(M,NU(NY,NX),NY,NX)*VLNH4(NU(NY,NX),NY,NX)
  VOLPMB(NU(NY,NX),NY,NX)=VOLPM(M,NU(NY,NX),NY,NX)*VLNHB(NU(NY,NX),NY,NX)
  FLVM(NU(NY,NX),NY,NX)=FLPM(M,NU(NY,NX),NY,NX)*XNPT
  FLQM(3,NU(NY,NX),NY,NX)=(FLWM(M,3,NU(NY,NX),NY,NX)+FLWHM(M,3,NU(NY,NX),NY,NX))*XNPT
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE RESIDUE
!
!     VOLT,DLYR,AREA=litter volume, thickness, area
!     VOLWM=micropore water-filled porosity from watsub.f
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
  call LitterAtmosExchange(M,NY,NX)

  call SoilAtmosExchange(M,NY,NX)

!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
  call ConvectiveSurfaceSoluteFlux(M,NY,NX,FLWRM1)
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  call LitterSurfSoilExchange(M,NY,NX,FLWRM1,SDifFlx)

!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES
!
  call TotalPoreFluxAdjacentCell(NY,NX,SDifFlx)

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

  subroutine ConvectiveSurfaceSoluteFlux(M,NY,NX,FLWRM1)
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8),intent(out) :: FLWRM1
  REAL(R8) :: VFLW
  integer :: K

  FLWRM1=FLWRM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!     VOLWM=litter water volume
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
    IF(VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VOLWM(M,0,NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,0,NY,NX))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,0,NY,NX))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,0,NY,NX))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,0,NY,NX))
    ENDDO
    RFLCOS=VFLW*AZMAX1(trc_solml2(idg_CO2,0,NY,NX))
    RFLCHS=VFLW*AZMAX1(trc_solml2(idg_CH4,0,NY,NX))
    RFLOXS=VFLW*AZMAX1(trc_solml2(idg_O2,0,NY,NX))
    RFLNGS=VFLW*AZMAX1(trc_solml2(idg_N2,0,NY,NX))
    RFLN2S=VFLW*AZMAX1(trc_solml2(idg_N2O,0,NY,NX))
    RFLHGS=VFLW*AZMAX1(trc_solml2(idg_H2,0,NY,NX))
    RFLNH4=VFLW*AZMAX1(trc_solml2(ids_NH4,0,NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNH3=VFLW*AZMAX1(trc_solml2(idg_NH3,0,NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNO3=VFLW*AZMAX1(trc_solml2(ids_NO3,0,NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLNO2=VFLW*AZMAX1(trc_solml2(ids_NO2,0,NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLP14=VFLW*AZMAX1(trc_solml2(ids_H1PO4,0,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLPO4=VFLW*AZMAX1(trc_solml2(ids_H2PO4,0,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLN4B=VFLW*AZMAX1(trc_solml2(ids_NH4,0,NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLN3B=VFLW*AZMAX1(trc_solml2(idg_NH3,0,NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLNOB=VFLW*AZMAX1(trc_solml2(ids_NO3,0,NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLN2B=VFLW*AZMAX1(trc_solml2(ids_NO2,0,NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLP1B=VFLW*AZMAX1(trc_solml2(ids_H1PO4,0,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
    RFLPOB=VFLW*AZMAX1(trc_solml2(ids_H2PO4,0,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE
!
!     VOLWM=micropore water-filled porosity from watsub.f
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
    IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VOLWM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,NU(NY,NX),NY,NX))
    ENDDO
    RFLCOS=VFLW*AZMAX1(trc_solml2(idg_CO2,NU(NY,NX),NY,NX))
    RFLCHS=VFLW*AZMAX1(trc_solml2(idg_CH4,NU(NY,NX),NY,NX))
    RFLOXS=VFLW*AZMAX1(trc_solml2(idg_O2,NU(NY,NX),NY,NX))
    RFLNGS=VFLW*AZMAX1(trc_solml2(idg_N2,NU(NY,NX),NY,NX))
    RFLN2S=VFLW*AZMAX1(trc_solml2(idg_N2O,NU(NY,NX),NY,NX))
    RFLHGS=VFLW*AZMAX1(trc_solml2(idg_H2,NU(NY,NX),NY,NX))
    RFLNH4=VFLW*AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX))
    RFLNH3=VFLW*AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX))
    RFLNO3=VFLW*AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX))
    RFLNO2=VFLW*AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX))
    RFLP14=VFLW*AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX))
    RFLPO4=VFLW*AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX))
    RFLN4B=VFLW*AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))
    RFLN3B=VFLW*AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX))
    RFLNOB=VFLW*AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX))
    RFLN2B=VFLW*AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX))
    RFLP1B=VFLW*AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX))
    RFLPOB=VFLW*AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX))
  ENDIF
  end subroutine ConvectiveSurfaceSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine LitterSurfSoilExchange(M,NY,NX,FLWRM1,SDifFlx)
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8), intent(in) :: FLWRM1
  real(r8), intent(out):: SDifFlx(ids_beg:ids_end)
  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: DIFOC,DIFON,DIFOP,DIFOA

  real(r8) :: DIFOC0,DIFON0,DIFOP0,DIFOA0
  real(r8) :: DIFOC1,DIFON1,DIFOP1,DIFOA1
  real(r8) :: DISPN,DIF0,DIF1
  real(r8) :: SDifc(ids_beg:ids_end)
  integer  :: K,NTS
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VOLWM=micropore water-filled porosity from watsub.f
!
  IF((VOLT(0,NY,NX).GT.ZEROS2(NY,NX).AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND.(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN
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
    SDifFlx(idg_CO2)=SDifc(idg_CO2)*(CCO2S1-CCO2S2)
    SDifFlx(idg_CH4)=SDifc(idg_CH4)*(CCH4S1-CCH4S2)
    SDifFlx(idg_O2) =SDifc(idg_O2)*(COXYS1-COXYS2)
    SDifFlx(idg_N2) =SDifc(idg_N2)*(CZ2GS1-CZ2GS2)
    SDifFlx(idg_N2O)=SDifc(idg_N2O)*(CZ2OS1-CZ2OS2)
    SDifFlx(idg_H2)=SDifc(idg_H2)*(CH2GS1-CH2GS2)
    SDifFlx(ids_NH4)=SDifc(ids_NH4)*(CNH4S1-CNH4S2)*VLNH4(NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3)=SDifc(idg_NH3)*(CNH3S1-CNH3S2)*VLNH4(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO3)=SDifc(ids_NO3)*(CNO3S1-CNO3S2)*VLNO3(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO2)=SDifc(ids_NO2)*(CNO2S1-CNO2S2)*VLNO3(NU(NY,NX),NY,NX)
    SDifFlx(ids_H1PO4)=SDifc(ids_H1PO4)*(CP14S1-CP14S2)*VLPO4(NU(NY,NX),NY,NX)
    SDifFlx(ids_H2PO4)=SDifc(ids_H2PO4)*(CPO4S1-CPO4S2)*VLPO4(NU(NY,NX),NY,NX)
    SDifFlx(ids_NH4B)=SDifc(ids_NH4B)*(CNH4S1-CNH4B2)*VLNHB(NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3B)=SDifc(idg_NH3B)*(CNH3S1-CNH3B2)*VLNHB(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO3B)=SDifc(ids_NO3B)*(CNO3S1-CNO3B2)*VLNOB(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO2B)=SDifc(ids_NO2B)*(CNO2S1-CNO2B2)*VLNOB(NU(NY,NX),NY,NX)
    SDifFlx(ids_H1PO4B)=SDifc(ids_H1PO4B)*(CP14S1-CP14B2)*VLPOB(NU(NY,NX),NY,NX)
    SDifFlx(ids_H2PO4B)=SDifc(ids_H2PO4B)*(CPO4S1-CPO4B2)*VLPOB(NU(NY,NX),NY,NX)
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

  subroutine TotalPoreFluxAdjacentCell(NY,NX,SDifFlx)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: NY, NX
  real(r8), intent(in) :: SDifFlx(ids_beg:ids_end)
  integer :: K
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
  RPoreSolFlx(idg_CO2,3,0,NY,NX)=RCOFL0(NY,NX)+RCOFLS0-RFLCOS-SDifFlx(idg_CO2)
  RPoreSolFlx(idg_CH4,3,0,NY,NX)=RCHFL0(NY,NX)+RCHFLS0-RFLCHS-SDifFlx(idg_CH4)
  RPoreSolFlx(idg_O2,3,0,NY,NX)=ROXFL0(NY,NX)+ROXFLS0-RFLOXS-SDifFlx(idg_O2)
  RPoreSolFlx(idg_N2,3,0,NY,NX)=RNGFL0(NY,NX)+RNGFLS0-RFLNGS-SDifFlx(idg_N2)
  RPoreSolFlx(idg_N2O,3,0,NY,NX)=RN2FL0(NY,NX)+RN2FLS0-RFLN2S-SDifFlx(idg_N2O)
  RPoreSolFlx(idg_H2,3,0,NY,NX)=RHGFL0(NY,NX)-RFLHGS-SDifFlx(idg_H2)
  RPoreSolFlx(ids_NH4,3,0,NY,NX)=RN4FL0(NY,NX)+RN4FLW0-RFLNH4-SDifFlx(ids_NH4)-RFLN4B-SDifFlx(ids_NH4B)
  RPoreSolFlx(idg_NH3,3,0,NY,NX)=RN3FL0(NY,NX)+RN3FLW0-RFLNH3-SDifFlx(idg_NH3)-RFLN3B-SDifFlx(idg_NH3B)
  RPoreSolFlx(ids_NO3,3,0,NY,NX)=RNOFL0(NY,NX)+RNOFLW0-RFLNO3-SDifFlx(ids_NO3)-RFLNOB-SDifFlx(ids_NO3B)
  RPoreSolFlx(ids_NO2,3,0,NY,NX)=RNXFL0(NY,NX)-RFLNO2-SDifFlx(ids_NO2)-RFLN2B-SDifFlx(ids_NO2B)
  RPoreSolFlx(ids_H1PO4,3,0,NY,NX)=RH1PF0(NY,NX)+RH1PFS0-RFLP14-SDifFlx(ids_H1PO4)-RFLP1B-SDifFlx(ids_H1PO4B)
  RPoreSolFlx(ids_H2PO4,3,0,NY,NX)=RH2PF0(NY,NX)+RH2PFS0-RFLPO4-SDifFlx(ids_H2PO4)-RFLPOB-SDifFlx(ids_H2PO4B)
  RPoreSolFlx(idg_CO2,3,NU(NY,NX),NY,NX)=RCOFL1(NY,NX)+RCOFLS1+RFLCOS+SDifFlx(idg_CO2)

  RPoreSolFlx(idg_CH4,3,NU(NY,NX),NY,NX)=RCHFL1(NY,NX)+RCHFLS1+RFLCHS+SDifFlx(idg_CH4)
  RPoreSolFlx(idg_O2,3,NU(NY,NX),NY,NX)=ROXFL1(NY,NX)+ROXFLS1+RFLOXS+SDifFlx(idg_O2)
  RPoreSolFlx(idg_N2,3,NU(NY,NX),NY,NX)=RNGFL1(NY,NX)+RNGFLS1+RFLNGS+SDifFlx(idg_N2)
  RPoreSolFlx(idg_N2O,3,NU(NY,NX),NY,NX)=RN2FL1(NY,NX)+RN2FLS1+RFLN2S+SDifFlx(idg_N2O)
  RPoreSolFlx(idg_H2,3,NU(NY,NX),NY,NX)=RHGFL1(NY,NX)+RFLHGS+SDifFlx(idg_H2)
  RPoreSolFlx(ids_NH4,3,NU(NY,NX),NY,NX)=RN4FL1(NY,NX)+RN4FLW1+RFLNH4+SDifFlx(ids_NH4)
  RPoreSolFlx(idg_NH3,3,NU(NY,NX),NY,NX)=RN3FL1(NY,NX)+RN3FLW1+RFLNH3+SDifFlx(idg_NH3)
  RPoreSolFlx(ids_NO3,3,NU(NY,NX),NY,NX)=RNOFL1(NY,NX)+RNOFLW1+RFLNO3+SDifFlx(ids_NO3)
  RPoreSolFlx(ids_NO2,3,NU(NY,NX),NY,NX)=RNXFL1(NY,NX)+RFLNO2+SDifFlx(ids_NO2)
  RPoreSolFlx(ids_H1PO4,3,NU(NY,NX),NY,NX)=RH1PF1(NY,NX)+RH1PFS1+RFLP14+SDifFlx(ids_H1PO4)
  RPoreSolFlx(ids_H2PO4,3,NU(NY,NX),NY,NX)=RH2PF1(NY,NX)+RH2PFS1+RFLPO4+SDifFlx(ids_H2PO4)
  RPoreSolFlx(ids_NH4B,3,NU(NY,NX),NY,NX)=RN4FL2(NY,NX)+RN4FLB1+RFLN4B+SDifFlx(ids_NH4B)
  RPoreSolFlx(idg_NH3B,3,NU(NY,NX),NY,NX)=RN3FL2(NY,NX)+RN3FLB1+RFLN3B+SDifFlx(idg_NH3B)
  RPoreSolFlx(ids_NO3B,3,NU(NY,NX),NY,NX)=RNOFL2(NY,NX)+RNOFLB1+RFLNOB+SDifFlx(ids_NO3B)
  RPoreSolFlx(ids_NO2B,3,NU(NY,NX),NY,NX)=RNXFL2(NY,NX)+RFLN2B+SDifFlx(ids_NO2B)
  RPoreSolFlx(ids_H1PO4B,3,NU(NY,NX),NY,NX)=RH1BF2(NY,NX)+RH1BFB1+RFLP1B+SDifFlx(ids_H1PO4B)
  RPoreSolFlx(ids_H2PO4B,3,NU(NY,NX),NY,NX)=RH2BF2(NY,NX)+RH2BFB1+RFLPOB+SDifFlx(ids_H2PO4B)
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
  XCOFLS(3,0,NY,NX)=XCOFLS(3,0,NY,NX)+RCOFLS0-RFLCOS-SDifFlx(idg_CO2)
  XCHFLS(3,0,NY,NX)=XCHFLS(3,0,NY,NX)+RCHFLS0-RFLCHS-SDifFlx(idg_CH4)
  XOXFLS(3,0,NY,NX)=XOXFLS(3,0,NY,NX)+ROXFLS0-RFLOXS-SDifFlx(idg_O2)
  XNGFLS(3,0,NY,NX)=XNGFLS(3,0,NY,NX)+RNGFLS0-RFLNGS-SDifFlx(idg_N2)
  XN2FLS(3,0,NY,NX)=XN2FLS(3,0,NY,NX)+RN2FLS0-RFLN2S-SDifFlx(idg_N2O)
  XHGFLS(3,0,NY,NX)=XHGFLS(3,0,NY,NX)-RFLHGS-SDifFlx(idg_H2)
  XN4FLW(3,0,NY,NX)=XN4FLW(3,0,NY,NX)+RN4FLW0-RFLNH4-SDifFlx(ids_NH4)-RFLN4B-SDifFlx(ids_NH4B)
  XN3FLW(3,0,NY,NX)=XN3FLW(3,0,NY,NX)+RN3FLW0-RFLNH3-SDifFlx(idg_NH3)-RFLN3B-SDifFlx(idg_NH3B)
  XNOFLW(3,0,NY,NX)=XNOFLW(3,0,NY,NX)+RNOFLW0-RFLNO3-SDifFlx(ids_NO3)-RFLNOB-SDifFlx(ids_NO3B)
  XNXFLS(3,0,NY,NX)=XNXFLS(3,0,NY,NX)-RFLNO2-SDifFlx(ids_NO2)-RFLN2B-SDifFlx(ids_NO2B)
  XH1PFS(3,0,NY,NX)=XH1PFS(3,0,NY,NX)+RH1PFS0-RFLP14-SDifFlx(ids_H1PO4)-RFLP1B-SDifFlx(ids_H1PO4B)
  XH2PFS(3,0,NY,NX)=XH2PFS(3,0,NY,NX)+RH2PFS0-RFLPO4-SDifFlx(ids_H2PO4)-RFLPOB-SDifFlx(ids_H2PO4B)
  XCOFLS(3,NU(NY,NX),NY,NX)=XCOFLS(3,NU(NY,NX),NY,NX)+RCOFLS1+RFLCOS+SDifFlx(idg_CO2)
  XCHFLS(3,NU(NY,NX),NY,NX)=XCHFLS(3,NU(NY,NX),NY,NX)+RCHFLS1+RFLCHS+SDifFlx(idg_CH4)
  XOXFLS(3,NU(NY,NX),NY,NX)=XOXFLS(3,NU(NY,NX),NY,NX)+ROXFLS1+RFLOXS+SDifFlx(idg_O2)
  XNGFLS(3,NU(NY,NX),NY,NX)=XNGFLS(3,NU(NY,NX),NY,NX)+RNGFLS1+RFLNGS+SDifFlx(idg_N2)
  XN2FLS(3,NU(NY,NX),NY,NX)=XN2FLS(3,NU(NY,NX),NY,NX)+RN2FLS1+RFLN2S+SDifFlx(idg_N2O)
  XHGFLS(3,NU(NY,NX),NY,NX)=XHGFLS(3,NU(NY,NX),NY,NX)+RFLHGS+SDifFlx(idg_H2)
  XN4FLW(3,NU(NY,NX),NY,NX)=XN4FLW(3,NU(NY,NX),NY,NX)+RN4FLW1+RFLNH4+SDifFlx(ids_NH4)
  XN3FLW(3,NU(NY,NX),NY,NX)=XN3FLW(3,NU(NY,NX),NY,NX)+RN3FLW1+RFLNH3+SDifFlx(idg_NH3)
  XNOFLW(3,NU(NY,NX),NY,NX)=XNOFLW(3,NU(NY,NX),NY,NX)+RNOFLW1+RFLNO3+SDifFlx(ids_NO3)
  XNXFLS(3,NU(NY,NX),NY,NX)=XNXFLS(3,NU(NY,NX),NY,NX)+RFLNO2+SDifFlx(ids_NO2)
  XH1PFS(3,NU(NY,NX),NY,NX)=XH1PFS(3,NU(NY,NX),NY,NX)+RH1PFS1+RFLP14+SDifFlx(ids_H1PO4)
  XH2PFS(3,NU(NY,NX),NY,NX)=XH2PFS(3,NU(NY,NX),NY,NX)+RH2PFS1+RFLPO4+SDifFlx(ids_H2PO4)
  XN4FLB(3,NU(NY,NX),NY,NX)=XN4FLB(3,NU(NY,NX),NY,NX)+RN4FLB1+RFLN4B+SDifFlx(ids_NH4B)
  XN3FLB(3,NU(NY,NX),NY,NX)=XN3FLB(3,NU(NY,NX),NY,NX)+RN3FLB1+RFLN3B+SDifFlx(idg_NH3B)
  XNOFLB(3,NU(NY,NX),NY,NX)=XNOFLB(3,NU(NY,NX),NY,NX)+RNOFLB1+RFLNOB+SDifFlx(ids_NO3B)
  XNXFLB(3,NU(NY,NX),NY,NX)=XNXFLB(3,NU(NY,NX),NY,NX)+RFLN2B+SDifFlx(ids_NO2B)
  XH1BFB(3,NU(NY,NX),NY,NX)=XH1BFB(3,NU(NY,NX),NY,NX)+RH1BFB1+RFLP1B+SDifFlx(ids_H1PO4B)
  XH2BFB(3,NU(NY,NX),NY,NX)=XH2BFB(3,NU(NY,NX),NY,NX)+RH2BFB1+RFLPOB+SDifFlx(ids_H2PO4B)
  end subroutine TotalPoreFluxAdjacentCell
!------------------------------------------------------------------------------------------

  subroutine MacMicPoresTransfer(M,NY,NX)
  implicit none

  integer, intent(in) :: M, NY, NX
  integer :: K,NTS
  real(r8) :: VFLW
  real(r8) :: VOLWHS,VOLWT
  real(r8) :: SDifFlx(ids_beg:ids_end)
  real(r8) :: SAdvFlx(ids_beg:ids_end)
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
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
    IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,NU(NY,NX),NY,NX) &
      /VOLWHM(M,NU(NY,NX),NY,NX)))
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
    SAdvFlx(ids_NH4)=VFLW*AZMAX1(trc_soHml2(ids_NH4,NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3)=VFLW*AZMAX1(trc_soHml2(idg_NH3,NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3)=VFLW*AZMAX1(trc_soHml2(ids_NO3,NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2)=VFLW*AZMAX1(trc_soHml2(ids_NO2,NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4)=VFLW*AZMAX1(trc_soHml2(ids_H1PO4,NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4)=VFLW*AZMAX1(trc_soHml2(ids_H2PO4,NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NH4B)=VFLW*AZMAX1(trc_soHml2(ids_NH4B,NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3B)=VFLW*AZMAX1(trc_soHml2(idg_NH3B,NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3B)=VFLW*AZMAX1(trc_soHml2(ids_NO3B,NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2B)=VFLW*AZMAX1(trc_soHml2(ids_NO2B,NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4B)=VFLW*AZMAX1(trc_soHml2(ids_H1PO4B,NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4B)=VFLW*AZMAX1(trc_soHml2(ids_H2PO4B,NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FINHM(M,NU(NY,NX),NY,NX).LT.0.0_r8)THEN
    IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX)))
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
    SAdvFlx(ids_NH4)=VFLW*AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3)=VFLW*AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3)=VFLW*AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2)=VFLW*AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4)=VFLW*AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4)=VFLW*AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NH4B)=VFLW*AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3B)=VFLW*AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO3B)=VFLW*AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_NO2B)=VFLW*AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H1PO4B)=VFLW*AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
    SAdvFlx(ids_H2PO4B)=VFLW*AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
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
!     VOLWM,VOLWHM=micropore,macropore water volume
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
  IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLWHS=AMIN1(XFRS*VOLT(NU(NY,NX),NY,NX),VOLWHM(M,NU(NY,NX),NY,NX))
    VOLWT=VOLWM(M,NU(NY,NX),NY,NX)+VOLWHS
    DO  K=1,jcplx
      DFVOC(K)=XNPH*(AZMAX1(OQCH2(K,NU(NY,NX),NY,NX))*VOLWM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQC2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
      DFVON(K)=XNPH*(AZMAX1(OQNH2(K,NU(NY,NX),NY,NX))*VOLWM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQN2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
      DFVOP(K)=XNPH*(AZMAX1(OQPH2(K,NU(NY,NX),NY,NX))*VOLWM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQP2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
      DFVOA(K)=XNPH*(AZMAX1(OQAH2(K,NU(NY,NX),NY,NX))*VOLWM(M,NU(NY,NX),NY,NX) &
        -AZMAX1(OQA2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    ENDDO

    SDifFlx(idg_CO2)=XNPH*(AZMAX1(trc_soHml2(idg_CO2,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_CO2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    SDifFlx(idg_CH4)=XNPH*(AZMAX1(trc_soHml2(idg_CH4,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_CH4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    SDifFlx(idg_O2)=XNPH*(AZMAX1(trc_soHml2(idg_O2,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_O2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    SDifFlx(idg_N2)=XNPH*(AZMAX1(trc_soHml2(idg_N2,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_N2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    SDifFlx(idg_N2O)=XNPH*(AZMAX1(trc_soHml2(idg_N2O,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_N2O,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    SDifFlx(idg_H2)=XNPH*(AZMAX1(trc_soHml2(idg_H2,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_H2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    SDifFlx(ids_NH4)=XNPH*(AZMAX1(trc_soHml2(ids_NH4,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNH4(NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3)=XNPH*(AZMAX1(trc_soHml2(idg_NH3,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNH4(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO3)=XNPH*(AZMAX1(trc_soHml2(ids_NO3,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNO3(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO2)=XNPH*(AZMAX1(trc_soHml2(ids_NO2,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNO3(NU(NY,NX),NY,NX)
    SDifFlx(ids_H1PO4)=XNPH*(AZMAX1(trc_soHml2(ids_H1PO4,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPO4(NU(NY,NX),NY,NX)
    SDifFlx(ids_H2PO4)=XNPH*(AZMAX1(trc_soHml2(ids_H2PO4,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPO4(NU(NY,NX),NY,NX)
    SDifFlx(ids_NH4B)=XNPH*(AZMAX1(trc_soHml2(ids_NH4B,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNHB(NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3B)=XNPH*(AZMAX1(trc_soHml2(idg_NH3B,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNHB(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO3B)=XNPH*(AZMAX1(trc_soHml2(ids_NO3B,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNOB(NU(NY,NX),NY,NX)
    SDifFlx(ids_NO2B)=XNPH*(AZMAX1(trc_soHml2(ids_NO2B,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNOB(NU(NY,NX),NY,NX)
    SDifFlx(ids_H1PO4B)=XNPH*(AZMAX1(trc_soHml2(ids_H1PO4B,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPOB(NU(NY,NX),NY,NX)
    SDifFlx(ids_H2PO4B)=XNPH*(AZMAX1(trc_soHml2(ids_H2PO4B,NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX)-AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPOB(NU(NY,NX),NY,NX)
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
  XCOFXS(NU(NY,NX),NY,NX)=XCOFXS(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_CO2,NU(NY,NX),NY,NX)
  XCHFXS(NU(NY,NX),NY,NX)=XCHFXS(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_CH4,NU(NY,NX),NY,NX)
  XOXFXS(NU(NY,NX),NY,NX)=XOXFXS(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_O2,NU(NY,NX),NY,NX)
  XNGFXS(NU(NY,NX),NY,NX)=XNGFXS(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_N2,NU(NY,NX),NY,NX)
  XN2FXS(NU(NY,NX),NY,NX)=XN2FXS(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_N2O,NU(NY,NX),NY,NX)
  XHGFXS(NU(NY,NX),NY,NX)=XHGFXS(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_H2,NU(NY,NX),NY,NX)
  XN4FXW(NU(NY,NX),NY,NX)=XN4FXW(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_NH4,NU(NY,NX),NY,NX)
  XN3FXW(NU(NY,NX),NY,NX)=XN3FXW(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_NH3,NU(NY,NX),NY,NX)
  XNOFXW(NU(NY,NX),NY,NX)=XNOFXW(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_NO3,NU(NY,NX),NY,NX)
  XNXFXS(NU(NY,NX),NY,NX)=XNXFXS(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_NO2,NU(NY,NX),NY,NX)
  XH1PXS(NU(NY,NX),NY,NX)=XH1PXS(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_H1PO4,NU(NY,NX),NY,NX)
  XH2PXS(NU(NY,NX),NY,NX)=XH2PXS(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_H2PO4,NU(NY,NX),NY,NX)
  XN4FXB(NU(NY,NX),NY,NX)=XN4FXB(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_NH4B,NU(NY,NX),NY,NX)
  XN3FXB(NU(NY,NX),NY,NX)=XN3FXB(NU(NY,NX),NY,NX)+RporeSoXFlx(idg_NH3B,NU(NY,NX),NY,NX)
  XNOFXB(NU(NY,NX),NY,NX)=XNOFXB(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_NO3B,NU(NY,NX),NY,NX)
  XNXFXB(NU(NY,NX),NY,NX)=XNXFXB(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_NO2B,NU(NY,NX),NY,NX)
  XH1BXB(NU(NY,NX),NY,NX)=XH1BXB(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_H1PO4B,NU(NY,NX),NY,NX)
  XH2BXB(NU(NY,NX),NY,NX)=XH2BXB(NU(NY,NX),NY,NX)+RporeSoXFlx(ids_H2PO4B,NU(NY,NX),NY,NX)
!
  end subroutine MacMicPoresTransfer
!------------------------------------------------------------------------------------------

  subroutine OverlandFlowSnowdriftTransport(M,NY,NX,NHE,NHW,NVS,NVN)
  implicit none

  integer, intent(in) :: M, NY, NX, NHE, NHW, NVS, NVN
  real(r8) :: FQRM,VFLW
  integer :: K,N,NN,N1,N2,N4,N5,N4B,N5B
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VOLWM=litter water volume from watsub.f
!     *S2=litter solute content
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     X*QRS=accumulated hourly solute in runoff
!
  N1=NX
  N2=NY
  IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VOLWM(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,QRM(M,N2,N1)/VOLWM(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      RQROC0(K,N2,N1)=VFLW*AZMAX1(OQC2(K,0,N2,N1))
      RQRON0(K,N2,N1)=VFLW*AZMAX1(OQN2(K,0,N2,N1))
      RQROP0(K,N2,N1)=VFLW*AZMAX1(OQP2(K,0,N2,N1))
      RQROA0(K,N2,N1)=VFLW*AZMAX1(OQA2(K,0,N2,N1))
    ENDDO
    RQRCOS0(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_CO2,0,N2,N1))
    RQRCHS0(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_CH4,0,N2,N1))
    RQROXS0(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_O2,0,N2,N1))
    RQRNGS0(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_N2,0,N2,N1))
    RQRN2S0(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_N2O,0,N2,N1))
    RQRHGS0(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_H2,0,N2,N1))
    RQRNH40(N2,N1)=VFLW*AZMAX1(trc_solml2(ids_NH4,0,N2,N1))
    RQRNH30(N2,N1)=VFLW*AZMAX1(trc_solml2(idg_NH3,0,N2,N1))
    RQRNO30(N2,N1)=VFLW*AZMAX1(trc_solml2(ids_NO3,0,N2,N1))
    RQRNO20(N2,N1)=VFLW*AZMAX1(trc_solml2(ids_NO2,0,N2,N1))
    RQRH1P0(N2,N1)=VFLW*AZMAX1(trc_solml2(ids_H1PO4,0,N2,N1))
    RQRH2P0(N2,N1)=VFLW*AZMAX1(trc_solml2(ids_H2PO4,0,N2,N1))
  ELSE
    DO K=1,jcplx
      RQROC0(K,N2,N1)=0.0_r8
      RQRON0(K,N2,N1)=0.0_r8
      RQROP0(K,N2,N1)=0.0_r8
      RQROA0(K,N2,N1)=0.0_r8
    ENDDO
    RQRCOS0(N2,N1)=0.0_r8
    RQRCHS0(N2,N1)=0.0_r8
    RQROXS0(N2,N1)=0.0_r8
    RQRNGS0(N2,N1)=0.0_r8
    RQRN2S0(N2,N1)=0.0_r8
    RQRHGS0(N2,N1)=0.0_r8
    RQRNH40(N2,N1)=0.0_r8
    RQRNH30(N2,N1)=0.0_r8
    RQRNO30(N2,N1)=0.0_r8
    RQRNO20(N2,N1)=0.0_r8
    RQRH1P0(N2,N1)=0.0_r8
    RQRH2P0(N2,N1)=0.0_r8
  ENDIF
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO 4310 N=1,2
    DO 4305 NN=1,2
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
!     VOLWM=litter water volume from watsub.f
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
          RQRCOS(N,2,N5,N4)=RQRCOS0(N2,N1)*FQRM
          RQRCHS(N,2,N5,N4)=RQRCHS0(N2,N1)*FQRM
          RQROXS(N,2,N5,N4)=RQROXS0(N2,N1)*FQRM
          RQRNGS(N,2,N5,N4)=RQRNGS0(N2,N1)*FQRM
          RQRN2S(N,2,N5,N4)=RQRN2S0(N2,N1)*FQRM
          RQRHGS(N,2,N5,N4)=RQRHGS0(N2,N1)*FQRM
          RQRNH4(N,2,N5,N4)=RQRNH40(N2,N1)*FQRM
          RQRNH3(N,2,N5,N4)=RQRNH30(N2,N1)*FQRM
          RQRNO3(N,2,N5,N4)=RQRNO30(N2,N1)*FQRM
          RQRNO2(N,2,N5,N4)=RQRNO20(N2,N1)*FQRM
          RQRH1P(N,2,N5,N4)=RQRH1P0(N2,N1)*FQRM
          RQRH2P(N,2,N5,N4)=RQRH2P0(N2,N1)*FQRM
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
          XCOQRS(N,2,N5,N4)=XCOQRS(N,2,N5,N4)+RQRCOS(N,2,N5,N4)
          XCHQRS(N,2,N5,N4)=XCHQRS(N,2,N5,N4)+RQRCHS(N,2,N5,N4)
          XOXQRS(N,2,N5,N4)=XOXQRS(N,2,N5,N4)+RQROXS(N,2,N5,N4)
          XNGQRS(N,2,N5,N4)=XNGQRS(N,2,N5,N4)+RQRNGS(N,2,N5,N4)
          XN2QRS(N,2,N5,N4)=XN2QRS(N,2,N5,N4)+RQRN2S(N,2,N5,N4)
          XHGQRS(N,2,N5,N4)=XHGQRS(N,2,N5,N4)+RQRHGS(N,2,N5,N4)
          XN4QRW(N,2,N5,N4)=XN4QRW(N,2,N5,N4)+RQRNH4(N,2,N5,N4)
          XN3QRW(N,2,N5,N4)=XN3QRW(N,2,N5,N4)+RQRNH3(N,2,N5,N4)
          XNOQRW(N,2,N5,N4)=XNOQRW(N,2,N5,N4)+RQRNO3(N,2,N5,N4)
          XNXQRS(N,2,N5,N4)=XNXQRS(N,2,N5,N4)+RQRNO2(N,2,N5,N4)
          XP1QRW(N,2,N5,N4)=XP1QRW(N,2,N5,N4)+RQRH1P(N,2,N5,N4)
          XP4QRW(N,2,N5,N4)=XP4QRW(N,2,N5,N4)+RQRH2P(N,2,N5,N4)
        ELSE
          DO K=1,jcplx
            RQROC(K,N,2,N5,N4)=0.0_r8
            RQRON(K,N,2,N5,N4)=0.0_r8
            RQROP(K,N,2,N5,N4)=0.0_r8
            RQROA(K,N,2,N5,N4)=0.0_r8
          ENDDO
          RQRCOS(N,2,N5,N4)=0.0_r8
          RQRCHS(N,2,N5,N4)=0.0_r8
          RQROXS(N,2,N5,N4)=0.0_r8
          RQRNGS(N,2,N5,N4)=0.0_r8
          RQRN2S(N,2,N5,N4)=0.0_r8
          RQRHGS(N,2,N5,N4)=0.0_r8
          RQRNH4(N,2,N5,N4)=0.0_r8
          RQRNH3(N,2,N5,N4)=0.0_r8
          RQRNO3(N,2,N5,N4)=0.0_r8
          RQRNO2(N,2,N5,N4)=0.0_r8
          RQRH1P(N,2,N5,N4)=0.0_r8
          RQRH2P(N,2,N5,N4)=0.0_r8
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
            RQRCOS(N,1,N5B,N4B)=RQRCOS0(N2,N1)*FQRM
            RQRCHS(N,1,N5B,N4B)=RQRCHS0(N2,N1)*FQRM
            RQROXS(N,1,N5B,N4B)=RQROXS0(N2,N1)*FQRM
            RQRNGS(N,1,N5B,N4B)=RQRNGS0(N2,N1)*FQRM
            RQRN2S(N,1,N5B,N4B)=RQRN2S0(N2,N1)*FQRM
            RQRHGS(N,1,N5B,N4B)=RQRHGS0(N2,N1)*FQRM
            RQRNH4(N,1,N5B,N4B)=RQRNH40(N2,N1)*FQRM
            RQRNH3(N,1,N5B,N4B)=RQRNH30(N2,N1)*FQRM
            RQRNO3(N,1,N5B,N4B)=RQRNO30(N2,N1)*FQRM
            RQRNO2(N,1,N5B,N4B)=RQRNO20(N2,N1)*FQRM
            RQRH1P(N,1,N5B,N4B)=RQRH1P0(N2,N1)*FQRM
            RQRH2P(N,1,N5B,N4B)=RQRH2P0(N2,N1)*FQRM
            DO K=1,jcplx
              XOCQRS(K,N,1,N5B,N4B)=XOCQRS(K,N,1,N5B,N4B)+RQROC(K,N,1,N5B,N4B)
              XONQRS(K,N,1,N5B,N4B)=XONQRS(K,N,1,N5B,N4B)+RQRON(K,N,1,N5B,N4B)
              XOPQRS(K,N,1,N5B,N4B)=XOPQRS(K,N,1,N5B,N4B)+RQROP(K,N,1,N5B,N4B)
              XOAQRS(K,N,1,N5B,N4B)=XOAQRS(K,N,1,N5B,N4B)+RQROA(K,N,1,N5B,N4B)
            ENDDO
            XCOQRS(N,1,N5B,N4B)=XCOQRS(N,1,N5B,N4B)+RQRCOS(N,1,N5B,N4B)
            XCHQRS(N,1,N5B,N4B)=XCHQRS(N,1,N5B,N4B)+RQRCHS(N,1,N5B,N4B)
            XOXQRS(N,1,N5B,N4B)=XOXQRS(N,1,N5B,N4B)+RQROXS(N,1,N5B,N4B)
            XNGQRS(N,1,N5B,N4B)=XNGQRS(N,1,N5B,N4B)+RQRNGS(N,1,N5B,N4B)
            XN2QRS(N,1,N5B,N4B)=XN2QRS(N,1,N5B,N4B)+RQRN2S(N,1,N5B,N4B)
            XHGQRS(N,1,N5B,N4B)=XHGQRS(N,1,N5B,N4B)+RQRHGS(N,1,N5B,N4B)
            XN4QRW(N,1,N5B,N4B)=XN4QRW(N,1,N5B,N4B)+RQRNH4(N,1,N5B,N4B)
            XN3QRW(N,1,N5B,N4B)=XN3QRW(N,1,N5B,N4B)+RQRNH3(N,1,N5B,N4B)
            XNOQRW(N,1,N5B,N4B)=XNOQRW(N,1,N5B,N4B)+RQRNO3(N,1,N5B,N4B)
            XNXQRS(N,1,N5B,N4B)=XNXQRS(N,1,N5B,N4B)+RQRNO2(N,1,N5B,N4B)
            XP1QRW(N,1,N5B,N4B)=XP1QRW(N,1,N5B,N4B)+RQRH1P(N,1,N5B,N4B)
            XP4QRW(N,1,N5B,N4B)=XP4QRW(N,1,N5B,N4B)+RQRH2P(N,1,N5B,N4B)
          ELSE
            DO  K=1,jcplx
              RQROC(K,N,1,N5B,N4B)=0.0_r8
              RQRON(K,N,1,N5B,N4B)=0.0_r8
              RQROP(K,N,1,N5B,N4B)=0.0_r8
              RQROA(K,N,1,N5B,N4B)=0.0_r8
            ENDDO
            RQRCOS(N,1,N5B,N4B)=0.0_r8
            RQRCHS(N,1,N5B,N4B)=0.0_r8
            RQROXS(N,1,N5B,N4B)=0.0_r8
            RQRNGS(N,1,N5B,N4B)=0.0_r8
            RQRN2S(N,1,N5B,N4B)=0.0_r8
            RQRHGS(N,1,N5B,N4B)=0.0_r8
            RQRNH4(N,1,N5B,N4B)=0.0_r8
            RQRNH3(N,1,N5B,N4B)=0.0_r8
            RQRNO3(N,1,N5B,N4B)=0.0_r8
            RQRNO2(N,1,N5B,N4B)=0.0_r8
            RQRH1P(N,1,N5B,N4B)=0.0_r8
            RQRH2P(N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        DO K=1,jcplx
          RQROC(K,N,2,N5,N4)=0.0_r8
          RQRON(K,N,2,N5,N4)=0.0_r8
          RQROP(K,N,2,N5,N4)=0.0_r8
          RQROA(K,N,2,N5,N4)=0.0_r8
        ENDDO
        RQRCOS(N,2,N5,N4)=0.0_r8
        RQRCHS(N,2,N5,N4)=0.0_r8
        RQROXS(N,2,N5,N4)=0.0_r8
        RQRNGS(N,2,N5,N4)=0.0_r8
        RQRN2S(N,2,N5,N4)=0.0_r8
        RQRHGS(N,2,N5,N4)=0.0_r8
        RQRNH4(N,2,N5,N4)=0.0_r8
        RQRNH3(N,2,N5,N4)=0.0_r8
        RQRNO3(N,2,N5,N4)=0.0_r8
        RQRNO2(N,2,N5,N4)=0.0_r8
        RQRH1P(N,2,N5,N4)=0.0_r8
        RQRH2P(N,2,N5,N4)=0.0_r8
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          DO  K=1,jcplx
            RQROC(K,N,1,N5B,N4B)=0.0_r8
            RQRON(K,N,1,N5B,N4B)=0.0_r8
            RQROP(K,N,1,N5B,N4B)=0.0_r8
            RQROA(K,N,1,N5B,N4B)=0.0_r8
          ENDDO
          RQRCOS(N,1,N5B,N4B)=0.0_r8
          RQRCHS(N,1,N5B,N4B)=0.0_r8
          RQROXS(N,1,N5B,N4B)=0.0_r8
          RQRNGS(N,1,N5B,N4B)=0.0_r8
          RQRN2S(N,1,N5B,N4B)=0.0_r8
          RQRHGS(N,1,N5B,N4B)=0.0_r8
          RQRNH4(N,1,N5B,N4B)=0.0_r8
          RQRNH3(N,1,N5B,N4B)=0.0_r8
          RQRNO3(N,1,N5B,N4B)=0.0_r8
          RQRNO2(N,1,N5B,N4B)=0.0_r8
          RQRH1P(N,1,N5B,N4B)=0.0_r8
          RQRH2P(N,1,N5B,N4B)=0.0_r8
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
          RQSCOS(N,N5,N4)=0.0_r8
          RQSCHS(N,N5,N4)=0.0_r8
          RQSOXS(N,N5,N4)=0.0_r8
          RQSNGS(N,N5,N4)=0.0_r8
          RQSN2S(N,N5,N4)=0.0_r8
          RQSNH4(N,N5,N4)=0.0_r8
          RQSNH3(N,N5,N4)=0.0_r8
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
          RQSCOS(N,N5,N4)=VFLW*AZMAX1(CO2W2(1,N2,N1))
          RQSCHS(N,N5,N4)=VFLW*AZMAX1(CH4W2(1,N2,N1))
          RQSOXS(N,N5,N4)=VFLW*AZMAX1(OXYW2(1,N2,N1))
          RQSNGS(N,N5,N4)=VFLW*AZMAX1(ZNGW2(1,N2,N1))
          RQSN2S(N,N5,N4)=VFLW*AZMAX1(ZN2W2(1,N2,N1))
          RQSNH4(N,N5,N4)=VFLW*AZMAX1(ZN4W2(1,N2,N1))
          RQSNH3(N,N5,N4)=VFLW*AZMAX1(ZN3W2(1,N2,N1))
          RQSNO3(N,N5,N4)=VFLW*AZMAX1(ZNOW2(1,N2,N1))
          RQSH1P(N,N5,N4)=VFLW*AZMAX1(Z1PW2(1,N2,N1))
          RQSH2P(N,N5,N4)=VFLW*AZMAX1(ZHPW2(1,N2,N1))
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
          IF(VOLS(N5,N4).GT.ZEROS2(N5,N4))THEN
            VFLW=AZMIN1(AMAX1(-VFLWX,QSM(M,N,N5,N4)/VOLS(N5,N4)))
          ELSE
            VFLW=-VFLWX
          ENDIF
          RQSCOS(N,N5,N4)=VFLW*AZMAX1(CO2W2(1,N5,N4))
          RQSCHS(N,N5,N4)=VFLW*AZMAX1(CH4W2(1,N5,N4))
          RQSOXS(N,N5,N4)=VFLW*AZMAX1(OXYW2(1,N5,N4))
          RQSNGS(N,N5,N4)=VFLW*AZMAX1(ZNGW2(1,N5,N4))
          RQSN2S(N,N5,N4)=VFLW*AZMAX1(ZN2W2(1,N5,N4))
          RQSNH4(N,N5,N4)=VFLW*AZMAX1(ZN4W2(1,N5,N4))
          RQSNH3(N,N5,N4)=VFLW*AZMAX1(ZN3W2(1,N5,N4))
          RQSNO3(N,N5,N4)=VFLW*AZMAX1(ZNOW2(1,N5,N4))
          RQSH1P(N,N5,N4)=VFLW*AZMAX1(Z1PW2(1,N5,N4))
          RQSH2P(N,N5,N4)=VFLW*AZMAX1(ZHPW2(1,N5,N4))
        ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QSS=hourly solute in snow flux
!     RQS*=solute in snow flux
!
        XCOQSS(N,N5,N4)=XCOQSS(N,N5,N4)+RQSCOS(N,N5,N4)
        XCHQSS(N,N5,N4)=XCHQSS(N,N5,N4)+RQSCHS(N,N5,N4)
        XOXQSS(N,N5,N4)=XOXQSS(N,N5,N4)+RQSOXS(N,N5,N4)
        XNGQSS(N,N5,N4)=XNGQSS(N,N5,N4)+RQSNGS(N,N5,N4)
        XN2QSS(N,N5,N4)=XN2QSS(N,N5,N4)+RQSN2S(N,N5,N4)
        XN4QSS(N,N5,N4)=XN4QSS(N,N5,N4)+RQSNH4(N,N5,N4)
        XN3QSS(N,N5,N4)=XN3QSS(N,N5,N4)+RQSNH3(N,N5,N4)
        XNOQSS(N,N5,N4)=XNOQSS(N,N5,N4)+RQSNO3(N,N5,N4)
        XP1QSS(N,N5,N4)=XP1QSS(N,N5,N4)+RQSH1P(N,N5,N4)
        XP4QSS(N,N5,N4)=XP4QSS(N,N5,N4)+RQSH2P(N,N5,N4)
      ENDIF
4305  CONTINUE
4310  CONTINUE
  end subroutine OverlandFlowSnowdriftTransport
!------------------------------------------------------------------------------------------

  subroutine LitterGasVolatilDissol(M,NY,NX)
  implicit none

  integer, intent(in) :: M,NY, NX
  real(r8) :: VOLHGR(JY,JX),VOLN3R(JY,JX)
  real(r8) :: VOLN2R(JY,JX),VOLOXR(JY,JX)
  real(r8) :: VOLCOR(JY,JX),VOLCHR(JY,JX)
  real(r8) :: VOLNGR(JY,JX)
  real(r8) :: CNH3S0,CNH4S0,CO2G0,CH4G0
  real(r8) :: Z2GG0,Z2OG0,ZN3G0,H2GG0,OXYG0
!
!     VOLT=litter volume from hour1.f
!     VOLWM,VOLPM=micropore water volume, air volume from watsub.f
!     VOLW*=equivalent aqueous volume for gas
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
  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VOLPM(M,0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN

    VOLWCO(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_CO2,0,NY,NX)
    VOLWCH(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_CH4,0,NY,NX)
    VOLWOX(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_O2,0,NY,NX)
    VOLWNG(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_N2,0,NY,NX)
    VOLWN2(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_N2O,0,NY,NX)
    VOLWN3(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_NH3,0,NY,NX)
    VOLWHG(0,NY,NX)=VOLWM(M,0,NY,NX)*GSolbility(idg_H2,0,NY,NX)

    VOLWXA(0,NY,NX)=natomw*VOLWM(M,0,NY,NX)
    CO2G0=CCO2G(0,NY,NX)*VOLPM(M,0,NY,NX)
    CH4G0=CCH4G(0,NY,NX)*VOLPM(M,0,NY,NX)
    OXYG0=COXYG(0,NY,NX)*VOLPM(M,0,NY,NX)
    Z2GG0=CZ2GG(0,NY,NX)*VOLPM(M,0,NY,NX)
    Z2OG0=CZ2OG(0,NY,NX)*VOLPM(M,0,NY,NX)
    ZN3G0=CNH3G(0,NY,NX)*VOLPM(M,0,NY,NX)
    H2GG0=CH2GG(0,NY,NX)*VOLPM(M,0,NY,NX)
    VOLCOR(NY,NX)=VOLWCO(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLCHR(NY,NX)=VOLWCH(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLOXR(NY,NX)=VOLWOX(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLNGR(NY,NX)=VOLWNG(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLN2R(NY,NX)=VOLWN2(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLN3R(NY,NX)=VOLWN3(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLHGR(NY,NX)=VOLWHG(0,NY,NX)+VOLPM(M,0,NY,NX)

    RGasDSFlx(idg_CO2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),CO2G0)*VOLWCO(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_CO2,0,NY,NX)+RCODXR)*VOLPM(M,0,NY,NX))/VOLCOR(NY,NX)
    RGasDSFlx(idg_CH4,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),CH4G0)*VOLWCH(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_CH4,0,NY,NX)+RCHDXR)*VOLPM(M,0,NY,NX))/VOLCHR(NY,NX)
    RGasDSFlx(idg_O2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),OXYG0)*VOLWOX(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_O2,0,NY,NX)+ROXDXR)*VOLPM(M,0,NY,NX))/VOLOXR(NY,NX)
    RGasDSFlx(idg_N2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),Z2GG0)*VOLWNG(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_N2,0,NY,NX)+RNGDXR)*VOLPM(M,0,NY,NX))/VOLNGR(NY,NX)
    RGasDSFlx(idg_N2O,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),Z2OG0)*VOLWN2(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_N2O,0,NY,NX)+RN2DXR)*VOLPM(M,0,NY,NX))/VOLN2R(NY,NX)
    RGasDSFlx(idg_NH3,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),ZN3G0)*VOLWN3(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_NH3,0,NY,NX)+RN3DXR)*VOLPM(M,0,NY,NX))/VOLN3R(NY,NX)
    CNH3S0=AZMAX1((trc_solml2(idg_NH3,0,NY,NX)+RGasDSFlx(idg_NH3,0,NY,NX)))/VOLWXA(0,NY,NX)
    CNH4S0=AZMAX1(trc_solml2(ids_NH4,0,NY,NX))/VOLWXA(0,NY,NX)
    RGasDSFlx(idg_H2,0,NY,NX)=DFGS(M,0,NY,NX)*(AMAX1(ZEROS(NY,NX),H2GG0)*VOLWHG(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2(idg_H2,0,NY,NX)+RHGDXR)*VOLPM(M,0,NY,NX))/VOLHGR(NY,NX)

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly surface gas volatilization
!     R*DFG=surface gas volatilization
!
    XCODFG(0,NY,NX)=XCODFG(0,NY,NX)+RGasDSFlx(idg_CO2,0,NY,NX)
    XCHDFG(0,NY,NX)=XCHDFG(0,NY,NX)+RGasDSFlx(idg_CH4,0,NY,NX)
    XOXDFG(0,NY,NX)=XOXDFG(0,NY,NX)+RGasDSFlx(idg_O2,0,NY,NX)
    XNGDFG(0,NY,NX)=XNGDFG(0,NY,NX)+RGasDSFlx(idg_N2,0,NY,NX)
    XN2DFG(0,NY,NX)=XN2DFG(0,NY,NX)+RGasDSFlx(idg_N2O,0,NY,NX)
    XN3DFG(0,NY,NX)=XN3DFG(0,NY,NX)+RGasDSFlx(idg_NH3,0,NY,NX)
    XHGDFG(0,NY,NX)=XHGDFG(0,NY,NX)+RGasDSFlx(idg_H2,0,NY,NX)
  ELSE
    RGasDSFlx(idg_CO2,0,NY,NX)=0.0_r8
    RGasDSFlx(idg_CH4,0,NY,NX)=0.0_r8
    RGasDSFlx(idg_O2,0,NY,NX)=0.0_r8
    RGasDSFlx(idg_N2,0,NY,NX)=0.0_r8
    RGasDSFlx(idg_N2O,0,NY,NX)=0.0_r8
    RGasDSFlx(idg_NH3,0,NY,NX)=0.0_r8
    RGasDSFlx(idg_H2,0,NY,NX)=0.0_r8
  ENDIF
  end subroutine LitterGasVolatilDissol

!------------------------------------------------------------------------------------------

  subroutine SurfSoilFluxGasDifAdv(M,NY,NX,FLQM)
!
! DESCRIPTION:
! surface soil gaseous diffusion, advection, dissolution & volatilization
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8),intent(in) :: FLQM(3,JD,JV,JH)
  real(r8) :: VFLW
  real(r8) :: DFVHGG,DFVN3G,DFVN2G
  real(r8) :: DNH3GQ,DH2GGQ
!
!     THETPM=air-filled porosity from watsub.f
!     BKDS=bulk density
!
  IF(THETPM(M,NU(NY,NX),NY,NX).GT.THETX.AND.BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
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
    XCOFLG(3,NU(NY,NX),NY,NX)=XCOFLG(3,NU(NY,NX),NY,NX)+RCOFLG(3,NU(NY,NX),NY,NX)
    XCHFLG(3,NU(NY,NX),NY,NX)=XCHFLG(3,NU(NY,NX),NY,NX)+RCHFLG(3,NU(NY,NX),NY,NX)
    XOXFLG(3,NU(NY,NX),NY,NX)=XOXFLG(3,NU(NY,NX),NY,NX)+ROXFLG(3,NU(NY,NX),NY,NX)
    XNGFLG(3,NU(NY,NX),NY,NX)=XNGFLG(3,NU(NY,NX),NY,NX)+RNGFLG(3,NU(NY,NX),NY,NX)
    XN2FLG(3,NU(NY,NX),NY,NX)=XN2FLG(3,NU(NY,NX),NY,NX)+RN2FLG(3,NU(NY,NX),NY,NX)
    XN3FLG(3,NU(NY,NX),NY,NX)=XN3FLG(3,NU(NY,NX),NY,NX)+RN3FLG(3,NU(NY,NX),NY,NX)
    XHGFLG(3,NU(NY,NX),NY,NX)=XHGFLG(3,NU(NY,NX),NY,NX)+RHGFLG(3,NU(NY,NX),NY,NX)

    IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      call SurfSoilFluxDisolVapor(M,NY,NX)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      XCODFG(NU(NY,NX),NY,NX)=XCODFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_CO2,NU(NY,NX),NY,NX)
      XCHDFG(NU(NY,NX),NY,NX)=XCHDFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_CH4,NU(NY,NX),NY,NX)
      XOXDFG(NU(NY,NX),NY,NX)=XOXDFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_O2,NU(NY,NX),NY,NX)
      XNGDFG(NU(NY,NX),NY,NX)=XNGDFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_N2,NU(NY,NX),NY,NX)
      XN2DFG(NU(NY,NX),NY,NX)=XN2DFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_N2O,NU(NY,NX),NY,NX)
      XN3DFG(NU(NY,NX),NY,NX)=XN3DFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)
      XHGDFG(NU(NY,NX),NY,NX)=XHGDFG(NU(NY,NX),NY,NX)+RGasDSFlx(idg_H2,NU(NY,NX),NY,NX)
      XNBDFG(NU(NY,NX),NY,NX)=XNBDFG(NU(NY,NX),NY,NX)+RNBDFG(NU(NY,NX),NY,NX)
    ELSE
      RGasDSFlx(idg_CO2,NU(NY,NX),NY,NX)=0.0_r8
      RGasDSFlx(idg_CH4,NU(NY,NX),NY,NX)=0.0_r8
      RGasDSFlx(idg_O2,NU(NY,NX),NY,NX)=0.0_r8
      RGasDSFlx(idg_N2,NU(NY,NX),NY,NX)=0.0_r8
      RGasDSFlx(idg_N2O,NU(NY,NX),NY,NX)=0.0_r8
      RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)=0.0_r8
      RGasDSFlx(idg_H2,NU(NY,NX),NY,NX)=0.0_r8
      RNBDFG(NU(NY,NX),NY,NX)=0.0_r8
    ENDIF
  ELSE
    RCOFLG(3,NU(NY,NX),NY,NX)=0.0_r8
    RCHFLG(3,NU(NY,NX),NY,NX)=0.0_r8
    ROXFLG(3,NU(NY,NX),NY,NX)=0.0_r8
    RNGFLG(3,NU(NY,NX),NY,NX)=0.0_r8
    RN2FLG(3,NU(NY,NX),NY,NX)=0.0_r8
    RN3FLG(3,NU(NY,NX),NY,NX)=0.0_r8
    RHGFLG(3,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_CO2,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_CH4,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_O2,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_N2O,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_N2,NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)=0.0_r8
    RNBDFG(NU(NY,NX),NY,NX)=0.0_r8
    RGasDSFlx(idg_H2,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  end subroutine SurfSoilFluxGasDifAdv

!------------------------------------------------------------------------------------------

  subroutine LitterAtmosExchange(M,NY,NX)
  implicit none
  integer, intent(in) :: NY,NX,M
  integer :: K

  real(r8) :: CCH4GQ,CCO2GQ,CH2GGQ,COXYGQ,CZ2GGQ,CZ2OGQ,CZN3GQ
  real(r8) :: DFGSCH,DFGSCO,DFGSHL,DFGSN2,DFGSN3,DFGSNG,DFGSOX
  real(r8) :: DLYR0,TORT0

  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX).AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX)) !vertical layer thickness
    TORT0=TORT(M,0,NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR0)*CVRD(NY,NX)

    DFGSCO=SolDifcc(idg_CO2,0,NY,NX)*TORT0
    DFGSCH=SolDifcc(idg_CH4,0,NY,NX)*TORT0
    DFGSOX=SolDifcc(idg_O2,0,NY,NX)*TORT0
    DFGSNG=SolDifcc(idg_N2,0,NY,NX)*TORT0
    DFGSN2=SolDifcc(idg_NH3,0,NY,NX)*TORT0
    DFGSN3=SolDifcc(idg_N2O,0,NY,NX)*TORT0
    DFGSHL=SolDifcc(idg_H2,0,NY,NX)*TORT0
    DO  K=1,jcplx
      COQC1(K)=AZMAX1(OQC2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
      COQN1(K)=AZMAX1(OQN2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
      COQP1(K)=AZMAX1(OQP2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
      COQA1(K)=AZMAX1(OQA2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
    ENDDO


    CCO2S1=AZMAX1(trc_solml2(idg_CO2,0,NY,NX)/VOLWM(M,0,NY,NX))
    CCH4S1=AZMAX1(trc_solml2(idg_CH4,0,NY,NX)/VOLWM(M,0,NY,NX))
    COXYS1=AZMAX1(trc_solml2(idg_O2,0,NY,NX)/VOLWM(M,0,NY,NX))
    CZ2GS1=AZMAX1(trc_solml2(idg_N2,0,NY,NX)/VOLWM(M,0,NY,NX))
    CZ2OS1=AZMAX1(trc_solml2(idg_N2O,0,NY,NX)/VOLWM(M,0,NY,NX))
    CH2GS1=AZMAX1(trc_solml2(idg_H2,0,NY,NX)/VOLWM(M,0,NY,NX))
    CNH4S1=AZMAX1(trc_solml2(ids_NH4,0,NY,NX)/VOLWM(M,0,NY,NX))
    CNH3S1=AZMAX1(trc_solml2(idg_NH3,0,NY,NX)/VOLWM(M,0,NY,NX))
    CNO3S1=AZMAX1(trc_solml2(ids_NO3,0,NY,NX)/VOLWM(M,0,NY,NX))
    CNO2S1=AZMAX1(trc_solml2(ids_NO2,0,NY,NX)/VOLWM(M,0,NY,NX))
    CP14S1=AZMAX1(trc_solml2(ids_H1PO4,0,NY,NX)/VOLWM(M,0,NY,NX))
    CPO4S1=AZMAX1(trc_solml2(ids_H2PO4,0,NY,NX)/VOLWM(M,0,NY,NX))
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
    CCO2GQ=(PARR(NY,NX)*AtmGgms(idg_CO2,NY,NX)*GSolbility(idg_CO2,0,NY,NX)+DFGSCO*CCO2S1)/(DFGSCO+PARR(NY,NX))
    CCH4GQ=(PARR(NY,NX)*AtmGgms(idg_CH4,NY,NX)*GSolbility(idg_CH4,0,NY,NX)+DFGSCH*CCH4S1)/(DFGSCH+PARR(NY,NX))
    COXYGQ=(PARR(NY,NX)*AtmGgms(idg_O2,NY,NX)*GSolbility(idg_O2,0,NY,NX)+DFGSOX*COXYS1)/(DFGSOX+PARR(NY,NX))
    CZ2GGQ=(PARR(NY,NX)*AtmGgms(idg_N2,NY,NX)*GSolbility(idg_N2,0,NY,NX)+DFGSNG*CZ2GS1)/(DFGSNG+PARR(NY,NX))
    CZ2OGQ=(PARR(NY,NX)*AtmGgms(idg_N2O,NY,NX)*GSolbility(idg_N2O,0,NY,NX)+DFGSN2*CZ2OS1)/(DFGSN2+PARR(NY,NX))
    CZN3GQ=(PARR(NY,NX)*AtmGgms(idg_NH3,NY,NX)*GSolbility(idg_NH3,0,NY,NX)+DFGSN3*CNH3S1)/(DFGSN3+PARR(NY,NX))
    CH2GGQ=(PARR(NY,NX)*AtmGgms(idg_H2,NY,NX)*GSolbility(idg_H2,0,NY,NX)+DFGSHL*CH2GS1)/(DFGSHL+PARR(NY,NX))
!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND RESIDUE SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!
!     R*DFR,X*DFR=current,hourly gas exchange between atmosphere and surface litter water
!     gas code: *CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at litter surface
!     VOLWM=litter water volume
!     DFGS*=effective solute diffusivity
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     R*DXR=R*DFR for gas flux calculations
!
    RCODFR(NY,NX)=(CCO2GQ-CCO2S1)*AMIN1(VOLWM(M,0,NY,NX),DFGSCO)
    RCHDFR(NY,NX)=(CCH4GQ-CCH4S1)*AMIN1(VOLWM(M,0,NY,NX),DFGSCH)
    ROXDFR(NY,NX)=(COXYGQ-COXYS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSOX)
    RNGDFR(NY,NX)=(CZ2GGQ-CZ2GS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSNG)
    RN2DFR(NY,NX)=(CZ2OGQ-CZ2OS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSN2)
    RN3DFR(NY,NX)=(CZN3GQ-CNH3S1)*AMIN1(VOLWM(M,0,NY,NX),DFGSN3)
    RHGDFR(NY,NX)=(CH2GGQ-CH2GS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSHL)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
    XCODFR(NY,NX)=XCODFR(NY,NX)+RCODFR(NY,NX)
    XCHDFR(NY,NX)=XCHDFR(NY,NX)+RCHDFR(NY,NX)
    XOXDFR(NY,NX)=XOXDFR(NY,NX)+ROXDFR(NY,NX)
    XNGDFR(NY,NX)=XNGDFR(NY,NX)+RNGDFR(NY,NX)
    XN2DFR(NY,NX)=XN2DFR(NY,NX)+RN2DFR(NY,NX)
    XN3DFR(NY,NX)=XN3DFR(NY,NX)+RN3DFR(NY,NX)
    XHGDFR(NY,NX)=XHGDFR(NY,NX)+RHGDFR(NY,NX)

  ELSE
    RCODFR(NY,NX)=0.0_r8
    RCHDFR(NY,NX)=0.0_r8
    ROXDFR(NY,NX)=0.0_r8
    RNGDFR(NY,NX)=0.0_r8
    RN2DFR(NY,NX)=0.0_r8
    RN3DFR(NY,NX)=0.0_r8
    RHGDFR(NY,NX)=0.0_r8
  ENDIF
  RCODXR=RCODFR(NY,NX)*XNPT
  RCHDXR=RCHDFR(NY,NX)*XNPT
  ROXDXR=ROXDFR(NY,NX)*XNPT
  RNGDXR=RNGDFR(NY,NX)*XNPT
  RN2DXR=RN2DFR(NY,NX)*XNPT
  RN3DXR=RN3DFR(NY,NX)*XNPT
  RHGDXR=RHGDFR(NY,NX)*XNPT
  end subroutine LitterAtmosExchange
!------------------------------------------------------------------------------------------
  subroutine SoilAtmosExchange(M,NY,NX)

  implicit none
  integer, intent(in) :: M,NY,NX

  real(r8) :: DLYR1,TORT1,VOLWOA,VOLWOB,VOLWPA,VOLWPB
  real(r8) :: DFGS(idg_beg:idg_end)
  real(r8) :: trc_gsolc(idg_beg:idg_end)
  real(r8) :: trc_gsolc2
  integer  :: K,NTG
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE SOIL LAYER
!
!     VOLWM=micropore water-filled porosity from watsub.f
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
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLWOA=VOLWM(M,NU(NY,NX),NY,NX)*VLNO3(NU(NY,NX),NY,NX)
    VOLWOB=VOLWM(M,NU(NY,NX),NY,NX)*VLNOB(NU(NY,NX),NY,NX)
    VOLWPA=VOLWM(M,NU(NY,NX),NY,NX)*VLPO4(NU(NY,NX),NY,NX)
    VOLWPB=VOLWM(M,NU(NY,NX),NY,NX)*VLPOB(NU(NY,NX),NY,NX)
    DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
    TORT1=TORT(M,NU(NY,NX),NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR1)

    DO  K=1,jcplx
      COQC2(K)=AZMAX1(OQC2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      COQN2(K)=AZMAX1(OQN2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      COQP2(K)=AZMAX1(OQP2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      COQA2(K)=AZMAX1(OQA2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    ENDDO

    CCO2S2=AZMAX1(trc_solml2(idg_CO2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CCH4S2=AZMAX1(trc_solml2(idg_CH4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    COXYS2=AZMAX1(trc_solml2(idg_O2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CZ2GS2=AZMAX1(trc_solml2(idg_N2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CZ2OS2=AZMAX1(trc_solml2(idg_N2O,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CH2GS2=AZMAX1(trc_solml2(idg_H2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    IF(VOLWMA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      CNH3S2=AZMAX1(trc_solml2(idg_NH3,NU(NY,NX),NY,NX)/VOLWMA(NU(NY,NX),NY,NX))
      CNH4S2=AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX)/VOLWMA(NU(NY,NX),NY,NX))
    ELSE
      CNH3S2=0.0_r8
      CNH4S2=0.0_r8
    ENDIF
    IF(VOLWOA.GT.ZEROS2(NY,NX))THEN
      CNO3S2=AZMAX1(trc_solml2(ids_NO3,NU(NY,NX),NY,NX)/VOLWOA)
      CNO2S2=AZMAX1(trc_solml2(ids_NO2,NU(NY,NX),NY,NX)/VOLWOA)
    ELSE
      CNO3S2=0.0_r8
      CNO2S2=0.0_r8
    ENDIF
    IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
      CP14S2=AZMAX1(trc_solml2(ids_H1PO4,NU(NY,NX),NY,NX)/VOLWPA)
      CPO4S2=AZMAX1(trc_solml2(ids_H2PO4,NU(NY,NX),NY,NX)/VOLWPA)
    ELSE
      CP14S2=0.0_r8
      CPO4S2=0.0_r8
    ENDIF
    IF(VOLWMB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      CNH3B2=AZMAX1(trc_solml2(idg_NH3B,NU(NY,NX),NY,NX)/VOLWMB(NU(NY,NX),NY,NX))
      CNH4B2=AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX)/VOLWMB(NU(NY,NX),NY,NX))
    ELSE
      CNH3B2=CNH3S2
      CNH4B2=CNH4S2
    ENDIF
    IF(VOLWOB.GT.ZEROS2(NY,NX))THEN
      CNO3B2=AZMAX1(trc_solml2(ids_NO3B,NU(NY,NX),NY,NX)/VOLWOB)
      CNO2B2=AZMAX1(trc_solml2(ids_NO2B,NU(NY,NX),NY,NX)/VOLWOB)
    ELSE
      CNO3B2=CNO3S2
      CNO2B2=CNO2S2
    ENDIF
    IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
      CP14B2=AZMAX1(trc_solml2(ids_H1PO4B,NU(NY,NX),NY,NX)/VOLWPB)
      CPO4B2=AZMAX1(trc_solml2(ids_H2PO4B,NU(NY,NX),NY,NX)/VOLWPB)
    ELSE
      CP14B2=CP14S2
      CPO4B2=CPO4S2
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

    DO NTG=idg_beg,idg_end
      DFGS(NTG)=SolDifcc(NTG,NU(NY,NX),NY,NX)*TORT1
      trc_gsolc2=AZMAX1(trc_solml2(NTG,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      trc_gsolc(NTG)=(PARG(M,NY,NX)*AtmGgms(NTG,NY,NX)*GSolbility(NTG,NU(NY,NX),NY,NX) &
        +DFGS(NTG)*trc_gsolc2)/(DFGS(NTG)+PARG(M,NY,NX))
    ENDDO
!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND SOIL SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!     R*DFS, dissolution/volatilization
!     R*DFS,X*DFS=current,cumulative gas exchange between atmosphere and soil surface water
!     gas code:*CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at soil surface
!     VOLWM=litter water volume
!     DFGS*=effective solute diffusivity
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     R*DXS=R*DFS for gas flux calculations
!
    RGasSSVol(idg_CO2,NY,NX)=(trc_gsolc(idg_CO2)-CCO2S2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGS(idg_CO2))
    RGasSSVol(idg_CH4,NY,NX)=(trc_gsolc(idg_CH4)-CCH4S2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGS(idg_CH4))
    RGasSSVol(idg_O2,NY,NX)=(trc_gsolc(idg_O2)-COXYS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGS(idg_O2))
    RGasSSVol(idg_N2,NY,NX)=(trc_gsolc(idg_N2)-CZ2GS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGS(idg_N2))
    RGasSSVol(idg_N2O,NY,NX)=(trc_gsolc(idg_N2O)-CZ2OS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGS(idg_N2O))
    RGasSSVol(idg_NH3,NY,NX)=(trc_gsolc(idg_NH3)-CNH3S2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX)*VLNH4(NU(NY,NX),NY,NX),DFGS(idg_NH3))
    RGasSSVol(idg_NH3B,NY,NX)=(trc_gsolc(idg_NH3B)-CNH3B2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX)*VLNHB(NU(NY,NX),NY,NX),DFGS(idg_NH3B))
    RGasSSVol(idg_H2,NY,NX)=(trc_gsolc(idg_H2)-CH2GS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGS(idg_H2))
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
    XCODFS(NY,NX)=XCODFS(NY,NX)+RGasSSVol(idg_CO2,NY,NX)
    XCHDFS(NY,NX)=XCHDFS(NY,NX)+RGasSSVol(idg_CH4,NY,NX)
    XOXDFS(NY,NX)=XOXDFS(NY,NX)+RGasSSVol(idg_O2,NY,NX)
    XNGDFS(NY,NX)=XNGDFS(NY,NX)+RGasSSVol(idg_N2,NY,NX)
    XN2DFS(NY,NX)=XN2DFS(NY,NX)+RGasSSVol(idg_N2O,NY,NX)
    XN3DFS(NY,NX)=XN3DFS(NY,NX)+RGasSSVol(idg_NH3,NY,NX)
    XNBDFS(NY,NX)=XNBDFS(NY,NX)+RGasSSVol(idg_NH3B,NY,NX)
    XHGDFS(NY,NX)=XHGDFS(NY,NX)+RGasSSVol(idg_H2,NY,NX)

  ELSE
    RGasSSVol(idg_CO2,NY,NX)=0.0_r8
    RGasSSVol(idg_CH4,NY,NX)=0.0_r8
    RGasSSVol(idg_O2,NY,NX)=0.0_r8
    RGasSSVol(idg_N2,NY,NX)=0.0_r8
    RGasSSVol(idg_N2O,NY,NX)=0.0_r8
    RGasSSVol(idg_NH3,NY,NX)=0.0_r8
    RGasSSVol(idg_NH3B,NY,NX)=0.0_r8
    RGasSSVol(idg_H2,NY,NX)=0.0_r8
  ENDIF
  RCODXS=RGasSSVol(idg_CO2,NY,NX)*XNPT
  RCHDXS=RGasSSVol(idg_CH4,NY,NX)*XNPT
  ROXDXS=RGasSSVol(idg_O2,NY,NX)*XNPT
  RNGDXS=RGasSSVol(idg_N2,NY,NX)*XNPT
  RN2DXS=RGasSSVol(idg_N2O,NY,NX)*XNPT
  RN3DXS=RGasSSVol(idg_NH3,NY,NX)*XNPT
  RNBDXS=RGasSSVol(idg_NH3B,NY,NX)*XNPT
  RHGDXS=RGasSSVol(idg_H2,NY,NX)*XNPT
  end subroutine SoilAtmosExchange
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSnowpackDisch(M,NY,NX)

  implicit none

  integer, intent(in) :: M, NY, NX
  integer :: ICHKL,L,L2
  real(r8) :: VFLWS,VFLWW,VFLWR
  real(r8) :: VFLWNOB,VFLWNO3,VFLWNHB
  real(r8) :: VFLWNH4,VFLWPO4,VFLWPOB
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
        RCOBLS(L2,NY,NX)=CO2W2(L,NY,NX)*VFLWW
        RCHBLS(L2,NY,NX)=CH4W2(L,NY,NX)*VFLWW
        ROXBLS(L2,NY,NX)=OXYW2(L,NY,NX)*VFLWW
        RNGBLS(L2,NY,NX)=ZNGW2(L,NY,NX)*VFLWW
        RN2BLS(L2,NY,NX)=ZN2W2(L,NY,NX)*VFLWW
        RN4BLW(L2,NY,NX)=ZN4W2(L,NY,NX)*VFLWW
        RN3BLW(L2,NY,NX)=ZN3W2(L,NY,NX)*VFLWW
        RNOBLW(L2,NY,NX)=ZNOW2(L,NY,NX)*VFLWW
        RH1PBS(L2,NY,NX)=Z1PW2(L,NY,NX)*VFLWW
        RH2PBS(L2,NY,NX)=ZHPW2(L,NY,NX)*VFLWW
        XCOBLS(L2,NY,NX)=XCOBLS(L2,NY,NX)+RCOBLS(L2,NY,NX)
        XCHBLS(L2,NY,NX)=XCHBLS(L2,NY,NX)+RCHBLS(L2,NY,NX)
        XOXBLS(L2,NY,NX)=XOXBLS(L2,NY,NX)+ROXBLS(L2,NY,NX)
        XNGBLS(L2,NY,NX)=XNGBLS(L2,NY,NX)+RNGBLS(L2,NY,NX)
        XN2BLS(L2,NY,NX)=XN2BLS(L2,NY,NX)+RN2BLS(L2,NY,NX)
        XN4BLW(L2,NY,NX)=XN4BLW(L2,NY,NX)+RN4BLW(L2,NY,NX)
        XN3BLW(L2,NY,NX)=XN3BLW(L2,NY,NX)+RN3BLW(L2,NY,NX)
        XNOBLW(L2,NY,NX)=XNOBLW(L2,NY,NX)+RNOBLW(L2,NY,NX)
        XH1PBS(L2,NY,NX)=XH1PBS(L2,NY,NX)+RH1PBS(L2,NY,NX)
        XH2PBS(L2,NY,NX)=XH2PBS(L2,NY,NX)+RH2PBS(L2,NY,NX)
      ELSE
        IF(L.LT.JS)THEN
          RCOBLS(L2,NY,NX)=0.0_r8
          RCHBLS(L2,NY,NX)=0.0_r8
          ROXBLS(L2,NY,NX)=0.0_r8
          RNGBLS(L2,NY,NX)=0.0_r8
          RN2BLS(L2,NY,NX)=0.0_r8
          RN4BLW(L2,NY,NX)=0.0_r8
          RN3BLW(L2,NY,NX)=0.0_r8
          RNOBLW(L2,NY,NX)=0.0_r8
          RH1PBS(L2,NY,NX)=0.0_r8
          RH2PBS(L2,NY,NX)=0.0_r8
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
            VFLWR=AZMAX1(AMIN1(1.0,FLQRM(M,NY,NX)/VOLWSL(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0,(FLQSM(M,NY,NX)+FLQHM(M,NY,NX))/VOLWSL(L,NY,NX)))
          ELSE
            VFLWR=CVRD(NY,NX)
            VFLWS=BARE(NY,NX)
          ENDIF
          VFLWNH4=VFLWS*VLNH4(NU(NY,NX),NY,NX)
          VFLWNHB=VFLWS*VLNHB(NU(NY,NX),NY,NX)
          VFLWNO3=VFLWS*VLNO3(NU(NY,NX),NY,NX)
          VFLWNOB=VFLWS*VLNOB(NU(NY,NX),NY,NX)
          VFLWPO4=VFLWS*VLPO4(NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*VLPOB(NU(NY,NX),NY,NX)
          RCOFLS0=CO2W2(L,NY,NX)*VFLWR
          RCHFLS0=CH4W2(L,NY,NX)*VFLWR
          ROXFLS0=OXYW2(L,NY,NX)*VFLWR
          RNGFLS0=ZNGW2(L,NY,NX)*VFLWR
          RN2FLS0=ZN2W2(L,NY,NX)*VFLWR
          RN4FLW0=ZN4W2(L,NY,NX)*VFLWR
          RN3FLW0=ZN3W2(L,NY,NX)*VFLWR
          RNOFLW0=ZNOW2(L,NY,NX)*VFLWR
          RH1PFS0=Z1PW2(L,NY,NX)*VFLWR
          RH2PFS0=ZHPW2(L,NY,NX)*VFLWR
          RCOFLS1=CO2W2(L,NY,NX)*VFLWS
          RCHFLS1=CH4W2(L,NY,NX)*VFLWS
          ROXFLS1=OXYW2(L,NY,NX)*VFLWS
          RNGFLS1=ZNGW2(L,NY,NX)*VFLWS
          RN2FLS1=ZN2W2(L,NY,NX)*VFLWS
          RN4FLW1=ZN4W2(L,NY,NX)*VFLWNH4
          RN3FLW1=ZN3W2(L,NY,NX)*VFLWNH4
          RNOFLW1=ZNOW2(L,NY,NX)*VFLWNO3
          RH1PFS1=Z1PW2(L,NY,NX)*VFLWPO4
          RH2PFS1=ZHPW2(L,NY,NX)*VFLWPO4
          RN4FLB1=ZN4W2(L,NY,NX)*VFLWNHB
          RN3FLB1=ZN3W2(L,NY,NX)*VFLWNHB
          RNOFLB1=ZNOW2(L,NY,NX)*VFLWNOB
          RH1BFB1=Z1PW2(L,NY,NX)*VFLWPOB
          RH2BFB1=ZHPW2(L,NY,NX)*VFLWPOB
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      RCOFLS0=0.0_r8
      RCHFLS0=0.0_r8
      ROXFLS0=0.0_r8
      RNGFLS0=0.0_r8
      RN2FLS0=0.0_r8
      RN4FLW0=0.0_r8
      RN3FLW0=0.0_r8
      RNOFLW0=0.0_r8
      RH1PFS0=0.0_r8
      RH2PFS0=0.0_r8
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
  real(r8) :: RFLCOG,RFLCHG,RFLOXG,RFLNGG
  real(r8) :: RFLN2G,RFLN3G,RFLH2G

!     FLQM=total water flux into soil micropore+macropore from watsub.f
!     VOLPM=air-filled porosity
!     RFL*G=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     *G2=gaseous content
!     C*E=atmospheric gas concentration from hour1.f
!

  IF(FLQM(3,NU(NY,NX),NY,NX).GT.0.0_r8)THEN
    IF(VOLPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQM(3,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    RFLCOG=VFLW*AZMAX1(trc_gasml2(idg_CO2,NU(NY,NX),NY,NX))
    RFLCHG=VFLW*AZMAX1(trc_gasml2(idg_CH4,NU(NY,NX),NY,NX))
    RFLOXG=VFLW*AZMAX1(trc_gasml2(idg_O2,NU(NY,NX),NY,NX))
    RFLNGG=VFLW*AZMAX1(trc_gasml2(idg_N2,NU(NY,NX),NY,NX))
    RFLN2G=VFLW*AZMAX1(trc_gasml2(idg_N2O,NU(NY,NX),NY,NX))
    RFLN3G=VFLW*AZMAX1(trc_gasml2(idg_NH3,NU(NY,NX),NY,NX))
    RFLH2G=VFLW*AZMAX1(trc_gasml2(idg_H2,NU(NY,NX),NY,NX))
  ELSE
    RFLCOG=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_CO2,NY,NX)
    RFLCHG=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_CH4,NY,NX)
    RFLOXG=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_O2,NY,NX)
    RFLNGG=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_N2,NY,NX)
    RFLN2G=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_N2O,NY,NX)
    RFLN3G=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_NH3,NY,NX)
    RFLH2G=-FLQM(3,NU(NY,NX),NY,NX)*AtmGgms(idg_H2,NY,NX)
  ENDIF
!     TOTAL SOIL GAS FLUX + CONVECTIVE FLUX
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux

  RCOFLG(3,NU(NY,NX),NY,NX)=RCOFLG(3,NU(NY,NX),NY,NX)+RFLCOG
  RCHFLG(3,NU(NY,NX),NY,NX)=RCHFLG(3,NU(NY,NX),NY,NX)+RFLCHG
  ROXFLG(3,NU(NY,NX),NY,NX)=ROXFLG(3,NU(NY,NX),NY,NX)+RFLOXG
  RNGFLG(3,NU(NY,NX),NY,NX)=RNGFLG(3,NU(NY,NX),NY,NX)+RFLNGG
  RN2FLG(3,NU(NY,NX),NY,NX)=RN2FLG(3,NU(NY,NX),NY,NX)+RFLN2G
  RN3FLG(3,NU(NY,NX),NY,NX)=RN3FLG(3,NU(NY,NX),NY,NX)+RFLN3G
  RHGFLG(3,NU(NY,NX),NY,NX)=RHGFLG(3,NU(NY,NX),NY,NX)+RFLH2G
  end subroutine SurfSoillAdvFlux

! ----------------------------------------------------------------------
  subroutine SurfSoilDifFlux(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: DFVCOG,DFVCHG,DFVOXG
  real(r8) :: DFVNGG,DFVN2G,DFVN3G
  real(r8) :: DFVHGG,DFLG2
  real(r8) :: CCO2G2,COXYG2,CZ2GG2,CZ2OG2,CNH3G2,CCH4G2,CH2GG2
  real(r8) :: DCO2GQ,DCH4GQ,DOXYGQ,DZ2GGQ,DZ2OGQ,DNH3GQ,DH2GGQ
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

    DO NTG=idg_beg,idg_end
      DifuscG(NTG,3,NU(NY,NX),NY,NX)=DFLG2*GasDifcc(NTG,NU(NY,NX),NY,NX)
    ENDDO
!
!     SURFACE GAS CONCENTRATIONS
!
!     C*G2=gaseous concentration
!     *G2=gaseous content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     VOLPM=air-filled porosity
!
    CCO2G2=AZMAX1(trc_gasml2(idg_CO2,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
    CCH4G2=AZMAX1(trc_gasml2(idg_CH4,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
    COXYG2=AZMAX1(trc_gasml2(idg_O2,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
    CZ2GG2=AZMAX1(trc_gasml2(idg_N2,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
    CZ2OG2=AZMAX1(trc_gasml2(idg_N2O,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
    CNH3G2=AZMAX1(trc_gasml2(idg_NH3,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
    CH2GG2=AZMAX1(trc_gasml2(idg_H2,NU(NY,NX),NY,NX)/VOLPM(M,NU(NY,NX),NY,NX))
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
    DCO2GQ=DifuscG(idg_CO2,3,NU(NY,NX),NY,NX)*PARGCO(NY,NX)/(DifuscG(idg_CO2,3,NU(NY,NX),NY,NX)+PARGCO(NY,NX))
    DCH4GQ=DifuscG(idg_CH4,3,NU(NY,NX),NY,NX)*PARGCH(NY,NX)/(DifuscG(idg_CH4,3,NU(NY,NX),NY,NX)+PARGCH(NY,NX))
    DOXYGQ=DifuscG(idg_O2,3,NU(NY,NX),NY,NX)*PARGOX(NY,NX)/(DifuscG(idg_O2,3,NU(NY,NX),NY,NX)+PARGOX(NY,NX))
    DZ2GGQ=DifuscG(idg_N2,3,NU(NY,NX),NY,NX)*PARGNG(NY,NX)/(DifuscG(idg_N2,3,NU(NY,NX),NY,NX)+PARGNG(NY,NX))
    DZ2OGQ=DifuscG(idg_N2O,3,NU(NY,NX),NY,NX)*PARGN2(NY,NX)/(DifuscG(idg_N2O,3,NU(NY,NX),NY,NX)+PARGN2(NY,NX))
    DNH3GQ=DifuscG(idg_NH3,3,NU(NY,NX),NY,NX)*PARGN3(NY,NX)/(DifuscG(idg_NH3,3,NU(NY,NX),NY,NX)+PARGN3(NY,NX))
    DH2GGQ=DifuscG(idg_H2,3,NU(NY,NX),NY,NX)*PARGH2(NY,NX)/(DifuscG(idg_H2,3,NU(NY,NX),NY,NX)+PARGH2(NY,NX))

    DFVCOG=DCO2GQ*(AtmGgms(idg_CO2,NY,NX)-CCO2G2)
    DFVCHG=DCH4GQ*(AtmGgms(idg_CH4,NY,NX)-CCH4G2)
    DFVOXG=DOXYGQ*(AtmGgms(idg_O2,NY,NX)-COXYG2)
    DFVNGG=DZ2GGQ*(AtmGgms(idg_N2,NY,NX)-CZ2GG2)
    DFVN2G=DZ2OGQ*(AtmGgms(idg_N2O,NY,NX)-CZ2OG2)
    DFVN3G=DNH3GQ*(AtmGgms(idg_NH3,NY,NX)-CNH3G2)
    DFVHGG=DH2GGQ*(AtmGgms(idg_H2,NY,NX)-CH2GG2)

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux

    RCOFLG(3,NU(NY,NX),NY,NX)=DFVCOG
    RCHFLG(3,NU(NY,NX),NY,NX)=DFVCHG
    ROXFLG(3,NU(NY,NX),NY,NX)=DFVOXG
    RNGFLG(3,NU(NY,NX),NY,NX)=DFVNGG
    RN2FLG(3,NU(NY,NX),NY,NX)=DFVN2G
    RN3FLG(3,NU(NY,NX),NY,NX)=DFVN3G
    RHGFLG(3,NU(NY,NX),NY,NX)=DFVHGG
  end subroutine SurfSoilDifFlux
! ----------------------------------------------------------------------

  subroutine SurfSoilFluxDisolVapor(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: CNH4S0,CNH3S0
  real(r8) :: CNH3B0,CNH4B0
  real(r8) :: VOLHGT(JY,JX),VOLNBT(JY,JX)
  real(r8) :: VOLN3T(JY,JX),VOLN2T(JY,JX)
  real(r8) :: VOLCOT(JY,JX),VOLCHT(JY,JX)
  real(r8) :: VOLOXT(JY,JX),VOLNGT(JY,JX)

!     VOLWM=micropore water-filled porosity from watsub.f
!     VOLW*=equivalent aqueous volume for gas
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     S*L=solubility of gas in water from hour1.f
!     VOLPM=air-filled porosity
!     R*DFG=water-air gas flux
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     R*DXS=gas exchange between atmosphere and soil surface water
!

  VOLWCO(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*GSolbility(idg_CO2,NU(NY,NX),NY,NX)
  VOLWCH(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*GSolbility(idg_CH4,NU(NY,NX),NY,NX)
  VOLWOX(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*GSolbility(idg_O2,NU(NY,NX),NY,NX)
  VOLWNG(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*GSolbility(idg_N2,NU(NY,NX),NY,NX)
  VOLWN2(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*GSolbility(idg_N2O,NU(NY,NX),NY,NX)
  VOLWN3(NU(NY,NX),NY,NX)=VOLWMA(NU(NY,NX),NY,NX)*GSolbility(idg_NH3,NU(NY,NX),NY,NX)
  VOLWNB(NU(NY,NX),NY,NX)=VOLWMB(NU(NY,NX),NY,NX)*GSolbility(idg_NH3,NU(NY,NX),NY,NX)
  VOLWHG(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*GSolbility(idg_H2,NU(NY,NX),NY,NX)

  VOLCOT(NY,NX)=VOLWCO(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
  VOLCHT(NY,NX)=VOLWCH(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
  VOLOXT(NY,NX)=VOLWOX(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
  VOLNGT(NY,NX)=VOLWNG(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
  VOLN2T(NY,NX)=VOLWN2(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
  VOLN3T(NY,NX)=VOLWN3(NU(NY,NX),NY,NX)+VOLPMA(NU(NY,NX),NY,NX)
  VOLNBT(NY,NX)=VOLWNB(NU(NY,NX),NY,NX)+VOLPMB(NU(NY,NX),NY,NX)
  VOLHGT(NY,NX)=VOLWHG(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
  RGasDSFlx(idg_CO2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_CO2,NU(NY,NX),NY,NX)) &
    *VOLWCO(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_CO2,NU(NY,NX),NY,NX)+RCODXS) &
    *VOLPM(M,NU(NY,NX),NY,NX))/VOLCOT(NY,NX)
  RGasDSFlx(idg_CH4,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_CH4,NU(NY,NX),NY,NX)) &
    *VOLWCH(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_CH4,NU(NY,NX),NY,NX)+RCHDXS) &
    *VOLPM(M,NU(NY,NX),NY,NX))/VOLCHT(NY,NX)
  RGasDSFlx(idg_O2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_O2,NU(NY,NX),NY,NX)) &
    *VOLWOX(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_O2,NU(NY,NX),NY,NX)+ROXDXS) &
    *VOLPM(M,NU(NY,NX),NY,NX))/VOLOXT(NY,NX)
  RGasDSFlx(idg_N2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_N2,NU(NY,NX),NY,NX)) &
    *VOLWNG(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_N2,NU(NY,NX),NY,NX)+RNGDXS) &
  *VOLPM(M,NU(NY,NX),NY,NX))/VOLNGT(NY,NX)
  RGasDSFlx(idg_N2O,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_N2O,NU(NY,NX),NY,NX)) &
    *VOLWN2(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_N2O,NU(NY,NX),NY,NX)+RN2DXS) &
    *VOLPM(M,NU(NY,NX),NY,NX))/VOLN2T(NY,NX)
  IF(VOLN3T(NY,NX).GT.ZEROS2(NY,NX).AND.VOLWXA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_NH3,NU(NY,NX),NY,NX)) &
      *VOLWN3(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,trc_solml2(idg_NH3,NU(NY,NX),NY,NX)+RN3DXS) &
      *VOLPMA(NU(NY,NX),NY,NX))/VOLN3T(NY,NX)
    CNH3S0=AZMAX1((trc_solml2(idg_NH3,NU(NY,NX),NY,NX) &
      +RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX))/VOLWXA(NU(NY,NX),NY,NX))
    CNH4S0=AZMAX1(trc_solml2(ids_NH4,NU(NY,NX),NY,NX))/VOLWXA(NU(NY,NX),NY,NX)
  ELSE
    RGasDSFlx(idg_NH3,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  IF(VOLNBT(NY,NX).GT.ZEROS2(NY,NX).AND.VOLWXB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    RNBDFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_NH3,NU(NY,NX),NY,NX)) &
      *VOLWNB(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,trc_solml2(idg_NH3B,NU(NY,NX),NY,NX)+RNBDXS) &
      *VOLPMB(NU(NY,NX),NY,NX))/VOLNBT(NY,NX)
    CNH3B0=AZMAX1((trc_solml2(idg_NH3B,NU(NY,NX),NY,NX) &
      +RNBDFG(NU(NY,NX),NY,NX))/VOLWXB(NU(NY,NX),NY,NX))
    CNH4B0=AZMAX1(trc_solml2(ids_NH4B,NU(NY,NX),NY,NX))/VOLWXB(NU(NY,NX),NY,NX)
  ELSE
    RNBDFG(NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  RGasDSFlx(idg_H2,NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
    *(AMAX1(ZEROS(NY,NX),trc_gasml2(idg_H2,NU(NY,NX),NY,NX)) &
    *VOLWHG(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
    ,trc_solml2(idg_H2,NU(NY,NX),NY,NX)+RHGDXS) &
    *VOLPM(M,NU(NY,NX),NY,NX))/VOLHGT(NY,NX)
  end subroutine SurfSoilFluxDisolVapor
! ----------------------------------------------------------------------
end module SurfaceFluxMod
