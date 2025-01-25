module SurfaceFluxMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: AZMAX1, AZMIN1, isclose
  use EcoSIMCtrlMod, only: PrintInfo
  use GridDataType  
  use TranspNoSaltDataMod
  use SoilBGCDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SurfLitterDataType
  use AqueChemDatatype
  use SnowDataType
  use EcoSimConst
  use TracerIDMod
  use ChemTranspDataType
  use ClimForcDataType
  use EcoSIMSolverPar
  use LandSurfDataType
  use SurfSoilDataType
  use SoilPropertyDataType
implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: SoluteFluxSurfaceM
  public :: LitterGasVolatilDissolMM
  public :: SurfSoilFluxGasDifAdvMM
  public :: SnowSoluteDischargeM
contains


!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSurfaceM(I,J,M,NY,NX,NHE,NHW,NVS,NVN,&
    WaterFlow2Soil,trcg_FloSno2LitR,trcn_FloSno2LitR,RGas_Dif_Atm2Soil_FlxMM,&
    RGas_Dif_Atm2Litr_FlxMM,RGasAtmDisol2SoilM,RGasAtmDisol2LitrM)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX,M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8),intent(in) :: trcg_FloSno2LitR(idg_beg:idg_NH3)
  real(r8),intent(in) :: trcn_FloSno2LitR(ids_nut_beg:ids_nuts_end)
  real(r8),intent(inout) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8),intent(out) :: RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)
  real(r8),intent(out) :: RGas_Dif_Atm2Litr_FlxMM(idg_beg:idg_end-1)
  real(r8),intent(out) :: RGasAtmDisol2SoilM(idg_beg:idg_end)
  real(r8),intent(out) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3)
  character(len=*), parameter :: subname='SoluteFluxSurfaceM'

  real(r8) :: FLWRM1
  real(r8) :: trcs_cl1(ids_beg:ids_end)
  real(r8) :: trcs_cl2(ids_beg:ids_end)
  real(r8) :: RFLs_adv(ids_beg:ids_end)
  real(r8) :: SDifFlx(ids_beg:ids_end)


!     VLWatMicPM,VLWatMacPM,VLsoiAirPM,ReductVLsoiAirPM=micropore,macropore water volume, air volume and change in air volume
!     WaterFlow2MicPM_3D,WaterFlow2MacPM_3D=water flux into soil micropore,macropore from watsub.f
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     CumReductVLsoiAirPM,FLM=air,water flux in gas flux calculations
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!
  call PrintInfo('beg '//subname)

  VLWatMicPMA_vr(NU(NY,NX),NY,NX)      = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
  VLWatMicPMB_vr(NU(NY,NX),NY,NX)      = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
  VLWatMicPXA(NU(NY,NX),NY,NX)         = natomw*VLWatMicPMA_vr(NU(NY,NX),NY,NX)
  VLWatMicPXB(NU(NY,NX),NY,NX)         = natomw*VLWatMicPMB_vr(NU(NY,NX),NY,NX)
  VLsoiAirPMA(NU(NY,NX),NY,NX)         = VLsoiAirPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
  VLsoiAirPMB(NU(NY,NX),NY,NX)         = VLsoiAirPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
  CumReductVLsoiAirPM(NU(NY,NX),NY,NX) = ReductVLsoiAirPM(M,NU(NY,NX),NY,NX)*dt_GasCyc
  WaterFlow2Soil(3,NU(NY,NX),NY,NX)    = (WaterFlow2MicPM_3D(M,3,NU(NY,NX),NY,NX)+WaterFlow2MacPM_3D(M,3,NU(NY,NX),NY,NX))*dt_GasCyc
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
!     DiffusivitySolutEff*=effective solute diffusivity
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in litter
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in litter
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in litter
!     C*1=solute concentration in litter
!
  call LitterAtmosExchangeM(I,J,M,NY,NX,trcs_cl1,RGas_Dif_Atm2Litr_FlxMM,RGasAtmDisol2LitrM)

  call SoilAtmosExchangeM(I,J,M,NY,NX,trcs_cl2,RGas_Dif_Atm2Soil_FlxMM,RGasAtmDisol2SoilM)

!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
  call ConvectiveSurfaceSoluteFlux(I,J,M,NY,NX,FLWRM1,RFLs_adv)
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  call LitterSurfSoilExchange(I,J,M,NY,NX,FLWRM1,trcs_cl1,trcs_cl2,SDifFlx)

!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES
!
  call TotalPoreFluxAdjacentCell(I,J,M,NY,NX,SDifFlx,RFLs_adv,trcg_FloSno2LitR,trcn_FloSno2LitR)

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
  call OverLandTracerTranspt(M,NY,NX,NHE,NHW,NVS,NVN)

  call PrintInfo('end '//subname)
  end subroutine SoluteFluxSurfaceM

!------------------------------------------------------------------------------------------

  subroutine ConvectiveSurfaceSoluteFlux(I,J,M,NY,NX,FLWRM1,RFLs_adv)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, NY, NX
  real(r8),intent(out) :: RFLs_adv(ids_beg:ids_end)
  real(r8),intent(out) :: FLWRM1
  REAL(R8) :: VFLW
  integer :: K,nnut,idom,idg,nnut1

  FLWRM1=WatFLo2LitrM(M,NY,NX)
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
    IF(VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VLWatMicPM_vr(M,0,NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2(idom,K,0,NY,NX))
      ENDDO
    ENDDO

    DO idg=idg_beg,idg_end
      RFLs_adv(idg)=VFLW*AZMAX1(trc_solml2_vr(idg,0,NY,NX))
    enddo

    !from NH4 to H2PO4
    DO nnut=ids_nut_beg,ids_nuts_end
      RFLs_adv(nnut)=VFLW*AZMAX1(trc_solml2_vr(nnut,0,NY,NX))*trcs_VLN_vr(nnut,NU(NY,NX),NY,NX)
    ENDDO

    !from NH4B to ids_H2PO4B
    DO nnut=0,ids_nuts     
      nnut1=ids_NH4B+nnut
      RFLs_adv(nnut1)=VFLW*AZMAX1(trc_solml2_vr(nnut+ids_NH4,0,NY,NX))*trcs_VLN_vr(nnut1,NU(NY,NX),NY,NX)
    ENDDO

    RFLs_adv(idg_NH3B)=VFLW*AZMAX1(trc_solml2_vr(idg_NH3,0,NY,NX))*trcs_VLN_vr(idg_NH3B,NU(NY,NX),NY,NX)

!    write(117,*)I+J/24.,M,'NO3adv',RFLs_adv(ids_NO3),VFLW,AZMAX1(trc_solml2_vr(ids_NO3,0,NY,NX)),trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)
!    write(117,*)I+J/24.,M,'NO3Badv',RFLs_adv(ids_NO3B),VFLW,AZMAX1(trc_solml2_vr(ids_NO3B,0,NY,NX)),trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)

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
    IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2(idom,K,NU(NY,NX),NY,NX))
      ENDDO
    ENDDO

    DO nnut=ids_beg,ids_end
      RFLs_adv(nnut)=VFLW*AZMAX1(trc_solml2_vr(nnut,NU(NY,NX),NY,NX))
    ENDDO
  ENDIF
  end subroutine ConvectiveSurfaceSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine LitterSurfSoilExchange(I,J,M,NY,NX,FLWRM1,trcs_cl1,trcs_cl2,SDifFlx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, NY, NX
  real(r8), intent(in) :: FLWRM1
  real(r8), intent(in) :: trcs_cl1(ids_beg:ids_end)
  real(r8), intent(in) :: trcs_cl2(ids_beg:ids_end)
  real(r8), intent(out):: SDifFlx(ids_beg:ids_end)

  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: DIFDOM(idom_beg:idom_end)

  real(r8) :: DIFDOM0
  real(r8) :: DIFDOM1
  real(r8) :: DISPN,DIF0,DIF1
  real(r8) :: SDifc(ids_beg:ids_end)
  integer  :: K,nnut,idg,idom
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!
  IF((VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND.(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN
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
    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX))
    TORT0 = TortMicPM_vr(M,0,NY,NX)/DLYR0*FracSurfByLitR_col(NY,NX)
    DLYR1 = AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))
    TORT1 = TortMicPM_vr(M,NU(NY,NX),NY,NX)/DLYR1
    DISPN = DISP(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))


    DO idom=idom_beg,idom_end
      DIFDOM0      = (DOMdiffusivity2_vr(idom,0,NY,NX)*TORT0+DISPN)
      DIFDOM1      = (DOMdiffusivity2_vr(idom,NU(NY,NX),NY,NX)*TORT1+DISPN)
      DIFDOM(idom) = DIFDOM0*DIFDOM1/(DIFDOM0+DIFDOM1)*AREA(3,NU(NY,NX),NY,NX)
    ENDDO
    DO nnut=ids_beg,ids_end
      DIF0        = (SoluteDifusvtytscal_vr(nnut,0,NY,NX)*TORT0+DISPN)
      DIF1        = (SoluteDifusvtytscal_vr(nnut,NU(NY,NX),NY,NX)*TORT1+DISPN)
      SDifc(nnut) = DIF0*DIF1/(DIF0+DIF1)*AREA(3,NU(NY,NX),NY,NX)
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
      DO idom=idom_beg,idom_end
        Difus_Micp_flx_DOM(idom,K)=DIFDOM(idom)*(CDOM_MicP1(idom,K)-CDOM_MicP2(idom,K))
      ENDDO
    ENDDO
!exclude NH3B and NH3
    DO idg=idg_beg,idg_end-2
        SDifFlx(idg)=SDifc(idg)*(trcs_cl1(idg)-trcs_cl2(idg))
    ENDDO

    DO nnut=ids_nuts_beg,ids_nuts_end
      SDifFlx(nnut)=SDifc(nnut)*(trcs_cl1(nnut)-trcs_cl2(nnut))*trcs_VLN_vr(nnut,NU(NY,NX),NY,NX)
    ENDDO
  ELSE
    DO  K=1,jcplx
      Difus_Micp_flx_DOM(idom_beg:idom_end,K)=0.0_r8
    ENDDO
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF
  end subroutine LitterSurfSoilExchange


!------------------------------------------------------------------------------------------

  subroutine TotalPoreFluxAdjacentCell(I,J,M,NY,NX,SDifFlx,RFLs_adv,trcg_FloSno2LitR,trcn_FloSno2LitR)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY, NX
  real(r8), intent(in) :: SDifFlx(ids_beg:ids_end)
  real(r8), intent(in) :: RFLs_adv(ids_beg:ids_end)
  real(r8), intent(in) :: trcg_FloSno2LitR(idg_beg:idg_NH3)
  real(r8), intent(in) :: trcn_FloSno2LitR(ids_nut_beg:ids_nuts_end)
  integer :: K,ntg,idom,nts

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
  DO K=1,micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)=RDOMFL0(idom,K,NY,NX)-DOM_Adv2MicP_flx(idom,K)-Difus_Micp_flx_DOM(idom,K)
      DOM_MicpTranspFlxM_3D(idom,K,3,NU(NY,NX),NY,NX)=RDOMFL1(idom,K,NY,NX)+DOM_Adv2MicP_flx(idom,K)+Difus_Micp_flx_DOM(idom,K)
    ENDDO
  ENDDO

  DO ntg=idg_beg,idg_NH3
    R3PoreSolFlx_3D(ntg,3,0,NY,NX)=trcg_RFL0(ntg,NY,NX)+trcg_FloSno2LitR(ntg)-RFLs_adv(ntg)-SDifFlx(ntg)
  ENDDO

  do nts=ids_nut_beg,ids_nuts_end
    R3PoreSolFlx_3D(nts,3,0,NY,NX)=trcn_RFL0(nts,NY,NX)+trcn_FloSno2LitR(nts)-RFLs_adv(nts)-SDifFlx(nts)
  enddo
  R3PoreSolFlx_3D(idg_NH3,3,0,NY,NX)=R3PoreSolFlx_3D(idg_NH3,3,0,NY,NX)-RFLs_adv(idg_NH3B)-SDifFlx(idg_NH3B)
  R3PoreSolFlx_3D(ids_NH4,3,0,NY,NX)=R3PoreSolFlx_3D(ids_NH4,3,0,NY,NX)-RFLs_adv(ids_NH4B)-SDifFlx(ids_NH4B)
  R3PoreSolFlx_3D(ids_NO3,3,0,NY,NX)=R3PoreSolFlx_3D(ids_NO3,3,0,NY,NX)-RFLs_adv(ids_NO3B)-SDifFlx(ids_NO3B)
  R3PoreSolFlx_3D(ids_NO2,3,0,NY,NX)=R3PoreSolFlx_3D(ids_NO2,3,0,NY,NX)-RFLs_adv(ids_NO2B)-SDifFlx(ids_NO2B)
  R3PoreSolFlx_3D(ids_H1PO4,3,0,NY,NX)=R3PoreSolFlx_3D(ids_H1PO4,3,0,NY,NX)-RFLs_adv(ids_H1PO4B)-SDifFlx(ids_H1PO4B)
  R3PoreSolFlx_3D(ids_H2PO4,3,0,NY,NX)=R3PoreSolFlx_3D(ids_H2PO4,3,0,NY,NX)-RFLs_adv(ids_H2PO4B)-SDifFlx(ids_H2PO4B)

  DO ntg=idg_beg,idg_NH3
    R3PoreSolFlx_3D(ntg,3,NU(NY,NX),NY,NX)=trcs_RFL1(ntg,NY,NX)+trcg_VFloSnow(ntg)+RFLs_adv(ntg)+SDifFlx(ntg)
  ENDDO

  do nts=ids_nut_beg,ids_nuts_end
    R3PoreSolFlx_3D(nts,3,NU(NY,NX),NY,NX)=trcs_RFL1(nts,NY,NX)+trcn_soil_VFloSnow(nts)+RFLs_adv(nts)+SDifFlx(nts)
  enddo

  do nts=ids_nutb_beg,ids_nutb_end
    R3PoreSolFlx_3D(nts,3,NU(NY,NX),NY,NX)=trcs_RFL1(nts,NY,NX)+trcn_band_VFloSnow(nts)+RFLs_adv(nts)+SDifFlx(nts)
  enddo
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  DO K=1,micpar%NumOfLitrCmplxs
    do idom=idom_beg,idom_end
      DOM_MicpTransp_3D(idom,K,3,0,NY,NX)=DOM_MicpTransp_3D(idom,K,3,0,NY,NX) &
        -DOM_Adv2MicP_flx(idom,K)-Difus_Micp_flx_DOM(idom,K)
      DOM_MicpTransp_3D(idom,K,3,NU(NY,NX),NY,NX)=DOM_MicpTransp_3D(idom,K,3,NU(NY,NX),NY,NX) &
        +DOM_Adv2MicP_flx(idom,K)+Difus_Micp_flx_DOM(idom,K)
    ENDDO
  ENDDO

  do ntg=idg_beg,idg_NH3
    trcs_TransptMicP_3D(ntg,3,0,NY,NX)=trcs_TransptMicP_3D(ntg,3,0,NY,NX)+trcg_FloSno2LitR(ntg)-RFLs_adv(ntg)-SDifFlx(ntg)
  ENDDO
!  write(116,*)I+J/24.,M,'tp30',trcs_TransptMicP_3D(ids_NO3,3,0,NY,NX)

  do nts=ids_nut_beg,ids_nuts_end
    trcs_TransptMicP_3D(nts,3,0,NY,NX)=trcs_TransptMicP_3D(nts,3,0,NY,NX)+trcn_FloSno2LitR(nts)-RFLs_adv(nts)-SDifFlx(nts)
  enddo

!  write(116,*)I+J/24.,M,'tp31',trcs_TransptMicP_3D(ids_NO3,3,0,NY,NX),trcn_FloSno2LitR(ids_NO3),-RFLs_adv(ids_NO3),-SDifFlx(ids_NO3)

  trcs_TransptMicP_3D(idg_NH3,3,0,NY,NX)=trcs_TransptMicP_3D(idg_NH3,3,0,NY,NX)-RFLs_adv(idg_NH3B)-SDifFlx(idg_NH3B)  
  trcs_TransptMicP_3D(ids_NH4,3,0,NY,NX)=trcs_TransptMicP_3D(ids_NH4,3,0,NY,NX)-RFLs_adv(ids_NH4B)-SDifFlx(ids_NH4B)
  trcs_TransptMicP_3D(ids_NO3,3,0,NY,NX)=trcs_TransptMicP_3D(ids_NO3,3,0,NY,NX)-RFLs_adv(ids_NO3B)-SDifFlx(ids_NO3B)  
  trcs_TransptMicP_3D(ids_NO2,3,0,NY,NX)=trcs_TransptMicP_3D(ids_NO2,3,0,NY,NX)-RFLs_adv(ids_NO2B)-SDifFlx(ids_NO2B)  
  trcs_TransptMicP_3D(ids_H1PO4,3,0,NY,NX)=trcs_TransptMicP_3D(ids_H1PO4,3,0,NY,NX)-RFLs_adv(ids_H1PO4B)-SDifFlx(ids_H1PO4B)
  trcs_TransptMicP_3D(ids_H2PO4,3,0,NY,NX)=trcs_TransptMicP_3D(ids_H2PO4,3,0,NY,NX)-RFLs_adv(ids_H2PO4B)-SDifFlx(ids_H2PO4B)

!  write(116,*)I+J/24.,M,'tp32',trcs_TransptMicP_3D(ids_NO3,3,0,NY,NX),-RFLs_adv(ids_NO3B),-SDifFlx(ids_NO3B)  

  do ntg=idg_beg,idg_NH3
    trcs_TransptMicP_3D(ntg,3,NU(NY,NX),NY,NX)=trcs_TransptMicP_3D(ntg,3,NU(NY,NX),NY,NX)+trcg_VFloSnow(ntg)+RFLs_adv(ntg)+SDifFlx(ntg)
  enddo

  do nts=ids_nut_beg,ids_nuts_end
    trcs_TransptMicP_3D(nts,3,NU(NY,NX),NY,NX)=trcs_TransptMicP_3D(nts,3,NU(NY,NX),NY,NX)+trcn_soil_VFloSnow(nts)+RFLs_adv(nts)+SDifFlx(nts)
  enddo

  do nts=ids_nutb_beg,ids_nutb_end
    trcs_TransptMicP_3D(nts,3,NU(NY,NX),NY,NX)=trcs_TransptMicP_3D(nts,3,NU(NY,NX),NY,NX)+trcn_band_VFloSnow(nts)+RFLs_adv(nts)+SDifFlx(nts)
  enddo
  end subroutine TotalPoreFluxAdjacentCell
!------------------------------------------------------------------------------------------

  subroutine MacMicPoresTransfer(M,NY,NX)
  implicit none

  integer, intent(in) :: M, NY, NX
  integer :: K,nnut,idg,nsol,idom
  real(r8) :: VFLW
  real(r8) :: VLWatMacPS,VOLWT
  real(r8) :: SDifFlx(ids_beg:ids_end)
  real(r8) :: SAdvFlx(ids_beg:ids_end)
!
!     FWatExMacP2MicPM_vr=macro-micropore water transfer from watsub.f
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
  IF(FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX).GT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX)/VLWatMacPM_vr(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MacP2(idom,K,NU(NY,NX),NY,NX))
      enddo
    ENDDO

    do idg=idg_beg,idg_NH3-1
      SAdvFlx(idg)=VFLW*AZMAX1(trc_soHml2_vr(idg,NU(NY,NX),NY,NX))
    enddo
    SAdvFlx(idg_NH3)=VFLW*AZMAX1(trc_soHml2_vr(idg_NH3,NU(NY,NX),NY,NX))*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3B)=VFLW*AZMAX1(trc_soHml2_vr(idg_NH3B,NU(NY,NX),NY,NX))*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)

    do nsol=ids_beg,ids_end
      SAdvFlx(nsol)=VFLW*AZMAX1(trc_soHml2_vr(nsol,NU(NY,NX),NY,NX))*trcs_VLN_vr(nsol,NU(NY,NX),NY,NX)
    enddo 
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX).LT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2(idom,K,NU(NY,NX),NY,NX))
      enddo
    ENDDO
    do idg=idg_beg,idg_NH3-1
      SAdvFlx(idg)=VFLW*AZMAX1(trc_solml2_vr(idg,NU(NY,NX),NY,NX))
    enddo
    SAdvFlx(idg_NH3)  = VFLW*AZMAX1(trc_solml2_vr(idg_NH3,NU(NY,NX),NY,NX))*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    SAdvFlx(idg_NH3B) = VFLW*AZMAX1(trc_solml2_vr(idg_NH3B,NU(NY,NX),NY,NX))*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)

    do nsol=ids_beg,ids_end
      SAdvFlx(nsol)=VFLW*AZMAX1(trc_solml2_vr(nsol,NU(NY,NX),NY,NX))*trcs_VLN_vr(nsol,NU(NY,NX),NY,NX)
    enddo
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    DO  K=1,jcplx      
      DOM_Adv2MicP_flx(idom_beg:idom_end,K)=0.0_r8
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
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *H2,*2=macropore,micropore solute content
!
  IF(VLWatMacPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMacPS = AMIN1(XFRS*VGeomLayer_vr(NU(NY,NX),NY,NX),VLWatMacPM_vr(M,NU(NY,NX),NY,NX))
    VOLWT      = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)+VLWatMacPS
    DO K = 1, jcplx
      !diffusion flux in micropore
      do idom=idom_beg,idom_end
        Difus_Micp_flx_DOM(idom,K)=dts_HeatWatTP*(AZMAX1(DOM_MacP2(idom,K,NU(NY,NX),NY,NX))*VLWatMicPM_vr(M,NU(NY,NX),NY,NX) &
          -AZMAX1(DOM_MicP2(idom,K,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
      enddo
    ENDDO
!    diffusion flux in macropores
    do idg=idg_beg,idg_NH3-1
      SDifFlx(idg)=dts_HeatWatTP*(AZMAX1(trc_soHml2_vr(idg,NU(NY,NX),NY,NX))*VLWatMicPM_vr(M,NU(NY,NX),NY,NX) &
        -AZMAX1(trc_solml2_vr(idg,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT
    enddo  
    SDifFlx(idg_NH3)=dts_HeatWatTP*(AZMAX1(trc_soHml2_vr(idg_NH3,NU(NY,NX),NY,NX))*VLWatMicPM_vr(M,NU(NY,NX),NY,NX) &
      -AZMAX1(trc_solml2_vr(idg_NH3,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    SDifFlx(idg_NH3B)=dts_HeatWatTP*(AZMAX1(trc_soHml2_vr(idg_NH3B,NU(NY,NX),NY,NX))*VLWatMicPM_vr(M,NU(NY,NX),NY,NX) &
      -AZMAX1(trc_solml2_vr(idg_NH3B,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
      *trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)

    !diffusion flux from macropre to micropore
    do nsol=ids_beg,ids_end
      SDifFlx(nsol)=dts_HeatWatTP*(AZMAX1(trc_soHml2_vr(nsol,NU(NY,NX),NY,NX))*VLWatMicPM_vr(M,NU(NY,NX),NY,NX) &
        -AZMAX1(trc_solml2_vr(nsol,NU(NY,NX),NY,NX))*VLWatMacPS)/VOLWT &
        *trcs_VLN_vr(nsol,NU(NY,NX),NY,NX)
    enddo    

  ELSE
    DO  K=1,jcplx
      Difus_Micp_flx_DOM(idom_beg:idom_end,K)=0.0_r8
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
    do idom=idom_beg,idom_end
      DOM_XPoreTransp_flx(idom,K,NU(NY,NX),NY,NX)=DOM_Adv2MicP_flx(idom,K)+Difus_Micp_flx_DOM(idom,K)
    enddo
  ENDDO

  DO nnut=ids_beg,ids_end
    RMac2MicSolFlx_vr(nnut,NU(NY,NX),NY,NX)=SAdvFlx(nnut)+SDifFlx(nnut)
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
    do idom=idom_beg,idom_end
      DOM_PoreTranspFlx(idom,K,NU(NY,NX),NY,NX)=DOM_PoreTranspFlx(idom,K,NU(NY,NX),NY,NX)+DOM_XPoreTransp_flx(idom,K,NU(NY,NX),NY,NX)
    enddo
  ENDDO

  DO nnut=ids_beg,ids_end
    trcs_Mac2MicXfer_vr(nnut,NU(NY,NX),NY,NX)=trcs_Mac2MicXfer_vr(nnut,NU(NY,NX),NY,NX)+RMac2MicSolFlx_vr(nnut,NU(NY,NX),NY,NX)
  ENDDO

!
  end subroutine MacMicPoresTransfer
!------------------------------------------------------------------------------------------

  subroutine OverLandTracerTranspt(M,NY,NX,NHE,NHW,NVS,NVN)
  implicit none

  integer, intent(in) :: M, NY, NX, NHE, NHW, NVS, NVN
  real(r8) :: FQRM,VFLW
  integer :: K,N,NN,N1,N2,N4,N5,N4B,N5B,idg,nnut
  integer :: idom
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
  N1=NX;N2=NY
  IF(WatFlux4ErosionM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VLWatMicPM_vr(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,WatFlux4ErosionM_2DH(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
      VFLW=AMIN1(VFLWX,WatFlux4ErosionM_2DH(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        dom_FloXSurRunoff(idom,K,N2,N1)=VFLW*AZMAX1(DOM_MicP2(idom,K,0,N2,N1))
      enddo
    ENDDO
    DO idg=idg_beg,idg_NH3
      trcg_FloXSurRunoff(idg,N2,N1)=VFLW*AZMAX1(trc_solml2_vr(idg,0,N2,N1))
    ENDDO

    DO nnut=ids_nut_beg,ids_nuts_end
      trcn_FloXSurRunoff(nnut,N2,N1)=VFLW*AZMAX1(trc_solml2_vr(nnut,0,N2,N1))
    ENDDO
  ELSE
    DO K=1,jcplx
      dom_FloXSurRunoff(idom_beg:idom_end,K,N2,N1)=0.0_r8
    ENDDO
    trcg_FloXSurRunoff(idg_beg:idg_NH3,N2,N1)          = 0.0_r8
    trcn_FloXSurRunoff(ids_nut_beg:ids_nuts_end,N2,N1) = 0.0_r8
  ENDIF
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  D4310: DO N=1,2
    D4305: DO NN=1,2
      IF(N.EQ.1)THEN
        !east-west
        IF((NX.EQ.NHE.AND.NN.EQ.1) .OR. (NX.EQ.NHW.AND.NN.EQ.2))THEN
          cycle
        ELSE
          N4  = NX+1
          N5  = NY
          N4B = NX-1
          N5B = NY
        ENDIF
      ELSEIF(N.EQ.2)THEN
        !south-north
        IF((NY.EQ.NVS.AND.NN.EQ.1) .OR. (NY.EQ.NVN.AND.NN.EQ.2))THEN
          cycle
        ELSE
          N4  = NX
          N5  = NY+1
          N4B = NX
          N5B = NY-1
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
      IF(WatFlux4ErosionM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.1)THEN
          FQRM=QflxSurfRunoffM_2DH(M,N,2,N5,N4)/WatFlux4ErosionM_2DH(M,N2,N1)
          DO  K=1,jcplx
            do idom=idom_beg,idom_end
              dom_2DFloXSurRunoffM(idom,K,N,2,N5,N4)=dom_FloXSurRunoff(idom,K,N2,N1)*FQRM
            enddo
          ENDDO

          DO idg=idg_beg,idg_NH3
            trcg_2DFloXSurRunoffM(idg,N,2,N5,N4)=trcg_FloXSurRunoff(idg,N2,N1)*FQRM
          ENDDO

          DO nnut=ids_nut_beg,ids_nuts_end
            trcn_2DFloXSurRunoffM(nnut,N,2,N5,N4)=trcn_FloXSurRunoff(nnut,N2,N1)*FQRM
          ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
          DO  K=1,jcplx
            do idom=idom_beg,idom_end
              dom_2DFloXSurRunoff(idom,K,N,2,N5,N4)=dom_2DFloXSurRunoff(idom,K,N,2,N5,N4)+dom_2DFloXSurRunoffM(idom,K,N,2,N5,N4)
            enddo
          ENDDO
          DO idg=idg_beg,idg_NH3
            trcg_FloXSurRunoff_2D(idg,N,2,N5,N4)=trcg_FloXSurRunoff_2D(idg,N,2,N5,N4)+trcg_2DFloXSurRunoffM(idg,N,2,N5,N4)
          ENDDO

          DO nnut=ids_nut_beg,ids_nuts_end
            trcn_FloXSurRunoff_2D(nnut,N,2,N5,N4)=trcn_FloXSurRunoff_2D(nnut,N,2,N5,N4)+trcn_2DFloXSurRunoffM(nnut,N,2,N5,N4)
          ENDDO
        ELSE
          DO K=1,jcplx
            dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,2,N5,N4)=0.0_r8
          ENDDO
          trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,N,2,N5,N4)=0.0_r8

          trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,2,N5,N4)=0.0_r8
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            FQRM=QflxSurfRunoffM_2DH(M,N,1,N5B,N4B)/WatFlux4ErosionM_2DH(M,N2,N1)
            DO  K=1,jcplx
              do idom=idom_beg,idom_end
                dom_2DFloXSurRunoffM(idom,K,N,1,N5B,N4B)=dom_FloXSurRunoff(idom,K,N2,N1)*FQRM
              enddo
            ENDDO

            DO idg=idg_beg,idg_NH3
              trcg_2DFloXSurRunoffM(idg,N,1,N5B,N4B)=trcg_FloXSurRunoff(idg,N2,N1)*FQRM
            ENDDO

            DO nnut=ids_nut_beg,ids_nuts_end
              trcn_2DFloXSurRunoffM(nnut,N,1,N5B,N4B)=trcn_FloXSurRunoff(nnut,N2,N1)*FQRM
            ENDDO

            DO K=1,jcplx
              do idom=idom_beg,idom_end
                dom_2DFloXSurRunoff(idom,K,N,1,N5B,N4B)=dom_2DFloXSurRunoff(idom,K,N,1,N5B,N4B)+dom_2DFloXSurRunoffM(idom,K,N,1,N5B,N4B)
              enddo
            ENDDO

            DO idg=idg_beg,idg_NH3
              trcg_FloXSurRunoff_2D(idg,N,1,N5B,N4B)=trcg_FloXSurRunoff_2D(idg,N,1,N5B,N4B)+trcg_2DFloXSurRunoffM(idg,N,1,N5B,N4B)
            ENDDO

            DO nnut=ids_nut_beg,ids_nuts_end
              trcn_FloXSurRunoff_2D(nnut,N,1,N5B,N4B)=trcn_FloXSurRunoff_2D(nnut,N,1,N5B,N4B)+trcn_2DFloXSurRunoffM(nnut,N,1,N5B,N4B)
            ENDDO
          ELSE
            DO  K=1,jcplx
              dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,1,N5B,N4B)=0.0_r8
            ENDDO
            trcg_2DFloXSurRunoffM(idg_beg:idg_end,N,1,N5B,N4B)          = 0.0_r8
            trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,1,N5B,N4B) = 0.0_r8
          ENDIF
        ENDIF
      ELSE
        DO K=1,jcplx
          dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,2,N5,N4)=0.0_r8
        ENDDO
        trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,N,2,N5,N4)          = 0.0_r8
        trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,2,N5,N4) = 0.0_r8

        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          DO  K=1,jcplx
            dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,1,N5B,N4B)=0.0_r8
          ENDDO
          trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,N,1,N5B,N4B)          = 0.0_r8
          trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,1,N5B,N4B) = 0.0_r8
        ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     subroutine SnowdriftTransport(M)
!
!     DrySnoFlxBySnoRedistM=snow transfer from watsub.f
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
        IF(ABS(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
          trcg_2DSnowDrift(idg_beg:idg_NH3,N,N5,N4)          = 0.0_r8
          trcn_2DSnowDrift(ids_nut_beg:ids_nuts_end,N,N5,N4) = 0.0_r8
          VFLW                                               = 0._r8
    !
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
        ELSEIF(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4).GT.ZEROS2(N2,N1))THEN
          IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(N2,N1))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)/VcumSnoDWI_col(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
        ELSEIF(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
          IF(VcumSnoDWI_col(N5,N4).GT.ZEROS2(N5,N4))THEN
            VFLW=AZMIN1(AMAX1(-VFLWX,DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)/VcumSnoDWI_col(N5,N4)))
          ELSE
            VFLW=-VFLWX
          ENDIF
        ENDIF

        if(.not.isclose(VFLW,0._r8))then
          DO idg=idg_beg,idg_NH3
            trcg_2DSnowDrift(idg,N,N5,N4)  = VFLW*AZMAX1(trcg_solsml2_snvr(idg,1,N5,N4))
            trcg_FloXSnow_2DH(idg,N,N5,N4) = trcg_FloXSnow_2DH(idg,N,N5,N4)+trcg_2DSnowDrift(idg,N,N5,N4)
          ENDDO

          DO nnut=ids_nut_beg,ids_nuts_end
            trcn_2DSnowDrift(nnut,N,N5,N4)  = VFLW*AZMAX1(trcn_solsml2_snvr(nnut,1,N5,N4))
            trcn_FloXSnow_2DH(nnut,N,N5,N4) = trcn_FloXSnow_2DH(nnut,N,N5,N4)+trcn_2DSnowDrift(nnut,N,N5,N4)
          ENDDO
        endif  

      ENDIF
    ENDDO D4305
  ENDDO D4310
  end subroutine OverLandTracerTranspt
!------------------------------------------------------------------------------------------

  subroutine LitterGasVolatilDissolMM(I,J,M,NY,NX,RGas_Dif_Atm2Litr_FlxMM)
  implicit none

  integer, intent(in) :: I,J, M,NY, NX
  real(r8), intent(in) :: RGas_Dif_Atm2Litr_FlxMM(idg_beg:idg_end-1)
  real(r8) :: VOLGas
  real(r8) :: trc_gascl
  integer  :: idg
!
!     VOLT=litter volume from hour1.f
!     VLWatMicPM,VLsoiAirPM=micropore water volume, air volume from watsub.f
!     VLWatMicP*=equivalent aqueous volume for gas
!     S*L=solubility of gas in water from hour1.f
!     C*G=gas concentration
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     R*DFG=surface gas volatilization
!     *S2=litter solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     R*DXR=gas exchange between atmosphere and surface litter water for gas flux calculations
!

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLsoiAirPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN

    VLWatMicPXA(0,NY,NX)=natomw*VLWatMicPM_vr(M,0,NY,NX)

    do idg=idg_beg,idg_NH3
      trcg_VLWatMicP_vr(idg,0,NY,NX)  = VLWatMicPM_vr(M,0,NY,NX)*GasSolbility_vr(idg,0,NY,NX)
      trc_gascl                       = trc_gascl_vr(idg,0,NY,NX)*VLsoiAirPM_vr(M,0,NY,NX)
      VOLGas                          = trcg_VLWatMicP_vr(idg,0,NY,NX)+VLsoiAirPM_vr(M,0,NY,NX)
      RGas_Disol_FlxM_vr(idg,0,NY,NX) = DiffusivitySolutEffM_vr(M,0,NY,NX) &
        *(AMAX1(ZEROS(NY,NX),trc_gascl)*trcg_VLWatMicP_vr(idg,0,NY,NX) &
        - AMAX1(ZEROS(NY,NX),trc_solml2_vr(idg,0,NY,NX)+RGas_Dif_Atm2Litr_FlxMM(idg))*VLsoiAirPM_vr(M,0,NY,NX))/VOLGas
      Gas_Disol_Flx_vr(idg,0,NY,NX) = Gas_Disol_Flx_vr(idg,0,NY,NX)+RGas_Disol_FlxM_vr(idg,0,NY,NX)
    ENDDO
!
  ELSE
    RGas_Disol_FlxM_vr(idg_beg:idg_end,0,NY,NX)=0.0_r8
  ENDIF
  end subroutine LitterGasVolatilDissolMM

!------------------------------------------------------------------------------------------

  subroutine SurfSoilFluxGasDifAdvMM(M,NY,NX,WaterFlow2Soil,RGas_Dif_Atm2Soil_FlxMM)
!
! DESCRIPTION:
! surface soil gaseous diffusion, advection, dissolution & volatilization
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8),intent(in) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8),intent(in) :: RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)
  real(r8) :: VFLW
  integer :: idg
!
!     AirFilledSoilPoreM_vr=air-filled porosity from watsub.f
!     SoilBulkDensity_vr=bulk density
!
  IF(AirFilledSoilPoreM_vr(M,NU(NY,NX),NY,NX).GT.THETX .AND. SoilBulkDensity_vr(NU(NY,NX),NY,NX).GT.ZERO)THEN
!
    call SurfSoilDifusFluxMM(M,NY,NX)

!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SURFACE WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE SOIL SURFACE
!     OR THE ATMOSPHERE DEPENDING ON WATER FLUX DIRECTION
!
    call SurfSoillAdvFlux(M,NY,NX,WaterFlow2Soil)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLG=hourly convective+diffusive gas flux
!
    DO idg=idg_beg,idg_NH3
      Gas_3DAdvDif_Flx_vr(idg,3,NU(NY,NX),NY,NX)=Gas_3DAdvDif_Flx_vr(idg,3,NU(NY,NX),NY,NX) &
        +RGasADFlx_3D(idg,3,NU(NY,NX),NY,NX)
    ENDDO

    IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      call SurfSoilFluxDisolVapor(M,NY,NX,RGas_Dif_Atm2Soil_FlxMM)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      DO idg=idg_beg,idg_end
        Gas_Disol_Flx_vr(idg,NU(NY,NX),NY,NX)=Gas_Disol_Flx_vr(idg,NU(NY,NX),NY,NX)+RGas_Disol_FlxM_vr(idg,NU(NY,NX),NY,NX)
      ENDDO
    ELSE
      RGas_Disol_FlxM_vr(idg_beg:idg_end,NU(NY,NX),NY,NX)=0.0_r8
    ENDIF
  ELSE
    RGasADFlx_3D(idg_beg:idg_NH3,3,NU(NY,NX),NY,NX)    = 0.0_r8
    RGas_Disol_FlxM_vr(idg_beg:idg_end,NU(NY,NX),NY,NX) = 0.0_r8
  ENDIF
  end subroutine SurfSoilFluxGasDifAdvMM

!------------------------------------------------------------------------------------------

  subroutine LitterAtmosExchangeM(I,J,M,NY,NX,trcs_cl1,RGas_Dif_Atm2Litr_FlxMM,RGasAtmDisol2LitrM)
  implicit none
  integer, intent(in) :: I,J,NY,NX,M
  real(r8),intent(out) :: trcs_cl1(ids_beg:ids_end)
  real(r8),intent(out) :: RGas_Dif_Atm2Litr_FlxMM(idg_beg:idg_end-1)
  real(r8),intent(out) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3)
  integer :: K,ids,idg,idom

  real(r8) :: trc_gasq(idg_beg:idg_NH3)
  real(r8) :: DFGcc(idg_beg:idg_NH3)
  real(r8) :: DLYR0,TORT0

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX)) !vertical layer thickness
    TORT0 = TortMicPM_vr(M,0,NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR0)*FracSurfByLitR_col(NY,NX)

    DO idg=idg_beg,idg_NH3
      DFGcc(idg)=SoluteDifusvtytscal_vr(idg,0,NY,NX)*TORT0
    ENDDO

    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MicP1(idom,K)=AZMAX1(DOM_MicP2(idom,K,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
      enddo
    ENDDO

    !temporary solute concentration
    DO ids=ids_beg,ids_end
      if(ids==idg_O2)then
        trcs_cl1(ids)=AZMAX1(trc_solml2_vr(ids,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
      else
        trcs_cl1(ids)=AZMAX1(trc_solml2_vr(ids,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
      endif
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
!     DiffusivitySolutEff*=effective solute diffusivity
!
    DO idg=idg_beg,idg_NH3
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND RESIDUE SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!     RGasAtmDisol2LitrM: dissolution into water
      !equivalent dissolved gas concentration
      trc_gasq(idg)=(PARR_col(NY,NX)*GasSolbility_vr(idg,0,NY,NX)*AtmGasCgperm3(idg,NY,NX) &
        +DFGcc(idg)*trcs_cl1(idg))/(DFGcc(idg)+PARR_col(NY,NX))

      RGasAtmDisol2LitrM(idg)=(trc_gasq(idg)-trcs_cl1(idg))*AMIN1(VLWatMicPM_vr(M,0,NY,NX),DFGcc(idg))
!
!     ACCUMULATE HOURLY FLUXES, dissolve into water
!      
      trcg_DisolEvap_Atm2Litr_flx(idg,NY,NX)=trcg_DisolEvap_Atm2Litr_flx(idg,NY,NX)+RGasAtmDisol2LitrM(idg)

      RGas_Dif_Atm2Litr_FlxMM(idg)=RGasAtmDisol2LitrM(idg)*dt_GasCyc
    ENDDO

  ELSE
    RGasAtmDisol2LitrM(idg_beg:idg_NH3)=0.0_r8
  ENDIF

  end subroutine LitterAtmosExchangeM
!------------------------------------------------------------------------------------------
  subroutine SoilAtmosExchangeM(I,J,M,NY,NX,trcs_cl2,RGas_Dif_Atm2Soil_FlxMM,RGasAtmDisol2SoilM)

  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcs_cl2(ids_beg:ids_end)
  real(r8), intent(out) :: RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)   !diffusive gas flux from atmosphere to soil
  real(r8), intent(out) :: RGasAtmDisol2SoilM(idg_beg:idg_end)        !gas dissolve from atmosphere to top soil
  character(len=*), parameter :: subname='SoilAtmosExchangeM'

  real(r8) :: DLYR1,TORT1,VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: DiffusivitySolutEff
  real(r8) :: trc_gsolc,trc_gsolc2
  integer  :: K,idg,idom
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE SOIL LAYER
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     TORT=tortuosity from hour1.f
!     *SGL*=solute diffusivity
!     DiffusivitySolutEff*=effective solute diffusivity
!
  call PrintInfo('beg '//subname)
  IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMicPOA = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)
    VLWatMicPOB = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)
    VLWatMicPPA = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    VLWatMicPPB = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
    DLYR1       = AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))
    TORT1       = TortMicPM_vr(M,NU(NY,NX),NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR1)

    D8910: DO  K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MicP2(idom,K)=AZMAX1(DOM_MicP2(idom,K,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX))
      enddo
    ENDDO D8910

   !excldue  NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      trcs_cl2(idg)=trc_solml2_vr(idg,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)
    ENDDO

    IF(VLWatMicPMA_vr(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl2(idg_NH3) = AZMAX1(trc_solml2_vr(idg_NH3,NU(NY,NX),NY,NX)/VLWatMicPMA_vr(NU(NY,NX),NY,NX))
      trcs_cl2(ids_NH4) = AZMAX1(trc_solml2_vr(ids_NH4,NU(NY,NX),NY,NX)/VLWatMicPMA_vr(NU(NY,NX),NY,NX))
    ELSE
      trcs_cl2(idg_NH3) = 0.0_r8
      trcs_cl2(ids_NH4) = 0.0_r8
    ENDIF

    IF(VLWatMicPOA.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_NO3) = AZMAX1(trc_solml2_vr(ids_NO3,NU(NY,NX),NY,NX)/VLWatMicPOA)
      trcs_cl2(ids_NO2) = AZMAX1(trc_solml2_vr(ids_NO2,NU(NY,NX),NY,NX)/VLWatMicPOA)
    ELSE
      trcs_cl2(ids_NO3) = 0.0_r8
      trcs_cl2(ids_NO2) = 0.0_r8
    ENDIF
    IF(VLWatMicPPA.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_H1PO4) = AZMAX1(trc_solml2_vr(ids_H1PO4,NU(NY,NX),NY,NX)/VLWatMicPPA)
      trcs_cl2(ids_H2PO4) = AZMAX1(trc_solml2_vr(ids_H2PO4,NU(NY,NX),NY,NX)/VLWatMicPPA)
    ELSE
      trcs_cl2(ids_H1PO4) = 0.0_r8
      trcs_cl2(ids_H2PO4) = 0.0_r8
    ENDIF
    IF(VLWatMicPMB_vr(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl2(idg_NH3B) = AZMAX1(trc_solml2_vr(idg_NH3B,NU(NY,NX),NY,NX)/VLWatMicPMB_vr(NU(NY,NX),NY,NX))
      trcs_cl2(ids_NH4B) = AZMAX1(trc_solml2_vr(ids_NH4B,NU(NY,NX),NY,NX)/VLWatMicPMB_vr(NU(NY,NX),NY,NX))
    ELSE
      trcs_cl2(idg_NH3B) = trcs_cl2(idg_NH3)
      trcs_cl2(ids_NH4B) = trcs_cl2(ids_NH4)
    ENDIF

    IF(VLWatMicPOB.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_NO3B) = AZMAX1(trc_solml2_vr(ids_NO3B,NU(NY,NX),NY,NX)/VLWatMicPOB)
      trcs_cl2(ids_NO2)  = AZMAX1(trc_solml2_vr(ids_NO2B,NU(NY,NX),NY,NX)/VLWatMicPOB)
    ELSE
      trcs_cl2(ids_NO3B) = trcs_cl2(ids_NO3)
      trcs_cl2(ids_NO2)  = trcs_cl2(ids_NO2)
    ENDIF
    IF(VLWatMicPPB.GT.ZEROS2(NY,NX))THEN
      trcs_cl2(ids_H1PO4B) = AZMAX1(trc_solml2_vr(ids_H1PO4B,NU(NY,NX),NY,NX)/VLWatMicPPB)
      trcs_cl2(ids_H2PO4B) = AZMAX1(trc_solml2_vr(ids_H2PO4B,NU(NY,NX),NY,NX)/VLWatMicPPB)
    ELSE
      trcs_cl2(ids_H1PO4B) = trcs_cl2(ids_H1PO4)
      trcs_cl2(ids_H2PO4B) = trcs_cl2(ids_H2PO4)
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
!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND SOIL SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!include NH3B
    DO idg=idg_beg,idg_end
      DiffusivitySolutEff = SoluteDifusvtytscal_vr(idg,NU(NY,NX),NY,NX)*TORT1
      trc_gsolc           = (CondGasXSnowM_col(M,NY,NX)*AtmGasCgperm3(idg,NY,NX)*GasSolbility_vr(idg,NU(NY,NX),NY,NX) &
        +DiffusivitySolutEff*trcs_cl2(idg))/(DiffusivitySolutEff+CondGasXSnowM_col(M,NY,NX))

      RGasAtmDisol2SoilM(idg)=(trc_gsolc-trcs_cl2(idg)) &
        *AMIN1(VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(idg,NU(NY,NX),NY,NX),DiffusivitySolutEff)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!     seven gas species plus aqueous NH3 in band

      trcg_DisolEvap_Atm2Soil_flx(idg,NY,NX) = trcg_DisolEvap_Atm2Soil_flx(idg,NY,NX)+RGasAtmDisol2SoilM(idg)    

      RGas_Dif_Atm2Soil_FlxMM(idg)=RGasAtmDisol2SoilM(idg)*dt_GasCyc
    ENDDO

  ELSE
    RGasAtmDisol2SoilM(idg_beg:idg_end)=0.0_r8
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine SoilAtmosExchangeM
!------------------------------------------------------------------------------------------

  subroutine SnowSoluteDischargeM(M,NY,NX,trcg_FloSno2LitR,trcn_FloSno2LitR)

  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8), intent(out) :: trcg_FloSno2LitR(idg_beg:idg_NH3)
  real(r8), intent(out) :: trcn_FloSno2LitR(ids_nut_beg:ids_nuts_end)
  integer :: ICHKL  !flag for dealing with the bottom snow layer
  integer :: L,L2
  real(r8) :: VFLWS,VFLWW,VFLWR
  real(r8) :: VFLWNOB,VFLWNO3,VFLWNHB
  real(r8) :: VFLWNH4,VFLWPO4,VFLWPOB
  integer  :: idg,nnut
!
!     VLSnowHeatCapM,VLHeatCapSnowMin_col=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     WatFlowInSnowM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  ICHKL=0
  DO  L=1,JS
    IF(VLSnowHeatCapM_snvr(M,L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      L2=MIN(JS,L+1)
      IF(L.LT.JS .AND. VLSnowHeatCapM_snvr(M,L2,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
        !water flow from L to L2
        IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0_r8,WatFlowInSnowM_snvr(M,L2,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
          !no water is layer L, all water is gone to L2
        ELSE          
          VFLWW=1.0_r8
        ENDIF
        DO idg=idg_beg,idg_NH3
          trcg_advW_snvr(idg,L2,NY,NX) = trcg_solsml2_snvr(idg,L,NY,NX)*VFLWW
          trcVolatile_Xbndl_flx_snvr(idg,L2,NY,NX) = trcVolatile_Xbndl_flx_snvr(idg,L2,NY,NX)+trcg_advW_snvr(idg,L2,NY,NX)
        ENDDO

        DO nnut=ids_nut_beg,ids_nuts_end
          trcn_RBLS(nnut,L2,NY,NX)=trcn_solsml2_snvr(nnut,L,NY,NX)*VFLWW
          trcn_Xbndl_flx(nnut,L2,NY,NX)=trcn_Xbndl_flx(nnut,L2,NY,NX)+trcn_RBLS(nnut,L2,NY,NX)
        ENDDO
      ELSE
        IF(L.LT.JS)THEN
          trcg_advW_snvr(idg_beg:idg_NH3,L2,NY,NX)=0.0_r8
          trcn_RBLS(ids_nut_beg:ids_nuts_end,L2,NY,NX)=0.0_r8
        ENDIF
!
!     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
!
!     VOLWSL=snowpack water content
!     WatFlowSno2LitRM_col,WatFlowSno2MicPM_col,WatFlowSno2MacPM_col=total water flux to litter,soil micropore,macropore
!     CVRD,BARE=litter cover fraction,1-CVRD
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     R*S0=solute flux to surface litter
!     R*S1,R*B1=solute flux to soil surface non-band,band
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
        IF(ICHKL.EQ.0)THEN
          IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,WatFlowSno2LitRM_col(M,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX))/VLWatSnow_snvr(L,NY,NX)))
          ELSE
            VFLWR=FracSurfByLitR_col(NY,NX)
            VFLWS=FracSurfBareSoil_col(NY,NX)
          ENDIF
          VFLWNH4=VFLWS*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
          VFLWNHB=VFLWS*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
          VFLWNO3=VFLWS*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)
          VFLWNOB=VFLWS*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)
          VFLWPO4=VFLWS*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

          DO idg=idg_beg,idg_NH3
            trcg_FloSno2LitR(idg)=trcg_solsml2_snvr(idg,L,NY,NX)*VFLWR
          ENDDO

          DO nnut=ids_nut_beg,ids_nuts_end
            trcn_FloSno2LitR(nnut)=trcn_solsml2_snvr(nnut,L,NY,NX)*VFLWR
          ENDDO

          trcg_VFloSnow(idg_CO2)  = trcg_solsml2_snvr(idg_CO2,L,NY,NX)*VFLWS
          trcg_VFloSnow(idg_CH4)  = trcg_solsml2_snvr(idg_CH4,L,NY,NX)*VFLWS
          trcg_VFloSnow(idg_O2)   = trcg_solsml2_snvr(idg_O2,L,NY,NX)*VFLWS
          trcg_VFloSnow(idg_N2)   = trcg_solsml2_snvr(idg_N2,L,NY,NX)*VFLWS
          trcg_VFloSnow(idg_N2O)  = trcg_solsml2_snvr(idg_N2O,L,NY,NX)*VFLWS
          trcg_VFloSnow(idg_NH3)  = trcg_solsml2_snvr(idg_NH3,L,NY,NX)*VFLWNH4
          trcg_VFloSnow(idg_NH3B) = trcg_solsml2_snvr(idg_NH3,L,NY,NX)*VFLWNHB

          trcn_soil_VFloSnow(ids_NH4)    = trcn_solsml2_snvr(ids_NH4,L,NY,NX)*VFLWNH4
          trcn_soil_VFloSnow(ids_NO3)    = trcn_solsml2_snvr(ids_NO3,L,NY,NX)*VFLWNO3
          trcn_soil_VFloSnow(ids_H1PO4)  = trcn_solsml2_snvr(ids_H1PO4,L,NY,NX)*VFLWPO4
          trcn_soil_VFloSnow(ids_H2PO4)  = trcn_solsml2_snvr(ids_H2PO4,L,NY,NX)*VFLWPO4
          trcn_band_VFloSnow(ids_NH4B)   = trcn_solsml2_snvr(ids_NH4,L,NY,NX)*VFLWNHB
          trcn_band_VFloSnow(ids_NO3B)   = trcn_solsml2_snvr(ids_NO3,L,NY,NX)*VFLWNOB
          trcn_band_VFloSnow(ids_H1PO4B) = trcn_solsml2_snvr(ids_H1PO4,L,NY,NX)*VFLWPOB
          trcn_band_VFloSnow(ids_H2PO4B) = trcn_solsml2_snvr(ids_H2PO4,L,NY,NX)*VFLWPOB
          ICHKL                          = 1
        ENDIF
      ENDIF
    ELSE
      trcg_FloSno2LitR(idg_beg:idg_NH3)          = 0.0_r8
      trcn_FloSno2LitR(ids_nut_beg:ids_nuts_end) = 0.0_r8
      trcg_VFloSnow(idg_beg:idg_end)             = 0.0_r8

      trcn_soil_VFloSnow(ids_NH4)   = 0.0_r8
      trcn_soil_VFloSnow(ids_NO3)   = 0.0_r8
      trcn_soil_VFloSnow(ids_H1PO4) = 0.0_r8
      trcn_soil_VFloSnow(ids_H2PO4) = 0.0_r8
      trcn_band_VFloSnow(ids_NH4B)  = 0.0_r8

      trcn_band_VFloSnow(ids_NO3B)   = 0.0_r8
      trcn_band_VFloSnow(ids_H1PO4B) = 0.0_r8
      trcn_band_VFloSnow(ids_H2PO4B) = 0.0_r8
    ENDIF
  ENDDO
  end subroutine SnowSoluteDischargeM

! ----------------------------------------------------------------------
  subroutine SurfSoillAdvFlux(M,NY,NX,WaterFlow2Soil)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8) :: VFLW
  real(r8) :: RFLg_ADV(idg_beg:idg_NH3)
  integer :: idg

!     WaterFlow2Soil=total water flux into soil micropore+macropore from watsub.f
!     VLsoiAirPM=air-filled porosity
!     RFL*G=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     *G2=gaseous content
!     C*E=atmospheric gas concentration from hour1.f
!

  IF(WaterFlow2Soil(3,NU(NY,NX),NY,NX).GT.0.0_r8)THEN
    IF(VLsoiAirPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,WaterFlow2Soil(3,NU(NY,NX),NY,NX)/VLsoiAirPM_vr(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO idg=idg_beg,idg_NH3
      RFLg_ADV(idg)=VFLW*AZMAX1(trc_gasml2_vr(idg,NU(NY,NX),NY,NX))
    ENDDO
  ELSE
    DO idg=idg_beg,idg_NH3
      RFLg_ADV(idg)=-WaterFlow2Soil(3,NU(NY,NX),NY,NX)*AtmGasCgperm3(idg,NY,NX)
    ENDDO
  ENDIF
!     TOTAL SOIL GAS FLUX + CONVECTIVE FLUX
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux
  DO idg=idg_beg,idg_NH3
    RGasADFlx_3D(idg,3,NU(NY,NX),NY,NX)=RGasADFlx_3D(idg,3,NU(NY,NX),NY,NX)+RFLg_ADV(idg)
  ENDDO
  end subroutine SurfSoillAdvFlux

! ----------------------------------------------------------------------
  subroutine SurfSoilDifusFluxMM(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: DFV_g
  real(r8) :: DFLG2
  real(r8) :: trcg_cl2
  real(r8) :: DGQ_cef
  integer  :: idg
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
  DFLG2=AZMAX1(AirFilledSoilPoreM_vr(M,NU(NY,NX),NY,NX))*POROQ &
    *AirFilledSoilPoreM_vr(M,NU(NY,NX),NY,NX)/POROS_vr(NU(NY,NX),NY,NX) &
    *AREA(3,NU(NY,NX),NY,NX)/AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))

  DO idg=idg_beg,idg_NH3
    DifuscG_vr(idg,3,NU(NY,NX),NY,NX)=DFLG2*GasDifc_vrc(idg,NU(NY,NX),NY,NX)
!
!     SURFACE GAS CONCENTRATIONS
!
!     C*G2=gaseous concentration
!     *G2=gaseous content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     VLsoiAirPM=air-filled porosity
!
    trcg_cl2=AZMAX1(trc_gasml2_vr(idg,NU(NY,NX),NY,NX)/VLsoiAirPM_vr(M,NU(NY,NX),NY,NX))
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
    DGQ_cef=DifuscG_vr(idg,3,NU(NY,NX),NY,NX)*PARG_cef(idg,NY,NX) &
      /(DifuscG_vr(idg,3,NU(NY,NX),NY,NX)+PARG_cef(idg,NY,NX))

    DFV_g=DGQ_cef*(AtmGasCgperm3(idg,NY,NX)-trcg_cl2)

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux

    RGasADFlx_3D(idg,3,NU(NY,NX),NY,NX)=DFV_g
  ENDDO

  end subroutine SurfSoilDifusFluxMM
! ----------------------------------------------------------------------

  subroutine SurfSoilFluxDisolVapor(M,NY,NX,RGas_Dif_Atm2Soil_FlxMM)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)
  real(r8) :: trcg_VOLT(idg_beg:idg_end,JY,JX)
  integer :: idg

!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     VLWatMicP*=equivalent aqueous volume for gas
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     S*L=solubility of gas in water from hour1.f
!     VLsoiAirPM=air-filled porosity
!     R*DFG=water-air gas flux
!     DiffusivitySolutEff=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     R*DXS=gas exchange between atmosphere and soil surface water
!
  do idg=idg_beg,idg_NH3-1
    trcg_VLWatMicP_vr(idg,NU(NY,NX),NY,NX) = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*GasSolbility_vr(idg,NU(NY,NX),NY,NX)
    trcg_VOLT(idg,NY,NX)                   = trcg_VLWatMicP_vr(idg,NU(NY,NX),NY,NX)+VLsoiAirPM_vr(M,NU(NY,NX),NY,NX)
    RGas_Disol_FlxM_vr(idg,NU(NY,NX),NY,NX) = DiffusivitySolutEffM_vr(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2_vr(idg,NU(NY,NX),NY,NX))*trcg_VLWatMicP_vr(idg,NU(NY,NX),NY,NX) &
       -AMAX1(ZEROS(NY,NX),trc_solml2_vr(idg,NU(NY,NX),NY,NX)+RGas_Dif_Atm2Soil_FlxMM(idg))*VLsoiAirPM_vr(M,NU(NY,NX),NY,NX)) &
       /trcg_VOLT(idg,NY,NX)
  enddo

  trcg_VLWatMicP_vr(idg_NH3,NU(NY,NX),NY,NX)  = VLWatMicPMA_vr(NU(NY,NX),NY,NX)*GasSolbility_vr(idg_NH3,NU(NY,NX),NY,NX)
  trcg_VLWatMicP_vr(idg_NH3B,NU(NY,NX),NY,NX) = VLWatMicPMB_vr(NU(NY,NX),NY,NX)*GasSolbility_vr(idg_NH3,NU(NY,NX),NY,NX)
  trcg_VOLT(idg_NH3,NY,NX)                    = trcg_VLWatMicP_vr(idg_NH3,NU(NY,NX),NY,NX)+VLsoiAirPMA(NU(NY,NX),NY,NX)
  trcg_VOLT(idg_NH3B,NY,NX)                   = trcg_VLWatMicP_vr(idg_NH3B,NU(NY,NX),NY,NX)+VLsoiAirPMB(NU(NY,NX),NY,NX)

  IF(trcg_VOLT(idg_NH3,NY,NX).GT.ZEROS2(NY,NX) .AND. VLWatMicPXA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    RGas_Disol_FlxM_vr(idg_NH3,NU(NY,NX),NY,NX)=DiffusivitySolutEffM_vr(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2_vr(idg_NH3,NU(NY,NX),NY,NX))*trcg_VLWatMicP_vr(idg_NH3,NU(NY,NX),NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2_vr(idg_NH3,NU(NY,NX),NY,NX)+RGas_Dif_Atm2Soil_FlxMM(idg_NH3))*VLsoiAirPMA(NU(NY,NX),NY,NX)) &
      /trcg_VOLT(idg_NH3,NY,NX)
  ELSE
    RGas_Disol_FlxM_vr(idg_NH3,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF

  IF(trcg_VOLT(idg_NH3B,NY,NX).GT.ZEROS2(NY,NX) .AND. VLWatMicPXB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    RGas_Disol_FlxM_vr(idg_NH3B,NU(NY,NX),NY,NX)=DiffusivitySolutEffM_vr(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),trc_gasml2_vr(idg_NH3,NU(NY,NX),NY,NX))*trcg_VLWatMicP_vr(idg_NH3B,NU(NY,NX),NY,NX) &
      -AMAX1(ZEROS(NY,NX),trc_solml2_vr(idg_NH3B,NU(NY,NX),NY,NX)+RGas_Dif_Atm2Soil_FlxMM(idg_NH3B))*VLsoiAirPMB(NU(NY,NX),NY,NX)) &
      /trcg_VOLT(idg_NH3B,NY,NX)
  ELSE
    RGas_Disol_FlxM_vr(idg_NH3B,NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  end subroutine SurfSoilFluxDisolVapor
! ----------------------------------------------------------------------
end module SurfaceFluxMod
