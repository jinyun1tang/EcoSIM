module IngridTranspMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use TranspSaltDataMod
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use GridDataType
  use SurfLitterDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SnowDataType
  use ChemTranspDataType
  use EcoSIMSolverPar
  use minimathmod, only : AZMAX1,AZMIN1
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), PARAMETER :: XFRS=0.05_r8
  public :: SaltModelSoluteHydroFlux
  contains


!------------------------------------------------------------------------------------------

  subroutine SaltModelSoluteHydroFlux(I,M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,M,NHW,NHE,NVN,NVS

  integer :: NY,NX,N1,N2
  real(r8) :: FLWRM1
  real(r8) :: trcSaltSnoFlo2Soil(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  real(r8) :: trcSaltSnoFlo2LitR(idsalt_beg:idsalt_end)
  real(r8) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
!     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     ADD SINKS FROM SOLUTE.F
      call SoluteSinksInSoil(NY,NX)
!
!     SOLUTE FLUXES FROM MELTING SNOWPACK TO
!     SOIL SURFACE FROM SNOWMELT IN 'WATSUB' AND
!     CONCENTRATIONS IN SNOWPACK
      call SoluteFluxInSnowpack(M,NY,NX,trcSaltSnoFlo2Soil,trcSaltSnoFlo2LitR)
!
!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
      FLWRM1=WatFLo2LitrM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f

      IF(FLWRM1.GT.0.0_r8)THEN
        call Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1,trcSalt_RFL)
!
      ELSE
        call Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1,trcSalt_RFL)
      ENDIF
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!
      IF((VGeomLayer_vr(0,NY,NX).GT.ZEROS(NY,NX).AND.VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
        .AND.(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN

        call TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1,trcSalt_flx_diffus)
      ELSE
        trcSalt_flx_diffus(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES

      call TopsoilResidueFluxAdvPlusDifus(NY,NX,trcSaltSnoFlo2Soil,&
        trcSalt_RFL,trcSaltSnoFlo2LitR,trcSalt_flx_diffus)
!
      call AccumHourlyTopsoilReisdueFlux(NY,NX,trcSalt_RFL,trcSalt_flx_diffus,&
        trcSaltSnoFlo2LitR,trcSaltSnoFlo2Soil)
!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE IN SOIL
!     SURFACE LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FWatExMacP2MicPM=macro-micropore water transfer from watsub.f
!     VLWatMicPM,VLWatMacPM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
      IF(FWatExMacP2MicPM(M,NU(NY,NX),NY,NX).GT.0.0)THEN
        call MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
      ELSEIF(FWatExMacP2MicPM(M,NU(NY,NX),NY,NX).LT.0.0)THEN
        call MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
      ELSE
        trcSalt_RFL(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
      IF(VLWatMacPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
        call MacMicPoreSoluteDifusExchange(M,NY,NX,trcSalt_flx_diffus)
      ELSE
        trcSalt_flx_diffus(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!      call MacMicPoreFluxAdvPlusDifus(NY,NX,trcSalt_flx_diffus,trcSalt_RFL)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
      call AccumHourlyMicMacPoreFlux(NY,NX)
!
!     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
!     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE LAYER
!
!     QRM=runoff from watsub.f
!     RQR*0=solute in runoff
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLWatMicPM=litter water volume from watsub.f
!     *S2=litter solute content
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     X*QRS=accumulated hourly solute in runoff
!
!
      N1=NX;N2=NY

      call SoluteFluxBySurfaceOutflow(M,N1,N2)
!
      call UpdateSoluteInSurfNeighbors(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
      call UpdateSoluteInSubsurfNeighbors(M,NY,NX,NHW,NHE,NVN,NVS)

    ENDDO
  ENDDO
  end subroutine SaltModelSoluteHydroFlux

!------------------------------------------------------------------------------------------

  subroutine SoluteSinksInSoil(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,nsalts
!     begin_execution
!
!     RZ*2=solute flux at time step for flux calculations
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO L=NU(NY,NX),NL(NY,NX)
    DO nsalts=idsalt_beg,idsaltb_end
      trcSalt_solml2(nsalts,L,NY,NX)=trcSalt_solml2(nsalts,L,NY,NX)-trcSalt_solml2R(nsalts,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine SoluteSinksInSoil
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxInSnowpack(M,NY,NX,trcSaltSnoFlo2Soil,trcSaltSnoFlo2LitR)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L,ICHKL,L2,nsalts,ids
  real(r8) :: VFLWW,VFLWR,VFLWS,VFLWPO4,VFLWPOB
  real(r8), intent(out) :: trcSaltSnoFlo2Soil(idsalt_beg:idsaltb_end)
  real(r8), intent(out) :: trcSaltSnoFlo2LitR(idsalt_beg:idsalt_end)
  
!     begin_execution
!
!     VLSnowHeatCapM,VLHeatCapSnowMin_col=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     WatFlowInSnowM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack
! VFLWW=fraction flowed into next snow layer
  ICHKL=0
  DO L=1,JS
    !
    IF(VLSnowHeatCapM_snvr(M,L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VLSnowHeatCapM_snvr(M,L2,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
        IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0_r8,WatFlowInSnowM_snvr(M,L2,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
        ELSE
          VFLWW=1.0_r8
        ENDIF

        DO nsalts=idsalt_beg,idsalt_end
          trcSaltAdv2SowLay(nsalts,L2,NY,NX)=trc_Saltml2_snvr(nsalts,L,NY,NX)*VFLWW
          trcSaltFlo2SnowLay(nsalts,L2,NY,NX)=trcSaltFlo2SnowLay(nsalts,L2,NY,NX)+trcSaltAdv2SowLay(nsalts,L2,NY,NX)
        ENDDO

      ELSE
        !bottom snow layer or insignificant snow layer
        IF(L.LT.JS)THEN
          !insignificant snow layer
          trcSaltAdv2SowLay(idsalt_beg:idsalt_end,L2,NY,NX)=0.0_r8
        ENDIF
!
!     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
!
!     VOLWSL=snowpack water content
!     WatFlowSno2LitRM,WatFlowSno2MicPM,WatFlowSno2MacPM=total water flux to litter,soil micropore,macropore
!     CVRD,BARE=litter cover fraction,1-CVRD
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     R*S0=solute flux to surface litter
!     R*S1,R*B1=solute flux to soil surface non-band,band
!
        IF(ICHKL.EQ.0)THEN
          !flow into soil
          IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,WatFlowSno2LitRM(M,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(WatFlowSno2MicPM(M,NY,NX)+WatFlowSno2MacPM(M,NY,NX))/VLWatSnow_snvr(L,NY,NX)))
          ELSE
            VFLWR=FracSurfByLitR_col(NY,NX)
            VFLWS=FracSurfBareSoil_col(NY,NX)
          ENDIF
          VFLWPO4=VFLWS*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

          DO nsalts=idsalt_beg,idsalt_end
            trcSaltSnoFlo2LitR(nsalts)=trc_Saltml2_snvr(nsalts,L,NY,NX)*VFLWR
          ENDDO

          DO nsalts=idsalt_beg,idsalt_KSO4
            trcSaltSnoFlo2Soil(nsalts)=trc_Saltml2_snvr(nsalts,L,NY,NX)*VFLWS
          ENDDO

          DO nsalts=idsalt_H0PO4,idsalt_MgHPO4
            ids=nsalts-idsalt_H0PO4+idsalt_H0PO4B  
            trcSaltSnoFlo2Soil(nsalts)=trc_Saltml2_snvr(nsalts,L,NY,NX)*VFLWPO4
            trcSaltSnoFlo2Soil(ids)=trc_Saltml2_snvr(ids,L,NY,NX)*VFLWPOB
          ENDDO
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      trcSaltSnoFlo2LitR(idsalt_beg:idsalt_end)=0.0_r8
      trcSaltSnoFlo2Soil(idsalt_beg:idsaltb_end)=0.0_r8
    ENDIF
  ENDDO
  end subroutine SoluteFluxInSnowpack
!------------------------------------------------------------------------------------------

  subroutine Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1,trcSalt_RFL)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW
  integer :: nsalts,ids

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!     VLWatMicPM=litter water volume
!     RFL*=soil-litter convective solute flux
!     Z*2=litter solute content
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  IF(VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VLWatMicPM_vr(M,0,NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF
  DO nsalts=idsalt_beg,idsalt_KSO4
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,0,NY,NX))
  ENDDO

  DO nsalts=idsalt_H0PO4,idsalt_MgHPO4
    ids=nsalts-idsalt_H0PO4+idsalt_H0PO4B  
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,0,NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcSalt_RFL(ids)=VFLW*AZMAX1(trcSalt_solml2(ids,0,NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine Residue2TopsoilSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1,trcSalt_RFL)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  integer :: nsalts
  real(r8) :: VFLW

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE
!
!     VLWatMicPM=litter water volume
!     RFL*=soil-litter convective solute flux
!     Z*2=soil solute content
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
  IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX))
  ENDDO
  end subroutine Topsoil2ResidueSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1,trcSalt_flx_diffus)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)    !tracer solute concentration in source cell
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)   !tracer solutte concentration in destination cell
  real(r8) :: VOLWPB,VOLWPA,TORT0,TORT1
  real(r8) :: DLYR0
  integer :: nsalts
!     begin_execution
!
!     MICROPORE CONCENTRATIONS FROM WATER IN RESIDUE AND SOIL SURFACE
!
!     C*1,C*2=solute concentration in litter, soil
!     Z*1,Z*2=solute content in litter, soil
!
  VOLWPA=VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  VOLWPB=VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

  DO nsalts=idsalt_beg,idsalt_end
    trcSalt_conc_src(nsalts)=AZMAX1(trcSalt_solml2(nsalts,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
  ENDDO

  DO nsalts=idsalt_beg,idsalt_KSO4
    trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX))
  ENDDO

  IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX)/VOLWPA)
    ENDDO
  ELSE
    trcSalt_conc_dest(idsalt_psoiL_beg:idsalt_psoil_end)=0.0_r8
  ENDIF
  
  IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX)/VOLWPB)
    ENDDO
  ELSE
    DO nsalts=0,idsalt_nuts
      trcSalt_conc_dest(idsalt_H0PO4B+nsalts)=trcSalt_conc_dest(idsalt_H0PO4+nsalts)
    ENDDO
  ENDIF
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
!
!     DLYR0,DLYR1=litter, soil surface thickness
!     TORT=tortuosity from hour1.f
!     CVRD,BARE=litter cover fraction,1-CVRD
!     DISP=dispersivity parameter
!     FLWRM=litter-soil water flux from watsub.f
!     DIF*=aqueous diffusivity-dispersivity in soil surface
!     *SGL*=solute diffusivity from hour1.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!
  DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX))
  TORT0=TortMicPM_vr(M,0,NY,NX)*FracSurfByLitR_col(NY,NX)
  DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
  TORT1=TortMicPM_vr(M,NU(NY,NX),NY,NX)
  TORTL=AMIN1(1.0_r8,(TORT0+TORT1)/(DLYR0+DLYR1))
  DISPN=DISP(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))

  DIFPO=(POSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)

  DIFAL=(AquaIonDifusivty2_vr(idsalt_Al,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance((/idsalt_Al,idsalt_AlOH,idsalt_AlOH2,idsalt_AlOH3,idsalt_AlOH4,idsalt_AlSO4/))=DIFAL

  DIFFE=(AquaIonDifusivty2_vr(idsalt_Fe,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance((/idsalt_Fe,idsalt_FeOH,idsalt_FeOH2,idsalt_FeOH3,idsalt_FeOH4,idsalt_FeSO4/))=DIFFE

  trcSaltDiffConductance(idsalt_Hp)=(AquaIonDifusivty2_vr(idsalt_Hp,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)

  DIFCA=(AquaIonDifusivty2_vr(idsalt_Ca,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance((/idsalt_Ca,idsalt_CaOH,idsalt_CO3,idsalt_HCO3,idsalt_CaSO4/))=DIFCA

  DIFMG=(AquaIonDifusivty2_vr(idsalt_Mg,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance((/idsalt_Mg,idsalt_MgOH2,idsalt_MgCO3,idsalt_MgHCO3,idsalt_SO4/))=DIFMG

  DIFNA=(AquaIonDifusivty2_vr(idsalt_Na,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance((/idsalt_Na,idsalt_NaCO3,idsalt_NaSO4/))=DIFNA

  DIFKA=(AquaIonDifusivty2_vr(idsalt_K,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance((/idsalt_K,idsalt_KSO4/))=DIFKA

  trcSaltDiffConductance(idsalt_OH)=(AquaIonDifusivty2_vr(idsalt_OH,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance(idsalt_SO4)=(AquaIonDifusivty2_vr(idsalt_SO4,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance(idsalt_Cl)=(AquaIonDifusivty2_vr(idsalt_Cl,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance(idsalt_CO3)=(AquaIonDifusivty2_vr(idsalt_CO3,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcSaltDiffConductance(idsalt_HCO3)=(AquaIonDifusivty2_vr(idsalt_HCO3,NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!     DFV*S,DFV*B=diffusive solute flux between litter and soil surface in non-band,band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     C*1,C*2=solute concentration in litter,soil surface
!

  DO nsalts=idsalt_beg,idsalt_KSO4
    trcSalt_flx_diffus(nsalts)=trcSaltDiffConductance(nsalts)*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))
  ENDDO

  DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_flx_diffus(nsalts)=DIFPO*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO nsalts=0,idsalt_nuts
    trcSalt_flx_diffus(idsalt_H0PO4B+nsalts)=DIFPO*(trcSalt_conc_src(idsalt_H0PO4+nsalts)&
      -trcSalt_conc_dest(idsalt_H0PO4B+nsalts))*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine TopsoilResidueSolutedifusExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueFluxAdvPlusDifus(NY,NX,trcSaltSnoFlo2Soil,&
    trcSalt_RFL,trcSaltSnoFlo2LitR,trcSalt_flx_diffus)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: trcSaltSnoFlo2Soil(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSaltSnoFlo2LitR(idsalt_beg:idsalt_end)
  real(r8), intent(in) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
  integer :: nsalts
!     begin_execution
!
!     R*FLS=convective + diffusive solute flux between litter, soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into soil in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     R*FL0,R*FL1=convective flux into surface litter, soil surface
!     RFL*=convective flux between surface litter and soil surface
!     DFV*=diffusive solute flux between litter and soil surface
!
  DO nsalts=idsalt_beg,idsalt_end
    trcSalt3DFlo2CellM(nsalts,3,0,NY,NX)=trcSalt_RFL0(nsalts,NY,NX)+trcSaltSnoFlo2LitR(nsalts) &
      -trcSalt_RFL(nsalts)-trcSalt_flx_diffus(nsalts)
  ENDDO

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt3DFlo2CellM(nsalts,3,NU(NY,NX),NY,NX)=trcSalt_RFL1(nsalts,NY,NX) &
      +trcSaltSnoFlo2Soil(nsalts)+trcSalt_RFL(nsalts)+trcSalt_flx_diffus(nsalts)
  ENDDO

  end subroutine TopsoilResidueFluxAdvPlusDifus
!------------------------------------------------------------------------------------------

  subroutine AccumHourlyTopsoilReisdueFlux(NY,NX,trcSalt_RFL,&
    trcSalt_flx_diffus,trcSaltSnoFlo2LitR,trcSaltSnoFlo2Soil)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8), intent(in) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSaltSnoFlo2LitR(idsalt_beg:idsalt_end)
  real(r8), intent(in) :: trcSaltSnoFlo2Soil(idsalt_beg:idsaltb_end)
  integer :: nsalts
!     begin_execution
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  DO nsalts=idsalt_beg,idsalt_end
    trcSalt3DFlo2Cell(nsalts,3,0,NY,NX)=trcSalt3DFlo2Cell(nsalts,3,0,NY,NX)+trcSaltSnoFlo2LitR(nsalts) &
      -trcSalt_RFL(nsalts)-trcSalt_flx_diffus(nsalts)
  ENDDO

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt3DFlo2Cell(nsalts,3,NU(NY,NX),NY,NX)=trcSalt3DFlo2Cell(nsalts,3,NU(NY,NX),NY,NX) &
      +trcSaltSnoFlo2Soil(nsalts)+trcSalt_RFL(nsalts)+trcSalt_flx_diffus(nsalts)
  ENDDO
  end subroutine AccumHourlyTopsoilReisdueFlux
!------------------------------------------------------------------------------------------

  subroutine MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW
  integer :: nsalts
!     begin_execution

  IF(VLWatMacPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM(M,NU(NY,NX),NY,NX)/VLWatMacPM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF

  DO nsalts=idsalt_beg,idsalt_KSO4
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,NU(NY,NX),NY,NX))
  ENDDO

  DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO nsalts=idsalt_pband_beg,idsalt_pband_end
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFXS(nsalts,NU(NY,NX),NY,NX)=trcSalt_RFL(nsalts)
  ENDDO
  end subroutine MacToMicPoreSoluteAdvExchange
!------------------------------------------------------------------------------------------

  subroutine MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  integer :: nsalts
  real(r8) :: VFLW
!     begin_execution
  IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM(M,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO nsalts=idsalt_beg,idsalt_KSO4
    trcSalt_RFL(idsalt_Al)=VFLW*AZMAX1(trcSalt_solml2(idsalt_Al,NU(NY,NX),NY,NX))
  ENDDO

  DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO nsalts=idsalt_pband_beg,idsalt_pband_end
    trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFXS(nsalts,NU(NY,NX),NY,NX)=trcSalt_RFL(nsalts)
  ENDDO
  end subroutine MicToMacPoreSoluteAdvExchange

!------------------------------------------------------------------------------------------

  subroutine MacMicPoreSoluteDifusExchange(M,NY,NX,trcSalt_flx_diffus)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
  real(r8) :: VOLWHS,VOLWT,VLWatMicPMNU
  integer :: nsalts
!     begin_execution
!
!     VLWatMicPM,VLWatMacPM=micropore,macropore water volume
!     XFRS*VOLT=maximum macropore volume for solute transfer
!     DFV*=diffusive macropore-micropore solute transfer
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *H2,*2=macropore,micropore solute content
  VLWatMicPMNU=VLWatMicPM_vr(M,NU(NY,NX),NY,NX)
  VOLWHS=AMIN1(XFRS*VGeomLayer_vr(NU(NY,NX),NY,NX),VLWatMacPM(M,NU(NY,NX),NY,NX))
  VOLWT=VLWatMicPM_vr(M,NU(NY,NX),NY,NX)+VOLWHS

  DO nsalts=idsalt_beg,idsalt_KSO4
    trcSalt_flx_diffus(nsalts)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2(nsalts,NU(NY,NX),NY,NX))*VLWatMicPMNU &
      -AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  ENDDO

  DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_flx_diffus(nsalts)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2(nsalts,NU(NY,NX),NY,NX))*VLWatMicPMNU &
      -AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO nsalts=idsalt_pband_beg,idsalt_pband_end
    trcSalt_flx_diffus(nsalts)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2(nsalts,NU(NY,NX),NY,NX))*VLWatMicPMNU &
      -AZMAX1(trcSalt_solml2(nsalts,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO


!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFXS(nsalts,NU(NY,NX),NY,NX)=trcSalt_RFXS(nsalts,NU(NY,NX),NY,NX)+trcSalt_flx_diffus(nsalts)
  ENDDO

  end subroutine MacMicPoreSoluteDifusExchange

!------------------------------------------------------------------------------------------

  subroutine AccumHourlyMicMacPoreFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  INTEGER :: nsalts
!     begin_execution
!
!     X*FXS=hourly convective + diffusive solute flux between macropores and micropores
!     R*FXS=total convective + diffusive solute flux between macropores and micropores
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_XFXS(nsalts,NU(NY,NX),NY,NX)=trcSalt_XFXS(nsalts,NU(NY,NX),NY,NX)+trcSalt_RFXS(nsalts,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine AccumHourlyMicMacPoreFlux
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxBySurfaceOutflow(M,N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2
  real(r8) :: VFLW
  INTEGER :: nsalts
!     begin_execution

  IF(WatFlux4ErosionM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VLWatMicPM_vr(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,WatFlux4ErosionM_2DH(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO nsalts=idsalt_beg,idsalt_end
      trcSalt_RQR0(nsalts,N2,N1)=VFLW*AZMAX1(trcSalt_solml2(nsalts,0,N2,N1))
    ENDDO
  else
    trcSalt_RQR0(idsalt_beg:idsalt_end,N2,N1)=0.0_r8
  endif  
  end subroutine SoluteFluxBySurfaceOutflow


!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSurfNeighbors(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX,NHW,NHE,NVN,NVS
  integer :: N,NN,N4,N5,N4B,N5B,nsalts
  real(r8) :: FQRM,VFLW
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO N=1,2
    DO  NN=1,2
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
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
      IF(WatFlux4ErosionM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.1)THEN
          FQRM=QflxSurfRunoffM_2DH(M,N,2,N5,N4)/WatFlux4ErosionM_2DH(M,N2,N1)

          DO nsalts=idsalt_beg,idsalt_end
            trcSalt_RQR(nsalts,N,2,N5,N4)=trcSalt_RQR0(nsalts,N2,N1)*FQRM
          ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
          DO nsalts=idsalt_beg,idsalt_end
            trc_salt_rof_bounds(nsalts,N,2,N5,N4)=trc_salt_rof_bounds(nsalts,N,2,N5,N4)+trcSalt_RQR(nsalts,N,2,N5,N4)
          ENDDO
        ELSE
          trcSalt_RQR(idsalt_beg:idsalt_end,N,2,N5,N4)=0.0_r8
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            FQRM=QflxSurfRunoffM_2DH(M,N,1,N5B,N4B)/WatFlux4ErosionM_2DH(M,N2,N1)
            DO nsalts=idsalt_beg,idsalt_end
              trcSalt_RQR(nsalts,N,1,N5B,N4B)=trcSalt_RQR0(nsalts,N2,N1)*FQRM
            ENDDO
      !
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
            DO nsalts=idsalt_beg,idsalt_end
              trc_salt_rof_bounds(nsalts,N,1,N5B,N4B)=trc_salt_rof_bounds(nsalts,N,1,N5B,N4B)+trcSalt_RQR(nsalts,N,1,N5B,N4B)
            ENDDO
          ELSE
            trcSalt_RQR(idsalt_beg:idsalt_end,N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        trcSalt_RQR(idsalt_beg:idsalt_end,N,2,N5,N4)=0.0_r8

        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          trcSalt_RQR(idsalt_beg:idsalt_end,N,1,N5B,N4B)=0.0_r8
        ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     DrySnoFlxBySnoRedistM=snow transfer from watsub.f
!     RQS*=solute flux in snow transfer
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VOLS=volume of snowpack from watsub.f
!     *W2=solute content of snowpack
!
      IF(NN.EQ.1)THEN
!
!     IF NO SNOW DRIFT THEN NO TRANSPORT
!
        IF(ABS(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
          trcSalt_RQ(idsalt_beg:idsaltb_end,N,N5,N4)=0.0_r8
!
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
        ELSEIF(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4).GT.0.0)THEN
          IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(NY,NX))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)/VcumSnoDWI_col(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF

          DO nsalts=idsalt_beg,idsaltb_end
            trcSalt_RQ(nsalts,N,N5,N4)=VFLW*AZMAX1(trc_Saltml2_snvr(nsalts,1,N2,N1))
          ENDDO
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
        ELSEIF(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
          IF(VcumSnoDWI_col(N5,N4).GT.ZEROS2(N5,N4))THEN
            VFLW=AZMIN1(AMAX1(-VFLWX,DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)/VcumSnoDWI_col(N5,N4)))
          ELSE
            VFLW=-VFLWX
          ENDIF
          DO nsalts=idsalt_beg,idsaltb_end
            trcSalt_RQ(nsalts,N,N5,N4)=VFLW*AZMAX1(trc_Saltml2_snvr(nsalts,1,N5,N4))
          ENDDO
        ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQS*=hourly solute in snow transfer
!     RQS*=solute in snow transfer
!
        DO nsalts=idsalt_beg,idsalt_end
          trcSalt_XQS(nsalts,N,N5,N4)=trcSalt_XQS(nsalts,N,N5,N4)+trcSalt_RQ(nsalts,N,N5,N4)
        ENDDO
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSurfNeighbors
!------------------------------------------------------------------------------------------
   subroutine SoluteAdvMicropore(M,N,N1,N2,N3,N4,N5,N6)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8) :: VFLW
   integer :: nsalts
   real(r8) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     FLWM=water flux through soil micropore from watsub.f
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     RFL*S=solute diffusive flux through micropore
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *S2,*B2=micropore solute content in non-band,band
!
  IF(WaterFlow2MicPM_3D(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
    IF(VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,N6,N5,N4)/VLWatMicPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO nsalts=idsalt_beg,idsalt_KSO4
        trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N3,N2,N1))
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N3,N2,N1))*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N3,N2,N1))*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    ENDDO
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS TO CURRENT FROM
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
  ELSE
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MicPM_3D(M,N,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ENDIF

!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt3DFlo2CellM(nsalts,N,N6,N5,N4)=trcSalt_RFL(nsalts)
  ENDDO
  end subroutine SoluteAdvMicropore
!------------------------------------------------------------------------------------------
  subroutine SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  REAL(R8) :: VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
  integer :: nsalts
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  IF(THETW1(N3,N2,N1).GT.THETY_vr(N3,N2,N1) &
    .AND.THETW1(N6,N5,N4).GT.THETY_vr(N6,N5,N4) &
    .AND.VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!
!     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!     C*1,C*2=solute concentration in source,destination layer
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *S2,*B2=soil solute content in non-band,band
!
    VLWPA1=VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    VLWPB1=VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    VLWPA2=VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    VLWPB2=VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)

    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_conc_src(nsalts)=AZMAX1(trcSalt_solml2(nsalts,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VLWPA1.GT.ZEROS(N2,N1))THEN
      DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_src(nsalts)=AZMAX1(trcSalt_solml2(nsalts,N3,N2,N1)/VLWPA1)
      ENDDO
    ELSE
      trcSalt_conc_src(idsalt_psoil_beg:idsalt_psoil_end)=0.0_r8
    ENDIF

    IF(VLWPB1.GT.ZEROS(N2,N1))THEN
      DO nsalts=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_src(nsalts)=AZMAX1(trcSalt_solml2(nsalts,N3,N2,N1)/VLWPB1)
      ENDDO
    ELSE
      DO nsalts=0,idsalt_nuts
        trcSalt_conc_src(idsalt_H0PO4B+nsalts)=trcSalt_conc_src(idsalt_H0PO4+nsalts)
      ENDDO
    ENDIF

    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO

    IF(VLWPA2.GT.ZEROS(N5,N4))THEN
      DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      ENDDO
    ELSE
      trcSalt_conc_dest(idsalt_psoil_beg:idsalt_psoil_end)=0.0_r8
    ENDIF
    IF(VLWPB2.GT.ZEROS(N5,N4))THEN
      DO nsalts=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      ENDDO
    ELSE
      DO nsalts=0,idsalt_nuts
        trcSalt_conc_dest(idsalt_H0PO4B+nsalts)=trcSalt_conc_dest(idsalt_H0PO4+nsalts)
      ENDDO
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!
!     DLYR=soil layer thickness
!     TORT=micropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWM=water flux through soil micropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through micropore
!     *SGL2=solute diffusivity from hour1.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     XDPTH=cross-sectional area/distance between layers
!     C*1,C*2=micropore solute concentration in source,destination layer
!     DFV*=diffusive solute transfer through soil micropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!

    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TortMicPM_vr(M,N3,N2,N1)*DLYR1+TortMicPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)

    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFAL=(AquaIonDifusivty2_vr(idsalt_Al,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Al,idsalt_AlOH,idsalt_AlOH2,idsalt_AlOH3,idsalt_AlOH4,idsalt_AlSO4/))=DIFAL

    DIFFE=(AquaIonDifusivty2_vr(idsalt_Fe,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Fe,idsalt_FeOH,idsalt_FeOH2,idsalt_FeOH3,idsalt_FeOH4,idsalt_FeSO4/))=DIFFE

    trcSaltDiffConductance(idsalt_Hp)=(AquaIonDifusivty2_vr(idsalt_Hp,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFCA=(AquaIonDifusivty2_vr(idsalt_Ca,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Ca,idsalt_CaOH,idsalt_CO3,idsalt_HCO3,idsalt_CaSO4/))=DIFCA

    DIFMG=(AquaIonDifusivty2_vr(idsalt_Mg,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Mg,idsalt_MgOH2,idsalt_MgCO3,idsalt_MgHCO3,idsalt_SO4/))=DIFMG

    DIFNA=(AquaIonDifusivty2_vr(idsalt_Na,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Na,idsalt_NaCO3,idsalt_NaSO4/))=DIFNA

    DIFKA=(AquaIonDifusivty2_vr(idsalt_K,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_K,idsalt_KSO4/))=DIFKA

    trcSaltDiffConductance(idsalt_OH)=(AquaIonDifusivty2_vr(idsalt_OH,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_SO4)=(AquaIonDifusivty2_vr(idsalt_SO4,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_Cl)=(AquaIonDifusivty2_vr(idsalt_Cl,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_CO3)=(AquaIonDifusivty2_vr(idsalt_CO3,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_HCO3)=(AquaIonDifusivty2_vr(idsalt_HCO3,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!

    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_flx_diffus(nsalts)=trcSaltDiffConductance(nsalts)*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoiL_end
      trcSalt_flx_diffus(nsalts)=DIFPO*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_flx_diffus(nsalts)=DIFPO*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_flx_diffus(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt3DFlo2CellM(nsalts,N,N6,N5,N4)=trcSalt3DFlo2CellM(nsalts,N,N6,N5,N4)+trcSalt_flx_diffus(nsalts)
  ENDDO

  end subroutine SoluteDifsMicropore
!----------------------------------------------------------------------
  subroutine SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  integer :: nsalts
  real(r8) :: trcSalt_RFH(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW
!     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
!
  IF(WaterFlow2MacPM_3D(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     VOLWAH=macropore porosity
!     RFH*=solute diffusive flux through macropore
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *SH2,*BH2=macropore solute content in non-band,band
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    IF(VLWatMacPM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE IN VERTICAL FLUX
!
    IF(N.EQ.3.AND.VLMacP_vr(N6,N5,N4).GT.VLWatMacPM(M,N6,N5,N4))THEN
      DO nsalts=idsalt_beg,idsalt_KSO4
        trcSalt_RFH(nsalts)=VFLW*AZMAX1((trcSalt_soHml2(nsalts,N3,N2,N1) &
          -AZMIN1(trcSalt_RFXS(nsalts,NU(N2,N1),N2,N1))))
      ENDDO

      DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_RFH(nsalts)=VFLW*AZMAX1((trcSalt_soHml2(nsalts,N3,N2,N1) &
          -AZMIN1(trcSalt_RFXS(nsalts,NU(N2,N1),N2,N1))))*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
      ENDDO

      DO nsalts=idsalt_pband_beg,idsalt_pband_end
        trcSalt_RFH(nsalts)=VFLW*AZMAX1((trcSalt_soHml2(nsalts,N3,N2,N1) &
          -AZMIN1(trcSalt_RFXS(nsalts,NU(N2,N1),N2,N1))))*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
      ENDDO
!
!     OTHERWISE
!
    ELSE

      DO nsalts=idsalt_beg,idsalt_KSO4
        trcSalt_RFH(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N3,N2,N1))
      ENDDO

      DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_RFH(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N3,N2,N1)) &
          *trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
      ENDDO

      DO nsalts=idsalt_pband_beg,idsalt_pband_end
        trcSalt_RFH(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N3,N2,N1)) &
          *trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
      ENDDO
    ENDIF
  ELSEIF(WaterFlow2MacPM_3D(M,N,N6,N5,N4).LT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
    IF(VLWatMacPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_RFH(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))
    ENDDO
    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFH(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4)) &
        *trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFH(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4)) &
        *trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE FLUX
!
  ELSE
    trcSalt_RFH(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFHS(nsalts,N,N6,N5,N4)=trcSalt_RFH(nsalts)
  ENDDO

  end subroutine SoluteAdvMacropore

!----------------------------------------------------------------------
  subroutine SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_DFH(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)
  integer  :: nsalts

!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume

  IF(VLWatMacPM(M,N3,N2,N1).GT.THETY_vr(N3,N2,N1)*VLMacP_vr(N3,N2,N1) &
    .AND.VLWatMacPM(M,N6,N5,N4).GT.THETY_vr(N6,N5,N4)*VLMacP_vr(N6,N5,N4))THEN
!
!     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
!
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     *H2=macropore solute content
!     VLWatMacPM=macropore water content
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
    DO nsalts=idsalt_beg,idsaltb_end
      trcSalt_conc_src(nsalts)=AZMAX1(trcSalt_soHml2(nsalts,N3,N2,N1)/VLWatMacPM(M,N3,N2,N1))
      trcSalt_conc_dest(nsalts)=AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4))
    ENDDO
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     DLYR=soil layer thickness
!     TortMacPM=macropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through macropore
!     *SGL2=solute diffusivity from hour1.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     XDPTH=cross-sectional area/distance between layers
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     DFH*=diffusive solute transfer through soil macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TortMacPM(M,N3,N2,N1)*DLYR1+TortMacPM(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFAL=(AquaIonDifusivty2_vr(idsalt_Al,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Al,idsalt_AlOH,idsalt_AlOH2,idsalt_AlOH3,idsalt_AlOH4,idsalt_AlSO4/))=DIFAL

    DIFFE=(AquaIonDifusivty2_vr(idsalt_Fe,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Fe,idsalt_FeOH,idsalt_FeOH2,idsalt_FeOH3,idsalt_FeOH4,idsalt_FeSO4/))=DIFFE

    trcSaltDiffConductance(idsalt_Hp)=(AquaIonDifusivty2_vr(idsalt_Hp,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFCA=(AquaIonDifusivty2_vr(idsalt_Ca,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Ca,idsalt_CaOH,idsalt_CO3,idsalt_HCO3,idsalt_CaSO4/))=DIFCA

    DIFMG=(AquaIonDifusivty2_vr(idsalt_Mg,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Mg,idsalt_MgOH2,idsalt_MgCO3,idsalt_MgHCO3,idsalt_SO4/))=DIFMG

    DIFNA=(AquaIonDifusivty2_vr(idsalt_Na,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_Na,idsalt_NaCO3,idsalt_NaSO4/))=DIFNA

    DIFKA=(AquaIonDifusivty2_vr(idsalt_K,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance((/idsalt_K,idsalt_KSO4/))=DIFKA

    trcSaltDiffConductance(idsalt_OH)=(AquaIonDifusivty2_vr(idsalt_OH,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_SO4)=(AquaIonDifusivty2_vr(idsalt_SO4,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_Cl)=(AquaIonDifusivty2_vr(idsalt_Cl,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_CO3)=(AquaIonDifusivty2_vr(idsalt_CO3,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcSaltDiffConductance(idsalt_HCO3)=(AquaIonDifusivty2_vr(idsalt_HCO3,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!

    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_DFH(nsalts)=trcSaltDiffConductance(nsalts)*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_DFH(nsalts)=DIFPO*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_DFH(nsalts)=DIFPO*(trcSalt_conc_src(nsalts)-trcSalt_conc_dest(nsalts))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_DFH(idsalt_beg:idsaltb_end)=0._r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFHS(nsalts,N,N6,N5,N4)=trcSalt_RFHS(nsalts,N,N6,N5,N4)+trcSalt_DFH(nsalts)
  ENDDO
  end subroutine SoluteDifsMacropore
!----------------------------------------------------------------------
  subroutine SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
  integer :: nsalts
  real(r8) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FWatExMacP2MicPM=macro-micropore water transfer from watsub.f
!     VLWatMicPM,VLWatMacPM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FWatExMacP2MicPM(M,N6,N5,N4).GT.0.0)THEN
    IF(VLWatMacPM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM(M,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FWatExMacP2MicPM(M,N6,N5,N4).LT.0.0)THEN
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM(M,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFL(nsalts)=VFLW*AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    trcSalt_RFL(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFXS(nsalts,N6,N5,N4)=trcSalt_RFL(nsalts)
  ENDDO
  end subroutine SoluteAdvExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6
  real(r8) :: VOLWHS,VOLWT
  integer :: nsalts
  real(r8) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
!     VLWatMicPM,VLWatMacPM=micropore,macropore water-filled porosity from watsub.f
!     DFV*S,DFV*B=diffusive solute flux between macro- and micropore in non-band,band
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *2,*H2=solute content of micropores,macropores
!
  IF(VLWatMacPM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
    VOLWHS=AMIN1(XFRS*VGeomLayer_vr(N6,N5,N4),VLWatMacPM(M,N6,N5,N4))
    VOLWT=VLWatMicPM_vr(M,N6,N5,N4)+VOLWHS
    DO nsalts=idsalt_beg,idsalt_KSO4
      trcSalt_flx_diffus(nsalts)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*VOLWHS)/VOLWT
    ENDDO

    DO nsalts=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_flx_diffus(nsalts)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO nsalts=idsalt_pband_beg,idsalt_pband_end
      trcSalt_flx_diffus(nsalts)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2(nsalts,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2(nsalts,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_flx_diffus(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFXS(nsalts,N6,N5,N4)=trcSalt_RFXS(nsalts,N6,N5,N4)+trcSalt_flx_diffus(nsalts)
  ENDDO
  end subroutine SoluteDifsExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  integer :: nsalts
!     SOLUTE TRANSPORT IN MICROPORES

  call SoluteAdvMicropore(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     SOLUTE TRANSPORT IN MACROPORES
  call SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT
!     FLUXES = CONVECTIVE + DIFFUSIVE FLUXES
!
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux through micropores
!     X*FLW,X*FLB= hourly convective + diffusive solute flux through micropores in non-band,band
!     X*FHS=hourly convective + diffusive solute flux through macropores
!     X*FHW,X*FHB= hourly convective + diffusive solute flux through macropores in non-band,band
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,X*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,X*FHB=convective + diffusive solute flux through macropores in non-band,band
!

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt3DFlo2Cell(nsalts,N,N6,N5,N4)=trcSalt3DFlo2Cell(nsalts,N,N6,N5,N4) &
      +trcSalt3DFlo2CellM(nsalts,N,N6,N5,N4)
    trcSalt_XFHS(nsalts,N,N6,N5,N4)=trcSalt_XFHS(nsalts,N,N6,N5,N4) &
      +trcSalt_RFHS(nsalts,N,N6,N5,N4)
  ENDDO
  end subroutine SoluteAdvDifsMicMacpore
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  integer :: nsalts

  call SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)


!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS,R*FXB=total convective + diffusive solute flux between macro- and micropore in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     RFL*=convective flux between macro- and micropore
!     DFV*=diffusive solute flux between macro- and micropore
!

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FXS,X*FXB= hourly convective + diffusive solute flux between macro- and micropore in non-band,band
!     R*FXS,R*FXB=convective + diffusive solute flux between macro- and micropore in non-band,band
!
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_XFXS(nsalts,N6,N5,N4)=trcSalt_XFXS(nsalts,N6,N5,N4)+trcSalt_RFXS(nsalts,N6,N5,N4)
  ENDDO

  end subroutine SoluteAdvDifsExchMicMacpore
!----------------------------------------------------------------------
  subroutine UpdateSoluteInSubsurfNeighbors(M,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX,NHW,NHE,NVN,NVS
  integer :: IFLGB,L,N1,N2,N3,N4,N5,N6
  integer :: LL,N
  real(r8) :: VOLWHS,VOLWT
  real(r8) :: VFLW,VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: THETW1(JZ,JY,JX)

!     begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  IFLGB=0
  DO  L=1,NL(NY,NX)

    N1=NX;N2=NY;N3=L
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
    DO  N=FlowDirIndicator(N2,N1),3

      IF(N.EQ.1)THEN
        !west-east
        IF(NX.EQ.NHE)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N6=L
        ENDIF
      ELSEIF(N.EQ.2)THEN
        !north-south
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N6=L
        ENDIF
      ELSEIF(N.EQ.3)THEN
        !vertical
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4=NX
          N5=NY
          N6=L+1
        ENDIF
      ENDIF

      DO LL=N6,NL(NY,NX)
        IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          exit
        ENDIF
      enddo
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
!     WATER CONTENTS AND WATER FLUXES 'FLQM' FROM 'WATSUB'
!
!     VLSoilPoreMicP_vr,VLSoilMicP=soil volume excluding rock, macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     VLWatMicPM,VLWatMacPM=micropore,macropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!     FLPM=change in air volume
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
        IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4).AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
          THETW1(N3,N2,N1)=AZMAX1(VLWatMicPM_vr(M,N3,N2,N1)/VLSoilMicP_vr(N3,N2,N1))
          THETW1(N6,N5,N4)=AZMAX1(VLWatMicPM_vr(M,N6,N5,N4)/VLSoilMicP_vr(N6,N5,N4))
!
          call SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
          IF(N.EQ.3)THEN
             call SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
          ENDIF
        ELSEIF(N.NE.3)THEN
          THETW1(N3,N2,N1)=0.0_r8
          THETW1(N6,N5,N4)=0.0_r8
          trcSalt3DFlo2CellM(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8

          trcSalt_RFHS(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
        ENDIF
      ELSE
        THETW1(N3,N2,N1)=0.0_r8
        THETW1(N6,N5,N4)=0.0_r8

        trcSalt3DFlo2CellM(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
        trcSalt_RFHS(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSubsurfNeighbors
!------------------------------------------------------------------------------------------

  subroutine MacMicPoreFluxAdvPlusDifus(NY,NX,trcSalt_flx_diffus,trcSalt_RFL)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: trcSalt_flx_diffus(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_RFL(idsalt_beg:idsaltb_end)
  integer :: nsalts
!     begin_execution
!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_RFXS(nsalts,NU(NY,NX),NY,NX)=trcSalt_RFL(nsalts)+trcSalt_flx_diffus(nsalts)
  ENDDO
  end subroutine MacMicPoreFluxAdvPlusDifus
end module IngridTranspMod
