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
  public :: SaltHydroTranpModelM
  contains


!------------------------------------------------------------------------------------------

  subroutine SaltHydroTranpModelM(I,M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,M,NHW,NHE,NVN,NVS

  integer :: NY,NX,N1,N2
  real(r8) :: FLWRM1
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)
!     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     ADD SINKS FROM SOLUTE.F
      call SaltGeochemSinksM(NY,NX)
!
!     SOLUTE FLUXES FROM MELTING SNOWPACK TO
!     SOIL SURFACE FROM SNOWMELT IN 'WATSUB' AND
!     CONCENTRATIONS IN SNOWPACK
!
      call SnowSaltTranptM(M,NY,NX)
!
!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
      FLWRM1=WatFLo2LitrM(M,NY,NX)
!

      IF(FLWRM1.GT.0.0_r8)THEN
        call Litr2SoilAdvTransptM(M,NY,NX,FLWRM1,trcSalt_Adv2MicP_flx)
!
      ELSE
        call Soil2LitrAdvTransptM(M,NY,NX,FLWRM1,trcSalt_Adv2MicP_flx)
      ENDIF
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!
      IF((VGeomLayer_vr(0,NY,NX).GT.ZEROS(NY,NX) .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
        .AND. (VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN

        call Litr2SoilDifusTransptM(M,NY,NX,FLWRM1,trcSalt_flx_diffusM)
      ELSE
        trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES

      call LitrSoilAdvDifusTransptM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM)
!
      call AccumLitrSoilSaltTranspFluxM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM)
!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE IN SOIL
!     SURFACE LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!
!     MACROPORE TO MICROPORE TRANSFER
!
      IF(FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX).GT.0.0)THEN
        call MacPoreSaltAdvExM(M,NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
      ELSEIF(FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX).LT.0.0)THEN
        call MicPoreSaltAdvExM(M,NY,NX)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
      ELSE
        trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
      IF(VLWatMacPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
        call MacMicPoreSaltDifusExM(M,NY,NX,trcSalt_Difus_Mac2Micp_flxM)
      ELSE
        trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
!
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
      call AccumMicMacPoreSaltTranspM(NY,NX)
!
!     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
!     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE LAYER
!
!
      N1=NX;N2=NY

      call SaltFlowThruLitrRofM(M,N1,N2)
!
      call SaltFlowThruSnowRofM(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
      call SaltExchXGridsM(M,NY,NX,NHW,NHE,NVN,NVS)

    ENDDO
  ENDDO
  end subroutine SaltHydroTranpModelM

!------------------------------------------------------------------------------------------

  subroutine SaltGeochemSinksM(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,idsalt
!     begin_execution
!

  DO L=NU(NY,NX),NL(NY,NX)
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_solml2_vr(idsalt,L,NY,NX)=trcSalt_solml2_vr(idsalt,L,NY,NX)-trcSalt_RGeoChem_flxM_vr(idsalt,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine SaltGeochemSinksM
!------------------------------------------------------------------------------------------

  subroutine SnowSaltTranptM(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX

  integer :: L,L2,idsalt,ids
  logical :: ICHKL
  real(r8) :: VFLWW   !fraction flowed into next snow layer
  real(r8) :: VFLWR,VFLWS,VFLWPO4,VFLWPOB
  
!     begin_execution
!
!     VLSnowHeatCapM,VLHeatCapSnowMin_col=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     WatFlowInSnowM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack

  ICHKL=.false.
  DO L=1,JS
    !
    IF(VLSnowHeatCapM_snvr(M,L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      L2=MIN(JS,L+1)

      !Not bottom layer, next layer has significant snow
      IF(L.LT.JS .AND. VLSnowHeatCapM_snvr(M,L2,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
        IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0_r8,WatFlowInSnowM_snvr(M,L2,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
        ELSE
          VFLWW=1.0_r8
        ENDIF

        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_AquaAdv_flxM_snvr(idsalt,L2,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWW
          trcSalt_AquaAdv_flx_snvr(idsalt,L2,NY,NX)  = trcSalt_AquaAdv_flx_snvr(idsalt,L2,NY,NX)+trcSalt_AquaAdv_flxM_snvr(idsalt,L2,NY,NX)
        ENDDO

        !bottom snow layer or L2 has insignificant snow layer
      ELSE

        !Insignificant snow layer
        IF(L.LT.JS)THEN
          trcSalt_AquaAdv_flxM_snvr(idsalt_beg:idsalt_end,L2,NY,NX)=0.0_r8
        ENDIF
!
        ! SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
        ! WatFlowSno2LitRM_col,WatFlowSno2MicPM_col,WatFlowSno2MacPM_col=total water flux to litter,soil micropore,macropore
        ! CVRD,BARE=litter cover fraction,1-CVRD
        ! VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
        ! VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
        IF(.not.ICHKL)THEN
          !flow into soil
          IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,WatFlowSno2LitRM_col(M,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX))/VLWatSnow_snvr(L,NY,NX)))
          ELSE
            VFLWR=FracSurfByLitR_col(NY,NX)
            VFLWS=FracSurfBareSoil_col(NY,NX)
          ENDIF
          VFLWPO4=VFLWS*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWR
            trcSalt_AquaADV_Snow2Litr_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Litr_flx(idsalt,NY,NX)+trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX)
          ENDDO

          DO idsalt=idsalt_beg,idsalt_KSO4
            trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWS
            trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)+trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)
          ENDDO

          DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
            ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
            !soil
            trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWPO4
            trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)+trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)
            !band
            trcSalt_AquaADV_Snow2Soil_flxM(ids,NY,NX) = trcSalt_ml2_snvr(ids,L,NY,NX)*VFLWPOB
            trcSalt_AquaADV_Snow2Soil_flx(ids,NY,NX)  = trcSalt_AquaADV_Snow2Soil_flx(ids,NY,NX)+trcSalt_AquaADV_Snow2Soil_flxM(ids,NY,NX)
          ENDDO
          ICHKL=.true.
        ENDIF
      ENDIF
    ELSE
      trcSalt_AquaADV_Snow2Litr_flxM(idsalt_beg:idsalt_end,NY,NX)  = 0.0_r8
      trcSalt_AquaADV_Snow2Soil_flxM(idsalt_beg:idsaltb_end,NY,NX) = 0.0_r8
    ENDIF
  ENDDO
  end subroutine SnowSaltTranptM
!------------------------------------------------------------------------------------------

  subroutine Litr2SoilAdvTransptM(M,NY,NX,FLWRM1,trcSalt_Adv2MicP_flx)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW
  integer :: idsalt,ids

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!
  IF(VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VLWatMicPM_vr(M,0,NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF
  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))
  ENDDO

  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
    trcSalt_Adv2MicP_flx(idsalt) = VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcSalt_Adv2MicP_flx(ids)    = VFLW*AZMAX1(trcSalt_solml2_vr(ids,0,NY,NX))*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine Litr2SoilAdvTransptM
!------------------------------------------------------------------------------------------

  subroutine Soil2LitrAdvTransptM(M,NY,NX,FLWRM1,trcSalt_Adv2MicP_flx)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  integer :: idsalt
  real(r8) :: VFLW

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE

  IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX))
  ENDDO
  end subroutine Soil2LitrAdvTransptM
!------------------------------------------------------------------------------------------

  subroutine Litr2SoilDifusTransptM(M,NY,NX,FLWRM1,trcSalt_flx_diffusM)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)    !tracer solute concentration in source cell
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)   !tracer solutte concentration in destination cell
  real(r8) :: VOLWPB,VOLWPA,TORT0,TORT1
  real(r8) :: DLYR0
  integer :: idsalt
!     begin_execution
!
!     MICROPORE CONCENTRATIONS FROM WATER IN RESIDUE AND SOIL SURFACE
!
!
  VOLWPA=VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  VOLWPB=VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
  ENDDO

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX))
  ENDDO

  IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX)/VOLWPA)
    ENDDO
  ELSE
    trcSalt_conc_dest(idsalt_psoiL_beg:idsalt_psoil_end)=0.0_r8
  ENDIF
  
  IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX)/VOLWPB)
    ENDDO
  ELSE
    DO idsalt=0,idsalt_nuts
      trcSalt_conc_dest(idsalt_H0PO4B+idsalt)=trcSalt_conc_dest(idsalt_H0PO4+idsalt)
    ENDDO
  ENDIF
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
!
!
  DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX))
  TORT0 = TortMicPM_vr(M,0,NY,NX)*FracSurfByLitR_col(NY,NX)
  DLYR1 = AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))
  TORT1 = TortMicPM_vr(M,NU(NY,NX),NY,NX)
  TORTL = AMIN1(1.0_r8,(TORT0+TORT1)/(DLYR0+DLYR1))
  DISPN = DISP_3D(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))

  DIFPO=(SoluteDifusivitytscaledM_vr(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)

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
!

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_flx_diffusM(idsalt)=trcSaltDiffConductance(idsalt)*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))
  ENDDO

  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_flx_diffusM(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO idsalt=0,idsalt_nuts
    trcSalt_flx_diffusM(idsalt_H0PO4B+idsalt)=DIFPO*(trcSalt_conc_src(idsalt_H0PO4+idsalt)&
      -trcSalt_conc_dest(idsalt_H0PO4B+idsalt))*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine Litr2SoilDifusTransptM
!------------------------------------------------------------------------------------------

  subroutine LitrSoilAdvDifusTransptM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  integer :: idsalt
!     begin_execution
!
!
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)=trcSalt_Precip2LitrM(idsalt,NY,NX)+trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX) &
      -trcSalt_Adv2MicP_flx(idsalt)-trcSalt_flx_diffusM(idsalt)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MicpTranspFlxM_3D(idsalt,3,NU(NY,NX),NY,NX)=trcSalt_Precip2MicpM(idsalt,NY,NX) &
      +trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)+trcSalt_Adv2MicP_flx(idsalt)+trcSalt_flx_diffusM(idsalt)
  ENDDO

  end subroutine LitrSoilAdvDifusTransptM
!------------------------------------------------------------------------------------------

  subroutine AccumLitrSoilSaltTranspFluxM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8), intent(in) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)

  integer :: idsalt
!     begin_execution
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_TransptMicP_3D(idsalt,3,0,NY,NX)=trcSalt_TransptMicP_3D(idsalt,3,0,NY,NX)+trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX) &
      -trcSalt_Adv2MicP_flx(idsalt)-trcSalt_flx_diffusM(idsalt)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_TransptMicP_3D(idsalt,3,NU(NY,NX),NY,NX)=trcSalt_TransptMicP_3D(idsalt,3,NU(NY,NX),NY,NX) &
      +trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)+trcSalt_Adv2MicP_flx(idsalt)+trcSalt_flx_diffusM(idsalt)
  ENDDO
  end subroutine AccumLitrSoilSaltTranspFluxM
!------------------------------------------------------------------------------------------

  subroutine MacPoreSaltAdvExM(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcSalt_Adv2MacP_flx(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW
  integer :: idsalt
!     begin_execution

  IF(VLWatMacPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX)/VLWatMacPM_vr(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Adv2MacP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,NU(NY,NX),NY,NX))
  ENDDO

  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_Adv2MacP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_pband_beg,idsalt_pband_end
    trcSalt_Adv2MacP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(NY,NX),NY,NX)=trcSalt_Adv2MacP_flx(idsalt)
  ENDDO
  end subroutine MacPoreSaltAdvExM
!------------------------------------------------------------------------------------------

  subroutine MicPoreSaltAdvExM(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  integer :: idsalt
  real(r8) :: VFLW
!     begin_execution
  IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Adv2MicP_flx(idsalt_Al)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt_Al,NU(NY,NX),NY,NX))
  ENDDO

  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_pband_beg,idsalt_pband_end
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(NY,NX),NY,NX)=trcSalt_Adv2MicP_flx(idsalt)
  ENDDO
  end subroutine MicPoreSaltAdvExM

!------------------------------------------------------------------------------------------

  subroutine MacMicPoreSaltDifusExM(M,NY,NX,trcSalt_Difus_Mac2Micp_flxM)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)
  real(r8) :: VOLWHS,VOLWT,VLWatMicPMNU
  integer :: idsalt
!     begin_execution
!

  VLWatMicPMNU = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)
  VOLWHS       = AMIN1(XFRS*VGeomLayer_vr(NU(NY,NX),NY,NX),VLWatMacPM_vr(M,NU(NY,NX),NY,NX))
  VOLWT        = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)+VOLWHS

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,NU(NY,NX),NY,NX))*VLWatMicPMNU &
      -AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  ENDDO

  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,NU(NY,NX),NY,NX))*VLWatMicPMNU &
      -AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_pband_beg,idsalt_pband_end
    trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,NU(NY,NX),NY,NX))*VLWatMicPMNU &
      -AZMAX1(trcSalt_solml2_vr(idsalt,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

!
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(NY,NX),NY,NX)=trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(NY,NX),NY,NX)+trcSalt_Difus_Mac2Micp_flxM(idsalt)
  ENDDO

  end subroutine MacMicPoreSaltDifusExM

!------------------------------------------------------------------------------------------

  subroutine AccumMicMacPoreSaltTranspM(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  INTEGER :: idsalt
!     begin_execution
!
!
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flx_vr(idsalt,NU(NY,NX),NY,NX)=trcSalt_Mac2MicPore_flx_vr(idsalt,NU(NY,NX),NY,NX)+trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine AccumMicMacPoreSaltTranspM
!------------------------------------------------------------------------------------------

  subroutine SaltFlowThruLitrRofM(M,N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2
  real(r8) :: VFLW
  INTEGER :: idsalt
!     begin_execution

  IF(SurfRunoffWatFluxM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VLWatMicPM_vr(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,SurfRunoffWatFluxM_2DH(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_FloXSurRof_flxM(idsalt,N2,N1)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,N2,N1))
    ENDDO
  else
    trcSalt_FloXSurRof_flxM(idsalt_beg:idsalt_end,N2,N1)=0.0_r8
  endif  
  end subroutine SaltFlowThruLitrRofM


!------------------------------------------------------------------------------------------

  subroutine SaltFlowThruSnowRofM(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX,NHW,NHE,NVN,NVS
  integer :: N,NN,N4,N5,N4B,N5B,idsalt
  real(r8) :: FQRM,VFLW
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO N=1,2
    DO  NN=1,2
      IF(N.EQ.iEastWestDirection)THEN
        IF(NX.EQ.NHE .AND. NN.EQ.iOutflow .OR. NX.EQ.NHW .AND. NN.EQ.iInflow)THEN
          cycle
        ELSE
          N4  = NX+1  !eastward
          N5  = NY
          N4B = NX-1  !westward
          N5B = NY
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        IF(NY.EQ.NVS .AND. NN.EQ.iOutflow .OR. NY.EQ.NVN .AND. NN.EQ.iInflow)THEN
          cycle
        ELSE
          N4  = NX
          N5  = NY+1   !southward
          N4B = NX
          N5B = NY-1   !northward
        ENDIF
      ENDIF
!
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
      IF(SurfRunoffWatFluxM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.iOutflow)THEN
          FQRM=QflxSurfRunoffM_2DH(M,N,2,N5,N4)/SurfRunoffWatFluxM_2DH(M,N2,N1)

          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_FloXSurRof_flxM_2DH(idsalt,N,2,N5,N4)=trcSalt_FloXSurRof_flxM(idsalt,N2,N1)*FQRM
          ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_FloXSurRunoff_2D(idsalt,N,2,N5,N4)=trcSalt_FloXSurRunoff_2D(idsalt,N,2,N5,N4)+trcSalt_FloXSurRof_flxM_2DH(idsalt,N,2,N5,N4)
          ENDDO
        ELSE
          trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,2,N5,N4)=0.0_r8
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.iInflow)THEN
          !legitimate inner grid
          IF(N4B.GT.0 .AND. N5B.GT.0)THEN
            FQRM=QflxSurfRunoffM_2DH(M,N,1,N5B,N4B)/SurfRunoffWatFluxM_2DH(M,N2,N1)
            DO idsalt=idsalt_beg,idsalt_end
              trcSalt_FloXSurRof_flxM_2DH(idsalt,N,1,N5B,N4B)=trcSalt_FloXSurRof_flxM(idsalt,N2,N1)*FQRM
            ENDDO
      !
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
            DO idsalt=idsalt_beg,idsalt_end
              trcSalt_FloXSurRunoff_2D(idsalt,N,1,N5B,N4B)=trcSalt_FloXSurRunoff_2D(idsalt,N,1,N5B,N4B)+trcSalt_FloXSurRof_flxM_2DH(idsalt,N,1,N5B,N4B)
            ENDDO
          ELSE
            trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,2,N5,N4)=0.0_r8

        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,1,N5B,N4B)=0.0_r8
        ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     DrySnoFlxBySnoRedistM=snow transfer from watsub.f

!     VOLS=volume of snowpack from watsub.f
!     *W2=solute content of snowpack
!
      IF(NN.EQ.iOutflow)THEN
!
!     IF NO SNOW DRIFT THEN NO TRANSPORT
!
        IF(ABS(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
          trcSalt_SnowDrift_flxM_2D(idsalt_beg:idsaltb_end,N,N5,N4)=0.0_r8
!
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
        ELSEIF(DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4).GT.0.0)THEN
          IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(NY,NX))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxBySnoRedistM_2DH(M,N,N5,N4)/VcumSnoDWI_col(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF

          DO idsalt=idsalt_beg,idsaltb_end
            trcSalt_SnowDrift_flxM_2D(idsalt,N,N5,N4)=VFLW*AZMAX1(trcSalt_ml2_snvr(idsalt,1,N2,N1))
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
          DO idsalt=idsalt_beg,idsaltb_end
            trcSalt_SnowDrift_flxM_2D(idsalt,N,N5,N4)=VFLW*AZMAX1(trcSalt_ml2_snvr(idsalt,1,N5,N4))
          ENDDO
        ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_FloXSnow_2DH(idsalt,N,N5,N4)=trcSalt_FloXSnow_2DH(idsalt,N,N5,N4)+trcSalt_SnowDrift_flxM_2D(idsalt,N,N5,N4)
        ENDDO
      ENDIF
    enddo
  enddo
  end subroutine SaltFlowThruSnowRofM
!------------------------------------------------------------------------------------------
   subroutine SoluteAdvMicroporeM(M,N,N1,N2,N3,N4,N5,N6)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8) :: VFLW
   integer :: idsalt
   real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
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

    DO idsalt=idsalt_beg,idsalt_KSO4
        trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1))*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1))*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
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
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ENDIF

!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MicpTranspFlxM_3D(idsalt,N,N6,N5,N4)=trcSalt_Adv2MicP_flx(idsalt)
  ENDDO
  end subroutine SoluteAdvMicroporeM
!------------------------------------------------------------------------------------------
  subroutine SoluteDifsMicroporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N3
  integer, intent(in) :: N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  REAL(R8) :: VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  integer :: idsalt
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  IF(THETW1(N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1) &
    .AND.THETW1(N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4) &
    .AND.VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!
!     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!
    VLWPA1=VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    VLWPB1=VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    VLWPA2=VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    VLWPB2=VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VLWPA1.GT.ZEROS(N2,N1))THEN
      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1)/VLWPA1)
      ENDDO
    ELSE
      trcSalt_conc_src(idsalt_psoil_beg:idsalt_psoil_end)=0.0_r8
    ENDIF

    IF(VLWPB1.GT.ZEROS(N2,N1))THEN
      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1)/VLWPB1)
      ENDDO
    ELSE
      DO idsalt=0,idsalt_nuts
        trcSalt_conc_src(idsalt_H0PO4B+idsalt)=trcSalt_conc_src(idsalt_H0PO4+idsalt)
      ENDDO
    ENDIF

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO

    IF(VLWPA2.GT.ZEROS(N5,N4))THEN
      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      ENDDO
    ELSE
      trcSalt_conc_dest(idsalt_psoil_beg:idsalt_psoil_end)=0.0_r8
    ENDIF
    IF(VLWPB2.GT.ZEROS(N5,N4))THEN
      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      ENDDO
    ELSE
      DO idsalt=0,idsalt_nuts
        trcSalt_conc_dest(idsalt_H0PO4B+idsalt)=trcSalt_conc_dest(idsalt_H0PO4+idsalt)
      ENDDO
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!

    DLYR1=AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL=(TortMicPM_vr(M,N3,N2,N1)*DLYR1+TortMicPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)

    DISPN=DISP_3D(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(SoluteDifusivitytscaledM_vr(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

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

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_flx_diffusM(idsalt)=trcSaltDiffConductance(idsalt)*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoiL_end
      trcSalt_flx_diffusM(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_flx_diffusM(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MicpTranspFlxM_3D(idsalt,N,N6,N5,N4)=trcSalt_MicpTranspFlxM_3D(idsalt,N,N6,N5,N4)+trcSalt_flx_diffusM(idsalt)
  ENDDO

  end subroutine SoluteDifsMicroporeM
!----------------------------------------------------------------------
  subroutine SoluteAdvMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  integer :: idsalt
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
    IF(VLWatMacPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE IN VERTICAL FLUX
!
    IF(N.EQ.3.AND.VLMacP_vr(N6,N5,N4).GT.VLWatMacPM_vr(M,N6,N5,N4))THEN
      DO idsalt=idsalt_beg,idsalt_KSO4
        trcSalt_RFH(idsalt)=VFLW*AZMAX1((trcSalt_soHml2_vr(idsalt,N3,N2,N1) &
          -AZMIN1(trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(N2,N1),N2,N1))))
      ENDDO

      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1((trcSalt_soHml2_vr(idsalt,N3,N2,N1) &
          -AZMIN1(trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(N2,N1),N2,N1))))*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
      ENDDO

      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1((trcSalt_soHml2_vr(idsalt,N3,N2,N1) &
          -AZMIN1(trcSalt_Mac2MicPore_flxM_vr(idsalt,NU(N2,N1),N2,N1))))*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
      ENDDO
!
!     OTHERWISE
!
    ELSE

      DO idsalt=idsalt_beg,idsalt_KSO4
        trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1))
      ENDDO

      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1)) &
          *trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
      ENDDO

      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1)) &
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
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4)) &
        *trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4)) &
        *trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE FLUX
!
  ELSE
    trcSalt_RFH(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
!

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MacpTranspFlxM_3D(idsalt,N,N6,N5,N4)=trcSalt_RFH(idsalt)
  ENDDO

  end subroutine SoluteAdvMacroporeM

!----------------------------------------------------------------------
  subroutine SoluteDifsMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_MacPore_Diffus(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)
  integer  :: idsalt

!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume

  IF(VLWatMacPM_vr(M,N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1)*VLMacP_vr(N3,N2,N1) &
    .AND.VLWatMacPM_vr(M,N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4)*VLMacP_vr(N6,N5,N4))THEN
!
!     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
!
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
    ENDDO
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     XDPTH=cross-sectional area/distance between layers
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     DFH*=diffusive solute transfer through soil macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMacPM_vr(M,N3,N2,N1)*DLYR1+TortMacPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN = DISP_3D(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO = (SoluteDifusivitytscaledM_vr(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

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

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_MacPore_Diffus(idsalt)=trcSaltDiffConductance(idsalt)*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_MacPore_Diffus(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_MacPore_Diffus(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_MacPore_Diffus(idsalt_beg:idsaltb_end)=0._r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MacpTranspFlxM_3D(idsalt,N,N6,N5,N4)=trcSalt_MacpTranspFlxM_3D(idsalt,N,N6,N5,N4)+trcSalt_MacPore_Diffus(idsalt)
  ENDDO
  end subroutine SoluteDifsMacroporeM
!----------------------------------------------------------------------
  subroutine SaltAdvXMicMacporesM(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
  integer :: idsalt
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FWatExMacP2MicPM_vr(M,N6,N5,N4).GT.0.0)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FWatExMacP2MicPM_vr(M,N6,N5,N4).LT.0.0)THEN
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)=trcSalt_Adv2MicP_flx(idsalt)
  ENDDO
  end subroutine SaltAdvXMicMacporesM
!----------------------------------------------------------------------
  subroutine SoluteDifsExchMicMacporeM(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6
  real(r8) :: VOLWHS,VOLWT
  integer :: idsalt
  real(r8) :: trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
  IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
    VOLWHS=AMIN1(XFRS*VGeomLayer_vr(N6,N5,N4),VLWatMacPM_vr(M,N6,N5,N4))
    VOLWT=VLWatMicPM_vr(M,N6,N5,N4)+VOLWHS
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*VOLWHS)/VOLWT
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)=trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)+trcSalt_Difus_Mac2Micp_flxM(idsalt)
  ENDDO
  end subroutine SoluteDifsExchMicMacporeM
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsMicMacporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  integer :: idsalt
!     SOLUTE TRANSPORT IN MICROPORES

  call SoluteAdvMicroporeM(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMicroporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     SOLUTE TRANSPORT IN MACROPORES
  call SoluteAdvMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
  !   TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT
  !   FLUXES = CONVECTIVE + DIFFUSIVE FLUXES
  !   ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
  
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_TransptMicP_3D(idsalt,N,N6,N5,N4)=trcSalt_TransptMicP_3D(idsalt,N,N6,N5,N4) &
      +trcSalt_MicpTranspFlxM_3D(idsalt,N,N6,N5,N4)
    trcSalt_TransptMacP_3D(idsalt,N,N6,N5,N4)=trcSalt_TransptMacP_3D(idsalt,N,N6,N5,N4) &
      +trcSalt_MacpTranspFlxM_3D(idsalt,N,N6,N5,N4)
  ENDDO
  end subroutine SoluteAdvDifsMicMacporeM
!----------------------------------------------------------------------
  subroutine SaltAdvDifuXMicMacporesM(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  integer :: idsalt

  call SaltAdvXMicMacporesM(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsExchMicMacporeM(M,N,NY,NX,N1,N2,N3,N4,N5,N6)

!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
  !
  !     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
  !
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flx_vr(idsalt,N6,N5,N4)=trcSalt_Mac2MicPore_flx_vr(idsalt,N6,N5,N4)+trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)
  ENDDO

  end subroutine SaltAdvDifuXMicMacporesM
!----------------------------------------------------------------------
  subroutine SaltExchXGridsM(M,NY,NX,NHW,NHE,NVN,NVS)
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

      IF(N.EQ.iEastWestDirection)THEN
        !west-east
        IF(NX.EQ.NHE)THEN
          cycle
        ELSE
          N4 = NX+1
          N5 = NY
          N6 = L
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !north-south
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4 = NX
          N5 = NY+1
          N6 = L
        ENDIF
      ELSEIF(N.EQ.iVerticalDirection)THEN
        !vertical
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4 = NX
          N5 = NY
          N6 = L+1
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
!
      IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
        IF(N3.GE.NUM(N2,N1) .AND. N6.GE.NUM(N5,N4) .AND. N3.LE.NL(N2,N1) .AND. N6.LE.NL(N5,N4))THEN
          THETW1(N3,N2,N1)=AZMAX1(VLWatMicPM_vr(M,N3,N2,N1)/VLSoilMicP_vr(N3,N2,N1))
          THETW1(N6,N5,N4)=AZMAX1(VLWatMicPM_vr(M,N6,N5,N4)/VLSoilMicP_vr(N6,N5,N4))
!
          call SoluteAdvDifsMicMacporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)
!
          ! MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
          ! LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
          ! FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
          IF(N.EQ.iVerticalDirection)THEN
             call SaltAdvDifuXMicMacporesM(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
          ENDIF
        ELSEIF(N.NE.iVerticalDirection)THEN
          THETW1(N3,N2,N1)=0.0_r8
          THETW1(N6,N5,N4)=0.0_r8
          trcSalt_MicpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
          trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
        ENDIF
      ELSE
        THETW1(N3,N2,N1) = 0.0_r8
        THETW1(N6,N5,N4) = 0.0_r8

        trcSalt_MicpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
        trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,N6,N5,N4)=0.0_r8
      ENDIF
    enddo
  enddo
  end subroutine SaltExchXGridsM

end module IngridTranspMod
