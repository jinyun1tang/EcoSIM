module GeochemAPI

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoluteChemDataType, only : solute_flx_type, chem_var_type
  use SoluteMod
  use AqueChemDatatype
  use SoilPropertyDataType
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use SOMDataType
  use EcoSIMSolverPar
  use SoilWaterDataType
  use GridDataType
  implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: soluteModel
  contains


  subroutine soluteModel(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOLUTE TRANSFORMATIONS
!     FROM THERMODYNAMIC EQUILIBRIA
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW, NHE, NVN, NVS

  ! declaration of local variables
  integer :: L,NY,NX,NPI
  real(r8) :: RHP1,RHP2,RN3S,RN4S
  type(solute_flx_type) :: solflx
  type(chem_var_type) :: chemvar

!     begin_execution
  NPI=INT(NPH/2)
  DO   NX=NHW,NHE
    DO   NY=NVN,NVS
      DO   L=NU(NY,NX),NL(NY,NX)
        IF(VLSoilPoreMicP(L,NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
!
!     WATER VOLUME IN NON-BAND AND BAND SOIL ZONES
!
!     VLWatMicPM=soil water volume
!     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
!     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
!     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
!     SoilMicPMassLayer=soil mass
!
          chemvar%VLWatMicPNH=VLWatMicPM(NPH,L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)
          chemvar%VLWatMicPNB=VLWatMicPM(NPH,L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)
          chemvar%VLWatMicPNO=VLWatMicPM(NPH,L,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)
          chemvar%VLWatMicPNZ=VLWatMicPM(NPH,L,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)
          chemvar%VLWatMicPPO=VLWatMicPM(NPH,L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
          chemvar%VLWatMicPPB=VLWatMicPM(NPH,L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
          IF(SoilMicPMassLayer(L,NY,NX).GT.ZEROS(NY,NX))THEN
            chemvar%SoilMicPMassLayerX=SoilMicPMassLayer(L,NY,NX)
            chemvar%BKVLNH=SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)
            chemvar%BKVLNB=SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)
            chemvar%BKVLNO=SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)
            chemvar%BKVLNZ=SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)
            chemvar%BKVLPO=SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
            chemvar%BKVLPB=SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
          ELSE
            chemvar%SoilMicPMassLayerX=VLMicP(L,NY,NX)
            chemvar%BKVLNH=chemvar%VLWatMicPNH
            chemvar%BKVLNB=chemvar%VLWatMicPNB
            chemvar%BKVLNO=chemvar%VLWatMicPNO
            chemvar%BKVLNZ=chemvar%VLWatMicPNZ
            chemvar%BKVLPO=chemvar%VLWatMicPPO
            chemvar%BKVLPB=chemvar%VLWatMicPPB
          ENDIF

          call UpdateSoilFertlizer(L,NY,NX,chemvar)

          call GeochemAPISend(L,NY,NX,chemvar,solflx)

          call GeoChemEquilibria(chemvar, solflx)

          call GeochemAPIRecv(L,NY,NX,solflx)

          call UpdateFertilizerBand(L,NY,NX)
        ENDIF
      ENDDO
!
!     SURFACE RESIDUE
!
      call UpdateSoluteInSurfaceResidue(NX,NY)
!
    ENDDO
  ENDDO
  end subroutine soluteModel

!------------------------------------------------------------------------

  subroutine GeochemAPISend(L,NY,NX,chemvar,solflx)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(solute_flx_type), intent(inout) :: solflx
  type(chem_var_type), intent(inout) :: chemvar

  chemvar%XCEC=trcx_solml(idx_CEC,L,NY,NX)
  chemvar%PH=PH(L,NY,NX)
  chemvar%CAL=CAL(L,NY,NX)
  chemvar%CFE=CFE(L,NY,NX)
  chemvar%VLWatMicPM=VLWatMicPM(NPH,L,NY,NX)
  chemvar%ZMG=trcSalt_solml(idsalt_Mg,L,NY,NX)
  chemvar%ZNA=trcSalt_solml(idsalt_Na,L,NY,NX)
  chemvar%ZKA=trcSalt_solml(idsalt_K,L,NY,NX)
  chemvar%CCA=CCA(L,NY,NX)
  chemvar%SoilMicPMassLayer=SoilMicPMassLayer(L,NY,NX)
  chemvar%XAEC=trcx_solml(idx_AEC,L,NY,NX)
  chemvar%VLNH4=trcs_VLN(ids_NH4,L,NY,NX)
  chemvar%GKC4=GKC4(L,NY,NX)
  chemvar%VLNHB=trcs_VLN(ids_NH4B,L,NY,NX)
  chemvar%GKCA=GKCA(L,NY,NX)
  chemvar%GKCH=GKCH(L,NY,NX)
  chemvar%GKCM=GKCM(L,NY,NX)
  chemvar%GKCK=GKCK(L,NY,NX)
  chemvar%GKCN=GKCN(L,NY,NX)
  chemvar%ZNO3S=trc_solml(ids_NO3,L,NY,NX)
  chemvar%ZNO3B=trc_solml(ids_NO3B,L,NY,NX)
  chemvar%XZHYS=XZHYS(L,NY,NX)
  chemvar%ZHY=trcSalt_solml(idsalt_Hp,L,NY,NX)
  chemvar%ZOH=trcSalt_solml(idsalt_OH,L,NY,NX)
  chemvar%ZAL=trcSalt_solml(idsalt_Al,L,NY,NX)
  chemvar%ZFE=trcSalt_solml(idsalt_Fe,L,NY,NX)
  chemvar%ZCA=trcSalt_solml(idsalt_Ca,L,NY,NX)
  chemvar%ZSO4=trcSalt_solml(idsalt_SO4,L,NY,NX)
  chemvar%ZCL=trcSalt_solml(idsalt_Cl,L,NY,NX)
  chemvar%ZCO3=trcSalt_solml(idsalt_CO3,L,NY,NX)
  chemvar%ZHCO3=trcSalt_solml(idsalt_HCO3,L,NY,NX)
  chemvar%CO2S=trc_solml(idg_CO2,L,NY,NX)
  chemvar%ZALOH1=trcSalt_solml(idsalt_AlOH,L,NY,NX)
  chemvar%ZALOH2=trcSalt_solml(idsalt_AlOH2,L,NY,NX)
  chemvar%ZALOH3=trcSalt_solml(idsalt_AlOH3,L,NY,NX)
  chemvar%ZALOH4=trcSalt_solml(idsalt_AlOH4,L,NY,NX)
  chemvar%ZALS=trcSalt_solml(idsalt_AlSO4,L,NY,NX)
  chemvar%ZFEOH1=trcSalt_solml(idsalt_FeOH,L,NY,NX)
  chemvar%ZFEOH2=trcSalt_solml(idsalt_FeOH2,L,NY,NX)
  chemvar%ZFEOH3=trcSalt_solml(idsalt_FeOH3,L,NY,NX)
  chemvar%ZFEOH4=trcSalt_solml(idsalt_FeOH4,L,NY,NX)
  chemvar%ZFES=trcSalt_solml(idsalt_FeSO4,L,NY,NX)
  chemvar%ZCAO=trcSalt_solml(idsalt_CaOH2,L,NY,NX)
  chemvar%ZCAC=trcSalt_solml(idsalt_CaCO3,L,NY,NX)
  chemvar%ZCAH=trcSalt_solml(idsalt_CaHCO3,L,NY,NX)
  chemvar%ZCAS=trcSalt_solml(idsalt_CaSO4,L,NY,NX)
  chemvar%ZMGO=trcSalt_solml(idsalt_MgOH2,L,NY,NX)
  chemvar%ZMGC=trcSalt_solml(idsalt_MgCO3,L,NY,NX)
  chemvar%ZMGH=trcSalt_solml(idsalt_MgHCO3,L,NY,NX)
  chemvar%ZMGS=trcSalt_solml(idsalt_MgSO4,L,NY,NX)
  chemvar%ZNAC=trcSalt_solml(idsalt_NaCO3,L,NY,NX)
  chemvar%ZNAS=trcSalt_solml(idsalt_NaSO4,L,NY,NX)
  chemvar%ZKAS=trcSalt_solml(idsalt_KSO4,L,NY,NX)
  chemvar%H0PO4=trcSalt_solml(idsalt_H0PO4,L,NY,NX)
  chemvar%H3PO4=trcSalt_solml(idsalt_H3PO4,L,NY,NX)
  chemvar%ZFE1P=trcSalt_solml(idsalt_FeHPO4,L,NY,NX)
  chemvar%ZFE2P=trcSalt_solml(idsalt_FeH2PO4,L,NY,NX)
  chemvar%ZCA0P=trcSalt_solml(idsalt_CaPO4,L,NY,NX)
  chemvar%ZCA1P=trcSalt_solml(idsalt_CaHPO4,L,NY,NX)
  chemvar%ZCA2P=trcSalt_solml(idsalt_CaH2PO4,L,NY,NX)
  chemvar%ZMG1P=trcSalt_solml(idsalt_MgHPO4,L,NY,NX)
  chemvar%H0POB=trcSalt_solml(idsalt_H0PO4B,L,NY,NX)
  chemvar%H3POB=trcSalt_solml(idsalt_H3PO4B,L,NY,NX)
  chemvar%ZFE1PB=trcSalt_solml(idsalt_FeHPO4B,L,NY,NX)
  chemvar%ZFE2PB=trcSalt_solml(idsalt_FeH2PO4B,L,NY,NX)
  chemvar%ZCA0PB=trcSalt_solml(idsalt_CaPO4B,L,NY,NX)
  chemvar%ZCA1PB=trcSalt_solml(idsalt_CaHPO4B,L,NY,NX)
  chemvar%ZCA2PB=trcSalt_solml(idsalt_CaH2PO4B,L,NY,NX)
  chemvar%ZMG1PB=trcSalt_solml(idsalt_MgHPO4B,L,NY,NX)
  chemvar%XHY=trcx_solml(idx_Hp,L,NY,NX)
  chemvar%XAL=trcx_solml(idx_Al,L,NY,NX)
  chemvar%XFE=trcx_solml(idx_Fe,L,NY,NX)
  chemvar%XCA=trcx_solml(idx_Ca,L,NY,NX)
  chemvar%XMG=trcx_solml(idx_Mg,L,NY,NX)
  chemvar%XNA=trcx_solml(idx_Na,L,NY,NX)
  chemvar%XKA=trcx_solml(idx_K,L,NY,NX)
  chemvar%XHC=trcx_solml(idx_COOH,L,NY,NX)
  chemvar%XALO2=trcx_solml(idx_AlOH2,L,NY,NX)
  chemvar%XFEO2=trcx_solml(idx_FeOH2,L,NY,NX)
  chemvar%ORGC=ORGC(L,NY,NX)
  chemvar%PALOH=trcp_salml(idsp_AlOH3,L,NY,NX)
  chemvar%PFEOH=trcp_salml(idsp_FeOH3,L,NY,NX)
  chemvar%PCACO=trcp_salml(idsp_CaCO3,L,NY,NX)
  chemvar%PCASO=trcp_salml(idsp_CaSO4,L,NY,NX)
  chemvar%VLPOB=trcs_VLN(ids_H1PO4B,L,NY,NX)
  chemvar%VLPO4=trcs_VLN(ids_H1PO4,L,NY,NX)
  chemvar%VLNOB=trcs_VLN(ids_NO3B,L,NY,NX)
  chemvar%VLNO3=trcs_VLN(ids_NO3,L,NY,NX)

!  solflx%TR_CaCO3_precip_soil=trcp_TR(idsp_CaCO3,L,NY,NX)
!  solflx%TR_NaCO3_soil=trcSalt_TR(idsalt_NaCO3,L,NY,NX)
!  solflx%TR_MgCO3_soil=trcSalt_TR(idsalt_MgCO3,L,NY,NX)
!  solflx%TR_CaCO3_soil=trcSalt_TR(idsalt_CaCO3,L,NY,NX)
!  solflx%TR_MgHCO3_soil=trcSalt_TR(idsalt_MgHCO3,L,NY,NX)
!  solflx%TR_CaHCO3_soil=trcSalt_TR(idsalt_CaHCO3,L,NY,NX)
!  solflx%TRHCO=trcSalt_TR(idsalt_HCO3,L,NY,NX)
!  solflx%TR_HCO3_sorbed_soil=TR_HCO3_sorbed_soil(L,NY,NX)
!  solflx%TR_CaH2PO4_soil=trcSalt_TR(idsalt_CaH2PO4,L,NY,NX)
!  solflx%TR_FeH2PO4_soil=trcSalt_TR(idsalt_FeH2PO4,L,NY,NX)
!  solflx%TR_CaH2PO4_band_soil=trcSalt_TR(idsalt_CaH2PO4B,L,NY,NX)
!  solflx%TR_FeH2PO4_band_soil=trcSalt_TR(idsalt_FeH2PO4B,L,NY,NX)
!  solflx%TR_MgHPO4_soil=trcSalt_TR(idsalt_MgHPO4,L,NY,NX)
!  solflx%TR_CaHPO4_soil=trcSalt_TR(idsalt_CaHPO4,L,NY,NX)
!  solflx%TR_FeHPO4_soil=trcSalt_TR(idsalt_FeHPO4,L,NY,NX)
!  solflx%TR_MgHPO4_band_soil=trcSalt_TR(idsalt_MgHPO4B,L,NY,NX)
!  solflx%TR_CaHPO4_band_soil=trcSalt_TR(idsalt_CaHPO4B,L,NY,NX)
!  solflx%TR_FeHPO4_band_soil=trcSalt_TR(idsalt_FeHPO4B,L,NY,NX)
!  solflx%TR_CaPO4_band_soil=trcSalt_TR(idsalt_CaPO4B,L,NY,NX)
!  solflx%TR_PO4_band_soil=trcSalt_TR(idsalt_H0PO4B,L,NY,NX)
!  solflx%TR_CaPO4_soil=trcSalt_TR(idsalt_CaPO4,L,NY,NX)
!  solflx%TR_PO4_soil=trcSalt_TR(idsalt_H0PO4,L,NY,NX)
!  solflx%TR_FePO4_precip_soil=trcp_TR(idsp_FePO4,L,NY,NX)
!  solflx%TR_AlPO4_precip_soil=trcp_TR(idsp_AlPO4,L,NY,NX)
!  solflx%TR_FePO4_precip_band_soil=trcp_TR(idsp_FePO4B,L,NY,NX)
!  solflx%TR_AlPO4_precip_band_soil=trcp_TR(idsp_AlPO4B,L,NY,NX)
!  solflx%TR_CaH4P2O8_precip_soil=trcp_TR(idsp_CaH2PO4,L,NY,NX)
!  solflx%TR_CaHPO4_precip_soil=trcp_TR(idsp_CaHPO4,L,NY,NX)
!  solflx%TR_CaH4P2O8_precip_band_soil=trcp_TR(idsp_CaH2PO4B,L,NY,NX)
!  solflx%TR_CaHPO4_precip_band_soil=trcp_TR(idsp_CaHPO4B,L,NY,NX)
!  solflx%TR_apatite_precip_band_soil=trcp_TR(idsp_HAB,L,NY,NX)
!  solflx%TR_apatite_precip_soil=trcp_TR(idsp_HA,L,NY,NX)
!  solflx%TR_FeOH_soil=trcSalt_TR(idsalt_FeOH,L,NY,NX)
!  solflx%TR_Al_3p_soil=trcSalt_TR(idsalt_Al,L,NY,NX)
!  solflx%TR_AlOH_soil=trcSalt_TR(idsalt_AlOH,L,NY,NX)
!  solflx%TR_FeO2H2_soil=trcSalt_TR(idsalt_FeOH2,L,NY,NX)
!  solflx%TR_AlO2H2_soil=trcSalt_TR(idsalt_AlOH2,L,NY,NX)
!  solflx%TR_FeO3H3_soil=trcSalt_TR(idsalt_FeOH3,L,NY,NX)
!  solflx%TR_AlO3H3_soil=trcSalt_TR(idsalt_AlOH3,L,NY,NX)
!  solflx%TR_FeO4H4_soil=trcSalt_TR(idsalt_FeOH4,L,NY,NX)
!  solflx%TR_AlO4H4_soil=trcSalt_TR(idsalt_AlOH4,L,NY,NX)
!  solflx%TR_MgOH_soil=trcSalt_TR(idsalt_MgOH2,L,NY,NX)
!  solflx%TR_CaOH_soil=trcSalt_TR(idsalt_CaOH2,L,NY,NX)
!  solflx%TR_FeOH3_precip_soil=trcp_TR(idsp_FeOH3,L,NY,NX)
!  solflx%TR_AlOH3_precip_soil=trcp_TR(idsp_AlOH3,L,NY,NX)
!  solflx%TR_FeO2H2_sorbed_soil=TR_FeO2H2_sorbed_soil(L,NY,NX)
!  solflx%TR_AlSO4_soil =trcSalt_TR(idsalt_AlSO4,L,NY,NX)
!  solflx%TR_RHPO4_sorbed_band_soil =trcx_TR(idx_HPO4B,L,NY,NX)
!  solflx%TR_RH2PO4_sorbed_band_soil =trcx_TR(idx_H2PO4B,L,NY,NX)
!  solflx%TR_RO_sorbed_band_soil =trcx_TR(idx_OHeB,L,NY,NX)
!  solflx%TR_ROH_sorbed_band_soil =trcx_TR(idx_OHB,L,NY,NX)
!  solflx%TR_ROH2_sorbed_band_soil =trcx_TR(idx_OHpB,L,NY,NX)
!  solflx%TR_Ca_2p_soil  =trcSalt_TR(idsalt_Ca,L,NY,NX)
!  solflx%TR_CaSO4_soil =trcSalt_TR(idsalt_CaSO4,L,NY,NX)
!  solflx%TR_CaSO4_precip_soil=trcp_TR(idsp_CaSO4,L,NY,NX)
!  solflx%TR_CO2_aqu_soil =TR_CO2_aqu_soil(L,NY,NX)
!  solflx%TR_CO3_2e_soil =trcSalt_TR(idsalt_CO3,L,NY,NX)
!  solflx%TR_Fe_3p_soil  =trcSalt_TR(idsalt_Fe,L,NY,NX)
!  solflx%TR_FeSO4_soil =trcSalt_TR(idsalt_FeSO4,L,NY,NX)
!  solflx%TR_H1PO4_band_soil =TR_H1PO4_band_soil(L,NY,NX)
!  solflx%TR_H2PO4_band_soil =TR_H2PO4_band_soil(L,NY,NX)
!  solflx%TR_H1PO4_soil =TR_H1PO4_soil(L,NY,NX)
!  solflx%TR_H2PO4_soil =TR_H2PO4_soil(L,NY,NX)
!  solflx%TR_H3PO4_band_soil =trcSalt_TR(idsalt_H3PO4B,L,NY,NX)
!  solflx%TR_H3PO4_sorbed_soil =trcSalt_TR(idsalt_H3PO4,L,NY,NX)
!  solflx%TR_H_p_soil  =trcSalt_TR(idsalt_Hp,L,NY,NX)
!  solflx%TR_K_1p_soil  =trcSalt_TR(idsalt_K,L,NY,NX)
!  solflx%TR_KSO4_soil =trcSalt_TR(idsalt_KSO4,L,NY,NX)
!  solflx%TR_Mg_2p_soil  =trcSalt_TR(idsalt_Mg,L,NY,NX)
!  solflx%TR_MgSO4_soil =trcSalt_TR(idsalt_MgSO4,L,NY,NX)
!  solflx%TR_NH3_band_soil =TR_NH3_band_soil(L,NY,NX)
!  solflx%TR_NH3_soil =TR_NH3_soil(L,NY,NX)
!  solflx%TR_NH4_band_soil =TR_NH4_band_soil(L,NY,NX)
!  solflx%TR_NH4_soil =TR_NH4_soil(L,NY,NX)
!  solflx%TR_Na_p_soil  =trcSalt_TR(idsalt_Na,L,NY,NX)
!  solflx%TR_NaSO4_soil =trcSalt_TR(idsalt_NaSO4,L,NY,NX)
!  solflx%TR_OH_1e_soil  =trcSalt_TR(idsalt_OH,L,NY,NX)
!  solflx%TR_SO4_2e_soil =trcSalt_TR(idsalt_SO4,L,NY,NX)
!  solflx%TR_RHPO4_sorbed_soil =trcx_TR(idx_HPO4,L,NY,NX)
!  solflx%TR_RH2PO4_sorbed_soil =trcx_TR(idx_H2PO4,L,NY,NX)
!  solflx%TR_Al_sorbed_soil =TR_Al_sorbed_soil(L,NY,NX)
!  solflx%TR_AlO2H2_sorbed_soil=TR_AlO2H2_sorbed_soil(L,NY,NX)
!  solflx%TR_Ca_sorbed_soil =TR_Ca_sorbed_soil(L,NY,NX)
!  solflx%TR_Fe_sorbed_soil =TR_Fe_sorbed_soil(L,NY,NX)
!  solflx%TR_RO_sorbed_soil =trcx_TR(idx_OHe,L,NY,NX)
!  solflx%TR_ROH_sorbed_soil =trcx_TR(idx_OH,L,NY,NX)
!  solflx%TR_ROH2_sorbed_soil =trcx_TR(idx_OHp,L,NY,NX)
!  solflx%TR_H_p_sorbed_soil =TR_H_p_sorbed_soil(L,NY,NX)
!  solflx%TR_K_sorbed_soil =TR_K_sorbed_soil(L,NY,NX)
!  solflx%TR_Mg_sorbed_soil =TR_Mg_sorbed_soil(L,NY,NX)
!  solflx%TR_NH4_sorbed_soil =trcx_TR(idx_NH4,L,NY,NX)
!  solflx%TR_Na_sorbed_soil =TR_Na_sorbed_soil(L,NY,NX)
!  solflx%TR_NH4_sorbed_band_soil =trcx_TR(idx_NH4B,L,NY,NX)
!  solflx%TBCO2 =TBCO2(L,NY,NX)
!  solflx%TBION =TBION(L,NY,NX)
!  solflx%TRH2O =TRH2O(L,NY,NX)
!  Integers

  end subroutine GeochemAPISend

!-----------------------------------------------------------------------

  subroutine GeochemAPIRecv(L,NY,NX,solflx)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(solute_flx_type), intent(in) :: solflx

  TR_NH4_soil(L,NY,NX)=solflx%TR_NH4_soil
  TR_NH4_band_soil(L,NY,NX)=solflx%TR_NH4_band_soil
  TR_NH3_soil(L,NY,NX)=solflx%TR_NH3_soil
  TR_NH3_band_soil(L,NY,NX)=solflx%TR_NH3_band_soil
  TR_H1PO4_soil(L,NY,NX)=solflx%TR_H1PO4_soil
  TR_H2PO4_soil(L,NY,NX)=solflx%TR_H2PO4_soil
  TR_H1PO4_band_soil(L,NY,NX)=solflx%TR_H1PO4_band_soil
  TR_H2PO4_band_soil(L,NY,NX)=solflx%TR_H2PO4_band_soil
  trcx_TR(idx_NH4,L,NY,NX)=solflx%TR_NH4_sorbed_soil
  trcx_TR(idx_NH4B,L,NY,NX)=solflx%TR_NH4_sorbed_band_soil
  trcx_TR(idx_OH,L,NY,NX)=solflx%TR_ROH_sorbed_soil
  trcx_TR(idx_OHp,L,NY,NX)=solflx%TR_ROH2_sorbed_soil
  trcx_TR(idx_HPO4,L,NY,NX)=solflx%TR_RHPO4_sorbed_soil
  trcx_TR(idx_H2PO4,L,NY,NX)=solflx%TR_RH2PO4_sorbed_soil
  trcx_TR(idx_OHB,L,NY,NX)=solflx%TR_ROH_sorbed_band_soil
  trcx_TR(idx_OHpB,L,NY,NX)=solflx%TR_ROH2_sorbed_band_soil
  trcx_TR(idx_HPO4B,L,NY,NX)=solflx%TR_RHPO4_sorbed_band_soil
  trcx_TR(idx_H2PO4B,L,NY,NX)=solflx%TR_RH2PO4_sorbed_band_soil
  trcp_TR(idsp_AlPO4,L,NY,NX)=solflx%TR_AlPO4_precip_soil
  trcp_TR(idsp_FePO4,L,NY,NX)=solflx%TR_FePO4_precip_soil
  trcp_TR(idsp_CaHPO4,L,NY,NX)=solflx%TR_CaHPO4_precip_soil
  trcp_TR(idsp_HA,L,NY,NX)=solflx%TR_apatite_precip_soil
  trcp_TR(idsp_CaH2PO4,L,NY,NX)=solflx%TR_CaH4P2O8_precip_soil
  trcp_TR(idsp_AlPO4B,L,NY,NX)=solflx%TR_AlPO4_precip_band_soil
  trcp_TR(idsp_FePO4B,L,NY,NX)=solflx%TR_FePO4_precip_band_soil
  trcp_TR(idsp_CaHPO4B,L,NY,NX)=solflx%TR_CaHPO4_precip_band_soil
  trcp_TR(idsp_HAB,L,NY,NX)=solflx%TR_apatite_precip_band_soil
  trcp_TR(idsp_CaH2PO4B,L,NY,NX)=solflx%TR_CaH4P2O8_precip_band_soil
  trcSalt_TR(idsalt_Al,L,NY,NX)=solflx%TR_Al_3p_soil
  trcSalt_TR(idsalt_Fe,L,NY,NX)=solflx%TR_Fe_3p_soil
  trcSalt_TR(idsalt_Hp,L,NY,NX)=solflx%TR_H_p_soil
  trcSalt_TR(idsalt_Ca,L,NY,NX)=solflx%TR_Ca_2p_soil
  trcSalt_TR(idsalt_Mg,L,NY,NX)=solflx%TR_Mg_2p_soil
  trcSalt_TR(idsalt_Na,L,NY,NX)=solflx%TR_Na_p_soil
  trcSalt_TR(idsalt_K,L,NY,NX)=solflx%TR_K_1p_soil
  trcSalt_TR(idsalt_OH,L,NY,NX)=solflx%TR_OH_1e_soil
  trcSalt_TR(idsalt_SO4,L,NY,NX)=solflx%TR_SO4_2e_soil
  trcSalt_TR(idsalt_CO3,L,NY,NX)=solflx%TR_CO3_2e_soil
  trcSalt_TR(idsalt_HCO3,L,NY,NX)=solflx%TRHCO
  TR_CO2_aqu_soil(L,NY,NX)=solflx%TR_CO2_aqu_soil
  trcSalt_TR(idsalt_AlOH,L,NY,NX)=solflx%TR_AlOH_soil
  trcSalt_TR(idsalt_AlOH2,L,NY,NX)=solflx%TR_AlO2H2_soil
  trcSalt_TR(idsalt_AlOH3,L,NY,NX)=solflx%TR_AlO3H3_soil
  trcSalt_TR(idsalt_AlOH4,L,NY,NX)=solflx%TR_AlO4H4_soil
  trcSalt_TR(idsalt_AlSO4,L,NY,NX)=solflx%TR_AlSO4_soil
  trcSalt_TR(idsalt_FeOH,L,NY,NX)=solflx%TR_FeOH_soil
  trcSalt_TR(idsalt_FeOH2,L,NY,NX)=solflx%TR_FeO2H2_soil
  trcSalt_TR(idsalt_FeOH3,L,NY,NX)=solflx%TR_FeO3H3_soil
  trcSalt_TR(idsalt_FeOH4,L,NY,NX)=solflx%TR_FeO4H4_soil
  trcSalt_TR(idsalt_FeSO4,L,NY,NX)=solflx%TR_FeSO4_soil
  trcSalt_TR(idsalt_CaOH2,L,NY,NX)=solflx%TR_CaOH_soil
  trcSalt_TR(idsalt_CaCO3,L,NY,NX)=solflx%TR_CaCO3_soil
  trcSalt_TR(idsalt_CaHCO3,L,NY,NX)=solflx%TR_CaHCO3_soil
  trcSalt_TR(idsalt_CaSO4,L,NY,NX)=solflx%TR_CaSO4_soil
  trcSalt_TR(idsalt_MgOH2,L,NY,NX)=solflx%TR_MgOH_soil
  trcSalt_TR(idsalt_MgCO3,L,NY,NX)=solflx%TR_MgCO3_soil
  trcSalt_TR(idsalt_MgHCO3,L,NY,NX)=solflx%TR_MgHCO3_soil
  trcSalt_TR(idsalt_MgSO4,L,NY,NX)=solflx%TR_MgSO4_soil
  trcSalt_TR(idsalt_NaCO3,L,NY,NX)=solflx%TR_NaCO3_soil
  trcSalt_TR(idsalt_NaSO4,L,NY,NX)=solflx%TR_NaSO4_soil
  trcSalt_TR(idsalt_KSO4,L,NY,NX)=solflx%TR_KSO4_soil
  trcSalt_TR(idsalt_H0PO4,L,NY,NX)=solflx%TR_PO4_soil
  trcSalt_TR(idsalt_H3PO4,L,NY,NX)=solflx%TR_H3PO4_sorbed_soil
  trcSalt_TR(idsalt_FeHPO4,L,NY,NX)=solflx%TR_FeHPO4_soil
  trcSalt_TR(idsalt_FeH2PO4,L,NY,NX)=solflx%TR_FeH2PO4_soil
  trcSalt_TR(idsalt_CaPO4,L,NY,NX)=solflx%TR_CaPO4_soil
  trcSalt_TR(idsalt_CaHPO4,L,NY,NX)=solflx%TR_CaHPO4_soil
  trcSalt_TR(idsalt_CaH2PO4,L,NY,NX)=solflx%TR_CaH2PO4_soil
  trcSalt_TR(idsalt_MgHPO4,L,NY,NX)=solflx%TR_MgHPO4_soil
  trcSalt_TR(idsalt_H0PO4B,L,NY,NX)=solflx%TR_PO4_band_soil
  trcSalt_TR(idsalt_H3PO4B,L,NY,NX)=solflx%TR_H3PO4_band_soil
  trcSalt_TR(idsalt_FeHPO4B,L,NY,NX)=solflx%TR_FeHPO4_band_soil
  trcSalt_TR(idsalt_FeH2PO4B,L,NY,NX)=solflx%TR_FeH2PO4_band_soil
  trcSalt_TR(idsalt_CaPO4B,L,NY,NX)=solflx%TR_CaPO4_band_soil
  trcSalt_TR(idsalt_CaHPO4B,L,NY,NX)=solflx%TR_CaHPO4_band_soil
  trcSalt_TR(idsalt_CaH2PO4B,L,NY,NX)=solflx%TR_CaH2PO4_band_soil
  trcSalt_TR(idsalt_MgHPO4B,L,NY,NX)=solflx%TR_MgHPO4_band_soil
  TR_H_p_sorbed_soil(L,NY,NX)=solflx%TR_H_p_sorbed_soil
  TR_Al_sorbed_soil(L,NY,NX)=solflx%TR_Al_sorbed_soil
  TR_Fe_sorbed_soil(L,NY,NX)=solflx%TR_Fe_sorbed_soil
  TR_Ca_sorbed_soil(L,NY,NX)=solflx%TR_Ca_sorbed_soil
  TR_Mg_sorbed_soil(L,NY,NX)=solflx%TR_Mg_sorbed_soil
  TR_Na_sorbed_soil(L,NY,NX)=solflx%TR_Na_sorbed_soil
  TR_K_sorbed_soil(L,NY,NX)=solflx%TR_K_sorbed_soil
  TR_HCO3_sorbed_soil(L,NY,NX)=solflx%TR_HCO3_sorbed_soil
  TR_AlO2H2_sorbed_soil(L,NY,NX)=solflx%TR_AlO2H2_sorbed_soil
  TR_FeO2H2_sorbed_soil(L,NY,NX)=solflx%TR_FeO2H2_sorbed_soil
  trcx_TR(idx_OHe,L,NY,NX)=solflx%TR_RO_sorbed_soil
  trcx_TR(idx_OHeB,L,NY,NX)=solflx%TR_RO_sorbed_band_soil
  trcp_TR(idsp_AlOH3,L,NY,NX)=solflx%TR_AlOH3_precip_soil
  trcp_TR(idsp_FeOH3,L,NY,NX)=solflx%TR_FeOH3_precip_soil
  trcp_TR(idsp_CaCO3,L,NY,NX)=solflx%TR_CaCO3_precip_soil
  trcp_TR(idsp_CaSO4,L,NY,NX)=solflx%TR_CaSO4_precip_soil
  TRH2O(L,NY,NX)=solflx%TRH2O
  TBION(L,NY,NX)=solflx%TBION
  TBCO2(L,NY,NX)=solflx%TBCO2
  end subroutine GeochemAPIRecv

end module GeochemAPI
