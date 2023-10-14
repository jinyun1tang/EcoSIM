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

!  solflx%TRCACO=trcp_TR(idsp_CaCO3,L,NY,NX)
!  solflx%TRNAC=trcSalt_TR(idsalt_NaCO3,L,NY,NX)
!  solflx%TRMGC=trcSalt_TR(idsalt_MgCO3,L,NY,NX)
!  solflx%TRCAC=trcSalt_TR(idsalt_CaCO3,L,NY,NX)
!  solflx%TRMGH=trcSalt_TR(idsalt_MgHCO3,L,NY,NX)
!  solflx%TRCAH=trcSalt_TR(idsalt_CaHCO3,L,NY,NX)
!  solflx%TRHCO=trcSalt_TR(idsalt_HCO3,L,NY,NX)
!  solflx%TRXHC=TRXHC(L,NY,NX)
!  solflx%TRC2P=trcSalt_TR(idsalt_CaH2PO4,L,NY,NX)
!  solflx%TRF2P=trcSalt_TR(idsalt_FeH2PO4,L,NY,NX)
!  solflx%TRC2B=trcSalt_TR(idsalt_CaH2PO4B,L,NY,NX)
!  solflx%TRF2B=trcSalt_TR(idsalt_FeH2PO4B,L,NY,NX)
!  solflx%TRM1P=trcSalt_TR(idsalt_MgHPO4,L,NY,NX)
!  solflx%TRC1P=trcSalt_TR(idsalt_CaHPO4,L,NY,NX)
!  solflx%TRF1P=trcSalt_TR(idsalt_FeHPO4,L,NY,NX)
!  solflx%TRM1B=trcSalt_TR(idsalt_MgHPO4B,L,NY,NX)
!  solflx%TRC1B=trcSalt_TR(idsalt_CaHPO4B,L,NY,NX)
!  solflx%TRF1B=trcSalt_TR(idsalt_FeHPO4B,L,NY,NX)
!  solflx%TRC0B=trcSalt_TR(idsalt_CaPO4B,L,NY,NX)
!  solflx%TRH0B=trcSalt_TR(idsalt_H0PO4B,L,NY,NX)
!  solflx%TRC0P=trcSalt_TR(idsalt_CaPO4,L,NY,NX)
!  solflx%TRH0P=trcSalt_TR(idsalt_H0PO4,L,NY,NX)
!  solflx%TRFEPO=trcp_TR(idsp_FePO4,L,NY,NX)
!  solflx%TR_AlPO4=trcp_TR(idsp_AlPO4,L,NY,NX)
!  solflx%TRFEPB=trcp_TR(idsp_FePO4B,L,NY,NX)
!  solflx%TRALPB=trcp_TR(idsp_AlPO4B,L,NY,NX)
!  solflx%TRCAPM=trcp_TR(idsp_CaH2PO4,L,NY,NX)
!  solflx%TRCAPD=trcp_TR(idsp_CaHPO4,L,NY,NX)
!  solflx%TRCPMB=trcp_TR(idsp_CaH2PO4B,L,NY,NX)
!  solflx%TRCPDB=trcp_TR(idsp_CaHPO4B,L,NY,NX)
!  solflx%TRCPHB=trcp_TR(idsp_HAB,L,NY,NX)
!  solflx%TRCAPH=trcp_TR(idsp_HA,L,NY,NX)
!  solflx%TRFE1=trcSalt_TR(idsalt_FeOH,L,NY,NX)
!  solflx%TRAL=trcSalt_TR(idsalt_Al,L,NY,NX)
!  solflx%TRAL1=trcSalt_TR(idsalt_AlOH,L,NY,NX)
!  solflx%TRFE2=trcSalt_TR(idsalt_FeOH2,L,NY,NX)
!  solflx%TRAL2=trcSalt_TR(idsalt_AlOH2,L,NY,NX)
!  solflx%TRFE3=trcSalt_TR(idsalt_FeOH3,L,NY,NX)
!  solflx%TRAL3=trcSalt_TR(idsalt_AlOH3,L,NY,NX)
!  solflx%TRFE4=trcSalt_TR(idsalt_FeOH4,L,NY,NX)
!  solflx%TRAL4=trcSalt_TR(idsalt_AlOH4,L,NY,NX)
!  solflx%TRMGO=trcSalt_TR(idsalt_MgOH2,L,NY,NX)
!  solflx%TRCAO=trcSalt_TR(idsalt_CaOH2,L,NY,NX)
!  solflx%TRFEOH=trcp_TR(idsp_FeOH3,L,NY,NX)
!  solflx%TRALOH=trcp_TR(idsp_AlOH3,L,NY,NX)
!  solflx%TRXFE2=TRXFE2(L,NY,NX)
!  solflx%TRALS =trcSalt_TR(idsalt_AlSO4,L,NY,NX)
!  solflx%TRB1P =trcx_TR(idx_HPO4B,L,NY,NX)
!  solflx%TRB2P =trcx_TR(idx_H2PO4B,L,NY,NX)
!  solflx%TRBH0 =trcx_TR(idx_OHeB,L,NY,NX)
!  solflx%TRBH1 =trcx_TR(idx_OHB,L,NY,NX)
!  solflx%TRBH2 =trcx_TR(idx_OHpB,L,NY,NX)
!  solflx%TRCA  =trcSalt_TR(idsalt_Ca,L,NY,NX)
!  solflx%TRCAS =trcSalt_TR(idsalt_CaSO4,L,NY,NX)
!  solflx%TRCASO=trcp_TR(idsp_CaSO4,L,NY,NX)
!  solflx%TRCO2 =TRCO2(L,NY,NX)
!  solflx%TRCO3 =trcSalt_TR(idsalt_CO3,L,NY,NX)
!  solflx%TRFE  =trcSalt_TR(idsalt_Fe,L,NY,NX)
!  solflx%TRFES =trcSalt_TR(idsalt_FeSO4,L,NY,NX)
!  solflx%TRH1B =TRH1B(L,NY,NX)
!  solflx%TRH2B =TRH2B(L,NY,NX)
!  solflx%TRH1P =TRH1P(L,NY,NX)
!  solflx%TRH2P =TRH2P(L,NY,NX)
!  solflx%TRH3B =trcSalt_TR(idsalt_H3PO4B,L,NY,NX)
!  solflx%TRH3P =trcSalt_TR(idsalt_H3PO4,L,NY,NX)
!  solflx%TRHY  =trcSalt_TR(idsalt_Hp,L,NY,NX)
!  solflx%TRKA  =trcSalt_TR(idsalt_K,L,NY,NX)
!  solflx%TRKAS =trcSalt_TR(idsalt_KSO4,L,NY,NX)
!  solflx%TRMG  =trcSalt_TR(idsalt_Mg,L,NY,NX)
!  solflx%TRMGS =trcSalt_TR(idsalt_MgSO4,L,NY,NX)
!  solflx%TRN3B =TRN3B(L,NY,NX)
!  solflx%TRN3S =TRN3S(L,NY,NX)
!  solflx%TRN4B =TRN4B(L,NY,NX)
!  solflx%TRN4S =TRN4S(L,NY,NX)
!  solflx%TRNA  =trcSalt_TR(idsalt_Na,L,NY,NX)
!  solflx%TRNAS =trcSalt_TR(idsalt_NaSO4,L,NY,NX)
!  solflx%TROH  =trcSalt_TR(idsalt_OH,L,NY,NX)
!  solflx%TRSO4 =trcSalt_TR(idsalt_SO4,L,NY,NX)
!  solflx%TRX1P =trcx_TR(idx_HPO4,L,NY,NX)
!  solflx%TRX2P =trcx_TR(idx_H2PO4,L,NY,NX)
!  solflx%TRXAL =TRXAL(L,NY,NX)
!  solflx%TRXAL2=TRXAL2(L,NY,NX)
!  solflx%TRXCA =TRXCA(L,NY,NX)
!  solflx%TRXFE =TRXFE(L,NY,NX)
!  solflx%TRXH0 =trcx_TR(idx_OHe,L,NY,NX)
!  solflx%TRXH1 =trcx_TR(idx_OH,L,NY,NX)
!  solflx%TRXH2 =trcx_TR(idx_OHp,L,NY,NX)
!  solflx%TRXHY =TRXHY(L,NY,NX)
!  solflx%TRXKA =TRXKA(L,NY,NX)
!  solflx%TRXMG =TRXMG(L,NY,NX)
!  solflx%TRXN4 =trcx_TR(idx_NH4,L,NY,NX)
!  solflx%TRXNA =TRXNA(L,NY,NX)
!  solflx%TRXNB =trcx_TR(idx_NH4B,L,NY,NX)
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

  TRN4S(L,NY,NX)=solflx%TRN4S
  TRN4B(L,NY,NX)=solflx%TRN4B
  TRN3S(L,NY,NX)=solflx%TRN3S
  TRN3B(L,NY,NX)=solflx%TRN3B
  TRH1P(L,NY,NX)=solflx%TRH1P
  TRH2P(L,NY,NX)=solflx%TRH2P
  TRH1B(L,NY,NX)=solflx%TRH1B
  TRH2B(L,NY,NX)=solflx%TRH2B
  trcx_TR(idx_NH4,L,NY,NX)=solflx%TRXN4
  trcx_TR(idx_NH4B,L,NY,NX)=solflx%TRXNB
  trcx_TR(idx_OH,L,NY,NX)=solflx%TRXH1
  trcx_TR(idx_OHp,L,NY,NX)=solflx%TRXH2
  trcx_TR(idx_HPO4,L,NY,NX)=solflx%TRX1P
  trcx_TR(idx_H2PO4,L,NY,NX)=solflx%TRX2P
  trcx_TR(idx_OHB,L,NY,NX)=solflx%TRBH1
  trcx_TR(idx_OHpB,L,NY,NX)=solflx%TRBH2
  trcx_TR(idx_HPO4B,L,NY,NX)=solflx%TRB1P
  trcx_TR(idx_H2PO4B,L,NY,NX)=solflx%TRB2P
  trcp_TR(idsp_AlPO4,L,NY,NX)=solflx%TR_AlPO4
  trcp_TR(idsp_FePO4,L,NY,NX)=solflx%TRFEPO
  trcp_TR(idsp_CaHPO4,L,NY,NX)=solflx%TRCAPD
  trcp_TR(idsp_HA,L,NY,NX)=solflx%TRCAPH
  trcp_TR(idsp_CaH2PO4,L,NY,NX)=solflx%TRCAPM
  trcp_TR(idsp_AlPO4B,L,NY,NX)=solflx%TRALPB
  trcp_TR(idsp_FePO4B,L,NY,NX)=solflx%TRFEPB
  trcp_TR(idsp_CaHPO4B,L,NY,NX)=solflx%TRCPDB
  trcp_TR(idsp_HAB,L,NY,NX)=solflx%TRCPHB
  trcp_TR(idsp_CaH2PO4B,L,NY,NX)=solflx%TRCPMB
  trcSalt_TR(idsalt_Al,L,NY,NX)=solflx%TRAL
  trcSalt_TR(idsalt_Fe,L,NY,NX)=solflx%TRFE
  trcSalt_TR(idsalt_Hp,L,NY,NX)=solflx%TRHY
  trcSalt_TR(idsalt_Ca,L,NY,NX)=solflx%TRCA
  trcSalt_TR(idsalt_Mg,L,NY,NX)=solflx%TRMG
  trcSalt_TR(idsalt_Na,L,NY,NX)=solflx%TRNA
  trcSalt_TR(idsalt_K,L,NY,NX)=solflx%TRKA
  trcSalt_TR(idsalt_OH,L,NY,NX)=solflx%TROH
  trcSalt_TR(idsalt_SO4,L,NY,NX)=solflx%TRSO4
  trcSalt_TR(idsalt_CO3,L,NY,NX)=solflx%TRCO3
  trcSalt_TR(idsalt_HCO3,L,NY,NX)=solflx%TRHCO
  TRCO2(L,NY,NX)=solflx%TRCO2
  trcSalt_TR(idsalt_AlOH,L,NY,NX)=solflx%TRAL1
  trcSalt_TR(idsalt_AlOH2,L,NY,NX)=solflx%TRAL2
  trcSalt_TR(idsalt_AlOH3,L,NY,NX)=solflx%TRAL3
  trcSalt_TR(idsalt_AlOH4,L,NY,NX)=solflx%TRAL4
  trcSalt_TR(idsalt_AlSO4,L,NY,NX)=solflx%TRALS
  trcSalt_TR(idsalt_FeOH,L,NY,NX)=solflx%TRFE1
  trcSalt_TR(idsalt_FeOH2,L,NY,NX)=solflx%TRFE2
  trcSalt_TR(idsalt_FeOH3,L,NY,NX)=solflx%TRFE3
  trcSalt_TR(idsalt_FeOH4,L,NY,NX)=solflx%TRFE4
  trcSalt_TR(idsalt_FeSO4,L,NY,NX)=solflx%TRFES
  trcSalt_TR(idsalt_CaOH2,L,NY,NX)=solflx%TRCAO
  trcSalt_TR(idsalt_CaCO3,L,NY,NX)=solflx%TRCAC
  trcSalt_TR(idsalt_CaHCO3,L,NY,NX)=solflx%TRCAH
  trcSalt_TR(idsalt_CaSO4,L,NY,NX)=solflx%TRCAS
  trcSalt_TR(idsalt_MgOH2,L,NY,NX)=solflx%TRMGO
  trcSalt_TR(idsalt_MgCO3,L,NY,NX)=solflx%TRMGC
  trcSalt_TR(idsalt_MgHCO3,L,NY,NX)=solflx%TRMGH
  trcSalt_TR(idsalt_MgSO4,L,NY,NX)=solflx%TRMGS
  trcSalt_TR(idsalt_NaCO3,L,NY,NX)=solflx%TRNAC
  trcSalt_TR(idsalt_NaSO4,L,NY,NX)=solflx%TRNAS
  trcSalt_TR(idsalt_KSO4,L,NY,NX)=solflx%TRKAS
  trcSalt_TR(idsalt_H0PO4,L,NY,NX)=solflx%TRH0P
  trcSalt_TR(idsalt_H3PO4,L,NY,NX)=solflx%TRH3P
  trcSalt_TR(idsalt_FeHPO4,L,NY,NX)=solflx%TRF1P
  trcSalt_TR(idsalt_FeH2PO4,L,NY,NX)=solflx%TRF2P
  trcSalt_TR(idsalt_CaPO4,L,NY,NX)=solflx%TRC0P
  trcSalt_TR(idsalt_CaHPO4,L,NY,NX)=solflx%TRC1P
  trcSalt_TR(idsalt_CaH2PO4,L,NY,NX)=solflx%TRC2P
  trcSalt_TR(idsalt_MgHPO4,L,NY,NX)=solflx%TRM1P
  trcSalt_TR(idsalt_H0PO4B,L,NY,NX)=solflx%TRH0B
  trcSalt_TR(idsalt_H3PO4B,L,NY,NX)=solflx%TRH3B
  trcSalt_TR(idsalt_FeHPO4B,L,NY,NX)=solflx%TRF1B
  trcSalt_TR(idsalt_FeH2PO4B,L,NY,NX)=solflx%TRF2B
  trcSalt_TR(idsalt_CaPO4B,L,NY,NX)=solflx%TRC0B
  trcSalt_TR(idsalt_CaHPO4B,L,NY,NX)=solflx%TRC1B
  trcSalt_TR(idsalt_CaH2PO4B,L,NY,NX)=solflx%TRC2B
  trcSalt_TR(idsalt_MgHPO4B,L,NY,NX)=solflx%TRM1B
  TRXHY(L,NY,NX)=solflx%TRXHY
  TRXAL(L,NY,NX)=solflx%TRXAL
  TRXFE(L,NY,NX)=solflx%TRXFE
  TRXCA(L,NY,NX)=solflx%TRXCA
  TRXMG(L,NY,NX)=solflx%TRXMG
  TRXNA(L,NY,NX)=solflx%TRXNA
  TRXKA(L,NY,NX)=solflx%TRXKA
  TRXHC(L,NY,NX)=solflx%TRXHC
  TRXAL2(L,NY,NX)=solflx%TRXAL2
  TRXFE2(L,NY,NX)=solflx%TRXFE2
  trcx_TR(idx_OHe,L,NY,NX)=solflx%TRXH0
  trcx_TR(idx_OHeB,L,NY,NX)=solflx%TRBH0
  trcp_TR(idsp_AlOH3,L,NY,NX)=solflx%TRALOH
  trcp_TR(idsp_FeOH3,L,NY,NX)=solflx%TRFEOH
  trcp_TR(idsp_CaCO3,L,NY,NX)=solflx%TRCACO
  trcp_TR(idsp_CaSO4,L,NY,NX)=solflx%TRCASO
  TRH2O(L,NY,NX)=solflx%TRH2O
  TBION(L,NY,NX)=solflx%TBION
  TBCO2(L,NY,NX)=solflx%TBCO2
  end subroutine GeochemAPIRecv

end module GeochemAPI
