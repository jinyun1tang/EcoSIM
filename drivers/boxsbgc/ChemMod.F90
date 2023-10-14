module ChemMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ChemIDMod
  use MicForcTypeMod      , only : micforctype
implicit none

  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: RunModel_nosalt
  contains
! ----------------------------------------------------------------------

  subroutine RunModel_nosalt(forc,micfor,nvars,ystates0l, ystatesfl, err_status)
  use ChemEquilibriaMod
  use ModelStatusType     , only : model_status_type
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  use ForcTypeMod         , only : forc_type
  implicit none
  type(forc_type), intent(in) :: forc
  type(micforctype), intent(in)    :: micfor
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  real(r8), intent(inout) :: ystatesfl(nvars)

  type(model_status_type), intent(out) :: err_status

  type(solute_flx_type) :: solflx
  type(chem_var_type)   :: chemvar

  integer :: jj

  call err_status%reset()
!  print*,'SetChemVar'
  call SetChemVar(forc,micfor,nvars, ystates0l, chemvar)
!  print*,'NoSaltChemEquilibria'
  call NoSaltChemEquilibria(chemvar,solflx)
!  print*,'RetrieveYstatef'
  call RetrieveYstatef(nvars, ystates0l, ystatesfl, chemvar, solflx)

  end subroutine RunModel_nosalt


! ----------------------------------------------------------------------

  subroutine RetrieveYstatef(nvars, ystates0l, ystatesfl, chemvar, solflx)
!
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  type(chem_var_type)  , intent(in) :: chemvar
  type(solute_flx_type), intent(in) :: solflx
  real(r8), intent(inout) :: ystatesfl(nvars)

  !ZNH4S=ZNH4S+TRN4S*Natomw
  ystatesfl(cid_NH4_1p_conc)=ystates0l(cid_NH4_1p_conc)+solflx%TRN4S/chemvar%VLWatMicPNH
  ystatesfl(fid_TRN4S)=solflx%TRN4S

  !ZNH3S=ZNH3S+TRN3S*Natomw
  ystatesfl(cid_NH3_aqu_conc)=ystates0l(cid_NH3_aqu_conc)+solflx%TRN3S/chemvar%VLWatMicPNH
  ystatesfl(fid_TRN3S)=solflx%TRN3S

  !XN4  =XN4+TRXN4
  ystatesfl(cid_XNH4_conc)=ystates0l(cid_XNH4_conc)+solflx%TRXN4/chemvar%VLWatMicPNH
  ystatesfl(fid_TRXN4)=solflx%TRXN4

  if(chemvar%VLWatMicPNB>0._r8)THEN
    !ZNH4B=ZNH4B+TRN4B*Natomw
    ystatesfl(cid_NH4_1p_Bconc)=ystates0l(cid_NH4_1p_Bconc)+solflx%TRN3B/chemvar%VLWatMicPNB

    !ZNH3B=ZNH3B+TRN3B*Natomw
    ystatesfl(cid_NH3_aqu_Bconc)=ystates0l(cid_NH3_aqu_Bconc)+solflx%TRN3B/chemvar%VLWatMicPNB

    !XNB  = XNB+TRXNB
    ystatesfl(cid_XNH4_Bconc)=ystates0l(cid_XNH4_Bconc)+solflx%TRXNB/chemvar%VLWatMicPNB

  else
    ystatesfl(cid_NH4_1p_Bconc)=ystates0l(cid_NH4_1p_Bconc)
    ystatesfl(cid_NH3_aqu_Bconc)=ystates0l(cid_NH3_aqu_Bconc)
    ystatesfl(cid_XNH4_Bconc)=ystates0l(cid_XNH4_Bconc)
  endif
  ystatesfl(fid_TRN3B)=solflx%TRN3B
  ystatesfl(fid_TRN3B)=solflx%TRN3B
  ystatesfl(fid_TRXNB)=solflx%TRXNB

  !H1PO4=H1PO4+TRH1P*Patomw
  ystatesfl(cid_H1PO4_2e_conc)=ystates0l(cid_H1PO4_2e_conc)+solflx%TRH1P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRH1P)=solflx%TRH1P

  !H2PO4=H2PO4+TRH2P*Patomw
  ystatesfl(cid_H2PO4_1e_conc)=ystates0l(cid_H2PO4_1e_conc)+solflx%TRH2P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRH2P)=solflx%TRH2P

  !XOH1 =XOH1+TRXH1
  ystatesfl(cid_XROH1_conc)=ystates0l(cid_XROH1_conc)+solflx%TRXH1/chemvar%VLWatMicPPO
  ystatesfl(fid_TRXH1)=solflx%TRXH1

  !XOH2 =XOH2+TRXH2
  ystatesfl(cid_XROH2_conc)=ystates0l(cid_XROH2_conc)+solflx%TRXH2/chemvar%VLWatMicPPO
  ystatesfl(fid_TRXH2)=solflx%TRXH2

  !XH1P =XH1P+TRX1P
  ystatesfl(cid_XHPO4_conc)=ystates0l(cid_XHPO4_conc)+solflx%TRX1P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRX1P)=solflx%TRX1P

  !XH2P =XH2P+TRX2P
  ystatesfl(cid_XH2PO4_conc)=ystates0l(cid_XH2PO4_conc)+solflx%TRX2P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRX2P)=solflx%TRX2P

  !PALPO=PALPO+TR_AlPO4
  ystatesfl(cid_Precp_AlPO4_conc)=ystates0l(cid_Precp_AlPO4_conc)+solflx%TR_AlPO4/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_AlPO4)=solflx%TR_AlPO4

  !PFEPO=PFEPO+TRFEPO
  ystatesfl(cid_Precp_FePO4_conc)=ystates0l(cid_Precp_FePO4_conc)+solflx%TRFEPO/chemvar%VLWatMicPPO
  ystatesfl(fid_TRFEPO)=solflx%TRFEPO

  !PCAPD=PCAPD+TRCAPD
  ystatesfl(cid_Precp_CaHPO4_conc)=ystates0l(cid_Precp_CaHPO4_conc)+solflx%TRCAPD/chemvar%VLWatMicPPO
  ystatesfl(fid_TRCAPD)=solflx%TRCAPD

  !PCAPH=PCAPH+TRCAPH
  ystatesfl(cid_Precp_Ca5P3O12O3H3_conc)=ystates0l(cid_Precp_Ca5P3O12O3H3_conc)+solflx%TRCAPH/chemvar%VLWatMicPPO
  ystatesfl(fid_TRCAPH)=solflx%TRCAPH

  !PCAPM=PCAPM+TRCAPM
  ystatesfl(cid_Precp_CaH4P2O8_conc)=ystates0l(cid_Precp_CaH4P2O8_conc)+solflx%TRCAPM/chemvar%VLWatMicPPO
  ystatesfl(fid_TRCAPM)=solflx%TRCAPM


  if(chemvar%VLWatMicPPB>0._r8)then
    !H1POB=H1POB+TRH1B*Patomw
    ystatesfl(cid_H1PO4_2e_Bconc)=ystates0l(cid_H1PO4_2e_Bconc)+solflx%TRH1B/chemvar%VLWatMicPPB

    !H2POB=H2POB+TRH2B*Patomw
    ystatesfl(cid_H2PO4_1e_Bconc)=ystates0l(cid_H2PO4_1e_Bconc)+solflx%TRH2B/chemvar%VLWatMicPPB

    !XOH1B=XOH1B+TRBH1
    ystatesfl(cid_XROH_Bconc)=ystates0l(cid_XROH_Bconc)+solflx%TRBH1/chemvar%VLWatMicPPB

    !XOH2B=XOH2B+TRBH2
    ystatesfl(cid_XROH2_Bconc)=ystates0l(cid_XROH2_Bconc)+solflx%TRBH2/chemvar%VLWatMicPPB

    !XHPO4_Bconc=XHPO4_Bconc+TRB1P
    ystatesfl(cid_XHPO4_Bconc)=ystates0l(cid_XHPO4_Bconc)+solflx%TRB1P/chemvar%VLWatMicPPB

    !XH2PB=XH2PB+TRB2P
    ystatesfl(cid_XH2PO4_Bconc)=ystates0l(cid_XH2PO4_Bconc)+solflx%TRB2P/chemvar%VLWatMicPPB

    !PALPB=PALPB+TRALPB
    ystatesfl(cid_PrecpB_AlPO4_conc)=ystates0l(cid_PrecpB_AlPO4_conc)+solflx%TRALPB/chemvar%VLWatMicPPB

    !PFEPB=PFEPB+TRFEPB
    ystatesfl(cid_PrecpB_FePO4_con)=ystates0l(cid_PrecpB_FePO4_con)+solflx%TRFEPB/chemvar%VLWatMicPPB

    !PCPDB=PCPDB+TRCPDB
    ystatesfl(cid_PrecpB_CaHPO4_conc)=ystates0l(cid_PrecpB_CaHPO4_conc)+solflx%TRCPDB/chemvar%VLWatMicPPB

    !PCPHB=PCPHB+TRCPHB
    ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)=ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)+solflx%TRCPHB/chemvar%VLWatMicPPB

    !PCPMB=PCPMB+TRCPMB
    ystatesfl(cid_PrecpB_CaH2PO4_con)=ystates0l(cid_PrecpB_CaH2PO4_con)+solflx%TRCPMB/chemvar%VLWatMicPPB

  else
    ystatesfl(cid_H1PO4_2e_Bconc)=ystates0l(cid_H1PO4_2e_Bconc)
    ystatesfl(cid_H2PO4_1e_Bconc)=ystates0l(cid_H2PO4_1e_Bconc)
    ystatesfl(cid_XROH_Bconc)=ystates0l(cid_XROH_Bconc)
    ystatesfl(cid_XROH2_Bconc)=ystates0l(cid_XROH2_Bconc)
    ystatesfl(cid_XHPO4_Bconc)=ystates0l(cid_XHPO4_Bconc)
    ystatesfl(cid_XH2PO4_Bconc)=ystates0l(cid_XH2PO4_Bconc)
    ystatesfl(cid_PrecpB_AlPO4_conc)=ystates0l(cid_PrecpB_AlPO4_conc)
    ystatesfl(cid_PrecpB_FePO4_con)=ystates0l(cid_PrecpB_FePO4_con)
    ystatesfl(cid_PrecpB_CaHPO4_conc)=ystates0l(cid_PrecpB_CaHPO4_conc)
    ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)=ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)
    ystatesfl(cid_PrecpB_CaH2PO4_con)=ystates0l(cid_PrecpB_CaH2PO4_con)
  endif
  ystatesfl(fid_TRH1B)=solflx%TRH1B
  ystatesfl(fid_TRH2B)=solflx%TRH2B
  ystatesfl(fid_TRBH1)=solflx%TRBH1
  ystatesfl(fid_TRBH2)=solflx%TRBH2
  ystatesfl(fid_TRB1P)=solflx%TRB1P
  ystatesfl(fid_TRB2P)=solflx%TRB2P
  ystatesfl(fid_TRALPB)=solflx%TRALPB
  ystatesfl(fid_TRFEPB)=solflx%TRFEPB
  ystatesfl(fid_TRCPDB)=solflx%TRCPDB
  ystatesfl(fid_TRCPHB)=solflx%TRCPHB
  ystatesfl(fid_TRCPMB)=solflx%TRCPMB
  end subroutine RetrieveYstatef
! ----------------------------------------------------------------------
  subroutine SetChemVar(forc,micfor,nvars, ystates0l, chemvar)
!
  use SoluteChemDataType, only : chem_var_type
  use ForcTypeMod         , only : forc_type
  implicit none
  type(micforctype), intent(in)    :: micfor
  type(forc_type), intent(in) :: forc
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  type(chem_var_type), intent(out)   :: chemvar

! Gapon selectivity coefficient
  real(r8), parameter :: GKCA=0.25_r8   !for Ca-Al
  real(r8), parameter :: GKCH=0.25_r8   !for Ca-H
  real(r8), parameter :: GKCM=0.60_r8   !for Ca-Mg
  real(r8), parameter :: GKCK=3.0_r8    !for Ca-K
  real(r8), parameter :: GKCN=0.16_r8   !for Ca-Na
  real(r8), parameter :: GKC4=2.5e-2_r8 !for Ca-NH4

  chemvar%GKC4 = GKC4
  chemvar%GKCA = GKCA
  chemvar%GKCH = GKCH
  chemvar%GKCM = GKCM
  chemvar%GKCK = GKCK
  chemvar%GKCN = GKCN

  chemvar%SoilMicPMassLayerX   =forc%SoilMicPMassLayer
  chemvar%SoilMicPMassLayer    =forc%SoilMicPMassLayer
  chemvar%VLWatMicPM   =forc%VLWatMicP
  chemvar%VLWatMicPPO  =forc%VLWatMicP*forc%VLPO4
  chemvar%VLWatMicPPB  =forc%VLWatMicP*forc%VLPOB
  chemvar%VLWatMicPNH  =forc%VLWatMicP*forc%VLNH4
  chemvar%VLWatMicPNB  =forc%VLWatMicP*forc%VLNHB
  chemvar%XCEC    =forc%XCEC
  chemvar%PH      =forc%pH
  chemvar%VLNH4   =forc%VLNH4
  chemvar%VLNHB   =forc%VLNHB

  chemvar%CAL     =forc%CAL
  chemvar%CFE     =forc%CFE
  chemvar%CCA     =forc%CCA
  chemvar%XAEC    =forc%XAEC

  chemvar%ZMG     =forc%ZMG
  chemvar%ZNA     =forc%ZNA
  chemvar%ZKA     =forc%ZKA
  chemvar%CO2S    =ystates0l(cid_CO2S)

! those below are variable

  chemvar%XROH1_conc   =ystates0l(cid_XROH1_conc)
  chemvar%XNH4_conc    =ystates0l(cid_XNH4_conc)
  chemvar%XNH4_Bconc    =ystates0l(cid_XNH4_Bconc)
  chemvar%H1PO4_2e_Bconc   =ystates0l(cid_H1PO4_2e_Bconc)  !H1POB
  chemvar%H1PO4_2e_conc   =ystates0l(cid_H1PO4_2e_conc)  !H1PO4
  chemvar%H2PO4_1e_conc   =ystates0l(cid_H2PO4_1e_conc)  !H2PO4
  chemvar%H2PO4_1e_Bconc   =ystates0l(cid_H2PO4_1e_Bconc)  !H2POB
  chemvar%XHPO4_Bconc   =ystates0l(cid_XHPO4_Bconc)
  chemvar%XH2PO4_Bconc   =ystates0l(cid_XH2PO4_Bconc)
  chemvar%XROH_Bconc   =ystates0l(cid_XROH_Bconc)
  chemvar%XHPO4_conc   =ystates0l(cid_XHPO4_conc)
  chemvar%XROH2_Bconc   =ystates0l(cid_XROH2_Bconc)
  chemvar%XH2PO4_conc   =ystates0l(cid_XH2PO4_conc)
  chemvar%XROH2_conc   =ystates0l(cid_XROH2_conc)
  chemvar%NH3_aqu_conc    =ystates0l(cid_NH3_aqu_conc)    !ZNH3S
  chemvar%NH3_aqu_Bconc    =ystates0l(cid_NH3_aqu_Bconc)    !ZNH3B
  chemvar%NH4_1p_conc    =ystates0l(cid_NH4_1p_conc)    !ZNH4S
  chemvar%NH4_1p_Bconc    =ystates0l(cid_NH4_1p_Bconc)    !ZNH4B
  chemvar%Precp_AlPO4_conc  =ystates0l(cid_Precp_AlPO4_conc)
  chemvar%PrecpB_AlPO4_conc  =ystates0l(cid_PrecpB_AlPO4_conc)
  chemvar%Precp_CaHPO4_conc  =ystates0l(cid_Precp_CaHPO4_conc)
  chemvar%PrecpB_CaHPO4_conc  =ystates0l(cid_PrecpB_CaHPO4_conc)
  chemvar%Precp_Ca5P3O12O3H3_conc  =ystates0l(cid_Precp_Ca5P3O12O3H3_conc)
  chemvar%PrecpB_Ca5P3O12O3H3_conc  =ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)
  chemvar%Precp_CaH4P2O8_conc  =ystates0l(cid_Precp_CaH4P2O8_conc)
  chemvar%PrecpB_CaH2PO4_con  =ystates0l(cid_PrecpB_CaH2PO4_con)
  chemvar%Precp_FePO4_conc  =ystates0l(cid_Precp_FePO4_conc)
  chemvar%PrecpB_FePO4_con  =ystates0l(cid_PrecpB_FePO4_con)

  end subroutine SetChemVar

end module ChemMod
