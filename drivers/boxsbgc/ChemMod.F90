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

  !ZNH4S=ZNH4S+TR_NH4_soil*Natomw
  ystatesfl(cid_NH4_1p_conc)=ystates0l(cid_NH4_1p_conc)+solflx%TR_NH4_soil/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH4_soil)=solflx%TR_NH4_soil

  !ZNH3S=ZNH3S+TR_NH3_soil_vr*Natomw
  ystatesfl(cid_NH3_aqu_conc)=ystates0l(cid_NH3_aqu_conc)+solflx%TR_NH3_soil_vr/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH3_soil_vr)=solflx%TR_NH3_soil_vr

  !XN4  =XN4+TR_NH4_sorbed_soil
  ystatesfl(cid_XNH4_conc)=ystates0l(cid_XNH4_conc)+solflx%TR_NH4_sorbed_soil/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH4_sorbed_soil)=solflx%TR_NH4_sorbed_soil

  if(chemvar%VLWatMicPNB>0._r8)THEN
    !ZNH4B=ZNH4B+TR_NH4_band_soil*Natomw
    ystatesfl(cid_NH4_1p_band_conc)=ystates0l(cid_NH4_1p_band_conc)+solflx%TR_NH3_band_soil/chemvar%VLWatMicPNB

    !ZNH3B=ZNH3B+TR_NH3_band_soil*Natomw
    ystatesfl(cid_NH3_aqu_band_conc)=ystates0l(cid_NH3_aqu_band_conc)+solflx%TR_NH3_band_soil/chemvar%VLWatMicPNB

    !XNB  = XNB+TR_NH4_sorbed_band_soil
    ystatesfl(cid_XNH4_band_conc)=ystates0l(cid_XNH4_band_conc)+solflx%TR_NH4_sorbed_band_soil/chemvar%VLWatMicPNB

  else
    ystatesfl(cid_NH4_1p_band_conc)=ystates0l(cid_NH4_1p_band_conc)
    ystatesfl(cid_NH3_aqu_band_conc)=ystates0l(cid_NH3_aqu_band_conc)
    ystatesfl(cid_XNH4_band_conc)=ystates0l(cid_XNH4_band_conc)
  endif
  ystatesfl(fid_TR_NH3_band_soil)=solflx%TR_NH3_band_soil
  ystatesfl(fid_TR_NH3_band_soil)=solflx%TR_NH3_band_soil
  ystatesfl(fid_TR_NH4_sorbed_band_soil)=solflx%TR_NH4_sorbed_band_soil

  !H1PO4=H1PO4+TR_H1PO4_soil*Patomw
  ystatesfl(cid_H1PO4_2e_conc)=ystates0l(cid_H1PO4_2e_conc)+solflx%TR_H1PO4_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_H1PO4_soil)=solflx%TR_H1PO4_soil

  !H2PO4=H2PO4+TR_H2PO4_soil*Patomw
  ystatesfl(cid_H2PO4_1e_conc)=ystates0l(cid_H2PO4_1e_conc)+solflx%TR_H2PO4_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_H2PO4_soil)=solflx%TR_H2PO4_soil

  !XOH1 =XOH1+TR_ROH_sorbed_soil
  ystatesfl(cid_XROH1_conc)=ystates0l(cid_XROH1_conc)+solflx%TR_ROH_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_ROH_sorbed_soil)=solflx%TR_ROH_sorbed_soil

  !XOH2 =XOH2+TR_ROH2_sorbed_soil
  ystatesfl(cid_XROH2_conc)=ystates0l(cid_XROH2_conc)+solflx%TR_ROH2_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_ROH2_sorbed_soil)=solflx%TR_ROH2_sorbed_soil

  !XH1P =XH1P+TR_RHPO4_sorbed_soil
  ystatesfl(cid_XHPO4_conc)=ystates0l(cid_XHPO4_conc)+solflx%TR_RHPO4_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_RHPO4_sorbed_soil)=solflx%TR_RHPO4_sorbed_soil

  !XH2P =XH2P+TR_RH2PO4_sorbed_soil
  ystatesfl(cid_XH2PO4_conc)=ystates0l(cid_XH2PO4_conc)+solflx%TR_RH2PO4_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_RH2PO4_sorbed_soil)=solflx%TR_RH2PO4_sorbed_soil

  !PALPO=PALPO+TR_AlPO4_precip_soil
  ystatesfl(cid_Precp_AlPO4_conc)=ystates0l(cid_Precp_AlPO4_conc)+solflx%TR_AlPO4_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_AlPO4_precip_soil)=solflx%TR_AlPO4_precip_soil

  !PFEPO=PFEPO+TR_FePO4_precip_soil
  ystatesfl(cid_Precp_FePO4_conc)=ystates0l(cid_Precp_FePO4_conc)+solflx%TR_FePO4_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_FePO4_precip_soil)=solflx%TR_FePO4_precip_soil

  !PCAPD=PCAPD+TR_CaHPO4_precip_soil
  ystatesfl(cid_Precp_CaHPO4_conc)=ystates0l(cid_Precp_CaHPO4_conc)+solflx%TR_CaHPO4_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_CaHPO4_precip_soil)=solflx%TR_CaHPO4_precip_soil

  !PCAPH=PCAPH+TR_apatite_precip_soil
  ystatesfl(cid_Precp_Ca5P3O12O3H3_conc)=ystates0l(cid_Precp_Ca5P3O12O3H3_conc)+solflx%TR_apatite_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_apatite_precip_soil)=solflx%TR_apatite_precip_soil

  !PCAPM=PCAPM+TR_CaH4P2O8_precip_soil
  ystatesfl(cid_Precp_CaH4P2O8_conc)=ystates0l(cid_Precp_CaH4P2O8_conc)+solflx%TR_CaH4P2O8_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_CaH4P2O8_precip_soil)=solflx%TR_CaH4P2O8_precip_soil


  if(chemvar%VLWatMicPPB>0._r8)then
    !H1POB=H1POB+TR_H1PO4_band_soil*Patomw
    ystatesfl(cid_H1PO4_2e_band_conc)=ystates0l(cid_H1PO4_2e_band_conc)+solflx%TR_H1PO4_band_soil/chemvar%VLWatMicPPB

    !H2POB=H2POB+TR_H2PO4_band_soil*Patomw
    ystatesfl(cid_H2PO4_1e_band_conc)=ystates0l(cid_H2PO4_1e_band_conc)+solflx%TR_H2PO4_band_soil/chemvar%VLWatMicPPB

    !XOH1B=XOH1B+TR_ROH_sorbed_band_soil
    ystatesfl(cid_XROH_band_conc)=ystates0l(cid_XROH_band_conc)+solflx%TR_ROH_sorbed_band_soil/chemvar%VLWatMicPPB

    !XOH2B=XOH2B+TR_ROH2_sorbed_band_soil
    ystatesfl(cid_XROH2_band_conc)=ystates0l(cid_XROH2_band_conc)+solflx%TR_ROH2_sorbed_band_soil/chemvar%VLWatMicPPB

    !XHPO4_band_conc=XHPO4_band_conc+TR_RHPO4_sorbed_band_soil
    ystatesfl(cid_XHPO4_band_conc)=ystates0l(cid_XHPO4_band_conc)+solflx%TR_RHPO4_sorbed_band_soil/chemvar%VLWatMicPPB

    !XH2PB=XH2PB+TR_RH2PO4_sorbed_band_soil
    ystatesfl(cid_XH2PO4_band_conc)=ystates0l(cid_XH2PO4_band_conc)+solflx%TR_RH2PO4_sorbed_band_soil/chemvar%VLWatMicPPB

    !PALPB=PALPB+TR_AlPO4_precip_band_soil
    ystatesfl(cid_PrecpB_AlPO4_conc)=ystates0l(cid_PrecpB_AlPO4_conc)+solflx%TR_AlPO4_precip_band_soil/chemvar%VLWatMicPPB

    !PFEPB=PFEPB+TR_FePO4_precip_band_soil
    ystatesfl(cid_PrecpB_FePO4_con)=ystates0l(cid_PrecpB_FePO4_con)+solflx%TR_FePO4_precip_band_soil/chemvar%VLWatMicPPB

    !PCPDB=PCPDB+TR_CaHPO4_precip_band_soil
    ystatesfl(cid_PrecpB_CaHPO4_conc)=ystates0l(cid_PrecpB_CaHPO4_conc)+solflx%TR_CaHPO4_precip_band_soil/chemvar%VLWatMicPPB

    !PCPHB=PCPHB+TR_apatite_precip_band_soil
    ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)=ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)+solflx%TR_apatite_precip_band_soil/chemvar%VLWatMicPPB

    !PCPMB=PCPMB+TR_CaH4P2O8_precip_band_soil
    ystatesfl(cid_PrecpB_CaH4P2O8_conc)=ystates0l(cid_PrecpB_CaH4P2O8_conc)+solflx%TR_CaH4P2O8_precip_band_soil/chemvar%VLWatMicPPB

  else
    ystatesfl(cid_H1PO4_2e_band_conc)=ystates0l(cid_H1PO4_2e_band_conc)
    ystatesfl(cid_H2PO4_1e_band_conc)=ystates0l(cid_H2PO4_1e_band_conc)
    ystatesfl(cid_XROH_band_conc)=ystates0l(cid_XROH_band_conc)
    ystatesfl(cid_XROH2_band_conc)=ystates0l(cid_XROH2_band_conc)
    ystatesfl(cid_XHPO4_band_conc)=ystates0l(cid_XHPO4_band_conc)
    ystatesfl(cid_XH2PO4_band_conc)=ystates0l(cid_XH2PO4_band_conc)
    ystatesfl(cid_PrecpB_AlPO4_conc)=ystates0l(cid_PrecpB_AlPO4_conc)
    ystatesfl(cid_PrecpB_FePO4_con)=ystates0l(cid_PrecpB_FePO4_con)
    ystatesfl(cid_PrecpB_CaHPO4_conc)=ystates0l(cid_PrecpB_CaHPO4_conc)
    ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)=ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)
    ystatesfl(cid_PrecpB_CaH4P2O8_conc)=ystates0l(cid_PrecpB_CaH4P2O8_conc)
  endif
  ystatesfl(fid_TR_H1PO4_band_soil)=solflx%TR_H1PO4_band_soil
  ystatesfl(fid_TR_H2PO4_band_soil)=solflx%TR_H2PO4_band_soil
  ystatesfl(fid_TR_ROH_sorbed_band_soil)=solflx%TR_ROH_sorbed_band_soil
  ystatesfl(fid_TR_ROH2_sorbed_band_soil)=solflx%TR_ROH2_sorbed_band_soil
  ystatesfl(fid_TR_RHPO4_sorbed_band_soil)=solflx%TR_RHPO4_sorbed_band_soil
  ystatesfl(fid_TR_RH2PO4_sorbed_band_soil)=solflx%TR_RH2PO4_sorbed_band_soil
  ystatesfl(fid_TR_AlPO4_precip_band_soil)=solflx%TR_AlPO4_precip_band_soil
  ystatesfl(fid_TR_FePO4_precip_band_soil)=solflx%TR_FePO4_precip_band_soil
  ystatesfl(fid_TR_CaHPO4_precip_band_soil)=solflx%TR_CaHPO4_precip_band_soil
  ystatesfl(fid_TR_apatite_precip_band_soil)=solflx%TR_apatite_precip_band_soil
  ystatesfl(fid_TR_CaH4P2O8_precip_band_soil)=solflx%TR_CaH4P2O8_precip_band_soil
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

  chemvar%SoilMicPMassLayerX = forc%SoilMicPMassLayer
  chemvar%SoilMicPMassLayer  = forc%SoilMicPMassLayer
  chemvar%VLWatMicPM         = forc%VLWatMicP
  chemvar%VLWatMicPPO        = forc%VLWatMicP*forc%VLPO4
  chemvar%VLWatMicPPB        = forc%VLWatMicP*forc%VLPOB
  chemvar%VLWatMicPNH        = forc%VLWatMicP*forc%VLNH4
  chemvar%VLWatMicPNB        = forc%VLWatMicP*forc%VLNHB
  chemvar%XCEC               = forc%XCEC
  chemvar%PH                 = forc%pH
  chemvar%VLNH4              = forc%VLNH4
  chemvar%VLNHB              = forc%VLNHB

  chemvar%CAL     =forc%CAL
  chemvar%CFE     =forc%CFE
  chemvar%CCA     =forc%CCA
  chemvar%XAEC    =forc%XAEC

  chemvar%ZMG     =forc%ZMG
  chemvar%ZNA     =forc%ZNA
  chemvar%ZKA     =forc%ZKA
  chemvar%CO2S    =ystates0l(cid_CO2S)

! those below are variable

  chemvar%XROH1_conc               = ystates0l(cid_XROH1_conc)
  chemvar%XNH4_conc                = ystates0l(cid_XNH4_conc)
  chemvar%XNH4_band_conc           = ystates0l(cid_XNH4_band_conc)
  chemvar%H1PO4_2e_band_conc       = ystates0l(cid_H1PO4_2e_band_conc)  !H1POB
  chemvar%H1PO4_2e_conc            = ystates0l(cid_H1PO4_2e_conc)  !H1PO4
  chemvar%H2PO4_1e_conc            = ystates0l(cid_H2PO4_1e_conc)  !H2PO4
  chemvar%H2PO4_1e_band_conc       = ystates0l(cid_H2PO4_1e_band_conc)  !H2POB
  chemvar%XHPO4_band_conc          = ystates0l(cid_XHPO4_band_conc)
  chemvar%XH2PO4_band_conc         = ystates0l(cid_XH2PO4_band_conc)
  chemvar%XROH_band_conc           = ystates0l(cid_XROH_band_conc)
  chemvar%XHPO4_conc               = ystates0l(cid_XHPO4_conc)
  chemvar%XROH2_band_conc          = ystates0l(cid_XROH2_band_conc)
  chemvar%XH2PO4_conc              = ystates0l(cid_XH2PO4_conc)
  chemvar%XROH2_conc               = ystates0l(cid_XROH2_conc)
  chemvar%NH3_aqu_conc             = ystates0l(cid_NH3_aqu_conc)    !ZNH3S
  chemvar%NH3_aqu_band_conc        = ystates0l(cid_NH3_aqu_band_conc)    !ZNH3B
  chemvar%NH4_1p_conc              = ystates0l(cid_NH4_1p_conc)    !ZNH4S
  chemvar%NH4_1p_band_conc         = ystates0l(cid_NH4_1p_band_conc)    !ZNH4B
  chemvar%Precp_AlPO4_conc         = ystates0l(cid_Precp_AlPO4_conc)
  chemvar%PrecpB_AlPO4_conc        = ystates0l(cid_PrecpB_AlPO4_conc)
  chemvar%Precp_CaHPO4_conc        = ystates0l(cid_Precp_CaHPO4_conc)
  chemvar%PrecpB_CaHPO4_conc       = ystates0l(cid_PrecpB_CaHPO4_conc)
  chemvar%Precp_Ca5P3O12O3H3_conc  = ystates0l(cid_Precp_Ca5P3O12O3H3_conc)
  chemvar%PrecpB_Ca5P3O12O3H3_conc = ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)
  chemvar%Precp_CaH4P2O8_conc      = ystates0l(cid_Precp_CaH4P2O8_conc)
  chemvar%PrecpB_CaH4P2O8_conc     = ystates0l(cid_PrecpB_CaH4P2O8_conc)
  chemvar%Precp_FePO4_conc         = ystates0l(cid_Precp_FePO4_conc)
  chemvar%PrecpB_FePO4_con         = ystates0l(cid_PrecpB_FePO4_con)

  end subroutine SetChemVar

end module ChemMod
