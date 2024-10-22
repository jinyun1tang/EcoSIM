module AquachemMod
  use data_kind_mod  , only : r8 => DAT_KIND_R8
  use MiniMathMod    , only : addone
  use ModelStatusType, only : model_status_type
  use EcosimConst
  use ChemIDMod
  use AquaSaltChemMod
implicit none
  private
  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: getvarlist
  public :: initmodel
  public :: getvarllen
  public :: runchem
contains
  function getvarllen(salton)result(nvars)
  implicit none
  logical, intent(in) :: salton
  integer :: nvars

  if (salton)then
    call Init_geochem_salt(nvars)
  else
    call Init_geochem_nosalt(nvars)
  endif
  end function getvarllen


! ----------------------------------------------------------------------

  subroutine Init_geochem_nosalt(nvars)

  implicit none
  integer, intent(out) :: nvars
  integer :: itemp
  itemp=0

  cid_H1PO4_2e_conc            = addone(itemp)
  cid_H1PO4_2e_band_conc       = addone(itemp)
  cid_H2PO4_1e_conc            = addone(itemp)
  cid_H2PO4_1e_band_conc       = addone(itemp)
  cid_NH3_aqu_conc             = addone(itemp)
  cid_NH3_aqu_band_conc        = addone(itemp)
  cid_NH4_1p_conc              = addone(itemp)
  cid_NH4_1p_band_conc         = addone(itemp)
  cid_XNH4_conc                = addone(itemp)
  cid_XNH4_band_conc           = addone(itemp)
  cid_XHPO4_band_conc          = addone(itemp)
  cid_XH2PO4_band_conc         = addone(itemp)
  cid_XROH_band_conc           = addone(itemp)
  cid_XHPO4_conc               = addone(itemp)
  cid_XROH2_band_conc          = addone(itemp)
  cid_XH2PO4_conc              = addone(itemp)
  cid_XROH1_conc               = addone(itemp)
  cid_XROH2_conc               = addone(itemp)
  cid_Precp_AlPO4_conc         = addone(itemp)
  cid_PrecpB_AlPO4_conc        = addone(itemp)
  cid_Precp_CaHPO4_conc        = addone(itemp)
  cid_PrecpB_CaHPO4_conc       = addone(itemp)
  cid_Precp_Ca5P3O12O3H3_conc  = addone(itemp)
  cid_PrecpB_Ca5P3O12O3H3_conc = addone(itemp)
  cid_Precp_CaH4P2O8_conc      = addone(itemp)
  cid_PrecpB_CaH4P2O8_conc     = addone(itemp)
  cid_Precp_FePO4_conc         = addone(itemp)
  cid_PrecpB_FePO4_con         = addone(itemp)

  fid_TR_NH4_soil                  = addone(itemp)
  fid_TR_NH4_band_soil             = addone(itemp)
  fid_TR_NH3_soil_vr               = addone(itemp)
  fid_TR_NH3_band_soil             = addone(itemp)
  fid_TR_H1PO4_soil                = addone(itemp)
  fid_TR_H2PO4_soil                = addone(itemp)
  fid_TR_H1PO4_band_soil           = addone(itemp)
  fid_TR_H2PO4_band_soil           = addone(itemp)
  fid_TR_NH4_sorbed_soil           = addone(itemp)
  fid_TR_NH4_sorbed_band_soil      = addone(itemp)
  fid_TR_ROH_sorbed_soil           = addone(itemp)
  fid_TR_ROH2_sorbed_soil          = addone(itemp)
  fid_TR_RHPO4_sorbed_soil         = addone(itemp)
  fid_TR_RH2PO4_sorbed_soil        = addone(itemp)
  fid_TR_ROH_sorbed_band_soil      = addone(itemp)
  fid_TR_ROH2_sorbed_band_soil     = addone(itemp)
  fid_TR_RHPO4_sorbed_band_soil    = addone(itemp)
  fid_TR_RH2PO4_sorbed_band_soil   = addone(itemp)
  fid_TR_AlPO4_precip_soil         = addone(itemp)
  fid_TR_FePO4_precip_soil         = addone(itemp)
  fid_TR_CaHPO4_precip_soil        = addone(itemp)
  fid_TR_apatite_precip_soil       = addone(itemp)
  fid_TR_CaH4P2O8_precip_soil      = addone(itemp)
  fid_TR_AlPO4_precip_band_soil    = addone(itemp)
  fid_TR_FePO4_precip_band_soil    = addone(itemp)
  fid_TR_CaHPO4_precip_band_soil   = addone(itemp)
  fid_TR_apatite_precip_band_soil  = addone(itemp)
  fid_TR_CaH4P2O8_precip_band_soil = addone(itemp)
  fid_TR_Al_3p_soil                = addone(itemp)
  nvars                            = itemp
  end subroutine Init_geochem_nosalt


!------------------------------------------------------------------
  subroutine getvarlist(nvars, varl, varlnml,unitl, vartypes, salton)

  use bhistMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  use ChemIDMod
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)
  logical, intent(in) :: salton

  if(salton)then
    call getvarlist_salt(nvars, varl, varlnml, unitl, vartypes)
  else
    call getvarlist_nosalt(nvars, varl, varlnml, unitl, vartypes)
  endif
  end subroutine getvarlist
! ----------------------------------------------------------------------

  subroutine initmodel(nvars, ystatesfl, salton, err_status)

  implicit none
  integer, intent(in) :: nvars
  logical, intent(in) :: salton
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()
  if (salton)then
    call initmodel_salt(nvars, ystatesfl,err_status)
  else
    call initmodel_nosalt(nvars, ystatesfl,err_status)
  endif

  end subroutine initmodel


! ----------------------------------------------------------------------

  subroutine initmodel_nosalt(nvars, ystatesfl,err_status)
!
! DESCRIPTION:
!
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status

  !zero out fluxes
  ystatesfl(fid_TR_NH4_soil:fid_TR_Al_3p_soil)=0._r8

!  ystatesfl(cid_CO2S)    =2.31E-002

  ystatesfl(cid_H1PO4_2e_conc)   =7.6E-3_r8
  ystatesfl(cid_H1PO4_2e_band_conc)   =7.6E-3_r8
  ystatesfl(cid_H2PO4_1e_conc)   =0.56_r8
  ystatesfl(cid_H2PO4_1e_band_conc)   =0.56_r8
  ystatesfl(cid_NH3_aqu_conc)    =1.39E-004_r8
  ystatesfl(cid_NH3_aqu_band_conc)    =1.95E-004_r8
  ystatesfl(cid_NH4_1p_conc)    =1.21_r8
  ystatesfl(cid_NH4_1p_band_conc)    =372.2_r8
  ystatesfl(cid_XNH4_conc)    =0.32_r8
  ystatesfl(cid_XNH4_band_conc)    =0.32_r8

  ystatesfl(cid_XHPO4_band_conc)   =0.0_r8
  ystatesfl(cid_XH2PO4_band_conc)   =0.0_r8
  ystatesfl(cid_XROH_band_conc)   =0.0_r8
  ystatesfl(cid_XHPO4_conc)   =0.0_r8
  ystatesfl(cid_XROH2_band_conc)   =0.0_r8
  ystatesfl(cid_XH2PO4_conc)   =0.0_r8
  ystatesfl(cid_XROH1_conc)   =0.0_r8
  ystatesfl(cid_XROH2_conc)   =0.0_r8

  ystatesfl(cid_Precp_AlPO4_conc)  =0.623_r8
  ystatesfl(cid_PrecpB_AlPO4_conc)  =0.623_r8
  ystatesfl(cid_Precp_CaHPO4_conc)  =0.0_r8
  ystatesfl(cid_PrecpB_CaHPO4_conc)  =0.0_r8
  ystatesfl(cid_Precp_Ca5P3O12O3H3_conc)  =0.0_r8
  ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)  =0.0_r8
  ystatesfl(cid_Precp_CaH4P2O8_conc)  =0.0_r8
  ystatesfl(cid_PrecpB_CaH4P2O8_conc)  =600.89_r8
  ystatesfl(cid_Precp_FePO4_conc)  =0.0_r8
  ystatesfl(cid_PrecpB_FePO4_con)  =0.0_r8

  end subroutine initmodel_nosalt
! ----------------------------------------------------------------------

  subroutine Runchem(nvars,ystates0l, ystatesfl, err_status, salton)

  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  real(r8), intent(out) :: ystatesfl(nvars)
  logical,  intent(in)  :: salton
  type(model_status_type), intent(out) :: err_status

  integer :: jj

  if(salton)then
    call RunModel_salt(nvars,ystates0l, ystatesfl, err_status)
  else
    call RunModel_nosalt(nvars,ystates0l, ystatesfl, err_status)
  endif


  end subroutine Runchem
! ----------------------------------------------------------------------

  subroutine RunModel_nosalt(nvars,ystates0l, ystatesfl, err_status)
  use ChemEquilibriaMod
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  real(r8), intent(out) :: ystatesfl(nvars)

  type(model_status_type), intent(out) :: err_status

  type(solute_flx_type) :: solflx
  type(chem_var_type)   :: chemvar

  integer :: jj

  call err_status%reset()

  call SetChemVar(nvars, ystates0l, chemvar)

  call NoSaltChemEquilibria(chemvar,solflx)

  call RetrieveYstatef(nvars, ystates0l, ystatesfl, chemvar, solflx)

  end subroutine RunModel_nosalt

! ----------------------------------------------------------------------
  subroutine SetChemVar(nvars, ystates0l, chemvar)
!
  use SoluteChemDataType, only : chem_var_type
  implicit none
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

  chemvar%SoilMicPMassLayerX   =4.08E-002_r8
  chemvar%SoilMicPMassLayer    =4.08E-002_r8
  chemvar%VLWatMicPM   =1.54E-002_r8
  chemvar%VLWatMicPPO  =1.53E-002_r8
  chemvar%VLWatMicPPB  =5.06E-005_r8
  chemvar%VLWatMicPNH  =1.51E-002_r8
  chemvar%VLWatMicPNB  =3.16E-004_r8
  chemvar%XCEC    =0.724_r8
  chemvar%PH      =5.33_r8
  chemvar%VLNH4   =0.98_r8
  chemvar%VLNHB   =0.02_r8

  chemvar%CAL     =37.04_r8
  chemvar%CFE     =0.0_r8
  chemvar%CCA     =2.5_r8
  chemvar%XAEC    =0.116_r8

  chemvar%ZMG     =0._r8
  chemvar%ZNA     =0._r8
  chemvar%ZKA     =0._r8
  chemvar%CO2S    =2.31E-002_r8

! those below are variable

  chemvar%XROH1_conc   =ystates0l(cid_XROH1_conc)
  chemvar%XNH4_conc    =ystates0l(cid_XNH4_conc)
  chemvar%XNH4_band_conc    =ystates0l(cid_XNH4_band_conc)
  chemvar%H1PO4_2e_band_conc   =ystates0l(cid_H1PO4_2e_band_conc)
  chemvar%H1PO4_2e_conc   =ystates0l(cid_H1PO4_2e_conc)
  chemvar%H2PO4_1e_conc   =ystates0l(cid_H2PO4_1e_conc)
  chemvar%H2PO4_1e_band_conc   =ystates0l(cid_H2PO4_1e_band_conc)
  chemvar%XHPO4_band_conc   =ystates0l(cid_XHPO4_band_conc)
  chemvar%XH2PO4_band_conc   =ystates0l(cid_XH2PO4_band_conc)
  chemvar%XROH_band_conc   =ystates0l(cid_XROH_band_conc)
  chemvar%XHPO4_conc   =ystates0l(cid_XHPO4_conc)
  chemvar%XROH2_band_conc   =ystates0l(cid_XROH2_band_conc)
  chemvar%XH2PO4_conc   =ystates0l(cid_XH2PO4_conc)
  chemvar%XROH2_conc   =ystates0l(cid_XROH2_conc)
  chemvar%NH3_aqu_conc    =ystates0l(cid_NH3_aqu_conc)
  chemvar%NH3_aqu_band_conc    =ystates0l(cid_NH3_aqu_band_conc)
  chemvar%NH4_1p_conc    =ystates0l(cid_NH4_1p_conc)
  chemvar%NH4_1p_band_conc    =ystates0l(cid_NH4_1p_band_conc)
  chemvar%Precp_AlPO4_conc  =ystates0l(cid_Precp_AlPO4_conc)
  chemvar%PrecpB_AlPO4_conc  =ystates0l(cid_PrecpB_AlPO4_conc)
  chemvar%Precp_CaHPO4_conc  =ystates0l(cid_Precp_CaHPO4_conc)
  chemvar%PrecpB_CaHPO4_conc  =ystates0l(cid_PrecpB_CaHPO4_conc)
  chemvar%Precp_Ca5P3O12O3H3_conc  =ystates0l(cid_Precp_Ca5P3O12O3H3_conc)
  chemvar%PrecpB_Ca5P3O12O3H3_conc  =ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)
  chemvar%Precp_CaH4P2O8_conc  =ystates0l(cid_Precp_CaH4P2O8_conc)
  chemvar%PrecpB_CaH4P2O8_conc  =ystates0l(cid_PrecpB_CaH4P2O8_conc)
  chemvar%Precp_FePO4_conc  =ystates0l(cid_Precp_FePO4_conc)
  chemvar%PrecpB_FePO4_con  =ystates0l(cid_PrecpB_FePO4_con)

  end subroutine SetChemVar

! ----------------------------------------------------------------------

  subroutine RetrieveYstatef(nvars, ystates0l, ystatesfl, chemvar, solflx)
!
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  type(chem_var_type)  , intent(in) :: chemvar
  type(solute_flx_type), intent(in) :: solflx
  real(r8), intent(out) :: ystatesfl(nvars)

  !ZNH4S=ZNH4S+TR_NH4_soil*Natomw
  ystatesfl(cid_NH4_1p_conc)=ystates0l(cid_NH4_1p_conc)+solflx%TR_NH4_soil/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH4_soil)=solflx%TR_NH4_soil

  !ZNH3S=ZNH3S+TR_NH3_soil_vr*Natomw
  ystatesfl(cid_NH3_aqu_conc)=ystates0l(cid_NH3_aqu_conc)+solflx%TR_NH3_soil_vr/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH3_soil_vr)=solflx%TR_NH3_soil_vr

  !XN4  =XN4+TR_NH4_sorbed_soil
  ystatesfl(cid_XNH4_conc)=ystates0l(cid_XNH4_conc)+solflx%TR_NH4_sorbed_soil/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH4_sorbed_soil)=solflx%TR_NH4_sorbed_soil

  !ZNH4B=ZNH4B+TR_NH4_band_soil*Natomw
  ystatesfl(cid_NH4_1p_band_conc)=ystates0l(cid_NH4_1p_band_conc)+solflx%TR_NH3_band_soil/chemvar%VLWatMicPNB
  ystatesfl(fid_TR_NH3_band_soil)=solflx%TR_NH3_band_soil

  !ZNH3B=ZNH3B+TR_NH3_band_soil*Natomw
  ystatesfl(cid_NH3_aqu_band_conc)=ystates0l(cid_NH3_aqu_band_conc)+solflx%TR_NH3_band_soil/chemvar%VLWatMicPNB
  ystatesfl(fid_TR_NH3_band_soil)=solflx%TR_NH3_band_soil

  !XNB  = XNB+TR_NH4_sorbed_band_soil
  ystatesfl(cid_XNH4_band_conc)=ystates0l(cid_XNH4_band_conc)+solflx%TR_NH4_sorbed_band_soil/chemvar%VLWatMicPNB
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

  !H1POB=H1POB+TR_H1PO4_band_soil*Patomw
  ystatesfl(cid_H1PO4_2e_band_conc)=ystates0l(cid_H1PO4_2e_band_conc)+solflx%TR_H1PO4_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_H1PO4_band_soil)=solflx%TR_H1PO4_band_soil

  !H2POB=H2POB+TR_H2PO4_band_soil*Patomw
  ystatesfl(cid_H2PO4_1e_band_conc)=ystates0l(cid_H2PO4_1e_band_conc)+solflx%TR_H2PO4_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_H2PO4_band_soil)=solflx%TR_H2PO4_band_soil

  !XOH1B=XOH1B+TR_ROH_sorbed_band_soil
  ystatesfl(cid_XROH_band_conc)=ystates0l(cid_XROH_band_conc)+solflx%TR_ROH_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_ROH_sorbed_band_soil)=solflx%TR_ROH_sorbed_band_soil

  !XOH2B=XOH2B+TR_ROH2_sorbed_band_soil
  ystatesfl(cid_XROH2_band_conc)=ystates0l(cid_XROH2_band_conc)+solflx%TR_ROH2_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_ROH2_sorbed_band_soil)=solflx%TR_ROH2_sorbed_band_soil

  !XHPO4_band_conc=XHPO4_band_conc+TR_RHPO4_sorbed_band_soil
  ystatesfl(cid_XHPO4_band_conc)=ystates0l(cid_XHPO4_band_conc)+solflx%TR_RHPO4_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_RHPO4_sorbed_band_soil)=solflx%TR_RHPO4_sorbed_band_soil

  !XH2PB=XH2PB+TR_RH2PO4_sorbed_band_soil
  ystatesfl(cid_XH2PO4_band_conc)=ystates0l(cid_XH2PO4_band_conc)+solflx%TR_RH2PO4_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_RH2PO4_sorbed_band_soil)=solflx%TR_RH2PO4_sorbed_band_soil

  !PALPB=PALPB+TR_AlPO4_precip_band_soil
  ystatesfl(cid_PrecpB_AlPO4_conc)=ystates0l(cid_PrecpB_AlPO4_conc)+solflx%TR_AlPO4_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_AlPO4_precip_band_soil)=solflx%TR_AlPO4_precip_band_soil

  !PFEPB=PFEPB+TR_FePO4_precip_band_soil
  ystatesfl(cid_PrecpB_FePO4_con)=ystates0l(cid_PrecpB_FePO4_con)+solflx%TR_FePO4_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_FePO4_precip_band_soil)=solflx%TR_FePO4_precip_band_soil

  !PCPDB=PCPDB+TR_CaHPO4_precip_band_soil
  ystatesfl(cid_PrecpB_CaHPO4_conc)=ystates0l(cid_PrecpB_CaHPO4_conc)+solflx%TR_CaHPO4_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_CaHPO4_precip_band_soil)=solflx%TR_CaHPO4_precip_band_soil

  !PCPHB=PCPHB+TR_apatite_precip_band_soil
  ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)=ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)+solflx%TR_apatite_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_apatite_precip_band_soil)=solflx%TR_apatite_precip_band_soil

  !PCPMB=PCPMB+TR_CaH4P2O8_precip_band_soil
  ystatesfl(cid_PrecpB_CaH4P2O8_conc)=ystates0l(cid_PrecpB_CaH4P2O8_conc)+solflx%TR_CaH4P2O8_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_CaH4P2O8_precip_band_soil)=solflx%TR_CaH4P2O8_precip_band_soil
  end subroutine RetrieveYstatef

end module AquachemMod
