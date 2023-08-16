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

  cid_CH1P1  =addone(itemp)
  cid_CH1PB  =addone(itemp)
  cid_CH2P1  =addone(itemp)
  cid_CH2PB  =addone(itemp)
  cid_CN31   =addone(itemp)
  cid_CN3B   =addone(itemp)
  cid_CN41   =addone(itemp)
  cid_CN4B   =addone(itemp)
  cid_XN41   =addone(itemp)
  cid_XN4B   =addone(itemp)
  cid_X1P1B  =addone(itemp)
  cid_X2P1B  =addone(itemp)
  cid_XH11B  =addone(itemp)
  cid_XH1P1  =addone(itemp)
  cid_XH21B  =addone(itemp)
  cid_XH2P1  =addone(itemp)
  cid_XOH11  =addone(itemp)
  cid_XOH21  =addone(itemp)
  cid_PALPO1 =addone(itemp)
  cid_PALPOB =addone(itemp)
  cid_PCAPD1 =addone(itemp)
  cid_PCAPDB =addone(itemp)
  cid_PCAPH1 =addone(itemp)
  cid_PCAPHB =addone(itemp)
  cid_PCAPM1 =addone(itemp)
  cid_PCAPMB =addone(itemp)
  cid_PFEPO1 =addone(itemp)
  cid_PFEPOB =addone(itemp)

  fid_TRN4S = addone(itemp)
  fid_TRN4B = addone(itemp)
  fid_TRN3S = addone(itemp)
  fid_TRN3B = addone(itemp)
  fid_TRH1P = addone(itemp)
  fid_TRH2P = addone(itemp)
  fid_TRH1B = addone(itemp)
  fid_TRH2B = addone(itemp)
  fid_TRXN4 = addone(itemp)
  fid_TRXNB = addone(itemp)
  fid_TRXH1 = addone(itemp)
  fid_TRXH2 = addone(itemp)
  fid_TRX1P = addone(itemp)
  fid_TRX2P = addone(itemp)
  fid_TRBH1 = addone(itemp)
  fid_TRBH2 = addone(itemp)
  fid_TRB1P = addone(itemp)
  fid_TRB2P = addone(itemp)
  fid_TRALPO= addone(itemp)
  fid_TRFEPO= addone(itemp)
  fid_TRCAPD= addone(itemp)
  fid_TRCAPH= addone(itemp)
  fid_TRCAPM= addone(itemp)
  fid_TRALPB= addone(itemp)
  fid_TRFEPB= addone(itemp)
  fid_TRCPDB= addone(itemp)
  fid_TRCPHB= addone(itemp)
  fid_TRCPMB= addone(itemp)
  fid_TRAL  = addone(itemp)
  nvars = itemp
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
  ystatesfl(fid_TRN4S:fid_TRAL)=0._r8

!  ystatesfl(cid_CO2S)    =2.31E-002

  ystatesfl(cid_CH1P1)   =7.6E-3_r8
  ystatesfl(cid_CH1PB)   =7.6E-3_r8
  ystatesfl(cid_CH2P1)   =0.56_r8
  ystatesfl(cid_CH2PB)   =0.56_r8
  ystatesfl(cid_CN31)    =1.39E-004_r8
  ystatesfl(cid_CN3B)    =1.95E-004_r8
  ystatesfl(cid_CN41)    =1.21_r8
  ystatesfl(cid_CN4B)    =372.2_r8
  ystatesfl(cid_XN41)    =0.32_r8
  ystatesfl(cid_XN4B)    =0.32_r8

  ystatesfl(cid_X1P1B)   =0.0_r8
  ystatesfl(cid_X2P1B)   =0.0_r8
  ystatesfl(cid_XH11B)   =0.0_r8
  ystatesfl(cid_XH1P1)   =0.0_r8
  ystatesfl(cid_XH21B)   =0.0_r8
  ystatesfl(cid_XH2P1)   =0.0_r8
  ystatesfl(cid_XOH11)   =0.0_r8
  ystatesfl(cid_XOH21)   =0.0_r8

  ystatesfl(cid_PALPO1)  =0.623_r8
  ystatesfl(cid_PALPOB)  =0.623_r8
  ystatesfl(cid_PCAPD1)  =0.0_r8
  ystatesfl(cid_PCAPDB)  =0.0_r8
  ystatesfl(cid_PCAPH1)  =0.0_r8
  ystatesfl(cid_PCAPHB)  =0.0_r8
  ystatesfl(cid_PCAPM1)  =0.0_r8
  ystatesfl(cid_PCAPMB)  =600.89_r8
  ystatesfl(cid_PFEPO1)  =0.0_r8
  ystatesfl(cid_PFEPOB)  =0.0_r8

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

  chemvar%XOH11   =ystates0l(cid_XOH11)
  chemvar%XN41    =ystates0l(cid_XN41)
  chemvar%XN4B    =ystates0l(cid_XN4B)
  chemvar%CH1PB   =ystates0l(cid_CH1PB)
  chemvar%CH1P1   =ystates0l(cid_CH1P1)
  chemvar%CH2P1   =ystates0l(cid_CH2P1)
  chemvar%CH2PB   =ystates0l(cid_CH2PB)
  chemvar%X1P1B   =ystates0l(cid_X1P1B)
  chemvar%X2P1B   =ystates0l(cid_X2P1B)
  chemvar%XH11B   =ystates0l(cid_XH11B)
  chemvar%XH1P1   =ystates0l(cid_XH1P1)
  chemvar%XH21B   =ystates0l(cid_XH21B)
  chemvar%XH2P1   =ystates0l(cid_XH2P1)
  chemvar%XOH21   =ystates0l(cid_XOH21)
  chemvar%CN31    =ystates0l(cid_CN31)
  chemvar%CN3B    =ystates0l(cid_CN3B)
  chemvar%CN41    =ystates0l(cid_CN41)
  chemvar%CN4B    =ystates0l(cid_CN4B)
  chemvar%PALPO1  =ystates0l(cid_PALPO1)
  chemvar%PALPOB  =ystates0l(cid_PALPOB)
  chemvar%PCAPD1  =ystates0l(cid_PCAPD1)
  chemvar%PCAPDB  =ystates0l(cid_PCAPDB)
  chemvar%PCAPH1  =ystates0l(cid_PCAPH1)
  chemvar%PCAPHB  =ystates0l(cid_PCAPHB)
  chemvar%PCAPM1  =ystates0l(cid_PCAPM1)
  chemvar%PCAPMB  =ystates0l(cid_PCAPMB)
  chemvar%PFEPO1  =ystates0l(cid_PFEPO1)
  chemvar%PFEPOB  =ystates0l(cid_PFEPOB)

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

  !ZNH4S=ZNH4S+TRN4S*Natomw
  ystatesfl(cid_CN41)=ystates0l(cid_CN41)+solflx%TRN4S/chemvar%VLWatMicPNH
  ystatesfl(fid_TRN4S)=solflx%TRN4S

  !ZNH3S=ZNH3S+TRN3S*Natomw
  ystatesfl(cid_CN31)=ystates0l(cid_CN31)+solflx%TRN3S/chemvar%VLWatMicPNH
  ystatesfl(fid_TRN3S)=solflx%TRN3S

  !XN4  =XN4+TRXN4
  ystatesfl(cid_XN41)=ystates0l(cid_XN41)+solflx%TRXN4/chemvar%VLWatMicPNH
  ystatesfl(fid_TRXN4)=solflx%TRXN4

  !ZNH4B=ZNH4B+TRN4B*Natomw
  ystatesfl(cid_CN4B)=ystates0l(cid_CN4B)+solflx%TRN3B/chemvar%VLWatMicPNB
  ystatesfl(fid_TRN3B)=solflx%TRN3B

  !ZNH3B=ZNH3B+TRN3B*Natomw
  ystatesfl(cid_CN3B)=ystates0l(cid_CN3B)+solflx%TRN3B/chemvar%VLWatMicPNB
  ystatesfl(fid_TRN3B)=solflx%TRN3B

  !XNB  = XNB+TRXNB
  ystatesfl(cid_XN4B)=ystates0l(cid_XN4B)+solflx%TRXNB/chemvar%VLWatMicPNB
  ystatesfl(fid_TRXNB)=solflx%TRXNB

  !H1PO4=H1PO4+TRH1P*Patomw
  ystatesfl(cid_CH1P1)=ystates0l(cid_CH1P1)+solflx%TRH1P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRH1P)=solflx%TRH1P

  !H2PO4=H2PO4+TRH2P*Patomw
  ystatesfl(cid_CH2P1)=ystates0l(cid_CH2P1)+solflx%TRH2P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRH2P)=solflx%TRH2P

  !XOH1 =XOH1+TRXH1
  ystatesfl(cid_XOH11)=ystates0l(cid_XOH11)+solflx%TRXH1/chemvar%VLWatMicPPO
  ystatesfl(fid_TRXH1)=solflx%TRXH1

  !XOH2 =XOH2+TRXH2
  ystatesfl(cid_XOH21)=ystates0l(cid_XOH21)+solflx%TRXH2/chemvar%VLWatMicPPO
  ystatesfl(fid_TRXH2)=solflx%TRXH2

  !XH1P =XH1P+TRX1P
  ystatesfl(cid_XH1P1)=ystates0l(cid_XH1P1)+solflx%TRX1P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRX1P)=solflx%TRX1P

  !XH2P =XH2P+TRX2P
  ystatesfl(cid_XH2P1)=ystates0l(cid_XH2P1)+solflx%TRX2P/chemvar%VLWatMicPPO
  ystatesfl(fid_TRX2P)=solflx%TRX2P

  !PALPO=PALPO+TRALPO
  ystatesfl(cid_PALPO1)=ystates0l(cid_PALPO1)+solflx%TRALPO/chemvar%VLWatMicPPO
  ystatesfl(fid_TRALPO)=solflx%TRALPO

  !PFEPO=PFEPO+TRFEPO
  ystatesfl(cid_PFEPO1)=ystates0l(cid_PFEPO1)+solflx%TRFEPO/chemvar%VLWatMicPPO
  ystatesfl(fid_TRFEPO)=solflx%TRFEPO

  !PCAPD=PCAPD+TRCAPD
  ystatesfl(cid_PCAPD1)=ystates0l(cid_PCAPD1)+solflx%TRCAPD/chemvar%VLWatMicPPO
  ystatesfl(fid_TRCAPD)=solflx%TRCAPD

  !PCAPH=PCAPH+TRCAPH
  ystatesfl(cid_PCAPH1)=ystates0l(cid_PCAPH1)+solflx%TRCAPH/chemvar%VLWatMicPPO
  ystatesfl(fid_TRCAPH)=solflx%TRCAPH

  !PCAPM=PCAPM+TRCAPM
  ystatesfl(cid_PCAPM1)=ystates0l(cid_PCAPM1)+solflx%TRCAPM/chemvar%VLWatMicPPO
  ystatesfl(fid_TRCAPM)=solflx%TRCAPM

  !H1POB=H1POB+TRH1B*Patomw
  ystatesfl(cid_CH1PB)=ystates0l(cid_CH1PB)+solflx%TRH1B/chemvar%VLWatMicPPB
  ystatesfl(fid_TRH1B)=solflx%TRH1B

  !H2POB=H2POB+TRH2B*Patomw
  ystatesfl(cid_CH2PB)=ystates0l(cid_CH2PB)+solflx%TRH2B/chemvar%VLWatMicPPB
  ystatesfl(fid_TRH2B)=solflx%TRH2B

  !XOH1B=XOH1B+TRBH1
  ystatesfl(cid_XH11B)=ystates0l(cid_XH11B)+solflx%TRBH1/chemvar%VLWatMicPPB
  ystatesfl(fid_TRBH1)=solflx%TRBH1

  !XOH2B=XOH2B+TRBH2
  ystatesfl(cid_XH21B)=ystates0l(cid_XH21B)+solflx%TRBH2/chemvar%VLWatMicPPB
  ystatesfl(fid_TRBH2)=solflx%TRBH2

  !XH1PB=XH1PB+TRB1P
  ystatesfl(cid_X1P1B)=ystates0l(cid_X1P1B)+solflx%TRB1P/chemvar%VLWatMicPPB
  ystatesfl(fid_TRB1P)=solflx%TRB1P

  !XH2PB=XH2PB+TRB2P
  ystatesfl(cid_X2P1B)=ystates0l(cid_X2P1B)+solflx%TRB2P/chemvar%VLWatMicPPB
  ystatesfl(fid_TRB2P)=solflx%TRB2P

  !PALPB=PALPB+TRALPB
  ystatesfl(cid_PALPOB)=ystates0l(cid_PALPOB)+solflx%TRALPB/chemvar%VLWatMicPPB
  ystatesfl(fid_TRALPB)=solflx%TRALPB

  !PFEPB=PFEPB+TRFEPB
  ystatesfl(cid_PFEPOB)=ystates0l(cid_PFEPOB)+solflx%TRFEPB/chemvar%VLWatMicPPB
  ystatesfl(fid_TRFEPB)=solflx%TRFEPB

  !PCPDB=PCPDB+TRCPDB
  ystatesfl(cid_PCAPDB)=ystates0l(cid_PCAPDB)+solflx%TRCPDB/chemvar%VLWatMicPPB
  ystatesfl(fid_TRCPDB)=solflx%TRCPDB

  !PCPHB=PCPHB+TRCPHB
  ystatesfl(cid_PCAPHB)=ystates0l(cid_PCAPHB)+solflx%TRCPHB/chemvar%VLWatMicPPB
  ystatesfl(fid_TRCPHB)=solflx%TRCPHB

  !PCPMB=PCPMB+TRCPMB
  ystatesfl(cid_PCAPMB)=ystates0l(cid_PCAPMB)+solflx%TRCPMB/chemvar%VLWatMicPPB
  ystatesfl(fid_TRCPMB)=solflx%TRCPMB
  end subroutine RetrieveYstatef

end module AquachemMod
