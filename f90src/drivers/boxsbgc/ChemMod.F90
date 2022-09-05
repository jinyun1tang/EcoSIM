module ChemMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ChemIDMod
  use MicForcTypeMod      , only : micforctype
implicit none

  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: RunModel_nosalt
  contains
! ----------------------------------------------------------------------

  subroutine RunModel_nosalt(micfor,nvars,ystates0l, ystatesfl, err_status)
  use ChemEquilibriaMod
  use ModelStatusType     , only : model_status_type
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  implicit none
  type(micforctype), intent(in)    :: micfor
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  real(r8), intent(inout) :: ystatesfl(nvars)

  type(model_status_type), intent(out) :: err_status

  type(solute_flx_type) :: solflx
  type(chem_var_type)   :: chemvar

  integer :: jj

  call err_status%reset()

  call SetChemVar(micfor,nvars, ystates0l, chemvar)

  call NoSaltChemEquilibria(chemvar,solflx)

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
  ystatesfl(cid_CN41)=ystates0l(cid_CN41)+solflx%TRN4S/chemvar%VOLWNH
  ystatesfl(fid_TRN4S)=solflx%TRN4S

  !ZNH3S=ZNH3S+TRN3S*Natomw
  ystatesfl(cid_CN31)=ystates0l(cid_CN31)+solflx%TRN3S/chemvar%VOLWNH
  ystatesfl(fid_TRN3S)=solflx%TRN3S

  !XN4  =XN4+TRXN4
  ystatesfl(cid_XN41)=ystates0l(cid_XN41)+solflx%TRXN4/chemvar%VOLWNH
  ystatesfl(fid_TRXN4)=solflx%TRXN4

  !ZNH4B=ZNH4B+TRN4B*Natomw
  ystatesfl(cid_CN4B)=ystates0l(cid_CN4B)+solflx%TRN3B/chemvar%VOLWNB
  ystatesfl(fid_TRN3B)=solflx%TRN3B

  !ZNH3B=ZNH3B+TRN3B*Natomw
  ystatesfl(cid_CN3B)=ystates0l(cid_CN3B)+solflx%TRN3B/chemvar%VOLWNB
  ystatesfl(fid_TRN3B)=solflx%TRN3B

  !XNB  = XNB+TRXNB
  ystatesfl(cid_XN4B)=ystates0l(cid_XN4B)+solflx%TRXNB/chemvar%VOLWNB
  ystatesfl(fid_TRXNB)=solflx%TRXNB

  !H1PO4=H1PO4+TRH1P*Patomw
  ystatesfl(cid_CH1P1)=ystates0l(cid_CH1P1)+solflx%TRH1P/chemvar%VOLWPO
  ystatesfl(fid_TRH1P)=solflx%TRH1P

  !H2PO4=H2PO4+TRH2P*Patomw
  ystatesfl(cid_CH2P1)=ystates0l(cid_CH2P1)+solflx%TRH2P/chemvar%VOLWPO
  ystatesfl(fid_TRH2P)=solflx%TRH2P

  !XOH1 =XOH1+TRXH1
  ystatesfl(cid_XOH11)=ystates0l(cid_XOH11)+solflx%TRXH1/chemvar%VOLWPO
  ystatesfl(fid_TRXH1)=solflx%TRXH1

  !XOH2 =XOH2+TRXH2
  ystatesfl(cid_XOH21)=ystates0l(cid_XOH21)+solflx%TRXH2/chemvar%VOLWPO
  ystatesfl(fid_TRXH2)=solflx%TRXH2

  !XH1P =XH1P+TRX1P
  ystatesfl(cid_XH1P1)=ystates0l(cid_XH1P1)+solflx%TRX1P/chemvar%VOLWPO
  ystatesfl(fid_TRX1P)=solflx%TRX1P

  !XH2P =XH2P+TRX2P
  ystatesfl(cid_XH2P1)=ystates0l(cid_XH2P1)+solflx%TRX2P/chemvar%VOLWPO
  ystatesfl(fid_TRX2P)=solflx%TRX2P

  !PALPO=PALPO+TRALPO
  ystatesfl(cid_PALPO1)=ystates0l(cid_PALPO1)+solflx%TRALPO/chemvar%VOLWPO
  ystatesfl(fid_TRALPO)=solflx%TRALPO

  !PFEPO=PFEPO+TRFEPO
  ystatesfl(cid_PFEPO1)=ystates0l(cid_PFEPO1)+solflx%TRFEPO/chemvar%VOLWPO
  ystatesfl(fid_TRFEPO)=solflx%TRFEPO

  !PCAPD=PCAPD+TRCAPD
  ystatesfl(cid_PCAPD1)=ystates0l(cid_PCAPD1)+solflx%TRCAPD/chemvar%VOLWPO
  ystatesfl(fid_TRCAPD)=solflx%TRCAPD

  !PCAPH=PCAPH+TRCAPH
  ystatesfl(cid_PCAPH1)=ystates0l(cid_PCAPH1)+solflx%TRCAPH/chemvar%VOLWPO
  ystatesfl(fid_TRCAPH)=solflx%TRCAPH

  !PCAPM=PCAPM+TRCAPM
  ystatesfl(cid_PCAPM1)=ystates0l(cid_PCAPM1)+solflx%TRCAPM/chemvar%VOLWPO
  ystatesfl(fid_TRCAPM)=solflx%TRCAPM

  !H1POB=H1POB+TRH1B*Patomw
  ystatesfl(cid_CH1PB)=ystates0l(cid_CH1PB)+solflx%TRH1B/chemvar%VOLWPB
  ystatesfl(fid_TRH1B)=solflx%TRH1B

  !H2POB=H2POB+TRH2B*Patomw
  ystatesfl(cid_CH2PB)=ystates0l(cid_CH2PB)+solflx%TRH2B/chemvar%VOLWPB
  ystatesfl(fid_TRH2B)=solflx%TRH2B

  !XOH1B=XOH1B+TRBH1
  ystatesfl(cid_XH11B)=ystates0l(cid_XH11B)+solflx%TRBH1/chemvar%VOLWPB
  ystatesfl(fid_TRBH1)=solflx%TRBH1

  !XOH2B=XOH2B+TRBH2
  ystatesfl(cid_XH21B)=ystates0l(cid_XH21B)+solflx%TRBH2/chemvar%VOLWPB
  ystatesfl(fid_TRBH2)=solflx%TRBH2

  !XH1PB=XH1PB+TRB1P
  ystatesfl(cid_X1P1B)=ystates0l(cid_X1P1B)+solflx%TRB1P/chemvar%VOLWPB
  ystatesfl(fid_TRB1P)=solflx%TRB1P

  !XH2PB=XH2PB+TRB2P
  ystatesfl(cid_X2P1B)=ystates0l(cid_X2P1B)+solflx%TRB2P/chemvar%VOLWPB
  ystatesfl(fid_TRB2P)=solflx%TRB2P

  !PALPB=PALPB+TRALPB
  ystatesfl(cid_PALPOB)=ystates0l(cid_PALPOB)+solflx%TRALPB/chemvar%VOLWPB
  ystatesfl(fid_TRALPB)=solflx%TRALPB

  !PFEPB=PFEPB+TRFEPB
  ystatesfl(cid_PFEPOB)=ystates0l(cid_PFEPOB)+solflx%TRFEPB/chemvar%VOLWPB
  ystatesfl(fid_TRFEPB)=solflx%TRFEPB

  !PCPDB=PCPDB+TRCPDB
  ystatesfl(cid_PCAPDB)=ystates0l(cid_PCAPDB)+solflx%TRCPDB/chemvar%VOLWPB
  ystatesfl(fid_TRCPDB)=solflx%TRCPDB

  !PCPHB=PCPHB+TRCPHB
  ystatesfl(cid_PCAPHB)=ystates0l(cid_PCAPHB)+solflx%TRCPHB/chemvar%VOLWPB
  ystatesfl(fid_TRCPHB)=solflx%TRCPHB

  !PCPMB=PCPMB+TRCPMB
  ystatesfl(cid_PCAPMB)=ystates0l(cid_PCAPMB)+solflx%TRCPMB/chemvar%VOLWPB
  ystatesfl(fid_TRCPMB)=solflx%TRCPMB
  end subroutine RetrieveYstatef
! ----------------------------------------------------------------------
  subroutine SetChemVar(micfor,nvars, ystates0l, chemvar)
!
  use SoluteChemDataType, only : chem_var_type
  implicit none
  type(micforctype), intent(in)    :: micfor
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

  chemvar%BKVLX   =4.08E-002_r8
  chemvar%BKVL    =4.08E-002_r8
  chemvar%VOLWM   =1.54E-002_r8
  chemvar%VOLWPO  =1.53E-002_r8
  chemvar%VOLWPB  =5.06E-005_r8
  chemvar%VOLWNH  =1.51E-002_r8
  chemvar%VOLWNB  =3.16E-004_r8
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

end module ChemMod
