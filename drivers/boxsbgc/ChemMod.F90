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
  ystatesfl(cid_CN41)=ystates0l(cid_CN41)+solflx%TRN4S/chemvar%VLWatMicPNH
  ystatesfl(fid_TRN4S)=solflx%TRN4S

  !ZNH3S=ZNH3S+TRN3S*Natomw
  ystatesfl(cid_CN31)=ystates0l(cid_CN31)+solflx%TRN3S/chemvar%VLWatMicPNH
  ystatesfl(fid_TRN3S)=solflx%TRN3S

  !XN4  =XN4+TRXN4
  ystatesfl(cid_XN41)=ystates0l(cid_XN41)+solflx%TRXN4/chemvar%VLWatMicPNH
  ystatesfl(fid_TRXN4)=solflx%TRXN4

  if(chemvar%VLWatMicPNB>0._r8)THEN
    !ZNH4B=ZNH4B+TRN4B*Natomw
    ystatesfl(cid_CN4B)=ystates0l(cid_CN4B)+solflx%TRN3B/chemvar%VLWatMicPNB

    !ZNH3B=ZNH3B+TRN3B*Natomw
    ystatesfl(cid_CN3B)=ystates0l(cid_CN3B)+solflx%TRN3B/chemvar%VLWatMicPNB

    !XNB  = XNB+TRXNB
    ystatesfl(cid_XN4B)=ystates0l(cid_XN4B)+solflx%TRXNB/chemvar%VLWatMicPNB

  else
    ystatesfl(cid_CN4B)=ystates0l(cid_CN4B)
    ystatesfl(cid_CN3B)=ystates0l(cid_CN3B)
    ystatesfl(cid_XN4B)=ystates0l(cid_XN4B)
  endif
  ystatesfl(fid_TRN3B)=solflx%TRN3B
  ystatesfl(fid_TRN3B)=solflx%TRN3B
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


  if(chemvar%VLWatMicPPB>0._r8)then
    !H1POB=H1POB+TRH1B*Patomw
    ystatesfl(cid_CH1PB)=ystates0l(cid_CH1PB)+solflx%TRH1B/chemvar%VLWatMicPPB

    !H2POB=H2POB+TRH2B*Patomw
    ystatesfl(cid_CH2PB)=ystates0l(cid_CH2PB)+solflx%TRH2B/chemvar%VLWatMicPPB

    !XOH1B=XOH1B+TRBH1
    ystatesfl(cid_XH11B)=ystates0l(cid_XH11B)+solflx%TRBH1/chemvar%VLWatMicPPB

    !XOH2B=XOH2B+TRBH2
    ystatesfl(cid_XH21B)=ystates0l(cid_XH21B)+solflx%TRBH2/chemvar%VLWatMicPPB

    !XH1PB=XH1PB+TRB1P
    ystatesfl(cid_X1P1B)=ystates0l(cid_X1P1B)+solflx%TRB1P/chemvar%VLWatMicPPB

    !XH2PB=XH2PB+TRB2P
    ystatesfl(cid_X2P1B)=ystates0l(cid_X2P1B)+solflx%TRB2P/chemvar%VLWatMicPPB

    !PALPB=PALPB+TRALPB
    ystatesfl(cid_PALPOB)=ystates0l(cid_PALPOB)+solflx%TRALPB/chemvar%VLWatMicPPB

    !PFEPB=PFEPB+TRFEPB
    ystatesfl(cid_PFEPOB)=ystates0l(cid_PFEPOB)+solflx%TRFEPB/chemvar%VLWatMicPPB

    !PCPDB=PCPDB+TRCPDB
    ystatesfl(cid_PCAPDB)=ystates0l(cid_PCAPDB)+solflx%TRCPDB/chemvar%VLWatMicPPB

    !PCPHB=PCPHB+TRCPHB
    ystatesfl(cid_PCAPHB)=ystates0l(cid_PCAPHB)+solflx%TRCPHB/chemvar%VLWatMicPPB

    !PCPMB=PCPMB+TRCPMB
    ystatesfl(cid_PCAPMB)=ystates0l(cid_PCAPMB)+solflx%TRCPMB/chemvar%VLWatMicPPB

  else
    ystatesfl(cid_CH1PB)=ystates0l(cid_CH1PB)
    ystatesfl(cid_CH2PB)=ystates0l(cid_CH2PB)
    ystatesfl(cid_XH11B)=ystates0l(cid_XH11B)
    ystatesfl(cid_XH21B)=ystates0l(cid_XH21B)
    ystatesfl(cid_X1P1B)=ystates0l(cid_X1P1B)
    ystatesfl(cid_X2P1B)=ystates0l(cid_X2P1B)
    ystatesfl(cid_PALPOB)=ystates0l(cid_PALPOB)
    ystatesfl(cid_PFEPOB)=ystates0l(cid_PFEPOB)
    ystatesfl(cid_PCAPDB)=ystates0l(cid_PCAPDB)
    ystatesfl(cid_PCAPHB)=ystates0l(cid_PCAPHB)
    ystatesfl(cid_PCAPMB)=ystates0l(cid_PCAPMB)
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

  chemvar%XOH11   =ystates0l(cid_XOH11)
  chemvar%XN41    =ystates0l(cid_XN41)
  chemvar%XN4B    =ystates0l(cid_XN4B)
  chemvar%CH1PB   =ystates0l(cid_CH1PB)  !H1POB
  chemvar%CH1P1   =ystates0l(cid_CH1P1)  !H1PO4
  chemvar%CH2P1   =ystates0l(cid_CH2P1)  !H2PO4
  chemvar%CH2PB   =ystates0l(cid_CH2PB)  !H2POB
  chemvar%X1P1B   =ystates0l(cid_X1P1B)
  chemvar%X2P1B   =ystates0l(cid_X2P1B)
  chemvar%XH11B   =ystates0l(cid_XH11B)
  chemvar%XH1P1   =ystates0l(cid_XH1P1)
  chemvar%XH21B   =ystates0l(cid_XH21B)
  chemvar%XH2P1   =ystates0l(cid_XH2P1)
  chemvar%XOH21   =ystates0l(cid_XOH21)
  chemvar%CN31    =ystates0l(cid_CN31)    !ZNH3S
  chemvar%CN3B    =ystates0l(cid_CN3B)    !ZNH3B
  chemvar%CN41    =ystates0l(cid_CN41)    !ZNH4S
  chemvar%CN4B    =ystates0l(cid_CN4B)    !ZNH4B
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
