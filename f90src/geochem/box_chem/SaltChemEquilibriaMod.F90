module SaltChemEquilibriaMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  use SoluteParMod
  use EcosimConst
  use EcoSIMSolverPar
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8), parameter :: ZERO=1.0E-15_r8
  real(r8), parameter :: ZEROS = 1.0E-015_r8
  real(r8), parameter :: ZEROS2= 1.0E-008_r8

  real(r8) :: RAL,RAL1,RAL2,RAL3,RAL4,RALS
  real(r8) :: RB1P,RB2P,RBH0,RBH1,RBH2,RC0B,RC0P
  real(r8) :: RC1B,RC1P,RC2P,RCA,RCAC,RCAH,RCAS,RMGO
  real(r8) :: RC2B,RCAO,RCO2,RCO3,RF1B,RF2B,RFE
  real(r8) :: RFEO1,RFEO2,RFEO3,RFEO4,RFES
  real(r8) :: RFE1,RFE2,RFE3,RFE4,RF1P,RF2P
  real(r8) :: RHB0,RHB3,RHB2,RHCO,RHCO3,RHB1,RHHY,RXH0
  real(r8) :: RHOH,RHP0,RHP1,RHP2,RHP3,RHY,RKA,RKAS
  real(r8) :: RM1B,RH1B,RM1P,RMG
  real(r8) :: RMGC,RMGH,RMGS,RH3B
  real(r8) :: XHY1,XAL1,XFE1
  real(r8) :: XCA1,XMG1,XNA1,XKA1,XHC1,XALO21,XFEO21,XCOOH,XCOO
  real(r8) :: RN3B,RN3S,RN4B,RN4S,ROH
  real(r8) :: RNA,RNAC,RNAS,RPALBX,RPALOX
  real(r8) :: RHAL1,RHALO1,RHALO2,RHALO3,RHALO4,RHFE1
  real(r8) :: R1,PCASO1,PCACO1,PFEOH1,PALOH1
  real(r8) :: CH0PB,CH3PB,CF1PB,CF2PB,CC0PB,CC1PB,CC2PB,CM1PB
  real(r8) :: CH0P1,CH3P1,CF1P1,CF2P1,CC0P1,CC1P1,CC2P1,CM1P1
  real(r8) :: CCAS1,CMGO1,CMGC1,CMGH1,CMGS1,CNAC1,CNAS1,CKAS1
  real(r8) :: CFEO1,CFEO2,CFEO3,CFEO4,CFES1,CCAO1,CCAC1,CCAH1
  real(r8) :: CSO41,CCL1,CHCO31,CALO1,CALO2,CALO3,CALO4,CALS1
  real(r8) :: CNO1,CNOB,CAL1,CCEC,CCO21,CMG1,CNA1,CFE1,CHY1
  real(r8) :: CCO31,CKA1,RPCACX,RH2P,RNH4
  real(r8) :: RHCACO,RPCASO,RHA0P1,RHA1P1,RHA2P1,RHA3P1
  real(r8) :: RHFEO1,RHFEO2,RHFEO3,RHFEO4,RPFEOX,RHCAC3,RHCACH
  real(r8) :: RHA4P1,RHA0P2,RHA1P2,RHA2P2,RHA3P2,RHA4P2,RHF0P1
  real(r8) :: RHF1P1,RHF2P1,RHF3P1,RHF4P1,RHF0P2,RHF1P2,RHF2P2
  real(r8) :: RHF3P2,RHF4P2,RPCAD1,RHCAD2,RHCAH1,RHCAH2,RHA0B1
  real(r8) :: RHA1B1,RHA2B1,RHA3B1,RHA4B1,RHA0B2,RHA1B2,RHA2B2
  real(r8) :: RHA3B2,RHA4B2,RHF0B1,RHF1B1,RHF2B1,RHF3B1,RHF4B1
  real(r8) :: RHF0B2,RHF1B2,RHF2B2,RHF3B2,RHF4B2,RPCDB1,RHCDB2
  real(r8) :: RHCHB1,RHCHB2,RXOH2,RXOH1,RXO2B,RXO1B,RXHY,RXAL
  real(r8) :: ASO41,AOH1,RH1P,RH3P,RPALPX,RPCADX,RPCAHX,RPCAMX
  real(r8) :: ACO21,AH0P1,AH1P1,AH2P1,AH3P1,AF1P1,AF2P1
  real(r8) :: AC0P1,AC1P1,AC2P1,AM1P1,AH0PB,AH1PB,AH2PB,AH3PB
  real(r8) :: A2,AHY1,AAL1,AALO1,AALO2,AALO3,AALO4
  real(r8) :: AFE1,AFEO1,AFEO2,AFEO3,AFEO4,ACA1,ACO31,AHCO31
  real(r8) :: RPCDBX,RPCHBX,RPCMBX,RPFEBX,RSO4,RX1P,RX2P
  real(r8) :: RXFE,RXCA,RXMG,RXNA,RXKA,RXHC,RXALO2,RXFEO2,RCO2Q
  real(r8) :: RXH1,RXH2,RXNB,CCA1,COH1,RXH1P,RXH2B,RPFEPX
  real(r8) :: RALO1,RALO2,RALO3,RALO4,RNHB,RXH1B,RXN4
  real(r8) :: AF1PB,AF2PB,AC0PB,AC1PB,AC2PB,AM1PB,AN41,AN4B
  real(r8) :: AN31,AN3B,AMG1,ANA1,AKA1,AALS1,AFES1,ACAO1,ACAC1
  real(r8) :: RH2B,RXH2P,RYH2B,RYH2P,AKAS1,AALX,AFEX,ACAX,AMGX
  real(r8) :: ACAS1,ACAH1,AMGO1,AMGC1,AMGH1,AMGS1,ANAC1,ANAS1
  real(r8) :: VOLWBK

  real(r8), pointer :: VOLWM
  real(r8), pointer :: VOLWNH
  real(r8), pointer :: VOLWNB
  real(r8), pointer :: VOLWPB
  real(r8), pointer :: VOLWPO
  real(r8), pointer :: VLPOB
  real(r8), pointer :: VLPO4
  real(r8), pointer :: VLNHB
  real(r8), pointer :: VLNH4
  real(r8), pointer :: VLNOB
  real(r8), pointer :: VLNO3
  real(r8), pointer :: BKVL
  real(r8), pointer :: BKVLX
  real(r8), pointer :: VOLWNO
  real(r8), pointer :: VOLWNZ
  real(r8), pointer :: GKC4     !Ca-NH4 Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCA     !Ca-Al Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCH     !Ca-H Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCM     !Ca-Mg Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCK     !Ca-K Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCN     !Ca-Na Gapon selectivity coefficient, [-]
  real(r8), pointer :: XCEC     !cation exchange capacity, [mol d-2]
  real(r8), pointer :: XAEC     !anion exchange capacity, [mol d-2]
  real(r8), pointer :: ORGC     !total soil organic C [g d-2]
  real(r8), pointer :: PH       !pH value of the system
  real(r8), pointer :: XZHYS    !total H+ production, [flux]

  real(r8), pointer :: CO2S     !aqueous CO2  micropore	[g d-2]
  real(r8), pointer :: CH1P1    !soil aqueous HPO4 content micropore non-band, [mol m-3]
  real(r8), pointer :: CH1PB    !soil aqueous HPO4 content micropore band, [mol m-3]
  real(r8), pointer :: CH2P1    !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  real(r8), pointer :: CH2PB    !soil aqueous H2PO4 content micropore  band, [mol m-3]
  real(r8), pointer :: CN31     !soil NH3 concentration in non-band soil, [mol m-3]
  real(r8), pointer :: CN3B     !soil NH3 concentration in band soil, [mol m-3]
  real(r8), pointer :: CN41     !soil NH4 concentration in non-band soil, [mol m-3]
  real(r8), pointer :: CN4B     !soil NH4 concentration in band soil, [mol m-3]
  real(r8), pointer :: ZNO3S    !NO3 mass non-band micropore, [g d-2]
  real(r8), pointer :: ZNO3B    !NO3 mass band micropore, [g d-2]
  real(r8), pointer :: ZHY      !soil aqueous H content micropore, [mol d-2]
  real(r8), pointer :: ZOH      !soil aqueous OH content micropore, [mol d-2]
  real(r8), pointer :: ZAL      !soil aqueous Al content micropore, [mol d-2]
  real(r8), pointer :: ZFE      !soil aqueous Fe content micropore, [mol d-2]
  real(r8), pointer :: ZCA      !soil aqueous Ca content micropore, [mol d-2]
  real(r8), pointer :: ZMG      !soil aqueous Mg content micropore, [mol d-2]
  real(r8), pointer :: ZNA      !soil aqueous Na content micropore, [mol d-2]
  real(r8), pointer :: ZKA      !soil aqueous K content micropore, [mol d-2]
  real(r8), pointer :: ZSO4     !soil aqueous SO4 content micropore, [mol d-2]
  real(r8), pointer :: ZCL      !soil aqueous Cl content micropore, [mol d-2]
  real(r8), pointer :: ZCO3     !soil aqueous CO3 content micropore, [mol d-2]
  real(r8), pointer :: ZHCO3    !soil aqueous HCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZALOH1   !soil aqueous AlOH content micropore, [mol d-2]
  real(r8), pointer :: ZALOH2   !soil aqueous AlOH2 content micropore, [mol d-2]
  real(r8), pointer :: ZALOH3   !soil aqueous AlOH3 content micropore, [mol d-2]
  real(r8), pointer :: ZALOH4   !soil aqueous AlOH4 content micropore, [mol d-2]
  real(r8), pointer :: ZALS     !soil aqueous AlSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH1   !soil aqueous FeOH content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH2   !soil aqueous FeOH2 content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH3   !soil aqueous FeOH3 content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH4   !soil aqueous FeOH4 content micropore, [mol d-2]
  real(r8), pointer :: ZFES     !soil aqueous FeSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZCAO     !soil aqueous CaOH2 content micropore, [mol d-2]
  real(r8), pointer :: ZCAC     !soil aqueous CACO3 content micropore, [mol d-2]
  real(r8), pointer :: ZCAH     !soil aqueous CaHCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZCAS     !soil aqueous CaSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZMGO     !soil aqueous MgOH content micropore, [mol d-2]
  real(r8), pointer :: ZMGC     !soil aqueous MgCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZMGH     !soil aqueous MgHCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZMGS     !soil aqueous MgSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZNAC     !soil aqueous NaCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZNAS     !soil aqueous NaSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZKAS     !soil aqueous KSO4 content micropore, [mol d-2]
  real(r8), pointer :: H0PO4    !soil aqueous PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: H3PO4    !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZFE1P    !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZFE2P    !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZCA0P    !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZCA1P    !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZCA2P    !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZMG1P    !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  real(r8), pointer :: H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  real(r8), pointer :: ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZCA2PB   !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  real(r8), pointer :: PALPO1   !precipitated AlPO4 non-band, [mol m-3]
  real(r8), pointer :: PALPOB   !precipitated AlPO4 band soil, [mol m-3]
  real(r8), pointer :: PCAPD1   !precipitated CaHPO4 non-band soil, [mol m-3]
  real(r8), pointer :: PCAPDB   !precipitated CaHPO4 band soil, [mol m-3]
  real(r8), pointer :: PCAPH1   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  real(r8), pointer :: PCAPHB   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  real(r8), pointer :: PCAPM1   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  real(r8), pointer :: PCAPMB   !precipitated CaH2PO4 band soil, [mol m-3]
  real(r8), pointer :: PFEPO1   !precipitated FePO4 non-band soil, [mol m-3]
  real(r8), pointer :: PFEPOB   !precipitated FePO4 band soil, [mol m-3]
  real(r8), pointer :: PALOH    !precipitated Al(OH)3, [mol d-2]
  real(r8), pointer :: PFEOH    !precipitated Fe(OH)3, [mol d-2]
  real(r8), pointer :: PCACO    !precipitated CaCO3, [mol d-2]
  real(r8), pointer :: PCASO    !precipitated CaSO4, [mol d-2]
  real(r8), pointer :: XHY      !exchangeable H(+) , [mol d-2]
  real(r8), pointer :: XAL      !exchangeable Al, [mol d-2]
  real(r8), pointer :: XFE      !exchangeable Fe, [mol d-2]
  real(r8), pointer :: XCA      !exchangeable Ca, [mol d-2]
  real(r8), pointer :: XMG      !exchangeable Mg, [mol d-2]
  real(r8), pointer :: XNA      !exchangeable Na, [mol d-2]
  real(r8), pointer :: XKA      !exchangeable K, [mol d-2]
  real(r8), pointer :: XHC      !exchangeable COOH , [mol d-2]
  real(r8), pointer :: XOH11    !exchangeable OH  non-band, [mol d-2]
  real(r8), pointer :: XALO2    !exchangeable AlOH2 , [mol d-2]
  real(r8), pointer :: XFEO2    !exchangeable Fe(OH)2, [mol d-2]
  real(r8), pointer :: XN41     !exchangeable NH4 non-band soil, [mol d-2]
  real(r8), pointer :: XN4B     !exchangeable NH4 band soil, [mol d-2]
  real(r8), pointer :: XH01B    !exchangeable OH- band, [mol d-2]
  real(r8), pointer :: XOH01    !exchangeable OH- non-band, [mol d-2]
  real(r8), pointer :: X1P1B    !exchangeable HPO4 concentration band-soil, [mol m-3]
  real(r8), pointer :: X2P1B    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  real(r8), pointer :: XH11B    !exchangeable OH band-soil, [mol m-3]
  real(r8), pointer :: XH1P1    !exchangeable HPO4  non-band, [mol m-3]
  real(r8), pointer :: XH21B    !exchangeable OH2 band-soil, [mol m-3]
  real(r8), pointer :: XOH21    !exchangeable OH2  non-band soil, [mol m-3]
  real(r8), pointer :: XH2P1    !exchangeable H2PO4  non-band soil, [mol m-3]

! fluxes
  real(r8), pointer :: TRCACO   !total precipitated CaCO3 transformation, [mol d-2 h-1]
  real(r8), pointer :: TRNAC
  real(r8), pointer :: TRMGC
  real(r8), pointer :: TRCAC
  real(r8), pointer :: TRMGH
  real(r8), pointer :: TRCAH
  real(r8), pointer :: TRHCO
  real(r8), pointer :: TRXHC
  real(r8), pointer :: TRC2P
  real(r8), pointer :: TRF2P
  real(r8), pointer :: TRC2B
  real(r8), pointer :: TRF2B
  real(r8), pointer :: TRM1P
  real(r8), pointer :: TRC1P
  real(r8), pointer :: TRF1P
  real(r8), pointer :: TRM1B
  real(r8), pointer :: TRC1B
  real(r8), pointer :: TRF1B
  real(r8), pointer :: TRC0B
  real(r8), pointer :: TRH0B
  real(r8), pointer :: TRC0P
  real(r8), pointer :: TRH0P
  real(r8), pointer :: TRFEPO
  real(r8), pointer :: TRALPO
  real(r8), pointer :: TRFEPB
  real(r8), pointer :: TRALPB
  real(r8), pointer :: TRCAPM
  real(r8), pointer :: TRCAPD
  real(r8), pointer :: TRCPMB
  real(r8), pointer :: TRCPDB
  real(r8), pointer :: TRCPHB
  real(r8), pointer :: TRCAPH
  real(r8), pointer :: TRFE1
  real(r8), pointer :: TRAL1
  real(r8), pointer :: TRFE2
  real(r8), pointer :: TRAL2
  real(r8), pointer :: TRFE3
  real(r8), pointer :: TRAL3
  real(r8), pointer :: TRFE4
  real(r8), pointer :: TRAL4
  real(r8), pointer :: TRMGO
  real(r8), pointer :: TRCAO
  real(r8), pointer :: TRFEOH
  real(r8), pointer :: TRALOH
  real(r8), pointer :: TRXFE2
  real(r8), pointer :: TRAL
  real(r8), pointer :: TRALS
  real(r8), pointer :: TRB1P
  real(r8), pointer :: TRB2P
  real(r8), pointer :: TRBH0
  real(r8), pointer :: TRBH1
  real(r8), pointer :: TRBH2
  real(r8), pointer :: TRCA
  real(r8), pointer :: TRCAS
  real(r8), pointer :: TRCASO
  real(r8), pointer :: TRCO2
  real(r8), pointer :: TRCO3
  real(r8), pointer :: TRFE
  real(r8), pointer :: TRFES
  real(r8), pointer :: TRH1B
  real(r8), pointer :: TRH2B
  real(r8), pointer :: TRH3B
  real(r8), pointer :: TRH1P
  real(r8), pointer :: TRH2P
  real(r8), pointer :: TRH3P
  real(r8), pointer :: TRHY
  real(r8), pointer :: TRKA
  real(r8), pointer :: TRKAS
  real(r8), pointer :: TRMG
  real(r8), pointer :: TRMGS
  real(r8), pointer :: TRN3B
  real(r8), pointer :: TRN3S
  real(r8), pointer :: TRN4B
  real(r8), pointer :: TRN4S
  real(r8), pointer :: TRNA
  real(r8), pointer :: TRNAS
  real(r8), pointer :: TROH
  real(r8), pointer :: TRSO4
  real(r8), pointer :: TRX1P
  real(r8), pointer :: TRX2P
  real(r8), pointer :: TRXAL
  real(r8), pointer :: TRXAL2   !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  real(r8), pointer :: TRXCA
  real(r8), pointer :: TRXFE
  real(r8), pointer :: TRXH0
  real(r8), pointer :: TRXH1
  real(r8), pointer :: TRXH2
  real(r8), pointer :: TRXHY
  real(r8), pointer :: TRXKA
  real(r8), pointer :: TRXMG
  real(r8), pointer :: TRXN4
  real(r8), pointer :: TRXNA
  real(r8), pointer :: TRXNB
  real(r8), pointer :: TBCO2
  real(r8), pointer :: TBION
  real(r8), pointer :: TRH2O
  public :: SaltChemEquilibria
  contains

!--------------------------------------------------------------------------
  subroutine SaltChemEquilibria(chemvar,solflx)

  implicit none
  type(solute_flx_type),target, intent(inout) :: solflx
  type(chem_var_type),target, intent(inout) :: chemvar

!     begin_execution

  BKVLX  => chemvar%BKVLX
  VOLWNZ => chemvar%VOLWNZ
  VOLWNO => chemvar%VOLWNO
  VOLWNB => chemvar%VOLWNB
  VOLWNH => chemvar%VOLWNH
  VOLWPB => chemvar%VOLWPB
  VOLWPO => chemvar%VOLWPO
  XOH11  => chemvar%XOH11
  XOH21  => chemvar%XOH21
  XN41   => chemvar%XN41
  XN4B   => chemvar%XN4B
  CH1P1  => chemvar%CH1P1
  CH1PB  => chemvar%CH1PB
  CH2P1  => chemvar%CH2P1
  CH2PB  => chemvar%CH2PB
  X1P1B  => chemvar%X1P1B
  X2P1B  => chemvar%X2P1B
  XH11B  => chemvar%XH11B
  XH1P1  => chemvar%XH1P1
  XH21B  => chemvar%XH21B
  XH2P1  => chemvar%XH2P1
  CN31   => chemvar%CN31
  CN3B   => chemvar%CN3B
  CN41   => chemvar%CN41
  CN4B   => chemvar%CN4B
  PALPO1 => chemvar%PALPO1
  PALPOB => chemvar%PALPOB
  PCAPD1 => chemvar%PCAPD1
  PCAPDB => chemvar%PCAPDB
  PCAPH1 => chemvar%PCAPH1
  PCAPHB => chemvar%PCAPHB
  PCAPM1 => chemvar%PCAPM1
  PCAPMB => chemvar%PCAPMB
  PFEPO1 => chemvar%PFEPO1
  PFEPOB => chemvar%PFEPOB
  ZNO3S  => chemvar%ZNO3S
  ZNO3B  => chemvar%ZNO3B
  VOLWM  => chemvar%VOLWM
  XZHYS  => chemvar%XZHYS
  ZHY    => chemvar%ZHY
  XCEC   => chemvar%XCEC
  ZOH    => chemvar%ZOH
  ZAL    => chemvar%ZAL
  ZFE    => chemvar%ZFE
  ZCA    => chemvar%ZCA
  ZMG    => chemvar%ZMG
  ZNA    => chemvar%ZNA
  ZKA    => chemvar%ZKA
  ZSO4   => chemvar%ZSO4
  ZCL    => chemvar%ZCL
  ZCO3   => chemvar%ZCO3
  ZHCO3  => chemvar%ZHCO3
  CO2S   => chemvar%CO2S
  ZALOH1 => chemvar%ZALOH1
  ZALOH2 => chemvar%ZALOH2
  ZALOH3 => chemvar%ZALOH3
  ZALOH4 => chemvar%ZALOH4
  ZALS   => chemvar%ZALS
  ZFEOH1 => chemvar%ZFEOH1
  ZFEOH2 => chemvar%ZFEOH2
  ZFEOH3 => chemvar%ZFEOH3
  ZFEOH4 => chemvar%ZFEOH4
  ZFES   => chemvar%ZFES
  ZCAO   => chemvar%ZCAO
  ZCAC   => chemvar%ZCAC
  ZCAH   => chemvar%ZCAH
  ZCAS   => chemvar%ZCAS
  ZMGO   => chemvar%ZMGO
  ZMGC   => chemvar%ZMGC
  ZMGH   => chemvar%ZMGH
  ZMGS   => chemvar%ZMGS
  ZNAC   => chemvar%ZNAC
  ZNAS   => chemvar%ZNAS
  ZKAS   => chemvar%ZKAS
  H0PO4  => chemvar%H0PO4
  H3PO4  => chemvar%H3PO4
  ZFE1P  => chemvar%ZFE1P
  ZFE2P  => chemvar%ZFE2P
  ZCA0P  => chemvar%ZCA0P
  ZCA1P  => chemvar%ZCA1P
  ZCA2P  => chemvar%ZCA2P
  ZMG1P  => chemvar%ZMG1P
  H0POB  => chemvar%H0POB
  H3POB  => chemvar%H3POB
  ZFE1PB => chemvar%ZFE1PB
  ZFE2PB => chemvar%ZFE2PB
  ZCA0PB => chemvar%ZCA0PB
  ZCA1PB => chemvar%ZCA1PB
  ZCA2PB => chemvar%ZCA2PB
  ZMG1PB => chemvar%ZMG1PB
  XHY    => chemvar%XHY
  XAL    => chemvar%XAL
  XFE    => chemvar%XFE
  XCA    => chemvar%XCA
  XMG    => chemvar%XMG
  XNA    => chemvar%XNA
  XKA    => chemvar%XKA
  XHC    => chemvar%XHC
  XALO2  => chemvar%XALO2
  XFEO2  => chemvar%XFEO2
  ORGC   => chemvar%ORGC
  PALOH  => chemvar%PALOH
  PFEOH  => chemvar%PFEOH
  PCACO  => chemvar%PCACO
  PCASO  => chemvar%PCASO
  PH     => chemvar%PH
  VLPOB  => chemvar%VLPOB
  VLPO4  => chemvar%VLPO4
  VLNHB  => chemvar%VLNHB
  VLNH4  => chemvar%VLNH4
  VLNOB  => chemvar%VLNOB
  VLNO3  => chemvar%VLNO3
  BKVL   => chemvar%BKVL
  XAEC   => chemvar%XAEC
  GKC4   => chemvar%GKC4
  GKCA   => chemvar%GKCA
  GKCH   => chemvar%GKCH
  GKCM   => chemvar%GKCM
  GKCK   => chemvar%GKCK
  GKCN   => chemvar%GKCN
  XH01B  => chemvar%XH01B
  XOH01  => chemvar%XOH01

  TRCACO => solflx%TRCACO
  TRNAC  => solflx%TRNAC
  TRMGC  => solflx%TRMGC
  TRCAC  => solflx%TRCAC
  TRMGH  => solflx%TRMGH
  TRCAH  => solflx%TRCAH
  TRHCO  => solflx%TRHCO
  TRXHC  => solflx%TRXHC
  TRC2P  => solflx%TRC2P
  TRF2P  => solflx%TRF2P
  TRC2B  => solflx%TRC2B
  TRF2B  => solflx%TRF2B
  TRM1P  => solflx%TRM1P
  TRC1P  => solflx%TRC1P
  TRF1P  => solflx%TRF1P
  TRM1B  => solflx%TRM1B
  TRC1B  => solflx%TRC1B
  TRF1B  => solflx%TRF1B
  TRC0B  => solflx%TRC0B
  TRH0B  => solflx%TRH0B
  TRC0P  => solflx%TRC0P
  TRH0P  => solflx%TRH0P
  TRFEPO => solflx%TRFEPO
  TRALPO => solflx%TRALPO
  TRFEPB => solflx%TRFEPB
  TRALPB => solflx%TRALPB
  TRCAPM => solflx%TRCAPM
  TRCAPD => solflx%TRCAPD
  TRCPMB => solflx%TRCPMB
  TRCPDB => solflx%TRCPDB
  TRCPHB => solflx%TRCPHB
  TRCAPH => solflx%TRCAPH
  TRFE1  => solflx%TRFE1
  TRAL1  => solflx%TRAL1
  TRFE2  => solflx%TRFE2
  TRAL2  => solflx%TRAL2
  TRFE3  => solflx%TRFE3
  TRAL3  => solflx%TRAL3
  TRFE4  => solflx%TRFE4
  TRAL4  => solflx%TRAL4
  TRMGO  => solflx%TRMGO
  TRCAO  => solflx%TRCAO
  TRFEOH => solflx%TRFEOH
  TRALOH => solflx%TRALOH
  TRXFE2 => solflx%TRXFE2
  TRAL   => solflx%TRAL
  TRALS  => solflx%TRALS
  TRB1P  => solflx%TRB1P
  TRB2P  => solflx%TRB2P
  TRBH0  => solflx%TRBH0
  TRBH1  => solflx%TRBH1
  TRBH2  => solflx%TRBH2
  TRCA   => solflx%TRCA
  TRCAS  => solflx%TRCAS
  TRCASO => solflx%TRCASO
  TRCO2  => solflx%TRCO2
  TRCO3  => solflx%TRCO3
  TRFE   => solflx%TRFE
  TRFES  => solflx%TRFES
  TRH1B  => solflx%TRH1B
  TRH2B  => solflx%TRH2B
  TRH3B  => solflx%TRH3B
  TRH1P  => solflx%TRH1P
  TRH2P  => solflx%TRH2P
  TRH3P  => solflx%TRH3P
  TRHY   => solflx%TRHY
  TRKA   => solflx%TRKA
  TRKAS  => solflx%TRKAS
  TRMG   => solflx%TRMG
  TRMGS  => solflx%TRMGS
  TRN3B  => solflx%TRN3B
  TRN3S  => solflx%TRN3S
  TRN4B  => solflx%TRN4B
  TRN4S  => solflx%TRN4S
  TRNA   => solflx%TRNA
  TRNAS  => solflx%TRNAS
  TROH   => solflx%TROH
  TRSO4  => solflx%TRSO4
  TRX1P  => solflx%TRX1P
  TRX2P  => solflx%TRX2P
  TRXAL  => solflx%TRXAL
  TRXAL2 => solflx%TRXAL2
  TRXCA  => solflx%TRXCA
  TRXFE  => solflx%TRXFE
  TRXH0  => solflx%TRXH0
  TRXH1  => solflx%TRXH1
  TRXH2  => solflx%TRXH2
  TRXHY  => solflx%TRXHY
  TRXKA  => solflx%TRXKA
  TRXMG  => solflx%TRXMG
  TRXN4  => solflx%TRXN4
  TRXNA  => solflx%TRXNA
  TRXNB  => solflx%TRXNB
  TBCO2  => solflx%TBCO2
  TBION  => solflx%TBION
  TRH2O  => solflx%TRH2O
!
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     C*,X*,Z*=soluble,exchangeable concentration, mass
!
  call PrepIonConcentrations
!
!     CONVERGENCE TOWARDS SOLUTE EQILIBRIA
!
  call SolveChemEquilibria
!

  call SummarizeIonFluxes

  end subroutine SaltChemEquilibria

!------------------------------------------------------------------------
  subroutine PrepIonConcentrations
  implicit none
  real(r8) :: VOLWPX
!     begin_execution
!     SOLUBLE NO3 CONCENTRATIONS
!     IN NON-BAND AND BAND SOIL ZONES
!
!     VOLWNO,VOLWNZ=soil water volume in NO3 non-band,band
!     ZNO3S,ZNO3B=NO3 mass in non-band,band
!     CNO1,CNOB=NO3 concentrations in non-band,band
!
  IF(VOLWNO.GT.ZEROS2)THEN
    CNO1=AMAX1(0._r8,ZNO3S/(natomw*VOLWNO))
  ELSE
    CNO1=0._r8
  ENDIF

  IF(VOLWNZ.GT.ZEROS2)THEN
    CNOB=AMAX1(0._r8,ZNO3B/(natomw*VOLWNZ))
  ELSE
    CNOB=0._r8
  ENDIF
!
!     H CONCENTRATION
!
!     XZHYS=total H+ production from nitro.f
!
  CHY1=AMAX1(0._r8,ZHY+XZHYS)/VOLWM
!
!     SOLUTE ION AND ION PAIR CONCENTRATIONS
!
!     CCEC,XCEC=cation exchange concentration,capacity
!     BKVLX=soil mass
!     VOLWM=soil water volume
!
  IF(BKVLX.GT.ZEROS)THEN
    CCEC=AMAX1(ZERO,XCEC/BKVLX)
  ELSE
    CCEC=ZERO
  ENDIF
  IF(VOLWM.GT.ZEROS2)THEN
    COH1=AMAX1(0._r8,ZOH/VOLWM)
    CAL1=AMAX1(0._r8,ZAL/VOLWM)
    CFE1=AMAX1(0._r8,ZFE/VOLWM)
    CCA1=AMAX1(0._r8,ZCA/VOLWM)
    CMG1=AMAX1(0._r8,ZMG/VOLWM)
    CNA1=AMAX1(0._r8,ZNA/VOLWM)
    CKA1=AMAX1(0._r8,ZKA/VOLWM)
    CSO41=AMAX1(0._r8,ZSO4/VOLWM)
    CCL1=AMAX1(0._r8,ZCL/VOLWM)
    CCO31=AMAX1(0._r8,ZCO3/VOLWM)
    CHCO31=AMAX1(0._r8,ZHCO3/VOLWM)
    CCO21=AMAX1(0._r8,CO2S/(catomw*VOLWM))
    CALO1=AMAX1(0._r8,ZALOH1/VOLWM)
    CALO2=AMAX1(0._r8,ZALOH2/VOLWM)
    CALO3=AMAX1(0._r8,ZALOH3/VOLWM)
    CALO4=AMAX1(0._r8,ZALOH4/VOLWM)
    CALS1=AMAX1(0._r8,ZALS/VOLWM)
    CFEO1=AMAX1(0._r8,ZFEOH1/VOLWM)
    CFEO2=AMAX1(0._r8,ZFEOH2/VOLWM)
    CFEO3=AMAX1(0._r8,ZFEOH3/VOLWM)
    CFEO4=AMAX1(0._r8,ZFEOH4/VOLWM)
    CFES1=AMAX1(0._r8,ZFES/VOLWM)
    CCAO1=AMAX1(0._r8,ZCAO/VOLWM)
    CCAC1=AMAX1(0._r8,ZCAC/VOLWM)
    CCAH1=AMAX1(0._r8,ZCAH/VOLWM)
    CCAS1=AMAX1(0._r8,ZCAS/VOLWM)
    CMGO1=AMAX1(0._r8,ZMGO/VOLWM)
    CMGC1=AMAX1(0._r8,ZMGC/VOLWM)
    CMGH1=AMAX1(0._r8,ZMGH/VOLWM)
    CMGS1=AMAX1(0._r8,ZMGS/VOLWM)
    CNAC1=AMAX1(0._r8,ZNAC/VOLWM)
    CNAS1=AMAX1(0._r8,ZNAS/VOLWM)
    CKAS1=AMAX1(0._r8,ZKAS/VOLWM)
  ELSE
    COH1=0._r8
    CAL1=0._r8
    CFE1=0._r8
    CCA1=0._r8
    CMG1=0._r8
    CNA1=0._r8
    CKA1=0._r8
    CSO41=0._r8
    CCL1=0._r8
    CCO31=0._r8
    CHCO31=0._r8
    CCO21=0._r8
    CALO1=0._r8
    CALO2=0._r8
    CALO3=0._r8
    CALO4=0._r8
    CALS1=0._r8
    CFEO1=0._r8
    CFEO2=0._r8
    CFEO3=0._r8
    CFEO4=0._r8
    CFES1=0._r8
    CCAO1=0._r8
    CCAC1=0._r8
    CCAH1=0._r8
    CCAS1=0._r8
    CMGO1=0._r8
    CMGC1=0._r8
    CMGH1=0._r8
    CMGS1=0._r8
    CNAC1=0._r8
    CNAS1=0._r8
    CKAS1=0._r8
  ENDIF
!
!     PO4 CONCENTRATIONS IN NON-BAND AND BAND SOIL ZONES
!
!     VOLWPO,VOLWPB=water volume in PO4 non-band,band
!
  IF(VOLWPO.GT.ZEROS2)THEN
    VOLWPX=patomw*VOLWPO
    CH0P1=AMAX1(0._r8,H0PO4/VOLWPO)
    CH3P1=AMAX1(0._r8,H3PO4/VOLWPO)
    CF1P1=AMAX1(0._r8,ZFE1P/VOLWPO)
    CF2P1=AMAX1(0._r8,ZFE2P/VOLWPO)
    CC0P1=AMAX1(0._r8,ZCA0P/VOLWPO)
    CC1P1=AMAX1(0._r8,ZCA1P/VOLWPO)
    CC2P1=AMAX1(0._r8,ZCA2P/VOLWPO)
    CM1P1=AMAX1(0._r8,ZMG1P/VOLWPO)
  ELSE
    CH0P1=0._r8
    CH3P1=0._r8
    CF1P1=0._r8
    CF2P1=0._r8
    CC0P1=0._r8
    CC1P1=0._r8
    CC2P1=0._r8
    CM1P1=0._r8
  ENDIF
  IF(VOLWPB.GT.ZEROS2)THEN
    CH0PB=AMAX1(0._r8,H0POB/VOLWPB)
    CH3PB=AMAX1(0._r8,H3POB/VOLWPB)
    CF1PB=AMAX1(0._r8,ZFE1PB/VOLWPB)
    CF2PB=AMAX1(0._r8,ZFE2PB/VOLWPB)
    CC0PB=AMAX1(0._r8,ZCA0PB/VOLWPB)
    CC1PB=AMAX1(0._r8,ZCA1PB/VOLWPB)
    CC2PB=AMAX1(0._r8,ZCA2PB/VOLWPB)
    CM1PB=AMAX1(0._r8,ZMG1PB/VOLWPB)
  ELSE
    CH0PB=0._r8
    CH3PB=0._r8
    CF1PB=0._r8
    CF2PB=0._r8
    CC0PB=0._r8
    CC1PB=0._r8
    CC2PB=0._r8
    CM1PB=0._r8
  ENDIF
!
!     EXCHANGEABLE ION CONCENTRATIONS
!
  IF(BKVLX.GT.ZEROS)THEN
    XHY1=AMAX1(0._r8,XHY/BKVLX)
    XAL1=AMAX1(0._r8,XAL/BKVLX)
    XFE1=AMAX1(0._r8,XFE/BKVLX)
    XCA1=AMAX1(0._r8,XCA/BKVLX)
    XMG1=AMAX1(0._r8,XMG/BKVLX)
    XNA1=AMAX1(0._r8,XNA/BKVLX)
    XKA1=AMAX1(0._r8,XKA/BKVLX)
    XHC1=AMAX1(0._r8,XHC/BKVLX)
    XALO21=AMAX1(0._r8,XALO2/BKVLX)
    XFEO21=AMAX1(0._r8,XFEO2/BKVLX)
    XCOOH=AMAX1(0._r8,COOH*ORGC/BKVLX)
!
!     PRECIPITATE CONCENTRATIONS
!
    PALOH1=AMAX1(0._r8,PALOH/BKVLX)
    PFEOH1=AMAX1(0._r8,PFEOH/BKVLX)
    PCACO1=AMAX1(0._r8,PCACO/BKVLX)
    PCASO1=AMAX1(0._r8,PCASO/BKVLX)
  ELSE
    XHY1=0._r8
    XAL1=0._r8
    XFE1=0._r8
    XCA1=0._r8
    XMG1=0._r8
    XNA1=0._r8
    XKA1=0._r8
    XHC1=0._r8
    XALO21=0._r8
    XFEO21=0._r8
    XCOOH=0._r8
    PALOH1=0._r8
    PFEOH1=0._r8
    PCACO1=0._r8
    PCASO1=0._r8
  ENDIF
  end subroutine PrepIonConcentrations
!--------------------------------------------------------------------------
  subroutine SolveChemEquilibria

  implicit none
  real(r8) :: SP,SPX,S0,S1,P1,P2,PX
  integer :: M,NR1,NP2
!     begin_execution

  DO 1000 M=1,MRXN
!
!     SOLUTE CONCENTRATIONS
!
    call GetSoluteConcentrations
!
    call IonStrengthActivity
!
!     ALUMINUM HYDROXIDE (GIBBSITE)
!
    PX=AMAX1(AAL1,AALO1,AALO2,AALO3,AALO4)
    IF(isclose(PX,AAL1))THEN
      R1=AHY1
      P1=AAL1
      P2=AOH1
      NR1=3
      NP2=0
      SP=SHALO
    ELSEIF(isclose(PX,AALO1))THEN
      R1=AHY1
      P1=AALO1
      P2=AOH1
      NR1=2
      NP2=0
      SP=SHAL1
    ELSEIF(isclose(PX,AALO2))THEN
      R1=AHY1
      P1=AALO2
      P2=AOH1
      NR1=1
      NP2=0
      SP=SHAL2
    ELSEIF(isclose(PX,AALO3))THEN
      R1=AHY1
      P1=AALO3
      P2=AOH1
      NR1=0
      NP2=0
      SP=SPAL3
    ELSEIF(isclose(PX,AALO4))THEN
      R1=AOH1
      P1=AALO4
      P2=AHY1
      NR1=0
      NP2=1
      SP=SHAL4
    ENDIF
    RHAL1=0._r8
    RHALO1=0._r8
    RHALO2=0._r8
    RHALO3=0._r8
    RHALO4=0._r8
!     R1=AMAX1(ZERO,R1)
!     P1=AMAX1(ZERO,P1)
!     P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPALOX=AMAX1(-PALOH1,TPDX*(P1-SPX))
    IF(isclose(PX,AAL1))THEN
      RHAL1=RPALOX
    ELSEIF(isclose(PX,AALO1))THEN
      RHALO1=RPALOX
    ELSEIF(isclose(PX,AALO2))THEN
      RHALO2=RPALOX
    ELSEIF(isclose(PX,AALO3))THEN
      RHALO3=RPALOX
    ELSEIF(isclose(PX,AALO4))THEN
      RHALO4=RPALOX
    ENDIF
!
!     IRON HYDROXIDE
!
    PX=AMAX1(AFE1,AFEO1,AFEO2,AFEO3,AFEO4)
    IF(isclose(PX,AFE1))THEN
      R1=AHY1
      P1=AFE1
      P2=AOH1
      NR1=3
      NP2=0
      SP=SHFEO
    ELSEIF(isclose(PX,AFEO1))THEN
      R1=AHY1
      P1=AFEO1
      P2=AOH1
      NR1=2
      NP2=0
      SP=SHFE1
    ELSEIF(isclose(PX,AFEO2))THEN
      R1=AHY1
      P1=AFEO2
      P2=AOH1
      NR1=1
      NP2=0
      SP=SHFE2
    ELSEIF(isclose(PX,AFEO3))THEN
      R1=AHY1
      P1=AFEO3
      P2=AOH1
      NR1=0
      NP2=0
      SP=SPFE3
    ELSEIF(isclose(PX,AFEO4))THEN
      R1=AOH1
      P1=AFEO4
      P2=AHY1
      NR1=0
      NP2=1
      SP=SHFE4
    ENDIF
    RHFE1=0._r8
    RHFEO1=0._r8
    RHFEO2=0._r8
    RHFEO3=0._r8
    RHFEO4=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPFEOX=AMAX1(-PFEOH1,TPDX*(P1-SPX))
    IF(isclose(PX,AFE1))THEN
      RHFE1=RPFEOX
    ELSEIF(isclose(PX,AFEO1))THEN
      RHFEO1=RPFEOX
    ELSEIF(isclose(PX,AFEO2))THEN
      RHFEO2=RPFEOX
    ELSEIF(isclose(PX,AFEO3))THEN
      RHFEO3=RPFEOX
    ELSEIF(isclose(PX,AFEO4))THEN
      RHFEO4=RPFEOX
    ENDIF
!     IF(I.EQ.180.AND.J.EQ.12)THEN
!     WRITE(*,1112)'FEOH',I,J,L,M,PFEOH1,AFE1,AFEO1,AFEO2,AFEO3,AFEO4
!    2,AOH1,R1,P1,P2,SP,SPX,RPFEOX,RHFE1,RHFEO1,RHFEO2,RHFEO3,RHFEO4
!    3,AFE1*AOH1**3,SPFEO
!     ENDIF
!
!     CALCITE
!
    PX=AMAX1(ACO31,AHCO31,ACO21)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,ACO31))THEN
      P2=ACO31
      NR1=0
      SP=SPCAC
    ELSEIF(isclose(PX,AHCO31))THEN
      P2=AHCO31
      NR1=1
      SP=SHCAC1
    ELSEIF(isclose(PX,ACO21))THEN
      P2=ACO21
      NR1=2
      SP=SHCAC2
    ENDIF
    RHCAC3=0._r8
    RHCACH=0._r8
    RHCACO=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2._r8-4.0_r8*(P1*P2-SPX))
    RPCACX=AMAX1(-PCACO1,TPDX*(S0-SQRT(S1)))
    IF(isclose(PX,ACO31))THEN
      RHCAC3=RPCACX
    ELSEIF(isclose(PX,AHCO31))THEN
      RHCACH=RPCACX
    ELSEIF(isclose(PX,ACO21))THEN
      RHCACO=RPCACX
    ENDIF
!
!     GYPSUM
!
    P1=ACA1
    P2=ASO41
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SPCAS
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2-4.0*(P1*P2-SPX))
    RPCASO=AMAX1(-PCASO1,TPDX*(S0-SQRT(S1)))
!     IF((M/10)*10.EQ.M)THEN
!     WRITE(*,1112)'CALC',I,J,L,M,PCASO1,ACO31,AHCO31,ACO21,CHY1
!    2,COH1,R1,P1,P2,P3,SP,Z,TX,RPCACX,RHCAC3,RHCACH,RHCACO
!    3,CCA1*A2*CCO3*A2,SPCAC
!     ENDIF
!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN NON-BAND SOIL ZONE
!
    call PhospPrecipDissolNonBand
!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN BAND SOIL ZONE
!
    call PhospPrecipDissolBand
!
    call PhospAnionExchNoBand
!
    call PhospAnionExchBand
!
!     CATION EXCHANGE FROM GAPON SELECTIVITY COEFFICIENTS
!     FOR CA-NH4, CA-H, CA-AL, CA-MG, CA-NA, CA-K
!
    call CationExchange
!
!     SOLUTE DISSOCIATION REACTIONS
!
    call SoluteDissociation
!
!     TOTAL ION FLUXES FOR CURRENT ITERATION
!     FROM ALL REACTIONS ABOVE
!
    call UpdateIonFluxCurentIter

    call UpdateIonConcCurrentIter
!
    call AccumulateIonFlux
!
!     GO TO NEXT ITERATION
!
1000  CONTINUE
!     CONVERGENCE ITERATIONS COMPLETED
!
!     IF(J.EQ.24)THEN
!     WRITE(*,1119)'GAPON',I,J,L,M,CH0P1,CAL1,CFE1,CH0P1*A3*CAL1*A3
!    2,SPALP,CH0P1*A3*CFE1*A3,SPFEP
!    6,SPOH2,XOH11*CHY1*A1/XOH21,SPOH1,XOH01*CHY1*A1/XOH11
!    7,SPH2P,XOH21*CH2P1*A1/XH2P1,SXH2P,XOH11*CH2P1/(XH2P1*COH1)
!    8,SPH1P,XOH11*CH1P1*A2/(XH1P1*COH1*A1)
!    9,COH1*A1,CHY1*A1
!1119  FORMAT(A8,4I4,24E11.3)
!     WRITE(*,1119)'CATION',I,J,L,M,CCEC,XN41+XHY1+3*XAL1+2*(XCA1+XMG1)
!    2+XNA1+XKA1,XN41,XHY1,XAL1,XCA1,XMG1,XNA1,XKA1,CN41,CHY1,CAL1,CCA1
!    2,CMG1,CNA1,CKA1,(CCA1*A2)**0.5*XN41/(CN41*A1*XCA1*2)
!    3,(CCA1*A2)**0.5*XHY1/(CHY1*A1*XCA1*2)
!    2,(CCA1*A2)**0.5*XAL1*3/((CAL1*A3)**0.333*XCA1*2)
!    3,(CCA1*A2)**0.5*XMG1*2/((CMG1*A2)**0.5*XCA1*2)
!    3,(CCA1*A2)**0.5*XNA1/(CNA1*A1*XCA1*2)
!    5,(CCA1*A2)**0.5*XKA1/(CKA1*A1*XCA1*2)
!    6,CHY1*A1*XCOO/XHC1,CALO2*A1*XCOO/XALO21
!     ENDIF
  end subroutine SolveChemEquilibria
!------------------------------------------------------------------------------------------

  subroutine PhospPrecipDissolNonBand

  implicit none
  real(r8) :: SP,SPX,S0,S1,P1,P2,P3,PX,PY
  integer :: NR1,NP3
  IF(VOLWPO.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
    PX=AMAX1(AAL1,AALO1,AALO2,AALO3,AALO4)
    PY=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P3=AHY1
    IF(isclose(PY,AH1P1))THEN
      P2=AH1P1
      IF(isclose(PX,AAL1))THEN
        P1=AAL1
        NR1=1
        NP3=0
        SP=SHA0P1
      ELSEIF(isclose(PX,AALO1))THEN
        P1=AALO1
        NR1=0
        NP3=0
        SP=SPA1P1
      ELSEIF(isclose(PX,AALO2))THEN
        P1=AALO2
        NR1=0
        NP3=1
        SP=SHA2P1
      ELSEIF(isclose(PX,AALO3))THEN
        P1=AALO3
        NR1=0
        NP3=2
        SP=SHA3P1
      ELSEIF(isclose(PX,AALO4))THEN
        P1=AALO4
        NR1=0
        NP3=3
      SP=SHA4P1
      ENDIF
    ELSE
      P2=AH2P1
      IF(isclose(PX,AAL1))THEN
        P1=AAL1
        NR1=2
        NP3=0
        SP=SHA0P2
      ELSEIF(isclose(PX,AALO1))THEN
        P1=AALO1
        NR1=1
        NP3=0
        SP=SHA1P2
      ELSEIF(isclose(PX,AALO2))THEN
        P1=AALO2
        NR1=0
        NP3=0
        SP=SPA2P2
      ELSEIF(isclose(PX,AALO3))THEN
        P1=AALO3
        NR1=0
        NP3=1
        SP=SHA3P2
      ELSEIF(isclose(PX,AALO4))THEN
        P1=AALO4
        NR1=0
        NP3=2
        SP=SHA4P2
      ENDIF
    ENDIF
    RHA0P1=0._r8
    RHA1P1=0._r8
    RHA2P1=0._r8
    RHA3P1=0._r8
    RHA4P1=0._r8
    RHA0P2=0._r8
    RHA1P2=0._r8
    RHA2P2=0._r8
    RHA3P2=0._r8
    RHA4P2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2-4.0*(P1*P2-SPX))
    RPALPX=AMAX1(-PALPO1,TPDX*(S0-SQRT(S1)))
    IF(isclose(PY,AH1P1))THEN
      IF(isclose(PX,AAL1))THEN
        RHA0P1=RPALPX
      ELSEIF(isclose(PX,AALO1))THEN
        RHA1P1=RPALPX
      ELSEIF(isclose(PX,AALO2))THEN
        RHA2P1=RPALPX
      ELSEIF(isclose(PX,AALO3))THEN
        RHA3P1=RPALPX
      ELSEIF(isclose(PX,AALO4))THEN
        RHA4P1=RPALPX
      ENDIF
    ELSE
      IF(isclose(PX,AAL1))THEN
        RHA0P2=RPALPX
      ELSEIF(isclose(PX,AALO1))THEN
        RHA1P2=RPALPX
      ELSEIF(isclose(PX,AALO2))THEN
        RHA2P2=RPALPX
      ELSEIF(isclose(PX,AALO3))THEN
        RHA3P2=RPALPX
      ELSEIF(isclose(PX,AALO4))THEN
        RHA4P2=RPALPX
      ENDIF
    ENDIF
!     IF(I.EQ.180.AND.J.EQ.12)THEN
!     WRITE(*,1112)'ALPO4',I,J,L,M,PALPO1,AAL1,AALO1,AALO2,AALO3,AALO4
!    2,AH0P1,AH1P1,AH2P1,AHY1,AOH1,RPALPX,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    3,RHA4P1,RHA0P2,RHA1P2,RHA2P2,RHA3P2,RHA4P2,SP,SPX,AAL1*AH0P1
!    4,SPALP,CH0P1,CH1P1,CH2P1
!     ENDIF
!1112  FORMAT(A8,4I5,80E12.4)
!     ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
    PX=AMAX1(AFE1,AFEO1,AFEO2,AFEO3,AFEO4)
    PY=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P3=AHY1
    IF(isclose(PY,AH1P1))THEN
      P2=AH1P1
      IF(isclose(PX,AFE1))THEN
        P1=AFE1
        NR1=1
        NP3=0
        SP=SHF0P1
      ELSEIF(isclose(PX,AFEO1))THEN
        P1=AFEO1
        NR1=0
        NP3=0
        SP=SPF1P1
      ELSEIF(isclose(PX,AFEO2))THEN
        P1=AFEO2
        NR1=0
        NP3=1
        SP=SHF2P1
      ELSEIF(isclose(PX,AFEO3))THEN
        P1=AFEO3
        NR1=0
        NP3=2
        SP=SHF3P1
      ELSEIF(isclose(PX,AFEO4))THEN
        P1=AFEO4
        NR1=0
        NP3=3
        SP=SHF4P1
      ENDIF
    ELSE
      P2=AH2P1
      IF(isclose(PX,AFE1))THEN
        P1=AFE1
        NR1=2
        NP3=0
        SP=SHF0P2
      ELSEIF(isclose(PX,AFEO1))THEN
        P1=AFEO1
        NR1=1
        NP3=0
        SP=SHF1P2
      ELSEIF(isclose(PX,AFEO2))THEN
        P1=AFEO2
        NR1=0
        NP3=0
        SP=SPF2P2
      ELSEIF(isclose(PX,AFEO3))THEN
        P1=AFEO3
        NR1=0
        NP3=1
        SP=SHF3P2
      ELSEIF(isclose(PX,AFEO4))THEN
        P1=AFEO4
        NR1=0
        NP3=2
        SP=SHF4P2
      ENDIF
    ENDIF
    RHF0P1=0._r8
    RHF1P1=0._r8
    RHF2P1=0._r8
    RHF3P1=0._r8
    RHF4P1=0._r8
    RHF0P2=0._r8
    RHF1P2=0._r8
    RHF2P2=0._r8
    RHF3P2=0._r8
    RHF4P2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2-4.0*(P1*P2-SPX))
    RPFEPX=AMAX1(-PFEPO1,TPDX*(S0-SQRT(S1)))
    IF(isclose(PY,AH1P1))THEN
      IF(isclose(PX,AFE1))THEN
        RHF0P1=RPFEPX
      ELSEIF(isclose(PX,AFEO1))THEN
        RHF1P1=RPFEPX
      ELSEIF(isclose(PX,AFEO2))THEN
        RHF2P1=RPFEPX
      ELSEIF(isclose(PX,AFEO3))THEN
        RHF3P1=RPFEPX
      ELSEIF(isclose(PX,AFEO4))THEN
        RHF4P1=RPFEPX
      ENDIF
    ELSE
      IF(isclose(PX,AFE1))THEN
        RHF0P2=RPFEPX
      ELSEIF(isclose(PX,AFEO1))THEN
        RHF1P2=RPFEPX
      ELSEIF(isclose(PX,AFEO2))THEN
        RHF2P2=RPFEPX
      ELSEIF(isclose(PX,AFEO3))THEN
        RHF3P2=RPFEPX
      ELSEIF(isclose(PX,AFEO4))THEN
        RHF4P2=RPFEPX
      ENDIF
    ENDIF
!     IF(I.EQ.180.AND.J.EQ.12)THEN
!     WRITE(*,1112)'FEPO4',I,J,L,M,PFEPO1,AFE1,AFEO1,AFEO2,AFEO3,AFEO4
!    2,AH0P1,AH1P1,AH2P1,AHY1,AOH1,RPFEPX,RHF0P1,RHF1P1,RHF2P1,RHF3P1
!    3,RHF4P1,RHF0P2,RHF1P2,RHF2P2,RHF3P2,RHF4P2,SP,SPX,AFE1*AH0P1
!    4,SPFEP
!     ENDIF
!
!     DICALCIUM PHOSPHATE
!
    PX=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,AH1P1))THEN
      P2=AH1P1
      NR1=0
      SP=SPCAD
    ELSEIF(isclose(PX,AH2P1))THEN
      P2=AH2P1
      NR1=1
      SP=SHCAD2
    ENDIF
    RPCAD1=0._r8
    RHCAD2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2-4.0*(P1*P2-SPX))
    RPCADX=AMAX1(-PCAPD1,TPDX*(S0-SQRT(S1)))
    IF(isclose(PX,AH1P1))THEN
      RPCAD1=RPCADX
    ELSEIF(isclose(PX,AH2P1))THEN
      RHCAD2=RPCADX
    ENDIF
!     IF((M/10)*10.EQ.M)THEN
!     WRITE(*,1112)'CAPO4',I,J,L,M,PCAPM1,PCAPD1,CCA1
!    2,CH1P1,CH2P1,CHY1,COH1,RPCADX,RPCAD1,RHCAD2,R1,P1,P2,P3
!    3,SP,Z,FX,Y,X,TX,A2,CCA1*A2*CH1P1*A2,SPCAD
!     ENDIF
!
!     HYDROXYAPATITE
!
    PX=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,AH1P1))THEN
      P2=AH1P1
      NR1=4
      SP=SHCAH1
    ELSEIF(isclose(PX,AH2P1))THEN
      P2=AH2P1
      NR1=7
      SP=SHCAH2
    ENDIF
    RHCAH1=0._r8
    RHCAH2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=(SP*R1**NR1/P1**5)**0.333
    RPCAHX=AMAX1(-PCAPH1,TPDX*(P2-SPX))
    IF(isclose(PX,AH1P1))THEN
      RHCAH1=RPCAHX
    ELSEIF(isclose(PX,AH2P1))THEN
      RHCAH2=RPCAHX
    ENDIF
!     IF((I/10)*10.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,1112)'A1',I,L,K,M,A1,A2,A3,FSTR2,CSTR1
!    2,CSTR2,CC3,CA3,CC2,CA2,CC1,CA1,VOLWM
!     WRITE(*,1112)'APATITE',I,J,L,M,PCAPH1,ACA1,XCA1
!    2,AH0P1,AH1P1,AH2P1,AHY1,AOH1,RPCAHX,RHCAH1,RHCAH2
!    3,SP,SPX,ACA1**5*AH0P1**3*AOH1,SPCAH,SHCAH1,SHCAH2
!    3,CH0P1,CH1P1,CH2P1,XOH01,XOH11,XOH21,XH1P1,XH2P1
!    4,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    2,RHA4P1,RHF0P1,RHF1P1,RHF2P1
!    3,RHF3P1,RHF4P1,RPCAD1,3.0*RHCAH1
!    4,RXH1P,RH1P,RH2P,RF1P,RC1P,RM1P
!    5,RHA0P2,RHA1P2,RHA2P2,RHA3P2
!    2,RHA4P2,RHF0P2,RHF1P2,RHF2P2
!    3,RHF3P2,RHF4P2,RHCAD2,3.0*RHCAH2
!    4,RXH2P,RYH2P,RH2P,RH3P,RF2P,RC2P,RH3P
!     ENDIF
!
!     MONOCALCIUM PHOSPHATE
!
    P1=ACA1
    P2=AH2P1
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SPCAM
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2-4.0*(P1*P2-SPX))
    RPCAMX=AMAX1(-PCAPM1,TPDX*(S0-SQRT(S1)))
  ELSE
    RPALPX=0._r8
    RPFEPX=0._r8
    RPCADX=0._r8
    RPCAHX=0._r8
    RHA0P1=0._r8
    RHA1P1=0._r8
    RHA2P1=0._r8
    RHA3P1=0._r8
    RHA4P1=0._r8
    RHA0P2=0._r8
    RHA1P2=0._r8
    RHA2P2=0._r8
    RHA3P2=0._r8
    RHA4P2=0._r8
    RHF0P1=0._r8
    RHF1P1=0._r8
    RHF2P1=0._r8
    RHF3P1=0._r8
    RHF4P1=0._r8
    RHF0P2=0._r8
    RHF1P2=0._r8
    RHF2P2=0._r8
    RHF3P2=0._r8
    RHF4P2=0._r8
    RPCAD1=0._r8
    RHCAD2=0._r8
    RHCAH1=0._r8
    RHCAH2=0._r8
    RPCAMX=0._r8
  ENDIF
  end subroutine PhospPrecipDissolNonBand

!----------------------------------------------------------------------------
  subroutine SummarizeIonFluxes

  implicit none
!     begin_execution
!     CONVERT TOTAL ION FLUXES FROM CHANGES IN CONCENTRATION
!     TO CHANGES IN MASS PER UNIT AREA FOR USE IN 'REDIST'
!
  TRN4S=TRN4S*VOLWNH
  TRN4B=TRN4B*VOLWNB
  TRN3S=TRN3S*VOLWNH
  TRN3B=TRN3B*VOLWNB
  TRAL=TRAL*VOLWM
  TRFE=TRFE*VOLWM
  TRHY=TRHY*VOLWM
  TRCA=TRCA*VOLWM
  TRMG=TRMG*VOLWM
  TRNA=TRNA*VOLWM
  TRKA=TRKA*VOLWM
  TROH=TROH*VOLWM
  TRSO4=TRSO4*VOLWM
  TRCO3=TRCO3*VOLWM
  TRHCO=TRHCO*VOLWM
  TRCO2=TRCO2*VOLWM
  TRAL1=TRAL1*VOLWM
  TRAL2=TRAL2*VOLWM
  TRAL3=TRAL3*VOLWM
  TRAL4=TRAL4*VOLWM
  TRALS=TRALS*VOLWM
  TRFE1=TRFE1*VOLWM
  TRFE2=TRFE2*VOLWM
  TRFE3=TRFE3*VOLWM
  TRFE4=TRFE4*VOLWM
  TRFES=TRFES*VOLWM
  TRCAO=TRCAO*VOLWM
  TRCAC=TRCAC*VOLWM
  TRCAH=TRCAH*VOLWM
  TRCAS=TRCAS*VOLWM
  TRMGO=TRMGO*VOLWM
  TRMGC=TRMGC*VOLWM
  TRMGH=TRMGH*VOLWM
  TRMGS=TRMGS*VOLWM
  TRNAC=TRNAC*VOLWM
  TRNAS=TRNAS*VOLWM
  TRKAS=TRKAS*VOLWM
  TRH0P=TRH0P*VOLWPO
  TRH1P=TRH1P*VOLWPO
  TRH2P=TRH2P*VOLWPO
  TRH3P=TRH3P*VOLWPO
  TRF1P=TRF1P*VOLWPO
  TRF2P=TRF2P*VOLWPO
  TRC0P=TRC0P*VOLWPO
  TRC1P=TRC1P*VOLWPO
  TRC2P=TRC2P*VOLWPO
  TRM1P=TRM1P*VOLWPO
  TRH0B=TRH0B*VOLWPB
  TRH1B=TRH1B*VOLWPB
  TRH2B=TRH2B*VOLWPB
  TRH3B=TRH3B*VOLWPB
  TRF1B=TRF1B*VOLWPB
  TRF2B=TRF2B*VOLWPB
  TRC0B=TRC0B*VOLWPB
  TRC1B=TRC1B*VOLWPB
  TRC2B=TRC2B*VOLWPB
  TRM1B=TRM1B*VOLWPB
  TRXN4=TRXN4*VOLWNH
  TRXNB=TRXNB*VOLWNB
  TRXHY=TRXHY*VOLWM
  TRXAL=TRXAL*VOLWM
  TRXFE=TRXFE*VOLWM
  TRXCA=TRXCA*VOLWM
  TRXMG=TRXMG*VOLWM
  TRXNA=TRXNA*VOLWM
  TRXKA=TRXKA*VOLWM
  TRXHC=TRXHC*VOLWM
  TRXAL2=TRXAL2*VOLWM
  TRXFE2=TRXFE2*VOLWM
  TRXH0=TRXH0*VOLWPO
  TRXH1=TRXH1*VOLWPO
  TRXH2=TRXH2*VOLWPO
  TRX1P=TRX1P*VOLWPO
  TRX2P=TRX2P*VOLWPO
  TRBH0=TRBH0*VOLWPB
  TRBH1=TRBH1*VOLWPB
  TRBH2=TRBH2*VOLWPB
  TRB1P=TRB1P*VOLWPB
  TRB2P=TRB2P*VOLWPB
  TRALOH=TRALOH*VOLWM
  TRFEOH=TRFEOH*VOLWM
  TRCACO=TRCACO*VOLWM
  TRCASO=TRCASO*VOLWM
  TRALPO=TRALPO*VOLWPO
  TRFEPO=TRFEPO*VOLWPO
  TRCAPD=TRCAPD*VOLWPO
  TRCAPH=TRCAPH*VOLWPO
  TRCAPM=TRCAPM*VOLWPO
  TRALPB=TRALPB*VOLWPB
  TRFEPB=TRFEPB*VOLWPB
  TRCPDB=TRCPDB*VOLWPB
  TRCPHB=TRCPHB*VOLWPB
  TRCPMB=TRCPMB*VOLWPB
!
!     BOUNDARY SALT FLUXES FOR C, H, OH, P, AL+FE, CA, OH
!     USED TO CHECK MATERIAL BALANCES IN REDIST.F
!
!     TBCO2=CO2 net change from all solute equilibria
!     TRH2O=H2O net change from all solute equilibria
!     TBION=total solute net change from all solute equilibria
!
  TBCO2=TRCO3+TRCAC+TRMGC+TRNAC+TRCACO &
    +2.0*(TRHCO+TRCAH+TRMGH)
  TRH2O=TRHY+TROH+TRXHY+TRXHC
  TBION=4.0*(TRH3P+TRH3B) &
    +3.0*(TRF2P+TRC2P &
    +TRF2B+TRC2B) &
    +2.0*(TRF1P+TRC1P+TRM1P &
    +TRF1B+TRC1B+TRM1B) &
    +TRH0P+TRC0P+TRH0B+TRC0B &
    -(TRALPO+TRFEPO &
    +TRALPB+TRFEPB) &
    -(TRCAPD+TRCAPM &
    +TRCPDB+TRCPMB &
    +5.0*(TRCAPH+TRCPHB)) &
    +TRAL1+TRFE1 &
    +2.0*(TRAL2+TRFE2) &
    +3.0*(TRAL3+TRFE3) &
    +4.0*(TRAL4+TRFE4) &
    +TRCAO+TRMGO &
    +3.0*(TRALOH+TRFEOH)
!     IF(L.EQ.11)THEN
!     WRITE(*,1111)'TRCO2',I,J,L,M,TRCO2,TRCO3
!    2,TRHCO,TRCAC,TRMGC
!    2,TRNAC,TRCAH
!    2,TRMGH,TRCACO,VOLWM,RCO2
!    3,RHCO,RHCACH,RCO2Q,RCAH,RMGH,RHCO3,AHY1,AHCO31,ACO21,DPCO2
!     WRITE(*,1111)'TBION',I,J,L,M,TBION
!     ENDIF
  end subroutine SummarizeIonFluxes
!------------------------------------------------------------------------------------------

  subroutine GetSoluteConcentrations
  implicit none
!     begin_execution

  CN41=AMAX1(ZERO,CN41)
  CN4B=AMAX1(ZERO,CN4B)
  CN31=AMAX1(ZERO,CN31)
  CN3B=AMAX1(ZERO,CN3B)
  CHY1=AMAX1(ZERO,10.0**(-(PH-3.0)))
  COH1=AMAX1(ZERO,DPH2O/CHY1)
  CCO31=AMAX1(ZERO,CCO31)
  CAL1=AMAX1(ZERO,CAL1)
  CFE1=AMAX1(ZERO,CFE1)
  CCA1=AMAX1(ZERO,AMIN1(CCAMX,CCA1))
  CMG1=AMAX1(ZERO,CMG1)
  CNA1=AMAX1(ZERO,CNA1)
  CKA1=AMAX1(ZERO,CKA1)
  CSO41=AMAX1(ZERO,CSO41)
  CHCO31=AMAX1(ZERO,CHCO31)
  CCO21=AMAX1(ZERO,CCO21)
  CALO1=AMAX1(ZERO,CALO1)
  CALO2=AMAX1(ZERO,CALO2)
  CALO3=AMAX1(ZERO,CALO3)
  CALO4=AMAX1(ZERO,CALO4)
  CALS1=AMAX1(ZERO,CALS1)
  CFEO1=AMAX1(ZERO,CFEO1)
  CFEO2=AMAX1(ZERO,CFEO2)
  CFEO3=AMAX1(ZERO,CFEO3)
  CFEO4=AMAX1(ZERO,CFEO4)
  CFES1=AMAX1(ZERO,CFES1)
  CCAO1=AMAX1(ZERO,CCAO1)
  CCAC1=AMAX1(ZERO,CCAC1)
  CCAH1=AMAX1(ZERO,CCAH1)
  CCAS1=AMAX1(ZERO,CCAS1)
  CMGO1=AMAX1(ZERO,CMGO1)
  CMGC1=AMAX1(ZERO,CMGC1)
  CMGH1=AMAX1(ZERO,CMGH1)
  CMGS1=AMAX1(ZERO,CMGS1)
  CNAC1=AMAX1(ZERO,CNAC1)
  CNAS1=AMAX1(ZERO,CNAS1)
  CKAS1=AMAX1(ZERO,CKAS1)
  CH0P1=AMAX1(ZERO,CH0P1)
  CH1P1=AMAX1(ZERO,CH1P1)
  CH2P1=AMAX1(ZERO,CH2P1)
  CH3P1=AMAX1(ZERO,CH3P1)
  CF1P1=AMAX1(ZERO,CF1P1)
  CF2P1=AMAX1(ZERO,CF2P1)
  CC0P1=AMAX1(ZERO,CC0P1)
  CC1P1=AMAX1(ZERO,CC1P1)
  CC2P1=AMAX1(ZERO,CC2P1)
  CM1P1=AMAX1(ZERO,CM1P1)
  CH0PB=AMAX1(ZERO,CH0PB)
  CH1PB=AMAX1(ZERO,CH1PB)
  CH2PB=AMAX1(ZERO,CH2PB)
  CH3PB=AMAX1(ZERO,CH3PB)
  CF1PB=AMAX1(ZERO,CF1PB)
  CF2PB=AMAX1(ZERO,CF2PB)
  CC0PB=AMAX1(ZERO,CC0PB)
  CC1PB=AMAX1(ZERO,CC1PB)
  CC2PB=AMAX1(ZERO,CC2PB)
  CM1PB=AMAX1(ZERO,CM1PB)
  XCOO=AMAX1(0._r8,XCOOH-XHC1-XALO21-XFEO21)
  end subroutine GetSoluteConcentrations
!------------------------------------------------------------------------------------------

  subroutine IonStrengthActivity
  implicit none
  real(r8) :: CC3,CA3,CC2,CA2,CC1,CA1,CSTR1,CSTR2
  real(r8) :: A1,A3,FSTR2
!     begin_execution
!     IONIC STRENGTH FROM SUMS OF ION CONCENTRATIONS
!
!     CC3,CA3,CC2,CA2,CC1,CA1=total tri-,di-,univalent cations C,anions A
!     CSTR1=ion strength
!
  CC3=CAL1+CFE1
  CA3=CH0P1*VLPO4+CH0PB*VLPOB
  CC2=CCA1+CMG1+CALO1+CFEO1+CF2P1*VLPO4+CF2PB*VLPOB
  CA2=CSO41+CCO31+CH1P1*VLPO4+CH1PB*VLPOB
  CC1=CN41*VLNH4+CN4B*VLNHB+CHY1+CNA1+CKA1 &
    +CALO2+CFEO2+CALS1+CFES1+CCAO1+CCAH1+CMGO1+CMGH1 &
    +(CF1P1+CC2P1)*VLPO4+(CF1PB+CC2PB)*VLPOB
  CA1=CNO1*VLNO3+CNOB*VLNOB+COH1+CHCO31+CCL1 &
    +CALO4+CFEO4+CNAC1+CNAS1+CKAS1+(CH2P1+CC0P1)*VLPO4 &
    +(CH2PB+CC0PB)*VLPOB
  CSTR1=AMAX1(0._r8,0.5E-03*(9.0*(CC3+CA3)+4.0*(CC2+CA2)+CC1+CA1))
  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0+CSTR2)
!
!     ACTIVITY COEFFICIENTS CALCULATED FROM ION STRENGTH
!
  A1=AMIN1(1.0,10.0**(-0.509*1.0*FSTR2+0.20*CSTR2))
  A2=AMIN1(1.0,10.0**(-0.509*4.0*FSTR2+0.20*CSTR2))
  A3=AMIN1(1.0,10.0**(-0.509*9.0*FSTR2+0.20*CSTR2))
!
!     PRECIPITATION-DISSOLUTION REACTIONS IN NON-BAND, BAND
!     CALCULATED FROM ACTIVITIES OF REACTANTS AND PRODUCTS THROUGH
!     CONVERGENCE SOLUTIONS FOR THEIR EQUILIBRIUM CONCENTRATIONS USING
!     SOLUTE FORMS CURRENTLY AT HIGHEST CONCENTRATIONS
!
!     for all reactions:
!     A*=ion activity
!     PX,PY=solute forms with greatest activity
!     R*,P*=reactant,product
!     NR*,NP*=reactant,product stoichiometry
!     SP=solubility product of PX from parameters above
!     SPX=equilibrium product concentration
!     R*X=precipitation-dissolution rate
!
  AHY1=CHY1*A1
  AOH1=COH1*A1
  AAL1=CAL1*A3
  AALO1=CALO1*A2
  AALO2=CALO2*A1
  AALO3=CALO3
  AALO4=CALO4*A1
  AFE1=CFE1*A3
  AFEO1=CFEO1*A2
  AFEO2=CFEO2*A1
  AFEO3=CFEO3
  AFEO4=CFEO4*A1
  ACA1=CCA1*A2
  ACO31=CCO31*A2
  AHCO31=CHCO31*A1
  ACO21=CCO21*A0
  ASO41=CSO41*A2
  AH0P1=CH0P1*A3
  AH1P1=CH1P1*A2
  AH2P1=CH2P1*A1
  AH3P1=CH3P1*A0
  AF1P1=CF1P1*A2
  AF2P1=CF2P1*A1
  AC0P1=CC0P1*A1
  AC1P1=CC1P1*A0
  AC2P1=CC2P1*A1
  AM1P1=CM1P1*A0
  AH0PB=CH0PB*A3
  AH1PB=CH1PB*A2
  AH2PB=CH2PB*A1
  AH3PB=CH3PB*A0
  AF1PB=CF1PB*A2
  AF2PB=CF2PB*A1
  AC0PB=CC0PB*A1
  AC1PB=CC1PB*A0
  AC2PB=CC2PB*A1
  AM1PB=CM1PB*A0
  AN41=CN41*A1
  AN4B=CN4B*A1
  AN31=CN31*A0
  AN3B=CN3B*A0
  AMG1=CMG1*A2
  ANA1=CNA1*A1
  AKA1=CKA1*A1
  AALS1=CALS1*A1
  AFES1=CFES1*A1
  ACAO1=CCAO1*A1
  ACAC1=CCAC1*A0
  ACAS1=CCAS1*A0
  ACAH1=CCAH1*A1
  AMGO1=CMGO1*A1
  AMGC1=CMGC1*A0
  AMGH1=CMGH1*A1
  AMGS1=CMGS1*A0
  ANAC1=CNAC1*A1
  ANAS1=CNAS1*A1
  AKAS1=CKAS1*A1
  end subroutine IonStrengthActivity
!------------------------------------------------------------------------------------------

  subroutine PhospPrecipDissolBand

  implicit none
  real(r8) :: SP,SPX,S0,S1,P1,P2,P3,PX,PY
  integer :: NR1,NP3
!     begin_execution

  IF(VOLWPB.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
    PX=AMAX1(AAL1,AALO1,AALO2,AALO3,AALO4)
    PY=AMAX1(AH1PB,AH2PB)
    R1=AHY1
    P3=AHY1
    IF(isclose(PY,AH1PB))THEN
      P2=AH1PB
      IF(isclose(PX,AAL1))THEN
        P1=AAL1
        NR1=1
        NP3=0
        SP=SHA0P1
      ELSEIF(isclose(PX,AALO1))THEN
        P1=AALO1
        NR1=0
        NP3=0
        SP=SPA1P1
      ELSEIF(isclose(PX,AALO2))THEN
        P1=AALO2
        NR1=0
        NP3=1
        SP=SHA2P1
      ELSEIF(isclose(PX,AALO3))THEN
        P1=AALO3
        NR1=0
        NP3=2
        SP=SHA3P1
      ELSEIF(isclose(PX,AALO4))THEN
        P1=AALO4
        NR1=0
        NP3=3
        SP=SHA4P1
      ENDIF
    ELSE
      P2=AH2PB
      IF(isclose(PX,AAL1))THEN
        P1=AAL1
        NR1=2
        NP3=0
        SP=SHA0P2
      ELSEIF(isclose(PX,AALO1))THEN
        P1=AALO1
        NR1=1
        NP3=0
        SP=SHA1P2
      ELSEIF(isclose(PX,AALO2))THEN
        P1=AALO2
        NR1=0
        NP3=0
        SP=SPA2P2
      ELSEIF(isclose(PX,AALO3))THEN
        P1=AALO3
        NR1=0
        NP3=1
        SP=SHA3P2
      ELSEIF(isclose(PX,AALO4))THEN
        P1=AALO4
        NR1=0
        NP3=2
        SP=SHA4P2
      ENDIF
    ENDIF
    RHA0B1=0._r8
    RHA1B1=0._r8
    RHA2B1=0._r8
    RHA3B1=0._r8
    RHA4B1=0._r8
    RHA0B2=0._r8
    RHA1B2=0._r8
    RHA2B2=0._r8
    RHA3B2=0._r8
    RHA4B2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2-4.0*(P1*P2-SPX))
    RPALBX=AMAX1(-PALPOB,TPDX*(S0-SQRT(S1)))
    IF(isclose(PY,AH1PB))THEN
      IF(isclose(PX,AAL1))THEN
        RHA0B1=RPALBX
      ELSEIF(isclose(PX,AALO1))THEN
        RHA1B1=RPALBX
      ELSEIF(isclose(PX,AALO2))THEN
        RHA2B1=RPALBX
      ELSEIF(isclose(PX,AALO3))THEN
        RHA3B1=RPALBX
      ELSEIF(isclose(PX,AALO4))THEN
        RHA4B1=RPALBX
      ENDIF
    ELSE
      IF(isclose(PX,AAL1))THEN
        RHA0B2=RPALBX
      ELSEIF(isclose(PX,AALO1))THEN
        RHA1B2=RPALBX
      ELSEIF(isclose(PX,AALO2))THEN
        RHA2B2=RPALBX
      ELSEIF(isclose(PX,AALO3))THEN
        RHA3B2=RPALBX
      ELSEIF(isclose(PX,AALO4))THEN
        RHA4B2=RPALBX
      ENDIF
    ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
    PX=AMAX1(AFE1,AFEO1,AFEO2,AFEO3,AFEO4)
    PY=AMAX1(AH1PB,AH2PB)
    R1=AHY1
    P3=AHY1
    IF(isclose(PY,AH1PB))THEN
      P2=AH1PB
      IF(isclose(PX,AFE1))THEN
        P1=AFE1
        NR1=1
        NP3=0
        SP=SHF0P1
      ELSEIF(isclose(PX,AFEO1))THEN
        P1=AFEO1
        NR1=0
        NP3=0
        SP=SPF1P1
      ELSEIF(isclose(PX,AFEO2))THEN
        P1=AFEO2
        NR1=0
        NP3=1
        SP=SHF2P1
      ELSEIF(isclose(PX,AFEO3))THEN
        P1=AFEO3
        NR1=0
        NP3=2
        SP=SHF3P1
      ELSEIF(isclose(PX,AFEO4))THEN
        P1=AFEO4
        NR1=0
        NP3=3
        SP=SHF4P1
      ENDIF
    ELSE
      P2=AH2PB
      IF(isclose(PX,AFE1))THEN
        P1=AFE1
        NR1=2
        NP3=0
        SP=SHF0P2
      ELSEIF(isclose(PX,AFEO1))THEN
        P1=AFEO1
        NR1=1
        NP3=0
        SP=SHF1P2
      ELSEIF(isclose(PX,AFEO2))THEN
        P1=AFEO2
        NR1=0
        NP3=0
        SP=SPF2P2
      ELSEIF(isclose(PX,AFEO3))THEN
        P1=AFEO3
        NR1=0
        NP3=1
        SP=SHF3P2
      ELSEIF(isclose(PX,AFEO4))THEN
        P1=AFEO4
        NR1=0
        NP3=2
        SP=SHF4P2
      ENDIF
    ENDIF
    RHF0B1=0._r8
    RHF1B1=0._r8
    RHF2B1=0._r8
    RHF3B1=0._r8
    RHF4B1=0._r8
    RHF0B2=0._r8
    RHF1B2=0._r8
    RHF2B2=0._r8
    RHF3B2=0._r8
    RHF4B2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2._r8-4.0_r8*(P1*P2-SPX))
    RPFEBX=AMAX1(-PFEPOB,TPDX*(S0-SQRT(S1)))
    IF(isclose(PY,AH1PB))THEN
      IF(isclose(PX,AFE1))THEN
        RHF0B1=RPFEBX
      ELSEIF(isclose(PX,AFEO1))THEN
        RHF1B1=RPFEBX
      ELSEIF(isclose(PX,AFEO2))THEN
        RHF2B1=RPFEBX
      ELSEIF(isclose(PX,AFEO3))THEN
        RHF3B1=RPFEBX
      ELSEIF(isclose(PX,AFEO4))THEN
        RHF4B1=RPFEBX
      ENDIF
    ELSE
      IF(isclose(PX,AFE1))THEN
        RHF0B2=RPFEBX
      ELSEIF(isclose(PX,AFEO1))THEN
        RHF1B2=RPFEBX
      ELSEIF(isclose(PX,AFEO2))THEN
        RHF2B2=RPFEBX
      ELSEIF(isclose(PX,AFEO3))THEN
        RHF3B2=RPFEBX
      ELSEIF(isclose(PX,AFEO4))THEN
        RHF4B2=RPFEBX
      ENDIF
    ENDIF
!
!     DICALCIUM PHOSPHATE
!
    PX=AMAX1(AH1PB,AH2PB)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,AH1PB))THEN
      P2=AH1PB
      NR1=0
      SP=SPCAD
    ELSEIF(isclose(PX,AH2PB))THEN
      P2=AH2PB
      NR1=1
      SP=SHCAD2
    ENDIF
    RPCDB1=0._r8
    RHCDB2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2._r8-4.0_r8*(P1*P2-SPX))
    RPCDBX=AMAX1(-PCAPDB,TPDX*(S0-SQRT(S1)))
    IF(isclose(PX,AH1PB))THEN
      RPCDB1=RPCDBX
    ELSEIF(isclose(PX,AH2PB))THEN
      RHCDB2=RPCDBX
    ENDIF
!
!     HYDROXYAPATITE
!
    PX=AMAX1(AH1PB,AH2PB)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,AH1PB))THEN
      P2=AH1PB
      NR1=4
      SP=SHCAH1
    ELSEIF(isclose(PX,AH2PB))THEN
      P2=AH2PB
      NR1=7
      SP=SHCAH2
    ENDIF
    RHCHB1=0._r8
    RHCHB2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=(SP*R1**NR1/P1**5._r8)**0.333_r8
    RPCHBX=AMAX1(-PCAPHB,TPDX*(P2-SPX))
    IF(isclose(PX,AH1PB))THEN
      RHCHB1=RPCHBX
    ELSEIF(isclose(PX,AH2PB))THEN
      RHCHB2=RPCHBX
    ENDIF
!
!     MONOCALCIUM PHOSPHATE
!
    P1=ACA1
    P2=AH2PB
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SPCAM
    S0=P1+P2
    S1=AMAX1(0._r8,S0**2._r8-4.0_r8*(P1*P2-SPX))
    RPCMBX=AMAX1(-PCAPMB,TPDX*(S0-SQRT(S1)))
  ELSE
    RPALBX=0._r8
    RPFEBX=0._r8
    RPCDBX=0._r8
    RPCHBX=0._r8
    RPCMBX=0._r8
    RHA0B1=0._r8
    RHA1B1=0._r8
    RHA2B1=0._r8
    RHA3B1=0._r8
    RHA4B1=0._r8
    RHA0B2=0._r8
    RHA1B2=0._r8
    RHA2B2=0._r8
    RHA3B2=0._r8
    RHA4B2=0._r8
    RHF0B1=0._r8
    RHF1B1=0._r8
    RHF2B1=0._r8
    RHF3B1=0._r8
    RHF4B1=0._r8
    RHF0B2=0._r8
    RHF1B2=0._r8
    RHF2B2=0._r8
    RHF3B2=0._r8
    RHF4B2=0._r8
    RPCDB1=0._r8
    RHCDB2=0._r8
    RHCHB1=0._r8
    RHCHB2=0._r8
  ENDIF
  end subroutine PhospPrecipDissolBand
!------------------------------------------------------------------------------------------

  subroutine PhospAnionExchNoBand
  implicit none
  real(r8) :: SPH1P,SPH2P
!     begin_execution
!     PHOSPHORUS ANION EXCHANGE IN NON-BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!
!     BKVL=soil mass
!     VOLWM=soil water volume
!     XAEC=anion exchange capacity
!     VOLWPO=soil water volume in non-band
!     TADAX=adsorption rate constant
!     RXOH2,RXOH1=OH2,OH exchange with R-OH2,R-OH in non-band
!     SXOH2,SXOH1=equilibrium constant for OH2,OH exchange with R-OH2,R-OH
!
  IF(VOLWM.GT.ZEROS2)THEN
    VOLWBK=AMIN1(1.0,BKVL/VOLWM)
  ELSE
    VOLWBK=1._r8
  ENDIF
  IF(VOLWPO.GT.ZEROS2.AND.XAEC.GT.ZEROS)THEN
    RXOH2=TADAX*(XOH11*AHY1-SXOH2*XOH21)/(XOH11+SXOH2)*VOLWBK
    RXOH1=TADAX*(XOH01*AHY1-SXOH1*XOH11)/(XOH01+SXOH1)*VOLWBK
!
!     H2PO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     RXH2P,RYH2P=H2PO4 exchange with R-OH2,R-OH in non-band
!
    SPH2P=SXH2P*DPH2O
    RXH2P=TADAX*(XOH21*AH2P1-SPH2P*XH2P1)/(XOH21+SPH2P)*VOLWBK
    RYH2P=TADAX*(XOH11*AH2P1-SXH2P*XH2P1*AOH1)/(XOH11+SXH2P)*VOLWBK
!
!     HPO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     RXH1P=HPO4 exchange with R-OH in non-band
!
    SPH1P=SXH1P*DPH2O/DPH2P
    RXH1P=TADAX*(XOH11*AH1P1-SPH1P*XH1P1)/(XOH11+SPH1P)*VOLWBK
  ELSE
    RXOH2=0._r8
    RXOH1=0._r8
    RXH2P=0._r8
    RYH2P=0._r8
    RXH1P=0._r8
  ENDIF
  end subroutine PhospAnionExchNoBand
!------------------------------------------------------------------------------------------

  subroutine PhospAnionExchBand
  implicit none
  real(r8) :: SPH1P,SPH2P
!     begin_execution
!     PHOSPHORUS ANION EXCHANGE IN BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!
  IF(VOLWPB.GT.ZEROS2.AND.XAEC.GT.ZEROS)THEN
!
!     RXO2B,RXO1B=OH2,OH exchange with R-OH2,R-OH in band
!     SXOH2,SXOH1=equilibrium constant for OH2,OH exchange with R-OH2,R-OH
!
    RXO2B=TADAX*(XH11B*AHY1-SXOH2*XH21B)/(XH11B+SXOH2)*VOLWBK
    RXO1B=TADAX*(XH01B*AHY1-SXOH1*XH11B)/(XH01B+SXOH1)*VOLWBK
!
!     H2PO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     RXH2B,RYH2B=H2PO4 exchange with R-OH2,R-OH in band
!
    SPH2P=SXH2P*DPH2O
    RXH2B=TADAX*(XH21B*AH2PB-SPH2P*X2P1B)/(XH21B+SPH2P)*VOLWBK
    RYH2B=TADAX*(XH11B*AH2PB-SXH2P*X2P1B*AOH1)/(XH11B+SXH2P)*VOLWBK
!
!     HPO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     RXH1B=HPO4 exchange with R-OH in band
!
    SPH1P=SXH1P*DPH2O/DPH2P
    RXH1B=TADAX*(XH11B*AH1PB-SPH1P*X1P1B)/(XH11B+SPH1P)*VOLWBK
!     WRITE(*,2226)'RXH1B',I,J,L,M,RXH1B,XH11B,XH21B,CH1PB
!    2,SPH1P,X1P1B,RYH2B,CH2PB,SXH2P,X2P1B,COH1,CHY1,ROH
!2226  FORMAT(A8,4I4,20E12.4)
  ELSE
    RXO2B=0._r8
    RXO1B=0._r8
    RXH2B=0._r8
    RYH2B=0._r8
    RXH1B=0._r8
  ENDIF
  end subroutine PhospAnionExchBand
!------------------------------------------------------------------------------------------

  subroutine CationExchange
  implicit none
  real(r8) :: FX,XALQ,XCAQ,XCAXl,XKAQ,XCAX,XFEQ
  real(r8) :: XHYQ,XMGQ,XN4Q,XNBQ,XTLQ,XNAQ
!     begin_execution

  IF(XCEC.GT.ZEROS)THEN
!
!     CATION CONCENTRATIONS
!
!     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS
!     AND CATION CONCENTRATIONS
!
!     CCEC,XCEC=cation exchange concentration,capacity
!     XCAX=equilibrium R-Ca concentration
!     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
!     Ca-NH4,Ca-H,Ca-Al,Ca-Mg,Ca-Na,Ca-K
!     X*Q=equilibrium exchangeable concentrations
!     XTLQ=total equilibrium exchangeable concentration
!
    AALX=AAL1**0.333_r8
    AFEX=AFE1**0.333_r8
    ACAX=ACA1**0.500_r8
    AMGX=AMG1**0.500_r8
    XCAX=CCEC/(1.0+GKC4*AN41/ACAX*VLNH4 &
      +GKC4*AN4B/ACAX*VLNHB &
      +GKCH*AHY1/ACAX+GKCA*AALX/ACAX &
      +GKCA*AFEX/ACAX+GKCM*AMGX/ACAX &
      +GKCN*ANA1/ACAX+GKCK*AKA1/ACAX)
    XN4Q=XCAX*AN41*GKC4
    XNBQ=XCAX*AN4B*GKC4
    XHYQ=XCAX*AHY1*GKCH
    XALQ=XCAX*AALX*GKCA
    XFEQ=XCAX*AFEX*GKCA
    XCAQ=XCAX*ACAX
    XMGQ=XCAX*AMGX*GKCM
    XNAQ=XCAX*ANA1*GKCN
    XKAQ=XCAX*AKA1*GKCK
    XTLQ=XN4Q*VLNH4+XNBQ*VLNHB+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CCEC/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XN4Q=FX*XN4Q
    XNBQ=FX*XNBQ
    XHYQ=FX*XHYQ
    XALQ=FX*XALQ/3.0_r8
    XFEQ=FX*XFEQ/3.0_r8
    XCAQ=FX*XCAQ/2.0_r8
    XMGQ=FX*XMGQ/2.0_r8
    XNAQ=FX*XNAQ
    XKAQ=FX*XKAQ
!
!     NH4 EXCHANGE IN NON-BAND AND BAND SOIL ZONES
!
!     RXN4,RXNB=NH4 adsorption in non-band,band
!     TADCX=adsorption rate constant
!
    RXN4=TADCX*AMAX1(AMIN1((XN4Q-XN41)*AN41/XN4Q,CN41),-XN41)
    RXNB=TADCX*AMAX1(AMIN1((XNBQ-XN4B)*AN4B/XNBQ,CN4B),-XN4B)
!
!     H,AL,FE,CA,MG,NA,K EXCHANGE
!
!     RX*=ion adsorption
!
    RXHY=TADCX*AMIN1((XHYQ-XHY1)*AHY1/XHYQ,CHY1)
    RXAL=TADCX*AMIN1((XALQ-XAL1)*AALX/XALQ,CAL1)
    RXFE=TADCX*AMIN1((XFEQ-XFE1)*AFEX/XFEQ,CFE1)
    RXCA=TADCX*AMIN1((XCAQ-XCA1)*ACAX/XCAQ,CCA1)
    RXMG=TADCX*AMIN1((XMGQ-XMG1)*AMGX/XMGQ,CMG1)
    RXNA=TADCX*AMIN1((XNAQ-XNA1)*ANA1/XNAQ,CNA1)
    RXKA=TADCX*AMIN1((XKAQ-XKA1)*AKA1/XKAQ,CKA1)
!     IF(I.EQ.256.AND.L.EQ.1)THEN
!     WRITE(*,1112)'RXAL',I,J,L,M,RXAL,TADCX,XALQ,XAL1,CAL1
!    2,AAL1,AALX,CCEC,XCAX,GKCA,FX
!     ENDIF
  ELSE
    RXN4=0._r8
    RXNB=0._r8
    RXHY=0._r8
    RXAL=0._r8
    RXFE=0._r8
    RXCA=0._r8
    RXMG=0._r8
    RXNA=0._r8
    RXKA=0._r8
  ENDIF
  end subroutine CationExchange
!------------------------------------------------------------------------------------------

  subroutine SoluteDissociation
  implicit none
  real(r8) :: S0,S1
!     begin_execution
!     for all reactions:
!     S0,S1=equilibrium constant,equilibrium solute concentration**2
!     TSLX=dissociation rate constant
!
!     DISSOCIATION OF CARBOXYL RADICALS
!     AND ADSORPTION OF AL AND FE OH2
!
!     RXHC=COOH-COO+H dissociation
!     RXALO2,RXFLO2=Al(OH2)2+,FeOH2+ adsorption
! COOH <-> COO(-) + H(+)
  S0=AHY1+XCOO+DPCOH
  S1=AMAX1(0._r8,S0**2-4.0_r8*(AHY1*XCOO-DPCOH*XHC1))
  RXHC=TADCX*(S0-SQRT(S1))
  S0=AALO2+XCOO+DPALO
  S1=AMAX1(0._r8,S0**2-4.0_r8*(AALO2*XCOO-DPALO*XALO21))
  RXALO2=TADAX*(S0-SQRT(S1))
  S0=AFEO2+XCOO+DPFEO
  S1=AMAX1(0._r8,S0**2-4.0_r8*(AFEO2*XCOO-DPFEO*XFEO21))
  RXFEO2=TADAX*(S0-SQRT(S1))
!
!     RNH4,RNHB-NH4-NH3+H dissociation in non-band,band
!     DPN4=NH4 dissociation constant
! NH4(+) <-> NH3 + H(+)
  IF(VOLWNH.GT.ZEROS2)THEN
    RNH4=TSLX*(AHY1*AN31-DPN4*AN41)/(DPN4+AHY1)
  ELSE
    RNH4=0._r8
  ENDIF
  IF(VOLWNB.GT.ZEROS2)THEN
    RNHB=TSLX*(AHY1*AN3B-DPN4*AN4B)/(DPN4+AHY1)
  ELSE
    RNHB=0._r8
  ENDIF
!
! RCO2Q=CO2-HCO3+H dissociation
! H2CO3 <-> HCO3(-) + H(+)

  S0=AHY1+AHCO31+DPCO2
  S1=AMAX1(0._r8,S0**2-4.0_r8*(AHY1*AHCO31-DPCO2*ACO21))
  RCO2Q=TSLX*(S0-SQRT(S1))
!
!     RHCO3=HCO3-CO3+H dissociation
! HCO3(-) <-> CO3(--)+H(+)
  S0=AHY1+ACO31+DPHCO
  S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AHY1*ACO31-DPHCO*AHCO31))
  RHCO3=TSLX*(S0-SQRT(S1))
!
!     RALO1=ALOH(++) <-> AL(+++)+OH(-) dissociation
!
  RALO1=TSLX*(AAL1*AOH1-DPAL1*AALO1)/(AOH1+DPAL1)
!
!     RALO2=ALOH2-ALOH+OH dissociation
! Al(OH)2(+) <-> Al(OH)(++)+OH(-)
  RALO2=TSLX*(AALO1*AOH1-DPAL2*AALO2)/(AOH1+DPAL2)
!
!     RALO3=ALOH3 <-> ALOH2+OH dissociation
! Al(OH)3 <-> Al(OH)2(+)+OH(-)
  RALO3=TSLX*(AALO2*AOH1-DPAL3*AALO3)/(AOH1+DPAL3)
!
!     RALO4=ALOH4 <-> ALOH3+OH dissociation
! Al(OH)4(-) <-> Al(OH)3+OH(-)
  RALO4=TSLX*(AALO3*AOH1-DPAL4*AALO4)/(AOH1+DPAL4)
!
!     RALS=ALSO4 <-> AL+SO4 dissociation
! AlSO4(+) <-> Al(3+) + SO4(2-)
  S0=AAL1+ASO41+DPALS
  S1=AMAX1(0._r8,S0**2-4.0*(AAL1*ASO41-DPALS*AALS1))
  RALS=TSLX*(S0-SQRT(S1))
!
!     RFEO1=FEOH <-> FE+OH dissociation
! Fe(OH)(++) <-> Fe(+++)+OH(-)
  RFEO1=TSLX*(AFE1*AOH1-DPFE1*AFEO1)/(AOH1+DPFE1)
!
!     RFEO2=FEOH2 <-> FEOH+OH dissociation
! Fe(OH)2(+) <-> Fe(OH)(++)+OH(-)
  RFEO2=TSLX*(AFEO1*AOH1-DPFE2*AFEO2)/(AOH1+DPFE2)
!
!     RFEO3=FEOH3 <-> FEOH2+OH dissociation
! Fe(OH)3 <-> Fe(OH)2(+)+OH(-)
  RFEO3=TSLX*(AFEO2*AOH1-DPFE3*AFEO3)/(AOH1+DPFE3)
!
!     RFE04=FEOH4 <-> FEOH3+OH dissociation
! Fe(OH)4(-) <-> Fe(OH)3 + OH(-)
  RFEO4=TSLX*(AFEO3*AOH1-DPFE4*AFEO4)/(AOH1+DPFE4)
!
!     RFES <-> FE+SO4 dissociation
! FeSO4(+) <-> Fe(+++)+SO4(--)
  S0=AFE1+ASO41+DPFES
  S1=AMAX1(0._r8,S0**2-4.0*(AFE1*ASO41-DPFES*AFES1))
  RFES=TSLX*(S0-SQRT(S1))
!
!     RCAO=CAOH <-> CA+OH dissociation
! Ca(OH)(+) <-> Ca(++)+OH(-)
  RCAO=TSLX*(ACA1*AOH1-DPCAO*ACAO1)/(AOH1+DPCAO)
!
!     RCAC=CACO3 <-> CA+CO3 dissociation
! CaCO3 <-> Ca(++)+CO3(--)
  S0=ACA1+ACO31+DPCAC
  S1=AMAX1(0._r8,S0**2-4.0*(ACA1*ACO31-DPCAC*ACAC1))
  RCAC=TSLX*(S0-SQRT(S1))
!
!     RCAH=CAHCO3<->CA+HCO3 dissociation
! CaHCO3(+) <-> Ca(++) + HCO3(-)
  S0=ACA1+AHCO31+DPCAH
  S1=AMAX1(0._r8,S0**2-4.0*(ACA1*AHCO31-DPCAH*ACAH1))
  RCAH=TSLX*(S0-SQRT(S1))
!
!     RCAS=CASO4<->CA+SO4 dissociation
! CaSO4 <-> Ca(++) + SO4(--)
  S0=ACA1+ASO41+DPCAS
  S1=AMAX1(0._r8,S0**2-4.0*(ACA1*ASO41-DPCAS*ACAS1))
  RCAS=TSLX*(S0-SQRT(S1))
!
!     RMGO=MGOH<->MG+OH dissociation
! MgOH(+) <-> Mg(++)+OH(-)
  RMGO=TSLX*(AMG1*AOH1-DPMGO*AMGO1)/(AOH1+DPMGO)
!
!     RMGC=MGCO3<->MG+CO3 dissociation
! MgCO3 <-> Mg(++)+CO3(--)
  S0=AMG1+ACO31+DPMGC
  S1=AMAX1(0._r8,S0**2-4.0*(AMG1*ACO31-DPMGC*AMGC1))
  RMGC=TSLX*(S0-SQRT(S1))
!
!     RMGH=MGHCO3<->MG+HCO3 dissociation
! MgHCO3(+) <-> Mg(++) + HCO3(-)
  S0=AMG1+AHCO31+DPMGH
  S1=AMAX1(0._r8,S0**2-4.0*(AMG1*AHCO31-DPMGH*AMGH1))
  RMGH=TSLX*(S0-SQRT(S1))
!
!     RMGS=MGSO4<->MG+SO4 dissociation
! MgSO4 <-> Mg(++) + SO4(--)
  S0=AMG1+ASO41+DPMGS
  S1=AMAX1(0._r8,S0**2-4.0*(AMG1*ASO41-DPMGS*AMGS1))
  RMGS=TSLX*(S0-SQRT(S1))
!
!     RNAC=NACO3<->NA+CO3 dissociation
! NaCO3(-) <-> Na(+) +CO3(--)
  S0=ANA1+ACO31+DPNAC
  S1=AMAX1(0._r8,S0**2-4.0*(ANA1*ACO31-DPNAC*ANAC1))
  RNAC=TSLX*(S0-SQRT(S1))
!
!     RNAS=NASO4<->NA+SO4 dissociation
! NaSO4(-) <-> Na(+)+SO4(--)
  S0=ANA1+ASO41+DPNAS
  S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(ANA1*ASO41-DPNAS*ANAS1))
  RNAS=TSLX*(S0-SQRT(S1))
!
!     RKAS=KSO4<->K+SO4 dissociation
! KSO4(-) <-> K(+)+SO4(--)
  S0=AKA1+ASO41+DPKAS
  S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AKA1*ASO41-DPKAS*AKAS1))
  RKAS=TSLX*(S0-SQRT(S1))
!
!     PHOSPHORUS IN NON-BAND SOIL ZONE
!
  IF(VOLWPO.GT.ZEROS2)THEN
!
!     RH1P=HPO4<->H+PO4 dissociation in non-band
! HPO4(--) <-> H(+)+PO4(---)
    RH1P=TSLX*(AH0P1*AHY1-DPH1P*AH1P1)/(DPH1P+AHY1)
!
!     RH2P=H2PO4 <-> H+HPO4 dissociation in non-band
! H2PO4(-) <-> H(+) + HPO4(--)
    RH2P=TSLX*(AH1P1*AHY1-DPH2P*AH2P1)/(DPH2P+AHY1)
!
!     RH3P=H3PO4 <-> H+H2PO4 dissociation in non-band
! H3PO4 <-> H(+)+H2PO4(-)
    RH3P=TSLX*(AH2P1*AHY1-DPH3P*AH3P1)/(DPH3P+AHY1)
!
!     RF1P=FEHPO4 <-> FE+HPO4 dissociation in non-band
! FeHPO4(+) <-> Fe(+++)+HPO4(--)
    S0=AFE1+AH1P1+DPF1P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AFE1*AH1P1-DPF1P*AF1P1))
    RF1P=TSLX*(S0-SQRT(S1))
!
!     RF2P=FEH2PO4 <-> FE+H2PO4 dissociation in non-band
! FeH2PO4(++) <-> Fe(+++)+H2PO4(-)
    S0=AFE1+AH2P1+DPF2P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AFE1*AH2P1-DPF2P*AF2P1))
    RF2P=TSLX*(S0-SQRT(S1))
!
!     RC0P=CAPO4 <-> CA+PO4 dissociation in non-band
! CaPO4(-) <-> Ca(++)+PO4(---)
    S0=ACA1+AH0P1+DPC0P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(ACA1*AH0P1-DPC0P*AC0P1))
    RC0P=TSLX*(S0-SQRT(S1))
!
!     RC1P=CAHPO4-CA+HPO4 dissociation in non-band
! CaHPO4 <-> Ca(++)+HPO4(--)

    S0=ACA1+AH1P1+DPC1P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(ACA1*AH1P1-DPC1P*AC1P1))
    RC1P=TSLX*(S0-SQRT(S1))
!
!     RC2P=CAH2PO4-CA+H2PO4 dissociation in non-band
! CaH2PO4(+) <-> Ca(++)+H2PO4(-)
    S0=ACA1+AH2P1+DPC2P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(ACA1*AH2P1-DPC2P*AC2P1))
    RC2P=TSLX*(S0-SQRT(S1))
!
!     RM1P=MGHPO4-MG+HPO4 dissociation in non-band
! MgHPO4 <-> Mg(++)+HPO4(--)
    S0=AMG1+AH1P1+DPM1P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AMG1*AH1P1-DPM1P*AM1P1))
    RM1P=TSLX*(S0-SQRT(S1))
  ELSE
    RH1P=0._r8
    RH2P=0._r8
    RH3P=0._r8
    RF1P=0._r8
    RF2P=0._r8
    RC0P=0._r8
    RC1P=0._r8
    RC2P=0._r8
    RM1P=0._r8
  ENDIF
!
!     PHOSPHORUS IN BAND SOIL ZONE
!
  IF(VOLWPB.GT.ZEROS2)THEN
!
!     RH1B=HPO4-H+PO4 dissociation in band
! HPO4(--) <-> H(+)+PO4(---)
    RH1B=TSLX*(AH0PB*AHY1-DPH1P*AH1PB)/(AHY1+DPH1P)
!
!     RH2B=H2PO4-H+HPO4 dissociation in band
! H2PO4(-) <-> H(+)+HPO4(--)
    RH2B=TSLX*(AH1PB*AHY1-DPH2P*AH2PB)/(AHY1+DPH2P)
!
!     RH3B=H3PO4-H+H2PO4 dissociation in band
! H3PO4 <-> H(+) + H2PO4(-)
    RH3B=TSLX*(AH2PB*AHY1-DPH3P*AH3PB)/(AHY1+DPH3P)
!
!     RF1B=FEHPO4-FE+HPO4 dissociation in band
! FeHPO4(+) <-> Fe(+++)+HPO4(--)
    S0=AFE1+AH1PB+DPF1P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AFE1*AH1PB-DPF1P*AF1PB))
    RF1B=TSLX*(S0-SQRT(S1))
!
!     RF2B=FEH2PO4-FE+H2PO4 dissociation in band
! FeH2PO4(+) <-> Fe(+++)+H2PO4(-)
    S0=AFE1+AH2PB+DPF2P
    S1=AMAX1(0._r8,S0**2_r8-4.0_r8*(AFE1*AH2PB-DPF2P*AF2PB))
    RF2B=TSLX*(S0-SQRT(S1))
!
!     RC0B=CAPO4-CA+PO4 dissociation in band
! CaPO4(-) <-> Ca(++)+PO4(---)
    S0=ACA1+AH0PB+DPC0P
    S1=AMAX1(0._r8,S0**2-4.0*(ACA1*AH0PB-DPC0P*AC0PB))
    RC0B=TSLX*(S0-SQRT(S1))
!
!     RC1B=CAHPO4-CA+HPO4 dissociation in band
! CaHPO4 <-> Ca(++)+HPO4(--)
    S0=ACA1+AH1PB+DPC1P
    S1=AMAX1(0._r8,S0**2-4.0*(ACA1*AH1PB-DPC1P*AC1PB))
    RC1B=TSLX*(S0-SQRT(S1))
!
!     RC2B=CAH2PO4-CA+H2PO4 dissociation in band
! CaH2PO4(+) <-> Ca(++)+H2PO4(-)
    S0=ACA1+AH2PB+DPC2P
    S1=AMAX1(0._r8,S0**2-4.0*(ACA1*AH2PB-DPC2P*AC2PB))
    RC2B=TSLX*(S0-SQRT(S1))
!
!     RM1B=MGHPO4-MG+HPO4 dissociation in band
! MgHPO4 <-> Mg(++)+HPO4(--)
    S0=AMG1+AH1PB+DPM1P
    S1=AMAX1(0._r8,S0**2-4.0*(AMG1*AH1PB-DPM1P*AM1PB))
    RM1B=TSLX*(S0-SQRT(S1))
  ELSE
    RH1B=0._r8
    RH2B=0._r8
    RH3B=0._r8
    RF1B=0._r8
    RF2B=0._r8
    RC0B=0._r8
    RC1B=0._r8
    RC2B=0._r8
    RM1B=0._r8
  ENDIF
  end subroutine SoluteDissociation
!------------------------------------------------------------------------------------------

  subroutine UpdateIonFluxCurentIter
  implicit none
!     begin_execution
!     RN4S,RN4B=net NH4 flux in non-band,band
!     RN3S,RN3B=net NH3 flux in non-band,band
!     RAL,RFE,RHY,RCA,RMG,RNA,RKA,ROH=net Al,Fe,H,Ca,Mg,Na,K,OH flux
!     RSO4,RCO3,RHCO,RCO2=net SO4,CO3,HCO3,CO2 flux
!     RAL1,RAL2,RAL3,RAL4,RALS=net AlOH,AlOH2,AlOH3,AlOH4,AlSO4
!     RFE1,RFE2,RFE3,RFE4,RFES=net FeOH,FeOH2,FeOH3,FeOH4,FeSO4
!     RHP0,RHP1,RHP2,RHP3=net PO4,HPO4,H2PO4,H3PO4 flux in non-band
!     RXH0,RXH1,RXH2,RX1P,RX2P=net R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 in non-band
!     RHB0,RHB1,RHB2,RHB3=net PO4,HPO4,H2PO4,H3PO4 flux in band
!     RBH0,RBH1,RBH2,RB1P,RB2P=net R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 in band
!
  RN4S=RNH4-RXN4
  RN4B=RNHB-RXNB
  RN3S=-RNH4
  RN3B=-RNHB
  RAL=-RHAL1-RXAL-RALO1-RALS-(RHA0P1+RHA0P2)*VLPO4 &
    -(RHA0B1+RHA0B2)*VLPOB
  RFE=-RHFE1-RXFE-RFEO1-RFES-(RHF0P1+RHF0P2+RF1P+RF2P)*VLPO4 &
    -(RHF0B1+RHF0B2+RF1B+RF2B)*VLPOB
  RHY=-RNH4*VLNH4-RNHB*VLNHB &
    -RXHY-RXHC+2.0*(RHALO1+RHFEO1+RHCACO &
    +(RHA0P2+RHF0P2-RHA3P1-RHA4P2-RHF3P1-RHF4P2)*VLPO4 &
    +(RHA0B2+RHF0B2-RHA3B1-RHA4B2-RHF3B1-RHF4B2)*VLPOB) &
    +3.0*(RHAL1+RHFE1-(RHA4P1+RHF4P1)*VLPO4 &
    -(RHF4B1+RHA4B1)*VLPOB)+4.0*(RHCAH1*VLPO4+RHCHB1*VLPOB) &
    +7.0*(RHCAH2*VLPO4+RHCHB2*VLPOB) &
    +RHALO2+RHFEO2-RHALO4-RHFEO4+RHCACH-RCO2Q-RHCO3 &
    +(RHA0P1-RHA2P1+RHA1P2-RHA3P2+RHF0P1-RHF2P1+RHF1P2-RHF3P2 &
    +RHCAD2-RXOH2-RXOH1-RH1P-RH2P-RH3P)*VLPO4 &
    +(RHA0B1-RHA2B1+RHA1B2-RHA3B2+RHF0B1-RHF2B1+RHF1B2-RHF3B2 &
    +RHCDB2-RXO2B-RXO1B-RH1B-RH2B-RH3B)*VLPOB
  RCA=-RPCACX-RPCASO-RXCA-RCAO-RCAC-RCAH-RCAS &
    -(RPCADX+RPCAMX+RC0P+RC1P+RC2P)*VLPO4 &
    -(RPCDBX+RPCMBX+RC0B+RC1B+RC2B)*VLPOB &
    -5.0*(RPCAHX*VLPO4+RPCHBX*VLPOB)
  RMG=-RXMG-RMGO-RMGC-RMGH-RMGS-RM1P*VLPO4-RM1B*VLPOB
  RNA=-RXNA-RNAC-RNAS
  RKA=-RXKA-RKAS
  ROH=-RCAO-RMGO-RALO1-RALO2-RALO3-RALO4-RFEO1-RFEO2-RFEO3-RFEO4 &
    -(-RYH2P-RXH1P)*VLPO4-(-RYH2B-RXH1B)*VLPOB
  RSO4=-RPCASO-RALS-RFES-RCAS-RMGS-RNAS-RKAS
  RCO3=-RHCAC3-RHCO3-RCAC-RMGC-RNAC
  RHCO=-RHCACH-RCO2Q-RCAH-RMGH+RHCO3
  RCO2=-RHCACO+RCO2Q
  RAL1=-RHALO1+RALO1-RALO2-(RHA1P1+RHA1P2)*VLPO4-(RHA1B1+RHA1B2)*VLPOB
  RAL2=-RHALO2+RALO2-RALO3-(RHA2P1+RHA2P2)*VLPO4-(RHA2B1+RHA2B2)*VLPOB
  RAL3=-RHALO3+RALO3-RALO4-(RHA3P1+RHA3P2)*VLPO4-(RHA3B1+RHA3B2)*VLPOB
  RAL4=-RHALO4+RALO4-(RHA4P1+RHA4P2)*VLPO4-(RHA4B1+RHA4B2)*VLPOB
  RFE1=-RHFEO1+RFEO1-RFEO2-(RHF1P1+RHF1P2)*VLPO4-(RHF1B1+RHF1B2)*VLPOB
  RFE2=-RHFEO2+RFEO2-RFEO3-(RHF2P1+RHF2P2)*VLPO4-(RHF2B1+RHF2B2)*VLPOB
  RFE3=-RHFEO3+RFEO3-RFEO4-(RHF3P1+RHF3P2)*VLPO4-(RHF3B1+RHF3B2)*VLPOB
  RFE4=-RHFEO4+RFEO4-(RHF4P1+RHF4P2)*VLPO4-(RHF4B1+RHF4B2)*VLPOB
  RHP0=-RH1P-RC0P
  RHP1=-RHA0P1-RHA1P1-RHA2P1-RHA3P1 &
    -RHA4P1-RHF0P1-RHF1P1-RHF2P1 &
    -RHF3P1-RHF4P1-RPCAD1-3.0*RHCAH1-RXH1P &
    +RH1P-RH2P-RF1P-RC1P-RM1P
  RHP2=-RHA0P2-RHA1P2-RHA2P2-RHA3P2 &
    -RHA4P2-RHF0P2-RHF1P2-RHF2P2 &
    -RHF3P2-RHF4P2-RHCAD2-3.0*RHCAH2 &
    -2.0*RPCAMX-RXH2P-RYH2P+RH2P-RH3P-RF2P-RC2P
  RHP3=RH3P
  RXH0=-RXOH1
  RXH1=RXOH1-RXOH2-RYH2P-RXH1P
  RXH2=RXOH2-RXH2P
  RX1P=RXH1P
  RX2P=RXH2P+RYH2P
  RHB0=-RH1B-RC0B
  RHB1=-RHA0B1-RHA1B1-RHA2B1-RHA3B1 &
    -RHA4B1-RHF0B1-RHF1B1-RHF2B1 &
    -RHF3B1-RHF4B1-RPCDB1-3.0*RHCHB1-RXH1B &
    +RH1B-RH2B-RF1B-RC1B-RM1B
  RHB2=-RHA0B2-RHA1B2-RHA2B2-RHA3B2 &
    -RHA4B2-RHF0B2-RHF1B2-RHF2B2 &
    -RHF3B2-RHF4B2-RHCDB2-3.0*RHCHB2 &
    -2.0*RPCMBX-RXH2B-RYH2B+RH2B-RH3B-RF2B-RC2B
  RHB3=RH3B
  RBH0=-RXO1B
  RBH1=RXO1B-RXO2B-RYH2B-RXH1B
  RBH2=RXO2B-RXH2B
  RB1P=RXH1B
  RB2P=RXH2B+RYH2B
!
  end subroutine UpdateIonFluxCurentIter
!------------------------------------------------------------------------------------------

  subroutine UpdateIonConcCurrentIter
  implicit none
  real(r8) :: CHY2,COH2
!     begin_execution
!     UPDATE ION CONCENTRATIONS FOR CURRENT ITERATION
!     FROM TOTAL ION FLUXES
!
  CN41=CN41+RN4S
  CN4B=CN4B+RN4B
  CN31=CN31+RN3S
  CN3B=CN3B+RN3B
  CAL1=CAL1+RAL
  CFE1=CFE1+RFE
  CHY1=CHY1+RHY
  CCA1=CCA1+RCA
  CMG1=CMG1+RMG
  CNA1=CNA1+RNA
  CKA1=CKA1+RKA
  COH1=COH1+ROH
  CSO41=CSO41+RSO4
  CCO31=CCO31+RCO3
  CHCO31=CHCO31+RHCO
  CCO21=CCO21+RCO2
  CALO1=CALO1+RAL1
  CALO2=CALO2+RAL2
  CALO3=CALO3+RAL3
  CALO4=CALO4+RAL4
  CALS1=CALS1+RALS
  CFEO1=CFEO1+RFE1
  CFEO2=CFEO2+RFE2
  CFEO3=CFEO3+RFE3
  CFEO4=CFEO4+RFE4
  CFES1=CFES1+RFES
  CCAO1=CCAO1+RCAO
  CCAC1=CCAC1+RCAC
  CCAH1=CCAH1+RCAH
  CCAS1=CCAS1+RCAS
  CMGO1=CMGO1+RMGO
  CMGC1=CMGC1+RMGC
  CMGH1=CMGH1+RMGH
  CMGS1=CMGS1+RMGS
  CNAC1=CNAC1+RNAC
  CNAS1=CNAS1+RNAS
  CKAS1=CKAS1+RKAS
  CH0P1=CH0P1+RHP0
  CH1P1=CH1P1+RHP1
  CH2P1=CH2P1+RHP2
  CH3P1=CH3P1+RHP3
  CF1P1=CF1P1+RF1P
  CF2P1=CF2P1+RF2P
  CC0P1=CC0P1+RC0P
  CC1P1=CC1P1+RC1P
  CC2P1=CC2P1+RC2P
  CM1P1=CM1P1+RM1P
  CH0PB=CH0PB+RHB0
  CH1PB=CH1PB+RHB1
  CH2PB=CH2PB+RHB2
  CH3PB=CH3PB+RHB3
  CF1PB=CF1PB+RF1B
  CF2PB=CF2PB+RF2B
  CC0PB=CC0PB+RC0B
  CC1PB=CC1PB+RC1B
  CC2PB=CC2PB+RC2B
  CM1PB=CM1PB+RM1B
!
!     RHHY,RHOH=H2O-H+OH equilibration
!
  CHY2=10.0**(-PH)*1.0E+03
  COH2=DPH2O/CHY2
  RHHY=CHY2-CHY1
  RHOH=COH2-COH1
  CHY1=CHY1+RHHY
  COH1=COH1+RHOH

!
!     UPDATE EXCHANGEABLE ION CONCENTRATIONS IN CURRENT
!     ITERATION FROM TOTAL ION FLUXES
!
  XN41=XN41+RXN4
  XN4B=XN4B+RXNB
  XHY1=XHY1+RXHY
  XAL1=XAL1+RXAL
  XFE1=XFE1+RXFE
  XCA1=XCA1+RXCA
  XMG1=XMG1+RXMG
  XNA1=XNA1+RXNA
  XKA1=XKA1+RXKA
  XHC1=XHC1+RXHC
  XALO21=XALO21+RXALO2
  XFEO21=XFEO21+RXFEO2
  XOH01=XOH01+RXH0
  XOH11=XOH11+RXH1
  XOH21=XOH21+RXH2
  XH1P1=XH1P1+RX1P
  XH2P1=XH2P1+RX2P
  XH01B=XH01B+RBH0
  XH11B=XH11B+RBH1
  XH21B=XH21B+RBH2
  X1P1B=X1P1B+RB1P
  X2P1B=X2P1B+RB2P
!
!     UPDATE PRECIPITATE CONCENTRATIONS IN CURRENT
!     ITERATION FROM TOTAL ION FLUXES
!
  PALOH1=PALOH1+RPALOX       !amorphous Al(OH)3 precipitation
  PFEOH1=PFEOH1+RPFEOX       !Fe(OH)3 precipitation
  PCACO1=PCACO1+RPCACX       !Calcite CaCO3 precipitation
  PCASO1=PCASO1+RPCASO       !Gypsum CaSO4 precipitation
  PALPO1=PALPO1+RPALPX       !Variscite AlPO4 precipitation non-band
  PFEPO1=PFEPO1+RPFEPX       !FePO4 precipitation
  PCAPD1=PCAPD1+RPCADX       !CaHPO4 precpitation
  PCAPH1=PCAPH1+RPCAHX       !Ca5(PO4)3OH (hydroxyapatite) precipitation
  PCAPM1=PCAPM1+RPCAMX       !Ca(H2PO4)2 precipitation non-band
  PALPOB=PALPOB+RPALBX       !AlPO4 precpitation band
  PFEPOB=PFEPOB+RPFEBX       !FePO4 precpitation band
  PCAPDB=PCAPDB+RPCDBX       !CaHPO4 precipitation band
  PCAPHB=PCAPHB+RPCHBX       !Ca5(PO4)3OH hydroxyapatite precpitation band
  PCAPMB=PCAPMB+RPCMBX       !Ca(H2PO4)2 precipitation band
  end subroutine UpdateIonConcCurrentIter
!------------------------------------------------------------------------------------------

  subroutine AccumulateIonFlux
  implicit none
!     begin_execution
!     ACCUMULATE TOTAL ION FLUXES FOR ALL ITERATIONS
!
!     TRN4S,TRN4B=total NH4 flux in non-band,band
!     TRN3S,TRN3B=total NH3 flux in non-band,band
!     TRAL,TRFE,TRHY,TRCA,TRMG,TRNA,TRKA,TROH=totalAl,Fe,H,Ca,Mg,Na,K,OH flux
!     TRSO4,TRCO3,TRHCO,TRCO2=total SO4,CO3,HCO3,CO2 flux
!     TRAL1,TRAL2,TRAL3,TRAL4,TRALS=total AlOH,AlOH2,AlOH3,AlOH4,AlSO4
!     TRFE1,TRFE2,TRFE3,TRFE4,TRFES=total FeOH,FeOH2,FeOH3,FeOH4,FeSO4
!     TRCAO,TRCAC,TRCAH,TRCAS=total CaOH,CaCO3,CaHCO3,CaSO4 flux
!     TRMGO,TRMGC,TRMGH,TRMGS=total MgOH,MgCO3,MgHCO3,MgSO4 flux
!     TRNAC,TRNAS,TRKAS=total NaCO3,NaSO4,KSO4 flux
!     TRH0P,TRH1P,TRH2P,TRH3P=net PO4,HPO4,H2PO4,H3PO4 flux in non-band
!     TRF1P,TRF2P,TRC0P,TRC1P,TRC2P,TRM1P
!     =total FeHPO4,FeH2PO4,CaPO4,CaHPO4,CaH2PO4,MgHPO4 in non-band
!     TRH0B,TRH1B,TRH2B,TRH3B=net PO4,HPO4,H2PO4,H3PO4 flux in band
!     TRF1B,TRF2B,TRC0B,TRC1B,TRC2B,TRM1B
!     =total FeHPO4,FeH2PO4,CaPO4,CaHPO4,CaH2PO4,MgHPO4 in band
!     TRNX4,TRNXB=total NH4 adsorption in non-band,band
!     TRXHY,TRXAL,TRXFE,TRXCA,TRXMG,TRXNA,TRXKA=total H,Al,Fe,Ca,Mg,Na,K adsorpn
!     TRXHC,TRXAL2,TRXFE2=total HCO3,AlOH2,FeOH2 adsorption
!     TRXH0,TRXH1,TRXH2,TRX1P,TRX2P
!     =total R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in non-band
!     TRBH0,TRBH1,TRBH2,TRB1P,TRB2P
!     =total R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in band
!     TRALOH,TRFEOH,TRCACO,TRCASO=total AlOH3,FeOH3,CaCO3,CaSO4 precipitation
!     TRALPO,TRFEPO,TRCAPD,TRCAPH,TRCAPM
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation in non-band
!     TRALPB,TRFEPB,TRCPDB,TRCPHB,TRCPMB
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation in band
!
  TRN4S=TRN4S+RN4S    !net NH4 flux in non-band
  TRN4B=TRN4B+RN4B    !net NH4 flux in band
  TRN3S=TRN3S+RN3S    !net NH3 flux in non-band
  TRN3B=TRN3B+RN3B    !net NH3 flux in band
  TRAL=TRAL+RAL       !total Al flux
  TRFE=TRFE+RFE       !total Fe(3+) flux
  TRHY=TRHY+RHY+RHHY  !total H(+) flux
  TRCA=TRCA+RCA       !total Ca(2+) flux
  TRMG=TRMG+RMG       !total Mg(2+) flux
  TRNA=TRNA+RNA       !total Na(+) flux
  TRKA=TRKA+RKA       !total K(+) flux
  TROH=TROH+ROH+RHOH  !total OH(-) flux
  TRSO4=TRSO4+RSO4    !total SO4(2-) flux
  TRCO3=TRCO3+RCO3    !total CO3(2-) flux
  TRHCO=TRHCO+RHCO    !total HCO3(2-) flux
  TRCO2=TRCO2+RCO2    !total CO2 flux due to dissociation
  TRAL1=TRAL1+RAL1    !total Al(OH)(2+) flux
  TRAL2=TRAL2+RAL2    !total Al(OH)2(+) flux
  TRAL3=TRAL3+RAL3    !total Al(OH)3 flux
  TRAL4=TRAL4+RAL4    !total Al(OH)4(-) flux
  TRALS=TRALS+RALS    !total Al(SO4)(+) flux
  TRFE1=TRFE1+RFE1    !total Fe(OH)(2+) flux
  TRFE2=TRFE2+RFE2    !total Fe(OH)2(-) flux
  TRFE3=TRFE3+RFE3    !total Fe(OH)3 flux
  TRFE4=TRFE4+RFE4    !total Fe(OH4) flux
  TRFES=TRFES+RFES    !total FeSO4(+) flux
  TRCAO=TRCAO+RCAO    !total Ca(OH)(+) flux
  TRCAC=TRCAC+RCAC    !total CaCO3 flux
  TRCAH=TRCAH+RCAH    !CaHCO3 flux
  TRCAS=TRCAS+RCAS    !CaSO4 flux
  TRMGO=TRMGO+RMGO    !MgOH  (+)
  TRMGC=TRMGC+RMGC    !MgCO3
  TRMGH=TRMGH+RMGH    !MgHCO3 (+)
  TRMGS=TRMGS+RMGS    !MgSO4
  TRNAC=TRNAC+RNAC    !NaCO3(-)
  TRNAS=TRNAS+RNAS    !NaSO4(-)
  TRKAS=TRKAS+RKAS    !KSO4(-)
  TRH0P=TRH0P+RHP0    !PO4 (3-)
  TRH1P=TRH1P+RHP1    !HPO4 (2-)
  TRH2P=TRH2P+RHP2    !H2PO4 (-)
  TRH3P=TRH3P+RHP3    !H3PO4
  TRF1P=TRF1P+RF1P    !FeHPO4
  TRF2P=TRF2P+RF2P    !FeH2PO4
  TRC0P=TRC0P+RC0P
  TRC1P=TRC1P+RC1P
  TRC2P=TRC2P+RC2P
  TRM1P=TRM1P+RM1P
  TRH0B=TRH0B+RHB0
  TRH1B=TRH1B+RHB1
  TRH2B=TRH2B+RHB2
  TRH3B=TRH3B+RHB3
  TRF1B=TRF1B+RF1B
  TRF2B=TRF2B+RF2B
  TRC0B=TRC0B+RC0B
  TRC1B=TRC1B+RC1B
  TRC2B=TRC2B+RC2B
  TRM1B=TRM1B+RM1B
  TRXN4=TRXN4+RXN4
  TRXNB=TRXNB+RXNB
  TRXHY=TRXHY+RXHY
  TRXAL=TRXAL+RXAL
  TRXFE=TRXFE+RXFE
  TRXCA=TRXCA+RXCA
  TRXMG=TRXMG+RXMG
  TRXNA=TRXNA+RXNA
  TRXKA=TRXKA+RXKA
  TRXHC=TRXHC+RXHC
  TRXAL2=TRXAL2+RXALO2
  TRXFE2=TRXFE2+RXFEO2
  TRXH0=TRXH0+RXH0
  TRXH1=TRXH1+RXH1
  TRXH2=TRXH2+RXH2
  TRX1P=TRX1P+RX1P
  TRX2P=TRX2P+RX2P
  TRBH0=TRBH0+RBH0
  TRBH1=TRBH1+RBH1
  TRBH2=TRBH2+RBH2
  TRB1P=TRB1P+RB1P
  TRB2P=TRB2P+RB2P
  TRALOH=TRALOH+RPALOX
  TRFEOH=TRFEOH+RPFEOX
  TRCACO=TRCACO+RPCACX
  TRCASO=TRCASO+RPCASO
  TRALPO=TRALPO+RPALPX
  TRFEPO=TRFEPO+RPFEPX
  TRCAPD=TRCAPD+RPCADX
  TRCAPH=TRCAPH+RPCAHX
  TRCAPM=TRCAPM+RPCAMX
  TRALPB=TRALPB+RPALBX
  TRFEPB=TRFEPB+RPFEBX
  TRCPDB=TRCPDB+RPCDBX
  TRCPHB=TRCPHB+RPCHBX
  TRCPMB=TRCPMB+RPCMBX
  end subroutine AccumulateIonFlux

end module SaltChemEquilibriaMod
