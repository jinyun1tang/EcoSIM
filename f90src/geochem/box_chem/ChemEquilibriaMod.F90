module ChemEquilibriaMod
!
! Description
! code to solve equilibrium chemistry, without the salt equilbrium
! the solver assumes prescirbed pH, and therefore [OH(-)]
! Note:
! the model does not include H3PO4, which seems problematic, July 28, 2022
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoluteParMod
  use EcosimConst
  use SoluteChemDataType
  use EcoSIMCtrlDataType
  use EcoSIMSolverPar
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__
  public ::  NoSaltChemEquilibria
  contains

!------------------------------------------------------------------
  subroutine NoSaltChemEquilibria(chemvar,solflx)
!
!  Solve for chemical concentrations using input pH
!     IF NO SALTS IS SELECTED IN SITE FILE THEN A SUBSET
!     OF THE EQUILIBRIA REACTIONS ARE SOLVED: MOSTLY THOSE
!     FOR PHOSPHORUS AND CO-REACTANTS
! Reference:
! Grant R.F. and Heaney D.J., 1997, Inorganic phosphorus transformation
! and transport in soils: mathematical modeling in ecosys, Soil Sci. Soc. Am. J.
! 61: 752-764.
  implicit none
  type(solute_flx_type),target, intent(inout) :: solflx
  type(chem_var_type), target, intent(inout) :: chemvar

  real(r8) :: CAL1,CALX,CCA1
  real(r8) :: CCAX,CCEC
  real(r8) :: CCO21,CCO31
  real(r8) :: CFEX,CFE1
  real(r8) :: CH2PA,CH2PD,CH2PF,CH2PH,CH2PM
  real(r8) :: CHY1,CNA1,CKA1,CMG1,CMGX
  real(r8) :: COH1,DP,FX
  real(r8) :: RB1P,RB2P,RBH1,RBH2,RH2B
  real(r8) :: RHB1,RHB2,RHP1,RHP2
  real(r8) :: RN3B,RN3S,RN4B,RN4S
  real(r8) :: RNHB,RPALBX,RPCDBX,RPCHBX
  real(r8) :: RPCMBX,RPFEBX,RX1P,RX2P
  real(r8) :: RXH1,RXH2,RXH1B,RXH2B,RXH1P,RXH2P,RYH2P
  real(r8) :: RXNB,RYH2B
  real(r8) :: S0,S1,SPH1P,SPH2P
  real(r8) :: XALQ,XCAQ,XCAX,XFEQ,XHYQ,XKAQ,XMGQ
  real(r8) :: XN4Q,XTLQ,XNAQ,XNBQ
  real(r8) :: RH2P,RNH4,RPALPX,RPCADX
  real(r8) :: RPCAHX,RPCAMX,RPFEPX,RXN4
  real(r8) :: VLWatMicPBK
  real(r8), parameter :: ZEROS = 1.0E-015_r8
  real(r8), parameter :: ZEROS2= 1.0E-008_r8

  real(r8), pointer :: SoilMicPMassLayerX
  real(r8), pointer :: VLWatMicPM
  real(r8), pointer :: SoilMicPMassLayer
  real(r8), pointer :: GKCA
  real(r8), pointer :: GKCH
  real(r8), pointer :: GKCM
  real(r8), pointer :: GKCK
  real(r8), pointer :: GKCN
  real(r8), pointer :: GKC4
  real(r8), pointer :: VLWatMicPNB
  real(r8), pointer :: VLWatMicPNH
  real(r8), pointer :: VLWatMicPPB
  real(r8), pointer :: VLWatMicPPO
  real(r8), pointer :: VLNH4
  real(r8), pointer :: VLNHB

  real(r8), pointer :: XCEC
  real(r8), pointer :: PH
  real(r8), pointer :: XAEC

  real(r8), pointer :: CAL
  real(r8), pointer :: CFE
  real(r8), pointer :: ZMG
  real(r8), pointer :: ZNA
  real(r8), pointer :: ZKA
  real(r8), pointer :: CO2S
  real(r8), pointer :: CCA
  real(r8), pointer :: CH1P1
  real(r8), pointer :: CH1PB
  real(r8), pointer :: CH2P1
  real(r8), pointer :: CH2PB
  real(r8), pointer :: CN31
  real(r8), pointer :: CN3B
  real(r8), pointer :: CN41
  real(r8), pointer :: CN4B
  real(r8), pointer :: XN41
  real(r8), pointer :: XN4B
  real(r8), pointer :: X1P1B
  real(r8), pointer :: X2P1B
  real(r8), pointer :: XH11B
  real(r8), pointer :: XH1P1
  real(r8), pointer :: XH21B
  real(r8), pointer :: XH2P1
  real(r8), pointer :: XOH11
  real(r8), pointer :: XOH21
  real(r8), pointer :: PALPO1
  real(r8), pointer :: PALPOB
  real(r8), pointer :: PCAPD1
  real(r8), pointer :: PCAPDB
  real(r8), pointer :: PCAPH1
  real(r8), pointer :: PCAPHB
  real(r8), pointer :: PCAPM1
  real(r8), pointer :: PCAPMB
  real(r8), pointer :: PFEPO1
  real(r8), pointer :: PFEPOB


  real(r8), pointer :: TRN4S
  real(r8), pointer :: TRN4B
  real(r8), pointer :: TRN3S
  real(r8), pointer :: TRN3B
  real(r8), pointer :: TRH1P
  real(r8), pointer :: TRH2P
  real(r8), pointer :: TRH1B
  real(r8), pointer :: TRH2B
  real(r8), pointer :: TRXN4
  real(r8), pointer :: TRXNB
  real(r8), pointer :: TRXH1
  real(r8), pointer :: TRXH2
  real(r8), pointer :: TRX1P
  real(r8), pointer :: TRX2P
  real(r8), pointer :: TRBH1
  real(r8), pointer :: TRBH2
  real(r8), pointer :: TRB1P
  real(r8), pointer :: TRB2P
  real(r8), pointer :: TRALPO
  real(r8), pointer :: TRFEPO
  real(r8), pointer :: TRCAPD
  real(r8), pointer :: TRCAPH
  real(r8), pointer :: TRCAPM
  real(r8), pointer :: TRALPB
  real(r8), pointer :: TRFEPB
  real(r8), pointer :: TRCPDB
  real(r8), pointer :: TRCPHB
  real(r8), pointer :: TRCPMB
  real(r8), pointer :: TRAL

!     begin_execution
  call solflx%SetZero()

  VLWatMicPNB => chemvar%VLWatMicPNB
  VLWatMicPNH => chemvar%VLWatMicPNH
  VLWatMicPPB => chemvar%VLWatMicPPB
  VLWatMicPPO => chemvar%VLWatMicPPO
  XOH11 => chemvar%XOH11
  XN41  => chemvar%XN41
  XN4B  => chemvar%XN4B
  CH1PB => chemvar%CH1PB
  CH1P1 => chemvar%CH1P1
  CH2P1 => chemvar%CH2P1
  CH2PB => chemvar%CH2PB
  X1P1B => chemvar%X1P1B
  X2P1B => chemvar%X2P1B
  XH11B => chemvar%XH11B
  XH1P1 => chemvar%XH1P1
  XH21B => chemvar%XH21B
  XH2P1 => chemvar%XH2P1
  XOH21 => chemvar%XOH21
  CN31  => chemvar%CN31
  CN3B  => chemvar%CN3B
  CN41  => chemvar%CN41
  CN4B  => chemvar%CN4B
  PALPO1=> chemvar%PALPO1
  PALPOB=> chemvar%PALPOB
  PCAPD1=> chemvar%PCAPD1
  PCAPDB=> chemvar%PCAPDB
  PCAPH1=> chemvar%PCAPH1
  PCAPHB=> chemvar%PCAPHB
  PCAPM1=> chemvar%PCAPM1
  PCAPMB=> chemvar%PCAPMB
  PFEPO1=> chemvar%PFEPO1
  PFEPOB=> chemvar%PFEPOB
  SoilMicPMassLayerX =>  chemvar%SoilMicPMassLayerX
  XCEC  =>  chemvar%XCEC
  PH    =>  chemvar%PH
  CAL   =>  chemvar%CAL
  CFE   =>  chemvar%CFE
  VLWatMicPM =>  chemvar%VLWatMicPM
  ZMG   =>  chemvar%ZMG
  ZNA   =>  chemvar%ZNA
  ZKA   =>  chemvar%ZKA
  CO2S  =>  chemvar%CO2S
  CCA   =>  chemvar%CCA
  SoilMicPMassLayer  =>  chemvar%SoilMicPMassLayer
  XAEC  =>  chemvar%XAEC
  VLNH4 =>  chemvar%VLNH4
  GKC4  =>  chemvar%GKC4
  VLNHB =>  chemvar%VLNHB
  GKCA  =>  chemvar%GKCA
  GKCH  =>  chemvar%GKCH
  GKCM  =>  chemvar%GKCM
  GKCK  =>  chemvar%GKCK
  GKCN  =>  chemvar%GKCN
  TRN4S  => solflx%TRN4S
  TRN4B  => solflx%TRN4B
  TRN3S  => solflx%TRN3S
  TRN3B  => solflx%TRN3B
  TRH1P  => solflx%TRH1P
  TRH2P  => solflx%TRH2P
  TRH1B  => solflx%TRH1B
  TRH2B  => solflx%TRH2B
  TRXN4  => solflx%TRXN4
  TRXNB  => solflx%TRXNB
  TRXH1  => solflx%TRXH1
  TRXH2  => solflx%TRXH2
  TRX1P  => solflx%TRX1P
  TRX2P  => solflx%TRX2P
  TRBH1  => solflx%TRBH1
  TRBH2  => solflx%TRBH2
  TRB1P  => solflx%TRB1P
  TRB2P  => solflx%TRB2P
  TRAL   => solflx%TRAL
  TRALPO => solflx%TRALPO
  TRFEPO => solflx%TRFEPO
  TRCAPD => solflx%TRCAPD
  TRCAPH => solflx%TRCAPH
  TRCAPM => solflx%TRCAPM
  TRALPB => solflx%TRALPB
  TRFEPB => solflx%TRFEPB
  TRCPDB => solflx%TRCPDB
  TRCPHB => solflx%TRCPHB
  TRCPMB => solflx%TRCPMB

!
!     PRECIPITATION-DISSOLUTION CALCULATED FROM ACTIVITIES
!     OF REACTANTS AND PRODUCTS THROUGH SOLUTIONS
!     FOR THEIR EQUILIBRIUM CONSTANTS USING CURRENT
!     ION CONCENTRATION
!
!     CCEC,XCEC=cation exchange concentration [mol/m3],capacity
!     BLVLX=soil mass
!     C*,Z*=solute concentration, mass
!     DP*=dissociation constant from PARAMETER above
!     SP*=solubility product from PARAMETER above
!     C*<0.0=solve for C* from equilibrium with other solutes
!  call PrintChemVar(chemvar)

  IF(SoilMicPMassLayerX.GT.ZEROS)THEN
    CCEC=AMAX1(ZERO,XCEC/SoilMicPMassLayerX)
  ELSE
    CCEC=ZERO
  ENDIF
  IF(VLWatMicPM.GT.ZEROS2)THEN
  ! aqueous Mg(2+)
    CMG1=AMAX1(ZERO,ZMG/VLWatMicPM)
  ! aqueous Na(+)
    CNA1=AMAX1(ZERO,ZNA/VLWatMicPM)
  ! aqueous K(-)
    CKA1=AMAX1(ZERO,ZKA/VLWatMicPM)

! aqueous CO2 (H2CO3), mol/m3
    CCO21=AMAX1(ZERO,CO2S/(Catomw*VLWatMicPM))

    VLWatMicPBK=AMIN1(1.0,SoilMicPMassLayer/VLWatMicPM)
  ELSE
    CMG1=0._r8
    CNA1=0._r8
    CKA1=0._r8
    CCO21=0._r8
    VLWatMicPBK=1._r8
  ENDIF


![H(+)] is prescribed
  CHY1=AMAX1(ZERO,10.0**(-(PH-3.0)))
! [OH(-)]
  COH1=AMAX1(ZERO,DPH2O/CHY1)
! Aluminum,[Al(3+)], Al(3+)+3OH(-)<->Al(OH)3
  IF(CAL.LT.0.0)THEN
    CAL1=AMAX1(ZERO,SPALO/COH1**3)
  ELSE
    CAL1=AMAX1(ZERO,AMIN1(CAL,SPALO/COH1**3))
  ENDIF
! Fe(3+) + 3OH(-)<->Fe(OH)3
  IF(CFE.LT.0.0)THEN
    CFE1=AMAX1(ZERO,SPFEO/COH1**3)
  ELSE
    CFE1=AMAX1(ZERO,AMIN1(CFE,SPFEO/COH1**3))
  ENDIF
!
!     CA CONCENTRATION FROM CURRENT CO2 CONCENTRATION
! CO3(2-)+2H(+) <->H2CO3
  CCO31=AMAX1(ZERO,CCO21*DPCO3/CHY1**2)
! [Ca(2+)], Ca(2+)+CO3(2-)<->CaCO3
  IF(CCA.LT.0.0_r8)THEN
    CCA1=AMAX1(ZERO,AMIN1(CCAMX,SPCAC/CCO31))
  ELSE
    CCA1=AMAX1(ZERO,AMIN1(CCA,SPCAC/CCO31))
  ENDIF

!
!     PHOSPHORUS TRANSFORMATIONS IN NON-BAND SOIL ZONE
!
  IF(VLWatMicPPO.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE, AlPO4)
!
!     CH2PA,CH2P1=equilibrium, current H2PO4 concentration in non-band
!     SYA0P2=solubility product derived from SPALO
!     RPALPX=H2PO4 dissolution from AlPO4 in non-band
!
!   PALP01: precipitated Al(PO4)
!   AlPO4+2H2O <-> 2OH(-)+H2PO4(-)+Al(3+)
    CH2PA=SYA0P2/(CAL1*COH1**2)
    RPALPX=AMAX1(-PALPO1,TPD*(CH2P1-CH2PA))
!
!     IRON PHOSPHATE (STRENGITE)
!
!     CH2PF,CH2P1=equilibrium,current H2PO4 concentration in non-band
!     SYF0P2=solubility product derived from SPALO
!     RPFEPX=H2PO4 dissolution from FePO4 in non-band
!
!   PFEP01: precipitated FePO4
!   FePO4+2H2O <-> Fe(3+)+2OH(-) + H2PO4(-)
    CH2PF=SYF0P2/(CFE1*COH1**2)
    RPFEPX=AMAX1(-PFEPO1,TPD*(CH2P1-CH2PF))
!
!     DICALCIUM PHOSPHATE (CaHPO4)
!
!     CH2PD,CH2P1=equilibrium,current H2PO4 concentration in non-band
!     SYCAD2=solubility product derived from SPALO
!     RPCADX=H2PO4 dissolution from CaHPO4 in non-band
!  PCAPD1: precipitated CaHPO4
!  CaHPO4 + H2O <-> Ca(2+) + OH(-)+H2PO4(-)
    CH2PD=SYCAD2/(CCA1*COH1)
    RPCADX=AMAX1(-PCAPD1,TPD*(CH2P1-CH2PD))
!
!     print*,'HYDROXYAPATITE'
!
!     CH2PH,CH2P1=equilibrium,current H2PO4 concentration in non-band
!     SYCAH2=solubility product derived from SPALO
!     RPCAHX=H2PO4 dissolution from apatite in non-band
!
!  PCAPH1: precipitated Ca5(PO4)3OH
!  Ca5(PO4)3OH+6H2O<->5Ca(2+) +7OH(-)+3H2PO4(-)
!
    CH2PH=(SYCAH2/(CCA1**5*COH1**7))**0.333_r8

!    print*,'dissocilation rate',PCAPH1,TPD,CH2P1,CH2PH
    RPCAHX=AMAX1(-PCAPH1,TPD*(CH2P1-CH2PH))
!
!     print*,'MONOCALCIUM PHOSPHATE'
!
!     CH2PM,CH2P1=equilibrium,current H2PO4 concentration in non-band
!     SPCAM=solubility product for Ca(H2PO4)2
!     RPCAMX=H2PO4 dissolution from Ca(H2PO4)2 in non-band
!   PCAPM1: precipitated Ca(H2PO4)2
!   Ca(H2PO4)2 <-> Ca(2+)+2H2PO4(-)
    CH2PM=SQRT(SPCAM/CCA1)
!   dissociation rate
    RPCAMX=AMAX1(-PCAPM1*SPPO4,TPD*(CH2P1-CH2PM))
!
!     PHOSPHORUS ANION EXCHANGE IN NON-BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!
!    print*,'XAEC: anaion exchange capacity'
    IF(XAEC.GT.ZEROS)THEN
!     Anion adsorption equilibria
!     H2PO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     RXH2P,RYH2P=H2PO4 exchange with R-OH2,R-OH in non-band
!     TADA: rate constant for adsorption equilibria
!     XOH21: X-OH2(+)
!     XOH11: X-OH
!     X-H2PO4+H2O   <-> X-OH2(+)+H2PO4(-)
!     X-H2PO4+OH(-) <-> X-OH+H2PO4(-)
      SPH2P=SXH2P*DPH2O
      RXH2P=TADA*(XOH21*CH2P1-SPH2P*XH2P1)/(XOH21+SPH2P)*VLWatMicPBK
      RYH2P=TADA*(XOH11*CH2P1-SXH2P*COH1*XH2P1)/(XOH11+SXH2P)*VLWatMicPBK
!
!     HPO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     RXH1P=HPO4 exchange with R-OH in non-band
!     X-HPO4(-)+ OH(-)<-> X-OH+HPO4(2-)
      SPH1P=SXH1P*DPH2O/DPH2P
      RXH1P=TADA*(XOH11*CH1P1-SPH1P*XH1P1)/(XOH11+SPH1P)*VLWatMicPBK
    ELSE
      RXH2P=0._r8
      RYH2P=0._r8
      RXH1P=0._r8
    ENDIF
!
!    print*,'H2PO4(-) <-> H(+)+HPO4(2-)'
!
!     DPH2P=dissociation constant
!     S1=equilibrium concentration in non-band
!     RH2P=H2PO4- <-> H(+) + HPO4(--) dissociation in non-band
!   The problem is equivalent to solving a quadratic equation, with S1 being the \Delta
! of the following equation
!  x^2+(A+B+K)x+AB-KC=0
!  corresponding to
!  (AB) <-> A+B,  (A+x)(B+x)/((AB)-x)=K
    DP=DPH2P
    S0=CH1P1+CHY1+DP
    S1=AMAX1(0._r8,S0**2-4.0_r8*(CH1P1*CHY1-DP*CH2P1))
    RH2P=TSL*(S0-SQRT(S1))
  ELSE
    RPALPX=0._r8
    RPFEPX=0._r8
    RPCADX=0._r8
    RPCAHX=0._r8
    RPCAMX=0._r8
    RXH2P=0._r8
    RYH2P=0._r8
    RXH1P=0._r8
    RH2P=0._r8
  ENDIF

!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN BAND SOIL ZONE
!
  IF(VLWatMicPPB.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
!     CH2PA,CH2PB=equilibrium,current H2PO4 concentration in band
!     SYA0P2=solubility product derived from SPALO
!     RPALBX=H2PO4 dissolution from AlPO4 in band
! Al(3+) and OH(-) are shared between band and non-band soils
!
    CH2PA=SYA0P2/(CAL1*COH1**2)
    RPALBX=AMAX1(-PALPOB,TPD*(CH2PB-CH2PA))
!
!     IRON PHOSPHATE (STRENGITE)
!
!     CH2PF,CH2PB=equilibrium,current H2PO4 concentration in band
!     SYF0P2=solubility product derived from SPALO
!     RPFEBX=H2PO4 dissolution from FePO4 in band
! Fe(3+) is shared between band and non-band soils
    CH2PF=SYF0P2/(CFE1*COH1**2._r8)
    RPFEBX=AMAX1(-PFEPOB,TPD*(CH2PB-CH2PF))
!
!     print*,'DICALCIUM PHOSPHATE'
!
!     CH2PD,CH2PB=equilibrium,current H2PO4 concentration in band
!     SYCAD2=solubility product derived from SPALO
!     RPCDBX=H2PO4 dissolution from CaHPO4 in band
!
    CH2PD=SYCAD2/(CCA1*COH1)
    RPCDBX=AMAX1(-PCAPDB,TPD*(CH2PB-CH2PD))
!
!     HYDROXYAPATITE
!
!     CH2PH,CH2PB=equilibrium,current H2PO4 concentration in band
!     SYCAH2=solubility product derived from SPALO
!     RPCHBX=H2PO4 dissolution from apatite in band
!
    CH2PH=(SYCAH2/(CCA1**5*COH1**7))**0.333_r8
    RPCHBX=AMAX1(-PCAPHB,TPD*(CH2PB-CH2PH))
!
!     print*,'MONOCALCIUM PHOSPHATE'
!
!     CH2PM,CH2PB=equilibrium,current H2PO4 concentration in band
!     SPCAM=solubility product for Ca(H2PO4)2
!     RPCMBX=H2PO4 dissolution from Ca(H2PO4)2 in band
!
    CH2PM=SQRT(SPCAM/CCA1)
    RPCMBX=AMAX1(-PCAPMB*SPPO4,TPD*(CH2PB-CH2PM))
!     IF(I.GT.315)THEN
!     WRITE(*,1117)'RPPOB',I,J,L,RPCMBX,CH2PM,CH2PB,SPCAM,CCA1
!    2,PCAPMB,SPPO4,TPD,PCPMB,BKVLPB
!     ENDIF
!
!     PHOSPHORUS ANION EXCHANGE IN BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!  Anion exchange concentration is also shared among band and non-band soils

    IF(XAEC.GT.ZEROS)THEN
!
!     H2PO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     RXH2B,RYH2B=H2PO4 exchange with R-OH2,R-OH in band
!     R-H2PO4 <-> R-OH2+H2PO4(-)
!   CH2PB: H2PO4(-)
!   X2P1B: R-H2PO4
!   XH21B: R-OH2
      SPH2P=SXH2P*DPH2O
      RXH2B=TADA*(XH21B*CH2PB-SPH2P*X2P1B)/(XH21B+SPH2P)*VLWatMicPBK
      RYH2B=TADA*(XH11B*CH2PB-SXH2P*X2P1B*COH1)/(XH11B+SXH2P)*VLWatMicPBK
!
!     HPO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     RXH1B=HPO4 exchange with R-OH in band
!
      SPH1P=SXH1P*DPH2O/DPH2P
      RXH1B=TADA*(XH11B*CH1PB-SPH1P*X1P1B)/(XH11B+SPH1P)*VLWatMicPBK
    ELSE
      RXH2B=0._r8
      RYH2B=0._r8
      RXH1B=0._r8
    ENDIF
!
!     H2PO4-H+HPO4
!
!     DPH2P=dissociation constant
!     S1=eqilibriunm concentration in band
!     RH2B=H2PO4-H+HPO4 dissociation in band
!
    DP=DPH2P
    S0=CH1PB+CHY1+DP
    S1=AMAX1(0._r8,S0**2-4.0*(CH1PB*CHY1-DP*CH2PB))
    RH2B=TSLX*(S0-SQRT(S1))
  ELSE
    RPALBX=0._r8
    RPFEBX=0._r8
    RPCDBX=0._r8
    RPCHBX=0._r8
    RPCMBX=0._r8
    RXH2B=0._r8
    RYH2B=0._r8
    RXH1B=0._r8
    RH2B=0._r8
  ENDIF
!
!     CATION EXCHANGE FROM GAPON SELECTIVITY COEFFICIENTS
!     FOR CA-NH4, CA-H, CA-AL
! cation-exchange capacity
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
!     CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
!     X*Q=equilibrium exchangeable concentrations
!     XTLQ=total equilibrium exchangeable concentration
!   non-band X-NH4 + H(+) <-> X-H + NH4(+)
    CN41=AMAX1(ZERO,CN41)

!   band X-NH4 <-> X-H + NH4(+)
    CN4B=AMAX1(ZERO,CN4B)

!   X-Al <-> X-H + (1/3)Al(3+)
    CALX=AMAX1(ZERO,CAL1)**0.333

!   X-Fe <-> X-H + (1/3)Fe(3+)
    CFEX=AMAX1(ZERO,CFE1)**0.333

!   X-Ca <-> X-H + (1/2)Ca(2+)
    CCAX=AMAX1(ZERO,CCA1)**0.500

!   X-Mg <-> X-H + (1/2)Mg(2+)
    CMGX=AMAX1(ZERO,CMG1)**0.500

!   X-Na <-> X-H + Na(+)
    CNA1=AMAX1(ZERO,CNA1)

!   X-K  <-> X-H + K(+)
    CKA1=AMAX1(ZERO,CKA1)
!
!   EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
!   CONCENTRATIONS
!
!   The following equation can be found in Robbins et al. (1980)
!   XCAX designates X-Ca
!   Given a general reaction
!   X-H + 1/n A(n+) <-> X-A(n+) + H(+)
!   X-A(n+)/[(A)^(1/n)*(X-H)]=GK
!   X-A(n+)=GK*(A)^(1/n)*(X-H)
    XCAX=CCEC/(1.0+GKC4*CN41/CCAX*VLNH4 &
      +GKC4*CN4B/CCAX*VLNHB &
      +GKCH*CHY1/CCAX+GKCA*CALX/CCAX &
      +GKCA*CFEX/CCAX+GKCM*CMGX/CCAX &
      +GKCN*CNA1/CCAX+GKCK*CKA1/CCAX)
!   X-H+NH4(+) <-> X-NH4 + H(+)
!   XN4Q: equilibrium exchangeable NH4(+) non-band soil
    XN4Q=XCAX*CN41*GKC4

!   XNBq: equilibrium exchangeable NH4(+) band soil
    XNBQ=XCAX*CN4B*GKC4
!   XHYQ: equilibrium exchangeable H(+) in soil
    XHYQ=XCAX*CHY1*GKCH

!   3X-H+Al(3+) <-> X-Al+ 3H(+)
!   XALQ: equilibrium exchangeable Al(3+) in soil
    XALQ=XCAX*CALX*GKCA

!   3X-H+Fe(3+) <-> X-Fe + 3(H+)
!   XFEQ: equilibrium exchangeable Fe(3+) soil
    XFEQ=XCAX*CFEX*GKCA

    XCAQ=XCAX*CCAX
    XMGQ=XCAX*CMGX*GKCM
    XNAQ=XCAX*CNA1*GKCN
    XKAQ=XCAX*CKA1*GKCK
    XTLQ=XN4Q*VLNH4+XNBQ*VLNHB+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CCEC/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XN4Q=FX*XN4Q
    XNBQ=FX*XNBQ
!
!    print*,'NH4 EXCHANGE IN NON-BAND AND BAND SOIL ZONES'
!
!     RXN4,RXNB=NH4 adsorption in non-band,band
!     TADC=adsorption rate constant
!  if RNX4>0, XN41 will be increased
    RXN4=TADC*AMAX1(AMIN1((XN4Q-XN41)*CN41/XN4Q,CN41),-XN41)
    RXNB=TADC*AMAX1(AMIN1((XNBQ-XN4B)*CN4B/XNBQ,CN4B),-XN4B)
  ELSE
    RXN4=0._r8
    RXNB=0._r8
  ENDIF
!
!     NH4-NH3+H IN NON-BAND AND BAND SOIL ZONES
!
!     RNH4,RNHB=NH4-NH3+H dissociation in non-band,band
!     DPN4=NH4 dissociation constant
! NH4(+) <-> NH3 + H(+)
! K=[NH3]*[H(+)]/[NH4(+)]
! H(+): CHY1
! CN31: NH3
! CN41: NH4(+)
! DPN4: K
! RNH4: if RNH4 > 0., then the equilibrium should move to association, and NH4(+) will be increased
  IF(VLWatMicPNH.GT.ZEROS2)THEN
    RNH4=(CHY1*CN31-DPN4*CN41)/(DPN4+CHY1)
  ELSE
    RNH4=0._r8
  ENDIF
  IF(VLWatMicPNB.GT.ZEROS2)THEN
    RNHB=(CHY1*CN3B-DPN4*CN4B)/(DPN4+CHY1)
  ELSE
    RNHB=0._r8
  ENDIF
!
!     TOTAL ION FLUXES FOR ALL REACTIONS ABOVE
!
!     RN4S,RN4B=net NH4 flux in non-band,band
!     RN3S,RN3B=net NH3 flux in non-band,band
!     RAL,RFE,RHY,RCA,RMG,RNA,RKA,ROH=net Al,Fe,H,Ca,Mg,Na,K,OH flux
!     RHP1,RHP2=net HPO4,H2PO4 flux in non-band
!     RXH1,RXH2,RX1P,RX2P=net R-OH,R-OH2,R-HPO4,R-H2PO4 in non-band
!     RHB1,RHB2=net HPO4,H2PO4 flux in band
!     RBH1,RBH2,RB1P,RB2P=net R-OH,R-OH2,R-HPO4,R-H2PO4 in band
!
  RN4S=RNH4-RXN4
  RN4B=RNHB-RXNB
  RN3S=-RNH4
  RN3B=-RNHB
  RHP1=-RH2P-RXH1P
  RHP2=RH2P-RXH2P-RYH2P-RPALPX-RPFEPX-RPCADX-2.0_r8*RPCAMX-3.0_r8*RPCAHX
  RHB1=-RH2B-RXH1B
  RHB2=RH2B-RXH2B-RYH2B-RPALBX-RPFEBX-RPCDBX-2.0_r8*RPCMBX-3.0_r8*RPCHBX
  RXH1=-RYH2P-RXH1P
  RXH2=-RXH2P
  RX1P=RXH1P
  RX2P=RXH2P+RYH2P
  RBH1=-RYH2B-RXH1B
  RBH2=-RXH2B
  RB1P=RXH1B
  RB2P=RXH2B+RYH2B
!
!     CONVERT TOTAL ION FLUXES FROM CHANGES IN CONCENTRATION
!     TO CHANGES IN MASS PER UNIT AREA FOR USE IN 'REDIST'
!
!     TRN4S,TRN4B=total NH4 flux in non-band,band
!     TRN3S,TRN3B=total NH3 flux in non-band,band
!     TRH1P,TRH2P=net HPO4,H2PO4 flux in non-band
!     TRH1B,TRH2B=net HPO4,H2PO4 flux in band
!     TRNX4,TRNXB=total NH4 adsorption in non-band,band
!     TRXH1,TRXH2,TRX1P,TRX2P
!     =total R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in non-band
!     TRBH1,TRBH2,TRB1P,TRB2P
!     =total R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in band
!     TRALPO,TRFEPO,TRCAPD,TRCAPH,TRCAPM
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation in non-band
!     TRALPB,TRFEPB,TRCPDB,TRCPHB,TRCPMB
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation in band
!
  TRN4S=TRN4S+RN4S*VLWatMicPNH
  TRN3S=TRN3S+RN3S*VLWatMicPNH
  TRXN4=TRXN4+RXN4*VLWatMicPNH
  TRN4B=TRN4B+RN4B*VLWatMicPNB
  TRN3B=TRN3B+RN3B*VLWatMicPNB
  TRXNB=TRXNB+RXNB*VLWatMicPNB
  TRH1P=TRH1P+RHP1*VLWatMicPPO
  TRH2P=TRH2P+RHP2*VLWatMicPPO
  TRXH1=TRXH1+RXH1*VLWatMicPPO
  TRXH2=TRXH2+RXH2*VLWatMicPPO
  TRX1P=TRX1P+RX1P*VLWatMicPPO
  TRX2P=TRX2P+RX2P*VLWatMicPPO
  TRALPO=TRALPO+RPALPX*VLWatMicPPO
  TRFEPO=TRFEPO+RPFEPX*VLWatMicPPO
  TRCAPD=TRCAPD+RPCADX*VLWatMicPPO
  TRCAPH=TRCAPH+RPCAHX*VLWatMicPPO
  TRCAPM=TRCAPM+RPCAMX*VLWatMicPPO
  TRH1B=TRH1B+RHB1*VLWatMicPPB
  TRH2B=TRH2B+RHB2*VLWatMicPPB
  TRBH1=TRBH1+RBH1*VLWatMicPPB
  TRBH2=TRBH2+RBH2*VLWatMicPPB
  TRB1P=TRB1P+RB1P*VLWatMicPPB
  TRB2P=TRB2P+RB2P*VLWatMicPPB
  TRALPB=TRALPB+RPALBX*VLWatMicPPB
  TRFEPB=TRFEPB+RPFEBX*VLWatMicPPB
  TRCPDB=TRCPDB+RPCDBX*VLWatMicPPB
  TRCPHB=TRCPHB+RPCHBX*VLWatMicPPB
  TRCPMB=TRCPMB+RPCMBX*VLWatMicPPB

  end subroutine NoSaltChemEquilibria

!--------------------------------------------------------------------------
  subroutine PrintChemVar(chemvar)
  implicit none
  type(chem_var_type), target, intent(in) :: chemvar

  write(*,*)'XCEC'  ,chemvar%XCEC
  write(*,*)'SoilMicPMassLayerX' ,chemvar%SoilMicPMassLayerX
  write(*,*)'SoilMicPMassLayer'  ,chemvar%SoilMicPMassLayer
  write(*,*)'VLWatMicPM' ,chemvar%VLWatMicPM
  write(*,*)'VLWatMicPPO',chemvar%VLWatMicPPO
  write(*,*)'VLWatMicPPB',chemvar%VLWatMicPPB
  write(*,*)'VLWatMicPNH',chemvar%VLWatMicPNH
  write(*,*)'VLWatMicPNB',chemvar%VLWatMicPNB
  write(*,*)'VLNH4' ,chemvar%VLNH4
  write(*,*)'VLNHB' ,chemvar%VLNHB
  write(*,*)'ZMG'   ,chemvar%ZMG
  write(*,*)'ZNA'   ,chemvar%ZNA
  write(*,*)'ZKA'   ,chemvar%ZKA
  write(*,*)'XAEC'  ,chemvar%XAEC
  write(*,*)'CO2S'  ,chemvar%CO2S
  write(*,*)'PH'    ,chemvar%PH
  write(*,*)'CHY1'  ,10.0**(-(chemvar%PH-3.0))
  write(*,*)'CAL'   ,chemvar%CAL
  write(*,*)'CFE'   ,chemvar%CFE
  write(*,*)'CCA'   ,chemvar%CCA
  write(*,*)'CH2P1' ,chemvar%CH2P1
  write(*,*)'PALPO1',chemvar%PALPO1
  write(*,*)'PFEPO1',chemvar%PFEPO1
  write(*,*)'PCAPD1',chemvar%PCAPD1
  write(*,*)'PCAPH1',chemvar%PCAPH1
  write(*,*)'PCAPM1',chemvar%PCAPM1
  write(*,*)'XOH21' ,chemvar%XOH21
  write(*,*)'XH2P1' ,chemvar%XH2P1
  write(*,*)'XOH11' ,chemvar%XOH11
  write(*,*)'XH1P1' ,chemvar%XH1P1
  write(*,*)'CH1PB' ,chemvar%CH1PB
  write(*,*)'CH1P1' ,chemvar%CH1P1
  write(*,*)'CH2PB' ,chemvar%CH2PB
  write(*,*)'PALPOB',chemvar%PALPOB
  write(*,*)'PFEPOB',chemvar%PFEPOB
  write(*,*)'PCAPDB',chemvar%PCAPDB
  write(*,*)'PCAPHB',chemvar%PCAPHB
  write(*,*)'PCAPMB',chemvar%PCAPMB
  write(*,*)'XH21B' ,chemvar%XH21B
  write(*,*)'X2P1B' ,chemvar%X2P1B
  write(*,*)'XH11B' ,chemvar%XH11B
  write(*,*)'X1P1B' ,chemvar%X1P1B
  write(*,*)'CH1PB' ,chemvar%CH1PB
  write(*,*)'CN31'  ,chemvar%CN31
  write(*,*)'CN3B'  ,chemvar%CN3B
  write(*,*)'CN41'  ,chemvar%CN41
  write(*,*)'CN4B'  ,chemvar%CN4B
  write(*,*)'XN41'  ,chemvar%XN41
  write(*,*)'XN4B'  ,chemvar%XN4B
  IF(chemvar%VLWatMicPNH>0. .AND. chemvar%VLWatMicPNB>0. &
    .AND. chemvar%VLWatMicPPO>0. .AND. chemvar%VLWatMicPPB>0.)PAUSE
  end subroutine PrintChemVar

end module ChemEquilibriaMod
