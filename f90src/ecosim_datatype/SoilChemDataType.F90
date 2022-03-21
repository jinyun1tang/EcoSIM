module SoilChemDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none

  save

  real(r8),allocatable :: CZ2GS(:,:,:)
  real(r8),allocatable :: CNH4S(:,:,:)
  real(r8),allocatable :: CNH3S(:,:,:)
  real(r8),allocatable :: CNO3S(:,:,:)
  real(r8),allocatable :: CPO4S(:,:,:)
  real(r8),allocatable :: CNH4B(:,:,:)
  real(r8),allocatable :: CNH3B(:,:,:)
  real(r8),allocatable :: CNO3B(:,:,:)
  real(r8),allocatable :: CPO4B(:,:,:)
  real(r8),allocatable :: CNO2S(:,:,:)
  real(r8),allocatable :: CNH3G(:,:,:)
  real(r8),allocatable :: CZ2GG(:,:,:)
  real(r8),allocatable :: CZ2OG(:,:,:)
  real(r8),allocatable :: CZ2OS(:,:,:)
  real(r8),allocatable :: OXYG(:,:,:)
  real(r8),allocatable :: OXYS(:,:,:)
  real(r8),allocatable :: OXYSH(:,:,:)
  real(r8),allocatable :: CO2G(:,:,:)
  real(r8),allocatable :: CO2S(:,:,:)
  real(r8),allocatable :: CO2SH(:,:,:)
  real(r8),allocatable :: CH4G(:,:,:)
  real(r8),allocatable :: CH4S(:,:,:)
  real(r8),allocatable :: CH4SH(:,:,:)
  real(r8),allocatable :: COXYG(:,:,:)
  real(r8),allocatable :: CCH4G(:,:,:)
  real(r8),allocatable :: COXYS(:,:,:)
  real(r8),allocatable :: CCO2G(:,:,:)
  real(r8),allocatable :: CCO2S(:,:,:)
  real(r8),allocatable :: CCH4S(:,:,:)
  real(r8),allocatable :: CH1P4(:,:,:)
  real(r8),allocatable :: CH1P4B(:,:,:)
  real(r8),allocatable :: CNO2B(:,:,:)
  real(r8),allocatable :: H2GS(:,:,:)
  real(r8),allocatable :: CH2GS(:,:,:)
  real(r8),allocatable :: CH2P4(:,:,:)
  real(r8),allocatable :: CH2P4B(:,:,:)
  real(r8),allocatable :: H2GSH(:,:,:)
  real(r8),allocatable :: H2GG(:,:,:)
  real(r8),allocatable :: CH2GG(:,:,:)
  real(r8),allocatable :: PH(:,:,:)
  real(r8),allocatable :: CEC(:,:,:)
  real(r8),allocatable :: AEC(:,:,:)

  real(r8) :: CNH4(JZ,JY,JX)                    !soil NH4 content, [mg kg-1]
  real(r8) :: CNO3(JZ,JY,JX)                    !soil NO3 content, [mg kg-1]
  real(r8) :: CPO4(JZ,JY,JX)                    !soil PO4 content, [mg kg-1]
  real(r8) :: CAL(JZ,JY,JX)                     !soil Al content, [mg kg-1]
  real(r8) :: CFE(JZ,JY,JX)                     !soil Fe content, [mg kg-1]
  real(r8) :: CCA(JZ,JY,JX)                     !soil Ca content, [mg kg-1]
  real(r8) :: CMG(JZ,JY,JX)                     !soil Mg content, [mg kg-1]
  real(r8) :: CNA(JZ,JY,JX)                     !soil Na content, [mg kg-1]
  real(r8) :: CKA(JZ,JY,JX)                     !soil K content, [mg kg-1]
  real(r8) :: CSO4(JZ,JY,JX)                    !soil SO4 content, [mg kg-1]
  real(r8) :: CCL(JZ,JY,JX)                     !soil Cl content, [mg kg-1]
  real(r8) :: CALOH(JZ,JY,JX)                   !soil AlOH3 content, [mg kg-1]
  real(r8) :: CFEOH(JZ,JY,JX)                   !soil FeOH3 content, [mg kg-1]
  real(r8) :: CCACO(JZ,JY,JX)                   !soil CaCO3 content, [mg kg-1]
  real(r8) :: CCASO(JZ,JY,JX)                   !soil CaSO4 content, [mg kg-1]
  real(r8) :: CALPO(JZ,JY,JX)                   !soil AlPO4 content, [mg kg-1]
  real(r8) :: CFEPO(JZ,JY,JX)                   !soil FePO4 content, [mg kg-1]
  real(r8) :: CCAPD(JZ,JY,JX)                   !soil CaHPO4 content, [mg kg-1]
  real(r8) :: CCAPH(JZ,JY,JX)                   !soil apatite content, [mg kg-1]
  real(r8) :: GKC4(JZ,JY,JX)                    !Ca-NH4 Gapon selectivity coefficient, [-]
  real(r8) :: GKCH(JZ,JY,JX)                    !Ca-H Gapon selectivity coefficient, [-]
  real(r8) :: GKCA(JZ,JY,JX)                    !Ca-Al Gapon selectivity coefficient, [-]
  real(r8) :: GKCM(JZ,JY,JX)                    !Ca-Mg Gapon selectivity coefficient, [-]
  real(r8) :: GKCN(JZ,JY,JX)                    !Ca-Na Gapon selectivity coefficient, [-]
  real(r8) :: GKCK(JZ,JY,JX)                    !Ca-K Gapon selectivity coefficient, [-]

  real(r8) :: ZLSGL(0:JZ,JY,JX)                 !aqueous N2 diffusivity, [m2 h-1]
  real(r8) :: ZHSGL(JZ,JY,JX)                   !aqueous NH4 diffusivity, [m2 h-1]
  real(r8) :: ZNSGL(0:JZ,JY,JX)                 !aqueous NH3 diffusivity, [m2 h-1]
  real(r8) :: ZOSGL(0:JZ,JY,JX)                 !aqueous NO3 diffusivity, [m2 h-1]
  real(r8) :: POSGL(0:JZ,JY,JX)                 !aqueous PO4 diffusivity, [m2 h-1]
  real(r8) :: OCSGL(0:JZ,JY,JX)                 !aqueous DOC diffusivity, [m2 h-1]
  real(r8) :: ONSGL(0:JZ,JY,JX)                 !aqueous DON diffusivity, [m2 h-1]
  real(r8) :: OPSGL(0:JZ,JY,JX)                 !aqueous DOP diffusivity, [m2 h-1]
  real(r8) :: OASGL(0:JZ,JY,JX)                 !aqueous acetate diffusivity, [m2 h-1]
  real(r8) :: Z2SGL(JZ,JY,JX)                   !gaseous N2O diffusivity, [m2 h-1]
  real(r8) :: ZVSGL(0:JZ,JY,JX)                 !aqueous N2O diffusivity, [m2 h-1]
  real(r8) :: WGSGL(JZ,JY,JX)                   !water vapor diffusivity, [m2 h-1]
  real(r8) :: WGSGW(JS,JY,JX)                   !water vapor diffusivity, [m2 h-1]
  real(r8) :: WGSGR(JY,JX)                      !water vapor diffusivity, [m2 h-1]
  real(r8) :: WGSGA(JY,JX)                      !water vapor diffusivity, [m2 h-1]
  real(r8) :: SCO2L(0:JZ,JY,JX)                 !solubility of CO2, [m3 m-3]
  real(r8) :: SOXYL(0:JZ,JY,JX)                 !solubility of O2, [m3 m-3]
  real(r8) :: SCH4L(0:JZ,JY,JX)                 !solubility of CH4, [m3 m-3]
  real(r8) :: SN2OL(0:JZ,JY,JX)                 !solubility of N2O, [m3 m-3]
  real(r8) :: SN2GL(0:JZ,JY,JX)                 !solubility of N2, [m3 m-3]
  real(r8) :: SNH3L(0:JZ,JY,JX)                 !solubility of NH3, [m3 m-3]
  real(r8) :: SH2GL(0:JZ,JY,JX)                 !solubility of H2, [m3 m-3]
  real(r8) :: HGSGL(JZ,JY,JX)                   !gaseous H2 diffusivity, [m2 h-1]
  real(r8) :: HLSGL(0:JZ,JY,JX)                 !aqueous H2 diffusivity, [m2 h-1]

  real(r8) :: H2PO4(0:JZ,JY,JX)                 !PO4 non-band micropore, [g d-2]
  real(r8) :: ZNH4B(0:JZ,JY,JX)                 !NH4 band micropore, [g d-2]
  real(r8) :: ZNH3B(0:JZ,JY,JX)                 !NH3 band micropore, [g d-2]
  real(r8) :: ZNO3B(0:JZ,JY,JX)                 !NO3 band micropore, [g d-2]
  real(r8) :: H2POB(0:JZ,JY,JX)                 !PO4 band micropore, [g d-2]
  real(r8) :: ZNO2S(0:JZ,JY,JX)                 !NO2  non-band micropore, [g d-2]
  real(r8) :: ZNH3G(JZ,JY,JX)                   !gaseous NH3, [g d-2]
  real(r8) :: Z2GG(JZ,JY,JX)                    !gaseous N2, [g d-2]
  real(r8) :: Z2GS(0:JZ,JY,JX)                  !aqueous N2 micropore, [g d-2]
  real(r8) :: Z2OG(JZ,JY,JX)                    !gaseous N2O, [g d-2]
  real(r8) :: Z2OS(0:JZ,JY,JX)                  !aqueous N2O micropore, [g d-2]
  real(r8) :: ZNH4SH(JZ,JY,JX)                  !NH4 non-band macropore, [g d-2]
  real(r8) :: ZNH3SH(JZ,JY,JX)                  !NH3 non-band macropore, [g d-2]
  real(r8) :: ZNO3SH(JZ,JY,JX)                  !NO3 non-band macropore, [g d-2]
  real(r8) :: H2PO4H(JZ,JY,JX)                  !PO4 non-band macropore, [g d-2]
  real(r8) :: ZNH4BH(JZ,JY,JX)                  !NH4 band macropore, [g d-2]
  real(r8) :: ZNH3BH(JZ,JY,JX)                  !NH3 band macropore, [g d-2]
  real(r8) :: ZNO3BH(JZ,JY,JX)                  !NO3 band macropore, [g d-2]
  real(r8) :: H2POBH(JZ,JY,JX)                  !PO4 band macropore, [g d-2]
  real(r8) :: ZNO2SH(JZ,JY,JX)                  !NO2  non-band macropore, [g d-2]
  real(r8) :: Z2GSH(JZ,JY,JX)                   !aqueous N2 macropore, [g d-2]
  real(r8) :: Z2OSH(JZ,JY,JX)                   !aqueous N2O macropore, [g d-2]
  real(r8) :: ZNO2BH(JZ,JY,JX)                  !NO2 band macropore, [g d-2]
  real(r8) :: ZNO2B(0:JZ,JY,JX)                 !NO2  band micropore, [g d-2]
  real(r8) :: ZNH4S(0:JZ,JY,JX)                 !NH4 non-band micropore, [g d-2]
  real(r8) :: ZNH3S(0:JZ,JY,JX)                 !NH3 non-band micropore, [g d-2]
  real(r8) :: ZNO3S(0:JZ,JY,JX)                 !NO3 non-band micropore, [g d-2]
  real(r8) :: H1PO4(0:JZ,JY,JX)                 !soil aqueous HPO4 content micropore non-band, [mol d-2]
  real(r8) :: H1POB(0:JZ,JY,JX)                 !soil aqueous HPO4 content micropore band, [mol d-2]
  real(r8) :: H1PO4H(JZ,JY,JX)                  !soil aqueous HPO4 content non-band macropore, [mol d-2]
  real(r8) :: H1POBH(JZ,JY,JX)                  !soil aqueous HPO4 content band macropore, [mol d-2]
  real(r8) :: ZNFNI(0:JZ,JY,JX)                 !current nitrification inhibition activity
  real(r8) :: ZNFN0(0:JZ,JY,JX)                 !initial nitrification inhibition activity
  real(r8) :: ZNHUI(0:JZ,JY,JX)                 !current inhibition activity
  real(r8) :: ZNHU0(0:JZ,JY,JX)                 !urea hydrolysis inhibition activity

  real(r8) :: XCODFG(0:JZ,JY,JX)                !soil CO2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XCHDFG(0:JZ,JY,JX)                !soil CH4 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XOXDFG(0:JZ,JY,JX)                !soil O2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XNGDFG(0:JZ,JY,JX)                !soil N2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XN2DFG(0:JZ,JY,JX)                !soil N2O dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XN3DFG(0:JZ,JY,JX)                !soil NH3 dissolution (+ve) - volatilization (-ve) non-band, [g d-2 h-1]
  real(r8) :: XNBDFG(0:JZ,JY,JX)                !soil NH3 dissolution (+ve) - volatilization (-ve) band, [g d-2 h-1]
  real(r8) :: XHGDFG(0:JZ,JY,JX)                !soil H2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: RCO2F(0:JZ,JY,JX)                 !net gaseous CO2 flux, [g d-2 h-1]
  real(r8) :: RCH4L(0:JZ,JY,JX)                 !net aqueous CH4 flux, [g d-2 h-1]
  real(r8) :: ROXYF(0:JZ,JY,JX)                 !net gaseous O2 flux, [g d-2 h-1]
  real(r8) :: ROXYL(0:JZ,JY,JX)                 !net aqueous O2 flux, [g d-2 h-1]
  real(r8) :: RCH4F(0:JZ,JY,JX)                 !net gaseous CH4 flux, [g d-2 h-1]
  real(r8) :: ZAL(0:JZ,JY,JX)                   !soil aqueous Al content micropore, [mol d-2]
  real(r8) :: ZFE(0:JZ,JY,JX)                   !soil aqueous Fe content micropore, [mol d-2]
  real(r8) :: ZHY(0:JZ,JY,JX)                   !soil aqueous H content micropore, [mol d-2]
  real(r8) :: ZCA(0:JZ,JY,JX)                   !soil aqueous Ca content micropore, [mol d-2]
  real(r8) :: ZMG(0:JZ,JY,JX)                   !soil aqueous Mg content micropore, [mol d-2]
  real(r8) :: ZNA(0:JZ,JY,JX)                   !soil aqueous Na content micropore, [mol d-2]
  real(r8) :: ZKA(0:JZ,JY,JX)                   !soil aqueous K content micropore, [mol d-2]
  real(r8) :: ZOH(0:JZ,JY,JX)                   !soil aqueous OH content micropore, [mol d-2]
  real(r8) :: ZSO4(0:JZ,JY,JX)                  !soil aqueous SO4 content micropore, [mol d-2]
  real(r8) :: ZCL(0:JZ,JY,JX)                   !soil aqueous Cl content micropore, [mol d-2]
  real(r8) :: ZCO3(0:JZ,JY,JX)                  !soil aqueous CO3 content micropore, [mol d-2]
  real(r8) :: ZHCO3(0:JZ,JY,JX)                 !soil aqueous HCO3 content micropore, [mol d-2]
  real(r8) :: ZALOH1(0:JZ,JY,JX)                !soil aqueous AlOH content micropore, [mol d-2]
  real(r8) :: ZALOH2(0:JZ,JY,JX)                !soil aqueous AlOH2 content micropore, [mol d-2]
  real(r8) :: ZALOH3(0:JZ,JY,JX)                !soil aqueous AlOH3 content micropore, [mol d-2]
  real(r8) :: ZALOH4(0:JZ,JY,JX)                !soil aqueous AlOH4 content micropore, [mol d-2]
  real(r8) :: ZALS(0:JZ,JY,JX)                  !soil aqueous AlSO4 content micropore, [mol d-2]
  real(r8) :: ZFEOH1(0:JZ,JY,JX)                !soil aqueous FeOH content micropore, [mol d-2]
  real(r8) :: ZFEOH2(0:JZ,JY,JX)                !soil aqueous FeOH2 content micropore, [mol d-2]
  real(r8) :: ZFEOH3(0:JZ,JY,JX)                !soil aqueous FeOH3 content micropore, [mol d-2]
  real(r8) :: ZFEOH4(0:JZ,JY,JX)                !soil aqueous FeOH4 content micropore, [mol d-2]
  real(r8) :: ZFES(0:JZ,JY,JX)                  !soil aqueous FeSO4 content micropore, [mol d-2]
  real(r8) :: ZCAO(0:JZ,JY,JX)                  !soil aqueous CaOH2 content micropore, [mol d-2]
  real(r8) :: ZCAC(0:JZ,JY,JX)                  !soil aqueous CACO3 content micropore, [mol d-2]
  real(r8) :: ZCAH(0:JZ,JY,JX)                  !soil aqueous CaHCO3 content micropore, [mol d-2]
  real(r8) :: ZCAS(0:JZ,JY,JX)                  !soil aqueous CaSO4 content micropore, [mol d-2]
  real(r8) :: ZMGO(0:JZ,JY,JX)                  !soil aqueous MgOH content micropore, [mol d-2]
  real(r8) :: ZMGC(0:JZ,JY,JX)                  !soil aqueous MgCO3 content micropore, [mol d-2]
  real(r8) :: ZMGH(0:JZ,JY,JX)                  !soil aqueous MgHCO3 content micropore, [mol d-2]
  real(r8) :: ZMGS(0:JZ,JY,JX)                  !soil aqueous MgSO4 content micropore, [mol d-2]
  real(r8) :: ZNAC(0:JZ,JY,JX)                  !soil aqueous NaCO3 content micropore, [mol d-2]
  real(r8) :: ZNAS(0:JZ,JY,JX)                  !soil aqueous NaSO4 content micropore, [mol d-2]
  real(r8) :: ZKAS(0:JZ,JY,JX)                  !soil aqueous KSO4 content micropore, [mol d-2]
  real(r8) :: H0PO4(0:JZ,JY,JX)                 !soil aqueous PO4 content micropore non-band, [mol d-2]
  real(r8) :: H3PO4(0:JZ,JY,JX)                 !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  real(r8) :: ZFE1P(0:JZ,JY,JX)                 !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  real(r8) :: ZFE2P(0:JZ,JY,JX)                 !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  real(r8) :: ZCA0P(0:JZ,JY,JX)                 !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  real(r8) :: ZCA1P(0:JZ,JY,JX)                 !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  real(r8) :: ZCA2P(0:JZ,JY,JX)                 !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  real(r8) :: ZMG1P(0:JZ,JY,JX)                 !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  real(r8) :: H0POB(JZ,JY,JX)                   !soil aqueous PO4 content micropore band, [mol d-2]
  real(r8) :: H3POB(JZ,JY,JX)                   !soil aqueous H3PO4 content micropore band, [mol d-2]
  real(r8) :: ZFE1PB(JZ,JY,JX)                  !soil aqueous FeHPO4 content micropore band, [mol d-2]
  real(r8) :: ZFE2PB(JZ,JY,JX)                  !soil aqueous FeH2PO4 content micropore band, [mol d-2]
  real(r8) :: ZCA0PB(JZ,JY,JX)                  !soil aqueous CaPO4 content micropore band, [mol d-2]
  real(r8) :: ZCA1PB(JZ,JY,JX)                  !soil aqueous CaHPO4 content micropore band, [mol d-2]
  real(r8) :: ZCA2PB(JZ,JY,JX)                  !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  real(r8) :: ZMG1PB(JZ,JY,JX)                  !soil aqueous MgHPO4 content micropore band, [mol d-2]
  real(r8) :: XN4(0:JZ,JY,JX)                   !exchangeable NH4 non-band, [mol d-2]
  real(r8) :: XNB(0:JZ,JY,JX)                   !exchangeable NH4 band, [mol d-2]
  real(r8) :: XHY(JZ,JY,JX)                     !exchangeable H , [mol d-2]
  real(r8) :: XAL(JZ,JY,JX)                     !exchangeable Al, [mol d-2]
  real(r8) :: XCA(JZ,JY,JX)                     !exchangeable Ca, [mol d-2]
  real(r8) :: XMG(JZ,JY,JX)                     !exchangeable Mg , [mol d-2]
  real(r8) :: XNA(JZ,JY,JX)                     !exchangeable Na, [mol d-2]
  real(r8) :: XKA(JZ,JY,JX)                     !exchangeable K, [mol d-2]
  real(r8) :: XHC(JZ,JY,JX)                     !exchangeable COOH , [mol d-2]
  real(r8) :: XALO2(JZ,JY,JX)                   !exchangeable AlOH2 , [mol d-2]
  real(r8) :: XOH0(0:JZ,JY,JX)                  !exchangeable OH- non-band, [mol d-2]
  real(r8) :: XOH1(0:JZ,JY,JX)                  !exchangeable OH  non-band, [mol d-2]
  real(r8) :: XOH2(0:JZ,JY,JX)                  !exchangeable OH2  non-band, [mol d-2]
  real(r8) :: XH1P(0:JZ,JY,JX)                  !exchangeable HPO4  non-band, [mol d-2]
  real(r8) :: XH2P(0:JZ,JY,JX)                  !exchangeable H2PO4  non-band, [mol d-2]
  real(r8) :: XOH0B(0:JZ,JY,JX)                 !exchangeable OH- band, [mol d-2]
  real(r8) :: XFE(JZ,JY,JX)
  real(r8) :: XFEO2(JZ,JY,JX)

  real(r8) :: XOH1B(0:JZ,JY,JX)                 !exchangeable OH  band, [mol d-2]
  real(r8) :: XOH2B(0:JZ,JY,JX)                 !exchangeable OH2  band, [mol d-2]
  real(r8) :: XH1PB(0:JZ,JY,JX)                 !exchangeable HPO4  band, [mol d-2]
  real(r8) :: XH2PB(0:JZ,JY,JX)                 !exchangeable H2PO4  band, [mol d-2]
  real(r8) :: PALOH(JZ,JY,JX)                   !precipitated AlOH3, [mol d-2]
  real(r8) :: PFEOH(JZ,JY,JX)                   !precipitated FeOH3, [mol d-2]
  real(r8) :: PCACO(JZ,JY,JX)                   !precipitated CaCO3, [mol d-2]
  real(r8) :: PCASO(JZ,JY,JX)                   !precipitated CaSO4, [mol d-2]
  real(r8) :: PALPO(0:JZ,JY,JX)                 !precipitated AlPO4 non-band, [mol d-2]
  real(r8) :: PFEPO(0:JZ,JY,JX)                 !precipitated FePO4 non-band, [mol d-2]
  real(r8) :: PCAPD(0:JZ,JY,JX)                 !precipitated CaHPO4 non-band, [mol d-2]
  real(r8) :: PCAPH(0:JZ,JY,JX)                 !precipitated hydroxyapatite non-band, [mol d-2]
  real(r8) :: PCAPM(0:JZ,JY,JX)                 !precipitated CaH2PO4 non-band, [mol d-2]
  real(r8) :: PALPB(0:JZ,JY,JX)                 !precipitated AlPO4 band, [mol d-2]
  real(r8) :: PFEPB(0:JZ,JY,JX)                 !precipitated FePO4 band, [mol d-2]
  real(r8) :: PCPDB(0:JZ,JY,JX)                 !precipitated CaHPO4 band, [mol d-2]
  real(r8) :: PCPMB(0:JZ,JY,JX)                 !precipitated CaH2PO4 , [mol d-2]
  real(r8) :: ECND(JZ,JY,JX)                    !electrical conductivity , [dS m-1]
  real(r8) :: CSTR(JZ,JY,JX)                    !solution ion strength, [mol m-3]
  real(r8) :: CION(JZ,JY,JX)                    !solution ion concentratiom, [mol m-3]
  real(r8) :: XCEC(JZ,JY,JX)                    !cation exchange capacity, [mol d-2]
  real(r8) :: XAEC(JZ,JY,JX)                    !anion exchange capacity, [mol d-2]
  real(r8) :: ALSGL(JZ,JY,JX)                   !aqueous Al diffusivity, [m2 h-1]
  real(r8) :: FESGL(JZ,JY,JX)                   !aqueous Fe diffusivity, [m2 h-1]
  real(r8) :: HYSGL(JZ,JY,JX)                   !aqueous H diffusivity, [m2 h-1]
  real(r8) :: CASGL(JZ,JY,JX)                   !aqueous Ca diffusivity, [m2 h-1]
  real(r8) :: GMSGL(JZ,JY,JX)                   !aqueous Mg diffusivity, [m2 h-1]
  real(r8) :: ANSGL(JZ,JY,JX)                   !aqueous Na diffusivity, [m2 h-1]
  real(r8) :: AKSGL(JZ,JY,JX)                   !aqueous K diffusivity, [m2 h-1]
  real(r8) :: OHSGL(JZ,JY,JX)                   !aqueous OH diffusivity, [m2 h-1]
  real(r8) :: C3SGL(JZ,JY,JX)                   !aqueous CO3 diffusivity, [m2 h-1]
  real(r8) :: HCSGL(JZ,JY,JX)                   !aqueous HCO3 diffusivity, [m2 h-1]
  real(r8) :: SOSGL(JZ,JY,JX)                   !aqueous SO4 diffusivity, [m2 h-1]
  real(r8) :: CLSXL(JZ,JY,JX)                   !aqueous Cl diffusivity, [m2 h-1]
  real(r8) :: ZALH(JZ,JY,JX)                    !aqueous Al content macropore, [mol d-2]
  real(r8) :: ZFEH(JZ,JY,JX)                    !aqueous Fe content macropore, [mol d-2]
  real(r8) :: ZHYH(JZ,JY,JX)                    !aqueous H content macropore, [mol d-2]
  real(r8) :: ZCCH(JZ,JY,JX)                    !aqueous CO3 content macropore, [mol d-2]
  real(r8) :: ZMAH(JZ,JY,JX)                    !aqueous Mg content macropore, [mol d-2]
  real(r8) :: ZNAH(JZ,JY,JX)                    !aqueous Na content macropore, [mol d-2]
  real(r8) :: ZKAH(JZ,JY,JX)                    !aqueous K content macropore, [mol d-2]
  real(r8) :: ZOHH(JZ,JY,JX)                    !aqueous OH content macropore, [mol d-2]
  real(r8) :: ZSO4H(JZ,JY,JX)                   !aqueous SO4 content macropore, [mol d-2]
  real(r8) :: ZCLH(JZ,JY,JX)                    !aqueous Cl content macropore, [mol d-2]
  real(r8) :: ZCO3H(JZ,JY,JX)                   !aqueous CO3 content macropore, [mol d-2]
  real(r8) :: ZHCO3H(JZ,JY,JX)                  !aqueous HCO3 content macropore, [mol d-2]
  real(r8) :: ZALO1H(JZ,JY,JX)                  !aqueous AlOH content macropore, [mol d-2]
  real(r8) :: ZALO2H(JZ,JY,JX)                  !aqueous AlOH2 content macropore, [mol d-2]
  real(r8) :: ZALO3H(JZ,JY,JX)                  !aqueous AlOH3 content macropore, [mol d-2]
  real(r8) :: ZALO4H(JZ,JY,JX)                  !aqueous AlOH4 content macropore, [mol d-2]
  real(r8) :: ZALSH(JZ,JY,JX)                   !aqueous AlSO4 content macropore, [mol d-2]
  real(r8) :: ZFEO1H(JZ,JY,JX)                  !aqueous FeOH content macropore, [mol d-2]
  real(r8) :: ZFEO2H(JZ,JY,JX)                  !aqueous FeOH2 content macropore, [mol d-2]
  real(r8) :: ZFEO3H(JZ,JY,JX)                  !aqueous FeOH3 content macropore, [mol d-2]
  real(r8) :: ZFEO4H(JZ,JY,JX)                  !aqueous FeOH4 content macropore, [mol d-2]
  real(r8) :: ZFESH(JZ,JY,JX)                   !aqueous FeSO4 content macropore, [mol d-2]
  real(r8) :: ZCAOH(JZ,JY,JX)                   !aqueous CaOH content macropore, [mol d-2]
  real(r8) :: ZCACH(JZ,JY,JX)                   !aqueous CaCO3 content macropore, [mol d-2]
  real(r8) :: ZCAHH(JZ,JY,JX)                   !aqueous CaHCO3 content macropore, [mol d-2]
  real(r8) :: ZCASH(JZ,JY,JX)                   !aqueous CaSO4 content macropore, [mol d-2]
  real(r8) :: PCPHB(0:JZ,JY,JX)                 !precipitated hydroxyapatite band, [mol d-2]
  real(r8) :: ZMGOH(JZ,JY,JX)                   !aqueous MgOH content macropore, [mol d-2]
  real(r8) :: ZMGCH(JZ,JY,JX)                   !aqueous MgCO3 content macropore, [mol d-2]
  real(r8) :: ZMGHH(JZ,JY,JX)                   !aqueous MgHCO3 content macropore, [mol d-2]
  real(r8) :: ZMGSH(JZ,JY,JX)                   !aqueous MgSO4 content macropore, [mol d-2]
  real(r8) :: ZNACH(JZ,JY,JX)                   !aqueous NaCO3 content macropore, [mol d-2]
  real(r8) :: ZNASH(JZ,JY,JX)                   !aqueous NaSO4 content macropore, [mol d-2]
  real(r8) :: ZKASH(JZ,JY,JX)                   !aqueous KSO4 content macropore, [mol d-2]
  real(r8) :: H0PO4H(JZ,JY,JX)                  !soil aqueous PO4 content non-band macropore , [mol d-2]
  real(r8) :: H3PO4H(JZ,JY,JX)                  !soil aqueous H3PO4 content non-band macropore, [mol d-2]
  real(r8) :: ZFE1PH(JZ,JY,JX)                  !soil aqueous FeHPO4 content non-band macropore, [mol d-2]
  real(r8) :: ZFE2PH(JZ,JY,JX)                  !soil aqueous FeH2PO4 content non-band macropore, [mol d-2]
  real(r8) :: ZCA0PH(JZ,JY,JX)                  !soil aqueous CaPO4 content non-band macropore, [mol d-2]
  real(r8) :: ZCA1PH(JZ,JY,JX)                  !soil aqueous CaHPO4 content non-band macropore, [mol d-2]
  real(r8) :: ZCA2PH(JZ,JY,JX)                  !soil aqueous CaH2PO4 content non-band macropore, [mol d-2]
  real(r8) :: ZMG1PH(JZ,JY,JX)                  !soil aqueous MgHPO4 content non-band macropore, [mol d-2]
  real(r8) :: H0POBH(JZ,JY,JX)                  !soil aqueous PO4 content band macropore , [mol d-2]
  real(r8) :: H3POBH(JZ,JY,JX)                  !soil aqueous H3PO4 content band macropore, [mol d-2]
  real(r8) :: ZFE1BH(JZ,JY,JX)                  !soil aqueous FeHPO4 content band macropore, [mol d-2]
  real(r8) :: ZFE2BH(JZ,JY,JX)                  !soil aqueous FeH2PO4 content band macropore, [mol d-2]
  real(r8) :: ZCA0BH(JZ,JY,JX)                  !soil aqueous CaPO4 content band macropore, [mol d-2]
  real(r8) :: ZCA1BH(JZ,JY,JX)                  !soil aqueous CaHPO4 content band macropore, [mol d-2]
  real(r8) :: ZCA2BH(JZ,JY,JX)                  !soil aqueous CaH2PO4 content band macropore, [mol d-2]
  real(r8) :: ZMG1BH(JZ,JY,JX)                  !soil aqueous MgHPO4 content band macropore, [mol d-2]

  real(r8) :: XQRAL(2,2,JV,JH)                  !total Al in runoff, [mol d-2 h-1]
  real(r8) :: XQRFE(2,2,JV,JH)                  !total Fe in runoff, [mol d-2 h-1]
  real(r8) :: XQRHY(2,2,JV,JH)                  !total H in runoff, [mol d-2 h-1]
  real(r8) :: XQRCA(2,2,JV,JH)                  !total Ca in runoff, [mol d-2 h-1]
  real(r8) :: XQRMG(2,2,JV,JH)                  !total Mg in runoff, [mol d-2 h-1]
  real(r8) :: XQRNA(2,2,JV,JH)                  !total Na in runoff, [mol d-2 h-1]
  real(r8) :: XQRKA(2,2,JV,JH)                  !total K in runoff, [mol d-2 h-1]
  real(r8) :: XQROH(2,2,JV,JH)                  !total OH in runoff, [mol d-2 h-1]
  real(r8) :: XQRSO(2,2,JV,JH)                  !total SO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRCL(2,2,JV,JH)                  !total Cl in runoff, [mol d-2 h-1]
  real(r8) :: XQRC3(2,2,JV,JH)                  !total CO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRHC(2,2,JV,JH)                  !total HCO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRAL1(2,2,JV,JH)                 !total AlOH in runoff, [mol d-2 h-1]
  real(r8) :: XQRAL2(2,2,JV,JH)                 !total AlOH2 in runoff, [mol d-2 h-1]
  real(r8) :: XQRAL3(2,2,JV,JH)                 !total AlOH3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRAL4(2,2,JV,JH)                 !total AlOH4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRALS(2,2,JV,JH)                 !total AlSO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRFE1(2,2,JV,JH)                 !total FeOH in runoff, [mol d-2 h-1]
  real(r8) :: XQRFE2(2,2,JV,JH)                 !total FeOH2 in runoff, [mol d-2 h-1]
  real(r8) :: XQRFE3(2,2,JV,JH)                 !total FeOH3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRFE4(2,2,JV,JH)                 !total FeOH4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRFES(2,2,JV,JH)                 !total FeSO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRCAO(2,2,JV,JH)                 !total CaOH in runoff, [mol d-2 h-1]
  real(r8) :: XQRCAC(2,2,JV,JH)                 !total CaCO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRCAH(2,2,JV,JH)                 !total CaHCO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRCAS(2,2,JV,JH)                 !total CaSO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRMGO(2,2,JV,JH)                 !total MgOH in runoff, [mol d-2 h-1]
  real(r8) :: XQRMGC(2,2,JV,JH)                 !total MgCO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRMGH(2,2,JV,JH)                 !total MgHCO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRMGS(2,2,JV,JH)                 !total MgSO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRNAC(2,2,JV,JH)                 !total NaCO3 in runoff, [mol d-2 h-1]
  real(r8) :: XQRNAS(2,2,JV,JH)                 !total NaSO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRKAS(2,2,JV,JH)                 !total KSO4 in runoff, [mol d-2 h-1]
  real(r8) :: XQRH0P(2,2,JV,JH)                 !total PO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRH3P(2,2,JV,JH)                 !total H3PO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRF1P(2,2,JV,JH)                 !total FeHPO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRF2P(2,2,JV,JH)                 !total FeH2PO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRC0P(2,2,JV,JH)                 !total CaPO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRC1P(2,2,JV,JH)                 !total CaHPO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRC2P(2,2,JV,JH)                 !total CaH2PO4 in runoff non-band, [mol d-2 h-1]
  real(r8) :: XQRM1P(2,2,JV,JH)                 !total MgHPO4 in runoff non-band, [mol d-2 h-1]

  private :: InitAllocate

  contains


  subroutine InitSoilChemData

  implicit none

  call InitAllocate

  end subroutine InitSoilChemData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  allocate(CZ2GS(0:JZ,JY,JX));CZ2GS(0:JZ,JY,JX)=0._r8
  allocate(CNH4S(0:JZ,JY,JX));CNH4S(0:JZ,JY,JX)=0._r8
  allocate(CNH3S(0:JZ,JY,JX));CNH3S(0:JZ,JY,JX)=0._r8
  allocate(CNO3S(0:JZ,JY,JX));CNO3S(0:JZ,JY,JX)=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(CNH4B(0:JZ,JY,JX));CNH4B(0:JZ,JY,JX)=0._r8
  allocate(CNH3B(0:JZ,JY,JX));CNH3B(0:JZ,JY,JX)=0._r8
  allocate(CNO3B(0:JZ,JY,JX));CNO3B(0:JZ,JY,JX)=0._r8
  allocate(CPO4B(0:JZ,JY,JX));CPO4B(0:JZ,JY,JX)=0._r8
  allocate(CNO2S(0:JZ,JY,JX));CNO2S(0:JZ,JY,JX)=0._r8
  allocate(CNH3G(0:JZ,JY,JX));CNH3G(0:JZ,JY,JX)=0._r8
  allocate(CZ2GG(0:JZ,JY,JX));CZ2GG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OG(0:JZ,JY,JX));CZ2OG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OS(0:JZ,JY,JX));CZ2OS(0:JZ,JY,JX)=0._r8
  allocate(OXYG(JZ,JY,JX));OXYG(JZ,JY,JX)=0._r8
  allocate(OXYS(0:JZ,JY,JX));OXYS(0:JZ,JY,JX)=0._r8
  allocate(OXYSH(JZ,JY,JX));OXYSH(JZ,JY,JX)=0._r8
  allocate(CO2G(JZ,JY,JX));CO2G(JZ,JY,JX)=0._r8
  allocate(CO2S(0:JZ,JY,JX));CO2S(0:JZ,JY,JX)=0._r8
  allocate(CO2SH(JZ,JY,JX));CO2SH(JZ,JY,JX)=0._r8
  allocate(CH4G(JZ,JY,JX));CH4G(JZ,JY,JX)=0._r8
  allocate(CH4S(0:JZ,JY,JX));CH4S(0:JZ,JY,JX)=0._r8
  allocate(CH4SH(JZ,JY,JX));CH4SH(JZ,JY,JX)=0._r8
  allocate(COXYG(0:JZ,JY,JX));COXYG(0:JZ,JY,JX)=0._r8
  allocate(CCH4G(0:JZ,JY,JX));CCH4G(0:JZ,JY,JX)=0._r8
  allocate(COXYS(0:JZ,JY,JX));COXYS(0:JZ,JY,JX)=0._r8
  allocate(CCO2G(0:JZ,JY,JX));CCO2G(0:JZ,JY,JX)=0._r8
  allocate(CCO2S(0:JZ,JY,JX));CCO2S(0:JZ,JY,JX)=0._r8
  allocate(CCH4S(0:JZ,JY,JX));CCH4S(0:JZ,JY,JX)=0._r8
  allocate(CH1P4(0:JZ,JY,JX));CH1P4(0:JZ,JY,JX)=0._r8
  allocate(CH1P4B(0:JZ,JY,JX));CH1P4B(0:JZ,JY,JX)=0._r8
  allocate(CNO2B(0:JZ,JY,JX));CNO2B(0:JZ,JY,JX)=0._r8
  allocate(H2GS(0:JZ,JY,JX));H2GS(0:JZ,JY,JX)=0._r8
  allocate(CH2GS(0:JZ,JY,JX));CH2GS(0:JZ,JY,JX)=0._r8
  allocate(CH2P4(0:JZ,JY,JX));CH2P4(0:JZ,JY,JX)=0._r8
  allocate(CH2P4B(0:JZ,JY,JX));CH2P4B(0:JZ,JY,JX)=0._r8
  allocate(H2GSH(JZ,JY,JX));H2GSH(JZ,JY,JX)=0._r8
  allocate(H2GG(JZ,JY,JX));H2GG(JZ,JY,JX)=0._r8
  allocate(CH2GG(0:JZ,JY,JX));CH2GG(0:JZ,JY,JX)=0._r8
  allocate(PH(0:JZ,JY,JX));PH(0:JZ,JY,JX)=0._r8
  allocate(CEC(JZ,JY,JX));CEC(JZ,JY,JX)=0._r8
  allocate(AEC(JZ,JY,JX));AEC(JZ,JY,JX)=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilChemData

  implicit none

  if(allocated(CZ2GS))deallocate(CZ2GS)
  if(allocated(CNH4S))deallocate(CNH4S)
  if(allocated(CNH3S))deallocate(CNH3S)
  if(allocated(CNO3S))deallocate(CNO3S)
  if(allocated(CPO4S))deallocate(CPO4S)
  if(allocated(CNH4B))deallocate(CNH4B)
  if(allocated(CNH3B))deallocate(CNH3B)
  if(allocated(CNO3B))deallocate(CNO3B)
  if(allocated(CPO4B))deallocate(CPO4B)
  if(allocated(CNO2S))deallocate(CNO2S)
  if(allocated(CNH3G))deallocate(CNH3G)
  if(allocated(CZ2GG))deallocate(CZ2GG)
  if(allocated(CZ2OG))deallocate(CZ2OG)
  if(allocated(CZ2OS))deallocate(CZ2OS)
  if(allocated(OXYG))deallocate(OXYG)
  if(allocated(OXYS))deallocate(OXYS)
  if(allocated(OXYSH))deallocate(OXYSH)
  if(allocated(CO2G))deallocate(CO2G)
  if(allocated(CO2S))deallocate(CO2S)
  if(allocated(CO2SH))deallocate(CO2SH)
  if(allocated(CH4G))deallocate(CH4G)
  if(allocated(CH4S))deallocate(CH4S)
  if(allocated(CH4SH))deallocate(CH4SH)
  if(allocated(COXYG))deallocate(COXYG)
  if(allocated(CCH4G))deallocate(CCH4G)
  if(allocated(COXYS))deallocate(COXYS)
  if(allocated(CCO2G))deallocate(CCO2G)
  if(allocated(CCO2S))deallocate(CCO2S)
  if(allocated(CCH4S))deallocate(CCH4S)
  if(allocated(CH1P4))deallocate(CH1P4)
  if(allocated(CH1P4B))deallocate(CH1P4B)
  if(allocated(CNO2B))deallocate(CNO2B)
  if(allocated(H2GS))deallocate(H2GS)
  if(allocated(CH2GS))deallocate(CH2GS)
  if(allocated(CH2P4))deallocate(CH2P4)
  if(allocated(CH2P4B))deallocate(CH2P4B)
  if(allocated(H2GSH))deallocate(H2GSH)
  if(allocated(H2GG))deallocate(H2GG)
  if(allocated(CH2GG))deallocate(CH2GG)
  if(allocated(PH))deallocate(PH)
  if(allocated(CEC))deallocate(CEC)
  if(allocated(AEC))deallocate(AEC)
  end subroutine DestructSoilChemData
end module SoilChemDataType
