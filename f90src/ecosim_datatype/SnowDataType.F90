module SnowDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  REAL(R8) :: ALBS(JY,JX)                       !snowpack albedo
  real(r8) :: DENS0(JY,JX)                      !snowpack density, [Mg m-3]
  real(r8) :: VHCPWM(60,JS,JY,JX)               !volumetric heat capacity of snowpack
  real(r8) :: FLQWM(60,JS,JY,JX)                !snowpack water flux
  real(r8) :: TCW(JS,JY,JX)                     !snow temperature, [oC]
  real(r8) :: TKW(JS,JY,JX)                     !snow temperature, [K]
  real(r8) :: VHCPW(JS,JY,JX)                   !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: VOLSSL(JS,JY,JX)                  !snow water equivalent volume in snowpack layer
  real(r8) :: VOLWSL(JS,JY,JX)                  !snow water volume in snowpack layer
  real(r8) :: VOLISL(JS,JY,JX)                  !snow ice volume in snowpack layer
  real(r8) :: VOLSL(JS,JY,JX)                   !snow volume in snowpack layer
  real(r8) :: DENSS(JS,JY,JX)                   !snowpack density, [Mg m-3]
  real(r8) :: DLYRS(JS,JY,JX)                   !snowpack layer depth
  real(r8) :: XFLWW(JS,JY,JX)                   !hourly snow water transfer
  real(r8) :: XFLWS(JS,JY,JX)                   !hourly snow transfer
  real(r8) :: XFLWI(JS,JY,JX)                   !hourly snow ice transfer
  real(r8) :: XHFLWW(JS,JY,JX)                  !hourly convective heat flux from water transfer
  real(r8) :: XWFLXS(JS,JY,JX)                  !hourly convective heat flux from snow transfer
  real(r8) :: XWFLXI(JS,JY,JX)                  !hourly convective heat flux from ice transfer
  real(r8) :: CDPTHS(0:JS,JY,JX)                !cumulative depth to bottom of snowpack layer
  real(r8) :: VOLSI(JS,JY,JX)                   !Initial snowpack volume, [m3 d-2]
  real(r8) :: DPTHS(JY,JX)                      !snowpack depth, [m]
  real(r8) :: VOLSS(JY,JX)                      !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: VOLWS(JY,JX)                      !water volume in snowpack, [m3 d-2]
  real(r8) :: VOLIS(JY,JX)                      !ice volume in snowpack, [m3 d-2]
  real(r8) :: VOLS(JY,JX)                       !snowpack volume, [m3 d-2]
  real(r8) :: VHCPWX(JY,JX)                     !snowpack heat capacity from previous time step, [MJ d-2 K-1]
  real(r8) :: FLSW(JS,JY,JX)                    !water from snowpack to soil micropores
  real(r8) :: FLSWH(JS,JY,JX)                   !water from snowpack to soil macropores
  real(r8) :: HFLSW(JS,JY,JX)                   !convective heat from snowpack to soil
  real(r8) :: FLSWR(JS,JY,JX)                   !water flux from snowpack to litter
  real(r8) :: HFLSWR(JS,JY,JX)                  !convective heat flux from snowpack to litter
  real(r8) :: QS(2,JV,JH)                       !snowpack runoff snow, [m3 d-2 h-1]
  real(r8) :: QW(2,JV,JH)                       !snowpack runoff water, [m3 d-2 h-1]
  real(r8) :: QI(2,JV,JH)                       !snowpack runoff ice, [m3 d-2 h-1]
  real(r8) :: HQS(2,JV,JH)                      !snowpack runoff heat, [MJ d-2 h-1]
  real(r8) :: QSM(60,2,JV,JH)                   !runoff snow flux, [m3 d-2 t-1]

  real(r8) :: XCOQSS(2,JV,JH)                   !snowpack runoff CO2 flux, [g d-2 h-1]
  real(r8) :: XCHQSS(2,JV,JH)                   !snowpack runoff CH4 flux, [g d-2 h-1]
  real(r8) :: XOXQSS(2,JV,JH)                   !snowpack runoff O2 flux, [g d-2 h-1]
  real(r8) :: XNGQSS(2,JV,JH)                   !snowpack runoff N2 flux, [g d-2 h-1]
  real(r8) :: XN2QSS(2,JV,JH)                   !snowpack runoff N2O flux, [g d-2 h-1]
  real(r8) :: XN4QSS(2,JV,JH)                   !snowpack runoff NH4 flux, [g d-2 h-1]
  real(r8) :: XN3QSS(2,JV,JH)                   !snowpack runoff NH3 flux, [g d-2 h-1]
  real(r8) :: XNOQSS(2,JV,JH)                   !snowpack runoff NO3 flux, [g d-2 h-1]
  real(r8) :: XP4QSS(2,JV,JH)                   !snowpack runoff PO4 flux, [g d-2 h-1]
  real(r8) :: XP1QSS(2,JV,JH)                   !snowpack runoff HPO4 flux, [g d-2 h-1]

  real(r8) :: CO2W(JS,JY,JX)                    !snowpack CO2, [mol d-2]
  real(r8) :: CH4W(JS,JY,JX)                    !snowpack CH4, [mol d-2]
  real(r8) :: OXYW(JS,JY,JX)                    !snowpack O2, [mol d-2]
  real(r8) :: ZN2W(JS,JY,JX)                    !snowpack N2O, [mol d-2]
  real(r8) :: ZNGW(JS,JY,JX)                    !snowpack N2, [mol d-2]
  real(r8) :: ZN4W(JS,JY,JX)                    !snowpack NH4, [mol d-2]
  real(r8) :: ZN3W(JS,JY,JX)                    !snowpack NH3, [mol d-2]
  real(r8) :: ZNOW(JS,JY,JX)                    !snowpack NO3, [mol d-2]
  real(r8) :: Z1PW(JS,JY,JX)                    !snowpack HPO4,[mol d-2]
  real(r8) :: ZHPW(JS,JY,JX)                    !snowpack H2PO4, [mol d-2]
  real(r8) :: ZALW(JS,JY,JX)                    !snowpack Al, [mol d-2]
  real(r8) :: ZFEW(JS,JY,JX)                    !snowpack Fe, [mol d-2]
  real(r8) :: ZHYW(JS,JY,JX)                    !snowpack H, [mol d-2]
  real(r8) :: ZCAW(JS,JY,JX)                    !snowpack Ca, [mol d-2]
  real(r8) :: ZMGW(JS,JY,JX)                    !snowpack Mg, [mol d-2]
  real(r8) :: ZNAW(JS,JY,JX)                    !snowpack Na, [mol d-2]
  real(r8) :: ZKAW(JS,JY,JX)                    !snowpack K, [mol d-2]
  real(r8) :: ZOHW(JS,JY,JX)                    !snowpack OH, [mol d-2]
  real(r8) :: ZSO4W(JS,JY,JX)                   !snowpack SO4, [mol d-2]
  real(r8) :: ZCLW(JS,JY,JX)                    !snowpack Cl, [mol d-2]
  real(r8) :: ZCO3W(JS,JY,JX)                   !snowpack CO3, [mol d-2]
  real(r8) :: ZHCO3W(JS,JY,JX)                  !snowpack HCO3, [mol d-2]
  real(r8) :: ZALH1W(JS,JY,JX)                  !snowpack AlOH, [mol d-2]
  real(r8) :: ZALH2W(JS,JY,JX)                  !snowpack AlOH2, [mol d-2]
  real(r8) :: ZALH3W(JS,JY,JX)                  !snowpack AlOH3, [mol d-2]
  real(r8) :: ZALH4W(JS,JY,JX)                  !snowpack AlOH4, [mol d-2]
  real(r8) :: ZALSW(JS,JY,JX)                   !snowpack AlSO4, [mol d-2]
  real(r8) :: ZFEH1W(JS,JY,JX)                  !snowpack FeOH, [mol d-2]
  real(r8) :: ZFEH2W(JS,JY,JX)                  !snowpack FeOH2, [mol d-2]
  real(r8) :: ZFEH3W(JS,JY,JX)                  !snowpack FeOH3, [mol d-2]
  real(r8) :: ZFEH4W(JS,JY,JX)                  !snowpack F3OH4, [mol d-2]
  real(r8) :: ZFESW(JS,JY,JX)                   !snowpack FeSO4, [mol d-2]
  real(r8) :: ZCAOW(JS,JY,JX)                   !snowpack CaOH2, [mol d-2]
  real(r8) :: ZCACW(JS,JY,JX)                   !snowpack CaCO3, [mol d-2]
  real(r8) :: ZCAHW(JS,JY,JX)                   !snowpack CaHCO3, [mol d-2]
  real(r8) :: ZCASW(JS,JY,JX)                   !snowpack CaSO4, [mol d-2]
  real(r8) :: ZMGOW(JS,JY,JX)                   !snowpack MgOH2, [mol d-2]
  real(r8) :: ZMGCW(JS,JY,JX)                   !snowpack MgCO3, [mol d-2]
  real(r8) :: ZMGHW(JS,JY,JX)                   !snowpack MgHCO3, [mol d-2]
  real(r8) :: ZMGSW(JS,JY,JX)                   !snowpack MgSO4, [mol d-2]
  real(r8) :: ZNACW(JS,JY,JX)                   !snowpack NaCO3, [mol d-2]
  real(r8) :: ZNASW(JS,JY,JX)                   !snowpack NaSO4, [mol d-2]
  real(r8) :: ZKASW(JS,JY,JX)                   !snowpack KSO4, [mol d-2]
  real(r8) :: H0PO4W(JS,JY,JX)                  !snowpack PO4, [mol d-2]
  real(r8) :: H3PO4W(JS,JY,JX)                  !snowpack H3PO4, [mol d-2]
  real(r8) :: ZFE1PW(JS,JY,JX)                  !snowpack FeHPO4, [mol d-2]
  real(r8) :: ZFE2PW(JS,JY,JX)                  !snowpack FeH2PO4, [mol d-2]
  real(r8) :: ZCA0PW(JS,JY,JX)                  !snowpack CaPO4, [mol d-2]
  real(r8) :: ZCA1PW(JS,JY,JX)                  !snowpack CaHPO4, [mol d-2]
  real(r8) :: ZCA2PW(JS,JY,JX)                  !snowpack CaH2PO4, [mol d-2]
  real(r8) :: ZMG1PW(JS,JY,JX)                  !snowpack MgHPO4, [mol d-2]

  real(r8) :: XQSAL(2,JV,JH)                    !total Al in snow drift, [mol d-2 h-1]
  real(r8) :: XQSFE(2,JV,JH)                    !total Fe in snow drift, [mol d-2 h-1]
  real(r8) :: XQSHY(2,JV,JH)                    !total H in snow drift, [mol d-2 h-1]
  real(r8) :: XQSCA(2,JV,JH)                    !total Ca in snow drift, [mol d-2 h-1]
  real(r8) :: XQSMG(2,JV,JH)                    !total Mg in snow drift, [mol d-2 h-1]
  real(r8) :: XQSNA(2,JV,JH)                    !total Na in snow drift, [mol d-2 h-1]
  real(r8) :: XQSKA(2,JV,JH)                    !total K in snow drift, [mol d-2 h-1]
  real(r8) :: XQSOH(2,JV,JH)                    !total OH in snow drift, [mol d-2 h-1]
  real(r8) :: XQSSO(2,JV,JH)                    !total SO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSCL(2,JV,JH)                    !total Cl in snow drift, [mol d-2 h-1]
  real(r8) :: XQSC3(2,JV,JH)                    !total CO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSHC(2,JV,JH)                    !total HCO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSAL1(2,JV,JH)                   !total AlOH in snow drift, [mol d-2 h-1]
  real(r8) :: XQSAL2(2,JV,JH)                   !total AlOH2 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSAL3(2,JV,JH)                   !total AlOH3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSAL4(2,JV,JH)                   !total AlOH4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSALS(2,JV,JH)                   !total AlSO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSFE1(2,JV,JH)                   !total FeOH in snow drift, [mol d-2 h-1]
  real(r8) :: XQSFE2(2,JV,JH)                   !total FeOH2 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSFE3(2,JV,JH)                   !total FeOH3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSFE4(2,JV,JH)                   !total FeOH4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSFES(2,JV,JH)                   !total FeSO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSCAO(2,JV,JH)                   !total CaOH in snow drift, [mol d-2 h-1]
  real(r8) :: XQSCAC(2,JV,JH)                   !total CaCO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSCAH(2,JV,JH)                   !total CaHCO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSCAS(2,JV,JH)                   !total CaSO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSMGO(2,JV,JH)                   !total MgOH in snow drift, [mol d-2 h-1]
  real(r8) :: XQSMGC(2,JV,JH)                   !total MgCO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSMGH(2,JV,JH)                   !total MgHCO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSMGS(2,JV,JH)                   !total MgSO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSNAC(2,JV,JH)                   !total NaCO3 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSNAS(2,JV,JH)                   !total NaSO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSKAS(2,JV,JH)                   !total KSO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSH0P(2,JV,JH)                   !total PO4 in snow drif, [mol d-2 h-1]
  real(r8) :: XQSH1P(2,JV,JH)                   !total HPO4 in snow drift , [mol d-2 h-1]
  real(r8) :: XQSH3P(2,JV,JH)                   !total H3PO4 in snow drift , [mol d-2 h-1]
  real(r8) :: XQSF1P(2,JV,JH)                   !total FeHPO4 in snow drift , [mol d-2 h-1]
  real(r8) :: XQSF2P(2,JV,JH)                   !total FeH2PO4 in snow drift , [mol d-2 h-1]
  real(r8) :: XQSC0P(2,JV,JH)                   !total CaPO4 in snow drift , [mol d-2 h-1]
  real(r8) :: XQSC1P(2,JV,JH)                   !total CaHPO4 in snow drift, [mol d-2 h-1]
  real(r8) :: XQSC2P(2,JV,JH)                   !total CaH2PO4 in snow drift , [mol d-2 h-1]
  real(r8) :: XQSM1P(2,JV,JH)                   !total MgHPO4 in snow drift , [mol d-2 h-1]

end module SnowDataType
