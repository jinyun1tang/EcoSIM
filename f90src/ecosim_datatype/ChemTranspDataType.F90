module ChemTranspDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: TFND(0:JZ,JY,JX)                  !temperature effect on diffusivity
  real(r8) :: DISP(3,JD,JV,JH)                  !aqueous dispersivity
  real(r8) :: CGSGL(JZ,JY,JX)                   !gaseous CO2 diffusivity	[m2 h-1]
  real(r8) :: CLSGL(0:JZ,JY,JX)                 !aqueous CO2 diffusivity	[m2 h-1]
  real(r8) :: OGSGL(JZ,JY,JX)                   !gaseous O2 diffusivity	m2 h-1
  real(r8) :: OLSGL(0:JZ,JY,JX)                 !aqueous CO2 diffusivity	m2 h-1
  real(r8) :: ZGSGL(JZ,JY,JX)                   !gaseous N2 diffusivity	m2 h-1
  real(r8) :: CHSGL(JZ,JY,JX)                   !gaseous CH4 diffusivity	m2 h-1
  real(r8) :: CQSGL(0:JZ,JY,JX)                 !aqueous CH4 diffusivity	m2 h-1
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



  real(r8) :: XCOQRS(2,2,JV,JH)                 !surface runoff CO2 flux, [g d-2 h-1]
  real(r8) :: XCHQRS(2,2,JV,JH)                 !surface brunoff CH4 flux, [g d-2 h-1]
  real(r8) :: XOXQRS(2,2,JV,JH)                 !surface runoff O2 flux, [g d-2 h-1]
  real(r8) :: XNGQRS(2,2,JV,JH)                 !surface runoff N2 flux, [g d-2 h-1]
  real(r8) :: XN2QRS(2,2,JV,JH)                 !surface runoff N2O flux, [g d-2 h-1]
  real(r8) :: XHGQRS(2,2,JV,JH)                 !surface runoff H2 flux, [g d-2 h-1]
  real(r8) :: XN4QRW(2,2,JV,JH)                 !surface runoff NH4 flux non-band, [g d-2 h-1]
  real(r8) :: XN3QRW(2,2,JV,JH)                 !surface runoff NH3 flux non-band, [g d-2 h-1]
  real(r8) :: XNOQRW(2,2,JV,JH)                 !surface runoff NO3 flux non-band, [g d-2 h-1]
  real(r8) :: XNXQRS(2,2,JV,JH)                 !surface runoff NO2 flux, [g d-2 h-1]
  real(r8) :: XP4QRW(2,2,JV,JH)                 !surface runoff PO4 flux, [g d-2 h-1]
  real(r8) :: XP1QRW(2,2,JV,JH)                 !surface runoff HPO4 flux, [g d-2 h-1]
  real(r8) :: XOCQRS(0:jcplx1,2,2,JV,JH)             !surface runoff DOC flux, [g d-2 h-1]
  real(r8) :: XONQRS(0:jcplx1,2,2,JV,JH)             !surface runoff DON flux, [g d-2 h-1]
  real(r8) :: XOPQRS(0:jcplx1,2,2,JV,JH)             !surface runoff DOP flux, [g d-2 h-1]
  real(r8) :: XOAQRS(0:jcplx1,2,2,JV,JH)             !surface runoff acetate flux, [g d-2 h-1]
  private :: InitAllocate

  contains


  subroutine InitChemTranspData

  implicit none

  call InitAllocate

  end subroutine InitChemTranspData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructChemTranspData

  implicit none

  end subroutine DestructChemTranspData
end module ChemTranspDataType
