module ChemTranspDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8),allocatable ::  TFND(:,:,:)                        !temperature effect on diffusivity
  real(r8),allocatable ::  DISP(:,:,:,:)                      !aqueous dispersivity
  real(r8),allocatable ::  CGSGL(:,:,:)                       !gaseous CO2 diffusivity	[m2 h-1]
  real(r8),allocatable ::  CLSGL(:,:,:)                       !aqueous CO2 diffusivity	[m2 h-1]
  real(r8),allocatable ::  OGSGL(:,:,:)                       !gaseous O2 diffusivity	m2 h-1
  real(r8),allocatable ::  OLSGL(:,:,:)                       !aqueous CO2 diffusivity	m2 h-1
  real(r8),allocatable ::  ZGSGL(:,:,:)                       !gaseous N2 diffusivity	m2 h-1
  real(r8),allocatable ::  CHSGL(:,:,:)                       !gaseous CH4 diffusivity	m2 h-1
  real(r8),allocatable ::  CQSGL(:,:,:)                       !aqueous CH4 diffusivity	m2 h-1
  real(r8),allocatable ::  ZLSGL(:,:,:)                       !aqueous N2 diffusivity, [m2 h-1]
  real(r8),allocatable ::  ZHSGL(:,:,:)                       !aqueous NH4 diffusivity, [m2 h-1]
  real(r8),allocatable ::  ZNSGL(:,:,:)                       !aqueous NH3 diffusivity, [m2 h-1]
  real(r8),allocatable ::  ZOSGL(:,:,:)                       !aqueous NO3 diffusivity, [m2 h-1]
  real(r8),allocatable ::  POSGL(:,:,:)                       !aqueous PO4 diffusivity, [m2 h-1]
  real(r8),allocatable ::  OCSGL(:,:,:)                       !aqueous DOC diffusivity, [m2 h-1]
  real(r8),allocatable ::  ONSGL(:,:,:)                       !aqueous DON diffusivity, [m2 h-1]
  real(r8),allocatable ::  OPSGL(:,:,:)                       !aqueous DOP diffusivity, [m2 h-1]
  real(r8),allocatable ::  OASGL(:,:,:)                       !aqueous acetate diffusivity, [m2 h-1]
  real(r8),allocatable ::  Z2SGL(:,:,:)                       !gaseous N2O diffusivity, [m2 h-1]
  real(r8),allocatable ::  ZVSGL(:,:,:)                       !aqueous N2O diffusivity, [m2 h-1]
  real(r8),allocatable ::  WGSGL(:,:,:)                       !water vapor diffusivity, [m2 h-1]
  real(r8),allocatable ::  WGSGW(:,:,:)                       !water vapor diffusivity, [m2 h-1]
  real(r8),allocatable ::  WGSGR(:,:)                         !water vapor diffusivity, [m2 h-1]
  real(r8),allocatable ::  WGSGA(:,:)                         !water vapor diffusivity, [m2 h-1]
  real(r8),allocatable ::  SCO2L(:,:,:)                       !solubility of CO2, [m3 m-3]
  real(r8),allocatable ::  SOXYL(:,:,:)                       !solubility of O2, [m3 m-3]
  real(r8),allocatable ::  SCH4L(:,:,:)                       !solubility of CH4, [m3 m-3]
  real(r8),allocatable ::  SN2OL(:,:,:)                       !solubility of N2O, [m3 m-3]
  real(r8),allocatable ::  SN2GL(:,:,:)                       !solubility of N2, [m3 m-3]
  real(r8),allocatable ::  SNH3L(:,:,:)                       !solubility of NH3, [m3 m-3]
  real(r8),allocatable ::  SH2GL(:,:,:)                       !solubility of H2, [m3 m-3]
  real(r8),allocatable ::  HGSGL(:,:,:)                       !gaseous H2 diffusivity, [m2 h-1]
  real(r8),allocatable ::  HLSGL(:,:,:)                       !aqueous H2 diffusivity, [m2 h-1]
  real(r8),allocatable ::  XCODFG(:,:,:)                      !soil CO2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XCHDFG(:,:,:)                      !soil CH4 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XOXDFG(:,:,:)                      !soil O2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XNGDFG(:,:,:)                      !soil N2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XN2DFG(:,:,:)                      !soil N2O dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XN3DFG(:,:,:)                      !soil NH3 dissolution (+ve) - volatilization (-ve) non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNBDFG(:,:,:)                      !soil NH3 dissolution (+ve) - volatilization (-ve) band, [g d-2 h-1]
  real(r8),allocatable ::  XHGDFG(:,:,:)                      !soil H2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  RCO2F(:,:,:)                       !net gaseous CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  RCH4L(:,:,:)                       !net aqueous CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  ROXYF(:,:,:)                       !net gaseous O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  ROXYL(:,:,:)                       !net aqueous O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  RCH4F(:,:,:)                       !net gaseous CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  ALSGL(:,:,:)                       !aqueous Al diffusivity, [m2 h-1]
  real(r8),allocatable ::  FESGL(:,:,:)                       !aqueous Fe diffusivity, [m2 h-1]
  real(r8),allocatable ::  HYSGL(:,:,:)                       !aqueous H diffusivity, [m2 h-1]
  real(r8),allocatable ::  CASGL(:,:,:)                       !aqueous Ca diffusivity, [m2 h-1]
  real(r8),allocatable ::  GMSGL(:,:,:)                       !aqueous Mg diffusivity, [m2 h-1]
  real(r8),allocatable ::  ANSGL(:,:,:)                       !aqueous Na diffusivity, [m2 h-1]
  real(r8),allocatable ::  AKSGL(:,:,:)                       !aqueous K diffusivity, [m2 h-1]
  real(r8),allocatable ::  OHSGL(:,:,:)                       !aqueous OH diffusivity, [m2 h-1]
  real(r8),allocatable ::  C3SGL(:,:,:)                       !aqueous CO3 diffusivity, [m2 h-1]
  real(r8),allocatable ::  HCSGL(:,:,:)                       !aqueous HCO3 diffusivity, [m2 h-1]
  real(r8),allocatable ::  SOSGL(:,:,:)                       !aqueous SO4 diffusivity, [m2 h-1]
  real(r8),allocatable ::  CLSXL(:,:,:)                       !aqueous Cl diffusivity, [m2 h-1]
  real(r8),allocatable ::  XQRAL(:,:,:,:)                     !total Al in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRFE(:,:,:,:)                     !total Fe in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRHY(:,:,:,:)                     !total H in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRCA(:,:,:,:)                     !total Ca in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRMG(:,:,:,:)                     !total Mg in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRNA(:,:,:,:)                     !total Na in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRKA(:,:,:,:)                     !total K in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQROH(:,:,:,:)                     !total OH in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRSO(:,:,:,:)                     !total SO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRCL(:,:,:,:)                     !total Cl in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRC3(:,:,:,:)                     !total CO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRHC(:,:,:,:)                     !total HCO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRAL1(:,:,:,:)                    !total AlOH in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRAL2(:,:,:,:)                    !total AlOH2 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRAL3(:,:,:,:)                    !total AlOH3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRAL4(:,:,:,:)                    !total AlOH4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRALS(:,:,:,:)                    !total AlSO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRFE1(:,:,:,:)                    !total FeOH in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRFE2(:,:,:,:)                    !total FeOH2 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRFE3(:,:,:,:)                    !total FeOH3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRFE4(:,:,:,:)                    !total FeOH4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRFES(:,:,:,:)                    !total FeSO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRCAO(:,:,:,:)                    !total CaOH in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRCAC(:,:,:,:)                    !total CaCO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRCAH(:,:,:,:)                    !total CaHCO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRCAS(:,:,:,:)                    !total CaSO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRMGO(:,:,:,:)                    !total MgOH in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRMGC(:,:,:,:)                    !total MgCO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRMGH(:,:,:,:)                    !total MgHCO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRMGS(:,:,:,:)                    !total MgSO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRNAC(:,:,:,:)                    !total NaCO3 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRNAS(:,:,:,:)                    !total NaSO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRKAS(:,:,:,:)                    !total KSO4 in runoff, [mol d-2 h-1]
  real(r8),allocatable ::  XQRH0P(:,:,:,:)                    !total PO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRH3P(:,:,:,:)                    !total H3PO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRF1P(:,:,:,:)                    !total FeHPO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRF2P(:,:,:,:)                    !total FeH2PO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRC0P(:,:,:,:)                    !total CaPO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRC1P(:,:,:,:)                    !total CaHPO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRC2P(:,:,:,:)                    !total CaH2PO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XQRM1P(:,:,:,:)                    !total MgHPO4 in runoff non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XCOQRS(:,:,:,:)                    !surface runoff CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XCHQRS(:,:,:,:)                    !surface brunoff CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  XOXQRS(:,:,:,:)                    !surface runoff O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XNGQRS(:,:,:,:)                    !surface runoff N2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XN2QRS(:,:,:,:)                    !surface runoff N2O flux, [g d-2 h-1]
  real(r8),allocatable ::  XHGQRS(:,:,:,:)                    !surface runoff H2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XN4QRW(:,:,:,:)                    !surface runoff NH4 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  XN3QRW(:,:,:,:)                    !surface runoff NH3 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNOQRW(:,:,:,:)                    !surface runoff NO3 flux non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNXQRS(:,:,:,:)                    !surface runoff NO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XP4QRW(:,:,:,:)                    !surface runoff PO4 flux, [g d-2 h-1]
  real(r8),allocatable ::  XP1QRW(:,:,:,:)                    !surface runoff HPO4 flux, [g d-2 h-1]
  real(r8),allocatable ::  XOCQRS(:,:,:,:,:)                  !surface runoff DOC flux, [g d-2 h-1]
  real(r8),allocatable ::  XONQRS(:,:,:,:,:)                  !surface runoff DON flux, [g d-2 h-1]
  real(r8),allocatable ::  XOPQRS(:,:,:,:,:)                  !surface runoff DOP flux, [g d-2 h-1]
  real(r8),allocatable ::  XOAQRS(:,:,:,:,:)                  !surface runoff acetate flux, [g d-2 h-1]

  private :: InitAllocate

  contains


  subroutine InitChemTranspData

  implicit none

  call InitAllocate

  end subroutine InitChemTranspData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none
  allocate(TFND(0:JZ,JY,JX));   TFND=0._r8
  allocate(DISP(3,JD,JV,JH));   DISP=0._r8
  allocate(CGSGL(JZ,JY,JX));    CGSGL=0._r8
  allocate(CLSGL(0:JZ,JY,JX));  CLSGL=0._r8
  allocate(OGSGL(JZ,JY,JX));    OGSGL=0._r8
  allocate(OLSGL(0:JZ,JY,JX));  OLSGL=0._r8
  allocate(ZGSGL(JZ,JY,JX));    ZGSGL=0._r8
  allocate(CHSGL(JZ,JY,JX));    CHSGL=0._r8
  allocate(CQSGL(0:JZ,JY,JX));  CQSGL=0._r8
  allocate(ZLSGL(0:JZ,JY,JX));  ZLSGL=0._r8
  allocate(ZHSGL(JZ,JY,JX));    ZHSGL=0._r8
  allocate(ZNSGL(0:JZ,JY,JX));  ZNSGL=0._r8
  allocate(ZOSGL(0:JZ,JY,JX));  ZOSGL=0._r8
  allocate(POSGL(0:JZ,JY,JX));  POSGL=0._r8
  allocate(OCSGL(0:JZ,JY,JX));  OCSGL=0._r8
  allocate(ONSGL(0:JZ,JY,JX));  ONSGL=0._r8
  allocate(OPSGL(0:JZ,JY,JX));  OPSGL=0._r8
  allocate(OASGL(0:JZ,JY,JX));  OASGL=0._r8
  allocate(Z2SGL(JZ,JY,JX));    Z2SGL=0._r8
  allocate(ZVSGL(0:JZ,JY,JX));  ZVSGL=0._r8
  allocate(WGSGL(JZ,JY,JX));    WGSGL=0._r8
  allocate(WGSGW(JS,JY,JX));    WGSGW=0._r8
  allocate(WGSGR(JY,JX));       WGSGR=0._r8
  allocate(WGSGA(JY,JX));       WGSGA=0._r8
  allocate(SCO2L(0:JZ,JY,JX));  SCO2L=0._r8
  allocate(SOXYL(0:JZ,JY,JX));  SOXYL=0._r8
  allocate(SCH4L(0:JZ,JY,JX));  SCH4L=0._r8
  allocate(SN2OL(0:JZ,JY,JX));  SN2OL=0._r8
  allocate(SN2GL(0:JZ,JY,JX));  SN2GL=0._r8
  allocate(SNH3L(0:JZ,JY,JX));  SNH3L=0._r8
  allocate(SH2GL(0:JZ,JY,JX));  SH2GL=0._r8
  allocate(HGSGL(JZ,JY,JX));    HGSGL=0._r8
  allocate(HLSGL(0:JZ,JY,JX));  HLSGL=0._r8
  allocate(XCODFG(0:JZ,JY,JX)); XCODFG=0._r8
  allocate(XCHDFG(0:JZ,JY,JX)); XCHDFG=0._r8
  allocate(XOXDFG(0:JZ,JY,JX)); XOXDFG=0._r8
  allocate(XNGDFG(0:JZ,JY,JX)); XNGDFG=0._r8
  allocate(XN2DFG(0:JZ,JY,JX)); XN2DFG=0._r8
  allocate(XN3DFG(0:JZ,JY,JX)); XN3DFG=0._r8
  allocate(XNBDFG(0:JZ,JY,JX)); XNBDFG=0._r8
  allocate(XHGDFG(0:JZ,JY,JX)); XHGDFG=0._r8
  allocate(RCO2F(0:JZ,JY,JX));  RCO2F=0._r8
  allocate(RCH4L(0:JZ,JY,JX));  RCH4L=0._r8
  allocate(ROXYF(0:JZ,JY,JX));  ROXYF=0._r8
  allocate(ROXYL(0:JZ,JY,JX));  ROXYL=0._r8
  allocate(RCH4F(0:JZ,JY,JX));  RCH4F=0._r8
  allocate(ALSGL(JZ,JY,JX));    ALSGL=0._r8
  allocate(FESGL(JZ,JY,JX));    FESGL=0._r8
  allocate(HYSGL(JZ,JY,JX));    HYSGL=0._r8
  allocate(CASGL(JZ,JY,JX));    CASGL=0._r8
  allocate(GMSGL(JZ,JY,JX));    GMSGL=0._r8
  allocate(ANSGL(JZ,JY,JX));    ANSGL=0._r8
  allocate(AKSGL(JZ,JY,JX));    AKSGL=0._r8
  allocate(OHSGL(JZ,JY,JX));    OHSGL=0._r8
  allocate(C3SGL(JZ,JY,JX));    C3SGL=0._r8
  allocate(HCSGL(JZ,JY,JX));    HCSGL=0._r8
  allocate(SOSGL(JZ,JY,JX));    SOSGL=0._r8
  allocate(CLSXL(JZ,JY,JX));    CLSXL=0._r8
  allocate(XQRAL(2,2,JV,JH));   XQRAL=0._r8
  allocate(XQRFE(2,2,JV,JH));   XQRFE=0._r8
  allocate(XQRHY(2,2,JV,JH));   XQRHY=0._r8
  allocate(XQRCA(2,2,JV,JH));   XQRCA=0._r8
  allocate(XQRMG(2,2,JV,JH));   XQRMG=0._r8
  allocate(XQRNA(2,2,JV,JH));   XQRNA=0._r8
  allocate(XQRKA(2,2,JV,JH));   XQRKA=0._r8
  allocate(XQROH(2,2,JV,JH));   XQROH=0._r8
  allocate(XQRSO(2,2,JV,JH));   XQRSO=0._r8
  allocate(XQRCL(2,2,JV,JH));   XQRCL=0._r8
  allocate(XQRC3(2,2,JV,JH));   XQRC3=0._r8
  allocate(XQRHC(2,2,JV,JH));   XQRHC=0._r8
  allocate(XQRAL1(2,2,JV,JH));  XQRAL1=0._r8
  allocate(XQRAL2(2,2,JV,JH));  XQRAL2=0._r8
  allocate(XQRAL3(2,2,JV,JH));  XQRAL3=0._r8
  allocate(XQRAL4(2,2,JV,JH));  XQRAL4=0._r8
  allocate(XQRALS(2,2,JV,JH));  XQRALS=0._r8
  allocate(XQRFE1(2,2,JV,JH));  XQRFE1=0._r8
  allocate(XQRFE2(2,2,JV,JH));  XQRFE2=0._r8
  allocate(XQRFE3(2,2,JV,JH));  XQRFE3=0._r8
  allocate(XQRFE4(2,2,JV,JH));  XQRFE4=0._r8
  allocate(XQRFES(2,2,JV,JH));  XQRFES=0._r8
  allocate(XQRCAO(2,2,JV,JH));  XQRCAO=0._r8
  allocate(XQRCAC(2,2,JV,JH));  XQRCAC=0._r8
  allocate(XQRCAH(2,2,JV,JH));  XQRCAH=0._r8
  allocate(XQRCAS(2,2,JV,JH));  XQRCAS=0._r8
  allocate(XQRMGO(2,2,JV,JH));  XQRMGO=0._r8
  allocate(XQRMGC(2,2,JV,JH));  XQRMGC=0._r8
  allocate(XQRMGH(2,2,JV,JH));  XQRMGH=0._r8
  allocate(XQRMGS(2,2,JV,JH));  XQRMGS=0._r8
  allocate(XQRNAC(2,2,JV,JH));  XQRNAC=0._r8
  allocate(XQRNAS(2,2,JV,JH));  XQRNAS=0._r8
  allocate(XQRKAS(2,2,JV,JH));  XQRKAS=0._r8
  allocate(XQRH0P(2,2,JV,JH));  XQRH0P=0._r8
  allocate(XQRH3P(2,2,JV,JH));  XQRH3P=0._r8
  allocate(XQRF1P(2,2,JV,JH));  XQRF1P=0._r8
  allocate(XQRF2P(2,2,JV,JH));  XQRF2P=0._r8
  allocate(XQRC0P(2,2,JV,JH));  XQRC0P=0._r8
  allocate(XQRC1P(2,2,JV,JH));  XQRC1P=0._r8
  allocate(XQRC2P(2,2,JV,JH));  XQRC2P=0._r8
  allocate(XQRM1P(2,2,JV,JH));  XQRM1P=0._r8
  allocate(XCOQRS(2,2,JV,JH));  XCOQRS=0._r8
  allocate(XCHQRS(2,2,JV,JH));  XCHQRS=0._r8
  allocate(XOXQRS(2,2,JV,JH));  XOXQRS=0._r8
  allocate(XNGQRS(2,2,JV,JH));  XNGQRS=0._r8
  allocate(XN2QRS(2,2,JV,JH));  XN2QRS=0._r8
  allocate(XHGQRS(2,2,JV,JH));  XHGQRS=0._r8
  allocate(XN4QRW(2,2,JV,JH));  XN4QRW=0._r8
  allocate(XN3QRW(2,2,JV,JH));  XN3QRW=0._r8
  allocate(XNOQRW(2,2,JV,JH));  XNOQRW=0._r8
  allocate(XNXQRS(2,2,JV,JH));  XNXQRS=0._r8
  allocate(XP4QRW(2,2,JV,JH));  XP4QRW=0._r8
  allocate(XP1QRW(2,2,JV,JH));  XP1QRW=0._r8
  allocate(XOCQRS(0:jcplx1,2,2,JV,JH));XOCQRS=0._r8
  allocate(XONQRS(0:jcplx1,2,2,JV,JH));XONQRS=0._r8
  allocate(XOPQRS(0:jcplx1,2,2,JV,JH));XOPQRS=0._r8
  allocate(XOAQRS(0:jcplx1,2,2,JV,JH));XOAQRS=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructChemTranspData
  use abortutils, only : destroy

  implicit none
  call destroy(TFND)
  call destroy(DISP)
  call destroy(CGSGL)
  call destroy(CLSGL)
  call destroy(OGSGL)
  call destroy(OLSGL)
  call destroy(ZGSGL)
  call destroy(CHSGL)
  call destroy(CQSGL)
  call destroy(ZLSGL)
  call destroy(ZHSGL)
  call destroy(ZNSGL)
  call destroy(ZOSGL)
  call destroy(POSGL)
  call destroy(OCSGL)
  call destroy(ONSGL)
  call destroy(OPSGL)
  call destroy(OASGL)
  call destroy(Z2SGL)
  call destroy(ZVSGL)
  call destroy(WGSGL)
  call destroy(WGSGW)
  call destroy(WGSGR)
  call destroy(WGSGA)
  call destroy(SCO2L)
  call destroy(SOXYL)
  call destroy(SCH4L)
  call destroy(SN2OL)
  call destroy(SN2GL)
  call destroy(SNH3L)
  call destroy(SH2GL)
  call destroy(HGSGL)
  call destroy(HLSGL)
  call destroy(XCODFG)
  call destroy(XCHDFG)
  call destroy(XOXDFG)
  call destroy(XNGDFG)
  call destroy(XN2DFG)
  call destroy(XN3DFG)
  call destroy(XNBDFG)
  call destroy(XHGDFG)
  call destroy(RCO2F)
  call destroy(RCH4L)
  call destroy(ROXYF)
  call destroy(ROXYL)
  call destroy(RCH4F)
  call destroy(ALSGL)
  call destroy(FESGL)
  call destroy(HYSGL)
  call destroy(CASGL)
  call destroy(GMSGL)
  call destroy(ANSGL)
  call destroy(AKSGL)
  call destroy(OHSGL)
  call destroy(C3SGL)
  call destroy(HCSGL)
  call destroy(SOSGL)
  call destroy(CLSXL)
  call destroy(XQRAL)
  call destroy(XQRFE)
  call destroy(XQRHY)
  call destroy(XQRCA)
  call destroy(XQRMG)
  call destroy(XQRNA)
  call destroy(XQRKA)
  call destroy(XQROH)
  call destroy(XQRSO)
  call destroy(XQRCL)
  call destroy(XQRC3)
  call destroy(XQRHC)
  call destroy(XQRAL1)
  call destroy(XQRAL2)
  call destroy(XQRAL3)
  call destroy(XQRAL4)
  call destroy(XQRALS)
  call destroy(XQRFE1)
  call destroy(XQRFE2)
  call destroy(XQRFE3)
  call destroy(XQRFE4)
  call destroy(XQRFES)
  call destroy(XQRCAO)
  call destroy(XQRCAC)
  call destroy(XQRCAH)
  call destroy(XQRCAS)
  call destroy(XQRMGO)
  call destroy(XQRMGC)
  call destroy(XQRMGH)
  call destroy(XQRMGS)
  call destroy(XQRNAC)
  call destroy(XQRNAS)
  call destroy(XQRKAS)
  call destroy(XQRH0P)
  call destroy(XQRH3P)
  call destroy(XQRF1P)
  call destroy(XQRF2P)
  call destroy(XQRC0P)
  call destroy(XQRC1P)
  call destroy(XQRC2P)
  call destroy(XQRM1P)
  call destroy(XCOQRS)
  call destroy(XCHQRS)
  call destroy(XOXQRS)
  call destroy(XNGQRS)
  call destroy(XN2QRS)
  call destroy(XHGQRS)
  call destroy(XN4QRW)
  call destroy(XN3QRW)
  call destroy(XNOQRW)
  call destroy(XNXQRS)
  call destroy(XP4QRW)
  call destroy(XP1QRW)
  call destroy(XOCQRS)
  call destroy(XONQRS)
  call destroy(XOPQRS)
  call destroy(XOAQRS)

  end subroutine DestructChemTranspData
end module ChemTranspDataType
