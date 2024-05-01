module ChemTranspDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc
  use TracerIDMod
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  TScal4Difsvity_vr(:,:,:)                        !temperature effect on diffusivity
  real(r8),target,allocatable ::  DISP(:,:,:,:)                      !aqueous dispersivity

  real(r8),target,allocatable ::  Gas_Disol_Flx_vr(:,:,:,:)                 !Gas dissolution flux
  real(r8),target,allocatable ::  GasDifc_vr(:,:,:,:)                   !gaseous diffusivity [m2 h-1]
  real(r8),target,allocatable ::  SolDifc_vr(:,:,:,:)
  real(r8),target,allocatable ::  CGSGL(:,:,:)                       !gaseous CO2 diffusivity	[m2 h-1]
  real(r8),target,allocatable ::  CLSGL(:,:,:)                       !aqueous CO2 diffusivity	[m2 h-1]
  real(r8),target,allocatable ::  OGSGL(:,:,:)                       !gaseous O2 diffusivity	m2 h-1
  real(r8),target,allocatable ::  O2AquaDiffusvity(:,:,:)                       !aqueous CO2 diffusivity	m2 h-1
  real(r8),target,allocatable ::  ZGSGL(:,:,:)                       !gaseous N2 diffusivity	m2 h-1
  real(r8),target,allocatable ::  CHSGL(:,:,:)                       !gaseous CH4 diffusivity	m2 h-1
  real(r8),target,allocatable ::  CQSGL(:,:,:)                       !aqueous CH4 diffusivity	m2 h-1
  real(r8),target,allocatable ::  ZLSGL(:,:,:)                       !aqueous N2 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  ZHSGL(:,:,:)                       !aqueous NH4 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  ZNSGL(:,:,:)                       !aqueous NH3 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  ZOSGL(:,:,:)                       !aqueous NO3 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  POSGL(:,:,:)                       !aqueous PO4 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  OCSGL(:,:,:)                       !aqueous DOC diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  ONSGL(:,:,:)                       !aqueous DON diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  OPSGL(:,:,:)                       !aqueous DOP diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  OASGL(:,:,:)                       !aqueous acetate diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  Z2SGL(:,:,:)                       !gaseous N2O diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  ZVSGL(:,:,:)                       !aqueous N2O diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  WGSGL(:,:,:)                       !water vapor diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  H2OVapDifscSno(:,:,:)                       !water vapor diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  VaporDiffusivityLitR(:,:)                         !water vapor diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  WGSGA(:,:)                         !water vapor diffusivity, [m2 h-1]

  real(r8),target,allocatable ::  GasSolbility_vr(:,:,:,:)                !solubility of gases
  real(r8),target,allocatable ::  HGSGL(:,:,:)                       !gaseous H2 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  HLSGL(:,:,:)                       !aqueous H2 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  XCODFG(:,:,:)                      !soil CO2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XCHDFG(:,:,:)                      !soil CH4 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XOXDFG(:,:,:)                      !soil O2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XNGDFG(:,:,:)                      !soil N2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XN2DFG(:,:,:)                      !soil N2O dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XN3DFG(:,:,:)                      !soil NH3 dissolution (+ve) - volatilization (-ve) non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNBDFG(:,:,:)                      !soil NH3 dissolution (+ve) - volatilization (-ve) band, [g d-2 h-1]
  real(r8),target,allocatable ::  XHGDFG(:,:,:)                      !soil H2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2F(:,:,:)                       !net gaseous CO2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RCH4L(:,:,:)                       !net aqueous CH4 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RO2GasXchangePrev_vr(:,:,:)                       !net gaseous O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RO2AquaXchangePrev_vr(:,:,:)                       !net aqueous O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RCH4F(:,:,:)                       !net gaseous CH4 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  ALSGL(:,:,:)                       !aqueous Al diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  FESGL(:,:,:)                       !aqueous Fe diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  HYSGL(:,:,:)                       !aqueous H diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  CASGL(:,:,:)                       !aqueous Ca diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  GMSGL(:,:,:)                       !aqueous Mg diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  ANSGL(:,:,:)                       !aqueous Na diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  AKSGL(:,:,:)                       !aqueous K diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  OHSGL(:,:,:)                       !aqueous OH diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  C3SGL(:,:,:)                       !aqueous CO3 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  HCSGL(:,:,:)                       !aqueous HCO3 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  SOSGL(:,:,:)                       !aqueous SO4 diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  CLSXL(:,:,:)                       !aqueous Cl diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  trc_salt_rof_bounds(:,:,:,:,:)                     !total Al in runoff, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcg_2DFloXSurRunoff(:,:,:,:,:)                    !surface runoff gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcn_2DFloXSurRunoff(:,:,:,:,:)                    !surface runoff nutrient flux, [g d-2 h-1]
  real(r8),target,allocatable ::  dom_2DFloXSurRunoff(:,:,:,:,:,:)                  !surface runoff DOC flux, [g d-2 h-1]

  private :: InitAllocate

  contains


  subroutine InitChemTranspData

  implicit none

  call InitAllocate

  end subroutine InitChemTranspData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none
  allocate(TScal4Difsvity_vr(0:JZ,JY,JX));   TScal4Difsvity_vr=0._r8
  allocate(DISP(3,JD,JV,JH));   DISP=0._r8

  allocate(GasDifc_vr(idg_beg:idg_end,JZ,JY,JX));GasDifc_vr=0._r8
  allocate(SolDifc_vr(ids_beg:ids_end,0:JZ,JY,JX));SolDifc_vr=0._r8
  allocate(CGSGL(JZ,JY,JX));    CGSGL=0._r8
  allocate(CLSGL(0:JZ,JY,JX));  CLSGL=0._r8
  allocate(OGSGL(JZ,JY,JX));    OGSGL=0._r8
  allocate(O2AquaDiffusvity(0:JZ,JY,JX));  O2AquaDiffusvity=0._r8
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
  allocate(H2OVapDifscSno(JS,JY,JX));    H2OVapDifscSno=0._r8
  allocate(VaporDiffusivityLitR(JY,JX));       VaporDiffusivityLitR=0._r8
  allocate(WGSGA(JY,JX));       WGSGA=0._r8
  allocate(GasSolbility_vr(idg_beg:idg_end,0:JZ,JY,JX)); GasSolbility_vr=0._r8

  allocate(Gas_Disol_Flx_vr(idg_beg:idg_end,0:JZ,JY,JX)); Gas_Disol_Flx_vr=0._r8
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
  allocate(RO2GasXchangePrev_vr(0:JZ,JY,JX));  RO2GasXchangePrev_vr=0._r8
  allocate(RO2AquaXchangePrev_vr(0:JZ,JY,JX));  RO2AquaXchangePrev_vr=0._r8
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
  allocate(trc_salt_rof_bounds(idsalt_beg:idsalt_end,2,2,JV,JH));   trc_salt_rof_bounds=0._r8
  allocate(trcg_2DFloXSurRunoff(idg_beg:idg_end-1,2,2,JV,JH));  trcg_2DFloXSurRunoff=0._r8
  allocate(trcn_2DFloXSurRunoff(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_2DFloXSurRunoff=0._r8
  allocate(dom_2DFloXSurRunoff(idom_beg:idom_end,1:jcplx,2,2,JV,JH));dom_2DFloXSurRunoff=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructChemTranspData
  use abortutils, only : destroy

  implicit none

  call destroy(trc_salt_rof_bounds)
  call destroy(GasDifc_vr)
  call destroy(SolDifc_vr)
  call destroy(TScal4Difsvity_vr)
  call destroy(DISP)
  call destroy(CGSGL)
  call destroy(CLSGL)
  call destroy(OGSGL)
  call destroy(O2AquaDiffusvity)
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
  call destroy(H2OVapDifscSno)
  call destroy(VaporDiffusivityLitR)
  call destroy(WGSGA)

  call destroy(trcn_2DFloXSurRunoff)
  call destroy(trcg_2DFloXSurRunoff)
  call destroy(GasSolbility_vr)

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
  call destroy(RO2GasXchangePrev_vr)
  call destroy(RO2AquaXchangePrev_vr)
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
  call destroy(dom_2DFloXSurRunoff)
  end subroutine DestructChemTranspData
end module ChemTranspDataType
