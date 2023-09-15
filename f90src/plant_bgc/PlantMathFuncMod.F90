module PlantMathFuncMod
!
!DESCRIPTION
! code for small functions used by plant processes
use data_kind_mod, only : r8 => DAT_KIND_R8
use EcoSimConst
use MiniMathMod
implicit none
  
  
contains

  pure function get_FDM(PSICanP)result(FDMP)
  !
  !compute the Ratio of leaf+sheath dry mass to symplasmic water (g g–1)
  !as a function of absolute value of leaf water potential (MPa)
  
  implicit none
  real(r8), intent(in) :: PSICanP   !canopy water potential, MPa
  real(r8) :: APSILT
  real(r8) :: FDMP

  APSILT=ABS(PSICanP)
  FDMP=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)

  end function get_FDM

!--------------------------------------------------------------------------------
  subroutine update_osmo_turg_pressure(PSIO,CCPOLT,OSMO,TKP,PSIOsmo,PSITurg,FDMP1)
  !
  !DESCRIPTION
  !update the osmotic and turgor pressure of a plant organ
  implicit none
  real(r8), intent(in) :: PSIO   !plant orgran pressure, MPa
  real(r8), intent(in) :: CCPOLT !total organ dry mass, C+N+P
  real(r8), intent(in) :: OSMO   !canopy osmotic potential when canopy water potential = 0 MPa
  real(r8), intent(in) :: TKP    !organ temperature, Kelvin
  real(r8), intent(out) :: PSIOsmo  !osmotic pressure of the organ, MPa
  real(r8), intent(out) :: PSITurg  !turgor pressure of the organ, MPa
  real(r8), optional, intent(out) :: FDMP1

  real(r8) :: OSWT
  real(r8) :: FDMP

  FDMP=get_FDM(PSIO)
  if(present(fdmp1))FDMP1=FDMP
  OSWT=36.0_r8+840.0_r8*AZMAX1(CCPOLT)
  PSIOsmo=FDMP/0.16_r8*OSMO-RGAS*TKP*FDMP*CCPOLT/OSWT
  PSITurg=AZMAX1(PSIO-PSIOsmo)

  end subroutine update_osmo_turg_pressure
!--------------------------------------------------------------------------------

  subroutine calc_seed_geometry(SeedCMass,SeedVolume,SeedLength,SeedArea)
  !
  !DESCRIPTION
  !assuming the seed is spherical, compute its volume, diameter(=length), and surface area
  !     SeedVolume,SeedLength,SeedArea=seed volume(m3),length(m),area(m2)
  !     SeedCMass=seed C mass (g) from PFT file
  !

  implicit none
  real(r8), intent(in)  :: SeedCMass   !carbon mass per seed
  real(r8), intent(out) :: SeedVolume,SeedLength,SeedArea

  SeedVolume=SeedCMass*5.0E-06_r8
  SeedLength=2.0_r8*(0.75_r8*SeedVolume/PICON)**0.33_r8
  SeedArea=4.0_r8*PICON*(SeedLength/2.0_r8)**2_r8

  end subroutine calc_seed_geometry
!--------------------------------------------------------------------------------
  pure function calc_root_grow_tempf(TKSO)result(fT_root)
  !
  !DESCRIPTION
  !compute the temperature dependence for plant root growth
  implicit none
  real(r8), intent(in) :: TKSO   !apparent temperature felt by the root
  real(r8) :: fT_root
  real(r8) :: RTK,STK,ACTV

  RTK=RGAS*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  FT_ROOT=EXP(25.229_r8-62500._r8/RTK)/ACTV

  end function calc_root_grow_tempf
!--------------------------------------------------------------------------------
  pure function calc_canopy_grow_tempf(TKGO)result(fT_canp)
  !
  !DESCRIPTION
  !compute the temperature dependence for plant canopy growth
  implicit none
  real(r8), intent(in) :: TKGO   !apparent temperature felt by the canopy
  real(r8) :: fT_canp
  real(r8) :: RTK,STK,ACTV

  RTK=RGAS*TKGO
  STK=710.0_r8*TKGO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  FT_canp=EXP(25.229_r8-62500._r8/RTK)/ACTV

  end function calc_canopy_grow_tempf

end module PlantMathFuncMod
