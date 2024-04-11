module UptakePars
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  use ElmIDMod
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
!
!     MaxIterNum=max number of cycles in convergence soln for water uptake
!     DIFFX,DIFFY=acceptance criteria in convergence soln
!     FMN=min PFT:total population ratio
!     MinCanopyBndlResist_pft,RACX=min,max canopy boundary layer resistance (h m-1)
!     RZ=surface resistance to evaporation (h m-1)
!     EMODW=wood modulus of elasticity (MPa)
!     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
!     SNH3X=NH3 solubility at 25 oC (g m-3 water/(g m-3 air))
!     EMMC=canopy emissivity
!     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition on root,myco N,P uptake(g g-1)
!     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation (h-1)
!

  integer :: MaxIterNum
  real(r8) :: DIFFX
  real(r8) :: DIFFY
  real(r8) :: FMN
  real(r8) :: MinCanopyBndlResist_pft
  real(r8) :: RACX
  real(r8) :: RZ
  real(r8) :: DSTK
  real(r8) :: VSTK
  real(r8) :: SNH3X
  real(r8) :: EMMC
  real(r8) :: EMODW
  real(r8) :: ZCKI
  real(r8) :: PCKI
  real(r8) :: ZPKI
  real(r8) :: PZKI
  real(r8) :: FEXUDE(NumPlantChemElms)

  contains
  subroutine InitUptakePars
  implicit none

  MaxIterNum=200
  DIFFX=1.0E-09_r8
  DIFFY=0.5E-02_r8
  FMN=ppmc
  MinCanopyBndlResist_pft=0.00139_r8
  RACX=0.0278_r8
  RZ=0.0139_r8
  DSTK=0.225_r8
  VSTK=ppmc/DSTK
  SNH3X=2.852E+02_r8
  EMMC=0.97_r8
  EMODW=50.0_r8
  ZCKI=0.5E-01_r8
  PCKI=0.5E-02_r8
  ZPKI=ZCKI/PCKI
  PZKI=PCKI/ZCKI
  FEXUDE(ielmc)=0.5E-03_r8
  FEXUDE(ielmn)=1.0E-02_r8
  FEXUDE(ielmp)=1.0E-02_r8
  
  end subroutine InitUptakePars
end module UptakePars
