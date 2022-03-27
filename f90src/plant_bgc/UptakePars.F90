module UptakePars
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
!
!     MXN=max number of cycles in convergence soln for water uptake
!     DIFFX,DIFFY=acceptance criteria in convergence soln
!     FMN=min PFT:total population ratio
!     RACM,RACX=min,max canopy boundary layer resistance (h m-1)
!     RZ=surface resistance to evaporation (h m-1)
!     EMODW=wood modulus of elasticity (MPa)
!     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
!     SNH3X=NH3 solubility at 25 oC (g m-3 water/(g m-3 air))
!     EMMC=canopy emissivity
!     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition on root,myco N,P uptake(g g-1)
!     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation (h-1)
!

  integer :: MXN
  real(r8) :: DIFFX
  real(r8) :: DIFFY
  real(r8) :: FMN
  real(r8) :: RACM
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
  real(r8) :: FEXUC
  real(r8) :: FEXUN
  real(r8) :: FEXUP

  contains
  subroutine InitUptakePars
  implicit none

  MXN=200
  DIFFX=1.0E-09
  DIFFY=0.5E-02
  FMN=1.0E-06
  RACM=0.00139
  RACX=0.0278
  RZ=0.0139
  DSTK=0.225
  VSTK=1.0E-06/DSTK
  SNH3X=2.852E+02
  EMMC=0.97
  EMODW=50.0
  ZCKI=0.5E-01
  PCKI=0.5E-02
  ZPKI=ZCKI/PCKI
  PZKI=PCKI/ZCKI
  FEXUC=0.5E-03
  FEXUN=1.0E-02
  FEXUP=1.0E-02
  end subroutine InitUptakePars
end module UptakePars
