module PrescribePhenolMod
  use data_kind_mod, only : r8 => DAT_KIND_R8

implicit none

  !pft information from CLM4.5
  !pft                                             ra             rb           ztop(m)      zbot(m)
  !Needleleaf evergreen tree/NET Temperate         7.0            2.0          17           8.5
  !Needleleaf evergreen tree/NET Boreal            7.0            2.0          17           8.5
  !Needleleaf deciduous tree/NDT Boreal            7.0            2.0          14           7
  !Broadleave evergreen tree/BET tropical          7.0            1.0          35           1
  !Broadleave evergreen tree/BET temperate         7.0            1.0          35           1
  !Broadleave deciduous tree/BDT tropical          6.0            2.0          18           10
  !Broadleave deciduous tree/BDT temperate         6.0            2.0          20           11.5
  !Broadleave deciduous tree/BDT Boreal            6.0            2.0          20           11.5
  !Broadleave evergreen shrub/BES temperate        7.0            1.5          0.5          0.1
  !Broadleave deciduous shrub/BDS temperate        7.0            1.5          0.5          0.1
  !Broadleave deciduous shrub/BDS Boreal           7.0            1.5          0.5          0.1
  !C3 Arctic grass                                11.0            2.0          0.5          0.01
  !C3 Grass                                       11.0            2.0          0.5          0.01
  !C4 grass                                       11.0            2.0          0.5          0.01
  !C3 Crop Rainfed                                6.0             3.0          0.5          0.01
  !C3 Crop irrigated                              6.0             3.0          0.5          0.01
  !Corn Rainfed                                   6.0             3.0          0.5          0.01
  !Corn irrigated                                 6.0             3.0          0.5          0.01
  !Rainfed temperate cereals                      6.0             3.0          0.5          0.01
  !irrigated temperate cereals                    6.0             3.0          0.5          0.01
  !Rainfed winter cereals                         6.0             3.0          0.5          0.01
  !irrigated winter cereals                       6.0             3.0          0.5          0.01
  !Rainfed soybean                                6.0             3.0          0.5          0.01
  !irrigated soybean                              6.0             3.0          0.5          0.01
  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: GetRootProfile
  contains

  subroutine GetRootProfile(ra,rb,NL,ZH,RootProfl)
  implicit none
  real(r8), intent(in) :: ra
  real(r8), intent(in) :: rb
  integer,  intent(in) :: NL
  real(r8), intent(in) :: Zh(0:NL)
  real(r8), intent(out) :: RootProfl(1:NL)
  real(r8) :: rz(0:NL-1)
  integer :: L

  DO L=0,NL
    rz(L)=exp(-ra*zh(L))+exp(-rb*zh(L))
  enddo
  DO L=1,NL-1
    rootProfl(L)=0.5_r8*(rz(L-1)-rz(L))
  ENDDO
  rootProfl(NL)=0.5_r8*rz(NL-1)

  end subroutine GetRootProfile

end module PrescribePhenolMod
