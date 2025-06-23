module PrescribePhenolMod
  use data_kind_mod, only : r8 => DAT_KIND_R8

  use PlantAPIData
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
  public :: SetCanopyProfile
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

!------------------------------------------------------------------------

  subroutine get_LAI_profile(zc,NZ)
  !Description:
  !Distribute leaf area in through the canopy
  !  
  implicit none
  real(r8), intent(in) :: ZC  !plant height
  integer, intent(in) :: NZ   !plant index

  end subroutine get_LAI_profile

!------------------------------------------------------------------------

  subroutine get_StemArea_profile(NZ)
  !
  !Description:
  !Distribute steam area in through the canopy
  !
  implicit none
  integer, intent(in) :: NZ   !plant index


  end subroutine get_StemArea_profile

!------------------------------------------------------------------------

  function TreeStem_diameter_taperEq(h,p_stump_rate,p_stemp_rate,top_ht_rate,h_dbh,dbh,h_cb,h_top)result(d)
 
  ! Estimate diameter at height h using segment-specific taper rates as in Larsen (2017)
  implicit none
 
  real(r8), intent(in) :: h             ! - height at which to estimate diameter
  real(r8), intent(in) :: p_stump_rate  ! - stump taper rate [unitless]
  real(r8), intent(in) :: p_stemp_rate  ! - stem taper rate  [uitless]
  real(r8), intent(in) :: top_ht_rate   ! - crown base height/crown length [unitless]
  real(r8), intent(in) :: h_dbh    ! - breast height [m]
  real(r8), intent(in) :: dbh      ! - diameter at breast height [m]
  real(r8), intent(in) :: h_cb     ! - height to crown base [m]
  real(r8), intent(in) :: h_top    ! - total tree height (real)
  !Output
  real(r8) ::   d                  ! - estimated diameter at height h [m]
  
  if(h<=h_dbh)then
    !below breast height
    d=dbh+p_stump_rate*(h-h_dbh)
  elseif(h<=h_cb)then
    !below crown base 
    d=dbh+p_stemp_rate*(h-h_dbh)  
  else
    !above crown base
    d=(h_top-h)*top_ht_rate
  endif

  end function TreeStem_diameter_taperEq
!------------------------------------------------------------------------

  subroutine SetCanopyProfile(I,J,DepthSurfWatIce,LeafAreaZsec_pft,StemAreaZsec_pft)
  !
  !Description:
  !distribute leaf area and stem area within the crown defined by top and bottom heights.
  !For simplicity, unifrom distribution is adopted
  !
  implicit none
  integer , intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce !water+ice depth at surface
  real(r8),intent(out) :: LeafAreaZsec_pft(NumOfLeafZenithSectors1,NumOfCanopyLayers1,JP1)  !leaf area in different zenith sectors
  real(r8),intent(out) :: StemAreaZsec_pft(NumOfLeafZenithSectors1,NumOfCanopyLayers1,JP1)    !stem area in different zenith sectors

  StemAreaZsec_pft=0._r8
  LeafAreaZsec_pft=0._r8

  end subroutine SetCanopyProfile
end module PrescribePhenolMod
