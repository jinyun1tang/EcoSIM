module SoilPhysDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: SLOPE(0:3,JY,JX)                 !slope	in three directions [o]
  real(r8) :: FC(0:JZ,JY,JX)                   !water contents at field capacity
  real(r8) :: WP(0:JZ,JY,JX)                   !water contents at wilting point
  real(r8) :: SCNV(0:JZ,JY,JX)                 !soil vertical saturated hydraulic conductivity [mm h-1]
  real(r8) :: SCNH(JZ,JY,JX)                   !soil horizontal saturated hydraulic conductivity, [mm h-1]

  real(r8) :: PSIFC(JY,JX)                     !water potentials at field capacity, [MPa]
  real(r8) :: PSIWP(JY,JX)                     !water potentials at wilting point [MPa]
  real(r8) :: THW(JZ,JY,JX)                    !initial soil water content
  real(r8) :: THI(JZ,JY,JX)                    !initial ice content
  REAL(R8) :: ALBX(JY,JX)                      !Surface albedo
  real(r8) :: PSL(0:JZ,JY,JX)                  !log soil porosity	-
  real(r8) :: FCL(0:JZ,JY,JX)                  !log water content at field capacity
  real(r8) :: WPL(0:JZ,JY,JX)                  !log water content at wilting point
  real(r8) :: PSD(0:JZ,JY,JX)                  !log soil porosity - log water content at field capacity
  real(r8) :: FCD(0:JZ,JY,JX)                  !log water content at field capacity

  real(r8) :: SRP(0:JZ,JY,JX)                  !shape parameter for water desorption
  real(r8) :: FSLOPE(2,JY,JX)                  !fraction of slope in 1 and 2
  REAL(R8) :: VOLAI(0:JZ,JY,JX)                !initial total soil micropore porosity	m3 d-2
  REAL(R8) :: PSIMS(JY,JX)                     !log water potential at saturation	MPa
  REAL(R8) :: PSIMX(JY,JX)                     !log water potential at field capacity	-
  REAL(R8) :: PSIMN(JY,JX)                     !log water potential at wilting point
  REAL(R8) :: PSISD(JY,JX)                     !log water potential at field capacity 	-
  REAL(R8) :: PSIMD(JY,JX)                     !log water potential at saturation - log water potential at field capacity
  real(r8) :: VHCM(0:JZ,JY,JX)                 !soil solid heat capacity [MPa m-3 K-1]
  real(r8) :: DPTHA(JY,JX)                     !active layer depth, [n]


end module SoilPhysDataType
