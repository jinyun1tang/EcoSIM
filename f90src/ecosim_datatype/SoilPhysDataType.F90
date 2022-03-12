module SoilPhysDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
implicit none

  save
  real(r8) :: CORGCI(JZ,JY,JX)
  real(r8) :: POROSI(JZ,JY,JX)
  real(r8) :: FHOLI(JZ,JY,JX)
  real(r8) :: SLOPE(0:3,JY,JX)
  real(r8) :: CSAND(JZ,JY,JX)   !soil sand content [kg Mg-1]
  real(r8) :: CSILT(JZ,JY,JX)   !soil silt content [kg Mg-1]
  real(r8) :: CCLAY(JZ,JY,JX)   !soil clay content [kg Mg-1]
  REAL(R8) :: ALBS(JY,JX)       !snowpack albedo
  REAL(R8) :: ROCK(JZ,JY,JX)    !Rock fraction
  real(r8) :: TSED(JY,JX)       !erosion rate [Mg d-2 t-1]
  real(r8) :: BKDS(0:JZ,JY,JX)  !soil bulk density, [Mg m-3]
  real(r8) :: FC(0:JZ,JY,JX)    !water contents at field capacity
  real(r8) :: WP(0:JZ,JY,JX)    !water contents at wilting point
  real(r8) :: SCNV(0:JZ,JY,JX)  !soil vertical saturated hydraulic conductivity [mm h-1]
  real(r8) :: SCNH(JZ,JY,JX)    !soil horizontal saturated hydraulic conductivity, [mm h-1]
  real(r8) :: BKDSI(JZ,JY,JX)   !initial bulk density [Mg m-3,0=water]
  real(r8) :: FMPR(0:JZ,JY,JX)  !micropore fraction
  real(r8) :: FHOL(JZ,JY,JX)    !macropore fraction
  real(r8) :: PHOL(JZ,JY,JX)    !path length between macopores
  real(r8) :: HRAD(JZ,JY,JX)    !radius of macropores
  real(r8) :: PSIFC(JY,JX)      !water potentials at field capacity, [MPa]
  real(r8) :: PSIWP(JY,JX)      !water potentials at wilting point [MPa]
  real(r8) :: THW(JZ,JY,JX)     !initial soil water content
  real(r8) :: THI(JZ,JY,JX)     !initial ice content
  integer :: NHOL(JZ,JY,JX)     !number of macropores
  REAL(R8) :: ALBX(JY,JX)       !Surface albedo
  real(r8) :: POROS(0:JZ,JY,JX) !
  real(r8) :: PSL(0:JZ,JY,JX)   !
  real(r8) :: FCL(0:JZ,JY,JX)   !
  real(r8) :: WPL(0:JZ,JY,JX)   !
  real(r8) :: PSD(0:JZ,JY,JX)   !
  real(r8) :: FCD(0:JZ,JY,JX)   !
  real(r8) :: VOLX(0:JZ,JY,JX)  !
  real(r8) :: VOLY(0:JZ,JY,JX)  !
  real(r8) :: BKVL(0:JZ,JY,JX)  !
  real(r8) :: SRP(0:JZ,JY,JX)   !
  real(r8) :: POROS0(JY,JX)     !
  real(r8) :: FSLOPE(2,JY,JX)
  real(r8) :: BKVLNM(JY,JX)
  real(r8) :: BKVLNU(JY,JX)
  REAL(R8) :: VOLAI(0:JZ,JY,JX)
  REAL(R8) :: PSIMS(JY,JX)
  REAL(R8) :: PSIMX(JY,JX)
  REAL(R8) :: PSIMN(JY,JX)
  REAL(R8) :: PSISD(JY,JX)
  REAL(R8) :: PSIMD(JY,JX)
  real(r8) :: SAND(JZ,JY,JX)
  real(r8) :: SILT(JZ,JY,JX)
  real(r8) :: CLAY(JZ,JY,JX) 
end module SoilPhysDataType
