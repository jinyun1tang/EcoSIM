module SoilHeatDatatype


  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: TKSZ(366,24,JZ)
  real(r8) :: TKS(0:JZ,JY,JX)
  real(r8) :: THAW(JZ,JY,JX)     !hourly accumulated freeze-thaw flux in micropores
  real(r8) :: HTHAW(JZ,JY,JX)    !hourly accumulated freeze-thaw latent heat flux
  real(r8) :: THAWH(JZ,JY,JX)    !hourly accumulated freeze-thaw flux in macropores
  real(r8) :: XTHAWW(JS,JY,JX)   !hourly accumulated latent heat flux from freeze-thaw
  real(r8) :: TSMX(0:JZ,JY,JX)   !daily maximum soil temperature [oC]
  real(r8) :: TSMN(0:JZ,JY,JX)   !daily minimum soil temperature [oC]
  real(r8) :: VHCP(0:JZ,JY,JX)   !soil heat capacity [MJ m-3 K-1]
  real(r8) :: TCS(0:JZ,JY,JX)    !soil temperature [oC]
  real(r8) :: STC(JZ,JY,JX)      !numerator for soil solid thermal conductivity [MJ m h-1 K-1]
  real(r8) :: DTC(JZ,JY,JX)      !denominator for soil solid thermal conductivity
end module SoilHeatDatatype
