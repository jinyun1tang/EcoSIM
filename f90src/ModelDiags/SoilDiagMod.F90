module SoilDiagMod
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use CanopyDataType, only: QVegET_col
  use GridDataType,   only: NU, NL
  use EcoSimConst,    only: DENSICE
  use abortutils,     only: endrun
  use SoilBGCDataType
  use GridDataType
  use SurfLitterDataType
  use CanopyDataType
  use SnowDataType
  use BalanceCheckDataType
  use SoilWaterDataType
  use ClimForcDataType
  use SurfSoilDataType
  use EcosimBGCFluxType
  use SoilHeatDataType
  use PlantDataRateType
  use EcoSimSumDataType
implicit none
  public :: DiagSoilGasPressure
  contains
end module SoilDiagMod
