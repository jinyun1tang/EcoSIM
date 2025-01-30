module SoilDiagsMod
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
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: DiagSoilGasPressure
  contains

  subroutine DiagSoilGasPressure(I,J,NHW,NHE,NVN,NVS)  

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: idg, NY,NX, L


  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)    
        do idg=idg_beg,idg_end

        enddo
      ENDDO
    ENDDO
  ENDDO    
  end subroutine DiagSoilGasPressure
end module SoilDiagsMod
