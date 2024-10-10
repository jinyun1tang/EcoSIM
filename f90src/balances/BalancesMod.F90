module BalancesMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use CanopyDataType, only : QvET_col  
  use BalanceCheckDataType
  use SoilWaterDataType
  use ClimForcDataType
  use SurfSoilDataType
  use EcosimBGCFluxType
  use SoilHeatDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: BegCheckBalances
  public :: EndCheckBalances
contains

  subroutine BegCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX
  
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      WaterErr_col(NY,NX) = WatMass_col(NY,NX)
      HeatErr_col(NY,NX)  = HeatStore_col(NY,NX)
    ENDDO
  ENDDO  
  end subroutine BegCheckBalances 

!------------------------------------------------------------------------------------------

  subroutine EndCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      WaterErr_col(NY,NX) = WaterErr_col(NY,NX)-WatMass_col(NY,NX)+(RainFalPrec(NY,NX)+SnoFalPrec_col(NY,NX)) &
        +VapXAir2GSurf_col(NY,NX)+QvET_col(NY,NX)-QDischar_col(NY,NX)+QRunSurf_col(NY,NX)-QDrain_col(NY,NX)
      HeatErr_col(NY,NX) = HeatErr_col(NY,NX)- HeatStore_col(NY,NX)+Eco_NetRad_col(NY,NX)+Eco_Heat_Latent_col(NY,NX) &
        +Eco_Heat_Sens_col(NY,NX)+HeatRunSurf_col(NY,NX)-HeatDrain_col(NY,NX)

!      write(114,*)I+J/24., WaterErr_col(NY,NX),WatMass_col(NY,NX),(RainFalPrec(NY,NX)+SnoFalPrec_col(NY,NX)), &
!        VapXAir2GSurf_col(NY,NX)+QvET_col(NY,NX),-QDischar_col(NY,NX)-QDrain_col(NY,NX)+QRunSurf_col(NY,NX)
    ENDDO
  ENDDO
  end subroutine EndCheckBalances 

end module BalancesMod
