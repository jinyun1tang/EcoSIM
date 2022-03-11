module nitroMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils  , only : endrun
  use minimathmod, only : safe_adb
  use MicrobialDataType
  use NitroPars
  use SOMDataType
  use SoilChemDataType
  use FertilizerDataType
  use NitroDiagTypes
  use nitro1LayerMod
  use NitroDisturbMod
  use GridDataType
  implicit none

  private

  include "blk8a.h"

  character(len=*), parameter :: mod_filename = __FILE__

!
  public :: nitro, initNitro

  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro

  implicit none

  call initNitro1Layer

  end subroutine initNitro

!------------------------------------------------------------------------------------------

  SUBROUTINE nitro(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOIL BIOLOGICAL TRANSFORMATIONS
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: L,NX,NY

!   begin_execution

  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
!
!       VOLWZ=water volume used to calculate aqueous microbial
!       concentrations that drive microbial density effects on
!       decomposition
!
      DO 998 L=0,NL(NY,NX)
        call SoilBGCOneLayer(I,J,L,NY,NX)
998   CONTINUE
!
!       SOC LOSS IF FIRE OR REMOVAL EVENT IS ENTERED IN DISTURBANCE FILE
!
      call SOMRemovalByDisturbance(I,J,NY,NX)
9990  CONTINUE
9995  CONTINUE
  RETURN
  END subroutine nitro

end module nitroMod
