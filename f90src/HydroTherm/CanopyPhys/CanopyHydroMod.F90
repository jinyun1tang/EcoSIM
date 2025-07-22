module CanopyHydroMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use abortutils,        only: endrun
  use DebugToolMod,      only: PrintInfo
  use MiniMathMod,       only: AZMAX1
  use PlantMgmtDataType, only: NP_col
  USE SoilWaterDataType
  use PlantTraitDataType
  use CanopyDataType
  use ClimForcDataType
implicit none

  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  public :: CanopyInterceptPrecp

  real(r8), parameter :: FoliarWatRetcap(0:3)=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)

  contains
!------------------------------------------------------------------------------------------

  subroutine CanopyInterceptPrecp(NY,NX)
  !
  !DESCRIPTION
  !precipitation intercepation by canopy
  implicit none
  integer, intent(in) :: NY,NX
  integer :: NZ
  real(r8) :: CanopyWatHeldCap  !maximum precipitation holding capacity by canopy (leaf+stem) [m3 H2O]
  real(r8) :: prec2canopy_pft   !precipiation onto canopy [m H2O/h]
!
!     CANOPY RETENTION OF PRECIPITATION
!
!     FoliarWatRetcap=foliar surface water retention capacity
!     CanopyLeafArea_pft,CanopyStemArea_pft=leaf,stalk area of PFT
!     FLWC,TFLWC=water retention of PFT,combined canopy
!     PRECA=precipitation+irrigation
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     VOLWC=canopy surface water retention
!
!     Warning: No snofall intercepation is considered at the moment.

  DO  NZ=1,NP_col(NY,NX)
    CanopyWatHeldCap                 = FoliarWatRetcap(iPlantRootProfile_pft(NZ,NY,NX)) &
      *(CanopyLeafArea_pft(NZ,NY,NX)+CanopyStemArea_pft(NZ,NY,NX))
      
    prec2canopy_pft                  = PrecRainAndIrrig_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX)
    PrecIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(prec2canopy_pft,CanopyWatHeldCap-WatHeldOnCanopy_pft(NZ,NY,NX)))
    Prec2Canopy_col(NY,NX)           = Prec2Canopy_col(NY,NX)+prec2canopy_pft
    PrecIntceptByCanopy_col(NY,NX)   = PrecIntceptByCanopy_col(NY,NX)+PrecIntcptByCanopy_pft(NZ,NY,NX)
  ENDDO
  RainPrecThrufall_col(NY,NX) = PrecRainAndIrrig_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX)

  end subroutine CanopyInterceptPrecp

end module CanopyHydroMod