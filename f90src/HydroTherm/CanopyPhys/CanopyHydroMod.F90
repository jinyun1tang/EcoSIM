module CanopyHydroMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use abortutils,        only: endrun
  use DebugToolMod,      only: PrintInfo
  use MiniMathMod,       only: AZMAX1
  use PlantMgmtDataType, only: NP_col
  USE EcoSIMCtrlDataType,ONLY: ZEROS
  USE SoilWaterDataType
  use PlantTraitDataType
  use CanopyDataType
  use ClimForcDataType
implicit none

  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  public :: CanopyInterceptPrecip

  real(r8), parameter :: FoliarWatRetcap(0:3)=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)                !bryophyte, grasses, deep roots, very deep roots
  real(r8), parameter :: FoliarSnowRetcap(0:4)=real((/1.0E-04,2.0E-04,3.0e-04,2.E-03,7.0E-03/),r8)        !bryophyte, grasses, shrubs, deciduous trees, confiers
  real(r8), parameter :: BulkFactor4Snow(0:4)=real((/0.03E+03,0.02E+03,0.07E+03,0.04E+03,0.1E+03/),r8)    !bulking factor for snow covered canopy radiation calculation [m2/m3 SWE] bryophyte, grasses, shrubs, deciduous trees, confiers
  contains
!------------------------------------------------------------------------------------------

  subroutine CanopyInterceptPrecip(NY,NX)
  !
  !DESCRIPTION
  !precipitation intercepation by canopy
  implicit none
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname='CanopyInterceptPrecip'
  integer :: NZ
  real(r8) :: CanopyWatHeldCap,CanopySnowHeldCap  !maximum rain/snow holding capacity by canopy (leaf+stem) [m3 H2O]~ [ton H2O]
  real(r8) :: rain2canopy_pft,snow2canopy_pft     !rain/snowfall onto canopy [m H2O/h]
  real(r8) :: SLArea                              !stem+leaf surface area [m2]
  real(r8) :: SnowWindUnload,SnowTempUnload,CanopySnowUnload_col,SnowUnload_pft
  real(r8) :: ENGYS
!
!     CANOPY RETENTION OF PRECIPITATION
!
!     FoliarWatRetcap=foliar surface water retention capacity
!     CanopyLeafArea_pft,CanopyStemSurfArea_pft=leaf,stalk area of PFT
!     FLWC,TFLWC=water retention of PFT,combined canopy
!     PRECA=precipitation+irrigation
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     VOLWC=canopy surface water retention
!
!     Warning: No snofall intercepation is considered at the moment.
  call PrintInfo('beg '//subname)
  CanopySnowUnload_col       = 0._r8
  BulkFactor4Snow_col(NY,NX) = 0.0_r8
  fSnowCanopy_col(NY,NX)     = 0._r8
  ENGYS                      = 0._R8
  DO  NZ=1,NP_col(NY,NX)
    if(LeafStalkArea_pft(NZ,NY,NX).GT.0._r8)then
      CanopyWatHeldCap  = FoliarWatRetcap(iPlantRootProfile_pft(NZ,NY,NX))*LeafStalkArea_pft(NZ,NY,NX)
      CanopySnowHeldCap = FoliarSnowRetcap(iPlantSnowIntercepType_pft(NZ,NY,NX))*LeafStalkArea_pft(NZ,NY,NX)

      snow2canopy_pft                  = SnoFalPrec_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX)
      rain2canopy_pft                  = PrecRainAndIrrig_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX)
      SnowIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(snow2canopy_pft,CanopySnowHeldCap-SnowOnCanopy_pft(NZ,NY,NX)))
      RainIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(rain2canopy_pft,CanopyWatHeldCap-WatHeldOnCanopy_pft(NZ,NY,NX)))
      Rain2Canopy_col(NY,NX)           = Rain2Canopy_col(NY,NX)+rain2canopy_pft
      RainIntceptByCanopy_col(NY,NX)   = RainIntceptByCanopy_col(NY,NX)+RainIntcptByCanopy_pft(NZ,NY,NX)
      SnowIntceptByCanopy_col(NY,NX)    = SnowIntceptByCanopy_col(NY,NX)+SnowIntcptByCanopy_pft(NZ,NY,NX)

      !wind-induced unloading
      SnowWindUnload= WindSpeedAtm_col(NY,NX)*SnowOnCanopy_pft(NZ,NY,NX)/1.56E5_r8
      !temperature-induced unloading
      SnowTempUnload = SnowOnCanopy_pft(NZ,NY,NX)*AZMAX1(TKCanopy_pft(NZ,NY,NX)-270._r8)/(1.87E5_r8*3600._r8)

      !update by unloading
      SnowUnload_pft                = AMIN1(SnowWindUnload+SnowTempUnload,SnowOnCanopy_pft(NZ,NY,NX))
      SnowOnCanopy_pft(NZ,NY,NX)    = SnowOnCanopy_pft(NZ,NY,NX)-SnowUnload_pft
      ENGYS                         = ENGYS+SnowUnload_pft*TKCanopy_pft(NZ,NY,NX)
      CanopySnowUnload_col          = CanopySnowUnload_col+SnowUnload_pft
      fSnowCanopy_pft(NZ,NY,NX)     = AZMAX1((SnowOnCanopy_pft(NZ,NY,NX)+SnowIntcptByCanopy_pft(NZ,NY,NX))/CanopySnowHeldCap)**0.15_r8
      BulkFactor4Snow_pft(NZ,NY,NX) = BulkFactor4Snow(iPlantSnowIntercepType_pft(NZ,NY,NX))
      BulkFactor4Snow_col(NY,NX)    = BulkFactor4Snow_col(NY,NX)+BulkFactor4Snow(iPlantSnowIntercepType_pft(NZ,NY,NX))*LeafStalkArea_pft(NZ,NY,NX)
      fSnowCanopy_col(NY,NX)        = fSnowCanopy_col(NY,NX)+fSnowCanopy_pft(NZ,NY,NX)*LeafStalkArea_pft(NZ,NY,NX)
    else
      SnowIntcptByCanopy_pft(NZ,NY,NX) = 0._r8
      RainIntcptByCanopy_pft(NZ,NY,NX) = 0._r8
      fSnowCanopy_pft(NZ,NY,NX)        = 0._r8
      BulkFactor4Snow_pft(NZ,NY,NX)    = 0._R8
    endif
  ENDDO
  RainPrecThrufall_col(NY,NX) = PrecRainAndIrrig_col(NY,NX)-RainIntceptByCanopy_col(NY,NX)
  SnowPrecThrufall_col(NY,NX) = SnoFalPrec_col(NY,NX)-SnowIntceptByCanopy_col(NY,NX)+CanopySnowUnload_col
  IF(ENGYS.GT.ZEROS(NY,NX))THEN
    TKSnowThrufall_col(NY,NX) = (TairK_col(NY,NX)*(SnoFalPrec_col(NY,NX)-SnowIntceptByCanopy_col(NY,NX))+ENGYS)/SnowPrecThrufall_col(NY,NX)
  else
    TKSnowThrufall_col(NY,NX) = TairK_col(NY,NX)
  endif
  if(LeafStalkArea_col(NY,NX).GT.ZEROS(NY,NX))THEN
    BulkFactor4Snow_col(NY,NX) = BulkFactor4Snow_col(NY,NX)/LeafStalkArea_col(NY,NX)
    fSnowCanopy_col(NY,NX)     = fSnowCanopy_col(NY,NX) /LeafStalkArea_col(NY,NX)
  ENDIF  
  call PrintInfo('end '//subname)
  end subroutine CanopyInterceptPrecip

end module CanopyHydroMod
