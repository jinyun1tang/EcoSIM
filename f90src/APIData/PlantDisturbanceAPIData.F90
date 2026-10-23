module PlantDisturbanceAPIData
  ! Owns the disturbance API type, instance and lifecycle procedures.
  use PlantAPICommonData
  implicit none
  private
  public :: plt_disturb_init, plt_disturb_destroy

  type, public :: plant_disturb_type

  real(r8) :: XCORP                     !factor for surface litter incorporation and soil mixing,[-]
  real(r8) :: DCORP                     !soil mixing fraction with tillage, [-]
  integer  :: IYTYP                     !fertilizer release type from fertilizer input file,[-]
  integer  :: iSoilDisturbType_col      !soil disturbance type, [-]

  integer,  pointer :: iHarvestDay_pft(:)       => null()  !day of harvest, [-]
  integer,  pointer :: iYearPlanting_pft(:)     => null()  !year of planting,[-]
  integer,  pointer :: iPlantingDay_pft(:)      => null()  !day of planting,[-]
  integer,  pointer :: iPlantingYear_pft(:)     => null()  !year of planting,[-]
  integer,  pointer :: iHarvestYear_pft(:)      => null()  !year of harvest,[-]
  integer,  pointer :: iDayPlanting_pft(:)      => null()  !day of planting,[-]
  integer,  pointer :: iDayPlantHarvest_pft(:)  => null()  !day of harvest,[-]
  integer,  pointer :: iYearPlantHarvest_pft(:) => null()  !year of harvest,[-]
  integer,  pointer :: iHarvstType_pft(:)       => null()  !type of harvest,[-]
  integer,  pointer :: jHarvstType_pft(:)           => null()  !flag for stand replacing disturbance,[-]

  real(r8), pointer :: EcoHavstElmnt_CumYr_col(:)   => null()  !ecosystem harvest element,                                    [gC d-2]
  real(r8), pointer :: O2ByFire_CumYr_pft(:)        => null()  !plant O2 uptake from fire,                                    [g d-2 ]
  real(r8), pointer :: FERT(:)                      => null()  !fertilizer application,                                       [g m-2]
  real(r8), pointer :: FracBiomHarvsted(:,:,:)      => null()  !harvest efficiency,                                           [-]
  real(r8), pointer :: CanopyCutProxy_pft(:)   => null()  !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  real(r8), pointer :: FireLossE_pft(:,:)        => null() !plant element lost by fire, [g d-2 h-1]
  real(r8), pointer :: RootLost2Fire_pft(:,:)     => null() !plant root biomass lost by fire, [g d-2 h-1]
  real(r8), pointer :: THIN_pft(:)                  => null()  !thinning of plant population,                                 [-]
  real(r8), pointer :: CH4ByFire_CumYr_pft(:)       => null()  !plant CH4 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: CO2ByFire_CumYr_pft(:)       => null()  !plant CO2 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: N2ObyFire_CumYr_pft(:)       => null()  !plant N2O emission from fire,                                 [g d-2 ]
  real(r8), pointer :: NH3byFire_CumYr_pft(:)       => null()  !plant NH3 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: PO4byFire_CumYr_pft(:)       => null()  !plant PO4 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: EcoHavstElmntCum_pft(:,:)    => null()  !total plant element harvest,                                  [gC d-2 ]
  real(r8), pointer :: EcoHavstElmnt_CumYr_pft(:,:) => null()  !plant element harvest,                                        [g d-2 ]
  real(r8), pointer :: PlantElmDistLoss_pft(:,:)    => null()  !plant element loss to disturbance,                            [g d-2 h-1]
  contains
    procedure, public :: Init    =>  plt_disturb_init
    procedure, public :: Destroy => plt_disturb_destroy
  end type plant_disturb_type


  type(plant_disturb_type)  , public, target :: plt_distb     !plant disturbance type

contains

  subroutine plt_disturb_init(this)

  implicit none
  class(plant_disturb_type) :: this

  allocate(this%THIN_pft(JP1));this%THIN_pft=spval
  allocate(this%EcoHavstElmnt_CumYr_col(NumPlantChemElms));this%EcoHavstElmnt_CumYr_col=spval
  allocate(this%EcoHavstElmntCum_pft(NumPlantChemElms,JP1));this%EcoHavstElmntCum_pft=spval
  allocate(this%EcoHavstElmnt_CumYr_pft(NumPlantChemElms,JP1));this%EcoHavstElmnt_CumYr_pft=spval
  allocate(this%CH4ByFire_CumYr_pft(JP1));this%CH4ByFire_CumYr_pft=spval
  allocate(this%CO2ByFire_CumYr_pft(JP1));this%CO2ByFire_CumYr_pft=spval
  allocate(this%FireLossE_pft(NumPlantChemElms,JP1)); this%FireLossE_pft=spval
  allocate(this%RootLost2Fire_pft(NumPlantChemElms,JP1)); this%RootLost2Fire_pft=spval
  allocate(this%PlantElmDistLoss_pft(NumPlantChemElms,JP1)); this%PlantElmDistLoss_pft=spval
  allocate(this%N2ObyFire_CumYr_pft(JP1));this%N2ObyFire_CumYr_pft=spval
  allocate(this%NH3byFire_CumYr_pft(JP1));this%NH3byFire_CumYr_pft=spval
  allocate(this%PO4byFire_CumYr_pft(JP1));this%PO4byFire_CumYr_pft=spval

  allocate(this%iHarvestDay_pft(JP1)); this%iHarvestDay_pft=0
  allocate(this%O2ByFire_CumYr_pft(JP1));this%O2ByFire_CumYr_pft=spval
  allocate(this%FracBiomHarvsted(1:2,1:4,JP1));this%FracBiomHarvsted=spval
  allocate(this%CanopyCutProxy_pft(JP1)); this%CanopyCutProxy_pft=spval
  allocate(this%iYearPlantHarvest_pft(JP1));this%iYearPlantHarvest_pft=0
  allocate(this%FERT(1:20));this%FERT=spval
  allocate(this%iYearPlanting_pft(JP1));this%iYearPlanting_pft=0
  allocate(this%iPlantingDay_pft(JP1));this%iPlantingDay_pft=0
  allocate(this%iPlantingYear_pft(JP1));this%iPlantingYear_pft=0
  allocate(this%iHarvestYear_pft(JP1));this%iHarvestYear_pft=0
  allocate(this%iDayPlanting_pft(JP1));this%iDayPlanting_pft=0
  allocate(this%iDayPlantHarvest_pft(JP1));this%iDayPlantHarvest_pft=0
  allocate(this%iHarvstType_pft(JP1));this%iHarvstType_pft=-1
  allocate(this%jHarvstType_pft(JP1));this%jHarvstType_pft=0

  end subroutine plt_disturb_init

  subroutine plt_disturb_destroy(this)
  implicit none
  class(plant_disturb_type) :: this


!  if(allocated(THIN_pft))deallocate(THIN_pft)
!  if(allocated(EcoHavstElmnt_CumYr_pft))deallocate(EcoHavstElmnt_CumYr_pft)
!  if(allocated(EcoHavstElmntCum_pft))deallocate(EcoHavstElmntCum_pft)
!  if(allocated(CH4ByFire_CumYr_pft))deallocate(CH4ByFire_CumYr_pft)
!  if(allocated(CO2ByFire_CumYr_pft))deallocate(CO2ByFire_CumYr_pft)
!  if(allocated(N2ObyFire_CumYr_pft))deallocate(N2ObyFire_CumYr_pft)
!  if(allocated(NH3byFire_CumYr_pft))deallocate(NH3byFire_CumYr_pft)
!  if(allocated(PO4byFire_CumYr_pft))deallocate(PO4byFire_CumYr_pft)

!  if(allocated(iHarvestDay_pft))deallocate(iHarvestDay_pft)
!  if(allocated(O2ByFire_CumYr_pft))deallocate(O2ByFire_CumYr_pft)
!  if(allocated(FracBiomHarvsted))deallocate(FracBiomHarvsted)
!  if(allocated(HVST))deallocate(HVST)
!  if(allocated(iYearPlantHarvest_pft))deallocate(iYearPlantHarvest_pft)
!  if(allocated(iDayPlanting_pft))deallocate(iDayPlanting_pft)
!  if(allocated(iDayPlantHarvest_pft))deallocate(iDayPlantHarvest_pft)
!  if(allocated(iHarvestYear_pft))deallocate(iHarvestYear_pft)
!  if(allocated(iPlantingDay_pft))deallocate(iPlantingDay_pft)
!  if(allocated(iYearPlanting_pft))deallocate(iYearPlanting_pft)
!  if(allocated(FERT))deallocate(FERT)
!  if(allocated(iPlantingYear_pft))deallocate(iPlantingYear_pft)
!  if(allocated(iHarvstType_pft))deallocate(iHarvstType_pft)
!  if(allocated(jHarvstType_pft))deallocate(jHarvstType_pft)

  end subroutine plt_disturb_destroy
end module PlantDisturbanceAPIData
