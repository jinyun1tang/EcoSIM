module FlagDataType
  use GridConsts
  use ElmIDMod  
implicit none

  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  integer  ::  ISALTG                              !salt option
  integer  ::  IERSNG                              !erosion option
  integer  ::  ICLM                                !changes to weather data (0=none,1=step,2=transient)
  integer  ::  IMNG                                !flag for land management
  integer  ::  IWTHR(2)                            !weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene

  integer,target,allocatable ::  IYTYP(:,:,:,:)                      !fertilizer release type from fertilizer input file
  integer,target,allocatable ::  ITILL(:,:,:)                        !soil disturbance type, [-]
  integer,target,allocatable ::  IETYP(:,:)                          !Koppen climate zone
  integer,target,allocatable ::  IFLGV(:,:)                          !flag for irrigation criterion,0=SWC,1=canopy water potential
  integer,target,allocatable ::  IFLGS(:,:)                          !disturbance flag
  integer,target,allocatable ::  IFNHB(:,:)                          !banded NH4 fertilizer flag
  integer,target,allocatable ::  IFNOB(:,:)                          !banded NO3 fertilizer flag
  integer,target,allocatable ::  IFPOB(:,:)                          !banded H2PO4 fertilizer flag
  integer,target,allocatable ::  ISOIL(:,:,:,:)                      !flag for calculating FC(1),WP(2),SCNV(3),SCNH(4)
  integer,target,allocatable ::  ISOILR(:,:)                         !natural(0),reconstructed(1) soil profile

  integer,target,allocatable ::  IUTYP(:,:)                          !urea hydrolysis inhibitor type (1=no,2=yes)
  integer,target,allocatable ::  ITILL1(:,:)                         !soil disturbance type, [-]
  integer,target,allocatable ::  IsPlantActive(:,:,:)                        ! flag for living pft
  integer,target,allocatable ::  doInitPlant(:,:,:)                        !PFT initialization flag:0=no,1=yes
  integer,target,allocatable ::  iPlantPhotosynthesisType(:,:,:)                        !plant photosynthetic type (C3 or C4)
  integer,target,allocatable ::  iPlantMorphologyType(:,:,:)                        !plant growth type (vascular, non-vascular)
  integer,target,allocatable ::  iPlantPhenologyPattern(:,:,:)                        !plant growth habit (annual or perennial)
  integer,target,allocatable ::  iPlantDevelopPattern(:,:,:)                        !plant growth habit (determinate or indeterminate)
  integer,target,allocatable ::  iPlantNfixType(:,:,:)                        !N2 fixation type
  integer,target,allocatable ::  iPlantPhenologyType(:,:,:)                        !climate signal for phenological progress none, temperature, water stress)
  integer,target,allocatable ::  iPlantPhotoperiodType(:,:,:)                        !photoperiod type (neutral, long day, short day)
  integer,target,allocatable ::  iPlantTurnoverPattern(:,:,:)                        !phenologically-driven above-ground turnover (all, foliar only, none)
  integer,target,allocatable ::  iPlantGrainType(:,:,:)                        !grain type (below or above-ground), e.g. potato and onion are below
  integer,target,allocatable ::  MY(:,:,:)                           !mycorrhizal type (no or yes)
  integer,target,allocatable ::  IDWaterTable(:,:)                   !water table flag from site file
!----------------------------------------------------------------------

contains
  subroutine InitFlagData

  implicit none
  allocate(IYTYP(0:2,366,JY,JX));IYTYP=0
  allocate(ITILL(366,JY,JX));   ITILL=0
  allocate(IETYP(JY,JX));       IETYP=0
  allocate(IFLGV(JY,JX));       IFLGV=0
  allocate(IFLGS(JY,JX));       IFLGS=0
  allocate(IFNHB(JY,JX));       IFNHB=0
  allocate(IFNOB(JY,JX));       IFNOB=0
  allocate(IFPOB(JY,JX));       IFPOB=0
  allocate(ISOIL(4,JZ,JY,JX));  ISOIL=0
  allocate(ISOILR(JY,JX));      ISOILR=0
  allocate(IUTYP(JY,JX));       IUTYP=0
  allocate(ITILL1(JY,JX));      ITILL1=0
  allocate(IsPlantActive(JP,JY,JX));    IsPlantActive=0
  allocate(doInitPlant(JP,JY,JX));    doInitPlant=ifalse
  allocate(iPlantPhotosynthesisType(JP,JY,JX));    iPlantPhotosynthesisType=0
  allocate(iPlantMorphologyType(JP,JY,JX));    iPlantMorphologyType=0
  allocate(iPlantPhenologyPattern(JP,JY,JX));    iPlantPhenologyPattern=0
  allocate(iPlantDevelopPattern(JP,JY,JX));    iPlantDevelopPattern=0
  allocate(iPlantNfixType(JP,JY,JX));    iPlantNfixType=0
  allocate(iPlantPhenologyType(JP,JY,JX));    iPlantPhenologyType=0
  allocate(iPlantPhotoperiodType(JP,JY,JX));    iPlantPhotoperiodType=0
  allocate(iPlantTurnoverPattern(JP,JY,JX));    iPlantTurnoverPattern=0
  allocate(iPlantGrainType(JP,JY,JX));    iPlantGrainType=0
  allocate(MY(JP,JY,JX));       MY=0
  allocate(IDWaterTable(JY,JX));       IDWaterTable=0
  end subroutine InitFlagData

!----------------------------------------------------------------------
  subroutine DestructFlagData
  use abortutils, only : destroy
  implicit none
  call destroy(IYTYP)
  call destroy(ITILL)
  call destroy(IETYP)
  call destroy(IFLGV)
  call destroy(IFLGS)
  call destroy(IFNHB)
  call destroy(IFNOB)
  call destroy(IFPOB)
  call destroy(ISOIL)
  call destroy(ISOILR)
  call destroy(IUTYP)
  call destroy(ITILL1)
  call destroy(IsPlantActive)
  call destroy(doInitPlant)
  call destroy(iPlantPhotosynthesisType)
  call destroy(iPlantMorphologyType)
  call destroy(iPlantPhenologyPattern)
  call destroy(iPlantDevelopPattern)
  call destroy(iPlantNfixType)
  call destroy(iPlantPhenologyType)
  call destroy(iPlantPhotoperiodType)
  call destroy(iPlantTurnoverPattern)
  call destroy(iPlantGrainType)
  call destroy(MY)
  call destroy(IDWaterTable)
  end subroutine DestructFlagData

end module FlagDataType
