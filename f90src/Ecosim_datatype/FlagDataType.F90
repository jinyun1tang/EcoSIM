module FlagDataType
  use GridConsts
  use ElmIDMod  
implicit none

  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  integer  ::  ICLM                                !changes to weather data (0=none,1=step,2=transient)
  integer  ::  IMNG                                !flag for land management
  integer  ::  IWTHR                               !weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene

  integer,target,allocatable ::  IYTYP(:,:,:,:)                      !fertilizer release type from fertilizer input file
  integer,target,allocatable ::  iSoilDisturbType_col(:,:,:)                        !soil disturbance type, [-]
  integer,target,allocatable ::  KoppenClimZone_col(:,:)                          !Koppen climate zone
  integer,target,allocatable ::  IFLGV(:,:)                          !flag for irrigation criterion,0=SWC,1=canopy water potential
  integer,target,allocatable ::  iResetSoilProf_col(:,:)                          !disturbance flag
  integer,target,allocatable ::  IFNHB(:,:)                          !banded NH4 fertilizer flag
  integer,target,allocatable ::  IFNOB(:,:)                          !banded NO3 fertilizer flag
  integer,target,allocatable ::  IFPOB(:,:)                          !banded H2PO4 fertilizer flag
  integer,target,allocatable ::  ISOIL(:,:,:,:)                      !flag for calculating FC(1),WP(2),SCNV(3),SCNH(4)
  integer,target,allocatable ::  ISOILR(:,:)                         !natural(0),reconstructed(1) soil profile

  integer,target,allocatable ::  IUTYP(:,:)                          !urea hydrolysis inhibitor type (1=no,2=yes)
  integer,target,allocatable ::  ITILL1(:,:)                         !soil disturbance type, [-]
  integer,target,allocatable ::  IsPlantActive_pft(:,:,:)                        ! flag for living pft
  integer,target,allocatable ::  doInitPlant_pft(:,:,:)                        !PFT initialization flag:0=no,1=yes
  integer,target,allocatable ::  iPlantPhotosynthesisType(:,:,:)                        !plant photosynthetic type (C3 or C4)
  integer,target,allocatable ::  iPlantRootProfile_pft(:,:,:)                        !plant growth type (vascular, non-vascular)
  integer,target,allocatable ::  iPlantPhenolPattern_pft(:,:,:)                        !plant growth habit (annual or perennial)
  integer,target,allocatable ::  iPlantDevelopPattern_pft(:,:,:)                        !plant growth habit (determinate or indeterminate)
  integer,target,allocatable ::  iPlantNfixType_pft(:,:,:)                        !N2 fixation type
  integer,target,allocatable ::  iPlantPhenolType_pft(:,:,:)                        !climate signal for phenological progress none, temperature, water stress)
  integer,target,allocatable ::  iPlantPhotoperiodType_pft(:,:,:)                        !photoperiod type (neutral, long day, short day)
  integer,target,allocatable ::  iPlantTurnoverPattern_pft(:,:,:)                        !phenologically-driven above-ground turnover (all, foliar only, none)
  integer,target,allocatable ::  iPlantGrainType_pft(:,:,:)                        !grain type (below or above-ground), e.g. potato and onion are below
  integer,target,allocatable ::  MY(:,:,:)                           !mycorrhizal type (no or yes)
  integer,target,allocatable ::  IDWaterTable_col(:,:)                   !water table flag from site file
!----------------------------------------------------------------------

contains
  subroutine InitFlagData

  implicit none
  allocate(IYTYP(0:2,366,JY,JX));IYTYP=0
  allocate(iSoilDisturbType_col(366,JY,JX));   iSoilDisturbType_col=0
  allocate(KoppenClimZone_col(JY,JX));       KoppenClimZone_col=0
  allocate(IFLGV(JY,JX));       IFLGV=0
  allocate(iResetSoilProf_col(JY,JX));       iResetSoilProf_col=itrue
  allocate(IFNHB(JY,JX));       IFNHB=0
  allocate(IFNOB(JY,JX));       IFNOB=0
  allocate(IFPOB(JY,JX));       IFPOB=0
  allocate(ISOIL(4,JZ,JY,JX));  ISOIL=isoi_unset     !soil properties unset by default
  allocate(ISOILR(JY,JX));      ISOILR=0
  allocate(IUTYP(JY,JX));       IUTYP=0
  allocate(ITILL1(JY,JX));      ITILL1=0
  allocate(IsPlantActive_pft(JP,JY,JX));    IsPlantActive_pft=iDormant
  allocate(doInitPlant_pft(JP,JY,JX));    doInitPlant_pft=ifalse
  allocate(iPlantPhotosynthesisType(JP,JY,JX));    iPlantPhotosynthesisType=0
  allocate(iPlantRootProfile_pft(JP,JY,JX));    iPlantRootProfile_pft=0
  allocate(iPlantPhenolPattern_pft(JP,JY,JX));    iPlantPhenolPattern_pft=0
  allocate(iPlantDevelopPattern_pft(JP,JY,JX));    iPlantDevelopPattern_pft=0
  allocate(iPlantNfixType_pft(JP,JY,JX));    iPlantNfixType_pft=0
  allocate(iPlantPhenolType_pft(JP,JY,JX));    iPlantPhenolType_pft=0
  allocate(iPlantPhotoperiodType_pft(JP,JY,JX));    iPlantPhotoperiodType_pft=0
  allocate(iPlantTurnoverPattern_pft(JP,JY,JX));    iPlantTurnoverPattern_pft=0
  allocate(iPlantGrainType_pft(JP,JY,JX));    iPlantGrainType_pft=0
  allocate(MY(JP,JY,JX));       MY=0
  allocate(IDWaterTable_col(JY,JX));       IDWaterTable_col=0
  end subroutine InitFlagData

!----------------------------------------------------------------------
  subroutine DestructFlagData
  use abortutils, only : destroy
  implicit none
  call destroy(IYTYP)
  call destroy(iSoilDisturbType_col)
  call destroy(KoppenClimZone_col)
  call destroy(IFLGV)
  call destroy(iResetSoilProf_col)
  call destroy(IFNHB)
  call destroy(IFNOB)
  call destroy(IFPOB)
  call destroy(ISOIL)
  call destroy(ISOILR)
  call destroy(IUTYP)
  call destroy(ITILL1)
  call destroy(IsPlantActive_pft)
  call destroy(doInitPlant_pft)
  call destroy(iPlantPhotosynthesisType)
  call destroy(iPlantRootProfile_pft)
  call destroy(iPlantPhenolPattern_pft)
  call destroy(iPlantDevelopPattern_pft)
  call destroy(iPlantNfixType_pft)
  call destroy(iPlantPhenolType_pft)
  call destroy(iPlantPhotoperiodType_pft)
  call destroy(iPlantTurnoverPattern_pft)
  call destroy(iPlantGrainType_pft)
  call destroy(MY)
  call destroy(IDWaterTable_col)
  end subroutine DestructFlagData

end module FlagDataType
