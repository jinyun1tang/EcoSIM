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

  integer,target,allocatable ::  IYTYP(:,:,:,:)                      !fertilizer release type from fertilizer input file,[-]
  integer,target,allocatable ::  iSoilDisturbType_col(:,:,:)         !soil disturbance type, [-]
  integer,target,allocatable ::  KoppenClimZone_col(:,:)             !Koppen climate zone,[-]
  integer,target,allocatable ::  IFLGV_col(:,:)                      !flag for irrigation criterion,0=SWC,1=canopy water potential,[-]
  integer,target,allocatable ::  iResetSoilProf_col(:,:)             !disturbance flag,[-]
  integer,target,allocatable ::  IFNHB_col(:,:)                      !banded NH4 fertilizer flag,[-]
  integer,target,allocatable ::  IFNOB_col(:,:)                      !banded NO3 fertilizer flag,[-]
  integer,target,allocatable ::  IFPOB_col(:,:)                      !banded H2PO4 fertilizer flag,[-]
  integer,target,allocatable ::  ISOIL_vr(:,:,:,:)                   !flag for calculating FC(1),WP(2),SCNV(3),SCNH(4),[-]
  integer,target,allocatable ::  ISOILR_col(:,:)                     !natural(0),reconstructed(1) soil profile,[-]

  integer,target,allocatable ::  iUreaHydInhibitorType_col(:,:)        !urea hydrolysis inhibitor type (1=no,2=yes),[-]
  integer,target,allocatable ::  ITILL1_col(:,:)                       !soil disturbance type, [-]
  integer,target,allocatable ::  IsPlantActive_pft(:,:,:)              ! flag for living pft,[-]
  integer,target,allocatable ::  doInitPlant_pft(:,:,:)                !PFT initialization flag:0=no,1=yes,[-]
  integer,target,allocatable ::  IDWaterTable_col(:,:)                 !water table flag from site file,[-]
!----------------------------------------------------------------------

contains
  subroutine InitFlagData

  implicit none
  allocate(IYTYP(0:2,366,JY,JX));IYTYP=0
  allocate(iSoilDisturbType_col(366,JY,JX));   iSoilDisturbType_col=0
  allocate(KoppenClimZone_col(JY,JX));       KoppenClimZone_col=0
  allocate(IFLGV_col(JY,JX));       IFLGV_col=0
  allocate(iResetSoilProf_col(JY,JX));       iResetSoilProf_col=itrue
  allocate(IFNHB_col(JY,JX));       IFNHB_col=0
  allocate(IFNOB_col(JY,JX));       IFNOB_col=0
  allocate(IFPOB_col(JY,JX));       IFPOB_col=0
  allocate(ISOIL_vr(4,JZ,JY,JX));  ISOIL_vr=isoi_unset     !soil properties unset by default
  allocate(ISOILR_col(JY,JX));      ISOILR_col=0
  allocate(iUreaHydInhibitorType_col(JY,JX));       iUreaHydInhibitorType_col=0
  allocate(ITILL1_col(JY,JX));      ITILL1_col=0
  allocate(IsPlantActive_pft(JP,JY,JX));    IsPlantActive_pft=iDormant
  allocate(doInitPlant_pft(JP,JY,JX));    doInitPlant_pft=ifalse
  allocate(IDWaterTable_col(JY,JX));       IDWaterTable_col=0
  end subroutine InitFlagData

!----------------------------------------------------------------------
  subroutine DestructFlagData
  use abortutils, only : destroy
  implicit none
  call destroy(IYTYP)
  call destroy(iSoilDisturbType_col)
  call destroy(KoppenClimZone_col)
  call destroy(IFLGV_col)
  call destroy(iResetSoilProf_col)
  call destroy(IFNHB_col)
  call destroy(IFNOB_col)
  call destroy(IFPOB_col)
  call destroy(ISOIL_vr)
  call destroy(ISOILR_col)
  call destroy(iUreaHydInhibitorType_col)
  call destroy(ITILL1_col)
  call destroy(IsPlantActive_pft)
  call destroy(doInitPlant_pft)
  call destroy(IDWaterTable_col)
  end subroutine DestructFlagData

end module FlagDataType
