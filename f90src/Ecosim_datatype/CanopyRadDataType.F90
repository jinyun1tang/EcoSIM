module CanopyRadDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GridConsts
  implicit none
  public
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable :: SineLeafAngle(:)                   !sine of leaf angle,	[-]
  real(r8),target,allocatable :: CosineLeafAngle(:)                 !cosine of leaf angle,	[-]
  real(r8),target,allocatable :: OMEGA(:,:,:)                       !sine of indirect sky radiation on leaf surface, [-]
  real(r8),target,allocatable :: OMEGX(:,:,:)                       !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation, [-]
  integer,target,allocatable :: iScatteringDiffus(:,:,:)            !flag for calculating backscattering of radiation in canopy, [-]
  real(r8),target,allocatable :: RadDifPAR_zsec(:,:,:,:,:,:)        !diffuse incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: RadTotPAR_zsec(:,:,:,:,:,:)           !direct incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: LeafAngleClass_pft(:,:,:,:)        !fractionction of leaves in different angle classes, [-]
  real(r8),target,allocatable :: LeafAreaZsec_brch(:,:,:,:,:,:,:)   !leaf surface area, [m2 d-2]
  real(r8),target,allocatable :: LeafAreaSunlit_zsec(:,:,:,:,:,:,:)  !leaf irradiated surface area, [m2 d-2]
  real(r8),target,allocatable :: StemAreaZsec_brch(:,:,:,:,:,:)     !stem surface area, [m2 d-2]
  real(r8),target,allocatable :: RadSW_Canopy_col(:,:)              !canopy intercepted shortwave radiation, [MJ d-2]
  real(r8) :: TotSineSkyAngles_grd
  real(r8) :: dangle
  private :: InitAllocate

  contains

  subroutine  InitCanopyRad

  implicit none

  real(r8) :: da
  real(r8) :: aa
  integer :: N
! NumLeafZenithSectors: number of leaf inclination groups, 90 deg into NumLeafZenithSectors groups

  call InitAllocate

  dangle=PICON2h/real(NumLeafZenithSectors,r8)         !the angle section width

  DO N = 1, NumLeafZenithSectors
    aa                 = real(N-0.5,r8)*dangle
    SineLeafAngle(N)   = sin(aa)
    CosineLeafAngle(N) = cos(aa)
  ENDDO
  TotSineSkyAngles_grd = 0._r8
  end subroutine InitCanopyRad

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  !JLS: number of leaf sectors divide the whole horizontal circle
  integer :: ncols
  ncols = bounds%ncols

  allocate(SineLeafAngle(NumLeafZenithSectors))
  allocate(CosineLeafAngle(NumLeafZenithSectors))
  allocate(OMEGA(NumOfSkyAzimuthSects,NumLeafZenithSectors,NumOfLeafAzimuthSectors));OMEGA=0._r8
  allocate(OMEGX(NumOfSkyAzimuthSects,NumLeafZenithSectors,NumOfLeafAzimuthSectors));OMEGX=0._r8
  allocate(iScatteringDiffus(NumOfSkyAzimuthSects,NumLeafZenithSectors,NumOfLeafAzimuthSectors))
  allocate(LeafAngleClass_pft(NumLeafZenithSectors,JP,JY,JX));LeafAngleClass_pft=0._r8
  allocate(LeafAreaZsec_brch(NumLeafZenithSectors,NumCanopyLayers,MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAreaZsec_brch=0._r8
  allocate(LeafAreaSunlit_zsec(NumLeafZenithSectors,NumCanopyLayers,MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAreaSunlit_zsec=0._r8
  allocate(RadTotPAR_zsec(NumLeafZenithSectors,NumOfSkyAzimuthSects,NumCanopyLayers,JP,JY,JX));RadTotPAR_zsec=0._r8
  allocate(RadDifPAR_zsec(NumLeafZenithSectors,NumOfSkyAzimuthSects,NumCanopyLayers,JP,JY,JX));RadDifPAR_zsec=0._r8
  allocate(StemAreaZsec_brch(NumLeafZenithSectors,NumCanopyLayers,MaxNumBranches,JP,JY,JX));StemAreaZsec_brch=0._r8
  allocate(RadSW_Canopy_col(JY,JX)); RadSW_Canopy_col=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructCanopyRad
  use abortutils, only : destroy
  use abortutils, only : destroy
  implicit none
  call destroy(SineLeafAngle)
  call destroy(CosineLeafAngle)
  call destroy(OMEGA)
  call destroy(OMEGX)
  call destroy(iScatteringDiffus)
  call destroy(LeafAngleClass_pft)
  call destroy(LeafAreaZsec_brch)
  call destroy(LeafAreaSunlit_zsec)
  call destroy(RadTotPAR_zsec)
  call destroy(RadDifPAR_zsec)
  call destroy(StemAreaZsec_brch)
  call destroy(RadSW_Canopy_col)
  end subroutine DestructCanopyRad
end module CanopyRadDataType
