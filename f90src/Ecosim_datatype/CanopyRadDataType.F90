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
  real(r8),target,allocatable :: RadPAR_zsec(:,:,:,:,:,:)           !direct incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: LeafAngleClass_pft(:,:,:,:)        !fractionction of leaves in different angle classes, [-]
  real(r8),target,allocatable :: LeafAreaZsec_brch(:,:,:,:,:,:,:)   !leaf surface area, [m2 d-2]
  real(r8),target,allocatable :: LeafAUnshaded_zsec(:,:,:,:,:,:,:)  !leaf irradiated surface area, [m2 d-2]
  real(r8),target,allocatable :: StemAreaZsec_brch(:,:,:,:,:,:)     !stem surface area, [m2 d-2]

  real(r8) :: TotSineSkyAngles_grd
  real(r8) :: dangle
  private :: InitAllocate

  contains

  subroutine  InitCanopyRad

  implicit none

  real(r8) :: da
  real(r8) :: aa
  integer :: N
! NumOfLeafZenithSectors: number of leaf inclination groups, 90 deg into NumOfLeafZenithSectors groups

  call InitAllocate

  dangle=PICON2h/real(NumOfLeafZenithSectors,r8)         !the angle section width

  DO N = 1, NumOfLeafZenithSectors
    aa=real(N-0.5,r8)*dangle
    SineLeafAngle(N)=sin(aa)
    CosineLeafAngle(N)=cos(aa)
  ENDDO
  TotSineSkyAngles_grd = 0._r8
  end subroutine InitCanopyRad

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  !JLS: number of leaf sectors divide the whole horizontal circle
  integer :: ncols
  ncols = bounds%ncols

  allocate(SineLeafAngle(NumOfLeafZenithSectors))
  allocate(CosineLeafAngle(NumOfLeafZenithSectors))
  allocate(OMEGA(NumOfSkyAzimuthSects,NumOfLeafZenithSectors,NumOfLeafAzimuthSectors));OMEGA=0._r8
  allocate(OMEGX(NumOfSkyAzimuthSects,NumOfLeafZenithSectors,NumOfLeafAzimuthSectors));OMEGX=0._r8
  allocate(iScatteringDiffus(NumOfSkyAzimuthSects,NumOfLeafZenithSectors,NumOfLeafAzimuthSectors))
  allocate(LeafAngleClass_pft(NumOfLeafZenithSectors,JP,JY,JX));LeafAngleClass_pft=0._r8
  allocate(LeafAreaZsec_brch(NumOfLeafZenithSectors,NumOfCanopyLayers,MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAreaZsec_brch=0._r8
  allocate(LeafAUnshaded_zsec(NumOfLeafZenithSectors,NumOfCanopyLayers,MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAUnshaded_zsec=0._r8
  allocate(RadPAR_zsec(NumOfLeafZenithSectors,NumOfSkyAzimuthSects,NumOfCanopyLayers,JP,JY,JX));RadPAR_zsec=0._r8
  allocate(RadDifPAR_zsec(NumOfLeafZenithSectors,NumOfSkyAzimuthSects,NumOfCanopyLayers,JP,JY,JX));RadDifPAR_zsec=0._r8
  allocate(StemAreaZsec_brch(NumOfLeafZenithSectors,NumOfCanopyLayers,MaxNumBranches,JP,JY,JX));StemAreaZsec_brch=0._r8

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
  call destroy(LeafAUnshaded_zsec)
  call destroy(RadPAR_zsec)
  call destroy(RadDifPAR_zsec)
  call destroy(StemAreaZsec_brch)
  end subroutine DestructCanopyRad
end module CanopyRadDataType
