module GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: TAREA               !total area of landscape	[m2]
  real(r8),target,allocatable ::  CumDepz2LayBottom_vr(:,:,:)                       !depth to bottom of soil layer [m]
  real(r8),target,allocatable ::  DLYR_3D(:,:,:,:)                      !thickness of soil layer [m]
  real(r8),target,allocatable ::  DLYRI_3D(:,:,:,:)                     !thickness of soil layer in 3 directions [m]
  real(r8),target,allocatable ::  XDPTH(:,:,:,:)                     !cross-sectional area / distance between adjacent grid cells [m]
  real(r8),target,allocatable ::  SoilDepthMidLay_vr(:,:,:)                        !depth to middle of soil layer [m]
  real(r8),target,allocatable ::  CumSoilThickness_vr(:,:,:)                      !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8),target,allocatable ::  CumSoilThickMidL_vr(:,:,:)                       !depth to middle of soil layer from  surface of grid cell [m]
  real(r8),target,allocatable ::  AREA(:,:,:,:)                      !cross-sectional area  [m2 d-2]
  real(r8),target,allocatable ::  DIST(:,:,:,:)                      !distance between adjacent layers:1=EW,2=NS,3=vertical [m]
  integer,target,allocatable ::  NU(:,:)                             !soil surface layer number
  integer,target,allocatable ::  NUI(:,:)                            !initial soil surface layer number
  integer,target,allocatable ::  MaxNumRootLays(:,:)                             !maximum root layer number
  integer,target,allocatable ::  NK(:,:)                             !additional soil lower boundary layers
  integer,target,allocatable ::  NLI(:,:)                            !initial lowest soil layer number
  integer,target,allocatable ::  NL(:,:)                             !lowest soil layer number
  integer,target,allocatable ::  NUM(:,:)                            !new surface layer number
  real(r8),target,allocatable ::  CumLitRDepzInit_col(:,:)               !initial position of the bottom of liter layer [m]
  REAL(R8),target,allocatable ::  ALAT(:,:)                          !latitude	[degrees]
  real(r8),target,allocatable ::  DH(:,:)                            !EW width of the grid cells, [m]
  real(r8),target,allocatable ::  DV(:,:)                            !NS width of the grid cells, [m]
  integer,target,allocatable ::  FlowDirIndicator(:,:)               !dimension of low
  integer,target,allocatable ::  LSG(:,:,:)                          !match PFT from different scenarios
!----------------------------------------------------------------------

contains

  subroutine InitGridData

  implicit none
  allocate(CumDepz2LayBottom_vr(0:JZ,JY,JX));  CumDepz2LayBottom_vr=0._r8
  allocate(DLYR_3D(3,0:JZ,JY,JX)); DLYR_3D=0._r8
  allocate(DLYRI_3D(3,0:JZ,JY,JX));DLYRI_3D=0._r8
  allocate(XDPTH(3,JZ,JY,JX));  XDPTH=0._r8
  allocate(SoilDepthMidLay_vr(JZ,JY,JX));     SoilDepthMidLay_vr=0._r8
  allocate(CumSoilThickness_vr(0:JZ,JY,JX)); CumSoilThickness_vr=0._r8
  allocate(CumSoilThickMidL_vr(JZ,JY,JX));    CumSoilThickMidL_vr=0._r8
  allocate(AREA(3,0:JZ,JY,JX)); AREA=0._r8
  allocate(DIST(3,JD,JV,JH));   DIST=0._r8
  allocate(NU(JY,JX));          NU=0
  allocate(NUI(JY,JX));         NUI=0
  allocate(MaxNumRootLays(JY,JX));          MaxNumRootLays=0
  allocate(NK(JY,JX));          NK=0
  allocate(NLI(JV,JH));         NLI=0
  allocate(NL(JV,JH));          NL=0
  allocate(NUM(JY,JX));         NUM=0
  allocate(CumLitRDepzInit_col(JY,JX));      CumLitRDepzInit_col=0._r8
  allocate(ALAT(JY,JX));        ALAT=0._r8
  allocate(DH(JY,JX));          DH=0._r8
  allocate(DV(JY,JX));          DV=0._r8
  allocate(FlowDirIndicator(JY,JX));         FlowDirIndicator=3   !vertical by default
  allocate(LSG(JZ,JY,JX));      LSG=0

  end subroutine InitGridData

!----------------------------------------------------------------------
  subroutine DestructGridData
  use abortutils, only : destroy
  implicit none

  call destroy(CumDepz2LayBottom_vr)
  call destroy(DLYR_3D)
  call destroy(DLYRI_3D)
  call destroy(XDPTH)
  call destroy(SoilDepthMidLay_vr)
  call destroy(CumSoilThickness_vr)
  call destroy(CumSoilThickMidL_vr)
  call destroy(AREA)
  call destroy(DIST)
  call destroy(NU)
  call destroy(NUI)
  call destroy(MaxNumRootLays)
  call destroy(NK)
  call destroy(NLI)
  call destroy(NL)
  call destroy(NUM)
  call destroy(CumLitRDepzInit_col)
  call destroy(ALAT)
  call destroy(DH)
  call destroy(DV)
  call destroy(FlowDirIndicator)
  call destroy(LSG)
  end subroutine DestructGridData

end module GridDataType
