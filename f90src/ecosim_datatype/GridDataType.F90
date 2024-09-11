module GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: TAREA               !total area of landscape	[m2]
  real(r8),target,allocatable ::  CumDepz2LayerBot_vr(:,:,:)                       !depth to bottom of soil layer [m]
  real(r8),target,allocatable ::  DLYR(:,:,:,:)                      !thickness of soil layer [m]
  real(r8),target,allocatable ::  DLYRI_3D(:,:,:,:)                     !thickness of soil layer in 3 directions [m]
  real(r8),target,allocatable ::  XDPTH(:,:,:,:)                     !cross-sectional area / distance between adjacent grid cells [m]
  real(r8),target,allocatable ::  SoiDepthMidLay_vr(:,:,:)                        !depth to middle of soil layer [m]
  real(r8),target,allocatable ::  CumSoilThickness_vr(:,:,:)                      !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8),target,allocatable ::  DPTHZ_vr(:,:,:)                       !depth to middle of soil layer from  surface of grid cell [m]
  real(r8),target,allocatable ::  AREA(:,:,:,:)                      !cross-sectional area  [m2 d-2]
  real(r8),target,allocatable ::  DIST(:,:,:,:)                      !distance between adjacent layers:1=EW,2=NS,3=vertical [m]
  integer,target,allocatable ::  NU(:,:)                             !soil surface layer number
  integer,target,allocatable ::  NUI(:,:)                            !initial soil surface layer number
  integer,target,allocatable ::  MaxNumRootLays(:,:)                             !maximum root layer number
  integer,target,allocatable ::  NK(:,:)                             !additional soil lower boundary layers
  integer,target,allocatable ::  NLI(:,:)                            !initial lowest soil layer number
  integer,target,allocatable ::  NL(:,:)                             !lowest soil layer number
  integer,target,allocatable ::  NUM(:,:)                            !new surface layer number
  real(r8),target,allocatable ::  CumSoilDeptht0(:,:)                !initial depth from surface to bottom of soil layer [m]
  REAL(R8),target,allocatable ::  ALAT(:,:)                          !latitude	[degrees]
  real(r8),target,allocatable ::  DH(:,:)                            !number of EW grid cells, [-]
  real(r8),target,allocatable ::  DV(:,:)                            !number of EW grid cells, [-]
  integer,target,allocatable ::  FlowDirIndicator(:,:)                            !dimension of low
  integer,target,allocatable ::  LSG(:,:,:)                          !match PFT from different scenarios
  integer,target,allocatable ::  NP(:,:)                             !number of plant species
  integer,target,allocatable ::  NP0(:,:)                            !intitial number of plant species
!----------------------------------------------------------------------

contains

  subroutine InitGridData

  implicit none
  allocate(CumDepz2LayerBot_vr(0:JZ,JY,JX));  CumDepz2LayerBot_vr=0._r8
  allocate(DLYR(3,0:JZ,JY,JX)); DLYR=0._r8
  allocate(DLYRI_3D(3,0:JZ,JY,JX));DLYRI_3D=0._r8
  allocate(XDPTH(3,JZ,JY,JX));  XDPTH=0._r8
  allocate(SoiDepthMidLay_vr(JZ,JY,JX));     SoiDepthMidLay_vr=0._r8
  allocate(CumSoilThickness_vr(0:JZ,JY,JX)); CumSoilThickness_vr=0._r8
  allocate(DPTHZ_vr(JZ,JY,JX));    DPTHZ_vr=0._r8
  allocate(AREA(3,0:JZ,JY,JX)); AREA=0._r8
  allocate(DIST(3,JD,JV,JH));   DIST=0._r8
  allocate(NU(JY,JX));          NU=0
  allocate(NUI(JY,JX));         NUI=0
  allocate(MaxNumRootLays(JY,JX));          MaxNumRootLays=0
  allocate(NK(JY,JX));          NK=0
  allocate(NLI(JV,JH));         NLI=0
  allocate(NL(JV,JH));          NL=0
  allocate(NUM(JY,JX));         NUM=0
  allocate(CumSoilDeptht0(JY,JX));      CumSoilDeptht0=0._r8
  allocate(ALAT(JY,JX));        ALAT=0._r8
  allocate(DH(JY,JX));          DH=0._r8
  allocate(DV(JY,JX));          DV=0._r8
  allocate(FlowDirIndicator(JY,JX));         FlowDirIndicator=3   !vertical by default
  allocate(LSG(JZ,JY,JX));      LSG=0
  allocate(NP(JY,JX));          NP=0
  allocate(NP0(JY,JX));         NP0=0
  end subroutine InitGridData

!----------------------------------------------------------------------
  subroutine DestructGridData
  use abortutils, only : destroy
  implicit none

  call destroy(CumDepz2LayerBot_vr)
  call destroy(DLYR)
  call destroy(DLYRI_3D)
  call destroy(XDPTH)
  call destroy(SoiDepthMidLay_vr)
  call destroy(CumSoilThickness_vr)
  call destroy(DPTHZ_vr)
  call destroy(AREA)
  call destroy(DIST)
  call destroy(NU)
  call destroy(NUI)
  call destroy(MaxNumRootLays)
  call destroy(NK)
  call destroy(NLI)
  call destroy(NL)
  call destroy(NUM)
  call destroy(CumSoilDeptht0)
  call destroy(ALAT)
  call destroy(DH)
  call destroy(DV)
  call destroy(FlowDirIndicator)
  call destroy(LSG)
  call destroy(NP)
  call destroy(NP0)
  end subroutine DestructGridData

end module GridDataType
