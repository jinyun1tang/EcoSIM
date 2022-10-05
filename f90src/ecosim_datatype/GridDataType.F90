module GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__


  integer :: NPX                  !number of E-W grid cells
  integer :: NPY                  !number of N-S grid cells
  real(r8) :: TAREA               !total area of landscape	[m2]
  real(r8),allocatable ::  CDPTH(:,:,:)                       !depth to bottom of soil layer [m]
  real(r8),allocatable ::  DLYR(:,:,:,:)                      !thickness of soil layer [m]
  real(r8),allocatable ::  DLYRI(:,:,:,:)                     !thickness of soil layer [m]
  real(r8),allocatable ::  XDPTH(:,:,:,:)                     !cross-sectional area / distance between adjacent grid cells [m]
  real(r8),allocatable ::  DPTH(:,:,:)                        !depth to middle of soil layer [m]
  real(r8),allocatable ::  CDPTHZ(:,:,:)                      !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8),allocatable ::  DPTHZ(:,:,:)                       !depth to middle of soil layer from  surface of grid cell [m]
  real(r8),allocatable ::  AREA(:,:,:,:)                      !cross-sectional area  [m2 d-2]
  real(r8),allocatable ::  DIST(:,:,:,:)                      !distance between adjacent layers:1=EW,2=NS,3=vertical [m]
  integer,allocatable ::  NU(:,:)                             !soil surface layer number
  integer,allocatable ::  NUI(:,:)                            !initial soil surface layer number
  integer,allocatable ::  NJ(:,:)                             !maximum root layer number
  integer,allocatable ::  NK(:,:)                             !additional soil lower boundary layers
  integer,allocatable ::  NLI(:,:)                            !initial lowest soil layer number
  integer,allocatable ::  NL(:,:)                             !lowest soil layer number
  integer,allocatable ::  NUM(:,:)                            !new surface layer number
  real(r8),allocatable ::  CDPTHI(:,:)                        !initial depth to bottom of soil layer [m]
  REAL(R8),allocatable ::  ALAT(:,:)                          !latitude	[degrees]
  real(r8),allocatable ::  DH(:,:)                            !number of EW grid cells, [-]
  real(r8),allocatable ::  DV(:,:)                            !number of EW grid cells, [-]
  integer,allocatable ::  NCN(:,:)                            !number of dimensions for grid cell connections
  integer,allocatable ::  LSG(:,:,:)                          !match PFT from different scenarios
  integer,allocatable ::  NP(:,:)                             !number of plant species
  integer,allocatable ::  NP0(:,:)                            !intitial number of plant species
!----------------------------------------------------------------------

contains

  subroutine InitGridData

  implicit none
  allocate(CDPTH(0:JZ,JY,JX));  CDPTH=0._r8
  allocate(DLYR(3,0:JZ,JY,JX)); DLYR=0._r8
  allocate(DLYRI(3,0:JZ,JY,JX));DLYRI=0._r8
  allocate(XDPTH(3,JZ,JY,JX));  XDPTH=0._r8
  allocate(DPTH(JZ,JY,JX));     DPTH=0._r8
  allocate(CDPTHZ(0:JZ,JY,JX)); CDPTHZ=0._r8
  allocate(DPTHZ(JZ,JY,JX));    DPTHZ=0._r8
  allocate(AREA(3,0:JZ,JY,JX)); AREA=0._r8
  allocate(DIST(3,JD,JV,JH));   DIST=0._r8
  allocate(NU(JY,JX));          NU=0
  allocate(NUI(JY,JX));         NUI=0
  allocate(NJ(JY,JX));          NJ=0
  allocate(NK(JY,JX));          NK=0
  allocate(NLI(JV,JH));         NLI=0
  allocate(NL(JV,JH));          NL=0
  allocate(NUM(JY,JX));         NUM=0
  allocate(CDPTHI(JY,JX));      CDPTHI=0._r8
  allocate(ALAT(JY,JX));        ALAT=0._r8
  allocate(DH(JY,JX));          DH=0._r8
  allocate(DV(JY,JX));          DV=0._r8
  allocate(NCN(JY,JX));         NCN=0
  allocate(LSG(JZ,JY,JX));      LSG=0
  allocate(NP(JY,JX));          NP=0
  allocate(NP0(JY,JX));         NP0=0
  end subroutine InitGridData

!----------------------------------------------------------------------
  subroutine DestructGridData
  use abortutils, only : destroy
  implicit none

  call destroy(CDPTH)
  call destroy(DLYR)
  call destroy(DLYRI)
  call destroy(XDPTH)
  call destroy(DPTH)
  call destroy(CDPTHZ)
  call destroy(DPTHZ)
  call destroy(AREA)
  call destroy(DIST)
  call destroy(NU)
  call destroy(NUI)
  call destroy(NJ)
  call destroy(NK)
  call destroy(NLI)
  call destroy(NL)
  call destroy(NUM)
  call destroy(CDPTHI)
  call destroy(ALAT)
  call destroy(DH)
  call destroy(DV)
  call destroy(NCN)
  call destroy(LSG)
  call destroy(NP)
  call destroy(NP0)
  end subroutine DestructGridData

end module GridDataType
