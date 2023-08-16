module SoilPropertyDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  "SoilPropertyDataType.F90"

   real(r8) ,target,allocatable ::  CORGCI(:,:,:)                    !soil organic C content   [g kg-1]
   real(r8) ,target,allocatable ::  POROSI(:,:,:)                    !soil porosity            [m3 m-3]
   real(r8) ,target,allocatable ::  SoilFracAsMacPt0(:,:,:)                     !soil macropore fraction
   real(r8) ,target,allocatable ::  CSAND(:,:,:)                     !soil sand content [kg Mg-1]
   real(r8) ,target,allocatable ::  CSILT(:,:,:)                     !soil silt content [kg Mg-1]
   real(r8) ,target,allocatable ::  CCLAY(:,:,:)                     !soil clay content [kg Mg-1]
   real(r8) ,target,allocatable ::  ROCK(:,:,:)                      !Rock fraction
   real(r8) ,target,allocatable ::  SoiBulkDensityt0(:,:,:)                     !initial bulk density [Mg m-3,0=water]
   real(r8) ,target,allocatable ::  FracSoiAsMicP(:,:,:)            !micropore fraction
   real(r8) ,target,allocatable ::  SoilFracAsMacP(:,:,:)                      !macropore fraction
   real(r8) ,target,allocatable ::  PathLenMacP(:,:,:)                      !path length between macopores
   real(r8) ,target,allocatable ::  MacPRadius(:,:,:)                      !radius of macropores
   real(r8) ,target,allocatable ::  SoiBulkDensity(:,:,:)                      !soil bulk density, [Mg m-3]
   integer  ,target,allocatable ::  MacPNumLayer(:,:,:)                 !number of macropores
   real(r8) ,target,allocatable ::  POROS(:,:,:)                     !soil porosity
   real(r8) ,target,allocatable ::  VLSoilPoreMicP(:,:,:)            !micropore volume of soil layer	m3 d-2
   real(r8) ,target,allocatable ::  VLSoilMicP(:,:,:)                      !micropore volume
   real(r8) ,target,allocatable ::  SoilMicPMassLayer(:,:,:)                      !mass of soil layer	Mg d-2
   real(r8) ,target,allocatable ::  SoilMicPMassLayerMn(:,:)                      !minimum soil layer mass
   real(r8) ,target,allocatable ::  SoilMicPMassLayerMX(:,:)                      !maximum soil layer mass
   real(r8) ,target,allocatable ::  SAND(:,:,:)                      !soil sand content	Mg d-2
   real(r8) ,target,allocatable ::  SILT(:,:,:)                      !soil silt content	Mg d-2
   real(r8) ,target,allocatable ::  CLAY(:,:,:)                      !soil clay content	Mg d-2
   real(r8) ,target,allocatable ::  VLMicP(:,:,:)                    !total micropore volume in layer
   real(r8) ,target,allocatable ::  VLMacP(:,:,:)                    !total macropore volume in layer
   real(r8) ,target,allocatable ::  VGeomLayer(:,:,:)                      !soil volume including  macropores+rock [m3 d-2]
   real(r8) ,target,allocatable ::  VGeomLayert0(:,:,:)                     !initial soil volume including  macropores+rock [m3 d-2]
  private :: InitAllocate

contains


  subroutine InitSoilProperty

  implicit none

  call InitAllocate

  end subroutine InitSoilProperty

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(CORGCI(JZ,JY,JX));    CORGCI=0._r8
  allocate(POROSI(0:JZ,JY,JX));  POROSI=0._r8
  allocate(SoilFracAsMacPt0(JZ,JY,JX));     SoilFracAsMacPt0=0._r8
  allocate(CSAND(JZ,JY,JX));     CSAND=0._r8
  allocate(CSILT(JZ,JY,JX));     CSILT=0._r8
  allocate(CCLAY(JZ,JY,JX));     CCLAY=0._r8
  allocate(ROCK(JZ,JY,JX));      ROCK=0._r8
  allocate(SoiBulkDensityt0(JZ,JY,JX));     SoiBulkDensityt0=0._r8
  allocate(FracSoiAsMicP(0:JZ,JY,JX));    FracSoiAsMicP=0._r8
  allocate(SoilFracAsMacP(JZ,JY,JX));      SoilFracAsMacP=0._r8
  allocate(PathLenMacP(JZ,JY,JX));      PathLenMacP=0._r8
  allocate(MacPRadius(JZ,JY,JX));      MacPRadius=0._r8
  allocate(SoiBulkDensity(0:JZ,JY,JX));    SoiBulkDensity=0._r8
  allocate(MacPNumLayer(JZ,JY,JX));      MacPNumLayer=0
  allocate(POROS(0:JZ,JY,JX));   POROS=0._r8
  allocate(VLSoilPoreMicP(0:JZ,JY,JX));    VLSoilPoreMicP=0._r8
  allocate(VLSoilMicP(0:JZ,JY,JX));    VLSoilMicP=0._r8
  allocate(SoilMicPMassLayer(0:JZ,JY,JX));    SoilMicPMassLayer=0._r8
  allocate(SoilMicPMassLayerMn(JY,JX));       SoilMicPMassLayerMn=0._r8
  allocate(SoilMicPMassLayerMX(JY,JX));       SoilMicPMassLayerMX=0._r8
  allocate(SAND(JZ,JY,JX));      SAND=0._r8
  allocate(SILT(JZ,JY,JX));      SILT=0._r8
  allocate(CLAY(JZ,JY,JX));      CLAY=0._r8
  allocate(VLMicP(0:JZ,JY,JX));    VLMicP=0._r8
  allocate(VLMacP(JZ,JY,JX));     VLMacP=0._r8
  allocate(VGeomLayer(0:JZ,JY,JX));    VGeomLayer=0._r8
  allocate(VGeomLayert0(0:JZ,JY,JX));   VGeomLayert0=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilProperty

  use abortutils, only : destroy
  implicit none
  call destroy(CORGCI)
  call destroy(POROSI)
  call destroy(SoilFracAsMacPt0)
  call destroy(CSAND)
  call destroy(CSILT)
  call destroy(CCLAY)
  call destroy(ROCK)
  call destroy(SoiBulkDensityt0)
  call destroy(FracSoiAsMicP)
  call destroy(SoilFracAsMacP)
  call destroy(PathLenMacP)
  call destroy(MacPRadius)
  call destroy(SoiBulkDensity)
  call destroy(MacPNumLayer)
  call destroy(POROS)
  call destroy(VLSoilPoreMicP)
  call destroy(VLSoilMicP)
  call destroy(SoilMicPMassLayer)
  call destroy(SoilMicPMassLayerMn)
  call destroy(SoilMicPMassLayerMX)
  call destroy(SAND)
  call destroy(SILT)
  call destroy(CLAY)
  call destroy(VLMicP)
  call destroy(VLMacP)
  call destroy(VGeomLayer)
  call destroy(VGeomLayert0)
  end subroutine DestructSoilProperty

end module SoilPropertyDataType
