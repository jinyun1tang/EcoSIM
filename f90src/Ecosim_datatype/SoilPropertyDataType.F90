module SoilPropertyDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

   real(r8) ,target,allocatable ::  CORGCI_vr(:,:,:)                    !soil organic C content   [g kg-1]
   real(r8) ,target,allocatable ::  POROSI_vr(:,:,:)                    !soil porosity            [m3 m-3]
   real(r8) ,target,allocatable ::  SoilFracAsMacPt0_vr(:,:,:)          !soil macropore fraction
   real(r8) ,target,allocatable ::  CSAND_vr(:,:,:)                     !soil sand content [kg Mg-1]
   real(r8) ,target,allocatable ::  CSILT(:,:,:)                        !soil silt content [kg Mg-1]
   real(r8) ,target,allocatable ::  CCLAY_vr(:,:,:)                     !soil clay content [kg Mg-1]
   real(r8) ,target,allocatable ::  ROCK_vr(:,:,:)                      !Rock fraction, [0-1]
   real(r8) ,target,allocatable ::  SoiBulkDensityt0_vr(:,:,:)          !initial bulk density [Mg m-3,0=water]
   real(r8) ,target,allocatable ::  FracSoiAsMicP_vr(:,:,:)             !micropore fraction, [0-1]
   real(r8) ,target,allocatable ::  SoilFracAsMacP_vr(:,:,:)            !macropore fraction, [0-1]
   real(r8) ,target,allocatable ::  PathLenMacP(:,:,:)                  !path length between macopores
   real(r8) ,target,allocatable ::  MacPRadius(:,:,:)                   !radius of macropores
   real(r8) ,target,allocatable ::  SoilBulkDensity_vr(:,:,:)           !soil bulk density, [Mg m-3]
   integer  ,target,allocatable ::  MacPNumLayer(:,:,:)                 !number of macropores
   real(r8) ,target,allocatable ::  POROS_vr(:,:,:)                     !soil porosity [m3 m-3]
   real(r8) ,target,allocatable ::  VLSoilPoreMicP_vr(:,:,:)            !Volume of soil occupied by micropores	[m3 d-2]
   real(r8) ,target,allocatable ::  VLSoilMicP_vr(:,:,:)                !volume of micropores [m3 d-2]
   real(r8) ,target,allocatable ::  VLSoilMicPMass_vr(:,:,:)            !mass of soil layer	[Mg d-2]
   real(r8) ,target,allocatable ::  SoilMicPMassLayerMn(:,:)            !minimum soil layer mass [Mg d-2]
   real(r8) ,target,allocatable ::  SoilMicPMassLayerMX(:,:)            !maximum soil layer mass [Mg d-2]
   real(r8) ,target,allocatable ::  SAND(:,:,:)                         !soil sand content	[Mg d-2]
   real(r8) ,target,allocatable ::  SILT(:,:,:)                         !soil silt content	[Mg d-2]
   real(r8) ,target,allocatable ::  CLAY(:,:,:)                         !soil clay content	[Mg d-2]
   real(r8) ,target,allocatable ::  VLMicP_vr(:,:,:)                    !total micropore volume in layer [m3 d-2]
   real(r8) ,target,allocatable ::  VLMacP_vr(:,:,:)                    !total macropore volume in layer [m3 d-2]
   real(r8) ,target,allocatable ::  VGeomLayer_vr(:,:,:)                !soil volume including  macropores+rock [m3 d-2]
   real(r8) ,target,allocatable ::  VGeomLayert0_vr(:,:,:)              !initial soil volume including  macropores+rock [m3 d-2]
   real(r8) ,target,allocatable ::  VOLTX_vr(:,:,:)                     !maximum soil pore (mac+mic) volume allowed [m3 d-2]  
  private :: InitAllocate

contains


  subroutine InitSoilProperty

  implicit none

  call InitAllocate

  end subroutine InitSoilProperty

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(CORGCI_vr(JZ,JY,JX));    CORGCI_vr=0._r8
  allocate(POROSI_vr(0:JZ,JY,JX));  POROSI_vr=0._r8
  allocate(SoilFracAsMacPt0_vr(JZ,JY,JX));     SoilFracAsMacPt0_vr=0._r8
  allocate(CSAND_vr(JZ,JY,JX));     CSAND_vr=0._r8
  allocate(CSILT(JZ,JY,JX));     CSILT=0._r8
  allocate(CCLAY_vr(JZ,JY,JX));     CCLAY_vr=0._r8
  allocate(ROCK_vr(JZ,JY,JX));      ROCK_vr=0._r8
  allocate(SoiBulkDensityt0_vr(JZ,JY,JX));     SoiBulkDensityt0_vr=0._r8
  allocate(FracSoiAsMicP_vr(0:JZ,JY,JX));    FracSoiAsMicP_vr=0._r8
  allocate(SoilFracAsMacP_vr(JZ,JY,JX));      SoilFracAsMacP_vr=0._r8
  allocate(PathLenMacP(JZ,JY,JX));      PathLenMacP=0._r8
  allocate(MacPRadius(JZ,JY,JX));      MacPRadius=0._r8
  allocate(SoilBulkDensity_vr(0:JZ,JY,JX));    SoilBulkDensity_vr=0._r8
  allocate(MacPNumLayer(JZ,JY,JX));      MacPNumLayer=0
  allocate(POROS_vr(0:JZ,JY,JX));   POROS_vr=0._r8
  allocate(VLSoilPoreMicP_vr(0:JZ,JY,JX));    VLSoilPoreMicP_vr=0._r8
  allocate(VLSoilMicP_vr(0:JZ,JY,JX));    VLSoilMicP_vr=0._r8
  allocate(VLSoilMicPMass_vr(0:JZ,JY,JX));    VLSoilMicPMass_vr=0._r8
  allocate(SoilMicPMassLayerMn(JY,JX));       SoilMicPMassLayerMn=0._r8
  allocate(SoilMicPMassLayerMX(JY,JX));       SoilMicPMassLayerMX=0._r8
  allocate(SAND(JZ,JY,JX));      SAND=0._r8
  allocate(SILT(JZ,JY,JX));      SILT=0._r8
  allocate(CLAY(JZ,JY,JX));      CLAY=0._r8
  allocate(VLMicP_vr(0:JZ,JY,JX));    VLMicP_vr=0._r8
  allocate(VLMacP_vr(JZ,JY,JX));     VLMacP_vr=0._r8
  allocate(VGeomLayer_vr(0:JZ,JY,JX));    VGeomLayer_vr=0._r8
  allocate(VGeomLayert0_vr(0:JZ,JY,JX));   VGeomLayert0_vr=0._r8
  allocate(VOLTX_vr(JZ,JY,JX));  VOLTX_vr=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilProperty

  use abortutils, only : destroy
  implicit none
  call destroy(CORGCI_vr)
  call destroy(POROSI_vr)
  call destroy(SoilFracAsMacPt0_vr)
  call destroy(CSAND_vr)
  call destroy(CSILT)
  call destroy(CCLAY_vr)
  call destroy(ROCK_vr)
  call destroy(SoiBulkDensityt0_vr)
  call destroy(FracSoiAsMicP_vr)
  call destroy(SoilFracAsMacP_vr)
  call destroy(PathLenMacP)
  call destroy(MacPRadius)
  call destroy(SoilBulkDensity_vr)
  call destroy(MacPNumLayer)
  call destroy(POROS_vr)
  call destroy(VLSoilPoreMicP_vr)
  call destroy(VLSoilMicP_vr)
  call destroy(VLSoilMicPMass_vr)
  call destroy(SoilMicPMassLayerMn)
  call destroy(SoilMicPMassLayerMX)
  call destroy(SAND)
  call destroy(SILT)
  call destroy(CLAY)
  call destroy(VLMicP_vr)
  call destroy(VLMacP_vr)
  call destroy(VGeomLayer_vr)
  call destroy(VGeomLayert0_vr)
  call destroy(VOLTX_vr)
  end subroutine DestructSoilProperty

end module SoilPropertyDataType
