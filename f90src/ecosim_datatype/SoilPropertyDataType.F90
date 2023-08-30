module SoilPropertyDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

   real(r8) ,target,allocatable ::  CORGCI(:,:,:)                    !soil organic C content   [g kg-1]
   real(r8) ,target,allocatable ::  POROSI(:,:,:)                    !soil porosity            [m3 m-3]
   real(r8) ,target,allocatable ::  FHOLI(:,:,:)                     !soil macropore fraction
   real(r8) ,target,allocatable ::  CSAND(:,:,:)                     !soil sand content [kg Mg-1]
   real(r8) ,target,allocatable ::  CSILT(:,:,:)                     !soil silt content [kg Mg-1]
   real(r8) ,target,allocatable ::  CCLAY(:,:,:)                     !soil clay content [kg Mg-1]
   real(r8) ,target,allocatable ::  ROCK(:,:,:)                      !Rock fraction
   real(r8) ,target,allocatable ::  BKDSI(:,:,:)                     !initial bulk density [Mg m-3,0=water]
   real(r8) ,target,allocatable ::  FMPR(:,:,:)                      !micropore fraction
   real(r8) ,target,allocatable ::  FHOL(:,:,:)                      !macropore fraction
   real(r8) ,target,allocatable ::  PHOL(:,:,:)                      !path length between macopores
   real(r8) ,target,allocatable ::  HRAD(:,:,:)                      !radius of macropores
   real(r8) ,target,allocatable ::  BKDS(:,:,:)                      !soil bulk density, [Mg m-3]
   integer  ,target,allocatable ::  NHOL(:,:,:)                      !number of macropores
   real(r8) ,target,allocatable ::  POROS(:,:,:)                     !soil porosity
   real(r8) ,target,allocatable ::  VSoilPoreMicP(:,:,:)                      !volume of soil layer	m3 d-2
   real(r8) ,target,allocatable ::  VOLY(:,:,:)                      !micropore volume
   real(r8) ,target,allocatable ::  BKVL(:,:,:)                      !mass of soil layer	Mg d-2
   real(r8) ,target,allocatable ::  BKVLNM(:,:)                      !minimum soil layer mass
   real(r8) ,target,allocatable ::  BKVLNU(:,:)                      !maximum soil layer mass
   real(r8) ,target,allocatable ::  SAND(:,:,:)                      !soil sand content	Mg d-2
   real(r8) ,target,allocatable ::  SILT(:,:,:)                      !soil silt content	Mg d-2
   real(r8) ,target,allocatable ::  CLAY(:,:,:)                      !soil clay content	Mg d-2
   real(r8) ,target,allocatable ::  VMicP(:,:,:)                      !total volume in micropores
   real(r8) ,target,allocatable ::  VAirMacP(:,:,:)                     !total volume in macropores
   real(r8) ,target,allocatable ::  VOLT(:,:,:)                      !soil volume including  macropores+rock [m3 d-2]
   real(r8) ,target,allocatable ::  VOLTI(:,:,:)                     !initial soil volume including  macropores+rock [m3 d-2]
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
  allocate(FHOLI(JZ,JY,JX));     FHOLI=0._r8
  allocate(CSAND(JZ,JY,JX));     CSAND=0._r8
  allocate(CSILT(JZ,JY,JX));     CSILT=0._r8
  allocate(CCLAY(JZ,JY,JX));     CCLAY=0._r8
  allocate(ROCK(JZ,JY,JX));      ROCK=0._r8
  allocate(BKDSI(JZ,JY,JX));     BKDSI=0._r8
  allocate(FMPR(0:JZ,JY,JX));    FMPR=0._r8
  allocate(FHOL(JZ,JY,JX));      FHOL=0._r8
  allocate(PHOL(JZ,JY,JX));      PHOL=0._r8
  allocate(HRAD(JZ,JY,JX));      HRAD=0._r8
  allocate(BKDS(0:JZ,JY,JX));    BKDS=0._r8
  allocate(NHOL(JZ,JY,JX));      NHOL=0
  allocate(POROS(0:JZ,JY,JX));   POROS=0._r8
  allocate(VSoilPoreMicP(0:JZ,JY,JX));    VSoilPoreMicP=0._r8
  allocate(VOLY(0:JZ,JY,JX));    VOLY=0._r8
  allocate(BKVL(0:JZ,JY,JX));    BKVL=0._r8
  allocate(BKVLNM(JY,JX));       BKVLNM=0._r8
  allocate(BKVLNU(JY,JX));       BKVLNU=0._r8
  allocate(SAND(JZ,JY,JX));      SAND=0._r8
  allocate(SILT(JZ,JY,JX));      SILT=0._r8
  allocate(CLAY(JZ,JY,JX));      CLAY=0._r8
  allocate(VMicP(0:JZ,JY,JX));    VMicP=0._r8
  allocate(VAirMacP(JZ,JY,JX));     VAirMacP=0._r8
  allocate(VOLT(0:JZ,JY,JX));    VOLT=0._r8
  allocate(VOLTI(0:JZ,JY,JX));   VOLTI=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilProperty

  use abortutils, only : destroy
  implicit none
  call destroy(CORGCI)
  call destroy(POROSI)
  call destroy(FHOLI)
  call destroy(CSAND)
  call destroy(CSILT)
  call destroy(CCLAY)
  call destroy(ROCK)
  call destroy(BKDSI)
  call destroy(FMPR)
  call destroy(FHOL)
  call destroy(PHOL)
  call destroy(HRAD)
  call destroy(BKDS)
  call destroy(NHOL)
  call destroy(POROS)
  call destroy(VSoilPoreMicP)
  call destroy(VOLY)
  call destroy(BKVL)
  call destroy(BKVLNM)
  call destroy(BKVLNU)
  call destroy(SAND)
  call destroy(SILT)
  call destroy(CLAY)
  call destroy(VMicP)
  call destroy(VAirMacP)
  call destroy(VOLT)
  call destroy(VOLTI)
  end subroutine DestructSoilProperty

end module SoilPropertyDataType
