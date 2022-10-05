module SoilPropertyDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

   real(r8) ,allocatable ::  CORGCI(:,:,:)                    !soil organic C content   [g kg-1]
   real(r8) ,allocatable ::  POROSI(:,:,:)                    !soil porosity            [m3 m-3]
   real(r8) ,allocatable ::  FHOLI(:,:,:)                     !soil macropore fraction
   real(r8) ,allocatable ::  CSAND(:,:,:)                     !soil sand content [kg Mg-1]
   real(r8) ,allocatable ::  CSILT(:,:,:)                     !soil silt content [kg Mg-1]
   real(r8) ,allocatable ::  CCLAY(:,:,:)                     !soil clay content [kg Mg-1]
   real(r8) ,allocatable ::  ROCK(:,:,:)                      !Rock fraction
   real(r8) ,allocatable ::  BKDSI(:,:,:)                     !initial bulk density [Mg m-3,0=water]
   real(r8) ,allocatable ::  FMPR(:,:,:)                      !micropore fraction
   real(r8) ,allocatable ::  FHOL(:,:,:)                      !macropore fraction
   real(r8) ,allocatable ::  PHOL(:,:,:)                      !path length between macopores
   real(r8) ,allocatable ::  HRAD(:,:,:)                      !radius of macropores
   real(r8) ,allocatable ::  BKDS(:,:,:)                      !soil bulk density, [Mg m-3]
   integer  ,allocatable ::  NHOL(:,:,:)                      !number of macropores
   real(r8) ,allocatable ::  POROS(:,:,:)                     !soil porosity
   real(r8) ,allocatable ::  VOLX(:,:,:)                      !volume of soil layer	m3 d-2
   real(r8) ,allocatable ::  VOLY(:,:,:)                      !micropore volume
   real(r8) ,allocatable ::  BKVL(:,:,:)                      !mass of soil layer	Mg d-2
   real(r8) ,allocatable ::  BKVLNM(:,:)                      !minimum soil layer mass
   real(r8) ,allocatable ::  BKVLNU(:,:)                      !maximum soil layer mass
   real(r8) ,allocatable ::  SAND(:,:,:)                      !soil sand content	Mg d-2
   real(r8) ,allocatable ::  SILT(:,:,:)                      !soil silt content	Mg d-2
   real(r8) ,allocatable ::  CLAY(:,:,:)                      !soil clay content	Mg d-2
   real(r8) ,allocatable ::  VOLA(:,:,:)                      !total volume in micropores
   real(r8) ,allocatable ::  VOLAH(:,:,:)                     !total volume in macropores
   real(r8) ,allocatable ::  VOLT(:,:,:)                      !soil volume including  macropores+rock [m3 d-2]
   real(r8) ,allocatable ::  VOLTI(:,:,:)                     !initial soil volume including  macropores+rock [m3 d-2]
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
  allocate(VOLX(0:JZ,JY,JX));    VOLX=0._r8
  allocate(VOLY(0:JZ,JY,JX));    VOLY=0._r8
  allocate(BKVL(0:JZ,JY,JX));    BKVL=0._r8
  allocate(BKVLNM(JY,JX));       BKVLNM=0._r8
  allocate(BKVLNU(JY,JX));       BKVLNU=0._r8
  allocate(SAND(JZ,JY,JX));      SAND=0._r8
  allocate(SILT(JZ,JY,JX));      SILT=0._r8
  allocate(CLAY(JZ,JY,JX));      CLAY=0._r8
  allocate(VOLA(0:JZ,JY,JX));    VOLA=0._r8
  allocate(VOLAH(JZ,JY,JX));     VOLAH=0._r8
  allocate(VOLT(0:JZ,JY,JX));    VOLT=0._r8
  allocate(VOLTI(0:JZ,JY,JX));   VOLTI=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilProperty
  if (allocated(CORGCI))   deallocate(CORGCI)
  if (allocated(POROSI))   deallocate(POROSI)
  if (allocated(FHOLI))    deallocate(FHOLI)
  if (allocated(CSAND))    deallocate(CSAND)
  if (allocated(CSILT))    deallocate(CSILT)
  if (allocated(CCLAY))    deallocate(CCLAY)
  if (allocated(ROCK))     deallocate(ROCK)
  if (allocated(BKDSI))    deallocate(BKDSI)
  if (allocated(FMPR))     deallocate(FMPR)
  if (allocated(FHOL))     deallocate(FHOL)
  if (allocated(PHOL))     deallocate(PHOL)
  if (allocated(HRAD))     deallocate(HRAD)
  if (allocated(BKDS))     deallocate(BKDS)
  if (allocated(NHOL))     deallocate(NHOL)
  if (allocated(POROS))    deallocate(POROS)
  if (allocated(VOLX))     deallocate(VOLX)
  if (allocated(VOLY))     deallocate(VOLY)
  if (allocated(BKVL))     deallocate(BKVL)
  if (allocated(BKVLNM))   deallocate(BKVLNM)
  if (allocated(BKVLNU))   deallocate(BKVLNU)
  if (allocated(SAND))     deallocate(SAND)
  if (allocated(SILT))     deallocate(SILT)
  if (allocated(CLAY))     deallocate(CLAY)
  if (allocated(VOLA))     deallocate(VOLA)
  if (allocated(VOLAH))    deallocate(VOLAH)
  if (allocated(VOLT))     deallocate(VOLT)
  if (allocated(VOLTI))    deallocate(VOLTI)
  end subroutine DestructSoilProperty

end module SoilPropertyDataType
