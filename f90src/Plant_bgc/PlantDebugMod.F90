module PlantDebugMod

  use TracerIDMod
  use PlantAPIData

implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: PrintRootTracer
  contains

  subroutine PrintRootTracer(I,J,NZ,header)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), intent(in) :: header
  associate(                                     &
    NU              => plt_site%NU,              &
    NK              => plt_site%NK,              &
    MY_pft          => plt_morph%MY_pft,         &
    trcg_rootml_pvr => plt_rbgc%trcg_rootml_pvr, &
    trcs_rootml_pvr => plt_rbgc%trcs_rootml_pvr  &
  )
  return
  write(223,*)NZ,MY_pft(NZ),NU,NK
  write(223,*)I*1000+J,trim(header),sum(trcg_rootml_pvr(idg_CO2,1:MY_pft(NZ),NU:NK,NZ))+sum(trcg_rootml_pvr(idg_CO2,1:MY_pft(NZ),NU:NK,NZ))
  end associate
  end subroutine PrintRootTracer
end module PlantDebugMod
