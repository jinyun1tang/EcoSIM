module PlantDebugMod

  use TracerIDMod
  use PlantAPIData

implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: PrintRootTracer
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine PrintRootTracer(I,J,NZ,header)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), intent(in) :: header
  associate(                                      &
    NU              => plt_site%NU               ,& !input  :current soil surface layer number, [-]
    NK              => plt_site%NK               ,& !input  :current hydrologically active layer, [-]
    Myco_pft        => plt_morph%Myco_pft        ,& !input  :mycorrhizal type (no or yes),[-]
    trcg_rootml_pvr => plt_rbgc%trcg_rootml_pvr   & !input  :root gas content, [g d-2]
  )
  return
  write(223,*)NZ,Myco_pft(NZ),NU,NK
  write(223,*)I*1000+J,trim(header),sum(trcg_rootml_pvr(idg_CO2,1:Myco_pft(NZ),NU:NK,NZ))+sum(trcg_rootml_pvr(idg_CO2,1:Myco_pft(NZ),NU:NK,NZ))
  end associate
  end subroutine PrintRootTracer
  ![tail]
end module PlantDebugMod
