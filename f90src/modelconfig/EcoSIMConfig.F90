module EcoSIMConfig

implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  logical :: is_restart_run=.false.
  logical :: is_first_year=.false.
  logical :: transport_on=.true.
  logical :: column_mode=.false.
  logical :: do_instequil=.false.
  integer :: IFLGW                 !flag for raising Z0G with vegn

  integer, parameter :: ndbiomcpc = 2 !# of microbial residue components
  integer, parameter :: nlbiomcpc = 3 !# of living biomass components
  integer, parameter :: jskenc    = 4 !# of kinetic components of the substrates
  integer, parameter :: jcplxc    = 5 !# of microbe-substrate complexes
  integer, parameter :: jcplx1c   = jcplxc-1
  integer, parameter :: NFGsc     = 7 !# of microbial functional groups in each complex
  character(len=36):: case_name
end module EcoSIMConfig
