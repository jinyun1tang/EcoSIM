module FlagDataType
  use GridConsts
implicit none

  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  integer :: IETYP(JY,JX)          !Koppen climate zone
  integer :: IFLGV(JY,JX)          !flag for irrigation criterion,0=SWC,1=canopy water potential
  integer :: IFLGS(JY,JX)          !disturbance flag
  integer :: IFNHB(JY,JX)          !banded NH4 fertilizer flag
  integer :: IFNOB(JY,JX)          !banded NO3 fertilizer flag
  integer :: IFPOB(JY,JX)          !banded H2PO4 fertilizer flag
  integer :: ISALTG                !salt option
  integer :: IERSNG                !erosion option
  integer :: ICLM                  !changes to weather data (0=none,1=step,2=transient)
  integer :: IMNG                  !flag for land management
  integer :: IFLGW                 !flag for raising Z0G with vegn
  integer :: ISOIL(4,JZ,JY,JX)     !flag for calculating FC(1),WP(2),SCNV(3),SCNH(4)
  integer :: ISOILR(JY,JX)         !natural(0),reconstructed(1) soil profile
  integer :: IWTHR(2)              !weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene
  integer :: IUTYP(JY,JX)          !urea hydrolysis inhibitor type (1=no,2=yes)
  integer :: IYTYP(0:2,366,JY,JX)  !fertilizer release type from fertilizer input file
  integer :: ITILL(366,JY,JX)      !soil disturbance type, [-]

  integer :: IFLGC(JP,JY,JX)       ! flag for living pft
  integer :: IFLGI(JP,JY,JX)       !PFT initialization flag:0=no,1=yes
  integer :: ICTYP(JP,JY,JX)       !plant photosynthetic type (C3 or C4)
  integer :: IGTYP(JP,JY,JX)       !plant growth type (vascular, non-vascular)
  integer :: ISTYP(JP,JY,JX)       !plant growth habit (annual or perennial)
  integer :: IDTYP(JP,JY,JX)       !plant growth habit (determinate or indeterminate)
  integer :: INTYP(JP,JY,JX)       !N2 fixation type
  integer :: IWTYP(JP,JY,JX)       !climate signal for phenological progress none, temperature, water stress)
  integer :: IPTYP(JP,JY,JX)       !photoperiod type (neutral, long day, short day)
  integer :: IBTYP(JP,JY,JX)       !phenologically-driven above-ground turnover (all, foliar only, none)
  integer :: IRTYP(JP,JY,JX)       !grain type (below or above-ground)
  integer :: MY(JP,JY,JX)          !mycorrhizal type (no or yes)
  integer :: IDTBL(JY,JX)         !water table flag from site file


end module FlagDataType
