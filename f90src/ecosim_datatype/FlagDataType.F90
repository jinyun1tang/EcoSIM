module FlagDataType
  use GridDataType
implicit none

  public
  save

  integer :: IETYP(JY,JX)         !Koppen climate zone
  integer :: IFLGV(JY,JX)         !flag for irrigation criterion,0=SWC,1=canopy water potential
  integer :: IFLGS(JY,JX)         !disturbance flag
  integer :: IFNHB(JY,JX)         !banded NH4 fertilizer flag
  integer :: IFNOB(JY,JX)         !banded NO3 fertilizer flag
  integer :: IFPOB(JY,JX)         !banded H2PO4 fertilizer flag
  integer :: ISALTG               !salt option
  integer :: IERSNG               !erosion option
  integer :: ICLM                 !changes to weather data (0=none,1=step,2=transient)
  integer :: IMNG                 !flag for land management
  integer :: IFLGW                !flag for raising Z0G with vegn
  integer :: ISOIL(4,JZ,JY,JX)    !flag for calculating FC(1),WP(2),SCNV(3),SCNH(4)
  integer :: ISOILR(JY,JX)        !natural(0),reconstructed(1) soil profile
  integer :: IWTHR(2)             !weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene
  integer :: IUTYP(JY,JX)         !urea hydrolysis inhibitor type (1=no,2=yes)
  integer :: IYTYP(0:2,366,JY,JX) !fertilizer release type from fertilizer input file


  integer :: IFLGC(JP,JY,JX)      ! flag for living pft
  integer :: IFLGI(JP,JY,JX)      !PFT initialization flag:0=no,1=yes
  integer :: ICTYP(JP,JY,JX)
  integer :: IGTYP(JP,JY,JX)
  integer :: ISTYP(JP,JY,JX)
  integer :: IDTYP(JP,JY,JX)
  integer :: INTYP(JP,JY,JX)
  integer :: IWTYP(JP,JY,JX)
  integer :: IPTYP(JP,JY,JX)
  integer :: IBTYP(JP,JY,JX)
  integer :: IRTYP(JP,JY,JX)
  integer :: MY(JP,JY,JX)



end module FlagDataType
