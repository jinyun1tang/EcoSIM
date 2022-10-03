module FlagDataType
  use GridConsts
implicit none

  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__
  integer  ::  ISALTG                              !salt option
  integer  ::  IERSNG                              !erosion option
  integer  ::  ICLM                                !changes to weather data (0=none,1=step,2=transient)
  integer  ::  IMNG                                !flag for land management

  integer,allocatable ::  IYTYP(:,:,:,:)                      !fertilizer release type from fertilizer input file
  integer,allocatable ::  ITILL(:,:,:)                        !soil disturbance type, [-]
  integer,allocatable ::  IETYP(:,:)                          !Koppen climate zone
  integer,allocatable ::  IFLGV(:,:)                          !flag for irrigation criterion,0=SWC,1=canopy water potential
  integer,allocatable ::  IFLGS(:,:)                          !disturbance flag
  integer,allocatable ::  IFNHB(:,:)                          !banded NH4 fertilizer flag
  integer,allocatable ::  IFNOB(:,:)                          !banded NO3 fertilizer flag
  integer,allocatable ::  IFPOB(:,:)                          !banded H2PO4 fertilizer flag
  integer,allocatable ::  ISOIL(:,:,:,:)                      !flag for calculating FC(1),WP(2),SCNV(3),SCNH(4)
  integer,allocatable ::  ISOILR(:,:)                         !natural(0),reconstructed(1) soil profile
  integer,allocatable ::  IWTHR(:)                            !weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene
  integer,allocatable ::  IUTYP(:,:)                          !urea hydrolysis inhibitor type (1=no,2=yes)
  integer,allocatable ::  ITILL1(:,:)                         !soil disturbance type, [-]
  integer,allocatable ::  IFLGC(:,:,:)                        ! flag for living pft
  integer,allocatable ::  IFLGI(:,:,:)                        !PFT initialization flag:0=no,1=yes
  integer,allocatable ::  ICTYP(:,:,:)                        !plant photosynthetic type (C3 or C4)
  integer,allocatable ::  IGTYP(:,:,:)                        !plant growth type (vascular, non-vascular)
  integer,allocatable ::  ISTYP(:,:,:)                        !plant growth habit (annual or perennial)
  integer,allocatable ::  IDTYP(:,:,:)                        !plant growth habit (determinate or indeterminate)
  integer,allocatable ::  INTYP(:,:,:)                        !N2 fixation type
  integer,allocatable ::  IWTYP(:,:,:)                        !climate signal for phenological progress none, temperature, water stress)
  integer,allocatable ::  IPTYP(:,:,:)                        !photoperiod type (neutral, long day, short day)
  integer,allocatable ::  IBTYP(:,:,:)                        !phenologically-driven above-ground turnover (all, foliar only, none)
  integer,allocatable ::  IRTYP(:,:,:)                        !grain type (below or above-ground)
  integer,allocatable ::  MY(:,:,:)                           !mycorrhizal type (no or yes)
  integer,allocatable ::  IDTBL(:,:)                          !water table flag from site file
!----------------------------------------------------------------------

contains
  subroutine InitFlagData

  implicit none
  allocate(IYTYP(0:2,366,JY,JX));IYTYP=0
  allocate(ITILL(366,JY,JX));   ITILL=0
  allocate(IETYP(JY,JX));       IETYP=0
  allocate(IFLGV(JY,JX));       IFLGV=0
  allocate(IFLGS(JY,JX));       IFLGS=0
  allocate(IFNHB(JY,JX));       IFNHB=0
  allocate(IFNOB(JY,JX));       IFNOB=0
  allocate(IFPOB(JY,JX));       IFPOB=0
  allocate(ISOIL(4,JZ,JY,JX));  ISOIL=0
  allocate(ISOILR(JY,JX));      ISOILR=0
  allocate(IWTHR(2));           IWTHR=0
  allocate(IUTYP(JY,JX));       IUTYP=0
  allocate(ITILL1(JY,JX));      ITILL1=0
  allocate(IFLGC(JP,JY,JX));    IFLGC=0
  allocate(IFLGI(JP,JY,JX));    IFLGI=0
  allocate(ICTYP(JP,JY,JX));    ICTYP=0
  allocate(IGTYP(JP,JY,JX));    IGTYP=0
  allocate(ISTYP(JP,JY,JX));    ISTYP=0
  allocate(IDTYP(JP,JY,JX));    IDTYP=0
  allocate(INTYP(JP,JY,JX));    INTYP=0
  allocate(IWTYP(JP,JY,JX));    IWTYP=0
  allocate(IPTYP(JP,JY,JX));    IPTYP=0
  allocate(IBTYP(JP,JY,JX));    IBTYP=0
  allocate(IRTYP(JP,JY,JX));    IRTYP=0
  allocate(MY(JP,JY,JX));       MY=0
  allocate(IDTBL(JY,JX));       IDTBL=0
  end subroutine InitFlagData

!----------------------------------------------------------------------
  subroutine DestructFlagData
  use abortutils, only : destroy
  implicit none
  call destroy(IYTYP)
  call destroy(ITILL)
  call destroy(IETYP)
  call destroy(IFLGV)
  call destroy(IFLGS)
  call destroy(IFNHB)
  call destroy(IFNOB)
  call destroy(IFPOB)
  call destroy(ISOIL)
  call destroy(ISOILR)
  call destroy(IWTHR)
  call destroy(IUTYP)
  call destroy(ITILL1)
  call destroy(IFLGC)
  call destroy(IFLGI)
  call destroy(ICTYP)
  call destroy(IGTYP)
  call destroy(ISTYP)
  call destroy(IDTYP)
  call destroy(INTYP)
  call destroy(IWTYP)
  call destroy(IPTYP)
  call destroy(IBTYP)
  call destroy(IRTYP)
  call destroy(MY)
  call destroy(IDTBL)
  end subroutine DestructFlagData

end module FlagDataType
