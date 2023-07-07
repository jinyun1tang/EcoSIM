Module SharedDataMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  implicit none
  character(len=*), private, parameter :: mod_filename=__FILE__
!figure out the grid for ATS
  integer :: JZSOI   !number of soil layers
  integer :: JSNO    !number of snow layers
  
! temporary data holder in ecosim
  real(r8) :: atm_n2, atm_o2,atm_co2,atm_ch4,atm_N2o,atm_H2,atm_NH3
  real(r8), allocatable :: a_csand(:,:)   !sand mass fraction
  real(r8), allocatable :: a_CSILT(:,:)   !silt mass fraction
  real(r8), allocatable :: a_BKDSI(:,:)   !bulk density
  real(r8), allocatable :: a_CDPTH(:,:)   !dpeth (from surfce to bottom)
  real(r8), allocatable :: a_FC(:,:)      !field capacity
  real(r8), allocatable :: a_WP(:,:)      !wilting point
  real(r8), allocatable :: a_FHOL(:,:)    !macropore fraction
  real(r8), allocatable :: a_ROCK(:,:)    !mass fraction as rock
  real(r8), allocatable :: a_CORGC(:,:)   !organic carbon content
  real(r8), allocatable :: a_CORGN(:,:)   !organic nitrogen content
  real(r8), allocatable :: a_CORGP(:,:)   !organic phosphorus content
  real(r8), allocatable :: a_poros(:,:)
!  real(r8), allocatable ::a_CORGR(:,:)  !organic nitrogen  content
  real(r8), allocatable :: a_ASP(:)
  real(r8), allocatable :: a_ALT(:)
  real(r8), allocatable :: a_ATKA(:)

  real(r8), allocatable :: tairc(:)  !air temperature oC
  real(r8), allocatable :: uwind(:)  !wind speed, m/s
  real(r8), allocatable :: prec(:)   !precipitation, mm H2O/hr
  real(r8), allocatable :: srad(:)   !solar radiation,
  real(r8), allocatable :: vpa(:)    !vapor pressure deficit
  real(r8), allocatable :: a_AREA3(:)  
  integer,  allocatable :: a_NU(:)
  integer,  allocatable :: a_NL(:)
  integer,  allocatable :: a_NJ(:)
  integer :: NYS                     !total number of columns
  contains

  subroutine InitSharedData(JZs,ncol)
  use GridConsts, only : JX,JY,JZ
  implicit none
  integer, intent(in) :: JZs    !number of vertical layers
  integer, intent(in) :: NCOL  !NUMBER of cols
  !set # of soil layers
  JZSOI=JZs
  
  JX=1;JY=ncol;jz=jzs
  allocate(a_csand(JZSOI,ncol))  
  allocate(a_CSILT(jzsoi,ncol))   !silt mass fraction
  allocate(a_BKDSI(jzsoi,ncol))   !bulk density
  allocate(a_CDPTH(jzsoi,ncol))   !dpeth (from surfce to bottom)
  allocate(a_FC(jzsoi,ncol))      !field capacity
  allocate(a_WP(jzsoi,ncol))      !wilting point
  allocate(a_FHOL(jzsoi,ncol))    !macropore fraction
  allocate(a_ROCK(jzsoi,ncol))    !mass fraction as rock
  allocate(a_CORGC(jzsoi,ncol))   !organic carbon content
  allocate(a_CORGN(jzsoi,ncol))   !organic nitrogen content
  allocate(a_CORGP(jzsoi,ncol))   !organic phosphorus content
  allocate(a_poros(jzsoi,ncol))
  allocate(a_AREA3(ncol))
  allocate(a_NU(ncol))
  allocate(a_NL(ncol))
  allocate(a_ASP(ncol))
  allocate(a_ALT(ncol))
  allocate(tairc(1:ncol))
  allocate(uwind(1:ncol))
  allocate(prec(1:ncol))
  allocate(srad(1:ncol))
  allocate(vpa(1:ncol))  
  allocate(a_ATKA(1:ncol))
  end subroutine InitSharedData

!------------------------------------------------------------------------------------------

  subroutine DestroySharedData()
  implicit none

  call destroy(a_csand)  
  call destroy(a_CSILT)
  call destroy(a_BKDSI)
  call destroy(a_CDPTH)
  call destroy(a_FC)
  call destroy(a_WP)
  call destroy(a_FHOL)
  call destroy(a_ROCK)
  call destroy(a_CORGC)
  call destroy(a_CORGN)
  call destroy(a_CORGP)
  call destroy(a_poros)
  call destroy(a_AREA3)
  call destroy(a_ASP)
  call destroy(a_ALT)
  call destroy(tairc)
  call destroy(uwind)
  call destroy(prec)
  call destroy(srad)
  call destroy(vpa)  
  call destroy(a_ATKA)
  end subroutine DestroySharedData

end module SharedDataMod