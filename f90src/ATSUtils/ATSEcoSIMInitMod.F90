module ATSEcoSIMInitMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use ClimForcDataType
  use SoilPropertyDataType
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: Init_EcoSIM_Soil
  contains

  subroutine Init_EcoSIM_Soil(NYS)
  use EcoSimConst
  use GridMod           , only : SetMesh
  use InitAllocMod
  use StartsMod, only : startsim
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS
  integer, intent(in) :: NYS
  real(r8) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)

  NHW=1;NHE=1;NVN=1;NVS=NYS

  write(*,*) "In Init_EcoSIM_Soil"
  write(*,*) "Setting Mesh..."

  !call SetMesh(NHW,NVN,NHE,NVS)

  call SetMeshATS(NHW,NVN,NHE,NVS)

  write(*,*) "Finished mesh, Allocating ..."

  call InitAlloc(NOMicrobeGuilds=1)

  NX=1

  write(*,*) "Finished Allocate, beginning loop"

  do NY=1,NYS

    NU(NY,NX)=a_NU(NY)
    NL(NY,NX)=a_NL(NY)
    AREA(3,0,NY,NX)=a_AREA3(NY)
    AREA(3,NU(NY,NX),NY,NX)=a_AREA3(NY)
    ASP(NY,NX)=a_ASP(NY)
    TairKClimMean(NY,NX)=a_ATKA(NY)
    CO2E(NY,NX)=atm_co2
    CH4E(NY,NX)=atm_ch4
    OXYE(NY,NX)=atm_o2
    Z2GE(NY,NX)=atm_n2
    Z2OE(NY,NX)=atm_n2o
    ZNH3E(NY,NX)=atm_nh3
    H2GE(NY,NX)=atm_H2

    DO L=NU(NY,NX),NL(NY,NX)
      FieldCapacity(L,NY,NX)=a_FC(L,ny)
      WiltPoint(L,NY,NX)=a_WP(L,NY)
      CumDepth2LayerBottom(L,NY,NX)=a_CumDepth2LayerBottom(L,NY)
      SoiBulkDensityt0(L,NY,NX)=a_BKDSI(L,NY)
      CORGC(L,NY,NX)=a_CORGC(L,NY)
      CORGN(L,NY,NX)=a_CORGN(L,NY)
      CORGP(L,NY,NX)=a_CORGP(L,NY)
    ENDDO
  ENDDO

  write(*,*) "Finished loop, starting simulation..."
  call startsim(NHW,NHE,NVN,NVS)

  write(*,*) "Finished Subroutine"
  end subroutine Init_EcoSIM_Soil


end module ATSEcoSIMInitMod
