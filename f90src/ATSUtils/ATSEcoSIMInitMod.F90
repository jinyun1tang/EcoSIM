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
  use GridMod           , only : SetMeshATS
  use InitAllocMod
  use StartsMod, only : startsim
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS
  integer, intent(in) :: NYS
  real(r8) :: YSIN(NumOfSkyAzimuSects),YCOS(NumOfSkyAzimuSects),SkyAzimuthAngle(NumOfSkyAzimuSects)

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
    !write(*,*) "For loop: ", NY, " out of ", NYS

    !write(*,*) "For variable NU:"
    !write(*,*) "NU(NY,NX) = ", NU(NY,NX)
    !write(*,*) "a_NU(NY)  = ", a_NU(NY)
    !NU(NY,NX)=a_NU(NY)

    !write(*,*) "For variable NL:"
    !write(*,*) "NL(NY,NX) = ", NL(NY,NX)
    !write(*,*) "a_NL(NY)  = ", a_NL(NY)
    !NL(NY,NX)=a_NL(NY)
  end do

  do NY=1,NYS
    !write(*,*) "For loop: ", NY, " out of ", NYS 
    
    !write(*,*) "For variable NU:"
    !write(*,*) "NU(NY,NX) = ", NU(NY,NX)
    !write(*,*) "a_NU(NY)  = ", a_NU(NY)
    NU(NY,NX)=a_NU(NY)
    
    !write(*,*) "For variable NL:"
    !write(*,*) "NL(NY,NX) = ", NL(NY,NX)
    !write(*,*) "a_NL(NY)  = ", a_NL(NY)    
    NL(NY,NX)=a_NL(NY)
    
    !write(*,*) "For variable AREA:"
    !write(*,*) "AREA(3,0,NY,NX) = ", AREA(3,0,NY,NX)
    !write(*,*) "a_NU(NY)  = ", a_AREA3(NY)


    AREA(3,0,NY,NX)=a_AREA3(NY)

    !write(*,*) "For variable AREA:"
    !write(*,*) "AREA(3,NU(NY,NX),NY,NX) = ", AREA(3,NU(NY,NX),NY,NX)
    !write(*,*) "a_NU(NY)  = ", a_AREA3(NY)

    AREA(3,NU(NY,NX),NY,NX)=a_AREA3(NY)

    !write(*,*) "For variable ASP:"
    !write(*,*) "ASP(NY,NX) = ", ASP(NY,NX)
    !write(*,*) "a_ASP(NY)  = ", a_ASP(NY)

    ASP(NY,NX)=a_ASP(NY)

    !write(*,*) "For variable Tair:"
    !write(*,*) "TairKClimMean(NY,NX) = ", TairKClimMean(NY,NX)
    !write(*,*) "a_ATKA(NY)  = ", a_ATKA(NY)

    TairKClimMean(NY,NX)=a_ATKA(NY)
    CO2E(NY,NX)=atm_co2
    CH4E(NY,NX)=atm_ch4
    OXYE(NY,NX)=atm_o2
    Z2GE(NY,NX)=atm_n2
    Z2OE(NY,NX)=atm_n2o
    ZNH3E(NY,NX)=atm_nh3
    H2GE(NY,NX)=atm_H2
  
    DO L=NU(NY,NX),NL(NY,NX)
      !write(*,*) "On Loop L: ", L, ", From ", NU(NY,NX), " to ", NL(NY,NX)

      !write(*,*) "For variable FC:"
      !write(*,*) "FieldCapacity(L,NY,NX) = ", FieldCapacity(L,NY,NX)
      !write(*,*) "a_FC(L,NY)  = ", a_FC(L,NY)

      FieldCapacity(L,NY,NX)=a_FC(L,ny)

      !write(*,*) "For variable WP:"
      !write(*,*) "WiltPoint(L,NY,NX) = ", WiltPoint(L,NY,NX)
      !write(*,*) "a_WP(L,NY)  = ", a_WP(L,NY)

      WiltPoint(L,NY,NX)=a_WP(L,NY)

      !write(*,*) "For variable CD2LB:"
      !write(*,*) "CumDepth2LayerBottom(L,NY,NX) = ", CumDepth2LayerBottom(L,NY,NX)
      !write(*,*) "a_CumDepth2LayerBottom(L,NY)  = ", a_CumDepth2LayerBottom(L,NY)

      CumDepth2LayerBottom(L,NY,NX)=a_CumDepth2LayerBottom(L,NY)

      !write(*,*) "For variable SBD:"
      !write(*,*) "SoiBulkDensityt0(L,NY,NX) = ", SoiBulkDensityt0(L,NY,NX)
      !write(*,*) "a_BKDSI(L,NY)  = ", a_BKDSI(L,NY)

      SoiBulkDensityt0(L,NY,NX)=a_BKDSI(L,NY)

      !write(*,*) "For variable CORGC:"
      !write(*,*) "CORGC(L,NY,NX) = ", CORGC(L,NY,NX)
      !write(*,*) "a_CORGC(L,NY)  = ", a_CORGC(L,NY)

      CORGC(L,NY,NX)=a_CORGC(L,NY)

      !write(*,*) "For variable CORGN:"
      !write(*,*) "CORGN(L,NY,NX) = ", CORGN(L,NY,NX)
      !write(*,*) "a_CORGN(L,NY)  = ", a_CORGN(L,NY)

      CORGN(L,NY,NX)=a_CORGN(L,NY)

      !write(*,*) "For variable CORGP:"
      !write(*,*) "CORGP(L,NY,NX) = ", CORGP(L,NY,NX)
      !write(*,*) "a_CORGP(L,NY)  = ", a_CORGP(L,NY)

      CORGP(L,NY,NX)=a_CORGP(L,NY)
    ENDDO
  ENDDO

  write(*,*) "Finished loop, starting simulation..."
  !call startsim(NHW,NHE,NVN,NVS)

  write(*,*) "Finished Subroutine"
  end subroutine Init_EcoSIM_Soil


end module ATSEcoSIMInitMod
