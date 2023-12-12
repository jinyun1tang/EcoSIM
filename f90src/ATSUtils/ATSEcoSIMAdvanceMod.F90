module ATSEcoSIMAdvanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use GridConsts
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use CanopyDataType, only: SWRadOnGrnd
  !use PlantAPIData, only: CO2E, CH4E, OXYE, Z2GE, Z2OE, ZNH3E, &
  !    H2GE
  use ClimForcDataType, only : LWRadSky, TairK, &
      VPA, WindSpeedAtm, RainH  
  use SoilPropertyDataType
  use HydroThermData, only : PSISM1, TKSoi1, VLHeatCapacity, &
      SoilFracAsMicP, VLWatMicP1, VLiceMicP1 !need the only as some vars are double defined
  use EcoSIMSolverPar, only : NPH
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: RunEcoSIMSurfaceBalance
  contains

  subroutine RunEcoSIMSurfaceBalance(NYS)
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use SurfPhysMod       , only : RunSurfacePhysModel, StageSurfacePhysModel
  use StartsMod         , only : set_ecosim_solver
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS, I, J, M, heat_vec_size
  integer, intent(in) :: NYS
  real(r8) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)
  real(r8) :: ResistanceLitRLay(JY,JX)
  real(r8) :: KSatReductByRainKineticEnergyS(JY,JX)
  real(r8) :: HeatFlux2Ground(JY,JX)
  real(r8) :: TopLayWatVol(JY,JX)

  NHW=1;NHE=1;NVN=1;NVS=NYS

  write(*,*) "In Run Surface Balance"

  call SetMeshATS(NHW,NVN,NHE,NVS)

  !call set_ecosim_solver(1, 1, 1, 1)
  
  NX=1

  write(*,*) "Checking dimensions:"
  write(*, *) "Dimensions of TKSoi1:", SIZE(TKSoi1, 1), SIZE(TKSoi1, 2), SIZE(TKSoi1, 3)
  write(*, *) "Dimensions of a_TEMP:", SIZE(a_TEMP, 1), SIZE(a_TEMP, 2)

  write(*,*) "Starting loop: "

  do NY=1,NYS
    NU(NY,NX)=a_NU(NY)
    NL(NY,NX)=a_NL(NY)
    AREA(3,0,NY,NX)=a_AREA3(NY)
    AREA(3,NU(NY,NX),NY,NX)=a_AREA3(NY)
    ASP(NY,NX)=a_ASP(NY)
    !TairKClimMean(NY,NX)=a_ATKA(NY)
    !CO2E(NY,NX)=atm_co2
    !CH4E(NY,NX)=atm_ch4
    !OXYE(NY,NX)=atm_o2
    !Z2GE(NY,NX)=atm_n2
    !Z2OE(NY,NX)=atm_n2o
    !ZNH3E(NY,NX)=atm_nh3
    !H2GE(NY,NX)=atm_H2
    TairK(NY,NX)=tairc(NY)
    VPA(NY,NX) = vpair(NY)
    WindSpeedAtm(NY,NX) = uwind(NY)
    SWRadOnGrnd(NY,NX) = swrad(NY)
    LWRadSky(NY,NX) = sunrad(NY)
    !RainH(NY,NX) = prec
    DO L=NU(NY,NX),NL(NY,NX)
      !FieldCapacity(L,NY,NX)=a_FC(L,ny)
      !WiltPoint(L,NY,NX)=a_WP(L,NY)
      CumDepth2LayerBottom(L,NY,NX)=a_CumDepth2LayerBottom(L,NY)
      SoiBulkDensityt0(L,NY,NX)=a_BKDSI(L,NY)
      CORGC(L,NY,NX)=a_CORGC(L,NY)
      CORGN(L,NY,NX)=a_CORGN(L,NY)
      CORGP(L,NY,NX)=a_CORGP(L,NY)
      VLWatMicP1(L,NY,NX)=a_WC(L,NY)
      VLiceMicP1(L,NY,NX)=0.0
      TKSoi1(L,NY,NX) = a_TEMP(L,NY)
      VLHeatCapacity(L,NY,NX) = heat_capacity
      SoilFracAsMicP(L,NY,NX) = 1.0
      PSISM1(L,NY,NX) = a_MATP(L,NY)
    ENDDO
  ENDDO

  PSIAtFldCapacity = pressure_at_field_capacity
  PSIAtWiltPoint = pressure_at_wilting_point

  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  DO M=1,NPH
    call RunSurfacePhysModel(M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
      KSatReductByRainKineticEnergyS,HeatFlux2Ground,TopLayWatVol)
  ENDDO

  do NY=1,NYS
    write(*,*) "Col: ", NY , "Heat Flux: ", HeatFlux2Ground(NY,1)
  ENDDO

  DO NY=1,NYS
    !for every column send the top layer to the transfer var
    surf_e_source(NY) = HeatFlux2Ground(NY,1) / (dts_HeatWatTP*3600._r8)
  ENDDO

  write(*,*) "Finished Subroutine RunEcoSIMSurfaceBalance"
  end subroutine RunEcoSIMSurfaceBalance

end module ATSEcoSIMAdvanceMod
