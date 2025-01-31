module SoilDiagsMod
  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use CanopyDataType,     only: QVegET_col
  use GridDataType,       only: NU, NL
  use EcoSimConst,        only: RGasC
  use abortutils,         only: endrun
  use TracerPropMod,      only: MolecularWeight
  use ChemTranspDataType, only: GasSolbility_vr
  use SoilBGCDataType
  use GridDataType
  use SurfLitterDataType
  use CanopyDataType
  use SnowDataType
  use BalanceCheckDataType
  use SoilWaterDataType
  use ClimForcDataType
  use SurfSoilDataType
  use EcosimBGCFluxType
  use SoilHeatDataType
  use PlantDataRateType
  use EcoSimSumDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: DiagSoilGasPressure
  contains

  subroutine DiagSoilGasPressure(I,J,NHW,NHE,NVN,NVS)  

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: idg, NY,NX, L
  real(r8) :: GasMassSolubility(idg_beg:idg_end)    !g m3 /(mol d-2)
  real(r8) :: GasPres(idg_beg:idg_end)              !gas pressure for each tracer [Pa]

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)    
        do idg=idg_beg,idg_NH3
          GasMassSolubility(idg) =MolecularWeight(idg)*GasSolbility_vr(idg,L,NY,NX)*VLWatMicP_vr(L,NY,NX)  !conver into carbon g C/mol
        enddo  
        GasMassSolubility(idg_NH3B)=GasMassSolubility(idg_NH3)

        Soil_Gas_pressure_vr(L,NY,NX)=0._r8
        do idg=idg_beg,idg_end
          GasPres(idg)=trcs_solml_vr(idg,L,NY,NX)*RGasC*TKS_vr(L,NY,NX)/GasMassSolubility(idg)
          Soil_Gas_pressure_vr(L,NY,NX)=Soil_Gas_pressure_vr(L,NY,NX)+GasPres(idg)
        enddo
        CO2_Gas_Frac_vr(L,NY,NX) = GasPres(idg_CO2)/Soil_Gas_pressure_vr(L,NY,NX)*1.E6
        CH4_Gas_Frac_vr(L,NY,NX) = GasPres(idg_CH4)/Soil_Gas_pressure_vr(L,NY,NX)*1.E6
        Ar_Gas_Frac_vr(L,NY,NX) = GasPres(idg_Ar)/Soil_Gas_pressure_vr(L,NY,NX)*1.E6
      ENDDO
    ENDDO
  ENDDO    
  end subroutine DiagSoilGasPressure
end module SoilDiagsMod
