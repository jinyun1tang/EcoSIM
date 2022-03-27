module EcoSimSumDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8) :: TLH2G    !total soil N2	g d-2
  real(r8) :: H2GIN    !total soil H2 g d-2
  real(r8) :: H2GOU    !cumulative H2 loss through lateral and lower boundaries
  real(r8) :: TION     !total soil ion content	mol d-2
  real(r8) :: TIONIN   !total surface ion flux	mol d-2
  real(r8) :: TIONOU   !total subsurface ion flux	mol d-2
  real(r8) :: TSEDSO   !total soil sediment	Mg d-2
  real(r8) :: TSEDOU   !total sediment subsurface flux	Mg d-2
  real(r8) :: VOLWSO   !total soil water content	m3 d-2
  real(r8) :: HEATSO   !total soil heat content	MJ d-2
  real(r8) :: OXYGSO   !total soil O2 content	g d-2
  real(r8) :: TLRSDC   !total soil litter C content	g d-2
  real(r8) :: TLRSDN   !total soil Litter N content	g d-2
  real(r8) :: TLRSDP   !total soil Litter P content	g d-2
  real(r8) :: TLORGC   !total soil POM + humus C content	g d-2
  real(r8) :: TLORGN   !total soil POM + humus N content	g d-2
  real(r8) :: TLORGP   !total soil POM + humus P content	g d-2
  real(r8) :: TLNH4    !total soil NH4 content	g d-2
  real(r8) :: TLNO3    !total soil NO3 content	g d-2
  real(r8) :: TLPO4    !total soil PO4 content	g d-2
  real(r8) :: TBALC    !total plant C balance	g d-2
  real(r8) :: TBALN    !total plant N balance	g d-2
  real(r8) :: TBALP    !total plant P balance	g d-2
  real(r8) :: CRAIN    !total precipitation	m3 d-2
  real(r8) :: HEATIN   !total surface heat flux	MJ d-2
  real(r8) :: OXYGIN   !total surface O2 flux	g d-2
  real(r8) :: TORGF    !total organic C amendment	g d-2
  real(r8) :: TORGN    !total organic N amendment	g d-2
  real(r8) :: TORGP    !total organic P amendment	g d-2
  real(r8) :: CO2GIN   !total surface CO2 flux	g d-2
  real(r8) :: ZN2GIN   !total surface N2 flux	g d-2
  real(r8) :: VOLWOU   !total subsurface water flux	m3 d-2
  real(r8) :: CEVAP    !total evaporation	m3 d-2
  real(r8) :: CRUN     !total surface runoff	m3 d-2
  real(r8) :: HEATOU   !total subsurface heat flux	MJ d-2
  real(r8) :: OXYGOU   !total subsurface O2 flux	g d-2
  real(r8) :: TCOU     !total subsurface C flux	g d-2
  real(r8) :: TZOU     !total subsurface N flux	g d-2
  real(r8) :: TPOU     !total subsurface P flux	g d-2
  real(r8) :: TZIN     !total surface N flux	g d-2
  real(r8) :: TPIN     !total surface P flux	g d-2
  real(r8) :: XCSN     !total litterfall C	g d-2
  real(r8) :: XZSN     !total litterfall N	g d-2
  real(r8) :: XPSN     !total litterfall P	g d-2
  real(r8) :: TLCO2G   !total soil CO2	g d-2
  real(r8) :: TLN2G    !total soil N2	g d-2

end module EcoSimSumDataType
