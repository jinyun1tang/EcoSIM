module PlantDataRateType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: RDFOMC(2,0:4,JZ,JP,JY,JX)         !root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8) :: RDFOMN(2,0:4,JZ,JP,JY,JX)         !root uptake (+ve) - exudation (-ve) of DON, [g d-2 h-1]
  real(r8) :: RDFOMP(2,0:4,JZ,JP,JY,JX)         !root uptake (+ve) - exudation (-ve) of DOP, [g d-2 h-1]
  real(r8) :: RUPNH4(2,JZ,JP,JY,JX)             !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8) :: RUPNHB(2,JZ,JP,JY,JX)             !root uptake of NH4 band, [g d-2 h-1]
  real(r8) :: RUPNO3(2,JZ,JP,JY,JX)             !root uptake of NO3 non-band, [g d-2 h-1]
  real(r8) :: RUPNOB(2,JZ,JP,JY,JX)             !root uptake of NO3 band, [g d-2 h-1]
  real(r8) :: RUPH2P(2,JZ,JP,JY,JX)             !root uptake of PO4 non-band, [g d-2 h-1]
  real(r8) :: RUPH2B(2,JZ,JP,JY,JX)             !root uptake of PO4 band, [g d-2 h-1]
  real(r8) :: RUPNF(JZ,JP,JY,JX)                !root N2 fixation, [g d-2 h-1]
  real(r8) :: RUPHGS(2,JZ,JP,JY,JX)             !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8) :: RUPH1P(2,JZ,JP,JY,JX)
  real(r8) :: RUPH1B(2,JZ,JP,JY,JX)
  real(r8) :: RUPP2P(2,JZ,JP,JY,JX)
  real(r8) :: RUPP2B(2,JZ,JP,JY,JX)
  real(r8) :: RUPP1P(2,JZ,JP,JY,JX)
  real(r8) :: RUPP1B(2,JZ,JP,JY,JX)
  real(r8) :: RCZLX(JC,JP,JY,JX)                !N translocated from leaf during senescence, [g d-2 h-1]
  real(r8) :: RCPLX(JC,JP,JY,JX)                !P translocated from leaf during senescence, [g d-2 h-1]
  real(r8) :: RCCLX(JC,JP,JY,JX)                !C translocated from leaf during senescence, [g d-2 h-1]
  real(r8) :: RCZSX(JC,JP,JY,JX)                !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8) :: RCPSX(JC,JP,JY,JX)                !P translocated from sheath during senescence, [g d-2 h-1]
  real(r8) :: RCCSX(JC,JP,JY,JX)                !C translocated from sheath during senescence, [g d-2 h-1]

  real(r8) :: CARBN(JP,JY,JX)                   !total gross CO2 fixation, [g d-2 ]
  real(r8) :: TCSNC(JP,JY,JX)                   !total plant C litterfall , [g d-2 ]
  real(r8) :: TZSNC(JP,JY,JX)                   !total plant N litterfall , [g d-2 ]
  real(r8) :: TPSNC(JP,JY,JX)                   !total plant P litterfall , [g d-2 ]
  real(r8) :: TZUPFX(JP,JY,JX)                  !total plant N2 fixation, [g d-2 ]
  real(r8) :: TCO2T(JP,JY,JX)                   !total plant respiration, [g d-2 ]
  real(r8) :: BALC(JP,JY,JX)                    !plant C balance, [g d-2]
  real(r8) :: BALN(JP,JY,JX)                    !plant N balance, [g d-2]
  real(r8) :: BALP(JP,JY,JX)                    !plant P balance, [g d-2]
  real(r8) :: HCSNC(JP,JY,JX)                   !plant C litterfall, [g d-2 h-1]
  real(r8) :: HZSNC(JP,JY,JX)                   !plant N litterfall, [g d-2 h-1]
  real(r8) :: HPSNC(JP,JY,JX)                   !plant C litterfall, [g d-2 h-1]
  real(r8) :: ZNPP(JP,JY,JX)                    !total net primary productivity, [g d-2]
  real(r8) :: CTRAN(JP,JY,JX)                   !total transpiration, [m d-2]

  real(r8) :: HVSTC(JP,JY,JX)                   !plant C harvest, [g d-2 ]
  real(r8) :: HVSTN(JP,JY,JX)                   !plant N harvest, [g d-2 ]
  real(r8) :: HVSTP(JP,JY,JX)                   !plant P harvest, [g d-2 ]
  real(r8) :: THVSTC(JP,JY,JX)                  !total plant C harvest, [g d-2 ]
  real(r8) :: THVSTN(JP,JY,JX)                  !total plant N harvest, [g d-2 ]
  real(r8) :: THVSTP(JP,JY,JX)                  !total plant P harvest, [g d-2 ]
  real(r8) :: TCO2A(JP,JY,JX)                   !total autotrophic respiration, [g d-2 ]
  real(r8) :: VCO2F(JP,JY,JX)                   !plant CO2 emission from fire, [g d-2 ]
  real(r8) :: VCH4F(JP,JY,JX)                   !plant CH4 emission from fire, [g d-2 ]
  real(r8) :: VOXYF(JP,JY,JX)                   !plant O2 uptake from fire, [g d-2 ]
  real(r8) :: VNH3F(JP,JY,JX)                   !plant NH3 emission from fire, [g d-2 ]
  real(r8) :: VN2OF(JP,JY,JX)                   !plant N2O emission from fire, [g d-2 ]
  real(r8) :: VPO4F(JP,JY,JX)                   !plant PO4 emission from fire, [g d-2 ]

end module PlantDataRateType
