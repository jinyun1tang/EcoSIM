module CanopyDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save

  real(r8) :: RSMN(JP,JY,JX)                    !canopy minimum stomatal resistance, [s m-1]
  real(r8) :: RAC(JY,JX)                        !canopy boundary layer resistance, [m h-1]
  real(r8) :: TVOLWC(JY,JX)                     !canopy surface water content, [m3 d-2]
  real(r8) :: TFLWCI(JY,JX)                     !net ice transfer to canopy, [MJ d-2 t-1]
  real(r8) :: TFLWC(JY,JX)                      !net water transfer to canopy, [MJ d-2 t-1]
  real(r8) :: RAD1(JP,JY,JX)                    !canopy net radiation , [MJ d-2 h-1]
  real(r8) :: THRM1(JP,JY,JX)                   !canopy longwave radiation , [MJ d-2 h-1]
  real(r8) :: EFLXC(JP,JY,JX)                   !canopy latent heat flux, [MJ d-2 h-1]
  real(r8) :: SFLXC(JP,JY,JX)                   !canopy sensible heat flux, [MJ d-2 h-1]
  real(r8) :: HFLXC(JP,JY,JX)                   !canopy storage heat flux, [MJ d-2 h-1]
  real(r8) :: ENGYX(JP,JY,JX)                   !canopy heat storage from previous time step, [MJ d-2]
  real(r8) :: VHCPC(JP,JY,JX)                   !canopy heat capacity, [MJ d-2 K-1]
  real(r8) :: PSILT(JP,JY,JX)                   !canopy total water potential , [Mpa]
  real(r8) :: PSILG(JP,JY,JX)                   !canopy turgor water potential, [Mpa]
  real(r8) :: PSILO(JP,JY,JX)                   !canopy osmotic water potential, [Mpa]
  real(r8) :: RC(JP,JY,JX)                      !canopy stomatal resistance, [h m-1]
  real(r8) :: RA(JP,JY,JX)                      !canopy boundary layer resistance, [h m-1]
  real(r8) :: EP(JP,JY,JX)                      !canopy transpiration, [m2 d-2 h-1]
  real(r8) :: EVAPC(JP,JY,JX)                   !canopy evaporation, [m2 d-2 h-1]
  real(r8) :: VOLWP(JP,JY,JX)                   !canopy water content, [m3 d-2]
  real(r8) :: RNH3C(JP,JY,JX)                   !canopy NH3 flux, [g d-2 h-1]
  real(r8) :: RSETC(JP,JY,JX)                   !effect of canopy C status on seed set , []
  real(r8) :: RSETN(JP,JY,JX)                   !effect of canopy N status on seed set , []
  real(r8) :: RSETP(JP,JY,JX)                   !effect of canopy P status on seed set , []
  real(r8) :: TNH3C(JP,JY,JX)                   !total canopy NH3 flux, [g d-2 ]
  real(r8) :: TZSN0(JP,JY,JX)                   !total surface litterfall N, [g d-2]
  real(r8) :: TPSN0(JP,JY,JX)                   !total surface litterfall P, [g d-2]
  real(r8) :: TCSN0(JP,JY,JX)                   !total surface litterfall C, [g d-2]
  real(r8) :: ARLFC(JY,JX)                      !total canopy leaf area, [m2 d-2]
  real(r8) :: ARSTC(JY,JX)                      !total canopy stem area, [m2 d-2]
  real(r8) :: TEVAPP(JY,JX)                     !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: TEVAPC(JY,JX)                     !total canopy evaporation, [m3 d-2]
  real(r8) :: TENGYC(JY,JX)                     !total canopy heat content, [MJ  d-2]
  real(r8) :: THFLXC(JY,JX)                     !total canopy heat flux, [MJ  d-2]
  real(r8) :: TVOLWP(JY,JX)                     !total canopy water content, [m3 d-2]
  real(r8) :: THRMC(JY,JX)                      !total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: TCNET(JY,JX)                      !total canopy net CO2 exchange, [g d-2 h-1]
  real(r8) :: ARLSS(JY,JX)                      !stalk area of combined,each PFT canopy
  real(r8) :: WGLFT(JC,JY,JX)                   !total leaf mass, [g d-2]
  real(r8) :: ARLFT(JC,JY,JX)                   !total leaf area, [m2 d-2]
  real(r8) :: ARSTT(JC,JY,JX)                   !total stem area, [m2 d-2]

end module CanopyDataType
