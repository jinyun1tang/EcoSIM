module PlantMngmtDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__



  real(r8) :: THIN(05,366,JY,JX)                !thinning of plant population, [-]
  real(r8) :: EHVST(2,4,05,366,JY,JX)           !harvest efficiency, [-]
  real(r8) :: HVST(05,366,JY,JX)                !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  integer  :: IYR0(JP,JY,JX)                    !year of planting, [-]
  integer  :: IYRH(JP,JY,JX)                    !year of harvest, [-]
  integer  :: IDAY0(JP,JY,JX)                   !day of planting, [-]
  integer  :: IDAYH(JP,JY,JX)                   !day of harvest, [-]
  integer  :: IHVST(05,366,JY,JX)               !type of harvest, [-]
  integer  :: JHVST(05,366,JY,JX)               !flag for stand replacing disturbance, [-]
  integer  :: IDTH(JP,JY,JX)                    !flag for species death, [-]
  integer  :: IYRX(JP,JY,JX)                    !alternate year of planting, [-]
  integer  :: IDAYX(JP,JY,JX)                   !alternate day of planting, [-]
  integer  :: IYRY(JP,JY,JX)                    !alternate year of harvest, [-]
  integer  :: IDAYY(JP,JY,JX)                   !alternate day of harvest, [-]

  real(r8) :: UCO2F(JY,JX)                      !total CO2 flux from fire, [g d-2]
  real(r8) :: UCH4F(JY,JX)                      !total CH4 flux from fire, [g d-2]
  real(r8) :: UOXYF(JY,JX)                      !total O2 flux from fire, [g d-2]
  real(r8) :: UNH3F(JY,JX)                      !total NH3 flux from fire, [g d-2]
  real(r8) :: UN2OF(JY,JX)                      !total N2O flux from fire, [g d-2]
  real(r8) :: UPO4F(JY,JX)                      !total PO4 flux from fire, [g d-2]

end module PlantMngmtDataType
