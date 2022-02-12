module InitEcoSIM


  implicit none
  private
  public ::   InitModules

  contains


  subroutine InitModules(nmicbguilds)
  use NitroMod, only : InitNitro
  use RedistMod, only : InitRedist
  use MicrobialDataType, only : InitMicrobialData
  implicit  none
  integer, intent(in) :: nmicbguilds

  call InitMicrobialData(nmicbguilds)

  call initNitro

  call InitRedist
  end subroutine InitModules
end module InitEcoSIM
