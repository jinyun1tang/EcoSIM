module CyanoBacteriaMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use abortutils,           only: endrun,   destroy
  use EcoSIMCtrlMod,        only: etimer
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: DebugPrint,PrintInfo
  use minimathmod            
  use ElmIDMod
  use TracerIDMod
  use MicAutoCPLXMod
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod
  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains
!------------------------------------------------------------------------------------------
  subroutine CyanoBacteriaCatabolism(I,J,N,RMOMK,TOMEAutoKC,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !

  implicit none
  integer, intent(in) :: I,J, N
  real(r8), intent(in) :: RMOMK(2)        
  real(r8), intent(in) :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag    
  character(len=*), parameter :: subname='CyanoBacteriaCatabolism'

  call PrintInfo('beg '//subname)

  !photosynthesis
  !CO2 + H2O + light -> CH2O + O2
  ! PhotoCyano = Pmax_cyano * BiomassC * fPAR * fMoist * fTemp * fPH

  call PrintInfo('end '//subname)

  end subroutine CyanoBacteriaCatabolism

  end module CyanoBacteriaMod