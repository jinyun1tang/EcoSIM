! EcoSIM_ZM: local policy for float32 history serialization, not model arithmetic.
module HistoryOutputPrecisionMod
  use, intrinsic :: ieee_exceptions, only : ieee_underflow, ieee_support_halting, &
       ieee_support_flag, ieee_get_halting_mode, ieee_set_halting_mode, &
       ieee_get_flag, ieee_set_flag
  implicit none
  private
  public :: history_output_environment, begin_history_output, end_history_output

  type :: history_output_environment
    private
    logical :: active = .false.
    logical :: halting = .false.
    logical :: flag = .false.
  end type history_output_environment

contains

  subroutine begin_history_output(environment, single_precision)
    type(history_output_environment), intent(inout) :: environment
    logical, intent(in) :: single_precision
    if (environment%active) error stop 'Nested history output environment'
    if (.not. single_precision) return
    if (.not. ieee_support_halting(ieee_underflow) .or. &
        .not. ieee_support_flag(ieee_underflow)) then
      error stop 'Float32 history output requires IEEE underflow control'
    end if
    call ieee_get_halting_mode(ieee_underflow, environment%halting)
    call ieee_get_flag(ieee_underflow, environment%flag)
    ! A finite double below the float32 range may round to a subnormal or zero
    ! inside NetCDF. Permit that storage conversion, without clipping model
    ! state or changing invalid/divide-by-zero/overflow/inexact trap settings.
    call ieee_set_halting_mode(ieee_underflow, .false.)
    environment%active = .true.
  end subroutine begin_history_output

  subroutine end_history_output(environment)
    type(history_output_environment), intent(inout) :: environment
    if (.not. environment%active) return
    ! Clear the transient output exception before re-enabling its trap, then
    ! restore the incoming flag. Some targets clear a flag when changing its
    ! halting mode, so setting the saved flag first would lose caller state.
    ! No other exception flag or trap setting is changed here.
    call ieee_set_flag(ieee_underflow, .false.)
    call ieee_set_halting_mode(ieee_underflow, environment%halting)
    call ieee_set_flag(ieee_underflow, environment%flag)
    environment%active = .false.
  end subroutine end_history_output

end module HistoryOutputPrecisionMod
