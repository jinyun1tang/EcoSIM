program history_output_precision_probe
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: ieee_exceptions
  use HistoryOutputPrecisionMod
  use netcdf
  implicit none
  type(history_output_environment) :: environment
  integer :: nc, dx, dy, v1, v2, vd
  real(real64), volatile :: values(6), values2(6,1), recovered(6), double_values(6)
  real(real32), volatile :: physical_conversion
  logical :: halt, flag, fatal_modes(3), after_modes(3), expected_halt, expected_flag
  character(len=32) :: mode
  call get_command_argument(1, mode)
  values = [2.5006010422333257e-49_real64, -2.5006010422333257e-49_real64, &
            1.e-40_real64, 3.2344201e-6_real64, 1.5_real64, -2.25_real64]
  values2(:,1) = values
  call ieee_get_halting_mode(ieee_invalid, fatal_modes(1))
  call ieee_get_halting_mode(ieee_divide_by_zero, fatal_modes(2))
  call ieee_get_halting_mode(ieee_overflow, fatal_modes(3))
  if (.not. all(fatal_modes)) error stop 'Fatal traps were not enabled by compiler flags'
  call ieee_get_halting_mode(ieee_underflow, halt)
  if (.not. halt) error stop 'Underflow trap was not enabled'
  expected_halt = .true.
  expected_flag = .false.
  if (trim(mode) == 'existing-flag') then
    expected_halt = .false.
    expected_flag = .true.
    call ieee_set_halting_mode(ieee_underflow, expected_halt)
    call ieee_set_flag(ieee_underflow, expected_flag)
  end if
  call check(nf90_create('precision_probe.nc', nf90_clobber, nc))
  call check(nf90_def_dim(nc, 'x', 6, dx))
  call check(nf90_def_dim(nc, 'y', 1, dy))
  call check(nf90_def_var(nc, 'float_rank1', nf90_float, [dx], v1))
  call check(nf90_def_var(nc, 'float_rank2', nf90_float, [dx,dy], v2))
  call check(nf90_def_var(nc, 'double_rank1', nf90_double, [dx], vd))
  call check(nf90_enddef(nc))
  if (trim(mode) == 'existing-flag') call ieee_set_flag(ieee_underflow, .true.)
  if (trim(mode) /= 'unguarded') call begin_history_output(environment, .true.)
  if (trim(mode) == 'overflow') values(1) = 1.e100_real64
  call check(nf90_put_var(nc, v1, values))
  call check(nf90_put_var(nc, v2, values2))
  if (trim(mode) /= 'unguarded') call end_history_output(environment)
  call ieee_get_halting_mode(ieee_underflow, halt)
  call ieee_get_flag(ieee_underflow, flag)
  if ((halt .neqv. expected_halt) .or. (flag .neqv. expected_flag)) &
    error stop 'Underflow environment was not restored'
  call ieee_get_halting_mode(ieee_invalid, after_modes(1))
  call ieee_get_halting_mode(ieee_divide_by_zero, after_modes(2))
  call ieee_get_halting_mode(ieee_overflow, after_modes(3))
  if (any(after_modes .neqv. fatal_modes)) error stop 'Fatal trap modes changed'
  if (trim(mode) == 'physics') then
    physical_conversion = real(values(1), real32)
    print *, physical_conversion
    error stop 'Physics underflow did not trap'
  end if
  call begin_history_output(environment, .false.)
  call check(nf90_put_var(nc, vd, values))
  call end_history_output(environment)
  call check(nf90_get_var(nc, v1, recovered))
  if (any(recovered(1:2) /= 0._real64)) error stop 'Tiny float32 values did not round to zero'
  if (recovered(3) <= 0._real64) error stop 'Representable subnormal was flushed to zero'
  if (any(recovered(5:6) /= values(5:6))) error stop 'Ordinary float32 values changed'
  call check(nf90_get_var(nc, v2, values2))
  if (any(values2(:,1) /= recovered)) error stop 'Rank-two writer differs'
  call check(nf90_get_var(nc, vd, double_values))
  if (any(double_values /= values)) error stop 'Double history values changed'
  call check(nf90_close(nc))
  ! The harness began with compiler underflow traps enabled. On some targets
  ! main's IEEE finalization restores that initial trap mode; do not carry our
  ! deliberately injected test flag into finalization. Its preservation has
  ! already been asserted immediately after the guarded write above.
  if (trim(mode) == 'existing-flag') call ieee_set_flag(ieee_underflow, .false.)
  print *, 'PASS: rank1/rank2 float output, double values and IEEE traps preserved'
contains
  subroutine check(code)
    integer, intent(in) :: code
    if (code /= nf90_noerr) then
      print *, trim(nf90_strerror(code))
      error stop 'NetCDF probe failed'
    end if
  end subroutine check
end program history_output_precision_probe
