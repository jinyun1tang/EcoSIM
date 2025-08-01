module ATSUtilsMod
  use ChemTranspDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilBGCDataType
  use GridConsts
  use GridDataType
  use SoilWaterDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: ForceGasAquaEquil, ComputeDatefromATS
  ! Add this interface block
  contains

  subroutine PhasePartition(VG,VW,SL,CW,CG)
!
! Let Cg and CW be the gaseous and aqueous concentration, then
! Ct=CG*VG+CW*VW=CG*VG+CG*SL*VW=CG*(VG+SL*VW)
! VW=CG*SL
! VG: VOLP
! VW: VOLW
! SL: solubility
  implicit none
  real(r8), intent(in) :: VG,VW,SL
  real(r8), intent(inout) :: CW,CG
  real(r8) :: CT,VT

  CT=CG*VG+CW*VW
  VT=VG+SL*VW
  if(VT>0._r8)then
    CG=CT/VT
  ELSE
    CG=0._r8
  ENDIF
  CW=CG*SL

  end subroutine PhasePartition
!------------------------------------------------------------------------------------------
  subroutine ForceGasAquaEquil(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,NTG

  DO L =0,NL_col(NY,NX)
    DO NTG=idg_beg,idg_NH3
      call PhasePartition(VLsoiAirP_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX),GasSolbility_vr(NTG,L,NY,NX),&
        trc_solcl_vr(NTG,L,NY,NX),trcg_gascl_vr(NTG,L,NY,NX))
    ENDDO
  ENDDO

  end subroutine ForceGasAquaEquil

!------------------------------------------------------------------------------------------
  subroutine ComputeDatefromATS(current_day, current_year, current_month, day_of_month, total_days_in_month)
    !This takes the current day (0-364) and current year
    !and returns the current month, day of that month and total days in that month
    implicit none
    integer, intent(in) :: current_day, current_year
    integer, intent(out) :: current_month, day_of_month, total_days_in_month

    integer, parameter :: days_in_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    integer :: cumulative_days, i, days_this_month
    logical :: is_leap

    ! Check if leap year
    !is_leap = (mod(current_year, 4) == 0 .and. mod(current_year, 100) /= 0) .or. (mod(current_year, 400) == 0)

    ! Find month and day
    cumulative_days = 0
    do i = 1, 12
      days_this_month = days_in_month(i)
      if (is_leap .and. i == 2) days_this_month = 29

      if (current_day < cumulative_days + days_this_month) then
        current_month = i
        day_of_month = current_day - cumulative_days + 1
        total_days_in_month = days_this_month
        exit
      endif
      cumulative_days = cumulative_days + days_this_month
    enddo

  end subroutine ComputeDatefromATS
end module ATSUtilsMod
