module ATSUtilsMod
  use ChemTranspDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use SoilBGCDataType
  use GridConsts
  use GridDataType
  use SoilWaterDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__
  public :: ForceGasAquaEquil
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

  integer :: L

  DO L =0,NL(NY,NX)

    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SCO2L(L,NY,NX),CCO2S(L,NY,NX),CCO2G(L,NY,NX))
    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SCH4L(L,NY,NX),CCH4S(L,NY,NX),CCH4G(L,NY,NX))
    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SOXYL(L,NY,NX),COXYS(L,NY,NX),COXYG(L,NY,NX))
    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SN2GL(L,NY,NX),CZ2GS(L,NY,NX),CZ2GG(L,NY,NX))
    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SN2OL(L,NY,NX),CZ2OS(L,NY,NX),CZ2OG(L,NY,NX))
    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SNH3L(L,NY,NX),CNH3S(L,NY,NX),CNH3G(L,NY,NX))
    call PhasePartition(VOLP(L,NY,NX),VOLW(L,NY,NX),SH2GL(L,NY,NX),CH2GS(L,NY,NX),CH2GG(L,NY,NX))

  ENDDO

  end subroutine ForceGasAquaEquil

end module ATSUtilsMod
