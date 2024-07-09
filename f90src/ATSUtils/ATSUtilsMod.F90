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

  integer :: L,NTG

  DO L =0,NL(NY,NX)
    DO NTG=idg_beg,idg_end-1
      call PhasePartition(VLsoiAirP_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX),GasSolbility_vr(NTG,L,NY,NX),&
        trc_solcl_vr(NTG,L,NY,NX),trc_gascl_vr(NTG,L,NY,NX))
    ENDDO
  ENDDO

  end subroutine ForceGasAquaEquil

end module ATSUtilsMod
