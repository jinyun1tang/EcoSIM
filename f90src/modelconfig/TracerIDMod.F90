module TracerIDMod

  use MiniMathMod, only : addone
implicit none
  save
  CHARACTER(LEN=*), private, PARAMETER :: MOD_FILENAME=__FILE__

  integer, parameter :: idg_CO2=1
  integer, parameter :: idg_CH4=2
  integer, parameter :: idg_O2 =3
  integer, parameter :: idg_N2 =4
  integer, parameter :: idg_N2O=5
  integer, parameter :: idg_H2 =6
  integer, parameter :: idg_NH3=7
  integer  :: idg_NH3B
  integer :: ids_NH4,ids_NH4B
  integer :: ids_NO3,ids_NO3B
  integer :: ids_NO2,ids_NO2B
  integer :: ids_H1PO4,ids_H1PO4B
  integer :: ids_H2PO4,ids_H2PO4B

  integer :: ids_DOC,ids_DON,ids_DOP,ids_ODA
  integer :: idg_beg,idg_end
  integer :: ids_beg,ids_end
  integer :: ids_nut_beg,ids_nuts_beg,ids_nuts_end

! salt tracers
  integer :: idsa_beg,idsa_end

  integer :: idsa_Al      ! Al
  integer :: idsa_Fe      ! Fe
  integer :: idsa_Hp      ! H+
  integer :: idsa_Ca      ! Ca
  integer :: idsa_Mg      ! Mg
  integer :: idsa_Na      ! Na
  integer :: idsa_K       ! K
  integer :: idsa_OH      ! OH
  integer :: idsa_SO4     ! SO4
  integer :: idsa_Cl      ! Cl
  integer :: idsa_CO3     ! CO3,
  integer :: idsa_HCO3    ! HCO3
  integer :: idsa_AlOH    ! AlOH
  integer :: idsa_AlOH2   ! AlOH2
  integer :: idsa_AlOH3   ! AlOH3
  integer :: idsa_AlOH4   ! AlOH4
  integer :: idsa_AlSO4   ! AlSO4
  integer :: idsa_FeOH    ! FeOH
  integer :: idsa_FeOH2   ! FeOH2
  integer :: idsa_FeOH3   ! FeOH3
  integer :: idsa_FeOH4   ! FeOH4
  integer :: idsa_FeSO4   ! FeSO4
  integer :: idsa_CaOH2   ! CaOH2
  integer :: idsa_CaCO3   ! CaCO3
  integer :: idsa_CaHCO3  ! CaHCO3
  integer :: idsa_CaSO4   ! CaSO4
  integer :: idsa_MgOH2   ! MgOH2
  integer :: idsa_MgCO3   ! MgCO3
  integer :: idsa_MgHCO3  ! MgHCO3
  integer :: idsa_MgSO4   ! MgSO4
  integer :: idsa_NaCO3   ! NaCO3
  integer :: idsa_NaSO4   ! NaSO4
  integer :: idsa_KSO4    ! KSO4
  integer :: idsa_H0PO4   ! PO4
  integer :: idsa_H3PO4   ! H3PO4
  integer :: idsa_FeHPO4  ! FeHPO4
  integer :: idsa_FeH2PO4 ! FeH2PO4
  integer :: idsa_CaPO4   ! CaPO4
  integer :: idsa_CaHPO4  ! CaHPO4
  integer :: idsa_CaH2PO4 ! CaH2PO4
  integer :: idsa_MgHPO4 ! MgHPO4
!band salt
  integer :: idsa_H0PO4B   ! PO4
  integer :: idsa_H3PO4B   ! H3PO4
  integer :: idsa_FeHPO4B  ! FeHPO4
  integer :: idsa_FeH2PO4B ! FeH2PO4
  integer :: idsa_CaPO4B   ! CaPO4
  integer :: idsa_CaHPO4B  ! CaHPO4
  integer :: idsa_CaH2PO4B ! CaH2PO4
  integer :: idsa_MgHPO4B ! MgHPO4
  integer :: idsab_beg
  integer :: idsab_end

  contains

  subroutine InitTracerIDs(lsalt_model)
  implicit none
  logical, intent(in) :: lsalt_model

  idg_beg=1; ids_beg=1
! for better array manipulation of land-atmosphere exchange,
! banded NH3 is considered as (potential) gas too. However, there
! is no banded gas NH3 concentration.

  idg_NH3B=idg_NH3+1
  idg_end=idg_NH3B;


  ids_nuts_beg=idg_NH3;  !the first nutrient tracer, including band
  ids_end=idg_end   !initalize the solute counter
  ids_NH4B=addone(ids_end);ids_NO3B=addone(ids_end);
  ids_NO2B=addone(ids_end);ids_H1PO4B=addone(ids_end)
  ids_H2PO4B=addone(ids_end)

  ids_NH4=addone(ids_end);
  ids_NO3=addone(ids_end);
  ids_NO2=addone(ids_end);
  ids_H1PO4=addone(ids_end);
  ids_H2PO4=addone(ids_end);

  ids_nut_beg=ids_NH4;  !the first non-band non-gaseous nutrient tracer
  ids_nuts_end=ids_H2PO4;!the last non-band nutrient tracer
  if(lsalt_model)then
    idsa_beg=1;idsa_end=0
    idsa_Al=addone(idsa_end)
    idsa_Fe=addone(idsa_end)
    idsa_Hp=addone(idsa_end)
    idsa_Ca=addone(idsa_end)
    idsa_Mg=addone(idsa_end)
    idsa_Na=addone(idsa_end)
    idsa_K=addone(idsa_end)
    idsa_OH=addone(idsa_end)
    idsa_SO4=addone(idsa_end)
    idsa_Cl=addone(idsa_end)
    idsa_CO3=addone(idsa_end)
    idsa_HCO3=addone(idsa_end)
    idsa_AlOH=addone(idsa_end)
    idsa_AlOH2=addone(idsa_end)
    idsa_AlOH3=addone(idsa_end)
    idsa_AlOH4=addone(idsa_end)
    idsa_AlSO4=addone(idsa_end)
    idsa_FeOH =addone(idsa_end)
    idsa_FeOH2=addone(idsa_end)
    idsa_FeOH3=addone(idsa_end)
    idsa_FeOH4=addone(idsa_end)
    idsa_FeSO4=addone(idsa_end)
    idsa_CaOH2=addone(idsa_end)
    idsa_CaCO3=addone(idsa_end)
    idsa_CaHCO3=addone(idsa_end)
    idsa_CaSO4 =addone(idsa_end)
    idsa_MgOH2 =addone(idsa_end)
    idsa_MgCO3 =addone(idsa_end)
    idsa_MgHCO3=addone(idsa_end)
    idsa_MgSO4 =addone(idsa_end)
    idsa_NaCO3 =addone(idsa_end)
    idsa_NaSO4 =addone(idsa_end)
    idsa_KSO4  =addone(idsa_end)
    idsa_H0PO4 =addone(idsa_end)
    idsa_H3PO4 =addone(idsa_end)
    idsa_FeHPO4=addone(idsa_end)
    idsa_FeH2PO4=addone(idsa_end)
    idsa_CaPO4 =addone(idsa_end)
    idsa_CaHPO4=addone(idsa_end)
    idsa_CaH2PO4=addone(idsa_end)
    idsa_MgHPO4=addone(idsa_end)
    idsa_end=idsa_MgHPO4
    idsab_beg=idsa_end+1
    idsab_end=idsab_beg
    idsa_H0PO4B=addone(idsab_end)   ! PO4
    idsa_H3PO4B=addone(idsab_end)  ! H3PO4
    idsa_FeHPO4B=addone(idsab_end)  ! FeHPO4
    idsa_FeH2PO4B=addone(idsab_end) ! FeH2PO4
    idsa_CaPO4B=addone(idsab_end)   ! CaPO4
    idsa_CaHPO4B=addone(idsab_end)  ! CaHPO4
    idsa_CaH2PO4B=addone(idsab_end) ! CaH2PO4
    idsa_MgHPO4B =addone(idsab_end)! MgHPO4
    idsab_end=idsa_MgHPO4B
  else

  endif
  end subroutine InitTracerIDs
end module TracerIDMod
