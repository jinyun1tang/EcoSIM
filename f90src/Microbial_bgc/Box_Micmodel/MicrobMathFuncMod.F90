module MicrobMathFuncMod
  use EcosimConst
  use data_kind_mod, only : r8 => DAT_KIND_R8  
implicit none

  character(len=*), parameter, private :: mod_filename=&
  __FILE__

  real(r8), parameter, private :: ZERO=1.0E-15_r8
  contains
!------------------------------------------------------------------------

  subroutine MicrobPhysTempFun(TKSO, TSensGrowth, TSensMaintR)
  !
  !the physiological temperature dependence of microbes
  implicit none
  real(r8), intent(in) :: TKSO
  real(r8), intent(out):: TSensGrowth !temperature sensitivity for growth respiration
  real(r8), intent(out):: TSensMaintR !temperature sensitivity for maintenance respiration
  real(r8) :: RTK,STK
  real(r8) :: ACTV,ACTVM

  RTK=RGASC*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TSensGrowth=EXP(25.229_r8-62500._r8/RTK)/ACTV
  ACTVM=1+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TSensMaintR=EXP(25.214_r8-62500._r8/RTK)/ACTVM

  end subroutine MicrobPhysTempFun
!------------------------------------------------------------------------


  pure function TranspBasedsubstrateUptake(S_conc,diffusc, KM, V_max, zeros)result(uptake)
  !
  !transport based substrate uptake
  implicit none
  real(r8), intent(in) :: S_conc   !substrate concentration
  real(r8), intent(in) :: diffusc !diffusion coefficient
  real(r8), intent(in) :: KM      !half saturation parameter
  real(r8), intent(in) :: V_max   !maximum uptake rate
  real(r8), optional, intent(in) :: zeros !threshold for active uptake
  real(r8) :: uptake
  real(r8) :: X,B,C
  real(r8) :: zero1

  if(present(zeros))then
    zero1=zeros
  else
    zero1=zero
  endif  

  !obtain uptake flux
  X=diffusc*S_conc
  IF(X.GT.ZERO1)THEN
    B=-V_max-diffusc*KM-X
    C=X*V_max
    uptake=(-B-SQRT(B*B-4.0_r8*C))/2.0_r8
  ELSE
    uptake=0.0_r8
  ENDIF
  end function TranspBasedsubstrateUptake
!------------------------------------------------------------------------

end module MicrobMathFuncMod
