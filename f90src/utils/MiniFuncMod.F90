module MiniFuncMod

  use data_kind_mod     , only : r8 => SHR_KIND_R8
implicit none

  character(len=*),private, parameter :: mod_filename = __FILE__
  contains
!------------------------------------------------------------------------------------------

  function FilmThickness(PSISM, is_top_layer)result(ans)

  implicit none
  real(r8), intent(in) :: PSISM
  logical, optional, intent(in) :: is_top_layer
  real(r8) :: ans

  if(PRESENT(is_top_layer) .and. is_top_layer)then
    ans=AMAX1(1.0E-06_r8,EXP(-13.650_r8-0.857_r8*LOG(-PSISM)))
  else
    ans=AMAX1(1.0E-06_r8,EXP(-13.833_r8-0.857_r8*LOG(-PSISM)))
  endif
  end function FilmThickness
!------------------------------------------------------------------------------------------


  function TEFAQUDIF(TK)result(ans)
!
! temperautre effect on aqueous diffusivity
  implicit none
  real(r8), intent(in) :: TK !temperature, [Kelvin]
  real(r8) :: ans

  ans =(TK/298.15_r8)**6._r8

  end function TEFAQUDIF
!------------------------------------------------------------------------------------------


  function TEFGASDIF(TK)result(ans)
!
! temperautre effect on gaseous diffusivity
  implicit none
  real(r8), intent(in) :: TK !temperature, [Kelvin]
  real(r8) :: ans

  ans =(TK/298.15_r8)**1.75_r8

  end function TEFGASDIF
!------------------------------------------------------------------------------------------

  function TortMicporew(THETWT)result(ans)
!
! aqueous tortuosity in micropore
  implicit none
  real(r8), intent(in) :: THETWT     !relative saturation of micropore
  real(r8) :: ans
  ans =0.7_r8*THETWT**2._r8
  end function TortMicporew
!------------------------------------------------------------------------------------------
  function TortMacporew(THETWH)result(ans)
  implicit none
  real(r8), intent(in) :: THETWH   !realative saturation of macropore

  real(r8) :: ans

  ans=AMIN1(1.0_r8,2.8_r8*THETWH**3_r8)

  end function TortMacporew

!------------------------------------------------------------------------------------------
  function fDFGS(scalar,THETWA,Z3SR,is_litter)result(ans)
!
! compute coefficient for air-water gas transfer
  implicit none
  real(r8), intent(in) :: scalar               !rate scalar including temperature effect
  real(r8), intent(in) :: THETWA               !relative saturation
  real(r8), intent(in) :: Z3SR                 !saturation threshold
  logical, optional, intent(in) :: is_litter   !indicator if it is litter layer
  real(r8) :: ans
  real(r8), parameter :: Z1S=0.010_r8
  real(r8), parameter :: Z2SW=12.0_r8
  real(r8), parameter :: Z2SD=12.0_r8
  real(r8), parameter :: Z3SX=0.50_r8
  real(r8), parameter :: Z1R=0.01_r8
  real(r8), parameter :: Z2RW=12.0_r8
  real(r8), parameter :: Z2RD=12.0_r8
  real(r8), parameter :: Z3R=0.50_r8
  real(r8) :: Z3S

! Z1S,Z2SW,Z2SD,Z3SX=parameters for air-water gas transfers in soil
! Z1R,Z2RW,Z2RD,Z3RX=parameters for air-water gas transfers in litter

  if(present(is_litter) .and. is_litter)then
    IF(THETWA.GT.Z3R)THEN
      ans=AMAX1(0.0_r8,scalar/((Z1R**(-1._r8))*EXP(Z2RW*(THETWA-Z3R))))
    ELSE
      ans=AMIN1(1.0_r8,scalar/((Z1R**(-1._r8))*EXP(Z2RD*(THETWA-Z3R))))
    ENDIF
  else
    Z3S=AMAX1(Z3SX,Z3SR)
    IF(THETWA.GT.Z3S)THEN
      ans=AMAX1(0.0_r8,scalar/((Z1S**(-1._r8))*EXP(Z2SW*(THETWA-Z3S))))
    ELSE
      ans=AMIN1(1.0_r8,scalar/((Z1S**(-1._r8))*EXP(Z2SD*(THETWA-Z3S))))
    ENDIF
  endif
  end function fDFGS
!------------------------------------------------------------------------------------------
  function fOFFSET(atcs)result(ans)
  implicit none
  real(r8), intent(in) :: atcs   !mean annual soil temperature, [oC]
  real(r8) :: ans

  ans=0.333_r8*(12.5_r8-AMAX1(0.0_r8,AMIN1(25.0_r8,ATCS)))

  end function fOFFSET
end module MiniFuncMod
