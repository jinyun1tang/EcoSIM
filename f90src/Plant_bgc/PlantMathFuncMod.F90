module PlantMathFuncMod
!
!DESCRIPTION
  ! code for small functions used by plant processes
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: endrun
  use EcoSimConst
  use MiniMathMod
  use ElmIDMod
implicit none
  character(len=*), parameter, private :: mod_filename=&
  __FILE__

  type, public  :: PlantSoluteUptakeConfig_type
    real(r8) :: SolAdvFlx  
    real(r8) :: SolDifusFlx    
    real(r8) :: UptakeRateMax   
    real(r8) :: O2Stress        
    real(r8) :: PlantPopulation 
    real(r8) :: CAvailStress    
    real(r8) :: SoluteMassMax   
    real(r8) :: SoluteConc      
    real(r8) :: SoluteKM      
    real(r8) :: SoluteConcMin       
  end type PlantSoluteUptakeConfig_type
contains

  function get_FDM(PSIOrgan,FDMP0)result(FDMP)
  !
  !compute the Ratio of leaf+sheath dry mass to symplasmic water (g g–1)
  !as a function of absolute value of leaf water potential (MPa)
  
  implicit none
  real(r8), intent(in) :: PSIOrgan !canopy water potential, MPa
  real(r8), optional, intent(out) ::  FDMP0
  
  real(r8) :: APSILT    !abosolute value of psi
  real(r8) :: FDMP      !=dry matter/water

  APSILT = ABS(PSIOrgan)
  FDMP   = 0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
  if(present(FDMP0))FDMP0=0.16_r8

  end function get_FDM

!--------------------------------------------------------------------------------
  subroutine update_osmo_turg_pressure(PSIOrgan,CCPOLT,OSMO,TKP,PSIOsmo,PSITurg,FDMP1)
  !
  !DESCRIPTION
  !update the osmotic and turgor pressure of a plant organ
  implicit none
  real(r8), intent(in) :: PSIOrgan    !plant orgran pressure, [MPa]
  real(r8), intent(in) :: CCPOLT      !total organ non-structrual elemental concentration, g/g
  real(r8), intent(in) :: OSMO        !Organ osmotic potential when water potential = 0 MPa
  real(r8), intent(in) :: TKP         !organ temperature, Kelvin
  real(r8), intent(out) :: PSIOsmo    !osmotic pressure of the organ, MPa
  real(r8), intent(out) :: PSITurg    !turgor pressure of the organ, MPa
  real(r8), optional, intent(out) :: FDMP1

  real(r8) :: OSWT   !average molecular weight of CCPOLT 
  real(r8) :: FDMP   !Ratio of leaf+sheath dry mass to symplasmic water (g/g)
  real(r8) :: FDMP0  !FDMP at zero canopy water potential. 
  
  FDMP    = get_FDM(PSIOrgan,FDMP0)
  OSWT    = 36.0_r8+840.0_r8*AZMAX1(CCPOLT)
  PSIOsmo = FDMP*(OSMO/FDMP0-RGASC*TKP*CCPOLT/OSWT)
  PSITurg = AZMAX1(PSIOrgan-PSIOsmo)

  if(present(fdmp1))FDMP1=FDMP
  end subroutine update_osmo_turg_pressure
!--------------------------------------------------------------------------------

  function get_zero_turg_ccpolt(PSIOrg,OSMO,TKP)result(CCPOLT)
  implicit none
  real(r8), intent(in) :: PSIOrg      !plant orgran pressure, [MPa]
  real(r8), intent(in) :: OSMO        !Organ osmotic potential [MPa], when water potential = 0 MPa
  real(r8), intent(in) :: TKP         !organ temperature, [Kelvin]

  real(r8) :: CCPOLT      !total organ non-structrual elemental concentration, g/g
  !--- local variables ---
  real(r8) :: FDMP,FDMP0,ratio

  FDMP  = get_FDM(PSIOrg,FDMP0)

  ratio = RGASC*TKP/(OSMO/FDMP0-PSIOrg/FDMP)

  CCPOLT = 36._r8/(ratio-840._r8)
  end function get_zero_turg_ccpolt
!--------------------------------------------------------------------------------

  subroutine calc_seed_geometry(SeedCMass,rwidth2lenSeed,SeedVolumeMean,SeedLengthMean,SeedArea)
  !
  !DESCRIPTION
  !assuming the seed is spherical, compute its volume, diameter(=length), and surface area
  !     SeedVolumeMean,SeedLengthMean,SeedAreaMean=seed volume(m3),length(m),AREA_3D(m2)
  !     SeedCMass=seed C mass (g) from PFT file
  !

  implicit none
  real(r8), intent(in)  :: SeedCMass        !carbon mass per seed, gC/seed
  real(r8), intent(in)  :: rwidth2lenSeed   !seed width to length ratio
  real(r8), intent(out) :: SeedVolumeMean,SeedLengthMean,SeedArea
  real(r8) :: seedHalfLength
  real(r8), parameter :: pp=1.6075_r8

  SeedVolumeMean = SeedCMass*5.0E-06_r8
  seedHalfLength=(0.75_r8*SeedVolumeMean/PICON)**0.33_r8*rwidth2lenSeed**(-0.667_r8) !assume prolate
  SeedLengthMean = 2.0_r8*seedHalfLength  
  SeedArea       = 4.0_r8*PICON*seedHalfLength**2*((2._r8*rwidth2lenSeed**pp+rwidth2lenSeed**(2._r8*pp))/3._r8)**(1._r8/pp)

  end subroutine calc_seed_geometry
!--------------------------------------------------------------------------------
  pure function calc_root_grow_tempf(TKSO)result(fT_root)
  !
  !DESCRIPTION
  !compute the temperature dependence for plant root growth
  implicit none
  real(r8), intent(in) :: TKSO   !apparent temperature felt by the root
  real(r8) :: fT_root
  real(r8) :: RTK,STK,ACTV

  RTK=RGASC*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  FT_ROOT=EXP(25.229_r8-62500._r8/RTK)/ACTV

  end function calc_root_grow_tempf
!--------------------------------------------------------------------------------
  pure function calc_canopy_grow_tempf(TKGO)result(fT_canp)
  !
  !DESCRIPTION
  !compute the temperature dependence for plant canopy growth
  implicit none
  real(r8), intent(in) :: TKGO   !apparent temperature felt by the canopy
  real(r8) :: fT_canp
  real(r8) :: RTK,STK,ACTV

  RTK     = RGASC*TKGO
  STK     = 710.0_r8*TKGO
  ACTV    = 1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  FT_canp = EXP(25.229_r8-62500._r8/RTK)/ACTV

  end function calc_canopy_grow_tempf

!--------------------------------------------------------------------------------
  pure function calc_leave_grow_tempf(TKCO)result(TFNP)
  !DESCRIPTION
  !compute the temperature dependence for leave growth
  implicit none
  real(r8), intent(in) :: TKCO   !apparent temperature felt by the leave
  real(r8) :: TFNP
  real(r8) :: RTK,STK,ACTV
  
  RTK  = RGASC*TKCO
  STK  = 710.0_r8*TKCO
  ACTV = 1+EXP((197500_r8-STK)/RTK)+EXP((STK-218500._r8)/RTK)
  TFNP = EXP(24.269_r8-60000._r8/RTK)/ACTV
  end function calc_leave_grow_tempf

!--------------------------------------------------------------------------------
  pure function calc_plant_maint_tempf(TKCM)result(TFN5)
  implicit none
  real(r8), intent(in) :: TKCM

  real(r8) :: TFN5
  real(r8) :: RTK,STK,ACTVM

  RTK   = RGASC*TKCM
  STK   = 710.0_r8*TKCM
  ACTVM = 1._r8+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TFN5  = EXP(25.214_r8-62500._r8/RTK)/ACTVM
  END function calc_plant_maint_tempf
!--------------------------------------------------------------------------------

  subroutine SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig, PltUptake_Ol, &
    PltUptake_Sl, PltUptake_OSl, PltUptake_OSCl,ldebug)
  !
  !DESCRIPTION
  !solve for substrate uptake rate as a function of solute concentration
  !
  !Q^2−(v+X-Y+DK)Q+(X−Y)v=0
  !Q is uptake rate
  !v is maximum uptake rate
  !K is affinity parameter
  !X=(q+D)C, with C as micropore solute concentration
  !Y=D*Cm, with Cm being the minimum concentration for uptake

  implicit none
  type(PlantSoluteUptakeConfig_type), intent(in) :: PlantSoluteUptakeConfig
  real(r8), intent(out) :: PltUptake_Ol     !oxygen limited but solute or carbon unlimited
  real(r8), intent(out) :: PltUptake_Sl     !oxygen and carbon unlimited but solute limited uptake
  real(r8), intent(out) :: PltUptake_OSl    !oxygen and solute limited, but not carbon limited
  real(r8), intent(out) :: PltUptake_OSCl   !oxygen, solute and carbon limited uptake
  logical, optional, intent(in) :: ldebug
  real(r8) :: UptakeRateMax_Ol   !oxygen limited maximum uptake rate
  real(r8) :: X, Y, B, C, BP, CP, delta
  real(r8) :: Uptake, Uptake_Ol
  logical :: lldebug
  associate(                                                    &
  SolAdvFlx       => PlantSoluteUptakeConfig%SolAdvFlx        , &
  SolDifusFlx     => PlantSoluteUptakeConfig%SolDifusFlx      , &
  UptakeRateMax   => PlantSoluteUptakeConfig%UptakeRateMax    , &
  O2Stress        => PlantSoluteUptakeConfig%O2Stress         , &
  PlantPopulation => PlantSoluteUptakeConfig%PlantPopulation  , &
  CAvailStress    => PlantSoluteUptakeConfig%CAvailStress     , &
  SoluteMassMax   => PlantSoluteUptakeConfig%SoluteMassMax    , &
  SoluteConc      => PlantSoluteUptakeConfig%SoluteConc       , &
  SoluteKM        => PlantSoluteUptakeConfig%SoluteKM         , &
  SoluteConcMin   => PlantSoluteUptakeConfig%SoluteConcMin      &
  )
  lldebug=.false.
  if(present(ldebug))lldebug=ldebug

  UptakeRateMax_Ol=UptakeRateMax*O2Stress  


  X=(SolDifusFlx+SolAdvFlx)*SoluteConc
  Y=SolDifusFlx*SoluteConcMin

  !Oxygen limited but not solute or carbon limited uptake
  ! u^2+Bu+C=0., it requires when C=0, delta=1, u=0
  B     = -AZMAX1(UptakeRateMax_Ol+X-Y+SolDifusFlx*SoluteKM)
  C     = AZMAX1(X-Y)*UptakeRateMax_Ol
  delta = B*B-4.0_r8*C

  if(delta<0._r8)then
    Uptake_Ol=0._r8
  else
    Uptake_Ol=AZMAX1(-B-SQRT(delta))/2.0_r8
  endif

  !Oxygen, and carbon unlimited solute uptake
  BP    = -AZMAX1(UptakeRateMax+X-Y+SolDifusFlx*SoluteKM)
  CP    = AZMAX1(X-Y)*UptakeRateMax
  delta = BP*BP-4.0_r8*CP
  if(delta<0._r8)then
    Uptake=0._r8
  else
    Uptake=AZMAX1(-BP-SQRT(delta))/2.0_r8
  endif
  if(lldebug)write(115,*)'delta2',delta,Uptake,'BP=',BP,CP

  !oxygen and solute limited but carbon unlimited
  PltUptake_Ol=AZMAX1(Uptake_Ol*PlantPopulation)

  !oxygen and solute limited, but not carbon limited
  PltUptake_OSl=AMIN1(SoluteMassMax,PltUptake_Ol)

  !oxygen and carbon unlimited but solute limited uptake
  PltUptake_Sl=AMIN1(SoluteMassMax,Uptake*PlantPopulation)

  !oxygen, solute and carbon limited uptake
  PltUptake_OSCl=PltUptake_OSl/CAvailStress

  end associate
  end subroutine SoluteUptakeByPlantRoots
!--------------------------------------------------------------------------------

  pure function is_plant_woody_vascular(iPlantRootProfile_pft)result(ans)
!
! currently, there are three plant growth types defined as
! iplt_bryophyte=0
! iplt_grasslike=1
!  iplt_treelike=2

  implicit none
  integer, intent(in) :: iPlantRootProfile_pft
  logical :: ans

  ans=iPlantRootProfile_pft > 1
  end function is_plant_woody_vascular

!--------------------------------------------------------------------------------

  pure function is_root_shallow(iPlantRootProfile_pft)result(ans)
!
! currently, there are three plant growth types defined as
! iplt_bryophyte=0
! iplt_grasslike=1
!  iplt_treelike=2
! only bryophyte is considered as shallow roots
  implicit none
  integer, intent(in) :: iPlantRootProfile_pft
  logical :: ans

  ans=(iPlantRootProfile_pft == iplt_bryophyte)
  end function is_root_shallow

!--------------------------------------------------------------------------------
  pure function is_root_N2fix(iPlantNfixType_pft)result(yesno)
  implicit none
  integer, intent(in) :: iPlantNfixType_pft

  logical :: yesno

  yesno=iPlantNfixType_pft.GE.in2fixtyp_root_fast.AND.iPlantNfixType_pft.LE.in2fixtyp_root_slow

  end function is_root_N2fix
!--------------------------------------------------------------------------------
  pure function is_canopy_N2fix(iPlantNfixType_pft)result(yesno)
  implicit none
  integer, intent(in) :: iPlantNfixType_pft

  logical :: yesno

  yesno=iPlantNfixType_pft.GE.in2fixtyp_canopy_fast.AND.iPlantNfixType_pft.LE.in2fixtyp_canopy_slow

  end function is_canopy_N2fix
!--------------------------------------------------------------------------------
  pure function is_plant_N2fix(iPlantNfixType_pft)result(yesno)
  implicit none
  integer, intent(in) :: iPlantNfixType_pft

  logical :: yesno

  yesno=iPlantNfixType_pft.NE.iN2fixtyp_none

  end function is_plant_N2fix
!--------------------------------------------------------------------------------

  pure function fRespWatSens(WFN,iPlantRootProfile)result(ans)

  implicit none
  real(r8), intent(in) :: WFN               !turgor based leaf/root elongation
  integer, intent(in) :: iPlantRootProfile
  real(r8) :: ans

  IF(is_root_shallow(iPlantRootProfile))THEN
    ans=WFN**0.10_r8
  ELSE
    ans=WFN**0.25_r8
  ENDIF

  end function fRespWatSens

!----------------------------------------------------------------------------------------------------

  subroutine ExchFluxLimiter(fromState,toState,XFRE)
  implicit none
  real(r8),intent(in) :: fromState
  real(r8),intent(in) :: toState
  real(r8), intent(inout) :: XFRE

  IF(XFRE>0._r8)then
    XFRE=AMIN1(fromState*0.9999_r8,XFRE)
  ELSE  
    XFRE=AMAX1(-toState*0.9999_r8,XFRE)
  ENDIF
  end subroutine ExchFluxLimiter
!----------------------------------------------------------------------------------------------------

  subroutine advect_remap_mass_loss(n, xr, c, ur, dt, c_new, xL, lost_mass)
    !--------------------------------------------------------------------
    ! Conservative 1D advect-remap with proportional mass loss at right boundary
    ! Workspace arrays must be preallocated by the caller (no allocate/deallocate here).
    !
    ! Inputs:
    !   n   - number of cells
    !   xr  - right-edge locations of each cell (size n), strictly increasing
    !   c   - cell-average concentration in each cell (size n)
    !   ur  - velocities at right-edge of each cell (size n), nonnegative
    !   dt  - time step (positive)
    !   xL  - left boundary (default 0.0 in caller)
    !
    ! Outputs:
    !   c_new    - updated cell-average concentration (size n)
    !   lost_mass- total mass lost this step
    !
    ! Workspace (caller provides arrays):
    !   xE, xE_star : real(r8), size n+1
    !   dx, m, m_kept, frac_kept : real(r8), size n
    !   M_star, M_on_fixed : real(r8), size n+1
    !--------------------------------------------------------------------
    implicit none

    integer, intent(in) :: n
    real(r8), intent(in) :: xr(n), c(n), ur(n), dt 
    real(r8), optional, intent(in):: xL
    real(r8), intent(out) :: c_new(n)
    real(r8), optional, intent(out) :: lost_mass

    ! workspace arrays provided by caller (no allocate/deallocate here)
    real(r8)  :: xE(n+1), xE_star(n+1)
    real(r8)  :: dx(n), m(n), m_kept(n), frac_kept(n)
    real(r8)  :: M_star(n+1), M_on_fixed(n+1)
    
    ! local scalars
    integer :: i
    real(r8) :: xR_most
    real(r8) :: total_initial, total_kept
    real(r8) :: a, b, w_star, inside_left, inside_right, overlap, frac
    real(r8) :: dt_res,dt_loc
    logical :: lhalf
    real(r8), parameter :: tiny = 1.0d-14

    ! Basic checks (lightweight)
    if (n <= 0) then
       if(present(lost_mass))lost_mass = 0._r8
       return
    end if
    if (dt <= 0._r8) then
       call endrun('Error: dt must be positive.  in '//trim(mod_filename)//' at line',__LINE__)                    
    end if

    do i = 2, n
       if (xr(i) <= xr(i-1)) then
          call endrun('Error: xr must be strictly increasing.  in '//trim(mod_filename)//' at line',__LINE__)                              
       end if
    end do
    do i = 1, n
       if (ur(i) < 0._r8) then
          call endrun('Error: ur must be nonnegative.  in '//trim(mod_filename)//' at line',__LINE__)                    
       end if
    end do

    ! rightmost fixed boundary
    xR_most = xr(n)

    ! fixed edges: xE(1) = xL, xE(2:n+1) = xr(1:n)
    if(present(xL))then
      xE(1) = xL
    else
      xE(1) = 0._r8
    endif
    do i = 1, n
       xE(i+1) = xr(i)
    end do

    ! compute dx and masses
    total_initial = 0._r8
    do i = 1, n
       dx(i) = xE(i+1) - xE(i)
       if (dx(i) <= 0._r8) then
          call endrun('Error: non-positive cell width encountered.  in '//trim(mod_filename)//' at line',__LINE__)                    
       end if
    end do
    m = c * dx
    total_initial = total_initial + sum(m)

  dt_res=dt; dt_loc=dt  
  DO
    ! compute moved edges (left edge fixed velocity = 0)

    do
      lhalf=.false.    
      xE_star(1) = xE(1) + 0._r8 * dt_loc      
      do i = 1, n
        xE_star(i+1) = xE(i+1) + ur(i) * dt_loc
        if (xE_star(i+1) + tiny < xE_star(i) .and. xE_star(i+1)<xR_most) then
          dt_loc=dt_loc*0.5_r8
          lhalf=.true.
          exit
        endif
      end do
      if(.not.lhalf)exit 
    enddo

    ! monotonicity check (advected edges should be increasing)
    do i = 2, n+1
      if (xE_star(i) + tiny < xE_star(i-1) .and. xE_star(i)<xR_most) then
        write(*,*)'ur',ur
        write(*,*)'xR_most',xR_most,dt_loc
        write(*,*)'xE_star',xE_star(2:n+1) 
        call endrun('Error: advected edges are not monotone increasing. Reduce dt.  in '//trim(mod_filename)//' at line',__LINE__)                    
      end if
    end do

    ! compute kept mass per moved cell via overlap with [xL, xR_most]
    total_kept = 0._r8
    do i = 1, n
       a = xE_star(i)
       b = xE_star(i+1)
       w_star = b - a
       if (w_star <= 0._r8) then
          frac_kept(i) = 0._r8
          m_kept(i) = 0._r8
       else
          inside_left = max(a, xL)
          inside_right = min(b, xR_most)
          overlap = inside_right - inside_left
          if (overlap < 0._r8) overlap = 0._r8
          frac = overlap / w_star
          if (frac < 0._r8) frac = 0._r8
          if (frac > 1.0d0) frac = 1.0d0
          frac_kept(i) = frac
          m_kept(i) = m(i) * frac
       end if
       total_kept = total_kept + m_kept(i)
    end do
    if(present(lost_mass))lost_mass = total_initial - total_kept

    ! build cumulative M_star on moved mesh edges
    M_star(1) = 0._r8
    do i = 1, n
       M_star(i+1) = M_star(i) + m_kept(i)
    end do

    ! interpolate M_star back to fixed edges xE -> M_on_fixed
    call interp_linear_clamped(n+1, xE_star, M_star, n+1, xE, M_on_fixed, 0._r8, total_kept)

    ! new cell masses and concentrations
    do i = 1, n
       m(i) = M_on_fixed(i+1) - M_on_fixed(i)    ! reuse m for new masses
       c_new(i) = m(i) / dx(i)
    end do
    dt_res=dt_res-dt_loc
    if(dt_res<dt*1.e-2_r8)exit
    dt_loc=dt_res
  enddo  

  end subroutine advect_remap_mass_loss

  !--------------------------------------------------------------------
  ! Simple piecewise-linear interpolation with clamped extrapolation:
  !   Inputs:
  !     nx  - length of x array
  !     x   - array of x nodes, length nx (must be nondecreasing)
  !     y   - array of y nodes, length nx
  !     nq  - number of query points (size of xq)
  !     xq  - query points (length nq)
  !   Outputs:
  !     yq  - interpolated values at xq (length nq)
  !   Extrapolation:
  !     xq <= x(1) -> y_left
  !     xq >= x(nx) -> y_right
  !--------------------------------------------------------------------
  subroutine interp_linear_clamped(nx, x, y, nq, xq, yq, y_left, y_right)
    implicit none
    integer, intent(in) :: nx, nq
    real(r8), intent(in) :: x(nx), y(nx), xq(nq)
    real(r8), intent(out) :: yq(nq)
    real(r8), intent(in) :: y_left, y_right

    integer :: iq, k
    real(r8) :: t

    do iq = 1, nq
       if (xq(iq) <= x(1)) then
          yq(iq) = y_left
       else if (xq(iq) >= x(nx)) then
          yq(iq) = y_right
       else
          ! find k s.t. x(k) <= xq < x(k+1)
          k = 1
          do while (k < nx - 1 .and. xq(iq) >= x(k+1))
             k = k + 1
          end do
          if (x(k+1) == x(k)) then
             yq(iq) = y(k)
          else
             t = (xq(iq) - x(k)) / (x(k+1) - x(k))
             yq(iq) = (1.0d0 - t) * y(k) + t * y(k+1)
          end if
       end if
    end do
  end subroutine interp_linear_clamped

end module PlantMathFuncMod
