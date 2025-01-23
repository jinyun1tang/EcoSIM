module PlantMathFuncMod
!
!DESCRIPTION
  ! code for small functions used by plant processes
  use data_kind_mod, only : r8 => DAT_KIND_R8
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

  pure function get_FDM(PSICanopy)result(FDMP)
  !
  !compute the Ratio of leaf+sheath dry mass to symplasmic water (g g–1)
  !as a function of absolute value of leaf water potential (MPa)
  
  implicit none
  real(r8), intent(in) :: PSICanopy !canopy water potential, MPa
  
  real(r8) :: APSILT    !abosolute value of psi
  real(r8) :: FDMP      !=dry matter/water

  APSILT = ABS(PSICanopy)
  FDMP   = 0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)

  end function get_FDM

!--------------------------------------------------------------------------------
  subroutine update_osmo_turg_pressure(PSICanopy,CCPOLT,OSMO,TKP,PSIOsmo,PSITurg,FDMP1)
  !
  !DESCRIPTION
  !update the osmotic and turgor pressure of a plant organ
  implicit none
  real(r8), intent(in) :: PSICanopy   !plant orgran pressure, MPa
  real(r8), intent(in) :: CCPOLT !total organ dry mass, C+N+P
  real(r8), intent(in) :: OSMO   !canopy osmotic potential when canopy water potential = 0 MPa
  real(r8), intent(in) :: TKP    !organ temperature, Kelvin
  real(r8), intent(out) :: PSIOsmo  !osmotic pressure of the organ, MPa
  real(r8), intent(out) :: PSITurg  !turgor pressure of the organ, MPa
  real(r8), optional, intent(out) :: FDMP1

  real(r8) :: OSWT
  real(r8) :: FDMP

  FDMP=get_FDM(PSICanopy)
  if(present(fdmp1))FDMP1=FDMP
  OSWT    = 36.0_r8+840.0_r8*AZMAX1(CCPOLT)
  PSIOsmo = FDMP/0.16_r8*OSMO-RGASC*TKP*FDMP*CCPOLT/OSWT
  PSITurg = AZMAX1(PSICanopy-PSIOsmo)

  end subroutine update_osmo_turg_pressure
!--------------------------------------------------------------------------------

  subroutine calc_seed_geometry(SeedCMass,SeedVolumeMean,SeedLengthMean,SeedArea)
  !
  !DESCRIPTION
  !assuming the seed is spherical, compute its volume, diameter(=length), and surface area
  !     SeedVolumeMean,SeedLengthMean,SeedAreaMean=seed volume(m3),length(m),area(m2)
  !     SeedCMass=seed C mass (g) from PFT file
  !

  implicit none
  real(r8), intent(in)  :: SeedCMass   !carbon mass per seed
  real(r8), intent(out) :: SeedVolumeMean,SeedLengthMean,SeedArea

  SeedVolumeMean=SeedCMass*5.0E-06_r8
  SeedLengthMean=2.0_r8*(0.75_r8*SeedVolumeMean/PICON)**0.33_r8
  SeedArea=4.0_r8*PICON*(SeedLengthMean/2.0_r8)**2_r8

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
  B     = -(UptakeRateMax_Ol+X-Y+SolDifusFlx*SoluteKM)
  C     = AZMAX1(X-Y)*UptakeRateMax_Ol
  delta = B*B-4.0_r8*C

  if(delta<0._r8)then
    Uptake_Ol=0._r8
  else
    Uptake_Ol=AZMAX1(-B-SQRT(delta))/2.0_r8
  endif

  !Oxygen, and carbon unlimited solute uptake
  BP    = -(UptakeRateMax+X-Y+SolDifusFlx*SoluteKM)
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

  pure function is_plant_treelike(iPlantRootProfile_pft)result(ans)
!
! currently, there are three plant growth types defined as
! iplt_bryophyte=0
! iplt_grasslike=1
!  iplt_treelike=2

  implicit none
  integer, intent(in) :: iPlantRootProfile_pft
  logical :: ans

  ans=iPlantRootProfile_pft > 1
  end function is_plant_treelike

!--------------------------------------------------------------------------------

  pure function is_root_shallow(iPlantRootProfile_pft)result(ans)
!
! currently, there are three plant growth types defined as
! iplt_bryophyte=0
! iplt_grasslike=1
!  iplt_treelike=2

  implicit none
  integer, intent(in) :: iPlantRootProfile_pft
  logical :: ans

  ans=(iPlantRootProfile_pft == iplt_bryophyte)
  end function is_root_shallow
!--------------------------------------------------------------------------------
  pure function pMOD(a,b)result(c)
  !
  ! compute c=MOD(A,B)
  ! if C==0 and A/=0 then C=B
  implicit none
  integer, intent(in) :: a,B
  integer :: C
  
  c=mod(a,b)
  if(c==0 .and. a/=0)c=b
  end function pMOD
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
  real(r8), intent(in) :: WFN   !turgor based leaf/root expansion 
  integer, intent(in) :: iPlantRootProfile
  real(r8) :: ans

  IF(is_root_shallow(iPlantRootProfile))THEN
    ans=WFN**0.10_r8
  ELSE
    ans=WFN**0.25_r8
  ENDIF

  end function fRespWatSens
  
end module PlantMathFuncMod
