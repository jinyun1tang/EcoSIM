program EcoATSTest
  use ATSCPLMod
  use BGCContainers_module
  use SharedDataMod
  use GridDataType
  implicit none
  
  ! Declare variables
  type (BGCState) :: state
  type (BGCProperties) :: props
  type (BGCSizes) :: sizes
  integer :: NY, NX, L
  integer :: ncells_per_col_, ncol
   
  NX = 1
  NYS = 1
  ncells_per_col_ = 100
  ncol = 1
  
  ! Initialize sizes structure
  sizes%num_components = 1
  sizes%ncells_per_col_ = 100
  sizes%num_columns = 1


  call Init_ATSEcoSIM_driver()

  call Init_EcoSIM(sizes)

  ! Call Run_EcoSIM_one_step to advance EcoSIM
  call Run_EcoSIM_one_step(sizes)

end program EcoATSTest

! ----------------------------------------------------------------------------

subroutine Init_ATSEcoSIM_driver()
  use BGCContainers_module
  use SharedDataMod
  use GridDataType
  implicit none
  
  ! Declare variables
  type (BGCState) :: state
  type (BGCProperties) :: props
  type (BGCSizes) :: sizes
  integer :: NY, NX, L
  integer :: ncells_per_col_, ncol
  NX = 1
  NYS = 1
  ncells_per_col_ = 100
  ncol = 1

  ! Initialize sizes structure
  sizes%num_components = 1
  sizes%ncells_per_col_ = 100
  sizes%num_columns = 1

  ! needed in starts
  pressure_at_field_capacity = -0.001
  pressure_at_wilting_point = -0.001

  !need to allocate these because c_f_pointer sets them in the actual coupler
  allocate(a_ASP(ncells_per_col_))
  allocate(tairc(1:ncells_per_col_))
  allocate(vpair(1:ncells_per_col_))
  allocate(uwind(1:ncells_per_col_))
  allocate(swrad(1:ncells_per_col_))
  allocate(sunrad(1:ncells_per_col_))
  allocate(p_rain(1:ncells_per_col_))
  allocate(surf_e_source(1:ncells_per_col_))
  allocate(surf_w_source(1:ncells_per_col_))
  allocate(surf_snow_depth(1:ncells_per_col_))
  allocate(a_TEMP(ncells_per_col_, ncol))
  allocate(a_CumDepz2LayerBot_vr(ncells_per_col_, ncol))
  allocate(a_AREA3(ncells_per_col_, ncol))
  allocate(a_BKDSI(ncells_per_col_, ncol))
  allocate(a_WC(ncells_per_col_, ncol))
  allocate(a_MATP(ncells_per_col_, ncol))
  allocate(a_PORO(ncells_per_col_, ncol))

  do NY=1,NYS
    do L=1,ncells_per_col_
      a_AREA3(L,NY) = 0.1
      !DH(NY,NX) = 0.316229
      !DV(NY,NX) = 0.316229
      a_CumDepz2LayerBot_vr(L,NY) = 0.0
      a_BKDSI(L,NY) = 26000.0/1.0e3_r8
      !a_CORGC(L,NY) = 0.0
      !a_CORGN(L,NY) = 0.0
      !a_CORGP(L,NY) = 0.0
      a_WC(L,NY) = 2500.0
      a_TEMP(L,NY) = 265.0
      a_MATP(L,NY) = 100.0
      a_PORO(L,NY) = 0.5
    enddo
    a_ASP(NY) = 0.0
    tairc(NY) = 265.0
    vpair(NY) = 736.3/1.0e6_r8
    uwind(NY) = 1.0*3600.0_r8
    swrad(NY) = 400.0*0.0036_r8
    sunrad(NY) =  219.78*0.0036_r8
    
    p_rain(NY) = 0.0
    !p_rain(NY) = 3.e-8*1000.0_r8*3600.0_r8
  enddo

end subroutine Init_ATSEcoSIM_driver

