Module SharedDataMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  implicit none
  character(len=*), private, parameter :: mod_filename= &
  __FILE__
!figure out the grid for ATS
!  integer, intent(in) :: jzsoi            !number of soil layers
!  integer, intent(in) :: js               !number of snow layers
!  integer, intent(in) :: jzsoi            !number of soil layers
!  integer, intent(in) :: js               !number of snow layers

! temporary data holder in ecosim
  real(r8) :: atm_n2, atm_o2,atm_co2,atm_ch4,atm_N2o,atm_H2,atm_NH3
  real(r8) :: heat_capacity, pressure_at_field_capacity, pressure_at_wilting_point
  real(r8), allocatable :: a_csand(:,:)   !sand mass fraction
  real(r8), allocatable :: a_CSILT(:,:)   !silt mass fraction
  real(r8), allocatable :: a_BKDSI(:,:)   !bulk density
  real(r8), allocatable :: a_CumDepz2LayBottom_vr(:,:)   !dpeth (from surfce to bottom)
  real(r8), allocatable :: a_Volume(:,:)   !volume
  real(r8), allocatable :: a_dz(:,:)      !distance between layers  
  real(r8), allocatable :: a_AreaZ(:,:)   !Area normal to z axis
  real(r8), allocatable :: a_FC(:,:)      !field capacity
  real(r8), allocatable :: a_WP(:,:)      !wilting point
  real(r8), allocatable :: a_FHOL(:,:)    !macropore fraction
  real(r8), allocatable :: a_ROCK(:,:)    !mass fraction as rock
  real(r8), allocatable :: a_CORGC(:,:)   !organic carbon content
  real(r8), allocatable :: a_CORGN(:,:)   !organic nitrogen content
  real(r8), allocatable :: a_CORGP(:,:)   !organic phosphorus content
  real(r8), allocatable :: a_PORO(:,:)    !Porosity
  real(r8), allocatable :: a_MATP(:,:)    !Matric Pressure
!  real(r8), allocatable ::a_CORGR(:,:)   !organic nitrogen  content
  real(r8), allocatable :: a_ASP(:)       !Aspect
  real(r8), allocatable :: a_ALT(:)       !Altitude
  real(r8), allocatable :: a_ATKA(:)
  real(r8), allocatable :: a_WC(:,:)      !Soil water content
  real(r8), allocatable :: a_LSAT(:,:)    !liquid saturation
  real(r8), allocatable :: a_RELPERM(:,:) !relative_permeability
  real(r8), allocatable :: a_HCOND(:,:)   !hydraulic conductivity
  real(r8), allocatable :: a_TEMP(:,:)    !temperature
  real(r8), allocatable :: a_TEMP_pad(:,:)    !temperature
  real(r8), allocatable :: a_SSES(:,:)    !subsurface energy source
  real(r8), allocatable :: a_SSWS(:,:)    !subsurface water source
  real(r8), allocatable :: a_AREA3(:,:) 
  real(r8), allocatable :: a_AREA3_pad(:,:)

  real(r8), allocatable :: tairc(:)       !air temperature oC
  real(r8), allocatable :: uwind(:)       !wind speed, m/s
  real(r8), allocatable :: p_rain(:)      !precipitation, mm H2O/hr
  real(r8), allocatable :: p_snow(:)      !precipitation snow, mm H2O/hr 
  real(r8), allocatable :: sunrad(:)      !solar radiation,
  real(r8), allocatable :: swrad(:)       !shortwave radiation,
  real(r8), allocatable :: vpair(:)       !vapor pressure deficit
  real(r8), allocatable :: surf_e_source(:) !surface energy source
  real(r8), allocatable :: surf_w_source(:) !surface water source
  real(r8), allocatable :: surf_snow_depth(:) !snow depth source
  !real(r8), allocatable :: a_AREA3(:)
  integer,  allocatable :: a_NU(:)        !upper soil layer index
  integer,  allocatable :: a_NL(:)        !lower soil layer index
  integer,  allocatable :: a_NJ(:)
  integer,  allocatable :: a_MaxNumRootLays(:)
  integer :: NYS, I                       !total number of columns
  contains

  subroutine InitSharedData(ncells_per_col_,ncol)
    use GridConsts, only : JX,JY,JZ
    implicit none
    integer, intent(in) :: ncells_per_col_   !number of vertical layers
    integer, intent(in) :: ncol  !NUMBER of cols
    !set # of soil layers
    !JZSOI=JZs
    write(*,*) "(InitSharedData) JX, JY, JZ, N_cells, N_cols: ", JX, JY, JZ, &
            ncells_per_col_, ncol
    
    JX=1;JY=ncol
    JZ = ncells_per_col_
   
    !write(*,*) "In Shared data after setting: "
    !write(*,*) "JX=", JX, ", JY=", JY, ", JZ=", JZ

    allocate(a_csand(ncells_per_col_,ncol))
    allocate(a_CSILT(ncells_per_col_,ncol))   !silt mass fraction
    !allocate(a_AreaZ(ncells_per_col_,ncol))   !actually need to allocate area
    !allocate(a_BKDSI(ncells_per_col_,ncol))   !bulk density
    !allocate(a_CumDepz2LayBottom_vr(ncells_per_col_,ncol))   !dpeth (from surfce to bottom)
    !allocate(a_FC(ncells_per_col_,ncol))      !field capacity
    !allocate(a_WP(ncells_per_col_,ncol))      !wilting point
    allocate(a_FHOL(ncells_per_col_,ncol))    !macropore fraction
    allocate(a_ROCK(ncells_per_col_,ncol))    !mass fraction as rock
    allocate(a_CORGC(ncells_per_col_,ncol))   !organic carbon content
    allocate(a_CORGN(ncells_per_col_,ncol))   !organic nitrogen content
    allocate(a_CORGP(ncells_per_col_,ncol))   !organic phosphorus content
    !allocate(a_PORO(ncells_per_col_,ncol))
    !allocate(a_AREA3(ncells_per_col_))
    allocate(a_NU(ncol))
    allocate(a_NL(ncol))
    !allocate(a_ASP(ncells_per_col_))
    allocate(a_ALT(ncells_per_col_))
    !allocate(tairc(1:ncells_per_col_))
    !allocate(uwind(1:ncells_per_col_))
    !allocate(prec(1:ncells_per_col_))
    !allocate(sunrad(1:ncells_per_col_))
    !allocate(vpair(1:ncells_per_col_))
    allocate(a_ATKA(1:ncells_per_col_))
 
    !allocate(a_BKDSI(ncells_per_col_, ncol))        ! bulk density
    !allocate(a_CumDepth2LayerBottom(ncells_per_col_, ncol))  ! depth (from surface to bottom)
    !allocate(a_Volume(ncells_per_col_, ncol))       ! volume
    !allocate(a_dz(ncells_per_col_, ncol))           ! distance between layers
    !allocate(a_AreaZ(ncells_per_col_, ncol))        ! area normal to z-axis
    !allocate(a_FC(ncells_per_col_, ncol))           ! field capacity
    !allocate(a_WP(ncells_per_col_, ncol))           ! wilting point
    !allocate(a_PORO(ncells_per_col_, ncol))         ! porosity
    !allocate(a_MATP(ncells_per_col_, ncol))         ! matric pressure
    !allocate(a_WC(ncells_per_col_, ncol))           ! soil water content
    !allocate(a_LSAT(ncells_per_col_, ncol))         ! liquid saturation
    !allocate(a_RELPERM(ncells_per_col_, ncol))      ! relative permeability
    !allocate(a_HCOND(ncells_per_col_, ncol))        ! hydraulic conductivity
    !allocate(a_TEMP(ncells_per_col_, ncol))         ! temperature
    !allocate(a_SSES(ncells_per_col_, ncol))         ! subsurface energy source
    !allocate(a_SSWS(ncells_per_col_, ncol))         ! subsurface water source
    !allocate(a_AREA3(ncells_per_col_, ncol)) 
    !allocate(a_AREA3_pad(ncells_per_col_, ncol))
    !allocate(a_TEMP_pad(ncells_per_col_, ncol))

    !allocate(a_ASP(ncells_per_col_))                ! aspect
    !allocate(tairc(ncells_per_col_))                ! air temperature
    !allocate(uwind(ncells_per_col_))                ! wind speed
    !allocate(p_rain(ncells_per_col_))               ! precipitation (rain)
    !allocate(p_snow(ncells_per_col_))               ! precipitation (snow)
    !allocate(sunrad(ncells_per_col_))               ! solar radiation
    !allocate(swrad(ncells_per_col_))                ! shortwave radiation
    !allocate(vpair(ncells_per_col_))                ! vapor pressure deficit
    !allocate(surf_e_source(ncells_per_col_))        ! surface energy source
    !allocate(surf_w_source(ncells_per_col_))        ! surface water source
    !allocate(surf_snow_depth(ncells_per_col_))        ! surface snow depth

    a_NU=1
    a_NL=ncells_per_col_

    write(*,*) "(InitSharedData f) JX, JY, JZ, N_cells, N_cols: ", JX, JY, JZ, &
        ncells_per_col_, ncol
  end subroutine InitSharedData

!------------------------------------------------------------------------------------------

  subroutine DestroySharedData()
  implicit none

  call destroy(a_csand)
  call destroy(a_CSILT)
  call destroy(a_BKDSI)
  call destroy(a_MATP)
  call destroy(a_AreaZ)
  call destroy(a_CumDepz2LayBottom_vr)
  call destroy(a_FC)
  call destroy(a_WP)
  call destroy(a_FHOL)
  call destroy(a_ROCK)
  call destroy(a_CORGC)
  call destroy(a_CORGN)
  call destroy(a_CORGP)
  call destroy(a_PORO)
  call destroy(a_AREA3)
  call destroy(a_ASP)
  call destroy(a_ALT)
  call destroy(tairc)
  call destroy(uwind)
  call destroy(p_rain)
  call destroy(p_snow)
  call destroy(sunrad)
  call destroy(swrad)
  call destroy(vpair)
  call destroy(a_ATKA)
  call destroy(a_AREA3_pad)
  call destroy(a_TEMP_pad)
  end subroutine DestroySharedData

end module SharedDataMod
