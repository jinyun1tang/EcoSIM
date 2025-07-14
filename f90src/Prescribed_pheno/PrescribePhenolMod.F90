module PrescribePhenolMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSimConst,   only: PICON2h
  use EcoSIMCtrlMod, only: etimer
  use EcoSIMCtrlDataType, only : ZEROS
  use ElmIDMod
  use GridDataType
  use RootDataType
  use PlantAPIData
  use ClimForcDataType
  use PlantTraitDataType
  use PlantMgmtDataType
  use CanopyDataType
  use CanopyRadDataType
implicit none

  !pft information from CLM4.5
  !pft                                                 ra             rb           ztop(m)      zbot(m)
  ! 1. Needleleaf evergreen tree/NET Temperate         7.0            2.0          17           8.5
  ! 2. Needleleaf evergreen tree/NET Boreal            7.0            2.0          17           8.5
  ! 3. Needleleaf deciduous tree/NDT Boreal            7.0            2.0          14           7
  ! 4. Broadleave evergreen tree/BET tropical          7.0            1.0          35           1
  ! 5. Broadleave evergreen tree/BET temperate         7.0            1.0          35           1
  ! 6. Broadleave deciduous tree/BDT tropical          6.0            2.0          18           10
  ! 7. Broadleave deciduous tree/BDT temperate         6.0            2.0          20           11.5
  ! 8. Broadleave deciduous tree/BDT Boreal            6.0            2.0          20           11.5
  ! 9. Broadleave evergreen shrub/BES temperate        7.0            1.5          0.5          0.1
  !10. Broadleave deciduous shrub/BDS temperate        7.0            1.5          0.5          0.1
  !11. Broadleave deciduous shrub/BDS Boreal           7.0            1.5          0.5          0.1
  !12. C3 Arctic grass                                11.0            2.0          0.5          0.01
  !13. C3 Grass                                       11.0            2.0          0.5          0.01
  !14. C4 grass                                       11.0            2.0          0.5          0.01
  !15. C3 Crop Rainfed                                6.0             3.0          0.5          0.01
  !C3 Crop irrigated                              6.0             3.0          0.5          0.01
  !Corn Rainfed                                   6.0             3.0          0.5          0.01
  !Corn irrigated                                 6.0             3.0          0.5          0.01
  !Rainfed temperate cereals                      6.0             3.0          0.5          0.01
  !irrigated temperate cereals                    6.0             3.0          0.5          0.01
  !Rainfed winter cereals                         6.0             3.0          0.5          0.01
  !irrigated winter cereals                       6.0             3.0          0.5          0.01
  !Rainfed soybean                                6.0             3.0          0.5          0.01
  !irrigated soybean                              6.0             3.0          0.5          0.01
  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: GetRootProfile
  public :: SetCanopyProfile
  public :: PrescribePhenologyInterp
  contains

  subroutine GetRootProfile(ra,rb,NL,ZH,RootProfl)
  implicit none
  real(r8), intent(in) :: ra
  real(r8), intent(in) :: rb
  integer,  intent(in) :: NL
  real(r8), intent(in) :: Zh(0:NL)
  real(r8), intent(out) :: RootProfl(1:NL)
  real(r8) :: rz(0:NL-1)
  integer :: L

  DO L=0,NL
    rz(L)=exp(-ra*zh(L))+exp(-rb*zh(L))
  enddo
  DO L=1,NL-1
    rootProfl(L)=0.5_r8*(rz(L-1)-rz(L))
  ENDDO
  rootProfl(NL)=0.5_r8*rz(NL-1)

  end subroutine GetRootProfile

!------------------------------------------------------------------------

  subroutine get_LAI_profile(zc,NZ)
  !Description:
  !Distribute leaf area in through the canopy
  !  
  implicit none
  real(r8), intent(in) :: ZC  !plant height
  integer, intent(in) :: NZ   !plant index

  end subroutine get_LAI_profile

!------------------------------------------------------------------------

  subroutine get_StemArea_profile(NZ)
  !
  !Description:
  !Distribute steam area in through the canopy
  !
  implicit none
  integer, intent(in) :: NZ   !plant index


  end subroutine get_StemArea_profile

!------------------------------------------------------------------------

  function TreeStem_diameter_taperEq(h,p_stump_rate,p_stemp_rate,top_ht_rate,h_dbh,dbh,h_cb,h_top)result(d)
 
  ! Estimate diameter at height h using segment-specific taper rates as in Larsen (2017)
  implicit none
 
  real(r8), intent(in) :: h             ! - height at which to estimate diameter
  real(r8), intent(in) :: p_stump_rate  ! - stump taper rate [unitless]
  real(r8), intent(in) :: p_stemp_rate  ! - stem taper rate  [uitless]
  real(r8), intent(in) :: top_ht_rate   ! - crown base height/crown length [unitless]
  real(r8), intent(in) :: h_dbh    ! - breast height [m]
  real(r8), intent(in) :: dbh      ! - diameter at breast height [m]
  real(r8), intent(in) :: h_cb     ! - height to crown base [m]
  real(r8), intent(in) :: h_top    ! - total tree height (real)
  !Output
  real(r8) ::   d                  ! - estimated diameter at height h [m]
  
  if(h<=h_dbh)then
    !below breast height
    d=dbh+p_stump_rate*(h-h_dbh)
  elseif(h<=h_cb)then
    !below crown base 
    d=dbh+p_stemp_rate*(h-h_dbh)  
  else
    !above crown base
    d=(h_top-h)*top_ht_rate
  endif

  end function TreeStem_diameter_taperEq
!------------------------------------------------------------------------

  subroutine SetCanopyProfile(I,J,LeafAreaZsec_lpft,StemAreaZsec_lpft)
  !
  !Description:
  !distribute leaf area and stem area within the crown defined by top and bottom heights.
  !For simplicity, unifrom distribution is adopted
  !
  implicit none
  integer , intent(in) :: I,J
  real(r8),intent(out) :: LeafAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)     !leaf area in different zenith sectors in different canopy layers
  real(r8),intent(out) :: StemAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)     !stem area in different zenith sectors in different canopy layers
  real(r8) :: dangle
  integer  :: N,L,NZ,NB

  associate(                                                 &
    CanopyStemAreaZ_pft   => plt_morph%CanopyStemAreaZ_pft  ,& !input  :plant canopy layer stem area, [m2 d-2]
    CanopyLeafAreaZ_pft   => plt_morph%CanopyLeafAreaZ_pft  ,& !inoput :total leaf area, [m2 d-2]
    LeafStalkArea_pft     => plt_morph%LeafStalkArea_pft    ,& !inoput :plant leaf+stem/stalk area, [m2 d-2]
    LeafStalkArea_col     => plt_morph%LeafStalkArea_col    ,& !inoput :stalk area of combined, each PFT canopy,[m^2 d-2]
    StemArea_col          => plt_morph%StemArea_col         ,& !input  :grid canopy stem area, [m2 d-2]
    LeafAngleClass_pft    => plt_morph%LeafAngleClass_pft   ,& !input  :fractionction of leaves in different angle classes, [-]    
    CanopyLeafArea_col    => plt_morph%CanopyLeafArea_col   ,& !input  :grid canopy leaf area, [m2 d-2]
    LeafAreaZsec_brch     => plt_morph%LeafAreaZsec_brch    ,& !input  :leaf surface area, [m2 d-2]    
    NP                    => plt_site%NP                     & !input  :current number of plant species,[-]
  )

  
  StemAreaZsec_lpft  = 0._r8
  LeafAreaZsec_lpft  = 0._r8
  dangle            = PICON2h/real(NumLeafZenithSectors1,r8)         !the angle section width
  LeafStalkArea_col = StemArea_col+CanopyLeafArea_col
  DO NZ=1,NP
    LeafStalkArea_pft(NZ)=0._r8
    NB=1    
    DO L=1,NumCanopyLayers1
      DO N=1,NumLeafZenithSectors1
        LeafAreaZsec_lpft(N,L,NZ)=LeafAngleClass_pft(N,NZ)*CanopyLeafAreaZ_pft(L,NZ)        
      ENDDO      
      StemAreaZsec_lpft(1:NumLeafZenithSectors1,L,NZ)=CanopyStemAreaZ_pft(L,NZ)/real(NumLeafZenithSectors1,kind=r8)      
      LeafStalkArea_pft(NZ)=LeafStalkArea_pft(NZ)+CanopyLeafAreaZ_pft(L,NZ)+CanopyStemAreaZ_pft(L,NZ)         

      !Assuming uniform azimuth desitribution for leaves  
      DO N=1,NumLeafZenithSectors1
        LeafAreaZsec_brch(N,L,L,NB,NZ)=LeafAreaZsec_lpft(N,L,NZ)/real(NumOfLeafAzimuthSectors1,r8)
      ENDDO  
    ENDDO

  ENDDO
  
  end associate
  end subroutine SetCanopyProfile

!------------------------------------------------------------------------------------------     
  subroutine PrescribePhenologyInterp(I, NHW, NHE, NVN, NVS)

  implicit none
  integer, intent(in) :: I
  Integer, intent(in) ::  NHW, NHE, NVN, NVS
  real(r8) :: t
  integer :: it(2)
  integer :: months(2)
  integer :: kmo,dofmon,ndaysmon,NY,NX,NZ,L
  real(r8) :: timwt(2)
  real(r8) :: ZL1(0:NumCanopyLayers)
  real(r8) :: AreaInterval,AreaL
  real(r8) :: ARX  !interval canopy area: leaf+stem
  real(r8) :: DZL  !canopy interval height 
  !============
  !example inputs
  real(r8) :: lai(12)=(/1.1852, 1.1821, 1.1554, 1.2433, 1.2922, 1.3341, 1.2296, 1.4118, 1.4343, 1.3941, 1.2721, 1.2218/)
  real(r8) :: sai(12)=(/0.3190, 0.3058, 0.3058, 0.3032, 0.3058, 0.3117, 0.3433, 0.3032, 0.3084, 0.3292, 0.3656, 0.3249/)
  integer  :: irootType=1
  REAL(R8) :: PerPlantRootC_vr(1:JZ)
  REAL(R8) :: PerPlantRootLen_vr(1:JZ)
  !============

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        !FOR test only
        tlai_mon_pft(:,NZ,NY,NX)       = LAI
        tsai_mon_pft(:,NZ,NY,NX)       = SAI
        height_top_mon_pft(:,NZ,NY,NX) = 17._r8
        LeafAngleClass_pft(:,NZ,NY,NX) = 1./real(NumLeafZenithSectors,kind=r8)
      ENDDO
    ENDDO
  ENDDO    

  ndaysmon = etimer%get_curr_mon_days()
  dofmon   = etimer%get_curr_dom()
  kmo      = etimer%get_curr_mon()
  t = (dofmon-0.5_r8) / ndaysmon
  it(1) = t + 0.5_r8
  it(2) = it(1) + 1
  months(1) = kmo + it(1) - 1
  months(2) = kmo + it(2) - 1
  if (months(1) <  1) months(1) = 12
  if (months(2) > 12) months(2) = 1

  timwt(1) = (it(1)+0.5_r8) - t
  timwt(2) = 1._r8-timwt(1)
  
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      CanopyLeafArea_col(NY,NX) = 0._r8
      StemArea_col(NY,NX)       = 0._r8
      CanopyHeight_col(NY,NX)   = 0.0_r8
      NP_col(NY,NX)=1
      DO NZ=1,NP_col(NY,NX)
        tlai_day_pft(NZ,NY,NX)     = timwt(1)*tlai_mon_pft(months(1),NZ,NY,NX)+timwt(2)*tlai_mon_pft(months(2),NZ,NY,NX)
        tsai_day_pft(NZ,NY,NX)     = timwt(1)*tsai_mon_pft(months(1),NZ,NY,NX)+timwt(2)*tsai_mon_pft(months(2),NZ,NY,NX)
        CanopyHeight_pft(NZ,NY,NX) = timwt(1)*height_top_mon_pft(months(1),NZ,NY,NX)+timwt(2)*height_top_mon_pft(months(2),NZ,NY,NX)
        CanopyLeafArea_col(NY,NX)  = CanopyLeafArea_col(NY,NX)+tlai_day_pft(NZ,NY,NX)
        StemArea_col(NY,NX)        = StemArea_col(NY,NX)+tsai_day_pft(NZ,NY,NX)
        CanopyHeight_col(NY,NX)    = AMAX1(CanopyHeight_col(NY,NX),CanopyHeight_pft(NZ,NY,NX))
        CanopyLeafAreaZ_pft(1:NumCanopyLayers,NZ,NY,NX)=tlai_day_pft(NZ,NY,NX)/real(NumCanopyLayers,kind=r8)
        CanopyStemAreaZ_pft(1:NumCanopyLayers,NZ,NY,NX)=tsai_day_pft(NZ,NY,NX)/real(NumCanopyLayers,kind=r8)

      ENDDO

      !set vertical desitribution of LAI and             
      CanopyLeafAareZ_col(1:NumCanopyLayers,NY,NX)=CanopyLeafArea_col(NY,NX)/real(NumCanopyLayers,kind=r8)
      CanopyStemAareZ_col(1:NumCanopyLayers,NY,NX)=StemArea_col(NY,NX)/real(NumCanopyLayers,kind=r8)

      !divide canopy height      
      CanopyHeightZ_col(NumCanopyLayers,NY,NX) = CanopyHeight_col(NY,NX)+0.01_r8
      ZL1(NumCanopyLayers)               = CanopyHeightZ_col(NumCanopyLayers,NY,NX)
      ZL1(0)                                = 0.0_r8
      
      !for simplicity, right now unifrom division is used
      !divide total are into NumCanopyLayers1, from top to bottom
      AreaInterval=(CanopyLeafArea_col(NY,NX)+StemArea_col(NY,NX))/NumCanopyLayers  
      DZL=CanopyHeightZ_col(NumCanopyLayers,NY,NX)/real(NumCanopyLayers,kind=r8)

      IF(AreaInterval.GT.ZEROS(NY,NX))THEN
         DO L=NumCanopyLayers,1,-1
           CanopyHeightZ_col(L-1,NY,NX)=CanopyHeightZ_col(L,NY,NX)-DZL
         ENDDO
      ENDIF

      !temporary set TEST values, assuming trees      
      NZ=1
      PlantPopulation_pft(NZ,NY,NX)=0.6_r8
      !
      call SetRootProfileZ(irootType,NL_col(NY,NX),CumDepz2LayBottom_vr(1:NL_col(NY,NX),NY,NX),PerPlantRootC_vr(1:NL_col(NY,NX)),PerPlantRootLen_vr(1:NL_col(NY,NX)))
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        RootLenDensPerPlant_pvr(ipltroot,L,NZ,NY,NX) = PerPlantRootLen_vr(L)/DLYR_3D(3,L,NY,NX)
        PopuRootMycoC_pvr(ipltroot,L,NZ,NY,NX)       = PerPlantRootC_vr(L)*PlantPopulation_pft(NZ,NY,NX)
      ENDDO
    ENDDO
  ENDDO  




  end subroutine PrescribePhenologyInterp
!------------------------------------------------------------------------------------------     
  subroutine SetRootProfileZ(irootType,NL,cdepthz,PerPlantRootC_vr,PerPlantRootLen_vr)
  !|Biome | Total fine root biomass (kg m$^{-2}$) | Live fine root biomass (kg m$^{-2}$) | Live fine root length (km m$^{-2}$) | Live fine root area index (m2 m$^{-2}$) |$\beta$|
  !|----------|----------|----------|----------|----------|----------|
  !|1. Boreal forest|0.60 | 0.23 | 2.6 | 4.6 | 0.943|
  !|2. Desert | 0.27|0.13|4.0 | 5.5 | 0.970|
  !|3. Sclerophyllous shrubs and tress |0.52|0.28|8.4|11.6|0.950|
  !|4. Temperate confierous forest|0.82|0.50|6.1|11.0|0.980|
  !|5. Temperature deciduous forest|0.78|0.44|5.4|9.8|0.967|
  !|6. Temperate grassland|1.51|0.95|112|79.1|0.943|
  !|7. Tropical deciduous forest|0.57|0.28|3.5|6.3|0.982|
  !|8. Tropical evergreen forest|0.57|0.33|4.1|7.4|0.972|
  !|9. Tropical grassland/savanna|0.99|0.51|60.4|42.5|0.972|
  !|10. Tundra|0.96|0.34|7.4|5.2|0.909|  
  implicit none
  integer, intent(in) :: iRootType
  integer, intent(in) :: NL
  real(r8),intent(in) :: cdepthz(1:NL)        !cumulative depth [m]
  real(r8),intent(out) :: PerPlantRootC_vr(1:NL)   !fine root C in each layer [gC m-2]
  real(r8),intent(out) :: PerPlantRootLen_vr(1:NL)    !root length in each layer [gC m-2]
  real(r8), parameter :: beta(10)=real((/0.943,0.970,0.950,0.980,0.967,0.943,0.982,0.972,0.972,0.909/),kind=r8)
  real(r8), parameter :: totfrootC(10)=(/0.6,0.27,0.52,0.82,0.78,1.51,0.57,0.57,0.99,0.96/)*0.488_r8 !total fine root C
  real(r8), parameter :: frootLen(10)=(/2.6,4.0,8.4,6.1,5.4,112.,3.5,4.1,60.4,7.4/)*1.e3_r8  !total fine root length
  real(r8), parameter :: PltPopDef(10)=(/0.6,1.,1.0,0.6,0.6,40.,0.6,0.6,40.,40./) !default plant population [1/m2]
  integer :: L 
  real(r8) :: CumRootFrac_vr(1:NL)  !Cumfraction of root in soil layers 
  real(r8) :: RootFrac_vr(1:NL)     !fraction of root in soil layers 

  CumRootFrac_vr(1)=1._r8-beta(irootType)**(cdepthz(1)*100._r8)
  RootFrac_vr(1)=CumRootFrac_vr(1)
  DO L=2,NL
    CumRootFrac_vr(L)=1._r8-beta(irootType)**(cdepthz(L)*100._r8)
!    if(CumRootFrac_vr(L-1)<0.99_r8)
      RootFrac_vr(L) = CumRootFrac_vr(L)-CumRootFrac_vr(L-1)
!    else
!      RootFrac_vr(L) = 1._r8-CumRootFrac_vr(L-1)
!    endif
  enddo
  
  DO L=1,NL
    PerPlantRootLen_vr(L)=frootLen(irootType)*RootFrac_vr(L)/PltPopDef(irootType)
    PerPlantRootC_vr(L)=RootFrac_vr(L)*totfrootC(irootType)/PltPopDef(irootType)
  ENDDO

  end subroutine SetRootProfileZ
end module PrescribePhenolMod
