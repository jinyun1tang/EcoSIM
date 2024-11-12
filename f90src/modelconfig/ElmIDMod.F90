module ElmIDMod
!
! DESCRIPTION:
!  Chemical element ids
implicit none
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
  integer, parameter :: ielmc=1    !carbon element
  integer, parameter :: ielmn=2    !nitrogen element
  integer, parameter :: ielmp=3    !phosphorus element
  integer, parameter :: NumPlantChemElms=3   !totally three elements
! erosion model options
  integer, parameter :: ieros_noaction       = -1
  integer, parameter :: ieros_frzthawelv     = 0
  integer, parameter :: ieros_frzthaweros    = 1
  integer, parameter :: ieros_frzthawsom     = 2
  integer, parameter :: ieros_frzthawsomeros = 3

! water flux direction
  integer, parameter :: iEastWestDirection=1   !east-west direction
  integer, parameter :: iNorthSouthDirection=2   !north-south direction
  integer, parameter :: iVerticalDirection=3   !vertical direction

! soil properties
  integer, parameter :: isoi_fc    = 1   !field capacity
  integer, parameter :: isoi_wp    = 2   !wilting point
  integer, parameter :: isoi_scnv  = 3  !vertical hydraulic conductivity
  integer, parameter :: isoi_scnh  = 4  !horizontal hydraulic conductivity
  integer, parameter :: isoi_set   = 0
  integer, parameter :: isoi_unset = 1
! plant harvest
  integer, parameter :: iplthvst_leaf=1 !leaf
  integer, parameter :: iplthvst_finenonleaf=2 !fine non-leaf
  integer, parameter :: iplthvst_woody=3 !woody
  integer, parameter :: iplthvst_stdead=4 !standing dead
! photosynthesis
  integer, parameter :: ic4_photo=4
  integer, parameter :: ic3_photo=3
! fertilizer
  integer, parameter :: ifert_nh4               = 1
  integer, parameter :: ifert_nh3               = 2
  integer, parameter :: ifert_urea              = 3
  integer, parameter :: ifert_no3               = 4
  integer, parameter :: ifert_nh4_band          = 5
  integer, parameter :: ifert_nh3_band          = 6
  integer, parameter :: ifert_urea_band         = 7
  integer, parameter :: ifert_no3_band          = 8
  integer, parameter :: ifert_P_Ca_H2PO4_2      = 9
  integer, parameter :: ifert_P_Ca_H2PO4_2_band = 10
  integer, parameter :: ifert_P_apatite         = 11
  integer, parameter :: ifert_Ca_lime           = 12
  integer, parameter :: ifert_Ca_gypsum         = 13
  integer, parameter :: ifert_plant_resC        = 14
  integer, parameter :: ifert_plant_resN        = 15
  integer, parameter :: ifert_plant_resP        = 16
  integer, parameter :: ifert_plant_manuC       = 17
  integer, parameter :: ifert_plant_manuN       = 18
  integer, parameter :: ifert_plant_manuP       = 19
! root order
  integer, parameter :: iroot_1st =1
  integer, parameter :: iroot_2nd =2
! root profile
  integer, parameter :: irootp_shallow=0
! mycorrhizae
  integer, parameter :: imycorr_none=1
  integer, parameter :: imycorr_arbu=2   !arbuscular

  integer, parameter :: ipltroot=1
  integer, parameter :: imycorrhz=2

  integer, parameter :: iliving_branch=0

  integer, parameter :: itrue=1
  integer, parameter :: ifalse=0

  integer, parameter :: iEnable  = 0
  integer, parameter :: iDisable = 1

  integer, parameter :: iLive=0  !plant organ alive
  integer, parameter :: iDead=1  !plant organ dead

  integer, parameter :: iActive  = 1  !plant active
  integer, parameter :: iDormant = 0  !plant dormant/dead

  integer, parameter :: jharvtyp_noaction  = 0
  integer, parameter :: jharvtyp_terminate = 1
  integer, parameter :: jharvtyp_tmareseed = 2

  integer, parameter :: ihav_pft = 1
  integer, parameter :: ihav_eco = 2
  integer, parameter :: iharvtyp_none    = 0
  integer, parameter :: iharvtyp_grain   = 1
  integer, parameter :: iharvtyp_allabv  = 2
  integer, parameter :: iharvtyp_pruning = 3
  integer, parameter :: iharvtyp_grazing = 4
  integer, parameter :: iharvtyp_fire    = 5
  integer, parameter :: iharvtyp_herbivo = 6

  integer, parameter :: iplt_annual    = 0
  integer, parameter :: iplt_perennial = 1

  integer, parameter :: iplt_bryophyte=0
  integer, parameter :: iplt_grasslike=1
  integer, parameter :: iplt_treelike =2
  
  integer, parameter :: igraintyp_abvgrnd=0
  integer, parameter :: igraintyp_blwgrnd=1

  integer, parameter :: ideterminate          = 0
  integer, parameter :: inondeterminate       = 1
  integer, parameter :: ipltcal_Planting      = 0
  integer, parameter :: ipltcal_Emerge        = 1
  integer, parameter :: ipltcal_InitFloral    = 2
  integer, parameter :: ipltcal_Jointing      = 3
  integer, parameter :: ipltcal_Elongation    = 4
  integer, parameter :: ipltcal_Heading       = 5
  integer, parameter :: ipltcal_Anthesis      = 6
  integer, parameter :: ipltcal_BeginSeedFill = 7
  integer, parameter :: ipltcal_SetSeedNumber = 8
  integer, parameter :: ipltcal_SetSeedMass   = 9
  integer, parameter :: ipltcal_EndSeedFill   = 10

  integer, parameter :: iphotop_neutral=0
  integer, parameter :: iphotop_short  =1
  integer, parameter :: iphotop_long   =2

  integer, parameter :: iphenotyp_evgreen       = 0
  integer, parameter :: iphenotyp_coldecid      = 1
  integer, parameter :: iphenotyp_drouhtdecidu  = 2
  integer, parameter :: iphenotyp_coldroutdecid = 3

  integer, parameter :: in2fixtyp_none          = 0
  integer, parameter :: in2fixtyp_root_fast     = 1
  integer, parameter :: in2fixtyp_root_medium   = 2
  integer, parameter :: in2fixtyp_root_slow     = 3
  integer, parameter :: in2fixtyp_canopy_fast   = 4
  integer, parameter :: in2fixtyp_canopy_medium = 5
  integer, parameter :: in2fixtyp_canopy_slow   = 6

  integer, parameter :: ithermozone_arcboreal = 1
  integer, parameter :: ithermozone_cooltempr = 2
  integer, parameter :: ithermozone_warmtempr = 3
  integer, parameter :: ithermozone_subtropic = 4
  integer, parameter :: ithermozone_tropical  = 5
end module ElmIDMod
