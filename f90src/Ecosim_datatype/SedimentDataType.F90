module SedimentDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  TSED(:,:)                         !erosion rate [Mg d-2 t-1]
  real(r8),target,allocatable ::  SoilDetachability4Erosion1(:,:)                         !soil detachment, [-]
  real(r8),target,allocatable ::  SoilDetachability4Erosion2(:,:)                         !soil detachability
  real(r8),target,allocatable ::  CER(:,:)                          !soil detachment/deposition, [-]
  real(r8),target,allocatable ::  XER(:,:)                          !soil detachment/deposition, [-]
  real(r8),target,allocatable ::  ParticleDensitySurfLay(:,:)       !particle density of surface layer
  real(r8),target,allocatable ::  VLS(:,:)                          !hourly sinking rate
  real(r8),target,allocatable ::  SED(:,:)                          !sediment transport, [Mg d-2 h-1]
  real(r8),target,allocatable ::  XSANER(:,:,:,:)                   !total sand erosion , [Mg d-2 h-1]
  real(r8),target,allocatable ::  XSILER(:,:,:,:)                   !total silt erosion , [Mg d-2 h-1]
  real(r8),target,allocatable ::  XCLAER(:,:,:,:)                   !total clay erosion , [Mg d-2 h-1]
  real(r8),target,allocatable ::  XNH4ER(:,:,:,:)                   !total NH4 fertilizer erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XNH3ER(:,:,:,:)                   !total NH3 fertilizer erosion non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNHUER(:,:,:,:)                   !total urea fertilizer erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XNO3ER(:,:,:,:)                   !total NO3 fertilizer erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XNH4EB(:,:,:,:)                   !total NH4 fertilizer erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XNH3EB(:,:,:,:)                   !total NH3 fertilizer erosion band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNHUEB(:,:,:,:)                   !total urea fertilizer erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XNO3EB(:,:,:,:)                   !total NO3 fertilizer erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  trcx_XER(:,:,:,:,:)               !total adsorbed sediment erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  trcp_ER(:,:,:,:,:)                !total adsorbed ALOH3  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  cumSedErosion(:,:,:,:)            !sediment erosion, [Mg d-2 h-1]
  real(r8),target,allocatable ::  ORMER(:,:,:,:,:,:,:)                !microbial residue C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OHMER(:,:,:,:,:,:)                  !adsorbed C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OSMER(:,:,:,:,:,:,:)                !humus C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OSAER(:,:,:,:,:,:)                !colonized humus C  erosion , [g d-2 h-1]
  private :: InitAllocate

  contains
!----------------------------------------------------------------------

  subroutine InitSedimentData

  implicit none

  call InitAllocate
  end subroutine InitSedimentData

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(TSED(JY,JX));        TSED=0._r8
  allocate(SoilDetachability4Erosion1(JY,JX));        SoilDetachability4Erosion1=0._r8
  allocate(SoilDetachability4Erosion2(JY,JX));        SoilDetachability4Erosion2=0._r8
  allocate(CER(JY,JX));         CER=0._r8
  allocate(XER(JY,JX));         XER=0._r8
  allocate(ParticleDensitySurfLay(JY,JX));      ParticleDensitySurfLay=0._r8
  allocate(VLS(JY,JX));         VLS=0._r8
  allocate(SED(JY,JX));         SED=0._r8
  allocate(XSANER(2,2,JV,JH));  XSANER=0._r8
  allocate(XSILER(2,2,JV,JH));  XSILER=0._r8
  allocate(XCLAER(2,2,JV,JH));  XCLAER=0._r8
  allocate(XNH4ER(2,2,JV,JH));  XNH4ER=0._r8
  allocate(XNH3ER(2,2,JV,JH));  XNH3ER=0._r8
  allocate(XNHUER(2,2,JV,JH));  XNHUER=0._r8
  allocate(XNO3ER(2,2,JV,JH));  XNO3ER=0._r8
  allocate(XNH4EB(2,2,JV,JH));  XNH4EB=0._r8
  allocate(XNH3EB(2,2,JV,JH));  XNH3EB=0._r8
  allocate(XNHUEB(2,2,JV,JH));  XNHUEB=0._r8
  allocate(XNO3EB(2,2,JV,JH));  XNO3EB=0._r8
  allocate(trcx_XER(idx_beg:idx_end,2,2,JV,JH));   trcx_XER=0._r8
  allocate(trcp_ER(idsp_beg:idsp_end,2,2,JV,JH));  trcp_ER=0._r8
  allocate(cumSedErosion(2,2,JV,JH));  cumSedErosion=0._r8
  allocate(ORMER(NumPlantChemElms,ndbiomcp,1:jcplx,2,2,JV,JH));ORMER=0._r8
  allocate(OHMER(idom_beg:idom_end,1:jcplx,2,2,JV,JH));OHMER=0._r8
  allocate(OSMER(NumPlantChemElms,jsken,1:jcplx,2,2,JV,JH));OSMER=0._r8
  allocate(OSAER(jsken,1:jcplx,2,2,JV,JH));OSAER=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSedimentData

  use abortutils, only : destroy
  implicit none

  call destroy(trcp_ER)
  call destroy(trcx_XER)
  call destroy(TSED)
  call destroy(SoilDetachability4Erosion1)
  call destroy(SoilDetachability4Erosion2)
  call destroy(CER)
  call destroy(XER)
  call destroy(ParticleDensitySurfLay)
  call destroy(VLS)
  call destroy(SED)
  call destroy(XSANER)
  call destroy(XSILER)
  call destroy(XCLAER)
  call destroy(XNH4ER)
  call destroy(XNH3ER)
  call destroy(XNHUER)
  call destroy(XNO3ER)
  call destroy(XNH4EB)
  call destroy(XNH3EB)
  call destroy(XNHUEB)
  call destroy(XNO3EB)
  call destroy(cumSedErosion)
  call destroy(ORMER)
  call destroy(OHMER)
  call destroy(OSMER)
  call destroy(OSAER)
  end subroutine DestructSedimentData



end module SedimentDataType
