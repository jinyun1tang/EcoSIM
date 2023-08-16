module SedimentDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,ndbiomcp=>ndbiomcpc,jsken=>jskenc
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  TSED(:,:)                         !erosion rate [Mg d-2 t-1]
  real(r8),target,allocatable ::  DETS(:,:)                         !soil detachment, [-]
  real(r8),target,allocatable ::  DETE(:,:)                         !soil detachability
  real(r8),target,allocatable ::  CER(:,:)                          !soil detachment/deposition, [-]
  real(r8),target,allocatable ::  XER(:,:)                          !soil detachment/deposition, [-]
  real(r8),target,allocatable ::  PTDSNU(:,:)                       !particle density of surface layer
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
  real(r8),target,allocatable ::  trcp_ER(:,:,:,:,:)                   !total adsorbed ALOH3  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XSEDER(:,:,:,:)                   !sediment erosion, [Mg d-2 h-1]
  real(r8),target,allocatable ::  ORCER(:,:,:,:,:,:)                !microbial residue C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  ORNER(:,:,:,:,:,:)                !microbial residue N  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  ORPER(:,:,:,:,:,:)                !microbial residue P  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OHCER(:,:,:,:,:)                  !adsorbed C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OHNER(:,:,:,:,:)                  !adsorbed N  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OHPER(:,:,:,:,:)                  !adsorbed P  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OHAER(:,:,:,:,:)                  !adsorbed acetate  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OSCER(:,:,:,:,:,:)                !humus C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OSAER(:,:,:,:,:,:)                !colonized humus C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OSNER(:,:,:,:,:,:)                !humus N  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  OSPER(:,:,:,:,:,:)                !humus P  erosion , [g d-2 h-1]
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
  allocate(DETS(JY,JX));        DETS=0._r8
  allocate(DETE(JY,JX));        DETE=0._r8
  allocate(CER(JY,JX));         CER=0._r8
  allocate(XER(JY,JX));         XER=0._r8
  allocate(PTDSNU(JY,JX));      PTDSNU=0._r8
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
  allocate(XSEDER(2,2,JV,JH));  XSEDER=0._r8
  allocate(ORCER(ndbiomcp,1:jcplx,2,2,JV,JH));ORCER=0._r8
  allocate(ORNER(ndbiomcp,1:jcplx,2,2,JV,JH));ORNER=0._r8
  allocate(ORPER(ndbiomcp,1:jcplx,2,2,JV,JH));ORPER=0._r8
  allocate(OHCER(1:jcplx,2,2,JV,JH));OHCER=0._r8
  allocate(OHNER(1:jcplx,2,2,JV,JH));OHNER=0._r8
  allocate(OHPER(1:jcplx,2,2,JV,JH));OHPER=0._r8
  allocate(OHAER(1:jcplx,2,2,JV,JH));OHAER=0._r8
  allocate(OSCER(jsken,1:jcplx,2,2,JV,JH));OSCER=0._r8
  allocate(OSAER(jsken,1:jcplx,2,2,JV,JH));OSAER=0._r8
  allocate(OSNER(jsken,1:jcplx,2,2,JV,JH));OSNER=0._r8
  allocate(OSPER(jsken,1:jcplx,2,2,JV,JH));OSPER=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSedimentData

  use abortutils, only : destroy
  implicit none

  call destroy(trcp_ER)
  call destroy(trcx_XER)
  call destroy(TSED)
  call destroy(DETS)
  call destroy(DETE)
  call destroy(CER)
  call destroy(XER)
  call destroy(PTDSNU)
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
  call destroy(XSEDER)
  call destroy(ORCER)
  call destroy(ORNER)
  call destroy(ORPER)
  call destroy(OHCER)
  call destroy(OHNER)
  call destroy(OHPER)
  call destroy(OHAER)
  call destroy(OSCER)
  call destroy(OSAER)
  call destroy(OSNER)
  call destroy(OSPER)
  end subroutine DestructSedimentData



end module SedimentDataType
