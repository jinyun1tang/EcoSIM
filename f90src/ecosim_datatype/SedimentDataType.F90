module SedimentDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc,ndbiomcp=>ndbiomcpc,jsken=>jskenc
implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

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
  real(r8),target,allocatable ::  XN4ER(:,:,:,:)                    !total adsorbed NH4 erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XNBER(:,:,:,:)                    !total adsorbed NH4 erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XH1PEB(:,:,:,:)                   !total adsorbed HPO4  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XH2PEB(:,:,:,:)                   !total adsorbed H2PO4  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XCECER(:,:,:,:)                   !total CEC erosion , [cmol kg-1 d-2 h-1]
  real(r8),target,allocatable ::  XAECER(:,:,:,:)                   !total AEC erosion , [cmol kg-1 d-2 h-1]
  real(r8),target,allocatable ::  PALOER(:,:,:,:)                   !total adsorbed ALOH3  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XHYER(:,:,:,:)                    !total adsorbed H  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XALER(:,:,:,:)                    !total adsorbed Al  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XCAER(:,:,:,:)                    !total adsorbed Ca  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XMGER(:,:,:,:)                    !total adsorbed Mg  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XNAER(:,:,:,:)                    !total adsorbed Na  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XKAER(:,:,:,:)                    !total adsorbed K  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XHCER(:,:,:,:)                    !total adsorbed COOH  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XAL2ER(:,:,:,:)                   !total adsorbed AlOH2  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  XOH0ER(:,:,:,:)                   !total adsorbed OH-  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XOH1ER(:,:,:,:)                   !total adsorbed OH  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XOH2ER(:,:,:,:)                   !total adsorbed OH2  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XH1PER(:,:,:,:)                   !total adsorbed HPO4  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XH2PER(:,:,:,:)                   !total adsorbed H2PO4  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  XOH0EB(:,:,:,:)                   !total adsorbed OH-  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XOH1EB(:,:,:,:)                   !total adsorbed OH  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XOH2EB(:,:,:,:)                   !total adsorbed OH2  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  PFEOER(:,:,:,:)                   !total adsorbed FeOH3  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  PCACER(:,:,:,:)                   !total adsorbed CaCO3  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  PCASER(:,:,:,:)                   !total adsorbed CaSO4  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  PALPER(:,:,:,:)                   !total adsorbed AlPO4  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  PFEPER(:,:,:,:)                   !total adsorbed FePO4  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  PCPDER(:,:,:,:)                   !total adsorbed CaHPO4  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  PCPHER(:,:,:,:)                   !total adsorbed apatite  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  PCPMER(:,:,:,:)                   !total adsorbed CaH2PO4  erosion non-band , [g d-2 h-1]
  real(r8),target,allocatable ::  PALPEB(:,:,:,:)                   !total adsorbed AlPO4  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  PFEPEB(:,:,:,:)                   !total adsorbed FePO4  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  PCPDEB(:,:,:,:)                   !total adsorbed CaHPO4  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  PCPHEB(:,:,:,:)                   !total adsorbed apatite  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  PCPMEB(:,:,:,:)                   !total adsorbed CaH2PO4  erosion band , [g d-2 h-1]
  real(r8),target,allocatable ::  XFE2ER(:,:,:,:)                   !sedimentation rate of Fe(OH)2
  real(r8),target,allocatable ::  XSEDER(:,:,:,:)                   !sediment erosion, [Mg d-2 h-1]
  real(r8),target,allocatable ::  XFEER(:,:,:,:)                    !erosion loss rate of Fe
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
  allocate(XN4ER(2,2,JV,JH));   XN4ER=0._r8
  allocate(XNBER(2,2,JV,JH));   XNBER=0._r8
  allocate(XH1PEB(2,2,JV,JH));  XH1PEB=0._r8
  allocate(XH2PEB(2,2,JV,JH));  XH2PEB=0._r8
  allocate(XCECER(2,2,JV,JH));  XCECER=0._r8
  allocate(XAECER(2,2,JV,JH));  XAECER=0._r8
  allocate(PALOER(2,2,JV,JH));  PALOER=0._r8
  allocate(XHYER(2,2,JV,JH));   XHYER=0._r8
  allocate(XALER(2,2,JV,JH));   XALER=0._r8
  allocate(XCAER(2,2,JV,JH));   XCAER=0._r8
  allocate(XMGER(2,2,JV,JH));   XMGER=0._r8
  allocate(XNAER(2,2,JV,JH));   XNAER=0._r8
  allocate(XKAER(2,2,JV,JH));   XKAER=0._r8
  allocate(XHCER(2,2,JV,JH));   XHCER=0._r8
  allocate(XAL2ER(2,2,JV,JH));  XAL2ER=0._r8
  allocate(XOH0ER(2,2,JV,JH));  XOH0ER=0._r8
  allocate(XOH1ER(2,2,JV,JH));  XOH1ER=0._r8
  allocate(XOH2ER(2,2,JV,JH));  XOH2ER=0._r8
  allocate(XH1PER(2,2,JV,JH));  XH1PER=0._r8
  allocate(XH2PER(2,2,JV,JH));  XH2PER=0._r8
  allocate(XOH0EB(2,2,JV,JH));  XOH0EB=0._r8
  allocate(XOH1EB(2,2,JV,JH));  XOH1EB=0._r8
  allocate(XOH2EB(2,2,JV,JH));  XOH2EB=0._r8
  allocate(PFEOER(2,2,JV,JH));  PFEOER=0._r8
  allocate(PCACER(2,2,JV,JH));  PCACER=0._r8
  allocate(PCASER(2,2,JV,JH));  PCASER=0._r8
  allocate(PALPER(2,2,JV,JH));  PALPER=0._r8
  allocate(PFEPER(2,2,JV,JH));  PFEPER=0._r8
  allocate(PCPDER(2,2,JV,JH));  PCPDER=0._r8
  allocate(PCPHER(2,2,JV,JH));  PCPHER=0._r8
  allocate(PCPMER(2,2,JV,JH));  PCPMER=0._r8
  allocate(PALPEB(2,2,JV,JH));  PALPEB=0._r8
  allocate(PFEPEB(2,2,JV,JH));  PFEPEB=0._r8
  allocate(PCPDEB(2,2,JV,JH));  PCPDEB=0._r8
  allocate(PCPHEB(2,2,JV,JH));  PCPHEB=0._r8
  allocate(PCPMEB(2,2,JV,JH));  PCPMEB=0._r8
  allocate(XFE2ER(2,2,JV,JH));  XFE2ER=0._r8
  allocate(XSEDER(2,2,JV,JH));  XSEDER=0._r8
  allocate(XFEER(2,2,JV,JH));   XFEER=0._r8
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
  if (allocated(TSED))     deallocate(TSED)
  if (allocated(DETS))     deallocate(DETS)
  if (allocated(DETE))     deallocate(DETE)
  if (allocated(CER))      deallocate(CER)
  if (allocated(XER))      deallocate(XER)
  if (allocated(PTDSNU))   deallocate(PTDSNU)
  if (allocated(VLS))      deallocate(VLS)
  if (allocated(SED))      deallocate(SED)
  if (allocated(XSANER))   deallocate(XSANER)
  if (allocated(XSILER))   deallocate(XSILER)
  if (allocated(XCLAER))   deallocate(XCLAER)
  if (allocated(XNH4ER))   deallocate(XNH4ER)
  if (allocated(XNH3ER))   deallocate(XNH3ER)
  if (allocated(XNHUER))   deallocate(XNHUER)
  if (allocated(XNO3ER))   deallocate(XNO3ER)
  if (allocated(XNH4EB))   deallocate(XNH4EB)
  if (allocated(XNH3EB))   deallocate(XNH3EB)
  if (allocated(XNHUEB))   deallocate(XNHUEB)
  if (allocated(XNO3EB))   deallocate(XNO3EB)
  if (allocated(XN4ER))    deallocate(XN4ER)
  if (allocated(XNBER))    deallocate(XNBER)
  if (allocated(XH1PEB))   deallocate(XH1PEB)
  if (allocated(XH2PEB))   deallocate(XH2PEB)
  if (allocated(XCECER))   deallocate(XCECER)
  if (allocated(XAECER))   deallocate(XAECER)
  if (allocated(PALOER))   deallocate(PALOER)
  if (allocated(XHYER))    deallocate(XHYER)
  if (allocated(XALER))    deallocate(XALER)
  if (allocated(XCAER))    deallocate(XCAER)
  if (allocated(XMGER))    deallocate(XMGER)
  if (allocated(XNAER))    deallocate(XNAER)
  if (allocated(XKAER))    deallocate(XKAER)
  if (allocated(XHCER))    deallocate(XHCER)
  if (allocated(XAL2ER))   deallocate(XAL2ER)
  if (allocated(XOH0ER))   deallocate(XOH0ER)
  if (allocated(XOH1ER))   deallocate(XOH1ER)
  if (allocated(XOH2ER))   deallocate(XOH2ER)
  if (allocated(XH1PER))   deallocate(XH1PER)
  if (allocated(XH2PER))   deallocate(XH2PER)
  if (allocated(XOH0EB))   deallocate(XOH0EB)
  if (allocated(XOH1EB))   deallocate(XOH1EB)
  if (allocated(XOH2EB))   deallocate(XOH2EB)
  if (allocated(PFEOER))   deallocate(PFEOER)
  if (allocated(PCACER))   deallocate(PCACER)
  if (allocated(PCASER))   deallocate(PCASER)
  if (allocated(PALPER))   deallocate(PALPER)
  if (allocated(PFEPER))   deallocate(PFEPER)
  if (allocated(PCPDER))   deallocate(PCPDER)
  if (allocated(PCPHER))   deallocate(PCPHER)
  if (allocated(PCPMER))   deallocate(PCPMER)
  if (allocated(PALPEB))   deallocate(PALPEB)
  if (allocated(PFEPEB))   deallocate(PFEPEB)
  if (allocated(PCPDEB))   deallocate(PCPDEB)
  if (allocated(PCPHEB))   deallocate(PCPHEB)
  if (allocated(PCPMEB))   deallocate(PCPMEB)
  if (allocated(XFE2ER))   deallocate(XFE2ER)
  if (allocated(XSEDER))   deallocate(XSEDER)
  if (allocated(XFEER))    deallocate(XFEER)
  if (allocated(ORCER))    deallocate(ORCER)
  if (allocated(ORNER))    deallocate(ORNER)
  if (allocated(ORPER))    deallocate(ORPER)
  if (allocated(OHCER))    deallocate(OHCER)
  if (allocated(OHNER))    deallocate(OHNER)
  if (allocated(OHPER))    deallocate(OHPER)
  if (allocated(OHAER))    deallocate(OHAER)
  if (allocated(OSCER))    deallocate(OSCER)
  if (allocated(OSAER))    deallocate(OSAER)
  if (allocated(OSNER))    deallocate(OSNER)
  if (allocated(OSPER))    deallocate(OSPER)
  end subroutine DestructSedimentData



end module SedimentDataType
