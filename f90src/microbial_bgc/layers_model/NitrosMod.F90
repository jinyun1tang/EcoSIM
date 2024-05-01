module nitrosMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils  , only : endrun
  use minimathmod, only : safe_adb,AZMAX1
  use EcoSiMParDataMod, only : micpar
  use MicrobialDataType
  use NitroPars
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use NitroDiagTypes
  use NitroDisturbMod
  use GridConsts
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use SoilPhysDataType
  use SoilPropertyDataType
  use GridDataType
  use MicBGCMod
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

!
  public :: initNitro
  public :: VerticalLitterMixLvsLL
  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro

  implicit none

  call initNitro1Layer

  end subroutine initNitro

!------------------------------------------------------------------------------------------

  subroutine VerticalLitterMixLvsLL(I,J,L,NY,NX)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: FOSCXS,FOSCXD
  real(r8) :: ORGRL,ORGRLL
  real(r8) :: OSCXD
  integer :: LL,LN
!     begin_execution

  IF(FOSCZ0.GT.ZERO)THEN
!     ORGR=total litter C
!     FOSCZ0=rate constant for mixing surface litter
!     FOSCXS=mixing fraction for surface litter
!     TOQCK=total active biomass respiration activity
!     TFNX=temperature function
!     VLSoilPoreMicP_vr=soil layer volume
!     OSCXD=mixing required for equilibrating litter concentration
!     FOSCXD=mixing fraction for equilibrating subsurface litter
!     FOSCXS=mixing fraction for subsurface litter
!
    IF(L.LT.NL(NY,NX))THEN
!     get mixing rate
      IF(L.EQ.0)THEN
        LL=NU(NY,NX)
        IF(ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXS=AMIN1(1.0_r8,FOSCZ0/ORGR(L,NY,NX)*TOQCK(L,NY,NX))
        ELSE
          FOSCXS=0.0_r8
        ENDIF
      ELSE
        D1100: DO LN=L+1,NL(NY,NX)
          IF(VLSoilPoreMicP_vr(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
            LL=LN
            exit
          ENDIF
        ENDDO D1100
        ORGRL=AZMAX1(ORGR(L,NY,NX))
        ORGRLL=AZMAX1(ORGR(LL,NY,NX))
        OSCXD=(ORGRL*VGeomLayer(LL,NY,NX)-ORGRLL*VGeomLayer(L,NY,NX))/(VGeomLayer(L,NY,NX)+VGeomLayer(LL,NY,NX))
        IF(OSCXD.GT.0.0_r8.AND.ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXD=OSCXD/ORGR(L,NY,NX)
        ELSEIF(OSCXD.LT.0.0_r8.AND.ORGR(LL,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXD=OSCXD/ORGR(LL,NY,NX)
        ELSE
          FOSCXD=0.0_r8
        ENDIF
        IF(VGeomLayer(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FOSCXS=FOSCZL*FOSCXD*TOQCK(L,NY,NX)/VGeomLayer(L,NY,NX)
        ELSE
          FOSCXS=0.0_r8
        ENDIF
      ENDIF

!     apply mixing
      call ApplyVerticalMix(FOSCXS,L,LL,NY,NX)
    ENDIF

  ENDIF
  end subroutine VerticalLitterMixLvsLL
!------------------------------------------------------------------------------------------

  subroutine ApplyVerticalMix(FOSCXS,L,LL,NY,NX)

  implicit none
  real(r8), intent(in) :: FOSCXS
  integer, intent(in) :: L,LL,NY,NX

  real(r8) :: OMCXS,OMNXS,OMPXS
  real(r8) :: ORCXS,ORNXS,ORPXS
  real(r8) :: OQCXS,OQCHXS,OHCXS,OQAXS
  real(r8) :: OQAHXS,OHAXS,OQNXS,OQNHXS
  real(r8) :: OHNXS,OQPXS,OQPHXS,OHPXS
  real(r8) :: OSCXS,OSAXS,OSNXS,OSPXS
  integer :: K,M,N,NGL,MID
  
!     begin_execution
  IF(FOSCXS.GT.ZERO)THEN
    D7971: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      D7961: DO N=1,NumMicbFunGroups
        DO NGL=JGnio(N),JGnfo(N)
          D7962: DO M=1,micpar%nlbiomcp
            MID=micpar%get_micb_id(M,NGL)
            IF(FOSCXS.GT.0.0)THEN
              OMCXS=FOSCXS*AZMAX1(OMEheter(ielmc,MID,K,L,NY,NX))
              OMNXS=FOSCXS*AZMAX1(OMEheter(ielmn,MID,K,L,NY,NX))
              OMPXS=FOSCXS*AZMAX1(OMEheter(ielmp,MID,K,L,NY,NX))
            ELSE
              OMCXS=FOSCXS*AZMAX1(OMEheter(ielmc,MID,K,LL,NY,NX))
              OMNXS=FOSCXS*AZMAX1(OMEheter(ielmn,MID,K,LL,NY,NX))
              OMPXS=FOSCXS*AZMAX1(OMEheter(ielmp,MID,K,LL,NY,NX))
            ENDIF
            OMEheter(ielmc,MID,K,L,NY,NX)=OMEheter(ielmc,MID,K,L,NY,NX)-OMCXS
            OMEheter(ielmn,MID,K,L,NY,NX)=OMEheter(ielmn,MID,K,L,NY,NX)-OMNXS
            OMEheter(ielmp,MID,K,L,NY,NX)=OMEheter(ielmp,MID,K,L,NY,NX)-OMPXS
            OMEheter(ielmc,MID,K,LL,NY,NX)=OMEheter(ielmc,MID,K,LL,NY,NX)+OMCXS
            OMEheter(ielmn,MID,K,LL,NY,NX)=OMEheter(ielmn,MID,K,LL,NY,NX)+OMNXS
            OMEheter(ielmp,MID,K,LL,NY,NX)=OMEheter(ielmp,MID,K,LL,NY,NX)+OMPXS
          ENDDO D7962
        ENDDO
      ENDDO D7961
    ENDDO D7971

    D7901: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      D7941: DO M=1,micpar%ndbiomcp
        IF(FOSCXS.GT.0.0_r8)THEN
          ORCXS=FOSCXS*AZMAX1(OMBioResdu_vr(ielmc,M,K,L,NY,NX))
          ORNXS=FOSCXS*AZMAX1(OMBioResdu_vr(ielmn,M,K,L,NY,NX))
          ORPXS=FOSCXS*AZMAX1(OMBioResdu_vr(ielmp,M,K,L,NY,NX))
        ELSE
          ORCXS=FOSCXS*AZMAX1(OMBioResdu_vr(ielmc,M,K,LL,NY,NX))
          ORNXS=FOSCXS*AZMAX1(OMBioResdu_vr(ielmn,M,K,LL,NY,NX))
          ORPXS=FOSCXS*AZMAX1(OMBioResdu_vr(ielmp,M,K,LL,NY,NX))
        ENDIF
        OMBioResdu_vr(ielmc,M,K,L,NY,NX)=OMBioResdu_vr(ielmc,M,K,L,NY,NX)-ORCXS
        OMBioResdu_vr(ielmn,M,K,L,NY,NX)=OMBioResdu_vr(ielmn,M,K,L,NY,NX)-ORNXS
        OMBioResdu_vr(ielmp,M,K,L,NY,NX)=OMBioResdu_vr(ielmp,M,K,L,NY,NX)-ORPXS
        OMBioResdu_vr(ielmc,M,K,LL,NY,NX)=OMBioResdu_vr(ielmc,M,K,LL,NY,NX)+ORCXS
        OMBioResdu_vr(ielmn,M,K,LL,NY,NX)=OMBioResdu_vr(ielmn,M,K,LL,NY,NX)+ORNXS
        OMBioResdu_vr(ielmp,M,K,LL,NY,NX)=OMBioResdu_vr(ielmp,M,K,LL,NY,NX)+ORPXS
      ENDDO D7941
      IF(FOSCXS.GT.0.0_r8)THEN
        OQCXS=FOSCXS*AZMAX1(DOM(idom_doc,K,L,NY,NX))
        OQCHXS=FOSCXS*AZMAX1(DOM_Macp(idom_doc,K,L,NY,NX))
        OHCXS=FOSCXS*AZMAX1(SorbedOM_vr(ielmc,K,L,NY,NX))
        OQAXS=FOSCXS*AZMAX1(DOM(idom_acetate,K,L,NY,NX))
        OQAHXS=FOSCXS*AZMAX1(DOM_Macp(idom_acetate,K,L,NY,NX))
        OHAXS=FOSCXS*AZMAX1(SorbedOM_vr(idom_acetate,K,L,NY,NX))
        OQNXS=FOSCXS*AZMAX1(DOM(idom_don,K,L,NY,NX))
        OQNHXS=FOSCXS*AZMAX1(DOM_Macp(idom_don,K,L,NY,NX))
        OHNXS=FOSCXS*AZMAX1(SorbedOM_vr(ielmn,K,L,NY,NX))
        OQPXS=FOSCXS*AZMAX1(DOM(idom_dop,K,L,NY,NX))
        OQPHXS=FOSCXS*AZMAX1(DOM_Macp(idom_dop,K,L,NY,NX))
        OHPXS=FOSCXS*AZMAX1(SorbedOM_vr(ielmp,K,L,NY,NX))
      ELSE
        OQCXS=FOSCXS*AZMAX1(DOM(idom_doc,K,LL,NY,NX))
        OQCHXS=FOSCXS*AZMAX1(DOM_Macp(idom_doc,K,LL,NY,NX))
        OHCXS=FOSCXS*AZMAX1(SorbedOM_vr(ielmc,K,LL,NY,NX))
        OQAXS=FOSCXS*AZMAX1(DOM(idom_acetate,K,LL,NY,NX))
        OQAHXS=FOSCXS*AZMAX1(DOM_Macp(idom_acetate,K,LL,NY,NX))
        OHAXS=FOSCXS*AZMAX1(SorbedOM_vr(idom_acetate,K,LL,NY,NX))
        OQNXS=FOSCXS*AZMAX1(DOM(idom_don,K,LL,NY,NX))
        OQNHXS=FOSCXS*AZMAX1(DOM_Macp(idom_don,K,LL,NY,NX))
        OHNXS=FOSCXS*AZMAX1(SorbedOM_vr(ielmn,K,LL,NY,NX))
        OQPXS=FOSCXS*AZMAX1(DOM(idom_dop,K,LL,NY,NX))
        OQPHXS=FOSCXS*AZMAX1(DOM_Macp(idom_dop,K,LL,NY,NX))
        OHPXS=FOSCXS*AZMAX1(SorbedOM_vr(ielmp,K,LL,NY,NX))
      ENDIF
      DOM(idom_doc,K,L,NY,NX)=DOM(idom_doc,K,L,NY,NX)-OQCXS
      DOM_Macp(idom_doc,K,L,NY,NX)=DOM_Macp(idom_doc,K,L,NY,NX)-OQCHXS
      SorbedOM_vr(ielmc,K,L,NY,NX)=SorbedOM_vr(ielmc,K,L,NY,NX)-OHCXS
      DOM(idom_acetate,K,L,NY,NX)=DOM(idom_acetate,K,L,NY,NX)-OQAXS
      DOM_Macp(idom_acetate,K,L,NY,NX)=DOM_Macp(idom_acetate,K,L,NY,NX)-OQAHXS
      SorbedOM_vr(idom_acetate,K,L,NY,NX)=SorbedOM_vr(idom_acetate,K,L,NY,NX)-OHAXS
      DOM(idom_don,K,L,NY,NX)=DOM(idom_don,K,L,NY,NX)-OQNXS
      DOM_Macp(idom_don,K,L,NY,NX)=DOM_Macp(idom_don,K,L,NY,NX)-OQNHXS
      SorbedOM_vr(ielmn,K,L,NY,NX)=SorbedOM_vr(ielmn,K,L,NY,NX)-OHNXS
      DOM(idom_dop,K,L,NY,NX)=DOM(idom_dop,K,L,NY,NX)-OQPXS
      DOM_Macp(idom_dop,K,L,NY,NX)=DOM_Macp(idom_dop,K,L,NY,NX)-OQPHXS
      SorbedOM_vr(ielmp,K,L,NY,NX)=SorbedOM_vr(ielmp,K,L,NY,NX)-OHPXS
      DOM(idom_doc,K,LL,NY,NX)=DOM(idom_doc,K,LL,NY,NX)+OQCXS
      DOM_Macp(idom_doc,K,LL,NY,NX)=DOM_Macp(idom_doc,K,LL,NY,NX)+OQCHXS
      SorbedOM_vr(ielmc,K,LL,NY,NX)=SorbedOM_vr(ielmc,K,LL,NY,NX)+OHCXS
      DOM(idom_acetate,K,LL,NY,NX)=DOM(idom_acetate,K,LL,NY,NX)+OQAXS
      DOM_Macp(idom_acetate,K,LL,NY,NX)=DOM_Macp(idom_acetate,K,LL,NY,NX)+OQAHXS
      SorbedOM_vr(idom_acetate,K,LL,NY,NX)=SorbedOM_vr(idom_acetate,K,LL,NY,NX)+OHAXS
      DOM(idom_don,K,LL,NY,NX)=DOM(idom_don,K,LL,NY,NX)+OQNXS
      DOM_Macp(idom_don,K,LL,NY,NX)=DOM_Macp(idom_don,K,LL,NY,NX)+OQNHXS
      SorbedOM_vr(ielmn,K,LL,NY,NX)=SorbedOM_vr(ielmn,K,LL,NY,NX)+OHNXS
      DOM(idom_dop,K,LL,NY,NX)=DOM(idom_dop,K,LL,NY,NX)+OQPXS
      DOM_Macp(idom_dop,K,LL,NY,NX)=DOM_Macp(idom_dop,K,LL,NY,NX)+OQPHXS
      SorbedOM_vr(ielmp,K,LL,NY,NX)=SorbedOM_vr(ielmp,K,LL,NY,NX)+OHPXS
      D7931: DO M=1,jsken
        IF(FOSCXS.GT.0.0_r8)THEN
          OSCXS=FOSCXS*AZMAX1(SolidOM_vr(ielmc,M,K,L,NY,NX))
          OSAXS=FOSCXS*AZMAX1(SolidOMAct_vr(M,K,L,NY,NX))
          OSNXS=FOSCXS*AZMAX1(SolidOM_vr(ielmn,M,K,L,NY,NX))
          OSPXS=FOSCXS*AZMAX1(SolidOM_vr(ielmp,M,K,L,NY,NX))
        ELSE
          OSCXS=FOSCXS*AZMAX1(SolidOM_vr(ielmc,M,K,LL,NY,NX))
          OSNXS=FOSCXS*AZMAX1(SolidOM_vr(ielmn,M,K,LL,NY,NX))
          OSPXS=FOSCXS*AZMAX1(SolidOM_vr(ielmp,M,K,LL,NY,NX))
          OSAXS=FOSCXS*AZMAX1(SolidOMAct_vr(M,K,LL,NY,NX))          
        ENDIF
        SolidOM_vr(ielmc,M,K,L,NY,NX)=SolidOM_vr(ielmc,M,K,L,NY,NX)-OSCXS
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)-OSAXS
        SolidOM_vr(ielmn,M,K,L,NY,NX)=SolidOM_vr(ielmn,M,K,L,NY,NX)-OSNXS
        SolidOM_vr(ielmp,M,K,L,NY,NX)=SolidOM_vr(ielmp,M,K,L,NY,NX)-OSPXS
        SolidOM_vr(ielmc,M,K,LL,NY,NX)=SolidOM_vr(ielmc,M,K,LL,NY,NX)+OSCXS
        SolidOMAct_vr(M,K,LL,NY,NX)=SolidOMAct_vr(M,K,LL,NY,NX)+OSAXS
        SolidOM_vr(ielmn,M,K,LL,NY,NX)=SolidOM_vr(ielmn,M,K,LL,NY,NX)+OSNXS
        SolidOM_vr(ielmp,M,K,LL,NY,NX)=SolidOM_vr(ielmp,M,K,LL,NY,NX)+OSPXS
      ENDDO D7931
    ENDDO D7901
  ENDIF
  end subroutine ApplyVerticalMix

end module nitrosMod
