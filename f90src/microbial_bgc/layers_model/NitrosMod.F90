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
  use EcoSIMConfig , only : ndbiomcp => NumDeadMicrbCompts  
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

!
  public :: initNitro
  public :: DownwardMixOM
  public :: sumORGMLayL
  public :: sumLitrOMLayL
  public :: sumSurfOMCK
  public :: sumMicBiomLayL
  public :: sumHumOMLayL
  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro

  implicit none

  call initNitro1Layer

  end subroutine initNitro

!------------------------------------------------------------------------------------------

  subroutine DownwardMixOM(I,J,L,NY,NX)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: FracLitrMix,FOSCXD
  real(r8) :: ORGRL,ORGRLL
  real(r8) :: OSCXD
  integer :: LL,LN
!     begin_execution

  IF(FOSCZ0.GT.ZERO)THEN
!     OMLitrC_vr=total litter C
!     FOSCZ0=rate constant for mixing surface litter
!     FracLitrMix=mixing fraction for surface litter
!     TMicHeterAct_vr=total active biomass respiration activity
!     VLSoilPoreMicP_vr=soil layer volume
!     OSCXD=mixing required for equilibrating litter concentration
!     FOSCXD=mixing fraction for equilibrating subsurface litter
!     FracLitrMix=mixing fraction for subsurface litter
!
    IF(L.LT.NL(NY,NX))THEN
!     get mixing rate
      IF(L.EQ.0)THEN
        LL=NU(NY,NX)
        IF(OMLitrC_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FracLitrMix=AMIN1(1.0_r8,FOSCZ0/OMLitrC_vr(L,NY,NX)*TMicHeterAct_vr(L,NY,NX))
        ELSE
          FracLitrMix=0.0_r8
        ENDIF
      ELSE
        D1100: DO LN=L+1,NL(NY,NX)
          IF(VLSoilPoreMicP_vr(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
            LL=LN
            exit
          ENDIF
        ENDDO D1100
        ORGRL=AZMAX1(OMLitrC_vr(L,NY,NX))
        ORGRLL=AZMAX1(OMLitrC_vr(LL,NY,NX))
        OSCXD=(ORGRL*VGeomLayer(LL,NY,NX)-ORGRLL*VGeomLayer(L,NY,NX))/(VGeomLayer(L,NY,NX)+VGeomLayer(LL,NY,NX))
        IF(OSCXD.GT.0.0_r8.AND.OMLitrC_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXD=OSCXD/OMLitrC_vr(L,NY,NX)
        ELSEIF(OSCXD.LT.0.0_r8.AND.OMLitrC_vr(LL,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXD=OSCXD/OMLitrC_vr(LL,NY,NX)
        ELSE
          FOSCXD=0.0_r8
        ENDIF
        IF(VGeomLayer(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FracLitrMix=FOSCZL*FOSCXD*TMicHeterAct_vr(L,NY,NX)/VGeomLayer(L,NY,NX)
        ELSE
          FracLitrMix=0.0_r8
        ENDIF
      ENDIF

!     apply mixing
      call ApplyVerticalMix(FracLitrMix,L,LL,NY,NX)
    ENDIF

  ENDIF
  end subroutine DownwardMixOM
!------------------------------------------------------------------------------------------

  subroutine ApplyVerticalMix(FracLitrMix,L,LL,NY,NX)

  implicit none
  real(r8), intent(in) :: FracLitrMix
  integer, intent(in) :: L,LL,NY,NX

  real(r8) :: OMEXS,ORMXS
  real(r8) :: OQMXS,OQMHXS,OHMXS
  real(r8) :: OSMXS,OSAXS
  integer :: K,M,N,NGL,MID,NE,L1
  
!     begin_execution
  IF(FracLitrMix.GT.ZERO)THEN
    !mix microbial biomass
    D7971: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      D7961: DO N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          D7962: DO M=1,micpar%nlbiomcp
            MID=micpar%get_micb_id(M,NGL)            
            IF(FracLitrMix.GT.0.0_r8)THEN
              L1=L
            ELSE
              L1=LL
            ENDIF
            DO NE=1,NumPlantChemElms            
              OMEXS=FracLitrMix*AZMAX1(mBOMHeter_vr(NE,MID,K,L1,NY,NX))
              mBOMHeter_vr(NE,MID,K,L,NY,NX)=mBOMHeter_vr(NE,MID,K,L,NY,NX)-OMEXS
              mBOMHeter_vr(NE,MID,K,LL,NY,NX)=mBOMHeter_vr(NE,MID,K,LL,NY,NX)+OMEXS
            ENDDO
          ENDDO D7962
        ENDDO
      ENDDO D7961
    ENDDO D7971

    D7901: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      !mix fine litter
      D7941: DO M=1,micpar%ndbiomcp
        IF(FracLitrMix.GT.0.0_r8)THEN
          L1=L
        ELSE
          L1=LL
        ENDIF

        DO NE=1,NumPlantChemElms
          ORMXS=FracLitrMix*AZMAX1(OMBioResdu_vr(NE,M,K,L1,NY,NX))        
          OMBioResdu_vr(NE,M,K,L,NY,NX)=OMBioResdu_vr(NE,M,K,L,NY,NX)-ORMXS
          OMBioResdu_vr(NE,M,K,LL,NY,NX)=OMBioResdu_vr(NE,M,K,LL,NY,NX)+ORMXS
        ENDDO
      ENDDO D7941
      !mix dissolved organic matter
      IF(FracLitrMix.GT.0.0_r8)THEN
        L1=L
      ELSE
        L1=LL
      ENDIF      
      DO NE=idom_beg,idom_end
        OQMXS=FracLitrMix*AZMAX1(DOM_vr(NE,K,L1,NY,NX))
        OQMHXS=FracLitrMix*AZMAX1(DOM_MacP_vr(NE,K,L1,NY,NX))
        OHMXS=FracLitrMix*AZMAX1(SorbedOM_vr(NE,K,L1,NY,NX))

        DOM_vr(NE,K,L,NY,NX)=DOM_vr(NE,K,L,NY,NX)-OQMXS
        DOM_MacP_vr(NE,K,L,NY,NX)=DOM_MacP_vr(NE,K,L,NY,NX)-OQMHXS
        SorbedOM_vr(NE,K,L,NY,NX)=SorbedOM_vr(NE,K,L,NY,NX)-OHMXS

        DOM_vr(NE,K,LL,NY,NX)=DOM_vr(NE,K,LL,NY,NX)+OQMXS
        DOM_MacP_vr(NE,K,LL,NY,NX)=DOM_MacP_vr(NE,K,LL,NY,NX)+OQMHXS
        SorbedOM_vr(NE,K,LL,NY,NX)=SorbedOM_vr(NE,K,LL,NY,NX)+OHMXS
      ENDDO

      !mix solid organic matter
      D7931: DO M=1,jsken
        IF(FracLitrMix.GT.0.0_r8)THEN
          L1=L
        ELSE
          L1=LL
        ENDIF
        DO NE=1,NumPlantChemElms
          OSMXS=FracLitrMix*AZMAX1(SolidOM_vr(NE,M,K,L1,NY,NX))
          SolidOM_vr(NE,M,K,L,NY,NX)=SolidOM_vr(NE,M,K,L,NY,NX)-OSMXS
          SolidOM_vr(NE,M,K,LL,NY,NX)=SolidOM_vr(NE,M,K,LL,NY,NX)+OSMXS
        ENDDO
        OSAXS=FracLitrMix*AZMAX1(SolidOMAct_vr(M,K,L1,NY,NX))
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)-OSAXS
        SolidOMAct_vr(M,K,LL,NY,NX)=SolidOMAct_vr(M,K,LL,NY,NX)+OSAXS
      ENDDO D7931
    ENDDO D7901
  ENDIF
  end subroutine ApplyVerticalMix

!------------------------------------------------------------------------------------------
  subroutine sumORGMLayL(L,NY,NX,ORGM,conly)
  !
  !sum up organic matter in layer L
  !including live microbial biomass, microbial residue, sorbed dom+acetate, and solid SOM
  
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  logical, optional, intent(in) :: conly
  integer :: K,N,NGL,M,MID,NE,jcplx1,nelms
  logical :: conly_loc

  if(present(conly))then
    conly_loc=conly
  else
    conly_loc=.false.
  endif

  ORGM=0._r8
  !sumup heterotrophic microbes
  if(L==0)then
    jcplx1=micpar%NumOfLitrCmplxs
  else
    jcplx1=jcplx
  endif
  if(conly_loc)then
    nelms=1
  else
    nelms=NumPlantChemElms
  endif  
  !add autotrophic microbes
  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,nelms
          ORGM(NE)=ORGM(NE)+mBOMAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  DO K=1,jcplx1
    !add heterotrophic microbes
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,nelms
            ORGM(NE)=ORGM(NE)+mBOMHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo

  !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,nelms
        ORGM(NE)=ORGM(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
      ENDDO    
    ENDDO

    !add dom
    DO NE=1,nelms
      ORGM(NE)=ORGM(NE)+DOM_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    ORGM(ielmc)=ORGM(ielmc)+DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,nelms
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO    

  end subroutine sumORGMLayL
!------------------------------------------------------------------------------------------

  subroutine sumLitrOMLayL(L,NY,NX,ORGM)
  !
  !sum up litter OM (DOM + litter OM) in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer :: K,N,NGL,M,MID,NE

  ORGM=0._r8

  !add autotrophic microbes
  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          ORGM(NE)=ORGM(NE)+mBOMAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo


  DO K=1,micpar%NumOfLitrCmplxs
    !add live heterotrophic microbes
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            ORGM(NE)=ORGM(NE)+mBOMHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo

   !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
      ENDDO    
    ENDDO

    !add dom
    DO NE=1,NumPlantChemElms
      ORGM(NE)=ORGM(NE)+DOM_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    ORGM(ielmc)=ORGM(ielmc)+DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO    

  end subroutine sumLitrOMLayL

!------------------------------------------------------------------------------------------

  subroutine sumHumOMLayL(L,NY,NX,ORGM)
  !
  !sum up litter OM in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer :: K,N,NGL,M,MID,NE

  ORGM=0._r8

  !add autotrophic microbes
  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          ORGM(NE)=ORGM(NE)+mBOMAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  DO  K=micpar%NumOfLitrCmplxs+1,jcplx
   !sumup heterotrophic microbes
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            ORGM(NE)=ORGM(NE)+mBOMHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo
  !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
      ENDDO    
    ENDDO

    !add dom
    DO NE=1,NumPlantChemElms
      ORGM(NE)=ORGM(NE)+DOM_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    !add acetate
    ORGM(ielmc)=ORGM(ielmc)+DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO    

  end subroutine sumHumOMLayL
  

!------------------------------------------------------------------------------------------

  subroutine sumMicBiomLayL(L,NY,NX,ORGM)
  !
  !sum up litter OM in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer :: K,N,NGL,M,MID,NE,jcplx1

  ORGM=0._r8

  if(L==0)then
    jcplx1=micpar%NumOfLitrCmplxs
  else
    jcplx1=jcplx
  endif

  !add autotrophic microbes
  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          ORGM(NE)=ORGM(NE)+mBOMAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  DO K=1,jcplx1
    !add heterotrophic microbes
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            ORGM(NE)=ORGM(NE)+mBOMHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo
  enddo  
  end subroutine sumMicBiomLayL
!------------------------------------------------------------------------------------------

  subroutine sumSurfOMCK(NY,NX,SOMHeterK,SOMAutor)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(out) :: SOMHeterK(1:micpar%NumOfLitrCmplxs)
  real(r8), intent(out) :: SOMAutor
  integer :: K,N,NGL,M,MID,NE,L

  SOMHeterK=0._r8
  SOMAutor=0._r8
  L=0

  DO  N=1,NumMicbFunGrupsPerCmplx
    do NGL=JGniA(n),JGnfA(n)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,ielmc
          SOMAutor=SOMAutor+mBOMAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  enddo

  DO K=1,micpar%NumOfLitrCmplxs
    DO  N=1,NumMicbFunGrupsPerCmplx
      do NGL=JGnio(n),JGnfo(n)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,ielmc
            SOMHeterK(K)=SOMHeterK(K)+mBOMHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo  
    enddo    

    DO  M=1,ndbiomcp
      DO NE=1,ielmc
        SOMHeterK(K)=SOMHeterK(K)+OMBioResdu_vr(NE,M,K,L,NY,NX)        
      ENDDO
    ENDDO  

    !add dom
    DO NE=1,ielmc
      SOMHeterK(K)=SOMHeterK(K)+DOM_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    !add acetate
    SOMHeterK(K)=SOMHeterK(K)+DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    DO M=1,jsken
      DO NE=1,ielmc
        SOMHeterK(K)=SOMHeterK(K)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO
    ENDDO  
  ENDDO
  
  end subroutine sumSurfOMCK
end module nitrosMod
