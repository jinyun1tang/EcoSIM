module SoilBGCNLayMod
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
  use MicrobeDiagTypes
  use SoilDisturbMod
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
  public :: SumMicbGroup
  public :: sumDOML
  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro

  implicit none

  call initNitro1Layer

  end subroutine initNitro

!------------------------------------------------------------------------------------------

  subroutine DownwardMixOM(I,J,L,NY,NX,FracLitrMix)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8), intent(out) :: FracLitrMix
  real(r8) :: FOSCXD
  real(r8) :: ORGRL,ORGRLL
  real(r8) :: OSCXD
  integer :: LL,LN,K
  real(r8) :: DOC_s,DOC_u,ActL,ActLL,ActD
!     begin_execution

  IF(FOSCZ0.GT.ZERO)THEN
!     OMLitrC_vr=total litter C
!     FOSCZ0=rate constant for mixing surface litter
!     FracLitrMix=mixing fraction for surface litter
!     TMicHeterActivity_vr=total active biomass respiration activity
!     VLSoilPoreMicP_vr=soil layer volume
!     OSCXD=mixing required for equilibrating litter concentration
!     FOSCXD=mixing fraction for equilibrating subsurface litter
!     FracLitrMix=mixing fraction for subsurface litter
!
    IF(L.LT.NL(NY,NX))THEN
!     get mixing rate
      IF(L.EQ.0)THEN
        LL=NU(NY,NX)
        DOC_s=0._r8;DOC_u=0._r8
        DO K=1,micpar%NumOfLitrCmplxs
          if(.not.micpar%is_finelitter(K))cycle
          DOC_s=DOC_s+DOM_vr(idom_doc,K,L,NY,NX)
          DOC_u=DOC_u+DOM_vr(idom_doc,K,LL,NY,NX)
        ENDDO
        ActL  = DOC_s
        ActLL = DOC_u
        ActD  = (ActL*OMLitrC_vr(LL,NY,NX)-ActLL*OMLitrC_vr(L,NY,NX))/(OMLitrC_vr(L,NY,NX)+OMLitrC_vr(LL,NY,NX))

        IF(OMLitrC_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FracLitrMix=FOSCZ0/OMLitrC_vr(L,NY,NX)*TMicHeterActivity_vr(L,NY,NX)
!          if(ActD>0._r8)then
!            FracLitrMix=AMIN1(1.0_r8,FracLitrMix*ActD/DOC_s)
!          elseif(actD<0._r8)then
!            FracLitrMix=AMIN1(1.0_r8,FracLitrMix*ActD/DOC_u)
!          endif
        ELSE
          FracLitrMix=0.0_r8
        ENDIF
!        write(113,*)I+J/24.,ActD,DOM_vr(idom_doc,micpar%k_fine_litr,L,NY,NX),FracLitrMix

      ELSE
        D1100: DO LN=L+1,NL(NY,NX)
          IF(VLSoilPoreMicP_vr(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
            LL=LN
            exit
          ENDIF
        ENDDO D1100
        ORGRL  = AZMAX1(OMLitrC_vr(L,NY,NX))
        ORGRLL = AZMAX1(OMLitrC_vr(LL,NY,NX))
        OSCXD  = (ORGRL*VGeomLayer_vr(LL,NY,NX)-ORGRLL*VGeomLayer_vr(L,NY,NX))/(VGeomLayer_vr(L,NY,NX)+VGeomLayer_vr(LL,NY,NX))

        IF(OSCXD.GT.0.0_r8 .AND. OMLitrC_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
         !mass gradient pointing from LL to L
          FOSCXD=OSCXD/OMLitrC_vr(L,NY,NX)
        ELSEIF(OSCXD.LT.0.0_r8 .AND. OMLitrC_vr(LL,NY,NX).GT.ZEROS(NY,NX))THEN
          !mass gradient pointing from layer L to LL
          FOSCXD=OSCXD/OMLitrC_vr(LL,NY,NX)
        ELSE
          FOSCXD=0.0_r8
        ENDIF

        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FracLitrMix=FOSCZL*FOSCXD*TMicHeterActivity_vr(L,NY,NX)/VGeomLayer_vr(L,NY,NX)          
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
  integer, intent(in) :: NY,NX        !horizontal location of the grid
  integer, intent(in) :: L            !source grid
  integer, intent(in) :: LL           !destination grid

  real(r8) :: OMEXS,ORMXS
  real(r8) :: OQMXS,OQMHXS,OHMXS
  real(r8) :: OSMXS,OSAXS
  integer :: K,M,N,NGL,MID,NE,L1
  
!     begin_execution
  IF(FracLitrMix.GT.0.0_r8)THEN
    L1=L
  ELSE
    L1=LL
  ENDIF
  !how about heat, water flux and capacity?
! only downward mixing is considered
! upward mixing is yet to be considered, even though FracLitrMix could be negative

  IF(FracLitrMix.GT.ZERO)THEN
    !mix microbial biomass
    D7971: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      D7961: DO N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          D7962: DO M=1,micpar%nlbiomcp
            MID=micpar%get_micb_id(M,NGL)            
            DO NE=1,NumPlantChemElms            
              OMEXS                             = FracLitrMix*AZMAX1(mBiomeHeter_vr(NE,MID,K,L1,NY,NX))
              mBiomeHeter_vr(NE,MID,K,L,NY,NX)  = mBiomeHeter_vr(NE,MID,K,L,NY,NX)-OMEXS
              mBiomeHeter_vr(NE,MID,K,LL,NY,NX) = mBiomeHeter_vr(NE,MID,K,LL,NY,NX)+OMEXS
            ENDDO
          ENDDO D7962
        ENDDO
      ENDDO D7961
    ENDDO D7971

    !mix microbial residual
    D7901: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      !mix fine litter
      D7941: DO M=1,micpar%ndbiomcp

        DO NE=1,NumPlantChemElms
          ORMXS                          = FracLitrMix*AZMAX1(OMBioResdu_vr(NE,M,K,L1,NY,NX))
          OMBioResdu_vr(NE,M,K,L,NY,NX)  = OMBioResdu_vr(NE,M,K,L,NY,NX)-ORMXS
          OMBioResdu_vr(NE,M,K,LL,NY,NX) = OMBioResdu_vr(NE,M,K,LL,NY,NX)+ORMXS
        ENDDO
      ENDDO D7941

      !mix dissolved organic matter

      DO NE=idom_beg,idom_end
        OQMXS  = FracLitrMix*AZMAX1(DOM_vr(NE,K,L1,NY,NX))
        OQMHXS = FracLitrMix*AZMAX1(DOM_MacP_vr(NE,K,L1,NY,NX))
        OHMXS  = FracLitrMix*AZMAX1(SorbedOM_vr(NE,K,L1,NY,NX))

        DOM_vr(NE,K,L,NY,NX)      = DOM_vr(NE,K,L,NY,NX)-OQMXS
        DOM_MacP_vr(NE,K,L,NY,NX) = DOM_MacP_vr(NE,K,L,NY,NX)-OQMHXS
        SorbedOM_vr(NE,K,L,NY,NX) = SorbedOM_vr(NE,K,L,NY,NX)-OHMXS

        DOM_vr(NE,K,LL,NY,NX)      = DOM_vr(NE,K,LL,NY,NX)+OQMXS
        DOM_MacP_vr(NE,K,LL,NY,NX) = DOM_MacP_vr(NE,K,LL,NY,NX)+OQMHXS
        SorbedOM_vr(NE,K,LL,NY,NX) = SorbedOM_vr(NE,K,LL,NY,NX)+OHMXS
      ENDDO

      !mix solid organic matter
      D7931: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          OSMXS                       = FracLitrMix*AZMAX1(SolidOM_vr(NE,M,K,L1,NY,NX))
          SolidOM_vr(NE,M,K,L,NY,NX)  = SolidOM_vr(NE,M,K,L,NY,NX)-OSMXS
          SolidOM_vr(NE,M,K,LL,NY,NX) = SolidOM_vr(NE,M,K,LL,NY,NX)+OSMXS
        ENDDO
        OSAXS                       = FracLitrMix*AZMAX1(SolidOMAct_vr(M,K,L1,NY,NX))
        SolidOMAct_vr(M,K,L,NY,NX)  = SolidOMAct_vr(M,K,L,NY,NX)-OSAXS
        SolidOMAct_vr(M,K,LL,NY,NX) = SolidOMAct_vr(M,K,LL,NY,NX)+OSAXS
      ENDDO D7931
    ENDDO D7901
  ENDIF
  end subroutine ApplyVerticalMix

!------------------------------------------------------------------------------------------
  subroutine sumORGMLayL(L,NY,NX,ORGM,conly,info)
  !
  !sum up organic matter in layer L
  !including live microbial biomass, microbial residue, sorbed dom+acetate, and solid SOM
  
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  logical, optional, intent(in) :: conly
  character(len=*), optional, intent(in) :: info
  integer :: K,N,NGL,M,MID,NE,jcplx1,nelms,idom
  logical :: conly_loc
  real(r8) :: orgm0(1:NumPlantChemElms)

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
          ORGM(NE)=ORGM(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  DK100: DO K=1,jcplx1
    !add heterotrophic microbes
    DC100: DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGnio(N),JGnfo(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,nelms
            ORGM(NE)=ORGM(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo DC100
  ENDDO DK100
!  write(*,*)'OMC',ORGM(ielmc)

  DK200: DO K=1,jcplx1
    !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,nelms
        ORGM(NE)=ORGM(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
        if(ORGM(NE)<0._r8)then
        print*,'orgm3',ORGM
        stop
        endif
      ENDDO    
    ENDDO

    !add dom
    DO idom=idom_beg,idom_end
      if(abs(DOM_vr(idom,K,L,NY,NX))<1.e-12_r8)DOM_vr(idom,K,L,NY,NX)=0._r8    
      if(abs(SorbedOM_vr(idom,K,L,NY,NX))<1.e-12_r8)SorbedOM_vr(idom,K,L,NY,NX)=0._r8
      if(abs(DOM_MacP_vr(idom,K,L,NY,NX))<1.e-12_r8)DOM_MacP_vr(idom,K,L,NY,NX)=0._r8
    ENDDO  

    DO NE=1,nelms      
      ORGM(NE)=ORGM(NE)+DOM_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
        if(ORGM(NE)<0._r8)then
        print*,'orgm2',NE,ORGM,DOM_vr(NE,K,L,NY,NX),DOM_MacP_vr(NE,K,L,NY,NX),SorbedOM_vr(NE,K,L,NY,NX)
        stop
        endif
    ENDDO

    ORGM(ielmc)=ORGM(ielmc)+DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    if(ORGM(ielmc)<0._r8)then
    print*,'orgmxxx',ORGM,DOM_vr(idom_acetate,K,L,NY,NX),DOM_MacP_vr(idom_acetate,K,L,NY,NX),SorbedOM_vr(idom_acetate,K,L,NY,NX)    
    print*,'L=',L
    if(present(info))print*,info
    stop
    endif

    !add solid organic matter, litter, manure etc
    DO  M=1,jsken
      DO NE=1,nelms
        orgm0(NE)=ORGM(NE)
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
        if(ORGM(NE)<0._r8)then
        print*,K,NE,nelms,trim(info)//' orgm1 in '//trim(mod_filename)
        print*,'lay ',L,ORGM,SolidOM_vr(:,M,K,L,NY,NX),orgm0
        stop
        endif
      ENDDO  
    ENDDO  
  ENDDO DK200   

  end subroutine sumORGMLayL
!------------------------------------------------------------------------------------------

  subroutine sumLitrOMLayL(L,NY,NX,ORGM,I,J)
  !
  !sum up litter OM (DOM + litter OM) in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer , optional, intent(in) :: I,J  
  integer :: K,N,NGL,M,MID,NE,idom
  real(r8) :: DOM_micp(idom_beg:idom_end)
  real(r8) :: DOM_macp(idom_beg:idom_end) 


  ORGM=0._r8
  DOM_micp=0._r8
  DOM_macp=0._r8
  !add autotrophic microbes
  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          ORGM(NE)=ORGM(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
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
            ORGM(NE)=ORGM(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
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
    DO idom=idom_beg,idom_end
      DOM_micp(idom) = DOM_micp(idom) + DOM_vr(NE,K,L,NY,NX)
      DOM_macp(idom) = DOM_macp(idom) + DOM_MacP_vr(NE,K,L,NY,NX)
    ENDDO

    !add dom
    DO NE=1,NumPlantChemElms
      ORGM(NE)=ORGM(NE)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    ORGM(ielmc)=ORGM(ielmc)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO    
  ORGM=ORGM+DOM_micp(1:NumPlantChemElms) + DOM_macp(1:NumPlantChemElms)

  ORGM(ielmc)=ORGM(ielmc)+DOM_micp(idom_acetate) + DOM_macp(idom_acetate)
  end subroutine sumLitrOMLayL

!------------------------------------------------------------------------------------------

  subroutine sumHumOMLayL(L,NY,NX,ORGM)
  !
  !sum up litter OM in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer :: K,N,NGL,M,MID,NE,idom
  real(r8) :: DOM_micp(idom_beg:idom_end)
  real(r8) :: DOM_macp(idom_beg:idom_end) 

  ORGM=0._r8
  DOM_micp=0._r8
  DOM_macp=0._r8

  !add autotrophic microbes
  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          ORGM(NE)=ORGM(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
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
            ORGM(NE)=ORGM(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
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
    DO idom=idom_beg,idom_end
      DOM_micp(idom) = DOM_micp(idom) + DOM_vr(NE,K,L,NY,NX)
      DOM_macp(idom) = DOM_macp(idom) + DOM_MacP_vr(NE,K,L,NY,NX)
    ENDDO
    DO NE=1,NumPlantChemElms
      ORGM(NE)=ORGM(NE)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    !add acetate
    ORGM(ielmc)=ORGM(ielmc)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO    

  ORGM=ORGM+DOM_micp(1:NumPlantChemElms) + DOM_macp(1:NumPlantChemElms)

  ORGM(ielmc)=ORGM(ielmc) +DOM_micp(idom_acetate) + DOM_macp(idom_acetate)

  
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
          ORGM(NE)=ORGM(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
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
            ORGM(NE)=ORGM(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
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
        NE=ielmc
        SOMAutor=SOMAutor+mBiomeAutor_vr(NE,MID,L,NY,NX)        
      ENDDO
    ENDDO
  enddo

  DO K=1,micpar%NumOfLitrCmplxs
    DO  N=1,NumMicbFunGrupsPerCmplx
      do NGL=JGnio(n),JGnfo(n)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          NE=ielmc
          SOMHeterK(K)=SOMHeterK(K)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)          
        enddo
      enddo  
    enddo    

    DO  M=1,ndbiomcp
      NE=ielmc
      SOMHeterK(K)=SOMHeterK(K)+OMBioResdu_vr(NE,M,K,L,NY,NX)              
    ENDDO  

    !add dom
    NE=ielmc
    SOMHeterK(K)=SOMHeterK(K)+DOM_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
    
    !add acetate
    SOMHeterK(K)=SOMHeterK(K)+DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    DO M=1,jsken
      NE=ielmc
      SOMHeterK(K)=SOMHeterK(K)+SolidOM_vr(NE,M,K,L,NY,NX)      
    ENDDO  
  ENDDO
  
  end subroutine sumSurfOMCK
!------------------------------------------------------------------------------------------

  subroutine SumMicbGroup(L,NY,NX,igroup,MicbE,isauto)
  !
  !Description
  !
  implicit none
  integer,  intent(in) :: L,NY,NX
  integer,  intent(in) :: igroup
  real(r8), intent(out):: MicbE(1:NumPlantChemElms)
  logical, optional, intent(in) :: isauto
  logical :: isauto_loc
  integer :: K,NE,M,MID,NGL
  if(present(isauto))then
    isauto_loc=isauto
  else
    isauto_loc=.false.
  endif  

  micBE=0._r8
  if(isauto_loc)then
    if(igroup /= micpar%mid_AmmoniaOxidBacter        .and. &
       igroup /= micpar%mid_NitriteOxidBacter        .and. & 
       igroup /= micpar%mid_AerobicMethanotrofBacter .and. &
       igroup /= micpar%mid_H2GenoMethanogArchea) then
      call endrun('undefined autotroph group in '//trim(mod_filename),__LINE__)
    endif

    do NGL=JGniA(igroup),JGnfA(igroup)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          micBE(NE)=micBE(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO

  else
    if(igroup /= micpar%mid_Aerob_HeteroBacter  .and. &
       igroup /= micpar%mid_Facult_DenitBacter  .and. &
       igroup /= micpar%mid_Aerob_Fungi         .and. &
       igroup /= micpar%mid_fermentor           .and. &
       igroup /= micpar%mid_AcetoMethanogArchea .and. &
       igroup /= micpar%mid_aerob_N2Fixer       .and. &
       igroup /= micpar%mid_Anaerob_N2Fixer) then
      call endrun('undefined heterotroph group in '//trim(mod_filename),__LINE__)
    endif
    
    DO K=1,micpar%NumOfLitrCmplxs    
      do NGL=JGnio(igroup),JGnfo(igroup)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            micBE(NE)=micBE(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo  
    enddo    

  endif
  end subroutine SumMicbGroup

!------------------------------------------------------------------------------------------
  subroutine sumDOML(L,NY,NX,DOM)

  implicit none
  integer,  intent(in) :: L,NY,NX
  real(r8), intent(out):: DOM(idom_beg:idom_end)

  integer :: idom, K
  real(r8) :: DOM_micp(idom_beg:idom_end)
  real(r8) :: DOM_macp(idom_beg:idom_end) 
  
  DOM      = 0._r8
  DOM_micp = 0._r8
  DOM_macp = 0._r8
  DO K = 1, micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      DOM_micp(idom) = DOM_micp(idom)+DOM_vr(idom,K,L,NY,NX)
      DOM_macp(idom) = DOM_macp(idom)+DOM_MacP_vr(idom,K,L,NY,NX)
    ENDDO
  ENDDO
  DOM=DOM_micp+DOM_macp

  end subroutine sumDOML
end module SoilBGCNLayMod
