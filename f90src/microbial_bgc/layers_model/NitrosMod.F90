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

  character(len=*), parameter :: mod_filename = __FILE__

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
!     VSoilPoreMicP=soil layer volume
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
          IF(VSoilPoreMicP(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
            LL=LN
            exit
          ENDIF
        ENDDO D1100
        ORGRL=AZMAX1(ORGR(L,NY,NX))
        ORGRLL=AZMAX1(ORGR(LL,NY,NX))
        OSCXD=(ORGRL*VOLT(LL,NY,NX)-ORGRLL*VOLT(L,NY,NX))/(VOLT(L,NY,NX)+VOLT(LL,NY,NX))
        IF(OSCXD.GT.0.0_r8.AND.ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXD=OSCXD/ORGR(L,NY,NX)
        ELSEIF(OSCXD.LT.0.0_r8.AND.ORGR(LL,NY,NX).GT.ZEROS(NY,NX))THEN
          FOSCXD=OSCXD/ORGR(LL,NY,NX)
        ELSE
          FOSCXD=0.0_r8
        ENDIF
        IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FOSCXS=FOSCZL*FOSCXD*TOQCK(L,NY,NX)/VOLT(L,NY,NX)
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
  integer :: K,M,N,NGL
!     begin_execution
  IF(FOSCXS.GT.ZERO)THEN
    D7971: DO K=1,micpar%n_litrsfk
      if(.not.micpar%is_finelitter(K))cycle
      D7961: DO N=1,NFGs
        DO NGL=JGnio(N),JGnfo(N)
          D7962: DO M=1,micpar%nlbiomcp
            IF(FOSCXS.GT.0.0)THEN
              OMCXS=FOSCXS*AZMAX1(OMC(M,NGL,K,L,NY,NX))
              OMNXS=FOSCXS*AZMAX1(OMN(M,NGL,K,L,NY,NX))
              OMPXS=FOSCXS*AZMAX1(OMP(M,NGL,K,L,NY,NX))
            ELSE
              OMCXS=FOSCXS*AZMAX1(OMC(M,NGL,K,LL,NY,NX))
              OMNXS=FOSCXS*AZMAX1(OMN(M,NGL,K,LL,NY,NX))
              OMPXS=FOSCXS*AZMAX1(OMP(M,NGL,K,LL,NY,NX))
            ENDIF
            OMC(M,NGL,K,L,NY,NX)=OMC(M,NGL,K,L,NY,NX)-OMCXS
            OMN(M,NGL,K,L,NY,NX)=OMN(M,NGL,K,L,NY,NX)-OMNXS
            OMP(M,NGL,K,L,NY,NX)=OMP(M,NGL,K,L,NY,NX)-OMPXS
            OMC(M,NGL,K,LL,NY,NX)=OMC(M,NGL,K,LL,NY,NX)+OMCXS
            OMN(M,NGL,K,LL,NY,NX)=OMN(M,NGL,K,LL,NY,NX)+OMNXS
            OMP(M,NGL,K,LL,NY,NX)=OMP(M,NGL,K,LL,NY,NX)+OMPXS
          ENDDO D7962
        ENDDO
      ENDDO D7961
    ENDDO D7971

    D7901: DO K=1,micpar%n_litrsfk
      if(.not.micpar%is_finelitter(K))cycle
      D7941: DO M=1,micpar%ndbiomcp
        IF(FOSCXS.GT.0.0_r8)THEN
          ORCXS=FOSCXS*AZMAX1(ORC(M,K,L,NY,NX))
          ORNXS=FOSCXS*AZMAX1(ORN(M,K,L,NY,NX))
          ORPXS=FOSCXS*AZMAX1(ORP(M,K,L,NY,NX))
        ELSE
          ORCXS=FOSCXS*AZMAX1(ORC(M,K,LL,NY,NX))
          ORNXS=FOSCXS*AZMAX1(ORN(M,K,LL,NY,NX))
          ORPXS=FOSCXS*AZMAX1(ORP(M,K,LL,NY,NX))
        ENDIF
        ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-ORCXS
        ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ORNXS
        ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-ORPXS
        ORC(M,K,LL,NY,NX)=ORC(M,K,LL,NY,NX)+ORCXS
        ORN(M,K,LL,NY,NX)=ORN(M,K,LL,NY,NX)+ORNXS
        ORP(M,K,LL,NY,NX)=ORP(M,K,LL,NY,NX)+ORPXS
      ENDDO D7941
      IF(FOSCXS.GT.0.0_r8)THEN
        OQCXS=FOSCXS*AZMAX1(OQC(K,L,NY,NX))
        OQCHXS=FOSCXS*AZMAX1(OQCH(K,L,NY,NX))
        OHCXS=FOSCXS*AZMAX1(OHC(K,L,NY,NX))
        OQAXS=FOSCXS*AZMAX1(OQA(K,L,NY,NX))
        OQAHXS=FOSCXS*AZMAX1(OQAH(K,L,NY,NX))
        OHAXS=FOSCXS*AZMAX1(OHA(K,L,NY,NX))
        OQNXS=FOSCXS*AZMAX1(OQN(K,L,NY,NX))
        OQNHXS=FOSCXS*AZMAX1(OQNH(K,L,NY,NX))
        OHNXS=FOSCXS*AZMAX1(OHN(K,L,NY,NX))
        OQPXS=FOSCXS*AZMAX1(OQP(K,L,NY,NX))
        OQPHXS=FOSCXS*AZMAX1(OQPH(K,L,NY,NX))
        OHPXS=FOSCXS*AZMAX1(OHP(K,L,NY,NX))
      ELSE
        OQCXS=FOSCXS*AZMAX1(OQC(K,LL,NY,NX))
        OQCHXS=FOSCXS*AZMAX1(OQCH(K,LL,NY,NX))
        OHCXS=FOSCXS*AZMAX1(OHC(K,LL,NY,NX))
        OQAXS=FOSCXS*AZMAX1(OQA(K,LL,NY,NX))
        OQAHXS=FOSCXS*AZMAX1(OQAH(K,LL,NY,NX))
        OHAXS=FOSCXS*AZMAX1(OHA(K,LL,NY,NX))
        OQNXS=FOSCXS*AZMAX1(OQN(K,LL,NY,NX))
        OQNHXS=FOSCXS*AZMAX1(OQNH(K,LL,NY,NX))
        OHNXS=FOSCXS*AZMAX1(OHN(K,LL,NY,NX))
        OQPXS=FOSCXS*AZMAX1(OQP(K,LL,NY,NX))
        OQPHXS=FOSCXS*AZMAX1(OQPH(K,LL,NY,NX))
        OHPXS=FOSCXS*AZMAX1(OHP(K,LL,NY,NX))
      ENDIF
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OQCXS
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OQCHXS
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OHCXS
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OQAXS
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OQAHXS
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OHAXS
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-OQNXS
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-OQNHXS
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-OHNXS
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OQPXS
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OQPHXS
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OHPXS
      OQC(K,LL,NY,NX)=OQC(K,LL,NY,NX)+OQCXS
      OQCH(K,LL,NY,NX)=OQCH(K,LL,NY,NX)+OQCHXS
      OHC(K,LL,NY,NX)=OHC(K,LL,NY,NX)+OHCXS
      OQA(K,LL,NY,NX)=OQA(K,LL,NY,NX)+OQAXS
      OQAH(K,LL,NY,NX)=OQAH(K,LL,NY,NX)+OQAHXS
      OHA(K,LL,NY,NX)=OHA(K,LL,NY,NX)+OHAXS
      OQN(K,LL,NY,NX)=OQN(K,LL,NY,NX)+OQNXS
      OQNH(K,LL,NY,NX)=OQNH(K,LL,NY,NX)+OQNHXS
      OHN(K,LL,NY,NX)=OHN(K,LL,NY,NX)+OHNXS
      OQP(K,LL,NY,NX)=OQP(K,LL,NY,NX)+OQPXS
      OQPH(K,LL,NY,NX)=OQPH(K,LL,NY,NX)+OQPHXS
      OHP(K,LL,NY,NX)=OHP(K,LL,NY,NX)+OHPXS
      D7931: DO M=1,jsken
        IF(FOSCXS.GT.0.0_r8)THEN
          OSCXS=FOSCXS*AZMAX1(OSC(M,K,L,NY,NX))
          OSAXS=FOSCXS*AZMAX1(OSA(M,K,L,NY,NX))
          OSNXS=FOSCXS*AZMAX1(OSN(M,K,L,NY,NX))
          OSPXS=FOSCXS*AZMAX1(OSP(M,K,L,NY,NX))
        ELSE
          OSCXS=FOSCXS*AZMAX1(OSC(M,K,LL,NY,NX))
          OSAXS=FOSCXS*AZMAX1(OSA(M,K,LL,NY,NX))
          OSNXS=FOSCXS*AZMAX1(OSN(M,K,LL,NY,NX))
          OSPXS=FOSCXS*AZMAX1(OSP(M,K,LL,NY,NX))
        ENDIF
        OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OSCXS
        OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OSAXS
        OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-OSNXS
        OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OSPXS
        OSC(M,K,LL,NY,NX)=OSC(M,K,LL,NY,NX)+OSCXS
        OSA(M,K,LL,NY,NX)=OSA(M,K,LL,NY,NX)+OSAXS
        OSN(M,K,LL,NY,NX)=OSN(M,K,LL,NY,NX)+OSNXS
        OSP(M,K,LL,NY,NX)=OSP(M,K,LL,NY,NX)+OSPXS
      ENDDO D7931
    ENDDO D7901
  ENDIF
  end subroutine ApplyVerticalMix

end module nitrosMod
