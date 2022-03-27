module nitro1LayerMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils  , only : endrun
  use minimathmod, only : safe_adb
  use MicrobialDataType
  use NitroPars
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use NitroDiagTypes
  use GridConsts
  use GridDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use EcosimConst
  use SurfLitterDataType
  use RootDataType
  use EcosysBGCFluxType
  use SoilPropertyDataType
  use IrrigationDataType
  use PlantDataRateType
  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = __FILE__

!
  public :: initNitro1Layer, SoilBGCOneLayer

  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro1Layer

  implicit none

  call initNitroPars

  end subroutine initNitro1Layer

!------------------------------------------------------------------------------------------

  subroutine SoilBGCOneLayer(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX

  integer :: LL,K,KL,NGL
  integer :: M,N
  type(NitroMicDiagType) :: nmicdiag
  type(NitroAQMFluxDiagType) :: naqfdiag
  type(NitroMicStateType) :: nmics
  type(NitroMicFluxType) :: nmicf
  type(NitroOMcplxFluxType) :: ncplxf
  type(NitroOMcplxStateType) :: ncplxs

! begin_execution
  call nmicf%Init(JG,jcplx)
  call nmics%Init(JG,jcplx)
  call ncplxf%Init()
  call ncplxs%Init()
  call naqfdiag%ZeroOut()


  IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
!
      !write(*,*)'StageBGCEnvironCondition'
    call StageBGCEnvironCondition(L,NY,NX,KL,nmicdiag,nmics,ncplxs)
!
      !write(*,*)'MicrobialCatabolism'
    call MicrobialCatabolism(L,NY,NX,nmicdiag,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
!
      !write(*,*)'ChemoDenitrification'
    call ChemoDenitrification(L,NY,NX,nmicdiag,naqfdiag)
!
!     DECOMPOSITION
!
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
    DO 1870 K=0,KL
      DO 1875 N=1,7
        DO NGL=1,JG
          ncplxf%ROQCK(K)=ncplxf%ROQCK(K)+nmicf%ROQCD(NGL,N,K)
        enddo
1875  CONTINUE
1870  CONTINUE
        !
        !write(*,*)'PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C'
    call OMTransferForPriming(KL,L,NY,NX,nmicf,nmics,ncplxf,ncplxs)
        !
        !     TRANSFER ALL PRIMING AMONG ALL K
        !
        !     TOQCK=total respiration of DOC+DOA in soil layer
        !     ROQCK=total respiration of DOC+DOA in substrate complex
        !     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
        !     OMC,OMN,OMP=microbial C,N,P
        !
    DO 1790 K=0,KL
!
          !write(*,*)'DECOMPOSITION OF ORGANIC SUBSTRATES'
!
      call SolidOMDecomposition(K,L,NY,NX,nmicdiag,ncplxf,ncplxs)
!
          !write(*,*)'DOC ADSORPTION - DESORPTION'
!
      call DOMSorption(K,L,NY,NX,nmicf,ncplxf,ncplxs)

1790  CONTINUE
        !write(*,*)'RedistDecompositionProduct'
      call RedistDecompositionProduct(KL,L,NY,NX,nmicdiag,nmicf,ncplxf,ncplxs)
!
        !write(*,*)'MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION'

      call MicrobialAnabolicUpdate(L,NY,NX,nmicf)
!
        !write(*,*)'MICROBIAL COLONIZATION OF NEW LITTER'
!
      call MicrobialLitterColonization(KL,L,NY,NX,ncplxf,ncplxs)
!
!     AGGREGATE ALL TRANSFORMATIONS CALCULATED ABOVE FOR EACH N,K
!
      call AggregateTransformations(L,NY,NX,nmicdiag,naqfdiag,nmicf,ncplxf)
  ELSE
    RCO2O(L,NY,NX)=0.0_r8
    RCH4O(L,NY,NX)=0.0_r8
    RH2GO(L,NY,NX)=0.0_r8
    RUPOXO(L,NY,NX)=0.0_r8
    RN2G(L,NY,NX)=0.0_r8
    RN2O(L,NY,NX)=0.0_r8
    XNH4S(L,NY,NX)=0.0_r8
    XNO3S(L,NY,NX)=0.0_r8
    XNO2S(L,NY,NX)=0.0_r8
    XH2PS(L,NY,NX)=0.0_r8
    XH1PS(L,NY,NX)=0.0_r8
    XNH4B(L,NY,NX)=0.0_r8
    XNO3B(L,NY,NX)=0.0_r8
    XNO2B(L,NY,NX)=0.0_r8
    XH2BS(L,NY,NX)=0.0_r8
    XH1BS(L,NY,NX)=0.0_r8
    XN2GS(L,NY,NX)=0.0_r8
  ENDIF
!

  call nmics%destroy()
  call nmicf%destroy()
  call ncplxf%destroy()
  call ncplxs%destroy()
  end subroutine SoilBGCOneLayer
!------------------------------------------------------------------------------------------

  subroutine StageBGCEnvironCondition(L,NY,NX,KL,nmicdiag,nmics,ncplxs)
  implicit none
  integer, intent(in) :: L,NY,NX
  integer, intent(out) :: KL
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroOMcplxStateType),intent(inout):: ncplxs
  integer  :: K
  integer  :: M,N,NGL
  real(r8) :: ACTV,ACTVM
  real(r8) :: ORGCL
  real(r8) :: RTK,STK,TKSO
  real(r8) :: TOSC,TOSA,TOHC
  real(r8) :: TSRH
!     begin_execution
  associate(                   &
    CNOMA  => nmics%CNOMA,     &
    CPOMA  => nmics%CPOMA,     &
    OMA    => nmics%OMA  ,     &
    OMC2   => nmics%OMC2 ,     &
    OMN2   => nmics%OMN2 ,     &
    FOM2   => nmics%FOM2 ,     &
    FCN    => nmics%FCN  ,     &
    FCP    => nmics%FCP  ,     &
    FCNP   => nmics%FCNP ,     &
    OSRH  => ncplxs%OSRH,      &
    TOMK  => ncplxs%TOMK,      &
    TONK  => ncplxs%TONK,      &
    TOPK  => ncplxs%TOPK,      &
    FOCA  => ncplxs%FOCA,      &
    FOAA  => ncplxs%FOAA,      &
    CNQ  => ncplxs%CNQ,        &
    CPQ  => ncplxs%CPQ,        &
    ORCT  => ncplxs%ORCT,      &
    OSCT  => ncplxs%OSCT,      &
    OSAT  => ncplxs%OSAT,      &
    TONX  => ncplxs%TONX,      &
    TOPX  => ncplxs%TOPX,      &
    TORC  =>  nmicdiag%TORC,   &
    TOMA  =>  nmicdiag%TOMA,   &
    TOMN  =>  nmicdiag%TOMN,   &
    XCO2  =>  nmicdiag%XCO2,   &
    TFNX  =>  nmicdiag%TFNX,   &
    TFNY  =>  nmicdiag%TFNY,   &
    THETR => nmicdiag%THETR,   &
    THETZ => nmicdiag%THETZ,   &
    VOLWZ => nmicdiag%VOLWZ,   &
    ZNH4T  =>  nmicdiag%ZNH4T, &
    ZNO3T  =>  nmicdiag%ZNO3T, &
    ZNO2T  =>  nmicdiag%ZNO2T, &
    H2P4T  =>  nmicdiag%H2P4T, &
    H1P4T  =>  nmicdiag%H1P4T  &
  )

! get KL, the number of mic-om complexes

!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH OFFSET FOR THERMAL ADAPTATION
  IF(L.EQ.0)THEN
    ! surface litter layer
    KL=2
    !write(*,*)'VOLR(NY,NX)=',VOLR(NY,NX)
    IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETR=VOLW(0,NY,NX)/VOLR(NY,NX)
      THETZ=AMAX1(0.0_r8,THETR-THETY(L,NY,NX))
      VOLWZ=THETZ*VOLR(NY,NX)
    ELSE
      VOLWZ=0.0_r8
    ENDIF
  ELSE
!     non-surface layer
    KL=4
    THETZ=AMAX1(0.0,(AMIN1(AMAX1(0.5*POROS(L,NY,NX),FC(L,NY,NX)) &
      ,THETW(L,NY,NX))-THETY(L,NY,NX)))
    VOLWZ=THETZ*VOLY(L,NY,NX)

  ENDIF


!     TKS=soil temperature
!     OFFSET=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TFNX,TFNY=temperature function for growth,maintenance respiration
!
  TKSO=TKS(L,NY,NX)+OFFSET(NY,NX)
  RTK=8.3143*TKSO
  STK=710.0*TKSO
  ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
  TFNX=EXP(25.229-62500/RTK)/ACTV
  ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
  TFNY=EXP(25.214-62500/RTK)/ACTVM
!
!     OXYI=inhibition of fermenters by O2
!     ORGCL=SOC used to calculate microbial concentration
!
!  OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS(L,NY,NX)+2.5)))
!  ORGCL=AMIN1(1.0E+05*BKVL(L,NY,NX),ORGC(L,NY,NX))
!
!     TOTAL MINERAL NH4, NO3 AND PO4
!
!     allocate NH4, NO3, HPO4, H2PO4 to non-band and band fractions
!
  if(L==0)then
    ZNH4T=AMAX1(0.0,ZNH4S(NU(NY,NX),NY,NX))+AMAX1(0.0,ZNH4B(NU(NY,NX),NY,NX))
    ZNO3T=AMAX1(0.0,ZNO3S(NU(NY,NX),NY,NX))+AMAX1(0.0,ZNO3B(NU(NY,NX),NY,NX))
    H1P4T=AMAX1(0.0,H1PO4(NU(NY,NX),NY,NX))+AMAX1(0.0,H1POB(NU(NY,NX),NY,NX))
    H2P4T=AMAX1(0.0,H2PO4(NU(NY,NX),NY,NX))+AMAX1(0.0,H2POB(NU(NY,NX),NY,NX))
    ZNO2T=AMAX1(0.0,ZNO2S(NU(NY,NX),NY,NX))+AMAX1(0.0,ZNO2B(NU(NY,NX),NY,NX))
  else
    ZNH4T=AMAX1(0.0,ZNH4S(L,NY,NX))+AMAX1(0.0,ZNH4B(L,NY,NX))
    ZNO3T=AMAX1(0.0,ZNO3S(L,NY,NX))+AMAX1(0.0,ZNO3B(L,NY,NX))
    H1P4T=AMAX1(0.0,H1PO4(L,NY,NX))+AMAX1(0.0,H1POB(L,NY,NX))
    ZNO2T=AMAX1(0.0,ZNO2S(L,NY,NX))+AMAX1(0.0,ZNO2B(L,NY,NX))
    H2P4T=AMAX1(0.0,H2PO4(L,NY,NX))+AMAX1(0.0,H2POB(L,NY,NX))
  endif
!
!     CCO2S=aqueous CO2 concentration
!
  XCO2=CCO2S(L,NY,NX)/(CCO2S(L,NY,NX)+CCKM)
!
!     TOTAL SUBSTRATE
!
!     TOSC=total SOC, TOSA=total colonized SOC
!     TORC=total microbial residue, TOHC=total adsorbed C
!     in each K:
!     OSCT=total SOC n each K, OSAT=total colonized SOC
!     ORCT=total microbial residue, OHCT=total adsorbed C
!
  TOSC=0.0_r8
  TOSA=0.0_r8
  TORC=0.0_r8
  TOHC=0.0_r8
!
!     TOTAL SOLID SUBSTRATE
!
  DO 870 K=0,KL
    OSCT(K)=0.0_r8
    OSAT(K)=0.0_r8
    DO 865 M=1,jsken
      OSCT(K)=OSCT(K)+OSC(M,K,L,NY,NX)
      OSAT(K)=OSAT(K)+OSA(M,K,L,NY,NX)
865 CONTINUE
    TOSC=TOSC+OSCT(K)
    TOSA=TOSA+OSAT(K)
870   CONTINUE
!
!     TOTAL BIORESIDUE
!
  DO 880 K=0,KL
    ORCT(K)=0.0_r8
    DO 875 M=1,2
      ORCT(K)=ORCT(K)+ORC(M,K,L,NY,NX)
!     IF(L.EQ.4.AND.K.EQ.2)THEN
!     WRITE(*,876)'ORCT',I,J,NX,NY,L,K,M,ORCT(K)
!    2,ORC(M,K,L,NY,NX)
!876   FORMAT(A8,7I4,60E12.4)
!     ENDIF
875 CONTINUE
    TORC=TORC+ORCT(K)
!
!     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
!
!     OSRH=total SOC
!
    TOHC=TOHC+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
880 CONTINUE
  DO 860 K=0,KL
    OSRH(K)=OSAT(K)+ORCT(K)+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
!     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
!     WRITE(*,861)'OSRH',I,J,NX,NY,L,K,OSRH(K),OSCT(K)
!    2,OSAT(K),ORCT(K),OHC(K,L,NY,NX),OHA(K,L,NY,NX)
!861   FORMAT(A8,6I4,20E12.4)
!     ENDIF
860 CONTINUE
  TSRH=TOSA+TORC+TOHC
!
!     C:N AND C:P RATIOS OF TOTAL BIOMASS
!     CNOMA,CPOMA=N,P contents of active biomass OMA
!     FCN,FCP=effects of N,P limitations on biomass activity
!
  TOMA=0.0_r8
  TOMN=0.0_r8
  DO 890 K=0,5
    IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 895 N=1,7
        IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
          DO NGL=1,JG
            IF(OMC(1,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
              CNOMA(NGL,N,K)=AMAX1(0.0,OMN(1,NGL,N,K,L,NY,NX)/OMC(1,NGL,N,K,L,NY,NX))
              CPOMA(NGL,N,K)=AMAX1(0.0,OMP(1,NGL,N,K,L,NY,NX)/OMC(1,NGL,N,K,L,NY,NX))
            ELSE
              CNOMA(NGL,N,K)=CNOMC(1,NGL,N,K)
              CPOMA(NGL,N,K)=CPOMC(1,NGL,N,K)
            ENDIF
            OMA(NGL,N,K)=AMAX1(0.0,OMC(1,NGL,N,K,L,NY,NX)/FL(1))
            FCN(NGL,N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CNOMA(NGL,N,K)/CNOMC(1,NGL,N,K))))
            FCP(NGL,N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CPOMA(NGL,N,K)/CPOMC(1,NGL,N,K))))
            FCNP(NGL,N,K)=AMIN1(FCN(NGL,N,K),FCP(NGL,N,K))
          ENDDO

!
  !       TOTAL BIOMASS
  !       OMC2=active biomass in recalcitrant fraction
!
          IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
            DO NGL=1,JG
              TOMA=TOMA+OMA(NGL,N,K)
            ENDDO
          ENDIF
          IF((K.LE.4.AND.N.EQ.2).OR.(K.EQ.5.AND.N.EQ.1))THEN
            DO NGL=1,JG
              TOMN=TOMN+OMA(NGL,N,K)
            ENDDO
          ENDIF
          DO NGL=1,JG
            OMC2(NGL,N,K)=AMAX1(0.0,AMIN1(OMA(NGL,N,K)*FL(2),OMC(2,NGL,N,K,L,NY,NX)))
            IF(OMC(2,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
              FOM2(NGL,N,K)=AMAX1(0.0,OMC2(NGL,N,K)/OMC(2,NGL,N,K,L,NY,NX))
              OMN2(NGL,N,K)=AMAX1(0.0,FOM2(NGL,N,K)*OMN(2,NGL,N,K,L,NY,NX))
            ELSE
              FOM2(NGL,N,K)=0.0_r8
              OMN2(NGL,N,K)=0.0_r8
            ENDIF
          ENDDO
        ENDIF
895   CONTINUE
    ENDIF
890 CONTINUE
  DO 690 K=0,KL
    TOMK(K)=0.0_r8
    TONK(K)=0.0_r8
    TOPK(K)=0.0_r8
    TONX(K)=0.0_r8
    TOPX(K)=0.0_r8
    DO 685 N=1,7
      DO NGL=1,JG
        TOMK(K)=TOMK(K)+OMA(NGL,N,K)
        TONK(K)=TONK(K)+OMA(NGL,N,K)*CNOMA(NGL,N,K)
        TOPK(K)=TOPK(K)+OMA(NGL,N,K)*CPOMA(NGL,N,K)
        TONX(K)=TONX(K)+OMA(NGL,N,K)*CNOMC(1,NGL,N,K)
        TOPX(K)=TOPX(K)+OMA(NGL,N,K)*CPOMC(1,NGL,N,K)
      ENDDO
685 CONTINUE
690 CONTINUE
!
!     FOSRH=fraction of total SOC in each substrate complex K
!
  DO 790 K=0,KL
    IF(TSRH.GT.ZEROS(NY,NX))THEN
      FOSRH(K,L,NY,NX)=OSRH(K)/TSRH
    ELSE
      FOSRH(K,L,NY,NX)=1.0
    ENDIF
    !
    !     DOC CONCENTRATIONS
    !
    !     COQC,COQA=aqueous DOC,acetate concentrations
    !     VOLWM=soil water content, FOSRH=fraction of total SOC
    !     occupied by each substrate complex K
    !
    IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(FOSRH(K,L,NY,NX).GT.ZERO)THEN
        COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX)/(VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)))
        COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX)/(VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)))
      ELSE
        COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX)/VOLWM(NPH,L,NY,NX))
        COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX)/VOLWM(NPH,L,NY,NX))
      ENDIF
    ELSE
      COQC(K,L,NY,NX)=0.0_r8
      COQA(K,L,NY,NX)=0.0_r8
    ENDIF
!
!     CNQ,CPQ=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
!
    IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNQ(K)=AMAX1(0.0,OQN(K,L,NY,NX)/OQC(K,L,NY,NX))
      CPQ(K)=AMAX1(0.0,OQP(K,L,NY,NX)/OQC(K,L,NY,NX))
    ELSE
      CNQ(K)=0.0_r8
      CPQ(K)=0.0_r8
    ENDIF
    IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX).AND.OQA(K,L,NY,NX) &
      .GT.ZEROS(NY,NX))THEN
      FOCA(K)=OQC(K,L,NY,NX)/(OQC(K,L,NY,NX)+OQA(K,L,NY,NX))
      FOAA(K)=1.0-FOCA(K)
    ELSEIF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOCA(K)=1.0
      FOAA(K)=0.0_r8
    ELSE
      FOCA(K)=0.0_r8
      FOAA(K)=1.0
    ENDIF
790 CONTINUE
!
  end associate
  end subroutine StageBGCEnvironCondition
!------------------------------------------------------------------------------------------

  subroutine MicrobialCatabolism(L,NY,NX,nmicdiag,naqfdiag, nmicf, nmics,ncplxf,ncplxs)
  !
  !  Description:
  !
  implicit none
  integer, intent(in) :: L,NY,NX
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout):: ncplxs
  integer :: K,M,N,NGL
  real(r8) :: TOMCNK(2)
  REAL(R8) :: OXKX
! begin_execution
  associate(                   &
    TFNR => nmics%TFNR,        &
    TFNG => nmics%TFNG,        &
    OMA  => nmics%OMA,         &
    TCGOQC  => ncplxf%TCGOQC,  &
    TCGOAC  => ncplxf%TCGOAC,  &
    TCGOMN  => ncplxf%TCGOMN,  &
    TCGOMP  => ncplxf%TCGOMP,  &
    XCO2  =>  nmicdiag%XCO2,   &
    TFNX  =>  nmicdiag%TFNX,   &
    TOMA  =>  nmicdiag%TOMA,   &
    TOMN  =>  nmicdiag%TOMN,   &
    RH2GZ  =>  nmicdiag%RH2GZ, &
    WFNG  =>  nmicdiag%WFNG,   &
    TFNY  =>  nmicdiag%TFNY,   &
    ZNH4T  =>  nmicdiag%ZNH4T, &
    ZNO3T  =>  nmicdiag%ZNO3T, &
    ZNO2T  =>  nmicdiag%ZNO2T, &
    H2P4T  =>  nmicdiag%H2P4T, &
    H1P4T  =>  nmicdiag%H1P4T  &
  )
  RH2GZ=0.0_r8
  DO 760 K=0,5
    IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      TCGOQC(K)=0.0_r8
      TCGOAC(K)=0.0_r8
      TCGOMN(K)=0.0_r8
      TCGOMP(K)=0.0_r8
      DO 750 N=1,7
        IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
          TOMCNK(:)=0.0
          DO NGL=1,JG
            DO M=1,2
              TOMCNK(M)=TOMCNK(M)+OMC(M,NGL,N,K,L,NY,NX)
            ENDDO
          ENDDO

          DO NGL=1,JG
            IF(K.LE.4)THEN
              IF(N.EQ.3)THEN
!               WFNG=water potential (PSISM) effect on microbial respiration
!               OXKX=Km for O2 uptake
!               OXKM=Km for heterotrophic O2 uptake set in starts.f
!               TFNG=combined temp and water stress effect on growth respiration
!               TFNR=temperature effect on maintenance respiration
                WFNG=EXP(0.1*PSISM(L,NY,NX))
              ELSE
                WFNG=EXP(0.2*PSISM(L,NY,NX))
              ENDIF
              OXKX=OXKM
            ELSE
              WFNG=EXP(0.2*PSISM(L,NY,NX))
              OXKX=OXKA
            ENDIF
            TFNG(NGL,N,K)=TFNX*WFNG
            TFNR(NGL,N,K)=TFNY
            IF(OMA(NGL,N,K).GT.0.0)THEN
              call ActiveMicrobes(NGL,N,K,L,NY,NX,XCO2,TFNX,WFNG,TOMCNK,&
                OXKX,TOMA,TOMN,RH2GZ,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,&
                naqfdiag,nmicf,nmics,ncplxf,ncplxs)
            ENDIF
          ENDDO
        ENDIF
750   CONTINUE
    ENDIF
760 CONTINUE
  end associate
  end subroutine MicrobialCatabolism
!------------------------------------------------------------------------------------------
  subroutine ActiveMicrobes(NGL,N,K,L,NY,NX,XCO2,TFNX,WFNG,TOMCNK,OXKX,TOMA,TOMN,RH2GZ,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,naqfdiag,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in):: OXKX,tomcnk(2),WFNG,TFNX,XCO2
  real(r8), intent(in) :: TOMA,TOMN
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  real(r8), intent(out):: RH2GZ
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  integer  :: M
  real(r8) :: COMC
  real(r8) :: ECHZ
  real(r8) :: ORGCL
  real(r8) :: FOXYX
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX
  real(r8) :: FOQC,FOQA
  real(r8) :: FGOCP,FGOAP
  real(r8) :: RGOCP
  real(r8) :: RGOMP
  real(r8) :: RVOXP
  real(r8) :: RVOXPA
  real(r8) :: RVOXPB
  real(r8) :: RGOMT
  real(r8) :: RXOMT
  real(r8) :: RMOMT
  real(r8) :: SPOMK(2)
  real(r8) :: RMOMK(2)
  REAL(R8) :: VOLWZ
! begin_execution
  associate(                  &
    FOMA => nmics%FOMA,       &
    FOMN => nmics%FOMN,       &
    FOMK => nmics%FOMK,       &
    OMA => nmics%OMA,         &
    RUPOX  => nmicf%RUPOX,    &
    RGOMO  => nmicf%RGOMO,    &
    ROXYM => nmicf%ROXYM ,    &
    ROXYP=> nmicf%ROXYP  ,    &
    ROXYO => nmicf%ROXYO ,    &
    RDNO3 => nmicf%RDNO3 ,    &
    RDNOB => nmicf%RDNOB ,    &
    RDNO2 => nmicf%RDNO2 ,    &
    RDN2B => nmicf%RDN2B ,    &
    RDN2O => nmicf%RDN2O ,    &
    RGOMD  => nmicf%RGOMD,    &
    RH2GX  => nmicf%RH2GX,    &
    RGOMY  => nmicf%RGOMY,    &
    RCO2X  => nmicf%RCO2X,    &
    RCH3X  => nmicf%RCH3X,    &
    RCH4X  => nmicf%RCH4X,    &
    TOMK  => ncplxs%TOMK      &
  )
! FOMA,FOMN=fraction of total active biomass C,N in each N and K

  IF(TOMA.GT.ZEROS(NY,NX))THEN
    FOMA(NGL,N,K)=OMA(NGL,N,K)/TOMA
  ELSE
    FOMA(NGL,N,K)=1.0
  ENDIF
  IF(TOMN.GT.ZEROS(NY,NX))THEN
    FOMN(NGL,N,K)=OMA(NGL,N,K)/TOMN
  ELSE
    FOMN(NGL,N,K)=1.0
  ENDIF
  IF(TOMK(K).GT.ZEROS(NY,NX))THEN
    FOMK(NGL,N,K)=OMA(NGL,N,K)/TOMK(K)
  ELSE
    FOMK(NGL,N,K)=1.0
  ENDIF
!
  !     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
  !
  !     COMC=microbial C concentration relative to substrate
  !     SPOMK=effect of microbial C concentration on microbial decay
  !     RMOMK=effect of microbial C concentration on maintenance respn
  !
  ORGCL=AMIN1(1.0E+05*BKVL(L,NY,NX),ORGC(L,NY,NX))
  IF(ORGCL.GT.ZEROS(NY,NX))THEN
    DO 765 M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
765 CONTINUE
  ELSE
    DO 770 M=1,2
      SPOMK(M)=1.0
      RMOMK(M)=1.0
770 CONTINUE
  ENDIF
!
! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
            !write(*,*)'SubstrateCompetitionFactors'
  call SubstrateCompetitionFactors(NGL,N,K,L,NY,NX,FOXYX,FNH4X,&
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA,&
    naqfdiag,nmicf,nmics)
!
! HETEROTROPHIC BIOMASS RESPIRATION
          !
  IF(K.LE.4)THEN
!
!   RESPIRATION BY HETEROTROPHIC AEROBES:
!   N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!   (6)N2 FIXERS
!
    IF(N.LE.3.OR.N.EQ.6)THEN
!     write(*,*)'AerobicHeterotrophCatabolism'
      call AerobicHeterotrophCatabolism(NGL,N,K,L,NY,NX,TFNX,WFNG,FOQC,FOQA,&
        ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,nmicf,nmics,ncplxs)
!     RESPIRATION BY HETEROTROPHIC ANAEROBES:
!     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
!
!     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND
!     ACETATE CONCENTRATION
!
!     GH2F=energy yield of acetotrophic methanogenesis per g C
!     GHAX=H2 effect on energy yield of fermentation
!     GOAX=acetate effect on energy yield of fermentation
!     ECHZ=growth respiration efficiency of fermentation
!
    ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
!     write(*,*)'AnaerobCatabolism'
      call AnaerobCatabolism(NGL,N,K,L,NY,NX,TFNX,WFNG,FOQC,ECHZ,FGOCP,FGOAP,RGOMP,&
        naqfdiag,nmicf,nmics)
!     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
!
!     GOMX=acetate effect on energy yield
!     ECHZ=growth respiration efficiency of aceto. methanogenesis
!
    ELSEIF(N.EQ.5)THEN
!     write(*,*)'AcetoMethanogenCatabolism'
      call AcetoMethanogenCatabolism(NGL,N,K,L,NY,NX,TFNX,WFNG,FOQA,ECHZ,&
        FGOCP,FGOAP,RGOMP,naqfdiag,nmicf,nmics)
    ENDIF
!
!     RESPIRATION RATES BY AUTOTROPHS 'RGOMP' FROM SPECIFIC
!     OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR. N=(1) NH4 OXIDIZERS (2) NO2
!     OXIDIZERS,(3) CH4 OXIDIZERS, (5) H2TROPHIC METHANOGENS
!
  ELSEIF(K.EQ.5)THEN
!
!   NH3 OXIDIZERS
!
    IF(N.EQ.1)THEN
      call NH3OxidizerCatabolism(NGL,N,K,L,NY,NX,XCO2,VOLWZ,TFNX,ECHZ,RGOMP,RVOXP,&
        RVOXPA,RVOXPB,naqfdiag,nmicf,nmics)
  !     NO2 OXIDIZERS
  !
    ELSEIF(N.EQ.2)THEN
!     write(*,*)'NO2OxidizerCatabolism'
      call NO2OxidizerCatabolism(NGL,N,K,L,NY,NX,XCO2,ECHZ,RGOMP,RVOXP,RVOXPA,RVOXPB,&
        naqfdiag,nmicf,nmics)
  !     H2TROPHIC METHANOGENS
  !
    ELSEIF(N.EQ.5)THEN

      call H2MethanogensCatabolism(NGL,N,K,L,NY,NX,ECHZ,RGOMP,XCO2,naqfdiag,nmicf,nmics)
  !     METHANOTROPHS
!
    ELSEIF(N.EQ.3)THEN
      call MethanotrophCatabolism(NGL,N,K,L,NY,NX,ECHZ,RGOMP,&
        RVOXP,RVOXPA,RVOXPB,naqfdiag,nmicf,nmics)
    ELSE
      RGOMP=0.0_r8
      ROXYM(NGL,N,K)=0.0_r8
      ROXYP(NGL,N,K)=0.0_r8
      ROXYS(NGL,N,K,L,NY,NX)=0.0_r8
    ENDIF
  ELSE
    RGOMP=0.0_r8
    ROXYM(NGL,N,K)=0.0_r8
    ROXYP(NGL,N,K)=0.0_r8
    ROXYS(NGL,N,K,L,NY,NX)=0.0_r8
  ENDIF
!
!  write(*,*)'O2 UPTAKE BY AEROBES'
!
! RUPOX, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
! RUPMX=O2-unlimited rate of O2 uptake
! FOXYX=fraction of O2 uptake by N,K relative to total
! XNPG=1/(NPH*NPT)
! ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
! OLSGL=aqueous O2 diffusivity
! OXYG,OXYS=gaseous, aqueous O2 amounts
! FLQRQ,FLQRI=surface water flux from precipitation, irrigation
! COXR,COXQ=O2 concentration in FLQRQ,FLQRI
!
  RUPOX(NGL,N,K)=0.0_r8
  IF(N.LE.3.OR.N.EQ.6)THEN
!  N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
!    (6)N2 FIXERS
!   write(*,*)'AerobsO2Uptake'
    call AerobsO2Uptake(NGL,N,K,L,NY,NX,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,&
      nmicf,nmics)
  ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
    RGOMO(NGL,N,K)=RGOMP
    RCO2X(NGL,N,K)=0.333*RGOMO(NGL,N,K)
    RCH3X(NGL,N,K)=0.667*RGOMO(NGL,N,K)
    RCH4X(NGL,N,K)=0.0_r8
    ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
    IF(K.LE.4)THEN
      RH2GX(NGL,N,K)=0.111*RGOMO(NGL,N,K)
    ELSE
      RH2GX(NGL,N,K)=0.0_r8
    ENDIF
  ELSEIF(N.EQ.5)THEN
    RGOMO(NGL,N,K)=RGOMP
    IF(K.LE.4)THEN
      RCO2X(NGL,N,K)=0.50_r8*RGOMO(NGL,N,K)
      RCH3X(NGL,N,K)=0.0_r8
      RCH4X(NGL,N,K)=0.50_r8*RGOMO(NGL,N,K)
      ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
      RH2GX(NGL,N,K)=0.0_r8
    ELSEIF(K.EQ.5)THEN
      RCO2X(NGL,N,K)=0.0_r8
      RCH3X(NGL,N,K)=0.0_r8
      RCH4X(NGL,N,K)=RGOMO(NGL,N,K)
      ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
      RH2GX(NGL,N,K)=0.0_r8
      RH2GZ=0.667_r8*RGOMO(NGL,N,K)
    ENDIF
  ENDIF
!
!  write(*,*)'HETEROTROPHIC DENITRIFICATION'
!
  IF(K.LE.4.AND.N.EQ.2.AND.ROXYM(NGL,N,K).GT.0.0 &
    .AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN
    call HeteroDenitrificCatabolism(NGL,N,K,L,NY,NX,FOQC,RGOCP,&
      VOLWZ,naqfdiag,nmicf,nmics,ncplxs)
!     AUTOTROPHIC DENITRIFICATION
!
  ELSEIF(K.EQ.5.AND.N.EQ.1.AND.ROXYM(NGL,N,K).GT.0.0 &
    .AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN

    call AutotrophDenitrificCatabolism(NGL,N,K,L,NY,NX,XCO2,VOLWZ,naqfdiag,nmicf,nmics)
  ELSE
    RDNO3(NGL,N,K)=0.0_r8
    RDNOB(NGL,N,K)=0.0_r8
    RDNO2(NGL,N,K)=0.0_r8
    RDN2B(NGL,N,K)=0.0_r8
    RDN2O(NGL,N,K)=0.0_r8
    RGOMY(NGL,N,K)=0.0_r8
    RGOMD(NGL,N,K)=0.0_r8
  ENDIF
!
!     BIOMASS DECOMPOSITION AND MINERALIZATION
!
  call BiomassMineralization(NGL,N,K,L,NY,NX,FNH4X, &
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,nmicf,nmics)
!
  call GatherMicrobialRespiration(NGL,N,K,L,NY,NX,RMOMK,RGOMT,RXOMT,RMOMT,nmicf,nmics)
!
  call GetMicrobialAnabolismFlux(NGL,N,K,L,NY,NX,ECHZ,FGOCP,&
    FGOAP,RGOMT,RXOMT,RMOMT,spomk,rmomk,nmicf,nmics,ncplxf,ncplxs)
  end associate
  end subroutine ActiveMicrobes
!------------------------------------------------------------------------------------------

  subroutine ChemoDenitrification(L,NY,NX,nmicdiag,naqfdiag)
  implicit none
  integer, intent(in) :: L,NY,NX
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroAQMFluxDiagType), INTENT(INOUT):: naqfdiag
  real(r8) :: chy1,CHNO2,CHNOB
  real(r8) :: FNO3S,FNO3B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: VMXC4S,VMXC4B
!     begin_execution
  associate(                   &
    TFNX  =>  nmicdiag%TFNX,   &
    RCNO2  =>  nmicdiag%RCNO2, &
    RCNOB  =>  nmicdiag%RCNOB, &
    RCN2O  =>  nmicdiag%RCN2O, &
    RCN2B  =>  nmicdiag%RCN2B, &
    RCNO3  =>  nmicdiag%RCNO3, &
    RCN3B  =>  nmicdiag%RCN3B, &
    RCOQN  =>  nmicdiag%RCOQN  &
  )
!
!     FNO2,FNB2=fraction of total NO2 demand in non-band,band
!     VMXC4S,VMXC4B=substrate-unlimited NO2 reduction in non-band,band
!     CHNO2,CHNOB=nitrous acid concentration in non-band,band
!     VOLWM=soil water content
!     FNO3S,FNO3B=fractions of NO2 in non-band,band
!     TFNX=temperature stress function
!     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
!     RCN2O,RCN2B=N2O production from nitrous acid reduction in non-band,band
!     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
!     RCOQN=DON production from nitrous acid reduction
!     RVMXC,RVMBC=demand for NO2 reduction in non-band,band
!     nitrous acid concn CHNO2
  CHY1=AMAX1(ZERO,10.0**(-(PH(L,NY,NX)-3.0)))
  CHNO2=CNO2S(L,NY,NX)*CHY1/0.5
  CHNOB=CNO2B(L,NY,NX)*CHY1/0.5

  IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO2=AMAX1(FMN,RVMXC(L,NY,NX)/RNO2Y(L,NY,NX))
  ELSE
    FNO2=FMN*VLNO3(L,NY,NX)
  ENDIF
  IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB2=AMAX1(FMN,RVMBC(L,NY,NX)/RN2BY(L,NY,NX))
  ELSE
    FNB2=FMN*VLNOB(L,NY,NX)
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
  FNO3S=VLNO3(L,NY,NX)
  FNO3B=VLNOB(L,NY,NX)
  VMXC4S=7.5E-02*CHNO2*VOLWM(NPH,L,NY,NX)*FNO3S*TFNX
  VMXC4B=7.5E-02*CHNOB*VOLWM(NPH,L,NY,NX)*FNO3B*TFNX
  RCNO2=AMAX1(0.0,AMIN1(ZNO2S(L,NY,NX)*FNO2,VMXC4S))
  RCNOB=AMAX1(0.0,AMIN1(ZNO2B(L,NY,NX)*FNB2,VMXC4B))
  RCN2O=0.10*RCNO2
  RCN2B=0.10*RCNOB
  RCNO3=0.80*RCNO2
  RCN3B=0.80*RCNOB
  RCOQN=0.10*(RCNO2+RCNOB)
  RVMXC(L,NY,NX)=VMXC4S
  RVMBC(L,NY,NX)=VMXC4B

  end associate
  end subroutine ChemoDenitrification
!------------------------------------------------------------------------------------------

  subroutine OMTransferForPriming(KL,L,NY,NX,nmicf,nmics,ncplxf,ncplxs)

  implicit none
  integer, intent(in) :: KL,L,NY,NX
  type(NitroMicFLuxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroOMcplxFluxType), intent(inout):: ncplxf
  type(NitroOMcplxStateType),intent(inout):: ncplxs
  integer  :: K,M,N,KK,NGL
  real(r8) :: OSRT
  real(r8) :: XFRK,XFRC
  real(r8) :: XFRN,XFRP,XFRA
  real(r8) :: XFMC,XFMN,XFMP
!     begin_execution
  associate(                 &
    TFNG  => nmics%TFNG,     &
    XOMCZ  => nmicf%XOMCZ,   &
    XOMNZ  => nmicf%XOMNZ,   &
    XOMPZ  => nmicf%XOMPZ,   &
    ROQCK  => ncplxf%ROQCK,  &
    XOQCK  => ncplxf%XOQCK,  &
    XOQCZ  => ncplxf%XOQCZ,  &
    XOQNZ  => ncplxf%XOQNZ,  &
    XOQPZ  => ncplxf%XOQPZ,  &
    XOQAZ  => ncplxf%XOQAZ,  &
    OSRH  => ncplxs%OSRH     &
  )
!
!     OSRH=total SOC in each K
!     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
!     between each K and KK, FPRIM=priming transfer rate constant
!     TFND=temperature effect on priming transfers
!     ROQCK,OQC,OQN,OQP=respiration,DOC,DON,DOP
!     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
!
  DO 795 K=0,KL
    IF(K.LE.KL-1)THEN
      DO 800 KK=K+1,KL
        OSRT=OSRH(K)+OSRH(KK)
        IF(OSRH(K).GT.ZEROS(NY,NX).AND.OSRH(KK).GT.ZEROS(NY,NX))THEN
          XFRK=FPRIM*TFND(L,NY,NX)*(ROQCK(K)*OSRH(KK)-ROQCK(KK)*OSRH(K))/OSRT
          XFRC=FPRIM*TFND(L,NY,NX)*(OQC(K,L,NY,NX)*OSRH(KK)-OQC(KK,L,NY,NX)*OSRH(K))/OSRT
          XFRN=FPRIM*TFND(L,NY,NX)*(OQN(K,L,NY,NX)*OSRH(KK)-OQN(KK,L,NY,NX)*OSRH(K))/OSRT
          XFRP=FPRIM*TFND(L,NY,NX)*(OQP(K,L,NY,NX)*OSRH(KK)-OQP(KK,L,NY,NX)*OSRH(K))/OSRT
          XFRA=FPRIM*TFND(L,NY,NX)*(OQA(K,L,NY,NX)*OSRH(KK)-OQA(KK,L,NY,NX)*OSRH(K))/OSRT
          IF(ROQCK(K)+XOQCK(K)-XFRK.GT.0.0 &
            .AND.ROQCK(KK)+XOQCK(KK)+XFRK.GT.0.0)THEN
            XOQCK(K)=XOQCK(K)-XFRK
            XOQCK(KK)=XOQCK(KK)+XFRK
!     IF(I.EQ.116)THEN
!     WRITE(*,4442)'XOQCK',I,J,NX,NY,L,K,KK,XFRC,ROQCK(K)
!    2,OSRH(K),ROQCK(KK),OSRH(KK),XOQCK(K),XOQCK(KK)
!    3,OQC(K,L,NY,NX),OQC(KK,L,NY,NX)
!4442  FORMAT(A8,7I4,12E12.4)
!     ENDIF
          ENDIF
          IF(OQC(K,L,NY,NX)+XOQCZ(K)-XFRC.GT.0.0 &
            .AND.OQC(KK,L,NY,NX)+XOQCZ(KK)+XFRC.GT.0.0)THEN
            XOQCZ(K)=XOQCZ(K)-XFRC
            XOQCZ(KK)=XOQCZ(KK)+XFRC
    !     IF(I.EQ.116)THEN
!     WRITE(*,4442)'XOQCZ',I,J,NX,NY,L,K,KK,XFRC,OQC(K,L,NY,NX)
!    2,OSRH(K),OQC(KK,L,NY,NX),OSRH(KK),XOQCZ(K),XOQCZ(KK)
!    3,OQC(K,L,NY,NX),OQC(KK,8,NY,NX)
!     ENDIF
          ENDIF
          IF(OQN(K,L,NY,NX)+XOQNZ(K)-XFRN.GT.0.0 &
            .AND.OQN(KK,L,NY,NX)+XOQNZ(KK)+XFRN.GT.0.0)THEN
            XOQNZ(K)=XOQNZ(K)-XFRN
            XOQNZ(KK)=XOQNZ(KK)+XFRN
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
!     WRITE(*,4442)'XOQNZ',I,J,NX,NY,L,K,KK,XFRN,OQN(K,L,NY,NX)
!    2,OSRH(K),OQN(KK,L,NY,NX),OSRH(KK),XOQNZ(K),XOQNZ(KK)
!     ENDIF
          ENDIF
          IF(OQP(K,L,NY,NX)+XOQPZ(K)-XFRP.GT.0.0 &
            .AND.OQP(KK,L,NY,NX)+XOQPZ(KK)+XFRP.GT.0.0)THEN
            XOQPZ(K)=XOQPZ(K)-XFRP
            XOQPZ(KK)=XOQPZ(KK)+XFRP
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
!     WRITE(*,4442)'XOQPZ',I,J,NX,NY,L,K,KK,XFRP,OQP(K,L,NY,NX)
!    2,OSRH(K),OQP(KK,L,NY,NX),OSRH(KK),XOQPZ(K),XOQPZ(KK)
!     ENDIF
          ENDIF
          IF(OQA(K,L,NY,NX)+XOQAZ(K)-XFRA.GT.0.0 &
            .AND.OQA(KK,L,NY,NX)+XOQAZ(KK)+XFRA.GT.0.0)THEN
            XOQAZ(K)=XOQAZ(K)-XFRA
            XOQAZ(KK)=XOQAZ(KK)+XFRA
!     IF((I/1)*1.EQ.I.AND.L.EQ.3.AND.K.EQ.1)THEN
!     WRITE(*,4442)'XOQAZ',I,J,NX,NY,L,K,KK,XFRA,OQA(K,L,NY,NX)
!    2,OSRH(K),OQA(KK,L,NY,NX),OSRH(KK),XOQAZ(K),XOQAZ(KK)
!     ENDIF
          ENDIF
!
!     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
!
!     XFMC,XFMN,XFMP=transfer of microbial C,N,P
!     between each K and KK, FPRIMM=priming transfer rate constant
!     TFNG=temperature+water effect
!     OMC,OMN,OMP=microbial C,N,P
!     OSRH=total SOC in each K
!     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
!
          DO 850 N=1,7
            DO  M=1,3
              DO NGL=1,JG
                XFMC=FPRIMM*TFNG(NGL,N,K)*(OMC(M,NGL,N,K,L,NY,NX)*OSRH(KK) &
                  -OMC(M,NGL,N,KK,L,NY,NX)*OSRH(K))/OSRT
                XFMN=FPRIMM*TFNG(NGL,N,K)*(OMN(M,NGL,N,K,L,NY,NX)*OSRH(KK) &
                  -OMN(M,NGL,N,KK,L,NY,NX)*OSRH(K))/OSRT
                XFMP=FPRIMM*TFNG(NGL,N,K)*(OMP(M,NGL,N,K,L,NY,NX)*OSRH(KK) &
                  -OMP(M,NGL,N,KK,L,NY,NX)*OSRH(K))/OSRT
                IF(OMC(M,NGL,N,K,L,NY,NX)+XOMCZ(M,NGL,N,K)-XFMC.GT.0.0 &
                  .AND.OMC(M,NGL,N,KK,L,NY,NX)+XOMCZ(M,NGL,N,KK)+XFMC.GT.0.0)THEN
                  XOMCZ(M,NGL,N,K)=XOMCZ(M,NGL,N,K)-XFMC
                  XOMCZ(M,NGL,N,KK)=XOMCZ(M,NGL,N,KK)+XFMC
                ENDIF
                IF(OMN(M,NGL,N,K,L,NY,NX)+XOMNZ(M,NGL,N,K)-XFMN.GT.0.0 &
                  .AND.OMN(M,NGL,N,KK,L,NY,NX)+XOMNZ(M,NGL,N,KK)+XFMN.GT.0.0)THEN
                  XOMNZ(M,NGL,N,K)=XOMNZ(M,NGL,N,K)-XFMN
                  XOMNZ(M,NGL,N,KK)=XOMNZ(M,NGL,N,KK)+XFMN
                ENDIF
                IF(OMP(M,NGL,N,K,L,NY,NX)+XOMPZ(M,NGL,N,K)-XFMP.GT.0.0 &
                  .AND.OMP(M,NGL,N,KK,L,NY,NX)+XOMPZ(M,NGL,N,KK)+XFMP.GT.0.0)THEN
                  XOMPZ(M,NGL,N,K)=XOMPZ(M,NGL,N,K)-XFMP
                  XOMPZ(M,NGL,N,KK)=XOMPZ(M,NGL,N,KK)+XFMP
                ENDIF
              enddo
            enddo
850       CONTINUE
        ENDIF
800   CONTINUE
    ENDIF
795   CONTINUE
!
!     TRANSFER ALL PRIMING AMONG ALL K
!
!     TOQCK=total respiration of DOC+DOA in soil layer
!     ROQCK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
  TOQCK(L,NY,NX)=0.0_r8
  DO 840 K=0,KL
    ROQCK(K)=ROQCK(K)+XOQCK(K)
    TOQCK(L,NY,NX)=TOQCK(L,NY,NX)+ROQCK(K)
    OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+XOQCZ(K)
    OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+XOQNZ(K)
    OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+XOQPZ(K)
    OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+XOQAZ(K)
    DO  N=1,7
      DO  M=1,3
        do NGL=1,JG
          OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)+XOMCZ(M,NGL,N,K)
          OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)+XOMNZ(M,NGL,N,K)
          OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)+XOMPZ(M,NGL,N,K)
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,5559)'XOM',I,J,NX,NY,L,K,N,M,OMC(M,NGL,N,K,L,NY,NX)
!    2,OMN(M,NGL,N,K,L,NY,NX),OMP(M,NGL,N,K,L,NY,NX)
!    3,XOMCZ(M,NGL,N,K),XOMNZ(M,NGL,N,K),XOMPZ(M,NGL,N,K)
!5559  FORMAT(A8,8I4,12E12.4)
!     ENDIF
        enddo
      enddo
    enddo
840   CONTINUE
  end associate
  end subroutine OMTransferForPriming
!------------------------------------------------------------------------------------------

  subroutine DOMSorption(K,L,NY,NX,nmicf,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: K,L,NY,NX
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  real(r8) :: AECX
  real(r8) :: OQCX
  real(r8) :: OQNX,OQPX,OQAX
  real(r8) :: OHCX,OHNX,OHPX,OHAX
  real(r8) :: VOLXX
  real(r8) :: VOLXW,VOLCX,VOLCW,VOLAX,VOLAW
!     begin_execution
  associate(                   &
    CGOMN  => nmicf%CGOMN,     &
    CGOMP  => nmicf%CGOMP,     &
    CGOQC  => nmicf%CGOQC,     &
    CGOAC  => nmicf%CGOAC,     &
    CSORP  => ncplxf%CSORP,    &
    ZSORP  => ncplxf%ZSORP,    &
    PSORP  => ncplxf%PSORP,    &
    CSORPA  => ncplxf%CSORPA,  &
    TCGOQC  => ncplxf%TCGOQC,  &
    TCGOAC  => ncplxf%TCGOAC,  &
    TCGOMN  => ncplxf%TCGOMN,  &
    TCGOMP  => ncplxf%TCGOMP,  &
    OSRH  => ncplxs%OSRH,      &
    FOCA  => ncplxs%FOCA,      &
    FOAA  => ncplxs%FOAA       &
  )
!     VOLWM=soil water content, FOSRH=fraction of total SOC
!     AEC,AECX=anion exchange capacity
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     TCGOQC,TCGOMN,TCGOMP,TCGOAC=total uptake of DOC,DON,DOP,acetate
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     TSORP,HSORP=sorption rate constant and coefficient for OHC
!     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
  IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
    IF(L.EQ.0)THEN
      AECX=0.5E+03
    ELSE
      AECX=AEC(L,NY,NX)
    ENDIF
    OQCX=AMAX1(ZEROS(NY,NX),OQC(K,L,NY,NX)-TCGOQC(K))
    OQNX=AMAX1(ZEROS(NY,NX),OQN(K,L,NY,NX)-TCGOMN(K))
    OQPX=AMAX1(ZEROS(NY,NX),OQP(K,L,NY,NX)-TCGOMP(K))
    OQAX=AMAX1(ZEROS(NY,NX),OQA(K,L,NY,NX)-TCGOAC(K))
    OHCX=AMAX1(ZEROS(NY,NX),OHC(K,L,NY,NX))
    OHNX=AMAX1(ZEROS(NY,NX),OHN(K,L,NY,NX))
    OHPX=AMAX1(ZEROS(NY,NX),OHP(K,L,NY,NX))
    OHAX=AMAX1(ZEROS(NY,NX),OHA(K,L,NY,NX))
    VOLXX=BKVL(L,NY,NX)*AECX*HSORP*FOSRH(K,L,NY,NX)
    VOLXW=VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)
    IF(FOCA(K).GT.ZERO)THEN
      VOLCX=FOCA(K)*VOLXX
      VOLCW=FOCA(K)*VOLXW
      CSORP(K)=TSORP*(OQCX*VOLCX-OHCX*VOLCW)/(VOLCX+VOLCW)
    ELSE
      CSORP(K)=TSORP*(OQCX*VOLXX-OHCX*VOLXW)/(VOLXX+VOLXW)
    ENDIF
    IF(FOAA(K).GT.ZERO)THEN
      VOLAX=FOAA(K)*VOLXX
      VOLAW=FOAA(K)*VOLXW
      CSORPA(K)=TSORP*(OQAX*VOLAX-OHAX*VOLAW)/(VOLAX+VOLAW)
    ELSE
      CSORPA(K)=TSORP*(OQAX*VOLXX-OHAX*VOLXW)/(VOLXX+VOLXW)
    ENDIF
    ZSORP(K)=TSORP*(OQNX*VOLXX-OHNX*VOLXW)/(VOLXX+VOLXW)
    PSORP(K)=TSORP*(OQPX*VOLXX-OHPX*VOLXW)/(VOLXX+VOLXW)
  ELSE
    CSORP(K)=0.0_r8
    CSORPA(K)=0.0_r8
    ZSORP(K)=0.0_r8
    PSORP(K)=0.0_r8
  ENDIF
  end associate
  end subroutine DOMSorption
!------------------------------------------------------------------------------------------

  subroutine SolidOMDecomposition(K,L,NY,NX,nmicdiag,ncplxf,ncplxs)

  implicit none
  integer, intent(in) :: K,L,NY,NX
  type(NitroMicDiagType), intent(inout) :: nmicdiag
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  integer  :: M
  real(r8) :: CNOMX,CPOMX
  real(r8) :: COQCK,COSC
  real(r8) :: CPR,CNR
  real(r8) :: DCKD
  real(r8) :: DFNS
  real(r8) :: OQCI
  real(r8) :: RHOSCM
  real(r8) :: FCNK(0:jcplx1),FCPK(0:jcplx1)
  real(r8) :: CNS(4,0:jcplx1),CPS(4,0:jcplx1)

!     begin_execution
  associate(                 &
    RDOSC  => ncplxf%RDOSC,  &
    RDOSN  => ncplxf%RDOSN,  &
    RDOSP  => ncplxf%RDOSP,  &
    RHOSC  => ncplxf%RHOSC,  &
    RHOSN  => ncplxf%RHOSN,  &
    RHOSP  => ncplxf%RHOSP,  &
    RCOSC  => ncplxf%RCOSC,  &
    RCOSN  => ncplxf%RCOSN,  &
    RCOSP  => ncplxf%RCOSP,  &
    RDORC  => ncplxf%RDORC,  &
    RDORN  => ncplxf%RDORN,  &
    RDORP  => ncplxf%RDORP,  &
    RDOHC  => ncplxf%RDOHC,  &
    RDOHN  => ncplxf%RDOHN,  &
    RDOHP  => ncplxf%RDOHP,  &
    RDOHA  => ncplxf%RDOHA,  &
    ROQCK  => ncplxf%ROQCK,  &
    OSRH  => ncplxs%OSRH,    &
    TOMK  => ncplxs%TOMK,    &
    TONK  => ncplxs%TONK,    &
    TOPK  => ncplxs%TOPK,    &
    CNH  => ncplxs%CNH,      &
    CPH  => ncplxs%CPH,      &
    TONX  => ncplxs%TONX,    &
    TOPX  => ncplxs%TOPX,    &
    VOLWZ =>  nmicdiag%VOLWZ,&
    TFNX  =>  nmicdiag%TFNX  &
  )
!     FCPK=N,P limitation to microbial activity in each K
!     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
!     COQCK=aqueous concentration of microbial activity
!     DCKD=Km for decomposition of SOC at current COQCK
!     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
!     DCKI=inhibition of decomposition by microbial concentration
!     OSRH=total SOC
!     COSC=concentration of total SOC
!     BKVL,VOLX=mass, volume of soil layer
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     OQKI=DOC product inhibition constant for decomposition
!
  IF(TOMK(K).GT.ZEROS(NY,NX))THEN
    CNOMX=TONK(K)/TONX(K)
    CPOMX=TOPK(K)/TOPX(K)
    FCNK(K)=AMIN1(1.0,AMAX1(0.50,CNOMX))
    FCPK(K)=AMIN1(1.0,AMAX1(0.50,CPOMX))
  ELSE
    FCNK(K)=1.0
    FCPK(K)=1.0
  ENDIF
!
!     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
!     CONSTANT FOR DECOMPOSITION
!
  IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
    COQCK=AMIN1(0.1E+06,ROQCK(K)/VOLWZ)
  ELSE
    COQCK=0.1E+06
  ENDIF
  IF(L.EQ.0)THEN
    DCKD=DCKM0*(1.0+COQCK/DCKI)
  ELSE
    DCKD=DCKML*(1.0+COQCK/DCKI)
  ENDIF
  IF(OSRH(K).GT.ZEROS(NY,NX))THEN
    IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      COSC=OSRH(K)/BKVL(L,NY,NX)
    ELSE
      COSC=OSRH(K)/VOLY(L,NY,NX)
    ENDIF
    DFNS=COSC/(COSC+DCKD)
    OQCI=1.0/(1.0+COQC(K,L,NY,NX)/OQKI)
!     IF(L.EQ.0.AND.J.EQ.15)THEN
!     WRITE(*,4242)'COSC',I,J,L,K,DFNS,COSC,COQCK,DCKD,OSRH(K)
!    2,OSAT(K),OSCT(K),ORCT(K),OHC(K,L,NY,NX),BKVL(L,NY,NX),ROQCK(K)
!    3,VOLWZ,VOLWRX(NY,NX),VOLW(0,NY,NX),FCR(NY,NX)
!    4,THETY(L,NY,NX)
!4242  FORMAT(A8,4I4,30E12.4)
!     ENDIF
!
!     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES 'RDOS*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     CNS,CPS=N:C,P:C ratios of SOC
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OSA,OSN,OSP=active biomass C,N,P
!     SPOSC=specific decomposition rate constant
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
    DO 785 M=1,jsken
      IF(OSC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
        CNS(M,K)=AMAX1(0.0,OSN(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
        CPS(M,K)=AMAX1(0.0,OSP(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
        RDOSC(M,K)=AMAX1(0.0,AMIN1(0.5*OSA(M,K,L,NY,NX) &
          ,SPOSC(M,K)*ROQCK(K)*DFNS*OQCI*TFNX*OSA(M,K,L,NY,NX)/OSRH(K)))
        RDOSN(M,K)=AMAX1(0.0,AMIN1(OSN(M,K,L,NY,NX),CNS(M,K)*RDOSC(M,K)))/FCNK(K)
        RDOSP(M,K)=AMAX1(0.0,AMIN1(OSP(M,K,L,NY,NX),CPS(M,K)*RDOSC(M,K)))/FCPK(K)

      ELSE
        CNS(M,K)=CNOSC(M,K,L,NY,NX)
        CPS(M,K)=CPOSC(M,K,L,NY,NX)
        RDOSC(M,K)=0.0_r8
        RDOSN(M,K)=0.0_r8
        RDOSP(M,K)=0.0_r8
      ENDIF
785 CONTINUE
!
!     HUMIFICATION OF DECOMPOSED RESIDUE LIGNIN WITH PROTEIN,
!     CH2O AND CELLULOSE 'RHOS*' WITH REMAINDER 'RCOS*' TO DOC,DON,DOP
!
!     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP
!     CNRH,CPRH=N:C,P:C in POC
!     EPOC=fraction of RDOSC allocated to POC from hour1.f
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
    IF(K.LE.2)THEN
      RHOSC(4,K)=AMAX1(0.0,AMIN1(RDOSN(4,K)/CNRH(3) &
        ,RDOSP(4,K)/CPRH(3),EPOC(L,NY,NX)*RDOSC(4,K)))
      RHOSCM=0.10*RHOSC(4,K)
      RHOSC(1,K)=AMAX1(0.0,AMIN1(RDOSC(1,K),RDOSN(1,K)/CNRH(3) &
        ,RDOSP(1,K)/CPRH(3),RHOSCM))
      RHOSC(2,K)=AMAX1(0.0,AMIN1(RDOSC(2,K),RDOSN(2,K)/CNRH(3) &
        ,RDOSP(2,K)/CPRH(3),RHOSCM))
      RHOSC(3,K)=AMAX1(0.0,AMIN1(RDOSC(3,K),RDOSN(3,K)/CNRH(3) &
        ,RDOSP(3,K)/CPRH(3),RHOSCM-RHOSC(2,K)))
      DO 805 M=1,jsken
        RHOSN(M,K)=AMIN1(RDOSN(M,K),RHOSC(M,K)*CNRH(3))
        RHOSP(M,K)=AMIN1(RDOSP(M,K),RHOSC(M,K)*CPRH(3))
        RCOSC(M,K)=RDOSC(M,K)-RHOSC(M,K)
        RCOSN(M,K)=RDOSN(M,K)-RHOSN(M,K)
        RCOSP(M,K)=RDOSP(M,K)-RHOSP(M,K)
805   CONTINUE
    ELSE
      DO 810 M=1,jsken
        RHOSC(M,K)=0.0_r8
        RHOSN(M,K)=0.0_r8
        RHOSP(M,K)=0.0_r8
        RCOSC(M,K)=RDOSC(M,K)
        RCOSN(M,K)=RDOSN(M,K)
        RCOSP(M,K)=RDOSP(M,K)
810   CONTINUE
    ENDIF
  ELSE
    DO 780 M=1,jsken
      RDOSC(M,K)=0.0_r8
      RDOSN(M,K)=0.0_r8
      RDOSP(M,K)=0.0_r8
      RHOSC(M,K)=0.0_r8
      RHOSN(M,K)=0.0_r8
      RHOSP(M,K)=0.0_r8
      RCOSC(M,K)=0.0_r8
      RCOSN(M,K)=0.0_r8
      RCOSP(M,K)=0.0_r8
780 CONTINUE
  ENDIF
!
!     C, N, P DECOMPOSITION RATE OF BIORESIDUE 'RDOR*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     CNR,CPR=N:C,P:C ratios of microbial residue
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     SPORC=specific decomposition rate constant for microbial residue
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
  IF(OSRH(K).GT.ZEROS(NY,NX))THEN
    DO 775 M=1,2
      IF(ORC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
        CNR=AMAX1(0.0,ORN(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
        CPR=AMAX1(0.0,ORP(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
        RDORC(M,K)=AMAX1(0.0,AMIN1(ORC(M,K,L,NY,NX) &
          ,SPORC(M)*ROQCK(K)*DFNS*OQCI*TFNX*ORC(M,K,L,NY,NX)/OSRH(K)))
    !    3*AMIN1(FCNK(K),FCPK(K))
        RDORN(M,K)=AMAX1(0.0,AMIN1(ORN(M,K,L,NY,NX),CNR*RDORC(M,K)))/FCNK(K)
        RDORP(M,K)=AMAX1(0.0,AMIN1(ORP(M,K,L,NY,NX),CPR*RDORC(M,K)))/FCPK(K)
      ELSE
        RDORC(M,K)=0.0_r8
        RDORN(M,K)=0.0_r8
        RDORP(M,K)=0.0_r8
      ENDIF
775 CONTINUE
  ELSE
    DO 776 M=1,2
      RDORC(M,K)=0.0_r8
      RDORN(M,K)=0.0_r8
      RDORP(M,K)=0.0_r8
776 CONTINUE
  ENDIF
!
!     C, N, P DECOMPOSITION RATE OF SORBED SUBSTRATES 'RDOH*' FROM
!     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
!     TEMPERATURE, SUBSTRATE C:N, C:P
!
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     CNH,CPH=N:C,P:C ratios of adsorbed C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     SPOHC=specific decomposition rate constant for adsorbed C
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     TFNX=temperature stress effect
!     OSRH=total SOC
!     FCNK,FCPK=N,P limitation to microbial activity in each K
!
  IF(OSRH(K).GT.ZEROS(NY,NX))THEN
    IF(OHC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNH(K)=AMAX1(0.0,OHN(K,L,NY,NX)/OHC(K,L,NY,NX))
      CPH(K)=AMAX1(0.0,OHP(K,L,NY,NX)/OHC(K,L,NY,NX))
      RDOHC(K)=AMAX1(0.0,AMIN1(OHC(K,L,NY,NX) &
        ,SPOHC*ROQCK(K)*DFNS*OQCI*TFNX*OHC(K,L,NY,NX)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
      RDOHN(K)=AMAX1(0.0,AMIN1(OHN(K,L,NY,NX),CNH(K)*RDOHC(K)))/FCNK(K)
      RDOHP(K)=AMAX1(0.0,AMIN1(OHP(K,L,NY,NX),CPH(K)*RDOHC(K)))/FCPK(K)
      RDOHA(K)=AMAX1(0.0,AMIN1(OHA(K,L,NY,NX) &
        ,SPOHA*ROQCK(K)*DFNS*TFNX*OHA(K,L,NY,NX)/OSRH(K)))
!    3*AMIN1(FCNK(K),FCPK(K))
    ELSE
      CNH(K)=0.0_r8
      CPH(K)=0.0_r8
      RDOHC(K)=0.0_r8
      RDOHN(K)=0.0_r8
      RDOHP(K)=0.0_r8
      RDOHA(K)=0.0_r8
    ENDIF
  ELSE
    CNH(K)=0.0_r8
    CPH(K)=0.0_r8
    RDOHC(K)=0.0_r8
    RDOHN(K)=0.0_r8
    RDOHP(K)=0.0_r8
    RDOHA(K)=0.0_r8
  ENDIF
  end associate
  end subroutine SolidOMDecomposition
!------------------------------------------------------------------------------------------

  subroutine RedistDecompositionProduct(KL,L,NY,NX,nmicdiag,nmicf,ncplxf,ncplxs)

  implicit none
  integer, intent(in) :: KL,L,NY,NX
  type(NitroMicDiagType), intent(in) :: nmicdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout):: ncplxs
  integer :: K,M,N,NGL
  real(r8) :: FORC(0:jcplx)
!     begin_execution
  associate(                   &
    CGOMN  => nmicf%CGOMN,     &
    CGOMP  => nmicf%CGOMP,     &
    CGOQC  => nmicf%CGOQC,     &
    CGOAC  => nmicf%CGOAC,     &
    RCOMC  => nmicf%RCOMC,     &
    RCOMN  => nmicf%RCOMN,     &
    RCOMP  => nmicf%RCOMP,     &
    RCMMC  => nmicf%RCMMC,     &
    RCMMN  => nmicf%RCMMN,     &
    RCMMP  => nmicf%RCMMP,     &
    RCCMC  => nmicf%RCCMC,     &
    RCCMN  => nmicf%RCCMN,     &
    RCCMP  => nmicf%RCCMP,     &
    RCH3X  => nmicf%RCH3X,     &
    RDOSC  => ncplxf%RDOSC,  &
    RDOSN  => ncplxf%RDOSN,  &
    RDOSP  => ncplxf%RDOSP,  &
    RHOSC  => ncplxf%RHOSC,  &
    RHOSN  => ncplxf%RHOSN,  &
    RHOSP  => ncplxf%RHOSP,  &
    RCOSC  => ncplxf%RCOSC,  &
    RCOSN  => ncplxf%RCOSN,  &
    RCOSP  => ncplxf%RCOSP,  &
    RDORC  => ncplxf%RDORC,  &
    RDORN  => ncplxf%RDORN,  &
    RDORP  => ncplxf%RDORP,  &
    RDOHC  => ncplxf%RDOHC,  &
    RDOHN  => ncplxf%RDOHN,  &
    RDOHP  => ncplxf%RDOHP,  &
    RDOHA  => ncplxf%RDOHA,  &
    CSORP  => ncplxf%CSORP,  &
    ZSORP  => ncplxf%ZSORP,  &
    PSORP  => ncplxf%PSORP,  &
    CSORPA  => ncplxf%CSORPA,&
    ORCT  => ncplxs%ORCT    ,&
    RCOQN =>  nmicdiag%RCOQN,&
    TORC  =>  nmicdiag%TORC  &
  )
!
!     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
!     HETEROTROPHIC SUBSTRATE-MICROBE COMPLEXES
!
!     FORC=fraction of total microbial residue
!     ORCT=microbial residue
!     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
  DO 1690 K=0,KL
    IF(TORC.GT.ZEROS(NY,NX))THEN
      FORC(K)=ORCT(K)/TORC
    ELSE
      IF(K.EQ.3)THEN
        FORC(K)=1.0
      ELSE
        FORC(K)=0.0_r8
      ENDIF
    ENDIF
    DO 1685 N=1,7
      DO 1680 M=1,2
        DO NGL=1,JG
          RCCMC(M,NGL,N,K)=(RCOMC(M,NGL,N,5)+RCMMC(M,NGL,N,5))*FORC(K)
          RCCMN(M,NGL,N,K)=(RCOMN(M,NGL,N,5)+RCMMN(M,NGL,N,5))*FORC(K)
          RCCMP(M,NGL,N,K)=(RCOMP(M,NGL,N,5)+RCMMP(M,NGL,N,5))*FORC(K)
        ENDDO
!     IF(L.EQ.0)THEN
!     WRITE(*,8821)'RCCMC',I,J,L,K,N,M,RCCMC(M,NGL,N,K)
!    2,RCOMC(M,N,5),RCMMC(M,N,5),FORC(K)
!     ENDIF
1680  CONTINUE
1685  CONTINUE
1690  CONTINUE
!
!   REDISTRIBUTE C,N AND P TRANSFORMATIONS AMONG STATE
!   VARIABLES IN SUBSTRATE-MICROBE COMPLEXES
!

  DO 590 K=0,KL
    DO 580 M=1,jsken
!
!     SUBSTRATE DECOMPOSITION PRODUCTS
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
!     OQC,OQN,OQP,OQA=DOC,DON,DOP
!     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
!
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-RDOSC(M,K)
!     OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-RDOSC(M,K)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-RDOSN(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-RDOSP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RCOSC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RCOSN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RCOSP(M,K)
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.K.EQ.4)THEN
!     WRITE(*,4444)'RDOSC',I,J,NX,NY,L,K,M,OSC(M,K,L,NY,NX)
!    2,RDOSC(M,K)
!     ENDIF
!
!     LIGNIFICATION PRODUCTS
!
!       RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
!
      IF(L.NE.0)THEN
        OSC(1,3,L,NY,NX)=OSC(1,3,L,NY,NX)+RHOSC(M,K)
    !     OSA(1,3,L,NY,NX)=OSA(1,3,L,NY,NX)+RHOSC(M,K)
        OSN(1,3,L,NY,NX)=OSN(1,3,L,NY,NX)+RHOSN(M,K)
        OSP(1,3,L,NY,NX)=OSP(1,3,L,NY,NX)+RHOSP(M,K)
      ELSE
        OSC(1,3,NU(NY,NX),NY,NX)=OSC(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
    !     OSA(1,3,NU(NY,NX),NY,NX)=OSA(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
        OSN(1,3,NU(NY,NX),NY,NX)=OSN(1,3,NU(NY,NX),NY,NX)+RHOSN(M,K)
        OSP(1,3,NU(NY,NX),NY,NX)=OSP(1,3,NU(NY,NX),NY,NX)+RHOSP(M,K)
      ENDIF
580 CONTINUE
!
!     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
!     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
!     RCOQN=DON production from nitrous acid reduction
!
    DO 575 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-RDORC(M,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-RDORN(M,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-RDORP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RDORC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RDORN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RDORP(M,K)
575 CONTINUE
    OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RDOHC(K)
    OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RDOHN(K)+RCOQN*FORC(K)
    OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RDOHP(K)
    OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+RDOHA(K)
    OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-RDOHC(K)
    OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-RDOHN(K)
    OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-RDOHP(K)
    OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-RDOHA(K)
!
!     MICROBIAL UPTAKE OF DISSOLVED C, N, P
!
!     CGOQC,CGOAC,CGOMN,CGOMP=DOC,acetate,DON,DOP uptake
!     RCH3X=acetate production from fermentation
!
    DO 570 N=1,7
      DO NGL=1,JG
        OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CGOQC(NGL,N,K)
        OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-CGOMN(NGL,N,K)
        OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-CGOMP(NGL,N,K)
        OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CGOAC(NGL,N,K)+RCH3X(NGL,N,K)
!
!     MICROBIAL DECOMPOSITION PRODUCTS
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
        DO 565 M=1,2
          ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+RCOMC(M,NGL,N,K)+RCCMC(M,NGL,N,K) &
            +RCMMC(M,NGL,N,K)
          ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+RCOMN(M,NGL,N,K)+RCCMN(M,NGL,N,K) &
            +RCMMN(M,NGL,N,K)
          ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+RCOMP(M,NGL,N,K)+RCCMP(M,NGL,N,K) &
            +RCMMP(M,NGL,N,K)
!     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4.AND.K.EQ.2)THEN
!     WRITE(*,8821)'ORC',I,J,L,K,N,M,ORC(M,K,L,NY,NX)
!    2,RCOMC(M,NGL,N,K),RCCMC(M,NGL,N,K),RCMMC(M,NGL,N,K),RDORC(M,K)
!     WRITE(*,8821)'ORP',I,J,L,K,N,M,ORP(M,K,L,NY,NX)
!    2,RCOMP(M,NGL,N,K),RCCMP(M,NGL,N,K),RCMMP(M,NGL,N,K),RDORP(M,K)
!8821  FORMAT(A8,6I4,40E12.4)
!     ENDIF
565     CONTINUE
      enddo
570 CONTINUE
!
!     SORPTION PRODUCTS
!
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
    OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CSORP(K)
    OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ZSORP(K)
    OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-PSORP(K)
    OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CSORPA(K)
    OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+CSORP(K)
    OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+ZSORP(K)
    OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+PSORP(K)
    OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+CSORPA(K)

590 CONTINUE
  end associate
  end subroutine RedistDecompositionProduct
!------------------------------------------------------------------------------------------

  subroutine MicrobialAnabolicUpdate(L,NY,NX,nmicf)
  implicit none
  integer, intent(in) :: L,NY,NX
  type(NitroMicFluxType), intent(inout) :: nmicf
  integer  :: K,M,N,NGL
  real(r8) ::CGROMC
!     begin_execution
  associate(                  &
    CGOMC    => nmicf%CGOMC,  &
    CGOMN  => nmicf%CGOMN,    &
    CGOMP  => nmicf%CGOMP,    &
    CGOMS  => nmicf%CGOMS,    &
    CGONS  => nmicf%CGONS,    &
    CGOPS  => nmicf%CGOPS,    &
    RGN2F  => nmicf%RGN2F,    &
    RGOMO  => nmicf%RGOMO,    &
    RGOMD  => nmicf%RGOMD,    &
    RINO3  => nmicf%RINO3,    &
    RCO2X  => nmicf%RCO2X,    &
    RIPO4  => nmicf%RIPO4,    &
    RINB4 => nmicf%RINB4 ,    &
    RINB3 => nmicf%RINB3 ,    &
    RIPOB => nmicf%RIPOB ,    &
    RHOMC => nmicf%RHOMC ,    &
    RHOMN => nmicf%RHOMN ,    &
    RHOMP => nmicf%RHOMP ,    &
    RHMMC  => nmicf%RHMMC,    &
    RHMMN  => nmicf%RHMMN,    &
    RHMMP  => nmicf%RHMMP,    &
    RN2FX  => nmicf%RN2FX,    &
    RXOMC  => nmicf%RXOMC,    &
    RXOMN  => nmicf%RXOMN,    &
    RXOMP  => nmicf%RXOMP,    &
    R3OMC  => nmicf%R3OMC,    &
    R3OMN  => nmicf%R3OMN,    &
    R3OMP  => nmicf%R3OMP,    &
    RXMMC  => nmicf%RXMMC,    &
    RXMMN  => nmicf%RXMMN,    &
    RXMMP  => nmicf%RXMMP,    &
    R3MMC  => nmicf%R3MMC,    &
    R3MMN  => nmicf%R3MMN,    &
    R3MMP  => nmicf%R3MMP,    &
    RINH4R  => nmicf%RINH4R,  &
    RINO3R  => nmicf%RINO3R,  &
    RIPO4R  => nmicf%RIPO4R,  &
    RIP14   => nmicf%RIP14 ,  &
    RIP1B   => nmicf%RIP1B ,  &
    RIP14R  => nmicf%RIP14R,  &
    RINH4   => nmicf%RINH4    &
  )
!
!     OMC,OMN,OMP=microbial C,N,P
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
!     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
!

  DO 550 K=0,5
    IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 545 N=1,7
        IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
          DO NGL=1,JG
            DO 540 M=1,2
              OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)+CGOMS(M,NGL,N,K) &
                -RXOMC(M,NGL,N,K)-RXMMC(M,NGL,N,K)
              OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)+CGONS(M,NGL,N,K) &
                -RXOMN(M,NGL,N,K)-RXMMN(M,NGL,N,K)
              OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)+CGOPS(M,NGL,N,K) &
                -RXOMP(M,NGL,N,K)-RXMMP(M,NGL,N,K)

!
!     HUMIFICATION PRODUCTS
!
!     CFOMC=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!
              IF(L.NE.0)THEN
                OSC(1,4,L,NY,NX)=OSC(1,4,L,NY,NX)+CFOMC(1,L,NY,NX) &
                  *(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
!     OSA(1,4,L,NY,NX)=OSA(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
!    2*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
                OSN(1,4,L,NY,NX)=OSN(1,4,L,NY,NX)+CFOMC(1,L,NY,NX) &
                  *(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
                OSP(1,4,L,NY,NX)=OSP(1,4,L,NY,NX)+CFOMC(1,L,NY,NX) &
                  *(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
                OSC(2,4,L,NY,NX)=OSC(2,4,L,NY,NX)+CFOMC(2,L,NY,NX) &
                  *(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
!     OSA(2,4,L,NY,NX)=OSA(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
!    2*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
                OSN(2,4,L,NY,NX)=OSN(2,4,L,NY,NX)+CFOMC(2,L,NY,NX) &
                  *(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
                OSP(2,4,L,NY,NX)=OSP(2,4,L,NY,NX)+CFOMC(2,L,NY,NX) &
                  *(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
!     IF((I/10)*10.EQ.I.AND.J.EQ.24)THEN
!     WRITE(*,4445)'RHOMC',I,J,NX,NY,L,K,M,N,OSC(1,4,L,NY,NX)
!    2,OSC(2,4,L,NY,NX),CFOMC(1,L,NY,NX),CFOMC(2,L,NY,NX)
!    3,RHOMC(M,NGL,N,K),RHMMC(M,NGL,N,K)
!4445  FORMAT(A8,8I4,40E12.4)
!     ENDIF
              ELSE
                OSC(1,4,NU(NY,NX),NY,NX)=OSC(1,4,NU(NY,NX),NY,NX) &
                  +CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
!     OSA(1,4,NU(NY,NX),NY,NX)=OSA(1,4,NU(NY,NX),NY,NX)
!    2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
                OSN(1,4,NU(NY,NX),NY,NX)=OSN(1,4,NU(NY,NX),NY,NX) &
                  +CFOMC(1,NU(NY,NX),NY,NX)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
                OSP(1,4,NU(NY,NX),NY,NX)=OSP(1,4,NU(NY,NX),NY,NX) &
                  +CFOMC(1,NU(NY,NX),NY,NX)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
                OSC(2,4,NU(NY,NX),NY,NX)=OSC(2,4,NU(NY,NX),NY,NX) &
                  +CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
!     OSA(2,4,NU(NY,NX),NY,NX)=OSA(2,4,NU(NY,NX),NY,NX)
!    2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
                OSN(2,4,NU(NY,NX),NY,NX)=OSN(2,4,NU(NY,NX),NY,NX) &
                  +CFOMC(2,NU(NY,NX),NY,NX)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
                OSP(2,4,NU(NY,NX),NY,NX)=OSP(2,4,NU(NY,NX),NY,NX) &
                  +CFOMC(2,NU(NY,NX),NY,NX)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
              ENDIF
540         CONTINUE

!
!     INPUTS TO NONSTRUCTURAL POOLS
!
!     CGOMC=total DOC+acetate uptake
!     RGOMO=total respiration
!     RGOMD=respiration for denitrifcation
!     RGN2F=respiration for N2 fixation
!     RCO2X=total CO2 emission
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
!     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!     CGOMN,CGOMP=DON, DOP uptake
!     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
!     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
!     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
!     RINH4R,RINO3R =substrate-limited NH4,NO3 mineraln-immobiln
!     RIPO4R,RIP14R=substrate-limited H2PO4,HPO4 mineraln-immobiln
!
            CGROMC=CGOMC(NGL,N,K)-RGOMO(NGL,N,K)-RGOMD(NGL,N,K)-RGN2F(NGL,N,K)
            RCO2X(NGL,N,K)=RCO2X(NGL,N,K)+RGN2F(NGL,N,K)
            DO 555 M=1,2
              OMC(3,NGL,N,K,L,NY,NX)=OMC(3,NGL,N,K,L,NY,NX)-CGOMS(M,NGL,N,K)+R3OMC(M,NGL,N,K)
              OMN(3,NGL,N,K,L,NY,NX)=OMN(3,NGL,N,K,L,NY,NX)-CGONS(M,NGL,N,K)+R3OMN(M,NGL,N,K)+R3MMN(M,NGL,N,K)
              OMP(3,NGL,N,K,L,NY,NX)=OMP(3,NGL,N,K,L,NY,NX)-CGOPS(M,NGL,N,K)+R3OMP(M,NGL,N,K)+R3MMP(M,NGL,N,K)
              RCO2X(NGL,N,K)=RCO2X(NGL,N,K)+R3MMC(M,NGL,N,K)
555         CONTINUE
            OMC(3,NGL,N,K,L,NY,NX)=OMC(3,NGL,N,K,L,NY,NX)+CGROMC
            OMN(3,NGL,N,K,L,NY,NX)=OMN(3,NGL,N,K,L,NY,NX)+CGOMN(NGL,N,K) &
              +RINH4(NGL,N,K)+RINB4(NGL,N,K)+RINO3(NGL,N,K)+RINB3(NGL,N,K)+RN2FX(NGL,N,K)
            OMP(3,NGL,N,K,L,NY,NX)=OMP(3,NGL,N,K,L,NY,NX)+CGOMP(NGL,N,K) &
              +RIPO4(NGL,N,K)+RIPOB(NGL,N,K)+RIP14(NGL,N,K)+RIP1B(NGL,N,K)
            IF(L.EQ.0)THEN
              OMN(3,NGL,N,K,L,NY,NX)=OMN(3,NGL,N,K,L,NY,NX)+RINH4R(NGL,N,K)+RINO3R(NGL,N,K)
              OMP(3,NGL,N,K,L,NY,NX)=OMP(3,NGL,N,K,L,NY,NX)+RIPO4R(NGL,N,K)+RIP14R(NGL,N,K)
            ENDIF
          enddo

        ENDIF
545   CONTINUE
    ENDIF
550 CONTINUE
  end associate
  end subroutine MicrobialAnabolicUpdate
!------------------------------------------------------------------------------------------

  subroutine MicrobialLitterColonization(KL,L,NY,NX,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: KL,L,NY,NX
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  integer  :: K,M
  real(r8) :: DOSAK
!     begin_execution
  associate(                 &
    ROQCK  => ncplxf%ROQCK,  &
    OSCT  => ncplxs%OSCT,    &
    OSAT  => ncplxs%OSAT     &
  )
!     OSCT,OSAT,OSCX=total,colonized,uncolonized SOC
!     OSA,OSC=colonized,total litter
!     DOSA=rate constant for litter colonization
!     ROQCK=total respiration of DOC+DOA used to represent microbial activity
!
  DO 475 K=0,KL
    OSCT(K)=0.0_r8
    OSAT(K)=0.0_r8
    DO  M=1,jsken
      OSCT(K)=OSCT(K)+OSC(M,K,L,NY,NX)
      OSAT(K)=OSAT(K)+OSA(M,K,L,NY,NX)
    enddo
475 CONTINUE
  DO 480 K=0,KL
    IF(OSCT(K).GT.ZEROS(NY,NX))THEN
      DOSAK=DOSA(K)*AMAX1(0.0,ROQCK(K))
      DO 485 M=1,jsken
        OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX) &
          ,OSA(M,K,L,NY,NX)+DOSAK*OSC(M,K,L,NY,NX)/OSCT(K))
485   CONTINUE
    ELSE
      DO 490 M=1,jsken
        OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX),OSA(M,K,L,NY,NX))
490   CONTINUE
    ENDIF

480 CONTINUE
  end associate
  end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

  subroutine AggregateTransformations(L,NY,NX,nmicdiag,naqfdiag,nmicf,ncplxf)
  implicit none
  integer, intent(in) :: L,NY,NX
  type(NitroMicDiagType),intent(in) :: nmicdiag
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  integer  :: K,M,N,NGL
!     begin_execution
  associate(                   &
    CGOMC  => nmicf%CGOMC,     &
    CGOMN  => nmicf%CGOMN,     &
    CGOMP  => nmicf%CGOMP,     &
    CGOQC  => nmicf%CGOQC,     &
    CGOAC  => nmicf%CGOAC,     &
    RUPOX  => nmicf%RUPOX,     &
    RDNO3 => nmicf%RDNO3 ,     &
    RDNOB => nmicf%RDNOB ,     &
    RDNO2 => nmicf%RDNO2 ,     &
    RDN2B => nmicf%RDN2B ,     &
    RDN2O => nmicf%RDN2O ,     &
    RGOMD  => nmicf%RGOMD,     &
    RINH4 => nmicf%RINH4 ,     &
    RINO3  => nmicf%RINO3,     &
    RIPO4  => nmicf%RIPO4,     &
    RINB4 => nmicf%RINB4 ,     &
    RINB3 => nmicf%RINB3 ,     &
    RIPOB => nmicf%RIPOB ,     &
    RH2GX  => nmicf%RH2GX,     &
    RINH4R  => nmicf%RINH4R,   &
    RINO3R  => nmicf%RINO3R,   &
    RIPO4R  => nmicf%RIPO4R,   &
    RCO2X  => nmicf%RCO2X  ,   &
    RCH3X  => nmicf%RCH3X ,    &
    RCH4X  => nmicf%RCH4X ,    &
    RVOXA  => nmicf%RVOXA,     &
    RVOXB  => nmicf%RVOXB,     &
    RIP14  => nmicf%RIP14,     &
    RIP1B  => nmicf%RIP1B,     &
    RIP14R  => nmicf%RIP14R,   &
    RN2FX  => nmicf%RN2FX  ,   &
    RCOSC  => ncplxf%RCOSC,    &
    RCOSN  => ncplxf%RCOSN,    &
    RCOSP  => ncplxf%RCOSP,    &
    RDORC  => ncplxf%RDORC,    &
    RDORN  => ncplxf%RDORN,    &
    RDORP  => ncplxf%RDORP,    &
    RDOHC  => ncplxf%RDOHC,    &
    RDOHN  => ncplxf%RDOHN,    &
    RDOHP  => ncplxf%RDOHP,    &
    RDOHA  => ncplxf%RDOHA,    &
    CSORP  => ncplxf%CSORP,    &
    ZSORP  => ncplxf%ZSORP,    &
    PSORP  => ncplxf%PSORP,    &
    CSORPA  => ncplxf%CSORPA,  &
    TFNX  =>  nmicdiag%TFNX,   &
    RH2GZ  =>  nmicdiag%RH2GZ, &
    RCNO2  =>  nmicdiag%RCNO2, &
    RCNOB  =>  nmicdiag%RCNOB, &
    RCN2O  =>  nmicdiag%RCN2O, &
    RCN2B  =>  nmicdiag%RCN2B, &
    RCNO3  =>  nmicdiag%RCNO3, &
    RCN3B  =>  nmicdiag%RCN3B, &
    VOLWZ  =>  nmicdiag%VOLWZ  &
  )
  DO 650 K=0,5
    IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 640 N=1,7
        IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
          DO NGL=1,JG
            naqfdiag%TRINH=naqfdiag%TRINH+RINH4(NGL,N,K)
            naqfdiag%TRINO=naqfdiag%TRINO+RINO3(NGL,N,K)
            naqfdiag%TRIPO=naqfdiag%TRIPO+RIPO4(NGL,N,K)
            naqfdiag%TRIP1=naqfdiag%TRIP1+RIP14(NGL,N,K)
            naqfdiag%TRINB=naqfdiag%TRINB+RINB4(NGL,N,K)
            naqfdiag%TRIOB=naqfdiag%TRIOB+RINB3(NGL,N,K)
            naqfdiag%TRIPB=naqfdiag%TRIPB+RIPOB(NGL,N,K)
            naqfdiag%TRIB1=naqfdiag%TRIB1+RIP1B(NGL,N,K)
            naqfdiag%TRN2F=naqfdiag%TRN2F+RN2FX(NGL,N,K)
            IF(L.EQ.NU(NY,NX))THEN
              naqfdiag%TRINH=naqfdiag%TRINH+RINH4R(NGL,N,K)
              naqfdiag%TRINO=naqfdiag%TRINO+RINO3R(NGL,N,K)
              naqfdiag%TRIPO=naqfdiag%TRIPO+RIPO4R(NGL,N,K)
              naqfdiag%TRIP1=naqfdiag%TRIP1+RIP14R(NGL,N,K)
            ENDIF
            naqfdiag%TRGOM=naqfdiag%TRGOM+RCO2X(NGL,N,K)
            naqfdiag%TRGOC=naqfdiag%TRGOC+RCH4X(NGL,N,K)
            naqfdiag%TRGOD=naqfdiag%TRGOD+RGOMD(NGL,N,K)
            naqfdiag%TUPOX=naqfdiag%TUPOX+RUPOX(NGL,N,K)
            naqfdiag%TRDN3=naqfdiag%TRDN3+RDNO3(NGL,N,K)
            naqfdiag%TRDNB=naqfdiag%TRDNB+RDNOB(NGL,N,K)
            naqfdiag%TRDN2=naqfdiag%TRDN2+RDNO2(NGL,N,K)
            naqfdiag%TRD2B=naqfdiag%TRD2B+RDN2B(NGL,N,K)
            naqfdiag%TRDNO=naqfdiag%TRDNO+RDN2O(NGL,N,K)
            naqfdiag%TRGOH=naqfdiag%TRGOH+RH2GX(NGL,N,K)
          ENDDO

        ENDIF
640   CONTINUE
    ENDIF
650 CONTINUE
  DO 645 N=1,7
    IF(N.LE.3.OR.N.EQ.5)THEN
      IF(N.NE.3)THEN
        DO NGL=1,JG
          naqfdiag%TRGOA=naqfdiag%TRGOA+CGOMC(NGL,N,5)
        ENDDO
      ENDIF
    ENDIF
645 CONTINUE
!
!     ALLOCATE AGGREGATED TRANSFORMATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2O=net CO2 uptake
!     TRGOA=total CO2 uptake by autotrophs
!     TRGOM total CO2 emission by heterotrophs reducing O2
!     TRGOD=total CO2 emission by denitrifiers reducing NOx
!     RVOXA(3)=CH4 oxidation
!     RCH4O=net CH4 uptake
!     CGOMC=total CH4 uptake by autotrophs
!     TRGOC=total CH4 emission
!     RH2GO=net H2 uptake
!     RH2GZ,TRGOH=total H2 uptake, emission
!     RUPOXO,TUPOX=total O2 uptake
!     RN2G=total N2 production
!     TRDNO=total N2O reduction
!     RN2O=total N2O uptake
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RCN2O,RCN2B=nitrous acid reduction in non-band,band
!
  RCO2O(L,NY,NX)=naqfdiag%TRGOA-naqfdiag%TRGOM-naqfdiag%TRGOD
  RCH4O(L,NY,NX)=-naqfdiag%TRGOC
  DO NGL=1,JG
    RCO2O(L,NY,NX)=RCO2O(L,NY,NX)-RVOXA(NGL,3)
    RCH4O(L,NY,NX)=RCH4O(L,NY,NX)+RVOXA(NGL,3)+CGOMC(NGL,3,5)
  ENDDO
  RH2GO(L,NY,NX)=RH2GZ-naqfdiag%TRGOH
  RUPOXO(L,NY,NX)=naqfdiag%TUPOX
  RN2G(L,NY,NX)=-naqfdiag%TRDNO
  RN2O(L,NY,NX)=-naqfdiag%TRDN2-naqfdiag%TRD2B-RCN2O-RCN2B+naqfdiag%TRDNO
!
!     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate
!
  DO 655 K=0,4
    DO 660 M=1,jsken
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RCOSC(M,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RCOSN(M,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RCOSP(M,K)
660 CONTINUE
    DO 665 M=1,2
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RDORC(M,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RDORN(M,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RDORP(M,K)
665 CONTINUE
    XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RDOHC(K)
    XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RDOHN(K)
    XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RDOHP(K)
    XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)+RDOHA(K)
    DO 670 N=1,7
      DO NGL=1,JG
        XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CGOQC(NGL,N,K)
        XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-CGOMN(NGL,N,K)
        XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-CGOMP(NGL,N,K)
        XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CGOAC(NGL,N,K)+RCH3X(NGL,N,K)
      ENDDO
670 CONTINUE
    XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CSORP(K)
    XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-ZSORP(K)
    XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-PSORP(K)
    XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CSORPA(K)
655 CONTINUE
!
!     XNH4S,XNH4B=net change in NH4 in band,non-band
!     TRINH,TRINB=total NH4 mineraln-immobn in non-band,band
!     RVOXA(1),RVOXB(1)=total NH4 oxidation in non-band,band
!     XNO3S,XNO3B=net change in NO3 in band,non-band
!     TRINO,TRIOB=total NO3 immobn in non-band,band
!     RVOXA(2),RVOXB(2)=total NO2 oxidation in non-band,band
!     TRDN3,TRDNB=total NO3 reduction in non-band,band
!     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
!     XNO2S,XNO2B=net change in NO3 in band,non-band
!     TRDN2,TRD2B=total NO2 reduction in non-band,band
!     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
!     XH2PS,XH2BS=net change in H2PO4 in band,non-band
!     TRIPO,TRIPB=total H2PO4 mineraln-immobn in non-band,band
!     XH1PS,XH1BS=net change in HPO4 in band,non-band
!     TRIP1,TRIB1=total HPO4 mineraln-immobn in non-band,band
!     XN2GS=total N2 fixation
!     XZHYS=total H+ production
!     TRN2F=total N2 fixation
!
  XNH4S(L,NY,NX)=-naqfdiag%TRINH
  XNO3S(L,NY,NX)=-naqfdiag%TRINO-naqfdiag%TRDN3+RCNO3
  XNO2S(L,NY,NX)=+naqfdiag%TRDN3-naqfdiag%TRDN2-RCNO2
  XH2PS(L,NY,NX)=-naqfdiag%TRIPO
  XH1PS(L,NY,NX)=-naqfdiag%TRIP1
  XNH4B(L,NY,NX)=-naqfdiag%TRINB
  XNO3B(L,NY,NX)=-naqfdiag%TRIOB-naqfdiag%TRDNB+RCN3B
  XNO2B(L,NY,NX)=naqfdiag%TRDNB-naqfdiag%TRD2B-RCNOB
  DO NGL=1,JG
    XNH4S(L,NY,NX)=XNH4S(L,NY,NX)-RVOXA(NGL,1)
    XNO3S(L,NY,NX)=XNO3S(L,NY,NX)+RVOXA(NGL,2)
    XNO2S(L,NY,NX)=XNO2S(L,NY,NX)+RVOXA(NGL,1)-RVOXA(NGL,2)
    XNH4B(L,NY,NX)=XNH4B(L,NY,NX)-RVOXB(NGL,1)
    XNO3B(L,NY,NX)=XNO3B(L,NY,NX)+RVOXB(NGL,2)
    XNO2B(L,NY,NX)=XNO2B(L,NY,NX)+RVOXB(NGL,1)-RVOXB(NGL,2)
  ENDDO

  XH2BS(L,NY,NX)=-naqfdiag%TRIPB
  XH1BS(L,NY,NX)=-naqfdiag%TRIB1
  XN2GS(L,NY,NX)=naqfdiag%TRN2F
  TFNQ(L,NY,NX)=TFNX
  VOLQ(L,NY,NX)=VOLWZ
  end associate
  end subroutine AggregateTransformations
!------------------------------------------------------------------------------------------

  subroutine SubstrateCompetitionFactors(NGL,N,K,L,NY,NX,FOXYX,&
    FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQC,FOQA,&
    naqfdiag,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(out):: FOXYX,FNH4X
  real(r8),intent(out) :: FNB3X,FNB4X,FNO3X
  real(r8),intent(out) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8),intent(out) :: FOQC,FOQA
  type(NitroAQMFluxDiagType),INTENT(INOUT)::  naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
! begin_execution
  associate(                  &
    FOMA    => nmics%FOMA,    &
    FOMK    => nmics%FOMK,    &
    FNH4XR  =>nmicf%FNH4XR,   &
    FNO3XR  =>nmicf%FNO3XR,   &
    FP14XR  => nmicf%FP14XR,  &
    FPO4XR  => nmicf%FPO4XR   &
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  IF(ROXYY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FOXYX=AMAX1(FMN,ROXYS(NGL,N,K,L,NY,NX)/ROXYY(L,NY,NX))
  ELSE
    FOXYX=AMAX1(FMN,FOMA(NGL,N,K))
  ENDIF
  IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNH4X=AMAX1(FMN,RINHO(NGL,N,K,L,NY,NX)/RNH4Y(L,NY,NX))
  ELSE
    FNH4X=AMAX1(FMN,FOMA(NGL,N,K)*VLNH4(L,NY,NX))
  ENDIF
  IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB4X=AMAX1(FMN,RINHB(NGL,N,K,L,NY,NX)/RNHBY(L,NY,NX))
  ELSE
    FNB4X=AMAX1(FMN,FOMA(NGL,N,K)*VLNHB(L,NY,NX))
  ENDIF
  IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO3X=AMAX1(FMN,RINOO(NGL,N,K,L,NY,NX)/RNO3Y(L,NY,NX))
  ELSE
    FNO3X=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3(L,NY,NX))
  ENDIF
  IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB3X=AMAX1(FMN,RINOB(NGL,N,K,L,NY,NX)/RN3BY(L,NY,NX))
  ELSE
    FNB3X=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB(L,NY,NX))
  ENDIF
  IF(RPO4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FPO4X=AMAX1(FMN,RIPOO(NGL,N,K,L,NY,NX)/RPO4Y(L,NY,NX))
  ELSE
    FPO4X=AMAX1(FMN,FOMA(NGL,N,K)*VLPO4(L,NY,NX))
  ENDIF
  IF(RPOBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FPOBX=AMAX1(FMN,RIPBO(NGL,N,K,L,NY,NX)/RPOBY(L,NY,NX))
  ELSE
    FPOBX=AMAX1(FMN,FOMA(NGL,N,K)*VLPOB(L,NY,NX))
  ENDIF
  IF(RP14Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FP14X=AMAX1(FMN,RIPO1(NGL,N,K,L,NY,NX)/RP14Y(L,NY,NX))
  ELSE
    FP14X=AMAX1(FMN,FOMA(NGL,N,K)*VLPO4(L,NY,NX))
  ENDIF
  IF(RP1BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FP1BX=AMAX1(FMN,RIPB1(NGL,N,K,L,NY,NX)/RP1BY(L,NY,NX))
  ELSE
    FP1BX=AMAX1(FMN,FOMA(NGL,N,K)*VLPOB(L,NY,NX))
  ENDIF
  IF(K.LE.4)THEN
    IF(ROQCY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQC=AMAX1(FMN,ROQCS(NGL,N,K,L,NY,NX)/ROQCY(K,L,NY,NX))
    ELSE
      FOQC=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    naqfdiag%TFOQC=naqfdiag%TFOQC+FOQC
    IF(ROQAY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQA=AMAX1(FMN,ROQAS(NGL,N,K,L,NY,NX)/ROQAY(K,L,NY,NX))
    ELSE
      FOQA=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    naqfdiag%TFOQA=naqfdiag%TFOQA+FOQA
  ENDIF
  naqfdiag%TFOXYX=naqfdiag%TFOXYX+FOXYX
  naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4X
  naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3X
  naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4X
  naqfdiag%TFP14X=naqfdiag%TFP14X+FP14X
  naqfdiag%TFNH4B=naqfdiag%TFNH4B+FNB4X
  naqfdiag%TFNO3B=naqfdiag%TFNO3B+FNB3X
  naqfdiag%TFPO4B=naqfdiag%TFPO4B+FPOBX
  naqfdiag%TFP14B=naqfdiag%TFP14B+FP1BX
!
! FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
! MICROBIAL POPULATIONS IN SURFACE RESIDUE
! F*=fraction of substrate uptake relative to total uptake from
! previous hour in surface litter, labels as for soil layers above
!
  IF(L.EQ.0)THEN
    IF(RNH4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4XR(NGL,N,K)=AMAX1(FMN,RINHOR(NGL,N,K,NY,NX)/RNH4Y(NU(NY,NX),NY,NX))
    ELSE
      FNH4XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    IF(RNO3Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3XR(NGL,N,K)=AMAX1(FMN,RINOOR(NGL,N,K,NY,NX)/RNO3Y(NU(NY,NX),NY,NX))
    ELSE
      FNO3XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    IF(RPO4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4XR(NGL,N,K)=AMAX1(FMN,RIPOOR(NGL,N,K,NY,NX)/RPO4Y(NU(NY,NX),NY,NX))
    ELSE
      FPO4XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
    IF(RP14Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FP14XR(NGL,N,K)=AMAX1(FMN,RIPO1R(NGL,N,K,NY,NX)/RP14Y(NU(NY,NX),NY,NX))
    ELSE
      FP14XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
    ENDIF
  ENDIF
  IF(L.EQ.NU(NY,NX).AND.K.NE.3.AND.K.NE.4 &
    .AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
    naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4XR(NGL,N,K)
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3XR(NGL,N,K)
    naqfdiag%TFPO4X=naqfdiag%TFPO4X+FPO4XR(NGL,N,K)
    naqfdiag%TFP14X=naqfdiag%TFP14X+FP14XR(NGL,N,K)
  ENDIF
  end associate
  end subroutine SubstrateCompetitionFactors
!------------------------------------------------------------------------------------------

  subroutine NH3OxidizerCatabolism(NGL,N,K,L,NY,NX,XCO2,VOLWZ,TFNX,ECHZ,RGOMP,RVOXP,RVOXPA,&
    RVOXPB,naqfdiag,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  REAL(R8), INTENT(IN) :: VOLWZ,TFNX
  REAL(R8), INTENT(IN) :: XCO2
  real(r8),intent(out) :: ECHZ
  real(r8),intent(out) :: RGOMP,RVOXP
  real(r8),intent(out) :: RVOXPA,RVOXPB
  type(NitroAQMFluxDiagType),INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: FNH4S,FNHBS
  real(r8) :: FNH4,FNB4
  real(r8) :: FCN4S,FCN4B
  real(r8) :: RNNH4,RNNHB
  real(r8) :: VMXX,VMX4S
  real(r8) :: VMX4B
  real(r8) :: ZNFN4S,ZNFN4B
  real(r8) :: FSBST
  real(r8) :: VMXA

!     begin_execution
  associate(             &
    TFNG => nmics%TFNG,  &
    FCNP => nmics%FCNP,  &
    FOMA => nmics%FOMA,  &
    OMA  => nmics%OMA ,  &
    ROXYM=> nmicf%ROXYM, &
    ROXYP=> nmicf%ROXYP  &
  )
!
!     FACTOR TO REGULATE COMPETITION FOR NH4 AMONG DIFFERENT
!     MICROBIAL AND ROOT POPULATIONS FNH4
!
!     FNH4,FNB4=frac of total biol demand for NH4 in non-band, band
!
  FNH4S=VLNH4(L,NY,NX)
  FNHBS=VLNHB(L,NY,NX)
  IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNH4=AMAX1(FMN,RVMX4(NGL,N,K,L,NY,NX)/RNH4Y(L,NY,NX))
  ELSE
    FNH4=AMAX1(FMN,VLNH4(L,NY,NX)*FOMA(NGL,N,K))
  ENDIF
  IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB4=AMAX1(FMN,RVMB4(NGL,N,K,L,NY,NX)/RNHBY(L,NY,NX))
  ELSE
    FNB4=AMAX1(FMN,VLNHB(L,NY,NX)*FOMA(NGL,N,K))
  ENDIF
  naqfdiag%TFNH4X=naqfdiag%TFNH4X+FNH4
  naqfdiag%TFNH4B=naqfdiag%TFNH4B+FNB4
!
!     NITRIFICATION INHIBITION
!
!     ZNFN0=inhibition when fertilizer added
!     ZNFNI=reduction in inhibition since fertilizer added
!     CNH4S,CNH4B=NH4 concentrations in non-band, band
!     TFNX=temperature effect
!     RNFNI=rate constant for inhibition decline
!     ZHKI=inhibition from high CNH4
!     ZNFN4S,ZNFN4B=inhibition in non-band, band
!
  IF(ZNFN0(L,NY,NX).GT.ZEROS(NY,NX))THEN
    ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*(1.0-RNFNI*TFNX)
    ZNFN4S=ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX)/(1.0+CNH4S(L,NY,NX)/ZHKI)
    ZNFN4B=ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX)/(1.0+CNH4B(L,NY,NX)/ZHKI)
  ELSE
    ZNFN4S=1.0
    ZNFN4B=1.0
  ENDIF
!
!     NH3 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     NH3 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXX=potential NH3 oxidation, VMXH=specific oxidation
!     TFNG=temperature+water limitation, FCNP=N,P limitation
!     XCO2=aqueous CO2 limitation, OMA=active biomass
!     VMXA= non-substrate limited NH3 oxidation
!     VHKI=nonlinear increase in VMXA with VMXH
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     CNH4S,CNH4B=NH4 concentration in non-band, band
!     ZHKM=Km for NH4 uptake
!     FNH4,FNB4=fractions of total NH4 demand in non-band, band
!     ZNH4S,ZNH4B=NH4 amount in non-band, band
!     RNNH4,RNNHB=NH3 oxidation in non-band, band
!     RGOMP=O2-unlimited respiration
!     ECNH=efficiency CO2 conversion to biomass
!     RVMX4,RVMXB=nitrifier demand for NH4 in non-band, band
!
  ECHZ=EO2X
  VMXX=VMXH*TFNG(NGL,N,K)*FCNP(NGL,N,K)*XCO2*OMA(NGL,N,K)
  IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
    VMXA=VMXX/(1.0+VMXX/(VHKI*VOLWZ))
  ELSE
    VMXA=0.0_r8
  ENDIF
  FCN4S=FNH4S*CNH4S(L,NY,NX)/(CNH4S(L,NY,NX)+ZHKM)
  FCN4B=FNHBS*CNH4B(L,NY,NX)/(CNH4B(L,NY,NX)+ZHKM)
  FSBST=FCN4S+FCN4B
  VMX4S=VMXA*FCN4S
  VMX4B=VMXA*FCN4B
  RNNH4=AMAX1(0.0,AMIN1(VMX4S,FNH4*ZNH4S(L,NY,NX)))*ZNFN4S
  RNNHB=AMAX1(0.0,AMIN1(VMX4B,FNB4*ZNH4B(L,NY,NX)))*ZNFN4B
  RVOXP=RNNH4+RNNHB
  RVOXPA=RNNH4
  RVOXPB=RNNHB
  RGOMP=AMAX1(0.0,RVOXP*ECNH*ECHZ)
  RVMX4(NGL,N,K,L,NY,NX)=VMX4S
  RVMB4(NGL,N,K,L,NY,NX)=VMX4B
!
!     O2 DEMAND FROM NH3 OXIDATION
!
!     ROXYM=O2 demand from respiration by nitrifiers
!     ROXYP,ROXYM=O2 demand from respiration + NH3 oxidation
!
  ROXYM(NGL,N,K)=2.667*RGOMP
  ROXYP(NGL,N,K)=ROXYM(NGL,N,K)+3.429*RVOXP
  ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
!
  end associate
  end subroutine NH3OxidizerCatabolism
!------------------------------------------------------------------------------------------

  subroutine NO2OxidizerCatabolism(NGL,N,K,L,NY,NX,XCO2,ECHZ,RGOMP,RVOXP,&
    RVOXPA,RVOXPB,naqfdiag,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  REAL(r8), intent(in) :: XCO2
  real(r8), intent(out) :: ECHZ,RGOMP,RVOXP
  real(r8), intent(out) :: RVOXPA,RVOXPB
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: FNH4S,FNHBS
  REAL(R8) :: fno2,FNB2
  real(r8) :: FCN2S,FCN2B
  real(r8) :: RNNO2,RNNOB
  real(r8) :: VMX2S,VMX2B
  real(r8) :: FSBST
  real(r8) :: VMXA

!     begin_execution
  associate(               &
    TFNG => nmics%TFNG,    &
    FCNP => nmics%FCNP,    &
    FOMN => nmics%FOMN,    &
    OMA  => nmics%OMA ,    &
    ROXYM => nmicf%ROXYM,  &
    ROXYP => nmicf%ROXYP   &
  )
!     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
!     MICROBIAL POPULATIONS
!
!     FNO2=fraction of total biological demand for NO2 in non-band, band
!
  FNH4S=VLNH4(L,NY,NX)
  FNHBS=VLNHB(L,NY,NX)
  IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO2=AMAX1(FMN,RVMX2(NGL,N,K,L,NY,NX)/RNO2Y(L,NY,NX))
  ELSE
    FNO2=AMAX1(FMN,FOMN(NGL,N,K)*VLNO3(L,NY,NX))
  ENDIF
  IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB2=AMAX1(FMN,RVMB2(NGL,N,K,L,NY,NX)/RN2BY(L,NY,NX))
  ELSE
    FNB2=AMAX1(FMN,FOMN(NGL,N,K)*VLNOB(L,NY,NX))
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     NO2 CONCENTRATIONS
!
!     ECHZ=growth respiration efficiency
!     VMXA= non-substrate limited NH3 oxidation
!     VMXN=specific oxidation
!     TFNG=temperature+water limitation, FCNP=N,P limitation
!     XCO2=aqueous CO2 limitation, OMA=active biomass
!     OMA=active biomass
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     CNO2S,CNO2B=NO2 concentration in non-band, band
!     ZNKM=Km for NO2 uptake
!     FNO2,FNB2=fractions of total NO2 demand in non-band, band
!     ZNO2S,ZNO2B=NO2 amount in non-band, band
!     RNNO2,RNNOB=NO2 oxidation in non-band, band
!     RGOMP=O2-unlimited respiration
!     ECNO=efficiency CO2 conversion to biomass
!     RVMX2,RVMB2=nitrifier demand for NO2 in non-band, band
!
  ECHZ=EO2X
  VMXA=TFNG(NGL,N,K)*FCNP(NGL,N,K)*XCO2*OMA(NGL,N,K)*VMXN
  FCN2S=FNH4S*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+ZNKM)
  FCN2B=FNHBS*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+ZNKM)
  FSBST=FCN2S+FCN2B
  VMX2S=VMXA*FCN2S
  VMX2B=VMXA*FCN2B
  RNNO2=AMAX1(0.0,AMIN1(VMX2S,FNO2*ZNO2S(L,NY,NX)))
  RNNOB=AMAX1(0.0,AMIN1(VMX2B,FNB2*ZNO2B(L,NY,NX)))
  RVOXP=RNNO2+RNNOB
  RVOXPA=RNNO2
  RVOXPB=RNNOB
  RGOMP=AMAX1(0.0,RVOXP*ECNO*ECHZ)
  RVMX2(NGL,N,K,L,NY,NX)=VMX2S
  RVMB2(NGL,N,K,L,NY,NX)=VMX2B
!
!     O2 DEMAND FROM NO2 OXIDATION
!
!     ROXYM=O2 demand from respiration by nitrifiers
!     ROXYP,ROXYM=O2 demand from respiration + NO2 oxidation
!
  ROXYM(NGL,N,K)=2.667*RGOMP
  ROXYP(NGL,N,K)=ROXYM(NGL,N,K)+1.143*RVOXP
  ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)

  end associate
  end subroutine NO2OxidizerCatabolism
!------------------------------------------------------------------------------------------

  subroutine H2MethanogensCatabolism(NGL,N,K,L,NY,NX,ECHZ,RGOMP,XCO2,naqfdiag,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(out) :: ECHZ,RGOMP
  REAL(R8), INTENT(IN) :: XCO2
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: GH2X,GH2H
  real(r8) :: H2GSX
  real(r8) :: FSBST
  real(r8) :: VMXA

  associate(                  &
    TFNG => nmics%TFNG,       &
    FCNP => nmics%FCNP,       &
    OMA  => nmics%OMA ,       &
    ROXYM => nmicf%ROXYM,     &
    ROXYP=> nmicf%ROXYP       &
  )
!     begin_execution
!
!     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
!
!     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
!     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
!     VMXA=substrate-unlimited H2 oxidation rate
!     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (TRH2G)
!     CH2GS=H2 concentration, H2KM=Km for H2 uptake
!     RGOMP=H2 oxidation, ROXY*=O2 demand
!
!     and energy yield of hydrogenotrophic
!     methanogenesis GH2X at ambient H2 concentration CH2GS

  GH2X=8.3143E-03*TKS(L,NY,NX)*LOG((AMAX1(1.0E-03,CH2GS(L,NY,NX))/H2KI)**4)
  GH2H=GH2X/12.0
  ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0/(1.0+AMAX1(0.0,(GCOX+GH2H))/EOMH)))
  VMXA=TFNG(NGL,N,K)*FCNP(NGL,N,K)*XCO2*OMA(NGL,N,K)*VMXC
  H2GSX=H2GS(L,NY,NX)+0.111_r8*naqfdiag%TRH2G
  FSBST=CH2GS(L,NY,NX)/(CH2GS(L,NY,NX)+H2KM)
  RGOMP=AMAX1(0.0,AMIN1(1.5*H2GSX,VMXA*FSBST))
  ROXYM(NGL,N,K)=0.0_r8
  ROXYP(NGL,N,K)=0.0_r8
  ROXYS(NGL,N,K,L,NY,NX)=0.0_r8
  naqfdiag%TCH4A=naqfdiag%TCH4A+RGOMP

!
  end associate
  end subroutine H2MethanogensCatabolism
!------------------------------------------------------------------------------------------

  subroutine MethanotrophCatabolism(NGL,N,K,L,NY,NX,ECHZ,RGOMP,&
    RVOXP,RVOXPA,RVOXPB,naqfdiag,nmicf,nmics)

  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(out) :: ECHZ,RGOMP,RVOXP
  real(r8), intent(out) :: RVOXPA,RVOXPB
  type(NitroAQMFluxDiagType), intent(in) :: naqfdiag
  type(NitroMicFLuxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  integer  :: M,MM
  real(r8) :: CH4G1,CH4S1,CCH4S1
  real(r8) :: RCH4L1,RCH4F1,RCH4S1
  real(r8) :: RGOMP1,RCHDF
  real(r8) :: VMXA1,VOLWCH
  real(r8) :: FSBST
  real(r8) :: RVOXP1
  real(r8) :: VMXA
  REAL(R8) :: VOLWPM

  associate(             &
    TFNG => nmics%TFNG,  &
    FCNP => nmics%FCNP,  &
    OMA  => nmics%OMA ,  &
    ROXYM=> nmicf%ROXYM, &
    ROXYP=> nmicf%ROXYP  &
  )
!     begin_execution
!
!     CH4 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     CH4 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMXA=potential oxidation
!     TFNG=temperature+water effect,FCNP=N,P limitation
!     OMA=active biomass,VMX4=specific respiration rate
!     RCH4L=total aqueous CH4 exchange from previous hour
!     RCH4F=total gaseous CH4 exchange from previous hour
!     TCH4H+TCH4A=total CH4 generated from methanogenesis
!     XNPG=1.0/(NPH*NPT)
!     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
!     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
!     VOLPM,VOLWM=air,water-filled porosity
!     SCH4L=CH4 aqueous solubility
!     CCK4=Km for CH4 uptake
!     ECHO=efficiency CO2 conversion to biomass
!     RGOMP1=substrate-limited CH4 oxidation
!     RCHDF=gaseous-aqueous CH4 exchange
!     DFGS=rate constant for gaseous-aqueous exchange
!
  ECHZ=EH4X
  VMXA=TFNG(NGL,N,K)*FCNP(NGL,N,K)*OMA(NGL,N,K)*VMX4
  RCH4L1=RCH4L(L,NY,NX)*XNPG
  RCH4F1=RCH4F(L,NY,NX)*XNPG
  RCH4S1=(naqfdiag%TCH4H+naqfdiag%TCH4A)*XNPG
  IF(L.EQ.0)THEN
    CH4G1=CCH4E(NY,NX)*VOLPM(1,L,NY,NX)
  ELSE
    CH4G1=CCH4G(L,NY,NX)*VOLPM(1,L,NY,NX)
  ENDIF
  CH4S1=CH4S(L,NY,NX)
  VMXA1=VMXA*XNPG
  RVOXP=0.0_r8
  RGOMP=0.0_r8
!
!     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
!     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
!
  DO 320 M=1,NPH
    IF(VOLWM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLWCH=VOLWM(M,L,NY,NX)*SCH4L(L,NY,NX)
      VOLWPM=VOLWCH+VOLPM(M,L,NY,NX)
      DO 325 MM=1,NPT
        CH4G1=CH4G1+RCH4F1
        CH4S1=CH4S1+RCH4L1+RCH4S1
        CCH4S1=AMAX1(0.0,safe_adb(CH4S1,VOLWM(M,L,NY,NX)))
        FSBST=CCH4S1/(CCH4S1+CCK4)
        RVOXP1=AMIN1(AMAX1(0.0,CH4S1)/(1.0+ECHO*ECHZ),VMXA1*FSBST)
        RGOMP1=RVOXP1*ECHO*ECHZ
        CH4S1=CH4S1-RVOXP1-RGOMP1
        IF(THETPM(M,L,NY,NX).GT.THETX)THEN
          RCHDF=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),CH4G1)*VOLWCH &
            -CH4S1*VOLPM(M,L,NY,NX))/VOLWPM
        ELSE
          RCHDF=0.0_r8
        ENDIF
        CH4G1=CH4G1-RCHDF
        CH4S1=CH4S1+RCHDF
        RVOXP=RVOXP+RVOXP1
        RGOMP=RGOMP+RGOMP1

325   CONTINUE
    ENDIF
320   CONTINUE
  RVOXPA=RVOXP
  RVOXPB=0.0_r8
!
!     O2 DEMAND FROM CH4 OXIDATION
!
!     ROXYM=O2 demand from respiration
!     ROXYP=O2 demand from respiration + CH4 oxidation
!
  ROXYM(NGL,N,K)=2.667*RGOMP
  ROXYP(NGL,N,K)=ROXYM(NGL,N,K)+4.00*RVOXP
  ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
  end associate
  end subroutine MethanotrophCatabolism
!------------------------------------------------------------------------------------------

  subroutine AcetoMethanogenCatabolism(NGL,N,K,L,NY,NX,TFNX,WFNG,FOQA,ECHZ,&
    FGOCP,FGOAP,RGOMP,naqfdiag,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in) :: FOQA,TFNX,WFNG
  real(r8), intent(out) :: ECHZ
  REAL(R8), intent(out) :: FGOCP,FGOAP,RGOMP
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  real(r8) :: FSBST
  real(r8) :: GOMX,GOMM
  real(r8) :: RGOGY,RGOGZ,RGOGX

! begin_execution
  associate(                  &
    FCNP  => nmics%FCNP,      &
    OMA   => nmics%OMA ,      &
    ROXYP => nmicf%ROXYP,     &
    ROXYM => nmicf%ROXYM,     &
    ROQCD => nmicf%ROQCD      &
  )
  GOMX=8.3143E-03*TKS(L,NY,NX)*LOG((AMAX1(ZERO,COQA(K,L,NY,NX))/OAKI))
  GOMM=GOMX/24.0
  ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0/(1.0+AMAX1(0.0,(GC4X+GOMM))/EOMH)))
!
!     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     COQA=DOA concentration
!     OQKAM=Km for acetate uptake,FCNP=N,P limitation
!     VMXM=specific respiration rate
!     WFNG=water stress effect, OMA=active biomass
!     TFNX=temp stress effect, FOQA= acetate limitation
!     RGOGX=substrate-limited respiration of acetate
!     RGOGX=competition-limited respiration of acetate
!     OQA=acetate, FOQA=fraction of biological demand for acetate
!     RGOMP=O2-unlimited respiration of acetate
!     ROXY*=O2 demand, ROQCS,ROQCA=DOC, acetate demand
!     ROQCD=microbial respiration used to represent microbial activity
!
  FSBST=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKAM)
  RGOGY=AMAX1(0.0,FCNP(NGL,N,K)*VMXM*WFNG*OMA(NGL,N,K))
  RGOGZ=RGOGY*FSBST*TFNX
  RGOGX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*ECHZ)
  RGOMP=AMIN1(RGOGX,RGOGZ)
  FGOCP=0.0_r8
  FGOAP=1.0
  ROXYM(NGL,N,K)=0.0_r8
  ROXYP(NGL,N,K)=0.0_r8
  ROXYS(NGL,N,K,L,NY,NX)=0.0_r8
  ROQCS(NGL,N,K,L,NY,NX)=0.0_r8
  ROQAS(NGL,N,K,L,NY,NX)=RGOGZ
  ROQCD(NGL,N,K)=0.0_r8
  naqfdiag%TCH4H=naqfdiag%TCH4H+0.5*RGOMP
  end associate
  end subroutine AcetoMethanogenCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobicHeterotrophCatabolism(NGL,N,K,L,NY,NX,TFNX,WFNG,FOQC,FOQA,&
    ECHZ,FGOCP,FGOAP,RGOCP,RGOMP,nmicf,nmics,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  REAL(R8), INTENT(IN) :: FOQC,FOQA,WFNG,TFNX
  real(r8), intent(out) :: ECHZ,RGOMP
  REAL(R8), INTENT(OUT) :: FGOCP,FGOAP
  real(r8), intent(out) :: RGOCP
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroOMcplxStateType),intent(inout):: ncplxs
  real(r8) :: EO2Q
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX
  real(r8) :: RGOAP
  real(r8) :: ROXYSX,ROQCSX,ROQASX
  real(r8) :: FSBST
!     begin_execution
  associate(                  &
    OMA   => nmics%OMA  ,     &
    FCNP  => nmics%FCNP ,     &
    ROXYM => nmicf%ROXYM,     &
    ROXYP=> nmicf%ROXYP ,     &
    ROQCD => nmicf%ROQCD,     &
    FOCA  => ncplxs%FOCA,     &
    FOAA  => ncplxs%FOAA      &
  )
!     ENERGY YIELDS OF O2 REDOX REACTIONS
!     E* = growth respiration efficiency calculated in PARAMETERS
!
  IF(N.EQ.1)THEN
    EO2Q=EO2X
  ELSEIF(N.EQ.2)THEN
    EO2Q=EO2D
  ELSEIF(N.EQ.3)THEN
    EO2Q=EO2G
  ELSEIF(N.EQ.6)THEN
    EO2Q=ENFX
  ENDIF
!
! O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
! 'RGO*Z'FROM SPECIFIC RESPIRATION RATE, ACTIVE BIOMASS, DOC OR
! ACETATE CONCENTRATION,MICROBIAL C:N:P FACTOR, AND TEMPERATURE
! FOLLOWED BY POTENTIAL RESPIRATION RATES 'RGO*P' WITH UNLIMITED
! SUBSTRATE USED FOR MICROBIAL COMPETITION FACTOR

! COQC,COQA=DOC,DOA concentration, FOCA,FOAA=DOC,DOA vs DOC+DOA
! FCNP=N,P limitation,VMXO=specific respiration rate
! WFNG=water stress effect, OMA=active biomass
! TFNX=temp stress effect,FOQC,FOQA=OQC,OQA limitation
! RGOMP=O2-unlimited respiration of DOC+DOA
! RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
!
  FSBSTC=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)
  FSBSTA=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKA)
  FSBST=FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
  RGOCY=AMAX1(0.0,FCNP(NGL,N,K)*VMXO*WFNG*OMA(NGL,N,K))
  RGOCZ=RGOCY*FSBSTC*FOCA(K)*TFNX
  RGOAZ=RGOCY*FSBSTA*FOAA(K)*TFNX
  RGOCX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*EO2Q)
  RGOAX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*EO2A)
  RGOCP=AMIN1(RGOCX,RGOCZ)
  RGOAP=AMIN1(RGOAX,RGOAZ)
  RGOMP=RGOCP+RGOAP
  IF(RGOMP.GT.ZEROS(NY,NX))THEN
    FGOCP=RGOCP/RGOMP
    FGOAP=RGOAP/RGOMP
  ELSE
    FGOCP=1.0
    FGOAP=0.0_r8
  ENDIF
!
! ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
! BY HETEROTROPHIC AEROBES

! ECHZ=growth respiration yield
! ROXYM,ROXYP,ROXYS=O2 demand from DOC,DOA oxidation
! ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
! ROQCD=microbial respiration used to represent microbial activity
!
  ECHZ=EO2Q*FGOCP+EO2A*FGOAP
  ROXYM(NGL,N,K)=2.667*RGOMP
  ROXYP(NGL,N,K)=ROXYM(NGL,N,K)
  ROXYSX=ROXYS(NGL,N,K,L,NY,NX)
  ROQCSX=ROQCS(NGL,N,K,L,NY,NX)
  ROQASX=ROQAS(NGL,N,K,L,NY,NX)
  ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
  ROQCS(NGL,N,K,L,NY,NX)=RGOCZ
  ROQAS(NGL,N,K,L,NY,NX)=RGOAZ
  ROQCD(NGL,N,K)=RGOCY
!
  end associate
  end subroutine AerobicHeterotrophCatabolism
!------------------------------------------------------------------------------------------

  subroutine AnaerobCatabolism(NGL,N,K,L,NY,NX,TFNX,WFNG,FOQC,ECHZ,FGOCP,FGOAP,RGOMP,&
    naqfdiag,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  REAL(R8), INTENT(IN) :: FOQC,WFNG,TFNX
  real(r8), intent(out) :: ECHZ,RGOMP
  real(r8), intent(out) :: FGOCP,FGOAP
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType), intent(inout):: nmics
  real(r8) :: GH2X,GH2F
  real(r8) :: GOAX,GOAF
  real(r8) :: GHAX
  REAL(R8) :: oxyi
  real(r8) :: RGOFX,RGOFY,RGOFZ
  real(r8) :: FSBST
!     begin_execution
  associate(                    &
    FCNP   => nmics%FCNP,       &
    OMA    => nmics%OMA ,       &
    ROXYM  => nmicf%ROXYM,      &
    ROXYP  => nmicf%ROXYP  ,    &
    ROQCD  => nmicf%ROQCD       &
  )
  GH2X=8.3143E-03*TKS(L,NY,NX)*LOG((AMAX1(1.0E-03,CH2GS(L,NY,NX))/H2KI)**4)
  GH2F=GH2X/72.0
  GOAX=8.3143E-03*TKS(L,NY,NX)*LOG((AMAX1(ZERO,COQA(K,L,NY,NX))/OAKI)**2)
  GOAF=GOAX/72.0
  GHAX=GH2F+GOAF
  IF(N.EQ.4)THEN
    ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0/(1.0+AMAX1(0.0,(GCHX-GHAX))/EOMF)))
  ELSE
    ECHZ=AMAX1(ENFX,AMIN1(1.0,1.0/(1.0+AMAX1(0.0,(GCHX-GHAX))/EOMN)))
  ENDIF
!
!     RESPIRATION RATES BY HETEROTROPHIC ANAEROBES 'RGOMP' FROM
!     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
!     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
!     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
!     MICROBIAL COMPETITION FACTOR
!
!     OXYI=O2 inhibition of fermentation
!     FCNP=N,P limitation on respiration
!     VMXF=maximum respiration rate by fermenters
!     WFNG=water stress effect on respiration
!     OMA=active fermenter biomass
!     TFNX=temp stress effect, FOQC=OQC limitation
!     RFOMP=O2-unlimited respiration of DOC
!     ROQCD=microbial respiration used to represent microbial activity
!
  OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS(L,NY,NX)+2.5)))
  FSBST=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)*OXYI
  RGOFY=AMAX1(0.0,FCNP(NGL,N,K)*VMXF*WFNG*OMA(NGL,N,K))
  RGOFZ=RGOFY*FSBST*TFNX
  RGOFX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*ECHZ)
  RGOMP=AMIN1(RGOFX,RGOFZ)
  FGOCP=1.0
  FGOAP=0.0_r8
  ROXYM(NGL,N,K)=0.0_r8
  ROXYP(NGL,N,K)=0.0_r8
  ROXYS(NGL,N,K,L,NY,NX)=0.0_r8
  ROQCS(NGL,N,K,L,NY,NX)=RGOFZ
  ROQAS(NGL,N,K,L,NY,NX)=0.0_r8
  ROQCD(NGL,N,K)=RGOFY
  naqfdiag%TRH2G=naqfdiag%TRH2G+RGOMP
!
  end associate
  end subroutine AnaerobCatabolism
!------------------------------------------------------------------------------------------

  subroutine HeteroDenitrificCatabolism(NGL,N,K,L,NY,NX,FOQC,RGOCP,&
    VOLWZ,naqfdiag,nmicf,nmics,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  REAL(R8), INTENT(IN) :: FOQC,RGOCP
  real(r8), intent(in) :: VOLWZ
  type(NitroAQMFluxDiagType), intent(inout) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxStateType),intent(inout) :: ncplxs
  real(r8) :: FNO3S,FNO3B
  real(r8) :: FNO2S,FNO2B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: FNO3,FNB3
  real(r8) :: FVMXDX
  REAL(R8) :: FN2O
  real(r8) :: OQCZ3
  real(r8) :: OQCD3
  real(r8) :: OQCD3S
  real(r8) :: OQCD3B
  real(r8) :: OQCZ2
  real(r8) :: OQCD2
  real(r8) :: OQCD2S
  real(r8) :: OQCD2B
  real(r8) :: OQCZ1
  real(r8) :: OQCD1
  real(r8) :: ROXYD
  real(r8) :: RDNO3X,RDNOBX,RDNOX
  real(r8) :: RDNOT
  real(r8) :: RGOM3X,RGOMD3
  real(r8) :: RDNO2X,RDN2X,RDN2T,RGOM2X,RGOMD2,RDN2OX,RGOM1X
  real(r8) :: RGOMD1
  real(r8) :: VMXD3
  real(r8) :: VMXDXS
  real(r8) :: VMXDXB
  real(r8) :: VMXDXT
  real(r8) :: VMXD3S,VMXD3B,VMXD2,VMXD2S,VMXD2B,VMXD1
  real(r8) :: VMXD1S
  real(r8) :: ZNO3SX,ZNO3BX
  real(r8) :: ZNO2SX,ZNO2BX
  real(r8) :: Z2OSX

! begin_execution
  associate(                  &
    WFN    => nmics%WFN,      &
    FOMA   => nmics%FOMA,     &
    ROXYM => nmicf%ROXYM,     &
    ROXYO => nmicf%ROXYO,     &
    RDNO3 => nmicf%RDNO3,     &
    RDNOB => nmicf%RDNOB,     &
    RDNO2 => nmicf%RDNO2,     &
    RDN2B => nmicf%RDN2B,     &
    RDN2O => nmicf%RDN2O,     &
    RGOMD  => nmicf%RGOMD,    &
    RGOMY  => nmicf%RGOMY,    &
    OSRH  => ncplxs%OSRH      &
  )
!
! FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
! AND ROOT POPULATIONS
!
! FNO3,FNB3=fraction of total biological demand for NO3
!

  FNO3S=VLNO3(L,NY,NX)
  FNO3B=VLNOB(L,NY,NX)
  IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO3=AMAX1(FMN,RVMX3(NGL,N,K,L,NY,NX)/RNO3Y(L,NY,NX))
  ELSE
    FNO3=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3(L,NY,NX))
  ENDIF
  IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB3=AMAX1(FMN,RVMB3(NGL,N,K,L,NY,NX)/RN3BY(L,NY,NX))
  ELSE
    FNB3=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB(L,NY,NX))
  ENDIF
  naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3
  naqfdiag%TFNO3B=naqfdiag%TFNO3B+FNB3
!
!     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
!
!     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
!     VMXD3=demand for NO3-N reduction
!     VMXDXS,VMXDXB=maximum NO3 reduction in non-band, band
!     FNO3S,FNO3B=fractions of total NO3 in non-band, band
!     CNO3S,CNO3B=NO3 concentrations in non-band, band
!     Z3KM,Z2KM=Km for NO3, NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD3S,VMXD3B=substrate-unlimited NO3 reduction in non-band,band
!     OQCD3S,OQCD3B=DOC limitation to NO3 reduction in non-band, band
!     RDNO3,RDNOB=substrate-limited NO3 reduction in non-band,band
!     RGOM3X,RGOMD3=substrate-unltd,-ltd respn from NO3 reduction
!     RVMX3,RVMB3=demand for NO3 reduction in non-band,band
!
  ROXYD=AMAX1(0.0,ROXYM(NGL,N,K)-ROXYO(NGL,N,K))
  VMXD3=0.875*ROXYD
  IF(CNO3S(L,NY,NX).GT.ZERO)THEN
    VMXDXS=FNO3S*VMXD3*CNO3S(L,NY,NX)/(CNO3S(L,NY,NX)+Z3KM) &
      /(1.0+(CNO2S(L,NY,NX)*Z3KM)/(CNO3S(L,NY,NX)*Z2KM))
  ELSE
    VMXDXS=0.0_r8
  ENDIF
  IF(CNO3B(L,NY,NX).GT.ZERO)THEN
    VMXDXB=FNO3B*VMXD3*CNO3B(L,NY,NX)/(CNO3B(L,NY,NX)+Z3KM) &
      /(1.0+(CNO2B(L,NY,NX)*Z3KM)/(CNO3B(L,NY,NX)*Z2KM))
  ELSE
    VMXDXB=0.0_r8
  ENDIF
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
    FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD3S=VMXDXS*FVMXDX
  VMXD3B=VMXDXB*FVMXDX
  OQCZ3=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC-RGOCP*WFN(NGL,N,K))
  OQCD3=OQCZ3/ECN3
  OQCD3S=OQCD3*FNO3S
  OQCD3B=OQCD3*FNO3B
  ZNO3SX=ZNO3S(L,NY,NX)*FNO3
  ZNO3BX=ZNO3B(L,NY,NX)*FNB3
  RDNO3X=AMAX1(0.0,AMIN1(ZNO3SX,VMXD3S))
  RDNOBX=AMAX1(0.0,AMIN1(ZNO3BX,VMXD3B))
  RDNO3(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD3S,OQCD3S,ZNO3SX))
  RDNOB(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD3B,OQCD3B,ZNO3BX))
  RDNOX=RDNO3X+RDNOBX
  RDNOT=RDNO3(NGL,N,K)+RDNOB(NGL,N,K)
  RGOM3X=ECN3*RDNOX
  RGOMD3=ECN3*RDNOT
  RVMX3(NGL,N,K,L,NY,NX)=VMXD3S
  RVMB3(NGL,N,K,L,NY,NX)=VMXD3B
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
  FNO2S=VLNO3(L,NY,NX)
  FNO2B=VLNOB(L,NY,NX)
  IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO2=AMAX1(FMN,RVMX2(NGL,N,K,L,NY,NX)/RNO2Y(L,NY,NX))
  ELSE
    FNO2=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3(L,NY,NX))
  ENDIF
  IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB2=AMAX1(FMN,RVMB2(NGL,N,K,L,NY,NX)/RN2BY(L,NY,NX))
  ELSE
    FNB2=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB(L,NY,NX))
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2 AND NO3 IN BAND AND NON-BAND SOIL ZONES
!
!     VMXD2=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM,Z1KM=Km for NO2, N2O uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD2S,VMXD2B=substrate-unlimited NO2 reduction in non-band,band
!     OQCD2S,OQCD2B=DOC limitation to NO2 reduction in non-band, band
!     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
!     RGOM2X,RGOMD2=substrate-unltd,-ltd respn from NO2 reduction
!
  VMXD2=VMXD3-RDNOT
  IF(CNO2S(L,NY,NX).GT.ZERO)THEN
    VMXDXS=FNO2S*VMXD2*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM) &
      /(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2S(L,NY,NX)*Z1KM))
  ELSE
    VMXDXS=0.0_r8
  ENDIF
  IF(CNO2B(L,NY,NX).GT.ZERO)THEN
    VMXDXB=FNO2B*VMXD2*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM) &
      /(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2B(L,NY,NX)*Z1KM))
  ELSE
    VMXDXB=0.0_r8
  ENDIF
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
    FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD2S=VMXDXS*FVMXDX
  VMXD2B=VMXDXB*FVMXDX
  OQCZ2=AMAX1(0.0,OQCZ3-RGOMD3)
  OQCD2=OQCZ2/ECN2
  OQCD2S=OQCD2*FNO3S
  OQCD2B=OQCD2*FNO3B
  ZNO2SX=(ZNO2S(L,NY,NX)+RDNO3(NGL,N,K))*FNO2
  ZNO2BX=(ZNO2B(L,NY,NX)+RDNOB(NGL,N,K))*FNB2
  RDNO2X=AMAX1(0.0,AMIN1(ZNO2SX,VMXD2S))
  RDNOBX=AMAX1(0.0,AMIN1(ZNO2BX,VMXD2B))
  RDNO2(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD2S,OQCD2S,ZNO2SX))
  RDN2B(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD2B,OQCD2B,ZNO2BX))
  RDN2X=RDNO2X+RDNOBX
  RDN2T=RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
  RGOM2X=ECN2*RDN2X
  RGOMD2=ECN2*RDN2T
  RVMX2(NGL,N,K,L,NY,NX)=VMXD2S
  RVMB2(NGL,N,K,L,NY,NX)=VMXD2B
!
!     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
!     AND ROOT POPULATIONS
!
!     FN2O=fraction of total biological demand for N2O
!
  IF(RN2OY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FN2O=AMAX1(FMN,RVMX1(NGL,N,K,L,NY,NX)/RN2OY(L,NY,NX))
  ELSE
    FN2O=AMAX1(FMN,FOMA(NGL,N,K))
  ENDIF
  naqfdiag%TFN2OX=naqfdiag%TFN2OX+FN2O
!
!     N2O REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS N2O
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2, NO3 AND NO2 IN BAND AND NON-BAND SOIL ZONES
!
!     VMXD1=demand for N2O-N reduction
!     VMXDXS=maximum N2O reduction
!     CZ2OS=N2O concentrations
!     Z1KM=Km for N2O uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD1S=substrate-unlimited N2O reduction
!     OQCD1=DOC limitation to N2O reduction
!     RDN2O=substrate-limited N2O reduction
!     RGOM1X,RGOMD1=substrate-unltd,-ltd  respn from N2O reduction
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NOx reduction
!     RVMX1=demand for N2O reduction
!
  VMXD1=(VMXD2-RDN2T)*2.0
  VMXDXS=VMXD1*CZ2OS(L,NY,NX)/(CZ2OS(L,NY,NX)+Z1KM)
  IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
    FVMXDX=1.0/(1.0+VMXDXS/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD1S=VMXDXS*FVMXDX
  OQCZ1=AMAX1(0.0,OQCZ2-RGOMD2)
  OQCD1=OQCZ1/ECN1
  Z2OSX=(Z2OS(L,NY,NX)+RDN2T)*FN2O
  RDN2OX=AMAX1(0.0,AMIN1(Z2OSX,VMXD1S))
  RDN2O(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD1S,OQCD1,Z2OSX))
  RGOM1X=ECN1*RDN2OX
  RGOMD1=ECN1*RDN2O(NGL,N,K)
  RGOMY(NGL,N,K)=RGOM3X+RGOM2X+RGOM1X
  RGOMD(NGL,N,K)=RGOMD3+RGOMD2+RGOMD1
  RVMX1(NGL,N,K,L,NY,NX)=VMXD1S
  end associate
  end subroutine HeteroDenitrificCatabolism
!------------------------------------------------------------------------------------------

  subroutine AutotrophDenitrificCatabolism(NGL,N,K,L,NY,NX,XCO2,VOLWZ,naqfdiag,nmicf,nmics)

  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in) :: VOLWZ,XCO2
  type(NitroAQMFluxDiagType), INTENT(INOUT) :: naqfdiag
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  real(r8) :: FNO2S,FNO2B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: FVMXDX
  real(r8) :: ROXYD,RDNOT
  real(r8) :: VMXDXS
  real(r8) :: VMXDXB
  real(r8) :: VMXDXT
  real(r8) :: VMXD4
  real(r8) :: VMXD4S
  real(r8) :: VMXD4B
  real(r8) :: ZNO2SX,ZNO2BX

!     begin_execution
  associate(                  &
    FOMN => nmics%FOMN,       &
    ROXYM  => nmicf%ROXYM,    &
    ROXYO => nmicf%ROXYO ,    &
    RDNO3 => nmicf%RDNO3 ,    &
    RDNOB => nmicf%RDNOB ,    &
    RDNO2 => nmicf%RDNO2 ,    &
    RDN2B => nmicf%RDN2B ,    &
    RDN2O => nmicf%RDN2O ,    &
    RGOMD  => nmicf%RGOMD,    &
    RGOMY  => nmicf%RGOMY,    &
    RVOXA  => nmicf%RVOXA,    &
    RVOXB  => nmicf%RVOXB     &
  )
!
!     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
!     POPULATIONS
!
!     FNO2,FNB2=fraction of total biological demand for NO2
!
  IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNO2=AMAX1(FMN,RVMX2(NGL,N,K,L,NY,NX)/RNO2Y(L,NY,NX))
  ELSE
    FNO2=AMAX1(FMN,FOMN(NGL,N,K)*VLNO3(L,NY,NX))
  ENDIF
  IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
    FNB2=AMAX1(FMN,RVMB2(NGL,N,K,L,NY,NX)/RN2BY(L,NY,NX))
  ELSE
    FNB2=AMAX1(FMN,FOMN(NGL,N,K)*VLNOB(L,NY,NX))
  ENDIF
  naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
!
!     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
!     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
!     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
!     NOT ACCEPTED BY O2
!
!     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
!     VMXD4=demand for NO2-N reduction
!     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
!     FNO2S,FNO2B=fractions of total NO2 in non-band, band
!     CNO2S,CNO2B=NO2 concentrations in non-band, band
!     Z2KM=Km for NO2 uptake
!     FVMXDX=nonlinear effect of product inhibition for NOx reduction
!     VMKI=product inhibition for NOx reduction
!     VMXD4S,VMXD4B=substrate-unlimited NO2 reduction in non-band,band
!     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
!     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NO2 reduction
!     ECNO=efficiency CO2 conversion to biomass
!     ECHZ=growth respiration efficiency
!     RVOXA,RVOXB=total O2-limited (1)NH4,(2)NO2,(3)CH4 oxidation
!
  FNO2S=VLNO3(L,NY,NX)
  FNO2B=VLNOB(L,NY,NX)
  ROXYD=AMAX1(0.0,ROXYM(NGL,N,K)-ROXYO(NGL,N,K))
  VMXD4=0.875*ROXYD*XCO2
  VMXDXS=FNO2S*VMXD4*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM)
  VMXDXB=FNO2B*VMXD4*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM)
  VMXDXT=VMXDXS+VMXDXB
  IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
    FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ))
  ELSE
    FVMXDX=0.0_r8
  ENDIF
  VMXD4S=VMXDXS*FVMXDX
  VMXD4B=VMXDXB*FVMXDX
  ZNO2SX=ZNO2S(L,NY,NX)+RVOXA(NGL,1)
  ZNO2BX=ZNO2B(L,NY,NX)+RVOXB(NGL,1)
  RDNO2(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD4S,ZNO2SX))
  RDN2B(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD4B,ZNO2BX))
  RDNOT=RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
  RGOMY(NGL,N,K)=0.0_r8
  RGOMD(NGL,N,K)=RDNOT*ECNO*ENOX
  RDNO3(NGL,N,K)=0.0_r8
  RDNOB(NGL,N,K)=0.0_r8
  RDN2O(NGL,N,K)=0.0_r8
  RVMX2(NGL,N,K,L,NY,NX)=VMXD4S
  RVMB2(NGL,N,K,L,NY,NX)=VMXD4B
  RVOXA(NGL,N)=RVOXA(NGL,N)+0.333*RDNO2(NGL,N,K)
  RVOXB(NGL,N)=RVOXB(NGL,N)+0.333*RDN2B(NGL,N,K)
!     TRN2ON(NY,NX)=TRN2ON(NY,NX)+RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
  end associate
  end subroutine AutotrophDenitrificCatabolism
!------------------------------------------------------------------------------------------

  subroutine AerobsO2Uptake(NGL,N,K,L,NY,NX,FOXYX,OXKX,RGOMP,RVOXP,RVOXPA,RVOXPB,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in) :: OXKX,FOXYX,RGOMP,RVOXP
  real(r8), intent(in) :: RVOXPA
  real(r8), intent(in) :: RVOXPB
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType),intent(inout) :: nmics
  integer  :: M,MX
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,OLSGL1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMX,ROXYFX
  real(r8) :: ROXYLX
  real(r8) :: RRADO,RMPOX,ROXDFQ
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X
  real(r8) :: VOLWPM

  ! begin_execution
  associate(                  &
    WFN    => nmics%WFN,      &
    OMA    => nmics%OMA,      &
    RUPOX  => nmicf%RUPOX,    &
    RGOMO  => nmicf%RGOMO,    &
    ROXYM => nmicf%ROXYM,     &
    ROXYP=> nmicf%ROXYP ,     &
    ROXYO => nmicf%ROXYO,     &
    RH2GX  => nmicf%RH2GX,    &
    ROQCD  => nmicf%ROQCD,    &
    RCO2X  => nmicf%RCO2X,    &
    RCH3X  => nmicf%RCH3X,    &
    RCH4X  => nmicf%RCH4X,    &
    RVOXA  => nmicf%RVOXA,    &
    RVOXB  => nmicf%RVOXB     &
  )
  IF(ROXYP(NGL,N,K).GT.ZEROS(NY,NX).AND.FOXYX.GT.ZERO)THEN
    IF(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMX=ROXYP(NGL,N,K)*XNPG
      ROXYFX=ROXYF(L,NY,NX)*XNPG*FOXYX
      OLSGL1=OLSGL(L,NY,NX)*XNPG
      IF(L.NE.0)THEN
        OXYG1=OXYG(L,NY,NX)*FOXYX
        ROXYLX=ROXYL(L,NY,NX)*XNPG*FOXYX
      ELSE
        OXYG1=COXYG(L,NY,NX)*VOLPM(1,L,NY,NX)*FOXYX
        ROXYLX=(ROXYL(L,NY,NX)+FLQRQ(NY,NX)*COXR(NY,NX) &
          +FLQRI(NY,NX)*COXQ(NY,NX))*XNPG*FOXYX
      ENDIF
      OXYS1=OXYS(L,NY,NX)*FOXYX
!
      !write(*,*)'O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP'
!     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
!
      DO 420 M=1,NPH
        !
        !     ACTUAL REDUCTION OF AQUEOUS BY AEROBES CALCULATED
        !     FROM MASS FLOW PLUS DIFFUSION = ACTIVE UPTAKE
        !     COUPLED WITH DISSOLUTION OF GASEOUS O2 DURING REDUCTION
        !     OF AQUEOUS O2 FROM DISSOLUTION RATE CONSTANT 'DFGS'
        !     CALCULATED IN 'WATSUB'
        !
        !     VOLWM,VOLPM,VOLX=water, air and total volumes
        !     ORAD=microbial radius,FILM=water film thickness
        !     DIFOX=aqueous O2 diffusion, TORT=tortuosity
        !     BIOS=microbial number, OMA=active biomass
        !     SOXYL=O2 solubility, OXKX=Km for O2 uptake
        !     OXYS,COXYS=aqueous O2 amount, concentration
        !     OXYG,COXYG=gaseous O2 amount, concentration
        !     RMPOX,ROXSK=O2 uptake
        !
        !write(*,*)'VOLY(L,NY,NX)=',VOLY(L,NY,NX)
        THETW1=AMAX1(0.0_r8,safe_adb(VOLWM(M,L,NY,NX),VOLY(L,NY,NX)))
        RRADO=ORAD*(FILM(M,L,NY,NX)+ORAD)/FILM(M,L,NY,NX)
        DIFOX=TORT(M,L,NY,NX)*OLSGL1*12.57_r8*BIOS*OMA(NGL,N,K)*RRADO
        VOLWOX=VOLWM(M,L,NY,NX)*SOXYL(L,NY,NX)
        VOLPOX=VOLPM(M,L,NY,NX)
        VOLWPM=VOLWOX+VOLPOX
        !write(*,*)'VOLWPM=',VOLWPM
        DO 425 MX=1,NPT
          OXYG1=OXYG1+ROXYFX
          OXYS1=OXYS1+ROXYLX
          COXYS1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX) &
            ,AMAX1(0.0,safe_adb(OXYS1,(VOLWM(M,L,NY,NX)*FOXYX))))
          X=DIFOX*COXYS1
          IF(X.GT.ZEROS(NY,NX).AND.OXYS1.GT.ZEROS(NY,NX))THEN
            B=-RUPMX-DIFOX*OXKX-X
            C=X*RUPMX
            RMPOX=(-B-SQRT(B*B-4.0*C))/2.0_r8
          ELSE
            RMPOX=0.0_r8
          ENDIF
          OXYS1=OXYS1-RMPOX
          IF(THETPM(M,L,NY,NX).GT.THETX.AND.VOLPOX.GT.ZEROS(NY,NX))THEN
            ROXDFQ=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
          ELSE
            ROXDFQ=0.0_r8
          ENDIF
          OXYG1=OXYG1-ROXDFQ
          OXYS1=OXYS1+ROXDFQ
          RUPOX(NGL,N,K)=RUPOX(NGL,N,K)+RMPOX
          ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RMPOX

425     CONTINUE
420   CONTINUE
      !write(*,*)'420'
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
      !
      !     WFN=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RVMX2,RVMB2=NH3,NO2 oxidation in non-band, band
      !
      WFN(NGL,N,K)=AMIN1(1.0,AMAX1(0.0,RUPOX(NGL,N,K)/ROXYP(NGL,N,K)))
!     IF(K.LE.4)THEN
!       ROQCS(NGL,N,K,L,NY,NX)=ROQCS(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
!       ROQAS(NGL,N,K,L,NY,NX)=ROQAS(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
!       ROQCD(NGL,N,K)=ROQCD(NGL,N,K)*WFN(NGL,N,K)
!     ENDIF
      IF(K.EQ.5)THEN
        IF(N.EQ.1)THEN
          RVMX4(NGL,N,K,L,NY,NX)=RVMX4(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
          RVMB4(NGL,N,K,L,NY,NX)=RVMB4(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
        ELSEIF(N.EQ.2)THEN
          RVMX2(NGL,N,K,L,NY,NX)=RVMX2(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
          RVMB2(NGL,N,K,L,NY,NX)=RVMB2(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
        ENDIF
      ENDIF
    ELSE
      RUPOX(NGL,N,K)=ROXYP(NGL,N,K)
      WFN(NGL,N,K)=1.0
    ENDIF
  ELSE
    RUPOX(NGL,N,K)=0.0_r8
    WFN(NGL,N,K)=1.0
  ENDIF
  !write(*,*)'RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2'
  !
  !     RGOMO,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2X,RCH3X,RCH4X,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
  !     ROXYO=O2-limited O2 uptake
  !     RVOXA,RVOXB=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !
  RGOMO(NGL,N,K)=RGOMP*WFN(NGL,N,K)
  RCO2X(NGL,N,K)=RGOMO(NGL,N,K)
  RCH3X(NGL,N,K)=0.0_r8
  RCH4X(NGL,N,K)=0.0_r8
  ROXYO(NGL,N,K)=ROXYM(NGL,N,K)*WFN(NGL,N,K)
  RH2GX(NGL,N,K)=0.0_r8
  IF(K.EQ.5)THEN
    RVOXA(NGL,N)=RVOXPA*WFN(NGL,N,K)
    RVOXB(NGL,N)=RVOXPB*WFN(NGL,N,K)
  ENDIF
  !write(*,*)'finish AerobsO2Uptake'
  end associate
  end subroutine AerobsO2Uptake
!------------------------------------------------------------------------------------------

  subroutine BiomassMineralization(NGL,N,K,L,NY,NX,FNH4X,&
    FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in) :: FNH4X
  real(r8), intent(in) :: FNB3X,FNB4X,FNO3X
  real(r8), intent(in) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroMicStateType),intent(inout) :: nmics
  real(r8) :: CNH4X,CNH4Y,CNO3X,CNO3Y
  real(r8) :: CH2PX,CH2PY
  real(r8) :: CH1PX,CH1PY
  real(r8) :: FNH4S,FNHBS
  real(r8) :: FNO3S,FNO3B
  real(r8) :: FH1PS,FH1PB
  real(r8) :: FH2PS,FH2PB
  real(r8) :: H2POM,H2PBM
  real(r8) :: H1POM,H1PBM
  real(r8) :: H1P4M,H2P4M
  real(r8) :: RINHP
  real(r8) :: RINHX,RINOP,RINOX,RIPOP,RIPOX,RIP1P
  real(r8) :: RIP1X,RINHPR,RINOPR,RIPOPR,RIP1PR
  real(r8) :: ZNH4M,ZNHBM
  real(r8) :: ZNO3M
  real(r8) :: ZNOBM

!     begin_execution
  associate(                 &
   TFNG    => nmics%TFNG  ,  &
   OMA     => nmics%OMA   ,  &
   FNH4XR  => nmicf%FNH4XR,  &
   FNO3XR  => nmicf%FNO3XR,  &
   FPO4XR  => nmicf%FPO4XR,  &
   FP14XR  => nmicf%FP14XR,  &
   RINH4 => nmicf%RINH4   ,  &
   RINO3  => nmicf%RINO3  ,  &
   RIPO4  => nmicf%RIPO4  ,  &
   RINB4 => nmicf%RINB4   ,  &
   RINB3 => nmicf%RINB3   ,  &
   RIPOB => nmicf%RIPOB   ,  &
   RINH4R  => nmicf%RINH4R,  &
   RINO3R  => nmicf%RINO3R,  &
   RIPO4R  => nmicf%RIPO4R,  &
   RIP14  => nmicf%RIP14,    &
   RIP1B  => nmicf%RIP1B,    &
   RIP14R  => nmicf%RIP14R   &
  )
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
!     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMN=microbial nonstructural C,N
!     CNOMC=maximum microbial N:C ratio
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     RINHX=microbially limited NH4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHO,RINHB=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
!     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  FNH4S=VLNH4(L,NY,NX)
  FNHBS=VLNHB(L,NY,NX)
  RINHP=(OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K)-OMN(3,NGL,N,K,L,NY,NX))
  IF(RINHP.GT.0.0)THEN
    CNH4X=AMAX1(0.0_r8,CNH4S(L,NY,NX)-Z4MN)
    CNH4Y=AMAX1(0.0_r8,CNH4B(L,NY,NX)-Z4MN)
    RINHX=AMIN1(RINHP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*Z4MX)
    RINHO(NGL,N,K,L,NY,NX)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RINHB(NGL,N,K,L,NY,NX)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    ZNH4M=Z4MN*VOLW(L,NY,NX)*FNH4S
    ZNHBM=Z4MN*VOLW(L,NY,NX)*FNHBS
    RINH4(NGL,N,K)=AMIN1(FNH4X*AMAX1(0.0,(ZNH4S(L,NY,NX)-ZNH4M)) &
      ,RINHO(NGL,N,K,L,NY,NX))
    RINB4(NGL,N,K)=AMIN1(FNB4X*AMAX1(0.0,(ZNH4B(L,NY,NX)-ZNHBM)) &
      ,RINHB(NGL,N,K,L,NY,NX))
  ELSE
    RINHO(NGL,N,K,L,NY,NX)=0.0_r8
    RINHB(NGL,N,K,L,NY,NX)=0.0_r8
    RINH4(NGL,N,K)=RINHP*FNH4S
    RINB4(NGL,N,K)=RINHP*FNHBS
  ENDIF
  TRINH4(NY,NX)=TRINH4(NY,NX)+(RINH4(NGL,N,K)+RINB4(NGL,N,K))
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
!     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RINOP=NO3 immobilization (+ve) demand
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     min NO3 concentration and Km for NO3 uptake
!     RINOX=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RINOO,RINOB=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
!     TRINH4=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
  FNO3S=VLNO3(L,NY,NX)
  FNO3B=VLNOB(L,NY,NX)
  RINOP=AMAX1(0.0,RINHP-RINH4(NGL,N,K)-RINB4(NGL,N,K))
  IF(RINOP.GT.0.0)THEN
    CNO3X=AMAX1(0.0,CNO3S(L,NY,NX)-ZOMN)
    CNO3Y=AMAX1(0.0,CNO3B(L,NY,NX)-ZOMN)
    RINOX=AMIN1(RINOP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*ZOMX)
    RINOO(NGL,N,K,L,NY,NX)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RINOB(NGL,N,K,L,NY,NX)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M=ZOMN*VOLW(L,NY,NX)*FNO3S
    ZNOBM=ZOMN*VOLW(L,NY,NX)*FNO3B
    RINO3(NGL,N,K)=AMIN1(FNO3X*AMAX1(0.0_r8,(ZNO3S(L,NY,NX)-ZNO3M)) &
      ,RINOO(NGL,N,K,L,NY,NX))
    RINB3(NGL,N,K)=AMIN1(FNB3X*AMAX1(0.0_r8,(ZNO3B(L,NY,NX)-ZNOBM)) &
      ,RINOB(NGL,N,K,L,NY,NX))
  ELSE
    RINOO(NGL,N,K,L,NY,NX)=0.0_r8
    RINOB(NGL,N,K,L,NY,NX)=0.0_r8
    RINO3(NGL,N,K)=RINOP*FNO3S
    RINB3(NGL,N,K)=RINOP*FNO3B
  ENDIF
  TRINH4(NY,NX)=TRINH4(NY,NX)+(RINO3(NGL,N,K)+RINB3(NGL,N,K))
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMP=microbial nonstructural C,P
!     CPOMC=maximum microbial P:C ratio
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     min H2PO4 concentration and Km for H2PO4 uptake
!     RIPOX=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RIPOO,RIPBO=substrate-unlimited H2PO4 mineraln-immobiln
!     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
!     VOLW=water content
!     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
!     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
  FH2PS=VLPO4(L,NY,NX)
  FH2PB=VLPOB(L,NY,NX)
  RIPOP=(OMC(3,NGL,N,K,L,NY,NX)*CPOMC(3,NGL,N,K)-OMP(3,NGL,N,K,L,NY,NX))
  IF(RIPOP.GT.0.0)THEN
    CH2PX=AMAX1(0.0,CH2P4(L,NY,NX)-HPMN)
    CH2PY=AMAX1(0.0,CH2P4B(L,NY,NX)-HPMN)
    RIPOX=AMIN1(RIPOP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX)
    RIPOO(NGL,N,K,L,NY,NX)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RIPBO(NGL,N,K,L,NY,NX)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM=HPMN*VOLW(L,NY,NX)*FH2PS
    H2PBM=HPMN*VOLW(L,NY,NX)*FH2PB
    RIPO4(NGL,N,K)=AMIN1(FPO4X*AMAX1(0.0,(H2PO4(L,NY,NX)-H2POM)) &
      ,RIPOO(NGL,N,K,L,NY,NX))
    RIPOB(NGL,N,K)=AMIN1(FPOBX*AMAX1(0.0,(H2POB(L,NY,NX)-H2PBM)) &
      ,RIPBO(NGL,N,K,L,NY,NX))
  ELSE
    RIPOO(NGL,N,K,L,NY,NX)=0.0_r8
    RIPBO(NGL,N,K,L,NY,NX)=0.0_r8
    RIPO4(NGL,N,K)=RIPOP*FH2PS
    RIPOB(NGL,N,K)=RIPOP*FH2PB
  ENDIF
  TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIPO4(NGL,N,K)+RIPOB(NGL,N,K))
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RIP1P=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     min HPO4 concentration and Km for HPO4 uptake
!     RIP1X=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RIPO1,RIPB1=substrate-unlimited HPO4 mineraln-immobiln
!     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
!     VOLW=water content
!     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
!     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
!     TRIPO4=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
  FH1PS=VLPO4(L,NY,NX)
  FH1PB=VLPOB(L,NY,NX)
  RIP1P=0.1*AMAX1(0.0,RIPOP-RIPO4(NGL,N,K)-RIPOB(NGL,N,K))
  IF(RIP1P.GT.0.0)THEN
    CH1PX=AMAX1(0.0,CH1P4(L,NY,NX)-HPMN)
    CH1PY=AMAX1(0.0,CH1P4B(L,NY,NX)-HPMN)
    RIP1X=AMIN1(RIP1P,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX)
    RIPO1(NGL,N,K,L,NY,NX)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RIPB1(NGL,N,K,L,NY,NX)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM=HPMN*VOLW(L,NY,NX)*FH1PS
    H1PBM=HPMN*VOLW(L,NY,NX)*FH1PB
    RIP14(NGL,N,K)=AMIN1(FP14X*AMAX1(0.0,(H1PO4(L,NY,NX)-H1POM)) &
      ,RIPO1(NGL,N,K,L,NY,NX))
    RIP1B(NGL,N,K)=AMIN1(FP1BX*AMAX1(0.0,(H1POB(L,NY,NX)-H1PBM)) &
      ,RIPB1(NGL,N,K,L,NY,NX))
  ELSE
    RIPO1(NGL,N,K,L,NY,NX)=0.0_r8
    RIPB1(NGL,N,K,L,NY,NX)=0.0_r8
    RIP14(NGL,N,K)=RIP1P*FH1PS
    RIP1B(NGL,N,K)=RIP1P*FH1PB
  ENDIF
  TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIP14(NGL,N,K)+RIP1B(NGL,N,K))
!
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RINHPR=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RINHOR=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M=NH4 not available for uptake
!     FNH4XR=fractions of biological NH4 demand
!     RINH4R=substrate-limited NH4 mineraln-immobiln
!     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  IF(L.EQ.0)THEN
    RINHPR=RINHP-RINH4(NGL,N,K)-RINO3(NGL,N,K)
    IF(RINHPR.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(NU(NY,NX),NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(NU(NY,NX),NY,NX)-Z4MN)
      RINHOR(NGL,N,K,NY,NX)=AMIN1(RINHPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*Z4MX) &
        *(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VOLW(NU(NY,NX),NY,NX)
      RINH4R(NGL,N,K)=AMIN1(FNH4XR(NGL,N,K)*AMAX1(0.0 &
        ,(ZNH4T-ZNH4M)),RINHOR(NGL,N,K,NY,NX))
    ELSE
      RINHOR(NGL,N,K,NY,NX)=0.0_r8
      RINH4R(NGL,N,K)=RINHPR
    ENDIF
    TRINH4(NY,NX)=TRINH4(NY,NX)+RINH4R(NGL,N,K)
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RINOPR=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     minimum NO3 concentration and Km for NO3 uptake
!     RINOOR=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RINO3R=substrate-unlimited NO3 immobiln
!     VOLW=water content
!     ZNO3M=NO3 not available for uptake
!     FNO3XR=fraction of biological NO3 demand
!     RINO3R=substrate-limited NO3 immobiln
!     TRINH4=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
    RINOPR=AMAX1(0.0,RINHPR-RINH4R(NGL,N,K))
    IF(RINOPR.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(NU(NY,NX),NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(NU(NY,NX),NY,NX)-ZOMN)
      RINOOR(NGL,N,K,NY,NX)=AMAX1(RINOPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*ZOMX) &
        *(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VOLW(NU(NY,NX),NY,NX)
      RINO3R(NGL,N,K)=AMIN1(FNO3XR(NGL,N,K)*AMAX1(0.0 &
        ,(ZNO3T-ZNO3M)),RINOOR(NGL,N,K,NY,NX))
    ELSE
      RINOOR(NGL,N,K,NY,NX)=0.0_r8
      RINO3R(NGL,N,K)=RINOPR
    ENDIF
    TRINH4(NY,NX)=TRINH4(NY,NX)+RINO3R(NGL,N,K)
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RIPOPR=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     minimum H2PO4 concentration and Km for H2PO4 uptake
!     RIPOOR=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RIPOOR=substrate-unlimited H2PO4 mineraln-immobiln
!     VOLW=water content
!     H2P4M=H2PO4 not available for uptake
!     FPO4XR=fractions of biological H2PO4 demand
!     RIPO4R=substrate-limited H2PO4 mineraln-immobiln
!     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    RIPOPR=RIPOP-RIPO4(NGL,N,K)
    IF(RIPOPR.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(NU(NY,NX),NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPOOR(NGL,N,K,NY,NX)=AMIN1(RIPOPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX) &
        *(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M=HPMN*VOLW(NU(NY,NX),NY,NX)
      RIPO4R(NGL,N,K)=AMIN1(FPO4XR(NGL,N,K)*AMAX1(0.0 &
        ,(H2P4T-H2P4M)),RIPOOR(NGL,N,K,NY,NX))
    ELSE
      RIPOOR(NGL,N,K,NY,NX)=0.0_r8
      RIPO4R(NGL,N,K)=RIPOPR
    ENDIF
    TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIPO4R(NGL,N,K)
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RIP1PR=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     minimum HPO4 concentration and Km for HPO4 uptake
!     RIPO1R=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RIPO1R=substrate-unlimited HPO4 mineraln-immobiln
!     VOLW=water content
!     H1P4M=HPO4 not available for uptake
!     FP14XR=fraction of biological HPO4 demand
!     RIP14R=substrate-limited HPO4 minereraln-immobiln
!     TRIPO4=total HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4(L,NY,NX)
    FH1PB=VLPOB(L,NY,NX)
    RIP1PR=0.1*AMAX1(0.0,RIPOPR-RIPO4R(NGL,N,K))
    IF(RIP1PR.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(NU(NY,NX),NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPO1R(NGL,N,K,NY,NX)=AMIN1(RIP1PR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX) &
        *(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M=HPMN*VOLW(NU(NY,NX),NY,NX)
      RIP14R(NGL,N,K)=AMIN1(FP14XR(NGL,N,K)*AMAX1(0.0 &
        ,(H1P4T-H1P4M)),RIPO1R(NGL,N,K,NY,NX))
    ELSE
      RIPO1R(NGL,N,K,NY,NX)=0.0_r8
      RIP14R(NGL,N,K)=RIP1PR
    ENDIF
    TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIP14R(NGL,N,K)
  ENDIF
  end associate
  end subroutine BiomassMineralization
!------------------------------------------------------------------------------------------

  subroutine GatherMicrobialRespiration(NGL,N,K,L,NY,NX,RMOMK,RGOMT,RXOMT,RMOMT,&
    nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(out) :: RGOMT,RXOMT
  real(r8), intent(out) :: RMOMT
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P
!     begin_execution
  associate(                  &
    TFNR  => nmics%TFNR,      &
    OMN2  => nmics%OMN2,      &
    RGN2F => nmicf%RGN2F,     &
    RGOMO  => nmicf%RGOMO,    &
    RMOMC => nmicf%RMOMC ,    &
    RINH4 => nmicf%RINH4 ,    &
    RINO3  => nmicf%RINO3,    &
    RN2FX  => nmicf%RN2FX     &
  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TFNR=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!     RMOMK=effect of low microbial C concentration on mntc respn
!
  FPH=1.0+AMAX1(0.0,0.25*(6.5-PH(L,NY,NX)))
  RMOMX=RMOM*TFNR(NGL,N,K)*FPH
  RMOMC(1,NGL,N,K)=OMN(1,NGL,N,K,L,NY,NX)*RMOMX*RMOMK(1)
  RMOMC(2,NGL,N,K)=OMN2(NGL,N,K)*RMOMX*RMOMK(2)
!
!     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
!
!     RMOMT=total maintenance respiration
!     RGOMT=growth respiration
!     RXOMT=senescence respiration
!
  RMOMT=RMOMC(1,NGL,N,K)+RMOMC(2,NGL,N,K)
  RGOMT=AMAX1(0.0,RGOMO(NGL,N,K)-RMOMT)
  RXOMT=AMAX1(0.0,RMOMT-RGOMO(NGL,N,K))
!
!     N2 FIXATION: N=(6) AEROBIC, (7) ANAEROBIC
!     FROM GROWTH RESPIRATION, FIXATION ENERGY REQUIREMENT,
!     MICROBIAL N REQUIREMENT IN LABILE (1) AND
!     RESISTANT (2) FRACTIONS
!
!     RGN2P=respiration to meet N2 fixation demand
!     OMC,OMN=microbial nonstructural C,N
!     CNOMC=maximum microbial N:C ratio
!     EN2F=N2 fixation yield per unit nonstructural C
!     RGOMT=growth respiration
!     RGN2F=respiration for N2 fixation
!     CZ2GS=aqueous N2 concentration
!     ZFKM=Km for N2 uptake
!     OMGR*OMC(3,NGL,N,K,L,NY,NX)=nonstructural C limitation to RGN2F
!     RN2FX=N2 fixation rate
!
  IF(K.LE.4.AND.(N.EQ.6.OR.N.EQ.7))THEN
    RGN2P=AMAX1(0.0,OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K) &
      -OMN(3,NGL,N,K,L,NY,NX))/EN2F(N)
    IF(RGOMT.GT.ZEROS(NY,NX))THEN
      RGN2F(NGL,N,K)=AMIN1(RGOMT*RGN2P/(RGOMT+RGN2P) &
        *CZ2GS(L,NY,NX)/(CZ2GS(L,NY,NX)+ZFKM) &
        ,OMGR*OMC(3,NGL,N,K,L,NY,NX))
    ELSE
      RGN2F(NGL,N,K)=0.0_r8
    ENDIF
    RN2FX(NGL,N,K)=RGN2F(NGL,N,K)*EN2F(N)
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,5566)'N2 FIX',I,J,NX,NY,L,K,N,RN2FX(NGL,N,K),EN2F(N)
!    2,OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K),OMN(3,NGL,N,K,L,NY,NX)
!    3,RINH4(NGL,N,K),RINO3(NGL,N,K),RGN2P,RGN2F(NGL,N,K),RGOMT
!    4,CZ2GS(L,NY,NX)
!5566  FORMAT(A8,7I4,30E12.4)
!     ENDIF
  ELSE
    RN2FX(NGL,N,K)=0.0_r8
    RGN2F(NGL,N,K)=0.0_r8
  ENDIF
  end associate
  end subroutine GatherMicrobialRespiration
!------------------------------------------------------------------------------------------

  subroutine GetMicrobialAnabolismFlux(NGL,N,K,L,NY,NX,ECHZ,FGOCP,FGOAP,&
    RGOMT,RXOMT,RMOMT,spomk,rmomk,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: NGL,N,K,L,NY,NX
  real(r8), intent(in) :: ECHZ
  real(r8), intent(in) :: FGOCP,FGOAP,RGOMT,RXOMT,RMOMT
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(NitroMicStateType), intent(inout) :: nmics
  type(NitroMicFluxType), intent(inout) :: nmicf
  type(NitroOMcplxFluxType), intent(inout) :: ncplxf
  type(NitroOMcplxStateType), intent(inout) :: ncplxs
  integer :: M
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: CXC
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FRM
!     begin_execution
  associate(                  &
    CNOMA  => nmics%CNOMA,    &
    CPOMA  => nmics%CPOMA,    &
    TFNG   => nmics%TFNG ,    &
    FCN    => nmics%FCN  ,    &
    FCP    => nmics%FCP  ,    &
    FOMK   => nmics%FOMK ,    &
    CGOMC  => nmicf%CGOMC,    &
    CGOMN  => nmicf%CGOMN,    &
    CGOMP  => nmicf%CGOMP,    &
    CGOQC  => nmicf%CGOQC,    &
    CGOAC  => nmicf%CGOAC,    &
    CGOMS  => nmicf%CGOMS,    &
    CGONS  => nmicf%CGONS,    &
    CGOPS  => nmicf%CGOPS,    &
    RGOMO  => nmicf%RGOMO,    &
    RGOMD  => nmicf%RGOMD,    &
    RMOMC => nmicf%RMOMC ,    &
    RDOMC  => nmicf%RDOMC,    &
    RDOMN => nmicf%RDOMN ,    &
    RDOMP => nmicf%RDOMP ,    &
    RHOMC => nmicf%RHOMC ,    &
    RHOMN => nmicf%RHOMN ,    &
    RHOMP => nmicf%RHOMP ,    &
    RCOMC  => nmicf%RCOMC,    &
    RCOMN  => nmicf%RCOMN,    &
    RCOMP  => nmicf%RCOMP,    &
    RDMMC  => nmicf%RDMMC,    &
    RHMMC  => nmicf%RHMMC,    &
    RCMMC  => nmicf%RCMMC,    &
    RDMMN  => nmicf%RDMMN,    &
    RHMMN  => nmicf%RHMMN,    &
    RCMMN  => nmicf%RCMMN,    &
    RDMMP  => nmicf%RDMMP,    &
    RHMMP  => nmicf%RHMMP,    &
    RCMMP  => nmicf%RCMMP,    &
    RXOMC  => nmicf%RXOMC,    &
    RXOMN  => nmicf%RXOMN,    &
    RXOMP  => nmicf%RXOMP,    &
    R3OMC  => nmicf%R3OMC,    &
    R3OMN  => nmicf%R3OMN,    &
    R3OMP  => nmicf%R3OMP,    &
    RXMMC  => nmicf%RXMMC,    &
    RXMMN  => nmicf%RXMMN,    &
    RXMMP  => nmicf%RXMMP,    &
    R3MMC  => nmicf%R3MMC,    &
    R3MMN  => nmicf%R3MMN,    &
    R3MMP  => nmicf%R3MMP,    &
    RCO2X  => nmicf%RCO2X,    &
    RGN2F  => nmicf%RGN2F,    &
    TCGOQC  => ncplxf%TCGOQC, &
    TCGOAC  => ncplxf%TCGOAC, &
    TCGOMN  => ncplxf%TCGOMN, &
    TCGOMP  => ncplxf%TCGOMP, &
    CNQ  => ncplxs%CNQ,       &
    CPQ  => ncplxs%CPQ        &
  )

!     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
!     FROM O2, NOX AND C REDUCTION
!
!     CGOMX=DOC+acetate uptake from aerobic growth respiration
!     CGOMD=DOC+acetate uptake from denitrifier growth respiration
!     RMOMT=maintenance respiration
!     RGOMO=total respiration
!     RGOMD=respiration for denitrifcation
!     RGN2F=respiration for N2 fixation
!     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
!     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake(heterotrophs
!     CGOMC=total CO2,CH4 uptake (autotrophs)
!     CGOMN,CGOMP=DON, DOP uptake
!     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
!     OQN,OPQ=DON,DOP
!     FOMK=faction of OMA in total OMA
!     CNQ,CPQ=DON/DOC, DOP/DOC
!     FCN,FCP=limitation from N,P
!

  CGOMX=AMIN1(RMOMT,RGOMO(NGL,N,K))+RGN2F(NGL,N,K) &
    +(RGOMT-RGN2F(NGL,N,K))/ECHZ
  CGOMD=RGOMD(NGL,N,K)/ENOX
  CGOMC(NGL,N,K)=CGOMX+CGOMD
  IF(K.LE.4)THEN
    CGOQC(NGL,N,K)=CGOMX*FGOCP+CGOMD
    CGOAC(NGL,N,K)=CGOMX*FGOAP
    CGOXC=CGOQC(NGL,N,K)+CGOAC(NGL,N,K)
    CGOMN(NGL,N,K)=AMAX1(0.0,AMIN1(OQN(K,L,NY,NX)*FOMK(NGL,N,K) &
      ,CGOXC*CNQ(K)/FCN(NGL,N,K)))
    CGOMP(NGL,N,K)=AMAX1(0.0,AMIN1(OQP(K,L,NY,NX)*FOMK(NGL,N,K) &
      ,CGOXC*CPQ(K)/FCP(NGL,N,K)))
  ELSE
    CGOQC(NGL,N,K)=CGOMX+CGOMD
    CGOAC(NGL,N,K)=0.0_r8
    CGOMN(NGL,N,K)=0.0_r8
    CGOMP(NGL,N,K)=0.0_r8
  ENDIF
  TCGOQC(K)=TCGOQC(K)+CGOQC(NGL,N,K)
  TCGOAC(K)=TCGOAC(K)+CGOAC(NGL,N,K)
  TCGOMN(K)=TCGOMN(K)+CGOMN(NGL,N,K)
  TCGOMP(K)=TCGOMP(K)+CGOMP(NGL,N,K)
!
!     TRANSFER UPTAKEN C,N,P FROM STORAGE TO ACTIVE BIOMASS
!
!     OMC,OMN,OMP=nonstructural C,N,P
!     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
!     CNOMC,CPOMC=maximum microbial N:C, P:C ratios
!     RCCC,RCCN,RCCP=C,N,P recycling fractions
!     RCCZ,RCCY=min, max C recycling fractions
!     RCCX,RCCQ=max N,P recycling fractions
!
  IF(OMC(3,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX) &
    .AND.OMC(1,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
    CCC=AMAX1(0.0,AMIN1(1.0 &
      ,OMN(3,NGL,N,K,L,NY,NX)/(OMN(3,NGL,N,K,L,NY,NX) &
      +OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K)) &
      ,OMP(3,NGL,N,K,L,NY,NX)/(OMP(3,NGL,N,K,L,NY,NX) &
      +OMC(3,NGL,N,K,L,NY,NX)*CPOMC(3,NGL,N,K))))
    CXC=OMC(3,NGL,N,K,L,NY,NX)/OMC(1,NGL,N,K,L,NY,NX)
    C3C=1.0/(1.0+CXC/CKC)
    CNC=AMAX1(0.0,AMIN1(1.0 &
      ,OMC(3,NGL,N,K,L,NY,NX)/(OMC(3,NGL,N,K,L,NY,NX) &
      +OMN(3,NGL,N,K,L,NY,NX)/CNOMC(3,NGL,N,K))))
    CPC=AMAX1(0.0,AMIN1(1.0 &
      ,OMC(3,NGL,N,K,L,NY,NX)/(OMC(3,NGL,N,K,L,NY,NX) &
      +OMP(3,NGL,N,K,L,NY,NX)/CPOMC(3,NGL,N,K))))
    RCCC=RCCZ+AMAX1(CCC,C3C)*RCCY
    RCCN=CNC*RCCX
    RCCP=CPC*RCCQ
  ELSE
    RCCC=RCCZ
    RCCN=0.0_r8
    RCCP=0.0_r8
  ENDIF
!     IF((I/120)*120.EQ.I.AND.J.EQ.24)THEN
!     WRITE(*,5555)'RCCC',I,J,NX,NY,L,K,N,RCCC,RCCN,RCCP
!    2,OMC(3,NGL,N,K,L,NY,NX),OMN(3,NGL,N,K,L,NY,NX),OMP(3,NGL,N,K,L,NY,NX)
!    3,CCC,C3C,CNC,CPC
!     ENDIF
!
!     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
!
!     CGOMZ=transfer from nonstructural to structural microbial C
!     TFNG=temperature+water stress function
!     OMGR=rate constant for transferring nonstructural to structural C
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     FL=partitioning between labile and resistant microbial components
!     OMC,OMN,OMP=nonstructural microbial C,N,P
!
  CGOMZ=TFNG(NGL,N,K)*OMGR*AMAX1(0.0,OMC(3,NGL,N,K,L,NY,NX))
  DO 745 M=1,2
    CGOMS(M,NGL,N,K)=FL(M)*CGOMZ
    IF(OMC(3,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CGONS(M,NGL,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMN(3,NGL,N,K,L,NY,NX)) &
        ,CGOMS(M,NGL,N,K)*OMN(3,NGL,N,K,L,NY,NX)/OMC(3,NGL,N,K,L,NY,NX))
      CGOPS(M,NGL,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMP(3,NGL,N,K,L,NY,NX)) &
        ,CGOMS(M,NGL,N,K)*OMP(3,NGL,N,K,L,NY,NX)/OMC(3,NGL,N,K,L,NY,NX))
    ELSE
      CGONS(M,NGL,N,K)=0.0_r8
      CGOPS(M,NGL,N,K)=0.0_r8
    ENDIF
!
!     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
!     RATE, TEMPERATURE
!
!     SPOMX=rate constant for microbial decomposition
!     SPOMC=basal decomposition rate
!     SPOMK=effect of low microbial C concentration on microbial decay
!     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
!     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall
!     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
!
    SPOMX=SQRT(TFNG(NGL,N,K))*SPOMC(M)*SPOMK(M)
    RXOMC(M,NGL,N,K)=AMAX1(0.0,OMC(M,NGL,N,K,L,NY,NX)*SPOMX)
    RXOMN(M,NGL,N,K)=AMAX1(0.0,OMN(M,NGL,N,K,L,NY,NX)*SPOMX)
    RXOMP(M,NGL,N,K)=AMAX1(0.0,OMP(M,NGL,N,K,L,NY,NX)*SPOMX)
    RDOMC(M,NGL,N,K)=RXOMC(M,NGL,N,K)*(1.0-RCCC)
    RDOMN(M,NGL,N,K)=RXOMN(M,NGL,N,K)*(1.0-RCCC)*(1.0-RCCN)
    RDOMP(M,NGL,N,K)=RXOMP(M,NGL,N,K)*(1.0-RCCC)*(1.0-RCCP)
    R3OMC(M,NGL,N,K)=RXOMC(M,NGL,N,K)-RDOMC(M,NGL,N,K)
    R3OMN(M,NGL,N,K)=RXOMN(M,NGL,N,K)-RDOMN(M,NGL,N,K)
    R3OMP(M,NGL,N,K)=RXOMP(M,NGL,N,K)-RDOMP(M,NGL,N,K)
!
!     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
!     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
!
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
!     EHUM=humus transfer fraction from hour1.f
!     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
!
    RHOMC(M,NGL,N,K)=AMAX1(0.0,RDOMC(M,NGL,N,K)*EHUM(L,NY,NX))
    RHOMN(M,NGL,N,K)=AMAX1(0.0,RDOMN(M,NGL,N,K)*EHUM(L,NY,NX))
    RHOMP(M,NGL,N,K)=AMAX1(0.0,RDOMP(M,NGL,N,K)*EHUM(L,NY,NX))
!     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
!     WRITE(*,8821)'RHOMC',I,J,L,K,N,M
!    2,RXOMC(M,NGL,N,K),RXOMN(M,NGL,N,K),RXOMP(M,NGL,N,K)
!    2,RDOMC(M,NGL,N,K),RDOMN(M,NGL,N,K),RDOMP(M,NGL,N,K)
!    2,R3OMC(M,NGL,N,K),R3OMN(M,NGL,N,K),R3OMP(M,NGL,N,K)
!    2,RHOMC(M,NGL,N,K),RHOMN(M,NGL,N,K),RHOMP(M,NGL,N,K)
!    4,OMC(M,NGL,N,K,L,NY,NX),OMN(M,NGL,N,K,L,NY,NX)
!    5,OMP(M,NGL,N,K,L,NY,NX)
!    4,OMC(3,NGL,N,K,L,NY,NX),OMN(3,NGL,N,K,L,NY,NX)
!    5,OMP(3,NGL,N,K,L,NY,NX)
!    6,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
!    2,SPOMX,RCCC,RCCN,RCCP
!     ENDIF
!
!     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
!
    RCOMC(M,NGL,N,K)=RDOMC(M,NGL,N,K)-RHOMC(M,NGL,N,K)
    RCOMN(M,NGL,N,K)=RDOMN(M,NGL,N,K)-RHOMN(M,NGL,N,K)
    RCOMP(M,NGL,N,K)=RDOMP(M,NGL,N,K)-RHOMP(M,NGL,N,K)
745 CONTINUE
!
!     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
!     EXCEEDS UPTAKE
!
!     OMC,OMN,OMP=microbial C,N,P
!     RMOMT=total maintenance respiration
!     RXOMT=senescence respiration
!     RCCC=C recycling fraction
!     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
!     RMOMC=maintenance respiration
!     CNOMA,CPOMA=N:C,P:C ratios of active biomass
!     RDMMC,RDMMN,RDMMP=microbial C,N,P litterfall from senescence
!     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
!
  IF(RXOMT.GT.ZEROS(NY,NX).AND.RMOMT.GT.ZEROS(NY,NX) &
    .AND.RCCC.GT.ZERO)THEN
    FRM=RXOMT/RMOMT
    DO 730 M=1,2
      RXMMC(M,NGL,N,K)=AMIN1(OMC(M,NGL,N,K,L,NY,NX) &
        ,AMAX1(0.0,FRM*RMOMC(M,NGL,N,K)/RCCC))
      RXMMN(M,NGL,N,K)=AMIN1(OMN(M,NGL,N,K,L,NY,NX) &
        ,AMAX1(0.0,RXMMC(M,NGL,N,K)*CNOMA(NGL,N,K)))
      RXMMP(M,NGL,N,K)=AMIN1(OMP(M,NGL,N,K,L,NY,NX) &
        ,AMAX1(0.0,RXMMC(M,NGL,N,K)*CPOMA(NGL,N,K)))
      RDMMC(M,NGL,N,K)=RXMMC(M,NGL,N,K)*(1.0-RCCC)
      RDMMN(M,NGL,N,K)=RXMMN(M,NGL,N,K)*(1.0-RCCN)*(1.0-RCCC)
      RDMMP(M,NGL,N,K)=RXMMP(M,NGL,N,K)*(1.0-RCCP)*(1.0-RCCC)
      R3MMC(M,NGL,N,K)=RXMMC(M,NGL,N,K)-RDMMC(M,NGL,N,K)
      R3MMN(M,NGL,N,K)=RXMMN(M,NGL,N,K)-RDMMN(M,NGL,N,K)
      R3MMP(M,NGL,N,K)=RXMMP(M,NGL,N,K)-RDMMP(M,NGL,N,K)
!
!     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
!     PRODUCTS
!
!     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
!     EHUM=humus transfer fraction
!     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
!
      RHMMC(M,NGL,N,K)=AMAX1(0.0,RDMMC(M,NGL,N,K)*EHUM(L,NY,NX))
      RHMMN(M,NGL,N,K)=AMAX1(0.0,RDMMN(M,NGL,N,K)*EHUM(L,NY,NX))
      RHMMP(M,NGL,N,K)=AMAX1(0.0,RDMMP(M,NGL,N,K)*EHUM(L,NY,NX))
      RCMMC(M,NGL,N,K)=RDMMC(M,NGL,N,K)-RHMMC(M,NGL,N,K)
      RCMMN(M,NGL,N,K)=RDMMN(M,NGL,N,K)-RHMMN(M,NGL,N,K)
      RCMMP(M,NGL,N,K)=RDMMP(M,NGL,N,K)-RHMMP(M,NGL,N,K)

730 CONTINUE
  ELSE
    DO 720 M=1,2
      RXMMC(M,NGL,N,K)=0.0_r8
      RXMMN(M,NGL,N,K)=0.0_r8
      RXMMP(M,NGL,N,K)=0.0_r8
      RDMMC(M,NGL,N,K)=0.0_r8
      RDMMN(M,NGL,N,K)=0.0_r8
      RDMMP(M,NGL,N,K)=0.0_r8
      R3MMC(M,NGL,N,K)=0.0_r8
      R3MMN(M,NGL,N,K)=0.0_r8
      R3MMP(M,NGL,N,K)=0.0_r8
      RHMMC(M,NGL,N,K)=0.0_r8
      RHMMN(M,NGL,N,K)=0.0_r8
      RHMMP(M,NGL,N,K)=0.0_r8
      RCMMC(M,NGL,N,K)=0.0_r8
      RCMMN(M,NGL,N,K)=0.0_r8
      RCMMP(M,NGL,N,K)=0.0_r8
720 CONTINUE
  ENDIF
  end associate
  end subroutine GetMicrobialAnabolismFlux

end module nitro1LayerMod
