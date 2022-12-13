module TrnsfrsMod
!!
! Description:
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb,AZMAX1
  use SOMDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use EcoSIMSolverPar
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use FertilizerDataType
  use SoilWaterDataType
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use ChemTranspDataType
  use LandSurfDataType
  use SoilBGCDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use IrrigationDataType
  use PlantDataRateType
  use GridDataType
  use IngridTranspMod
  use TrnsfrsDataMod
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8) :: RCHQF,RCHGFU,RCHGFT
  real(r8) :: XN
  real(r8),allocatable ::  trcsa_TQR(:,:,:)                         !

  real(r8),allocatable ::  trcsa_TQ(:,:,:)                         !
  real(r8),allocatable ::  TALFLS(:,:,:)                      !
  real(r8),allocatable ::  TFEFLS(:,:,:)                      !
  real(r8),allocatable ::  TCAFLS(:,:,:)                      !
  real(r8),allocatable ::  THYFLS(:,:,:)                      !
  real(r8),allocatable ::  TMGFLS(:,:,:)                      !
  real(r8),allocatable ::  TNAFLS(:,:,:)                      !
  real(r8),allocatable ::  TKAFLS(:,:,:)                      !
  real(r8),allocatable ::  TOHFLS(:,:,:)                      !
  real(r8),allocatable ::  TSOFLS(:,:,:)                      !
  real(r8),allocatable ::  TCLFLS(:,:,:)                      !
  real(r8),allocatable ::  TC3FLS(:,:,:)                      !
  real(r8),allocatable ::  THCFLS(:,:,:)                      !
  real(r8),allocatable ::  TAL1FS(:,:,:)                      !
  real(r8),allocatable ::  TAL2FS(:,:,:)                      !
  real(r8),allocatable ::  TAL3FS(:,:,:)                      !
  real(r8),allocatable ::  TAL4FS(:,:,:)                      !
  real(r8),allocatable ::  TALSFS(:,:,:)                      !
  real(r8),allocatable ::  TFE1FS(:,:,:)                      !
  real(r8),allocatable ::  TFE2FS(:,:,:)                      !
  real(r8),allocatable ::  TFE3FS(:,:,:)                      !
  real(r8),allocatable ::  TFE4FS(:,:,:)                      !
  real(r8),allocatable ::  TFESFS(:,:,:)                      !
  real(r8),allocatable ::  TCAOFS(:,:,:)                      !
  real(r8),allocatable ::  TCACFS(:,:,:)                      !
  real(r8),allocatable ::  TCAHFS(:,:,:)                      !
  real(r8),allocatable ::  TCASFS(:,:,:)                      !
  real(r8),allocatable ::  TMGOFS(:,:,:)                      !
  real(r8),allocatable ::  TMGCFS(:,:,:)                      !
  real(r8),allocatable ::  TMGHFS(:,:,:)                      !
  real(r8),allocatable ::  TMGSFS(:,:,:)                      !
  real(r8),allocatable ::  TNACFS(:,:,:)                      !
  real(r8),allocatable ::  TNASFS(:,:,:)                      !
  real(r8),allocatable ::  TKASFS(:,:,:)                      !
  real(r8),allocatable ::  TH0PFS(:,:,:)                      !
  real(r8),allocatable ::  TH3PFS(:,:,:)                      !
  real(r8),allocatable ::  TF1PFS(:,:,:)                      !
  real(r8),allocatable ::  TF2PFS(:,:,:)                      !
  real(r8),allocatable ::  TC0PFS(:,:,:)                      !
  real(r8),allocatable ::  TC1PFS(:,:,:)                      !
  real(r8),allocatable ::  TC2PFS(:,:,:)                      !
  real(r8),allocatable ::  TM1PFS(:,:,:)                      !
  real(r8),allocatable ::  TH0BFB(:,:,:)                      !
  real(r8),allocatable ::  TH3BFB(:,:,:)                      !
  real(r8),allocatable ::  TF1BFB(:,:,:)                      !
  real(r8),allocatable ::  TF2BFB(:,:,:)                      !
  real(r8),allocatable ::  TC0BFB(:,:,:)                      !
  real(r8),allocatable ::  TC1BFB(:,:,:)                      !
  real(r8),allocatable ::  TC2BFB(:,:,:)                      !
  real(r8),allocatable ::  TM1BFB(:,:,:)                      !
  real(r8),allocatable ::  TALFHS(:,:,:)                      !
  real(r8),allocatable ::  TFEFHS(:,:,:)                      !
  real(r8),allocatable ::  THYFHS(:,:,:)                      !
  real(r8),allocatable ::  TCAFHS(:,:,:)                      !
  real(r8),allocatable ::  TMGFHS(:,:,:)                      !
  real(r8),allocatable ::  TNAFHS(:,:,:)                      !
  real(r8),allocatable ::  TKAFHS(:,:,:)                      !
  real(r8),allocatable ::  TOHFHS(:,:,:)                      !
  real(r8),allocatable ::  TSOFHS(:,:,:)                      !
  real(r8),allocatable ::  TCLFHS(:,:,:)                      !
  real(r8),allocatable ::  TC3FHS(:,:,:)                      !
  real(r8),allocatable ::  THCFHS(:,:,:)                      !
  real(r8),allocatable ::  TAL1HS(:,:,:)                      !
  real(r8),allocatable ::  TAL2HS(:,:,:)                      !
  real(r8),allocatable ::  TAL3HS(:,:,:)                      !
  real(r8),allocatable ::  TAL4HS(:,:,:)                      !
  real(r8),allocatable ::  TALSHS(:,:,:)                      !
  real(r8),allocatable ::  TFE1HS(:,:,:)                      !
  real(r8),allocatable ::  TFE2HS(:,:,:)                      !
  real(r8),allocatable ::  TFE3HS(:,:,:)                      !
  real(r8),allocatable ::  TFE4HS(:,:,:)                      !
  real(r8),allocatable ::  TFESHS(:,:,:)                      !
  real(r8),allocatable ::  TCAOHS(:,:,:)                      !
  real(r8),allocatable ::  TCACHS(:,:,:)                      !
  real(r8),allocatable ::  TCAHHS(:,:,:)                      !
  real(r8),allocatable ::  TCASHS(:,:,:)                      !
  real(r8),allocatable ::  TMGOHS(:,:,:)                      !
  real(r8),allocatable ::  TMGCHS(:,:,:)                      !
  real(r8),allocatable ::  TMGHHS(:,:,:)                      !
  real(r8),allocatable ::  TMGSHS(:,:,:)                      !
  real(r8),allocatable ::  TNACHS(:,:,:)                      !
  real(r8),allocatable ::  TNASHS(:,:,:)                      !
  real(r8),allocatable ::  TKASHS(:,:,:)                      !
  real(r8),allocatable ::  TH0PHS(:,:,:)                      !
  real(r8),allocatable ::  TH3PHS(:,:,:)                      !
  real(r8),allocatable ::  TF1PHS(:,:,:)                      !
  real(r8),allocatable ::  TF2PHS(:,:,:)                      !
  real(r8),allocatable ::  TC0PHS(:,:,:)                      !
  real(r8),allocatable ::  TC1PHS(:,:,:)                      !
  real(r8),allocatable ::  TC2PHS(:,:,:)                      !
  real(r8),allocatable ::  TM1PHS(:,:,:)                      !
  real(r8),allocatable ::  TH0BHB(:,:,:)                      !
  real(r8),allocatable ::  TH3BHB(:,:,:)                      !
  real(r8),allocatable ::  TF1BHB(:,:,:)                      !
  real(r8),allocatable ::  TF2BHB(:,:,:)                      !
  real(r8),allocatable ::  TC0BHB(:,:,:)                      !
  real(r8),allocatable ::  TC1BHB(:,:,:)                      !
  real(r8),allocatable ::  TC2BHB(:,:,:)                      !
  real(r8),allocatable ::  TM1BHB(:,:,:)                      !
  real(r8),allocatable ::  trcsa_RFLZ(:,:,:,:)                      !
!----------------------------------------------------------------------



  public :: trnsfrs
  public :: initTrnsfrs
  public :: destructTrnsfrs
  contains
!----------------------------------------------------------------------
  subroutine DestructTrnsfrs
  use abortutils, only : destroy
  implicit none

  call destroy(trcsa_TQR)


  call destroy(TALFLS)
  call destroy(TFEFLS)
  call destroy(TCAFLS)
  call destroy(THYFLS)
  call destroy(TMGFLS)
  call destroy(TNAFLS)
  call destroy(TKAFLS)
  call destroy(TOHFLS)
  call destroy(TSOFLS)
  call destroy(TCLFLS)
  call destroy(TC3FLS)
  call destroy(THCFLS)
  call destroy(TAL1FS)
  call destroy(TAL2FS)
  call destroy(TAL3FS)
  call destroy(TAL4FS)
  call destroy(TALSFS)
  call destroy(TFE1FS)
  call destroy(TFE2FS)
  call destroy(TFE3FS)
  call destroy(TFE4FS)
  call destroy(TFESFS)
  call destroy(TCAOFS)
  call destroy(TCACFS)
  call destroy(TCAHFS)
  call destroy(TCASFS)
  call destroy(TMGOFS)
  call destroy(TMGCFS)
  call destroy(TMGHFS)
  call destroy(TMGSFS)
  call destroy(TNACFS)
  call destroy(TNASFS)
  call destroy(TKASFS)
  call destroy(TH0PFS)
  call destroy(TH3PFS)
  call destroy(TF1PFS)
  call destroy(TF2PFS)
  call destroy(TC0PFS)
  call destroy(TC1PFS)
  call destroy(TC2PFS)
  call destroy(TM1PFS)
  call destroy(TH0BFB)
  call destroy(TH3BFB)
  call destroy(TF1BFB)
  call destroy(TF2BFB)
  call destroy(TC0BFB)
  call destroy(TC1BFB)
  call destroy(TC2BFB)
  call destroy(TM1BFB)
  call destroy(TALFHS)
  call destroy(TFEFHS)
  call destroy(THYFHS)
  call destroy(TCAFHS)
  call destroy(TMGFHS)
  call destroy(TNAFHS)
  call destroy(TKAFHS)
  call destroy(TOHFHS)
  call destroy(TSOFHS)
  call destroy(TCLFHS)
  call destroy(TC3FHS)
  call destroy(THCFHS)
  call destroy(TAL1HS)
  call destroy(TAL2HS)
  call destroy(TAL3HS)
  call destroy(TAL4HS)
  call destroy(TALSHS)
  call destroy(TFE1HS)
  call destroy(TFE2HS)
  call destroy(TFE3HS)
  call destroy(TFE4HS)
  call destroy(TFESHS)
  call destroy(TCAOHS)
  call destroy(TCACHS)
  call destroy(TCAHHS)
  call destroy(TCASHS)
  call destroy(TMGOHS)
  call destroy(TMGCHS)
  call destroy(TMGHHS)
  call destroy(TMGSHS)
  call destroy(TNACHS)
  call destroy(TNASHS)
  call destroy(TKASHS)
  call destroy(TH0PHS)
  call destroy(TH3PHS)
  call destroy(TF1PHS)
  call destroy(TF2PHS)
  call destroy(TC0PHS)
  call destroy(TC1PHS)
  call destroy(TC2PHS)
  call destroy(TM1PHS)
  call destroy(TH0BHB)
  call destroy(TH3BHB)
  call destroy(TF1BHB)
  call destroy(TF2BHB)
  call destroy(TC0BHB)
  call destroy(TC1BHB)
  call destroy(TC2BHB)
  call destroy(TM1BHB)

  end subroutine DestructTrnsfrs

!----------------------------------------------------------------------

  subroutine initTrnsfrs()
  implicit none

  allocate(trcsa_TQR(idsa_beg:idsa_end,JY,JX));       trcsa_TQR=0._r8

  allocate(TALFLS(JZ,JY,JX));   TALFLS=0._r8
  allocate(TFEFLS(JZ,JY,JX));   TFEFLS=0._r8
  allocate(TCAFLS(JZ,JY,JX));   TCAFLS=0._r8
  allocate(THYFLS(JZ,JY,JX));   THYFLS=0._r8
  allocate(TMGFLS(JZ,JY,JX));   TMGFLS=0._r8
  allocate(TNAFLS(JZ,JY,JX));   TNAFLS=0._r8
  allocate(TKAFLS(JZ,JY,JX));   TKAFLS=0._r8
  allocate(TOHFLS(JZ,JY,JX));   TOHFLS=0._r8
  allocate(TSOFLS(JZ,JY,JX));   TSOFLS=0._r8
  allocate(TCLFLS(JZ,JY,JX));   TCLFLS=0._r8
  allocate(TC3FLS(JZ,JY,JX));   TC3FLS=0._r8
  allocate(THCFLS(JZ,JY,JX));   THCFLS=0._r8
  allocate(TAL1FS(JZ,JY,JX));   TAL1FS=0._r8
  allocate(TAL2FS(JZ,JY,JX));   TAL2FS=0._r8
  allocate(TAL3FS(JZ,JY,JX));   TAL3FS=0._r8
  allocate(TAL4FS(JZ,JY,JX));   TAL4FS=0._r8
  allocate(TALSFS(JZ,JY,JX));   TALSFS=0._r8
  allocate(TFE1FS(JZ,JY,JX));   TFE1FS=0._r8
  allocate(TFE2FS(JZ,JY,JX));   TFE2FS=0._r8
  allocate(TFE3FS(JZ,JY,JX));   TFE3FS=0._r8
  allocate(TFE4FS(JZ,JY,JX));   TFE4FS=0._r8
  allocate(TFESFS(JZ,JY,JX));   TFESFS=0._r8
  allocate(TCAOFS(JZ,JY,JX));   TCAOFS=0._r8
  allocate(TCACFS(JZ,JY,JX));   TCACFS=0._r8
  allocate(TCAHFS(JZ,JY,JX));   TCAHFS=0._r8
  allocate(TCASFS(JZ,JY,JX));   TCASFS=0._r8
  allocate(TMGOFS(JZ,JY,JX));   TMGOFS=0._r8
  allocate(TMGCFS(JZ,JY,JX));   TMGCFS=0._r8
  allocate(TMGHFS(JZ,JY,JX));   TMGHFS=0._r8
  allocate(TMGSFS(JZ,JY,JX));   TMGSFS=0._r8
  allocate(TNACFS(JZ,JY,JX));   TNACFS=0._r8
  allocate(TNASFS(JZ,JY,JX));   TNASFS=0._r8
  allocate(TKASFS(JZ,JY,JX));   TKASFS=0._r8
  allocate(TH0PFS(JZ,JY,JX));   TH0PFS=0._r8
  allocate(TH3PFS(JZ,JY,JX));   TH3PFS=0._r8
  allocate(TF1PFS(JZ,JY,JX));   TF1PFS=0._r8
  allocate(TF2PFS(JZ,JY,JX));   TF2PFS=0._r8
  allocate(TC0PFS(JZ,JY,JX));   TC0PFS=0._r8
  allocate(TC1PFS(JZ,JY,JX));   TC1PFS=0._r8
  allocate(TC2PFS(JZ,JY,JX));   TC2PFS=0._r8
  allocate(TM1PFS(JZ,JY,JX));   TM1PFS=0._r8
  allocate(TH0BFB(JZ,JY,JX));   TH0BFB=0._r8
  allocate(TH3BFB(JZ,JY,JX));   TH3BFB=0._r8
  allocate(TF1BFB(JZ,JY,JX));   TF1BFB=0._r8
  allocate(TF2BFB(JZ,JY,JX));   TF2BFB=0._r8
  allocate(TC0BFB(JZ,JY,JX));   TC0BFB=0._r8
  allocate(TC1BFB(JZ,JY,JX));   TC1BFB=0._r8
  allocate(TC2BFB(JZ,JY,JX));   TC2BFB=0._r8
  allocate(TM1BFB(JZ,JY,JX));   TM1BFB=0._r8
  allocate(TALFHS(JZ,JY,JX));   TALFHS=0._r8
  allocate(TFEFHS(JZ,JY,JX));   TFEFHS=0._r8
  allocate(THYFHS(JZ,JY,JX));   THYFHS=0._r8
  allocate(TCAFHS(JZ,JY,JX));   TCAFHS=0._r8
  allocate(TMGFHS(JZ,JY,JX));   TMGFHS=0._r8
  allocate(TNAFHS(JZ,JY,JX));   TNAFHS=0._r8
  allocate(TKAFHS(JZ,JY,JX));   TKAFHS=0._r8
  allocate(TOHFHS(JZ,JY,JX));   TOHFHS=0._r8
  allocate(TSOFHS(JZ,JY,JX));   TSOFHS=0._r8
  allocate(TCLFHS(JZ,JY,JX));   TCLFHS=0._r8
  allocate(TC3FHS(JZ,JY,JX));   TC3FHS=0._r8
  allocate(THCFHS(JZ,JY,JX));   THCFHS=0._r8
  allocate(TAL1HS(JZ,JY,JX));   TAL1HS=0._r8
  allocate(TAL2HS(JZ,JY,JX));   TAL2HS=0._r8
  allocate(TAL3HS(JZ,JY,JX));   TAL3HS=0._r8
  allocate(TAL4HS(JZ,JY,JX));   TAL4HS=0._r8
  allocate(TALSHS(JZ,JY,JX));   TALSHS=0._r8
  allocate(TFE1HS(JZ,JY,JX));   TFE1HS=0._r8
  allocate(TFE2HS(JZ,JY,JX));   TFE2HS=0._r8
  allocate(TFE3HS(JZ,JY,JX));   TFE3HS=0._r8
  allocate(TFE4HS(JZ,JY,JX));   TFE4HS=0._r8
  allocate(TFESHS(JZ,JY,JX));   TFESHS=0._r8
  allocate(TCAOHS(JZ,JY,JX));   TCAOHS=0._r8
  allocate(TCACHS(JZ,JY,JX));   TCACHS=0._r8
  allocate(TCAHHS(JZ,JY,JX));   TCAHHS=0._r8
  allocate(TCASHS(JZ,JY,JX));   TCASHS=0._r8
  allocate(TMGOHS(JZ,JY,JX));   TMGOHS=0._r8
  allocate(TMGCHS(JZ,JY,JX));   TMGCHS=0._r8
  allocate(TMGHHS(JZ,JY,JX));   TMGHHS=0._r8
  allocate(TMGSHS(JZ,JY,JX));   TMGSHS=0._r8
  allocate(TNACHS(JZ,JY,JX));   TNACHS=0._r8
  allocate(TNASHS(JZ,JY,JX));   TNASHS=0._r8
  allocate(TKASHS(JZ,JY,JX));   TKASHS=0._r8
  allocate(TH0PHS(JZ,JY,JX));   TH0PHS=0._r8
  allocate(TH3PHS(JZ,JY,JX));   TH3PHS=0._r8
  allocate(TF1PHS(JZ,JY,JX));   TF1PHS=0._r8
  allocate(TF2PHS(JZ,JY,JX));   TF2PHS=0._r8
  allocate(TC0PHS(JZ,JY,JX));   TC0PHS=0._r8
  allocate(TC1PHS(JZ,JY,JX));   TC1PHS=0._r8
  allocate(TC2PHS(JZ,JY,JX));   TC2PHS=0._r8
  allocate(TM1PHS(JZ,JY,JX));   TM1PHS=0._r8
  allocate(TH0BHB(JZ,JY,JX));   TH0BHB=0._r8
  allocate(TH3BHB(JZ,JY,JX));   TH3BHB=0._r8
  allocate(TF1BHB(JZ,JY,JX));   TF1BHB=0._r8
  allocate(TF2BHB(JZ,JY,JX));   TF2BHB=0._r8
  allocate(TC0BHB(JZ,JY,JX));   TC0BHB=0._r8
  allocate(TC1BHB(JZ,JY,JX));   TC1BHB=0._r8
  allocate(TC2BHB(JZ,JY,JX));   TC2BHB=0._r8
  allocate(TM1BHB(JZ,JY,JX));   TM1BHB=0._r8
  allocate(trcsa_RFLZ(idsa_beg:idsab_end,JZ,JY,JX));   trcsa_RFLZ=0._r8
  end subroutine InitTrnsfrs
!----------------------------------------------------------------------
  SUBROUTINE trnsfrs(I,J,NHW,NHE,NVN,NVS)
!
!     Description:
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     SALT SOLUTES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NX,NY,L,M
!     execution begins here

  IF(ISALTG.EQ.0)RETURN
!
!     TIME STEPS FOR SOLUTE FLUX CALCULATIONS
!
  call SaltModelSoluteFlux(I,NHW,NHE,NVN,NVS)
!
!     TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
!
  D30: DO M=1,NPH
!
    call InitFluxArrays(NHW,NHE,NVN,NVS)

    call SaltModelSoluteHydroFlux(I,M,NHW,NHE,NVN,NVS)
!
!     BOUNDARY SOLUTE AND GAS FLUXES
!
    call SaltModelInternalFlux(M,NHW,NHE,NVN,NVS)
!
!     UPDATE STATE VARIABLES FROM TOTAL FLUXES CALCULATED ABOVE
!
    D9695: DO NX=NHW,NHE
      D9690: DO NY=NVN,NVS
!
!     STATE VARIABLES FOR SOLUTES IN MICROPORES AND MACROPORES IN
!     SOIL SURFACE LAYER FROM OVERLAND FLOW
!
        call UpdateSoluteInSnow(NY,NX)
!
!     STATE VARIABLES FOR SOLUTES IN SURFACE RESIDUE FROM OVERLAND
!     FLOW AND SURFACE FLUX
        call UpdateSoluteInResidue(NY,NX)
!
!     STATE VARIABLES FOR GASES AND FOR SOLUTES IN MICROPORES AND
!     MACROPORES IN SOIL LAYERS FROM SUBSURFACE FLOW, EQUILIBRIUM
!     REACTIONS IN SOLUTE
        call UpdateSoluteInMicMacpores(NY,NX)
      ENDDO D9690
    ENDDO D9695
  ENDDO D30
  RETURN

  END subroutine trnsfrs
!------------------------------------------------------------------------------------------

  subroutine ZeroAtmosSoluteFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: NTSA
!     begin_execution
  DO NTSA=idsa_beg,idsa_end
    trcsa_XBLS(NTSA,1,NY,NX)=0.0
    trcsa_XFLS(NTSA,3,0,NY,NX)=0.0
  ENDDO

  DO NTSA=idsa_beg,idsab_end
    trcsa_XFLS(NTSA,3,NU(NY,NX),NY,NX)=0.0
  ENDDO
  end subroutine ZeroAtmosSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToSnowpack(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: NTSA

!     begin_execution

  trcsa_XBLS(idsa_Al,1,NY,NX)=FLQGQ(NY,NX)*CALR(NY,NX)+FLQGI(NY,NX)*CALQ(I,NY,NX)
  trcsa_XBLS(idsa_Fe,1,NY,NX)=FLQGQ(NY,NX)*CFER(NY,NX)+FLQGI(NY,NX)*CFEQ(I,NY,NX)
  trcsa_XBLS(idsa_Hp,1,NY,NX)=FLQGQ(NY,NX)*CHYR(NY,NX)+FLQGI(NY,NX)*CHYQ(I,NY,NX)
  trcsa_XBLS(idsa_Ca,1,NY,NX)=FLQGQ(NY,NX)*CCAR(NY,NX)+FLQGI(NY,NX)*CCAQ(I,NY,NX)
  trcsa_XBLS(idsa_Mg,1,NY,NX)=FLQGQ(NY,NX)*CMGR(NY,NX)+FLQGI(NY,NX)*CMGQ(I,NY,NX)
  trcsa_XBLS(idsa_Na,1,NY,NX)=FLQGQ(NY,NX)*CNAR(NY,NX)+FLQGI(NY,NX)*CNAQ(I,NY,NX)
  trcsa_XBLS(idsa_K,1,NY,NX)=FLQGQ(NY,NX)*CKAR(NY,NX)+FLQGI(NY,NX)*CKAQ(I,NY,NX)
  trcsa_XBLS(idsa_OH,1,NY,NX)=FLQGQ(NY,NX)*COHR(NY,NX)+FLQGI(NY,NX)*COHQ(I,NY,NX)
  trcsa_XBLS(idsa_SO4,1,NY,NX)=FLQGQ(NY,NX)*CSOR(NY,NX)+FLQGI(NY,NX)*CSOQ(I,NY,NX)
  trcsa_XBLS(idsa_Cl,1,NY,NX)=FLQGQ(NY,NX)*CCLR(NY,NX)+FLQGI(NY,NX)*CCLQ(I,NY,NX)
  trcsa_XBLS(idsa_CO3,1,NY,NX)=FLQGQ(NY,NX)*CC3R(NY,NX)+FLQGI(NY,NX)*CC3Q(I,NY,NX)
  trcsa_XBLS(idsa_HCO3,1,NY,NX)=FLQGQ(NY,NX)*CHCR(NY,NX)+FLQGI(NY,NX)*CHCQ(I,NY,NX)
  trcsa_XBLS(idsa_AlOH,1,NY,NX)=FLQGQ(NY,NX)*CAL1R(NY,NX)+FLQGI(NY,NX)*CAL1Q(I,NY,NX)
  trcsa_XBLS(idsa_AlOH2,1,NY,NX)=FLQGQ(NY,NX)*CAL2R(NY,NX)+FLQGI(NY,NX)*CAL2Q(I,NY,NX)
  trcsa_XBLS(idsa_AlOH3,1,NY,NX)=FLQGQ(NY,NX)*CAL3R(NY,NX)+FLQGI(NY,NX)*CAL3Q(I,NY,NX)
  trcsa_XBLS(idsa_AlOH4,1,NY,NX)=FLQGQ(NY,NX)*CAL4R(NY,NX)+FLQGI(NY,NX)*CAL4Q(I,NY,NX)
  trcsa_XBLS(idsa_AlSO4,1,NY,NX)=FLQGQ(NY,NX)*CALSR(NY,NX)+FLQGI(NY,NX)*CALSQ(I,NY,NX)
  trcsa_XBLS(idsa_FeOH,1,NY,NX)=FLQGQ(NY,NX)*CFE1R(NY,NX)+FLQGI(NY,NX)*CFE1Q(I,NY,NX)
  trcsa_XBLS(idsa_FeOH2,1,NY,NX)=FLQGQ(NY,NX)*CFE2R(NY,NX)+FLQGI(NY,NX)*CFE2Q(I,NY,NX)
  trcsa_XBLS(idsa_FeOH3,1,NY,NX)=FLQGQ(NY,NX)*CFE3R(NY,NX)+FLQGI(NY,NX)*CFE3Q(I,NY,NX)
  trcsa_XBLS(idsa_FeOH4,1,NY,NX)=FLQGQ(NY,NX)*CFE4R(NY,NX)+FLQGI(NY,NX)*CFE4Q(I,NY,NX)
  trcsa_XBLS(idsa_FeSO4,1,NY,NX)=FLQGQ(NY,NX)*CFESR(NY,NX)+FLQGI(NY,NX)*CFESQ(I,NY,NX)
  trcsa_XBLS(idsa_CaOH2,1,NY,NX)=FLQGQ(NY,NX)*CCAOR(NY,NX)+FLQGI(NY,NX)*CCAOQ(I,NY,NX)
  trcsa_XBLS(idsa_CaCO3,1,NY,NX)=FLQGQ(NY,NX)*CCACR(NY,NX)+FLQGI(NY,NX)*CCACQ(I,NY,NX)
  trcsa_XBLS(idsa_CaHCO3,1,NY,NX)=FLQGQ(NY,NX)*CCAHR(NY,NX)+FLQGI(NY,NX)*CCAHQ(I,NY,NX)
  trcsa_XBLS(idsa_CaSO4,1,NY,NX)=FLQGQ(NY,NX)*CCASR(NY,NX)+FLQGI(NY,NX)*CCASQ(I,NY,NX)
  trcsa_XBLS(idsa_MgOH2,1,NY,NX)=FLQGQ(NY,NX)*CMGOR(NY,NX)+FLQGI(NY,NX)*CMGOQ(I,NY,NX)
  trcsa_XBLS(idsa_MgCO3,1,NY,NX)=FLQGQ(NY,NX)*CMGCR(NY,NX)+FLQGI(NY,NX)*CMGCQ(I,NY,NX)
  trcsa_XBLS(idsa_MgHCO3,1,NY,NX)=FLQGQ(NY,NX)*CMGHR(NY,NX)+FLQGI(NY,NX)*CMGHQ(I,NY,NX)
  trcsa_XBLS(idsa_MgSO4,1,NY,NX)=FLQGQ(NY,NX)*CMGSR(NY,NX)+FLQGI(NY,NX)*CMGSQ(I,NY,NX)
  trcsa_XBLS(idsa_NaCO3,1,NY,NX)=FLQGQ(NY,NX)*CNACR(NY,NX)+FLQGI(NY,NX)*CNACQ(I,NY,NX)
  trcsa_XBLS(idsa_NaSO4,1,NY,NX)=FLQGQ(NY,NX)*CNASR(NY,NX)+FLQGI(NY,NX)*CNASQ(I,NY,NX)
  trcsa_XBLS(idsa_KSO4,1,NY,NX)=FLQGQ(NY,NX)*CKASR(NY,NX)+FLQGI(NY,NX)*CKASQ(I,NY,NX)
  trcsa_XBLS(idsa_H0PO4,1,NY,NX)=FLQGQ(NY,NX)*CH0PR(NY,NX)+FLQGI(NY,NX)*CH0PQ(I,NY,NX)
  trcsa_XBLS(idsa_H3PO4,1,NY,NX)=FLQGQ(NY,NX)*CH3PR(NY,NX)+FLQGI(NY,NX)*CH3PQ(I,NY,NX)
  trcsa_XBLS(idsa_FeHPO4,1,NY,NX)=FLQGQ(NY,NX)*CF1PR(NY,NX)+FLQGI(NY,NX)*CF1PQ(I,NY,NX)
  trcsa_XBLS(idsa_FeH2PO4,1,NY,NX)=FLQGQ(NY,NX)*CF2PR(NY,NX)+FLQGI(NY,NX)*CF2PQ(I,NY,NX)
  trcsa_XBLS(idsa_CaPO4,1,NY,NX)=FLQGQ(NY,NX)*CC0PR(NY,NX)+FLQGI(NY,NX)*CC0PQ(I,NY,NX)
  trcsa_XBLS(idsa_CaHPO4,1,NY,NX)=FLQGQ(NY,NX)*CC1PR(NY,NX)+FLQGI(NY,NX)*CC1PQ(I,NY,NX)
  trcsa_XBLS(idsa_CaH2PO4,1,NY,NX)=FLQGQ(NY,NX)*CC2PR(NY,NX)+FLQGI(NY,NX)*CC2PQ(I,NY,NX)
  trcsa_XBLS(idsa_MgHPO4,1,NY,NX)=FLQGQ(NY,NX)*CM1PR(NY,NX)+FLQGI(NY,NX)*CM1PQ(I,NY,NX)
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ARE ZERO IF SNOWPACK IS PRESENT
!
!     X*FLS,X*FLB=hourly solute flux to micropores in non-band,band
!
  DO NTSA=idsa_beg,idsa_end
    trcsa_XFLS(NTSA,3,0,NY,NX)=0.0
  ENDDO

  DO NTSA=idsa_beg,idsab_end
    trcsa_XFLS(NTSA,3,NU(NY,NX),NY,NX)=0.0
  ENDDO
  end subroutine AtmosSoluteFluxToSnowpack
!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToTopsoil(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: NTSA
!     begin_execution
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION IS ZERO IF SNOWPACK IS ABSENT
!
!     PRECQ,PRECI=snow+rain,irrigation
!     X*BLS=hourly solute flux to snowpack
!     X*FLS,X*FLB=hourly solute flux to surface litter,soil surface micropore non-band,band
!     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to soil surface from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO NTSA=idsa_beg,idsa_end
    trcsa_XBLS(NTSA,1,NY,NX)=0.0_r8
  ENDDO

  trcsa_XFLS(idsa_Al,3,0,NY,NX)=FLQRQ(NY,NX)*CALR(NY,NX)+FLQRI(NY,NX)*CALQ(I,NY,NX)
  trcsa_XFLS(idsa_Fe,3,0,NY,NX)=FLQRQ(NY,NX)*CFER(NY,NX)+FLQRI(NY,NX)*CFEQ(I,NY,NX)
  trcsa_XFLS(idsa_Hp,3,0,NY,NX)=FLQRQ(NY,NX)*CHYR(NY,NX)+FLQRI(NY,NX)*CHYQ(I,NY,NX)
  trcsa_XFLS(idsa_Ca,3,0,NY,NX)=FLQRQ(NY,NX)*CCAR(NY,NX)+FLQRI(NY,NX)*CCAQ(I,NY,NX)
  trcsa_XFLS(idsa_Mg,3,0,NY,NX)=FLQRQ(NY,NX)*CMGR(NY,NX)+FLQRI(NY,NX)*CMGQ(I,NY,NX)
  trcsa_XFLS(idsa_Na,3,0,NY,NX)=FLQRQ(NY,NX)*CNAR(NY,NX)+FLQRI(NY,NX)*CNAQ(I,NY,NX)
  trcsa_XFLS(idsa_K,3,0,NY,NX)=FLQRQ(NY,NX)*CKAR(NY,NX)+FLQRI(NY,NX)*CKAQ(I,NY,NX)
  trcsa_XFLS(idsa_OH,3,0,NY,NX)=FLQRQ(NY,NX)*COHR(NY,NX)+FLQRI(NY,NX)*COHQ(I,NY,NX)
  trcsa_XFLS(idsa_SO4,3,0,NY,NX)=FLQRQ(NY,NX)*CSOR(NY,NX)+FLQRI(NY,NX)*CSOQ(I,NY,NX)
  trcsa_XFLS(idsa_Cl,3,0,NY,NX)=FLQRQ(NY,NX)*CCLR(NY,NX)+FLQRI(NY,NX)*CCLQ(I,NY,NX)
  trcsa_XFLS(idsa_CO3,3,0,NY,NX)=FLQRQ(NY,NX)*CC3R(NY,NX)+FLQRI(NY,NX)*CC3Q(I,NY,NX)
  trcsa_XFLS(idsa_HCO3,3,0,NY,NX)=FLQRQ(NY,NX)*CHCR(NY,NX)+FLQRI(NY,NX)*CHCQ(I,NY,NX)
  trcsa_XFLS(idsa_AlOH,3,0,NY,NX)=FLQRQ(NY,NX)*CAL1R(NY,NX)+FLQRI(NY,NX)*CAL1Q(I,NY,NX)
  trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)=FLQRQ(NY,NX)*CAL2R(NY,NX)+FLQRI(NY,NX)*CAL2Q(I,NY,NX)
  trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)=FLQRQ(NY,NX)*CAL3R(NY,NX)+FLQRI(NY,NX)*CAL3Q(I,NY,NX)
  trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)=FLQRQ(NY,NX)*CAL4R(NY,NX)+FLQRI(NY,NX)*CAL4Q(I,NY,NX)
  trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)=FLQRQ(NY,NX)*CALSR(NY,NX)+FLQRI(NY,NX)*CALSQ(I,NY,NX)
  trcsa_XFLS(idsa_FeOH,3,0,NY,NX)=FLQRQ(NY,NX)*CFE1R(NY,NX)+FLQRI(NY,NX)*CFE1Q(I,NY,NX)
  trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)=FLQRQ(NY,NX)*CFE2R(NY,NX)+FLQRI(NY,NX)*CFE2Q(I,NY,NX)
  trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)=FLQRQ(NY,NX)*CFE3R(NY,NX)+FLQRI(NY,NX)*CFE3Q(I,NY,NX)
  trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)=FLQRQ(NY,NX)*CFE4R(NY,NX)+FLQRI(NY,NX)*CFE4Q(I,NY,NX)
  trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)=FLQRQ(NY,NX)*CFESR(NY,NX)+FLQRI(NY,NX)*CFESQ(I,NY,NX)
  trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)=FLQRQ(NY,NX)*CCAOR(NY,NX)+FLQRI(NY,NX)*CCAOQ(I,NY,NX)
  trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)=FLQRQ(NY,NX)*CCACR(NY,NX)+FLQRI(NY,NX)*CCACQ(I,NY,NX)
  trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)=FLQRQ(NY,NX)*CCAHR(NY,NX)+FLQRI(NY,NX)*CCAHQ(I,NY,NX)
  trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)=FLQRQ(NY,NX)*CCASR(NY,NX)+FLQRI(NY,NX)*CCASQ(I,NY,NX)
  trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)=FLQRQ(NY,NX)*CMGOR(NY,NX)+FLQRI(NY,NX)*CMGOQ(I,NY,NX)
  trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)=FLQRQ(NY,NX)*CMGCR(NY,NX)+FLQRI(NY,NX)*CMGCQ(I,NY,NX)
  trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)=FLQRQ(NY,NX)*CMGHR(NY,NX)+FLQRI(NY,NX)*CMGHQ(I,NY,NX)
  trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)=FLQRQ(NY,NX)*CMGSR(NY,NX)+FLQRI(NY,NX)*CMGSQ(I,NY,NX)
  trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)=FLQRQ(NY,NX)*CNACR(NY,NX)+FLQRI(NY,NX)*CNACQ(I,NY,NX)
  trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)=FLQRQ(NY,NX)*CNASR(NY,NX)+FLQRI(NY,NX)*CNASQ(I,NY,NX)
  trcsa_XFLS(idsa_KSO4,3,0,NY,NX)=FLQRQ(NY,NX)*CKASR(NY,NX)+FLQRI(NY,NX)*CKASQ(I,NY,NX)
  trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)=FLQRQ(NY,NX)*CH0PR(NY,NX)+FLQRI(NY,NX)*CH0PQ(I,NY,NX)
  trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)=FLQRQ(NY,NX)*CH3PR(NY,NX)+FLQRI(NY,NX)*CH3PQ(I,NY,NX)
  trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)=FLQRQ(NY,NX)*CF1PR(NY,NX)+FLQRI(NY,NX)*CF1PQ(I,NY,NX)
  trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)=FLQRQ(NY,NX)*CF2PR(NY,NX)+FLQRI(NY,NX)*CF2PQ(I,NY,NX)
  trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)=FLQRQ(NY,NX)*CC0PR(NY,NX)+FLQRI(NY,NX)*CC0PQ(I,NY,NX)
  trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)=FLQRQ(NY,NX)*CC1PR(NY,NX)+FLQRI(NY,NX)*CC1PQ(I,NY,NX)
  trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)=FLQRQ(NY,NX)*CC2PR(NY,NX)+FLQRI(NY,NX)*CC2PQ(I,NY,NX)
  trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)=FLQRQ(NY,NX)*CM1PR(NY,NX)+FLQRI(NY,NX)*CM1PQ(I,NY,NX)
  trcsa_XFLS(idsa_Al,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CALR(NY,NX)+FLQGI(NY,NX)*CALQ(I,NY,NX)
  trcsa_XFLS(idsa_Fe,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFER(NY,NX)+FLQGI(NY,NX)*CFEQ(I,NY,NX)
  trcsa_XFLS(idsa_Hp,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CHYR(NY,NX)+FLQGI(NY,NX)*CHYQ(I,NY,NX)
  trcsa_XFLS(idsa_Ca,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCAR(NY,NX)+FLQGI(NY,NX)*CCAQ(I,NY,NX)
  trcsa_XFLS(idsa_Mg,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGR(NY,NX)+FLQGI(NY,NX)*CMGQ(I,NY,NX)
  trcsa_XFLS(idsa_Na,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNAR(NY,NX)+FLQGI(NY,NX)*CNAQ(I,NY,NX)
  trcsa_XFLS(idsa_K,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CKAR(NY,NX)+FLQGI(NY,NX)*CKAQ(I,NY,NX)
  trcsa_XFLS(idsa_OH,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*COHR(NY,NX)+FLQGI(NY,NX)*COHQ(I,NY,NX)
  trcsa_XFLS(idsa_SO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CSOR(NY,NX)+FLQGI(NY,NX)*CSOQ(I,NY,NX)
  trcsa_XFLS(idsa_Cl,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCLR(NY,NX)+FLQGI(NY,NX)*CCLQ(I,NY,NX)
  trcsa_XFLS(idsa_CO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CC3R(NY,NX)+FLQGI(NY,NX)*CC3Q(I,NY,NX)
  trcsa_XFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CHCR(NY,NX)+FLQGI(NY,NX)*CHCQ(I,NY,NX)
  trcsa_XFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL1R(NY,NX)+FLQGI(NY,NX)*CAL1Q(I,NY,NX)
  trcsa_XFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL2R(NY,NX)+FLQGI(NY,NX)*CAL2Q(I,NY,NX)
  trcsa_XFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL3R(NY,NX)+FLQGI(NY,NX)*CAL3Q(I,NY,NX)
  trcsa_XFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL4R(NY,NX)+FLQGI(NY,NX)*CAL4Q(I,NY,NX)
  trcsa_XFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CALSR(NY,NX)+FLQGI(NY,NX)*CALSQ(I,NY,NX)
  trcsa_XFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE1R(NY,NX)+FLQGI(NY,NX)*CFE1Q(I,NY,NX)
  trcsa_XFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE2R(NY,NX)+FLQGI(NY,NX)*CFE2Q(I,NY,NX)
  trcsa_XFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE3R(NY,NX)+FLQGI(NY,NX)*CFE3Q(I,NY,NX)
  trcsa_XFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE4R(NY,NX)+FLQGI(NY,NX)*CFE4Q(I,NY,NX)
  trcsa_XFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFESR(NY,NX)+FLQGI(NY,NX)*CFESQ(I,NY,NX)
  trcsa_XFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCAOR(NY,NX)+FLQGI(NY,NX)*CCAOQ(I,NY,NX)
  trcsa_XFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCACR(NY,NX)+FLQGI(NY,NX)*CCACQ(I,NY,NX)
  trcsa_XFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCAHR(NY,NX)+FLQGI(NY,NX)*CCAHQ(I,NY,NX)
  trcsa_XFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCASR(NY,NX)+FLQGI(NY,NX)*CCASQ(I,NY,NX)
  trcsa_XFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGOR(NY,NX)+FLQGI(NY,NX)*CMGOQ(I,NY,NX)
  trcsa_XFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGCR(NY,NX)+FLQGI(NY,NX)*CMGCQ(I,NY,NX)
  trcsa_XFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGHR(NY,NX)+FLQGI(NY,NX)*CMGHQ(I,NY,NX)
  trcsa_XFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGSR(NY,NX)+FLQGI(NY,NX)*CMGSQ(I,NY,NX)
  trcsa_XFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNACR(NY,NX)+FLQGI(NY,NX)*CNACQ(I,NY,NX)
  trcsa_XFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNASR(NY,NX)+FLQGI(NY,NX)*CNASQ(I,NY,NX)
  trcsa_XFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CKASR(NY,NX)+FLQGI(NY,NX)*CKASQ(I,NY,NX)
  trcsa_XFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH0PR(NY,NX)+FLQGI(NY,NX)*CH0PQ(I,NY,NX)) &
    *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH3PR(NY,NX)+FLQGI(NY,NX) &
    *CH3PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF1PR(NY,NX)+FLQGI(NY,NX) &
    *CF1PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF2PR(NY,NX)+FLQGI(NY,NX) &
    *CF2PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC0PR(NY,NX)+FLQGI(NY,NX) &
    *CC0PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC1PR(NY,NX)+FLQGI(NY,NX) &
    *CC1PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC2PR(NY,NX)+FLQGI(NY,NX) &
    *CC2PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CM1PR(NY,NX)+FLQGI(NY,NX) &
    *CM1PQ(I,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_H0PO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH0PR(NY,NX)+FLQGI(NY,NX) &
    *CH0PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_H3PO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH3PR(NY,NX)+FLQGI(NY,NX) &
    *CH3PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_FeHPO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF1PR(NY,NX)+FLQGI(NY,NX) &
    *CF1PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_FeH2PO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF2PR(NY,NX)+FLQGI(NY,NX) &
    *CF2PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_CaPO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC0PR(NY,NX)+FLQGI(NY,NX) &
    *CC0PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_CaHPO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC1PR(NY,NX)+FLQGI(NY,NX) &
    *CC1PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_CaH2PO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC2PR(NY,NX)+FLQGI(NY,NX) &
    *CC2PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  trcsa_XFLS(idsa_MgHPO4B,3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CM1PR(NY,NX)+FLQGI(NY,NX) &
    *CM1PQ(I,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  end subroutine AtmosSoluteFluxToTopsoil
!------------------------------------------------------------------------------------------

  subroutine InitSolutesInSnowpack(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,NTA
!     begin_execution
!
!     Z*W=solute content in snowpacl
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-C
  D20: DO L=1,JS
    DO NTA=idsa_beg,idsa_end
      trcs_solsml2(NTA,L,NY,NX)=trcs_solsml(NTA,L,NY,NX)
    ENDDO
  ENDDO D20
  end subroutine InitSolutesInSnowpack
!------------------------------------------------------------------------------------------

  subroutine GetSubHourFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
!
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     R*BLS,R*FL0,R*FL1,R*FL2=solute flux to snowpack,surface litter,soil surface non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  RALBLS(1,NY,NX)=trcsa_XBLS(idsa_Al,1,NY,NX)*XNPH
  RFEBLS(1,NY,NX)=trcsa_XBLS(idsa_Fe,1,NY,NX)*XNPH
  RHYBLS(1,NY,NX)=trcsa_XBLS(idsa_Hp,1,NY,NX)*XNPH
  RCABLS(1,NY,NX)=trcsa_XBLS(idsa_Ca,1,NY,NX)*XNPH
  RMGBLS(1,NY,NX)=trcsa_XBLS(idsa_Mg,1,NY,NX)*XNPH
  RNABLS(1,NY,NX)=trcsa_XBLS(idsa_Na,1,NY,NX)*XNPH
  RKABLS(1,NY,NX)=trcsa_XBLS(idsa_K,1,NY,NX)*XNPH
  ROHBLS(1,NY,NX)=trcsa_XBLS(idsa_OH,1,NY,NX)*XNPH
  RSOBLS(1,NY,NX)=trcsa_XBLS(idsa_SO4,1,NY,NX)*XNPH
  RCLBLS(1,NY,NX)=trcsa_XBLS(idsa_Cl,1,NY,NX)*XNPH
  RC3BLS(1,NY,NX)=trcsa_XBLS(idsa_CO3,1,NY,NX)*XNPH
  RHCBLS(1,NY,NX)=trcsa_XBLS(idsa_HCO3,1,NY,NX)*XNPH
  RAL1BS(1,NY,NX)=trcsa_XBLS(idsa_AlOH,1,NY,NX)*XNPH
  RAL2BS(1,NY,NX)=trcsa_XBLS(idsa_AlOH2,1,NY,NX)*XNPH
  RAL3BS(1,NY,NX)=trcsa_XBLS(idsa_AlOH3,1,NY,NX)*XNPH
  RAL4BS(1,NY,NX)=trcsa_XBLS(idsa_AlOH4,1,NY,NX)*XNPH
  RALSBS(1,NY,NX)=trcsa_XBLS(idsa_AlSO4,1,NY,NX)*XNPH
  RFE1BS(1,NY,NX)=trcsa_XBLS(idsa_FeOH,1,NY,NX)*XNPH
  RFE2BS(1,NY,NX)=trcsa_XBLS(idsa_FeOH2,1,NY,NX)*XNPH
  RFE3BS(1,NY,NX)=trcsa_XBLS(idsa_FeOH3,1,NY,NX)*XNPH
  RFE4BS(1,NY,NX)=trcsa_XBLS(idsa_FeOH4,1,NY,NX)*XNPH
  RFESBS(1,NY,NX)=trcsa_XBLS(idsa_FeSO4,1,NY,NX)*XNPH
  RCAOBS(1,NY,NX)=trcsa_XBLS(idsa_CaOH2,1,NY,NX)*XNPH
  RCACBS(1,NY,NX)=trcsa_XBLS(idsa_CaCO3,1,NY,NX)*XNPH
  RCAHBS(1,NY,NX)=trcsa_XBLS(idsa_CaHCO3,1,NY,NX)*XNPH
  RCASBS(1,NY,NX)=trcsa_XBLS(idsa_CaSO4,1,NY,NX)*XNPH
  RMGOBS(1,NY,NX)=trcsa_XBLS(idsa_MgOH2,1,NY,NX)*XNPH
  RMGCBS(1,NY,NX)=trcsa_XBLS(idsa_MgCO3,1,NY,NX)*XNPH
  RMGHBS(1,NY,NX)=trcsa_XBLS(idsa_MgHCO3,1,NY,NX)*XNPH
  RMGSBS(1,NY,NX)=trcsa_XBLS(idsa_MgSO4,1,NY,NX)*XNPH
  RNACBS(1,NY,NX)=trcsa_XBLS(idsa_NaCO3,1,NY,NX)*XNPH
  RNASBS(1,NY,NX)=trcsa_XBLS(idsa_NaSO4,1,NY,NX)*XNPH
  RKASBS(1,NY,NX)=trcsa_XBLS(idsa_KSO4,1,NY,NX)*XNPH
  RH0PBS(1,NY,NX)=trcsa_XBLS(idsa_H0PO4,1,NY,NX)*XNPH
  RH3PBS(1,NY,NX)=trcsa_XBLS(idsa_H3PO4,1,NY,NX)*XNPH
  RF1PBS(1,NY,NX)=trcsa_XBLS(idsa_FeHPO4,1,NY,NX)*XNPH
  RF2PBS(1,NY,NX)=trcsa_XBLS(idsa_FeH2PO4,1,NY,NX)*XNPH
  RC0PBS(1,NY,NX)=trcsa_XBLS(idsa_CaPO4,1,NY,NX)*XNPH
  RC1PBS(1,NY,NX)=trcsa_XBLS(idsa_CaHPO4,1,NY,NX)*XNPH
  RC2PBS(1,NY,NX)=trcsa_XBLS(idsa_CaH2PO4,1,NY,NX)*XNPH
  RM1PBS(1,NY,NX)=trcsa_XBLS(idsa_MgHPO4,1,NY,NX)*XNPH
  RALFL0(NY,NX)=trcsa_XFLS(idsa_Al,3,0,NY,NX)*XNPH
  RFEFL0(NY,NX)=trcsa_XFLS(idsa_Fe,3,0,NY,NX)*XNPH
  RHYFL0(NY,NX)=trcsa_XFLS(idsa_Hp,3,0,NY,NX)*XNPH
  RCAFL0(NY,NX)=trcsa_XFLS(idsa_Ca,3,0,NY,NX)*XNPH
  RMGFL0(NY,NX)=trcsa_XFLS(idsa_Mg,3,0,NY,NX)*XNPH
  RNAFL0(NY,NX)=trcsa_XFLS(idsa_Na,3,0,NY,NX)*XNPH
  RKAFL0(NY,NX)=trcsa_XFLS(idsa_K,3,0,NY,NX)*XNPH
  ROHFL0(NY,NX)=trcsa_XFLS(idsa_OH,3,0,NY,NX)*XNPH
  RSOFL0(NY,NX)=trcsa_XFLS(idsa_SO4,3,0,NY,NX)*XNPH
  RCLFL0(NY,NX)=trcsa_XFLS(idsa_Cl,3,0,NY,NX)*XNPH
  RC3FL0(NY,NX)=trcsa_XFLS(idsa_CO3,3,0,NY,NX)*XNPH
  RHCFL0(NY,NX)=trcsa_XFLS(idsa_HCO3,3,0,NY,NX)*XNPH
  RAL1F0(NY,NX)=trcsa_XFLS(idsa_AlOH,3,0,NY,NX)*XNPH
  RAL2F0(NY,NX)=trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)*XNPH
  RAL3F0(NY,NX)=trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)*XNPH
  RAL4F0(NY,NX)=trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)*XNPH
  RALSF0(NY,NX)=trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)*XNPH
  RFE1F0(NY,NX)=trcsa_XFLS(idsa_FeOH,3,0,NY,NX)*XNPH
  RFE2F0(NY,NX)=trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)*XNPH
  RFE3F0(NY,NX)=trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)*XNPH
  RFE4F0(NY,NX)=trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)*XNPH
  RFESF0(NY,NX)=trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)*XNPH
  RCAOF0(NY,NX)=trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)*XNPH
  RCACF0(NY,NX)=trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)*XNPH
  RCAHF0(NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)*XNPH
  RCASF0(NY,NX)=trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)*XNPH
  RMGOF0(NY,NX)=trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)*XNPH
  RMGCF0(NY,NX)=trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)*XNPH
  RMGHF0(NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)*XNPH
  RMGSF0(NY,NX)=trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)*XNPH
  RNACF0(NY,NX)=trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)*XNPH
  RNASF0(NY,NX)=trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)*XNPH
  RKASF0(NY,NX)=trcsa_XFLS(idsa_KSO4,3,0,NY,NX)*XNPH
  RH0PF0(NY,NX)=trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)*XNPH
  RH3PF0(NY,NX)=trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)*XNPH
  RF1PF0(NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)*XNPH
  RF2PF0(NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)*XNPH
  RC0PF0(NY,NX)=trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)*XNPH
  RC1PF0(NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)*XNPH
  RC2PF0(NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)*XNPH
  RM1PF0(NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)*XNPH

  RALFL1(NY,NX)=trcsa_XFLS(idsa_Al,3,NU(NY,NX),NY,NX)*XNPH
  RFEFL1(NY,NX)=trcsa_XFLS(idsa_Fe,3,NU(NY,NX),NY,NX)*XNPH
  RHYFL1(NY,NX)=trcsa_XFLS(idsa_Hp,3,NU(NY,NX),NY,NX)*XNPH
  RCAFL1(NY,NX)=trcsa_XFLS(idsa_Ca,3,NU(NY,NX),NY,NX)*XNPH
  RMGFL1(NY,NX)=trcsa_XFLS(idsa_Mg,3,NU(NY,NX),NY,NX)*XNPH
  RNAFL1(NY,NX)=trcsa_XFLS(idsa_Na,3,NU(NY,NX),NY,NX)*XNPH
  RKAFL1(NY,NX)=trcsa_XFLS(idsa_K,3,NU(NY,NX),NY,NX)*XNPH
  ROHFL1(NY,NX)=trcsa_XFLS(idsa_OH,3,NU(NY,NX),NY,NX)*XNPH
  RSOFL1(NY,NX)=trcsa_XFLS(idsa_SO4,3,NU(NY,NX),NY,NX)*XNPH
  RCLFL1(NY,NX)=trcsa_XFLS(idsa_Cl,3,NU(NY,NX),NY,NX)*XNPH
  RC3FL1(NY,NX)=trcsa_XFLS(idsa_CO3,3,NU(NY,NX),NY,NX)*XNPH
  RHCFL1(NY,NX)=trcsa_XFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)*XNPH
  RAL1F1(NY,NX)=trcsa_XFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)*XNPH
  RAL2F1(NY,NX)=trcsa_XFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)*XNPH
  RAL3F1(NY,NX)=trcsa_XFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)*XNPH
  RAL4F1(NY,NX)=trcsa_XFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)*XNPH
  RALSF1(NY,NX)=trcsa_XFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)*XNPH
  RFE1F1(NY,NX)=trcsa_XFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)*XNPH
  RFE2F1(NY,NX)=trcsa_XFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)*XNPH
  RFE3F1(NY,NX)=trcsa_XFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)*XNPH
  RFE4F1(NY,NX)=trcsa_XFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)*XNPH
  RFESF1(NY,NX)=trcsa_XFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)*XNPH
  RCAOF1(NY,NX)=trcsa_XFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)*XNPH
  RCACF1(NY,NX)=trcsa_XFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)*XNPH
  RCAHF1(NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)*XNPH
  RCASF1(NY,NX)=trcsa_XFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)*XNPH
  RMGOF1(NY,NX)=trcsa_XFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)*XNPH
  RMGCF1(NY,NX)=trcsa_XFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)*XNPH
  RMGHF1(NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)*XNPH
  RMGSF1(NY,NX)=trcsa_XFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)*XNPH
  RNACF1(NY,NX)=trcsa_XFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)*XNPH
  RNASF1(NY,NX)=trcsa_XFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)*XNPH
  RKASF1(NY,NX)=trcsa_XFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)*XNPH
  RH0PF1(NY,NX)=trcsa_XFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)*XNPH
  RH3PF1(NY,NX)=trcsa_XFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)*XNPH
  RF1PF1(NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)*XNPH
  RF2PF1(NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)*XNPH
  RC0PF1(NY,NX)=trcsa_XFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)*XNPH
  RC1PF1(NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)*XNPH
  RC2PF1(NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)*XNPH
  RM1PF1(NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)*XNPH
  RH0BF2(NY,NX)=trcsa_XFLS(idsa_H0PO4B,3,NU(NY,NX),NY,NX)*XNPH
  RH3BF2(NY,NX)=trcsa_XFLS(idsa_H3PO4B,3,NU(NY,NX),NY,NX)*XNPH
  RF1BF2(NY,NX)=trcsa_XFLS(idsa_FeHPO4B,3,NU(NY,NX),NY,NX)*XNPH
  RF2BF2(NY,NX)=trcsa_XFLS(idsa_FeH2PO4B,3,NU(NY,NX),NY,NX)*XNPH
  RC0BF2(NY,NX)=trcsa_XFLS(idsa_CaPO4B,3,NU(NY,NX),NY,NX)*XNPH
  RC1BF2(NY,NX)=trcsa_XFLS(idsa_CaHPO4B,3,NU(NY,NX),NY,NX)*XNPH
  RC2BF2(NY,NX)=trcsa_XFLS(idsa_CaH2PO4B,3,NU(NY,NX),NY,NX)*XNPH
  RM1BF2(NY,NX)=trcsa_XFLS(idsa_MgHPO4B,3,NU(NY,NX),NY,NX)*XNPH
  end subroutine GetSubHourFlux
!------------------------------------------------------------------------------------------

  subroutine GetSubHourlyFluxByLayer(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: L,NTSA
  real(r8) :: FLWU(JZ,JY,JX)


!     begin_execution
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     RZ*2=solute flux at time step for flux calculations
!     TR*=solute transformations from solute.f
!
  D10: DO L=NU(NY,NX),NL(NY,NX)
    trcsa_solml2R(idsa_Al,L,NY,NX)=-trcsa_TR(idsa_Al,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_Fe,L,NY,NX)=-trcsa_TR(idsa_Fe,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_Hp,L,NY,NX)=-(trcsa_TR(idsa_Hp,L,NY,NX)+XZHYS(L,NY,NX))*XNPH
    trcsa_solml2R(idsa_Ca,L,NY,NX)=-trcsa_TR(idsa_Ca,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_Mg,L,NY,NX)=-trcsa_TR(idsa_Mg,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_Na,L,NY,NX)=-trcsa_TR(idsa_Na,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_K,L,NY,NX)=-trcsa_TR(idsa_K,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_OH,L,NY,NX)=-trcsa_TR(idsa_OH,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_SO4,L,NY,NX)=-trcsa_TR(idsa_SO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_Cl,L,NY,NX)=0.0_r8
    trcsa_solml2R(idsa_CO3,L,NY,NX)=-trcsa_TR(idsa_CO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_HCO3,L,NY,NX)=-trcsa_TR(idsa_HCO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_AlOH,L,NY,NX)=-trcsa_TR(idsa_AlOH,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_AlOH2,L,NY,NX)=-trcsa_TR(idsa_AlOH2,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_AlOH3,L,NY,NX)=-trcsa_TR(idsa_AlOH3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_AlOH4,L,NY,NX)=-trcsa_TR(idsa_AlOH4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_AlSO4,L,NY,NX)=-trcsa_TR(idsa_AlSO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeOH,L,NY,NX)=-trcsa_TR(idsa_FeOH,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeOH2,L,NY,NX)=-trcsa_TR(idsa_FeOH2,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeOH3,L,NY,NX)=-trcsa_TR(idsa_FeOH3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeOH4,L,NY,NX)=-trcsa_TR(idsa_FeOH4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeSO4,L,NY,NX)=-trcsa_TR(idsa_FeSO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaOH2,L,NY,NX)=-trcsa_TR(idsa_CaOH2,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaCO3,L,NY,NX)=-trcsa_TR(idsa_CaCO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaHCO3,L,NY,NX)=-trcsa_TR(idsa_CaHCO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaSO4,L,NY,NX)=-trcsa_TR(idsa_CaSO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_MgOH2,L,NY,NX)=-trcsa_TR(idsa_MgOH2,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_MgCO3,L,NY,NX)=-trcsa_TR(idsa_MgCO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_MgHCO3,L,NY,NX)=-trcsa_TR(idsa_MgHCO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_MgSO4,L,NY,NX)=-trcsa_TR(idsa_MgSO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_NaCO3,L,NY,NX)=-trcsa_TR(idsa_NaCO3,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_NaSO4,L,NY,NX)=-trcsa_TR(idsa_NaSO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_KSO4,L,NY,NX)=-trcsa_TR(idsa_KSO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_H0PO4,L,NY,NX)=-trcsa_TR(idsa_H0PO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_H3PO4,L,NY,NX)=-trcsa_TR(idsa_H3PO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeHPO4,L,NY,NX)=-trcsa_TR(idsa_FeHPO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeH2PO4,L,NY,NX)=-trcsa_TR(idsa_FeH2PO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaPO4,L,NY,NX)=-trcsa_TR(idsa_CaPO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaHPO4,L,NY,NX)=-trcsa_TR(idsa_CaHPO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaH2PO4,L,NY,NX)=-trcsa_TR(idsa_CaH2PO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_MgHPO4,L,NY,NX)=-trcsa_TR(idsa_MgHPO4,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_H0PO4B,L,NY,NX)=-trcsa_TR(idsa_H0PO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_H3PO4B,L,NY,NX)=-trcsa_TR(idsa_H3PO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeHPO4B,L,NY,NX)=-trcsa_TR(idsa_FeHPO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_FeH2PO4B,L,NY,NX)=-trcsa_TR(idsa_FeH2PO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaPO4B,L,NY,NX)=-trcsa_TR(idsa_CaPO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaHPO4B,L,NY,NX)=-trcsa_TR(idsa_CaHPO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_CaH2PO4B,L,NY,NX)=-trcsa_TR(idsa_CaH2PO4B,L,NY,NX)*XNPH
    trcsa_solml2R(idsa_MgHPO4B,L,NY,NX)=-trcsa_TR(idsa_MgHPO4B,L,NY,NX)*XNPH
!
!     SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     FLU=subsurface water flux from watsub.f
!     R*FLU,R*FBU=subsurface solute flux in non-band,band
!     C*Q=irrigation solute concentrations
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    FLWU(L,NY,NX)=TUPWTR(L,NY,NX)*XNPH
    trcsa_RFLU(idsa_Al,L,NY,NX)=FLU(L,NY,NX)*CALQ(I,NY,NX)
    trcsa_RFLU(idsa_Fe,L,NY,NX)=FLU(L,NY,NX)*CFEQ(I,NY,NX)
    trcsa_RFLU(idsa_Hp,L,NY,NX)=FLU(L,NY,NX)*CHYQ(I,NY,NX)
    trcsa_RFLU(idsa_Ca,L,NY,NX)=FLU(L,NY,NX)*CCAQ(I,NY,NX)
    trcsa_RFLU(idsa_Mg,L,NY,NX)=FLU(L,NY,NX)*CMGQ(I,NY,NX)
    trcsa_RFLU(idsa_Na,L,NY,NX)=FLU(L,NY,NX)*CNAQ(I,NY,NX)
    trcsa_RFLU(idsa_K,L,NY,NX)=FLU(L,NY,NX)*CKAQ(I,NY,NX)
    trcsa_RFLU(idsa_OH,L,NY,NX)=FLU(L,NY,NX)*COHQ(I,NY,NX)
    trcsa_RFLU(idsa_SO4,L,NY,NX)=FLU(L,NY,NX)*CSOQ(I,NY,NX)
    trcsa_RFLU(idsa_Cl,L,NY,NX)=FLU(L,NY,NX)*CCLQ(I,NY,NX)
    trcsa_RFLU(idsa_CO3,L,NY,NX)=FLU(L,NY,NX)*CC3Q(I,NY,NX)
    trcsa_RFLU(idsa_HCO3,L,NY,NX)=FLU(L,NY,NX)*CHCQ(I,NY,NX)
    trcsa_RFLU(idsa_AlOH,L,NY,NX)=FLU(L,NY,NX)*CAL1Q(I,NY,NX)
    trcsa_RFLU(idsa_AlOH2,L,NY,NX)=FLU(L,NY,NX)*CAL2Q(I,NY,NX)
    trcsa_RFLU(idsa_AlOH3,L,NY,NX)=FLU(L,NY,NX)*CAL3Q(I,NY,NX)
    trcsa_RFLU(idsa_AlOH4,L,NY,NX)=FLU(L,NY,NX)*CAL4Q(I,NY,NX)
    trcsa_RFLU(idsa_AlSO4,L,NY,NX)=FLU(L,NY,NX)*CALSQ(I,NY,NX)
    trcsa_RFLU(idsa_FeOH,L,NY,NX)=FLU(L,NY,NX)*CFE1Q(I,NY,NX)
    trcsa_RFLU(idsa_FeOH2,L,NY,NX)=FLU(L,NY,NX)*CFE2Q(I,NY,NX)
    trcsa_RFLU(idsa_FeOH3,L,NY,NX)=FLU(L,NY,NX)*CFE3Q(I,NY,NX)
    trcsa_RFLU(idsa_FeOH4,L,NY,NX)=FLU(L,NY,NX)*CFE4Q(I,NY,NX)
    trcsa_RFLU(idsa_FeSO4,L,NY,NX)=FLU(L,NY,NX)*CFESQ(I,NY,NX)
    trcsa_RFLU(idsa_CaOH2,L,NY,NX)=FLU(L,NY,NX)*CCAOQ(I,NY,NX)
    trcsa_RFLU(idsa_CaCO3,L,NY,NX)=FLU(L,NY,NX)*CCACQ(I,NY,NX)
    trcsa_RFLU(idsa_CaHCO3,L,NY,NX)=FLU(L,NY,NX)*CCAHQ(I,NY,NX)
    trcsa_RFLU(idsa_CaSO4,L,NY,NX)=FLU(L,NY,NX)*CCASQ(I,NY,NX)
    trcsa_RFLU(idsa_MgOH2,L,NY,NX)=FLU(L,NY,NX)*CMGOQ(I,NY,NX)
    trcsa_RFLU(idsa_MgCO3,L,NY,NX)=FLU(L,NY,NX)*CMGCQ(I,NY,NX)
    trcsa_RFLU(idsa_MgHCO3,L,NY,NX)=FLU(L,NY,NX)*CMGHQ(I,NY,NX)
    trcsa_RFLU(idsa_MgSO4,L,NY,NX)=FLU(L,NY,NX)*CMGSQ(I,NY,NX)
    trcsa_RFLU(idsa_NaCO3,L,NY,NX)=FLU(L,NY,NX)*CNACQ(I,NY,NX)
    trcsa_RFLU(idsa_NaSO4,L,NY,NX)=FLU(L,NY,NX)*CNASQ(I,NY,NX)
    trcsa_RFLU(idsa_KSO4,L,NY,NX)=FLU(L,NY,NX)*CKASQ(I,NY,NX)
    trcsa_RFLU(idsa_H0PO4,L,NY,NX)=FLU(L,NY,NX)*CH0PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_H3PO4,L,NY,NX)=FLU(L,NY,NX)*CH3PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_FeHPO4,L,NY,NX)=FLU(L,NY,NX)*CF1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_FeH2PO4,L,NY,NX)=FLU(L,NY,NX)*CF2PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_CaPO4,L,NY,NX)=FLU(L,NY,NX)*CC0PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_CaHPO4,L,NY,NX)=FLU(L,NY,NX)*CC1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_CaH2PO4,L,NY,NX)=FLU(L,NY,NX)*CC2PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_MgHPO4,L,NY,NX)=FLU(L,NY,NX)*CM1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_RFLU(idsa_H0PO4B,L,NY,NX)=FLU(L,NY,NX)*CH0PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_H3PO4B,L,NY,NX)=FLU(L,NY,NX)*CH3PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_FeHPO4B,L,NY,NX)=FLU(L,NY,NX)*CF1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_FeH2PO4B,L,NY,NX)=FLU(L,NY,NX)*CF2PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_CaPO4B,L,NY,NX)=FLU(L,NY,NX)*CC0PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_CaHPO4B,L,NY,NX)=FLU(L,NY,NX)*CC1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_CaH2PO4B,L,NY,NX)=FLU(L,NY,NX)*CC2PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_RFLU(idsa_MgHPO4B,L,NY,NX)=FLU(L,NY,NX)*CM1PQ(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!
    DO NTSA=idsa_beg,idsab_end
      trcsa_RFLZ(NTSA,L,NY,NX)=trcsa_RFLU(NTSA,L,NY,NX)*XNPH
    ENDDO
!
!     SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:PO=PO4,AL=Al,FE=Fe,HY=H,CA=Ca,GM=Mg,AN=Na,AK=KOH=OH
!                :SO=SO4,CL=Cl,C3=CO3,HC=HCO3
    POSGL2(L,NY,NX)=SolDifc(ids_H1PO4,L,NY,NX)*XNPH
    ALSGL2(L,NY,NX)=ALSGL(L,NY,NX)*XNPH
    FESGL2(L,NY,NX)=FESGL(L,NY,NX)*XNPH
    HYSGL2(L,NY,NX)=HYSGL(L,NY,NX)*XNPH
    CASGL2(L,NY,NX)=CASGL(L,NY,NX)*XNPH
    GMSGL2(L,NY,NX)=GMSGL(L,NY,NX)*XNPH
    ANSGL2(L,NY,NX)=ANSGL(L,NY,NX)*XNPH
    AKSGL2(L,NY,NX)=AKSGL(L,NY,NX)*XNPH
    OHSGL2(L,NY,NX)=OHSGL(L,NY,NX)*XNPH
    SOSGL2(L,NY,NX)=SOSGL(L,NY,NX)*XNPH
    CLSXL2(L,NY,NX)=CLSXL(L,NY,NX)*XNPH
    C3SGL2(L,NY,NX)=C3SGL(L,NY,NX)*XNPH
    HCSGL2(L,NY,NX)=HCSGL(L,NY,NX)*XNPH
!
!     STATE VARIABLES FOR SOLUTES USED IN 'TRNSFRS'
!     TO STORE SUB-HOURLY CHANGES DURING FLUX CALCULATIONS
!     INCLUDING TRANSFORMATIONS FROM REACTIONS IN SOLUTE.F
!
!     Z*,Z*2=soil solute contents
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
    trcsa_solml2(idsa_Al,L,NY,NX)=trcsa_solml(idsa_Al,L,NY,NX)
    trcsa_solml2(idsa_Fe,L,NY,NX)=trcsa_solml(idsa_Fe,L,NY,NX)
    trcsa_solml2(idsa_Hp,L,NY,NX)=trcsa_solml(idsa_Hp,L,NY,NX)
    trcsa_solml2(idsa_Ca,L,NY,NX)=trcsa_solml(idsa_Ca,L,NY,NX)
    trcsa_solml2(idsa_Mg,L,NY,NX)=trcsa_solml(idsa_Mg,L,NY,NX)
    trcsa_solml2(idsa_Na,L,NY,NX)=trcsa_solml(idsa_Na,L,NY,NX)
    trcsa_solml2(idsa_K,L,NY,NX)=trcsa_solml(idsa_K,L,NY,NX)
    trcsa_solml2(idsa_OH,L,NY,NX)=trcsa_solml(idsa_OH,L,NY,NX)
    trcsa_solml2(idsa_SO4,L,NY,NX)=trcsa_solml(idsa_SO4,L,NY,NX)
    trcsa_solml2(idsa_Cl,L,NY,NX)=trcsa_solml(idsa_Cl,L,NY,NX)
    trcsa_solml2(idsa_CO3,L,NY,NX)=trcsa_solml(idsa_CO3,L,NY,NX)
    trcsa_solml2(idsa_HCO3,L,NY,NX)=trcsa_solml(idsa_HCO3,L,NY,NX)
    trcsa_solml2(idsa_AlOH,L,NY,NX)=trcsa_solml(idsa_AlOH,L,NY,NX)
    trcsa_solml2(idsa_AlOH2,L,NY,NX)=trcsa_solml(idsa_AlOH2,L,NY,NX)
    trcsa_solml2(idsa_AlOH3,L,NY,NX)=trcsa_solml(idsa_AlOH3,L,NY,NX)
    trcsa_solml2(idsa_AlOH4,L,NY,NX)=trcsa_solml(idsa_AlOH4,L,NY,NX)
    trcsa_solml2(idsa_AlSO4,L,NY,NX)=trcsa_solml(idsa_AlSO4,L,NY,NX)
    trcsa_solml2(idsa_FeOH,L,NY,NX)=trcsa_solml(idsa_FeOH,L,NY,NX)
    trcsa_solml2(idsa_FeOH2,L,NY,NX)=trcsa_solml(idsa_FeOH2,L,NY,NX)
    trcsa_solml2(idsa_FeOH3,L,NY,NX)=trcsa_solml(idsa_FeOH3,L,NY,NX)
    trcsa_solml2(idsa_FeOH4,L,NY,NX)=trcsa_solml(idsa_FeOH4,L,NY,NX)
    trcsa_solml2(idsa_FeSO4,L,NY,NX)=trcsa_solml(idsa_FeSO4,L,NY,NX)
    trcsa_solml2(idsa_CaOH2,L,NY,NX)=trcsa_solml(idsa_CaOH2,L,NY,NX)
    trcsa_solml2(idsa_CaCO3,L,NY,NX)=trcsa_solml(idsa_CaCO3,L,NY,NX)
    trcsa_solml2(idsa_CaHCO3,L,NY,NX)=trcsa_solml(idsa_CaHCO3,L,NY,NX)
    trcsa_solml2(idsa_CaSO4,L,NY,NX)=trcsa_solml(idsa_CaSO4,L,NY,NX)
    trcsa_solml2(idsa_MgOH2,L,NY,NX)=trcsa_solml(idsa_MgOH2,L,NY,NX)
    trcsa_solml2(idsa_MgCO3,L,NY,NX)=trcsa_solml(idsa_MgCO3,L,NY,NX)
    trcsa_solml2(idsa_MgHCO3,L,NY,NX)=trcsa_solml(idsa_MgHCO3,L,NY,NX)
    trcsa_solml2(idsa_MgSO4,L,NY,NX)=trcsa_solml(idsa_MgSO4,L,NY,NX)
    trcsa_solml2(idsa_NaCO3,L,NY,NX)=trcsa_solml(idsa_NaCO3,L,NY,NX)
    trcsa_solml2(idsa_NaSO4,L,NY,NX)=trcsa_solml(idsa_NaSO4,L,NY,NX)
    trcsa_solml2(idsa_KSO4,L,NY,NX)=trcsa_solml(idsa_KSO4,L,NY,NX)
    trcsa_solml2(idsa_H0PO4,L,NY,NX)=trcsa_solml(idsa_H0PO4,L,NY,NX)
    trcsa_solml2(idsa_H3PO4,L,NY,NX)=trcsa_solml(idsa_H3PO4,L,NY,NX)
    trcsa_solml2(idsa_FeHPO4,L,NY,NX)=trcsa_solml(idsa_FeHPO4,L,NY,NX)
    trcsa_solml2(idsa_FeH2PO4,L,NY,NX)=trcsa_solml(idsa_FeH2PO4,L,NY,NX)
    trcsa_solml2(idsa_CaPO4,L,NY,NX)=trcsa_solml(idsa_CaPO4,L,NY,NX)
    trcsa_solml2(idsa_CaHPO4,L,NY,NX)=trcsa_solml(idsa_CaHPO4,L,NY,NX)
    trcsa_solml2(idsa_CaH2PO4,L,NY,NX)=trcsa_solml(idsa_CaH2PO4,L,NY,NX)
    trcsa_solml2(idsa_MgHPO4,L,NY,NX)=trcsa_solml(idsa_MgHPO4,L,NY,NX)
    trcsa_solml2(idsa_H0PO4B,L,NY,NX)=trcsa_solml(idsa_H0PO4B,L,NY,NX)
    trcsa_solml2(idsa_H3PO4B,L,NY,NX)=trcsa_solml(idsa_H3PO4B,L,NY,NX)
    trcsa_solml2(idsa_FeHPO4B,L,NY,NX)=trcsa_solml(idsa_FeHPO4B,L,NY,NX)
    trcsa_solml2(idsa_FeH2PO4B,L,NY,NX)=trcsa_solml(idsa_FeH2PO4B,L,NY,NX)
    trcsa_solml2(idsa_CaPO4B,L,NY,NX)=trcsa_solml(idsa_CaPO4B,L,NY,NX)
    trcsa_solml2(idsa_CaHPO4B,L,NY,NX)=trcsa_solml(idsa_CaHPO4B,L,NY,NX)
    trcsa_solml2(idsa_CaH2PO4B,L,NY,NX)=trcsa_solml(idsa_CaH2PO4B,L,NY,NX)
    trcsa_solml2(idsa_MgHPO4B,L,NY,NX)=trcsa_solml(idsa_MgHPO4B,L,NY,NX)

    ZALH2(L,NY,NX)=trcsa_soHml(idsa_Al,L,NY,NX)
    ZFEH2(L,NY,NX)=trcsa_soHml(idsa_Fe,L,NY,NX)
    ZHYH2(L,NY,NX)=trcsa_soHml(idsa_Hp,L,NY,NX)
    ZCCH2(L,NY,NX)=trcsa_soHml(idsa_Ca,L,NY,NX)
    ZMAH2(L,NY,NX)=trcsa_soHml(idsa_Mg,L,NY,NX)
    ZNAH2(L,NY,NX)=trcsa_soHml(idsa_Na,L,NY,NX)
    ZKAH2(L,NY,NX)=trcsa_soHml(idsa_K,L,NY,NX)
    ZOHH2(L,NY,NX)=trcsa_soHml(idsa_OH,L,NY,NX)
    ZSO4H2(L,NY,NX)=trcsa_soHml(idsa_SO4,L,NY,NX)
    ZCLH2(L,NY,NX)=trcsa_soHml(idsa_Cl,L,NY,NX)
    ZCO3H2(L,NY,NX)=trcsa_soHml(idsa_CO3,L,NY,NX)
    ZHCOH2(L,NY,NX)=trcsa_soHml(idsa_HCO3,L,NY,NX)
    ZAL1H2(L,NY,NX)=trcsa_soHml(idsa_AlOH,L,NY,NX)
    ZAL2H2(L,NY,NX)=trcsa_soHml(idsa_AlOH2,L,NY,NX)
    ZAL3H2(L,NY,NX)=trcsa_soHml(idsa_AlOH3,L,NY,NX)
    ZAL4H2(L,NY,NX)=trcsa_soHml(idsa_AlOH4,L,NY,NX)
    ZALSH2(L,NY,NX)=trcsa_soHml(idsa_AlSO4,L,NY,NX)
    ZFE1H2(L,NY,NX)=trcsa_soHml(idsa_FeOH,L,NY,NX)
    ZFE2H2(L,NY,NX)=trcsa_soHml(idsa_FeOH2,L,NY,NX)
    ZFE3H2(L,NY,NX)=trcsa_soHml(idsa_FeOH3,L,NY,NX)
    ZFE4H2(L,NY,NX)=trcsa_soHml(idsa_FeOH4,L,NY,NX)
    ZFESH2(L,NY,NX)=trcsa_soHml(idsa_FeSO4,L,NY,NX)
    ZCAOH2(L,NY,NX)=trcsa_soHml(idsa_CaOH2,L,NY,NX)
    ZCACH2(L,NY,NX)=trcsa_soHml(idsa_CaCO3,L,NY,NX)
    ZCAHH2(L,NY,NX)=trcsa_soHml(idsa_CaHCO3,L,NY,NX)
    ZCASH2(L,NY,NX)=trcsa_soHml(idsa_CaSO4,L,NY,NX)
    ZMGOH2(L,NY,NX)=trcsa_soHml(idsa_MgOH2,L,NY,NX)
    ZMGCH2(L,NY,NX)=trcsa_soHml(idsa_MgCO3,L,NY,NX)
    ZMGHH2(L,NY,NX)=trcsa_soHml(idsa_MgHCO3,L,NY,NX)
    ZMGSH2(L,NY,NX)=trcsa_soHml(idsa_MgSO4,L,NY,NX)
    ZNACH2(L,NY,NX)=trcsa_soHml(idsa_NaCO3,L,NY,NX)
    ZNASH2(L,NY,NX)=trcsa_soHml(idsa_NaSO4,L,NY,NX)
    ZKASH2(L,NY,NX)=trcsa_soHml(idsa_KSO4,L,NY,NX)
    H0P4H2(L,NY,NX)=trcsa_soHml(idsa_H0PO4,L,NY,NX)
    H3P4H2(L,NY,NX)=trcsa_soHml(idsa_H3PO4,L,NY,NX)
    ZF1PH2(L,NY,NX)=trcsa_soHml(idsa_FeHPO4,L,NY,NX)
    ZF2PH2(L,NY,NX)=trcsa_soHml(idsa_FeH2PO4,L,NY,NX)
    ZC0PH2(L,NY,NX)=trcsa_soHml(idsa_CaPO4,L,NY,NX)
    ZC1PH2(L,NY,NX)=trcsa_soHml(idsa_CaHPO4,L,NY,NX)
    ZC2PH2(L,NY,NX)=trcsa_soHml(idsa_CaH2PO4,L,NY,NX)
    ZM1PH2(L,NY,NX)=trcsa_soHml(idsa_MgHPO4,L,NY,NX)
    H0PBH2(L,NY,NX)=trcsa_soHml(idsa_H0PO4B,L,NY,NX)
    H3PBH2(L,NY,NX)=trcsa_soHml(idsa_H3PO4B,L,NY,NX)
    ZF1BH2(L,NY,NX)=trcsa_soHml(idsa_FeHPO4B,L,NY,NX)
    ZF2BH2(L,NY,NX)=trcsa_soHml(idsa_FeH2PO4B,L,NY,NX)
    ZC0BH2(L,NY,NX)=trcsa_soHml(idsa_CaPO4B,L,NY,NX)
    ZC1BH2(L,NY,NX)=trcsa_soHml(idsa_CaHPO4B,L,NY,NX)
    ZC2BH2(L,NY,NX)=trcsa_soHml(idsa_CaH2PO4B,L,NY,NX)
    ZM1BH2(L,NY,NX)=trcsa_soHml(idsa_MgHPO4B,L,NY,NX)
  ENDDO D10
  end subroutine GetSubHourlyFluxByLayer
!------------------------------------------------------------------------------------------

      subroutine SaltModelSoluteFlux(I,NHW,NHE,NVN,NVS)
!
!     Description:
!
      implicit none
      integer, intent(in) :: I,NHW,NHE,NVN,NVS

      integer :: NY,NX
!     begin_execution

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
!     PRECW,PRECR=snow,rain
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     X*BLS=hourly solute flux to snowpack
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!
      IF(PRECW(NY,NX).GT.0.0.OR.(PRECR(NY,NX).GT.0.0 &
      .AND.VHCPWM(1,1,NY,NX).GT.VHCPWX(NY,NX)))THEN
!
      call AtmosSoluteFluxToSnowpack(I,NY,NX)
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
      ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0) &
        .AND.VHCPWM(1,1,NY,NX).LE.VHCPWX(NY,NX))THEN
!
      call AtmosSoluteFluxToTopsoil(I,NY,NX)
!
!     NO SOLUTE FLUXES FROM ATMOSPHERE
!
      ELSE
      call ZeroAtmosSoluteFlux(NY,NX)
      ENDIF
!
!     GAS AND SOLUTE FLUXES AT SUB-HOURLY FLUX TIME STEP
!     ENTERED IN SITE FILE
      call GetSubHourFlux(NY,NX)
!
!     INITIAL SOLUTES IN SNOWPACK
      call InitSolutesInSnowpack(NY,NX)
!
!     SOLUTE FLUXES FROM SOLUTE.F
      call GetSubHourlyFluxByLayer(I,NY,NX)

    ENDDO D9990
  ENDDO D9995
  end subroutine SaltModelSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine InitFluxAccumulatorsInSnowpack(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,NTSA

  D9855: DO L=1,JS
    DO NTSA=idsa_beg,idsa_end
      trcsa_TBLS(NTSA,L,NY,NX)=0.0
    ENDDO
  ENDDO D9855
  end subroutine InitFluxAccumulatorsInSnowpack

!------------------------------------------------------------------------------------------

  subroutine BoundarySnowFlux(N,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,M5,M4
!     begin_execution
  trcsa_RQ(idsa_Al,N,M5,M4)=0.0
  trcsa_RQ(idsa_Fe,N,M5,M4)=0.0
  trcsa_RQ(idsa_Hp,N,M5,M4)=0.0
  trcsa_RQ(idsa_Ca,N,M5,M4)=0.0
  trcsa_RQ(idsa_Mg,N,M5,M4)=0.0
  trcsa_RQ(idsa_Na,N,M5,M4)=0.0
  trcsa_RQ(idsa_K,N,M5,M4)=0.0
  trcsa_RQ(idsa_OH,N,M5,M4)=0.0
  trcsa_RQ(idsa_SO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_Cl,N,M5,M4)=0.0
  trcsa_RQ(idsa_CO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_HCO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_AlOH,N,M5,M4)=0.0
  trcsa_RQ(idsa_AlOH2,N,M5,M4)=0.0
  trcsa_RQ(idsa_AlOH3,N,M5,M4)=0.0
  trcsa_RQ(idsa_AlOH4,N,M5,M4)=0.0
  trcsa_RQ(idsa_AlSO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeOH,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeOH2,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeOH3,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeOH4,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeSO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaOH2,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaCO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaHCO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaSO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_MgOH2,N,M5,M4)=0.0
  trcsa_RQ(idsa_MgCO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_MgHCO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_MgSO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_NaCO3,N,M5,M4)=0.0
  trcsa_RQ(idsa_NaSO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_KSO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_H0PO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_H3PO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeHPO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_FeH2PO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaPO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaHPO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_CaH2PO4,N,M5,M4)=0.0
  trcsa_RQ(idsa_MgHPO4,N,M5,M4)=0.0
  end subroutine BoundarySnowFlux
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteFluxFromRecharge(N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,NN,M5,M4
!     begin_execution

  RQRAL(N,NN,M5,M4)=0.0
  RQRFE(N,NN,M5,M4)=0.0
  RQRHY(N,NN,M5,M4)=0.0
  RQRCA(N,NN,M5,M4)=0.0
  RQRMG(N,NN,M5,M4)=0.0
  RQRNA(N,NN,M5,M4)=0.0
  RQRKA(N,NN,M5,M4)=0.0
  RQROH(N,NN,M5,M4)=0.0
  RQRSO(N,NN,M5,M4)=0.0
  RQRCL(N,NN,M5,M4)=0.0
  RQRC3(N,NN,M5,M4)=0.0
  RQRHC(N,NN,M5,M4)=0.0
  RQRAL1(N,NN,M5,M4)=0.0
  RQRAL2(N,NN,M5,M4)=0.0
  RQRAL3(N,NN,M5,M4)=0.0
  RQRAL4(N,NN,M5,M4)=0.0
  RQRALS(N,NN,M5,M4)=0.0
  RQRFE1(N,NN,M5,M4)=0.0
  RQRFE2(N,NN,M5,M4)=0.0
  RQRFE3(N,NN,M5,M4)=0.0
  RQRFE4(N,NN,M5,M4)=0.0
  RQRFES(N,NN,M5,M4)=0.0
  RQRCAO(N,NN,M5,M4)=0.0
  RQRCAC(N,NN,M5,M4)=0.0
  RQRCAH(N,NN,M5,M4)=0.0
  RQRCAS(N,NN,M5,M4)=0.0
  RQRMGO(N,NN,M5,M4)=0.0
  RQRMGC(N,NN,M5,M4)=0.0
  RQRMGH(N,NN,M5,M4)=0.0
  RQRMGS(N,NN,M5,M4)=0.0
  RQRNAC(N,NN,M5,M4)=0.0
  RQRNAS(N,NN,M5,M4)=0.0
  RQRKAS(N,NN,M5,M4)=0.0
  RQRH0P(N,NN,M5,M4)=0.0
  RQRH3P(N,NN,M5,M4)=0.0
  RQRF1P(N,NN,M5,M4)=0.0
  RQRF2P(N,NN,M5,M4)=0.0
  RQRC0P(N,NN,M5,M4)=0.0
  RQRC1P(N,NN,M5,M4)=0.0
  RQRC2P(N,NN,M5,M4)=0.0
  RQRM1P(N,NN,M5,M4)=0.0
  end subroutine ZeroSoluteFluxFromRecharge
!------------------------------------------------------------------------------------------

  subroutine SoluteExportThruBoundary(N1,N2,M,N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N1,N2,M,N,NN,M5,M4
real(r8) :: FQRM
!     begin_execution

  FQRM=QRMN(M,N,NN,M5,M4)/QRM(M,N2,N1)
  RQRAL(N,NN,M5,M4)=RQRAL0(N2,N1)*FQRM
  RQRFE(N,NN,M5,M4)=RQRFE0(N2,N1)*FQRM
  RQRHY(N,NN,M5,M4)=RQRHY0(N2,N1)*FQRM
  RQRCA(N,NN,M5,M4)=RQRCA0(N2,N1)*FQRM
  RQRMG(N,NN,M5,M4)=RQRMG0(N2,N1)*FQRM
  RQRNA(N,NN,M5,M4)=RQRNA0(N2,N1)*FQRM
  RQRKA(N,NN,M5,M4)=RQRKA0(N2,N1)*FQRM
  RQROH(N,NN,M5,M4)=RQROH0(N2,N1)*FQRM
  RQRSO(N,NN,M5,M4)=RQRSO0(N2,N1)*FQRM
  RQRCL(N,NN,M5,M4)=RQRCL0(N2,N1)*FQRM
  RQRC3(N,NN,M5,M4)=RQRC30(N2,N1)*FQRM
  RQRHC(N,NN,M5,M4)=RQRHC0(N2,N1)*FQRM
  RQRAL1(N,NN,M5,M4)=RQRAL10(N2,N1)*FQRM
  RQRAL2(N,NN,M5,M4)=RQRAL20(N2,N1)*FQRM
  RQRAL3(N,NN,M5,M4)=RQRAL30(N2,N1)*FQRM
  RQRAL4(N,NN,M5,M4)=RQRAL40(N2,N1)*FQRM
  RQRALS(N,NN,M5,M4)=RQRALS0(N2,N1)*FQRM
  RQRFE1(N,NN,M5,M4)=RQRFE10(N2,N1)*FQRM
  RQRFE2(N,NN,M5,M4)=RQRFE20(N2,N1)*FQRM
  RQRFE3(N,NN,M5,M4)=RQRFE30(N2,N1)*FQRM
  RQRFE4(N,NN,M5,M4)=RQRFE40(N2,N1)*FQRM
  RQRFES(N,NN,M5,M4)=RQRFES0(N2,N1)*FQRM
  RQRCAO(N,NN,M5,M4)=RQRCAO0(N2,N1)*FQRM
  RQRCAC(N,NN,M5,M4)=RQRCAC0(N2,N1)*FQRM
  RQRCAH(N,NN,M5,M4)=RQRCAH0(N2,N1)*FQRM
  RQRCAS(N,NN,M5,M4)=RQRCAS0(N2,N1)*FQRM
  RQRMGO(N,NN,M5,M4)=RQRMGO0(N2,N1)*FQRM
  RQRMGC(N,NN,M5,M4)=RQRMGC0(N2,N1)*FQRM
  RQRMGH(N,NN,M5,M4)=RQRMGH0(N2,N1)*FQRM
  RQRMGS(N,NN,M5,M4)=RQRMGS0(N2,N1)*FQRM
  RQRNAC(N,NN,M5,M4)=RQRNAC0(N2,N1)*FQRM
  RQRNAS(N,NN,M5,M4)=RQRNAS0(N2,N1)*FQRM
  RQRKAS(N,NN,M5,M4)=RQRKAS0(N2,N1)*FQRM
  RQRH0P(N,NN,M5,M4)=RQRH0P0(N2,N1)*FQRM
  RQRH3P(N,NN,M5,M4)=RQRH3P0(N2,N1)*FQRM
  RQRF1P(N,NN,M5,M4)=RQRF1P0(N2,N1)*FQRM
  RQRF2P(N,NN,M5,M4)=RQRF2P0(N2,N1)*FQRM
  RQRC0P(N,NN,M5,M4)=RQRC0P0(N2,N1)*FQRM
  RQRC1P(N,NN,M5,M4)=RQRC1P0(N2,N1)*FQRM
  RQRC2P(N,NN,M5,M4)=RQRC2P0(N2,N1)*FQRM
  RQRM1P(N,NN,M5,M4)=RQRM1P0(N2,N1)*FQRM
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QRS=hourly solute in runoff
!     RQR*=solute in runoff
!
  XQRAL(N,NN,M5,M4)=XQRAL(N,NN,M5,M4)+RQRAL(N,NN,M5,M4)
  XQRFE(N,NN,M5,M4)=XQRFE(N,NN,M5,M4)+RQRFE(N,NN,M5,M4)
  XQRHY(N,NN,M5,M4)=XQRHY(N,NN,M5,M4)+RQRHY(N,NN,M5,M4)
  XQRCA(N,NN,M5,M4)=XQRCA(N,NN,M5,M4)+RQRCA(N,NN,M5,M4)
  XQRMG(N,NN,M5,M4)=XQRMG(N,NN,M5,M4)+RQRMG(N,NN,M5,M4)
  XQRNA(N,NN,M5,M4)=XQRNA(N,NN,M5,M4)+RQRNA(N,NN,M5,M4)
  XQRKA(N,NN,M5,M4)=XQRKA(N,NN,M5,M4)+RQRKA(N,NN,M5,M4)
  XQROH(N,NN,M5,M4)=XQROH(N,NN,M5,M4)+RQROH(N,NN,M5,M4)
  XQRSO(N,NN,M5,M4)=XQRSO(N,NN,M5,M4)+RQRSO(N,NN,M5,M4)
  XQRCL(N,NN,M5,M4)=XQRCL(N,NN,M5,M4)+RQRCL(N,NN,M5,M4)
  XQRC3(N,NN,M5,M4)=XQRC3(N,NN,M5,M4)+RQRC3(N,NN,M5,M4)
  XQRHC(N,NN,M5,M4)=XQRHC(N,NN,M5,M4)+RQRHC(N,NN,M5,M4)
  XQRAL1(N,NN,M5,M4)=XQRAL1(N,NN,M5,M4)+RQRAL1(N,NN,M5,M4)
  XQRAL2(N,NN,M5,M4)=XQRAL2(N,NN,M5,M4)+RQRAL2(N,NN,M5,M4)
  XQRAL3(N,NN,M5,M4)=XQRAL3(N,NN,M5,M4)+RQRAL3(N,NN,M5,M4)
  XQRAL4(N,NN,M5,M4)=XQRAL4(N,NN,M5,M4)+RQRAL4(N,NN,M5,M4)
  XQRALS(N,NN,M5,M4)=XQRALS(N,NN,M5,M4)+RQRALS(N,NN,M5,M4)
  XQRFE1(N,NN,M5,M4)=XQRFE1(N,NN,M5,M4)+RQRFE1(N,NN,M5,M4)
  XQRFE2(N,NN,M5,M4)=XQRFE2(N,NN,M5,M4)+RQRFE2(N,NN,M5,M4)
  XQRFE3(N,NN,M5,M4)=XQRFE3(N,NN,M5,M4)+RQRFE3(N,NN,M5,M4)
  XQRFE4(N,NN,M5,M4)=XQRFE4(N,NN,M5,M4)+RQRFE4(N,NN,M5,M4)
  XQRFES(N,NN,M5,M4)=XQRFES(N,NN,M5,M4)+RQRFES(N,NN,M5,M4)
  XQRCAO(N,NN,M5,M4)=XQRCAO(N,NN,M5,M4)+RQRCAO(N,NN,M5,M4)
  XQRCAC(N,NN,M5,M4)=XQRCAC(N,NN,M5,M4)+RQRCAC(N,NN,M5,M4)
  XQRCAH(N,NN,M5,M4)=XQRCAH(N,NN,M5,M4)+RQRCAH(N,NN,M5,M4)
  XQRCAS(N,NN,M5,M4)=XQRCAS(N,NN,M5,M4)+RQRCAS(N,NN,M5,M4)
  XQRMGO(N,NN,M5,M4)=XQRMGO(N,NN,M5,M4)+RQRMGO(N,NN,M5,M4)
  XQRMGC(N,NN,M5,M4)=XQRMGC(N,NN,M5,M4)+RQRMGC(N,NN,M5,M4)
  XQRMGH(N,NN,M5,M4)=XQRMGH(N,NN,M5,M4)+RQRMGH(N,NN,M5,M4)
  XQRMGS(N,NN,M5,M4)=XQRMGS(N,NN,M5,M4)+RQRMGS(N,NN,M5,M4)
  XQRNAC(N,NN,M5,M4)=XQRNAC(N,NN,M5,M4)+RQRNAC(N,NN,M5,M4)
  XQRNAS(N,NN,M5,M4)=XQRNAS(N,NN,M5,M4)+RQRNAS(N,NN,M5,M4)
  XQRKAS(N,NN,M5,M4)=XQRKAS(N,NN,M5,M4)+RQRKAS(N,NN,M5,M4)
  XQRH0P(N,NN,M5,M4)=XQRH0P(N,NN,M5,M4)+RQRH0P(N,NN,M5,M4)
  XQRH3P(N,NN,M5,M4)=XQRH3P(N,NN,M5,M4)+RQRH3P(N,NN,M5,M4)
  XQRF1P(N,NN,M5,M4)=XQRF1P(N,NN,M5,M4)+RQRF1P(N,NN,M5,M4)
  XQRF2P(N,NN,M5,M4)=XQRF2P(N,NN,M5,M4)+RQRF2P(N,NN,M5,M4)
  XQRC0P(N,NN,M5,M4)=XQRC0P(N,NN,M5,M4)+RQRC0P(N,NN,M5,M4)
  XQRC1P(N,NN,M5,M4)=XQRC1P(N,NN,M5,M4)+RQRC1P(N,NN,M5,M4)
  XQRC2P(N,NN,M5,M4)=XQRC2P(N,NN,M5,M4)+RQRC2P(N,NN,M5,M4)
  XQRM1P(N,NN,M5,M4)=XQRM1P(N,NN,M5,M4)+RQRM1P(N,NN,M5,M4)
  end subroutine SoluteExportThruBoundary
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,NN,M5,M4
!     begin_execution

  RQRAL(N,NN,M5,M4)=0.0
  RQRFE(N,NN,M5,M4)=0.0
  RQRHY(N,NN,M5,M4)=0.0
  RQRCA(N,NN,M5,M4)=0.0
  RQRMG(N,NN,M5,M4)=0.0
  RQRNA(N,NN,M5,M4)=0.0
  RQRKA(N,NN,M5,M4)=0.0
  RQROH(N,NN,M5,M4)=0.0
  RQRSO(N,NN,M5,M4)=0.0
  RQRCL(N,NN,M5,M4)=0.0
  RQRC3(N,NN,M5,M4)=0.0
  RQRHC(N,NN,M5,M4)=0.0
  RQRAL1(N,NN,M5,M4)=0.0
  RQRAL2(N,NN,M5,M4)=0.0
  RQRAL3(N,NN,M5,M4)=0.0
  RQRAL4(N,NN,M5,M4)=0.0
  RQRALS(N,NN,M5,M4)=0.0
  RQRFE1(N,NN,M5,M4)=0.0
  RQRFE2(N,NN,M5,M4)=0.0
  RQRFE3(N,NN,M5,M4)=0.0
  RQRFE4(N,NN,M5,M4)=0.0
  RQRFES(N,NN,M5,M4)=0.0
  RQRCAO(N,NN,M5,M4)=0.0
  RQRCAC(N,NN,M5,M4)=0.0
  RQRCAH(N,NN,M5,M4)=0.0
  RQRCAS(N,NN,M5,M4)=0.0
  RQRMGO(N,NN,M5,M4)=0.0
  RQRMGC(N,NN,M5,M4)=0.0
  RQRMGH(N,NN,M5,M4)=0.0
  RQRMGS(N,NN,M5,M4)=0.0
  RQRNAC(N,NN,M5,M4)=0.0
  RQRNAS(N,NN,M5,M4)=0.0
  RQRKAS(N,NN,M5,M4)=0.0
  RQRH0P(N,NN,M5,M4)=0.0
  RQRH3P(N,NN,M5,M4)=0.0
  RQRF1P(N,NN,M5,M4)=0.0
  RQRF2P(N,NN,M5,M4)=0.0
  RQRC0P(N,NN,M5,M4)=0.0
  RQRC1P(N,NN,M5,M4)=0.0
  RQRC2P(N,NN,M5,M4)=0.0
  RQRM1P(N,NN,M5,M4)=0.0
  end subroutine ZeroSoluteInfluxThruBoundary
!------------------------------------------------------------------------------------------

  subroutine SoluteGainSubsurfMicropore(M3,M2,M1,M,N,M6,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M3,M2,M1,M,N,M6,M5,M4
!     begin_execution
  trcsa_RFLS(idsa_Al,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CALU(M3,M2,M1)
  trcsa_RFLS(idsa_Fe,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFEU(M3,M2,M1)
  trcsa_RFLS(idsa_Hp,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CHYU(M3,M2,M1)
  trcsa_RFLS(idsa_Ca,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCAU(M3,M2,M1)
  trcsa_RFLS(idsa_Mg,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGU(M3,M2,M1)
  trcsa_RFLS(idsa_Na,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNAU(M3,M2,M1)
  trcsa_RFLS(idsa_K,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CKAU(M3,M2,M1)
  trcsa_RFLS(idsa_OH,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*COHU(M3,M2,M1)
  trcsa_RFLS(idsa_SO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CSOU(M3,M2,M1)
  trcsa_RFLS(idsa_Cl,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCLU(M3,M2,M1)
  trcsa_RFLS(idsa_CO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC3U(M3,M2,M1)
  trcsa_RFLS(idsa_HCO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CHCU(M3,M2,M1)
  trcsa_RFLS(idsa_AlOH,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL1U(M3,M2,M1)
  trcsa_RFLS(idsa_AlOH2,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL2U(M3,M2,M1)
  trcsa_RFLS(idsa_AlOH3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL3U(M3,M2,M1)
  trcsa_RFLS(idsa_AlOH4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL4U(M3,M2,M1)
  trcsa_RFLS(idsa_AlSO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CALSU(M3,M2,M1)
  trcsa_RFLS(idsa_FeOH,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE1U(M3,M2,M1)
  trcsa_RFLS(idsa_FeOH2,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE2U(M3,M2,M1)
  trcsa_RFLS(idsa_FeOH3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE3U(M3,M2,M1)
  trcsa_RFLS(idsa_FeOH4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE4U(M3,M2,M1)
  trcsa_RFLS(idsa_FeSO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFESU(M3,M2,M1)
  trcsa_RFLS(idsa_CaOH2,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCAOU(M3,M2,M1)
  trcsa_RFLS(idsa_CaCO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCACU(M3,M2,M1)
  trcsa_RFLS(idsa_CaHCO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCAHU(M3,M2,M1)
  trcsa_RFLS(idsa_CaSO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCASU(M3,M2,M1)
  trcsa_RFLS(idsa_MgOH2,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGOU(M3,M2,M1)
  trcsa_RFLS(idsa_MgCO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGCU(M3,M2,M1)
  trcsa_RFLS(idsa_MgHCO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGHU(M3,M2,M1)
  trcsa_RFLS(idsa_MgSO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGSU(M3,M2,M1)
  trcsa_RFLS(idsa_NaCO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNACU(M3,M2,M1)
  trcsa_RFLS(idsa_NaSO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNASU(M3,M2,M1)
  trcsa_RFLS(idsa_KSO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CKASU(M3,M2,M1)
  trcsa_RFLS(idsa_H0PO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH0PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_H3PO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH3PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_FeHPO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_FeH2PO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF2PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_CaPO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC0PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_CaHPO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_CaH2PO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC2PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_MgHPO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CM1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
  trcsa_RFLS(idsa_H0PO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH0PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_H3PO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH3PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_FeHPO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_FeH2PO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF2PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_CaPO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC0PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_CaHPO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_CaH2PO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC2PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  trcsa_RFLS(idsa_MgHPO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CM1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  end subroutine SoluteGainSubsurfMicropore
!------------------------------------------------------------------------------------------

      subroutine SoluteLossSubsurfMicropore(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
      implicit none
      integer, intent(in) :: M,N,M1,M2,M3,M4,M5,M6
  real(r8) :: VFLW
!     begin_execution

      IF(VOLWM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWM(M,N,M6,M5,M4) &
      /VOLWM(M,M3,M2,M1)))
      ELSE
      VFLW=0.0
      ENDIF
      trcsa_RFLS(idsa_Al,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Al,M3,M2,M1))
      trcsa_RFLS(idsa_Fe,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,M3,M2,M1))
      trcsa_RFLS(idsa_Hp,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,M3,M2,M1))
      trcsa_RFLS(idsa_Ca,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,M3,M2,M1))
      trcsa_RFLS(idsa_Mg,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,M3,M2,M1))
      trcsa_RFLS(idsa_Na,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Na,M3,M2,M1))
      trcsa_RFLS(idsa_K,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_K,M3,M2,M1))
      trcsa_RFLS(idsa_OH,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_OH,M3,M2,M1))
      trcsa_RFLS(idsa_SO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,M3,M2,M1))
      trcsa_RFLS(idsa_Cl,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,M3,M2,M1))
      trcsa_RFLS(idsa_CO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,M3,M2,M1))
      trcsa_RFLS(idsa_HCO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,M3,M2,M1))
      trcsa_RFLS(idsa_AlOH,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,M3,M2,M1))
      trcsa_RFLS(idsa_AlOH2,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,M3,M2,M1))
      trcsa_RFLS(idsa_AlOH3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,M3,M2,M1))
      trcsa_RFLS(idsa_AlOH4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,M3,M2,M1))
      trcsa_RFLS(idsa_AlSO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,M3,M2,M1))
      trcsa_RFLS(idsa_FeOH,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,M3,M2,M1))
      trcsa_RFLS(idsa_FeOH2,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,M3,M2,M1))
      trcsa_RFLS(idsa_FeOH3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,M3,M2,M1))
      trcsa_RFLS(idsa_FeOH4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,M3,M2,M1))
      trcsa_RFLS(idsa_FeSO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,M3,M2,M1))
      trcsa_RFLS(idsa_CaOH2,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,M3,M2,M1))
      trcsa_RFLS(idsa_CaCO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,M3,M2,M1))
      trcsa_RFLS(idsa_CaHCO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,M3,M2,M1))
      trcsa_RFLS(idsa_CaSO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,M3,M2,M1))
      trcsa_RFLS(idsa_MgOH2,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,M3,M2,M1))
      trcsa_RFLS(idsa_MgCO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,M3,M2,M1))
      trcsa_RFLS(idsa_MgHCO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,M3,M2,M1))
      trcsa_RFLS(idsa_MgSO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,M3,M2,M1))
      trcsa_RFLS(idsa_NaCO3,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,M3,M2,M1))
      trcsa_RFLS(idsa_NaSO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,M3,M2,M1))
      trcsa_RFLS(idsa_KSO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,M3,M2,M1))
      trcsa_RFLS(idsa_H0PO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_H3PO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_FeHPO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_FeH2PO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_CaPO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_CaHPO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_CaH2PO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_MgHPO4,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
      trcsa_RFLS(idsa_H0PO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_H3PO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_FeHPO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_FeH2PO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_CaPO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_CaHPO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_CaH2PO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      trcsa_RFLS(idsa_MgHPO4B,N,M6,M5,M4)=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
      end subroutine SoluteLossSubsurfMicropore
!------------------------------------------------------------------------------------------

  subroutine SoluteLossSubsurfMacropore(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N,M1,M2,M3,M4,M5,M6
  real(r8) :: VFLW
!     begin_execution

  IF(VOLWHM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
    VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWHM(M,N,M6,M5,M4)/VOLWHM(M,M3,M2,M1)))
  ELSE
    VFLW=0.0
  ENDIF
  RALFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZALH2(M3,M2,M1))
  RFEFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZFEH2(M3,M2,M1))
  RHYFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZHYH2(M3,M2,M1))
  RCAFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCCH2(M3,M2,M1))
  RMGFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZMAH2(M3,M2,M1))
  RNAFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZNAH2(M3,M2,M1))
  RKAFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZKAH2(M3,M2,M1))
  ROHFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZOHH2(M3,M2,M1))
  RSOFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZSO4H2(M3,M2,M1))
  RCLFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCLH2(M3,M2,M1))
  RC3FHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCO3H2(M3,M2,M1))
  RHCFHS(N,M6,M5,M4)=VFLW*AZMAX1(ZHCOH2(M3,M2,M1))
  RAL1HS(N,M6,M5,M4)=VFLW*AZMAX1(ZAL1H2(M3,M2,M1))
  RAL2HS(N,M6,M5,M4)=VFLW*AZMAX1(ZAL2H2(M3,M2,M1))
  RAL3HS(N,M6,M5,M4)=VFLW*AZMAX1(ZAL3H2(M3,M2,M1))
  RAL4HS(N,M6,M5,M4)=VFLW*AZMAX1(ZAL4H2(M3,M2,M1))
  RALSHS(N,M6,M5,M4)=VFLW*AZMAX1(ZALSH2(M3,M2,M1))
  RFE1HS(N,M6,M5,M4)=VFLW*AZMAX1(ZFE1H2(M3,M2,M1))
  RFE2HS(N,M6,M5,M4)=VFLW*AZMAX1(ZFE2H2(M3,M2,M1))
  RFE3HS(N,M6,M5,M4)=VFLW*AZMAX1(ZFE3H2(M3,M2,M1))
  RFE4HS(N,M6,M5,M4)=VFLW*AZMAX1(ZFE4H2(M3,M2,M1))
  RFESHS(N,M6,M5,M4)=VFLW*AZMAX1(ZFESH2(M3,M2,M1))
  RCAOHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCAOH2(M3,M2,M1))
  RCACHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCACH2(M3,M2,M1))
  RCAHHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCAHH2(M3,M2,M1))
  RCASHS(N,M6,M5,M4)=VFLW*AZMAX1(ZCASH2(M3,M2,M1))
  RMGOHS(N,M6,M5,M4)=VFLW*AZMAX1(ZMGOH2(M3,M2,M1))
  RMGCHS(N,M6,M5,M4)=VFLW*AZMAX1(ZMGCH2(M3,M2,M1))
  RMGHHS(N,M6,M5,M4)=VFLW*AZMAX1(ZMGHH2(M3,M2,M1))
  RMGSHS(N,M6,M5,M4)=VFLW*AZMAX1(ZMGSH2(M3,M2,M1))
  RNACHS(N,M6,M5,M4)=VFLW*AZMAX1(ZNACH2(M3,M2,M1))
  RNASHS(N,M6,M5,M4)=VFLW*AZMAX1(ZNASH2(M3,M2,M1))
  RKASHS(N,M6,M5,M4)=VFLW*AZMAX1(ZKASH2(M3,M2,M1))
  RH0PHS(N,M6,M5,M4)=VFLW*AZMAX1(H0P4H2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RH3PHS(N,M6,M5,M4)=VFLW*AZMAX1(H3P4H2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RF1PHS(N,M6,M5,M4)=VFLW*AZMAX1(ZF1PH2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RF2PHS(N,M6,M5,M4)=VFLW*AZMAX1(ZF2PH2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RC0PHS(N,M6,M5,M4)=VFLW*AZMAX1(ZC0PH2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RC1PHS(N,M6,M5,M4)=VFLW*AZMAX1(ZC1PH2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RC2PHS(N,M6,M5,M4)=VFLW*AZMAX1(ZC2PH2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RM1PHS(N,M6,M5,M4)=VFLW*AZMAX1(ZM1PH2(M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RH0BHB(N,M6,M5,M4)=VFLW*AZMAX1(H0PBH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RH3BHB(N,M6,M5,M4)=VFLW*AZMAX1(H3PBH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RF1BHB(N,M6,M5,M4)=VFLW*AZMAX1(ZF1BH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RF2BHB(N,M6,M5,M4)=VFLW*AZMAX1(ZF2BH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RC0BHB(N,M6,M5,M4)=VFLW*AZMAX1(ZC0BH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RC1BHB(N,M6,M5,M4)=VFLW*AZMAX1(ZC1BH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RC2BHB(N,M6,M5,M4)=VFLW*AZMAX1(ZC2BH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RM1BHB(N,M6,M5,M4)=VFLW*AZMAX1(ZM1BH2(M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  end subroutine SoluteLossSubsurfMacropore
!------------------------------------------------------------------------------------------

  subroutine ZeroSolueGainSubsurfMacropore(N,M6,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,M6,M5,M4
!     begin_execution
  RALFHS(N,M6,M5,M4)=0.0
  RFEFHS(N,M6,M5,M4)=0.0
  RHYFHS(N,M6,M5,M4)=0.0
  RCAFHS(N,M6,M5,M4)=0.0
  RMGFHS(N,M6,M5,M4)=0.0
  RNAFHS(N,M6,M5,M4)=0.0
  RKAFHS(N,M6,M5,M4)=0.0
  ROHFHS(N,M6,M5,M4)=0.0
  RSOFHS(N,M6,M5,M4)=0.0
  RCLFHS(N,M6,M5,M4)=0.0
  RC3FHS(N,M6,M5,M4)=0.0
  RHCFHS(N,M6,M5,M4)=0.0
  RAL1HS(N,M6,M5,M4)=0.0
  RAL2HS(N,M6,M5,M4)=0.0
  RAL3HS(N,M6,M5,M4)=0.0
  RAL4HS(N,M6,M5,M4)=0.0
  RALSHS(N,M6,M5,M4)=0.0
  RFE1HS(N,M6,M5,M4)=0.0
  RFE2HS(N,M6,M5,M4)=0.0
  RFE3HS(N,M6,M5,M4)=0.0
  RFE4HS(N,M6,M5,M4)=0.0
  RFESHS(N,M6,M5,M4)=0.0
  RCAOHS(N,M6,M5,M4)=0.0
  RCACHS(N,M6,M5,M4)=0.0
  RCAHHS(N,M6,M5,M4)=0.0
  RCASHS(N,M6,M5,M4)=0.0
  RMGOHS(N,M6,M5,M4)=0.0
  RMGCHS(N,M6,M5,M4)=0.0
  RMGHHS(N,M6,M5,M4)=0.0
  RMGSHS(N,M6,M5,M4)=0.0
  RNACHS(N,M6,M5,M4)=0.0
  RNASHS(N,M6,M5,M4)=0.0
  RKASHS(N,M6,M5,M4)=0.0
  RH0PHS(N,M6,M5,M4)=0.0
  RH3PHS(N,M6,M5,M4)=0.0
  RF1PHS(N,M6,M5,M4)=0.0
  RF2PHS(N,M6,M5,M4)=0.0
  RC0PHS(N,M6,M5,M4)=0.0
  RC1PHS(N,M6,M5,M4)=0.0
  RC2PHS(N,M6,M5,M4)=0.0
  RM1PHS(N,M6,M5,M4)=0.0
  RH0BHB(N,M6,M5,M4)=0.0
  RH3BHB(N,M6,M5,M4)=0.0
  RF1BHB(N,M6,M5,M4)=0.0
  RF2BHB(N,M6,M5,M4)=0.0
  RC0BHB(N,M6,M5,M4)=0.0
  RC1BHB(N,M6,M5,M4)=0.0
  RC2BHB(N,M6,M5,M4)=0.0
  RM1BHB(N,M6,M5,M4)=0.0
  end subroutine ZeroSolueGainSubsurfMacropore
!------------------------------------------------------------------------------------------

  subroutine AccumFluxMacMicPores(N,M6,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,M6,M5,M4
!     begin_execution
!
!     X*FLS,X*FLW,X*FLB=hourly solute flux in non-band,band micropores
!     X*FHS,X*FHW,X*FHB=hourly solute flux in non-band,band macropores
!     R*FLS,R*FLW,R*FLB=solute flux in non-band,band micropores
!     R*FHS,R*FHW,R*FHB=solute flux in non-band,band macropores
!
  trcsa_XFLS(idsa_Al,N,M6,M5,M4)=trcsa_XFLS(idsa_Al,N,M6,M5,M4)+trcsa_RFLS(idsa_Al,N,M6,M5,M4)
  trcsa_XFLS(idsa_Fe,N,M6,M5,M4)=trcsa_XFLS(idsa_Fe,N,M6,M5,M4)+trcsa_RFLS(idsa_Fe,N,M6,M5,M4)
  trcsa_XFLS(idsa_Hp,N,M6,M5,M4)=trcsa_XFLS(idsa_Hp,N,M6,M5,M4)+trcsa_RFLS(idsa_Hp,N,M6,M5,M4)
  trcsa_XFLS(idsa_Ca,N,M6,M5,M4)=trcsa_XFLS(idsa_Ca,N,M6,M5,M4)+trcsa_RFLS(idsa_Ca,N,M6,M5,M4)
  trcsa_XFLS(idsa_Mg,N,M6,M5,M4)=trcsa_XFLS(idsa_Mg,N,M6,M5,M4)+trcsa_RFLS(idsa_Mg,N,M6,M5,M4)
  trcsa_XFLS(idsa_Na,N,M6,M5,M4)=trcsa_XFLS(idsa_Na,N,M6,M5,M4)+trcsa_RFLS(idsa_Na,N,M6,M5,M4)
  trcsa_XFLS(idsa_K,N,M6,M5,M4)=trcsa_XFLS(idsa_K,N,M6,M5,M4)+trcsa_RFLS(idsa_K,N,M6,M5,M4)
  trcsa_XFLS(idsa_OH,N,M6,M5,M4)=trcsa_XFLS(idsa_OH,N,M6,M5,M4)+trcsa_RFLS(idsa_OH,N,M6,M5,M4)
  trcsa_XFLS(idsa_SO4,N,M6,M5,M4)=trcsa_XFLS(idsa_SO4,N,M6,M5,M4)+trcsa_RFLS(idsa_SO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_Cl,N,M6,M5,M4)=trcsa_XFLS(idsa_Cl,N,M6,M5,M4)+trcsa_RFLS(idsa_Cl,N,M6,M5,M4)
  trcsa_XFLS(idsa_CO3,N,M6,M5,M4)=trcsa_XFLS(idsa_CO3,N,M6,M5,M4)+trcsa_RFLS(idsa_CO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_HCO3,N,M6,M5,M4)=trcsa_XFLS(idsa_HCO3,N,M6,M5,M4)+trcsa_RFLS(idsa_HCO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_AlOH,N,M6,M5,M4)=trcsa_XFLS(idsa_AlOH,N,M6,M5,M4)+trcsa_RFLS(idsa_AlOH,N,M6,M5,M4)
  trcsa_XFLS(idsa_AlOH2,N,M6,M5,M4)=trcsa_XFLS(idsa_AlOH2,N,M6,M5,M4)+trcsa_RFLS(idsa_AlOH2,N,M6,M5,M4)
  trcsa_XFLS(idsa_AlOH3,N,M6,M5,M4)=trcsa_XFLS(idsa_AlOH3,N,M6,M5,M4)+trcsa_RFLS(idsa_AlOH3,N,M6,M5,M4)
  trcsa_XFLS(idsa_AlOH4,N,M6,M5,M4)=trcsa_XFLS(idsa_AlOH4,N,M6,M5,M4)+trcsa_RFLS(idsa_AlOH4,N,M6,M5,M4)
  trcsa_XFLS(idsa_AlSO4,N,M6,M5,M4)=trcsa_XFLS(idsa_AlSO4,N,M6,M5,M4)+trcsa_RFLS(idsa_AlSO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeOH,N,M6,M5,M4)=trcsa_XFLS(idsa_FeOH,N,M6,M5,M4)+trcsa_RFLS(idsa_FeOH,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeOH2,N,M6,M5,M4)=trcsa_XFLS(idsa_FeOH2,N,M6,M5,M4)+trcsa_RFLS(idsa_FeOH2,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeOH3,N,M6,M5,M4)=trcsa_XFLS(idsa_FeOH3,N,M6,M5,M4)+trcsa_RFLS(idsa_FeOH3,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeOH4,N,M6,M5,M4)=trcsa_XFLS(idsa_FeOH4,N,M6,M5,M4)+trcsa_RFLS(idsa_FeOH4,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeSO4,N,M6,M5,M4)=trcsa_XFLS(idsa_FeSO4,N,M6,M5,M4)+trcsa_RFLS(idsa_FeSO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaOH2,N,M6,M5,M4)=trcsa_XFLS(idsa_CaOH2,N,M6,M5,M4)+trcsa_RFLS(idsa_CaOH2,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaCO3,N,M6,M5,M4)=trcsa_XFLS(idsa_CaCO3,N,M6,M5,M4)+trcsa_RFLS(idsa_CaCO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaHCO3,N,M6,M5,M4)=trcsa_XFLS(idsa_CaHCO3,N,M6,M5,M4)+trcsa_RFLS(idsa_CaHCO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaSO4,N,M6,M5,M4)=trcsa_XFLS(idsa_CaSO4,N,M6,M5,M4)+trcsa_RFLS(idsa_CaSO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_MgOH2,N,M6,M5,M4)=trcsa_XFLS(idsa_MgOH2,N,M6,M5,M4)+trcsa_RFLS(idsa_MgOH2,N,M6,M5,M4)
  trcsa_XFLS(idsa_MgCO3,N,M6,M5,M4)=trcsa_XFLS(idsa_MgCO3,N,M6,M5,M4)+trcsa_RFLS(idsa_MgCO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_MgHCO3,N,M6,M5,M4)=trcsa_XFLS(idsa_MgHCO3,N,M6,M5,M4)+trcsa_RFLS(idsa_MgHCO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_MgSO4,N,M6,M5,M4)=trcsa_XFLS(idsa_MgSO4,N,M6,M5,M4)+trcsa_RFLS(idsa_MgSO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_NaCO3,N,M6,M5,M4)=trcsa_XFLS(idsa_NaCO3,N,M6,M5,M4)+trcsa_RFLS(idsa_NaCO3,N,M6,M5,M4)
  trcsa_XFLS(idsa_NaSO4,N,M6,M5,M4)=trcsa_XFLS(idsa_NaSO4,N,M6,M5,M4)+trcsa_RFLS(idsa_NaSO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_KSO4,N,M6,M5,M4)=trcsa_XFLS(idsa_KSO4,N,M6,M5,M4)+trcsa_RFLS(idsa_KSO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_H0PO4,N,M6,M5,M4)=trcsa_XFLS(idsa_H0PO4,N,M6,M5,M4)+trcsa_RFLS(idsa_H0PO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_H3PO4,N,M6,M5,M4)=trcsa_XFLS(idsa_H3PO4,N,M6,M5,M4)+trcsa_RFLS(idsa_H3PO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeHPO4,N,M6,M5,M4)=trcsa_XFLS(idsa_FeHPO4,N,M6,M5,M4)+trcsa_RFLS(idsa_FeHPO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeH2PO4,N,M6,M5,M4)=trcsa_XFLS(idsa_FeH2PO4,N,M6,M5,M4)+trcsa_RFLS(idsa_FeH2PO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaPO4,N,M6,M5,M4)=trcsa_XFLS(idsa_CaPO4,N,M6,M5,M4)+trcsa_RFLS(idsa_CaPO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaHPO4,N,M6,M5,M4)=trcsa_XFLS(idsa_CaHPO4,N,M6,M5,M4)+trcsa_RFLS(idsa_CaHPO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaH2PO4,N,M6,M5,M4)=trcsa_XFLS(idsa_CaH2PO4,N,M6,M5,M4)+trcsa_RFLS(idsa_CaH2PO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_MgHPO4,N,M6,M5,M4)=trcsa_XFLS(idsa_MgHPO4,N,M6,M5,M4)+trcsa_RFLS(idsa_MgHPO4,N,M6,M5,M4)
  trcsa_XFLS(idsa_H0PO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_H0PO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_H0PO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_H3PO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_H3PO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_H3PO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeHPO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_FeHPO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_FeHPO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_FeH2PO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_FeH2PO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_FeH2PO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaPO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_CaPO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_CaPO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaHPO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_CaHPO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_CaHPO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_CaH2PO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_CaH2PO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_CaH2PO4B,N,M6,M5,M4)
  trcsa_XFLS(idsa_MgHPO4B,N,M6,M5,M4)=trcsa_XFLS(idsa_MgHPO4B,N,M6,M5,M4)+trcsa_RFLS(idsa_MgHPO4B,N,M6,M5,M4)
  trcsa_XFHS(idsa_Al,N,M6,M5,M4)=trcsa_XFHS(idsa_Al,N,M6,M5,M4)+RALFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_Fe,N,M6,M5,M4)=trcsa_XFHS(idsa_Fe,N,M6,M5,M4)+RFEFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_Hp,N,M6,M5,M4)=trcsa_XFHS(idsa_Hp,N,M6,M5,M4)+RHYFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_Ca,N,M6,M5,M4)=trcsa_XFHS(idsa_Ca,N,M6,M5,M4)+RCAFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_Mg,N,M6,M5,M4)=trcsa_XFHS(idsa_Mg,N,M6,M5,M4)+RMGFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_Na,N,M6,M5,M4)=trcsa_XFHS(idsa_Na,N,M6,M5,M4)+RNAFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_K,N,M6,M5,M4)=trcsa_XFHS(idsa_K,N,M6,M5,M4)+RKAFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_OH,N,M6,M5,M4)=trcsa_XFHS(idsa_OH,N,M6,M5,M4)+ROHFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_SO4,N,M6,M5,M4)=trcsa_XFHS(idsa_SO4,N,M6,M5,M4)+RSOFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_Cl,N,M6,M5,M4)=trcsa_XFHS(idsa_Cl,N,M6,M5,M4)+RCLFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CO3,N,M6,M5,M4)=trcsa_XFHS(idsa_CO3,N,M6,M5,M4)+RC3FHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_HCO3,N,M6,M5,M4)=trcsa_XFHS(idsa_HCO3,N,M6,M5,M4)+RHCFHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_AlOH,N,M6,M5,M4)=trcsa_XFHS(idsa_AlOH,N,M6,M5,M4)+RAL1HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_AlOH2,N,M6,M5,M4)=trcsa_XFHS(idsa_AlOH2,N,M6,M5,M4)+RAL2HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_AlOH3,N,M6,M5,M4)=trcsa_XFHS(idsa_AlOH3,N,M6,M5,M4)+RAL3HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_AlOH4,N,M6,M5,M4)=trcsa_XFHS(idsa_AlOH4,N,M6,M5,M4)+RAL4HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_AlSO4,N,M6,M5,M4)=trcsa_XFHS(idsa_AlSO4,N,M6,M5,M4)+RALSHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeOH,N,M6,M5,M4)=trcsa_XFHS(idsa_FeOH,N,M6,M5,M4)+RFE1HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeOH2,N,M6,M5,M4)=trcsa_XFHS(idsa_FeOH2,N,M6,M5,M4)+RFE2HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeOH3,N,M6,M5,M4)=trcsa_XFHS(idsa_FeOH3,N,M6,M5,M4)+RFE3HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeOH4,N,M6,M5,M4)=trcsa_XFHS(idsa_FeOH4,N,M6,M5,M4)+RFE4HS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeSO4,N,M6,M5,M4)=trcsa_XFHS(idsa_FeSO4,N,M6,M5,M4)+RFESHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaOH2,N,M6,M5,M4)=trcsa_XFHS(idsa_CaOH2,N,M6,M5,M4)+RCAOHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaCO3,N,M6,M5,M4)=trcsa_XFHS(idsa_CaCO3,N,M6,M5,M4)+RCACHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaHCO3,N,M6,M5,M4)=trcsa_XFHS(idsa_CaHCO3,N,M6,M5,M4)+RCAHHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaSO4,N,M6,M5,M4)=trcsa_XFHS(idsa_CaSO4,N,M6,M5,M4)+RCASHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_MgOH2,N,M6,M5,M4)=trcsa_XFHS(idsa_MgOH2,N,M6,M5,M4)+RMGOHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_MgCO3,N,M6,M5,M4)=trcsa_XFHS(idsa_MgCO3,N,M6,M5,M4)+RMGCHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_MgHCO3,N,M6,M5,M4)=trcsa_XFHS(idsa_MgHCO3,N,M6,M5,M4)+RMGHHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_MgSO4,N,M6,M5,M4)=trcsa_XFHS(idsa_MgSO4,N,M6,M5,M4)+RMGSHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_NaCO3,N,M6,M5,M4)=trcsa_XFHS(idsa_NaCO3,N,M6,M5,M4)+RNACHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_NaSO4,N,M6,M5,M4)=trcsa_XFHS(idsa_NaSO4,N,M6,M5,M4)+RNASHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_KSO4,N,M6,M5,M4)=trcsa_XFHS(idsa_KSO4,N,M6,M5,M4)+RKASHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_H0PO4,N,M6,M5,M4)=trcsa_XFHS(idsa_H0PO4,N,M6,M5,M4)+RH0PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_H3PO4,N,M6,M5,M4)=trcsa_XFHS(idsa_H3PO4,N,M6,M5,M4)+RH3PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeHPO4,N,M6,M5,M4)=trcsa_XFHS(idsa_FeHPO4,N,M6,M5,M4)+RF1PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeH2PO4,N,M6,M5,M4)=trcsa_XFHS(idsa_FeH2PO4,N,M6,M5,M4)+RF2PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaPO4,N,M6,M5,M4)=trcsa_XFHS(idsa_CaPO4,N,M6,M5,M4)+RC0PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaHPO4,N,M6,M5,M4)=trcsa_XFHS(idsa_CaHPO4,N,M6,M5,M4)+RC1PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaH2PO4,N,M6,M5,M4)=trcsa_XFHS(idsa_CaH2PO4,N,M6,M5,M4)+RC2PHS(N,M6,M5,M4)
  trcsa_XFHS(idsa_MgHPO4,N,M6,M5,M4)=trcsa_XFHS(idsa_MgHPO4,N,M6,M5,M4)+RM1PHS(N,M6,M5,M4)

  trcsa_XFHS(idsa_H0PO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_H0PO4B,N,M6,M5,M4)+RH0BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_H3PO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_H3PO4B,N,M6,M5,M4)+RH3BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeHPO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_FeHPO4B,N,M6,M5,M4)+RF1BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_FeH2PO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_FeH2PO4B,N,M6,M5,M4)+RF2BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaPO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_CaPO4B,N,M6,M5,M4)+RC0BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaHPO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_CaHPO4B,N,M6,M5,M4)+RC1BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_CaH2PO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_CaH2PO4B,N,M6,M5,M4)+RC2BHB(N,M6,M5,M4)
  trcsa_XFHS(idsa_MgHPO4B,N,M6,M5,M4)=trcsa_XFHS(idsa_MgHPO4B,N,M6,M5,M4)+RM1BHB(N,M6,M5,M4)
  end subroutine AccumFluxMacMicPores
!------------------------------------------------------------------------------------------

      subroutine NetOverloadFluxInWater(M,N,N1,N2,N4,N5,N4B,N5B)
!
!     Description:
!
      integer, intent(in) :: M,N,N1,N2,N4,N5,N4B,N5B
      integer :: NN
!     begin_execution
!
!     TQR*=net overland solute flux
!     RQR*=overland solute flux
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
  D1202: DO NN=1,2
    trcsa_TQR(idsa_Al,N2,N1)=trcsa_TQR(idsa_Al,N2,N1)+RQRAL(N,NN,N2,N1)
    trcsa_TQR(idsa_Fe,N2,N1)=trcsa_TQR(idsa_Fe,N2,N1)+RQRFE(N,NN,N2,N1)
    trcsa_TQR(idsa_Hp,N2,N1)=trcsa_TQR(idsa_Hp,N2,N1)+RQRHY(N,NN,N2,N1)
    trcsa_TQR(idsa_Ca,N2,N1)=trcsa_TQR(idsa_Ca,N2,N1)+RQRCA(N,NN,N2,N1)
    trcsa_TQR(idsa_Mg,N2,N1)=trcsa_TQR(idsa_Mg,N2,N1)+RQRMG(N,NN,N2,N1)
    trcsa_TQR(idsa_Na,N2,N1)=trcsa_TQR(idsa_Na,N2,N1)+RQRNA(N,NN,N2,N1)
    trcsa_TQR(idsa_K,N2,N1)=trcsa_TQR(idsa_K,N2,N1)+RQRKA(N,NN,N2,N1)
    trcsa_TQR(idsa_OH,N2,N1)=trcsa_TQR(idsa_OH,N2,N1)+RQROH(N,NN,N2,N1)
    trcsa_TQR(idsa_SO4,N2,N1)=trcsa_TQR(idsa_SO4,N2,N1)+RQRSO(N,NN,N2,N1)
    trcsa_TQR(idsa_Cl,N2,N1)=trcsa_TQR(idsa_Cl,N2,N1)+RQRCL(N,NN,N2,N1)
    trcsa_TQR(idsa_CO3,N2,N1)=trcsa_TQR(idsa_CO3,N2,N1)+RQRC3(N,NN,N2,N1)
    trcsa_TQR(idsa_HCO3,N2,N1)=trcsa_TQR(idsa_HCO3,N2,N1)+RQRHC(N,NN,N2,N1)
    trcsa_TQR(idsa_AlOH,N2,N1)=trcsa_TQR(idsa_AlOH,N2,N1)+RQRAL1(N,NN,N2,N1)
    trcsa_TQR(idsa_AlOH2,N2,N1)=trcsa_TQR(idsa_AlOH2,N2,N1)+RQRAL2(N,NN,N2,N1)
    trcsa_TQR(idsa_AlOH3,N2,N1)=trcsa_TQR(idsa_AlOH3,N2,N1)+RQRAL3(N,NN,N2,N1)
    trcsa_TQR(idsa_AlOH4,N2,N1)=trcsa_TQR(idsa_AlOH4,N2,N1)+RQRAL4(N,NN,N2,N1)
    trcsa_TQR(idsa_AlSO4,N2,N1)=trcsa_TQR(idsa_AlSO4,N2,N1)+RQRALS(N,NN,N2,N1)
    trcsa_TQR(idsa_FeOH,N2,N1)=trcsa_TQR(idsa_FeOH,N2,N1)+RQRFE1(N,NN,N2,N1)
    trcsa_TQR(idsa_FeOH2,N2,N1)=trcsa_TQR(idsa_FeOH2,N2,N1)+RQRFE2(N,NN,N2,N1)
    trcsa_TQR(idsa_FeOH3,N2,N1)=trcsa_TQR(idsa_FeOH3,N2,N1)+RQRFE3(N,NN,N2,N1)
    trcsa_TQR(idsa_FeOH4,N2,N1)=trcsa_TQR(idsa_FeOH4,N2,N1)+RQRFE4(N,NN,N2,N1)
    trcsa_TQR(idsa_FeSO4,N2,N1)=trcsa_TQR(idsa_FeSO4,N2,N1)+RQRFES(N,NN,N2,N1)
    trcsa_TQR(idsa_CaOH2,N2,N1)=trcsa_TQR(idsa_CaOH2,N2,N1)+RQRCAO(N,NN,N2,N1)
    trcsa_TQR(idsa_CaCO3,N2,N1)=trcsa_TQR(idsa_CaCO3,N2,N1)+RQRCAC(N,NN,N2,N1)
    trcsa_TQR(idsa_CaHCO3,N2,N1)=trcsa_TQR(idsa_CaHCO3,N2,N1)+RQRCAH(N,NN,N2,N1)
    trcsa_TQR(idsa_CaSO4,N2,N1)=trcsa_TQR(idsa_CaSO4,N2,N1)+RQRCAS(N,NN,N2,N1)
    trcsa_TQR(idsa_MgOH2,N2,N1)=trcsa_TQR(idsa_MgOH2,N2,N1)+RQRMGO(N,NN,N2,N1)
    trcsa_TQR(idsa_MgCO3,N2,N1)=trcsa_TQR(idsa_MgCO3,N2,N1)+RQRMGC(N,NN,N2,N1)
    trcsa_TQR(idsa_MgHCO3,N2,N1)=trcsa_TQR(idsa_MgHCO3,N2,N1)+RQRMGH(N,NN,N2,N1)
    trcsa_TQR(idsa_MgSO4,N2,N1)=trcsa_TQR(idsa_MgSO4,N2,N1)+RQRMGS(N,NN,N2,N1)
    trcsa_TQR(idsa_NaCO3,N2,N1)=trcsa_TQR(idsa_NaCO3,N2,N1)+RQRNAC(N,NN,N2,N1)
    trcsa_TQR(idsa_NaSO4,N2,N1)=trcsa_TQR(idsa_NaSO4,N2,N1)+RQRNAS(N,NN,N2,N1)
    trcsa_TQR(idsa_KSO4,N2,N1)=trcsa_TQR(idsa_KSO4,N2,N1)+RQRKAS(N,NN,N2,N1)
    trcsa_TQR(idsa_H0PO4,N2,N1)=trcsa_TQR(idsa_H0PO4,N2,N1)+RQRH0P(N,NN,N2,N1)
    trcsa_TQR(idsa_H3PO4,N2,N1)=trcsa_TQR(idsa_H3PO4,N2,N1)+RQRH3P(N,NN,N2,N1)
    trcsa_TQR(idsa_FeHPO4,N2,N1)=trcsa_TQR(idsa_FeHPO4,N2,N1)+RQRF1P(N,NN,N2,N1)
    trcsa_TQR(idsa_FeH2PO4,N2,N1)=trcsa_TQR(idsa_FeH2PO4,N2,N1)+RQRF2P(N,NN,N2,N1)
    trcsa_TQR(idsa_CaPO4,N2,N1)=trcsa_TQR(idsa_CaPO4,N2,N1)+RQRC0P(N,NN,N2,N1)
    trcsa_TQR(idsa_CaHPO4,N2,N1)=trcsa_TQR(idsa_CaHPO4,N2,N1)+RQRC1P(N,NN,N2,N1)
    trcsa_TQR(idsa_CaH2PO4,N2,N1)=trcsa_TQR(idsa_CaH2PO4,N2,N1)+RQRC2P(N,NN,N2,N1)
    trcsa_TQR(idsa_MgHPO4,N2,N1)=trcsa_TQR(idsa_MgHPO4,N2,N1)+RQRM1P(N,NN,N2,N1)
    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      trcsa_TQR(idsa_Al,N2,N1)=trcsa_TQR(idsa_Al,N2,N1)-RQRAL(N,NN,N5,N4)
      trcsa_TQR(idsa_Fe,N2,N1)=trcsa_TQR(idsa_Fe,N2,N1)-RQRFE(N,NN,N5,N4)
      trcsa_TQR(idsa_Hp,N2,N1)=trcsa_TQR(idsa_Hp,N2,N1)-RQRHY(N,NN,N5,N4)
      trcsa_TQR(idsa_Ca,N2,N1)=trcsa_TQR(idsa_Ca,N2,N1)-RQRCA(N,NN,N5,N4)
      trcsa_TQR(idsa_Mg,N2,N1)=trcsa_TQR(idsa_Mg,N2,N1)-RQRMG(N,NN,N5,N4)
      trcsa_TQR(idsa_Na,N2,N1)=trcsa_TQR(idsa_Na,N2,N1)-RQRNA(N,NN,N5,N4)
      trcsa_TQR(idsa_K,N2,N1)=trcsa_TQR(idsa_K,N2,N1)-RQRKA(N,NN,N5,N4)
      trcsa_TQR(idsa_OH,N2,N1)=trcsa_TQR(idsa_OH,N2,N1)-RQROH(N,NN,N5,N4)
      trcsa_TQR(idsa_SO4,N2,N1)=trcsa_TQR(idsa_SO4,N2,N1)-RQRSO(N,NN,N5,N4)
      trcsa_TQR(idsa_Cl,N2,N1)=trcsa_TQR(idsa_Cl,N2,N1)-RQRCL(N,NN,N5,N4)
      trcsa_TQR(idsa_CO3,N2,N1)=trcsa_TQR(idsa_CO3,N2,N1)-RQRC3(N,NN,N5,N4)
      trcsa_TQR(idsa_HCO3,N2,N1)=trcsa_TQR(idsa_HCO3,N2,N1)-RQRHC(N,NN,N5,N4)
      trcsa_TQR(idsa_AlOH,N2,N1)=trcsa_TQR(idsa_AlOH,N2,N1)-RQRAL1(N,NN,N5,N4)
      trcsa_TQR(idsa_AlOH2,N2,N1)=trcsa_TQR(idsa_AlOH2,N2,N1)-RQRAL2(N,NN,N5,N4)
      trcsa_TQR(idsa_AlOH3,N2,N1)=trcsa_TQR(idsa_AlOH3,N2,N1)-RQRAL3(N,NN,N5,N4)
      trcsa_TQR(idsa_AlOH4,N2,N1)=trcsa_TQR(idsa_AlOH4,N2,N1)-RQRAL4(N,NN,N5,N4)
      trcsa_TQR(idsa_AlSO4,N2,N1)=trcsa_TQR(idsa_AlSO4,N2,N1)-RQRALS(N,NN,N5,N4)
      trcsa_TQR(idsa_FeOH,N2,N1)=trcsa_TQR(idsa_FeOH,N2,N1)-RQRFE1(N,NN,N5,N4)
      trcsa_TQR(idsa_FeOH2,N2,N1)=trcsa_TQR(idsa_FeOH2,N2,N1)-RQRFE2(N,NN,N5,N4)
      trcsa_TQR(idsa_FeOH3,N2,N1)=trcsa_TQR(idsa_FeOH3,N2,N1)-RQRFE3(N,NN,N5,N4)
      trcsa_TQR(idsa_FeOH4,N2,N1)=trcsa_TQR(idsa_FeOH4,N2,N1)-RQRFE4(N,NN,N5,N4)
      trcsa_TQR(idsa_FeSO4,N2,N1)=trcsa_TQR(idsa_FeSO4,N2,N1)-RQRFES(N,NN,N5,N4)
      trcsa_TQR(idsa_CaOH2,N2,N1)=trcsa_TQR(idsa_CaOH2,N2,N1)-RQRCAO(N,NN,N5,N4)
      trcsa_TQR(idsa_CaCO3,N2,N1)=trcsa_TQR(idsa_CaCO3,N2,N1)-RQRCAC(N,NN,N5,N4)
      trcsa_TQR(idsa_CaHCO3,N2,N1)=trcsa_TQR(idsa_CaHCO3,N2,N1)-RQRCAH(N,NN,N5,N4)
      trcsa_TQR(idsa_CaSO4,N2,N1)=trcsa_TQR(idsa_CaSO4,N2,N1)-RQRCAS(N,NN,N5,N4)
      trcsa_TQR(idsa_MgOH2,N2,N1)=trcsa_TQR(idsa_MgOH2,N2,N1)-RQRMGO(N,NN,N5,N4)
      trcsa_TQR(idsa_MgCO3,N2,N1)=trcsa_TQR(idsa_MgCO3,N2,N1)-RQRMGC(N,NN,N5,N4)
      trcsa_TQR(idsa_MgHCO3,N2,N1)=trcsa_TQR(idsa_MgHCO3,N2,N1)-RQRMGH(N,NN,N5,N4)
      trcsa_TQR(idsa_MgSO4,N2,N1)=trcsa_TQR(idsa_MgSO4,N2,N1)-RQRMGS(N,NN,N5,N4)
      trcsa_TQR(idsa_NaCO3,N2,N1)=trcsa_TQR(idsa_NaCO3,N2,N1)-RQRNAC(N,NN,N5,N4)
      trcsa_TQR(idsa_NaSO4,N2,N1)=trcsa_TQR(idsa_NaSO4,N2,N1)-RQRNAS(N,NN,N5,N4)
      trcsa_TQR(idsa_KSO4,N2,N1)=trcsa_TQR(idsa_KSO4,N2,N1)-RQRKAS(N,NN,N5,N4)
      trcsa_TQR(idsa_H0PO4,N2,N1)=trcsa_TQR(idsa_H0PO4,N2,N1)-RQRH0P(N,NN,N5,N4)
      trcsa_TQR(idsa_H3PO4,N2,N1)=trcsa_TQR(idsa_H3PO4,N2,N1)-RQRH3P(N,NN,N5,N4)
      trcsa_TQR(idsa_FeHPO4,N2,N1)=trcsa_TQR(idsa_FeHPO4,N2,N1)-RQRF1P(N,NN,N5,N4)
      trcsa_TQR(idsa_FeH2PO4,N2,N1)=trcsa_TQR(idsa_FeH2PO4,N2,N1)-RQRF2P(N,NN,N5,N4)
      trcsa_TQR(idsa_CaPO4,N2,N1)=trcsa_TQR(idsa_CaPO4,N2,N1)-RQRC0P(N,NN,N5,N4)
      trcsa_TQR(idsa_CaHPO4,N2,N1)=trcsa_TQR(idsa_CaHPO4,N2,N1)-RQRC1P(N,NN,N5,N4)
      trcsa_TQR(idsa_CaH2PO4,N2,N1)=trcsa_TQR(idsa_CaH2PO4,N2,N1)-RQRC2P(N,NN,N5,N4)
      trcsa_TQR(idsa_MgHPO4,N2,N1)=trcsa_TQR(idsa_MgHPO4,N2,N1)-RQRM1P(N,NN,N5,N4)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      trcsa_TQR(idsa_Al,N2,N1)=trcsa_TQR(idsa_Al,N2,N1)-RQRAL(N,NN,N5B,N4B)
      trcsa_TQR(idsa_Fe,N2,N1)=trcsa_TQR(idsa_Fe,N2,N1)-RQRFE(N,NN,N5B,N4B)
      trcsa_TQR(idsa_Hp,N2,N1)=trcsa_TQR(idsa_Hp,N2,N1)-RQRHY(N,NN,N5B,N4B)
      trcsa_TQR(idsa_Ca,N2,N1)=trcsa_TQR(idsa_Ca,N2,N1)-RQRCA(N,NN,N5B,N4B)
      trcsa_TQR(idsa_Mg,N2,N1)=trcsa_TQR(idsa_Mg,N2,N1)-RQRMG(N,NN,N5B,N4B)
      trcsa_TQR(idsa_Na,N2,N1)=trcsa_TQR(idsa_Na,N2,N1)-RQRNA(N,NN,N5B,N4B)
      trcsa_TQR(idsa_K,N2,N1)=trcsa_TQR(idsa_K,N2,N1)-RQRKA(N,NN,N5B,N4B)
      trcsa_TQR(idsa_OH,N2,N1)=trcsa_TQR(idsa_OH,N2,N1)-RQROH(N,NN,N5B,N4B)
      trcsa_TQR(idsa_SO4,N2,N1)=trcsa_TQR(idsa_SO4,N2,N1)-RQRSO(N,NN,N5B,N4B)
      trcsa_TQR(idsa_Cl,N2,N1)=trcsa_TQR(idsa_Cl,N2,N1)-RQRCL(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CO3,N2,N1)=trcsa_TQR(idsa_CO3,N2,N1)-RQRC3(N,NN,N5B,N4B)
      trcsa_TQR(idsa_HCO3,N2,N1)=trcsa_TQR(idsa_HCO3,N2,N1)-RQRHC(N,NN,N5B,N4B)
      trcsa_TQR(idsa_AlOH,N2,N1)=trcsa_TQR(idsa_AlOH,N2,N1)-RQRAL1(N,NN,N5B,N4B)
      trcsa_TQR(idsa_AlOH2,N2,N1)=trcsa_TQR(idsa_AlOH2,N2,N1)-RQRAL2(N,NN,N5B,N4B)
      trcsa_TQR(idsa_AlOH3,N2,N1)=trcsa_TQR(idsa_AlOH3,N2,N1)-RQRAL3(N,NN,N5B,N4B)
      trcsa_TQR(idsa_AlOH4,N2,N1)=trcsa_TQR(idsa_AlOH4,N2,N1)-RQRAL4(N,NN,N5B,N4B)
      trcsa_TQR(idsa_AlSO4,N2,N1)=trcsa_TQR(idsa_AlSO4,N2,N1)-RQRALS(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeOH,N2,N1)=trcsa_TQR(idsa_FeOH,N2,N1)-RQRFE1(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeOH2,N2,N1)=trcsa_TQR(idsa_FeOH2,N2,N1)-RQRFE2(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeOH3,N2,N1)=trcsa_TQR(idsa_FeOH3,N2,N1)-RQRFE3(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeOH4,N2,N1)=trcsa_TQR(idsa_FeOH4,N2,N1)-RQRFE4(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeSO4,N2,N1)=trcsa_TQR(idsa_FeSO4,N2,N1)-RQRFES(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaOH2,N2,N1)=trcsa_TQR(idsa_CaOH2,N2,N1)-RQRCAO(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaCO3,N2,N1)=trcsa_TQR(idsa_CaCO3,N2,N1)-RQRCAC(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaHCO3,N2,N1)=trcsa_TQR(idsa_CaHCO3,N2,N1)-RQRCAH(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaSO4,N2,N1)=trcsa_TQR(idsa_CaSO4,N2,N1)-RQRCAS(N,NN,N5B,N4B)
      trcsa_TQR(idsa_MgOH2,N2,N1)=trcsa_TQR(idsa_MgOH2,N2,N1)-RQRMGO(N,NN,N5B,N4B)
      trcsa_TQR(idsa_MgCO3,N2,N1)=trcsa_TQR(idsa_MgCO3,N2,N1)-RQRMGC(N,NN,N5B,N4B)
      trcsa_TQR(idsa_MgHCO3,N2,N1)=trcsa_TQR(idsa_MgHCO3,N2,N1)-RQRMGH(N,NN,N5B,N4B)
      trcsa_TQR(idsa_MgSO4,N2,N1)=trcsa_TQR(idsa_MgSO4,N2,N1)-RQRMGS(N,NN,N5B,N4B)
      trcsa_TQR(idsa_NaCO3,N2,N1)=trcsa_TQR(idsa_NaCO3,N2,N1)-RQRNAC(N,NN,N5B,N4B)
      trcsa_TQR(idsa_NaSO4,N2,N1)=trcsa_TQR(idsa_NaSO4,N2,N1)-RQRNAS(N,NN,N5B,N4B)
      trcsa_TQR(idsa_KSO4,N2,N1)=trcsa_TQR(idsa_KSO4,N2,N1)-RQRKAS(N,NN,N5B,N4B)
      trcsa_TQR(idsa_H0PO4,N2,N1)=trcsa_TQR(idsa_H0PO4,N2,N1)-RQRH0P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_H3PO4,N2,N1)=trcsa_TQR(idsa_H3PO4,N2,N1)-RQRH3P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeHPO4,N2,N1)=trcsa_TQR(idsa_FeHPO4,N2,N1)-RQRF1P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_FeH2PO4,N2,N1)=trcsa_TQR(idsa_FeH2PO4,N2,N1)-RQRF2P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaPO4,N2,N1)=trcsa_TQR(idsa_CaPO4,N2,N1)-RQRC0P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaHPO4,N2,N1)=trcsa_TQR(idsa_CaHPO4,N2,N1)-RQRC1P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_CaH2PO4,N2,N1)=trcsa_TQR(idsa_CaH2PO4,N2,N1)-RQRC2P(N,NN,N5B,N4B)
      trcsa_TQR(idsa_MgHPO4,N2,N1)=trcsa_TQR(idsa_MgHPO4,N2,N1)-RQRM1P(N,NN,N5B,N4B)
    ENDIF
  ENDDO D1202
  end subroutine NetOverloadFluxInWater
!------------------------------------------------------------------------------------------

    subroutine NetOverloadFLuxInSnow(N,N1,N2,N4,N5)
!
!     Description:
!
    implicit none
    integer, intent(in) :: N,N1,N2,N4,N5
!     begin_execution
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
    trcsa_TQ(idsa_Al,N2,N1)=trcsa_TQ(idsa_Al,N2,N1)+trcsa_RQ(idsa_Al,N,N2,N1)-trcsa_RQ(idsa_Al,N,N5,N4)
    trcsa_TQ(idsa_Fe,N2,N1)=trcsa_TQ(idsa_Fe,N2,N1)+trcsa_RQ(idsa_Fe,N,N2,N1)-trcsa_RQ(idsa_Fe,N,N5,N4)
    trcsa_TQ(idsa_Hp,N2,N1)=trcsa_TQ(idsa_Hp,N2,N1)+trcsa_RQ(idsa_Hp,N,N2,N1)-trcsa_RQ(idsa_Hp,N,N5,N4)
    trcsa_TQ(idsa_Ca,N2,N1)=trcsa_TQ(idsa_Ca,N2,N1)+trcsa_RQ(idsa_Ca,N,N2,N1)-trcsa_RQ(idsa_Ca,N,N5,N4)
    trcsa_TQ(idsa_Mg,N2,N1)=trcsa_TQ(idsa_Mg,N2,N1)+trcsa_RQ(idsa_Mg,N,N2,N1)-trcsa_RQ(idsa_Mg,N,N5,N4)
    trcsa_TQ(idsa_Na,N2,N1)=trcsa_TQ(idsa_Na,N2,N1)+trcsa_RQ(idsa_Na,N,N2,N1)-trcsa_RQ(idsa_Na,N,N5,N4)
    trcsa_TQ(idsa_K,N2,N1)=trcsa_TQ(idsa_K,N2,N1)+trcsa_RQ(idsa_K,N,N2,N1)-trcsa_RQ(idsa_K,N,N5,N4)
    trcsa_TQ(idsa_OH,N2,N1)=trcsa_TQ(idsa_OH,N2,N1)+trcsa_RQ(idsa_OH,N,N2,N1)-trcsa_RQ(idsa_OH,N,N5,N4)
    trcsa_TQ(idsa_SO4,N2,N1)=trcsa_TQ(idsa_SO4,N2,N1)+trcsa_RQ(idsa_SO4,N,N2,N1)-trcsa_RQ(idsa_SO4,N,N5,N4)
    trcsa_TQ(idsa_Cl,N2,N1)=trcsa_TQ(idsa_Cl,N2,N1)+trcsa_RQ(idsa_Cl,N,N2,N1)-trcsa_RQ(idsa_Cl,N,N5,N4)
    trcsa_TQ(idsa_CO3,N2,N1)=trcsa_TQ(idsa_CO3,N2,N1)+trcsa_RQ(idsa_CO3,N,N2,N1)-trcsa_RQ(idsa_CO3,N,N5,N4)
    trcsa_TQ(idsa_HCO3,N2,N1)=trcsa_TQ(idsa_HCO3,N2,N1)+trcsa_RQ(idsa_HCO3,N,N2,N1)-trcsa_RQ(idsa_HCO3,N,N5,N4)
    trcsa_TQ(idsa_AlOH,N2,N1)=trcsa_TQ(idsa_AlOH,N2,N1)+trcsa_RQ(idsa_AlOH,N,N2,N1)-trcsa_RQ(idsa_AlOH,N,N5,N4)
    trcsa_TQ(idsa_AlOH2,N2,N1)=trcsa_TQ(idsa_AlOH2,N2,N1)+trcsa_RQ(idsa_AlOH2,N,N2,N1)-trcsa_RQ(idsa_AlOH2,N,N5,N4)
    trcsa_TQ(idsa_AlOH3,N2,N1)=trcsa_TQ(idsa_AlOH3,N2,N1)+trcsa_RQ(idsa_AlOH3,N,N2,N1)-trcsa_RQ(idsa_AlOH3,N,N5,N4)
    trcsa_TQ(idsa_AlOH4,N2,N1)=trcsa_TQ(idsa_AlOH4,N2,N1)+trcsa_RQ(idsa_AlOH4,N,N2,N1)-trcsa_RQ(idsa_AlOH4,N,N5,N4)
    trcsa_TQ(idsa_AlSO4,N2,N1)=trcsa_TQ(idsa_AlSO4,N2,N1)+trcsa_RQ(idsa_AlSO4,N,N2,N1)-trcsa_RQ(idsa_AlSO4,N,N5,N4)
    trcsa_TQ(idsa_FeOH,N2,N1)=trcsa_TQ(idsa_FeOH,N2,N1)+trcsa_RQ(idsa_FeOH,N,N2,N1)-trcsa_RQ(idsa_FeOH,N,N5,N4)
    trcsa_TQ(idsa_FeOH2,N2,N1)=trcsa_TQ(idsa_FeOH2,N2,N1)+trcsa_RQ(idsa_FeOH2,N,N2,N1)-trcsa_RQ(idsa_FeOH2,N,N5,N4)
    trcsa_TQ(idsa_FeOH3,N2,N1)=trcsa_TQ(idsa_FeOH3,N2,N1)+trcsa_RQ(idsa_FeOH3,N,N2,N1)-trcsa_RQ(idsa_FeOH3,N,N5,N4)
    trcsa_TQ(idsa_FeOH4,N2,N1)=trcsa_TQ(idsa_FeOH4,N2,N1)+trcsa_RQ(idsa_FeOH4,N,N2,N1)-trcsa_RQ(idsa_FeOH4,N,N5,N4)
    trcsa_TQ(idsa_FeSO4,N2,N1)=trcsa_TQ(idsa_FeSO4,N2,N1)+trcsa_RQ(idsa_FeSO4,N,N2,N1)-trcsa_RQ(idsa_FeSO4,N,N5,N4)
    trcsa_TQ(idsa_CaOH2,N2,N1)=trcsa_TQ(idsa_CaOH2,N2,N1)+trcsa_RQ(idsa_CaOH2,N,N2,N1)-trcsa_RQ(idsa_CaOH2,N,N5,N4)
    trcsa_TQ(idsa_CaCO3,N2,N1)=trcsa_TQ(idsa_CaCO3,N2,N1)+trcsa_RQ(idsa_CaCO3,N,N2,N1)-trcsa_RQ(idsa_CaCO3,N,N5,N4)
    trcsa_TQ(idsa_CaHCO3,N2,N1)=trcsa_TQ(idsa_CaHCO3,N2,N1)+trcsa_RQ(idsa_CaHCO3,N,N2,N1)-trcsa_RQ(idsa_CaHCO3,N,N5,N4)
    trcsa_TQ(idsa_CaSO4,N2,N1)=trcsa_TQ(idsa_CaSO4,N2,N1)+trcsa_RQ(idsa_CaSO4,N,N2,N1)-trcsa_RQ(idsa_CaSO4,N,N5,N4)
    trcsa_TQ(idsa_MgOH2,N2,N1)=trcsa_TQ(idsa_MgOH2,N2,N1)+trcsa_RQ(idsa_MgOH2,N,N2,N1)-trcsa_RQ(idsa_MgOH2,N,N5,N4)
    trcsa_TQ(idsa_MgCO3,N2,N1)=trcsa_TQ(idsa_MgCO3,N2,N1)+trcsa_RQ(idsa_MgCO3,N,N2,N1)-trcsa_RQ(idsa_MgCO3,N,N5,N4)
    trcsa_TQ(idsa_MgHCO3,N2,N1)=trcsa_TQ(idsa_MgHCO3,N2,N1)+trcsa_RQ(idsa_MgHCO3,N,N2,N1)-trcsa_RQ(idsa_MgHCO3,N,N5,N4)
    trcsa_TQ(idsa_MgSO4,N2,N1)=trcsa_TQ(idsa_MgSO4,N2,N1)+trcsa_RQ(idsa_MgSO4,N,N2,N1)-trcsa_RQ(idsa_MgSO4,N,N5,N4)
    trcsa_TQ(idsa_NaCO3,N2,N1)=trcsa_TQ(idsa_NaCO3,N2,N1)+trcsa_RQ(idsa_NaCO3,N,N2,N1)-trcsa_RQ(idsa_NaCO3,N,N5,N4)
    trcsa_TQ(idsa_NaSO4,N2,N1)=trcsa_TQ(idsa_NaSO4,N2,N1)+trcsa_RQ(idsa_NaSO4,N,N2,N1)-trcsa_RQ(idsa_NaSO4,N,N5,N4)
    trcsa_TQ(idsa_KSO4,N2,N1)=trcsa_TQ(idsa_KSO4,N2,N1)+trcsa_RQ(idsa_KSO4,N,N2,N1)-trcsa_RQ(idsa_KSO4,N,N5,N4)
    trcsa_TQ(idsa_H0PO4,N2,N1)=trcsa_TQ(idsa_H0PO4,N2,N1)+trcsa_RQ(idsa_H0PO4,N,N2,N1)-trcsa_RQ(idsa_H0PO4,N,N5,N4)
    trcsa_TQ(idsa_H3PO4,N2,N1)=trcsa_TQ(idsa_H3PO4,N2,N1)+trcsa_RQ(idsa_H3PO4,N,N2,N1)-trcsa_RQ(idsa_H3PO4,N,N5,N4)
    trcsa_TQ(idsa_FeHPO4,N2,N1)=trcsa_TQ(idsa_FeHPO4,N2,N1)+trcsa_RQ(idsa_FeHPO4,N,N2,N1)-trcsa_RQ(idsa_FeHPO4,N,N5,N4)
    trcsa_TQ(idsa_FeH2PO4,N2,N1)=trcsa_TQ(idsa_FeH2PO4,N2,N1)+trcsa_RQ(idsa_FeH2PO4,N,N2,N1)-trcsa_RQ(idsa_FeH2PO4,N,N5,N4)
    trcsa_TQ(idsa_CaPO4,N2,N1)=trcsa_TQ(idsa_CaPO4,N2,N1)+trcsa_RQ(idsa_CaPO4,N,N2,N1)-trcsa_RQ(idsa_CaPO4,N,N5,N4)
    trcsa_TQ(idsa_CaHPO4,N2,N1)=trcsa_TQ(idsa_CaHPO4,N2,N1)+trcsa_RQ(idsa_CaHPO4,N,N2,N1)-trcsa_RQ(idsa_CaHPO4,N,N5,N4)
    trcsa_TQ(idsa_CaH2PO4,N2,N1)=trcsa_TQ(idsa_CaH2PO4,N2,N1)+trcsa_RQ(idsa_CaH2PO4,N,N2,N1)-trcsa_RQ(idsa_CaH2PO4,N,N5,N4)
    trcsa_TQ(idsa_MgHPO4,N2,N1)=trcsa_TQ(idsa_MgHPO4,N2,N1)+trcsa_RQ(idsa_MgHPO4,N,N2,N1)-trcsa_RQ(idsa_MgHPO4,N,N5,N4)
    end subroutine NetOverloadFLuxInSnow
!------------------------------------------------------------------------------------------

  subroutine NetFluxInSnowpack(M,NY,NX,N1,N2)
!
!     Description:
!
      integer, intent(in) :: M,NY,NX,N1,N2
      integer :: LS,LS2
!     begin_execution

  D1205: DO LS=1,JS
    IF(VHCPWM(M,LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
  !
    IF(LS.LT.JS.AND.VHCPWM(M,LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
      trcsa_TBLS(idsa_Al,LS,NY,NX)=trcsa_TBLS(idsa_Al,LS,NY,NX)+RALBLS(LS,NY,NX)-RALBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_Fe,LS,NY,NX)=trcsa_TBLS(idsa_Fe,LS,NY,NX)+RFEBLS(LS,NY,NX)-RFEBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_Hp,LS,NY,NX)=trcsa_TBLS(idsa_Hp,LS,NY,NX)+RHYBLS(LS,NY,NX)-RHYBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_Ca,LS,NY,NX)=trcsa_TBLS(idsa_Ca,LS,NY,NX)+RCABLS(LS,NY,NX)-RCABLS(LS2,NY,NX)
      trcsa_TBLS(idsa_Mg,LS,NY,NX)=trcsa_TBLS(idsa_Mg,LS,NY,NX)+RMGBLS(LS,NY,NX)-RMGBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_Na,LS,NY,NX)=trcsa_TBLS(idsa_Na,LS,NY,NX)+RNABLS(LS,NY,NX)-RNABLS(LS2,NY,NX)
      trcsa_TBLS(idsa_K,LS,NY,NX)=trcsa_TBLS(idsa_K,LS,NY,NX)+RKABLS(LS,NY,NX)-RKABLS(LS2,NY,NX)
      trcsa_TBLS(idsa_OH,LS,NY,NX)=trcsa_TBLS(idsa_OH,LS,NY,NX)+ROHBLS(LS,NY,NX)-ROHBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_SO4,LS,NY,NX)=trcsa_TBLS(idsa_SO4,LS,NY,NX)+RSOBLS(LS,NY,NX)-RSOBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_Cl,LS,NY,NX)=trcsa_TBLS(idsa_Cl,LS,NY,NX)+RCLBLS(LS,NY,NX)-RCLBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_CO3,LS,NY,NX)=trcsa_TBLS(idsa_CO3,LS,NY,NX)+RC3BLS(LS,NY,NX)-RC3BLS(LS2,NY,NX)
      trcsa_TBLS(idsa_HCO3,LS,NY,NX)=trcsa_TBLS(idsa_HCO3,LS,NY,NX)+RHCBLS(LS,NY,NX)-RHCBLS(LS2,NY,NX)
      trcsa_TBLS(idsa_AlOH,LS,NY,NX)=trcsa_TBLS(idsa_AlOH,LS,NY,NX)+RAL1BS(LS,NY,NX)-RAL1BS(LS2,NY,NX)
      trcsa_TBLS(idsa_AlOH2,LS,NY,NX)=trcsa_TBLS(idsa_AlOH2,LS,NY,NX)+RAL2BS(LS,NY,NX)-RAL2BS(LS2,NY,NX)
      trcsa_TBLS(idsa_AlOH3,LS,NY,NX)=trcsa_TBLS(idsa_AlOH3,LS,NY,NX)+RAL3BS(LS,NY,NX)-RAL3BS(LS2,NY,NX)
      trcsa_TBLS(idsa_AlOH4,LS,NY,NX)=trcsa_TBLS(idsa_AlOH4,LS,NY,NX)+RAL4BS(LS,NY,NX)-RAL4BS(LS2,NY,NX)
      trcsa_TBLS(idsa_AlSO4,LS,NY,NX)=trcsa_TBLS(idsa_AlSO4,LS,NY,NX)+RALSBS(LS,NY,NX)-RALSBS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeOH,LS,NY,NX)=trcsa_TBLS(idsa_FeOH,LS,NY,NX)+RFE1BS(LS,NY,NX)-RFE1BS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeOH2,LS,NY,NX)=trcsa_TBLS(idsa_FeOH2,LS,NY,NX)+RFE2BS(LS,NY,NX)-RFE2BS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeOH3,LS,NY,NX)=trcsa_TBLS(idsa_FeOH3,LS,NY,NX)+RFE3BS(LS,NY,NX)-RFE3BS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeOH4,LS,NY,NX)=trcsa_TBLS(idsa_FeOH4,LS,NY,NX)+RFE4BS(LS,NY,NX)-RFE4BS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeSO4,LS,NY,NX)=trcsa_TBLS(idsa_FeSO4,LS,NY,NX)+RFESBS(LS,NY,NX)-RFESBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaOH2,LS,NY,NX)=trcsa_TBLS(idsa_CaOH2,LS,NY,NX)+RCAOBS(LS,NY,NX)-RCAOBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaCO3,LS,NY,NX)+RCACBS(LS,NY,NX)-RCACBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)+RCAHBS(LS,NY,NX)-RCAHBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaSO4,LS,NY,NX)=trcsa_TBLS(idsa_CaSO4,LS,NY,NX)+RCASBS(LS,NY,NX)-RCASBS(LS2,NY,NX)
      trcsa_TBLS(idsa_MgOH2,LS,NY,NX)=trcsa_TBLS(idsa_MgOH2,LS,NY,NX)+RMGOBS(LS,NY,NX)-RMGOBS(LS2,NY,NX)
      trcsa_TBLS(idsa_MgCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgCO3,LS,NY,NX)+RMGCBS(LS,NY,NX)-RMGCBS(LS2,NY,NX)
      trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)+RMGHBS(LS,NY,NX)-RMGHBS(LS2,NY,NX)
      trcsa_TBLS(idsa_MgSO4,LS,NY,NX)=trcsa_TBLS(idsa_MgSO4,LS,NY,NX)+RMGSBS(LS,NY,NX)-RMGSBS(LS2,NY,NX)
      trcsa_TBLS(idsa_NaCO3,LS,NY,NX)=trcsa_TBLS(idsa_NaCO3,LS,NY,NX)+RNACBS(LS,NY,NX)-RNACBS(LS2,NY,NX)
      trcsa_TBLS(idsa_NaSO4,LS,NY,NX)=trcsa_TBLS(idsa_NaSO4,LS,NY,NX)+RNASBS(LS,NY,NX)-RNASBS(LS2,NY,NX)
      trcsa_TBLS(idsa_KSO4,LS,NY,NX)=trcsa_TBLS(idsa_KSO4,LS,NY,NX)+RKASBS(LS,NY,NX)-RKASBS(LS2,NY,NX)
      trcsa_TBLS(idsa_H0PO4,LS,NY,NX)=trcsa_TBLS(idsa_H0PO4,LS,NY,NX)+RH0PBS(LS,NY,NX)-RH0PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_H3PO4,LS,NY,NX)=trcsa_TBLS(idsa_H3PO4,LS,NY,NX)+RH3PBS(LS,NY,NX)-RH3PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)=trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)+RF1PBS(LS,NY,NX)-RF1PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)+RF2PBS(LS,NY,NX)-RF2PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaPO4,LS,NY,NX)+RC0PBS(LS,NY,NX)-RC0PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)+RC1PBS(LS,NY,NX)-RC1PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)+RC2PBS(LS,NY,NX)-RC2PBS(LS2,NY,NX)
      trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)=trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)+RM1PBS(LS,NY,NX)-RM1PBS(LS2,NY,NX)
    ELSE
  !
  !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
  !
      trcsa_TBLS(idsa_Al,LS,NY,NX)=trcsa_TBLS(idsa_Al,LS,NY,NX)+RALBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Al,3,0,N2,N1)-trcsa_RFLS(idsa_Al,3,NUM(N2,N1),N2,N1)
  !    3-RALFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Fe,LS,NY,NX)=trcsa_TBLS(idsa_Fe,LS,NY,NX)+RFEBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Fe,3,0,N2,N1)-trcsa_RFLS(idsa_Fe,3,NUM(N2,N1),N2,N1)
  !    3-RFEFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Hp,LS,NY,NX)=trcsa_TBLS(idsa_Hp,LS,NY,NX)+RHYBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Hp,3,0,N2,N1)-trcsa_RFLS(idsa_Hp,3,NUM(N2,N1),N2,N1)
  !    3-RHYFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Ca,LS,NY,NX)=trcsa_TBLS(idsa_Ca,LS,NY,NX)+RCABLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Ca,3,0,N2,N1)-trcsa_RFLS(idsa_Ca,3,NUM(N2,N1),N2,N1)
  !    3-RCAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Mg,LS,NY,NX)=trcsa_TBLS(idsa_Mg,LS,NY,NX)+RMGBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Mg,3,0,N2,N1)-trcsa_RFLS(idsa_Mg,3,NUM(N2,N1),N2,N1)
  !    3-RMGFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Na,LS,NY,NX)=trcsa_TBLS(idsa_Na,LS,NY,NX)+RNABLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Na,3,0,N2,N1)-trcsa_RFLS(idsa_Na,3,NUM(N2,N1),N2,N1)
  !    3-RNAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_K,LS,NY,NX)=trcsa_TBLS(idsa_K,LS,NY,NX)+RKABLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_K,3,0,N2,N1)-trcsa_RFLS(idsa_K,3,NUM(N2,N1),N2,N1)
  !    3-RKAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_OH,LS,NY,NX)=trcsa_TBLS(idsa_OH,LS,NY,NX)+ROHBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_OH,3,0,N2,N1)-trcsa_RFLS(idsa_OH,3,NUM(N2,N1),N2,N1)
  !    3-ROHFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_SO4,LS,NY,NX)=trcsa_TBLS(idsa_SO4,LS,NY,NX)+RSOBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_SO4,3,0,N2,N1)-trcsa_RFLS(idsa_SO4,3,NUM(N2,N1),N2,N1)
  !    3-RSOFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Cl,LS,NY,NX)=trcsa_TBLS(idsa_Cl,LS,NY,NX)+RCLBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_Cl,3,0,N2,N1)-trcsa_RFLS(idsa_Cl,3,NUM(N2,N1),N2,N1)
  !    3-RCLFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CO3,LS,NY,NX)=trcsa_TBLS(idsa_CO3,LS,NY,NX)+RC3BLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CO3,3,0,N2,N1)-trcsa_RFLS(idsa_CO3,3,NUM(N2,N1),N2,N1)
  !    3-RC3FHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_HCO3,LS,NY,NX)=trcsa_TBLS(idsa_HCO3,LS,NY,NX)+RHCBLS(LS,NY,NX) &
        -trcsa_RFLS(idsa_HCO3,3,0,N2,N1)-trcsa_RFLS(idsa_HCO3,3,NUM(N2,N1),N2,N1)
  !    3-RHCFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH,LS,NY,NX)=trcsa_TBLS(idsa_AlOH,LS,NY,NX)+RAL1BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH,3,NUM(N2,N1),N2,N1)
  !    3-RAL1HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH2,LS,NY,NX)=trcsa_TBLS(idsa_AlOH2,LS,NY,NX)+RAL2BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH2,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH2,3,NUM(N2,N1),N2,N1)
  !    3-RAL2HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH3,LS,NY,NX)=trcsa_TBLS(idsa_AlOH3,LS,NY,NX)+RAL3BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH3,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH3,3,NUM(N2,N1),N2,N1)
  !    3-RAL3HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH4,LS,NY,NX)=trcsa_TBLS(idsa_AlOH4,LS,NY,NX)+RAL4BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH4,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH4,3,NUM(N2,N1),N2,N1)
  !    3-RAL4HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlSO4,LS,NY,NX)=trcsa_TBLS(idsa_AlSO4,LS,NY,NX)+RALSBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_AlSO4,3,0,N2,N1)-trcsa_RFLS(idsa_AlSO4,3,NUM(N2,N1),N2,N1)
  !    3-RALSHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH,LS,NY,NX)=trcsa_TBLS(idsa_FeOH,LS,NY,NX)+RFE1BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH,3,NUM(N2,N1),N2,N1)
  !    3-RFE1HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH2,LS,NY,NX)=trcsa_TBLS(idsa_FeOH2,LS,NY,NX)+RFE2BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH2,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH2,3,NUM(N2,N1),N2,N1)
  !    3-RFE2HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH3,LS,NY,NX)=trcsa_TBLS(idsa_FeOH3,LS,NY,NX)+RFE3BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH3,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH3,3,NUM(N2,N1),N2,N1)
  !    3-RFE3HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH4,LS,NY,NX)=trcsa_TBLS(idsa_FeOH4,LS,NY,NX)+RFE4BS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH4,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH4,3,NUM(N2,N1),N2,N1)
  !    3-RFE4HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeSO4,LS,NY,NX)=trcsa_TBLS(idsa_FeSO4,LS,NY,NX)+RFESBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeSO4,3,0,N2,N1)-trcsa_RFLS(idsa_FeSO4,3,NUM(N2,N1),N2,N1)
  !    3-RFESHS(3,NUM(N2,N1),N2,N1)
        trcsa_TBLS(idsa_CaOH2,LS,NY,NX)=trcsa_TBLS(idsa_CaOH2,LS,NY,NX)+RCAOBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaOH2,3,0,N2,N1)-trcsa_RFLS(idsa_CaOH2,3,NUM(N2,N1),N2,N1)
  !    3-RCAOHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaCO3,LS,NY,NX)+RCACBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaCO3,3,0,N2,N1)-trcsa_RFLS(idsa_CaCO3,3,NUM(N2,N1),N2,N1)
  !    3-RCACHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)+RCAHBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaHCO3,3,0,N2,N1)-trcsa_RFLS(idsa_CaHCO3,3,NUM(N2,N1),N2,N1)
  !    3-RCAHHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaSO4,LS,NY,NX)=trcsa_TBLS(idsa_CaSO4,LS,NY,NX)+RCASBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaSO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaSO4,3,NUM(N2,N1),N2,N1)
  !    3-RCASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgOH2,LS,NY,NX)=trcsa_TBLS(idsa_MgOH2,LS,NY,NX)+RMGOBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_MgOH2,3,0,N2,N1)-trcsa_RFLS(idsa_MgOH2,3,NUM(N2,N1),N2,N1)
  !    3-RMGOHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgCO3,LS,NY,NX)+RMGCBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_MgCO3,3,0,N2,N1)-trcsa_RFLS(idsa_MgCO3,3,NUM(N2,N1),N2,N1)
  !    3-RMGCHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)+RMGHBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_MgHCO3,3,0,N2,N1)-trcsa_RFLS(idsa_MgHCO3,3,NUM(N2,N1),N2,N1)
  !    3-RMGHHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgSO4,LS,NY,NX)=trcsa_TBLS(idsa_MgSO4,LS,NY,NX)+RMGSBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_MgSO4,3,0,N2,N1)-trcsa_RFLS(idsa_MgSO4,3,NUM(N2,N1),N2,N1)
  !    3-RMGSHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_NaCO3,LS,NY,NX)=trcsa_TBLS(idsa_NaCO3,LS,NY,NX)+RNACBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_NaCO3,3,0,N2,N1)-trcsa_RFLS(idsa_NaCO3,3,NUM(N2,N1),N2,N1)
  !    3-RNACHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_NaSO4,LS,NY,NX)=trcsa_TBLS(idsa_NaSO4,LS,NY,NX)+RNASBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_NaSO4,3,0,N2,N1)-trcsa_RFLS(idsa_NaSO4,3,NUM(N2,N1),N2,N1)
  !    3-RNASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_KSO4,LS,NY,NX)=trcsa_TBLS(idsa_KSO4,LS,NY,NX)+RKASBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_KSO4,3,0,N2,N1)-trcsa_RFLS(idsa_KSO4,3,NUM(N2,N1),N2,N1)
  !    3-RKASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_H0PO4,LS,NY,NX)=trcsa_TBLS(idsa_H0PO4,LS,NY,NX)+RH0PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_H0PO4,3,0,N2,N1)-trcsa_RFLS(idsa_H0PO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_H0PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RH0PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_H0PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_H3PO4,LS,NY,NX)=trcsa_TBLS(idsa_H3PO4,LS,NY,NX)+RH3PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_H3PO4,3,0,N2,N1)-RH3PHS(3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_H3PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RH3PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_H3PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)=trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)+RF1PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeHPO4,3,0,N2,N1)-trcsa_RFLS(idsa_FeHPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_FeHPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RF1PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_FeHPO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)+RF2PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_FeH2PO4,3,0,N2,N1)-trcsa_RFLS(idsa_FeH2PO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_FeH2PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RF2PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_FeH2PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaPO4,LS,NY,NX)+RC0PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaPO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_CaPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RC0PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_CaPO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)+RC1PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaHPO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaHPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_CaHPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RC1PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_CaHPO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)+RC2PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_CaH2PO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaH2PO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_CaH2PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RC2PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_CaH2PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)=trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)+RM1PBS(LS,NY,NX) &
        -trcsa_RFLS(idsa_MgHPO4,3,0,N2,N1)-trcsa_RFLS(idsa_MgHPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_MgHPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RM1PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_MgHPO4B,3,NUM(N2,N1),N2,N1)
      ENDIF
    ENDIF
  ENDDO D1205
  end subroutine NetFluxInSnowpack
!------------------------------------------------------------------------------------------

  subroutine TotFluxInMacMicPores(N,N1,N2,N3,N4,N5,NY,NX,N6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,NY,NX
  integer, intent(inout) :: N6
  integer :: LL
!     begin_execution
!
!     T*FLS=net convective + diffusive solute flux through micropores
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     T*FHS=net convective + diffusive solute flux through macropores
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band

  IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
    D1200: DO LL=N6,NL(NY,NX)
      IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
        N6=LL
        exit
      ENDIF
    ENDDO D1200

    IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      TALFLS(N3,N2,N1)=TALFLS(N3,N2,N1)+trcsa_RFLS(idsa_Al,N,N3,N2,N1)-trcsa_RFLS(idsa_Al,N,N6,N5,N4)
      TFEFLS(N3,N2,N1)=TFEFLS(N3,N2,N1)+trcsa_RFLS(idsa_Fe,N,N3,N2,N1)-trcsa_RFLS(idsa_Fe,N,N6,N5,N4)
      THYFLS(N3,N2,N1)=THYFLS(N3,N2,N1)+trcsa_RFLS(idsa_Hp,N,N3,N2,N1)-trcsa_RFLS(idsa_Hp,N,N6,N5,N4)
      TCAFLS(N3,N2,N1)=TCAFLS(N3,N2,N1)+trcsa_RFLS(idsa_Ca,N,N3,N2,N1)-trcsa_RFLS(idsa_Ca,N,N6,N5,N4)
      TMGFLS(N3,N2,N1)=TMGFLS(N3,N2,N1)+trcsa_RFLS(idsa_Mg,N,N3,N2,N1)-trcsa_RFLS(idsa_Mg,N,N6,N5,N4)
      TNAFLS(N3,N2,N1)=TNAFLS(N3,N2,N1)+trcsa_RFLS(idsa_Na,N,N3,N2,N1)-trcsa_RFLS(idsa_Na,N,N6,N5,N4)
      TKAFLS(N3,N2,N1)=TKAFLS(N3,N2,N1)+trcsa_RFLS(idsa_K,N,N3,N2,N1)-trcsa_RFLS(idsa_K,N,N6,N5,N4)
      TOHFLS(N3,N2,N1)=TOHFLS(N3,N2,N1)+trcsa_RFLS(idsa_OH,N,N3,N2,N1)-trcsa_RFLS(idsa_OH,N,N6,N5,N4)
      TSOFLS(N3,N2,N1)=TSOFLS(N3,N2,N1)+trcsa_RFLS(idsa_SO4,N,N3,N2,N1)-trcsa_RFLS(idsa_SO4,N,N6,N5,N4)
      TCLFLS(N3,N2,N1)=TCLFLS(N3,N2,N1)+trcsa_RFLS(idsa_Cl,N,N3,N2,N1)-trcsa_RFLS(idsa_Cl,N,N6,N5,N4)
      TC3FLS(N3,N2,N1)=TC3FLS(N3,N2,N1)+trcsa_RFLS(idsa_CO3,N,N3,N2,N1)-trcsa_RFLS(idsa_CO3,N,N6,N5,N4)
      THCFLS(N3,N2,N1)=THCFLS(N3,N2,N1)+trcsa_RFLS(idsa_HCO3,N,N3,N2,N1)-trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)
      TAL1FS(N3,N2,N1)=TAL1FS(N3,N2,N1)+trcsa_RFLS(idsa_AlOH,N,N3,N2,N1)-trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)
      TAL2FS(N3,N2,N1)=TAL2FS(N3,N2,N1)+trcsa_RFLS(idsa_AlOH2,N,N3,N2,N1)-trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)
      TAL3FS(N3,N2,N1)=TAL3FS(N3,N2,N1)+trcsa_RFLS(idsa_AlOH3,N,N3,N2,N1)-trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)
      TAL4FS(N3,N2,N1)=TAL4FS(N3,N2,N1)+trcsa_RFLS(idsa_AlOH4,N,N3,N2,N1)-trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)
      TALSFS(N3,N2,N1)=TALSFS(N3,N2,N1)+trcsa_RFLS(idsa_AlSO4,N,N3,N2,N1)-trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)
      TFE1FS(N3,N2,N1)=TFE1FS(N3,N2,N1)+trcsa_RFLS(idsa_FeOH,N,N3,N2,N1)-trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)
      TFE2FS(N3,N2,N1)=TFE2FS(N3,N2,N1)+trcsa_RFLS(idsa_FeOH2,N,N3,N2,N1)-trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)
      TFE3FS(N3,N2,N1)=TFE3FS(N3,N2,N1)+trcsa_RFLS(idsa_FeOH3,N,N3,N2,N1)-trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)
      TFE4FS(N3,N2,N1)=TFE4FS(N3,N2,N1)+trcsa_RFLS(idsa_FeOH4,N,N3,N2,N1)-trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)
      TFESFS(N3,N2,N1)=TFESFS(N3,N2,N1)+trcsa_RFLS(idsa_FeSO4,N,N3,N2,N1)-trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)
      TCAOFS(N3,N2,N1)=TCAOFS(N3,N2,N1)+trcsa_RFLS(idsa_CaOH2,N,N3,N2,N1)-trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)
      TCACFS(N3,N2,N1)=TCACFS(N3,N2,N1)+trcsa_RFLS(idsa_CaCO3,N,N3,N2,N1)-trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)
      TCAHFS(N3,N2,N1)=TCAHFS(N3,N2,N1)+trcsa_RFLS(idsa_CaHCO3,N,N3,N2,N1)-trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)
      TCASFS(N3,N2,N1)=TCASFS(N3,N2,N1)+trcsa_RFLS(idsa_CaSO4,N,N3,N2,N1)-trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)
      TMGOFS(N3,N2,N1)=TMGOFS(N3,N2,N1)+trcsa_RFLS(idsa_MgOH2,N,N3,N2,N1)-trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)
      TMGCFS(N3,N2,N1)=TMGCFS(N3,N2,N1)+trcsa_RFLS(idsa_MgCO3,N,N3,N2,N1)-trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)
      TMGHFS(N3,N2,N1)=TMGHFS(N3,N2,N1)+trcsa_RFLS(idsa_MgHCO3,N,N3,N2,N1)-trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)
      TMGSFS(N3,N2,N1)=TMGSFS(N3,N2,N1)+trcsa_RFLS(idsa_MgSO4,N,N3,N2,N1)-trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)
      TNACFS(N3,N2,N1)=TNACFS(N3,N2,N1)+trcsa_RFLS(idsa_NaCO3,N,N3,N2,N1)-trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)
      TNASFS(N3,N2,N1)=TNASFS(N3,N2,N1)+trcsa_RFLS(idsa_NaSO4,N,N3,N2,N1)-trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)
      TKASFS(N3,N2,N1)=TKASFS(N3,N2,N1)+trcsa_RFLS(idsa_KSO4,N,N3,N2,N1)-trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)
      TH0PFS(N3,N2,N1)=TH0PFS(N3,N2,N1)+trcsa_RFLS(idsa_H0PO4,N,N3,N2,N1)-trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)
      TH3PFS(N3,N2,N1)=TH3PFS(N3,N2,N1)+trcsa_RFLS(idsa_H3PO4,N,N3,N2,N1)-trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)
      TF1PFS(N3,N2,N1)=TF1PFS(N3,N2,N1)+trcsa_RFLS(idsa_FeHPO4,N,N3,N2,N1)-trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)
      TF2PFS(N3,N2,N1)=TF2PFS(N3,N2,N1)+trcsa_RFLS(idsa_FeH2PO4,N,N3,N2,N1)-trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)
      TC0PFS(N3,N2,N1)=TC0PFS(N3,N2,N1)+trcsa_RFLS(idsa_CaPO4,N,N3,N2,N1)-trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)
      TC1PFS(N3,N2,N1)=TC1PFS(N3,N2,N1)+trcsa_RFLS(idsa_CaHPO4,N,N3,N2,N1)-trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)
      TC2PFS(N3,N2,N1)=TC2PFS(N3,N2,N1)+trcsa_RFLS(idsa_CaH2PO4,N,N3,N2,N1)-trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)
      TM1PFS(N3,N2,N1)=TM1PFS(N3,N2,N1)+trcsa_RFLS(idsa_MgHPO4,N,N3,N2,N1)-trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)
      TH0BFB(N3,N2,N1)=TH0BFB(N3,N2,N1)+trcsa_RFLS(idsa_H0PO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)
      TH3BFB(N3,N2,N1)=TH3BFB(N3,N2,N1)+trcsa_RFLS(idsa_H3PO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)
      TF1BFB(N3,N2,N1)=TF1BFB(N3,N2,N1)+trcsa_RFLS(idsa_FeHPO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)
      TF2BFB(N3,N2,N1)=TF2BFB(N3,N2,N1)+trcsa_RFLS(idsa_FeH2PO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)
      TC0BFB(N3,N2,N1)=TC0BFB(N3,N2,N1)+trcsa_RFLS(idsa_CaPO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)
      TC1BFB(N3,N2,N1)=TC1BFB(N3,N2,N1)+trcsa_RFLS(idsa_CaHPO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)
      TC2BFB(N3,N2,N1)=TC2BFB(N3,N2,N1)+trcsa_RFLS(idsa_CaH2PO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)
      TM1BFB(N3,N2,N1)=TM1BFB(N3,N2,N1)+trcsa_RFLS(idsa_MgHPO4B,N,N3,N2,N1)-trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)
      TALFHS(N3,N2,N1)=TALFHS(N3,N2,N1)+RALFHS(N,N3,N2,N1)-RALFHS(N,N6,N5,N4)
      TFEFHS(N3,N2,N1)=TFEFHS(N3,N2,N1)+RFEFHS(N,N3,N2,N1)-RFEFHS(N,N6,N5,N4)
      THYFHS(N3,N2,N1)=THYFHS(N3,N2,N1)+RHYFHS(N,N3,N2,N1)-RHYFHS(N,N6,N5,N4)
      TCAFHS(N3,N2,N1)=TCAFHS(N3,N2,N1)+RCAFHS(N,N3,N2,N1)-RCAFHS(N,N6,N5,N4)
      TMGFHS(N3,N2,N1)=TMGFHS(N3,N2,N1)+RMGFHS(N,N3,N2,N1)-RMGFHS(N,N6,N5,N4)
      TNAFHS(N3,N2,N1)=TNAFHS(N3,N2,N1)+RNAFHS(N,N3,N2,N1)-RNAFHS(N,N6,N5,N4)
      TKAFHS(N3,N2,N1)=TKAFHS(N3,N2,N1)+RKAFHS(N,N3,N2,N1)-RKAFHS(N,N6,N5,N4)
      TOHFHS(N3,N2,N1)=TOHFHS(N3,N2,N1)+ROHFHS(N,N3,N2,N1)-ROHFHS(N,N6,N5,N4)
      TSOFHS(N3,N2,N1)=TSOFHS(N3,N2,N1)+RSOFHS(N,N3,N2,N1)-RSOFHS(N,N6,N5,N4)
      TCLFHS(N3,N2,N1)=TCLFHS(N3,N2,N1)+RCLFHS(N,N3,N2,N1)-RCLFHS(N,N6,N5,N4)
      TC3FHS(N3,N2,N1)=TC3FHS(N3,N2,N1)+RC3FHS(N,N3,N2,N1)-RC3FHS(N,N6,N5,N4)
      THCFHS(N3,N2,N1)=THCFHS(N3,N2,N1)+RHCFHS(N,N3,N2,N1)-RHCFHS(N,N6,N5,N4)
      TAL1HS(N3,N2,N1)=TAL1HS(N3,N2,N1)+RAL1HS(N,N3,N2,N1)-RAL1HS(N,N6,N5,N4)
      TAL2HS(N3,N2,N1)=TAL2HS(N3,N2,N1)+RAL2HS(N,N3,N2,N1)-RAL2HS(N,N6,N5,N4)
      TAL3HS(N3,N2,N1)=TAL3HS(N3,N2,N1)+RAL3HS(N,N3,N2,N1)-RAL3HS(N,N6,N5,N4)
      TAL4HS(N3,N2,N1)=TAL4HS(N3,N2,N1)+RAL4HS(N,N3,N2,N1)-RAL4HS(N,N6,N5,N4)
      TALSHS(N3,N2,N1)=TALSHS(N3,N2,N1)+RALSHS(N,N3,N2,N1)-RALSHS(N,N6,N5,N4)
      TFE1HS(N3,N2,N1)=TFE1HS(N3,N2,N1)+RFE1HS(N,N3,N2,N1)-RFE1HS(N,N6,N5,N4)
      TFE2HS(N3,N2,N1)=TFE2HS(N3,N2,N1)+RFE2HS(N,N3,N2,N1)-RFE2HS(N,N6,N5,N4)
      TFE3HS(N3,N2,N1)=TFE3HS(N3,N2,N1)+RFE3HS(N,N3,N2,N1)-RFE3HS(N,N6,N5,N4)
      TFE4HS(N3,N2,N1)=TFE4HS(N3,N2,N1)+RFE4HS(N,N3,N2,N1)-RFE4HS(N,N6,N5,N4)
      TFESHS(N3,N2,N1)=TFESHS(N3,N2,N1)+RFESHS(N,N3,N2,N1)-RFESHS(N,N6,N5,N4)
      TCAOHS(N3,N2,N1)=TCAOHS(N3,N2,N1)+RCAOHS(N,N3,N2,N1)-RCAOHS(N,N6,N5,N4)
      TCACHS(N3,N2,N1)=TCACHS(N3,N2,N1)+RCACHS(N,N3,N2,N1)-RCACHS(N,N6,N5,N4)
      TCAHHS(N3,N2,N1)=TCAHHS(N3,N2,N1)+RCAHHS(N,N3,N2,N1)-RCAHHS(N,N6,N5,N4)
      TCASHS(N3,N2,N1)=TCASHS(N3,N2,N1)+RCASHS(N,N3,N2,N1)-RCASHS(N,N6,N5,N4)
      TMGOHS(N3,N2,N1)=TMGOHS(N3,N2,N1)+RMGOHS(N,N3,N2,N1)-RMGOHS(N,N6,N5,N4)
      TMGCHS(N3,N2,N1)=TMGCHS(N3,N2,N1)+RMGCHS(N,N3,N2,N1)-RMGCHS(N,N6,N5,N4)
      TMGHHS(N3,N2,N1)=TMGHHS(N3,N2,N1)+RMGHHS(N,N3,N2,N1)-RMGHHS(N,N6,N5,N4)
      TMGSHS(N3,N2,N1)=TMGSHS(N3,N2,N1)+RMGSHS(N,N3,N2,N1)-RMGSHS(N,N6,N5,N4)
      TNACHS(N3,N2,N1)=TNACHS(N3,N2,N1)+RNACHS(N,N3,N2,N1)-RNACHS(N,N6,N5,N4)
      TNASHS(N3,N2,N1)=TNASHS(N3,N2,N1)+RNASHS(N,N3,N2,N1)-RNASHS(N,N6,N5,N4)
      TKASHS(N3,N2,N1)=TKASHS(N3,N2,N1)+RKASHS(N,N3,N2,N1)-RKASHS(N,N6,N5,N4)
      TH0PHS(N3,N2,N1)=TH0PHS(N3,N2,N1)+RH0PHS(N,N3,N2,N1)-RH0PHS(N,N6,N5,N4)
      TH3PHS(N3,N2,N1)=TH3PHS(N3,N2,N1)+RH3PHS(N,N3,N2,N1)-RH3PHS(N,N6,N5,N4)
      TF1PHS(N3,N2,N1)=TF1PHS(N3,N2,N1)+RF1PHS(N,N3,N2,N1)-RF1PHS(N,N6,N5,N4)
      TF2PHS(N3,N2,N1)=TF2PHS(N3,N2,N1)+RF2PHS(N,N3,N2,N1)-RF2PHS(N,N6,N5,N4)
      TC0PHS(N3,N2,N1)=TC0PHS(N3,N2,N1)+RC0PHS(N,N3,N2,N1)-RC0PHS(N,N6,N5,N4)
      TC1PHS(N3,N2,N1)=TC1PHS(N3,N2,N1)+RC1PHS(N,N3,N2,N1)-RC1PHS(N,N6,N5,N4)
      TC2PHS(N3,N2,N1)=TC2PHS(N3,N2,N1)+RC2PHS(N,N3,N2,N1)-RC2PHS(N,N6,N5,N4)
      TM1PHS(N3,N2,N1)=TM1PHS(N3,N2,N1)+RM1PHS(N,N3,N2,N1)-RM1PHS(N,N6,N5,N4)
      TH0BHB(N3,N2,N1)=TH0BHB(N3,N2,N1)+RH0BHB(N,N3,N2,N1)-RH0BHB(N,N6,N5,N4)
      TH3BHB(N3,N2,N1)=TH3BHB(N3,N2,N1)+RH3BHB(N,N3,N2,N1)-RH3BHB(N,N6,N5,N4)
      TF1BHB(N3,N2,N1)=TF1BHB(N3,N2,N1)+RF1BHB(N,N3,N2,N1)-RF1BHB(N,N6,N5,N4)
      TF2BHB(N3,N2,N1)=TF2BHB(N3,N2,N1)+RF2BHB(N,N3,N2,N1)-RF2BHB(N,N6,N5,N4)
      TC0BHB(N3,N2,N1)=TC0BHB(N3,N2,N1)+RC0BHB(N,N3,N2,N1)-RC0BHB(N,N6,N5,N4)
      TC1BHB(N3,N2,N1)=TC1BHB(N3,N2,N1)+RC1BHB(N,N3,N2,N1)-RC1BHB(N,N6,N5,N4)
      TC2BHB(N3,N2,N1)=TC2BHB(N3,N2,N1)+RC2BHB(N,N3,N2,N1)-RC2BHB(N,N6,N5,N4)
      TM1BHB(N3,N2,N1)=TM1BHB(N3,N2,N1)+RM1BHB(N,N3,N2,N1)-RM1BHB(N,N6,N5,N4)
    ELSE
      TALFLS(N3,N2,N1)=0.0
      TFEFLS(N3,N2,N1)=0.0
      THYFLS(N3,N2,N1)=0.0
      TCAFLS(N3,N2,N1)=0.0
      TMGFLS(N3,N2,N1)=0.0
      TNAFLS(N3,N2,N1)=0.0
      TKAFLS(N3,N2,N1)=0.0
      TOHFLS(N3,N2,N1)=0.0
      TSOFLS(N3,N2,N1)=0.0
      TCLFLS(N3,N2,N1)=0.0
      TC3FLS(N3,N2,N1)=0.0
      THCFLS(N3,N2,N1)=0.0
      TAL1FS(N3,N2,N1)=0.0
      TAL2FS(N3,N2,N1)=0.0
      TAL3FS(N3,N2,N1)=0.0
      TAL4FS(N3,N2,N1)=0.0
      TALSFS(N3,N2,N1)=0.0
      TFE1FS(N3,N2,N1)=0.0
      TFE2FS(N3,N2,N1)=0.0
      TFE3FS(N3,N2,N1)=0.0
      TFE4FS(N3,N2,N1)=0.0
      TFESFS(N3,N2,N1)=0.0
      TCAOFS(N3,N2,N1)=0.0
      TCACFS(N3,N2,N1)=0.0
      TCAHFS(N3,N2,N1)=0.0
      TCASFS(N3,N2,N1)=0.0
      TMGOFS(N3,N2,N1)=0.0
      TMGCFS(N3,N2,N1)=0.0
      TMGHFS(N3,N2,N1)=0.0
      TMGSFS(N3,N2,N1)=0.0
      TNACFS(N3,N2,N1)=0.0
      TNASFS(N3,N2,N1)=0.0
      TKASFS(N3,N2,N1)=0.0
      TH0PFS(N3,N2,N1)=0.0
      TH3PFS(N3,N2,N1)=0.0
      TF1PFS(N3,N2,N1)=0.0
      TF2PFS(N3,N2,N1)=0.0
      TC0PFS(N3,N2,N1)=0.0
      TC1PFS(N3,N2,N1)=0.0
      TC2PFS(N3,N2,N1)=0.0
      TM1PFS(N3,N2,N1)=0.0
      TH0BFB(N3,N2,N1)=0.0
      TH3BFB(N3,N2,N1)=0.0
      TF1BFB(N3,N2,N1)=0.0
      TF2BFB(N3,N2,N1)=0.0
      TC0BFB(N3,N2,N1)=0.0
      TC1BFB(N3,N2,N1)=0.0
      TC2BFB(N3,N2,N1)=0.0
      TM1BFB(N3,N2,N1)=0.0
      TALFHS(N3,N2,N1)=0.0
      TFEFHS(N3,N2,N1)=0.0
      THYFHS(N3,N2,N1)=0.0
      TCAFHS(N3,N2,N1)=0.0
      TMGFHS(N3,N2,N1)=0.0
      TNAFHS(N3,N2,N1)=0.0
      TKAFHS(N3,N2,N1)=0.0
      TOHFHS(N3,N2,N1)=0.0
      TSOFHS(N3,N2,N1)=0.0
      TCLFHS(N3,N2,N1)=0.0
      TC3FHS(N3,N2,N1)=0.0
      THCFHS(N3,N2,N1)=0.0
      TAL1HS(N3,N2,N1)=0.0
      TAL2HS(N3,N2,N1)=0.0
      TAL3HS(N3,N2,N1)=0.0
      TAL4HS(N3,N2,N1)=0.0
      TALSHS(N3,N2,N1)=0.0
      TFE1HS(N3,N2,N1)=0.0
      TFE2HS(N3,N2,N1)=0.0
      TFE3HS(N3,N2,N1)=0.0
      TFE4HS(N3,N2,N1)=0.0
      TFESHS(N3,N2,N1)=0.0
      TCAOHS(N3,N2,N1)=0.0
      TCACHS(N3,N2,N1)=0.0
      TCAHHS(N3,N2,N1)=0.0
      TCASHS(N3,N2,N1)=0.0
      TMGOHS(N3,N2,N1)=0.0
      TMGCHS(N3,N2,N1)=0.0
      TMGHHS(N3,N2,N1)=0.0
      TMGSHS(N3,N2,N1)=0.0
      TNACHS(N3,N2,N1)=0.0
      TNASHS(N3,N2,N1)=0.0
      TKASHS(N3,N2,N1)=0.0
      TH0PHS(N3,N2,N1)=0.0
      TH3PHS(N3,N2,N1)=0.0
      TF1PHS(N3,N2,N1)=0.0
      TF2PHS(N3,N2,N1)=0.0
      TC0PHS(N3,N2,N1)=0.0
      TC1PHS(N3,N2,N1)=0.0
      TC2PHS(N3,N2,N1)=0.0
      TM1PHS(N3,N2,N1)=0.0
      TH0BHB(N3,N2,N1)=0.0
      TH3BHB(N3,N2,N1)=0.0
      TF1BHB(N3,N2,N1)=0.0
      TF2BHB(N3,N2,N1)=0.0
      TC0BHB(N3,N2,N1)=0.0
      TC1BHB(N3,N2,N1)=0.0
      TC2BHB(N3,N2,N1)=0.0
      TM1BHB(N3,N2,N1)=0.0
    ENDIF
  ENDIF
  end subroutine TotFluxInMacMicPores
!------------------------------------------------------------------------------------------

      subroutine SaltModelInternalFlux(M,NHW,NHE,NVN,NVS)
!
!     Description:
!
      implicit none
      integer, intent(in) :: M,NHW,NHE,NVN,NVS
      integer :: NY,NX,L,N1,N2,N3,N4,N5,N6,N4B,N5B,N,NN
      integer :: M1,M2,M3,M4,M5,M6
!     begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     M6,M5,M4=L,NY,NX of destination grid cell
!

  D9595: DO NX=NHW,NHE
    D9590: DO NY=NVN,NVS
      D9585: DO L=NU(NY,NX),NL(NY,NX)
        N1=NX
        N2=NY
        N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        D9580: DO N=1,3
          D9575: DO NN=1,2
            IF(N.EQ.1)THEN
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.1)THEN
                IF(NX.EQ.NHE)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX+1
                  M5=NY
                  M6=L
                  XN=-1.0
                  RCHQF=RCHQE(M2,M1)
                  RCHGFU=RCHGEU(M2,M1)
                  RCHGFT=RCHGET(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                IF(NX.EQ.NHW)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RCHQW(M5,M4)
                  RCHGFU=RCHGWU(M5,M4)
                  RCHGFT=RCHGWT(M5,M4)
                ELSE
                  CYCLE
                ENDIF
              ENDIF
            ELSEIF(N.EQ.2)THEN
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.1)THEN
                IF(NY.EQ.NVS)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY+1
                  M6=L
                  XN=-1.0
                  RCHQF=RCHQS(M2,M1)
                  RCHGFU=RCHGSU(M2,M1)
                  RCHGFT=RCHGST(M2,M1)
                ELSE
                  CYCLE
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RCHQN(M5,M4)
                  RCHGFU=RCHGNU(M5,M4)
                  RCHGFT=RCHGNT(M5,M4)
                ELSE
                  CYCLE
                ENDIF
              ENDIF
            ELSEIF(N.EQ.3)THEN
              N1=NX
              N2=NY
              N3=L
              N4=NX
              N5=NY
              N6=L+1
              IF(NN.EQ.1)THEN
                IF(L.EQ.NL(NY,NX))THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L+1
                  XN=-1.0
                ELSE
                  CYCLE
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                CYCLE
              ENDIF
            ENDIF
!
!     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
!     RUNOFF IN WATSUB AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
!
!
!     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
!     RUNOFF IN 'WATSUB' AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
!
!     QRM =runoff from watsub.f
!     RQR*=solute in runoff
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
            IF(L.EQ.NUM(M2,M1).AND.N.NE.3)THEN
              IF(IRCHG(NN,N,N2,N1).EQ.0.OR.test_aeqb(RCHQF,0.0_r8) &
                .OR.QRM(M,N2,N1).LE.ZEROS(N2,N1))THEN
                call ZeroSoluteFluxFromRecharge(N,NN,M5,M4)
              ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
                IF((NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
                  .OR.(NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
                  call SoluteExportThruBoundary(N1,N2,M,N,NN,M5,M4)
!
!     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
                ELSEIF((NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
                  .OR.(NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
                  call ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
                ELSE
                  call ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
                ENDIF
              ENDIF
!
!     BOUNDARY SNOW FLUX
!
              IF(NN.EQ.1)THEN
                call BoundarySnowFlux(N,M5,M4)
              ENDIF
            ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MICROPORE WATER LOSS
!
!     FLWM=water flux through soil micropore from watsub.f
!     VOLWM=micropore water-filled porosity from watsub.f
!     R*FLS=convective solute flux through micropores
!     R*FLW,R*FLB=convective solute flux through micropores in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
            IF(VOLX(N3,N2,N1).GT.ZEROS(NY,NX))THEN

              IF(NCN(M2,M1).NE.3.OR.N.EQ.3)THEN
                IF(NN.EQ.1.AND.FLWM(M,N,M6,M5,M4).GT.0.0 &
                  .OR.NN.EQ.2.AND.FLWM(M,N,M6,M5,M4).LT.0.0)THEN

                  call SoluteLossSubsurfMicropore(M,N,M1,M2,M3,M4,M5,M6)

                ELSE
!     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
!
                  call SoluteGainSubsurfMicropore(M3,M2,M1,M,N,M6,M5,M4)
                ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
!
!     FLWHM=water flux through soil macropore from watsub.f
!     VOLWHM=macropore water-filled porosity from watsub.f
!     RFH*S=solute diffusive flux through macropore
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
                IF(NN.EQ.1.AND.FLWHM(M,N,M6,M5,M4).GT.0.0 &
                  .OR.NN.EQ.2.AND.FLWHM(M,N,M6,M5,M4).LT.0.0)THEN

                  call SoluteLossSubsurfMacropore(M,N,M1,M2,M3,M4,M5,M6)

                ELSE
!
!     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
!
                  call ZeroSolueGainSubsurfMacropore(N,M6,M5,M4)
                ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
                call AccumFluxMacMicPores(N,M6,M5,M4)
              ENDIF
            ENDIF
          ENDDO D9575
!
!     TOTAL SOLUTE FLUXES IN EACH GRID CELL
!
          IF(L.EQ.NUM(N2,N1))THEN
            IF(N.NE.3)THEN
!
!     NET OVERLAND SOLUTE FLUX IN WATER
!
              call NetOverloadFluxInWater(M,N,N1,N2,N4,N5,N4B,N5B)
!
!     NET OVERLAND SOLUTE FLUX IN SNOW
              call NetOverloadFLuxInSnow(N,N1,N2,N4,N5)
!
            ELSEIF(N.EQ.3)THEN
!     NET SOLUTE FLUX IN SNOWPACK
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     T*BLS=net solute flux in snowpack
!     R*BLS=solute flux in snowpack
              call NetFluxInSnowpack(M,NY,NX,N1,N2)
            ENDIF
          ENDIF
!
!     TOTAL SOLUTE FLUX IN MICROPORES AND MACROPORES
!
          call TotFluxInMacMicPores(N,N1,N2,N3,N4,N5,NY,NX,N6)
        ENDDO D9580
      ENDDO D9585
    ENDDO D9590
  ENDDO D9595
  end subroutine SaltModelInternalFlux
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSnow(NY,NX)
!
!     Description:
!
  integer, intent(in) :: NY,NX
  integer :: L,NTA
!     begin_execution
!     *W2=solute content of snowpack
!     TQS*=net overland solute flux in snow
!     T*BLS=net solute flux in snowpack
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO NTA=idsa_beg,idsa_end
    trcs_solsml2(NTA,1,NY,NX)=trcs_solsml2(NTA,1,NY,NX)+trcsa_TQ(NTA,NY,NX)
  ENDDO
  D9670: DO L=1,JS
    DO NTA=idsa_beg,idsa_end
      trcs_solsml2(NTA,L,NY,NX)=trcs_solsml2(NTA,L,NY,NX)+trcsa_TBLS(NTA,L,NY,NX)
    ENDDO
  ENDDO D9670
  end subroutine UpdateSoluteInSnow
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInResidue(NY,NX)
!
!     Description:
!
  integer, intent(in) :: NY,NX

  integer :: NTSA
!     begin_execution
!
!     *S2=litter solute content
!     R*DFR=gas exchange between atmosphere and surface litter water
!     R*DFS=gas exchange between atmosphere and soil surface water
!     R*FLS=convective + diffusive solute flux into litter,soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into litter from non-band,band
!     TQR*=net overland solute flux
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!

  trcsa_solml2(idsa_Al,0,NY,NX)=trcsa_solml2(idsa_Al,0,NY,NX)+trcsa_TQR(idsa_Al,NY,NX)+trcsa_RFLS(idsa_Al,3,0,NY,NX)
  trcsa_solml2(idsa_Fe,0,NY,NX)=trcsa_solml2(idsa_Fe,0,NY,NX)+trcsa_TQR(idsa_Fe,NY,NX)+trcsa_RFLS(idsa_Fe,3,0,NY,NX)
  trcsa_solml2(idsa_Hp,0,NY,NX)=trcsa_solml2(idsa_Hp,0,NY,NX)+trcsa_TQR(idsa_Hp,NY,NX)+trcsa_RFLS(idsa_Hp,3,0,NY,NX)
  trcsa_solml2(idsa_Ca,0,NY,NX)=trcsa_solml2(idsa_Ca,0,NY,NX)+trcsa_TQR(idsa_Ca,NY,NX)+trcsa_RFLS(idsa_Ca,3,0,NY,NX)
  trcsa_solml2(idsa_Mg,0,NY,NX)=trcsa_solml2(idsa_Mg,0,NY,NX)+trcsa_TQR(idsa_Mg,NY,NX)+trcsa_RFLS(idsa_Mg,3,0,NY,NX)
  trcsa_solml2(idsa_Na,0,NY,NX)=trcsa_solml2(idsa_Na,0,NY,NX)+trcsa_TQR(idsa_Na,NY,NX)+trcsa_RFLS(idsa_Na,3,0,NY,NX)
  trcsa_solml2(idsa_K,0,NY,NX)=trcsa_solml2(idsa_K,0,NY,NX)+trcsa_TQR(idsa_K,NY,NX)+trcsa_RFLS(idsa_K,3,0,NY,NX)
  trcsa_solml2(idsa_OH,0,NY,NX)=trcsa_solml2(idsa_OH,0,NY,NX)+trcsa_TQR(idsa_OH,NY,NX)+trcsa_RFLS(idsa_OH,3,0,NY,NX)
  trcsa_solml2(idsa_SO4,0,NY,NX)=trcsa_solml2(idsa_SO4,0,NY,NX)+trcsa_TQR(idsa_SO4,NY,NX)+trcsa_RFLS(idsa_SO4,3,0,NY,NX)
  trcsa_solml2(idsa_Cl,0,NY,NX)=trcsa_solml2(idsa_Cl,0,NY,NX)+trcsa_TQR(idsa_Cl,NY,NX)+trcsa_RFLS(idsa_Cl,3,0,NY,NX)
  trcsa_solml2(idsa_CO3,0,NY,NX)=trcsa_solml2(idsa_CO3,0,NY,NX)+trcsa_TQR(idsa_CO3,NY,NX)+trcsa_RFLS(idsa_CO3,3,0,NY,NX)
  trcsa_solml2(idsa_HCO3,0,NY,NX)=trcsa_solml2(idsa_HCO3,0,NY,NX)+trcsa_TQR(idsa_HCO3,NY,NX)+trcsa_RFLS(idsa_HCO3,3,0,NY,NX)
  trcsa_solml2(idsa_AlOH,0,NY,NX)=trcsa_solml2(idsa_AlOH,0,NY,NX)+trcsa_TQR(idsa_AlOH,NY,NX)+trcsa_RFLS(idsa_AlOH,3,0,NY,NX)
  trcsa_solml2(idsa_AlOH2,0,NY,NX)=trcsa_solml2(idsa_AlOH2,0,NY,NX)+trcsa_TQR(idsa_AlOH2,NY,NX)+trcsa_RFLS(idsa_AlOH2,3,0,NY,NX)
  trcsa_solml2(idsa_AlOH3,0,NY,NX)=trcsa_solml2(idsa_AlOH3,0,NY,NX)+trcsa_TQR(idsa_AlOH3,NY,NX)+trcsa_RFLS(idsa_AlOH3,3,0,NY,NX)
  trcsa_solml2(idsa_AlOH4,0,NY,NX)=trcsa_solml2(idsa_AlOH4,0,NY,NX)+trcsa_TQR(idsa_AlOH4,NY,NX)+trcsa_RFLS(idsa_AlOH4,3,0,NY,NX)
  trcsa_solml2(idsa_AlSO4,0,NY,NX)=trcsa_solml2(idsa_AlSO4,0,NY,NX)+trcsa_TQR(idsa_AlSO4,NY,NX)+trcsa_RFLS(idsa_AlSO4,3,0,NY,NX)
  trcsa_solml2(idsa_FeOH,0,NY,NX)=trcsa_solml2(idsa_FeOH,0,NY,NX)+trcsa_TQR(idsa_FeOH,NY,NX)+trcsa_RFLS(idsa_FeOH,3,0,NY,NX)
  trcsa_solml2(idsa_FeOH2,0,NY,NX)=trcsa_solml2(idsa_FeOH2,0,NY,NX)+trcsa_TQR(idsa_FeOH2,NY,NX)+trcsa_RFLS(idsa_FeOH2,3,0,NY,NX)
  trcsa_solml2(idsa_FeOH3,0,NY,NX)=trcsa_solml2(idsa_FeOH3,0,NY,NX)+trcsa_TQR(idsa_FeOH3,NY,NX)+trcsa_RFLS(idsa_FeOH3,3,0,NY,NX)
  trcsa_solml2(idsa_FeOH4,0,NY,NX)=trcsa_solml2(idsa_FeOH4,0,NY,NX)+trcsa_TQR(idsa_FeOH4,NY,NX)+trcsa_RFLS(idsa_FeOH4,3,0,NY,NX)
  trcsa_solml2(idsa_FeSO4,0,NY,NX)=trcsa_solml2(idsa_FeSO4,0,NY,NX)+trcsa_TQR(idsa_FeSO4,NY,NX)+trcsa_RFLS(idsa_FeSO4,3,0,NY,NX)
  trcsa_solml2(idsa_CaOH2,0,NY,NX)=trcsa_solml2(idsa_CaOH2,0,NY,NX)+trcsa_TQR(idsa_CaOH2,NY,NX)+trcsa_RFLS(idsa_CaOH2,3,0,NY,NX)
  trcsa_solml2(idsa_CaCO3,0,NY,NX)=trcsa_solml2(idsa_CaCO3,0,NY,NX)+trcsa_TQR(idsa_CaCO3,NY,NX)+trcsa_RFLS(idsa_CaCO3,3,0,NY,NX)
  trcsa_solml2(idsa_CaHCO3,0,NY,NX)=trcsa_solml2(idsa_CaHCO3,0,NY,NX)+trcsa_TQR(idsa_CaHCO3,NY,NX)+trcsa_RFLS(idsa_CaHCO3,3,0,NY,NX)
  trcsa_solml2(idsa_CaSO4,0,NY,NX)=trcsa_solml2(idsa_CaSO4,0,NY,NX)+trcsa_TQR(idsa_CaSO4,NY,NX)+trcsa_RFLS(idsa_CaSO4,3,0,NY,NX)
  trcsa_solml2(idsa_MgOH2,0,NY,NX)=trcsa_solml2(idsa_MgOH2,0,NY,NX)+trcsa_TQR(idsa_MgOH2,NY,NX)+trcsa_RFLS(idsa_MgOH2,3,0,NY,NX)
  trcsa_solml2(idsa_MgCO3,0,NY,NX)=trcsa_solml2(idsa_MgCO3,0,NY,NX)+trcsa_TQR(idsa_MgCO3,NY,NX)+trcsa_RFLS(idsa_MgCO3,3,0,NY,NX)
  trcsa_solml2(idsa_MgHCO3,0,NY,NX)=trcsa_solml2(idsa_MgHCO3,0,NY,NX)+trcsa_TQR(idsa_MgHCO3,NY,NX)+trcsa_RFLS(idsa_MgHCO3,3,0,NY,NX)
  trcsa_solml2(idsa_MgSO4,0,NY,NX)=trcsa_solml2(idsa_MgSO4,0,NY,NX)+trcsa_TQR(idsa_MgSO4,NY,NX)+trcsa_RFLS(idsa_MgSO4,3,0,NY,NX)
  trcsa_solml2(idsa_NaCO3,0,NY,NX)=trcsa_solml2(idsa_NaCO3,0,NY,NX)+trcsa_TQR(idsa_NaCO3,NY,NX)+trcsa_RFLS(idsa_NaCO3,3,0,NY,NX)
  trcsa_solml2(idsa_NaSO4,0,NY,NX)=trcsa_solml2(idsa_NaSO4,0,NY,NX)+trcsa_TQR(idsa_NaSO4,NY,NX)+trcsa_RFLS(idsa_NaSO4,3,0,NY,NX)
  trcsa_solml2(idsa_KSO4,0,NY,NX)=trcsa_solml2(idsa_KSO4,0,NY,NX)+trcsa_TQR(idsa_KSO4,NY,NX)+trcsa_RFLS(idsa_KSO4,3,0,NY,NX)
  trcsa_solml2(idsa_H0PO4,0,NY,NX)=trcsa_solml2(idsa_H0PO4,0,NY,NX)+trcsa_TQR(idsa_H0PO4,NY,NX)+trcsa_RFLS(idsa_H0PO4,3,0,NY,NX)
  trcsa_solml2(idsa_H3PO4,0,NY,NX)=trcsa_solml2(idsa_H3PO4,0,NY,NX)+trcsa_TQR(idsa_H3PO4,NY,NX)+trcsa_RFLS(idsa_H3PO4,3,0,NY,NX)
  trcsa_solml2(idsa_FeHPO4,0,NY,NX)=trcsa_solml2(idsa_FeHPO4,0,NY,NX)+trcsa_TQR(idsa_FeHPO4,NY,NX)+trcsa_RFLS(idsa_FeHPO4,3,0,NY,NX)
  trcsa_solml2(idsa_FeH2PO4,0,NY,NX)=trcsa_solml2(idsa_FeH2PO4,0,NY,NX)+trcsa_TQR(idsa_FeH2PO4,NY,NX)+trcsa_RFLS(idsa_FeH2PO4,3,0,NY,NX)
  trcsa_solml2(idsa_CaPO4,0,NY,NX)=trcsa_solml2(idsa_CaPO4,0,NY,NX)+trcsa_TQR(idsa_CaPO4,NY,NX)+trcsa_RFLS(idsa_CaPO4,3,0,NY,NX)
  trcsa_solml2(idsa_CaHPO4,0,NY,NX)=trcsa_solml2(idsa_CaHPO4,0,NY,NX)+trcsa_TQR(idsa_CaHPO4,NY,NX)+trcsa_RFLS(idsa_CaHPO4,3,0,NY,NX)
  trcsa_solml2(idsa_CaH2PO4,0,NY,NX)=trcsa_solml2(idsa_CaH2PO4,0,NY,NX)+trcsa_TQR(idsa_CaH2PO4,NY,NX)+trcsa_RFLS(idsa_CaH2PO4,3,0,NY,NX)
  trcsa_solml2(idsa_MgHPO4,0,NY,NX)=trcsa_solml2(idsa_MgHPO4,0,NY,NX)+trcsa_TQR(idsa_MgHPO4,NY,NX)+trcsa_RFLS(idsa_MgHPO4,3,0,NY,NX)
  end subroutine UpdateSoluteInResidue
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInMicMacpores(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
      integer :: L
!     begin_execution
!
!     *S2,*B2=micropore solute content in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     T*FLS=net convective + diffusive solute flux through micropores
!     T*FHS=net convective + diffusive solute flux through macropores
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!
  D9685: DO L=NU(NY,NX),NL(NY,NX)
    IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      trcsa_solml2(idsa_Al,L,NY,NX)=trcsa_solml2(idsa_Al,L,NY,NX)+TALFLS(L,NY,NX)+RALFXS(L,NY,NX)+trcsa_RFLZ(idsa_Al,L,NY,NX)
      trcsa_solml2(idsa_Fe,L,NY,NX)=trcsa_solml2(idsa_Fe,L,NY,NX)+TFEFLS(L,NY,NX)+RFEFXS(L,NY,NX)+trcsa_RFLZ(idsa_Fe,L,NY,NX)
      trcsa_solml2(idsa_Hp,L,NY,NX)=trcsa_solml2(idsa_Hp,L,NY,NX)+THYFLS(L,NY,NX)+RHYFXS(L,NY,NX)+trcsa_RFLZ(idsa_Hp,L,NY,NX)
      trcsa_solml2(idsa_Ca,L,NY,NX)=trcsa_solml2(idsa_Ca,L,NY,NX)+TCAFLS(L,NY,NX)+RCAFXS(L,NY,NX)+trcsa_RFLZ(idsa_Ca,L,NY,NX)
      trcsa_solml2(idsa_Mg,L,NY,NX)=trcsa_solml2(idsa_Mg,L,NY,NX)+TMGFLS(L,NY,NX)+RMGFXS(L,NY,NX)+trcsa_RFLZ(idsa_Mg,L,NY,NX)
      trcsa_solml2(idsa_Na,L,NY,NX)=trcsa_solml2(idsa_Na,L,NY,NX)+TNAFLS(L,NY,NX)+RNAFXS(L,NY,NX)+trcsa_RFLZ(idsa_Na,L,NY,NX)
      trcsa_solml2(idsa_K,L,NY,NX)=trcsa_solml2(idsa_K,L,NY,NX)+TKAFLS(L,NY,NX)+RKAFXS(L,NY,NX)+trcsa_RFLZ(idsa_K,L,NY,NX)
      trcsa_solml2(idsa_OH,L,NY,NX)=trcsa_solml2(idsa_OH,L,NY,NX)+TOHFLS(L,NY,NX)+ROHFXS(L,NY,NX)+trcsa_RFLZ(idsa_OH,L,NY,NX)
      trcsa_solml2(idsa_SO4,L,NY,NX)=trcsa_solml2(idsa_SO4,L,NY,NX)+TSOFLS(L,NY,NX)+RSOFXS(L,NY,NX)+trcsa_RFLZ(idsa_SO4,L,NY,NX)
      trcsa_solml2(idsa_Cl,L,NY,NX)=trcsa_solml2(idsa_Cl,L,NY,NX)+TCLFLS(L,NY,NX)+RCLFXS(L,NY,NX)+trcsa_RFLZ(idsa_Cl,L,NY,NX)
      trcsa_solml2(idsa_CO3,L,NY,NX)=trcsa_solml2(idsa_CO3,L,NY,NX)+TC3FLS(L,NY,NX)+RC3FXS(L,NY,NX)+trcsa_RFLZ(idsa_CO3,L,NY,NX)
      trcsa_solml2(idsa_HCO3,L,NY,NX)=trcsa_solml2(idsa_HCO3,L,NY,NX)+THCFLS(L,NY,NX)+RHCFXS(L,NY,NX)+trcsa_RFLZ(idsa_HCO3,L,NY,NX)
      trcsa_solml2(idsa_AlOH,L,NY,NX)=trcsa_solml2(idsa_AlOH,L,NY,NX)+TAL1FS(L,NY,NX)+RAL1XS(L,NY,NX)+trcsa_RFLZ(idsa_AlOH,L,NY,NX)
      trcsa_solml2(idsa_AlOH2,L,NY,NX)=trcsa_solml2(idsa_AlOH2,L,NY,NX)+TAL2FS(L,NY,NX)+RAL2XS(L,NY,NX)+trcsa_RFLZ(idsa_AlOH2,L,NY,NX)
      trcsa_solml2(idsa_AlOH3,L,NY,NX)=trcsa_solml2(idsa_AlOH3,L,NY,NX)+TAL3FS(L,NY,NX)+RAL3XS(L,NY,NX)+trcsa_RFLZ(idsa_AlOH3,L,NY,NX)
      trcsa_solml2(idsa_AlOH4,L,NY,NX)=trcsa_solml2(idsa_AlOH4,L,NY,NX)+TAL4FS(L,NY,NX)+RAL4XS(L,NY,NX)+trcsa_RFLZ(idsa_AlOH4,L,NY,NX)
      trcsa_solml2(idsa_AlSO4,L,NY,NX)=trcsa_solml2(idsa_AlSO4,L,NY,NX)+TALSFS(L,NY,NX)+RALSXS(L,NY,NX)+trcsa_RFLZ(idsa_AlSO4,L,NY,NX)
      trcsa_solml2(idsa_FeOH,L,NY,NX)=trcsa_solml2(idsa_FeOH,L,NY,NX)+TFE1FS(L,NY,NX)+RFE1XS(L,NY,NX)+trcsa_RFLZ(idsa_FeOH,L,NY,NX)
      trcsa_solml2(idsa_FeOH2,L,NY,NX)=trcsa_solml2(idsa_FeOH2,L,NY,NX)+TFE2FS(L,NY,NX)+RFE2XS(L,NY,NX)+trcsa_RFLZ(idsa_FeOH2,L,NY,NX)
      trcsa_solml2(idsa_FeOH3,L,NY,NX)=trcsa_solml2(idsa_FeOH3,L,NY,NX)+TFE3FS(L,NY,NX)+RFE3XS(L,NY,NX)+trcsa_RFLZ(idsa_FeOH3,L,NY,NX)
      trcsa_solml2(idsa_FeOH4,L,NY,NX)=trcsa_solml2(idsa_FeOH4,L,NY,NX)+TFE4FS(L,NY,NX)+RFE4XS(L,NY,NX)+trcsa_RFLZ(idsa_FeOH4,L,NY,NX)
      trcsa_solml2(idsa_FeSO4,L,NY,NX)=trcsa_solml2(idsa_FeSO4,L,NY,NX)+TFESFS(L,NY,NX)+RFESXS(L,NY,NX)+trcsa_RFLZ(idsa_FeSO4,L,NY,NX)
      trcsa_solml2(idsa_CaOH2,L,NY,NX)=trcsa_solml2(idsa_CaOH2,L,NY,NX)+TCAOFS(L,NY,NX)+RCAOXS(L,NY,NX)+trcsa_RFLZ(idsa_CaOH2,L,NY,NX)
      trcsa_solml2(idsa_CaCO3,L,NY,NX)=trcsa_solml2(idsa_CaCO3,L,NY,NX)+TCACFS(L,NY,NX)+RCACXS(L,NY,NX)+trcsa_RFLZ(idsa_CaCO3,L,NY,NX)
      trcsa_solml2(idsa_CaHCO3,L,NY,NX)=trcsa_solml2(idsa_CaHCO3,L,NY,NX)+TCAHFS(L,NY,NX)+RCAHXS(L,NY,NX)+trcsa_RFLZ(idsa_CaHCO3,L,NY,NX)
      trcsa_solml2(idsa_CaSO4,L,NY,NX)=trcsa_solml2(idsa_CaSO4,L,NY,NX)+TCASFS(L,NY,NX)+RCASXS(L,NY,NX)+trcsa_RFLZ(idsa_CaSO4,L,NY,NX)
      trcsa_solml2(idsa_MgOH2,L,NY,NX)=trcsa_solml2(idsa_MgOH2,L,NY,NX)+TMGOFS(L,NY,NX)+RMGOXS(L,NY,NX)+trcsa_RFLZ(idsa_MgOH2,L,NY,NX)
      trcsa_solml2(idsa_MgCO3,L,NY,NX)=trcsa_solml2(idsa_MgCO3,L,NY,NX)+TMGCFS(L,NY,NX)+RMGCXS(L,NY,NX)+trcsa_RFLZ(idsa_MgCO3,L,NY,NX)
      trcsa_solml2(idsa_MgHCO3,L,NY,NX)=trcsa_solml2(idsa_MgHCO3,L,NY,NX)+TMGHFS(L,NY,NX)+RMGHXS(L,NY,NX)+trcsa_RFLZ(idsa_MgHCO3,L,NY,NX)
      trcsa_solml2(idsa_MgSO4,L,NY,NX)=trcsa_solml2(idsa_MgSO4,L,NY,NX)+TMGSFS(L,NY,NX)+RMGSXS(L,NY,NX)+trcsa_RFLZ(idsa_MgSO4,L,NY,NX)
      trcsa_solml2(idsa_NaCO3,L,NY,NX)=trcsa_solml2(idsa_NaCO3,L,NY,NX)+TNACFS(L,NY,NX)+RNACXS(L,NY,NX)+trcsa_RFLZ(idsa_NaCO3,L,NY,NX)
      trcsa_solml2(idsa_NaSO4,L,NY,NX)=trcsa_solml2(idsa_NaSO4,L,NY,NX)+TNASFS(L,NY,NX)+RNASXS(L,NY,NX)+trcsa_RFLZ(idsa_NaSO4,L,NY,NX)
      trcsa_solml2(idsa_KSO4,L,NY,NX)=trcsa_solml2(idsa_KSO4,L,NY,NX)+TKASFS(L,NY,NX)+RKASXS(L,NY,NX)+trcsa_RFLZ(idsa_KSO4,L,NY,NX)
      trcsa_solml2(idsa_H0PO4,L,NY,NX)=trcsa_solml2(idsa_H0PO4,L,NY,NX)+TH0PFS(L,NY,NX)+RH0PXS(L,NY,NX)+trcsa_RFLZ(idsa_H0PO4,L,NY,NX)
      trcsa_solml2(idsa_H3PO4,L,NY,NX)=trcsa_solml2(idsa_H3PO4,L,NY,NX)+TH3PFS(L,NY,NX)+RH3PXS(L,NY,NX)+trcsa_RFLZ(idsa_H3PO4,L,NY,NX)
      trcsa_solml2(idsa_FeHPO4,L,NY,NX)=trcsa_solml2(idsa_FeHPO4,L,NY,NX)+TF1PFS(L,NY,NX)+RF1PXS(L,NY,NX)+trcsa_RFLZ(idsa_FeHPO4,L,NY,NX)
      trcsa_solml2(idsa_FeH2PO4,L,NY,NX)=trcsa_solml2(idsa_FeH2PO4,L,NY,NX)+TF2PFS(L,NY,NX)+RF2PXS(L,NY,NX)+trcsa_RFLZ(idsa_FeH2PO4,L,NY,NX)
      trcsa_solml2(idsa_CaPO4,L,NY,NX)=trcsa_solml2(idsa_CaPO4,L,NY,NX)+TC0PFS(L,NY,NX)+RC0PXS(L,NY,NX)+trcsa_RFLZ(idsa_CaPO4,L,NY,NX)
      trcsa_solml2(idsa_CaHPO4,L,NY,NX)=trcsa_solml2(idsa_CaHPO4,L,NY,NX)+TC1PFS(L,NY,NX)+RC1PXS(L,NY,NX)+trcsa_RFLZ(idsa_CaHPO4,L,NY,NX)
      trcsa_solml2(idsa_CaH2PO4,L,NY,NX)=trcsa_solml2(idsa_CaH2PO4,L,NY,NX)+TC2PFS(L,NY,NX)+RC2PXS(L,NY,NX)+trcsa_RFLZ(idsa_CaH2PO4,L,NY,NX)
      trcsa_solml2(idsa_MgHPO4,L,NY,NX)=trcsa_solml2(idsa_MgHPO4,L,NY,NX)+TM1PFS(L,NY,NX)+RM1PXS(L,NY,NX)+trcsa_RFLZ(idsa_MgHPO4,L,NY,NX)
      trcsa_solml2(idsa_H0PO4B,L,NY,NX)=trcsa_solml2(idsa_H0PO4B,L,NY,NX)+TH0BFB(L,NY,NX)+RH0BXB(L,NY,NX)+trcsa_RFLZ(idsa_H0PO4B,L,NY,NX)
      trcsa_solml2(idsa_H3PO4B,L,NY,NX)=trcsa_solml2(idsa_H3PO4B,L,NY,NX)+TH3BFB(L,NY,NX)+RH3BXB(L,NY,NX)+trcsa_RFLZ(idsa_H3PO4B,L,NY,NX)
      trcsa_solml2(idsa_FeHPO4B,L,NY,NX)=trcsa_solml2(idsa_FeHPO4B,L,NY,NX)+TF1BFB(L,NY,NX)+RF1BXB(L,NY,NX)+trcsa_RFLZ(idsa_FeHPO4B,L,NY,NX)
      trcsa_solml2(idsa_FeH2PO4B,L,NY,NX)=trcsa_solml2(idsa_FeH2PO4B,L,NY,NX)+TF2BFB(L,NY,NX)+RF2BXB(L,NY,NX)+trcsa_RFLZ(idsa_FeH2PO4B,L,NY,NX)
      trcsa_solml2(idsa_CaPO4B,L,NY,NX)=trcsa_solml2(idsa_CaPO4B,L,NY,NX)+TC0BFB(L,NY,NX)+RC0BXB(L,NY,NX)+trcsa_RFLZ(idsa_CaPO4B,L,NY,NX)
      trcsa_solml2(idsa_CaHPO4B,L,NY,NX)=trcsa_solml2(idsa_CaHPO4B,L,NY,NX)+TC1BFB(L,NY,NX)+RC1BXB(L,NY,NX)+trcsa_RFLZ(idsa_CaHPO4B,L,NY,NX)
      trcsa_solml2(idsa_CaH2PO4B,L,NY,NX)=trcsa_solml2(idsa_CaH2PO4B,L,NY,NX)+TC2BFB(L,NY,NX)+RC2BXB(L,NY,NX)+trcsa_RFLZ(idsa_CaH2PO4B,L,NY,NX)
      trcsa_solml2(idsa_MgHPO4B,L,NY,NX)=trcsa_solml2(idsa_MgHPO4B,L,NY,NX)+TM1BFB(L,NY,NX)+RM1BXB(L,NY,NX)+trcsa_RFLZ(idsa_MgHPO4B,L,NY,NX)
      ZALH2(L,NY,NX)=ZALH2(L,NY,NX)+TALFHS(L,NY,NX)-RALFXS(L,NY,NX)
      ZFEH2(L,NY,NX)=ZFEH2(L,NY,NX)+TFEFHS(L,NY,NX)-RFEFXS(L,NY,NX)
      ZHYH2(L,NY,NX)=ZHYH2(L,NY,NX)+THYFHS(L,NY,NX)-RHYFXS(L,NY,NX)
      ZCCH2(L,NY,NX)=ZCCH2(L,NY,NX)+TCAFHS(L,NY,NX)-RCAFXS(L,NY,NX)
      ZMAH2(L,NY,NX)=ZMAH2(L,NY,NX)+TMGFHS(L,NY,NX)-RMGFXS(L,NY,NX)
      ZNAH2(L,NY,NX)=ZNAH2(L,NY,NX)+TNAFHS(L,NY,NX)-RNAFXS(L,NY,NX)
      ZKAH2(L,NY,NX)=ZKAH2(L,NY,NX)+TKAFHS(L,NY,NX)-RKAFXS(L,NY,NX)
      ZOHH2(L,NY,NX)=ZOHH2(L,NY,NX)+TOHFHS(L,NY,NX)-ROHFXS(L,NY,NX)
      ZSO4H2(L,NY,NX)=ZSO4H2(L,NY,NX)+TSOFHS(L,NY,NX)-RSOFXS(L,NY,NX)
      ZCLH2(L,NY,NX)=ZCLH2(L,NY,NX)+TCLFHS(L,NY,NX)-RCLFXS(L,NY,NX)
      ZCO3H2(L,NY,NX)=ZCO3H2(L,NY,NX)+TC3FHS(L,NY,NX)-RC3FXS(L,NY,NX)
      ZHCOH2(L,NY,NX)=ZHCOH2(L,NY,NX)+THCFHS(L,NY,NX)-RHCFXS(L,NY,NX)
      ZAL1H2(L,NY,NX)=ZAL1H2(L,NY,NX)+TAL1HS(L,NY,NX)-RAL1XS(L,NY,NX)
      ZAL2H2(L,NY,NX)=ZAL2H2(L,NY,NX)+TAL2HS(L,NY,NX)-RAL2XS(L,NY,NX)
      ZAL3H2(L,NY,NX)=ZAL3H2(L,NY,NX)+TAL3HS(L,NY,NX)-RAL3XS(L,NY,NX)
      ZAL4H2(L,NY,NX)=ZAL4H2(L,NY,NX)+TAL4HS(L,NY,NX)-RAL4XS(L,NY,NX)
      ZALSH2(L,NY,NX)=ZALSH2(L,NY,NX)+TALSHS(L,NY,NX)-RALSXS(L,NY,NX)
      ZFE1H2(L,NY,NX)=ZFE1H2(L,NY,NX)+TFE1HS(L,NY,NX)-RFE1XS(L,NY,NX)
      ZFE2H2(L,NY,NX)=ZFE2H2(L,NY,NX)+TFE2HS(L,NY,NX)-RFE2XS(L,NY,NX)
      ZFE3H2(L,NY,NX)=ZFE3H2(L,NY,NX)+TFE3HS(L,NY,NX)-RFE3XS(L,NY,NX)
      ZFE4H2(L,NY,NX)=ZFE4H2(L,NY,NX)+TFE4HS(L,NY,NX)-RFE4XS(L,NY,NX)
      ZFESH2(L,NY,NX)=ZFESH2(L,NY,NX)+TFESHS(L,NY,NX)-RFESXS(L,NY,NX)
      ZCAOH2(L,NY,NX)=ZCAOH2(L,NY,NX)+TCAOHS(L,NY,NX)-RCAOXS(L,NY,NX)
      ZCACH2(L,NY,NX)=ZCACH2(L,NY,NX)+TCACHS(L,NY,NX)-RCACXS(L,NY,NX)
      ZCAHH2(L,NY,NX)=ZCAHH2(L,NY,NX)+TCAHHS(L,NY,NX)-RCAHXS(L,NY,NX)
      ZCASH2(L,NY,NX)=ZCASH2(L,NY,NX)+TCASHS(L,NY,NX)-RCASXS(L,NY,NX)
      ZMGOH2(L,NY,NX)=ZMGOH2(L,NY,NX)+TMGOHS(L,NY,NX)-RMGOXS(L,NY,NX)
      ZMGCH2(L,NY,NX)=ZMGCH2(L,NY,NX)+TMGCHS(L,NY,NX)-RMGCXS(L,NY,NX)
      ZMGHH2(L,NY,NX)=ZMGHH2(L,NY,NX)+TMGHHS(L,NY,NX)-RMGHXS(L,NY,NX)
      ZMGSH2(L,NY,NX)=ZMGSH2(L,NY,NX)+TMGSHS(L,NY,NX)-RMGSXS(L,NY,NX)
      ZNACH2(L,NY,NX)=ZNACH2(L,NY,NX)+TNACHS(L,NY,NX)-RNACXS(L,NY,NX)
      ZNASH2(L,NY,NX)=ZNASH2(L,NY,NX)+TNASHS(L,NY,NX)-RNASXS(L,NY,NX)
      ZKASH2(L,NY,NX)=ZKASH2(L,NY,NX)+TKASHS(L,NY,NX)-RKASXS(L,NY,NX)
      H0P4H2(L,NY,NX)=H0P4H2(L,NY,NX)+TH0PHS(L,NY,NX)-RH0PXS(L,NY,NX)
      H3P4H2(L,NY,NX)=H3P4H2(L,NY,NX)+TH3PHS(L,NY,NX)-RH3PXS(L,NY,NX)
      ZF1PH2(L,NY,NX)=ZF1PH2(L,NY,NX)+TF1PHS(L,NY,NX)-RF1PXS(L,NY,NX)
      ZF2PH2(L,NY,NX)=ZF2PH2(L,NY,NX)+TF2PHS(L,NY,NX)-RF2PXS(L,NY,NX)
      ZC0PH2(L,NY,NX)=ZC0PH2(L,NY,NX)+TC0PHS(L,NY,NX)-RC0PXS(L,NY,NX)
      ZC1PH2(L,NY,NX)=ZC1PH2(L,NY,NX)+TC1PHS(L,NY,NX)-RC1PXS(L,NY,NX)
      ZC2PH2(L,NY,NX)=ZC2PH2(L,NY,NX)+TC2PHS(L,NY,NX)-RC2PXS(L,NY,NX)
      ZM1PH2(L,NY,NX)=ZM1PH2(L,NY,NX)+TM1PHS(L,NY,NX)-RM1PXS(L,NY,NX)
      H0PBH2(L,NY,NX)=H0PBH2(L,NY,NX)+TH0BHB(L,NY,NX)-RH0BXB(L,NY,NX)
      H3PBH2(L,NY,NX)=H3PBH2(L,NY,NX)+TH3BHB(L,NY,NX)-RH3BXB(L,NY,NX)
      ZF1BH2(L,NY,NX)=ZF1BH2(L,NY,NX)+TF1BHB(L,NY,NX)-RF1BXB(L,NY,NX)
      ZF2BH2(L,NY,NX)=ZF2BH2(L,NY,NX)+TF2BHB(L,NY,NX)-RF2BXB(L,NY,NX)
      ZC0BH2(L,NY,NX)=ZC0BH2(L,NY,NX)+TC0BHB(L,NY,NX)-RC0BXB(L,NY,NX)
      ZC1BH2(L,NY,NX)=ZC1BH2(L,NY,NX)+TC1BHB(L,NY,NX)-RC1BXB(L,NY,NX)
      ZC2BH2(L,NY,NX)=ZC2BH2(L,NY,NX)+TC2BHB(L,NY,NX)-RC2BXB(L,NY,NX)
      ZM1BH2(L,NY,NX)=ZM1BH2(L,NY,NX)+TM1BHB(L,NY,NX)-RM1BXB(L,NY,NX)
!     IF(I.EQ.268.AND.L.EQ.1)THEN
!     WRITE(*,444)'ZOH2',I,J,M,NX,NY,L,trcsa_solml2(idsa_OH,L,NY,NX)
!    2,TOHFLS(L,NY,NX),ROHFXS(L,NY,NX),trcsa_RFLZ(idsa_OH,L,NY,NX)
!    3,trcsa_RFLS(idsa_OH,3,L-1,NY,NX),trcsa_RFLS(idsa_OH,3,L,NY,NX)
!     WRITE(*,444)'ZAL2',I,J,M,NX,NY,L,trcsa_solml2(idsa_Al,L,NY,NX)
!    2,TALFLS(L,NY,NX),RALFXS(L,NY,NX),trcsa_RFLZ(idsa_Al,L,NY,NX)
!    3,trcsa_solml2R(idsa_Al,L,NY,NX),trcsa_TR(idsa_Al,L,NY,NX)
!444   FORMAT(A8,6I4,20E12.4)
!     ENDIF
    ENDIF
  ENDDO D9685
  end subroutine UpdateSoluteInMicMacpores

!------------------------------------------------------------------------------------------

  subroutine InitFluxArrays(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     INITIALIZE SOLUTE RUNOFF NET FLUX ACCUMULATORS
      call InitFluxAccumlatorsInRunoff(NY,NX)
!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
      call InitFluxAccumulatorsInSnowpack(NY,NX)
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
      call InitFluxAccumulatorsInSoil(NY,NX)
    ENDDO
  ENDDO
  end subroutine InitFluxArrays


!------------------------------------------------------------------------------------------

  subroutine InitFluxAccumlatorsInRunoff(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
!
  trcsa_TQR(idsa_beg:idsa_end,NY,NX)=0.0

  trcsa_TQ(idsa_beg:idsa_end,NY,NX)=0.0
  end subroutine InitFluxAccumlatorsInRunoff
!------------------------------------------------------------------------------------------

  subroutine InitFluxAccumulatorsInSoil(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  D9885: DO L=NU(NY,NX),NL(NY,NX)
    TALFLS(L,NY,NX)=0.0
    TFEFLS(L,NY,NX)=0.0
    THYFLS(L,NY,NX)=0.0
    TCAFLS(L,NY,NX)=0.0
    TMGFLS(L,NY,NX)=0.0
    TNAFLS(L,NY,NX)=0.0
    TKAFLS(L,NY,NX)=0.0
    TOHFLS(L,NY,NX)=0.0
    TSOFLS(L,NY,NX)=0.0
    TCLFLS(L,NY,NX)=0.0
    TC3FLS(L,NY,NX)=0.0
    THCFLS(L,NY,NX)=0.0
    TAL1FS(L,NY,NX)=0.0
    TAL2FS(L,NY,NX)=0.0
    TAL3FS(L,NY,NX)=0.0
    TAL4FS(L,NY,NX)=0.0
    TALSFS(L,NY,NX)=0.0
    TFE1FS(L,NY,NX)=0.0
    TFE2FS(L,NY,NX)=0.0
    TFE3FS(L,NY,NX)=0.0
    TFE4FS(L,NY,NX)=0.0
    TFESFS(L,NY,NX)=0.0
    TCAOFS(L,NY,NX)=0.0
    TCACFS(L,NY,NX)=0.0
    TCAHFS(L,NY,NX)=0.0
    TCASFS(L,NY,NX)=0.0
    TMGOFS(L,NY,NX)=0.0
    TMGCFS(L,NY,NX)=0.0
    TMGHFS(L,NY,NX)=0.0
    TMGSFS(L,NY,NX)=0.0
    TNACFS(L,NY,NX)=0.0
    TNASFS(L,NY,NX)=0.0
    TKASFS(L,NY,NX)=0.0
    TH0PFS(L,NY,NX)=0.0
    TH3PFS(L,NY,NX)=0.0
    TF1PFS(L,NY,NX)=0.0
    TF2PFS(L,NY,NX)=0.0
    TC0PFS(L,NY,NX)=0.0
    TC1PFS(L,NY,NX)=0.0
    TC2PFS(L,NY,NX)=0.0
    TM1PFS(L,NY,NX)=0.0
    TH0BFB(L,NY,NX)=0.0
    TH3BFB(L,NY,NX)=0.0
    TF1BFB(L,NY,NX)=0.0
    TF2BFB(L,NY,NX)=0.0
    TC0BFB(L,NY,NX)=0.0
    TC1BFB(L,NY,NX)=0.0
    TC2BFB(L,NY,NX)=0.0
    TM1BFB(L,NY,NX)=0.0
    TALFHS(L,NY,NX)=0.0
    TFEFHS(L,NY,NX)=0.0
    THYFHS(L,NY,NX)=0.0
    TCAFHS(L,NY,NX)=0.0
    TMGFHS(L,NY,NX)=0.0
    TNAFHS(L,NY,NX)=0.0
    TKAFHS(L,NY,NX)=0.0
    TOHFHS(L,NY,NX)=0.0
    TSOFHS(L,NY,NX)=0.0
    TCLFHS(L,NY,NX)=0.0
    TC3FHS(L,NY,NX)=0.0
    THCFHS(L,NY,NX)=0.0
    TAL1HS(L,NY,NX)=0.0
    TAL2HS(L,NY,NX)=0.0
    TAL3HS(L,NY,NX)=0.0
    TAL4HS(L,NY,NX)=0.0
    TALSHS(L,NY,NX)=0.0
    TFE1HS(L,NY,NX)=0.0
    TFE2HS(L,NY,NX)=0.0
    TFE3HS(L,NY,NX)=0.0
    TFE4HS(L,NY,NX)=0.0
    TFESHS(L,NY,NX)=0.0
    TCAOHS(L,NY,NX)=0.0
    TCACHS(L,NY,NX)=0.0
    TCAHHS(L,NY,NX)=0.0
    TCASHS(L,NY,NX)=0.0
    TMGOHS(L,NY,NX)=0.0
    TMGCHS(L,NY,NX)=0.0
    TMGHHS(L,NY,NX)=0.0
    TMGSHS(L,NY,NX)=0.0
    TNACHS(L,NY,NX)=0.0
    TNASHS(L,NY,NX)=0.0
    TKASHS(L,NY,NX)=0.0
    TH0PHS(L,NY,NX)=0.0
    TH3PHS(L,NY,NX)=0.0
    TF1PHS(L,NY,NX)=0.0
    TF2PHS(L,NY,NX)=0.0
    TC0PHS(L,NY,NX)=0.0
    TC1PHS(L,NY,NX)=0.0
    TC2PHS(L,NY,NX)=0.0
    TM1PHS(L,NY,NX)=0.0
    TH0BHB(L,NY,NX)=0.0
    TH3BHB(L,NY,NX)=0.0
    TF1BHB(L,NY,NX)=0.0
    TF2BHB(L,NY,NX)=0.0
    TC0BHB(L,NY,NX)=0.0
    TC1BHB(L,NY,NX)=0.0
    TC2BHB(L,NY,NX)=0.0
    TM1BHB(L,NY,NX)=0.0
  ENDDO D9885
  end subroutine InitFluxAccumulatorsInSoil

end module TrnsfrsMod
