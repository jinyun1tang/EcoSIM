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
  real(r8),allocatable ::  trcsa_TFLS(:,:,:,:)                      !
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
  call destroy(trcsa_TFLS)

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

  allocate(trcsa_TFLS(idsa_beg:idsab_end,JZ,JY,JX));   trcsa_TFLS=0._r8
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

  integer :: NTSA
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
  DO NTSA=idsa_beg,idsa_end
    trcsa_RBLS(NTSA,1,NY,NX)=trcsa_XBLS(NTSA,1,NY,NX)*XNPH
  ENDDO

  trcsa_RFL0(idsa_Al,NY,NX)=trcsa_XFLS(idsa_Al,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_Fe,NY,NX)=trcsa_XFLS(idsa_Fe,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_Hp,NY,NX)=trcsa_XFLS(idsa_Hp,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_Ca,NY,NX)=trcsa_XFLS(idsa_Ca,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_Mg,NY,NX)=trcsa_XFLS(idsa_Mg,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_Na,NY,NX)=trcsa_XFLS(idsa_Na,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_K,NY,NX)=trcsa_XFLS(idsa_K,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_OH,NY,NX)=trcsa_XFLS(idsa_OH,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_SO4,NY,NX)=trcsa_XFLS(idsa_SO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_Cl,NY,NX)=trcsa_XFLS(idsa_Cl,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CO3,NY,NX)=trcsa_XFLS(idsa_CO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_HCO3,NY,NX)=trcsa_XFLS(idsa_HCO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_AlOH,NY,NX)=trcsa_XFLS(idsa_AlOH,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_AlOH2,NY,NX)=trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_AlOH3,NY,NX)=trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_AlOH4,NY,NX)=trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_AlSO4,NY,NX)=trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeOH,NY,NX)=trcsa_XFLS(idsa_FeOH,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeOH2,NY,NX)=trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeOH3,NY,NX)=trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeOH4,NY,NX)=trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeSO4,NY,NX)=trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaOH2,NY,NX)=trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaCO3,NY,NX)=trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaHCO3,NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaSO4,NY,NX)=trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_MgOH2,NY,NX)=trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_MgCO3,NY,NX)=trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_MgHCO3,NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_MgSO4,NY,NX)=trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_NaCO3,NY,NX)=trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_NaSO4,NY,NX)=trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_KSO4,NY,NX)=trcsa_XFLS(idsa_KSO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_H0PO4,NY,NX)=trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_H3PO4,NY,NX)=trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeHPO4,NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_FeH2PO4,NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaPO4,NY,NX)=trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaHPO4,NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_CaH2PO4,NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)*XNPH
  trcsa_RFL0(idsa_MgHPO4,NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)*XNPH

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
    DO NTSA=idsa_beg,idsab_end
      trcsa_solml2R(NTSA,L,NY,NX)=-trcsa_TR(NTSA,L,NY,NX)*XNPH
    ENDDO

    trcsa_solml2R(idsa_Hp,L,NY,NX)=trcsa_solml2R(idsa_Hp,L,NY,NX)-(XZHYS(L,NY,NX))*XNPH
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
    DO NTSA=idsa_beg,idsab_end
      trcsa_solml2(NTSA,L,NY,NX)=trcsa_solml(NTSA,L,NY,NX)
    ENDDO

    trcsa_soHml2(idsa_Al,L,NY,NX)=trcsa_soHml(idsa_Al,L,NY,NX)
    trcsa_soHml2(idsa_Fe,L,NY,NX)=trcsa_soHml(idsa_Fe,L,NY,NX)
    trcsa_soHml2(idsa_Hp,L,NY,NX)=trcsa_soHml(idsa_Hp,L,NY,NX)
    trcsa_soHml2(idsa_Ca,L,NY,NX)=trcsa_soHml(idsa_Ca,L,NY,NX)
    trcsa_soHml2(idsa_Mg,L,NY,NX)=trcsa_soHml(idsa_Mg,L,NY,NX)
    trcsa_soHml2(idsa_Na,L,NY,NX)=trcsa_soHml(idsa_Na,L,NY,NX)
    trcsa_soHml2(idsa_K,L,NY,NX)=trcsa_soHml(idsa_K,L,NY,NX)
    trcsa_soHml2(idsa_OH,L,NY,NX)=trcsa_soHml(idsa_OH,L,NY,NX)
    trcsa_soHml2(idsa_SO4,L,NY,NX)=trcsa_soHml(idsa_SO4,L,NY,NX)
    trcsa_soHml2(idsa_Cl,L,NY,NX)=trcsa_soHml(idsa_Cl,L,NY,NX)
    trcsa_soHml2(idsa_CO3,L,NY,NX)=trcsa_soHml(idsa_CO3,L,NY,NX)
    trcsa_soHml2(idsa_HCO3,L,NY,NX)=trcsa_soHml(idsa_HCO3,L,NY,NX)
    trcsa_soHml2(idsa_AlOH,L,NY,NX)=trcsa_soHml(idsa_AlOH,L,NY,NX)
    trcsa_soHml2(idsa_AlOH2,L,NY,NX)=trcsa_soHml(idsa_AlOH2,L,NY,NX)
    trcsa_soHml2(idsa_AlOH3,L,NY,NX)=trcsa_soHml(idsa_AlOH3,L,NY,NX)
    trcsa_soHml2(idsa_AlOH4,L,NY,NX)=trcsa_soHml(idsa_AlOH4,L,NY,NX)
    trcsa_soHml2(idsa_AlSO4,L,NY,NX)=trcsa_soHml(idsa_AlSO4,L,NY,NX)
    trcsa_soHml2(idsa_FeOH,L,NY,NX)=trcsa_soHml(idsa_FeOH,L,NY,NX)
    trcsa_soHml2(idsa_FeOH2,L,NY,NX)=trcsa_soHml(idsa_FeOH2,L,NY,NX)
    trcsa_soHml2(idsa_FeOH3,L,NY,NX)=trcsa_soHml(idsa_FeOH3,L,NY,NX)
    trcsa_soHml2(idsa_FeOH4,L,NY,NX)=trcsa_soHml(idsa_FeOH4,L,NY,NX)
    trcsa_soHml2(idsa_FeSO4,L,NY,NX)=trcsa_soHml(idsa_FeSO4,L,NY,NX)
    trcsa_soHml2(idsa_CaOH2,L,NY,NX)=trcsa_soHml(idsa_CaOH2,L,NY,NX)
    trcsa_soHml2(idsa_CaCO3,L,NY,NX)=trcsa_soHml(idsa_CaCO3,L,NY,NX)
    trcsa_soHml2(idsa_CaHCO3,L,NY,NX)=trcsa_soHml(idsa_CaHCO3,L,NY,NX)
    trcsa_soHml2(idsa_CaSO4,L,NY,NX)=trcsa_soHml(idsa_CaSO4,L,NY,NX)
    trcsa_soHml2(idsa_MgOH2,L,NY,NX)=trcsa_soHml(idsa_MgOH2,L,NY,NX)
    trcsa_soHml2(idsa_MgCO3,L,NY,NX)=trcsa_soHml(idsa_MgCO3,L,NY,NX)
    trcsa_soHml2(idsa_MgHCO3,L,NY,NX)=trcsa_soHml(idsa_MgHCO3,L,NY,NX)
    trcsa_soHml2(idsa_MgSO4,L,NY,NX)=trcsa_soHml(idsa_MgSO4,L,NY,NX)
    trcsa_soHml2(idsa_NaCO3,L,NY,NX)=trcsa_soHml(idsa_NaCO3,L,NY,NX)
    trcsa_soHml2(idsa_NaSO4,L,NY,NX)=trcsa_soHml(idsa_NaSO4,L,NY,NX)
    trcsa_soHml2(idsa_KSO4,L,NY,NX)=trcsa_soHml(idsa_KSO4,L,NY,NX)
    trcsa_soHml2(idsa_H0PO4,L,NY,NX)=trcsa_soHml(idsa_H0PO4,L,NY,NX)
    trcsa_soHml2(idsa_H3PO4,L,NY,NX)=trcsa_soHml(idsa_H3PO4,L,NY,NX)
    trcsa_soHml2(idsa_FeHPO4,L,NY,NX)=trcsa_soHml(idsa_FeHPO4,L,NY,NX)
    trcsa_soHml2(idsa_FeH2PO4,L,NY,NX)=trcsa_soHml(idsa_FeH2PO4,L,NY,NX)
    trcsa_soHml2(idsa_CaPO4,L,NY,NX)=trcsa_soHml(idsa_CaPO4,L,NY,NX)
    trcsa_soHml2(idsa_CaHPO4,L,NY,NX)=trcsa_soHml(idsa_CaHPO4,L,NY,NX)
    trcsa_soHml2(idsa_CaH2PO4,L,NY,NX)=trcsa_soHml(idsa_CaH2PO4,L,NY,NX)
    trcsa_soHml2(idsa_MgHPO4,L,NY,NX)=trcsa_soHml(idsa_MgHPO4,L,NY,NX)
    trcsa_soHml2(idsa_H0PO4B,L,NY,NX)=trcsa_soHml(idsa_H0PO4B,L,NY,NX)
    trcsa_soHml2(idsa_H3PO4B,L,NY,NX)=trcsa_soHml(idsa_H3PO4B,L,NY,NX)
    trcsa_soHml2(idsa_FeHPO4B,L,NY,NX)=trcsa_soHml(idsa_FeHPO4B,L,NY,NX)
    trcsa_soHml2(idsa_FeH2PO4B,L,NY,NX)=trcsa_soHml(idsa_FeH2PO4B,L,NY,NX)
    trcsa_soHml2(idsa_CaPO4B,L,NY,NX)=trcsa_soHml(idsa_CaPO4B,L,NY,NX)
    trcsa_soHml2(idsa_CaHPO4B,L,NY,NX)=trcsa_soHml(idsa_CaHPO4B,L,NY,NX)
    trcsa_soHml2(idsa_CaH2PO4B,L,NY,NX)=trcsa_soHml(idsa_CaH2PO4B,L,NY,NX)
    trcsa_soHml2(idsa_MgHPO4B,L,NY,NX)=trcsa_soHml(idsa_MgHPO4B,L,NY,NX)
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

  subroutine ZeroBoundarySnowFlux(N,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,M5,M4
!     begin_execution
  trcsa_RQ(idsa_beg:idsa_end,N,M5,M4)=0.0

  end subroutine ZeroBoundarySnowFlux
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteFluxFromRecharge(N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,NN,M5,M4
!     begin_execution

  trcsa_RQR(idsa_beg:idsa_end,N,NN,M5,M4)=0.0

  end subroutine ZeroSoluteFluxFromRecharge
!------------------------------------------------------------------------------------------

  subroutine SoluteExportThruBoundary(N1,N2,M,N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N1,N2,M,N,NN,M5,M4
  real(r8) :: FQRM
  integer :: NTSA
!     begin_execution

  FQRM=QRMN(M,N,NN,M5,M4)/QRM(M,N2,N1)
  trcsa_RQR(idsa_Al,N,NN,M5,M4)=trcsa_RQR0(idsa_Al,N2,N1)*FQRM
  trcsa_RQR(idsa_Fe,N,NN,M5,M4)=trcsa_RQR0(idsa_Fe,N2,N1)*FQRM
  trcsa_RQR(idsa_Hp,N,NN,M5,M4)=trcsa_RQR0(idsa_Hp,N2,N1)*FQRM
  trcsa_RQR(idsa_Ca,N,NN,M5,M4)=trcsa_RQR0(idsa_Ca,N2,N1)*FQRM
  trcsa_RQR(idsa_Mg,N,NN,M5,M4)=trcsa_RQR0(idsa_Mg,N2,N1)*FQRM
  trcsa_RQR(idsa_Na,N,NN,M5,M4)=trcsa_RQR0(idsa_Na,N2,N1)*FQRM
  trcsa_RQR(idsa_K,N,NN,M5,M4)=trcsa_RQR0(idsa_K,N2,N1)*FQRM
  trcsa_RQR(idsa_OH,N,NN,M5,M4)=trcsa_RQR0(idsa_OH,N2,N1)*FQRM
  trcsa_RQR(idsa_SO4,N,NN,M5,M4)=trcsa_RQR0(idsa_SO4,N2,N1)*FQRM
  trcsa_RQR(idsa_Cl,N,NN,M5,M4)=trcsa_RQR0(idsa_Cl,N2,N1)*FQRM
  trcsa_RQR(idsa_CO3,N,NN,M5,M4)=trcsa_RQR0(idsa_CO3,N2,N1)*FQRM
  trcsa_RQR(idsa_HCO3,N,NN,M5,M4)=trcsa_RQR0(idsa_HCO3,N2,N1)*FQRM
  trcsa_RQR(idsa_AlOH,N,NN,M5,M4)=trcsa_RQR0(idsa_AlOH,N2,N1)*FQRM
  trcsa_RQR(idsa_AlOH2,N,NN,M5,M4)=trcsa_RQR0(idsa_AlOH2,N2,N1)*FQRM
  trcsa_RQR(idsa_AlOH3,N,NN,M5,M4)=trcsa_RQR0(idsa_AlOH3,N2,N1)*FQRM
  trcsa_RQR(idsa_AlOH4,N,NN,M5,M4)=trcsa_RQR0(idsa_AlOH4,N2,N1)*FQRM
  trcsa_RQR(idsa_AlSO4,N,NN,M5,M4)=trcsa_RQR0(idsa_AlSO4,N2,N1)*FQRM
  trcsa_RQR(idsa_FeOH,N,NN,M5,M4)=trcsa_RQR0(idsa_FeOH,N2,N1)*FQRM
  trcsa_RQR(idsa_FeOH2,N,NN,M5,M4)=trcsa_RQR0(idsa_FeOH2,N2,N1)*FQRM
  trcsa_RQR(idsa_FeOH3,N,NN,M5,M4)=trcsa_RQR0(idsa_FeOH3,N2,N1)*FQRM
  trcsa_RQR(idsa_FeOH4,N,NN,M5,M4)=trcsa_RQR0(idsa_FeOH4,N2,N1)*FQRM
  trcsa_RQR(idsa_FeSO4,N,NN,M5,M4)=trcsa_RQR0(idsa_FeSO4,N2,N1)*FQRM
  trcsa_RQR(idsa_CaOH2,N,NN,M5,M4)=trcsa_RQR0(idsa_CaOH2,N2,N1)*FQRM
  trcsa_RQR(idsa_CaCO3,N,NN,M5,M4)=trcsa_RQR0(idsa_CaCO3,N2,N1)*FQRM
  trcsa_RQR(idsa_CaHCO3,N,NN,M5,M4)=trcsa_RQR0(idsa_CaHCO3,N2,N1)*FQRM
  trcsa_RQR(idsa_CaSO4,N,NN,M5,M4)=trcsa_RQR0(idsa_CaSO4,N2,N1)*FQRM
  trcsa_RQR(idsa_MgOH2,N,NN,M5,M4)=trcsa_RQR0(idsa_MgOH2,N2,N1)*FQRM
  trcsa_RQR(idsa_MgCO3,N,NN,M5,M4)=trcsa_RQR0(idsa_MgCO3,N2,N1)*FQRM
  trcsa_RQR(idsa_MgHCO3,N,NN,M5,M4)=trcsa_RQR0(idsa_MgHCO3,N2,N1)*FQRM
  trcsa_RQR(idsa_MgSO4,N,NN,M5,M4)=trcsa_RQR0(idsa_MgSO4,N2,N1)*FQRM
  trcsa_RQR(idsa_NaCO3,N,NN,M5,M4)=trcsa_RQR0(idsa_NaCO3,N2,N1)*FQRM
  trcsa_RQR(idsa_NaSO4,N,NN,M5,M4)=trcsa_RQR0(idsa_NaSO4,N2,N1)*FQRM
  trcsa_RQR(idsa_KSO4,N,NN,M5,M4)=trcsa_RQR0(idsa_KSO4,N2,N1)*FQRM
  trcsa_RQR(idsa_H0PO4,N,NN,M5,M4)=trcsa_RQR0(idsa_H0PO4,N2,N1)*FQRM
  trcsa_RQR(idsa_H3PO4,N,NN,M5,M4)=trcsa_RQR0(idsa_H3PO4,N2,N1)*FQRM
  trcsa_RQR(idsa_FeHPO4,N,NN,M5,M4)=trcsa_RQR0(idsa_FeHPO4,N2,N1)*FQRM
  trcsa_RQR(idsa_FeH2PO4,N,NN,M5,M4)=trcsa_RQR0(idsa_FeH2PO4,N2,N1)*FQRM
  trcsa_RQR(idsa_CaPO4,N,NN,M5,M4)=trcsa_RQR0(idsa_CaPO4,N2,N1)*FQRM
  trcsa_RQR(idsa_CaHPO4,N,NN,M5,M4)=trcsa_RQR0(idsa_CaHPO4,N2,N1)*FQRM
  trcsa_RQR(idsa_CaH2PO4,N,NN,M5,M4)=trcsa_RQR0(idsa_CaH2PO4,N2,N1)*FQRM
  trcsa_RQR(idsa_MgHPO4,N,NN,M5,M4)=trcsa_RQR0(idsa_MgHPO4,N2,N1)*FQRM
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QRS=hourly solute in runoff
!     RQR*=solute in runoff
!
  DO NTSA=idsa_beg,idsa_end
    trcsa_XQR(NTSA,N,NN,M5,M4)=trcsa_XQR(NTSA,N,NN,M5,M4)+trcsa_RQR(NTSA,N,NN,M5,M4)
  ENDDO
  end subroutine SoluteExportThruBoundary
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,NN,M5,M4
!     begin_execution

  trcsa_RQR(idsa_beg:idsa_end,N,NN,M5,M4)=0.0

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
  RALFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Al,M3,M2,M1))
  RFEFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Fe,M3,M2,M1))
  RHYFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Hp,M3,M2,M1))
  RCAFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Ca,M3,M2,M1))
  RMGFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Mg,M3,M2,M1))
  RNAFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Na,M3,M2,M1))
  RKAFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_K,M3,M2,M1))
  ROHFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_OH,M3,M2,M1))
  RSOFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_SO4,M3,M2,M1))
  RCLFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_Cl,M3,M2,M1))
  RC3FHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CO3,M3,M2,M1))
  RHCFHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_HCO3,M3,M2,M1))
  RAL1HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH,M3,M2,M1))
  RAL2HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH2,M3,M2,M1))
  RAL3HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH3,M3,M2,M1))
  RAL4HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH4,M3,M2,M1))
  RALSHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_AlSO4,M3,M2,M1))
  RFE1HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH,M3,M2,M1))
  RFE2HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH2,M3,M2,M1))
  RFE3HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH3,M3,M2,M1))
  RFE4HS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH4,M3,M2,M1))
  RFESHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeSO4,M3,M2,M1))
  RCAOHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaOH2,M3,M2,M1))
  RCACHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaCO3,M3,M2,M1))
  RCAHHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHCO3,M3,M2,M1))
  RCASHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaSO4,M3,M2,M1))
  RMGOHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_MgOH2,M3,M2,M1))
  RMGCHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_MgCO3,M3,M2,M1))
  RMGHHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHCO3,M3,M2,M1))
  RMGSHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_MgSO4,M3,M2,M1))
  RNACHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_NaCO3,M3,M2,M1))
  RNASHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_NaSO4,M3,M2,M1))
  RKASHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_KSO4,M3,M2,M1))
  RH0PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RH3PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RF1PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RF2PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RC0PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RC1PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RC2PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RM1PHS(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4,M3,M2,M1))*trcs_VLN(ids_H1PO4,M3,M2,M1)
  RH0BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RH3BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RF1BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RF2BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RC0BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RC1BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RC2BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
  RM1BHB(N,M6,M5,M4)=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4B,M3,M2,M1))*trcs_VLN(ids_H1PO4B,M3,M2,M1)
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
      integer :: NN,NTSA
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
    DO NTSA=idsa_beg,idsa_end
      trcsa_TQR(NTSA,N2,N1)=trcsa_TQR(NTSA,N2,N1)+trcsa_RQR(NTSA,N,NN,N2,N1)
    ENDDO

    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO NTSA=idsa_beg,idsa_end
        trcsa_TQR(NTSA,N2,N1)=trcsa_TQR(NTSA,N2,N1)-trcsa_RQR(NTSA,N,NN,N5,N4)
      ENDDO

    ENDIF
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO NTSA=idsa_beg,idsa_end
        trcsa_TQR(NTSA,N2,N1)=trcsa_TQR(NTSA,N2,N1)-trcsa_RQR(NTSA,N,NN,N5B,N4B)
      ENDDO
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

    integer :: NTSA
!     begin_execution
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
    DO NTSA=idsa_beg,idsa_end
      trcsa_TQ(NTSA,N2,N1)=trcsa_TQ(NTSA,N2,N1)+trcsa_RQ(NTSA,N,N2,N1)-trcsa_RQ(NTSA,N,N5,N4)
    ENDDO
    end subroutine NetOverloadFLuxInSnow
!------------------------------------------------------------------------------------------

  subroutine NetFluxInSnowpack(M,NY,NX,N1,N2)
!
!     Description:
!
  integer, intent(in) :: M,NY,NX,N1,N2
  integer :: LS,LS2,NTSA
!     begin_execution

  D1205: DO LS=1,JS
    IF(VHCPWM(M,LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
  !
    IF(LS.LT.JS.AND.VHCPWM(M,LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
      DO NTSA=idsa_beg,idsa_end
        trcsa_TBLS(NTSA,LS,NY,NX)=trcsa_TBLS(NTSA,LS,NY,NX)+trcsa_RBLS(NTSA,LS,NY,NX)-trcsa_RBLS(NTSA,LS2,NY,NX)
      ENDDO
    ELSE
  !
  !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
  !
      trcsa_TBLS(idsa_Al,LS,NY,NX)=trcsa_TBLS(idsa_Al,LS,NY,NX)+trcsa_RBLS(idsa_Al,LS,NY,NX) &
        -trcsa_RFLS(idsa_Al,3,0,N2,N1)-trcsa_RFLS(idsa_Al,3,NUM(N2,N1),N2,N1)
  !    3-RALFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Fe,LS,NY,NX)=trcsa_TBLS(idsa_Fe,LS,NY,NX)+trcsa_RBLS(idsa_Fe,LS,NY,NX) &
        -trcsa_RFLS(idsa_Fe,3,0,N2,N1)-trcsa_RFLS(idsa_Fe,3,NUM(N2,N1),N2,N1)
  !    3-RFEFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Hp,LS,NY,NX)=trcsa_TBLS(idsa_Hp,LS,NY,NX)+trcsa_RBLS(idsa_Hp,LS,NY,NX) &
        -trcsa_RFLS(idsa_Hp,3,0,N2,N1)-trcsa_RFLS(idsa_Hp,3,NUM(N2,N1),N2,N1)
  !    3-RHYFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Ca,LS,NY,NX)=trcsa_TBLS(idsa_Ca,LS,NY,NX)+trcsa_RBLS(idsa_Ca,LS,NY,NX) &
        -trcsa_RFLS(idsa_Ca,3,0,N2,N1)-trcsa_RFLS(idsa_Ca,3,NUM(N2,N1),N2,N1)
  !    3-RCAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Mg,LS,NY,NX)=trcsa_TBLS(idsa_Mg,LS,NY,NX)+trcsa_RBLS(idsa_Mg,LS,NY,NX) &
        -trcsa_RFLS(idsa_Mg,3,0,N2,N1)-trcsa_RFLS(idsa_Mg,3,NUM(N2,N1),N2,N1)
  !    3-RMGFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Na,LS,NY,NX)=trcsa_TBLS(idsa_Na,LS,NY,NX)+trcsa_RBLS(idsa_Na,LS,NY,NX) &
        -trcsa_RFLS(idsa_Na,3,0,N2,N1)-trcsa_RFLS(idsa_Na,3,NUM(N2,N1),N2,N1)
  !    3-RNAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_K,LS,NY,NX)=trcsa_TBLS(idsa_K,LS,NY,NX)+trcsa_RBLS(idsa_K,LS,NY,NX) &
        -trcsa_RFLS(idsa_K,3,0,N2,N1)-trcsa_RFLS(idsa_K,3,NUM(N2,N1),N2,N1)
  !    3-RKAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_OH,LS,NY,NX)=trcsa_TBLS(idsa_OH,LS,NY,NX)+trcsa_RBLS(idsa_OH,LS,NY,NX) &
        -trcsa_RFLS(idsa_OH,3,0,N2,N1)-trcsa_RFLS(idsa_OH,3,NUM(N2,N1),N2,N1)
  !    3-ROHFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_SO4,LS,NY,NX)=trcsa_TBLS(idsa_SO4,LS,NY,NX)+trcsa_RBLS(idsa_SO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_SO4,3,0,N2,N1)-trcsa_RFLS(idsa_SO4,3,NUM(N2,N1),N2,N1)
  !    3-RSOFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Cl,LS,NY,NX)=trcsa_TBLS(idsa_Cl,LS,NY,NX)+trcsa_RBLS(idsa_Cl,LS,NY,NX) &
        -trcsa_RFLS(idsa_Cl,3,0,N2,N1)-trcsa_RFLS(idsa_Cl,3,NUM(N2,N1),N2,N1)
  !    3-RCLFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CO3,LS,NY,NX)=trcsa_TBLS(idsa_CO3,LS,NY,NX)+trcsa_RBLS(idsa_CO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_CO3,3,0,N2,N1)-trcsa_RFLS(idsa_CO3,3,NUM(N2,N1),N2,N1)
  !    3-RC3FHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_HCO3,LS,NY,NX)=trcsa_TBLS(idsa_HCO3,LS,NY,NX)+trcsa_RBLS(idsa_HCO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_HCO3,3,0,N2,N1)-trcsa_RFLS(idsa_HCO3,3,NUM(N2,N1),N2,N1)
  !    3-RHCFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH,LS,NY,NX)=trcsa_TBLS(idsa_AlOH,LS,NY,NX)+trcsa_RBLS(idsa_AlOH,LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH,3,NUM(N2,N1),N2,N1)
  !    3-RAL1HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH2,LS,NY,NX)=trcsa_TBLS(idsa_AlOH2,LS,NY,NX)+trcsa_RBLS(idsa_AlOH2,LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH2,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH2,3,NUM(N2,N1),N2,N1)
  !    3-RAL2HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH3,LS,NY,NX)=trcsa_TBLS(idsa_AlOH3,LS,NY,NX)+trcsa_RBLS(idsa_AlOH3,LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH3,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH3,3,NUM(N2,N1),N2,N1)
  !    3-RAL3HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH4,LS,NY,NX)=trcsa_TBLS(idsa_AlOH4,LS,NY,NX)+trcsa_RBLS(idsa_AlOH4,LS,NY,NX) &
        -trcsa_RFLS(idsa_AlOH4,3,0,N2,N1)-trcsa_RFLS(idsa_AlOH4,3,NUM(N2,N1),N2,N1)
  !    3-RAL4HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlSO4,LS,NY,NX)=trcsa_TBLS(idsa_AlSO4,LS,NY,NX)+trcsa_RBLS(idsa_AlSO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_AlSO4,3,0,N2,N1)-trcsa_RFLS(idsa_AlSO4,3,NUM(N2,N1),N2,N1)
  !    3-RALSHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH,LS,NY,NX)=trcsa_TBLS(idsa_FeOH,LS,NY,NX)+trcsa_RBLS(idsa_FeOH,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH,3,NUM(N2,N1),N2,N1)
  !    3-RFE1HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH2,LS,NY,NX)=trcsa_TBLS(idsa_FeOH2,LS,NY,NX)+trcsa_RBLS(idsa_FeOH2,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH2,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH2,3,NUM(N2,N1),N2,N1)
  !    3-RFE2HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH3,LS,NY,NX)=trcsa_TBLS(idsa_FeOH3,LS,NY,NX)+trcsa_RBLS(idsa_FeOH3,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH3,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH3,3,NUM(N2,N1),N2,N1)
  !    3-RFE3HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH4,LS,NY,NX)=trcsa_TBLS(idsa_FeOH4,LS,NY,NX)+trcsa_RBLS(idsa_FeOH4,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeOH4,3,0,N2,N1)-trcsa_RFLS(idsa_FeOH4,3,NUM(N2,N1),N2,N1)
  !    3-RFE4HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeSO4,LS,NY,NX)=trcsa_TBLS(idsa_FeSO4,LS,NY,NX)+trcsa_RBLS(idsa_FeSO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeSO4,3,0,N2,N1)-trcsa_RFLS(idsa_FeSO4,3,NUM(N2,N1),N2,N1)
  !    3-RFESHS(3,NUM(N2,N1),N2,N1)
        trcsa_TBLS(idsa_CaOH2,LS,NY,NX)=trcsa_TBLS(idsa_CaOH2,LS,NY,NX)+trcsa_RBLS(idsa_CaOH2,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaOH2,3,0,N2,N1)-trcsa_RFLS(idsa_CaOH2,3,NUM(N2,N1),N2,N1)
  !    3-RCAOHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaCO3,LS,NY,NX)+trcsa_RBLS(idsa_CaCO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaCO3,3,0,N2,N1)-trcsa_RFLS(idsa_CaCO3,3,NUM(N2,N1),N2,N1)
  !    3-RCACHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)+trcsa_RBLS(idsa_CaHCO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaHCO3,3,0,N2,N1)-trcsa_RFLS(idsa_CaHCO3,3,NUM(N2,N1),N2,N1)
  !    3-RCAHHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaSO4,LS,NY,NX)=trcsa_TBLS(idsa_CaSO4,LS,NY,NX)+trcsa_RBLS(idsa_CaSO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaSO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaSO4,3,NUM(N2,N1),N2,N1)
  !    3-RCASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgOH2,LS,NY,NX)=trcsa_TBLS(idsa_MgOH2,LS,NY,NX)+trcsa_RBLS(idsa_MgOH2,LS,NY,NX) &
        -trcsa_RFLS(idsa_MgOH2,3,0,N2,N1)-trcsa_RFLS(idsa_MgOH2,3,NUM(N2,N1),N2,N1)
  !    3-RMGOHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgCO3,LS,NY,NX)+trcsa_RBLS(idsa_MgCO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_MgCO3,3,0,N2,N1)-trcsa_RFLS(idsa_MgCO3,3,NUM(N2,N1),N2,N1)
  !    3-RMGCHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)+trcsa_RBLS(idsa_MgHCO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_MgHCO3,3,0,N2,N1)-trcsa_RFLS(idsa_MgHCO3,3,NUM(N2,N1),N2,N1)
  !    3-RMGHHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgSO4,LS,NY,NX)=trcsa_TBLS(idsa_MgSO4,LS,NY,NX)+trcsa_RBLS(idsa_MgSO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_MgSO4,3,0,N2,N1)-trcsa_RFLS(idsa_MgSO4,3,NUM(N2,N1),N2,N1)
  !    3-RMGSHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_NaCO3,LS,NY,NX)=trcsa_TBLS(idsa_NaCO3,LS,NY,NX)+trcsa_RBLS(idsa_NaCO3,LS,NY,NX) &
        -trcsa_RFLS(idsa_NaCO3,3,0,N2,N1)-trcsa_RFLS(idsa_NaCO3,3,NUM(N2,N1),N2,N1)
  !    3-RNACHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_NaSO4,LS,NY,NX)=trcsa_TBLS(idsa_NaSO4,LS,NY,NX)+trcsa_RBLS(idsa_NaSO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_NaSO4,3,0,N2,N1)-trcsa_RFLS(idsa_NaSO4,3,NUM(N2,N1),N2,N1)
  !    3-RNASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_KSO4,LS,NY,NX)=trcsa_TBLS(idsa_KSO4,LS,NY,NX)+trcsa_RBLS(idsa_KSO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_KSO4,3,0,N2,N1)-trcsa_RFLS(idsa_KSO4,3,NUM(N2,N1),N2,N1)
  !    3-RKASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_H0PO4,LS,NY,NX)=trcsa_TBLS(idsa_H0PO4,LS,NY,NX)+trcsa_RBLS(idsa_H0PO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_H0PO4,3,0,N2,N1)-trcsa_RFLS(idsa_H0PO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_H0PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RH0PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_H0PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_H3PO4,LS,NY,NX)=trcsa_TBLS(idsa_H3PO4,LS,NY,NX)+trcsa_RBLS(idsa_H3PO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_H3PO4,3,0,N2,N1)-RH3PHS(3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_H3PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RH3PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_H3PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)=trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)+trcsa_RBLS(idsa_FeHPO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeHPO4,3,0,N2,N1)-trcsa_RFLS(idsa_FeHPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_FeHPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RF1PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_FeHPO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)+trcsa_RBLS(idsa_FeH2PO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_FeH2PO4,3,0,N2,N1)-trcsa_RFLS(idsa_FeH2PO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_FeH2PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RF2PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_FeH2PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaPO4,LS,NY,NX)+trcsa_RBLS(idsa_CaPO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaPO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_CaPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RC0PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_CaPO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)+trcsa_RBLS(idsa_CaHPO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaHPO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaHPO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_CaHPO4B,3,NUM(N2,N1),N2,N1)
  !    3-RC1PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_CaHPO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)+trcsa_RBLS(idsa_CaH2PO4,LS,NY,NX) &
        -trcsa_RFLS(idsa_CaH2PO4,3,0,N2,N1)-trcsa_RFLS(idsa_CaH2PO4,3,NUM(N2,N1),N2,N1) &
        -trcsa_RFLS(idsa_CaH2PO4B,3,NUM(N2,N1),N2,N1)
  !    3-RC2PHS(3,NUM(N2,N1),N2,N1)-trcsa_RFLS(idsa_CaH2PO4B,3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)=trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)+trcsa_RBLS(idsa_MgHPO4,LS,NY,NX) &
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
  integer :: LL,NTSA
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
      DO NTSA=idsa_beg,idsab_end
        trcsa_TFLS(NTSA,N3,N2,N1)=trcsa_TFLS(NTSA,N3,N2,N1)+trcsa_RFLS(NTSA,N,N3,N2,N1)-trcsa_RFLS(NTSA,N,N6,N5,N4)
      ENDDO

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
      trcsa_TFLS(idsa_beg:idsab_end,N3,N2,N1)=0.0
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
                call ZeroBoundarySnowFlux(N,M5,M4)
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
  DO NTSA=idsa_beg,idsa_end
    trcsa_solml2(NTSA,0,NY,NX)=trcsa_solml2(NTSA,0,NY,NX)+trcsa_TQR(NTSA,NY,NX)+trcsa_RFLS(NTSA,3,0,NY,NX)
  ENDDO
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
      trcsa_solml2(idsa_Al,L,NY,NX)=trcsa_solml2(idsa_Al,L,NY,NX)+trcsa_TFLS(idsa_Al,L,NY,NX)+trcsa_RFXS(idsa_Al,L,NY,NX)+trcsa_RFLZ(idsa_Al,L,NY,NX)
      trcsa_solml2(idsa_Fe,L,NY,NX)=trcsa_solml2(idsa_Fe,L,NY,NX)+trcsa_TFLS(idsa_Fe,L,NY,NX)+trcsa_RFXS(idsa_Fe,L,NY,NX)+trcsa_RFLZ(idsa_Fe,L,NY,NX)
      trcsa_solml2(idsa_Hp,L,NY,NX)=trcsa_solml2(idsa_Hp,L,NY,NX)+trcsa_TFLS(idsa_Hp,L,NY,NX)+trcsa_RFXS(idsa_Hp,L,NY,NX)+trcsa_RFLZ(idsa_Hp,L,NY,NX)
      trcsa_solml2(idsa_Ca,L,NY,NX)=trcsa_solml2(idsa_Ca,L,NY,NX)+trcsa_TFLS(idsa_Ca,L,NY,NX)+trcsa_RFXS(idsa_Ca,L,NY,NX)+trcsa_RFLZ(idsa_Ca,L,NY,NX)
      trcsa_solml2(idsa_Mg,L,NY,NX)=trcsa_solml2(idsa_Mg,L,NY,NX)+trcsa_TFLS(idsa_Mg,L,NY,NX)+trcsa_RFXS(idsa_Mg,L,NY,NX)+trcsa_RFLZ(idsa_Mg,L,NY,NX)
      trcsa_solml2(idsa_Na,L,NY,NX)=trcsa_solml2(idsa_Na,L,NY,NX)+trcsa_TFLS(idsa_Na,L,NY,NX)+trcsa_RFXS(idsa_Na,L,NY,NX)+trcsa_RFLZ(idsa_Na,L,NY,NX)
      trcsa_solml2(idsa_K,L,NY,NX)=trcsa_solml2(idsa_K,L,NY,NX)+trcsa_TFLS(idsa_K,L,NY,NX)+trcsa_RFXS(idsa_K,L,NY,NX)+trcsa_RFLZ(idsa_K,L,NY,NX)
      trcsa_solml2(idsa_OH,L,NY,NX)=trcsa_solml2(idsa_OH,L,NY,NX)+trcsa_TFLS(idsa_OH,L,NY,NX)+trcsa_RFXS(idsa_OH,L,NY,NX)+trcsa_RFLZ(idsa_OH,L,NY,NX)
      trcsa_solml2(idsa_SO4,L,NY,NX)=trcsa_solml2(idsa_SO4,L,NY,NX)+trcsa_TFLS(idsa_SO4,L,NY,NX)+trcsa_RFXS(idsa_SO4,L,NY,NX)+trcsa_RFLZ(idsa_SO4,L,NY,NX)
      trcsa_solml2(idsa_Cl,L,NY,NX)=trcsa_solml2(idsa_Cl,L,NY,NX)+trcsa_TFLS(idsa_Cl,L,NY,NX)+trcsa_RFXS(idsa_Cl,L,NY,NX)+trcsa_RFLZ(idsa_Cl,L,NY,NX)
      trcsa_solml2(idsa_CO3,L,NY,NX)=trcsa_solml2(idsa_CO3,L,NY,NX)+trcsa_TFLS(idsa_CO3,L,NY,NX)+trcsa_RFXS(idsa_CO3,L,NY,NX)+trcsa_RFLZ(idsa_CO3,L,NY,NX)
      trcsa_solml2(idsa_HCO3,L,NY,NX)=trcsa_solml2(idsa_HCO3,L,NY,NX)+trcsa_TFLS(idsa_HCO3,L,NY,NX)+trcsa_RFXS(idsa_HCO3,L,NY,NX)+trcsa_RFLZ(idsa_HCO3,L,NY,NX)
      trcsa_solml2(idsa_AlOH,L,NY,NX)=trcsa_solml2(idsa_AlOH,L,NY,NX)+trcsa_TFLS(idsa_AlOH,L,NY,NX)+trcsa_RFXS(idsa_AlOH,L,NY,NX)+trcsa_RFLZ(idsa_AlOH,L,NY,NX)
      trcsa_solml2(idsa_AlOH2,L,NY,NX)=trcsa_solml2(idsa_AlOH2,L,NY,NX)+trcsa_TFLS(idsa_AlOH2,L,NY,NX)+trcsa_RFXS(idsa_AlOH2,L,NY,NX)+trcsa_RFLZ(idsa_AlOH2,L,NY,NX)
      trcsa_solml2(idsa_AlOH3,L,NY,NX)=trcsa_solml2(idsa_AlOH3,L,NY,NX)+trcsa_TFLS(idsa_AlOH3,L,NY,NX)+trcsa_RFXS(idsa_AlOH3,L,NY,NX)+trcsa_RFLZ(idsa_AlOH3,L,NY,NX)
      trcsa_solml2(idsa_AlOH4,L,NY,NX)=trcsa_solml2(idsa_AlOH4,L,NY,NX)+trcsa_TFLS(idsa_AlOH4,L,NY,NX)+trcsa_RFXS(idsa_AlOH4,L,NY,NX)+trcsa_RFLZ(idsa_AlOH4,L,NY,NX)
      trcsa_solml2(idsa_AlSO4,L,NY,NX)=trcsa_solml2(idsa_AlSO4,L,NY,NX)+trcsa_TFLS(idsa_AlSO4,L,NY,NX)+trcsa_RFXS(idsa_AlSO4,L,NY,NX)+trcsa_RFLZ(idsa_AlSO4,L,NY,NX)
      trcsa_solml2(idsa_FeOH,L,NY,NX)=trcsa_solml2(idsa_FeOH,L,NY,NX)+trcsa_TFLS(idsa_FeOH,L,NY,NX)+trcsa_RFXS(idsa_FeOH,L,NY,NX)+trcsa_RFLZ(idsa_FeOH,L,NY,NX)
      trcsa_solml2(idsa_FeOH2,L,NY,NX)=trcsa_solml2(idsa_FeOH2,L,NY,NX)+trcsa_TFLS(idsa_FeOH2,L,NY,NX)+trcsa_RFXS(idsa_FeOH2,L,NY,NX)+trcsa_RFLZ(idsa_FeOH2,L,NY,NX)
      trcsa_solml2(idsa_FeOH3,L,NY,NX)=trcsa_solml2(idsa_FeOH3,L,NY,NX)+trcsa_TFLS(idsa_FeOH3,L,NY,NX)+trcsa_RFXS(idsa_FeOH3,L,NY,NX)+trcsa_RFLZ(idsa_FeOH3,L,NY,NX)
      trcsa_solml2(idsa_FeOH4,L,NY,NX)=trcsa_solml2(idsa_FeOH4,L,NY,NX)+trcsa_TFLS(idsa_FeOH4,L,NY,NX)+trcsa_RFXS(idsa_FeOH4,L,NY,NX)+trcsa_RFLZ(idsa_FeOH4,L,NY,NX)
      trcsa_solml2(idsa_FeSO4,L,NY,NX)=trcsa_solml2(idsa_FeSO4,L,NY,NX)+trcsa_TFLS(idsa_FeSO4,L,NY,NX)+trcsa_RFXS(idsa_FeSO4,L,NY,NX)+trcsa_RFLZ(idsa_FeSO4,L,NY,NX)
      trcsa_solml2(idsa_CaOH2,L,NY,NX)=trcsa_solml2(idsa_CaOH2,L,NY,NX)+trcsa_TFLS(idsa_CaOH2,L,NY,NX)+trcsa_RFXS(idsa_CaOH2,L,NY,NX)+trcsa_RFLZ(idsa_CaOH2,L,NY,NX)
      trcsa_solml2(idsa_CaCO3,L,NY,NX)=trcsa_solml2(idsa_CaCO3,L,NY,NX)+trcsa_TFLS(idsa_CaCO3,L,NY,NX)+trcsa_RFXS(idsa_CaCO3,L,NY,NX)+trcsa_RFLZ(idsa_CaCO3,L,NY,NX)
      trcsa_solml2(idsa_CaHCO3,L,NY,NX)=trcsa_solml2(idsa_CaHCO3,L,NY,NX)+trcsa_TFLS(idsa_CaHCO3,L,NY,NX)+trcsa_RFXS(idsa_CaHCO3,L,NY,NX)+trcsa_RFLZ(idsa_CaHCO3,L,NY,NX)
      trcsa_solml2(idsa_CaSO4,L,NY,NX)=trcsa_solml2(idsa_CaSO4,L,NY,NX)+trcsa_TFLS(idsa_CaSO4,L,NY,NX)+trcsa_RFXS(idsa_CaSO4,L,NY,NX)+trcsa_RFLZ(idsa_CaSO4,L,NY,NX)
      trcsa_solml2(idsa_MgOH2,L,NY,NX)=trcsa_solml2(idsa_MgOH2,L,NY,NX)+trcsa_TFLS(idsa_MgOH2,L,NY,NX)+trcsa_RFXS(idsa_MgOH2,L,NY,NX)+trcsa_RFLZ(idsa_MgOH2,L,NY,NX)
      trcsa_solml2(idsa_MgCO3,L,NY,NX)=trcsa_solml2(idsa_MgCO3,L,NY,NX)+trcsa_TFLS(idsa_MgCO3,L,NY,NX)+trcsa_RFXS(idsa_MgCO3,L,NY,NX)+trcsa_RFLZ(idsa_MgCO3,L,NY,NX)
      trcsa_solml2(idsa_MgHCO3,L,NY,NX)=trcsa_solml2(idsa_MgHCO3,L,NY,NX)+trcsa_TFLS(idsa_MgHCO3,L,NY,NX)+trcsa_RFXS(idsa_MgHCO3,L,NY,NX)+trcsa_RFLZ(idsa_MgHCO3,L,NY,NX)
      trcsa_solml2(idsa_MgSO4,L,NY,NX)=trcsa_solml2(idsa_MgSO4,L,NY,NX)+trcsa_TFLS(idsa_MgSO4,L,NY,NX)+trcsa_RFXS(idsa_MgSO4,L,NY,NX)+trcsa_RFLZ(idsa_MgSO4,L,NY,NX)
      trcsa_solml2(idsa_NaCO3,L,NY,NX)=trcsa_solml2(idsa_NaCO3,L,NY,NX)+trcsa_TFLS(idsa_NaCO3,L,NY,NX)+trcsa_RFXS(idsa_NaCO3,L,NY,NX)+trcsa_RFLZ(idsa_NaCO3,L,NY,NX)
      trcsa_solml2(idsa_NaSO4,L,NY,NX)=trcsa_solml2(idsa_NaSO4,L,NY,NX)+trcsa_TFLS(idsa_NaSO4,L,NY,NX)+trcsa_RFXS(idsa_NaSO4,L,NY,NX)+trcsa_RFLZ(idsa_NaSO4,L,NY,NX)
      trcsa_solml2(idsa_KSO4,L,NY,NX)=trcsa_solml2(idsa_KSO4,L,NY,NX)+trcsa_TFLS(idsa_KSO4,L,NY,NX)+trcsa_RFXS(idsa_KSO4,L,NY,NX)+trcsa_RFLZ(idsa_KSO4,L,NY,NX)
      trcsa_solml2(idsa_H0PO4,L,NY,NX)=trcsa_solml2(idsa_H0PO4,L,NY,NX)+trcsa_TFLS(idsa_H0PO4,L,NY,NX)+trcsa_RFXS(idsa_H0PO4,L,NY,NX)+trcsa_RFLZ(idsa_H0PO4,L,NY,NX)
      trcsa_solml2(idsa_H3PO4,L,NY,NX)=trcsa_solml2(idsa_H3PO4,L,NY,NX)+trcsa_TFLS(idsa_H3PO4,L,NY,NX)+trcsa_RFXS(idsa_H3PO4,L,NY,NX)+trcsa_RFLZ(idsa_H3PO4,L,NY,NX)
      trcsa_solml2(idsa_FeHPO4,L,NY,NX)=trcsa_solml2(idsa_FeHPO4,L,NY,NX)+trcsa_TFLS(idsa_FeHPO4,L,NY,NX)+trcsa_RFXS(idsa_FeHPO4,L,NY,NX)+trcsa_RFLZ(idsa_FeHPO4,L,NY,NX)
      trcsa_solml2(idsa_FeH2PO4,L,NY,NX)=trcsa_solml2(idsa_FeH2PO4,L,NY,NX)+trcsa_TFLS(idsa_FeH2PO4,L,NY,NX)+trcsa_RFXS(idsa_FeH2PO4,L,NY,NX)+trcsa_RFLZ(idsa_FeH2PO4,L,NY,NX)
      trcsa_solml2(idsa_CaPO4,L,NY,NX)=trcsa_solml2(idsa_CaPO4,L,NY,NX)+trcsa_TFLS(idsa_CaPO4,L,NY,NX)+trcsa_RFXS(idsa_CaPO4,L,NY,NX)+trcsa_RFLZ(idsa_CaPO4,L,NY,NX)
      trcsa_solml2(idsa_CaHPO4,L,NY,NX)=trcsa_solml2(idsa_CaHPO4,L,NY,NX)+trcsa_TFLS(idsa_CaHPO4,L,NY,NX)+trcsa_RFXS(idsa_CaHPO4,L,NY,NX)+trcsa_RFLZ(idsa_CaHPO4,L,NY,NX)
      trcsa_solml2(idsa_CaH2PO4,L,NY,NX)=trcsa_solml2(idsa_CaH2PO4,L,NY,NX)+trcsa_TFLS(idsa_CaH2PO4,L,NY,NX)+trcsa_RFXS(idsa_CaH2PO4,L,NY,NX)+trcsa_RFLZ(idsa_CaH2PO4,L,NY,NX)
      trcsa_solml2(idsa_MgHPO4,L,NY,NX)=trcsa_solml2(idsa_MgHPO4,L,NY,NX)+trcsa_TFLS(idsa_MgHPO4,L,NY,NX)+trcsa_RFXS(idsa_MgHPO4,L,NY,NX)+trcsa_RFLZ(idsa_MgHPO4,L,NY,NX)
      trcsa_solml2(idsa_H0PO4B,L,NY,NX)=trcsa_solml2(idsa_H0PO4B,L,NY,NX)+trcsa_TFLS(idsa_H0PO4B,L,NY,NX)+trcsa_RFXS(idsa_H0PO4B,L,NY,NX)+trcsa_RFLZ(idsa_H0PO4B,L,NY,NX)
      trcsa_solml2(idsa_H3PO4B,L,NY,NX)=trcsa_solml2(idsa_H3PO4B,L,NY,NX)+trcsa_TFLS(idsa_H3PO4B,L,NY,NX)+trcsa_RFXS(idsa_H3PO4B,L,NY,NX)+trcsa_RFLZ(idsa_H3PO4B,L,NY,NX)
      trcsa_solml2(idsa_FeHPO4B,L,NY,NX)=trcsa_solml2(idsa_FeHPO4B,L,NY,NX)+trcsa_TFLS(idsa_FeHPO4B,L,NY,NX)+trcsa_RFXS(idsa_FeHPO4B,L,NY,NX)+trcsa_RFLZ(idsa_FeHPO4B,L,NY,NX)
      trcsa_solml2(idsa_FeH2PO4B,L,NY,NX)=trcsa_solml2(idsa_FeH2PO4B,L,NY,NX)+trcsa_TFLS(idsa_FeH2PO4B,L,NY,NX)+trcsa_RFXS(idsa_FeH2PO4B,L,NY,NX)+trcsa_RFLZ(idsa_FeH2PO4B,L,NY,NX)
      trcsa_solml2(idsa_CaPO4B,L,NY,NX)=trcsa_solml2(idsa_CaPO4B,L,NY,NX)+trcsa_TFLS(idsa_CaPO4B,L,NY,NX)+trcsa_RFXS(idsa_CaPO4B,L,NY,NX)+trcsa_RFLZ(idsa_CaPO4B,L,NY,NX)
      trcsa_solml2(idsa_CaHPO4B,L,NY,NX)=trcsa_solml2(idsa_CaHPO4B,L,NY,NX)+trcsa_TFLS(idsa_CaHPO4B,L,NY,NX)+trcsa_RFXS(idsa_CaHPO4B,L,NY,NX)+trcsa_RFLZ(idsa_CaHPO4B,L,NY,NX)
      trcsa_solml2(idsa_CaH2PO4B,L,NY,NX)=trcsa_solml2(idsa_CaH2PO4B,L,NY,NX)+trcsa_TFLS(idsa_CaH2PO4B,L,NY,NX)+trcsa_RFXS(idsa_CaH2PO4B,L,NY,NX)+trcsa_RFLZ(idsa_CaH2PO4B,L,NY,NX)
      trcsa_solml2(idsa_MgHPO4B,L,NY,NX)=trcsa_solml2(idsa_MgHPO4B,L,NY,NX)+trcsa_TFLS(idsa_MgHPO4B,L,NY,NX)+trcsa_RFXS(idsa_MgHPO4B,L,NY,NX)+trcsa_RFLZ(idsa_MgHPO4B,L,NY,NX)

      trcsa_soHml2(idsa_Al,L,NY,NX)=trcsa_soHml2(idsa_Al,L,NY,NX)+TALFHS(L,NY,NX)-trcsa_RFXS(idsa_Al,L,NY,NX)
      trcsa_soHml2(idsa_Fe,L,NY,NX)=trcsa_soHml2(idsa_Fe,L,NY,NX)+TFEFHS(L,NY,NX)-trcsa_RFXS(idsa_Fe,L,NY,NX)
      trcsa_soHml2(idsa_Hp,L,NY,NX)=trcsa_soHml2(idsa_Hp,L,NY,NX)+THYFHS(L,NY,NX)-trcsa_RFXS(idsa_Hp,L,NY,NX)
      trcsa_soHml2(idsa_Ca,L,NY,NX)=trcsa_soHml2(idsa_Ca,L,NY,NX)+TCAFHS(L,NY,NX)-trcsa_RFXS(idsa_Ca,L,NY,NX)
      trcsa_soHml2(idsa_Mg,L,NY,NX)=trcsa_soHml2(idsa_Mg,L,NY,NX)+TMGFHS(L,NY,NX)-trcsa_RFXS(idsa_Mg,L,NY,NX)
      trcsa_soHml2(idsa_Na,L,NY,NX)=trcsa_soHml2(idsa_Na,L,NY,NX)+TNAFHS(L,NY,NX)-trcsa_RFXS(idsa_Na,L,NY,NX)
      trcsa_soHml2(idsa_K,L,NY,NX)=trcsa_soHml2(idsa_K,L,NY,NX)+TKAFHS(L,NY,NX)-trcsa_RFXS(idsa_K,L,NY,NX)
      trcsa_soHml2(idsa_OH,L,NY,NX)=trcsa_soHml2(idsa_OH,L,NY,NX)+TOHFHS(L,NY,NX)-trcsa_RFXS(idsa_OH,L,NY,NX)
      trcsa_soHml2(idsa_SO4,L,NY,NX)=trcsa_soHml2(idsa_SO4,L,NY,NX)+TSOFHS(L,NY,NX)-trcsa_RFXS(idsa_SO4,L,NY,NX)
      trcsa_soHml2(idsa_Cl,L,NY,NX)=trcsa_soHml2(idsa_Cl,L,NY,NX)+TCLFHS(L,NY,NX)-trcsa_RFXS(idsa_Cl,L,NY,NX)
      trcsa_soHml2(idsa_CO3,L,NY,NX)=trcsa_soHml2(idsa_CO3,L,NY,NX)+TC3FHS(L,NY,NX)-trcsa_RFXS(idsa_CO3,L,NY,NX)
      trcsa_soHml2(idsa_HCO3,L,NY,NX)=trcsa_soHml2(idsa_HCO3,L,NY,NX)+THCFHS(L,NY,NX)-trcsa_RFXS(idsa_HCO3,L,NY,NX)
      trcsa_soHml2(idsa_AlOH,L,NY,NX)=trcsa_soHml2(idsa_AlOH,L,NY,NX)+TAL1HS(L,NY,NX)-trcsa_RFXS(idsa_AlOH,L,NY,NX)
      trcsa_soHml2(idsa_AlOH2,L,NY,NX)=trcsa_soHml2(idsa_AlOH2,L,NY,NX)+TAL2HS(L,NY,NX)-trcsa_RFXS(idsa_AlOH2,L,NY,NX)
      trcsa_soHml2(idsa_AlOH3,L,NY,NX)=trcsa_soHml2(idsa_AlOH3,L,NY,NX)+TAL3HS(L,NY,NX)-trcsa_RFXS(idsa_AlOH3,L,NY,NX)
      trcsa_soHml2(idsa_AlOH4,L,NY,NX)=trcsa_soHml2(idsa_AlOH4,L,NY,NX)+TAL4HS(L,NY,NX)-trcsa_RFXS(idsa_AlOH4,L,NY,NX)
      trcsa_soHml2(idsa_AlSO4,L,NY,NX)=trcsa_soHml2(idsa_AlSO4,L,NY,NX)+TALSHS(L,NY,NX)-trcsa_RFXS(idsa_AlSO4,L,NY,NX)
      trcsa_soHml2(idsa_FeOH,L,NY,NX)=trcsa_soHml2(idsa_FeOH,L,NY,NX)+TFE1HS(L,NY,NX)-trcsa_RFXS(idsa_FeOH,L,NY,NX)
      trcsa_soHml2(idsa_FeOH2,L,NY,NX)=trcsa_soHml2(idsa_FeOH2,L,NY,NX)+TFE2HS(L,NY,NX)-trcsa_RFXS(idsa_FeOH2,L,NY,NX)
      trcsa_soHml2(idsa_FeOH3,L,NY,NX)=trcsa_soHml2(idsa_FeOH3,L,NY,NX)+TFE3HS(L,NY,NX)-trcsa_RFXS(idsa_FeOH3,L,NY,NX)
      trcsa_soHml2(idsa_FeOH4,L,NY,NX)=trcsa_soHml2(idsa_FeOH4,L,NY,NX)+TFE4HS(L,NY,NX)-trcsa_RFXS(idsa_FeOH4,L,NY,NX)
      trcsa_soHml2(idsa_FeSO4,L,NY,NX)=trcsa_soHml2(idsa_FeSO4,L,NY,NX)+TFESHS(L,NY,NX)-trcsa_RFXS(idsa_FeSO4,L,NY,NX)
      trcsa_soHml2(idsa_CaOH2,L,NY,NX)=trcsa_soHml2(idsa_CaOH2,L,NY,NX)+TCAOHS(L,NY,NX)-trcsa_RFXS(idsa_CaOH2,L,NY,NX)
      trcsa_soHml2(idsa_CaCO3,L,NY,NX)=trcsa_soHml2(idsa_CaCO3,L,NY,NX)+TCACHS(L,NY,NX)-trcsa_RFXS(idsa_CaCO3,L,NY,NX)
      trcsa_soHml2(idsa_CaHCO3,L,NY,NX)=trcsa_soHml2(idsa_CaHCO3,L,NY,NX)+TCAHHS(L,NY,NX)-trcsa_RFXS(idsa_CaHCO3,L,NY,NX)
      trcsa_soHml2(idsa_CaSO4,L,NY,NX)=trcsa_soHml2(idsa_CaSO4,L,NY,NX)+TCASHS(L,NY,NX)-trcsa_RFXS(idsa_CaSO4,L,NY,NX)
      trcsa_soHml2(idsa_MgOH2,L,NY,NX)=trcsa_soHml2(idsa_MgOH2,L,NY,NX)+TMGOHS(L,NY,NX)-trcsa_RFXS(idsa_MgOH2,L,NY,NX)
      trcsa_soHml2(idsa_MgCO3,L,NY,NX)=trcsa_soHml2(idsa_MgCO3,L,NY,NX)+TMGCHS(L,NY,NX)-trcsa_RFXS(idsa_MgCO3,L,NY,NX)
      trcsa_soHml2(idsa_MgHCO3,L,NY,NX)=trcsa_soHml2(idsa_MgHCO3,L,NY,NX)+TMGHHS(L,NY,NX)-trcsa_RFXS(idsa_MgHCO3,L,NY,NX)
      trcsa_soHml2(idsa_MgSO4,L,NY,NX)=trcsa_soHml2(idsa_MgSO4,L,NY,NX)+TMGSHS(L,NY,NX)-trcsa_RFXS(idsa_MgSO4,L,NY,NX)
      trcsa_soHml2(idsa_NaCO3,L,NY,NX)=trcsa_soHml2(idsa_NaCO3,L,NY,NX)+TNACHS(L,NY,NX)-trcsa_RFXS(idsa_NaCO3,L,NY,NX)
      trcsa_soHml2(idsa_NaSO4,L,NY,NX)=trcsa_soHml2(idsa_NaSO4,L,NY,NX)+TNASHS(L,NY,NX)-trcsa_RFXS(idsa_NaSO4,L,NY,NX)
      trcsa_soHml2(idsa_KSO4,L,NY,NX)=trcsa_soHml2(idsa_KSO4,L,NY,NX)+TKASHS(L,NY,NX)-trcsa_RFXS(idsa_KSO4,L,NY,NX)
      trcsa_soHml2(idsa_H0PO4,L,NY,NX)=trcsa_soHml2(idsa_H0PO4,L,NY,NX)+TH0PHS(L,NY,NX)-trcsa_RFXS(idsa_H0PO4,L,NY,NX)
      trcsa_soHml2(idsa_H3PO4,L,NY,NX)=trcsa_soHml2(idsa_H3PO4,L,NY,NX)+TH3PHS(L,NY,NX)-trcsa_RFXS(idsa_H3PO4,L,NY,NX)
      trcsa_soHml2(idsa_FeHPO4,L,NY,NX)=trcsa_soHml2(idsa_FeHPO4,L,NY,NX)+TF1PHS(L,NY,NX)-trcsa_RFXS(idsa_FeHPO4,L,NY,NX)
      trcsa_soHml2(idsa_FeH2PO4,L,NY,NX)=trcsa_soHml2(idsa_FeH2PO4,L,NY,NX)+TF2PHS(L,NY,NX)-trcsa_RFXS(idsa_FeH2PO4,L,NY,NX)
      trcsa_soHml2(idsa_CaPO4,L,NY,NX)=trcsa_soHml2(idsa_CaPO4,L,NY,NX)+TC0PHS(L,NY,NX)-trcsa_RFXS(idsa_CaPO4,L,NY,NX)
      trcsa_soHml2(idsa_CaHPO4,L,NY,NX)=trcsa_soHml2(idsa_CaHPO4,L,NY,NX)+TC1PHS(L,NY,NX)-trcsa_RFXS(idsa_CaHPO4,L,NY,NX)
      trcsa_soHml2(idsa_CaH2PO4,L,NY,NX)=trcsa_soHml2(idsa_CaH2PO4,L,NY,NX)+TC2PHS(L,NY,NX)-trcsa_RFXS(idsa_CaH2PO4,L,NY,NX)
      trcsa_soHml2(idsa_MgHPO4,L,NY,NX)=trcsa_soHml2(idsa_MgHPO4,L,NY,NX)+TM1PHS(L,NY,NX)-trcsa_RFXS(idsa_MgHPO4,L,NY,NX)
      trcsa_soHml2(idsa_H0PO4B,L,NY,NX)=trcsa_soHml2(idsa_H0PO4B,L,NY,NX)+TH0BHB(L,NY,NX)-trcsa_RFXS(idsa_H0PO4B,L,NY,NX)
      trcsa_soHml2(idsa_H3PO4B,L,NY,NX)=trcsa_soHml2(idsa_H3PO4B,L,NY,NX)+TH3BHB(L,NY,NX)-trcsa_RFXS(idsa_H3PO4B,L,NY,NX)
      trcsa_soHml2(idsa_FeHPO4B,L,NY,NX)=trcsa_soHml2(idsa_FeHPO4B,L,NY,NX)+TF1BHB(L,NY,NX)-trcsa_RFXS(idsa_FeHPO4B,L,NY,NX)
      trcsa_soHml2(idsa_FeH2PO4B,L,NY,NX)=trcsa_soHml2(idsa_FeH2PO4B,L,NY,NX)+TF2BHB(L,NY,NX)-trcsa_RFXS(idsa_FeH2PO4B,L,NY,NX)
      trcsa_soHml2(idsa_CaPO4B,L,NY,NX)=trcsa_soHml2(idsa_CaPO4B,L,NY,NX)+TC0BHB(L,NY,NX)-trcsa_RFXS(idsa_CaPO4B,L,NY,NX)
      trcsa_soHml2(idsa_CaHPO4B,L,NY,NX)=trcsa_soHml2(idsa_CaHPO4B,L,NY,NX)+TC1BHB(L,NY,NX)-trcsa_RFXS(idsa_CaHPO4B,L,NY,NX)
      trcsa_soHml2(idsa_CaH2PO4B,L,NY,NX)=trcsa_soHml2(idsa_CaH2PO4B,L,NY,NX)+TC2BHB(L,NY,NX)-trcsa_RFXS(idsa_CaH2PO4B,L,NY,NX)
      trcsa_soHml2(idsa_MgHPO4B,L,NY,NX)=trcsa_soHml2(idsa_MgHPO4B,L,NY,NX)+TM1BHB(L,NY,NX)-trcsa_RFXS(idsa_MgHPO4B,L,NY,NX)

!     IF(I.EQ.268.AND.L.EQ.1)THEN
!     WRITE(*,444)'ZOH2',I,J,M,NX,NY,L,trcsa_solml2(idsa_OH,L,NY,NX)
!    2,trcsa_TFLS(idsa_OH,L,NY,NX),trcsa_RFXS(idsa_OH,L,NY,NX),trcsa_RFLZ(idsa_OH,L,NY,NX)
!    3,trcsa_RFLS(idsa_OH,3,L-1,NY,NX),trcsa_RFLS(idsa_OH,3,L,NY,NX)
!     WRITE(*,444)'ZAL2',I,J,M,NX,NY,L,trcsa_solml2(idsa_Al,L,NY,NX)
!    2,trcsa_TFLS(idsa_Al,L,NY,NX),trcsa_RFXS(idsa_Al,L,NY,NX),trcsa_RFLZ(idsa_Al,L,NY,NX)
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
    trcsa_TFLS(idsa_beg:idsab_end,L,NY,NX)=0.0
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
