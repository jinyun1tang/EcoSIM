module TrnsfrsMod
!!
! Description:
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
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
  real(r8),allocatable ::  TQRAL(:,:)                         !
  real(r8),allocatable ::  TQRFE(:,:)                         !
  real(r8),allocatable ::  TQRHY(:,:)                         !
  real(r8),allocatable ::  TQRCA(:,:)                         !
  real(r8),allocatable ::  TQRMG(:,:)                         !
  real(r8),allocatable ::  TQRNA(:,:)                         !
  real(r8),allocatable ::  TQRKA(:,:)                         !
  real(r8),allocatable ::  TQROH(:,:)                         !
  real(r8),allocatable ::  TQRSO(:,:)                         !
  real(r8),allocatable ::  TQRCL(:,:)                         !
  real(r8),allocatable ::  TQRC3(:,:)                         !
  real(r8),allocatable ::  TQRHC(:,:)                         !
  real(r8),allocatable ::  TQRAL1(:,:)                        !
  real(r8),allocatable ::  TQRAL2(:,:)                        !
  real(r8),allocatable ::  TQRAL3(:,:)                        !
  real(r8),allocatable ::  TQRAL4(:,:)                        !
  real(r8),allocatable ::  TQRALS(:,:)                        !
  real(r8),allocatable ::  TQRFE1(:,:)                        !
  real(r8),allocatable ::  TQRFE2(:,:)                        !
  real(r8),allocatable ::  TQRFE3(:,:)                        !
  real(r8),allocatable ::  TQRFE4(:,:)                        !
  real(r8),allocatable ::  TQRFES(:,:)                        !
  real(r8),allocatable ::  TQRCAO(:,:)                        !
  real(r8),allocatable ::  TQRCAC(:,:)                        !
  real(r8),allocatable ::  TQRCAH(:,:)                        !
  real(r8),allocatable ::  TQRCAS(:,:)                        !
  real(r8),allocatable ::  TQRMGO(:,:)                        !
  real(r8),allocatable ::  TQRMGC(:,:)                        !
  real(r8),allocatable ::  TQRMGH(:,:)                        !
  real(r8),allocatable ::  TQRMGS(:,:)                        !
  real(r8),allocatable ::  TQRNAC(:,:)                        !
  real(r8),allocatable ::  TQRNAS(:,:)                        !
  real(r8),allocatable ::  TQRKAS(:,:)                        !
  real(r8),allocatable ::  TQRH0P(:,:)                        !
  real(r8),allocatable ::  TQRH3P(:,:)                        !
  real(r8),allocatable ::  TQRF1P(:,:)                        !
  real(r8),allocatable ::  TQRF2P(:,:)                        !
  real(r8),allocatable ::  TQRC0P(:,:)                        !
  real(r8),allocatable ::  TQRC1P(:,:)                        !
  real(r8),allocatable ::  TQRC2P(:,:)                        !
  real(r8),allocatable ::  TQRM1P(:,:)                        !
  real(r8),allocatable ::  TQSAL(:,:)                         !
  real(r8),allocatable ::  TQSFE(:,:)                         !
  real(r8),allocatable ::  TQSHY(:,:)                         !
  real(r8),allocatable ::  TQSCA(:,:)                         !
  real(r8),allocatable ::  TQSMG(:,:)                         !
  real(r8),allocatable ::  TQSNA(:,:)                         !
  real(r8),allocatable ::  TQSKA(:,:)                         !
  real(r8),allocatable ::  TQSOH(:,:)                         !
  real(r8),allocatable ::  TQSSO(:,:)                         !
  real(r8),allocatable ::  TQSCL(:,:)                         !
  real(r8),allocatable ::  TQSC3(:,:)                         !
  real(r8),allocatable ::  TQSHC(:,:)                         !
  real(r8),allocatable ::  TQSAL1(:,:)                        !
  real(r8),allocatable ::  TQSAL2(:,:)                        !
  real(r8),allocatable ::  TQSAL3(:,:)                        !
  real(r8),allocatable ::  TQSAL4(:,:)                        !
  real(r8),allocatable ::  TQSALS(:,:)                        !
  real(r8),allocatable ::  TQSFE1(:,:)                        !
  real(r8),allocatable ::  TQSFE2(:,:)                        !
  real(r8),allocatable ::  TQSFE3(:,:)                        !
  real(r8),allocatable ::  TQSFE4(:,:)                        !
  real(r8),allocatable ::  TQSFES(:,:)                        !
  real(r8),allocatable ::  TQSCAO(:,:)                        !
  real(r8),allocatable ::  TQSCAC(:,:)                        !
  real(r8),allocatable ::  TQSCAH(:,:)                        !
  real(r8),allocatable ::  TQSCAS(:,:)                        !
  real(r8),allocatable ::  TQSMGO(:,:)                        !
  real(r8),allocatable ::  TQSMGC(:,:)                        !
  real(r8),allocatable ::  TQSMGH(:,:)                        !
  real(r8),allocatable ::  TQSMGS(:,:)                        !
  real(r8),allocatable ::  TQSNAC(:,:)                        !
  real(r8),allocatable ::  TQSNAS(:,:)                        !
  real(r8),allocatable ::  TQSKAS(:,:)                        !
  real(r8),allocatable ::  TQSH0P(:,:)                        !
  real(r8),allocatable ::  TQSH3P(:,:)                        !
  real(r8),allocatable ::  TQSF1P(:,:)                        !
  real(r8),allocatable ::  TQSF2P(:,:)                        !
  real(r8),allocatable ::  TQSC0P(:,:)                        !
  real(r8),allocatable ::  TQSC1P(:,:)                        !
  real(r8),allocatable ::  TQSC2P(:,:)                        !
  real(r8),allocatable ::  TQSM1P(:,:)                        !
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
  real(r8),allocatable ::  RALFLZ(:,:,:)                      !
  real(r8),allocatable ::  RFEFLZ(:,:,:)                      !
  real(r8),allocatable ::  RHYFLZ(:,:,:)                      !
  real(r8),allocatable ::  RCAFLZ(:,:,:)                      !
  real(r8),allocatable ::  RMGFLZ(:,:,:)                      !
  real(r8),allocatable ::  RNAFLZ(:,:,:)                      !
  real(r8),allocatable ::  RKAFLZ(:,:,:)                      !
  real(r8),allocatable ::  ROHFLZ(:,:,:)                      !
  real(r8),allocatable ::  RSOFLZ(:,:,:)                      !
  real(r8),allocatable ::  RCLFLZ(:,:,:)                      !
  real(r8),allocatable ::  RC3FLZ(:,:,:)                      !
  real(r8),allocatable ::  RHCFLZ(:,:,:)                      !
  real(r8),allocatable ::  RAL1FZ(:,:,:)                      !
  real(r8),allocatable ::  RAL2FZ(:,:,:)                      !
  real(r8),allocatable ::  RAL3FZ(:,:,:)                      !
  real(r8),allocatable ::  RAL4FZ(:,:,:)                      !
  real(r8),allocatable ::  RALSFZ(:,:,:)                      !
  real(r8),allocatable ::  RCAOFZ(:,:,:)                      !
  real(r8),allocatable ::  RFE1FZ(:,:,:)                      !
  real(r8),allocatable ::  RFE2FZ(:,:,:)                      !
  real(r8),allocatable ::  RFE3FZ(:,:,:)                      !
  real(r8),allocatable ::  RFE4FZ(:,:,:)                      !
  real(r8),allocatable ::  RFESFZ(:,:,:)                      !
  real(r8),allocatable ::  RCACFZ(:,:,:)                      !
  real(r8),allocatable ::  RCAHFZ(:,:,:)                      !
  real(r8),allocatable ::  RCASFZ(:,:,:)                      !
  real(r8),allocatable ::  RMGOFZ(:,:,:)                      !
  real(r8),allocatable ::  RMGCFZ(:,:,:)                      !
  real(r8),allocatable ::  RMGHFZ(:,:,:)                      !
  real(r8),allocatable ::  RMGSFZ(:,:,:)                      !
  real(r8),allocatable ::  RNACFZ(:,:,:)                      !
  real(r8),allocatable ::  RNASFZ(:,:,:)                      !
  real(r8),allocatable ::  RKASFZ(:,:,:)                      !
  real(r8),allocatable ::  RH0PFZ(:,:,:)                      !
  real(r8),allocatable ::  RH3PFZ(:,:,:)                      !
  real(r8),allocatable ::  RF1PFZ(:,:,:)                      !
  real(r8),allocatable ::  RF2PFZ(:,:,:)                      !
  real(r8),allocatable ::  RC0PFZ(:,:,:)                      !
  real(r8),allocatable ::  RC1PFZ(:,:,:)                      !
  real(r8),allocatable ::  RC2PFZ(:,:,:)                      !
  real(r8),allocatable ::  RM1PFZ(:,:,:)                      !
  real(r8),allocatable ::  RH0BBZ(:,:,:)                      !
  real(r8),allocatable ::  RH3BBZ(:,:,:)                      !
  real(r8),allocatable ::  RF1BBZ(:,:,:)                      !
  real(r8),allocatable ::  RF2BBZ(:,:,:)                      !
  real(r8),allocatable ::  RC0BBZ(:,:,:)                      !
  real(r8),allocatable ::  RC1BBZ(:,:,:)                      !
  real(r8),allocatable ::  RC2BBZ(:,:,:)                      !
  real(r8),allocatable ::  RM1BBZ(:,:,:)                      !
!----------------------------------------------------------------------



  public :: trnsfrs
  public :: initTrnsfrs
  public :: destructTrnsfrs
  contains
!----------------------------------------------------------------------
  subroutine DestructTrnsfrs
  use abortutils, only : destroy
  implicit none

  call destroy(TQRAL)
  call destroy(TQRFE)
  call destroy(TQRHY)
  call destroy(TQRCA)
  call destroy(TQRMG)
  call destroy(TQRNA)
  call destroy(TQRKA)
  call destroy(TQROH)
  call destroy(TQRSO)
  call destroy(TQRCL)
  call destroy(TQRC3)
  call destroy(TQRHC)
  call destroy(TQRAL1)
  call destroy(TQRAL2)
  call destroy(TQRAL3)
  call destroy(TQRAL4)
  call destroy(TQRALS)
  call destroy(TQRFE1)
  call destroy(TQRFE2)
  call destroy(TQRFE3)
  call destroy(TQRFE4)
  call destroy(TQRFES)
  call destroy(TQRCAO)
  call destroy(TQRCAC)
  call destroy(TQRCAH)
  call destroy(TQRCAS)
  call destroy(TQRMGO)
  call destroy(TQRMGC)
  call destroy(TQRMGH)
  call destroy(TQRMGS)
  call destroy(TQRNAC)
  call destroy(TQRNAS)
  call destroy(TQRKAS)
  call destroy(TQRH0P)
  call destroy(TQRH3P)
  call destroy(TQRF1P)
  call destroy(TQRF2P)
  call destroy(TQRC0P)
  call destroy(TQRC1P)
  call destroy(TQRC2P)
  call destroy(TQRM1P)
  call destroy(TQSAL)
  call destroy(TQSFE)
  call destroy(TQSHY)
  call destroy(TQSCA)
  call destroy(TQSMG)
  call destroy(TQSNA)
  call destroy(TQSKA)
  call destroy(TQSOH)
  call destroy(TQSSO)
  call destroy(TQSCL)
  call destroy(TQSC3)
  call destroy(TQSHC)
  call destroy(TQSAL1)
  call destroy(TQSAL2)
  call destroy(TQSAL3)
  call destroy(TQSAL4)
  call destroy(TQSALS)
  call destroy(TQSFE1)
  call destroy(TQSFE2)
  call destroy(TQSFE3)
  call destroy(TQSFE4)
  call destroy(TQSFES)
  call destroy(TQSCAO)
  call destroy(TQSCAC)
  call destroy(TQSCAH)
  call destroy(TQSCAS)
  call destroy(TQSMGO)
  call destroy(TQSMGC)
  call destroy(TQSMGH)
  call destroy(TQSMGS)
  call destroy(TQSNAC)
  call destroy(TQSNAS)
  call destroy(TQSKAS)
  call destroy(TQSH0P)
  call destroy(TQSH3P)
  call destroy(TQSF1P)
  call destroy(TQSF2P)
  call destroy(TQSC0P)
  call destroy(TQSC1P)
  call destroy(TQSC2P)
  call destroy(TQSM1P)
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
  call destroy(RALFLZ)
  call destroy(RFEFLZ)
  call destroy(RHYFLZ)
  call destroy(RCAFLZ)
  call destroy(RMGFLZ)
  call destroy(RNAFLZ)
  call destroy(RKAFLZ)
  call destroy(ROHFLZ)
  call destroy(RSOFLZ)
  call destroy(RCLFLZ)
  call destroy(RC3FLZ)
  call destroy(RHCFLZ)
  call destroy(RAL1FZ)
  call destroy(RAL2FZ)
  call destroy(RAL3FZ)
  call destroy(RAL4FZ)
  call destroy(RALSFZ)
  call destroy(RCAOFZ)
  call destroy(RFE1FZ)
  call destroy(RFE2FZ)
  call destroy(RFE3FZ)
  call destroy(RFE4FZ)
  call destroy(RFESFZ)
  call destroy(RCACFZ)
  call destroy(RCAHFZ)
  call destroy(RCASFZ)
  call destroy(RMGOFZ)
  call destroy(RMGCFZ)
  call destroy(RMGHFZ)
  call destroy(RMGSFZ)
  call destroy(RNACFZ)
  call destroy(RNASFZ)
  call destroy(RKASFZ)
  call destroy(RH0PFZ)
  call destroy(RH3PFZ)
  call destroy(RF1PFZ)
  call destroy(RF2PFZ)
  call destroy(RC0PFZ)
  call destroy(RC1PFZ)
  call destroy(RC2PFZ)
  call destroy(RM1PFZ)
  call destroy(RH0BBZ)
  call destroy(RH3BBZ)
  call destroy(RF1BBZ)
  call destroy(RF2BBZ)
  call destroy(RC0BBZ)
  call destroy(RC1BBZ)
  call destroy(RC2BBZ)
  call destroy(RM1BBZ)

  end subroutine DestructTrnsfrs

!----------------------------------------------------------------------

  subroutine initTrnsfrs()
  implicit none

  allocate(TQRAL(JY,JX));       TQRAL=0._r8
  allocate(TQRFE(JY,JX));       TQRFE=0._r8
  allocate(TQRHY(JY,JX));       TQRHY=0._r8
  allocate(TQRCA(JY,JX));       TQRCA=0._r8
  allocate(TQRMG(JY,JX));       TQRMG=0._r8
  allocate(TQRNA(JY,JX));       TQRNA=0._r8
  allocate(TQRKA(JY,JX));       TQRKA=0._r8
  allocate(TQROH(JY,JX));       TQROH=0._r8
  allocate(TQRSO(JY,JX));       TQRSO=0._r8
  allocate(TQRCL(JY,JX));       TQRCL=0._r8
  allocate(TQRC3(JY,JX));       TQRC3=0._r8
  allocate(TQRHC(JY,JX));       TQRHC=0._r8
  allocate(TQRAL1(JY,JX));      TQRAL1=0._r8
  allocate(TQRAL2(JY,JX));      TQRAL2=0._r8
  allocate(TQRAL3(JY,JX));      TQRAL3=0._r8
  allocate(TQRAL4(JY,JX));      TQRAL4=0._r8
  allocate(TQRALS(JY,JX));      TQRALS=0._r8
  allocate(TQRFE1(JY,JX));      TQRFE1=0._r8
  allocate(TQRFE2(JY,JX));      TQRFE2=0._r8
  allocate(TQRFE3(JY,JX));      TQRFE3=0._r8
  allocate(TQRFE4(JY,JX));      TQRFE4=0._r8
  allocate(TQRFES(JY,JX));      TQRFES=0._r8
  allocate(TQRCAO(JY,JX));      TQRCAO=0._r8
  allocate(TQRCAC(JY,JX));      TQRCAC=0._r8
  allocate(TQRCAH(JY,JX));      TQRCAH=0._r8
  allocate(TQRCAS(JY,JX));      TQRCAS=0._r8
  allocate(TQRMGO(JY,JX));      TQRMGO=0._r8
  allocate(TQRMGC(JY,JX));      TQRMGC=0._r8
  allocate(TQRMGH(JY,JX));      TQRMGH=0._r8
  allocate(TQRMGS(JY,JX));      TQRMGS=0._r8
  allocate(TQRNAC(JY,JX));      TQRNAC=0._r8
  allocate(TQRNAS(JY,JX));      TQRNAS=0._r8
  allocate(TQRKAS(JY,JX));      TQRKAS=0._r8
  allocate(TQRH0P(JY,JX));      TQRH0P=0._r8
  allocate(TQRH3P(JY,JX));      TQRH3P=0._r8
  allocate(TQRF1P(JY,JX));      TQRF1P=0._r8
  allocate(TQRF2P(JY,JX));      TQRF2P=0._r8
  allocate(TQRC0P(JY,JX));      TQRC0P=0._r8
  allocate(TQRC1P(JY,JX));      TQRC1P=0._r8
  allocate(TQRC2P(JY,JX));      TQRC2P=0._r8
  allocate(TQRM1P(JY,JX));      TQRM1P=0._r8
  allocate(TQSAL(JY,JX));       TQSAL=0._r8
  allocate(TQSFE(JY,JX));       TQSFE=0._r8
  allocate(TQSHY(JY,JX));       TQSHY=0._r8
  allocate(TQSCA(JY,JX));       TQSCA=0._r8
  allocate(TQSMG(JY,JX));       TQSMG=0._r8
  allocate(TQSNA(JY,JX));       TQSNA=0._r8
  allocate(TQSKA(JY,JX));       TQSKA=0._r8
  allocate(TQSOH(JY,JX));       TQSOH=0._r8
  allocate(TQSSO(JY,JX));       TQSSO=0._r8
  allocate(TQSCL(JY,JX));       TQSCL=0._r8
  allocate(TQSC3(JY,JX));       TQSC3=0._r8
  allocate(TQSHC(JY,JX));       TQSHC=0._r8
  allocate(TQSAL1(JY,JX));      TQSAL1=0._r8
  allocate(TQSAL2(JY,JX));      TQSAL2=0._r8
  allocate(TQSAL3(JY,JX));      TQSAL3=0._r8
  allocate(TQSAL4(JY,JX));      TQSAL4=0._r8
  allocate(TQSALS(JY,JX));      TQSALS=0._r8
  allocate(TQSFE1(JY,JX));      TQSFE1=0._r8
  allocate(TQSFE2(JY,JX));      TQSFE2=0._r8
  allocate(TQSFE3(JY,JX));      TQSFE3=0._r8
  allocate(TQSFE4(JY,JX));      TQSFE4=0._r8
  allocate(TQSFES(JY,JX));      TQSFES=0._r8
  allocate(TQSCAO(JY,JX));      TQSCAO=0._r8
  allocate(TQSCAC(JY,JX));      TQSCAC=0._r8
  allocate(TQSCAH(JY,JX));      TQSCAH=0._r8
  allocate(TQSCAS(JY,JX));      TQSCAS=0._r8
  allocate(TQSMGO(JY,JX));      TQSMGO=0._r8
  allocate(TQSMGC(JY,JX));      TQSMGC=0._r8
  allocate(TQSMGH(JY,JX));      TQSMGH=0._r8
  allocate(TQSMGS(JY,JX));      TQSMGS=0._r8
  allocate(TQSNAC(JY,JX));      TQSNAC=0._r8
  allocate(TQSNAS(JY,JX));      TQSNAS=0._r8
  allocate(TQSKAS(JY,JX));      TQSKAS=0._r8
  allocate(TQSH0P(JY,JX));      TQSH0P=0._r8
  allocate(TQSH3P(JY,JX));      TQSH3P=0._r8
  allocate(TQSF1P(JY,JX));      TQSF1P=0._r8
  allocate(TQSF2P(JY,JX));      TQSF2P=0._r8
  allocate(TQSC0P(JY,JX));      TQSC0P=0._r8
  allocate(TQSC1P(JY,JX));      TQSC1P=0._r8
  allocate(TQSC2P(JY,JX));      TQSC2P=0._r8
  allocate(TQSM1P(JY,JX));      TQSM1P=0._r8
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
  allocate(RALFLZ(JZ,JY,JX));   RALFLZ=0._r8
  allocate(RFEFLZ(JZ,JY,JX));   RFEFLZ=0._r8
  allocate(RHYFLZ(JZ,JY,JX));   RHYFLZ=0._r8
  allocate(RCAFLZ(JZ,JY,JX));   RCAFLZ=0._r8
  allocate(RMGFLZ(JZ,JY,JX));   RMGFLZ=0._r8
  allocate(RNAFLZ(JZ,JY,JX));   RNAFLZ=0._r8
  allocate(RKAFLZ(JZ,JY,JX));   RKAFLZ=0._r8
  allocate(ROHFLZ(JZ,JY,JX));   ROHFLZ=0._r8
  allocate(RSOFLZ(JZ,JY,JX));   RSOFLZ=0._r8
  allocate(RCLFLZ(JZ,JY,JX));   RCLFLZ=0._r8
  allocate(RC3FLZ(JZ,JY,JX));   RC3FLZ=0._r8
  allocate(RHCFLZ(JZ,JY,JX));   RHCFLZ=0._r8
  allocate(RAL1FZ(JZ,JY,JX));   RAL1FZ=0._r8
  allocate(RAL2FZ(JZ,JY,JX));   RAL2FZ=0._r8
  allocate(RAL3FZ(JZ,JY,JX));   RAL3FZ=0._r8
  allocate(RAL4FZ(JZ,JY,JX));   RAL4FZ=0._r8
  allocate(RALSFZ(JZ,JY,JX));   RALSFZ=0._r8
  allocate(RCAOFZ(JZ,JY,JX));   RCAOFZ=0._r8
  allocate(RFE1FZ(JZ,JY,JX));   RFE1FZ=0._r8
  allocate(RFE2FZ(JZ,JY,JX));   RFE2FZ=0._r8
  allocate(RFE3FZ(JZ,JY,JX));   RFE3FZ=0._r8
  allocate(RFE4FZ(JZ,JY,JX));   RFE4FZ=0._r8
  allocate(RFESFZ(JZ,JY,JX));   RFESFZ=0._r8
  allocate(RCACFZ(JZ,JY,JX));   RCACFZ=0._r8
  allocate(RCAHFZ(JZ,JY,JX));   RCAHFZ=0._r8
  allocate(RCASFZ(JZ,JY,JX));   RCASFZ=0._r8
  allocate(RMGOFZ(JZ,JY,JX));   RMGOFZ=0._r8
  allocate(RMGCFZ(JZ,JY,JX));   RMGCFZ=0._r8
  allocate(RMGHFZ(JZ,JY,JX));   RMGHFZ=0._r8
  allocate(RMGSFZ(JZ,JY,JX));   RMGSFZ=0._r8
  allocate(RNACFZ(JZ,JY,JX));   RNACFZ=0._r8
  allocate(RNASFZ(JZ,JY,JX));   RNASFZ=0._r8
  allocate(RKASFZ(JZ,JY,JX));   RKASFZ=0._r8
  allocate(RH0PFZ(JZ,JY,JX));   RH0PFZ=0._r8
  allocate(RH3PFZ(JZ,JY,JX));   RH3PFZ=0._r8
  allocate(RF1PFZ(JZ,JY,JX));   RF1PFZ=0._r8
  allocate(RF2PFZ(JZ,JY,JX));   RF2PFZ=0._r8
  allocate(RC0PFZ(JZ,JY,JX));   RC0PFZ=0._r8
  allocate(RC1PFZ(JZ,JY,JX));   RC1PFZ=0._r8
  allocate(RC2PFZ(JZ,JY,JX));   RC2PFZ=0._r8
  allocate(RM1PFZ(JZ,JY,JX));   RM1PFZ=0._r8
  allocate(RH0BBZ(JZ,JY,JX));   RH0BBZ=0._r8
  allocate(RH3BBZ(JZ,JY,JX));   RH3BBZ=0._r8
  allocate(RF1BBZ(JZ,JY,JX));   RF1BBZ=0._r8
  allocate(RF2BBZ(JZ,JY,JX));   RF2BBZ=0._r8
  allocate(RC0BBZ(JZ,JY,JX));   RC0BBZ=0._r8
  allocate(RC1BBZ(JZ,JY,JX));   RC1BBZ=0._r8
  allocate(RC2BBZ(JZ,JY,JX));   RC2BBZ=0._r8
  allocate(RM1BBZ(JZ,JY,JX));   RM1BBZ=0._r8
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
  DO 30 M=1,NPH
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
      DO 9695 NX=NHW,NHE
        DO 9690 NY=NVN,NVS
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
9690    CONTINUE
9695  CONTINUE
30  CONTINUE
  RETURN

  END subroutine trnsfrs
!------------------------------------------------------------------------------------------

  subroutine ZeroAtmosSoluteFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
  XALBLS(1,NY,NX)=0.0
  XFEBLS(1,NY,NX)=0.0
  XHYBLS(1,NY,NX)=0.0
  XCABLS(1,NY,NX)=0.0
  XMGBLS(1,NY,NX)=0.0
  XNABLS(1,NY,NX)=0.0
  XKABLS(1,NY,NX)=0.0
  XOHBLS(1,NY,NX)=0.0
  XSOBLS(1,NY,NX)=0.0
  XCLBLS(1,NY,NX)=0.0
  XC3BLS(1,NY,NX)=0.0
  XHCBLS(1,NY,NX)=0.0
  XAL1BS(1,NY,NX)=0.0
  XAL2BS(1,NY,NX)=0.0
  XAL3BS(1,NY,NX)=0.0
  XAL4BS(1,NY,NX)=0.0
  XALSBS(1,NY,NX)=0.0
  XFE1BS(1,NY,NX)=0.0
  XFE2BS(1,NY,NX)=0.0
  XFE3BS(1,NY,NX)=0.0
  XFE4BS(1,NY,NX)=0.0
  XFESBS(1,NY,NX)=0.0
  XCAOBS(1,NY,NX)=0.0
  XCACBS(1,NY,NX)=0.0
  XCAHBS(1,NY,NX)=0.0
  XCASBS(1,NY,NX)=0.0
  XMGOBS(1,NY,NX)=0.0
  XMGCBS(1,NY,NX)=0.0
  XMGHBS(1,NY,NX)=0.0
  XMGSBS(1,NY,NX)=0.0
  XNACBS(1,NY,NX)=0.0
  XNASBS(1,NY,NX)=0.0
  XKASBS(1,NY,NX)=0.0
  XH0PBS(1,NY,NX)=0.0
  XH3PBS(1,NY,NX)=0.0
  XF1PBS(1,NY,NX)=0.0
  XF2PBS(1,NY,NX)=0.0
  XC0PBS(1,NY,NX)=0.0
  XC1PBS(1,NY,NX)=0.0
  XC2PBS(1,NY,NX)=0.0
  XM1PBS(1,NY,NX)=0.0
  XALFLS(3,0,NY,NX)=0.0
  XFEFLS(3,0,NY,NX)=0.0
  XHYFLS(3,0,NY,NX)=0.0
  XCAFLS(3,0,NY,NX)=0.0
  XMGFLS(3,0,NY,NX)=0.0
  XNAFLS(3,0,NY,NX)=0.0
  XKAFLS(3,0,NY,NX)=0.0
  XOHFLS(3,0,NY,NX)=0.0
  XSOFLS(3,0,NY,NX)=0.0
  XCLFLS(3,0,NY,NX)=0.0
  XC3FLS(3,0,NY,NX)=0.0
  XHCFLS(3,0,NY,NX)=0.0
  XAL1FS(3,0,NY,NX)=0.0
  XAL2FS(3,0,NY,NX)=0.0
  XAL3FS(3,0,NY,NX)=0.0
  XAL4FS(3,0,NY,NX)=0.0
  XALSFS(3,0,NY,NX)=0.0
  XFE1FS(3,0,NY,NX)=0.0
  XFE2FS(3,0,NY,NX)=0.0
  XFE3FS(3,0,NY,NX)=0.0
  XFE4FS(3,0,NY,NX)=0.0
  XFESFS(3,0,NY,NX)=0.0
  XCAOFS(3,0,NY,NX)=0.0
  XCACFS(3,0,NY,NX)=0.0
  XCAHFS(3,0,NY,NX)=0.0
  XCASFS(3,0,NY,NX)=0.0
  XMGOFS(3,0,NY,NX)=0.0
  XMGCFS(3,0,NY,NX)=0.0
  XMGHFS(3,0,NY,NX)=0.0
  XMGSFS(3,0,NY,NX)=0.0
  XNACFS(3,0,NY,NX)=0.0
  XNASFS(3,0,NY,NX)=0.0
  XKASFS(3,0,NY,NX)=0.0
  XH0PFS(3,0,NY,NX)=0.0
  XH3PFS(3,0,NY,NX)=0.0
  XF1PFS(3,0,NY,NX)=0.0
  XF2PFS(3,0,NY,NX)=0.0
  XC0PFS(3,0,NY,NX)=0.0
  XC1PFS(3,0,NY,NX)=0.0
  XC2PFS(3,0,NY,NX)=0.0
  XM1PFS(3,0,NY,NX)=0.0
  XALFLS(3,NU(NY,NX),NY,NX)=0.0
  XFEFLS(3,NU(NY,NX),NY,NX)=0.0
  XHYFLS(3,NU(NY,NX),NY,NX)=0.0
  XCAFLS(3,NU(NY,NX),NY,NX)=0.0
  XMGFLS(3,NU(NY,NX),NY,NX)=0.0
  XNAFLS(3,NU(NY,NX),NY,NX)=0.0
  XKAFLS(3,NU(NY,NX),NY,NX)=0.0
  XOHFLS(3,NU(NY,NX),NY,NX)=0.0
  XSOFLS(3,NU(NY,NX),NY,NX)=0.0
  XCLFLS(3,NU(NY,NX),NY,NX)=0.0
  XC3FLS(3,NU(NY,NX),NY,NX)=0.0
  XHCFLS(3,NU(NY,NX),NY,NX)=0.0
  XAL1FS(3,NU(NY,NX),NY,NX)=0.0
  XAL2FS(3,NU(NY,NX),NY,NX)=0.0
  XAL3FS(3,NU(NY,NX),NY,NX)=0.0
  XAL4FS(3,NU(NY,NX),NY,NX)=0.0
  XALSFS(3,NU(NY,NX),NY,NX)=0.0
  XFE1FS(3,NU(NY,NX),NY,NX)=0.0
  XFE2FS(3,NU(NY,NX),NY,NX)=0.0
  XFE3FS(3,NU(NY,NX),NY,NX)=0.0
  XFE4FS(3,NU(NY,NX),NY,NX)=0.0
  XFESFS(3,NU(NY,NX),NY,NX)=0.0
  XCAOFS(3,NU(NY,NX),NY,NX)=0.0
  XCACFS(3,NU(NY,NX),NY,NX)=0.0
  XCAHFS(3,NU(NY,NX),NY,NX)=0.0
  XCASFS(3,NU(NY,NX),NY,NX)=0.0
  XMGOFS(3,NU(NY,NX),NY,NX)=0.0
  XMGCFS(3,NU(NY,NX),NY,NX)=0.0
  XMGHFS(3,NU(NY,NX),NY,NX)=0.0
  XMGSFS(3,NU(NY,NX),NY,NX)=0.0
  XNACFS(3,NU(NY,NX),NY,NX)=0.0
  XNASFS(3,NU(NY,NX),NY,NX)=0.0
  XKASFS(3,NU(NY,NX),NY,NX)=0.0
  XH0PFS(3,NU(NY,NX),NY,NX)=0.0
  XH3PFS(3,NU(NY,NX),NY,NX)=0.0
  XF1PFS(3,NU(NY,NX),NY,NX)=0.0
  XF2PFS(3,NU(NY,NX),NY,NX)=0.0
  XC0PFS(3,NU(NY,NX),NY,NX)=0.0
  XC1PFS(3,NU(NY,NX),NY,NX)=0.0
  XC2PFS(3,NU(NY,NX),NY,NX)=0.0
  XM1PFS(3,NU(NY,NX),NY,NX)=0.0
  XH0BFB(3,NU(NY,NX),NY,NX)=0.0
  XH3BFB(3,NU(NY,NX),NY,NX)=0.0
  XF1BFB(3,NU(NY,NX),NY,NX)=0.0
  XF2BFB(3,NU(NY,NX),NY,NX)=0.0
  XC0BFB(3,NU(NY,NX),NY,NX)=0.0
  XC1BFB(3,NU(NY,NX),NY,NX)=0.0
  XC2BFB(3,NU(NY,NX),NY,NX)=0.0
  XM1BFB(3,NU(NY,NX),NY,NX)=0.0
  end subroutine ZeroAtmosSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToSnowpack(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
!     begin_execution

  XALBLS(1,NY,NX)=FLQGQ(NY,NX)*CALR(NY,NX)+FLQGI(NY,NX)*CALQ(I,NY,NX)
  XFEBLS(1,NY,NX)=FLQGQ(NY,NX)*CFER(NY,NX)+FLQGI(NY,NX)*CFEQ(I,NY,NX)
  XHYBLS(1,NY,NX)=FLQGQ(NY,NX)*CHYR(NY,NX)+FLQGI(NY,NX)*CHYQ(I,NY,NX)
  XCABLS(1,NY,NX)=FLQGQ(NY,NX)*CCAR(NY,NX)+FLQGI(NY,NX)*CCAQ(I,NY,NX)
  XMGBLS(1,NY,NX)=FLQGQ(NY,NX)*CMGR(NY,NX)+FLQGI(NY,NX)*CMGQ(I,NY,NX)
  XNABLS(1,NY,NX)=FLQGQ(NY,NX)*CNAR(NY,NX)+FLQGI(NY,NX)*CNAQ(I,NY,NX)
  XKABLS(1,NY,NX)=FLQGQ(NY,NX)*CKAR(NY,NX)+FLQGI(NY,NX)*CKAQ(I,NY,NX)
  XOHBLS(1,NY,NX)=FLQGQ(NY,NX)*COHR(NY,NX)+FLQGI(NY,NX)*COHQ(I,NY,NX)
  XSOBLS(1,NY,NX)=FLQGQ(NY,NX)*CSOR(NY,NX)+FLQGI(NY,NX)*CSOQ(I,NY,NX)
  XCLBLS(1,NY,NX)=FLQGQ(NY,NX)*CCLR(NY,NX)+FLQGI(NY,NX)*CCLQ(I,NY,NX)
  XC3BLS(1,NY,NX)=FLQGQ(NY,NX)*CC3R(NY,NX)+FLQGI(NY,NX)*CC3Q(I,NY,NX)
  XHCBLS(1,NY,NX)=FLQGQ(NY,NX)*CHCR(NY,NX)+FLQGI(NY,NX)*CHCQ(I,NY,NX)
  XAL1BS(1,NY,NX)=FLQGQ(NY,NX)*CAL1R(NY,NX)+FLQGI(NY,NX)*CAL1Q(I,NY,NX)
  XAL2BS(1,NY,NX)=FLQGQ(NY,NX)*CAL2R(NY,NX)+FLQGI(NY,NX)*CAL2Q(I,NY,NX)
  XAL3BS(1,NY,NX)=FLQGQ(NY,NX)*CAL3R(NY,NX)+FLQGI(NY,NX)*CAL3Q(I,NY,NX)
  XAL4BS(1,NY,NX)=FLQGQ(NY,NX)*CAL4R(NY,NX)+FLQGI(NY,NX)*CAL4Q(I,NY,NX)
  XALSBS(1,NY,NX)=FLQGQ(NY,NX)*CALSR(NY,NX)+FLQGI(NY,NX)*CALSQ(I,NY,NX)
  XFE1BS(1,NY,NX)=FLQGQ(NY,NX)*CFE1R(NY,NX)+FLQGI(NY,NX)*CFE1Q(I,NY,NX)
  XFE2BS(1,NY,NX)=FLQGQ(NY,NX)*CFE2R(NY,NX)+FLQGI(NY,NX)*CFE2Q(I,NY,NX)
  XFE3BS(1,NY,NX)=FLQGQ(NY,NX)*CFE3R(NY,NX)+FLQGI(NY,NX)*CFE3Q(I,NY,NX)
  XFE4BS(1,NY,NX)=FLQGQ(NY,NX)*CFE4R(NY,NX)+FLQGI(NY,NX)*CFE4Q(I,NY,NX)
  XFESBS(1,NY,NX)=FLQGQ(NY,NX)*CFESR(NY,NX)+FLQGI(NY,NX)*CFESQ(I,NY,NX)
  XCAOBS(1,NY,NX)=FLQGQ(NY,NX)*CCAOR(NY,NX)+FLQGI(NY,NX)*CCAOQ(I,NY,NX)
  XCACBS(1,NY,NX)=FLQGQ(NY,NX)*CCACR(NY,NX)+FLQGI(NY,NX)*CCACQ(I,NY,NX)
  XCAHBS(1,NY,NX)=FLQGQ(NY,NX)*CCAHR(NY,NX)+FLQGI(NY,NX)*CCAHQ(I,NY,NX)
  XCASBS(1,NY,NX)=FLQGQ(NY,NX)*CCASR(NY,NX)+FLQGI(NY,NX)*CCASQ(I,NY,NX)
  XMGOBS(1,NY,NX)=FLQGQ(NY,NX)*CMGOR(NY,NX)+FLQGI(NY,NX)*CMGOQ(I,NY,NX)
  XMGCBS(1,NY,NX)=FLQGQ(NY,NX)*CMGCR(NY,NX)+FLQGI(NY,NX)*CMGCQ(I,NY,NX)
  XMGHBS(1,NY,NX)=FLQGQ(NY,NX)*CMGHR(NY,NX)+FLQGI(NY,NX)*CMGHQ(I,NY,NX)
  XMGSBS(1,NY,NX)=FLQGQ(NY,NX)*CMGSR(NY,NX)+FLQGI(NY,NX)*CMGSQ(I,NY,NX)
  XNACBS(1,NY,NX)=FLQGQ(NY,NX)*CNACR(NY,NX)+FLQGI(NY,NX)*CNACQ(I,NY,NX)
  XNASBS(1,NY,NX)=FLQGQ(NY,NX)*CNASR(NY,NX)+FLQGI(NY,NX)*CNASQ(I,NY,NX)
  XKASBS(1,NY,NX)=FLQGQ(NY,NX)*CKASR(NY,NX)+FLQGI(NY,NX)*CKASQ(I,NY,NX)
  XH0PBS(1,NY,NX)=FLQGQ(NY,NX)*CH0PR(NY,NX)+FLQGI(NY,NX)*CH0PQ(I,NY,NX)
  XH3PBS(1,NY,NX)=FLQGQ(NY,NX)*CH3PR(NY,NX)+FLQGI(NY,NX)*CH3PQ(I,NY,NX)
  XF1PBS(1,NY,NX)=FLQGQ(NY,NX)*CF1PR(NY,NX)+FLQGI(NY,NX)*CF1PQ(I,NY,NX)
  XF2PBS(1,NY,NX)=FLQGQ(NY,NX)*CF2PR(NY,NX)+FLQGI(NY,NX)*CF2PQ(I,NY,NX)
  XC0PBS(1,NY,NX)=FLQGQ(NY,NX)*CC0PR(NY,NX)+FLQGI(NY,NX)*CC0PQ(I,NY,NX)
  XC1PBS(1,NY,NX)=FLQGQ(NY,NX)*CC1PR(NY,NX)+FLQGI(NY,NX)*CC1PQ(I,NY,NX)
  XC2PBS(1,NY,NX)=FLQGQ(NY,NX)*CC2PR(NY,NX)+FLQGI(NY,NX)*CC2PQ(I,NY,NX)
  XM1PBS(1,NY,NX)=FLQGQ(NY,NX)*CM1PR(NY,NX)+FLQGI(NY,NX)*CM1PQ(I,NY,NX)
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ARE ZERO IF SNOWPACK IS PRESENT
!
!     X*FLS,X*FLB=hourly solute flux to micropores in non-band,band
!
  XALFLS(3,0,NY,NX)=0.0
  XFEFLS(3,0,NY,NX)=0.0
  XHYFLS(3,0,NY,NX)=0.0
  XCAFLS(3,0,NY,NX)=0.0
  XMGFLS(3,0,NY,NX)=0.0
  XNAFLS(3,0,NY,NX)=0.0
  XKAFLS(3,0,NY,NX)=0.0
  XOHFLS(3,0,NY,NX)=0.0
  XSOFLS(3,0,NY,NX)=0.0
  XCLFLS(3,0,NY,NX)=0.0
  XC3FLS(3,0,NY,NX)=0.0
  XHCFLS(3,0,NY,NX)=0.0
  XAL1FS(3,0,NY,NX)=0.0
  XAL2FS(3,0,NY,NX)=0.0
  XAL3FS(3,0,NY,NX)=0.0
  XAL4FS(3,0,NY,NX)=0.0
  XALSFS(3,0,NY,NX)=0.0
  XFE1FS(3,0,NY,NX)=0.0
  XFE2FS(3,0,NY,NX)=0.0
  XFE3FS(3,0,NY,NX)=0.0
  XFE4FS(3,0,NY,NX)=0.0
  XFESFS(3,0,NY,NX)=0.0
  XCAOFS(3,0,NY,NX)=0.0
  XCACFS(3,0,NY,NX)=0.0
  XCAHFS(3,0,NY,NX)=0.0
  XCASFS(3,0,NY,NX)=0.0
  XMGOFS(3,0,NY,NX)=0.0
  XMGCFS(3,0,NY,NX)=0.0
  XMGHFS(3,0,NY,NX)=0.0
  XMGSFS(3,0,NY,NX)=0.0
  XNACFS(3,0,NY,NX)=0.0
  XNASFS(3,0,NY,NX)=0.0
  XKASFS(3,0,NY,NX)=0.0
  XH0PFS(3,0,NY,NX)=0.0
  XH3PFS(3,0,NY,NX)=0.0
  XF1PFS(3,0,NY,NX)=0.0
  XF2PFS(3,0,NY,NX)=0.0
  XC0PFS(3,0,NY,NX)=0.0
  XC1PFS(3,0,NY,NX)=0.0
  XC2PFS(3,0,NY,NX)=0.0
  XM1PFS(3,0,NY,NX)=0.0
  XALFLS(3,NU(NY,NX),NY,NX)=0.0
  XFEFLS(3,NU(NY,NX),NY,NX)=0.0
  XHYFLS(3,NU(NY,NX),NY,NX)=0.0
  XCAFLS(3,NU(NY,NX),NY,NX)=0.0
  XMGFLS(3,NU(NY,NX),NY,NX)=0.0
  XNAFLS(3,NU(NY,NX),NY,NX)=0.0
  XKAFLS(3,NU(NY,NX),NY,NX)=0.0
  XOHFLS(3,NU(NY,NX),NY,NX)=0.0
  XSOFLS(3,NU(NY,NX),NY,NX)=0.0
  XCLFLS(3,NU(NY,NX),NY,NX)=0.0
  XC3FLS(3,NU(NY,NX),NY,NX)=0.0
  XHCFLS(3,NU(NY,NX),NY,NX)=0.0
  XAL1FS(3,NU(NY,NX),NY,NX)=0.0
  XAL2FS(3,NU(NY,NX),NY,NX)=0.0
  XAL3FS(3,NU(NY,NX),NY,NX)=0.0
  XAL4FS(3,NU(NY,NX),NY,NX)=0.0
  XALSFS(3,NU(NY,NX),NY,NX)=0.0
  XFE1FS(3,NU(NY,NX),NY,NX)=0.0
  XFE2FS(3,NU(NY,NX),NY,NX)=0.0
  XFE3FS(3,NU(NY,NX),NY,NX)=0.0
  XFE4FS(3,NU(NY,NX),NY,NX)=0.0
  XFESFS(3,NU(NY,NX),NY,NX)=0.0
  XCAOFS(3,NU(NY,NX),NY,NX)=0.0
  XCACFS(3,NU(NY,NX),NY,NX)=0.0
  XCAHFS(3,NU(NY,NX),NY,NX)=0.0
  XCASFS(3,NU(NY,NX),NY,NX)=0.0
  XMGOFS(3,NU(NY,NX),NY,NX)=0.0
  XMGCFS(3,NU(NY,NX),NY,NX)=0.0
  XMGHFS(3,NU(NY,NX),NY,NX)=0.0
  XMGSFS(3,NU(NY,NX),NY,NX)=0.0
  XNACFS(3,NU(NY,NX),NY,NX)=0.0
  XNASFS(3,NU(NY,NX),NY,NX)=0.0
  XKASFS(3,NU(NY,NX),NY,NX)=0.0
  XH0PFS(3,NU(NY,NX),NY,NX)=0.0
  XH3PFS(3,NU(NY,NX),NY,NX)=0.0
  XF1PFS(3,NU(NY,NX),NY,NX)=0.0
  XF2PFS(3,NU(NY,NX),NY,NX)=0.0
  XC0PFS(3,NU(NY,NX),NY,NX)=0.0
  XC1PFS(3,NU(NY,NX),NY,NX)=0.0
  XC2PFS(3,NU(NY,NX),NY,NX)=0.0
  XM1PFS(3,NU(NY,NX),NY,NX)=0.0
  XH0BFB(3,NU(NY,NX),NY,NX)=0.0
  XH3BFB(3,NU(NY,NX),NY,NX)=0.0
  XF1BFB(3,NU(NY,NX),NY,NX)=0.0
  XF2BFB(3,NU(NY,NX),NY,NX)=0.0
  XC0BFB(3,NU(NY,NX),NY,NX)=0.0
  XC1BFB(3,NU(NY,NX),NY,NX)=0.0
  XC2BFB(3,NU(NY,NX),NY,NX)=0.0
  XM1BFB(3,NU(NY,NX),NY,NX)=0.0
  end subroutine AtmosSoluteFluxToSnowpack
!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToTopsoil(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
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
  XALBLS(1,NY,NX)=0.0
  XFEBLS(1,NY,NX)=0.0
  XHYBLS(1,NY,NX)=0.0
  XCABLS(1,NY,NX)=0.0
  XMGBLS(1,NY,NX)=0.0
  XNABLS(1,NY,NX)=0.0
  XKABLS(1,NY,NX)=0.0
  XOHBLS(1,NY,NX)=0.0
  XSOBLS(1,NY,NX)=0.0
  XCLBLS(1,NY,NX)=0.0
  XC3BLS(1,NY,NX)=0.0
  XHCBLS(1,NY,NX)=0.0
  XAL1BS(1,NY,NX)=0.0
  XAL2BS(1,NY,NX)=0.0
  XAL3BS(1,NY,NX)=0.0
  XAL4BS(1,NY,NX)=0.0
  XALSBS(1,NY,NX)=0.0
  XFE1BS(1,NY,NX)=0.0
  XFE2BS(1,NY,NX)=0.0
  XFE3BS(1,NY,NX)=0.0
  XFE4BS(1,NY,NX)=0.0
  XFESBS(1,NY,NX)=0.0
  XCAOBS(1,NY,NX)=0.0
  XCACBS(1,NY,NX)=0.0
  XCAHBS(1,NY,NX)=0.0
  XCASBS(1,NY,NX)=0.0
  XMGOBS(1,NY,NX)=0.0
  XMGCBS(1,NY,NX)=0.0
  XMGHBS(1,NY,NX)=0.0
  XMGSBS(1,NY,NX)=0.0
  XNACBS(1,NY,NX)=0.0
  XNASBS(1,NY,NX)=0.0
  XKASBS(1,NY,NX)=0.0
  XH0PBS(1,NY,NX)=0.0
  XH3PBS(1,NY,NX)=0.0
  XF1PBS(1,NY,NX)=0.0
  XF2PBS(1,NY,NX)=0.0
  XC0PBS(1,NY,NX)=0.0
  XC1PBS(1,NY,NX)=0.0
  XC2PBS(1,NY,NX)=0.0
  XM1PBS(1,NY,NX)=0.0
  XALFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CALR(NY,NX)+FLQRI(NY,NX)*CALQ(I,NY,NX)
  XFEFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CFER(NY,NX)+FLQRI(NY,NX)*CFEQ(I,NY,NX)
  XHYFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CHYR(NY,NX)+FLQRI(NY,NX)*CHYQ(I,NY,NX)
  XCAFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CCAR(NY,NX)+FLQRI(NY,NX)*CCAQ(I,NY,NX)
  XMGFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CMGR(NY,NX)+FLQRI(NY,NX)*CMGQ(I,NY,NX)
  XNAFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CNAR(NY,NX)+FLQRI(NY,NX)*CNAQ(I,NY,NX)
  XKAFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CKAR(NY,NX)+FLQRI(NY,NX)*CKAQ(I,NY,NX)
  XOHFLS(3,0,NY,NX)=FLQRQ(NY,NX)*COHR(NY,NX)+FLQRI(NY,NX)*COHQ(I,NY,NX)
  XSOFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CSOR(NY,NX)+FLQRI(NY,NX)*CSOQ(I,NY,NX)
  XCLFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CCLR(NY,NX)+FLQRI(NY,NX)*CCLQ(I,NY,NX)
  XC3FLS(3,0,NY,NX)=FLQRQ(NY,NX)*CC3R(NY,NX)+FLQRI(NY,NX)*CC3Q(I,NY,NX)
  XHCFLS(3,0,NY,NX)=FLQRQ(NY,NX)*CHCR(NY,NX)+FLQRI(NY,NX)*CHCQ(I,NY,NX)
  XAL1FS(3,0,NY,NX)=FLQRQ(NY,NX)*CAL1R(NY,NX)+FLQRI(NY,NX)*CAL1Q(I,NY,NX)
  XAL2FS(3,0,NY,NX)=FLQRQ(NY,NX)*CAL2R(NY,NX)+FLQRI(NY,NX)*CAL2Q(I,NY,NX)
  XAL3FS(3,0,NY,NX)=FLQRQ(NY,NX)*CAL3R(NY,NX)+FLQRI(NY,NX)*CAL3Q(I,NY,NX)
  XAL4FS(3,0,NY,NX)=FLQRQ(NY,NX)*CAL4R(NY,NX)+FLQRI(NY,NX)*CAL4Q(I,NY,NX)
  XALSFS(3,0,NY,NX)=FLQRQ(NY,NX)*CALSR(NY,NX)+FLQRI(NY,NX)*CALSQ(I,NY,NX)
  XFE1FS(3,0,NY,NX)=FLQRQ(NY,NX)*CFE1R(NY,NX)+FLQRI(NY,NX)*CFE1Q(I,NY,NX)
  XFE2FS(3,0,NY,NX)=FLQRQ(NY,NX)*CFE2R(NY,NX)+FLQRI(NY,NX)*CFE2Q(I,NY,NX)
  XFE3FS(3,0,NY,NX)=FLQRQ(NY,NX)*CFE3R(NY,NX)+FLQRI(NY,NX)*CFE3Q(I,NY,NX)
  XFE4FS(3,0,NY,NX)=FLQRQ(NY,NX)*CFE4R(NY,NX)+FLQRI(NY,NX)*CFE4Q(I,NY,NX)
  XFESFS(3,0,NY,NX)=FLQRQ(NY,NX)*CFESR(NY,NX)+FLQRI(NY,NX)*CFESQ(I,NY,NX)
  XCAOFS(3,0,NY,NX)=FLQRQ(NY,NX)*CCAOR(NY,NX)+FLQRI(NY,NX)*CCAOQ(I,NY,NX)
  XCACFS(3,0,NY,NX)=FLQRQ(NY,NX)*CCACR(NY,NX)+FLQRI(NY,NX)*CCACQ(I,NY,NX)
  XCAHFS(3,0,NY,NX)=FLQRQ(NY,NX)*CCAHR(NY,NX)+FLQRI(NY,NX)*CCAHQ(I,NY,NX)
  XCASFS(3,0,NY,NX)=FLQRQ(NY,NX)*CCASR(NY,NX)+FLQRI(NY,NX)*CCASQ(I,NY,NX)
  XMGOFS(3,0,NY,NX)=FLQRQ(NY,NX)*CMGOR(NY,NX)+FLQRI(NY,NX)*CMGOQ(I,NY,NX)
  XMGCFS(3,0,NY,NX)=FLQRQ(NY,NX)*CMGCR(NY,NX)+FLQRI(NY,NX)*CMGCQ(I,NY,NX)
  XMGHFS(3,0,NY,NX)=FLQRQ(NY,NX)*CMGHR(NY,NX)+FLQRI(NY,NX)*CMGHQ(I,NY,NX)
  XMGSFS(3,0,NY,NX)=FLQRQ(NY,NX)*CMGSR(NY,NX)+FLQRI(NY,NX)*CMGSQ(I,NY,NX)
  XNACFS(3,0,NY,NX)=FLQRQ(NY,NX)*CNACR(NY,NX)+FLQRI(NY,NX)*CNACQ(I,NY,NX)
  XNASFS(3,0,NY,NX)=FLQRQ(NY,NX)*CNASR(NY,NX)+FLQRI(NY,NX)*CNASQ(I,NY,NX)
  XKASFS(3,0,NY,NX)=FLQRQ(NY,NX)*CKASR(NY,NX)+FLQRI(NY,NX)*CKASQ(I,NY,NX)
  XH0PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CH0PR(NY,NX)+FLQRI(NY,NX)*CH0PQ(I,NY,NX)
  XH3PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CH3PR(NY,NX)+FLQRI(NY,NX)*CH3PQ(I,NY,NX)
  XF1PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CF1PR(NY,NX)+FLQRI(NY,NX)*CF1PQ(I,NY,NX)
  XF2PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CF2PR(NY,NX)+FLQRI(NY,NX)*CF2PQ(I,NY,NX)
  XC0PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CC0PR(NY,NX)+FLQRI(NY,NX)*CC0PQ(I,NY,NX)
  XC1PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CC1PR(NY,NX)+FLQRI(NY,NX)*CC1PQ(I,NY,NX)
  XC2PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CC2PR(NY,NX)+FLQRI(NY,NX)*CC2PQ(I,NY,NX)
  XM1PFS(3,0,NY,NX)=FLQRQ(NY,NX)*CM1PR(NY,NX)+FLQRI(NY,NX)*CM1PQ(I,NY,NX)
  XALFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CALR(NY,NX)+FLQGI(NY,NX)*CALQ(I,NY,NX)
  XFEFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFER(NY,NX)+FLQGI(NY,NX)*CFEQ(I,NY,NX)
  XHYFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CHYR(NY,NX)+FLQGI(NY,NX)*CHYQ(I,NY,NX)
  XCAFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCAR(NY,NX)+FLQGI(NY,NX)*CCAQ(I,NY,NX)
  XMGFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGR(NY,NX)+FLQGI(NY,NX)*CMGQ(I,NY,NX)
  XNAFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNAR(NY,NX)+FLQGI(NY,NX)*CNAQ(I,NY,NX)
  XKAFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CKAR(NY,NX)+FLQGI(NY,NX)*CKAQ(I,NY,NX)
  XOHFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*COHR(NY,NX)+FLQGI(NY,NX)*COHQ(I,NY,NX)
  XSOFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CSOR(NY,NX)+FLQGI(NY,NX)*CSOQ(I,NY,NX)
  XCLFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCLR(NY,NX)+FLQGI(NY,NX)*CCLQ(I,NY,NX)
  XC3FLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CC3R(NY,NX)+FLQGI(NY,NX)*CC3Q(I,NY,NX)
  XHCFLS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CHCR(NY,NX)+FLQGI(NY,NX)*CHCQ(I,NY,NX)
  XAL1FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL1R(NY,NX)+FLQGI(NY,NX)*CAL1Q(I,NY,NX)
  XAL2FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL2R(NY,NX)+FLQGI(NY,NX)*CAL2Q(I,NY,NX)
  XAL3FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL3R(NY,NX)+FLQGI(NY,NX)*CAL3Q(I,NY,NX)
  XAL4FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CAL4R(NY,NX)+FLQGI(NY,NX)*CAL4Q(I,NY,NX)
  XALSFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CALSR(NY,NX)+FLQGI(NY,NX)*CALSQ(I,NY,NX)
  XFE1FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE1R(NY,NX)+FLQGI(NY,NX)*CFE1Q(I,NY,NX)
  XFE2FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE2R(NY,NX)+FLQGI(NY,NX)*CFE2Q(I,NY,NX)
  XFE3FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE3R(NY,NX)+FLQGI(NY,NX)*CFE3Q(I,NY,NX)
  XFE4FS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFE4R(NY,NX)+FLQGI(NY,NX)*CFE4Q(I,NY,NX)
  XFESFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CFESR(NY,NX)+FLQGI(NY,NX)*CFESQ(I,NY,NX)
  XCAOFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCAOR(NY,NX)+FLQGI(NY,NX)*CCAOQ(I,NY,NX)
  XCACFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCACR(NY,NX)+FLQGI(NY,NX)*CCACQ(I,NY,NX)
  XCAHFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCAHR(NY,NX)+FLQGI(NY,NX)*CCAHQ(I,NY,NX)
  XCASFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CCASR(NY,NX)+FLQGI(NY,NX)*CCASQ(I,NY,NX)
  XMGOFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGOR(NY,NX)+FLQGI(NY,NX)*CMGOQ(I,NY,NX)
  XMGCFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGCR(NY,NX)+FLQGI(NY,NX)*CMGCQ(I,NY,NX)
  XMGHFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGHR(NY,NX)+FLQGI(NY,NX)*CMGHQ(I,NY,NX)
  XMGSFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CMGSR(NY,NX)+FLQGI(NY,NX)*CMGSQ(I,NY,NX)
  XNACFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNACR(NY,NX)+FLQGI(NY,NX)*CNACQ(I,NY,NX)
  XNASFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CNASR(NY,NX)+FLQGI(NY,NX)*CNASQ(I,NY,NX)
  XKASFS(3,NU(NY,NX),NY,NX)=FLQGQ(NY,NX)*CKASR(NY,NX)+FLQGI(NY,NX)*CKASQ(I,NY,NX)
  XH0PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH0PR(NY,NX)+FLQGI(NY,NX)*CH0PQ(I,NY,NX)) &
    *VLPO4(NU(NY,NX),NY,NX)
  XH3PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH3PR(NY,NX)+FLQGI(NY,NX) &
    *CH3PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XF1PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF1PR(NY,NX)+FLQGI(NY,NX) &
    *CF1PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XF2PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF2PR(NY,NX)+FLQGI(NY,NX) &
    *CF2PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XC0PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC0PR(NY,NX)+FLQGI(NY,NX) &
    *CC0PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XC1PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC1PR(NY,NX)+FLQGI(NY,NX) &
    *CC1PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XC2PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC2PR(NY,NX)+FLQGI(NY,NX) &
    *CC2PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XM1PFS(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CM1PR(NY,NX)+FLQGI(NY,NX) &
    *CM1PQ(I,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
  XH0BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH0PR(NY,NX)+FLQGI(NY,NX) &
    *CH0PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XH3BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CH3PR(NY,NX)+FLQGI(NY,NX) &
    *CH3PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XF1BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF1PR(NY,NX)+FLQGI(NY,NX) &
    *CF1PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XF2BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CF2PR(NY,NX)+FLQGI(NY,NX) &
    *CF2PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XC0BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC0PR(NY,NX)+FLQGI(NY,NX) &
    *CC0PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XC1BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC1PR(NY,NX)+FLQGI(NY,NX) &
    *CC1PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XC2BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CC2PR(NY,NX)+FLQGI(NY,NX) &
    *CC2PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  XM1BFB(3,NU(NY,NX),NY,NX)=(FLQGQ(NY,NX)*CM1PR(NY,NX)+FLQGI(NY,NX) &
    *CM1PQ(I,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
  end subroutine AtmosSoluteFluxToTopsoil
!------------------------------------------------------------------------------------------

  subroutine InitSolutesInSnowpack(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
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
  DO 20 L=1,JS
    ZALW2(L,NY,NX)=ZALW(L,NY,NX)
    ZFEW2(L,NY,NX)=ZFEW(L,NY,NX)
    ZHYW2(L,NY,NX)=ZHYW(L,NY,NX)
    ZCAW2(L,NY,NX)=ZCAW(L,NY,NX)
    ZMGW2(L,NY,NX)=ZMGW(L,NY,NX)
    ZNAW2(L,NY,NX)=ZNAW(L,NY,NX)
    ZKAW2(L,NY,NX)=ZKAW(L,NY,NX)
    ZOHW2(L,NY,NX)=ZOHW(L,NY,NX)
    ZSO4W2(L,NY,NX)=ZSO4W(L,NY,NX)
    ZCLW2(L,NY,NX)=ZCLW(L,NY,NX)
    ZCO3W2(L,NY,NX)=ZCO3W(L,NY,NX)
    ZHCO3W2(L,NY,NX)=ZHCO3W(L,NY,NX)
    ZALH1W2(L,NY,NX)=ZALH1W(L,NY,NX)
    ZALH2W2(L,NY,NX)=ZALH2W(L,NY,NX)
    ZALH3W2(L,NY,NX)=ZALH3W(L,NY,NX)
    ZALH4W2(L,NY,NX)=ZALH4W(L,NY,NX)
    ZALSW2(L,NY,NX)=ZALSW(L,NY,NX)
    ZFEH1W2(L,NY,NX)=ZFEH1W(L,NY,NX)
    ZFEH2W2(L,NY,NX)=ZFEH2W(L,NY,NX)
    ZFEH3W2(L,NY,NX)=ZFEH3W(L,NY,NX)
    ZFEH4W2(L,NY,NX)=ZFEH4W(L,NY,NX)
    ZFESW2(L,NY,NX)=ZFESW(L,NY,NX)
    ZCAOW2(L,NY,NX)=ZCAOW(L,NY,NX)
    ZCACW2(L,NY,NX)=ZCACW(L,NY,NX)
    ZCAHW2(L,NY,NX)=ZCAHW(L,NY,NX)
    ZCASW2(L,NY,NX)=ZCASW(L,NY,NX)
    ZMGOW2(L,NY,NX)=ZMGOW(L,NY,NX)
    ZMGCW2(L,NY,NX)=ZMGCW(L,NY,NX)
    ZMGHW2(L,NY,NX)=ZMGHW(L,NY,NX)
    ZMGSW2(L,NY,NX)=ZMGSW(L,NY,NX)
    ZNACW2(L,NY,NX)=ZNACW(L,NY,NX)
    ZNASW2(L,NY,NX)=ZNASW(L,NY,NX)
    ZKASW2(L,NY,NX)=ZKASW(L,NY,NX)
    H0PO4W2(L,NY,NX)=H0PO4W(L,NY,NX)
    H3PO4W2(L,NY,NX)=H3PO4W(L,NY,NX)
    ZFE1PW2(L,NY,NX)=ZFE1PW(L,NY,NX)
    ZFE2PW2(L,NY,NX)=ZFE2PW(L,NY,NX)
    ZCA0PW2(L,NY,NX)=ZCA0PW(L,NY,NX)
    ZCA1PW2(L,NY,NX)=ZCA1PW(L,NY,NX)
    ZCA2PW2(L,NY,NX)=ZCA2PW(L,NY,NX)
    ZMG1PW2(L,NY,NX)=ZMG1PW(L,NY,NX)
20  CONTINUE
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
  RALBLS(1,NY,NX)=XALBLS(1,NY,NX)*XNPH
  RFEBLS(1,NY,NX)=XFEBLS(1,NY,NX)*XNPH
  RHYBLS(1,NY,NX)=XHYBLS(1,NY,NX)*XNPH
  RCABLS(1,NY,NX)=XCABLS(1,NY,NX)*XNPH
  RMGBLS(1,NY,NX)=XMGBLS(1,NY,NX)*XNPH
  RNABLS(1,NY,NX)=XNABLS(1,NY,NX)*XNPH
  RKABLS(1,NY,NX)=XKABLS(1,NY,NX)*XNPH
  ROHBLS(1,NY,NX)=XOHBLS(1,NY,NX)*XNPH
  RSOBLS(1,NY,NX)=XSOBLS(1,NY,NX)*XNPH
  RCLBLS(1,NY,NX)=XCLBLS(1,NY,NX)*XNPH
  RC3BLS(1,NY,NX)=XC3BLS(1,NY,NX)*XNPH
  RHCBLS(1,NY,NX)=XHCBLS(1,NY,NX)*XNPH
  RAL1BS(1,NY,NX)=XAL1BS(1,NY,NX)*XNPH
  RAL2BS(1,NY,NX)=XAL2BS(1,NY,NX)*XNPH
  RAL3BS(1,NY,NX)=XAL3BS(1,NY,NX)*XNPH
  RAL4BS(1,NY,NX)=XAL4BS(1,NY,NX)*XNPH
  RALSBS(1,NY,NX)=XALSBS(1,NY,NX)*XNPH
  RFE1BS(1,NY,NX)=XFE1BS(1,NY,NX)*XNPH
  RFE2BS(1,NY,NX)=XFE2BS(1,NY,NX)*XNPH
  RFE3BS(1,NY,NX)=XFE3BS(1,NY,NX)*XNPH
  RFE4BS(1,NY,NX)=XFE4BS(1,NY,NX)*XNPH
  RFESBS(1,NY,NX)=XFESBS(1,NY,NX)*XNPH
  RCAOBS(1,NY,NX)=XCAOBS(1,NY,NX)*XNPH
  RCACBS(1,NY,NX)=XCACBS(1,NY,NX)*XNPH
  RCAHBS(1,NY,NX)=XCAHBS(1,NY,NX)*XNPH
  RCASBS(1,NY,NX)=XCASBS(1,NY,NX)*XNPH
  RMGOBS(1,NY,NX)=XMGOBS(1,NY,NX)*XNPH
  RMGCBS(1,NY,NX)=XMGCBS(1,NY,NX)*XNPH
  RMGHBS(1,NY,NX)=XMGHBS(1,NY,NX)*XNPH
  RMGSBS(1,NY,NX)=XMGSBS(1,NY,NX)*XNPH
  RNACBS(1,NY,NX)=XNACBS(1,NY,NX)*XNPH
  RNASBS(1,NY,NX)=XNASBS(1,NY,NX)*XNPH
  RKASBS(1,NY,NX)=XKASBS(1,NY,NX)*XNPH
  RH0PBS(1,NY,NX)=XH0PBS(1,NY,NX)*XNPH
  RH3PBS(1,NY,NX)=XH3PBS(1,NY,NX)*XNPH
  RF1PBS(1,NY,NX)=XF1PBS(1,NY,NX)*XNPH
  RF2PBS(1,NY,NX)=XF2PBS(1,NY,NX)*XNPH
  RC0PBS(1,NY,NX)=XC0PBS(1,NY,NX)*XNPH
  RC1PBS(1,NY,NX)=XC1PBS(1,NY,NX)*XNPH
  RC2PBS(1,NY,NX)=XC2PBS(1,NY,NX)*XNPH
  RM1PBS(1,NY,NX)=XM1PBS(1,NY,NX)*XNPH
  RALFL0(NY,NX)=XALFLS(3,0,NY,NX)*XNPH
  RFEFL0(NY,NX)=XFEFLS(3,0,NY,NX)*XNPH
  RHYFL0(NY,NX)=XHYFLS(3,0,NY,NX)*XNPH
  RCAFL0(NY,NX)=XCAFLS(3,0,NY,NX)*XNPH
  RMGFL0(NY,NX)=XMGFLS(3,0,NY,NX)*XNPH
  RNAFL0(NY,NX)=XNAFLS(3,0,NY,NX)*XNPH
  RKAFL0(NY,NX)=XKAFLS(3,0,NY,NX)*XNPH
  ROHFL0(NY,NX)=XOHFLS(3,0,NY,NX)*XNPH
  RSOFL0(NY,NX)=XSOFLS(3,0,NY,NX)*XNPH
  RCLFL0(NY,NX)=XCLFLS(3,0,NY,NX)*XNPH
  RC3FL0(NY,NX)=XC3FLS(3,0,NY,NX)*XNPH
  RHCFL0(NY,NX)=XHCFLS(3,0,NY,NX)*XNPH
  RAL1F0(NY,NX)=XAL1FS(3,0,NY,NX)*XNPH
  RAL2F0(NY,NX)=XAL2FS(3,0,NY,NX)*XNPH
  RAL3F0(NY,NX)=XAL3FS(3,0,NY,NX)*XNPH
  RAL4F0(NY,NX)=XAL4FS(3,0,NY,NX)*XNPH
  RALSF0(NY,NX)=XALSFS(3,0,NY,NX)*XNPH
  RFE1F0(NY,NX)=XFE1FS(3,0,NY,NX)*XNPH
  RFE2F0(NY,NX)=XFE2FS(3,0,NY,NX)*XNPH
  RFE3F0(NY,NX)=XFE3FS(3,0,NY,NX)*XNPH
  RFE4F0(NY,NX)=XFE4FS(3,0,NY,NX)*XNPH
  RFESF0(NY,NX)=XFESFS(3,0,NY,NX)*XNPH
  RCAOF0(NY,NX)=XCAOFS(3,0,NY,NX)*XNPH
  RCACF0(NY,NX)=XCACFS(3,0,NY,NX)*XNPH
  RCAHF0(NY,NX)=XCAHFS(3,0,NY,NX)*XNPH
  RCASF0(NY,NX)=XCASFS(3,0,NY,NX)*XNPH
  RMGOF0(NY,NX)=XMGOFS(3,0,NY,NX)*XNPH
  RMGCF0(NY,NX)=XMGCFS(3,0,NY,NX)*XNPH
  RMGHF0(NY,NX)=XMGHFS(3,0,NY,NX)*XNPH
  RMGSF0(NY,NX)=XMGSFS(3,0,NY,NX)*XNPH
  RNACF0(NY,NX)=XNACFS(3,0,NY,NX)*XNPH
  RNASF0(NY,NX)=XNASFS(3,0,NY,NX)*XNPH
  RKASF0(NY,NX)=XKASFS(3,0,NY,NX)*XNPH
  RH0PF0(NY,NX)=XH0PFS(3,0,NY,NX)*XNPH
  RH3PF0(NY,NX)=XH3PFS(3,0,NY,NX)*XNPH
  RF1PF0(NY,NX)=XF1PFS(3,0,NY,NX)*XNPH
  RF2PF0(NY,NX)=XF2PFS(3,0,NY,NX)*XNPH
  RC0PF0(NY,NX)=XC0PFS(3,0,NY,NX)*XNPH
  RC1PF0(NY,NX)=XC1PFS(3,0,NY,NX)*XNPH
  RC2PF0(NY,NX)=XC2PFS(3,0,NY,NX)*XNPH
  RM1PF0(NY,NX)=XM1PFS(3,0,NY,NX)*XNPH
  RALFL1(NY,NX)=XALFLS(3,NU(NY,NX),NY,NX)*XNPH
  RFEFL1(NY,NX)=XFEFLS(3,NU(NY,NX),NY,NX)*XNPH
  RHYFL1(NY,NX)=XHYFLS(3,NU(NY,NX),NY,NX)*XNPH
  RCAFL1(NY,NX)=XCAFLS(3,NU(NY,NX),NY,NX)*XNPH
  RMGFL1(NY,NX)=XMGFLS(3,NU(NY,NX),NY,NX)*XNPH
  RNAFL1(NY,NX)=XNAFLS(3,NU(NY,NX),NY,NX)*XNPH
  RKAFL1(NY,NX)=XKAFLS(3,NU(NY,NX),NY,NX)*XNPH
  ROHFL1(NY,NX)=XOHFLS(3,NU(NY,NX),NY,NX)*XNPH
  RSOFL1(NY,NX)=XSOFLS(3,NU(NY,NX),NY,NX)*XNPH
  RCLFL1(NY,NX)=XCLFLS(3,NU(NY,NX),NY,NX)*XNPH
  RC3FL1(NY,NX)=XC3FLS(3,NU(NY,NX),NY,NX)*XNPH
  RHCFL1(NY,NX)=XHCFLS(3,NU(NY,NX),NY,NX)*XNPH
  RAL1F1(NY,NX)=XAL1FS(3,NU(NY,NX),NY,NX)*XNPH
  RAL2F1(NY,NX)=XAL2FS(3,NU(NY,NX),NY,NX)*XNPH
  RAL3F1(NY,NX)=XAL3FS(3,NU(NY,NX),NY,NX)*XNPH
  RAL4F1(NY,NX)=XAL4FS(3,NU(NY,NX),NY,NX)*XNPH
  RALSF1(NY,NX)=XALSFS(3,NU(NY,NX),NY,NX)*XNPH
  RFE1F1(NY,NX)=XFE1FS(3,NU(NY,NX),NY,NX)*XNPH
  RFE2F1(NY,NX)=XFE2FS(3,NU(NY,NX),NY,NX)*XNPH
  RFE3F1(NY,NX)=XFE3FS(3,NU(NY,NX),NY,NX)*XNPH
  RFE4F1(NY,NX)=XFE4FS(3,NU(NY,NX),NY,NX)*XNPH
  RFESF1(NY,NX)=XFESFS(3,NU(NY,NX),NY,NX)*XNPH
  RCAOF1(NY,NX)=XCAOFS(3,NU(NY,NX),NY,NX)*XNPH
  RCACF1(NY,NX)=XCACFS(3,NU(NY,NX),NY,NX)*XNPH
  RCAHF1(NY,NX)=XCAHFS(3,NU(NY,NX),NY,NX)*XNPH
  RCASF1(NY,NX)=XCASFS(3,NU(NY,NX),NY,NX)*XNPH
  RMGOF1(NY,NX)=XMGOFS(3,NU(NY,NX),NY,NX)*XNPH
  RMGCF1(NY,NX)=XMGCFS(3,NU(NY,NX),NY,NX)*XNPH
  RMGHF1(NY,NX)=XMGHFS(3,NU(NY,NX),NY,NX)*XNPH
  RMGSF1(NY,NX)=XMGSFS(3,NU(NY,NX),NY,NX)*XNPH
  RNACF1(NY,NX)=XNACFS(3,NU(NY,NX),NY,NX)*XNPH
  RNASF1(NY,NX)=XNASFS(3,NU(NY,NX),NY,NX)*XNPH
  RKASF1(NY,NX)=XKASFS(3,NU(NY,NX),NY,NX)*XNPH
  RH0PF1(NY,NX)=XH0PFS(3,NU(NY,NX),NY,NX)*XNPH
  RH3PF1(NY,NX)=XH3PFS(3,NU(NY,NX),NY,NX)*XNPH
  RF1PF1(NY,NX)=XF1PFS(3,NU(NY,NX),NY,NX)*XNPH
  RF2PF1(NY,NX)=XF2PFS(3,NU(NY,NX),NY,NX)*XNPH
  RC0PF1(NY,NX)=XC0PFS(3,NU(NY,NX),NY,NX)*XNPH
  RC1PF1(NY,NX)=XC1PFS(3,NU(NY,NX),NY,NX)*XNPH
  RC2PF1(NY,NX)=XC2PFS(3,NU(NY,NX),NY,NX)*XNPH
  RM1PF1(NY,NX)=XM1PFS(3,NU(NY,NX),NY,NX)*XNPH
  RH0BF2(NY,NX)=XH0BFB(3,NU(NY,NX),NY,NX)*XNPH
  RH3BF2(NY,NX)=XH3BFB(3,NU(NY,NX),NY,NX)*XNPH
  RF1BF2(NY,NX)=XF1BFB(3,NU(NY,NX),NY,NX)*XNPH
  RF2BF2(NY,NX)=XF2BFB(3,NU(NY,NX),NY,NX)*XNPH
  RC0BF2(NY,NX)=XC0BFB(3,NU(NY,NX),NY,NX)*XNPH
  RC1BF2(NY,NX)=XC1BFB(3,NU(NY,NX),NY,NX)*XNPH
  RC2BF2(NY,NX)=XC2BFB(3,NU(NY,NX),NY,NX)*XNPH
  RM1BF2(NY,NX)=XM1BFB(3,NU(NY,NX),NY,NX)*XNPH
  end subroutine GetSubHourFlux
!------------------------------------------------------------------------------------------

  subroutine GetSubHourlyFluxByLayer(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: L
  real(r8) :: FLWU(JZ,JY,JX)


!     begin_execution
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     RZ*2=solute flux at time step for flux calculations
!     TR*=solute transformations from solute.f
!
      DO 10 L=NU(NY,NX),NL(NY,NX)
      RZAL2(L,NY,NX)=-TRAL(L,NY,NX)*XNPH
      RZFE2(L,NY,NX)=-TRFE(L,NY,NX)*XNPH
      RZHY2(L,NY,NX)=-(TRHY(L,NY,NX)+XZHYS(L,NY,NX))*XNPH
      RZCA2(L,NY,NX)=-TRCA(L,NY,NX)*XNPH
      RZMG2(L,NY,NX)=-TRMG(L,NY,NX)*XNPH
      RZNA2(L,NY,NX)=-TRNA(L,NY,NX)*XNPH
      RZKA2(L,NY,NX)=-TRKA(L,NY,NX)*XNPH
      RZOH2(L,NY,NX)=-TROH(L,NY,NX)*XNPH
      RZSO42(L,NY,NX)=-TRSO4(L,NY,NX)*XNPH
      RZCL2(L,NY,NX)=0.0
      RZCO32(L,NY,NX)=-TRCO3(L,NY,NX)*XNPH
      RZHCO32(L,NY,NX)=-TRHCO(L,NY,NX)*XNPH
      RZAL12(L,NY,NX)=-TRAL1(L,NY,NX)*XNPH
      RZAL22(L,NY,NX)=-TRAL2(L,NY,NX)*XNPH
      RZAL32(L,NY,NX)=-TRAL3(L,NY,NX)*XNPH
      RZAL42(L,NY,NX)=-TRAL4(L,NY,NX)*XNPH
      RZALS2(L,NY,NX)=-TRALS(L,NY,NX)*XNPH
      RZFE12(L,NY,NX)=-TRFE1(L,NY,NX)*XNPH
      RZFE22(L,NY,NX)=-TRFE2(L,NY,NX)*XNPH
      RZFE32(L,NY,NX)=-TRFE3(L,NY,NX)*XNPH
      RZFE42(L,NY,NX)=-TRFE4(L,NY,NX)*XNPH
      RZFES2(L,NY,NX)=-TRFES(L,NY,NX)*XNPH
      RZCAO2(L,NY,NX)=-TRCAO(L,NY,NX)*XNPH
      RZCAC2(L,NY,NX)=-TRCAC(L,NY,NX)*XNPH
      RZCAH2(L,NY,NX)=-TRCAH(L,NY,NX)*XNPH
      RZCAS2(L,NY,NX)=-TRCAS(L,NY,NX)*XNPH
      RZMGO2(L,NY,NX)=-TRMGO(L,NY,NX)*XNPH
      RZMGC2(L,NY,NX)=-TRMGC(L,NY,NX)*XNPH
      RZMGH2(L,NY,NX)=-TRMGH(L,NY,NX)*XNPH
      RZMGS2(L,NY,NX)=-TRMGS(L,NY,NX)*XNPH
      RZNAC2(L,NY,NX)=-TRNAC(L,NY,NX)*XNPH
      RZNAS2(L,NY,NX)=-TRNAS(L,NY,NX)*XNPH
      RZKAS2(L,NY,NX)=-TRKAS(L,NY,NX)*XNPH
      RH0PO42(L,NY,NX)=-TRH0P(L,NY,NX)*XNPH
      RH3PO42(L,NY,NX)=-TRH3P(L,NY,NX)*XNPH
      RZFE1P2(L,NY,NX)=-TRF1P(L,NY,NX)*XNPH
      RZFE2P2(L,NY,NX)=-TRF2P(L,NY,NX)*XNPH
      RZCA0P2(L,NY,NX)=-TRC0P(L,NY,NX)*XNPH
      RZCA1P2(L,NY,NX)=-TRC1P(L,NY,NX)*XNPH
      RZCA2P2(L,NY,NX)=-TRC2P(L,NY,NX)*XNPH
      RZMG1P2(L,NY,NX)=-TRM1P(L,NY,NX)*XNPH
      RH0POB2(L,NY,NX)=-TRH0B(L,NY,NX)*XNPH
      RH3POB2(L,NY,NX)=-TRH3B(L,NY,NX)*XNPH
      RZF1PB2(L,NY,NX)=-TRF1B(L,NY,NX)*XNPH
      RZF2PB2(L,NY,NX)=-TRF2B(L,NY,NX)*XNPH
      RZC0PB2(L,NY,NX)=-TRC0B(L,NY,NX)*XNPH
      RZC1PB2(L,NY,NX)=-TRC1B(L,NY,NX)*XNPH
      RZC2PB2(L,NY,NX)=-TRC2B(L,NY,NX)*XNPH
      RZM1PB2(L,NY,NX)=-TRM1B(L,NY,NX)*XNPH
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
      RALFLU(L,NY,NX)=FLU(L,NY,NX)*CALQ(I,NY,NX)
      RFEFLU(L,NY,NX)=FLU(L,NY,NX)*CFEQ(I,NY,NX)
      RHYFLU(L,NY,NX)=FLU(L,NY,NX)*CHYQ(I,NY,NX)
      RCAFLU(L,NY,NX)=FLU(L,NY,NX)*CCAQ(I,NY,NX)
      RMGFLU(L,NY,NX)=FLU(L,NY,NX)*CMGQ(I,NY,NX)
      RNAFLU(L,NY,NX)=FLU(L,NY,NX)*CNAQ(I,NY,NX)
      RKAFLU(L,NY,NX)=FLU(L,NY,NX)*CKAQ(I,NY,NX)
      ROHFLU(L,NY,NX)=FLU(L,NY,NX)*COHQ(I,NY,NX)
      RSOFLU(L,NY,NX)=FLU(L,NY,NX)*CSOQ(I,NY,NX)
      RCLFLU(L,NY,NX)=FLU(L,NY,NX)*CCLQ(I,NY,NX)
      RC3FLU(L,NY,NX)=FLU(L,NY,NX)*CC3Q(I,NY,NX)
      RHCFLU(L,NY,NX)=FLU(L,NY,NX)*CHCQ(I,NY,NX)
      RAL1FU(L,NY,NX)=FLU(L,NY,NX)*CAL1Q(I,NY,NX)
      RAL2FU(L,NY,NX)=FLU(L,NY,NX)*CAL2Q(I,NY,NX)
      RAL3FU(L,NY,NX)=FLU(L,NY,NX)*CAL3Q(I,NY,NX)
      RAL4FU(L,NY,NX)=FLU(L,NY,NX)*CAL4Q(I,NY,NX)
      RALSFU(L,NY,NX)=FLU(L,NY,NX)*CALSQ(I,NY,NX)
      RFE1FU(L,NY,NX)=FLU(L,NY,NX)*CFE1Q(I,NY,NX)
      RFE2FU(L,NY,NX)=FLU(L,NY,NX)*CFE2Q(I,NY,NX)
      RFE3FU(L,NY,NX)=FLU(L,NY,NX)*CFE3Q(I,NY,NX)
      RFE4FU(L,NY,NX)=FLU(L,NY,NX)*CFE4Q(I,NY,NX)
      RFESFU(L,NY,NX)=FLU(L,NY,NX)*CFESQ(I,NY,NX)
      RCAOFU(L,NY,NX)=FLU(L,NY,NX)*CCAOQ(I,NY,NX)
      RCACFU(L,NY,NX)=FLU(L,NY,NX)*CCACQ(I,NY,NX)
      RCAHFU(L,NY,NX)=FLU(L,NY,NX)*CCAHQ(I,NY,NX)
      RCASFU(L,NY,NX)=FLU(L,NY,NX)*CCASQ(I,NY,NX)
      RMGOFU(L,NY,NX)=FLU(L,NY,NX)*CMGOQ(I,NY,NX)
      RMGCFU(L,NY,NX)=FLU(L,NY,NX)*CMGCQ(I,NY,NX)
      RMGHFU(L,NY,NX)=FLU(L,NY,NX)*CMGHQ(I,NY,NX)
      RMGSFU(L,NY,NX)=FLU(L,NY,NX)*CMGSQ(I,NY,NX)
      RNACFU(L,NY,NX)=FLU(L,NY,NX)*CNACQ(I,NY,NX)
      RNASFU(L,NY,NX)=FLU(L,NY,NX)*CNASQ(I,NY,NX)
      RKASFU(L,NY,NX)=FLU(L,NY,NX)*CKASQ(I,NY,NX)
      RH0PFU(L,NY,NX)=FLU(L,NY,NX)*CH0PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RH3PFU(L,NY,NX)=FLU(L,NY,NX)*CH3PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RF1PFU(L,NY,NX)=FLU(L,NY,NX)*CF1PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RF2PFU(L,NY,NX)=FLU(L,NY,NX)*CF2PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RC0PFU(L,NY,NX)=FLU(L,NY,NX)*CC0PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RC1PFU(L,NY,NX)=FLU(L,NY,NX)*CC1PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RC2PFU(L,NY,NX)=FLU(L,NY,NX)*CC2PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RM1PFU(L,NY,NX)=FLU(L,NY,NX)*CM1PQ(I,NY,NX)*VLPO4(L,NY,NX)
      RH0BBU(L,NY,NX)=FLU(L,NY,NX)*CH0PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RH3BBU(L,NY,NX)=FLU(L,NY,NX)*CH3PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RF1BBU(L,NY,NX)=FLU(L,NY,NX)*CF1PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RF2BBU(L,NY,NX)=FLU(L,NY,NX)*CF2PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RC0BBU(L,NY,NX)=FLU(L,NY,NX)*CC0PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RC1BBU(L,NY,NX)=FLU(L,NY,NX)*CC1PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RC2BBU(L,NY,NX)=FLU(L,NY,NX)*CC2PQ(I,NY,NX)*VLPOB(L,NY,NX)
      RM1BBU(L,NY,NX)=FLU(L,NY,NX)*CM1PQ(I,NY,NX)*VLPOB(L,NY,NX)
!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!
      RALFLZ(L,NY,NX)=RALFLU(L,NY,NX)*XNPH
      RFEFLZ(L,NY,NX)=RFEFLU(L,NY,NX)*XNPH
      RHYFLZ(L,NY,NX)=RHYFLU(L,NY,NX)*XNPH
      RCAFLZ(L,NY,NX)=RCAFLU(L,NY,NX)*XNPH
      RMGFLZ(L,NY,NX)=RMGFLU(L,NY,NX)*XNPH
      RNAFLZ(L,NY,NX)=RNAFLU(L,NY,NX)*XNPH
      RKAFLZ(L,NY,NX)=RKAFLU(L,NY,NX)*XNPH
      ROHFLZ(L,NY,NX)=ROHFLU(L,NY,NX)*XNPH
      RSOFLZ(L,NY,NX)=RSOFLU(L,NY,NX)*XNPH
      RCLFLZ(L,NY,NX)=RCLFLU(L,NY,NX)*XNPH
      RC3FLZ(L,NY,NX)=RC3FLU(L,NY,NX)*XNPH
      RHCFLZ(L,NY,NX)=RHCFLU(L,NY,NX)*XNPH
      RAL1FZ(L,NY,NX)=RAL1FU(L,NY,NX)*XNPH
      RAL2FZ(L,NY,NX)=RAL2FU(L,NY,NX)*XNPH
      RAL3FZ(L,NY,NX)=RAL3FU(L,NY,NX)*XNPH
      RAL4FZ(L,NY,NX)=RAL4FU(L,NY,NX)*XNPH
      RALSFZ(L,NY,NX)=RALSFU(L,NY,NX)*XNPH
      RFE1FZ(L,NY,NX)=RFE1FU(L,NY,NX)*XNPH
      RFE2FZ(L,NY,NX)=RFE2FU(L,NY,NX)*XNPH
      RFE3FZ(L,NY,NX)=RFE3FU(L,NY,NX)*XNPH
      RFE4FZ(L,NY,NX)=RFE4FU(L,NY,NX)*XNPH
      RFESFZ(L,NY,NX)=RFESFU(L,NY,NX)*XNPH
      RCAOFZ(L,NY,NX)=RCAOFU(L,NY,NX)*XNPH
      RCACFZ(L,NY,NX)=RCACFU(L,NY,NX)*XNPH
      RCAHFZ(L,NY,NX)=RCAHFU(L,NY,NX)*XNPH
      RCASFZ(L,NY,NX)=RCASFU(L,NY,NX)*XNPH
      RMGOFZ(L,NY,NX)=RMGOFU(L,NY,NX)*XNPH
      RMGCFZ(L,NY,NX)=RMGCFU(L,NY,NX)*XNPH
      RMGHFZ(L,NY,NX)=RMGHFU(L,NY,NX)*XNPH
      RMGSFZ(L,NY,NX)=RMGSFU(L,NY,NX)*XNPH
      RNACFZ(L,NY,NX)=RNACFU(L,NY,NX)*XNPH
      RNASFZ(L,NY,NX)=RNASFU(L,NY,NX)*XNPH
      RKASFZ(L,NY,NX)=RKASFU(L,NY,NX)*XNPH
      RH0PFZ(L,NY,NX)=RH0PFU(L,NY,NX)*XNPH
      RH3PFZ(L,NY,NX)=RH3PFU(L,NY,NX)*XNPH
      RF1PFZ(L,NY,NX)=RF1PFU(L,NY,NX)*XNPH
      RF2PFZ(L,NY,NX)=RF2PFU(L,NY,NX)*XNPH
      RC0PFZ(L,NY,NX)=RC0PFU(L,NY,NX)*XNPH
      RC1PFZ(L,NY,NX)=RC1PFU(L,NY,NX)*XNPH
      RC2PFZ(L,NY,NX)=RC2PFU(L,NY,NX)*XNPH
      RM1PFZ(L,NY,NX)=RM1PFU(L,NY,NX)*XNPH
      RH0BBZ(L,NY,NX)=RH0BBU(L,NY,NX)*XNPH
      RH3BBZ(L,NY,NX)=RH3BBU(L,NY,NX)*XNPH
      RF1BBZ(L,NY,NX)=RF1BBU(L,NY,NX)*XNPH
      RF2BBZ(L,NY,NX)=RF2BBU(L,NY,NX)*XNPH
      RC0BBZ(L,NY,NX)=RC0BBU(L,NY,NX)*XNPH
      RC1BBZ(L,NY,NX)=RC1BBU(L,NY,NX)*XNPH
      RC2BBZ(L,NY,NX)=RC2BBU(L,NY,NX)*XNPH
      RM1BBZ(L,NY,NX)=RM1BBU(L,NY,NX)*XNPH
!
!     SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
!
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:PO=PO4,AL=Al,FE=Fe,HY=H,CA=Ca,GM=Mg,AN=Na,AK=KOH=OH
!                :SO=SO4,CL=Cl,C3=CO3,HC=HCO3
      POSGL2(L,NY,NX)=POSGL(L,NY,NX)*XNPH
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
      ZAL2(L,NY,NX)=ZAL(L,NY,NX)
      ZFE2(L,NY,NX)=ZFE(L,NY,NX)
      ZHY2(L,NY,NX)=ZHY(L,NY,NX)
      ZCA2(L,NY,NX)=ZCA(L,NY,NX)
      ZMG2(L,NY,NX)=ZMG(L,NY,NX)
      ZNA2(L,NY,NX)=ZNA(L,NY,NX)
      ZKA2(L,NY,NX)=ZKA(L,NY,NX)
      ZOH2(L,NY,NX)=ZOH(L,NY,NX)
      ZSO42(L,NY,NX)=ZSO4(L,NY,NX)
      ZCL2(L,NY,NX)=ZCL(L,NY,NX)
      ZCO32(L,NY,NX)=ZCO3(L,NY,NX)
      ZHCO32(L,NY,NX)=ZHCO3(L,NY,NX)
      ZAL12(L,NY,NX)=ZALOH1(L,NY,NX)
      ZAL22(L,NY,NX)=ZALOH2(L,NY,NX)
      ZAL32(L,NY,NX)=ZALOH3(L,NY,NX)
      ZAL42(L,NY,NX)=ZALOH4(L,NY,NX)
      ZALS2(L,NY,NX)=ZALS(L,NY,NX)
      ZFE12(L,NY,NX)=ZFEOH1(L,NY,NX)
      ZFE22(L,NY,NX)=ZFEOH2(L,NY,NX)
      ZFE32(L,NY,NX)=ZFEOH3(L,NY,NX)
      ZFE42(L,NY,NX)=ZFEOH4(L,NY,NX)
      ZFES2(L,NY,NX)=ZFES(L,NY,NX)
      ZCAO2(L,NY,NX)=ZCAO(L,NY,NX)
      ZCAC2(L,NY,NX)=ZCAC(L,NY,NX)
      ZCAH2(L,NY,NX)=ZCAH(L,NY,NX)
      ZCAS2(L,NY,NX)=ZCAS(L,NY,NX)
      ZMGO2(L,NY,NX)=ZMGO(L,NY,NX)
      ZMGC2(L,NY,NX)=ZMGC(L,NY,NX)
      ZMGH2(L,NY,NX)=ZMGH(L,NY,NX)
      ZMGS2(L,NY,NX)=ZMGS(L,NY,NX)
      ZNAC2(L,NY,NX)=ZNAC(L,NY,NX)
      ZNAS2(L,NY,NX)=ZNAS(L,NY,NX)
      ZKAS2(L,NY,NX)=ZKAS(L,NY,NX)
      H0PO42(L,NY,NX)=H0PO4(L,NY,NX)
      H3PO42(L,NY,NX)=H3PO4(L,NY,NX)
      ZFE1P2(L,NY,NX)=ZFE1P(L,NY,NX)
      ZFE2P2(L,NY,NX)=ZFE2P(L,NY,NX)
      ZCA0P2(L,NY,NX)=ZCA0P(L,NY,NX)
      ZCA1P2(L,NY,NX)=ZCA1P(L,NY,NX)
      ZCA2P2(L,NY,NX)=ZCA2P(L,NY,NX)
      ZMG1P2(L,NY,NX)=ZMG1P(L,NY,NX)
      H0POB2(L,NY,NX)=H0POB(L,NY,NX)
      H3POB2(L,NY,NX)=H3POB(L,NY,NX)
      ZF1PB2(L,NY,NX)=ZFE1PB(L,NY,NX)
      ZF2PB2(L,NY,NX)=ZFE2PB(L,NY,NX)
      ZC0PB2(L,NY,NX)=ZCA0PB(L,NY,NX)
      ZC1PB2(L,NY,NX)=ZCA1PB(L,NY,NX)
      ZC2PB2(L,NY,NX)=ZCA2PB(L,NY,NX)
      ZM1PB2(L,NY,NX)=ZMG1PB(L,NY,NX)
      ZALH2(L,NY,NX)=ZALH(L,NY,NX)
      ZFEH2(L,NY,NX)=ZFEH(L,NY,NX)
      ZHYH2(L,NY,NX)=ZHYH(L,NY,NX)
      ZCCH2(L,NY,NX)=ZCCH(L,NY,NX)
      ZMAH2(L,NY,NX)=ZMAH(L,NY,NX)
      ZNAH2(L,NY,NX)=ZNAH(L,NY,NX)
      ZKAH2(L,NY,NX)=ZKAH(L,NY,NX)
      ZOHH2(L,NY,NX)=ZOHH(L,NY,NX)
      ZSO4H2(L,NY,NX)=ZSO4H(L,NY,NX)
      ZCLH2(L,NY,NX)=ZCLH(L,NY,NX)
      ZCO3H2(L,NY,NX)=ZCO3H(L,NY,NX)
      ZHCOH2(L,NY,NX)=ZHCO3H(L,NY,NX)
      ZAL1H2(L,NY,NX)=ZALO1H(L,NY,NX)
      ZAL2H2(L,NY,NX)=ZALO2H(L,NY,NX)
      ZAL3H2(L,NY,NX)=ZALO3H(L,NY,NX)
      ZAL4H2(L,NY,NX)=ZALO4H(L,NY,NX)
      ZALSH2(L,NY,NX)=ZALSH(L,NY,NX)
      ZFE1H2(L,NY,NX)=ZFEO1H(L,NY,NX)
      ZFE2H2(L,NY,NX)=ZFEO2H(L,NY,NX)
      ZFE3H2(L,NY,NX)=ZFEO3H(L,NY,NX)
      ZFE4H2(L,NY,NX)=ZFEO4H(L,NY,NX)
      ZFESH2(L,NY,NX)=ZFESH(L,NY,NX)
      ZCAOH2(L,NY,NX)=ZCAOH(L,NY,NX)
      ZCACH2(L,NY,NX)=ZCACH(L,NY,NX)
      ZCAHH2(L,NY,NX)=ZCAHH(L,NY,NX)
      ZCASH2(L,NY,NX)=ZCASH(L,NY,NX)
      ZMGOH2(L,NY,NX)=ZMGOH(L,NY,NX)
      ZMGCH2(L,NY,NX)=ZMGCH(L,NY,NX)
      ZMGHH2(L,NY,NX)=ZMGHH(L,NY,NX)
      ZMGSH2(L,NY,NX)=ZMGSH(L,NY,NX)
      ZNACH2(L,NY,NX)=ZNACH(L,NY,NX)
      ZNASH2(L,NY,NX)=ZNASH(L,NY,NX)
      ZKASH2(L,NY,NX)=ZKASH(L,NY,NX)
      H0P4H2(L,NY,NX)=H0PO4H(L,NY,NX)
      H3P4H2(L,NY,NX)=H3PO4H(L,NY,NX)
      ZF1PH2(L,NY,NX)=ZFE1PH(L,NY,NX)
      ZF2PH2(L,NY,NX)=ZFE2PH(L,NY,NX)
      ZC0PH2(L,NY,NX)=ZCA0PH(L,NY,NX)
      ZC1PH2(L,NY,NX)=ZCA1PH(L,NY,NX)
      ZC2PH2(L,NY,NX)=ZCA2PH(L,NY,NX)
      ZM1PH2(L,NY,NX)=ZMG1PH(L,NY,NX)
      H0PBH2(L,NY,NX)=H0POBH(L,NY,NX)
      H3PBH2(L,NY,NX)=H3POBH(L,NY,NX)
      ZF1BH2(L,NY,NX)=ZFE1BH(L,NY,NX)
      ZF2BH2(L,NY,NX)=ZFE2BH(L,NY,NX)
      ZC0BH2(L,NY,NX)=ZCA0BH(L,NY,NX)
      ZC1BH2(L,NY,NX)=ZCA1BH(L,NY,NX)
      ZC2BH2(L,NY,NX)=ZCA2BH(L,NY,NX)
      ZM1BH2(L,NY,NX)=ZMG1BH(L,NY,NX)
10    CONTINUE
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

      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
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

9990  CONTINUE
9995  CONTINUE
      end subroutine SaltModelSoluteFlux
!------------------------------------------------------------------------------------------

      subroutine InitFluxAccumulatorsInSnowpack(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
      integer :: L

  DO 9855 L=1,JS
      TALBLS(L,NY,NX)=0.0
      TFEBLS(L,NY,NX)=0.0
      THYBLS(L,NY,NX)=0.0
      TCABLS(L,NY,NX)=0.0
      TMGBLS(L,NY,NX)=0.0
      TNABLS(L,NY,NX)=0.0
      TKABLS(L,NY,NX)=0.0
      TOHBLS(L,NY,NX)=0.0
      TSOBLS(L,NY,NX)=0.0
      TCLBLS(L,NY,NX)=0.0
      TC3BLS(L,NY,NX)=0.0
      THCBLS(L,NY,NX)=0.0
      TAL1BS(L,NY,NX)=0.0
      TAL2BS(L,NY,NX)=0.0
      TAL3BS(L,NY,NX)=0.0
      TAL4BS(L,NY,NX)=0.0
      TALSBS(L,NY,NX)=0.0
      TFE1BS(L,NY,NX)=0.0
      TFE2BS(L,NY,NX)=0.0
      TFE3BS(L,NY,NX)=0.0
      TFE4BS(L,NY,NX)=0.0
      TFESBS(L,NY,NX)=0.0
      TCAOBS(L,NY,NX)=0.0
      TCACBS(L,NY,NX)=0.0
      TCAHBS(L,NY,NX)=0.0
      TCASBS(L,NY,NX)=0.0
      TMGOBS(L,NY,NX)=0.0
      TMGCBS(L,NY,NX)=0.0
      TMGHBS(L,NY,NX)=0.0
      TMGSBS(L,NY,NX)=0.0
      TNACBS(L,NY,NX)=0.0
      TNASBS(L,NY,NX)=0.0
      TKASBS(L,NY,NX)=0.0
      TH0PBS(L,NY,NX)=0.0
      TH3PBS(L,NY,NX)=0.0
      TF1PBS(L,NY,NX)=0.0
      TF2PBS(L,NY,NX)=0.0
      TC0PBS(L,NY,NX)=0.0
      TC1PBS(L,NY,NX)=0.0
      TC2PBS(L,NY,NX)=0.0
      TM1PBS(L,NY,NX)=0.0
9855  CONTINUE
  end subroutine InitFluxAccumulatorsInSnowpack

!------------------------------------------------------------------------------------------

      subroutine BoundarySnowFlux(N,M5,M4)
!
!     Description:
!
      implicit none
      integer, intent(in) :: N,M5,M4
!     begin_execution
      RQSAL(N,M5,M4)=0.0
      RQSFE(N,M5,M4)=0.0
      RQSHY(N,M5,M4)=0.0
      RQSCA(N,M5,M4)=0.0
      RQSMG(N,M5,M4)=0.0
      RQSNA(N,M5,M4)=0.0
      RQSKA(N,M5,M4)=0.0
      RQSOH(N,M5,M4)=0.0
      RQSSO(N,M5,M4)=0.0
      RQSCL(N,M5,M4)=0.0
      RQSC3(N,M5,M4)=0.0
      RQSHC(N,M5,M4)=0.0
      RQSAL1(N,M5,M4)=0.0
      RQSAL2(N,M5,M4)=0.0
      RQSAL3(N,M5,M4)=0.0
      RQSAL4(N,M5,M4)=0.0
      RQSALS(N,M5,M4)=0.0
      RQSFE1(N,M5,M4)=0.0
      RQSFE2(N,M5,M4)=0.0
      RQSFE3(N,M5,M4)=0.0
      RQSFE4(N,M5,M4)=0.0
      RQSFES(N,M5,M4)=0.0
      RQSCAO(N,M5,M4)=0.0
      RQSCAC(N,M5,M4)=0.0
      RQSCAH(N,M5,M4)=0.0
      RQSCAS(N,M5,M4)=0.0
      RQSMGO(N,M5,M4)=0.0
      RQSMGC(N,M5,M4)=0.0
      RQSMGH(N,M5,M4)=0.0
      RQSMGS(N,M5,M4)=0.0
      RQSNAC(N,M5,M4)=0.0
      RQSNAS(N,M5,M4)=0.0
      RQSKAS(N,M5,M4)=0.0
      RQSH0P(N,M5,M4)=0.0
      RQSH3P(N,M5,M4)=0.0
      RQSF1P(N,M5,M4)=0.0
      RQSF2P(N,M5,M4)=0.0
      RQSC0P(N,M5,M4)=0.0
      RQSC1P(N,M5,M4)=0.0
      RQSC2P(N,M5,M4)=0.0
      RQSM1P(N,M5,M4)=0.0
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
      RALFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CALU(M3,M2,M1)
      RFEFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFEU(M3,M2,M1)
      RHYFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CHYU(M3,M2,M1)
      RCAFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCAU(M3,M2,M1)
      RMGFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGU(M3,M2,M1)
      RNAFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNAU(M3,M2,M1)
      RKAFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CKAU(M3,M2,M1)
      ROHFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*COHU(M3,M2,M1)
      RSOFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CSOU(M3,M2,M1)
      RCLFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCLU(M3,M2,M1)
      RC3FLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC3U(M3,M2,M1)
      RHCFLS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CHCU(M3,M2,M1)
      RAL1FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL1U(M3,M2,M1)
      RAL2FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL2U(M3,M2,M1)
      RAL3FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL3U(M3,M2,M1)
      RAL4FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CAL4U(M3,M2,M1)
      RALSFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CALSU(M3,M2,M1)
      RFE1FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE1U(M3,M2,M1)
      RFE2FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE2U(M3,M2,M1)
      RFE3FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE3U(M3,M2,M1)
      RFE4FS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFE4U(M3,M2,M1)
      RFESFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CFESU(M3,M2,M1)
      RCAOFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCAOU(M3,M2,M1)
      RCACFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCACU(M3,M2,M1)
      RCAHFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCAHU(M3,M2,M1)
      RCASFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CCASU(M3,M2,M1)
      RMGOFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGOU(M3,M2,M1)
      RMGCFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGCU(M3,M2,M1)
      RMGHFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGHU(M3,M2,M1)
      RMGSFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CMGSU(M3,M2,M1)
      RNACFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNACU(M3,M2,M1)
      RNASFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNASU(M3,M2,M1)
      RKASFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CKASU(M3,M2,M1)
      RH0PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH0PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RH3PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH3PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RF1PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF1PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RF2PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF2PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RC0PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC0PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RC1PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC1PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RC2PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC2PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RM1PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CM1PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RH0BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH0PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RH3BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH3PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RF1BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF1PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RF2BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CF2PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RC0BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC0PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RC1BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC1PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RC2BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CC2PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RM1BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CM1PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
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
      RALFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL2(M3,M2,M1))
      RFEFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE2(M3,M2,M1))
      RHYFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZHY2(M3,M2,M1))
      RCAFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCA2(M3,M2,M1))
      RMGFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMG2(M3,M2,M1))
      RNAFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNA2(M3,M2,M1))
      RKAFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZKA2(M3,M2,M1))
      ROHFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZOH2(M3,M2,M1))
      RSOFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZSO42(M3,M2,M1))
      RCLFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCL2(M3,M2,M1))
      RC3FLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCO32(M3,M2,M1))
      RHCFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZHCO32(M3,M2,M1))
      RAL1FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL12(M3,M2,M1))
      RAL2FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL22(M3,M2,M1))
      RAL3FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL32(M3,M2,M1))
      RAL4FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL42(M3,M2,M1))
      RALSFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZALS2(M3,M2,M1))
      RFE1FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE12(M3,M2,M1))
      RFE2FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE22(M3,M2,M1))
      RFE3FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE32(M3,M2,M1))
      RFE4FS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE42(M3,M2,M1))
      RFESFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFES2(M3,M2,M1))
      RCAOFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCAO2(M3,M2,M1))
      RCACFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCAC2(M3,M2,M1))
      RCAHFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCAH2(M3,M2,M1))
      RCASFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCAS2(M3,M2,M1))
      RMGOFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGO2(M3,M2,M1))
      RMGCFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGC2(M3,M2,M1))
      RMGHFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGH2(M3,M2,M1))
      RMGSFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGS2(M3,M2,M1))
      RNACFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNAC2(M3,M2,M1))
      RNASFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNAS2(M3,M2,M1))
      RKASFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZKAS2(M3,M2,M1))
      RH0PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H0PO42(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RH3PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H3PO42(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RF1PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE1P2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RF2PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE2P2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RC0PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCA0P2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RC1PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCA1P2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RC2PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCA2P2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RM1PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMG1P2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RH0BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H0POB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RH3BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H3POB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RF1BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZF1PB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RF2BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZF2PB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RC0BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC0PB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RC1BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC1PB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RC2BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC2PB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RM1BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZM1PB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
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
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWHM(M,N,M6,M5,M4) &
      /VOLWHM(M,M3,M2,M1)))
      ELSE
      VFLW=0.0
      ENDIF
      RALFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZALH2(M3,M2,M1))
      RFEFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFEH2(M3,M2,M1))
      RHYFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZHYH2(M3,M2,M1))
      RCAFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCCH2(M3,M2,M1))
      RMGFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMAH2(M3,M2,M1))
      RNAFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNAH2(M3,M2,M1))
      RKAFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZKAH2(M3,M2,M1))
      ROHFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZOHH2(M3,M2,M1))
      RSOFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZSO4H2(M3,M2,M1))
      RCLFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCLH2(M3,M2,M1))
      RC3FHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCO3H2(M3,M2,M1))
      RHCFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZHCOH2(M3,M2,M1))
      RAL1HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL1H2(M3,M2,M1))
      RAL2HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL2H2(M3,M2,M1))
      RAL3HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL3H2(M3,M2,M1))
      RAL4HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZAL4H2(M3,M2,M1))
      RALSHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZALSH2(M3,M2,M1))
      RFE1HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE1H2(M3,M2,M1))
      RFE2HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE2H2(M3,M2,M1))
      RFE3HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE3H2(M3,M2,M1))
      RFE4HS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFE4H2(M3,M2,M1))
      RFESHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZFESH2(M3,M2,M1))
      RCAOHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCAOH2(M3,M2,M1))
      RCACHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCACH2(M3,M2,M1))
      RCAHHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCAHH2(M3,M2,M1))
      RCASHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZCASH2(M3,M2,M1))
      RMGOHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGOH2(M3,M2,M1))
      RMGCHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGCH2(M3,M2,M1))
      RMGHHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGHH2(M3,M2,M1))
      RMGSHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZMGSH2(M3,M2,M1))
      RNACHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNACH2(M3,M2,M1))
      RNASHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNASH2(M3,M2,M1))
      RKASHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZKASH2(M3,M2,M1))
      RH0PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H0P4H2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RH3PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H3P4H2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RF1PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZF1PH2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RF2PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZF2PH2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RC0PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC0PH2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RC1PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC1PH2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RC2PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC2PH2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RM1PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZM1PH2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RH0BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H0PBH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RH3BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H3PBH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RF1BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZF1BH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RF2BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZF2BH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RC0BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC0BH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RC1BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC1BH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RC2BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZC2BH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RM1BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZM1BH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
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
      XALFLS(N,M6,M5,M4)=XALFLS(N,M6,M5,M4)+RALFLS(N,M6,M5,M4)
      XFEFLS(N,M6,M5,M4)=XFEFLS(N,M6,M5,M4)+RFEFLS(N,M6,M5,M4)
      XHYFLS(N,M6,M5,M4)=XHYFLS(N,M6,M5,M4)+RHYFLS(N,M6,M5,M4)
      XCAFLS(N,M6,M5,M4)=XCAFLS(N,M6,M5,M4)+RCAFLS(N,M6,M5,M4)
      XMGFLS(N,M6,M5,M4)=XMGFLS(N,M6,M5,M4)+RMGFLS(N,M6,M5,M4)
      XNAFLS(N,M6,M5,M4)=XNAFLS(N,M6,M5,M4)+RNAFLS(N,M6,M5,M4)
      XKAFLS(N,M6,M5,M4)=XKAFLS(N,M6,M5,M4)+RKAFLS(N,M6,M5,M4)
      XOHFLS(N,M6,M5,M4)=XOHFLS(N,M6,M5,M4)+ROHFLS(N,M6,M5,M4)
      XSOFLS(N,M6,M5,M4)=XSOFLS(N,M6,M5,M4)+RSOFLS(N,M6,M5,M4)
      XCLFLS(N,M6,M5,M4)=XCLFLS(N,M6,M5,M4)+RCLFLS(N,M6,M5,M4)
      XC3FLS(N,M6,M5,M4)=XC3FLS(N,M6,M5,M4)+RC3FLS(N,M6,M5,M4)
      XHCFLS(N,M6,M5,M4)=XHCFLS(N,M6,M5,M4)+RHCFLS(N,M6,M5,M4)
      XAL1FS(N,M6,M5,M4)=XAL1FS(N,M6,M5,M4)+RAL1FS(N,M6,M5,M4)
      XAL2FS(N,M6,M5,M4)=XAL2FS(N,M6,M5,M4)+RAL2FS(N,M6,M5,M4)
      XAL3FS(N,M6,M5,M4)=XAL3FS(N,M6,M5,M4)+RAL3FS(N,M6,M5,M4)
      XAL4FS(N,M6,M5,M4)=XAL4FS(N,M6,M5,M4)+RAL4FS(N,M6,M5,M4)
      XALSFS(N,M6,M5,M4)=XALSFS(N,M6,M5,M4)+RALSFS(N,M6,M5,M4)
      XFE1FS(N,M6,M5,M4)=XFE1FS(N,M6,M5,M4)+RFE1FS(N,M6,M5,M4)
      XFE2FS(N,M6,M5,M4)=XFE2FS(N,M6,M5,M4)+RFE2FS(N,M6,M5,M4)
      XFE3FS(N,M6,M5,M4)=XFE3FS(N,M6,M5,M4)+RFE3FS(N,M6,M5,M4)
      XFE4FS(N,M6,M5,M4)=XFE4FS(N,M6,M5,M4)+RFE4FS(N,M6,M5,M4)
      XFESFS(N,M6,M5,M4)=XFESFS(N,M6,M5,M4)+RFESFS(N,M6,M5,M4)
      XCAOFS(N,M6,M5,M4)=XCAOFS(N,M6,M5,M4)+RCAOFS(N,M6,M5,M4)
      XCACFS(N,M6,M5,M4)=XCACFS(N,M6,M5,M4)+RCACFS(N,M6,M5,M4)
      XCAHFS(N,M6,M5,M4)=XCAHFS(N,M6,M5,M4)+RCAHFS(N,M6,M5,M4)
      XCASFS(N,M6,M5,M4)=XCASFS(N,M6,M5,M4)+RCASFS(N,M6,M5,M4)
      XMGOFS(N,M6,M5,M4)=XMGOFS(N,M6,M5,M4)+RMGOFS(N,M6,M5,M4)
      XMGCFS(N,M6,M5,M4)=XMGCFS(N,M6,M5,M4)+RMGCFS(N,M6,M5,M4)
      XMGHFS(N,M6,M5,M4)=XMGHFS(N,M6,M5,M4)+RMGHFS(N,M6,M5,M4)
      XMGSFS(N,M6,M5,M4)=XMGSFS(N,M6,M5,M4)+RMGSFS(N,M6,M5,M4)
      XNACFS(N,M6,M5,M4)=XNACFS(N,M6,M5,M4)+RNACFS(N,M6,M5,M4)
      XNASFS(N,M6,M5,M4)=XNASFS(N,M6,M5,M4)+RNASFS(N,M6,M5,M4)
      XKASFS(N,M6,M5,M4)=XKASFS(N,M6,M5,M4)+RKASFS(N,M6,M5,M4)
      XH0PFS(N,M6,M5,M4)=XH0PFS(N,M6,M5,M4)+RH0PFS(N,M6,M5,M4)
      XH3PFS(N,M6,M5,M4)=XH3PFS(N,M6,M5,M4)+RH3PFS(N,M6,M5,M4)
      XF1PFS(N,M6,M5,M4)=XF1PFS(N,M6,M5,M4)+RF1PFS(N,M6,M5,M4)
      XF2PFS(N,M6,M5,M4)=XF2PFS(N,M6,M5,M4)+RF2PFS(N,M6,M5,M4)
      XC0PFS(N,M6,M5,M4)=XC0PFS(N,M6,M5,M4)+RC0PFS(N,M6,M5,M4)
      XC1PFS(N,M6,M5,M4)=XC1PFS(N,M6,M5,M4)+RC1PFS(N,M6,M5,M4)
      XC2PFS(N,M6,M5,M4)=XC2PFS(N,M6,M5,M4)+RC2PFS(N,M6,M5,M4)
      XM1PFS(N,M6,M5,M4)=XM1PFS(N,M6,M5,M4)+RM1PFS(N,M6,M5,M4)
      XH0BFB(N,M6,M5,M4)=XH0BFB(N,M6,M5,M4)+RH0BFB(N,M6,M5,M4)
      XH3BFB(N,M6,M5,M4)=XH3BFB(N,M6,M5,M4)+RH3BFB(N,M6,M5,M4)
      XF1BFB(N,M6,M5,M4)=XF1BFB(N,M6,M5,M4)+RF1BFB(N,M6,M5,M4)
      XF2BFB(N,M6,M5,M4)=XF2BFB(N,M6,M5,M4)+RF2BFB(N,M6,M5,M4)
      XC0BFB(N,M6,M5,M4)=XC0BFB(N,M6,M5,M4)+RC0BFB(N,M6,M5,M4)
      XC1BFB(N,M6,M5,M4)=XC1BFB(N,M6,M5,M4)+RC1BFB(N,M6,M5,M4)
      XC2BFB(N,M6,M5,M4)=XC2BFB(N,M6,M5,M4)+RC2BFB(N,M6,M5,M4)
      XM1BFB(N,M6,M5,M4)=XM1BFB(N,M6,M5,M4)+RM1BFB(N,M6,M5,M4)
      XALFHS(N,M6,M5,M4)=XALFHS(N,M6,M5,M4)+RALFHS(N,M6,M5,M4)
      XFEFHS(N,M6,M5,M4)=XFEFHS(N,M6,M5,M4)+RFEFHS(N,M6,M5,M4)
      XHYFHS(N,M6,M5,M4)=XHYFHS(N,M6,M5,M4)+RHYFHS(N,M6,M5,M4)
      XCAFHS(N,M6,M5,M4)=XCAFHS(N,M6,M5,M4)+RCAFHS(N,M6,M5,M4)
      XMGFHS(N,M6,M5,M4)=XMGFHS(N,M6,M5,M4)+RMGFHS(N,M6,M5,M4)
      XNAFHS(N,M6,M5,M4)=XNAFHS(N,M6,M5,M4)+RNAFHS(N,M6,M5,M4)
      XKAFHS(N,M6,M5,M4)=XKAFHS(N,M6,M5,M4)+RKAFHS(N,M6,M5,M4)
      XOHFHS(N,M6,M5,M4)=XOHFHS(N,M6,M5,M4)+ROHFHS(N,M6,M5,M4)
      XSOFHS(N,M6,M5,M4)=XSOFHS(N,M6,M5,M4)+RSOFHS(N,M6,M5,M4)
      XCLFHS(N,M6,M5,M4)=XCLFHS(N,M6,M5,M4)+RCLFHS(N,M6,M5,M4)
      XC3FHS(N,M6,M5,M4)=XC3FHS(N,M6,M5,M4)+RC3FHS(N,M6,M5,M4)
      XHCFHS(N,M6,M5,M4)=XHCFHS(N,M6,M5,M4)+RHCFHS(N,M6,M5,M4)
      XAL1HS(N,M6,M5,M4)=XAL1HS(N,M6,M5,M4)+RAL1HS(N,M6,M5,M4)
      XAL2HS(N,M6,M5,M4)=XAL2HS(N,M6,M5,M4)+RAL2HS(N,M6,M5,M4)
      XAL3HS(N,M6,M5,M4)=XAL3HS(N,M6,M5,M4)+RAL3HS(N,M6,M5,M4)
      XAL4HS(N,M6,M5,M4)=XAL4HS(N,M6,M5,M4)+RAL4HS(N,M6,M5,M4)
      XALSHS(N,M6,M5,M4)=XALSHS(N,M6,M5,M4)+RALSHS(N,M6,M5,M4)
      XFE1HS(N,M6,M5,M4)=XFE1HS(N,M6,M5,M4)+RFE1HS(N,M6,M5,M4)
      XFE2HS(N,M6,M5,M4)=XFE2HS(N,M6,M5,M4)+RFE2HS(N,M6,M5,M4)
      XFE3HS(N,M6,M5,M4)=XFE3HS(N,M6,M5,M4)+RFE3HS(N,M6,M5,M4)
      XFE4HS(N,M6,M5,M4)=XFE4HS(N,M6,M5,M4)+RFE4HS(N,M6,M5,M4)
      XFESHS(N,M6,M5,M4)=XFESHS(N,M6,M5,M4)+RFESHS(N,M6,M5,M4)
      XCAOHS(N,M6,M5,M4)=XCAOHS(N,M6,M5,M4)+RCAOHS(N,M6,M5,M4)
      XCACHS(N,M6,M5,M4)=XCACHS(N,M6,M5,M4)+RCACHS(N,M6,M5,M4)
      XCAHHS(N,M6,M5,M4)=XCAHHS(N,M6,M5,M4)+RCAHHS(N,M6,M5,M4)
      XCASHS(N,M6,M5,M4)=XCASHS(N,M6,M5,M4)+RCASHS(N,M6,M5,M4)
      XMGOHS(N,M6,M5,M4)=XMGOHS(N,M6,M5,M4)+RMGOHS(N,M6,M5,M4)
      XMGCHS(N,M6,M5,M4)=XMGCHS(N,M6,M5,M4)+RMGCHS(N,M6,M5,M4)
      XMGHHS(N,M6,M5,M4)=XMGHHS(N,M6,M5,M4)+RMGHHS(N,M6,M5,M4)
      XMGSHS(N,M6,M5,M4)=XMGSHS(N,M6,M5,M4)+RMGSHS(N,M6,M5,M4)
      XNACHS(N,M6,M5,M4)=XNACHS(N,M6,M5,M4)+RNACHS(N,M6,M5,M4)
      XNASHS(N,M6,M5,M4)=XNASHS(N,M6,M5,M4)+RNASHS(N,M6,M5,M4)
      XKASHS(N,M6,M5,M4)=XKASHS(N,M6,M5,M4)+RKASHS(N,M6,M5,M4)
      XH0PHS(N,M6,M5,M4)=XH0PHS(N,M6,M5,M4)+RH0PHS(N,M6,M5,M4)
      XH3PHS(N,M6,M5,M4)=XH3PHS(N,M6,M5,M4)+RH3PHS(N,M6,M5,M4)
      XF1PHS(N,M6,M5,M4)=XF1PHS(N,M6,M5,M4)+RF1PHS(N,M6,M5,M4)
      XF2PHS(N,M6,M5,M4)=XF2PHS(N,M6,M5,M4)+RF2PHS(N,M6,M5,M4)
      XC0PHS(N,M6,M5,M4)=XC0PHS(N,M6,M5,M4)+RC0PHS(N,M6,M5,M4)
      XC1PHS(N,M6,M5,M4)=XC1PHS(N,M6,M5,M4)+RC1PHS(N,M6,M5,M4)
      XC2PHS(N,M6,M5,M4)=XC2PHS(N,M6,M5,M4)+RC2PHS(N,M6,M5,M4)
      XM1PHS(N,M6,M5,M4)=XM1PHS(N,M6,M5,M4)+RM1PHS(N,M6,M5,M4)
      XH0BHB(N,M6,M5,M4)=XH0BHB(N,M6,M5,M4)+RH0BHB(N,M6,M5,M4)
      XH3BHB(N,M6,M5,M4)=XH3BHB(N,M6,M5,M4)+RH3BHB(N,M6,M5,M4)
      XF1BHB(N,M6,M5,M4)=XF1BHB(N,M6,M5,M4)+RF1BHB(N,M6,M5,M4)
      XF2BHB(N,M6,M5,M4)=XF2BHB(N,M6,M5,M4)+RF2BHB(N,M6,M5,M4)
      XC0BHB(N,M6,M5,M4)=XC0BHB(N,M6,M5,M4)+RC0BHB(N,M6,M5,M4)
      XC1BHB(N,M6,M5,M4)=XC1BHB(N,M6,M5,M4)+RC1BHB(N,M6,M5,M4)
      XC2BHB(N,M6,M5,M4)=XC2BHB(N,M6,M5,M4)+RC2BHB(N,M6,M5,M4)
      XM1BHB(N,M6,M5,M4)=XM1BHB(N,M6,M5,M4)+RM1BHB(N,M6,M5,M4)
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
      DO 1202 NN=1,2
      TQRAL(N2,N1)=TQRAL(N2,N1)+RQRAL(N,NN,N2,N1)
      TQRFE(N2,N1)=TQRFE(N2,N1)+RQRFE(N,NN,N2,N1)
      TQRHY(N2,N1)=TQRHY(N2,N1)+RQRHY(N,NN,N2,N1)
      TQRCA(N2,N1)=TQRCA(N2,N1)+RQRCA(N,NN,N2,N1)
      TQRMG(N2,N1)=TQRMG(N2,N1)+RQRMG(N,NN,N2,N1)
      TQRNA(N2,N1)=TQRNA(N2,N1)+RQRNA(N,NN,N2,N1)
      TQRKA(N2,N1)=TQRKA(N2,N1)+RQRKA(N,NN,N2,N1)
      TQROH(N2,N1)=TQROH(N2,N1)+RQROH(N,NN,N2,N1)
      TQRSO(N2,N1)=TQRSO(N2,N1)+RQRSO(N,NN,N2,N1)
      TQRCL(N2,N1)=TQRCL(N2,N1)+RQRCL(N,NN,N2,N1)
      TQRC3(N2,N1)=TQRC3(N2,N1)+RQRC3(N,NN,N2,N1)
      TQRHC(N2,N1)=TQRHC(N2,N1)+RQRHC(N,NN,N2,N1)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)+RQRAL1(N,NN,N2,N1)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)+RQRAL2(N,NN,N2,N1)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)+RQRAL3(N,NN,N2,N1)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)+RQRAL4(N,NN,N2,N1)
      TQRALS(N2,N1)=TQRALS(N2,N1)+RQRALS(N,NN,N2,N1)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)+RQRFE1(N,NN,N2,N1)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)+RQRFE2(N,NN,N2,N1)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)+RQRFE3(N,NN,N2,N1)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)+RQRFE4(N,NN,N2,N1)
      TQRFES(N2,N1)=TQRFES(N2,N1)+RQRFES(N,NN,N2,N1)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)+RQRCAO(N,NN,N2,N1)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)+RQRCAC(N,NN,N2,N1)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)+RQRCAH(N,NN,N2,N1)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)+RQRCAS(N,NN,N2,N1)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)+RQRMGO(N,NN,N2,N1)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)+RQRMGC(N,NN,N2,N1)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)+RQRMGH(N,NN,N2,N1)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)+RQRMGS(N,NN,N2,N1)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)+RQRNAC(N,NN,N2,N1)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)+RQRNAS(N,NN,N2,N1)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)+RQRKAS(N,NN,N2,N1)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)+RQRH0P(N,NN,N2,N1)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)+RQRH3P(N,NN,N2,N1)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)+RQRF1P(N,NN,N2,N1)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)+RQRF2P(N,NN,N2,N1)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)+RQRC0P(N,NN,N2,N1)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)+RQRC1P(N,NN,N2,N1)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)+RQRC2P(N,NN,N2,N1)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)+RQRM1P(N,NN,N2,N1)
      IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      TQRAL(N2,N1)=TQRAL(N2,N1)-RQRAL(N,NN,N5,N4)
      TQRFE(N2,N1)=TQRFE(N2,N1)-RQRFE(N,NN,N5,N4)
      TQRHY(N2,N1)=TQRHY(N2,N1)-RQRHY(N,NN,N5,N4)
      TQRCA(N2,N1)=TQRCA(N2,N1)-RQRCA(N,NN,N5,N4)
      TQRMG(N2,N1)=TQRMG(N2,N1)-RQRMG(N,NN,N5,N4)
      TQRNA(N2,N1)=TQRNA(N2,N1)-RQRNA(N,NN,N5,N4)
      TQRKA(N2,N1)=TQRKA(N2,N1)-RQRKA(N,NN,N5,N4)
      TQROH(N2,N1)=TQROH(N2,N1)-RQROH(N,NN,N5,N4)
      TQRSO(N2,N1)=TQRSO(N2,N1)-RQRSO(N,NN,N5,N4)
      TQRCL(N2,N1)=TQRCL(N2,N1)-RQRCL(N,NN,N5,N4)
      TQRC3(N2,N1)=TQRC3(N2,N1)-RQRC3(N,NN,N5,N4)
      TQRHC(N2,N1)=TQRHC(N2,N1)-RQRHC(N,NN,N5,N4)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)-RQRAL1(N,NN,N5,N4)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)-RQRAL2(N,NN,N5,N4)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)-RQRAL3(N,NN,N5,N4)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)-RQRAL4(N,NN,N5,N4)
      TQRALS(N2,N1)=TQRALS(N2,N1)-RQRALS(N,NN,N5,N4)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)-RQRFE1(N,NN,N5,N4)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)-RQRFE2(N,NN,N5,N4)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)-RQRFE3(N,NN,N5,N4)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)-RQRFE4(N,NN,N5,N4)
      TQRFES(N2,N1)=TQRFES(N2,N1)-RQRFES(N,NN,N5,N4)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)-RQRCAO(N,NN,N5,N4)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)-RQRCAC(N,NN,N5,N4)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)-RQRCAH(N,NN,N5,N4)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)-RQRCAS(N,NN,N5,N4)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)-RQRMGO(N,NN,N5,N4)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)-RQRMGC(N,NN,N5,N4)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)-RQRMGH(N,NN,N5,N4)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)-RQRMGS(N,NN,N5,N4)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)-RQRNAC(N,NN,N5,N4)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)-RQRNAS(N,NN,N5,N4)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)-RQRKAS(N,NN,N5,N4)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)-RQRH0P(N,NN,N5,N4)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)-RQRH3P(N,NN,N5,N4)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)-RQRF1P(N,NN,N5,N4)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)-RQRF2P(N,NN,N5,N4)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)-RQRC0P(N,NN,N5,N4)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)-RQRC1P(N,NN,N5,N4)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)-RQRC2P(N,NN,N5,N4)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)-RQRM1P(N,NN,N5,N4)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQRAL(N2,N1)=TQRAL(N2,N1)-RQRAL(N,NN,N5B,N4B)
      TQRFE(N2,N1)=TQRFE(N2,N1)-RQRFE(N,NN,N5B,N4B)
      TQRHY(N2,N1)=TQRHY(N2,N1)-RQRHY(N,NN,N5B,N4B)
      TQRCA(N2,N1)=TQRCA(N2,N1)-RQRCA(N,NN,N5B,N4B)
      TQRMG(N2,N1)=TQRMG(N2,N1)-RQRMG(N,NN,N5B,N4B)
      TQRNA(N2,N1)=TQRNA(N2,N1)-RQRNA(N,NN,N5B,N4B)
      TQRKA(N2,N1)=TQRKA(N2,N1)-RQRKA(N,NN,N5B,N4B)
      TQROH(N2,N1)=TQROH(N2,N1)-RQROH(N,NN,N5B,N4B)
      TQRSO(N2,N1)=TQRSO(N2,N1)-RQRSO(N,NN,N5B,N4B)
      TQRCL(N2,N1)=TQRCL(N2,N1)-RQRCL(N,NN,N5B,N4B)
      TQRC3(N2,N1)=TQRC3(N2,N1)-RQRC3(N,NN,N5B,N4B)
      TQRHC(N2,N1)=TQRHC(N2,N1)-RQRHC(N,NN,N5B,N4B)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)-RQRAL1(N,NN,N5B,N4B)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)-RQRAL2(N,NN,N5B,N4B)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)-RQRAL3(N,NN,N5B,N4B)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)-RQRAL4(N,NN,N5B,N4B)
      TQRALS(N2,N1)=TQRALS(N2,N1)-RQRALS(N,NN,N5B,N4B)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)-RQRFE1(N,NN,N5B,N4B)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)-RQRFE2(N,NN,N5B,N4B)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)-RQRFE3(N,NN,N5B,N4B)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)-RQRFE4(N,NN,N5B,N4B)
      TQRFES(N2,N1)=TQRFES(N2,N1)-RQRFES(N,NN,N5B,N4B)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)-RQRCAO(N,NN,N5B,N4B)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)-RQRCAC(N,NN,N5B,N4B)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)-RQRCAH(N,NN,N5B,N4B)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)-RQRCAS(N,NN,N5B,N4B)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)-RQRMGO(N,NN,N5B,N4B)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)-RQRMGC(N,NN,N5B,N4B)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)-RQRMGH(N,NN,N5B,N4B)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)-RQRMGS(N,NN,N5B,N4B)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)-RQRNAC(N,NN,N5B,N4B)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)-RQRNAS(N,NN,N5B,N4B)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)-RQRKAS(N,NN,N5B,N4B)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)-RQRH0P(N,NN,N5B,N4B)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)-RQRH3P(N,NN,N5B,N4B)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)-RQRF1P(N,NN,N5B,N4B)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)-RQRF2P(N,NN,N5B,N4B)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)-RQRC0P(N,NN,N5B,N4B)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)-RQRC1P(N,NN,N5B,N4B)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)-RQRC2P(N,NN,N5B,N4B)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)-RQRM1P(N,NN,N5B,N4B)
      ENDIF
1202  CONTINUE
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
      TQSAL(N2,N1)=TQSAL(N2,N1)+RQSAL(N,N2,N1)-RQSAL(N,N5,N4)
      TQSFE(N2,N1)=TQSFE(N2,N1)+RQSFE(N,N2,N1)-RQSFE(N,N5,N4)
      TQSHY(N2,N1)=TQSHY(N2,N1)+RQSHY(N,N2,N1)-RQSHY(N,N5,N4)
      TQSCA(N2,N1)=TQSCA(N2,N1)+RQSCA(N,N2,N1)-RQSCA(N,N5,N4)
      TQSMG(N2,N1)=TQSMG(N2,N1)+RQSMG(N,N2,N1)-RQSMG(N,N5,N4)
      TQSNA(N2,N1)=TQSNA(N2,N1)+RQSNA(N,N2,N1)-RQSNA(N,N5,N4)
      TQSKA(N2,N1)=TQSKA(N2,N1)+RQSKA(N,N2,N1)-RQSKA(N,N5,N4)
      TQSOH(N2,N1)=TQSOH(N2,N1)+RQSOH(N,N2,N1)-RQSOH(N,N5,N4)
      TQSSO(N2,N1)=TQSSO(N2,N1)+RQSSO(N,N2,N1)-RQSSO(N,N5,N4)
      TQSCL(N2,N1)=TQSCL(N2,N1)+RQSCL(N,N2,N1)-RQSCL(N,N5,N4)
      TQSC3(N2,N1)=TQSC3(N2,N1)+RQSC3(N,N2,N1)-RQSC3(N,N5,N4)
      TQSHC(N2,N1)=TQSHC(N2,N1)+RQSHC(N,N2,N1)-RQSHC(N,N5,N4)
      TQSAL1(N2,N1)=TQSAL1(N2,N1)+RQSAL1(N,N2,N1)-RQSAL1(N,N5,N4)
      TQSAL2(N2,N1)=TQSAL2(N2,N1)+RQSAL2(N,N2,N1)-RQSAL2(N,N5,N4)
      TQSAL3(N2,N1)=TQSAL3(N2,N1)+RQSAL3(N,N2,N1)-RQSAL3(N,N5,N4)
      TQSAL4(N2,N1)=TQSAL4(N2,N1)+RQSAL4(N,N2,N1)-RQSAL4(N,N5,N4)
      TQSALS(N2,N1)=TQSALS(N2,N1)+RQSALS(N,N2,N1)-RQSALS(N,N5,N4)
      TQSFE1(N2,N1)=TQSFE1(N2,N1)+RQSFE1(N,N2,N1)-RQSFE1(N,N5,N4)
      TQSFE2(N2,N1)=TQSFE2(N2,N1)+RQSFE2(N,N2,N1)-RQSFE2(N,N5,N4)
      TQSFE3(N2,N1)=TQSFE3(N2,N1)+RQSFE3(N,N2,N1)-RQSFE3(N,N5,N4)
      TQSFE4(N2,N1)=TQSFE4(N2,N1)+RQSFE4(N,N2,N1)-RQSFE4(N,N5,N4)
      TQSFES(N2,N1)=TQSFES(N2,N1)+RQSFES(N,N2,N1)-RQSFES(N,N5,N4)
      TQSCAO(N2,N1)=TQSCAO(N2,N1)+RQSCAO(N,N2,N1)-RQSCAO(N,N5,N4)
      TQSCAC(N2,N1)=TQSCAC(N2,N1)+RQSCAC(N,N2,N1)-RQSCAC(N,N5,N4)
      TQSCAH(N2,N1)=TQSCAH(N2,N1)+RQSCAH(N,N2,N1)-RQSCAH(N,N5,N4)
      TQSCAS(N2,N1)=TQSCAS(N2,N1)+RQSCAS(N,N2,N1)-RQSCAS(N,N5,N4)
      TQSMGO(N2,N1)=TQSMGO(N2,N1)+RQSMGO(N,N2,N1)-RQSMGO(N,N5,N4)
      TQSMGC(N2,N1)=TQSMGC(N2,N1)+RQSMGC(N,N2,N1)-RQSMGC(N,N5,N4)
      TQSMGH(N2,N1)=TQSMGH(N2,N1)+RQSMGH(N,N2,N1)-RQSMGH(N,N5,N4)
      TQSMGS(N2,N1)=TQSMGS(N2,N1)+RQSMGS(N,N2,N1)-RQSMGS(N,N5,N4)
      TQSNAC(N2,N1)=TQSNAC(N2,N1)+RQSNAC(N,N2,N1)-RQSNAC(N,N5,N4)
      TQSNAS(N2,N1)=TQSNAS(N2,N1)+RQSNAS(N,N2,N1)-RQSNAS(N,N5,N4)
      TQSKAS(N2,N1)=TQSKAS(N2,N1)+RQSKAS(N,N2,N1)-RQSKAS(N,N5,N4)
      TQSH0P(N2,N1)=TQSH0P(N2,N1)+RQSH0P(N,N2,N1)-RQSH0P(N,N5,N4)
      TQSH3P(N2,N1)=TQSH3P(N2,N1)+RQSH3P(N,N2,N1)-RQSH3P(N,N5,N4)
      TQSF1P(N2,N1)=TQSF1P(N2,N1)+RQSF1P(N,N2,N1)-RQSF1P(N,N5,N4)
      TQSF2P(N2,N1)=TQSF2P(N2,N1)+RQSF2P(N,N2,N1)-RQSF2P(N,N5,N4)
      TQSC0P(N2,N1)=TQSC0P(N2,N1)+RQSC0P(N,N2,N1)-RQSC0P(N,N5,N4)
      TQSC1P(N2,N1)=TQSC1P(N2,N1)+RQSC1P(N,N2,N1)-RQSC1P(N,N5,N4)
      TQSC2P(N2,N1)=TQSC2P(N2,N1)+RQSC2P(N,N2,N1)-RQSC2P(N,N5,N4)
      TQSM1P(N2,N1)=TQSM1P(N2,N1)+RQSM1P(N,N2,N1)-RQSM1P(N,N5,N4)
      end subroutine NetOverloadFLuxInSnow
!------------------------------------------------------------------------------------------

  subroutine NetFluxInSnowpack(M,NY,NX,N1,N2)
!
!     Description:
!
      integer, intent(in) :: M,NY,NX,N1,N2
      integer :: LS,LS2
!     begin_execution

      DO 1205 LS=1,JS
      IF(VHCPWM(M,LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
  IF(LS.LT.JS.AND.VHCPWM(M,LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
    TALBLS(LS,NY,NX)=TALBLS(LS,NY,NX)+RALBLS(LS,NY,NX) &
      -RALBLS(LS2,NY,NX)
    TFEBLS(LS,NY,NX)=TFEBLS(LS,NY,NX)+RFEBLS(LS,NY,NX) &
      -RFEBLS(LS2,NY,NX)
    THYBLS(LS,NY,NX)=THYBLS(LS,NY,NX)+RHYBLS(LS,NY,NX) &
      -RHYBLS(LS2,NY,NX)
    TCABLS(LS,NY,NX)=TCABLS(LS,NY,NX)+RCABLS(LS,NY,NX) &
      -RCABLS(LS2,NY,NX)
    TMGBLS(LS,NY,NX)=TMGBLS(LS,NY,NX)+RMGBLS(LS,NY,NX) &
      -RMGBLS(LS2,NY,NX)
    TNABLS(LS,NY,NX)=TNABLS(LS,NY,NX)+RNABLS(LS,NY,NX) &
      -RNABLS(LS2,NY,NX)
      TKABLS(LS,NY,NX)=TKABLS(LS,NY,NX)+RKABLS(LS,NY,NX) &
      -RKABLS(LS2,NY,NX)
      TOHBLS(LS,NY,NX)=TOHBLS(LS,NY,NX)+ROHBLS(LS,NY,NX) &
      -ROHBLS(LS2,NY,NX)
      TSOBLS(LS,NY,NX)=TSOBLS(LS,NY,NX)+RSOBLS(LS,NY,NX) &
      -RSOBLS(LS2,NY,NX)
      TCLBLS(LS,NY,NX)=TCLBLS(LS,NY,NX)+RCLBLS(LS,NY,NX) &
      -RCLBLS(LS2,NY,NX)
      TC3BLS(LS,NY,NX)=TC3BLS(LS,NY,NX)+RC3BLS(LS,NY,NX) &
      -RC3BLS(LS2,NY,NX)
      THCBLS(LS,NY,NX)=THCBLS(LS,NY,NX)+RHCBLS(LS,NY,NX) &
      -RHCBLS(LS2,NY,NX)
      TAL1BS(LS,NY,NX)=TAL1BS(LS,NY,NX)+RAL1BS(LS,NY,NX) &
      -RAL1BS(LS2,NY,NX)
      TAL2BS(LS,NY,NX)=TAL2BS(LS,NY,NX)+RAL2BS(LS,NY,NX) &
      -RAL2BS(LS2,NY,NX)
      TAL3BS(LS,NY,NX)=TAL3BS(LS,NY,NX)+RAL3BS(LS,NY,NX) &
      -RAL3BS(LS2,NY,NX)
      TAL4BS(LS,NY,NX)=TAL4BS(LS,NY,NX)+RAL4BS(LS,NY,NX) &
      -RAL4BS(LS2,NY,NX)
      TALSBS(LS,NY,NX)=TALSBS(LS,NY,NX)+RALSBS(LS,NY,NX) &
      -RALSBS(LS2,NY,NX)
      TFE1BS(LS,NY,NX)=TFE1BS(LS,NY,NX)+RFE1BS(LS,NY,NX) &
      -RFE1BS(LS2,NY,NX)
      TFE2BS(LS,NY,NX)=TFE2BS(LS,NY,NX)+RFE2BS(LS,NY,NX) &
      -RFE2BS(LS2,NY,NX)
      TFE3BS(LS,NY,NX)=TFE3BS(LS,NY,NX)+RFE3BS(LS,NY,NX) &
      -RFE3BS(LS2,NY,NX)
      TFE4BS(LS,NY,NX)=TFE4BS(LS,NY,NX)+RFE4BS(LS,NY,NX) &
      -RFE4BS(LS2,NY,NX)
      TFESBS(LS,NY,NX)=TFESBS(LS,NY,NX)+RFESBS(LS,NY,NX) &
      -RFESBS(LS2,NY,NX)
      TCAOBS(LS,NY,NX)=TCAOBS(LS,NY,NX)+RCAOBS(LS,NY,NX) &
      -RCAOBS(LS2,NY,NX)
      TCACBS(LS,NY,NX)=TCACBS(LS,NY,NX)+RCACBS(LS,NY,NX) &
      -RCACBS(LS2,NY,NX)
      TCAHBS(LS,NY,NX)=TCAHBS(LS,NY,NX)+RCAHBS(LS,NY,NX) &
      -RCAHBS(LS2,NY,NX)
      TCASBS(LS,NY,NX)=TCASBS(LS,NY,NX)+RCASBS(LS,NY,NX) &
      -RCASBS(LS2,NY,NX)
      TMGOBS(LS,NY,NX)=TMGOBS(LS,NY,NX)+RMGOBS(LS,NY,NX) &
      -RMGOBS(LS2,NY,NX)
      TMGCBS(LS,NY,NX)=TMGCBS(LS,NY,NX)+RMGCBS(LS,NY,NX) &
      -RMGCBS(LS2,NY,NX)
      TMGHBS(LS,NY,NX)=TMGHBS(LS,NY,NX)+RMGHBS(LS,NY,NX) &
      -RMGHBS(LS2,NY,NX)
      TMGSBS(LS,NY,NX)=TMGSBS(LS,NY,NX)+RMGSBS(LS,NY,NX) &
      -RMGSBS(LS2,NY,NX)
      TNACBS(LS,NY,NX)=TNACBS(LS,NY,NX)+RNACBS(LS,NY,NX) &
      -RNACBS(LS2,NY,NX)
      TNASBS(LS,NY,NX)=TNASBS(LS,NY,NX)+RNASBS(LS,NY,NX) &
      -RNASBS(LS2,NY,NX)
      TKASBS(LS,NY,NX)=TKASBS(LS,NY,NX)+RKASBS(LS,NY,NX) &
      -RKASBS(LS2,NY,NX)
      TH0PBS(LS,NY,NX)=TH0PBS(LS,NY,NX)+RH0PBS(LS,NY,NX) &
      -RH0PBS(LS2,NY,NX)
      TH3PBS(LS,NY,NX)=TH3PBS(LS,NY,NX)+RH3PBS(LS,NY,NX) &
      -RH3PBS(LS2,NY,NX)
      TF1PBS(LS,NY,NX)=TF1PBS(LS,NY,NX)+RF1PBS(LS,NY,NX) &
      -RF1PBS(LS2,NY,NX)
      TF2PBS(LS,NY,NX)=TF2PBS(LS,NY,NX)+RF2PBS(LS,NY,NX) &
      -RF2PBS(LS2,NY,NX)
      TC0PBS(LS,NY,NX)=TC0PBS(LS,NY,NX)+RC0PBS(LS,NY,NX) &
      -RC0PBS(LS2,NY,NX)
      TC1PBS(LS,NY,NX)=TC1PBS(LS,NY,NX)+RC1PBS(LS,NY,NX) &
      -RC1PBS(LS2,NY,NX)
      TC2PBS(LS,NY,NX)=TC2PBS(LS,NY,NX)+RC2PBS(LS,NY,NX) &
      -RC2PBS(LS2,NY,NX)
      TM1PBS(LS,NY,NX)=TM1PBS(LS,NY,NX)+RM1PBS(LS,NY,NX) &
      -RM1PBS(LS2,NY,NX)
      ELSE
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
      TALBLS(LS,NY,NX)=TALBLS(LS,NY,NX)+RALBLS(LS,NY,NX) &
      -RALFLS(3,0,N2,N1)-RALFLS(3,NUM(N2,N1),N2,N1)
!    3-RALFHS(3,NUM(N2,N1),N2,N1)
      TFEBLS(LS,NY,NX)=TFEBLS(LS,NY,NX)+RFEBLS(LS,NY,NX) &
      -RFEFLS(3,0,N2,N1)-RFEFLS(3,NUM(N2,N1),N2,N1)
!    3-RFEFHS(3,NUM(N2,N1),N2,N1)
      THYBLS(LS,NY,NX)=THYBLS(LS,NY,NX)+RHYBLS(LS,NY,NX) &
      -RHYFLS(3,0,N2,N1)-RHYFLS(3,NUM(N2,N1),N2,N1)
!    3-RHYFHS(3,NUM(N2,N1),N2,N1)
      TCABLS(LS,NY,NX)=TCABLS(LS,NY,NX)+RCABLS(LS,NY,NX) &
      -RCAFLS(3,0,N2,N1)-RCAFLS(3,NUM(N2,N1),N2,N1)
!    3-RCAFHS(3,NUM(N2,N1),N2,N1)
      TMGBLS(LS,NY,NX)=TMGBLS(LS,NY,NX)+RMGBLS(LS,NY,NX) &
      -RMGFLS(3,0,N2,N1)-RMGFLS(3,NUM(N2,N1),N2,N1)
!    3-RMGFHS(3,NUM(N2,N1),N2,N1)
      TNABLS(LS,NY,NX)=TNABLS(LS,NY,NX)+RNABLS(LS,NY,NX) &
      -RNAFLS(3,0,N2,N1)-RNAFLS(3,NUM(N2,N1),N2,N1)
!    3-RNAFHS(3,NUM(N2,N1),N2,N1)
      TKABLS(LS,NY,NX)=TKABLS(LS,NY,NX)+RKABLS(LS,NY,NX) &
      -RKAFLS(3,0,N2,N1)-RKAFLS(3,NUM(N2,N1),N2,N1)
!    3-RKAFHS(3,NUM(N2,N1),N2,N1)
      TOHBLS(LS,NY,NX)=TOHBLS(LS,NY,NX)+ROHBLS(LS,NY,NX) &
      -ROHFLS(3,0,N2,N1)-ROHFLS(3,NUM(N2,N1),N2,N1)
!    3-ROHFHS(3,NUM(N2,N1),N2,N1)
      TSOBLS(LS,NY,NX)=TSOBLS(LS,NY,NX)+RSOBLS(LS,NY,NX) &
      -RSOFLS(3,0,N2,N1)-RSOFLS(3,NUM(N2,N1),N2,N1)
!    3-RSOFHS(3,NUM(N2,N1),N2,N1)
      TCLBLS(LS,NY,NX)=TCLBLS(LS,NY,NX)+RCLBLS(LS,NY,NX) &
      -RCLFLS(3,0,N2,N1)-RCLFLS(3,NUM(N2,N1),N2,N1)
!    3-RCLFHS(3,NUM(N2,N1),N2,N1)
      TC3BLS(LS,NY,NX)=TC3BLS(LS,NY,NX)+RC3BLS(LS,NY,NX) &
      -RC3FLS(3,0,N2,N1)-RC3FLS(3,NUM(N2,N1),N2,N1)
!    3-RC3FHS(3,NUM(N2,N1),N2,N1)
      THCBLS(LS,NY,NX)=THCBLS(LS,NY,NX)+RHCBLS(LS,NY,NX) &
      -RHCFLS(3,0,N2,N1)-RHCFLS(3,NUM(N2,N1),N2,N1)
!    3-RHCFHS(3,NUM(N2,N1),N2,N1)
      TAL1BS(LS,NY,NX)=TAL1BS(LS,NY,NX)+RAL1BS(LS,NY,NX) &
      -RAL1FS(3,0,N2,N1)-RAL1FS(3,NUM(N2,N1),N2,N1)
!    3-RAL1HS(3,NUM(N2,N1),N2,N1)
      TAL2BS(LS,NY,NX)=TAL2BS(LS,NY,NX)+RAL2BS(LS,NY,NX) &
      -RAL2FS(3,0,N2,N1)-RAL2FS(3,NUM(N2,N1),N2,N1)
!    3-RAL2HS(3,NUM(N2,N1),N2,N1)
      TAL3BS(LS,NY,NX)=TAL3BS(LS,NY,NX)+RAL3BS(LS,NY,NX) &
      -RAL3FS(3,0,N2,N1)-RAL3FS(3,NUM(N2,N1),N2,N1)
!    3-RAL3HS(3,NUM(N2,N1),N2,N1)
      TAL4BS(LS,NY,NX)=TAL4BS(LS,NY,NX)+RAL4BS(LS,NY,NX) &
      -RAL4FS(3,0,N2,N1)-RAL4FS(3,NUM(N2,N1),N2,N1)
!    3-RAL4HS(3,NUM(N2,N1),N2,N1)
      TALSBS(LS,NY,NX)=TALSBS(LS,NY,NX)+RALSBS(LS,NY,NX) &
      -RALSFS(3,0,N2,N1)-RALSFS(3,NUM(N2,N1),N2,N1)
!    3-RALSHS(3,NUM(N2,N1),N2,N1)
      TFE1BS(LS,NY,NX)=TFE1BS(LS,NY,NX)+RFE1BS(LS,NY,NX) &
      -RFE1FS(3,0,N2,N1)-RFE1FS(3,NUM(N2,N1),N2,N1)
!    3-RFE1HS(3,NUM(N2,N1),N2,N1)
      TFE2BS(LS,NY,NX)=TFE2BS(LS,NY,NX)+RFE2BS(LS,NY,NX) &
      -RFE2FS(3,0,N2,N1)-RFE2FS(3,NUM(N2,N1),N2,N1)
!    3-RFE2HS(3,NUM(N2,N1),N2,N1)
      TFE3BS(LS,NY,NX)=TFE3BS(LS,NY,NX)+RFE3BS(LS,NY,NX) &
      -RFE3FS(3,0,N2,N1)-RFE3FS(3,NUM(N2,N1),N2,N1)
!    3-RFE3HS(3,NUM(N2,N1),N2,N1)
      TFE4BS(LS,NY,NX)=TFE4BS(LS,NY,NX)+RFE4BS(LS,NY,NX) &
      -RFE4FS(3,0,N2,N1)-RFE4FS(3,NUM(N2,N1),N2,N1)
!    3-RFE4HS(3,NUM(N2,N1),N2,N1)
      TFESBS(LS,NY,NX)=TFESBS(LS,NY,NX)+RFESBS(LS,NY,NX) &
      -RFESFS(3,0,N2,N1)-RFESFS(3,NUM(N2,N1),N2,N1)
!    3-RFESHS(3,NUM(N2,N1),N2,N1)
      TCAOBS(LS,NY,NX)=TCAOBS(LS,NY,NX)+RCAOBS(LS,NY,NX) &
      -RCAOFS(3,0,N2,N1)-RCAOFS(3,NUM(N2,N1),N2,N1)
!    3-RCAOHS(3,NUM(N2,N1),N2,N1)
      TCACBS(LS,NY,NX)=TCACBS(LS,NY,NX)+RCACBS(LS,NY,NX) &
      -RCACFS(3,0,N2,N1)-RCACFS(3,NUM(N2,N1),N2,N1)
!    3-RCACHS(3,NUM(N2,N1),N2,N1)
      TCAHBS(LS,NY,NX)=TCAHBS(LS,NY,NX)+RCAHBS(LS,NY,NX) &
      -RCAHFS(3,0,N2,N1)-RCAHFS(3,NUM(N2,N1),N2,N1)
!    3-RCAHHS(3,NUM(N2,N1),N2,N1)
      TCASBS(LS,NY,NX)=TCASBS(LS,NY,NX)+RCASBS(LS,NY,NX) &
      -RCASFS(3,0,N2,N1)-RCASFS(3,NUM(N2,N1),N2,N1)
!    3-RCASHS(3,NUM(N2,N1),N2,N1)
      TMGOBS(LS,NY,NX)=TMGOBS(LS,NY,NX)+RMGOBS(LS,NY,NX) &
      -RMGOFS(3,0,N2,N1)-RMGOFS(3,NUM(N2,N1),N2,N1)
!    3-RMGOHS(3,NUM(N2,N1),N2,N1)
      TMGCBS(LS,NY,NX)=TMGCBS(LS,NY,NX)+RMGCBS(LS,NY,NX) &
      -RMGCFS(3,0,N2,N1)-RMGCFS(3,NUM(N2,N1),N2,N1)
!    3-RMGCHS(3,NUM(N2,N1),N2,N1)
      TMGHBS(LS,NY,NX)=TMGHBS(LS,NY,NX)+RMGHBS(LS,NY,NX) &
      -RMGHFS(3,0,N2,N1)-RMGHFS(3,NUM(N2,N1),N2,N1)
!    3-RMGHHS(3,NUM(N2,N1),N2,N1)
      TMGSBS(LS,NY,NX)=TMGSBS(LS,NY,NX)+RMGSBS(LS,NY,NX) &
      -RMGSFS(3,0,N2,N1)-RMGSFS(3,NUM(N2,N1),N2,N1)
!    3-RMGSHS(3,NUM(N2,N1),N2,N1)
      TNACBS(LS,NY,NX)=TNACBS(LS,NY,NX)+RNACBS(LS,NY,NX) &
      -RNACFS(3,0,N2,N1)-RNACFS(3,NUM(N2,N1),N2,N1)
!    3-RNACHS(3,NUM(N2,N1),N2,N1)
      TNASBS(LS,NY,NX)=TNASBS(LS,NY,NX)+RNASBS(LS,NY,NX) &
      -RNASFS(3,0,N2,N1)-RNASFS(3,NUM(N2,N1),N2,N1)
!    3-RNASHS(3,NUM(N2,N1),N2,N1)
      TKASBS(LS,NY,NX)=TKASBS(LS,NY,NX)+RKASBS(LS,NY,NX) &
      -RKASFS(3,0,N2,N1)-RKASFS(3,NUM(N2,N1),N2,N1)
!    3-RKASHS(3,NUM(N2,N1),N2,N1)
      TH0PBS(LS,NY,NX)=TH0PBS(LS,NY,NX)+RH0PBS(LS,NY,NX) &
      -RH0PFS(3,0,N2,N1)-RH0PFS(3,NUM(N2,N1),N2,N1) &
      -RH0BFB(3,NUM(N2,N1),N2,N1)
!    3-RH0PHS(3,NUM(N2,N1),N2,N1)-RH0BFB(3,NUM(N2,N1),N2,N1)
      TH3PBS(LS,NY,NX)=TH3PBS(LS,NY,NX)+RH3PBS(LS,NY,NX) &
      -RH3PFS(3,0,N2,N1)-RH3PHS(3,NUM(N2,N1),N2,N1) &
      -RH3BFB(3,NUM(N2,N1),N2,N1)
!    3-RH3PHS(3,NUM(N2,N1),N2,N1)-RH3BFB(3,NUM(N2,N1),N2,N1)
      TF1PBS(LS,NY,NX)=TF1PBS(LS,NY,NX)+RF1PBS(LS,NY,NX) &
      -RF1PFS(3,0,N2,N1)-RF1PFS(3,NUM(N2,N1),N2,N1) &
      -RF1BFB(3,NUM(N2,N1),N2,N1)
!    3-RF1PHS(3,NUM(N2,N1),N2,N1)-RF1BFB(3,NUM(N2,N1),N2,N1)
      TF2PBS(LS,NY,NX)=TF2PBS(LS,NY,NX)+RF2PBS(LS,NY,NX) &
      -RF2PFS(3,0,N2,N1)-RF2PFS(3,NUM(N2,N1),N2,N1) &
      -RF2BFB(3,NUM(N2,N1),N2,N1)
!    3-RF2PHS(3,NUM(N2,N1),N2,N1)-RF2BFB(3,NUM(N2,N1),N2,N1)
      TC0PBS(LS,NY,NX)=TC0PBS(LS,NY,NX)+RC0PBS(LS,NY,NX) &
      -RC0PFS(3,0,N2,N1)-RC0PFS(3,NUM(N2,N1),N2,N1) &
      -RC0BFB(3,NUM(N2,N1),N2,N1)
!    3-RC0PHS(3,NUM(N2,N1),N2,N1)-RC0BFB(3,NUM(N2,N1),N2,N1)
      TC1PBS(LS,NY,NX)=TC1PBS(LS,NY,NX)+RC1PBS(LS,NY,NX) &
      -RC1PFS(3,0,N2,N1)-RC1PFS(3,NUM(N2,N1),N2,N1) &
      -RC1BFB(3,NUM(N2,N1),N2,N1)
!    3-RC1PHS(3,NUM(N2,N1),N2,N1)-RC1BFB(3,NUM(N2,N1),N2,N1)
      TC2PBS(LS,NY,NX)=TC2PBS(LS,NY,NX)+RC2PBS(LS,NY,NX) &
      -RC2PFS(3,0,N2,N1)-RC2PFS(3,NUM(N2,N1),N2,N1) &
      -RC2BFB(3,NUM(N2,N1),N2,N1)
!    3-RC2PHS(3,NUM(N2,N1),N2,N1)-RC2BFB(3,NUM(N2,N1),N2,N1)
      TM1PBS(LS,NY,NX)=TM1PBS(LS,NY,NX)+RM1PBS(LS,NY,NX) &
      -RM1PFS(3,0,N2,N1)-RM1PFS(3,NUM(N2,N1),N2,N1) &
      -RM1BFB(3,NUM(N2,N1),N2,N1)
!    3-RM1PHS(3,NUM(N2,N1),N2,N1)-RM1BFB(3,NUM(N2,N1),N2,N1)
      ENDIF
      ENDIF
1205  CONTINUE
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
      DO 1200 LL=N6,NL(NY,NX)
      IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      GO TO 1201
      ENDIF
1200  CONTINUE
1201  CONTINUE
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      TALFLS(N3,N2,N1)=TALFLS(N3,N2,N1)+RALFLS(N,N3,N2,N1) &
      -RALFLS(N,N6,N5,N4)
      TFEFLS(N3,N2,N1)=TFEFLS(N3,N2,N1)+RFEFLS(N,N3,N2,N1) &
      -RFEFLS(N,N6,N5,N4)
      THYFLS(N3,N2,N1)=THYFLS(N3,N2,N1)+RHYFLS(N,N3,N2,N1) &
      -RHYFLS(N,N6,N5,N4)
      TCAFLS(N3,N2,N1)=TCAFLS(N3,N2,N1)+RCAFLS(N,N3,N2,N1) &
      -RCAFLS(N,N6,N5,N4)
      TMGFLS(N3,N2,N1)=TMGFLS(N3,N2,N1)+RMGFLS(N,N3,N2,N1) &
      -RMGFLS(N,N6,N5,N4)
      TNAFLS(N3,N2,N1)=TNAFLS(N3,N2,N1)+RNAFLS(N,N3,N2,N1) &
      -RNAFLS(N,N6,N5,N4)
      TKAFLS(N3,N2,N1)=TKAFLS(N3,N2,N1)+RKAFLS(N,N3,N2,N1) &
      -RKAFLS(N,N6,N5,N4)
      TOHFLS(N3,N2,N1)=TOHFLS(N3,N2,N1)+ROHFLS(N,N3,N2,N1) &
      -ROHFLS(N,N6,N5,N4)
      TSOFLS(N3,N2,N1)=TSOFLS(N3,N2,N1)+RSOFLS(N,N3,N2,N1) &
      -RSOFLS(N,N6,N5,N4)
      TCLFLS(N3,N2,N1)=TCLFLS(N3,N2,N1)+RCLFLS(N,N3,N2,N1) &
      -RCLFLS(N,N6,N5,N4)
      TC3FLS(N3,N2,N1)=TC3FLS(N3,N2,N1)+RC3FLS(N,N3,N2,N1) &
      -RC3FLS(N,N6,N5,N4)
      THCFLS(N3,N2,N1)=THCFLS(N3,N2,N1)+RHCFLS(N,N3,N2,N1) &
      -RHCFLS(N,N6,N5,N4)
      TAL1FS(N3,N2,N1)=TAL1FS(N3,N2,N1)+RAL1FS(N,N3,N2,N1) &
      -RAL1FS(N,N6,N5,N4)
      TAL2FS(N3,N2,N1)=TAL2FS(N3,N2,N1)+RAL2FS(N,N3,N2,N1) &
      -RAL2FS(N,N6,N5,N4)
      TAL3FS(N3,N2,N1)=TAL3FS(N3,N2,N1)+RAL3FS(N,N3,N2,N1) &
      -RAL3FS(N,N6,N5,N4)
      TAL4FS(N3,N2,N1)=TAL4FS(N3,N2,N1)+RAL4FS(N,N3,N2,N1) &
      -RAL4FS(N,N6,N5,N4)
      TALSFS(N3,N2,N1)=TALSFS(N3,N2,N1)+RALSFS(N,N3,N2,N1) &
      -RALSFS(N,N6,N5,N4)
      TFE1FS(N3,N2,N1)=TFE1FS(N3,N2,N1)+RFE1FS(N,N3,N2,N1) &
      -RFE1FS(N,N6,N5,N4)
      TFE2FS(N3,N2,N1)=TFE2FS(N3,N2,N1)+RFE2FS(N,N3,N2,N1) &
      -RFE2FS(N,N6,N5,N4)
      TFE3FS(N3,N2,N1)=TFE3FS(N3,N2,N1)+RFE3FS(N,N3,N2,N1) &
      -RFE3FS(N,N6,N5,N4)
      TFE4FS(N3,N2,N1)=TFE4FS(N3,N2,N1)+RFE4FS(N,N3,N2,N1) &
      -RFE4FS(N,N6,N5,N4)
      TFESFS(N3,N2,N1)=TFESFS(N3,N2,N1)+RFESFS(N,N3,N2,N1) &
      -RFESFS(N,N6,N5,N4)
      TCAOFS(N3,N2,N1)=TCAOFS(N3,N2,N1)+RCAOFS(N,N3,N2,N1) &
      -RCAOFS(N,N6,N5,N4)
      TCACFS(N3,N2,N1)=TCACFS(N3,N2,N1)+RCACFS(N,N3,N2,N1) &
      -RCACFS(N,N6,N5,N4)
      TCAHFS(N3,N2,N1)=TCAHFS(N3,N2,N1)+RCAHFS(N,N3,N2,N1) &
      -RCAHFS(N,N6,N5,N4)
      TCASFS(N3,N2,N1)=TCASFS(N3,N2,N1)+RCASFS(N,N3,N2,N1) &
      -RCASFS(N,N6,N5,N4)
      TMGOFS(N3,N2,N1)=TMGOFS(N3,N2,N1)+RMGOFS(N,N3,N2,N1) &
      -RMGOFS(N,N6,N5,N4)
      TMGCFS(N3,N2,N1)=TMGCFS(N3,N2,N1)+RMGCFS(N,N3,N2,N1) &
      -RMGCFS(N,N6,N5,N4)
      TMGHFS(N3,N2,N1)=TMGHFS(N3,N2,N1)+RMGHFS(N,N3,N2,N1) &
      -RMGHFS(N,N6,N5,N4)
      TMGSFS(N3,N2,N1)=TMGSFS(N3,N2,N1)+RMGSFS(N,N3,N2,N1) &
      -RMGSFS(N,N6,N5,N4)
      TNACFS(N3,N2,N1)=TNACFS(N3,N2,N1)+RNACFS(N,N3,N2,N1) &
      -RNACFS(N,N6,N5,N4)
      TNASFS(N3,N2,N1)=TNASFS(N3,N2,N1)+RNASFS(N,N3,N2,N1) &
      -RNASFS(N,N6,N5,N4)
      TKASFS(N3,N2,N1)=TKASFS(N3,N2,N1)+RKASFS(N,N3,N2,N1) &
      -RKASFS(N,N6,N5,N4)
      TH0PFS(N3,N2,N1)=TH0PFS(N3,N2,N1)+RH0PFS(N,N3,N2,N1) &
      -RH0PFS(N,N6,N5,N4)
      TH3PFS(N3,N2,N1)=TH3PFS(N3,N2,N1)+RH3PFS(N,N3,N2,N1) &
      -RH3PFS(N,N6,N5,N4)
      TF1PFS(N3,N2,N1)=TF1PFS(N3,N2,N1)+RF1PFS(N,N3,N2,N1) &
      -RF1PFS(N,N6,N5,N4)
      TF2PFS(N3,N2,N1)=TF2PFS(N3,N2,N1)+RF2PFS(N,N3,N2,N1) &
      -RF2PFS(N,N6,N5,N4)
      TC0PFS(N3,N2,N1)=TC0PFS(N3,N2,N1)+RC0PFS(N,N3,N2,N1) &
      -RC0PFS(N,N6,N5,N4)
      TC1PFS(N3,N2,N1)=TC1PFS(N3,N2,N1)+RC1PFS(N,N3,N2,N1) &
      -RC1PFS(N,N6,N5,N4)
      TC2PFS(N3,N2,N1)=TC2PFS(N3,N2,N1)+RC2PFS(N,N3,N2,N1) &
      -RC2PFS(N,N6,N5,N4)
      TM1PFS(N3,N2,N1)=TM1PFS(N3,N2,N1)+RM1PFS(N,N3,N2,N1) &
      -RM1PFS(N,N6,N5,N4)
      TH0BFB(N3,N2,N1)=TH0BFB(N3,N2,N1)+RH0BFB(N,N3,N2,N1) &
      -RH0BFB(N,N6,N5,N4)
      TH3BFB(N3,N2,N1)=TH3BFB(N3,N2,N1)+RH3BFB(N,N3,N2,N1) &
      -RH3BFB(N,N6,N5,N4)
      TF1BFB(N3,N2,N1)=TF1BFB(N3,N2,N1)+RF1BFB(N,N3,N2,N1) &
      -RF1BFB(N,N6,N5,N4)
      TF2BFB(N3,N2,N1)=TF2BFB(N3,N2,N1)+RF2BFB(N,N3,N2,N1) &
      -RF2BFB(N,N6,N5,N4)
      TC0BFB(N3,N2,N1)=TC0BFB(N3,N2,N1)+RC0BFB(N,N3,N2,N1) &
      -RC0BFB(N,N6,N5,N4)
      TC1BFB(N3,N2,N1)=TC1BFB(N3,N2,N1)+RC1BFB(N,N3,N2,N1) &
      -RC1BFB(N,N6,N5,N4)
      TC2BFB(N3,N2,N1)=TC2BFB(N3,N2,N1)+RC2BFB(N,N3,N2,N1) &
      -RC2BFB(N,N6,N5,N4)
      TM1BFB(N3,N2,N1)=TM1BFB(N3,N2,N1)+RM1BFB(N,N3,N2,N1) &
      -RM1BFB(N,N6,N5,N4)
      TALFHS(N3,N2,N1)=TALFHS(N3,N2,N1)+RALFHS(N,N3,N2,N1) &
      -RALFHS(N,N6,N5,N4)
      TFEFHS(N3,N2,N1)=TFEFHS(N3,N2,N1)+RFEFHS(N,N3,N2,N1) &
      -RFEFHS(N,N6,N5,N4)
      THYFHS(N3,N2,N1)=THYFHS(N3,N2,N1)+RHYFHS(N,N3,N2,N1) &
      -RHYFHS(N,N6,N5,N4)
      TCAFHS(N3,N2,N1)=TCAFHS(N3,N2,N1)+RCAFHS(N,N3,N2,N1) &
      -RCAFHS(N,N6,N5,N4)
      TMGFHS(N3,N2,N1)=TMGFHS(N3,N2,N1)+RMGFHS(N,N3,N2,N1) &
      -RMGFHS(N,N6,N5,N4)
      TNAFHS(N3,N2,N1)=TNAFHS(N3,N2,N1)+RNAFHS(N,N3,N2,N1) &
      -RNAFHS(N,N6,N5,N4)
      TKAFHS(N3,N2,N1)=TKAFHS(N3,N2,N1)+RKAFHS(N,N3,N2,N1) &
      -RKAFHS(N,N6,N5,N4)
      TOHFHS(N3,N2,N1)=TOHFHS(N3,N2,N1)+ROHFHS(N,N3,N2,N1) &
      -ROHFHS(N,N6,N5,N4)
      TSOFHS(N3,N2,N1)=TSOFHS(N3,N2,N1)+RSOFHS(N,N3,N2,N1) &
      -RSOFHS(N,N6,N5,N4)
      TCLFHS(N3,N2,N1)=TCLFHS(N3,N2,N1)+RCLFHS(N,N3,N2,N1) &
      -RCLFHS(N,N6,N5,N4)
      TC3FHS(N3,N2,N1)=TC3FHS(N3,N2,N1)+RC3FHS(N,N3,N2,N1) &
      -RC3FHS(N,N6,N5,N4)
      THCFHS(N3,N2,N1)=THCFHS(N3,N2,N1)+RHCFHS(N,N3,N2,N1) &
      -RHCFHS(N,N6,N5,N4)
      TAL1HS(N3,N2,N1)=TAL1HS(N3,N2,N1)+RAL1HS(N,N3,N2,N1) &
      -RAL1HS(N,N6,N5,N4)
      TAL2HS(N3,N2,N1)=TAL2HS(N3,N2,N1)+RAL2HS(N,N3,N2,N1) &
      -RAL2HS(N,N6,N5,N4)
      TAL3HS(N3,N2,N1)=TAL3HS(N3,N2,N1)+RAL3HS(N,N3,N2,N1) &
      -RAL3HS(N,N6,N5,N4)
      TAL4HS(N3,N2,N1)=TAL4HS(N3,N2,N1)+RAL4HS(N,N3,N2,N1) &
      -RAL4HS(N,N6,N5,N4)
      TALSHS(N3,N2,N1)=TALSHS(N3,N2,N1)+RALSHS(N,N3,N2,N1) &
      -RALSHS(N,N6,N5,N4)
      TFE1HS(N3,N2,N1)=TFE1HS(N3,N2,N1)+RFE1HS(N,N3,N2,N1) &
      -RFE1HS(N,N6,N5,N4)
      TFE2HS(N3,N2,N1)=TFE2HS(N3,N2,N1)+RFE2HS(N,N3,N2,N1) &
      -RFE2HS(N,N6,N5,N4)
      TFE3HS(N3,N2,N1)=TFE3HS(N3,N2,N1)+RFE3HS(N,N3,N2,N1) &
      -RFE3HS(N,N6,N5,N4)
      TFE4HS(N3,N2,N1)=TFE4HS(N3,N2,N1)+RFE4HS(N,N3,N2,N1) &
      -RFE4HS(N,N6,N5,N4)
      TFESHS(N3,N2,N1)=TFESHS(N3,N2,N1)+RFESHS(N,N3,N2,N1) &
      -RFESHS(N,N6,N5,N4)
      TCAOHS(N3,N2,N1)=TCAOHS(N3,N2,N1)+RCAOHS(N,N3,N2,N1) &
      -RCAOHS(N,N6,N5,N4)
      TCACHS(N3,N2,N1)=TCACHS(N3,N2,N1)+RCACHS(N,N3,N2,N1) &
      -RCACHS(N,N6,N5,N4)
      TCAHHS(N3,N2,N1)=TCAHHS(N3,N2,N1)+RCAHHS(N,N3,N2,N1) &
      -RCAHHS(N,N6,N5,N4)
      TCASHS(N3,N2,N1)=TCASHS(N3,N2,N1)+RCASHS(N,N3,N2,N1) &
      -RCASHS(N,N6,N5,N4)
      TMGOHS(N3,N2,N1)=TMGOHS(N3,N2,N1)+RMGOHS(N,N3,N2,N1) &
      -RMGOHS(N,N6,N5,N4)
      TMGCHS(N3,N2,N1)=TMGCHS(N3,N2,N1)+RMGCHS(N,N3,N2,N1) &
      -RMGCHS(N,N6,N5,N4)
      TMGHHS(N3,N2,N1)=TMGHHS(N3,N2,N1)+RMGHHS(N,N3,N2,N1) &
      -RMGHHS(N,N6,N5,N4)
      TMGSHS(N3,N2,N1)=TMGSHS(N3,N2,N1)+RMGSHS(N,N3,N2,N1) &
      -RMGSHS(N,N6,N5,N4)
      TNACHS(N3,N2,N1)=TNACHS(N3,N2,N1)+RNACHS(N,N3,N2,N1) &
      -RNACHS(N,N6,N5,N4)
      TNASHS(N3,N2,N1)=TNASHS(N3,N2,N1)+RNASHS(N,N3,N2,N1) &
      -RNASHS(N,N6,N5,N4)
      TKASHS(N3,N2,N1)=TKASHS(N3,N2,N1)+RKASHS(N,N3,N2,N1) &
      -RKASHS(N,N6,N5,N4)
      TH0PHS(N3,N2,N1)=TH0PHS(N3,N2,N1)+RH0PHS(N,N3,N2,N1) &
      -RH0PHS(N,N6,N5,N4)
      TH3PHS(N3,N2,N1)=TH3PHS(N3,N2,N1)+RH3PHS(N,N3,N2,N1) &
      -RH3PHS(N,N6,N5,N4)
      TF1PHS(N3,N2,N1)=TF1PHS(N3,N2,N1)+RF1PHS(N,N3,N2,N1) &
      -RF1PHS(N,N6,N5,N4)
      TF2PHS(N3,N2,N1)=TF2PHS(N3,N2,N1)+RF2PHS(N,N3,N2,N1) &
      -RF2PHS(N,N6,N5,N4)
      TC0PHS(N3,N2,N1)=TC0PHS(N3,N2,N1)+RC0PHS(N,N3,N2,N1) &
      -RC0PHS(N,N6,N5,N4)
      TC1PHS(N3,N2,N1)=TC1PHS(N3,N2,N1)+RC1PHS(N,N3,N2,N1) &
      -RC1PHS(N,N6,N5,N4)
      TC2PHS(N3,N2,N1)=TC2PHS(N3,N2,N1)+RC2PHS(N,N3,N2,N1) &
      -RC2PHS(N,N6,N5,N4)
      TM1PHS(N3,N2,N1)=TM1PHS(N3,N2,N1)+RM1PHS(N,N3,N2,N1) &
      -RM1PHS(N,N6,N5,N4)
      TH0BHB(N3,N2,N1)=TH0BHB(N3,N2,N1)+RH0BHB(N,N3,N2,N1) &
      -RH0BHB(N,N6,N5,N4)
      TH3BHB(N3,N2,N1)=TH3BHB(N3,N2,N1)+RH3BHB(N,N3,N2,N1) &
      -RH3BHB(N,N6,N5,N4)
      TF1BHB(N3,N2,N1)=TF1BHB(N3,N2,N1)+RF1BHB(N,N3,N2,N1) &
      -RF1BHB(N,N6,N5,N4)
      TF2BHB(N3,N2,N1)=TF2BHB(N3,N2,N1)+RF2BHB(N,N3,N2,N1) &
      -RF2BHB(N,N6,N5,N4)
      TC0BHB(N3,N2,N1)=TC0BHB(N3,N2,N1)+RC0BHB(N,N3,N2,N1) &
      -RC0BHB(N,N6,N5,N4)
      TC1BHB(N3,N2,N1)=TC1BHB(N3,N2,N1)+RC1BHB(N,N3,N2,N1) &
      -RC1BHB(N,N6,N5,N4)
      TC2BHB(N3,N2,N1)=TC2BHB(N3,N2,N1)+RC2BHB(N,N3,N2,N1) &
      -RC2BHB(N,N6,N5,N4)
      TM1BHB(N3,N2,N1)=TM1BHB(N3,N2,N1)+RM1BHB(N,N3,N2,N1) &
      -RM1BHB(N,N6,N5,N4)
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

      DO 9595 NX=NHW,NHE
      DO 9590 NY=NVN,NVS
      DO 9585 L=NU(NY,NX),NL(NY,NX)
      N1=NX
      N2=NY
      N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
      DO 9580 N=1,3
      DO 9575 NN=1,2
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
      GO TO 9575
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
      GO TO 9575
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
      GO TO 9575
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
      GO TO 9575
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
      GO TO 9575
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      GO TO 9575
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
9575  CONTINUE
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
9580  CONTINUE
9585  CONTINUE
9590  CONTINUE
9595  CONTINUE
      end subroutine SaltModelInternalFlux
!------------------------------------------------------------------------------------------

      subroutine UpdateSoluteInSnow(NY,NX)
!
!     Description:
!
      integer, intent(in) :: NY,NX
      integer :: L
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
      ZALW2(1,NY,NX)=ZALW2(1,NY,NX)+TQSAL(NY,NX)
      ZFEW2(1,NY,NX)=ZFEW2(1,NY,NX)+TQSFE(NY,NX)
      ZHYW2(1,NY,NX)=ZHYW2(1,NY,NX)+TQSHY(NY,NX)
      ZCAW2(1,NY,NX)=ZCAW2(1,NY,NX)+TQSCA(NY,NX)
      ZMGW2(1,NY,NX)=ZMGW2(1,NY,NX)+TQSMG(NY,NX)
      ZNAW2(1,NY,NX)=ZNAW2(1,NY,NX)+TQSNA(NY,NX)
      ZKAW2(1,NY,NX)=ZKAW2(1,NY,NX)+TQSKA(NY,NX)
      ZOHW2(1,NY,NX)=ZOHW2(1,NY,NX)+TQSOH(NY,NX)
      ZSO4W2(1,NY,NX)=ZSO4W2(1,NY,NX)+TQSSO(NY,NX)
      ZCLW2(1,NY,NX)=ZCLW2(1,NY,NX)+TQSCL(NY,NX)
      ZCO3W2(1,NY,NX)=ZCO3W2(1,NY,NX)+TQSC3(NY,NX)
      ZHCO3W2(1,NY,NX)=ZHCO3W2(1,NY,NX)+TQSHC(NY,NX)
      ZALH1W2(1,NY,NX)=ZALH1W2(1,NY,NX)+TQSAL1(NY,NX)
      ZALH2W2(1,NY,NX)=ZALH2W2(1,NY,NX)+TQSAL2(NY,NX)
      ZALH3W2(1,NY,NX)=ZALH3W2(1,NY,NX)+TQSAL3(NY,NX)
      ZALH4W2(1,NY,NX)=ZALH4W2(1,NY,NX)+TQSAL4(NY,NX)
      ZALSW2(1,NY,NX)=ZALSW2(1,NY,NX)+TQSALS(NY,NX)
      ZFEH1W2(1,NY,NX)=ZFEH1W2(1,NY,NX)+TQSFE1(NY,NX)
      ZFEH2W2(1,NY,NX)=ZFEH2W2(1,NY,NX)+TQSFE2(NY,NX)
      ZFEH3W2(1,NY,NX)=ZFEH3W2(1,NY,NX)+TQSFE3(NY,NX)
      ZFEH4W2(1,NY,NX)=ZFEH4W2(1,NY,NX)+TQSFE4(NY,NX)
      ZFESW2(1,NY,NX)=ZFESW2(1,NY,NX)+TQSFES(NY,NX)
      ZCAOW2(1,NY,NX)=ZCAOW2(1,NY,NX)+TQSCAO(NY,NX)
      ZCACW2(1,NY,NX)=ZCACW2(1,NY,NX)+TQSCAC(NY,NX)
      ZCAHW2(1,NY,NX)=ZCAHW2(1,NY,NX)+TQSCAH(NY,NX)
      ZCASW2(1,NY,NX)=ZCASW2(1,NY,NX)+TQSCAS(NY,NX)
      ZMGOW2(1,NY,NX)=ZMGOW2(1,NY,NX)+TQSMGO(NY,NX)
      ZMGCW2(1,NY,NX)=ZMGCW2(1,NY,NX)+TQSMGC(NY,NX)
      ZMGHW2(1,NY,NX)=ZMGHW2(1,NY,NX)+TQSMGH(NY,NX)
      ZMGSW2(1,NY,NX)=ZMGSW2(1,NY,NX)+TQSMGS(NY,NX)
      ZNACW2(1,NY,NX)=ZNACW2(1,NY,NX)+TQSNAC(NY,NX)
      ZNASW2(1,NY,NX)=ZNASW2(1,NY,NX)+TQSNAS(NY,NX)
      ZKASW2(1,NY,NX)=ZKASW2(1,NY,NX)+TQSKAS(NY,NX)
      H0PO4W2(1,NY,NX)=H0PO4W2(1,NY,NX)+TQSH0P(NY,NX)
      H3PO4W2(1,NY,NX)=H3PO4W2(1,NY,NX)+TQSH3P(NY,NX)
      ZFE1PW2(1,NY,NX)=ZFE1PW2(1,NY,NX)+TQSF1P(NY,NX)
      ZFE2PW2(1,NY,NX)=ZFE2PW2(1,NY,NX)+TQSF2P(NY,NX)
      ZCA0PW2(1,NY,NX)=ZCA0PW2(1,NY,NX)+TQSC0P(NY,NX)
      ZCA1PW2(1,NY,NX)=ZCA1PW2(1,NY,NX)+TQSC1P(NY,NX)
      ZCA2PW2(1,NY,NX)=ZCA2PW2(1,NY,NX)+TQSC2P(NY,NX)
      ZMG1PW2(1,NY,NX)=ZMG1PW2(1,NY,NX)+TQSM1P(NY,NX)
      DO 9670 L=1,JS
      ZALW2(L,NY,NX)=ZALW2(L,NY,NX)+TALBLS(L,NY,NX)
      ZFEW2(L,NY,NX)=ZFEW2(L,NY,NX)+TFEBLS(L,NY,NX)
      ZHYW2(L,NY,NX)=ZHYW2(L,NY,NX)+THYBLS(L,NY,NX)
      ZCAW2(L,NY,NX)=ZCAW2(L,NY,NX)+TCABLS(L,NY,NX)
      ZMGW2(L,NY,NX)=ZMGW2(L,NY,NX)+TMGBLS(L,NY,NX)
      ZNAW2(L,NY,NX)=ZNAW2(L,NY,NX)+TNABLS(L,NY,NX)
      ZKAW2(L,NY,NX)=ZKAW2(L,NY,NX)+TKABLS(L,NY,NX)
      ZOHW2(L,NY,NX)=ZOHW2(L,NY,NX)+TOHBLS(L,NY,NX)
      ZSO4W2(L,NY,NX)=ZSO4W2(L,NY,NX)+TSOBLS(L,NY,NX)
      ZCLW2(L,NY,NX)=ZCLW2(L,NY,NX)+TCLBLS(L,NY,NX)
      ZCO3W2(L,NY,NX)=ZCO3W2(L,NY,NX)+TC3BLS(L,NY,NX)
      ZHCO3W2(L,NY,NX)=ZHCO3W2(L,NY,NX)+THCBLS(L,NY,NX)
      ZALH1W2(L,NY,NX)=ZALH1W2(L,NY,NX)+TAL1BS(L,NY,NX)
      ZALH2W2(L,NY,NX)=ZALH2W2(L,NY,NX)+TAL2BS(L,NY,NX)
      ZALH3W2(L,NY,NX)=ZALH3W2(L,NY,NX)+TAL3BS(L,NY,NX)
      ZALH4W2(L,NY,NX)=ZALH4W2(L,NY,NX)+TAL4BS(L,NY,NX)
      ZALSW2(L,NY,NX)=ZALSW2(L,NY,NX)+TALSBS(L,NY,NX)
      ZFEH1W2(L,NY,NX)=ZFEH1W2(L,NY,NX)+TFE1BS(L,NY,NX)
      ZFEH2W2(L,NY,NX)=ZFEH2W2(L,NY,NX)+TFE2BS(L,NY,NX)
      ZFEH3W2(L,NY,NX)=ZFEH3W2(L,NY,NX)+TFE3BS(L,NY,NX)
      ZFEH4W2(L,NY,NX)=ZFEH4W2(L,NY,NX)+TFE4BS(L,NY,NX)
      ZFESW2(L,NY,NX)=ZFESW2(L,NY,NX)+TFESBS(L,NY,NX)
      ZCAOW2(L,NY,NX)=ZCAOW2(L,NY,NX)+TCAOBS(L,NY,NX)
      ZCACW2(L,NY,NX)=ZCACW2(L,NY,NX)+TCACBS(L,NY,NX)
      ZCAHW2(L,NY,NX)=ZCAHW2(L,NY,NX)+TCAHBS(L,NY,NX)
      ZCASW2(L,NY,NX)=ZCASW2(L,NY,NX)+TCASBS(L,NY,NX)
      ZMGOW2(L,NY,NX)=ZMGOW2(L,NY,NX)+TMGOBS(L,NY,NX)
      ZMGCW2(L,NY,NX)=ZMGCW2(L,NY,NX)+TMGCBS(L,NY,NX)
      ZMGHW2(L,NY,NX)=ZMGHW2(L,NY,NX)+TMGHBS(L,NY,NX)
      ZMGSW2(L,NY,NX)=ZMGSW2(L,NY,NX)+TMGSBS(L,NY,NX)
      ZNACW2(L,NY,NX)=ZNACW2(L,NY,NX)+TNACBS(L,NY,NX)
      ZNASW2(L,NY,NX)=ZNASW2(L,NY,NX)+TNASBS(L,NY,NX)
      ZKASW2(L,NY,NX)=ZKASW2(L,NY,NX)+TKASBS(L,NY,NX)
      H0PO4W2(L,NY,NX)=H0PO4W2(L,NY,NX)+TH0PBS(L,NY,NX)
      H3PO4W2(L,NY,NX)=H3PO4W2(L,NY,NX)+TH3PBS(L,NY,NX)
      ZFE1PW2(L,NY,NX)=ZFE1PW2(L,NY,NX)+TF1PBS(L,NY,NX)
      ZFE2PW2(L,NY,NX)=ZFE2PW2(L,NY,NX)+TF2PBS(L,NY,NX)
      ZCA0PW2(L,NY,NX)=ZCA0PW2(L,NY,NX)+TC0PBS(L,NY,NX)
      ZCA1PW2(L,NY,NX)=ZCA1PW2(L,NY,NX)+TC1PBS(L,NY,NX)
      ZCA2PW2(L,NY,NX)=ZCA2PW2(L,NY,NX)+TC2PBS(L,NY,NX)
      ZMG1PW2(L,NY,NX)=ZMG1PW2(L,NY,NX)+TM1PBS(L,NY,NX)
9670  CONTINUE
      end subroutine UpdateSoluteInSnow
!------------------------------------------------------------------------------------------

      subroutine UpdateSoluteInResidue(NY,NX)
!
!     Description:
!
      integer, intent(in) :: NY,NX
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
      ZAL2(0,NY,NX)=ZAL2(0,NY,NX)+TQRAL(NY,NX)+RALFLS(3,0,NY,NX)
      ZFE2(0,NY,NX)=ZFE2(0,NY,NX)+TQRFE(NY,NX)+RFEFLS(3,0,NY,NX)
      ZHY2(0,NY,NX)=ZHY2(0,NY,NX)+TQRHY(NY,NX)+RHYFLS(3,0,NY,NX)
      ZCA2(0,NY,NX)=ZCA2(0,NY,NX)+TQRCA(NY,NX)+RCAFLS(3,0,NY,NX)
      ZMG2(0,NY,NX)=ZMG2(0,NY,NX)+TQRMG(NY,NX)+RMGFLS(3,0,NY,NX)
      ZNA2(0,NY,NX)=ZNA2(0,NY,NX)+TQRNA(NY,NX)+RNAFLS(3,0,NY,NX)
      ZKA2(0,NY,NX)=ZKA2(0,NY,NX)+TQRKA(NY,NX)+RKAFLS(3,0,NY,NX)
      ZOH2(0,NY,NX)=ZOH2(0,NY,NX)+TQROH(NY,NX)+ROHFLS(3,0,NY,NX)
      ZSO42(0,NY,NX)=ZSO42(0,NY,NX)+TQRSO(NY,NX)+RSOFLS(3,0,NY,NX)
      ZCL2(0,NY,NX)=ZCL2(0,NY,NX)+TQRCL(NY,NX)+RCLFLS(3,0,NY,NX)
      ZCO32(0,NY,NX)=ZCO32(0,NY,NX)+TQRC3(NY,NX)+RC3FLS(3,0,NY,NX)
      ZHCO32(0,NY,NX)=ZHCO32(0,NY,NX)+TQRHC(NY,NX)+RHCFLS(3,0,NY,NX)
      ZAL12(0,NY,NX)=ZAL12(0,NY,NX)+TQRAL1(NY,NX)+RAL1FS(3,0,NY,NX)
      ZAL22(0,NY,NX)=ZAL22(0,NY,NX)+TQRAL2(NY,NX)+RAL2FS(3,0,NY,NX)
      ZAL32(0,NY,NX)=ZAL32(0,NY,NX)+TQRAL3(NY,NX)+RAL3FS(3,0,NY,NX)
      ZAL42(0,NY,NX)=ZAL42(0,NY,NX)+TQRAL4(NY,NX)+RAL4FS(3,0,NY,NX)
      ZALS2(0,NY,NX)=ZALS2(0,NY,NX)+TQRALS(NY,NX)+RALSFS(3,0,NY,NX)
      ZFE12(0,NY,NX)=ZFE12(0,NY,NX)+TQRFE1(NY,NX)+RFE1FS(3,0,NY,NX)
      ZFE22(0,NY,NX)=ZFE22(0,NY,NX)+TQRFE2(NY,NX)+RFE2FS(3,0,NY,NX)
      ZFE32(0,NY,NX)=ZFE32(0,NY,NX)+TQRFE3(NY,NX)+RFE3FS(3,0,NY,NX)
      ZFE42(0,NY,NX)=ZFE42(0,NY,NX)+TQRFE4(NY,NX)+RFE4FS(3,0,NY,NX)
      ZFES2(0,NY,NX)=ZFES2(0,NY,NX)+TQRFES(NY,NX)+RFESFS(3,0,NY,NX)
      ZCAO2(0,NY,NX)=ZCAO2(0,NY,NX)+TQRCAO(NY,NX)+RCAOFS(3,0,NY,NX)
      ZCAC2(0,NY,NX)=ZCAC2(0,NY,NX)+TQRCAC(NY,NX)+RCACFS(3,0,NY,NX)
      ZCAH2(0,NY,NX)=ZCAH2(0,NY,NX)+TQRCAH(NY,NX)+RCAHFS(3,0,NY,NX)
      ZCAS2(0,NY,NX)=ZCAS2(0,NY,NX)+TQRCAS(NY,NX)+RCASFS(3,0,NY,NX)
      ZMGO2(0,NY,NX)=ZMGO2(0,NY,NX)+TQRMGO(NY,NX)+RMGOFS(3,0,NY,NX)
      ZMGC2(0,NY,NX)=ZMGC2(0,NY,NX)+TQRMGC(NY,NX)+RMGCFS(3,0,NY,NX)
      ZMGH2(0,NY,NX)=ZMGH2(0,NY,NX)+TQRMGH(NY,NX)+RMGHFS(3,0,NY,NX)
      ZMGS2(0,NY,NX)=ZMGS2(0,NY,NX)+TQRMGS(NY,NX)+RMGSFS(3,0,NY,NX)
      ZNAC2(0,NY,NX)=ZNAC2(0,NY,NX)+TQRNAC(NY,NX)+RNACFS(3,0,NY,NX)
      ZNAS2(0,NY,NX)=ZNAS2(0,NY,NX)+TQRNAS(NY,NX)+RNASFS(3,0,NY,NX)
      ZKAS2(0,NY,NX)=ZKAS2(0,NY,NX)+TQRKAS(NY,NX)+RKASFS(3,0,NY,NX)
      H0PO42(0,NY,NX)=H0PO42(0,NY,NX)+TQRH0P(NY,NX)+RH0PFS(3,0,NY,NX)
      H3PO42(0,NY,NX)=H3PO42(0,NY,NX)+TQRH3P(NY,NX)+RH3PFS(3,0,NY,NX)
      ZFE1P2(0,NY,NX)=ZFE1P2(0,NY,NX)+TQRF1P(NY,NX)+RF1PFS(3,0,NY,NX)
      ZFE2P2(0,NY,NX)=ZFE2P2(0,NY,NX)+TQRF2P(NY,NX)+RF2PFS(3,0,NY,NX)
      ZCA0P2(0,NY,NX)=ZCA0P2(0,NY,NX)+TQRC0P(NY,NX)+RC0PFS(3,0,NY,NX)
      ZCA1P2(0,NY,NX)=ZCA1P2(0,NY,NX)+TQRC1P(NY,NX)+RC1PFS(3,0,NY,NX)
      ZCA2P2(0,NY,NX)=ZCA2P2(0,NY,NX)+TQRC2P(NY,NX)+RC2PFS(3,0,NY,NX)
      ZMG1P2(0,NY,NX)=ZMG1P2(0,NY,NX)+TQRM1P(NY,NX)+RM1PFS(3,0,NY,NX)
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
      DO 9685 L=NU(NY,NX),NL(NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      ZAL2(L,NY,NX)=ZAL2(L,NY,NX)+TALFLS(L,NY,NX)+RALFXS(L,NY,NX) &
      +RALFLZ(L,NY,NX)
      ZFE2(L,NY,NX)=ZFE2(L,NY,NX)+TFEFLS(L,NY,NX)+RFEFXS(L,NY,NX) &
      +RFEFLZ(L,NY,NX)
      ZHY2(L,NY,NX)=ZHY2(L,NY,NX)+THYFLS(L,NY,NX)+RHYFXS(L,NY,NX) &
      +RHYFLZ(L,NY,NX)
      ZCA2(L,NY,NX)=ZCA2(L,NY,NX)+TCAFLS(L,NY,NX)+RCAFXS(L,NY,NX) &
      +RCAFLZ(L,NY,NX)
      ZMG2(L,NY,NX)=ZMG2(L,NY,NX)+TMGFLS(L,NY,NX)+RMGFXS(L,NY,NX) &
      +RMGFLZ(L,NY,NX)
      ZNA2(L,NY,NX)=ZNA2(L,NY,NX)+TNAFLS(L,NY,NX)+RNAFXS(L,NY,NX) &
      +RNAFLZ(L,NY,NX)
      ZKA2(L,NY,NX)=ZKA2(L,NY,NX)+TKAFLS(L,NY,NX)+RKAFXS(L,NY,NX) &
      +RKAFLZ(L,NY,NX)
      ZOH2(L,NY,NX)=ZOH2(L,NY,NX)+TOHFLS(L,NY,NX)+ROHFXS(L,NY,NX) &
      +ROHFLZ(L,NY,NX)
      ZSO42(L,NY,NX)=ZSO42(L,NY,NX)+TSOFLS(L,NY,NX)+RSOFXS(L,NY,NX) &
      +RSOFLZ(L,NY,NX)
      ZCL2(L,NY,NX)=ZCL2(L,NY,NX)+TCLFLS(L,NY,NX)+RCLFXS(L,NY,NX) &
      +RCLFLZ(L,NY,NX)
      ZCO32(L,NY,NX)=ZCO32(L,NY,NX)+TC3FLS(L,NY,NX)+RC3FXS(L,NY,NX) &
      +RC3FLZ(L,NY,NX)
      ZHCO32(L,NY,NX)=ZHCO32(L,NY,NX)+THCFLS(L,NY,NX)+RHCFXS(L,NY,NX) &
      +RHCFLZ(L,NY,NX)
      ZAL12(L,NY,NX)=ZAL12(L,NY,NX)+TAL1FS(L,NY,NX)+RAL1XS(L,NY,NX) &
      +RAL1FZ(L,NY,NX)
      ZAL22(L,NY,NX)=ZAL22(L,NY,NX)+TAL2FS(L,NY,NX)+RAL2XS(L,NY,NX) &
      +RAL2FZ(L,NY,NX)
      ZAL32(L,NY,NX)=ZAL32(L,NY,NX)+TAL3FS(L,NY,NX)+RAL3XS(L,NY,NX) &
      +RAL3FZ(L,NY,NX)
      ZAL42(L,NY,NX)=ZAL42(L,NY,NX)+TAL4FS(L,NY,NX)+RAL4XS(L,NY,NX) &
      +RAL4FZ(L,NY,NX)
      ZALS2(L,NY,NX)=ZALS2(L,NY,NX)+TALSFS(L,NY,NX)+RALSXS(L,NY,NX) &
      +RALSFZ(L,NY,NX)
      ZFE12(L,NY,NX)=ZFE12(L,NY,NX)+TFE1FS(L,NY,NX)+RFE1XS(L,NY,NX) &
      +RFE1FZ(L,NY,NX)
      ZFE22(L,NY,NX)=ZFE22(L,NY,NX)+TFE2FS(L,NY,NX)+RFE2XS(L,NY,NX) &
      +RFE2FZ(L,NY,NX)
      ZFE32(L,NY,NX)=ZFE32(L,NY,NX)+TFE3FS(L,NY,NX)+RFE3XS(L,NY,NX) &
      +RFE3FZ(L,NY,NX)
      ZFE42(L,NY,NX)=ZFE42(L,NY,NX)+TFE4FS(L,NY,NX)+RFE4XS(L,NY,NX) &
      +RFE4FZ(L,NY,NX)
      ZFES2(L,NY,NX)=ZFES2(L,NY,NX)+TFESFS(L,NY,NX)+RFESXS(L,NY,NX) &
      +RFESFZ(L,NY,NX)
      ZCAO2(L,NY,NX)=ZCAO2(L,NY,NX)+TCAOFS(L,NY,NX)+RCAOXS(L,NY,NX) &
      +RCAOFZ(L,NY,NX)
      ZCAC2(L,NY,NX)=ZCAC2(L,NY,NX)+TCACFS(L,NY,NX)+RCACXS(L,NY,NX) &
      +RCACFZ(L,NY,NX)
      ZCAH2(L,NY,NX)=ZCAH2(L,NY,NX)+TCAHFS(L,NY,NX)+RCAHXS(L,NY,NX) &
      +RCAHFZ(L,NY,NX)
      ZCAS2(L,NY,NX)=ZCAS2(L,NY,NX)+TCASFS(L,NY,NX)+RCASXS(L,NY,NX) &
      +RCASFZ(L,NY,NX)
      ZMGO2(L,NY,NX)=ZMGO2(L,NY,NX)+TMGOFS(L,NY,NX)+RMGOXS(L,NY,NX) &
      +RMGOFZ(L,NY,NX)
      ZMGC2(L,NY,NX)=ZMGC2(L,NY,NX)+TMGCFS(L,NY,NX)+RMGCXS(L,NY,NX) &
      +RMGCFZ(L,NY,NX)
      ZMGH2(L,NY,NX)=ZMGH2(L,NY,NX)+TMGHFS(L,NY,NX)+RMGHXS(L,NY,NX) &
      +RMGHFZ(L,NY,NX)
      ZMGS2(L,NY,NX)=ZMGS2(L,NY,NX)+TMGSFS(L,NY,NX)+RMGSXS(L,NY,NX) &
      +RMGSFZ(L,NY,NX)
      ZNAC2(L,NY,NX)=ZNAC2(L,NY,NX)+TNACFS(L,NY,NX)+RNACXS(L,NY,NX) &
      +RNACFZ(L,NY,NX)
      ZNAS2(L,NY,NX)=ZNAS2(L,NY,NX)+TNASFS(L,NY,NX)+RNASXS(L,NY,NX) &
      +RNASFZ(L,NY,NX)
      ZKAS2(L,NY,NX)=ZKAS2(L,NY,NX)+TKASFS(L,NY,NX)+RKASXS(L,NY,NX) &
      +RKASFZ(L,NY,NX)
      H0PO42(L,NY,NX)=H0PO42(L,NY,NX)+TH0PFS(L,NY,NX)+RH0PXS(L,NY,NX) &
      +RH0PFZ(L,NY,NX)
      H3PO42(L,NY,NX)=H3PO42(L,NY,NX)+TH3PFS(L,NY,NX)+RH3PXS(L,NY,NX) &
      +RH3PFZ(L,NY,NX)
      ZFE1P2(L,NY,NX)=ZFE1P2(L,NY,NX)+TF1PFS(L,NY,NX)+RF1PXS(L,NY,NX) &
      +RF1PFZ(L,NY,NX)
      ZFE2P2(L,NY,NX)=ZFE2P2(L,NY,NX)+TF2PFS(L,NY,NX)+RF2PXS(L,NY,NX) &
      +RF2PFZ(L,NY,NX)
      ZCA0P2(L,NY,NX)=ZCA0P2(L,NY,NX)+TC0PFS(L,NY,NX)+RC0PXS(L,NY,NX) &
      +RC0PFZ(L,NY,NX)
      ZCA1P2(L,NY,NX)=ZCA1P2(L,NY,NX)+TC1PFS(L,NY,NX)+RC1PXS(L,NY,NX) &
      +RC1PFZ(L,NY,NX)
      ZCA2P2(L,NY,NX)=ZCA2P2(L,NY,NX)+TC2PFS(L,NY,NX)+RC2PXS(L,NY,NX) &
      +RC2PFZ(L,NY,NX)
      ZMG1P2(L,NY,NX)=ZMG1P2(L,NY,NX)+TM1PFS(L,NY,NX)+RM1PXS(L,NY,NX) &
      +RM1PFZ(L,NY,NX)
      H0POB2(L,NY,NX)=H0POB2(L,NY,NX)+TH0BFB(L,NY,NX)+RH0BXB(L,NY,NX) &
      +RH0BBZ(L,NY,NX)
      H3POB2(L,NY,NX)=H3POB2(L,NY,NX)+TH3BFB(L,NY,NX)+RH3BXB(L,NY,NX) &
      +RH3BBZ(L,NY,NX)
      ZF1PB2(L,NY,NX)=ZF1PB2(L,NY,NX)+TF1BFB(L,NY,NX)+RF1BXB(L,NY,NX) &
      +RF1BBZ(L,NY,NX)
      ZF2PB2(L,NY,NX)=ZF2PB2(L,NY,NX)+TF2BFB(L,NY,NX)+RF2BXB(L,NY,NX) &
      +RF2BBZ(L,NY,NX)
      ZC0PB2(L,NY,NX)=ZC0PB2(L,NY,NX)+TC0BFB(L,NY,NX)+RC0BXB(L,NY,NX) &
      +RC0BBZ(L,NY,NX)
      ZC1PB2(L,NY,NX)=ZC1PB2(L,NY,NX)+TC1BFB(L,NY,NX)+RC1BXB(L,NY,NX) &
      +RC1BBZ(L,NY,NX)
      ZC2PB2(L,NY,NX)=ZC2PB2(L,NY,NX)+TC2BFB(L,NY,NX)+RC2BXB(L,NY,NX) &
      +RC2BBZ(L,NY,NX)
      ZM1PB2(L,NY,NX)=ZM1PB2(L,NY,NX)+TM1BFB(L,NY,NX)+RM1BXB(L,NY,NX) &
      +RM1BBZ(L,NY,NX)
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
!     WRITE(*,444)'ZOH2',I,J,M,NX,NY,L,ZOH2(L,NY,NX)
!    2,TOHFLS(L,NY,NX),ROHFXS(L,NY,NX),ROHFLZ(L,NY,NX)
!    3,ROHFLS(3,L-1,NY,NX),ROHFLS(3,L,NY,NX)
!     WRITE(*,444)'ZAL2',I,J,M,NX,NY,L,ZAL2(L,NY,NX)
!    2,TALFLS(L,NY,NX),RALFXS(L,NY,NX),RALFLZ(L,NY,NX)
!    3,RZAL2(L,NY,NX),TRAL(L,NY,NX)
!444   FORMAT(A8,6I4,20E12.4)
!     ENDIF
      ENDIF
9685  CONTINUE
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
  TQRAL(NY,NX)=0.0
  TQRFE(NY,NX)=0.0
  TQRHY(NY,NX)=0.0
  TQRCA(NY,NX)=0.0
  TQRMG(NY,NX)=0.0
  TQRNA(NY,NX)=0.0
  TQRKA(NY,NX)=0.0
  TQROH(NY,NX)=0.0
  TQRSO(NY,NX)=0.0
  TQRCL(NY,NX)=0.0
  TQRC3(NY,NX)=0.0
  TQRHC(NY,NX)=0.0
  TQRAL1(NY,NX)=0.0
  TQRAL2(NY,NX)=0.0
  TQRAL3(NY,NX)=0.0
  TQRAL4(NY,NX)=0.0
  TQRALS(NY,NX)=0.0
  TQRFE1(NY,NX)=0.0
  TQRFE2(NY,NX)=0.0
  TQRFE3(NY,NX)=0.0
  TQRFE4(NY,NX)=0.0
  TQRFES(NY,NX)=0.0
  TQRCAO(NY,NX)=0.0
  TQRCAC(NY,NX)=0.0
  TQRCAH(NY,NX)=0.0
  TQRCAS(NY,NX)=0.0
  TQRMGO(NY,NX)=0.0
  TQRMGC(NY,NX)=0.0
  TQRMGH(NY,NX)=0.0
  TQRMGS(NY,NX)=0.0
  TQRNAC(NY,NX)=0.0
  TQRNAS(NY,NX)=0.0
  TQRKAS(NY,NX)=0.0
  TQRH0P(NY,NX)=0.0
  TQRH3P(NY,NX)=0.0
  TQRF1P(NY,NX)=0.0
  TQRF2P(NY,NX)=0.0
  TQRC0P(NY,NX)=0.0
  TQRC1P(NY,NX)=0.0
  TQRC2P(NY,NX)=0.0
  TQRM1P(NY,NX)=0.0
  TQSAL(NY,NX)=0.0
  TQSFE(NY,NX)=0.0
  TQSHY(NY,NX)=0.0
  TQSCA(NY,NX)=0.0
  TQSMG(NY,NX)=0.0
  TQSNA(NY,NX)=0.0
  TQSKA(NY,NX)=0.0
  TQSOH(NY,NX)=0.0
  TQSSO(NY,NX)=0.0
  TQSCL(NY,NX)=0.0
  TQSC3(NY,NX)=0.0
  TQSHC(NY,NX)=0.0
  TQSAL1(NY,NX)=0.0
  TQSAL2(NY,NX)=0.0
  TQSAL3(NY,NX)=0.0
  TQSAL4(NY,NX)=0.0
  TQSALS(NY,NX)=0.0
  TQSFE1(NY,NX)=0.0
  TQSFE2(NY,NX)=0.0
  TQSFE3(NY,NX)=0.0
  TQSFE4(NY,NX)=0.0
  TQSFES(NY,NX)=0.0
  TQSCAO(NY,NX)=0.0
  TQSCAC(NY,NX)=0.0
  TQSCAH(NY,NX)=0.0
  TQSCAS(NY,NX)=0.0
  TQSMGO(NY,NX)=0.0
  TQSMGC(NY,NX)=0.0
  TQSMGH(NY,NX)=0.0
  TQSMGS(NY,NX)=0.0
  TQSNAC(NY,NX)=0.0
  TQSNAS(NY,NX)=0.0
  TQSKAS(NY,NX)=0.0
  TQSH0P(NY,NX)=0.0
  TQSH3P(NY,NX)=0.0
  TQSF1P(NY,NX)=0.0
  TQSF2P(NY,NX)=0.0
  TQSC0P(NY,NX)=0.0
  TQSC1P(NY,NX)=0.0
  TQSC2P(NY,NX)=0.0
  TQSM1P(NY,NX)=0.0
  end subroutine InitFluxAccumlatorsInRunoff
!------------------------------------------------------------------------------------------

      subroutine InitFluxAccumulatorsInSoil(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
      integer :: L

  DO 9885 L=NU(NY,NX),NL(NY,NX)
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
9885  CONTINUE
      end subroutine InitFluxAccumulatorsInSoil

end module TrnsfrsMod
