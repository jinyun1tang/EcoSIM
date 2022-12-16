module TrnsfrsDataMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use SoilPropertyDataType
  use TracerIDMod
  use AqueChemDatatype
implicit none
  public
  CHARACTER(LEN=*),private,PARAMETER :: MOD_FILENAME=__FILE__

  real(r8) :: RFLFE,RFLHY,RFLCA,RFLMG,RFLNA,RFLKA,RFLOH,RFLSO
  real(r8) :: RFLCL,RFLC3,RFLHC,RFLAL1,RFLAL2,RFLAL3,RFLAL4,RFLAL
  real(r8) :: RFLALS,RFLFE1,RFLFE2,RFLFE3,RFLFE4,RFLFES,RFLCAO
  real(r8) :: RFLCAC,RFLCAH,RFLCAS,RFLMGO,RFLMGC,RFLMGH,RFLMGS
  real(r8) :: RFLNAC,RFLNAS,RFLKAS,RFLH0P,RFLH3P,RFLF1P,RFLF2P
  real(r8) :: RFLC0P,RFLC1P,RFLC2P,RFLM1P,RFLH0B,RFLH3B,RFLF1B
  real(r8) :: RFLF2B,RFLC0B,RFLC1B,RFLC2B,RFLM1B
  real(r8) :: DFVCA,DFVMG,DFVNA,DFVKA,DFVOH,DFVSO,DFVCL,DFVC3
  real(r8) :: DFVHC,DFVAL1,DFVAL2,DFVAL3,DFVAL4,DFVALS,DFVFE1
  real(r8) :: DFVFE2,DFVFE3,DFVFE4,DFVFES,DFVCAO,DFVCAC,DFVCAH
  real(r8) :: DFVCAS,DFVMGO,DFVMGC,DFVMGH,DFVMGS,DFVNAC,DFVNAS
  real(r8) :: DFVKAS,DFVH0P,DFVH3P,DFVF1P,DFVF2P,DFVC0P,DFVC1P
  real(r8) :: DFVC2P,DFVM1P,DFVH0B,DFVH3B,DFVF1B,DFVF2B,DFVC0B
  real(r8) :: DFVC1B,DFVC2B,DFVM1B,DFVF22,DFVAL,DFVFE,DFVHY

  real(r8), PARAMETER :: VFLWX=0.5_r8
  real(r8),allocatable ::  ALSGL2(:,:,:)                      !
  real(r8),allocatable ::  FESGL2(:,:,:)                      !
  real(r8),allocatable ::  HYSGL2(:,:,:)                      !
  real(r8),allocatable ::  CASGL2(:,:,:)                      !
  real(r8),allocatable ::  GMSGL2(:,:,:)                      !
  real(r8),allocatable ::  ANSGL2(:,:,:)                      !
  real(r8),allocatable ::  AKSGL2(:,:,:)                      !
  real(r8),allocatable ::  OHSGL2(:,:,:)                      !
  real(r8),allocatable ::  C3SGL2(:,:,:)                      !
  real(r8),allocatable ::  HCSGL2(:,:,:)                      !
  real(r8),allocatable ::  SOSGL2(:,:,:)                      !
  real(r8),allocatable ::  CLSXL2(:,:,:)                      !
  real(r8),allocatable ::  trcsa_RQR(:,:,:,:,:)                     !
  real(r8),allocatable ::  trcsa_RQ(:,:,:,:)                       !
  real(r8),allocatable ::  trcsa_RQR0(:,:,:)                        !

  real(r8),allocatable ::  trcsa_solml2(:,:,:,:)              !
  real(r8),allocatable ::  ZFE2(:,:,:)                        !
  real(r8),allocatable ::  ZHCO32(:,:,:)                      !
  real(r8),allocatable ::  trcsa_solml2R(:,:,:,:)                       !
  real(r8),allocatable ::  trcsa_RBLS(:,:,:,:)                      !

  real(r8),allocatable ::  trcsa_TBLS(:,:,:,:)

  real(r8),allocatable ::  POSGL2(:,:,:)                      !
  real(r8),allocatable ::  trcsa_RFLS(:,:,:,:,:)                    !
  real(r8),allocatable ::  RALFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RFEFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RHYFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCAFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RMGFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RNAFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RKAFHS(:,:,:,:)                    !
  real(r8),allocatable ::  ROHFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RC3FHS(:,:,:,:)                    !
  real(r8),allocatable ::  RHCFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RSOFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCLFHS(:,:,:,:)                    !
  real(r8),allocatable ::  RAL1HS(:,:,:,:)                    !
  real(r8),allocatable ::  RAL2HS(:,:,:,:)                    !
  real(r8),allocatable ::  RAL3HS(:,:,:,:)                    !
  real(r8),allocatable ::  RAL4HS(:,:,:,:)                    !
  real(r8),allocatable ::  RALSHS(:,:,:,:)                    !
  real(r8),allocatable ::  RFE1HS(:,:,:,:)                    !
  real(r8),allocatable ::  RFE2HS(:,:,:,:)                    !
  real(r8),allocatable ::  RFE3HS(:,:,:,:)                    !
  real(r8),allocatable ::  RFE4HS(:,:,:,:)                    !
  real(r8),allocatable ::  RFESHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCAOHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCACHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCAHHS(:,:,:,:)                    !
  real(r8),allocatable ::  RCASHS(:,:,:,:)                    !
  real(r8),allocatable ::  RMGOHS(:,:,:,:)                    !
  real(r8),allocatable ::  RMGCHS(:,:,:,:)                    !
  real(r8),allocatable ::  RMGHHS(:,:,:,:)                    !
  real(r8),allocatable ::  RMGSHS(:,:,:,:)                    !
  real(r8),allocatable ::  RNACHS(:,:,:,:)                    !
  real(r8),allocatable ::  RNASHS(:,:,:,:)                    !
  real(r8),allocatable ::  RKASHS(:,:,:,:)                    !
  real(r8),allocatable ::  RH0PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RH3PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RF1PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RF2PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RC0PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RC1PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RC2PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RM1PHS(:,:,:,:)                    !
  real(r8),allocatable ::  RH0BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RH3BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RF1BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RF2BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RC0BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RC1BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RC2BHB(:,:,:,:)                    !
  real(r8),allocatable ::  RM1BHB(:,:,:,:)                    !
  real(r8),allocatable ::  trcsa_soHml2(:,:,:,:)                       !
  real(r8),allocatable ::  trcsa_RFXS(:,:,:,:)                      !
  real(r8), allocatable ::  trcs_solsml2(:,:,:,:)               ! snowpack salt dissolved tracers
  real(r8),allocatable ::  trcsa_RFL0(:,:,:)                        !
  real(r8),allocatable ::  RALFL1(:,:)                        !
  real(r8),allocatable ::  RFEFL1(:,:)                        !
  real(r8),allocatable ::  RHYFL1(:,:)                        !
  real(r8),allocatable ::  RCAFL1(:,:)                        !
  real(r8),allocatable ::  RMGFL1(:,:)                        !
  real(r8),allocatable ::  RNAFL1(:,:)                        !
  real(r8),allocatable ::  RKAFL1(:,:)                        !
  real(r8),allocatable ::  ROHFL1(:,:)                        !
  real(r8),allocatable ::  RSOFL1(:,:)                        !
  real(r8),allocatable ::  RCLFL1(:,:)                        !
  real(r8),allocatable ::  RC3FL1(:,:)                        !
  real(r8),allocatable ::  RHCFL1(:,:)                        !
  real(r8),allocatable ::  RAL1F1(:,:)                        !
  real(r8),allocatable ::  RAL2F1(:,:)                        !
  real(r8),allocatable ::  RAL3F1(:,:)                        !
  real(r8),allocatable ::  RAL4F1(:,:)                        !
  real(r8),allocatable ::  RALSF1(:,:)                        !
  real(r8),allocatable ::  RFE1F1(:,:)                        !
  real(r8),allocatable ::  RFE2F1(:,:)                        !
  real(r8),allocatable ::  RFE3F1(:,:)                        !
  real(r8),allocatable ::  RFE4F1(:,:)                        !
  real(r8),allocatable ::  RFESF1(:,:)                        !
  real(r8),allocatable ::  RCAOF1(:,:)                        !
  real(r8),allocatable ::  RCACF1(:,:)                        !
  real(r8),allocatable ::  RCAHF1(:,:)                        !
  real(r8),allocatable ::  RCASF1(:,:)                        !
  real(r8),allocatable ::  RMGOF1(:,:)                        !
  real(r8),allocatable ::  RMGCF1(:,:)                        !
  real(r8),allocatable ::  RMGHF1(:,:)                        !
  real(r8),allocatable ::  RMGSF1(:,:)                        !
  real(r8),allocatable ::  RNACF1(:,:)                        !
  real(r8),allocatable ::  RNASF1(:,:)                        !
  real(r8),allocatable ::  RKASF1(:,:)                        !
  real(r8),allocatable ::  RH0PF1(:,:)                        !
  real(r8),allocatable ::  RH3PF1(:,:)                        !
  real(r8),allocatable ::  RF1PF1(:,:)                        !
  real(r8),allocatable ::  RF2PF1(:,:)                        !
  real(r8),allocatable ::  RC0PF1(:,:)                        !
  real(r8),allocatable ::  RC1PF1(:,:)                        !
  real(r8),allocatable ::  RC2PF1(:,:)                        !
  real(r8),allocatable ::  RM1PF1(:,:)                        !
  real(r8),allocatable ::  RH0BF2(:,:)                        !
  real(r8),allocatable ::  RH3BF2(:,:)                        !
  real(r8),allocatable ::  RF1BF2(:,:)                        !
  real(r8),allocatable ::  RF2BF2(:,:)                        !
  real(r8),allocatable ::  RC0BF2(:,:)                        !
  real(r8),allocatable ::  RC1BF2(:,:)                        !
  real(r8),allocatable ::  RC2BF2(:,:)                        !
  real(r8),allocatable ::  RM1BF2(:,:)                        !
!----------------------------------------------------------------------

contains
  subroutine InitTrnsfrsData

  implicit none
  allocate(ALSGL2(JZ,JY,JX));   ALSGL2=0._r8
  allocate(FESGL2(JZ,JY,JX));   FESGL2=0._r8
  allocate(HYSGL2(JZ,JY,JX));   HYSGL2=0._r8
  allocate(CASGL2(JZ,JY,JX));   CASGL2=0._r8
  allocate(GMSGL2(JZ,JY,JX));   GMSGL2=0._r8
  allocate(ANSGL2(JZ,JY,JX));   ANSGL2=0._r8
  allocate(AKSGL2(JZ,JY,JX));   AKSGL2=0._r8
  allocate(OHSGL2(JZ,JY,JX));   OHSGL2=0._r8
  allocate(C3SGL2(JZ,JY,JX));   C3SGL2=0._r8
  allocate(HCSGL2(JZ,JY,JX));   HCSGL2=0._r8
  allocate(SOSGL2(JZ,JY,JX));   SOSGL2=0._r8
  allocate(CLSXL2(JZ,JY,JX));   CLSXL2=0._r8
  allocate(trcsa_RQR(idsa_beg:idsa_end,2,2,JV,JH));   trcsa_RQR=0._r8
  allocate(trcsa_RQ(idsa_beg:idsab_end,2,JV,JH));     trcsa_RQ=0._r8
  allocate(trcsa_RQR0(idsa_beg:idsa_end,JY,JX));     trcsa_RQR0=0._r8
  allocate(trcsa_RBLS(idsa_beg:idsa_end,JS,JY,JX));   trcsa_RBLS=0._r8
  allocate(trcsa_TBLS(idsa_beg:idsa_end,JS,JY,JX)); trcsa_TBLS=0._r8

  allocate(POSGL2(JZ,JY,JX));   POSGL2=0._r8
  allocate(trcsa_RFLS(idsa_beg:idsab_end,3,0:JD,JV,JH));trcsa_RFLS=0._r8
  allocate(RALFHS(3,JD,JV,JH)); RALFHS=0._r8
  allocate(RFEFHS(3,JD,JV,JH)); RFEFHS=0._r8
  allocate(RHYFHS(3,JD,JV,JH)); RHYFHS=0._r8
  allocate(RCAFHS(3,JD,JV,JH)); RCAFHS=0._r8
  allocate(RMGFHS(3,JD,JV,JH)); RMGFHS=0._r8
  allocate(RNAFHS(3,JD,JV,JH)); RNAFHS=0._r8
  allocate(RKAFHS(3,JD,JV,JH)); RKAFHS=0._r8
  allocate(ROHFHS(3,JD,JV,JH)); ROHFHS=0._r8
  allocate(RC3FHS(3,JD,JV,JH)); RC3FHS=0._r8
  allocate(RHCFHS(3,JD,JV,JH)); RHCFHS=0._r8
  allocate(RSOFHS(3,JD,JV,JH)); RSOFHS=0._r8
  allocate(RCLFHS(3,JD,JV,JH)); RCLFHS=0._r8
  allocate(RAL1HS(3,JD,JV,JH)); RAL1HS=0._r8
  allocate(RAL2HS(3,JD,JV,JH)); RAL2HS=0._r8
  allocate(RAL3HS(3,JD,JV,JH)); RAL3HS=0._r8
  allocate(RAL4HS(3,JD,JV,JH)); RAL4HS=0._r8
  allocate(RALSHS(3,JD,JV,JH)); RALSHS=0._r8
  allocate(RFE1HS(3,JD,JV,JH)); RFE1HS=0._r8
  allocate(RFE2HS(3,JD,JV,JH)); RFE2HS=0._r8
  allocate(RFE3HS(3,JD,JV,JH)); RFE3HS=0._r8
  allocate(RFE4HS(3,JD,JV,JH)); RFE4HS=0._r8
  allocate(RFESHS(3,JD,JV,JH)); RFESHS=0._r8
  allocate(RCAOHS(3,JD,JV,JH)); RCAOHS=0._r8
  allocate(RCACHS(3,JD,JV,JH)); RCACHS=0._r8
  allocate(RCAHHS(3,JD,JV,JH)); RCAHHS=0._r8
  allocate(RCASHS(3,JD,JV,JH)); RCASHS=0._r8
  allocate(RMGOHS(3,JD,JV,JH)); RMGOHS=0._r8
  allocate(RMGCHS(3,JD,JV,JH)); RMGCHS=0._r8
  allocate(RMGHHS(3,JD,JV,JH)); RMGHHS=0._r8
  allocate(RMGSHS(3,JD,JV,JH)); RMGSHS=0._r8
  allocate(RNACHS(3,JD,JV,JH)); RNACHS=0._r8
  allocate(RNASHS(3,JD,JV,JH)); RNASHS=0._r8
  allocate(RKASHS(3,JD,JV,JH)); RKASHS=0._r8
  allocate(RH0PHS(3,JD,JV,JH)); RH0PHS=0._r8
  allocate(RH3PHS(3,JD,JV,JH)); RH3PHS=0._r8
  allocate(RF1PHS(3,JD,JV,JH)); RF1PHS=0._r8
  allocate(RF2PHS(3,JD,JV,JH)); RF2PHS=0._r8
  allocate(RC0PHS(3,JD,JV,JH)); RC0PHS=0._r8
  allocate(RC1PHS(3,JD,JV,JH)); RC1PHS=0._r8
  allocate(RC2PHS(3,JD,JV,JH)); RC2PHS=0._r8
  allocate(RM1PHS(3,JD,JV,JH)); RM1PHS=0._r8
  allocate(RH0BHB(3,JD,JV,JH)); RH0BHB=0._r8
  allocate(RH3BHB(3,JD,JV,JH)); RH3BHB=0._r8
  allocate(RF1BHB(3,JD,JV,JH)); RF1BHB=0._r8
  allocate(RF2BHB(3,JD,JV,JH)); RF2BHB=0._r8
  allocate(RC0BHB(3,JD,JV,JH)); RC0BHB=0._r8
  allocate(RC1BHB(3,JD,JV,JH)); RC1BHB=0._r8
  allocate(RC2BHB(3,JD,JV,JH)); RC2BHB=0._r8
  allocate(RM1BHB(3,JD,JV,JH)); RM1BHB=0._r8
  allocate(trcsa_soHml2(idsa_beg:idsab_end,JZ,JY,JX));    trcsa_soHml2=0._r8
  allocate(trcsa_RFXS(idsa_beg:idsab_end,JZ,JY,JX));   trcsa_RFXS=0._r8

  allocate(trcs_solsml2(idsa_beg:idsa_end,JS,JY,JX)); trcs_solsml2=0._r8

  allocate(trcsa_RFL0(idsa_beg:idsa_end,JY,JX));      trcsa_RFL0=0._r8
  allocate(RALFL1(JY,JX));      RALFL1=0._r8
  allocate(RFEFL1(JY,JX));      RFEFL1=0._r8
  allocate(RHYFL1(JY,JX));      RHYFL1=0._r8
  allocate(RCAFL1(JY,JX));      RCAFL1=0._r8
  allocate(RMGFL1(JY,JX));      RMGFL1=0._r8
  allocate(RNAFL1(JY,JX));      RNAFL1=0._r8
  allocate(RKAFL1(JY,JX));      RKAFL1=0._r8
  allocate(ROHFL1(JY,JX));      ROHFL1=0._r8
  allocate(RSOFL1(JY,JX));      RSOFL1=0._r8
  allocate(RCLFL1(JY,JX));      RCLFL1=0._r8
  allocate(RC3FL1(JY,JX));      RC3FL1=0._r8
  allocate(RHCFL1(JY,JX));      RHCFL1=0._r8
  allocate(RAL1F1(JY,JX));      RAL1F1=0._r8
  allocate(RAL2F1(JY,JX));      RAL2F1=0._r8
  allocate(RAL3F1(JY,JX));      RAL3F1=0._r8
  allocate(RAL4F1(JY,JX));      RAL4F1=0._r8
  allocate(RALSF1(JY,JX));      RALSF1=0._r8
  allocate(RFE1F1(JY,JX));      RFE1F1=0._r8
  allocate(RFE2F1(JY,JX));      RFE2F1=0._r8
  allocate(RFE3F1(JY,JX));      RFE3F1=0._r8
  allocate(RFE4F1(JY,JX));      RFE4F1=0._r8
  allocate(RFESF1(JY,JX));      RFESF1=0._r8
  allocate(RCAOF1(JY,JX));      RCAOF1=0._r8
  allocate(RCACF1(JY,JX));      RCACF1=0._r8
  allocate(RCAHF1(JY,JX));      RCAHF1=0._r8
  allocate(RCASF1(JY,JX));      RCASF1=0._r8
  allocate(RMGOF1(JY,JX));      RMGOF1=0._r8
  allocate(RMGCF1(JY,JX));      RMGCF1=0._r8
  allocate(RMGHF1(JY,JX));      RMGHF1=0._r8
  allocate(RMGSF1(JY,JX));      RMGSF1=0._r8
  allocate(RNACF1(JY,JX));      RNACF1=0._r8
  allocate(RNASF1(JY,JX));      RNASF1=0._r8
  allocate(RKASF1(JY,JX));      RKASF1=0._r8
  allocate(RH0PF1(JY,JX));      RH0PF1=0._r8
  allocate(RH3PF1(JY,JX));      RH3PF1=0._r8
  allocate(RF1PF1(JY,JX));      RF1PF1=0._r8
  allocate(RF2PF1(JY,JX));      RF2PF1=0._r8
  allocate(RC0PF1(JY,JX));      RC0PF1=0._r8
  allocate(RC1PF1(JY,JX));      RC1PF1=0._r8
  allocate(RC2PF1(JY,JX));      RC2PF1=0._r8
  allocate(RM1PF1(JY,JX));      RM1PF1=0._r8
  allocate(RH0BF2(JY,JX));      RH0BF2=0._r8
  allocate(RH3BF2(JY,JX));      RH3BF2=0._r8
  allocate(RF1BF2(JY,JX));      RF1BF2=0._r8
  allocate(RF2BF2(JY,JX));      RF2BF2=0._r8
  allocate(RC0BF2(JY,JX));      RC0BF2=0._r8
  allocate(RC1BF2(JY,JX));      RC1BF2=0._r8
  allocate(RC2BF2(JY,JX));      RC2BF2=0._r8
  allocate(RM1BF2(JY,JX));      RM1BF2=0._r8
  end subroutine InitTrnsfrsData

!----------------------------------------------------------------------
  subroutine DestructTrnsfrsData
  use abortutils, only : destroy
  implicit none

  call destroy(trcsa_RFLS)
  call destroy(trcsa_solml2)
  call destroy(trcsa_solml2R)

  call destroy(ALSGL2)
  call destroy(FESGL2)
  call destroy(HYSGL2)
  call destroy(CASGL2)
  call destroy(GMSGL2)
  call destroy(ANSGL2)
  call destroy(AKSGL2)
  call destroy(OHSGL2)
  call destroy(C3SGL2)
  call destroy(HCSGL2)
  call destroy(SOSGL2)
  call destroy(CLSXL2)
  call destroy(trcsa_RQ)

  call destroy(POSGL2)
  call destroy(RALFHS)
  call destroy(RFEFHS)
  call destroy(RHYFHS)
  call destroy(RCAFHS)
  call destroy(RMGFHS)
  call destroy(RNAFHS)
  call destroy(RKAFHS)
  call destroy(ROHFHS)
  call destroy(RC3FHS)
  call destroy(RHCFHS)
  call destroy(RSOFHS)
  call destroy(RCLFHS)
  call destroy(RAL1HS)
  call destroy(RAL2HS)
  call destroy(RAL3HS)
  call destroy(RAL4HS)
  call destroy(RALSHS)
  call destroy(RFE1HS)
  call destroy(RFE2HS)
  call destroy(RFE3HS)
  call destroy(RFE4HS)
  call destroy(RFESHS)
  call destroy(RCAOHS)
  call destroy(RCACHS)
  call destroy(RCAHHS)
  call destroy(RCASHS)
  call destroy(RMGOHS)
  call destroy(RMGCHS)
  call destroy(RMGHHS)
  call destroy(RMGSHS)
  call destroy(RNACHS)
  call destroy(RNASHS)
  call destroy(RKASHS)
  call destroy(RH0PHS)
  call destroy(RH3PHS)
  call destroy(RF1PHS)
  call destroy(RF2PHS)
  call destroy(RC0PHS)
  call destroy(RC1PHS)
  call destroy(RC2PHS)
  call destroy(RM1PHS)
  call destroy(RH0BHB)
  call destroy(RH3BHB)
  call destroy(RF1BHB)
  call destroy(RF2BHB)
  call destroy(RC0BHB)
  call destroy(RC1BHB)
  call destroy(RC2BHB)
  call destroy(RM1BHB)
  call destroy(RALFL1)
  call destroy(RFEFL1)
  call destroy(RHYFL1)
  call destroy(RCAFL1)
  call destroy(RMGFL1)
  call destroy(RNAFL1)
  call destroy(RKAFL1)
  call destroy(ROHFL1)
  call destroy(RSOFL1)
  call destroy(RCLFL1)
  call destroy(RC3FL1)
  call destroy(RHCFL1)
  call destroy(RAL1F1)
  call destroy(RAL2F1)
  call destroy(RAL3F1)
  call destroy(RAL4F1)
  call destroy(RALSF1)
  call destroy(RFE1F1)
  call destroy(RFE2F1)
  call destroy(RFE3F1)
  call destroy(RFE4F1)
  call destroy(RFESF1)
  call destroy(RCAOF1)
  call destroy(RCACF1)
  call destroy(RCAHF1)
  call destroy(RCASF1)
  call destroy(RMGOF1)
  call destroy(RMGCF1)
  call destroy(RMGHF1)
  call destroy(RMGSF1)
  call destroy(RNACF1)
  call destroy(RNASF1)
  call destroy(RKASF1)
  call destroy(RH0PF1)
  call destroy(RH3PF1)
  call destroy(RF1PF1)
  call destroy(RF2PF1)
  call destroy(RC0PF1)
  call destroy(RC1PF1)
  call destroy(RC2PF1)
  call destroy(RM1PF1)
  call destroy(RH0BF2)
  call destroy(RH3BF2)
  call destroy(RF1BF2)
  call destroy(RF2BF2)
  call destroy(RC0BF2)
  call destroy(RC1BF2)
  call destroy(RC2BF2)
  call destroy(RM1BF2)
  call destroy(trcsa_RQR)
  call destroy(trcsa_RFXS)
  end subroutine DestructTrnsfrsData

end module TrnsfrsDataMod
