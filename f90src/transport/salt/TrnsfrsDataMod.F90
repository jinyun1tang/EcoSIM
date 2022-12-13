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
  real(r8),allocatable ::  RQRAL(:,:,:,:)                     !
  real(r8),allocatable ::  RQRFE(:,:,:,:)                     !
  real(r8),allocatable ::  RQRHY(:,:,:,:)                     !
  real(r8),allocatable ::  RQRCA(:,:,:,:)                     !
  real(r8),allocatable ::  RQRMG(:,:,:,:)                     !
  real(r8),allocatable ::  RQRNA(:,:,:,:)                     !
  real(r8),allocatable ::  RQRKA(:,:,:,:)                     !
  real(r8),allocatable ::  RQROH(:,:,:,:)                     !
  real(r8),allocatable ::  RQRSO(:,:,:,:)                     !
  real(r8),allocatable ::  RQRCL(:,:,:,:)                     !
  real(r8),allocatable ::  RQRC3(:,:,:,:)                     !
  real(r8),allocatable ::  RQRHC(:,:,:,:)                     !
  real(r8),allocatable ::  RQRAL1(:,:,:,:)                    !
  real(r8),allocatable ::  RQRAL2(:,:,:,:)                    !
  real(r8),allocatable ::  RQRAL3(:,:,:,:)                    !
  real(r8),allocatable ::  RQRAL4(:,:,:,:)                    !
  real(r8),allocatable ::  RQRALS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRFE1(:,:,:,:)                    !
  real(r8),allocatable ::  RQRFE2(:,:,:,:)                    !
  real(r8),allocatable ::  RQRFE3(:,:,:,:)                    !
  real(r8),allocatable ::  RQRFE4(:,:,:,:)                    !
  real(r8),allocatable ::  RQRFES(:,:,:,:)                    !
  real(r8),allocatable ::  RQRCAO(:,:,:,:)                    !
  real(r8),allocatable ::  RQRCAC(:,:,:,:)                    !
  real(r8),allocatable ::  RQRCAH(:,:,:,:)                    !
  real(r8),allocatable ::  RQRCAS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRMGO(:,:,:,:)                    !
  real(r8),allocatable ::  RQRMGC(:,:,:,:)                    !
  real(r8),allocatable ::  RQRMGH(:,:,:,:)                    !
  real(r8),allocatable ::  RQRMGS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNAC(:,:,:,:)                    !
  real(r8),allocatable ::  RQRNAS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRKAS(:,:,:,:)                    !
  real(r8),allocatable ::  RQRH0P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRH3P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRF1P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRF2P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRC0P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRC1P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRC2P(:,:,:,:)                    !
  real(r8),allocatable ::  RQRM1P(:,:,:,:)                    !
  real(r8),allocatable ::  trcsa_RQ(:,:,:,:)                       !
  real(r8),allocatable ::  RQRAL0(:,:)                        !
  real(r8),allocatable ::  RQRFE0(:,:)                        !
  real(r8),allocatable ::  RQRHY0(:,:)                        !
  real(r8),allocatable ::  RQRCA0(:,:)                        !
  real(r8),allocatable ::  RQRMG0(:,:)                        !
  real(r8),allocatable ::  RQRNA0(:,:)                        !
  real(r8),allocatable ::  RQRKA0(:,:)                        !
  real(r8),allocatable ::  RQROH0(:,:)                        !
  real(r8),allocatable ::  RQRSO0(:,:)                        !
  real(r8),allocatable ::  RQRCL0(:,:)                        !
  real(r8),allocatable ::  RQRC30(:,:)                        !
  real(r8),allocatable ::  RQRHC0(:,:)                        !
  real(r8),allocatable ::  RQRAL10(:,:)                       !
  real(r8),allocatable ::  RQRAL20(:,:)                       !
  real(r8),allocatable ::  RQRAL30(:,:)                       !
  real(r8),allocatable ::  RQRAL40(:,:)                       !
  real(r8),allocatable ::  RQRALS0(:,:)                       !
  real(r8),allocatable ::  RQRFE10(:,:)                       !
  real(r8),allocatable ::  RQRFE20(:,:)                       !
  real(r8),allocatable ::  RQRFE30(:,:)                       !
  real(r8),allocatable ::  RQRFE40(:,:)                       !
  real(r8),allocatable ::  RQRFES0(:,:)                       !
  real(r8),allocatable ::  RQRCAO0(:,:)                       !
  real(r8),allocatable ::  RQRCAC0(:,:)                       !
  real(r8),allocatable ::  RQRCAH0(:,:)                       !
  real(r8),allocatable ::  RQRCAS0(:,:)                       !
  real(r8),allocatable ::  RQRMGO0(:,:)                       !
  real(r8),allocatable ::  RQRMGC0(:,:)                       !
  real(r8),allocatable ::  RQRMGH0(:,:)                       !
  real(r8),allocatable ::  RQRMGS0(:,:)                       !
  real(r8),allocatable ::  RQRNAC0(:,:)                       !
  real(r8),allocatable ::  RQRNAS0(:,:)                       !
  real(r8),allocatable ::  RQRKAS0(:,:)                       !
  real(r8),allocatable ::  RQRH0P0(:,:)                       !
  real(r8),allocatable ::  RQRH3P0(:,:)                       !
  real(r8),allocatable ::  RQRF1P0(:,:)                       !
  real(r8),allocatable ::  RQRF2P0(:,:)                       !
  real(r8),allocatable ::  RQRC0P0(:,:)                       !
  real(r8),allocatable ::  RQRC1P0(:,:)                       !
  real(r8),allocatable ::  RQRC2P0(:,:)                       !
  real(r8),allocatable ::  RQRM1P0(:,:)                       !

  real(r8),allocatable ::  trcsa_solml2(:,:,:,:)              !
  real(r8),allocatable ::  ZFE2(:,:,:)                        !
  real(r8),allocatable ::  ZHCO32(:,:,:)                      !
  real(r8),allocatable ::  trcsa_solml2R(:,:,:,:)                       !
  real(r8),allocatable ::  RALBLS(:,:,:)                      !
  real(r8),allocatable ::  RFEBLS(:,:,:)                      !
  real(r8),allocatable ::  RHYBLS(:,:,:)                      !
  real(r8),allocatable ::  RCABLS(:,:,:)                      !
  real(r8),allocatable ::  RMGBLS(:,:,:)                      !
  real(r8),allocatable ::  RNABLS(:,:,:)                      !
  real(r8),allocatable ::  RKABLS(:,:,:)                      !
  real(r8),allocatable ::  ROHBLS(:,:,:)                      !
  real(r8),allocatable ::  RSOBLS(:,:,:)                      !
  real(r8),allocatable ::  RCLBLS(:,:,:)                      !
  real(r8),allocatable ::  RC3BLS(:,:,:)                      !
  real(r8),allocatable ::  RHCBLS(:,:,:)                      !
  real(r8),allocatable ::  RAL1BS(:,:,:)                      !
  real(r8),allocatable ::  RAL2BS(:,:,:)                      !
  real(r8),allocatable ::  RAL3BS(:,:,:)                      !
  real(r8),allocatable ::  RAL4BS(:,:,:)                      !
  real(r8),allocatable ::  RALSBS(:,:,:)                      !
  real(r8),allocatable ::  RFE1BS(:,:,:)                      !
  real(r8),allocatable ::  RFE2BS(:,:,:)                      !
  real(r8),allocatable ::  RFE3BS(:,:,:)                      !
  real(r8),allocatable ::  RFE4BS(:,:,:)                      !
  real(r8),allocatable ::  RFESBS(:,:,:)                      !
  real(r8),allocatable ::  RCAOBS(:,:,:)                      !
  real(r8),allocatable ::  RCACBS(:,:,:)                      !
  real(r8),allocatable ::  RCAHBS(:,:,:)                      !
  real(r8),allocatable ::  RCASBS(:,:,:)                      !
  real(r8),allocatable ::  RMGOBS(:,:,:)                      !
  real(r8),allocatable ::  RMGCBS(:,:,:)                      !
  real(r8),allocatable ::  RMGHBS(:,:,:)                      !
  real(r8),allocatable ::  RMGSBS(:,:,:)                      !
  real(r8),allocatable ::  RNACBS(:,:,:)                      !
  real(r8),allocatable ::  RNASBS(:,:,:)                      !
  real(r8),allocatable ::  RKASBS(:,:,:)                      !
  real(r8),allocatable ::  RH0PBS(:,:,:)                      !
  real(r8),allocatable ::  RH3PBS(:,:,:)                      !
  real(r8),allocatable ::  RF1PBS(:,:,:)                      !
  real(r8),allocatable ::  RF2PBS(:,:,:)                      !
  real(r8),allocatable ::  RC0PBS(:,:,:)                      !
  real(r8),allocatable ::  RC1PBS(:,:,:)                      !
  real(r8),allocatable ::  RC2PBS(:,:,:)                      !
  real(r8),allocatable ::  RM1PBS(:,:,:)                      !

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
  real(r8),allocatable ::  ZALH2(:,:,:)                       !
  real(r8),allocatable ::  ZFEH2(:,:,:)                       !
  real(r8),allocatable ::  ZHYH2(:,:,:)                       !
  real(r8),allocatable ::  ZCCH2(:,:,:)                       !
  real(r8),allocatable ::  ZMAH2(:,:,:)                       !
  real(r8),allocatable ::  ZNAH2(:,:,:)                       !
  real(r8),allocatable ::  ZKAH2(:,:,:)                       !
  real(r8),allocatable ::  ZOHH2(:,:,:)                       !
  real(r8),allocatable ::  ZSO4H2(:,:,:)                      !
  real(r8),allocatable ::  ZCLH2(:,:,:)                       !
  real(r8),allocatable ::  ZCO3H2(:,:,:)                      !
  real(r8),allocatable ::  ZHCOH2(:,:,:)                      !
  real(r8),allocatable ::  ZAL1H2(:,:,:)                      !
  real(r8),allocatable ::  ZAL2H2(:,:,:)                      !
  real(r8),allocatable ::  ZAL3H2(:,:,:)                      !
  real(r8),allocatable ::  ZAL4H2(:,:,:)                      !
  real(r8),allocatable ::  ZALSH2(:,:,:)                      !
  real(r8),allocatable ::  ZFE1H2(:,:,:)                      !
  real(r8),allocatable ::  ZFE2H2(:,:,:)                      !
  real(r8),allocatable ::  ZFE3H2(:,:,:)                      !
  real(r8),allocatable ::  ZFE4H2(:,:,:)                      !
  real(r8),allocatable ::  ZFESH2(:,:,:)                      !
  real(r8),allocatable ::  ZCAOH2(:,:,:)                      !
  real(r8),allocatable ::  ZCACH2(:,:,:)                      !
  real(r8),allocatable ::  ZCAHH2(:,:,:)                      !
  real(r8),allocatable ::  ZCASH2(:,:,:)                      !
  real(r8),allocatable ::  ZMGOH2(:,:,:)                      !
  real(r8),allocatable ::  ZMGCH2(:,:,:)                      !
  real(r8),allocatable ::  ZMGHH2(:,:,:)                      !
  real(r8),allocatable ::  ZMGSH2(:,:,:)                      !
  real(r8),allocatable ::  ZNACH2(:,:,:)                      !
  real(r8),allocatable ::  ZNASH2(:,:,:)                      !
  real(r8),allocatable ::  ZKASH2(:,:,:)                      !
  real(r8),allocatable ::  H0P4H2(:,:,:)                      !
  real(r8),allocatable ::  H3P4H2(:,:,:)                      !
  real(r8),allocatable ::  ZF1PH2(:,:,:)                      !
  real(r8),allocatable ::  ZF2PH2(:,:,:)                      !
  real(r8),allocatable ::  ZC0PH2(:,:,:)                      !
  real(r8),allocatable ::  ZC1PH2(:,:,:)                      !
  real(r8),allocatable ::  ZC2PH2(:,:,:)                      !
  real(r8),allocatable ::  ZM1PH2(:,:,:)                      !
  real(r8),allocatable ::  H0PBH2(:,:,:)                      !
  real(r8),allocatable ::  H3PBH2(:,:,:)                      !
  real(r8),allocatable ::  ZF1BH2(:,:,:)                      !
  real(r8),allocatable ::  ZF2BH2(:,:,:)                      !
  real(r8),allocatable ::  ZC0BH2(:,:,:)                      !
  real(r8),allocatable ::  ZC1BH2(:,:,:)                      !
  real(r8),allocatable ::  ZC2BH2(:,:,:)                      !
  real(r8),allocatable ::  ZM1BH2(:,:,:)                      !
  real(r8),allocatable ::  RH0PXS(:,:,:)                      !
  real(r8),allocatable ::  RH3PXS(:,:,:)                      !
  real(r8),allocatable :: RF1PXS(:,:,:)                       !
  real(r8),allocatable ::  RF2PXS(:,:,:)                      !
  real(r8),allocatable ::  RC0PXS(:,:,:)                      !
  real(r8),allocatable :: RC1PXS(:,:,:)                       !
  real(r8),allocatable ::  RC2PXS(:,:,:)                      !
  real(r8),allocatable ::  RM1PXS(:,:,:)                      !
  real(r8),allocatable :: RH0BXB(:,:,:)                       !
  real(r8),allocatable ::  RH3BXB(:,:,:)                      !
  real(r8),allocatable :: RF1BXB(:,:,:)                       !
  real(r8),allocatable ::  RF2BXB(:,:,:)                      !
  real(r8),allocatable ::  RC0BXB(:,:,:)                      !
  real(r8),allocatable :: RC1BXB(:,:,:)                       !
  real(r8),allocatable ::  RC2BXB(:,:,:)                      !
  real(r8),allocatable ::  RM1BXB(:,:,:)                      !
  real(r8),allocatable ::  RALFXS(:,:,:)                      !
  real(r8),allocatable :: RFEFXS(:,:,:)                       !
  real(r8),allocatable ::  RHYFXS(:,:,:)                      !
  real(r8),allocatable ::  RCAFXS(:,:,:)                      !
  real(r8),allocatable :: RMGFXS(:,:,:)                       !
  real(r8),allocatable ::  RNAFXS(:,:,:)                      !
  real(r8),allocatable ::  RKAFXS(:,:,:)                      !
  real(r8),allocatable :: ROHFXS(:,:,:)                       !
  real(r8),allocatable ::  RSOFXS(:,:,:)                      !
  real(r8),allocatable ::  RCLFXS(:,:,:)                      !
  real(r8),allocatable :: RC3FXS(:,:,:)                       !
  real(r8),allocatable ::  RHCFXS(:,:,:)                      !
  real(r8),allocatable ::  RAL1XS(:,:,:)                      !
  real(r8),allocatable :: RAL2XS(:,:,:)                       !
  real(r8),allocatable ::  RAL3XS(:,:,:)                      !
  real(r8),allocatable ::  RAL4XS(:,:,:)                      !
  real(r8),allocatable :: RALSXS(:,:,:)                       !
  real(r8),allocatable ::  RFE1XS(:,:,:)                      !
  real(r8),allocatable ::  RFE2XS(:,:,:)                      !
  real(r8),allocatable :: RFE3XS(:,:,:)                       !
  real(r8),allocatable ::  RFE4XS(:,:,:)                      !
  real(r8),allocatable ::  RFESXS(:,:,:)                      !
  real(r8),allocatable :: RCACXS(:,:,:)                       !
  real(r8),allocatable ::  RCAOXS(:,:,:)                      !
  real(r8),allocatable ::  RCAHXS(:,:,:)                      !
  real(r8),allocatable :: RCASXS(:,:,:)                       !
  real(r8),allocatable ::  RMGOXS(:,:,:)                      !
  real(r8),allocatable ::  RMGCXS(:,:,:)                      !
  real(r8),allocatable :: RMGHXS(:,:,:)                       !
  real(r8),allocatable ::  RMGSXS(:,:,:)                      !
  real(r8),allocatable ::  RNACXS(:,:,:)                      !
  real(r8),allocatable :: RNASXS(:,:,:)                       !
  real(r8),allocatable ::  RKASXS(:,:,:)                      !

  real(r8), allocatable ::  trcs_solsml2(:,:,:,:)               ! snowpack salt dissolved tracers


  real(r8),allocatable ::  RALFL0(:,:)                        !
  real(r8),allocatable ::  RFEFL0(:,:)                        !
  real(r8),allocatable ::  RHYFL0(:,:)                        !
  real(r8),allocatable ::  RCAFL0(:,:)                        !
  real(r8),allocatable ::  RMGFL0(:,:)                        !
  real(r8),allocatable ::  RNAFL0(:,:)                        !
  real(r8),allocatable ::  RKAFL0(:,:)                        !
  real(r8),allocatable ::  ROHFL0(:,:)                        !
  real(r8),allocatable ::  RSOFL0(:,:)                        !
  real(r8),allocatable ::  RCLFL0(:,:)                        !
  real(r8),allocatable ::  RC3FL0(:,:)                        !
  real(r8),allocatable ::  RHCFL0(:,:)                        !
  real(r8),allocatable ::  RAL1F0(:,:)                        !
  real(r8),allocatable ::  RAL2F0(:,:)                        !
  real(r8),allocatable ::  RAL3F0(:,:)                        !
  real(r8),allocatable ::  RAL4F0(:,:)                        !
  real(r8),allocatable ::  RALSF0(:,:)                        !
  real(r8),allocatable ::  RFE1F0(:,:)                        !
  real(r8),allocatable ::  RFE2F0(:,:)                        !
  real(r8),allocatable ::  RFE3F0(:,:)                        !
  real(r8),allocatable ::  RFE4F0(:,:)                        !
  real(r8),allocatable ::  RFESF0(:,:)                        !
  real(r8),allocatable ::  RCAOF0(:,:)                        !
  real(r8),allocatable ::  RCACF0(:,:)                        !
  real(r8),allocatable ::  RCAHF0(:,:)                        !
  real(r8),allocatable ::  RCASF0(:,:)                        !
  real(r8),allocatable ::  RMGOF0(:,:)                        !
  real(r8),allocatable ::  RMGCF0(:,:)                        !
  real(r8),allocatable ::  RMGHF0(:,:)                        !
  real(r8),allocatable ::  RMGSF0(:,:)                        !
  real(r8),allocatable ::  RNACF0(:,:)                        !
  real(r8),allocatable ::  RNASF0(:,:)                        !
  real(r8),allocatable ::  RKASF0(:,:)                        !
  real(r8),allocatable ::  RH0PF0(:,:)                        !
  real(r8),allocatable ::  RH3PF0(:,:)                        !
  real(r8),allocatable ::  RF1PF0(:,:)                        !
  real(r8),allocatable ::  RF2PF0(:,:)                        !
  real(r8),allocatable ::  RC0PF0(:,:)                        !
  real(r8),allocatable ::  RC1PF0(:,:)                        !
  real(r8),allocatable ::  RC2PF0(:,:)                        !
  real(r8),allocatable ::  RM1PF0(:,:)                        !
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
  allocate(RQRAL(2,2,JV,JH));   RQRAL=0._r8
  allocate(RQRFE(2,2,JV,JH));   RQRFE=0._r8
  allocate(RQRHY(2,2,JV,JH));   RQRHY=0._r8
  allocate(RQRCA(2,2,JV,JH));   RQRCA=0._r8
  allocate(RQRMG(2,2,JV,JH));   RQRMG=0._r8
  allocate(RQRNA(2,2,JV,JH));   RQRNA=0._r8
  allocate(RQRKA(2,2,JV,JH));   RQRKA=0._r8
  allocate(RQROH(2,2,JV,JH));   RQROH=0._r8
  allocate(RQRSO(2,2,JV,JH));   RQRSO=0._r8
  allocate(RQRCL(2,2,JV,JH));   RQRCL=0._r8
  allocate(RQRC3(2,2,JV,JH));   RQRC3=0._r8
  allocate(RQRHC(2,2,JV,JH));   RQRHC=0._r8
  allocate(RQRAL1(2,2,JV,JH));  RQRAL1=0._r8
  allocate(RQRAL2(2,2,JV,JH));  RQRAL2=0._r8
  allocate(RQRAL3(2,2,JV,JH));  RQRAL3=0._r8
  allocate(RQRAL4(2,2,JV,JH));  RQRAL4=0._r8
  allocate(RQRALS(2,2,JV,JH));  RQRALS=0._r8
  allocate(RQRFE1(2,2,JV,JH));  RQRFE1=0._r8
  allocate(RQRFE2(2,2,JV,JH));  RQRFE2=0._r8
  allocate(RQRFE3(2,2,JV,JH));  RQRFE3=0._r8
  allocate(RQRFE4(2,2,JV,JH));  RQRFE4=0._r8
  allocate(RQRFES(2,2,JV,JH));  RQRFES=0._r8
  allocate(RQRCAO(2,2,JV,JH));  RQRCAO=0._r8
  allocate(RQRCAC(2,2,JV,JH));  RQRCAC=0._r8
  allocate(RQRCAH(2,2,JV,JH));  RQRCAH=0._r8
  allocate(RQRCAS(2,2,JV,JH));  RQRCAS=0._r8
  allocate(RQRMGO(2,2,JV,JH));  RQRMGO=0._r8
  allocate(RQRMGC(2,2,JV,JH));  RQRMGC=0._r8
  allocate(RQRMGH(2,2,JV,JH));  RQRMGH=0._r8
  allocate(RQRMGS(2,2,JV,JH));  RQRMGS=0._r8
  allocate(RQRNAC(2,2,JV,JH));  RQRNAC=0._r8
  allocate(RQRNAS(2,2,JV,JH));  RQRNAS=0._r8
  allocate(RQRKAS(2,2,JV,JH));  RQRKAS=0._r8
  allocate(RQRH0P(2,2,JV,JH));  RQRH0P=0._r8
  allocate(RQRH3P(2,2,JV,JH));  RQRH3P=0._r8
  allocate(RQRF1P(2,2,JV,JH));  RQRF1P=0._r8
  allocate(RQRF2P(2,2,JV,JH));  RQRF2P=0._r8
  allocate(RQRC0P(2,2,JV,JH));  RQRC0P=0._r8
  allocate(RQRC1P(2,2,JV,JH));  RQRC1P=0._r8
  allocate(RQRC2P(2,2,JV,JH));  RQRC2P=0._r8
  allocate(RQRM1P(2,2,JV,JH));  RQRM1P=0._r8
  allocate(trcsa_RQ(idsa_beg:idsab_end,2,JV,JH));     trcsa_RQ=0._r8
  allocate(RQRAL0(JY,JX));      RQRAL0=0._r8
  allocate(RQRFE0(JY,JX));      RQRFE0=0._r8
  allocate(RQRHY0(JY,JX));      RQRHY0=0._r8
  allocate(RQRCA0(JY,JX));      RQRCA0=0._r8
  allocate(RQRMG0(JY,JX));      RQRMG0=0._r8
  allocate(RQRNA0(JY,JX));      RQRNA0=0._r8
  allocate(RQRKA0(JY,JX));      RQRKA0=0._r8
  allocate(RQROH0(JY,JX));      RQROH0=0._r8
  allocate(RQRSO0(JY,JX));      RQRSO0=0._r8
  allocate(RQRCL0(JY,JX));      RQRCL0=0._r8
  allocate(RQRC30(JY,JX));      RQRC30=0._r8
  allocate(RQRHC0(JY,JX));      RQRHC0=0._r8
  allocate(RQRAL10(JY,JX));     RQRAL10=0._r8
  allocate(RQRAL20(JY,JX));     RQRAL20=0._r8
  allocate(RQRAL30(JY,JX));     RQRAL30=0._r8
  allocate(RQRAL40(JY,JX));     RQRAL40=0._r8
  allocate(RQRALS0(JY,JX));     RQRALS0=0._r8
  allocate(RQRFE10(JY,JX));     RQRFE10=0._r8
  allocate(RQRFE20(JY,JX));     RQRFE20=0._r8
  allocate(RQRFE30(JY,JX));     RQRFE30=0._r8
  allocate(RQRFE40(JY,JX));     RQRFE40=0._r8
  allocate(RQRFES0(JY,JX));     RQRFES0=0._r8
  allocate(RQRCAO0(JY,JX));     RQRCAO0=0._r8
  allocate(RQRCAC0(JY,JX));     RQRCAC0=0._r8
  allocate(RQRCAH0(JY,JX));     RQRCAH0=0._r8
  allocate(RQRCAS0(JY,JX));     RQRCAS0=0._r8
  allocate(RQRMGO0(JY,JX));     RQRMGO0=0._r8
  allocate(RQRMGC0(JY,JX));     RQRMGC0=0._r8
  allocate(RQRMGH0(JY,JX));     RQRMGH0=0._r8
  allocate(RQRMGS0(JY,JX));     RQRMGS0=0._r8
  allocate(RQRNAC0(JY,JX));     RQRNAC0=0._r8
  allocate(RQRNAS0(JY,JX));     RQRNAS0=0._r8
  allocate(RQRKAS0(JY,JX));     RQRKAS0=0._r8
  allocate(RQRH0P0(JY,JX));     RQRH0P0=0._r8
  allocate(RQRH3P0(JY,JX));     RQRH3P0=0._r8
  allocate(RQRF1P0(JY,JX));     RQRF1P0=0._r8
  allocate(RQRF2P0(JY,JX));     RQRF2P0=0._r8
  allocate(RQRC0P0(JY,JX));     RQRC0P0=0._r8
  allocate(RQRC1P0(JY,JX));     RQRC1P0=0._r8
  allocate(RQRC2P0(JY,JX));     RQRC2P0=0._r8
  allocate(RQRM1P0(JY,JX));     RQRM1P0=0._r8
  allocate(RALBLS(JS,JY,JX));   RALBLS=0._r8
  allocate(RFEBLS(JS,JY,JX));   RFEBLS=0._r8
  allocate(RHYBLS(JS,JY,JX));   RHYBLS=0._r8
  allocate(RCABLS(JS,JY,JX));   RCABLS=0._r8
  allocate(RMGBLS(JS,JY,JX));   RMGBLS=0._r8
  allocate(RNABLS(JS,JY,JX));   RNABLS=0._r8
  allocate(RKABLS(JS,JY,JX));   RKABLS=0._r8
  allocate(ROHBLS(JS,JY,JX));   ROHBLS=0._r8
  allocate(RSOBLS(JS,JY,JX));   RSOBLS=0._r8
  allocate(RCLBLS(JS,JY,JX));   RCLBLS=0._r8
  allocate(RC3BLS(JS,JY,JX));   RC3BLS=0._r8
  allocate(RHCBLS(JS,JY,JX));   RHCBLS=0._r8
  allocate(RAL1BS(JS,JY,JX));   RAL1BS=0._r8
  allocate(RAL2BS(JS,JY,JX));   RAL2BS=0._r8
  allocate(RAL3BS(JS,JY,JX));   RAL3BS=0._r8
  allocate(RAL4BS(JS,JY,JX));   RAL4BS=0._r8
  allocate(RALSBS(JS,JY,JX));   RALSBS=0._r8
  allocate(RFE1BS(JS,JY,JX));   RFE1BS=0._r8
  allocate(RFE2BS(JS,JY,JX));   RFE2BS=0._r8
  allocate(RFE3BS(JS,JY,JX));   RFE3BS=0._r8
  allocate(RFE4BS(JS,JY,JX));   RFE4BS=0._r8
  allocate(RFESBS(JS,JY,JX));   RFESBS=0._r8
  allocate(RCAOBS(JS,JY,JX));   RCAOBS=0._r8
  allocate(RCACBS(JS,JY,JX));   RCACBS=0._r8
  allocate(RCAHBS(JS,JY,JX));   RCAHBS=0._r8
  allocate(RCASBS(JS,JY,JX));   RCASBS=0._r8
  allocate(RMGOBS(JS,JY,JX));   RMGOBS=0._r8
  allocate(RMGCBS(JS,JY,JX));   RMGCBS=0._r8
  allocate(RMGHBS(JS,JY,JX));   RMGHBS=0._r8
  allocate(RMGSBS(JS,JY,JX));   RMGSBS=0._r8
  allocate(RNACBS(JS,JY,JX));   RNACBS=0._r8
  allocate(RNASBS(JS,JY,JX));   RNASBS=0._r8
  allocate(RKASBS(JS,JY,JX));   RKASBS=0._r8
  allocate(RH0PBS(JS,JY,JX));   RH0PBS=0._r8
  allocate(RH3PBS(JS,JY,JX));   RH3PBS=0._r8
  allocate(RF1PBS(JS,JY,JX));   RF1PBS=0._r8
  allocate(RF2PBS(JS,JY,JX));   RF2PBS=0._r8
  allocate(RC0PBS(JS,JY,JX));   RC0PBS=0._r8
  allocate(RC1PBS(JS,JY,JX));   RC1PBS=0._r8
  allocate(RC2PBS(JS,JY,JX));   RC2PBS=0._r8
  allocate(RM1PBS(JS,JY,JX));   RM1PBS=0._r8
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
  allocate(ZALH2(JZ,JY,JX));    ZALH2=0._r8
  allocate(ZFEH2(JZ,JY,JX));    ZFEH2=0._r8
  allocate(ZHYH2(JZ,JY,JX));    ZHYH2=0._r8
  allocate(ZCCH2(JZ,JY,JX));    ZCCH2=0._r8
  allocate(ZMAH2(JZ,JY,JX));    ZMAH2=0._r8
  allocate(ZNAH2(JZ,JY,JX));    ZNAH2=0._r8
  allocate(ZKAH2(JZ,JY,JX));    ZKAH2=0._r8
  allocate(ZOHH2(JZ,JY,JX));    ZOHH2=0._r8
  allocate(ZSO4H2(JZ,JY,JX));   ZSO4H2=0._r8
  allocate(ZCLH2(JZ,JY,JX));    ZCLH2=0._r8
  allocate(ZCO3H2(JZ,JY,JX));   ZCO3H2=0._r8
  allocate(ZHCOH2(JZ,JY,JX));   ZHCOH2=0._r8
  allocate(ZAL1H2(JZ,JY,JX));   ZAL1H2=0._r8
  allocate(ZAL2H2(JZ,JY,JX));   ZAL2H2=0._r8
  allocate(ZAL3H2(JZ,JY,JX));   ZAL3H2=0._r8
  allocate(ZAL4H2(JZ,JY,JX));   ZAL4H2=0._r8
  allocate(ZALSH2(JZ,JY,JX));   ZALSH2=0._r8
  allocate(ZFE1H2(JZ,JY,JX));   ZFE1H2=0._r8
  allocate(ZFE2H2(JZ,JY,JX));   ZFE2H2=0._r8
  allocate(ZFE3H2(JZ,JY,JX));   ZFE3H2=0._r8
  allocate(ZFE4H2(JZ,JY,JX));   ZFE4H2=0._r8
  allocate(ZFESH2(JZ,JY,JX));   ZFESH2=0._r8
  allocate(ZCAOH2(JZ,JY,JX));   ZCAOH2=0._r8
  allocate(ZCACH2(JZ,JY,JX));   ZCACH2=0._r8
  allocate(ZCAHH2(JZ,JY,JX));   ZCAHH2=0._r8
  allocate(ZCASH2(JZ,JY,JX));   ZCASH2=0._r8
  allocate(ZMGOH2(JZ,JY,JX));   ZMGOH2=0._r8
  allocate(ZMGCH2(JZ,JY,JX));   ZMGCH2=0._r8
  allocate(ZMGHH2(JZ,JY,JX));   ZMGHH2=0._r8
  allocate(ZMGSH2(JZ,JY,JX));   ZMGSH2=0._r8
  allocate(ZNACH2(JZ,JY,JX));   ZNACH2=0._r8
  allocate(ZNASH2(JZ,JY,JX));   ZNASH2=0._r8
  allocate(ZKASH2(JZ,JY,JX));   ZKASH2=0._r8
  allocate(H0P4H2(JZ,JY,JX));   H0P4H2=0._r8
  allocate(H3P4H2(JZ,JY,JX));   H3P4H2=0._r8
  allocate(ZF1PH2(JZ,JY,JX));   ZF1PH2=0._r8
  allocate(ZF2PH2(JZ,JY,JX));   ZF2PH2=0._r8
  allocate(ZC0PH2(JZ,JY,JX));   ZC0PH2=0._r8
  allocate(ZC1PH2(JZ,JY,JX));   ZC1PH2=0._r8
  allocate(ZC2PH2(JZ,JY,JX));   ZC2PH2=0._r8
  allocate(ZM1PH2(JZ,JY,JX));   ZM1PH2=0._r8
  allocate(H0PBH2(JZ,JY,JX));   H0PBH2=0._r8
  allocate(H3PBH2(JZ,JY,JX));   H3PBH2=0._r8
  allocate(ZF1BH2(JZ,JY,JX));   ZF1BH2=0._r8
  allocate(ZF2BH2(JZ,JY,JX));   ZF2BH2=0._r8
  allocate(ZC0BH2(JZ,JY,JX));   ZC0BH2=0._r8
  allocate(ZC1BH2(JZ,JY,JX));   ZC1BH2=0._r8
  allocate(ZC2BH2(JZ,JY,JX));   ZC2BH2=0._r8
  allocate(ZM1BH2(JZ,JY,JX));   ZM1BH2=0._r8
  allocate(RH0PXS(JZ,JY,JX));   RH0PXS=0._r8
  allocate(RH3PXS(JZ,JY,JX));   RH3PXS=0._r8
  allocate(RF1PXS(JZ,JY,JX));   RF1PXS=0._r8
  allocate(RF2PXS(JZ,JY,JX));   RF2PXS=0._r8
  allocate(RC0PXS(JZ,JY,JX));   RC0PXS=0._r8
  allocate(RC1PXS(JZ,JY,JX));   RC1PXS=0._r8
  allocate(RC2PXS(JZ,JY,JX));   RC2PXS=0._r8
  allocate(RM1PXS(JZ,JY,JX));   RM1PXS=0._r8
  allocate(RH0BXB(JZ,JY,JX));   RH0BXB=0._r8
  allocate(RH3BXB(JZ,JY,JX));   RH3BXB=0._r8
  allocate(RF1BXB(JZ,JY,JX));   RF1BXB=0._r8
  allocate(RF2BXB(JZ,JY,JX));   RF2BXB=0._r8
  allocate(RC0BXB(JZ,JY,JX));   RC0BXB=0._r8
  allocate(RC1BXB(JZ,JY,JX));   RC1BXB=0._r8
  allocate(RC2BXB(JZ,JY,JX));   RC2BXB=0._r8
  allocate(RM1BXB(JZ,JY,JX));   RM1BXB=0._r8
  allocate(RALFXS(JZ,JY,JX));   RALFXS=0._r8
  allocate(RFEFXS(JZ,JY,JX));   RFEFXS=0._r8
  allocate(RHYFXS(JZ,JY,JX));   RHYFXS=0._r8
  allocate(RCAFXS(JZ,JY,JX));   RCAFXS=0._r8
  allocate(RMGFXS(JZ,JY,JX));   RMGFXS=0._r8
  allocate(RNAFXS(JZ,JY,JX));   RNAFXS=0._r8
  allocate(RKAFXS(JZ,JY,JX));   RKAFXS=0._r8
  allocate(ROHFXS(JZ,JY,JX));   ROHFXS=0._r8
  allocate(RSOFXS(JZ,JY,JX));   RSOFXS=0._r8
  allocate(RCLFXS(JZ,JY,JX));   RCLFXS=0._r8
  allocate(RC3FXS(JZ,JY,JX));   RC3FXS=0._r8
  allocate(RHCFXS(JZ,JY,JX));   RHCFXS=0._r8
  allocate(RAL1XS(JZ,JY,JX));   RAL1XS=0._r8
  allocate(RAL2XS(JZ,JY,JX));   RAL2XS=0._r8
  allocate(RAL3XS(JZ,JY,JX));   RAL3XS=0._r8
  allocate(RAL4XS(JZ,JY,JX));   RAL4XS=0._r8
  allocate(RALSXS(JZ,JY,JX));   RALSXS=0._r8
  allocate(RFE1XS(JZ,JY,JX));   RFE1XS=0._r8
  allocate(RFE2XS(JZ,JY,JX));   RFE2XS=0._r8
  allocate(RFE3XS(JZ,JY,JX));   RFE3XS=0._r8
  allocate(RFE4XS(JZ,JY,JX));   RFE4XS=0._r8
  allocate(RFESXS(JZ,JY,JX));   RFESXS=0._r8
  allocate(RCACXS(JZ,JY,JX));   RCACXS=0._r8
  allocate(RCAOXS(JZ,JY,JX));   RCAOXS=0._r8
  allocate(RCAHXS(JZ,JY,JX));   RCAHXS=0._r8
  allocate(RCASXS(JZ,JY,JX));   RCASXS=0._r8
  allocate(RMGOXS(JZ,JY,JX));   RMGOXS=0._r8
  allocate(RMGCXS(JZ,JY,JX));   RMGCXS=0._r8
  allocate(RMGHXS(JZ,JY,JX));   RMGHXS=0._r8
  allocate(RMGSXS(JZ,JY,JX));   RMGSXS=0._r8
  allocate(RNACXS(JZ,JY,JX));   RNACXS=0._r8
  allocate(RNASXS(JZ,JY,JX));   RNASXS=0._r8
  allocate(RKASXS(JZ,JY,JX));   RKASXS=0._r8

  allocate(trcs_solsml2(idsa_beg:idsa_end,JS,JY,JX)); trcs_solsml2=0._r8

  allocate(RALFL0(JY,JX));      RALFL0=0._r8
  allocate(RFEFL0(JY,JX));      RFEFL0=0._r8
  allocate(RHYFL0(JY,JX));      RHYFL0=0._r8
  allocate(RCAFL0(JY,JX));      RCAFL0=0._r8
  allocate(RMGFL0(JY,JX));      RMGFL0=0._r8
  allocate(RNAFL0(JY,JX));      RNAFL0=0._r8
  allocate(RKAFL0(JY,JX));      RKAFL0=0._r8
  allocate(ROHFL0(JY,JX));      ROHFL0=0._r8
  allocate(RSOFL0(JY,JX));      RSOFL0=0._r8
  allocate(RCLFL0(JY,JX));      RCLFL0=0._r8
  allocate(RC3FL0(JY,JX));      RC3FL0=0._r8
  allocate(RHCFL0(JY,JX));      RHCFL0=0._r8
  allocate(RAL1F0(JY,JX));      RAL1F0=0._r8
  allocate(RAL2F0(JY,JX));      RAL2F0=0._r8
  allocate(RAL3F0(JY,JX));      RAL3F0=0._r8
  allocate(RAL4F0(JY,JX));      RAL4F0=0._r8
  allocate(RALSF0(JY,JX));      RALSF0=0._r8
  allocate(RFE1F0(JY,JX));      RFE1F0=0._r8
  allocate(RFE2F0(JY,JX));      RFE2F0=0._r8
  allocate(RFE3F0(JY,JX));      RFE3F0=0._r8
  allocate(RFE4F0(JY,JX));      RFE4F0=0._r8
  allocate(RFESF0(JY,JX));      RFESF0=0._r8
  allocate(RCAOF0(JY,JX));      RCAOF0=0._r8
  allocate(RCACF0(JY,JX));      RCACF0=0._r8
  allocate(RCAHF0(JY,JX));      RCAHF0=0._r8
  allocate(RCASF0(JY,JX));      RCASF0=0._r8
  allocate(RMGOF0(JY,JX));      RMGOF0=0._r8
  allocate(RMGCF0(JY,JX));      RMGCF0=0._r8
  allocate(RMGHF0(JY,JX));      RMGHF0=0._r8
  allocate(RMGSF0(JY,JX));      RMGSF0=0._r8
  allocate(RNACF0(JY,JX));      RNACF0=0._r8
  allocate(RNASF0(JY,JX));      RNASF0=0._r8
  allocate(RKASF0(JY,JX));      RKASF0=0._r8
  allocate(RH0PF0(JY,JX));      RH0PF0=0._r8
  allocate(RH3PF0(JY,JX));      RH3PF0=0._r8
  allocate(RF1PF0(JY,JX));      RF1PF0=0._r8
  allocate(RF2PF0(JY,JX));      RF2PF0=0._r8
  allocate(RC0PF0(JY,JX));      RC0PF0=0._r8
  allocate(RC1PF0(JY,JX));      RC1PF0=0._r8
  allocate(RC2PF0(JY,JX));      RC2PF0=0._r8
  allocate(RM1PF0(JY,JX));      RM1PF0=0._r8
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
  call destroy(RQRAL)
  call destroy(RQRFE)
  call destroy(RQRHY)
  call destroy(RQRCA)
  call destroy(RQRMG)
  call destroy(RQRNA)
  call destroy(RQRKA)
  call destroy(RQROH)
  call destroy(RQRSO)
  call destroy(RQRCL)
  call destroy(RQRC3)
  call destroy(RQRHC)
  call destroy(RQRAL1)
  call destroy(RQRAL2)
  call destroy(RQRAL3)
  call destroy(RQRAL4)
  call destroy(RQRALS)
  call destroy(RQRFE1)
  call destroy(RQRFE2)
  call destroy(RQRFE3)
  call destroy(RQRFE4)
  call destroy(RQRFES)
  call destroy(RQRCAO)
  call destroy(RQRCAC)
  call destroy(RQRCAH)
  call destroy(RQRCAS)
  call destroy(RQRMGO)
  call destroy(RQRMGC)
  call destroy(RQRMGH)
  call destroy(RQRMGS)
  call destroy(RQRNAC)
  call destroy(RQRNAS)
  call destroy(RQRKAS)
  call destroy(RQRH0P)
  call destroy(RQRH3P)
  call destroy(RQRF1P)
  call destroy(RQRF2P)
  call destroy(RQRC0P)
  call destroy(RQRC1P)
  call destroy(RQRC2P)
  call destroy(RQRM1P)
  call destroy(trcsa_RQ)
  call destroy(RQRAL0)
  call destroy(RQRFE0)
  call destroy(RQRHY0)
  call destroy(RQRCA0)
  call destroy(RQRMG0)
  call destroy(RQRNA0)
  call destroy(RQRKA0)
  call destroy(RQROH0)
  call destroy(RQRSO0)
  call destroy(RQRCL0)
  call destroy(RQRC30)
  call destroy(RQRHC0)
  call destroy(RQRAL10)
  call destroy(RQRAL20)
  call destroy(RQRAL30)
  call destroy(RQRAL40)
  call destroy(RQRALS0)
  call destroy(RQRFE10)
  call destroy(RQRFE20)
  call destroy(RQRFE30)
  call destroy(RQRFE40)
  call destroy(RQRFES0)
  call destroy(RQRCAO0)
  call destroy(RQRCAC0)
  call destroy(RQRCAH0)
  call destroy(RQRCAS0)
  call destroy(RQRMGO0)
  call destroy(RQRMGC0)
  call destroy(RQRMGH0)
  call destroy(RQRMGS0)
  call destroy(RQRNAC0)
  call destroy(RQRNAS0)
  call destroy(RQRKAS0)
  call destroy(RQRH0P0)
  call destroy(RQRH3P0)
  call destroy(RQRF1P0)
  call destroy(RQRF2P0)
  call destroy(RQRC0P0)
  call destroy(RQRC1P0)
  call destroy(RQRC2P0)
  call destroy(RQRM1P0)
  call destroy(RALBLS)
  call destroy(RFEBLS)
  call destroy(RHYBLS)
  call destroy(RCABLS)
  call destroy(RMGBLS)
  call destroy(RNABLS)
  call destroy(RKABLS)
  call destroy(ROHBLS)
  call destroy(RSOBLS)
  call destroy(RCLBLS)
  call destroy(RC3BLS)
  call destroy(RHCBLS)
  call destroy(RAL1BS)
  call destroy(RAL2BS)
  call destroy(RAL3BS)
  call destroy(RAL4BS)
  call destroy(RALSBS)
  call destroy(RFE1BS)
  call destroy(RFE2BS)
  call destroy(RFE3BS)
  call destroy(RFE4BS)
  call destroy(RFESBS)
  call destroy(RCAOBS)
  call destroy(RCACBS)
  call destroy(RCAHBS)
  call destroy(RCASBS)
  call destroy(RMGOBS)
  call destroy(RMGCBS)
  call destroy(RMGHBS)
  call destroy(RMGSBS)
  call destroy(RNACBS)
  call destroy(RNASBS)
  call destroy(RKASBS)
  call destroy(RH0PBS)
  call destroy(RH3PBS)
  call destroy(RF1PBS)
  call destroy(RF2PBS)
  call destroy(RC0PBS)
  call destroy(RC1PBS)
  call destroy(RC2PBS)
  call destroy(RM1PBS)

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
  call destroy(ZALH2)
  call destroy(ZFEH2)
  call destroy(ZHYH2)
  call destroy(ZCCH2)
  call destroy(ZMAH2)
  call destroy(ZNAH2)
  call destroy(ZKAH2)
  call destroy(ZOHH2)
  call destroy(ZSO4H2)
  call destroy(ZCLH2)
  call destroy(ZCO3H2)
  call destroy(ZHCOH2)
  call destroy(ZAL1H2)
  call destroy(ZAL2H2)
  call destroy(ZAL3H2)
  call destroy(ZAL4H2)
  call destroy(ZALSH2)
  call destroy(ZFE1H2)
  call destroy(ZFE2H2)
  call destroy(ZFE3H2)
  call destroy(ZFE4H2)
  call destroy(ZFESH2)
  call destroy(ZCAOH2)
  call destroy(ZCACH2)
  call destroy(ZCAHH2)
  call destroy(ZCASH2)
  call destroy(ZMGOH2)
  call destroy(ZMGCH2)
  call destroy(ZMGHH2)
  call destroy(ZMGSH2)
  call destroy(ZNACH2)
  call destroy(ZNASH2)
  call destroy(ZKASH2)
  call destroy(H0P4H2)
  call destroy(H3P4H2)
  call destroy(ZF1PH2)
  call destroy(ZF2PH2)
  call destroy(ZC0PH2)
  call destroy(ZC1PH2)
  call destroy(ZC2PH2)
  call destroy(ZM1PH2)
  call destroy(H0PBH2)
  call destroy(H3PBH2)
  call destroy(ZF1BH2)
  call destroy(ZF2BH2)
  call destroy(ZC0BH2)
  call destroy(ZC1BH2)
  call destroy(ZC2BH2)
  call destroy(ZM1BH2)
  call destroy(RH0PXS)
  call destroy(RH3PXS)
  call destroy(RF1PXS)
  call destroy(RF2PXS)
  call destroy(RC0PXS)
  call destroy(RC1PXS)
  call destroy(RC2PXS)
  call destroy(RM1PXS)
  call destroy(RH0BXB)
  call destroy(RH3BXB)
  call destroy(RF1BXB)
  call destroy(RF2BXB)
  call destroy(RC0BXB)
  call destroy(RC1BXB)
  call destroy(RC2BXB)
  call destroy(RM1BXB)
  call destroy(RALFXS)
  call destroy(RFEFXS)
  call destroy(RHYFXS)
  call destroy(RCAFXS)
  call destroy(RMGFXS)
  call destroy(RNAFXS)
  call destroy(RKAFXS)
  call destroy(ROHFXS)
  call destroy(RSOFXS)
  call destroy(RCLFXS)
  call destroy(RC3FXS)
  call destroy(RHCFXS)
  call destroy(RAL1XS)
  call destroy(RAL2XS)
  call destroy(RAL3XS)
  call destroy(RAL4XS)
  call destroy(RALSXS)
  call destroy(RFE1XS)
  call destroy(RFE2XS)
  call destroy(RFE3XS)
  call destroy(RFE4XS)
  call destroy(RFESXS)
  call destroy(RCACXS)
  call destroy(RCAOXS)
  call destroy(RCAHXS)
  call destroy(RCASXS)
  call destroy(RMGOXS)
  call destroy(RMGCXS)
  call destroy(RMGHXS)
  call destroy(RMGSXS)
  call destroy(RNACXS)
  call destroy(RNASXS)
  call destroy(RKASXS)

  call destroy(RALFL0)
  call destroy(RFEFL0)
  call destroy(RHYFL0)
  call destroy(RCAFL0)
  call destroy(RMGFL0)
  call destroy(RNAFL0)
  call destroy(RKAFL0)
  call destroy(ROHFL0)
  call destroy(RSOFL0)
  call destroy(RCLFL0)
  call destroy(RC3FL0)
  call destroy(RHCFL0)
  call destroy(RAL1F0)
  call destroy(RAL2F0)
  call destroy(RAL3F0)
  call destroy(RAL4F0)
  call destroy(RALSF0)
  call destroy(RFE1F0)
  call destroy(RFE2F0)
  call destroy(RFE3F0)
  call destroy(RFE4F0)
  call destroy(RFESF0)
  call destroy(RCAOF0)
  call destroy(RCACF0)
  call destroy(RCAHF0)
  call destroy(RCASF0)
  call destroy(RMGOF0)
  call destroy(RMGCF0)
  call destroy(RMGHF0)
  call destroy(RMGSF0)
  call destroy(RNACF0)
  call destroy(RNASF0)
  call destroy(RKASF0)
  call destroy(RH0PF0)
  call destroy(RH3PF0)
  call destroy(RF1PF0)
  call destroy(RF2PF0)
  call destroy(RC0PF0)
  call destroy(RC1PF0)
  call destroy(RC2PF0)
  call destroy(RM1PF0)
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
  end subroutine DestructTrnsfrsData

end module TrnsfrsDataMod
