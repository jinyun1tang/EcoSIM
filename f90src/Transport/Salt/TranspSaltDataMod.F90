module TranspSaltDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use SoilPropertyDataType
  use TracerIDMod
  use AqueChemDatatype
implicit none
  public
  CHARACTER(LEN=*),private,PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), PARAMETER :: VFLWX=0.5_r8
  real(r8),allocatable ::  AquaIonDifusivty2_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcSalt_RQR(:,:,:,:,:)                     !
  real(r8),allocatable ::  trcSalt_RQ(:,:,:,:)                       !
  real(r8),allocatable ::  trcSalt_RQR0(:,:,:)                        !

  real(r8),allocatable ::  trcSalt_solml2(:,:,:,:)              !
  real(r8),allocatable ::  ZFE2(:,:,:)                        !
  real(r8),allocatable ::  ZHCO32(:,:,:)                      !
  real(r8),allocatable ::  trcSalt_solml2R(:,:,:,:)                       !
  real(r8),allocatable ::  trcSaltAdv2SowLay(:,:,:,:)                      !

  real(r8),allocatable ::  trcSalt_TBLS(:,:,:,:)

  real(r8),allocatable ::  POSGL2(:,:,:)                      !
  real(r8),allocatable ::  trcSalt3DFlo2CellM(:,:,:,:,:)                    !
  real(r8),allocatable ::  trcSalt_RFHS(:,:,:,:,:)                    !
  real(r8),allocatable ::  trcSalt_soHml2(:,:,:,:)                       !
  real(r8),allocatable ::  trcSalt_RFXS(:,:,:,:)                      !
  real(r8), allocatable ::  trc_Saltml2_snvr(:,:,:,:)               ! snowpack salt dissolved tracers
  real(r8),allocatable ::  trcSalt_RFL0(:,:,:)                        !
  real(r8),allocatable ::  trcSalt_RFL1(:,:,:)                        !
!----------------------------------------------------------------------

contains
  subroutine InitTranspSaltData

  implicit none
  allocate(AquaIonDifusivty2_vr(idsalt_beg:idsalt_mend,JZ,JY,JX));   AquaIonDifusivty2_vr=0._r8

  allocate(trcSalt_RQR(idsalt_beg:idsalt_end,2,2,JV,JH));   trcSalt_RQR=0._r8
  allocate(trcSalt_RQ(idsalt_beg:idsaltb_end,2,JV,JH));     trcSalt_RQ=0._r8
  allocate(trcSalt_RQR0(idsalt_beg:idsalt_end,JY,JX));     trcSalt_RQR0=0._r8
  allocate(trcSaltAdv2SowLay(idsalt_beg:idsalt_end,JS,JY,JX));   trcSaltAdv2SowLay=0._r8
  allocate(trcSalt_TBLS(idsalt_beg:idsalt_end,JS,JY,JX)); trcSalt_TBLS=0._r8

  allocate(POSGL2(JZ,JY,JX));   POSGL2=0._r8
  allocate(trcSalt3DFlo2CellM(idsalt_beg:idsaltb_end,3,0:JD,JV,JH));trcSalt3DFlo2CellM=0._r8
  allocate(trcSalt_RFHS(idsalt_beg:idsaltb_end,3,JD,JV,JH)); trcSalt_RFHS=0._r8
  allocate(trcSalt_soHml2(idsalt_beg:idsaltb_end,JZ,JY,JX));    trcSalt_soHml2=0._r8
  allocate(trcSalt_RFXS(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_RFXS=0._r8

  allocate(trc_Saltml2_snvr(idsalt_beg:idsalt_end,JS,JY,JX)); trc_Saltml2_snvr=0._r8

  allocate(trcSalt_RFL0(idsalt_beg:idsalt_end,JY,JX));      trcSalt_RFL0=0._r8
  allocate(trcSalt_RFL1(idsalt_beg:idsaltb_end,JY,JX));     trcSalt_RFL1=0._r8
  end subroutine InitTranspSaltData

!----------------------------------------------------------------------
  subroutine DestructTranspSaltData
  use abortutils, only : destroy
  implicit none

  call destroy(trcSalt3DFlo2CellM)
  call destroy(trcSalt_solml2)
  call destroy(trcSalt_solml2R)

  call destroy(trcSalt_RQ)

  call destroy(POSGL2)


  call destroy(trcSalt_RQR)
  call destroy(trcSalt_RFXS)
  end subroutine DestructTranspSaltData

end module TranspSaltDataMod
