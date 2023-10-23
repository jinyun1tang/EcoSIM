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
  real(r8), allocatable ::  trcSalt_sosml2(:,:,:,:)               ! snowpack salt dissolved tracers
  real(r8),allocatable ::  trcSalt_RFL0(:,:,:)                        !
  real(r8),allocatable ::  trcSalt_RFL1(:,:,:)                        !
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

  allocate(trcSalt_sosml2(idsalt_beg:idsalt_end,JS,JY,JX)); trcSalt_sosml2=0._r8

  allocate(trcSalt_RFL0(idsalt_beg:idsalt_end,JY,JX));      trcSalt_RFL0=0._r8
  allocate(trcSalt_RFL1(idsalt_beg:idsaltb_end,JY,JX));     trcSalt_RFL1=0._r8
  end subroutine InitTrnsfrsData

!----------------------------------------------------------------------
  subroutine DestructTrnsfrsData
  use abortutils, only : destroy
  implicit none

  call destroy(trcSalt3DFlo2CellM)
  call destroy(trcSalt_solml2)
  call destroy(trcSalt_solml2R)

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
  call destroy(trcSalt_RQ)

  call destroy(POSGL2)


  call destroy(trcSalt_RQR)
  call destroy(trcSalt_RFXS)
  end subroutine DestructTrnsfrsData

end module TranspSaltDataMod
