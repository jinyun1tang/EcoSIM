module TrnsfrsDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use SoilPropertyDataType
  use TracerIDMod
  use AqueChemDatatype
implicit none
  public
  CHARACTER(LEN=*),private,PARAMETER :: MOD_FILENAME=__FILE__

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
  real(r8),allocatable ::  trcsa_RFHS(:,:,:,:,:)                    !
  real(r8),allocatable ::  trcsa_soHml2(:,:,:,:)                       !
  real(r8),allocatable ::  trcsa_RFXS(:,:,:,:)                      !
  real(r8), allocatable ::  trcsa_sosml2(:,:,:,:)               ! snowpack salt dissolved tracers
  real(r8),allocatable ::  trcsa_RFL0(:,:,:)                        !
  real(r8),allocatable ::  trcsa_RFL1(:,:,:)                        !
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
  allocate(trcsa_RFHS(idsa_beg:idsab_end,3,JD,JV,JH)); trcsa_RFHS=0._r8
  allocate(trcsa_soHml2(idsa_beg:idsab_end,JZ,JY,JX));    trcsa_soHml2=0._r8
  allocate(trcsa_RFXS(idsa_beg:idsab_end,JZ,JY,JX));   trcsa_RFXS=0._r8

  allocate(trcsa_sosml2(idsa_beg:idsa_end,JS,JY,JX)); trcsa_sosml2=0._r8

  allocate(trcsa_RFL0(idsa_beg:idsa_end,JY,JX));      trcsa_RFL0=0._r8
  allocate(trcsa_RFL1(idsa_beg:idsab_end,JY,JX));     trcsa_RFL1=0._r8
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


  call destroy(trcsa_RQR)
  call destroy(trcsa_RFXS)
  end subroutine DestructTrnsfrsData

end module TrnsfrsDataMod
