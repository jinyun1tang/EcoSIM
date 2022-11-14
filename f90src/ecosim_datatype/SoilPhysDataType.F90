module SoilPhysDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  SLOPE(:,:,:)                       !slope	in three directions [o]
  real(r8),target,allocatable ::  FC(:,:,:)                          !water contents at field capacity
  real(r8),target,allocatable ::  WP(:,:,:)                          !water contents at wilting point
  real(r8),target,allocatable ::  SCNV(:,:,:)                        !soil vertical saturated hydraulic conductivity [mm h-1]
  real(r8),target,allocatable ::  SCNH(:,:,:)                        !soil horizontal saturated hydraulic conductivity, [mm h-1]
  real(r8),target,allocatable ::  PSIFC(:,:)                         !water potentials at field capacity, [MPa]
  real(r8),target,allocatable ::  PSIWP(:,:)                         !water potentials at wilting point [MPa]
  real(r8),target,allocatable ::  THW(:,:,:)                         !initial soil water content
  real(r8),target,allocatable ::  THI(:,:,:)                         !initial ice content
  REAL(R8),target,allocatable ::  ALBX(:,:)                          !Surface albedo
  real(r8),target,allocatable ::  PSL(:,:,:)                         !log soil porosity	-
  real(r8),target,allocatable ::  FCL(:,:,:)                         !log water content at field capacity
  real(r8),target,allocatable ::  WPL(:,:,:)                         !log water content at wilting point
  real(r8),target,allocatable ::  PSD(:,:,:)                         !log soil porosity - log water content at field capacity
  real(r8),target,allocatable ::  FCD(:,:,:)                         !log water content at field capacity
  real(r8),target,allocatable ::  SRP(:,:,:)                         !shape parameter for water desorption
  real(r8),target,allocatable ::  FSLOPE(:,:,:)                      !fraction of slope in 1 and 2
  REAL(R8),target,allocatable ::  VOLAI(:,:,:)                       !initial total soil micropore porosity	m3 d-2
  REAL(R8),target,allocatable ::  PSIMS(:,:)                         !log water potential at saturation	MPa
  REAL(R8),target,allocatable ::  PSIMX(:,:)                         !log water potential at field capacity	-
  REAL(R8),target,allocatable ::  PSIMN(:,:)                         !log water potential at wilting point
  REAL(R8),target,allocatable ::  PSISD(:,:)                         !log water potential at field capacity 	-
  REAL(R8),target,allocatable ::  PSIMD(:,:)                         !log water potential at saturation - log water potential at field capacity
  real(r8),target,allocatable ::  VHCM(:,:,:)                        !soil solid heat capacity [MPa m-3 K-1]
  real(r8),target,allocatable ::  DPTHA(:,:)                         !active layer depth, [n]
!----------------------------------------------------------------------

contains
  subroutine InitSoilPhysData

  implicit none
  allocate(SLOPE(0:3,JY,JX));   SLOPE=0._r8
  allocate(FC(0:JZ,JY,JX));     FC=0._r8
  allocate(WP(0:JZ,JY,JX));     WP=0._r8
  allocate(SCNV(0:JZ,JY,JX));   SCNV=0._r8
  allocate(SCNH(JZ,JY,JX));     SCNH=0._r8
  allocate(PSIFC(JY,JX));       PSIFC=0._r8
  allocate(PSIWP(JY,JX));       PSIWP=0._r8
  allocate(THW(JZ,JY,JX));      THW=0._r8
  allocate(THI(JZ,JY,JX));      THI=0._r8
  allocate(ALBX(JY,JX));        ALBX=0._r8
  allocate(PSL(0:JZ,JY,JX));    PSL=0._r8
  allocate(FCL(0:JZ,JY,JX));    FCL=0._r8
  allocate(WPL(0:JZ,JY,JX));    WPL=0._r8
  allocate(PSD(0:JZ,JY,JX));    PSD=0._r8
  allocate(FCD(0:JZ,JY,JX));    FCD=0._r8
  allocate(SRP(0:JZ,JY,JX));    SRP=0._r8
  allocate(FSLOPE(2,JY,JX));    FSLOPE=0._r8
  allocate(VOLAI(0:JZ,JY,JX));  VOLAI=0._r8
  allocate(PSIMS(JY,JX));       PSIMS=0._r8
  allocate(PSIMX(JY,JX));       PSIMX=0._r8
  allocate(PSIMN(JY,JX));       PSIMN=0._r8
  allocate(PSISD(JY,JX));       PSISD=0._r8
  allocate(PSIMD(JY,JX));       PSIMD=0._r8
  allocate(VHCM(0:JZ,JY,JX));   VHCM=0._r8
  allocate(DPTHA(JY,JX));       DPTHA=0._r8
  end subroutine InitSoilPhysData

!----------------------------------------------------------------------
  subroutine DestructSoilPhysData
  use abortutils, only : destroy
  implicit none
  call destroy(SLOPE)
  call destroy(FC)
  call destroy(WP)
  call destroy(SCNV)
  call destroy(SCNH)
  call destroy(PSIFC)
  call destroy(PSIWP)
  call destroy(THW)
  call destroy(THI)
  call destroy(ALBX)
  call destroy(PSL)
  call destroy(FCL)
  call destroy(WPL)
  call destroy(PSD)
  call destroy(FCD)
  call destroy(SRP)
  call destroy(FSLOPE)
  call destroy(VOLAI)
  call destroy(PSIMS)
  call destroy(PSIMX)
  call destroy(PSIMN)
  call destroy(PSISD)
  call destroy(PSIMD)
  call destroy(VHCM)
  call destroy(DPTHA)
  end subroutine DestructSoilPhysData

end module SoilPhysDataType
