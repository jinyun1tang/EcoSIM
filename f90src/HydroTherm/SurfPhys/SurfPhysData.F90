module SurfPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  save
  character(len=*), private, parameter :: mod_filename=&
  __FILE__

  real(r8),allocatable ::  XVLMobileWaterLitR(:,:)                         !
  real(r8),allocatable ::  XVLMobileWatMicP(:,:)                         !
  real(r8),allocatable ::  XVLiceMicP_col(:,:)                         !
  real(r8),allocatable ::  LWEmscefLitR(:,:)                         !
  real(r8),allocatable ::  VLPoreLitR(:,:)                        !
  real(r8),allocatable ::  LWRad2LitR(:,:)                         !
  real(r8),allocatable ::  LWEmscefSoil(:,:)                         !
  real(r8),allocatable ::  RadSWonLitR(:,:)                         !
  real(r8),allocatable ::  LWRad2Grnd(:,:)                         !
  real(r8),allocatable ::  RadSWonSoi(:,:)                         !
  real(r8),allocatable ::  RAGR(:,:)                          !
  real(r8),allocatable ::  PAREG(:,:)                         !
  real(r8),allocatable ::  PARER(:,:)                         !
  real(r8),allocatable ::  RARG(:,:)                          !
  real(r8),allocatable ::  PARSR(:,:)                         !
  real(r8),allocatable ::  VapDiffusResistanceLitR(:,:)                           !
  real(r8),allocatable ::  WatFLow2LitR_col(:,:)                         !
  real(r8),allocatable ::  HeatFLoByWat2LitRi(:,:)                        !  
  real(r8),allocatable ::  PARSG(:,:)                         !
  real(r8),allocatable ::  VapXAir2LitR(:,:)                  !water vapor flux from canopy air to litr
!  real(r8),allocatable ::  VapXAir2TopLay(:,:)                !water vapor flux from canopy air to top layer of soi/lake
  real(r8),allocatable ::  LitrIceFlxThaw(:,:)                !Water flux from ice thaw in litter 
  real(r8),allocatable ::  LitrIceHeatFlxFrez(:,:)            !Heat associated with ice freeze in litter (>0 freeze) 
  real(r8),allocatable ::  RAGS(:,:)                          !    
  real(r8),allocatable ::  CVRDW(:,:)                         !
  real(r8),allocatable ::  Prec2SoiMacP1(:,:)                          !
  real(r8),allocatable ::  PRECM(:,:)                         !
  real(r8),allocatable ::  Prec2SoiMicP1(:,:)                          !
  real(r8),allocatable ::  PrecHeat2SoiMicP1(:,:)                        !
  real(r8),allocatable ::  PrecHeat2LitR1(:,:)                        !
  real(r8),allocatable ::  Prec2LitR1(:,:)                          !
  real(r8),allocatable ::  BAREW(:,:)                         !
  real(r8),allocatable ::  HCNDR(:,:)                         !

  real(r8),allocatable :: watflw(:,:)
  real(r8),allocatable :: waticefl(:,:)

  public :: InitSurfPhysData,DestructSurfPhysData

  contains
!------------------------------------------------------------------------------------------  

  subroutine InitSurfPhysData  
  implicit none

  allocate(watflw(JY, JX))
  allocate(waticefl(JY,JX))

  allocate(XVLMobileWaterLitR(JY,JX));       XVLMobileWaterLitR=0._r8
  allocate(XVLMobileWatMicP(JY,JX));       XVLMobileWatMicP=0._r8
  allocate(XVLiceMicP_col(JY,JX));       XVLiceMicP_col=0._r8
  allocate(LWEmscefLitR(JY,JX));       LWEmscefLitR=0._r8
  allocate(VLPoreLitR(JY,JX));      VLPoreLitR=0._r8
  allocate(LWRad2LitR(JY,JX));       LWRad2LitR=0._r8
  allocate(LWEmscefSoil(JY,JX));       LWEmscefSoil=0._r8
  allocate(RadSWonLitR(JY,JX));       RadSWonLitR=0._r8    
  allocate(LWRad2Grnd(JY,JX));       LWRad2Grnd=0._r8
  allocate(RadSWonSoi(JY,JX));       RadSWonSoi=0._r8
  allocate(RAGR(JY,JX));        RAGR=0._r8
  allocate(PAREG(JY,JX));       PAREG=0._r8
  allocate(PARER(JY,JX));       PARER=0._r8
  allocate(RARG(JY,JX));        RARG=0._r8
  allocate(PARSR(JY,JX));       PARSR=0._r8
  allocate(VapDiffusResistanceLitR(JY,JX));         VapDiffusResistanceLitR=0._r8  
  allocate(WatFLow2LitR_col(JY,JX));       WatFLow2LitR_col=0._r8
  allocate(HeatFLoByWat2LitRi(JY,JX));      HeatFLoByWat2LitRi=0._r8  
  allocate(PARSG(JY,JX));       PARSG=0._r8  
  allocate(VapXAir2LitR(JY,JX));       VapXAir2LitR=0._r8  
!  allocate(VapXAir2TopLay(JY,JX));       VapXAir2TopLay=0._r8  
  allocate(LitrIceFlxThaw(JY,JX));       LitrIceFlxThaw=0._r8  
  allocate(RAGS(JY,JX));        RAGS=0._r8  
  allocate(LitrIceHeatFlxFrez(JY,JX));       LitrIceHeatFlxFrez=0._r8  
  allocate(CVRDW(JY,JX));       CVRDW=0._r8
  allocate(Prec2SoiMacP1(JY,JX));        Prec2SoiMacP1=0._r8
  allocate(PRECM(JY,JX));       PRECM=0._r8
  allocate(Prec2SoiMicP1(JY,JX));        Prec2SoiMicP1=0._r8
  allocate(PrecHeat2SoiMicP1(JY,JX));      PrecHeat2SoiMicP1=0._r8
  allocate(PrecHeat2LitR1(JY,JX));      PrecHeat2LitR1=0._r8  
  allocate(Prec2LitR1(JY,JX));        Prec2LitR1=0._r8  
  allocate(BAREW(JY,JX));       BAREW=0._r8  
  allocate(HCNDR(JY,JX));       HCNDR=0._r8  
  end subroutine InitSurfPhysData  
!------------------------------------------------------------------------------------------  

  subroutine DestructSurfPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(XVLMobileWaterLitR)
  call destroy(XVLMobileWatMicP)
  call destroy(XVLiceMicP_col)
  call destroy(LWEmscefLitR)
  call destroy(VLPoreLitR)
  call destroy(LWRad2LitR)
  call destroy(LWEmscefSoil)
  call destroy(RadSWonLitR)  
  call destroy(LWRad2Grnd)
  call destroy(RadSWonSoi)
  call destroy(RAGR)
  call destroy(PAREG)
  call destroy(PARER)
  call destroy(RARG)
  call destroy(PARSR)
  call destroy(VapDiffusResistanceLitR)  
  call destroy(WatFLow2LitR_col)
  call destroy(HeatFLoByWat2LitRi)  
  call destroy(PARSG)  
  call destroy(VapXAir2LitR)
!  call destroy(VapXAir2TopLay)
  call destroy(LitrIceFlxThaw)  
  call destroy(RAGS) 
  call destroy(LitrIceHeatFlxFrez)
  call destroy(CVRDW)
  call destroy(Prec2SoiMacP1)
  call destroy(PRECM)
  call destroy(Prec2SoiMicP1)
  call destroy(PrecHeat2SoiMicP1)
  call destroy(PrecHeat2LitR1)
  call destroy(Prec2LitR1)
  call destroy(BAREW)
  call destroy(HCNDR)  
  end subroutine DestructSurfPhysData

end module SurfPhysData
