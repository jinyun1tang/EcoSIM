module SurfPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  save
  character(len=*), private, parameter :: mod_filename=__FILE__

  real(r8),allocatable ::  XVGeomLayer(:,:)                         !
  real(r8),allocatable ::  XVLWatMicP(:,:)                         !
  real(r8),allocatable ::  XVLiceMicP(:,:)                         !
  real(r8),allocatable ::  THRMR(:,:)                         !
  real(r8),allocatable ::  VLPoreLitR(:,:)                        !
  real(r8),allocatable ::  LWRad2LitR(:,:)                         !
  real(r8),allocatable ::  THRMS(:,:)                         !
  real(r8),allocatable ::  RADXR(:,:)                         !
  real(r8),allocatable ::  LWRad2Grnd(:,:)                         !
  real(r8),allocatable ::  RADXG(:,:)                         !
  real(r8),allocatable ::  RAGR(:,:)                          !
  real(r8),allocatable ::  PAREG(:,:)                         !
  real(r8),allocatable ::  PARER(:,:)                         !
  real(r8),allocatable ::  RARG(:,:)                          !
  real(r8),allocatable ::  PARSR(:,:)                         !
  real(r8),allocatable ::  RAR(:,:)                           !
  real(r8),allocatable ::  FLWRL(:,:)                         !
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

  public :: InitSurfPhysData,DestructSurfPhysData

  contains
!------------------------------------------------------------------------------------------  

  subroutine InitSurfPhysData  
  implicit none

  allocate(XVGeomLayer(JY,JX));       XVGeomLayer=0._r8
  allocate(XVLWatMicP(JY,JX));       XVLWatMicP=0._r8
  allocate(XVLiceMicP(JY,JX));       XVLiceMicP=0._r8
  allocate(THRMR(JY,JX));       THRMR=0._r8
  allocate(VLPoreLitR(JY,JX));      VLPoreLitR=0._r8
  allocate(LWRad2LitR(JY,JX));       LWRad2LitR=0._r8
  allocate(THRMS(JY,JX));       THRMS=0._r8
  allocate(RADXR(JY,JX));       RADXR=0._r8    
  allocate(LWRad2Grnd(JY,JX));       LWRad2Grnd=0._r8
  allocate(RADXG(JY,JX));       RADXG=0._r8
  allocate(RAGR(JY,JX));        RAGR=0._r8
  allocate(PAREG(JY,JX));       PAREG=0._r8
  allocate(PARER(JY,JX));       PARER=0._r8
  allocate(RARG(JY,JX));        RARG=0._r8
  allocate(PARSR(JY,JX));       PARSR=0._r8
  allocate(RAR(JY,JX));         RAR=0._r8  
  allocate(FLWRL(JY,JX));       FLWRL=0._r8
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

  call destroy(XVGeomLayer)
  call destroy(XVLWatMicP)
  call destroy(XVLiceMicP)
  call destroy(THRMR)
  call destroy(VLPoreLitR)
  call destroy(LWRad2LitR)
  call destroy(THRMS)
  call destroy(RADXR)  
  call destroy(LWRad2Grnd)
  call destroy(RADXG)
  call destroy(RAGR)
  call destroy(PAREG)
  call destroy(PARER)
  call destroy(RARG)
  call destroy(PARSR)
  call destroy(RAR)  
  call destroy(FLWRL)
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
