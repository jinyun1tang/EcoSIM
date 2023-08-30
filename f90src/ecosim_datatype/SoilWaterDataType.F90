module SoilWaterDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  THETP(:,:,:)                      !air concentration [m3 m-3]
  real(r8),target,allocatable ::  VsoiP(:,:,:)                       !soil air content [m3 d-2]
  real(r8),target,allocatable ::  THETW(:,:,:)                      !volumetric water content [m3 m-3]
  real(r8),target,allocatable ::  THETI(:,:,:)                      !volumetric ice content [m3 m-3]
  real(r8),target,allocatable ::  THETWZ(:,:,:)                     !volumetric moblize water [m3 m-3]
  real(r8),target,allocatable ::  THETIZ(:,:,:)                     !volumetric mobile ice [m3 m-3]
  real(r8),target,allocatable ::  VWatMicP(:,:,:)                       !soil micropore water content [m3 d-2]
  real(r8),target,allocatable ::  ViceMicP(:,:,:)                       !soil micropore ice content   [m3 d-2]
  real(r8),target,allocatable ::  VWatMacP(:,:,:)                      !soil macropore water content [m3 d-2]
  real(r8),target,allocatable ::  PSISM(:,:,:)                      !soil micropore matric water potential [MPa]
  real(r8),target,allocatable ::  TotalSoilH2OPSIMPa(:,:,:)                      !soil micropore total water potential [MPa]
  real(r8),target,allocatable ::  VWatMicPX(:,:,:)                      !soil micropore water content before wetting front [m3 d-2]
  real(r8),target,allocatable ::  FINH(:,:,:)                       !soil macropore - micropore water transfer [m3 d-2 h-1]
  real(r8),target,allocatable ::  ViceMacP(:,:,:)                      !soil macropore ice content [m3 d-2]
  real(r8),target,allocatable ::  VWatMicPM(:,:,:,:)                    !soil micropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VWatMacPM(:,:,:,:)                   !soil macropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VsoiPM(:,:,:,:)                    !soil air content, [m3 d-2]
  real(r8),target,allocatable ::  FILM(:,:,:,:)                     !soil water film thickness , [m]
  real(r8),target,allocatable ::  DTBLG(:,:)                        !slope of water table relative to surface slope, [-]
  real(r8),target,allocatable ::  DTBLDI(:,:)                       !depth of artificial water table
  real(r8),target,allocatable ::  DTBLY(:,:)                        !artificial water table depth, [m]
  real(r8),target,allocatable ::  DTBLD(:,:)                        !depth of artificial water table adjusted for elevation
  real(r8),target,allocatable ::  DPTHT(:,:)                        !internal water table depth, [m]
  real(r8),target,allocatable ::  DTBLZ(:,:)                        !external water table depth, [m]
  real(r8),target,allocatable ::  DTBLX(:,:)                        !external water table depth, [m]
  real(r8),target,allocatable ::  DTBLI(:,:)                        !external water table depth, [m]
  real(r8),target,allocatable ::  ENGYPM(:,:,:)                     !total energy impact for erosion
  real(r8),target,allocatable ::  XVOLTM(:,:,:)                     !excess water+ice
  real(r8),target,allocatable ::  XVWatMicPM(:,:,:)                     !excess water
  real(r8),target,allocatable ::  XViceMicPM(:,:,:)                     !excess ice
  real(r8),target,allocatable ::  HCND(:,:,:,:,:)                   !saturated hydraulic conductivity
  real(r8),target,allocatable ::  CNDH(:,:,:)                       !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),target,allocatable ::  CNDU(:,:,:)                       !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8),target,allocatable ::  QRM(:,:,:)                        !runoff water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  QRV(:,:,:)                        !runoff velocity, [m t-1]
   integer,target,allocatable ::  IFLBM(:,:,:,:,:)                   !flag for directional surface runoff
   integer,target,allocatable ::  IRCHG(:,:,:,:)                     !enables or disables boundary water flux depending on aspect, [-]
   integer,target,allocatable ::  IFLBH(:,:,:,:)                     !flag for directional runoff, related to IFLBM
  real(r8),target,allocatable ::  RCHGNU(:,:)                       !northern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHGEU(:,:)                       !eastern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHGSU(:,:)                       !southern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHGWU(:,:)                       !western subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHGNT(:,:)                       !northern subsurface boundary water flux rate constant, [h-1]
  real(r8),target,allocatable ::  RCHGET(:,:)                       !eastern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RCHGST(:,:)                       !southern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RCHGWT(:,:)                       !western subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RCHQN(:,:)                        !northern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHQE(:,:)                        !eastern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHQS(:,:)                        !southern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHQW(:,:)                        !western surface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHGD(:,:)                        !lower subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  FLWM(:,:,:,:,:)                   !micropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  WaterFlowMacPi(:,:,:,:,:)                  !macropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  FLPM(:,:,:,:)                     !soil air flux, [g d-2 t-1]
  real(r8),target,allocatable ::  FINHM(:,:,:,:)                    !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8),target,allocatable ::  FLQSM(:,:,:)                      !meltwater flux into soil micropores
  real(r8),target,allocatable ::  FLQHM(:,:,:)                      !meltwater flux into soil macropores
  real(r8),target,allocatable ::  THETPM(:,:,:,:)                   !soil air-filled porosity, [m3 m-3]
  real(r8),target,allocatable ::  TORT(:,:,:,:)                     !soil tortuosity, []
  real(r8),target,allocatable ::  TORTH(:,:,:,:)                    !macropore tortuosity, []
  real(r8),target,allocatable ::  DFGS(:,:,:,:)                     !coefficient for dissolution - volatilization, []
  real(r8),target,allocatable ::  RSCS(:,:,:)                       !soil hydraulic resistance, [MPa h m-2]
  real(r8),target,allocatable ::  PSISE(:,:,:)                      !soil water potential at saturation, [Mpa]
  real(r8),target,allocatable ::  PSISA(:,:,:)                      !soil water potential at air entry, [Mpa]
  real(r8),target,allocatable ::  PSISO(:,:,:)                      !osmotic soil water potential , [Mpa]
  real(r8),target,allocatable ::  PSISH(:,:,:)                      !gravimetric soil water potential , [Mpa]
  real(r8),target,allocatable ::  THETY(:,:,:)                      !air-dry water content, [m3 m-3]
  real(r8),target,allocatable ::  THETS(:,:,:)                      !micropore class water content
  real(r8),target,allocatable ::  FLWX(:,:,:,:)                     !unsaturated water flux , [m3 d-2 h-1]
  real(r8),target,allocatable ::  UEVAP(:,:)                        !total evaporation, [m3 d-2]
  real(r8),target,allocatable ::  URAIN(:,:)                        !total precipitation, [m3 d-2]
  real(r8),target,allocatable ::  URUN(:,:)                         !total surface runoff, [m3 d-2]
  real(r8),target,allocatable ::  UVWatMicP(:,:)                        !total soil water content, [m3 d-2]
  real(r8),target,allocatable ::  UVOLO(:,:)                        !total subsurface water flux, [m3 d-2]
  real(r8),target,allocatable ::  UDRAIN(:,:)                       !total water drainage below root zone, [m3 d-2]
  real(r8),target,allocatable ::  QR(:,:,:,:)                       !soil surface runoff water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HQR(:,:,:,:)                      !soil surface runoff heat, [MJ d-2 h-1]
  real(r8),target,allocatable ::  WQRH(:,:)                         !runoff from surface water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HVOLO(:,:)                        !water discharge, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QRMN(:,:,:,:,:)                   !surface runoff,
  private :: InitAllocate
  contains

  subroutine InitSoilWater
  implicit none

  call InitAllocate
  end subroutine InitSoilWater

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(THETP(0:JZ,JY,JX));  THETP=0._r8
  allocate(VsoiP(0:JZ,JY,JX));   VsoiP=0._r8
  allocate(THETW(0:JZ,JY,JX));  THETW=0._r8
  allocate(THETI(0:JZ,JY,JX));  THETI=0._r8
  allocate(THETWZ(0:JZ,JY,JX)); THETWZ=0._r8
  allocate(THETIZ(0:JZ,JY,JX)); THETIZ=0._r8
  allocate(VWatMicP(0:JZ,JY,JX));   VWatMicP=0._r8
  allocate(ViceMicP(0:JZ,JY,JX));   ViceMicP=0._r8
  allocate(VWatMacP(JZ,JY,JX));    VWatMacP=0._r8
  allocate(PSISM(0:JZ,JY,JX));  PSISM=0._r8
  allocate(TotalSoilH2OPSIMPa(0:JZ,JY,JX));  TotalSoilH2OPSIMPa=0._r8
  allocate(VWatMicPX(0:JZ,JY,JX));  VWatMicPX=0._r8
  allocate(FINH(JZ,JY,JX));     FINH=0._r8
  allocate(ViceMacP(JZ,JY,JX));    ViceMacP=0._r8
  allocate(VWatMicPM(60,0:JZ,JY,JX));VWatMicPM=0._r8
  allocate(VWatMacPM(60,JZ,JY,JX));VWatMacPM=0._r8
  allocate(VsoiPM(60,0:JZ,JY,JX));VsoiPM=0._r8
  allocate(FILM(60,0:JZ,JY,JX));FILM=0._r8
  allocate(DTBLG(JY,JX));       DTBLG=0._r8
  allocate(DTBLDI(JY,JX));      DTBLDI=0._r8
  allocate(DTBLY(JY,JX));       DTBLY=0._r8
  allocate(DTBLD(JY,JX));       DTBLD=0._r8
  allocate(DPTHT(JY,JX));       DPTHT=0._r8
  allocate(DTBLZ(JY,JX));       DTBLZ=0._r8
  allocate(DTBLX(JY,JX));       DTBLX=0._r8
  allocate(DTBLI(JY,JX));       DTBLI=0._r8
  allocate(ENGYPM(60,JY,JX));   ENGYPM=0._r8
  allocate(XVOLTM(60,JY,JX));   XVOLTM=0._r8
  allocate(XVWatMicPM(60,JY,JX));   XVWatMicPM=0._r8
  allocate(XViceMicPM(60,JY,JX));   XViceMicPM=0._r8
  allocate(HCND(3,100,0:JZ,JY,JX));HCND=0._r8
  allocate(CNDH(JZ,JY,JX));     CNDH=0._r8
  allocate(CNDU(JZ,JY,JX));     CNDU=0._r8
  allocate(QRM(60,JV,JH));      QRM=0._r8
  allocate(QRV(60,JY,JX));      QRV=0._r8
  allocate(IFLBM(60,2,2,JY,JX));IFLBM=0
  allocate(IRCHG(2,2,JY,JX));   IRCHG=0
  allocate(IFLBH(2,2,JY,JX));   IFLBH=0
  allocate(RCHGNU(JY,JX));      RCHGNU=0._r8
  allocate(RCHGEU(JY,JX));      RCHGEU=0._r8
  allocate(RCHGSU(JY,JX));      RCHGSU=0._r8
  allocate(RCHGWU(JY,JX));      RCHGWU=0._r8
  allocate(RCHGNT(JY,JX));      RCHGNT=0._r8
  allocate(RCHGET(JY,JX));      RCHGET=0._r8
  allocate(RCHGST(JY,JX));      RCHGST=0._r8
  allocate(RCHGWT(JY,JX));      RCHGWT=0._r8
  allocate(RCHQN(JY,JX));       RCHQN=0._r8
  allocate(RCHQE(JY,JX));       RCHQE=0._r8
  allocate(RCHQS(JY,JX));       RCHQS=0._r8
  allocate(RCHQW(JY,JX));       RCHQW=0._r8
  allocate(RCHGD(JY,JX));       RCHGD=0._r8
  allocate(FLWM(60,3,JD,JV,JH));FLWM=0._r8
  allocate(WaterFlowMacPi(60,3,JD,JV,JH));WaterFlowMacPi=0._r8
  allocate(FLPM(60,JZ,JY,JX));  FLPM=0._r8
  allocate(FINHM(60,JZ,JY,JX)); FINHM=0._r8
  allocate(FLQSM(60,JY,JX));    FLQSM=0._r8
  allocate(FLQHM(60,JY,JX));    FLQHM=0._r8
  allocate(THETPM(60,0:JZ,JY,JX));THETPM=0._r8
  allocate(TORT(60,0:JZ,JY,JX));TORT=0._r8
  allocate(TORTH(60,JZ,JY,JX)); TORTH=0._r8
  allocate(DFGS(60,0:JZ,JY,JX));DFGS=0._r8
  allocate(RSCS(JZ,JY,JX));     RSCS=0._r8
  allocate(PSISE(0:JZ,JY,JX));  PSISE=0._r8
  allocate(PSISA(0:JZ,JY,JX));  PSISA=0._r8
  allocate(PSISO(0:JZ,JY,JX));  PSISO=0._r8
  allocate(PSISH(0:JZ,JY,JX));  PSISH=0._r8
  allocate(THETY(0:JZ,JY,JX));  THETY=0._r8
  allocate(THETS(0:JZ,JY,JX));  THETS=0._r8
  allocate(FLWX(3,JD,JV,JH));   FLWX=0._r8
  allocate(UEVAP(JY,JX));       UEVAP=0._r8
  allocate(URAIN(JY,JX));       URAIN=0._r8
  allocate(URUN(JY,JX));        URUN=0._r8
  allocate(UVWatMicP(JY,JX));       UVWatMicP=0._r8
  allocate(UVOLO(JY,JX));       UVOLO=0._r8
  allocate(UDRAIN(JY,JX));      UDRAIN=0._r8
  allocate(QR(2,2,JV,JH));      QR=0._r8
  allocate(HQR(2,2,JV,JH));     HQR=0._r8
  allocate(WQRH(JY,JX));        WQRH=0._r8
  allocate(HVOLO(JY,JX));       HVOLO=0._r8
  allocate(QRMN(60,2,2,JV,JH)); QRMN=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilWater
  use abortutils, only : destroy
  implicit none
  call destroy(THETP)
  call destroy(VsoiP)
  call destroy(THETW)
  call destroy(THETI)
  call destroy(THETWZ)
  call destroy(THETIZ)
  call destroy(VWatMicP)
  call destroy(ViceMicP)
  call destroy(VWatMacP)
  call destroy(PSISM)
  call destroy(TotalSoilH2OPSIMPa)
  call destroy(VWatMicPX)
  call destroy(FINH)
  call destroy(ViceMacP)
  call destroy(VWatMicPM)
  call destroy(VWatMacPM)
  call destroy(VsoiPM)
  call destroy(FILM)
  call destroy(DTBLG)
  call destroy(DTBLDI)
  call destroy(DTBLY)
  call destroy(DTBLD)
  call destroy(DPTHT)
  call destroy(DTBLZ)
  call destroy(DTBLX)
  call destroy(DTBLI)
  call destroy(ENGYPM)
  call destroy(XVOLTM)
  call destroy(XVWatMicPM)
  call destroy(XViceMicPM)
  call destroy(HCND)
  call destroy(CNDH)
  call destroy(CNDU)
  call destroy(QRM)
  call destroy(QRV)
  call destroy(IFLBM)
  call destroy(IRCHG)
  call destroy(IFLBH)
  call destroy(RCHGNU)
  call destroy(RCHGEU)
  call destroy(RCHGSU)
  call destroy(RCHGWU)
  call destroy(RCHGNT)
  call destroy(RCHGET)
  call destroy(RCHGST)
  call destroy(RCHGWT)
  call destroy(RCHQN)
  call destroy(RCHQE)
  call destroy(RCHQS)
  call destroy(RCHQW)
  call destroy(RCHGD)
  call destroy(FLWM)
  call destroy(WaterFlowMacPi)
  call destroy(FLPM)
  call destroy(FINHM)
  call destroy(FLQSM)
  call destroy(FLQHM)
  call destroy(THETPM)
  call destroy(TORT)
  call destroy(TORTH)
  call destroy(DFGS)
  call destroy(RSCS)
  call destroy(PSISE)
  call destroy(PSISA)
  call destroy(PSISO)
  call destroy(PSISH)
  call destroy(THETY)
  call destroy(THETS)
  call destroy(FLWX)
  call destroy(UEVAP)
  call destroy(URAIN)
  call destroy(URUN)
  call destroy(UVWatMicP)
  call destroy(UVOLO)
  call destroy(UDRAIN)
  call destroy(QR)
  call destroy(HQR)
  call destroy(WQRH)
  call destroy(HVOLO)
  call destroy(QRMN)
  end subroutine DestructSoilWater

end module SoilWaterDataType
