module SoilWaterDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  THETP(:,:,:)                      !air concentration [m3 m-3]
  real(r8),allocatable ::  VOLP(:,:,:)                       !soil air content [m3 d-2]
  real(r8),allocatable ::  THETW(:,:,:)                      !volumetric water content [m3 m-3]
  real(r8),allocatable ::  THETI(:,:,:)                      !volumetric ice content [m3 m-3]
  real(r8),allocatable ::  THETWZ(:,:,:)                     !volumetric moblize water [m3 m-3]
  real(r8),allocatable ::  THETIZ(:,:,:)                     !volumetric mobile ice [m3 m-3]
  real(r8),allocatable ::  VOLW(:,:,:)                       !soil micropore water content [m3 d-2]
  real(r8),allocatable ::  VOLI(:,:,:)                       !soil micropore ice content   [m3 d-2]
  real(r8),allocatable ::  VOLWH(:,:,:)                      !soil macropore water content [m3 d-2]
  real(r8),allocatable ::  PSISM(:,:,:)                      !soil micropore matric water potential [MPa]
  real(r8),allocatable ::  PSIST(:,:,:)                      !soil micropore total water potential [MPa]
  real(r8),allocatable ::  VOLWX(:,:,:)                      !soil micropore water content before wetting front [m3 d-2]
  real(r8),allocatable ::  FINH(:,:,:)                       !soil macropore - micropore water transfer [m3 d-2 h-1]
  real(r8),allocatable ::  VOLIH(:,:,:)                      !soil macropore ice content [m3 d-2]
  real(r8),allocatable ::  VOLWM(:,:,:,:)                    !soil micropore water content, [m3 d-2]
  real(r8),allocatable ::  VOLWHM(:,:,:,:)                   !soil macropore water content, [m3 d-2]
  real(r8),allocatable ::  VOLPM(:,:,:,:)                    !soil air content, [m3 d-2]
  real(r8),allocatable ::  FILM(:,:,:,:)                     !soil water film thickness , [m]
  real(r8),allocatable ::  DTBLG(:,:)                        !slope of water table relative to surface slope, [-]
  real(r8),allocatable ::  DTBLDI(:,:)                       !depth of artificial water table
  real(r8),allocatable ::  DTBLY(:,:)                        !artificial water table depth, [m]
  real(r8),allocatable ::  DTBLD(:,:)                        !depth of artificial water table adjusted for elevation
  real(r8),allocatable ::  DPTHT(:,:)                        !internal water table depth, [m]
  real(r8),allocatable ::  DTBLZ(:,:)                        !external water table depth, [m]
  real(r8),allocatable ::  DTBLX(:,:)                        !external water table depth, [m]
  real(r8),allocatable ::  DTBLI(:,:)                        !external water table depth, [m]
  real(r8),allocatable ::  ENGYPM(:,:,:)                     !total energy impact for erosion
  real(r8),allocatable ::  XVOLTM(:,:,:)                     !excess water+ice
  real(r8),allocatable ::  XVOLWM(:,:,:)                     !excess water
  real(r8),allocatable ::  XVOLIM(:,:,:)                     !excess ice
  real(r8),allocatable ::  HCND(:,:,:,:,:)                   !saturated hydraulic conductivity
  real(r8),allocatable ::  CNDH(:,:,:)                       !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),allocatable ::  CNDU(:,:,:)                       !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8),allocatable ::  QRM(:,:,:)                        !runoff water flux, [m3 d-2 t-1]
  real(r8),allocatable ::  QRV(:,:,:)                        !runoff velocity, [m t-1]
   integer,allocatable ::  IFLBM(:,:,:,:,:)                   !flag for directional surface runoff
   integer,allocatable ::  IRCHG(:,:,:,:)                     !enables or disables boundary water flux depending on aspect, [-]
   integer,allocatable ::  IFLBH(:,:,:,:)                     !flag for directional runoff, related to IFLBM
  real(r8),allocatable ::  RCHGNU(:,:)                       !northern subsurface boundary water flux , [-]
  real(r8),allocatable ::  RCHGEU(:,:)                       !eastern subsurface boundary water flux , [-]
  real(r8),allocatable ::  RCHGSU(:,:)                       !southern subsurface boundary water flux , [-]
  real(r8),allocatable ::  RCHGWU(:,:)                       !western subsurface boundary water flux , [-]
  real(r8),allocatable ::  RCHGNT(:,:)                       !northern subsurface boundary water flux rate constant, [h-1]
  real(r8),allocatable ::  RCHGET(:,:)                       !eastern subsurface boundary water flux  rate constant, [h-1]
  real(r8),allocatable ::  RCHGST(:,:)                       !southern subsurface boundary water flux  rate constant, [h-1]
  real(r8),allocatable ::  RCHGWT(:,:)                       !western subsurface boundary water flux  rate constant, [h-1]
  real(r8),allocatable ::  RCHQN(:,:)                        !northern surface boundary water flux , [-]
  real(r8),allocatable ::  RCHQE(:,:)                        !eastern surface boundary water flux , [-]
  real(r8),allocatable ::  RCHQS(:,:)                        !southern surface boundary water flux , [-]
  real(r8),allocatable ::  RCHQW(:,:)                        !western surface boundary water flux , [-]
  real(r8),allocatable ::  RCHGD(:,:)                        !lower subsurface boundary water flux , [-]
  real(r8),allocatable ::  FLWM(:,:,:,:,:)                   !micropore water flux, [m3 d-2 t-1]
  real(r8),allocatable ::  FLWHM(:,:,:,:,:)                  !macropore water flux, [m3 d-2 t-1]
  real(r8),allocatable ::  FLPM(:,:,:,:)                     !soil air flux, [g d-2 t-1]
  real(r8),allocatable ::  FINHM(:,:,:,:)                    !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8),allocatable ::  FLQSM(:,:,:)                      !meltwater flux into soil micropores
  real(r8),allocatable ::  FLQHM(:,:,:)                      !meltwater flux into soil macropores
  real(r8),allocatable ::  THETPM(:,:,:,:)                   !soil air-filled porosity, [m3 m-3]
  real(r8),allocatable ::  TORT(:,:,:,:)                     !soil tortuosity, []
  real(r8),allocatable ::  TORTH(:,:,:,:)                    !macropore tortuosity, []
  real(r8),allocatable ::  DFGS(:,:,:,:)                     !coefficient for dissolution - volatilization, []
  real(r8),allocatable ::  RSCS(:,:,:)                       !soil hydraulic resistance, [MPa h m-2]
  real(r8),allocatable ::  PSISE(:,:,:)                      !soil water potential at saturation, [Mpa]
  real(r8),allocatable ::  PSISA(:,:,:)                      !soil water potential at air entry, [Mpa]
  real(r8),allocatable ::  PSISO(:,:,:)                      !osmotic soil water potential , [Mpa]
  real(r8),allocatable ::  PSISH(:,:,:)                      !gravimetric soil water potential , [Mpa]
  real(r8),allocatable ::  THETY(:,:,:)                      !air-dry water content, [m3 m-3]
  real(r8),allocatable ::  THETS(:,:,:)                      !micropore class water content
  real(r8),allocatable ::  FLWX(:,:,:,:)                     !unsaturated water flux , [m3 d-2 h-1]
  real(r8),allocatable ::  UEVAP(:,:)                        !total evaporation, [m3 d-2]
  real(r8),allocatable ::  URAIN(:,:)                        !total precipitation, [m3 d-2]
  real(r8),allocatable ::  URUN(:,:)                         !total surface runoff, [m3 d-2]
  real(r8),allocatable ::  UVOLW(:,:)                        !total soil water content, [m3 d-2]
  real(r8),allocatable ::  UVOLO(:,:)                        !total subsurface water flux, [m3 d-2]
  real(r8),allocatable ::  UDRAIN(:,:)                       !total water drainage below root zone, [m3 d-2]
  real(r8),allocatable ::  QR(:,:,:,:)                       !soil surface runoff water, [m3 d-2 h-1]
  real(r8),allocatable ::  HQR(:,:,:,:)                      !soil surface runoff heat, [MJ d-2 h-1]
  real(r8),allocatable ::  WQRH(:,:)                         !runoff from surface water, [m3 d-2 h-1]
  real(r8),allocatable ::  HVOLO(:,:)                        !water discharge, [m3 d-2 h-1]
  real(r8),allocatable ::  QRMN(:,:,:,:,:)                   !surface runoff,
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
  allocate(VOLP(0:JZ,JY,JX));   VOLP=0._r8
  allocate(THETW(0:JZ,JY,JX));  THETW=0._r8
  allocate(THETI(0:JZ,JY,JX));  THETI=0._r8
  allocate(THETWZ(0:JZ,JY,JX)); THETWZ=0._r8
  allocate(THETIZ(0:JZ,JY,JX)); THETIZ=0._r8
  allocate(VOLW(0:JZ,JY,JX));   VOLW=0._r8
  allocate(VOLI(0:JZ,JY,JX));   VOLI=0._r8
  allocate(VOLWH(JZ,JY,JX));    VOLWH=0._r8
  allocate(PSISM(0:JZ,JY,JX));  PSISM=0._r8
  allocate(PSIST(0:JZ,JY,JX));  PSIST=0._r8
  allocate(VOLWX(0:JZ,JY,JX));  VOLWX=0._r8
  allocate(FINH(JZ,JY,JX));     FINH=0._r8
  allocate(VOLIH(JZ,JY,JX));    VOLIH=0._r8
  allocate(VOLWM(60,0:JZ,JY,JX));VOLWM=0._r8
  allocate(VOLWHM(60,JZ,JY,JX));VOLWHM=0._r8
  allocate(VOLPM(60,0:JZ,JY,JX));VOLPM=0._r8
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
  allocate(XVOLWM(60,JY,JX));   XVOLWM=0._r8
  allocate(XVOLIM(60,JY,JX));   XVOLIM=0._r8
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
  allocate(FLWHM(60,3,JD,JV,JH));FLWHM=0._r8
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
  allocate(UVOLW(JY,JX));       UVOLW=0._r8
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
  if (allocated(THETP))    deallocate(THETP)
  if (allocated(VOLP))     deallocate(VOLP)
  if (allocated(THETW))    deallocate(THETW)
  if (allocated(THETI))    deallocate(THETI)
  if (allocated(THETWZ))   deallocate(THETWZ)
  if (allocated(THETIZ))   deallocate(THETIZ)
  if (allocated(VOLW))     deallocate(VOLW)
  if (allocated(VOLI))     deallocate(VOLI)
  if (allocated(VOLWH))    deallocate(VOLWH)
  if (allocated(PSISM))    deallocate(PSISM)
  if (allocated(PSIST))    deallocate(PSIST)
  if (allocated(VOLWX))    deallocate(VOLWX)
  if (allocated(FINH))     deallocate(FINH)
  if (allocated(VOLIH))    deallocate(VOLIH)
  if (allocated(VOLWM))    deallocate(VOLWM)
  if (allocated(VOLWHM))   deallocate(VOLWHM)
  if (allocated(VOLPM))    deallocate(VOLPM)
  if (allocated(FILM))     deallocate(FILM)
  if (allocated(DTBLG))    deallocate(DTBLG)
  if (allocated(DTBLDI))   deallocate(DTBLDI)
  if (allocated(DTBLY))    deallocate(DTBLY)
  if (allocated(DTBLD))    deallocate(DTBLD)
  if (allocated(DPTHT))    deallocate(DPTHT)
  if (allocated(DTBLZ))    deallocate(DTBLZ)
  if (allocated(DTBLX))    deallocate(DTBLX)
  if (allocated(DTBLI))    deallocate(DTBLI)
  if (allocated(ENGYPM))   deallocate(ENGYPM)
  if (allocated(XVOLTM))   deallocate(XVOLTM)
  if (allocated(XVOLWM))   deallocate(XVOLWM)
  if (allocated(XVOLIM))   deallocate(XVOLIM)
  if (allocated(HCND))     deallocate(HCND)
  if (allocated(CNDH))     deallocate(CNDH)
  if (allocated(CNDU))     deallocate(CNDU)
  if (allocated(QRM))      deallocate(QRM)
  if (allocated(QRV))      deallocate(QRV)
  if (allocated(IFLBM))    deallocate(IFLBM)
  if (allocated(IRCHG))    deallocate(IRCHG)
  if (allocated(IFLBH))    deallocate(IFLBH)
  if (allocated(RCHGNU))   deallocate(RCHGNU)
  if (allocated(RCHGEU))   deallocate(RCHGEU)
  if (allocated(RCHGSU))   deallocate(RCHGSU)
  if (allocated(RCHGWU))   deallocate(RCHGWU)
  if (allocated(RCHGNT))   deallocate(RCHGNT)
  if (allocated(RCHGET))   deallocate(RCHGET)
  if (allocated(RCHGST))   deallocate(RCHGST)
  if (allocated(RCHGWT))   deallocate(RCHGWT)
  if (allocated(RCHQN))    deallocate(RCHQN)
  if (allocated(RCHQE))    deallocate(RCHQE)
  if (allocated(RCHQS))    deallocate(RCHQS)
  if (allocated(RCHQW))    deallocate(RCHQW)
  if (allocated(RCHGD))    deallocate(RCHGD)
  if (allocated(FLWM))     deallocate(FLWM)
  if (allocated(FLWHM))    deallocate(FLWHM)
  if (allocated(FLPM))     deallocate(FLPM)
  if (allocated(FINHM))    deallocate(FINHM)
  if (allocated(FLQSM))    deallocate(FLQSM)
  if (allocated(FLQHM))    deallocate(FLQHM)
  if (allocated(THETPM))   deallocate(THETPM)
  if (allocated(TORT))     deallocate(TORT)
  if (allocated(TORTH))    deallocate(TORTH)
  if (allocated(DFGS))     deallocate(DFGS)
  if (allocated(RSCS))     deallocate(RSCS)
  if (allocated(PSISE))    deallocate(PSISE)
  if (allocated(PSISA))    deallocate(PSISA)
  if (allocated(PSISO))    deallocate(PSISO)
  if (allocated(PSISH))    deallocate(PSISH)
  if (allocated(THETY))    deallocate(THETY)
  if (allocated(THETS))    deallocate(THETS)
  if (allocated(FLWX))     deallocate(FLWX)
  if (allocated(UEVAP))    deallocate(UEVAP)
  if (allocated(URAIN))    deallocate(URAIN)
  if (allocated(URUN))     deallocate(URUN)
  if (allocated(UVOLW))    deallocate(UVOLW)
  if (allocated(UVOLO))    deallocate(UVOLO)
  if (allocated(UDRAIN))   deallocate(UDRAIN)
  if (allocated(QR))       deallocate(QR)
  if (allocated(HQR))      deallocate(HQR)
  if (allocated(WQRH))     deallocate(WQRH)
  if (allocated(HVOLO))    deallocate(HVOLO)
  if (allocated(QRMN))     deallocate(QRMN)
  end subroutine DestructSoilWater

end module SoilWaterDataType
