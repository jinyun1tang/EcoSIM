module SoilWaterDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save

  real(r8) :: FSNW(JY,JX)                       !fraction of snow cover
  REAL(R8) :: THETP(0:JZ,JY,JX)       !air concentration [m3 m-3]
  real(r8) :: VOLP(0:JZ,JY,JX)        !soil air content [m3 d-2]

  real(r8) :: HCND(3,100,0:JZ,JY,JX)  !saturated hydraulic conductivity
  real(r8) :: THETW(0:JZ,JY,JX)       !volumetric water content [m3 m-3]
  real(r8) :: THETI(0:JZ,JY,JX)       !volumetric ice content [m3 m-3]
  real(r8) :: THETWZ(0:JZ,JY,JX)
  real(r8) :: THETIZ(0:JZ,JY,JX)
  real(r8) :: VOLW(0:JZ,JY,JX)        !soil micropore water content [m3 d-2]
  real(r8) :: VOLI(0:JZ,JY,JX)        !soil micropore ice content   [m3 d-2]
  real(r8) :: VOLWH(JZ,JY,JX)         !soil macropore water content [m3 d-2]
  real(r8) :: PSISM(0:JZ,JY,JX)       !soil micropore matric water potential [MPa]
  real(r8) :: PSIST(0:JZ,JY,JX)       !soil micropore total water potential [MPa]
  real(r8) :: CNDH(JZ,JY,JX)          !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8) :: CNDU(JZ,JY,JX)          !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8) :: VOLWX(0:JZ,JY,JX)       !soil micropore water content before wetting front [m3 d-2]
  real(r8) :: FINH(JZ,JY,JX)          !soil macropore - micropore water transfer [m3 d-2 h-1]
  real(r8) :: VOLIH(JZ,JY,JX)         !soil macropore ice content [m3 d-2]
  real(r8) :: VOLWM(60,0:JZ,JY,JX)              !soil micropore water content, [m3 d-2]
  real(r8) :: VOLWHM(60,JZ,JY,JX)               !soil macropore water content, [m3 d-2]
  real(r8) :: VOLPM(60,0:JZ,JY,JX)              !soil air content, [m3 d-2]
  real(r8) :: FILM(60,0:JZ,JY,JX)               !soil water film thickness , [m]

  real(r8) :: QRM(60,JV,JH)                     !runoff water flux, [m3 d-2 t-1]
  real(r8) :: QRV(60,JY,JX)                     !runoff velocity, [m t-1]
  integer  :: IFLBM(60,2,2,JY,JX)
  integer  :: IRCHG(2,2,JY,JX)                  !enables or disables boundary water flux depending on aspect, [-]
  integer  :: IFLBH(2,2,JY,JX)
  real(r8) :: RCHGNU(JY,JX)                     !northern subsurface boundary water flux , [-]
  real(r8) :: RCHGEU(JY,JX)                     !eastern subsurface boundary water flux , [-]
  real(r8) :: RCHGSU(JY,JX)                     !southern subsurface boundary water flux , [-]
  real(r8) :: RCHGWU(JY,JX)                     !western subsurface boundary water flux , [-]
  real(r8) :: RCHGNT(JY,JX)                     !northern subsurface boundary water flux rate constant, [h-1]
  real(r8) :: RCHGET(JY,JX)                     !eastern subsurface boundary water flux  rate constant, [h-1]
  real(r8) :: RCHGST(JY,JX)                     !southern subsurface boundary water flux  rate constant, [h-1]
  real(r8) :: RCHGWT(JY,JX)                     !western subsurface boundary water flux  rate constant, [h-1]
  real(r8) :: RCHQN(JY,JX)                      !northern surface boundary water flux , [-]
  real(r8) :: RCHQE(JY,JX)                      !eastern surface boundary water flux , [-]
  real(r8) :: RCHQS(JY,JX)                      !southern surface boundary water flux , [-]
  real(r8) :: RCHQW(JY,JX)                      !western surface boundary water flux , [-]
  real(r8) :: RCHGD(JY,JX)                      !lower subsurface boundary water flux , [-]
  real(r8) :: DTBLG(JY,JX)                      !slope of water table relative to surface slope, [-]
  real(r8) :: DTBLDI(JY,JX)                     !depth of artificial water table
  real(r8) :: DTBLY(JY,JX)                      !artificial water table depth, [m]
  real(r8) :: DTBLD(JY,JX)                      !depth of artificial water table adjusted for elevation
  real(r8) :: DPTHT(JY,JX)                      !internal water table depth, [m]
  real(r8) :: DTBLZ(JY,JX)                      !external water table depth, [m]
  real(r8) :: DTBLX(JY,JX)                      !external water table depth, [m]
  real(r8) :: DTBLI(JY,JX)                      !external water table depth, [m]

  real(r8) :: FLWM(60,3,JD,JV,JH)               !micropore water flux, [m3 d-2 t-1]
  real(r8) :: FLWHM(60,3,JD,JV,JH)              !macropore water flux, [m3 d-2 t-1]
  real(r8) :: QRMN(60,2,2,JV,JH)                !runoff
  real(r8) :: QSM(60,2,JV,JH)                   !runoff snow flux, [m3 d-2 t-1]
  real(r8) :: ENGYPM(60,JY,JX)                  !total energy impact for erosion
  real(r8) :: XVOLTM(60,JY,JX)                  !excess water+ice
  real(r8) :: XVOLWM(60,JY,JX)                  !excess water
  real(r8) :: XVOLIM(60,JY,JX)                  !excess ice
  real(r8) :: FLPM(60,JZ,JY,JX)                 !soil air flux, [g d-2 t-1]
  real(r8) :: FINHM(60,JZ,JY,JX)                !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8) :: FLQSM(60,JY,JX)                   !meltwater flux into soil micropores
  real(r8) :: FLQHM(60,JY,JX)                   !meltwater flux into soil macropores

  real(r8) :: ROXSK(60,0:JZ,JY,JX)              !total O2 sink, [g d-2 t-1]
  real(r8) :: THETPM(60,0:JZ,JY,JX)             !soil air-filled porosity, [m3 m-3]
  real(r8) :: TORT(60,0:JZ,JY,JX)               !soil tortuosity, []
  real(r8) :: TORTH(60,JZ,JY,JX)                !macropore tortuosity, []
  real(r8) :: DFGS(60,0:JZ,JY,JX)               ! coefficient for dissolution - volatilization, []
  real(r8) :: RSCS(JZ,JY,JX)                    !soil hydraulic resistance, [MPa h m-2]

  real(r8) :: PSISE(0:JZ,JY,JX)                 !soil water potential at saturation, [Mpa]
  real(r8) :: PSISA(0:JZ,JY,JX)                 !soil water potential at air entry, [Mpa]
  real(r8) :: PSISO(0:JZ,JY,JX)                 !osmotic soil water potential , [Mpa]
  real(r8) :: PSISH(0:JZ,JY,JX)                 !gravimetric soil water potential , [Mpa]
  real(r8) :: THETY(0:JZ,JY,JX)                 !air-dry water content, [m3 m-3]
  real(r8) :: THETS(0:JZ,JY,JX)                 !micropore class water content
  real(r8) :: FLWX(3,JD,JV,JH)                  !unsaturated water flux , [m3 d-2 h-1]
  real(r8) :: UEVAP(JY,JX)                      !total evaporation, [m3 d-2]
  real(r8) :: URAIN(JY,JX)                      !total precipitation, [m3 d-2]
  real(r8) :: URUN(JY,JX)                       !total surface runoff, [m3 d-2]
  real(r8) :: UVOLW(JY,JX)                      !total soil water content, [m3 d-2]
  real(r8) :: UVOLO(JY,JX)                      !total subsurface water flux, [m3 d-2]
  real(r8) :: UDRAIN(JY,JX)                     !total water drainage below root zone, [m3 d-2]
  real(r8) :: WQRH(JY,JX)                       !water flux
  real(r8) :: HVOLO(JY,JX)                      !water discharge, [m3 d-2 h-1]

end module SoilWaterDataType
