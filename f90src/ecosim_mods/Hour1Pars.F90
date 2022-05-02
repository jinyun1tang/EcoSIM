module Hour1Pars
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
!
!     *SG=diffusivity (m2 h-1):CG=CO2g,CL=CO2s,CH=CH4g,CQ=CH4s,OG=O2g
!     OL=O2s,ZG=N2g,ZL=N2s,Z2=N2Og,ZV=N2Os,ZH=NH3g,ZN=NH3s,ZO=NO3
!     PO=H2PO4,OC=DOC,ON=DON,OP=DOP,OA=acetate,WG=H2Og,AL=Al,FE=Fe
!     HY=H,CA=Ca,GM=Mg,AN=Na,AK=K,OH=OH,C3=CO3,HC=HCO3,SO=SO4,CL=Cl
!     HG=H2g,HL=H2s
!
!
!     SC*X=solubility (g m-3/g m-3):CO2=CO2,CH4=CH4,OXY=O2,N2G=N2g
!     N2O=N2O,NH3=NH3,H2G=H2
!     AC*X=activity (g m-3):CO2=CO2,CH4=CH4,OXY=O2,N2G=N2g
!     N2O=N2O,NH3=NH3,H2G=H2
!
!     VISCW=water viscosity (Mg m-1 s)
!     BKDSX=maximm soil bulk density
!     THETPW=minimum air-filled porosity for saturation (m3 m-3)
!

  real(r8) :: CGSG
  real(r8) :: CLSG
  real(r8) :: CHSG
  real(r8) :: CQSG
  real(r8) :: OGSG
  real(r8) :: OLSG
  real(r8) :: ZGSG
  real(r8) :: ZLSG
  real(r8) :: Z2SG
  real(r8) :: ZVSG
  real(r8) :: ZHSG
  real(r8) :: ZNSG
  real(r8) :: ZOSG
  real(r8) :: POSG
  real(r8) :: OCSG
  real(r8) :: ONSG
  real(r8) :: OPSG
  real(r8) :: OASG
  real(r8) :: WGSG
  real(r8) :: ALSG
  real(r8) :: FESG
  real(r8) :: HYSG
  real(r8) :: CASG
  real(r8) :: GMSG
  real(r8) :: ANSG
  real(r8) :: AKSG
  real(r8) :: OHSG
  real(r8) :: C3SG
  real(r8) :: HCSG
  real(r8) :: SOSG
  real(r8) :: CLSX
  real(r8) :: HGSG
  real(r8) :: HLSG
  real(r8) :: SCO2X
  real(r8) :: SCH4X
  real(r8) :: SOXYX
  real(r8) :: SN2GX
  real(r8) :: SN2OX
  real(r8) :: SNH3X
  real(r8) :: SH2GX
  real(r8) :: ACO2X
  real(r8) :: ACH4X
  real(r8) :: AOXYX
  real(r8) :: AN2GX
  real(r8) :: AN2OX
  real(r8) :: ANH3X
  real(r8) :: AH2GX
  real(r8) :: VISCW
  real(r8) :: BKDSX
  real(r8) :: FORGW
  real(r8) :: DTHETW
  real(r8) :: THETPW
  real(r8) :: THETWP
  real(r8) :: RAM
  real(r8) :: XVOLWC(0:3),THETRX(0:2)
!
!     XVOLWC=foliar water retention capacity (m3 m-2)
!     THETRX=litter water retention capacity (m3 g C-1)
!

  contains

  subroutine initHour1Pars
  implicit none
  CGSG=4.68E-02
  CLSG=4.25E-06
  CHSG=7.80E-02
  CQSG=7.08E-06
  OGSG=6.43E-02
  OLSG=8.57E-06
  ZGSG=5.57E-02
  ZLSG=7.34E-06
  Z2SG=5.57E-02
  ZVSG=5.72E-06
  ZHSG=6.67E-02
  ZNSG=4.00E-06
  ZOSG=6.00E-06
  POSG=3.00E-06
  OCSG=1.0E-08
  ONSG=1.0E-08
  OPSG=1.0E-08
  OASG=3.64E-06
  WGSG=7.70E-02
  ALSG=5.0E-06
  FESG=5.0E-06
  HYSG=5.0E-06
  CASG=5.0E-06
  GMSG=5.0E-06
  ANSG=5.0E-06
  AKSG=5.0E-06
  OHSG=5.0E-06
  C3SG=5.0E-06
  HCSG=5.0E-06
  SOSG=5.0E-06
  CLSX=5.0E-06
  HGSG=5.57E-02
  HLSG=7.34E-06
  SCO2X=7.391E-01
  SCH4X=3.156E-02
  SOXYX=2.925E-02
  SN2GX=1.510E-02
  SN2OX=5.241E-01
  SNH3X=2.852E+02
  SH2GX=3.156E-02
  ACO2X=0.14
  ACH4X=0.14
  AOXYX=0.31
  AN2GX=0.23
  AN2OX=0.23
  ANH3X=0.07
  AH2GX=0.14
  VISCW=1.0E-06
  BKDSX=1.89
  FORGW=0.25E+06
  DTHETW=1.0E-06
  THETPW=0.01
  THETWP=1.0-THETPW

  XVOLWC=(/5.0E-04,2.5E-04,2.5E-04,2.5E-04/)
  THETRX=(/4.0E-06,8.0E-06,8.0E-06/)

  end subroutine initHour1Pars

end module Hour1Pars
