module SurfLitterPhysMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use MiniMathMod
  use HydroThermData
  use SurfSoilDataType
  use GridDataType
  use SurfLitterDataType
  use ChemTranspDataType
  use SoilWaterDataType
  use SoilPropertyDataType
  use SurfPhysData
  use SoilPhysDataType
  use EcoSIMCtrlDataType
  use SOMDataType
  use LandSurfDataType
  use EcoSimConst
  use EcoSIMSolverPar
  use abortutils
  use PhysPars
  use UnitMod, only : units  
  use SoilPhysParaMod
  USE SoilHeatDataType
  USE ClimForcDataType
  use SnowPhysData
implicit none
  character(len=*), parameter, private :: mod_filename=&
  __FILE__
  public :: SRFLitterEnergyBalance
  public :: UpdateLitRPhys
  public :: UpdateLitRB4RunoffM
  public :: UpdateLitRAftRunoff
  contains


!------------------------------------------------------------------------------------------
  subroutine SRFLitterEnergyBalance(M,NY,NX,PSISV1,Prec2LitR2,PrecHeat2LitR2,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR)
!!
! Description:
  implicit none
  integer , intent(in) :: M       !soil heat-water iteration id
  integer , intent(in) :: NY,NX
  real(r8), intent(in) :: PSISV1
  real(r8), intent(in) :: Prec2LitR2,PrecHeat2LitR2
  real(r8), intent(out) :: HeatSensLitR2Soi1
  real(r8), intent(out) :: HeatSensVapLitR2Soi1
  real(r8), intent(out) :: EvapLitR2Soi1
  real(r8), intent(out) :: TotHeatAir2LitR
  real(r8) :: HeatSensAir2LitR
  real(r8) :: HeatSensEvapAir2LitR
  real(r8) :: LWRadLitR  
  real(r8) :: VWatLitr2,ATCNDR,ATCNVR,Radt2LitR
  real(r8) :: CNVR,CNV1
  real(r8) :: TCNDR
  real(r8) :: TCND1
  real(r8) :: DTKX
  real(r8) :: LitRAlbedo  
  real(r8) :: Radnet2LitR,LatentHeatAir2LitR

! begin_execution
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!

  VapXAir2LitR(NY,NX)  = 0.0_r8
  Radnet2LitR          = 0.0_r8
  LatentHeatAir2LitR   = 0.0_r8
  HeatSensEvapAir2LitR = 0.0_r8
  HeatSensAir2LitR     = 0.0_r8
  TotHeatAir2LitR      = 0.0_r8

  EvapLitR2Soi1        = 0.0_r8
  HeatSensVapLitR2Soi1 = 0.0_r8
  HeatSensLitR2Soi1    = 0.0_r8
  LWRadLitR            = 0.0_r8

  IF(VHeatCapacity1_vr(0,NY,NX).LE.VHeatCapLitRMin_col(NY,NX))THEN
    TKSoi1_vr(0,NY,NX)=TKSoi1_vr(NUM(NY,NX),NY,NX)
    RETURN
  ENDIF  
!
! NET RADIATION AT RESIDUE SURFACE
!
! LitRAlbedo=litter albedo
! VLWatMicP1,VLiceMicP1=water,ice volume in litter
! RadSWonLitR,LWRad2LitR=incoming shortwave,longwave radiation

!
! litter layer
! albedo
  LitRAlbedo=(0.20_r8*VLSoilMicPMass_vr(0,NY,NX)+0.06_r8*VLWatMicP1_vr(0,NY,NX)+0.30_r8 &
    *VLiceMicP1_vr(0,NY,NX))/(VLSoilMicPMass_vr(0,NY,NX)+VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX))

  !radiation incident on litter layer  
  Radt2LitR=(1.0_r8-LitRAlbedo)*RadSWonLitR(NY,NX)+LWRad2LitR(NY,NX)  

  !
  ! THERMAL CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
  !
  ! CNVR,CNV1=litter,soil vapor conductivity
  ! THETPM=litter air concentration
  ! POROS,POROQ=litter porosity, tortuosity
  ! VaporDiffusivityLitR_col,WVapDifusvitySoil_vr=litter,soil vapor diffusivity
  ! CVRD=litter cover fraction
  ! ATCNVR=litter-soil vapor conductance
  ! DLYRR,DLYR=litter,soil depths
  ! THETRR=dry litter concentration
  ! DTH*,RYL*,DNU*,TRB*=turbulence effects on thermal conductivity
  ! WTHET0,WTHET1=multiplier for air concn in thermal conductivity
  ! TCNDW*,TCNDA*=thermal conductivity of water,air
  ! TCNDR,TCND1=litter,soil thermal conductivity
  ! ATCNDR=litter-soil thermal conductance
  !
  CNVR=VaporDiffusivityLitR_col(NY,NX)*THETPM(M,0,NY,NX)*POROQ*THETPM(M,0,NY,NX)/POROS_vr(0,NY,NX)
  CNV1=WVapDifusvitySoil_vr(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ &
    *THETPM(M,NUM(NY,NX),NY,NX)/POROS_vr(NUM(NY,NX),NY,NX)

  IF(FracSurfByLitR_col(NY,NX).GT.ZERO)THEN
    !there is litter layer
    IF(CNVR.GT.ZERO.AND.CNV1.GT.ZERO)THEN
      ATCNVR=2.0_r8*CNVR*CNV1/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV1*DLYRR_COL(NY,NX))
    ELSE
      !below is a numerical hack
      ATCNVR=2.0_r8*CNVR/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR_COL(NY,NX))*FracSurfByLitR_col(NY,NX)
    ENDIF
  ELSE
    ATCNVR=0.0_r8
  ENDIF


  DTKX=ABS(TKSoi1_vr(0,NY,NX)-TKSoi1_vr(NUM(NY,NX),NY,NX))*ppmc

  call CalcLitRThermConductivity(NY,NX,DTKX,TCNDR)

  call CalcSoilThermConductivity(NX,NY,NUM(NY,NX),DTKX,TCND1)

  ATCNDR=2.0_r8*TCNDR*TCND1/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)+TCND1*DLYRR_COL(NY,NX))
!
! SMALLER TIME STEP FOR SOLVING SURFACE RESIDUE ENERGY EXCHANGE
!
  call SurfLitterIterateNN(M,NY,NX,ATCNDR,ATCNVR,PSISV1,Radt2LitR,Prec2LitR2,PrecHeat2LitR2,&
    EvapLitR2Soi1,HeatSensAir2LitR,HeatSensEvapAir2LitR,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,LatentHeatAir2LitR,LWRadLitR,Radnet2LitR,TotHeatAir2LitR)

  HeatByRadiation_col(NY,NX)=HeatByRadiation_col(NY,NX)+Radnet2LitR
  HeatNet2Surf_col(NY,NX)=HeatNet2Surf_col(NY,NX)+Radnet2LitR+HeatSensEvapAir2LitR+HeatSensAir2LitR
  HeatEvapAir2Surf_col(NY,NX)=HeatEvapAir2Surf_col(NY,NX)+LatentHeatAir2LitR
  HeatSensVapAir2Surf_col(NY,NX)=HeatSensVapAir2Surf_col(NY,NX)+HeatSensEvapAir2LitR
  LWRadBySurf_col(NY,NX)=LWRadBySurf_col(NY,NX)+LWRadLitR

  end subroutine SRFLitterEnergyBalance
!------------------------------------------------------------------------------------------

  subroutine CalcLitRThermConductivity(NY,NX,DTKX,TCNDR)

  use PhysPars
  implicit none
  integer , intent(in) :: NY,NX
  real(r8), intent(in) :: DTKX    !absolute temperautre gradient between litter and soil
  real(r8), intent(out) :: TCNDR
  real(r8) :: HeatDiffusByWat0,HeatDiffusByAir0
  real(r8) :: RYLXW0,RYLXA0,RYLNA0,RYLNW0
  real(r8) :: XNUSW0,XNUSA0
  real(r8) :: TCNDW0,TCNDA0,WTHET0,THETRR

  !THETRR=litter in relative volume
  !FracSoilAsAirt=relative volume as air  
  THETRR=AZMAX1(1.0_r8-FracSoiPAsAir_vr(0,NY,NX)-FracSoiPAsWat_vr(0,NY,NX)-FracSoiPAsIce_vr(0,NY,NX))

  HeatDiffusByWat0=AZMAX1(FracSoiPAsWat_vr(0,NY,NX)-TRBW)**3._r8
  HeatDiffusByAir0=AZMAX1(FracSoiPAsAir_vr(0,NY,NX)-TRBA)**3._r8
  RYLXW0=DTKX*HeatDiffusByWat0
  RYLXA0=DTKX*HeatDiffusByAir0
  RYLNW0=AMIN1(1.0E+04_r8,RYLXW*RYLXW0)
  RYLNA0=AMIN1(1.0E+04_r8,RYLXA*RYLXA0)
  XNUSW0=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW0**0.25_r8/DNUSW)
  XNUSA0=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA0**0.25_r8/DNUSA)
  TCNDW0=2.067E-03_r8*XNUSW0
  TCNDA0=9.050E-05_r8*XNUSA0  
  WTHET0=1.467_r8-0.467_r8*FracSoilAsAirt(0,NY,NX)
  TCNDR=(0.779_r8*THETRR*9.050E-04_r8+0.622_r8*FracSoiPAsWat_vr(0,NY,NX)*TCNDW0 &
    +0.380_r8*FracSoiPAsIce_vr(0,NY,NX)*7.844E-03_r8+WTHET0*FracSoiPAsAir_vr(0,NY,NX)*TCNDA0) &
    /(0.779_r8*THETRR+0.622_r8*FracSoiPAsWat_vr(0,NY,NX)+0.380_r8*FracSoiPAsIce_vr(0,NY,NX)+WTHET0*FracSoiPAsAir_vr(0,NY,NX))
  end subroutine CalcLitRThermConductivity

!------------------------------------------------------------------------------------------

  subroutine SurfLitterIterateNN(M,NY,NX,ATCNDR,ATCNVR,PSISV1,Radt2LitR,Prec2LitR2,PrecHeat2LitR2,&
    EvapLitR2Soi1,HeatSensAir2LitR,HeatSensEvapAir2LitR,HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,&
    LatentHeatAir2LitR,LWRadLitR,Radnet2LitR,TotHeatAir2LitR)

  implicit none
  integer, intent(in) :: M     !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: ATCNDR,ATCNVR,PSISV1,Radt2LitR
  real(r8), intent(in):: Prec2LitR2,PrecHeat2LitR2
  real(r8), intent(inout) :: EvapLitR2Soi1
  real(r8), intent(inout) :: HeatSensAir2LitR
  real(r8), intent(inout) :: HeatSensEvapAir2LitR
  real(r8), intent(inout) :: HeatSensLitR2Soi1
  real(r8), intent(inout) :: HeatSensVapLitR2Soi1
  real(r8), intent(inout) :: LatentHeatAir2LitR
  real(r8), intent(inout) :: LWRadLitR
  real(r8), intent(inout) :: Radnet2LitR
  real(r8), intent(inout) :: TotHeatAir2LitR
  integer  :: NN
  real(r8) :: tk1pre,VHCPRXX
  real(r8) :: RAa,VaporSoi1,VapLitR,VPY
  real(r8) :: TKXR,TK1X,TKY
  real(r8) :: FLVC,FLVX,ENGYR
  real(r8) :: HFLWC,HFLWX,EVAPR2,LatentHeatAir2LitR2
  real(r8) :: HeatSensEvapAir2LitR2
  real(r8) :: NetHeatAir2LitR2,LWRadLitR2,ThetaWLitR,HeatSensAir2LitR2,EvapLitR2Soi2,HeatSensVapLitR2Soi2
  real(r8) :: CdLitREvap,CdLitRHSens,RAGX,RI,Radnet2LitR2,TotHeatAir2LitR2,HeatSensLitR2Soi2
  real(r8) :: VLHeatCapcityLitR2,VLHeatCapacity2
  real(r8) :: VLWatMicP12,VWatLitr2
  real(r8) :: TKR1,TKS1
  real(r8) :: dt_litrHeat
! begin_execution

  VWatLitr2          = VLWatMicP1_vr(0,NY,NX)
  VLHeatCapcityLitR2 = VHeatCapacity1_vr(0,NY,NX)
  VLHeatCapacity2    = VHeatCapacity1_vr(NUM(NY,NX),NY,NX)
  VLWatMicP12        = VLWatMicP1_vr(NUM(NY,NX),NY,NX)
  TKR1               = TKSoi1_vr(0,NY,NX)
  TKS1               = TKSoi1_vr(NUM(NY,NX),NY,NX)

  !embedded iteration, local time step size
  dt_litrHeat=dts_HeatWatTP/real(NPR,kind=r8)      !time step for litter flux calculation
  D5000: DO NN=1,NPR
!    write(*,*)'TotHeatAir2LitR  x',NN,TotHeatAir2LitR
    IF(VLHeatCapcityLitR2.GT.VHeatCapLitRMin_col(NY,NX))THEN
      !
      ! AERODYNAMIC RESISTANCE ABOVE RESIDUE INCLUDING
      ! RESISTANCE IMPOSED BY PLANT CANOPY
      !
      ! RI=Richardsons number
      ! RIB=isothermal RI
      ! TKQ,TKR1=canopy air,litter temperature
      ! RZ=surface resistance to evaporation, given as a prescribed parameter
      ! RAGX,RA=litter blr
      ! RAG,RAGR=isothermal blr at ground surface
      ! CdLitREvap,CdLitRHSens=conductance for litter latent,sensible heat fluxes
      !
      RI=RichardsonNumber(RIB(NY,NX),TKQ(NY,NX),TKR1)
      RAGX=AMAX1(RAM,0.8_r8*RAGR(NY,NX),AMIN1(1.2_r8*RAGR(NY,NX),RARG(NY,NX)/(1.0_r8-10.0_r8*RI)))
      RAGR(NY,NX)=RAGX
      RAa=RAGX
      CdLitREvap=PARER(NY,NX)/(RAa+RZ)
      CdLitRHSens=PARSR(NY,NX)/RAa
!
!     NET RADIATION AT RESIDUE SURFACE
!
!     LWRadLitR2=longwave radiation emitted by litter
!     Radnet2LitR2=litter net radiation,after taking out outgoing radiation from litter 
!     ThetaWLitR=litter water content
!     VWatLitRHoldCapcity=maximum water retention by litter
!     PSISM1=litter matric water potential
!
      LWRadLitR2=LWEmscefLitR_col(NY,NX)*TKR1**4._r8/real(NPR,kind=r8)
      Radnet2LitR2=Radt2LitR-LWRadLitR2
!      write(*,*)'Radnet2LitR2',Radt2LitR,LWEmscefLitR_col(NY,NX),TKR1
      IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS2(NY,NX))THEN
        ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VWatLitr2)/VLitR_col(NY,NX)
      ELSE
        ThetaWLitR=POROS0(NY,NX)
      ENDIF

      IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX).AND.VLWatMicP1_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
        ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatMicP1_vr(0,NY,NX))/VLitR_col(NY,NX)
        IF(ThetaWLitR.LT.FieldCapacity_vr(0,NY,NX))THEN
          PSISM1_vr(0,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX) &
            +((LOGFldCapacity_vr(0,NY,NX)-LOG(ThetaWLitR))/FCD(0,NY,NX)*LOGPSIMND(NY,NX))))
        ELSEIF(ThetaWLitR.LT.POROS0(NY,NX))THEN
          PSISM1_vr(0,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(0,NY,NX)-LOG(ThetaWLitR)) &
            /PSD(0,NY,NX))**SRP(0,NY,NX)*LOGPSIMXD(NY,NX)))
        ELSE
          ThetaWLitR=POROS0(NY,NX)
          PSISM1_vr(0,NY,NX)=PSISE_vr(0,NY,NX)
        ENDIF
      ELSE
        ThetaWLitR=POROS0(NY,NX)
        PSISM1_vr(0,NY,NX)=PSISE_vr(0,NY,NX)
      ENDIF
!
!     VAPOR FLUX AT RESIDUE SURFACE
!
!     VapLitR,VaporSoi1,VPQ_col=vapor pressure in litter,soil,canopy air, m3/m-3
!     TKS1=soil temperature
!     EVAPR2=negative of litter evaporation,<0 into atmosphere
!     LatentHeatAir2LitR2=litter latent heat flux
!     VAP=latent heat of evaporation
!     HeatSensEvapAir2LitR2=convective heat of evaporation flux
!
!     in litter      
      VapLitR=vapsat(TKR1)*EXP(18.0_r8*PSISM1_vr(0,NY,NX)/(RGASC*TKR1))   
      if(abs(VapLitR)>1.e20_r8)then
        write(*,*)'TKR1=',TKR1,TKSoi1_vr(0,NY,NX),TKSoi1_vr(NUM(NY,NX),NY,NX)
        write(*,*)'PSISM1_vr(0,NY,NX)=',PSISM1_vr(0,NY,NX)
        call endrun(trim(mod_filename)//'at line',__LINE__)
      endif

      VaporSoi1 = vapsat(TKS1)*EXP(18.0_r8*PSISV1/(RGASC*TKS1))    !in soil, ton/m3
      EVAPR2    = AMAX1(-AZMAX1(VWatLitr2*dts_wat),CdLitREvap*(VPQ_col(NY,NX)-VapLitR)) ![ton hr/m]
      LatentHeatAir2LitR2   = EVAPR2*EvapLHTC          !latent energy flux, kJ/kg *10^3 kg hr/m = MJ hr/m
      HeatSensEvapAir2LitR2 = EVAPR2*cpw*TKR1        !mass energy flux
      !
      ! SOLVE FOR RESIDUE TO SOIL SURFACE HEAT FLUXES
      !
      ! FLVC,FLVX=vapor unconstrained,vapor constrained vapor flux
      ! dt_litrHeat=time step for litter flux calculations
      ! VPY=equilibrium vapor concentration
      ! VsoiPM=litter,soil air filled porosity
      ! EvapLitR2Soi2=litter soil vapor flux
      ! HeatSensVapLitR2Soi2=convective heat of litter soil vapor flux
      ! TKXR,TK1X=interim calculation of litter,soil temperatures
      ! TKY=equilibrium litter-soil temperature
      ! HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat
      ! HeatSensLitR2Soi2=litter-soil heat flux
      ! THETPM: air-filled porosity

      IF(THETPM(M,0,NY,NX).GT.THETX.AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
        FLVC=ATCNVR*(VapLitR-VaporSoi1)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree(NY,NX)*FracSurfByLitR_col(NY,NX)*dt_litrHeat
        if(abs(FLVC)>1.e20_r8)then
          write(*,*)'ATCNVR*(VapLitR-VaporSoi1)=',ATCNVR,VapLitR,VaporSoi1
          write(*,*)'FracSurfSnoFree(NY,NX),FracSurfByLitR_col(NY,NX)=',FracSurfSnoFree(NY,NX),FracSurfByLitR_col(NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif

        VPY=(VapLitR*VLsoiAirPM(M,0,NY,NX)+VaporSoi1*VLsoiAirPM(M,NUM(NY,NX),NY,NX)) &
          /(VLsoiAirPM(M,0,NY,NX)+VLsoiAirPM(M,NUM(NY,NX),NY,NX))
        FLVX=(VapLitR-VPY)*VLsoiAirPM(M,0,NY,NX)*XNPB
        if(abs(FLVX)>1.e20_r8)then
          write(*,*)'(VapLitR-VPY)*VLsoiAirPM(M,0,NY,NX)',VapLitR,VPY,VLsoiAirPM(M,0,NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif

        IF(FLVC.GE.0.0_r8)THEN
          !from litter to soil
          EvapLitR2Soi2=AZMAX1(AMIN1(FLVC,FLVX))
          if(abs(EvapLitR2Soi2)>1.e20_r8)then
            write(*,*)'FLVC,FLVX',FLVC,FLVX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          HeatSensVapLitR2Soi2=(cpw*TKR1+EvapLHTC)*EvapLitR2Soi2
        ELSE
          !from soil to litter
          EvapLitR2Soi2=AZMIN1(AMAX1(FLVC,FLVX))
          if(abs(EvapLitR2Soi2)>1.e20_r8)then
            write(*,*)'FLVC,FLVX',FLVC,FLVX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          HeatSensVapLitR2Soi2=(cpw*TKS1+EvapLHTC)*EvapLitR2Soi2
        ENDIF
      ELSE
        EvapLitR2Soi2        = 0.0_r8
        HeatSensVapLitR2Soi2 = 0.0_r8
      ENDIF

      TKXR  = TKR1-HeatSensVapLitR2Soi2/VLHeatCapcityLitR2         !update litter layer temperature
      TK1X  = TKS1+HeatSensVapLitR2Soi2/VLHeatCapacity2            !update soil layer temperature
      TKY   = (TKXR*VLHeatCapcityLitR2+TK1X*VLHeatCapacity2)/(VLHeatCapcityLitR2+VLHeatCapacity2)   !equilibrium temperature
      HFLWX = (TKXR-TKY)*VLHeatCapcityLitR2*XNPB                   !sensible heat flux > 0 into soil
      HFLWC = ATCNDR*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree(NY,NX)*FracSurfByLitR_col(NY,NX)*dt_litrHeat
      IF(HFLWC.GE.0.0_r8)THEN
        HeatSensLitR2Soi2=AZMAX1(AMIN1(HFLWX,HFLWC))
      ELSE
        HeatSensLitR2Soi2=AZMIN1(AMAX1(HFLWX,HFLWC))
      ENDIF
!
!     SOLVE FOR RESIDUE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
!     HeatSensAir2LitR2,Radnet2LitR2,LatentHeatAir2LitR2=litter sensible,net radn,latent heat fluxes
!     NetHeatAir2LitR2,TotHeatAir2LitR2=storage,total litter heat flux
!
      HeatSensAir2LitR2 = CdLitRHSens*(TKQ(NY,NX)-TKR1)    !sensible heat flux between canopy air and litter surface
      NetHeatAir2LitR2  = Radnet2LitR2+LatentHeatAir2LitR2+HeatSensAir2LitR2      !
      TotHeatAir2LitR2  = NetHeatAir2LitR2+HeatSensEvapAir2LitR2
!
!     AGGREGATE WATER AND ENERGY FLUXES FROM TIME STEP FOR LITTER
!     CALCULATIONS TO THAT FOR SOIL PROFILE
!
      VapXAir2LitR(NY,NX)  = VapXAir2LitR(NY,NX)+EVAPR2
      Radnet2LitR          = Radnet2LitR+Radnet2LitR2
      LatentHeatAir2LitR   = LatentHeatAir2LitR+LatentHeatAir2LitR2
      HeatSensEvapAir2LitR = HeatSensEvapAir2LitR+HeatSensEvapAir2LitR2
      HeatSensAir2LitR     = HeatSensAir2LitR+HeatSensAir2LitR2
      TotHeatAir2LitR      = TotHeatAir2LitR+TotHeatAir2LitR2

      ! write(*,*)'.........TotHeatAir2LitR',TotHeatAir2LitR,TotHeatAir2LitR2,TKQ(NY,NX),TKR1
      ! write(*,*)'NetHeatAir2LitR2+HeatSensEvapAir2LitR2,NN',NetHeatAir2LitR2,HeatSensEvapAir2LitR2,NN
      ! write(*,*)'Radnet2LitR2+LatentHeatAir2LitR2+HeatSensAir2LitR2',Radnet2LitR2,LatentHeatAir2LitR2,HeatSensAir2LitR2
      ! write(*,*)'CdLitRHSens*(TKQ(NY,NX)-TKR1)',CdLitRHSens,TKQ(NY,NX),TKR1
      ! call print_info(trim(mod_filename)//' at line',__LINE__)
      EvapLitR2Soi1        = EvapLitR2Soi1+EvapLitR2Soi2
      HeatSensVapLitR2Soi1 = HeatSensVapLitR2Soi1+HeatSensVapLitR2Soi2
      HeatSensLitR2Soi1    = HeatSensLitR2Soi1+HeatSensLitR2Soi2
      LWRadLitR            = LWRadLitR+LWRadLitR2
    ELSE
      !not significant litter layer heat capacity
      EVAPR2                = 0.0_r8
      Radnet2LitR2          = 0.0_r8
      LatentHeatAir2LitR2   = 0.0_r8
      HeatSensEvapAir2LitR2 = 0.0_r8
      HeatSensAir2LitR2     = 0.0_r8
      TotHeatAir2LitR2      = 0.0_r8
      EvapLitR2Soi2         = 0.0_r8
      HeatSensVapLitR2Soi2  = 0.0_r8
      HeatSensLitR2Soi2     = 0.0_r8
      LWRadLitR2            = 0.0_r8
    ENDIF
    !VWatLitr2,VLWatMicP12=water in litter, topsoil layer 
    VWatLitr2   = VWatLitr2+Prec2LitR2+EVAPR2-EvapLitR2Soi2
    VLWatMicP12 = VLWatMicP12+EvapLitR2Soi2
    ENGYR       = VLHeatCapcityLitR2*TKR1
    VHCPRXX     = VLHeatCapcityLitR2
    
    ! VLHeatCapcityLitR2: heat capacity, kJ/kg/Kelvin
    VLHeatCapcityLitR2 = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VWatLitr2+cpi*VLiceMicP1_vr(0,NY,NX)
    VLHeatCapacity2    = VLHeatCapacity2+cpw*EvapLitR2Soi2
    tk1pre             = TKR1
    TKR1               = (ENGYR+TotHeatAir2LitR2+PrecHeat2LitR2-HeatSensVapLitR2Soi2-HeatSensLitR2Soi2)/VLHeatCapcityLitR2
    TKS1               = TKS1+(HeatSensVapLitR2Soi2+HeatSensLitR2Soi2)/VLHeatCapacity2
!    print*,'tkr1',NN,TKR1,TKS1,TKQ(NY,NX),TairK_col(NY,NX),TKSoi1_vr(0,NY,NX),TKSoi1_vr(NUM(NY,NX),NY,NX)
!    print*,'cplitr',NN,VLHeatCapcityLitR2,SoilOrgM_vr(ielmc,0,NY,NX),VWatLitr2,VLiceMicP1_vr(0,NY,NX)
  ENDDO D5000
  end subroutine SurfLitterIterateNN

!------------------------------------------------------------------------------------------
  subroutine UpdateLitRPhys(I,J,NY,NX,WatInByRunoff,HEATIN_lndByRunoff,HeatStore_lnd,HEATIN_lnd)
  !
  !Description
  !Update Litter physical variables
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: NY,NX
  real(r8), intent(in) :: WatInByRunoff
  real(r8), intent(in) :: HEATIN_lndByRunoff
  real(r8), intent(inout) :: HeatStore_lnd
  real(r8), intent(inout) :: HEATIN_lnd

  real(r8) :: VHeatCapacityLitrX  !old litr heat capacity
  real(r8) :: VHeatCapacityLitR   !current litr heat capacity
  real(r8) :: dVHeatCapacityLitr  !change in heat capacity
  real(r8) :: tkspre,ENGYR,VLWatMicPr,VLiceMicPr
  real(r8) :: ENGYZ,HeatByLitrMassChange,HS
  integer :: LS

  ! CALCULATE SURFACE RESIDUE TEMPERATURE FROM ITS CHANGE
  ! IN HEAT STORAGE
  !
  VHeatCapacityLitrX = VHeatCapacity_vr(0,NY,NX)
  VHeatCapacityLitR  = cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)+cpo*SoilOrgM_vr(ielmc,0,NY,NX)
  VLWatMicPr         = VLWatMicP_vr(0,NY,NX)
  VLiceMicPr         = VLiceMicP_vr(0,NY,NX)
  if(VLWatMicP_vr(0,NY,NX)<0._r8 .or. VLiceMicP_vr(0,NY,NX)<0._r8 .or. VHeatCapacityLitR<0._r8)then
    write(*,*)'negative litr water',VLWatMicP_vr(0,NY,NX),VLiceMicP_vr(0,NY,NX),VHeatCapacityLitR,SoilOrgM_vr(ielmc,0,NY,NX) 
    call endrun(trim(mod_filename)//' at line',__LINE__)    
  endif
  !TairK: air temperature in kelvin, HeatByLitrMassChange represents increase heat in litr
  dVHeatCapacityLitr   = VHeatCapacityLitR-VHeatCapacityLitrX
  HeatByLitrMassChange = dVHeatCapacityLitr*TairK_col(NY,NX)
  ENGYZ                = VHeatCapacityLitrX*TKS_vr(0,NY,NX)

  !update water, ice content and heat capacity of residue
  VLWatMicP_vr(0,NY,NX)     = AZMAX1(VLWatMicP_vr(0,NY,NX)+WatFLo2Litr(NY,NX)+TLitrIceFlxThaw(NY,NX)+WatInByRunoff)
  VLiceMicP_vr(0,NY,NX)     = AZMAX1(VLiceMicP_vr(0,NY,NX)-TLitrIceFlxThaw(NY,NX)/DENSICE)
  VHeatCapacity_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)

  IF(VHeatCapacity_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    !when there are still significant heat capacity of the residual layer
    tkspre=TKS_vr(0,NY,NX)
    TKS_vr(0,NY,NX)=(ENGYZ+HeatFLo2LitrByWat(NY,NX)+TLitrIceHeatFlxFrez(NY,NX)+HeatByLitrMassChange &
      +HEATIN_lndByRunoff)/VHeatCapacity_vr(0,NY,NX)
    HS                   = TKS_vr(0,NY,NX)*VHeatCapacity_vr(0,NY,NX)
    HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+HS
    HeatStore_lnd        = HeatStore_lnd+HS
    if(TKS_vr(0,NY,NX)<100._r8 .or. TKS_vr(0,NY,NX)>360._r8)then
      write(*,*)I,J,NY,NX,TKS_vr(0,NY,NX),tkspre
      write(*,*)'WatFLo2Litr, WatInByRunoff=',WatFLo2Litr(NY,NX),WatInByRunoff
      write(*,*)'wat flo2litr icethaw runoff',VLWatMicPr,VLWatMicP_vr(0,NY,NX),WatFLo2Litr(NY,NX),TLitrIceFlxThaw(NY,NX),WatInByRunoff
      write(*,*)'ice',VLiceMicPr,VLiceMicP_vr(0,NY,NX)
      write(*,*)'engy',ENGYZ,HeatFLo2LitrByWat(NY,NX),TLitrIceHeatFlxFrez(NY,NX),HeatByLitrMassChange, &
        HEATIN_lndByRunoff,VHeatCapacity_vr(0,NY,NX)        
      write(*,*)'vhc',VHeatCapacityLitrX,VHeatCapacityLitR,dVHeatCapacityLitR,TairK_col(NY,NX),SoilOrgM_vr(ielmc,0,NY,NX)   
      write(*,*)'tengz',ENGYZ/VHeatCapacity_vr(0,NY,NX),HeatFLo2LitrByWat(NY,NX)/VHeatCapacity_vr(0,NY,NX),&
        TLitrIceHeatFlxFrez(NY,NX)/VHeatCapacity_vr(0,NY,NX),HeatByLitrMassChange/VHeatCapacity_vr(0,NY,NX), &
        HEATIN_lndByRunoff/VHeatCapacity_vr(0,NY,NX)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif  
    HEATIN_lnd=HEATIN_lnd+HeatByLitrMassChange
    Ls=NUM(NY,NX)
    !if(curday>=175)write(*,*)'at line',__LINE__,TKS_vr(0,NY,NX),TKS_vr(Ls,ny,nx),tkspre
!    if(abs(VHeatCapacity_vr(0,NY,NX)/VHeatCapacityLitrX-1._r8)>0.025_r8.or. &
!      abs(TKS_vr(0,NY,NX)/tkspre-1._r8)>0.025_r8)then
!      TKS_vr(0,NY,NX)=TKS_vr(NUM(NY,NX),NY,NX)
!    endif
  ELSE
    HEATIN_lnd=HEATIN_lnd+HeatByLitrMassChange+(TKS_vr(NUM(NY,NX),NY,NX)-TKS_vr(0,NY,NX))*VHeatCapacity_vr(0,NY,NX)
    TKS_vr(0,NY,NX)=TKS_vr(NUM(NY,NX),NY,NX)
  ENDIF

  TCS(0,NY,NX)=units%Kelvin2Celcius(TKS_vr(0,NY,NX))
    
  ENGYR=VHeatCapacity_vr(0,NY,NX)*TKS_vr(0,NY,NX)

  HeatStore_lnd=HeatStore_lnd+ENGYR
  HEATIN_lnd=HEATIN_lnd+TLitrIceHeatFlxFrez(NY,NX)

  end subroutine UpdateLitRPhys

!------------------------------------------------------------------------------------------
  subroutine UpdateLitRB4RunoffM(I,J,M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX,I,J
  real(r8) :: VOLIRZ,ENGYR,VLHeatCapLitRPre
  real(r8) :: TK0Prev,TVWatIceLitR,VWatLitrZ
  real(r8) :: VLWatMicP10,VLiceMicP10
  ! SURFACE RESIDUE WATER AND TEMPERATURE
  !
  ! XVOLT,XVOLW=free water+ice,water in litter layer
  ! VOLWM,VsoiPM=surface water,air content for use in TranspNoSalt.f
  ! VWatLitRHoldCapcity=maximum water retention by litter
  ! VHeatCapacity1_vr=volumetric heat capacity of litter
  ! VOLA1,VLWatMicP1,VLiceMicP1,VOLP1=pore,water,ice,air volumes of litter
  ! VWatLitRHoldCapcity=maximum water retention by litter
  ! LitrIceHeatFlxFrez,LitrIceFlxThaw=litter water,latent heat flux from freeze-thaw
  ! VLitR=dry litter volume
  ! THETWX,FracSoiPAsIce,FracSoiPAsAir=water,ice,air concentrations
  ! VHeatCapacity1_vr=volumetric heat capacity of litter
  ! TK1=litter temperature
  ! HFLWRL,LitrIceHeatFlxFrez,cumHeatFlx2LitRByRunoff=litter total cond+conv,latent,runoff heat flux
  VLWatMicP10                   = VLWatMicP1_vr(0,NY,NX)
  VLiceMicP10                   = VLiceMicP1_vr(0,NY,NX)
  VLWatMicP1_vr(0,NY,NX)        = AZMAX1(VLWatMicP1_vr(0,NY,NX)+WatFLow2LitR_col(NY,NX)+LitrIceFlxThaw(NY,NX))
  VLiceMicP1_vr(0,NY,NX)        = AZMAX1(VLiceMicP1_vr(0,NY,NX)-LitrIceFlxThaw(NY,NX)/DENSICE)
  VLairMicP1_vr(0,NY,NX)        = AZMAX1(VLPoreLitR(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))
  VLWatMicPM_vr(M+1,0,NY,NX)    = VLWatMicP1_vr(0,NY,NX)
  VLsoiAirPM(M+1,0,NY,NX)       = VLairMicP1_vr(0,NY,NX)
  TVWatIceLitR                  = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX) = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ               = VLWatMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    VOLIRZ                  = VLiceMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    XVLMobileWatMicP(NY,NX) = AZMAX1(VLWatMicP1_vr(0,NY,NX)-VWatLitrZ)
    XVLiceMicP_col(NY,NX)   = AZMAX1(VLiceMicP1_vr(0,NY,NX)-VOLIRZ)
  ELSE
    XVLMobileWatMicP(NY,NX) = 0.0_r8
    XVLiceMicP_col(NY,NX)   = 0.0_r8
  ENDIF
  XVLMobileWaterLitRM(M+1,NY,NX) = XVLMobileWaterLitR_col(NY,NX)
  XVLMobileWatMicPM(M+1,NY,NX)   = XVLMobileWatMicP(NY,NX)
  XVLiceMicPM(M+1,NY,NX)         = XVLiceMicP_col(NY,NX)
  IF(VLitR_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat_vr(0,NY,NX)=AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsIce_vr(0,NY,NX)=AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsAir_vr(0,NY,NX)=AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)) &
      *AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/MaxVLWatByLitR_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)=0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)=0.0_r8
    FracSoiPAsAir_vr(0,NY,NX)=1.0_r8
  ENDIF
  THETPM(M+1,0,NY,NX) = FracSoiPAsAir_vr(0,NY,NX)
  VLHeatCapLitRPre    = VHeatCapacity1_vr(0,NY,NX)                !heat capacity
  TK0Prev             = TKSoi1_vr(0,NY,NX)                                 !residual temperature
  ENGYR               = VHeatCapacity1_vr(0,NY,NX)*TKSoi1_vr(0,NY,NX)  !initial energy content
  VHeatCapacity1_vr(0,NY,NX)=cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1_vr(0,NY,NX)+cpi*VLiceMicP1_vr(0,NY,NX)  !update heat capacity

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    TKSoi1_vr(0,NY,NX)=(ENGYR+HeatFLoByWat2LitRi_col(NY,NX)+LitrIceHeatFlxFrez(NY,NX))/VHeatCapacity1_vr(0,NY,NX)
!    if(I==280)print*,I+J/24.,J,'cpksoi0',VHeatCapacity1_vr(0,NY,NX),SoilOrgM_vr(ielmc,0,NY,NX),VLWatMicP1_vr(0,NY,NX),VLiceMicP1_vr(0,NY,NX)
    if(TKSoi1_vr(0,NY,NX)<100._r8 .or. TKSoi1_vr(0,NY,NX)>400._r8)then
      write(*,*)'weird litter temp UpdateLitRB4RunoffM=',I,J,M,NY,NX,TKSoi1_vr(0,NY,NX),TK0Prev,VHeatCapacity1_vr(0,NY,NX)
      write(*,*)'engy',ENGYR/VHeatCapacity1_vr(0,NY,NX),HeatFLoByWat2LitRi_col(NY,NX)/VHeatCapacity1_vr(0,NY,NX),&
        LitrIceHeatFlxFrez(NY,NX)/VHeatCapacity1_vr(0,NY,NX)
      write(*,*)'cpo',cpo*SoilOrgM_vr(ielmc,0,NY,NX),cpw*VLWatMicP1_vr(0,NY,NX),cpi*VLiceMicP1_vr(0,NY,NX),VLHeatCapLitRPre,VHeatCapLitRMin_col(NY,NX) 
      write(*,*)'cpw',cpw*VLWatMicP10,VLWatMicP10,VLiceMicP10,VLWatMicP1_vr(0,NY,NX),VLiceMicP1_vr(0,NY,NX)
      write(*,*)'watflw',ENGYR,WatFLow2LitR_col(NY,NX),HeatFLoByWat2LitRi_col(NY,NX),HeatFLoByWat2LitRi_col(NY,NX)/(WatFLow2LitR_col(NY,NX)*cpw)
      write(*,*)'vlwat',VLWatMicP10,WatFLow2LitR_col(NY,NX),cumWatFlx2LitRByRunoff(NY,NX),cumHeatFlx2LitRByRunoff(NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
!    IF(ABS(VHeatCapacity1_vr(0,NY,NX)/VLHeatCapLitRPre-1._r8)>0.025_r8.or. &
!      abs(TKSoi1_vr(0,NY,NX)/TK0Prev-1._r8)>0.025_r8)THEN
!      TKSoi1_vr(0,NY,NX)=TKSoi1_vr(NUM(NY,NX),NY,NX)
!    ENDIF
  ELSE
    TKSoi1_vr(0,NY,NX)=TKSoi1_vr(NUM(NY,NX),NY,NX)
  ENDIF
!  if(NY==1 .AND. NX==1)write(191,*)'M  UpdateLitRB4RunoffM',M,VLWatMicP_vr(0,NY,NX),VLWatMicP10,VLWatMicP1_vr(0,NY,NX),&
!    watflw(NY,NX),waticefl(NY,NX)

  watflw(NY, NX)=watflw(NY,NX)+WatFLow2LitR_col(NY,NX)
  waticefl(NY,NX)=waticefl(NY,NX)+LitrIceFlxThaw(NY,NX)
!  if(NY==1 .AND. NX==1)write(*,*)'afM  UpdateLitRB4RunoffM',M,watflw(NY,NX),waticefl(NY,NX),VLWatMicP10,WatFLow2LitR_col(NY,NX)

  end subroutine UpdateLitRB4RunoffM

!------------------------------------------------------------------------------------------
  subroutine UpdateLitRAftRunoff(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VOLIRZ,ENGYR,VLHeatCapLitRPre
  real(r8) :: TK0Prev,TVWatIceLitR,VWatLitrZ
  real(r8) :: VLWatMicP10,VLiceMicP10
  ! SURFACE RESIDUE WATER AND TEMPERATURE
  !
  ! XVOLT,XVOLW=free water+ice,water in litter layer
  ! VOLWM,VsoiPM=surface water,air content for use in TranspNoSalt.f
  ! VWatLitRHoldCapcity=maximum water retention by litter
  ! VHeatCapacity1_vr=volumetric heat capacity of litter
  ! VOLA1,VLWatMicP1,VLiceMicP1,VOLP1=pore,water,ice,air volumes of litter
  ! VWatLitRHoldCapcity=maximum water retention by litter
  ! LitrIceHeatFlxFrez,LitrIceFlxThaw=litter water,latent heat flux from freeze-thaw
  ! VLitR=dry litter volume
  ! THETWX,FracSoiPAsIce,FracSoiPAsAir=water,ice,air concentrations
  ! VHeatCapacity1_vr=volumetric heat capacity of litter
  ! TK1=litter temperature
  ! HFLWRL,LitrIceHeatFlxFrez,cumHeatFlx2LitRByRunoff=litter total cond+conv,latent,runoff heat flux
  VLWatMicP10 = VLWatMicP1_vr(0,NY,NX)
  VLiceMicP10 = VLiceMicP1_vr(0,NY,NX)

  VLWatMicP1_vr(0,NY,NX)     = AZMAX1(VLWatMicP1_vr(0,NY,NX)+cumWatFlx2LitRByRunoff(NY,NX))
  VLairMicP1_vr(0,NY,NX)     = AZMAX1(VLPoreLitR(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))
  VLWatMicPM_vr(M+1,0,NY,NX) = VLWatMicP1_vr(0,NY,NX)
  VLsoiAirPM(M+1,0,NY,NX)    = VLairMicP1_vr(0,NY,NX)
  TVWatIceLitR               = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX)  = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ               = VLWatMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    VOLIRZ                  = VLiceMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    XVLMobileWatMicP(NY,NX) = AZMAX1(VLWatMicP1_vr(0,NY,NX)-VWatLitrZ)
    XVLiceMicP_col(NY,NX)   = AZMAX1(VLiceMicP1_vr(0,NY,NX)-VOLIRZ)
  ELSE
    XVLMobileWatMicP(NY,NX) = 0.0_r8
    XVLiceMicP_col(NY,NX)   = 0.0_r8
  ENDIF
  XVLMobileWaterLitRM(M+1,NY,NX) = XVLMobileWaterLitR_col(NY,NX)
  XVLMobileWatMicPM(M+1,NY,NX)   = XVLMobileWatMicP(NY,NX)
  XVLiceMicPM(M+1,NY,NX)         = XVLiceMicP_col(NY,NX)
  
  IF(VLitR_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat_vr(0,NY,NX)=AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsIce_vr(0,NY,NX)=AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsAir_vr(0,NY,NX)=AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)) &
      *AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/MaxVLWatByLitR_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)=0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)=0.0_r8
    FracSoiPAsAir_vr(0,NY,NX)=1.0_r8
  ENDIF
  THETPM(M+1,0,NY,NX)        = FracSoiPAsAir_vr(0,NY,NX)
  VLHeatCapLitRPre           = VHeatCapacity1_vr(0,NY,NX)                !heat capacity
  TK0Prev                    = TKSoi1_vr(0,NY,NX)                                 !residual temperature
  ENGYR                      = VHeatCapacity1_vr(0,NY,NX)*TKSoi1_vr(0,NY,NX)  !initial energy content
  VHeatCapacity1_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1_vr(0,NY,NX)+cpi*VLiceMicP1_vr(0,NY,NX)  !update heat capacity

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    TKSoi1_vr(0,NY,NX)=(ENGYR+cumHeatFlx2LitRByRunoff(NY,NX))/VHeatCapacity1_vr(0,NY,NX)

    if(TKSoi1_vr(0,NY,NX)<100._r8 .or. TKSoi1_vr(0,NY,NX)>400._r8)then
!      if(TKSoi1_vr(0,NY,NX)-TK0Prev>20._r8 .and. HeatFLoByWat2LitRi_col(NY,NX)<0._r8)then
!        dEvap=AMIN1((TKSoi1_vr(0,NY,NX)-TK0Prev)*VHeatCapacity1_vr(0,NY,NX)/EvapLHTC,VLWatMicP1_vr(0,NY,NX))
!        VLWatMicP1_vr(0,NY,NX)=VLWatMicP1_vr(0,NY,NX)-dEvap

!        VHeatCapacity1_vr(0,NY,NX)=cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1_vr(0,NY,NX)+cpi*VLiceMicP1_vr(0,NY,NX)
!      endif
      write(*,*)'weird litter temp UpdateLitRAftRunoff  =',M,NY,NX,TKSoi1_vr(0,NY,NX),TK0Prev,VHeatCapacity1_vr(0,NY,NX)
      write(*,*)ENGYR/VHeatCapacity1_vr(0,NY,NX),cumHeatFlx2LitRByRunoff(NY,NX)/VHeatCapacity1_vr(0,NY,NX)
      write(*,*)cpo*SoilOrgM_vr(ielmc,0,NY,NX),cpw*VLWatMicP1_vr(0,NY,NX),cpi*VLiceMicP1_vr(0,NY,NX),VLHeatCapLitRPre,VHeatCapLitRMin_col(NY,NX) 
      write(*,*)cpw*VLWatMicP10,VLWatMicP10,VLiceMicP10,VLWatMicP1_vr(0,NY,NX),VLiceMicP1_vr(0,NY,NX)
      write(*,*)VLWatMicP10,cumWatFlx2LitRByRunoff(NY,NX),cumHeatFlx2LitRByRunoff(NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
!    IF(ABS(VHeatCapacity1_vr(0,NY,NX)/VLHeatCapLitRPre-1._r8)>0.025_r8.or. &
!      abs(TKSoi1_vr(0,NY,NX)/TK0Prev-1._r8)>0.025_r8)THEN
!      TKSoi1_vr(0,NY,NX)=TKSoi1_vr(NUM(NY,NX),NY,NX)
!    ENDIF
  ELSE
    TKSoi1_vr(0,NY,NX)=TKSoi1_vr(NUM(NY,NX),NY,NX)
  ENDIF

  end subroutine UpdateLitRAftRunoff  
end module SurfLitterPhysMod