module SnowBalMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : endrun
  USE SnowDataType
  use GridConsts
  use GridDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use SoilHeatDataType
  use ClimForcDataType
  use TFlxTypeMod
  use EcosimConst
  USE SoilPhysDataType
  USE EcoSimSumDataType
  USE FlagDataType
implicit none

  private

  character(len=*), parameter :: mod_filename = __FILE__

  public :: SnowDynUpdate
  contains

  subroutine SnowDynUpdate(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: VOLSWI,ENGYW
  integer :: L
      !
      !     CALCULATE SNOWPACK TEMPERATURE FROM ITS CHANGE
      !     IN HEAT STORAGE
      !
  VOLSS(NY,NX)=0.0_r8
  VOLWS(NY,NX)=0.0_r8
  VOLIS(NY,NX)=0.0_r8
  VOLS(NY,NX)=0.0_r8
  DPTHS(NY,NX)=0.0_r8
  VOLSWI=0.0_r8
  DO 9780 L=1,JS
    call UpdateSnowLayers(L,NY,NX,VOLSWI)

    call UpdateSoluteInSnow(L,NY,NX)
9780  CONTINUE
!
!     SNOW RUNOFF
!
  VOLSSL(1,NY,NX)=VOLSSL(1,NY,NX)+TQS(NY,NX)
  VOLWSL(1,NY,NX)=VOLWSL(1,NY,NX)+TQW(NY,NX)
  VOLISL(1,NY,NX)=VOLISL(1,NY,NX)+TQI(NY,NX)
  ENGYW=VHCPW(1,NY,NX)*TKW(1,NY,NX)
  VHCPW(1,NY,NX)=cps*VOLSSL(1,NY,NX)+cpw*VOLWSL(1,NY,NX) &
    +cpi*VOLISL(1,NY,NX)
  IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
    TKW(1,NY,NX)=(ENGYW+THQS(NY,NX))/VHCPW(1,NY,NX)
  ELSE
    TKW(1,NY,NX)=TKA(NY,NX)
  ENDIF
!
! IF SNOWPACK DISAPPEARS
  call SnowpackDisapper(NY,NX)

  end subroutine SnowDynUpdate

!------------------------------------------------------------------------------------------

  subroutine SnowpackDisapper(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
  real(r8) :: tksx
  real(r8) :: FLWS,FLWW,FLWI
  real(r8) :: HFLWS,ENGYS,ENGY2,ENGY1,ENGY
!     begin_execution
!
  IF(VHCPW(1,NY,NX).GT.0.0_r8.AND.VHCPW(1,NY,NX).LE.VHCPWX(NY,NX) &
    .AND.TKA(NY,NX).GT.TFICE)THEN
    ENGYS=TKW(1,NY,NX)*VHCPW(1,NY,NX)
    ENGY1=TKS(NUM(NY,NX),NY,NX)*VHCP(NUM(NY,NX),NY,NX)
    FLWS=VOLSSL(1,NY,NX)
    FLWW=VOLWSL(1,NY,NX)
    FLWI=VOLISL(1,NY,NX)
    HFLWS=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TKW(1,NY,NX)
    VOLSSL(1,NY,NX)=0.0_r8
    VOLWSL(1,NY,NX)=0.0_r8
    VOLISL(1,NY,NX)=0.0_r8
    VHCPW(1,NY,NX)=0.0_r8
    VOLSS(NY,NX)=0.0_r8
    VOLWS(NY,NX)=0.0_r8
    VOLIS(NY,NX)=0.0_r8
    VOLS(NY,NX)=0.0_r8
    DPTHS(NY,NX)=0.0_r8
    DO 9770 L=1,JS
      DENSS(L,NY,NX)=DENS0(NY,NX)
9770  CONTINUE
    VOLW(0,NY,NX)=VOLW(0,NY,NX)+FLWW
    VOLI(0,NY,NX)=VOLI(0,NY,NX)+FLWI+FLWS/DENSI
    ENGY=VHCP(NUM(NY,NX),NY,NX)*TKS(NUM(NY,NX),NY,NX)
    VHCP(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX) &
      +cpw*(VOLW(NUM(NY,NX),NY,NX)+VOLWH(NUM(NY,NX),NY,NX)) &
      +cpi*(VOLI(NUM(NY,NX),NY,NX)+VOLIH(NUM(NY,NX),NY,NX))
    TKSX=TKS(NUM(NY,NX),NY,NX)
    IF(VHCP(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      TKS(NUM(NY,NX),NY,NX)=(ENGY+HFLWS)/VHCP(NUM(NY,NX),NY,NX)
    ELSE
      TKS(NUM(NY,NX),NY,NX)=TKA(NY,NX)
    ENDIF
    ENGY2=VHCP(NUM(NY,NX),NY,NX)*TKS(NUM(NY,NX),NY,NX)
!     WRITE(*,2222)'SNW',I,J,IYRC,NX,NY,HEATIN
!    2,(cps-cpi)*FLWS/DENSI*TKSX,HFLWS,ENGY
!    3,VHCP(NUM(NY,NX),NY,NX)*TKS(NUM(NY,NX),NY,NX),TKSX
!    3,TKS(NUM(NY,NX),NY,NX),FLWS,DENSI,TKW(1,NY,NX),VHCPWZ(1,NY,NX)
!    4,ENGY2,ENGYS,ENGY1
!    5,VOLW(NUM(NY,NX),NY,NX),VOLI(NUM(NY,NX),NY,NX)
  ENDIF
  end subroutine SnowpackDisapper

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowLayers(L,NY,NX,VOLSWI)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(inout) :: VOLSWI
  real(r8) :: TKWX,VHCPWZ(JZ,JY,JX)
  real(r8) :: ENGYW
  real(r8) :: VOLSF,TCASF
  real(r8) :: DENSF,DDENS2,DDENS1
  real(r8) :: CVISC,DENSX
! begin_execution
!
! ADD CHANGES IN SNOW, WATER AND ICE
!
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
! XWFLXS,XWFLXI=freeze-thaw flux from watsub.f
!
!
  VOLSSL(L,NY,NX)=VOLSSL(L,NY,NX)+TFLWS(L,NY,NX)-XWFLXS(L,NY,NX)
  VOLWSL(L,NY,NX)=VOLWSL(L,NY,NX)+TFLWW(L,NY,NX)+XWFLXS(L,NY,NX) &
      +XWFLXI(L,NY,NX)
  if(VOLWSL(L,NY,NX)<0._r8)then
    if(L>1)then
      VOLWSL(L-1,NY,NX)=VOLWSL(L-1,NY,NX)+VOLWSL(L,NY,NX)
    elseif(L==1)then
      VOLW(0,NY,NX)=VOLW(0,NY,NX)+VOLWSL(L,NY,NX)
    endif
    VOLWSL(L,NY,NX)=0._r8
  endif
  VOLISL(L,NY,NX)=VOLISL(L,NY,NX)+TFLWI(L,NY,NX)-XWFLXI(L,NY,NX)/DENSI
  if(VOLWSL(L,NY,NX)<0._r8)then
    if(L>1)then
      VOLISL(L-1,NY,NX)=VOLISL(L-1,NY,NX)+VOLISL(L,NY,NX)
    elseif(L==1)then
      VOLI(0,NY,NX)=VOLI(0,NY,NX)+VOLISL(L,NY,NX)
    endif
    VOLISL(L,NY,NX)=0._r8
  endif
!  if(VOLISL(L,NY,NX)/=VOLISL(L,NY,NX))then
!    call print_info('VOLISL(L,NY,NX)=VOLISL(L,NY,NX)', &
!      (/padr('TFLWI',10),padr('XWFLXI',10)/),(/TFLWI(L,NY,NX) &
!      ,XWFLXI(L,NY,NX)/))
!  endif
!
! ACCUMULATE SNOW MASS FOR CALCULATING COMPRESSION
!
! VOLWSI=accumulated water equivalent volume in snowpack
! XFLWS=snow transfer from watsub.f
! VOLSF=snowfall volume
! DENSS=snow density in layer
!
  IF(L.EQ.1)THEN
    VOLSWI=VOLSWI+0.5*(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
      +VOLISL(L,NY,NX)*DENSI)
    if(VOLSWI<0._r8)then
      write(*,*)'VOLSWI=',VOLSWI,VOLSSL(L,NY,NX),VOLWSL(L,NY,NX), &
        VOLISL(L,NY,NX)*DENSI
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
!
!   RESET SNOW SURFACE DENSITY FOR SNOWFALL
!
    IF(XFLWS(L,NY,NX).GT.0.0)THEN
      DENSX=DENSS(L,NY,NX)
      TCASF=AMAX1(-15.0,AMIN1(2.0,TCA(NY,NX)))
      DENSF=0.05+1.7E-03*(TCASF+15.0)**1.5
      VOLSF=XFLWS(L,NY,NX)/DENSF &
        +(VOLSSL(L,NY,NX)-XFLWS(L,NY,NX))/DENSS(L,NY,NX)
      DENSS(L,NY,NX)=VOLSSL(L,NY,NX)/VOLSF
      if(DENSS(L,NY,NX)<0._r8)write(*,*)'VOLSSL=',VOLSSL(L,NY,NX),VOLSF
    ENDIF
  ELSE
    VOLSWI=VOLSWI+0.5*(VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX) &
      +VOLISL(L-1,NY,NX)*DENSI+VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
      +VOLISL(L,NY,NX)*DENSI)
    if(VOLSWI<0._r8)then
      write(*,*)'iVOLSWI=',VOLSWI,VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX) &
        +VOLISL(L-1,NY,NX)*DENSI,VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
        +VOLISL(L,NY,NX)*DENSI
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDIF
!
! SNOWPACK COMPRESSION
!
! DDENS1,DDENS2=temperature, compression effect on snow density
! DENSS=snow density in layer
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! VOLSL=snowpack layer volume
! DLYRS=snowpack layer depth
! CDPTHS=cumulative depth to bottom of snowpack layer
! VHCPW=snowpack layer heat capacity
! TKW,TCW=snowpack layer temperature K,oC
! THFLWW=convective heat fluxes of snow,water,ice in snowpack
! XTHAWW=latent heat flux from freeze-thaw from watsub.f
! HEATIN=cumulative net surface heat transfer
! VOLSS,VOLWS,VOLIS=total snow water equivalent, water, ice content of snowpack
! VOLS,DPTHS=total snowpack volume, depth
!
  IF(DENSS(L,NY,NX).LT.0.25)THEN
    DDENS1=DENSS(L,NY,NX)*1.0E-05*EXP(0.04*TCW(L,NY,NX))
  ELSE
    DDENS1=0.0_r8
  ENDIF
  CVISC=0.25_r8*EXP(-0.08_r8*TCW(L,NY,NX)+23.0_r8*DENSS(L,NY,NX))
  !if(curday>=83)write(*,*)'CVISC=',CVISC,TCW(L,NY,NX),DENSS(L,NY,NX)
  DDENS2=DENSS(L,NY,NX)*VOLSWI/(AREA(3,NU(NY,NX),NY,NX)*CVISC)
  if(DDENS2<0._r8)write(*,*)'DDENS2=',DENSS(L,NY,NX),VOLSWI,CVISC
  DENSS(L,NY,NX)=DENSS(L,NY,NX)+DDENS1+DDENS2
  if(DENSS(L,NY,NX)<0._r8)write(*,*)'DDENS1=',DENSS(L,NY,NX),DDENS1,DDENS2
  IF(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX) &
    .GT.ZEROS2(NY,NX))THEN
    VOLSL(L,NY,NX)=VOLSSL(L,NY,NX)/DENSS(L,NY,NX) &
      +VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)
    DLYRS(L,NY,NX)=AMAX1(0.0,VOLSL(L,NY,NX)) &
      /AREA(3,NU(NY,NX),NY,NX)
    CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)+DLYRS(L,NY,NX)
    VHCPWZ(L,NY,NX)=VHCPW(L,NY,NX)
    TKWX=TKW(L,NY,NX)
    ENGYW=VHCPW(L,NY,NX)*TKW(L,NY,NX)
    VHCPW(L,NY,NX)=cps*VOLSSL(L,NY,NX)+cpw*VOLWSL(L,NY,NX) &
      +cpi*VOLISL(L,NY,NX)
    IF(VHCPWZ(L,NY,NX).GT.VHCPWX(NY,NX) &
      .AND.VHCPW(L,NY,NX).GT.ZEROS(NY,NX))THEN
      TKW(L,NY,NX)=(ENGYW+THFLWW(L,NY,NX)+XTHAWW(L,NY,NX))/VHCPW(L,NY,NX)
    ELSE
      IF(L.EQ.1)THEN
        TKW(L,NY,NX)=TKA(NY,NX)
      ELSE
        TKW(L,NY,NX)=TKW(L-1,NY,NX)
      ENDIF
      IF(VHCPWZ(L,NY,NX).GT.ZEROS(NY,NX))THEN
        HEATIN=HEATIN+(TKW(L,NY,NX)-TKWX)*VHCPW(L,NY,NX)
      ENDIF
    ENDIF
    VOLSS(NY,NX)=VOLSS(NY,NX)+VOLSSL(L,NY,NX)
    VOLWS(NY,NX)=VOLWS(NY,NX)+VOLWSL(L,NY,NX)
    VOLIS(NY,NX)=VOLIS(NY,NX)+VOLISL(L,NY,NX)
    VOLS(NY,NX)=VOLS(NY,NX)+VOLSL(L,NY,NX)
    DPTHS(NY,NX)=DPTHS(NY,NX)+DLYRS(L,NY,NX)
!     IF(L.EQ.5)THEN
!     WRITE(*,7753)'VOLSS',I,J,NX,NY,L
!    3,VOLSSL(L,NY,NX),VOLWSL(L,NY,NX),VOLISL(L,NY,NX),VOLSL(L,NY,NX)
!    3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
!    3,XFLWS(L+1,NY,NX),XFLWW(L+1,NY,NX),XFLWI(L+1,NY,NX)
!    4,TFLWS(L,NY,NX),TFLWW(L,NY,NX),TFLWI(L,NY,NX)
!    5,XWFLXS(L,NY,NX),XWFLXI(L,NY,NX)
!    3,VOLSS(NY,NX),VOLWS(NY,NX),VOLIS(NY,NX),VOLS(NY,NX),DPTHS(NY,NX)
!    5,DLYRS(L,NY,NX),DENSS(L,NY,NX),DDENS1,DDENS2,TCW(L,NY,NX)
!    6,CVISC,VOLSWI,TCA(NY,NX)
!     WRITE(*,7753)'TKW',I,J,NX,NY,L,TKW(L,NY,NX),TCW(L,NY,NX)
!    3,THFLWW(L,NY,NX),XTHAWW(L,NY,NX),XHFLWW(L,NY,NX)
!    4,XHFLWW(L+1,NY,NX),HFLSWR(L,NY,NX),HFLSW(L,NY,NX)
!    4,VHCPWX(NY,NX),VHCPW(L,NY,NX),VHCPWZ(L,NY,NX)
!7753  FORMAT(A8,5I4,100E14.6)
!     ENDIF
  ELSE
    VOLSSL(L,NY,NX)=0.0_r8
    VOLWSL(L,NY,NX)=0.0_r8
    VOLISL(L,NY,NX)=0.0_r8
    VOLSL(L,NY,NX)=0.0_r8
    DLYRS(L,NY,NX)=0.0_r8
    CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)
    VHCPW(L,NY,NX)=0.0_r8
    IF(L.EQ.1)THEN
      TKW(L,NY,NX)=TKA(NY,NX)
    ELSE
      TKW(L,NY,NX)=TKW(L-1,NY,NX)
    ENDIF
  ENDIF
  TCW(L,NY,NX)=TKW(L,NY,NX)-TC2K
  end subroutine UpdateSnowLayers
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSnow(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
  !     begin_execution
  !     SNOWPACK SOLUTE CONTENT
  !
  !     *W2=solute content of snowpack
  !     T*BLS=net solute flux in snowpack
  !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
  !             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
  !
  CO2W(L,NY,NX)=CO2W(L,NY,NX)+TCOBLS(L,NY,NX)
  CH4W(L,NY,NX)=CH4W(L,NY,NX)+TCHBLS(L,NY,NX)
  OXYW(L,NY,NX)=OXYW(L,NY,NX)+TOXBLS(L,NY,NX)
  ZNGW(L,NY,NX)=ZNGW(L,NY,NX)+TNGBLS(L,NY,NX)
  ZN2W(L,NY,NX)=ZN2W(L,NY,NX)+TN2BLS(L,NY,NX)
  ZN4W(L,NY,NX)=ZN4W(L,NY,NX)+TN4BLW(L,NY,NX)
  ZN3W(L,NY,NX)=ZN3W(L,NY,NX)+TN3BLW(L,NY,NX)
  ZNOW(L,NY,NX)=ZNOW(L,NY,NX)+TNOBLW(L,NY,NX)
  Z1PW(L,NY,NX)=Z1PW(L,NY,NX)+TH1PBS(L,NY,NX)
  ZHPW(L,NY,NX)=ZHPW(L,NY,NX)+TH2PBS(L,NY,NX)
  !
  !     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
  !          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
  !          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
  !          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
  !          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
  !          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
  !          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
  !     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
  !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
  !
  IF(ISALTG.NE.0)THEN
    ZALW(L,NY,NX)=ZALW(L,NY,NX)+TALBLS(L,NY,NX)
    ZFEW(L,NY,NX)=ZFEW(L,NY,NX)+TFEBLS(L,NY,NX)
    ZHYW(L,NY,NX)=ZHYW(L,NY,NX)+THYBLS(L,NY,NX)
    ZCAW(L,NY,NX)=ZCAW(L,NY,NX)+TCABLS(L,NY,NX)
    ZMGW(L,NY,NX)=ZMGW(L,NY,NX)+TMGBLS(L,NY,NX)
    ZNAW(L,NY,NX)=ZNAW(L,NY,NX)+TNABLS(L,NY,NX)
    ZKAW(L,NY,NX)=ZKAW(L,NY,NX)+TKABLS(L,NY,NX)
    ZOHW(L,NY,NX)=ZOHW(L,NY,NX)+TOHBLS(L,NY,NX)
    ZSO4W(L,NY,NX)=ZSO4W(L,NY,NX)+TSOBLS(L,NY,NX)
    ZCLW(L,NY,NX)=ZCLW(L,NY,NX)+TCLBLS(L,NY,NX)
    ZCO3W(L,NY,NX)=ZCO3W(L,NY,NX)+TC3BLS(L,NY,NX)
    ZHCO3W(L,NY,NX)=ZHCO3W(L,NY,NX)+THCBLS(L,NY,NX)
    ZALH1W(L,NY,NX)=ZALH1W(L,NY,NX)+TAL1BS(L,NY,NX)
    ZALH2W(L,NY,NX)=ZALH2W(L,NY,NX)+TAL2BS(L,NY,NX)
    ZALH3W(L,NY,NX)=ZALH3W(L,NY,NX)+TAL3BS(L,NY,NX)
    ZALH4W(L,NY,NX)=ZALH4W(L,NY,NX)+TAL4BS(L,NY,NX)
    ZALSW(L,NY,NX)=ZALSW(L,NY,NX)+TALSBS(L,NY,NX)
    ZFEH1W(L,NY,NX)=ZFEH1W(L,NY,NX)+TFE1BS(L,NY,NX)
    ZFEH2W(L,NY,NX)=ZFEH2W(L,NY,NX)+TFE2BS(L,NY,NX)
    ZFEH3W(L,NY,NX)=ZFEH3W(L,NY,NX)+TFE3BS(L,NY,NX)
    ZFEH4W(L,NY,NX)=ZFEH4W(L,NY,NX)+TFE4BS(L,NY,NX)
    ZFESW(L,NY,NX)=ZFESW(L,NY,NX)+TFESBS(L,NY,NX)
    ZCAOW(L,NY,NX)=ZCAOW(L,NY,NX)+TCAOBS(L,NY,NX)
    ZCACW(L,NY,NX)=ZCACW(L,NY,NX)+TCACBS(L,NY,NX)
    ZCAHW(L,NY,NX)=ZCAHW(L,NY,NX)+TCAHBS(L,NY,NX)
    ZCASW(L,NY,NX)=ZCASW(L,NY,NX)+TCASBS(L,NY,NX)
    ZMGOW(L,NY,NX)=ZMGOW(L,NY,NX)+TMGOBS(L,NY,NX)
    ZMGCW(L,NY,NX)=ZMGCW(L,NY,NX)+TMGCBS(L,NY,NX)
    ZMGHW(L,NY,NX)=ZMGHW(L,NY,NX)+TMGHBS(L,NY,NX)
    ZMGSW(L,NY,NX)=ZMGSW(L,NY,NX)+TMGSBS(L,NY,NX)
    ZNACW(L,NY,NX)=ZNACW(L,NY,NX)+TNACBS(L,NY,NX)
    ZNASW(L,NY,NX)=ZNASW(L,NY,NX)+TNASBS(L,NY,NX)
    ZKASW(L,NY,NX)=ZKASW(L,NY,NX)+TKASBS(L,NY,NX)
    H0PO4W(L,NY,NX)=H0PO4W(L,NY,NX)+TH0PBS(L,NY,NX)
    H3PO4W(L,NY,NX)=H3PO4W(L,NY,NX)+TH3PBS(L,NY,NX)
    ZFE1PW(L,NY,NX)=ZFE1PW(L,NY,NX)+TF1PBS(L,NY,NX)
    ZFE2PW(L,NY,NX)=ZFE2PW(L,NY,NX)+TF2PBS(L,NY,NX)
    ZCA0PW(L,NY,NX)=ZCA0PW(L,NY,NX)+TC0PBS(L,NY,NX)
    ZCA1PW(L,NY,NX)=ZCA1PW(L,NY,NX)+TC1PBS(L,NY,NX)
    ZCA2PW(L,NY,NX)=ZCA2PW(L,NY,NX)+TC2PBS(L,NY,NX)
    ZMG1PW(L,NY,NX)=ZMG1PW(L,NY,NX)+TM1PBS(L,NY,NX)
  ENDIF

  end subroutine UpdateSoluteInSnow

end module SnowBalMod
